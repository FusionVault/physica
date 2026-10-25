//! Potential fields on a level grid, in the wavenumber domain: upward continuation, reduction to
//! the pole and pseudo-gravity (Blakely, *Potential Theory in Gravity and Magnetic
//! Applications*, 1995, ch. 12).
//!
//! With `z` down and the sources below the grid, a field harmonic above them has a spectrum that
//! falls as `e^(−|k|h)` with height `h`. A derivative along a unit vector `û` multiplies it by
//! `|k|Θ_u`, with `Θ_u = û_z + i(û_x k_x + û_y k_y)/|k|`. The transforms run on
//! [`crate::fft::fft`], which treats the grid as periodic. So an anomaly should fall off well
//! inside the grid, or be padded first.
//!
//! ```
//! use geopotential::planar::Grid;
//!
//! // A point mass 50 m below the centre of a 128 × 128 grid at 10 m, seen 30 m higher.
//! let grid = Grid { nx: 128, ny: 128, dx: 10.0, dy: 10.0 };
//! let gz = |depth: f64| grid.points().map(|(x, y)| depth / ((x - 640.0).powi(2) + (y - 640.0).powi(2) + depth * depth).powf(1.5)).collect::<Vec<_>>();
//! let up = grid.upward_continue(&gz(50.0), 30.0);
//! let centre = 64 * 128 + 64;
//! assert!((up[centre] / gz(80.0)[centre] - 1.0).abs() < 0.03); // periodic wrap-around costs ~2%
//! ```

use crate::fft::fft;
use crate::linalg::{cholesky, cholesky_solve};
use metrology::{G, MU_0};

/// A level grid of `nx × ny` values at spacings `dx` and `dy` (m). Value `i·ny + j` is at
/// `(i·dx, j·dy)`, with `x` north and `y` east. Both sizes must be powers of two.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Grid {
    pub nx: usize,
    pub ny: usize,
    pub dx: f64,
    pub dy: f64,
}

impl Grid {
    /// The grid's points `(x, y)`, in the order of its values.
    pub fn points(&self) -> impl Iterator<Item = (f64, f64)> + '_ {
        (0..self.nx).flat_map(move |i| (0..self.ny).map(move |j| (i as f64 * self.dx, j as f64 * self.dy)))
    }

    /// Multiplies the spectrum of `values` by `kernel(k_x, k_y) = (re, im)`, with wavenumbers in
    /// rad/m, and returns the real part of the result.
    pub fn filter(&self, values: &[f64], kernel: impl Fn(f64, f64) -> (f64, f64)) -> Vec<f64> {
        assert_eq!(values.len(), self.nx * self.ny, "one value per grid point");
        let (mut re, mut im) = (values.to_vec(), vec![0.0; values.len()]);
        self.fft2(&mut re, &mut im, false);
        for i in 0..self.nx {
            let kx = wavenumber(i, self.nx, self.dx);
            for j in 0..self.ny {
                let (a, b) = kernel(kx, wavenumber(j, self.ny, self.dy));
                let k = i * self.ny + j;
                (re[k], im[k]) = (re[k] * a - im[k] * b, re[k] * b + im[k] * a);
            }
        }
        self.fft2(&mut re, &mut im, true);
        re
    }

    /// The field `h` metres higher: its spectrum times `e^(−|k|h)`. A negative `h` continues it
    /// downward, which amplifies short wavelengths and soon noise.
    pub fn upward_continue(&self, values: &[f64], h: f64) -> Vec<f64> {
        self.filter(values, |kx, ky| ((-kx.hypot(ky) * h).exp(), 0.0))
    }

    /// A total-field anomaly reduced to the pole: as if the magnetisation `m` and the ambient
    /// field `f` (unit vectors, `z` down) were both vertical. The spectrum is divided by
    /// `Θ_m Θ_f` and the mean set to zero. The filter's gain grows as `1/sin²I` toward the
    /// magnetic equator, so it amplifies noise at low inclination.
    pub fn reduce_to_pole(&self, values: &[f64], m: [f64; 3], f: [f64; 3]) -> Vec<f64> {
        self.filter(values, |kx, ky| pole_kernel(kx, ky, m, f, 1.0))
    }

    /// Baranov's pseudo-gravity: the vertical attraction (m/s², positive down) that the body
    /// behind a total-field anomaly would have if its density contrast were
    /// `density_per_magnetisation` (kg/m³ per A/m) times its magnetisation. By Poisson's relation
    /// the anomaly reduced to the pole is `(μ₀/4π)(M/GΔρ)∂g_z/∂z`, so the spectrum is reduced
    /// and divided by `|k|`. The mean is set to zero.
    pub fn pseudo_gravity(&self, values: &[f64], m: [f64; 3], f: [f64; 3], density_per_magnetisation: f64) -> Vec<f64> {
        let scale = G * density_per_magnetisation * 4.0 * core::f64::consts::PI / MU_0;
        self.filter(values, |kx, ky| pole_kernel(kx, ky, m, f, scale / kx.hypot(ky)))
    }

    /// The 2-D FFT, rows then columns.
    fn fft2(&self, re: &mut [f64], im: &mut [f64], inverse: bool) {
        let (nx, ny) = (self.nx, self.ny);
        for i in 0..nx {
            fft(&mut re[i * ny..(i + 1) * ny], &mut im[i * ny..(i + 1) * ny], inverse);
        }
        let (mut cr, mut ci) = (vec![0.0; nx], vec![0.0; nx]);
        for j in 0..ny {
            for i in 0..nx {
                (cr[i], ci[i]) = (re[i * ny + j], im[i * ny + j]);
            }
            fft(&mut cr, &mut ci, inverse);
            for i in 0..nx {
                (re[i * ny + j], im[i * ny + j]) = (cr[i], ci[i]);
            }
        }
    }
}

/// Euler deconvolution of a profile (Thompson 1982; Reid & Thurston 2014): the source position
/// `[x₀, z₀]` and background `B` that best satisfy Euler's homogeneity equation
/// `(x − x₀)∂f/∂x + (z − z₀)∂f/∂z = −N(f − B)` at the points `(x, z)`, `z` down, by least
/// squares, returned as `[x₀, z₀, B]`.
///
/// The structural index `N` is the field's rate of fall-off from its source: 3 for a sphere's
/// (a dipole's) magnetic anomaly, 2 for a horizontal cylinder's magnetic anomaly or a sphere's
/// gravity, 1 for a cylinder's gravity. `fz` can come from the grid as
/// `filter(values, |kx, ky| (kx.hypot(ky), 0.0))`.
///
/// # Panics
/// If the slices differ in length, or the points cannot fix three unknowns.
pub fn euler(x: &[f64], z: &[f64], f: &[f64], fx: &[f64], fz: &[f64], n: f64) -> [f64; 3] {
    let len = x.len();
    assert!([z.len(), f.len(), fx.len(), fz.len()].iter().all(|&l| l == len), "one value per point");
    // Rows [f_x, f_z, N]·[x₀, z₀, B] = x f_x + z f_z + N f, solved by normal equations whose
    // columns are scaled to unit norm first (f_x and N differ by orders of magnitude).
    let row = |i: usize| [fx[i], fz[i], n];
    let scale: [f64; 3] = core::array::from_fn(|c| 1.0 / (0..len).map(|i| row(i)[c].powi(2)).sum::<f64>().sqrt());
    let (mut normal, mut rhs) = ([0.0; 9], [0.0; 3]);
    for i in 0..len {
        let (a, b) = (row(i), x[i] * fx[i] + z[i] * fz[i] + n * f[i]);
        for r in 0..3 {
            rhs[r] += a[r] * scale[r] * b;
            for c in 0..=r {
                normal[r * 3 + c] += a[r] * scale[r] * a[c] * scale[c];
            }
        }
    }
    cholesky(&mut normal, 3);
    let y = cholesky_solve(&normal, 3, &rhs);
    core::array::from_fn(|c| y[c] * scale[c])
}

/// The FFT's wavenumber (rad/m) at index `i` of `n` samples `d` apart, negative above `n/2`.
fn wavenumber(i: usize, n: usize, d: f64) -> f64 {
    let i = if i < n / 2 { i as f64 } else { i as f64 - n as f64 };
    core::f64::consts::TAU * i / (n as f64 * d)
}

/// `scale/(Θ_m Θ_f)`, and zero at `k = 0`.
fn pole_kernel(kx: f64, ky: f64, m: [f64; 3], f: [f64; 3], scale: f64) -> (f64, f64) {
    let k = kx.hypot(ky);
    if k == 0.0 {
        return (0.0, 0.0);
    }
    let theta = |u: [f64; 3]| (u[2], (u[0] * kx + u[1] * ky) / k);
    let ((a, b), (c, d)) = (theta(m), theta(f));
    let (re, im) = (a * c - b * d, a * d + b * c);
    let norm = scale / (re * re + im * im);
    (re * norm, -im * norm)
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::f64::consts::PI;

    const GRID: Grid = Grid { nx: 128, ny: 128, dx: 10.0, dy: 10.0 };

    /// `f̂·∇∇(1/r)·m̂` of a dipole at depth `depth` below the grid's centre: a total-field anomaly
    /// up to its constant.
    fn dipole(depth: f64, m: [f64; 3], f: [f64; 3]) -> Vec<f64> {
        GRID.points()
            .map(|(x, y)| {
                let d = [x - 640.0, y - 640.0, -depth];
                let r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
                let (dm, df, mf) = ((0..3).map(|i| d[i] * m[i]).sum::<f64>(), (0..3).map(|i| d[i] * f[i]).sum::<f64>(), (0..3).map(|i| m[i] * f[i]).sum::<f64>());
                (3.0 * dm * df - r2 * mf) / (r2 * r2 * r2.sqrt())
            })
            .collect()
    }

    fn rms_relative(a: &[f64], b: &[f64]) -> f64 {
        (a.iter().zip(b).map(|(x, y)| (x - y).powi(2)).sum::<f64>() / b.iter().map(|y| y * y).sum::<f64>()).sqrt()
    }

    fn unit(inclination: f64, declination: f64) -> [f64; 3] {
        let (i, d) = (inclination.to_radians(), declination.to_radians());
        [i.cos() * d.cos(), i.cos() * d.sin(), i.sin()]
    }

    #[test]
    fn upward_continuation() {
        // A pole anomaly 50 m down, continued 30 m: 6.9×10⁻⁴ rms against the anomaly 80 m down,
        // the notes' value, all of it from the periodic grid.
        let pole = [0.0, 0.0, 1.0];
        let up = GRID.upward_continue(&dipole(50.0, pole, pole), 30.0);
        let err = rms_relative(&up, &dipole(80.0, pole, pole));
        assert!((err / 6.9e-4 - 1.0).abs() < 0.02, "{err:e}");
        // Up and back down is the identity, to round-off times the gain e^(|k|h) of the way down:
        // 9 at 5 m on this grid, 6×10⁵ at 30 m.
        let field = dipole(50.0, pole, pole);
        for (h, tol) in [(5.0, 1e-14), (30.0, 2e-11)] {
            let err = rms_relative(&GRID.upward_continue(&GRID.upward_continue(&field, h), -h), &field);
            assert!(err < tol, "{h} m: {err:e}");
        }
    }

    #[test]
    fn reduction_to_the_pole() {
        // At I = 45° the reduced anomaly is the pole anomaly to 1.1% rms (the notes' value); with
        // Θ's imaginary part of the wrong sign it would be 168% off.
        let (tilted, pole) = (unit(45.0, 0.0), [0.0, 0.0, 1.0]);
        let reduced = GRID.reduce_to_pole(&dipole(50.0, tilted, tilted), tilted, tilted);
        let err = rms_relative(&reduced, &dipole(50.0, pole, pole));
        assert!((err / 0.011 - 1.0).abs() < 0.1, "{err}");
        let flipped = [-tilted[0], -tilted[1], tilted[2]];
        assert!(rms_relative(&GRID.reduce_to_pole(&dipole(50.0, tilted, tilted), flipped, flipped), &dipole(50.0, pole, pole)) > 1.0);
        // At the pole the filter only removes the mean.
        let field = dipole(50.0, pole, pole);
        let mean = field.iter().sum::<f64>() / field.len() as f64;
        let same = GRID.reduce_to_pole(&field, pole, pole);
        assert!(same.iter().zip(&field).all(|(a, b)| (a - (b - mean)).abs() < 1e-14 * field[64 * 128 + 64]));
    }

    /// Value and derivatives of `field(x, z)` along `z = 0`, by fourth-order 5-point stencils.
    fn sampled(field: impl Fn(f64, f64) -> f64, xs: &[f64]) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let h = 1e-2;
        let d = |g: &dyn Fn(f64) -> f64| (g(-2.0 * h) - 8.0 * g(-h) + 8.0 * g(h) - g(2.0 * h)) / (12.0 * h);
        let f = xs.iter().map(|&x| field(x, 0.0)).collect();
        let fx = xs.iter().map(|&x| d(&|e| field(x + e, 0.0))).collect();
        let fz = xs.iter().map(|&x| d(&|e| field(x, e))).collect();
        (f, fx, fz)
    }

    #[test]
    fn euler_deconvolution() {
        // Sources at x₀ = 3 m, z₀ = 20 m (z down), seen along z = 0. Each field is homogeneous of
        // degree −N in (x − x₀, z − z₀), the structural indices of Reid & Thurston (2014).
        let (x0, z0) = (3.0, 20.0);
        let tilted = unit(60.0, 0.0);
        let dipole = move |x: f64, z: f64| {
            let d = [x - x0, 0.0, z - z0];
            let r2 = d[0] * d[0] + d[2] * d[2];
            let dm = d[0] * tilted[0] + d[2] * tilted[2];
            (3.0 * dm * dm - r2) / (r2 * r2 * r2.sqrt())
        };
        let line_dipole = move |x: f64, z: f64| {
            let (dx, dz) = (x - x0, z - z0);
            let r2 = dx * dx + dz * dz;
            ((dz * dz - dx * dx) * (tilted[2].powi(2) - tilted[0].powi(2)) + 4.0 * dx * dz * tilted[0] * tilted[2]) / (r2 * r2)
        };
        let point_mass = move |x: f64, z: f64| (z0 - z) / ((x - x0).powi(2) + (z0 - z).powi(2)).powf(1.5);
        let line_mass = move |x: f64, z: f64| (z0 - z) / ((x - x0).powi(2) + (z0 - z).powi(2));
        let xs: Vec<f64> = (0..41).map(|i| 2.0 * i as f64 - 40.0).collect();
        let zs = vec![0.0; xs.len()];
        type Profile<'a> = &'a dyn Fn(f64, f64) -> f64;
        let cases: [(Profile, f64); 4] = [(&dipole, 3.0), (&line_dipole, 2.0), (&point_mass, 2.0), (&line_mass, 1.0)];
        for (k, (field, n)) in cases.into_iter().enumerate() {
            // A background of 1/10 of the peak, which the solve should return.
            let (f, fx, fz) = sampled(field, &xs);
            let background = 0.1 * f.iter().fold(0f64, |m, v| m.max(v.abs()));
            let f: Vec<f64> = f.iter().map(|v| v + background).collect();
            let [x, z, b] = euler(&xs, &zs, &f, &fx, &fz, n);
            // 2×10⁻¹¹ m at worst, set by the stencils (the notes' second-order ones gave 1.2×10⁻⁹).
            assert!((x - x0).abs() < 1e-10 && (z - z0).abs() < 1e-10, "case {k}: {x} {z}");
            assert!((b / background - 1.0).abs() < 1e-11, "case {k}: background {b:e}");
        }
        // The wrong index biases the depth: a dipole read as a point mass looks shallower.
        let (f, fx, fz) = sampled(dipole, &xs);
        assert!(euler(&xs, &zs, &f, &fx, &fz, 2.0)[1] < 0.9 * z0);
    }

    #[test]
    fn pseudo_gravity_is_poisson() {
        // A 10 m sphere 50 m down, at 1 A/m and 2700 kg/m³: its pseudo-gravity is its gravity,
        // peak minus corner, to −0.13% (the notes' value; the periodic grid again).
        let (a, depth, rho) = (10.0f64, 50.0, 2700.0);
        let volume = 4.0 / 3.0 * PI * a.powi(3);
        let pole = [0.0, 0.0, 1.0];
        let anomaly: Vec<f64> = dipole(depth, pole, pole).iter().map(|t| MU_0 / (4.0 * PI) * volume * t).collect();
        let gz: Vec<f64> = GRID.points().map(|(x, y)| G * rho * volume * depth / ((x - 640.0).powi(2) + (y - 640.0).powi(2) + depth * depth).powf(1.5)).collect();
        let pseudo = GRID.pseudo_gravity(&anomaly, pole, pole, rho);
        let centre = 64 * 128 + 64;
        let err = (pseudo[centre] - pseudo[0]) / (gz[centre] - gz[0]) - 1.0;
        assert!((err / -1.3e-3 - 1.0).abs() < 0.05, "{err:e}");
    }
}
