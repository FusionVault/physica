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
use crate::linalg::symmetric_eigen;
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

    /// `values` on a grid twice the size, the data in its middle and each edge value repeated
    /// outward over the margin, with a cosine roll-off to zero across the margin if `taper`. The
    /// FFT treats a grid as periodic, so an anomaly that reaches an edge wraps round; padding this
    /// way cut the error of continuing a dipole 300 m deep by 30 m on a 128 × 128 grid from
    /// 1.8×10⁻² to 1.1×10⁻³ rms, where padding with zeros helps little (a step at the edge).
    /// Filter the padded grid, then take [`unpadded`](Self::unpadded).
    pub fn padded(&self, values: &[f64], taper: bool) -> (Grid, Vec<f64>) {
        assert_eq!(values.len(), self.nx * self.ny, "one value per grid point");
        let big = Grid { nx: 2 * self.nx, ny: 2 * self.ny, ..*self };
        let (px, py) = (self.nx / 2, self.ny / 2);
        let window = |i: usize, n: usize, pad: usize| -> f64 {
            if !taper || (pad..pad + n).contains(&i) {
                1.0
            } else {
                let out = if i < pad { (pad - i) as f64 } else { (i + 1 - pad - n) as f64 };
                0.5 + 0.5 * (core::f64::consts::PI * out / pad as f64).cos()
            }
        };
        let out = (0..big.nx)
            .flat_map(|i| {
                (0..big.ny).map(move |j| {
                    let (si, sj) = (i.saturating_sub(px).min(self.nx - 1), j.saturating_sub(py).min(self.ny - 1));
                    values[si * self.ny + sj] * window(i, self.nx, px) * window(j, self.ny, py)
                })
            })
            .collect();
        (big, out)
    }

    /// The middle of `values` on the grid [`padded`](Self::padded) made, back on this grid.
    pub fn unpadded(&self, values: &[f64]) -> Vec<f64> {
        assert_eq!(values.len(), 4 * self.nx * self.ny, "values on the padded grid");
        let (px, py) = (self.nx / 2, self.ny / 2);
        (0..self.nx).flat_map(|i| (0..self.ny).map(move |j| values[(i + px) * 2 * self.ny + j + py])).collect()
    }

    /// The field `h` metres higher: its spectrum times `e^(−|k|h)`. A negative `h` continues it
    /// downward, which amplifies short wavelengths and soon noise.
    pub fn upward_continue(&self, values: &[f64], h: f64) -> Vec<f64> {
        self.filter(values, |kx, ky| ((-kx.hypot(ky) * h).exp(), 0.0))
    }

    /// A total-field anomaly reduced to the pole: as if the magnetisation `m` and the ambient
    /// field `f` (unit vectors, `z` down) were both vertical. The spectrum is divided by
    /// `Θ_m Θ_f` and the mean set to zero. The filter's gain grows as `1/sin²I` toward the
    /// magnetic equator, so it amplifies noise at low inclination; see
    /// [`reduce_to_pole_damped`](Self::reduce_to_pole_damped) and
    /// [`reduce_to_equator`](Self::reduce_to_equator).
    pub fn reduce_to_pole(&self, values: &[f64], m: [f64; 3], f: [f64; 3]) -> Vec<f64> {
        self.redirect(values, (m, f), POLE, 0.0)
    }

    /// Reduction to the pole with Wiener damping `epsilon`: 0.02–0.03 was best for 1% noise at
    /// 10° inclination in the tests, where it brought the error from 0.84 to 0.38 rms; choosing it
    /// by cross-validation is usual. See [`redirect`](Self::redirect).
    pub fn reduce_to_pole_damped(&self, values: &[f64], m: [f64; 3], f: [f64; 3], epsilon: f64) -> Vec<f64> {
        self.redirect(values, (m, f), POLE, epsilon)
    }

    /// A total-field anomaly reduced to the equator: as if the magnetisation and the ambient
    /// field were horizontal, at their own declinations. The kernel `Θ_m′ Θ_f′/(Θ_m Θ_f)` then has
    /// a gain of at most 1 at every inclination, where reduction to the pole's reaches `1/sin²I`,
    /// so it is stable at low latitudes; the anomaly keeps the equator's shape, a low over the
    /// body between two highs.
    pub fn reduce_to_equator(&self, values: &[f64], m: [f64; 3], f: [f64; 3]) -> Vec<f64> {
        let flat = |u: [f64; 3]| {
            let h = u[0].hypot(u[1]);
            [u[0] / h, u[1] / h, 0.0]
        };
        self.redirect(values, (m, f), (flat(m), flat(f)), 0.0)
    }

    /// A total-field anomaly as it would be with the magnetisation and ambient field along
    /// `to = (m′, f′)` instead of `from = (m, f)` (unit vectors, `z` down): the spectrum times
    /// `Θ_m′ Θ_f′ · conj(t)/(|t|² + ε)`, `t = Θ_m Θ_f`, with the mean set to zero. With `epsilon`
    /// 0 that is the exact `Θ_m′Θ_f′/t`; a positive `epsilon` is Wiener damping, which leaves the
    /// kernel where `|t|² ≫ ε` and sends it to zero where the exact gain would exceed `1/√ε`.
    /// `|t| ≤ 1`, so `epsilon` is a fraction of the gain's floor.
    pub fn redirect(&self, values: &[f64], from: ([f64; 3], [f64; 3]), to: ([f64; 3], [f64; 3]), epsilon: f64) -> Vec<f64> {
        self.filter(values, |kx, ky| direction_kernel(kx, ky, from, to, 1.0, epsilon))
    }

    /// Baranov's pseudo-gravity: the vertical attraction (m/s², positive down) that the body
    /// behind a total-field anomaly would have if its density contrast were
    /// `density_per_magnetisation` (kg/m³ per A/m) times its magnetisation. By Poisson's relation
    /// the anomaly reduced to the pole is `(μ₀/4π)(M/GΔρ)∂g_z/∂z`, so the spectrum is reduced
    /// and divided by `|k|`. The mean is set to zero.
    pub fn pseudo_gravity(&self, values: &[f64], m: [f64; 3], f: [f64; 3], density_per_magnetisation: f64) -> Vec<f64> {
        let scale = G * density_per_magnetisation * 4.0 * core::f64::consts::PI / MU_0;
        self.filter(values, |kx, ky| direction_kernel(kx, ky, (m, f), POLE, scale / kx.hypot(ky), 0.0))
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
/// If the slices differ in length.
pub fn euler(x: &[f64], z: &[f64], f: &[f64], fx: &[f64], fz: &[f64], n: f64) -> [f64; 3] {
    let len = x.len();
    assert!([z.len(), f.len(), fx.len(), fz.len()].iter().all(|&l| l == len), "one value per point");
    least_squares(len, |i| [fx[i], fz[i], n], |i| x[i] * fx[i] + z[i] * fz[i] + n * f[i]).0
}

/// Euler deconvolution of gridded data (Reid et al. 1990): the source `[x₀, y₀, z₀]` and
/// background `B` that best satisfy `(x − x₀)∂f/∂x + (y − y₀)∂f/∂y + (z − z₀)∂f/∂z = −N(f − B)`
/// at `points` `(x, y, z)`, `z` down, by least squares, with their covariance
/// `σ̂₀²(AᵀA)⁻¹`, `σ̂₀² = ‖residual‖²/(N_points − 4)`, which only holds if the derivatives are
/// exact and is what ranks windows. A two-dimensional source (a long cylinder or contact)
/// leaves the position along its strike undetermined: the solution is then the one nearest
/// the origin along that direction, with an infinite variance there.
///
/// # Panics
/// If the slices differ in length.
pub fn euler3(points: &[[f64; 3]], f: &[f64], fx: &[f64], fy: &[f64], fz: &[f64], n: f64) -> ([f64; 4], [[f64; 4]; 4]) {
    let len = points.len();
    assert!([f.len(), fx.len(), fy.len(), fz.len()].iter().all(|&l| l == len), "one value per point");
    least_squares(len, |i| [fx[i], fy[i], fz[i], n], |i| points[i][0] * fx[i] + points[i][1] * fy[i] + points[i][2] * fz[i] + n * f[i])
}

/// Euler deconvolution in moving windows of `size × size` grid points, every `step` points,
/// over `values` observed on `grid` at `z = 0` with the derivatives `fx`, `fy`, `fz`: the
/// solutions of [`euler3`] whose depth variance is among the smallest `keep` fraction, best
/// first, each as `[x₀, y₀, z₀, B]`. Reid et al.'s (1990) way of mapping many sources; windows
/// without a source give deep, scattered, high-variance solutions that the ranking drops.
pub fn euler_windows(grid: &Grid, values: &[f64], fx: &[f64], fy: &[f64], fz: &[f64], n: f64, size: usize, step: usize, keep: f64) -> Vec<[f64; 4]> {
    let mut solutions = Vec::new();
    for i0 in (0..=grid.nx - size).step_by(step.max(1)) {
        for j0 in (0..=grid.ny - size).step_by(step.max(1)) {
            let index: Vec<usize> = (i0..i0 + size).flat_map(|i| (j0..j0 + size).map(move |j| i * grid.ny + j)).collect();
            let points: Vec<[f64; 3]> = index.iter().map(|&k| [(k / grid.ny) as f64 * grid.dx, (k % grid.ny) as f64 * grid.dy, 0.0]).collect();
            let pick = |v: &[f64]| index.iter().map(|&k| v[k]).collect::<Vec<_>>();
            let (p, cov) = euler3(&points, &pick(values), &pick(fx), &pick(fy), &pick(fz), n);
            solutions.push((cov[2][2], p));
        }
    }
    solutions.sort_by(|a, b| a.0.total_cmp(&b.0));
    solutions.truncate(((solutions.len() as f64 * keep).round() as usize).max(1));
    solutions.into_iter().map(|(_, p)| p).collect()
}

/// Euler inversion ([Uieda et al. 2025](https://doi.org/10.1093/gji/ggaf114), Algorithm 1): the
/// source `[x₀, y₀, z₀]` and background `B` fitted together with predicted data that satisfy
/// Euler's equation exactly, by Gauss–Newton on the merit `‖W^½(d_obs − d)‖ + ν‖Euler(d, p)‖`
/// from the deconvolution's solution, with the paper's weights (1, 0.1, 0.1, 0.025) on `f` and
/// its three derivatives, `ν = 0.1`, at most 20 iterations, stopping at a 10% change of the
/// merit or when a step would raise it. It returns the solution and the weighted rms misfit
/// `√(Σ W(d_obs − d)²/4N)`, which is smallest at the right structural index: run it for `N` of
/// 1, 2 and 3 and keep the best (their Algorithm 2). Noisy derivatives bias the deconvolution's
/// depth badly (−1553 m for a source at −3000 in the paper's test); the inversion gave −2663.
///
/// # Panics
/// If the slices differ in length.
pub fn euler_inversion(points: &[[f64; 3]], f: &[f64], fx: &[f64], fy: &[f64], fz: &[f64], n: f64) -> ([f64; 4], f64) {
    let len = points.len();
    assert!([f.len(), fx.len(), fy.len(), fz.len()].iter().all(|&l| l == len), "one value per point");
    const WEIGHTS: [f64; 4] = [1.0, 0.1, 0.1, 0.025];
    const NU: f64 = 0.1;
    let observed: Vec<[f64; 4]> = (0..len).map(|i| [f[i], fx[i], fy[i], fz[i]]).collect();
    let mut predicted: Vec<[f64; 4]> = observed.iter().map(|d| d.map(|x| 0.9 * x)).collect();
    let mut p = euler3(points, f, fx, fy, fz, n).0;
    let euler = |d: &[f64; 4], q: &[f64; 4], r: &[f64; 3]| (r[0] - q[0]) * d[1] + (r[1] - q[1]) * d[2] + (r[2] - q[2]) * d[3] + n * (d[0] - q[3]);
    let merit = |pred: &[[f64; 4]], q: &[f64; 4]| {
        let misfit: f64 = observed.iter().zip(pred).map(|(o, d)| (0..4).map(|k| (WEIGHTS[k] * (o[k] - d[k])).powi(2)).sum::<f64>()).sum::<f64>().sqrt();
        let residual: f64 = pred.iter().zip(points).map(|(d, r)| euler(d, q, r).powi(2)).sum::<f64>().sqrt();
        misfit + NU * residual
    };
    let mut m0 = merit(&predicted, &p);
    for _ in 0..20 {
        // Each point's Euler equation is linear in its data with gradient B_k = [N, x − x₀, y − y₀,
        // z − z₀] and in the parameters with gradient A = −[fx, fy, fz, N]; eliminating the data
        // update gives a 4×4 system for the parameter step, then the data step (Uieda et al.'s
        // Eqs. 17–19).
        let (mut normal, mut rhs) = ([[0.0; 4]; 4], [0.0; 4]);
        let rows: Vec<([f64; 4], [f64; 4], f64, f64)> = predicted
            .iter()
            .zip(&observed)
            .zip(points)
            .map(|((d, o), r)| {
                let a = [-d[1], -d[2], -d[3], -n];
                let b = [n, r[0] - p[0], r[1] - p[1], r[2] - p[2]];
                let q_inv = 1.0 / (0..4).map(|k| b[k] * b[k] / WEIGHTS[k]).sum::<f64>();
                let br: f64 = (0..4).map(|k| b[k] * (o[k] - d[k])).sum();
                (a, b, q_inv, euler(d, &p, r) + br)
            })
            .collect();
        for (a, _, q_inv, e) in &rows {
            for i in 0..4 {
                rhs[i] -= a[i] * q_inv * e;
                for j in 0..4 {
                    normal[i][j] += a[i] * q_inv * a[j];
                }
            }
        }
        let dp = pseudo_solve(&normal, &rhs);
        let next: Vec<[f64; 4]> = rows
            .iter()
            .zip(&predicted)
            .zip(&observed)
            .map(|(((a, b, q_inv, e), d), o)| {
                let lambda = q_inv * (e + (0..4).map(|i| a[i] * dp[i]).sum::<f64>());
                core::array::from_fn(|k| o[k] - b[k] * lambda / WEIGHTS[k] + (d[k] - d[k]))
            })
            .collect();
        let q: [f64; 4] = core::array::from_fn(|i| p[i] + dp[i]);
        let m1 = merit(&next, &q);
        if m1 > m0 {
            break;
        }
        (p, predicted) = (q, next);
        if (m0 - m1).abs() / m0 < 0.1 {
            break;
        }
        m0 = m1;
    }
    let wrmse = (observed.iter().zip(&predicted).map(|(o, d)| (0..4).map(|k| WEIGHTS[k] * (o[k] - d[k]).powi(2)).sum::<f64>()).sum::<f64>() / (4 * len) as f64).sqrt();
    (p, wrmse)
}

/// The least-squares solution of `row(i)·y = rhs(i)` over `len` rows, and its covariance
/// `σ̂₀²(AᵀA)⁻¹`. The columns are scaled to unit norm first (a derivative and a structural index
/// differ by orders of magnitude), and the normal matrix is inverted through its eigenvalues,
/// dropping those below `10⁻¹²` of the largest, so a rank-deficient fit returns the minimum-norm
/// solution with infinite variance along the missing directions.
fn least_squares<const K: usize>(len: usize, row: impl Fn(usize) -> [f64; K], rhs: impl Fn(usize) -> f64) -> ([f64; K], [[f64; K]; K]) {
    let scale: [f64; K] = core::array::from_fn(|c| 1.0 / (0..len).map(|i| row(i)[c].powi(2)).sum::<f64>().sqrt());
    let (mut normal, mut b) = ([[0.0; K]; K], [0.0; K]);
    for i in 0..len {
        let (a, r) = (row(i), rhs(i));
        for p in 0..K {
            b[p] += a[p] * scale[p] * r;
            for q in 0..K {
                normal[p][q] += a[p] * scale[p] * a[q] * scale[q];
            }
        }
    }
    let (inverse, undetermined, rank) = pseudo_inverse(&normal);
    let y: [f64; K] = core::array::from_fn(|p| (0..K).map(|q| inverse[p][q] * b[q]).sum::<f64>() * scale[p]);
    let residual: f64 = (0..len).map(|i| (rhs(i) - (0..K).map(|c| row(i)[c] * y[c]).sum::<f64>()).powi(2)).sum();
    let sigma2 = if len > rank { residual / (len - rank) as f64 } else { 0.0 };
    (y, core::array::from_fn(|p| core::array::from_fn(|q| if undetermined[p][q] { f64::INFINITY } else { sigma2 * inverse[p][q] * scale[p] * scale[q] })))
}

/// The pseudo-inverse of a symmetric matrix through its eigenvalues, which entries the null
/// space (eigenvalues below `10⁻¹²` of the largest) reaches, and the rank.
fn pseudo_inverse<const K: usize>(a: &[[f64; K]; K]) -> ([[f64; K]; K], [[bool; K]; K], usize) {
    let flat: Vec<f64> = a.iter().flatten().copied().collect();
    let (values, vectors) = symmetric_eigen(&flat, K);
    let cut = 1e-12 * values.iter().fold(0f64, |m, v| m.max(v.abs()));
    let rank = values.iter().filter(|&&v| v > cut).count();
    let inv = core::array::from_fn(|p| core::array::from_fn(|q| (0..K).filter(|&k| values[k] > cut).map(|k| vectors[p * K + k] * vectors[q * K + k] / values[k]).sum::<f64>()));
    let null = core::array::from_fn(|p| core::array::from_fn(|q| (0..K).any(|k| values[k] <= cut && (vectors[p * K + k] * vectors[q * K + k]).abs() > 1e-12)));
    (inv, null, rank)
}

/// `A⁺b` for a symmetric `A`, the minimum-norm solution when `A` is singular.
fn pseudo_solve<const K: usize>(a: &[[f64; K]; K], b: &[f64; K]) -> [f64; K] {
    let (inverse, _, _) = pseudo_inverse(a);
    core::array::from_fn(|p| (0..K).map(|q| inverse[p][q] * b[q]).sum())
}

/// The FFT's wavenumber (rad/m) at index `i` of `n` samples `d` apart, negative above `n/2`.
fn wavenumber(i: usize, n: usize, d: f64) -> f64 {
    let i = if i < n / 2 { i as f64 } else { i as f64 - n as f64 };
    core::f64::consts::TAU * i / (n as f64 * d)
}

/// Both directions vertical: the pole.
const POLE: ([f64; 3], [f64; 3]) = ([0.0, 0.0, 1.0], [0.0, 0.0, 1.0]);

/// `scale · Θ_m′Θ_f′ · conj(t)/(|t|² + ε)` with `t = Θ_m Θ_f`: the spectrum of an anomaly seen
/// with magnetisation and field along `from = (m, f)` turned into the one along `to`, damped by
/// `epsilon` (0 for the exact `1/t`), and zero at `k = 0`.
fn direction_kernel(kx: f64, ky: f64, from: ([f64; 3], [f64; 3]), to: ([f64; 3], [f64; 3]), scale: f64, epsilon: f64) -> (f64, f64) {
    let k = kx.hypot(ky);
    if k == 0.0 {
        return (0.0, 0.0);
    }
    let theta = |u: [f64; 3]| (u[2], (u[0] * kx + u[1] * ky) / k);
    let product = |(a, b): (f64, f64), (c, d): (f64, f64)| (a * c - b * d, a * d + b * c);
    let (re, im) = product(theta(from.0), theta(from.1));
    let (nr, ni) = product(theta(to.0), theta(to.1));
    let norm = scale / (re * re + im * im + epsilon);
    product((nr, ni), (re * norm, -im * norm))
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
    fn padding() {
        // A dipole 300 m deep reaches the grid's edges: continued 30 m, the periodic wrap costs
        // 1.8×10⁻² rms; the edge-extended, tapered pad to 256 × 256 leaves 1.1×10⁻³ (the notes'
        // values). For the 50 m dipole, which decays inside the grid, 6.9×10⁻⁴ becomes 2.5×10⁻⁴.
        let pole = [0.0, 0.0, 1.0];
        let continued = |depth: f64, pad: bool| {
            let field = dipole(depth, pole, pole);
            let up = if pad {
                let (big, values) = GRID.padded(&field, true);
                GRID.unpadded(&big.upward_continue(&values, 30.0))
            } else {
                GRID.upward_continue(&field, 30.0)
            };
            rms_relative(&up, &dipole(depth + 30.0, pole, pole))
        };
        assert!((continued(300.0, false) / 1.8e-2 - 1.0).abs() < 0.1 && (continued(300.0, true) / 1.1e-3 - 1.0).abs() < 0.1, "{} {}", continued(300.0, false), continued(300.0, true));
        assert!((continued(50.0, true) / 2.5e-4 - 1.0).abs() < 0.1);
        // The pad repeats the edge values and the taper ends at zero; unpadded is the inverse.
        let field = dipole(300.0, pole, pole);
        let (big, values) = GRID.padded(&field, false);
        assert_eq!((big.nx, big.ny, values.len()), (256, 256, 65536));
        assert_eq!(values[0], field[0]);
        assert_eq!(values[64 * 256 + 300 - 64], field[0 * 128 + 127]);
        assert!(GRID.unpadded(&values) == field && GRID.unpadded(&GRID.padded(&field, true).1) == field);
        let tapered = GRID.padded(&field, true).1;
        assert!(tapered[0].abs() < 1e-12 * field[64 * 128 + 64] && tapered[64 * 256 + 64] == field[0]);
    }

    #[test]
    fn low_inclination() {
        // A dipole 100 m deep on a 64 × 64 grid at 25 m, induced at I = 10°, padded. Reduced to
        // the equator the clean data give the true equatorial anomaly to 1%; reduced to the pole
        // they amplify 1% noise to an error of order 1, which Wiener damping at ε = 0.025 cuts by
        // more than half; the damped kernel is the plain one as ε → 0.
        let grid = Grid { nx: 64, ny: 64, dx: 25.0, dy: 25.0 };
        let field = |m: [f64; 3], f: [f64; 3]| -> Vec<f64> {
            grid.points()
                .map(|(x, y)| {
                    let d = [x - 800.0, y - 800.0, -100.0];
                    let r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
                    let (dm, df, mf) = ((0..3).map(|i| d[i] * m[i]).sum::<f64>(), (0..3).map(|i| d[i] * f[i]).sum::<f64>(), (0..3).map(|i| m[i] * f[i]).sum::<f64>());
                    (3.0 * dm * df - r2 * mf) / (r2 * r2 * r2.sqrt())
                })
                .collect()
        };
        let apply = |values: &[f64], f: &dyn Fn(&Grid, &[f64]) -> Vec<f64>| {
            let (big, padded) = grid.padded(values, true);
            grid.unpadded(&f(&big, &padded))
        };
        let (low, equator, pole) = (unit(10.0, 0.0), unit(0.0, 0.0), [0.0, 0.0, 1.0]);
        let (data, at_equator, at_pole) = (field(low, low), field(equator, equator), field(pole, pole));
        let rte = apply(&data, &|g, v| g.reduce_to_equator(v, low, low));
        assert!(rms_relative(&rte, &at_equator) < 0.01, "{}", rms_relative(&rte, &at_equator));
        // Round trip: the pole anomaly turned back to I = 10° is the data less its mean, to the
        // 10⁻⁵ the Nyquist rows carry (a complex kernel there has no real inverse transform). The
        // equator cannot be turned back from: Θ_m′Θ_f′ vanishes across the declination.
        let back = grid.redirect(&grid.reduce_to_pole(&data, low, low), POLE, (low, low), 0.0);
        let mean = data.iter().sum::<f64>() / data.len() as f64;
        let round_trip = rms_relative(&back.iter().map(|v| v + mean).collect::<Vec<_>>(), &data);
        assert!(round_trip < 1e-4, "{round_trip:e}");
        // Noise of 1% of the peak, from a fixed generator (sums of twelve uniforms, so Gaussian).
        let mut state = 7u64;
        let peak = data.iter().fold(0f64, |m, v| m.max(v.abs()));
        let noisy: Vec<f64> = data
            .iter()
            .map(|v| {
                let g: f64 = (0..12)
                    .map(|_| {
                        state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                        (state >> 11) as f64 / (1u64 << 53) as f64
                    })
                    .sum::<f64>()
                    - 6.0;
                v + 0.01 * peak * g
            })
            .collect();
        let plain = rms_relative(&apply(&noisy, &|g, v| g.reduce_to_pole(v, low, low)), &at_pole);
        let damped = rms_relative(&apply(&noisy, &|g, v| g.reduce_to_pole_damped(v, low, low, 0.025)), &at_pole);
        assert!(plain > 0.5 && damped < 0.6 && damped < 0.6 * plain, "plain {plain} damped {damped}");
        let (tiny, exact) = (grid.reduce_to_pole_damped(&data, low, low, 1e-18), grid.reduce_to_pole(&data, low, low));
        assert!(rms_relative(&tiny, &exact) < 1e-12);
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
    fn euler_in_three_dimensions() {
        // The notes' noise-free test: a grid of 20 × 16 km at 250 m, 500 m above the ground
        // (z = −500, z down), sources 2 km deep at (10, 8) km, derivatives by central differences
        // of 1 m. Deconvolution at the true index is exact to the stencil (10⁻³ m); the inversion's
        // weighted misfit picks the index. A line mass (SI 1) along azimuth 30° leaves its
        // position along strike undetermined: the error across strike is what is checked.
        let (x0, y0, z0) = (10e3, 8e3, 2e3);
        let tilted = unit(60.0, 0.0);
        let dipole = move |r: [f64; 3]| {
            let d = [r[0] - x0, r[1] - y0, r[2] - z0];
            let r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
            let dm: f64 = (0..3).map(|i| d[i] * tilted[i]).sum();
            1e11 * (3.0 * dm * dm - r2) / (r2 * r2 * r2.sqrt())
        };
        let point_mass = move |r: [f64; 3]| 1e13 * (z0 - r[2]) / ((r[0] - x0).powi(2) + (r[1] - y0).powi(2) + (z0 - r[2]).powi(2)).powf(1.5);
        let (ca, sa) = (30f64.to_radians().cos(), 30f64.to_radians().sin());
        let line_mass = move |r: [f64; 3]| 1e9 * (z0 - r[2]) / (((r[0] - x0) * ca - (r[1] - y0) * sa).powi(2) + (z0 - r[2]).powi(2));
        let points: Vec<[f64; 3]> = (0..81).flat_map(|i| (0..65).map(move |j| [i as f64 * 250.0, j as f64 * 250.0, -500.0])).collect();
        type Source<'a> = &'a dyn Fn([f64; 3]) -> f64;
        let cases: [(Source, f64, bool); 3] = [(&dipole, 3.0, false), (&point_mass, 2.0, false), (&line_mass, 1.0, true)];
        for (field, n, two_d) in cases {
            let sample = |shift: [f64; 3]| -> Vec<f64> { points.iter().map(|r| field([r[0] + shift[0], r[1] + shift[1], r[2] + shift[2]]) + 50.0).collect() };
            let f = sample([0.0; 3]);
            let d = |axis: usize| -> Vec<f64> {
                let (mut up, mut down) = ([0.0; 3], [0.0; 3]);
                (up[axis], down[axis]) = (1.0, -1.0);
                sample(up).iter().zip(sample(down)).map(|(a, b)| (a - b) / 2.0).collect()
            };
            let (fx, fy, fz) = (d(0), d(1), d(2));
            let (p, cov) = euler3(&points, &f, &fx, &fy, &fz, n);
            let across = if two_d { ((p[0] - x0) * ca - (p[1] - y0) * sa).abs() } else { (p[0] - x0).hypot(p[1] - y0) };
            assert!(across < 2e-3 && (p[2] - z0).abs() < 2e-3 && (p[3] - 50.0).abs() < 1e-2, "N = {n}: {p:?}");
            assert!(cov[2][2].is_finite() && cov[2][2] >= 0.0);
            if two_d {
                // Along strike the solution is the minimum-norm one and the variance infinite.
                let along = (p[0] - x0) * sa + (p[1] - y0) * ca;
                assert!(along.abs() > 100.0 && (cov[0][0].is_infinite() || cov[1][1].is_infinite()), "{p:?} {cov:?}");
            }
            let fits = [1.0, 2.0, 3.0].map(|si| euler_inversion(&points, &f, &fx, &fy, &fz, si).1);
            let best = (0..3).min_by(|&a, &b| fits[a].total_cmp(&fits[b])).unwrap();
            assert!((best + 1) as f64 == n, "N = {n}: misfits {fits:?}");
            let (q, _) = euler_inversion(&points, &f, &fx, &fy, &fz, n);
            assert!((q[2] - z0).abs() < 1.0, "inversion at N = {n}: {q:?}");
        }
        // Windows over the dipole on the grid itself, with derivatives from the spectrum: the
        // best-ranked window puts the source within 2% of its depth and 100 m of its position.
        let grid = Grid { nx: 128, ny: 64, dx: 250.0, dy: 250.0 };
        let field: Vec<f64> = grid.points().map(|(x, y)| dipole([x, y, 0.0])).collect();
        let (big, padded) = grid.padded(&field, true);
        let spectral = |kernel: &dyn Fn(f64, f64) -> (f64, f64)| grid.unpadded(&big.filter(&padded, kernel));
        let (fx, fy, fz) = (spectral(&|kx, _| (0.0, kx)), spectral(&|_, ky| (0.0, ky)), spectral(&|kx, ky| (kx.hypot(ky), 0.0)));
        let kept = euler_windows(&grid, &field, &fx, &fy, &fz, 3.0, 12, 4, 0.1);
        let best = kept[0];
        assert!((best[2] / z0 - 1.0).abs() < 0.02 && (best[0] - x0).hypot(best[1] - y0) < 100.0, "{best:?}");
        assert!(kept.len() >= 10 && kept.iter().all(|p| p[2] > 0.0));
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
