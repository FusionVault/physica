//! A minimal closed-loop simulation of a GRACE-type satellite pair.
//!
//! Two satellites fly the same circular orbit a fixed distance apart over a rotating Earth. The
//! observable is the line-of-sight difference of their gravitational accelerations,
//! `e_AB · (∇V(r_B) − ∇V(r_A))` (the "acceleration approach"; it is what Zingerle et al. 2024 call
//! product noise is measured against). Noise with a chosen ASD is added, the series is whitened by
//! the same ASD, and the spherical-harmonic coefficients from degree 2 are estimated by least
//! squares. The errors are the estimate minus the truth.
//!
//! Not modelled: orbit perturbations and dynamic orbit determination, non-gravitational forces,
//! time-variable gravity and therefore temporal aliasing, attitude errors, regularization.

use geopotential::{Coefficients, Legendre};

use crate::noise::{colored_noise, whiten};

/// The pair's orbit.
#[derive(Clone, Copy, Debug)]
pub struct Pair {
    /// Altitude above the reference radius, m.
    pub altitude: f64,
    /// Inclination, radians.
    pub inclination: f64,
    /// Distance between the satellites along the orbit, m.
    pub separation: f64,
}

/// Earth rotation rate, rad/s (IERS).
pub const EARTH_ROTATION: f64 = 7.292_115e-5;

impl Pair {
    /// Earth-fixed positions of the leading (A) and trailing (B) satellite `t` seconds after
    /// A crosses the equator northbound at longitude 0.
    pub fn positions(&self, gm: f64, radius: f64, t: f64) -> [[f64; 3]; 2] {
        let a = radius + self.altitude;
        let n = (gm / a.powi(3)).sqrt();
        let theta = EARTH_ROTATION * t;
        let (si, ci) = self.inclination.sin_cos();
        let at = |u: f64| {
            let (su, cu) = u.sin_cos();
            let (x, y, z) = (a * cu, a * su * ci, a * su * si);
            let (st, ct) = theta.sin_cos();
            [x * ct + y * st, -x * st + y * ct, z]
        };
        [at(n * t), at(n * t - self.separation / a)]
    }
}

/// Unknowns: `(n, m, is_sine)` for `n = 2..=nmax`, cosine then sine.
fn unknowns(nmax: usize) -> Vec<(usize, usize, bool)> {
    (2..=nmax).flat_map(|n| (0..=n).flat_map(move |m| [(n, m, false), (n, m, true)].into_iter().filter(move |u| !(u.1 == 0 && u.2)))).collect()
}

/// Gravitational acceleration (ECEF) of each unknown's unit coefficient at `x`, accumulated into
/// `out` with weight `w`: `out[k] += w · e·∇Y_k(x)`.
#[allow(clippy::too_many_arguments)]
fn add_los_partials(leg: &mut Legendre, gm: f64, radius: f64, x: [f64; 3], e: [f64; 3], w: f64, cols: &[(usize, usize, bool)], out: &mut [f64]) {
    let r = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
    let (lat, lon) = ((x[2] / r).asin(), x[1].atan2(x[0]));
    leg.compute_with_derivative(lat);
    let ((sp, cp), (sl, cl)) = (lat.sin_cos(), lon.sin_cos());
    // e projected on the local north, east and up directions.
    let e_n = -sp * cl * e[0] - sp * sl * e[1] + cp * e[2];
    let e_e = -sl * e[0] + cl * e[1];
    let e_u = cp * cl * e[0] + cp * sl * e[1] + sp * e[2];
    let q = radius / r;
    let mut qn = vec![0.0; leg.max_degree() + 1];
    qn.iter_mut().enumerate().for_each(|(n, v)| *v = gm / r * q.powi(n as i32));
    for (k, &(n, m, sine)) in cols.iter().enumerate() {
        let (s, c) = (m as f64 * lon).sin_cos();
        let (trig, dtrig) = if sine { (s, m as f64 * c) } else { (c, -(m as f64) * s) };
        let v = qn[n];
        let up = -((n + 1) as f64) / r * v * leg.p(n, m) * trig;
        let north = v / r * leg.dp(n, m) * trig;
        let east = if m == 0 { 0.0 } else { v * leg.p(n, m) * dtrig / (r * cp) };
        out[k] += w * (e_n * north + e_e * east + e_u * up);
    }
}

/// Result of one simulated solution.
pub struct Solution {
    /// Estimate minus truth (degrees 0 and 1 are zero).
    pub errors: Coefficients,
    /// Formal standard deviations from the whitened normal equations.
    pub formal: Coefficients,
}

/// Simulate `duration` seconds of observations every `dt` seconds with noise of ASD `asd`
/// (m s⁻² Hz⁻¹ᐟ²), add the signal of `truth` (if any), and estimate degrees 2..=`nmax`.
///
/// Memory is `unknowns × samples` f64 values (about 1.2 GB for degree 40, one month at 30 s).
pub fn simulate(
    pair: &Pair,
    gm: f64,
    radius: f64,
    nmax: usize,
    duration: f64,
    dt: f64,
    asd: impl Fn(f64) -> f64 + Copy,
    truth: Option<&Coefficients>,
    seed: u64,
) -> Solution {
    let cols = unknowns(nmax);
    let (k, samples) = (cols.len(), (duration / dt) as usize);
    let mut a = vec![0.0; k * samples]; // column-major
    let mut y = colored_noise(asd, dt, samples, seed);
    let mut leg = Legendre::new(nmax);
    let mut row = vec![0.0; k];
    for t in 0..samples {
        let [ra, rb] = pair.positions(gm, radius, t as f64 * dt);
        let d = [rb[0] - ra[0], rb[1] - ra[1], rb[2] - ra[2]];
        let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        let e = [d[0] / len, d[1] / len, d[2] / len];
        row.iter_mut().for_each(|v| *v = 0.0);
        add_los_partials(&mut leg, gm, radius, rb, e, 1.0, &cols, &mut row);
        add_los_partials(&mut leg, gm, radius, ra, e, -1.0, &cols, &mut row);
        if let Some(tr) = truth {
            y[t] += cols.iter().zip(&row).map(|(&(n, m, s), r)| r * if n > tr.max_degree() { 0.0 } else if s { tr.s(n, m) } else { tr.c(n, m) }).sum::<f64>();
        }
        for (j, v) in row.iter().enumerate() {
            a[j * samples + t] = *v;
        }
    }
    // Whiten observations and every column by the noise ASD.
    whiten(&mut y, dt, asd);
    let mut aty = vec![0.0; k];
    for j in 0..k {
        let col = &mut a[j * samples..(j + 1) * samples];
        whiten(col, dt, asd);
        aty[j] = col.iter().zip(&y).map(|(c, v)| c * v).sum();
    }
    // Normal equations N = AᵀA (lower triangle), Cholesky N = LLᵀ.
    let mut nm = vec![0.0; k * k];
    for i in 0..k {
        let ci = &a[i * samples..(i + 1) * samples];
        for j in 0..=i {
            let cj = &a[j * samples..(j + 1) * samples];
            nm[i * k + j] = ci.iter().zip(cj).map(|(p, q)| p * q).sum();
        }
    }
    drop(a);
    cholesky(&mut nm, k);
    let x = cholesky_solve(&nm, k, &aty);
    let inv_diag = cholesky_inverse_diagonal(&nm, k);
    // Whitened noise has per-sample variance 1/(2·dt); scale the formal errors accordingly.
    let s0 = (1.0 / (2.0 * dt)).sqrt();
    let mut errors = Coefficients::zeros(gm, radius, nmax);
    let mut formal = Coefficients::zeros(gm, radius, nmax);
    for (i, &(n, m, sine)) in cols.iter().enumerate() {
        let tv = truth.filter(|tr| n <= tr.max_degree()).map_or(0.0, |tr| if sine { tr.s(n, m) } else { tr.c(n, m) });
        let (e, f) = (x[i] - tv, s0 * inv_diag[i].sqrt());
        let (ec, es, fc, fs) = (errors.c(n, m), errors.s(n, m), formal.c(n, m), formal.s(n, m));
        if sine {
            errors.set(n, m, ec, e);
            formal.set(n, m, fc, f);
        } else {
            errors.set(n, m, e, es);
            formal.set(n, m, f, fs);
        }
    }
    Solution { errors, formal }
}

/// In-place Cholesky of a symmetric positive-definite matrix stored in the lower triangle.
fn cholesky(a: &mut [f64], n: usize) {
    for j in 0..n {
        let mut d = a[j * n + j];
        for p in 0..j {
            d -= a[j * n + p] * a[j * n + p];
        }
        assert!(d > 0.0, "normal matrix not positive definite at {j}");
        let d = d.sqrt();
        a[j * n + j] = d;
        for i in j + 1..n {
            let mut s = a[i * n + j];
            for p in 0..j {
                s -= a[i * n + p] * a[j * n + p];
            }
            a[i * n + j] = s / d;
        }
    }
}

fn cholesky_solve(l: &[f64], n: usize, b: &[f64]) -> Vec<f64> {
    let mut z = b.to_vec();
    for i in 0..n {
        for p in 0..i {
            z[i] -= l[i * n + p] * z[p];
        }
        z[i] /= l[i * n + i];
    }
    for i in (0..n).rev() {
        for p in i + 1..n {
            z[i] -= l[p * n + i] * z[p];
        }
        z[i] /= l[i * n + i];
    }
    z
}

/// Diagonal of N⁻¹ = L⁻ᵀL⁻¹: the squared column norms of L⁻¹.
fn cholesky_inverse_diagonal(l: &[f64], n: usize) -> Vec<f64> {
    let mut diag = vec![0.0; n];
    let mut col = vec![0.0; n];
    for j in 0..n {
        // Solve L·c = e_j.
        col.iter_mut().for_each(|v| *v = 0.0);
        col[j] = 1.0 / l[j * n + j];
        for i in j + 1..n {
            let mut s = 0.0;
            for p in j..i {
                s -= l[i * n + p] * col[p];
            }
            col[i] = s / l[i * n + i];
        }
        for (i, c) in col.iter().enumerate().skip(j) {
            diag[i] += c * c;
        }
    }
    diag
}

#[cfg(test)]
mod tests {
    use super::*;

    const GM: f64 = 3.986004415e14;
    const R: f64 = 6378136.3;
    const GRACE: Pair = Pair { altitude: 463e3, inclination: 89f64.to_radians(), separation: 200e3 };

    #[test]
    fn orbit_geometry() {
        let [a, b] = GRACE.positions(GM, R, 1234.0);
        let r = |p: [f64; 3]| (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        assert!((r(a) - (R + 463e3)).abs() < 1e-6 && (r(b) - (R + 463e3)).abs() < 1e-6);
        let chord = r([a[0] - b[0], a[1] - b[1], a[2] - b[2]]);
        assert!((chord - 200e3).abs() < 30.0, "{chord}"); // the chord is ~7 m shorter than the 200 km arc
    }

    #[test]
    fn los_partials_match_finite_differences_of_the_field() {
        // One coefficient's LOS acceleration difference against geopotential's own gradient.
        let nmax = 6;
        let cols = unknowns(nmax);
        let mut leg = Legendre::new(nmax);
        let [ra, rb] = GRACE.positions(GM, R, 777.0);
        let d = [rb[0] - ra[0], rb[1] - ra[1], rb[2] - ra[2]];
        let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        let e = [d[0] / len, d[1] / len, d[2] / len];
        let mut row = vec![0.0; cols.len()];
        add_los_partials(&mut leg, GM, R, rb, e, 1.0, &cols, &mut row);
        add_los_partials(&mut leg, GM, R, ra, e, -1.0, &cols, &mut row);
        for (k, &(n, m, sine)) in cols.iter().enumerate().step_by(5) {
            let mut c = Coefficients::zeros(GM, R, nmax);
            c.set(n, m, if sine { 0.0 } else { 1.0 }, if sine { 1.0 } else { 0.0 });
            let v = |p: [f64; 3]| {
                let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
                c.partials(r, (p[2] / r).asin(), p[1].atan2(p[0])).v
            };
            let h = 0.5;
            let grad_e = |p: [f64; 3]| (v([p[0] + h * e[0], p[1] + h * e[1], p[2] + h * e[2]]) - v([p[0] - h * e[0], p[1] - h * e[1], p[2] - h * e[2]])) / (2.0 * h);
            let fd = grad_e(rb) - grad_e(ra);
            assert!((row[k] - fd).abs() < 1e-6 * row[k].abs().max(1e-3), "({n},{m},{sine}): {} vs {fd}", row[k]);
        }
    }

    #[test]
    fn noise_free_closed_loop_recovers_the_truth() {
        let mut truth = Coefficients::zeros(GM, R, 8);
        let mut s = 1u64;
        for n in 2..=8 {
            for m in 0..=n {
                s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                let v = 1e-6 / (n * n) as f64 * (((s >> 33) as f64 / (1u64 << 31) as f64) - 1.0);
                truth.set(n, m, v, if m == 0 { 0.0 } else { -v });
            }
        }
        // Noise far below the signal: the estimate must reproduce the truth closely.
        let sol = simulate(&GRACE, GM, R, 8, 2.0 * 86_400.0, 30.0, |_| 1e-16, Some(&truth), 1);
        let worst = (2..=8).flat_map(|n| (0..=n).map(move |m| (n, m))).map(|(n, m)| sol.errors.c(n, m).abs().max(sol.errors.s(n, m).abs())).fold(0.0, f64::max);
        assert!(worst < 1e-12, "worst error {worst:e}");
    }

    #[test]
    fn empirical_errors_agree_with_formal_errors() {
        // White noise: the RMS of (estimate − truth)/formal over all coefficients is ~1.
        let sol = simulate(&GRACE, GM, R, 10, 4.0 * 86_400.0, 30.0, |_| 1e-10, None, 42);
        let mut sum = 0.0;
        let mut count = 0.0;
        for n in 2..=10 {
            for m in 0..=n {
                sum += (sol.errors.c(n, m) / sol.formal.c(n, m)).powi(2);
                count += 1.0;
                if m > 0 {
                    sum += (sol.errors.s(n, m) / sol.formal.s(n, m)).powi(2);
                    count += 1.0;
                }
            }
        }
        let ratio = (sum / count).sqrt();
        assert!((ratio - 1.0).abs() < 0.15, "RMS(error/formal) = {ratio}");
    }
}
