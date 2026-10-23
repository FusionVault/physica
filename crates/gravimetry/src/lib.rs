//! GRACE-style processing of monthly gravity fields into water-storage anomalies.
//!
//! The chain is: replace `C̄20`/`C̄30` and add degree 1 (from `gravity-field`'s TN-14/TN-13
//! readers), [`anomaly`] against a baseline [`mean`], [`smooth`] with [`gaussian_weights`], then
//! synthesize [`geopotential::Functional::WaterHeight`] with [`love_numbers`].
//!
//! ```
//! use geopotential::{Coefficients, Functional};
//!
//! let (gm, r) = (3.986004415e14, 6378136.3);
//! let mut month = Coefficients::zeros(gm, r, 60);
//! month.set(2, 0, 1e-10, 0.0);
//! let baseline = gravimetry::mean(&[Coefficients::zeros(gm, r, 60)]).unwrap();
//! let w = gravimetry::gaussian_weights(300e3, r, 60);
//! let dewh = gravimetry::smooth(&gravimetry::anomaly(&month, &baseline), &w);
//! let love = gravimetry::love_numbers(60).unwrap();
//! let pole = dewh.grid(&Functional::WaterHeight { love: &love }, r, &[1.5707963], &[0.0])[0];
//! assert!(pole > 0.01 && pole < 0.02); // ~19 mm of water, slightly smoothed
//! ```

use geopotential::Coefficients;

/// Elastic load Love numbers `k'_n` (PREM, Han & Wahr 1995) as tabulated by Wahr et al. (1998,
/// JGR 103, Table 1), linearly interpolated between tabulated degrees (error < 0.05%), for
/// `n = 0..=nmax`. `k'_1 = 0.027` is the centre-of-figure value. `None` beyond degree 200, the end
/// of the table.
pub fn love_numbers(nmax: usize) -> Option<Vec<f64>> {
    const TABLE: [(usize, f64); 21] = [
        (0, 0.0), (1, 0.027), (2, -0.303), (3, -0.194), (4, -0.132), (5, -0.104), (6, -0.089),
        (7, -0.081), (8, -0.076), (9, -0.072), (10, -0.069), (12, -0.064), (15, -0.058),
        (20, -0.051), (30, -0.040), (40, -0.033), (50, -0.027), (70, -0.020), (100, -0.014),
        (150, -0.010), (200, -0.007),
    ];
    if nmax > 200 {
        return None;
    }
    let k = (0..=nmax).map(|n| {
        let hi = TABLE.iter().position(|&(d, _)| d >= n).expect("n ≤ 200");
        let (d1, k1) = TABLE[hi];
        if d1 == n {
            return k1;
        }
        let (d0, k0) = TABLE[hi - 1];
        k0 + (k1 - k0) * (n - d0) as f64 / (d1 - d0) as f64
    });
    Some(k.collect())
}

/// Gaussian smoothing weights for averaging radius `radius` (m; the distance at which the
/// spatial kernel falls to half), on a sphere of radius `earth_radius`, for `n = 0..=nmax`.
/// Jekeli (1981) / Wahr et al. (1998, eqs. 32–34), normalized so `W_0 = 1`. The recursion is
/// unstable once the weights are tiny, so they are held at 1e-10 from there on, as gravity-toolkit
/// does.
pub fn gaussian_weights(radius: f64, earth_radius: f64, nmax: usize) -> Vec<f64> {
    let b = core::f64::consts::LN_2 / (1.0 - (radius / earth_radius).cos());
    let e = (-2.0 * b).exp();
    let mut w = vec![1.0];
    if nmax >= 1 {
        w.push((1.0 + e) / (1.0 - e) - 1.0 / b);
    }
    for n in 1..nmax {
        let next = -((2 * n + 1) as f64) / b * w[n] + w[n - 1];
        let floor = w[n] <= 1e-10 || !(1e-10..=w[n]).contains(&next);
        w.push(if floor { 1e-10 } else { next });
    }
    w
}

/// Multiply every degree `n` by `weights[n]` (e.g. [`gaussian_weights`]).
///
/// # Panics
/// If `weights` is shorter than `max_degree() + 1`.
pub fn smooth(field: &Coefficients, weights: &[f64]) -> Coefficients {
    let mut out = field.clone();
    for (n, &w) in weights[..=field.max_degree()].iter().enumerate() {
        for m in 0..=n {
            out.set(n, m, field.c(n, m) * w, field.s(n, m) * w);
        }
    }
    out
}

/// The mean of several fields, in the first field's `GM` and radius (others are rescaled), to the
/// lowest degree among them. `None` if `fields` is empty.
pub fn mean(fields: &[Coefficients]) -> Option<Coefficients> {
    let first = fields.first()?;
    let nmax = fields.iter().map(Coefficients::max_degree).min()?;
    let mut sum = Coefficients::zeros(first.gm, first.radius, nmax);
    sum.tide_system = first.tide_system;
    let k = 1.0 / fields.len() as f64;
    for f in fields {
        let f = f.rescaled(first.gm, first.radius);
        for n in 0..=nmax {
            for m in 0..=n {
                sum.set(n, m, sum.c(n, m) + k * f.c(n, m), sum.s(n, m) + k * f.s(n, m));
            }
        }
    }
    Some(sum)
}

/// `field − reference`, with the reference rescaled to the field's `GM` and radius, to the lower
/// of the two degrees.
pub fn anomaly(field: &Coefficients, reference: &Coefficients) -> Coefficients {
    let r = reference.rescaled(field.gm, field.radius);
    let nmax = field.max_degree().min(r.max_degree());
    let mut out = Coefficients::zeros(field.gm, field.radius, nmax);
    out.tide_system = field.tide_system;
    for n in 0..=nmax {
        for m in 0..=n {
            out.set(n, m, field.c(n, m) - r.c(n, m), field.s(n, m) - r.s(n, m));
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn love_numbers_hit_the_table_and_interpolate() {
        let k = love_numbers(200).unwrap();
        assert_eq!((k[0], k[1], k[2], k[50], k[200]), (0.0, 0.027, -0.303, -0.027, -0.007));
        assert!((k[11] - -0.0665).abs() < 1e-12); // halfway between n = 10 and 12
        assert!((k[60] - -0.0235).abs() < 1e-12);
        assert!(love_numbers(201).is_none());
    }

    #[test]
    fn gaussian_weights_behave() {
        let r = 6_371_000.0;
        let w = gaussian_weights(300e3, r, 120);
        assert_eq!(w[0], 1.0);
        assert!(w.windows(2).all(|p| p[1] <= p[0]), "monotone");
        // Small-angle limit: W_n ≈ exp(−n(n+1)/(2b)), b = ln2/(1 − cos(r/a)) (Wahr et al. 1998).
        let b = core::f64::consts::LN_2 / (1.0 - (300e3 / r).cos());
        for n in [1, 5, 10, 20, 40] {
            let approx = (-((n * (n + 1)) as f64) / (2.0 * b)).exp();
            assert!((w[n] - approx).abs() < 0.01, "n={n}: {} vs {approx}", w[n]);
        }
        // Half weight where n(n+1) = 2b·ln2, degree ≈ 29 for 300 km.
        let half = w.iter().position(|&x| x < 0.5).unwrap();
        assert!((28..=30).contains(&half), "{half}");
        // A 1000 km filter decays to the floor without going negative or oscillating.
        let wide = gaussian_weights(1000e3, r, 120);
        assert!(wide.iter().all(|&x| x >= 1e-10) && wide[120] == 1e-10);
    }

    #[test]
    fn mean_and_anomaly_rescale_and_truncate() {
        let mut a = Coefficients::zeros(4e14, 6.4e6, 3);
        a.set(2, 1, 2.0e-9, -1.0e-9);
        let b = a.rescaled(4.1e14, 6.5e6); // same field, other convention, as in another file
        let m = mean(&[a.clone(), b]).unwrap();
        assert!((m.c(2, 1) - 2.0e-9).abs() < 1e-24 && (m.s(2, 1) + 1.0e-9).abs() < 1e-24);
        let d = anomaly(&a, &Coefficients::zeros(4e14, 6.4e6, 2));
        assert_eq!((d.max_degree(), d.c(2, 1)), (2, 2.0e-9));
        assert!(mean(&[]).is_none());
    }
}
