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

/// Remove north–south stripes (correlated errors) the way Swenson & Wahr (2006, GRL 33, L08402)
/// do: for each order `m ≥ 5`, fit a moving quadratic across degree, separately for even and odd
/// degrees, over `2nₛ+1` coefficients with `nₛ = max(round(15·e^(−m/10)), 2)`, and subtract the
/// fit. Sequences too short for one window are removed entirely. Apply it to an anomaly
/// (typically truncated to degree 60), before Gaussian smoothing.
///
/// Ported from gravity-toolkit's `destripe_harmonics` (MIT, Copyright (c) 2020 Tyler C
/// Sutterley; after I. Velicogna's `remove_errors.f`), default settings.
pub fn destripe(field: &Coefficients) -> Coefficients {
    let nmax = field.max_degree();
    let mut out = field.clone();
    for m in 5..=nmax {
        let ns = ((15.0 * (-(m as f64) / 10.0).exp()).round() as usize).max(2);
        // Window sums of ll⁰, ll², ll⁴ for ll = −nₛ..=nₛ (odd powers vanish by symmetry).
        let (s0, s2, s4) = (-(ns as i64)..=ns as i64).fold((0.0, 0.0, 0.0), |(a, b, c), ll| {
            let l2 = (ll * ll) as f64;
            (a + 1.0, b + l2, c + l2 * l2)
        });
        let det = s0 * s4 - s2 * s2;
        for parity in 0..2 {
            let degrees: Vec<usize> = (m..=nmax).filter(|l| l % 2 == parity).collect();
            for part in 0..2 {
                let x: Vec<f64> = degrees.iter().map(|&l| if part == 0 { field.c(l, m) } else { field.s(l, m) }).collect();
                let fit = smooth_part(&x, ns, s2, (s0, s4, det));
                for (&l, f) in degrees.iter().zip(fit) {
                    let (c, s) = (out.c(l, m), out.s(l, m));
                    if part == 0 { out.set(l, m, c - f, s) } else { out.set(l, m, c, s - f) }
                }
            }
        }
    }
    out
}

/// The moving-quadratic smooth of `x` (see [`destripe`]); `x` itself when it is shorter than one
/// window.
fn smooth_part(x: &[f64], ns: usize, s2: f64, (s0, s4, det): (f64, f64, f64)) -> Vec<f64> {
    if x.len() < 2 * ns + 1 {
        return x.to_vec();
    }
    let mut out = vec![0.0; x.len()];
    let last = x.len() - 1 - ns;
    for l in ns..=last {
        let (mut r0, mut r1, mut r2) = (0.0, 0.0, 0.0);
        for ll in -(ns as i64)..=ns as i64 {
            let v = x[(l as i64 + ll) as usize];
            let llf = ll as f64;
            (r0, r1, r2) = (r0 + v, r1 + v * llf, r2 + v * llf * llf);
        }
        // Normal equations of a + b·ll + c·ll² on a symmetric window.
        let (a, b, c) = ((s4 * r0 - s2 * r2) / det, r1 / s2, (s0 * r2 - s2 * r0) / det);
        let poly = |ll: f64| a + b * ll + c * ll * ll;
        out[l] = a;
        // The first and last windows also cover the ends of the sequence.
        if l == ns {
            (0..ns).for_each(|k| out[k] = poly(k as f64 - ns as f64));
        }
        if l == last {
            (1..=ns).for_each(|k| out[l + k] = poly(k as f64));
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
    fn destripe_removes_smooth_patterns_and_keeps_spikes() {
        let nmax = 60;
        let mut f = Coefficients::zeros(4e14, 6.4e6, nmax);
        for n in 0..=nmax {
            for m in 0..=n {
                let q = 1e-10 * (1.0 + 0.1 * n as f64 + 0.01 * (n * n) as f64);
                f.set(n, m, q, if m == 0 { 0.0 } else { -q });
            }
        }
        let d = destripe(&f);
        let worst = (0..=nmax).flat_map(|n| (5..=n).map(move |m| (n, m)))
            .map(|(n, m)| d.c(n, m).abs().max(d.s(n, m).abs()))
            .fold(0.0, f64::max);
        assert!(worst < 1e-20, "quadratic left {worst:e}");
        assert_eq!((d.c(40, 4), d.s(40, 4)), (f.c(40, 4), f.s(40, 4)), "orders < 5 untouched");
        // A lone coefficient is mostly signal to the filter: at m = 5 (19-point window) 88% stays.
        let mut spike = Coefficients::zeros(4e14, 6.4e6, nmax);
        spike.set(30, 5, 1e-10, 0.0);
        let kept = destripe(&spike).c(30, 5) / 1e-10;
        assert!((kept - (1.0 - 30666.0 / 257754.0)).abs() < 1e-12, "{kept}");
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
