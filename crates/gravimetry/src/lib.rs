#![doc = include_str!("../README.md")]

use geopotential::{Coefficients, Legendre};

pub mod atmosphere;
pub mod atom;
pub mod coseismic;
pub mod gradiometry;
pub mod kalman;
pub mod local;
pub mod missions;
pub mod noise;
pub mod sim;
pub mod technotes;

/// Area of the world ocean, m²: 3.618×10⁸ km², so 361.8 Gt of water raise sea level 1 mm.
pub const OCEAN_AREA: f64 = 3.618e14;

/// Global mean sea-level rise (m) from `mass` (kg) of water added to the ocean and spread evenly:
/// `M/(ρ_w A)`. Ignores self-attraction and loading (sea-level fingerprints).
pub fn sea_level_equivalent(mass: f64) -> f64 {
    mass / (geopotential::RHO_WATER * OCEAN_AREA)
}

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
/// Jekeli (1981) / Wahr et al. (1998, eqs. 32–34), normalized so `W_0 = 1`.
///
/// The weights are ratios of modified spherical Bessel functions, `W_n = i_n(b)/i_0(b)` with
/// `b = ln 2/(1 − cos(r/a))`: the decaying solution of their three-term recurrence. Run forward,
/// as Wahr et al. print it, that recurrence loses every digit once the weights are small (at
/// 1000 km it turns negative by degree 60). So the ratios `ρ_n = W_n/W_{n−1} =
/// 1/((2n+1)/b + ρ_{n+1})` are run backward instead (Gautschi, SIAM Rev. 9, 24, 1967), from a
/// depth where the starting guess no longer matters: an error there shrinks like
/// `exp(−(N² − n²)/b)`, so `N² ≥ nmax² + 40b` leaves it below 10⁻¹⁷. Every weight is then exact
/// to about 10⁻¹⁴, however small: W₉₆ at 1000 km is 1.39×10⁻³¹. For radii so small that
/// `b > 10⁶(nmax+1)²` the weights are `exp(−n(n+1)/2b)` to 10⁻¹⁴, which avoids millions of steps.
pub fn gaussian_weights(radius: f64, earth_radius: f64, nmax: usize) -> Vec<f64> {
    let b = core::f64::consts::LN_2 / (1.0 - (radius / earth_radius).cos());
    let top = ((nmax + 1) * (nmax + 1)) as f64;
    if b >= 1e6 * top {
        // Also radius 0, where b is infinite and every weight is 1.
        return (0..=nmax).map(|n| (-((n * (n + 1)) as f64) / (2.0 * b)).exp()).collect();
    }
    let start = ((nmax * nmax) as f64 + 40.0 * b).sqrt().ceil() as usize + 10;
    // Amos's estimate of i_N/i_{N−1} (Math. Comp. 28, 239, 1974) as the starting guess.
    let mut rho = b / (start as f64 + ((start * start) as f64 + b * b).sqrt());
    let mut ratios = vec![0.0; nmax + 1];
    for n in (1..=start).rev() {
        rho = 1.0 / ((2 * n + 1) as f64 / b + rho);
        if n <= nmax {
            ratios[n] = rho;
        }
    }
    let mut w = Vec::with_capacity(nmax + 1);
    w.push(1.0);
    for n in 1..=nmax {
        w.push(w[n - 1] * ratios[n]);
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

/// Spherical-harmonic coefficients, to degree `nmax`, of point masses `(lat, lon, kg)` (radians)
/// on the sphere of radius `radius`, in the convention of a field with `gm`:
/// `C̄nm + iS̄nm = (m_k / M) · P̄nm(sin φ_k) · e^{imλ_k} / (2n + 1)` summed over the masses,
/// with `M = gm / G`. Degree 0 holds the net mass; it is kept so the caller can check it.
///
/// ```
/// use gravimetry::mass_coefficients;
///
/// // One gigatonne at the North Pole.
/// let pole = (90f64.to_radians(), 0.0, 1e12);
/// let field = mass_coefficients(&[pole], 3.986004415e14, 6_378_136.3, 10);
///
/// // A mass on the axis changes only the zonal (m = 0) terms.
/// assert!(field.c(2, 0) > 0.0);
/// assert!(field.c(2, 1).abs() < 1e-25);
/// ```
pub fn mass_coefficients(points: &[(f64, f64, f64)], gm: f64, radius: f64, nmax: usize) -> Coefficients {
    let earth_mass = gm / metrology::G;
    let mut out = Coefficients::zeros(gm, radius, nmax);
    let mut leg = Legendre::new(nmax);
    let (mut c, mut s) = (vec![0.0; (nmax + 1) * (nmax + 2) / 2], vec![0.0; (nmax + 1) * (nmax + 2) / 2]);
    for &(lat, lon, kg) in points {
        leg.compute(lat);
        let k = kg / earth_mass;
        for m in 0..=nmax {
            let (sm, cm) = (m as f64 * lon).sin_cos();
            for n in m..=nmax {
                let v = k * leg.p(n, m) / (2 * n + 1) as f64;
                let i = n * (n + 1) / 2 + m;
                c[i] += v * cm;
                s[i] += v * sm;
            }
        }
    }
    for n in 0..=nmax {
        for m in 0..=n {
            let i = n * (n + 1) / 2 + m;
            out.set(n, m, c[i], s[i]);
        }
    }
    out
}

/// Degree weights `R²/(n(n+1))` (zero at degree 0) that turn the rate of change of a field into
/// its minimum-norm mass-transport potential φ, the solution of `∇²φ = −dσ/dt` on the sphere.
/// Apply them with [`smooth`] to a field of yearly change, synthesize
/// [`geopotential::Functional::WaterHeight`] and multiply by the water density: that is φ in
/// kg/yr, and its horizontal gradient is the smallest mass flux (kg per metre per year) that
/// could carry the change. It points toward where mass is gathering.
pub fn transport_weights(radius: f64, nmax: usize) -> Vec<f64> {
    (0..=nmax).map(|n| if n == 0 { 0.0 } else { radius * radius / (n * (n + 1)) as f64 }).collect()
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
        // Against the Legendre transform of the Gaussian kernel, by quadrature (and, with
        // a = 6 378 136.3 m, the value the cross-domain report quotes: 0.475 974 043 972 641).
        assert!((w[30] / 0.475_182_474_894_050 - 1.0).abs() < 1e-12, "{}", w[30]);
        let report = gaussian_weights(300e3, 6_378_136.3, 30);
        assert!((report[30] / 0.475_974_043_972_641 - 1.0).abs() < 1e-12);
        // A 1000 km filter decays without a floor, positive and monotone, to the values of the
        // positive-term series of i_n(b)/i_0(b): W₉₆ = 1.391 224 300 882 06×10⁻³¹,
        // W₁₂₀ = 2.702 118 043 068 33×10⁻⁴⁶.
        let wide = gaussian_weights(1000e3, r, 120);
        assert!(wide.windows(2).all(|p| p[1] < p[0] && p[1] > 0.0));
        assert!((wide[96] / 1.391_224_300_882_055e-31 - 1.0).abs() < 1e-12, "{}", wide[96]);
        assert!((wide[120] / 2.702_118_043_068_331e-46 - 1.0).abs() < 1e-12, "{}", wide[120]);
        // A 20 km filter needs a deep start (2 383 steps): a start at nmax + 400 is 2×10⁻⁸ off.
        let narrow = gaussian_weights(20e3, r, 60);
        assert!((narrow[60] / 0.987_075_319_875_3 - 1.0).abs() < 1e-12, "{}", narrow[60]);
        // No smoothing at radius 0, and the large-b form joins the recurrence smoothly.
        assert!(gaussian_weights(0.0, r, 10).iter().all(|&x| x == 1.0));
        // For degree 60 the switch at b = 10⁶·61² falls at a 123 m radius. Just above it, where the
        // recurrence still runs (123.1 m, 2.4×10⁵ steps), it agrees with exp(−n(n+1)/2b).
        let rad = 123.1;
        let b = core::f64::consts::LN_2 / (1.0 - (rad / r).cos());
        assert!(b < 1e6 * 61.0 * 61.0);
        let w = gaussian_weights(rad, r, 60);
        for n in [1, 30, 60] {
            let asymptotic = (-((n * (n + 1)) as f64) / (2.0 * b)).exp();
            assert!((w[n] / asymptotic - 1.0).abs() < 1e-13, "n={n}: {} vs {asymptotic}", w[n]);
        }
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

    #[test]
    fn point_masses_match_the_addition_theorem() {
        // A mass's potential from its coefficients equals G·m / distance outside the sphere.
        let (gm, r) = (3.986004415e14, 6378136.3);
        let (lat0, lon0, kg) = (0.4, 1.1, 1e13);
        let c = mass_coefficients(&[(lat0, lon0, kg)], gm, r, 200);
        let (rr, lat, lon) = (r + 400e3, 0.45, 1.13);
        let v = c.partials(rr, lat, lon).v;
        let p = |la: f64, lo: f64, rad: f64| [rad * la.cos() * lo.cos(), rad * la.cos() * lo.sin(), rad * la.sin()];
        let (a, b) = (p(lat0, lon0, r), p(lat, lon, rr));
        let dist = ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt();
        let exact = metrology::G * kg / dist;
        assert!((v - exact).abs() < 1e-3 * exact, "{v:e} vs {exact:e}");
        assert!((c.c(0, 0) / (kg * metrology::G / gm) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn transport_weights_invert_the_laplacian() {
        let w = transport_weights(2.0, 3);
        assert_eq!(w, vec![0.0, 2.0, 4.0 / 6.0, 4.0 / 12.0]);
    }
}
