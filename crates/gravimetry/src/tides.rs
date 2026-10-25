//! The Earth tide at a gravimeter: Longman's lunar and solar tide, the Love-number factors that
//! scale a rigid Earth's tide to the observed one, and the permanent tide that separates the
//! mean-tide, zero-tide and tide-free conventions.
//!
//! ```
//! use geopotential::julian_epoch;
//! use gravimetry::tides::{delta, longman, permanent};
//!
//! // State College, Pennsylvania, at 2015-04-23 00:00 UT: 3.5 µGal of lunar and solar tide
//! // (LongmanTide's test vector), with Longman's own Love numbers.
//! let (lat, lon) = (40.7914f64.to_radians(), 282.1414f64.to_radians());
//! let (moon, sun, total) = longman(lat, lon, 370.0, julian_epoch(2015, 4, 23, 0, 0));
//! assert!((total / 1e-8 - 3.5347).abs() < 1e-4 && (moon + sun - total).abs() < 1e-18);
//! assert!((delta(2, 0.612, 0.303) - 1.1575).abs() < 1e-15);
//! // The permanent tide lowers gravity 30.49 µGal at the equator (Mäkinen 2021, Eq. 20).
//! assert!((permanent(0.0, 0.0).1 / 1e-8 + 30.49).abs() < 0.01);
//! ```

use core::f64::consts::TAU;

use geopotential::geocentric;
use metrology::WGS84_A;

/// The time average of the tide-generating potential in Mäkinen's (2021, Eq. 18) form
/// `W_T = A(r/a)²(sin²ψ − 1/3)`, m²/s²: −2.9166, from the IERS Conventions 2010, which Mäkinen
/// proposes for the International Height Reference Frame. Spectral expansions from Doodson's to
/// Kudryavtsev's agree to 0.0002 (his Table 1).
pub const PERMANENT_A: f64 = -2.9166;

/// The gravimetric factor `δₙ = 1 + (2/n)hₙ − ((n+1)/n)kₙ` of degree `n` for Love numbers `h`
/// and `k`: the observed tide in gravity over a rigid Earth's. It sums the direct pull `nW/r` of
/// a tide potential `W ∝ rⁿ`, the free-air change `2hW/r` of the surface lifted `hW/g`, and the
/// pull `−(n + 1)kW/r` of the Earth's own redistributed mass. With the IERS Conventions' nominal
/// `h₂ = 0.6078` and `k₂₀ = 0.30190` it is 1.155, the conventional 1.16.
pub fn delta(n: u32, h: f64, k: f64) -> f64 {
    let n = f64::from(n);
    1.0 + 2.0 * h / n - (n + 1.0) * k / n
}

/// The diminishing factor `γ = 1 + k − h` for Love numbers `h` and `k` of one degree: the
/// observed tide in tilt, or in a potential difference measured between points that rise with
/// the ground, over a rigid Earth's. About 0.69 at degree 2 (Mäkinen 2021, Eq. 52).
pub fn gamma(h: f64, k: f64) -> f64 {
    1.0 + k - h
}

/// The lunar and solar tide in gravity at geodetic latitude `lat` and east longitude `lon` (rad),
/// height `h` (m) and Julian epoch `epoch` (UT, as from [`geopotential::julian_epoch`]), m/s²:
/// `(moon, sun, total)`, positive upward. A tide that pulls up lowers gravity, so it is also the
/// correction to add to a reading.
///
/// [Longman's (1959)](https://doi.org/10.1029/JZ064i012p02351) closed form from the mean
/// orbital elements: the lunar tide to degree 3, the solar to degree 2, both times
/// [`delta`]`(2, 0.612, 0.303)` = 1.1575, with Longman's own cgs constants and Love numbers, as
/// [LongmanTide](https://github.com/jrleeman/LongmanTide) (MIT) codes it. Its masses make the
/// lunar tide 0.09% and the solar 0.2% larger than today's `GM`s, and its degree-2 factor also
/// scales the degree-3 term, whose own `δ₃` is 1.07. Averaged over the 18.61-year nodal cycle it
/// is `δ₂` times the permanent tide of [`permanent`], to 0.2 µGal.
pub fn longman(lat: f64, lon: f64, h: f64, epoch: f64) -> (f64, f64, f64) {
    // Longman's constants, cgs: G, the lunar and solar masses, the lunar eccentricity, the ratio
    // of the mean motions of the sun and moon, the mean distances of the moon and sun, the
    // Earth's equatorial radius, the lunar orbit's inclination to the ecliptic and the obliquity.
    let (mu, moon, sun) = (6.673e-8, 7.3537e25, 1.993e33);
    let (e, m, c, c1, a) = (0.05490, 0.074804, 3.84402e10, 1.495e13, 6.378270e8);
    let (i, omega): (f64, f64) = (0.08979719, 23.452f64.to_radians());
    // Julian centuries from 1899 December 31, 12h, J2000.0 less one century, and the hour angle
    // of the mean sun west of the station, from the days since J2000.0, 12h.
    let tc = (epoch - 2000.0) / 100.0 + 1.0;
    let day = (epoch - 2000.0) * 365.25;
    let t = TAU * (day - day.floor()) + lon;
    let (tc2, tc3) = (tc * tc, tc * tc * tc);

    // The moon: mean longitude, perigee, the sun's mean longitude and the ascending node.
    let s = 4.72000889397 + 8399.70927456 * tc + 3.45575191895e-05 * tc2 + 3.49065850399e-08 * tc3;
    let p = 5.83515162814 + 71.0180412089 * tc + 0.000180108282532 * tc2 + 1.74532925199e-07 * tc3;
    let hs = 4.88162798259 + 628.331950894 * tc + 5.23598775598e-06 * tc2;
    let node = 4.52360161181 - 33.757146295 * tc + 3.6264063347e-05 * tc2 + 3.39369576777e-08 * tc3;
    let incl = (omega.cos() * i.cos() - omega.sin() * i.sin() * node.cos()).acos();
    let nu = (i.sin() * node.sin() / incl.sin()).asin();
    let chi = t + hs - nu;
    let (cos_alpha, sin_alpha) = (node.cos() * nu.cos() + node.sin() * nu.sin() * omega.cos(), omega.sin() * node.sin() / incl.sin());
    let sigma = s - (node - 2.0 * (sin_alpha / (1.0 + cos_alpha)).atan());
    let l = sigma + 2.0 * e * (s - p).sin() + 1.25 * e * e * (2.0 * (s - p)).sin()
        + 3.75 * m * e * (s - 2.0 * hs + p).sin() + 11.0 / 8.0 * m * m * (2.0 * (s - hs)).sin();

    // The sun: perigee, the Earth's eccentricity and the sun's longitude.
    let p1 = 4.90822941839 + 0.0300025492114 * tc + 7.85398163397e-06 * tc2 + 5.3329504922e-08 * tc3;
    let e1 = 0.01675104 - 0.00004180 * tc - 0.000000126 * tc2;
    let chi1 = t + hs;
    let l1 = hs + 2.0 * e1 * (hs - p1).sin();

    // Zenith angles, the station's radius and the two distances.
    let (sl, cl) = lat.sin_cos();
    let cos_theta = sl * incl.sin() * l.sin() + cl * ((0.5 * incl).cos().powi(2) * (l - chi).cos() + (0.5 * incl).sin().powi(2) * (l + chi).cos());
    let cos_phi = sl * omega.sin() * l1.sin() + cl * ((0.5 * omega).cos().powi(2) * (l1 - chi1).cos() + (0.5 * omega).sin().powi(2) * (l1 + chi1).cos());
    let r = (1.0 / (1.0 + 0.006738 * sl * sl)).sqrt() * a + h * 100.0;
    let (ap, ap1) = (1.0 / (c * (1.0 - e * e)), 1.0 / (c1 * (1.0 - e1 * e1)));
    let d = 1.0 / (1.0 / c + ap * e * (s - p).cos() + ap * e * e * (2.0 * (s - p)).cos()
        + 15.0 / 8.0 * ap * m * e * (s - 2.0 * hs + p).cos() + ap * m * m * (2.0 * (s - hs)).cos());
    let big_d = 1.0 / (1.0 / c1 + ap1 * e1 * (hs - p1).cos());

    let gm = mu * moon * r / (d * d * d) * (3.0 * cos_theta * cos_theta - 1.0)
        + 1.5 * mu * moon * r * r / (d * d * d * d) * (5.0 * cos_theta.powi(3) - 3.0 * cos_theta);
    let gs = mu * sun * r / (big_d * big_d * big_d) * (3.0 * cos_phi * cos_phi - 1.0);
    // Gal to m/s², times the gravimetric factor.
    let k = 1e-2 * delta(2, 0.612, 0.303);
    (gm * k, gs * k, (gm + gs) * k)
}

/// The permanent tide at geodetic latitude `lat` (rad) and ellipsoidal height `h` (m):
/// `(W_T, g_T)`, the time-averaged tide-generating potential (m²/s²) and its part of gravity along
/// the ellipsoidal normal (m/s²), from [`PERMANENT_A`] (Mäkinen 2021, Eq. 18).
///
/// [Mäkinen (2021)](https://doi.org/10.1007/s00190-021-01541-5): `g_T` is what to add to
/// zero-tide gravity (the IAG convention) for mean-tide gravity, −30.49 µGal at the equator and
/// +60.77 at the poles (his Eq. 20); `W_T/γ`, with `γ` normal gravity, is the height of the
/// mean-tide geoid above the zero-tide one, +99.4 mm at the equator and −196.4 mm at the poles
/// (Eq. 21). In Cartesian coordinates `W_T = (A/3a²)(2z² − p²)`, so this is exact for any height,
/// and his polynomials in `sin²φ` (Eqs. 19–21) agree with it to their printed digits. A tide-free
/// value differs from the zero-tide one by the Earth's response: `k` times `W_T`, and `δ − 1`
/// times `g_T` for the [`delta`] its tide correction used.
///
/// Mean-tide gravity cannot be taken from a field whose `C̄20` was moved to
/// [`TideSystem::MeanTide`](geopotential::TideSystem::MeanTide): there `W_T` is a degree-2
/// term of the Earth's own field, falling as `r⁻³` instead of growing as `r²`, which matches it on
/// the sphere `r = R` but gives +45.7 µGal at the equator instead of −30.5.
pub fn permanent(lat: f64, h: f64) -> (f64, f64) {
    let (r, psi) = geocentric(lat, h);
    let (p, z) = (r * psi.cos(), r * psi.sin());
    let k = PERMANENT_A / (3.0 * WGS84_A * WGS84_A);
    // g_T = −∂W_T/∂h, along the normal (cos φ, sin φ) in the meridian plane.
    (k * (2.0 * z * z - p * p), 2.0 * k * (p * lat.cos() - 2.0 * z * lat.sin()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use geopotential::normal::GRS80;
    use geopotential::{Coefficients, TideSystem, julian_epoch};
    use metrology::GM_EARTH;

    const MICROGAL: f64 = crate::local::MICROGAL;

    #[test]
    fn longmans_tide() {
        // LongmanTide's test vector (tests/test_basics.py), mGal. The Julian epoch holds the time
        // to about 4 µs, which is 3×10⁻¹¹ mGal.
        let (lat, lon) = (40.7914f64.to_radians(), 282.1414f64.to_radians());
        let (moon, sun, total) = longman(lat, lon, 370.0, julian_epoch(2015, 4, 23, 0, 0));
        for (got, want) in [(moon, 0.0324029651226), (sun, -0.0288682178454), (total, 0.00353474727722)] {
            assert!((got / 1e-5 - want).abs() < 1e-10, "{} against {want}", got / 1e-5);
        }
    }

    #[test]
    fn averaged_over_the_nodal_cycle_it_is_the_permanent_tide() {
        // The mean of the routine tide is not zero (Honkasalo 1964): over 18.61 years Longman's
        // tide is δ₂ times −g_T, Mäkinen's average from the spectral expansions, to 0.2 µGal (0.19
        // at 45°).
        let step = 1.37 / 24.0 / 365.25; // years, incommensurate with the day and month
        let n = (18.61 / step) as usize;
        for deg in [0.0, 45.0, 90.0] {
            let lat = f64::to_radians(deg);
            let mean = (0..n).map(|k| longman(lat, 0.0, 0.0, 2000.0 + k as f64 * step).2).sum::<f64>() / n as f64;
            let want = -delta(2, 0.612, 0.303) * permanent(lat, 0.0).1;
            assert!((mean - want).abs() < 0.25 * MICROGAL, "{deg}°: {} against {} µGal", mean / MICROGAL, want / MICROGAL);
        }
    }

    #[test]
    fn love_number_factors() {
        // Longman's own Love numbers, and the IERS Conventions 2010 nominal values: h₂ = 0.6078
        // (ch. 7), k₂₀ = 0.30190 (Table 6.3), h₃ = 0.292 (ch. 7), k₃ = 0.093 (Table 6.3).
        assert!((delta(2, 0.612, 0.303) - 1.1575).abs() < 1e-15);
        assert!((delta(2, 0.6078, 0.30190) - 1.15495).abs() < 1e-14);
        assert!((delta(2, 0.6078, 0.30190) - 1.16).abs() < 0.01); // the conventional factor
        assert!((delta(3, 0.292, 0.093) - 1.070_666_666_666_667).abs() < 1e-15);
        // Mäkinen's Eq. 52, 1 + k − h ≈ 0.69, with k ≈ 0.30 and the nominal h.
        assert!((gamma(0.6078, 0.30190) - 0.6941).abs() < 1e-12 && (gamma(0.6078, 0.30190) - 0.69).abs() < 0.005);
        // A rigid Earth (h = k = 0) shows the tide as it is; a fluid one (h = 1 + k) shows no tilt.
        assert_eq!((delta(2, 0.0, 0.0), gamma(0.0, 0.0), gamma(1.25, 0.25)), (1.0, 1.0, 0.0));
    }

    #[test]
    fn makinens_permanent_tide() {
        // Eq. 18's A follows from the IERS H₀ = −0.31460 m (Eq. 16): A = (3/2)H₀√(5/4π)g_e(a/R_e)².
        let a_iers = 1.5 * -0.31460 * (5.0 / (4.0 * core::f64::consts::PI)).sqrt() * 9.79828685 * (WGS84_A / 6378136.55f64).powi(2);
        assert!((a_iers - PERMANENT_A).abs() < 5e-5, "{a_iers}");
        // Eqs. 19–21 in sin²φ, to their printed digits, with GRS80 normal gravity.
        for deg in (0..=90).step_by(5) {
            let lat = f64::to_radians(f64::from(deg));
            let s2 = lat.sin().powi(2);
            let (w, g) = permanent(lat, 0.0);
            assert!((w - (0.9722 - 2.8841 * s2 - 0.0195 * s2 * s2)).abs() < 5e-5, "{deg}°: W_T {w}");
            assert!((g / MICROGAL - (-30.49 + 90.95 * s2 + 0.31 * s2 * s2)).abs() < 0.005, "{deg}°: g_T {}", g / MICROGAL);
            assert!((w / GRS80.gamma(lat) * 1e3 - (99.40 - 295.41 * s2 - 0.42 * s2 * s2)).abs() < 0.005, "{deg}°: H_T");
        }
        // The spherical form −2A/a·(sin²φ − 1/3) at the equator: −30.485 µGal; the mean-tide geoid
        // 9.9 cm above the zero-tide one there.
        assert!((permanent(0.0, 0.0).1 / MICROGAL + 30.485).abs() < 1e-3);
        assert!((permanent(0.0, 0.0).0 / GRS80.gamma(0.0) - 0.0994).abs() < 1e-4);
        // Eq. 19's height factor 1 + 2h/a, to 4×10⁻⁶ at 1 km.
        let lat = 45f64.to_radians();
        assert!((permanent(lat, 1e3).0 / permanent(lat, 0.0).0 - (1.0 + 2e3 / WGS84_A)).abs() < 5e-6);
    }

    #[test]
    fn the_mean_tide_c20_is_right_only_on_its_sphere() {
        // Moving C̄20 to the mean-tide system adds W_T as part of the Earth's own field. At r = R on
        // the equator its potential is Mäkinen's to 2×10⁻⁴ (EGM's −3.11080×10⁻⁸/√5 against
        // IERS's A₀H₀ = −1.39141×10⁻⁸); its gravity is 3W/r where the tide's is −2W/r.
        let mut zero = Coefficients::zeros(GM_EARTH, WGS84_A, 2);
        zero.tide_system = TideSystem::ZeroTide;
        let mean = zero.to_tide_system(TideSystem::MeanTide).unwrap();
        let v = |c: &Coefficients, r: f64| c.partials(r, 0.0, 0.0).v;
        let added = v(&mean, WGS84_A) - v(&zero, WGS84_A);
        let (w, g) = permanent(0.0, 0.0);
        assert!((added / w - 1.0).abs() < 2e-4, "{}", added / w);
        let dr = 1.0;
        let pull = -((v(&mean, WGS84_A + dr) - v(&zero, WGS84_A + dr)) - added) / dr;
        assert!((pull / g + 1.5).abs() < 1e-3 && (pull / MICROGAL - 45.7).abs() < 0.1, "{}", pull / MICROGAL);
    }
}
