//! Clocks as gravimeters. A clock higher in the Earth's field runs faster by `gΔh/c²`, 1.1×10⁻¹⁶
//! per metre, so comparing two clocks measures their geopotential difference, and an optical
//! clock at 10⁻¹⁸ resolves a centimetre of height (chronometric levelling). The same physics
//! sets the rates of the relativistic time scales and the corrections a GNSS receiver applies.
//!
//! `W` is the geodesist's geopotential, largest at the geoid and smaller higher up, so the higher
//! clock has the smaller `W` and the higher frequency.
//!
//! ```
//! use geopotential::normal::GRS80;
//! use gravimetry::clocks::{frequency_shift, height_difference};
//!
//! // Two optical clocks compared over a fibre differ by 3.2 parts in 10¹⁷: the faster one is
//! // 29 cm higher.
//! let g = GRS80.gamma(47f64.to_radians());
//! assert!((height_difference(3.2e-17, g) - 0.2933).abs() < 1e-4);
//! // The Tokyo Skytree: 452.6 m of height is 4.93×10⁻¹⁴, 21.18 Hz on strontium's clock line.
//! let shift = frequency_shift(-9.797_248 * 452.596);
//! assert!((shift * 429_228_004_229_873.0 - 21.18).abs() < 0.01);
//! ```

use metrology::{C, EARTH_ROTATION, GM_EARTH};

const C2: f64 = C * C;

/// The fractional frequency `Δf/f` of a clock whose geopotential is `dw = W_clock − W_ref`
/// (m²/s²) above a reference clock's: `−ΔW/c²`. Against Terrestrial Time, `dw = W −`
/// [`metrology::W0_IAU`]: zero on the geoid, −1.1×10⁻¹⁶ a metre below it.
pub fn frequency_shift(dw: f64) -> f64 {
    -dw / C2
}

/// The geopotential difference `W_B − W_A` (m²/s²) that a measured `(f_B − f_A)/f` means:
/// `−c²Δf/f`. 10⁻¹⁶ is 9 m²/s².
pub fn potential_difference(df: f64) -> f64 {
    -C2 * df
}

/// The height `h_B − h_A` (m) of clock B above clock A from their `(f_B − f_A)/f` and the local
/// gravity `g` (m/s²): `c²(Δf/f)/g`. A clock's uncertainty propagates the same way: 10⁻¹⁸ is 0.92 cm.
pub fn height_difference(df: f64, g: f64) -> f64 {
    C2 * df / g
}

/// The rate of a moving clock against Geocentric Coordinate Time, `dτ/dTCG − 1 = −(U + v²/2)/c²`,
/// for gravitational potential `u` (m²/s², positive) and coordinate speed `v` (m/s): the
/// gravitational redshift and the second-order Doppler shift. For an Earth-fixed clock `v = ωr cos φ`
/// makes `v²/2` the centrifugal potential, so `U + v²/2 = W`, and a geoid clock's rate is
/// `−`[`metrology::L_G`].
pub fn rate_tcg(u: f64, v: f64) -> f64 {
    -(u + 0.5 * v * v) / C2
}

/// Proper time minus coordinate time (s) a clock accumulates along samples of `(t, U, v)` (s,
/// m²/s², m/s), the trapezoidal integral of [`rate_tcg`]; zero for fewer than two samples.
pub fn proper_time_offset(samples: &[(f64, f64, f64)]) -> f64 {
    samples.windows(2).map(|w| 0.5 * (rate_tcg(w[0].1, w[0].2) + rate_tcg(w[1].1, w[1].2)) * (w[1].0 - w[0].0)).sum()
}

/// The fractional offset GPS satellite clocks are built with, so that they read 10.23 MHz on
/// orbit against a geoid clock: −4.4647×10⁻¹⁰, the net of the gravitational blueshift and the
/// orbital Doppler redshift (IS-GPS-200, §3.3.1.1).
pub const GPS_CLOCK_RATE_OFFSET: f64 = -4.4647e-10;

/// `F = −2√μ/c²`, s/√m, the coefficient of the GNSS eccentricity correction (IS-GPS-200,
/// §20.3.3.3.3.1), with its `μ` the GRS80 value [`metrology::GRS80_GM`].
pub const GNSS_ECCENTRICITY_F: f64 = -4.442_807_633e-10;

/// The periodic relativistic clock correction `Δt_r = F e √a sin E` (s) for a satellite with
/// semi-major axis `a` (m), eccentricity `e` and eccentric anomaly `big_e` (rad): up to ±23 ns for
/// GPS, added to the satellite clock reading.
pub fn eccentricity_correction(a: f64, e: f64, big_e: f64) -> f64 {
    GNSS_ECCENTRICITY_F * e * a.sqrt() * big_e.sin()
}

/// The Sagnac correction `(ω/c²)(x_e y_r − y_e x_r)` (s) for a signal from `emitter` to
/// `receiver`, both Earth-fixed `[x, y, z]` (m), in the rotating frame: up to 133 ns for a ground
/// receiver, with the sign of the east–west geometry. GPS's ω differs from the IERS value used
/// here by 2×10⁻⁹.
pub fn sagnac(emitter: [f64; 3], receiver: [f64; 3]) -> f64 {
    EARTH_ROTATION / C2 * (emitter[0] * receiver[1] - emitter[1] * receiver[0])
}

/// The Shapiro delay (s) of light past a mass `gm` (m³/s²) between points `r1` and `r2` (m) from
/// it, `r12` apart: `(2GM/c³) ln((r₁ + r₂ + r₁₂)/(r₁ + r₂ − r₁₂))`. 262 µs for Cassini's round trip
/// past the Sun, tens of picoseconds between a GNSS satellite and the ground.
pub fn shapiro_delay(gm: f64, r1: f64, r2: f64, r12: f64) -> f64 {
    2.0 * gm / (C2 * C) * ((r1 + r2 + r12) / (r1 + r2 - r12)).ln()
}

/// The one-way coordinate time (s) a signal takes from `emitter` to `receiver` (Earth-fixed
/// `[x, y, z]`, m): the light time plus the Earth's [`shapiro_delay`] and the [`sagnac`] term. In a
/// two-way link the first two cancel by reciprocity and the Sagnac term, which changes sign with
/// the direction, is the residual.
pub fn one_way_transfer(emitter: [f64; 3], receiver: [f64; 3]) -> f64 {
    let norm = |v: [f64; 3]| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    let d = norm([receiver[0] - emitter[0], receiver[1] - emitter[1], receiver[2] - emitter[2]]);
    d / C + shapiro_delay(GM_EARTH, norm(emitter), norm(receiver), d) + sagnac(emitter, receiver)
}

#[cfg(test)]
mod tests {
    use super::*;
    use geopotential::normal::GRS80;
    use metrology::{AU, G0, GM_SUN, GRS80_GM, L_G, W0_IAU, W0_IHRS};

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    #[test]
    fn chronometric_levelling() {
        // The headline numbers: 1 cm is 1.09×10⁻¹⁸ at g₀, and 10⁻¹⁸ is 0.92 cm; the three
        // conversions invert one another.
        assert!(close(frequency_shift(-G0 * 0.01), 1.0912e-18, 1e-21) && close(height_difference(1e-18, G0), 0.009_17, 1e-5));
        assert!(close(frequency_shift(potential_difference(3.2e-17)), 3.2e-17, 1e-30));
        // Terrestrial Time: a geoid clock runs at the TT rate, and 1 m below it 1.09×10⁻¹⁶ slow;
        // against TCG the geoid clock's rate is −L_G by definition. The IHRS W₀ sits 26 cm higher.
        assert_eq!(frequency_shift(W0_IAU - W0_IAU), 0.0);
        assert!(close(frequency_shift(9.8), -1.09e-16, 1e-18) && close(rate_tcg(W0_IAU, 0.0), -L_G, 1e-19));
        assert!(close(height_difference(frequency_shift(W0_IHRS - W0_IAU), G0), 0.265, 1e-3));
    }

    #[test]
    fn the_tokyo_skytree() {
        // Takamoto et al. (2020; Shinkai, Takamoto & Katori 2025, arXiv:2502.06104): two Sr clocks
        // Δh = 452.596(13) m apart, g = 9.797248(24) m/s², Δν ≈ 21.18 Hz on 429 228 004 229 873 Hz.
        // GRS80 normal gravity at 35.71°N, 230 m up, is that g to 2 parts in 10⁶.
        let (dh, g, nu) = (452.596, 9.797_248, 429_228_004_229_873.0);
        assert!(close(frequency_shift(-g * dh) * nu, 21.18, 0.01));
        assert!(close(height_difference(frequency_shift(-g * dh), g), dh, 1e-9));
        let gamma = GRS80.gamma_at(35.71f64.to_radians(), 230.0);
        assert!(close(gamma, g, 2e-5), "{gamma}");
    }

    #[test]
    fn gradients_across_a_sample_and_the_suns_redshift() {
        // Bothwell et al. (2022, arXiv:2109.12238): 1.09×10⁻¹⁹ per mm at JILA's g = 9.796,
        // measured 0.98(23). Zheng et al. (2023, arXiv:2207.07145): 10.9×10⁻¹⁹ per cm at
        // g = 9.803, measured 12.4 ± 0.7 ± 2.5, 0.57σ away.
        let jila = frequency_shift(-9.796 * 1e-3);
        assert!(close(jila, 1.09e-19, 5e-22) && (jila - 0.98e-19).abs() < 2.3e-20);
        let wisconsin = frequency_shift(-9.803 * 1e-2);
        assert!(close(wisconsin, 10.9e-19, 5e-21) && (((12.4e-19 - wisconsin) / 0.7f64.hypot(2.5) / 1e-19) - 0.575).abs() < 0.001);
        // González Hernández et al. (2020, arXiv:2009.10558): the Sun's lines, as a speed, are
        // 636.31 m/s at the photosphere, 633.35 at 1 AU and 633.10 at the Earth on the observing
        // date (147.75 million km, less the Earth's own 0.21; the paper calls that distance
        // correction "+4 cm/s", but it lowers the value); measured 638 ± 6.
        let (r_sun, r_earth) = (6.957e8, 6.3781e6);
        let speed = |w_clock: f64, w_ref: f64| -C * frequency_shift(w_clock - w_ref);
        assert!(close(speed(GM_SUN / r_sun, 0.0), 636.31, 0.005) && close(speed(GM_SUN / r_sun, GM_SUN / AU), 633.35, 0.005));
        let seen = speed(GM_SUN / r_sun, GM_SUN / 147.75e9 + GM_EARTH / r_earth);
        assert!(close(seen, 633.10, 0.01) && (638.0 - seen).abs() < 6.0);
    }

    #[test]
    fn gnss_corrections() {
        // IS-GPS-200: F from the GRS80 μ to its ten digits, and the clock offset of a satellite at
        // a = 26 561.75 km (√a = 5153.8 m^½) against a geoid clock: 4.4647×10⁻¹⁰, to 3×10⁻¹⁴.
        assert!(close(GNSS_ECCENTRICITY_F, -2.0 * GRS80_GM.sqrt() / C2, 1e-19));
        let a = 5153.8f64.powi(2);
        let orbit = rate_tcg(GM_EARTH / a, (GM_EARTH / a).sqrt()) - rate_tcg(W0_IAU, 0.0);
        assert!(close(orbit, -GPS_CLOCK_RATE_OFFSET, 3e-14), "{orbit:e}");
        // e = 0.01 at E = 90°: −22.9 ns. Sagnac for a satellite 20 000 km east of the receiver: 133 ns
        // is its ceiling; the Earth's Shapiro delay to the ground is tens of picoseconds.
        assert!(close(eccentricity_correction(a, 0.01, core::f64::consts::FRAC_PI_2), -22.9e-9, 0.05e-9));
        let (sat, rcv) = ([6_378_137.0, 20_000e3, 0.0], [6_378_137.0, 0.0, 0.0]);
        let s = sagnac(sat, rcv);
        assert!(s < -10e-9 && s > -133e-9 && sagnac(rcv, sat) == -s);
        let light = 20_000e3 / C;
        let shapiro = one_way_transfer(sat, rcv) - light - s;
        assert!(shapiro > 1e-12 && shapiro < 1e-10, "{shapiro:e}");
    }

    #[test]
    fn cassini_and_proper_time() {
        // Cassini's 2002 conjunction (Will 2014, arXiv:1403.7377): the spacecraft at 8.43 AU, the
        // ray 1.6 R☉ from the Sun; twice the one-way delay is 262.4 µs, and Will's Eq. 63 shorthand
        // ½(1 + γ)[240 − 20 ln(d²/r)] µs gives 263.8.
        let b: f64 = 1.6 * 6.957e8;
        let (r1, r2) = (AU, 8.43 * AU);
        let r12 = (r1 * r1 - b * b).sqrt() + (r2 * r2 - b * b).sqrt();
        let round_trip = 2.0 * shapiro_delay(GM_SUN, r1, r2, r12) * 1e6;
        assert!(close(round_trip, 262.42, 0.01) && close(round_trip / (240.0 - 20.0 * (1.6f64.powi(2) / 8.43).ln()), 1.0, 0.006));
        // A GPS clock over 60 s: −2.5×10⁻¹⁰ × 60 s = −15 ns against TCG; constant rate, so the
        // trapezoid is exact.
        let a = 26_561.75e3;
        let (u, v) = (GM_EARTH / a, (GM_EARTH / a).sqrt());
        let offset = proper_time_offset(&[(0.0, u, v), (20.0, u, v), (60.0, u, v)]);
        assert!(close(offset, 60.0 * rate_tcg(u, v), 1e-24) && close(offset, -1.503e-8, 0.002e-8));
        assert_eq!(proper_time_offset(&[(0.0, u, v)]), 0.0);
    }
}
