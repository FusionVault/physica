//! Waves in the atmosphere that a gravity satellite can feel: sound, buoyancy ("gravity") waves
//! and the Lamb wave, and the mass a pressure change puts above a square metre of ground.
//!
//! The isothermal acoustic-gravity dispersion relation is
//! `k_z² = (ω² − ω_a²)/c² − k_h²(1 − N²/ω²)`, with the acoustic cut-off `ω_a = γg/(2c)` and the
//! buoyancy (Brunt–Väisälä) frequency `N = √(γ−1)·g/c`. A Lamb wave runs along the ground at the
//! speed of sound and decays upward. The 2022 Hunga Tonga eruption sent one around the world at
//! 318 m/s ([Wright et al. 2022](https://doi.org/10.1038/s41586-022-05012-5)).
//!
//! ```
//! use gravimetry::atmosphere::{acoustic_cutoff, sound_speed};
//!
//! let c = sound_speed(250.0); // 317.0 m/s
//! let period_minutes = 2.0 * std::f64::consts::PI / acoustic_cutoff(c) / 60.0;
//! assert!((period_minutes - 4.83).abs() < 0.01);
//! ```

use metrology::R_GAS;

/// Ratio of specific heats of dry air.
pub const GAMMA: f64 = 1.4;
/// Molar mass of dry air, kg/mol.
pub const MOLAR_MASS: f64 = 0.028_964_4;
/// Standard gravity, [`metrology::G0`] (this module once carried 9.81, 0.03% more).
pub const G0: f64 = metrology::G0;

/// Speed of sound in dry air at temperature `t` (K): `c = √(γRT/M)`, m/s.
pub fn sound_speed(t: f64) -> f64 {
    (GAMMA * R_GAS * t / MOLAR_MASS).sqrt()
}

/// Acoustic cut-off frequency `ω_a = γg/(2c)`, rad/s.
pub fn acoustic_cutoff(c: f64) -> f64 {
    GAMMA * G0 / (2.0 * c)
}

/// Buoyancy (Brunt–Väisälä) frequency of an isothermal atmosphere, `N = √(γ−1)·g/c`, rad/s.
pub fn buoyancy_frequency(c: f64) -> f64 {
    (GAMMA - 1.0).sqrt() * G0 / c
}

/// Squared vertical wavenumber of a wave of angular frequency `omega` and horizontal wavenumber
/// `kh` in an isothermal atmosphere with sound speed `c`. Positive: the wave propagates upward;
/// negative: it is trapped (evanescent).
pub fn vertical_wavenumber_sq(omega: f64, kh: f64, c: f64) -> f64 {
    let (wa, n) = (acoustic_cutoff(c), buoyancy_frequency(c));
    (omega * omega - wa * wa) / (c * c) - kh * kh * (1.0 - n * n / (omega * omega))
}

/// How much a wave's particle velocity grows climbing to altitude `z` (m) in an atmosphere of
/// scale height `scale_height` (m): conserving the energy flux `½ρv²c` through exponentially
/// thinning air gives `exp(z/2H)`, a factor of 1.3×10³ by 100 km for H = 7 km.
///
/// `None` above `ceiling`, where the formula stops applying: molecular viscosity and thermal
/// conduction dissipate the wave instead, so extrapolating it to 300 km (a factor of 2×10⁹) is
/// meaningless. The turbopause is near 105 km and dissipation grows quickly above it, so a
/// ceiling of 110–150 km is usual; the caller states its own.
pub fn wave_growth(z: f64, scale_height: f64, ceiling: f64) -> Option<f64> {
    (z <= ceiling).then(|| (z / (2.0 * scale_height)).exp())
}

/// Mass per square metre above the ground that a surface-pressure change `dp` (Pa) represents:
/// `Δp/g`, kg/m². A 1 hPa pulse is 10 kg/m², the same as 1 cm of water.
pub fn pressure_to_mass(dp: f64) -> f64 {
    dp / G0
}

/// Ocean mass that moves with a pressure pulse travelling at `speed` (m/s) over water `depth` (m)
/// deep, as a multiple of the pulse's own atmospheric mass: `gH/(U² − gH)`, the forced
/// shallow-water solution ([Proudman 1929](https://doi.org/10.1111/j.2044-8312.1929.tb05286.x)).
/// A slow pulse gives −1, the inverse barometer, which cancels the air's mass. A pulse faster than
/// ocean long waves (`U > √(gH)`, like the Tonga Lamb wave over the deep Pacific) raises the sea
/// instead, which adds to it. Near `U = √(gH)` the response is resonant and this formula fails.
pub fn forced_ocean_response(speed: f64, depth: f64) -> f64 {
    let gh = G0 * depth;
    gh / (speed * speed - gh)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn isothermal_values_at_250_kelvin() {
        let c = sound_speed(250.0);
        assert!((c - 316.97).abs() < 0.01, "{c}");
        assert!((acoustic_cutoff(c) - 0.02166).abs() < 1e-5);
        assert!((buoyancy_frequency(c) - 0.01957).abs() < 1e-5);
        // Wright et al.'s 318.2 m/s Lamb wave is sound in an isothermal atmosphere at about 252 K.
        assert!((sound_speed(252.0) / 318.2 - 1.0).abs() < 1e-3);
    }

    #[test]
    fn the_lamb_wave_is_trapped() {
        // ω = c·k_h: k_z² = (N² − ω_a²)/c² < 0 for air (γ = 1.4), at every wavelength.
        let c = sound_speed(250.0);
        for kh in [1e-6, 1e-5, 1e-4] {
            let kz2 = vertical_wavenumber_sq(c * kh, kh, c);
            let expect = (buoyancy_frequency(c).powi(2) - acoustic_cutoff(c).powi(2)) / (c * c);
            assert!(kz2 < 0.0 && (kz2 / expect - 1.0).abs() < 1e-9, "{kz2:e}");
        }
        // Sound well above the cut-off propagates upward.
        assert!(vertical_wavenumber_sq(0.1, 1e-5, c) > 0.0);
        assert!((pressure_to_mass(100.0) - 10.197).abs() < 0.001);
    }

    #[test]
    fn ocean_response_limits() {
        assert!((forced_ocean_response(1e-3, 4000.0) + 1.0).abs() < 1e-9); // the inverse barometer
        // Tonga's Lamb wave over 4 km and 5 km of ocean raises the sea by 0.63 and 0.94 of Δp/(ρg).
        assert!((forced_ocean_response(318.2, 4000.0) - 0.6324).abs() < 1e-4);
        assert!((forced_ocean_response(318.2, 5000.0) - 0.9390).abs() < 1e-4);
    }

    #[test]
    fn waves_grow_then_dissipate() {
        // Energy-flux conservation: a factor of 1.27e3 by 100 km for H = 7 km.
        assert!((wave_growth(100e3, 7e3, 150e3).unwrap() / 1265.0 - 1.0).abs() < 1e-3);
        // The formula is refused where dissipation takes over, instead of returning 2e9 at 300 km.
        assert!(wave_growth(300e3, 7e3, 150e3).is_none());
        assert!((wave_growth(0.0, 7e3, 150e3).unwrap() - 1.0).abs() < 1e-12);
        // The Hunga Tonga Lamb wave: 780 Pa gives 2.04 m/s of surface particle velocity.
        let v0: f64 = 780.0 / (1.2 * 318.0);
        assert!((v0 - 2.04).abs() < 0.01);
    }
}
