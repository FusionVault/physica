#![doc = include_str!("../README.md")]

// Exact SI defining constants (2019).

/// Speed of light in vacuum, m/s (exact).
pub const C: f64 = 299_792_458.0;
/// Planck constant, J s (exact).
pub const H: f64 = 6.626_070_15e-34;
/// Reduced Planck constant `h/2π`, J s.
pub const HBAR: f64 = 1.054_571_817e-34;
/// Elementary charge, C (exact).
pub const E_CHARGE: f64 = 1.602_176_634e-19;
/// Boltzmann constant, J/K (exact).
pub const K_B: f64 = 1.380_649e-23;
/// Avogadro constant, mol⁻¹ (exact).
pub const N_A: f64 = 6.022_140_76e23;
/// Molar gas constant `N_A k_B`, J mol⁻¹ K⁻¹.
pub const R_GAS: f64 = 8.314_462_618;

// CODATA 2018 measured constants.

/// Newtonian constant of gravitation, m³ kg⁻¹ s⁻² (uncertainty 22 ppm).
pub const G: f64 = 6.674_30e-11;
/// Fine-structure constant.
pub const ALPHA: f64 = 7.297_352_569_3e-3;
/// Vacuum magnetic permeability, N/A².
pub const MU_0: f64 = 1.256_637_062_12e-6;
/// Vacuum electric permittivity, F/m.
pub const EPS_0: f64 = 8.854_187_812_8e-12;
/// Electron mass, kg.
pub const M_E: f64 = 9.109_383_701_5e-31;
/// Proton mass, kg.
pub const M_P: f64 = 1.672_621_923_69e-27;
/// Atomic mass constant (one dalton), kg.
pub const AMU: f64 = 1.660_539_066_60e-27;
/// Mass of a hydrogen atom, kg.
pub const M_H: f64 = 1.007_825_032_23 * AMU;
/// Einstein's coupling constant `8πG/c⁴`, s² kg⁻¹ m⁻¹: the stiffness of spacetime in
/// `G_μν = κT_μν`. Its reciprocal, 4.8×10⁴² N, is why laboratory stress-energy curves spacetime
/// so little.
pub const EINSTEIN_KAPPA: f64 = 8.0 * core::f64::consts::PI * G / (C * C * C * C);
/// Bohr magneton `eħ/2m_e`, J/T.
pub const MU_B: f64 = 9.274_010_078_3e-24;
/// Proton gyromagnetic ratio, s⁻¹ T⁻¹ (γ_p/2π = 42.577 MHz/T, the MRI Larmor frequency).
pub const GAMMA_P: f64 = 2.675_221_874_4e8;
/// Planck length `√(ħG/c³)`, m.
pub const PLANCK_LENGTH: f64 = 1.616_255e-35;

// Energy units.

/// One electronvolt, J.
pub const EV: f64 = E_CHARGE;
/// One kiloelectronvolt, J.
pub const KEV: f64 = 1e3 * E_CHARGE;

// Earth and astronomy.

/// Standard acceleration of gravity, m/s² (exact, conventional).
pub const G0: f64 = 9.806_65;
/// Earth's gravitational parameter GM (IERS 2010, including the atmosphere), m³/s².
pub const GM_EARTH: f64 = 3.986_004_418e14;
/// WGS84 ellipsoid semi-major axis, m.
pub const WGS84_A: f64 = 6_378_137.0;
/// WGS84 ellipsoid flattening.
pub const WGS84_F: f64 = 1.0 / 298.257_223_563;
/// Nominal solar mass parameter GM☉ (IAU 2015 Resolution B3), m³/s².
pub const GM_SUN: f64 = 1.327_124_4e20;
/// Solar mass `GM☉/G`, kg.
pub const SOLAR_MASS: f64 = GM_SUN / G;
/// Astronomical unit, m (exact, IAU 2012).
pub const AU: f64 = 149_597_870_700.0;
/// Parsec, m.
pub const PARSEC: f64 = 3.085_677_581_491_367e16;
/// Megaparsec, m.
pub const MPC: f64 = 1e6 * PARSEC;
/// Julian year, s.
pub const YEAR: f64 = 365.25 * 86_400.0;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn derived_constants_agree() {
        assert!((HBAR / (H / (2.0 * core::f64::consts::PI)) - 1.0).abs() < 1e-9);
        assert!((R_GAS / (N_A * K_B) - 1.0).abs() < 1e-9);
        assert!((MU_0 * EPS_0 * C * C - 1.0).abs() < 1e-9);
        assert!((ALPHA / (E_CHARGE * E_CHARGE / (4.0 * core::f64::consts::PI * EPS_0 * HBAR * C)) - 1.0).abs() < 1e-9);
        assert!((PLANCK_LENGTH / (HBAR * G / C.powi(3)).sqrt() - 1.0).abs() < 1e-5);
        assert!((SOLAR_MASS / 1.988_4e30 - 1.0).abs() < 1e-4);
        assert!((MU_B / (E_CHARGE * HBAR / (2.0 * M_E)) - 1.0).abs() < 1e-9);
        assert!((EINSTEIN_KAPPA / 2.076_647e-43 - 1.0).abs() < 1e-6);
        assert!((PARSEC / (AU * 648_000.0 / core::f64::consts::PI) - 1.0).abs() < 1e-12);
    }
}
