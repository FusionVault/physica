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
/// Earth's nominal mean angular velocity, rad/s (IERS Conventions 2010, Table 1.1; also the GRS80
/// and WGS84 defining value). GPS uses 7.292 115 146 7×10⁻⁵ (IS-GPS-200).
pub const EARTH_ROTATION: f64 = 7.292_115e-5;
/// The Moon-to-Earth mass ratio (IERS Conventions 2010, Table 1.1).
pub const MOON_EARTH_MASS_RATIO: f64 = 0.012_300_037_1;
/// The Moon's gravitational parameter, m³/s²: [`GM_EARTH`] × [`MOON_EARTH_MASS_RATIO`],
/// 4.902 800 2×10¹² and, like `GM_EARTH`, TT-compatible. The JPL DE430 ephemeris gives
/// 4 902.800 066 km³/s² on its own time scale, 3×10⁻⁸ lower.
pub const GM_MOON: f64 = GM_EARTH * MOON_EARTH_MASS_RATIO;
/// WGS84 ellipsoid semi-major axis, m.
pub const WGS84_A: f64 = 6_378_137.0;
/// WGS84 ellipsoid flattening.
pub const WGS84_F: f64 = 1.0 / 298.257_223_563;
/// GRS80 flattening (Moritz 2000), derived there from its defining `J₂`; it differs from WGS84's
/// in the ninth digit. GRS80's semi-major axis and rotation rate are [`WGS84_A`] and
/// [`EARTH_ROTATION`].
pub const GRS80_F: f64 = 1.0 / 298.257_222_101;
/// GRS80's geocentric gravitational constant (Moritz 2000), m³/s²: the 1980 value, 1.5×10⁻⁷ above
/// [`GM_EARTH`], and the `μ` of the GPS interface specification.
pub const GRS80_GM: f64 = 3.986_005e14;
/// Earth's spin angular momentum `Iω`, kg m²/s: 5.86×10³³, the value the LAGEOS and LARES
/// frame-dragging analyses use (Ciufolini & Pavlis 2004).
pub const EARTH_ANGULAR_MOMENTUM: f64 = 5.86e33;
/// The geoid potential `W₀` that Terrestrial Time's rate was defined from (IAU 2000 Resolution
/// B1.9; IERS Conventions 2010, Table 1.1), m²/s². For time scales and clock rates.
pub const W0_IAU: f64 = 62_636_856.0;
/// The reference potential of the International Height Reference System (IAG 2015 Resolution 1),
/// m²/s²: 2.6 m²/s², about 26 cm, below [`W0_IAU`]. For heights.
pub const W0_IHRS: f64 = 62_636_853.4;
/// `L_G = 1 − dTT/dTCG`, the defining constant of Terrestrial Time (IAU 2000 Resolution B1.9):
/// [`W0_IAU`]`/c²` to its ten digits.
pub const L_G: f64 = 6.969_290_134e-10;
/// `L_B = 1 − dTDB/dTCB` (IAU 2006 Resolution B3).
pub const L_B: f64 = 1.550_519_768e-8;
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
        assert!((GM_MOON / 4.902_800_222e12 - 1.0).abs() < 1e-9);
        assert!((GM_MOON / 4.902_800_066e12 - 1.0).abs() < 5e-8); // JPL DE430, its own time scale
        assert!((1.0 / MOON_EARTH_MASS_RATIO - 81.300_57).abs() < 1e-4);
        // A sidereal day is 2π/ω: 86 164.1 s.
        assert!((2.0 * core::f64::consts::PI / EARTH_ROTATION - 86_164.1).abs() < 0.1);
        assert!((HBAR / (H / (2.0 * core::f64::consts::PI)) - 1.0).abs() < 1e-9);
        assert!((R_GAS / (N_A * K_B) - 1.0).abs() < 1e-9);
        assert!((MU_0 * EPS_0 * C * C - 1.0).abs() < 1e-9);
        assert!((ALPHA / (E_CHARGE * E_CHARGE / (4.0 * core::f64::consts::PI * EPS_0 * HBAR * C)) - 1.0).abs() < 1e-9);
        assert!((PLANCK_LENGTH / (HBAR * G / C.powi(3)).sqrt() - 1.0).abs() < 1e-5);
        assert!((SOLAR_MASS / 1.988_4e30 - 1.0).abs() < 1e-4);
        assert!((MU_B / (E_CHARGE * HBAR / (2.0 * M_E)) - 1.0).abs() < 1e-9);
        assert!((EINSTEIN_KAPPA / 2.076_647e-43 - 1.0).abs() < 1e-6);
        assert!((PARSEC / (AU * 648_000.0 / core::f64::consts::PI) - 1.0).abs() < 1e-12);
        assert!((L_G - W0_IAU / (C * C)).abs() < 1e-19 && (W0_IAU - W0_IHRS - 2.6).abs() < 1e-8);
        assert!((GRS80_GM / GM_EARTH - 1.0 - 1.46e-7).abs() < 1e-9 && GRS80_F > WGS84_F);
    }
}
