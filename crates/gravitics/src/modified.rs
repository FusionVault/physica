//! Alternatives to Newton and Einstein that make numbers: MOND and the radial acceleration
//! relation, short-range Yukawa forces, Einstein–Cartan torsion and Weber-like gravity.
//!
//! The radial acceleration relation is data: galaxies' observed acceleration follows their
//! baryons' with one scale, g† = 1.2×10⁻¹⁰ m/s² ([McGaugh, Lelli & Schombert 2016](https://arxiv.org/abs/1609.05917)).
//! Whether that means modified gravity is contested; wide binaries are its live test. Any
//! gravitational-strength Yukawa force must have a range below 38.6 µm
//! ([Lee et al. 2020](https://arxiv.org/abs/2002.11761)). Einstein–Cartan torsion matters only near
//! 10⁴⁷ g/cm³. Weber-like velocity-dependent gravity cannot give both Mercury's perihelion and
//! the right light bending.
//!
//! ```
//! use gravitics::modified::{G_DAGGER, rar};
//!
//! // Deep in the low-acceleration regime the RAR becomes √(g_bar g†).
//! let g = 1e-14;
//! assert!((rar(g, G_DAGGER) / (g * G_DAGGER).sqrt() - 1.0).abs() < 0.01);
//! ```

use core::f64::consts::PI;

use metrology::{C, E_CHARGE, G, HBAR};

/// The acceleration scale of the radial acceleration relation, m/s² (McGaugh et al. 2016;
/// MOND's a₀).
pub const G_DAGGER: f64 = 1.20e-10;

/// Limits on a Yukawa addition `V = V_N[1 + α e^(−r/λ)]`: range λ (m), the largest allowed
/// |α| there, and the source.
pub const YUKAWA_LIMITS: [(f64, f64, &str); 3] = [
    (38.6e-6, 1.0, "Lee et al. 2020, Eöt-Wash (95%)"),
    (56e-6, 1.0, "Kapner et al. 2007 (95%, superseded)"),
    (10e-6, 1.4e4, "Geraci et al. 2008"),
];

/// The radial acceleration relation `g_obs = g_bar/(1 − e^(−√(g_bar/g†)))` for baryonic
/// acceleration `g_bar` and scale `g_dagger` (m/s²); Newtonian well above g†.
pub fn rar(g_bar: f64, g_dagger: f64) -> f64 {
    g_bar / -(-(g_bar / g_dagger).sqrt()).exp_m1()
}

/// The baryonic Tully–Fisher velocity `(GMa₀)^(1/4)` (m/s) of a galaxy of baryonic mass `m` (kg).
pub fn btfr_velocity(m: f64, a0: f64) -> f64 {
    (G * m * a0).powf(0.25)
}

/// The MOND radius `√(GM/a₀)` (m), beyond which a mass `m` (kg) is in the low-acceleration regime:
/// about 7000 au for the Sun, where wide binaries test it.
pub fn mond_radius(m: f64, a0: f64) -> f64 {
    (G * m / a0).sqrt()
}

/// The force on a point mass with a Yukawa term, relative to Newton:
/// `1 + α(1 + r/λ)e^(−r/λ)` at separation `r`.
pub fn yukawa_factor(r: f64, alpha: f64, lambda: f64) -> f64 {
    1.0 + alpha * (1.0 + r / lambda) * (-r / lambda).exp()
}

/// The dark-energy length `(ħc/ρ)^(1/4)` (m) for vacuum energy density `rho` (J/m³): 85 µm for
/// the observed 3.8 keV/cm³, a natural place for new short-range physics.
pub fn dark_energy_length(rho: f64) -> f64 {
    (HBAR * C / rho).powf(0.25)
}

/// Mass (eV/c²) of the boson whose Yukawa range is `lambda` (m): `ħc/λ`.
pub fn boson_mass_ev(lambda: f64) -> f64 {
    HBAR * C / (lambda * E_CHARGE)
}

/// The density (kg/m³) at which Einstein–Cartan spin–torsion energy rivals mass energy for
/// fermions of mass `m` (kg): `m²c⁴/(8πGħ²)`.
pub fn einstein_cartan_density(m: f64) -> f64 {
    m * m * C.powi(4) / (8.0 * PI * G * HBAR * HBAR)
}

/// Perihelion advance of Weber-like gravity, `F = −(GMm/r²)(1 − μṙ²/c² + γ r r̈/c²)`, relative to
/// general relativity: `γ/6` for γ in units of 1/c², independent of μ. Matching Mercury needs
/// γ = 6 (Assis), and no choice also gives GR's light bending.
pub fn weber_perihelion_ratio(gamma: f64) -> f64 {
    gamma / 6.0
}

/// Fractional inertial-mass shift `qU/(3mc²)` that Weber–Assis electrodynamics predicts for a
/// charge `q` (C) of mass `m` (kg) inside a shell at potential `u` (V): 0.65% for an electron at
/// 10 kV, which experiment does not see.
pub fn weber_mass_shift(q: f64, u: f64, m: f64) -> f64 {
    q * u / (3.0 * m * C * C)
}

#[cfg(test)]
mod tests {
    use super::*;
    use metrology::{AMU, KEV, M_E, PARSEC, SOLAR_MASS};

    fn close(a: f64, b: f64, rel: f64) -> bool {
        (a / b - 1.0).abs() < rel
    }

    #[test]
    fn radial_acceleration_relation() {
        for (g, obs) in [(1e-12, 1.14621e-11), (1e-11, 3.98812e-11), (1e-10, 1.67048e-10), (1.2e-10, 1.89837e-10), (1e-9, 1.05905e-9), (1e-8, 1.00011e-8)] {
            assert!(close(rar(g, G_DAGGER), obs, 1e-5), "{g}");
        }
        assert!(close(btfr_velocity(1e11 * SOLAR_MASS, G_DAGGER), 199.767e3, 1e-5));
        let au = 1.495_978_707e11;
        assert!(close(mond_radius(2.0 * SOLAR_MASS, G_DAGGER) / au, 9941.57, 1e-5));
        assert!(close(mond_radius(SOLAR_MASS, G_DAGGER) / au, 7030.0, 1e-3));
        assert!(close(mond_radius(1e11 * SOLAR_MASS, G_DAGGER) / (1e3 * PARSEC), 10.7774, 1e-5));
    }

    #[test]
    fn short_range_forces() {
        assert!(close(yukawa_factor(1.0, 1.0, 1.0), 1.735759, 1e-6));
        assert!(close(dark_energy_length(3.8 * KEV * 1e6), 84.8889e-6, 1e-5)); // Lee et al.: ≈ 85 µm
        assert!(close(boson_mass_ev(YUKAWA_LIMITS[0].0), 5.11210e-3, 1e-5)); // ≈ 5.1 meV
    }

    #[test]
    fn torsion_and_weber() {
        // 3.6×10⁴⁷ g/cm³ for electrons, 1.2×10⁵⁴ for neutrons (Ruggiero & Tartaglia).
        assert!(close(einstein_cartan_density(M_E), 3.59305e50, 1e-5));
        assert!(close(einstein_cartan_density(1.008_664_915_95 * AMU), 1.21472e57, 1e-5));
        // Tisserand (γ = 2), Weber (1), Assis (6); a direct orbit integration gives 0.3333, 0.1667, 0.9998.
        assert!(close(weber_perihelion_ratio(2.0), 1.0 / 3.0, 1e-12) && close(weber_perihelion_ratio(6.0), 1.0, 1e-12));
        assert!(close(weber_mass_shift(E_CHARGE, 1e4, M_E), 6.52317e-3, 1e-5));
    }
}
