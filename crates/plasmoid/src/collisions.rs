//! Collisions in a plasma: Spitzer's resistivity, which the pinch and reconnection calculators
//! share.
//!
//! The NRL Plasma Formulary gives the resistivity across the field as
//! `η⊥ = 1.03×10⁻⁴ Z lnΛ T^(−3/2)` Ω m, with `T` the electron temperature in eV, and along the
//! field `η∥ = η⊥/1.96` for singly charged ions (Spitzer & Härm 1953). The coefficient is
//! Braginskii's `m_e/(n_e e² τ_e)` rounded: 1.0314×10⁻⁴ from the constants.
//!
//! ```
//! use plasmoid::collisions::{Resistivity, spitzer_eta};
//!
//! // A 1 keV hydrogen plasma with lnΛ = 10 conducts better than copper (1.7×10⁻⁸ Ω m)
//! // only along the field.
//! let across = spitzer_eta(1000.0, 1.0, 10.0, Resistivity::Perpendicular);
//! let along = spitzer_eta(1000.0, 1.0, 10.0, Resistivity::Parallel);
//! assert!((across / 3.257e-8 - 1.0).abs() < 1e-3 && along < across);
//! ```

/// The coefficient of Spitzer's resistivity across the field (Ω m, for `T` in eV), from the NRL
/// Plasma Formulary.
pub const SPITZER: f64 = 1.03e-4;

/// `η⊥/η∥` for singly charged ions (Spitzer & Härm 1953, via the NRL Plasma Formulary). It rises
/// with the ion charge, toward the Lorentz-gas limit `32/(3π)` = 3.40.
pub const PARALLEL_RATIO: f64 = 1.96;

/// Which Spitzer resistivity a calculation uses.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Resistivity {
    /// Across the field, `η⊥ = 1.03×10⁻⁴ Z lnΛ T^(−3/2)` Ω m.
    Perpendicular,
    /// Along the field, `η⊥/1.96`.
    Parallel,
}

/// Spitzer's resistivity (Ω m) at electron temperature `te` (eV), ion charge `z` and Coulomb
/// logarithm `ln_lambda`, across or along the field.
pub fn spitzer_eta(te: f64, z: f64, ln_lambda: f64, resistivity: Resistivity) -> f64 {
    let across = SPITZER * z * ln_lambda * te.powf(-1.5);
    match resistivity {
        Resistivity::Perpendicular => across,
        Resistivity::Parallel => across / PARALLEL_RATIO,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::f64::consts::PI;
    use metrology::{E_CHARGE as E, EPS_0, M_E};

    #[test]
    fn the_coefficient_is_braginskii_rounded() {
        // η⊥ = m_e/(n e² τ_e) with τ_e = 3√m_e T^(3/2)(4πε₀)²/(4√(2π) n Z e⁴ lnΛ), at 1 eV,
        // Z = 1, lnΛ = 1: 1.031 36×10⁻⁴ Ω m, which NRL rounds to 1.03×10⁻⁴.
        let braginskii = 4.0 * (2.0 * PI).sqrt() * M_E.sqrt() * E * E / (3.0 * (4.0 * PI * EPS_0).powi(2) * E.powf(1.5));
        assert!((braginskii / 1.031_362e-4 - 1.0).abs() < 1e-5, "{braginskii:e}");
        assert!((braginskii / spitzer_eta(1.0, 1.0, 1.0, Resistivity::Perpendicular) - 1.0).abs() < 2e-3);
    }

    #[test]
    fn scalings() {
        let at = |te: f64, z: f64, r| spitzer_eta(te, z, 10.0, r);
        assert!((at(1000.0, 1.0, Resistivity::Perpendicular) / 3.257_146e-8 - 1.0).abs() < 1e-6);
        // Hotter is cleaner as T^(−3/2); more charge is dirtier in proportion.
        assert!((at(4000.0, 1.0, Resistivity::Perpendicular) / at(1000.0, 1.0, Resistivity::Perpendicular) - 0.125).abs() < 1e-12);
        assert!((at(1000.0, 2.0, Resistivity::Perpendicular) / at(1000.0, 1.0, Resistivity::Perpendicular) - 2.0).abs() < 1e-12);
        assert!((at(1000.0, 1.0, Resistivity::Perpendicular) / at(1000.0, 1.0, Resistivity::Parallel) - 1.96).abs() < 1e-12);
    }
}
