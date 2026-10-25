//! Tides in a neutron-star inspiral: how much each star deforms in its companion's field, and the
//! one combination of the two that the waveform measures first (at 5PN order).
//!
//! A star of mass `M`, radius `R` and quadrupolar Love number `k₂` has the dimensionless tidal
//! deformability `Λ = (2/3)k₂(c²R/GM)⁵`. A binary's phase carries the mass-weighted
//! `Λ̃ = (16/13)[(m₁ + 12m₂)m₁⁴Λ₁ + (m₂ + 12m₁)m₂⁴Λ₂]/(m₁ + m₂)⁵`, which equals both when the masses
//! are equal ([LVC, PRX 9, 011001](https://doi.org/10.1103/PhysRevX.9.011001), Eq. 5).
//!
//! ```
//! use gravitational::tides::{GW170817_LOW_SPIN, lambda, lambda_tilde};
//! use metrology::SOLAR_MASS;
//!
//! // Two 1.36 M☉ stars of 12 km with k₂ = 0.1: Λ̃ = 508, inside GW170817's 70–720.
//! let l = lambda(0.1, 12e3, 1.36 * SOLAR_MASS);
//! let tilde = lambda_tilde(1.36, 1.36, l, l);
//! assert!((tilde - l).abs() < 1e-9 && (tilde - 508.0).abs() < 1.0);
//! assert!(tilde > GW170817_LOW_SPIN.0 && tilde < GW170817_LOW_SPIN.1);
//! ```

use metrology::{C, G};

/// GW170817's `Λ̃` with spins below 0.05: the 90% highest-posterior-density interval
/// 300 (+420/−230) (LVC, PRX 9, 011001, Table II, PhenomPNRT, flat prior in `Λ̃`).
pub const GW170817_LOW_SPIN: (f64, f64) = (70.0, 720.0);

/// GW170817's 90% upper limit on `Λ̃` with spins up to 0.89 (the same table).
pub const GW170817_HIGH_SPIN_MAX: f64 = 630.0;

/// The tidal deformability `Λ = (2/3)k₂(c²R/GM)⁵` of a star with Love number `k2`, radius `radius`
/// (m) and mass `mass` (kg).
pub fn lambda(k2: f64, radius: f64, mass: f64) -> f64 {
    2.0 / 3.0 * k2 * (C * C * radius / (G * mass)).powi(5)
}

/// The binary tidal deformability `Λ̃` of stars of masses `m1`, `m2` (any one unit) with
/// deformabilities `l1`, `l2`.
pub fn lambda_tilde(m1: f64, m2: f64, l1: f64, l2: f64) -> f64 {
    16.0 / 13.0 * ((m1 + 12.0 * m2) * m1.powi(4) * l1 + (m2 + 12.0 * m1) * m2.powi(4) * l2) / (m1 + m2).powi(5)
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::f64::consts::PI;
    use metrology::SOLAR_MASS;

    #[test]
    fn the_injections_of_table_vi() {
        // LVC's injections: (m₁, m₂), (Λ₁, Λ₂) → Λ̃, printed to the nearest unit.
        for (m1, m2, l1, l2, want) in [(1.38, 1.37, 275.0, 309.0, 292.0), (1.68, 1.13, 77.0, 973.0, 303.0), (1.38, 1.37, 1018.0, 1063.0, 1040.0)] {
            let got = lambda_tilde(m1, m2, l1, l2);
            assert!((got - want).abs() < 0.5, "({m1}, {m2}): {got}");
        }
        // Equal masses give Λ̃ = Λ, and the order of the stars does not matter.
        assert!((lambda_tilde(1.4, 1.4, 190.0, 190.0) - 190.0).abs() < 1e-12);
        assert!((lambda_tilde(1.4, 1.2, 190.0, 500.0) - lambda_tilde(1.2, 1.4, 500.0, 190.0)).abs() < 1e-12);
    }

    #[test]
    fn a_star_deforms_as_its_compactness_to_the_minus_fifth() {
        // A Newtonian n = 1 polytrope has k₂ = (15 − π²)/(2π²) = 0.25991; integrating the
        // Clairaut–Radau equation gives the same to 2×10⁻¹² (and 3/4 for uniform density). At
        // compactness C = GM/(c²R) its Λ is (2/3)k₂/C⁵, a hundred times larger at a C 2.5× smaller.
        let k2 = (15.0 - PI * PI) / (2.0 * PI * PI);
        assert!((k2 - 0.25991).abs() < 5e-6);
        let (m, r) = (1.4 * SOLAR_MASS, 12e3);
        let c = G * m / (C * C * r);
        assert!((lambda(k2, r, m) * c.powi(5) - 2.0 / 3.0 * k2).abs() < 1e-15);
        assert!((lambda(k2, 2.5 * r, m) / lambda(k2, r, m) - 2.5f64.powi(5)).abs() < 1e-9);
    }
}
