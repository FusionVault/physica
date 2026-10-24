//! Claims that a pulsed high-voltage or time-varying magnetic field makes a local gravitational
//! field, computed against the two classical forces that act in the same apparatus: the
//! dielectrophoretic (Kelvin) polarisation force and the angular momentum stored in crossed
//! electric and magnetic fields.
//!
//! The claim is in Pais's US Navy patents and, more recently, in a 2025 Chinese patent and in
//! replications by independent experimenters: fields changing fast enough should produce a
//! gravito-electric or "vortex gravitational" field. Two setups recur — a dielectric sheet
//! between high-voltage electrodes, and a suspended dielectric ball between offset coils in
//! vacuum.
//!
//! Both setups have a conventional explanation whose size this module computes. The reported
//! signature of the first, that the sheet deflects the same way when the polarity is reversed,
//! is exactly what [`dielectrophoretic_force`] predicts, because that force follows `∇E²` and so
//! does not care about the field's sign. The second setup stores angular momentum in its crossed
//! fields ([`field_angular_momentum`]), which is delivered as a torque when the field is switched.
//!
//! What would separate a real effect from either is [`mass_independence`]: the classical forces
//! scale with a body's volume and permittivity, while a gravitational field accelerates every
//! body alike. Swapping the test body for a denser one of the same shape therefore changes the
//! acceleration if the force is electrostatic, and leaves it unchanged if the field is
//! gravitational.
//!
//! ```
//! use gravitics::electrogravity::{dielectrophoretic_force, mass_independence};
//!
//! // A 2 cm square of 0.15 mm plastic (ε_r = 2.2) near a 40 kV electrode, where
//! // E ≈ 10⁷ V/m falls off over centimetres: ∇E² ≈ 2×10¹⁶ V²/m³.
//! let volume = 0.02 * 0.02 * 0.15e-3;
//! let f = dielectrophoretic_force(volume, 2.2, 2e16);
//! assert!((f * 1e3 - 6.4).abs() < 0.1); // 6.4 mN: the weight of 650 mg
//!
//! // Polyethylene against tungsten, same shape: a classical force gives 20× the acceleration
//! // on the lighter body, a gravitational field exactly the same on both.
//! let ratio = mass_independence(950.0, 19_300.0);
//! assert!((ratio - 20.3).abs() < 0.1);
//! ```

use metrology::EPS_0;

/// The dielectrophoretic (Kelvin) force (N) on a neutral dielectric of volume `volume` (m³) and
/// relative permittivity `eps_r` in a field whose squared magnitude has gradient `grad_e2`
/// (V²/m³): `F = ½Vε₀(ε_r − 1)∇(E²)`.
///
/// It points toward the strongest field, which is the sharpest electrode. Because it depends on
/// `E²`, **reversing the polarity leaves it unchanged** — the signature often reported as
/// evidence of a gravitational effect.
pub fn dielectrophoretic_force(volume: f64, eps_r: f64, grad_e2: f64) -> f64 {
    0.5 * volume * EPS_0 * (eps_r - 1.0) * grad_e2
}

/// The angular momentum (kg m²/s) stored in crossed static fields over a volume `volume` (m³):
/// `L = ε₀ r (E × B) V`, for field magnitudes `e` (V/m) and `b` (T) at lever arm `lever` (m),
/// taken perpendicular.
///
/// Switching the fields off transfers this to the apparatus as a mechanical impulse: the Feynman
/// disc paradox. A suspended body in such a field therefore turns when the field changes, with no
/// new physics.
pub fn field_angular_momentum(e: f64, b: f64, lever: f64, volume: f64) -> f64 {
    EPS_0 * lever * e * b * volume
}

/// How much more a classical polarisation force accelerates a light test body than a dense one of
/// the same shape: the density ratio `rho_dense/rho_light`.
///
/// A force set by volume and permittivity is the same on both, so the acceleration goes as `1/ρ`.
/// A gravitational field gives `a = E_g` whatever the body is made of, so the ratio is 1. Running
/// the same experiment with two bodies of equal shape and different density therefore separates
/// the two explanations, and it needs no absolute calibration.
pub fn mass_independence(rho_light: f64, rho_dense: f64) -> f64 {
    rho_dense / rho_light
}

/// The acceleration (m/s²) a force `force` (N) gives a body of volume `volume` (m³) and density
/// `rho` (kg/m³): `F/(ρV)`. Pair it with [`dielectrophoretic_force`] to turn a claimed
/// displacement into the acceleration a gravimeter would have to resolve.
pub fn acceleration(force: f64, volume: f64, rho: f64) -> f64 {
    force / (rho * volume)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, rel: f64) -> bool {
        (a / b - 1.0).abs() < rel
    }

    #[test]
    fn the_sheet_setup_is_dielectrophoresis() {
        let volume = 0.02 * 0.02 * 0.15e-3;
        let f = dielectrophoretic_force(volume, 2.2, 2e16);
        assert!(close(f, 6.375e-3, 1e-3), "{f:e}");
        // Reversing the polarity does not change it: E -> -E leaves grad(E^2) alone. This is the
        // reported behaviour of the claim's own experiment.
        assert!(dielectrophoretic_force(volume, 2.2, 2e16) == f);
        // It vanishes in a uniform field and for vacuum-like material.
        assert!(dielectrophoretic_force(volume, 2.2, 0.0) == 0.0);
        assert!(dielectrophoretic_force(volume, 1.0, 2e16) == 0.0);
        // 6.4 mN lifts 650 mg, far more than any claimed gravitational effect.
        assert!(close(f / metrology::G0 * 1e6, 650.0, 1e-2));
    }

    #[test]
    fn the_ball_setup_stores_field_angular_momentum() {
        // 10^6 V/m across 10^-3 T over 10^-4 m^3 at a 5 cm lever arm.
        let l = field_angular_momentum(1e6, 1e-3, 0.05, 1e-4);
        assert!(close(l, 4.427e-14, 1e-3), "{l:e}");
        // It scales with both fields, so it disappears if either is removed.
        assert!(field_angular_momentum(0.0, 1e-3, 0.05, 1e-4) == 0.0);
        assert!(close(field_angular_momentum(2e6, 1e-3, 0.05, 1e-4) / l, 2.0, 1e-12));
    }

    #[test]
    fn mass_independence_is_the_discriminator() {
        // Polyethylene (950) against tungsten (19 300) of the same shape.
        assert!(close(mass_independence(950.0, 19_300.0), 20.32, 1e-3));
        // The same force on both bodies: the acceleration falls by exactly that ratio.
        let (volume, force) = (1e-6, 1e-6);
        let light = acceleration(force, volume, 950.0);
        let dense = acceleration(force, volume, 19_300.0);
        assert!(close(light / dense, mass_independence(950.0, 19_300.0), 1e-12));
        // A gravitational field would give both the same acceleration, so the ratio would be 1.
        assert!(mass_independence(950.0, 950.0) == 1.0);
        // The polyethylene ball of the reported setup, if the claimed force were 1 mN on 1 g:
        assert!(close(acceleration(1e-3, 1e-6, 1000.0), 1.0, 1e-12));
    }
}
