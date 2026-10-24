//! Warp drives, wormholes, negative energy and Puthoff's polarizable vacuum: the closed forms
//! behind the DIA "DIRD" papers on spacetime engineering, from the mainstream literature.
//!
//! The formulas are standard general relativity and quantum field theory. A warp bubble needs a
//! large mass of negative energy, while quantum field theory allows negative energy only in tiny
//! amounts for short times (the quantum inequalities). Puthoff's polarizable-vacuum metric is a
//! genuine alternative: it matches general relativity at first order and differs at second order,
//! by 0.73 µas for light grazing the Sun, so it is testable in principle.
//!
//! ```
//! use gravitics::spacetime::alcubierre_energy;
//!
//! // A 100 m warp bubble at light speed with a 1 m wall needs half a Sun of negative energy.
//! let sun = 1.989e30;
//! println!("{:.2} solar masses", alcubierre_energy(100.0, 1.0, 1.0) / sun); // -0.56
//! ```

use core::f64::consts::PI;

use crate::{C, G, HBAR, K_B};

/// Total energy of an Alcubierre warp bubble of `radius` (m) and wall thickness `wall` (m) moving
/// at `speed` (in units of c), as a mass (kg): `−(v²/12)(R²/Δ)·c²/G`
/// ([Pfenning & Ford 1997](https://arxiv.org/abs/gr-qc/9702026), Eq. 28). Negative. (Their text's
/// worked numbers, "6.2×10⁶⁵ g" and "a quarter of a solar mass", disagree with their Eq. 28,
/// which gives 6.9×10⁶² kg and 0.56 solar masses for the same inputs.)
pub fn alcubierre_energy(radius: f64, wall: f64, speed: f64) -> f64 {
    -speed * speed / 12.0 * radius * radius / wall * C * C / G
}

/// The most negative energy density (J/m³) quantum field theory allows when averaged over a time
/// `tau` (s): `−3ħc/(32π²(cτ)⁴)` (the quantum inequality, Pfenning & Ford 1997 eq. 9).
pub fn quantum_inequality(tau: f64) -> f64 {
    -3.0 * HBAR * C / (32.0 * PI * PI * (C * tau).powi(4))
}

/// Casimir energy density (J/m³) between ideal plates a distance `d` (m) apart: `−π²ħc/(720d⁴)`.
pub fn casimir_energy_density(d: f64) -> f64 {
    -PI * PI * HBAR * C / (720.0 * d.powi(4))
}

/// Casimir pressure (Pa) between ideal plates `d` (m) apart: `−π²ħc/(240d⁴)` (attractive).
pub fn casimir_pressure(d: f64) -> f64 {
    -PI * PI * HBAR * C / (240.0 * d.powi(4))
}

/// Casimir force (N) between a sphere of `radius` (m) and a plate `d` (m) away, in the proximity
/// force approximation: `−π³ħcR/(360d³)`.
pub fn casimir_sphere_plate(radius: f64, d: f64) -> f64 {
    -PI.powi(3) * HBAR * C * radius / (360.0 * d.powi(3))
}

/// Radial tension (Pa) at the throat of a Morris–Thorne wormhole of throat radius `b0` (m):
/// `c⁴/(8πG b₀²)`. Holding the throat open needs this much negative energy.
pub fn wormhole_throat_tension(b0: f64) -> f64 {
    C.powi(4) / (8.0 * PI * G * b0 * b0)
}

/// Line tension (J/m) of the exotic edges of Visser's cubical wormhole: `−c⁴/(8G)`.
pub fn wormhole_edge_tension() -> f64 {
    -C.powi(4) / (8.0 * G)
}

/// Temperature (K) of a horizon of surface gravity `kappa` (s⁻¹): `ħκ/(2πk_B)`. With κ = a/c it is
/// the Unruh temperature, with c³/(4GM) Hawking's, with the Hubble rate Gibbons–Hawking's.
pub fn horizon_temperature(kappa: f64) -> f64 {
    HBAR * kappa / (2.0 * PI * K_B)
}

/// Unruh temperature (K) seen by an observer accelerating at `a` (m/s²).
pub fn unruh_temperature(a: f64) -> f64 {
    horizon_temperature(a / C)
}

/// Hawking temperature (K) of a black hole of mass `m` (kg).
pub fn hawking_temperature(m: f64) -> f64 {
    horizon_temperature(C.powi(3) / (4.0 * G * m))
}

/// Gibbons–Hawking temperature (K) of the cosmological horizon for Hubble rate `h` (s⁻¹).
pub fn gibbons_hawking_temperature(h: f64) -> f64 {
    horizon_temperature(h)
}

/// Evaporation time (s) of a black hole of mass `m` (kg) radiating photons only:
/// `5120πG²M³/(ħc⁴)`.
pub fn hawking_lifetime(m: f64) -> f64 {
    5120.0 * PI * G * G * m.powi(3) / (HBAR * C.powi(4))
}

/// Zero-point energy density (J/m³) of the electromagnetic field up to a cut-off `omega_c`
/// (rad/s): `ħω_c⁴/(8π²c³)`.
pub fn zero_point_energy_density(omega_c: f64) -> f64 {
    HBAR * omega_c.powi(4) / (8.0 * PI * PI * C.powi(3))
}

/// A static, spherical metric.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Metric {
    /// Schwarzschild.
    GeneralRelativity,
    /// Puthoff's polarizable vacuum, `K = exp(2GM/rc²)`
    /// ([Puthoff 2002](https://arxiv.org/abs/gr-qc/9909037)).
    PolarizableVacuum,
}

/// Deflection (rad) of light passing a `mass` (kg) at impact parameter `b` (m), to second order:
/// `4m/b + c₂(m/b)²` with `m = GM/c²`, where `c₂` is `15π/4` in general relativity and `4π` in the
/// polarizable vacuum.
pub fn light_bending(mass: f64, b: f64, metric: Metric) -> f64 {
    let x = G * mass / (C * C * b);
    let second = match metric {
        Metric::GeneralRelativity => 15.0 * PI / 4.0,
        Metric::PolarizableVacuum => 4.0 * PI,
    };
    4.0 * x + second * x * x
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::PLANCK_LENGTH;

    fn close(a: f64, b: f64, rel: f64) -> bool {
        (a / b - 1.0).abs() < rel
    }

    #[test]
    fn warp_bubbles_and_quantum_inequalities() {
        // R = 100 m with a wall of 100 Planck lengths (v = c): 6.9×10⁶² kg (the paper: 6.2×10⁶²).
        assert!(close(alcubierre_energy(100.0, 100.0 * PLANCK_LENGTH, 1.0), -6.94e62, 5e-3));
        // A 1 m wall: 0.56 solar masses.
        assert!(close(alcubierre_energy(100.0, 1.0, 1.0) / 1.989e30, -0.564, 5e-3));
        assert!(close(quantum_inequality(1e-9), -3.72e-26, 5e-3));
        assert!(close(quantum_inequality(1e-15), -3.72e-2, 5e-3));
    }

    #[test]
    fn casimir_closed_forms() {
        assert!(close(casimir_pressure(1e-6), -1.300e-3, 1e-3));
        // Pressure is −dE/dd of the energy per area u·d, and the sphere force is 2πR·u(d)·d.
        let d = 1e-6;
        let energy = |d: f64| casimir_energy_density(d) * d;
        let h = 1e-12;
        assert!(close(-(energy(d + h) - energy(d - h)) / (2.0 * h), casimir_pressure(d), 1e-6));
        assert!(close(casimir_sphere_plate(1e-4, d), 2.0 * PI * 1e-4 * energy(d), 1e-12));
    }

    #[test]
    fn wormholes() {
        assert!(close(wormhole_throat_tension(10.0), 4.82e40, 5e-3));
        assert!(close(wormhole_edge_tension(), -1.513e43, 1e-3));
    }

    #[test]
    fn horizons() {
        // Unruh: 4.055×10⁻²¹ K per m/s²; 1 K needs 2.466×10²⁰ m/s².
        assert!(close(unruh_temperature(1.0), 4.055e-21, 1e-3));
        assert!(close(1.0 / unruh_temperature(1.0), 2.466e20, 1e-3));
        // A solar-mass black hole: 61.7 nK, evaporating (photons only) in 2.1×10⁶⁷ years.
        let sun = 1.327_124_4e20 / G;
        assert!(close(hawking_temperature(sun), 6.17e-8, 2e-3));
        assert!(close(hawking_lifetime(sun) / 3.155_76e7, 2.095e67, 3e-3));
        // The cosmological horizon at H₀ = 67.4 km/s/Mpc: 2.7×10⁻³⁰ K.
        assert!(close(gibbons_hawking_temperature(67.4e3 / 3.085_677_581e22), 2.655e-30, 2e-3));
        // Zero-point energy to the electron's Compton frequency: 1.8×10²² J/m³.
        let compton = 9.109_383_701_5e-31 * C * C / HBAR;
        assert!(close(zero_point_energy_density(compton), 1.801e22, 2e-3));
    }

    #[test]
    fn schwinger_is_a_boltzmann_factor_at_twice_the_unruh_temperature() {
        // πE_S/E = m_ec²/(2k_B T_U) for an electron accelerated by E: exactly, at any field.
        use crate::strong_field::{E_SCHWINGER, schwinger_exponent};
        let (m_e, e) = (9.109_383_701_5e-31, 1.602_176_634e-19);
        for f in [1.0, 0.1, 1e-3] {
            let field = f * E_SCHWINGER;
            let t_u = unruh_temperature(e * field / m_e);
            assert!(close(schwinger_exponent(field), m_e * C * C / (2.0 * K_B * t_u), 1e-12));
        }
        // At the critical field itself, T_U = 9.44×10⁸ K (81.3 keV).
        assert!(close(unruh_temperature(e * E_SCHWINGER / m_e), 9.438e8, 1e-3));
    }

    #[test]
    fn polarizable_vacuum_against_general_relativity() {
        let (sun, limb) = (1.327_124_4e20 / G, 6.957e8);
        let arcsec = 180.0 / PI * 3600.0;
        let gr = light_bending(sun, limb, Metric::GeneralRelativity);
        assert!((gr * arcsec - 1.7512).abs() < 1e-4);
        // The two agree at first order and differ by (π/4)(m/b)² = 0.73 µas.
        let diff = light_bending(sun, limb, Metric::PolarizableVacuum) - gr;
        assert!(close(diff * arcsec * 1e6, 0.730, 5e-3));
    }
}
