//! Gravity and superconductors: Tajmar's gravitomagnetic London moment, Podkletnov's shielding
//! and impulse generator, and what general relativity predicts instead.
//!
//! **Status: refuted.**
//! - **Tajmar & de Matos (2006–2007)** reported frame dragging near spinning niobium rings of
//!   k ≈ 2.6×10⁻⁵ of the ring's rotation. General relativity gives `k = 4GM/(c²(R_o + R_i))`
//!   ≈ 7×10⁻²⁷, so the claim is 3.6×10²¹ times larger. Tajmar's own later work set a limit
//!   800× below the earlier signal and traced it to vibration.
//! - **Podkletnov (1997)** reported 0.3–2.1% weight loss above a spinning superconducting disc.
//!   Replications bound any effect below 2×10⁻⁸ (NASA Marshall), 10⁶ times smaller. A shield
//!   would have to cast a shadow that shrinks with height as `1 − cos³θ`; he reported none.
//! - **Podkletnov's impulse generator** reproduces from pendulum geometry alone, but the energy
//!   it implies needs 40× the charge his own Marx bank holds, and the swing's momentum is 10⁹
//!   times too large for any radiation of that energy.
//!
//! The legitimate physics is the London moment, `B = 2m_eω/e`, measured to ppm precision.
//!
//! ```
//! use gravity_claims::superconductors::{TAJMAR_CLAIM, ring_frame_dragging};
//!
//! // Tajmar's 0.349 kg niobium ring (radii 75 and 69 mm).
//! let gr = ring_frame_dragging(0.349, 0.075, 0.069);
//! println!("GR {gr:.1e}, claim {TAJMAR_CLAIM:.1e}, ratio {:.1e}", TAJMAR_CLAIM / gr); // 3.6e21
//! ```

use core::f64::consts::PI;

use crate::{C, E_CHARGE, G, G0, M_E, MU_0};

/// The frame-dragging ratio k (frame rotation over ring rotation) Tajmar & de Matos reported,
/// from a coupling of −9.64×10⁻⁸ s² at r = 3.6 cm.
pub const TAJMAR_CLAIM: f64 = 2.63e-5;

/// Podkletnov's largest reported weight loss, as a fraction (1997, while braking).
pub const PODKLETNOV_CLAIM: f64 = 2.1e-2;

/// Upper limits on the fractional weight change above a spinning superconductor.
pub const WEIGHT_BOUNDS: [(&str, f64); 4] = [
    ("Li et al. 1997, NASA Marshall", 2e-8),
    ("Tajmar et al. 2006", 5e-6),
    ("Hathaway et al. 2003", 1e-5),
    ("Woods et al. 2001", 4e-5),
];

/// The London moment: a spinning superconductor at `omega` (rad/s) carries a field
/// `B = 2m_eω/e`, T.
pub fn london_moment(omega: f64) -> f64 {
    2.0 * M_E * omega / E_CHARGE
}

/// General relativity's frame-dragging ratio inside a ring of `mass` (kg) and radii `r_out`,
/// `r_in` (m): `k = 4GM/(c²(R_o + R_i))`.
pub fn ring_frame_dragging(mass: f64, r_out: f64, r_in: f64) -> f64 {
    4.0 * G * mass / (C * C * (r_out + r_in))
}

/// The Cooper-pair mass density `ρ* = m_e²/(μ₀e²λ_L²)` (kg/m³) of a superconductor with London
/// depth `london_depth` (m). Tajmar & de Matos predicted `k = ρ*/ρ`.
pub fn cooper_pair_density(london_depth: f64) -> f64 {
    M_E * M_E / (MU_0 * E_CHARGE * E_CHARGE * london_depth * london_depth)
}

/// The gravitational analogue of the London depth, `c/√(4πGρ*)`, m: how far a "gravitational
/// Meissner effect" would reach. About 10¹⁴ m, so there is none in the lab.
pub fn gravitational_london_depth(rho_star: f64) -> f64 {
    C / (4.0 * PI * G * rho_star).sqrt()
}

/// The fraction of a uniform Earth's pull a disc of `radius` would block for a sample `height`
/// above it, if it were a perfect shield: `1 − cos³θ`, `θ = atan(radius/height)`.
pub fn shadow_fraction(radius: f64, height: f64) -> f64 {
    1.0 - (radius / height).atan().cos().powi(3)
}

/// A pendulum's response to a short push, from its swing.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Swing {
    /// Height gained, m.
    pub rise: f64,
    /// Energy given, J.
    pub energy: f64,
    /// Initial speed, m/s.
    pub speed: f64,
    /// Mean acceleration during the push, m/s².
    pub acceleration: f64,
}

/// What a pendulum of `length` (m) and `mass` (kg) swung `deflection` (m) sideways by a push
/// lasting `pulse` (s) received.
pub fn pendulum_impulse(length: f64, deflection: f64, mass: f64, pulse: f64) -> Swing {
    let rise = length - (length * length - deflection * deflection).sqrt();
    let speed = (2.0 * G0 * rise).sqrt();
    Swing { rise, energy: mass * G0 * rise, speed, acceleration: speed / pulse }
}

/// Energy (J) and charge (C) of a Marx generator of `stages` stages, each of `capacitance` (F)
/// charged to `voltage` (V).
pub fn marx_bank(capacitance: f64, stages: f64, voltage: f64) -> (f64, f64) {
    (stages * 0.5 * capacitance * voltage * voltage, capacitance * voltage)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, rel: f64) -> bool {
        (a / b - 1.0).abs() < rel
    }

    #[test]
    fn tajmar_against_general_relativity() {
        let gr = ring_frame_dragging(0.349, 0.075, 0.069);
        assert!(close(gr, 7.20e-27, 2e-3));
        assert!(close(TAJMAR_CLAIM / gr, 3.65e21, 5e-3));
        // Their predicted ρ*/ρ for niobium with λ_L = 27.6 nm.
        assert!(close(cooper_pair_density(27.6e-9) / 8570.0, 3.94e-6, 2e-3));
        // With λ_L = 39 nm the gravitational London depth is 8×10¹³ m.
        assert!(close(gravitational_london_depth(cooper_pair_density(39e-9)), 7.96e13, 5e-3));
        // The London moment, 1.137×10⁻¹¹ T s: 7.39 nT at 650 rad/s.
        assert!(close(london_moment(650.0), 7.39e-9, 2e-3));
    }

    #[test]
    fn podkletnov_against_replications_and_geometry() {
        let strictest = WEIGHT_BOUNDS.iter().map(|b| b.1).fold(f64::INFINITY, f64::min);
        assert!(close(PODKLETNOV_CLAIM / strictest, 1.05e6, 1e-9));
        // A perfect 275 mm shield blocks 0.31% of the pull 3 m above it, all of it right above.
        assert!(close(shadow_fraction(0.1375, 3.0), 3.14e-3, 2e-3));
        assert!((shadow_fraction(0.1375, 1e-9) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn impulse_generator_tables_reproduce() {
        // 18.5 g pendulum, 0.8 m long, swung 142 mm by a 0.1 ms pulse at 2 MV.
        let s = pendulum_impulse(0.8, 0.142, 0.0185, 1e-4);
        assert!(close(s.rise, 12.70e-3, 1e-3));
        assert!(close(s.energy, 2.305e-3, 1e-3));
        assert!(close(s.speed, 0.499, 1e-3));
        assert!(close(s.acceleration, 4.99e3, 1e-3)); // ~509 g, published "~10³ g"
        // Radiation carrying that energy would carry momentum E/c; the swing has 10⁹ times more.
        assert!(close(0.0185 * s.speed * C / s.energy, 1.2e9, 0.01));
        // The Marx bank: 20 stages of 25 nF at 100 kV hold 2.5 kJ and 2.5 mC, not 0.1 C.
        let (energy, charge) = marx_bank(25e-9, 20.0, 1e5);
        assert!(close(energy, 2500.0, 1e-12) && close(0.1 / charge, 40.0, 1e-12));
    }
}
