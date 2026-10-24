//! Laboratory gravity, sized honestly: the gravitational field of a light or particle beam, the
//! gravitational waves a spinning body radiates, and the energy a wave carries.
//!
//! A beam of light or matter really does carry gravity. Parallel beams do not attract each other
//! (Tolman, Ehrenfest & Podolsky 1931), but a test mass at rest near a beam of power P is pulled
//! by `a = 4GP/(c³ρ)` while the beam passes ([Rätzel, Wilkens & Menzel 2016](https://arxiv.org/abs/1511.01023)).
//! This is a near field confined to the beam's surroundings, not a projected ray. Radiated waves
//! from lab rotors are suppressed by G/c⁵ (the quadrupole formula).
//!
//! ```
//! use gwsens::lab;
//!
//! // A petawatt laser, felt 1 mm from the beam:
//! let a = lab::beam_acceleration(1e15, 1e-3);
//! assert!((a / 9.91e-18 - 1.0).abs() < 0.01); // 10⁻¹⁷ m/s²: no instrument can see it
//! ```

use crate::{C, G};
use core::f64::consts::PI;

/// Acceleration of a test mass at rest at distance `rho` (m) from a beam of power `power` (W),
/// while the beam passes: `a = 4GP/(c³ρ)` (twice the Newtonian line-mass value).
pub fn beam_acceleration(power: f64, rho: f64) -> f64 {
    4.0 * G * power / (C.powi(3) * rho)
}

/// Tidal acceleration gradient next to the beam, s⁻²: `−4GP/(c³ρ²)`. This is
/// [Spengler, Rätzel & Braun 2022](https://arxiv.org/abs/2104.09209)'s Eq. 1, `R₀ₓ₀ₓ ≃ −4GP/(c⁵ρ²)`
/// in m⁻² (with x⁰ = ct), times c².
pub fn beam_tidal(power: f64, rho: f64) -> f64 {
    -4.0 * G * power / (C.powi(3) * rho * rho)
}

/// Velocity kick from a pulse of energy `energy` (J) passing at distance `rho` (m): the
/// acceleration integrated over the pulse duration, `Δv = 4GE/(c³ρ)`.
pub fn beam_velocity_kick(energy: f64, rho: f64) -> f64 {
    4.0 * G * energy / (C.powi(3) * rho)
}

/// Gravitational-wave power from a rod of mass `m` (kg) and length `l` (m) spinning about its
/// middle at `omega` (rad/s): `P = (2/45) G m² l⁴ ω⁶/c⁵`, W.
pub fn spinning_rod_power(m: f64, l: f64, omega: f64) -> f64 {
    2.0 / 45.0 * G * m * m * l.powi(4) * omega.powi(6) / C.powi(5)
}

/// Gravitational-wave power from a rigid body spinning at `omega` about a principal axis with
/// the other two moments of inertia differing by `delta_i` (kg m²):
/// `P = (32/5)(G/c⁵)(ΔI)² ω⁶`, W.
pub fn rotor_power(delta_i: f64, omega: f64) -> f64 {
    32.0 / 5.0 * G / C.powi(5) * delta_i * delta_i * omega.powi(6)
}

/// Wavelength of the waves from a body spinning at `omega`: they come out at twice the spin
/// frequency, so `λ = πc/ω`, m.
pub fn rotor_wavelength(omega: f64) -> f64 {
    PI * C / omega
}

/// Energy flux of a plane wave of amplitude `h` at frequency `f`: `F = c³ω²h²/(32πG)`, W/m².
pub fn wave_flux(h: f64, f: f64) -> f64 {
    let w = 2.0 * PI * f;
    C.powi(3) * w * w * h * h / (32.0 * PI * G)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn beam_fixtures() {
        assert!((beam_acceleration(1e15, 1e-3) / 9.91e-18 - 1.0).abs() < 2e-3);
        // The LHC beam, 3.8×10¹² W, at 100 µm.
        assert!((beam_acceleration(3.8e12, 1e-4) / 3.77e-19 - 1.0).abs() < 2e-3);
        // A 100 J pulse at 1 mm: about 10⁻³⁰ m/s.
        assert!((beam_velocity_kick(100.0, 1e-3) / 9.9e-31 - 1.0).abs() < 0.01);
        // 1 PW at 1 mm: −1.10×10⁻³¹ m⁻², i.e. −9.91×10⁻¹⁵ s⁻²; the gradient of the acceleration.
        assert!((beam_tidal(1e15, 1e-3) / -9.91e-15 - 1.0).abs() < 2e-3);
        assert!((beam_tidal(1e15, 1e-3) * 1e-3 / -beam_acceleration(1e15, 1e-3) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn misner_thorne_wheeler_rod() {
        // 490 t, 20 m, spun to break-up at 28 rad/s: 2.27×10⁻²⁹ W at λ ≈ 3.36×10⁷ m.
        let p = spinning_rod_power(4.9e5, 20.0, 28.0);
        assert!((p / 2.27e-29 - 1.0).abs() < 0.01, "{p:e}");
        assert!((rotor_wavelength(28.0) / 3.36e7 - 1.0).abs() < 0.01);
        // The rod is a rotor with ΔI = m l²/12.
        let i = 4.9e5 * 400.0 / 12.0;
        assert!((rotor_power(i, 28.0) / p - 1.0).abs() < 1e-12);
    }

    #[test]
    fn flux_scales_as_h_squared() {
        assert!((wave_flux(1e-21, 100.0) / 1.5853e-3 - 1.0).abs() < 1e-4);
        assert!((wave_flux(2e-21, 100.0) / wave_flux(1e-21, 100.0) - 4.0).abs() < 1e-12);
    }
}
