//! Gravitational-wave memory: the permanent strain a burst leaves behind, from the waves' own
//! energy (the nonlinear Christodoulou memory), and how it shows in pulsar timing.
//!
//! The memory of a binary builds up over the inspiral and stays. At leading order it is
//! `(2ηMx/R) sin²Θ (17 + cos²Θ)/96` ([Favata 2010](https://arxiv.org/abs/1003.3486)), about a
//! twentieth of the peak strain. It has not been detected: LIGO needs to stack many events, and
//! a pulsar timing array sees it as a ramp in the residuals, a step in the pulsar's apparent
//! spin frequency.
//!
//! ```
//! use gravitational::memory::memory_0pn;
//! use metrology::{MPC, SOLAR_MASS};
//!
//! // GW150914 edge-on at 410 Mpc, evaluated at the innermost stable orbit (x = 1/6): 1.1×10⁻²².
//! let h = memory_0pn(36.0 * SOLAR_MASS, 29.0 * SOLAR_MASS, 410.0 * MPC, 1.0 / 6.0, core::f64::consts::FRAC_PI_2);
//! assert!((h / 1.10657e-22 - 1.0).abs() < 1e-5);
//! ```

use metrology::{C, G};

/// Leading-order (0PN) memory strain of a binary of masses `m1`, `m2` (kg) at distance `d` (m),
/// accumulated up to the post-Newtonian parameter `x` and seen at inclination `theta` (radians):
/// `(2ηMx/R) sin²Θ (17 + cos²Θ)/96`, with `M = G(m₁+m₂)/c²`. Zero face-on; a lower bound, as
/// the merger adds more.
pub fn memory_0pn(m1: f64, m2: f64, d: f64, x: f64, theta: f64) -> f64 {
    let m = G * (m1 + m2) / (C * C);
    let eta = m1 * m2 / (m1 + m2).powi(2);
    let (s, c) = theta.sin_cos();
    2.0 * eta * m * x / d * s * s * (17.0 + c * c) / 96.0
}

/// Timing residual (s) at time `t` of a pulsar a distance `l` (m) away from a memory burst of
/// strain `h0` that passed Earth at `t0`, for wave polarisation angle `psi` and angle `alpha`
/// between the pulsar and the wave's direction of travel: `B h₀[(t−t₀)Θ(t−t₀) − (t−tᵢ)Θ(t−tᵢ)]`,
/// with `B = ½ cos 2ψ (1 − cos α)` and the pulsar term at `tᵢ = t₀ + (l/c)(1 + cos α)`.
pub fn pta_memory_residual(t: f64, t0: f64, h0: f64, psi: f64, alpha: f64, l: f64) -> f64 {
    let b = 0.5 * (2.0 * psi).cos() * (1.0 - alpha.cos());
    let ramp = |since: f64| since.max(0.0);
    let ti = t0 + l / C * (1.0 + alpha.cos());
    b * h0 * (ramp(t - t0) - ramp(t - ti))
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::f64::consts::FRAC_PI_2;
    use metrology::{MPC, PARSEC, SOLAR_MASS, YEAR};

    #[test]
    fn binary_memory() {
        let h = memory_0pn(1e9 * SOLAR_MASS, 1e9 * SOLAR_MASS, 1e3 * MPC, 1.0 / 6.0, FRAC_PI_2);
        assert!((h / 1.41236e-15 - 1.0).abs() < 1e-5);
        assert_eq!(memory_0pn(36.0 * SOLAR_MASS, 29.0 * SOLAR_MASS, 410.0 * MPC, 1.0 / 6.0, 0.0), 0.0);
    }

    #[test]
    fn pulsar_ramp() {
        // B = ½, h₀ = 10⁻¹⁴, five years after the Earth term and before the pulsar term: 0.789 µs.
        let r = pta_memory_residual(5.0 * YEAR, 0.0, 1e-14, 0.0, FRAC_PI_2, 1e3 * PARSEC);
        assert!((r / 0.78894e-6 - 1.0).abs() < 1e-5);
        assert_eq!(pta_memory_residual(-1.0, 0.0, 1e-14, 0.0, FRAC_PI_2, 1e3 * PARSEC), 0.0);
        // Once the pulsar term arrives the residual stops growing.
        let l = 1e3 * PARSEC;
        let after = |dt: f64| pta_memory_residual(l / C + dt, 0.0, 1e-14, 0.0, FRAC_PI_2, l);
        assert!((after(10.0 * YEAR) - after(20.0 * YEAR)).abs() < 1e-18);
    }
}
