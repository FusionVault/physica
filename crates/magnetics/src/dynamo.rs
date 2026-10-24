//! How a moving conductor makes its own magnetic field: the mean-field αΩ dynamo's onset, and the
//! thresholds the laboratory dynamos actually reached.
//!
//! Theory gives thresholds, not predictions, and the experiments show why: VKS2 self-excited at a
//! magnetic Reynolds number of 31–35, where simulations had predicted 48–133, because its soft-iron
//! impellers changed the boundary ([Stefani, Gailitis & Gerbeth 2008](https://arxiv.org/abs/0807.0299)).
//!
//! ```
//! use magnetics::dynamo::{ALPHA_OMEGA_CRITICAL, alpha_omega_growth};
//!
//! // At the critical dynamo number αS/(η²k³) = 2 the plane wave neither grows nor decays, and
//! // its cycle frequency is ηk².
//! let (eta, k, shear): (f64, f64, f64) = (1.0, 1.0, 1414.2);
//! let alpha = ALPHA_OMEGA_CRITICAL * eta * eta * k.powi(3) / shear;
//! let (growth, frequency) = alpha_omega_growth(alpha, shear, eta, k);
//! assert!(growth.abs() < 1e-3 && (frequency - 1.0).abs() < 1e-3);
//! ```

/// The αΩ dynamo number at which a plane dynamo wave starts to grow.
pub const ALPHA_OMEGA_CRITICAL: f64 = 2.0;

/// Published self-excitation thresholds: dynamo or flow, the critical number, and what it is.
pub const LAB_THRESHOLDS: [(&str, f64, &str); 7] = [
    ("Ponomarenko flow", 17.7, "Rm, theory"),
    ("Riga, optimised profile (convective)", 14.3, "Rm, theory"),
    ("Riga, optimised profile (absolute)", 17.6, "Rm, theory"),
    ("Karlsruhe (predicted)", 8.12, "α-number C, mean-field theory; measured 8.4–9.3"),
    ("VKS2 (measured)", 31.0, "Rm, lower edge of 31–35; predicted 48–133"),
    ("Roberts flow", 3.49, "Rm, k_z-optimised"),
    ("Small-scale dynamo at Pm = 1", 35.0, "Rm, simulations"),
];

/// The αΩ dynamo number `αS/(η²k³)` for α effect `alpha` (m/s), shear `shear` (s⁻¹), turbulent
/// diffusivity `eta` (m²/s) and wavenumber `k` (m⁻¹).
pub fn dynamo_number(alpha: f64, shear: f64, eta: f64, k: f64) -> f64 {
    alpha * shear / (eta * eta * k.powi(3))
}

/// Growth rate and cycle frequency (s⁻¹) of a plane mean-field dynamo wave along the shear's
/// gradient: `λ = −ηk² + (α²k² − iαSk)^(1/2)` (the root with positive real part of the square root).
pub fn alpha_omega_growth(alpha: f64, shear: f64, eta: f64, k: f64) -> (f64, f64) {
    // √(a + ib) on the principal branch.
    let (a, b) = (alpha * alpha * k * k, -alpha * shear * k);
    let modulus = a.hypot(b);
    let (re, im) = (((modulus + a) / 2.0).sqrt(), ((modulus - a) / 2.0).sqrt().copysign(b));
    (re - eta * k * k, im.abs())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn onset_at_dynamo_number_two() {
        // αk/S = 10⁻⁴, η = k = 1: D = 2 is the onset, with frequency ηk² (Brandenburg & Subramanian).
        let shear = (2.0f64 / 1e-4).sqrt();
        let alpha = 1e-4 * shear;
        assert!((dynamo_number(alpha, shear, 1.0, 1.0) - 2.0).abs() < 1e-12);
        let (growth, frequency) = alpha_omega_growth(alpha, shear, 1.0, 1.0);
        assert!(growth.abs() < 2e-4 && (frequency - 1.0).abs() < 2e-4, "{growth} {frequency}");
        // Below onset it decays, above it grows.
        assert!(alpha_omega_growth(0.5 * alpha, shear, 1.0, 1.0).0 < 0.0 && alpha_omega_growth(2.0 * alpha, shear, 1.0, 1.0).0 > 0.0);
        // A pure α² dynamo (no shear) grows once α > ηk.
        assert!((alpha_omega_growth(2.0, 0.0, 1.0, 1.0).0 - 1.0).abs() < 1e-12);
    }
}
