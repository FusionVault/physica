//! When a magnetised plasma holds together: Kadomtsev's criteria for the sausage (m = 0) and kink
//! (m ≥ 1) modes of a pinch, magnetic Rayleigh–Taylor, and resistive tearing.
//!
//! The pinch criteria take the local logarithmic pressure gradient and beta, so they apply to any
//! profile. On the Bennett pinch (`p ∝ (1 + x²)⁻²`, local β = 1/x², the profile behind
//! `plasmoid::pinch`) the kink mode is marginal at `x = √3` and the sausage mode is unstable
//! everywhere unless γ ≥ 2. Tearing follows Furth–Killeen–Rosenbluth; once the layer would be
//! thinner than the Sweet–Parker sheet, plasmoids take over (`plasmoid::relaxation::PLASMOID_RATE`).
//!
//! ```
//! use magnetics::stability::kink_stable;
//!
//! // Bennett pinch at x = r/r_p = 1: d ln p/d ln r = −2, β = 1. Kink-stable for m = 1? No.
//! assert!(!kink_stable(-2.0, 1.0, 1.0));
//! ```

use core::f64::consts::PI;

use metrology::MU_0;

/// Furth–Killeen–Rosenbluth coefficient `[Γ(¼)/(2πΓ(¾))]^(4/5)`.
pub const FKR_COEFFICIENT: f64 = 0.547_436_725;

/// Kadomtsev's sausage (m = 0) criterion: stable while `−d ln p/d ln r ≤ 4γ/(2 + γβ)`.
/// `dlnp_dlnr` is the profile's logarithmic gradient (negative where pressure falls outward).
pub fn sausage_stable(dlnp_dlnr: f64, beta: f64, gamma: f64) -> bool {
    -dlnp_dlnr <= 4.0 * gamma / (2.0 + gamma * beta)
}

/// Kadomtsev's kink criterion for mode `m ≥ 1`: stable while `−d ln p/d ln r ≤ m²/β`.
pub fn kink_stable(dlnp_dlnr: f64, beta: f64, m: f64) -> bool {
    -dlnp_dlnr <= m * m / beta
}

/// The shortest magnetic Rayleigh–Taylor wavelength that grows (m), `4πB²cos²θ/(μ₀gΔρ)`, for a
/// heavy fluid `rho_h` above a light one `rho_l` (kg/m³) in gravity `g`, with the field `b` (T) at
/// angle `theta` (radians) to the wave vector. The fastest growth is at twice this.
pub fn rt_cutoff(b: f64, theta: f64, rho_h: f64, rho_l: f64, g: f64) -> f64 {
    4.0 * PI * b * b * theta.cos().powi(2) / (MU_0 * g * (rho_h - rho_l))
}

/// Tearing-mode growth rate in Alfvén times, `γτ_A = 0.547 (Δ′a)^(4/5) S^(−3/5)`, for stability
/// index `delta_prime_a` and Lundquist number `s`. `None` outside the constant-ψ regime
/// `Δ′a ≪ S^(1/3)` (taken as `Δ′a < S^(1/3)/10`), where the formula fails.
pub fn fkr_growth(delta_prime_a: f64, s: f64) -> Option<f64> {
    (delta_prime_a < s.cbrt() / 10.0).then(|| FKR_COEFFICIENT * delta_prime_a.powf(0.8) * s.powf(-0.6))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The Bennett pinch at x = r/r_p: (d ln p/d ln r, local β).
    fn bennett(x: f64) -> (f64, f64) {
        (-4.0 * x * x / (1.0 + x * x), 1.0 / (x * x))
    }

    #[test]
    fn bennett_pinch() {
        // Kink: marginal at x = √3, β = 1/3 (Kadomtsev; Crews & Turner).
        let (g, b) = bennett(3f64.sqrt());
        assert!((-g - 1.0 / b).abs() < 1e-12 && (b - 1.0 / 3.0).abs() < 1e-12);
        assert!(!kink_stable(bennett(1.5).0, bennett(1.5).1, 1.0) && kink_stable(bennett(2.0).0, bennett(2.0).1, 1.0));
        // Sausage: unstable at every radius for γ = 5/3; exactly marginal for γ = 2.
        for x in [0.5, 1.0, 3f64.sqrt(), 3.0] {
            let (g, b) = bennett(x);
            assert!(!sausage_stable(g, b, 5.0 / 3.0), "{x}");
            assert!((-g - 8.0 / (2.0 + 2.0 * b)).abs() < 1e-12);
        }
    }

    #[test]
    fn prominence_rayleigh_taylor() {
        // Khomenko et al. 2014: 10 G at 89° and 88°, ρ 2.9×10⁻¹⁰ over 3.7×10⁻¹² kg/m³, solar g:
        // 38 km and 155 km.
        let cut = |deg: f64| rt_cutoff(1e-3, deg.to_radians(), 2.9e-10, 3.7e-12, 274.0);
        assert!((cut(89.0) / 38.827e3 - 1.0).abs() < 1e-4 && (cut(88.0) / 155.262e3 - 1.0).abs() < 1e-4);
    }

    #[test]
    fn tearing() {
        // Γ(¼) = 3.625 609 908 2, Γ(¾) = 1.225 416 702 5.
        assert!((FKR_COEFFICIENT - (3.625_609_908_2 / (2.0 * PI * 1.225_416_702_5)).powf(0.8)).abs() < 1e-9);
        assert!((fkr_growth(1.0, 1e6).unwrap() / 1.375_099e-4 - 1.0).abs() < 1e-6);
        assert!((fkr_growth(1.0, 1e8).unwrap() / 8.676e-6 - 1.0).abs() < 1e-3);
        assert!(fkr_growth(100.0, 1e6).is_none());
    }
}
