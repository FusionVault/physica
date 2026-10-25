#![doc = include_str!("../README.md")]

pub mod curves;
pub mod detectors;
pub mod energy;
pub mod hf;
pub mod inspiral;
pub mod laboratory;
pub mod memory;
pub mod propagation;
pub mod quanta;
pub mod ringdown;
pub mod sgwb;
pub mod strain;
pub mod ultralight;

/// Hubble constant for h = 1 (100 km s⁻¹ Mpc⁻¹), s⁻¹. Multiply by the dimensionless h.
pub const H100: f64 = 1.0e5 / metrology::MPC;

/// `∫ g(f) df` from `lo` to `hi` (both > 0) by Simpson's rule on a logarithmic grid of `n`
/// intervals (`n` is made even). Accurate for the smooth power-law integrands used here.
pub fn integrate_log(lo: f64, hi: f64, n: usize, g: impl Fn(f64) -> f64) -> f64 {
    let n = (n.max(2) + 1) & !1;
    let (a, b) = (lo.ln(), hi.ln());
    let h = (b - a) / n as f64;
    let term = |k: usize| { let f = (a + k as f64 * h).exp(); g(f) * f };
    let inner: f64 = (1..n).map(|k| if k % 2 == 1 { 4.0 } else { 2.0 } * term(k)).sum();
    (term(0) + inner + term(n)) * h / 3.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn log_integral_is_exact_for_power_laws() {
        // ∫₁¹⁰ f⁻² df = 0.9
        assert!((integrate_log(1.0, 10.0, 200, |f| f.powi(-2)) - 0.9).abs() < 1e-9);
        assert!((H100 - 3.2408e-18).abs() < 1e-22);
    }
}
