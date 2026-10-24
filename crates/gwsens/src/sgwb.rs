//! Stochastic backgrounds and pulsar timing: power laws, power-law-integrated sensitivity curves
//! ([Thrane & Romano 2013](https://arxiv.org/abs/1310.5300)) and the Hellings–Downs curve.
//!
//! ```
//! use gwsens::sgwb::hellings_downs;
//!
//! // Pulsars 90° apart on the sky are anti-correlated.
//! assert!((hellings_downs(90f64.to_radians()) + 0.14486).abs() < 1e-5);
//! ```

use crate::{YEAR, integrate_log, strain::omega_from_sh};

/// The Hellings–Downs correlation between two pulsars an angle `psi` (radians) apart:
/// `(3/2) x ln x − x/4 + 1/2`, with `x = (1 − cos ψ)/2` (no pulsar-term auto-correlation).
pub fn hellings_downs(psi: f64) -> f64 {
    let x = (1.0 - psi.cos()) / 2.0;
    if x <= 0.0 { 0.5 } else { 1.5 * x * x.ln() - x / 4.0 + 0.5 }
}

/// Energy density of a pulsar-timing background `h_c = A (f/f_yr)^α`, with `α = (3 − γ)/2`
/// (γ = 13/3 for supermassive black-hole binaries), at `f`; `h` as in [`crate::strain`].
pub fn pta_omega(amplitude: f64, gamma: f64, f: f64, h: f64) -> f64 {
    let alpha = (3.0 - gamma) / 2.0;
    let hc = amplitude * (f * YEAR).powf(alpha);
    crate::strain::omega_from_hc(hc, f, h, crate::strain::Convention::OneSided)
}

/// The power-law-integrated sensitivity curve Ω_PI(f) at each frequency in `freqs`.
///
/// For each power-law index β in `betas`, the amplitude Ω_β is the one that reaches SNR
/// `snr` in observing time `t_obs` (s): `ρ = √(2T ∫ (Ω_β(f/f_ref)^β)²/Ω_eff² df)`, with
/// `Ω_eff(f)` the energy density of the network's effective noise `s_eff` (one-sided, Hz⁻¹;
/// see [`omega_from_sh`]). The curve is the envelope `max_β Ω_β (f/f_ref)^β`. The integral runs
/// over `band = (f_lo, f_hi)`; `h` is the Hubble parameter as elsewhere.
pub fn power_law_integrated(freqs: &[f64], s_eff: impl Fn(f64) -> f64, band: (f64, f64), t_obs: f64, snr: f64, h: f64, betas: &[f64]) -> Vec<f64> {
    let (f_lo, f_hi) = band;
    let f_ref = (f_lo * f_hi).sqrt();
    let omega_eff = |f: f64| omega_from_sh(s_eff(f), f, h);
    let amps: Vec<(f64, f64)> = betas.iter().map(|&beta| {
        let integral = integrate_log(f_lo, f_hi, 4000, |f| (f / f_ref).powf(2.0 * beta) / omega_eff(f).powi(2));
        (beta, snr / (2.0 * t_obs * integral).sqrt())
    }).collect();
    freqs.iter().map(|&f| amps.iter().map(|&(beta, a)| a * (f / f_ref).powf(beta)).fold(0.0, f64::max)).collect()
}

/// The integer power-law indices −8 … 8 used by Thrane & Romano.
pub fn standard_betas() -> Vec<f64> {
    (-8..=8).map(f64::from).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hellings_downs_values() {
        for (deg, expect) in [(10.0, 0.44250), (30.0, 0.21163), (60.0, -0.08236), (90.0, -0.14486), (120.0, -0.01114), (150.0, 0.16971), (180.0, 0.25)] {
            let v = hellings_downs(f64::to_radians(deg));
            assert!((v - expect).abs() < 2e-5, "{deg}: {v}");
        }
        assert_eq!(hellings_downs(0.0), 0.5);
    }

    #[test]
    fn nanograv_energy_density() {
        let v = pta_omega(2.4e-15, 13.0 / 3.0, 1.0 / YEAR, 1.0);
        assert!((v / 3.62e-9 - 1.0).abs() < 2e-3, "{v:e}");
    }

    #[test]
    fn pi_curve_touches_every_power_law_it_can_detect() {
        // For flat effective noise in Ω, a β = 0 background at the curve's minimum has SNR = 1.
        let s_eff = |f: f64| 1e-45 * (f / 100.0).powi(-3); // Ω_eff constant
        let freqs = [20.0, 100.0, 500.0];
        let pi = power_law_integrated(&freqs, s_eff, (20.0, 500.0), YEAR, 1.0, 1.0, &standard_betas());
        let omega_eff = omega_from_sh(s_eff(100.0), 100.0, 1.0);
        let flat = omega_eff / (2.0 * YEAR * (500.0 - 20.0)).sqrt();
        assert!(pi.iter().all(|&v| v >= flat * (1.0 - 1e-9)), "{pi:?} vs {flat:e}");
        assert!(pi[1] < pi[0] && pi[1] < pi[2]);
    }
}
