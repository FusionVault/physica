//! Converting between the three ways to state the size of a gravitational-wave background:
//! energy density Ω_GW, strain power spectral density S_h, and characteristic strain h_c.
//!
//! The standard one-sided convention ([Moore, Cole & Berry 2015](https://arxiv.org/abs/1408.0740);
//! [Thrane & Romano 2013](https://arxiv.org/abs/1310.5300)) is
//! `H₀² Ω_GW = (2π²/3) f³ S_h = (2π²/3) f² h_c²`, with `h_c² = f·S_h`. The 2025 Living Review on
//! high-frequency waves ([Aggarwal et al.](https://arxiv.org/abs/2501.11723)) uses two-sided
//! spectra, so its h_c is √2 smaller for the same Ω. Pass the [`Convention`] explicitly.
//!
//! ```
//! use gwsens::strain::{Convention, hc_from_omega};
//!
//! // The cosmological ceiling on any early-Universe background, at 1 GHz.
//! let hc = hc_from_omega(1.1e-6, 1e9, 1.0, Convention::OneSided);
//! assert!((hc - 1.325e-30).abs() < 0.01e-30);
//! ```

use crate::H100;
use core::f64::consts::PI;

/// Spectral convention.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Convention {
    /// One-sided spectra (Moore et al. 2015, Thrane & Romano 2013, LIGO and LISA usage).
    OneSided,
    /// Two-sided spectra of the 2025 high-frequency Living Review: h_c smaller by √2.
    TwoSidedLrr2025,
}

impl Convention {
    fn hc_factor(self) -> f64 {
        match self {
            Convention::OneSided => 1.0,
            Convention::TwoSidedLrr2025 => core::f64::consts::FRAC_1_SQRT_2,
        }
    }
}

/// Characteristic strain h_c at frequency `f` (Hz) for energy density `omega` (Ω_GW),
/// with Hubble parameter `h` (H₀ = h × 100 km/s/Mpc; pass 1.0 and Ω h² as `omega` to work in
/// the usual h-free combination).
pub fn hc_from_omega(omega: f64, f: f64, h: f64, convention: Convention) -> f64 {
    let h0 = h * H100;
    (3.0 * h0 * h0 * omega / (2.0 * PI * PI)).sqrt() / f * convention.hc_factor()
}

/// Energy density Ω_GW from characteristic strain (inverse of [`hc_from_omega`]).
pub fn omega_from_hc(hc: f64, f: f64, h: f64, convention: Convention) -> f64 {
    let h0 = h * H100;
    let hc = hc / convention.hc_factor();
    2.0 * PI * PI * f * f * hc * hc / (3.0 * h0 * h0)
}

/// Strain power spectral density S_h (Hz⁻¹) from characteristic strain: `S_h = h_c²/f`.
pub fn sh_from_hc(hc: f64, f: f64) -> f64 {
    hc * hc / f
}

/// Characteristic strain from S_h: `h_c = √(f·S_h)`.
pub fn hc_from_sh(sh: f64, f: f64) -> f64 {
    (f * sh).sqrt()
}

/// Energy density Ω_GW of a detector's noise S_n (or of a signal with spectrum S_h):
/// `Ω = (2π²/3H₀²) f³ S_h`, one-sided.
pub fn omega_from_sh(sh: f64, f: f64, h: f64) -> f64 {
    let h0 = h * H100;
    2.0 * PI * PI * f.powi(3) * sh / (3.0 * h0 * h0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn prefactors_and_round_trips() {
        // h_c = 1.2634×10⁻¹⁸ (1 Hz / f) √(h²Ω), one-sided; ÷√2 two-sided.
        assert!((hc_from_omega(1.0, 1.0, 1.0, Convention::OneSided) / 1.2634e-18 - 1.0).abs() < 1e-4);
        assert!((hc_from_omega(1.0, 1.0, 1.0, Convention::TwoSidedLrr2025) / 8.934e-19 - 1.0).abs() < 1e-4);
        for conv in [Convention::OneSided, Convention::TwoSidedLrr2025] {
            let hc = hc_from_omega(3.7e-9, 2.3e-8, 0.674, conv);
            assert!((omega_from_hc(hc, 2.3e-8, 0.674, conv) / 3.7e-9 - 1.0).abs() < 1e-12);
        }
        let sh = sh_from_hc(2.0e-20, 100.0);
        assert!((hc_from_sh(sh, 100.0) / 2.0e-20 - 1.0).abs() < 1e-12);
        assert!((omega_from_sh(sh, 100.0, 0.7) / omega_from_hc(2.0e-20, 100.0, 0.7, Convention::OneSided) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn nanograv_amplitude() {
        // NANOGrav 15-yr: A = 2.4×10⁻¹⁵ at f = 1/yr gives Ω h² = 3.62×10⁻⁹.
        let f = 1.0 / crate::YEAR;
        assert!((omega_from_hc(2.4e-15, f, 1.0, Convention::OneSided) / 3.62e-9 - 1.0).abs() < 2e-3);
    }
}
