//! Gravitational waves as seen by the sensors this crate models: long-baseline atom
//! interferometers, gradiometers, optical clocks and satellite pairs.
//!
//! A wave of strain h stretches a baseline L, so two free masses feel a tidal acceleration
//! `Δa = ½ω²hL`. A gravity gradiometer and a wave detector therefore measure the same thing; what
//! separates them is the baseline. A GOCE-class gradiometer reads as about 10⁻¹¹ strain per √Hz
//! at 0.1 Hz, while an atom-interferometer detector with a kilometre baseline reaches about 10⁻²¹.
//!
//! ```
//! use satgrav::gw::{AtomGwDetector, K_SR87};
//!
//! // AION-km: 2 km, T = 5 s, 40 000 photon kicks, 3×10⁻⁶ rad/√Hz (Badurina et al. 2020).
//! let aion_km = AtomGwDetector { baseline: 2000.0, t: 5.0, lmt: 40_000.0, phase_noise: 3e-6, k: K_SR87, diamonds: 1.0 };
//! assert!((aion_km.peak_strain() / 2.08e-21 - 1.0).abs() < 0.01);
//! ```

use core::f64::consts::PI;

/// Wavenumber `k = ω_A/c = 2π/λ` of the ⁸⁷Sr clock transition (698.4457 nm), used single-photon by
/// long-baseline detectors (AION, MAGIS-100, AICE), m⁻¹.
pub const K_SR87: f64 = 2.0 * PI / 698.445_7e-9;

const C: f64 = 299_792_458.0;

/// An atom-interferometer gravitational-wave detector: two atom clouds a `baseline` apart
/// interrogated by the same laser ([Graham et al. 2016](https://arxiv.org/abs/1606.01860)).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AtomGwDetector {
    /// Distance between the atom clouds, m.
    pub baseline: f64,
    /// Time between pulses T, s.
    pub t: f64,
    /// Large-momentum-transfer order n (photon kicks per pulse).
    pub lmt: f64,
    /// Phase noise, rad/√Hz (an input: squeezing and interleaving make it differ from 1/√N).
    pub phase_noise: f64,
    /// Laser wavenumber, m⁻¹ (e.g. [`K_SR87`]).
    pub k: f64,
    /// Number of diamond sequences Q in resonant mode (1 for broadband).
    pub diamonds: f64,
}

impl AtomGwDetector {
    /// Phase per unit strain, |Δφ/h| (Graham et al. 2016, Eq. 9):
    /// `k_eff L · sin(ωQT)/cos(ωT/2) · sinc(ωnL/2c) · sin(ωT/2 − ω(n−1)L/2c)`, `k_eff = n·k`,
    /// optimally oriented.
    pub fn phase_response(&self, f: f64) -> f64 {
        let w = 2.0 * PI * f;
        let (n, q, t, l) = (self.lmt, self.diamonds, self.t, self.baseline);
        let ratio = if (w * t / 2.0).cos().abs() < 1e-9 {
            -2.0 * q * (w * q * t).cos() / (w * t / 2.0).sin() // the limit at ωT = π(2j+1)
        } else {
            (w * q * t).sin() / (w * t / 2.0).cos()
        };
        let x = w * n * l / (2.0 * C);
        let sinc = if x.abs() < 1e-12 { 1.0 } else { x.sin() / x };
        (n * self.k * l * ratio * sinc * (w * t / 2.0 - w * (n - 1.0) * l / (2.0 * C)).sin()).abs()
    }

    /// Strain noise, Hz⁻¹ᐟ²: phase noise over the response.
    pub fn strain_asd(&self, f: f64) -> f64 {
        self.phase_noise / self.phase_response(f)
    }

    /// Best broadband strain noise, `δφ/(2nkL)`, at f = 1/(2T) (light travel time neglected).
    pub fn peak_strain(&self) -> f64 {
        self.phase_noise / (2.0 * self.lmt * self.k * self.baseline)
    }
}

/// Tidal acceleration a wave of strain `h` at frequency `f` puts across a `baseline` (m):
/// `Δa = ½ω²hL`, m/s².
pub fn tidal_acceleration(h: f64, f: f64, baseline: f64) -> f64 {
    0.5 * (2.0 * PI * f).powi(2) * h * baseline
}

/// The strain noise of a gravity gradiometer with gradient noise `gradient_asd` (s⁻² Hz⁻¹ᐟ²) at
/// frequency `f`: `h = 2√S_Γ/ω²`.
pub fn gradiometer_strain(gradient_asd: f64, f: f64) -> f64 {
    2.0 * gradient_asd / (2.0 * PI * f).powi(2)
}

/// Strain noise of two satellites whose separation `separation` (m) is measured with noise
/// `range_asd` (m Hz⁻¹ᐟ²): `h = x/x₀` (as in [Ross et al. 2020](https://arxiv.org/abs/2002.02044)
/// for GRACE-FO).
pub fn pair_strain(range_asd: f64, separation: f64) -> f64 {
    range_asd / separation
}

/// Shot-noise strain limit of an optical-clock detector
/// ([Kolkowitz et al. 2016](https://arxiv.org/abs/1606.01859)): `√Δ/(ν√(2πτN))`, Hz⁻¹ᐟ², for
/// clock frequency `nu` (Hz), linewidth `linewidth` (Hz), interrogation time `tau` (s) and
/// `atoms`.
pub fn clock_strain(nu: f64, linewidth: f64, tau: f64, atoms: f64) -> f64 {
    linewidth.sqrt() / (nu * (2.0 * PI * tau * atoms).sqrt())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn design(baseline: f64, t: f64, lmt: f64, phase_noise: f64) -> AtomGwDetector {
        AtomGwDetector { baseline, t, lmt, phase_noise, k: K_SR87, diamonds: 1.0 }
    }

    #[test]
    fn published_design_strains() {
        // Badurina et al. 2020 (AION) and Baynham et al. 2025 (AICE), δφ/(2nkL).
        assert!((design(2000.0, 5.0, 40_000.0, 3e-6).peak_strain() / 2.08e-21 - 1.0).abs() < 0.01);
        assert!((design(50.0, 1.4, 5494.0, 1e-5).peak_strain() / 2.02e-18 - 1.0).abs() < 0.01);
        assert!((design(100.0, 1.4, 1000.0, 1e-4).peak_strain() / 5.55e-17 - 1.0).abs() < 0.01);
        assert!((design(10.0, 1.4, 100.0, 1e-3).peak_strain() / 5.6e-14 - 1.0).abs() < 0.01);
    }

    #[test]
    fn response_limits() {
        // Without light-travel time, one diamond gives 2nkL·sin²(ωT/2): the tidal-acceleration
        // picture fed through the interferometer's transfer function.
        let d = AtomGwDetector { baseline: 1e-3, t: 1.0, lmt: 10.0, phase_noise: 1e-5, k: K_SR87, diamonds: 1.0 };
        for f in [0.05, 0.2, 0.37] {
            let expect = 2.0 * d.lmt * d.k * d.baseline * (PI * f * d.t).sin().powi(2);
            assert!((d.phase_response(f) / expect - 1.0).abs() < 1e-6, "{f}");
        }
        // At f = 1/(2T) the response peaks at 2Q·nkL.
        let q3 = AtomGwDetector { diamonds: 3.0, ..d };
        assert!((q3.phase_response(0.5) / (2.0 * 3.0 * d.lmt * d.k * d.baseline) - 1.0).abs() < 1e-6);
        assert!((d.strain_asd(0.5) / d.peak_strain() - 1.0).abs() < 1e-6);
    }

    #[test]
    fn gradiometers_clocks_and_pairs() {
        // A GOCE-class 10 mE/√Hz gradiometer read as a wave detector.
        assert!((gradiometer_strain(1e-11, 0.01) / 5.07e-9 - 1.0).abs() < 2e-3);
        assert!((gradiometer_strain(1e-11, 0.1) / 5.07e-11 - 1.0).abs() < 2e-3);
        // Kolkowitz et al.: Sr at 429.2 THz, 1 mHz linewidth, 7×10⁶ atoms.
        assert!((clock_strain(429.2e12, 1e-3, 1.0, 7e6) / 1.11e-20 - 1.0).abs() < 2e-3);
        // GRACE-FO: about 1 nm/√Hz of ranging over 220 km.
        assert!((pair_strain(1e-9, 2.2e5) - 4.5e-15).abs() < 1e-16);
        // Tidal acceleration and gradiometer strain are consistent: Δa/L = ½ω²h.
        let h = gradiometer_strain(1e-11, 0.1);
        assert!((tidal_acceleration(h, 0.1, 1.0) / 1e-11 - 1.0).abs() < 1e-12);
    }
}
