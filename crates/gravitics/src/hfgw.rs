//! High-frequency gravitational waves: the Li–Baker detector and generators, and the Pais patents.
//!
//! **The proposals.** Li, Baker and colleagues proposed detecting GHz waves by their conversion
//! to photons in a magnet threaded by a Gaussian laser beam. Their published "signal" is the
//! interference term between the reference beam and the converted photons, `2√(Ṅ₀Ṅ_GW)`, which
//! grows with the beam power. Their papers use a wave energy flux `c³ω²A²/(8πG)`, four times the
//! standard one for a single polarization, which makes their conversion rate eight times the
//! standard one. Pais's patents give radiated fluxes `f_G(σ²/ε₀)v` from charged spinning or
//! vibrating surfaces.
//!
//! **Standard physics.** [JASON (2008)](https://web.archive.org/web/2010/http://www.fas.org/irp/agency/dod/jason/gravwaves.pdf)
//! showed (§4.1) that the reference beam's own shot noise cancels the interference gain: the time
//! to see a wave is about `1/Ṅ_GW` whatever the beam power, and the sensitivity claim is "more than
//! 20 orders of magnitude" too optimistic. Pais's fluxes need surface fields millions of times
//! above the vacuum breakdown limit his own patent states.
//!
//! **Measured.** No Li–Baker detector was built. A Navy test of a Pais device (NAWCAD, 2016–2019)
//! did not show the effect.
//!
//! ```
//! use gravitics::hfgw::{interference_rate, shot_noise_time};
//!
//! // Li et al.'s parameters convert 2.56×10⁻¹⁸ photons/s in their convention. Their "signal" at
//! // x = 3.5 cm is the interference with 1.24×10²² background photons/s: 357/s (published 354/s).
//! assert!((interference_rate(1.24e22, 2.56e-18) / 354.0 - 1.0).abs() < 0.01);
//! // Against the background's shot noise, one signal photon takes 3.9×10¹⁷ s at their rate
//! // (3.1×10¹⁸ s at the standard rate), against the 3×10³ s claimed.
//! assert!(shot_noise_time(2.56e-18, 2.0) / 3e3 > 1e14);
//! ```

use core::f64::consts::PI;

use crate::{EPS_0, HBAR};

/// Which energy flux a strain amplitude carries.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Flux {
    /// `c³ω²h²/(32πG)`, one linear polarization.
    Standard,
    /// `c³ω²A²/(8πG)`, as in the Li–Baker papers: four times the standard flux.
    LiBaker,
}

/// Energy flux (W/m²) of a wave of strain amplitude `h` at frequency `f` (Hz).
pub fn wave_flux(h: f64, f: f64, convention: Flux) -> f64 {
    let standard = gwsens::lab::wave_flux(h, f);
    match convention {
        Flux::Standard => standard,
        Flux::LiBaker => 4.0 * standard,
    }
}

/// Strain amplitude of a wave carrying `flux` (W/m²) at `f` (Hz): the inverse of [`wave_flux`].
/// For the same flux, Baker's strain is half the standard one.
pub fn strain_from_flux(flux: f64, f: f64, convention: Flux) -> f64 {
    (flux / wave_flux(1.0, f, convention)).sqrt()
}

/// Gravitons per second through `area` (m²) of a wave with energy flux `flux` at `f`:
/// `flux·area/ħω`.
pub fn graviton_rate(flux: f64, area: f64, f: f64) -> f64 {
    flux * area / (HBAR * 2.0 * PI * f)
}

/// Rate (photons/s) of the interference term between a background of `n0` photons/s and `n_gw`
/// converted photons/s: `2√(n0·n_gw)`. This is the Li group's published "signal".
pub fn interference_rate(n0: f64, n_gw: f64) -> f64 {
    2.0 * (n0 * n_gw).sqrt()
}

/// Time (s) for the interference term to reach `snr` against the shot noise of the background
/// that makes it: `snr²/(4·n_gw)`. The background cancels (JASON 2008, §4.1).
pub fn shot_noise_time(n_gw: f64, snr: f64) -> f64 {
    snr * snr / (4.0 * n_gw)
}

/// The vacuum breakdown field Pais's own patent (US 10,135,366) states, V/m.
pub const VACUUM_BREAKDOWN: f64 = 1e9;

/// Pais's maximum radiated flux (W/m²), `f_G·(σ²/ε₀)·v` (US 10,144,532, eq. 1), for a surface
/// charge `sigma` (C/m²) moving at `v` (m/s); `f_g` is 1 for a disc.
pub fn pais_flux(f_g: f64, sigma: f64, v: f64) -> f64 {
    f_g * sigma * sigma / EPS_0 * v
}

/// The surface field `σ/ε₀` (V/m) of a surface charge `sigma` (C/m²).
pub fn surface_field(sigma: f64) -> f64 {
    sigma / EPS_0
}

/// Gravitational-wave power (W) from light of `power` (W) crossing a field `b` (T) over `l` (m),
/// returned as (Pais's use, SI). Pais applies JASON's formula `(4πG/c⁴)B²L²·P` (Eq. 3-14, written
/// in Heaviside–Lorentz form; the Gaussian factor is `GB²L²/c⁴`) to SI numbers; the SI
/// first-order result is [`gwsens::hf::gertsenshtein_probability`]·P. His US 10,322,827 gives an
/// input flux of 10³³ W/m² and no length: reproducing his "10¹⁰ W" needs 1 m and 1 m².
pub fn pais_conversion(b: f64, l: f64, power: f64) -> (f64, f64) {
    let si = gwsens::hf::gertsenshtein_probability(b, l) * power;
    (si * gwsens::MU_0, si)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, rel: f64) -> bool {
        (a / b - 1.0).abs() < rel
    }

    #[test]
    fn baker_strain_is_half_the_standard_one() {
        // Baker's FBAR case: 3.53 W/m² at 4.915 GHz gives 4.8×10⁻²⁸ (published 4.9×10⁻²⁸).
        assert!(close(strain_from_flux(3.53, 4.915e9, Flux::LiBaker), 4.80e-28, 2e-3));
        assert!(close(strain_from_flux(3.53, 4.915e9, Flux::Standard), 9.60e-28, 2e-3));
    }

    #[test]
    fn li_detector_numbers_reproduce() {
        // Graviton flux through the 5 cm waist for h = 10⁻³⁰ at 5 GHz: 3.77×10¹⁶ s⁻¹ published.
        let flux = wave_flux(1e-30, 5e9, Flux::LiBaker);
        assert!(close(graviton_rate(flux, PI * 0.05 * 0.05, 5e9), 3.77e16, 5e-3));
        // Table I's signal at x = 3.5 cm is the interference term (ratio 1.007).
        assert!(close(interference_rate(1.24e22, 2.56e-18), 3.54e2, 1e-2));
        // Li's rate is 8× the standard one for the same 3 T, 6 m, 0.01 m², 5 GHz and h = 10⁻³⁰.
        let standard = gwsens::hf::coherent_photon_rate(1.0, 0.01, 3.0, 6.0, 5e9, 1e-30);
        assert!(close(standard, 3.20e-19, 2e-3));
        assert!(close(2.56e-18 / standard, 8.0, 0.01));
        // JASON's shot-noise argument: one signal photon in 3.9×10¹⁷ s at Li's rate, 3.1×10¹⁸ s at
        // the standard rate; Li et al. claimed 3×10³ s.
        assert!(close(shot_noise_time(2.56e-18, 2.0), 3.9e17, 1e-2));
        assert!(close(shot_noise_time(standard, 2.0), 3.12e18, 1e-2));
    }

    #[test]
    fn pais_numbers_need_impossible_fields() {
        // US 10,144,532: σ = 5×10⁴ C/m² on a 2 m disc at 30,000 rpm gives ~10²⁴ W/m².
        assert!(close(pais_flux(1.0, 5e4, 2.0 * 3142.0), 1.77e24, 5e-3));
        assert!(close(surface_field(5e4) / VACUUM_BREAKDOWN, 5.65e6, 5e-3));
        // Capped at his own breakdown limit, the same disc gives 14 orders less.
        assert!(close(pais_flux(1.0, VACUUM_BREAKDOWN * EPS_0, 2.0 * 3142.0), 5.56e10, 5e-3));
        // 10¹⁰ T over 1 m on 10³³ W: "10¹⁰ W" as he uses it, 8.3×10¹⁵ W in consistent SI.
        let (as_used, si) = pais_conversion(1e10, 1.0, 1e33);
        assert!(close(as_used, 1.04e10, 1e-2));
        assert!(close(si, 8.3e15, 1e-2));
    }
}
