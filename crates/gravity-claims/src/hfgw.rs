//! The Li–Baker high-frequency gravitational-wave detector and the Pais patents.
//!
//! **Status: refuted.** [JASON (2008)](https://web.archive.org/web/2010/http://www.fas.org/irp/agency/dod/jason/gravwaves.pdf)
//! found the Li–Baker sensitivity claim too optimistic by more than 20 orders of magnitude. The
//! detector's published "signal" is the interference term between its reference beam and the
//! photons a wave converts, `2√(Ṅ₀Ṅ_GW)`. The reference beam's own shot noise cancels that gain,
//! so the time to see a wave is about `1/Ṅ_GW` whatever the beam power. The Li–Baker papers also
//! use a wave energy flux `c³ω²A²/(8πG)`, four times the standard one for a single polarization.
//!
//! The Pais patents' headline numbers come from `f_G(σ²/ε₀)v`, which reproduces them only with
//! surface fields millions of times above the vacuum breakdown limit his own patent states. A
//! Navy test (NAWCAD, 2016–2019, $508k) could not show the "Pais Effect".
//!
//! ```
//! use gravity_claims::hfgw::{interference_rate, shot_noise_time};
//!
//! // Li et al.'s parameters convert 2.56×10⁻¹⁸ photons/s. Their "signal" at x = 3.5 cm is the
//! // interference with 1.24×10²² background photons/s: 357/s (published 354/s).
//! assert!((interference_rate(1.24e22, 2.56e-18) / 354.0 - 1.0).abs() < 0.01);
//! // Read against the background's shot noise, one signal photon takes 3.9×10¹⁷ s, not 3×10³ s.
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
/// returned as (Pais's use, SI). Pais applies JASON's Gaussian-unit formula `(4πG/c⁴)B²L²·P`
/// to SI numbers; the SI first-order result is [`gwsens::hf::gertsenshtein_probability`]·P.
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
        // JASON: one signal photon in 3.9×10¹⁷ s; Li et al. claimed 3×10³ s.
        assert!(close(shot_noise_time(2.56e-18, 2.0), 3.9e17, 1e-2));
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
