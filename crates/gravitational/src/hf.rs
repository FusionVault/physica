//! High-frequency gravitational waves (above 10 kHz): the cosmological ceiling, conversion of
//! waves into photons in a magnetic field, resonant-cavity noise, source frequencies, and a
//! registry of published limits, each with its status.
//!
//! The main reference is the 2025 Living Review by
//! [Aggarwal et al.](https://arxiv.org/abs/2501.11723). Its key message is in [`CEILING_BBN_CMB`]:
//! any background from the early Universe is capped at Ω h² ≈ 10⁻⁶, about 10⁻³⁰ strain at 1 GHz,
//! while the best laboratory limits sit about 19–26 orders of magnitude above it.
//!
//! ```
//! use gravitational::hf;
//!
//! // Probability that a wave turns into a photon over 9.26 m of a 9 T magnet (the CAST helioscope).
//! let p = hf::gertsenshtein_probability(9.0, 9.26);
//! assert!((p / 5.74e-34 - 1.0).abs() < 0.01);
//! ```

use metrology::{C, E_CHARGE, EPS_0, EV, G, HBAR, K_B, M_E, MU_0};
use core::f64::consts::PI;

/// Big-bang nucleosynthesis and CMB bound on any background present before nucleosynthesis:
/// Ω_GW h² < 1.1 × 10⁻⁶ (ΔN_eff < 0.18; Review 2025 §4.2).
pub const CEILING_BBN_CMB: f64 = 1.1e-6;
/// CMB-only bound with homogeneous initial conditions: Ω_GW h² < 2.9 × 10⁻⁷.
pub const CEILING_CMB: f64 = 2.9e-7;
/// A future CMB measurement of N_eff to its theoretical precision: about 5.6 × 10⁻⁹.
pub const CEILING_FUTURE: f64 = 5.6e-9;

/// Ω_GW h² per unit ΔN_eff: 5.6 × 10⁻⁶.
pub fn ceiling_from_delta_neff(delta_neff: f64) -> f64 {
    5.6e-6 * delta_neff
}

/// Probability that a gravitational wave converts to a photon (inverse Gertsenshtein effect)
/// crossing a transverse magnetic field `b` (T) over length `l` (m), for `l` ≫ wavelength:
/// `P = 4πG B² L²/(μ₀ c⁴)`.
pub fn gertsenshtein_probability(b: f64, l: f64) -> f64 {
    4.0 * PI * G * b * b * l * l / (MU_0 * C.powi(4))
}

/// Photons per second from a flat stochastic background of characteristic strain `hc` over a
/// bandwidth `df` (Hz), for a magnet of field `b` (T), length `l` (m) and aperture `area` (m²)
/// with detection efficiency `eff`: `Ṅ = ε A B² L² h_c² 2πΔf/(4μ₀cħ)` (Ejlli et al. 2019, Eqs. 3–4).
/// G cancels: the conversion probability times the wave's energy flux. Like Ejlli et al., this
/// treats the isotropic background as a beam along the magnet; the 2025 Living Review's
/// field-of-view correction weakens the resulting limits to h_c ≈ 4×10⁻²³ (OSQAR II) and
/// 8×10⁻²⁶ (CAST).
pub fn stochastic_photon_rate(eff: f64, area: f64, b: f64, l: f64, hc: f64, df: f64) -> f64 {
    eff * area * b * b * l * l * hc * hc * 2.0 * PI * df / (4.0 * MU_0 * C * HBAR)
}

/// The characteristic strain that gives `rate` photons per second (inverse of
/// [`stochastic_photon_rate`]).
pub fn stochastic_strain_for_rate(rate: f64, eff: f64, area: f64, b: f64, l: f64, df: f64) -> f64 {
    (rate / stochastic_photon_rate(eff, area, b, l, 1.0, df)).sqrt()
}

/// Photons per second from a monochromatic wave of amplitude `h0` at frequency `f`:
/// `Ṅ = ε A B² L² ω h₀²/(8μ₀cħ)`.
pub fn coherent_photon_rate(eff: f64, area: f64, b: f64, l: f64, f: f64, h0: f64) -> f64 {
    eff * area * b * b * l * l * 2.0 * PI * f * h0 * h0 / (8.0 * MU_0 * C * HBAR)
}

/// One-sided strain noise of a resonant microwave cavity in a magnetic field, Hz⁻¹: twice the
/// Living Review 2025's two-sided Eq. 113 in SI, `S_h = 4 k_B T μ₀ c²/(ω³ η² Q B² L⁵)`, with system
/// temperature `t` (K), frequency `f`, coupling `eta`, quality factor `q`, field `b` (T) and
/// size `l` (m). One-sided like the rest of this crate, so it can go into
/// [`omega_from_sh`](crate::strain::omega_from_sh).
pub fn cavity_noise(t: f64, f: f64, eta: f64, q: f64, b: f64, l: f64) -> f64 {
    let w = 2.0 * PI * f;
    4.0 * K_B * t * MU_0 * C * C / (w.powi(3) * eta * eta * q * b * b * l.powi(5))
}

/// Gravitational-wave line from a superradiant cloud of bosons of mass `mass_ev` (eV) around a
/// black hole: `f = 2μc²/h`, Hz (binding corrections neglected).
pub fn superradiance_frequency(mass_ev: f64) -> f64 {
    2.0 * mass_ev * EV / (2.0 * PI * HBAR)
}

/// Plasma frequency `(1/2π)√(n_e e²/(ε₀m_e))` (Hz) of electron density `n_e` (m⁻³): 8.98 kHz at
/// 1 cm⁻³. A wave below it cannot convert to a propagating photon, so Gertsenshtein conversion in
/// the interstellar medium (≈ 0.03 cm⁻³, 1.56 kHz) is closed to the LIGO band.
pub fn plasma_frequency(n_e: f64) -> f64 {
    (n_e * E_CHARGE * E_CHARGE / (EPS_0 * M_E)).sqrt() / (2.0 * PI)
}

/// The electron density (m⁻³) at which [`plasma_frequency`] equals `f`: above it the plasma
/// reflects a wave at that frequency, below it refracts. 1.24×10¹⁸ m⁻³ for X-band (10 GHz).
pub fn critical_density(f: f64) -> f64 {
    (2.0 * PI * f).powi(2) * EPS_0 * M_E / (E_CHARGE * E_CHARGE)
}

/// Refractive index `√(1 − (f_p/f)²)` of an underdense plasma of electron density `n_e` (m⁻³) at
/// frequency `f` (Hz); `None` at or above [`critical_density`], where the wave does not
/// propagate. A graded index is what transformation optics bends rays with, and a static curved
/// spacetime acts on light as one. In isotropic coordinates, `ds² = −A c²dt² + B dx²`, the
/// flat-space index is `n = √(B/A)`, and near a mass `n − 1 ≈ 2GM/(rc²)`, which bends a ray
/// grazing the Sun by 1.75″. Taking `n = 1/√g₀₀` alone keeps only the time part and gives half
/// that, Einstein's 1911 value; it is right only together with the curved spatial metric.
///
/// ```
/// use metrology::{C, GM_SUN};
///
/// // First order, the deflection is twice the index excess at the impact parameter: 4GM/(c²b).
/// let limb = 6.957e8; // the Sun's radius, m
/// let arcsec = |rad: f64| rad.to_degrees() * 3600.0;
/// let full = arcsec(2.0 * 2.0 * GM_SUN / (C * C * limb)); // n − 1 = 2GM/(rc²)
/// let time_only = arcsec(2.0 * GM_SUN / (C * C * limb)); // n − 1 = GM/(rc²), from 1/√g₀₀
/// assert!((full - 1.7512).abs() < 1e-4 && (time_only - 0.8756).abs() < 1e-4);
/// ```
pub fn plasma_index(n_e: f64, f: f64) -> Option<f64> {
    let ratio = plasma_frequency(n_e) / f;
    (ratio < 1.0).then(|| (1.0 - ratio * ratio).sqrt())
}

/// Spectral shape of a first-order phase transition's sound waves, `S(x) = x³[7/(4 + 3x²)]^(7/2)`,
/// with x = f/f_peak (Review 2025 Eq. 101). S(1) = 1.
pub fn phase_transition_shape(x: f64) -> f64 {
    x.powi(3) * (7.0 / (4.0 + 3.0 * x * x)).powf(3.5)
}

/// How far a result is from being a real measurement.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Status {
    /// A limit published by the experiment from its own data.
    PublishedLimit,
    /// A limit derived later from an experiment built for something else (e.g. axion searches).
    Recast,
    /// A recast whose interpretation is contested in the literature.
    Contested,
    /// A claimed signal not accepted as gravitational waves.
    Disputed,
}

/// One entry of the high-frequency landscape.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Limit {
    /// Experiment or analysis.
    pub name: &'static str,
    /// Frequency range, Hz.
    pub band: (f64, f64),
    /// The quoted sensitivity or limit, in the unit of `quantity`.
    pub value: f64,
    /// What `value` is.
    pub quantity: &'static str,
    /// Its status.
    pub status: Status,
    /// Where it comes from.
    pub source: &'static str,
}

/// Published high-frequency limits and claims (from the 2025 Living Review and the papers it
/// cites; values as quoted there).
pub const REGISTRY: [Limit; 7] = [
    Limit { name: "QUEST (twin Michelsons, Cardiff)", band: (13e6, 80e6), value: 3e-20, quantity: "√S_h Hz⁻¹ᐟ² near 40 MHz", status: Status::PublishedLimit, source: "https://arxiv.org/abs/2410.09175" },
    Limit { name: "Fermilab Holometer", band: (1e6, 13e6), value: 1e-21, quantity: "√S_h Hz⁻¹ᐟ² (better than, 1–13 MHz)", status: Status::PublishedLimit, source: "https://arxiv.org/abs/1611.05560" },
    Limit { name: "CAPP-12T axion haloscope", band: (5.310e9, 5.312e9), value: 3.9e-21, quantity: "h₀ (90% CL)", status: Status::Recast, source: "https://arxiv.org/abs/2511.17817" },
    Limit { name: "Light shining through a wall (OSQAR, ALPS)", band: (2.7e14, 1.4e15), value: 6e-26, quantity: "h_c", status: Status::Contested, source: "https://arxiv.org/abs/1908.00232" },
    Limit { name: "Helioscopes, X-ray (CAST)", band: (5e18, 1.2e19), value: 5e-28, quantity: "h_c", status: Status::Contested, source: "https://arxiv.org/abs/1908.00232" },
    Limit { name: "MAGE quartz resonators", band: (4.993e6, 5.505e6), value: 1.35e-20, quantity: "√S_h Hz⁻¹ᐟ² measured at 4.993 MHz (1.10×10⁻²⁰ at 5.505 MHz)", status: Status::PublishedLimit, source: "https://arxiv.org/abs/2307.00715" },
    Limit { name: "2021 quartz-resonator rare events", band: (5.5e6, 5.5e6), value: 1e8, quantity: "Ω_GW needed to explain them", status: Status::Disputed, source: "https://arxiv.org/abs/2110.13319" },
];

#[cfg(test)]
mod tests {
    use super::*;
    use crate::strain::{Convention, hc_from_omega};

    #[test]
    fn ceiling_in_strain() {
        let hc = hc_from_omega(CEILING_BBN_CMB, 1e9, 1.0, Convention::OneSided);
        assert!((hc / 1.33e-30 - 1.0).abs() < 0.01, "{hc:e}");
        // Ω h² < 5.6×10⁻⁶ ΔN_eff: ΔN_eff = 0.2 gives the quoted 1.1×10⁻⁶ (1.12 before rounding).
        assert!((ceiling_from_delta_neff(0.2) / CEILING_BBN_CMB - 1.0).abs() < 0.02);
    }

    #[test]
    fn photon_conversion_fixtures() {
        assert!((gertsenshtein_probability(1.0, 1.0) / 8.263e-38 - 1.0).abs() < 1e-3);
        assert!((gertsenshtein_probability(5.3, 106.0) / 2.61e-32 - 1.0).abs() < 0.01);
        assert!((stochastic_photon_rate(1.0, 1.0, 1.0, 1.0, 1.0, 1.0) / 3.954e31 - 1.0).abs() < 1e-3);
        // One photon per second: Ejlli et al.'s published prefactor 1.6×10⁻¹⁶.
        assert!((stochastic_strain_for_rate(1.0, 1.0, 1.0, 1.0, 1.0, 1.0) / 1.590e-16 - 1.0).abs() < 1e-3);
        assert!((coherent_photon_rate(1.0, 1.0, 1.0, 1.0, 1.0, 1.0) / 1.977e31 - 1.0).abs() < 1e-3);
    }

    #[test]
    fn cavity_and_sources() {
        // ADMX-like cavity: 0.8 GHz, η = 0.1, Q = 8×10⁴, 7.5 T, 0.51 m, 0.6 K. The Review quotes
        // 9.7×10⁻²³ two-sided; one-sided is √2 larger.
        let asd = cavity_noise(0.6, 0.8e9, 0.1, 8e4, 7.5, 0.51).sqrt();
        assert!((asd / (9.74e-23 * 2f64.sqrt()) - 1.0).abs() < 0.01, "{asd:e}");
        assert!((superradiance_frequency(1.0) / 4.836e14 - 1.0).abs() < 1e-3);
        assert!((phase_transition_shape(1.0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn registry_values_are_the_published_ones() {
        let value = |name: &str| REGISTRY.iter().find(|l| l.name.starts_with(name)).map(|l| l.value);
        assert_eq!(value("QUEST"), Some(3e-20));
        assert_eq!(value("Fermilab Holometer"), Some(1e-21)); // "better than 10⁻²¹/√Hz"
        assert_eq!(value("CAPP-12T"), Some(3.9e-21));
        assert_eq!(value("MAGE"), Some(1.35e-20)); // measured, not the 6.6×10⁻²¹ prediction
        assert_eq!(value("Light shining"), Some(6e-26));
        assert_eq!(value("Helioscopes"), Some(5e-28));
    }

    #[test]
    fn registry_sits_far_above_the_ceiling() {
        // Every published √S_h limit, converted to Ω h², is ≥ 10¹⁰ times the ceiling.
        for l in REGISTRY.iter().filter(|l| l.quantity.starts_with("√S_h")) {
            let f = (l.band.0 * l.band.1).sqrt();
            let omega = crate::strain::omega_from_sh(l.value * l.value, f, 1.0);
            assert!(omega > 1e10 * CEILING_BBN_CMB, "{}: {omega:e}", l.name);
        }
    }

    #[test]
    fn plasma_cutoff() {
        assert!((plasma_frequency(1e6) / 8980.0 - 1.0).abs() < 1e-3);
        assert!((plasma_frequency(3e4) / 1555.15 - 1.0).abs() < 1e-5);
    }

    #[test]
    fn plasma_refraction() {
        // X-band: the critical density is 1.24e18 m^-3, where the index reaches zero.
        let nc = critical_density(10e9);
        assert!((nc / 1.240_44e18 - 1.0).abs() < 1e-5);
        assert!((plasma_frequency(nc) / 10e9 - 1.0).abs() < 1e-12);
        assert!(plasma_index(nc, 10e9).is_none() && plasma_index(2.0 * nc, 10e9).is_none());
        // Underdense: a usable graded index just below the cutoff.
        for (n_e, n) in [(1e17, 0.9588), (5e17, 0.7726), (1e18, 0.4403)] {
            assert!((plasma_index(n_e, 10e9).unwrap() - n).abs() < 1e-4, "{n_e}");
        }
    }
}
