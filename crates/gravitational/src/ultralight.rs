//! Ultralight bosonic dark matter as an oscillating signal: the Compton frequency of a very light
//! boson, and the detector band it falls in.
//!
//! A boson light enough to be a classical wave over galactic scales oscillates the local field at
//! its Compton frequency `f = mc²/h`. If it couples to the standard model, "constants" such as α
//! and the electron mass oscillate with it, which an atom interferometer or clock sees as a phase
//! modulation. The frequency band between LISA (below about 0.1 Hz) and the ground
//! interferometers (above about 10 Hz) is the gap mid-band detectors such as MAGIS-100 are built
//! for, and it maps to a definite mass window.
//!
//! ```
//! use gravitational::ultralight::{MIDBAND, compton_frequency, mass_from_frequency};
//!
//! // A 10⁻¹⁵ eV boson oscillates at 0.24 Hz, inside the mid-band.
//! assert!((compton_frequency(1e-15) - 0.2418).abs() < 1e-4);
//! let (lo, hi) = MIDBAND;
//! assert!((lo..hi).contains(&compton_frequency(1e-15)));
//! // The band's edges in mass: 4.1×10⁻¹⁷ to 4.1×10⁻¹⁴ eV.
//! assert!((mass_from_frequency(hi) / 4.136e-14 - 1.0).abs() < 1e-3);
//! ```

use metrology::{C, EV, H};

/// The frequency gap between LISA and the ground-based interferometers, Hz: what mid-band
/// detectors (MAGIS-100, AION, AEDGE) are designed to cover.
pub const MIDBAND: (f64, f64) = (0.01, 10.0);

/// Compton frequency `mc²/h` (Hz) of a boson of mass `m_ev` (eV/c²): the rate at which an
/// ultralight dark-matter field oscillates.
pub fn compton_frequency(m_ev: f64) -> f64 {
    m_ev * EV / H
}

/// The boson mass (eV/c²) whose Compton frequency is `f` (Hz): `hf/c²`, the inverse of
/// [`compton_frequency`].
pub fn mass_from_frequency(f: f64) -> f64 {
    H * f / EV
}

/// The de Broglie wavelength (m) of an ultralight boson of mass `m_ev` (eV/c²) moving at the
/// galactic virial speed `v` (m/s), `h/(mv)`: the field is coherent over this scale, so it
/// behaves as a classical wave only when it is much larger than the detector.
pub fn de_broglie_wavelength(m_ev: f64, v: f64) -> f64 {
    H / (m_ev * EV / (C * C) * v)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, rel: f64) -> bool {
        (a / b - 1.0).abs() < rel
    }

    #[test]
    fn the_midband_in_mass() {
        for (m, f) in [(1e-15, 0.241_8), (1e-14, 2.418), (1e-13, 24.18)] {
            assert!(close(compton_frequency(m), f, 1e-3), "{m}");
            assert!(close(mass_from_frequency(compton_frequency(m)), m, 1e-12));
        }
        // The LISA–LIGO gap maps to 4.1×10⁻¹⁷ … 4.1×10⁻¹⁴ eV.
        assert!(close(mass_from_frequency(MIDBAND.0), 4.136e-17, 1e-3));
        assert!(close(mass_from_frequency(MIDBAND.1), 4.136e-14, 1e-3));
    }

    #[test]
    fn coherence_over_the_galaxy() {
        // At 220 km/s a 10⁻¹⁵ eV field is coherent over 1.7×10¹² m, about 11 au: vastly larger
        // than any detector, so it acts as a uniform classical wave in the laboratory.
        let lambda = de_broglie_wavelength(1e-15, 220e3);
        assert!(close(lambda, 1.6895e12, 1e-3), "{lambda:e}");
        // Lighter means more coherent.
        assert!(de_broglie_wavelength(1e-16, 220e3) > lambda);
    }
}
