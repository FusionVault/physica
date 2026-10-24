//! The energy gravitational waves carry: the Planck luminosity `c⁵/G`, the strain a given
//! luminosity produces, the quadrupole flux of a circular binary and its multipole coefficients,
//! and the orbital decay of binary pulsars (Peters & Mathews 1963).
//!
//! `c⁵/G = 3.6×10⁵² W` is the natural unit of gravitational-wave power, not a bound: the
//! brightest events reach about a thousandth of it (GW150914 at 9.6×10⁻⁴, GW190521 at
//! 1.0×10⁻³). The pulsar formula is the oldest confirmed prediction here; B1913+16's decay
//! matches it to 0.2% ([Weisberg & Huang 2016](https://arxiv.org/abs/1606.02744)) and the
//! double pulsar's to 0.01% ([Kramer et al. 2021](https://arxiv.org/abs/2112.06795)).
//!
//! ```
//! use gravitational::energy::{PLANCK_LUMINOSITY, peters_pbdot};
//! use metrology::SOLAR_MASS;
//!
//! // GW150914 peaked at 3.5×10⁵⁶ erg/s: under a thousandth of c⁵/G.
//! assert!((3.5e49 / PLANCK_LUMINOSITY / 9.6465e-4 - 1.0).abs() < 1e-4);
//! // The Hulse–Taylor pulsar loses 2.40×10⁻¹² s of orbital period per second.
//! let pbdot = peters_pbdot(1.438 * SOLAR_MASS, 1.390 * SOLAR_MASS, 0.322997448918 * 86_400.0, 0.6171340);
//! assert!((pbdot / -2.40263e-12 - 1.0).abs() < 1e-3);
//! ```

use core::f64::consts::{PI, SQRT_2};

use metrology::{C, G};

/// The Planck luminosity `c⁵/G`, W. A gravitational-wave source radiating `ε c⁵/G` has
/// luminosity fraction `ε`.
pub const PLANCK_LUMINOSITY: f64 = C * C * C * C * C / G;

/// Amplitude `h₀` (sinusoidal, optimally oriented) at distance `r` (m) of a source radiating the
/// fraction `fraction` of [`PLANCK_LUMINOSITY`] at frequency `f` (Hz), from the flux
/// `c³ω²h₀²/(32πG)`: `h₀ = 2√2 √ε c/(2πf r)`. At one wavelength from a source at `c⁵/G` this is
/// 0.45: the strain cannot be of order one in the wave zone.
pub fn strain_from_luminosity(fraction: f64, f: f64, r: f64) -> f64 {
    2.0 * SQRT_2 * fraction.sqrt() * C / (2.0 * PI * f * r)
}

/// Quadrupole luminosity of a circular binary in units of [`PLANCK_LUMINOSITY`]:
/// `(32/5) η² x⁵`, with `eta` the symmetric mass ratio and `x = (GMω/c³)^(2/3)` the
/// post-Newtonian parameter (1/6 at the innermost stable orbit).
pub fn circular_flux(eta: f64, x: f64) -> f64 {
    32.0 / 5.0 * eta * eta * x.powi(5)
}

/// The energy-flux coefficients `(a_ℓ, b_ℓ)` of the mass and current multipoles of order `l ≥ 2`
/// ([Blanchet's Living Review](https://arxiv.org/abs/1310.1528), Eq. 68):
/// `a_ℓ = (ℓ+1)(ℓ+2)/((ℓ−1)ℓ ℓ!(2ℓ+1)!!)` and `b_ℓ = 4ℓ(ℓ+2)/((ℓ−1)(ℓ+1)!(2ℓ+1)!!)`.
/// `a₂ = 1/5` is the quadrupole formula's `G/(5c⁵)`.
pub fn energy_flux_coefficients(l: u32) -> (f64, f64) {
    let factorial = |n: u32| (1..=n).map(f64::from).product::<f64>();
    let double_factorial = |n: u32| (1..=n).rev().step_by(2).map(f64::from).product::<f64>();
    let l_f = f64::from(l);
    let a = (l_f + 1.0) * (l_f + 2.0) / ((l_f - 1.0) * l_f * factorial(l) * double_factorial(2 * l + 1));
    let b = 4.0 * l_f * (l_f + 2.0) / ((l_f - 1.0) * factorial(l + 1) * double_factorial(2 * l + 1));
    (a, b)
}

/// Peters' enhancement of the radiated power by eccentricity `e`:
/// `(1 + 73e²/24 + 37e⁴/96)/(1 − e²)^(7/2)`.
pub fn peters_enhancement(e: f64) -> f64 {
    let e2 = e * e;
    (1.0 + 73.0 / 24.0 * e2 + 37.0 / 96.0 * e2 * e2) / (1.0 - e2).powf(3.5)
}

/// Orbital-period derivative (dimensionless, s/s) of a binary of masses `mp` and `mc` (kg),
/// period `pb` (s) and eccentricity `e`, from gravitational radiation:
/// `−(192π/5)(G/c³)^(5/3)(P_b/2π)^(−5/3) f(e) m_p m_c (m_p+m_c)^(−1/3)`.
pub fn peters_pbdot(mp: f64, mc: f64, pb: f64, e: f64) -> f64 {
    -192.0 * PI / 5.0 * (G / C.powi(3)).powf(5.0 / 3.0) * (pb / (2.0 * PI)).powf(-5.0 / 3.0) * peters_enhancement(e) * mp * mc
        / (mp + mc).cbrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use metrology::{MPC, SOLAR_MASS};

    fn close(a: f64, b: f64, rel: f64) -> bool {
        (a / b - 1.0).abs() < rel
    }

    #[test]
    fn luminosity_and_strain() {
        assert!(close(PLANCK_LUMINOSITY, 3.62825e52, 1e-5) && close(PLANCK_LUMINOSITY / 4.0, 9.07064e51, 1e-5));
        assert!(close(3.7e49 / PLANCK_LUMINOSITY, 1.01977e-3, 1e-5)); // GW190521
        assert!(close(strain_from_luminosity(1e-3, 250.0, 410.0 * MPC), 1.34931e-21, 1e-5));
        assert!(close(strain_from_luminosity(1.0, 1.0, C), 0.450158, 1e-5));
        assert!(close(circular_flux(0.25, 1.0 / 6.0), 5.14403e-5, 1e-5));
    }

    #[test]
    fn multipole_coefficients() {
        for (l, a, b) in [(2, 1.0 / 5.0, 16.0 / 45.0), (3, 1.0 / 189.0, 1.0 / 84.0), (4, 1.0 / 9072.0, 4.0 / 14175.0)] {
            let (al, bl) = energy_flux_coefficients(l);
            assert!(close(al, a, 1e-12) && close(bl, b, 1e-12), "{l}: {al} {bl}");
        }
    }

    #[test]
    fn binary_pulsars() {
        assert!(close(peters_enhancement(0.6171340), 11.8568, 1e-5));
        let day = 86_400.0;
        // Masses quoted to 10⁻³: Weisberg & Huang's −2.40263(5)×10⁻¹² to that precision.
        let b1913 = peters_pbdot(1.438 * SOLAR_MASS, 1.390 * SOLAR_MASS, 0.322997448918 * day, 0.6171340);
        assert!(close(b1913, -2.40210e-12, 1e-5) && close(b1913, -2.40263e-12, 1e-3));
        let j0737 = peters_pbdot(1.338185 * SOLAR_MASS, 1.248868 * SOLAR_MASS, 0.1022515592973 * day, 0.087777023);
        assert!(close(j0737, -1.24781e-12, 1e-5));
    }
}
