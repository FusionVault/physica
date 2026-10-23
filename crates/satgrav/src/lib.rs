//! Instrument noise models for satellite gravimetry: the amplitude spectral densities (ASDs) of
//! accelerometers and inter-satellite laser ranging, and the line-of-sight noise of a GRACE-type
//! satellite pair.
//!
//! The models are equations 1–8 of Zingerle et al. (2024), "The Benefits of Future Quantum
//! Accelerometers for Satellite Gravimetry", *Earth and Space Science* 11(9),
//! doi:10.1029/2024EA003630 (CC BY 4.0). Frequencies are in Hz.
//!
//! ```
//! use satgrav::{Accelerometer, Ranging, pair_noise};
//!
//! // GRACE-FO class at 1 mHz: SuperSTAR accelerometers plus the GRACE-FO laser interferometer.
//! let n = pair_noise(Accelerometer::SuperStar, Ranging::GraceFo, 1e-3);
//! assert!((n - 3.5e-10).abs() < 0.05e-10);
//! ```

mod atom;
mod errors;

pub use atom::{AtomInterferometer, K_EFF_RB87};
pub use errors::{gradiometry_errors, radial_gradient_transfer, resolution};

use core::f64::consts::{SQRT_2, TAU};

/// An accelerometer noise model; [`asd`](Self::asd) is in m s⁻² Hz⁻¹ᐟ².
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Accelerometer {
    /// Scaled LISA-Pathfinder-type sensor, drag-free at 350 km (eq. 1).
    Alvarez350,
    /// The same sensor without drag-free control at 500 km (eq. 2).
    Alvarez500,
    /// ONERA "Dalin" concept (eq. 3).
    Dalin,
    /// GRACE / GRACE-FO SuperSTAR electrostatic accelerometer (eq. 4).
    SuperStar,
    /// NGGM-class MicroSTAR electrostatic accelerometer (eq. 5).
    MicroStar,
    /// A cold-atom interferometer with a white noise floor, m s⁻² Hz⁻¹ᐟ²: Zingerle's CAI11 is
    /// `1e-11` ("realistic"), CAI12 `1e-12` ("future").
    AtomInterferometer {
        /// White ASD level.
        white: f64,
    },
}

impl Accelerometer {
    /// Amplitude spectral density at frequency `f`, m s⁻² Hz⁻¹ᐟ².
    pub fn asd(&self, f: f64) -> f64 {
        match *self {
            Accelerometer::Alvarez350 => {
                3e-8 * f * f + 4e-13 * (1.0 + 7e-4 / f + (3e-4 / f).powi(2)).sqrt()
            }
            Accelerometer::Alvarez500 => 3e-8 * f * f + 5e-13 * (1.0 + (1.0 / f).powf(2.0 / 3.0)).sqrt(),
            Accelerometer::Dalin => 2.4e-13 * f.powf(-0.25) + 1e-18 / (f * f) + 8e-11 * f * f,
            Accelerometer::SuperStar => 1e-10 * (1.0 + 0.005 / f).sqrt(),
            Accelerometer::MicroStar => {
                1e-11 / SQRT_2
                    * ((1e-3 / f).powi(4) / ((1e-5 / f).powi(4) + 1.0) + 1.0 + (f / 0.1).powi(4)).sqrt()
            }
            Accelerometer::AtomInterferometer { white } => white,
        }
    }
}

/// An inter-satellite laser ranging noise model; [`asd`](Self::asd) is in m Hz⁻¹ᐟ².
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Ranging {
    /// The GRACE-FO laser ranging interferometer as flown (eq. 6).
    GraceFo,
    /// Projected next-generation LRI, 2030 (eq. 7).
    Lri2030,
    /// Projected next-generation LRI, 2033 (eq. 8).
    Lri2033,
}

impl Ranging {
    /// Ranging ASD at frequency `f`, m Hz⁻¹ᐟ².
    pub fn asd(&self, f: f64) -> f64 {
        let (white, drift) = match self {
            Ranging::GraceFo => (2.2e-10, 1e-12),
            Ranging::Lri2030 => (2.2e-10, 1e-13),
            Ranging::Lri2033 => (1.1e-10, 5e-14),
        };
        white / f.sqrt() + drift / (f * f)
    }

    /// The ranging noise expressed as range acceleration, `(2πf)²·asd(f)`, m s⁻² Hz⁻¹ᐟ².
    pub fn as_acceleration(&self, f: f64) -> f64 {
        (TAU * f).powi(2) * self.asd(f)
    }
}

/// Line-of-sight noise of a satellite pair, m s⁻² Hz⁻¹ᐟ²: the two accelerometers (×√2) and
/// the ranging (as acceleration) added in quadrature. This is the "product noise" of
/// Zingerle et al. (2024), Fig. 1b.
pub fn pair_noise(accelerometer: Accelerometer, ranging: Ranging, f: f64) -> f64 {
    (SQRT_2 * accelerometer.asd(f)).hypot(ranging.as_acceleration(f))
}

/// Noise density of two accelerometers combined by inverse-variance weighting at each frequency,
/// `1/s² = 1/s_a² + 1/s_b²`: a stand-in for an electrostatic + cold-atom hybrid. Zingerle et al.
/// (2024) "optimally combine" the two but do not publish how; real hybrids (Zahzam et al. 2022)
/// use bias loops or Kalman filters and also suffer from the CAI's dead time.
pub fn hybrid_asd(a: Accelerometer, b: Accelerometer, f: f64) -> f64 {
    let (x, y) = (a.asd(f), b.asd(f));
    x * y / x.hypot(y)
}

/// Spherical-harmonic degree that a signal at frequency `f` maps to along an orbit of period
/// `orbit_period` seconds, `n ≈ f·T` (one cycle per revolution is degree 1).
pub fn degree_at_frequency(f: f64, orbit_period: f64) -> f64 {
    f * orbit_period
}

/// Period of a circular orbit at altitude `h` metres above a sphere of radius `radius` with
/// gravitational constant `gm`, seconds.
pub fn orbit_period(gm: f64, radius: f64, h: f64) -> f64 {
    TAU * ((radius + h).powi(3) / gm).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, rel: f64) -> bool {
        (a - b).abs() <= rel * b.abs()
    }

    /// Values read off / recomputed from Zingerle et al. (2024) Fig. 1.
    #[test]
    fn zingerle_figure_1_checkpoints() {
        assert!(close(Ranging::GraceFo.as_acceleration(0.1), 3.1e-10, 0.02));
        assert!(close(Ranging::GraceFo.as_acceleration(1e-3), 4.0e-11, 0.02));
        assert!(close(pair_noise(Accelerometer::SuperStar, Ranging::GraceFo, 1e-3), 3.5e-10, 0.02));
        let cai12 = Accelerometer::AtomInterferometer { white: 1e-12 };
        assert!(close(pair_noise(cai12, Ranging::Lri2033, 1e-3), 2.5e-12, 0.03));
    }

    #[test]
    fn white_and_limiting_behaviour() {
        let cai = Accelerometer::AtomInterferometer { white: 1e-11 };
        assert_eq!(cai.asd(1e-4), cai.asd(1e-1));
        // SuperSTAR flattens to 1e-10 above its 5 mHz corner.
        assert!(close(Accelerometer::SuperStar.asd(10.0), 1e-10, 1e-3));
        // MicroSTAR's floor is 1e-11/√2 between 1 mHz and 0.1 Hz.
        assert!(close(Accelerometer::MicroStar.asd(1e-2), 1e-11 / SQRT_2, 0.02));
    }

    #[test]
    fn hybrid_takes_the_better_sensor_per_frequency() {
        let cai = Accelerometer::AtomInterferometer { white: 1e-11 };
        assert!(close(hybrid_asd(cai, cai, 1e-3), 1e-11 / SQRT_2, 1e-12));
        // Below SuperSTAR's corner the CAI dominates; well above it neither does.
        assert!(close(hybrid_asd(Accelerometer::SuperStar, cai, 1e-4), 1e-11, 0.01));
        assert!(hybrid_asd(Accelerometer::MicroStar, cai, 1e-2) < Accelerometer::MicroStar.asd(1e-2));
    }

    #[test]
    fn grace_orbit_maps_revolution_frequency_to_degree_one() {
        let t = orbit_period(3.986004418e14, 6_371_000.0, 463e3);
        assert!(close(t, 5620.0, 0.002), "{t}");
        assert!(close(degree_at_frequency(1.0 / t, t), 1.0, 1e-12));
    }
}
