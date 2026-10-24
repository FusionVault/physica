//! Compact-binary inspirals at leading (Newtonian) order: chirp mass, time to merger, frequency
//! drift, the innermost stable circular orbit, matched-filter signal-to-noise ratio, and how far
//! away a detector can see a binary ([Sathyaprakash & Schutz 2009](https://arxiv.org/abs/0903.0338);
//! [Chen et al. 2021](https://arxiv.org/abs/1709.08079)).
//!
//! ```
//! use gwsens::{SOLAR_MASS, MPC, curves, inspiral};
//!
//! // Two 1.4 M☉ neutron stars and the Advanced LIGO design fit.
//! let (m1, m2) = (1.4 * SOLAR_MASS, 1.4 * SOLAR_MASS);
//! let horizon = inspiral::horizon(m1, m2, 8.0, curves::aligo, 10.0);
//! println!("visible to {:.0} Mpc", horizon / MPC);
//! assert!((400.0..600.0).contains(&(horizon / MPC)));
//! ```

use crate::{C, G, integrate_log};
use core::f64::consts::PI;

/// Chirp mass `(m₁m₂)^(3/5)/(m₁+m₂)^(1/5)`, kg.
pub fn chirp_mass(m1: f64, m2: f64) -> f64 {
    (m1 * m2).powf(0.6) / (m1 + m2).powf(0.2)
}

/// Gravitational-wave frequency at the innermost stable circular orbit of total mass `m` (kg):
/// `c³/(6^(3/2) π G m)`, Hz.
pub fn f_isco(m: f64) -> f64 {
    C.powi(3) / (6f64.powf(1.5) * PI * G * m)
}

/// Time from gravitational-wave frequency `f` to merger, leading order:
/// `τ₀ = (5/256)(πf)^(−8/3)(G M_c/c³)^(−5/3)`, s.
pub fn time_to_merger(f: f64, chirp: f64) -> f64 {
    5.0 / 256.0 * (PI * f).powf(-8.0 / 3.0) * (G * chirp / C.powi(3)).powf(-5.0 / 3.0)
}

/// Frequency drift `ḟ = (96/5) π^(8/3) (G M_c/c³)^(5/3) f^(11/3)`, Hz/s.
pub fn frequency_drift(f: f64, chirp: f64) -> f64 {
    96.0 / 5.0 * PI.powf(8.0 / 3.0) * (G * chirp / C.powi(3)).powf(5.0 / 3.0) * f.powf(11.0 / 3.0)
}

/// Fourier amplitude of an optimally placed and oriented inspiral at distance `d` (m):
/// `|h̃(f)| = √(5/24) π^(−2/3) (c/d) (G M_c/c³)^(5/6) f^(−7/6)`, Hz⁻¹.
pub fn amplitude(f: f64, chirp: f64, d: f64) -> f64 {
    (5.0f64 / 24.0).sqrt() * PI.powf(-2.0 / 3.0) * (C / d) * (G * chirp / C.powi(3)).powf(5.0 / 6.0) * f.powf(-7.0 / 6.0)
}

/// Optimal matched-filter SNR `ρ = √(4 ∫ |h̃|²/S_n df)` between `f_lo` and `f_hi`, for a signal
/// amplitude `h(f)` and one-sided noise `s_n(f)`.
pub fn snr(h: impl Fn(f64) -> f64, s_n: impl Fn(f64) -> f64, f_lo: f64, f_hi: f64) -> f64 {
    (4.0 * integrate_log(f_lo, f_hi, 4000, |f| { let a = h(f); a * a / s_n(f) })).sqrt()
}

/// Horizon distance (m): the distance at which an optimally placed and oriented binary of masses
/// `m1`, `m2` (kg) reaches SNR `threshold`, integrating from `f_lo` to the ISCO frequency.
/// The sky- and orientation-averaged range is [`HORIZON_TO_RANGE`] times smaller.
pub fn horizon(m1: f64, m2: f64, threshold: f64, s_n: impl Fn(f64) -> f64, f_lo: f64) -> f64 {
    let chirp = chirp_mass(m1, m2);
    snr(|f| amplitude(f, chirp, 1.0), s_n, f_lo, f_isco(m1 + m2)) / threshold
}

/// Horizon over range for a single interferometer (Chen et al. 2021): 2.264.
pub const HORIZON_TO_RANGE: f64 = 2.264;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{MPC, SOLAR_MASS, curves};

    #[test]
    fn binary_neutron_star_numbers() {
        let m = 1.4 * SOLAR_MASS;
        let mc = chirp_mass(m, m);
        assert!((mc / SOLAR_MASS / 1.21877 - 1.0).abs() < 1e-5, "{}", mc / SOLAR_MASS);
        assert!((time_to_merger(20.0, mc) / 157.9 - 1.0).abs() < 1e-3, "{}", time_to_merger(20.0, mc));
        assert!((frequency_drift(100.0, mc) / 17.36 - 1.0).abs() < 1e-3);
        assert!((f_isco(2.0 * m) / 1570.4 - 1.0).abs() < 1e-4);
        assert!((f_isco(SOLAR_MASS) / 4397.0 - 1.0).abs() < 1e-3);
    }

    #[test]
    fn snr_scales_inversely_with_distance() {
        let mc = chirp_mass(30.0 * SOLAR_MASS, 30.0 * SOLAR_MASS);
        let near = snr(|f| amplitude(f, mc, 400.0 * MPC), curves::aligo, 20.0, 70.0);
        let far = snr(|f| amplitude(f, mc, 800.0 * MPC), curves::aligo, 20.0, 70.0);
        assert!((near / far - 2.0).abs() < 1e-9);
    }

    #[test]
    fn horizon_with_the_analytic_fit() {
        // The design fit gives 496.5 Mpc for 1.4 + 1.4 M☉ (a tabulated O4 curve gives 449 Mpc in
        // Chen et al.; the fit is for quick estimates).
        let d = horizon(1.4 * SOLAR_MASS, 1.4 * SOLAR_MASS, 8.0, curves::aligo, 20.0) / MPC;
        assert!((d / 496.5 - 1.0).abs() < 1e-3, "{d}");
    }
}
