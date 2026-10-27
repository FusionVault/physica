//! Compact-binary inspirals at leading (Newtonian) order: chirp mass, time to merger, frequency
//! drift, the innermost stable circular orbit, matched-filter signal-to-noise ratio, and how far
//! away a detector can see a binary ([Sathyaprakash & Schutz 2009](https://arxiv.org/abs/0903.0338);
//! [Chen et al. 2021](https://arxiv.org/abs/1709.08079)).
//!
//! ```
//! use gravitational::{curves, inspiral};
//! use metrology::{MPC, SOLAR_MASS};
//!
//! // Two 1.4 M☉ neutron stars and the Advanced LIGO design fit.
//! let (m1, m2) = (1.4 * SOLAR_MASS, 1.4 * SOLAR_MASS);
//! let horizon = inspiral::horizon(m1, m2, 8.0, curves::aligo, 10.0);
//! println!("visible to {:.0} Mpc", horizon / MPC);
//! assert!((400.0..600.0).contains(&(horizon / MPC)));
//! ```

use crate::integrate_log;
use metrology::{C, G};
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

/// The innermost stable circular orbit of a Kerr black hole whose spin `a` (`cJ/GM²`, from −1
/// to 1) is signed along the orbit's angular momentum, negative for a retrograde orbit:
/// `(r, E, L)` in units of `GM/c²`, `mc²` and `Gm M/c` (Bardeen, Press & Teukolsky 1972, as
/// printed in [Hofmann, Barausse & Rezzolla 2016](https://arxiv.org/abs/1605.01938)):
/// `r = 3 + Z₂ − sign(a)√((3 − Z₁)(3 + Z₁ + 2Z₂))`, `E = √(1 − 2/(3r))`,
/// `L = (2/(3√3))(1 + 2√(3r − 2))`. At `a = 0` it is [`f_isco`]'s 6M; it reaches 1M prograde
/// and 9M retrograde at `a = 1`.
pub fn kerr_isco(a: f64) -> (f64, f64, f64) {
    let z1 = 1.0 + (1.0 - a * a).cbrt() * ((1.0 + a).cbrt() + (1.0 - a).cbrt());
    let z2 = (3.0 * a * a + z1 * z1).sqrt();
    let r = 3.0 + z2 - a.signum() * ((3.0 - z1) * (3.0 + z1 + 2.0 * z2)).sqrt();
    (r, (1.0 - 2.0 / (3.0 * r)).sqrt(), 2.0 / (3.0 * 3f64.sqrt()) * (1.0 + 2.0 * (3.0 * r - 2.0).sqrt()))
}

/// The gravitational-wave frequency (Hz, twice the orbital frequency: the `l = m = 2` harmonic) at
/// the Kerr ISCO of [`kerr_isco`] for total mass `m` (kg) and signed spin `a`: `MΩ =
/// 1/(r_ISCO^(3/2) + a)` ([Isoyama et al. 2014](https://arxiv.org/abs/1404.6133), Eq. 3, after
/// Bardeen, Press & Teukolsky 1972) and `f = Ω/π = c³/(πGm(r^(3/2) + a))`. `MΩ` is `6^(−3/2)`
/// without spin ([`f_isco`]), 1/2 prograde and 1/26 retrograde at `a = 1`. This is the test-mass
/// geodesic value; the self-force shift is O(μ/M) and not included.
pub fn kerr_f_isco(m: f64, a: f64) -> f64 {
    let (r, _, _) = kerr_isco(a);
    C.powi(3) / (PI * G * m * (r.powf(1.5) + a))
}

/// The radius (units of `GM/c²`) of the circular photon orbit in the equatorial plane of a Kerr
/// black hole of signed spin `a`, as in [`kerr_isco`]: `2(1 + cos((2/3)arccos(−a)))`
/// (Bardeen, Press & Teukolsky 1972). 3M without spin; 1M prograde and 4M retrograde at `a = 1`.
pub fn kerr_photon_orbit(a: f64) -> f64 {
    2.0 * (1.0 + (2.0 / 3.0 * (-a).acos()).cos())
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
/// `|h̃(f)| = √(5/24) π^(−2/3) (c/d) (G M_c/c³)^(5/6) f^(−7/6)`, Hz⁻¹. With a sky-averaged
/// sensitivity such as [`lisa`](crate::curves::lisa) (Robson et al.), multiply by 4/5 for the
/// inclination average (Robson Eq. 16), or the SNR comes out too high.
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
/// Zero when the binary merges below `f_lo`.
pub fn horizon(m1: f64, m2: f64, threshold: f64, s_n: impl Fn(f64) -> f64, f_lo: f64) -> f64 {
    let chirp = chirp_mass(m1, m2);
    let f_hi = f_isco(m1 + m2);
    if f_lo >= f_hi {
        return 0.0;
    }
    snr(|f| amplitude(f, chirp, 1.0), s_n, f_lo, f_hi) / threshold
}

/// Energy spectrum of the inspiral, `dE/df = (π^(2/3)/3) G^(2/3) M_c^(5/3) f^(−1/3)`, J/Hz.
/// Divided by `4πd²` it is the spectral fluence at distance `d` for isotropic emission.
pub fn energy_spectrum(chirp: f64, f: f64) -> f64 {
    PI.powf(2.0 / 3.0) / 3.0 * G.powf(2.0 / 3.0) * chirp.powf(5.0 / 3.0) * f.powf(-1.0 / 3.0)
}

/// Horizon over range for a single interferometer (Chen et al. 2021): 2.264.
pub const HORIZON_TO_RANGE: f64 = 2.264;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::curves;
    use metrology::{MPC, SOLAR_MASS};

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
    fn kerr_orbits() {
        // The closed forms: 6M, √(8/9), 2√3 without spin; 1M, 1/√3, 2/√3 prograde and 9M,
        // √(25/27), 22/(3√3) retrograde at a = 1.
        let close = |a: f64, b: f64| (a - b).abs() < 1e-14 * b.abs().max(1.0);
        let cases = [(0.0, 6.0, (8.0f64 / 9.0).sqrt(), 2.0 * 3f64.sqrt()), (1.0, 1.0, 1.0 / 3f64.sqrt(), 2.0 / 3f64.sqrt()), (-1.0, 9.0, (25.0f64 / 27.0).sqrt(), 22.0 / (3.0 * 3f64.sqrt()))];
        for (a, r, e, l) in cases {
            let (ri, ei, li) = kerr_isco(a);
            assert!(close(ri, r) && close(ei, e) && close(li, l), "a = {a}: {ri} {ei} {li}");
        }
        for a in [-0.9, -0.5, 0.3, 0.7, 0.99] {
            // Marginal stability, r² − 6r + 8a√r − 3a² = 0, and the energy and angular momentum of
            // a circular orbit at that radius (Bardeen, Press & Teukolsky 1972, Eqs. 2.12–2.13).
            let (r, e, l) = kerr_isco(a);
            assert!((r * r - 6.0 * r + 8.0 * a * r.sqrt() - 3.0 * a * a).abs() < 1e-13 * r * r, "a = {a}");
            let d = r.powf(0.75) * (r.powf(1.5) - 3.0 * r.sqrt() + 2.0 * a).sqrt();
            assert!(close(e, (r.powf(1.5) - 2.0 * r.sqrt() + a) / d) && close(l, (r * r - 2.0 * a * r.sqrt() + a * a) / d), "a = {a}");
            // The photon orbit is where that denominator vanishes: r^(3/2) − 3r^(1/2) + 2a = 0.
            let p = kerr_photon_orbit(a);
            assert!((p.powf(1.5) - 3.0 * p.sqrt() + 2.0 * a).abs() < 1e-14, "photon orbit at a = {a}");
        }
        assert!(close(kerr_photon_orbit(0.0), 3.0) && close(kerr_photon_orbit(1.0), 1.0) && close(kerr_photon_orbit(-1.0), 4.0));
        // f_isco is the a = 0 case of kerr_f_isco.
        assert!(close(f_isco(SOLAR_MASS), kerr_f_isco(SOLAR_MASS, 0.0)));
    }

    #[test]
    fn kerr_isco_frequency() {
        // Isoyama et al. 2014, Eq. 3: MΩ = 6^(−3/2), 1/2 and 1/26 at a = 0, +1 and −1.
        let m_omega = |a: f64| kerr_f_isco(SOLAR_MASS, a) * PI * G * SOLAR_MASS / C.powi(3);
        for (a, exact) in [(0.0, 6f64.powf(-1.5)), (1.0, 0.5), (-1.0, 1.0 / 26.0)] {
            assert!((m_omega(a) / exact - 1.0).abs() < 1e-15, "a = {a}: {}", m_omega(a));
        }
        // The research script's a = 0.9: r = 2.3208830418, MΩ = 0.225441708835.
        assert!((m_omega(0.9) - 0.225441708835).abs() < 1e-12, "{}", m_omega(0.9));
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

    #[test]
    fn heavy_binaries_that_merge_below_the_band_have_no_horizon() {
        // 500 + 500 M☉ merges at 4.4 Hz, below a 10 Hz cut: horizon 0, not NaN.
        let m = 500.0 * SOLAR_MASS;
        assert_eq!(horizon(m, m, 8.0, crate::curves::aligo, 10.0), 0.0);
    }

    #[test]
    fn gw170817_spectral_fluence() {
        let fluence = energy_spectrum(1.188 * SOLAR_MASS, 900.0) / (4.0 * PI * (40.0 * MPC).powi(2));
        assert!((fluence / 2.66691e-7 - 1.0).abs() < 1e-5);
    }
}
