//! How gravitational waves travel: their speed, a graviton mass, gravitational lensing, and
//! absorption by matter.
//!
//! GW170817 arrived 1.74 s before its gamma-ray burst after 40 Mpc, which pins the speed of
//! gravity to c within −3×10⁻¹⁵ … +7×10⁻¹⁶ ([LVK, Fermi & INTEGRAL 2017](https://arxiv.org/abs/1710.05834))
//! and removed most relativistic dark-energy and MOND theories in one event. A massive graviton
//! would make waves dispersive, `1 − v/c = ½(c/(λ_g f))²` with `λ_g = h/(m_g c)`
//! ([Will 2014](https://arxiv.org/abs/1403.7377)). Matter is nearly transparent: a viscous
//! fluid damps waves at `16πGη/c²` ([Baym, Patil & Pethick 2017](https://arxiv.org/abs/1707.05192)),
//! an absorption length of 2.6×10¹¹ Gpc in water.
//!
//! ```
//! use gravitational::propagation::{GRAVITON_MASS_GW, compton_wavelength, speed_bounds};
//! use metrology::MPC;
//!
//! let (upper, lower) = speed_bounds(1.74, 26.0 * MPC, 10.0);
//! assert!((upper / 6.502e-16 - 1.0).abs() < 1e-3 && (lower / -3.0866e-15 - 1.0).abs() < 1e-3);
//! // The graviton's Compton wavelength is over 6×10¹³ km.
//! assert!(compton_wavelength(GRAVITON_MASS_GW) > 6.4e16);
//! ```

use core::f64::consts::PI;

use metrology::{C, EV, G, H};

/// Upper bound (90%) on the graviton mass from gravitational-wave dispersion, eV/c²: GWTC-4.0,
/// 91 events ([arXiv:2603.19020](https://arxiv.org/abs/2603.19020)). Earlier quotes differ by
/// catalogue version: GWTC-3 gave 1.27×10⁻²³ (v1) and 2.42×10⁻²³ (v3).
pub const GRAVITON_MASS_GW: f64 = 1.92e-23;
/// Upper bound (99.7%) on the graviton mass from planetary ephemerides (Mariani et al. 2023, as
/// cited in [arXiv:2603.19020](https://arxiv.org/abs/2603.19020)), eV/c².
pub const GRAVITON_MASS_EPHEMERIS: f64 = 1.01e-24;
/// The fraction of primordial tensor power that free-streaming neutrinos remove from modes
/// entering the horizon in the radiation era ([Weinberg 2004](https://arxiv.org/abs/astro-ph/0306304));
/// the amplitude falls by `√(1 − 0.356) = 0.80`.
pub const NEUTRINO_DAMPING: f64 = 0.356;

/// Bounds on `(v_GW − c)/c` from a gravitational wave arriving `delay` seconds before its
/// electromagnetic counterpart from distance `d` (m), when the light may have left up to
/// `max_emission` seconds after the waves: `(delay·c/d, (delay − max_emission)·c/d)`.
pub fn speed_bounds(delay: f64, d: f64, max_emission: f64) -> (f64, f64) {
    (delay * C / d, (delay - max_emission) * C / d)
}

/// Compton wavelength `λ_g = hc/(m_g c²)` (m) of a graviton of mass `m_ev` (eV/c²).
pub fn compton_wavelength(m_ev: f64) -> f64 {
    H * C / (m_ev * EV)
}

/// The speed deficit `1 − v/c = ½(m_g c²/(hf))²` of a wave of frequency `f` (Hz) for graviton
/// mass `m_ev` (eV/c²).
pub fn speed_deficit(m_ev: f64, f: f64) -> f64 {
    0.5 * (m_ev * EV / (H * f)).powi(2)
}

/// Extra travel time (s) over distance `d` (m) of a wave of frequency `f` for graviton mass
/// `m_ev`: `(d/c)(1 − v/c)`, ignoring cosmological expansion.
pub fn massive_graviton_delay(m_ev: f64, f: f64, d: f64) -> f64 {
    d / C * speed_deficit(m_ev, f)
}

/// The wave-optics lensing parameter `w = 8πGM_Lz f/c³` for a point lens of redshifted mass
/// `m_lens` (kg) at frequency `f` ([Takahashi & Nakamura 2003](https://arxiv.org/abs/astro-ph/0305055)).
/// Geometrical optics holds for `w ≫ 1`; below 1 the wave diffracts round the lens.
pub fn lens_parameter(m_lens: f64, f: f64) -> f64 {
    8.0 * PI * G * m_lens * f / C.powi(3)
}

/// The two images of a point lens in geometrical optics, at source offset `y` (in Einstein radii):
/// magnifications `μ± = ½ ± (y²+2)/(2y√(y²+4))` and their time delay (s) for redshifted lens mass
/// `m_lens` (kg), `(4GM/c³)[y√(y²+4)/2 + ln((√(y²+4)+y)/(√(y²+4)−y))]`.
pub fn point_lens_geometric(y: f64, m_lens: f64) -> (f64, f64, f64) {
    let s = (y * y + 4.0).sqrt();
    let mu = (y * y + 2.0) / (2.0 * y * s);
    let delay = 4.0 * G * m_lens / C.powi(3) * (y * s / 2.0 + ((s + y) / (s - y)).ln());
    (0.5 + mu, 0.5 - mu, delay)
}

/// Peak wave-optics amplification `|F(y = 0)|² = πw/(1 − e^(−πw))` of a perfectly aligned point
/// lens with parameter `w` ([`lens_parameter`]).
pub fn wave_optics_peak(w: f64) -> f64 {
    PI * w / -(-PI * w).exp_m1()
}

/// Damping rate (s⁻¹) of a wave's energy in a fluid of shear viscosity `eta` (Pa s),
/// `γ = 16πGη/c²`. `None` when `omega·tau > 1`: the formula needs the wave's angular frequency
/// `omega` to be slow next to the fluid's relaxation time `tau` (s), so it does not apply to,
/// say, a mantle viscosity at hertz frequencies.
pub fn viscous_damping(eta: f64, omega: f64, tau: f64) -> Option<f64> {
    (omega * tau <= 1.0).then(|| 16.0 * PI * G * eta / (C * C))
}

#[cfg(test)]
mod tests {
    use super::*;
    use metrology::{MPC, SOLAR_MASS};

    fn close(a: f64, b: f64, rel: f64) -> bool {
        (a / b - 1.0).abs() < rel
    }

    #[test]
    fn graviton_mass() {
        assert!(close(compton_wavelength(GRAVITON_MASS_GW), 6.45751e16, 1e-5));
        assert!(close(compton_wavelength(2.42e-23), 5.12331e16, 1e-5));
        assert!(close(compton_wavelength(GRAVITON_MASS_EPHEMERIS), 1.22757e18, 1e-5));
        // INPOP19a: 3.16×10⁻²³ eV is λ_g > 3.93×10¹³ km (Bernus et al. 2020).
        assert!(close(compton_wavelength(3.16e-23), 3.93e16, 1e-2));
        assert!(close(speed_deficit(GRAVITON_MASS_GW, 100.0), 1.07766e-21, 1e-5));
        assert!(close(massive_graviton_delay(GRAVITON_MASS_GW, 100.0, 1e3 * MPC), 1.10920e-4, 1e-5));
        // The ephemeris Yukawa factor e^(−r/λ_g) at 1 au: 1 − 2.3166×10⁻⁶.
        let au = 1.495_978_707e11;
        assert!(((-au / compton_wavelength(GRAVITON_MASS_GW)).exp_m1() / -2.3166e-6 - 1.0).abs() < 1e-4);
    }

    #[test]
    fn lensing() {
        assert!(close(lens_parameter(100.0 * SOLAR_MASS, 100.0), 1.23791, 1e-5));
        assert!(close(1.0 / lens_parameter(SOLAR_MASS, 100.0), 80.7813, 1e-5));
        let (plus, minus, delay) = point_lens_geometric(1.0, 1e8 * SOLAR_MASS);
        assert!(close(plus, 1.17082, 1e-5) && close(minus, -0.17082, 1e-4) && close(delay, 4098.91, 1e-5));
        assert!(close(wave_optics_peak(1.0), 3.28348, 1e-5));
    }

    #[test]
    fn absorption() {
        let water = viscous_damping(1e-3, 1e3, 1e-12).unwrap();
        assert!(close(water, 3.7328e-29, 1e-4));
        assert!(close(C / water / (1e3 * MPC), 2.60277e11, 1e-5));
        assert!(close(viscous_damping(1.8e-5, 1e3, 1e-10).unwrap(), 6.71903e-31, 1e-5));
        // Mantle: η = 10²¹ Pa s relaxes over ~10¹⁰ s; at 1 Hz the formula must refuse.
        assert!(viscous_damping(1e21, 2.0 * PI, 1e10).is_none());
        assert!(close((1.0 - NEUTRINO_DAMPING).sqrt(), 0.802496, 1e-6));
    }
}
