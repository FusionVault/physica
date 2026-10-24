//! Gravitons: how many quanta a classical wave holds, and whether a resonant bar could absorb
//! one ([Tobar et al. 2024](https://arxiv.org/abs/2308.15440); [Carney, Domcke & Rodd 2024](https://arxiv.org/abs/2308.12988)).
//!
//! A detectable wave is overwhelmingly classical: GW150914's strain at 150 Hz carries about
//! 4×10³⁶ gravitons per reduced wavelength squared per radian of phase. Yet a bar can still
//! register single quanta, because each absorbed graviton is one jump of a phonon mode.
//! GW170817 would have deposited about 0.19 ħω in an optimally oriented LSU-class bar
//! ([Aguiar 2011](https://arxiv.org/abs/1009.1138)'s cross-section): the quantum regime, with a
//! bar of the right mass, is in reach. Whether such a click would prove the graviton exists is
//! contested.
//!
//! ```
//! use gravitational::quanta::graviton_count;
//!
//! assert!((graviton_count(1e-21, 150.0) / 3.85283e36 - 1.0).abs() < 1e-5);
//! ```

use core::f64::consts::PI;

use metrology::{C, G, HBAR};

/// Gravitons per reduced wavelength squared per radian of phase in a wave of amplitude `h0` at
/// frequency `f`: the flux `c³ω²h₀²/(32πG)` times `(c/ω)²/ω` over `ħω`, i.e. `h₀²c⁵/(32πGħω²)`.
pub fn graviton_count(h0: f64, f: f64) -> f64 {
    let w = 2.0 * PI * f;
    h0 * h0 * C.powi(5) / (32.0 * PI * G * HBAR * w * w)
}

/// Mean number of gravitons in one cubic wavelength of a wave of amplitude `h` at frequency `f`:
/// the energy density `c²h²ω²/(32πG)` over `ħω`, times `(c/f)³`. Below one the wave is quantum;
/// that needs, for instance, `h = 10⁻³⁰` above 4.6×10¹² Hz.
pub fn occupation_per_cell(h: f64, f: f64) -> f64 {
    let w = 2.0 * PI * f;
    C * C * h * h * w * w / (32.0 * PI * G) / (HBAR * w) * (C / f).powi(3)
}

/// Rate (Hz) at which a wave of amplitude `h` stimulates single-graviton absorption in mode `l` of
/// a bar of mass `m` (kg) and sound speed `vs` (m/s), on resonance: `v_s²Mh²/(4l²π³ħ)`.
pub fn stimulated_rate(vs: f64, m: f64, h: f64, l: f64) -> f64 {
    vs * vs * m * h * h / (4.0 * l * l * PI.powi(3) * HBAR)
}

/// Spontaneous graviton emission rate (Hz) of a bar of density `rho` (kg/m³), sound speed `vs`,
/// radius `r` and length `l` (m) in its first excited mode: `8πGρv_s⁴R²/(Lc⁵)`.
pub fn spontaneous_rate(rho: f64, vs: f64, r: f64, l: f64) -> f64 {
    8.0 * PI * G * rho * vs.powi(4) * r * r / (l * C.powi(5))
}

/// The bar mass (kg) that absorbs on average one graviton from an inspiral of chirp mass `chirp`
/// (kg) whose amplitude is `h0` at frequency `f`, for sound speed `vs`:
/// `(24π²/5)(ħ/(h₀²v_s²))(GM_c/2c³)^(5/3) ω^(8/3)`.
pub fn optimal_bar_mass(vs: f64, h0: f64, chirp: f64, f: f64) -> f64 {
    let w = 2.0 * PI * f;
    24.0 * PI * PI / 5.0 * HBAR / (h0 * h0 * vs * vs) * (G * chirp / (2.0 * C.powi(3))).powf(5.0 / 3.0) * w.powf(8.0 / 3.0)
}

/// Frequency-integrated absorption cross-section (m² Hz) of mode `n` of a cylindrical bar of mass
/// `m` (kg), sound speed `vs`, length `l`, radius `r` (m) and Poisson ratio `poisson`, for a wave
/// arriving at angle `theta` to the axis ([Aguiar 2011](https://arxiv.org/abs/1009.1138), Eq. 4.1):
/// `(8/π)(M/n²)(G/c)(v_s/c)² sin⁴θ [1 + σ(1 − 2σ)(nπR/L)²/2]`. Times a spectral fluence
/// (J m⁻² Hz⁻¹, e.g. [`crate::inspiral::energy_spectrum`]) it gives the energy absorbed.
pub fn bar_cross_section(m: f64, vs: f64, l: f64, r: f64, poisson: f64, n: f64, theta: f64) -> f64 {
    let correction = 1.0 + poisson * (1.0 - 2.0 * poisson) * (n * PI * r / l).powi(2) / 2.0;
    8.0 / PI * m / (n * n) * (G / C) * (vs / C).powi(2) * theta.sin().powi(4) * correction
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::inspiral::energy_spectrum;
    use core::f64::consts::FRAC_PI_2;
    use metrology::{MPC, SOLAR_MASS};

    fn close(a: f64, b: f64, rel: f64) -> bool {
        (a / b - 1.0).abs() < rel
    }

    #[test]
    fn counting_gravitons() {
        // Tobar's ≈ 4×10³⁶ and Carney's ≈ 2×10³⁵.
        assert!(close(graviton_count(1e-21, 150.0), 4e36, 0.05));
        assert!(close(occupation_per_cell(1e-22, 1e3), 2.15031e35, 1e-5));
        assert!(close(occupation_per_cell(7e-28, 1e18), 1.05365e-5, 1e-5));
        assert!(close(occupation_per_cell(1e-30, 4.63715e12), 1.0, 1e-5));
    }

    #[test]
    fn tobar_rates_and_masses() {
        assert!(close(stimulated_rate(5400.0, 1800.0, 5e-22, 1.0), 1.00326, 1e-5)); // paper ≈ 1 Hz
        assert!(close(spontaneous_rate(8570.0, 5e3, 0.5, 1.0), 9.27562e-34, 1e-5)); // paper ≈ 10⁻³³
        let be = optimal_bar_mass(13e3, 2e-22, 1.19 * SOLAR_MASS, 100.0);
        let al = optimal_bar_mass(5.4e3, 2e-22, 1.19 * SOLAR_MASS, 150.0);
        assert!(close(be, 12.8458, 1e-5) && close(al, 219.501, 1e-5));
        assert!(close(be, 15.0, 0.25) && close(al, 250.0, 0.25)); // the paper's rounded values
    }

    #[test]
    fn gw170817_in_a_bar() {
        let sigma = bar_cross_section(2300.0, 5400.0, 3.0, 0.6, 0.345, 1.0, FRAC_PI_2);
        assert!(close(sigma, 4.31988e-25, 1e-5)); // Aguiar: 4×10⁻²⁵ m² Hz
        let fluence = energy_spectrum(1.188 * SOLAR_MASS, 900.0) / (4.0 * PI * (40.0 * MPC).powi(2));
        let e = sigma * fluence;
        assert!(close(e, 1.15207e-31, 1e-5));
        assert!(close(e / (HBAR * 2.0 * PI * 900.0), 0.193189, 1e-5));
        assert!(close(e / (C * 2300.0), 1.67083e-43, 1e-5));
    }
}
