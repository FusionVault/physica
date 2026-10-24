//! Cold-atom (light-pulse) interferometers as space accelerometers and gradiometers, at the
//! quantum-projection-noise limit.
//!
//! In microgravity the Mach–Zehnder π/2–π–π/2 sequence uses double diffraction, so the phase is
//! `Δφ = 2·k_eff·a·T²` (Lévèque et al. 2021, arXiv:2011.03382; Trimeche et al. 2019,
//! arXiv:1903.09828). Per-shot phase noise is `1/(C·√N)`, and the white noise density is the
//! per-shot noise times `√T_c`.

/// Two-photon Raman wave number `k_eff = 4π/λ` for ⁸⁷Rb on the D2 line (780.241 nm), m⁻¹.
pub const K_EFF_RB87: f64 = 4.0 * core::f64::consts::PI / 780.241e-9;

/// Interferometer phase per unit acceleration, rad per m/s²: `k_eff T²` (Kasevich–Chu), for
/// effective wavenumber `k_eff` (m⁻¹) and pulse separation `t` (s). Rosi et al.'s 0.16 s with the
/// doubled ⁸⁷Rb wavenumber gives 4.12×10⁵ rad per m/s², so Earth's gravity is 4.0×10⁶ rad. The
/// double-diffraction sequence of [`AtomInterferometer::phase`] carries a further factor of 2.
pub fn phase_sensitivity(k_eff: f64, t: f64) -> f64 {
    k_eff * t * t
}

/// Rate of gravitational Aharonov–Bohm phase accumulation (rad/s) for an atom of mass `m` (kg) at
/// distance `r` (m) from a source mass `source` (kg): `mGM/(ħr)`. The phase follows the
/// potential, not the force, which is what
/// [Overstreet et al. 2022](https://doi.org/10.1126/science.abl7152) measured: 1.63 rad/s for
/// ⁸⁷Rb at 7 cm from 1.25 kg.
pub fn aharonov_bohm_rate(m: f64, source: f64, r: f64) -> f64 {
    m * metrology::G * source / (metrology::HBAR * r)
}

/// One atom interferometer (or each of a gradiometer's two).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AtomInterferometer {
    /// Time between pulses `T`, s (the interferometer lasts `2T`).
    pub t: f64,
    /// Fringe contrast `C`, 0–1.
    pub contrast: f64,
    /// Atoms detected per shot `N`.
    pub atoms: f64,
    /// Cycle time `T_c`, s: `2T + T_prep` when shots run one after another, `T_prep` when
    /// interleaved.
    pub cycle: f64,
}

impl AtomInterferometer {
    /// Interferometer phase for acceleration `a` (m/s²), radians: `2·k_eff·a·T²`.
    pub fn phase(&self, a: f64) -> f64 {
        2.0 * K_EFF_RB87 * a * self.t * self.t
    }

    /// Quantum-projection phase noise per shot, `1/(C·√N)`, radians.
    pub fn phase_noise(&self) -> f64 {
        1.0 / (self.contrast * self.atoms.sqrt())
    }

    /// Acceleration noise per shot, m/s².
    pub fn acceleration_noise(&self) -> f64 {
        self.phase_noise() / (2.0 * K_EFF_RB87 * self.t * self.t)
    }

    /// White acceleration noise density, m s⁻² Hz⁻¹ᐟ²; use it as
    /// [`Accelerometer::AtomInterferometer`](crate::missions::Accelerometer::AtomInterferometer).
    pub fn acceleration_asd(&self) -> f64 {
        self.acceleration_noise() * self.cycle.sqrt()
    }

    /// Gravity-gradient noise per shot of a gradiometer made of two such interferometers a
    /// `baseline` (m) apart sharing the same lasers, s⁻²: `√2·σ_φ/(2·k_eff·D·T²)`.
    pub fn gradient_noise(&self, baseline: f64) -> f64 {
        std::f64::consts::SQRT_2 * self.phase_noise() / (2.0 * K_EFF_RB87 * baseline * self.t * self.t)
    }

    /// White gravity-gradient noise density of that gradiometer, s⁻² Hz⁻¹ᐟ².
    pub fn gradient_asd(&self, baseline: f64) -> f64 {
        self.gradient_noise(baseline) * self.cycle.sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, rel: f64) -> bool {
        (a - b).abs() <= rel * b.abs()
    }

    /// Published instrument parameter sets, recomputed from first principles.
    #[test]
    fn published_instruments() {
        // CARIOQA (Lévèque et al. 2022): T = 1 s, C = 0.8, N = 1e5, T_c = 3.5 s. The paper quotes
        // 3e-10 (conservative) and 8e-11 with N = 1e6; QPN alone gives 2.3e-10 and 7.3e-11.
        let carioqa = AtomInterferometer { t: 1.0, contrast: 0.8, atoms: 1e5, cycle: 3.5 };
        assert!(close(carioqa.phase_noise(), 3.95e-3, 0.01));
        assert!(close(carioqa.acceleration_asd(), 2.3e-10, 0.01));
        assert!(close(AtomInterferometer { atoms: 1e6, ..carioqa }.acceleration_asd(), 7.3e-11, 0.01));
        // GRICE (Lévèque et al. 2021): T = 0.5 s, C = 0.43, N = 5e5, T_c = 2 s → 6e-10.
        let grice = AtomInterferometer { t: 0.5, contrast: 0.43, atoms: 5e5, cycle: 2.0 };
        assert!(close(grice.acceleration_asd(), 5.8e-10, 0.01));
        // Trimeche et al. 2019 gradiometer: T = 5 s, N = 1e6, ideal contrast, D = 0.5 m → 3.5 mE
        // per shot, and a 1087 rad differential phase for T_zz = 2.7e-6 s⁻².
        let tri = AtomInterferometer { t: 5.0, contrast: 1.0, atoms: 1e6, cycle: 1.0 };
        assert!(close(tri.gradient_noise(0.5), 3.5e-12, 0.01));
        assert!(close(tri.phase(2.7e-6 * 0.5), 1087.0, 0.001));
    }

    #[test]
    fn phase_and_potential() {
        // Rosi et al. 2014: k_eff T^2 = 4.123e5 rad per m/s^2, so g0 is 4.04e6 rad.
        let s = phase_sensitivity(K_EFF_RB87, 0.160);
        assert!((s / 4.123_07e5 - 1.0).abs() < 1e-5 && (s * metrology::G0 / 4.043e6 - 1.0).abs() < 1e-3);
        // Overstreet et al. 2022: 1.63 rad/s for Rb-87 at 7 cm from 1.25 kg.
        let rb = 86.909 * metrology::AMU;
        assert!((aharonov_bohm_rate(rb, 1.25, 0.07) / 1.631 - 1.0).abs() < 1e-3);
        // It falls as 1/r, unlike a force.
        assert!((aharonov_bohm_rate(rb, 1.25, 0.14) / aharonov_bohm_rate(rb, 1.25, 0.07) - 0.5).abs() < 1e-12);
    }
}
