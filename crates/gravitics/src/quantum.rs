//! Tabletop tests of quantum gravity: entanglement through gravity (Bose–Marletto–Vedral) and
//! the radiation predicted by gravitational collapse of the wavefunction (Diósi–Penrose).
//!
//! In the BMV proposal two 10⁻¹⁴ kg microdiamonds, each in a superposition of two places, become
//! entangled only if gravity can carry quantum information; the phases are `Gm²τ/(ħr)`
//! ([Bose et al. 2017](https://arxiv.org/abs/1707.06050)). The parameter-free Diósi–Penrose
//! model predicted spontaneous X-rays from germanium at about ten per day; an underground
//! detector saw none, which refutes that version ([Donadi et al. 2021](https://arxiv.org/abs/2111.13490)).
//!
//! ```
//! use gravitics::quantum::bmv_phases;
//!
//! // Bose et al.: 10⁻¹⁴ kg, 450 µm apart, split by 250 µm, held 2.5 s.
//! let (near, closer, farther) = bmv_phases(1e-14, 450e-6, 250e-6, 2.5);
//! assert!(((closer - near) - 0.440).abs() < 1e-3 && ((farther - near) + 0.126).abs() < 1e-3);
//! ```

use core::f64::consts::PI;

use metrology::{C, E_CHARGE, EPS_0, G, HBAR};

/// Gravitational phases `Gm²τ/(ħr)` (radians) between two masses `m` (kg) held for `tau` (s) at
/// separations `d`, `d − dx` and `d + dx` (m), the three distances of the BMV interferometer.
pub fn bmv_phases(m: f64, d: f64, dx: f64, tau: f64) -> (f64, f64, f64) {
    let phase = |r: f64| G * m * m * tau / (HBAR * r);
    (phase(d), phase(d - dx), phase(d + dx))
}

/// Spontaneous photons per second between energies `e_lo` and `e_hi` (any common unit) that the
/// parameter-free Diósi–Penrose model predicts from `n_atoms` atoms of nuclear charge `z`, for
/// mass-density smearing radius `r0` (m), at unit detection efficiency:
/// `2Ge²Z²N/(3π^(3/2)ε₀c³R₀³) ln(e_hi/e_lo)` (Donadi et al. 2021).
pub fn dp_photon_rate(r0: f64, z: f64, n_atoms: f64, e_lo: f64, e_hi: f64) -> f64 {
    2.0 * G * E_CHARGE * E_CHARGE * z * z * n_atoms / (3.0 * PI.powf(1.5) * EPS_0 * C.powi(3) * r0.powi(3)) * (e_hi / e_lo).ln()
}

#[cfg(test)]
mod tests {
    use super::*;
    use metrology::N_A;

    #[test]
    fn bose_marletto_vedral() {
        let (p0, p_rl, p_lr) = bmv_phases(1e-14, 450e-6, 250e-6, 2.5);
        assert!((p0 / 0.351607 - 1.0).abs() < 1e-5 && (p_rl / 0.791115 - 1.0).abs() < 1e-5 && (p_lr / 0.226033 - 1.0).abs() < 1e-5);
    }

    #[test]
    fn diosi_penrose_germanium() {
        // 375 cm³ of germanium (5.323 g/cm³, 72.630 g/mol), 1–3.8 MeV, per day.
        let atoms = 375.0 * 5.323 / 72.630 * N_A;
        assert!((atoms / 1.655e25 - 1.0).abs() < 1e-3);
        let day = 86_400.0;
        assert!((dp_photon_rate(0.54e-10, 32.0, atoms, 1.0, 3.8) * day / 10.6742 - 1.0).abs() < 1e-5);
        assert!((dp_photon_rate(0.05e-10, 32.0, atoms, 1.0, 3.8) * day / 1.34464e4 - 1.0).abs() < 1e-5);
    }
}
