//! Electricity without a steam cycle: direct conversion of charged fusion products, inductive
//! recovery in pulsed machines, and the plant arithmetic that decides whether any of it nets power.
//!
//! A direct converter is not a heat engine, so Carnot does not limit it; the spread of the ions'
//! energies against the number of collector stages does. Measured efficiencies are 48–65%; the
//! 70–90% figures for multi-stage, travelling-wave and inverse-cyclotron converters are projections.
//!
//! ```
//! use plasmoid::conversion::collector_efficiency;
//!
//! // An ideal one-stage collector on a thermal ion flux returns 42%; eight stages 82%.
//! assert!((collector_efficiency(1, 0.0) - 0.419).abs() < 0.003);
//! assert!((collector_efficiency(8, 0.0) - 0.824).abs() < 0.003);
//! ```

use core::f64::consts::PI;

use crate::aneutronic::Q_PB11;
use metrology::KEV;

/// Best efficiency of an ideal `stages`-stage electrostatic collector on an ion flux
/// `f(E) ∝ (E − s)e^(−(E − s))` (energies in units of the temperature, `shift` = s the potential
/// drop the ions fall through first): each ion is caught by the highest stage it can reach and
/// returns that stage's potential. The stage potentials are optimised by coordinate search.
pub fn collector_efficiency(stages: usize, shift: f64) -> f64 {
    let (n, e_max) = (800usize, 40.0);
    let h = (e_max - shift) / n as f64;
    let flux: Vec<(f64, f64)> = (0..n).map(|i| {
        let e = shift + (i as f64 + 0.5) * h;
        (e, (e - shift) * (-(e - shift)).exp())
    }).collect();
    let total: f64 = flux.iter().map(|(e, f)| e * f).sum();
    let eff = |v: &[f64]| flux.iter().map(|(e, f)| f * v.iter().copied().filter(|x| x <= e).fold(0.0, f64::max)).sum::<f64>() / total;
    let mut v: Vec<f64> = (0..stages).map(|k| shift + (k + 1) as f64).collect();
    let mut step = 0.5;
    while step > 0.01 {
        let mut improved = false;
        for k in 0..stages {
            for d in [-step, step] {
                let mut w = v.clone();
                w[k] = (w[k] + d).max(0.0);
                w.sort_by(f64::total_cmp);
                if eff(&w) > eff(&v) + 1e-9 {
                    v = w;
                    improved = true;
                }
            }
        }
        if !improved {
            step /= 2.0;
        }
    }
    eff(&v)
}

/// Engineering gain of a steady magnetic-confinement plant: electricity out over electricity in,
/// `η_c η_d (Q_sci + 1)`, for conversion efficiency `eta_c` (heat or charged power to
/// electricity) and heating-system efficiency `eta_d`.
pub fn q_eng_magnetic(q_sci: f64, eta_c: f64, eta_d: f64) -> f64 {
    eta_c * eta_d * (q_sci + 1.0)
}

/// Engineering gain of an inertial plant: `η_c η_d G` for target gain `gain` and driver
/// efficiency `eta_d`.
pub fn q_eng_inertial(gain: f64, eta_c: f64, eta_d: f64) -> f64 {
    eta_c * eta_d * gain
}

/// The scientific gain a pulsed machine with inductive energy recovery needs to break even:
/// `(1 − η_rec)/(η_dc f_ch)` for recovery efficiency `eta_rec`, direct-conversion efficiency
/// `eta_dc` and charged-particle fraction `f_ch`.
pub fn pulsed_breakeven(eta_rec: f64, eta_dc: f64, f_ch: f64) -> f64 {
    (1.0 - eta_rec) / (eta_dc * f_ch)
}

/// Fusion energy out over laser energy in for a p–¹¹B shot that reports `alphas_per_sr` alphas
/// per steradian from `laser_energy` (J), assuming isotropic emission (three alphas per fusion).
pub fn laser_pb11_gain(alphas_per_sr: f64, laser_energy: f64) -> f64 {
    4.0 * PI * alphas_per_sr / 3.0 * Q_PB11 * 1e3 * KEV / laser_energy
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, rel: f64) -> bool {
        (a / b - 1.0).abs() < rel
    }

    #[test]
    fn collectors() {
        for (stages, shift, expect) in [(1, 0.0, 0.419), (2, 0.0, 0.593), (8, 0.0, 0.824), (8, 5.0, 0.951)] {
            let e = collector_efficiency(stages, shift);
            assert!((e - expect).abs() < 0.003, "{stages} {shift}: {e}");
        }
    }

    #[test]
    fn plants() {
        // A steam cycle (40%) and 70% heating need Q_sci = 2.57 to break even; 50% heating, 4.
        assert!(close(1.0 / (0.4 * 0.7) - 1.0, 2.571, 1e-3));
        assert!((q_eng_magnetic(1.0 / (0.4 * 0.5) - 1.0, 0.4, 0.5) - 1.0).abs() < 1e-12);
        // NIF 2025 (8.6 MJ from 2.08 MJ) against ~400 MJ at the wall: Q_eng ≈ 0.009 at 40%.
        assert!(close(0.4 * 8.6 / 400.0, 0.0086, 1e-9));
        assert!(close(q_eng_inertial(140.0, 0.4, 0.07), 3.92, 1e-3));
        // Pulsed inductive recovery: Q_sci 0.058 at 95% recovery, 0.47 at 60%.
        assert!(close(pulsed_breakeven(0.95, 0.9, 0.95), 0.05848, 1e-3));
        assert!(close(pulsed_breakeven(0.6, 0.9, 0.95), 0.4678, 1e-3));
    }

    #[test]
    fn laser_proton_boron() {
        // Margarone 2020: 5×10⁹ α/sr from 1.4 kJ gives 2.1×10⁻⁵ of the laser energy back.
        assert!(close(laser_pb11_gain(5e9, 1.4e3), 2.081e-5, 2e-3));
        assert!(close(laser_pb11_gain(1.7e9, 550.0), 1.801e-5, 2e-3));
    }
}
