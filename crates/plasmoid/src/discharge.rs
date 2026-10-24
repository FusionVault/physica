//! Gas discharges: when a gas gap breaks down, how far a molecule flies between collisions, and
//! when a wire goes into corona. Together they decide whether a high-voltage rig "in vacuum" is
//! free of ions, or full of them.
//!
//! Townsend breakdown follows Paschen's law, a function of pressure × gap alone:
//! `V_b = B·pd / (ln(A·pd) − ln ln(1 + 1/γ))`, where `α/p = A exp(−Bp/E)` is the gas's first
//! ionisation coefficient and γ the electrode's secondary-electron yield. The curve has a
//! minimum, about 300 V for air near 1 Pa·m, and no breakdown to its left, where an electron
//! crosses the gap without ionising enough. A rough vacuum of 1–100 Pa in a 10 cm vessel sits
//! right at the minimum, so a few hundred volts light it.
//!
//! ```
//! use plasmoid::discharge::{AIR, AIR_DIAMETER, mean_free_path, paschen};
//!
//! // A 10 cm tank at 10 Pa breaks down at about 310 V; at 1 Pa the gap is too thin to.
//! assert!((paschen(AIR, 0.01, 10.0 * 0.1).unwrap() - 307.0).abs() < 1.0);
//! assert!(paschen(AIR, 0.01, 1.0 * 0.1).is_none());
//!
//! // Below about 0.07 Pa a molecule crosses the tank without meeting another.
//! assert!(mean_free_path(0.066, 293.15, AIR_DIAMETER) > 0.1);
//! ```

use core::f64::consts::{E, PI};

use metrology::K_B;

/// A gas's Townsend constants: `a` (1/(Pa m)) and `b` (V/(Pa m)) in `α/p = a exp(−bp/E)`.
#[derive(Clone, Copy, Debug)]
pub struct Townsend {
    pub a: f64,
    pub b: f64,
}

const TORR_CM: f64 = 133.322_368 * 0.01;

/// Air: A = 15 /(cm Torr), B = 365 V/(cm Torr), fitted for E/p from 100 to 800 V/(cm Torr)
/// (Raizer, *Gas Discharge Physics*, 1991, table 4.1). For a 10 cm gap that covers pressures up
/// to about 150 Pa.
pub const AIR: Townsend = Townsend { a: 15.0 / TORR_CM, b: 365.0 / TORR_CM };

/// Air's Lennard-Jones collision diameter (m), 3.711 Å (Hirschfelder, Curtiss & Bird, *Molecular
/// Theory of Gases and Liquids*, 1954).
pub const AIR_DIAMETER: f64 = 3.711e-10;

/// The Townsend breakdown voltage (V) of a gap at pressure × distance `pd` (Pa m), for a cathode
/// of secondary yield `gamma` (about 0.01 for metals in air). `None` to the left of the curve,
/// where the gap cannot break down at any voltage.
pub fn paschen(gas: Townsend, gamma: f64, pd: f64) -> Option<f64> {
    let denominator = (gas.a * pd).ln() - (1.0 + 1.0 / gamma).ln().ln();
    (denominator > 0.0).then(|| gas.b * pd / denominator)
}

/// The minimum of the Paschen curve: `(pd, V)` in (Pa m, V), at `pd = e ln(1 + 1/γ)/A`.
pub fn paschen_minimum(gas: Townsend, gamma: f64) -> (f64, f64) {
    let pd = E * (1.0 + 1.0 / gamma).ln() / gas.a;
    (pd, gas.b * pd)
}

/// The mean free path (m) of a gas molecule of collision diameter `diameter` (m) at `pressure`
/// (Pa) and `temperature` (K): `k_B T/(√2 π d² p)`.
pub fn mean_free_path(pressure: f64, temperature: f64, diameter: f64) -> f64 {
    K_B * temperature / (2f64.sqrt() * PI * diameter * diameter * pressure)
}

/// The surface field (V/m) at which a smooth wire of radius `radius` (m) goes into corona in air
/// of relative density `relative_density` (1 at 20 °C and 1 atm): Peek's law
/// `E_c = 3.0 MV/m · δ(1 + 0.0301/√(δr))`, with `r` in metres (Peek, *Dielectric Phenomena in
/// High-Voltage Engineering*, 1929). Thin wires need more field, but reach it at lower voltage.
pub fn peek_onset(radius: f64, relative_density: f64) -> f64 {
    3.0e6 * relative_density * (1.0 + 0.0301 / (relative_density * radius).sqrt())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, rel: f64) -> bool {
        (a / b - 1.0).abs() < rel
    }

    #[test]
    fn the_paschen_minimum_of_air() {
        // 305 V at 1.115 Pa m (0.836 Torr cm), close to the measured minimum of air near 330 V.
        let (pd, v) = paschen_minimum(AIR, 0.01);
        assert!(close(pd, 1.1150, 1e-4) && close(v, 305.27, 1e-4), "{pd} {v}");
        assert!(close(pd / TORR_CM, 0.8363, 1e-3));
        // It is the minimum: both sides are higher.
        let at = |x: f64| paschen(AIR, 0.01, x).unwrap();
        assert!(close(at(pd), v, 1e-12) && at(0.9 * pd) > v && at(1.1 * pd) > v);
        // Left of where the denominator vanishes, no breakdown at any voltage.
        let edge = (1.0 + 1.0 / 0.01f64).ln() / AIR.a;
        assert!(paschen(AIR, 0.01, 0.99 * edge).is_none() && paschen(AIR, 0.01, 1.01 * edge).is_some());
    }

    #[test]
    fn a_ten_centimetre_tank() {
        // The patent's tank: 857 V at 100 Pa, 307 V at 10 Pa, nothing at 1 Pa.
        assert!(close(paschen(AIR, 0.01, 100.0 * 0.1).unwrap(), 857.23, 1e-4));
        assert!(close(paschen(AIR, 0.01, 10.0 * 0.1).unwrap(), 307.24, 1e-4));
        assert!(paschen(AIR, 0.01, 1.0 * 0.1).is_none());
    }

    #[test]
    fn mean_free_path_of_air() {
        // The vacuum-handbook rule for air at 20 °C: λ·p ≈ 6.6×10⁻³ m Pa.
        let lambda = mean_free_path(1.0, 293.15, AIR_DIAMETER);
        assert!(close(lambda, 6.6e-3, 0.02), "{lambda}");
        // At atmospheric pressure, 66 nm.
        assert!(close(mean_free_path(101_325.0, 293.15, AIR_DIAMETER) / lambda, 1.0 / 101_325.0, 1e-12));
    }

    #[test]
    fn peek_corona_onset() {
        // A 1 cm wire: 39.0 kV/cm. The patent's 0.75 mm copper: 63 kV/cm.
        assert!(close(peek_onset(0.01, 1.0), 3.903e6, 1e-3));
        assert!(close(peek_onset(0.75e-3, 1.0), 6.297e6, 1e-3));
        // Thinner air goes into corona sooner.
        assert!(peek_onset(0.01, 0.5) < peek_onset(0.01, 1.0));
    }
}
