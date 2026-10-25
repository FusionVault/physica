//! Distances in a flat ΛCDM universe: what a lensed quasar's time delays, a standard siren or a
//! supernova measure.
//!
//! `D_M(z) = (c/H₀)∫₀^z dz′/E(z′)` with `E(z) = √(Ω_m(1 + z)³ + 1 − Ω_m)`, radiation neglected
//! (its `Ω_r ≈ 9×10⁻⁵` moves `D_M` by about 10⁻⁴ at `z = 2`). In a flat universe the angular
//! diameter distance between redshifts `z₁ < z₂` is `(D_M(z₂) − D_M(z₁))/(1 + z₂)`, and a lens at
//! `z_d` in front of a source at `z_s` has the time-delay distance `(1 + z_d)D_d D_s/D_ds`, which
//! scales as `1/H₀` ([Wong et al. 2020](https://arxiv.org/abs/1907.04869)).
//!
//! ```
//! use gravitational::H100;
//! use gravitational::cosmology::{angular_diameter_distance, time_delay_distance};
//! use metrology::MPC;
//!
//! // RXJ1131−1231 (z_d = 0.295, z_s = 0.654) at H₀ = 70 km/s/Mpc and Ω_m = 0.3: 2334 Mpc,
//! // against the measured 2096 (+98/−83), so this lens alone says H₀ = 70 × 2334/2096 = 78.
//! let ddt = time_delay_distance(0.295, 0.654, 0.7 * H100, 0.3) / MPC;
//! assert!((ddt - 2334.06).abs() < 0.01);
//! assert!(angular_diameter_distance(0.0, 0.295, 0.7 * H100, 0.3) / MPC > 800.0);
//! ```

use metrology::C;

use crate::integrate_log;

/// The transverse comoving distance (m) to redshift `z` in a flat ΛCDM universe with Hubble
/// constant `h0` (s⁻¹, e.g. `0.7 * H100`) and matter density `omega_m`.
pub fn comoving_distance(z: f64, h0: f64, omega_m: f64) -> f64 {
    // Over x = 1 + z, so the log-grid Simpson rule starts above zero.
    C / h0 * integrate_log(1.0, 1.0 + z, 2000, |x| 1.0 / (omega_m * x * x * x + 1.0 - omega_m).sqrt())
}

/// The angular diameter distance (m) of an object at redshift `z2` seen from redshift `z1 < z2`
/// (`z1 = 0` for the observer), in the same universe.
pub fn angular_diameter_distance(z1: f64, z2: f64, h0: f64, omega_m: f64) -> f64 {
    (comoving_distance(z2, h0, omega_m) - comoving_distance(z1, h0, omega_m)) / (1.0 + z2)
}

/// The time-delay distance (m) of a lens at redshift `zd` and a source at `zs`:
/// `(1 + z_d)D_d D_s/D_ds`.
pub fn time_delay_distance(zd: f64, zs: f64, h0: f64, omega_m: f64) -> f64 {
    let (dd, ds, dds) = (angular_diameter_distance(0.0, zd, h0, omega_m), angular_diameter_distance(0.0, zs, h0, omega_m), angular_diameter_distance(zd, zs, h0, omega_m));
    (1.0 + zd) * dd * ds / dds
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::H100;
    use metrology::MPC;

    #[test]
    fn the_two_closed_forms() {
        // Einstein–de Sitter (Ω_m = 1): D_M = (2c/H₀)(1 − 1/√(1 + z)); de Sitter (Ω_m = 0): cz/H₀.
        let h0 = 0.7 * H100;
        for z in [0.1, 1.0, 3.0] {
            let eds = 2.0 * C / h0 * (1.0 - 1.0 / (1.0f64 + z).sqrt());
            assert!((comoving_distance(z, h0, 1.0) / eds - 1.0).abs() < 1e-13, "z = {z}");
            assert!((comoving_distance(z, h0, 0.0) / (C * z / h0) - 1.0).abs() < 1e-14, "z = {z}");
        }
        assert_eq!(comoving_distance(0.0, h0, 0.3), 0.0);
    }

    #[test]
    fn h0licow() {
        // Wong et al. (2020), Tables 1–2: the six lenses' (z_d, z_s) and median D_Δt (Mpc) with
        // its 16th–84th percentiles. D_Δt ∝ 1/H₀, so each gives H₀ = 70 D_Δt(70)/D_Δt; their
        // inverse-variance mean at Ω_m = 0.3 is the paper's flat-ΛCDM 73.3 +1.7/−1.8 to 0.1.
        let lenses = [
            (0.6304, 1.394, 5156.0, 296.0, 236.0),
            (0.295, 0.654, 2096.0, 98.0, 83.0),
            (0.4546, 1.693, 2707.0, 183.0, 168.0),
            (0.745, 1.789, 5769.0, 589.0, 471.0),
            (0.6575, 1.662, 4784.0, 399.0, 248.0),
            (0.311, 1.722, 1470.0, 137.0, 127.0),
        ];
        let (mut sum, mut weights) = (0.0, 0.0);
        for (zd, zs, ddt, up, down) in lenses {
            let h = 70.0 * time_delay_distance(zd, zs, 0.7 * H100, 0.3) / MPC / ddt;
            let sigma = h * 0.5 * (up + down) / ddt;
            sum += h / (sigma * sigma);
            weights += 1.0 / (sigma * sigma);
        }
        let (h0, sigma) = (sum / weights, weights.powf(-0.5));
        assert!((h0 - 73.3).abs() < 0.1 && (sigma - 1.75).abs() < 0.1, "{h0} ± {sigma}");
    }
}
