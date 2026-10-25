//! The black hole a merger leaves: its spin, the energy radiated, and its horizon area, for
//! binaries with spins aligned or anti-aligned with the orbit.
//!
//! The spin is the fit of [Hofmann, Barausse & Rezzolla 2016](https://arxiv.org/abs/1605.01938)
//! (HBR16, Eq. 7) with the six coefficients `n_M = 1, n_J = 2` that the authors find
//! avoids overfitting (their χ²_red is 0.985, against 0.69 for 20 coefficients). It is exact in
//! the test-particle limit and gives the equal-mass, non-spinning `0.68646` of Scheel et al. by
//! construction. The radiated energy is [Barausse, Morozova & Rezzolla 2012](https://arxiv.org/abs/1206.3803)
//! (BMR12, their final expression), exact in both limits too.
//!
//! ```
//! use gravitational::remnant::{final_spin, horizon_area, radiated_energy};
//! use metrology::SOLAR_MASS;
//!
//! // GW150914: 36 + 29 M☉ with spins near zero left 62 ± 4 M☉ spinning at 0.67 (+0.05/−0.07).
//! let (m1, m2) = (36.0, 29.0);
//! let spin = final_spin(m1, m2, 0.0, 0.0);
//! let mass = (m1 + m2) * (1.0 - radiated_energy(m1, m2, 0.0, 0.0));
//! assert!((0.60..0.72).contains(&spin) && (mass - 62.0).abs() < 4.0);
//! // Hawking's area theorem: the remnant's horizon is larger than the two it swallowed.
//! let area = |m: f64, chi: f64| horizon_area(m * SOLAR_MASS, chi);
//! assert!(area(mass, spin) > area(m1, 0.0) + area(m2, 0.0));
//! ```

use core::f64::consts::PI;

use metrology::{C, G};

use crate::inspiral::kerr_isco;

/// HBR16's `k_ij` for `n_M = 1, n_J = 2` (their Table I, top block) and `ξ`. `k₀₀` follows from
/// the equal-mass constraint `√3/2 + k₀₀/16 + k₁₀/64 = 0.68646`.
const K01: f64 = -1.2019;
const K02: f64 = -1.20764;
const K10: f64 = 3.79245;
const K11: f64 = 1.18385;
const K12: f64 = 4.90494;
const XI: f64 = 0.41616;
/// Scheel et al. (2009)'s equal-mass, non-spinning final spin, which HBR16 impose.
const EQUAL_MASS_SPIN: f64 = 0.68646;

/// BMR12's `p₀` and `p₁`, fitted to equal-mass binaries with aligned spins.
const P0: f64 = 0.04827;
const P1: f64 = 0.01707;

/// The symmetric mass ratio `ν` and the total spin `a_tot = (a₁m₁² + a₂m₂²)/(m₁ + m₂)²`.
fn ratios(m1: f64, m2: f64, a1: f64, a2: f64) -> (f64, f64) {
    let m = m1 + m2;
    (m1 * m2 / (m * m), (a1 * m1 * m1 + a2 * m2 * m2) / (m * m))
}

/// The spin `cJ/GM²` of the black hole left by two of masses `m1`, `m2` (any unit) with spins
/// `a1`, `a2` projected on the orbital angular momentum (−1 to 1) (HBR16 Eq. 7):
/// `a_tot + ν[L_ISCO(a_eff) − 2a_tot(E_ISCO(a_eff) − 1)] + Σ k_ij ν^(2+i) a_eff^j`, with
/// `a_eff = a_tot + ξν(a₁ + a₂)`.
pub fn final_spin(m1: f64, m2: f64, a1: f64, a2: f64) -> f64 {
    let (nu, a_tot) = ratios(m1, m2, a1, a2);
    let a_eff = a_tot + XI * nu * (a1 + a2);
    let (_, e, l) = kerr_isco(a_eff);
    let k00 = 16.0 * (EQUAL_MASS_SPIN - 3f64.sqrt() / 2.0) - K10 / 4.0;
    let fit = nu * nu * (k00 + K01 * a_eff + K02 * a_eff * a_eff) + nu.powi(3) * (K10 + K11 * a_eff + K12 * a_eff * a_eff);
    a_tot + nu * (l - 2.0 * a_tot * (e - 1.0)) + fit
}

/// The fraction of the binary's mass radiated in gravitational waves, for the same arguments as
/// [`final_spin`] (BMR12): `[1 − E_ISCO(ã)]ν + 4ν²[4p₀ + 16p₁ã(ã + 1) + E_ISCO(ã) − 1]`, with
/// `ã = a_tot`. 4.83% for equal masses without spin, 9.95% with both spins maximal and aligned.
pub fn radiated_energy(m1: f64, m2: f64, a1: f64, a2: f64) -> f64 {
    let (nu, a) = ratios(m1, m2, a1, a2);
    let (_, e, _) = kerr_isco(a);
    (1.0 - e) * nu + 4.0 * nu * nu * (4.0 * P0 + 16.0 * P1 * a * (a + 1.0) + e - 1.0)
}

/// The horizon area (m²) of a Kerr black hole of mass `mass` (kg) and spin `chi`:
/// `8π(GM/c²)²(1 + √(1 − χ²))`, from `16π(GM/c²)²` without spin down to half that at `χ = 1`.
pub fn horizon_area(mass: f64, chi: f64) -> f64 {
    let rg = G * mass / (C * C);
    8.0 * PI * rg * rg * (1.0 + (1.0 - chi * chi).sqrt())
}

#[cfg(test)]
mod tests {
    use super::*;
    use metrology::SOLAR_MASS;

    #[test]
    fn the_final_spin() {
        // Equal masses without spin: Scheel et al.'s 0.68646, through the derived k₀₀, which HBR16
        // quote as "≈ −3.82".
        assert!((final_spin(1.0, 1.0, 0.0, 0.0) - 0.68646).abs() < 1e-12);
        let k00 = 16.0 * (EQUAL_MASS_SPIN - 3f64.sqrt() / 2.0) - K10 / 4.0;
        assert!((k00 + 3.82).abs() < 0.005, "{k00}");
        // The test-particle limit (HBR16 Eq. 1): a_fin − a₁ → ν(L_ISCO(a₁) − 2a₁E_ISCO(a₁)).
        for a1 in [-0.9, 0.0, 0.5, 0.95] {
            let q = 1e-6;
            let nu = q / (1.0 + q) / (1.0 + q);
            let (_, e, l) = kerr_isco(a1);
            let slope = (final_spin(1.0, q, a1, 0.0) - a1) / nu;
            assert!((slope - (l - 2.0 * a1 * e)).abs() < 1e-4, "a₁ = {a1}: {slope}");
        }
        // Either body can be first, and aligned spins spin the remnant up.
        assert!((final_spin(29.0, 36.0, 0.3, -0.2) - final_spin(36.0, 29.0, -0.2, 0.3)).abs() < 1e-15);
        assert!(final_spin(1.0, 1.0, 0.9, 0.9) > final_spin(1.0, 1.0, 0.0, 0.0));
        assert!(final_spin(1.0, 1.0, 1.0, 1.0) < 1.0);
    }

    #[test]
    fn the_radiated_energy() {
        // Equal masses: p₀ = 4.827% without spin and p₀ + 3p₁ = 9.95% at a₁ = a₂ = 1 (BMR12).
        assert!((radiated_energy(1.0, 1.0, 0.0, 0.0) - 0.04827).abs() < 1e-15);
        assert!((radiated_energy(1.0, 1.0, 1.0, 1.0) - 0.0995).abs() < 5e-5);
        // The test-particle limit: 1 − E_ISCO per unit ν, 1 − √(8/9) without spin.
        let nu = 1e-6 / (1.0 + 1e-6f64).powi(2);
        assert!((radiated_energy(1.0, 1e-6, 0.0, 0.0) / nu - (1.0 - (8.0f64 / 9.0).sqrt())).abs() < 1e-5);
    }

    #[test]
    fn the_horizon_and_the_area_theorem() {
        let rg = G * SOLAR_MASS / (C * C);
        assert!((horizon_area(SOLAR_MASS, 0.0) / (16.0 * PI * rg * rg) - 1.0).abs() < 1e-15);
        assert!((horizon_area(SOLAR_MASS, 1.0) / (8.0 * PI * rg * rg) - 1.0).abs() < 1e-15);
        // An equal-mass, non-spinning merger grows the horizon area by 56%.
        let (spin, mass) = (final_spin(1.0, 1.0, 0.0, 0.0), 2.0 * (1.0 - radiated_energy(1.0, 1.0, 0.0, 0.0)));
        let growth = horizon_area(mass * SOLAR_MASS, spin) / (2.0 * horizon_area(SOLAR_MASS, 0.0));
        assert!((growth - 1.5644).abs() < 1e-3, "{growth}");
    }
}
