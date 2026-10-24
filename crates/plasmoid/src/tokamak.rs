//! The engineering scalings a tokamak design is checked against: the edge safety factor (and the
//! Kruskal–Shafranov kink limit behind it), the Greenwald density limit, the Troyon beta limit,
//! and the IPB98(y,2) energy-confinement law.
//!
//! These are empirical fits to many machines, not first-principles MHD (that is the `magnetics`
//! crate). They are what reactor designs such as ITER and SPARC quote, so a design point can be
//! checked in a few lines.
//!
//! ```
//! use plasmoid::tokamak::ITER;
//!
//! // ITER's edge safety factor: 3.0, above the kink limit of 2 with margin.
//! assert!((ITER.q95() - 3.0).abs() < 0.01);
//! // At 85% of the Greenwald limit, 1.01×10²⁰ m⁻³.
//! assert!((0.85 * ITER.greenwald() / 1.01e20 - 1.0).abs() < 0.01);
//! ```

use core::f64::consts::PI;

/// A tokamak's shape, field and current.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Tokamak {
    /// Minor radius a, m.
    pub minor_radius: f64,
    /// Major radius R, m.
    pub major_radius: f64,
    /// Toroidal field on axis, T.
    pub field: f64,
    /// Plasma current, MA.
    pub current: f64,
    /// Elongation κ.
    pub elongation: f64,
    /// Triangularity δ.
    pub triangularity: f64,
}

/// ITER's design point (ITER Physics Basis): a = 2.0 m, R = 6.2 m, 5.3 T, 15 MA, κ = 1.7, δ = 0.33.
pub const ITER: Tokamak = Tokamak { minor_radius: 2.0, major_radius: 6.2, field: 5.3, current: 15.0, elongation: 1.7, triangularity: 0.33 };
/// SPARC's design point (Creely et al. 2020): a = 0.57 m, R = 1.85 m, 12.2 T, 8.7 MA, κ = 1.75, δ = 0.45.
pub const SPARC: Tokamak = Tokamak { minor_radius: 0.57, major_radius: 1.85, field: 12.2, current: 8.7, elongation: 1.75, triangularity: 0.45 };

impl Tokamak {
    /// Inverse aspect ratio `a/R`.
    pub fn inverse_aspect_ratio(&self) -> f64 {
        self.minor_radius / self.major_radius
    }

    /// The cylindrical safety factor `5a²B/(RI)` (I in MA): the field-line twist the
    /// Kruskal–Shafranov kink limit requires to exceed 1.
    pub fn q_cylindrical(&self) -> f64 {
        5.0 * self.minor_radius.powi(2) * self.field / (self.major_radius * self.current)
    }

    /// Uckan's shaped safety factor `q* = q_cyl [1 + κ²(1 + 2δ² − 1.2δ³)]/2`.
    pub fn q_star(&self) -> f64 {
        let (k, d) = (self.elongation, self.triangularity);
        self.q_cylindrical() * (1.0 + k * k * (1.0 + 2.0 * d * d - 1.2 * d.powi(3))) / 2.0
    }

    /// The ITER Physics Basis edge safety factor `q95 = q* (1.17 − 0.65ε)/(1 − ε²)²`.
    pub fn q95(&self) -> f64 {
        let e = self.inverse_aspect_ratio();
        self.q_star() * (1.17 - 0.65 * e) / (1.0 - e * e).powi(2)
    }

    /// The Greenwald density limit `I/(πa²)` × 10²⁰ m⁻³ (I in MA).
    pub fn greenwald(&self) -> f64 {
        self.current / (PI * self.minor_radius.powi(2)) * 1e20
    }

    /// The beta (as a fraction) allowed by normalised beta `beta_n` (% m T/MA): `β_N I/(aB)` %.
    /// Troyon's limit is β_N ≈ 2.8.
    pub fn troyon_beta(&self, beta_n: f64) -> f64 {
        beta_n * self.current / (self.minor_radius * self.field) / 100.0
    }

    /// The IPB98(y,2) energy-confinement time (s) at line-averaged density `n19` (10¹⁹ m⁻³), net
    /// heating power `power` (MW), ion mass `mass` (u), area elongation `kappa_a` and
    /// confinement-enhancement factor `h`:
    /// `0.0562 H I^0.93 B^0.15 n^0.41 P^−0.69 R^1.97 κ_a^0.78 ε^0.58 M^0.19`. Scatter ±14% (1σ).
    pub fn ipb98(&self, n19: f64, power: f64, mass: f64, kappa_a: f64, h: f64) -> f64 {
        0.0562 * h * self.current.powf(0.93) * self.field.powf(0.15) * n19.powf(0.41) * power.powf(-0.69) * self.major_radius.powf(1.97)
            * kappa_a.powf(0.78) * self.inverse_aspect_ratio().powf(0.58) * mass.powf(0.19)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, rel: f64) -> bool {
        (a / b - 1.0).abs() < rel
    }

    #[test]
    fn safety_factors() {
        // ITER 3.0; SPARC q* 3.05 and q95 ≈ 3.6 (Creely et al. 2020).
        assert!(close(ITER.q95(), 2.996_339, 1e-6));
        assert!(close(SPARC.q_star(), 3.058_679, 1e-6) && close(SPARC.q95(), 3.620_936, 1e-6));
        assert!(ITER.q_cylindrical() > 1.0 && SPARC.q_cylindrical() > 1.0);
    }

    #[test]
    fn density_and_beta_limits() {
        assert!(close(ITER.greenwald(), 1.193_662e20, 1e-6) && close(0.85 * ITER.greenwald(), 1.01e20, 5e-3));
        assert!(close(SPARC.greenwald(), 8.523_533e20, 1e-6) && close(0.37 * SPARC.greenwald(), 3.15e20, 2e-3));
        assert!(close(ITER.troyon_beta(1.8), 0.025_472, 1e-4) && close(SPARC.troyon_beta(1.0), 0.012_511, 1e-4));
    }

    #[test]
    fn confinement() {
        // No source tabulates the net heating power, so test the inverse: ITER's published 3.7 s
        // needs 83.8 MW, SPARC's 0.77 s needs 29.0 MW.
        assert!(close(ITER.ipb98(10.1, 83.815, 2.5, 1.7, 1.0), 3.7, 1e-4));
        let n19 = 0.37 * SPARC.greenwald() / 1e19;
        assert!(close(SPARC.ipb98(n19, 29.036, 2.5, 1.75, 1.0), 0.77, 1e-4));
        // τ falls as P^−0.69.
        assert!(close(ITER.ipb98(10.1, 2.0 * 83.815, 2.5, 1.7, 1.0) / 3.7, 2f64.powf(-0.69), 1e-4));
    }
}
