//! Plasmoids: how magnetised plasma organises itself. Taylor relaxation, spheromaks,
//! field-reversed configurations, magnetic reconnection, and the energy a free plasmoid can hold.
//!
//! A turbulent magnetised plasma sheds energy faster than magnetic helicity, so it relaxes to the
//! minimum-energy state at fixed helicity: `∇×B = λB` (Woltjer, Taylor). That state, a spheromak
//! or reversed-field pinch, is what the plasma chooses by itself. Reconnection, by which field
//! lines break and rejoin, runs at the slow Sweet–Parker rate `S^(−1/2)` in small systems but at
//! about 0.01 once the current sheet breaks into plasmoids.
//!
//! ```
//! use plasmoid::relaxation::{SPHERE_EIGENVALUE, frc_beta};
//!
//! // A spheromak in a spherical flux conserver of radius a has λa = 4.4934.
//! assert!((SPHERE_EIGENVALUE - 4.4934).abs() < 1e-4);
//! // An FRC filling half its wall radius holds plasma at ⟨β⟩ = 0.875.
//! assert!((frc_beta(0.5) - 0.875).abs() < 1e-12);
//! ```

use core::f64::consts::PI;

use crate::collisions::{Resistivity, spitzer_eta};
use metrology::{M_P, MU_0 as MU0};

/// First zero of J₁, the cylinder's radial eigenvalue.
pub const J11: f64 = 3.831_705_970_207_512;
/// First zero of J₀: the Bessel-function model's field reverses at μa = 2.4048.
pub const J01: f64 = 2.404_825_557_695_773;
/// λa of a spheromak in a spherical flux conserver of radius a (first zero of j₁).
pub const SPHERE_EIGENVALUE: f64 = 4.493_409_457_909_064;

/// Bessel function Jₙ(x) by Bessel's integral.
pub fn bessel_j(n: u32, x: f64) -> f64 {
    let steps = 2000;
    let h = PI / steps as f64;
    let f = |t: f64| (f64::from(n) * t - x * t.sin()).cos();
    let inner: f64 = (1..steps).map(|k| if k % 2 == 1 { 4.0 } else { 2.0 } * f(k as f64 * h)).sum();
    (f(0.0) + inner + f(PI)) * h / 3.0 / PI
}

/// The Bessel-function model of a relaxed cylinder at `mu_a` = μa: the pinch parameter
/// Θ = B_θ(a)/⟨B_z⟩ and the reversal parameter F = B_z(a)/⟨B_z⟩.
pub fn bessel_model(mu_a: f64) -> (f64, f64) {
    (mu_a / 2.0, mu_a * bessel_j(0, mu_a) / (2.0 * bessel_j(1, mu_a)))
}

/// λa of the lowest axisymmetric Taylor state (a spheromak) in a cylinder of radius a and length
/// h: `√(j₁₁² + (πa/h)²)`.
pub fn spheromak_cylinder(h_over_a: f64) -> f64 {
    (J11 * J11 + (PI / h_over_a).powi(2)).sqrt()
}

/// Magnetic energy (J) of a Taylor state with eigenvalue `lambda` (m⁻¹) and helicity `helicity`
/// (Wb²): `λK/(2μ₀)`.
pub fn taylor_energy(lambda: f64, helicity: f64) -> f64 {
    lambda * helicity / (2.0 * MU0)
}

/// Average beta of a field-reversed configuration whose separatrix fills `xs` of the wall
/// radius: `1 − x_s²/2` (Tuszewski 1988, Eq. 10).
pub fn frc_beta(xs: f64) -> f64 {
    1.0 - xs * xs / 2.0
}

/// Bounds on an FRC's s-parameter (ion gyroradii across the plasma) for `xs` and wall radius over
/// the ion gyroradius `r_over_rho`: `S*x_s²/(4√2) < s < S*x_s/4`.
pub fn frc_s_bounds(xs: f64, r_over_rho: f64) -> (f64, f64) {
    (r_over_rho * xs * xs / (4.0 * 2f64.sqrt()), r_over_rho * xs / 4.0)
}

/// Lundquist number S and Alfvén speed (m/s) of a current sheet of length `length` (m) in plasma
/// of ion density `n` (m⁻³), electron temperature `te` (eV), field `b` (T), Coulomb logarithm
/// `ln_lambda` and ion mass number `a`, with Spitzer's parallel resistivity. The Alfvén speed counts
/// the ions' mass `a·m_p` only, 2.7×10⁻⁴ above `magnetics::mhd::alfven_speed` with the electrons'
/// mass included for hydrogen (`physica` checks the gap).
pub fn lundquist(length: f64, n: f64, te: f64, b: f64, ln_lambda: f64, a: f64) -> (f64, f64) {
    let eta = spitzer_eta(te, 1.0, ln_lambda, Resistivity::Parallel);
    let va = b / (MU0 * n * a * M_P).sqrt();
    (MU0 * length * va / eta, va)
}

/// Sweet–Parker reconnection rate `S^(−1/2)` (inflow over Alfvén speed).
pub fn sweet_parker_rate(s: f64) -> f64 {
    s.powf(-0.5)
}

/// The fast reconnection rate once the sheet breaks into plasmoids (S above ~10⁴), independent
/// of S (Loureiro et al. 2007; Bhattacharjee et al. 2009; Uzdensky et al. 2010).
pub const PLASMOID_RATE: f64 = 0.01;

/// The most energy (J) a free plasmoid of volume `volume` (m³) can hold against an ambient
/// pressure `p0` (Pa) by the virial theorem: `3p₀V`. For ball lightning in air, 0.3 MJ/m³.
pub fn virial_energy_cap(volume: f64, p0: f64) -> f64 {
    3.0 * p0 * volume
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, rel: f64) -> bool {
        (a / b - 1.0).abs() < rel
    }

    #[test]
    fn taylor_states() {
        assert!(bessel_j(0, J01).abs() < 1e-9 && bessel_j(1, J11).abs() < 1e-9);
        assert!(close(bessel_j(0, 1.0), 0.765_197_686_6, 1e-9));
        // Field reversal at μa = 2.4048, Θ = 1.2024; the helical state appears at μa = 3.11.
        let (theta, f) = bessel_model(J01);
        assert!(close(theta, 1.2024, 1e-4) && f.abs() < 1e-9);
        assert!(close(bessel_model(3.11).0, 1.555, 1e-3));
        assert!(close(spheromak_cylinder(1.0), 4.9550, 1e-4) && close(spheromak_cylinder(1.67), 4.2686, 1e-4));
        assert!(close(taylor_energy(10.0, 1e-3), 10.0 * 1e-3 / (2.0 * MU0), 1e-12));
    }

    #[test]
    fn field_reversed_configurations() {
        assert!(close(frc_beta(0.5), 0.875, 1e-12));
        let (lo, hi) = frc_s_bounds(0.5, 10.0);
        assert!(close(lo, 0.25 / 2.0 * 10.0 / (2.0 * 2f64.sqrt()), 1e-12) && close(hi, 1.25, 1e-12));
    }

    #[test]
    fn reconnection_rates() {
        // Ji & Daughton 2011: MRX ~3×10³, the corona ~10¹³, the magnetotail ~4×10¹⁵.
        let (mrx, va) = lundquist(0.2, 1e19, 10.0, 0.03, 10.0, 1.0);
        assert!(close(mrx, 3.13e3, 0.01) && close(va, 2.069e5, 0.01));
        assert!(close(lundquist(5e6, 1e15, 200.0, 2e-3, 20.0, 1.0).0, 2.33e13, 0.01));
        assert!(close(lundquist(3e8, 3e5, 600.0, 2e-8, 30.0, 1.0).0, 2.80e15, 0.01));
        // Bhattacharjee et al. 2009 Table 1.
        assert!(close(sweet_parker_rate(3.14e4), 5.6e-3, 0.01) && close(sweet_parker_rate(6.28e5), 1.26e-3, 0.01));
        // The corona: Sweet–Parker takes 200 days, the plasmoid rate 6 minutes.
        let t_a = 5e6 / lundquist(5e6, 1e15, 200.0, 2e-3, 20.0, 1.0).1;
        assert!(close(t_a / sweet_parker_rate(2.33e13) / 86_400.0, 203.0, 0.02));
        assert!(close(t_a / PLASMOID_RATE / 60.0, 6.04, 0.02));
    }

    #[test]
    fn plasmoid_energy_cap() {
        // A 20 cm ball in air holds at most 1.27 kJ; reported ball lightning is 1–10 MJ/m³.
        let v = PI * 0.2f64.powi(3) / 6.0;
        assert!(close(virial_energy_cap(v, 101_325.0), 1.273e3, 1e-3));
    }
}
