//! Gravitational lensing by galaxies, stars and dark-matter haloes: the Einstein radius of a
//! singular isothermal sphere, point-lens microlensing and the mass it weighs, the NFW halo's
//! surface density, and the Hubble constant from lensed quasars.
//!
//! ```
//! use gravitational::lensing::{KAPPA, lens_mass, microlensing};
//! use metrology::SOLAR_MASS;
//!
//! // OGLE-2011-BLG-0462, an isolated black hole (Sahu et al. 2025): θ_E = 5.496 mas and
//! // π_E = 0.0944 weigh 7.15 M☉.
//! let mas = core::f64::consts::PI / 180.0 / 3600.0 * 1e-3;
//! assert!((lens_mass(5.496 * mas, 0.0944) / SOLAR_MASS - 7.15).abs() < 0.01);
//! assert!((KAPPA * SOLAR_MASS / mas - 8.144).abs() < 1e-3); // mas per solar mass
//! // At one Einstein radius the source is magnified 3/√5 and shifted θ_E/3.
//! assert_eq!(microlensing(1.0), (3.0 / 5f64.sqrt(), 1.0 / 3.0));
//! ```

use core::f64::consts::PI;

use metrology::{AU, C, G};

/// `κ = 4G/(c² AU)` (rad/kg), so that a lens of mass `M` has `θ_E = κMπ_rel`: the deflection of
/// light grazing a mass at 1 AU, per unit mass. 8.144 mas per solar mass.
pub const KAPPA: f64 = 4.0 * G / (C * C * AU);

/// The Einstein radius (rad) of a singular isothermal sphere with velocity dispersion `sigma`
/// (m/s), for angular-diameter distances `dls` (lens to source) and `ds` (to the source):
/// `θ_E = 4π(σ/c)² D_ls/D_s`. Only the ratio of the distances enters, so `H₀` cancels.
pub fn sis(sigma: f64, dls: f64, ds: f64) -> f64 {
    4.0 * PI * (sigma / C).powi(2) * dls / ds
}

/// A point lens at separation `u` (in Einstein radii) from the source: the total magnification
/// `(u² + 2)/(u√(u² + 4))` (Paczyński 1986) and the shift of the light's centroid, in Einstein
/// radii, `u/(u² + 2)`, largest at `u = √2`.
pub fn microlensing(u: f64) -> (f64, f64) {
    ((u * u + 2.0) / (u * (u * u + 4.0).sqrt()), u / (u * u + 2.0))
}

/// The mass (kg) of a lens from its angular Einstein radius `theta_e` (rad) and microlensing
/// parallax `pi_e`: `θ_E/(κπ_E)`. For OGLE-2011-BLG-0462, Sahu et al. (2025) find 7.15 ± 0.83 M☉
/// and tabulate the history: Mróz et al. (2022) 7.88 ± 0.82, and Lam et al.'s 2022 2–4 M☉,
/// superseded by their 2023 6.0 (+1.2/−1.0).
pub fn lens_mass(theta_e: f64, pi_e: f64) -> f64 {
    theta_e / (KAPPA * pi_e)
}

/// An NFW halo's surface density and its mean inside radius `x` (in units of the scale radius
/// `r_s`), both in units of `2r_sδ_cρ_c` ([Wright & Brainerd 2000](https://arxiv.org/abs/astro-ph/9908213)).
/// The tangential shear is their difference over `Σ_c`. With `u = x² − 1` and
/// `F = arctan(√u)/√u` (`artanh(√−u)/√−u` inside `r_s`), the density is `(1 − F)/u` and the
/// mean `(2/x²)(F + ln(x/2))`. Near `x = 1` the density is summed as `Σ(−u)ᵏ/(2k + 3)`, which
/// the closed form would lose to cancellation; at `x = 1` it is 1/3.
pub fn nfw(x: f64) -> (f64, f64) {
    let u = x * x - 1.0;
    let f = if u > 0.0 { u.sqrt().atan() / u.sqrt() } else if u < 0.0 { (-u).sqrt().atanh() / (-u).sqrt() } else { 1.0 };
    let density = if u.abs() < 0.1 {
        // The series of (1 − F)/u: at |u| < 0.1 twenty terms are past round-off.
        (0..20).map(|k| (-u).powi(k) / (2 * k + 3) as f64).sum()
    } else {
        (1.0 - f) / u
    };
    (density, 2.0 / (x * x) * (f + (x / 2.0).ln()))
}

/// How a lensing `H₀` handles the mass-sheet degeneracy (Treu et al., as TDCOSMO 2025 uses
/// the terms).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Assumption {
    /// Power-law or composite mass profiles break the degeneracy by assumption, for about 2%.
    Assertive,
    /// The degeneracy is left free and constrained by stellar kinematics, for about 5%.
    Conservative,
}

/// One time-delay measurement of the Hubble constant.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Hubble {
    /// The analysis.
    pub name: &'static str,
    /// The median and its 84th and 16th percentiles' distances, km/s/Mpc.
    pub value: (f64, f64, f64),
    /// What breaks the mass-sheet degeneracy.
    pub assumption: Assumption,
    /// What else it uses, in flat ΛCDM.
    pub data: &'static str,
    /// Where it comes from.
    pub source: &'static str,
}

/// Time-delay `H₀` in flat ΛCDM, oldest first. TDCOSMO 2025 reanalyses the same six lenses with
/// new kinematics and calls its earlier kinematics superseded.
pub const H0_REGISTRY: [Hubble; 2] = [
    Hubble { name: "H0LiCOW (six lenses)", value: (73.3, 1.7, 1.8), assumption: Assumption::Assertive, data: "lenses only", source: "https://arxiv.org/abs/1907.04869" },
    Hubble { name: "TDCOSMO 2025 (the same six, new kinematics)", value: (71.6, 3.9, 3.3), assumption: Assumption::Conservative, data: "lenses with Pantheon+ supernovae", source: "https://arxiv.org/abs/2506.03023" },
];

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{H100, cosmology::angular_diameter_distance};
    use metrology::SOLAR_MASS;

    #[test]
    fn ngc_6505() {
        // The Euclid ring (O'Riordan et al. 2025): θ_E = 2.500″, z_l = 0.042, z_s = 0.4058. An
        // isothermal sphere needs σ = 313 km/s, against the DESI σ* = 303 ± 15 km/s, for any H₀.
        let (zl, zs, theta) = (0.042, 0.4058, 2.5 / 206_264.806_247);
        let (dls, ds) = (angular_diameter_distance(zl, zs, 0.7 * H100, 0.3), angular_diameter_distance(0.0, zs, 0.7 * H100, 0.3));
        let sigma = C * (theta * ds / (4.0 * PI * dls)).sqrt();
        assert!((sis(sigma, dls, ds) / theta - 1.0).abs() < 1e-14);
        assert!((sigma / 1e3 - 313.0).abs() < 1.0 && (sigma / 1e3 - 303.0).abs() < 15.0, "{}", sigma / 1e3);
        let other = angular_diameter_distance(zl, zs, 0.6 * H100, 0.3) / angular_diameter_distance(0.0, zs, 0.6 * H100, 0.3);
        assert!((other / (dls / ds) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn microlensing_and_the_lens_mass() {
        // Paczyński: A = 3/√5 at u = 1, A → 1/u as u → 0, and A → 1 far away; the centroid
        // shift peaks at θ_E/(2√2) at u = √2.
        assert!((microlensing(1e-6).0 * 1e-6 - 1.0).abs() < 1e-11 && (microlensing(1e3).0 - 1.0) < 1e-11);
        let peak = microlensing(2f64.sqrt()).1;
        assert!((peak - 1.0 / (2.0 * 2f64.sqrt())).abs() < 1e-15 && microlensing(1.3).1 < peak && microlensing(1.5).1 < peak);
        // OGLE-2011-BLG-0462: 7.15 ± 0.83 M☉ (Sahu et al. 2025). Its κ is the light deflection
        // 4GM☉/(c² AU) of a solar mass at 1 AU: 8.144 mas.
        let mas = PI / 180.0 / 3600.0 * 1e-3;
        assert!((lens_mass(5.496 * mas, 0.0944) / SOLAR_MASS - 7.149).abs() < 1e-3);
        assert!((KAPPA * SOLAR_MASS / mas - 8.144).abs() < 5e-4);
    }

    #[test]
    fn an_nfw_halo() {
        // 1/3 at x = 1 with slope −2/5; the series and the closed form meet at |u| = 0.1.
        assert!((nfw(1.0).0 - 1.0 / 3.0).abs() < 1e-16);
        let slope = (nfw(1.0 + 1e-6).0 - nfw(1.0 - 1e-6).0) / 2e-6;
        assert!((slope + 0.4).abs() < 1e-8);
        for edge in [1.1f64.sqrt(), 0.9f64.sqrt()] {
            let (a, b) = (nfw(edge * (1.0 - 1e-15)).0, nfw(edge * (1.0 + 1e-15)).0);
            assert!((a - b).abs() < 1e-14, "{edge}");
        }
        // The mean inside x is (2/x²)∫₀ˣ x′Σ dx′: Wright & Brainerd's two formulas agree.
        for x in [0.3, 1.0, 2.5] {
            let n = 4000;
            let h = x / n as f64;
            // Midpoint rule on x′Σ(x′), whose log singularity at 0 is integrable.
            let integral: f64 = (0..n).map(|i| { let s = (i as f64 + 0.5) * h; s * nfw(s).0 }).sum::<f64>() * h;
            assert!((2.0 / (x * x) * integral / nfw(x).1 - 1.0).abs() < 1e-5, "x = {x}");
        }
        assert!((nfw(1.0).1 - 2.0 * (1.0 + 0.5f64.ln())).abs() < 1e-15);
    }

    #[test]
    fn the_h0_registry() {
        // H0LiCOW is recovered from its six distances by cosmology::time_delay_distance; the
        // conservative reanalysis is wider and lower.
        let [h0licow, tdcosmo] = H0_REGISTRY;
        assert!(tdcosmo.value.0 < h0licow.value.0 && tdcosmo.value.1 > h0licow.value.1);
        assert_eq!((h0licow.assumption, tdcosmo.assumption), (Assumption::Assertive, Assumption::Conservative));
    }
}
