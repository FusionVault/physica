//! Gravitational lensing by galaxies, stars and dark-matter haloes: the Einstein radius of a
//! singular isothermal sphere, point-lens microlensing of point and finite sources and the mass
//! it weighs, the critical density, the NFW halo's normalisation, concentration and surface
//! density, and the Hubble constant from lensed quasars beside the other methods.
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

use metrology::{AU, C, G, SOLAR_MASS};

use crate::integrate_log;

/// `κ = 4G/(c² AU)` (rad/kg), so that a lens of mass `M` has `θ_E² = κMπ_rel`: the deflection of
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

/// A uniform disc of radius `rho` at separation `u` from a point lens, both in Einstein radii:
/// the exact magnification of [Lee et al. (2009)](https://arxiv.org/abs/0901.1316),
/// `(1/πρ²)∫₀^π [u₂√(u₂² + 4) − u₁√(u₁² + 4)]dϑ`, with `u₂` and `u₁` the far and near limb
/// along each direction `ϑ` from the lens (`u₁ = 0` when the lens is inside the disc). Their
/// substitution `u sin ϑ = ρ sin t` makes `t` the angle about the disc's centre and joins the two
/// limbs into one integral without endpoint singularities,
/// `(1/πρ)∫₀^π (ρ + u cos t)√(r² + 4)/r dt` with `r² = u² + ρ² + 2uρ cos t` (this joined form
/// is ours; it matches their integral to 10⁻¹³). At `u = 0` it is
/// `√(1 + 4/ρ²)`, with the lens on the limb at `ρ = 1` it is `1 + 2/π`, and far from the disc it
/// exceeds [`microlensing`] by `4ρ²/u⁶`. Gould's `A(u)B₀(u/ρ)` and Yoo et al.'s second order in `ρ`
/// ([2004](https://arxiv.org/abs/astro-ph/0309302)) are approximations to it: Gould's is 1.5% low
/// at `u = ρ/2` for `ρ = 0.5`. The result loses about `u/ρ` in relative precision to
/// cancellation, so for `u ≫ ρ` use [`microlensing`].
pub fn finite_source(u: f64, rho: f64) -> f64 {
    // Over w = π − t, shifted by δ so the log-grid Simpson rule crowds its points at the far limb,
    // where the integrand turns over within |u − ρ|/ρ when the lens is near the limb.
    let delta = ((u - rho).abs() / rho).clamp(1e-9, 1.0);
    let integrand = |v: f64| {
        let s2 = (0.5 * (v - delta)).sin().powi(2);
        // ρ + u cos t and r², written without cancellation near t = π.
        let (near, r2) = (rho - u + 2.0 * u * s2, (u - rho).powi(2) + 4.0 * u * rho * s2);
        if r2 == 0.0 { 0.0 } else { near * (r2 + 4.0).sqrt() / r2.sqrt() }
    };
    integrate_log(delta, PI + delta, 16000, integrand) / (PI * rho)
}

/// The critical surface density (kg/m²) of a lens at angular-diameter distance `dd` in front of
/// a source at `ds`, `dds` apart: `Σ_crit = c²D_s/(4πG D_d D_ds)`, so that `κ = Σ/Σ_crit` and
/// an Einstein ring encloses `πR_E²Σ_crit` ([Wright & Brainerd 2000](https://arxiv.org/abs/astro-ph/9908213)).
pub fn sigma_crit(dd: f64, ds: f64, dds: f64) -> f64 {
    C * C * ds / (4.0 * PI * G * dd * dds)
}

/// The characteristic overdensity of an NFW halo of concentration `c = r₂₀₀/r_s`:
/// `δ_c = (200/3)c³/[ln(1 + c) − c/(1 + c)]`, so that `ρ(r) = δ_cρ_c/[(r/r_s)(1 + r/r_s)²]` holds
/// `M₂₀₀ = (800π/3)ρ_c r₂₀₀³` inside `r₂₀₀`, with `ρ_c = 3H²(z)/8πG` ([Wright & Brainerd 2000](https://arxiv.org/abs/astro-ph/9908213)).
/// The surface densities of [`nfw`] are in units of `2r_sδ_cρ_c`.
pub fn nfw_delta_c(c: f64) -> f64 {
    200.0 / 3.0 * c.powi(3) / (c.ln_1p() - c / (1.0 + c))
}

/// The concentration `c₂₀₀` of a halo of mass `m200` (kg) at redshift `z`, for a Hubble constant
/// `h` in units of 100 km/s/Mpc: `log₁₀c = a + b log₁₀(M/10¹²h⁻¹M☉)` with `b = −0.101 + 0.026z`
/// and `a = 0.520 + 0.385 exp(−0.617z^1.21)`, fitted to Planck-cosmology simulations
/// ([Dutton & Macciò 2014](https://arxiv.org/abs/1402.7073)). Duffy et al.'s WMAP5 fit gives 6.05
/// rather than 8.04 for a 10¹² h⁻¹M☉ halo today, a measure of how uncertain concentrations are.
pub fn concentration_dutton(m200: f64, h: f64, z: f64) -> f64 {
    let (a, b) = (0.520 + 0.385 * (-0.617 * z.powf(1.21)).exp(), -0.101 + 0.026 * z);
    10f64.powf(a + b * (m200 * h / (1e12 * SOLAR_MASS)).log10())
}

/// The mass (kg) of a lens from its angular Einstein radius `theta_e` (rad) and microlensing
/// parallax `pi_e`: `θ_E/(κπ_E)`. For OGLE-2011-BLG-0462, Sahu et al. (2025) find 7.15 ± 0.83 M☉;
/// [`OGLE_2011_BLG_0462`] holds the earlier measurements.
pub fn lens_mass(theta_e: f64, pi_e: f64) -> f64 {
    theta_e / (KAPPA * pi_e)
}

/// One published weighing of a microlens.
#[non_exhaustive]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Weighing {
    /// The analysis.
    pub name: &'static str,
    /// The mass, M☉: the median and its 84th and 16th percentiles' distances.
    pub mass: (f64, f64, f64),
    /// The angular Einstein radius, mas (median).
    pub theta_e: f64,
    /// The microlensing parallax (median).
    pub pi_e: f64,
    /// Whether its authors' later analysis replaces it, as Sahu et al. (2025) say of Lam et al. 2022.
    pub superseded: bool,
    /// Where it comes from.
    pub source: &'static str,
}

/// The mass history of OGLE-2011-BLG-0462, the first isolated black hole, oldest first, from
/// Sahu et al. (2025), Table 6, which says Lam et al.'s 2023 reanalysis supersedes their 2022
/// one. That table gives Mróz et al. π_E = 0.095, their photometry-only parallax; the 7.88 M☉
/// comes from their joint fit's (π_E,N, π_E,E) = (0.073, −0.050), so 0.0885 here. Lam et al.'s
/// medians come from separate marginals, so their `θ_E/(κπ_E)` misses the mass by 2–5%.
pub const OGLE_2011_BLG_0462: &[Weighing] = &[
    Weighing { name: "Sahu et al. 2022", mass: (7.1, 1.3, 1.3), theta_e: 5.18, pi_e: 0.0894, superseded: false, source: "https://arxiv.org/abs/2201.13296" },
    Weighing { name: "Lam et al. 2022 (default weight)", mass: (3.79, 0.62, 0.57), theta_e: 3.89, pi_e: 0.12, superseded: true, source: "https://arxiv.org/abs/2202.01903" },
    Weighing { name: "Lam et al. 2022 (equal weight)", mass: (2.15, 0.67, 0.54), theta_e: 4.13, pi_e: 0.24, superseded: true, source: "https://arxiv.org/abs/2202.01903" },
    Weighing { name: "Mróz, Udalski & Gould 2022", mass: (7.88, 0.82, 0.82), theta_e: 5.68, pi_e: 0.0885, superseded: false, source: "https://arxiv.org/abs/2207.10729" },
    Weighing { name: "Lam & Lu 2023", mass: (6.03, 1.19, 1.04), theta_e: 4.79, pi_e: 0.10, superseded: false, source: "https://arxiv.org/abs/2308.03302" },
    Weighing { name: "Sahu et al. 2025", mass: (7.15, 0.83, 0.83), theta_e: 5.496, pi_e: 0.0944, superseded: false, source: "https://arxiv.org/abs/2503.07820" },
];

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

/// Time-delay `H₀` in flat ΛCDM, oldest first. TDCOSMO 2025 adds DES J0408 and WGD 2038 to the six
/// H0LiCOW lenses, with new kinematics, and calls its earlier kinematics superseded.
/// [`H0_MEASUREMENTS`] sets these beside the other methods.
pub const H0_REGISTRY: [Hubble; 2] = [
    Hubble { name: "H0LiCOW (six lenses)", value: (73.3, 1.7, 1.8), assumption: Assumption::Assertive, data: "lenses only", source: "https://arxiv.org/abs/1907.04869" },
    Hubble { name: "TDCOSMO 2025 (eight lenses, new kinematics)", value: (71.6, 3.9, 3.3), assumption: Assumption::Conservative, data: "lenses with Pantheon+ supernovae", source: "https://arxiv.org/abs/2506.03023" },
];

/// How a value of `H₀` was measured.
#[non_exhaustive]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Method {
    /// The CMB's acoustic scale, through a cosmological model.
    Cmb,
    /// Baryon acoustic oscillations, with the sound horizon from big-bang nucleosynthesis.
    Bao,
    /// Geometric anchors, then Cepheids or red giants, then type Ia supernovae.
    DistanceLadder,
    /// Lensed quasars' time delays, and how the mass-sheet degeneracy is broken.
    TimeDelay(Assumption),
    /// A gravitational-wave standard siren with its host galaxy's redshift.
    Siren,
}

/// One measurement of the Hubble constant, by any method.
#[non_exhaustive]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Measurement {
    /// The analysis.
    pub name: &'static str,
    /// The median and its 84th and 16th percentiles' distances, km/s/Mpc.
    pub value: (f64, f64, f64),
    /// How it was measured.
    pub method: Method,
    /// What it uses, and other forms of the result.
    pub data: &'static str,
    /// Where it comes from.
    pub source: &'static str,
}

const fn time_delay(h: Hubble) -> Measurement {
    Measurement { name: h.name, value: h.value, method: Method::TimeDelay(h.assumption), data: h.data, source: h.source }
}

/// `H₀` by every method, oldest first, all as the median and its 16th–84th percentiles. The
/// early-universe values (CMB, BAO) sit near 68 and the distance ladder near 73; SH0ES and Planck
/// differ by 4.85σ, which Riess et al. round to 5σ.
///
/// ```
/// use gravitational::lensing::H0_MEASUREMENTS;
///
/// let value = |name: &str| H0_MEASUREMENTS.iter().find(|m| m.name.starts_with(name)).unwrap().value;
/// let (shoes, planck) = (value("SH0ES"), value("Planck"));
/// let tension = (shoes.0 - planck.0) / shoes.2.hypot(planck.1);
/// assert!((tension - 4.85).abs() < 0.005);
/// ```
pub const H0_MEASUREMENTS: &[Measurement] = &[
    Measurement { name: "GW170817 and NGC 4993", value: (74.0, 16.0, 8.0), method: Method::Siren, data: "one siren; the abstract's 70.0 (+12.0/−8.0) is the maximum a posteriori and its narrowest 68% interval", source: "https://arxiv.org/abs/1710.05835" },
    Measurement { name: "Planck 2018 (TT,TE,EE+lowE+lensing)", value: (67.36, 0.54, 0.54), method: Method::Cmb, data: "flat ΛCDM; 67.66 ± 0.42 with BAO", source: "https://arxiv.org/abs/1807.06209" },
    time_delay(H0_REGISTRY[0]),
    Measurement { name: "SH0ES 2022", value: (73.04, 1.04, 1.04), method: Method::DistanceLadder, data: "Cepheids and type Ia supernovae, systematics included", source: "https://arxiv.org/abs/2112.04510" },
    Measurement { name: "CCHP 2025 (tip of the red-giant branch)", value: (70.39, 1.94, 1.94), method: Method::DistanceLadder, data: "24 supernova calibrators from HST and JWST; ±1.22 (stat), ±1.33 (sys) and ±0.70 (σ_SN) in quadrature", source: "https://arxiv.org/abs/2408.06153" },
    Measurement { name: "DESI DR2 (BAO + BBN)", value: (68.51, 0.58, 0.58), method: Method::Bao, data: "flat ΛCDM; the conclusions print 68.50", source: "https://arxiv.org/abs/2503.14738" },
    time_delay(H0_REGISTRY[1]),
];

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{H100, cosmology::{PLANCK18, angular_diameter_distance}};
    use metrology::MPC;

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
    fn finite_sources() {
        // Identities: A(0) = √(1 + 4/ρ²); with the lens on the limb, (1/π)[2/ρ + ((1 + ρ²)/ρ²)
        // (π/2 + arcsin((ρ² − 1)/(ρ² + 1)))], which is 1 + 2/π at ρ = 1 (Witt & Mao's case, as
        // recalled in the research notes and confirmed there against Lee et al.'s integral).
        for rho in [0.01, 0.05, 0.5, 1.0] {
            assert!((finite_source(0.0, rho) / (1.0 + 4.0 / (rho * rho)).sqrt() - 1.0).abs() < 1e-14, "{rho}");
            let limb = (2.0 / rho + (1.0 + rho * rho) / (rho * rho) * (PI / 2.0 + ((rho * rho - 1.0) / (rho * rho + 1.0)).asin())) / PI;
            assert!((finite_source(rho, rho) / limb - 1.0).abs() < 1e-12, "{rho}");
        }
        assert!((finite_source(1.0, 1.0) / (1.0 + 2.0 / PI) - 1.0).abs() < 1e-12);
        // Lee et al.'s Fig. 2 peak, A(u = 0.1; ρ = 0.5) = 4.086267 (research script, scipy quad of
        // their integral as printed), and the same quadrature either side of the limb.
        for (u, rho, a) in [(0.1, 0.5, 4.086_266_585_622_305), (0.4995, 0.5, 2.759_105_373_407_014), (0.5005, 0.5, 2.739_050_255_247_673), (1.0, 0.5, 1.381_003_927_897_093), (0.01, 0.05, 39.609_820_689_795_04)] {
            assert!((finite_source(u, rho) / a - 1.0).abs() < 1e-13, "{u} {rho}");
        }
        // Far away the disc averages A = 1 + 2/u⁴ − 8/u⁶ + …, so it exceeds Paczyński by
        // (ρ²/8)∇²A = 4ρ²/u⁶(1 − 9/u² + …) (ours).
        let (u, rho) = (30.0f64, 0.3f64);
        let excess = finite_source(u, rho) - microlensing(u).0;
        assert!((excess / (4.0 * rho * rho / u.powi(6)) - 1.0).abs() < 0.02, "{excess:e}");
    }

    #[test]
    fn nfw_normalisation() {
        // δ_c(4, 5, 10) = 5271.147547, 8694.810191, 44778.661516 (research script).
        for (c, dc) in [(4.0, 5271.147547), (5.0, 8694.810191), (10.0, 44778.661516)] {
            assert!((nfw_delta_c(c) - dc).abs() < 1e-6, "{c}");
        }
        // Dutton & Macciò at z = 0 is a pure power law: 10^0.905 = 8.03526 for 10¹² h⁻¹M☉.
        let h = 0.6766;
        assert!((concentration_dutton(1e12 / h * SOLAR_MASS, h, 0.0) - 10f64.powf(0.905)).abs() < 1e-13);
        assert!((concentration_dutton(1e14 / h * SOLAR_MASS, h, 0.5) - 4.35381).abs() < 5e-6);
        // A 10¹⁴ h⁻¹M☉ halo at z_d = 0.3 lensing a source at z_s = 1 in Planck18 (research script):
        // ρ_c = 1.180190×10⁻²⁶ kg/m³, c₂₀₀ = 4.6457, δ_c = 7361.493, r_s = 216.089 kpc,
        // Σ_crit = 5.749602 kg/m², κ_s = 2r_sδ_cρ_c/Σ_crit = 0.20151, and δ_c holds M₂₀₀ to 10⁻¹⁴.
        let (p, zd, zs) = (PLANCK18, 0.3, 1.0);
        let rho_c = 3.0 * (p.h0 * p.e(zd)).powi(2) / (8.0 * PI * G);
        let m200 = 1e14 / h * SOLAR_MASS;
        let r200 = (3.0 * m200 / (800.0 * PI * rho_c)).cbrt();
        let c = concentration_dutton(m200, h, zd);
        let (dc, rs) = (nfw_delta_c(c), r200 / c);
        assert!((rho_c / 1.180190e-26 - 1.0).abs() < 1e-6 && (c - 4.6457).abs() < 5e-5 && (dc - 7361.493).abs() < 5e-4);
        assert!((rs / (MPC / 1e3) - 216.089).abs() < 5e-4);
        assert!((4.0 * PI * rho_c * dc * rs.powi(3) * (c.ln_1p() - c / (1.0 + c)) / m200 - 1.0).abs() < 1e-14);
        let sigma = sigma_crit(p.angular_diameter_distance(0.0, zd), p.angular_diameter_distance(0.0, zs), p.angular_diameter_distance(zd, zs));
        assert!((sigma - 5.749602).abs() < 5e-7 && (2.0 * rs * dc * rho_c / sigma - 0.20151).abs() < 5e-6, "{sigma}");
    }

    #[test]
    fn ogle_2011_blg_0462() {
        // Where the published θ_E and π_E are a pair, θ_E/(κπ_E) gives back the published mass:
        // 7.1, 7.88 and 7.15 M☉.
        let mas = PI / 180.0 / 3600.0 * 1e-3;
        for w in OGLE_2011_BLG_0462.iter().filter(|w| !w.name.starts_with("Lam")) {
            let m = lens_mass(w.theta_e * mas, w.pi_e) / SOLAR_MASS;
            assert!((m - w.mass.0).abs() < 0.02, "{}: {m}", w.name);
        }
        for w in OGLE_2011_BLG_0462.iter().filter(|w| w.name.starts_with("Lam")) {
            let m = lens_mass(w.theta_e * mas, w.pi_e) / SOLAR_MASS;
            assert!((m / w.mass.0 - 1.0).abs() < 0.06, "{}: {m}", w.name);
        }
        let current: Vec<_> = OGLE_2011_BLG_0462.iter().filter(|w| !w.superseded).map(|w| w.mass.0).collect();
        assert!(current.iter().all(|&m| m > 5.0) && OGLE_2011_BLG_0462.iter().filter(|w| w.superseded).all(|w| w.mass.0 < 4.0));
    }

    #[test]
    fn the_h0_registry() {
        // H0LiCOW is recovered from its six distances by cosmology::time_delay_distance; the
        // conservative reanalysis is wider and lower.
        let [h0licow, tdcosmo] = H0_REGISTRY;
        assert!(tdcosmo.value.0 < h0licow.value.0 && tdcosmo.value.1 > h0licow.value.1);
        assert_eq!((h0licow.assumption, tdcosmo.assumption), (Assumption::Assertive, Assumption::Conservative));
        // The wider registry carries both rows unchanged, and every method.
        let lenses: Vec<_> = H0_MEASUREMENTS.iter().filter(|m| matches!(m.method, Method::TimeDelay(_))).collect();
        assert!(lenses.len() == 2 && lenses[0].value == h0licow.value && lenses[1].method == Method::TimeDelay(Assumption::Conservative));
        for method in [Method::Cmb, Method::Bao, Method::DistanceLadder, Method::Siren] {
            assert!(H0_MEASUREMENTS.iter().any(|m| m.method == method), "{method:?}");
        }
        // CCHP's total is its three errors in quadrature: 1.94.
        assert!((1.22f64.hypot(1.33).hypot(0.70) - 1.94).abs() < 0.005);
    }
}
