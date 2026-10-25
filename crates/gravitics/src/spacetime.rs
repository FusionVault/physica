//! Warp drives, wormholes, negative energy and Puthoff's polarizable vacuum: the closed forms
//! behind the DIA "DIRD" papers on spacetime engineering, from the mainstream literature.
//!
//! The formulas are standard general relativity and quantum field theory. A warp bubble needs a
//! large mass of negative energy, while quantum field theory allows negative energy only in tiny
//! amounts for short times (the quantum inequalities). Puthoff's polarizable-vacuum metric is a
//! genuine alternative: it matches general relativity at first order and differs at second order,
//! by 0.73 µas for light grazing the Sun, so it is testable in principle.
//!
//! ```
//! use gravitics::spacetime::alcubierre_energy;
//! use metrology::SOLAR_MASS;
//!
//! // A 100 m warp bubble at light speed with a 1 m wall needs half a Sun of negative energy.
//! println!("{:.2} solar masses", alcubierre_energy(100.0, 1.0, 1.0) / SOLAR_MASS); // -0.56
//! ```

use core::f64::consts::PI;

use metrology::{C, G, HBAR, K_B};

/// Total energy of an Alcubierre warp bubble of `radius` (m) and wall thickness `wall` (m) moving
/// at `speed` (in units of c), as a mass (kg): `−(v²/12)(R²/Δ)·c²/G`
/// ([Pfenning & Ford 1997](https://arxiv.org/abs/gr-qc/9702026), Eq. 28). Negative. (Their text's
/// worked numbers, "6.2×10⁶⁵ g" and "a quarter of a solar mass", disagree with their Eq. 28,
/// which gives 6.9×10⁶² kg and 0.56 solar masses for the same inputs.)
pub fn alcubierre_energy(radius: f64, wall: f64, speed: f64) -> f64 {
    -speed * speed / 12.0 * radius * radius / wall * C * C / G
}

/// Forward's "antigravity" in general relativity: the acceleration `4Gm a_f/(c²r)` (m/s²) induced
/// at distance `r` (m) by mass `m` (kg) accelerating at `a_flow` (m/s²) ([Forward 1963](https://pubs.aip.org/aapt/ajp/article/31/3/166/1046748/Guidelines-to-Antigravity)).
/// The coupling G/c² = 7.4×10⁻²⁸ m/kg makes it real and useless: one g needs m/r ≈ 3×10²¹ kg/m
/// even at a_f = 10⁶ m/s².
pub fn induced_acceleration(m: f64, a_flow: f64, r: f64) -> f64 {
    4.0 * G * m * a_flow / (C * C * r)
}

/// The most negative energy density (J/m³) quantum field theory allows when averaged over a time
/// `tau` (s): `−3ħc/(32π²(cτ)⁴)` (the quantum inequality, Pfenning & Ford 1997 eq. 9).
pub fn quantum_inequality(tau: f64) -> f64 {
    -3.0 * HBAR * C / (32.0 * PI * PI * (C * tau).powi(4))
}

/// Casimir energy density (J/m³) between ideal plates a distance `d` (m) apart: `−π²ħc/(720d⁴)`.
pub fn casimir_energy_density(d: f64) -> f64 {
    -PI * PI * HBAR * C / (720.0 * d.powi(4))
}

/// Casimir pressure (Pa) between ideal plates `d` (m) apart: `−π²ħc/(240d⁴)` (attractive).
pub fn casimir_pressure(d: f64) -> f64 {
    -PI * PI * HBAR * C / (240.0 * d.powi(4))
}

/// Casimir force (N) between a sphere of `radius` (m) and a plate `d` (m) away, in the proximity
/// force approximation: `−π³ħcR/(360d³)`.
pub fn casimir_sphere_plate(radius: f64, d: f64) -> f64 {
    -PI.powi(3) * HBAR * C * radius / (360.0 * d.powi(3))
}

/// Radial tension (Pa) at the throat of a Morris–Thorne wormhole of throat radius `b0` (m):
/// `c⁴/(8πG b₀²)`. Holding the throat open needs this much negative energy.
pub fn wormhole_throat_tension(b0: f64) -> f64 {
    C.powi(4) / (8.0 * PI * G * b0 * b0)
}

/// Line tension (J/m) of the exotic edges of Visser's cubical wormhole: `−c⁴/(8G)`.
pub fn wormhole_edge_tension() -> f64 {
    -C.powi(4) / (8.0 * G)
}

/// Temperature (K) of a horizon of surface gravity `kappa` (s⁻¹): `ħκ/(2πk_B)`. With κ = a/c it is
/// the Unruh temperature, with c³/(4GM) Hawking's, with the Hubble rate Gibbons–Hawking's.
///
/// An optical pulse in a fibre makes a horizon for light whose `κ` is the steepness `α` of its
/// refractive-index front ([Philbin et al. 2008](https://arxiv.org/abs/0711.4796)): a shock as
/// steep as twice the carrier frequency, 8×10¹⁴ Hz, "reaches 10³ K".
///
/// ```
/// use gravitics::spacetime::horizon_temperature;
///
/// assert!((horizon_temperature(8e14) - 972.5).abs() < 0.1);
/// ```
pub fn horizon_temperature(kappa: f64) -> f64 {
    HBAR * kappa / (2.0 * PI * K_B)
}

/// Unruh temperature (K) seen by an observer accelerating at `a` (m/s²).
pub fn unruh_temperature(a: f64) -> f64 {
    horizon_temperature(a / C)
}

/// Hawking temperature (K) of a black hole of mass `m` (kg).
pub fn hawking_temperature(m: f64) -> f64 {
    horizon_temperature(C.powi(3) / (4.0 * G * m))
}

/// Gibbons–Hawking temperature (K) of the cosmological horizon for Hubble rate `h` (s⁻¹).
pub fn gibbons_hawking_temperature(h: f64) -> f64 {
    horizon_temperature(h)
}

/// Evaporation time (s) of a black hole of mass `m` (kg) radiating photons only:
/// `5120πG²M³/(ħc⁴)`.
pub fn hawking_lifetime(m: f64) -> f64 {
    5120.0 * PI * G * G * m.powi(3) / (HBAR * C.powi(4))
}

/// Zero-point energy density (J/m³) of the electromagnetic field up to a cut-off `omega_c`
/// (rad/s): `ħω_c⁴/(8π²c³)`.
pub fn zero_point_energy_density(omega_c: f64) -> f64 {
    HBAR * omega_c.powi(4) / (8.0 * PI * PI * C.powi(3))
}

/// A static, spherical metric.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Metric {
    /// Schwarzschild.
    GeneralRelativity,
    /// Puthoff's polarizable vacuum, `K = exp(2GM/rc²)`
    /// ([Puthoff 2002](https://arxiv.org/abs/gr-qc/9909037)).
    PolarizableVacuum,
}

/// Deflection (rad) of light passing a `mass` (kg) at impact parameter `b` (m), to second order:
/// `4m/b + c₂(m/b)²` with `m = GM/c²`, where `c₂` is `15π/4` in general relativity and `4π` in the
/// polarizable vacuum.
pub fn light_bending(mass: f64, b: f64, metric: Metric) -> f64 {
    let x = G * mass / (C * C * b);
    let second = match metric {
        Metric::GeneralRelativity => 15.0 * PI / 4.0,
        Metric::PolarizableVacuum => 4.0 * PI,
    };
    4.0 * x + second * x * x
}

/// The exact deflection (rad) of light passing a `mass` (kg) at impact parameter `b` (m), at any
/// strength of field, or `None` when the light is captured: at or inside the critical impact
/// parameter, `3√3 GM/c²` in general relativity and `2e GM/c²` in the polarizable vacuum.
///
/// It is [Darwin's](https://doi.org/10.1098/rspa.1959.0015) integral `α = 2∫₀^{u₀} du/√F(u) − π` over `u = GM/(c²r)`, with
/// `F = (c²b/GM)⁻² − u²(1 − 2u)` in Schwarzschild's radius or `e^{4u}(c²b/GM)⁻² − u²` in
/// Puthoff's isotropic one, and the turning point `u₀` the root of `F`. Writing
/// `u = u₀(1 − t²)` divides that root out exactly, and panels graded toward `t = 0` follow the
/// integrand's peak as the ray nears the photon sphere, where the deflection grows as
/// `−ln(b/b_c − 1) + ln[216(7 − 4√3)] − π` ([Bozza 2002](https://arxiv.org/abs/gr-qc/0208075)).
/// Far out it is [`light_bending`] plus `(128/3)(GM/c²b)³`.
///
/// The photon ring's subrings are images whose light went round `n` half-orbits more. Each
/// half-orbit adds π to the deflection, so with Bozza's slope of 1 each subring sits `e^π` = 23.1
/// times closer to `b_c`: the Lyapunov exponent `γ = π` of a Schwarzschild black hole
/// ([Johnson et al. 2020](https://arxiv.org/abs/1907.04329)).
///
/// ```
/// use gravitics::spacetime::{Metric, exact_deflection};
/// use metrology::{C, G, SOLAR_MASS};
///
/// let rg = G * SOLAR_MASS / (C * C);
/// let critical = 27f64.sqrt() * rg;
/// // The impact parameter at which light turns by `angle`, by bisection.
/// let image = |angle: f64| {
///     let (mut lo, mut hi) = (critical, 2.0 * critical);
///     for _ in 0..200 {
///         let mid = 0.5 * (lo + hi);
///         if exact_deflection(SOLAR_MASS, mid, Metric::GeneralRelativity).map_or(true, |a| a > angle) { lo = mid } else { hi = mid }
///     }
///     0.5 * (lo + hi)
/// };
/// let ratio = (image(4.0 * std::f64::consts::PI) - critical) / (image(5.0 * std::f64::consts::PI) - critical);
/// assert!((ratio / std::f64::consts::PI.exp() - 1.0).abs() < 2e-5); // corrections of order e^(−4π)
/// ```
///
/// It is accurate to round-off except
/// close to `b_c`, where the deflection itself is ill-conditioned: `dα/d ln b = −1/(b/b_c − 1)`,
/// so the rounding of `b` alone moves it by `10⁻¹⁶/(b/b_c − 1)`.
pub fn exact_deflection(mass: f64, b: f64, metric: Metric) -> Option<f64> {
    let b = b * C * C / (G * mass);
    // u²(1 − 2u) or u²e^(−4u) reaches its maximum 1/b_c² at the photon sphere.
    let (g, photon): (fn(f64) -> f64, f64) = match metric {
        Metric::GeneralRelativity => (|u| u * u * (1.0 - 2.0 * u), 1.0 / 3.0),
        Metric::PolarizableVacuum => (|u| u * u * (-4.0 * u).exp(), 0.5),
    };
    if b * b * g(photon) <= 1.0 {
        return None;
    }
    let (mut lo, mut hi) = (0.0, photon);
    while hi - lo > f64::EPSILON * hi {
        let mid = 0.5 * (lo + hi);
        if b * b * g(mid) < 1.0 { lo = mid } else { hi = mid }
    }
    let u0 = 0.5 * (lo + hi);
    // (F(u) − F(u₀))/(u₀ − u), with the divided difference of e^{4u} taken by expm1.
    let q = |u: f64| match metric {
        Metric::GeneralRelativity => (u + u0) - 2.0 * (u * u + u * u0 + u0 * u0),
        Metric::PolarizableVacuum => {
            let d = u - u0;
            let slope = if d == 0.0 { 4.0 } else { (4.0 * d).exp_m1() / d };
            (u + u0) - (4.0 * u0).exp() * slope / (b * b)
        }
    };
    // Panels halving toward t = 0 keep the integrand's near-singularity, at the peak's width, a
    // panel away from each panel, so 16 points reach round-off for any width down to 2⁻⁴⁰.
    let (x, w) = geopotential::gauss_legendre(16);
    let mut sum = 0.0;
    let mut upper = 1.0;
    for k in 1..=41 {
        let lower = if k == 41 { 0.0 } else { 0.5f64.powi(k) };
        let (mid, half) = (0.5 * (upper + lower), 0.5 * (upper - lower));
        sum += half * x.iter().zip(&w).map(|(x, w)| w * 2.0 * u0.sqrt() / q(u0 * (1.0 - (mid + half * x).powi(2))).sqrt()).sum::<f64>();
        upper = lower;
    }
    Some(2.0 * sum - PI)
}

/// The deflection (rad) of light by a body's whole gravity field `field`, in the weak-field Born
/// approximation: the transverse attraction integrated along the undeflected straight ray,
/// `α = ((1 + γ)/c²)∫∇⊥U dl`, with `γ` the PPN parameter (1 in general relativity). The ray
/// passes through `point` with direction `dir`, both in the body-fixed frame of the field
/// (metres, z along the spin axis). The result is the deflection vector, perpendicular to the
/// ray and pointing toward the mass. For a ray in the equator of an oblate body it is
/// `(4GM/c²b)(1 + J₂R²/b²)`.
pub fn multipole_deflection(field: &geopotential::Coefficients, point: [f64; 3], dir: [f64; 3], gamma: f64) -> [f64; 3] {
    let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let n = dir.map(|x| x / dot(dir, dir).sqrt());
    let along = dot(point, n);
    let closest: [f64; 3] = core::array::from_fn(|i| point[i] - along * n[i]);
    let b = dot(closest, closest).sqrt();
    // s = b tan θ maps the line onto (−π/2, π/2), where the integrand is smooth: b sec²θ ∇⊥U.
    // An even number of nodes keeps them off θ = 0, which may lie on the spin axis.
    let (x, w) = geopotential::gauss_legendre(128);
    let mut sum = [0.0; 3];
    for (x, w) in x.iter().zip(&w) {
        let theta = core::f64::consts::FRAC_PI_2 * x;
        let s = b * theta.tan();
        let p: [f64; 3] = core::array::from_fn(|i| closest[i] + s * n[i]);
        let r = dot(p, p).sqrt();
        let (lat, lon) = ((p[2] / r).asin(), p[1].atan2(p[0]));
        let [north, east, up] = field.partials(r, lat, lon).acceleration_neu(r, lat);
        let ((sf, cf), (sl, cl)) = (lat.sin_cos(), lon.sin_cos());
        let g = [-sf * cl * north - sl * east + cf * cl * up, -sf * sl * north + cl * east + cf * sl * up, cf * north + sf * up];
        let gn = dot(g, n);
        let weight = core::f64::consts::FRAC_PI_2 * w * b / (theta.cos() * theta.cos());
        for i in 0..3 {
            sum[i] += weight * (g[i] - gn * n[i]);
        }
    }
    sum.map(|x| (1.0 + gamma) * x / (C * C))
}

/// The angular diameter (rad) of the shadow of a non-spinning black hole of mass `mass` (kg) at
/// distance `distance` (m): `6√3 GM/(c²D)`, the critical impact parameter `3√3 GM/c²` seen from
/// `D`. Spin and inclination make a Kerr shadow smaller, by up to about 7.5%
/// ([EHT Sgr A* Paper VI](https://arxiv.org/abs/2311.09484)).
pub fn shadow_diameter(mass: f64, distance: f64) -> f64 {
    6.0 * 3f64.sqrt() * G * mass / (C * C * distance)
}

/// The fractional deviation `δ = d_sh/(6√3 θ_g) − 1` of a shadow of angular diameter `diameter`
/// from a Schwarzschild black hole of angular gravitational radius `theta_g = GM/(c²D)` (both in
/// the same angular unit), as the EHT defines it. Kerr black holes give `δ` from about −0.075 to 0.
pub fn shadow_deviation(diameter: f64, theta_g: f64) -> f64 {
    diameter / (6.0 * 3f64.sqrt() * theta_g) - 1.0
}

/// The shadow radius, in units of `GM/c²`, that a Kerr black hole of spin `a` (−1 to 1) shows an
/// observer on its spin axis. Its edge is the spherical photon orbit with no angular momentum
/// about the axis, the largest root of `r³ − 3r² + a²r + a² = 0`, where Teo's constants
/// (Gen. Rel. Grav. 35, 1909, 2003) give the radius `√(η + a²) = √(8r³/(r² − 1))`. `√27` without
/// spin and `2(1 + √2)` at `a = 1`.
pub fn kerr_pole_on_shadow(a: f64) -> f64 {
    // Newton from r = 3 descends monotonically: the cubic is rising and convex for r > 1.
    let mut r = 3.0;
    for _ in 0..100 {
        let step = (r * r * r - 3.0 * r * r + a * a * r + a * a) / (3.0 * r * r - 6.0 * r + a * a);
        r -= step;
        if step.abs() <= 1e-16 * r {
            break;
        }
    }
    (8.0 * r * r * r / (r * r - 1.0)).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use metrology::PLANCK_LENGTH;

    fn close(a: f64, b: f64, rel: f64) -> bool {
        (a / b - 1.0).abs() < rel
    }

    #[test]
    fn warp_bubbles_and_quantum_inequalities() {
        // R = 100 m with a wall of 100 Planck lengths (v = c): 6.9×10⁶² kg (the paper: 6.2×10⁶²).
        assert!(close(alcubierre_energy(100.0, 100.0 * PLANCK_LENGTH, 1.0), -6.94e62, 5e-3));
        // A 1 m wall: 0.56 solar masses.
        assert!(close(alcubierre_energy(100.0, 1.0, 1.0) / 1.989e30, -0.564, 5e-3));
        assert!(close(quantum_inequality(1e-9), -3.72e-26, 5e-3));
        assert!(close(quantum_inequality(1e-15), -3.72e-2, 5e-3));
    }

    #[test]
    fn casimir_closed_forms() {
        assert!(close(casimir_pressure(1e-6), -1.300e-3, 1e-3));
        // Pressure is −dE/dd of the energy per area u·d, and the sphere force is 2πR·u(d)·d.
        let d = 1e-6;
        let energy = |d: f64| casimir_energy_density(d) * d;
        let h = 1e-12;
        assert!(close(-(energy(d + h) - energy(d - h)) / (2.0 * h), casimir_pressure(d), 1e-6));
        assert!(close(casimir_sphere_plate(1e-4, d), 2.0 * PI * 1e-4 * energy(d), 1e-12));
    }

    #[test]
    fn wormholes() {
        assert!(close(wormhole_throat_tension(10.0), 4.82e40, 5e-3));
        assert!(close(wormhole_edge_tension(), -1.513e43, 1e-3));
    }

    #[test]
    fn horizons() {
        // Unruh: 4.055×10⁻²¹ K per m/s²; 1 K needs 2.466×10²⁰ m/s².
        assert!(close(unruh_temperature(1.0), 4.055e-21, 1e-3));
        assert!(close(1.0 / unruh_temperature(1.0), 2.466e20, 1e-3));
        // A solar-mass black hole: 61.7 nK, evaporating (photons only) in 2.1×10⁶⁷ years.
        let sun = metrology::SOLAR_MASS;
        assert!(close(hawking_temperature(sun), 6.17e-8, 2e-3));
        assert!(close(hawking_lifetime(sun) / 3.155_76e7, 2.095e67, 3e-3));
        // The cosmological horizon at H₀ = 67.4 km/s/Mpc: 2.7×10⁻³⁰ K.
        assert!(close(gibbons_hawking_temperature(67.4e3 / 3.085_677_581e22), 2.655e-30, 2e-3));
        // Zero-point energy to the electron's Compton frequency: 1.8×10²² J/m³.
        let compton = metrology::M_E * C * C / HBAR;
        assert!(close(zero_point_energy_density(compton), 1.801e22, 2e-3));
    }

    #[test]
    fn schwinger_is_a_boltzmann_factor_at_twice_the_unruh_temperature() {
        // πE_S/E = m_ec²/(2k_B T_U) for an electron accelerated by E: exactly, at any field.
        use crate::strong_field::{E_SCHWINGER, schwinger_exponent};
        let (m_e, e) = (metrology::M_E, metrology::E_CHARGE);
        for f in [1.0, 0.1, 1e-3] {
            let field = f * E_SCHWINGER;
            let t_u = unruh_temperature(e * field / m_e);
            assert!(close(schwinger_exponent(field), m_e * C * C / (2.0 * K_B * t_u), 1e-12));
        }
        // At the critical field itself, T_U = 9.44×10⁸ K (81.3 keV).
        assert!(close(unruh_temperature(e * E_SCHWINGER / m_e), 9.438e8, 1e-3));
    }

    #[test]
    fn polarizable_vacuum_against_general_relativity() {
        let (sun, limb) = (metrology::SOLAR_MASS, 6.957e8);
        let arcsec = 180.0 / PI * 3600.0;
        let gr = light_bending(sun, limb, Metric::GeneralRelativity);
        assert!((gr * arcsec - 1.7512).abs() < 1e-4);
        // The two agree at first order and differ by (π/4)(m/b)² = 0.73 µas.
        let diff = light_bending(sun, limb, Metric::PolarizableVacuum) - gr;
        assert!(close(diff * arcsec * 1e6, 0.730, 5e-3));
    }

    #[test]
    fn forward_antigravity() {
        assert!((G / (C * C) / 7.426e-28 - 1.0).abs() < 1e-3);
        assert!((induced_acceleration(1e4, 1e4, 1.0) / 2.97046e-19 - 1.0).abs() < 1e-5);
        assert!((9.8 / induced_acceleration(1.0, 1e6, 1.0) / 3.29915e21 - 1.0).abs() < 1e-5);
    }

    #[test]
    fn the_shadow_of_sgr_a_star() {
        // EHT Paper VI's stellar-orbit priors: VLTI M = 4.297×10⁶ M☉ at 8277 pc and Keck
        // 3.951×10⁶ M☉ at 7935 pc give θ_g = 5.125 and 4.92 µas, so 6√3θ_g = 53.3 and 51.1 µas.
        let uas = PI / 180.0 / 3600.0 * 1e-6;
        let pc = metrology::PARSEC;
        let vlti = shadow_diameter(4.297e6 * metrology::SOLAR_MASS, 8277.0 * pc) / uas;
        let keck = shadow_diameter(3.951e6 * metrology::SOLAR_MASS, 7935.0 * pc) / uas;
        assert!((vlti / (6.0 * 3f64.sqrt() * 5.125) - 1.0).abs() < 1e-3 && (keck / (6.0 * 3f64.sqrt() * 4.92) - 1.0).abs() < 2e-3);
        assert!((vlti - 53.3).abs() < 0.05 && (keck - 51.1).abs() < 0.1, "{vlti} {keck}");
        // Its Table 1 shadow, 48.9 µas, against each prior: the fiducial δ = −0.08 and −0.04.
        assert!((shadow_deviation(48.9, 5.125) + 0.08).abs() < 0.005 && (shadow_deviation(48.9, 4.92) + 0.04).abs() < 0.005);
    }

    #[test]
    fn a_kerr_shadow_seen_from_the_pole() {
        assert!(close(kerr_pole_on_shadow(0.0), 27f64.sqrt(), 1e-15) && close(kerr_pole_on_shadow(1.0), 2.0 * (1.0 + 2f64.sqrt()), 1e-15));
        // Against Teo's η = r³(4a² − r(r − 3)²)/(a²(r − 1)²) at the same orbit, and inside the
        // EHT's Kerr range: from a = 0 to 1 the pole-on shadow shrinks by 7.1%.
        for a in [0.1, 0.5, 0.94, -0.7] {
            // The orbit by bisection on [1 + √2, 3], where the cubic changes sign.
            let cubic = |r: f64| r * r * r - 3.0 * r * r + a * a * r + a * a;
            let (mut lo, mut hi) = (1.0 + 2f64.sqrt(), 3.0);
            for _ in 0..200 {
                let mid = 0.5 * (lo + hi);
                if cubic(mid) > 0.0 { hi = mid } else { lo = mid }
            }
            let r = 0.5 * (lo + hi);
            let eta = r.powi(3) * (4.0 * a * a - r * (r - 3.0).powi(2)) / (a * a * (r - 1.0).powi(2));
            let shadow = kerr_pole_on_shadow(a);
            assert!(close(shadow, (eta + a * a).sqrt(), 1e-13), "a = {a}");
            assert!(shadow < 27f64.sqrt() && shadow > 2.0 * (1.0 + 2f64.sqrt()));
        }
        assert!((kerr_pole_on_shadow(1.0) / 27f64.sqrt() - 1.0 + 0.0708).abs() < 1e-3);
    }

    #[test]
    fn exact_deflection_against_darwin_bozza_and_the_series() {
        let m = metrology::SOLAR_MASS;
        let rg = G * m / (C * C);
        let at = |b: f64| exact_deflection(m, b * rg, Metric::GeneralRelativity).unwrap();
        // Darwin's closed form, as Iyer & Petters (2007) print it, from the closest approach r₀:
        // 4√(r₀/Q)[K(k) − F(Ψ, k)] − π, Q = √((r₀ − 2)(r₀ + 6)), k² = (Q − r₀ + 6)/(2Q),
        // sin²Ψ = (Q + 2 − r₀)/(Q + 6 − r₀); F through Carlson's R_F.
        // (Closer to the photon sphere α is ill-conditioned in b: the round trip r₀ → b → r₀
        // alone costs 10⁻¹⁶/(b/b_c − 1). Bozza's limit below covers that regime.)
        for r0 in [3.05f64, 3.5, 4.0, 6.0, 10.0, 100.0] {
            let q = ((r0 - 2.0) * (r0 + 6.0)).sqrt();
            let k2 = (q - r0 + 6.0) / (2.0 * q);
            let s2 = (q + 2.0 - r0) / (q + 6.0 - r0);
            let (big_k, _) = geopotential::special::complete_elliptic(k2);
            let f = s2.sqrt() * geopotential::special::carlson_rf(1.0 - s2, 1.0 - k2 * s2, 1.0);
            let darwin = 4.0 * (r0 / q).sqrt() * (big_k - f) - PI;
            let b = (r0.powi(3) / (r0 - 2.0)).sqrt();
            // Darwin's form subtracts π from nearly π and K from F near the photon sphere: 10⁻¹³.
            assert!((at(b) - darwin).abs() < 2e-13, "r₀ = {r0}: {} vs {darwin}", at(b));
        }
        // Bozza's strong-field limit, b̄ = ln[216(7 − 4√3)] − π = −0.4002, approached linearly:
        // 2.5×10⁻², 3.1×10⁻³, 3.8×10⁻⁴ at b/b_c − 1 = 10⁻², 10⁻³, 10⁻⁴.
        let bbar = (216.0 * (7.0 - 4.0 * 3f64.sqrt())).ln() - PI;
        assert!((bbar + 0.4002).abs() < 5e-5);
        for (eps, gap) in [(1e-2, 2.5e-2), (1e-3, 3.1e-3), (1e-4, 3.8e-4)] {
            let d = at(27f64.sqrt() * (1.0 + eps)) + eps.ln() - bbar;
            assert!((d / gap - 1.0).abs() < 0.03, "{eps}: {d:e}");
        }
        // Iyer & Petters' series in h = GM/(c²r₀) at h = 10⁻³, to its h⁵ term.
        let h = 1e-3;
        let series = 4.0 * h + (-4.0 + 15.0 * PI / 4.0) * h * h + (122.0 / 3.0 - 7.5 * PI) * h.powi(3) + (-130.0 + 3465.0 * PI / 64.0) * h.powi(4)
            + (778.3 - 3465.0 * PI / 16.0) * h.powi(5);
        assert!((at((1e9f64 / 998.0).sqrt()) - series).abs() < 1e-15, "{:e}", at((1e9f64 / 998.0).sqrt()) - series);
        // Far out, light_bending plus (128/3)(GM/c²b)³: the 2.7×10⁻⁸ beyond second order at 2×10⁴.
        let b = 2e4;
        let rest = at(b) - light_bending(m, b * rg, Metric::GeneralRelativity);
        assert!((rest / (128.0 / 3.0 / b.powi(3)) - 1.0).abs() < 2e-3, "{rest:e}");
        assert!(((rest / at(b)) / 2.67e-8 - 1.0).abs() < 2e-3);
        // Inside the photon sphere's critical impact parameter the light is captured.
        assert!(exact_deflection(m, 27f64.sqrt() * rg * (1.0 - 1e-12), Metric::GeneralRelativity).is_none());
    }

    #[test]
    fn exact_deflection_in_the_polarizable_vacuum() {
        let m = metrology::SOLAR_MASS;
        let rg = G * m / (C * C);
        // Far out it is light_bending's 4/b + 4π/b²; its photon sphere is at isotropic r = 2M, so
        // light is captured inside b = 2e, and it bends more than in general relativity.
        let b = 2e4;
        let pv = exact_deflection(m, b * rg, Metric::PolarizableVacuum).unwrap();
        assert!((pv / light_bending(m, b * rg, Metric::PolarizableVacuum) - 1.0).abs() < 1e-7);
        let critical = 2.0 * core::f64::consts::E;
        assert!(exact_deflection(m, critical * rg * (1.0 - 1e-9), Metric::PolarizableVacuum).is_none());
        assert!(exact_deflection(m, critical * rg * (1.0 + 1e-6), Metric::PolarizableVacuum).unwrap() > 10.0);
        let gr = exact_deflection(m, 10.0 * rg, Metric::GeneralRelativity).unwrap();
        assert!(exact_deflection(m, 10.0 * rg, Metric::PolarizableVacuum).unwrap() > gr);
    }

    #[test]
    fn jupiters_oblateness_in_the_deflection() {
        // Jupiter: the IAU 2015 nominal GM and equatorial radius, and Juno's J₂ = 14696.572×10⁻⁶
        // (Iess et al. 2018, as tabulated in arXiv:2310.02834).
        let (gm, radius, j2) = (1.266_865_3e17, 7.1492e7, 14_696.572e-6);
        let mut jupiter = geopotential::Coefficients::zeros(gm, radius, 2);
        jupiter.set(0, 0, 1.0, 0.0);
        jupiter.set(2, 0, -j2 / 5f64.sqrt(), 0.0);
        let mas = PI / 180.0 / 3600.0 * 1e-3;
        for (b, lon) in [(radius, 0.0), (2.0 * radius, 0.7), (radius * 1.1, -2.0)] {
            // An equatorial ray, at any longitude: (4GM/c²b)(1 + J₂R²/b²).
            let point = [b * f64::cos(lon), b * f64::sin(lon), 0.0];
            let dir = [-f64::sin(lon), f64::cos(lon), 0.0];
            let alpha = multipole_deflection(&jupiter, point, dir, 1.0);
            let want = 4.0 * gm / (C * C * b) * (1.0 + j2 * radius * radius / (b * b));
            let size = (alpha[0] * alpha[0] + alpha[1] * alpha[1] + alpha[2] * alpha[2]).sqrt();
            assert!(close(size, want, 1e-13) && alpha[2].abs() < 1e-15 * size, "b = {b}");
            // It points at the planet.
            assert!(close(-(alpha[0] * point[0] + alpha[1] * point[1]) / (b * size), 1.0, 1e-13));
        }
        // At the limb: 16.27 mas, of which J₂ adds 239 µas.
        let limb = 4.0 * gm / (C * C * radius);
        assert!((limb / mas - 16.27).abs() < 0.005 && (limb * j2 / mas * 1e3 - 239.0).abs() < 0.5);
        // A point mass, off the axis in any direction, gives the first order of light_bending.
        let mut point = geopotential::Coefficients::zeros(gm, radius, 0);
        point.set(0, 0, 1.0, 0.0);
        let alpha = multipole_deflection(&point, [1e8, 2e8, -0.5e8], [0.3, -0.2, 0.9], 1.0);
        let n = [0.3, -0.2, 0.9].map(|x: f64| x / 0.94f64.sqrt());
        let p = [1e8, 2e8, -0.5e8];
        let along = p[0] * n[0] + p[1] * n[1] + p[2] * n[2];
        let b = (0..3).map(|i| (p[i] - along * n[i]).powi(2)).sum::<f64>().sqrt();
        let size = (alpha[0] * alpha[0] + alpha[1] * alpha[1] + alpha[2] * alpha[2]).sqrt();
        assert!(close(size, 4.0 * gm / (C * C * b), 1e-13));
        // γ scales it as (1 + γ)/2.
        let half = multipole_deflection(&point, [1e8, 2e8, -0.5e8], [0.3, -0.2, 0.9], 0.0);
        assert!(close(half[0], alpha[0] / 2.0, 1e-14));
    }
}
