//! Warp drives, wormholes, negative energy and Puthoff's polarizable vacuum: the closed forms
//! behind the DIA "DIRD" papers on spacetime engineering, from the mainstream literature.
//!
//! The formulas are standard general relativity and quantum field theory. A warp bubble needs a
//! large mass of negative energy, while quantum field theory allows negative energy only in tiny
//! amounts for short times (the quantum inequalities).
//!
//! Puthoff's polarizable vacuum (PV) is two things, and they fare differently. Its static
//! geometry, the exponential metric `K = exp(2GM/rc²)`, matches general relativity at first
//! order and differs at second (0.73 µas for light grazing the Sun) and in the strong field,
//! where its shadow is 4.6% larger, 1–2σ from the EHT's Sgr A* results: **speculative**. Its
//! dynamics, Puthoff's scalar theory, radiates only 2/3 of general relativity's power from a
//! binary ([`pv_radiated_power_ratio`]), about 200σ from the Hulse–Taylor pulsar: **excluded**
//! ([Ibison 2003](https://arxiv.org/abs/astro-ph/0302273)).
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
/// at `speed` (in units of c), as a mass (kg), for Pfenning & Ford's linear wall:
/// `−(v²/12)(R²/Δ + Δ/12)·c²/G` ([Pfenning & Ford 1997](https://arxiv.org/abs/gr-qc/9702026),
/// Eq. 28). Negative. (Their text's worked numbers, "6.2×10⁶⁵ g" and "a quarter of a solar mass",
/// disagree with their Eq. 28, which gives 6.9×10⁶² kg and 0.56 solar masses for the same
/// inputs.) Alcubierre's own tanh wall of steepness σ costs exactly `−(v²/36)R²σ·c²/G` for
/// `σR ≫ 1`, 2/3 of this with `Δ ≈ 2/σ`; Lobo & Visser's estimate `−v²R²σ` is 36 times that.
pub fn alcubierre_energy(radius: f64, wall: f64, speed: f64) -> f64 {
    -speed * speed / 12.0 * (radius * radius / wall + wall / 12.0) * C * C / G
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

/// Hawking temperature (K) of a Kerr black hole of mass `m` (kg) and spin `chi` = `cJ/(GM²)`
/// (−1 to 1), from its surface gravity `κ = c³√(1 − χ²)/(2GM(1 + √(1 − χ²)))`
/// ([Townsend](https://arxiv.org/abs/gr-qc/9707012), ch. 4; T = κ/2π in
/// [Wald](https://arxiv.org/abs/gr-qc/9912119)). Spin cools it by the factor `2√(1 − χ²)/(1 + √(1 − χ²))`:
/// 3/4 at χ = 0.8, to zero for an extremal hole.
pub fn kerr_hawking_temperature(m: f64, chi: f64) -> f64 {
    let s = (1.0 - chi * chi).sqrt();
    horizon_temperature(C.powi(3) * s / (2.0 * G * m * (1.0 + s)))
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
    /// ([Puthoff 2002](https://arxiv.org/abs/gr-qc/9909037)): the exponential geometry,
    /// **speculative**. The PV theory that produces it is excluded by binary-pulsar radiation
    /// ([`pv_radiated_power_ratio`]).
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

/// The orbital decay of a binary of eccentricity `e` in Puthoff's PV theory, as a fraction of
/// general relativity's: `½(51e⁴ + 396e² + 128)/(37e⁴ + 292e² + 96)`, exactly 2/3 for a circular
/// orbit ([Ibison 2003](https://arxiv.org/abs/astro-ph/0302273), Eq. 65). Its scalar waves carry
/// a third less power than GR's tensor waves, and the binary pulsars agree with GR, so the PV
/// theory is excluded; the exponential geometry it produces is a separate question.
///
/// ```
/// use gravitics::spacetime::pv_radiated_power_ratio;
///
/// // Hulse–Taylor, e = 0.6171340: the observed decay is 0.9983 ± 0.0016 of GR's (Weisberg &
/// // Huang 2016), while PV predicts 0.6732, about 200σ away.
/// let pv = pv_radiated_power_ratio(0.6171340);
/// assert!((pv - 0.67321).abs() < 1e-5 && (0.9983 - pv) / 0.0016 > 200.0);
/// ```
pub fn pv_radiated_power_ratio(e: f64) -> f64 {
    let e2 = e * e;
    0.5 * (51.0 * e2 * e2 + 396.0 * e2 + 128.0) / (37.0 * e2 * e2 + 292.0 * e2 + 96.0)
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
/// `D`. Spin and inclination make a Kerr shadow smaller: its mean, median and area-equivalent
/// radii by at most `2(1 + √2)/(3√3) − 1` = −7.08% (`a = 1` seen from the pole, see
/// [`kerr_pole_on_shadow`]), its horizontal width by up to 13.4%. The EHT quotes about 7.5%
/// ([EHT Sgr A* Paper VI](https://arxiv.org/abs/2311.09484)).
pub fn shadow_diameter(mass: f64, distance: f64) -> f64 {
    6.0 * 3f64.sqrt() * G * mass / (C * C * distance)
}

/// The fractional deviation `δ = d_sh/(6√3 θ_g) − 1` of a shadow of angular diameter `diameter`
/// from a Schwarzschild black hole of angular gravitational radius `theta_g = GM/(c²D)` (both in
/// the same angular unit), as the EHT defines it. For the mean radius a Kerr black hole gives `δ`
/// from exactly −0.0708 to 0; the EHT uses −0.075 as its published Kerr bound.
pub fn shadow_deviation(diameter: f64, theta_g: f64) -> f64 {
    diameter / (6.0 * 3f64.sqrt() * theta_g) - 1.0
}

/// The shadow radius, in units of `GM/c²`, that a Kerr black hole of spin `a` (−1 to 1) shows an
/// observer on its spin axis. Its edge is the spherical photon orbit with no angular momentum
/// about the axis, the largest root of `r³ − 3r² + a²r + a² = 0`, where Teo's constants
/// (Gen. Rel. Grav. 35, 1909, 2003) give the radius `√(η + a²) = √(8r³/(r² − 1))`. `√27` without
/// spin and `2(1 + √2)` at `a = 1`.
pub fn kerr_pole_on_shadow(a: f64) -> f64 {
    let r = polar_photon_orbit(a);
    (8.0 * r * r * r / (r * r - 1.0)).sqrt()
}

/// The spherical photon orbit with `ξ = 0`, the largest root of `r³ − 3r² + a²r + a² = 0`.
fn polar_photon_orbit(a: f64) -> f64 {
    // Newton from r = 3 descends monotonically: the cubic is rising and convex for r > 1.
    let mut r = 3.0;
    for _ in 0..100 {
        let step = (r * r * r - 3.0 * r * r + a * a * r + a * a) / (3.0 * r * r - 6.0 * r + a * a);
        r -= step;
        if step.abs() <= 1e-16 * r {
            break;
        }
    }
    r
}

/// Bardeen's image point `(α, β²)` of the spherical photon orbit at radius `r` (`GM/c²` = 1):
/// `ξ = L_z/E` and `η = Q/E²` give `α = −ξ/sin i` and `β² = η + a²cos²i − ξ²cot²i`.
fn bardeen(a: f64, inclination: f64, r: f64) -> (f64, f64) {
    let xi = (r * r * (3.0 - r) - a * a * (r + 1.0)) / (a * (r - 1.0));
    let eta = r.powi(3) * (4.0 * a * a - r * (r - 3.0).powi(2)) / (a * a * (r - 1.0).powi(2));
    let (s, c) = inclination.sin_cos();
    (-xi / s, eta + a * a * c * c - xi * xi * c * c / (s * s))
}

/// A point `[α, β]` (β ≥ 0, units of `GM/c²`) on the edge of the shadow of a Kerr black hole of
/// spin `a` (nonzero, |a| < 1) seen from `inclination` (rad) to its spin axis, cast by the
/// spherical photon orbit of Boyer–Lindquist radius `r`, or `None` where that orbit's light does
/// not reach the observer (`β² < 0`). Bardeen's curve (1973, as restated by
/// [Cunha & Herdeiro 2018](https://arxiv.org/abs/1801.00860) and
/// [Johannsen 2013](https://arxiv.org/abs/1501.02814)):
///
/// `ξ = [r²(3 − r) − a²(r + 1)]/[a(r − 1)]`, `η = r³[4a² − r(r − 3)²]/[a²(r − 1)²]`,
/// `α = −ξ/sin i`, `β = √(η + a²cos²i − ξ²cot²i)`.
///
/// The image is centred on the black hole, with `β` along the projected spin axis; for `a > 0`
/// the prograde side is at negative `α`. Edge-on, `r = 3` gives `(2a, 3√3)` at every spin, the
/// top of the shadow, so its vertical diameter is `6√3` ([Cunha & Herdeiro 2018]).
///
/// [Cunha & Herdeiro 2018]: https://arxiv.org/abs/1801.00860
pub fn kerr_shadow_edge(a: f64, inclination: f64, r: f64) -> Option<[f64; 2]> {
    let (alpha, beta2) = bardeen(a, inclination, r);
    (beta2 >= 0.0).then(|| [alpha, beta2.sqrt()])
}

/// The mean radius `⟨R⟩` and the displacement `D` (units of `GM/c²`) of the shadow of a Kerr
/// black hole of spin `a` (|a| < 1) seen from `inclination` (rad) to its spin axis, as
/// [Johannsen & Psaltis 2010](https://arxiv.org/abs/1005.1931) define them on
/// [`kerr_shadow_edge`]: `D = |α_max + α_min|/2`, and `⟨R⟩ = ∫R dθ/∫dθ` with
/// `R = √((α − D)² + β²)` measured from the displaced centre but the angle `tan θ = β/α` measured
/// about the black hole. (The angle about the centre gives 5.031316 instead of 5.027489 at
/// `a = 0.9`, `i = 90°`.) The mean diameter is `2⟨R⟩`.
///
/// The edge runs over `r` between the roots of `β² = 0`, bracketed by the polar orbit of
/// [`kerr_pole_on_shadow`] and the equatorial ones `2{1 + cos[(2/3)arccos(∓a)]}` (Bardeen, Press
/// & Teukolsky 1972), which are the ends edge-on. Samples clustered as `r₁ + (r₂ − r₁)(1 − cos t)/2`
/// follow β's square-root ends, and the trapezoid rule over the 10⁵ used converges to `10⁻¹⁰`
/// (`2×10⁻¹⁰` at `a = 0.998`).
/// Seen from the pole it is [`kerr_pole_on_shadow`] and `D = 0`.
///
/// ```
/// use gravitics::spacetime::kerr_shadow_size;
///
/// // a = 0.9 edge-on: the mean diameter against Johannsen 2013's fit L_K, accurate to 0.09%.
/// let (radius, d) = kerr_shadow_size(0.9, std::f64::consts::FRAC_PI_2);
/// assert!((2.0 * radius - 10.054979).abs() < 1e-6 && (2.0 * radius / 10.056937 - 1.0).abs() < 9e-4);
/// assert!((d - 1.993949).abs() < 1e-6);
/// ```
pub fn kerr_shadow_size(a: f64, inclination: f64) -> (f64, f64) {
    let a = a.abs();
    if a == 0.0 {
        return (27f64.sqrt(), 0.0);
    }
    if inclination.sin() == 0.0 {
        return (kerr_pole_on_shadow(a), 0.0);
    }
    let polar = polar_photon_orbit(a);
    // Bisect from the polar orbit, where β² > 0 at every inclination, toward an equatorial one.
    let end = |sign: f64| {
        let (mut inside, mut outside) = (polar, 2.0 * (1.0 + (2.0 / 3.0 * (sign * a).acos()).cos()));
        for _ in 0..100 {
            let mid = 0.5 * (inside + outside);
            if bardeen(a, inclination, mid).1 > 0.0 { inside = mid } else { outside = mid }
        }
        inside
    };
    let (r1, r2) = (end(-1.0), end(1.0));
    let n = 100_000;
    let edge: Vec<[f64; 2]> = (0..=n)
        .map(|k| {
            let (alpha, beta2) = bardeen(a, inclination, r1 + (r2 - r1) * (1.0 - (PI * k as f64 / n as f64).cos()) / 2.0);
            // β is zero at the ends by definition; rounding r there leaves β ~ 10⁻⁷, a sliver of angle.
            [alpha, if k == 0 || k == n { 0.0 } else { beta2.max(0.0).sqrt() }]
        })
        .collect();
    let (lo, hi) = edge.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), p| (lo.min(p[0]), hi.max(p[0])));
    let d = 0.5 * (lo + hi);
    let radius = |p: [f64; 2]| (p[0] - d).hypot(p[1]);
    let angle = |p: [f64; 2]| p[1].atan2(p[0]);
    let sum: f64 = edge.windows(2).map(|w| 0.5 * (radius(w[0]) + radius(w[1])) * (angle(w[1]) - angle(w[0]))).sum();
    // The upper half of the edge runs from θ = π to 0.
    (-sum / PI, d)
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
    fn a_kerr_shadow_at_any_inclination() {
        let edge_on = core::f64::consts::FRAC_PI_2;
        // Edge-on, r = 3 gives ξ = −2a and η = 27 at every spin, the top of the shadow at
        // (2a, 3√3) where η is stationary: the vertical diameter is always 6√3.
        for a in [0.1, 0.5, 0.9, 0.998, -0.7] {
            let [alpha, beta] = kerr_shadow_edge(a, edge_on, 3.0).unwrap();
            assert!((alpha - 2.0 * a).abs() < 1e-14 && close(beta, 27f64.sqrt(), 1e-15), "a = {a}");
            assert!([2.95, 2.999, 3.001, 3.05].iter().all(|&r| kerr_shadow_edge(a, edge_on, r).is_none_or(|p| p[1] < beta)));
        }
        // Its ends edge-on are the equatorial photon orbits: at a = 0.9, −2.844421403 and 6.832319230.
        let orbit = |sign: f64| 2.0 * (1.0 + (2.0 / 3.0 * (sign * 0.9f64).acos()).cos());
        assert!((bardeen(0.9, edge_on, orbit(-1.0)).0 + 2.844_421_403).abs() < 1e-9);
        assert!((bardeen(0.9, edge_on, orbit(1.0)).0 - 6.832_319_230).abs() < 1e-9);
        // Seen from the pole it is kerr_pole_on_shadow, 4.91606338 at a = 0.9, and its limit.
        assert!((kerr_pole_on_shadow(0.9) - 4.916_063_38).abs() < 1e-8);
        for a in [0.3, 0.9, 0.999] {
            assert_eq!(kerr_shadow_size(a, 0.0), (kerr_pole_on_shadow(a), 0.0));
            let (radius, d) = kerr_shadow_size(a, 1e-5);
            assert!((radius - kerr_pole_on_shadow(a)).abs() < 1e-9 && d.abs() < 1e-4, "a = {a}: {radius}");
        }
        // ⟨R⟩ and D against the research's tables, recomputed here with an independent quadrature
        // over θ (a root in r for each θ): in units of GM/c².
        for (a, degrees, radius, d) in [
            (0.5, 90.0, 5.157_024_802_42, 1.020_944_533_00),
            (0.9, 17.0, 4.924_470_288_51, 0.600_227_062_12),
            (0.9, 60.0, 4.997_506_316_25, 1.746_750_502_48),
            (0.9, 90.0, 5.027_489_401_78, 1.993_948_913_49),
            (0.998, 90.0, 4.925_689_542_49, 2.442_889_238_42),
        ] {
            let (r, dd) = kerr_shadow_size(a, f64::to_radians(degrees));
            assert!((r - radius).abs() < 1e-9 && (dd - d).abs() < 1e-9, "a = {a}, i = {degrees}: {r} {dd}");
            // Inside the Kerr range, 2(1 + √2) to √27; spin's sign only mirrors it.
            assert!(r > 2.0 * (1.0 + 2f64.sqrt()) && r < 27f64.sqrt() && kerr_shadow_size(-a, f64::to_radians(degrees)) == (r, dd));
        }
        assert_eq!(kerr_shadow_size(0.0, 1.0), (27f64.sqrt(), 0.0));
    }

    #[test]
    fn a_spinning_black_hole_is_colder() {
        let sun = metrology::SOLAR_MASS;
        let cold = |chi: f64| kerr_hawking_temperature(sun, chi);
        assert!(close(cold(0.0), hawking_temperature(sun), 1e-15));
        // κ(χ)/κ(0) = 2√(1 − χ²)/(1 + √(1 − χ²)): exactly 8/9, 3/4 and 0 (Townsend, ch. 4).
        assert!(close(cold(0.6) / cold(0.0), 8.0 / 9.0, 1e-15) && close(cold(0.8) / cold(0.0), 0.75, 1e-15) && cold(1.0) == 0.0);
        // With SOLAR_MASS = GM☉/G: 61.70260 nK, and 46.27695 nK at χ = 0.8.
        assert!(close(cold(0.0), 6.170_260e-8, 1e-6) && close(cold(0.8), 4.627_695e-8, 1e-6));
        // The first law at fixed J, in G = c = M = 1: dM/dA = κ/8π with the horizon area
        // A = 8π(M² + √(M⁴ − J²)), 0.0082883527 at χ = 0.7.
        let kappa = 2.0 * PI * K_B * cold(0.7) / HBAR * G * sun / C.powi(3);
        let area = |m: f64| 8.0 * PI * (m * m + (m.powi(4) - 0.49).sqrt());
        let h = 1e-6;
        assert!(close(2.0 * h / (area(1.0 + h) - area(1.0 - h)), kappa / (8.0 * PI), 1e-9));
        assert!((kappa / (8.0 * PI) - 0.008_288_352_7).abs() < 1e-10);
    }

    #[test]
    fn exact_deflection_against_darwin_bozza_and_the_series() {
        let m = metrology::SOLAR_MASS;
        let rg = G * m / (C * C);
        let at = |b: f64| exact_deflection(m, b * rg, Metric::GeneralRelativity).unwrap();
        // Darwin's closed form, as Iyer & Petters (2007) print it, from the closest approach r₀:
        // 4√(r₀/Q)[K(k) − F(Ψ, k)] − π, Q = √((r₀ − 2)(r₀ + 6)), k² = (Q − r₀ + 6)/(2Q),
        // sin²Ψ = (Q + 2 − r₀)/(Q + 6 − r₀).
        // (Closer to the photon sphere α is ill-conditioned in b: the round trip r₀ → b → r₀
        // alone costs 10⁻¹⁶/(b/b_c − 1). Bozza's limit below covers that regime.)
        for r0 in [3.05f64, 3.5, 4.0, 6.0, 10.0, 100.0] {
            let q = ((r0 - 2.0) * (r0 + 6.0)).sqrt();
            let k2 = (q - r0 + 6.0) / (2.0 * q);
            let s2 = (q + 2.0 - r0) / (q + 6.0 - r0);
            let (big_k, _) = geopotential::special::complete_elliptic(k2);
            let f = geopotential::special::elliptic_f(s2.sqrt().asin(), k2);
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
        // Against an independent Python quadrature (the frontier-gaps research's PV notes): at
        // b = 6M, 10M and b_c(1 + 10⁻³), both metrics.
        let at = |b: f64, metric| exact_deflection(m, b * rg, metric).unwrap();
        for (b, gr, pv) in [(6.0, 1.719388310, 1.999746113), (10.0, 0.590395788, 0.610437437)] {
            assert!((at(b, Metric::GeneralRelativity) - gr).abs() < 1e-9 && (at(b, Metric::PolarizableVacuum) - pv).abs() < 1e-9, "b = {b}");
        }
        assert!((at(27f64.sqrt() * 1.001, Metric::GeneralRelativity) - 6.5106448).abs() < 1e-7);
        assert!((at(critical * 1.001, Metric::PolarizableVacuum) - 6.4454619).abs() < 1e-7);
        // Bozza's strong-field offset for PV, b̄ = −0.46581756 (two methods in the notes agree to
        // 2×10⁻⁸), against GR's −0.40023: α + ln ε at ε = 10⁻⁶.
        let eps = 1e-6;
        assert!((at(critical * (1.0 + eps), Metric::PolarizableVacuum) + eps.ln() + 0.46581756).abs() < 1e-5);
    }

    #[test]
    fn jupiters_oblateness_in_the_deflection() {
        let (gm, radius, j2) = (DUR2020_GM, DUR2020_RADIUS, -DUR2020_C[0]);
        let field = jupiter(2);
        let mas = PI / 180.0 / 3600.0 * 1e-3;
        for (b, lon) in [(radius, 0.0), (2.0 * radius, 0.7), (radius * 1.1, -2.0)] {
            // An equatorial ray, at any longitude: (4GM/c²b)(1 + J₂R²/b²).
            let point = [b * f64::cos(lon), b * f64::sin(lon), 0.0];
            let dir = [-f64::sin(lon), f64::cos(lon), 0.0];
            let alpha = multipole_deflection(&field, point, dir, 1.0);
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

    // Jupiter from Juno: Durante et al. 2020's PDS product (jnogrv_0002, SHADR
    // GRV_JUGR_SHA_JUP_U_DUR2020_20200203_V01): GM, the reference radius, and the unnormalised
    // zonal C_n0 = −J_n for n = 2 to 12. Its J₂ = 14696.5735×10⁻⁶ ± 0.0017 (3σ) includes a
    // 6.72×10⁻⁸ permanent tide (14696.5063 without it); Iess et al. 2018 give 14696.572 ± 0.014.
    const DUR2020_GM: f64 = 1.266_865_341e17;
    const DUR2020_RADIUS: f64 = 7.1492e7;
    const DUR2020_C: [f64; 11] = [-1.469_657_35e-2, 4.5e-8, 5.866_085e-4, 7.23e-8, -3.420_07e-5, -1.2e-7, 2.422e-6, 1.13e-7, -1.81e-7, -1.6e-8, -6.2e-8];

    /// The DUR2020 zonal field to degree `nmax`, normalised as `C̄_n0 = C_n0/√(2n + 1)`.
    fn jupiter(nmax: usize) -> geopotential::Coefficients {
        let mut field = geopotential::Coefficients::zeros(DUR2020_GM, DUR2020_RADIUS, nmax);
        field.set(0, 0, 1.0, 0.0);
        for n in 2..=nmax {
            field.set(n, 0, DUR2020_C[n - 2] / ((2 * n + 1) as f64).sqrt(), 0.0);
        }
        field
    }

    #[test]
    fn jupiters_quadrupole_at_any_orientation() {
        // Crosta & Mignard 2006 (astro-ph/0512359) Eqs. defqu1–defqu2: with n from the ray toward
        // the planet, t along it and m = t × n, α = Φ₁n + Φ₂m where
        // Φ₁ = (4GM/c²b)[1 + J₂(R/b)²(1 − 2(n·z)² − (t·z)²)] and Φ₂ = (8GM/c²b)J₂(R/b)²(m·z)(n·z).
        let field = jupiter(2);
        let (gm, radius, j2) = (DUR2020_GM, DUR2020_RADIUS, -DUR2020_C[0]);
        let uas = PI / 180.0 / 3600.0 * 1e-6;
        let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        let (s30, c30, s45) = (0.5, 0.75f64.sqrt(), 0.5f64.sqrt());
        // (point, direction, |α| and its quadrupole part in µas, from the research's Born
        // integral, which matched Crosta's forms to 2×10⁻¹²): over the pole, the quadrupole
        // subtracts; along the spin axis it vanishes; a ray at 30° to the equator, closest at 45°
        // latitude, has sin²ϑ = 0.75 to the axis, so 4GMJ₂R²sin²ϑ/(c²b³).
        for (point, dir, size, quad) in [
            ([0.0, 0.0, radius], [0.0, 1.0, 0.0], 16_028.272_6, 239.074_26),
            ([radius, 0.0, 0.0], [0.0, 0.0, 1.0], 16_267.346_9, 0.0),
            ([radius * s45, 0.0, radius * s45], [0.0, c30, s30], 17_423.170_1, 219.069_54),
        ] {
            let along = dot(point, dir);
            let closest: [f64; 3] = core::array::from_fn(|i| point[i] - along * dir[i]);
            let b = dot(closest, closest).sqrt();
            let n = closest.map(|x| -x / b);
            let m = [dir[1] * n[2] - dir[2] * n[1], dir[2] * n[0] - dir[0] * n[2], dir[0] * n[1] - dir[1] * n[0]];
            let mono = 4.0 * gm / (C * C * b);
            let q = j2 * (radius / b).powi(2);
            let phi1 = mono * (1.0 + q * (1.0 - 2.0 * n[2] * n[2] - dir[2] * dir[2]));
            let phi2 = 2.0 * mono * q * m[2] * n[2];
            let alpha = multipole_deflection(&field, point, dir, 1.0);
            let crosta: [f64; 3] = core::array::from_fn(|i| phi1 * n[i] + phi2 * m[i]);
            let gap: [f64; 3] = core::array::from_fn(|i| alpha[i] - crosta[i]);
            assert!(dot(gap, gap).sqrt() < 1e-11 * mono * q, "{dir:?}");
            assert!((dot(alpha, alpha).sqrt() / uas - size).abs() < 1e-4, "{dir:?}");
            let rest: [f64; 3] = core::array::from_fn(|i| alpha[i] - mono * n[i]);
            assert!((dot(rest, rest).sqrt() / uas - quad).abs() < 1e-5, "{dir:?}");
        }
        // Over the pole and along the axis these are 4GM/(c²R)(1 − J₂) and 4GM/(c²R) exactly.
        let limb = 4.0 * gm / (C * C * radius);
        let over = multipole_deflection(&field, [0.0, 0.0, radius], [0.0, 1.0, 0.0], 1.0);
        assert!(close(-over[2], limb * (1.0 - j2), 1e-12) && over[0].abs() + over[1].abs() < 1e-15 * limb);
        let axial = multipole_deflection(&field, [radius, 0.0, 0.0], [0.0, 0.0, 1.0], 1.0);
        assert!(close(-axial[0], limb, 1e-12) && axial[1].abs() + axial[2].abs() < 1e-15 * limb);
    }

    #[test]
    fn jupiters_full_zonal_field_on_an_equatorial_ray() {
        // For an equatorial ray the zonal field gives, exactly (the research's series, from
        // P_n(0) and ∫(b² + l²)^(−(n+1)/2) dl), α = (4GM/c²b)[1 + Σ(−1)^(k+1) J_2k (R/b)^2k] in the
        // plane and (4GM/c²b) Σ(−1)^(k+1) J_(2k+1) (R/b)^(2k+1) north: 16516.561365 and +0.004495 µas
        // at the limb, 8163.860339 and −0.000018 at 2R; J₄ alone is Crosta & Mignard's "10 µas".
        let field = jupiter(12);
        let uas = PI / 180.0 / 3600.0 * 1e-6;
        let j = |n: usize| -DUR2020_C[n - 2];
        for (scale, plane, north) in [(1.0, 16_516.561_365, 0.004_495), (2.0, 8_163.860_339, -0.000_018)] {
            let b = scale * DUR2020_RADIUS;
            let x = DUR2020_RADIUS / b;
            let mono = 4.0 * DUR2020_GM / (C * C * b);
            let even = 1.0 + (1..=6).map(|k| if k % 2 == 1 { 1.0 } else { -1.0 } * j(2 * k) * x.powi(2 * k as i32)).sum::<f64>();
            let odd = (1..=5).map(|k| if k % 2 == 1 { 1.0 } else { -1.0 } * j(2 * k + 1) * x.powi(2 * k as i32 + 1)).sum::<f64>();
            let alpha = multipole_deflection(&field, [b, 0.0, 0.0], [0.0, 1.0, 0.0], 1.0);
            assert!(close(-alpha[0], mono * even, 1e-12) && (alpha[2] - mono * odd).abs() < 1e-12 * mono, "b = {scale}R");
            assert!((-alpha[0] / uas - plane).abs() < 1e-6 && (alpha[2] / uas - north).abs() < 1e-6, "b = {scale}R");
        }
        assert!((DUR2020_GM / (C * C * DUR2020_RADIUS) * 4.0 * j(4) / uas + 9.5426).abs() < 1e-4);
    }
}
