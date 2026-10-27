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

/// The slope `f′(r)` of Alcubierre's shape function `f = [tanh σ(r + R) − tanh σ(r − R)]/(2 tanh σR)`.
fn alcubierre_slope(r: f64, radius: f64, sigma: f64) -> f64 {
    let sech2 = |x: f64| 1.0 / x.cosh().powi(2);
    sigma * (sech2(sigma * (r + radius)) - sech2(sigma * (r - radius))) / (2.0 * (sigma * radius).tanh())
}

/// The energy density (J/m³) that Eulerian observers see at `point` (m from the centre, x along
/// the motion) of Alcubierre's warp bubble of `radius` (m) and wall steepness `sigma` (m⁻¹) moving
/// at `speed` (in units of c): `−(c⁴/32πG)v²f′²(y² + z²)/r²`, with his tanh shape function `f`
/// ([Alcubierre 1994](https://arxiv.org/abs/gr-qc/0009013)). It is negative everywhere and
/// deepest where the wall crosses the plane through the centre across the motion, `−(c⁴/32πG)v²σ²/4`
/// at `r = R` for `σR ≫ 1`. [`stress_energy`] reproduces it from the metric.
///
/// Warp Factory's test case (v = 0.1, R = 300 m, σ = 0.015 m⁻¹) gives −6.775×10³⁵ J/m³ there,
/// 4.0 times less than the colour bar of its figure, about −2.7×10³⁶
/// ([Helmerich et al. 2024](https://arxiv.org/abs/2404.03095)), as if that was drawn for v = 0.2
/// or σ = 0.03.
///
/// White et al. 2021 ([EPJ C 81, 677](https://doi.org/10.1140/epjc/s10052-021-09484-z)) compared
/// the *shape* of a Casimir cavity's negative energy with a 1 µm bubble, with no correlation
/// coefficient and no data (**speculative**). Since the density goes as `v²f′²`, a shape carries no
/// speed; matching magnitudes instead, the cavity buys 1.7×10⁻²² m/s:
///
/// ```
/// use gravitics::spacetime::{alcubierre_energy_density, casimir_energy_density};
/// use metrology::C;
///
/// // A bubble 1 µm across (R = 0.5 µm, a 0.5 µm wall: σ = 4×10⁶ m⁻¹) at light speed, across its
/// // side wall: −5.2×10⁵⁴ J/m³ at the deepest.
/// let peak = (1..=2000).map(|k| alcubierre_energy_density(0.5e-6, 4e6, 1.0, [0.0, k as f64 * 5e-10, 0.0])).fold(0.0, f64::min);
/// assert!((peak / -5.1677e54 - 1.0).abs() < 1e-4);
/// // The 4 µm cavity's −1.7×10⁻⁶ J/m³ sets v² = u/peak.
/// let speed = (casimir_energy_density(4e-6) / peak).sqrt() * C;
/// assert!((speed / 1.7159e-22 - 1.0).abs() < 1e-4);
/// ```
pub fn alcubierre_energy_density(radius: f64, sigma: f64, speed: f64, point: [f64; 3]) -> f64 {
    let [x, y, z] = point;
    let r2 = x * x + y * y + z * z;
    if r2 == 0.0 {
        return 0.0;
    }
    -C.powi(4) / (32.0 * PI * G) * speed * speed * alcubierre_slope(r2.sqrt(), radius, sigma).powi(2) * (y * y + z * z) / r2
}

/// The total energy (as a mass, kg) of Alcubierre's bubble of `radius` (m) with his own tanh wall
/// of steepness `sigma` (m⁻¹), moving at `speed` (in units of c): the volume integral of
/// [`alcubierre_energy_density`], `−(v²/12)(c²/G)∫₀^∞ r²f′² dr`, by Gauss–Legendre panels a wall
/// wide. Negative. For `σR ≫ 1` it is `−(v²/36)R²σ·c²/G`, 2/3 of [`alcubierre_energy`] under
/// Pfenning & Ford's matching `Δ ≈ 2/σ`. Over all walls it is least, `0.6675R(v²/12)c²/G`, at
/// `σR = 1.229`, which gives up the flat interior; [`alcubierre_minimum_energy`] keeps one.
///
/// White's "Warp Field Mechanics 101" (2011, [NTRS 20110015936](https://ntrs.nasa.gov/citations/20110015936))
/// proposed seeing a "boost" of 1.0000001 over 1 cm in an interferometer. Computed both ways:
/// - in his own model, `γ = cosh(½ln|1 − v²|)` = 1.0000001 needs v = 0.0299c, and a 0.5 cm
///   bubble at that speed (wall ≈ radius) needs −3.4×10²⁰ kg;
/// - in general relativity, a generous capacitor (10⁶ V/m filling a 0.1 m sphere) perturbs the
///   metric by `2GM/(c²r)` = 3×10⁻⁴⁵, 38 orders below the target (**fringe** as a test). Lee &
///   Cleaver 2014 ([arXiv:1407.7772](https://arxiv.org/abs/1407.7772)) reach the same verdict,
///   though their 1.92×10⁻²³ m Schwarzschild radius is 10²¹–10²³ too large.
///
/// ```
/// use gravitics::spacetime::alcubierre_tanh_energy;
/// use metrology::{C, EPS_0, G};
///
/// // γ = (u + 1/u)/2 with u = √(1 − v²).
/// let boost: f64 = 1.0 + 1e-7;
/// let u = boost - (boost * boost - 1.0).sqrt();
/// let v = (1.0 - u * u).sqrt();
/// assert!((v - 0.029900).abs() < 1e-6);
/// assert!((alcubierre_tanh_energy(0.005, 200.0, v) / -3.4340e20 - 1.0).abs() < 1e-4);
/// // The capacitor: ε₀E²/2 = 4.4 J/m³ in a 0.1 m sphere.
/// let mass = EPS_0 * 1e12 / 2.0 * 4.0 / 3.0 * std::f64::consts::PI * 0.1f64.powi(3) / (C * C);
/// let perturbation = 2.0 * G * mass / (C * C * 0.1);
/// assert!((perturbation / 3.0645e-45 - 1.0).abs() < 1e-4 && (1e-7 / perturbation).log10() > 37.5);
/// ```
pub fn alcubierre_tanh_energy(radius: f64, sigma: f64, speed: f64) -> f64 {
    // f′ is below e⁻⁸⁰ of its peak beyond 40 wall widths 1/σ of the wall.
    let (lo, hi) = ((radius - 40.0 / sigma).max(0.0), radius + 40.0 / sigma);
    let panels = ((hi - lo) * sigma).ceil();
    let half = 0.5 * (hi - lo) / panels;
    let (x, w) = geopotential::gauss_legendre(16);
    let integral: f64 = (0..panels as usize)
        .map(|k| {
            let mid = lo + (2 * k + 1) as f64 * half;
            half * x.iter().zip(&w).map(|(x, w)| w * (mid + half * x).powi(2) * alcubierre_slope(mid + half * x, radius, sigma).powi(2)).sum::<f64>()
        })
        .sum();
    -speed * speed / 12.0 * integral * C * C / G
}

/// The least energy (as a mass, kg) that any bubble of Alcubierre's form moving at `speed` (in
/// units of c) needs to keep a flat interior of radius `r0` (m): `−(v²/12)r₀c²/G`. The energy is
/// `−(v²/12)(c²/G)∫r²f′² dr` for any shape `f`, which with `f = 1` inside `r₀` and 0 far away is a
/// sphere's electrostatic capacitance, least for `f = min(r₀/r, 1)`
/// ([Bobrick & Martire 2021](https://arxiv.org/abs/2102.06824), appendix). A 100 m interior at
/// light speed needs 1.122×10²⁸ kg, 5.9 Jupiter masses, and Alcubierre's own σR = 8 profile 2.68
/// times that, their "factor of three".
///
/// White's 2013 chart of the exotic mass a 10 m bubble needs at 10c
/// ([NTRS 20130011213](https://ntrs.nasa.gov/citations/20130011213)) cannot be reproduced: over all
/// thicknesses of Alcubierre's wall the energy has a floor of 3.7×10²⁸ kg, and the chart's thick
/// end, about 10²⁰ kg, sits 10^8.6 below it (**fringe**).
///
/// ```
/// use gravitics::spacetime::{alcubierre_minimum_energy, alcubierre_tanh_energy};
///
/// // R = 5 m at v = 10: σR from 0.5 to 3 finds the floor at σR ≈ 1.23.
/// let floor = (50..=300).map(|k| alcubierre_tanh_energy(5.0, k as f64 / 500.0, 10.0)).fold(f64::NEG_INFINITY, f64::max);
/// assert!((floor / -3.7451e28 - 1.0).abs() < 1e-4 && (floor / -1e20).log10() > 8.5);
/// // A flat 5 m interior costs more, 5.6×10²⁸ kg.
/// assert!((alcubierre_minimum_energy(5.0, 10.0) / -5.6108e28 - 1.0).abs() < 1e-4);
/// ```
pub fn alcubierre_minimum_energy(r0: f64, speed: f64) -> f64 {
    -speed * speed / 12.0 * r0 * C * C / G
}

/// `(tr K)²` and `tr K²` for the symmetric part `K` of a shift's Jacobian.
fn shift_invariants(jacobian: [[f64; 3]; 3]) -> (f64, f64) {
    let trace = jacobian[0][0] + jacobian[1][1] + jacobian[2][2];
    let square = (0..9).map(|n| (0.5 * (jacobian[n / 3][n % 3] + jacobian[n % 3][n / 3])).powi(2)).sum();
    (trace * trace, square)
}

/// The Eulerian energy density (J/m³) of a warp drive of Natário's class, flat slices and unit
/// lapse, `ds² = −c²dt² + |dx − v c dt|²`, from the Jacobian `jacobian[i][j] = ∂v_i/∂x_j` (m⁻¹,
/// with `v` in units of c) of its shift: `(c⁴/16πG)(K² − tr K²)`, `K` its symmetric part
/// ([Santiago, Schuster & Visser 2022](https://arxiv.org/abs/2105.03079)). Its volume integral is
/// `−(c⁴/32πG)∫|∇×v|² d³x ≤ 0`, so a drive without vorticity that is positive somewhere is negative
/// elsewhere, and every superluminal drive of this class violates the null energy condition
/// ([`shift_drive_null_sum`]).
///
/// The "positive-energy" drives rest on the Eulerian observer alone. Lentz 2021's positivity is
/// disputed as a sign error ([Celmaster & Rubin 2025](https://arxiv.org/abs/2511.18251), preprint):
/// **excluded** as to positivity. Fell & Heisenberg 2021 concede violations for other observers,
/// and their 9.25×10⁴³ J is 5.2×10⁻⁴ of the Sun's rest energy, 3.3 orders below it rather than
/// their "four" (**speculative**).
pub fn shift_drive_density(jacobian: [[f64; 3]; 3]) -> f64 {
    let (k2, square) = shift_invariants(jacobian);
    C.powi(4) / (16.0 * PI * G) * (k2 - square)
}

/// `ρ + p̄` (J/m³), the Eulerian density plus the mean principal pressure, of the drives of
/// [`shift_drive_density`]: `(c⁴/24πG)(−2ℒₙK + K² − 3 tr K²)`, with `K` the symmetric part of
/// `jacobian` and `lie_trace` = `ℒₙK = (∂/∂(ct) + v·∇)(∇·v)` (m⁻²)
/// ([Santiago, Schuster & Visser 2022](https://arxiv.org/abs/2105.03079)). It is the mean of
/// `T(k, k)` over the null directions `k = n + e`, so where it is negative the null energy condition
/// fails. The sign of `K` matters here: `K_ij = +∂₍ᵢv_ⱼ₎` reproduces the Einstein tensor and the
/// opposite sign does not. At Alcubierre's side wall (v = 1, σR = 8, G = c = 1) it is −0.37136.
pub fn shift_drive_null_sum(jacobian: [[f64; 3]; 3], lie_trace: f64) -> f64 {
    let (k2, square) = shift_invariants(jacobian);
    C.powi(4) / (24.0 * PI * G) * (-2.0 * lie_trace + k2 - 3.0 * square)
}

/// The surface energy density σ (J/m², negative) and surface pressure `p` (N/m, positive, a
/// surface tension of −p) of Visser's thin-shell wormhole, two Schwarzschild exteriors of `mass`
/// (kg) joined at `radius` (m): `σ = −(c⁴/2πGa)√(1 − 2GM/ac²)` and
/// `p = (c⁴/4πGa)(1 − GM/ac²)/√(1 − 2GM/ac²)` ([Visser 1989](https://arxiv.org/abs/0809.0927);
/// [Poisson & Visser 1995](https://arxiv.org/abs/gr-qc/9506083)). `σ + p < 0`, so the throat
/// violates the null energy condition. A 1 m throat holds −2.7×10²⁷ kg, the same order as the
/// twelve edges of Visser's 1 m cube (−2.0×10²⁷ kg, [`wormhole_edge_tension`]).
pub fn thin_shell_wormhole(radius: f64, mass: f64) -> (f64, f64) {
    let x = G * mass / (radius * C * C);
    // c⁴/(2πGa) is 4a times the Morris–Thorne throat tension at a.
    let scale = 4.0 * radius * wormhole_throat_tension(radius);
    (-scale * (1.0 - 2.0 * x).sqrt(), 0.5 * scale * (1.0 - x) / (1.0 - 2.0 * x).sqrt())
}

/// The largest throat (m) of a Morris–Thorne wormhole (`b = r₀²/r`) that quantum fields can hold
/// open, when the quantum inequality samples a time `f r₀/c`, a fraction `f` of the throat:
/// `√(3/4π) l_P/f²`, 7.9×10⁻³² m at f = 0.01 ([Ford & Roman 1996](https://arxiv.org/abs/gr-qc/9510071)
/// §IV, whose `l_P/(2f²)` rounds 3/4π). The throat's density, −[`wormhole_throat_tension`], may not
/// fall below [`quantum_inequality`].
pub fn wormhole_throat_bound(f: f64) -> f64 {
    // −τ(r₀) = −τ(1 m)/r₀² against QI(f r₀/c) = QI(f/c)/r₀⁴, r₀ in metres.
    (quantum_inequality(f / C) / -wormhole_throat_tension(1.0)).sqrt()
}

/// The thickest band (m) of exotic matter that quantum fields allow at the throat of radius `r0`
/// (m) of Ford & Roman's "absurdly benign" wormhole, `b = r₀[1 − (r − r₀)/a₀]²` out to `r₀ + a₀`,
/// when the quantum inequality samples a fraction `f` of the band: `(3r₀l_P²/(8πf⁴))^{1/3}`,
/// 1.46×10⁻²¹ m for a 1 m throat at f = 0.01 ([Ford & Roman 1996](https://arxiv.org/abs/gr-qc/9510071);
/// their "0.2 fermi" for a 1 light-year throat is 0.31 fm).
pub fn wormhole_band_bound(r0: f64, f: f64) -> f64 {
    // The density at the throat, b′/(8πr₀²), is 2r₀/a₀ times −τ(r₀) = −τ(1 m)/r₀².
    (quantum_inequality(f / C) * r0 / (-2.0 * wormhole_throat_tension(1.0))).cbrt()
}

/// The most negative energy density (J/m³) a massless scalar field allows when averaged with a
/// Gaussian weight `e^{−t²/τ²}/(√π τ)` of width `tau` (s): `−3ħc/(64π²(cτ)⁴)`, half of
/// [`quantum_inequality`]'s Lorentzian bound. It is Fewster & Eveson's bound
/// `−(ħ/16π²c³)∫((√w)″)² dt` for any weight `w` ([Fewster & Eveson 1998](https://arxiv.org/abs/gr-qc/9805024)),
/// which for the Lorentzian weight gives `27/(2048π²)`, 9/64 of Ford & Roman's constant.
pub fn quantum_inequality_gaussian(tau: f64) -> f64 {
    0.5 * quantum_inequality(tau)
}

/// `x` moved by `dx` along coordinate `a`.
fn nudge(x: [f64; 4], a: usize, dx: f64) -> [f64; 4] {
    let mut y = x;
    y[a] += dx;
    y
}

/// The fourth-order central difference from the values at −2h, −h, h and 2h.
fn central(v: [f64; 4], h: f64) -> f64 {
    (8.0 * (v[2] - v[1]) - (v[3] - v[0])) / (12.0 * h)
}

/// The inverse of a symmetric 4×4 matrix, a column per right-hand side.
fn inverse(g: [[f64; 4]; 4]) -> [[f64; 4]; 4] {
    core::array::from_fn(|j| {
        let rows = (0..4).map(|i| [&g[i][..], &[if i == j { 1.0 } else { 0.0 }]].concat()).collect();
        let column = geopotential::linalg::gauss(rows);
        core::array::from_fn(|i| column[i])
    })
}

/// `Γ^λ_μν = ½g^λσ(∂_μg_σν + ∂_νg_σμ − ∂_σg_μν)`, indexed `[λ][μ][ν]`.
fn christoffel(metric: &impl Fn([f64; 4]) -> [[f64; 4]; 4], x: [f64; 4], h: f64) -> [[[f64; 4]; 4]; 4] {
    let inv = inverse(metric(x));
    let dg: [[[f64; 4]; 4]; 4] = core::array::from_fn(|a| {
        let v = [-2.0, -1.0, 1.0, 2.0].map(|k| metric(nudge(x, a, k * h)));
        core::array::from_fn(|m| core::array::from_fn(|n| central(v.map(|g| g[m][n]), h)))
    });
    core::array::from_fn(|l| core::array::from_fn(|m| core::array::from_fn(|n| 0.5 * (0..4).map(|s| inv[l][s] * (dg[m][s][n] + dg[n][s][m] - dg[s][m][n])).sum::<f64>())))
}

/// The stress-energy `T_μν` (J/m³, coordinate components) that the spacetime `metric` needs at the
/// event `x`: Einstein's equation `T = G/κ` (`κ = 8πG/c⁴`), with the Einstein tensor `G_μν` from
/// nested fourth-order central differences of step `step` (m), the recipe of Warp Factory
/// ([Helmerich et al. 2024](https://arxiv.org/abs/2404.03095)). Coordinates are in metres with
/// `x[0] = ct`, and `metric` gives the dimensionless `g_μν` (signature −+++). [`energy_conditions`]
/// reads it.
///
/// The error relative to the curvature is the truncation `(h/L)⁴`, with `L` the metric's shortest
/// length scale, plus the rounding `10⁻¹⁶(L/h)²`, so `h ≈ 10⁻³L` leaves 10⁻¹⁰ to 10⁻⁹. A vacuum
/// shows that residual: 3×10⁻⁹ of the curvature for Schwarzschild at `r = 3GM/c²` with
/// `h = 10⁻³GM/c²`. It reproduces Alcubierre's density
/// ([`alcubierre_energy_density`]) and Santiago, Schuster & Visser's `ρ + p̄`
/// ([`shift_drive_null_sum`]) from the metric alone.
///
/// Everett & Roman's Krasnikov tube ([gr-qc/9702049](https://arxiv.org/abs/gr-qc/9702049)),
/// `ds² = −dt² + (1 − k)dx dt + k dx² + dρ² + ρ²dφ²` (**speculative**) needs negative energy in
/// its wall: their Eq. T for the static observer gives `T_tt ε² = −0.04222 c⁴/G` at n = 10,
/// δ = 0.01, and their fit −0.04217, for the wall `k = 1 − (1 − δ/2)θ_ε(nε − ρ)` with
/// `θ_ε(ξ) = ½{tanh[2(2ξ/ε − 1)] + 1}`:
///
/// ```
/// use gravitics::spacetime::stress_energy;
/// use metrology::{C, G};
///
/// // ε = 1 m; the inner edge of the wall, ρ = (n − 1)ε.
/// let (n, delta) = (10.0, 0.01);
/// let k = |rho: f64| 1.0 - (1.0 - delta / 2.0) * 0.5 * ((2.0 * (2.0 * (n - rho) - 1.0)).tanh() + 1.0);
/// let tube = |x: [f64; 4]| {
///     let k = k(x[2].hypot(x[3]));
///     [[-1.0, (1.0 - k) / 2.0, 0.0, 0.0], [(1.0 - k) / 2.0, k, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]]
/// };
/// let t = stress_energy(tube, [0.0, 0.0, n - 1.0, 0.0], 1e-3)[0][0] * G / C.powi(4);
/// assert!((t + 0.0422211).abs() < 1e-7 && (t / -0.04217 - 1.0).abs() < 2e-3);
/// ```
pub fn stress_energy(metric: impl Fn([f64; 4]) -> [[f64; 4]; 4], x: [f64; 4], step: f64) -> [[f64; 4]; 4] {
    let gamma = christoffel(&metric, x, step);
    // dgamma[a][λ][μ][ν] = ∂_a Γ^λ_μν.
    let dgamma: [[[[f64; 4]; 4]; 4]; 4] = core::array::from_fn(|a| {
        let v = [-2.0, -1.0, 1.0, 2.0].map(|k| christoffel(&metric, nudge(x, a, k * step), step));
        core::array::from_fn(|l| core::array::from_fn(|m| core::array::from_fn(|n| central(v.map(|c| c[l][m][n]), step))))
    });
    // R_μν = ∂_λΓ^λ_μν − ∂_νΓ^λ_μλ + Γ^λ_λσΓ^σ_μν − Γ^λ_νσΓ^σ_μλ.
    let ricci: [[f64; 4]; 4] = core::array::from_fn(|m| {
        core::array::from_fn(|n| (0..4).map(|l| dgamma[l][l][m][n] - dgamma[n][l][m][l] + (0..4).map(|s| gamma[l][l][s] * gamma[s][m][n] - gamma[l][n][s] * gamma[s][m][l]).sum::<f64>()).sum())
    });
    let g = metric(x);
    let inv = inverse(g);
    let scalar: f64 = (0..16).map(|k| inv[k / 4][k % 4] * ricci[k / 4][k % 4]).sum();
    core::array::from_fn(|m| core::array::from_fn(|n| (ricci[m][n] - 0.5 * g[m][n] * scalar) / metrology::EINSTEIN_KAPPA))
}

/// What a stress-energy tensor means for the energy conditions, from [`energy_conditions`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EnergyConditions {
    /// `T_ab` (J/m³) in the orthonormal frame of the Eulerian observer, who moves normal to the
    /// slices of constant `x⁰`, with spatial axes from the coordinate axes by Gram–Schmidt:
    /// `[0][0]` is the energy density that observer sees, `[i][j]` the stresses, and `[0][i]` minus
    /// the momentum density times c.
    pub eulerian: [[f64; 4]; 4],
    /// The least `T_ab kᵃkᵇ` (J/m³) over the null vectors `k = (1, n)` of that frame, found exactly
    /// as a trust-region problem in the stress's eigenbasis. The null energy condition holds where
    /// it is ≥ 0; that sign does not depend on how `k` is scaled, its size does. Minima over
    /// timelike observers are not given: wherever this is negative they are unbounded, growing as γ².
    pub null: f64,
    /// For Hawking–Ellis type I stress-energy, `[ρ, p₁, p₂, p₃]` (J/m³): the density in its rest
    /// frame and the principal pressures in ascending order, from which the null (`ρ + p₁ ≥ 0`),
    /// weak, strong and dominant conditions follow exactly. `None` for types II–IV (a null or a
    /// complex eigenvector), which need [`EnergyConditions::null`] instead. Eigenvalues closer than
    /// 10⁻⁴ of the largest component count as one, which absorbs errors in `T` up to 10⁻⁸.
    pub type_one: Option<[f64; 4]>,
}

/// The energy conditions for the stress-energy `stress` (coordinate `T_μν`, J/m³, as from
/// [`stress_energy`]) at an event where the metric is `metric`: the components the Eulerian
/// observer measures, the null energy condition's minimum, and the Hawking–Ellis type I
/// decomposition. Le 2026's warpax finds that Eulerian checks alone miss 15–28% of violations
/// ([arXiv:2602.18023](https://arxiv.org/abs/2602.18023)); the null minimum and the eigenvalues do
/// not depend on an observer.
///
/// At Alcubierre's side wall (v = 1, σR = 8, G = c = 1) the Eulerian density is −1/(2π) while the
/// null minimum is −5/(2π), along the motion; the stress-energy is type I with rest-frame
/// `ρ = −(√15 − 2)/(4π)` and `ρ + p₁ = −√15/(2π)`. On the wall's front and back, and on Natário's
/// zero-expansion wall, it is type IV.
pub fn energy_conditions(metric: [[f64; 4]; 4], stress: [[f64; 4]; 4]) -> EnergyConditions {
    let inv = inverse(metric);
    let dot = |u: [f64; 4], v: [f64; 4]| (0..16).map(|k| metric[k / 4][k % 4] * u[k / 4] * v[k % 4]).sum::<f64>();
    // The Eulerian observer nᵘ = −αg^{μ0}; the coordinate axes ∂ᵢ are normal to it already.
    let alpha = (-1.0 / inv[0][0]).sqrt();
    let mut frame = [inv[0].map(|g| -alpha * g); 4];
    for i in 1..4 {
        let mut u = [0.0; 4];
        u[i] = 1.0;
        for e in &frame[1..i] {
            let p = dot(u, *e);
            u.iter_mut().zip(e).for_each(|(x, e)| *x -= p * e);
        }
        let norm = dot(u, u).sqrt();
        frame[i] = u.map(|x| x / norm);
    }
    let eulerian = core::array::from_fn(|a| core::array::from_fn(|b| (0..16).map(|k| frame[a][k / 4] * frame[b][k % 4] * stress[k / 4][k % 4]).sum()));
    EnergyConditions { eulerian, null: null_minimum(eulerian), type_one: hawking_ellis_type_one(eulerian) }
}

/// `min ρ + 2f·n + nᵀSn` over unit `n` for the orthonormal `T`. In the eigenbasis of `S`
/// (eigenvalues λ, `g = Qᵀf`) the minimiser is `n = −(Λ − μ)⁻¹g` with the `μ ≤ λ_min` that makes
/// `|n| = 1`, and the minimum is `ρ + μ − Σg²/(λ − μ)`; when `g` has no part along `λ_min` and
/// `|n| < 1` even there (the "hard case"), `μ = λ_min`.
fn null_minimum(t: [[f64; 4]; 4]) -> f64 {
    let stress: Vec<f64> = (0..9).map(|k| t[k / 3 + 1][k % 3 + 1]).collect();
    let (values, vectors) = geopotential::linalg::symmetric_eigen(&stress, 3);
    let parts: Vec<(f64, f64)> = (0..3).map(|k| (values[k], (0..3).map(|i| vectors[i * 3 + k] * t[0][i + 1]).sum::<f64>())).filter(|p| p.1 != 0.0).collect();
    let low = values.iter().copied().fold(f64::INFINITY, f64::min);
    let size = parts.iter().map(|p| p.1 * p.1).sum::<f64>().sqrt();
    // |n(μ)|² rises with μ, from at most 1 at λ_min − |g|.
    let length = |mu: f64| parts.iter().map(|(l, g)| (g / (l - mu)).powi(2)).sum::<f64>();
    let (mut lo, mut hi) = (low - size, low);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if length(mid) > 1.0 { hi = mid } else { lo = mid }
    }
    t[0][0] + lo - parts.iter().map(|(l, g)| g * g / (l - lo)).sum::<f64>()
}

/// Hawking–Ellis type I: `[ρ, p₁, p₂, p₃]` from the eigenvalues of `Tᵃ_b = η^{ac}T_cb`, if they are
/// real and one eigenspace holds a timelike eigenvector. That eigenspace is the one onto which the
/// Eulerian observer `e₀` projects as timelike: the spectral projections `Π(M − λ_d)/(λ − λ_d)`
/// (over the other eigenvalues) of a diagonalisable `M` are η-orthogonal and sum to the identity,
/// so their η-norms sum to −1.
fn hawking_ellis_type_one(t: [[f64; 4]; 4]) -> Option<[f64; 4]> {
    let m: [[f64; 4]; 4] = core::array::from_fn(|a| t[a].map(|x| if a == 0 { -x } else { x }));
    let scale = m.iter().flatten().fold(0.0f64, |big, x| big.max(x.abs()));
    if scale == 0.0 {
        return Some([0.0; 4]);
    }
    // Errors split a repeated eigenvalue, a defective one (types II and III) by the square root
    // of the error, so this tolerance absorbs errors up to 10⁻⁸ of the largest component.
    let tol = 1e-4 * scale;
    let values = eigenvalues(m)?;
    if values.iter().any(|v| v.1.abs() > tol) {
        return None;
    }
    let mut lambda = values.map(|v| v.0);
    lambda.sort_by(f64::total_cmp);
    let mut bounds = vec![0];
    bounds.extend((1..4).filter(|&k| lambda[k] - lambda[k - 1] > tol));
    bounds.push(4);
    let means: Vec<f64> = bounds.windows(2).map(|w| lambda[w[0]..w[1]].iter().sum::<f64>() / (w[1] - w[0]) as f64).collect();
    let apply = |v: [f64; 4], shift: f64| -> [f64; 4] { core::array::from_fn(|a| (0..4).map(|b| m[a][b] * v[b]).sum::<f64>() - shift * v[a]) };
    for (c, &mean) in means.iter().enumerate() {
        let others = means.iter().enumerate().filter(|&(d, _)| d != c);
        let v = others.clone().fold([1.0, 0.0, 0.0, 0.0], |v, (_, &shift)| apply(v, shift)).map(|x| x / others.clone().map(|(_, d)| mean - d).product::<f64>());
        let size2: f64 = v.iter().map(|x| x * x).sum();
        // The timelike projection has η-norm −γ² ≤ −1, the others ≥ 0 up to the errors.
        let norm = size2 - 2.0 * v[0] * v[0];
        if norm > -0.5 {
            continue;
        }
        // A defective eigenvalue's projection is no eigenvector.
        if apply(v, mean).iter().map(|x| x * x).sum::<f64>() > (10.0 * tol).powi(2) * size2 {
            return None;
        }
        // The eigenvalue v carries, vᵀTv/vᵀηv, picks ρ out of a cluster.
        let carried = (0..16).map(|k| t[k / 4][k % 4] * v[k / 4] * v[k % 4]).sum::<f64>() / norm;
        let k = (bounds[c]..bounds[c + 1]).min_by(|&i, &j| (lambda[i] - carried).abs().total_cmp(&(lambda[j] - carried).abs()))?;
        let p: Vec<f64> = (0..4).filter(|&i| i != k).map(|i| lambda[i]).collect();
        return Some([-lambda[k], p[0], p[1], p[2]]);
    }
    None
}

/// The eigenvalues `(re, im)` of a real 4×4 matrix, by Householder reduction to Hessenberg form
/// and the Francis double-shift QR iteration (EISPACK's `hqr`, as in Press et al., *Numerical
/// Recipes*, 2nd ed., §11.6), or `None` if that does not converge. Indices run from 1, as in `hqr`.
#[allow(clippy::needless_range_loop)]
fn eigenvalues(matrix: [[f64; 4]; 4]) -> Option<[(f64, f64); 4]> {
    const N: usize = 4;
    let sign = |a: f64, b: f64| if b >= 0.0 { a.abs() } else { -a.abs() };
    let mut a = [[0.0; N + 1]; N + 1];
    for i in 0..N {
        a[i + 1][1..].copy_from_slice(&matrix[i]);
    }
    for k in 1..=N - 2 {
        let norm = (k + 1..=N).map(|i| a[i][k] * a[i][k]).sum::<f64>().sqrt();
        if norm == 0.0 {
            continue;
        }
        let mut v = [0.0; N + 1];
        v[k + 1] = a[k + 1][k] + sign(norm, a[k + 1][k]);
        for i in k + 2..=N {
            v[i] = a[i][k];
        }
        let vv: f64 = v.iter().map(|x| x * x).sum();
        for j in 1..=N {
            let s = 2.0 * (k + 1..=N).map(|i| v[i] * a[i][j]).sum::<f64>() / vv;
            for i in k + 1..=N {
                a[i][j] -= s * v[i];
            }
        }
        for i in 1..=N {
            let s = 2.0 * (k + 1..=N).map(|j| a[i][j] * v[j]).sum::<f64>() / vv;
            for j in k + 1..=N {
                a[i][j] -= s * v[j];
            }
        }
    }
    for i in 3..=N {
        for j in 1..=i - 2 {
            a[i][j] = 0.0;
        }
    }
    let (mut wr, mut wi) = ([0.0; N + 1], [0.0; N + 1]);
    let mut anorm = 0.0;
    for i in 1..=N {
        for j in (i - 1).max(1)..=N {
            anorm += a[i][j].abs();
        }
    }
    let mut nn = N;
    let mut t = 0.0;
    while nn >= 1 {
        let mut its = 0;
        loop {
            let mut l = nn;
            while l >= 2 {
                let mut s = a[l - 1][l - 1].abs() + a[l][l].abs();
                if s == 0.0 {
                    s = anorm;
                }
                if a[l][l - 1].abs() + s == s {
                    a[l][l - 1] = 0.0;
                    break;
                }
                l -= 1;
            }
            let mut x = a[nn][nn];
            if l == nn {
                wr[nn] = x + t;
                wi[nn] = 0.0;
                nn -= 1;
                break;
            }
            let mut y = a[nn - 1][nn - 1];
            let mut w = a[nn][nn - 1] * a[nn - 1][nn];
            if l == nn - 1 {
                let p = 0.5 * (y - x);
                let q = p * p + w;
                let z = q.abs().sqrt();
                x += t;
                if q >= 0.0 {
                    let z = p + sign(z, p);
                    wr[nn - 1] = x + z;
                    wr[nn] = if z != 0.0 { x - w / z } else { x + z };
                    wi[nn - 1] = 0.0;
                    wi[nn] = 0.0;
                } else {
                    wr[nn - 1] = x + p;
                    wr[nn] = x + p;
                    wi[nn - 1] = -z;
                    wi[nn] = z;
                }
                nn -= 2;
                break;
            }
            if its == 30 {
                return None;
            }
            if its == 10 || its == 20 {
                t += x;
                for i in 1..=nn {
                    a[i][i] -= x;
                }
                let s = a[nn][nn - 1].abs() + a[nn - 1][nn - 2].abs();
                x = 0.75 * s;
                y = x;
                w = -0.4375 * s * s;
            }
            its += 1;
            let mut m = nn - 2;
            let (mut p, mut q, mut r, mut z);
            loop {
                z = a[m][m];
                r = x - z;
                let s = y - z;
                p = (r * s - w) / a[m + 1][m] + a[m][m + 1];
                q = a[m + 1][m + 1] - z - r - s;
                r = a[m + 2][m + 1];
                let s = p.abs() + q.abs() + r.abs();
                p /= s;
                q /= s;
                r /= s;
                if m == l {
                    break;
                }
                let u = a[m][m - 1].abs() * (q.abs() + r.abs());
                let v = p.abs() * (a[m - 1][m - 1].abs() + z.abs() + a[m + 1][m + 1].abs());
                if u + v == v {
                    break;
                }
                m -= 1;
            }
            for i in m + 2..=nn {
                a[i][i - 2] = 0.0;
                if i != m + 2 {
                    a[i][i - 3] = 0.0;
                }
            }
            for k in m..nn {
                if k != m {
                    p = a[k][k - 1];
                    q = a[k + 1][k - 1];
                    r = if k != nn - 1 { a[k + 2][k - 1] } else { 0.0 };
                    x = p.abs() + q.abs() + r.abs();
                    if x != 0.0 {
                        p /= x;
                        q /= x;
                        r /= x;
                    }
                }
                let s = sign((p * p + q * q + r * r).sqrt(), p);
                if s == 0.0 {
                    continue;
                }
                if k == m {
                    if l != m {
                        a[k][k - 1] = -a[k][k - 1];
                    }
                } else {
                    a[k][k - 1] = -s * x;
                }
                p += s;
                x = p / s;
                y = q / s;
                z = r / s;
                q /= p;
                r /= p;
                for j in k..=nn {
                    p = a[k][j] + q * a[k + 1][j];
                    if k != nn - 1 {
                        p += r * a[k + 2][j];
                        a[k + 2][j] -= p * z;
                    }
                    a[k + 1][j] -= p * y;
                    a[k][j] -= p * x;
                }
                for i in l..=nn.min(k + 3) {
                    p = x * a[i][k] + y * a[i][k + 1];
                    if k != nn - 1 {
                        p += z * a[i][k + 2];
                        a[i][k + 2] -= p * r;
                    }
                    a[i][k + 1] -= p * q;
                    a[i][k] -= p;
                }
            }
        }
    }
    Some(core::array::from_fn(|i| (wr[i + 1], wi[i + 1])))
}

/// The Regge–Wheeler potential `V` (units of `(c²/GM)²`) of a test field of spin `s` (0, 1 or 2)
/// and multipole `l` on the static `metric`, at radius `r` (units of `GM/c²`), each metric in its
/// own radius as in [`exact_deflection`]. The field obeys `∂²ψ/∂t² = ∂²ψ/∂r*² − Vψ`.
/// - Schwarzschild, areal `r`: `V = (1 − 2/r)[ℓ(ℓ + 1)/r² + 2(1 − s²)/r³]`, with `dr* = dr/(1 − 2/r)`
///   (Regge & Wheeler 1957 for `s = 2`).
/// - Puthoff's isotropic `r`: `V = e^{−4/r}[ℓ(ℓ + 1) + (s² − 1)((1 − 1/r)² − 1)]/r²`, with
///   `dr* = e^{2/r}dr` ([Boonserm et al. 2018](https://arxiv.org/abs/1805.03781) §8). This puts
///   general relativity's perturbation equations on the exponential geometry, their "inverse
///   Cowling approximation"; PV theory itself radiates a single scalar (**speculative**).
///
/// Each is a single barrier, so neither has echoes. For `ℓ = s = 2` it peaks at 0.15129 at
/// `r = (9 + √17)/4` = 3.2808 (GR) and at 0.13188 at isotropic `r = 2.4376` (PV); for spin 1 it
/// is the photon's, peaking at the photon sphere at `ℓ(ℓ + 1)/b_c²`. Evolving `ψ` in time and
/// fitting its ringing gives `Mω = 0.3448 − 0.0851i` on the PV geometry, against Schwarzschild's
/// 0.3737 − 0.0890i ([Berti, Cardoso & Starinets 2009](https://arxiv.org/abs/0905.2975)), 7.7% lower.
pub fn regge_wheeler(metric: Metric, s: u32, l: u32, r: f64) -> f64 {
    let (ell, spin) = (f64::from(l * (l + 1)), f64::from(s * s));
    match metric {
        Metric::GeneralRelativity => (1.0 - 2.0 / r) * (ell / (r * r) + 2.0 * (1.0 - spin) / r.powi(3)),
        Metric::PolarizableVacuum => (-4.0 / r).exp() * (ell + (spin - 1.0) * ((1.0 - 1.0 / r).powi(2) - 1.0)) / (r * r),
    }
}

/// The dielectric function `K` of Puthoff's polarizable vacuum around a static `mass` (kg) of
/// `charge` (C), at isotropic radius `r` (m), with `g₀₀ = 1/K`: `√K = cosh(k/r) + (a/k)sinh(k/r)`
/// with `a = GM/c²`, `b² = GQ²/(4πε₀c⁴)` and `k = √(a² − b²)`, or the cos and sin of
/// `k = √(b² − a²)` when `b > a` ([Puthoff 2002](https://arxiv.org/abs/gr-qc/9909037), eq. 65)
/// (**speculative**). It matches Reissner–Nordström to second order,
/// `g₀₀ = 1 − 2x + (2 + b²/a²)x²` with `x = a/r`. `K` drops below 1, the faster clocks Puthoff's
/// craft needs, only when `Q/M` exceeds [`crate::gem::extremal_charge_to_mass`] (`b > a`): inside
/// `r = 1.654a` for `b = 2a`, where the solution ends at `√K = 0`, `r = 0.827a`.
pub fn charged_pv_dielectric(mass: f64, charge: f64, r: f64) -> f64 {
    let a = G * mass / (C * C);
    let b2 = G * charge * charge / (4.0 * PI * metrology::EPS_0 * C.powi(4));
    // (k/r)², negative when b > a.
    let q = (a * a - b2) / (r * r);
    let (even, odd) = if q > 0.0 {
        (q.sqrt().cosh(), q.sqrt().sinh() / q.sqrt())
    } else if q < 0.0 {
        ((-q).sqrt().cos(), (-q).sqrt().sin() / (-q).sqrt())
    } else {
        (1.0, 1.0)
    };
    (even + a / r * odd).powi(2)
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

impl Metric {
    /// The time component `g₀₀` at `x = GM/(c²r)`, with `r` the isotropic radius:
    /// `((1 − x/2)/(1 + x/2))²` in general relativity and `e^{−2x}` (that is, `1/K`) in the
    /// polarizable vacuum ([Puthoff 2002](https://arxiv.org/abs/gr-qc/9909037), eqs 42 and 45).
    /// Both are `1 − 2x + 2x²`, so they agree at first post-Newtonian order, and part at `x³`:
    /// `−3/2` against `−4/3`. A negative `x` is a negative mass.
    pub fn g00(self, x: f64) -> f64 {
        match self {
            Metric::GeneralRelativity => ((1.0 - x / 2.0) / (1.0 + x / 2.0)).powi(2),
            Metric::PolarizableVacuum => (-2.0 * x).exp(),
        }
    }

    /// The spatial component `−g_ii` at `x = GM/(c²r)`, isotropic: `(1 + x/2)⁴` in general
    /// relativity and `e^{2x}` (`K`) in the polarizable vacuum (Puthoff 2002 eqs 43 and 46). Both are
    /// `1 + 2x`, with `x²` coefficients `3/2` and `2`. In the polarizable vacuum
    /// `g₀₀·spatial = 1`, so a faster clock always comes with longer rulers.
    pub fn spatial(self, x: f64) -> f64 {
        match self {
            Metric::GeneralRelativity => (1.0 + x / 2.0).powi(4),
            Metric::PolarizableVacuum => (2.0 * x).exp(),
        }
    }

    /// The critical impact parameter in units of `GM/c²`, the smallest `r√(spatial/g₀₀)`: `3√3`
    /// in general relativity and `2e` in the polarizable vacuum, whose photon sphere is at
    /// isotropic `r = 2GM/c²` ([Boonserm et al. 2018](https://arxiv.org/abs/1805.03781), §7). Light
    /// passing closer is captured ([`exact_deflection`]).
    pub fn critical_impact(self) -> f64 {
        match self {
            Metric::GeneralRelativity => 27f64.sqrt(),
            Metric::PolarizableVacuum => 2.0 * core::f64::consts::E,
        }
    }

    /// The angular diameter (rad) of the shadow of a static `mass` (kg) at `distance` (m): twice
    /// the critical impact parameter seen from `D`, which is [`shadow_diameter`] in general
    /// relativity. The polarizable vacuum's is `2e/√27 − 1` = 4.63% larger, `δ = +0.0463` in
    /// [`shadow_deviation`]'s terms; the literature's "5%" rounds this. Spin would shrink a
    /// Kerr shadow, but no rotating-PV shadow has been computed.
    ///
    /// ```
    /// use gravitics::spacetime::{Metric, shadow_deviation};
    /// use metrology::{C, G, PARSEC, SOLAR_MASS};
    ///
    /// let uas = std::f64::consts::PI / 180.0 / 3600.0 * 1e-6;
    /// // Sgr A* with the EHT's VLTI prior, 4.297×10⁶ M☉ at 8277 pc: 53.254 µas against 55.718.
    /// let (mass, distance) = (4.297e6 * SOLAR_MASS, 8277.0 * PARSEC);
    /// let gr = Metric::GeneralRelativity.shadow_diameter(mass, distance) / uas;
    /// let pv = Metric::PolarizableVacuum.shadow_diameter(mass, distance) / uas;
    /// assert!((gr - 53.254).abs() < 5e-4 && (pv - 55.718).abs() < 5e-4);
    /// let theta_g = G * mass / (C * C * distance) / uas;
    /// assert!((shadow_deviation(pv, theta_g) - 0.046267).abs() < 5e-7);
    /// // M87*: δ = −0.01 ± 0.17 with the stellar-dynamics prior (EHT M87* Paper VI) puts PV at
    /// // 0.33σ and Schwarzschild at 0.06σ, so M87* does not constrain it.
    /// let sigma = |metric: Metric| (metric.critical_impact() / 27f64.sqrt() - 1.0 + 0.01) / 0.17;
    /// assert!((sigma(Metric::PolarizableVacuum) - 0.33).abs() < 5e-3 && (sigma(Metric::GeneralRelativity) - 0.06).abs() < 5e-3);
    /// ```
    ///
    /// The EHT's own δ comparisons are in [EHT Sgr A* Paper VI](https://arxiv.org/abs/2311.09484)
    /// and [EHT M87* Paper VI](https://arxiv.org/abs/1906.11243).
    pub fn shadow_diameter(self, mass: f64, distance: f64) -> f64 {
        shadow_diameter(mass, distance) * (self.critical_impact() / 27f64.sqrt())
    }
}

/// Bozza's strong-deflection coefficients `(ā, b̄)`: near the critical impact parameter the
/// deflection is `α ≈ −ā ln(b/b_c − 1) + b̄` ([Bozza 2002](https://arxiv.org/abs/gr-qc/0208075)).
/// General relativity has `(1, ln[216(7 − 4√3)] − π)` = `(1, −0.40023004)` (Bozza §3.1). The
/// polarizable vacuum has `ā = 1` exactly too, so its photon subrings also shrink by `e^π`, but
/// `b̄ = −0.46581756`: our value, as no published one could be read. Only `b_c` and `b̄` tell the
/// geometries apart; `b̄` moves the subrings' offsets from the critical curve by `e^{Δb̄}` = 0.936.
///
/// It evaluates Bozza's `ā = R/(2√β_m)` and `b̄ = −π + b_R + ā ln(2β_m/y_m)` from the photon
/// sphere's `y_m = A_m`, `C_m`, `C′_m`, `C″_m` and `A″_m`, with `b_R` his regular integral over
/// `z = (A − y_m)/(1 − y_m)`.
pub fn strong_deflection(metric: Metric) -> (f64, f64) {
    let e = core::f64::consts::E;
    // Areal x = 3 for GR, isotropic x = 2 for PV.
    let (ym, c, c1, c2, a2) = match metric {
        Metric::GeneralRelativity => (1.0 / 3.0, 9.0, 6.0, 2.0, -4.0 / 27.0),
        Metric::PolarizableVacuum => (1.0 / e, 4.0 * e, 2.0 * e, e, -1.0 / (4.0 * e)),
    };
    let beta = c * (1.0 - ym).powi(2) * (c2 * ym - c * a2) / (2.0 * ym * ym * c1 * c1);
    // A·B = 1 and C·A′ = 2 in both metrics, so Bozza's R(z, x_m) is the constant (1 − y_m)√C_m.
    let r = (1.0 - ym) * c.sqrt();
    // y_m − A·C_m/C, which has a double root at z = 0, written without cancellation there.
    let radicand = |z: f64| match metric {
        Metric::GeneralRelativity => z * z * (1.0 - 2.0 * z / 3.0),
        Metric::PolarizableVacuum => {
            let s = (e - 1.0) * z;
            let q = (1.0 + s) * s.ln_1p() - s;
            q * (2.0 - q) / e
        }
    };
    // The subtracted integrand is analytic on [0, 1]; 32 nodes reach round-off.
    let (x, w) = geopotential::gauss_legendre(32);
    let b_r = r * x.iter().zip(&w).map(|(x, w)| {
        let z = 0.5 * (1.0 + x);
        0.5 * w * (1.0 / radicand(z).sqrt() - 1.0 / (beta.sqrt() * z))
    }).sum::<f64>();
    let abar = r / (2.0 * beta.sqrt());
    (abar, -PI + b_r + abar * (2.0 * beta / ym).ln())
}

/// The innermost stable circular orbit `(areal radius, GMΩ/c³, E)`, with the radius in units of
/// `GM/c²` and `E` the energy per unit rest energy, so `1 − E` is the accretion efficiency.
/// General relativity: `(6, 6^{−3/2}, √(8/9))`. The polarizable vacuum: isotropic `r = (3 + √5)GM/c²`,
/// areal `r e^{1/r}` = 6.337940, `MΩ = e^{−2/r}/(r√(r − 1))` = 0.0633326 and
/// `E = e^{−1/r}√((r − 1)/(r − 2))` = 0.9452132 ([Boonserm et al. 2018](https://arxiv.org/abs/1805.03781), §7;
/// [Makukov & Mychelkin 2018](https://arxiv.org/abs/1809.05290)). PV's `Ω` is 6.92% lower and its
/// efficiency 5.479% against 5.719%, which a slightly retrograde Kerr spin would mimic:
/// **speculative**, not a clean test without an independent spin.
///
/// ```
/// use gravitics::spacetime::{Metric, isco};
/// use metrology::{C, G, SOLAR_MASS};
///
/// // The gravitational-wave frequency Ω/π at the ISCO of 10 M☉: 439.7 Hz against 409.3 Hz.
/// let f = |metric| isco(metric).1 / std::f64::consts::PI * C.powi(3) / (G * 10.0 * SOLAR_MASS);
/// assert!((f(Metric::GeneralRelativity) - 439.7).abs() < 0.05 && (f(Metric::PolarizableVacuum) - 409.3).abs() < 0.05);
/// ```
pub fn isco(metric: Metric) -> (f64, f64, f64) {
    match metric {
        Metric::GeneralRelativity => {
            let r = 6.0f64;
            (r, r.powf(-1.5), (1.0 - 2.0 / r) / (1.0 - 3.0 / r).sqrt())
        }
        Metric::PolarizableVacuum => {
            let r = 3.0 + 5f64.sqrt();
            (r * (1.0 / r).exp(), (-2.0 / r).exp() / (r * (r - 1.0).sqrt()), (-1.0 / r).exp() * ((r - 1.0) / (r - 2.0)).sqrt())
        }
    }
}

/// The periastron ratio `K = Ω_φ/Ω_r` of a near-circular orbit at radius `r` (units of `GM/c²`;
/// areal in general relativity, isotropic in the polarizable vacuum, as in [`exact_deflection`]):
/// `(1 − 6/r)^{−1/2}` and `r/√(r² − 6r + 4)`, infinite at the [`isco`]. The orbit advances by
/// `2π(K − 1)` per radial period. In terms of the frequency, `x = (GMΩ/c³)^{2/3}` with
/// `MΩ = r^{−3/2}` or `e^{−2/r}/(r√(r − 1))`, they are `1 + 3x + 13.5x² + 67.5x³` and
/// `1 + 3x + 14.5x² + 76.5x³`: PV adds `2πx²` per orbit, 0.35 µas per century for Mercury.
/// Ours, from the research notes' PV tests; **speculative** and untested.
pub fn periastron_ratio(metric: Metric, r: f64) -> f64 {
    let c = match metric {
        Metric::GeneralRelativity => 0.0,
        Metric::PolarizableVacuum => 4.0,
    };
    r / (r * r - 6.0 * r + c).sqrt()
}

/// The redshift `z` of light from a static emitter at areal radius `R`, with `u = GM/(c²R)`:
/// `1/√(1 − 2u) − 1` in general relativity, and `e^w − 1` in the polarizable vacuum, where
/// `w e^{−w} = u` gives its isotropic `GM/(c²r)` on the principal branch (up to the throat,
/// `u = 1/e`; beyond it, NaN). They part at `u³/6`: 0.2346620 against 0.2371177 at `u = 0.172`, a
/// 1% difference at a neutron-star surface, though PV has no published stellar interior. Ours,
/// from the research notes; **speculative**.
pub fn surface_redshift(metric: Metric, u: f64) -> f64 {
    let g00 = match metric {
        Metric::GeneralRelativity => 1.0 - 2.0 * u,
        Metric::PolarizableVacuum => {
            if u > (-1.0f64).exp() {
                return f64::NAN;
            }
            let (mut lo, mut hi) = (0.0, 1.0);
            while hi - lo > f64::EPSILON * hi {
                let mid = 0.5 * (lo + hi);
                if mid * (-mid).exp() < u { lo = mid } else { hi = mid }
            }
            metric.g00(0.5 * (lo + hi))
        }
    };
    1.0 / redshift_factor(g00) - 1.0
}

/// The second post-Newtonian Shapiro delay (s) of the polarizable vacuum minus general
/// relativity's, one way past a `mass` (kg) at impact parameter `b` (m) with the emitter and
/// receiver fixed: `(π/4)(GM/c²)²/(bc)`. That is 8.21 ps at the solar limb, and in Cassini's
/// Doppler observable it mimics `δγ = (π/8)GM/(bc²)` = 5.2×10⁻⁷ at 1.6 R☉, against Cassini's
/// σ = 2.3×10⁻⁵. Ours, derived in the research notes; untested.
pub fn shapiro_2pn_difference(mass: f64, b: f64) -> f64 {
    let m = G * mass / (C * C);
    PI / 4.0 * m * m / (b * C)
}

/// The rate of a static clock against a remote one, `√g₀₀`: Puthoff's Table 1 rows 1–3, where a
/// remote observer sees local intervals as `Δt/√g₀₀`, frequencies as `ω√g₀₀` and energies as
/// `E√g₀₀` ([Puthoff 2010](https://arxiv.org/abs/1204.2184)). His §3.1 example, ten seconds
/// remote judged as five locally, is `g₀₀ = 0.25`. Between two heights it gives the gravitational
/// redshift; for a `g₀₀` that close to 1 use [`clock_shift`], since an `f64` holds `g₀₀` itself only
/// to the nearest step of 2.2×10⁻¹⁶. **Established.**
pub fn redshift_factor(g00: f64) -> f64 {
    g00.sqrt()
}

/// The fractional clock shift `√g₀₀ − 1` from the metric's excess `g₀₀ − 1`, computed as
/// `expm1(ln_1p(excess)/2)` so weak fields keep every digit: Pound and Rebka's 22.5 m tower,
/// `excess = 2gh/c²`, gives `gh/c²` = 2.456×10⁻¹⁵ exactly, where [`redshift_factor`] of
/// `1 + excess` is 0.55% off.
pub fn clock_shift(excess: f64) -> f64 {
    (0.5 * excess.ln_1p()).exp_m1()
}

/// Puthoff's "effective mass" factor `spatial/√g₀₀` (Table 1 row 6 of
/// [Puthoff 2010](https://arxiv.org/abs/1204.2184); `K^{3/2}` in the polarizable vacuum,
/// [Puthoff 2002](https://arxiv.org/abs/gr-qc/9909037) eq. 9), with `spatial = −g₁₁`: the
/// coefficient of `½mv²` in the Lagrangian `−mc²√(g₀₀ − spatial·v²/c²)` for the coordinate
/// velocity `v`. It is bookkeeping, **not an invariant**: it depends on the coordinates, and it
/// rises near a mass, `1 + 3x + 4x²` in general relativity and `1 + 3x + 4.5x²` in PV, while the
/// conserved energy of a body at rest, `mc²√g₀₀`, falls.
pub fn effective_mass_factor(g00: f64, spatial: f64) -> f64 {
    spatial / g00.sqrt()
}

/// The mass (kg, negative) whose exterior has `g₀₀ = 1 + excess` at areal `radius` (m): Puthoff's
/// "spacetime-engineered" metric (`g₀₀ > 1`, faster clocks) is exactly the exterior of a negative
/// mass. General relativity: `−excess·Rc²/(2G)`. The polarizable vacuum: `−yRe^y c²/G` with
/// `y = ln(1 + excess)/2`, the same in the weak field and less in the strong: 3 times less for a
/// clock running 17.6 times fast. Ours: the paper gives no number. With ordinary matter and fields
/// `√g₀₀` is subharmonic and cannot exceed 1, so the region needs a strong-energy-condition
/// violation. The mathematics is established; as an explanation of a craft it is **fringe**.
///
/// ```
/// use gravitics::spacetime::{Metric, engineered_mass};
/// use metrology::{C, G};
///
/// // A clock speed-up of 10⁻¹⁸ at 10 m, the reach of optical clocks, needs 6.7×10⁹ kg of
/// // negative mass, and it pulls outward with c²δ/(2R) = 4.5 mm/s²...
/// let mass = engineered_mass(1e-18, 10.0, Metric::GeneralRelativity);
/// assert!((mass / -6.733e9 - 1.0).abs() < 1e-4);
/// let pull = |excess: f64| G * -engineered_mass(excess, 10.0, Metric::PolarizableVacuum) / 100.0;
/// assert!((pull(1e-18) / (C * C * 1e-18 / 20.0) - 1.0).abs() < 1e-12 && (pull(1e-18) - 4.5e-3).abs() < 1e-5);
/// // ...so a 1 µGal gravimeter at 10 m resolves δ = 2gR/c² = 2.2×10⁻²⁴, a million times less.
/// let delta = 2.0 * 1e-8 * 10.0 / (C * C);
/// assert!((delta - 2.2e-24).abs() < 5e-26 && (pull(delta) / 1e-8 - 1.0).abs() < 1e-12);
/// ```
///
/// Hiding the region from gravimeters would need a ± mass double layer of dipole density
/// `δc²/(8πG)`, 5.4×10⁷ kg/m for `δ = 10⁻¹⁸`, negative inside (ours, derived in the notes, not
/// found in a source).
pub fn engineered_mass(excess: f64, radius: f64, metric: Metric) -> f64 {
    match metric {
        Metric::GeneralRelativity => -excess * radius * C * C / (2.0 * G),
        Metric::PolarizableVacuum => {
            let y = excess.ln_1p() / 2.0;
            -y * radius * y.exp() * C * C / G
        }
    }
}

/// The extra transit time (s) of light across a region of `length` `L` (m) filled with
/// electromagnetic energy density `u` (J/m³, `B²/2μ₀` or `ε₀E²/2`), Puthoff, Davis & Maccone's
/// Levi-Civita effect, the field route Puthoff 2010 cites: `(L/c)·4πGuL²/(3c⁴)`
/// ([Puthoff, Davis & Maccone 2005](https://arxiv.org/abs/physics/0403064), eqs 8 and 15). At the
/// centre `K − 1 = πGuL²/c⁴`, which is `3cτ/(4L)` (eqs 6 and 13): 1.03×10⁻³⁸ for 1 T over 1 m,
/// with a delay of 4.6×10⁻⁴⁷ s. The sign is `K > 1`: field energy slows light like mass, the
/// opposite of what the engineered metric needs. Their length scale
/// `a = c²/(B√(4πG/μ₀))`, with `K − 1 = L²/(8a²)`, prints as 3.484×10¹⁸/B m, which takes
/// `c = 3×10⁸` m/s; the exact `c` gives 3.4789×10¹⁸. **Established** at this order, where PV and
/// general relativity agree.
pub fn levi_civita_delay(energy_density: f64, length: f64) -> f64 {
    length / C * 4.0 * PI * G * energy_density * length * length / (3.0 * C.powi(4))
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

    /// The Taylor coefficient of `f` at 0 that follows the `known` ones, by one Richardson step
    /// from `h` to `h/2`, so its error is of order `h²` times the coefficient two places on.
    fn next_coefficient(f: impl Fn(f64) -> f64, known: &[f64], h: f64) -> f64 {
        let rest = |x: f64| (f(x) - known.iter().rev().fold(0.0, |sum, c| sum * x + c)) / x.powi(known.len() as i32);
        2.0 * rest(h / 2.0) - rest(h)
    }

    #[test]
    fn exact_deflection_carries_both_weak_field_series() {
        // PV's 4π, 48 and 64π against GR's 15π/4, 128/3 and 3465π/64 (the research's series,
        // confirmed without cancellation in PR #47), from exact_deflection at b = 500 to 2000 GM/c².
        let m = metrology::SOLAR_MASS;
        let rg = G * m / (C * C);
        for (metric, [c2, c3, c4]) in [
            (Metric::GeneralRelativity, [15.0 * PI / 4.0, 128.0 / 3.0, 3465.0 * PI / 64.0]),
            (Metric::PolarizableVacuum, [4.0 * PI, 48.0, 64.0 * PI]),
        ] {
            let alpha = |x: f64| exact_deflection(m, rg / x, metric).unwrap();
            // Each estimate is off by about c_(k+2)h²/2: 10⁻⁴ for c₂ and 4×10⁻⁴ for c₃ at
            // h = 10⁻³. For c₄, α's absolute round-off of 4×10⁻¹⁶ (it is 2∫ − π) over h⁴ needs
            // h = 2×10⁻³, where the remainder is 8×10⁻³. The gaps PV − GR are 0.79, 5.3 and 31.
            assert!((next_coefficient(alpha, &[0.0, 4.0], 1e-3) - c2).abs() < 2e-4, "{metric:?}");
            assert!((next_coefficient(alpha, &[0.0, 4.0, c2], 1e-3) - c3).abs() < 1e-3, "{metric:?}");
            assert!((next_coefficient(alpha, &[0.0, 4.0, c2, c3], 2e-3) - c4).abs() < 2e-2, "{metric:?}");
        }
    }

    #[test]
    fn the_metrics_in_isotropic_coordinates() {
        let (gr, pv) = (Metric::GeneralRelativity, Metric::PolarizableVacuum);
        // Both g₀₀ = 1 − 2x + 2x², then −3/2 against −4/3; both spatial = 1 + 2x, then 3/2 against 2
        // (Puthoff 2002 eqs 42–46; Ibison 2006 eq. 3).
        for (metric, cubic, square) in [(gr, -1.5, 1.5), (pv, -4.0 / 3.0, 2.0)] {
            assert!((next_coefficient(|x| metric.g00(x), &[1.0, -2.0, 2.0], 1e-3) - cubic).abs() < 1e-5, "{metric:?}");
            assert!((next_coefficient(|x| metric.spatial(x), &[1.0, 2.0], 1e-3) - square).abs() < 1e-5, "{metric:?}");
        }
        // Isotropic Schwarzschild is areal 1 − 2GM/(c²R) with R = r(1 + x/2)².
        for x in [0.01, 0.3, 1.0] {
            assert!(close(gr.g00(x), 1.0 - 2.0 * x / (1.0 + x / 2.0).powi(2), 1e-14), "x = {x}");
            assert!(close(pv.g00(x) * pv.spatial(x), 1.0, 1e-15) && close(pv.g00(-x), pv.spatial(x), 1e-15));
        }
        // The critical impact parameter is the least r√(spatial/g₀₀), by golden-section search
        // over the isotropic radius; and exact_deflection captures light inside it.
        for metric in [gr, pv] {
            let b = |r: f64| r * (metric.spatial(1.0 / r) / metric.g00(1.0 / r)).sqrt();
            let (mut lo, mut hi) = (1.0, 10.0);
            for _ in 0..100 {
                let (a, c) = (hi - 0.618 * (hi - lo), lo + 0.618 * (hi - lo));
                if b(a) < b(c) { hi = c } else { lo = a }
            }
            assert!(close(b(0.5 * (lo + hi)), metric.critical_impact(), 1e-14), "{metric:?}");
            let rg = G * metrology::SOLAR_MASS / (C * C);
            let at = |scale: f64| exact_deflection(metrology::SOLAR_MASS, scale * metric.critical_impact() * rg, metric);
            assert!(at(1.0 - 1e-9).is_none() && at(1.0 + 1e-9).is_some(), "{metric:?}");
        }
        assert!(close(pv.critical_impact() / gr.critical_impact() - 1.0, 0.046_267, 1e-5));
        // The GR shadow is shadow_diameter's, bit for bit.
        let (mass, distance) = (4.297e6 * metrology::SOLAR_MASS, 8277.0 * metrology::PARSEC);
        assert_eq!(gr.shadow_diameter(mass, distance), shadow_diameter(mass, distance));
    }

    #[test]
    fn strong_deflection_both_ways() {
        // Bozza §3.1: ā = 1 and b̄ = ln[216(7 − 4√3)] − π = −0.40023004 for Schwarzschild.
        let (a, b) = strong_deflection(Metric::GeneralRelativity);
        let exact = (216.0 * (7.0 - 4.0 * 3f64.sqrt())).ln() - PI;
        assert!((a - 1.0).abs() < 1e-15 && (b - exact).abs() < 1e-14 && (b + 0.400_230_04).abs() < 5e-9, "{a} {b}");
        // PV: ā = 1 and b̄ = −0.46581756 (the research's PV notes, two methods to 2×10⁻⁸). A
        // long-double recomputation of Bozza's integral gives −0.465817563343; the notes'
        // ten-digit −0.4658175637 lost 4×10⁻¹⁰ to cancellation near the photon sphere.
        let (a, b) = strong_deflection(Metric::PolarizableVacuum);
        assert!((a - 1.0).abs() < 1e-15 && (b + 0.465_817_56).abs() < 5e-9, "{a} {b}");
        // Against exact_deflection at b_c(1 + 10⁻⁶): 13.4152855578 and 13.3496988228 (the notes),
        // with α + ln ε within 6×10⁻⁶ of b̄, the O(ε ln ε) remainder.
        let m = metrology::SOLAR_MASS;
        let rg = G * m / (C * C);
        let eps = 1e-6;
        for (metric, alpha) in [(Metric::GeneralRelativity, 13.415_285_557_8), (Metric::PolarizableVacuum, 13.349_698_822_8)] {
            let at = exact_deflection(m, metric.critical_impact() * (1.0 + eps) * rg, metric).unwrap();
            assert!((at - alpha).abs() < 5e-11, "{metric:?}");
            assert!((at + eps.ln() - strong_deflection(metric).1).abs() < 6e-6, "{metric:?}");
        }
    }

    #[test]
    fn the_innermost_stable_orbit_both_ways() {
        // GR: 6GM/c², MΩ = 6^(−3/2) = 0.0680414, E = √(8/9) = 0.9428090.
        let (r, omega, energy) = isco(Metric::GeneralRelativity);
        assert!(r == 6.0 && close(omega, 6f64.powf(-1.5), 1e-15) && close(energy, (8.0f64 / 9.0).sqrt(), 1e-15));
        assert!((omega - 0.068_041_4).abs() < 5e-8 && (energy - 0.942_809_0).abs() < 5e-8);
        // PV at isotropic (3 + √5)GM/c²: areal 6.337940 and L = 3.523216 (Boonserm et al. §7, who
        // print 6.337940263 and 3.523216438, 2×10⁻⁹ and 10⁻⁹ below the closed forms' 6.3379402649
        // and 3.5232164391), MΩ = 0.0633326 (Makukov & Mychelkin), E = 0.9452132.
        let (r, omega, energy) = isco(Metric::PolarizableVacuum);
        assert!((r - 6.337_940).abs() < 5e-7 && (omega - 0.063_332_6).abs() < 5e-8 && (energy - 0.945_213_2).abs() < 5e-8);
        let iso = 3.0 + 5f64.sqrt();
        let pv = Metric::PolarizableVacuum;
        // L = C Ω u^t with u^t = E/g₀₀, and C = r²·spatial.
        let l = iso * iso * pv.spatial(1.0 / iso) * omega * energy / pv.g00(1.0 / iso);
        assert!((l - 3.523_216).abs() < 5e-7, "{l}");
        // Ω is 6.92% lower and the efficiency 5.479% against 5.719%.
        let gr = isco(Metric::GeneralRelativity);
        assert!((omega / gr.1 - 1.0 + 0.069_20).abs() < 5e-6);
        assert!(((1.0 - energy) * 100.0 - 5.479).abs() < 5e-4 && ((1.0 - gr.2) * 100.0 - 5.719).abs() < 5e-4);
        // It is where the circular orbits' E(r) = e^(−1/r)√((r − 1)/(r − 2)) is least, and where
        // the periastron ratio diverges.
        let e = |r: f64| (-1.0 / r).exp() * ((r - 1.0) / (r - 2.0)).sqrt();
        assert!(close(e(iso), energy, 1e-15) && e(iso - 1e-3) > energy && e(iso + 1e-3) > energy);
        assert!(periastron_ratio(pv, iso * (1.0 + 1e-12)) > 1e5 && periastron_ratio(Metric::GeneralRelativity, 6.0 * (1.0 + 1e-12)) > 1e5);
    }

    #[test]
    fn periastron_ratio_both_ways() {
        // The research's numerical Ω_φ/Ω_r (pv_metric_tests.py §5): GR at areal 7, 10 and 20,
        // where it is (1 − 6x)^(−1/2); PV at isotropic 6, 9 and 19.
        for (metric, r, k) in [
            (Metric::GeneralRelativity, 7.0, 2.645_75),
            (Metric::GeneralRelativity, 10.0, 1.581_14),
            (Metric::GeneralRelativity, 20.0, 1.195_23),
            (Metric::PolarizableVacuum, 6.0, 3.0),
            (Metric::PolarizableVacuum, 9.0, 1.616_45),
            (Metric::PolarizableVacuum, 19.0, 1.199_27),
        ] {
            assert!((periastron_ratio(metric, r) - k).abs() < 5e-6, "{metric:?} r = {r}");
        }
        // In x = (MΩ)^(2/3): 1 + 3x + 13.5x² + 67.5x³ against 1 + 3x + 14.5x² + 76.5x³. The x³
        // coefficient by extrapolating (K − 1 − 3x − c₂x²)/x³ linearly in x from r = 10³ and 2×10³.
        for (metric, c2, c3) in [(Metric::GeneralRelativity, 13.5, 67.5), (Metric::PolarizableVacuum, 14.5, 76.5)] {
            let point = |r: f64| {
                let omega = match metric {
                    Metric::GeneralRelativity => r.powf(-1.5),
                    Metric::PolarizableVacuum => (-2.0 / r).exp() / (r * (r - 1.0).sqrt()),
                };
                let x = omega.powf(2.0 / 3.0);
                (x, (periastron_ratio(metric, r) - 1.0 - 3.0 * x - c2 * x * x) / x.powi(3))
            };
            let ((x1, d1), (x2, d2)) = (point(1e3), point(2e3));
            // The remainder is of order c₅x₁x₂, about 10⁻³.
            assert!(((d2 * x1 - d1 * x2) / (x1 - x2) - c3).abs() < 5e-3, "{metric:?}");
        }
    }

    #[test]
    fn surface_redshift_both_ways() {
        // At GM/(c²R) = 0.172: 0.2346620 against 0.2371177 (the research's PV notes, §8).
        let (gr, pv) = (Metric::GeneralRelativity, Metric::PolarizableVacuum);
        assert!((surface_redshift(gr, 0.172) - 0.234_662_0).abs() < 5e-8 && (surface_redshift(pv, 0.172) - 0.237_117_7).abs() < 5e-8);
        // At S2's pericentre, 1/1400, they differ by 6.096×10⁻¹¹, close to u³/6.
        let u = 1.0 / 1400.0;
        assert!(close(surface_redshift(pv, u) - surface_redshift(gr, u), 6.096e-11, 1e-4) && close(u.powi(3) / 6.0, 6.074e-11, 1e-3));
        // PV's throat is at areal eGM/c², where the isotropic radius is GM/c² and z = e − 1.
        assert!(close(surface_redshift(pv, (-1.0f64).exp()), core::f64::consts::E - 1.0, 1e-7) && surface_redshift(pv, 0.4).is_nan());
    }

    #[test]
    fn the_second_order_shapiro_difference() {
        // (π/4)(GM/c²)²/(bc): 8.211 ps at the solar limb and 5.132 ps at Cassini's 1.6 R☉.
        let (sun, limb) = (metrology::SOLAR_MASS, 6.957e8);
        assert!((shapiro_2pn_difference(sun, limb) * 1e12 - 8.211).abs() < 5e-4);
        assert!((shapiro_2pn_difference(sun, 1.6 * limb) * 1e12 - 5.132).abs() < 5e-4);
    }

    #[test]
    fn puthoffs_table_one() {
        // §3.1: ten seconds remote are five locally where g₀₀ = 0.25.
        assert_eq!(10.0 * redshift_factor(0.25), 5.0);
        // Pound and Rebka's gh/c² = 2.456×10⁻¹⁵ over 22.5 m (g = 9.81 m/s²), the bottom clock as
        // reference. f64 holds 1 + 2gh/c² only to the nearest of its 22 ulps above 1, so
        // redshift_factor reproduces the shift to 0.55%; clock_shift, from g₀₀ − 1, is exact.
        let excess = 2.0 * 9.81 * 22.5 / (C * C);
        assert!(((redshift_factor(1.0 + excess) - 1.0) / 2.456e-15 - 1.0).abs() < 6e-3);
        assert!(close(clock_shift(excess), excess / 2.0, 1e-14) && (clock_shift(excess) / 2.456e-15 - 1.0).abs() < 2e-4);
        assert!(close(clock_shift(-0.75), redshift_factor(0.25) - 1.0, 1e-15));
        // At the Sun's surface, −GM/(c²R) = −2.1225×10⁻⁶ to first order.
        let u = metrology::GM_SUN / (C * C * 6.957e8);
        assert!(close(redshift_factor(1.0 - 2.0 * u) - 1.0, -2.1225e-6, 5e-5));
        // The mass factor is K^(3/2) in PV, 8 at K = 4, and the ½v² coefficient of
        // −√(g₀₀ − spatial·v²): 1.912366 at (0.7, 1.6) (the research's §B).
        assert_eq!(effective_mass_factor(0.25, 4.0), 8.0);
        assert!((effective_mass_factor(0.7, 1.6) - 1.912_366).abs() < 5e-7);
        let lagrangian = |v: f64| -(0.7 - 1.6 * v * v).sqrt();
        let h = 1e-4;
        assert!(close((lagrangian(h) - 2.0 * lagrangian(0.0) + lagrangian(-h)) / (h * h), effective_mass_factor(0.7, 1.6), 1e-7));
        // Near a mass it rises as 1 + 3x + 4x² in GR and 1 + 3x + 4.5x² in PV (the research's §C).
        for (metric, square) in [(Metric::GeneralRelativity, 4.0), (Metric::PolarizableVacuum, 4.5)] {
            let mass = |x: f64| effective_mass_factor(metric.g00(x), metric.spatial(x));
            assert!((next_coefficient(mass, &[1.0, 3.0], 1e-3) - square).abs() < 1e-5, "{metric:?}");
        }
    }

    #[test]
    fn negative_mass_for_an_engineered_metric() {
        let (gr, pv) = (Metric::GeneralRelativity, Metric::PolarizableVacuum);
        // 10⁻¹⁸ at 10 m: 6.733×10⁹ kg in both (puthoff_metric_checks §G).
        assert!(close(engineered_mass(1e-18, 10.0, gr), -6.733e9, 1e-4) && close(engineered_mass(1e-18, 10.0, pv), engineered_mass(1e-18, 10.0, gr), 1e-15));
        // A 300 K body glowing at 550 nm, √g₀₀ = 17.56: 2.070×10³⁰ kg in GR, 6.777×10²⁹ in PV.
        let excess = (2.897_771_955e-3 / 300.0 / 550e-9f64).powi(2) - 1.0;
        assert!(close(engineered_mass(excess, 10.0, gr), -2.070e30, 1e-4) && close(engineered_mass(excess, 10.0, pv), -6.777e29, 1e-4));
        // Back again: the PV mass at areal R = r e^(−x), x = G|M|/(c²r), so x e^x = G|M|/(c²R), and
        // there g₀₀ = e^(2x) = 1 + excess.
        let target = -G * engineered_mass(excess, 10.0, pv) / (C * C * 10.0);
        let (mut lo, mut hi) = (0.0f64, 10.0);
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if mid * mid.exp() < target { lo = mid } else { hi = mid }
        }
        assert!(close(pv.g00(-0.5 * (lo + hi)), 1.0 + excess, 1e-13));
    }

    #[test]
    fn the_levi_civita_effect() {
        let mu0 = metrology::MU_0;
        // 1 T over 1 m: a delay of 4.594×10⁻⁴⁷ s and K − 1 = 3cτ/(4L) = 1.033×10⁻³⁸; a 10⁹ T
        // pulsar field over 10 km, 1.033×10⁻¹² (puthoff_metric_checks §H).
        let tau = levi_civita_delay(1.0 / (2.0 * mu0), 1.0);
        assert!(close(tau, 4.594e-47, 1.1e-4) && close(3.0 * C * tau / 4.0, 1.033e-38, 5e-4));
        assert!(close(3.0 * C * levi_civita_delay(1e18 / (2.0 * mu0), 1e4) / 4e4, 1.033e-12, 5e-4));
        // Puthoff, Davis & Maccone eq. B-2: a = c²/(B√(4πG/μ₀)), with K − 1 = L²/(8a²). Their
        // 3.484×10¹⁸ T·m is 3.4837×10¹⁸ with c = 3×10⁸; the exact c gives 3.4789×10¹⁸.
        let a = |c: f64| c * c / (4.0 * PI * G / mu0).sqrt();
        assert!((a(C) / 3.4789e18 - 1.0).abs() < 1.5e-5 && (a(3e8) / 3.4837e18 - 1.0).abs() < 1.5e-5);
        assert!(close(3.0 * C * tau / 4.0, 1.0 / (8.0 * a(C) * a(C)), 1e-13));
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

#[cfg(test)]
mod warp_tests {
    use super::*;
    use metrology::{EINSTEIN_KAPPA, HBAR, SOLAR_MASS, YEAR};

    fn close(a: f64, b: f64, rel: f64) -> bool {
        (a / b - 1.0).abs() < rel
    }

    /// J/m³ to G = c = 1 with lengths in metres: G/c⁴.
    const GEOMETRIC: f64 = EINSTEIN_KAPPA / (8.0 * PI);

    /// The shape function of Alcubierre's tanh wall.
    fn alcubierre_shape(r: f64, radius: f64, sigma: f64) -> f64 {
        ((sigma * (r + radius)).tanh() - (sigma * (r - radius)).tanh()) / (2.0 * (sigma * radius).tanh())
    }

    /// Alcubierre's shift `v = (v f(r_s), 0, 0)` at the event `x = (ct, x, y, z)`, the bubble
    /// centred at `x = v ct`.
    fn alcubierre_shift(v: f64, radius: f64, sigma: f64, x: [f64; 4]) -> [f64; 3] {
        let r = ((x[1] - v * x[0]).powi(2) + x[2] * x[2] + x[3] * x[3]).sqrt();
        [v * alcubierre_shape(r, radius, sigma), 0.0, 0.0]
    }

    /// Natário's zero-expansion shift along x: `−2vf cosθ e_r + v(2f + rf′) sinθ e_θ`, with
    /// `f = (1 − f_Alcubierre)/2` (gr-qc/0110086).
    fn natario_shift(v: f64, radius: f64, sigma: f64, x: [f64; 4]) -> [f64; 3] {
        let [_, x, y, z] = x;
        let (r, rho) = ((x * x + y * y + z * z).sqrt(), y.hypot(z));
        let (f, slope) = (0.5 * (1.0 - alcubierre_shape(r, radius, sigma)), -0.5 * alcubierre_slope(r, radius, sigma));
        let (c, s) = (x / r, rho / r);
        let (radial, polar) = ([x / r, y / r, z / r], [-s, c * y / rho, c * z / rho]);
        core::array::from_fn(|i| -2.0 * v * f * c * radial[i] + v * (2.0 * f + r * slope) * s * polar[i])
    }

    /// The metric `−dt² + |dx − v dt|²` of a shift `v` (Santiago, Schuster & Visser's convention).
    fn shift_metric(v: [f64; 3]) -> [[f64; 4]; 4] {
        core::array::from_fn(|a| {
            core::array::from_fn(|b| match (a, b) {
                (0, 0) => -1.0 + v[0] * v[0] + v[1] * v[1] + v[2] * v[2],
                (0, i) | (i, 0) => -v[i - 1],
                (i, j) => f64::from(u8::from(i == j)),
            })
        })
    }

    /// The Jacobian `∂v_i/∂x_j` of a shift and `ℒₙ(∇·v) = (∂₀ + v·∇)(∇·v)`, by nested central differences.
    fn shift_derivatives(shift: &impl Fn([f64; 4]) -> [f64; 3], x: [f64; 4], h: f64) -> ([[f64; 3]; 3], f64) {
        let jacobian = |y: [f64; 4]| -> [[f64; 3]; 3] { core::array::from_fn(|i| core::array::from_fn(|j| central([-2.0, -1.0, 1.0, 2.0].map(|k| shift(nudge(y, j + 1, k * h))[i]), h))) };
        let divergence = |y: [f64; 4]| {
            let j = jacobian(y);
            j[0][0] + j[1][1] + j[2][2]
        };
        let v = shift(x);
        let lie = (0..4).map(|a| if a == 0 { 1.0 } else { v[a - 1] } * central([-2.0, -1.0, 1.0, 2.0].map(|k| divergence(nudge(x, a, k * h))), h)).sum();
        (jacobian(x), lie)
    }

    /// The energy conditions of a metric at `x`, in G = c = 1.
    fn conditions(metric: impl Fn([f64; 4]) -> [[f64; 4]; 4], x: [f64; 4], h: f64) -> EnergyConditions {
        let c = energy_conditions(metric(x), stress_energy(&metric, x, h));
        EnergyConditions {
            eulerian: c.eulerian.map(|row| row.map(|t| t * GEOMETRIC)),
            null: c.null * GEOMETRIC,
            type_one: c.type_one.map(|t| t.map(|t| t * GEOMETRIC)),
        }
    }

    #[test]
    fn alcubierre_density_and_the_capacitance_floor() {
        // At the side wall (0, R, 0) with σR = 8 and v = 1: −1/(2π), since f′(R) = −σ/2 up to
        // tanh 8 = 1 − 2.3×10⁻⁷ (exactly −tanh⁴16/(2π tanh²8)).
        let side = alcubierre_energy_density(1.0, 8.0, 1.0, [0.0, 1.0, 0.0]) * GEOMETRIC;
        assert!(close(side, -1.0 / (2.0 * PI), 4.6e-7) && !close(side, -1.0 / (2.0 * PI), 4.4e-7));
        assert!(close(side, -16f64.tanh().powi(4) / (2.0 * PI * 8f64.tanh().powi(2)), 1e-14));
        // Warp Factory's case: −6.7751×10³⁵ J/m³ at (0, 300 m, 0). Zero on the axis and at the centre.
        assert!(close(alcubierre_energy_density(300.0, 0.015, 0.1, [0.0, 300.0, 0.0]), -6.7751e35, 1e-4));
        assert!(alcubierre_energy_density(1.0, 8.0, 1.0, [0.0; 3]) == 0.0 && alcubierre_energy_density(1.0, 8.0, 1.0, [1.0, 0.0, 0.0]) == 0.0);
        // Bobrick & Martire's floor for a 100 m interior at light speed: 1.122×10²⁸ kg, 5.9 M_J.
        assert!(close(alcubierre_minimum_energy(100.0, 1.0), -1.122_16e28, 1e-5));
        assert!(close(alcubierre_minimum_energy(100.0, 1.0) / -1.898_13e27, 5.912, 1e-3));
    }

    #[test]
    fn alcubierre_total_energy_of_the_tanh_wall() {
        // The density's volume integral, with ∫sin²θ dΩ = 8π/3, by panels of its own.
        let (x, w) = geopotential::gauss_legendre(32);
        let radial: f64 = (0..300)
            .map(|k| {
                x.iter()
                    .zip(&w)
                    .map(|(x, w)| {
                        let r = 0.02 * (k as f64 + 0.5 + 0.5 * x);
                        0.01 * w * r * r * alcubierre_energy_density(1.0, 8.0, 1.0, [0.0, r, 0.0])
                    })
                    .sum::<f64>()
            })
            .sum();
        let total = alcubierre_tanh_energy(1.0, 8.0, 1.0);
        assert!(close(8.0 * PI / 3.0 * radial / (C * C), total, 1e-12));
        // Alcubierre's σR = 8 profile needs 2.680104 times the capacitance floor (BM's "factor of three").
        assert!(close(total / alcubierre_minimum_energy(1.0, 1.0), 2.680_104, 1e-6));
        // Thin walls: −(v²/36)R²σ c²/G (0.3762 M☉ at R = 100 m, σ = 2 m⁻¹, up to 8×10⁻⁶), 2/3 of
        // Pfenning & Ford's Eq. 28 with their Δ = (1 + tanh²σR)²/(2σ tanh σR), and 1/36 of Lobo & Visser's "−v²R²σ".
        let thin = alcubierre_tanh_energy(100.0, 2.0, 1.0);
        assert!(close(thin / (-1e4 * 2.0 / 36.0 * C * C / G), 1.000_008_06, 1e-8));
        assert!(close(thin / -SOLAR_MASS, 0.3762, 1e-4));
        let delta = (1.0 + 200f64.tanh().powi(2)).powi(2) / (2.0 * 2.0 * 200f64.tanh());
        assert!((thin / alcubierre_energy(100.0, delta, 1.0) - 2.0 / 3.0).abs() < 1e-6);
        assert!(close(-1e4 * 20.0 * C * C / G / alcubierre_tanh_energy(100.0, 20.0, 1.0), 36.0, 1e-6));
        // Over all walls the least is 0.667477 R (v²/12)c²/G at σR = 1.2290.
        let (best, at) = (1000..=1500).map(|k| (alcubierre_tanh_energy(1.0, k as f64 / 1000.0, 1.0), k as f64 / 1000.0)).fold((f64::NEG_INFINITY, 0.0), |a, b| if b.0 > a.0 { b } else { a });
        assert!(close(best / alcubierre_minimum_energy(1.0, 1.0), 0.667_477, 1e-6) && (at - 1.229).abs() < 1.5e-3);
    }

    #[test]
    fn shift_drives_at_alcubierre_and_natario_walls() {
        // Alcubierre's side wall, v = 1, R = 1, σ = 8: ρ = −1/(2π) and ρ + p̄ = −7/(6π) = −0.37136,
        // both up to tanh 8 = 1 − 2.3×10⁻⁷. Analytically ∂v_x/∂y = vf′(R) and ℒₙK = v²(f − 1)f′/R there.
        let shift = |x: [f64; 4]| alcubierre_shift(1.0, 1.0, 8.0, x);
        let wall = [0.0, 0.0, 1.0, 0.0];
        let (jacobian, lie) = shift_derivatives(&shift, wall, 1e-3);
        let (f, slope) = (alcubierre_shape(1.0, 1.0, 8.0), alcubierre_slope(1.0, 1.0, 8.0));
        assert!((jacobian[0][1] - slope).abs() < 1e-8 && (lie - (f - 1.0) * slope).abs() < 1e-7);
        assert!(jacobian.iter().flatten().filter(|j| j.abs() > 1e-9).count() == 1);
        let rho = shift_drive_density(jacobian) * GEOMETRIC;
        let null = shift_drive_null_sum(jacobian, lie) * GEOMETRIC;
        // Fourth-order differences of step h leave (hσ)⁴ = 4×10⁻⁹ in the slope.
        assert!((rho - alcubierre_energy_density(1.0, 8.0, 1.0, [0.0, 1.0, 0.0]) * GEOMETRIC).abs() < 1e-8);
        assert!(close(null, -7.0 / (6.0 * PI), 1e-6) && (null + 0.37136).abs() < 5e-6);
        // The other sign of K gives −0.26526, which the Einstein tensor contradicts (below).
        assert!((shift_drive_null_sum(jacobian.map(|r| r.map(|j| -j)), -lie) * GEOMETRIC + 0.26526).abs() < 5e-6);
        // Natário (v = 1, R = 1, σ = 4) at (0.6, 0.85, 0): −0.0405272, his
        // −(v²/8π)[3f′²cos²θ + (f′ + rf″/2)²sin²θ], with no expansion.
        let shift = |x: [f64; 4]| natario_shift(1.0, 1.0, 4.0, x);
        let point = [0.0, 0.6, 0.85, 0.0];
        let (jacobian, _) = shift_derivatives(&shift, point, 1e-3);
        let r = 0.6f64.hypot(0.85);
        let slope = -0.5 * alcubierre_slope(r, 1.0, 4.0);
        let curve = -0.5 * (alcubierre_slope(r + 1e-4, 1.0, 4.0) - alcubierre_slope(r - 1e-4, 1.0, 4.0)) / 2e-4;
        let cos2 = 0.36 / (r * r);
        let natario = -1.0 / (8.0 * PI) * (3.0 * slope * slope * cos2 + (slope + r * curve / 2.0).powi(2) * (1.0 - cos2));
        let rho = shift_drive_density(jacobian) * GEOMETRIC;
        assert!(close(rho, natario, 1e-7) && (rho + 0.040_527_2).abs() < 1e-7);
        assert!((jacobian[0][0] + jacobian[1][1] + jacobian[2][2]).abs() < 1e-9);
    }

    #[test]
    fn wormhole_shells_and_quantum_inequality_bounds() {
        // A 1 m shell around nothing: −1.926×10⁴³ J/m², −2.693×10²⁷ kg in all, pressure +c⁴/(4πG).
        let (sigma, p) = thin_shell_wormhole(1.0, 0.0);
        assert!(close(sigma, -1.926_18e43, 1e-5) && close(4.0 * PI * sigma / (C * C), -2.693_18e27, 1e-5));
        assert!(close(p, 9.630_91e42, 1e-5) && sigma + p < 0.0);
        assert!(close(4.0 * PI * thin_shell_wormhole(1.0, 0.1 * C * C / G).0 / (C * C), -2.408_85e27, 1e-5));
        // Static shells keep σ′ = −(2/a)(σ + p) (Poisson & Visser), here at a = 4GM/c².
        let mass = C * C / G;
        let h = 1e-4;
        let slope = (thin_shell_wormhole(4.0 + h, mass).0 - thin_shell_wormhole(4.0 - h, mass).0) / (2.0 * h);
        let (sigma, p) = thin_shell_wormhole(4.0, mass);
        assert!(close(slope, -0.5 * (sigma + p), 1e-8));
        // Ford & Roman's bounds at f = 0.01, against the closed forms in l_P = √(ħG/c³).
        let planck = (HBAR * G / C.powi(3)).sqrt();
        assert!(close(wormhole_throat_bound(0.01), 7.8971e-32, 1e-4) && close(wormhole_throat_bound(0.01), (0.75 / PI).sqrt() * planck / 1e-4, 1e-13));
        let band = |r0: f64| (3.0 * r0 * planck * planck / (8.0 * PI * 1e-8)).cbrt();
        assert!(close(wormhole_band_bound(1.0, 0.01), 1.460_94e-21, 1e-5) && close(wormhole_band_bound(1.0, 0.01), band(1.0), 1e-13));
        // FR's "0.2 fermi" for a light-year throat is 0.31 fm; 1.43×10⁻¹⁴ m for 10⁵ ly.
        assert!(close(wormhole_band_bound(C * YEAR, 0.01), 3.0899e-16, 1e-4) && close(wormhole_band_bound(1e5 * C * YEAR, 0.01), 1.4342e-14, 1e-4));
        // The Gaussian bound at 1 ns: −1.859×10⁻²⁶ J/m³. Fewster & Eveson's (1/16π²)∫((√w)″)² by
        // quadrature: 3/(64π²) for the Gaussian, half of Ford & Roman's 3/(32π²), and 27/(2048π²),
        // 9/64 of it, for the Lorentzian (t = tan θ).
        assert!(close(quantum_inequality_gaussian(1e-9), -1.858_90e-26, 1e-5));
        let (x, w) = geopotential::gauss_legendre(64);
        let gauss = (0..24)
            .map(|k| {
                x.iter()
                    .zip(&w)
                    .map(|(x, w)| {
                        let t = -12.0 + k as f64 + 0.5 + 0.5 * x;
                        0.5 * w * ((t * t - 1.0) * (-t * t / 2.0).exp()).powi(2)
                    })
                    .sum::<f64>()
            })
            .sum::<f64>()
            / PI.sqrt();
        let lorentz = x
            .iter()
            .zip(&w)
            .map(|(x, w)| {
                let (s, c) = (PI / 2.0 * x).sin_cos();
                PI / 2.0 * w * (2.0 * s * s - c * c).powi(2) * c.powi(4)
            })
            .sum::<f64>()
            / PI;
        let ford_roman = 3.0 / (32.0 * PI * PI);
        assert!(close(gauss / (16.0 * PI * PI) / ford_roman, quantum_inequality_gaussian(1.0) / quantum_inequality(1.0), 1e-13));
        assert!(close(lorentz / (16.0 * PI * PI) / ford_roman, 9.0 / 64.0, 1e-13));
    }

    #[test]
    fn eigenvalues_of_general_four_by_four_matrices() {
        // B = [[2, −3], [3, 2]] ⊕ [[1, 5], [0, −4]] turned by a reflection Q = I − 2uuᵀ: 2 ± 3i, 1, −4.
        let u = [0.5, -0.5, 0.1, 0.7];
        let q: [[f64; 4]; 4] = core::array::from_fn(|i| core::array::from_fn(|j| if i == j { 1.0 } else { 0.0 } - 2.0 * u[i] * u[j]));
        let b = [[2.0, -3.0, 0.0, 0.0], [3.0, 2.0, 0.0, 0.0], [0.0, 0.0, 1.0, 5.0], [0.0, 0.0, 0.0, -4.0]];
        let m: [[f64; 4]; 4] = core::array::from_fn(|i| core::array::from_fn(|j| (0..16).map(|k| q[i][k / 4] * b[k / 4][k % 4] * q[j][k % 4]).sum()));
        let mut values = eigenvalues(m).unwrap();
        values.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
        for (got, want) in values.iter().zip([(-4.0, 0.0), (1.0, 0.0), (2.0, -3.0), (2.0, 3.0)]) {
            assert!((got.0 - want.0).abs() < 1e-13 && (got.1 - want.1).abs() < 1e-13, "{values:?}");
        }
        // A matrix already triangular, and the zero matrix.
        let mut values = eigenvalues([[3.0, 1.0, 4.0, 1.0], [0.0, -5.0, 9.0, 2.0], [0.0, 0.0, 6.0, 5.0], [0.0, 0.0, 0.0, 0.5]]).unwrap().map(|v| v.0);
        values.sort_by(f64::total_cmp);
        assert_eq!(values, [-5.0, 0.5, 3.0, 6.0]);
        assert_eq!(eigenvalues([[0.0; 4]; 4]).unwrap(), [(0.0, 0.0); 4]);
    }

    #[test]
    fn hawking_ellis_types_of_known_stress_energy() {
        // An anisotropic fluid ρ = 1, p = (−3, 0.5, 2), boosted by β = 0.5 along (1, 2, 2)/3: type I,
        // and the null condition fails (ρ + p₁ = −2).
        let n = [1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0];
        let (beta, gamma) = (0.5, 1.0 / 0.75f64.sqrt());
        let boost: [[f64; 4]; 4] = core::array::from_fn(|a| {
            core::array::from_fn(|b| match (a, b) {
                (0, 0) => gamma,
                (0, i) | (i, 0) => gamma * beta * n[i - 1],
                (i, j) => f64::from(u8::from(i == j)) + (gamma - 1.0) * n[i - 1] * n[j - 1],
            })
        });
        let rest = [1.0, -3.0, 0.5, 2.0];
        let t: [[f64; 4]; 4] = core::array::from_fn(|a| core::array::from_fn(|b| (0..4).map(|c| boost[c][a] * rest[c] * boost[c][b]).sum()));
        let found = hawking_ellis_type_one(t).unwrap();
        assert!(found.iter().zip(rest).all(|(a, b)| (a - b).abs() < 1e-12), "{found:?}");
        assert!(null_minimum(t) < 0.0);
        // Vacuum energy T = Λ·diag(1, −1, −1, −1): ρ = Λ, p = −Λ, all degenerate.
        assert_eq!(hawking_ellis_type_one([[2.0, 0.0, 0.0, 0.0], [0.0, -2.0, 0.0, 0.0], [0.0, 0.0, -2.0, 0.0], [0.0, 0.0, 0.0, -2.0]]), Some([2.0, -2.0, -2.0, -2.0]));
        // Null dust T = k_a k_b, k = (1, 0.6, 0.8, 0): type II, and the null condition holds at 0.
        let k = [-1.0, 0.6, 0.8, 0.0];
        let dust: [[f64; 4]; 4] = core::array::from_fn(|a| core::array::from_fn(|b| k[a] * k[b]));
        assert!(hawking_ellis_type_one(dust).is_none() && null_minimum(dust).abs() < 1e-15);
        // A flux larger than (ρ + p)/2: complex eigenvalues ±i, type IV.
        assert!(hawking_ellis_type_one([[0.0, 1.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]]).is_none());
        // The exact null minimum against a fine sampling of directions.
        let t = [[0.3, -0.4, 0.2, 0.1], [-0.4, 1.0, 0.3, -0.2], [0.2, 0.3, -0.5, 0.6], [0.1, -0.2, 0.6, 0.2]];
        let sampled = (0..400 * 400)
            .map(|k| {
                let (z, phi) = (-1.0 + (k / 400) as f64 / 199.5, PI * (k % 400) as f64 / 200.0);
                let s = (1.0 - z * z).max(0.0).sqrt();
                let e = [s * phi.cos(), s * phi.sin(), z];
                t[0][0] + 2.0 * (0..3).map(|i| t[0][i + 1] * e[i]).sum::<f64>() + (0..9).map(|i| t[i / 3 + 1][i % 3 + 1] * e[i / 3] * e[i % 3]).sum::<f64>()
            })
            .fold(f64::INFINITY, f64::min);
        let exact = null_minimum(t);
        assert!(exact <= sampled && sampled - exact < 1e-3, "{exact} {sampled}");
    }

    #[test]
    fn stress_energy_of_flat_space_and_of_schwarzschild() {
        // Minkowski exactly, and in Rindler's accelerated coordinates (a = 1 m⁻¹): no stress.
        let flat = stress_energy(|_| [[-1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]], [0.3, 1.0, -2.0, 0.5], 1e-3);
        assert_eq!(flat, [[0.0; 4]; 4]);
        let rindler = |x: [f64; 4]| [[-(1.0 + x[1]).powi(2), 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]];
        let t = conditions(rindler, [0.0, 0.5, 0.0, 0.0], 1e-3);
        assert!(t.eulerian.iter().flatten().all(|t| t.abs() < 1e-11) && t.null.abs() < 1e-11);
        // Schwarzschild in isotropic coordinates, GM/c² = 1 m, at r = 3 m and 6 m: a residual of
        // order the rounding 10⁻¹⁶/h², against curvature 1/r².
        let schwarzschild = |x: [f64; 4]| {
            let a = 0.5 / (x[1] * x[1] + x[2] * x[2] + x[3] * x[3]).sqrt();
            let s = (1.0 + a).powi(4);
            [[-((1.0 - a) / (1.0 + a)).powi(2), 0.0, 0.0, 0.0], [0.0, s, 0.0, 0.0], [0.0, 0.0, s, 0.0], [0.0, 0.0, 0.0, s]]
        };
        for r in [3.0, 6.0] {
            let t = conditions(schwarzschild, [0.0, r / 3f64.sqrt(), r / 3f64.sqrt(), r / 3f64.sqrt()], 1e-3);
            assert!(t.eulerian.iter().flatten().all(|t| t.abs() < 1e-10), "{:?}", t.eulerian);
        }
    }

    #[test]
    fn energy_conditions_of_warp_drives_from_their_metrics() {
        // Alcubierre (v = 1, R = 1, σ = 8) at the side wall (0, R, 0), t = 0: the density, ρ + p̄ of
        // Santiago, Schuster & Visser, the null minimum −5/(2π) along the motion, and type I with
        // ρ = −(√15 − 2)/(4π), p = (−(2 + √15), −3, 1)/(4π): closed forms up to tanh 8 = 1 − 2.3×10⁻⁷.
        let alcubierre = |x: [f64; 4]| shift_metric(alcubierre_shift(1.0, 1.0, 8.0, x));
        let c = conditions(alcubierre, [0.0, 0.0, 1.0, 0.0], 1e-3);
        let rho = c.eulerian[0][0];
        // Against the closed form, to the truncation (hσ)⁴ = 4×10⁻⁹ (of a curvature of order 1).
        assert!((rho - alcubierre_energy_density(1.0, 8.0, 1.0, [0.0, 1.0, 0.0]) * GEOMETRIC).abs() < 1e-8);
        let (jacobian, lie) = shift_derivatives(&|x| alcubierre_shift(1.0, 1.0, 8.0, x), [0.0, 0.0, 1.0, 0.0], 1e-3);
        let mean = rho + (c.eulerian[1][1] + c.eulerian[2][2] + c.eulerian[3][3]) / 3.0;
        assert!(close(mean, shift_drive_null_sum(jacobian, lie) * GEOMETRIC, 1e-7) && (mean + 0.37136).abs() < 5e-6);
        assert!(close(c.null, -5.0 / (2.0 * PI), 1e-6));
        let root = 15f64.sqrt();
        let want = [-(root - 2.0), -(2.0 + root), -3.0, 1.0].map(|x| x / (4.0 * PI));
        let found = c.type_one.unwrap();
        assert!(found.iter().zip(want).all(|(a, b)| close(*a, b, 1e-6)), "{found:?}");
        assert!(close(found[0] + found[1], -root / (2.0 * PI), 1e-6) && (found[0] + found[1] + 0.616_40).abs() < 5e-6);
        // At 45° and on the front of the wall it is type IV; there the null minima are −0.762533 and
        // −0.676539 (a numerical minimisation over directions of the research's evaluator).
        #[allow(clippy::approx_constant)] // the research's point, not 1/√2
        let diagonal = 0.7071;
        for (point, null) in [([0.0, diagonal, diagonal, 0.0], -0.762_533), ([0.0, 1.0, 1e-3, 0.0], -0.676_539)] {
            let c = conditions(alcubierre, point, 1e-3);
            let analytic = alcubierre_energy_density(1.0, 8.0, 1.0, [point[1], point[2], 0.0]) * GEOMETRIC;
            assert!((c.eulerian[0][0] - analytic).abs() < 1e-8);
            assert!(c.type_one.is_none() && (c.null - null).abs() < 2e-6, "{point:?}: {}", c.null);
        }
        // Warp Factory's comoving bubble in SI (v = 0.1, R = 300 m, σ = 0.015 m⁻¹, g_0x = v(1 − f)):
        // −6.7751×10³⁵ J/m³ at (0, 300 m, 0).
        let comoving = |x: [f64; 4]| {
            let b = 0.1 * (1.0 - alcubierre_shape((x[1] * x[1] + x[2] * x[2] + x[3] * x[3]).sqrt(), 300.0, 0.015));
            [[-1.0 + b * b, b, 0.0, 0.0], [b, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]]
        };
        let x = [0.0, 0.0, 300.0, 0.0];
        let c = energy_conditions(comoving(x), stress_energy(comoving, x, 0.5));
        assert!(close(c.eulerian[0][0], alcubierre_energy_density(300.0, 0.015, 0.1, [0.0, 300.0, 0.0]), 1e-7));
        // Natário (v = 1, R = 1, σ = 4) at (0.6, 0.85, 0): his −0.0405272, type IV, null minimum −1.225409.
        let natario = |x: [f64; 4]| shift_metric(natario_shift(1.0, 1.0, 4.0, x));
        let point = [0.0, 0.6, 0.85, 0.0];
        let c = conditions(natario, point, 1e-3);
        let (jacobian, _) = shift_derivatives(&|x| natario_shift(1.0, 1.0, 4.0, x), point, 1e-3);
        assert!(close(c.eulerian[0][0], shift_drive_density(jacobian) * GEOMETRIC, 1e-7) && (c.eulerian[0][0] + 0.040_527_2).abs() < 1e-7);
        assert!(c.type_one.is_none() && (c.null + 1.225_409).abs() < 2e-6);
    }

    #[test]
    fn energy_conditions_of_a_morris_thorne_throat() {
        // b = r₀²/r with r₀ = 1 m, Φ = 0, at r = 1.2: ρ = p_r = −1/(8πr⁴) and p_t = +1/(8πr⁴) (Ford &
        // Roman 1996), type I; the null minimum is ρ + p_r, radially.
        let wormhole = |x: [f64; 4]| {
            let r2 = x[1] * x[1] + x[2] * x[2] + x[3] * x[3];
            let extra = 1.0 / (1.0 - 1.0 / r2) - 1.0;
            core::array::from_fn(|a| {
                core::array::from_fn(|b| match (a, b) {
                    (0, 0) => -1.0,
                    (0, _) | (_, 0) => 0.0,
                    (i, j) => f64::from(u8::from(i == j)) + extra * x[i] * x[j] / r2,
                })
            })
        };
        let r = 1.2f64;
        let c = conditions(wormhole, [0.0, r / 2f64.sqrt(), r / 2f64.sqrt(), 0.0], 1e-3);
        let a = 1.0 / (8.0 * PI * r.powi(4));
        assert!(close(c.eulerian[0][0], -a, 1e-8) && close(c.null, -2.0 * a, 1e-8) && (c.null + 0.038_376_5).abs() < 1e-7);
        let found = c.type_one.unwrap();
        assert!(found.iter().zip([-a, -a, a, a]).all(|(x, y)| close(*x, y, 1e-6)), "{found:?}");
    }

    #[test]
    fn regge_wheeler_potentials() {
        let (gr, pv) = (Metric::GeneralRelativity, Metric::PolarizableVacuum);
        // The peak, by golden section.
        let peak = |metric, s, l, lo: f64, hi: f64| {
            let (mut lo, mut hi) = (lo, hi);
            for _ in 0..200 {
                let (a, b) = (hi - 0.618_034 * (hi - lo), lo + 0.618_034 * (hi - lo));
                if regge_wheeler(metric, s, l, a) < regge_wheeler(metric, s, l, b) { lo = a } else { hi = b }
            }
            (0.5 * (lo + hi), regge_wheeler(metric, s, l, 0.5 * (lo + hi)))
        };
        // Schwarzschild ℓ = s = 2: the root (9 + √17)/4 of dV/dr, V = 0.1512867.
        let (at, top) = peak(gr, 2, 2, 2.5, 5.0);
        assert!((at - (9.0 + 17f64.sqrt()) / 4.0).abs() < 1e-7 && (top - 0.151_286_70).abs() < 1e-8);
        // Against the research's grid values (0.15129, 0.13188, 0.24712, 0.22869) and our peaks; its
        // PV ℓ = s = 2 radius, 2.4365, is a grid's 2.43756.
        for (metric, s, l, lo, value, radius) in [
            (gr, 2, 2, 2.5, 0.151_286_699_6, 3.280_776_458),
            (pv, 2, 2, 1.5, 0.131_876_761_2, 2.437_564_895),
            (gr, 0, 2, 2.5, 0.247_116_781_7, 2.951_714_821),
            (pv, 0, 2, 1.2, 0.228_686_877_1, 1.928_849_147),
        ] {
            let (at, top) = peak(metric, s, l, lo, 5.0);
            assert!((top - value).abs() < 1e-9 && (at - radius).abs() < 1e-6, "{metric:?} {s} {l}: {top} at {at}");
        }
        assert!((peak(pv, 2, 2, 1.5, 5.0).1 - 0.13188).abs() < 5e-6);
        // Spin 1 peaks at the photon spheres, ℓ(ℓ + 1)/27 at r = 3 and ℓ(ℓ + 1)/(2e)² at r = 2 (Boonserm §8).
        for l in 1..4 {
            let ell = f64::from(l * (l + 1));
            assert!(close(regge_wheeler(gr, 1, l, 3.0), ell / 27.0, 1e-15) && close(regge_wheeler(pv, 1, l, 2.0), ell / (4.0 * core::f64::consts::E.powi(2)), 1e-15));
            assert!((peak(gr, 1, l, 2.5, 5.0).0 - 3.0).abs() < 1e-7 && (peak(pv, 1, l, 1.2, 5.0).0 - 2.0).abs() < 1e-7);
        }
        // PV spin 0, ℓ = 0: a peak at (1 + 1/√3) and a trough at (1 − 1/√3) (Boonserm §8).
        let slope = |r: f64| (regge_wheeler(pv, 0, 0, r + 1e-6) - regge_wheeler(pv, 0, 0, r - 1e-6)) / 2e-6;
        for r in [1.0 + 1.0 / 3f64.sqrt(), 1.0 - 1.0 / 3f64.sqrt()] {
            assert!(slope(r).abs() < 1e-8);
        }
        assert!(regge_wheeler(pv, 0, 0, 1.0 + 1.0 / 3f64.sqrt()) > 0.0 && regge_wheeler(pv, 0, 0, 1.0 - 1.0 / 3f64.sqrt()) < 0.0);
        // Far out the two agree at the same areal radius R = r e^{1/r} to first order: 5×10⁻⁵ at R = 100, 5×10⁻⁷ at 1000.
        for (areal, gap) in [(100.0f64, 5.1364e-5), (1000.0, 5.0134e-7)] {
            let (mut lo, mut hi) = (1.0, areal);
            for _ in 0..200 {
                let mid = 0.5 * (lo + hi);
                if mid * (1.0 / mid).exp() < areal { lo = mid } else { hi = mid }
            }
            let relative = regge_wheeler(pv, 2, 2, lo) / regge_wheeler(gr, 2, 2, areal) - 1.0;
            assert!(close(-relative, gap, 1e-3), "{areal}: {relative:e}");
        }
    }

    #[test]
    fn the_charged_polarizable_vacuum() {
        // GM/c² = a = 1 m; the charge for b = ratio·a is ratio times the extremal charge.
        let mass = C * C / G;
        let extremal = crate::gem::extremal_charge_to_mass() * mass;
        let k = |ratio: f64, r: f64| charged_pv_dielectric(mass, ratio * extremal, r);
        // Uncharged, the exponential metric; extremal, √K = 1 + a/r, which the cosh and cos forms
        // approach linearly from either side, so their mean is second order.
        assert!(close(k(0.0, 3.0), (2.0f64 / 3.0).exp(), 1e-14) && close(k(1.0, 5.0), 1.44, 1e-12));
        let (below, above) = (k(1.0 - 1e-6, 5.0) - 1.44, k(1.0 + 1e-6, 5.0) - 1.44);
        assert!(below > 1e-8 && above < -1e-8 && (below + above).abs() < 1e-12, "{below:e} {above:e}");
        for ratio in [0.5, 2.0] {
            // Reissner–Nordström in isotropic r, R = r + a + (a² − b²)/(4r): g₀₀ = 1 − 2a/R + b²/R², to O(x³).
            let r = 1e3;
            let areal = r + 1.0 + (1.0 - ratio * ratio) / (4.0 * r);
            let rn = 1.0 - 2.0 / areal + ratio * ratio / (areal * areal);
            assert!((1.0 / k(ratio, r) - rn).abs() < 1e-9, "{ratio}");
            // Puthoff's field equation (64b), √K″ + (2/r)√K′ = (√K′² − b²/r⁴)/√K, at r = 7a.
            let y = |r: f64| k(ratio, r).sqrt();
            let (r, h) = (7.0, 1e-3);
            let (d1, d2) = ((y(r + h) - y(r - h)) / (2.0 * h), (y(r + h) - 2.0 * y(r) + y(r - h)) / (h * h));
            assert!((d2 + 2.0 / r * d1 - (d1 * d1 - ratio * ratio / r.powi(4)) / y(r)).abs() < 1e-9, "{ratio}");
        }
        // Under the extremal ratio K > 1 everywhere; at b = 2a, √K = 0 at √3/(π − arctan √3) =
        // 0.826993a and K = 1 at 1.653987a.
        assert!((1..2000).all(|i| k(0.9, 0.005 * i as f64) > 1.0));
        let zero = 3f64.sqrt() / (PI - 3f64.sqrt().atan());
        assert!((zero - 0.826_993_343).abs() < 1e-9 && k(2.0, zero) < 1e-28);
        let (mut lo, mut hi) = (zero * 1.0001, 100.0);
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if k(2.0, mid) < 1.0 { lo = mid } else { hi = mid }
        }
        assert!((lo - 1.653_986_686).abs() < 1e-8);
    }
}
