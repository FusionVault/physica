//! Tides in a neutron-star inspiral: how much each star deforms in its companion's field, and the
//! one combination of the two that the waveform measures first (at 5PN order).
//!
//! A star of mass `M`, radius `R` and quadrupolar Love number `k₂` has the dimensionless tidal
//! deformability `Λ = (2/3)k₂(c²R/GM)⁵`. A binary's phase carries the mass-weighted
//! `Λ̃ = (16/13)[(m₁ + 12m₂)m₁⁴Λ₁ + (m₂ + 12m₁)m₂⁴Λ₂]/(m₁ + m₂)⁵`, which equals both when the masses
//! are equal ([LVC, PRX 9, 011001](https://doi.org/10.1103/PhysRevX.9.011001), Eq. 5).
//!
//! `k₂` itself comes from the star's structure: [`PiecewisePolytrope::star`] integrates the TOV
//! equations together with Postnikov et al.'s Riccati equation for `y = rH′/H`
//! ([arXiv:1004.5098](https://arxiv.org/abs/1004.5098)), from `y = 2` at the centre, and
//! [`love_number`] turns the surface `y` into `k₂`. The equation of state is one of
//! [Read et al.'s](https://arxiv.org/abs/0812.2163) four-parameter piecewise polytropes, such as SLy
//! (`R₁.₄ = 11.706` km, `k₂ = 0.07617`, `Λ₁.₄ = 295.6`) or AP4 (`Λ₁.₄ = 247.6`). These values are
//! our own, from Read's fits rather than the tabulated equations of state; Hinderer et al.'s
//! tabulated SLy has `k₂ = 0.0763` ([arXiv:0911.3535](https://arxiv.org/abs/0911.3535)).
//!
//! ```
//! use gravitational::tides::{GW170817_LOW_SPIN, PiecewisePolytrope, lambda, lambda_tilde};
//! use metrology::SOLAR_MASS;
//!
//! // Two 1.36 M☉ stars of 12 km with k₂ = 0.1: Λ̃ = 508, inside GW170817's 70–720.
//! let l = lambda(0.1, 12e3, 1.36 * SOLAR_MASS);
//! let tilde = lambda_tilde(1.36, 1.36, l, l);
//! assert!((tilde - l).abs() < 1e-9 && (tilde - 508.0).abs() < 1.0);
//! assert!(tilde > GW170817_LOW_SPIN.0 && tilde < GW170817_LOW_SPIN.1);
//!
//! // A 1.4 M☉ star of SLy: radius 11.706 km, Λ = 295.6.
//! let (m, r, k2) = PiecewisePolytrope::sly().star_with_mass(1.4 * SOLAR_MASS).unwrap();
//! assert!((r / 1e3 - 11.706).abs() < 1e-3 && (lambda(k2, r, m) - 295.6).abs() < 0.1);
//! ```

use core::f64::consts::{FRAC_PI_2, PI};

use metrology::{C, G};

/// GW170817's `Λ̃` with spins below 0.05: the 90% highest-posterior-density interval
/// 300 (+420/−230) (LVC, PRX 9, 011001, Table II, PhenomPNRT, flat prior in `Λ̃`).
pub const GW170817_LOW_SPIN: (f64, f64) = (70.0, 720.0);

/// GW170817's 90% upper limit on `Λ̃` with spins up to 0.89 (the same table).
pub const GW170817_HIGH_SPIN_MAX: f64 = 630.0;

/// The tidal deformability `Λ = (2/3)k₂(c²R/GM)⁵` of a star with Love number `k2`, radius `radius`
/// (m) and mass `mass` (kg).
pub fn lambda(k2: f64, radius: f64, mass: f64) -> f64 {
    2.0 / 3.0 * k2 * (C * C * radius / (G * mass)).powi(5)
}

/// The binary tidal deformability `Λ̃` of stars of masses `m1`, `m2` (any one unit) with
/// deformabilities `l1`, `l2`.
pub fn lambda_tilde(m1: f64, m2: f64, l1: f64, l2: f64) -> f64 {
    16.0 / 13.0 * ((m1 + 12.0 * m2) * m1.powi(4) * l1 + (m2 + 12.0 * m1) * m2.powi(4) * l2) / (m1 + m2).powi(5)
}

/// The quadrupolar Love number `k₂` of a star of compactness `C = GM/(c²R)` whose tidal
/// perturbation has `y = RH′(R)/H(R)` just outside its surface: Hinderer's Eq. 23 as corrected in
/// her erratum ([arXiv:0711.2420](https://arxiv.org/abs/0711.2420)),
/// `k₂ = (8C⁵/5)(1−2C)²u / {2C[6 − 3y + 3C(5y−8)] + 4C³[13 − 11y + C(3y−2) + 2C²(1+y)] + 3(1−2C)²u ln(1−2C)}`
/// with `u = 2 − y + 2C(y−1)`.
///
/// Below `C = 0.02` the braces cancel catastrophically (relative error ~1 at `C ≤ 10⁻⁴`), so the
/// series `k₂ = ((1−2C)²/2)(b₀ + b₁C + b₂C² + b₃C³)` is used, with `b₀ = (2−y)/(y+3)`,
/// `b₁ = (y² + 6y − 6)/(y+3)²`, `b₂ = (y³ + 34y² − 8y + 12)/(7(y+3)³)` and
/// `b₃ = (y⁴ + 62y³ + 84y² + 48y + 36)/(7(y+3)⁴)`: Postnikov et al.'s series (their Eq. 10) with its
/// misprints corrected (they print `−6y` in `b₁`, label `b₃` as the `C⁴` term and give a wrong
/// fifth term). The two agree to 3×10⁻⁸ at the switch (relatively 10⁻⁸–10⁻⁷ for `y ≤ 1`, worse
/// as `y → 2`, where `k₂` itself is of order `C`). At `C → 0` this is the Newtonian
/// `(2 − y)/(2(y + 3))`.
pub fn love_number(compactness: f64, y: f64) -> f64 {
    let (c, z) = (compactness, y + 3.0);
    if c < 0.02 {
        let b1 = (y * y + 6.0 * y - 6.0) / (z * z);
        let b2 = (y.powi(3) + 34.0 * y * y - 8.0 * y + 12.0) / (7.0 * z.powi(3));
        let b3 = (y.powi(4) + 62.0 * y.powi(3) + 84.0 * y * y + 48.0 * y + 36.0) / (7.0 * z.powi(4));
        return (1.0 - 2.0 * c).powi(2) / 2.0 * ((2.0 - y) / z + c * (b1 + c * (b2 + c * b3)));
    }
    let u = 2.0 - y + 2.0 * c * (y - 1.0);
    let braces = 2.0 * c * (6.0 - 3.0 * y + 3.0 * c * (5.0 * y - 8.0))
        + 4.0 * c.powi(3) * (13.0 - 11.0 * y + c * (3.0 * y - 2.0) + 2.0 * c * c * (1.0 + y))
        + 3.0 * (1.0 - 2.0 * c).powi(2) * u * (-2.0 * c).ln_1p();
    8.0 / 5.0 * c.powi(5) * (1.0 - 2.0 * c).powi(2) * u / braces
}

/// A piecewise-polytropic equation of state ([Read, Lackey, Owen & Friedman 2009](https://arxiv.org/abs/0812.2163)):
/// above the density `ρᵢ` the pressure is `p/c² = Kᵢρ^Γᵢ` and the energy density
/// `ε/c² = (1 + aᵢ)ρ + Kᵢρ^Γᵢ/(Γᵢ − 1)`, with `aᵢ` set by continuity of `ε` (so the density has no
/// jumps), and `ρ` the rest-mass density.
#[derive(Clone, Debug)]
pub struct PiecewisePolytrope {
    /// Each piece's `[ρᵢ (kg/m³), Kᵢ (SI, c = 1 form), Γᵢ, aᵢ, specific enthalpy ln((ε + p)/ρc²) at ρᵢ]`.
    pieces: Vec<[f64; 5]>,
}

/// Read et al.'s four-piece SLy crust (their Table II): `(ρᵢ, Kᵢ, Γᵢ)` with `ρᵢ` the lower bound in
/// g/cm³ and `Kᵢ` in the c = 1 form `p/c² = Kρ^Γ` (g/cm³), despite the table caption's dyn/cm².
const SLY_CRUST: [(f64, f64, f64); 4] =
    [(0.0, 6.80110e-9, 1.58425), (2.44034e7, 1.06186e-6, 1.28733), (3.78358e11, 5.32697e1, 0.62223), (2.62780e12, 3.99874e-8, 1.35692)];

/// RK4 steps across the whole star (more in proportion where a piece boundary splits it).
const STEPS: f64 = 2000.0;

impl PiecewisePolytrope {
    /// Read et al.'s parametrisation: the SLy crust, then three core pieces with exponents `gamma`
    /// meeting at 10^14.7 and 10^15 g/cm³, where the pressure is `10^log_p1` dyn/cm² at the first.
    /// The first core piece is extended down to where it meets the crust.
    pub fn read(log_p1: f64, gamma: [f64; 3]) -> Self {
        let p1 = 0.1 * 10f64.powf(log_p1) / (C * C);
        let (rho1, rho2) = (10f64.powf(17.7), 1e18);
        let k1 = p1 / rho1.powf(gamma[0]);
        let k2 = p1 / rho1.powf(gamma[1]);
        // ρ in g/cm³ to kg/m³: K scales by 1000^(1 − Γ).
        let mut pieces: Vec<_> = SLY_CRUST.iter().map(|&(rho, k, g)| (rho * 1e3, k * 1e3f64.powf(1.0 - g), g)).collect();
        let (_, kc, gc) = pieces[3];
        pieces.extend([((kc / k1).powf(1.0 / (gamma[0] - gc)), k1, gamma[0]), (rho1, k2, gamma[1]), (rho2, k2 * rho2.powf(gamma[1] - gamma[2]), gamma[2])]);
        Self::new(&pieces)
    }

    /// SLy (Read et al., Table III: log p₁ = 34.384, Γ = 3.005, 2.988, 2.851).
    pub fn sly() -> Self {
        Self::read(34.384, [3.005, 2.988, 2.851])
    }

    /// AP4 (Read et al., Table III: log p₁ = 34.269, Γ = 2.830, 3.445, 3.348).
    pub fn ap4() -> Self {
        Self::read(34.269, [2.830, 3.445, 3.348])
    }

    /// Pieces `(ρᵢ, Kᵢ, Γᵢ)` in SI, the first starting at `ρ = 0`.
    fn new(pieces: &[(f64, f64, f64)]) -> Self {
        let mut out: Vec<[f64; 5]> = Vec::new();
        for &(rho, k, g) in pieces {
            let a = out.last().map_or(0.0, |&[_, kp, gp, ap, _]| ap + kp * rho.powf(gp - 1.0) / (gp - 1.0) - k * rho.powf(g - 1.0) / (g - 1.0));
            out.push([rho, k, g, a, (a + k * g / (g - 1.0) * rho.powf(g - 1.0)).ln_1p()]);
        }
        Self { pieces: out }
    }

    /// `[ε/c², p/c², (ε + p)dε/dp/c²]` (kg/m³) at specific enthalpy `h` on the piece that holds
    /// enthalpy `inside`, where `ρ` follows in closed form from `e^h = 1 + aᵢ + KᵢΓᵢρ^(Γᵢ−1)/(Γᵢ − 1)`.
    fn at_enthalpy(&self, h: f64, inside: f64) -> [f64; 3] {
        let &[_, k, g, a, _] = self.pieces.iter().rev().find(|p| p[4] <= inside).unwrap_or(&self.pieces[0]);
        let x = (h.exp_m1() - a) * (g - 1.0) / g; // Kρ^(Γ−1)
        if x <= 0.0 {
            return [0.0; 3];
        }
        let rho = (x / k).powf(1.0 / (g - 1.0));
        let (p, e) = (x * rho, (1.0 + a) * rho + x * rho / (g - 1.0));
        [e, p, (e + p) * (1.0 + a + g * x / (g - 1.0)) / (g * x)]
    }

    /// The star of central rest-mass density `central_density` (kg/m³): its mass (kg), radius (m)
    /// and Love number `k₂`. [`lambda`] then gives its tidal deformability.
    pub fn star(&self, central_density: f64) -> (f64, f64, f64) {
        let &[_, k, g, a, _] = self.pieces.iter().rev().find(|p| p[0] <= central_density).unwrap_or(&self.pieces[0]);
        let hc = (a + k * g / (g - 1.0) * central_density.powf(g - 1.0)).ln_1p();
        let geo = G / (C * C);
        let breaks: Vec<f64> = self.pieces.iter().map(|p| p[4]).collect();
        let (m, r, y) = tov_love(hc, &breaks, |h, inside| self.at_enthalpy(h, inside).map(|v| v * geo));
        (m / geo, r, love_number(m / r, y))
    }

    /// The star of mass `mass` (kg) on the stable branch, as in [`star`](Self::star), searched
    /// upward from a central density of 3×10¹⁷ kg/m³ (about nuclear density; about 0.2 M☉ for SLy).
    /// `None` if it is lighter than that star or heavier than the maximum mass.
    pub fn star_with_mass(&self, mass: f64) -> Option<(f64, f64, f64)> {
        stable_root(|x| self.star(x.exp()), 3e17f64.ln(), |s| s.0 - mass)
    }
}

/// The TOV equations with Postnikov's `r y′ + y² + y e^λ[1 + 4πr²(p − ε)] + r²Q = 0`,
/// `Q = 4πe^λ(5ε + 9p + (ε + p)dε/dp) − 6e^λ/r² − ν′²`, in geometric units (`G = c = 1`), from the
/// centre at specific enthalpy `hc` to the surface at `h = 0` (Lindblom 1992). The variable is θ with
/// `h = h_c cos²θ`: `r ∝ sin θ` near the centre, and an `n < 1` polytrope's divergent `dε/dp` at
/// the surface becomes finite. `eos(h, inside)` gives `[ε, p, (ε + p)dε/dp]` at `h` in the form that
/// holds at enthalpy `inside`; `breaks` are the enthalpies where that form changes, which RK4 steps
/// onto, taking each segment's form from its midpoint. Returns `(M, R, y)`, with `y` reduced by
/// `4πR³ε₋/M` for a finite surface density `ε₋` (Damour & Nagar 2009; Hinderer et al. 2010).
/// Read's piecewise polytropes have no density jumps, so only the uniform-density test uses it.
fn tov_love(hc: f64, breaks: &[f64], eos: impl Fn(f64, f64) -> [f64; 3]) -> (f64, f64, f64) {
    // With μ = m/r³ as the second variable, every derivative stays regular at the centre.
    let rhs = |t: f64, [r, mu, y]: [f64; 3], inside: f64| {
        let [e, p, w] = eos(hc * t.cos().powi(2), inside);
        let el = 1.0 / (1.0 - 2.0 * mu * r * r);
        let dr = hc * (2.0 * t).sin() / (el * r * (mu + 4.0 * PI * p));
        let nu = 2.0 * el * (mu + 4.0 * PI * p) * r;
        let q = 4.0 * PI * el * (5.0 * e + 9.0 * p + w) - 6.0 * el / (r * r) - nu * nu;
        [dr, (4.0 * PI * e - 3.0 * mu) * dr / r, -(y * y + y * el * (1.0 + 4.0 * PI * r * r * (p - e)) + r * r * q) / r * dr]
    };
    // Start at θ₀ with r² = 3(h_c − h)/(2π(ε_c + 3p_c)), μ = 4πε_c/3 and y = 2 (errors ~θ₀²).
    let t0: f64 = 1e-6;
    let [e, p, _] = eos(hc, hc);
    let mut u = [(3.0 * hc * t0.sin().powi(2) / (2.0 * PI * (e + 3.0 * p))).sqrt(), 4.0 / 3.0 * PI * e, 2.0];
    let mut edges: Vec<f64> = breaks.iter().filter(|&&b| b > 0.0 && b < hc).map(|b| (b / hc).sqrt().acos()).collect();
    edges.extend([t0, FRAC_PI_2]);
    edges.sort_by(f64::total_cmp);
    for w in edges.windows(2) {
        let n = (STEPS * (w[1] - w[0]) / FRAC_PI_2).ceil().max(8.0);
        let dt = (w[1] - w[0]) / n;
        let inside = hc * (0.5 * (w[0] + w[1])).cos().powi(2);
        for i in 0..n as usize {
            let t = w[0] + i as f64 * dt;
            let add = |u: [f64; 3], d: [f64; 3], f: f64| [u[0] + f * d[0], u[1] + f * d[1], u[2] + f * d[2]];
            let k1 = rhs(t, u, inside);
            let k2 = rhs(t + dt / 2.0, add(u, k1, dt / 2.0), inside);
            let k3 = rhs(t + dt / 2.0, add(u, k2, dt / 2.0), inside);
            let k4 = rhs(t + dt, add(u, k3, dt), inside);
            u = core::array::from_fn(|j| u[j] + dt / 6.0 * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]));
        }
    }
    let [r, mu, y] = u;
    (mu * r.powi(3), r, y - 4.0 * PI * eos(0.0, 0.0)[0] / mu)
}

/// The star `star(x)` where `f` first turns non-negative as `x` rises from `x0` in steps of 0.2, by
/// bisection; `None` if the mass (`.0`) peaks first.
fn stable_root(star: impl Fn(f64) -> (f64, f64, f64), x0: f64, f: impl Fn((f64, f64, f64)) -> f64) -> Option<(f64, f64, f64)> {
    let (mut lo, mut last) = (x0, star(x0));
    if f(last) >= 0.0 {
        return None;
    }
    let mut hi = lo;
    loop {
        hi += 0.2;
        let s = star(hi);
        if f(s) >= 0.0 {
            break;
        }
        if s.0 < last.0 {
            return None;
        }
        (lo, last) = (hi, s);
    }
    for _ in 0..50 {
        let mid = 0.5 * (lo + hi);
        if f(star(mid)) < 0.0 { lo = mid } else { hi = mid }
    }
    Some(star(hi))
}

#[cfg(test)]
mod tests {
    use super::*;
    use metrology::SOLAR_MASS;

    /// `k₂` of the star of compactness `c` on the stable branch of `star`, which returns
    /// `(M, R, y)` in geometric units at central enthalpy `e^x`.
    fn k2_at(c: f64, star: impl Fn(f64) -> (f64, f64, f64)) -> f64 {
        let (m, r, y) = stable_root(star, 1e-6f64.ln(), |(m, r, _)| m / r - c).unwrap();
        love_number(m / r, y)
    }

    /// An "e-polytrope" `p = ε^Γ` (geometric units): `h = (Γ/(Γ−1)) ln(1 + ε^(Γ−1))`.
    fn e_polytrope(gamma: f64) -> impl Fn(f64) -> (f64, f64, f64) {
        move |x| {
            tov_love(x.exp(), &[], |h, _| {
                let e = ((gamma - 1.0) * h / gamma).exp_m1().max(0.0).powf(1.0 / (gamma - 1.0));
                if e == 0.0 { [0.0; 3] } else { [e, e.powf(gamma), (e.powf(2.0 - gamma) + e) / gamma] }
            })
        }
    }

    #[test]
    fn love_numbers_in_closed_form() {
        // Damour & Nagar (2009): the incompressible star at the Buchdahl limit C = 4/9 has y = 2
        // outside, and k₂ = 4096/(10935(308 − 81 ln 3)).
        let dn09 = 4096.0 / (10935.0 * (308.0 - 81.0 * 3f64.ln()));
        assert!((love_number(4.0 / 9.0, 2.0) / dn09 - 1.0).abs() < 1e-14 && (dn09 - 0.001_710_300_596_149_0).abs() < 5e-17);
        // The Newtonian limits: 3/4 for uniform density (y = −1) and (15 − π²)/(2π²) for an n = 1
        // polytrope (y = (π² − 9)/3), both from Postnikov et al.
        assert!((love_number(0.0, -1.0) - 0.75).abs() < 1e-15);
        assert!((love_number(0.0, (PI * PI - 9.0) / 3.0) - (15.0 - PI * PI) / (2.0 * PI * PI)).abs() < 1e-15);
        // The series and Eq. 23 meet at C = 0.02, where Eq. 23 is still good to 10⁻⁹.
        for y in [-1.0, 0.0, 0.29, 1.0, 2.0] {
            assert!((love_number(0.02 - 1e-15, y) - love_number(0.02, y)).abs() < 3e-8, "{y}");
        }
    }

    #[test]
    fn an_incompressible_star() {
        // Uniform ε = 1: h = ln(1 + p), no dε/dp term inside, and the surface jump y → y − 3. The
        // Schwarzschild interior gives p_c = (1 − q)/(3q − 1), q = √(1 − 2C), and the TOV
        // integration must land on that C. k₂ as recomputed in the research (gap5_love_tov.out).
        for (c, want) in [(1e-8f64, 0.75), (0.05, 0.556_257), (0.2, 0.183_019), (0.4, 0.010_347)] {
            let q = (1.0 - 2.0 * c).sqrt();
            let (m, r, y) = tov_love(((1.0 - q) / (3.0 * q - 1.0)).ln_1p(), &[], |h, _| [1.0, h.exp_m1(), 0.0]);
            assert!((m / r / c - 1.0).abs() < 1e-10, "{c}: {}", m / r);
            assert!((love_number(m / r, y) - want).abs() < 6e-7, "{c}: {}", love_number(m / r, y));
        }
    }

    #[test]
    fn polytropes_and_their_two_conventions() {
        // DN09's fit for the Γ = 2 "µ-polytrope" (p = Kρ^Γ in the rest-mass density, Table I):
        // k₂ = k₂ᴺ(1 − 2C)²ΣaₙCⁿ, within 4×10⁻⁴ from C = 0.05 to 0.2.
        let mu = PiecewisePolytrope::new(&[(0.0, 1e-19, 2.0)]);
        let a = [0.9991, -2.9287, -1.1373, 14.0013, -50.9711];
        for c in [0.05f64, 0.1, 0.15, 0.2] {
            let fit = (15.0 - PI * PI) / (2.0 * PI * PI) * (1.0 - 2.0 * c).powi(2) * a.iter().rev().fold(0.0, |s, an| s * c + an);
            let (_, _, k2) = stable_root(|x| mu.star(x.exp()), 1e16f64.ln(), |(m, r, _)| G * m / (C * C * r) - c).unwrap();
            assert!((k2 / fit - 1.0).abs() < 6e-4, "{c}: {k2} against {fit}");
        }
        // Hinderer's polytropes (erratum, Table 1) are "e-polytropes", p = Kε^Γ: her n = 1 rows
        // 0.122, 0.0776 and 0.0459 at C = 0.1, 0.15 and 0.2 agree to 1%, while the µ-polytrope is
        // 5–20% lower.
        for (c, hinderer) in [(0.1, 0.122), (0.15, 0.0776), (0.2, 0.0459)] {
            assert!((k2_at(c, e_polytrope(2.0)) / hinderer - 1.0).abs() < 0.01, "{c}");
        }
        // Her n = 0.5 rows (0.251, 0.173, 0.095) disagree by up to 13%; this enthalpy integration,
        // where the n < 1 surface is regular, matches the research's radial one to 10⁻⁴.
        for (c, ours, hinderer) in [(0.1, 0.22842, 0.251), (0.15, 0.15314, 0.173), (0.2, 0.09676, 0.095)] {
            let k2 = k2_at(c, e_polytrope(3.0));
            assert!((k2 - ours).abs() < 1e-5 && (hinderer / k2 - 1.0).abs() > 0.015, "{c}: {k2}");
        }
    }

    #[test]
    fn read_et_al_piecewise_polytropes() {
        // Read et al.'s fits, Table III: SLy's tabulated R₁.₄ = 11.736 km less the fit's 0.21%, so
        // 11.711 km, and M_max = 2.049 plus 0.02%, 2.0494; AP4's 11.428 km less 0.90%, 11.325 km,
        // and 2.213 less 1.08%, 2.1891. Hinderer et al. (2010, Table I) give k₂ = 0.0763 for
        // tabulated SLy. R, k₂ and Λ₁.₄ are ours: the research script gives the same
        // (11.70603 km, 0.0761665, 295.572 with its 1.98847×10³⁰ kg) once its surface cut at
        // p = 10⁻¹⁴p_c, which dropped the outer 0.95 m and printed 11.705 km and 0.07620, is lowered.
        for (eos, radius, read, k2, our_lambda, m_max) in [
            (PiecewisePolytrope::sly(), 11.706, 11.711, 0.07617, 295.6, 2.0494),
            (PiecewisePolytrope::ap4(), 11.320, 11.325, 0.07545, 247.6, 2.1891),
        ] {
            let (m, r, k) = eos.star_with_mass(1.4 * SOLAR_MASS).unwrap();
            assert!((r / 1e3 - radius).abs() < 5e-4 && (r / 1e3 / read - 1.0).abs() < 1e-3, "{r}");
            assert!((k - k2).abs() < 5e-6, "{k}");
            assert!((lambda(k, r, m) - our_lambda).abs() < 0.05, "{}", lambda(k, r, m));
            // The maximum mass, by golden section in ln ρ_c.
            let (mut a, mut b, g) = (1.5e18f64.ln(), 6e18f64.ln(), (5f64.sqrt() - 1.0) / 2.0);
            for _ in 0..60 {
                let (x1, x2) = (b - g * (b - a), a + g * (b - a));
                if eos.star(x1.exp()).0 > eos.star(x2.exp()).0 { b = x2 } else { a = x1 }
            }
            let most = eos.star(a.exp()).0 / SOLAR_MASS;
            assert!((most / m_max - 1.0).abs() < 1e-3, "{most}");
            assert!(eos.star_with_mass(1.001 * most * SOLAR_MASS).is_none());
        }
        let (_, _, k) = PiecewisePolytrope::sly().star_with_mass(1.4 * SOLAR_MASS).unwrap();
        assert!((k / 0.0763 - 1.0).abs() < 2e-3);
    }

    #[test]
    fn the_injections_of_table_vi() {
        // LVC's injections: (m₁, m₂), (Λ₁, Λ₂) → Λ̃, printed to the nearest unit.
        for (m1, m2, l1, l2, want) in [(1.38, 1.37, 275.0, 309.0, 292.0), (1.68, 1.13, 77.0, 973.0, 303.0), (1.38, 1.37, 1018.0, 1063.0, 1040.0)] {
            let got = lambda_tilde(m1, m2, l1, l2);
            assert!((got - want).abs() < 0.5, "({m1}, {m2}): {got}");
        }
        // Equal masses give Λ̃ = Λ, and the order of the stars does not matter.
        assert!((lambda_tilde(1.4, 1.4, 190.0, 190.0) - 190.0).abs() < 1e-12);
        assert!((lambda_tilde(1.4, 1.2, 190.0, 500.0) - lambda_tilde(1.2, 1.4, 500.0, 190.0)).abs() < 1e-12);
    }

    #[test]
    fn a_star_deforms_as_its_compactness_to_the_minus_fifth() {
        // A Newtonian n = 1 polytrope has k₂ = (15 − π²)/(2π²) = 0.25991; integrating the
        // Clairaut–Radau equation gives the same to 2×10⁻¹² (and 3/4 for uniform density). At
        // compactness C = GM/(c²R) its Λ is (2/3)k₂/C⁵, a hundred times larger at a C 2.5× smaller.
        let k2 = (15.0 - PI * PI) / (2.0 * PI * PI);
        assert!((k2 - 0.25991).abs() < 5e-6);
        let (m, r) = (1.4 * SOLAR_MASS, 12e3);
        let c = G * m / (C * C * r);
        assert!((lambda(k2, r, m) * c.powi(5) - 2.0 / 3.0 * k2).abs() < 1e-15);
        assert!((lambda(k2, 2.5 * r, m) / lambda(k2, r, m) - 2.5f64.powi(5)).abs() < 1e-9);
    }
}
