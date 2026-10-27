//! Incompressible Navier–Stokes and MHD on the flat torus, and the analytic pieces of OpenAI's
//! preprint *Velocity-Field Detection of Computation in Forced Navier–Stokes Flows* (27 September
//! 2026).
//!
//! The equations are `∂ₜu + (u·∇)u = −∇p + ν∆u + f` with `∇·u = 0`, on the torus of period `L`
//! in two or three dimensions ([`Torus`]), and with a magnetic field, incompressible resistive
//! MHD ([`Torus::mhd`]) and the kinematic induction equation ([`Torus::induction`]). The solvers
//! are pseudo-spectral. Derivatives come from [`geopotential::fft::fftn`] and products are taken
//! on the grid. The 2/3 rule removes aliasing, and the Leray projection `û − k(k·û)/|k|²` takes
//! the place of the pressure. What the solvers evolve is therefore the Galerkin truncation: the
//! PDE's exact properties (positivity, maximum principles) hold for it only to the truncation
//! error, while its quadratic invariants are exact. Time steps are fourth-order exponential time
//! differencing (ETDRK4), which carries diffusion exactly and a constant forcing too, so a
//! steadily forced state stays put. Diagnostics on fields given on the grid ([`Torus::advect`],
//! [`Torus::residual`], [`Torus::forcing`], [`Torus::pressure`]) take their products on a grid
//! twice as fine, so they are alias-free for any field below the Nyquist mode. Norms and inner
//! products are box averages, the paper's normalised measure.
//!
//! The preprint's device (its §2, eqs. 5–6) is this. Take a divergence-free planar field
//! `a(t, Y)` and a scalar `w(t, Y)`. The z-independent velocity `U = (a, w)` is divergence-free,
//! and `U·∇ = a·∇_Y`. Let the horizontal force be `F = ∂ₜa + (a·∇)a − ν∆a` ([`Torus::forcing`]),
//! and let `w` solve `∂ₜw + a·∇w = ν∆w + h` ([`Torus::advection_diffusion`]) with `w(0) = 0`.
//! Then `U` with `p = 0` solves Navier–Stokes under the force `(F, h)`, from rest if `a(0) = 0`.
//! The vertical velocity is a passive scalar. The paper steers it with a planar machine, so
//! that `u₃ > ½` somewhere in a fixed band exactly when a Turing machine halts (Theorem 1.1).
//! That event is therefore undecidable (Corollary 6.1), and no program can compute it. What this
//! module holds is the equations, the exact-solution family, the kernels, profiles and bounds
//! the proofs use, and solvers that reproduce them. The construction itself is the example
//! `navier_stokes_computation`. Status: the theorems are **mainstream mathematics (preprint)**,
//! not peer reviewed; the solvers and the exact solutions they are tested on are textbook.
//!
//! ```
//! use core::f64::consts::TAU;
//! use magnetics::navier_stokes::Torus;
//!
//! // Taylor–Green on the 2π torus: u = (sin x cos y, −cos x sin y) decays as e^(−2νt).
//! let torus = Torus { n: [8, 8], length: TAU };
//! let u = [torus.sample(|[x, y]| x.sin() * y.cos()), torus.sample(|[x, y]| -x.cos() * y.sin())];
//! let nu = 0.1;
//! let later = torus.navier_stokes(&u, nu, (0.0, 1.0), 10, |_| [vec![0.0; 64], vec![0.0; 64]]);
//! assert!((later[0][9] / u[0][9] - (-2.0 * nu).exp()).abs() < 1e-13);
//! ```
//!
//! Sources: OpenAI, [*Velocity-Field Detection of Computation in Forced Navier–Stokes Flows*](https://github.com/openai/math/blob/main/preprints/Velocity-Field-Detection-of-Computation-in-Forced-Navier-Stokes-Flows-September-27-2026/manuscript.pdf)
//! (27 September 2026, Apache-2.0): §2 eqs. 5–7, §3.1 eqs. 10–11, Lemmas 3.2–3.3, eq. 21, §5
//! eqs. 28, 35 and 37. The machinery it builds on: Moore,
//! [Phys. Rev. Lett. 64, 2354 (1990)](https://doi.org/10.1103/PhysRevLett.64.2354) and
//! [Nonlinearity 4, 199 (1991)](https://doi.org/10.1088/0951-7715/4/2/002); Bennett,
//! [IBM J. Res. Dev. 17, 525 (1973)](https://doi.org/10.1147/rd.176.0525); Cardona, Miranda &
//! Peralta-Salas, [IMRN 2022, 18092](https://doi.org/10.1093/imrn/rnab233); Dyhr,
//! González-Prieto, Miranda & Peralta-Salas, PNAS Nexus 5, pgag131 (2026). The test solutions:
//! Taylor & Green, [Proc. R. Soc. A 158, 499 (1937)](https://doi.org/10.1098/rspa.1937.0013);
//! Dombre et al., [J. Fluid Mech. 167, 353 (1986)](https://doi.org/10.1017/S0022112086002859)
//! (ABC flows). The numerics: Orszag, [J. Atmos. Sci. 28, 1074 (1971)](https://doi.org/10.1175/1520-0469(1971)028%3C1074:OTEOAI%3E2.0.CO;2)
//! (the 2/3 rule); Cox & Matthews, [J. Comput. Phys. 176, 430 (2002)](https://doi.org/10.1006/jcph.2002.6995)
//! (ETDRK4); Kassam & Trefethen, [SIAM J. Sci. Comput. 26, 1214 (2005)](https://doi.org/10.1137/S1064827502410633)
//! (its coefficients' cancellation, avoided here by series); Bernstein's inequality for the
//! gradient bound. MHD: Walén's exact Alfvén waves (Walén 1944, Ark. Mat. Astron. Fys. 30A, 15);
//! Biskamp, *Magnetohydrodynamic Turbulence* (Cambridge, 2003) for the ideal invariants.

use core::f64::consts::PI;

use geopotential::fft::{fftn, wavenumber};
use geopotential::linalg::symmetric_eigen;

/// The preprint's flat exponential `b(r) = e^(−1/r)` for `r > 0`, and 0 otherwise (§3.1). Every
/// derivative vanishes at 0.
pub fn flat_exp(r: f64) -> f64 {
    if r > 0.0 { (-1.0 / r).exp() } else { 0.0 }
}

/// The smooth step `θ(r) = b(r)/(b(r) + b(1 − r))`: 0 for `r ≤ 0`, 1 for `r ≥ 1`, flat at both
/// ends, with `θ(1 − r) = 1 − θ(r)` (§3.1). The denominator is at least `e^(−2)`.
pub fn smooth_step(r: f64) -> f64 {
    smooth_step_derivatives(r)[0]
}

/// `[θ, θ′, θ″]` at `r`, from `b′ = b/r²` and `b″ = b(1/r⁴ − 2/r³)`.
pub fn smooth_step_derivatives(r: f64) -> [f64; 3] {
    let (p, q) = (flat_exp(r), flat_exp(1.0 - r));
    if p == 0.0 || q == 0.0 {
        return [if p == 0.0 { 0.0 } else { 1.0 }, 0.0, 0.0];
    }
    let s = 1.0 - r;
    let (p1, q1) = (p / (r * r), -q / (s * s));
    let (p2, q2) = (p * (1.0 / r.powi(4) - 2.0 / r.powi(3)), q * (1.0 / s.powi(4) - 2.0 / s.powi(3)));
    let (sum, sum1) = (p + q, p1 + q1);
    let top = p1 * q - p * q1;
    [p / sum, top / (sum * sum), (p2 * q - p * q2) / (sum * sum) - 2.0 * top * sum1 / sum.powi(3)]
}

/// `[A, A′, A″]` of the plateau `A(s) = S(k(c + s))·S(k(c − s))` built from a step `S` given as
/// `[S, S′, S″]`.
fn plateau(step: impl Fn(f64) -> [f64; 3], s: f64, c: f64, k: f64) -> [f64; 3] {
    let ([u, u1, u2], [v, v1, v2]) = (step(k * (c + s)), step(k * (c - s)));
    [u * v, k * (u1 * v - u * v1), k * k * (u2 * v - 2.0 * u1 * v1 + u * v2)]
}

/// The preprint's bump `g(Z) = Πⱼ θ(2(Zⱼ + 1))θ(2(1 − Zⱼ))` (§3.1): in `[0, 1]`, 1 on
/// `|Z|∞ ≤ ½`, 0 outside `[−1, 1]²`. Its mass is `(3/2)² = 2.25` (eq. 11 bounds it by 4), and
/// `max‖∇g‖ = 4`, `max‖D²g‖ ≈ 40.46`, so the paper's unstated constant `C` is 41.
pub fn bump(z: [f64; 2]) -> f64 {
    bump_derivatives(z).0
}

/// The bump `g`, its gradient and its Hessian at `z`.
pub fn bump_derivatives(z: [f64; 2]) -> (f64, [f64; 2], [[f64; 2]; 2]) {
    let [x, y] = z.map(|s| plateau(smooth_step_derivatives, s, 1.0, 2.0));
    (x[0] * y[0], [x[1] * y[0], x[0] * y[1]], [[x[2] * y[0], x[1] * y[1]], [x[1] * y[1], x[0] * y[2]]])
}

/// The heat kernel of the unit two-torus at unit diffusivity, the periodised Gaussian
/// `K(τ, Z) = (4πτ)⁻¹ Σₖ e^(−|Z + k|²/4τ)` of Lemma 3.3. It factors into two circle kernels,
/// each summed over images for `τ < ¼` and as the Fourier series `1 + 2Σ e^(−4π²m²τ) cos 2πmz`
/// above, both to `e^(−40)` of their leading term. The lemma's eq. 14: `K ≤ H` when `τ ≥ 1` or
/// the torus distance of `Z` from 0 is at least `d`, for `(d, H) = (1/32, 10⁹)` or `(1/16, 10⁸)`.
/// Both are loose: the true supremum over `τ ≤ 1` at distance `≥ d` is `1/(πed²)` (120 and 30),
/// and over `τ ≥ 1` it is `K(1, 0) = 1 + 4e^(−4π²)`.
///
/// # Panics
/// If `tau` is not positive.
pub fn torus_heat_kernel(tau: f64, z: [f64; 2]) -> f64 {
    assert!(tau > 0.0, "the heat kernel needs a positive time");
    let circle = |x: f64| {
        let x = x - x.round();
        if tau < 0.25 {
            let k = (160.0 * tau).sqrt().ceil() as i32 + 1;
            (-k..=k).map(|j| (-(x + j as f64).powi(2) / (4.0 * tau)).exp()).sum::<f64>() / (4.0 * PI * tau).sqrt()
        } else {
            let m = (40.0 / (4.0 * PI * PI * tau)).sqrt().ceil() as i32 + 1;
            1.0 + 2.0 * (1..=m).map(|j| (-4.0 * PI * PI * (j * j) as f64 * tau).exp() * (2.0 * PI * j as f64 * x).cos()).sum::<f64>()
        }
    };
    circle(z[0]) * circle(z[1])
}

/// The heat kernel of the plane, `H_s(X) = (4πνs)⁻¹ e^(−|X|²/4νs)` (m⁻²), for diffusivity `nu`
/// (m²/s) after time `s` (s): the kernel of the preprint's mild solution, eq. 35.
pub fn plane_heat_kernel(nu: f64, s: f64, x: [f64; 2]) -> f64 {
    (-(x[0] * x[0] + x[1] * x[1]) / (4.0 * nu * s)).exp() / (4.0 * PI * nu * s)
}

/// `‖∂ₗH_s‖₁ = 1/√(πνs)` (m⁻¹), in closed form: `∂₁H = −X₁H/2νs`, and the mean of `|X₁|` under a
/// Gaussian of variance `2νs` is `√(4νs/π)`. This is the preprint's `C_ν s^(−1/2)` with
/// `C_ν = 1/√(πν)`, the bound that makes the transport term of eq. 35 a contraction over short
/// times.
pub fn plane_heat_gradient_l1(nu: f64, s: f64) -> f64 {
    1.0 / (PI * nu * s).sqrt()
}

/// The preprint's translation gate, eq. 28: `W = ∇⊥[χ(ξ/R)(G₁ξ₂ − G₂ξ₁)]` with `∇⊥ = (∂₂, −∂₁)`
/// and `χ(a) = Πₗ θ(aₗ + 3)θ(3 − aₗ)`, at offset `xi` from the gate's centre, for velocity `g`
/// and radius `r`. It is divergence-free as a perpendicular gradient, equals `G` where
/// `|ξ|∞ ≤ 2R` and vanishes where `|ξ|∞ ≥ 3R`.
pub fn gate(g: [f64; 2], xi: [f64; 2], r: f64) -> [f64; 2] {
    let [x, y] = xi.map(|s| plateau(smooth_step_derivatives, s / r, 3.0, 1.0));
    let (chi, phi) = (x[0] * y[0], g[0] * xi[1] - g[1] * xi[0]);
    [chi * g[0] + phi / r * x[0] * y[1], chi * g[1] - phi / r * x[1] * y[0]]
}

/// `[P, P′, P″]` of the quintic step `P(s) = 10s³ − 15s⁴ + 6s⁵` on `[0, 1]`, 0 to the left and 1
/// to the right: C² (§5.4). `|P″| ≤ 10/√3` at `s = (3 ∓ √3)/6`; the preprint's 360 bounds it.
pub fn quintic_step(s: f64) -> [f64; 3] {
    if s <= 0.0 {
        return [0.0; 3];
    }
    if s >= 1.0 {
        return [1.0, 0.0, 0.0];
    }
    [s.powi(3) * (10.0 - 15.0 * s + 6.0 * s * s), 30.0 * (s * (1.0 - s)).powi(2), 60.0 * s * (1.0 - s) * (1.0 - 2.0 * s)]
}

/// The preprint's moving cutoff `φ_R(a) = A∗(a₁/R)A∗(a₂/R)`, `A∗(s) = P(2(1 + s))P(2(1 − s))`,
/// and its Laplacian (eq. 37): `φ_R` is 1 on `|a|∞ ≤ R/2` and 0 on `|a|∞ ≥ R`, and
/// `‖∆φ_R‖∞ ≤ 2880/R²` in the paper (`43.49/R²` in fact). It is only C² (`P‴ = 60` at the
/// ends), so differentiate it analytically, as here, not spectrally.
pub fn cutoff(a: [f64; 2], r: f64) -> (f64, f64) {
    let [x, y] = a.map(|s| plateau(quintic_step, s / r, 1.0, 2.0));
    (x[0] * y[0], (x[2] * y[0] + x[0] * y[2]) / (r * r))
}

/// The preprint's viscosity rescaling, eq. 21. Given a solution `(a, w)` with force `(F, h)` at
/// unit viscosity, evaluated at time `νt`, it returns the solution `(νa, νa, w)` and force
/// `(ν²F, ν²F, νh)` at viscosity `nu` and time `t`. Every horizontal term gains `ν²` and every
/// vertical one `ν`. This holds for the z-independent family `(a, w)`, where `w` is passive, not
/// for a general 3-D flow (whose vertical velocity would scale as `ν` too). A time derivative maps
/// as a force does: `∂ₜ(νa(νt)) = ν²∂ₛa`, `∂ₜw(νt) = ν∂ₛw`.
pub fn rescale_viscosity(nu: f64, u: [f64; 3], f: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    ([nu * u[0], nu * u[1], u[2]], [nu * nu * f[0], nu * nu * f[1], nu * f[2]])
}

/// A periodic box of `n[0] × … × n[D−1]` points (each a power of two) and period `length` (m) on
/// every axis: the torus `T^D`. Value `i` is at grid index `i` read row-major, the last axis
/// fastest, as [`geopotential::planar::Grid`] stores its values. Fields are `Vec<f64>`, one value
/// per point, and vector fields `[Vec<f64>; D]`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Torus<const D: usize> {
    pub n: [usize; D],
    pub length: f64,
}

/// A spectrum's real and imaginary parts.
type Spectrum = (Vec<f64>, Vec<f64>);

/// A weighted transport term `c(w·g)` for [`combine`]: the weight, the vector field `w` and the
/// gradient `g` of one component.
type Transport<'a, V, const D: usize> = (f64, &'a [V], &'a [Vec<f64>; D]);

impl<const D: usize> Torus<D> {
    /// The number of grid points.
    pub fn size(&self) -> usize {
        self.n.iter().product()
    }

    /// The position (m) of grid point `i`.
    pub fn position(&self, i: usize) -> [f64; D] {
        let index = self.index(i);
        core::array::from_fn(|a| index[a] as f64 * self.length / self.n[a] as f64)
    }

    /// `f` at every grid point.
    pub fn sample(&self, f: impl Fn([f64; D]) -> f64) -> Vec<f64> {
        (0..self.size()).map(|i| f(self.position(i))).collect()
    }

    /// The box average of `v`.
    pub fn mean(&self, v: &[f64]) -> f64 {
        v.iter().sum::<f64>() / v.len() as f64
    }

    /// The box average of `u·v`, the paper's `⟨u, v⟩`.
    pub fn inner(&self, u: &[Vec<f64>; D], v: &[Vec<f64>; D]) -> f64 {
        (0..D).map(|a| u[a].iter().zip(&v[a]).map(|(x, y)| x * y).sum::<f64>()).sum::<f64>() / self.size() as f64
    }

    /// The gradient of `v`, spectrally.
    pub fn gradient(&self, v: &[f64]) -> [Vec<f64>; D] {
        let s = self.forward(v);
        core::array::from_fn(|a| self.inverse(self.derivative(&s, a)))
    }

    /// The Laplacian of `v`, spectrally.
    pub fn laplacian(&self, v: &[f64]) -> Vec<f64> {
        let (mut re, mut im) = self.forward(v);
        for i in 0..re.len() {
            let k2 = self.k2(i);
            (re[i], im[i]) = (-k2 * re[i], -k2 * im[i]);
        }
        self.inverse((re, im))
    }

    /// The divergence of `u`.
    pub fn divergence(&self, u: &[Vec<f64>; D]) -> Vec<f64> {
        let mut out = vec![0.0; self.size()];
        for (a, ua) in u.iter().enumerate() {
            let d = self.inverse(self.derivative(&self.forward(ua), a));
            out.iter_mut().zip(d).for_each(|(o, x)| *o += x);
        }
        out
    }

    /// `u·∇v`, derivatives spectral and the product taken on a grid twice as fine, then cut
    /// back: exactly the true product's modes `2|m| < n`, with no aliasing, for fields below the
    /// Nyquist mode. (A product on this grid would fold modes above `n/2` back onto kept ones.)
    pub fn advect(&self, u: &[Vec<f64>; D], v: &[f64]) -> Vec<f64> {
        let fine = Torus { n: self.n.map(|n| 2 * n), length: self.length };
        let s = self.forward(v);
        let mut product = vec![0.0; fine.size()];
        for (a, ua) in u.iter().enumerate() {
            let ua = fine.inverse(self.resample(&self.forward(ua), &fine));
            let d = fine.inverse(self.resample(&self.derivative(&s, a), &fine));
            product.iter_mut().zip(ua.iter().zip(d)).for_each(|(p, (x, y))| *p += x * y);
        }
        self.inverse(fine.resample(&fine.forward(&product), self))
    }

    /// `v` cut to the modes the 2/3 rule keeps (`3|m| < n` on every axis): the Galerkin
    /// projection `P_J` of Lemma 3.2, which the solvers apply to their data.
    pub fn truncate(&self, v: &[f64]) -> Vec<f64> {
        let s = self.flatten(vec![self.forward(v)]);
        let len = self.size();
        self.inverse((s[..len].to_vec(), s[len..].to_vec()))
    }

    /// The Leray projection of `u` onto divergence-free fields, `û − k(k·û)/|k|²`, which keeps
    /// the mean.
    pub fn project(&self, u: &[Vec<f64>; D]) -> [Vec<f64>; D] {
        let mut s: Vec<Spectrum> = u.iter().map(|v| self.forward(v)).collect();
        self.leray(&mut s);
        let mut s = s.into_iter();
        core::array::from_fn(|_| self.inverse(s.next().expect("one spectrum per component")))
    }

    /// The momentum residual `∂ₜu + (u·∇)u + ∇p − ν∆u − f` and the divergence of `u`, for fields
    /// and a time derivative `dudt` given on the grid. Both vanish, to the spectral accuracy of
    /// the grid, on a solution.
    pub fn residual(&self, u: &[Vec<f64>; D], dudt: &[Vec<f64>; D], p: &[f64], f: &[Vec<f64>; D], nu: f64) -> ([Vec<f64>; D], Vec<f64>) {
        let grad_p = self.gradient(p);
        let r = core::array::from_fn(|a| {
            let (conv, lap) = (self.advect(u, &u[a]), self.laplacian(&u[a]));
            (0..self.size()).map(|i| dudt[a][i] + conv[i] + grad_p[a][i] - nu * lap[i] - f[a][i]).collect()
        });
        (r, self.divergence(u))
    }

    /// The force `F = ∂ₜa + (a·∇)a − ν∆a` that makes a prescribed divergence-free field `a`,
    /// with time derivative `dadt`, a solution with zero pressure: the preprint's eq. 5. With a
    /// vertical source `h` alongside it, `(F, h)` forces `(a, w)` (see the module notes).
    pub fn forcing(&self, a: &[Vec<f64>; D], dadt: &[Vec<f64>; D], nu: f64) -> [Vec<f64>; D] {
        let zero = vec![0.0; self.size()];
        self.residual(a, dadt, &zero, &core::array::from_fn(|_| zero.clone()), nu).0
    }

    /// The pressure (mean zero) that keeps `u` divergence-free under the force `f`: the solution
    /// of `∆p = ∇·(f − (u·∇)u)`.
    pub fn pressure(&self, u: &[Vec<f64>; D], f: &[Vec<f64>; D]) -> Vec<f64> {
        let g: [Vec<f64>; D] = core::array::from_fn(|a| f[a].iter().zip(self.advect(u, &u[a])).map(|(x, c)| x - c).collect());
        self.inverse_laplacian(&self.divergence(&g))
    }

    /// `‖∇u‖²`, the box average of `Σ(∂ᵦuₐ)²`.
    pub fn gradient_norm2(&self, u: &[Vec<f64>; D]) -> f64 {
        u.iter().map(|ua| self.gradient(ua).iter().flatten().map(|x| x * x).sum::<f64>()).sum::<f64>() / self.size() as f64
    }

    /// An upper bound on `‖∇u‖∞,op`, the largest operator (spectral) norm of the velocity
    /// gradient anywhere, the rate in the uniqueness estimate; at most 8/7 of it. Each
    /// `ξ·∇u η` is a trigonometric sum with frequencies along any segment below
    /// `Σₐ kₐ|δₐ|`, so by Bernstein's inequality it falls by at most a fraction `θ²/2` between its
    /// maximum and the nearest point of a grid refined `r`-fold, `θ = Σₐ πmₐ/(rnₐ)` with `mₐ`
    /// the highest mode present. The bound is that grid's maximum over `1 − θ²/2`, with `r` the
    /// least power of two giving `θ ≤ ½`. Modes below 10⁻¹³ of the largest count as absent.
    pub fn gradient_sup(&self, u: &[Vec<f64>; D]) -> f64 {
        let spectra: Vec<Spectrum> = u.iter().map(|v| self.forward(v)).collect();
        let size = |(re, im): &Spectrum, i: usize| re[i].hypot(im[i]);
        let peak = spectra.iter().flat_map(|s| (0..self.size()).map(move |i| size(s, i))).fold(0.0, f64::max);
        let mut top = [0usize; D];
        for s in &spectra {
            for i in (0..self.size()).filter(|&i| size(s, i) > 1e-13 * peak) {
                let index = self.index(i);
                (0..D).for_each(|a| top[a] = top[a].max(index[a].min(self.n[a] - index[a])));
            }
        }
        let theta = |r: usize| (0..D).map(|a| PI * top[a] as f64 / (r * self.n[a]) as f64).sum::<f64>();
        let r = (0..).map(|p| 1usize << p).find(|&r| theta(r) <= 0.5).expect("θ halves with each doubling");
        let fine = Torus { n: self.n.map(|n| r * n), length: self.length };
        let grads: Vec<[Vec<f64>; D]> = spectra.iter().map(|s| core::array::from_fn(|a| fine.inverse(self.resample(&self.derivative(s, a), &fine)))).collect();
        let max = (0..fine.size())
            .map(|i| {
                let jtj: Vec<f64> = (0..D * D).map(|rc| (0..D).map(|a| grads[a][rc / D][i] * grads[a][rc % D][i]).sum()).collect();
                symmetric_eigen(&jtj, D).0.into_iter().fold(0f64, f64::max).sqrt()
            })
            .fold(0.0, f64::max);
        max / (1.0 - theta(r).powi(2) / 2.0)
    }

    /// The uniqueness estimate, eq. 7, for two solutions `u` and `v` under the same force. The
    /// difference `e = v − u` obeys `½ d‖e‖²/dt = −ν‖∇e‖² − ⟨(e·∇)u, e⟩` exactly on the torus,
    /// since `⟨(v·∇)e, e⟩ = 0` and `⟨∇q, e⟩ = 0`. So `½ d‖e‖²/dt + ν‖∇e‖² ≤ ‖∇u‖∞,op‖e‖²`, and
    /// Gronwall from `e(0) = 0` gives `e = 0`. It returns that rate and the right side
    /// `‖∇u‖∞,op‖e‖² − ν‖∇e‖²`, which bounds it.
    pub fn difference_energy_rate(&self, u: &[Vec<f64>; D], v: &[Vec<f64>; D], nu: f64) -> (f64, f64) {
        let e: [Vec<f64>; D] = core::array::from_fn(|a| v[a].iter().zip(&u[a]).map(|(x, y)| x - y).collect());
        let stretch: [Vec<f64>; D] = core::array::from_fn(|a| self.advect(&e, &u[a]));
        let (e2, de2) = (self.inner(&e, &e), self.gradient_norm2(&e));
        (-nu * de2 - self.inner(&stretch, &e), self.gradient_sup(u) * e2 - nu * de2)
    }

    /// Incompressible Navier–Stokes from `u` at `t0` to `t1`, in `steps` ETDRK4 steps at
    /// kinematic viscosity `nu` (m²/s), under the force `force(t)` given on
    /// the grid. The data are first cut to the modes the 2/3 rule keeps (`3|m| < n` on every axis)
    /// and projected onto divergence-free fields. Each step then takes the projected, dealiased
    /// `f − (u·∇)u`. A force's gradient part is balanced by the pressure, which
    /// [`pressure`](Self::pressure) returns. The result is the Galerkin truncation, whose
    /// properties hold to the truncation error. Diffusion is exact, so an explicit step need only
    /// resolve advection: `|u|·k_max·dt` below about 2.8.
    ///
    /// # Panics
    /// If `steps` is zero.
    pub fn navier_stokes(&self, u: &[Vec<f64>; D], nu: f64, (t0, t1): (f64, f64), steps: usize, force: impl Fn(f64) -> [Vec<f64>; D]) -> [Vec<f64>; D] {
        let mut start: Vec<Spectrum> = u.iter().map(|v| self.forward(v)).collect();
        self.leray(&mut start);
        let state = self.flatten(start);
        let tendency = |t: f64, s: &[f64]| {
            let (u, grads) = self.fields_and_gradients(&self.unflatten(s));
            let f = force(t);
            let mut rhs: Vec<Spectrum> = (0..D).map(|a| self.forward(&combine(&[(1.0, &f[a])], &[(-1.0, &u, &grads[a])]))).collect();
            self.leray(&mut rhs);
            self.flatten(rhs)
        };
        let end = self.integrate(state, &[nu; D], (t0, t1), steps, tendency);
        let mut fields = self.unflatten(&end).into_iter().map(|s| self.inverse(s));
        core::array::from_fn(|_| fields.next().expect("one spectrum per component"))
    }

    /// The advection–diffusion equation `∂ₜw + a·∇w = ν∆w + h` from `w` at `t0` to `t1`, with the
    /// velocity and source `fields(t) = (a, h)` given on the grid: Lemma 3.2's Galerkin scheme
    /// `∂ₜw_J = ν∆w_J + P_J(−a·∇w_J + h)`, `P_J` the 2/3-rule cut, in `steps` ETDRK4 steps. `a` should be divergence-free and, with `h`, resolved
    /// by the kept modes. Then the mean of `w` changes only by that of `h`. The PDE's maximum
    /// principle (`w ≥ 0` for non-negative data and source, a sup norm that never grows without a
    /// source, `‖w‖∞ ≤ A(t − t₀)` from zero under `|h| ≤ A`) holds to the truncation error.
    ///
    /// # Panics
    /// If `steps` is zero.
    pub fn advection_diffusion(&self, w: &[f64], nu: f64, (t0, t1): (f64, f64), steps: usize, fields: impl Fn(f64) -> ([Vec<f64>; D], Vec<f64>)) -> Vec<f64> {
        let state = self.flatten(vec![self.forward(w)]);
        let tendency = |t: f64, s: &[f64]| {
            let spectra = self.unflatten(s);
            let grad: [Vec<f64>; D] = core::array::from_fn(|a| self.inverse(self.derivative(&spectra[0], a)));
            let (a, h) = fields(t);
            self.flatten(vec![self.forward(&combine(&[(1.0, &h)], &[(-1.0, &a, &grad)]))])
        };
        let end = self.integrate(state, &[nu], (t0, t1), steps, tendency);
        self.inverse(self.unflatten(&end).remove(0))
    }

    /// Incompressible resistive MHD from `(u, b)` at `t0` to `t1`:
    /// `∂ₜu + (u·∇)u = −∇P + (b·∇)b + ν∆u + f` and `∂ₜb + (u·∇)b = (b·∇)u + η∆b`, with
    /// `∇·u = ∇·b = 0`, in `steps` ETDRK4 steps, `nu` the kinematic viscosity and `eta` the
    /// magnetic diffusivity (m²/s; [`crate::mhd::magnetic_diffusivity`] gives it from the
    /// conductivity). The field is in Alfvén units, `b = B/√(μ₀ρ)` (m/s, as
    /// [`crate::mhd::alfven_speed`] gives its size), so that `(b·∇)b = J×B/ρ + ∇(b²/2)` and `P` is
    /// the total pressure over the density. A uniform `B₀` is the mean of `b`. As in
    /// [`navier_stokes`](Self::navier_stokes), the data are cut by the 2/3 rule and projected,
    /// each tendency is dealiased and projected, and the result is the Galerkin truncation: its
    /// properties hold to the truncation error. Alfvén waves travel at `|b|`, so an explicit step
    /// needs `(|u| + |b|)·k_max·dt` below about 2.8. With `b = 0` it is exactly
    /// [`navier_stokes`](Self::navier_stokes).
    ///
    /// # Panics
    /// If `steps` is zero.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn mhd(&self, u: &[Vec<f64>; D], b: &[Vec<f64>; D], nu: f64, eta: f64, (t0, t1): (f64, f64), steps: usize, force: impl Fn(f64) -> [Vec<f64>; D]) -> ([Vec<f64>; D], [Vec<f64>; D]) {
        let mut start: Vec<Spectrum> = u.iter().chain(b).map(|v| self.forward(v)).collect();
        self.leray(&mut start[..D]);
        self.leray(&mut start[D..]);
        let state = self.flatten(start);
        let tendency = |t: f64, s: &[f64]| {
            let (fields, grads) = self.fields_and_gradients(&self.unflatten(s));
            let (u, b) = fields.split_at(D);
            let f = force(t);
            let mut rhs: Vec<Spectrum> = (0..D)
                .map(|a| self.forward(&combine(&[(1.0, &f[a])], &[(-1.0, u, &grads[a]), (1.0, b, &grads[D + a])])))
                .chain((0..D).map(|a| self.forward(&combine::<_, D>(&[], &[(1.0, b, &grads[a]), (-1.0, u, &grads[D + a])]))))
                .collect();
            self.leray(&mut rhs[..D]);
            self.leray(&mut rhs[D..]);
            self.flatten(rhs)
        };
        let diffusivity: Vec<f64> = [nu; D].into_iter().chain([eta; D]).collect();
        let end = self.integrate(state, &diffusivity, (t0, t1), steps, tendency);
        let mut fields = self.unflatten(&end).into_iter().map(|s| self.inverse(s));
        let mut next = || fields.next().expect("one spectrum per component");
        (core::array::from_fn(|_| next()), core::array::from_fn(|_| next()))
    }

    /// The kinematic induction equation `∂ₜb = ∇×(u×b) + η∆b = (b·∇)u − (u·∇)b + η∆b` from `b` at
    /// `t0` to `t1`, for a divergence-free velocity `velocity(t)` given on the grid and resolved by
    /// the kept modes: the field's growth or decay in a flow it does not act back on, the dynamo
    /// problem. Steps and Galerkin truncation as in [`mhd`](Self::mhd).
    ///
    /// # Panics
    /// If `steps` is zero.
    pub fn induction(&self, b: &[Vec<f64>; D], eta: f64, (t0, t1): (f64, f64), steps: usize, velocity: impl Fn(f64) -> [Vec<f64>; D]) -> [Vec<f64>; D] {
        let mut start: Vec<Spectrum> = b.iter().map(|v| self.forward(v)).collect();
        self.leray(&mut start);
        let state = self.flatten(start);
        let tendency = |t: f64, s: &[f64]| {
            let (b, grad_b) = self.fields_and_gradients(&self.unflatten(s));
            let u = velocity(t);
            let grad_u: Vec<[Vec<f64>; D]> = u.iter().map(|ua| self.gradient(ua)).collect();
            let mut rhs: Vec<Spectrum> = (0..D).map(|a| self.forward(&combine::<_, D>(&[], &[(1.0, &b, &grad_u[a]), (-1.0, &u, &grad_b[a])]))).collect();
            self.leray(&mut rhs);
            self.flatten(rhs)
        };
        let end = self.integrate(state, &[eta; D], (t0, t1), steps, tendency);
        let mut fields = self.unflatten(&end).into_iter().map(|s| self.inverse(s));
        core::array::from_fn(|_| fields.next().expect("one spectrum per component"))
    }

    /// Grid point `i`'s index on each axis.
    fn index(&self, mut i: usize) -> [usize; D] {
        let mut out = [0; D];
        for a in (0..D).rev() {
            out[a] = i % self.n[a];
            i /= self.n[a];
        }
        out
    }

    /// Grid point `i`'s index along `axis`.
    fn axis_index(&self, i: usize, axis: usize) -> usize {
        i / self.n[axis + 1..].iter().product::<usize>() % self.n[axis]
    }

    /// The wavenumber (rad/m) of mode `i` along `axis`.
    fn k(&self, i: usize, axis: usize) -> f64 {
        let n = self.n[axis];
        wavenumber(self.axis_index(i, axis), n, self.length / n as f64)
    }

    /// `|k|²` of mode `i`.
    fn k2(&self, i: usize) -> f64 {
        (0..D).map(|a| self.k(i, a).powi(2)).sum()
    }

    /// Whether the 2/3 rule keeps mode `i`: `3|m| < n` on every axis.
    fn kept(&self, i: usize) -> bool {
        let index = self.index(i);
        (0..D).all(|a| 3 * index[a].min(self.n[a] - index[a]) < self.n[a])
    }

    fn forward(&self, v: &[f64]) -> Spectrum {
        assert_eq!(v.len(), self.size(), "one value per grid point");
        let (mut re, mut im) = (v.to_vec(), vec![0.0; v.len()]);
        fftn(&mut re, &mut im, &self.n, false);
        (re, im)
    }

    /// The real part of the inverse transform.
    fn inverse(&self, (mut re, mut im): Spectrum) -> Vec<f64> {
        fftn(&mut re, &mut im, &self.n, true);
        re
    }

    /// The spectrum times `ik` along `axis`, with the Nyquist mode, whose derivative of a real
    /// field is not real, set to zero.
    fn derivative(&self, (re, im): &Spectrum, axis: usize) -> Spectrum {
        let (n, stride) = (self.n[axis], self.n[axis + 1..].iter().product::<usize>());
        let mut out = (vec![0.0; re.len()], vec![0.0; re.len()]);
        for i in 0..re.len() {
            let j = i / stride % n;
            if 2 * j != n {
                let k = wavenumber(j, n, self.length / n as f64);
                (out.0[i], out.1[i]) = (-k * im[i], k * re[i]);
            }
        }
        out
    }

    /// The Leray projection of a vector field's spectra, in place.
    fn leray(&self, s: &mut [Spectrum]) {
        for i in 0..self.size() {
            let k2 = self.k2(i);
            if k2 == 0.0 {
                continue;
            }
            let k: [f64; D] = core::array::from_fn(|a| self.k(i, a));
            let (kr, ki) = (0..D).fold((0.0, 0.0), |(r, m), a| (r + k[a] * s[a].0[i], m + k[a] * s[a].1[i]));
            for a in 0..D {
                s[a].0[i] -= k[a] * kr / k2;
                s[a].1[i] -= k[a] * ki / k2;
            }
        }
    }

    /// Spectra, cut by the 2/3 rule, laid end to end as `[re₀, im₀, re₁, im₁, …]`.
    fn flatten(&self, spectra: Vec<Spectrum>) -> Vec<f64> {
        let mask: Vec<f64> = (0..self.size()).map(|i| if self.kept(i) { 1.0 } else { 0.0 }).collect();
        spectra.into_iter().flat_map(|(re, im)| re.into_iter().chain(im)).zip(mask.iter().cycle()).map(|(x, m)| x * m).collect()
    }

    /// Spectra laid end to end by [`flatten`](Self::flatten), apart again.
    fn unflatten(&self, s: &[f64]) -> Vec<Spectrum> {
        s.chunks(2 * self.size()).map(|c| (c[..self.size()].to_vec(), c[self.size()..].to_vec())).collect()
    }

    /// The fields of `spectra` on the grid, and each one's gradient.
    fn fields_and_gradients(&self, spectra: &[Spectrum]) -> (Vec<Vec<f64>>, Vec<[Vec<f64>; D]>) {
        let fields = spectra.iter().map(|s| self.inverse(s.clone())).collect();
        let grads = spectra.iter().map(|s| core::array::from_fn(|a| self.inverse(self.derivative(s, a)))).collect();
        (fields, grads)
    }

    /// The spectrum of `v`, a field on this grid, on the grid `to` of the same period: every mode
    /// with `2|m| < n` on both grids, the Nyquist modes dropped. Up, it is the trigonometric
    /// interpolant; down, the cut to the coarse grid's modes.
    fn resample(&self, (re, im): &Spectrum, to: &Torus<D>) -> Spectrum {
        let scale = to.size() as f64 / self.size() as f64;
        let mut out = (vec![0.0; to.size()], vec![0.0; to.size()]);
        'mode: for i in 0..self.size() {
            let (index, mut j) = (self.index(i), 0);
            for (a, (&k, &n)) in index.iter().zip(&self.n).enumerate() {
                let m = if 2 * k < n { k as isize } else { k as isize - n as isize };
                if 2 * m.unsigned_abs() >= n.min(to.n[a]) {
                    continue 'mode;
                }
                j = j * to.n[a] + m.rem_euclid(to.n[a] as isize) as usize;
            }
            (out.0[j], out.1[j]) = (scale * re[i], scale * im[i]);
        }
        out
    }

    /// The mean-zero solution of `∆φ = v`.
    fn inverse_laplacian(&self, v: &[f64]) -> Vec<f64> {
        let (mut re, mut im) = self.forward(v);
        for i in 0..re.len() {
            let k2 = self.k2(i);
            let scale = if k2 > 0.0 { -1.0 / k2 } else { 0.0 };
            (re[i], im[i]) = (scale * re[i], scale * im[i]);
        }
        self.inverse((re, im))
    }

    /// Cox & Matthews' fourth-order exponential time differencing (ETDRK4) for
    /// `ŝ′ = −D|k|²ŝ + N(t, ŝ)` on flattened spectra, with one diffusivity `D` per spectrum. The
    /// linear part is exact, and so is a constant `N`, so a steadily forced state stays put.
    fn integrate(&self, mut s: Vec<f64>, diffusivity: &[f64], (t0, t1): (f64, f64), steps: usize, tendency: impl Fn(f64, &[f64]) -> Vec<f64>) -> Vec<f64> {
        assert!(steps > 0, "at least one step");
        let (len, dt) = (self.size(), (t1 - t0) / steps as f64);
        assert_eq!(s.len(), 2 * len * diffusivity.len(), "one diffusivity per spectrum");
        // Per entry: e^z, e^(z/2), (dt/2)φ₁(z/2) and dt·(f₁, f₂, f₃), z = −D|k|²dt.
        let w: Vec<[f64; 6]> = (0..s.len())
            .map(|i| {
                let z = -diffusivity[i / (2 * len)] * self.k2(i % len) * dt;
                let ([f1, f2, f3, _], half) = (etd_coefficients(z), etd_coefficients(z / 2.0)[3]);
                [z.exp(), (z / 2.0).exp(), dt / 2.0 * half, dt * f1, dt * f2, dt * f3]
            })
            .collect();
        for step in 0..steps {
            let t = t0 + step as f64 * dt;
            let n1 = tendency(t, &s);
            let a: Vec<f64> = (0..s.len()).map(|i| w[i][1] * s[i] + w[i][2] * n1[i]).collect();
            let n2 = tendency(t + 0.5 * dt, &a);
            let b: Vec<f64> = (0..s.len()).map(|i| w[i][1] * s[i] + w[i][2] * n2[i]).collect();
            let n3 = tendency(t + 0.5 * dt, &b);
            let c: Vec<f64> = (0..s.len()).map(|i| w[i][1] * a[i] + w[i][2] * (2.0 * n3[i] - n1[i])).collect();
            let n4 = tendency(t + dt, &c);
            s = (0..s.len()).map(|i| w[i][0] * s[i] + w[i][3] * n1[i] + 2.0 * w[i][4] * (n2[i] + n3[i]) + w[i][5] * n4[i]).collect();
        }
        s
    }
}

/// `Σ cᵢvᵢ + Σ cⱼ(wⱼ·gⱼ)` pointwise on the grid: weighted fields plus weighted transport terms,
/// each `w·g = Σₐ wₐgₐ` for a vector field `w` and the gradient `g` of one component.
fn combine<V: AsRef<[f64]>, const D: usize>(fields: &[(f64, &Vec<f64>)], transports: &[Transport<'_, V, D>]) -> Vec<f64> {
    let len = fields.first().map(|f| f.1.len()).unwrap_or_else(|| transports[0].2[0].len());
    (0..len)
        .map(|i| {
            let plain: f64 = fields.iter().map(|(c, v)| c * v[i]).sum();
            plain + transports.iter().map(|(c, w, g)| c * (0..D).map(|a| w[a].as_ref()[i] * g[a][i]).sum::<f64>()).sum::<f64>()
        })
        .collect()
}

/// ETDRK4's weights `[f₁, f₂, f₃, φ₁]` at `z ≤ 0`: `f₁ = (−4 − z + e^z(4 − 3z + z²))/z³`,
/// `f₂ = (2 + z + e^z(z − 2))/z³`, `f₃ = (−4 − 3z − z² + e^z(4 − z))/z³` and `φ₁ = (e^z − 1)/z`.
/// The closed forms cancel catastrophically near 0 (Kassam & Trefethen's point), so below
/// `|z| = 2` they come from their series `Σ zʲ(j + 1)²/(j + 3)!`, `Σ zʲ(j + 1)/(j + 3)!`,
/// `Σ zʲ(1 − j)/(j + 3)!` and `Σ zʲ/(j + 1)!`. Both branches hold to 6×10⁻¹⁶ of `φ₁` against
/// 60-digit arithmetic for `−10⁴ ≤ z ≤ −10⁻⁸`.
fn etd_coefficients(z: f64) -> [f64; 4] {
    if z.abs() >= 2.0 {
        let e = z.exp();
        let z3 = z * z * z;
        return [(-4.0 - z + e * (4.0 - 3.0 * z + z * z)) / z3, (2.0 + z + e * (z - 2.0)) / z3, (-4.0 - 3.0 * z - z * z + e * (4.0 - z)) / z3, z.exp_m1() / z];
    }
    let (mut out, mut power, mut fact3, mut fact1) = ([0.0; 4], 1.0, 6.0, 1.0);
    for j in 0..40 {
        let jf = j as f64;
        out[0] += power * (jf + 1.0).powi(2) / fact3;
        out[1] += power * (jf + 1.0) / fact3;
        out[2] += power * (1.0 - jf) / fact3;
        out[3] += power / fact1;
        power *= z;
        fact3 *= jf + 4.0;
        fact1 *= jf + 2.0;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::f64::consts::{E, TAU};

    fn close(a: f64, b: f64, rel: f64) -> bool {
        (a / b - 1.0).abs() < rel
    }

    fn max_abs(v: &[f64]) -> f64 {
        v.iter().fold(0.0, |m, x| m.max(x.abs()))
    }

    fn max_diff(a: &[f64], b: &[f64]) -> f64 {
        a.iter().zip(b).fold(0.0, |m, (x, y)| m.max((x - y).abs()))
    }

    #[test]
    fn smooth_step_and_bump() {
        // θ runs from 0 to 1 and is flat at both ends: θ(0.01) = e^(1/0.99 − 100) = 1.1×10⁻⁴³.
        assert_eq!([smooth_step(-1.0), smooth_step(0.0), smooth_step(1.0), smooth_step(2.0), smooth_step(0.5)], [0.0, 0.0, 1.0, 1.0, 0.5]);
        assert!(close(smooth_step(0.01), (1.0 / 0.99 - 100.0f64).exp(), 1e-13) && 1.0 - smooth_step(0.99) < 2e-43);
        // θ(1 − r) = 1 − θ(r), the derivatives match differences, and b(r) + b(1 − r) ≥ e⁻²
        // (its least value is 2e⁻² at r = ½).
        let h = 1e-5;
        for r in (1..100).map(|i| i as f64 / 100.0) {
            let [t, t1, t2] = smooth_step_derivatives(r);
            assert!((t + smooth_step(1.0 - r) - 1.0).abs() < 1e-15);
            assert!((t1 - (smooth_step(r + h) - smooth_step(r - h)) / (2.0 * h)).abs() < 1e-8 * (1.0 + t1.abs()));
            let d2 = (smooth_step_derivatives(r + h)[1] - smooth_step_derivatives(r - h)[1]) / (2.0 * h);
            assert!((t2 - d2).abs() < 1e-6 * (1.0 + t2.abs()), "{r}: {t2} {d2}");
            assert!(flat_exp(r) + flat_exp(1.0 - r) >= 2.0 / (E * E) - 1e-15);
        }
        // The bump is 1 on |Z|∞ ≤ ½, 0 outside [−1, 1]², between them in between.
        assert_eq!([bump([0.0, 0.0]), bump([0.5, -0.5]), bump([1.0, 0.2]), bump([0.3, -1.5])], [1.0, 1.0, 0.0, 0.0]);
        for z in [[0.7, 0.1], [-0.8, 0.6], [0.55, -0.9]] {
            let (g, grad, hess) = bump_derivatives(z);
            assert!(g > 0.0 && g < 1.0);
            let dx = (bump([z[0] + h, z[1]]) - bump([z[0] - h, z[1]])) / (2.0 * h);
            let dxy = (bump_derivatives([z[0], z[1] + h]).1[0] - bump_derivatives([z[0], z[1] - h]).1[0]) / (2.0 * h);
            assert!((grad[0] - dx).abs() < 1e-8 && (hess[0][1] - dxy).abs() < 1e-7 && hess[0][1] == hess[1][0]);
        }
        // Its steepest slope is 4 = 2θ′(½), and its mass is (3/2)² = 2.25, where the paper's
        // eq. 11 bounds it by 4 (per unit width²). The midpoint rule is spectral on a smooth
        // compactly supported function.
        assert!((bump_derivatives([-0.75, 0.0]).1[0] - 4.0).abs() < 1e-14 && smooth_step_derivatives(0.5)[1] == 2.0);
        let line: f64 = (0..400).map(|i| bump([-1.0 + (i as f64 + 0.5) / 200.0, 0.0])).sum::<f64>() / 200.0;
        assert!((line * line - 2.25).abs() < 1e-14, "{}", line * line);
    }

    #[test]
    fn torus_heat_kernel_mass_semigroup_and_eq_14() {
        // Unit mass by the midpoint rule (spectrally accurate for a smooth periodic function), on
        // both sides of the switch from images to Fourier series.
        let m = 64;
        let cell = |i: usize| (i as f64 + 0.5) / m as f64;
        let integral = |f: &dyn Fn([f64; 2]) -> f64| (0..m * m).map(|i| f([cell(i / m), cell(i % m)])).sum::<f64>() / (m * m) as f64;
        for tau in [0.003, 0.2, 0.3, 2.0] {
            assert!((integral(&|z| torus_heat_kernel(tau, z)) - 1.0).abs() < 1e-13, "{tau}");
        }
        // Semigroup: K(τ₁ + τ₂) = K(τ₁) ∗ K(τ₂).
        for z in [[0.1, 0.4], [0.5, 0.5]] {
            let conv = integral(&|y| torus_heat_kernel(0.01, [z[0] - y[0], z[1] - y[1]]) * torus_heat_kernel(0.25, y));
            assert!(close(conv, torus_heat_kernel(0.26, z), 1e-12));
        }
        // Near 0 at short times it is the plane kernel at unit diffusivity.
        assert!(close(torus_heat_kernel(1e-3, [0.01, 0.02]), plane_heat_kernel(1.0, 1e-3, [0.01, 0.02]), 1e-14));
        // Eq. 14. Over τ ≤ 1 at torus distance ≥ d the largest value is the main Gaussian's peak
        // 1/(πed²) at τ = d²/4, 120 at d = 1/32 and 30 at d = 1/16, against H = 10⁹ and 10⁸ and
        // the proof's 8×10⁴/(4πd²) = 6.5×10⁶ and 1.6×10⁶; for τ ≥ 1 it is K(1, 0) = 1 + 4e^(−4π²).
        for (d, big_h) in [(1.0 / 32.0, 1e9), (1.0 / 16.0, 1e8)] {
            let worst = (0..2000)
                .map(|i| d * d * 10f64.powf(-2.0 + 4.0 * i as f64 / 2000.0))
                .filter(|&tau| tau <= 1.0)
                .flat_map(|tau| [torus_heat_kernel(tau, [d, 0.0]), torus_heat_kernel(tau, [d / 2f64.sqrt(), d / 2f64.sqrt()])])
                .fold(0.0, f64::max);
            assert!(close(worst, 1.0 / (PI * E * d * d), 1e-5) && worst < big_h && worst < 8e4 / (4.0 * PI * d * d), "{d}: {worst}");
        }
        // (Below 10⁻¹⁶ from 1, so equal to it in f64.)
        assert!(close(torus_heat_kernel(1.0, [0.0, 0.0]), 1.0 + 4.0 * (-4.0 * PI * PI).exp(), 1e-15));
        assert!(torus_heat_kernel(5.0, [0.0, 0.0]) <= torus_heat_kernel(1.0, [0.0, 0.0]));
        // The proof's shell sum Σₖ e^(−|Z+k|²/8) = 8πK(2, Z) peaks at 8π(1 + 2e^(−8π²))² ≈ 25.1 at
        // Z = 0, against its bound 1 + 8·1056 < 10⁴.
        assert!(close(8.0 * PI * torus_heat_kernel(2.0, [0.0, 0.0]), 8.0 * PI, 1e-15));
    }

    #[test]
    fn plane_heat_kernel_derivative_norm() {
        // ‖∂₁H_s‖₁ by Gauss–Legendre quadrature of a central difference, X₁ over [0, 10σ] (twice,
        // by symmetry, so the kink at 0 is an end) and X₂ over ±10σ, against 1/√(πνs) = C_ν s^(−1/2).
        let (nu, s) = (0.7f64, 0.3);
        let sigma = (2.0 * nu * s).sqrt();
        let (x, w) = geopotential::gauss_legendre(80);
        let eps = 1e-6;
        let mut sum = 0.0;
        for (xi, wi) in x.iter().zip(&w) {
            for (yj, wj) in x.iter().zip(&w) {
                let p = [5.0 * sigma * (xi + 1.0), 10.0 * sigma * yj];
                let d = (plane_heat_kernel(nu, s, [p[0] + eps, p[1]]) - plane_heat_kernel(nu, s, [p[0] - eps, p[1]])) / (2.0 * eps);
                sum += 2.0 * wi * wj * 5.0 * sigma * 10.0 * sigma * d.abs();
            }
        }
        assert!(close(sum, plane_heat_gradient_l1(nu, s), 1e-9), "{sum}");
        assert!(close(plane_heat_gradient_l1(nu, 4.0 * s), plane_heat_gradient_l1(nu, s) / 2.0, 1e-15));
    }

    #[test]
    fn translation_gate() {
        // Divergence-free (by differences), G on |ξ|∞ ≤ 2R, zero beyond 3R.
        let (g, r, h) = ([0.3, -1.7], 5.0, 1e-4);
        for xi in [[11.0, 3.0], [-12.0, 13.5], [2.5, -14.0], [0.0, 12.0], [-13.0, -13.0]] {
            let div = (gate(g, [xi[0] + h, xi[1]], r)[0] - gate(g, [xi[0] - h, xi[1]], r)[0] + gate(g, [xi[0], xi[1] + h], r)[1] - gate(g, [xi[0], xi[1] - h], r)[1]) / (2.0 * h);
            assert!(div.abs() < 1e-8, "{xi:?}: {div}");
        }
        assert!(gate(g, [11.0, 3.0], r) != g && gate(g, [11.0, 3.0], r) != [0.0, 0.0]);
        for xi in [[0.0, 0.0], [10.0, -10.0], [-7.0, 9.9]] {
            assert_eq!(gate(g, xi, r), g);
        }
        for xi in [[15.0, 0.0], [3.0, -15.0], [40.0, 40.0]] {
            assert_eq!(gate(g, xi, r), [0.0, 0.0]);
        }
    }

    #[test]
    fn moving_cutoff_eq_37() {
        // P'' peaks at 10/√3 (paper: ≤ 360), A∗'' at 40/√3 (paper: ≤ 1440); R²‖∆φ_R‖∞ reaches
        // 43.49, against the paper's 2880 ≤ C∗ = 3000.
        let scan = |f: &dyn Fn(f64) -> f64| (0..=20000).map(|i| f(-1.2 + 2.4 * i as f64 / 20000.0).abs()).fold(0.0, f64::max);
        assert!(close(scan(&|s| quintic_step(s)[2]), 10.0 / 3f64.sqrt(), 1e-6));
        assert!(close(scan(&|s| plateau(quintic_step, s, 1.0, 2.0)[2]), 40.0 / 3f64.sqrt(), 1e-6));
        let r = 3.0;
        // R²∆φ_R = A∗″(s₁)A∗(s₂) + A∗(s₁)A∗″(s₂), scanned on a 4001² grid of (s₁, s₂).
        let profile: Vec<[f64; 3]> = (0..=4000).map(|i| plateau(quintic_step, -1.0 + i as f64 / 2000.0, 1.0, 2.0)).collect();
        let lap = profile.iter().flat_map(|p| profile.iter().map(move |q| (p[2] * q[0] + p[0] * q[2]).abs())).fold(0.0, f64::max);
        assert!((lap - 43.49).abs() < 0.01 && lap <= 2880.0, "{lap}");
        assert!(close(cutoff([2.3, 0.4], r).1 * r * r, plateau(quintic_step, 2.3 / r, 1.0, 2.0)[2] * plateau(quintic_step, 0.4 / r, 1.0, 2.0)[0], 1e-15));
        // 1 on |a|∞ ≤ R/2, 0 on |a|∞ ≥ R, in [0, 1], and the Laplacian matches differences.
        assert_eq!([cutoff([1.5, -1.5], r).0, cutoff([3.0, 0.0], r).0, cutoff([0.2, -3.5], r).0], [1.0, 0.0, 0.0]);
        let h = 1e-4;
        for a in [[2.0, 0.3], [-2.2, 2.6], [0.1, -1.9]] {
            let (phi, l) = cutoff(a, r);
            let fd = (cutoff([a[0] + h, a[1]], r).0 + cutoff([a[0] - h, a[1]], r).0 + cutoff([a[0], a[1] + h], r).0 + cutoff([a[0], a[1] - h], r).0 - 4.0 * phi) / (h * h);
            assert!((0.0..=1.0).contains(&phi) && (l - fd).abs() < 1e-5, "{a:?}: {l} {fd}");
        }
    }

    /// The 2π torus with `n` points a side in two dimensions.
    fn plane(n: usize) -> Torus<2> {
        Torus { n: [n, n], length: TAU }
    }

    #[test]
    fn taylor_green_decay() {
        // u = (sin x cos y, −cos x sin y) on the 2π torus (period 2π, not the paper's 1) decays
        // as e^(−2νt): its nonlinear term is a gradient, which the projection removes, and the
        // exponential step carries diffusion exactly. To 4×10⁻¹⁶ at 16², ν = 0.1, t = 2.
        let torus = plane(16);
        let tg = |t: f64| {
            let decay = (-0.2 * t).exp();
            [torus.sample(|[x, y]| x.sin() * y.cos() * decay), torus.sample(|[x, y]| -x.cos() * y.sin() * decay)]
        };
        let zero = || [vec![0.0; 256], vec![0.0; 256]];
        let u = torus.navier_stokes(&tg(0.0), 0.1, (0.0, 2.0), 8, |_| zero());
        let exact = tg(2.0);
        assert!(max_diff(&u[0], &exact[0]) < 1e-15 && max_diff(&u[1], &exact[1]) < 1e-15);
        // The pressure is p = (cos 2x + cos 2y)/4 times e^(−4νt).
        let p = torus.pressure(&exact, &zero());
        let decay = (-0.4f64 * 2.0).exp();
        assert!(max_diff(&p, &torus.sample(|[x, y]| ((2.0 * x).cos() + (2.0 * y).cos()) / 4.0 * decay)) < 1e-15);
    }

    /// The ABC (Arnold–Beltrami–Childress) flow with A = B = C = 1 on the 2π torus: ∇×u = u.
    fn abc(torus: &Torus<3>) -> [Vec<f64>; 3] {
        [torus.sample(|[_, y, z]| z.sin() + y.cos()), torus.sample(|[x, _, z]| x.sin() + z.cos()), torus.sample(|[x, y, _]| y.sin() + x.cos())]
    }

    #[test]
    fn forced_abc_flow_is_steady() {
        // Beltrami: ω × u = 0, so (u·∇)u = ∇|u|²/2, a pressure. Forced by f = νu = −ν∆u it stays
        // put; unforced it decays as e^(−νt). 8³ points, ν = 0.2, t = 5.
        let torus = Torus { n: [8, 8, 8], length: TAU };
        let (u0, nu) = (abc(&torus), 0.2);
        let force = u0.clone().map(|c| c.iter().map(|x| nu * x).collect::<Vec<_>>());
        // ETDRK4 integrates a constant forcing exactly, so the flow stays put to round-off even
        // in 5 steps of 1. (Fourth-order Runge–Kutta in the integrating factor drifted by
        // 2.7×10⁻⁹ in 20 steps.)
        for steps in [5, 20] {
            let steady = torus.navier_stokes(&u0, nu, (0.0, 5.0), steps, |_| force.clone());
            let drift = (0..3).map(|a| max_diff(&steady[a], &u0[a])).fold(0.0, f64::max);
            assert!(drift < 2e-15, "{steps} steps: {drift:e}");
        }
        let free = torus.navier_stokes(&u0, nu, (0.0, 5.0), 20, |_| core::array::from_fn(|_| vec![0.0; 512]));
        assert!((0..3).all(|a| max_diff(&free[a], &u0[a].iter().map(|x| x * (-1.0f64).exp()).collect::<Vec<_>>()) < 1e-13));
        // Its pressure is −|u|²/2 less its mean.
        let p = torus.pressure(&u0, &force);
        let ke: Vec<f64> = (0..512).map(|i| -(0..3).map(|a| u0[a][i] * u0[a][i]).sum::<f64>() / 2.0).collect();
        let mean = torus.mean(&ke);
        assert!(max_diff(&p, &ke.iter().map(|x| x - mean).collect::<Vec<_>>()) < 1e-14);
    }

    #[test]
    fn energy_identity() {
        // ½ d‖u‖²/dt = −ν‖∇u‖² + ⟨f, u⟩, on a genuinely nonlinear 3-D flow at 8³: the solver's
        // energy by central differences (h = 10⁻⁴) against the right side, to 10⁻⁸.
        let torus = Torus { n: [8, 8, 8], length: TAU };
        let u0 = torus.project(&[torus.sample(|[x, y, z]| (y + 2.0 * z).sin() + (2.0 * x).cos()), torus.sample(|[x, y, z]| (x - z).cos() + 0.5 * (2.0 * y).sin()), torus.sample(|[x, y, z]| (x + y).sin() * z.cos())]);
        let f = [torus.sample(|[_, y, _]| (2.0 * y).sin()), torus.sample(|[x, _, z]| (x + z).cos()), torus.sample(|[x, _, _]| x.sin())];
        let (nu, h) = (0.05, 1e-4);
        let run = |u: &[Vec<f64>; 3], t: (f64, f64), steps| torus.navier_stokes(u, nu, t, steps, |_| f.clone());
        let before = run(&u0, (0.0, 0.5 - h), 50);
        let now = run(&before, (0.5 - h, 0.5), 1);
        let after = run(&now, (0.5, 0.5 + h), 1);
        let rate = (torus.inner(&after, &after) - torus.inner(&before, &before)) / (4.0 * h);
        let budget = -nu * torus.gradient_norm2(&now) + torus.inner(&f, &now);
        assert!((rate - budget).abs() < 1e-8 * budget.abs(), "{rate} {budget}");
        // And the flow is not trivial: the nonlinear term moves energy between modes.
        assert!(torus.gradient_norm2(&now) > 0.1);
    }

    #[test]
    fn uniqueness_energy_estimate_eq_7() {
        // Two forced 2-D flows from data 10⁻³ apart: ½ d‖e‖²/dt (central differences) is the
        // identity's −ν‖∇e‖² − ⟨(e·∇)U, e⟩ and lies below ‖∇U‖∞‖e‖² − ν‖∇e‖².
        let torus = plane(16);
        let (nu, h) = (0.05, 1e-3);
        let u0 = torus.project(&[torus.sample(|[x, y]| (x + y).sin() + (2.0 * y).cos()), torus.sample(|[x, y]| (x - 2.0 * y).cos())]);
        let bump = torus.project(&[torus.sample(|[x, y]| (2.0 * x + y).cos()), torus.sample(|[x, _]| (3.0 * x).sin())]);
        let v0: [Vec<f64>; 2] = core::array::from_fn(|a| u0[a].iter().zip(&bump[a]).map(|(x, y)| x + 1e-3 * y).collect());
        let f = [torus.sample(|[_, y]| y.sin()), torus.sample(|[x, _]| (2.0 * x).cos())];
        let run = |u: &[Vec<f64>; 2], t: (f64, f64), steps| torus.navier_stokes(u, nu, t, steps, |_| f.clone());
        let (mut u, mut v, mut t) = (u0, v0, 0.0);
        for _ in 0..3 {
            let (ub, vb) = (run(&u, (t, t + 0.3 - h), 15), run(&v, (t, t + 0.3 - h), 15));
            let (un, vn) = (run(&ub, (t + 0.3 - h, t + 0.3), 1), run(&vb, (t + 0.3 - h, t + 0.3), 1));
            let (ua, va) = (run(&un, (t + 0.3, t + 0.3 + h), 1), run(&vn, (t + 0.3, t + 0.3 + h), 1));
            let e2 = |u: &[Vec<f64>; 2], v: &[Vec<f64>; 2]| (0..2).map(|a| u[a].iter().zip(&v[a]).map(|(x, y)| (x - y).powi(2)).sum::<f64>()).sum::<f64>() / 256.0;
            let rate = (e2(&ua, &va) - e2(&ub, &vb)) / (4.0 * h);
            let (identity, bound) = torus.difference_energy_rate(&un, &vn, nu);
            assert!((rate - identity).abs() < 1e-6 * identity.abs() && identity <= bound, "{rate:e} {identity:e} {bound:e}");
            (u, v, t) = (un, vn, t + 0.3);
        }
    }

    /// `v` on a 2-D grid repeated `nz` times along a third axis: a z-independent field on T³.
    fn lift(v: &[f64], nz: usize) -> Vec<f64> {
        v.iter().flat_map(|&x| core::iter::repeat_n(x, nz)).collect()
    }

    /// The manufactured pair of the reduction test on the 2π torus, at time `t` and viscosity
    /// `nu`: `a = sin t (cos y, sin x)`, `w = (1 − e^(−t))(1 + cos x cos y)`, and the source `h`
    /// that makes `w` solve the scalar equation, as `(a, ∂ₜa, w, ∂ₜw, h)`.
    #[allow(clippy::type_complexity)]
    fn manufactured(torus: &Torus<2>, t: f64, nu: f64) -> ([Vec<f64>; 2], [Vec<f64>; 2], Vec<f64>, Vec<f64>, Vec<f64>) {
        let a = [torus.sample(|[_, y]| t.sin() * y.cos()), torus.sample(|[x, _]| t.sin() * x.sin())];
        let dadt = [torus.sample(|[_, y]| t.cos() * y.cos()), torus.sample(|[x, _]| t.cos() * x.sin())];
        let w = torus.sample(|[x, y]| (1.0 - (-t).exp()) * (1.0 + x.cos() * y.cos()));
        let dwdt = torus.sample(|[x, y]| (-t).exp() * (1.0 + x.cos() * y.cos()));
        let (adv, lap) = (torus.advect(&a, &w), torus.laplacian(&w));
        let h = (0..w.len()).map(|i| dwdt[i] + adv[i] - nu * lap[i]).collect();
        (a, dadt, w, dwdt, h)
    }

    #[test]
    fn reduction_eqs_5_6() {
        // A planar a and source h; w from the scalar solver; then (a, w) as a 3-D flow on
        // 16 × 16 × 2 with the force (F, h) of eq. 5.
        let (t2, nu) = (plane(16), 0.3);
        let t3 = Torus { n: [16, 16, 2], length: TAU };
        let up = |v: &[f64]| lift(v, 2);
        // The scalar solver from w = 0 reproduces the manufactured w to 4.1×10⁻¹⁰ in steps of
        // 0.02, and 16 times better in steps of 0.01: the error is the stepper's fourth order.
        let scalar = |steps| {
            t2.advection_diffusion(&vec![0.0; 256], nu, (0.0, 1.0), steps, |t| {
                let (a, _, _, _, h) = manufactured(&t2, t, nu);
                (a, h)
            })
        };
        let (a, dadt, w_exact, dwdt, h) = manufactured(&t2, 1.0, nu);
        let (coarse, fine) = (max_diff(&scalar(50), &w_exact), max_diff(&scalar(100), &w_exact));
        assert!(coarse < 5e-10 && (14.0..18.0).contains(&(coarse / fine)), "{coarse:e} {fine:e}");
        // The full 3-D residual of U = (a, w), p = 0, under f = (F, h) is round-off, and so is the
        // pressure that the force and flow imply.
        let f = |a: &[Vec<f64>; 2], dadt: &[Vec<f64>; 2], h: &[f64]| {
            let big_f = t2.forcing(a, dadt, nu);
            [up(&big_f[0]), up(&big_f[1]), up(h)]
        };
        let u3 = [up(&a[0]), up(&a[1]), up(&w_exact)];
        let force = f(&a, &dadt, &h);
        let (r, div) = t3.residual(&u3, &[up(&dadt[0]), up(&dadt[1]), up(&dwdt)], &vec![0.0; 512], &force, nu);
        assert!(r.iter().all(|c| max_abs(c) < 1e-13) && max_abs(&div) < 1e-14);
        assert!(max_abs(&t3.pressure(&u3, &force)) < 1e-14);
        // And the 3-D solver itself, from rest, follows (a, w) with no pressure gradient to undo:
        // to 8.5×10⁻¹⁰ in steps of 0.02 (its w is carried by its own a, which has its own step
        // error), and again 16 times better in steps of 0.01.
        let full = |steps| {
            let u = t3.navier_stokes(&core::array::from_fn(|_| vec![0.0; 512]), nu, (0.0, 1.0), steps, |t| {
                let (a, dadt, _, _, h) = manufactured(&t2, t, nu);
                f(&a, &dadt, &h)
            });
            (0..3).map(|c| max_diff(&u[c], &u3[c])).fold(0.0, f64::max)
        };
        let (coarse, fine) = (full(50), full(100));
        assert!(coarse < 1e-9 && (14.0..18.0).contains(&(coarse / fine)), "{coarse:e} {fine:e}");
    }

    #[test]
    fn viscosity_rescaling_eq_21() {
        // The unit-viscosity pair rescaled to ν = 0.3 at t = 2 (unit time νt = 0.6) solves
        // Navier–Stokes at ν = 0.3 with the rescaled force.
        let (t2, nu, t) = (plane(16), 0.3, 2.0);
        let t3 = Torus { n: [16, 16, 2], length: TAU };
        let (a, dadt, w, dwdt, h) = manufactured(&t2, nu * t, 1.0);
        let big_f = t2.forcing(&a, &dadt, 1.0);
        let (mut u, mut f, mut dudt) = (core::array::from_fn(|_| vec![0.0; 256]), core::array::from_fn(|_| vec![0.0; 256]), core::array::from_fn(|_| vec![0.0; 256]));
        for i in 0..256 {
            let (ui, fi) = rescale_viscosity(nu, [a[0][i], a[1][i], w[i]], [big_f[0][i], big_f[1][i], h[i]]);
            let di = rescale_viscosity(nu, [0.0; 3], [dadt[0][i], dadt[1][i], dwdt[i]]).1;
            for c in 0..3 {
                (u[c][i], f[c][i], dudt[c][i]) = (ui[c], fi[c], di[c]);
            }
        }
        let up = |v: &[Vec<f64>; 3]| v.clone().map(|c| lift(&c, 2));
        let (r, _) = t3.residual(&up(&u), &up(&dudt), &vec![0.0; 512], &up(&f), nu);
        assert!(r.iter().all(|c| max_abs(c) < 1e-14));
        // At the wrong viscosity it does not.
        let (r, _) = t3.residual(&up(&u), &up(&dudt), &vec![0.0; 512], &up(&f), 1.0);
        assert!(max_abs(&r[2]) > 0.1);
    }

    /// A divergence-free cellular flow on the unit torus, `∇⊥ψ` for `ψ = cos(t) sin 2πx sin 2πy/2π`.
    fn cells(torus: &Torus<2>, t: f64) -> [Vec<f64>; 2] {
        let k = TAU;
        [torus.sample(|[x, y]| t.cos() * (k * x).sin() * (k * y).cos()), torus.sample(|[x, y]| -t.cos() * (k * x).cos() * (k * y).sin())]
    }

    #[test]
    fn lemma_3_2_maximum_principle_and_mass() {
        let torus = Torus { n: [32, 32], length: 1.0 };
        let (nu, k) = (0.02, TAU);
        let zero = vec![0.0; 1024];
        // A source in [0, A], A = 1, whose mean is ¼: from zero, w ≥ 0, ‖w‖∞ ≤ At and the mean
        // is exactly t/4.
        let source = torus.sample(|[x, y]| (1.0 + (k * x).cos()) * (1.0 + (k * y).cos()) / 4.0);
        let mut w = zero.clone();
        for step in 0..4 {
            let (t0, t1) = (step as f64 * 0.25, (step + 1) as f64 * 0.25);
            w = torus.advection_diffusion(&w, nu, (t0, t1), 10, |t| (cells(&torus, t), source.clone()));
            let min = w.iter().fold(f64::MAX, |m, &x| m.min(x));
            assert!(min >= 0.0 && max_abs(&w) <= t1, "{min} {}", max_abs(&w));
            assert!((torus.mean(&w) - t1 / 4.0).abs() < 1e-15);
        }
        // Without a source the sup norm never grows and the mean stays put.
        let w0 = torus.sample(|[x, y]| ((k * x).cos() + (k * y).sin()).exp());
        let (mut w, mut sup, mean) = (w0.clone(), max_abs(&w0), torus.mean(&w0));
        for step in 0..6 {
            w = torus.advection_diffusion(&w, nu, (step as f64 * 0.2, (step + 1) as f64 * 0.2), 10, |t| (cells(&torus, t), zero.clone()));
            assert!(max_abs(&w) <= sup + 1e-12 && (torus.mean(&w) - mean).abs() < 1e-14);
            sup = max_abs(&w);
        }
        assert!(sup < 0.9 * max_abs(&w0));
    }

    #[test]
    fn advection_diffusion_converges() {
        // Fourth order in time: halving the step cuts the error ~16 times. Spectral in space:
        // 16² is 10³ times worse than 32² against a 64² reference, on the points they share.
        let (nu, k) = (0.05, TAU);
        let w0 = |torus: &Torus<2>| torus.sample(|[x, y]| (0.5 * ((k * x).sin() + (k * y).cos())).exp());
        let solve = |n: usize, steps: usize| {
            let torus = Torus { n: [n, n], length: 1.0 };
            let zero = vec![0.0; n * n];
            torus.advection_diffusion(&w0(&torus), nu, (0.0, 0.5), steps, |t| (cells(&torus, t), zero.clone()))
        };
        let reference = solve(64, 40);
        let at = |w: &[f64], n: usize| -> f64 { (0..n * n).map(|i| (w[i] - reference[(i / n) * (64 / n) * 64 + (i % n) * (64 / n)]).abs()).fold(0.0, f64::max) };
        let (e16, e32) = (at(&solve(16, 40), 16), at(&solve(32, 40), 32));
        assert!(e16 > 1e3 * e32 && e32 < 1e-6, "{e16:e} {e32:e}");
        let (coarse, fine) = (solve(32, 5), solve(32, 10));
        let exact = solve(32, 80);
        let ratio = max_diff(&coarse, &exact) / max_diff(&fine, &exact);
        assert!((10.0..25.0).contains(&ratio), "{ratio}");
    }
    #[test]
    fn etd_weights_against_60_digit_values() {
        // [f₁, f₂, f₃, φ₁] from the closed forms in 60-digit arithmetic, on both sides of the
        // series' hand-over at |z| = 2; at 0 they are 1/6 and 1.
        assert_eq!(etd_coefficients(0.0), [1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 1.0]);
        for (z, exact) in [
            (-1e-3, [0.16650007497778274, 0.16658335832777876, 0.16666665833611052, 0.9995001666250083]),
            (-1.0, [0.056964470628461425, 0.10363832351432696, 0.16060279414278839, 0.6321205588285577]),
            (-1.999, [0.013189593548782956, 0.06769490085782291, 0.14851170036357098, 0.4324808973436456]),
            (-2.0, [0.01316325433592779, 0.06766764161830635, 0.14849853757254047, 0.43233235838169365]),
            (-30.0, [-0.0009629629629664079, 0.0010370370370371479, 0.03014814814814803, 0.03333333333333022]),
        ] {
            let got = etd_coefficients(z);
            assert!((0..4).all(|k| (got[k] - exact[k]).abs() < 1e-15 * exact[3]), "{z}: {got:?}");
            // A constant forcing is integrated exactly: f₁ + 4f₂ + f₃ = φ₁.
            assert!((got[0] + 4.0 * got[1] + got[2] - got[3]).abs() < 2e-16 * exact[3] * 4.0);
        }
    }

    #[test]
    fn padded_advection_is_alias_free() {
        // Fields filling the 8² grid's band (modes up to 3): u·∇v has modes up to 6, which a
        // product on the grid folds back onto kept modes. The padded product equals the true
        // product's 8² modes, here from a 32² grid that holds it exactly.
        let (coarse, fine) = (plane(8), plane(32));
        let u = |t: &Torus<2>| [t.sample(|[x, y]| (3.0 * y).sin() + (x + 2.0 * y).cos()), t.sample(|[x, y]| (3.0 * x).cos() - (2.0 * x - y).sin())];
        let v = |t: &Torus<2>| t.sample(|[x, y]| (3.0 * x + y).sin() + (2.0 * y).cos());
        let (uf, gf) = (u(&fine), fine.gradient(&v(&fine)));
        let exact_fine: Vec<f64> = (0..fine.size()).map(|i| uf[0][i] * gf[0][i] + uf[1][i] * gf[1][i]).collect();
        let exact = coarse.inverse(fine.resample(&fine.forward(&exact_fine), &coarse));
        let (uc, vc) = (u(&coarse), v(&coarse));
        assert!(max_diff(&coarse.advect(&uc, &vc), &exact) < 1e-14);
        let g = coarse.gradient(&vc);
        let naive: Vec<f64> = (0..coarse.size()).map(|i| uc[0][i] * g[0][i] + uc[1][i] * g[1][i]).collect();
        assert!(max_diff(&naive, &exact) > 1.0, "the grid product aliases");
    }

    #[test]
    fn gradient_sup_bounds_the_maximum_between_points() {
        // u = (0, sin(3x + 0.3)) on 8²: ‖∇u‖ = 3|cos(3x + 0.3)| peaks at 3 between the points,
        // where the grid sees only 3 cos 0.3 = 2.866. The bound lies in [3, 3·8/7].
        let torus = plane(8);
        let u = [vec![0.0; 64], torus.sample(|[x, _]| (3.0 * x + 0.3).sin())];
        let on_grid = torus.gradient(&u[1])[0].iter().fold(0f64, |m, x| m.max(x.abs()));
        let bound = torus.gradient_sup(&u);
        assert!((on_grid - 3.0 * 0.3f64.cos()).abs() < 1e-14 && (3.0..=3.0 * 8.0 / 7.0).contains(&bound), "{on_grid} {bound}");
        // And in three dimensions, u = (sin(2z + 0.4), 0, 0) on 8³: the bound is in [2, 2·8/7].
        let t3 = Torus { n: [8, 8, 8], length: TAU };
        let bound = t3.gradient_sup(&[t3.sample(|[_, _, z]| (2.0 * z + 0.4).sin()), vec![0.0; 512], vec![0.0; 512]]);
        assert!((2.0..=2.0 * 8.0 / 7.0).contains(&bound), "{bound}");
    }

    /// `∇⊥ψ = (∂ᵧψ, −∂ₓψ)` on a 2-D grid: a divergence-free field from its stream function.
    fn perp(torus: &Torus<2>, psi: &[f64]) -> [Vec<f64>; 2] {
        let [gx, gy] = torus.gradient(psi);
        [gy, gx.iter().map(|x| -x).collect()]
    }

    fn zero_force<const D: usize>(torus: &Torus<D>) -> impl Fn(f64) -> [Vec<f64>; D] + '_ {
        |_| core::array::from_fn(|_| vec![0.0; torus.size()])
    }

    #[test]
    fn mhd_without_a_field_is_navier_stokes() {
        let torus = plane(16);
        let u = perp(&torus, &torus.sample(|[x, y]| (x + 2.0 * y).sin() + 0.5 * (3.0 * x).cos() * y.sin()));
        let zero = [vec![0.0; 256], vec![0.0; 256]];
        let (um, bm) = torus.mhd(&u, &zero, 0.05, 0.1, (0.0, 1.0), 20, zero_force(&torus));
        let un = torus.navier_stokes(&u, 0.05, (0.0, 1.0), 20, zero_force(&torus));
        assert!((0..2).all(|a| max_diff(&um[a], &un[a]) == 0.0 && max_abs(&bm[a]) == 0.0));
    }

    #[test]
    fn alfven_waves_of_any_amplitude() {
        // About a uniform B₀, u = −b′ with ν = η is an exact nonlinear solution: the Elsasser
        // field z⁺ = u + b = B₀ is uniform, so the nonlinear terms cancel, and z⁻ = u − b is
        // carried along B₀ undistorted, each mode decaying as e^(−νk²t) (Walén). Here b′ is of
        // order 4 against |B₀| = 1.1, far from linear.
        let (torus, nu, b0) = (plane(16), 0.05, [1.0, 0.5]);
        let psi = |t: f64| {
            move |[x, y]: [f64; 2]| {
                let (x, y) = (x - b0[0] * t, y - b0[1] * t);
                0.8 * (-5.0 * nu * t).exp() * (x + 2.0 * y).sin() + 0.6 * (-10.0 * nu * t).exp() * (3.0 * x - y).cos()
            }
        };
        let wave = |t: f64| perp(&torus, &torus.sample(psi(t)));
        let b = |t: f64| wave(t).map(|c| c.iter().map(|x| x + 0.0).collect::<Vec<_>>());
        let with_mean = |w: [Vec<f64>; 2]| [w[0].iter().map(|x| x + b0[0]).collect::<Vec<_>>(), w[1].iter().map(|x| x + b0[1]).collect()];
        let u0 = wave(0.0).map(|c| c.iter().map(|x| -x).collect::<Vec<_>>());
        let run = |steps| torus.mhd(&u0, &with_mean(b(0.0)), nu, nu, (0.0, 1.0), steps, zero_force(&torus));
        let exact = with_mean(b(1.0));
        let error = |steps| {
            let (u, b) = run(steps);
            // z⁺ = B₀ to round-off at every step count: the cancellation is exact.
            assert!((0..2).all(|a| u[a].iter().zip(&b[a]).all(|(x, y)| (x + y - b0[a]).abs() < 1e-13)));
            (0..2).map(|a| max_diff(&b[a], &exact[a])).fold(0.0, f64::max)
        };
        let (coarse, fine) = (error(50), error(100));
        assert!(coarse < 1e-6 && (14.0..18.0).contains(&(coarse / fine)), "{coarse:e} {fine:e}");
    }

    #[test]
    fn free_magnetic_decay() {
        // With u = 0 and b = ∇⊥ψ, ψ a single-|k| shell (k² = 2), the Lorentz force is a gradient:
        // the flow stays at rest and b decays as e^(−ηk²t), to round-off.
        let (torus, eta) = (plane(16), 0.1);
        let psi = torus.sample(|[x, y]| x.sin() * y.sin() + (x + y).cos());
        let b0 = perp(&torus, &psi);
        let zero = [vec![0.0; 256], vec![0.0; 256]];
        let (u, b) = torus.mhd(&zero, &b0, 0.05, eta, (0.0, 2.0), 10, zero_force(&torus));
        let decay = (-2.0 * eta * 2.0).exp();
        assert!((0..2).all(|a| max_abs(&u[a]) < 1e-14 && max_diff(&b[a], &b0[a].iter().map(|x| x * decay).collect::<Vec<_>>()) < 1e-14));
    }

    #[test]
    fn beltrami_pair_decays_in_place() {
        // u = αA and b = βA with A the ABC flow (∇×A = A): the Lorentz and inertial forces are
        // both gradients and ∇×(u×b) = 0, so with ν = η each decays as e^(−νt), to round-off.
        let torus = Torus { n: [8, 8, 8], length: TAU };
        let (a, nu) = (abc(&torus), 0.1);
        let scaled = |c: f64| a.clone().map(|v| v.iter().map(|x| c * x).collect::<Vec<_>>());
        let (u, b) = torus.mhd(&scaled(0.5), &scaled(0.8), nu, nu, (0.0, 3.0), 12, zero_force(&torus));
        let decay = (-nu * 3.0f64).exp();
        let (ue, be) = (scaled(0.5 * decay), scaled(0.8 * decay));
        assert!((0..3).all(|c| max_diff(&u[c], &ue[c]) < 1e-14 && max_diff(&b[c], &be[c]) < 1e-14));
    }

    /// A nonlinear 3-D pair `(u, b)` on 8³ with modes up to 2, both divergence-free.
    fn turbulent_pair(torus: &Torus<3>) -> ([Vec<f64>; 3], [Vec<f64>; 3]) {
        let u = torus.project(&[torus.sample(|[x, y, z]| (y + 2.0 * z).sin() + (2.0 * x).cos()), torus.sample(|[x, y, z]| (x - z).cos() + 0.5 * (2.0 * y).sin()), torus.sample(|[x, y, z]| (x + y).sin() * z.cos())]);
        let b = torus.project(&[torus.sample(|[x, y, z]| 0.7 * (x + z).cos() + 0.3 * (2.0 * y).sin()), torus.sample(|[x, _, z]| (2.0 * z - x).sin()), torus.sample(|[x, y, _]| 0.6 * (x - 2.0 * y).cos() + 0.2)]);
        (u, b)
    }

    #[test]
    fn mhd_energy_identity_and_solenoidal_field() {
        // ½ d(‖u‖² + ‖b‖²)/dt = −ν‖∇u‖² − η‖∇b‖² + ⟨f, u⟩ (the Lorentz work cancels the
        // induction's), by central differences (h = 10⁻⁴) against the right side, to 10⁻⁸. ∇·b
        // stays at round-off.
        let torus = Torus { n: [8, 8, 8], length: TAU };
        let (u0, b0) = turbulent_pair(&torus);
        let f = [torus.sample(|[_, y, _]| (2.0 * y).sin()), torus.sample(|[x, _, z]| (x + z).cos()), torus.sample(|[x, _, _]| x.sin())];
        let (nu, eta, h) = (0.05, 0.08, 1e-4);
        let run = |(u, b): &([Vec<f64>; 3], [Vec<f64>; 3]), t: (f64, f64), steps| torus.mhd(u, b, nu, eta, t, steps, |_| f.clone());
        let before = run(&(u0, b0), (0.0, 0.5 - h), 50);
        let now = run(&before, (0.5 - h, 0.5), 1);
        let after = run(&now, (0.5, 0.5 + h), 1);
        let energy = |(u, b): &([Vec<f64>; 3], [Vec<f64>; 3])| torus.inner(u, u) + torus.inner(b, b);
        let rate = (energy(&after) - energy(&before)) / (4.0 * h);
        let budget = -nu * torus.gradient_norm2(&now.0) - eta * torus.gradient_norm2(&now.1) + torus.inner(&f, &now.0);
        assert!((rate - budget).abs() < 1e-8 * budget.abs(), "{rate} {budget}");
        assert!(max_abs(&torus.divergence(&now.1)) < 1e-13);
        // The field does work on the flow: the Lorentz term ⟨u, (b·∇)b⟩ is not small.
        let lorentz: [Vec<f64>; 3] = core::array::from_fn(|a| torus.advect(&now.1, &now.1[a]));
        assert!(torus.inner(&now.0, &lorentz).abs() > 0.05);
    }

    /// `∇×v` on T³.
    fn curl(torus: &Torus<3>, v: &[Vec<f64>; 3]) -> [Vec<f64>; 3] {
        let g: Vec<[Vec<f64>; 3]> = v.iter().map(|c| torus.gradient(c)).collect();
        core::array::from_fn(|a| {
            let (j, k) = ((a + 1) % 3, (a + 2) % 3);
            g[k][j].iter().zip(&g[j][k]).map(|(x, y)| x - y).collect()
        })
    }

    #[test]
    fn ideal_invariants() {
        // With ν = η = 0 the Galerkin system conserves the energy ½⟨u² + b²⟩, the cross-helicity
        // ⟨u·b⟩ and the magnetic helicity ⟨a·b⟩ (∇×a = b, ∇·a = 0) in 3-D, and ⟨A²⟩
        // (b = ∇⊥A) in 2-D, exactly; the steps keep them to the stepper's order or better (the
        // drift falls about 32-fold per halving here, fifth order).
        let t3 = Torus { n: [8, 8, 8], length: TAU };
        let (u0, b) = turbulent_pair(&t3);
        // A Beltrami part gives b magnetic helicity, and a part along u cross-helicity. The pair's
        // uniform part is removed: with a mean field B₀, ⟨a·b⟩ is not conserved on the torus
        // (B₀ has no periodic potential), and it drifts by 2% here.
        let a = abc(&t3);
        let b0: [Vec<f64>; 3] = core::array::from_fn(|c| {
            let mean = t3.mean(&b[c]);
            (0..512).map(|i| b[c][i] - mean + 0.4 * a[c][i] + 0.5 * u0[c][i]).collect()
        });
        let helicity = |b: &[Vec<f64>; 3]| {
            let a: [Vec<f64>; 3] = curl(&t3, b).map(|c| t3.inverse_laplacian(&c).iter().map(|x| -x).collect());
            t3.inner(&a, b)
        };
        let invariants3 = |(u, b): &([Vec<f64>; 3], [Vec<f64>; 3])| [t3.inner(u, u) + t3.inner(b, b), t3.inner(u, b), helicity(b)];
        let start = invariants3(&(u0.clone(), b0.clone()));
        let drift3 = |steps| {
            let end = invariants3(&t3.mhd(&u0, &b0, 0.0, 0.0, (0.0, 1.0), steps, zero_force(&t3)));
            (0..3).map(|k| ((end[k] - start[k]) / start[k]).abs()).fold(0.0, f64::max)
        };
        let (coarse, fine) = (drift3(20), drift3(40));
        assert!(coarse < 1e-5 && coarse / fine > 14.0 && start.iter().all(|x| x.abs() > 0.05), "{coarse:e} {fine:e} {start:?}");
        let t2 = plane(16);
        let (u, b) = (perp(&t2, &t2.sample(|[x, y]| (x + 2.0 * y).sin() + (3.0 * x).cos())), perp(&t2, &t2.sample(|[x, y]| 0.8 * (2.0 * x - y).cos() + (3.0 * y).sin() + 0.5 * (x + 2.0 * y).sin())));
        let flux = |b: &[Vec<f64>; 2]| {
            let [gx, gy] = [t2.gradient(&b[0]), t2.gradient(&b[1])];
            let j: Vec<f64> = gy[0].iter().zip(&gx[1]).map(|(dbydx, dbxdy)| dbydx - dbxdy).collect();
            t2.inverse_laplacian(&j.iter().map(|x| -x).collect::<Vec<_>>())
        };
        let a2 = |b: &[Vec<f64>; 2]| t2.mean(&flux(b).iter().map(|x| x * x).collect::<Vec<_>>());
        // The flux function is right: ∇⊥A gives back b.
        assert!((0..2).all(|c| max_diff(&perp(&t2, &flux(&b))[c], &b[c]) < 1e-14));
        let drift2 = |steps| {
            let (ue, be) = t2.mhd(&u, &b, 0.0, 0.0, (0.0, 1.0), steps, zero_force(&t2));
            [((t2.inner(&ue, &ue) + t2.inner(&be, &be)) / (t2.inner(&u, &u) + t2.inner(&b, &b)) - 1.0).abs(), (t2.inner(&ue, &be) / t2.inner(&u, &b) - 1.0).abs(), (a2(&be) / a2(&b) - 1.0).abs()]
                .into_iter()
                .fold(0.0, f64::max)
        };
        let (coarse, fine) = (drift2(80), drift2(160));
        assert!(coarse < 1e-5 && coarse / fine > 14.0, "{coarse:e} {fine:e}");
    }

    #[test]
    fn induction_matches_mhd_for_a_weak_field() {
        // A field too weak to act back: in a steadily forced ABC flow (f = νu holds it in place,
        // exactly, as above) the kinematic solver and the full one agree to the Lorentz term's
        // order, b², here 10⁻⁸ against b = 10⁻⁴. ∇·b stays at round-off.
        let torus = Torus { n: [8, 8, 8], length: TAU };
        let (u, nu, eta) = (abc(&torus), 0.2, 0.25);
        let (_, b0) = turbulent_pair(&torus);
        let b0 = b0.map(|c| c.iter().map(|x| 1e-4 * x).collect::<Vec<_>>());
        let force = u.clone().map(|c| c.iter().map(|x| nu * x).collect::<Vec<_>>());
        let (_, b_full) = torus.mhd(&u, &b0, nu, eta, (0.0, 2.0), 20, |_| force.clone());
        let b_kin = torus.induction(&b0, eta, (0.0, 2.0), 20, |_| u.clone());
        assert!((0..3).all(|c| max_diff(&b_full[c], &b_kin[c]) < 1e-10 && max_abs(&b_kin[c]) > 1e-5));
        assert!(max_abs(&torus.divergence(&b_kin)) < 1e-17);
    }
}
