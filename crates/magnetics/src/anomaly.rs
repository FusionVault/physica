//! Magnetic anomalies of uniformly magnetised bodies, from their gravity by Poisson's relation.
//!
//! A uniform magnetisation `M` (A/m) in a body whose unit-density Newtonian potential has the
//! Hessian `T` at a point makes the field `B = (μ₀/4π)(T − I tr T)·M` there (Poisson 1826;
//! Blakely, *Potential Theory in Gravity and Magnetic Applications*, 1995, ch. 5). Outside,
//! where `tr T = 0`, this is `μ₀H`. Inside, `tr T = −4π`, and the trace adds the magnetisation's
//! own `μ₀M`. So every body in [`geopotential::bodies`] has a magnetic anomaly too, and `−T/4π`
//! inside it is the demagnetising tensor: `1/3` on each axis at the centre of a cube or a sphere.
//!
//! Vectors are in any one right-handed frame; geophysics uses x north, y east, z down.
//!
//! ```
//! use geopotential::bodies::prism;
//! use magnetics::anomaly::{body, sphere, total_field};
//!
//! // A 5 m sphere magnetised at 1 A/m along a vertical field, its centre 20 m down: 13.09 nT
//! // straight above it (Blakely's 2μ₀m/4πz³).
//! let earth = [0.0, 0.0, 50_000e-9];
//! let b = sphere([0.0; 3], [0.0, 0.0, 20.0], 5.0, [0.0, 0.0, 1.0]);
//! assert!((total_field(earth, b) / 1e-9 - 13.09).abs() < 0.01);
//!
//! // Any prism works the same way through its gravity-gradient tensor.
//! let t = prism([0.0; 3], [-5.0, -5.0, 15.0], [5.0, 5.0, 25.0]).2;
//! assert!(total_field(earth, body(t, [0.0, 0.0, 1.0])) > 0.0);
//! ```

use core::f64::consts::PI;

use geopotential::linalg::gauss;
use geopotential::special::carlson_rd;
use metrology::MU_0;

/// `μ₀/4π` (T m/A).
const C_M: f64 = MU_0 / (4.0 * PI);

/// The field (T) at a point of a body magnetised uniformly at `m` (A/m), from the Hessian of the
/// body's unit-density potential there, inside or outside: `B = (μ₀/4π)(T − I tr T)·m`.
pub fn body(hessian: [[f64; 3]; 3], m: [f64; 3]) -> [f64; 3] {
    let trace = hessian[0][0] + hessian[1][1] + hessian[2][2];
    core::array::from_fn(|i| C_M * ((0..3).map(|j| hessian[i][j] * m[j]).sum::<f64>() - trace * m[i]))
}

/// A sphere of radius `a` (m) centred at `centre`, magnetised uniformly at `m` (A/m), seen from
/// `p`: outside, the dipole of moment `(4/3)πa³m`; inside, the uniform `(2/3)μ₀m`.
pub fn sphere(p: [f64; 3], centre: [f64; 3], a: f64, m: [f64; 3]) -> [f64; 3] {
    let d: [f64; 3] = core::array::from_fn(|i| p[i] - centre[i]);
    let r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
    let t = if r2 < a * a {
        core::array::from_fn(|i| core::array::from_fn(|j| if i == j { -4.0 * PI / 3.0 } else { 0.0 }))
    } else {
        let k = 4.0 * PI * a.powi(3) / (3.0 * r2 * r2 * r2.sqrt());
        core::array::from_fn(|i| core::array::from_fn(|j| k * (3.0 * d[i] * d[j] - if i == j { r2 } else { 0.0 })))
    };
    body(t, m)
}

/// An infinite cylinder of radius `a` (m) along the direction `axis` through `point`, magnetised
/// uniformly at `m` (A/m), seen from `p`. Outside, it is the two-dimensional dipole of the line
/// moment `πa²m`; inside, `μ₀m/2` across the axis and `μ₀m` along it.
pub fn cylinder(p: [f64; 3], point: [f64; 3], axis: [f64; 3], a: f64, m: [f64; 3]) -> [f64; 3] {
    let n = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
    let u = axis.map(|x| x / n);
    let d: [f64; 3] = core::array::from_fn(|i| p[i] - point[i]);
    let along = d[0] * u[0] + d[1] * u[1] + d[2] * u[2];
    let rho: [f64; 3] = core::array::from_fn(|i| d[i] - along * u[i]);
    let r2 = rho[0] * rho[0] + rho[1] * rho[1] + rho[2] * rho[2];
    // The projector across the axis, I − uuᵀ.
    let across = |i: usize, j: usize| if i == j { 1.0 } else { 0.0 } - u[i] * u[j];
    let t = if r2 < a * a {
        core::array::from_fn(|i| core::array::from_fn(|j| -2.0 * PI * across(i, j)))
    } else {
        let k = 2.0 * PI * a * a / (r2 * r2);
        core::array::from_fn(|i| core::array::from_fn(|j| k * (2.0 * rho[i] * rho[j] - r2 * across(i, j))))
    };
    body(t, m)
}

/// The demagnetising factors `(N_a, N_b, N_c)` of a uniformly magnetised ellipsoid with semi-axes
/// `a`, `b`, `c` along its own axes, in any order: `N_c = (abc/3)R_D(a², b², c²)` (Carlson's `R_D`
/// is the Landau–Lifshitz integral `(abc/2)∫ds/((s + c²)√((s + a²)(s + b²)(s + c²)))`), and the
/// others by permutation. They sum to 1, and are `1/3` each for a sphere, `0.1736` along the long
/// axis of a 2:1 prolate spheroid and `0.9845` across a 100:1 disc (Chen, Brug & Goldfarb 1991,
/// Table I; Osborn 1945). The field inside is `H = H₀ − N·M`.
pub fn ellipsoid_factors(a: f64, b: f64, c: f64) -> [f64; 3] {
    let factor = |x: f64, y: f64, z: f64| x * y * z / 3.0 * carlson_rd(x * x, y * y, z * z);
    [factor(b, c, a), factor(c, a, b), factor(a, b, c)]
}

/// The magnetisation (A/m) of a body of susceptibility `chi` (SI, isotropic) and remanence
/// `remanent` (A/m) in the field `h0 = F/μ₀` (A/m), given its demagnetising tensor `n` (diagonal
/// for an ellipsoid in its own axes, [`ellipsoid_factors`]): `M = (I + χN)⁻¹(χH₀ + M_r)`
/// (Clark, Saul & Emerson 1986, Eq. 38; Takahashi & Oliveira 2017, Eqs. 49–50). The remanence
/// is demagnetised too, so pass the intrinsic value, not a measurement on the body. For `χ ≪ 1`
/// it is the classical `χH₀ + M_r`; the error of using that for a sphere is `χ/3`. A prism's
/// field is not uniform, so this is exact for ellipsoids only.
pub fn magnetisation(chi: f64, h0: [f64; 3], remanent: [f64; 3], n: [[f64; 3]; 3]) -> [f64; 3] {
    let rows: Vec<Vec<f64>> = (0..3)
        .map(|i| (0..3).map(|j| if i == j { 1.0 } else { 0.0 } + chi * n[i][j]).chain([chi * h0[i] + remanent[i]]).collect())
        .collect();
    let m = gauss(rows);
    [m[0], m[1], m[2]]
}

/// The Königsberger ratio `Q = |M_r|/(χ|H₀|)` of remanent to induced magnetisation: above 1 the
/// remanence sets the anomaly's shape (Koenigsberger 1938).
pub fn koenigsberger(remanent: [f64; 3], chi: f64, h0: [f64; 3]) -> f64 {
    let norm = |v: [f64; 3]| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    norm(remanent) / (chi * norm(h0))
}

/// What a total-field magnetometer reads (T) when an anomaly `b` adds to the ambient field `f`:
/// `|F + B| − |F|`, computed as `(2F·B + B²)/(|F + B| + |F|)` without cancellation. It is `B·F̂`
/// to first order in `|B|/|F|`.
pub fn total_field(f: [f64; 3], b: [f64; 3]) -> f64 {
    let dot = |x: [f64; 3], y: [f64; 3]| x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
    let sum = [f[0] + b[0], f[1] + b[1], f[2] + b[2]];
    (2.0 * dot(f, b) + dot(b, b)) / (dot(sum, sum).sqrt() + dot(f, f).sqrt())
}

#[cfg(test)]
mod tests {
    use super::*;
    use geopotential::bodies::prism;

    const NT: f64 = 1e-9;

    #[test]
    fn a_sphere_at_the_pole_and_the_equator() {
        // m = (4/3)π 5³ × 1 A/m at 20 m: 2C_m m/z³ = 13.09 nT along the moment, −C_m m/z³ =
        // −6.545 nT across it.
        let m = 4.0 / 3.0 * PI * 125.0;
        let b = sphere([0.0; 3], [0.0, 0.0, 20.0], 5.0, [0.0, 0.0, 1.0]);
        assert!((b[2] / (2.0 * C_M * m / 8000.0) - 1.0).abs() < 1e-14 && (b[2] / NT - 13.09).abs() < 1e-3);
        let pole = total_field([0.0, 0.0, 50_000.0 * NT], b);
        assert!((pole - b[2]).abs() < 2e-3 * NT); // the second-order term B²/2F
        let b = sphere([0.0; 3], [0.0, 0.0, 20.0], 5.0, [1.0, 0.0, 0.0]);
        let equator = total_field([30_000.0 * NT, 0.0, 0.0], b);
        assert!((equator / NT + 6.545).abs() < 1e-3 && (equator / b[0] - 1.0).abs() < 1e-14);
        // Inside, (2/3)μ₀M.
        let inside = sphere([1.0, 2.0, 19.0], [0.0, 0.0, 20.0], 5.0, [0.0, 0.0, 1.0]);
        assert!((inside[2] - 2.0 / 3.0 * MU_0).abs() < 1e-20 && inside[0] == 0.0);
    }

    #[test]
    fn demagnetising_factors_of_ellipsoids() {
        // Chen, Brug & Goldfarb (1991), Table I: N along the axis of a spheroid of axis 2l and
        // diameter 2a, γ = l/a, to the printed digits; 1/3 for a sphere; the three sum to 1.
        for (gamma, want, digits) in [(0.01, 0.9845, 4), (0.1, 0.8608, 4), (0.5, 0.5272, 4), (1.0, 0.3333, 4), (2.0, 0.1736, 4), (5.0, 0.05582, 5), (10.0, 0.02029, 5), (100.0, 0.0004299, 7)] {
            let n = ellipsoid_factors(1.0, 1.0, gamma);
            assert!((n[2] - want).abs() < 0.5 * 10f64.powi(-digits) + 1e-12, "γ = {gamma}: {}", n[2]);
            assert!((n[0] + n[1] + n[2] - 1.0).abs() < 1e-14 && (n[0] - n[1]).abs() < 1e-15);
        }
        assert!((ellipsoid_factors(2.0, 2.0, 2.0)[1] - 1.0 / 3.0).abs() < 1e-15);
        // Takahashi & Oliveira (2017), Table 2: Warrego, a = 490.7, b = 69.7, c = 30 m, χ = 1.69.
        // Their Eq. 60 puts the susceptibility above which an 8% error needs demagnetisation at
        // ε/N_max = 0.116 SI; the order of the axes does not matter.
        let n = ellipsoid_factors(490.7, 69.7, 30.0);
        assert!((n[2] - 0.689_521).abs() < 1e-6 && (0.08 / n[2] - 0.116).abs() < 5e-4, "{n:?}");
        assert!((ellipsoid_factors(30.0, 490.7, 69.7)[0] - n[2]).abs() < 1e-15);
    }

    #[test]
    fn induced_and_remanent_magnetisation() {
        // A sphere: M = χH₀/(1 + χ/3), 0.75 of the classical value at χ = 1, 0.2308 at χ = 10; the
        // remanence is demagnetised by the same factor, and Q is the ratio of the two inputs.
        let n3 = [[1.0 / 3.0, 0.0, 0.0], [0.0, 1.0 / 3.0, 0.0], [0.0, 0.0, 1.0 / 3.0]];
        let h0 = [0.0, 0.0, 50e-6 / MU_0]; // 50 000 nT: 39.79 A/m
        assert!((h0[2] - 39.79).abs() < 0.01);
        for (chi, factor) in [(0.1, 0.967_742), (1.0, 0.75), (10.0, 0.230_769)] {
            let m = magnetisation(chi, h0, [0.0; 3], n3);
            assert!((m[2] / (chi * h0[2]) - factor).abs() < 1e-6 && m[0] == 0.0, "χ = {chi}: {m:?}");
        }
        let m = magnetisation(1.0, h0, [3.0, 0.0, 0.0], n3);
        assert!((m[0] - 2.25).abs() < 1e-12 && (m[2] - 0.75 * h0[2]).abs() < 1e-12);
        assert!((koenigsberger([3.0, 0.0, 0.0], 1.0, h0) - 3.0 / h0[2]).abs() < 1e-15);
        // Takahashi & Oliveira's confocal pair (Table 1): the (900, 500, 100) m ellipsoid at χ = 1.2
        // and the confocal (1676.31, 1500, 1417.74) m one at χ′ = 0.01415 have the same dipole
        // moment along a, so the same anomaly: a volume 79 times larger at a susceptibility 85
        // times smaller.
        let (a, b, c, chi, u) = (900.0f64, 500.0f64, 100.0f64, 1.2, 2e6);
        let (a2, b2, c2) = ((a * a + u).sqrt(), (b * b + u).sqrt(), (c * c + u).sqrt());
        assert!((a2 - 1676.31).abs() < 0.01 && (c2 - 1417.74).abs() < 0.01);
        let diag = |n: [f64; 3]| [[n[0], 0.0, 0.0], [0.0, n[1], 0.0], [0.0, 0.0, n[2]]];
        let h = [1.0, 0.0, 0.0];
        let moment = a * b * c * magnetisation(chi, h, [0.0; 3], diag(ellipsoid_factors(a, b, c)))[0];
        let chi2 = 0.014_15;
        let moment2 = a2 * b2 * c2 * magnetisation(chi2, h, [0.0; 3], diag(ellipsoid_factors(a2, b2, c2)))[0];
        assert!((moment2 / moment - 1.0).abs() < 1e-3 && (a2 * b2 * c2 / (a * b * c) - 79.2).abs() < 0.1 && (chi / chi2 - 84.8).abs() < 0.1);
        // At χ ≪ 1 the classical χH₀ + M_r returns; a cube's centre (−T/4π = I/3) gives the
        // sphere's factor, by the same route as the sphere.
        let weak = magnetisation(1e-4, h0, [0.1, 0.0, 0.0], n3);
        assert!((weak[2] / (1e-4 * h0[2]) - 1.0).abs() < 4e-5 && (weak[0] - 0.1).abs() < 1e-5);
        let t = prism([0.0; 3], [-0.5; 3], [0.5; 3]).2;
        let cube: [[f64; 3]; 3] = core::array::from_fn(|i| core::array::from_fn(|j| -t[i][j] / (4.0 * PI)));
        assert!((magnetisation(1.0, h0, [0.0; 3], cube)[2] / (0.75 * h0[2]) - 1.0).abs() < 1e-13);
    }

    #[test]
    fn a_cube_through_its_gravity() {
        // At a cube's centre H = −M/3, as in a sphere.
        let m = [0.3, -0.5, 0.8];
        let b = body(prism([0.0; 3], [-0.5; 3], [0.5; 3]).2, m);
        assert!((0..3).all(|i| (b[i] - 2.0 / 3.0 * MU_0 * m[i]).abs() < 1e-20));
        // Far off, a cube is a sphere of the same volume to (a/r)⁴, having no quadrupole:
        // 9×10⁻⁸ at 37 sides.
        let p = [12.0, -16.0, 30.0];
        let (bc, bs) = (body(prism(p, [-0.5; 3], [0.5; 3]).2, m), sphere(p, [0.0; 3], (3.0 / (4.0 * PI)).cbrt(), m));
        assert!((0..3).all(|i| (bc[i] / bs[i] - 1.0).abs() < 3e-7), "{bc:?} {bs:?}");
    }

    #[test]
    fn boundary_conditions_on_a_face() {
        // Across the top face of a prism the normal B and the tangential H = B/μ₀ − M are
        // continuous, which only holds if the trace term adds μ₀M exactly inside.
        let (min, max) = ([-1.0, -2.0, 0.0], [3.0, 1.0, 2.0]);
        let m = [0.4, 0.7, -0.6];
        let (p, h) = ([0.5, -0.3, 2.0], 1e-7);
        let out = body(prism([p[0], p[1], p[2] + h], min, max).2, m);
        let inn = body(prism([p[0], p[1], p[2] - h], min, max).2, m);
        assert!((out[2] - inn[2]).abs() < 1e-12, "B_z {} vs {}", out[2], inn[2]);
        for i in 0..2 {
            assert!((out[i] / MU_0 - (inn[i] / MU_0 - m[i])).abs() < 1e-6, "H_{i}");
        }
        // The jump is all in the tangential B: μ₀M_t.
        assert!((inn[0] - out[0] - MU_0 * m[0]).abs() < 1e-12);
    }

    #[test]
    fn a_cylinder() {
        // Far off, a long square prism of the same area: the square has no 2-D quadrupole, so
        // they differ by (a/r)⁴, 5×10⁻⁶ at 20 radii and 3×10⁻⁷ at 40.
        let (a, m) = (1.0, [0.6, 0.0, 0.8]);
        let s = (PI * a * a).sqrt() / 2.0;
        let p = [24.0, 0.0, -32.0];
        let bc = cylinder(p, [0.0; 3], [0.0, 1.0, 0.0], a, m);
        let bp = body(prism(p, [-s, -1e6, -s], [s, 1e6, s]).2, m);
        assert!((0..3).all(|i| (bc[i] - bp[i]).abs() < 1e-6 * bc[0].abs().max(bc[2].abs())), "{bc:?} {bp:?}");
        assert!(bc[1].abs() < 1e-30);
        // Inside, μ₀M/2 across the axis and μ₀M along it.
        let inside = cylinder([0.2, 5.0, 0.1], [0.0; 3], [0.0, 2.0, 0.0], a, [1.0, 1.0, 0.0]);
        assert!((inside[0] - MU_0 / 2.0).abs() < 1e-20 && (inside[1] - MU_0).abs() < 1e-20);
        // The 2-D dipole straight above: 2C_m(πa²M)/z² along a vertical moment.
        let b = cylinder([0.0; 3], [0.0, 0.0, 15.0], [0.0, 1.0, 0.0], 3.0, [0.0, 0.0, 1.0]);
        assert!((b[2] / (2.0 * C_M * PI * 9.0 / 225.0) - 1.0).abs() < 1e-14);
    }
}
