//! The Newtonian potential of uniform bodies, in closed form: what a survey's forward model sums
//! over a grid of prisms.
//!
//! Every function returns, at unit `Gρ`, the potential `V = ∫dV/r`, its gradient `∇V` (the
//! attraction, pointing at the body) and its Hessian (the gravity-gradient tensor). Multiply by
//! `Gρ` for SI units. The Hessian's trace is Poisson's `−4π` inside the body and 0 outside, and
//! on a face it is the limit from outside.
//!
//! ```
//! use core::f64::consts::PI;
//! use geopotential::bodies::prism;
//! use metrology::G;
//!
//! // The Bouguer gradient: a 1 m layer of 2670 kg/m³ rock, 4000 km wide, pulls 0.112 mGal on top.
//! let (_, g, _) = prism([0.0; 3], [-2e6, -2e6, -1.0], [2e6, 2e6, 0.0]);
//! assert!((-g[2] * G * 2670.0 / 1e-5 - 0.11197).abs() < 1e-5);
//! assert!((-g[2] / (2.0 * PI) - 1.0).abs() < 1e-6); // 2πGρ per metre
//! ```

/// Potential, gradient and Hessian, at unit `Gρ`.
pub type Newtonian = (f64, [f64; 3], [[f64; 3]; 3]);

/// Beyond this many body diagonals from its centre the corner sums lose more to cancellation
/// (2×10⁻¹⁰ for a 1×1×10 prism) than a 4×4×4 Gauss–Legendre cubature does (4×10⁻¹²).
// shortcut: the corner sums lose about 10⁻¹⁶ r³/(abc), so a 1×1×100 needle loses 2×10⁻⁸ near the
// hand-over; split long prisms into near-cubes if that matters.
const FAR: f64 = 10.0;

/// A rectangular prism with faces normal to the axes, from corner `min` to corner `max`, seen
/// from point `p`, by the corner sums of Nagy, Papp & Benedek (J. Geodesy 74, 552, 2000):
/// `V = |||XY ln(Z+R) + YZ ln(X+R) + ZX ln(Y+R) − (X²/2)atan(YZ/XR) − (Y²/2)atan(ZX/YR) −
/// (Z²/2)atan(XY/ZR)|||`, with `X = x − p_x` at each corner and `+` at `max`.
///
/// It is finite on faces, edges and their extensions, where only the off-diagonal gradients
/// along an edge are infinite, as they are physically. Far away it hands over to the cubature.
pub fn prism(p: [f64; 3], min: [f64; 3], max: [f64; 3]) -> Newtonian {
    let diagonal = (0..3).map(|a| (max[a] - min[a]).powi(2)).sum::<f64>().sqrt();
    let distance = (0..3).map(|a| (p[a] - 0.5 * (min[a] + max[a])).powi(2)).sum::<f64>().sqrt();
    if distance > FAR * diagonal {
        return cubature(p, min, max);
    }
    // A zero coefficient times the infinite log on an edge counts as zero.
    let times = |k: f64, v: f64| if k == 0.0 { 0.0 } else { k * v };
    let (mut v, mut g, mut t) = (0.0, [0.0; 3], [[0.0; 3]; 3]);
    for corner in 0..8u32 {
        let upper = [corner & 1 == 1, corner & 2 == 2, corner & 4 == 4];
        let s = if corner.count_ones() % 2 == 1 { 1.0 } else { -1.0 };
        // −0 at an upper face makes a point on the face take the outside limit.
        let x: [f64; 3] = core::array::from_fn(|a| {
            let x = if upper[a] { max[a] } else { min[a] } - p[a];
            if upper[a] && x == 0.0 { -0.0 } else { x }
        });
        let r = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
        let (mut ln, mut at) = ([0.0; 3], [0.0; 3]);
        for a in 0..3 {
            let (b, c) = ((a + 1) % 3, (a + 2) % 3);
            // ln(x + r) without cancellation. Beyond the prism both corners of a pair have x ≤ 0,
            // and −ln(r − x) has the same pair difference; between them x + r = ρ²/(r − x).
            ln[a] = if p[a] >= max[a] {
                -(r - x[a]).ln()
            } else if x[a] >= 0.0 {
                (x[a] + r).ln()
            } else {
                ((x[b] * x[b] + x[c] * x[c]) / (r - x[a])).ln()
            };
            let num = x[b] * x[c];
            at[a] = if num == 0.0 { 0.0 } else { (num / (x[a] * r)).atan() };
        }
        for a in 0..3 {
            let (b, c) = ((a + 1) % 3, (a + 2) % 3);
            v += s * (times(x[b] * x[c], ln[a]) - times(0.5 * x[a] * x[a], at[a]));
            g[a] -= s * (times(x[b], ln[c]) + times(x[c], ln[b]) - times(x[a], at[a]));
            t[a][a] -= s * at[a];
            t[b][c] += s * ln[a];
            t[c][b] = t[b][c];
        }
    }
    (v, g, t)
}

/// The same integrals by a 4×4×4 Gauss–Legendre cubature of `1/r`, which converges like
/// `(half-width/2r)⁸` and so is at round-off far from the body, where the corner sums cancel.
fn cubature(p: [f64; 3], min: [f64; 3], max: [f64; 3]) -> Newtonian {
    const NODES: [(f64, f64); 4] = [
        (-0.861_136_311_594_052_6, 0.347_854_845_137_453_9),
        (-0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
        (0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
        (0.861_136_311_594_052_6, 0.347_854_845_137_453_9),
    ];
    let half: [f64; 3] = core::array::from_fn(|a| 0.5 * (max[a] - min[a]));
    let (mut v, mut g, mut t) = (0.0, [0.0; 3], [[0.0; 3]; 3]);
    for &(u, wu) in &NODES {
        for &(w, ww) in &NODES {
            for &(z, wz) in &NODES {
                let weight = wu * ww * wz * half[0] * half[1] * half[2];
                let d: [f64; 3] = core::array::from_fn(|a| 0.5 * (min[a] + max[a]) + half[a] * [u, w, z][a] - p[a]);
                let r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
                let r = r2.sqrt();
                v += weight / r;
                for a in 0..3 {
                    g[a] += weight * d[a] / (r2 * r);
                    for b in 0..3 {
                        let delta = if a == b { r2 } else { 0.0 };
                        t[a][b] += weight * (3.0 * d[a] * d[b] - delta) / (r2 * r2 * r);
                    }
                }
            }
        }
    }
    (v, g, t)
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::f64::consts::PI;

    const MIN: [f64; 3] = [-0.5; 3];
    const MAX: [f64; 3] = [0.5; 3];

    fn trace(t: &[[f64; 3]; 3]) -> f64 {
        t[0][0] + t[1][1] + t[2][2]
    }

    #[test]
    fn a_unit_cube() {
        // At the centre, ∫dV/r = 3 ln(2 + √3) − π/2 in closed form, and by symmetry the Hessian
        // is a third of Poisson's −4π on each axis, with no attraction.
        let (v, g, t) = prism([0.0; 3], MIN, MAX);
        assert!((v - (3.0 * (2.0 + 3f64.sqrt()).ln() - PI / 2.0)).abs() < 1e-15, "{v}");
        assert!(g.iter().all(|x| x.abs() < 1e-15));
        for (a, row) in t.iter().enumerate() {
            for (b, &x) in row.iter().enumerate() {
                assert!((x - if a == b { -4.0 * PI / 3.0 } else { 0.0 }).abs() < 1e-14, "T[{a}][{b}] = {x}");
            }
        }
        // V scales as size²: at a corner of the unit cube it is half the centre value (a corner
        // of each of the centre's eight half-size cubes, each a quarter of the unit cube's).
        let (corner, g, _) = prism(MAX, MIN, MAX);
        assert!((corner - v / 2.0).abs() < 1e-15, "{corner}");
        assert!(g.iter().all(|x| x.is_finite()) && (g[0] - g[1]).abs() < 1e-15);
    }

    #[test]
    fn poisson_inside_laplace_outside() {
        let (min, max) = ([0.3, -2.0, 5.0], [1.3, -1.0, 15.0]);
        for p in [[0.8, -1.5, 10.0], [0.31, -1.99, 14.9], [1.0, -1.2, 5.5]] {
            assert!((trace(&prism(p, min, max).2) + 4.0 * PI).abs() < 1e-13, "inside {p:?}");
        }
        for p in [[2.0, 0.0, 0.0], [0.8, -1.5, 15.5], [-3.0, 4.0, 9.0], [0.8, -1.5, 15.0]] {
            assert!(trace(&prism(p, min, max).2).abs() < 1e-13, "outside or on the surface {p:?}");
        }
    }

    #[test]
    fn derivatives_match_finite_differences() {
        let (min, max) = ([0.3, -2.0, 5.0], [1.3, -1.0, 15.0]);
        for p in [[2.0, 0.5, 3.0], [0.8, -1.5, 10.0], [0.7, -3.0, 16.0]] {
            let (_, g, t) = prism(p, min, max);
            let h = 1e-4;
            for a in 0..3 {
                let (mut up, mut down) = (p, p);
                up[a] += h;
                down[a] -= h;
                let (vu, gu, _) = prism(up, min, max);
                let (vd, gd, _) = prism(down, min, max);
                assert!(((vu - vd) / (2.0 * h) - g[a]).abs() < 1e-8, "g[{a}] at {p:?}");
                for b in 0..3 {
                    assert!(((gu[b] - gd[b]) / (2.0 * h) - t[a][b]).abs() < 1e-7, "T[{a}][{b}] at {p:?}");
                }
            }
        }
    }

    #[test]
    fn edge_lines_and_faces_are_continuous() {
        // On an edge's extension, or a face's plane outside it, the field is analytic: the value
        // there matches a point 10⁻⁹ off it, with nothing infinite.
        let (min, max) = ([0.3, -2.0, 5.0], [1.3, -1.0, 15.0]);
        for p in [[3.0, -2.0, 5.0], [-1.0, -1.0, 15.0], [0.3, -2.0, 20.0], [1.3, 1.0, 10.0], [0.8, -1.5, 15.0]] {
            let (v, g, t) = prism(p, min, max);
            let (vn, gn, tn) = prism([p[0] + 1e-9, p[1] + 2e-9, p[2] + 3e-9], min, max);
            assert!((v - vn).abs() < 1e-7 && (0..3).all(|a| (g[a] - gn[a]).abs() < 1e-6), "{p:?}");
            let on_face = p[2] == 15.0 && p[0] > 0.3 && p[0] < 1.3;
            for a in 0..3 {
                for b in 0..3 {
                    // Crossing the top face from outside, T_zz jumps by 4π.
                    let tol = if on_face && a == 2 && b == 2 { 4.0 * PI + 1e-6 } else { 1e-6 };
                    assert!(t[a][b].is_finite() && (t[a][b] - tn[a][b]).abs() < tol, "T[{a}][{b}] at {p:?}");
                }
            }
        }
    }

    #[test]
    fn far_field() {
        // At the hand-over the corner sums and the cubature agree, so the field is continuous.
        let (min, max) = ([0.3, -2.0, 5.0], [1.3, -1.0, 15.0]);
        let diagonal = 102f64.sqrt();
        let centre = [0.8, -1.5, 10.0];
        let dir = [0.36, 0.48, 0.8];
        for k in [3.0, 9.99] {
            let p: [f64; 3] = core::array::from_fn(|a| centre[a] + k * diagonal * dir[a]);
            let ((v, g, t), (vc, gc, tc)) = (prism(p, min, max), cubature(p, min, max));
            let rel = if k < 5.0 { 1e-6 } else { 1e-9 };
            assert!((v / vc - 1.0).abs() < rel, "V at {k}");
            assert!((0..3).all(|a| (g[a] / gc[a] - 1.0).abs() < rel), "g at {k}");
            assert!((t[2][2] / tc[2][2] - 1.0).abs() < 10.0 * rel, "T at {k}");
        }
        // Far off, a cube is a point mass to its hexadecapole, (a/r)⁴.
        let p = [600.0, 0.0, 800.0];
        let (v, g, t) = prism(p, MIN, MAX);
        assert!((v * 1000.0 - 1.0).abs() < 1e-12 && (g[2] / (-800.0 / 1e9) - 1.0).abs() < 1e-12);
        assert!((t[2][2] / ((3.0 * 800.0 * 800.0 - 1e6) / 1e15) - 1.0).abs() < 1e-12);
    }
}
