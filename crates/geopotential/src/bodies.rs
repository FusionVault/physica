//! The Newtonian potential of uniform bodies, in closed form: what a survey's forward model sums
//! over a grid of prisms, or over the triangles of an asteroid's shape model.
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
    // crate::gauss_legendre(4), written out for the inner loop (a test holds them equal).
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

/// A closed polyhedron of triangles `faces`, each listing three indices into `verts` counter-
/// clockwise seen from outside, by the edge and face sums of Werner & Scheeres (Celest. Mech.
/// Dyn. Astron. 65, 313, 1997):
/// `V = ½Σ_e r_e·E_e·r_e L_e − ½Σ_f (n̂_f·r_f)² ω_f`, `∇V = −Σ_e E_e r_e L_e + Σ_f n̂_f(n̂_f·r_f)ω_f`,
/// Hessian `Σ_e E_e L_e − Σ_f n̂_f n̂_fᵀ ω_f`. `E_e = n̂_A n̂_eᴬᵀ + n̂_B n̂_eᴮᵀ` from the two faces
/// at an edge, `L_e = ln((r_a + r_b + e)/(r_a + r_b − e))`, and `ω_f` is the solid angle a face
/// subtends at `p`.
///
/// The trace is `−Σω_f`, so it doubles as an inside test: `−4π` inside, 0 outside. Each face's
/// own dyad `n̂_f n̂_eᵀ` annihilates the edge, so edges are summed face by face, with no pairing.
///
/// `V` and `∇V` are exact on the surface too. The Hessian jumps by `4π` across a face, where it
/// takes the limit from outside when `p` lies exactly in the face's plane; on an edge, including
/// the inner edges of a triangulated face, it is not defined and some entries are not finite.
// shortcut: far away the sums cancel like the prism's, losing 2×10⁻⁹ at 100 sizes and 7×10⁻⁷ at
// 1000 for a cube mesh, with no hand-over; add a cubature over tetrahedra when a caller needs it.
pub fn polyhedron(p: [f64; 3], verts: &[[f64; 3]], faces: &[[usize; 3]]) -> Newtonian {
    let sub = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let cross = |a: [f64; 3], b: [f64; 3]| [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
    let unit = |a: [f64; 3]| {
        let n = dot(a, a).sqrt();
        [a[0] / n, a[1] / n, a[2] / n]
    };
    let (mut v, mut g, mut t) = (0.0, [0.0; 3], [[0.0; 3]; 3]);
    for face in faces {
        let corners = face.map(|i| verts[i]);
        let n = unit(cross(sub(corners[1], corners[0]), sub(corners[2], corners[0])));
        let r = corners.map(|c| sub(c, p));
        let l = r.map(|r| dot(r, r).sqrt());
        // Van Oosterom & Strackee's solid angle of a triangle; −0 in its plane makes a point on
        // the face see it from outside (−2π), the prism's convention.
        let triple = dot(r[0], cross(r[1], r[2]));
        let triple = if triple == 0.0 { -0.0 } else { triple };
        let omega = 2.0 * triple.atan2(l[0] * l[1] * l[2] + l[0] * dot(r[1], r[2]) + l[1] * dot(r[2], r[0]) + l[2] * dot(r[0], r[1]));
        let rf = dot(n, r[0]);
        v -= 0.5 * rf * rf * omega;
        for a in 0..3 {
            g[a] += n[a] * rf * omega;
            for b in 0..3 {
                t[a][b] -= n[a] * n[b] * omega;
            }
        }
        for k in 0..3 {
            let (i, j) = (k, (k + 1) % 3);
            let edge = sub(corners[j], corners[i]);
            let ne = unit(cross(edge, n));
            let (sum, e) = (l[i] + l[j], dot(edge, edge).sqrt());
            // Near the edge itself r_a + r_b − e cancels; 2|r_a × r_b|²/((r_a r_b − r_a·r_b)(r_a + r_b + e))
            // is the same number without the cancellation.
            let ab = dot(r[i], r[j]);
            let gap = if ab < 0.0 {
                let c = cross(r[i], r[j]);
                2.0 * dot(c, c) / ((l[i] * l[j] - ab) * (sum + e))
            } else {
                sum - e
            };
            let le = ((sum + e) / gap).ln();
            for a in 0..3 {
                for b in 0..3 {
                    if n[a] * ne[b] != 0.0 {
                        t[a][b] += n[a] * ne[b] * le;
                    }
                }
            }
            // On the edge itself n̂_e·r = 0 exactly and only the Hessian is infinite.
            if le.is_finite() {
                let ner = dot(ne, r[i]);
                v += 0.5 * dot(n, r[i]) * ner * le;
                for a in 0..3 {
                    g[a] -= n[a] * ner * le;
                }
            }
        }
    }
    for (a, b) in [(0, 1), (0, 2), (1, 2)] {
        let mean = 0.5 * (t[a][b] + t[b][a]);
        (t[a][b], t[b][a]) = (mean, mean);
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

    /// The unit cube as 12 outward triangles.
    fn cube_mesh() -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
        let verts: Vec<[f64; 3]> = (0..8).map(|i| [i & 1, i >> 1 & 1, i >> 2 & 1].map(|b| b as f64 - 0.5)).collect();
        // Each face as a quad counter-clockwise from outside: −x, +x, −y, +y, −z, +z.
        let quads = [[0, 4, 6, 2], [1, 3, 7, 5], [0, 1, 5, 4], [2, 6, 7, 3], [0, 2, 3, 1], [4, 5, 7, 6]];
        (verts, quads.iter().flat_map(|&[a, b, c, d]| [[a, b, c], [a, c, d]]).collect())
    }

    #[test]
    fn a_triangulated_cube_is_the_prism() {
        let (verts, faces) = cube_mesh();
        // Both sums cancel as (r/size)³ away from the body. At [3, 2, −4] a 40-point Gauss–Legendre
        // reference puts the prism 4.6×10⁻¹⁴ off in V and the polyhedron 1.5×10⁻¹⁴, so agreement
        // is relative to the field's size.
        for p in [[0.9, -0.4, 1.3], [0.1, 0.2, -0.05], [3.0, 2.0, -4.0], [0.5, 0.2, 0.7], [0.0; 3]] {
            let ((v, g, t), (vp, gp, tp)) = (polyhedron(p, &verts, &faces), prism(p, MIN, MAX));
            let (gs, ts) = (gp.iter().fold(0f64, |m, x| m.max(x.abs())), tp.iter().flatten().fold(0f64, |m, x| m.max(x.abs())));
            assert!((v - vp).abs() < 1e-13 * vp, "V at {p:?}: {v} vs {vp}");
            for a in 0..3 {
                assert!((g[a] - gp[a]).abs() < 2e-13 * gs + 1e-15, "g at {p:?}");
                for b in 0..3 {
                    assert!((t[a][b] - tp[a][b]).abs() < 1e-14 * ts.max(1.0), "T[{a}][{b}] at {p:?}");
                }
            }
        }
        // Where the prism has no cancellation the two agree to round-off, as the notes found.
        let ((v, g, _), (vp, gp, _)) = (polyhedron([0.9, -0.4, 1.3], &verts, &faces), prism([0.9, -0.4, 1.3], MIN, MAX));
        assert!((v - vp).abs() < 3e-16 && (g[2] - gp[2]).abs() < 1e-15);
        // The trace is −Σω: an inside test.
        assert!((trace(&polyhedron([0.1, 0.2, -0.05], &verts, &faces).2) + 4.0 * PI).abs() < 1e-14);
        assert!(trace(&polyhedron([0.9, -0.4, 1.3], &verts, &faces).2).abs() < 1e-14);
        // On a face both take the outside limit of the Hessian.
        for p in [[0.5, 0.2, 0.1], [0.2, 0.1, 0.5], [-0.3, -0.5, 0.4]] {
            let (t, tp) = (polyhedron(p, &verts, &faces).2, prism(p, MIN, MAX).2);
            assert!((0..9).all(|k| (t[k / 3][k % 3] - tp[k / 3][k % 3]).abs() < 1e-14), "{p:?}");
        }
        // On an edge, a vertex and a face's inner diagonal, V and g are exact.
        for p in [[0.5, 0.5, 0.1], [0.5, 0.5, 0.5], [0.5, 0.2, 0.2]] {
            let (v, g, _) = polyhedron(p, &verts, &faces);
            let (vp, gp, _) = prism(p, MIN, MAX);
            assert!((v - vp).abs() < 1e-14 && (0..3).all(|a| (g[a] - gp[a]).abs() < 1e-14), "{p:?}");
        }
    }

    #[test]
    fn a_turned_cube() {
        // Turn the cube about an oblique axis: V is unchanged, ∇V turns with it and the Hessian
        // turns as RTRᵀ.
        let (verts, faces) = cube_mesh();
        let (c, s) = (0.6f64, 0.8f64);
        let turn = |x: [f64; 3]| [c * x[0] - s * x[1], s * c * x[0] + c * c * x[1] - s * x[2], s * s * x[0] + s * c * x[1] + c * x[2]];
        let r: [[f64; 3]; 3] = core::array::from_fn(|i| turn([0, 1, 2].map(|j| if j == i { 1.0 } else { 0.0 })));
        let turned: Vec<[f64; 3]> = verts.iter().map(|&x| turn(x)).collect();
        let p = [0.7, -1.1, 0.4];
        let (vp, gp, tp) = prism(p, MIN, MAX);
        let (v, g, t) = polyhedron(turn(p), &turned, &faces);
        assert!((v - vp).abs() < 1e-14);
        let gt = turn(gp);
        assert!((0..3).all(|a| (g[a] - gt[a]).abs() < 1e-14));
        for a in 0..3 {
            for b in 0..3 {
                // turn(e_i) is column i of R, so R_ai = r[i][a].
                let want: f64 = (0..3).flat_map(|i| (0..3).map(move |j| (i, j))).map(|(i, j)| r[i][a] * tp[i][j] * r[j][b]).sum();
                assert!((t[a][b] - want).abs() < 1e-14, "T[{a}][{b}]");
            }
        }
    }

    #[test]
    fn an_octahedron_far_off() {
        // |x| + |y| + |z| ≤ 1 has volume 4/3 and cubic symmetry, so no quadrupole: from 30 sizes
        // it is a point mass to (1/r)⁴.
        let verts = [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0], [0.0, 0.0, -1.0]];
        let faces = [[0, 2, 4], [2, 1, 4], [1, 3, 4], [3, 0, 4], [2, 0, 5], [1, 2, 5], [3, 1, 5], [0, 3, 5]];
        let p = [18.0, 24.0, 40.0]; // |p| = 50
        let (v, g, t) = polyhedron(p, &verts, &faces);
        let m = 4.0 / 3.0;
        assert!((v / (m / 50.0) - 1.0).abs() < 1e-6);
        assert!((g[2] / (-m * 40.0 / 125e3) - 1.0).abs() < 1e-6);
        assert!((t[0][1] / (3.0 * m * 18.0 * 24.0 / 50f64.powi(5)) - 1.0).abs() < 1e-5);
        assert!((trace(&polyhedron([0.2, -0.3, 0.1], &verts, &faces).2) + 4.0 * PI).abs() < 1e-14);
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
    fn the_cubature_is_gauss_legendre() {
        // The far-field nodes are the shared 4-point rule; the far field test pins their accuracy.
        let (x, w) = crate::gauss_legendre(4);
        let ours = cubature([30.0, 0.0, 0.0], MIN, MAX);
        let mut v = 0.0;
        for (&a, &wa) in x.iter().zip(&w) {
            for (&b, &wb) in x.iter().zip(&w) {
                for (&c, &wc) in x.iter().zip(&w) {
                    v += 0.125 * wa * wb * wc / ((30.0 - 0.5 * a).powi(2) + 0.25 * b * b + 0.25 * c * c).sqrt();
                }
            }
        }
        assert!((ours.0 - v).abs() < 1e-17, "{} vs {v}", ours.0);
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
