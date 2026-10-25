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

/// The relative error the hand-over between the corner sums and the cubature aims at. The
/// corner sums lose `ε·R²(π/2 + |ln R|)` in absolute terms, `R` being the farthest corner's
/// distance, so near a long side of a needle they lose far more than this.
const TOL: f64 = 1e-13;

/// A rectangular prism with faces normal to the axes, from corner `min` to corner `max`, seen
/// from point `p`: the corner sums of Nagy, Papp & Benedek (J. Geodesy 74, 552, 2000),
/// `V = |||XY ln(Z+R) + YZ ln(X+R) + ZX ln(Y+R) − (X²/2)atan(YZ/XR) − (Y²/2)atan(ZX/YR) −
/// (Z²/2)atan(XY/ZR)|||`, with `X = x − p_x` at each corner and `+` at `max`, inside and on the
/// body, and a Gauss–Legendre cubature of `1/r` outside it wherever that is affordable.
///
/// The corner sums are finite on faces, edges and their extensions, where only the off-diagonal
/// gradients along an edge are infinite, as they are physically. Away from the body they cancel,
/// losing about `ε r³/(abc)`, cubically with distance, as Fukushima (GJI 220, 610, 2019) found.
/// An `n`-point Gauss–Legendre rule on each axis integrates the moments to degree `2n − 1`
/// exactly, so the cubature is the multipole expansion to that degree with a convergent
/// remainder, `ρ⁻²ⁿ` for the Bernstein ellipse `ρ = (d + √(d² + h²))/h` of the half-width `h` at
/// distance `d` from the box; each axis gets the order that puts it at round-off, up to 16 and
/// 512 points in all. Nearer a long side, where no rule that size suffices, an elongated prism is
/// bisected along its longest side until its pieces are either affordable or short enough for
/// the corner sums to hold 10⁻¹³. Splitting only happens with `p` outside the box, so no split
/// plane passes through `p`. On a 1×1×100 needle this holds 10⁻¹⁴ from 0.01 off its side to 100
/// Brillouin radii, where the corner sums alone lose 5×10⁻⁸.
pub fn prism(p: [f64; 3], min: [f64; 3], max: [f64; 3]) -> Newtonian {
    let half: [f64; 3] = core::array::from_fn(|a| 0.5 * (max[a] - min[a]));
    let distance = (0..3).map(|a| (min[a] - p[a]).max(p[a] - max[a]).max(0.0).powi(2)).sum::<f64>().sqrt();
    if distance > 0.0 {
        let orders: [usize; 3] = core::array::from_fn(|a| {
            let rho = (distance + distance.hypot(half[a])) / half[a];
            (-(f64::EPSILON.ln()) / (2.0 * rho.ln())).ceil().max(1.0) as usize
        });
        if orders.iter().all(|&n| n <= 16) && orders.iter().product::<usize>() <= 512 {
            return cubature(p, min, max, orders);
        }
    }
    let sums = corner_sums(p, min, max);
    let longest = (0..3).max_by(|&a, &b| half[a].total_cmp(&half[b])).unwrap();
    let shortest = half.iter().fold(f64::INFINITY, |m, &h| m.min(h));
    if distance > 0.0 && half[longest] > 2.0 * shortest {
        let farthest = (0..8u32).map(|c| (0..3).map(|a| (if c >> a & 1 == 1 { max[a] } else { min[a] } - p[a]).powi(2)).sum::<f64>().sqrt()).fold(0.0, f64::max);
        if f64::EPSILON * farthest * farthest * (core::f64::consts::FRAC_PI_2 + farthest.ln().abs()) > TOL * sums.0.abs() {
            let mid = 0.5 * (min[longest] + max[longest]);
            let (mut lower_max, mut upper_min) = (max, min);
            (lower_max[longest], upper_min[longest]) = (mid, mid);
            let ((v1, g1, t1), (v2, g2, t2)) = (prism(p, min, lower_max), prism(p, upper_min, max));
            return (v1 + v2, core::array::from_fn(|a| g1[a] + g2[a]), core::array::from_fn(|a| core::array::from_fn(|b| t1[a][b] + t2[a][b])));
        }
    }
    sums
}

/// The corner sums of Nagy, Papp & Benedek, as [`prism`] describes them.
fn corner_sums(p: [f64; 3], min: [f64; 3], max: [f64; 3]) -> Newtonian {
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

/// The same integrals by a Gauss–Legendre cubature of `1/r` with `orders` points on the three
/// axes.
fn cubature(p: [f64; 3], min: [f64; 3], max: [f64; 3], orders: [usize; 3]) -> Newtonian {
    let rules = orders.map(crate::gauss_legendre);
    let half: [f64; 3] = core::array::from_fn(|a| 0.5 * (max[a] - min[a]));
    let (mut v, mut g, mut t) = (0.0, [0.0; 3], [[0.0; 3]; 3]);
    for (&u, &wu) in rules[0].0.iter().zip(&rules[0].1) {
        for (&w, &ww) in rules[1].0.iter().zip(&rules[1].1) {
            for (&z, &wz) in rules[2].0.iter().zip(&rules[2].1) {
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
/// own dyad `n̂_f n̂_eᵀ` annihilates the edge, so edges are summed face by face, with no pairing,
/// except where `p` lies on an edge, where the two faces' dyads are summed first so that a flat
/// face's inner diagonal cancels. Check a mesh once with [`check_mesh`].
///
/// `V` and `∇V` are exact on the surface too. The Hessian jumps by `4π` across a face, where it
/// takes the limit from outside when `p` lies exactly in the face's plane, including on the inner
/// edges of a triangulated face; on a real edge some entries are infinite, as they are
/// physically.
///
/// Far away the sums cancel like the prism's corner sums, losing 2×10⁻¹⁰ at 100 sizes and
/// 3×10⁻⁷ at 1000 for a cube mesh. Beyond [`POLYHEDRON_FAR`] circumradii of the vertex centroid
/// they hand over to a Gauss–Legendre cubature over the signed tetrahedra from that centroid to
/// each face, collapsed onto a cube by Duffy's map, 6³ points each: exact for the moments to
/// degree 9, so a multipole expansion to that degree with a convergent remainder, 10⁻¹³ or better
/// at the hand-over on a cube and on a 14 744-face asteroid.
pub fn polyhedron(p: [f64; 3], verts: &[[f64; 3]], faces: &[[usize; 3]]) -> Newtonian {
    debug_assert!(check_mesh(verts, faces).is_ok(), "{:?}", check_mesh(verts, faces));
    let mut used = vec![false; verts.len()];
    faces.iter().flatten().for_each(|&i| used[i] = true);
    let count = used.iter().filter(|&&u| u).count() as f64;
    let centroid: [f64; 3] = core::array::from_fn(|a| verts.iter().zip(&used).filter(|&(_, &u)| u).map(|(v, _)| v[a]).sum::<f64>() / count);
    let radius = verts.iter().zip(&used).filter(|&(_, &u)| u).map(|(v, _)| (0..3).map(|a| (v[a] - centroid[a]).powi(2)).sum::<f64>()).fold(0.0, f64::max).sqrt();
    if (0..3).map(|a| (p[a] - centroid[a]).powi(2)).sum::<f64>().sqrt() > POLYHEDRON_FAR * radius {
        tetrahedra(p, centroid, verts, faces)
    } else {
        edge_face_sums(p, verts, faces)
    }
}

/// Beyond this many circumradii from the vertex centroid [`polyhedron`] uses the tetrahedral
/// cubature: there a cube mesh's sums have lost 7×10⁻¹⁵ and the cubature holds 6×10⁻¹⁴ or better.
pub const POLYHEDRON_FAR: f64 = 5.0;

/// The cubature of [`polyhedron`]'s far field: 6-point Gauss–Legendre on each axis of the cube
/// that Duffy's map `x = A + u[(B − A) + v((C − B) + w(D − C))]`, with Jacobian `6·Vol·u²v`,
/// collapses the tetrahedron `ABCD` onto, for each face `BCD` from the `apex` `A`. The volumes are
/// signed, so the overlaps of a non-convex body cancel, and the faces' sums are added with
/// Neumaier's compensation, since they are large and of both signs.
fn tetrahedra(p: [f64; 3], apex: [f64; 3], verts: &[[f64; 3]], faces: &[[usize; 3]]) -> Newtonian {
    let (x, w) = crate::gauss_legendre(6);
    let rule: Vec<(f64, f64)> = x.iter().zip(&w).map(|(&x, &w)| (0.5 * (x + 1.0), 0.5 * w)).collect();
    let mut sums = [(0.0f64, 0.0f64); 13]; // (sum, compensation) for V, ∇V and the Hessian
    for face in faces {
        let [b, c, d] = face.map(|i| verts[i]);
        let (ab, bc, cd): ([f64; 3], [f64; 3], [f64; 3]) = (core::array::from_fn(|a| b[a] - apex[a]), core::array::from_fn(|a| c[a] - b[a]), core::array::from_fn(|a| d[a] - c[a]));
        let ac: [f64; 3] = core::array::from_fn(|a| c[a] - apex[a]);
        let ad: [f64; 3] = core::array::from_fn(|a| d[a] - apex[a]);
        let vol6 = ab[0] * (ac[1] * ad[2] - ac[2] * ad[1]) + ab[1] * (ac[2] * ad[0] - ac[0] * ad[2]) + ab[2] * (ac[0] * ad[1] - ac[1] * ad[0]);
        let mut terms = [0.0f64; 13];
        for &(u, wu) in &rule {
            for &(v, wv) in &rule {
                for &(t, wt) in &rule {
                    let weight = vol6 * u * u * v * wu * wv * wt;
                    let r: [f64; 3] = core::array::from_fn(|a| apex[a] + u * (ab[a] + v * (bc[a] + t * cd[a])) - p[a]);
                    let r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
                    let rn = r2.sqrt();
                    terms[0] += weight / rn;
                    for a in 0..3 {
                        terms[1 + a] += weight * r[a] / (r2 * rn);
                        for bb in 0..3 {
                            let delta = if a == bb { r2 } else { 0.0 };
                            terms[4 + 3 * a + bb] += weight * (3.0 * r[a] * r[bb] - delta) / (r2 * r2 * rn);
                        }
                    }
                }
            }
        }
        for (k, &x) in terms.iter().enumerate() {
            let (sum, carry) = &mut sums[k];
            let next = *sum + x;
            *carry += if sum.abs() >= x.abs() { (*sum - next) + x } else { (x - next) + *sum };
            *sum = next;
        }
    }
    let total = sums.map(|(s, c)| s + c);
    (total[0], [total[1], total[2], total[3]], core::array::from_fn(|a| core::array::from_fn(|b| total[4 + 3 * a + b])))
}

/// The edge and face sums of [`polyhedron`], at any distance.
fn edge_face_sums(p: [f64; 3], verts: &[[f64; 3]], faces: &[[usize; 3]]) -> Newtonian {
    let sub = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let cross = |a: [f64; 3], b: [f64; 3]| [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
    let unit = |a: [f64; 3]| {
        let n = dot(a, a).sqrt();
        [a[0] / n, a[1] / n, a[2] / n]
    };
    let (mut v, mut g, mut t) = (0.0, [0.0; 3], [[0.0; 3]; 3]);
    // Edges through p, keyed by their vertices, with their faces' summed dyads.
    type Dyads = ((usize, usize), [[f64; 3]; 3]);
    let mut on_edge: Vec<Dyads> = Vec::new();
    for face in faces {
        let corners = face.map(|i| verts[i]);
        let n = unit(cross(sub(corners[1], corners[0]), sub(corners[2], corners[0])));
        let r = corners.map(|c| sub(c, p));
        let l = r.map(|r| dot(r, r).sqrt());
        // Van Oosterom & Strackee's solid angle of a triangle. In its plane, the limit from outside
        // (the prism's convention) is minus the in-plane winding angle: −2π inside, −π on an edge,
        // minus the corner's angle at a vertex, 0 outside. Edges through p count for nothing.
        let triple = dot(r[0], cross(r[1], r[2]));
        let omega = if triple == 0.0 {
            -(0..3)
                .map(|k| {
                    let (a, b) = (r[k], r[(k + 1) % 3]);
                    let (sin, cos) = (dot(cross(a, b), n), dot(a, b));
                    if l[k] == 0.0 || l[(k + 1) % 3] == 0.0 || (sin == 0.0 && cos < 0.0) { 0.0 } else { sin.atan2(cos) }
                })
                .sum::<f64>()
        } else {
            2.0 * triple.atan2(l[0] * l[1] * l[2] + l[0] * dot(r[1], r[2]) + l[1] * dot(r[2], r[0]) + l[2] * dot(r[0], r[1]))
        };
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
            if le.is_finite() {
                for a in 0..3 {
                    for b in 0..3 {
                        t[a][b] += n[a] * ne[b] * le;
                    }
                }
            } else {
                let key = (face[i].min(face[j]), face[i].max(face[j]));
                let dyad: [[f64; 3]; 3] = core::array::from_fn(|a| core::array::from_fn(|b| n[a] * ne[b]));
                match on_edge.iter_mut().find(|(k, _)| *k == key) {
                    Some((_, sum)) => (0..9).for_each(|q| sum[q / 3][q % 3] += dyad[q / 3][q % 3]),
                    None => on_edge.push((key, dyad)),
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
    // On a real edge E_e L_e is infinite; between coplanar faces E_e vanishes and so does the term.
    for (_, e) in on_edge {
        for (a, row) in e.iter().enumerate() {
            for (b, &x) in row.iter().enumerate() {
                if x.abs() > 1e-12 {
                    t[a][b] += x * f64::INFINITY;
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

/// What [`check_mesh`] found wrong with a polyhedron's faces.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MeshError {
    /// A face with zero area, whose normal is undefined.
    Degenerate(usize),
    /// A directed edge used by two faces: one of them is wound the wrong way (or more than two
    /// faces meet at the edge).
    Flipped(usize, usize),
    /// An edge with no face on its other side: the mesh is not closed.
    Open(usize, usize),
    /// A closed mesh with no positive volume: it is wound clockwise, inside out.
    Inverted,
}

/// Checks that `faces` form a closed, consistently wound, outward-facing surface, as
/// [`polyhedron`] needs, and returns its volume. Every directed edge must appear exactly once and
/// its reverse once, and the signed volume `Σ v₀·(v₁ × v₂)/6` must be positive. Without the check
/// one flipped face goes unnoticed: on a cube it makes V 13.7% wrong inside.
pub fn check_mesh(verts: &[[f64; 3]], faces: &[[usize; 3]]) -> Result<f64, MeshError> {
    let mut edges: Vec<(usize, usize)> = faces.iter().flat_map(|f| [(f[0], f[1]), (f[1], f[2]), (f[2], f[0])]).collect();
    edges.sort_unstable();
    if let Some(w) = edges.windows(2).find(|w| w[0] == w[1]) {
        return Err(MeshError::Flipped(w[0].0, w[0].1));
    }
    if let Some(&(a, b)) = edges.iter().find(|&&(a, b)| edges.binary_search(&(b, a)).is_err()) {
        return Err(MeshError::Open(a, b));
    }
    let mut volume = 0.0;
    for (k, f) in faces.iter().enumerate() {
        let [a, b, c] = f.map(|i| verts[i]);
        let (u, w) = ([b[0] - a[0], b[1] - a[1], b[2] - a[2]], [c[0] - a[0], c[1] - a[1], c[2] - a[2]]);
        let n = [u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]];
        if n == [0.0; 3] {
            return Err(MeshError::Degenerate(k));
        }
        volume += (a[0] * n[0] + a[1] * n[1] + a[2] * n[2]) / 6.0;
    }
    if volume > 0.0 { Ok(volume) } else { Err(MeshError::Inverted) }
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
    fn checking_a_mesh() {
        let (verts, faces) = cube_mesh();
        assert!((check_mesh(&verts, &faces).unwrap() - 1.0).abs() < 1e-15);
        let inverted: Vec<[usize; 3]> = faces.iter().map(|&[a, b, c]| [a, c, b]).collect();
        assert_eq!(check_mesh(&verts, &inverted), Err(MeshError::Inverted));
        let mut flipped = faces.clone();
        flipped[5] = [flipped[5][0], flipped[5][2], flipped[5][1]];
        assert!(matches!(check_mesh(&verts, &flipped), Err(MeshError::Flipped(..))));
        assert!(matches!(check_mesh(&verts, &faces[1..]), Err(MeshError::Open(..))));
        let mut flat = faces.clone();
        flat[0] = [0, 0, 1];
        assert!(matches!(check_mesh(&verts, &flat), Err(MeshError::Degenerate(0)) | Err(MeshError::Flipped(..))));
        assert_eq!(check_mesh(&[[0.0; 3], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]], &[[0, 1, 2]]), Err(MeshError::Open(0, 1)));
    }

    #[test]
    fn a_face_diagonal_has_a_finite_hessian() {
        // On the triangulated +x face's inner diagonal y = z the two triangles' edge dyads cancel:
        // the Hessian is the face's outside limit, as the prism gives, not NaN.
        let (verts, faces) = cube_mesh();
        for p in [[0.5, 0.2, 0.2], [0.5, -0.3, -0.3]] {
            let (t, tp) = (polyhedron(p, &verts, &faces).2, prism(p, MIN, MAX).2);
            assert!((0..9).all(|k| (t[k / 3][k % 3] - tp[k / 3][k % 3]).abs() < 1e-14), "{p:?}: {t:?}");
        }
        // On a real edge only the edge's own off-diagonal entries are infinite.
        let t = polyhedron([0.5, 0.5, 0.1], &verts, &faces).2;
        assert!(t[0][1].is_infinite() && t.iter().flatten().all(|x| !x.is_nan()));
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
    fn the_polyhedrons_far_field() {
        // Across the hand-over at 5 circumradii the cube mesh's two methods agree with the prism
        // (now at round-off far off) and each other; at 100 and 1000 sizes the sums alone lose
        // 2×10⁻¹⁰ and 3×10⁻⁷.
        let (verts, faces) = cube_mesh();
        let dirs = [[0.36, 0.48, 0.8], [-0.6, 0.0, 0.8], [0.8, -0.6, 0.0], [1.0, 0.0, 0.0]];
        for k in [4.9, 5.1, 10.0, 100.0, 1000.0] {
            for dir in dirs {
                let p: [f64; 3] = dir.map(|d| d * k * 0.75f64.sqrt());
                let ((v, g, t), (vp, gp, tp)) = (polyhedron(p, &verts, &faces), prism(p, MIN, MAX));
                assert!((v / vp - 1.0).abs() < 1e-13, "V at {k}, {dir:?}: {v} vs {vp}");
                let (gs, ts) = (gp.iter().fold(0f64, |m, x| m.max(x.abs())), tp.iter().flatten().fold(0f64, |m, x| m.max(x.abs())));
                assert!((0..3).all(|a| (g[a] - gp[a]).abs() < 1e-13 * gs), "g at {k}, {dir:?}: {g:?} vs {gp:?}");
                assert!((0..9).all(|q| (t[q / 3][q % 3] - tp[q / 3][q % 3]).abs() < 1e-12 * ts), "T at {k}, {dir:?}");
            }
        }
        let p = [360.0, 480.0, 800.0];
        let loss = (edge_face_sums(p, &verts, &faces).0 / prism(p, MIN, MAX).0 - 1.0).abs();
        assert!(loss > 1e-11 && loss < 1e-6, "{loss:e}");
        // A dented cube (one corner pushed inside) is not convex: some of its tetrahedra from the
        // centroid have negative volume, and the signed cubature still equals the sums at 5 and 8
        // circumradii, where the sums hold 10⁻¹³. Far off it is a point mass of its volume.
        let mut dented = verts.clone();
        dented[7] = [-0.3, -0.3, -0.3];
        let volume = check_mesh(&dented, &faces).unwrap();
        assert!((volume - 0.2).abs() < 1e-12, "{volume}"); // the dent takes four fifths
        for k in [5.0, 8.0] {
            for dir in dirs {
                let p: [f64; 3] = dir.map(|d| d * k);
                let ((v, g, _), (vs, gs, _)) = (polyhedron(p, &dented, &faces), edge_face_sums(p, &dented, &faces));
                assert!((v / vs - 1.0).abs() < 1e-12 && (0..3).all(|a| (g[a] - gs[a]).abs() < 1e-11 * g[0].abs().max(g[2].abs())), "{k} {dir:?}: {v} vs {vs} {g:?} {gs:?}");
            }
        }
        // Far off it is a point mass of its volume to its dipole term (the dent moves the centre
        // of mass by a few tenths), and the two methods agree to the sums' loss there.
        let far = [600.0, 0.0, 800.0];
        let (v, vs) = (polyhedron(far, &dented, &faces).0, edge_face_sums(far, &dented, &faces).0);
        assert!((v * 1000.0 / volume - 1.0).abs() < 1e-3 && (v / vs - 1.0).abs() < 1e-6, "{v} {vs}");
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
        // A 4×4×4 rule against the shared nodes written out.
        let (x, w) = crate::gauss_legendre(4);
        let ours = cubature([30.0, 0.0, 0.0], MIN, MAX, [4; 3]);
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
    fn fukushimas_table_4() {
        // Fukushima (2019), Table 4: half-sides (2.5, 3.1, 5.0), the point (30, 40, 50) from the
        // centre, Gρ = 1. The corner sums are already 5×10⁻¹³ off here, at 11 Brillouin radii.
        let (min, max) = ([-2.5, -3.1, -5.0], [2.5, 3.1, 5.0]);
        let p = [30.0, 40.0, 50.0];
        let (v, g, t) = prism(p, min, max);
        assert!((v / 4.385_413_035_403_429 - 1.0).abs() < 1e-14, "{v}");
        for (got, want) in g.iter().zip([-2.638_376_163_531_892e-2, -3.515_460_560_153_808e-2, -4.380_838_400_289_37e-2]) {
            assert!((got / want - 1.0).abs() < 1e-14, "{got}");
        }
        assert!((t[0][0] / -4.019_417_444_564_296e-4 - 1.0).abs() < 1e-13, "{}", t[0][0]);
        assert!((corner_sums(p, min, max).0 / 4.385_413_035_403_429 - 1.0).abs() > 1e-13);
        // Order 0 is M/R, and the quadrupole coefficients are (2a² − b² − c²)/6 and its permutations:
        // the monopole plus quadrupole leaves O((h/R)⁴).
        let (m, r) = (310.0, 5000f64.sqrt());
        assert!((m / r - 4.384_062_043_356_595).abs() < 1e-15);
        let c = [-3.685, -2.005, 5.69];
        let quad = m / r + m * (c[0] * 900.0 + c[1] * 1600.0 + c[2] * 2500.0) / r.powi(5);
        assert!((quad / v - 1.0).abs() < 3e-5 && (quad / v - 1.0).abs() > 1e-8, "{:e}", quad / v - 1.0);
    }

    #[test]
    fn needles_and_plates() {
        // The notes' regression values (8- and 10-point split Gauss–Legendre, which agree to
        // 7×10⁻¹⁶), where the corner sums lose up to 4×10⁻¹¹.
        let needle = ([-0.5, -0.5, -50.0], [0.5, 0.5, 50.0]);
        let plate = ([-50.0, -50.0, -0.5], [50.0, 50.0, 0.5]);
        type Case = (([f64; 3], [f64; 3]), [f64; 3], f64, [f64; 3], f64, f64);
        let cases: [Case; 3] = [
            (needle, [300.0, 400.0, 1200.0], 7.695_259_993_072_591e-2, [-1.368_795_824_936_556_3e-5, -1.825_061_099_915_412_4e-5, -5.467_086_502_352_132e-5], 7.092_939_317_103_203e-8, 9.755_357_675_874_793e-9),
            (needle, [30.0, 40.0, 120.0], 7.998_939_225_898_392e-1, [-1.747_532_288_963_430_4e-3, -2.330_043_052_150_337e-3, -5.981_386_012_197_591e-3], 7.941_055_859_745_535e-5, 1.780_396_513_786_331e-5),
            (plate, [120.0, -90.0, 300.0], 2.966_034_694_225_653_5e1, [-3.076_382_506_458_516_5e-2, 2.307_149_084_728_148e-2, -7.861_952_504_954_749e-2], 3.667_048_272_270_566e-4, -7.013_391_218_529_682e-5),
        ];
        for ((min, max), p, v, g, tzz, txy) in cases {
            let (ov, og, ot) = prism(p, min, max);
            assert!((ov / v - 1.0).abs() < 3e-14, "V at {p:?}: {ov}");
            assert!((0..3).all(|a| (og[a] / g[a] - 1.0).abs() < 3e-14), "g at {p:?}: {og:?}");
            assert!((ot[2][2] / tzz - 1.0).abs() < 3e-14 && (ot[0][1] / txy - 1.0).abs() < 3e-14, "T at {p:?}");
        }
        // Near a long side the needle is bisected: the whole equals its halves at points off the
        // split plane, and the halves' corner sums are exact enough there.
        let (min, max) = needle;
        for p in [[0.8, 0.1, 3.0], [0.3, 0.5 + 0.01, -20.0], [1.5, 1.5, 49.0]] {
            let (v, g, t) = prism(p, min, max);
            let (v1, g1, t1) = prism(p, min, [0.5, 0.5, 0.0]);
            let (v2, g2, t2) = prism(p, [-0.5, -0.5, 0.0], max);
            assert!((v / (v1 + v2) - 1.0).abs() < 1e-13 && (0..3).all(|a| (g[a] / (g1[a] + g2[a]) - 1.0).abs() < 1e-12), "{p:?}");
            assert!((t[2][2] / (t1[2][2] + t2[2][2]) - 1.0).abs() < 1e-11, "{p:?}");
        }
        // A point on the needle's side, with no split possible, still gets the corner sums.
        assert!((prism([0.5, 0.0, 0.0], min, max).0 - corner_sums([0.5, 0.0, 0.0], min, max).0).abs() < 1e-15);
    }

    #[test]
    fn far_field() {
        // Across the hand-over the corner sums, the adaptive rule and a large fixed rule agree, so
        // the field is continuous to 10⁻¹³.
        let (min, max) = ([0.3, -2.0, 5.0], [1.3, -1.0, 15.0]);
        let diagonal = 102f64.sqrt();
        let centre = [0.8, -1.5, 10.0];
        let dir = [0.36, 0.48, 0.8];
        // The reference is a 12³ rule on each of the ten unit cubes the box splits into, which is
        // at round-off from 3.6 units away.
        let reference = |p: [f64; 3]| {
            (0..10).fold((0.0, [0.0; 3], [[0.0; 3]; 3]), |(v, g, t), k| {
                let (dv, dg, dt) = cubature(p, [0.3, -2.0, 5.0 + k as f64], [1.3, -1.0, 6.0 + k as f64], [12; 3]);
                (v + dv, core::array::from_fn(|a| g[a] + dg[a]), core::array::from_fn(|a| core::array::from_fn(|b| t[a][b] + dt[a][b])))
            })
        };
        for k in [0.7, 1.5, 3.0, 9.99, 30.0] {
            let p: [f64; 3] = core::array::from_fn(|a| centre[a] + k * diagonal * dir[a]);
            let ((v, g, t), (vc, gc, tc)) = (prism(p, min, max), reference(p));
            assert!((v / vc - 1.0).abs() < 1e-13, "V at {k}: {v} vs {vc}");
            assert!((0..3).all(|a| (g[a] / gc[a] - 1.0).abs() < 1e-13), "g at {k}");
            assert!((t[2][2] / tc[2][2] - 1.0).abs() < 1e-12, "T at {k}");
        }
        // Far off, a cube is a point mass to its hexadecapole, (a/r)⁴.
        let p = [600.0, 0.0, 800.0];
        let (v, g, t) = prism(p, MIN, MAX);
        assert!((v * 1000.0 - 1.0).abs() < 1e-12 && (g[2] / (-800.0 / 1e9) - 1.0).abs() < 1e-12);
        assert!((t[2][2] / ((3.0 * 800.0 * 800.0 - 1e6) / 1e15) - 1.0).abs() < 1e-12);
    }
}
