//! The electric field of high-voltage rigs: coil turns as charged rings about the z axis, and
//! wires, leads and gaps as straight charged segments anywhere, with the charge on each piece
//! solved from the potentials the rig applies.
//!
//! The potential and field of a ring come in closed form from the complete elliptic integrals
//! `K(m)` and `E(m)`, and those of a uniformly charged segment from logarithms. The charges follow
//! from a linear solve (collocation): each piece is held at its potential at one point (a
//! segment's midpoint) or on average around it (a ring, against the segments). A piece's own
//! potential is that of a thin conductor of wire radius `a`: `Q ln(8R/a)/(4π²ε₀R)` for a ring of
//! radius `R`, and `(Q/L)·2 asinh(L/2a)/(4πε₀)` for a segment of length `L`, both valid for
//! `a` much smaller than the piece. Ground is at infinity; dielectrics (insulation, glass) are
//! left out.
//!
//! ```
//! use gravitics::electrodes::{charges, field, stack, Ring};
//!
//! // US 12,417,870: two 12-turn coils of 3.5 cm radius above and below a 10 cm vacuum tank,
//! // copper 1.5 mm across, the upper coil at 50 kV and the lower one grounded.
//! let mut rings: Vec<Ring> = stack(0.035, 0.065, 0.01, 12);
//! rings.extend(stack(0.035, -0.065, -0.01, 12));
//! let volts: Vec<f64> = (0..24).map(|i| if i < 12 { 50e3 } else { 0.0 }).collect();
//! let q = charges(&rings, 0.75e-3, &volts);
//!
//! // At the ball, in the middle of the tank, the field points down at 2.5×10⁵ V/m.
//! let [e_rho, e_z] = field(&rings, &q, 0.0, 0.0);
//! assert!(e_rho == 0.0 && (e_z + 2.50e5).abs() < 1e3);
//! ```

use core::f64::consts::PI;

use metrology::EPS_0;

const K: f64 = 1.0 / (4.0 * PI * EPS_0);

/// One ring of charge about the z axis: its radius (m) and its height (m).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Ring {
    pub radius: f64,
    pub z: f64,
}

/// A straight piece of wire from `a` to `b` (m), its charge spread evenly along it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Segment {
    pub a: [f64; 3],
    pub b: [f64; 3],
}

/// A piece of an electrode.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Element {
    Ring(Ring),
    Segment(Segment),
}

/// `n` rings of radius `radius` (m), the first at height `z0` (m) and each next one `step` (m)
/// further along: the turns of a coil, or a straight wire when `radius` is the wire's own.
pub fn stack(radius: f64, z0: f64, step: f64, n: usize) -> Vec<Ring> {
    (0..n).map(|i| Ring { radius, z: z0 + step * i as f64 }).collect()
}

/// A wire along the polyline `points` (m), cut into segments no longer than `step` (m).
pub fn polyline(points: &[[f64; 3]], step: f64) -> Vec<Segment> {
    points
        .windows(2)
        .flat_map(|w| {
            let n = (dist(w[0], w[1]) / step).ceil().max(1.0) as usize;
            let at = move |i: usize| lerp(w[0], w[1], i as f64 / n as f64);
            (0..n).map(move |i| Segment { a: at(i), b: at(i + 1) })
        })
        .collect()
}

/// The charge (C) on each element when element `i` is held at `potentials[i]` (V) against ground
/// at infinity, for conductors of wire radius `wire` (m).
pub fn solve(elements: &[Element], wire: f64, potentials: &[f64]) -> Vec<f64> {
    assert_eq!(elements.len(), potentials.len(), "one potential per element");
    // Potential coefficients: volts at element i per coulomb on element j, with the applied
    // potential appended as the right-hand side.
    let rows = elements
        .iter()
        .zip(potentials)
        .map(|(ei, &v)| {
            let mut row: Vec<f64> = elements.iter().map(|ej| coefficient(ei, ej, wire)).collect();
            row.push(v);
            row
        })
        .collect();
    gauss(rows)
}

/// The potential (V) at the point `p` (m) of elements carrying `charges` (C).
pub fn potential_at(elements: &[Element], charges: &[f64], p: [f64; 3]) -> f64 {
    elements.iter().zip(charges).map(|(e, &q)| q * unit_potential(e, p)).sum()
}

/// The field `[E_x, E_y, E_z]` (V/m) at the point `p` (m) of elements carrying `charges` (C).
pub fn field_at(elements: &[Element], charges: &[f64], p: [f64; 3]) -> [f64; 3] {
    elements.iter().zip(charges).fold([0.0; 3], |acc, (e, &q)| {
        let f = unit_field(e, p);
        [acc[0] + q * f[0], acc[1] + q * f[1], acc[2] + q * f[2]]
    })
}

/// The charge (C) on each ring of a rig made only of rings: [`solve`] for rings.
pub fn charges(rings: &[Ring], wire: f64, potentials: &[f64]) -> Vec<f64> {
    solve(&as_elements(rings), wire, potentials)
}

/// The potential (V) at distance `rho` (m) from the axis and height `z` (m) of rings carrying
/// `charges` (C).
pub fn potential(rings: &[Ring], charges: &[f64], rho: f64, z: f64) -> f64 {
    potential_at(&as_elements(rings), charges, [rho, 0.0, z])
}

/// The field `[E_ρ, E_z]` (V/m) at distance `rho` (m) from the axis and height `z` (m) of rings
/// carrying `charges` (C). The field of rings has no azimuthal part.
pub fn field(rings: &[Ring], charges: &[f64], rho: f64, z: f64) -> [f64; 2] {
    let [e_rho, _, e_z] = field_at(&as_elements(rings), charges, [rho, 0.0, z]);
    [e_rho, e_z]
}

fn as_elements(rings: &[Ring]) -> Vec<Element> {
    rings.iter().map(|&r| Element::Ring(r)).collect()
}

// Volts at element i's collocation point per coulomb on element j.
fn coefficient(ei: &Element, ej: &Element, wire: f64) -> f64 {
    match (ei, ej) {
        _ if ei == ej => match ei {
            Element::Ring(r) => (8.0 * r.radius / wire).ln() / (4.0 * PI * PI * EPS_0 * r.radius),
            Element::Segment(s) => {
                let l = dist(s.a, s.b);
                K / l * 2.0 * (l / (2.0 * wire)).asinh()
            }
        },
        // A ring's potential is the same all the way round another ring, so one point does.
        (Element::Ring(r), Element::Ring(_)) => unit_potential(ej, [r.radius, 0.0, r.z]),
        (Element::Ring(r), Element::Segment(_)) => {
            let n = 16;
            (0..n)
                .map(|k| {
                    let t = 2.0 * PI * k as f64 / n as f64;
                    unit_potential(ej, [r.radius * t.cos(), r.radius * t.sin(), r.z])
                })
                .sum::<f64>()
                / n as f64
        }
        (Element::Segment(s), _) => unit_potential(ej, lerp(s.a, s.b, 0.5)),
    }
}

fn unit_potential(e: &Element, p: [f64; 3]) -> f64 {
    match e {
        Element::Ring(r) => ring_potential(1.0, *r, p[0].hypot(p[1]), p[2]),
        Element::Segment(s) => {
            let (l, r1, r2) = (dist(s.a, s.b), dist(p, s.a), dist(p, s.b));
            K / l * ((r1 + r2 + l) / (r1 + r2 - l)).ln()
        }
    }
}

fn unit_field(e: &Element, p: [f64; 3]) -> [f64; 3] {
    match e {
        Element::Ring(r) => {
            let rho = p[0].hypot(p[1]);
            let [e_rho, e_z] = ring_field(1.0, *r, rho, p[2]);
            if rho == 0.0 { [0.0, 0.0, e_z] } else { [e_rho * p[0] / rho, e_rho * p[1] / rho, e_z] }
        }
        Element::Segment(s) => segment_field(*s, p),
    }
}

// A uniform line charge of 1 C from a to b. With s measured along the segment from the foot of
// the perpendicular from p, the field is K/L·[(1/r₂ − 1/r₁) along it + (s₂/r₂ − s₁/r₁)/d across].
fn segment_field(s: Segment, p: [f64; 3]) -> [f64; 3] {
    let l = dist(s.a, s.b);
    let u = [(s.b[0] - s.a[0]) / l, (s.b[1] - s.a[1]) / l, (s.b[2] - s.a[2]) / l];
    let w = [p[0] - s.a[0], p[1] - s.a[1], p[2] - s.a[2]];
    let along = w[0] * u[0] + w[1] * u[1] + w[2] * u[2];
    let perp = [w[0] - along * u[0], w[1] - along * u[1], w[2] - along * u[2]];
    let d = (perp[0] * perp[0] + perp[1] * perp[1] + perp[2] * perp[2]).sqrt();
    let (s1, s2) = (-along, l - along);
    let (r1, r2) = (s1.hypot(d), s2.hypot(d));
    let par = K / l * (1.0 / r2 - 1.0 / r1);
    let across = if d > 0.0 { K / l * (s2 / r2 - s1 / r1) / d / d } else { 0.0 };
    [par * u[0] + across * perp[0], par * u[1] + across * perp[1], par * u[2] + across * perp[2]]
}

fn ring_potential(q: f64, ring: Ring, rho: f64, z: f64) -> f64 {
    let s = (ring.radius + rho).powi(2) + (z - ring.z).powi(2);
    let (k, _) = elliptic(4.0 * ring.radius * rho / s);
    q * K * 2.0 / PI * k / s.sqrt()
}

// Zypman, Am. J. Phys. 74, 295 (2006), checked against direct summation in the tests.
fn ring_field(q: f64, ring: Ring, rho: f64, z: f64) -> [f64; 2] {
    let (r, dz) = (ring.radius, z - ring.z);
    let s = (r + rho).powi(2) + dz * dz;
    let d = (r - rho).powi(2) + dz * dz;
    let (k, e) = elliptic(4.0 * r * rho / s);
    let c = q * K;
    let e_z = c * 2.0 / PI * dz * e / (s.sqrt() * d);
    let e_rho = if rho == 0.0 { 0.0 } else { c / (PI * rho * s.sqrt()) * (k - (r * r - rho * rho + dz * dz) / d * e) };
    [e_rho, e_z]
}

/// The complete elliptic integrals `(K(m), E(m))` for parameter `m = k²`, by the
/// arithmetic–geometric mean (Abramowitz & Stegun 17.6).
fn elliptic(m: f64) -> (f64, f64) {
    let (mut a, mut b, mut c) = (1.0, (1.0 - m).sqrt(), m.sqrt());
    let (mut sum, mut weight) = (0.5 * m, 0.5);
    for _ in 0..64 {
        if c <= f64::EPSILON * a {
            break;
        }
        (a, b, c) = (0.5 * (a + b), (a * b).sqrt(), 0.5 * (a - b));
        weight *= 2.0;
        sum += weight * c * c;
    }
    let k = PI / (2.0 * a);
    (k, k * (1.0 - sum))
}

// Gaussian elimination with partial pivoting on rows that carry their right-hand side last.
fn gauss(mut p: Vec<Vec<f64>>) -> Vec<f64> {
    let n = p.len();
    for col in 0..n {
        let pivot = (col..n).max_by(|&a, &b| p[a][col].abs().total_cmp(&p[b][col].abs())).unwrap();
        p.swap(col, pivot);
        let pivot_row = p[col].clone();
        for row in &mut p[col + 1..] {
            let f = row[col] / pivot_row[col];
            for (x, y) in row[col..].iter_mut().zip(&pivot_row[col..]) {
                *x -= f * y;
            }
        }
    }
    let mut q = vec![0.0; n];
    for row in (0..n).rev() {
        let s: f64 = (row + 1..n).map(|k| p[row][k] * q[k]).sum();
        q[row] = (p[row][n] - s) / p[row][row];
    }
    q
}

fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

fn lerp(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, rel: f64) -> bool {
        (a / b - 1.0).abs() < rel
    }

    fn patent_rig(top: f64, bottom: f64) -> (Vec<Ring>, Vec<f64>) {
        let mut rings = stack(0.035, 0.065, 0.01, 12);
        rings.extend(stack(0.035, -0.065, -0.01, 12));
        let volts: Vec<f64> = (0..24).map(|i| if i < 12 { top } else { bottom }).collect();
        let q = charges(&rings, 0.75e-3, &volts);
        (rings, q)
    }

    #[test]
    fn elliptic_integrals() {
        let (k, e) = elliptic(0.0);
        assert!(close(k, PI / 2.0, 1e-15) && close(e, PI / 2.0, 1e-15));
        // A&S Table 17.1 at m = 0.5.
        let (k, e) = elliptic(0.5);
        assert!(close(k, 1.854_074_677_301_372, 1e-14), "{k}");
        assert!(close(e, 1.350_643_881_047_675, 1e-14), "{e}");
    }

    #[test]
    fn one_ring_against_direct_summation() {
        let (ring, q, rho, z) = (Ring { radius: 0.04, z: 0.01 }, 1e-9, 0.023, 0.037);
        let n = 20_000;
        let (mut phi, mut er, mut ez) = (0.0, 0.0, 0.0);
        for i in 0..n {
            let t = 2.0 * PI * (i as f64 + 0.5) / n as f64;
            let (dx, dy, dz) = (rho - ring.radius * t.cos(), -ring.radius * t.sin(), z - ring.z);
            let d = (dx * dx + dy * dy + dz * dz).sqrt();
            let dq = q / n as f64 / (4.0 * PI * EPS_0);
            phi += dq / d;
            er += dq * dx / d.powi(3);
            ez += dq * dz / d.powi(3);
        }
        assert!(close(ring_potential(q, ring, rho, z), phi, 1e-12));
        let [a, b] = ring_field(q, ring, rho, z);
        assert!(close(a, er, 1e-12) && close(b, ez, 1e-12), "{a} {er} {b} {ez}");
        // On the axis the field is the textbook kQz/(R² + z²)^{3/2}.
        let [a, b] = ring_field(q, ring, 0.0, 0.05);
        let want = q / (4.0 * PI * EPS_0) * 0.04 / (0.04f64.powi(2) + 0.04f64.powi(2)).powf(1.5);
        assert!(a == 0.0 && close(b, want, 1e-14));
        // Far away the ring is a point charge.
        let [a, _] = ring_field(q, ring, 50.0, ring.z);
        assert!(close(a, q / (4.0 * PI * EPS_0 * 2500.0), 1e-6));
    }

    #[test]
    fn one_segment_against_direct_summation() {
        let seg = Segment { a: [0.01, -0.02, 0.03], b: [-0.04, 0.05, -0.01] };
        let p = [0.02, 0.01, -0.03];
        let (n, l) = (200_000, dist(seg.a, seg.b));
        let (mut phi, mut e) = (0.0, [0.0; 3]);
        for i in 0..n {
            let x = lerp(seg.a, seg.b, (i as f64 + 0.5) / n as f64);
            let w = [p[0] - x[0], p[1] - x[1], p[2] - x[2]];
            let d = dist(p, x);
            phi += K / n as f64 / d;
            for k in 0..3 {
                e[k] += K / n as f64 * w[k] / d.powi(3);
            }
        }
        let el = Element::Segment(seg);
        assert!(close(unit_potential(&el, p), phi, 1e-9), "{} {phi}", unit_potential(&el, p));
        let f = unit_field(&el, p);
        for k in 0..3 {
            assert!(close(f[k], e[k], 1e-8), "{k}: {} {}", f[k], e[k]);
        }
        // Far away it is a point charge; on its own axis beyond an end, K/(r₁r₂).
        assert!(close(unit_potential(&el, [100.0, 0.0, 0.0]), K / 100.0, 1e-3));
        let ax = Segment { a: [0.0, 0.0, 0.0], b: [0.0, 0.0, l] };
        assert!(close(segment_field(ax, [0.0, 0.0, 2.0 * l])[2], K / (2.0 * l * l), 1e-12));
    }

    #[test]
    fn a_wire_holds_its_potential() {
        // A 30 cm wire of 0.75 mm radius at 10 kV, in 1 cm segments.
        let segs = polyline(&[[0.0, 0.0, 0.0], [0.0, 0.0, 0.3]], 0.01);
        assert_eq!(segs.len(), 30);
        let els: Vec<Element> = segs.iter().map(|&s| Element::Segment(s)).collect();
        let q = solve(&els, 0.75e-3, &[10e3; 30]);
        // On its surface between two collocation points, and at its middle, the potential is back.
        let v = potential_at(&els, &q, [0.75e-3, 0.0, 0.150]);
        assert!(close(v, 10e3, 0.02), "{v}");
        // Charge crowds to the ends, as on any conductor.
        assert!(q[0] > q[15] && q[29] > q[15]);
        // The middle carries close to 2πε₀V/ln(L/a) per metre, the long-wire value.
        let lambda = q[15] / 0.01;
        assert!(close(lambda, 2.0 * PI * EPS_0 * 10e3 / (0.3f64 / 0.75e-3).ln(), 0.15), "{lambda:e}");
    }

    #[test]
    fn the_solved_charges_hold_the_potentials() {
        let (rings, q) = patent_rig(50e3, 0.0);
        // Halfway between two turns of the live coil, 0.75 mm off the wire's axis, the
        // potential is back at 50 kV to within the thin-ring approximation.
        let v = potential(&rings, &q, 0.035 + 0.75e-3, 0.065);
        assert!(close(v, 50e3, 0.01), "{v}");
        // The patent rig: 279 nC on the live coil, −58 nC induced on the grounded one, 7.0 mJ.
        let live: f64 = q[..12].iter().sum();
        let grounded: f64 = q[12..].iter().sum();
        assert!(close(live, 278.53e-9, 1e-3) && close(grounded, -57.54e-9, 1e-3), "{live:e} {grounded:e}");
        assert!(close(0.5 * live * 50e3, 6.963e-3, 1e-3));
    }

    #[test]
    fn the_field_at_the_ball() {
        let (rings, q) = patent_rig(50e3, 0.0);
        let [_, ez] = field(&rings, &q, 0.0, 0.0);
        assert!(close(ez, -2.4991e5, 1e-3), "{ez:e}");
        // The squared field grows toward the live coil: 1.3×10¹² V²/m³ at the ball.
        let e2 = |z: f64| {
            let [a, b] = field(&rings, &q, 0.0, z);
            a * a + b * b
        };
        let grad = (e2(1e-4) - e2(-1e-4)) / 2e-4;
        assert!(close(grad, 1.2753e12, 1e-3), "{grad:e}");
        // Driving the coils alone at ±25 kV makes the rig symmetric: the field at the ball is the
        // same size, and its gradient there vanishes. (The leads break this; see physica::rig.)
        let (rings, q) = patent_rig(25e3, -25e3);
        let e2 = |z: f64| {
            let [a, b] = field(&rings, &q, 0.0, z);
            a * a + b * b
        };
        assert!(((e2(1e-4) - e2(-1e-4)) / 2e-4).abs() < 1e-6 * grad);
    }
}
