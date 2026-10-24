//! The electric field of coaxial ring electrodes: coils, wires and gaps modelled as stacks of
//! charged rings, with the charge on each ring solved from the potentials the rig applies.
//!
//! A coil's turns are close to rings, and a straight wire is a stack of rings of its own radius,
//! so most high-voltage "gravity" rigs (two coils across a chamber, a wire with a gap) are
//! axially symmetric enough for this. The potential and field of one ring come in closed form
//! from the complete elliptic integrals `K(m)` and `E(m)`; the charges follow from a linear solve
//! with each ring's self-potential taken as that of a thin torus, `Q ln(8R/a)/(4π²ε₀R)` for wire
//! radius `a` (valid for `a ≪ R`).
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

/// One ring of charge: its radius (m) and its height on the common axis (m).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Ring {
    pub radius: f64,
    pub z: f64,
}

/// `n` rings of radius `radius` (m), the first at height `z0` (m) and each next one `step` (m)
/// further along: the turns of a coil, or a straight wire when `radius` is the wire's own.
pub fn stack(radius: f64, z0: f64, step: f64, n: usize) -> Vec<Ring> {
    (0..n).map(|i| Ring { radius, z: z0 + step * i as f64 }).collect()
}

/// The charge (C) on each ring when ring `i` is held at `potentials[i]` (V) against ground at
/// infinity, for conductors of wire radius `wire` (m).
pub fn charges(rings: &[Ring], wire: f64, potentials: &[f64]) -> Vec<f64> {
    let n = rings.len();
    assert_eq!(n, potentials.len(), "one potential per ring");
    // Potential coefficients P[i][j]: volts at ring i per coulomb on ring j. The last column
    // holds the right-hand side during elimination.
    let mut p: Vec<Vec<f64>> = rings
        .iter()
        .zip(potentials)
        .map(|(ri, &v)| {
            let mut row: Vec<f64> = rings
                .iter()
                .map(|rj| if ri == rj { (8.0 * ri.radius / wire).ln() / (4.0 * PI * PI * EPS_0 * ri.radius) } else { ring_potential(1.0, *rj, ri.radius, ri.z) })
                .collect();
            row.push(v);
            row
        })
        .collect();
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

/// The potential (V) at distance `rho` (m) from the axis and height `z` (m) of rings carrying
/// `charges` (C).
pub fn potential(rings: &[Ring], charges: &[f64], rho: f64, z: f64) -> f64 {
    rings.iter().zip(charges).map(|(&r, &q)| ring_potential(q, r, rho, z)).sum()
}

/// The field `[E_ρ, E_z]` (V/m) at distance `rho` (m) from the axis and height `z` (m) of rings
/// carrying `charges` (C). The field has no azimuthal part.
pub fn field(rings: &[Ring], charges: &[f64], rho: f64, z: f64) -> [f64; 2] {
    rings.iter().zip(charges).fold([0.0, 0.0], |[er, ez], (&r, &q)| {
        let [a, b] = ring_field(q, r, rho, z);
        [er + a, ez + b]
    })
}

fn ring_potential(q: f64, ring: Ring, rho: f64, z: f64) -> f64 {
    let s = (ring.radius + rho).powi(2) + (z - ring.z).powi(2);
    let (k, _) = elliptic(4.0 * ring.radius * rho / s);
    q / (4.0 * PI * EPS_0) * 2.0 / PI * k / s.sqrt()
}

// Zypman, Am. J. Phys. 74, 295 (2006), checked against direct summation in the tests.
fn ring_field(q: f64, ring: Ring, rho: f64, z: f64) -> [f64; 2] {
    let (r, dz) = (ring.radius, z - ring.z);
    let s = (r + rho).powi(2) + dz * dz;
    let d = (r - rho).powi(2) + dz * dz;
    let (k, e) = elliptic(4.0 * r * rho / s);
    let c = q / (4.0 * PI * EPS_0);
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
        // Driving the coils at ±25 kV instead makes the rig symmetric: the field at the ball is
        // the same size, and its gradient there vanishes.
        let (rings, q) = patent_rig(25e3, -25e3);
        let e2 = |z: f64| {
            let [a, b] = field(&rings, &q, 0.0, z);
            a * a + b * b
        };
        assert!(((e2(1e-4) - e2(-1e-4)) / 2e-4).abs() < 1e-6 * grad);
    }
}
