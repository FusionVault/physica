//! Special functions that potential fields share: the elliptic integrals of rings of mass, charge
//! and current, and the ring itself.
//!
//! `K(m)` and `E(m)` come from one arithmetic–geometric mean. Bulirsch's general complete
//! elliptic integral [`cel`] covers the combinations of them that would otherwise be computed as
//! a difference of two nearly equal numbers, the way a ring's field off its axis is. Carlson's
//! symmetric integrals `R_F`, `R_D`, `R_J` and `R_C` give the incomplete integrals and the third
//! kind. [`ring`] is the average of `1/r` over a ring, the kernel of every axisymmetric coil,
//! electrode and mass ring.
//!
//! ```
//! use geopotential::special::{cel, complete_elliptic};
//!
//! let (k, e) = complete_elliptic(0.5);
//! assert!((k - 1.854_074_677_301_372).abs() < 1e-14 && (e - 1.350_643_881_047_675).abs() < 1e-14);
//! let kc = 0.5f64.sqrt();
//! assert!((cel(kc, 1.0, 1.0, 1.0) - k).abs() < 1e-14); // cel(k_c, 1, 1, 1) = K
//! ```

use core::f64::consts::{FRAC_PI_2, PI};

/// The complete elliptic integrals of the first and second kind, `(K(m), E(m))`, for parameter
/// `m = k²` below 1, by the arithmetic–geometric mean (Abramowitz & Stegun 17.6).
pub fn complete_elliptic(m: f64) -> (f64, f64) {
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
    let k = FRAC_PI_2 / a;
    (k, k * (1.0 - sum))
}

/// Bulirsch's general complete elliptic integral
/// `cel(k_c, p, a, b) = ∫₀^{π/2} (a cos²φ + b sin²φ) / ((cos²φ + p sin²φ)√(cos²φ + k_c² sin²φ)) dφ`,
/// for `p > 0` (Bulirsch, Numer. Math. 13, 305, 1969). `cel(k_c, 1, 1, 1) = K`,
/// `cel(k_c, 1, 1, k_c²) = E` and `cel(k_c, k_c², 1, 1) = E/k_c²`, with `k_c² = 1 − m`. It is
/// computed directly rather than from `K` and `E`, so it keeps its accuracy where the
/// numerator nearly cancels.
///
/// # Panics
/// If `p ≤ 0`.
pub fn cel(kc: f64, p: f64, a: f64, b: f64) -> f64 {
    assert!(p > 0.0, "cel is implemented for p > 0");
    let (mut kc, mut e, mut m) = (kc.abs(), kc.abs(), 1.0);
    let mut p = p.sqrt();
    let (mut a, mut b) = (a, b / p);
    for _ in 0..64 {
        let f = a;
        a += b / p;
        let g = e / p;
        b = 2.0 * (b + f * g);
        p += g;
        let previous = m;
        m += kc;
        if (previous - kc).abs() <= previous * 1e-15 {
            break;
        }
        kc = 2.0 * e.sqrt();
        e = kc * m;
    }
    FRAC_PI_2 * (a * m + b) / (m * (m + p))
}

/// Carlson's relative error bound `r` for the duplication algorithms: round-off, so they are as
/// accurate as the arithmetic.
const CARLSON_R: f64 = f64::EPSILON;

/// Carlson's `R_F(x, y, z) = ½∫₀^∞ dt/√((t+x)(t+y)(t+z))`, for `x, y, z ≥ 0` with at most one of
/// them 0, by duplication (Carlson, Numer. Algorithms 10, 13, 1995). `K(m) = R_F(0, 1−m, 1)`,
/// and `F(φ, m) = sin φ R_F(cos²φ, 1 − m sin²φ, 1)`.
pub fn carlson_rf(x: f64, y: f64, z: f64) -> f64 {
    let a0 = (x + y + z) / 3.0;
    let q = (3.0 * CARLSON_R).powf(-1.0 / 6.0) * (a0 - x).abs().max((a0 - y).abs()).max((a0 - z).abs());
    let (mut v, mut a, mut f) = ([x, y, z], a0, 1.0);
    while f * q >= a.abs() {
        let [sx, sy, sz] = v.map(f64::sqrt);
        let l = sx * sy + sx * sz + sy * sz;
        (v, a, f) = (v.map(|w| (w + l) / 4.0), (a + l) / 4.0, f / 4.0);
    }
    let (dx, dy) = ((a0 - x) * f / a, (a0 - y) * f / a);
    let dz = -dx - dy;
    let (e2, e3) = (dx * dy - dz * dz, dx * dy * dz);
    (1.0 - e2 / 10.0 + e3 / 14.0 + e2 * e2 / 24.0 - 3.0 * e2 * e3 / 44.0) / a.sqrt()
}

/// Carlson's `R_D(x, y, z) = (3/2)∫₀^∞ dt/((t+z)√((t+x)(t+y)(t+z)))`, for `x, y ≥ 0` with at most
/// one of them 0 and `z > 0`. `E(m) = R_F(0, 1−m, 1) − (m/3)R_D(0, 1−m, 1)`.
pub fn carlson_rd(x: f64, y: f64, z: f64) -> f64 {
    let a0 = (x + y + 3.0 * z) / 5.0;
    let q = (CARLSON_R / 4.0).powf(-1.0 / 6.0) * (a0 - x).abs().max((a0 - y).abs()).max((a0 - z).abs());
    let (mut v, mut a, mut f, mut sum) = ([x, y, z], a0, 1.0, 0.0);
    while f * q >= a.abs() {
        let [sx, sy, sz] = v.map(f64::sqrt);
        let l = sx * sy + sx * sz + sy * sz;
        sum += f / (sz * (v[2] + l));
        (v, a, f) = (v.map(|w| (w + l) / 4.0), (a + l) / 4.0, f / 4.0);
    }
    let (dx, dy) = ((a0 - x) * f / a, (a0 - y) * f / a);
    let dz = -(dx + dy) / 3.0;
    let (e2, e3) = (dx * dy - 6.0 * dz * dz, (3.0 * dx * dy - 8.0 * dz * dz) * dz);
    let (e4, e5) = (3.0 * (dx * dy - dz * dz) * dz * dz, dx * dy * dz.powi(3));
    f / (a * a.sqrt()) * series(e2, e3, e4, e5) + 3.0 * sum
}

/// Carlson's `R_J(x, y, z, p) = (3/2)∫₀^∞ dt/((t+p)√((t+x)(t+y)(t+z)))`, for `x, y, z ≥ 0` with at
/// most one of them 0. For `p < 0` it is the Cauchy principal value, through Carlson's Eq. 2.23a
/// with `y` the middle argument; at `p = 0` it diverges. The complete integral of the third kind
/// is `Π(n, m) = R_F(0, 1−m, 1) + (n/3)R_J(0, 1−m, 1, 1−n)`, for `n > 1` too.
pub fn carlson_rj(x: f64, y: f64, z: f64, p: f64) -> f64 {
    if p == 0.0 {
        return f64::INFINITY;
    }
    if p < 0.0 {
        let mut v = [x, y, z];
        v.sort_by(f64::total_cmp);
        let ([x, y, z], q) = (v, -p);
        let pp = y + (z - y) * (y - x) / (y + q); // ≥ y > 0
        let tail = 3.0 * (x * y * z / (x * z + pp * q)).sqrt() * carlson_rc(x * z + pp * q, pp * q);
        return ((pp - y) * carlson_rj(x, y, z, pp) - 3.0 * carlson_rf(x, y, z) + tail) / (y + q);
    }
    let a0 = (x + y + z + 2.0 * p) / 5.0;
    let delta = (p - x) * (p - y) * (p - z);
    let q = (CARLSON_R / 4.0).powf(-1.0 / 6.0) * [x, y, z, p].iter().fold(0f64, |m, w| m.max((a0 - w).abs()));
    let (mut v, mut a, mut f, mut sum) = ([x, y, z, p], a0, 1.0, 0.0);
    while f * q >= a.abs() {
        let [sx, sy, sz, sp] = v.map(f64::sqrt);
        let l = sx * sy + sx * sz + sy * sz;
        let d = (sp + sx) * (sp + sy) * (sp + sz);
        sum += f / d * carlson_rc(1.0, 1.0 + f * f * f * delta / (d * d));
        (v, a, f) = (v.map(|w| (w + l) / 4.0), (a + l) / 4.0, f / 4.0);
    }
    let [dx, dy, dz] = [x, y, z].map(|w| (a0 - w) * f / a);
    let dp = -(dx + dy + dz) / 2.0;
    let e2 = dx * dy + dx * dz + dy * dz - 3.0 * dp * dp;
    let e3 = dx * dy * dz + 2.0 * e2 * dp + 4.0 * dp.powi(3);
    let e4 = (2.0 * dx * dy * dz + e2 * dp + 3.0 * dp.powi(3)) * dp;
    let e5 = dx * dy * dz * dp * dp;
    f / (a * a.sqrt()) * series(e2, e3, e4, e5) + 6.0 * sum
}

/// Carlson's `R_C(x, y) = R_F(x, y, y)`, for `x ≥ 0`: `ln` and `atan` in one function,
/// `R_C(0, 1/4) = π` and `R_C(9/4, 2) = ln 2`. For `y < 0` it is the Cauchy principal value
/// `√(x/(x − y))R_C(x − y, −y)` (Carlson's Eq. 2.13a); at `y = 0` it diverges.
pub fn carlson_rc(x: f64, y: f64) -> f64 {
    if y == 0.0 {
        return f64::INFINITY;
    }
    if y < 0.0 {
        return (x / (x - y)).sqrt() * carlson_rc(x - y, -y);
    }
    let a0 = (x + 2.0 * y) / 3.0;
    let q = (3.0 * CARLSON_R).powf(-1.0 / 8.0) * (a0 - x).abs();
    let (mut v, mut a, mut f) = ([x, y], a0, 1.0);
    while f * q >= a.abs() {
        let l = 2.0 * v[0].sqrt() * v[1].sqrt() + v[1];
        (v, a, f) = (v.map(|w| (w + l) / 4.0), (a + l) / 4.0, f / 4.0);
    }
    let s = (y - a0) * f / a;
    (1.0 + s * s * (3.0 / 10.0 + s * (1.0 / 7.0 + s * (3.0 / 8.0 + s * (9.0 / 22.0 + s * (159.0 / 208.0 + s * 9.0 / 8.0)))))) / a.sqrt()
}

/// `φ` reduced by the period `π` (DLMF 19.2.10): the integer `j` and `sin r`, `cos r` of
/// `r = φ − jπ`, `|r| ≤ π/2`, with `φ = ±π/2` itself left whole so that the complete integrals
/// come out directly.
fn reduce(phi: f64) -> (f64, f64, f64) {
    let q = phi / PI;
    let j = if (q - q.round()).abs() == 0.5 { q.trunc() } else { q.round() };
    let (s, c) = (phi - j * PI).sin_cos();
    (j, s, c)
}

/// `1 − m sin²φ` without cancellation for `0 < m ≤ 1`, as `cos²φ + (1 − m)sin²φ` (Boost's form,
/// DLMF 19.36(i)); near `m = 1`, `φ = π/2` the plain form loses 4×10⁻⁹ in `E`.
fn delta2(s: f64, c: f64, m: f64) -> f64 {
    if m > 0.0 && m <= 1.0 { c * c + (1.0 - m) * s * s } else { 1.0 - m * s * s }
}

/// Legendre's incomplete integral of the first kind, `F(φ|m) = ∫₀^φ dθ/√(1 − m sin²θ)`, as
/// `sin φ R_F(cos²φ, 1 − m sin²φ, 1)` (DLMF 19.25.5) with the period `π` added back as `2jK(m)`
/// (19.2.10). `F(φ|0) = φ`, `F(φ|1) = artanh(sin φ)`. For `m > 1` it is real only while
/// `m sin²φ ≤ 1` (DLMF 19.7.4), and NaN beyond that or beyond `|φ| = π/2`.
pub fn elliptic_f(phi: f64, m: f64) -> f64 {
    let (j, s, c) = reduce(phi);
    let d2 = delta2(s, c, m);
    if d2 < 0.0 || (j != 0.0 && m >= 1.0) {
        return f64::NAN;
    }
    s * carlson_rf(c * c, d2, 1.0) + if j == 0.0 { 0.0 } else { 2.0 * j * carlson_rf(0.0, 1.0 - m, 1.0) }
}

/// Legendre's incomplete integral of the second kind, `E(φ|m) = ∫₀^φ √(1 − m sin²θ) dθ`, by
/// whichever of DLMF 19.25.9–11 has no cancelling terms for the `m` at hand (their stated
/// ranges `m ≤ 0`, `0 ≤ m ≤ 1` and `1 ≤ m`), with `2jE(m)` for the periods. `E(φ|0) = φ`,
/// `E(φ|1) = sin φ`; the same NaN cases as [`elliptic_f`].
pub fn elliptic_e(phi: f64, m: f64) -> f64 {
    let (j, s, c) = reduce(phi);
    let (x, y) = (c * c, delta2(s, c, m));
    if y < 0.0 || (j != 0.0 && m > 1.0) {
        return f64::NAN;
    }
    let s3 = s * s * s;
    let value = if s == 0.0 {
        0.0
    } else if m <= 0.0 {
        s * carlson_rf(x, y, 1.0) - m / 3.0 * s3 * carlson_rd(x, y, 1.0)
    } else if m == 1.0 {
        s
    } else if m < 1.0 {
        (1.0 - m) * s * carlson_rf(x, y, 1.0) + m * (1.0 - m) / 3.0 * s3 * carlson_rd(x, 1.0, y) + m * s * c / y.sqrt()
    } else {
        -(1.0 - m) / 3.0 * s3 * carlson_rd(y, 1.0, x) + s * y.sqrt() / c
    };
    value + if j == 0.0 { 0.0 } else { 2.0 * j * (carlson_rf(0.0, 1.0 - m, 1.0) - m / 3.0 * carlson_rd(0.0, 1.0 - m, 1.0)) }
}

/// Legendre's incomplete integral of the third kind,
/// `Π(n; φ|m) = ∫₀^φ dθ/((1 − n sin²θ)√(1 − m sin²θ))` (DLMF 19.2.7; A&S 17.7.1, the same sign of
/// `n`, where Carlson and GSL use `1 + n sin²θ`), as `sin φ R_F + (n/3)sin³φ R_J(cos²φ, 1 − m sin²φ,
/// 1, 1 − n sin²φ)` (DLMF 19.25.14) with `2jΠ(n|m)` for the periods. Beyond its pole, `n sin²φ > 1`,
/// it is the Cauchy principal value, through `R_J`'s; at the pole it is infinite. `Π(0; φ|m) =
/// F(φ|m)`, `Π(1; φ|0) = tan φ`, and `Π(n; π/2|m)` is the complete integral.
pub fn elliptic_pi(n: f64, phi: f64, m: f64) -> f64 {
    let (j, s, c) = reduce(phi);
    let y = delta2(s, c, m);
    if y < 0.0 || (j != 0.0 && m >= 1.0) {
        return f64::NAN;
    }
    let value = if s == 0.0 { 0.0 } else { s * carlson_rf(c * c, y, 1.0) + n / 3.0 * s * s * s * carlson_rj(c * c, y, 1.0, 1.0 - n * s * s) };
    value + if j == 0.0 { 0.0 } else { 2.0 * j * (carlson_rf(0.0, 1.0 - m, 1.0) + n / 3.0 * carlson_rj(0.0, 1.0 - m, 1.0, 1.0 - n)) }
}

/// The degree-5 series that ends `R_D` and `R_J`.
fn series(e2: f64, e3: f64, e4: f64, e5: f64) -> f64 {
    1.0 - 3.0 * e2 / 14.0 + e3 / 6.0 + 9.0 * e2 * e2 / 88.0 - 3.0 * e4 / 22.0 - 9.0 * e2 * e3 / 52.0 + 3.0 * e5 / 26.0
}

/// The average of `1/|x − x′|` over a ring of radius `radius` in the plane `z = 0`, centred on the
/// axis, seen from `(ρ, z)`, and its gradient `(∂/∂ρ, ∂/∂z)`: `(2/π)K(m)/√s` with
/// `s = (radius + ρ)² + z²` and `m = 4 radius ρ/s`. A ring of charge `Q` has the potential
/// `Q/(4πε₀)` times it, a ring of mass `M` the potential `GM` times it.
///
/// `∂/∂z` is Zypman's closed form (Am. J. Phys. 74, 295, 2006). `∂/∂ρ` is
/// `−(2/(π s^{3/2}))·cel(k_c, k_c², ρ + R, ρ − R)`, since Zypman's `K − (…)E` form loses digits as
/// `ρ²/m` near the axis (5% at `ρ = 10⁻⁷R`) and far along it. Closer to the axis than
/// `3×10⁻⁴` of the distance to the ring, where even `cel`'s integrand cancels, the near-axis
/// series `−(ρ/2)∂²/∂z² + (ρ³/16)∂⁴/∂z⁴` of the on-axis `1/√(R² + z²)` is exact to 10⁻¹⁵.
pub fn ring(radius: f64, rho: f64, z: f64) -> (f64, [f64; 2]) {
    let r = radius;
    let s = (r + rho).powi(2) + z * z;
    let d = (r - rho).powi(2) + z * z;
    let (k, e) = complete_elliptic(4.0 * r * rho / s);
    let v = 2.0 / PI * k / s.sqrt();
    let dv_dz = -2.0 / PI * z * e / (s.sqrt() * d);
    let u = r * r + z * z;
    let dv_drho = if rho == 0.0 {
        0.0
    } else if rho * rho < 9e-8 * u {
        let (z2, r2) = (z * z, r * r);
        0.5 * rho * (r2 - 2.0 * z2) / u.powf(2.5) - rho.powi(3) / 16.0 * (-24.0 * z2 * z2 + 72.0 * z2 * r2 - 9.0 * r2 * r2) / u.powf(4.5)
    } else {
        let kc = (d / s).sqrt(); // √(1 − m), without the subtraction near the wire
        -2.0 / (PI * s.powf(1.5)) * cel(kc, kc * kc, rho + r, rho - r)
    };
    (v, [dv_drho, dv_dz])
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::f64::consts::PI;

    fn close(a: f64, b: f64, rel: f64) -> bool {
        (a / b - 1.0).abs() < rel
    }

    #[test]
    fn complete_elliptic_integrals() {
        let (k, e) = complete_elliptic(0.0);
        assert!(close(k, FRAC_PI_2, 1e-15) && close(e, FRAC_PI_2, 1e-15));
        // A&S Table 17.1 at m = 0.5.
        let (k, e) = complete_elliptic(0.5);
        assert!(close(k, 1.854_074_677_301_372, 1e-14) && close(e, 1.350_643_881_047_675, 1e-14));
        // Legendre's relation, E K' + E' K − K K' = π/2, at m = 0.3.
        let ((k, e), (k1, e1)) = (complete_elliptic(0.3), complete_elliptic(0.7));
        assert!(close(e * k1 + e1 * k - k * k1, FRAC_PI_2, 1e-14));
    }

    #[test]
    fn cel_reduces_to_k_and_e() {
        for m in [1e-9, 0.1, 0.5, 0.9, 0.999_999] {
            let (k, e) = complete_elliptic(m);
            let kc = (1.0 - m).sqrt();
            assert!(close(cel(kc, 1.0, 1.0, 1.0), k, 1e-14), "K at m = {m}");
            assert!(close(cel(kc, 1.0, 1.0, kc * kc), e, 1e-14), "E at m = {m}");
            // ∫(1 − m sin²φ)^(−3/2) dφ = E/(1 − m)
            assert!(close(cel(kc, kc * kc, 1.0, 1.0), e / (1.0 - m), 1e-13), "E/(1−m) at m = {m}");
        }
    }

    #[test]
    fn cel_keeps_a_cancelling_numerator() {
        // With a = 1 + δ and b = −1 + δ the numerator cos 2φ + δ nearly cancels:
        // ∫(cos 2φ + δ)dφ = δπ/2 at k_c = 1. With δ a power of two, so that 1 ± δ are exact, cel
        // returns it exactly, down to δ = 2⁻⁵⁰.
        for k in [10, 30, 50] {
            let delta = 2f64.powi(-k);
            assert!(close(cel(1.0, 1.0, 1.0 + delta, -1.0 + delta), delta * FRAC_PI_2, 1e-15), "δ = 2^-{k}");
        }
    }

    #[test]
    fn carlsons_table() {
        // Carlson (1995), section 3.
        assert!(close(carlson_rf(1.0, 2.0, 0.0), 1.311_028_777_146_1, 1e-13));
        assert!(close(carlson_rf(0.5, 1.0, 0.0), 1.854_074_677_301_4, 1e-13));
        assert!(close(carlson_rf(2.0, 3.0, 4.0), 0.584_082_841_677_15, 1e-13));
        assert!(close(carlson_rc(0.0, 0.25), PI, 1e-15));
        assert!(close(carlson_rc(2.25, 2.0), 2f64.ln(), 1e-15));
        assert!(close(carlson_rj(0.0, 1.0, 2.0, 3.0), 0.776_886_237_785_82, 1e-13));
        assert!(close(carlson_rj(2.0, 3.0, 4.0, 5.0), 0.142_975_796_671_57, 1e-13));
        assert!(close(carlson_rd(0.0, 2.0, 1.0), 1.797_210_352_103_4, 1e-13));
        assert!(close(carlson_rd(2.0, 3.0, 4.0), 0.165_105_272_942_61, 1e-13));
        // Homogeneity and the equal-argument values: R_F = x^(−1/2), R_D = R_J = x^(−3/2).
        assert!(close(carlson_rf(4.0, 4.0, 4.0), 0.5, 1e-15) && close(carlson_rd(4.0, 4.0, 4.0), 0.125, 1e-15));
        assert!(close(carlson_rj(4.0, 4.0, 4.0, 4.0), 0.125, 1e-15) && close(carlson_rj(2.0, 3.0, 4.0, 4.0), carlson_rd(2.0, 3.0, 4.0), 1e-14));
        assert!(close(carlson_rf(8.0, 12.0, 16.0), carlson_rf(2.0, 3.0, 4.0) / 2.0, 1e-15));
    }

    #[test]
    fn carlson_principal_values() {
        // Carlson (1995), section 3: R_C(1/4, −2) = ln2/3, R_J(2, 3, 4, −0.5) and R_J(2, 3, 4, −5).
        assert!(close(carlson_rc(0.25, -2.0), 2f64.ln() / 3.0, 1e-14));
        assert_eq!(carlson_rc(0.0, -2.0), 0.0);
        assert!(close(carlson_rj(2.0, 3.0, 4.0, -0.5), 0.247_238_197_030_52, 1e-13));
        assert!(close(carlson_rj(2.0, 3.0, 4.0, -5.0), -0.127_112_300_429_64, 1e-13));
        // The argument order does not matter, and the principal value changes sign at
        // p = −1.2521959347 (Carlson prints "approximately −1.2552"; scipy and QUADPACK agree).
        assert!(close(carlson_rj(4.0, 2.0, 3.0, -0.5), carlson_rj(2.0, 3.0, 4.0, -0.5), 1e-14));
        assert!(carlson_rj(2.0, 3.0, 4.0, -1.252_195_934_7).abs() < 1e-10);
        assert!(carlson_rj(2.0, 3.0, 4.0, -1.24) > 0.0 && carlson_rj(2.0, 3.0, 4.0, -1.26) < 0.0);
        assert!(carlson_rj(2.0, 3.0, 4.0, 0.0).is_infinite() && carlson_rc(1.0, 0.0).is_infinite());
    }

    #[test]
    fn carlson_against_the_agm_and_cel() {
        // K and E from Carlson against the arithmetic–geometric mean; Π from R_J against
        // Bulirsch's cel(k_c, 1 − n, 1, 1), three independent algorithms.
        for m in [1e-6, 0.1, 0.5, 0.9, 0.999] {
            let (k, e) = complete_elliptic(m);
            let (rf, rd) = (carlson_rf(0.0, 1.0 - m, 1.0), carlson_rd(0.0, 1.0 - m, 1.0));
            assert!(close(rf, k, 2e-15) && close(rf - m / 3.0 * rd, e, 2e-15), "m = {m}");
            let kc = (1.0 - m).sqrt();
            for n in [-2.0, 0.3, 0.8] {
                let pi3 = rf + n / 3.0 * carlson_rj(0.0, 1.0 - m, 1.0, 1.0 - n);
                assert!(close(pi3, cel(kc, 1.0 - n, 1.0, 1.0), 1e-14), "Π({n}, {m})");
            }
        }
        // R_C in its two closed forms.
        assert!(close(carlson_rc(1.0, 4.0), (3f64.sqrt()).atan() / 3f64.sqrt(), 1e-15));
        assert!(close(carlson_rc(4.0, 1.0), (2.0 + 3f64.sqrt()).ln() / 3f64.sqrt(), 1e-15));
    }

    #[test]
    fn legendres_incomplete_integrals_against_abramowitz_and_stegun() {
        // A&S Tables 17.5 and 17.6 (F and E to eight decimals, modular angle α, m = sin²α) and 17.9
        // (Π to five), each within its rounding.
        let d = f64::to_radians;
        for (phi, alpha, want) in [(30.0, 45.0, 0.535_622_73), (30.0, 60.0, 0.542_229_11), (30.0, 90.0, 0.549_306_14), (20.0, 40.0, 0.351_988_69)] {
            assert!((elliptic_f(d(phi), d(alpha).sin().powi(2)) - want).abs() < 5e-9, "F({phi}\\{alpha})");
        }
        for (phi, alpha, want) in [(30.0, 45.0, 0.512_049_32), (30.0, 60.0, 0.506_092_07), (90.0, 45.0, 1.350_643_88), (70.0, 60.0, 1.026_636_89), (75.0, 85.0, 0.969_922_12), (80.0, 88.0, 0.985_689_15), (85.0, 88.0, 0.997_483_92)] {
            assert!((elliptic_e(d(phi), d(alpha).sin().powi(2)) - want).abs() < 5e-9, "E({phi}\\{alpha})");
        }
        for (n, phi, alpha, want) in [(0.5, 60.0, 45.0, 1.38218), (0.3, 75.0, 60.0, 1.93879), (0.6, 45.0, 30.0, 0.91527), (0.5, 90.0, 45.0, 2.70129), (0.1, 90.0, 0.0, 1.65576), (0.0, 75.0, 90.0, 2.02759)] {
            assert!((elliptic_pi(n, d(phi), d(alpha).sin().powi(2)) - want).abs() < 5e-6, "Π({n};{phi}\\{alpha})");
        }
    }

    #[test]
    #[allow(clippy::approx_constant)] // φ = 1.5707 is the notes' case near π/2, not π/2
    fn legendres_integrals_special_values_and_periods() {
        // DLMF 19.6: F(φ|1) = gd⁻¹φ, E(φ|1) = sin φ, F(φ|0) = E(φ|0) = φ, Π(1;φ|0) = tan φ, and
        // Π(n;φ|0) = R_C(csc²φ − 1, csc²φ − n), a principal value for n = 3.
        assert!(close(elliptic_f(0.7, 1.0), 0.7f64.sin().atanh(), 1e-15) && close(elliptic_e(0.7, 1.0), 0.7f64.sin(), 1e-15));
        assert!(close(elliptic_f(1.3, 0.0), 1.3, 1e-15) && close(elliptic_e(1.3, 0.0), 1.3, 1e-15));
        assert!(close(elliptic_pi(1.0, 1.1, 0.0), 1.1f64.tan(), 1e-14));
        let csc2 = 1.0 / 1.1f64.sin().powi(2);
        assert!(close(elliptic_pi(3.0, 1.1, 0.0), carlson_rc(csc2 - 1.0, csc2 - 3.0), 1e-14));
        // The complete integrals at φ = π/2, and DLMF 19.2.10 with the complete ones beyond it;
        // F and E are odd.
        let m = 0.7;
        let (k, e) = complete_elliptic(m);
        assert!(close(elliptic_f(FRAC_PI_2, m), k, 1e-15) && close(elliptic_e(FRAC_PI_2, m), e, 1e-15));
        for (j, r) in [(1.0, 0.4), (2.0, -0.3), (-1.0, 1.0)] {
            assert!(close(elliptic_f(j * PI + r, m), 2.0 * j * k + elliptic_f(r, m), 1e-14), "F, j = {j}");
            assert!(close(elliptic_e(j * PI + r, m), 2.0 * j * e + elliptic_e(r, m), 1e-14), "E, j = {j}");
        }
        assert!((elliptic_f(-0.8, 0.6) + elliptic_f(0.8, 0.6)).abs() < 1e-16 && (elliptic_e(-0.8, 0.6) + elliptic_e(0.8, 0.6)).abs() < 1e-16);
        assert!(close(elliptic_pi(2.0, PI + 0.3, 0.5), 2.0 * elliptic_pi(2.0, FRAC_PI_2, 0.5) + elliptic_pi(2.0, 0.3, 0.5), 1e-14));
        // m > 1 by the reciprocal modulus (DLMF 19.7.4): F(φ|m) = F(β|1/m)/√m with sin β = √m sin φ,
        // E(φ|m) = √m[E(β|1/m) − (1 − 1/m)F(β|1/m)]; real only while m sin²φ ≤ 1.
        let (phi, m): (f64, f64) = (0.5, 2.5);
        let beta = (m.sqrt() * phi.sin()).asin();
        assert!(close(elliptic_f(phi, m), elliptic_f(beta, 1.0 / m) / m.sqrt(), 1e-14) && close(elliptic_f(phi, m), 0.570_176_528_700_536, 1e-14));
        assert!(close(elliptic_e(phi, m), m.sqrt() * (elliptic_e(beta, 1.0 / m) - (1.0 - 1.0 / m) * elliptic_f(beta, 1.0 / m)), 1e-14));
        assert!(close(elliptic_e(phi, m), 0.444_842_124_053_142, 1e-14) && elliptic_f(1.0, m).is_nan() && elliptic_e(1.0, m).is_nan());
        // E near m = 1 at φ = 1.5707: 19.25.10 with the plain Δ² = 1 − m sin²φ loses 4×10⁻⁹ there;
        // with cos²φ + (1 − m)sin²φ it agrees with 19.25.9 (which has no cancellation at this m) to
        // 10⁻¹⁴, the notes' quadrature check.
        let (phi, m): (f64, f64) = (1.5707, 1.0 - 1e-10);
        let (s, c) = phi.sin_cos();
        let (x, plain) = (c * c, 1.0 - m * s * s);
        let lossy = (1.0 - m) * s * carlson_rf(x, plain, 1.0) + m * (1.0 - m) / 3.0 * s * s * s * carlson_rd(x, 1.0, plain) + m * s * c / plain.sqrt();
        let other = s * carlson_rf(x, plain, 1.0) - m / 3.0 * s * s * s * carlson_rd(x, plain, 1.0);
        assert!(close(elliptic_e(phi, m), other, 1e-14) && (lossy / other - 1.0).abs() > 1e-10, "{} {lossy} {other}", elliptic_e(phi, m));
    }

    #[test]
    fn pi_beyond_its_pole_is_the_principal_value() {
        // n sin²φ > 1: DLMF 19.25.16 computes the same principal value with no negative R_J
        // argument, through ω² = m/n; the four cases of the notes agree to 10⁻¹⁴.
        for (n, phi, m, want) in [(2.0, 1.2, 0.5, 0.253_604_270_701_507), (4.0, 0.9, 0.3, 0.259_293_626_045_840), (1.5, 1.4, 0.9, -0.973_823_998_688_881), (10.0, 0.5, 0.2, 0.233_655_250_908_165)] {
            let c = 1.0 / f64::sin(phi).powi(2);
            let w2: f64 = m / n;
            let dlmf = -w2 / 3.0 * carlson_rj(c - 1.0, c - m, c, c - w2) + ((c - 1.0) * (c - m) / ((n - 1.0) * (1.0 - w2))).sqrt() * carlson_rc(c * (n - 1.0) * (1.0 - w2), (n - c) * (c - w2));
            let ours = elliptic_pi(n, phi, m);
            assert!(close(ours, dlmf, 1e-14) && close(ours, want, 1e-14), "Π({n};{phi}|{m}) = {ours}");
        }
        // The complete principal value, DLMF 19.6.5: Π(n|m) = K − Π(m/n|m) for n > 1.
        for (n, m) in [(2.0, 0.5), (5.0, 0.9)] {
            assert!(close(elliptic_pi(n, FRAC_PI_2, m), complete_elliptic(m).0 - elliptic_pi(m / n, FRAC_PI_2, m), 1e-13));
        }
        assert!(close(elliptic_pi(2.0, FRAC_PI_2, 0.5), -0.313_544_683_465_184, 1e-13));
        // At the pole itself it diverges, and Π(0; φ|m) = F(φ|m).
        assert!(elliptic_pi(1.0, FRAC_PI_2, 0.3).is_infinite() && close(elliptic_pi(0.0, 0.9, 0.4), elliptic_f(0.9, 0.4), 1e-15));
    }

    /// The ring average of 1/r and its gradient by direct summation over `n` points, compensated
    /// (Neumaier) so the near cancellation of the two sides costs only ε·R/ρ.
    fn summed_ring(radius: f64, rho: f64, z: f64, n: usize) -> [f64; 3] {
        let mut sums = [(0.0f64, 0.0f64); 3];
        for i in 0..n {
            let t = 2.0 * PI * (i as f64 + 0.5) / n as f64;
            let d2 = rho * rho + radius * radius - 2.0 * rho * radius * t.cos() + z * z;
            let terms = [1.0 / d2.sqrt(), -(rho - radius * t.cos()) / d2.powf(1.5), -z / d2.powf(1.5)];
            for (k, &x) in terms.iter().enumerate() {
                let x = x / n as f64;
                let (sum, carry) = &mut sums[k];
                let next = *sum + x;
                *carry += if sum.abs() >= x.abs() { (*sum - next) + x } else { (x - next) + *sum };
                *sum = next;
            }
        }
        sums.map(|(s, c)| s + c)
    }

    #[test]
    fn a_ring_against_direct_summation() {
        // Either side of the switch to the near-axis series, near the wire, in its plane and far
        // along the axis (100 R).
        let r = 0.04;
        for (rho, z) in [(0.023, 0.027), (r * 3.5e-4, 0.01), (r * 2.5e-4, 0.01), (r * 1e-3, 0.01), (1e-3, 4.0), (0.3, 4.0), (0.02, 0.0), (0.0399, 1e-4)] {
            let (v, [dr, dz]) = ring(r, rho, z);
            let [sv, sr, sz] = summed_ring(r, rho, z, 20_000);
            assert!(close(v, sv, 1e-12) && close(dr, sr, 1e-10), "ρ = {rho}, z = {z}: {dr:e} vs {sr:e}");
            assert!(if z == 0.0 { dz == 0.0 } else { close(dz, sz, 1e-12) }, "∂z at ρ = {rho}, z = {z}");
        }
        // Against the near-axis −(ρ/2)∂²/∂z² of the on-axis 1/√(R² + z²), where summation can't
        // resolve it: the old Zypman form was 1.4×10⁻⁶ off at ρ = 10⁻⁵R and 5% at 10⁻⁷R.
        for f in [1e-7, 1e-10] {
            let (rho, z) = (f * r, 0.01);
            let u = r * r + z * z;
            let near = 0.5 * rho * (r * r - 2.0 * z * z) / u.powf(2.5);
            assert!(close(ring(r, rho, z).1[0], near, 1e-12), "ρ = {f}R");
        }
        // On the axis, exactly 1/√(R² + z²) and −z/(R² + z²)^(3/2).
        let (v, [dr, dz]) = ring(2.0, 0.0, 1.5);
        assert!(close(v, 0.4, 1e-15) && dr == 0.0 && close(dz, -1.5 / 15.625, 1e-15));
    }
}
