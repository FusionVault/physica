//! Special functions that potential fields share: the complete elliptic integrals of rings of
//! mass, charge and current.
//!
//! `K(m)` and `E(m)` come from one arithmetic–geometric mean. Bulirsch's general complete
//! elliptic integral [`cel`] covers the combinations of them that would otherwise be computed as
//! a difference of two nearly equal numbers, the way a ring's field off its axis is.
//!
//! ```
//! use geopotential::special::{cel, complete_elliptic};
//!
//! let (k, e) = complete_elliptic(0.5);
//! assert!((k - 1.854_074_677_301_372).abs() < 1e-14 && (e - 1.350_643_881_047_675).abs() < 1e-14);
//! let kc = 0.5f64.sqrt();
//! assert!((cel(kc, 1.0, 1.0, 1.0) - k).abs() < 1e-14); // cel(k_c, 1, 1, 1) = K
//! ```

use core::f64::consts::FRAC_PI_2;

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
pub fn cel(kc: f64, p: f64, a: f64, b: f64) -> f64 {
    debug_assert!(p > 0.0, "cel is implemented for p > 0");
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

#[cfg(test)]
mod tests {
    use super::*;

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
}
