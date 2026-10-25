//! Slepian functions: the bandlimited functions on the sphere best concentrated inside a region,
//! here a polar cap (Simons, Dahlen & Wieczorek, SIAM Rev. 48, 504, 2006). They are the basis for
//! averaging a gravity field over a basin, or a magnetic field over a region, with the least
//! leakage from outside it.
//!
//! A function `g = Σ_l c_l P̄lm(cos θ) cos mλ` to degree `L` puts the fraction
//! `λ = cᵀDc / cᵀc` of its energy inside the cap `θ ≤ Θ`, where `D` is the cap's concentration
//! matrix. Its eigenvectors are the Slepian functions. The eigenvalues sum to the Shannon number
//! `(L + 1)²(1 − cos Θ)/2`, about how many of them are well concentrated.
//!
//! ```
//! use geopotential::slepian::polar_cap;
//!
//! // A 40° cap at degree 18: 42.2 functions' worth, and the best is concentrated to 1 − 10⁻⁷.
//! let cap = polar_cap(18, 40f64.to_radians());
//! let shannon: f64 = cap.iter().map(|s| if s.m == 0 { 1.0 } else { 2.0 } * s.concentration).sum();
//! assert!((shannon - 19.0 * 19.0 * (1.0 - 40f64.to_radians().cos()) / 2.0).abs() < 1e-12);
//! assert!(cap[0].concentration > 1.0 - 1e-7);
//! ```

use crate::legendre::{Legendre, gauss_legendre, idx};
use crate::linalg::symmetric_eigen;

/// One Slepian function of a polar cap.
#[derive(Clone, Debug, PartialEq)]
pub struct Slepian {
    /// Its order. An order above 0 stands for two functions, with `cos mλ` and `sin mλ`, which
    /// share the coefficients and the concentration.
    pub m: usize,
    /// The fraction of its energy inside the cap, between 0 and 1.
    pub concentration: f64,
    /// Its unit-norm coefficients of `P̄lm(cos θ)` (4π-normalized), for `l = m ..= L`.
    pub coefficients: Vec<f64>,
}

/// The Slepian functions to degree `lmax` of the cap within colatitude `theta` (radians) of the
/// pole, best concentrated first.
///
/// The concentration matrices are exact: `P̄lm P̄l′m` is a polynomial of degree `l + l′ ≤ 2L` in
/// `cos θ`, which `L + 1` Gauss–Legendre points over the cap integrate exactly. Their
/// eigenvalues cluster at 1 and 0, so their eigenvectors are poorly determined. The vectors come
/// instead from Grünbaum's tridiagonal matrix, which commutes with the concentration matrix and
/// has a simple spectrum (Simons et al. 2006, "Commuting tridiagonal matrix"), and the
/// concentrations are the Rayleigh quotients `cᵀDc`.
pub fn polar_cap(lmax: usize, theta: f64) -> Vec<Slepian> {
    let mut out = Vec::new();
    for m in 0..=lmax {
        let (d, g, k) = (concentration(m, lmax, theta), grunbaum(m, lmax, theta), lmax - m + 1);
        let (_, vectors) = symmetric_eigen(&g, k);
        for col in 0..k {
            let c: Vec<f64> = (0..k).map(|i| vectors[i * k + col]).collect();
            let concentration = (0..k).map(|i| c[i] * (0..k).map(|j| d[i * k + j] * c[j]).sum::<f64>()).sum();
            out.push(Slepian { m, concentration, coefficients: c });
        }
    }
    out.sort_by(|a, b| b.concentration.total_cmp(&a.concentration));
    out
}

/// The concentration matrix of order `m` for the cap, `D_ll′ = (1/4π)∫_cap P̄lm P̄l′m cos²mλ dΩ`,
/// over `l, l′ = m ..= L`.
fn concentration(m: usize, lmax: usize, theta: f64) -> Vec<f64> {
    let k = lmax - m + 1;
    let (x, w) = gauss_legendre(lmax + 1);
    let (x0, mut d) = (theta.cos(), vec![0.0; k * k]);
    let mut leg = Legendre::new(lmax);
    // ∫cos²mλ dλ/4π is 1/2 for m = 0 and 1/4 otherwise; the nodes map to [cos Θ, 1].
    let order = if m == 0 { 0.5 } else { 0.25 } * (1.0 - x0) / 2.0;
    for (xi, wi) in x.iter().zip(&w) {
        leg.compute((x0 + (1.0 - x0) * (xi + 1.0) / 2.0).asin());
        for i in 0..k {
            for j in 0..k {
                d[i * k + j] += order * wi * leg.p[idx(m + i, m)] * leg.p[idx(m + j, m)];
            }
        }
    }
    d
}

/// Grünbaum's tridiagonal matrix of order `m` for the cap: `G_ll = −l(l+1)cos Θ`,
/// `G_l,l+1 = [l(l+2) − L(L+2)]√(((l+1)² − m²)/((2l+1)(2l+3)))` (Simons et al. 2006).
fn grunbaum(m: usize, lmax: usize, theta: f64) -> Vec<f64> {
    let k = lmax - m + 1;
    let (lf, mf, cos) = (lmax as f64, m as f64, theta.cos());
    let mut g = vec![0.0; k * k];
    for i in 0..k {
        let l = (m + i) as f64;
        g[i * k + i] = -l * (l + 1.0) * cos;
        if i + 1 < k {
            let off = (l * (l + 2.0) - lf * (lf + 2.0)) * (((l + 1.0).powi(2) - mf * mf) / ((2.0 * l + 1.0) * (2.0 * l + 3.0))).sqrt();
            (g[i * k + i + 1], g[(i + 1) * k + i]) = (off, off);
        }
    }
    g
}

#[cfg(test)]
mod tests {
    use super::*;

    const L: usize = 18;

    fn theta() -> f64 {
        40f64.to_radians()
    }

    #[test]
    fn the_shannon_number() {
        // Σλ = (L + 1)²(1 − cos Θ)/2 = 42.229 for a 40° cap at L = 18, each order above 0 twice.
        let cap = polar_cap(L, theta());
        let sum: f64 = cap.iter().map(|s| if s.m == 0 { 1.0 } else { 2.0 } * s.concentration).sum();
        let shannon = (L as f64 + 1.0).powi(2) * (1.0 - theta().cos()) / 2.0;
        assert!((shannon - 42.229).abs() < 5e-4 && (sum / shannon - 1.0).abs() < 1e-13, "{sum}");
        // The notes' rank-ordered values, counting both functions of each order above 0:
        // λ₁…₃ round to 1.000000, λ₄₂ = 0.579, and 42 exceed one half.
        let mut all: Vec<f64> = cap.iter().flat_map(|s| vec![s.concentration; if s.m == 0 { 1 } else { 2 }]).collect();
        all.sort_by(|a, b| b.total_cmp(a));
        assert!(all[2] > 1.0 - 5e-7 && (all[41] - 0.579).abs() < 5e-4);
        assert_eq!(all.iter().filter(|&&l| l > 0.5).count(), 42);
        // On the whole sphere every function is fully concentrated.
        assert!(polar_cap(6, core::f64::consts::PI).iter().all(|s| (s.concentration - 1.0).abs() < 1e-13));
    }

    #[test]
    fn grunbaum_commutes_with_the_concentration_matrix() {
        // DG = GD, the identity that makes G's eigenvectors D's (Simons et al. 2006).
        for m in [0, 3, 11] {
            let (d, g, k) = (concentration(m, L, theta()), grunbaum(m, L, theta()), L - m + 1);
            let size = g.iter().fold(0f64, |a, x| a.max(x.abs()));
            for i in 0..k {
                for j in 0..k {
                    let (dg, gd): (f64, f64) = (0..k).map(|n| (d[i * k + n] * g[n * k + j], g[i * k + n] * d[n * k + j])).fold((0.0, 0.0), |s, x| (s.0 + x.0, s.1 + x.1));
                    assert!((dg - gd).abs() < 1e-14 * size, "m = {m}: [{i}][{j}] {dg} vs {gd}");
                }
            }
        }
    }

    #[test]
    fn each_function_is_concentrated_as_claimed() {
        // An eigenvector of D (Dc = λc), with its concentration found again by integrating g² over
        // the cap and the sphere with 200 points, far beyond the exact L + 1.
        let (x, w) = gauss_legendre(200);
        let mut leg = Legendre::new(L);
        let tables: Vec<(f64, Vec<f64>)> = x
            .iter()
            .zip(&w)
            .map(|(&xi, &wi)| {
                leg.compute(xi.asin());
                (wi, leg.p.clone())
            })
            .collect();
        let x0 = theta().cos();
        for s in polar_cap(L, theta()).iter().step_by(7) {
            let k = s.coefficients.len();
            let d = concentration(s.m, L, theta());
            for i in 0..k {
                let dc: f64 = (0..k).map(|j| d[i * k + j] * s.coefficients[j]).sum();
                assert!((dc - s.concentration * s.coefficients[i]).abs() < 1e-13);
            }
            // ∫g² dx over [cos Θ, 1] against [−1, 1], by mapping the 200 nodes onto each.
            let g = |p: &[f64]| (0..k).map(|i| s.coefficients[i] * p[idx(s.m + i, s.m)]).sum::<f64>();
            let whole: f64 = tables.iter().map(|(wi, p)| wi * g(p).powi(2)).sum();
            let inside: f64 = x
                .iter()
                .zip(&w)
                .map(|(&xi, &wi)| {
                    leg.compute((x0 + (1.0 - x0) * (xi + 1.0) / 2.0).asin());
                    wi * (1.0 - x0) / 2.0 * g(&leg.p).powi(2)
                })
                .sum();
            assert!((inside / whole - s.concentration).abs() < 1e-13, "m = {}: {} vs {}", s.m, inside / whole, s.concentration);
        }
    }
}
