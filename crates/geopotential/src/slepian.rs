//! Slepian functions: the bandlimited functions on the sphere best concentrated inside a region
//! (Simons, Dahlen & Wieczorek, SIAM Rev. 48, 504, 2006). They are the basis for averaging a
//! gravity field over a basin, or a magnetic field over a region, with the least leakage from
//! outside it.
//!
//! A function `g = Σ c_lm P̄lm(cos θ){cos, sin}mλ` to degree `L` puts the fraction
//! `λ = cᵀDc / cᵀc` of its energy inside a region `R`, where
//! `D = (1/4π)∫_R P̄lm P̄l′m′ T_m T_m′ dΩ` is the region's concentration matrix. Its eigenvectors
//! are the Slepian functions, and its eigenvalues sum to the Shannon number `(L + 1)²A/4π`,
//! about how many of them are well concentrated. For a [`polar_cap`] the matrix is one block
//! per order and Grünbaum's commuting matrix gives clean eigenvectors; any other [`region`] is
//! integrated the way Simons' `kernelc` does it, exactly in longitude and by Gauss–Legendre in
//! colatitude, here with the nodes clustered at the region's vertex colatitudes. The region's
//! first row of `D` is its [`basin_kernel`], which averages a field over it (Swenson & Wahr
//! 2002, Eqs. 17–18).
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

use core::f64::consts::{FRAC_PI_2, PI, TAU};

use crate::legendre::{Legendre, gauss_legendre, idx};
use crate::linalg::symmetric_eigen;
use crate::{Coefficients, Functional};

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

/// One Slepian function of any region: its concentration and its coefficients, a unit vector
/// in the crate's normalization (`(1/4π)∫g² dΩ = 1`), held as a field with unit `GM` and radius
/// so that it can be synthesized and compared with any [`Coefficients`] by their dot product.
#[derive(Clone, Debug, PartialEq)]
pub struct RegionSlepian {
    /// The fraction of its energy inside the region, between 0 and 1.
    pub concentration: f64,
    /// Its coefficients `C̄lm`, `S̄lm` to the degree the region was built for.
    pub coefficients: Coefficients,
}

/// A region as the longitude intervals it covers at a colatitude: at `θ` (radians) the function
/// returns each `(λ₁, λ₂)` with `λ₁ < λ₂` inside it, and nothing where the region is absent.
/// Longitudes may run beyond `[0, 2π)`: only the differences matter.
pub type Intervals<'a> = &'a dyn Fn(f64) -> Vec<(f64, f64)>;

/// The basis the general matrices are indexed by: `(l, m, sine)` for `l ≤ L`, `m ≤ l`, the sine
/// term for `m > 0` only.
fn basis(lmax: usize) -> Vec<(usize, usize, bool)> {
    (0..=lmax).flat_map(|l| (0..=l).flat_map(move |m| if m == 0 { vec![(l, m, false)] } else { vec![(l, m, false), (l, m, true)] })).collect()
}

/// The integrals of `cos fλ` and `sin fλ` over `intervals`, summed, for `f = 0 ..= 2L`, from which
/// every product of two longitude factors follows: `cos a cos b = ½[S(a−b) + S(a+b)]`,
/// `sin a sin b = ½[S(a−b) − S(a+b)]`, `sin a cos b = ½[C(a+b) + C(a−b)]`, with `S` even and `C`
/// odd in `f`.
fn longitude_integrals(lmax: usize, intervals: &[(f64, f64)]) -> (Vec<f64>, Vec<f64>) {
    let (mut s, mut c) = (vec![0.0; 2 * lmax + 1], vec![0.0; 2 * lmax + 1]);
    for &(lo, hi) in intervals {
        s[0] += hi - lo;
        for f in 1..=2 * lmax {
            let ff = f as f64;
            s[f] += ((ff * hi).sin() - (ff * lo).sin()) / ff;
            c[f] += ((ff * lo).cos() - (ff * hi).cos()) / ff;
        }
    }
    (s, c)
}

/// `∫ T_a T_b dλ` over the intervals behind `(s, c)`, for orders `a`, `b` and their sine flags.
fn longitude_product((s, c): (&[f64], &[f64]), (a, sa): (usize, bool), (b, sb): (usize, bool)) -> f64 {
    let even = |f: isize| s[f.unsigned_abs()];
    let odd = |f: isize| if f < 0 { -c[(-f) as usize] } else { c[f as usize] };
    let (sum, diff) = ((a + b) as isize, a as isize - b as isize);
    match (sa, sb) {
        (false, false) => 0.5 * (even(diff) + even(sum)),
        (true, true) => 0.5 * (even(diff) - even(sum)),
        (true, false) => 0.5 * (odd(sum) + odd(diff)),
        (false, true) => 0.5 * (odd(sum) - odd(diff)),
    }
}

/// The nodes and weights of the colatitude integration over `[a, b]`: `nodes` Gauss–Legendre
/// points in `t ∈ [0, π]` mapped as `θ = a + (b − a)(1 − cos t)/2`, which clusters them at both
/// ends, where a region's width goes like a square root; the weights carry `sin θ`.
fn colatitude_rule(a: f64, b: f64, nodes: usize) -> Vec<(f64, f64)> {
    let (x, w) = gauss_legendre(nodes);
    x.iter()
        .zip(&w)
        .map(|(&x, &w)| {
            let t = (x + 1.0) * FRAC_PI_2;
            let theta = a + (b - a) * (1.0 - t.cos()) / 2.0;
            (theta, w * FRAC_PI_2 * (b - a) * t.sin() / 2.0 * theta.sin())
        })
        .collect()
}

/// The concentration matrix `D` of a region to degree `lmax`, row-major over the basis `(l, m,
/// cos)`, `(l, m, sin)` for `l ≤ L`, `m ≤ l` (the sine for `m > 0` only), `l` outermost:
/// `D_ij = (1/4π)∫_R P̄_i P̄_j T_i T_j dΩ`, integrated exactly in longitude and with `nodes`
/// clustered Gauss–Legendre points in colatitude on each panel between consecutive `breaks`,
/// which should include every colatitude where the region's outline has a vertex or an extremum
/// (Simons' `kernelc`, with the panels added). Its trace is the Shannon number `(L + 1)²A/4π` and
/// its first entry the fractional area `A/4π`.
pub fn concentration_matrix(lmax: usize, intervals: Intervals, breaks: &[f64], nodes: usize) -> Vec<f64> {
    let basis = basis(lmax);
    let n = basis.len();
    let mut d = vec![0.0; n * n];
    let mut leg = Legendre::new(lmax);
    for panel in breaks.windows(2) {
        for (theta, weight) in colatitude_rule(panel[0], panel[1], nodes) {
            let at = intervals(theta);
            if at.is_empty() {
                continue;
            }
            let (s, c) = longitude_integrals(lmax, &at);
            leg.compute(FRAC_PI_2 - theta);
            let p: Vec<f64> = basis.iter().map(|&(l, m, _)| leg.p[idx(l, m)]).collect();
            for i in 0..n {
                for j in 0..=i {
                    let value = weight * p[i] * p[j] * longitude_product((&s, &c), (basis[i].1, basis[i].2), (basis[j].1, basis[j].2));
                    d[i * n + j] += value;
                    if i != j {
                        d[j * n + i] += value;
                    }
                }
            }
        }
    }
    d.iter_mut().for_each(|x| *x /= 4.0 * PI);
    d
}

/// The Slepian functions to degree `lmax` of a region ([`concentration_matrix`] has the
/// arguments), best concentrated first. Their eigenvalues cluster at 1 and 0, so vectors within
/// a cluster are determined only up to rotations among themselves; the sum of the top `N`
/// functions' projections, which is what a regional estimate uses, is not affected.
pub fn region(lmax: usize, intervals: Intervals, breaks: &[f64], nodes: usize) -> Vec<RegionSlepian> {
    let basis = basis(lmax);
    let n = basis.len();
    let (values, vectors) = symmetric_eigen(&concentration_matrix(lmax, intervals, breaks, nodes), n);
    let mut out: Vec<RegionSlepian> = (0..n)
        .map(|k| {
            let mut coefficients = Coefficients::zeros(1.0, 1.0, lmax);
            for (i, &(l, m, sine)) in basis.iter().enumerate() {
                let (c, s) = (coefficients.c(l, m), coefficients.s(l, m));
                if sine { coefficients.set(l, m, c, vectors[i * n + k]) } else { coefficients.set(l, m, vectors[i * n + k], s) }
            }
            RegionSlepian { concentration: values[k], coefficients }
        })
        .collect();
    out.sort_by(|a, b| b.concentration.total_cmp(&a.concentration));
    out
}

/// The longitude intervals of a spherical cap of angular `radius` centred at colatitude `colat`
/// and longitude `lon` (radians), and the breaks its integration needs: the cap's two colatitude
/// extremes and, when it holds a pole, the colatitude where it stops covering every longitude.
pub fn cap(colat: f64, lon: f64, radius: f64) -> (impl Fn(f64) -> Vec<(f64, f64)> + use<>, Vec<f64>) {
    let intervals = move |theta: f64| {
        let q = (radius.cos() - theta.cos() * colat.cos()) / (theta.sin() * colat.sin()).max(f64::MIN_POSITIVE);
        if q <= -1.0 {
            vec![(0.0, TAU)]
        } else if q >= 1.0 {
            vec![]
        } else {
            let half = q.acos();
            vec![(lon - half, lon + half)]
        }
    };
    let mut breaks = vec![(colat - radius).max(0.0), (colat - radius).abs(), (colat + radius).min(PI), (TAU - colat - radius).min(PI)];
    breaks.sort_by(f64::total_cmp);
    breaks.dedup();
    (intervals, breaks)
}

/// The longitude intervals of a polygon with `vertices` `(λ, φ)` (longitude and latitude,
/// radians), its edges straight in `(λ, φ)` as in Simons' `dphregion`, and the breaks at its
/// vertices' colatitudes. A vertical crossing of the polygon at a colatitude pairs the edges it
/// cuts from west to east, so the polygon may be concave; it must not cross the dateline or
/// hold a pole.
pub fn polygon(vertices: &[(f64, f64)]) -> (impl Fn(f64) -> Vec<(f64, f64)> + use<>, Vec<f64>) {
    let v: Vec<(f64, f64)> = vertices.iter().map(|&(lon, lat)| (lon, FRAC_PI_2 - lat)).collect();
    let mut breaks: Vec<f64> = v.iter().map(|&(_, c)| c).collect();
    breaks.sort_by(f64::total_cmp);
    breaks.dedup();
    let intervals = move |theta: f64| {
        let mut crossings: Vec<f64> = (0..v.len())
            .filter_map(|i| {
                let ((l0, c0), (l1, c1)) = (v[i], v[(i + 1) % v.len()]);
                ((c0 - theta) * (c1 - theta) < 0.0).then(|| l0 + (l1 - l0) * (theta - c0) / (c1 - c0))
            })
            .collect();
        crossings.sort_by(f64::total_cmp);
        crossings.chunks_exact(2).map(|pair| (pair[0], pair[1])).collect()
    };
    (intervals, breaks)
}

/// The basin kernel `ϑ̂_lm = (1/4π)∫_R P̄lm {cos, sin}mλ dΩ` of a region to degree `lmax`, as
/// coefficients with unit `GM` and radius: the first row of its [`concentration_matrix`], at a
/// cost of `n` per node instead of `n²`. `ϑ̂_00` is the fractional area. Swenson & Wahr's
/// (2002) `ϑ^c_lm` is `4πϑ̂` and their `Ω_region` is `4πϑ̂_00`.
pub fn basin_kernel(lmax: usize, intervals: Intervals, breaks: &[f64], nodes: usize) -> Coefficients {
    let basis = basis(lmax);
    let mut kernel = vec![0.0; basis.len()];
    let mut leg = Legendre::new(lmax);
    for panel in breaks.windows(2) {
        for (theta, weight) in colatitude_rule(panel[0], panel[1], nodes) {
            let at = intervals(theta);
            if at.is_empty() {
                continue;
            }
            let (s, c) = longitude_integrals(lmax, &at);
            leg.compute(FRAC_PI_2 - theta);
            for (k, &(l, m, sine)) in basis.iter().enumerate() {
                kernel[k] += weight * leg.p[idx(l, m)] * longitude_product((&s, &c), (m, sine), (0, false));
            }
        }
    }
    let mut out = Coefficients::zeros(1.0, 1.0, lmax);
    for (k, &(l, m, sine)) in basis.iter().enumerate() {
        let (c, s) = (out.c(l, m), out.s(l, m));
        let v = kernel[k] / (4.0 * PI);
        if sine { out.set(l, m, c, v) } else { out.set(l, m, v, s) }
    }
    out
}

/// The mean of a `functional` of `field` over the region of `kernel` at radius `r`:
/// `Σ_lm κ_l (C̄_lm ϑ̂^c_lm + S̄_lm ϑ̂^s_lm)/ϑ̂_00`, Swenson & Wahr's (2002) Eq. 17 in the
/// crate's normalization, to the lower of the two degrees. With
/// [`Functional::WaterHeight`] it is the basin's mean water height.
pub fn basin_average(field: &Coefficients, kernel: &Coefficients, functional: &Functional, r: f64) -> f64 {
    let nmax = field.max_degree().min(kernel.max_degree());
    let mut sum = 0.0;
    for n in 0..=nmax {
        let factor = functional.factor(n, field.gm, field.radius, r);
        for m in 0..=n {
            sum += factor * (field.c(n, m) * kernel.c(n, m) + field.s(n, m) * kernel.s(n, m));
        }
    }
    sum / kernel.c(0, 0)
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

    /// The polar cap's concentrations with their multiplicity, best first.
    fn polar_spectrum() -> Vec<f64> {
        let mut all: Vec<f64> = polar_cap(L, theta()).iter().flat_map(|s| vec![s.concentration; if s.m == 0 { 1 } else { 2 }]).collect();
        all.sort_by(|a, b| b.total_cmp(a));
        all
    }

    #[test]
    fn a_cap_anywhere_has_the_polar_caps_spectrum() {
        // The general integration of the polar cap itself reproduces the block method's
        // eigenvalues, and the same cap centred at colatitude 60°, longitude 30° has them too
        // (rotation invariance), with the trace the Shannon number: 48 nodes per panel reach
        // 10⁻¹³, where 24 leave 4×10⁻³ (the notes).
        let polar = polar_spectrum();
        let shannon = (L as f64 + 1.0).powi(2) * (1.0 - theta().cos()) / 2.0;
        let (at_pole, breaks) = cap(0.0, 0.0, theta());
        let (moved, moved_breaks) = cap(60f64.to_radians(), 30f64.to_radians(), theta());
        for (name, intervals, breaks) in [("polar", &at_pole as Intervals, &breaks[..]), ("moved", &moved as Intervals, &moved_breaks[..])] {
            let d = concentration_matrix(L, intervals, breaks, 48);
            let n = (L + 1) * (L + 1);
            let trace: f64 = (0..n).map(|i| d[i * n + i]).sum();
            assert!((trace - shannon).abs() < 1e-12, "{name}: {trace}");
            let functions = region(L, intervals, breaks, 48);
            for (k, (got, want)) in functions.iter().map(|s| s.concentration).zip(&polar).enumerate() {
                assert!((got - want).abs() < 1e-12, "{name} λ{k}: {got} vs {want}");
            }
            // Unit norm, and the concentration is the Rayleigh quotient of the matrix.
            let best = &functions[0];
            let norm: f64 = (0..=L).flat_map(|l| (0..=l).map(move |m| (l, m))).map(|(l, m)| best.coefficients.c(l, m).powi(2) + best.coefficients.s(l, m).powi(2)).sum();
            assert!((norm - 1.0).abs() < 1e-13 && best.concentration > 1.0 - 1e-6);
        }
        let coarse = region(L, &moved as Intervals, &moved_breaks, 24);
        let worst = coarse.iter().map(|s| s.concentration).zip(&polar).map(|(a, b)| (a - b).abs()).fold(0f64, f64::max);
        assert!(worst > 1e-4, "{worst:e}");
    }

    #[test]
    fn polygons_have_their_areas_and_shannon_numbers() {
        // A longitude–latitude box 10–50°E, 20°S–15°N: D₀₀ = A/4π exactly, and tr D = (L+1)²A/4π.
        let d = f64::to_radians;
        let (box_intervals, breaks) = polygon(&[(d(10.0), d(-20.0)), (d(50.0), d(-20.0)), (d(50.0), d(15.0)), (d(10.0), d(15.0))]);
        let area = (d(15.0).sin() - d(-20.0).sin()) * d(40.0);
        let n = 13 * 13;
        let m = concentration_matrix(12, &box_intervals as Intervals, &breaks, 32);
        assert!((m[0] - area / (4.0 * PI)).abs() < 1e-15, "{} vs {}", m[0], area / (4.0 * PI));
        let trace: f64 = (0..n).map(|i| m[i * n + i]).sum();
        assert!((trace / (169.0 * area / (4.0 * PI)) - 1.0).abs() < 1e-13);
        // A triangle (0,0), (40,10), (15,50) with edges straight in (λ, φ): its area by a fine
        // midpoint sum of its width times cos φ, 0.2601573187 of the sphere's... as 4πD₀₀.
        let (tri, breaks) = polygon(&[(0.0, 0.0), (d(40.0), d(10.0)), (d(15.0), d(50.0))]);
        let steps = 200_000;
        let fine: f64 = (0..steps)
            .map(|i| {
                let lat = d(50.0) * (i as f64 + 0.5) / steps as f64;
                tri(FRAC_PI_2 - lat).iter().map(|(a, b)| b - a).sum::<f64>() * lat.cos() * d(50.0) / steps as f64
            })
            .sum();
        let m = concentration_matrix(12, &tri as Intervals, &breaks, 32);
        assert!((m[0] * 4.0 * PI - fine).abs() < 1e-9 && (fine - 0.260_157_318_7).abs() < 1e-9, "{} vs {fine}", m[0] * 4.0 * PI);
        let trace: f64 = (0..n).map(|i| m[i * n + i]).sum();
        assert!((trace / (169.0 * fine / (4.0 * PI)) - 1.0).abs() < 1e-11, "{}", trace / (169.0 * fine / (4.0 * PI)) - 1.0);
        // Concave polygons work: an L shape's area is the two boxes'.
        let (ell, breaks) = polygon(&[(0.0, 0.0), (d(20.0), 0.0), (d(20.0), d(10.0)), (d(10.0), d(10.0)), (d(10.0), d(20.0)), (0.0, d(20.0))]);
        let want = d(20.0) * d(10.0).sin() + d(10.0) * (d(20.0).sin() - d(10.0).sin());
        assert!((basin_kernel(4, &ell as Intervals, &breaks, 32).c(0, 0) * 4.0 * PI - want).abs() < 1e-14);
    }

    #[test]
    fn the_basin_average_of_a_harmonic() {
        // The mean of P̄₂₀ over the 40° polar cap is √5 c(1 + c)/2 with c = cos 40° (the integral
        // of P₂ over the cap in closed form): 1.512553001033005. The kernel is D's first row.
        let (at_pole, breaks) = cap(0.0, 0.0, theta());
        let kernel = basin_kernel(L, &at_pole as Intervals, &breaks, 48);
        let c = theta().cos();
        assert!((kernel.c(0, 0) - (1.0 - c) / 2.0).abs() < 1e-16);
        let mut field = Coefficients::zeros(1.0, 1.0, L);
        field.set(2, 0, 1.0, 0.0);
        let mean = basin_average(&field, &kernel, &Functional::Potential, 1.0);
        assert!((mean - 5f64.sqrt() * c * (1.0 + c) / 2.0).abs() < 1e-14, "{mean}");
        let d = concentration_matrix(L, &at_pole as Intervals, &breaks, 48);
        let n = (L + 1) * (L + 1);
        for (k, &(l, m, sine)) in basis(L).iter().enumerate() {
            let want = if sine { kernel.s(l, m) } else { kernel.c(l, m) };
            assert!((d[k] - want).abs() < 1e-16 && (d[k * n] - want).abs() < 1e-16);
        }
        // A sine coefficient averages to zero over anything symmetric about λ = 0, and the mean
        // of a field over the whole sphere is its C̄₀₀.
        let (half, breaks) = polygon(&[(-1.0, -0.5), (1.0, -0.5), (1.0, 0.5), (-1.0, 0.5)]);
        let symmetric = basin_kernel(6, &half as Intervals, &breaks, 32);
        assert!((0..=6).all(|l| (1..=l).all(|m| symmetric.s(l, m).abs() < 1e-16)));
        let (all, breaks) = cap(0.0, 0.0, PI);
        let everywhere = basin_kernel(6, &all as Intervals, &breaks, 48);
        assert!((everywhere.c(0, 0) - 1.0).abs() < 1e-15 && (1..=6).all(|l| everywhere.c(l, 0).abs() < 1e-15));
    }
}
