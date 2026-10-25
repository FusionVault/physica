//! Fully normalized associated Legendre functions `P̄nm(sin φ)` (4π normalization, no
//! Condon–Shortley phase) and their latitude derivatives, stable to degree 2700 at every latitude.
//!
//! This is the modified forward column recursion of Holmes & Featherstone (2002, J. Geodesy 76,
//! 279–299). The recursion runs on `P̄nm / cos^m φ`, scaled by 1e-280 so it cannot overflow, and
//! `cos^m φ` is applied last, the way SHTOOLS' `PlmBar` does it. A value that underflows in that
//! last step is below 1e-308, which is negligible.

/// Highest degree the scaled recursion supports before `P̄nm / cos^m φ` can overflow near the
/// poles (Holmes & Featherstone 2002).
pub const MAX_DEGREE: usize = 2700;

const SCALE: f64 = 1e-280;

/// Index of `(n, m)` in triangular storage, `n(n+1)/2 + m`.
#[inline]
pub(crate) fn idx(n: usize, m: usize) -> usize {
    n * (n + 1) / 2 + m
}

/// A reusable table of `P̄nm` (and optionally `∂P̄nm/∂φ`) for all `0 ≤ m ≤ n ≤ nmax` at one
/// latitude. The recursion coefficients are computed once in [`Legendre::new`].
///
/// ```
/// use geopotential::Legendre;
///
/// let mut leg = Legendre::new(2);
/// leg.compute(0.3);
/// let t = 0.3f64.sin();
/// assert!((leg.p(2, 0) - 5f64.sqrt() / 2.0 * (3.0 * t * t - 1.0)).abs() < 1e-15);
/// ```
#[derive(Clone, Debug)]
pub struct Legendre {
    nmax: usize,
    a: Vec<f64>,
    b: Vec<f64>,
    pub(crate) p: Vec<f64>,
    pub(crate) dp: Vec<f64>,
    pub(crate) ddp: Vec<f64>,
}

impl Legendre {
    /// Prepare a table up to degree `nmax`.
    ///
    /// # Panics
    /// If `nmax > MAX_DEGREE`: the recursion would overflow near the poles and return NaN.
    // shortcut: no extended-range (X-number) arithmetic, so degree > 2700 (e.g. XGM2019e) is
    // refused; add Fukushima (2012) X-numbers when such models are needed.
    pub fn new(nmax: usize) -> Self {
        assert!(
            nmax <= MAX_DEGREE,
            "degree {nmax} exceeds the Holmes–Featherstone limit of {MAX_DEGREE}"
        );
        let len = idx(nmax, nmax) + 1;
        let (mut a, mut b) = (vec![0.0; len], vec![0.0; len]);
        for m in 0..=nmax {
            for n in m + 1..=nmax {
                let (nf, mf) = (n as f64, m as f64);
                let k = idx(n, m);
                a[k] = ((2.0 * nf - 1.0) * (2.0 * nf + 1.0) / ((nf - mf) * (nf + mf))).sqrt();
                if n > m + 1 {
                    b[k] = ((2.0 * nf + 1.0) * (nf + mf - 1.0) * (nf - mf - 1.0)
                        / ((nf - mf) * (nf + mf) * (2.0 * nf - 3.0)))
                        .sqrt();
                }
            }
        }
        Legendre {
            nmax,
            a,
            b,
            p: vec![0.0; len],
            dp: Vec::new(),
            ddp: Vec::new(),
        }
    }

    /// Maximum degree of the table.
    pub fn max_degree(&self) -> usize {
        self.nmax
    }

    /// Fill `P̄nm(sin φ)` at geocentric latitude `lat` (radians).
    pub fn compute(&mut self, lat: f64) {
        let (t, u) = lat.sin_cos();
        let mut pmm = SCALE; // P̄mm / u^m, scaled
        let mut rescale = 1.0 / SCALE; // u^m / SCALE
        for m in 0..=self.nmax {
            if m > 0 {
                let f = if m == 1 { 3.0 } else { (2 * m + 1) as f64 / (2 * m) as f64 };
                pmm *= f.sqrt();
                rescale *= u;
            }
            self.p[idx(m, m)] = pmm * rescale;
            let (mut q2, mut q1) = (0.0, pmm);
            for n in m + 1..=self.nmax {
                let k = idx(n, m);
                let q = self.a[k] * t * q1 - self.b[k] * q2;
                self.p[k] = q * rescale;
                (q2, q1) = (q1, q);
            }
        }
    }

    /// Like [`compute`](Self::compute), and also fill `∂P̄nm/∂φ`. The derivative is a
    /// combination of `P̄n,m±1`, with no `1/cos φ` factor, so it is exact at the poles.
    pub fn compute_with_derivative(&mut self, lat: f64) {
        self.compute(lat);
        self.dp.resize(self.p.len(), 0.0);
        derive(&self.p, &mut self.dp, self.nmax);
    }

    /// Like [`compute_with_derivative`](Self::compute_with_derivative), and also fill
    /// `∂²P̄nm/∂φ²`, by applying the same pole-safe operator to the first derivatives.
    pub fn compute_with_second_derivative(&mut self, lat: f64) {
        self.compute_with_derivative(lat);
        self.ddp.resize(self.p.len(), 0.0);
        derive(&self.dp, &mut self.ddp, self.nmax);
    }

    /// `P̄nm` from the last [`compute`](Self::compute).
    ///
    /// # Panics
    /// If `m > n` or `n > max_degree()`.
    pub fn p(&self, n: usize, m: usize) -> f64 {
        assert!(m <= n && n <= self.nmax, "P̄({n},{m}) outside table");
        self.p[idx(n, m)]
    }

    /// `∂P̄nm/∂φ` from the last [`compute_with_derivative`](Self::compute_with_derivative).
    ///
    /// # Panics
    /// If `m > n`, `n > max_degree()`, or no derivative has been computed.
    pub fn dp(&self, n: usize, m: usize) -> f64 {
        assert!(m <= n && n <= self.nmax, "∂P̄({n},{m}) outside table");
        self.dp[idx(n, m)]
    }

    /// `∂²P̄nm/∂φ²` from the last
    /// [`compute_with_second_derivative`](Self::compute_with_second_derivative).
    ///
    /// # Panics
    /// If `m > n`, `n > max_degree()`, or no second derivative has been computed.
    pub fn ddp(&self, n: usize, m: usize) -> f64 {
        assert!(m <= n && n <= self.nmax, "∂²P̄({n},{m}) outside table");
        self.ddp[idx(n, m)]
    }
}

/// `dst = ∂src/∂φ` for a table `src` of functions that transform like `P̄nm` under the latitude
/// derivative: `∂P̄nm/∂φ = −½[k_m·√((n+m)(n−m+1))·P̄n,m−1 − √((n−m)(n+m+1))·P̄n,m+1]`, with
/// `k_1 = √2`, `k_m = 1` otherwise, and `∂P̄n0/∂φ = √(n(n+1)/2)·P̄n1`. Being linear, it also maps
/// `∂P̄/∂φ` to `∂²P̄/∂φ²`.
fn derive(src: &[f64], dst: &mut [f64], nmax: usize) {
    dst[0] = 0.0;
    for n in 1..=nmax {
        let nf = n as f64;
        dst[idx(n, 0)] = (nf * (nf + 1.0) / 2.0).sqrt() * src[idx(n, 1)];
        for m in 1..=n {
            let mf = m as f64;
            let km = if m == 1 { std::f64::consts::SQRT_2 } else { 1.0 };
            let lower = km * ((nf + mf) * (nf - mf + 1.0)).sqrt() * src[idx(n, m - 1)];
            let upper = if m < n { ((nf - mf) * (nf + mf + 1.0)).sqrt() * src[idx(n, m + 1)] } else { 0.0 };
            dst[idx(n, m)] = -0.5 * (lower - upper);
        }
    }
}

/// The `n` nodes (ascending) and weights of Gauss–Legendre quadrature on `[−1, 1]`, exact for
/// polynomials of degree up to `2n − 1`: Newton's method on `P_n` from Tricomi's starting values,
/// with the weights `2/((1 − x²)P_n′²)` and the symmetry about 0 imposed.
pub fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
    let (mut x, mut w) = (vec![0.0; n], vec![0.0; n]);
    for i in 0..n.div_ceil(2) {
        let mut z = (core::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        let mut dp = 0.0;
        for _ in 0..100 {
            let (mut p0, mut p1) = (1.0, z);
            for k in 2..=n {
                (p0, p1) = (p1, ((2 * k - 1) as f64 * z * p1 - (k - 1) as f64 * p0) / k as f64);
            }
            if n == 1 {
                p0 = 1.0;
            }
            dp = n as f64 * (z * p1 - p0) / (z * z - 1.0);
            let step = p1 / dp;
            z -= step;
            if step.abs() <= 1e-16 {
                break;
            }
        }
        let weight = 2.0 / ((1.0 - z * z) * dp * dp);
        (x[i], x[n - 1 - i]) = (-z, z);
        (w[i], w[n - 1 - i]) = (weight, weight);
    }
    if n % 2 == 1 {
        x[n / 2] = 0.0;
    }
    (x, w)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gauss_legendre_quadrature() {
        for n in [1, 2, 5, 17, 64] {
            let (x, w) = gauss_legendre(n);
            // Exact for x^k up to k = 2n − 1: ∫x^k = 2/(k+1) for even k, 0 for odd.
            for k in 0..2 * n {
                let sum: f64 = x.iter().zip(&w).map(|(x, w)| w * x.powi(k as i32)).sum();
                let want = if k % 2 == 0 { 2.0 / (k + 1) as f64 } else { 0.0 };
                assert!((sum - want).abs() < 2e-15, "n = {n}, x^{k}: {sum}");
            }
            assert!(x.windows(2).all(|p| p[0] < p[1]));
        }
        // The 4-point rule in closed form: ±√(3/7 ∓ (2/7)√(6/5)), weights (18 ± √30)/36.
        let (x, w) = gauss_legendre(4);
        assert!((x[3] - (3.0 / 7.0 + 2.0 / 7.0 * 1.2f64.sqrt()).sqrt()).abs() < 1e-16);
        // Two ulp, as the closed forms themselves round.
        assert!((w[3] - (18.0 - 30f64.sqrt()) / 36.0).abs() < 4e-16 && (w[2] - (18.0 + 30f64.sqrt()) / 36.0).abs() < 4e-16);
    }

    const LATS_DEG: [f64; 7] = [-89.9999, -45.0, -3.0, 0.0, 30.0, 89.9, 89.9999];

    #[test]
    fn closed_forms_low_degree() {
        let mut leg = Legendre::new(2);
        for lat in LATS_DEG.map(f64::to_radians) {
            leg.compute(lat);
            let (t, u) = lat.sin_cos();
            let expect = [
                (0, 0, 1.0),
                (1, 0, 3f64.sqrt() * t),
                (1, 1, 3f64.sqrt() * u),
                (2, 0, 5f64.sqrt() / 2.0 * (3.0 * t * t - 1.0)),
                (2, 1, 15f64.sqrt() * t * u),
                (2, 2, 15f64.sqrt() / 2.0 * u * u),
            ];
            for (n, m, v) in expect {
                assert!((leg.p(n, m) - v).abs() < 1e-14, "P̄({n},{m}) at {lat}");
            }
        }
    }

    /// Addition theorem: Σ_m P̄nm² = 2n+1 for every n, including the near-pole latitudes where
    /// a naive recursion underflows.
    #[test]
    fn sum_rule_to_degree_2700() {
        let mut leg = Legendre::new(MAX_DEGREE);
        for lat in [0.0, 45.0, 89.9, 89.9999, -60.0].map(f64::to_radians) {
            leg.compute(lat);
            let mut worst = 0.0f64;
            for n in 0..=MAX_DEGREE {
                let sum: f64 = (0..=n).map(|m| leg.p(n, m).powi(2)).sum();
                worst = worst.max((sum / (2 * n + 1) as f64 - 1.0).abs());
            }
            assert!(worst < 1e-10, "lat {lat}: worst relative error {worst:e}");
        }
    }

    #[test]
    fn derivative_matches_finite_difference() {
        let (nmax, h) = (20, 1e-5);
        let (mut leg, mut lo, mut hi) = (Legendre::new(nmax), Legendre::new(nmax), Legendre::new(nmax));
        for lat in [-80.0, -10.0, 0.0, 33.0, 89.99].map(f64::to_radians) {
            leg.compute_with_derivative(lat);
            lo.compute(lat - h);
            hi.compute(lat + h);
            for n in 0..=nmax {
                for m in 0..=n {
                    let fd = (hi.p(n, m) - lo.p(n, m)) / (2.0 * h);
                    assert!((leg.dp(n, m) - fd).abs() < 1e-5, "∂P̄({n},{m}) at {lat}");
                }
            }
        }
    }

    /// Second derivatives against the Legendre equation, P̄'' = tanφ·P̄' − (n(n+1) − m²/cos²φ)·P̄
    /// in latitude, away from the poles (where that form is singular and ours is not).
    #[test]
    fn second_derivative_satisfies_the_legendre_equation() {
        let nmax = 40;
        let mut leg = Legendre::new(nmax);
        for lat in [-1.2, -0.3, 0.0, 0.7, 1.4f64] {
            leg.compute_with_second_derivative(lat);
            let (t, u) = (lat.tan(), lat.cos());
            for n in 0..=nmax {
                for m in 0..=n {
                    let rhs = t * leg.dp(n, m) - ((n * (n + 1)) as f64 - (m * m) as f64 / (u * u)) * leg.p(n, m);
                    let scale = 1.0 + ((n * (n + 1)) as f64) * leg.p(n, m).abs();
                    assert!((leg.ddp(n, m) - rhs).abs() < 1e-9 * scale, "({n},{m}) at {lat}: {} vs {rhs}", leg.ddp(n, m));
                }
            }
        }
    }

    /// At the north pole every derivative is finite: ∂P̄n1/∂φ = −√(n(n+1)(2n+1)/2), all others 0
    /// (to rounding, since cos(π/2) is 6e-17 in f64, not 0).
    #[test]
    fn derivative_exact_at_pole() {
        let nmax = 50;
        let mut leg = Legendre::new(nmax);
        leg.compute_with_derivative(std::f64::consts::FRAC_PI_2);
        for n in 1..=nmax {
            let nf = n as f64;
            let expect = -(nf * (nf + 1.0) * (2.0 * nf + 1.0) / 2.0).sqrt();
            assert!((leg.dp(n, 1) - expect).abs() < 1e-12 * expect.abs());
            for m in (0..=n).filter(|&m| m != 1) {
                assert!(leg.dp(n, m).abs() < 1e-10, "∂P̄({n},{m}) at pole");
            }
        }
    }
}
