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
        self.dp[0] = 0.0;
        for n in 1..=self.nmax {
            let nf = n as f64;
            self.dp[idx(n, 0)] = (nf * (nf + 1.0) / 2.0).sqrt() * self.p[idx(n, 1)];
            for m in 1..=n {
                let mf = m as f64;
                let km = if m == 1 { std::f64::consts::SQRT_2 } else { 1.0 };
                let lower = km * ((nf + mf) * (nf - mf + 1.0)).sqrt() * self.p[idx(n, m - 1)];
                let upper = if m < n {
                    ((nf - mf) * (nf + mf + 1.0)).sqrt() * self.p[idx(n, m + 1)]
                } else {
                    0.0
                };
                self.dp[idx(n, m)] = -0.5 * (lower - upper);
            }
        }
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
}

#[cfg(test)]
mod tests {
    use super::*;

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
