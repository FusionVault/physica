use crate::functional::Functional;
use crate::legendre::idx;

/// How the permanent tide is treated in `C̄20`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TideSystem {
    /// No permanent tide, direct or indirect (the "conventional tide-free" system).
    TideFree,
    /// Includes the permanent deformation (indirect) part, not the direct tidal potential.
    ZeroTide,
    /// Includes both the direct and the indirect permanent parts.
    MeanTide,
    /// Not stated; cannot be converted.
    Unknown,
}

impl TideSystem {
    /// `C̄20` in this system minus `C̄20` tide-free. Uses the EGM96/EGM2008 convention
    /// `ΔC̄20 = −3.11080e-8·k/√5` with `k = 0.3` (Rapp 1989), which gives the IERS 2010 value
    /// −4.1736e-9 for zero-tide minus tide-free.
    fn c20_offset(self) -> Option<f64> {
        let direct = -3.11080e-8 / 5f64.sqrt();
        let k = 0.3;
        match self {
            TideSystem::TideFree => Some(0.0),
            TideSystem::ZeroTide => Some(direct * k),
            TideSystem::MeanTide => Some(direct * (1.0 + k)),
            TideSystem::Unknown => None,
        }
    }
}

/// A spherical-harmonic gravity field: fully normalized (4π, no Condon–Shortley phase) Stokes
/// coefficients `C̄nm, S̄nm` up to a maximum degree, with the `GM` and reference radius they
/// belong to.
///
/// The same type holds a static model, one GRACE month, a difference of months, or formal
/// errors (store `σC, σS` as the coefficients and take [`degree_amplitudes`](Self::degree_amplitudes)).
#[derive(Clone, Debug, PartialEq)]
pub struct Coefficients {
    /// Geocentric gravitational constant `GM`, m³/s².
    pub gm: f64,
    /// Reference radius `R` of the expansion, m.
    pub radius: f64,
    /// Tide system of `C̄20`.
    pub tide_system: TideSystem,
    nmax: usize,
    pub(crate) c: Vec<f64>,
    pub(crate) s: Vec<f64>,
}

impl Coefficients {
    /// All-zero coefficients up to degree `nmax`, tide system [`TideSystem::Unknown`].
    pub fn zeros(gm: f64, radius: f64, nmax: usize) -> Self {
        let len = idx(nmax, nmax) + 1;
        Coefficients {
            gm,
            radius,
            tide_system: TideSystem::Unknown,
            nmax,
            c: vec![0.0; len],
            s: vec![0.0; len],
        }
    }

    /// Maximum degree.
    pub fn max_degree(&self) -> usize {
        self.nmax
    }

    /// `C̄nm`.
    ///
    /// # Panics
    /// If `m > n` or `n > max_degree()` (also for [`s`](Self::s) and [`set`](Self::set)).
    pub fn c(&self, n: usize, m: usize) -> f64 {
        self.c[self.checked(n, m)]
    }

    /// `S̄nm`.
    pub fn s(&self, n: usize, m: usize) -> f64 {
        self.s[self.checked(n, m)]
    }

    /// Set `C̄nm` and `S̄nm`.
    pub fn set(&mut self, n: usize, m: usize, c: f64, s: f64) {
        let k = self.checked(n, m);
        self.c[k] = c;
        self.s[k] = s;
    }

    fn checked(&self, n: usize, m: usize) -> usize {
        assert!(m <= n && n <= self.nmax, "coefficient ({n},{m}) outside degree {}", self.nmax);
        idx(n, m)
    }

    /// Degree variances `c_n = Σ_m (C̄nm² + S̄nm²)`, for `n = 0..=max_degree()`. The per-coefficient
    /// degree RMS is `√(c_n / (2n+1))`.
    pub fn degree_variances(&self) -> Vec<f64> {
        (0..=self.nmax)
            .map(|n| (0..=n).map(|m| self.c[idx(n, m)].powi(2) + self.s[idx(n, m)].powi(2)).sum())
            .collect()
    }

    /// Degree amplitudes `|κ_n|·√c_n` of a functional on the reference sphere (`r = R`), e.g.
    /// geoid-height signal or error per degree. The cumulative curve is the running
    /// root-sum-square of these.
    pub fn degree_amplitudes(&self, functional: &Functional) -> Vec<f64> {
        self.degree_variances()
            .into_iter()
            .enumerate()
            .map(|(n, cn)| functional.factor(n, self.gm, self.radius, self.radius).abs() * cn.sqrt())
            .collect()
    }

    /// The field cut to degree `nmax` (or unchanged if it is already at or below it).
    pub fn truncated(&self, nmax: usize) -> Coefficients {
        let nmax = nmax.min(self.nmax);
        let len = idx(nmax, nmax) + 1;
        Coefficients { nmax, c: self.c[..len].to_vec(), s: self.s[..len].to_vec(), ..self.clone() }
    }

    /// The same field expressed with another `GM` and reference radius:
    /// `C̄'nm = C̄nm · (GM/GM') · (R/R')^n`, so that the synthesized potential is unchanged.
    pub fn rescaled(&self, gm: f64, radius: f64) -> Coefficients {
        let mut out = self.clone();
        (out.gm, out.radius) = (gm, radius);
        for n in 0..=self.nmax {
            let f = self.gm / gm * (self.radius / radius).powi(n as i32);
            for m in 0..=n {
                out.c[idx(n, m)] *= f;
                out.s[idx(n, m)] *= f;
            }
        }
        out
    }

    /// The same field in another tide system (only `C̄20` changes). `None` if either system is
    /// [`TideSystem::Unknown`] or the field has no degree 2.
    pub fn to_tide_system(&self, target: TideSystem) -> Option<Coefficients> {
        let delta = target.c20_offset()? - self.tide_system.c20_offset()?;
        if self.nmax < 2 {
            return None;
        }
        let mut out = self.clone();
        out.c[idx(2, 0)] += delta;
        out.tide_system = target;
        Some(out)
    }
}

/// What a time-variable ICGEM record contributes to a coefficient at time `t`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum TermKind {
    /// `gfct`: the value at the reference epoch.
    Value,
    /// `trnd` (or legacy `dot`): a rate per year, times `t − epoch`.
    Trend,
    /// `acos`: amplitude × `cos(2π(t − epoch)/period)`, period in years.
    Cos {
        /// Period, years.
        period: f64,
    },
    /// `asin`: amplitude × `sin(2π(t − epoch)/period)`, period in years.
    Sin {
        /// Period, years.
        period: f64,
    },
}

/// One time-variable record of coefficient `(n, m)`. All times are decimal years.
#[derive(Clone, Debug, PartialEq)]
pub struct Term {
    /// Degree.
    pub n: usize,
    /// Order.
    pub m: usize,
    /// What the record contributes.
    pub kind: TermKind,
    /// Contribution to `C̄nm`.
    pub c: f64,
    /// Contribution to `S̄nm`.
    pub s: f64,
    /// Reference epoch `t0` of the record.
    pub epoch: f64,
    /// Interval in which the record applies: `t0..t1` for icgem2.0, unbounded
    /// (`f64::NEG_INFINITY..f64::INFINITY`) for icgem1.0.
    pub valid: core::ops::Range<f64>,
}

/// A time-variable model as ICGEM defines it (e.g. GOCO06s, EIGEN-6S4):
/// `G(t) = G(t0) + trnd·(t−t0) + Σ_k [asin_k·sin(2π(t−t0)/p_k) + acos_k·cos(2π(t−t0)/p_k)]`.
#[derive(Clone, Debug, PartialEq)]
pub struct TimeVariable {
    /// Static (`gfc`) part.
    pub base: Coefficients,
    /// Time-variable records (`gfct`, `trnd`, `acos`, `asin`).
    pub terms: Vec<Term>,
}

impl TimeVariable {
    /// Evaluate the model at time `t` (decimal years).
    ///
    /// # Panics
    /// If a term's `(n, m)` lies outside the base model's degree.
    pub fn at(&self, t: f64) -> Coefficients {
        let mut out = self.base.clone();
        for term in self.terms.iter().filter(|term| term.valid.contains(&t)) {
            let dt = t - term.epoch;
            let w = match term.kind {
                TermKind::Value => 1.0,
                TermKind::Trend => dt,
                TermKind::Cos { period } => (core::f64::consts::TAU * dt / period).cos(),
                TermKind::Sin { period } => (core::f64::consts::TAU * dt / period).sin(),
            };
            let k = out.checked(term.n, term.m);
            out.c[k] += w * term.c;
            out.s[k] += w * term.s;
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn iers_c20_tide_conversion() {
        let mut zt = Coefficients::zeros(3.986004415e14, 6378136.3, 2);
        zt.tide_system = TideSystem::ZeroTide;
        zt.set(2, 0, -0.48416948e-3, 0.0);
        let tf = zt.to_tide_system(TideSystem::TideFree).unwrap();
        // IERS 2010 §6.2.2: conventional zero-tide C̄20 is −0.48416531e-3 tide-free, an offset of
        // −4.1736e-9.
        assert!((tf.c(2, 0) - -0.48416531e-3).abs() < 5e-12);
        assert!((zt.c(2, 0) - tf.c(2, 0) - -4.1736e-9).abs() < 1e-13);
        let back = tf.to_tide_system(TideSystem::MeanTide).unwrap().to_tide_system(TideSystem::ZeroTide).unwrap();
        assert!((back.c(2, 0) - zt.c(2, 0)).abs() < 1e-18);
        let mut unknown = zt.clone();
        unknown.tide_system = TideSystem::Unknown;
        assert!(unknown.to_tide_system(TideSystem::TideFree).is_none());
    }

    #[test]
    fn degree_variances_sum_orders() {
        let mut c = Coefficients::zeros(1.0, 1.0, 3);
        c.set(2, 1, 3.0, 4.0);
        c.set(3, 3, 1.0, 0.0);
        assert_eq!(c.degree_variances(), vec![0.0, 0.0, 25.0, 1.0]);
        let t = c.truncated(2);
        assert_eq!((t.max_degree(), t.c(2, 1), t.s(2, 1)), (2, 3.0, 4.0));
        assert_eq!(c.truncated(9), c);
    }

    #[test]
    #[should_panic(expected = "outside degree")]
    fn order_above_degree_is_rejected() {
        // (1, 2) would otherwise alias (2, 0) in triangular storage.
        Coefficients::zeros(1.0, 1.0, 3).set(1, 2, 1.0, 0.0);
    }

    #[test]
    fn time_variable_evaluation() {
        let mut base = Coefficients::zeros(1.0, 1.0, 2);
        base.set(2, 0, -4.8e-4, 0.0);
        let interval = 2010.0..2020.0;
        let term = |kind, c, s| Term { n: 2, m: 1, kind, c, s, epoch: 2010.0, valid: interval.clone() };
        let model = TimeVariable {
            base,
            terms: vec![
                term(TermKind::Value, 1e-9, 2e-9),
                term(TermKind::Trend, 1e-11, -1e-11),
                term(TermKind::Cos { period: 1.0 }, 3e-11, 0.0),
                term(TermKind::Sin { period: 0.5 }, 0.0, 4e-11),
                Term { epoch: 2020.0, valid: 2020.0..2030.0, ..term(TermKind::Value, 9.0, 9.0) },
            ],
        };
        let g = model.at(2012.1);
        let dt: f64 = 2012.1 - 2010.0;
        let tau = core::f64::consts::TAU;
        let c21 = 1e-9 + 1e-11 * dt + 3e-11 * (tau * dt).cos();
        let s21 = 2e-9 - 1e-11 * dt + 4e-11 * (tau * dt / 0.5).sin();
        assert!((g.c(2, 1) - c21).abs() < 1e-24);
        assert!((g.s(2, 1) - s21).abs() < 1e-24);
        assert_eq!(g.c(2, 0), -4.8e-4);
        // Outside every interval only the static part remains.
        assert_eq!(model.at(2005.0).c(2, 1), 0.0);
        assert_eq!(model.at(2025.0).c(2, 1), 9.0);
    }
}
