//! Detector noise curves: one-sided strain power spectral densities S_n(f) in Hz⁻¹.
//!
//! - [`aligo`], [`einstein_telescope`]: the analytic fits of
//!   [Sathyaprakash & Schutz (2009)](https://arxiv.org/abs/0903.0338), Table 1. They are design-era
//!   fits, good for quick estimates; use a tabulated curve ([`Tabulated`]) for precise work.
//! - [`lisa`]: the [Robson, Cornish & Liu (2019)](https://arxiv.org/abs/1803.01944) curve, with
//!   or without the galactic confusion noise.
//!
//! Below a fit's low-frequency cut-off the noise is infinite, so integrals simply ignore it.
//!
//! ```
//! use gravitational::curves::{Confusion, lisa};
//!
//! // LISA at 3 mHz, instrument noise only: √S_n ≈ 1.74×10⁻²⁰ Hz⁻¹ᐟ².
//! assert!((lisa(3e-3, Confusion::None).sqrt() / 1.74e-20 - 1.0).abs() < 0.01);
//! ```

use core::f64::consts::PI;

/// Advanced LIGO design fit: `S0 [x^−4.14 − 5x^−2 + 111(1 − x² + x⁴/2)/(1 + x²/2)]`,
/// x = f/215 Hz, S0 = 10⁻⁴⁹ Hz⁻¹, cut off below 20 Hz.
pub fn aligo(f: f64) -> f64 {
    if f < 20.0 {
        return f64::INFINITY;
    }
    let x = f / 215.0;
    1e-49 * (x.powf(-4.14) - 5.0 / (x * x) + 111.0 * (1.0 - x * x + 0.5 * x.powi(4)) / (1.0 + 0.5 * x * x))
}

/// Einstein Telescope fit (ET-B era): x = f/200 Hz, S0 = 1.5×10⁻⁵² Hz⁻¹, cut off below 10 Hz.
pub fn einstein_telescope(f: f64) -> f64 {
    if f < 10.0 {
        return f64::INFINITY;
    }
    let x = f / 200.0;
    let poly = |c: &[f64]| c.iter().rev().fold(0.0, |acc, &k| acc * x + k);
    let num = poly(&[1.0, 31.0, -65.0, 52.0, -42.0, 10.0, 12.0]);
    let den = poly(&[1.0, 14.0, -37.0, 19.0, 27.0]);
    1.5e-52 * (x.powf(-4.1) + 186.0 * x.powf(-0.69) + 233.0 * num / den)
}

/// LISA's galactic confusion noise, by mission duration.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Confusion {
    /// Instrument noise only.
    None,
    /// After six months.
    HalfYear,
    /// After one year.
    OneYear,
    /// After two years.
    TwoYears,
    /// After four years.
    FourYears,
}

/// LISA sky-averaged sensitivity (Robson, Cornish & Liu 2019, Eqs. 10–14):
/// `S_n = (10/3L²)[P_OMS + 2(1 + cos²(f/f*))P_acc/(2πf)⁴](1 + 0.6(f/f*)²) + S_c`.
pub fn lisa(f: f64, confusion: Confusion) -> f64 {
    let l = 2.5e9;
    let fstar = metrology::C / (2.0 * PI * l);
    let p_oms = (1.5e-11f64).powi(2) * (1.0 + (2e-3 / f).powi(4));
    let p_acc = (3e-15f64).powi(2) * (1.0 + (0.4e-3 / f).powi(2)) * (1.0 + (f / 8e-3).powi(4));
    let r = f / fstar;
    let instrument = 10.0 / (3.0 * l * l) * (p_oms + 2.0 * (1.0 + r.cos().powi(2)) * p_acc / (2.0 * PI * f).powi(4)) * (1.0 + 0.6 * r * r);
    let fit = match confusion {
        Confusion::None => return instrument,
        Confusion::HalfYear => (0.133, 243.0, 482.0, 917.0, 0.00258),
        Confusion::OneYear => (0.171, 292.0, 1020.0, 1680.0, 0.00215),
        Confusion::TwoYears => (0.165, 299.0, 611.0, 1340.0, 0.00173),
        Confusion::FourYears => (0.138, -221.0, 521.0, 1680.0, 0.00113),
    };
    let (alpha, beta, kappa, gamma, fk) = fit;
    let sc = 9e-45 * f.powf(-7.0 / 3.0) * (-f.powf(alpha) + beta * f * (kappa * f).sin()).exp() * (1.0 + (gamma * (fk - f)).tanh());
    instrument + sc
}

/// A tabulated amplitude spectral density (e.g. an official LIGO or Cosmic Explorer file),
/// interpolated linearly in log–log. Outside the table the noise is infinite.
#[derive(Clone, Debug, PartialEq)]
pub struct Tabulated {
    log_f: Vec<f64>,
    log_asd: Vec<f64>,
}

impl Tabulated {
    /// From `(frequency Hz, ASD Hz⁻¹ᐟ²)` pairs in increasing frequency. `None` if fewer than
    /// two points, or if any value is not positive or frequencies do not increase.
    pub fn new(points: &[(f64, f64)]) -> Option<Self> {
        let ok = points.len() >= 2 && points.iter().all(|p| p.0 > 0.0 && p.1 > 0.0) && points.windows(2).all(|w| w[1].0 > w[0].0);
        ok.then(|| Tabulated { log_f: points.iter().map(|p| p.0.ln()).collect(), log_asd: points.iter().map(|p| p.1.ln()).collect() })
    }

    /// Read whitespace-separated "frequency ASD" lines; lines starting with `#` are skipped.
    pub fn parse(text: &str) -> Option<Self> {
        let pts: Vec<(f64, f64)> = text.lines().filter(|l| !l.trim_start().starts_with('#')).filter_map(|l| {
            let mut it = l.split_whitespace().map(str::parse::<f64>);
            Some((it.next()?.ok()?, it.next()?.ok()?))
        }).collect();
        Self::new(&pts)
    }

    /// The power spectral density S_n = ASD² at `f`.
    pub fn psd(&self, f: f64) -> f64 {
        let x = f.ln();
        let n = self.log_f.len();
        if x < self.log_f[0] || x > self.log_f[n - 1] {
            return f64::INFINITY;
        }
        let i = self.log_f.partition_point(|&v| v <= x).clamp(1, n - 1);
        let t = (x - self.log_f[i - 1]) / (self.log_f[i] - self.log_f[i - 1]);
        (2.0 * (self.log_asd[i - 1] + t * (self.log_asd[i] - self.log_asd[i - 1]))).exp()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lisa_instrument_values() {
        // Robson et al. Eq. 13, recomputed: √S_n at 0.1, 1, 3, 10, 100 mHz.
        for (f, expect) in [(1e-4, 4.60e-17), (1e-3, 1.28e-19), (3e-3, 1.74e-20), (1e-2, 1.20e-20), (1e-1, 4.61e-20)] {
            let v = lisa(f, Confusion::None).sqrt();
            assert!((v / expect - 1.0).abs() < 0.01, "{f}: {v:e} vs {expect:e}");
        }
        // Confusion noise only adds, and mostly below a few mHz: S_c(1 mHz) is 6.662×10⁻³⁸ after four
        // years and 1.6635×10⁻³⁷ after one (Robson Eq. 14, Table 1).
        assert!(((lisa(1e-3, Confusion::FourYears) - lisa(1e-3, Confusion::None)) / 6.662e-38 - 1.0).abs() < 2e-3);
        assert!(((lisa(1e-3, Confusion::OneYear) - lisa(1e-3, Confusion::None)) / 1.6635e-37 - 1.0).abs() < 2e-3);
        assert!((lisa(5e-2, Confusion::FourYears) / lisa(5e-2, Confusion::None) - 1.0).abs() < 1e-3);
    }

    #[test]
    fn ground_fits_have_the_right_shape() {
        // The aLIGO fit is 1.871×10⁻²⁴ Hz⁻¹ᐟ² at 200 Hz and bottoms out at 1.800×10⁻²⁴ near 228 Hz;
        // ET is about ten times better.
        let (a, e) = (aligo(200.0).sqrt(), einstein_telescope(200.0).sqrt());
        assert!((a / 1.871e-24 - 1.0).abs() < 1e-3, "{a:e}");
        let bottom = (100..400).map(|f| aligo(f64::from(f)).sqrt()).fold(f64::INFINITY, f64::min);
        assert!((bottom / 1.800e-24 - 1.0).abs() < 2e-3, "{bottom:e}");
        assert!(aligo(100.0) > aligo(200.0) && aligo(500.0) > aligo(200.0));
        assert!((e / 1.631e-25 - 1.0).abs() < 1e-3, "{e:e}");
        assert!(aligo(10.0).is_infinite() && einstein_telescope(5.0).is_infinite());
    }

    #[test]
    fn tabulated_interpolates_in_log_log() {
        let t = Tabulated::parse("# f asd\n10 1e-22\n1000 1e-24\n").unwrap();
        assert!((t.psd(100.0).sqrt() / 1e-23 - 1.0).abs() < 1e-12);
        assert!(t.psd(5.0).is_infinite());
        assert!(Tabulated::new(&[(10.0, 1.0)]).is_none());
    }
}
