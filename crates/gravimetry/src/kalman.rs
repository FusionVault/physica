//! A Kalman smoother for monthly series with gaps, such as one gravity coefficient over the
//! GRACE and GRACE-FO record.
//!
//! The series is a random signal plus noise. The signal is an Ornstein–Uhlenbeck process: it
//! wanders with a typical size `sd` and forgets its past over a correlation time `tau`. That is a
//! Gaussian process with a Matérn-½ covariance `sd² · exp(−|Δt| / tau)`. Each observation adds
//! its own known error. The filter runs forward through time, the Rauch–Tung–Striebel smoother
//! ([Rauch et al. 1965](https://doi.org/10.2514/3.3166)) runs back, and every time step gets the
//! best estimate of the signal and its uncertainty, including steps with no data.
//!
//! ```
//! use gravimetry::kalman::{Process, smooth};
//!
//! // Five months, the middle one missing, each observed with an error of 1.
//! let times = [0.0, 1.0 / 12.0, 2.0 / 12.0, 3.0 / 12.0, 4.0 / 12.0];
//! let observed = [Some((1.0, 1.0)), Some((2.0, 1.0)), None, Some((2.5, 1.0)), Some((2.0, 1.0))];
//!
//! let signal = Process { sd: 2.0, tau: 0.5 }; // varies by about 2, over about 6 months
//! let result = smooth(&times, &observed, signal);
//!
//! // The gap is filled from its neighbours, with a larger uncertainty than they have.
//! assert!(result.mean[2] > 1.5 && result.mean[2] < 2.5);
//! assert!(result.sd[2] > result.sd[1]);
//! ```

/// The signal: an Ornstein–Uhlenbeck process.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Process {
    /// Typical size (stationary standard deviation), in the units of the data.
    pub sd: f64,
    /// Correlation time, in the units of the times (years for monthly fields).
    pub tau: f64,
}

/// The smoothed signal at every time step.
#[derive(Clone, Debug, PartialEq)]
pub struct Smoothed {
    /// Best estimate of the signal.
    pub mean: Vec<f64>,
    /// Its standard deviation.
    pub sd: Vec<f64>,
    /// Log-likelihood of the observations under the model (for fitting `Process`).
    pub log_likelihood: f64,
}

/// Smooth a series. `times` must increase; `observed[k]` is `Some((value, error))` or `None` for
/// a step with no data (a gap, or a time where you only want the estimate).
///
/// # Panics
/// If `times` and `observed` differ in length.
pub fn smooth(times: &[f64], observed: &[Option<(f64, f64)>], process: Process) -> Smoothed {
    assert_eq!(times.len(), observed.len(), "one observation slot per time");
    let n = times.len();
    let var = process.sd * process.sd;
    let decay = |k: usize| (-(times[k] - times[k - 1]) / process.tau).exp();
    let (mut pred_x, mut pred_p) = (vec![0.0; n], vec![0.0; n]);
    let (mut x, mut p) = (vec![0.0; n], vec![0.0; n]);
    let mut log_likelihood = 0.0;
    for k in 0..n {
        // Predict from the previous step (or start from the stationary distribution).
        (pred_x[k], pred_p[k]) = if k == 0 { (0.0, var) } else {
            let f = decay(k);
            (f * x[k - 1], f * f * p[k - 1] + var * (1.0 - f * f))
        };
        // Update with this step's observation, if there is one.
        (x[k], p[k]) = match observed[k] {
            Some((y, error)) => {
                let s = pred_p[k] + error * error;
                let innovation = y - pred_x[k];
                log_likelihood -= 0.5 * ((2.0 * core::f64::consts::PI * s).ln() + innovation * innovation / s);
                let gain = pred_p[k] / s;
                (pred_x[k] + gain * innovation, (1.0 - gain) * pred_p[k])
            }
            None => (pred_x[k], pred_p[k]),
        };
    }
    // Run back in time, folding in what later observations say about earlier steps.
    let (mut mean, mut v) = (x.clone(), p.clone());
    for k in (0..n.saturating_sub(1)).rev() {
        let c = p[k] * decay(k + 1) / pred_p[k + 1];
        mean[k] = x[k] + c * (mean[k + 1] - pred_x[k + 1]);
        v[k] = p[k] + c * c * (v[k + 1] - pred_p[k + 1]);
    }
    Smoothed { mean, sd: v.iter().map(|v| v.max(0.0).sqrt()).collect(), log_likelihood }
}

/// One series to fit: its times and its observations (as for [`smooth`]).
pub type Series<'a> = (&'a [f64], &'a [Option<(f64, f64)>]);

/// The best-fitting signal and error scale for several series that share them, such as all the
/// orders of one spherical-harmonic degree.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Fit {
    /// The signal process.
    pub process: Process,
    /// Factor on the stated errors: above 1 means the data are noisier than their errors say.
    pub error_scale: f64,
    /// The summed log-likelihood at the best fit.
    pub log_likelihood: f64,
}

/// Find the [`Process`] and error scale with the highest summed likelihood over `series`, each
/// given as `(times, observed)`, by trying every combination of the candidate sizes, correlation
/// times and error scales.
///
/// # Panics
/// If any candidate list is empty.
pub fn fit(series: &[Series], sds: &[f64], taus: &[f64], error_scales: &[f64]) -> Fit {
    assert!(!sds.is_empty() && !taus.is_empty() && !error_scales.is_empty(), "candidate lists must not be empty");
    let mut best: Option<Fit> = None;
    for &error_scale in error_scales {
        let scaled: Vec<Vec<Option<(f64, f64)>>> = series.iter()
            .map(|(_, obs)| obs.iter().map(|o| o.map(|(y, e)| (y, e * error_scale))).collect()).collect();
        for &sd in sds {
            for &tau in taus {
                let process = Process { sd, tau };
                let log_likelihood = series.iter().zip(&scaled).map(|((t, _), obs)| smooth(t, obs, process).log_likelihood).sum();
                if best.is_none_or(|b| log_likelihood > b.log_likelihood) {
                    best = Some(Fit { process, error_scale, log_likelihood });
                }
            }
        }
    }
    best.expect("non-empty candidates")
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The same answer from the Gaussian-process formulas, by brute-force linear algebra.
    fn brute_force(times: &[f64], observed: &[Option<(f64, f64)>], p: Process) -> (Vec<f64>, Vec<f64>, f64) {
        let k = |a: f64, b: f64| p.sd * p.sd * (-(a - b).abs() / p.tau).exp();
        let idx: Vec<usize> = (0..times.len()).filter(|&i| observed[i].is_some()).collect();
        let m = idx.len();
        // A = K + Σ over the observed times; solve by Gauss–Jordan for A⁻¹.
        let mut a = vec![vec![0.0; 2 * m]; m];
        for (r, &i) in idx.iter().enumerate() {
            for (c, &j) in idx.iter().enumerate() {
                a[r][c] = k(times[i], times[j]) + if i == j { observed[i].unwrap().1.powi(2) } else { 0.0 };
            }
            a[r][m + r] = 1.0;
        }
        for c in 0..m {
            let piv = a[c][c];
            a[c].iter_mut().for_each(|v| *v /= piv);
            for r in (0..m).filter(|&r| r != c) {
                let f = a[r][c];
                let row = a[c].clone();
                a[r].iter_mut().zip(&row).for_each(|(v, w)| *v -= f * w);
            }
        }
        let inv: Vec<Vec<f64>> = a.iter().map(|r| r[m..].to_vec()).collect();
        let y: Vec<f64> = idx.iter().map(|&i| observed[i].unwrap().0).collect();
        let alpha: Vec<f64> = (0..m).map(|r| (0..m).map(|c| inv[r][c] * y[c]).sum()).collect();
        let (mut mean, mut sd) = (Vec::new(), Vec::new());
        for &t in times {
            let ks: Vec<f64> = idx.iter().map(|&j| k(t, times[j])).collect();
            mean.push(ks.iter().zip(&alpha).map(|(a, b)| a * b).sum());
            let quad: f64 = (0..m).map(|r| (0..m).map(|c| ks[r] * inv[r][c] * ks[c]).sum::<f64>()).sum();
            sd.push((p.sd * p.sd - quad).sqrt());
        }
        // log|A| from the pivots of a fresh Cholesky.
        let mut l = vec![vec![0.0; m]; m];
        let mut logdet = 0.0;
        for r in 0..m {
            for c in 0..=r {
                let aij = k(times[idx[r]], times[idx[c]]) + if r == c { observed[idx[r]].unwrap().1.powi(2) } else { 0.0 };
                let s: f64 = (0..c).map(|q| l[r][q] * l[c][q]).sum();
                if r == c { l[r][r] = (aij - s).sqrt(); logdet += 2.0 * l[r][r].ln(); } else { l[r][c] = (aij - s) / l[c][c]; }
            }
        }
        let quad: f64 = y.iter().zip(&alpha).map(|(a, b)| a * b).sum();
        (mean, sd, -0.5 * (quad + logdet + m as f64 * (2.0 * core::f64::consts::PI).ln()))
    }

    #[test]
    fn matches_the_gaussian_process_solution() {
        let times = [0.0, 0.05, 0.13, 0.2, 0.41, 0.43, 0.9, 1.0, 1.6];
        let observed = [Some((0.3, 0.5)), Some((1.1, 0.2)), None, Some((0.8, 0.4)), Some((-0.2, 0.3)), None, Some((-1.0, 0.6)), Some((-0.4, 0.1)), None];
        let p = Process { sd: 1.3, tau: 0.35 };
        let s = smooth(&times, &observed, p);
        let (mean, sd, ll) = brute_force(&times, &observed, p);
        for k in 0..times.len() {
            assert!((s.mean[k] - mean[k]).abs() < 1e-10, "mean {k}: {} vs {}", s.mean[k], mean[k]);
            assert!((s.sd[k] - sd[k]).abs() < 1e-10, "sd {k}: {} vs {}", s.sd[k], sd[k]);
        }
        assert!((s.log_likelihood - ll).abs() < 1e-10, "{} vs {ll}", s.log_likelihood);
    }

    #[test]
    fn fit_recovers_the_simulated_process() {
        // SplitMix64 normals for a simulated OU process with known parameters.
        let mut state = 7u64;
        let mut normal = || {
            let mut u = || { state = state.wrapping_add(0x9E37_79B9_7F4A_7C15); let mut z = state; z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9); z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB); ((z ^ (z >> 31)) >> 11) as f64 / (1u64 << 53) as f64 };
            (-2.0 * u().max(1e-300).ln()).sqrt() * (core::f64::consts::TAU * u()).cos()
        };
        let truth = Process { sd: 2.0, tau: 0.5 };
        let (error, scale) = (0.5, 1.5); // stated error 0.5, real noise 0.75
        let times: Vec<f64> = (0..300).map(|k| k as f64 / 12.0).collect();
        let mut series = Vec::new();
        for _ in 0..20 {
            let mut x = truth.sd * normal();
            let mut obs = Vec::new();
            for k in 0..times.len() {
                if k > 0 {
                    let f = (-1.0 / 12.0 / truth.tau).exp();
                    x = f * x + truth.sd * (1.0 - f * f).sqrt() * normal();
                }
                obs.push(Some((x + error * scale * normal(), error)));
            }
            series.push(obs);
        }
        let refs: Vec<Series> = series.iter().map(|o| (times.as_slice(), o.as_slice())).collect();
        let f = fit(&refs, &[1.0, 1.4, 2.0, 2.8, 4.0], &[0.125, 0.25, 0.5, 1.0, 2.0], &[0.7, 1.0, 1.5, 2.1]);
        assert_eq!((f.process.sd, f.process.tau, f.error_scale), (2.0, 0.5, 1.5), "{f:?}");
    }
}
