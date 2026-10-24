//! M3: a Kalman smoother for the whole COST-G record, and Q2: are COST-G's stated errors right?
//!
//! Every coefficient of degree 1–60 (anomaly against the 2004–2009 mean, after the TN-14/TN-13
//! replacements) gets a weighted least-squares fit of intercept, trend, annual, semi-annual and
//! 161-day terms. What is left is modelled as an Ornstein–Uhlenbeck signal plus the file's error
//! times a calibration factor κ. Signal size, correlation time and κ are fitted by maximum
//! likelihood per degree, shared by all its orders (`gravimetry::kalman::fit`).
//!
//! Q2 is κ itself (1 means the stated errors are right) and a cross-validation: every 10th
//! observed month is withheld from both the seasonal fit and the smoother and predicted; the
//! normalised prediction errors should have an rms of 1. The replaced C̄20 and C̄30 keep the
//! field file's errors (the technical-note readers give none); degree 1, which has none in the
//! files, gets a robust noise level from its month-to-month changes.

use std::fmt::Write;

use geopotential::Coefficients;
use gravimetry::kalman::{Process, Series, fit, smooth};
use geopotential::{julian_epoch, read_gfc};
use gravimetry::technotes::{read_tn13, read_tn14};

use crate::shared::{ols, quantile};
use crate::svg::{esc, open};
use crate::{Result, baseline_mean, calendar, corrected_month, gfc_files, span, write};

const NMAX: usize = 60;
const TAUS_MONTHS: [f64; 8] = [1.0, 2.0, 3.0, 6.0, 12.0, 24.0, 48.0, 96.0];
const SD_FACTORS: [f64; 12] = [0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.5, 4.0, 7.0, 12.0, 20.0, 35.0];
const ERROR_SCALES: [f64; 11] = [0.3, 0.45, 0.6, 0.8, 1.0, 1.25, 1.6, 2.0, 2.5, 3.2, 4.0];

/// One coefficient: degree, order, and whether it is the sine term.
type Term = (usize, usize, bool);

fn design(t: f64) -> Vec<f64> {
    let w = |p: f64| core::f64::consts::TAU * t / p;
    let alias = 161.0 / 365.25;
    vec![1.0, t - 2015.0, w(1.0).cos(), w(1.0).sin(), w(0.5).cos(), w(0.5).sin(), w(alias).cos(), w(alias).sin()]
}

/// Weighted seasonal fit on the steps that are observed and `used`; returns the fitted curve at
/// every step.
fn seasonal_fit(times: &[f64], obs: &[Option<(f64, f64)>], used: &dyn Fn(usize) -> bool) -> Result<Vec<f64>> {
    let (mut x, mut y) = (Vec::new(), Vec::new());
    for (k, o) in obs.iter().enumerate() {
        if let (Some((v, e)), true) = (o, used(k)) {
            x.push(design(times[k]).into_iter().map(|d| d / e).collect::<Vec<f64>>());
            y.push(v / e);
        }
    }
    let (b, _) = ols(&x, &y)?;
    Ok(times.iter().map(|&t| design(t).iter().zip(&b).map(|(d, b)| d * b).sum()).collect())
}

pub fn run(months_dir: &str, baseline_dir: &str, tn14: &str, tn13: &str, figures: &str) -> Result<String> {
    let tn14 = read_tn14(&std::fs::read_to_string(tn14)?)?;
    let tn13 = read_tn13(&std::fs::read_to_string(tn13)?)?;
    let (mean, _) = baseline_mean(baseline_dir, &tn14, &tn13)?;
    let (gm, radius) = (mean.gm, mean.radius);

    // Observed months: anomaly coefficients and their stated errors.
    let mut obs_t = Vec::new();
    let mut fields = Vec::new();
    let mut sigmas = Vec::new();
    for p in gfc_files(months_dir)? {
        let sigma = read_gfc(&std::fs::read_to_string(&p)?)?.sigma.ok_or("a month without errors")?;
        fields.push(gravimetry::anomaly(&corrected_month(&p, &tn14, &tn13)?, &mean).truncated(NMAX));
        sigmas.push(sigma.rescaled(gm, radius).truncated(NMAX));
        obs_t.push(span(&p)?.mid());
    }
    // Time steps: every observation, plus the middle of every calendar month with none.
    let (first, last) = (calendar(obs_t[0]), calendar(*obs_t.last().ok_or("no months")?));
    let mut calendar_months = Vec::new();
    let (mut y, mut m) = first;
    while (y, m) <= last {
        calendar_months.push((y, m));
        (y, m) = if m == 12 { (y + 1, 1) } else { (y, m + 1) };
    }
    let mut steps: Vec<(f64, Option<usize>)> = obs_t.iter().enumerate().map(|(i, &t)| (t, Some(i))).collect();
    for &(y, m) in &calendar_months {
        if !obs_t.iter().any(|&t| calendar(t) == (y, m)) {
            steps.push((julian_epoch(y, m, 16, 0, 0), None));
        }
    }
    steps.sort_by(|a, b| a.0.total_cmp(&b.0));
    let times: Vec<f64> = steps.iter().map(|s| s.0).collect();
    // The step shown for each calendar month: its first observation, or the inserted step.
    let frame_step: Vec<usize> = calendar_months.iter().map(|&cm| (0..steps.len()).find(|&k| calendar(times[k]) == cm).ok_or("no step for a month")).collect::<std::result::Result<_, _>>()?;
    eprintln!("kalman: {} observed months, {} time steps, {} calendar months", obs_t.len(), steps.len(), calendar_months.len());

    let terms: Vec<Term> = (1..=NMAX).flat_map(|n| (0..=n).flat_map(move |m| [(n, m, false), (n, m, true)].into_iter().filter(move |t| !(t.2 && t.1 == 0)))).collect();
    let get = |f: &Coefficients, (n, m, s): Term| if s { f.s(n, m) } else { f.c(n, m) };
    // Degree 1 (from TN-13) has no stated error in the files: give such coefficients a robust
    // noise level from their month-to-month changes, 1.4826 · MAD(Δy) / √2.
    let series = |term: Term| -> Vec<Option<(f64, f64)>> {
        let raw: Vec<Option<(f64, f64)>> = steps.iter().map(|&(_, i)| i.map(|i| (get(&fields[i], term), get(&sigmas[i], term)))).collect();
        if raw.iter().flatten().any(|o| o.1 > 1e-14) {
            return raw;
        }
        let obs: Vec<(f64, f64)> = steps.iter().zip(&raw).filter_map(|(s, o)| o.map(|o| (s.0, o.0))).collect();
        let d: Vec<f64> = obs.windows(2).filter(|w| w[1].0 - w[0].0 < 45.0 / 365.25).map(|w| w[1].1 - w[0].1).collect();
        let med = quantile(&d, 0.5);
        let mad = quantile(&d.iter().map(|x| (x - med).abs()).collect::<Vec<_>>(), 0.5);
        let nominal = 1.4826 * mad / 2f64.sqrt();
        raw.iter().map(|o| o.map(|(v, _)| (v, nominal))).collect()
    };

    // Seasonal fits, residual series, per-degree hyperparameters.
    let all = |_: usize| true;
    let mut fitted = Vec::with_capacity(terms.len());
    let mut resid = Vec::with_capacity(terms.len());
    for &term in &terms {
        let obs = series(term);
        let curve = seasonal_fit(&times, &obs, &all)?;
        resid.push(obs.iter().zip(&curve).map(|(o, c)| o.map(|(v, e)| (v - c, e))).collect::<Vec<_>>());
        fitted.push(curve);
    }
    let mut params = vec![None; NMAX + 1];
    let mut edges = 0;
    for (n, slot) in params.iter_mut().enumerate().skip(1) {
        let idx: Vec<usize> = (0..terms.len()).filter(|&k| terms[k].0 == n).collect();
        let errs: Vec<f64> = idx.iter().flat_map(|&k| resid[k].iter().flatten().map(|o| o.1)).collect();
        let base = quantile(&errs, 0.5);
        let sds: Vec<f64> = SD_FACTORS.iter().map(|f| f * base).collect();
        let taus: Vec<f64> = TAUS_MONTHS.iter().map(|m| m / 12.0).collect();
        let refs: Vec<Series> = idx.iter().map(|&k| (times.as_slice(), resid[k].as_slice())).collect();
        let f = fit(&refs, &sds, &taus, &ERROR_SCALES);
        let at_edge = [f.process.sd == sds[0], f.process.sd == sds[sds.len() - 1], f.process.tau == taus[0], f.process.tau == taus[taus.len() - 1], f.error_scale == ERROR_SCALES[0], f.error_scale == ERROR_SCALES[ERROR_SCALES.len() - 1]];
        edges += at_edge.iter().filter(|&&e| e).count();
        *slot = Some((f.process, f.error_scale, base));
    }
    let param = |n: usize| params[n].expect("fitted");

    // Cross-validation: withhold every 10th observed month from the seasonal fit and the smoother.
    let held: Vec<usize> = (0..steps.len()).filter(|&j| steps[j].1.is_some_and(|i| i % 10 == 5)).collect();
    let keep = |j: usize| !held.contains(&j);
    let (mut z2, mut count) = (vec![0.0; NMAX + 1], vec![0usize; NMAX + 1]);
    for &term in &terms {
        let obs = series(term);
        let curve = seasonal_fit(&times, &obs, &keep)?;
        let (process, kappa, _) = param(term.0);
        let masked: Vec<Option<(f64, f64)>> = obs.iter().zip(&curve).enumerate().map(|(j, (o, c))| if keep(j) { o.map(|(v, e)| (v - c, e * kappa)) } else { None }).collect();
        let s = smooth(&times, &masked, process);
        for &j in &held {
            let (v, e) = obs[j].expect("held months are observed");
            let z = (v - curve[j] - s.mean[j]) / (s.sd[j].powi(2) + (e * kappa).powi(2)).sqrt();
            z2[term.0] += z * z;
            count[term.0] += 1;
        }
    }
    let rms_z = |lo: usize, hi: usize| ((lo..=hi).map(|n| z2[n]).sum::<f64>() / (lo..=hi).map(|n| count[n]).sum::<usize>() as f64).sqrt();

    // Report.
    let mut rep = String::new();
    let observed = steps.iter().filter(|s| s.1.is_some()).count();
    let _ = writeln!(rep, "{observed} observed months, {} filled, {} calendar months {}-{:02} to {}-{:02}", calendar_months.len() - frame_step.iter().filter(|&&j| steps[j].1.is_some()).count(), calendar_months.len(), first.0, first.1, last.0, last.1);
    let kappas: Vec<f64> = (2..=NMAX).map(|n| param(n).1).collect();
    let _ = writeln!(rep, "error scale κ (1 = COST-G's errors are right), degrees 2–60: median {:.2}, range {:.2}–{:.2}", quantile(&kappas, 0.5), kappas.iter().fold(f64::MAX, |a, &b| a.min(b)), kappas.iter().fold(0.0, |a: f64, &b| a.max(b)));
    for (lo, hi) in [(2, 10), (11, 30), (31, 60)] {
        let k: Vec<f64> = (lo..=hi).map(|n| param(n).1).collect();
        let tau: Vec<f64> = (lo..=hi).map(|n| param(n).0.tau * 12.0).collect();
        let snr: Vec<f64> = (lo..=hi).map(|n| param(n).0.sd / param(n).2 / param(n).1).collect();
        let _ = writeln!(rep, "  degrees {lo}–{hi}: median κ {:.2}, correlation time {:.0} months, signal/noise {:.2}; cross-validation rms z {:.2}", quantile(&k, 0.5), quantile(&tau, 0.5), quantile(&snr, 0.5), rms_z(lo, hi));
    }
    let _ = writeln!(rep, "cross-validation (every 10th month withheld, {} months): rms z {:.2} over degrees 1–60 (1 = calibrated)", held.len(), rms_z(1, NMAX));
    let _ = writeln!(rep, "grid-edge hits in the hyperparameter search: {edges} of {}", 6 * NMAX);
    let mut table = String::from("degree,signal_sd,tau_months,error_scale,median_error,cv_rms_z\n");
    for n in 1..=NMAX {
        let (p, k, base) = param(n);
        let _ = writeln!(table, "{n},{:.4e},{:.0},{k:.2},{base:.4e},{:.3}", p.sd, p.tau * 12.0, (z2[n] / count[n].max(1) as f64).sqrt());
    }
    write(figures, "kalman_degrees.csv", &table)?;
    write(figures, "kalman_report.txt", &rep)?;
    write(figures, "kalman.svg", figure(&params, &z2, &count, &rep))?;
    Ok(rep)
}

type Params = Vec<Option<(Process, f64, f64)>>;

fn figure(params: &Params, z2: &[f64], count: &[usize], rep: &str) -> String {
    let (w, h) = (960.0, 560.0);
    let mut s = open(w, h, "Kalman smoother calibration");
    let _ = write!(s, r#"<text x="60" y="32" class="title">Are COST-G's stated errors right?</text>
<text x="60" y="54" class="subtitle">Per degree, from the Kalman smoother's maximum likelihood over 24 years of monthly fields (κ from a grid of candidate values), and from predicting withheld months.</text>
"#);
    let panel = |s: &mut String, title: &str, x0: f64, y0: f64, vals: &[(usize, f64)], lo: f64, hi: f64, log: bool, ticks: &[f64], cls: &str| {
        let (pw, ph) = (380.0, 300.0);
        let tr = |v: f64| if log { v.ln() } else { v };
        let x = |n: f64| x0 + pw * (n - 1.0) / 59.0;
        let y = |v: f64| y0 + ph * (tr(hi) - tr(v.clamp(lo, hi))) / (tr(hi) - tr(lo));
        let _ = writeln!(s, r#"<text x="{x0}" y="{:.1}" class="panel">{}</text>"#, y0 - 14.0, esc(title));
        for &t in ticks {
            let (yt, class) = (y(t), if t == 1.0 { "axis" } else { "grid" });
            let _ = writeln!(s, r#"<line x1="{x0}" x2="{:.1}" y1="{yt:.1}" y2="{yt:.1}" class="{class}"/><text x="{:.1}" y="{:.1}" class="tick" text-anchor="end">{t}</text>"#, x0 + pw, x0 - 6.0, yt + 4.0);
        }
        for n in [1, 10, 20, 30, 40, 50, 60] {
            let _ = writeln!(s, r#"<text x="{:.1}" y="{:.1}" class="tick" text-anchor="middle">{n}</text>"#, x(n as f64), y0 + ph + 16.0);
        }
        let _ = writeln!(s, r#"<text x="{:.1}" y="{:.1}" class="tick" text-anchor="middle">degree</text>"#, x0 + pw / 2.0, y0 + ph + 32.0);
        let d: String = vals.iter().enumerate().map(|(i, &(n, v))| format!("{}{:.1},{:.1}", if i == 0 { 'M' } else { 'L' }, x(n as f64), y(v))).collect();
        let _ = writeln!(s, r#"<path d="{d}" class="line {cls}"/>"#);
    };
    let kappa: Vec<(usize, f64)> = (1..params.len()).filter_map(|n| params[n].map(|p| (n, p.1))).collect();
    let cv: Vec<(usize, f64)> = (1..params.len()).map(|n| (n, (z2[n] / count[n].max(1) as f64).sqrt())).collect();
    panel(&mut s, "a. Error scale κ: the data's noise over COST-G's stated error", 70.0, 110.0, &kappa, 0.25, 4.0, true, &[0.25, 0.5, 1.0, 2.0, 4.0], "s1");
    panel(&mut s, "b. Withheld months: rms of normalised prediction errors", 540.0, 110.0, &cv, 0.5, 2.0, true, &[0.5, 1.0, 2.0], "s2");
    for (i, line) in rep.lines().filter(|l| l.starts_with("error scale") || l.starts_with("cross-validation")).enumerate() {
        let _ = writeln!(s, r#"<text x="60" y="{:.1}" class="subtitle">{}</text>"#, 490.0 + 20.0 * i as f64, esc(line));
    }
    s.push_str("</svg>\n");
    s
}
