//! Q0: does the noise of a COST-G month depend on solar activity, geomagnetic storms, the
//! accelerometer transplant, GRACE-FO's wide-pointing mode or how the centres were weighted?
//!
//! Noise metrics per month (fixed before looking):
//! 1. Open-ocean residual RMS, mm water: degree ≤ 60, 300 km Gaussian, no destriping; 2° cells
//!    within 60° of the equator and at least 3 cells from any land; each cell's trend, annual,
//!    semi-annual and 161-day terms removed over the whole record.
//! 2. Geoid amplitude of degrees 61–90 against the 2004–2009 mean, mm: almost all noise.
//! 3. |C₂₀(COST-G) − C₂₀(laser ranging)| ×10¹⁰, after removing each mission's median difference.
//!
//! Model: ln(metric) on mission-phase dummies (GRACE from 2011, GRACE transplant from September
//! 2016, GRACE-FO, wide-pointing from July 2023) and the standardised span means of F10.7 and
//! ln Ap, plus the effective number of centres 1/Σw² from the file header. Newey–West errors with
//! 12 lags; Benjamini–Hochberg over the 15 tests of interest. Output: a noise weight per month,
//! (median RMS / RMS)², from metric 1.

use std::fmt::Write;
use std::path::{Path, PathBuf};

use geopotential::Functional;
use gravity_field::{julian_epoch, read_gfc, read_tn13, read_tn14};

use crate::lab::{bh, land_fraction, normal_p, ols, ols_hac, quantile, ym};
use crate::svg::{esc, open};
use crate::{Result, baseline_mean, corrected_month, gfc_files, grid_axes, span};

const NAMES: [&str; 8] = ["intercept", "GRACE from 2011", "GRACE transplant", "GRACE-FO", "wide-pointing", "F10.7", "Ap", "centre evenness"];
/// Columns whose effects are tested (the others are nuisance).
const TESTED: [usize; 5] = [2, 4, 5, 6, 7];
const METRICS: [&str; 3] = ["ocean residual RMS", "degrees 61–90", "C₂₀ offset"];

pub struct Outcome {
    pub report: String,
    pub svg: String,
    pub csv: String,
    /// Metric, covariate, % effect with its 95% interval, p and q.
    pub effects: String,
}

struct Month {
    t: f64,
    mission: &'static str,
    metrics: [f64; 3],
    f107: f64,
    ap: f64,
    neff: f64,
}

/// Effective number of contributing centres, 1/Σw², from `CMMNT 0.20*…` header lines.
fn centre_evenness(path: &Path) -> Result<f64> {
    let text = std::fs::read_to_string(path)?;
    let w: Vec<f64> = text.lines().take_while(|l| !l.starts_with("end_of_head"))
        .filter_map(|l| l.strip_prefix("CMMNT ")?.split_once('*')?.0.trim().parse().ok()).collect();
    if w.is_empty() {
        return Err(format!("{}: no centre weights", path.display()).into());
    }
    Ok(1.0 / w.iter().map(|x| x * x).sum::<f64>())
}

/// Daily (epoch, F10.7 observed, Ap) from the GFZ `Kp_ap_Ap_SN_F107_since_1932.txt` file.
fn read_kp(path: &Path) -> Result<Vec<(f64, f64, f64)>> {
    let text = std::fs::read_to_string(path)?;
    Ok(text.lines().filter(|l| !l.starts_with('#')).filter_map(|l| {
        let f: Vec<&str> = l.split_whitespace().collect();
        let (y, m, d) = (f.first()?.parse().ok()?, f.get(1)?.parse().ok()?, f.get(2)?.parse().ok()?);
        let (ap, f107): (f64, f64) = (f.get(23)?.parse().ok()?, f.get(25)?.parse().ok()?);
        (ap >= 0.0 && f107 > 0.0).then(|| (julian_epoch(y, m, d, 12, 0), f107, ap))
    }).collect())
}

/// 2° cells within 60° of the equator with no land in the 7×7 block around them.
fn open_ocean(lats: &[f64], lons: &[f64]) -> Vec<usize> {
    let land = land_fraction(lats, lons);
    let (nlat, nlon) = (lats.len() as i64, lons.len() as i64);
    (0..lats.len() * lons.len()).filter(|&c| {
        let (i, j) = ((c / lons.len()) as i64, (c % lons.len()) as i64);
        lats[i as usize].to_degrees().abs() <= 60.0
            && (-3..=3).all(|di| (-3..=3).all(|dj| {
                let ii = (i + di).clamp(0, nlat - 1);
                land[(ii * nlon + (j + dj).rem_euclid(nlon)) as usize] == 0.0
            }))
    }).collect()
}

pub fn run(months_dir: &str, baseline_dir: &str, tn14: &str, tn13: &str, kp: &Path) -> Result<Outcome> {
    let tn14 = read_tn14(&std::fs::read_to_string(tn14)?)?;
    let tn13 = read_tn13(&std::fs::read_to_string(tn13)?)?;
    let (mean, _) = baseline_mean(baseline_dir, &tn14, &tn13)?;
    let daily = read_kp(kp)?;
    let (lats, lons) = grid_axes(2.0);
    let ocean = open_ocean(&lats, &lons);
    let weight: Vec<f64> = ocean.iter().map(|&c| lats[c / lons.len()].cos()).collect();
    let paths: Vec<PathBuf> = gfc_files(months_dir)?;
    let love = gravimetry::love_numbers(60).ok_or("degree above the Love-number table")?;
    let gauss = gravimetry::gaussian_weights(300e3, mean.radius, 60);

    let mut months = Vec::new();
    let mut ocean_series: Vec<Vec<f64>> = Vec::new(); // [month][ocean cell], mm
    let mut c20_raw = Vec::new();
    for p in &paths {
        let sp = span(p)?;
        let raw = read_gfc(&std::fs::read_to_string(p)?)?.model.base;
        let corrected = corrected_month(p, &tn14, &tn13)?;
        let anomaly = gravimetry::anomaly(&corrected, &mean);
        let dv = anomaly.degree_variances();
        let hideg = mean.radius * dv.get(61..=90).ok_or("month below degree 90")?.iter().sum::<f64>().sqrt() * 1e3;
        let smooth = gravimetry::smooth(&anomaly.truncated(60), &gauss);
        let grid = smooth.grid(&Functional::WaterHeight { love: &love }, mean.radius, &lats, &lons);
        ocean_series.push(ocean.iter().map(|&c| grid[c] * 1e3).collect());
        c20_raw.push((raw.c(2, 0) - corrected.c(2, 0)) * 1e10);
        let days: Vec<&(f64, f64, f64)> = daily.iter().filter(|d| (sp.start..sp.end).contains(&d.0)).collect();
        if days.is_empty() {
            return Err(format!("{}: no F10.7/Ap days in its span", p.display()).into());
        }
        let avg = |k: fn(&(f64, f64, f64)) -> f64| days.iter().map(|d| k(d)).sum::<f64>() / days.len() as f64;
        let mission = if p.to_string_lossy().contains("_GRFO_") { "GRACE-FO" } else { "GRACE" };
        months.push(Month { t: sp.mid(), mission, metrics: [0.0, hideg, 0.0], f107: avg(|d| d.1), ap: avg(|d| d.2), neff: centre_evenness(p)? });
    }
    // File names sort by start date; HAC needs time order.
    if months.windows(2).any(|w| w[1].t <= w[0].t) {
        return Err("months are not in time order".into());
    }
    eprintln!("q0: {} months, {} open-ocean cells", months.len(), ocean.len());

    // Metric 1: per-cell seasonal + trend residuals, area-weighted RMS per month.
    let season = |t: f64| {
        let w = |p: f64| core::f64::consts::TAU * t / p;
        let alias = 161.0 / 365.25;
        vec![1.0, t - 2015.0, w(1.0).cos(), w(1.0).sin(), w(0.5).cos(), w(0.5).sin(), w(alias).cos(), w(alias).sin()]
    };
    let xs: Vec<Vec<f64>> = months.iter().map(|m| season(m.t)).collect();
    let mut sq = vec![0.0; months.len()];
    for (k, w) in weight.iter().enumerate() {
        let y: Vec<f64> = ocean_series.iter().map(|g| g[k]).collect();
        let (_, resid) = ols(&xs, &y)?;
        resid.iter().zip(sq.iter_mut()).for_each(|(e, s)| *s += w * e * e);
    }
    let wsum: f64 = weight.iter().sum();
    sq.iter().zip(months.iter_mut()).for_each(|(s, m)| m.metrics[0] = (s / wsum).sqrt());
    // Metric 3: C20 offset from laser ranging, each mission's median removed.
    for mission in ["GRACE", "GRACE-FO"] {
        let idx: Vec<usize> = (0..months.len()).filter(|&i| months[i].mission == mission).collect();
        let med = quantile(&idx.iter().map(|&i| c20_raw[i]).collect::<Vec<_>>(), 0.5);
        idx.iter().for_each(|&i| months[i].metrics[2] = (c20_raw[i] - med).abs().max(1e-3));
    }

    // Design matrix.
    let z = |v: Vec<f64>| -> Vec<f64> {
        let m = v.iter().sum::<f64>() / v.len() as f64;
        let s = (v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (v.len() - 1) as f64).sqrt();
        v.iter().map(|x| (x - m) / s).collect()
    };
    let zf = z(months.iter().map(|m| m.f107).collect());
    let za = z(months.iter().map(|m| m.ap.ln()).collect());
    let zn = z(months.iter().map(|m| m.neff).collect());
    let (late, transplant, wpm) = (julian_epoch(2011, 1, 1, 0, 0), julian_epoch(2016, 9, 1, 0, 0), julian_epoch(2023, 7, 1, 0, 0));
    let x: Vec<Vec<f64>> = months.iter().enumerate().map(|(i, m)| {
        let grace = m.mission == "GRACE";
        let flag = |b: bool| if b { 1.0 } else { 0.0 };
        vec![1.0, flag(grace && m.t >= late && m.t < transplant), flag(grace && m.t >= transplant), flag(!grace), flag(!grace && m.t >= wpm), zf[i], za[i], zn[i]]
    }).collect();

    let mut fits = Vec::new();
    for k in 0..3 {
        let y: Vec<f64> = months.iter().map(|m| m.metrics[k].ln()).collect();
        fits.push(ols_hac(&x, &y, 12)?);
    }
    let p: Vec<f64> = fits.iter().flat_map(|(b, se, _)| TESTED.iter().map(|&c| normal_p(b[c] / se[c]))).collect();
    let q = bh(&p);

    let mut rep = String::new();
    let _ = writeln!(rep, "{} months ({} to {}), {} open-ocean cells; F10.7 and Ap from GFZ daily values over each month's span", months.len(), label(months[0].t), label(months[months.len() - 1].t), ocean.len());
    let _ = writeln!(rep, "effects as % change in the metric (dummies: against early GRACE; F10.7, Ap, evenness: per standard deviation); Newey–West errors, 12 lags");
    for (k, (b, se, resid)) in fits.iter().enumerate() {
        let y: Vec<f64> = months.iter().map(|m| m.metrics[k].ln()).collect();
        let my = y.iter().sum::<f64>() / y.len() as f64;
        let r2 = 1.0 - resid.iter().map(|e| e * e).sum::<f64>() / y.iter().map(|v| (v - my).powi(2)).sum::<f64>();
        let _ = writeln!(rep, "\n{} (median {:.2}{}; R² = {r2:.2})", METRICS[k], quantile(&months.iter().map(|m| m.metrics[k]).collect::<Vec<_>>(), 0.5), ["mm", "mm", "×10⁻¹⁰"][k]);
        for c in 1..NAMES.len() {
            let pct = |v: f64| 100.0 * (v.exp() - 1.0);
            let tested = TESTED.iter().position(|&t| t == c).map(|j| format!(", q = {:.3}", q[k * TESTED.len() + j])).unwrap_or_default();
            let _ = writeln!(rep, "  {:<18} {:>+7.1}%  [{:+.1}, {:+.1}]  p = {:.3}{tested}", NAMES[c], pct(b[c]), pct(b[c] - 1.96 * se[c]), pct(b[c] + 1.96 * se[c]), normal_p(b[c] / se[c]));
        }
    }
    let q_ref = &q;
    let flagged: Vec<String> = (0..3).flat_map(|k| TESTED.iter().enumerate().filter(move |&(j, _)| q_ref[k * TESTED.len() + j] < 0.1).map(move |(_, &c)| format!("{} → {}", NAMES[c], METRICS[k]))).collect();
    let _ = writeln!(rep, "\nq < 0.1: {}", if flagged.is_empty() { "none".into() } else { flagged.join("; ") });

    let med = quantile(&months.iter().map(|m| m.metrics[0]).collect::<Vec<_>>(), 0.5);
    let mut csv = String::from("epoch_year,mission,ocean_rms_mm,deg61_90_mm,c20_offset_e10,f107_sfu,ap,centre_evenness,ocean_rms_model_mm,noise_weight\n");
    let fitted: Vec<f64> = x.iter().map(|r| r.iter().zip(&fits[0].0).map(|(a, b)| a * b).sum::<f64>().exp()).collect();
    for (i, m) in months.iter().enumerate() {
        let _ = writeln!(csv, "{:.4},{},{:.3},{:.3},{:.3},{:.1},{:.1},{:.2},{:.3},{:.3}", m.t, m.mission, m.metrics[0], m.metrics[1], m.metrics[2], m.f107, m.ap, m.neff, fitted[i], (med / m.metrics[0]).powi(2));
    }
    let svg = figure(&months, &fitted, &rep);
    let mut effects = String::from("metric,covariate,effect_pct,lo_pct,hi_pct,p,q\n");
    for (k, (b, se, _)) in fits.iter().enumerate() {
        for c in 1..NAMES.len() {
            let pct = |v: f64| 100.0 * (v.exp() - 1.0);
            let q = TESTED.iter().position(|&t| t == c).map(|j| format!("{:.4}", q[k * TESTED.len() + j])).unwrap_or_default();
            let _ = writeln!(effects, "{},{},{:.1},{:.1},{:.1},{:.4},{q}", METRICS[k], NAMES[c], pct(b[c]), pct(b[c] - 1.96 * se[c]), pct(b[c] + 1.96 * se[c]), normal_p(b[c] / se[c]));
        }
    }
    Ok(Outcome { report: rep, svg, csv, effects })
}

fn label(t: f64) -> String {
    let (y, m) = ym(t);
    format!("{y}-{m:02}")
}

/// Panel title, values and an optional fitted line.
type Panel<'a> = (&'a str, Vec<f64>, Option<&'a [f64]>);

fn figure(months: &[Month], fitted: &[f64], rep: &str) -> String {
    let (w, left, right) = (960.0, 80.0, 930.0);
    let (x0, x1) = (2002.0, 2027.0);
    let x = |t: f64| left + (right - left) * (t - x0) / (x1 - x0);
    let panels: [Panel; 3] = [
        ("a. Open-ocean residual RMS (mm of water, log scale); blue: covariate model", months.iter().map(|m| m.metrics[0]).collect(), Some(fitted)),
        ("b. Degrees 61–90 against the 2004–2009 mean (mm of geoid, log scale)", months.iter().map(|m| m.metrics[1]).collect(), None),
        ("c. F10.7 solar radio flux, mean over each month's span (sfu)", months.iter().map(|m| m.f107).collect(), None),
    ];
    let (top, ph, gap) = (106.0, 120.0, 46.0);
    let h = top + 3.0 * (ph + gap) + 80.0;
    let mut s = open(w, h, "Noise in COST-G months");
    let _ = write!(s, r#"<text x="{left}" y="32" class="title">What makes a GRACE month noisy?</text>
<text x="{left}" y="54" class="subtitle">COST-G monthly fields, 2002–2026. Shading: GRACE single-accelerometer and transplant months (from September 2016)</text>
<text x="{left}" y="72" class="subtitle">and GRACE-FO wide-pointing mode (from July 2023).</text>
"#);
    for (k, (title, v, fit)) in panels.iter().enumerate() {
        let y0 = top + k as f64 * (ph + gap);
        // Noise metrics on a log axis (the regression is in logs; a few bad months span 10×);
        // F10.7 on a linear one.
        let log = k < 2;
        let all: Vec<f64> = v.iter().chain(fit.iter().flat_map(|f| f.iter())).copied().collect();
        let (lo, hi) = if log {
            (all.iter().fold(f64::MAX, |a, &b| a.min(b)) * 0.9, all.iter().fold(0.0, |a: f64, &b| a.max(b)) * 1.1)
        } else {
            (0.0, (all.iter().fold(0.0, |a: f64, &b| a.max(b)) / 50.0).ceil() * 50.0)
        };
        let tr = |val: f64| if log { val.ln() } else { val };
        let y = |val: f64| y0 + ph * (tr(hi) - tr(val)) / (tr(hi) - tr(lo));
        let _ = writeln!(s, r#"<text x="{left}" y="{:.1}" class="panel">{}</text>"#, y0 - 10.0, esc(title));
        for (a, b) in [(2016.67, 2017.5), (2023.5, 2026.6)] {
            let _ = writeln!(s, r#"<rect x="{:.1}" y="{y0}" width="{:.1}" height="{ph}" fill="var(--s4)" opacity="0.12"/>"#, x(a), x(b) - x(a));
        }
        let ticks: Vec<f64> = if log {
            (-2..4).flat_map(|e| [1.0, 2.0, 5.0].map(|m| m * 10f64.powi(e))).filter(|t| (lo..=hi).contains(t)).collect()
        } else {
            (0..=(hi / 50.0) as usize).map(|i| i as f64 * 50.0).collect()
        };
        for tick in ticks {
            let _ = writeln!(s, r#"<line x1="{left}" y1="{0:.1}" x2="{right}" y2="{0:.1}" class="grid"/><text x="{1:.1}" y="{2:.1}" class="tick" text-anchor="end">{3}</text>"#,
                y(tick), left - 6.0, y(tick) + 4.0, if tick >= 1.0 { format!("{tick:.0}") } else { format!("{tick}") });
        }
        let line = |vals: &[f64]| -> String {
            (0..months.len()).map(|i| format!("{}{:.1},{:.1}", if i == 0 || months[i].t - months[i - 1].t > 45.0 / 365.25 { 'M' } else { 'L' }, x(months[i].t), y(vals[i]))).collect()
        };
        let _ = writeln!(s, r#"<path d="{}" class="line {}" stroke-width="1.4"/>"#, line(v), if k == 2 { "s4" } else { "ink" });
        if let Some(f) = fit {
            let _ = writeln!(s, r#"<path d="{}" class="line s1" stroke-width="1.6" opacity="0.85"/>"#, line(f));
        }
    }
    let yb = top + 3.0 * (ph + gap) - gap + 18.0;
    for yr in (2004..=2026).step_by(2) {
        let _ = writeln!(s, r#"<text x="{:.1}" y="{yb:.1}" class="tick" text-anchor="middle">{yr}</text>"#, x(yr as f64));
    }
    let verdict = rep.lines().last().unwrap_or("").trim_start_matches("q < 0.1: ");
    let _ = writeln!(s, r#"<text x="{left}" y="{:.1}" class="subtitle">Effects with q &lt; 0.1 (Benjamini–Hochberg over 15 tests, Newey–West errors):</text>"#, yb + 30.0);
    let _ = writeln!(s, r#"<text x="{left}" y="{:.1}" class="subtitle">{}</text>"#, yb + 48.0, esc(verdict));
    s.push_str("</svg>\n");
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn open_ocean_mask_is_far_from_land() {
        let (lats, lons) = grid_axes(2.0);
        let cells = open_ocean(&lats, &lons);
        // Central Pacific in, Mediterranean out.
        let cell = |lat: f64, lon: f64| ((lat + 90.0) / 2.0) as usize * lons.len() + ((lon + 180.0) / 2.0) as usize;
        assert!(cells.contains(&cell(-10.0, -140.0)));
        assert!(!cells.contains(&cell(35.0, 18.0)));
        assert!((1500..6000).contains(&cells.len()), "{}", cells.len());
    }
}
