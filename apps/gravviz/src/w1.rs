//! W1: does global land water storage respond asymmetrically to El Niño and La Niña, and what does
//! that imply for 2026–27? Fits, tests, forecasts and draws the figure.
//!
//! S is COST-G land storage outside the ice sheets (km³, GRACE Tellus chain). S′ is S minus annual,
//! semi-annual and 161-day terms, fitted jointly with a linear trend. N is the Niño 3.4 anomaly,
//! N⁺ = max(N, 0) and N⁻ = min(N, 0). Two models, both fixed before any target month (December 2026,
//! March 2027, June 2027) exists:
//! - Rate model (as first announced): the central-difference rate of S′ against N⁺ and N⁻ at lags
//!   0–9 months. Its diagnostics (lag coefficients swinging ±400 between collinear lags; month-to-
//!   month noise of ~500 km³ amplified by differencing) led to:
//! - Level model (primary): S′ = a + b·(t − 2015) + γ⁺·N̄⁺ + γ⁻·N̄⁻ + e, with N̄ the mean Niño 3.4
//!   over lags 2–6 months and e an AR(1) process.
//!
//! After the last observed month, Niño 3.4 follows analogue years whose same-month anomaly exceeded
//! +1.4 °C, offset to the observed value. Parameter and noise uncertainty come from a moving-block
//! (12-month) bootstrap.

use std::collections::{BTreeSet, HashMap};
use std::fmt::Write;
use std::path::Path;

use crate::Result;
use crate::lab::{Month, Rng, block_resample, land_area_km2, land_volume_km3, months_between, ols, quantile, read_psl_index, ym};
use crate::svg::{esc, open};

type Cal = (i32, u32);
/// Coefficients, residuals and bootstrap coefficients.
type Fit = (Vec<f64>, Vec<f64>, Vec<Vec<f64>>);

const LAGS: usize = 10;
const TARGETS: [Cal; 3] = [(2026, 12), (2027, 3), (2027, 6)];
const BOOT: usize = 2000;
const DRAWS: usize = 4000;

/// Seasonal design row (no intercept or trend): annual, semi-annual and 161-day (S2 alias) terms.
fn seasonal(t: f64) -> [f64; 6] {
    let w = |p: f64| core::f64::consts::TAU * t / p;
    let alias = 161.0 / 365.25;
    [w(1.0).cos(), w(1.0).sin(), w(0.5).cos(), w(0.5).sin(), w(alias).cos(), w(alias).sin()]
}

/// The calendar month `k` months after `cm` (negative `k` goes back).
fn shift(cm: Cal, k: i32) -> Cal {
    let n = cm.0 * 12 + cm.1 as i32 - 1 + k;
    (n.div_euclid(12), n.rem_euclid(12) as u32 + 1)
}

/// Middle of a calendar month in years.
fn mid(cm: Cal) -> f64 {
    cm.0 as f64 + (cm.1 as f64 - 0.5) / 12.0
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

pub struct Outcome {
    pub report: String,
    pub svg: String,
    pub csv: String,
    pub forecast_csv: String,
    /// Month, observed Niño 3.4 and the level model's forecast quantiles, for the viewer.
    pub monthly_csv: String,
}

pub fn run(months: &[Month], lats: &[f64], lons: &[f64], nino_path: &Path) -> Result<Outcome> {
    let area = land_area_km2(lats, lons);
    let t: Vec<f64> = months.iter().map(|m| m.span.mid()).collect();
    let s: Vec<f64> = months.iter().map(|m| land_volume_km3(&m.cm, &area)).collect();
    let last = t.len() - 1;
    let consecutive = |i: usize| i > 0 && t[i] - t[i - 1] < 45.0 / 365.25;

    // S′: seasons fitted jointly with a trend, then removed (the trend stays).
    let x: Vec<Vec<f64>> = t.iter().map(|&ti| [vec![1.0, ti - 2015.0], seasonal(ti).to_vec()].concat()).collect();
    let (sb, _) = ols(&x, &s)?;
    let season = |ti: f64| dot(&seasonal(ti), &sb[2..]);
    let sp: Vec<f64> = t.iter().zip(&s).map(|(&ti, &si)| si - season(ti)).collect();

    // Niño 3.4: observed where available, else analogue year `ay`'s path offset to the last value.
    let nino: HashMap<Cal, f64> = read_psl_index(nino_path)?.into_iter().collect();
    let observed_to = *nino.keys().max().ok_or("empty Niño 3.4 file")?;
    let anchor = nino[&observed_to];
    let start = ym(t[last]);
    let analog_years: Vec<i32> = nino.iter()
        .filter(|&(&(y, m), &v)| m == observed_to.1 && v > 1.4 && y < start.0 && nino.contains_key(&(y + 1, 6)))
        .map(|(&(y, _), _)| y).collect::<BTreeSet<_>>().into_iter().collect();
    if analog_years.is_empty() {
        return Err("no Niño 3.4 analogue years".into());
    }
    let n34 = |cm: Cal, ay: i32| -> Option<f64> {
        if let Some(v) = nino.get(&cm) {
            return Some(*v);
        }
        let (k, base) = (months_between(observed_to, cm), (ay, observed_to.1));
        if k <= 0 {
            return None;
        }
        Some(nino.get(&shift(base, k))? - nino.get(&base)? + anchor)
    };

    // Rate model: central-difference rates of S′ (km³/month) on lagged N⁺ and N⁻.
    let rate_row = |cm: Cal, ay: i32| -> Option<Vec<f64>> {
        let mut row = vec![1.0];
        for positive in [true, false] {
            for k in 0..LAGS as i32 {
                let v = n34(shift(cm, -k), ay)?;
                row.push(if positive { v.max(0.0) } else { v.min(0.0) });
            }
        }
        Some(row)
    };
    let (mut rt, mut rx, mut ry) = (Vec::new(), Vec::new(), Vec::new());
    for i in 1..last {
        if let (true, true, Some(row)) = (consecutive(i), consecutive(i + 1), rate_row(ym(t[i]), 0)) {
            rt.push(t[i]);
            rx.push(row);
            ry.push((sp[i + 1] - sp[i - 1]) / (t[i + 1] - t[i - 1]) / 12.0);
        }
    }

    // Level model: S′ on trend and the lag 2–6 mean of N, split by sign.
    let level_row = |ti: f64, ay: i32| -> Option<Vec<f64>> {
        let cm = ym(ti);
        let nbar = (2..=6).map(|k| n34(shift(cm, -k), ay)).sum::<Option<f64>>()? / 5.0;
        Some(vec![1.0, ti - 2015.0, nbar.max(0.0), nbar.min(0.0)])
    };
    let lx: Vec<Vec<f64>> = t.iter().map(|&ti| level_row(ti, 0)).collect::<Option<_>>().ok_or("Niño 3.4 does not cover the record")?;

    let mut rng = Rng::new(20261024);
    let (rate_beta, rate_resid, rate_boot) = fit_boot(&rx, &ry, &mut rng)?;
    let (lev_beta, lev_resid, lev_boot) = fit_boot(&lx, &sp, &mut rng)?;

    // AR(1) for the level residuals, from consecutive months only.
    let pairs: Vec<usize> = (1..t.len()).filter(|&i| consecutive(i)).collect();
    let phi = pairs.iter().map(|&i| lev_resid[i] * lev_resid[i - 1]).sum::<f64>()
        / pairs.iter().map(|&i| lev_resid[i - 1].powi(2)).sum::<f64>();
    let innov: Vec<f64> = pairs.iter().map(|&i| lev_resid[i] - phi * lev_resid[i - 1]).collect();

    // Recovery criterion threshold: 5th percentile of 12-month residual means with full coverage.
    let means12: Vec<f64> = (0..t.len().saturating_sub(11))
        .filter(|&i| t[i + 11] - t[i] < 1.0)
        .map(|i| lev_resid[i..i + 12].iter().sum::<f64>() / 12.0).collect();
    let ratchet_threshold = quantile(&means12, 0.05);

    // Forecasts, month by month from the last field.
    let horizon = months_between(start, TARGETS[TARGETS.len() - 1]) as usize;
    let mut rate_paths = Vec::with_capacity(DRAWS);
    let mut lev_paths = Vec::with_capacity(DRAWS);
    for d in 0..DRAWS {
        let ay = analog_years[rng.below(analog_years.len())];
        let (rb, lb) = (&rate_boot[d % BOOT], &lev_boot[d % BOOT]);
        let e = block_resample(&rate_resid, 12, &mut rng);
        let mut level = sp[last];
        let mut e_lev = sp[last] - dot(&lx[last], lb);
        let (mut rp, mut lp) = (Vec::with_capacity(horizon), Vec::with_capacity(horizon));
        for h in 1..=horizon {
            let cm = shift(start, h as i32);
            level += dot(&rate_row(cm, ay).ok_or("analogue path too short")?, rb) + e[h - 1];
            rp.push(level);
            // shortcut: φ is not bootstrapped; widen the cone if φ's uncertainty ever matters.
            e_lev = phi * e_lev + innov[rng.below(innov.len())];
            lp.push(dot(&level_row(mid(cm), ay).ok_or("analogue path too short")?, lb) + e_lev);
        }
        rate_paths.push(rp);
        lev_paths.push(lp);
    }
    let at = |paths: &[Vec<f64>], h: usize| paths.iter().map(|p| p[h - 1]).collect::<Vec<f64>>();

    // Report.
    let ci = |v: &[f64]| format!("{:.0} [{:.0}, {:.0}]", quantile(v, 0.5), quantile(v, 0.05), quantile(v, 0.95));
    let mut rep = String::new();
    let w = &mut rep;
    let _ = writeln!(w, "fields: {} months, {} to {}; land area used: {:.1} million km²", t.len(), label(t[0]), label(t[last]), area.iter().sum::<f64>() / 1e6);
    let _ = writeln!(w, "season (frozen; S′ = S − these terms): annual cos/sin {:.3} {:.3}, semi-annual {:.3} {:.3}, 161-day {:.3} {:.3} km³", sb[2], sb[3], sb[4], sb[5], sb[6], sb[7]);
    let _ = writeln!(w, "trend in S: {:.1} km³/yr; annual amplitude {:.0} km³", sb[1], sb[2].hypot(sb[3]));
    let _ = writeln!(w, "Niño 3.4 observed to {}-{:02} ({:+.2} °C); analogue years (same month > +1.4 °C): {analog_years:?}", observed_to.0, observed_to.1, anchor);
    let _ = writeln!(w, "\nrate model ({} rate samples)", ry.len());
    let lags = |b: &[f64]| b.iter().map(|v| format!("{v:.0}")).collect::<Vec<_>>().join(" ");
    let _ = writeln!(w, "  intercept {:.1} km³/month; β⁺ lags 0–9: {}; β⁻ lags 0–9: {}", rate_beta[0], lags(&rate_beta[1..=LAGS]), lags(&rate_beta[LAGS + 1..]));
    let sum_p = |b: &[f64]| b[1..=LAGS].iter().sum::<f64>();
    let sum_m = |b: &[f64]| b[LAGS + 1..].iter().sum::<f64>();
    let _ = writeln!(w, "  Σβ⁺ {}; Σβ⁻ {} (km³/month per °C)", stat(&rate_beta, &rate_boot, sum_p), stat(&rate_beta, &rate_boot, sum_m));
    let _ = writeln!(w, "  asymmetry Σβ⁺ − Σβ⁻: {}", stat(&rate_beta, &rate_boot, |b| sum_p(b) - sum_m(b)));
    let _ = writeln!(w, "  fitted on GRACE, correlation with GRACE-FO rates: {}", cross_mission(&rt, &rx, &ry)?);
    let _ = writeln!(w, "\nlevel model, primary (frozen coefficients for scoring)");
    let _ = writeln!(w, "  S′ = {:.3} + {:.4}·(t − 2015) + {:.3}·N̄⁺ + {:.3}·N̄⁻ km³, N̄ = mean Niño 3.4 at lags 2–6 months", lev_beta[0], lev_beta[1], lev_beta[2], lev_beta[3]);
    let sd = (lev_resid.iter().map(|e| e * e).sum::<f64>() / (lev_resid.len() - 4) as f64).sqrt();
    let _ = writeln!(w, "  residual sd {sd:.0} km³, AR(1) φ = {phi:.3}");
    let _ = writeln!(w, "  γ⁺ (El Niño, km³ per °C): {}", stat(&lev_beta, &lev_boot, |b| b[2]));
    let _ = writeln!(w, "  γ⁻ (La Niña, km³ per °C): {}", stat(&lev_beta, &lev_boot, |b| b[3]));
    let _ = writeln!(w, "  asymmetry γ⁺ − γ⁻: {}", stat(&lev_beta, &lev_boot, |b| b[2] - b[3]));
    let _ = writeln!(w, "  fitted on GRACE, correlation with GRACE-FO S′: {}", cross_mission(&t, &lx, &sp)?);
    let _ = writeln!(w, "  ratchet criterion: mean level-model residual over 2027-07..2028-06 (observed Niño 3.4, frozen coefficients) below {ratchet_threshold:.0} km³ (5th percentile of {} full 12-month windows, 2002–2026)", means12.len());
    let _ = writeln!(w, "\nforecasts from S′ = {:.0} km³ at {} (median [5%, 95%], km³)", sp[last], label(t[last]));
    let mut fc = String::from("month,model,quantity,q05,q50,q95\n");
    for tg in TARGETS {
        let h = months_between(start, tg) as usize;
        let lev = at(&lev_paths, h);
        let full: Vec<f64> = lev.iter().map(|v| v + season(mid(tg))).collect();
        let rate = at(&rate_paths, h);
        let _ = writeln!(w, "  {}-{:02}: level model S′ {}, full anomaly incl. season {}; rate model S′ {}", tg.0, tg.1, ci(&lev), ci(&full), ci(&rate));
        for (model, quantity, v) in [("level", "s_prime", &lev), ("level", "s_full", &full), ("rate", "s_prime", &rate)] {
            let _ = writeln!(fc, "{}-{:02},{model},{quantity},{:.0},{:.0},{:.0}", tg.0, tg.1, quantile(v, 0.05), quantile(v, 0.5), quantile(v, 0.95));
        }
    }

    // Figure and data.
    let mut cone = vec![(t[last], sp[last], sp[last], sp[last])];
    cone.extend((1..=horizon).map(|h| {
        let v = at(&lev_paths, h);
        (mid(shift(start, h as i32)), quantile(&v, 0.05), quantile(&v, 0.5), quantile(&v, 0.95))
    }));
    let fit: Vec<f64> = lx.iter().map(|row| dot(row, &lev_beta)).collect();
    let caption = [
        format!("El Niño minus La Niña response, level model (km³ per °C): {}.", stat(&lev_beta, &lev_boot, |b| b[2] - b[3])),
        format!("The same for the rate model as first announced (km³/month per °C): {}.", stat(&rate_beta, &rate_boot, |b| sum_p(b) - sum_m(b))),
    ];
    let svg = figure(&t, &sp, &fit, &nino, &cone, &caption);
    let mut csv = String::from("epoch_year,mission,s_km3,s_prime_km3,level_model_km3\n");
    for i in 0..t.len() {
        let _ = writeln!(csv, "{:.4},{},{:.1},{:.1},{:.1}", t[i], months[i].mission, s[i], sp[i], fit[i]);
    }
    let mut monthly_csv = String::from("month,nino34,forecast_q05,forecast_q50,forecast_q95\n");
    let mut cm = (2002, 1);
    while cm <= TARGETS[TARGETS.len() - 1] {
        let n = nino.get(&cm).map(|v| format!("{v:.2}")).unwrap_or_default();
        let h = months_between(start, cm);
        let q = if h >= 1 && h as usize <= horizon {
            let v = at(&lev_paths, h as usize);
            format!("{:.0},{:.0},{:.0}", quantile(&v, 0.05), quantile(&v, 0.5), quantile(&v, 0.95))
        } else {
            ",,".into()
        };
        let _ = writeln!(monthly_csv, "{}-{:02},{n},{q}", cm.0, cm.1);
        cm = shift(cm, 1);
    }
    Ok(Outcome { report: rep, svg, csv, forecast_csv: fc, monthly_csv })
}

/// OLS fit, then a moving-block bootstrap (fitted values plus resampled residual blocks, refit).
fn fit_boot(x: &[Vec<f64>], y: &[f64], rng: &mut Rng) -> Result<Fit> {
    let (beta, resid) = ols(x, y)?;
    let fitted: Vec<f64> = y.iter().zip(&resid).map(|(y, e)| y - e).collect();
    let boot = (0..BOOT)
        .map(|_| {
            let e = block_resample(&resid, 12, rng);
            let yb: Vec<f64> = fitted.iter().zip(&e).map(|(f, e)| f + e).collect();
            ols(x, &yb).map(|(b, _)| b)
        })
        .collect::<Result<Vec<_>>>()?;
    Ok((beta, resid, boot))
}

/// "value, 90% CI [lo, hi], p = …" for a statistic of the coefficients; p is the two-sided
/// bootstrap probability of the other sign.
fn stat(beta: &[f64], boot: &[Vec<f64>], f: impl Fn(&[f64]) -> f64) -> String {
    let v: Vec<f64> = boot.iter().map(|b| f(b)).collect();
    let other = v.iter().filter(|&&x| x <= 0.0).count().min(v.iter().filter(|&&x| x >= 0.0).count());
    let p = if other == 0 { format!("p < {:.0e}", 1.0 / v.len() as f64) } else { format!("p = {:.3}", 2.0 * other as f64 / v.len() as f64) };
    format!("{:.1}, 90% CI [{:.1}, {:.1}], {p}", f(beta), quantile(&v, 0.05), quantile(&v, 0.95))
}

/// Fit on GRACE (before mid-2017), correlate the prediction with GRACE-FO.
fn cross_mission(t: &[f64], x: &[Vec<f64>], y: &[f64]) -> Result<String> {
    let (grace, fo): (Vec<usize>, Vec<usize>) = (0..t.len()).partition(|&i| t[i] < 2017.6);
    let (b, _) = ols(&grace.iter().map(|&i| x[i].clone()).collect::<Vec<_>>(), &grace.iter().map(|&i| y[i]).collect::<Vec<_>>())?;
    let pred: Vec<f64> = fo.iter().map(|&i| dot(&x[i], &b)).collect();
    let obs: Vec<f64> = fo.iter().map(|&i| y[i]).collect();
    Ok(format!("r = {:.2} over {} months", pearson(&pred, &obs), obs.len()))
}

fn label(t: f64) -> String {
    let (y, m) = ym(t);
    format!("{y}-{m:02}")
}

fn pearson(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len() as f64;
    let (ma, mb) = (a.iter().sum::<f64>() / n, b.iter().sum::<f64>() / n);
    let cov: f64 = a.iter().zip(b).map(|(x, y)| (x - ma) * (y - mb)).sum();
    let va: f64 = a.iter().map(|x| (x - ma).powi(2)).sum();
    let vb: f64 = b.iter().map(|y| (y - mb).powi(2)).sum();
    cov / (va * vb).sqrt()
}

fn figure(t: &[f64], sp: &[f64], fit: &[f64], nino: &HashMap<Cal, f64>, cone: &[(f64, f64, f64, f64)], caption: &[String]) -> String {
    let (w, h, left, right, top, bottom) = (960.0, 580.0, 72.0, 936.0, 110.0, 480.0);
    let (x0, x1) = (2002.0, 2027.75);
    let lo = sp.iter().chain(cone.iter().map(|c| &c.1)).fold(f64::MAX, |a, &b| a.min(b));
    let hi = sp.iter().chain(cone.iter().map(|c| &c.3)).fold(f64::MIN, |a, &b| a.max(b));
    let (lo, hi) = ((lo / 1000.0).floor() * 1000.0, (hi / 1000.0).ceil() * 1000.0);
    let x = |v: f64| left + (right - left) * (v - x0) / (x1 - x0);
    let y = |v: f64| top + (bottom - top) * (hi - v) / (hi - lo);
    let mut s = open(w, h, "Land water storage and El Niño");
    let _ = write!(s, r#"<text x="{left}" y="32" class="title">Land water storage and El Niño</text>
<text x="{left}" y="54" class="subtitle">Land outside the ice sheets, seasons removed: COST-G GRACE and GRACE-FO, km³ relative to 2004–2009. Shading: Niño 3.4 above +1 °C.</text>
<text x="{left}" y="72" class="subtitle">Black: observed. Blue: trend plus the Niño 3.4 level model (in-sample fit), then its forecast median and 90% range to June 2027,</text>
<text x="{left}" y="90" class="subtitle">committed as a git tag before those months were observed.</text>
"#);
    let mut shaded: Vec<(&Cal, &f64)> = nino.iter().collect();
    shaded.sort_by_key(|(cm, _)| **cm); // stable output across runs
    for (&(yr, mo), &v) in shaded {
        let tm = yr as f64 + (mo as f64 - 1.0) / 12.0;
        if v > 1.0 && (x0..x1).contains(&tm) {
            let _ = writeln!(s, r#"<rect x="{:.1}" y="{top}" width="{:.1}" height="{:.1}" fill="var(--s2)" opacity="0.15"/>"#, x(tm), x(tm + 1.0 / 12.0) - x(tm) + 0.3, bottom - top);
        }
    }
    let step = if hi - lo > 8000.0 { 2000.0 } else { 1000.0 };
    let mut v = lo;
    while v <= hi + 1.0 {
        let _ = writeln!(s, r#"<line x1="{left}" y1="{:.1}" x2="{right}" y2="{:.1}" class="{}"/><text x="{:.1}" y="{:.1}" class="tick" text-anchor="end">{:.0}</text>"#,
            y(v), y(v), if v == 0.0 { "axis" } else { "grid" }, left - 8.0, y(v) + 4.0, v);
        v += step;
    }
    for yr in (2004..=2026).step_by(2) {
        let _ = writeln!(s, r#"<text x="{:.1}" y="{:.1}" class="tick" text-anchor="middle">{yr}</text>"#, x(yr as f64), bottom + 18.0);
    }
    // Lines break at gaps longer than 45 days.
    let path = |vals: &[f64]| -> String {
        (0..t.len()).map(|i| format!("{}{:.1},{:.1}", if i == 0 || t[i] - t[i - 1] > 45.0 / 365.25 { 'M' } else { 'L' }, x(t[i]), y(vals[i]))).collect()
    };
    let _ = writeln!(s, r#"<path d="{}" class="line ink" stroke-width="1.5"><title>Observed S′</title></path>"#, path(sp));
    let _ = writeln!(s, r#"<path d="{}" class="line s1" stroke-width="1.5" opacity="0.8"><title>Level model fit</title></path>"#, path(fit));
    let band: String = cone.iter().map(|c| format!("{:.1},{:.1} ", x(c.0), y(c.3))).chain(cone.iter().rev().map(|c| format!("{:.1},{:.1} ", x(c.0), y(c.1)))).collect();
    let _ = writeln!(s, r#"<polygon points="{band}" fill="var(--s1)" opacity="0.25"><title>90% forecast range</title></polygon>"#);
    let med: String = cone.iter().enumerate().map(|(i, c)| format!("{}{:.1},{:.1}", if i == 0 { 'M' } else { 'L' }, x(c.0), y(c.2))).collect();
    let _ = writeln!(s, r#"<path d="{med}" class="line s1"><title>Forecast median</title></path>"#);
    for (i, line) in caption.iter().enumerate() {
        let _ = writeln!(s, r#"<text x="{left}" y="{:.1}" class="subtitle">{}</text>"#, bottom + 46.0 + 20.0 * i as f64, esc(line));
    }
    s.push_str("</svg>\n");
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn calendar_shifts() {
        assert_eq!(shift((2026, 7), 5), (2026, 12));
        assert_eq!(shift((2026, 7), 11), (2027, 6));
        assert_eq!(shift((2026, 2), -6), (2025, 8));
        assert_eq!(shift((2026, 1), -1), (2025, 12));
    }
}
