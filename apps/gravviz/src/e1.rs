//! E1 quick look: is the 29 July 2025 Kamchatka Mw 8.8 earthquake visible in COST-G months?
//!
//! Gravity anomaly (µGal) of every GRACE-FO month minus the 2004–2009 mean: degree ≤ 60, 300 km
//! Gaussian, no destriping (it removes north–south structure, which a subduction dipole has), 1°
//! grid. Every cell gets one least-squares fit: intercept, trend, annual, semi-annual and 161-day
//! terms, a step at the event (July 2025 weighted by the share of its span after the event) and a
//! post-event ramp. Statistic: the largest |step| within 500 km of the epicentre. Null: the same
//! statistic around 100 random ocean points at 40–65° latitude (GRACE noise depends on latitude),
//! more than 2,500 km from the epicentre.

use std::fmt::Write;

use geopotential::Functional;
use gravity_field::{julian_epoch, read_tn13, read_tn14};

use crate::lab::{Rng, land_rings, on_land, ols, quantile};
use crate::svg::{esc, open};
use crate::{Result, Span, baseline_mean, corrected_month, gfc_files, grid_axes, span};

const EPICENTRE: (f64, f64) = (52.5, 160.3); // USGS, rounded (lat, lon degrees)
const RADIUS_KM: f64 = 500.0;
const NULL_POINTS: usize = 100;

pub struct Outcome {
    pub report: String,
    pub svg: String,
    pub csv: String,
}

fn distance_km(a: (f64, f64), b: (f64, f64)) -> f64 {
    let (la1, lo1, la2, lo2) = (a.0.to_radians(), a.1.to_radians(), b.0.to_radians(), b.1.to_radians());
    let c = la1.sin() * la2.sin() + la1.cos() * la2.cos() * (lo1 - lo2).cos();
    6371.0 * c.clamp(-1.0, 1.0).acos()
}

/// Design row: intercept, trend, six seasonal terms, step, post-event ramp (years).
fn row(sp: Span, te: f64) -> Vec<f64> {
    let t = sp.mid();
    let w = |p: f64| core::f64::consts::TAU * t / p;
    let alias = 161.0 / 365.25;
    let step = ((sp.end - te) / (sp.end - sp.start)).clamp(0.0, 1.0);
    vec![1.0, t - 2022.0, w(1.0).cos(), w(1.0).sin(), w(0.5).cos(), w(0.5).sin(), w(alias).cos(), w(alias).sin(), step, (t - te).max(0.0) * step]
}

pub fn run(months_dir: &str, baseline_dir: &str, tn14: &str, tn13: &str) -> Result<Outcome> {
    let tn14 = read_tn14(&std::fs::read_to_string(tn14)?)?;
    let tn13 = read_tn13(&std::fs::read_to_string(tn13)?)?;
    let (mean, _) = baseline_mean(baseline_dir, &tn14, &tn13)?;
    let (lats, lons) = grid_axes(1.0);
    let te = julian_epoch(2025, 7, 29, 23, 25);
    let paths: Vec<_> = gfc_files(months_dir)?.into_iter().filter(|p| p.to_string_lossy().contains("_GRFO_")).collect();
    let mut spans = Vec::new();
    let mut grids = Vec::new();
    for p in &paths {
        let anomaly = gravimetry::anomaly(&corrected_month(p, &tn14, &tn13)?, &mean).truncated(60);
        let anomaly = gravimetry::smooth(&anomaly, &gravimetry::gaussian_weights(300e3, mean.radius, 60));
        grids.push(anomaly.grid(&Functional::GravityAnomaly, mean.radius, &lats, &lons).into_iter().map(|g| g * 1e8).collect::<Vec<f64>>());
        spans.push(span(p)?);
    }
    eprintln!("e1: {} GRACE-FO months", grids.len());

    // The step and ramp estimates are linear in the data: w·y, with w from one OLS per unit vector.
    let x: Vec<Vec<f64>> = spans.iter().map(|&s| row(s, te)).collect();
    let mut weights = vec![[0.0; 2]; x.len()];
    for (j, w) in weights.iter_mut().enumerate() {
        let unit: Vec<f64> = (0..x.len()).map(|i| if i == j { 1.0 } else { 0.0 }).collect();
        let (b, _) = ols(&x, &unit)?;
        *w = [b[8], b[9]];
    }
    let cells = lats.len() * lons.len();
    let estimate = |k: usize, c: usize| -> f64 { weights.iter().zip(&grids).map(|(w, g)| w[k] * g[c]).sum() };
    let step: Vec<f64> = (0..cells).map(|c| estimate(0, c)).collect();
    let coord = |c: usize| (lats[c / lons.len()].to_degrees(), lons[c % lons.len()].to_degrees());
    let disk = |centre: (f64, f64)| (0..cells).filter(move |&c| distance_km(coord(c), centre) < RADIUS_KM);
    let stat = |centre: (f64, f64)| disk(centre).map(|c| step[c].abs()).fold(0.0, f64::max);

    // Null centres, drawn before the Kamchatka statistic is computed.
    let rings = land_rings();
    let mut rng = Rng::new(20250729);
    let mut centres = Vec::new();
    while centres.len() < NULL_POINTS {
        let u = |rng: &mut Rng| (rng.next_u64() >> 11) as f64 / (1u64 << 53) as f64;
        let (s0, s1) = (40f64.to_radians().sin(), 65f64.to_radians().sin());
        let lat = (s0 + (s1 - s0) * u(&mut rng)).asin().to_degrees() * if u(&mut rng) < 0.5 { -1.0 } else { 1.0 };
        let lon = -180.0 + 360.0 * u(&mut rng);
        if !on_land(&rings, lon, lat) && distance_km((lat, lon), EPICENTRE) > 2500.0 {
            centres.push((lat, lon));
        }
    }
    let null: Vec<f64> = centres.iter().map(|&c| stat(c)).collect();
    let observed = stat(EPICENTRE);
    let rank_p = (1 + null.iter().filter(|&&v| v >= observed).count()) as f64 / (NULL_POINTS + 1) as f64;

    let near: Vec<usize> = disk(EPICENTRE).collect();
    let best = *near.iter().max_by(|&&a, &&b| step[a].abs().total_cmp(&step[b].abs())).ok_or("empty disk")?;
    let lo = near.iter().map(|&c| step[c]).fold(f64::MAX, f64::min);
    let hi = near.iter().map(|&c| step[c]).fold(f64::MIN, f64::max);
    let (blat, blon) = coord(best);
    let series: Vec<f64> = grids.iter().map(|g| g[best]).collect();
    let (beta, resid) = ols(&x, &series)?;
    let sd = (resid.iter().map(|e| e * e).sum::<f64>() / (resid.len() - beta.len()) as f64).sqrt();

    let mut rep = String::new();
    let _ = writeln!(rep, "GRACE-FO months: {} ({} after the event); event 2025-07-29 23:25 UTC at {:.1}°N {:.1}°E", spans.len(), spans.iter().filter(|s| s.end > te).count(), EPICENTRE.0, EPICENTRE.1);
    let _ = writeln!(rep, "largest |step| within {RADIUS_KM:.0} km: {:.2} µGal at {blat:.1}°N {blon:.1}°E (sign {:+.2}); step range in the disk {lo:.2} to {hi:.2} µGal", observed, step[best]);
    let _ = writeln!(rep, "at that cell: post-event ramp {:+.2} µGal/yr; residual sd {sd:.2} µGal", beta[9]);
    let _ = writeln!(rep, "null ({} ocean points, 40–65° latitude): median {:.2}, 95th percentile {:.2}, max {:.2} µGal", NULL_POINTS, quantile(&null, 0.5), quantile(&null, 0.95), null.iter().fold(0.0, |a: f64, &b| a.max(b)));
    let _ = writeln!(rep, "rank p = {rank_p:.3} ({})", if observed > quantile(&null, 0.95) { "above the null's 95th percentile: detected" } else { "not above the null's 95th percentile: no detection at this filtering" });

    let mut csv = String::from("lat,lon,step_ugal,ramp_ugal_per_yr\n");
    for c in (0..cells).filter(|&c| { let (la, lo) = coord(c); (40.0..66.0).contains(&la) && (140.0..180.0).contains(&lo) }) {
        let (la, lo) = coord(c);
        let _ = writeln!(csv, "{la:.1},{lo:.1},{:.3},{:.3}", step[c], estimate(1, c));
    }
    let fitted: Vec<f64> = x.iter().map(|r| r.iter().zip(&beta).map(|(a, b)| a * b).sum()).collect();
    let svg = figure(&lats, &lons, &step, &rings, &spans, &series, &fitted, te, (blat, blon), &null, observed, &rep);
    Ok(Outcome { report: rep, svg, csv })
}

#[allow(clippy::too_many_arguments)] // one-off figure; the pieces have nothing in common to bundle
fn figure(lats: &[f64], lons: &[f64], step: &[f64], rings: &[Vec<(f64, f64)>], spans: &[Span], series: &[f64], fitted: &[f64], te: f64, best: (f64, f64), null: &[f64], observed: f64, rep: &str) -> String {
    let mut s = open(1000.0, 600.0, "Kamchatka 2025 earthquake in GRACE-FO");
    let _ = write!(s, r#"<text x="40" y="32" class="title">Is the 2025 Kamchatka earthquake visible in GRACE-FO?</text>
<text x="40" y="54" class="subtitle">COST-G monthly gravity, degree ≤ 60, 300 km Gaussian, no destriping. Step at 29 July 2025 fitted with trend, seasons and a post-event ramp.</text>
"#);
    // a. Regional step map, 140–180°E, 40–66°N, equirectangular at 53°N.
    let (px, py) = (8.0, 8.0 / 53f64.to_radians().cos());
    let (mx, my) = (40.0, 100.0);
    let mapx = |lon: f64| mx + (lon - 140.0) * px;
    let mapy = |lat: f64| my + (66.0 - lat) * py;
    let limit = step.iter().enumerate().filter(|(c, _)| {
        let (la, lo) = (lats[c / lons.len()].to_degrees(), lons[c % lons.len()].to_degrees());
        (40.0..66.0).contains(&la) && (140.0..180.0).contains(&lo)
    }).map(|(_, v)| v.abs()).fold(0.5, f64::max);
    let _ = writeln!(s, r#"<text x="{mx}" y="{:.1}" class="panel">a. Step (µGal); colour limit ±{limit:.1}</text>"#, my - 12.0);
    for (i, &la) in lats.iter().enumerate() {
        let la = la.to_degrees();
        if !(40.0..66.0).contains(&la) { continue; }
        for (j, &lo) in lons.iter().enumerate() {
            let lo = lo.to_degrees();
            if !(140.0..180.0).contains(&lo) { continue; }
            let v = step[i * lons.len() + j];
            let _ = writeln!(s, r#"<rect x="{:.1}" y="{:.1}" width="{:.2}" height="{:.2}" fill="var({})" fill-opacity="{:.3}"/>"#,
                mapx(lo - 0.5), mapy(la + 0.5), px + 0.2, py + 0.2, if v >= 0.0 { "--s1" } else { "--s2" }, (v.abs() / limit).min(1.0));
        }
    }
    for ring in rings {
        let pts: Vec<String> = ring.iter().filter(|(lo, la)| (139.0..181.0).contains(lo) && (39.0..67.0).contains(la)).map(|&(lo, la)| format!("{:.1},{:.1}", mapx(lo), mapy(la))).collect();
        if pts.len() > 1 {
            let _ = writeln!(s, r#"<polyline points="{}" fill="none" stroke="var(--ink2)" stroke-width="0.8"/>"#, pts.join(" "));
        }
    }
    let _ = writeln!(s, r#"<rect x="{mx}" y="{my}" width="{:.1}" height="{:.1}" fill="none" class="axis"/>"#, 40.0 * px, 26.0 * py);
    let (ex, ey) = (mapx(EPICENTRE.1), mapy(EPICENTRE.0));
    let r_px = RADIUS_KM / 111.2 * py;
    let rx = RADIUS_KM / (111.2 * EPICENTRE.0.to_radians().cos()) * px;
    let _ = writeln!(s, r#"<ellipse cx="{ex:.1}" cy="{ey:.1}" rx="{rx:.1}" ry="{r_px:.1}" fill="none" stroke="var(--ink)" stroke-dasharray="4 3"/><text x="{:.1}" y="{:.1}" class="tick">500 km from the epicentre (★)</text><text x="{ex:.1}" y="{:.1}" text-anchor="middle" style="fill: var(--ink)">★</text>"#,
        ex - rx, ey - r_px - 6.0, ey + 4.0);
    let _ = writeln!(s, r#"<circle cx="{:.1}" cy="{:.1}" r="4" fill="none" stroke="var(--ink)" stroke-width="1.5"/>"#, mapx(best.1), mapy(best.0));
    for lo in (140..=180).step_by(10) {
        let _ = writeln!(s, r#"<text x="{:.1}" y="{:.1}" class="tick" text-anchor="middle">{lo}°E</text>"#, mapx(lo as f64), my + 26.0 * py + 16.0);
    }
    for la in (40..=65).step_by(5) {
        let _ = writeln!(s, r#"<text x="{:.1}" y="{:.1}" class="tick" text-anchor="end">{la}°N</text>"#, mx - 4.0, mapy(la as f64) + 4.0);
    }

    // b. Series at the circled cell.
    let (bx0, bx1, by0, by1) = (450.0, 960.0, 100.0, 300.0);
    let t: Vec<f64> = spans.iter().map(|s| s.mid()).collect();
    let (t0, t1) = (2018.4, 2026.7);
    let vmax = series.iter().chain(fitted).fold(0.0, |a: f64, &b| a.max(b.abs())).ceil();
    let bx = |v: f64| bx0 + (bx1 - bx0) * (v - t0) / (t1 - t0);
    let byy = |v: f64| by0 + (by1 - by0) * (vmax - v) / (2.0 * vmax);
    let _ = writeln!(s, r#"<text x="{bx0}" y="{:.1}" class="panel">b. Circled cell, {:.1}°N {:.1}°E (µGal)</text>"#, by0 - 12.0, best.0, best.1);
    let _ = writeln!(s, r#"<line x1="{bx0}" y1="{:.1}" x2="{bx1}" y2="{:.1}" class="axis"/><line x1="{:.1}" y1="{by0}" x2="{:.1}" y2="{by1}" class="grid" stroke-dasharray="3 3"/>"#, byy(0.0), byy(0.0), bx(te), bx(te));
    for v in [-vmax, vmax] {
        let _ = writeln!(s, r#"<text x="{:.1}" y="{:.1}" class="tick" text-anchor="end">{v:+.0}</text>"#, bx0 - 6.0, byy(v) + 4.0);
    }
    for yr in (2019..=2026).step_by(1) {
        let _ = writeln!(s, r#"<text x="{:.1}" y="{:.1}" class="tick" text-anchor="middle">{yr}</text>"#, bx(yr as f64), by1 + 16.0);
    }
    let line = |vals: &[f64]| -> String {
        (0..t.len()).map(|i| format!("{}{:.1},{:.1}", if i == 0 || t[i] - t[i - 1] > 45.0 / 365.25 { 'M' } else { 'L' }, bx(t[i]), byy(vals[i]))).collect()
    };
    let _ = writeln!(s, r#"<path d="{}" class="line ink" stroke-width="1.2"/><path d="{}" class="line s1" stroke-width="1.5"/>"#, line(series), line(fitted));
    let _ = writeln!(s, r#"<text x="{:.1}" y="{:.1}" class="tick" text-anchor="end">earthquake</text>"#, bx(te) - 4.0, by0 + 12.0);

    // c. Null histogram.
    let (hx0, hx1, hy0, hy1) = (450.0, 960.0, 380.0, 500.0);
    let top = null.iter().fold(observed, |a: f64, &b| a.max(b)) * 1.05;
    let bins = 25;
    let mut counts = vec![0usize; bins];
    for &v in null {
        counts[((v / top * bins as f64) as usize).min(bins - 1)] += 1;
    }
    let cmax = *counts.iter().max().unwrap_or(&1) as f64;
    let hx = |v: f64| hx0 + (hx1 - hx0) * v / top;
    let _ = writeln!(s, r#"<text x="{hx0}" y="{:.1}" class="panel">c. Largest |step| within 500 km (µGal): Kamchatka against 100 ocean points</text>"#, hy0 - 40.0);
    for (k, &n) in counts.iter().enumerate() {
        let h = (hy1 - hy0) * n as f64 / cmax;
        let _ = writeln!(s, r#"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="{h:.1}" fill="var(--muted)" opacity="0.6"/>"#, hx(k as f64 * top / bins as f64) + 0.5, hy1 - h, (hx1 - hx0) / bins as f64 - 1.0);
    }
    let q95 = quantile(null, 0.95);
    let _ = writeln!(s, r#"<line x1="{hx0}" y1="{hy1}" x2="{hx1}" y2="{hy1}" class="axis"/>"#);
    let _ = writeln!(s, r#"<line x1="{0:.1}" y1="{hy0}" x2="{0:.1}" y2="{hy1}" stroke="var(--ink2)" stroke-dasharray="4 3"/><text x="{0:.1}" y="{1:.1}" class="tick" text-anchor="middle">95%</text>"#, hx(q95), hy0 - 4.0);
    let _ = writeln!(s, r#"<line x1="{0:.1}" y1="{1:.1}" x2="{0:.1}" y2="{hy1}" class="s2" stroke-width="2.5"/><text x="{0:.1}" y="{2:.1}" class="tick" text-anchor="middle" style="fill: var(--s2)">Kamchatka {3:.1}</text>"#, hx(observed), hy0 - 14.0, hy0 - 18.0, observed);
    let ticks = if top > 10.0 { 2.0 } else { 1.0 };
    let mut v = 0.0;
    while v <= top {
        let _ = writeln!(s, r#"<text x="{:.1}" y="{:.1}" class="tick" text-anchor="middle">{v:.0}</text>"#, hx(v), hy1 + 16.0);
        v += ticks;
    }
    let verdict = rep.lines().last().unwrap_or("");
    let _ = write!(s, r#"<text x="40" y="560" class="subtitle">{}</text>
</svg>
"#, esc(verdict));
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn distances_and_mixed_month_step() {
        assert!((distance_km((0.0, 0.0), (0.0, 90.0)) - 10_007.5).abs() < 1.0);
        let te = julian_epoch(2025, 7, 29, 23, 25);
        let july = Span { start: julian_epoch(2025, 7, 1, 0, 0), end: julian_epoch(2025, 8, 1, 0, 0) };
        let r = row(july, te);
        assert!((r[8] - 2.02 / 31.0).abs() < 0.01, "{}", r[8]);
        let june = Span { start: julian_epoch(2025, 6, 1, 0, 0), end: julian_epoch(2025, 7, 1, 0, 0) };
        assert_eq!(row(june, te)[8..], [0.0, 0.0]);
    }
}
