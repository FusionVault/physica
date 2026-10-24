//! E1: the 29 July 2025 Kamchatka Mw 8.8 earthquake in COST-G monthly gravity, against a
//! forward model of the USGS finite-fault solution.
//!
//! Observations: gravity anomaly (µGal) of every GRACE-FO month minus the 2004–2009 mean, degree
//! ≤ 60, 300 km Gaussian, no destriping (it removes north–south structure, which a subduction
//! dipole has), 1° grid. Each cell gets one weighted least-squares fit (weights: the Q0 noise
//! weight of each month): intercept, trend, annual, semi-annual and 161-day terms, a step at the
//! event (July 2025 weighted by the share of its span after the event) and a post-event ramp.
//!
//! Tests, fixed before the first run of this version:
//! - S1, detection: the largest |step| within 500 km of the epicentre;
//! - S2, after-slip: the largest |ramp| within 500 km;
//! - S3, model match: α, the least-squares amplitude of the forward model in the step map within
//!   1,500 km (area-weighted), and the correlation between the two.
//!
//! The null for all three is the same statistic around 100 random ocean points at 40–65°
//! latitude, more than 2,500 km from the epicentre; for S3 the model template is shifted there.
//! Forward model: the USGS `coulomb.inp` subfaults in an Okubo (1992) half-space (2,700 kg/m³,
//! Poisson's ratio 0.25) with the ocean-layer correction (1,025 kg/m³) where the surface is sea,
//! on a 0.1° grid, expanded to degree 60 and filtered exactly like the data.

use std::fmt::Write;
use std::path::Path;

use geopotential::Functional;
use gravimetry::coseismic::{Fault, Medium, read_coulomb};
use geopotential::julian_epoch;
use gravimetry::technotes::{read_tn13, read_tn14};

use crate::shared::{Rng, land_rings, on_land, ols, quantile};
use crate::svg::{esc, open};
use crate::{Result, Span, baseline_mean, corrected_month, gfc_files, grid_axes, span};

const EPICENTRE: (f64, f64) = (52.4948, 160.2395); // USGS us6000qw60
const RADIUS_KM: f64 = 500.0;
const TEMPLATE_KM: f64 = 1500.0;
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

/// Q0's noise weight for each month, matched by mid-epoch.
fn noise_weights(path: &Path, spans: &[Span]) -> Result<Vec<f64>> {
    let text = std::fs::read_to_string(path)?;
    let mut lines = text.lines();
    let head: Vec<&str> = lines.next().ok_or("empty Q0 file")?.split(',').collect();
    let col = |name: &str| head.iter().position(|h| *h == name).ok_or(format!("no {name} in the Q0 file"));
    let (ct, cw) = (col("epoch_year")?, col("noise_weight")?);
    let table: Vec<(f64, f64)> = lines.filter_map(|l| { let f: Vec<&str> = l.split(',').collect(); Some((f.get(ct)?.parse().ok()?, f.get(cw)?.parse().ok()?)) }).collect();
    spans.iter().map(|s| {
        table.iter().find(|(t, _)| (t - s.mid()).abs() < 0.005).map(|(_, w)| *w).ok_or_else(|| format!("no Q0 weight for {:.3}", s.mid()).into())
    }).collect()
}

/// The forward model as a filtered gravity anomaly (µGal) on the analysis grid.
fn forward_model(faults: &[Fault], lats: &[f64], lons: &[f64], gm: f64, radius: f64) -> Vec<f64> {
    let rings = land_rings();
    let fine_lat: Vec<f64> = (0..170).map(|i| (44.05 + 0.1 * i as f64).to_radians()).collect();
    let fine_lon: Vec<f64> = (0..220).map(|j| (150.05 + 0.1 * j as f64).to_radians()).collect();
    let masses = gravimetry::coseismic::surface_masses(faults, &fine_lat, &fine_lon, &Medium::default(), |lat, lon| !on_land(&rings, lon.to_degrees(), lat.to_degrees()));
    let (net, gross) = masses.iter().fold((0.0, 0.0), |(n, g), m| (n + m.2, g + m.2.abs()));
    eprintln!("e1 model: net mass {:.2} Gt of gross {:.2} Gt (the ocean correction moves water away)", net / 1e12, gross / 1e12);
    let mut field = gravimetry::mass_coefficients(&masses, gm, radius, 60);
    for (n, m) in [(0, 0), (1, 0), (1, 1)] {
        field.set(n, m, 0.0, 0.0);
    }
    let field = gravimetry::smooth(&field, &gravimetry::gaussian_weights(300e3, radius, 60));
    field.grid(&Functional::GravityAnomaly, radius, lats, lons).into_iter().map(|g| g * 1e8).collect()
}

/// Our uplift against USGS's `surface_deformation.disp` (lon, lat, elevation, E, N, U):
/// (regression slope of ours on theirs, correlation, number of points).
fn uplift_check(faults: &[Fault], path: &Path) -> Result<(f64, f64, usize)> {
    let text = std::fs::read_to_string(path)?;
    let pts: Vec<(f64, f64, f64)> = text.lines().filter(|l| !l.starts_with('#')).filter_map(|l| {
        let f: Vec<f64> = l.split_whitespace().filter_map(|s| s.parse().ok()).collect();
        (f.len() >= 6).then(|| (f[1], f[0], f[5]))
    }).collect();
    let medium = Medium::default();
    let (mut sxy, mut sxx, mut syy, mut sx, mut sy) = (0.0, 0.0, 0.0, 0.0, 0.0);
    for &(lat, lon, usgs) in &pts {
        let ours: f64 = faults.iter().map(|f| {
            let e = 6_371_000.0 * f.lat.to_radians().cos() * (lon - f.lon).to_radians();
            f.surface_change(e, 6_371_000.0 * (lat - f.lat).to_radians(), &medium).uplift
        }).sum();
        (sxy, sxx, syy, sx, sy) = (sxy + usgs * ours, sxx + usgs * usgs, syy + ours * ours, sx + usgs, sy + ours);
    }
    let n = pts.len() as f64;
    let cov = sxy - sx * sy / n;
    Ok((cov / (sxx - sx * sx / n), cov / ((sxx - sx * sx / n) * (syy - sy * sy / n)).sqrt(), pts.len()))
}

pub fn run(months_dir: &str, baseline_dir: &str, tn14: &str, tn13: &str, usgs_dir: &Path, q0: &Path) -> Result<Outcome> {
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
    let w = noise_weights(q0, &spans)?;
    eprintln!("e1: {} GRACE-FO months, noise weights {:.2} to {:.2}", grids.len(), w.iter().fold(f64::MAX, |a, &b| a.min(b)), w.iter().fold(0.0, |a: f64, &b| a.max(b)));

    // Weighted least squares is OLS on rows scaled by √w. The step and ramp are linear in the
    // data, so one solve per unit vector gives their weights for every cell.
    let x: Vec<Vec<f64>> = spans.iter().map(|&s| row(s, te)).collect();
    let xw: Vec<Vec<f64>> = x.iter().zip(&w).map(|(r, wi)| r.iter().map(|v| v * wi.sqrt()).collect()).collect();
    let mut kernel = vec![[0.0; 2]; x.len()];
    for (j, k) in kernel.iter_mut().enumerate() {
        let unit: Vec<f64> = (0..x.len()).map(|i| if i == j { w[j].sqrt() } else { 0.0 }).collect();
        let (b, _) = ols(&xw, &unit)?;
        *k = [b[8], b[9]];
    }
    let cells = lats.len() * lons.len();
    let estimate = |k: usize, c: usize| -> f64 { kernel.iter().zip(&grids).map(|(w, g)| w[k] * g[c]).sum() };
    let step: Vec<f64> = (0..cells).map(|c| estimate(0, c)).collect();
    let ramp: Vec<f64> = (0..cells).map(|c| estimate(1, c)).collect();
    let coord = |c: usize| (lats[c / lons.len()].to_degrees(), lons[c % lons.len()].to_degrees());
    let cell_at = |lat: f64, lon: f64| -> usize {
        let i = ((lat + 89.5).round() as i64).clamp(0, lats.len() as i64 - 1) as usize;
        let j = ((lon + 179.5).round() as i64).rem_euclid(lons.len() as i64) as usize;
        i * lons.len() + j
    };
    let disk = |centre: (f64, f64)| (0..cells).filter(move |&c| distance_km(coord(c), centre) < RADIUS_KM);
    let max_abs = |map: &[f64], centre: (f64, f64)| disk(centre).map(|c| map[c].abs()).fold(0.0, f64::max);

    // Forward model and its template.
    let faults = read_coulomb(&std::fs::read_to_string(usgs_dir.join("coulomb.inp"))?)?;
    let model = forward_model(&faults, &lats, &lons, mean.gm, mean.radius);
    let template: Vec<(f64, f64, f64, f64)> = (0..cells).filter(|&c| distance_km(coord(c), EPICENTRE) < TEMPLATE_KM)
        .map(|c| { let (la, lo) = coord(c); (la - EPICENTRE.0, lo - EPICENTRE.1, model[c], la.to_radians().cos()) }).collect();
    // α and correlation of the step map around `centre` against the template shifted there.
    let project = |centre: (f64, f64)| -> (f64, f64) {
        let k = EPICENTRE.0.to_radians().cos() / centre.0.to_radians().cos();
        let (mut om, mut mm, mut oo) = (0.0, 0.0, 0.0);
        for &(dla, dlo, m, a) in &template {
            let o = step[cell_at(centre.0 + dla, centre.1 + dlo * k)];
            (om, mm, oo) = (om + a * o * m, mm + a * m * m, oo + a * o * o);
        }
        (om / mm, om / (mm * oo).sqrt())
    };

    // Null centres, drawn as in the quick look.
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
    let null_step: Vec<f64> = centres.iter().map(|&c| max_abs(&step, c)).collect();
    let null_ramp: Vec<f64> = centres.iter().map(|&c| max_abs(&ramp, c)).collect();
    let null_proj: Vec<(f64, f64)> = centres.iter().map(|&c| project(c)).collect();
    let (s1, s2) = (max_abs(&step, EPICENTRE), max_abs(&ramp, EPICENTRE));
    let (alpha, corr) = project(EPICENTRE);
    let rank = |obs: f64, null: &[f64]| (1 + null.iter().filter(|&&v| v >= obs).count()) as f64 / (null.len() + 1) as f64;
    let null_alpha: Vec<f64> = null_proj.iter().map(|p| p.0).collect();
    let alpha_sd = (null_alpha.iter().map(|a| a * a).sum::<f64>() / null_alpha.len() as f64).sqrt();
    let p_alpha = rank(alpha.abs(), &null_alpha.iter().map(|a| a.abs()).collect::<Vec<_>>());
    let p_corr = rank(corr, &null_proj.iter().map(|p| p.1).collect::<Vec<_>>());

    let near: Vec<usize> = disk(EPICENTRE).collect();
    let best = *near.iter().max_by(|&&a, &&b| step[a].abs().total_cmp(&step[b].abs())).ok_or("empty disk")?;
    let range = |map: &[f64], idx: &[usize]| idx.iter().fold((f64::MAX, f64::MIN), |(a, b), &c| (a.min(map[c]), b.max(map[c])));
    let (lo, hi) = range(&step, &near);
    let (mlo, mhi) = range(&model, &near);
    let (blat, blon) = coord(best);
    let series: Vec<f64> = grids.iter().map(|g| g[best]).collect();
    let yw: Vec<f64> = series.iter().zip(&w).map(|(v, wi)| v * wi.sqrt()).collect();
    let (beta, _) = ols(&xw, &yw)?;
    let resid: Vec<f64> = x.iter().zip(&series).map(|(r, y)| y - r.iter().zip(&beta).map(|(a, b)| a * b).sum::<f64>()).collect();
    let sd = (resid.iter().map(|e| e * e).sum::<f64>() / (resid.len() - beta.len()) as f64).sqrt();
    let (slope, r_uplift, n_uplift) = uplift_check(&faults, &usgs_dir.join("surface_deformation.disp"))?;
    let slip_area: f64 = faults.iter().map(|f| f.length * f.width * f.slip).sum();

    let mut rep = String::new();
    let _ = writeln!(rep, "GRACE-FO months: {} ({} after the event), weighted by Q0 noise; event 2025-07-29 23:25 UTC at {:.2}°N {:.2}°E", spans.len(), spans.iter().filter(|s| s.end > te).count(), EPICENTRE.0, EPICENTRE.1);
    let _ = writeln!(rep, "forward model: {} USGS subfaults, Σ area·slip = {:.3e} m³ (= {:.2e} N·m at μ = 40 GPa; USGS M0 2.25e22); our uplift vs USGS at {n_uplift} points: slope {slope:.2}, r = {r_uplift:.3}", faults.len(), slip_area, slip_area * 40e9);
    let _ = writeln!(rep, "S1 detection: largest |step| within {RADIUS_KM:.0} km {s1:.2} µGal at {blat:.1}°N {blon:.1}°E; null 95th percentile {:.2}, max {:.2}; rank p = {:.3}", quantile(&null_step, 0.95), null_step.iter().fold(0.0, |a: f64, &b| a.max(b)), rank(s1, &null_step));
    let _ = writeln!(rep, "  step range within {RADIUS_KM:.0} km: observed {lo:.2} to {hi:.2} µGal; model {mlo:.2} to {mhi:.2} µGal");
    let _ = writeln!(rep, "S2 after-slip: largest |ramp| within {RADIUS_KM:.0} km {s2:.2} µGal/yr (at the S1 cell {:+.2}); null 95th percentile {:.2}, max {:.2}; rank p = {:.3}", beta[9], quantile(&null_ramp, 0.95), null_ramp.iter().fold(0.0, |a: f64, &b| a.max(b)), rank(s2, &null_ramp));
    let _ = writeln!(rep, "S3 model match within {TEMPLATE_KM:.0} km: amplitude α = {alpha:.2} ± {alpha_sd:.2} (null rms), rank p = {p_alpha:.3}; correlation r = {corr:.2}, rank p = {p_corr:.3}");
    let _ = writeln!(rep, "  α > 1: more signal than the coseismic slip explains (early after-slip in the first months, or a stiffer or denser real Earth); α < 1: less.");
    let _ = writeln!(rep, "residual sd at the S1 cell {sd:.2} µGal");

    let in_region = |c: usize| { let (la, lo) = coord(c); (40.0..66.0).contains(&la) && (140.0..180.0).contains(&lo) };
    let mut csv = String::from("lat,lon,step_ugal,ramp_ugal_per_yr,model_ugal\n");
    for c in (0..cells).filter(|&c| in_region(c)) {
        let (la, lo) = coord(c);
        let _ = writeln!(csv, "{la:.1},{lo:.1},{:.3},{:.3},{:.3}", step[c], ramp[c], model[c]);
    }
    let fitted: Vec<f64> = x.iter().map(|r| r.iter().zip(&beta).map(|(a, b)| a * b).sum()).collect();
    let svg = figure(&Maps { lats: &lats, lons: &lons, step: &step, model: &model }, &rings, &spans, &series, &fitted, te, (blat, blon), &null_alpha, alpha, &rep);

    Ok(Outcome { report: rep, svg, csv })
}

struct Maps<'a> {
    lats: &'a [f64],
    lons: &'a [f64],
    step: &'a [f64],
    model: &'a [f64],
}

#[allow(clippy::too_many_arguments)] // one-off figure; the pieces have nothing in common to bundle
fn figure(maps: &Maps, rings: &[Vec<(f64, f64)>], spans: &[Span], series: &[f64], fitted: &[f64], te: f64, best: (f64, f64), null_alpha: &[f64], alpha: f64, rep: &str) -> String {
    let mut s = open(1000.0, 860.0, "Kamchatka 2025 earthquake in GRACE-FO");
    let _ = write!(s, r#"<text x="40" y="32" class="title">The 2025 Kamchatka earthquake in GRACE-FO, against the USGS fault model</text>
<text x="40" y="54" class="subtitle">COST-G monthly gravity, degree ≤ 60, 300 km Gaussian, no destriping. Step at 29 July 2025, fitted with trend, seasons and a post-event ramp.</text>
"#);
    let (px, py) = (8.0, 8.0 / 53f64.to_radians().cos());
    let limit = maps.lats.iter().enumerate().flat_map(|(i, la)| maps.lons.iter().enumerate().map(move |(j, lo)| (i, j, la.to_degrees(), lo.to_degrees())))
        .filter(|&(_, _, la, lo)| (40.0..66.0).contains(&la) && (140.0..180.0).contains(&lo))
        .map(|(i, j, _, _)| maps.step[i * maps.lons.len() + j].abs()).fold(0.5, f64::max);
    for (k, (title, map, mx)) in [("a. Observed step (µGal)", maps.step, 40.0), ("b. USGS fault model, same filter (µGal)", maps.model, 520.0)].into_iter().enumerate() {
        let my = 100.0;
        let mapx = |lon: f64| mx + (lon - 140.0) * px;
        let mapy = |lat: f64| my + (66.0 - lat) * py;
        let _ = writeln!(s, r#"<text x="{mx}" y="{:.1}" class="panel">{}</text>"#, my - 12.0, esc(title));
        for (i, &la) in maps.lats.iter().enumerate() {
            let la = la.to_degrees();
            if !(40.0..66.0).contains(&la) { continue; }
            for (j, &lo) in maps.lons.iter().enumerate() {
                let lo = lo.to_degrees();
                if !(140.0..180.0).contains(&lo) { continue; }
                let v = map[i * maps.lons.len() + j];
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
        let _ = writeln!(s, r#"<text x="{ex:.1}" y="{:.1}" text-anchor="middle" style="fill: var(--ink)">★</text>"#, ey + 4.0);
        if k == 0 {
            let _ = writeln!(s, r#"<circle cx="{:.1}" cy="{:.1}" r="4" fill="none" stroke="var(--ink)" stroke-width="1.5"/>"#, mapx(best.1), mapy(best.0));
        }
        for lo in (140..=180).step_by(10) {
            let _ = writeln!(s, r#"<text x="{:.1}" y="{:.1}" class="tick" text-anchor="middle">{lo}°E</text>"#, mapx(lo as f64), my + 26.0 * py + 16.0);
        }
        for la in (40..=65).step_by(5) {
            let _ = writeln!(s, r#"<text x="{:.1}" y="{:.1}" class="tick" text-anchor="end">{la}°N</text>"#, mx - 4.0, mapy(la as f64) + 4.0);
        }
    }
    let _ = writeln!(s, r#"<text x="40" y="490" class="subtitle">Colour limit ±{limit:.1} µGal on both maps; blue more mass, orange less. ★ epicentre; circle: the cell charted in c.</text>"#);

    // c. Series at the circled cell.
    let (bx0, bx1, by0, by1) = (60.0, 520.0, 560.0, 760.0);
    let t: Vec<f64> = spans.iter().map(|s| s.mid()).collect();
    let (t0, t1) = (2018.4, 2026.7);
    let vmax = series.iter().chain(fitted).fold(0.0, |a: f64, &b| a.max(b.abs())).ceil();
    let bx = |v: f64| bx0 + (bx1 - bx0) * (v - t0) / (t1 - t0);
    let byy = |v: f64| by0 + (by1 - by0) * (vmax - v) / (2.0 * vmax);
    let _ = writeln!(s, r#"<text x="{bx0}" y="{:.1}" class="panel">c. Circled cell, {:.1}°N {:.1}°E (µGal); blue: fit</text>"#, by0 - 14.0, best.0, best.1);
    let _ = writeln!(s, r#"<line x1="{bx0}" y1="{:.1}" x2="{bx1}" y2="{:.1}" class="axis"/><line x1="{:.1}" y1="{by0}" x2="{:.1}" y2="{by1}" class="grid" stroke-dasharray="3 3"/>"#, byy(0.0), byy(0.0), bx(te), bx(te));
    for v in [-vmax, vmax] {
        let _ = writeln!(s, r#"<text x="{:.1}" y="{:.1}" class="tick" text-anchor="end">{v:+.0}</text>"#, bx0 - 6.0, byy(v) + 4.0);
    }
    for yr in 2019..=2026 {
        let _ = writeln!(s, r#"<text x="{:.1}" y="{:.1}" class="tick" text-anchor="middle">{yr}</text>"#, bx(yr as f64), by1 + 16.0);
    }
    let line = |vals: &[f64]| -> String {
        (0..t.len()).map(|i| format!("{}{:.1},{:.1}", if i == 0 || t[i] - t[i - 1] > 45.0 / 365.25 { 'M' } else { 'L' }, bx(t[i]), byy(vals[i]))).collect()
    };
    let _ = writeln!(s, r#"<path d="{}" class="line ink" stroke-width="1.2"/><path d="{}" class="line s1" stroke-width="1.5"/>"#, line(series), line(fitted));

    // d. Null of α.
    let (hx0, hx1, hy0, hy1) = (580.0, 960.0, 580.0, 760.0);
    let top = null_alpha.iter().fold(alpha.abs(), |a: f64, &b| a.max(b.abs())) * 1.1;
    let bins = 24;
    let mut counts = vec![0usize; bins];
    for &v in null_alpha {
        counts[(((v + top) / (2.0 * top) * bins as f64) as usize).min(bins - 1)] += 1;
    }
    let cmax = *counts.iter().max().unwrap_or(&1) as f64;
    let hx = |v: f64| hx0 + (hx1 - hx0) * (v + top) / (2.0 * top);
    let _ = writeln!(s, r#"<text x="{hx0}" y="{:.1}" class="panel">d. Model amplitude α: Kamchatka against 100 ocean points</text>"#, hy0 - 34.0);
    for (k, &n) in counts.iter().enumerate() {
        let h = (hy1 - hy0) * n as f64 / cmax;
        let _ = writeln!(s, r#"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="{h:.1}" fill="var(--muted)" opacity="0.6"/>"#, hx(-top + k as f64 * 2.0 * top / bins as f64) + 0.5, hy1 - h, (hx1 - hx0) / bins as f64 - 1.0);
    }
    let _ = writeln!(s, r#"<line x1="{hx0}" y1="{hy1}" x2="{hx1}" y2="{hy1}" class="axis"/>"#);
    let _ = writeln!(s, r#"<line x1="{0:.1}" y1="{1:.1}" x2="{0:.1}" y2="{hy1}" class="s2" stroke-width="2.5"/><text x="{0:.1}" y="{2:.1}" class="tick" text-anchor="middle" style="fill: var(--s2)">Kamchatka α = {alpha:.2}</text>"#, hx(alpha), hy0 - 10.0, hy0 - 14.0);
    for v in [-top, 0.0, top] {
        let _ = writeln!(s, r#"<text x="{:.1}" y="{:.1}" class="tick" text-anchor="middle">{v:.1}</text>"#, hx(v), hy1 + 16.0);
    }
    for (i, line) in rep.lines().filter(|l| l.starts_with('S')).enumerate() {
        let _ = writeln!(s, r#"<text x="40" y="{:.1}" class="subtitle">{}</text>"#, 800.0 + 18.0 * i as f64, esc(line.split(';').next().unwrap_or(line)));
    }
    s.push_str("</svg>\n");
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
