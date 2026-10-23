//! Figures from the quantum-gravimeter crates.
//!
//! ```text
//! gravviz noise OUT_DIR
//! gravviz water-map MONTH.gfc BASELINE_DIR TN-14.txt TN-13.txt OUT_DIR
//! gravviz gradiometry STATIC_MODEL.gfc OUT_DIR
//! ```
//!
//! `noise` writes `noise_spectra.svg` and `.csv` (Zingerle et al. 2024, Fig. 1). `water-map`
//! writes `water_<yyyy-mm>.svg` and a 1° `.csv`: the month's equivalent water height minus the
//! mean of every `.gfc` in BASELINE_DIR, all from ICGEM GSM-2 files (`GSM-2_YYYYDDD-YYYYDDD_…`).
//! `gradiometry` writes `gradiometry_errors.svg` and `.csv`: gradiometer mission errors against
//! the static model's signal (e.g. GOCO06s).

mod gradiometry;
mod map;
mod noise;
mod svg;

use std::fmt::Write as _;
use std::path::Path;

use geopotential::{Coefficients, Functional};
use gravity_field::{Degree1, SlrC20C30, julian_epoch, read_gfc, read_tn13, read_tn14};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match args.iter().map(String::as_str).collect::<Vec<_>>().as_slice() {
        ["noise", out] => {
            let (svg, csv) = noise::figure();
            write(out, "noise_spectra.svg", &svg)?;
            write(out, "noise_spectra.csv", &csv)
        }
        ["water-map", month, baseline, tn14, tn13, out] => water_map(month, baseline, tn14, tn13, out),
        ["gradiometry", model, out] => {
            let gfc = read_gfc(&std::fs::read_to_string(model)?)?;
            let name = gfc.name.clone().unwrap_or_else(|| "model".into());
            let (svg, csv) = gradiometry::figure(&gfc.model.base, &name);
            write(out, "gradiometry_errors.svg", &svg)?;
            write(out, "gradiometry_errors.csv", &csv)
        }
        _ => Err("usage: gravviz noise OUT_DIR | gravviz water-map MONTH.gfc BASELINE_DIR TN-14.txt TN-13.txt OUT_DIR | gravviz gradiometry STATIC_MODEL.gfc OUT_DIR".into()),
    }
}

fn write(dir: &str, name: &str, content: &str) -> Result<()> {
    std::fs::create_dir_all(dir)?;
    let path = Path::new(dir).join(name);
    std::fs::write(&path, content)?;
    eprintln!("wrote {}", path.display());
    Ok(())
}

/// Start of the span in an ICGEM GSM-2 file name, `GSM-2_YYYYDDD-YYYYDDD_…`, as a Julian epoch,
/// plus its calendar (year, month) for labels.
fn span_start(path: &Path) -> Result<(f64, i32, u32)> {
    let name = path.file_name().and_then(|n| n.to_str()).unwrap_or_default();
    let digits = name.strip_prefix("GSM-2_").and_then(|r| r.get(..7)).ok_or_else(|| format!("{name}: not a GSM-2_YYYYDDD file name"))?;
    let (year, doy): (i32, u32) = (digits[..4].parse()?, digits[4..].parse()?);
    let t = julian_epoch(year, 1, 1, 0, 0) + f64::from(doy - 1) / 365.25;
    let month = (1..=12).rev().find(|&m| julian_epoch(year, m, 1, 0, 0) <= t + 1e-9).unwrap_or(1);
    Ok((t, year, month))
}

/// One monthly field with the standard replacements: `C̄20` (and `C̄30` from August 2016) from
/// TN-14 satellite laser ranging, degree 1 from TN-13.
fn corrected_month(path: &Path, tn14: &[SlrC20C30], tn13: &[Degree1]) -> Result<Coefficients> {
    let (t, ..) = span_start(path)?;
    let mut c = read_gfc(&std::fs::read_to_string(path)?)?.model.base;
    let near = |start: f64| (start - t).abs() < 10.0 / 365.25;
    let slr = tn14.iter().find(|r| near(r.start)).ok_or_else(|| format!("{}: no TN-14 month", path.display()))?;
    c.set(2, 0, slr.c20, 0.0);
    if let (Some(c30), true) = (slr.c30, t >= julian_epoch(2016, 8, 1, 0, 0)) {
        c.set(3, 0, c30, 0.0);
    }
    let d1 = tn13.iter().find(|r| near(r.start)).ok_or_else(|| format!("{}: no TN-13 month", path.display()))?;
    c.set(1, 0, d1.c10, 0.0);
    c.set(1, 1, d1.c11, d1.s11);
    Ok(c)
}

fn water_map(month: &str, baseline: &str, tn14: &str, tn13: &str, out: &str) -> Result<()> {
    let tn14 = read_tn14(&std::fs::read_to_string(tn14)?)?;
    let tn13 = read_tn13(&std::fs::read_to_string(tn13)?)?;
    let month_path = Path::new(month);
    let (_, year, mon) = span_start(month_path)?;
    let field = corrected_month(month_path, &tn14, &tn13)?;

    let mut paths: Vec<_> = std::fs::read_dir(baseline)?.filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.extension().is_some_and(|x| x == "gfc")).collect();
    paths.sort();
    let fields = paths.iter().map(|p| corrected_month(p, &tn14, &tn13)).collect::<Result<Vec<_>>>()?;
    let (first, last) = (span_start(&paths[0])?, span_start(paths.last().ok_or("empty baseline")?)?);
    let mean = gravimetry::mean(&fields).ok_or("empty baseline")?;
    eprintln!("baseline: {} months, {}-{:02} to {}-{:02}", fields.len(), first.1, first.2, last.1, last.2);

    // The GRACE Tellus chain (Landerer & Swenson 2012): degree 60, destriping, 300 km Gaussian.
    let nmax = 60;
    let anomaly = gravimetry::destripe(&gravimetry::anomaly(&field, &mean).truncated(nmax));
    let anomaly = gravimetry::smooth(&anomaly, &gravimetry::gaussian_weights(300e3, field.radius, nmax));
    let love = gravimetry::love_numbers(nmax).ok_or("degree above the Love-number table")?;
    let lats: Vec<f64> = (0..360).map(|i| (-89.75 + 0.5 * i as f64).to_radians()).collect();
    let lons: Vec<f64> = (0..720).map(|j| (-179.75 + 0.5 * j as f64).to_radians()).collect();
    let cm: Vec<f64> = anomaly.grid(&Functional::WaterHeight { love: &love }, field.radius, &lats, &lons)
        .into_iter().map(|m| m * 100.0).collect();

    // Colour limit: the area-weighted 95th percentile of |value|, rounded up to 5 cm.
    let mut weighted: Vec<(f64, f64)> = cm.iter().enumerate().map(|(k, v)| (v.abs(), lats[k / 720].cos())).collect();
    weighted.sort_by(|a, b| a.0.total_cmp(&b.0));
    let total: f64 = weighted.iter().map(|w| w.1).sum();
    let mut acc = 0.0;
    let p95 = weighted.iter().find(|w| { acc += w.1; acc >= 0.95 * total }).map_or(1.0, |w| w.0);
    let limit = (p95 / 5.0).ceil().max(1.0) * 5.0;
    let rms = (cm.iter().enumerate().map(|(k, v)| v * v * lats[k / 720].cos()).sum::<f64>() / total).sqrt();
    let (min, max) = cm.iter().fold((f64::MAX, f64::MIN), |(a, b), &v| (a.min(v), b.max(v)));
    eprintln!("EWH: area-weighted RMS {rms:.1} cm, min {min:.1}, max {max:.1}, 95th pct |EWH| {p95:.1} → colour limit ±{limit} cm");

    let coast_text = include_str!("../data/coastline_110m.txt");
    let coast: Vec<Vec<(f64, f64)>> = coast_text.split("\n\n").map(|seg| seg.lines().filter_map(|l| {
        let mut it = l.split_whitespace().map(str::parse::<f64>);
        Some((it.next()?.ok()?, it.next()?.ok()?))
    }).collect()).filter(|s: &Vec<_>| s.len() > 1).collect();

    let month_name = ["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"][mon as usize - 1];
    let title = format!("Water storage anomaly, {month_name} {year}");
    let sub1 = format!("Equivalent water height of the GRACE-FO {month_name} {year} field relative to the {}–{} mean ({} months); degree 60, destriped, 300 km Gaussian.", first.1, last.1, fields.len());
    let sub2 = "COST-G RL02.1 monthly fields (ICGEM); C₂₀ and C₃₀ from satellite laser ranging (TN-14), degree 1 from TN-13. Glacial isostatic adjustment not removed.";
    let svg = map::render(&map::Grid { lats: &lats, lons: &lons, values: &cm }, limit, &coast, &map::MapText {
        title: &title,
        subtitle: &[&sub1, sub2],
        unit: "cm",
        ends: ("less water", "more water"),
        footnote: "Equal Earth projection. Coastlines: Natural Earth (public domain). Data: COST-G (CC BY 4.0, doi:10.5880/COST-G.ICGEM_02_L2); NASA GRACE technical notes 13 and 14.",
    })?;
    let stem = format!("water_{year}-{mon:02}");
    write(out, &format!("{stem}.svg"), &svg)?;
    let mut csv = String::from("lat_deg,lon_deg,ewh_cm\n");
    for i in (0..360).step_by(2) {
        for j in (0..720).step_by(2) {
            let _ = writeln!(csv, "{:.2},{:.2},{:.2}", lats[i].to_degrees(), lons[j].to_degrees(), cm[i * 720 + j]);
        }
    }
    write(out, &format!("{stem}.csv"), &csv)
}
