//! Figures and viewer data from the quantum-gravimeter crates.
//!
//! ```text
//! gravviz noise OUT_DIR
//! gravviz geoid STATIC_MODEL.gfc OUT_DIR
//! gravviz water-map MONTH.gfc BASELINE_DIR TN-14.txt TN-13.txt OUT_DIR
//! gravviz series MONTHS_DIR BASELINE_DIR TN-14.txt TN-13.txt OUT_DIR
//! gravviz triangle MONTH.gfc OUT_DIR
//! gravviz gradiometry STATIC_MODEL.gfc OUT_DIR
//! gravviz sst OUT_DIR [REFERENCE_MONTH.gfc]
//! gravviz w1 MONTHS_DIR BASELINE_DIR TN-14.txt TN-13.txt NINO34.data OUT_DIR
//! gravviz e1 MONTHS_DIR BASELINE_DIR TN-14.txt TN-13.txt USGS_FAULT_DIR Q0_NOISE.csv OUT_DIR
//! gravviz screen W1_LAND_STORAGE.csv INDEX_DIR OUT_DIR
//! gravviz q0 MONTHS_DIR BASELINE_DIR TN-14.txt TN-13.txt KP_FILE OUT_DIR
//! gravviz kalman MONTHS_DIR BASELINE_DIR TN-14.txt TN-13.txt OUT_DIR VIEWER_DATA_DIR
//! gravviz m1 MONTHS_DIR BASELINE_DIR TN-14.txt TN-13.txt KALMAN_DEGREES.csv CATALOGUE.csv OUT_DIR VIEWER_DATA_DIR
//! ```
//!
//! `noise` writes `noise_spectra.svg` and `.csv` (Zingerle et al. 2024, Fig. 1). `geoid` maps a
//! static model's geoid against the GRS80 ellipsoid (`geoid_<model>.svg`). `water-map`
//! writes `water_<yyyy-mm>.svg` and a 1° `.csv`: the month's equivalent water height minus the
//! mean of every `.gfc` in BASELINE_DIR, all from ICGEM GSM-2 files (`GSM-2_YYYYDDD-YYYYDDD_…`).
//! `series` does the same for every month in MONTHS_DIR on a 2° grid and writes `tws_grids.json`
//! (zlib-compressed little-endian i16 millimetres, `[month][lat][lon]`, base64 in JSON because
//! artifact hosting serves no raw binary), `tws_index.json`, `flow_grids.json` (the minimum-norm
//! transport potential of each month's change, i16 × `flow_scale` Gt/yr) and the coastlines
//! (`coast.txt`) for the browser viewer in `viewer/`. `triangle` writes
//! `error_triangle_<yyyy-mm>.svg`: the month's calibrated coefficient errors as water height.
//! `gradiometry` writes `gradiometry_errors.svg` and `.csv`: gradiometer mission errors against
//! the static model's signal (e.g. GOCO06s). `sst` runs the closed-loop simulation of four
//! GRACE-type mission classes and writes `sst_closed_loop.svg` and `.csv`, with a real month's
//! calibrated errors as reference if given.
//!
//! The research tests: `w1` runs the pre-registered El Niño test on global land storage (NOAA PSL
//! Niño 3.4 anomaly file) and writes `w1_report.txt`, `w1_land_storage.svg` and `.csv`,
//! `w1_forecast.csv` and `w1_monthly.csv`. `e1` looks for the 2025 Kamchatka earthquake in
//! GRACE-FO gravity against a random-location null and the USGS fault model (`coulomb.inp` and
//! `surface_deformation.disp` in USGS_FAULT_DIR), weighting months by Q0, and writes `e1_kamchatka.svg`,
//! `e1_report.txt`, `e1_step.csv` and `e1_frames.json`. `screen` reruns the W1 level model with
//! every NOAA PSL index in INDEX_DIR, QBO and solar flux as negative controls, and writes
//! `screen_indices.svg`, `.txt` and `.csv`. `q0` regresses each month's noise on solar flux, Ap,
//! mission phase and centre weights (GFZ `Kp_ap_Ap_SN_F107_since_1932.txt`) and writes
//! `q0_noise.svg` and `.csv`, `q0_effects.csv` and `q0_report.txt`. `kalman` smooths every
//! coefficient through time (M3), tests COST-G's stated errors (Q2), and writes `kalman.svg`,
//! `kalman_report.txt` and `kalman_degrees.csv`, plus gap-free grids, uncertainty and flow for
//! the viewer. `m1` scans every month for sudden mass changes anywhere (a likelihood-ratio test
//! with Q2's error scaling), scores the detections against a USGS earthquake catalogue (CSV from
//! the FDSN service), and writes `m1_detector.svg`, `m1_report.txt` and the viewer's `m1_scan.json`.

mod claims;
mod e1;
mod gradiometry;
mod kalman;
mod lab;
mod m1;
mod map;
mod noise;
mod q0;
mod screen;
mod sst;
mod svg;
mod t1;
mod tonga;
mod triangle;
mod w1;

use std::fmt::Write as _;
use std::path::{Path, PathBuf};

use geopotential::{Coefficients, Functional};
use gravity_field::{Degree1, SlrC20C30, julian_epoch, read_gfc, read_tn13, read_tn14};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

const MONTH_NAMES: [&str; 12] =
    ["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"];

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match args.iter().map(String::as_str).collect::<Vec<_>>().as_slice() {
        ["noise", out] => {
            let (svg, csv) = noise::figure();
            write(out, "noise_spectra.svg", &svg)?;
            write(out, "noise_spectra.csv", &csv)
        }
        ["geoid", model, out] => geoid(model, out),
        ["water-map", month, baseline, tn14, tn13, out] => water_map(month, baseline, tn14, tn13, out),
        ["series", months, baseline, tn14, tn13, out] => series(months, baseline, tn14, tn13, out),
        ["triangle", month, out] => {
            let gfc = read_gfc(&std::fs::read_to_string(month)?)?;
            let sigma = gfc.sigma.ok_or("the file has no error columns")?;
            let (year, mon) = calendar(span(Path::new(month))?.mid());
            let label = format!("{} {year}", MONTH_NAMES[mon as usize - 1]);
            write(out, &format!("error_triangle_{year}-{mon:02}.svg"), &triangle::figure(&sigma, &label)?)
        }
        ["gradiometry", model, out] => {
            let gfc = read_gfc(&std::fs::read_to_string(model)?)?;
            let name = gfc.name.clone().unwrap_or_else(|| "model".into());
            let (svg, csv) = gradiometry::figure(&gfc.model.base, &name);
            write(out, "gradiometry_errors.svg", &svg)?;
            write(out, "gradiometry_errors.csv", &csv)
        }
        ["sst", out, rest @ ..] if rest.len() <= 1 => {
            let reference = match rest.first() {
                Some(path) => {
                    let sigma = read_gfc(&std::fs::read_to_string(path)?)?.sigma.ok_or("the reference file has no error columns")?;
                    let (year, mon) = calendar(span(Path::new(path))?.mid());
                    Some((sigma, format!("{} {year}", MONTH_NAMES[mon as usize - 1])))
                }
                None => None,
            };
            let (svg, csv) = sst::figure(reference.as_ref().map(|(s, l)| (s, l.as_str())));
            write(out, "sst_closed_loop.svg", &svg)?;
            write(out, "sst_closed_loop.csv", &csv)
        }
        ["w1", months, baseline, tn14, tn13, nino, out] => {
            let (lats, lons, fields) = lab::load_months(months, baseline, tn14, tn13)?;
            let w1::Outcome { report, svg, csv, forecast_csv, monthly_csv } = w1::run(&fields, &lats, &lons, Path::new(nino))?;
            print!("{report}");
            write(out, "w1_report.txt", &report)?;
            write(out, "w1_land_storage.svg", &svg)?;
            write(out, "w1_land_storage.csv", &csv)?;
            write(out, "w1_forecast.csv", &forecast_csv)?;
            write(out, "w1_monthly.csv", &monthly_csv)
        }
        ["e1", months, baseline, tn14, tn13, usgs, q0, out] => {
            let e1::Outcome { report, svg, csv, frames } = e1::run(months, baseline, tn14, tn13, Path::new(usgs), Path::new(q0))?;
            print!("{report}");
            write(out, "e1_report.txt", &report)?;
            write(out, "e1_kamchatka.svg", &svg)?;
            write(out, "e1_step.csv", &csv)?;
            write(out, "e1_frames.json", &frames)
        }
        ["screen", csv, indices, out] => {
            let screen::Outcome { report, svg, csv } = screen::run(Path::new(csv), Path::new(indices))?;
            print!("{report}");
            write(out, "screen_indices.txt", &report)?;
            write(out, "screen_indices.svg", &svg)?;
            write(out, "screen_indices.csv", &csv)
        }
        ["q0", months, baseline, tn14, tn13, kp, out] => {
            let q0::Outcome { report, svg, csv, effects } = q0::run(months, baseline, tn14, tn13, Path::new(kp))?;
            print!("{report}");
            write(out, "q0_report.txt", &report)?;
            write(out, "q0_noise.svg", &svg)?;
            write(out, "q0_noise.csv", &csv)?;
            write(out, "q0_effects.csv", &effects)
        }
        ["kalman", months, baseline, tn14, tn13, out, viewer] => {
            print!("{}", kalman::run(months, baseline, tn14, tn13, out, viewer)?);
            Ok(())
        }
        ["claims", out] => {
            for (name, svg) in claims::figures() {
                write(out, name, &svg)?;
            }
            Ok(())
        }
        ["t1-noise", l1b_dir, aod_dir, goco, costg, out, dates @ ..] if !dates.is_empty() => {
            let (report, svg) = t1::noise_study(Path::new(l1b_dir), Path::new(aod_dir), Path::new(goco), Path::new(costg), dates)?;
            print!("{report}");
            write(out, "t1_noise.txt", &report)?;
            write(out, "t1_noise.svg", &svg)
        }
        ["t1-test", l1b_dir, aod_dir, goco, costg, tides, out, rest @ ..] if !rest.is_empty() => {
            let unblind = rest[0] == "--unblind";
            t1::UNBLIND.store(unblind, std::sync::atomic::Ordering::Relaxed);
            let dates = if unblind { &rest[1..] } else { rest };
            let report = t1::test(Path::new(l1b_dir), Path::new(aod_dir), Path::new(goco), Path::new(costg), Path::new(tides), dates)?;
            print!("{report}");
            write(out, if unblind { "t1_result.txt" } else { "t1_null.txt" }, &report)
        }
        ["tonga", out] => {
            let (report, svg) = tonga::run()?;
            print!("{report}");
            write(out, "tonga_feasibility.txt", &report)?;
            write(out, "tonga_lamb_wave.svg", &svg)
        }
        ["m1", months, baseline, tn14, tn13, kalman, catalogue, out, viewer] => {
            print!("{}", m1::run(months, baseline, tn14, tn13, Path::new(kalman), Path::new(catalogue), out, viewer)?);
            Ok(())
        }
        _ => Err("usage: gravviz noise OUT_DIR | water-map MONTH.gfc BASELINE_DIR TN-14.txt TN-13.txt OUT_DIR | \
            series MONTHS_DIR BASELINE_DIR TN-14.txt TN-13.txt OUT_DIR | triangle MONTH.gfc OUT_DIR | \
            gradiometry STATIC_MODEL.gfc OUT_DIR | geoid STATIC_MODEL.gfc OUT_DIR | sst OUT_DIR [REFERENCE_MONTH.gfc] | \
            w1 MONTHS_DIR BASELINE_DIR TN-14.txt TN-13.txt NINO34.data OUT_DIR | \
            e1 MONTHS_DIR BASELINE_DIR TN-14.txt TN-13.txt USGS_FAULT_DIR Q0_NOISE.csv OUT_DIR | screen W1_LAND_STORAGE.csv INDEX_DIR OUT_DIR | \
            q0 MONTHS_DIR BASELINE_DIR TN-14.txt TN-13.txt KP_FILE OUT_DIR | \
            kalman MONTHS_DIR BASELINE_DIR TN-14.txt TN-13.txt OUT_DIR VIEWER_DATA_DIR | \
            m1 MONTHS_DIR BASELINE_DIR TN-14.txt TN-13.txt KALMAN_DEGREES.csv CATALOGUE.csv OUT_DIR VIEWER_DATA_DIR | tonga OUT_DIR | claims OUT_DIR | \
            t1-noise L1B_DIR AOD1B_DIR GOCO06s.gfc COSTG_2022-01.gfc OUT_DIR DATE... | \
            t1-test L1B_DIR AOD1B_DIR GOCO06s.gfc COSTG_2022-01.gfc GOT_TIDE_DIR OUT_DIR [--unblind] DATE...".into()),
    }
}

fn write(dir: &str, name: &str, content: impl AsRef<[u8]>) -> Result<()> {
    std::fs::create_dir_all(dir)?;
    let path = Path::new(dir).join(name);
    std::fs::write(&path, content)?;
    eprintln!("wrote {}", path.display());
    Ok(())
}

/// The span of an ICGEM GSM-2 file, from its name `GSM-2_YYYYDDD-YYYYDDD_…`, in Julian years.
#[derive(Clone, Copy)]
struct Span {
    start: f64,
    end: f64,
}

impl Span {
    fn mid(&self) -> f64 {
        (self.start + self.end) / 2.0
    }
}

fn span(path: &Path) -> Result<Span> {
    let name = path.file_name().and_then(|n| n.to_str()).unwrap_or_default();
    let bad = || format!("{name}: not a GSM-2_YYYYDDD-YYYYDDD file name");
    let digits = name.strip_prefix("GSM-2_").and_then(|r| r.get(..15)).ok_or_else(bad)?;
    let day = |s: &str| -> Result<f64> {
        let (year, doy): (i32, u32) = (s[..4].parse()?, s[4..].parse()?);
        Ok(julian_epoch(year, 1, 1, 0, 0) + f64::from(doy - 1) / 365.25)
    };
    Ok(Span { start: day(&digits[..7])?, end: day(&digits[8..15])? + 1.0 / 365.25 })
}

/// Calendar (year, month) containing Julian epoch `t`.
fn calendar(t: f64) -> (i32, u32) {
    let y = t.floor() as i32;
    let year = (y - 1..=y + 1).rev().find(|&y| julian_epoch(y, 1, 1, 0, 0) <= t).unwrap_or(y);
    let month = (1..=12).rev().find(|&m| julian_epoch(year, m, 1, 0, 0) <= t).unwrap_or(1);
    (year, month)
}

/// One monthly field with the standard replacements: `C̄20` (and `C̄30` from August 2016) from
/// TN-14 satellite laser ranging, degree 1 from TN-13.
fn corrected_month(path: &Path, tn14: &[SlrC20C30], tn13: &[Degree1]) -> Result<Coefficients> {
    let sp = span(path)?;
    let (t, mid) = (sp.start, sp.mid());
    let mut c = read_gfc(&std::fs::read_to_string(path)?)?.model.base;
    // Technical-note months can be offset from the field's (late-2016 GRACE months ran from the
    // 11th), so take the entry whose midpoint is nearest, within half a month.
    let nearest = |mids: &mut dyn Iterator<Item = (usize, f64)>| {
        mids.map(|(i, m)| (i, (m - mid).abs())).filter(|(_, d)| *d < 15.0 / 365.25).min_by(|a, b| a.1.total_cmp(&b.1)).map(|(i, _)| i)
    };
    let slr = nearest(&mut tn14.iter().map(|r| (r.start + r.end) / 2.0).enumerate())
        .map(|i| &tn14[i]).ok_or_else(|| format!("{}: no TN-14 month", path.display()))?;
    c.set(2, 0, slr.c20, 0.0);
    if let (Some(c30), true) = (slr.c30, t >= julian_epoch(2016, 8, 1, 0, 0)) {
        c.set(3, 0, c30, 0.0);
    }
    let d1 = nearest(&mut tn13.iter().map(|r| (r.start + r.end) / 2.0).enumerate())
        .map(|i| &tn13[i]).ok_or_else(|| format!("{}: no TN-13 month", path.display()))?;
    c.set(1, 0, d1.c10, 0.0);
    c.set(1, 1, d1.c11, d1.s11);
    Ok(c)
}

fn gfc_files(dir: &str) -> Result<Vec<PathBuf>> {
    let mut paths: Vec<_> = std::fs::read_dir(dir)?.filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.extension().is_some_and(|x| x == "gfc")).collect();
    paths.sort();
    Ok(paths)
}

/// The corrected mean of every month in `dir`, with the calendar months it spans and its size.
fn baseline_mean(dir: &str, tn14: &[SlrC20C30], tn13: &[Degree1]) -> Result<(Coefficients, String)> {
    let paths = gfc_files(dir)?;
    let fields = paths.iter().map(|p| corrected_month(p, tn14, tn13)).collect::<Result<Vec<_>>>()?;
    let mean = gravimetry::mean(&fields).ok_or("empty baseline")?;
    let first = calendar(span(&paths[0])?.mid());
    let last = calendar(span(paths.last().ok_or("empty baseline")?)?.mid());
    let label = if first.0 == last.0 { format!("{}", first.0) } else { format!("{}–{}", first.0, last.0) };
    eprintln!("baseline: {} months, {}-{:02} to {}-{:02}", fields.len(), first.0, first.1, last.0, last.1);
    Ok((mean, format!("{label} mean ({} months)", fields.len())))
}

/// Equivalent water height (cm) of `field − mean` on `lats × lons` (radians), through the GRACE
/// Tellus chain (Landerer & Swenson 2012): degree 60, destriping, 300 km Gaussian.
fn ewh_cm(field: &Coefficients, mean: &Coefficients, lats: &[f64], lons: &[f64]) -> Result<Vec<f64>> {
    let love = gravimetry::love_numbers(60).ok_or("degree above the Love-number table")?;
    Ok(filtered(field, mean).grid(&Functional::WaterHeight { love: &love }, field.radius, lats, lons).into_iter().map(|m| m * 100.0).collect())
}

/// `field − mean` through the Tellus filters: degree 60, destriping, 300 km Gaussian.
fn filtered(field: &Coefficients, mean: &Coefficients) -> Coefficients {
    let anomaly = gravimetry::destripe(&gravimetry::anomaly(field, mean).truncated(60));
    gravimetry::smooth(&anomaly, &gravimetry::gaussian_weights(300e3, field.radius, 60))
}

/// Minimum-norm mass transport (M2): the potential φ whose gradient is the smallest horizontal
/// mass flux that carries the observed change, ∇²φ = −σ̇ on the sphere, so φ_nm = R²σ̇_nm/(n(n+1)).
/// `rate` is the change of the filtered field per year; returns φ on the grid in Gt/yr (its
/// gradient divided by R is the flux in Gt per metre per year).
fn transport_potential(rate: &Coefficients, love: &[f64], lats: &[f64], lons: &[f64]) -> Vec<f64> {
    let r = rate.radius;
    let weights: Vec<f64> = (0..=rate.max_degree()).map(|n| if n == 0 { 0.0 } else { r * r / (n * (n + 1)) as f64 }).collect();
    // WaterHeight gives metres of water per year; × ρ_w = 1000 kg/m³ and ÷ 1e12 for Gt.
    gravimetry::smooth(rate, &weights).grid(&Functional::WaterHeight { love }, r, lats, lons).into_iter().map(|v| v * 1000.0 / 1e12).collect()
}

/// Cell centres of a regular `step`-degree grid, radians: (latitudes south to north, longitudes from −180°).
fn grid_axes(step: f64) -> (Vec<f64>, Vec<f64>) {
    let n = |span: f64| (span / step).round() as usize;
    let lats = (0..n(180.0)).map(|i| (-90.0 + step * (i as f64 + 0.5)).to_radians()).collect();
    let lons = (0..n(360.0)).map(|j| (-180.0 + step * (j as f64 + 0.5)).to_radians()).collect();
    (lats, lons)
}

fn read_technotes(tn14: &str, tn13: &str) -> Result<(Vec<SlrC20C30>, Vec<Degree1>)> {
    Ok((read_tn14(&std::fs::read_to_string(tn14)?)?, read_tn13(&std::fs::read_to_string(tn13)?)?))
}

/// The geoid of a static model against the GRS80 ellipsoid, on a 0.5° Equal Earth map.
fn geoid(model: &str, out: &str) -> Result<()> {
    let gfc = read_gfc(&std::fs::read_to_string(model)?)?;
    let name = gfc.name.clone().unwrap_or_else(|| "model".into());
    let mut field = gfc.model.base;
    // GRS80 normal field (Moritz 2000), fully normalized even zonals; degrees 0 and 1 dropped.
    for (n, c) in [(2, -4.841_667_749_85e-4), (4, 7.903_037_335_11e-7), (6, -1.687_249_611_51e-9), (8, 3.460_524_683_94e-12)] {
        field.set(n, 0, field.c(n, 0) - c, 0.0);
    }
    field.set(0, 0, 0.0, 0.0);
    field.set(1, 0, 0.0, 0.0);
    field.set(1, 1, 0.0, 0.0);
    let (lats, lons) = grid_axes(0.5);
    let n = field.grid(&Functional::GeoidHeight, field.radius, &lats, &lons);
    let (min, max) = n.iter().fold((f64::MAX, f64::MIN), |(a, b), &v| (a.min(v), b.max(v)));
    eprintln!("geoid: {min:.1} to {max:.1} m");
    let coast = coastlines();
    let title = format!("The geoid from {name}");
    let sub = format!("Height of the sea-level surface of Earth's gravity field above the GRS80 ellipsoid, from {name} to degree {}.", field.max_degree());
    let svg = map::render(&map::Grid { lats: &lats, lons: &lons, values: &n }, 100.0, &coast, &map::MapText {
        title: &title,
        subtitle: &[&sub, "Low: the Indian Ocean geoid low. High: New Guinea and the North Atlantic."],
        unit: "m",
        ends: ("below the ellipsoid", "above"),
        footnote: "Equal Earth projection. Coastlines: Natural Earth (public domain). Model: GOCO06s (CC BY 4.0, doi:10.5880/ICGEM.2019.002).",
    })?;
    write(out, &format!("geoid_{}.svg", name.to_lowercase()), &svg)
}

fn coastlines() -> Vec<Vec<(f64, f64)>> {
    include_str!("../data/coastline_110m.txt").split("\n\n").map(|seg| seg.lines().filter_map(|l| {
        let mut it = l.split_whitespace().map(str::parse::<f64>);
        Some((it.next()?.ok()?, it.next()?.ok()?))
    }).collect()).filter(|s: &Vec<_>| s.len() > 1).collect()
}

fn water_map(month: &str, baseline: &str, tn14: &str, tn13: &str, out: &str) -> Result<()> {
    let (tn14, tn13) = read_technotes(tn14, tn13)?;
    let month_path = Path::new(month);
    let (year, mon) = calendar(span(month_path)?.mid());
    let field = corrected_month(month_path, &tn14, &tn13)?;
    let (mean, baseline_label) = baseline_mean(baseline, &tn14, &tn13)?;
    let (lats, lons) = grid_axes(0.5);
    let cm = ewh_cm(&field, &mean, &lats, &lons)?;
    let nlon = lons.len();

    // Colour limit: the area-weighted 95th percentile of |value|, rounded up to 5 cm.
    let mut weighted: Vec<(f64, f64)> = cm.iter().enumerate().map(|(k, v)| (v.abs(), lats[k / nlon].cos())).collect();
    weighted.sort_by(|a, b| a.0.total_cmp(&b.0));
    let total: f64 = weighted.iter().map(|w| w.1).sum();
    let mut acc = 0.0;
    let p95 = weighted.iter().find(|w| { acc += w.1; acc >= 0.95 * total }).map_or(1.0, |w| w.0);
    let limit = (p95 / 5.0).ceil().max(1.0) * 5.0;
    let rms = (cm.iter().enumerate().map(|(k, v)| v * v * lats[k / nlon].cos()).sum::<f64>() / total).sqrt();
    let (min, max) = cm.iter().fold((f64::MAX, f64::MIN), |(a, b), &v| (a.min(v), b.max(v)));
    eprintln!("EWH: area-weighted RMS {rms:.1} cm, min {min:.1}, max {max:.1}, 95th pct |EWH| {p95:.1} → colour limit ±{limit} cm");

    let coast = coastlines();

    let month_name = MONTH_NAMES[mon as usize - 1];
    let title = format!("Water storage anomaly, {month_name} {year}");
    let sub1 = format!("Equivalent water height of the GRACE-FO {month_name} {year} field relative to the {baseline_label}; degree 60, destriped, 300 km Gaussian.");
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
    for i in (0..lats.len()).step_by(2) {
        for j in (0..nlon).step_by(2) {
            let _ = writeln!(csv, "{:.2},{:.2},{:.2}", lats[i].to_degrees(), lons[j].to_degrees(), cm[i * nlon + j]);
        }
    }
    write(out, &format!("{stem}.csv"), &csv)
}

/// Every month in `months` as a 2° grid of equivalent water height for the viewer.
fn series(months: &str, baseline: &str, tn14: &str, tn13: &str, out: &str) -> Result<()> {
    let (tn14, tn13) = read_technotes(tn14, tn13)?;
    let (mean, baseline_label) = baseline_mean(baseline, &tn14, &tn13)?;
    let (lats, lons) = grid_axes(2.0);
    let love = gravimetry::love_numbers(60).ok_or("degree above the Love-number table")?;
    let (mut bin, mut entries, mut fields) = (Vec::new(), Vec::new(), Vec::new());
    for path in gfc_files(months)? {
        let field = match corrected_month(&path, &tn14, &tn13) {
            Ok(f) => f,
            Err(e) => {
                eprintln!("skipping: {e}");
                continue;
            }
        };
        let sp = span(&path)?;
        let (year, mon) = calendar(sp.mid());
        let anomaly = filtered(&field, &mean);
        for v in anomaly.grid(&Functional::WaterHeight { love: &love }, field.radius, &lats, &lons) {
            bin.extend_from_slice(&encode_mm(v * 100.0).to_le_bytes());
        }
        entries.push((format!("{year}-{mon:02}"), sp.mid(), if path.to_string_lossy().contains("_GRFO_") { "GRACE-FO" } else { "GRACE" }));
        fields.push(anomaly);
    }

    // Transport potential from the change between neighbouring months (central where both lie
    // within 45 days, one-sided where only one does, none across longer gaps).
    let near = 45.0 / 365.25;
    let mut phi: Vec<Option<Vec<f64>>> = Vec::with_capacity(fields.len());
    for i in 0..fields.len() {
        let prev = (i > 0 && entries[i].1 - entries[i - 1].1 < near).then(|| i - 1);
        let next = (i + 1 < fields.len() && entries[i + 1].1 - entries[i].1 < near).then(|| i + 1);
        let (a, b) = match (prev, next) {
            (Some(p), Some(n)) => (p, n),
            (Some(p), None) => (p, i),
            (None, Some(n)) => (i, n),
            (None, None) => { phi.push(None); continue; }
        };
        let dt = entries[b].1 - entries[a].1;
        let mut rate = Coefficients::zeros(mean.gm, mean.radius, 60);
        for n in 0..=60 {
            for m in 0..=n {
                rate.set(n, m, (fields[b].c(n, m) - fields[a].c(n, m)) / dt, (fields[b].s(n, m) - fields[a].s(n, m)) / dt);
            }
        }
        phi.push(Some(transport_potential(&rate, &love, &lats, &lons)));
    }
    // One i16 scale for every month: twice the 99.9th percentile of |φ|, clamping the rare beyond.
    let mut mags: Vec<f64> = phi.iter().flatten().flatten().map(|v| v.abs()).collect();
    mags.sort_by(f64::total_cmp);
    let scale = 2.0 * mags.get(mags.len() * 999 / 1000).copied().unwrap_or(1.0) / 32767.0;
    let mut flow = Vec::new();
    for p in &phi {
        for k in 0..lats.len() * lons.len() {
            let v = p.as_ref().map_or(0.0, |p| p[k]);
            flow.extend_from_slice(&((v / scale).round().clamp(-32767.0, 32767.0) as i16).to_le_bytes());
        }
    }

    let months_json: Vec<String> = entries.iter().zip(&phi).map(|((label, t, mission), p)| {
        format!(r#"{{"label":"{label}","t":{t:.4},"mission":"{mission}","flow":{}}}"#, p.is_some())
    }).collect();
    let index = format!(
        r#"{{"source":"COST-G RL02.1 monthly fields (ICGEM, CC BY 4.0, doi:10.5880/COST-G.ICGEM_02_L2); C20/C30 from TN-14, degree 1 from TN-13",
"processing":"equivalent water height relative to the {baseline_label}; degree 60, destriped, 300 km Gaussian; GIA not removed",
"unit":"mm","lat0":{},"lon0":{},"step":2.0,"nlat":{},"nlon":{},
"flow":"minimum-norm transport potential φ (∇²φ = −dσ/dt) from the change between neighbouring months","flow_unit":"Gt/yr","flow_scale":{scale:e},
"months":[{}]}}
"#,
        lats[0].to_degrees(), lons[0].to_degrees(), lats.len(), lons.len(), months_json.join(",\n")
    );
    eprintln!("{} months", entries.len());
    let pack = |bytes: &[u8]| -> Result<String> {
        let mut z = flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::best());
        std::io::Write::write_all(&mut z, bytes)?;
        Ok(format!(r#"{{"encoding":"zlib+base64","data":"{}"}}"#, svg::base64(&z.finish()?)))
    };
    write(out, "tws_grids.json", pack(&bin)?)?;
    write(out, "flow_grids.json", pack(&flow)?)?;
    write(out, "coast.txt", include_str!("../data/coastline_110m.txt"))?;
    write(out, "tws_index.json", index)
}

fn encode_mm(cm: f64) -> i16 {
    (cm * 10.0).round().clamp(f64::from(i16::MIN), f64::from(i16::MAX)) as i16
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn spans_and_calendar_months() {
        let s = span(Path::new("GSM-2_2026182-2026212_GRFO_COSTG_BF01_0201.gfc")).unwrap();
        assert_eq!(s.start, julian_epoch(2026, 7, 1, 0, 0));
        assert_eq!(s.end, julian_epoch(2026, 8, 1, 0, 0));
        assert_eq!(calendar(s.mid()), (2026, 7));
        assert_eq!(calendar(julian_epoch(2004, 1, 1, 0, 0)), (2004, 1));
        assert_eq!(calendar(julian_epoch(2003, 12, 31, 23, 0)), (2003, 12));
        assert!(span(Path::new("ITSG_2026-07.gfc")).is_err());
    }

    /// The whole map chain on one known coefficient: ΔC̄20 = 1e-10 is R·ρe/(3ρw)·5/(1+k'₂)·√5·ΔC̄20
    /// of water at the pole, times the 300 km Gaussian weight of degree 2.
    #[test]
    fn map_chain_on_a_single_coefficient() {
        let (gm, r) = (3.986004415e14, 6378136.3);
        let mean = Coefficients::zeros(gm, r, 90);
        let mut month = mean.clone();
        month.set(2, 0, 1e-10, 0.0);
        let cm = ewh_cm(&month, &mean, &[90f64.to_radians()], &[0.0]).unwrap()[0];
        let w2 = gravimetry::gaussian_weights(300e3, r, 2)[2];
        let expect = r * 5517.0 / 3000.0 * 5.0 / (1.0 - 0.303) * 5f64.sqrt() * 1e-10 * w2 * 100.0;
        assert!((cm - expect).abs() < 1e-12 && (cm - 1.87).abs() < 0.01, "{cm} cm");
    }

    #[test]
    fn transport_potential_inverts_the_laplacian() {
        // ∇²φ = −σ̇: check with a finite-difference spherical Laplacian at 30°N, 20°E.
        let mut rate = Coefficients::zeros(3.986004415e14, 6_378_136.3, 12);
        for (n, m, c, sn) in [(2, 0, 3e-10, 0.0), (3, 1, -2e-10, 1e-10), (5, 4, 1e-10, 2e-10), (11, 7, -5e-11, 4e-11)] {
            rate.set(n, m, c, sn);
        }
        let love = gravimetry::love_numbers(12).unwrap();
        let (lat, lon, h) = (30f64.to_radians(), 20f64.to_radians(), 0.2f64.to_radians());
        let phi = |dlat: f64, dlon: f64| transport_potential(&rate, &love, &[lat + dlat], &[lon + dlon])[0];
        let r = rate.radius;
        let d_lat = (lat + h / 2.0).cos() * (phi(h, 0.0) - phi(0.0, 0.0)) - (lat - h / 2.0).cos() * (phi(0.0, 0.0) - phi(-h, 0.0));
        let laplacian = (d_lat / (h * h * lat.cos()) + (phi(0.0, h) - 2.0 * phi(0.0, 0.0) + phi(0.0, -h)) / (h * h * lat.cos().powi(2))) / (r * r);
        let ewh_rate = rate.grid(&Functional::WaterHeight { love: &love }, r, &[lat], &[lon])[0]; // m/yr
        let sigma_dot = ewh_rate * 1000.0 / 1e12; // Gt m⁻² yr⁻¹
        assert!((laplacian + sigma_dot).abs() < 1e-4 * sigma_dot.abs(), "{laplacian:e} vs {:e}", -sigma_dot);
    }

    #[test]
    fn grid_axes_and_encoding() {
        let (lats, lons) = grid_axes(2.0);
        assert_eq!((lats.len(), lons.len()), (90, 180));
        assert!((lats[0].to_degrees() + 89.0).abs() < 1e-12 && (lons[179].to_degrees() - 179.0).abs() < 1e-12);
        assert_eq!((encode_mm(-87.64), encode_mm(0.04), encode_mm(1e6)), (-876, 0, i16::MAX));
    }
}
