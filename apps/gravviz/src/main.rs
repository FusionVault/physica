//! Figures and viewer data from the quantum-gravimeter crates.
//!
//! ```text
//! gravviz noise OUT_DIR
//! gravviz water-map MONTH.gfc BASELINE_DIR TN-14.txt TN-13.txt OUT_DIR
//! gravviz series MONTHS_DIR BASELINE_DIR TN-14.txt TN-13.txt OUT_DIR
//! gravviz triangle MONTH.gfc OUT_DIR
//! gravviz gradiometry STATIC_MODEL.gfc OUT_DIR
//! gravviz sst OUT_DIR [REFERENCE_MONTH.gfc]
//! gravviz w1 MONTHS_DIR BASELINE_DIR TN-14.txt TN-13.txt NINO34.data OUT_DIR
//! gravviz e1 MONTHS_DIR BASELINE_DIR TN-14.txt TN-13.txt OUT_DIR
//! gravviz screen W1_LAND_STORAGE.csv INDEX_DIR OUT_DIR
//! ```
//!
//! `noise` writes `noise_spectra.svg` and `.csv` (Zingerle et al. 2024, Fig. 1). `water-map`
//! writes `water_<yyyy-mm>.svg` and a 1° `.csv`: the month's equivalent water height minus the
//! mean of every `.gfc` in BASELINE_DIR, all from ICGEM GSM-2 files (`GSM-2_YYYYDDD-YYYYDDD_…`).
//! `series` does the same for every month in MONTHS_DIR on a 2° grid and writes `tws_grids.json`
//! (zlib-compressed little-endian i16 millimetres, `[month][lat][lon]`, base64 in JSON because
//! artifact hosting serves no raw binary), `tws_index.json` and the coastlines (`coast.txt`) for the
//! browser viewer in `viewer/`. `triangle` writes `error_triangle_<yyyy-mm>.svg`: the month's calibrated
//! coefficient errors as water height. `gradiometry` writes `gradiometry_errors.svg` and `.csv`:
//! gradiometer mission errors against the static model's signal (e.g. GOCO06s). `sst` runs the
//! closed-loop simulation of four GRACE-type mission classes and writes `sst_closed_loop.svg` and
//! `.csv`, with a real month's calibrated errors as reference if given. `w1` runs the pre-registered
//! El Niño test on global land storage (NOAA PSL Niño 3.4 anomaly file) and writes `w1_report.txt`,
//! `w1_land_storage.svg` and `.csv`, and `w1_forecast.csv`. `e1` looks for the 2025 Kamchatka
//! earthquake in GRACE-FO gravity against a random-location null and writes `e1_kamchatka.svg`,
//! `e1_report.txt` and `e1_step.csv`. `screen` reruns the W1 level model with every NOAA PSL index in
//! INDEX_DIR, QBO and solar flux as negative controls, and writes `screen_indices.svg` and `.txt`.

mod e1;
mod gradiometry;
mod lab;
mod map;
mod noise;
mod screen;
mod sst;
mod svg;
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
            let w1::Outcome { report, svg, csv, forecast_csv } = w1::run(&fields, &lats, &lons, Path::new(nino))?;
            print!("{report}");
            write(out, "w1_report.txt", &report)?;
            write(out, "w1_land_storage.svg", &svg)?;
            write(out, "w1_land_storage.csv", &csv)?;
            write(out, "w1_forecast.csv", &forecast_csv)
        }
        ["e1", months, baseline, tn14, tn13, out] => {
            let e1::Outcome { report, svg, csv } = e1::run(months, baseline, tn14, tn13)?;
            print!("{report}");
            write(out, "e1_report.txt", &report)?;
            write(out, "e1_kamchatka.svg", &svg)?;
            write(out, "e1_step.csv", &csv)
        }
        ["screen", csv, indices, out] => {
            let screen::Outcome { report, svg } = screen::run(Path::new(csv), Path::new(indices))?;
            print!("{report}");
            write(out, "screen_indices.txt", &report)?;
            write(out, "screen_indices.svg", &svg)
        }
        _ => Err("usage: gravviz noise OUT_DIR | water-map MONTH.gfc BASELINE_DIR TN-14.txt TN-13.txt OUT_DIR | \
            series MONTHS_DIR BASELINE_DIR TN-14.txt TN-13.txt OUT_DIR | triangle MONTH.gfc OUT_DIR | \
            gradiometry STATIC_MODEL.gfc OUT_DIR | sst OUT_DIR [REFERENCE_MONTH.gfc] | \
            w1 MONTHS_DIR BASELINE_DIR TN-14.txt TN-13.txt NINO34.data OUT_DIR | \
            e1 MONTHS_DIR BASELINE_DIR TN-14.txt TN-13.txt OUT_DIR | screen W1_LAND_STORAGE.csv INDEX_DIR OUT_DIR".into()),
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
    let nmax = 60;
    let anomaly = gravimetry::destripe(&gravimetry::anomaly(field, mean).truncated(nmax));
    let anomaly = gravimetry::smooth(&anomaly, &gravimetry::gaussian_weights(300e3, field.radius, nmax));
    let love = gravimetry::love_numbers(nmax).ok_or("degree above the Love-number table")?;
    Ok(anomaly.grid(&Functional::WaterHeight { love: &love }, field.radius, lats, lons).into_iter().map(|m| m * 100.0).collect())
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

    let coast_text = include_str!("../data/coastline_110m.txt");
    let coast: Vec<Vec<(f64, f64)>> = coast_text.split("\n\n").map(|seg| seg.lines().filter_map(|l| {
        let mut it = l.split_whitespace().map(str::parse::<f64>);
        Some((it.next()?.ok()?, it.next()?.ok()?))
    }).collect()).filter(|s: &Vec<_>| s.len() > 1).collect();

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
    let (mut bin, mut entries) = (Vec::new(), Vec::new());
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
        for v in ewh_cm(&field, &mean, &lats, &lons)? {
            bin.extend_from_slice(&encode_mm(v).to_le_bytes());
        }
        entries.push(format!(r#"{{"label":"{year}-{mon:02}","t":{:.4},"mission":"{}"}}"#, sp.mid(),
            if path.to_string_lossy().contains("_GRFO_") { "GRACE-FO" } else { "GRACE" }));
    }
    let index = format!(
        r#"{{"source":"COST-G RL02.1 monthly fields (ICGEM, CC BY 4.0, doi:10.5880/COST-G.ICGEM_02_L2); C20/C30 from TN-14, degree 1 from TN-13",
"processing":"equivalent water height relative to the {baseline_label}; degree 60, destriped, 300 km Gaussian; GIA not removed",
"unit":"mm","lat0":{},"lon0":{},"step":2.0,"nlat":{},"nlon":{},
"months":[{}]}}
"#,
        lats[0].to_degrees(), lons[0].to_degrees(), lats.len(), lons.len(), entries.join(",\n")
    );
    eprintln!("{} months", entries.len());
    let mut z = flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::best());
    std::io::Write::write_all(&mut z, &bin)?;
    write(out, "tws_grids.json", format!(r#"{{"encoding":"zlib+base64","data":"{}"}}"#, svg::base64(&z.finish()?)))?;
    write(out, "coast.txt", include_str!("../data/coastline_110m.txt"))?;
    write(out, "tws_index.json", index)
}

/// Equivalent water height in cm → i16 millimetres, saturating (±32.7 m covers any real month).
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
    fn grid_axes_and_encoding() {
        let (lats, lons) = grid_axes(2.0);
        assert_eq!((lats.len(), lons.len()), (90, 180));
        assert!((lats[0].to_degrees() + 89.0).abs() < 1e-12 && (lons[179].to_degrees() - 179.0).abs() < 1e-12);
        assert_eq!((encode_mm(-87.64), encode_mm(0.04), encode_mm(1e6)), (-876, 0, i16::MAX));
    }
}
