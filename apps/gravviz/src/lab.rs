//! Shared pieces for the hypothesis tests: the monthly COST-G record as water-height grids, a land
//! mask, global land storage, climate-index readers and small statistics helpers.

use std::path::Path;

use gravity_field::{read_tn13, read_tn14};

use crate::{Result, Span, baseline_mean, calendar, corrected_month, ewh_cm, gfc_files, grid_axes, span};

/// One month of the record: its span and its 2° water-height grid (cm, GRACE Tellus chain).
pub struct Month {
    pub span: Span,
    pub mission: &'static str,
    pub cm: Vec<f64>,
}

/// Every month in `months_dir`, relative to the mean of `baseline_dir`, on the 2° grid.
pub fn load_months(months_dir: &str, baseline_dir: &str, tn14: &str, tn13: &str) -> Result<(Vec<f64>, Vec<f64>, Vec<Month>)> {
    let tn14 = read_tn14(&std::fs::read_to_string(tn14)?)?;
    let tn13 = read_tn13(&std::fs::read_to_string(tn13)?)?;
    let (mean, _) = baseline_mean(baseline_dir, &tn14, &tn13)?;
    let (lats, lons) = grid_axes(2.0);
    let mut out = Vec::new();
    for path in gfc_files(months_dir)? {
        let field = corrected_month(&path, &tn14, &tn13)?;
        let mission = if path.to_string_lossy().contains("_GRFO_") { "GRACE-FO" } else { "GRACE" };
        out.push(Month { span: span(&path)?, mission, cm: ewh_cm(&field, &mean, &lats, &lons)? });
    }
    out.sort_by(|a, b| a.span.mid().total_cmp(&b.span.mid()));
    Ok((lats, lons, out))
}

/// Land rings (lon, lat degrees) from the bundled Natural Earth 1:110m polygons.
pub fn land_rings() -> Vec<Vec<(f64, f64)>> {
    include_str!("../data/land_110m.txt")
        .split("\n\n")
        .map(|ring| ring.lines().filter_map(|l| {
            let mut it = l.split_whitespace().map(str::parse::<f64>);
            Some((it.next()?.ok()?, it.next()?.ok()?))
        }).collect::<Vec<_>>())
        .filter(|r| r.len() > 2)
        .collect()
}

/// Even–odd point-in-polygon over all rings (holes included).
pub fn on_land(rings: &[Vec<(f64, f64)>], lon: f64, lat: f64) -> bool {
    let mut inside = false;
    for ring in rings {
        let mut j = ring.len() - 1;
        for i in 0..ring.len() {
            let ((xi, yi), (xj, yj)) = (ring[i], ring[j]);
            if (yi > lat) != (yj > lat) && lon < (xj - xi) * (lat - yi) / (yj - yi) + xi {
                inside = !inside;
            }
            j = i;
        }
    }
    inside
}

/// Area (km²) of land each 2° cell contributes to "land outside the ice sheets": the land fraction
/// from 4×4 sub-samples, with Antarctica (south of 60°S) and a Greenland box (north of 58°N, 80°W
/// to 10°W, which also drops Baffin and Ellesmere Islands and Iceland, where Greenland's signal
/// leaks) set to zero.
pub fn land_area_km2(lats: &[f64], lons: &[f64]) -> Vec<f64> {
    let rings = land_rings();
    let r = 6371.0;
    let mut out = Vec::with_capacity(lats.len() * lons.len());
    for &lat in lats {
        let (la, step) = (lat.to_degrees(), 2.0_f64);
        let cell = r * r * step.to_radians() * ((la + 1.0).to_radians().sin() - (la - 1.0).to_radians().sin());
        for &lon in lons {
            let lo = lon.to_degrees();
            let excluded = la < -60.0 || (la > 58.0 && (-80.0..=-10.0).contains(&lo));
            let hits = (0..16).filter(|k| {
                let (dy, dx) = ((k / 4) as f64 * 0.5 - 0.75, (k % 4) as f64 * 0.5 - 0.75);
                on_land(&rings, lo + dx, la + dy)
            }).count();
            out.push(if excluded { 0.0 } else { cell * hits as f64 / 16.0 });
        }
    }
    out
}

/// Σ EWH · area over land, km³ (cm × km² → km³ is ×1e-5).
pub fn land_volume_km3(cm: &[f64], area_km2: &[f64]) -> f64 {
    cm.iter().zip(area_km2).map(|(c, a)| c * a).sum::<f64>() * 1e-5
}

/// A NOAA PSL climate-index file: a "first last" year line, then `year v1 … v12` rows, then a
/// missing-value line (−99.99, −9.90, −999, … depending on the index). Returns ((year, month),
/// value) for every non-missing month.
pub fn read_psl_index(path: &Path) -> Result<Vec<((i32, u32), f64)>> {
    let text = std::fs::read_to_string(path)?;
    let mut lines = text.lines();
    let first = lines.next().ok_or("empty index file")?;
    let years: Vec<i32> = first.split_whitespace().filter_map(|s| s.parse().ok()).collect();
    let [y0, y1] = years[..] else { return Err("bad index header".into()) };
    let mut table = Vec::new();
    let mut missing = None;
    for line in lines {
        let f: Vec<f64> = line.split_whitespace().filter_map(|s| s.parse().ok()).collect();
        if f.len() == 13 && f[0].fract() == 0.0 && (y0..=y1).contains(&(f[0] as i32)) {
            table.push(f);
        } else if let ([v], false) = (f.as_slice(), table.is_empty()) {
            let v = *v;
            missing = Some(v);
            break;
        }
    }
    let missing = missing.ok_or("index file has no missing-value line")?;
    Ok(table.iter().flat_map(|f| (1..13).filter(|&m| (f[m] - missing).abs() > 1e-6).map(|m| ((f[0] as i32, m as u32), f[m]))).collect())
}

/// Ordinary least squares by the normal equations (small, well-scaled designs only).
/// Returns the coefficients and the residuals.
pub fn ols(x: &[Vec<f64>], y: &[f64]) -> Result<(Vec<f64>, Vec<f64>)> {
    let p = x[0].len();
    let mut a = vec![0.0; p * (p + 1)];
    for (row, &yi) in x.iter().zip(y) {
        for i in 0..p {
            for j in 0..p {
                a[i * (p + 1) + j] += row[i] * row[j];
            }
            a[i * (p + 1) + p] += row[i] * yi;
        }
    }
    // Gauss–Jordan with partial pivoting on the augmented matrix.
    for c in 0..p {
        let piv = (c..p).max_by(|&i, &j| a[i * (p + 1) + c].abs().total_cmp(&a[j * (p + 1) + c].abs())).unwrap_or(c);
        if a[piv * (p + 1) + c].abs() < 1e-12 {
            return Err("singular design matrix".into());
        }
        for k in 0..=p {
            a.swap(c * (p + 1) + k, piv * (p + 1) + k);
        }
        for i in (0..p).filter(|&i| i != c) {
            let f = a[i * (p + 1) + c] / a[c * (p + 1) + c];
            for k in 0..=p {
                a[i * (p + 1) + k] -= f * a[c * (p + 1) + k];
            }
        }
    }
    let beta: Vec<f64> = (0..p).map(|i| a[i * (p + 1) + p] / a[i * (p + 1) + i]).collect();
    let resid = x.iter().zip(y).map(|(row, yi)| yi - row.iter().zip(&beta).map(|(a, b)| a * b).sum::<f64>()).collect();
    Ok((beta, resid))
}

/// A small deterministic random generator (SplitMix64) for bootstraps and ensembles.
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Self {
        Rng(seed)
    }
    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    pub fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }
}

/// Moving-block bootstrap resample of `x` with blocks of `block` consecutive values.
pub fn block_resample(x: &[f64], block: usize, rng: &mut Rng) -> Vec<f64> {
    let mut out = Vec::with_capacity(x.len());
    while out.len() < x.len() {
        let start = rng.below(x.len() - block + 1);
        out.extend_from_slice(&x[start..start + block]);
    }
    out.truncate(x.len());
    out
}

/// The `q`-quantile (0–1) of `v` (sorted copy, linear interpolation).
pub fn quantile(v: &[f64], q: f64) -> f64 {
    let mut s = v.to_vec();
    s.sort_by(f64::total_cmp);
    let pos = q * (s.len() - 1) as f64;
    let (i, f) = (pos.floor() as usize, pos.fract());
    if i + 1 < s.len() { s[i] * (1.0 - f) + s[i + 1] * f } else { s[i] }
}

/// Calendar (year, month) of a Julian epoch, re-exported for the analyses.
pub fn ym(t: f64) -> (i32, u32) {
    calendar(t)
}

/// Months between two calendar months.
pub fn months_between(a: (i32, u32), b: (i32, u32)) -> i32 {
    (b.0 - a.0) * 12 + b.1 as i32 - a.1 as i32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn land_mask_and_area_are_sane() {
        let rings = land_rings();
        assert!(on_land(&rings, 15.0, 50.0)); // central Europe
        assert!(on_land(&rings, -60.0, -10.0)); // Amazon
        assert!(!on_land(&rings, -150.0, 0.0)); // central Pacific
        assert!(!on_land(&rings, 60.0, -30.0)); // southern Indian Ocean
        let (lats, lons) = grid_axes(2.0);
        let total: f64 = land_area_km2(&lats, &lons).iter().sum();
        // Land is ~149 million km²; minus Antarctica (~14 M) and Greenland box (~4 M).
        assert!((1.15e8..1.40e8).contains(&total), "{total:e}");
    }

    #[test]
    fn ols_recovers_a_line() {
        let x: Vec<Vec<f64>> = (0..20).map(|i| vec![1.0, i as f64]).collect();
        let y: Vec<f64> = (0..20).map(|i| 3.0 - 0.5 * i as f64).collect();
        let (b, r) = ols(&x, &y).unwrap();
        assert!((b[0] - 3.0).abs() < 1e-10 && (b[1] + 0.5).abs() < 1e-10 && r.iter().all(|e| e.abs() < 1e-9));
    }

    #[test]
    fn psl_index_missing_values() {
        let path = std::env::temp_dir().join(format!("gravviz_psl_{}.data", std::process::id()));
        std::fs::write(&path, " 2024 2025\n 2024 1 2 3 4 5 6 7 8 9 10 11 12\n 2025 -1 -2 -9.90 -9.90 -9.90 -9.90 -9.90 -9.90 -9.90 -9.90 -9.90 -9.90\n -9.90\n PDO\n").unwrap();
        let rows = read_psl_index(&path).unwrap();
        std::fs::remove_file(&path).unwrap();
        assert_eq!(rows.len(), 14);
        assert_eq!(rows[13], ((2025, 2), -2.0));
    }

    #[test]
    fn quantiles_and_months() {
        assert_eq!(quantile(&[3.0, 1.0, 2.0], 0.5), 2.0);
        assert_eq!(months_between((2025, 11), (2026, 2)), 3);
    }
}
