//! Q1 negative controls and a climate-index screen: the W1 level model with each index in place
//! of Niño 3.4.
//!
//! For each index X, standardised over its record from 1950:
//! S′ = a + b·(t − 2015) + γ⁺·X̄⁺ + γ⁻·X̄⁻ + e, with X̄ the mean at lags 2–6 months, fitted on the
//! months X covers. Significance comes from 1,000 circular shifts of X by at least five years,
//! which keep its autocorrelation but break its alignment with S′. Benjamini–Hochberg q-values run
//! over every γ⁺, γ⁻ and γ⁺ − γ⁻. QBO and solar flux are the negative controls (Q1) and should be
//! null. The lags were chosen for Niño 3.4, so for the other indices this is a screen, not a test.

use std::collections::HashMap;
use std::fmt::Write;
use std::path::Path;

use crate::Result;
use crate::lab::{Rng, bh, ols, read_psl_index, ym};
use crate::svg::{esc, open};

type Cal = (i32, u32);

/// (file, label, negative control)
const INDICES: [(&str, &str, bool); 9] = [
    ("nina34.anom.data", "Niño 3.4", false),
    ("dmi.had.long.data", "Indian Ocean Dipole", false),
    ("pdo.data", "Pacific Decadal Osc.", false),
    ("amon.us.data", "Atlantic Multidecadal Osc.", false),
    ("nao.data", "North Atlantic Osc.", false),
    ("ao.data", "Arctic Osc.", false),
    ("aao.data", "Southern Annular Mode", false),
    ("qbo.data", "QBO (control)", true),
    ("solar.data", "Solar flux (control)", true),
];
const SHIFTS: usize = 1000;

struct Row {
    label: &'static str,
    control: bool,
    months: usize,
    /// (γ⁺, γ⁻, γ⁺ − γ⁻), km³ per standard deviation of the index.
    gamma: [f64; 3],
    p: [f64; 3],
    z: [f64; 3],
}

pub struct Outcome {
    pub report: String,
    pub svg: String,
    pub csv: String,
}

/// S′ from the W1 CSV: (epoch, km³).
fn read_s_prime(path: &Path) -> Result<Vec<(f64, f64)>> {
    let text = std::fs::read_to_string(path)?;
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().ok_or("empty CSV")?.split(',').collect();
    let col = |name: &str| header.iter().position(|h| *h == name).ok_or(format!("no {name} column"));
    let (ct, cs) = (col("epoch_year")?, col("s_prime_km3")?);
    lines.map(|l| {
        let f: Vec<&str> = l.split(',').collect();
        Ok((f[ct].parse()?, f[cs].parse()?))
    }).collect()
}

fn shift(cm: Cal, k: i32) -> Cal {
    let n = cm.0 * 12 + cm.1 as i32 - 1 + k;
    (n.div_euclid(12), n.rem_euclid(12) as u32 + 1)
}

pub fn run(w1_csv: &Path, index_dir: &Path) -> Result<Outcome> {
    let sp = read_s_prime(w1_csv)?;
    let mut rng = Rng::new(1950);
    let mut rows = Vec::new();
    for (file, label, control) in INDICES {
        // Standardised index as a contiguous monthly vector from its first month in or after 1950;
        // gaps inside the record become 0 (the mean).
        let raw: HashMap<Cal, f64> = read_psl_index(&index_dir.join(file))?.into_iter().filter(|((y, _), _)| *y >= 1950).collect();
        let first = *raw.keys().min().ok_or("empty index")?;
        let last = *raw.keys().max().ok_or("empty index")?;
        let n = (last.0 * 12 + last.1 as i32 - first.0 * 12 - first.1 as i32 + 1) as usize;
        let (mean, sd) = {
            let v: Vec<f64> = raw.values().copied().collect();
            let m = v.iter().sum::<f64>() / v.len() as f64;
            (m, (v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (v.len() - 1) as f64).sqrt())
        };
        let series: Vec<f64> = (0..n as i32).map(|k| raw.get(&shift(first, k)).map_or(0.0, |v| (v - mean) / sd)).collect();
        let pos = |cm: Cal| cm.0 * 12 + cm.1 as i32 - first.0 * 12 - first.1 as i32;

        // Months whose lags 2–6 are all observed, and their lag positions.
        let sample: Vec<(f64, f64, Vec<i32>)> = sp.iter().filter_map(|&(t, s)| {
            let cm = ym(t);
            let lags: Vec<i32> = (2..=6).map(|k| shift(cm, -k)).map(|c| raw.contains_key(&c).then(|| pos(c))).collect::<Option<_>>()?;
            Some((t, s, lags))
        }).collect();
        let y: Vec<f64> = sample.iter().map(|r| r.1).collect();
        let fit = |offset: usize| -> Result<[f64; 3]> {
            let x: Vec<Vec<f64>> = sample.iter().map(|(t, _, lags)| {
                let xbar = lags.iter().map(|&p| series[(p as usize + offset) % n]).sum::<f64>() / lags.len() as f64;
                vec![1.0, t - 2015.0, xbar.max(0.0), xbar.min(0.0)]
            }).collect();
            let (b, _) = ols(&x, &y)?;
            Ok([b[2], b[3], b[2] - b[3]])
        };
        let gamma = fit(0)?;
        let null: Vec<[f64; 3]> = (0..SHIFTS).map(|_| fit(60 + rng.below(n - 120))).collect::<Result<_>>()?;
        let mut p = [0.0; 3];
        let mut z = [0.0; 3];
        for k in 0..3 {
            let v: Vec<f64> = null.iter().map(|g| g[k]).collect();
            p[k] = (1 + v.iter().filter(|g| g.abs() >= gamma[k].abs()).count()) as f64 / (SHIFTS + 1) as f64;
            let m = v.iter().sum::<f64>() / v.len() as f64;
            let s = (v.iter().map(|g| (g - m).powi(2)).sum::<f64>() / (v.len() - 1) as f64).sqrt();
            z[k] = (gamma[k] - m) / s;
        }
        rows.push(Row { label, control, months: sample.len(), gamma, p, z });
    }

    // Benjamini–Hochberg over all tests.
    let all: Vec<f64> = rows.iter().flat_map(|r| r.p).collect();
    let q = bh(&all);

    let mut rep = String::from("index                       months   γ⁺ km³/σ  p      q      γ⁻ km³/σ  p      q      γ⁺−γ⁻   p      q\n");
    for (i, r) in rows.iter().enumerate() {
        let _ = write!(rep, "{:<27} {:>6}", r.label, r.months);
        for k in 0..3 {
            let _ = write!(rep, "  {:>8.0}  {:.3}  {:.3}", r.gamma[k], r.p[k], q[3 * i + k]);
        }
        rep.push('\n');
    }
    let _ = writeln!(rep, "p: two-sided, {SHIFTS} circular shifts of at least 5 years; q: Benjamini–Hochberg over all {} tests", all.len());
    let q_ref = &q;
    let flagged: Vec<String> = rows.iter().enumerate().flat_map(|(i, r)| (0..3).filter(move |&k| q_ref[3 * i + k] < 0.1).map(move |k| format!("{} {}", r.label, ["γ⁺", "γ⁻", "γ⁺−γ⁻"][k]))).collect();
    let _ = writeln!(rep, "q < 0.1: {}", if flagged.is_empty() { "none".into() } else { flagged.join("; ") });
    let controls_null = rows.iter().enumerate().filter(|(_, r)| r.control).all(|(i, _)| (0..3).all(|k| q[3 * i + k] >= 0.1));
    let _ = writeln!(rep, "negative controls (QBO, solar flux): {}", if controls_null { "null, as required" } else { "NOT null: the pipeline leaks" });
    let svg = figure(&rows, &q);
    let mut csv = String::from("index,control,months,parameter,gamma_km3_per_sd,z,p,q\n");
    for (i, r) in rows.iter().enumerate() {
        for (k, name) in ["positive", "negative", "asymmetry"].iter().enumerate() {
            let _ = writeln!(csv, "{},{},{},{name},{:.1},{:.3},{:.4},{:.4}", r.label, r.control, r.months, r.gamma[k], r.z[k], r.p[k], q[3 * i + k]);
        }
    }
    Ok(Outcome { report: rep, svg, csv })
}

fn figure(rows: &[Row], q: &[f64]) -> String {
    let (left, right, top, step) = (250.0, 900.0, 110.0, 38.0);
    let h = top + step * rows.len() as f64 + 90.0;
    let zmax = rows.iter().flat_map(|r| r.z).fold(4.0, |a: f64, b| a.max(b.abs())).ceil();
    let x = |z: f64| left + (right - left) * (z + zmax) / (2.0 * zmax);
    let mut s = open(940.0, h, "Land water storage against climate indices");
    let _ = write!(s, r#"<text x="20" y="32" class="title">Which climate patterns move global land water?</text>
<text x="20" y="54" class="subtitle">W1 level model with each index in place of Niño 3.4: how each phase moves land storage, as z against 1,000 time-shifted copies of the index.</text>
<text x="20" y="72" class="subtitle">Blue: positive phase. Orange: negative phase. Right: that phase adds water to land. Filled: q &lt; 0.1 over all 27 tests. Grey band: |z| &lt; 2.</text>
"#);
    let _ = writeln!(s, r#"<rect x="{:.1}" y="{top}" width="{:.1}" height="{:.1}" fill="var(--grid)" opacity="0.6"/>"#, x(-2.0), x(2.0) - x(-2.0), step * rows.len() as f64);
    let _ = writeln!(s, r#"<line x1="{0:.1}" y1="{top}" x2="{0:.1}" y2="{1:.1}" class="axis"/>"#, x(0.0), top + step * rows.len() as f64);
    let mut z = -zmax;
    while z <= zmax {
        let _ = writeln!(s, r#"<text x="{:.1}" y="{:.1}" class="tick" text-anchor="middle">{z:+.0}</text>"#, x(z), top + step * rows.len() as f64 + 18.0);
        z += 2.0;
    }
    let base = top + step * rows.len() as f64 + 38.0;
    let _ = writeln!(s, r#"<text x="{:.1}" y="{base:.1}" class="tick">← land loses water</text><text x="{:.1}" y="{base:.1}" class="tick" text-anchor="middle">z</text><text x="{:.1}" y="{base:.1}" class="tick" text-anchor="end">land gains water →</text>"#, x(-zmax), x(0.0), x(zmax));
    for (i, r) in rows.iter().enumerate() {
        let yy = top + step * (i as f64 + 0.5);
        let _ = writeln!(s, r#"<text x="{:.1}" y="{:.1}" text-anchor="end"{}>{}</text>"#, left - 14.0, yy + 4.0, if r.control { r#" class="tick""# } else { "" }, esc(r.label));
        // A phase's effect on storage per σ: γ⁺ for the positive phase, −γ⁻ for the negative one.
        for (k, sign, colour, phase) in [(0, 1.0, "--s1", "positive"), (1, -1.0, "--s2", "negative")] {
            let filled = q[3 * i + k] < 0.1;
            let dy = if k == 0 { -5.0 } else { 5.0 };
            let _ = writeln!(s, r#"<circle cx="{:.1}" cy="{:.1}" r="6" fill="{}" stroke="var({colour})" stroke-width="2"><title>{}, {phase} phase: {:+.0} km³ per σ, p = {:.3}, q = {:.3}</title></circle>"#,
                x((sign * r.z[k]).clamp(-zmax, zmax)), yy + dy, if filled { format!("var({colour})") } else { "var(--surface)".into() }, esc(r.label), sign * r.gamma[k], r.p[k], q[3 * i + k]);
        }
    }
    s.push_str("</svg>\n");
    s
}
