//! M1 and E2: a blind global scan of the COST-G record for sudden mass changes, scored against
//! the earthquake catalogue.
//!
//! For every candidate month τ (a step at the start of that month), each coefficient of degree
//! 2–60 gets a step estimate from the months within a year of τ: a least-squares fit of an
//! intercept, a local trend, annual and semi-annual terms and the step. Its variance under "no
//! step" comes from the M3 Kalman model of that degree: the natural signal (an Ornstein–Uhlenbeck
//! process) pushed through the step estimator, plus the files' errors scaled by Q2's κ. A
//! generalised likelihood-ratio test then looks for a
//! uniform spherical cap of mass (400 or 800 km radius, either sign, with elastic loading) centred
//! on every 2° cell:
//!
//! `Z = Σ g·s/v / √(Σ g²/v)`, with `g_nm = c_n · Ȳ_nm(centre)`,
//!
//! so the numerator at every centre is one spherical-harmonic synthesis, and the denominator one
//! synthesis of squared functions. Real errors are correlated (stripes), so each cell's Z is
//! divided by its robust spread across all months. The null is the month's global maximum when
//! the months' residuals are shifted circularly in time, which keeps each month's real noise
//! pattern and the record's temporal correlation but moves any step away from the epoch tested;
//! the threshold is the 99th percentile of a Gumbel law fitted to its quartiles (the shifted real
//! earthquakes contaminate its top tail). The catalogue (USGS, Mw ≥ 8) is used only to score.

use std::fmt::Write;
use std::path::Path;

use geopotential::{Coefficients, Legendre, RHO_EARTH};
use geopotential::read_gfc;
use gravimetry::technotes::{read_tn13, read_tn14};

use crate::shared::{Rng, land_rings, ols, quantile};
use crate::svg::{esc, open};
use crate::{Result, Span, baseline_mean, calendar, corrected_month, gfc_files, grid_axes, span, write};

const NMIN: usize = 2;
const NMAX: usize = 60;
const RADII_KM: [f64; 2] = [400.0, 800.0];
const WINDOW_YEARS: f64 = 1.0;
const SHUFFLES: usize = 8; // circular shifts for the null
const NULL_EPOCHS: usize = 40;

/// Coefficients as flat (C, S) arrays in triangular order.
#[derive(Clone)]
struct Flat {
    c: Vec<f64>,
    s: Vec<f64>,
}

fn idx(n: usize, m: usize) -> usize {
    n * (n + 1) / 2 + m
}

fn flatten(f: &Coefficients) -> Flat {
    let size = idx(NMAX, NMAX) + 1;
    let (mut c, mut s) = (vec![0.0; size], vec![0.0; size]);
    for n in 0..=NMAX {
        for m in 0..=n {
            c[idx(n, m)] = f.c(n, m);
            s[idx(n, m)] = f.s(n, m);
        }
    }
    Flat { c, s }
}

/// Legendre-polynomial coefficients `c_n` (Stokes per kg/m² of cap density) of a uniform cap of
/// angular radius `psi`, loading an elastic Earth: σ_n = ½[P_{n−1}(cos ψ) − P_{n+1}(cos ψ)] and
/// `c_n = 3 (1 + k'_n) σ_n / (R ρ_e (2n + 1)²)`.
fn cap_template(psi: f64, radius: f64, love: &[f64]) -> Vec<f64> {
    let x = psi.cos();
    let mut p = vec![1.0, x];
    for n in 1..=NMAX {
        p.push(((2 * n + 1) as f64 * x * p[n] - n as f64 * p[n - 1]) / (n + 1) as f64);
    }
    (0..=NMAX).map(|n| {
        if n < NMIN { return 0.0; }
        let sigma = 0.5 * (p[n - 1] - p[n + 1]);
        3.0 * (1.0 + love[n]) * sigma / (radius * RHO_EARTH * ((2 * n + 1) as f64).powi(2))
    }).collect()
}

/// `Σ_nm (a_nm cos mλ + b_nm sin mλ) P̄_nm(sin φ)` on the grid, and, if `squared`, the same with
/// P̄² cos² and P̄² sin² (for the GLRT denominator).
fn synthesize(a: &[f64], b: &[f64], squared: bool, lats: &[f64], lons: &[f64], leg: &mut Legendre) -> Vec<f64> {
    let mut out = Vec::with_capacity(lats.len() * lons.len());
    let trig: Vec<Vec<(f64, f64)>> = lons.iter().map(|&l| (0..=NMAX).map(|m| (m as f64 * l).sin_cos()).collect()).collect();
    for &lat in lats {
        leg.compute(lat);
        let (mut am, mut bm) = (vec![0.0; NMAX + 1], vec![0.0; NMAX + 1]);
        for n in NMIN..=NMAX {
            for m in 0..=n {
                let p = leg.p(n, m);
                let w = if squared { p * p } else { p };
                am[m] += a[idx(n, m)] * w;
                bm[m] += b[idx(n, m)] * w;
            }
        }
        for t in &trig {
            out.push((0..=NMAX).map(|m| {
                let (sm, cm) = t[m];
                if squared { am[m] * cm * cm + bm[m] * sm * sm } else { am[m] * cm + bm[m] * sm }
            }).sum());
        }
    }
    out
}

/// A catalogue earthquake.
struct Quake {
    t: f64,
    lat: f64,
    lon: f64,
    depth: f64,
    mag: f64,
    name: String,
}

fn read_catalogue(path: &Path) -> Result<Vec<Quake>> {
    let text = std::fs::read_to_string(path)?;
    text.lines().skip(1).map(|l| {
        let f: Vec<&str> = l.splitn(15, ',').collect();
        let date = f.first().ok_or("bad catalogue line")?;
        let (y, mo, d) = (date[0..4].parse::<i32>()?, date[5..7].parse::<u32>()?, date[8..10].parse::<u32>()?);
        let name = l.rsplit('"').nth(1).unwrap_or("").to_string();
        Ok(Quake { t: geopotential::julian_epoch(y, mo, d, 12, 0), lat: f[1].parse()?, lon: f[2].parse()?, depth: f[3].parse()?, mag: f[4].parse()?, name })
    }).collect()
}

fn distance_km(a: (f64, f64), b: (f64, f64)) -> f64 {
    let (la1, lo1, la2, lo2) = (a.0.to_radians(), a.1.to_radians(), b.0.to_radians(), b.1.to_radians());
    6371.0 * (la1.sin() * la2.sin() + la1.cos() * la2.cos() * (lo1 - lo2).cos()).clamp(-1.0, 1.0).acos()
}

/// One candidate month's scan: the signed Z map of the radius with the larger |Z| per cell,
/// before standardisation, for each radius.
struct Scanner<'a> {
    spans: &'a [Span],
    kappa: &'a [f64],
    /// The M3 signal per degree: (standard deviation, correlation time in years).
    signal: &'a [(f64, f64)],
    templates: Vec<Vec<f64>>,
    lats: &'a [f64],
    lons: &'a [f64],
}

impl Scanner<'_> {
    /// Step kernel for epoch `tau`: weights over months, or None if the window is too thin.
    fn kernel(&self, tau: f64) -> Result<Option<Vec<(usize, f64)>>> {
        let win: Vec<usize> = (0..self.spans.len()).filter(|&i| (self.spans[i].mid() - tau).abs() <= WINDOW_YEARS).collect();
        let before = win.iter().filter(|&&i| self.spans[i].mid() < tau).count();
        if before < 6 || win.len() - before < 6 {
            return Ok(None);
        }
        let row = |i: usize| {
            let sp = self.spans[i];
            let t = sp.mid();
            let w = |p: f64| core::f64::consts::TAU * t / p;
            let step = ((sp.end - tau) / (sp.end - sp.start)).clamp(0.0, 1.0);
            vec![1.0, t - tau, w(1.0).cos(), w(1.0).sin(), w(0.5).cos(), w(0.5).sin(), step]
        };
        let x: Vec<Vec<f64>> = win.iter().map(|&i| row(i)).collect();
        let mut k = Vec::with_capacity(win.len());
        for (j, &i) in win.iter().enumerate() {
            let unit: Vec<f64> = (0..win.len()).map(|r| if r == j { 1.0 } else { 0.0 }).collect();
            k.push((i, ols(&x, &unit)?.0[6]));
        }
        Ok(Some(k))
    }

    /// Z maps (one per radius) for data `d` and errors `sigma` with the step kernel `k`.
    fn scan(&self, k: &[(usize, f64)], d: &[Flat], sigma: &[Flat], leg: &mut Legendre) -> Vec<Vec<f64>> {
        let size = idx(NMAX, NMAX) + 1;
        let (mut sc, mut ss, mut vc, mut vs) = (vec![0.0; size], vec![0.0; size], vec![0.0; size], vec![0.0; size]);
        for &(i, w) in k {
            for j in 0..size {
                sc[j] += w * d[i].c[j];
                ss[j] += w * d[i].s[j];
                vc[j] += w * w * sigma[i].c[j] * sigma[i].c[j];
                vs[j] += w * w * sigma[i].s[j] * sigma[i].s[j];
            }
        }
        // The natural signal's share of the step variance, per degree: Σᵢⱼ kᵢkⱼ sd² e^(−|tᵢ−tⱼ|/τ).
        let natural: Vec<f64> = self.signal.iter().map(|&(sd, tau)| {
            k.iter().flat_map(|&(i, wi)| k.iter().map(move |&(j, wj)| (i, wi, j, wj)))
                .map(|(i, wi, j, wj)| wi * wj * (-(self.spans[i].mid() - self.spans[j].mid()).abs() / tau).exp())
                .sum::<f64>() * sd * sd
        }).collect();
        self.templates.iter().map(|tmpl| {
            let (mut ac, mut asn, mut bc, mut bs) = (vec![0.0; size], vec![0.0; size], vec![0.0; size], vec![0.0; size]);
            for (n, (&t, &kap)) in tmpl.iter().zip(self.kappa).enumerate().skip(NMIN) {
                let k2 = kap * kap;
                for m in 0..=n {
                    let j = idx(n, m);
                    let (v_c, v_s) = (natural[n] + k2 * vc[j], natural[n] + k2 * vs[j]);
                    if vc[j] > 0.0 {
                        ac[j] = t * sc[j] / v_c;
                        bc[j] = t * t / v_c;
                    }
                    if m > 0 && vs[j] > 0.0 {
                        asn[j] = t * ss[j] / v_s;
                        bs[j] = t * t / v_s;
                    }
                }
            }
            let num = synthesize(&ac, &asn, false, self.lats, self.lons, leg);
            let den = synthesize(&bc, &bs, true, self.lats, self.lons, leg);
            num.iter().zip(&den).map(|(a, b)| if *b > 0.0 { a / b.sqrt() } else { 0.0 }).collect()
        }).collect()
    }
}

#[allow(clippy::too_many_arguments)] // one argument per input file, as on the command line
pub fn run(months_dir: &str, baseline_dir: &str, tn14: &str, tn13: &str, kalman_csv: &Path, catalogue: &Path, figures: &str) -> Result<String> {
    let tn14 = read_tn14(&std::fs::read_to_string(tn14)?)?;
    let tn13 = read_tn13(&std::fs::read_to_string(tn13)?)?;
    let (mean, _) = baseline_mean(baseline_dir, &tn14, &tn13)?;
    let (mut spans, mut data, mut sigma) = (Vec::new(), Vec::new(), Vec::new());
    for p in gfc_files(months_dir)? {
        let s = read_gfc(&std::fs::read_to_string(&p)?)?.sigma.ok_or("a month without errors")?;
        data.push(flatten(&gravimetry::anomaly(&corrected_month(&p, &tn14, &tn13)?, &mean).truncated(NMAX)));
        sigma.push(flatten(&s.rescaled(mean.gm, mean.radius).truncated(NMAX)));
        spans.push(span(&p)?);
    }
    // The M3 Kalman model per degree: signal sd, correlation time (months) and κ (Q2).
    let mut kappa = vec![1.0; NMAX + 1];
    let mut signal = vec![(0.0, 1.0); NMAX + 1];
    for line in std::fs::read_to_string(kalman_csv)?.lines().skip(1) {
        let f: Vec<f64> = line.split(',').filter_map(|s| s.parse().ok()).collect();
        if let [n, sd, tau, k, ..] = f[..] {
            if (n as usize) <= NMAX {
                kappa[n as usize] = k;
                signal[n as usize] = (sd, tau / 12.0);
            }
        }
    }
    let quakes = read_catalogue(catalogue)?;
    let (lats, lons) = grid_axes(2.0);
    let cells = lats.len() * lons.len();
    let love = gravimetry::love_numbers(NMAX).ok_or("Love numbers")?;
    let scanner = Scanner {
        spans: &spans, kappa: &kappa, signal: &signal, lats: &lats, lons: &lons,
        templates: RADII_KM.iter().map(|&r| cap_template(r / 6371.0, mean.radius, &love)).collect(),
    };
    let mut leg = Legendre::new(NMAX);

    // Candidate epochs: the first day of every month.
    let (y0, m0) = calendar(spans[0].mid());
    let (y1, m1) = calendar(spans[spans.len() - 1].mid());
    let mut epochs = Vec::new();
    let (mut y, mut m) = (y0, m0);
    while (y, m) <= (y1, m1) {
        epochs.push(((y, m), geopotential::julian_epoch(y, m, 1, 0, 0)));
        (y, m) = if m == 12 { (y + 1, 1) } else { (y, m + 1) };
    }
    let mut kernels = Vec::new();
    for &(cm, tau) in &epochs {
        if let Some(k) = scanner.kernel(tau)? {
            kernels.push((cm, tau, k));
        }
    }
    eprintln!("m1: {} months, {} candidate epochs", spans.len(), kernels.len());

    // Residuals after a global seasonal + trend fit per coefficient: the scan's data, and the
    // null's raw material.
    let size = idx(NMAX, NMAX) + 1;
    let design = |t: f64| { let w = |p: f64| core::f64::consts::TAU * t / p; vec![1.0, t - 2015.0, w(1.0).cos(), w(1.0).sin(), w(0.5).cos(), w(0.5).sin()] };
    let x: Vec<Vec<f64>> = spans.iter().map(|s| design(s.mid())).collect();
    let mut resid: Vec<Flat> = data.clone();
    for j in 0..size {
        for part in 0..2 {
            let y: Vec<f64> = data.iter().map(|d| if part == 0 { d.c[j] } else { d.s[j] }).collect();
            if y.iter().all(|v| *v == 0.0) { continue; }
            let (_, e) = ols(&x, &y)?;
            for (i, e) in e.iter().enumerate() {
                if part == 0 { resid[i].c[j] = *e } else { resid[i].s[j] = *e }
            }
        }
    }
    let maps: Vec<Vec<Vec<f64>>> = kernels.iter().map(|(_, _, k)| scanner.scan(k, &resid, &sigma, &mut leg)).collect();
    // Each cell's robust spread across epochs, per radius.
    let scale: Vec<Vec<f64>> = (0..RADII_KM.len()).map(|r| (0..cells).map(|c| {
        let v: Vec<f64> = maps.iter().map(|m| m[r][c].abs()).collect();
        1.4826 * quantile(&v, 0.5)
    }).collect()).collect();
    // Standardised, combined over radii: the signed value with the largest |Z'|, and its radius.
    let combine = |m: &[Vec<f64>]| -> Vec<(f64, usize)> {
        (0..cells).map(|c| (0..RADII_KM.len()).map(|r| (m[r][c] / scale[r][c], r)).max_by(|a, b| a.0.abs().total_cmp(&b.0.abs())).unwrap_or((0.0, 0))).collect()
    };
    let zmaps: Vec<Vec<(f64, usize)>> = maps.iter().map(|m| combine(m)).collect();
    let peak = |z: &[(f64, usize)]| -> (usize, f64) {
        z.iter().enumerate().map(|(c, v)| (c, v.0.abs())).max_by(|a, b| a.1.total_cmp(&b.1)).unwrap_or((0, 0.0))
    };

    // Null: shift the residual months circularly by at least two years, keeping their order.
    let mut rng = Rng::new(2025);
    let mut null_max = Vec::new();
    let months = spans.len();
    for _ in 0..SHUFFLES {
        let shift = 24 + rng.below(months - 48);
        let (d, s): (Vec<Flat>, Vec<Flat>) = (0..months).map(|i| (resid[(i + shift) % months].clone(), sigma[(i + shift) % months].clone())).unzip();
        for e in 0..NULL_EPOCHS {
            let (_, _, k) = &kernels[e * kernels.len() / NULL_EPOCHS];
            null_max.push(peak(&combine(&scanner.scan(k, &d, &s, &mut leg))).1);
        }
    }
    // The shifted null still holds the record's real earthquakes, now at other times, so its top
    // tail is contaminated. Maxima of many cells follow a Gumbel law; fit it to the null's
    // quartiles, which a few contaminated values cannot move, and take its 99th percentile.
    let (q25, q75) = (quantile(&null_max, 0.25), quantile(&null_max, 0.75));
    let beta = (q75 - q25) / ((-(0.75f64).ln()).ln() - (-(0.25f64).ln()).ln()).abs();
    let mu = q25 + beta * (-(0.25f64).ln()).ln();
    let threshold = mu - beta * (-(0.99f64).ln()).ln();

    // Detections: epochs whose global maximum passes. The ±1-year window cannot time a step better
    // than about a year, so detections within a year and 1,000 km of each other are one event.
    let coord = |c: usize| (lats[c / lons.len()].to_degrees(), lons[c % lons.len()].to_degrees());
    let mut detections: Vec<(usize, usize, f64, usize)> = Vec::new(); // (epoch, cell, signed Z', radius)
    for (e, z) in zmaps.iter().enumerate() {
        let (c, v) = peak(z);
        if v < threshold { continue; }
        let tau = kernels[e].1;
        if let Some(prev) = detections.iter_mut().find(|d| (kernels[d.0].1 - tau).abs() <= WINDOW_YEARS && distance_km(coord(d.1), coord(c)) < 1000.0) {
            if v > prev.2.abs() { *prev = (e, c, z[c].0, z[c].1); }
            continue;
        }
        detections.push((e, c, z[c].0, z[c].1));
    }
    // Score each catalogue event: the largest |Z'| within 600 km and two months.
    let mut rep = String::new();
    let _ = writeln!(rep, "{} monthly fields, {} candidate epochs ({}-{:02} to {}-{:02}), caps of {:?} km, degrees {NMIN}–{NMAX}", spans.len(), kernels.len(), kernels[0].0.0, kernels[0].0.1, kernels[kernels.len() - 1].0.0, kernels[kernels.len() - 1].0.1, RADII_KM);
    let _ = writeln!(rep, "noise: the M3 model (natural signal + file errors × κ); median robust spread of Z per cell {:.2} (1 = the model is exact)", quantile(&scale[0], 0.5));
    let _ = writeln!(rep, "null: {} circularly shifted epoch maxima (median {:.2}, raw 99th percentile {:.2}); Gumbel fit to the quartiles gives the 1%-per-month threshold |Z'| = {threshold:.2}", null_max.len(), quantile(&null_max, 0.5), quantile(&null_max, 0.99));
    let _ = writeln!(rep, "\ncatalogue events (USGS Mw ≥ 8.0, {}):", quakes.len());
    let mut found = 0;
    let mut scored = 0;
    for q in &quakes {
        let near = kernels.iter().enumerate().filter(|(_, (_, tau, _))| (tau - q.t).abs() < 2.0 / 12.0)
            .flat_map(|(e, _)| (0..cells).filter(|&c| distance_km(coord(c), (q.lat, q.lon)) < 600.0).map(move |c| (e, c)))
            .map(|(e, c)| (e, c, zmaps[e][c].0.abs())).max_by(|a, b| a.2.total_cmp(&b.2));
        let Some((e, _, z)) = near else {
            let _ = writeln!(rep, "  {:<48} Mw {:.1} {:>4.0} km deep: no field window", q.name, q.mag, q.depth);
            continue;
        };
        scored += 1;
        let hit = z >= threshold;
        found += usize::from(hit);
        let _ = writeln!(rep, "  {:<48} Mw {:.1} {:>4.0} km deep: max |Z'| {z:5.2} at {}-{:02} {}", q.name, q.mag, q.depth, kernels[e].0.0, kernels[e].0.1, if hit { "DETECTED" } else { "" });
    }
    let _ = writeln!(rep, "recovered {found} of {scored} scorable events");
    let _ = writeln!(rep, "\nblind detections (global monthly maximum above the threshold), {}:", detections.len());
    let mut matched = 0;
    for &(e, c, z, r) in &detections {
        let (la, lo) = coord(c);
        let tau = kernels[e].1;
        let q = quakes.iter().find(|q| (q.t - tau).abs() < 2.5 / 12.0 && distance_km((q.lat, q.lon), (la, lo)) < 1000.0);
        matched += usize::from(q.is_some());
        let amplitude = z * scale[r][c]; // back to the unstandardised Z for the amplitude scale
        let _ = writeln!(rep, "  {}-{:02} at {la:5.1}°, {lo:6.1}°: Z' {z:+.2} ({:.0} km cap, raw Z {amplitude:+.1}) {}", kernels[e].0.0, kernels[e].0.1, RADII_KM[r], q.map_or("— not in the catalogue".to_string(), |q| format!("= {} (Mw {:.1})", q.name, q.mag)));
    }
    let years = (kernels[kernels.len() - 1].1 - kernels[0].1).max(1.0);
    let _ = writeln!(rep, "{matched} of {} detections match a catalogue earthquake; the other {} are {:.2} per year (hydrology, ice or noise)", detections.len(), detections.len() - matched, (detections.len() - matched) as f64 / years);

    write(figures, "m1_report.txt", &rep)?;
    let svg = figure(&kernels.iter().map(|k| k.1).collect::<Vec<_>>(), &zmaps.iter().map(|z| peak(z).1).collect::<Vec<_>>(), threshold, &quakes, &detections.iter().map(|&(e, c, z, _)| (kernels[e].1, coord(c), z)).collect::<Vec<_>>());
    write(figures, "m1_detector.svg", svg)?;
    Ok(rep)
}

fn figure(t: &[f64], peak: &[f64], threshold: f64, quakes: &[Quake], detections: &[(f64, (f64, f64), f64)]) -> String {
    let (w, h) = (960.0, 760.0);
    let mut s = open(w, h, "Blind detector of sudden mass changes");
    let _ = write!(s, r#"<text x="60" y="32" class="title">A blind search of 24 years of gravity for sudden mass changes</text>
<text x="60" y="54" class="subtitle">Each month: the strongest standardised step anywhere on Earth (likelihood-ratio scan for 400 and 800 km caps). Dashed: the false-alarm threshold</text>
<text x="60" y="72" class="subtitle">from time-shifted data. Orange: USGS earthquakes of magnitude 8 or more (taller for larger). The catalogue was not used in the search.</text>
"#);
    let (x0, x1, y0, y1) = (70.0, 930.0, 100.0, 330.0);
    let (t0, t1) = (2002.0, 2026.8);
    let top = peak.iter().fold(threshold, |a, &b| a.max(b)).ceil() + 1.0;
    let x = |v: f64| x0 + (x1 - x0) * (v - t0) / (t1 - t0);
    let y = |v: f64| y1 - (y1 - y0) * v / top;
    for q in quakes {
        let hq = (q.mag - 7.5) * 60.0;
        let _ = writeln!(s, r#"<line x1="{0:.1}" x2="{0:.1}" y1="{1:.1}" y2="{2:.1}" stroke="var(--s2)" stroke-width="2" opacity="0.55"><title>{3} Mw {4:.1}</title></line>"#, x(q.t), y1, y1 - hq, esc(&q.name), q.mag);
        if q.mag >= 8.5 {
            let _ = writeln!(s, r#"<text x="{:.1}" y="{:.1}" class="tick" text-anchor="middle" style="fill: var(--s2)">{:.1}</text>"#, x(q.t), y1 - hq - 4.0, q.mag);
        }
    }
    let mut v = 0.0;
    while v <= top {
        let _ = writeln!(s, r#"<line x1="{x0}" x2="{x1}" y1="{0:.1}" y2="{0:.1}" class="grid"/><text x="{1:.1}" y="{2:.1}" class="tick" text-anchor="end">{v:.0}</text>"#, y(v), x0 - 6.0, y(v) + 4.0);
        v += 2.0;
    }
    for yr in (2004..=2026).step_by(2) {
        let _ = writeln!(s, r#"<text x="{:.1}" y="{:.1}" class="tick" text-anchor="middle">{yr}</text>"#, x(yr as f64), y1 + 16.0);
    }
    let d: String = t.iter().zip(peak).enumerate().map(|(i, (&ti, &p))| format!("{}{:.1},{:.1}", if i == 0 || ti - t[i - 1] > 0.25 { 'M' } else { 'L' }, x(ti), y(p))).collect();
    let _ = writeln!(s, r#"<path d="{d}" class="line ink" stroke-width="1.4"/>"#);
    let _ = writeln!(s, r#"<line x1="{x0}" x2="{x1}" y1="{0:.1}" y2="{0:.1}" stroke="var(--ink)" stroke-dasharray="5 4"/><text x="{x1}" y="{1:.1}" class="tick" text-anchor="end">threshold {threshold:.1}</text>"#, y(threshold), y(threshold) - 5.0);
    let _ = writeln!(s, r#"<text x="{x0}" y="{:.1}" class="tick">largest |Z'| in the month</text>"#, y0 - 6.0);

    // Map of detections and catalogue events (equirectangular).
    let (mx0, mx1, my0, my1) = (70.0, 930.0, 380.0, 710.0);
    let mx = |lon: f64| mx0 + (mx1 - mx0) * (lon + 180.0) / 360.0;
    let my = |lat: f64| my0 + (my1 - my0) * (75.0 - lat) / 140.0;
    let _ = writeln!(s, r#"<text x="{mx0}" y="{:.1}" class="panel">Where: blind detections (blue rings, size by strength) and the catalogue (orange)</text>"#, my0 - 10.0);
    let _ = writeln!(s, r#"<rect x="{mx0}" y="{my0}" width="{:.1}" height="{:.1}" fill="none" class="axis"/>"#, mx1 - mx0, my1 - my0);
    for ring in land_rings() {
        let pts: Vec<String> = ring.iter().filter(|(_, la)| (-65.0..75.0).contains(la)).map(|&(lo, la)| format!("{:.1},{:.1}", mx(lo), my(la))).collect();
        if pts.len() > 2 {
            let _ = writeln!(s, r#"<polygon points="{}" fill="var(--grid)" stroke="none"/>"#, pts.join(" "));
        }
    }
    for q in quakes {
        let _ = writeln!(s, r#"<circle cx="{:.1}" cy="{:.1}" r="{:.1}" fill="var(--s2)" opacity="0.7"><title>{} Mw {:.1}</title></circle>"#, mx(q.lon), my(q.lat), (q.mag - 7.6) * 6.0, esc(&q.name), q.mag);
    }
    for &(tt, (la, lo), z) in detections {
        let (yy, mm) = calendar(tt);
        let _ = writeln!(s, r#"<circle cx="{:.1}" cy="{:.1}" r="{:.1}" fill="none" stroke="var(--s1)" stroke-width="2"><title>{yy}-{mm:02}: Z' {z:+.1}</title></circle>"#, mx(lo), my(la), 4.0 + z.abs());
    }
    s.push_str("</svg>\n");
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cap_template_matches_numerical_integration() {
        // σ_n = (2n+1)/2 ∫ P_n(x) dx over the cap, by Simpson's rule, against the closed form.
        let love = vec![0.0; NMAX + 1];
        let (radius, psi) = (6_371_000.0, 800.0 / 6371.0);
        let c = cap_template(psi, radius, &love);
        let legendre = |n: usize, x: f64| {
            let (mut p0, mut p1) = (1.0, x);
            for k in 1..n { (p0, p1) = (p1, ((2 * k + 1) as f64 * x * p1 - k as f64 * p0) / (k + 1) as f64); }
            if n == 0 { p0 } else { p1 }
        };
        let (a, b, steps) = (psi.cos(), 1.0, 4000);
        for n in [2, 7, 30, 60] {
            let h = (b - a) / steps as f64;
            let integral: f64 = (0..=steps).map(|k| {
                let w = if k == 0 || k == steps { 1.0 } else if k % 2 == 1 { 4.0 } else { 2.0 };
                w * legendre(n, a + k as f64 * h)
            }).sum::<f64>() * h / 3.0;
            let sigma = (2 * n + 1) as f64 / 2.0 * integral;
            let from_template = c[n] * radius * RHO_EARTH * ((2 * n + 1) as f64).powi(2) / 3.0;
            assert!((sigma - from_template).abs() < 1e-9, "n = {n}: {sigma} vs {from_template}");
        }
    }
}
