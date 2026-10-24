//! Closed-loop simulation of GRACE-type mission classes (Zingerle et al. 2024, Fig. 2a analogue:
//! product noise only, single pair, one month) against a real month's calibrated errors.

use std::fmt::Write;

use geopotential::{Coefficients, Functional};
use gravimetry::missions::{Accelerometer, Ranging, pair_noise};
use gravimetry::sim::{Pair, Setup, simulate};

use crate::svg::{esc, open, pow10};

const NMAX: usize = 40;
const DAYS: f64 = 31.0;
/// 2¹⁶ samples over the month (40.9 s; Nyquist ≈ degree 69), so noise whitening is exact.
const DT: f64 = DAYS * 86_400.0 / 65_536.0;

struct Curve {
    name: String,
    class: &'static str,
    mm: Vec<f64>,
    /// χ²/K of the errors against the formal covariance: 1 ± √(2/K) when the loop is consistent.
    consistency: Option<f64>,
}

fn ewh_mm(c: &Coefficients) -> Vec<f64> {
    let love = gravimetry::love_numbers(c.max_degree()).expect("degree ≤ 200");
    c.degree_amplitudes(&Functional::WaterHeight { love: &love }).into_iter().map(|m| m * 1e3).collect()
}

/// Run the four classes, then render the figure (SVG) and its CSV.
pub fn figure(reference: Option<(&Coefficients, &str)>) -> (String, String) {
    let (gm, radius) = (3.986004415e14, 6378136.3);
    let pair = Pair { altitude: 463e3, inclination: 89f64.to_radians(), separation: 200e3 };
    let cai = |white| Accelerometer::AtomInterferometer { white };
    let classes: [(&str, &str, Accelerometer, Ranging); 4] = [
        ("GRACE-FO class: SuperSTAR + LRI GRACE-FO", "ink", Accelerometer::SuperStar, Ranging::GraceFo),
        ("NGGM class: MicroSTAR + LRI 2030", "s1", Accelerometer::MicroStar, Ranging::Lri2030),
        ("CAI11 class: CAI 1e-11 + LRI 2033", "s2", cai(1e-11), Ranging::Lri2033),
        ("CAI12 class: CAI 1e-12 + LRI 2033", "s3", cai(1e-12), Ranging::Lri2033),
    ];
    let mut curves: Vec<Curve> = Vec::new();
    for (i, (name, class, acc, rng)) in classes.into_iter().enumerate() {
        let t0 = std::time::Instant::now();
        let setup = Setup { pair, gm, radius, nmax: NMAX, duration: DAYS * 86_400.0, dt: DT, seed: 100 + i as u64 };
        let sol = simulate(&setup, move |f| pair_noise(acc, rng, f), None);
        let consistency = sol.chi2_per_dof;
        eprintln!("{name}: {:.0} s, chi2/K = {consistency:.3}", t0.elapsed().as_secs_f64());
        curves.push(Curve { name: name.into(), class, mm: ewh_mm(&sol.errors), consistency: Some(consistency) });
    }
    if let Some((sigma, label)) = reference {
        let name = format!("Real GRACE-FO month: COST-G {label}, calibrated errors");
        curves.push(Curve { name, class: "muted", mm: ewh_mm(&sigma.truncated(NMAX)), consistency: None });
    }

    let (lo, hi) = curves.iter().flat_map(|c| c.mm[2..].iter()).fold((f64::MAX, f64::MIN), |(a, b), &v| (a.min(v), b.max(v)));
    let (e_lo, e_hi) = (lo.log10().floor() as i32, hi.log10().ceil() as i32);
    let (w, left, right, top, bottom) = (960.0, 92.0, 936.0, 164.0, 524.0);
    let x = |n: f64| left + (right - left) * (n - 2.0) / (NMAX as f64 - 2.0);
    let y = |v: f64| top + (bottom - top) * (f64::from(e_hi) - v.log10()) / f64::from(e_hi - e_lo);

    let mut s = open(w, 630.0, "Closed-loop errors of GRACE-type mission classes");
    let _ = write!(
        s,
        r#"<text x="{left}" y="32" class="title">Closed-loop errors of GRACE-type mission classes</text>
<text x="{left}" y="54" class="subtitle">Error per degree, mm of equivalent water height, of one month recovered from simulated noise only (single pair, 89°, 463 km, 200 km apart,</text>
<text x="{left}" y="72" class="subtitle">41 s sampling, degree {NMAX}, whitened least squares). Compare Zingerle et al. (2024), Fig. 2a; time-variable gravity and aliasing are not simulated.</text>
"#
    );
    for (i, c) in curves.iter().enumerate() {
        let (lx, ly) = (left + 430.0 * (i % 2) as f64, 98.0 + 20.0 * (i / 2) as f64);
        let label = if c.consistency.is_none() { &c.name } else { c.name.split(':').next().unwrap_or(&c.name) };
        let _ = writeln!(s, r#"<line x1="{lx}" y1="{ly}" x2="{}" y2="{ly}" class="line {}"/><text x="{}" y="{}">{}</text>"#, lx + 16.0, c.class, lx + 22.0, ly + 4.0, esc(label));
    }
    for e in e_lo..=e_hi {
        let yy = y(10f64.powi(e));
        let _ = writeln!(s, r#"<line x1="{left}" y1="{yy:.1}" x2="{right}" y2="{yy:.1}" class="grid"/><text x="{:.1}" y="{:.1}" class="tick" text-anchor="end">{} mm</text>"#, left - 8.0, yy + 4.0, pow10(e));
    }
    for n in (10..=NMAX).step_by(10).chain([2]) {
        let xx = x(n as f64);
        let _ = writeln!(s, r#"<line x1="{xx:.1}" y1="{top}" x2="{xx:.1}" y2="{bottom}" class="grid"/><text x="{xx:.1}" y="{}" class="tick" text-anchor="middle">{n}</text>"#, bottom + 18.0);
    }
    let _ = writeln!(s, r#"<line x1="{left}" y1="{bottom}" x2="{right}" y2="{bottom}" class="axis"/><text x="{}" y="{}" class="tick" text-anchor="middle">Spherical-harmonic degree</text>"#, (left + right) / 2.0, bottom + 38.0);
    for c in &curves {
        let d: String = (2..=NMAX).map(|n| format!("{}{:.1},{:.1}", if n == 2 { 'M' } else { 'L' }, x(n as f64), y(c.mm[n]))).collect();
        let _ = writeln!(s, r#"<path d="{d}" class="line {}"><title>{}</title></path>"#, c.class, esc(&c.name));
    }
    let unknowns = (NMAX + 1) * (NMAX + 1) - 4;
    let consist: Vec<String> = curves.iter().filter_map(|c| c.consistency.map(|r| format!("{r:.2}"))).collect();
    let _ = write!(
        s,
        r#"<text x="{left}" y="{}" class="subtitle">Consistency of actual and formal errors, χ²/K per class: {} (expected 1 ± {:.3} for {} unknowns).</text>
</svg>
"#,
        bottom + 66.0,
        consist.join(", "),
        (2.0 / unknowns as f64).sqrt(),
        unknowns
    );

    let mut csv = String::from("degree");
    for c in &curves {
        let _ = write!(csv, ",\"{} (mm)\"", c.name);
    }
    csv.push('\n');
    for n in 2..=NMAX {
        let _ = write!(csv, "{n}");
        for c in &curves {
            let _ = write!(csv, ",{:.4e}", c.mm[n]);
        }
        csv.push('\n');
    }
    (s, csv)
}
