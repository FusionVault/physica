//! Static-field geoid errors of gravity-gradiometer missions against a real gravity field's signal,
//! in the manner of Zingerle et al. (2024), Figs. 5 and 7.

use std::fmt::Write;

use geopotential::{Coefficients, Functional};
use satgrav::{AtomInterferometer, gradiometry_errors, resolution};

use crate::svg::{esc, open, pow10};

const NMAX: usize = 300;
const Y_MIN: f64 = 1e-5;
const Y_MAX: f64 = 1e2;
const MONTH: f64 = 365.25 / 12.0 * 86_400.0;

struct Curve {
    name: String,
    class: &'static str,
    values: Vec<f64>,
}

/// Render the figure (SVG) and its table view (CSV) for `signal` (e.g. GOCO06s).
pub fn figure(signal: &Coefficients, signal_name: &str) -> (String, String) {
    let (gm, radius) = (signal.gm, signal.radius);
    let geoid_error = |sigma: Vec<f64>| -> Vec<f64> {
        sigma.iter().enumerate().map(|(n, s)| radius * ((2 * n + 1) as f64).sqrt() * s).collect()
    };
    // Trimeche et al. (2019): T = 5 s, N = 1e6, D = 0.5 m, interleaved ~1 s cycle → ~3.5 mE/√Hz;
    // their mission simulation uses 5 mE/√Hz white, which is what is plotted.
    let cai = AtomInterferometer { t: 5.0, contrast: 1.0, atoms: 1e6, cycle: 1.0 };
    debug_assert!((cai.gradient_asd(0.5) - 3.5e-12).abs() < 0.1e-12);
    let goce_asd = |f: f64| 11e-12 * (1.0 + 5e-3 / f).sqrt();
    let signal_amp = signal.degree_amplitudes(&Functional::GeoidHeight);
    let curves = [
        Curve { name: format!("{signal_name} signal"), class: "ink", values: signal_amp[..=NMAX.min(signal.max_degree())].to_vec() },
        Curve { name: "GOCE-like: 11 mE/√Hz (rising below 5 mHz), 255 km, 47 months".into(), class: "s1",
            values: geoid_error(gradiometry_errors(gm, radius, 255e3, 47.0 * MONTH, NMAX, goce_asd)) },
        Curve { name: "Cold-atom gradiometer: 5 mE/√Hz, 239 km, 8 months".into(), class: "s2",
            values: geoid_error(gradiometry_errors(gm, radius, 239e3, 8.0 * MONTH, NMAX, |_| 5e-12)) },
        Curve { name: "Cold-atom gradiometer: 5 mE/√Hz, 239 km, 47 months".into(), class: "s3",
            values: geoid_error(gradiometry_errors(gm, radius, 239e3, 47.0 * MONTH, NMAX, |_| 5e-12)) },
    ];
    // Highest degree resolved: the first degree past 20 where the error exceeds the signal.
    let crossing = |c: &Curve| (20..=NMAX.min(signal.max_degree())).find(|&n| c.values[n] >= curves[0].values[n]);

    let (w, h, left, right, top, bottom) = (960.0, 640.0, 92.0, 936.0, 176.0, 540.0);
    let x = |n: f64| left + (right - left) * n / NMAX as f64;
    let y = |v: f64| top + (bottom - top) * (Y_MAX / v.clamp(1e-30, f64::MAX)).log10() / (Y_MAX / Y_MIN).log10();

    let mut s = open(w, h, "Static gravity-field errors of gradiometer missions");
    let _ = write!(
        s,
        r#"<text x="{left}" y="32" class="title">Static gravity-field errors of gradiometer missions</text>
<text x="{left}" y="54" class="subtitle">Geoid degree amplitudes, m: the {sig} signal against formal errors of a radial-gradient (V_zz) gradiometer, propagated in closed form</text>
<text x="{left}" y="72" class="subtitle">(uniform coverage, no polar gap, one tensor component). Compare Zingerle et al. (2024), Figs. 5 and 7, and Trimeche et al. (2019).</text>
<clipPath id="plot"><rect x="{left}" y="{top}" width="{pw}" height="{ph}"/></clipPath>
"#,
        sig = esc(signal_name),
        pw = right - left,
        ph = bottom - top
    );
    // Legend: two rows, short line keys, names in text ink.
    for (i, c) in curves.iter().enumerate() {
        let (lx, ly) = (left + 430.0 * (i % 2) as f64, 98.0 + 20.0 * (i / 2) as f64);
        let _ = writeln!(s, r#"<line x1="{lx}" y1="{ly}" x2="{}" y2="{ly}" class="line {}"/><text x="{}" y="{}">{}</text>"#, lx + 16.0, c.class, lx + 22.0, ly + 4.0, esc(&c.name));
    }
    // Top axis: half-wavelength resolution.
    let _ = writeln!(s, r#"<text x="{left}" y="146" class="tick">Resolution (half wavelength), km</text><line x1="{left}" y1="{top}" x2="{right}" y2="{top}" class="axis"/>"#);
    for n in [50, 100, 150, 200, 250, 300] {
        let xx = x(n as f64);
        let _ = writeln!(s, r#"<line x1="{xx:.1}" y1="{top}" x2="{xx:.1}" y2="{}" class="axis"/><text x="{xx:.1}" y="{}" class="tick" text-anchor="middle">{:.0}</text>"#, top - 4.0, top - 8.0, resolution(n, radius) / 1e3);
    }
    for e in -5..=2 {
        let yy = y(10f64.powi(e));
        let _ = writeln!(s, r#"<line x1="{left}" y1="{yy:.1}" x2="{right}" y2="{yy:.1}" class="grid"/><text x="{:.1}" y="{:.1}" class="tick" text-anchor="end">{}</text>"#, left - 8.0, yy + 4.0, pow10(e));
    }
    for n in (0..=NMAX).step_by(50) {
        let xx = x(n as f64);
        let _ = writeln!(s, r#"<line x1="{xx:.1}" y1="{top}" x2="{xx:.1}" y2="{bottom}" class="grid"/><text x="{xx:.1}" y="{}" class="tick" text-anchor="middle">{n}</text>"#, bottom + 18.0);
    }
    let _ = writeln!(s, r#"<line x1="{left}" y1="{bottom}" x2="{right}" y2="{bottom}" class="axis"/><text x="{}" y="{}" class="tick" text-anchor="middle">Spherical-harmonic degree</text>"#, (left + right) / 2.0, bottom + 38.0);
    for c in &curves {
        let d: String = c.values.iter().enumerate().skip(2)
            .map(|(n, &v)| format!("{}{:.1},{:.1}", if n == 2 { 'M' } else { 'L' }, x(n as f64), y(v)))
            .collect();
        let _ = writeln!(s, r#"<path d="{d}" class="line {}" clip-path="url(#plot)"><title>{}</title></path>"#, c.class, esc(&c.name));
    }
    // Where each error curve meets the signal: a dot with a surface ring, numbers in the footnote.
    let mut found = Vec::new();
    for c in &curves[1..] {
        if let Some(n) = crossing(c) {
            let _ = writeln!(s, r#"<circle cx="{:.1}" cy="{:.1}" r="5" fill="var(--{})" stroke="var(--surface)" stroke-width="2"><title>{}: degree {n}, {:.0} km</title></circle>"#,
                x(n as f64), y(curves[0].values[n]), c.class, esc(&c.name), resolution(n, radius) / 1e3);
            found.push(format!("degree {n} ({:.0} km)", resolution(n, radius) / 1e3));
        }
    }
    let _ = write!(
        s,
        r#"<text x="{left}" y="{}" class="subtitle">Error meets signal at {}: GOCE-like, cold-atom 8 months and 47 months.</text>
</svg>
"#,
        bottom + 70.0,
        found.join(", ")
    );

    let mut csv = String::from("degree,resolution_km");
    for c in &curves {
        let _ = write!(csv, ",\"{} (m)\"", c.name);
    }
    csv.push('\n');
    for n in 2..=NMAX.min(signal.max_degree()) {
        let _ = write!(csv, "{n},{:.1}", resolution(n, radius) / 1e3);
        for c in &curves {
            let _ = write!(csv, ",{:.4e}", c.values[n]);
        }
        csv.push('\n');
    }
    (s, csv)
}
