//! Figure 1 of Zingerle et al. (2024): accelerometer, laser-ranging and combined noise spectra of
//! GRACE-type missions, as three panels on one shared log-frequency scale.

use std::fmt::Write;

use satgrav::{Accelerometer, Ranging, pair_noise};

use crate::svg::{esc, open, pow10};

const F_MIN: f64 = 1.4e-4;
const F_MAX: f64 = 0.1;
const Y_MIN: f64 = 5e-13;
const Y_MAX: f64 = 1e-9;
/// Orbit period of the paper's 463 km GRACE-like orbit, s.
const T_ORBIT: f64 = 5620.0;
/// The paper draws the CAI accelerometers only up to this frequency.
const CAI_F_MAX: f64 = 0.03;

struct Series {
    name: &'static str,
    class: &'static str,
    values: Vec<f64>,
    f_max: f64,
    /// Direct label: frequency and vertical offset in px (negative = above the line).
    label: Option<(f64, f64)>,
}

const CAI11: Accelerometer = Accelerometer::AtomInterferometer { white: 1e-11 };
const CAI12: Accelerometer = Accelerometer::AtomInterferometer { white: 1e-12 };

/// The three panels' series, evaluated on `freqs`.
fn panels(freqs: &[f64]) -> [(&'static str, Vec<Series>); 3] {
    let acc = |name, class, a: Accelerometer, f_max| Series {
        name,
        class,
        values: freqs.iter().map(|&f| std::f64::consts::SQRT_2 * a.asd(f)).collect(),
        f_max,
        label: None, // seven converging lines: the legend carries identity
    };
    let lri = |name, class, r: Ranging| Series {
        name,
        class,
        values: freqs.iter().map(|&f| r.as_acceleration(f)).collect(),
        f_max: F_MAX,
        label: Some((F_MIN, -6.0)),
    };
    let pair = |name, class, a, r, label| Series {
        name,
        class,
        values: freqs.iter().map(|&f| pair_noise(a, r, f)).collect(),
        f_max: F_MAX,
        label: Some(label),
    };
    // Each instrument keeps its colour in every panel: SuperSTAR (today) is the ink reference,
    // MicroSTAR/NGGM, CAI11 and CAI12 take slots 1–3 (the three that separate across all pairs).
    [
        ("a · Accelerometers, ×√2 for the two satellites", vec![
            acc("SuperSTAR (GRACE-FO)", "ink", Accelerometer::SuperStar, F_MAX),
            acc("MicroSTAR (NGGM)", "s1", Accelerometer::MicroStar, F_MAX),
            acc("CAI11", "s2", CAI11, CAI_F_MAX),
            acc("CAI12", "s3", CAI12, CAI_F_MAX),
            acc("Alvarez@350", "s4", Accelerometer::Alvarez350, F_MAX),
            acc("Alvarez@500", "s5", Accelerometer::Alvarez500, F_MAX),
            acc("Dalin", "s6", Accelerometer::Dalin, F_MAX),
        ]),
        ("b · Laser ranging, as range acceleration (2πf)²·ASD", vec![
            lri("LRI GRACE-FO", "ink", Ranging::GraceFo),
            lri("LRI 2030", "ink2", Ranging::Lri2030),
            lri("LRI 2033", "muted", Ranging::Lri2033),
        ]),
        ("c · Combined line-of-sight noise per mission class", vec![
            pair("GRACE-FO class: SuperSTAR + LRI GRACE-FO", "ink", Accelerometer::SuperStar, Ranging::GraceFo, (1e-2, -8.0)),
            pair("NGGM class: MicroSTAR + LRI 2030", "s1", Accelerometer::MicroStar, Ranging::Lri2030, (F_MIN, 16.0)),
            pair("CAI11 class: CAI11 + LRI 2033", "s2", CAI11, Ranging::Lri2033, (F_MIN, -6.0)),
            pair("CAI12 class: CAI12 + LRI 2033", "s3", CAI12, Ranging::Lri2033, (F_MIN, -6.0)),
        ]),
    ]
}

/// Render the figure (SVG) and its table view (CSV).
pub fn figure() -> (String, String) {
    let freqs: Vec<f64> = (0..=400).map(|i| F_MIN * (F_MAX / F_MIN).powf(i as f64 / 400.0)).collect();
    let panels = panels(&freqs);

    let (w, left, right) = (960.0, 92.0, 936.0);
    let plot_h = 230.0;
    let x = |f: f64| left + (right - left) * (f / F_MIN).log10() / (F_MAX / F_MIN).log10();
    let y_in = |top: f64, v: f64| top + plot_h * (Y_MAX / v).log10() / (Y_MAX / Y_MIN).log10();

    let mut s = open(w, 1070.0, "Instrument noise of GRACE-type gravity missions");
    let _ = write!(
        s,
        r#"<text x="{left}" y="32" class="title">Instrument noise of GRACE-type gravity missions</text>
<text x="{left}" y="54" class="subtitle">Amplitude spectral density, m s⁻² Hz⁻¹ᐟ², along the inter-satellite line of sight. Reproduces Zingerle et al. (2024), Fig. 1, from its eqs. 1–8.</text>
"#
    );
    // Top axis: the spherical-harmonic degree a frequency maps to on a 463 km orbit.
    let axis_top = 104.0;
    let _ = writeln!(s, r#"<text x="{left}" y="82" class="tick">Corresponding spherical-harmonic degree, n ≈ f·T (463 km orbit, T = 5620 s)</text>"#);
    let _ = writeln!(s, r#"<line x1="{left}" y1="{axis_top}" x2="{right}" y2="{axis_top}" class="axis"/>"#);
    for n in [1.0, 10.0, 100.0, 500.0] {
        let xx = x(n / T_ORBIT);
        let _ = writeln!(s, r#"<line x1="{xx:.1}" y1="{axis_top}" x2="{xx:.1}" y2="{:.1}" class="axis"/><text x="{xx:.1}" y="{:.1}" class="tick" text-anchor="middle">{n}</text>"#, axis_top - 4.0, axis_top - 8.0);
    }

    let mut top = 120.0;
    for (title, series) in &panels {
        let _ = writeln!(s, r#"<text x="{left}" y="{:.1}" class="panel">{}</text>"#, top + 14.0, esc(title));
        // Legend row: a short line key in the series colour, the name in text ink.
        let mut lx = left;
        for se in series {
            let label = se.name.split(':').next().unwrap_or(se.name);
            let _ = writeln!(s, r#"<line x1="{lx:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" class="line {}"/><text x="{:.1}" y="{:.1}">{}</text>"#,
                top + 32.0, lx + 16.0, top + 32.0, se.class, lx + 22.0, top + 36.0, esc(label));
            lx += 22.0 + 7.0 * label.chars().count() as f64 + 20.0;
        }
        let plot_top = top + 48.0;
        let bottom = plot_top + plot_h;
        // Recessive frame: decade gridlines and the baseline.
        for e in -12..=-9 {
            let yy = y_in(plot_top, 10f64.powi(e));
            let _ = writeln!(s, r#"<line x1="{left}" y1="{yy:.1}" x2="{right}" y2="{yy:.1}" class="grid"/><text x="{:.1}" y="{:.1}" class="tick" text-anchor="end">{}</text>"#, left - 8.0, yy + 4.0, pow10(e));
        }
        for e in -3..=-1 {
            let xx = x(10f64.powi(e));
            let _ = writeln!(s, r#"<line x1="{xx:.1}" y1="{plot_top}" x2="{xx:.1}" y2="{bottom}" class="grid"/>"#);
        }
        let _ = writeln!(s, r#"<line x1="{left}" y1="{bottom}" x2="{right}" y2="{bottom}" class="axis"/>"#);
        for se in series {
            let d: Vec<String> = freqs.iter().zip(&se.values).filter(|(f, _)| **f <= se.f_max * 1.0001)
                .enumerate()
                .map(|(i, (&f, &v))| format!("{}{:.1},{:.1}", if i == 0 { 'M' } else { 'L' }, x(f), y_in(plot_top, v.clamp(Y_MIN, Y_MAX))))
                .collect();
            let _ = writeln!(s, r#"<path d="{}" class="line {}"><title>{}</title></path>"#, d.concat(), se.class, esc(se.name));
        }
        // Direct labels, each placed where its line is clear of the others.
        for se in series {
            let Some((f, dy)) = se.label else { continue };
            let i = freqs.iter().position(|&g| g >= f).unwrap_or(0);
            let label = se.name.split(':').next().unwrap_or(se.name);
            let _ = writeln!(s, r#"<text x="{:.1}" y="{:.1}">{}</text>"#, x(freqs[i]) + 8.0, y_in(plot_top, se.values[i]) + dy, esc(label));
        }
        top = bottom + 34.0;
    }
    let bottom = top - 34.0;
    for e in -3..=-1 {
        let _ = writeln!(s, r#"<text x="{:.1}" y="{:.1}" class="tick" text-anchor="middle">{}</text>"#, x(10f64.powi(e)), bottom + 18.0, pow10(e));
    }
    let _ = write!(
        s,
        r#"<text x="{:.1}" y="{:.1}" class="tick" text-anchor="middle">Frequency (Hz)</text>
<text x="{left}" y="{:.1}" class="subtitle">CAI accelerometers are drawn to 0.03 Hz, as in the paper. The CAI11/CAI12 classes use the CAI white level throughout: the paper's EA–CAI fusion is not published.</text>
</svg>
"#,
        (left + right) / 2.0,
        bottom + 38.0,
        bottom + 62.0
    );

    let mut csv = String::from("frequency_hz,sh_degree");
    for (_, series) in &panels {
        for se in series {
            let _ = write!(csv, ",\"{}\"", se.name);
        }
    }
    csv.push('\n');
    for (i, f) in freqs.iter().enumerate() {
        let _ = write!(csv, "{f:.6e},{:.2}", f * T_ORBIT);
        for (_, series) in &panels {
            for se in series {
                if *f <= se.f_max * 1.0001 {
                    let _ = write!(csv, ",{:.4e}", se.values[i]);
                } else {
                    csv.push(',');
                }
            }
        }
        csv.push('\n');
    }
    (s, csv)
}
