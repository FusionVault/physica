//! T1 feasibility: would the 2022 Hunga Tonga Lamb wave show in GRACE-FO's laser ranging?
//!
//! The eruption's Lamb wave carried a surface-pressure pulse around the world at 318.2 m/s
//! ([Wright et al. 2022](https://doi.org/10.1038/s41586-022-05012-5)); barometers far from the
//! volcano recorded 1.3–2 hPa (13–20 kg/m² of air; [arXiv 2303.04362](https://arxiv.org/abs/2303.04362),
//! [Frontiers 2022](https://doi.org/10.3389/feart.2022.931545)). GRACE-FO crosses such a front in
//! about a minute at 7 km/s over the ground, so the front barely moves while the pair flies through
//! it; [`template_moving`] includes its motion, which matters at shallow crossings. Locally the front is a straight strip of surface mass `σ(x)` (a Gaussian across
//! the front), which pulls each satellite horizontally with
//!
//! `g_x(x₀) = −2G ∫ σ(x)(x₀ − x)/((x₀ − x)² + h²) dx`.
//!
//! When the track crosses the front at angle φ, the range acceleration is the difference of the
//! along-track pulls, `ρ̈ = sin φ·[g_x(x_lead) − g_x(x_trail)]`, with `x = u·sin φ` along the
//! track. Its matched-filter SNR against the pair's instrument noise, `ρ² = 4∫|ρ̈(f)|²/S_n df`,
//! says whether one crossing can be seen. Everything is per hPa: SNR is linear in the pulse.

use std::fmt::Write;

use gravimetry::atmosphere::{forced_ocean_response, pressure_to_mass};
use gravimetry::missions::{Accelerometer, Ranging, pair_noise};
use gravimetry::noise::fft;

use crate::Result;
use crate::svg::{esc, open, pow10};

const G: f64 = 6.674_30e-11;
const GM: f64 = 3.986_004_418e14;
const R_EARTH: f64 = 6_371e3;
const LAMB_SPEED: f64 = 318.2;
const ALTITUDE: f64 = 490e3;
const SEPARATION: f64 = 220e3;
const DT: f64 = 2.0; // model sampling, s
const N: usize = 2048;
const WIDTHS_MIN: [f64; 4] = [1.25, 2.5, 5.0, 10.0]; // Gaussian σ of the pulse at a fixed place
const ANGLES_DEG: [f64; 4] = [90.0, 45.0, 20.0, 10.0];
const DETECT: f64 = 5.0;

// shortcut: flat Earth in ground coordinates (horizontal distances at 490 km are ~7% larger, and
// the line of sight tilts ~1° from the track); fine for feasibility, use the real orbit for T1.
fn ground_speed() -> f64 {
    (GM / (R_EARTH + ALTITUDE)).sqrt() * R_EARTH / (R_EARTH + ALTITUDE)
}

fn ground_separation() -> f64 {
    SEPARATION * R_EARTH / (R_EARTH + ALTITUDE)
}

/// Horizontal gravity (m/s², positive toward +x) at height `h` and offset `x0` from the centre of
/// an infinite strip of surface mass `sigma0·exp(−x²/2s²)` kg/m², by the trapezoid rule over ±8s.
fn strip_gx(x0: f64, sigma0: f64, s: f64, h: f64) -> f64 {
    const STEPS: usize = 400;
    let dx = 16.0 * s / STEPS as f64;
    let sum: f64 = (0..=STEPS).map(|i| {
        let x = -8.0 * s + i as f64 * dx;
        let r = x0 - x;
        let w = if i == 0 || i == STEPS { 0.5 } else { 1.0 };
        w * (-x * x / (2.0 * s * s)).exp() * r / (r * r + h * h)
    }).sum();
    -2.0 * G * sigma0 * sum * dx
}

/// Range acceleration (m/s²) for a 1 hPa pulse of time width `width` (s) crossed at `angle`
/// (rad) by a frozen front, sampled every `DT` with the crossing in the middle.
pub(crate) fn template(width: f64, angle: f64) -> Vec<f64> {
    template_moving(width, angle, ground_speed() * angle.sin())
}

/// As [`template`], for a pair whose midpoint closes on the front at `normal_speed` (m/s, its
/// speed across the front relative to the front's own motion). The pair's separation still
/// projects onto the front normal as `d·sin φ`.
pub(crate) fn template_moving(width: f64, angle: f64, normal_speed: f64) -> Vec<f64> {
    let (sigma0, s) = (pressure_to_mass(100.0), LAMB_SPEED * width);
    let (sin, d) = (angle.sin(), ground_separation());
    (0..N).map(|i| {
        let x = (i as f64 - (N / 2) as f64) * DT * normal_speed;
        sin * (strip_gx(x + d / 2.0 * sin, sigma0, s, ALTITUDE) - strip_gx(x - d / 2.0 * sin, sigma0, s, ALTITUDE))
    }).collect()
}

/// One-sided Fourier amplitudes `|ρ̈(f_k)|` (m s⁻² Hz⁻¹) of a series sampled every `DT`, k = 1…N/2−1.
pub(crate) fn spectrum(signal: &[f64]) -> Vec<(f64, f64)> {
    let (mut re, mut im) = (signal.to_vec(), vec![0.0; signal.len()]);
    fft(&mut re, &mut im, false);
    let df = 1.0 / (signal.len() as f64 * DT);
    (1..signal.len() / 2).map(|k| (k as f64 * df, DT * re[k].hypot(im[k]))).collect()
}

/// Matched-filter SNR against a one-sided noise ASD, and the band holding 90% of SNR².
fn snr(signal: &[f64], asd: impl Fn(f64) -> f64) -> (f64, (f64, f64)) {
    let spec = spectrum(signal);
    let df = spec[0].0;
    let terms: Vec<f64> = spec.iter().map(|&(f, a)| 4.0 * a * a / asd(f).powi(2) * df).collect();
    let total: f64 = terms.iter().sum();
    let mut acc = 0.0;
    let (mut lo, mut hi) = (spec[0].0, spec[spec.len() - 1].0);
    for (&(f, _), t) in spec.iter().zip(&terms) {
        let before = acc;
        acc += t;
        if before < 0.05 * total && acc >= 0.05 * total { lo = f; }
        if before < 0.95 * total && acc >= 0.95 * total { hi = f; }
    }
    (total.sqrt(), (lo, hi))
}

fn noise(f: f64) -> f64 {
    pair_noise(Accelerometer::SuperStar, Ranging::GraceFo, f)
}

pub fn run() -> Result<(String, String)> {
    let mut r = String::new();
    let ocean = 1.0 + forced_ocean_response(LAMB_SPEED, 4000.0);
    let antipode_h = std::f64::consts::PI * R_EARTH / LAMB_SPEED / 3600.0;
    writeln!(r, "T1 feasibility: the Hunga Tonga Lamb wave in GRACE-FO laser ranging")?;
    writeln!(r, "Pulse: Gaussian across the front, 1 hPa = {:.1} kg/m², moving at {LAMB_SPEED} m/s.", pressure_to_mass(100.0))?;
    writeln!(r, "Pair: {:.0} km altitude, {:.0} km apart, {:.2} km/s over the ground; noise SuperSTAR + GRACE-FO LRI (Zingerle et al. 2024).", ALTITUDE / 1e3, SEPARATION / 1e3, ground_speed() / 1e3)?;
    writeln!(r, "The front reaches Tonga's antipode after {antipode_h:.1} h and is back at Tonga after {:.1} h.\n", 2.0 * antipode_h)?;
    writeln!(r, "Peak range acceleration (nm/s²) and matched-filter SNR of one crossing, per hPa of pulse:")?;
    write!(r, "{:>22}", "pulse σ (min) →")?;
    for w in WIDTHS_MIN { write!(r, "{w:>16}")?; }
    writeln!(r)?;
    let mut bands = Vec::new();
    for a in ANGLES_DEG {
        write!(r, "{:>22}", format!("crossing at {a}°"))?;
        for w in WIDTHS_MIN {
            let t = template(w * 60.0, a.to_radians());
            let peak = t.iter().fold(0.0f64, |m, v| m.max(v.abs()));
            let (rho, band) = snr(&t, noise);
            if a == 90.0 { bands.push(band); }
            write!(r, "{:>16}", format!("{:.2} / {rho:.1}", peak * 1e9))?;
        }
        writeln!(r)?;
    }
    writeln!(r, "\n90% of SNR² lies in (mHz, perpendicular crossing):")?;
    for (w, (lo, hi)) in WIDTHS_MIN.iter().zip(&bands) {
        writeln!(r, "  σ = {w:>5} min: {:.1}–{:.1}", lo * 1e3, hi * 1e3)?;
    }
    writeln!(r, "\nReading it:")?;
    writeln!(r, "- Far from Tonga barometers recorded 1.3–2 hPa: multiply the SNRs by that.")?;
    writeln!(r, "- Over 4 km of ocean the sea rises under the pulse (it is faster than ocean long waves,")?;
    writeln!(r, "  Proudman 1929), adding its mass: ×{ocean:.2} (×{:.2} over 5 km).", 1.0 + forced_ocean_response(LAMB_SPEED, 5000.0))?;
    writeln!(r, "- SNR {DETECT} per crossing needs the real noise to be no worse than SNR/{DETECT} times nominal.")?;
    writeln!(r, "  GRACE-D's accelerometer data are transplanted from GRACE-C, and in 1–8 mHz (degrees ~7–43)")?;
    writeln!(r, "  weather and ocean signals are removed only by AOD1B, whose residuals there are unknown. The")?;
    writeln!(r, "  LRI itself resolves ~0.1 nm/s² (Ghobadi-Far et al. 2020); the margin is what must survive.")?;
    writeln!(r, "- Many crossings: each orbit cuts the expanding front about twice, ~30 a day.")?;
    Ok((r, figure()))
}

fn figure() -> String {
    let (w, left, right) = (960.0, 92.0, 936.0);
    let mut s = open(w, 700.0, "The Tonga Lamb wave as GRACE-FO would see it");
    let _ = write!(s, r#"<text x="{left}" y="32" class="title">The Tonga Lamb wave as GRACE-FO would see it</text>
<text x="{left}" y="54" class="subtitle">A 1 hPa pressure pulse crossed at right angles, for four pulse widths (σ at a fixed place). Flat-Earth forward model, 490 km, 220 km apart.</text>
"#);
    let series: Vec<(f64, &str, Vec<f64>)> = WIDTHS_MIN.iter().zip(["s1", "s2", "s3", "s4"])
        .map(|(&m, c)| (m, c, template(m * 60.0, std::f64::consts::FRAC_PI_2))).collect();

    // Panel 1: range acceleration against time from the crossing.
    let (top, h, t_max, y_max) = (110.0, 220.0, 400.0, 0.5);
    let x = |t: f64| left + (right - left) * (t + t_max) / (2.0 * t_max);
    let y = |v: f64| top + h / 2.0 - h / 2.0 * v / y_max;
    let _ = writeln!(s, r#"<text x="{left}" y="{:.1}" class="panel">Range acceleration, nm/s², against seconds from the crossing</text>"#, top - 24.0);
    legend(&mut s, left, top - 6.0, &series);
    for v in [-0.5, -0.25, 0.0, 0.25, 0.5] {
        let yy = y(v);
        let _ = writeln!(s, r#"<line x1="{left}" y1="{yy:.1}" x2="{right}" y2="{yy:.1}" class="{}"/><text x="{:.1}" y="{:.1}" class="tick" text-anchor="end">{v}</text>"#,
            if v == 0.0 { "axis" } else { "grid" }, left - 8.0, yy + 4.0);
    }
    for t in [-300.0, -150.0, 0.0, 150.0, 300.0] {
        let _ = writeln!(s, r#"<text x="{:.1}" y="{:.1}" class="tick" text-anchor="middle">{t}</text>"#, x(t), top + h + 18.0);
    }
    for (m, c, t) in &series {
        let d: String = t.iter().enumerate().map(|(i, v)| ((i as f64 - (N / 2) as f64) * DT, v))
            .filter(|(tt, _)| tt.abs() <= t_max)
            .enumerate().map(|(j, (tt, v))| format!("{}{:.1},{:.1}", if j == 0 { 'M' } else { 'L' }, x(tt), y(v * 1e9))).collect();
        let _ = writeln!(s, r#"<path d="{d}" class="line {c}"><title>σ = {m} min</title></path>"#);
    }

    // Panel 2: signal spectrum 2√f|ρ̈(f)| against the pair's noise ASD.
    let (top, h) = (430.0, 220.0);
    let (f_min, f_max, a_min, a_max) = (1e-4, 0.25, 1e-13, 1e-8);
    let x = |f: f64| left + (right - left) * (f / f_min).log10() / (f_max / f_min).log10();
    let y = |v: f64| top + h * (a_max / v.clamp(a_min, a_max)).log10() / (a_max / a_min).log10();
    let _ = writeln!(s, r#"<text x="{left}" y="{:.1}" class="panel">{}</text>"#, top - 24.0,
        esc("Signal 2√f·|ρ̈(f)| against the pair noise, m s⁻² Hz⁻¹ᐟ²: the gap between them is the SNR"));
    for e in -13..=-8 {
        let yy = y(10f64.powi(e));
        let _ = writeln!(s, r#"<line x1="{left}" y1="{yy:.1}" x2="{right}" y2="{yy:.1}" class="grid"/><text x="{:.1}" y="{:.1}" class="tick" text-anchor="end">{}</text>"#, left - 8.0, yy + 4.0, pow10(e));
    }
    for (f, label) in [(1e-4, "0.1 mHz"), (1e-3, "1 mHz"), (1e-2, "10 mHz"), (1e-1, "0.1 Hz")] {
        let _ = writeln!(s, r#"<text x="{:.1}" y="{:.1}" class="tick" text-anchor="middle">{label}</text>"#, x(f), top + h + 18.0);
    }
    let path = |pts: &mut dyn Iterator<Item = (f64, f64)>| -> String {
        pts.enumerate().map(|(j, (f, v))| format!("{}{:.1},{:.1}", if j == 0 { 'M' } else { 'L' }, x(f), y(v))).collect()
    };
    for (m, c, t) in &series {
        let d = path(&mut spectrum(t).into_iter().map(|(f, a)| (f, 2.0 * f.sqrt() * a)).take_while(|&(_, v)| v >= a_min));
        let _ = writeln!(s, r#"<path d="{d}" class="line {c}"><title>σ = {m} min</title></path>"#);
    }
    let d = path(&mut (0..=200).map(|i| { let f = f_min * (f_max / f_min).powf(i as f64 / 200.0); (f, noise(f)) }));
    let _ = writeln!(s, r#"<path d="{d}" class="line ink"><title>SuperSTAR + LRI pair noise</title></path>"#);
    let _ = writeln!(s, r#"<text x="{:.1}" y="{:.1}">pair noise</text>"#, x(0.12) - 40.0, y(noise(0.12)) - 10.0);
    s.push_str("</svg>\n");
    s
}

fn legend(s: &mut String, left: f64, y: f64, series: &[(f64, &str, Vec<f64>)]) {
    let mut lx = left;
    for (m, c, _) in series {
        let label = format!("σ = {m} min");
        let _ = writeln!(s, r#"<line x1="{lx:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" class="line {c}"/><text x="{:.1}" y="{y:.1}">{label}</text>"#,
            y - 4.0, lx + 16.0, y - 4.0, lx + 22.0);
        lx += 22.0 + 6.5 * label.chars().count() as f64 + 14.0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_narrow_strip_pulls_like_a_line_mass() {
        // s ≪ h: g_x → −2Gλ[K + (s²/2)K″] with λ = σ₀√(2π)s, K = x₀/(x₀² + h²) and
        // K″ = 2x₀(x₀² − 3h²)/(x₀² + h²)³ (the Gaussian's width correction).
        let (sigma0, s, h) = (20.0, 1e3, 490e3);
        let lambda = sigma0 * (2.0 * std::f64::consts::PI).sqrt() * s;
        for x0 in [50e3, 300e3, 1e6] {
            let q = x0 * x0 + h * h;
            let k = x0 / q + s * s / 2.0 * 2.0 * x0 * (x0 * x0 - 3.0 * h * h) / q.powi(3);
            let err = strip_gx(x0, sigma0, s, h) / (-2.0 * G * lambda * k) - 1.0;
            assert!(err.abs() < 1e-8, "{x0}: {err:e}");
        }
        assert!(strip_gx(0.0, sigma0, 50e3, h).abs() < 1e-20);
    }

    #[test]
    fn oblique_and_moving_crossings() {
        // A frozen front crossed at φ is the moving-front template at normal speed v·sin φ.
        let phi = 30f64.to_radians();
        let a = template(300.0, phi);
        let b = template_moving(300.0, phi, ground_speed() * phi.sin());
        assert!(a.iter().zip(&b).all(|(x, y)| (x - y).abs() < 1e-24));
        // The signal is even in time (direction of crossing does not matter), and a slower
        // closing speed stretches it: the same peak, wider.
        assert!((a[N / 2 + 40] - a[N / 2 - 40]).abs() < 1e-20);
        let slow = template_moving(300.0, phi, ground_speed() * phi.sin() / 2.0);
        assert!((slow[N / 2] / a[N / 2] - 1.0).abs() < 1e-9);
        assert!((slow[N / 2 + 40] - a[N / 2 + 20]).abs() < 1e-3 * a[N / 2].abs());
    }

    #[test]
    fn snr_of_a_pulse_in_white_noise() {
        // White one-sided PSD S₀: ρ² = 2∫ρ̈² dt/S₀. For the zero-mean pulse a(t/σ)e^(−t²/2σ²)
        // (like the template: the DC bin is left out) that is a²σ√π/S₀.
        let (a, sigma, s0) = (1e-9, 60.0, 1e-20);
        let pulse: Vec<f64> = (0..N).map(|i| { let t = (i as f64 - (N / 2) as f64) * DT; a * t / sigma * (-t * t / (2.0 * sigma * sigma)).exp() }).collect();
        let expect = (a * a * sigma * std::f64::consts::PI.sqrt() / s0).sqrt();
        let (rho, _) = snr(&pulse, |_| s0.sqrt());
        assert!((rho / expect - 1.0).abs() < 1e-3, "{rho} vs {expect}");
    }
}
