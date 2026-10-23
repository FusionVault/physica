//! A coefficient "triangle": one cell per (degree, order), sine coefficients to the left of the
//! zonal column and cosine to the right, as in Zingerle et al. (2024), Fig. 3a/c. Here it shows
//! the calibrated errors of a real monthly field, as water height per coefficient.

use std::fmt::Write;

use geopotential::{Coefficients, Functional};

use crate::svg::{esc, open, png_data_uri, pow10};

/// Sequential hue (blue): magnitude is its opacity over the surface, so small values recede
/// toward the surface in light and dark mode alike.
const HUE: [u8; 3] = [0x2a, 0x78, 0xd6];

pub fn figure(sigma: &Coefficients, month: &str) -> Result<String, Box<dyn std::error::Error>> {
    let nmax = sigma.max_degree();
    let love = gravimetry::love_numbers(nmax).ok_or("degree above the Love-number table")?;
    let ewh = Functional::WaterHeight { love: &love };
    let mm = |n: usize, v: f64| v * ewh.factor(n, sigma.gm, sigma.radius, sigma.radius) * 1e3;
    // One value per cell: x = nmax + m for C̄nm, nmax − m for S̄nm.
    let (cols, rows) = (2 * nmax + 1, nmax + 1);
    let mut cells = vec![f64::NAN; cols * rows];
    for n in 2..=nmax {
        for m in 0..=n {
            cells[n * cols + nmax + m] = mm(n, sigma.c(n, m));
            if m > 0 {
                cells[n * cols + nmax - m] = mm(n, sigma.s(n, m));
            }
        }
    }
    let finite: Vec<f64> = cells.iter().copied().filter(|v| v.is_finite() && *v > 0.0).collect();
    let (lo, hi) = finite.iter().fold((f64::MAX, f64::MIN), |(a, b), &v| (a.min(v), b.max(v)));
    let (e_lo, e_hi) = (lo.log10().floor() as i32, hi.log10().ceil() as i32);
    let alpha = |v: f64| 0.08 + 0.92 * ((v.log10() - f64::from(e_lo)) / f64::from(e_hi - e_lo)).clamp(0.0, 1.0);

    let px = 5u32;
    let (pw, ph) = (cols as u32 * px, rows as u32 * px);
    let mut rgba = vec![0u8; (pw * ph * 4) as usize];
    for (k, &v) in cells.iter().enumerate() {
        if !(v.is_finite() && v > 0.0) {
            continue;
        }
        let a = (alpha(v) * 255.0).round() as u8;
        let (cx, cy) = ((k % cols) as u32 * px, (k / cols) as u32 * px);
        for y in cy..cy + px {
            for x in cx..cx + px {
                let i = ((y * pw + x) * 4) as usize;
                rgba[i..i + 3].copy_from_slice(&HUE);
                rgba[i + 3] = a;
            }
        }
    }

    let (left, top) = (64.0, 132.0);
    let (w, h) = (left + f64::from(pw) + 40.0, top + f64::from(ph) + 120.0);
    let title = format!("Calibrated errors of the COST-G {month} field, per coefficient");
    let mut s = open(w, h, &title);
    let _ = write!(
        s,
        r#"<text x="{left}" y="32" class="title">{t}</text>
<text x="{left}" y="54" class="subtitle">Equivalent water height of each coefficient's error, mm (log scale). Sine coefficients left of the zonal column, cosine right.</text>
<text x="{left}" y="72" class="subtitle">Real data in the layout of Zingerle et al. (2024), Fig. 3a/c, which shows simulated NGGM and CAI11 errors.</text>
<image x="{left}" y="{top}" width="{pw}" height="{ph}" href="{img}" style="image-rendering: pixelated"/>
"#,
        t = esc(&title),
        img = png_data_uri(&rgba, pw, ph)?
    );
    let cell = f64::from(px);
    for m in [-(nmax as i64), -60, -30, 0, 30, 60, nmax as i64] {
        let xx = left + (nmax as f64 + m as f64 + 0.5) * cell;
        let _ = writeln!(s, r#"<text x="{xx:.1}" y="{}" class="tick" text-anchor="middle">{}</text>"#, top - 8.0, m.abs());
    }
    let _ = writeln!(s, r#"<text x="{left}" y="{}" class="tick">S̄nm ← order → C̄nm</text>"#, top - 30.0);
    for n in (0..=nmax).step_by(10) {
        let yy = top + (n as f64 + 0.5) * cell;
        let _ = writeln!(s, r#"<text x="{}" y="{:.1}" class="tick" text-anchor="end">{n}</text>"#, left - 8.0, yy + 4.0);
    }
    let _ = writeln!(s, r#"<text x="16" y="{}" class="tick" transform="rotate(-90 16 {})" text-anchor="middle">Degree</text>"#, top + f64::from(ph) / 2.0, top + f64::from(ph) / 2.0);

    // Colour bar: the same opacity ramp, with decade ticks.
    let (cb_y, cb_w, cb_x) = (top + f64::from(ph) + 36.0, 360.0, left + f64::from(pw) / 2.0 - 180.0);
    let hex = format!("#{:02x}{:02x}{:02x}", HUE[0], HUE[1], HUE[2]);
    let _ = writeln!(s, r#"<defs><linearGradient id="seq"><stop offset="0" stop-color="{hex}" stop-opacity="0.08"/><stop offset="1" stop-color="{hex}" stop-opacity="1"/></linearGradient></defs><rect x="{cb_x}" y="{cb_y}" width="{cb_w}" height="12" rx="2" fill="url(#seq)" stroke="var(--axis)"/>"#);
    for e in e_lo..=e_hi {
        let xx = cb_x + cb_w * f64::from(e - e_lo) / f64::from(e_hi - e_lo);
        let _ = writeln!(s, r#"<line x1="{xx:.1}" y1="{}" x2="{xx:.1}" y2="{}" class="axis"/><text x="{xx:.1}" y="{}" class="tick" text-anchor="middle">{} mm</text>"#, cb_y + 12.0, cb_y + 17.0, cb_y + 30.0, pow10(e));
    }
    let _ = write!(
        s,
        r#"<text x="{left}" y="{}" class="subtitle">Unsmoothed: errors grow with degree and toward high orders, the pattern that becomes north–south stripes on maps. Data: COST-G (CC BY 4.0).</text>
</svg>
"#,
        cb_y + 62.0
    );
    Ok(s)
}
