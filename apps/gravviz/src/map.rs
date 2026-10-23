//! A global map of a gridded field: Equal Earth projection, a blue↔red diverging raster whose
//! midpoint is the surface itself, vector coastlines and a colour bar.

use std::fmt::Write;

use crate::svg::{base64, esc, open};

const A1: f64 = 1.340264;
const A2: f64 = -0.081106;
const A3: f64 = 0.000893;
const A4: f64 = 0.003796;
const M: f64 = 0.866_025_403_784_438_6; // √3/2

/// Equal Earth (Šavrič, Patterson & Jenny 2018) on the unit sphere: (lon, lat) radians → (x, y).
pub fn forward(lon: f64, lat: f64) -> (f64, f64) {
    let t = (M * lat.sin()).asin();
    let (t2, t6) = (t * t, t.powi(6));
    let x = lon * t.cos() / (M * (A1 + 3.0 * A2 * t2 + t6 * (7.0 * A3 + 9.0 * A4 * t2)));
    (x, t * (A1 + A2 * t2 + t6 * (A3 + A4 * t2)))
}

/// Inverse of [`forward`]; `None` outside the map outline.
pub fn inverse(x: f64, y: f64) -> Option<(f64, f64)> {
    let mut t = y;
    for _ in 0..12 {
        let (t2, t6) = (t * t, t.powi(6));
        let fy = t * (A1 + A2 * t2 + t6 * (A3 + A4 * t2)) - y;
        t -= fy / (A1 + 3.0 * A2 * t2 + t6 * (7.0 * A3 + 9.0 * A4 * t2));
    }
    let (t2, t6) = (t * t, t.powi(6));
    let s = t.sin() / M;
    let lon = M * x * (A1 + 3.0 * A2 * t2 + t6 * (7.0 * A3 + 9.0 * A4 * t2)) / t.cos();
    (s.abs() <= 1.0 && lon.abs() <= std::f64::consts::PI).then(|| (lon, s.asin()))
}

/// A field on a regular grid of cell centres: `lats` ascending, `lons` ascending over 360°,
/// `values[i * lons.len() + j]`.
pub struct Grid<'a> {
    pub lats: &'a [f64],
    pub lons: &'a [f64],
    pub values: &'a [f64],
}

impl Grid<'_> {
    /// Bilinear sample at (lon, lat) in radians, wrapping in longitude and clamping at the poles.
    fn sample(&self, lon: f64, lat: f64) -> f64 {
        let (nlat, nlon) = (self.lats.len(), self.lons.len());
        let dlat = self.lats[1] - self.lats[0];
        let dlon = self.lons[1] - self.lons[0];
        let fi = ((lat - self.lats[0]) / dlat).clamp(0.0, (nlat - 1) as f64);
        let fj = (lon - self.lons[0]).rem_euclid(std::f64::consts::TAU) / dlon;
        let (i0, j0) = (fi.floor() as usize, fj.floor() as usize);
        let (i1, j1) = ((i0 + 1).min(nlat - 1), (j0 + 1) % nlon);
        let (u, v) = (fi - i0 as f64, fj - j0 as f64);
        let g = |i: usize, j: usize| self.values[i * nlon + j % nlon];
        (1.0 - u) * ((1.0 - v) * g(i0, j0) + v * g(i0, j1)) + u * ((1.0 - v) * g(i1, j0) + v * g(i1, j1))
    }
}

/// Diverging poles: positive values blue, negative red. The midpoint is transparent, so it reads
/// as "nothing" on both the light and the dark surface, and each arm gains contrast with |value|.
const POSITIVE: [u8; 3] = [0x2a, 0x78, 0xd6];
const NEGATIVE: [u8; 3] = [0xe3, 0x49, 0x48];

/// What to put on the map.
pub struct MapText<'a> {
    pub title: &'a str,
    pub subtitle: &'a [&'a str],
    /// Unit shown on the colour bar, e.g. "cm".
    pub unit: &'a str,
    /// Labels for the negative and positive ends of the colour bar.
    pub ends: (&'a str, &'a str),
    pub footnote: &'a str,
}

/// Render `grid` (values in display units) with colour limits ±`limit`, coastlines given as
/// (lon, lat)-degree segments.
pub fn render(grid: &Grid, limit: f64, coast: &[Vec<(f64, f64)>], text: &MapText) -> Result<String, png::EncodingError> {
    let (w, map_w) = (1200.0, 1120.0);
    let (x_max, y_max) = (forward(std::f64::consts::PI, 0.0).0, forward(0.0, std::f64::consts::FRAC_PI_2).1);
    let scale = map_w / (2.0 * x_max);
    let map_h = (2.0 * y_max * scale).ceil();
    let (map_x, map_y) = (40.0, 74.0 + 18.0 * text.subtitle.len() as f64);
    let to_px = |lon: f64, lat: f64| {
        let (x, y) = forward(lon.to_radians(), lat.to_radians());
        (map_x + map_w / 2.0 + x * scale, map_y + map_h / 2.0 - y * scale)
    };

    // Raster at the map's own pixel size.
    let (pw, ph) = (map_w as u32, map_h as u32);
    let mut rgba = vec![0u8; (pw * ph * 4) as usize];
    for py in 0..ph {
        for px in 0..pw {
            let x = (f64::from(px) + 0.5 - map_w / 2.0) / scale;
            let y = (map_h / 2.0 - f64::from(py) - 0.5) / scale;
            if let Some((lon, lat)) = inverse(x, y) {
                let t = (grid.sample(lon, lat) / limit).clamp(-1.0, 1.0);
                let rgb = if t >= 0.0 { POSITIVE } else { NEGATIVE };
                let k = ((py * pw + px) * 4) as usize;
                rgba[k..k + 3].copy_from_slice(&rgb);
                rgba[k + 3] = (t.abs() * 255.0).round() as u8;
            }
        }
    }
    let mut png_bytes = Vec::new();
    {
        let mut enc = png::Encoder::new(&mut png_bytes, pw, ph);
        enc.set_color(png::ColorType::Rgba);
        enc.set_depth(png::BitDepth::Eight);
        enc.write_header()?.write_image_data(&rgba)?;
    }

    let colorbar_y = map_y + map_h + 34.0;
    let mut s = open(w, colorbar_y + 96.0, text.title);
    let _ = writeln!(s, r#"<text x="{map_x}" y="32" class="title">{}</text>"#, esc(text.title));
    for (i, line) in text.subtitle.iter().enumerate() {
        let _ = writeln!(s, r#"<text x="{map_x}" y="{}" class="subtitle">{}</text>"#, 54.0 + 18.0 * i as f64, esc(line));
    }
    let _ = writeln!(s, r#"<image x="{map_x}" y="{map_y}" width="{map_w}" height="{map_h}" href="data:image/png;base64,{}"/>"#, base64(&png_bytes));

    // Coastlines, split where a segment crosses the antimeridian.
    let mut d = String::new();
    for seg in coast {
        let mut prev: Option<f64> = None;
        for &(lon, lat) in seg {
            let (px, py) = to_px(lon, lat);
            let cmd = if prev.is_some_and(|p| (p - lon).abs() < 180.0) { 'L' } else { 'M' };
            let _ = write!(d, "{cmd}{px:.1},{py:.1}");
            prev = Some(lon);
        }
    }
    let _ = writeln!(s, r#"<path d="{d}" fill="none" stroke="var(--ink2)" stroke-width="0.7" stroke-opacity="0.8"/>"#);
    // Map outline: the ±180° meridians.
    let mut outline = String::new();
    for (k, lon) in [-180.0, 180.0].into_iter().enumerate() {
        for i in 0..=90 {
            let lat = if k == 0 { -90.0 + 2.0 * i as f64 } else { 90.0 - 2.0 * i as f64 };
            let (px, py) = to_px(lon, lat);
            let _ = write!(outline, "{}{px:.1},{py:.1}", if k == 0 && i == 0 { 'M' } else { 'L' });
        }
    }
    let _ = writeln!(s, r#"<path d="{outline}Z" fill="none" class="axis"/>"#);

    // Colour bar: the same transparent-midpoint ramp as the raster.
    let (cb_w, cb_x) = (420.0, w / 2.0 - 210.0);
    let hex = |c: [u8; 3]| format!("#{:02x}{:02x}{:02x}", c[0], c[1], c[2]);
    let _ = write!(
        s,
        r#"<defs><linearGradient id="ramp">
<stop offset="0" stop-color="{neg}" stop-opacity="1"/><stop offset="0.5" stop-color="{neg}" stop-opacity="0"/>
<stop offset="0.5" stop-color="{pos}" stop-opacity="0"/><stop offset="1" stop-color="{pos}" stop-opacity="1"/>
</linearGradient></defs>
<rect x="{cb_x}" y="{colorbar_y}" width="{cb_w}" height="12" rx="2" fill="url(#ramp)" stroke="var(--axis)"/>
"#,
        neg = hex(NEGATIVE),
        pos = hex(POSITIVE)
    );
    for k in -2i32..=2 {
        let v = limit * f64::from(k) / 2.0;
        let xx = cb_x + cb_w * (0.5 + f64::from(k) / 4.0);
        let sign = if k > 0 { "+" } else if k < 0 { "−" } else { "" };
        let edge = if k.abs() == 2 { if k < 0 { "≤ " } else { "≥ " } } else { "" };
        let _ = writeln!(s, r#"<line x1="{xx:.1}" y1="{}" x2="{xx:.1}" y2="{}" class="axis"/><text x="{xx:.1}" y="{}" class="tick" text-anchor="middle">{edge}{sign}{} {}</text>"#,
            colorbar_y + 12.0, colorbar_y + 17.0, colorbar_y + 30.0, v.abs(), text.unit);
    }
    let _ = write!(
        s,
        r#"<text x="{:.1}" y="{:.1}" text-anchor="end">{}</text>
<text x="{:.1}" y="{:.1}">{}</text>
<text x="{map_x}" y="{:.1}" class="subtitle">{}</text>
</svg>
"#,
        cb_x - 10.0, colorbar_y + 10.5, esc(text.ends.0),
        cb_x + cb_w + 10.0, colorbar_y + 10.5, esc(text.ends.1),
        colorbar_y + 66.0, esc(text.footnote)
    );
    Ok(s)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn equal_earth_round_trips_and_has_the_published_extent() {
        let (x, y) = (forward(std::f64::consts::PI, 0.0).0, forward(0.0, std::f64::consts::FRAC_PI_2).1);
        assert!((x - 2.7066).abs() < 1e-4 && (y - 1.3173).abs() < 1e-4, "{x} {y}");
        for (lon, lat) in [(0.3, 0.2), (-2.9, 1.3), (1.0, -1.5), (3.1, 0.0)] {
            let (px, py) = forward(lon, lat);
            let (lo, la) = inverse(px, py).unwrap();
            assert!((lo - lon).abs() < 1e-9 && (la - lat).abs() < 1e-9);
        }
        assert!(inverse(2.8, 0.0).is_none());
    }
}
