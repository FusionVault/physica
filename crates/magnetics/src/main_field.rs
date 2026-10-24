//! Earth's main magnetic field from spherical-harmonic Gauss coefficients: the World Magnetic
//! Model (WMM2025, built in) or any model in its `.COF` format, evaluated with `geopotential`'s
//! gravity engine.
//!
//! The field is `B = −∇V` with `V = a Σ (a/r)^(n+1) Σ (gₙᵐ cos mλ + hₙᵐ sin mλ) Pₙᵐ(sin φ′)`,
//! `a = 6371.2 km` and Schmidt semi-normalised `Pₙᵐ`. Since Schmidt and full (4π) normalisation
//! differ by `√(2n+1)` and neither uses the Condon–Shortley phase, the model is a gravity field
//! with `GM = a²`, `R = a` and `C̄ₙₘ = gₙᵐ/√(2n+1)`: the same code evaluates both.
//!
//! ```
//! use magnetics::main_field::wmm2025;
//!
//! // Declination in Boulder, Colorado (40.0°N, 105.3°W, 1.6 km), mid-2026.
//! let e = wmm2025().elements(2026.5, 40f64.to_radians(), (-105.3f64).to_radians(), 1600.0);
//! println!("declination {:.2}°, inclination {:.2}°, {:.0} nT", e.declination.to_degrees(), e.inclination.to_degrees(), e.f);
//! assert!((7.0..9.0).contains(&e.declination.to_degrees()));
//! ```

use std::fmt;

use geopotential::{Coefficients, Field, geocentric};

/// The geomagnetic reference radius `a` of WMM and IGRF, m.
pub const REFERENCE_RADIUS: f64 = 6_371_200.0;

/// A main-field model: Gauss coefficients at an epoch and their linear secular variation.
#[derive(Clone, Debug, PartialEq)]
pub struct Model {
    /// The model's name, e.g. `WMM-2025`.
    pub name: String,
    /// The epoch of the coefficients, decimal years.
    pub epoch: f64,
    /// The field at the epoch, nT, as a gravity-engine [`Coefficients`] (see the module docs).
    pub main: Coefficients,
    /// Its rate of change, nT/yr, in the same form.
    pub secular: Coefficients,
}

/// The field at one point, in the local geodetic frame: north (X), east (Y) and down (Z), nT, and
/// the angles navigators use (radians), with the rates of each per year.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Elements {
    /// North component X, nT.
    pub x: f64,
    /// East component Y, nT.
    pub y: f64,
    /// Down component Z, nT.
    pub z: f64,
    /// Horizontal intensity H, nT.
    pub h: f64,
    /// Total intensity F, nT.
    pub f: f64,
    /// Inclination (dip) I, positive down, radians.
    pub inclination: f64,
    /// Declination D, east of true north, radians.
    pub declination: f64,
    /// Rates of X, Y, Z, H and F, nT/yr.
    pub rates: [f64; 5],
    /// Rate of inclination, radians/yr.
    pub inclination_rate: f64,
    /// Rate of declination, radians/yr.
    pub declination_rate: f64,
}

/// A main-field model evaluated at one date: a [`Field`] whose potential is in nT m and whose
/// vector is `B` in nT (north, east, up, geocentric).
#[derive(Clone, Debug, PartialEq)]
pub struct MainField(pub Coefficients);

impl Field for MainField {
    fn potential(&self, r: f64, lat: f64, lon: f64) -> f64 {
        self.0.potential(r, lat, lon)
    }

    fn vector(&self, r: f64, lat: f64, lon: f64) -> [f64; 3] {
        self.0.vector(r, lat, lon).map(|v| -v)
    }
}

impl Model {
    /// The field at decimal year `t`, extrapolated linearly from the epoch.
    pub fn at(&self, t: f64) -> MainField {
        let dt = t - self.epoch;
        let mut c = self.main.clone();
        for n in 1..=c.max_degree() {
            for m in 0..=n {
                c.set(n, m, c.c(n, m) + dt * self.secular.c(n, m), c.s(n, m) + dt * self.secular.s(n, m));
            }
        }
        MainField(c)
    }

    /// The field elements at decimal year `t`, geodetic latitude `lat` and longitude `lon`
    /// (radians) and height `height` (m) above the WGS84 ellipsoid: the WMM report's steps 1–5.
    pub fn elements(&self, t: f64, lat: f64, lon: f64, height: f64) -> Elements {
        let (r, lat_c) = geocentric(lat, height);
        // From geocentric (north, east, up) to geodetic (north, east, down).
        let psi = lat_c - lat;
        let frame = |[n, e, u]: [f64; 3]| [n * psi.cos() + u * psi.sin(), e, n * psi.sin() - u * psi.cos()];
        let [x, y, z] = frame(self.at(t).vector(r, lat_c, lon));
        let [xd, yd, zd] = frame(MainField(self.secular.clone()).vector(r, lat_c, lon));
        let h = x.hypot(y);
        let f = h.hypot(z);
        let hd = (x * xd + y * yd) / h;
        let fd = (x * xd + y * yd + z * zd) / f;
        Elements {
            x, y, z, h, f,
            inclination: z.atan2(h),
            declination: y.atan2(x),
            rates: [xd, yd, zd, hd, fd],
            inclination_rate: (h * zd - z * hd) / (f * f),
            declination_rate: (x * yd - y * xd) / (h * h),
        }
    }

    /// The dipole terms `(g₁₀, g₁₁, h₁₁)` at decimal year `t`, nT (Schmidt, as published).
    fn dipole(&self, t: f64) -> (f64, f64, f64) {
        let c = self.at(t).0;
        let k = 3f64.sqrt();
        (k * c.c(1, 0), k * c.c(1, 1), k * c.s(1, 1))
    }

    /// The dipole moment `4πa³/μ₀ · √(g₁₀² + g₁₁² + h₁₁²)` at decimal year `t`, A m².
    pub fn dipole_moment(&self, t: f64) -> f64 {
        let (g10, g11, h11) = self.dipole(t);
        4.0 * core::f64::consts::PI * REFERENCE_RADIUS.powi(3) / metrology::MU_0 * g10.hypot(g11).hypot(h11) * 1e-9
    }

    /// The north geomagnetic pole, where the dipole's axis meets the sphere: geocentric latitude
    /// and longitude (radians) at decimal year `t`. Not the dip pole, where the field is vertical.
    pub fn geomagnetic_pole(&self, t: f64) -> (f64, f64) {
        let (g10, g11, h11) = self.dipole(t);
        (core::f64::consts::FRAC_PI_2 - (-g10 / g10.hypot(g11).hypot(h11)).acos(), (-h11).atan2(-g11))
    }
}

/// Why a `.COF` file could not be read.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CofError {
    /// The 1-based line number.
    pub line: usize,
    /// What was wrong.
    pub reason: &'static str,
}

impl fmt::Display for CofError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "line {}: {}", self.line, self.reason)
    }
}

impl std::error::Error for CofError {}

/// Read a model in NOAA's `.COF` format: a header line `epoch name date`, then one line per
/// `n m g h ġ ḣ` (nT and nT/yr, Schmidt semi-normalised), ending at a line of 9s.
pub fn read_cof(text: &str) -> Result<Model, CofError> {
    let mut lines = text.lines().enumerate().filter(|(_, l)| !l.trim().is_empty());
    let (_, header) = lines.next().ok_or(CofError { line: 1, reason: "empty file" })?;
    let mut head = header.split_whitespace();
    let epoch = head.next().and_then(|v| v.parse().ok()).ok_or(CofError { line: 1, reason: "bad epoch" })?;
    let name = head.next().unwrap_or_default().to_string();
    let mut rows = Vec::new();
    for (i, line) in lines {
        if line.trim_start().starts_with("9999") {
            break;
        }
        let bad = CofError { line: i + 1, reason: "expected n m g h gdot hdot" };
        let v: Vec<f64> = line.split_whitespace().map(str::parse).collect::<Result<_, _>>().map_err(|_| bad.clone())?;
        let &[n, m, g, h, gd, hd] = v.as_slice() else { return Err(bad) };
        if n < 1.0 || m < 0.0 || m > n || n.fract() != 0.0 || m.fract() != 0.0 {
            return Err(CofError { line: i + 1, reason: "bad degree or order" });
        }
        rows.push((n as usize, m as usize, [g, h, gd, hd]));
    }
    let nmax = rows.iter().map(|r| r.0).max().ok_or(CofError { line: 1, reason: "no coefficients" })?;
    let gm = REFERENCE_RADIUS * REFERENCE_RADIUS;
    let (mut main, mut secular) = (Coefficients::zeros(gm, REFERENCE_RADIUS, nmax), Coefficients::zeros(gm, REFERENCE_RADIUS, nmax));
    for (n, m, [g, h, gd, hd]) in rows {
        let k = ((2 * n + 1) as f64).sqrt().recip();
        main.set(n, m, g * k, h * k);
        secular.set(n, m, gd * k, hd * k);
    }
    Ok(Model { name, epoch, main, secular })
}

/// The World Magnetic Model 2025 (NOAA NCEI and BGS, valid 2025.0–2030.0; public domain), from its
/// official `WMM.COF` file.
pub fn wmm2025() -> Model {
    read_cof(include_str!("../data/WMM2025.COF")).expect("the bundled WMM2025.COF parses")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn noaa_test_values() {
        // All of NOAA's WMM2025_TestValues.txt: D, I (2 decimals, degrees), H, X, Y, Z, F (nT) and
        // their rates. The 6-decimal values agree to 1e-3 nT (a 1e-8 relative difference in NOAA's
        // generator that the research notes could not trace).
        let wmm = wmm2025();
        let mut rows = 0;
        for line in include_str!("../tests/data/WMM2025_TestValues.txt").lines().filter(|l| !l.starts_with('#') && !l.trim().is_empty()) {
            let v: Vec<f64> = line.split_whitespace().map(|x| x.parse().unwrap()).collect();
            let e = wmm.elements(v[0], v[2].to_radians(), v[3].to_radians(), v[1] * 1e3);
            assert!((e.declination.to_degrees() - v[4]).abs() < 0.0051 && (e.inclination.to_degrees() - v[5]).abs() < 0.0051, "{line}");
            for (got, want) in [e.h, e.x, e.y, e.z, e.f].iter().zip(&v[6..11]) {
                assert!((got - want).abs() < 1e-3, "{line}: {got} vs {want}");
            }
            for (got, want) in [e.rates[3], e.rates[0], e.rates[1], e.rates[2], e.rates[4]].iter().zip(&v[13..18]) {
                assert!((got - want).abs() < 1e-3, "{line}: rate {got} vs {want}");
            }
            assert!((e.declination_rate.to_degrees() - v[11]).abs() < 1e-4 && (e.inclination_rate.to_degrees() - v[12]).abs() < 1e-4, "{line}");
            rows += 1;
        }
        assert_eq!(rows, 100);
    }

    #[test]
    fn report_table_6() {
        // WMM2025 report Table 6 (Chulliat et al. 2025), 2025.0 at sea level, 80°N 0°E.
        let e = wmm2025().elements(2025.0, 80f64.to_radians(), 0.0, 0.0);
        for (got, want) in [(e.x, 6521.6), (e.y, 145.9), (e.z, 54791.5), (e.h, 6523.2), (e.f, 55178.5)] {
            assert!((got - want).abs() < 0.051, "{got} vs {want}");
        }
        assert!((e.inclination.to_degrees() - 83.21).abs() < 0.0051 && (e.declination.to_degrees() - 1.28).abs() < 0.0051);
    }

    #[test]
    fn dipole_and_format() {
        // From g₁₀ = −29351.8, g₁₁ = −1410.8, h₁₁ = 4545.4 nT: 7.69×10²² A m² in 2025.
        let m = wmm2025().dipole_moment(2025.0);
        assert!((m / 7.6898e22 - 1.0).abs() < 1e-3, "{m:e}");
        // The north geomagnetic pole in 2025: 80.79°N, 72.76°W (geocentric).
        let (lat, lon) = wmm2025().geomagnetic_pole(2025.0);
        assert!((lat.to_degrees() - 80.790).abs() < 1e-3 && (lon.to_degrees() + 72.757).abs() < 1e-3);
        assert_eq!(wmm2025().name, "WMM-2025");
        assert_eq!(read_cof("2025.0 X\n1 2 0 0 0 0\n").unwrap_err().line, 2);
        assert!(read_cof("").is_err() && read_cof("2025.0 X\n9999\n").is_err());
    }
}
