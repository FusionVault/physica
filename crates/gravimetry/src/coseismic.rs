//! Gravity change from an earthquake: rectangular dislocations in an elastic half-space
//! (Okubo 1992), the USGS `coulomb.inp` finite-fault format, and the equivalent surface masses a
//! satellite sees.
//!
//! ```
//! use gravimetry::coseismic::{Fault, Medium};
//!
//! // A 40 × 20 km thrust, 10 m of slip, centroid 15 km deep.
//! let fault = Fault { lat: 0.0, lon: 0.0, depth: 15e3, strike: 0.0, dip: 20.0, rake: 90.0,
//!     length: 40e3, width: 20e3, slip: 10.0, opening: 0.0 };
//! let above = fault.surface_change(0.0, 0.0, &Medium::default());
//! assert!(above.uplift > 0.5);          // the hanging wall rises (metres)
//! assert!(above.gravity > 0.0);         // and the mass above the centroid grows (m/s²)
//! ```
//!
//! The formulas follow [Okubo (1992)](https://doi.org/10.1029/92JB00178) for gravity and
//! [Okada (1985)](https://doi.org/10.1785/BSSA0750041135) for uplift, ported from François
//! Beauducel's [`okubo92.m`](https://www.mathworks.com/matlabcentral/fileexchange/37158).

// Port of okubo92.m version 1.1.0, whose licence asks that this notice be kept:
//
// Copyright (c) 2013, François Beauducel, covered by BSD License. All rights reserved.
// Redistribution and use in source and binary forms, with or without modification, are
// permitted provided that the following conditions are met: Redistributions of source code
// must retain the above copyright notice, this list of conditions and the following
// disclaimer. Redistributions in binary form must reproduce the above copyright notice, this
// list of conditions and the following disclaimer in the documentation and/or other
// materials provided with the distribution. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT
// HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
// NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

use core::f64::consts::PI;

use metrology::G;

/// Normal free-air gravity gradient, s⁻².
pub const FREE_AIR: f64 = 3.086e-6;
const EARTH_RADIUS: f64 = 6_371_000.0;

/// A rectangular fault patch and its slip. Angles in degrees, lengths in metres.
#[derive(Clone, Debug, PartialEq)]
pub struct Fault {
    /// Latitude of the patch centroid.
    pub lat: f64,
    /// Longitude of the patch centroid.
    pub lon: f64,
    /// Depth of the centroid (positive down).
    pub depth: f64,
    /// Azimuth of the fault trace from north; the fault dips to the right of it.
    pub strike: f64,
    /// Dip from horizontal, 0–90°.
    pub dip: f64,
    /// Direction the hanging wall moves, from the strike direction: 90° is a thrust, −90° a
    /// normal fault, 0° left-lateral.
    pub rake: f64,
    /// Length along strike.
    pub length: f64,
    /// Width down dip.
    pub width: f64,
    /// Slip in the rake direction.
    pub slip: f64,
    /// Tensile opening.
    pub opening: f64,
}

/// The elastic half-space.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Medium {
    /// Rock density, kg/m³ (default 2,700).
    pub density: f64,
    /// Density of whatever fills an opening crack, kg/m³ (default: the rock's).
    pub fill_density: f64,
    /// Poisson's ratio (default 0.25, the Poisson solid of Okada's check list).
    pub poisson: f64,
    /// Seawater density for the ocean-layer correction, kg/m³ (default 1,025).
    pub water_density: f64,
}

impl Default for Medium {
    fn default() -> Self {
        Medium { density: 2700.0, fill_density: 2700.0, poisson: 0.25, water_density: 1025.0 }
    }
}

/// What a dislocation does at one surface point.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfaceChange {
    /// Vertical displacement of the surface, m (Okada 1985).
    pub uplift: f64,
    /// Gravity change from the moved mass, m/s², as felt at a point that does not move with the
    /// ground: what a satellite senses. Okubo's (1992) terms for dilatation, the raised or
    /// lowered surface layer, and crack filling.
    pub gravity: f64,
}

impl SurfaceChange {
    /// What a gravimeter standing on the moving ground records: [`gravity`](Self::gravity)
    /// minus the free-air effect of its own height change, `beta · uplift` (Okubo's Δg; use
    /// [`FREE_AIR`] for `beta`).
    pub fn ground_gravity(&self, beta: f64) -> f64 {
        self.gravity - beta * self.uplift
    }

    /// Equivalent surface density, kg/m², for degrees much lower than `R / depth`: the
    /// gravity change divided by 2πG. Under the sea, raising the seafloor also pushes away a
    /// water column as tall as the uplift, so `ocean = true` subtracts `water_density · uplift`.
    pub fn surface_density(&self, medium: &Medium, ocean: bool) -> f64 {
        self.gravity / (2.0 * PI * G) - if ocean { medium.water_density * self.uplift } else { 0.0 }
    }
}

impl Fault {
    /// The change at a surface point `east`, `north` metres from the centroid (flat Earth).
    pub fn surface_change(&self, east: f64, north: f64, medium: &Medium) -> SurfaceChange {
        let (dip, strike) = (self.dip.to_radians(), self.strike.to_radians());
        let (sd, cd) = dip.sin_cos();
        let (ss, cs) = strike.sin_cos();
        let (l, w) = (self.length, self.width);
        let u1 = self.rake.to_radians().cos() * self.slip;
        let u2 = self.rake.to_radians().sin() * self.slip;
        let u3 = self.opening;
        // Centroid-relative (east, north, depth) to Okada's frame: origin at the bottom edge,
        // x along strike, d the depth of the bottom edge.
        let d = self.depth + sd * w / 2.0;
        let ec = east + cs * cd * w / 2.0;
        let nc = north - ss * cd * w / 2.0;
        let x = cs * nc + ss * ec + l / 2.0;
        let y = ss * nc - cs * ec + cd * w;
        let p = y * cd + d * sd;
        let q = y * sd - d * cd;
        let k = Terms { q, sd, cd, nu: medium.poisson };
        // Chinnery's notation: f(x, p) − f(x, p − W) − f(x − L, p) + f(x − L, p − W).
        let chin = |f: fn(&Terms, f64, f64) -> f64| f(&k, x, p) - f(&k, x, p - w) - f(&k, x - l, p) + f(&k, x - l, p - w);
        let uplift = (u1 * chin(Terms::sh) + u2 * chin(Terms::dh) + u3 * chin(Terms::th)) / (2.0 * PI);
        let gravity = medium.density * G * (u1 * chin(Terms::sg) + u2 * chin(Terms::dg) + u3 * chin(Terms::tg))
            + (medium.fill_density - medium.density) * G * u3 * chin(Terms::cg);
        SurfaceChange { uplift, gravity }
    }

    /// Seismic moment `μ · length · width · slip`, N·m, for rigidity `mu` (Pa).
    pub fn moment(&self, mu: f64) -> f64 {
        mu * self.length * self.width * self.slip
    }
}

/// The parts of Okubo's (1992) integrands that do not change across the four Chinnery corners.
struct Terms {
    q: f64,
    sd: f64,
    cd: f64,
    nu: f64,
}

impl Terms {
    fn r(&self, xi: f64, eta: f64) -> f64 {
        (xi * xi + eta * eta + self.q * self.q).sqrt()
    }

    // Displacement, Okubo (1992) eqs. 58–64, which restate Okada (1985).
    fn sh(&self, xi: f64, eta: f64) -> f64 {
        let (r, q) = (self.r(xi, eta), self.q);
        let db = eta * self.sd - q * self.cd;
        -db * q / (r * (r + eta)) - q * self.sd / (r + eta) - self.i4(db, eta, r) * self.sd
    }

    fn dh(&self, xi: f64, eta: f64) -> f64 {
        let (r, q) = (self.r(xi, eta), self.q);
        let db = eta * self.sd - q * self.cd;
        -db * q / (r * (r + xi)) - self.sd * (xi * eta / (q * r)).atan() + self.i5(xi, eta, r, db) * self.sd * self.cd
    }

    fn th(&self, xi: f64, eta: f64) -> f64 {
        let (r, q) = (self.r(xi, eta), self.q);
        let db = eta * self.sd - q * self.cd;
        (eta * self.cd + q * self.sd) * q / (r * (r + xi)) + self.cd * (xi * q / (r * (r + eta)) - (xi * eta / (q * r)).atan())
            - self.i5(xi, eta, r, db) * self.sd * self.sd
    }

    fn i4(&self, db: f64, eta: f64, r: f64) -> f64 {
        if self.cd.abs() < 1e-12 {
            -(1.0 - 2.0 * self.nu) * self.q / (r + db)
        } else {
            (1.0 - 2.0 * self.nu) / self.cd * ((r + db).ln() - self.sd * (r + eta).ln())
        }
    }

    fn i5(&self, xi: f64, eta: f64, r: f64, db: f64) -> f64 {
        if self.cd.abs() < 1e-12 {
            -(1.0 - 2.0 * self.nu) * xi * self.sd / (r + db)
        } else if xi.abs() < 1e-9 * r {
            0.0 // Okada (1985): I5 = 0 where ξ = 0
        } else {
            (1.0 - 2.0 * self.nu) * 2.0 / self.cd * ((-self.q * self.cd + (1.0 + self.sd) * (r + eta)) / (xi * self.cd)).atan()
        }
    }

    // Gravity, Okubo (1992) eqs. 52–55 and 32.
    fn i2(&self, xi: f64, eta: f64, r: f64) -> f64 {
        ((r + xi + eta) / self.q).atan()
    }

    fn sg(&self, xi: f64, eta: f64) -> f64 {
        let (r, q) = (self.r(xi, eta), self.q);
        -q * self.sd / r + q * q * self.cd / (r * (r + eta))
    }

    fn dg(&self, xi: f64, eta: f64) -> f64 {
        let (r, q) = (self.r(xi, eta), self.q);
        let db = eta * self.sd - q * self.cd;
        2.0 * self.i2(xi, eta, r) * self.sd - q * db / (r * (r + xi))
    }

    fn tg(&self, xi: f64, eta: f64) -> f64 {
        let (r, q) = (self.r(xi, eta), self.q);
        let yb = eta * self.cd + q * self.sd;
        2.0 * self.i2(xi, eta, r) * self.cd + q * yb / (r * (r + xi)) + q * xi * self.cd / (r * (r + eta))
    }

    fn cg(&self, xi: f64, eta: f64) -> f64 {
        let r = self.r(xi, eta);
        2.0 * self.i2(xi, eta, r) * self.cd - self.sd * (r + xi).ln()
    }
}

/// Why a `coulomb.inp` file could not be read.
#[derive(Debug, Clone, PartialEq)]
pub enum CoulombError {
    /// No fault elements were found.
    NoElements,
    /// No `zero lon` / `zero lat` in the map-info block, so the elements cannot be placed.
    NoOrigin,
    /// A fault element line has fewer than 11 numbers or a bad one (line number, 1-based).
    BadElement(usize),
}

impl core::fmt::Display for CoulombError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            CoulombError::NoElements => write!(f, "no fault elements in the coulomb.inp file"),
            CoulombError::NoOrigin => write!(f, "no 'zero lon' and 'zero lat' in the coulomb.inp map info"),
            CoulombError::BadElement(line) => write!(f, "bad fault element on line {line}"),
        }
    }
}

impl std::error::Error for CoulombError {}

/// Read the fault elements of a Coulomb 3 input file in the "rake, net slip" form that USGS
/// publishes with every finite-fault model (`coulomb.inp`). Each element row gives the top edge
/// from start to finish in km east and north of the file's `zero lon`/`zero lat`, then kode, rake,
/// net slip (m), dip, and top and bottom depths (km). Elements are placed on the sphere with the
/// equirectangular projection about that origin.
pub fn read_coulomb(text: &str) -> Result<Vec<Fault>, CoulombError> {
    let value = |key: &str| text.lines().find(|l| l.contains(key)).and_then(|l| l.rsplit('=').next()?.trim().parse::<f64>().ok());
    let (lon0, lat0) = (value("zero lon").ok_or(CoulombError::NoOrigin)?, value("zero lat").ok_or(CoulombError::NoOrigin)?);
    let mut faults = Vec::new();
    let mut in_table = false;
    for (i, line) in text.lines().enumerate() {
        if line.trim_start().starts_with("xxx") {
            in_table = true;
            continue;
        }
        if !in_table {
            continue;
        }
        let f: Vec<&str> = line.split_whitespace().collect();
        if f.is_empty() || f[0].parse::<i64>().is_err() {
            if !faults.is_empty() {
                break;
            }
            continue;
        }
        let num = |k: usize| f.get(k).and_then(|s| s.parse::<f64>().ok()).ok_or(CoulombError::BadElement(i + 1));
        let (xs, ys, xf, yf) = (num(1)? * 1e3, num(2)? * 1e3, num(3)? * 1e3, num(4)? * 1e3);
        let (rake, slip, dip, top, bottom) = (num(6)?, num(7)?, num(8)?, num(9)? * 1e3, num(10)? * 1e3);
        let (dx, dy) = (xf - xs, yf - ys);
        let strike = dx.atan2(dy).to_degrees().rem_euclid(360.0);
        let (sd, cd) = dip.to_radians().sin_cos();
        let width = (bottom - top) / sd;
        // Centroid: middle of the top edge, then half the width down dip (to the right of strike).
        let (sx, sy) = (strike.to_radians().cos(), -strike.to_radians().sin()); // unit vector strike + 90°
        let (ex, ny) = ((xs + xf) / 2.0 + sx * cd * width / 2.0, (ys + yf) / 2.0 + sy * cd * width / 2.0);
        faults.push(Fault {
            lat: lat0 + (ny / EARTH_RADIUS).to_degrees(),
            lon: lon0 + (ex / (EARTH_RADIUS * lat0.to_radians().cos())).to_degrees(),
            depth: (top + bottom) / 2.0,
            strike, dip, rake, length: dx.hypot(dy), width, slip, opening: 0.0,
        });
    }
    if faults.is_empty() { Err(CoulombError::NoElements) } else { Ok(faults) }
}

/// The equivalent surface masses of a set of faults on a regular grid of cell centres `lats` ×
/// `lons` (radians), as `(lat, lon, kg)` per cell, ready for [`crate::mass_coefficients`]. Each
/// fault is evaluated in its own flat frame about its centroid. `is_ocean(lat, lon)` (radians)
/// switches on the ocean-layer correction for a cell.
///
/// This is a flat-Earth, non-gravitating half-space; compared with a layered spherical Earth it
/// is usually within 10–20% at GRACE resolution for great earthquakes.
pub fn surface_masses(faults: &[Fault], lats: &[f64], lons: &[f64], medium: &Medium, is_ocean: impl Fn(f64, f64) -> bool) -> Vec<(f64, f64, f64)> {
    let step = |v: &[f64]| if v.len() > 1 { (v[1] - v[0]).abs() } else { 0.0 };
    let (dlat, dlon) = (step(lats), step(lons));
    let mut out = Vec::with_capacity(lats.len() * lons.len());
    for &lat in lats {
        let area = EARTH_RADIUS * EARTH_RADIUS * dlon * ((lat + dlat / 2.0).sin() - (lat - dlat / 2.0).sin());
        for &lon in lons {
            let mut total = SurfaceChange { uplift: 0.0, gravity: 0.0 };
            for f in faults {
                let (flat, flon) = (f.lat.to_radians(), f.lon.to_radians());
                let dl = (lon - flon + PI).rem_euclid(2.0 * PI) - PI;
                let c = f.surface_change(EARTH_RADIUS * flat.cos() * dl, EARTH_RADIUS * (lat - flat), medium);
                total.uplift += c.uplift;
                total.gravity += c.gravity;
            }
            out.push((lat, lon, total.surface_density(medium, is_ocean(lat, lon)) * area));
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A fault in Okada's own frame: bottom edge at depth `d`, x from 0 to L along strike
    /// (north), y to the east of the trace origin. Our centroid frame is shifted accordingly.
    fn okada(x: f64, y: f64, d: f64, dip: f64, l: f64, w: f64, u: [f64; 3]) -> f64 {
        let (sd, cd) = dip.to_radians().sin_cos();
        let slip = u[0].hypot(u[1]);
        let rake = u[1].atan2(u[0]).to_degrees();
        let f = Fault { lat: 0.0, lon: 0.0, depth: d - sd * w / 2.0, strike: 0.0, dip, rake, length: l, width: w, slip, opening: u[2] };
        // Okada's x runs along strike (north here) from one end; his y points up dip (west, as
        // the fault dips east) from the bottom edge. The centroid is at x = L/2, y = W·cosδ/2.
        f.surface_change(cd * w / 2.0 - y, x - l / 2.0, &Medium::default()).uplift
    }

    #[test]
    fn okada_1985_check_list() {
        // Okada (1985) Table 2, case 2: x = 2, y = 3, d = 4, δ = 70°, L = 3, W = 2, U = 1.
        let uz = |u| okada(2.0, 3.0, 4.0, 70.0, 3.0, 2.0, u);
        for (u, expect) in [([1.0, 0.0, 0.0], -2.747e-3), ([0.0, 1.0, 0.0], -3.564e-2), ([0.0, 0.0, 1.0], 3.214e-3)] {
            // The table prints four significant digits.
            assert!((uz(u) / expect - 1.0).abs() < 2e-4, "{u:?}: {:e} vs {expect:e}", uz(u));
        }
    }

    #[test]
    fn shear_faults_conserve_mass() {
        // The equivalent mass of a buried shear dislocation integrates to zero.
        let fault = Fault { lat: 0.0, lon: 0.0, depth: 20e3, strike: 30.0, dip: 25.0, rake: 95.0, length: 60e3, width: 30e3, slip: 5.0, opening: 0.0 };
        let medium = Medium::default();
        let (mut net, mut gross) = (0.0, 0.0);
        let h = 4e3;
        for i in -250..250 {
            for j in -250..250 {
                let s = fault.surface_change(i as f64 * h, j as f64 * h, &medium).surface_density(&medium, false) * h * h;
                net += s;
                gross += s.abs();
            }
        }
        assert!(net.abs() < 0.02 * gross, "net {net:e} of gross {gross:e}");
    }

    #[test]
    fn strike_slip_is_antisymmetric() {
        let f = Fault { lat: 0.0, lon: 0.0, depth: 6e3, strike: 90.0, dip: 90.0, rake: 0.0, length: 10e3, width: 10e3, slip: 5.0, opening: 0.0 };
        let m = Medium::default();
        let a = f.surface_change(3e3, 4e3, &m);
        let b = f.surface_change(-3e3, 4e3, &m);
        assert!((a.gravity + b.gravity).abs() < 1e-12 && (a.uplift + b.uplift).abs() < 1e-9, "{a:?} {b:?}");
        // Okubo (1992) Fig. 4a: tens of µGal near a 10 km fault with 5 m of slip.
        let peak = (1..20).map(|k| f.surface_change(k as f64 * 500.0, 3e3, &m).ground_gravity(FREE_AIR).abs()).fold(0.0, f64::max);
        assert!((20e-8..200e-8).contains(&peak), "{peak:e}");
    }

    #[test]
    fn reads_a_usgs_coulomb_file() {
        let text = "header\n  #   X-start    Y-start     X-fin      Y-fin   Kode  rake     netslip   dip angle     top      bot\n\
            xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx\n\
              1   137.2731    35.6409   125.2368    19.6682 100   107.1532     2.7529     6.5000     7.3777     8.5097 FFM0\n\
            \n     Grid Parameters\n3 ---  zero lon =     160.1840000\n6 ---  zero lat =      52.5050000\n";
        let f = &read_coulomb(text).unwrap()[0];
        assert!((f.strike - 217.0).abs() < 0.01 && (f.length - 20e3).abs() < 1.0 && (f.width - 10e3).abs() < 5.0, "{f:?}");
        assert!((f.depth - 7943.7).abs() < 0.1 && (f.slip - 2.7529).abs() < 1e-9);
        // The centroid is half the width down dip (azimuth 307°) from the top edge's midpoint:
        // 127.288 km east and 30.644 km north of the origin.
        assert!((f.lon - 162.0646).abs() < 1e-3 && (f.lat - 52.7806).abs() < 1e-3, "{} {}", f.lat, f.lon);
        assert_eq!(read_coulomb("nothing"), Err(CoulombError::NoOrigin));
    }
}
