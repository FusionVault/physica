//! T1: the Hunga Tonga Lamb wave in GRACE-FO laser ranging (`research/tests/T1/prereg.md`).
//!
//! The residual range acceleration is the LRI's range acceleration minus a model of it at the
//! GNV1B orbits, in the Earth-fixed frame those orbits are given in:
//!
//! `ρ̈ = ê·Δg + ê·Δa_fictitious + ê·Δa_nongravitational + (|Δv|² − ρ̇²)/ρ`,
//!
//! where Δ is satellite D minus satellite C and ê points from C to D. The gravity is the COST-G
//! field of January 2022 with GOCO06s above it, AOD1B (GLO, linear between its 3-hourly fields,
//! degrees ≥ 2), solid-Earth tides (degrees 2 and 3) and the Sun and Moon. The non-gravitational
//! part is the accelerometer's x axis, which points at the other satellite. Then a per-revolution
//! fit (bias, trend, 1/rev, 2/rev) is removed and the result band-passed to 1–10 mHz.
//!
//! **Blinding.** Every record at or after the eruption is dropped on reading, until the pipeline is
//! frozen as `prereg/T1-frozen`.

use std::f64::consts::PI;
use std::fmt::Write;
use std::path::Path;

use geopotential::{Coefficients, Legendre};
use geopotential::{julian_epoch, read_gfc};
use gravimetry::missions::{Accelerometer, Ranging, pair_noise};

use crate::svg::{esc, open, pow10};
use crate::{Result, tonga};

/// The eruption, 2022-01-15 04:14:45 UTC, in GPS seconds past 2000-01-01 12:00 GPS (UTC + 18 s).
pub const ERUPTION: f64 = 695_492_103.0;

/// Lifted only by `t1-test --unblind`, after the pipeline is frozen.
pub static UNBLIND: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);

/// Records at or after this time are dropped on reading.
fn limit() -> f64 {
    if UNBLIND.load(std::sync::atomic::Ordering::Relaxed) { f64::INFINITY } else { ERUPTION }
}
const OMEGA: f64 = 7.292_115e-5;
const GM_E: f64 = 3.986_004_418e14;
const R_E: f64 = 6_378_136.3;
const GM_SUN: f64 = 1.327_124_400_18e20;
const GM_MOON: f64 = 4.902_800_1e12;
const DT: f64 = 2.0; // LRI1B sampling, s
const NMAX: usize = 180;
const BAND: (f64, f64) = (1e-3, 10e-3);

pub(crate) type V3 = [f64; 3];

fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn add(a: V3, b: V3) -> V3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
fn scale(a: V3, k: f64) -> V3 {
    [a[0] * k, a[1] * k, a[2] * k]
}
fn norm(a: V3) -> f64 {
    dot(a, a).sqrt()
}

// ---------------------------------------------------------------------------------------------
// Level-1B and AOD1B files

/// The numeric columns `cols` of every record of a GRACE-FO Level-1B ASCII file, from records
/// strictly before [`limit`].
fn l1b(path: &Path, cols: &[usize]) -> Result<Vec<Vec<f64>>> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let body = text.split_once("# End of YAML header").ok_or("no end of YAML header")?.1;
    let mut rows = Vec::new();
    for line in body.lines().filter(|l| !l.trim().is_empty()) {
        let f: Vec<&str> = line.split_whitespace().collect();
        let t: f64 = f[0].parse()?;
        if t >= limit() {
            continue; // blinding
        }
        rows.push(cols.iter().map(|&c| f[c].parse::<f64>()).collect::<std::result::Result<_, _>>()?);
    }
    Ok(rows)
}

/// One satellite's orbit and accelerometer, indexed by whole GPS second.
struct Satellite {
    states: std::collections::HashMap<i64, (V3, V3)>,
    accel_x: std::collections::HashMap<i64, f64>,
}

fn satellite(day: &Path, date: &str, id: char) -> Result<Satellite> {
    let gnv = l1b(&day.join(format!("GNV1B_{date}_{id}_04.txt")), &[0, 3, 4, 5, 9, 10, 11])?;
    let act = l1b(&day.join(format!("ACT1B_{date}_{id}_04.txt")), &[0, 2])?;
    Ok(Satellite {
        states: gnv.iter().map(|r| (r[0] as i64, ([r[1], r[2], r[3]], [r[4], r[5], r[6]]))).collect(),
        accel_x: act.iter().map(|r| (r[0] as i64, r[1])).collect(),
    })
}

/// LRI range acceleration with its light-time and antenna corrections, and the quality flag.
fn lri(day: &Path, date: &str) -> Result<Vec<(f64, f64, bool)>> {
    let path = day.join(format!("LRI1B_{date}_Y_04.txt"));
    let rows = l1b(&path, &[0, 3, 7, 10])?;
    let text = std::fs::read_to_string(&path)?;
    let body = text.split_once("# End of YAML header").ok_or("no end of YAML header")?.1;
    let flags: Vec<bool> = body
        .lines()
        .filter(|l| !l.trim().is_empty())
        .filter(|l| l.split_whitespace().next().and_then(|t| t.parse::<f64>().ok()).is_some_and(|t| t < limit()))
        .map(|l| l.split_whitespace().last().is_some_and(|q| q.ends_with('0'))) // bit 0: phase break
        .collect();
    Ok(rows.iter().zip(flags).map(|(r, ok)| (r[0], r[1] + r[2] + r[3], ok)).collect())
}

/// AOD1B GLO fields (degrees ≥ 2) of one day, as (GPS seconds, coefficients), before the eruption.
fn aod1b(path: &Path) -> Result<Vec<(f64, Coefficients)>> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let mut out = Vec::new();
    let mut lines = text.lines();
    while let Some(line) = lines.next() {
        let Some(rest) = line.strip_prefix("DATA SET") else { continue };
        let words: Vec<&str> = rest.split_whitespace().collect();
        let count: usize = words[1].parse()?;
        let kind = *words.last().ok_or("no type")?;
        let (date, time) = (words[4], words[5]);
        let block: Vec<&str> = lines.by_ref().take(count).collect();
        if kind != "glo" {
            continue;
        }
        let t = gps_seconds(date, time)?;
        if t >= limit() + 3.0 * 3600.0 {
            continue; // the field after the eruption is needed only to interpolate up to it
        }
        let mut c = Coefficients::zeros(GM_E, 6_378_136.6, NMAX);
        for l in block {
            let f: Vec<&str> = l.split_whitespace().collect();
            let (n, m): (usize, usize) = (f[0].parse()?, f[1].parse()?);
            if (2..=NMAX).contains(&n) {
                c.set(n, m, f[2].parse()?, f[3].parse()?);
            }
        }
        out.push((t, c));
    }
    Ok(out)
}

/// GPS seconds past 2000-01-01 12:00 GPS of a GPS date `YYYY-MM-DD` and time `hh:mm:ss`.
fn gps_seconds(date: &str, time: &str) -> Result<f64> {
    let d: Vec<i64> = date.split('-').map(str::parse).collect::<std::result::Result<_, _>>()?;
    let h: Vec<f64> = time.split(':').map(str::parse).collect::<std::result::Result<_, _>>()?;
    let days = days_from_civil(d[0], d[1], d[2]) - days_from_civil(2000, 1, 1);
    Ok(days as f64 * 86_400.0 - 43_200.0 + h[0] * 3600.0 + h[1] * 60.0 + h[2])
}

fn days_from_civil(y: i64, m: i64, d: i64) -> i64 {
    let y = if m <= 2 { y - 1 } else { y };
    let era = y.div_euclid(400);
    let yoe = y - era * 400;
    let doy = (153 * (m + if m > 2 { -3 } else { 9 }) + 2) / 5 + d - 1;
    era * 146_097 + yoe * 365 + yoe / 4 - yoe / 100 + doy
}

// ---------------------------------------------------------------------------------------------
// Forces

/// Gravitational acceleration (m/s², Earth-fixed) of `field` at the Earth-fixed position `x`.
fn gravity(field: &Coefficients, x: V3) -> V3 {
    let r = norm(x);
    let (lat, lon) = ((x[2] / r).asin(), x[1].atan2(x[0]));
    let [n, e, u] = field.partials(r, lat, lon).acceleration_neu(r, lat);
    let (sl, cl, so, co) = (lat.sin(), lat.cos(), lon.sin(), lon.cos());
    [-sl * co * n - so * e + cl * co * u, -sl * so * n + co * e + cl * so * u, cl * n + sl * u]
}

/// Sun and Moon positions (m, Earth-fixed) at GPS seconds `t`: the low-precision series of
/// Montenbruck & Gill (2000, §3.3.2), rotated by Greenwich mean sidereal time.
// shortcut: J2000 positions and GMST, without precession, nutation or polar motion (errors of a
// few 10⁻⁹ m/s² in the tidal pull, at 1–2 cycles per revolution, which the per-revolution fit
// and the band-pass remove).
fn sun_moon(t: f64) -> (V3, V3) {
    let tc = (t + 51.184) / (36_525.0 * 86_400.0); // TT centuries past J2000
    let (deg, sec) = (PI / 180.0, PI / (180.0 * 3600.0));
    let eps = 23.439_291_11 * deg;
    let ecliptic = |lam: f64, beta: f64, r: f64| -> V3 {
        let (x, y, z) = (r * lam.cos() * beta.cos(), r * lam.sin() * beta.cos(), r * beta.sin());
        [x, y * eps.cos() - z * eps.sin(), y * eps.sin() + z * eps.cos()]
    };
    let m = (357.5256 + 35_999.049 * tc) * deg;
    let sun = ecliptic(282.94 * deg + m + (6892.0 * m.sin() + 72.0 * (2.0 * m).sin()) * sec, 0.0,
        (149.619 - 2.499 * m.cos() - 0.021 * (2.0 * m).cos()) * 1e9);
    let l0 = (218.316_17 + 481_267.880_88 * tc - 1.3972 * tc) * deg;
    let l = (134.962_92 + 477_198.867_53 * tc) * deg;
    let lp = (357.525_43 + 35_999.049_44 * tc) * deg;
    let f = (93.272_83 + 483_202.018_73 * tc) * deg;
    let d = (297.850_27 + 445_267.111_35 * tc) * deg;
    let lam = l0
        + (22_640.0 * l.sin() + 769.0 * (2.0 * l).sin() - 4586.0 * (l - 2.0 * d).sin() + 2370.0 * (2.0 * d).sin()
            - 668.0 * lp.sin() - 412.0 * (2.0 * f).sin() - 212.0 * (2.0 * l - 2.0 * d).sin()
            - 206.0 * (l + lp - 2.0 * d).sin() + 192.0 * (l + 2.0 * d).sin() - 165.0 * (lp - 2.0 * d).sin()
            + 148.0 * (l - lp).sin() - 125.0 * d.sin() - 110.0 * (l + lp).sin() - 55.0 * (2.0 * f - 2.0 * d).sin())
            * sec;
    let beta = (18_520.0 * (f + lam - l0 + (412.0 * (2.0 * f).sin() + 541.0 * lp.sin()) * sec).sin()
        - 526.0 * (f - 2.0 * d).sin() + 44.0 * (l + f - 2.0 * d).sin() - 31.0 * (-l + f - 2.0 * d).sin()
        - 25.0 * (-2.0 * l + f).sin() - 23.0 * (lp + f - 2.0 * d).sin() + 21.0 * (-l + f).sin()
        + 11.0 * (-lp + f - 2.0 * d).sin())
        * sec;
    let r = (385_000.0 - 20_905.0 * l.cos() - 3699.0 * (2.0 * d - l).cos() - 2956.0 * (2.0 * d).cos()
        - 570.0 * (2.0 * l).cos() + 246.0 * (2.0 * l - 2.0 * d).cos() - 205.0 * (lp - 2.0 * d).cos()
        - 171.0 * (l + 2.0 * d).cos() - 152.0 * (l + lp - 2.0 * d).cos())
        * 1e3;
    let moon = ecliptic(lam, beta, r);
    let days = (t - 5.0) / 86_400.0; // UTC ≈ UT1 days past 2000-01-01 12:00 (5 leap seconds since)
    let gmst = (280.460_618_37 + 360.985_647_366_29 * days) * deg;
    let rot = |v: V3| [gmst.cos() * v[0] + gmst.sin() * v[1], -gmst.sin() * v[0] + gmst.cos() * v[1], v[2]];
    (rot(sun), rot(moon))
}

fn third_body(x: V3, body: V3, gm: f64) -> V3 {
    let d = sub(body, x);
    sub(scale(d, gm / norm(d).powi(3)), scale(body, gm / norm(body).powi(3)))
}

/// Solid-Earth tide coefficients of degrees 2 and 3 raised by the Sun and Moon (IERS Conventions
/// 2010, eq. 6.6, nominal elastic Love numbers).
fn solid_tide(sun: V3, moon: V3) -> Coefficients {
    const K: [[f64; 4]; 4] = [[0.0; 4], [0.0; 4], [0.295_25, 0.294_70, 0.298_01, 0.0], [0.093; 4]];
    let mut c = Coefficients::zeros(GM_E, R_E, 3);
    let mut leg = Legendre::new(3);
    for (body, gm) in [(sun, GM_SUN), (moon, GM_MOON)] {
        let r = norm(body);
        let (lat, lon) = ((body[2] / r).asin(), body[1].atan2(body[0]));
        leg.compute(lat);
        for (n, kn) in K.iter().enumerate().skip(2) {
            let k = gm / GM_E * (R_E / r).powi(n as i32 + 1) / (2 * n + 1) as f64;
            for (m, knm) in kn.iter().enumerate().take(n + 1) {
                let p = knm * k * leg.p(n, m);
                c.set(n, m, c.c(n, m) + p * (m as f64 * lon).cos(), c.s(n, m) + p * (m as f64 * lon).sin());
            }
        }
    }
    c
}

fn fictitious(x: V3, v: V3) -> V3 {
    // −2ω×v − ω×(ω×x) with ω = (0, 0, Ω).
    [2.0 * OMEGA * v[1] + OMEGA * OMEGA * x[0], -2.0 * OMEGA * v[0] + OMEGA * OMEGA * x[1], 0.0]
}

// ---------------------------------------------------------------------------------------------
// Ocean tides: GOT4.10c (Ray 2013), eight main constituents

const TIDE_NMAX: usize = 120;
const RHO_SEAWATER: f64 = 1025.0;

#[derive(Clone, Copy)]
enum Nodal {
    None,
    M2,
    K1,
    O1,
    K2,
}

/// Name, multipliers of (T, s, h, p) and the phase offset (degrees) of each constituent's
/// equilibrium argument (Schureman: K1 = T + h − 90°, O1 = T − 2s + h + 90°), and its nodal
/// correction. The phases were checked on pre-eruption residuals (`tide_check`).
const CONSTITUENTS: [(&str, [f64; 4], f64, Nodal); 8] = [
    ("m2", [2.0, -2.0, 2.0, 0.0], 0.0, Nodal::M2),
    ("s2", [2.0, 0.0, 0.0, 0.0], 0.0, Nodal::None),
    ("n2", [2.0, -3.0, 2.0, 1.0], 0.0, Nodal::M2),
    ("k2", [2.0, 0.0, 2.0, 0.0], 0.0, Nodal::K2),
    ("k1", [1.0, 0.0, 1.0, 0.0], -90.0, Nodal::K1),
    ("o1", [1.0, -2.0, 1.0, 0.0], 90.0, Nodal::O1),
    ("p1", [1.0, 0.0, -1.0, 0.0], 90.0, Nodal::None),
    ("q1", [1.0, -3.0, 1.0, 1.0], 90.0, Nodal::O1),
];

/// Mean solar hour angle T, Moon's mean longitude s, Sun's h, lunar perigee p and node N, degrees,
/// at GPS seconds `t`.
fn tidal_arguments(t: f64) -> [f64; 5] {
    let tc = (t + 51.184) / (36_525.0 * 86_400.0);
    let day = ((t - 5.0) / 86_400.0 + 0.5).rem_euclid(1.0); // UTC fraction of the day since midnight
    [
        180.0 + 360.0 * day,
        218.316_447_7 + 481_267.881_234_21 * tc,
        280.466_45 + 36_000.769_748_9 * tc,
        83.353_246_5 + 4_069.013_728_7 * tc,
        125.044_52 - 1_934.136_261 * tc,
    ]
}

/// Nodal factor f and phase u (degrees) for the lunar node `n` (degrees) (Schureman's series).
fn nodal(kind: Nodal, n: f64) -> (f64, f64) {
    let n = n.to_radians();
    let (c1, s1, c2, s2, c3, s3) = (n.cos(), n.sin(), (2.0 * n).cos(), (2.0 * n).sin(), (3.0 * n).cos(), (3.0 * n).sin());
    match kind {
        Nodal::None => (1.0, 0.0),
        Nodal::M2 => (1.0004 - 0.0373 * c1 + 0.0002 * c2, -2.14 * s1),
        Nodal::K1 => (1.0060 + 0.1150 * c1 - 0.0088 * c2 + 0.0006 * c3, -8.86 * s1 + 0.68 * s2 - 0.07 * s3),
        Nodal::O1 => (1.0089 + 0.1871 * c1 - 0.0147 * c2 + 0.0014 * c3, 10.80 * s1 - 1.34 * s2 + 0.19 * s3),
        Nodal::K2 => (1.0241 + 0.2863 * c1 + 0.0083 * c2 - 0.0015 * c3, -17.74 * s1 + 0.68 * s2 - 0.04 * s3),
    }
}

/// The amplitude (cm) and Greenwich phase-lag (degrees) grids of a GOT ASCII file: 361 × 720 at
/// 0.5°, south first, 999 on land, in fixed-width fields (11 per line).
fn got_grids(path: &Path) -> Result<(Vec<f64>, Vec<f64>)> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let mut grids = Vec::new();
    let mut lines = text.lines();
    // The format line, e.g. "(11F7.2)", gives the field width.
    while let Some(width) = lines.by_ref().find_map(|l| l.trim().strip_prefix("(11F").and_then(|r| r.split('.').next()).and_then(|w| w.parse::<usize>().ok())) {
        let mut values = Vec::with_capacity(361 * 720);
        for _ in 0..361 {
            let mut row = Vec::with_capacity(720);
            while row.len() < 720 {
                let line = lines.next().ok_or("short GOT grid")?;
                for k in 0..(line.len() / width).min(11) {
                    row.push(line[width * k..width * k + width].trim().parse::<f64>()?);
                }
            }
            values.extend(row);
        }
        grids.push(values);
    }
    let mut it = grids.into_iter();
    Ok((it.next().ok_or("no amplitude grid")?, it.next().ok_or("no phase grid")?))
}

/// 4π-normalised spherical-harmonic coefficients of a 0.5° global grid (361 × 720, south first),
/// by quadrature.
fn analyse_grid(grid: &[f64], nmax: usize) -> Coefficients {
    let mut out = Coefficients::zeros(GM_E, R_E, nmax);
    let mut leg = Legendre::new(nmax);
    let dl = 0.5f64.to_radians();
    let trig: Vec<Vec<(f64, f64)>> = (0..=nmax).map(|m| (0..720).map(|j| ((m as f64 * j as f64 * dl).cos(), (m as f64 * j as f64 * dl).sin())).collect()).collect();
    let (mut c, mut s) = (vec![0.0; (nmax + 1) * (nmax + 2) / 2], vec![0.0; (nmax + 1) * (nmax + 2) / 2]);
    for i in 0..361 {
        let lat = (-90.0 + 0.5 * i as f64).to_radians();
        let w = lat.cos() * dl * dl / (4.0 * PI);
        if w <= 0.0 {
            continue;
        }
        leg.compute(lat);
        let row = &grid[i * 720..(i + 1) * 720];
        for (m, tm) in trig.iter().enumerate() {
            let (cm, sm) = row.iter().zip(tm).fold((0.0, 0.0), |(a, b), (v, (cs, sn))| (a + v * cs, b + v * sn));
            for n in m..=nmax {
                let k = n * (n + 1) / 2 + m;
                c[k] += w * leg.p(n, m) * cm;
                s[k] += w * leg.p(n, m) * sm;
            }
        }
    }
    for n in 0..=nmax {
        for m in 0..=n {
            let k = n * (n + 1) / 2 + m;
            out.set(n, m, c[k], s[k]);
        }
    }
    out
}

/// The ocean tide as potential coefficients: for each constituent, the in-phase and quadrature
/// fields, so that `ΔC(t) = f·[P cos(V+u) + Q sin(V+u)]`.
struct OceanTide {
    parts: Vec<(usize, Coefficients, Coefficients)>, // index into CONSTITUENTS, P, Q
}

impl OceanTide {
    fn load(dir: &Path) -> Result<OceanTide> {
        let love = gravimetry::love_numbers(TIDE_NMAX).ok_or("no Love numbers")?;
        let mut parts = Vec::new();
        for (i, (name, ..)) in CONSTITUENTS.iter().enumerate() {
            let (amp, phase) = got_grids(&dir.join(format!("{name}.d")))?;
            let comp = |k: usize, f: fn(f64) -> f64| if amp[k] >= 999.0 { 0.0 } else { amp[k] * 0.01 * f(phase[k].to_radians()) };
            let a: Vec<f64> = (0..amp.len()).map(|k| comp(k, f64::cos)).collect();
            let b: Vec<f64> = (0..amp.len()).map(|k| comp(k, f64::sin)).collect();
            let (mut p, mut q) = (analyse_grid(&a, TIDE_NMAX), analyse_grid(&b, TIDE_NMAX));
            for (n, k_load) in love.iter().enumerate().take(TIDE_NMAX + 1) {
                let f = if n < 2 { 0.0 } else { 3.0 * RHO_SEAWATER * (1.0 + k_load) / (geopotential::RHO_EARTH * R_E * (2 * n + 1) as f64) };
                for m in 0..=n {
                    p.set(n, m, f * p.c(n, m), f * p.s(n, m));
                    q.set(n, m, f * q.c(n, m), f * q.s(n, m));
                }
            }
            parts.push((i, p, q));
        }
        Ok(OceanTide { parts })
    }

    /// Cosine and sine weights `f·cos(V+u)`, `f·sin(V+u)` of each part at time `t`.
    fn weights(&self, t: f64) -> Vec<(f64, f64)> {
        let a = tidal_arguments(t);
        self.parts.iter().map(|(i, ..)| {
            let (_, mult, offset, kind) = CONSTITUENTS[*i];
            let (f, u) = nodal(kind, a[4]);
            let v = (mult[0] * a[0] + mult[1] * a[1] + mult[2] * a[2] + mult[3] * a[3] + offset + u).to_radians();
            (f * v.cos(), f * v.sin())
        }).collect()
    }

    fn add_to(&self, t: f64, buf: &mut Coefficients) {
        for ((_, p, q), (wc, ws)) in self.parts.iter().zip(self.weights(t)) {
            for n in 2..=TIDE_NMAX {
                for m in 0..=n {
                    buf.set(n, m, buf.c(n, m) + wc * p.c(n, m) + ws * q.c(n, m), buf.s(n, m) + wc * p.s(n, m) + ws * q.s(n, m));
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------------------------
// Gravity model and dynamic orbits

/// The time-variable gravity model: the static field plus each AOD1B field, blended linearly in time.
struct Model {
    fields: Vec<(f64, Coefficients)>,
    tide: Option<OceanTide>,
}

impl Model {
    fn new(field: &Coefficients, aod: &[(f64, Coefficients)]) -> Model {
        let fields = aod.iter().map(|(t, a)| {
            let a = a.rescaled(field.gm, field.radius);
            let mut f = field.clone();
            for n in 2..=NMAX {
                for m in 0..=n {
                    f.set(n, m, f.c(n, m) + a.c(n, m), f.s(n, m) + a.s(n, m));
                }
            }
            (*t, f)
        });
        Model { fields: fields.collect(), tide: None }
    }

    fn covers(&self, t: f64) -> bool {
        self.fields.first().is_some_and(|f| f.0 <= t) && self.fields.last().is_some_and(|f| t < f.0)
    }

    /// The field at time `t`, written into `buf` (degree `NMAX`). `t` must be covered.
    fn at(&self, t: f64, buf: &mut Coefficients) {
        let j = self.fields.iter().rposition(|f| f.0 <= t).unwrap_or(0).min(self.fields.len() - 2);
        let (a, b) = (&self.fields[j], &self.fields[j + 1]);
        let w = (t - a.0) / (b.0 - a.0);
        for n in 0..=NMAX {
            for m in 0..=n {
                buf.set(n, m, (1.0 - w) * a.1.c(n, m) + w * b.1.c(n, m), (1.0 - w) * a.1.s(n, m) + w * b.1.s(n, m));
            }
        }
        if let Some(tide) = &self.tide {
            tide.add_to(t, buf);
        }
    }
}

/// Solid tides, Sun and Moon at time `t` and position `x`.
fn tides(t: f64, x: V3) -> V3 {
    let (sun, moon) = sun_moon(t);
    add(gravity(&solid_tide(sun, moon), x), add(third_body(x, sun, GM_SUN), third_body(x, moon, GM_MOON)))
}

/// Empirical accelerations in the radial, along-track and cross-track directions: a constant and
/// a once-per-revolution cosine and sine in each, in units of 10⁻⁸ m/s².
fn empirical(x: V3, v: V3, u: f64, e: &[f64]) -> V3 {
    let r = scale(x, 1.0 / norm(x));
    let h = [x[1] * v[2] - x[2] * v[1], x[2] * v[0] - x[0] * v[2], x[0] * v[1] - x[1] * v[0]];
    let n = scale(h, 1.0 / norm(h));
    let along = [n[1] * r[2] - n[2] * r[1], n[2] * r[0] - n[0] * r[2], n[0] * r[1] - n[1] * r[0]];
    let (cu, su) = (u.cos(), u.sin());
    let k = |i: usize| 1e-8 * (e[i] + e[i + 1] * cu + e[i + 2] * su);
    add(scale(r, k(0)), add(scale(along, k(3)), scale(n, k(6))))
}

/// Point mass and J₂ only: the cheap model the fit's partial derivatives come from.
fn j2_gravity(x: V3) -> V3 {
    const J2: f64 = 1.082_626_7e-3;
    let r2 = dot(x, x);
    let (k, z2, j) = (GM_E / (r2 * r2.sqrt()), x[2] * x[2] / r2, 1.5 * J2 * R_E * R_E / r2);
    [-k * x[0] * (1.0 + j * (1.0 - 5.0 * z2)), -k * x[1] * (1.0 + j * (1.0 - 5.0 * z2)), -k * x[2] * (1.0 + j * (3.0 - 5.0 * z2))]
}

/// Fourth-order Runge–Kutta for `ẍ = a(t, x, v)` in steps of `DT`: the states at `t0 + k·DT`,
/// k = 0…steps.
fn integrate(t0: f64, steps: usize, x0: V3, v0: V3, a: &mut dyn FnMut(f64, V3, V3) -> V3) -> Vec<(V3, V3)> {
    let mut out = Vec::with_capacity(steps + 1);
    let (mut x, mut v) = (x0, v0);
    out.push((x, v));
    for k in 0..steps {
        let t = t0 + k as f64 * DT;
        let h = DT;
        let a1 = a(t, x, v);
        let (x2, v2) = (add(x, scale(v, h / 2.0)), add(v, scale(a1, h / 2.0)));
        let a2 = a(t + h / 2.0, x2, v2);
        let (x3, v3) = (add(x, scale(v2, h / 2.0)), add(v, scale(a2, h / 2.0)));
        let a3 = a(t + h / 2.0, x3, v3);
        let (x4, v4) = (add(x, scale(v3, h)), add(v, scale(a3, h)));
        let a4 = a(t + h, x4, v4);
        x = add(x, scale(add(add(v, scale(v2, 2.0)), add(scale(v3, 2.0), v4)), h / 6.0));
        v = add(v, scale(add(add(a1, scale(a2, 2.0)), add(scale(a3, 2.0), a4)), h / 6.0));
        out.push((x, v));
    }
    out
}

/// An arc to fit: start time, number of `DT` steps, and the revolution's node time and period.
#[derive(Clone, Copy)]
struct Arc {
    t0: f64,
    steps: usize,
    node: f64,
    period: f64,
}

/// A dynamic orbit fitted to the GNV1B positions (every 10 s) over `arc`: initial state and nine
/// empirical accelerations, by Gauss–Newton with partials from the J₂ model. Returns the states
/// at every step and the rms fit to GNV1B, m.
fn fit_orbit(model: &Model, gnv: &States, arc: Arc, buf: &mut Coefficients) -> Option<(Vec<(V3, V3)>, f64)> {
    if !model.covers(arc.t0) || !model.covers(arc.t0 + arc.steps as f64 * DT) {
        return None;
    }
    let obs: Vec<(usize, V3)> = (0..=arc.steps).step_by(5).filter_map(|k| gnv.get(&((arc.t0 + k as f64 * DT) as i64)).map(|s| (k, s.0))).collect();
    if obs.len() < arc.steps / 10 {
        return None;
    }
    let &(x0, v0) = gnv.get(&(arc.t0 as i64))?;
    let mut p: Vec<f64> = x0.iter().chain(&v0).copied().chain([0.0; 9]).collect();
    let scales = [1.0, 1.0, 1.0, 1e-3, 1e-3, 1e-3, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
    let u = |t: f64| 2.0 * PI * (t - arc.node) / arc.period;
    let run = |p: &[f64], full: bool, buf: &mut Coefficients| -> Vec<(V3, V3)> {
        let e = p[6..].to_vec();
        integrate(arc.t0, arc.steps, [p[0], p[1], p[2]], [p[3], p[4], p[5]], &mut |t, x, v| {
            let g = if full {
                model.at(t, buf);
                add(gravity(buf, x), tides(t, x))
            } else {
                j2_gravity(x)
            };
            add(g, add(fictitious(x, v), empirical(x, v, u(t), &e)))
        })
    };
    let misfit = |states: &[(V3, V3)]| -> Vec<f64> { obs.iter().flat_map(|&(k, o)| sub(o, states[k].0)).collect() };
    for _ in 0..2 {
        let y = misfit(&run(&p, true, buf));
        let base = run(&p, false, buf);
        let columns: Vec<Vec<f64>> = (0..15).map(|i| {
            let mut q = p.clone();
            q[i] += scales[i];
            let moved = run(&q, false, buf);
            obs.iter().flat_map(|&(k, _)| sub(moved[k].0, base[k].0)).collect()
        }).collect();
        let rows: Vec<Vec<f64>> = (0..y.len()).map(|r| columns.iter().map(|c| c[r]).collect()).collect();
        let (dq, _) = crate::shared::ols(&rows, &y).ok()?;
        for i in 0..15 {
            p[i] += dq[i] * scales[i];
        }
    }
    let states = run(&p, true, buf);
    let y = misfit(&states);
    Some((states, (y.iter().map(|v| v * v).sum::<f64>() / y.len() as f64).sqrt()))
}

type States = std::collections::HashMap<i64, (V3, V3)>;

/// Dynamic orbits for every revolution (between consecutive `nodes`), each fitted over the
/// revolution plus 30 minutes either side and kept for the revolution itself. Runs in parallel.
fn fitted_orbits(model: &Model, gnv: &States, nodes: &[f64]) -> (States, Vec<f64>) {
    let template = model.fields[0].1.clone();
    let threads = std::thread::available_parallelism().map_or(4, |n| n.get().min(16));
    let windows: Vec<(f64, f64)> = nodes.windows(2).map(|w| (w[0], w[1])).filter(|w| w.1 - w.0 < 7000.0).collect();
    let chunk = windows.len().div_ceil(threads).max(1);
    let results: Vec<(States, Vec<f64>)> = std::thread::scope(|scope| {
        let handles: Vec<_> = windows.chunks(chunk).map(|part| {
            let mut buf = template.clone();
            scope.spawn(move || {
                let (mut states, mut fits) = (States::new(), Vec::new());
                for &(a, b) in part {
                    let t0 = 2.0 * ((a - 1800.0) / 2.0).floor();
                    let steps = ((b + 1800.0 - t0) / DT).ceil() as usize;
                    let Some((orbit, rms)) = fit_orbit(model, gnv, Arc { t0, steps, node: a, period: b - a }, &mut buf) else { continue };
                    fits.push(rms);
                    for (k, s) in orbit.into_iter().enumerate() {
                        let t = t0 + k as f64 * DT;
                        if t >= a && t < b {
                            states.insert(t as i64, s);
                        }
                    }
                }
                (states, fits)
            })
        }).collect();
        handles.into_iter().filter_map(|h| h.join().ok()).collect()
    });
    results.into_iter().fold((States::new(), Vec::new()), |(mut s, mut f), (s2, f2)| {
        s.extend(s2);
        f.extend(f2);
        (s, f)
    })
}

// ---------------------------------------------------------------------------------------------
// Residuals

/// One epoch's residual and the size of each model term along the line of sight.
#[derive(Clone, Copy, Default)]
struct Epoch {
    t: f64,
    residual: f64,
    terms: [f64; 5], // gravity (field + AOD1B), tides + Sun + Moon, rotating frame, accelerometers, velocity
    mid: V3,         // the pair's midpoint, Earth-fixed
    vel: V3,         // its velocity
}

struct Inputs<'a> {
    model: &'a Model,
    c: &'a States,
    d: &'a States,
    c_acc: &'a std::collections::HashMap<i64, f64>,
    d_acc: &'a std::collections::HashMap<i64, f64>,
}

fn epoch(inp: &Inputs, buf: &mut Coefficients, t: f64, observed: f64) -> Option<Epoch> {
    let s = t as i64;
    let (&(rc, vc), &(rd, vd)) = (inp.c.get(&s)?, inp.d.get(&s)?);
    let (axc, axd) = (*inp.c_acc.get(&s)?, *inp.d_acc.get(&s)?);
    if !inp.model.covers(t) {
        return None;
    }
    inp.model.at(t, buf);
    let dr = sub(rd, rc);
    let rho = norm(dr);
    let e = scale(dr, 1.0 / rho);
    let dv = sub(vd, vc);
    let along = |f: &dyn Fn(V3, V3) -> V3| dot(e, sub(f(rd, vd), f(rc, vc)));
    let terms = [
        along(&|x, _| gravity(buf, x)),
        along(&|x, _| tides(t, x)),
        along(&fictitious),
        -axd - axc, // each accelerometer's x axis points at the other satellite
        (dot(dv, dv) - dot(e, dv).powi(2)) / rho,
    ];
    Some(Epoch { t, residual: observed - terms.iter().sum::<f64>(), terms, mid: scale(add(rc, rd), 0.5), vel: scale(add(vc, vd), 0.5) })
}

/// Residuals at every good LRI epoch, computed in parallel.
fn residuals(inp: &Inputs, obs: &[(f64, f64, bool)]) -> Vec<Epoch> {
    let threads = std::thread::available_parallelism().map_or(4, |n| n.get().min(16));
    let chunk = obs.len().div_ceil(threads).max(1);
    let template = inp.model.fields[0].1.clone();
    let mut out: Vec<Epoch> = std::thread::scope(|scope| {
        let handles: Vec<_> = obs.chunks(chunk).map(|part| {
            let mut buf = template.clone();
            scope.spawn(move || part.iter().filter(|o| o.2).filter_map(|o| epoch(inp, &mut buf, o.0, o.1)).collect::<Vec<_>>())
        }).collect();
        handles.into_iter().flat_map(|h| h.join().unwrap_or_default()).collect()
    });
    out.sort_by(|a, b| a.t.total_cmp(&b.t));
    out
}

/// The epochs of each revolution (between consecutive ascending crossings in `nodes`), with the
/// revolution's start and period. `times` must be sorted.
fn revolutions(times: &[f64], nodes: &[f64]) -> Vec<(std::ops::Range<usize>, f64, f64)> {
    nodes.windows(2).filter_map(|w| {
        let (i, j) = (times.partition_point(|&t| t < w[0]), times.partition_point(|&t| t < w[1]));
        (j - i >= 50).then(|| (i..j, w[0], w[1] - w[0]))
    }).collect()
}

/// Removes, per revolution, a least-squares fit of a bias, a trend and 1/rev and 2/rev sines and
/// cosines from `values` (sampled at `times`).
fn detrend(times: &[f64], values: &mut [f64], revs: &[(std::ops::Range<usize>, f64, f64)]) {
    for (range, a, period) in revs {
        let basis = |t: f64| {
            let u = 2.0 * PI * (t - a) / period;
            vec![1.0, (t - a) / period, u.cos(), u.sin(), (2.0 * u).cos(), (2.0 * u).sin()]
        };
        let rows: Vec<Vec<f64>> = times[range.clone()].iter().map(|&t| basis(t)).collect();
        if values[range.clone()].iter().all(|v| *v == 0.0) {
            continue;
        }
        let Ok((beta, _)) = crate::shared::ols(&rows, &values[range.clone()]) else { continue };
        for (row, v) in rows.iter().zip(&mut values[range.clone()]) {
            *v -= row.iter().zip(&beta).map(|(x, b)| x * b).sum::<f64>();
        }
    }
}

fn remove_per_revolution(epochs: &mut [Epoch], nodes: &[f64]) {
    let times: Vec<f64> = epochs.iter().map(|e| e.t).collect();
    let mut values: Vec<f64> = epochs.iter().map(|e| e.residual).collect();
    detrend(&times, &mut values, &revolutions(&times, nodes));
    for (e, v) in epochs.iter_mut().zip(values) {
        e.residual = v;
    }
}

/// Times of the ascending equator crossings of an orbit.
fn ascending_nodes(sat: &Satellite) -> Vec<f64> {
    let mut keys: Vec<i64> = sat.states.keys().copied().collect();
    keys.sort_unstable();
    keys.windows(2)
        .filter(|w| w[1] == w[0] + 1)
        .filter_map(|w| {
            let (z0, z1) = (sat.states[&w[0]].0[2], sat.states[&w[1]].0[2]);
            (z0 < 0.0 && z1 >= 0.0).then(|| w[0] as f64 + z0 / (z0 - z1))
        })
        .collect()
}

// ---------------------------------------------------------------------------------------------
// Filtering and spectra

/// A biquad section `(b0, b1, b2, a1, a2)`, normalised by a0 (RBJ cookbook).
type Biquad = [f64; 5];

/// The four biquads of a 4th-order Butterworth band-pass (4th-order high-pass at `lo`, then
/// 4th-order low-pass at `hi`, Hz), for sampling interval `dt`.
fn butterworth_band(lo: f64, hi: f64, dt: f64) -> Vec<Biquad> {
    let qs = [1.0 / (2.0 * (PI / 8.0).cos()), 1.0 / (2.0 * (3.0 * PI / 8.0).cos())];
    let mut out = Vec::new();
    for (f, high) in [(lo, true), (hi, false)] {
        let w0 = 2.0 * PI * f * dt;
        for q in qs {
            let alpha = w0.sin() / (2.0 * q);
            let a0 = 1.0 + alpha;
            let (b0, b1) = if high { ((1.0 + w0.cos()) / 2.0, -(1.0 + w0.cos())) } else { ((1.0 - w0.cos()) / 2.0, 1.0 - w0.cos()) };
            out.push([b0 / a0, b1 / a0, b0 / a0, -2.0 * w0.cos() / a0, (1.0 - alpha) / a0]);
        }
    }
    out
}

fn run(sections: &[Biquad], x: &mut [f64]) {
    for s in sections {
        let (mut x1, mut x2, mut y1, mut y2) = (0.0, 0.0, 0.0, 0.0);
        for v in x.iter_mut() {
            let y = s[0] * *v + s[1] * x1 + s[2] * x2 - s[3] * y1 - s[4] * y2;
            (x2, x1, y2, y1) = (x1, *v, y1, y);
            *v = y;
        }
    }
}

/// Zero-phase filtering (forward, then backward), with odd reflection padding at both ends.
fn filtfilt(sections: &[Biquad], x: &[f64], pad: usize) -> Vec<f64> {
    let n = x.len();
    let pad = pad.min(n - 1);
    let mut y: Vec<f64> = (1..=pad).rev().map(|i| 2.0 * x[0] - x[i])
        .chain(x.iter().copied())
        .chain((n - 1 - pad..n - 1).rev().map(|i| 2.0 * x[n - 1] - x[i]))
        .collect();
    run(sections, &mut y);
    y.reverse();
    run(sections, &mut y);
    y.reverse();
    y[pad..pad + n].to_vec()
}

/// Contiguous runs (no step longer than 1.5 `DT`) of at least `min_len` samples.
fn segments(epochs: &[Epoch], min_len: usize) -> Vec<std::ops::Range<usize>> {
    let mut out = Vec::new();
    let mut start = 0;
    for i in 1..=epochs.len() {
        if i == epochs.len() || epochs[i].t - epochs[i - 1].t > 1.5 * DT {
            if i - start >= min_len {
                out.push(start..i);
            }
            start = i;
        }
    }
    out
}

/// One-sided amplitude spectral density by Welch's method (Hann window, segments of `len`
/// samples, half overlap) over several series, as (frequency, ASD).
fn welch(series: &[&[f64]], len: usize) -> Vec<(f64, f64)> {
    let w: Vec<f64> = (0..len).map(|i| 0.5 - 0.5 * (2.0 * PI * i as f64 / len as f64).cos()).collect();
    let norm2: f64 = w.iter().map(|v| v * v).sum();
    let mut acc = vec![0.0; len / 2];
    let mut count = 0usize;
    for s in series {
        let mut start = 0;
        while start + len <= s.len() {
            let mut re: Vec<f64> = s[start..start + len].iter().zip(&w).map(|(x, w)| x * w).collect();
            let mean = re.iter().sum::<f64>() / len as f64;
            re.iter_mut().zip(&w).for_each(|(x, w)| *x -= mean * w);
            let mut im = vec![0.0; len];
            gravimetry::noise::fft(&mut re, &mut im, false);
            for (k, a) in acc.iter_mut().enumerate() {
                *a += 2.0 * DT * (re[k] * re[k] + im[k] * im[k]) / norm2;
            }
            count += 1;
            start += len / 2;
        }
    }
    (1..len / 2).map(|k| (k as f64 / (len as f64 * DT), (acc[k] / count.max(1) as f64).sqrt())).collect()
}

fn rms(x: impl Iterator<Item = f64>) -> f64 {
    let (s, n) = x.fold((0.0, 0usize), |(s, n), v| (s + v * v, n + 1));
    (s / n.max(1) as f64).sqrt()
}

// ---------------------------------------------------------------------------------------------
// Stage 1

/// The static field: COST-G January 2022 to its maximum degree, GOCO06s above it to `NMAX`.
fn static_field(goco: &Path, costg: &Path) -> Result<(Coefficients, usize)> {
    let mut field = read_gfc(&std::fs::read_to_string(goco)?)?.model.at(julian_epoch(2022, 1, 15, 0, 0)).truncated(NMAX);
    let month = read_gfc(&std::fs::read_to_string(costg)?)?.model.base.rescaled(field.gm, field.radius);
    // T1_MONTHLY_MAX: how far up the monthly field replaces GOCO06s (blind tuning on pre-eruption data).
    let cap = std::env::var("T1_MONTHLY_MAX").ok().and_then(|v| v.parse().ok()).unwrap_or(90usize);
    let top = month.max_degree().min(NMAX).min(cap);
    for n in 0..=top {
        for m in 0..=n {
            field.set(n, m, month.c(n, m), month.s(n, m));
        }
    }
    Ok((field, top))
}

/// Every AOD1B GLO field in `dir`, in time order.
fn all_aod1b(dir: &Path) -> Result<Vec<(f64, Coefficients)>> {
    let mut files: Vec<_> = std::fs::read_dir(dir)?.filter_map(|e| e.ok().map(|e| e.path())).filter(|p| p.extension().is_some_and(|x| x == "asc")).collect();
    files.sort();
    let mut aod = Vec::new();
    for f in files {
        aod.extend(aod1b(&f)?);
    }
    aod.sort_by(|a, b| a.0.total_cmp(&b.0));
    Ok(aod)
}

/// LRI epochs: time, corrected range acceleration, and whether the record is good.
type Lri = Vec<(f64, f64, bool)>;

/// The two satellites and the LRI over several days, merged.
fn load_days(l1b_dir: &Path, dates: &[&str]) -> Result<(Satellite, Satellite, Lri)> {
    let empty = || Satellite { states: States::new(), accel_x: std::collections::HashMap::new() };
    let (mut c, mut d, mut obs) = (empty(), empty(), Vec::new());
    for date in dates {
        let day = l1b_dir.join(date);
        for (sat, id) in [(&mut c, 'C'), (&mut d, 'D')] {
            let s = satellite(&day, date, id)?;
            sat.states.extend(s.states);
            sat.accel_x.extend(s.accel_x);
        }
        obs.extend(lri(&day, date)?);
    }
    Ok((c, d, obs))
}

/// Per-revolution fit, band-pass and spectrum of a residual series: (band rms, ASD, hours used).
fn analyse(epochs: &mut [Epoch], nodes: &[f64]) -> (f64, Vec<(f64, f64)>, f64) {
    remove_per_revolution(epochs, nodes);
    let sections = butterworth_band(BAND.0, BAND.1, DT);
    let segs = segments(epochs, 4096);
    let filtered: Vec<Vec<f64>> = segs.iter().map(|r| {
        let x: Vec<f64> = epochs[r.clone()].iter().map(|e| e.residual).collect();
        filtfilt(&sections, &x, 1500)[750..x.len() - 750].to_vec()
    }).collect();
    let raw: Vec<Vec<f64>> = segs.iter().map(|r| epochs[r.clone()].iter().map(|e| e.residual).collect()).collect();
    let asd = welch(&raw.iter().map(Vec::as_slice).collect::<Vec<_>>(), 2048);
    let hours = filtered.iter().map(Vec::len).sum::<usize>() as f64 * DT / 3600.0;
    (rms(filtered.iter().flatten().copied()), asd, hours)
}

/// Checks the tide model's phases on the residual without it: regresses the band-passed residual on
/// each constituent's predicted signal and on the same signal shifted by 90°. Right phases give
/// coefficients near (1, 0).
fn tide_check(report: &mut String, model: &Model, tide: &OceanTide, c: &Satellite, d: &Satellite, obs: &[(f64, f64, bool)], nodes: &[f64]) -> Result<()> {
    let inp = Inputs { model, c: &c.states, d: &d.states, c_acc: &c.accel_x, d_acc: &d.accel_x };
    let epochs = residuals(&inp, obs);
    let times: Vec<f64> = epochs.iter().map(|e| e.t).collect();
    let threads = std::thread::available_parallelism().map_or(4, |n| n.get().min(16));
    let chunk = times.len().div_ceil(threads).max(1);
    let rows: Vec<Vec<f64>> = std::thread::scope(|scope| {
        let handles: Vec<_> = times.chunks(chunk).map(|part| scope.spawn(move || {
            part.iter().map(|&t| {
                let s = t as i64;
                let (rc, rd) = (c.states[&s].0, d.states[&s].0);
                let e = scale(sub(rd, rc), 1.0 / norm(sub(rd, rc)));
                let w = tide.weights(t);
                // Per species (diurnal, semidiurnal): predicted signal and its 90°-shifted copy.
                let mut row = vec![0.0; 4];
                for ((i, p, q), (wc, ws)) in tide.parts.iter().zip(&w) {
                    let gp = dot(e, sub(gravity(p, rd), gravity(p, rc)));
                    let gq = dot(e, sub(gravity(q, rd), gravity(q, rc)));
                    let species = if CONSTITUENTS[*i].1[0] == 1.0 { 0 } else { 2 };
                    row[species] += wc * gp + ws * gq;
                    row[species + 1] += -ws * gp + wc * gq;
                }
                row
            }).collect::<Vec<_>>()
        })).collect();
        handles.into_iter().flat_map(|h| h.join().unwrap_or_default()).collect()
    });
    let revs = revolutions(&times, nodes);
    let prep = |mut v: Vec<f64>| {
        detrend(&times, &mut v, &revs);
        band_pass(&times, &v)
    };
    let y = prep(epochs.iter().map(|e| e.residual).collect());
    let columns: Vec<Vec<f64>> = (0..4).map(|j| prep(rows.iter().map(|r| r[j]).collect())).collect();
    // In nm/s²: shared::ols has an absolute pivot threshold.
    let design: Vec<Vec<f64>> = (0..y.len()).map(|i| columns.iter().map(|col| 1e9 * col[i]).collect()).collect();
    let y_nm: Vec<f64> = y.iter().map(|v| 1e9 * v).collect();
    let (beta, resid) = crate::shared::ols(&design, &y_nm)?;
    let resid: Vec<f64> = resid.iter().map(|v| 1e-9 * v).collect();
    writeln!(report, "Tide phase check (no-tide residual on predicted / 90°-shifted signal; expect 1 / 0):")?;
    for (j, name) in ["diurnal", "semidiurnal"].iter().enumerate() {
        let size = rms(columns[2 * j].iter().copied());
        writeln!(report, "  {name:<12} in-phase {:>6.2}  quadrature {:>6.2}   (signal rms {size:.2e} m/s² in band)", beta[2 * j], beta[2 * j + 1])?;
    }
    writeln!(report, "  band rms: {:.3e} before, {:.3e} after the fit", rms(y.iter().copied()), rms(resid.iter().copied()))?;
    Ok(())
}

/// The blind noise study: residuals on the given pre-eruption days, with the GNV1B orbits as
/// given and with dynamic orbits fitted to them.
pub fn noise_study(l1b_dir: &Path, aod_dir: &Path, goco: &Path, costg: &Path, dates: &[&str]) -> Result<(String, String)> {
    let (field, top) = static_field(goco, costg)?;
    let mut model = Model::new(&field, &all_aod1b(aod_dir)?);
    let (c, d, obs) = load_days(l1b_dir, dates)?;
    let mut nodes = ascending_nodes(&c);
    nodes.sort_by(f64::total_cmp);
    let mut report = String::new();
    writeln!(report, "T1 stage 1 (blind): residual range acceleration before the eruption, {}", dates.join(", "))?;
    writeln!(report, "Model: COST-G January 2022 to degree {top}, GOCO06s to {NMAX}; AOD1B RL06 GLO; solid tides; Sun and Moon.")?;
    writeln!(report, "{} LRI epochs, {} with a phase break.", obs.len(), obs.iter().filter(|o| !o.2).count())?;
    if let Ok(dir) = std::env::var("T1_TIDE_DIR") {
        let tide = OceanTide::load(Path::new(&dir))?;
        writeln!(report, "Ocean tides: GOT4.10c, 8 constituents to degree {TIDE_NMAX}.")?;
        tide_check(&mut report, &model, &tide, &c, &d, &obs, &nodes)?;
        model.tide = Some(tide);
    }

    let quick = std::env::var("T1_QUICK").is_ok();
    let (fit_c, rms_c) = if quick { (States::new(), vec![0.0]) } else { fitted_orbits(&model, &c.states, &nodes) };
    let (fit_d, rms_d) = if quick { (States::new(), vec![0.0]) } else { fitted_orbits(&model, &d.states, &nodes) };
    let median = |v: &[f64]| crate::shared::quantile(v, 0.5);
    writeln!(report, "Dynamic orbits: {} and {} revolutions fitted; median rms fit to GNV1B {:.3} m (C), {:.3} m (D).",
        rms_c.len(), rms_d.len(), median(&rms_c), median(&rms_d))?;

    let names = ["gravity (field + AOD1B)", "tides, Sun and Moon", "rotating frame", "accelerometers", "velocity term"];
    let mut curves = Vec::new();
    for (label, sc, sd) in [("GNV1B orbits", &c.states, &d.states), ("fitted orbits", &fit_c, &fit_d)] {
        if sc.is_empty() {
            continue;
        }
        let inp = Inputs { model: &model, c: sc, d: sd, c_acc: &c.accel_x, d_acc: &d.accel_x };
        let mut epochs = residuals(&inp, &obs);
        writeln!(report, "\n{label}: {} residuals", epochs.len())?;
        let sections = butterworth_band(BAND.0, BAND.1, DT);
        for (k, name) in names.iter().enumerate() {
            let band: Vec<f64> = segments(&epochs, 4096).iter().flat_map(|r| {
                let x: Vec<f64> = epochs[r.clone()].iter().map(|e| e.terms[k]).collect();
                filtfilt(&sections, &x, 1500)[750..x.len() - 750].to_vec()
            }).collect();
            writeln!(report, "  {name:<26} rms {:>10.3e}   in 1–10 mHz {:>10.3e}", rms(epochs.iter().map(|e| e.terms[k])), rms(band.into_iter()))?;
        }
        writeln!(report, "  {:<26} rms {:>10.3e}", "residual", rms(epochs.iter().map(|e| e.residual)))?;
        let (band_rms, asd, hours) = analyse(&mut epochs, &nodes);
        writeln!(report, "  {:<26} rms {:>10.3e}   ({hours:.1} h)", "band-passed residual", band_rms)?;
        for f in [1e-3, 2e-3, 5e-3, 10e-3, 20e-3, 50e-3] {
            let a = asd.iter().min_by(|x, y| (x.0 - f).abs().total_cmp(&(y.0 - f).abs())).map_or(0.0, |x| x.1);
            let nominal = pair_noise(Accelerometer::SuperStar, Ranging::GraceFo, f);
            writeln!(report, "    {:>5.1} mHz  ASD {a:.2e}  nominal {nominal:.2e}  ratio {:.1}", f * 1e3, a / nominal)?;
        }
        curves.push((label, asd, band_rms));
    }
    Ok((report, figure(&curves)))
}

// ---------------------------------------------------------------------------------------------
// Stage 2: the pre-registered test

const LAMB_SPEED: f64 = 318.2;
const R_MEAN: f64 = 6_371e3;
const TONGA: (f64, f64) = (-20.546, -175.390);
const WIDTHS_MIN: [f64; 3] = [2.5, 5.0, 10.0];
const WINDOW: f64 = 86_400.0;
const AOD_STEP: f64 = 10_800.0;
const AOD_ORIGIN: f64 = 695_304_000.0; // 2022-01-13 00:00 GPS, an AOD1B epoch

fn unit_vector(lat_deg: f64, lon_deg: f64) -> V3 {
    let (la, lo) = (lat_deg.to_radians(), lon_deg.to_radians());
    [la.cos() * lo.cos(), la.cos() * lo.sin(), la.sin()]
}

/// The front's angular distance from its source at `t` for an origin at `t0`, folded into [0, π]
/// (past the antipode it comes back).
fn front(t: f64, t0: f64) -> f64 {
    let a = (LAMB_SPEED * (t - t0) / R_MEAN).rem_euclid(2.0 * PI);
    if a > PI { 2.0 * PI - a } else { a }
}

/// A crossing of a front by the pair: time, angular distance from the source, sin of the angle
/// between the ground track and the front, the pair's ground speed across the front relative to
/// the front's own motion, and a weight (1 for the real front).
#[derive(Clone, Copy, Debug)]
struct Crossing {
    t: f64,
    theta: f64,
    sin_phi: f64,
    normal_speed: f64,
    weight: f64,
}

/// Times in `[from, to)` at which the pair's midpoint is at angular distance `target(t)` from
/// `source`, outside 10° of the source and of its antipode.
fn crossings(epochs: &[Epoch], source: V3, from: f64, to: f64, target: &dyn Fn(f64) -> f64, weight: &dyn Fn(f64) -> f64, front_speed: &dyn Fn(f64) -> f64) -> Vec<Crossing> {
    let (i, j) = (epochs.partition_point(|e| e.t < from), epochs.partition_point(|e| e.t < to));
    let theta = |e: &Epoch| dot(scale(e.mid, 1.0 / norm(e.mid)), source).clamp(-1.0, 1.0).acos();
    let mut out = Vec::new();
    for k in i.max(1)..j {
        let (a, b) = (&epochs[k - 1], &epochs[k]);
        if b.t - a.t > 1.5 * DT {
            continue;
        }
        let (ga, gb) = (theta(a) - target(a.t), theta(b) - target(b.t));
        if ga.signum() == gb.signum() {
            continue;
        }
        let th = theta(b);
        if !(10f64.to_radians()..170f64.to_radians()).contains(&th) {
            continue;
        }
        let t = a.t + (b.t - a.t) * ga / (ga - gb);
        let p = scale(b.mid, 1.0 / norm(b.mid));
        let away = sub(scale(p, dot(p, source)), source); // grows the angular distance
        let track = sub(b.vel, scale(p, dot(b.vel, p)));
        let outward = dot(away, track) / (norm(away) * norm(track)); // signed sin φ
        let ground = norm(track) * R_MEAN / norm(b.mid);
        let normal_speed = (ground * outward - front_speed(t)).abs();
        out.push(Crossing { t, theta: th, sin_phi: outward.abs(), normal_speed, weight: weight(t) });
    }
    out
}

/// Crossings of the real front, and of AOD1B's 3-hourly copies of it, in the window after `t0`.
fn all_crossings(epochs: &[Epoch], source: V3, t0: f64, span: (f64, f64)) -> (Vec<Crossing>, Vec<Crossing>) {
    // Past the antipode the front comes back toward the source, so its outward speed changes sign.
    let outward_speed = |t: f64| if (LAMB_SPEED * (t - t0) / R_MEAN).rem_euclid(2.0 * PI) > PI { -LAMB_SPEED } else { LAMB_SPEED };
    let (w0, w1) = (t0 + span.0, t0 + span.1);
    let real = crossings(epochs, source, w0, w1, &|t| front(t, t0), &|_| 1.0, &outward_speed);
    let mut ghosts = Vec::new();
    let first = ((w0 - AOD_ORIGIN) / AOD_STEP).ceil() as i64;
    let last = ((w1 - AOD_ORIGIN) / AOD_STEP).ceil() as i64;
    for k in first..=last {
        let tj = AOD_ORIGIN + k as f64 * AOD_STEP;
        if tj < t0 {
            continue; // AOD1B can only hold the pulse after it exists
        }
        let theta_j = front(tj, t0);
        let (from, to) = ((tj - AOD_STEP).max(w0), (tj + AOD_STEP).min(w1));
        ghosts.extend(crossings(epochs, source, from, to, &|_| theta_j, &|t| 1.0 - (t - tj).abs() / AOD_STEP, &|_| 0.0));
    }
    (real, ghosts)
}

/// Template series (range acceleration per hPa at θ = 90°) of a set of crossings, at `times`.
fn template_series(times: &[f64], crossings: &[Crossing], width_min: f64, cache: &mut std::collections::HashMap<(u32, u32, u32), Vec<f64>>) -> Vec<f64> {
    let mut out = vec![0.0; times.len()];
    for c in crossings {
        let deg = c.sin_phi.clamp(0.01, 1.0).asin().to_degrees().round().max(1.0) as u32;
        let speed = (c.normal_speed / 50.0).round().max(1.0) as u32; // 50 m/s bins
        let tmpl = cache.entry(((width_min * 100.0) as u32, deg, speed))
            .or_insert_with(|| tonga::template_moving(width_min * 60.0, f64::from(deg).to_radians(), f64::from(speed) * 50.0));
        let amp = c.weight / c.theta.sin().sqrt();
        let half = (tmpl.len() / 2) as f64;
        let start = times.partition_point(|&t| t < c.t - (half - 2.0) * DT);
        for (i, &t) in times.iter().enumerate().skip(start) {
            let x = (t - c.t) / DT + half;
            if x >= (tmpl.len() - 2) as f64 {
                break;
            }
            let (k, f) = (x.floor() as usize, x.fract());
            out[i] += amp * (tmpl[k] * (1.0 - f) + tmpl[k + 1] * f);
        }
    }
    out
}

/// The 4th-order band-pass, run zero-phase over each contiguous run of `times`.
fn band_pass(times: &[f64], x: &[f64]) -> Vec<f64> {
    let sections = butterworth_band(BAND.0, BAND.1, DT);
    let mut out = vec![0.0; x.len()];
    let mut start = 0;
    for i in 1..=x.len() {
        if i == x.len() || times[i] - times[i - 1] > 1.5 * DT {
            if i - start > 64 {
                out[start..i].copy_from_slice(&filtfilt(&sections, &x[start..i], 1500));
            }
            start = i;
        }
    }
    out
}

/// One-sided PSD (Welch, Hann segments of `len` samples zero-padded to `nfft`, half overlap) of
/// the contiguous runs of a band-passed series, as (frequency, PSD).
fn psd(times: &[f64], x: &[f64], len: usize, nfft: usize) -> Vec<(f64, f64)> {
    let w: Vec<f64> = (0..len).map(|i| 0.5 - 0.5 * (2.0 * PI * i as f64 / len as f64).cos()).collect();
    let norm2: f64 = w.iter().map(|v| v * v).sum();
    let mut acc = vec![0.0; nfft / 2];
    let mut count = 0usize;
    let mut run_start = 0;
    for i in 1..=x.len() {
        if i == x.len() || times[i] - times[i - 1] > 1.5 * DT {
            let run = &x[run_start..i];
            let mut s = 0;
            while s + len <= run.len() {
                let mut re = vec![0.0; nfft];
                for k in 0..len {
                    re[k] = run[s + k] * w[k];
                }
                let mut im = vec![0.0; nfft];
                gravimetry::noise::fft(&mut re, &mut im, false);
                for (k, a) in acc.iter_mut().enumerate() {
                    *a += 2.0 * DT * (re[k] * re[k] + im[k] * im[k]) / norm2;
                }
                count += 1;
                s += len / 2;
            }
            run_start = i;
        }
    }
    (1..nfft / 2).map(|k| (k as f64 / (nfft as f64 * DT), acc[k] / count.max(1) as f64)).collect()
}

/// The noise-weighted inner product `4 Re Σ ã b̃* Δf / S(f)` over the band, summed over the
/// contiguous runs of `times` within `[from, to)`.
fn inner(times: &[f64], from: f64, to: f64, series: &[&[f64]], s: &[(f64, f64)]) -> Vec<Vec<f64>> {
    let n = series.len();
    let mut g = vec![vec![0.0; n]; n];
    let (i0, i1) = (times.partition_point(|&t| t < from), times.partition_point(|&t| t < to));
    let mut start = i0;
    for i in i0 + 1..=i1 {
        if i == i1 || times[i] - times[i - 1] > 1.5 * DT {
            let len = i - start;
            if len >= 256 {
                let nfft = len.next_power_of_two();
                let spectra: Vec<(Vec<f64>, Vec<f64>)> = series.iter().map(|x| {
                    let mut re = vec![0.0; nfft];
                    re[..len].copy_from_slice(&x[start..i]);
                    let mut im = vec![0.0; nfft];
                    gravimetry::noise::fft(&mut re, &mut im, false);
                    (re, im)
                }).collect();
                let df = 1.0 / (nfft as f64 * DT);
                for k in 1..nfft / 2 {
                    let f = k as f64 * df;
                    if !(0.7e-3..14e-3).contains(&f) {
                        continue;
                    }
                    let j = s.partition_point(|p| p.0 < f).min(s.len() - 1);
                    let weight = 4.0 * DT * DT * df / s[j].1;
                    for a in 0..n {
                        for b in a..n {
                            let v = weight * (spectra[a].0[k] * spectra[b].0[k] + spectra[a].1[k] * spectra[b].1[k]);
                            g[a][b] += v;
                            if a != b {
                                g[b][a] += v;
                            }
                        }
                    }
                }
            }
            start = i;
        }
    }
    g
}

/// The test statistic for one origin: for each width, the generalised least-squares amplitude of
/// the real-front template (with AOD1B's copy as a second template), its σ, and Z = Â/σ, over the
/// crossings between `span.0` and `span.1` seconds after `t0`. Returns the width with the largest
/// Z: width, Â, σ_A, Z, crossings, and AOD1B's amplitude B̂ with its σ (0, 0 when it has no
/// crossings).
struct Context<'a> {
    times: &'a [f64],
    residual: &'a [f64],
    revs: &'a [(std::ops::Range<usize>, f64, f64)],
    epochs: &'a [Epoch],
    psd: &'a [(f64, f64)],
}

type Stat = (f64, f64, f64, f64, usize, f64, f64);

fn statistic(ctx: &Context, source: V3, t0: f64, span: (f64, f64), cache: &mut std::collections::HashMap<(u32, u32, u32), Vec<f64>>) -> Option<Stat> {
    let (real, ghosts) = all_crossings(ctx.epochs, source, t0, span);
    if real.is_empty() {
        return None;
    }
    let mut best: Option<Stat> = None;
    for width in WIDTHS_MIN {
        let mut tr = template_series(ctx.times, &real, width, cache);
        let mut ta = template_series(ctx.times, &ghosts, width, cache);
        detrend(ctx.times, &mut tr, ctx.revs);
        detrend(ctx.times, &mut ta, ctx.revs);
        let (tr, ta) = (band_pass(ctx.times, &tr), band_pass(ctx.times, &ta));
        let neg: Vec<f64> = ta.iter().map(|v| -v).collect();
        let g = inner(ctx.times, t0 + span.0, t0 + span.1, &[&tr, &neg, ctx.residual], ctx.psd);
        let (a, sigma, b, sigma_b) = if g[1][1] > 1e-6 * g[0][0] {
            let det = g[0][0] * g[1][1] - g[0][1] * g[0][1];
            let a = (g[1][1] * g[0][2] - g[0][1] * g[1][2]) / det;
            let b = (g[0][0] * g[1][2] - g[0][1] * g[0][2]) / det;
            (a, (g[1][1] / det).sqrt(), b, (g[0][0] / det).sqrt())
        } else {
            (g[0][2] / g[0][0], (1.0 / g[0][0]).sqrt(), 0.0, 0.0)
        };
        let z = a / sigma;
        if best.is_none_or(|s| z > s.3) {
            best = Some((width, a, sigma, z, real.len(), b, sigma_b));
        }
    }
    best
}

/// Stage 2: the residual series over `dates`, the null from fake eruptions before the real one,
/// and (only with [`UNBLIND`]) the real statistic.
pub fn test(l1b_dir: &Path, aod_dir: &Path, goco: &Path, costg: &Path, tides: &Path, dates: &[&str]) -> Result<String> {
    let (field, top) = static_field(goco, costg)?;
    let mut model = Model::new(&field, &all_aod1b(aod_dir)?);
    model.tide = Some(OceanTide::load(tides)?);
    let (c, d, obs) = load_days(l1b_dir, dates)?;
    let mut nodes = ascending_nodes(&c);
    nodes.sort_by(f64::total_cmp);
    let (fit_c, _) = fitted_orbits(&model, &c.states, &nodes);
    let (fit_d, _) = fitted_orbits(&model, &d.states, &nodes);
    let inp = Inputs { model: &model, c: &fit_c, d: &fit_d, c_acc: &c.accel_x, d_acc: &d.accel_x };
    let epochs = residuals(&inp, &obs);
    let times: Vec<f64> = epochs.iter().map(|e| e.t).collect();
    let revs = revolutions(&times, &nodes);
    let mut residual: Vec<f64> = epochs.iter().map(|e| e.residual).collect();
    detrend(&times, &mut residual, &revs);
    let residual = band_pass(&times, &residual);
    // The noise spectrum comes from the residuals before the eruption only.
    let pre = times.partition_point(|&t| t < ERUPTION);
    let spectrum = psd(&times[..pre], &residual[..pre], 10_800, 16_384);
    let ctx = Context { times: &times, residual: &residual, revs: &revs, epochs: &epochs, psd: &spectrum };

    let mut report = String::new();
    writeln!(report, "T1 stage 2: the pre-registered test")?;
    writeln!(report, "Model: COST-G January 2022 to degree {top}, GOCO06s to {NMAX}; AOD1B RL06 GLO; GOT4.10c ocean tides to degree {TIDE_NMAX}; solid tides; Sun and Moon; fitted orbits; moving-front templates.")?;
    writeln!(report, "Residual epochs: {} ({} before the eruption).", epochs.len(), pre)?;

    let mut cache = std::collections::HashMap::new();
    let tonga = unit_vector(TONGA.0, TONGA.1);
    let mut null = Vec::new();
    let mut t0 = AOD_ORIGIN;
    while t0 + WINDOW <= ERUPTION {
        if let Some(r) = statistic(&ctx, tonga, t0, (0.0, WINDOW), &mut cache) {
            null.push(r);
        }
        t0 += 600.0;
    }
    let tonga_fakes = null.len();
    let rings = crate::shared::land_rings();
    let mut rng = crate::shared::Rng::new(2022);
    let mut uniform = || (rng.next_u64() >> 11) as f64 / (1u64 << 53) as f64;
    let mut sources = 0;
    while sources < 30 {
        let (lat, lon) = ((2.0 * uniform() - 1.0).asin().to_degrees(), 360.0 * uniform() - 180.0);
        if crate::shared::on_land(&rings, lon, lat) {
            continue;
        }
        let t0 = AOD_ORIGIN + uniform() * (ERUPTION - WINDOW - AOD_ORIGIN);
        if let Some(r) = statistic(&ctx, unit_vector(lat, lon), t0, (0.0, WINDOW), &mut cache) {
            null.push(r);
        }
        sources += 1;
    }
    let zs: Vec<f64> = null.iter().map(|r| r.3).collect();
    let threshold = crate::shared::quantile(&zs, 0.99);
    let sigmas: Vec<f64> = null.iter().map(|r| r.2).collect();
    writeln!(report, "Null: {} fake eruptions ({tonga_fakes} from Tonga every 10 min, {} from random ocean points).", null.len(), null.len() - tonga_fakes)?;
    writeln!(report, "  null Z: median {:.2}, 1st–99th percentile {:.2} to {:.2}; detection threshold Z > {threshold:.2}",
        crate::shared::quantile(&zs, 0.5), crate::shared::quantile(&zs, 0.01), threshold)?;
    writeln!(report, "  sensitivity: median σ of the amplitude {:.2} hPa (5–95%: {:.2}–{:.2}); crossings per window: median {}",
        crate::shared::quantile(&sigmas, 0.5), crate::shared::quantile(&sigmas, 0.05), crate::shared::quantile(&sigmas, 0.95),
        crate::shared::quantile(&null.iter().map(|r| r.4 as f64).collect::<Vec<_>>(), 0.5))?;
    if UNBLIND.load(std::sync::atomic::Ordering::Relaxed) {
        let (width, a, sigma, z, n, b, sigma_b) = statistic(&ctx, tonga, ERUPTION, (0.0, WINDOW), &mut cache).ok_or("no crossings after the eruption")?;
        let rank = zs.iter().filter(|&&v| v >= z).count();
        writeln!(report, "\nREAL: {n} crossings; best width σ = {width} min; amplitude {a:.2} ± {sigma:.2} hPa; Z = {z:.2}; p = {:.3} ({rank} of {} null values ≥ Z)",
            (rank + 1) as f64 / (zs.len() + 1) as f64, zs.len())?;
        writeln!(report, "Decision: {}", if z > threshold && a > 0.0 { "DETECTED" } else { "not detected" })?;
        writeln!(report, "95% upper limit on the amplitude: {:.2} hPa", a.max(0.0) + 1.645 * sigma)?;
        writeln!(report, "Intervals: A = {a:.2} hPa, 68% [{:.2}, {:.2}], 95% [{:.2}, {:.2}]; AOD1B's copy B = {b:.2} ± {sigma_b:.2} hPa",
            a - sigma, a + sigma, a - 1.96 * sigma, a + 1.96 * sigma)?;
        // Secondary (not part of the decision): the second passage, 24 h to the end of the LRI data.
        if let Some((w2, a2, s2, z2, n2, _, _)) = statistic(&ctx, tonga, ERUPTION, (WINDOW, 43.75 * 3600.0), &mut cache) {
            writeln!(report, "Secondary, second passage (24–43.75 h): {n2} crossings; width {w2} min; amplitude {a2:.2} ± {s2:.2} hPa; Z = {z2:.2}")?;
        }
    } else {
        writeln!(report, "\n(blind: the real window is not evaluated)")?;
    }
    Ok(report)
}

/// A labelled residual spectrum and its band-passed rms.
type Curve<'a> = (&'a str, Vec<(f64, f64)>, f64);

fn figure(curves: &[Curve]) -> String {
    let (w, h, left, right, top, bottom) = (960.0, 560.0, 92.0, 900.0, 110.0, 500.0);
    let (f_lo, f_hi, a_lo, a_hi) = (1e-4, 0.25, 1e-12, 1e-6);
    let x = |f: f64| left + (right - left) * (f / f_lo).log10() / (f_hi / f_lo).log10();
    let y = |a: f64| top + (bottom - top) * (a_hi / a.clamp(a_lo, a_hi)).log10() / (a_hi / a_lo).log10();
    let mut s = open(w, h, "T1 noise before the eruption");
    let _ = write!(s, r#"<text x="24" y="34" class="title">T1: GRACE-FO laser-ranging residuals before the eruption</text>
<text x="24" y="56" class="subtitle">Residual range acceleration on 13–14 January 2022 (blind), against the instrument noise and the Lamb-wave signal it must reveal. {}</text>
"#, esc(&curves.iter().map(|c| format!("{}: {:.1e} m/s² rms in band", c.0, c.2)).collect::<Vec<_>>().join("; ")));
    for e in -12..=-6 {
        let yy = y(10f64.powi(e));
        let _ = writeln!(s, r#"<line x1="{left}" y1="{yy:.1}" x2="{right}" y2="{yy:.1}" class="grid"/><text x="{:.1}" y="{:.1}" class="tick" text-anchor="end">{}</text>"#, left - 8.0, yy + 4.0, pow10(e));
    }
    for (f, label) in [(1e-4, "0.1 mHz"), (1e-3, "1 mHz"), (1e-2, "10 mHz"), (1e-1, "0.1 Hz")] {
        let xx = x(f);
        let _ = writeln!(s, r#"<line x1="{xx:.1}" y1="{top}" x2="{xx:.1}" y2="{bottom}" class="grid"/><text x="{xx:.1}" y="{:.1}" class="tick" text-anchor="middle">{label}</text>"#, bottom + 18.0);
    }
    let band = (x(BAND.0), x(BAND.1));
    let _ = writeln!(s, r#"<rect x="{:.1}" y="{top}" width="{:.1}" height="{:.1}" fill="var(--s3)" opacity="0.10"/><text x="{:.1}" y="{:.1}" class="tick">test band</text>"#, band.0, band.1 - band.0, bottom - top, band.0 + 6.0, top + 16.0);
    let path = |pts: &mut dyn Iterator<Item = (f64, f64)>| -> String {
        pts.filter(|p| p.0 >= f_lo && p.0 <= f_hi).enumerate().map(|(j, (f, a))| format!("{}{:.1},{:.1}", if j == 0 { 'M' } else { 'L' }, x(f), y(a))).collect()
    };
    let nominal = path(&mut (0..=200).map(|i| { let f = f_lo * (f_hi / f_lo).powf(i as f64 / 200.0); (f, pair_noise(Accelerometer::SuperStar, Ranging::GraceFo, f)) }));
    let lamb = path(&mut tonga::spectrum(&tonga::template(5.0 * 60.0, PI / 2.0)).into_iter().map(|(f, a)| (f, 2.0 * f.sqrt() * a)).filter(|p| p.1 > a_lo));
    let mut series: Vec<(String, &str, String)> = curves.iter().zip(["s4", "s1"]).map(|(c, class)| (path(&mut c.1.iter().copied()), class, format!("residual, {}", c.0))).collect();
    series.push((nominal, "ink", "instrument noise (nominal)".into()));
    series.push((lamb, "s2", "Lamb wave per crossing, 1 hPa".into()));
    for (d, class, label) in &series {
        let _ = writeln!(s, r#"<path d="{d}" class="line {class}"><title>{}</title></path>"#, esc(label));
    }
    let mut lx = 24.0;
    for (_, class, label) in &series {
        let label = label.as_str();
        let _ = writeln!(s, r#"<line x1="{lx:.1}" y1="84" x2="{:.1}" y2="84" class="line {class}"/><text x="{:.1}" y="88">{label}</text>"#, lx + 16.0, lx + 22.0);
        lx += 22.0 + 6.5 * label.chars().count() as f64 + 24.0;
    }
    s.push_str("</svg>\n");
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gps_time_of_the_eruption() {
        // 2022-01-15 04:14:45 UTC is 04:15:03 GPS.
        assert_eq!(gps_seconds("2022-01-15", "04:15:03").unwrap(), ERUPTION);
        assert_eq!(gps_seconds("2022-01-13", "00:00:00").unwrap(), 695_304_000.0);
    }

    #[test]
    fn band_pass_passes_the_band_only() {
        let sections = butterworth_band(1e-3, 10e-3, DT);
        let gain = |f: f64| {
            let x: Vec<f64> = (0..40_000).map(|i| (2.0 * PI * f * i as f64 * DT).sin()).collect();
            let y = filtfilt(&sections, &x, 1500);
            rms(y[10_000..30_000].iter().copied()) / rms(x[10_000..30_000].iter().copied())
        };
        assert!((gain(3e-3) - 1.0).abs() < 0.02, "{}", gain(3e-3));
        assert!(gain(0.2e-3) < 1e-4 && gain(50e-3) < 1e-4, "{} {}", gain(0.2e-3), gain(50e-3));
        // Zero-phase: each edge frequency passes at −6 dB (−3 dB each way).
        assert!((gain(1e-3) - 0.5).abs() < 0.03 && (gain(10e-3) - 0.5).abs() < 0.03, "{} {}", gain(1e-3), gain(10e-3));
    }

    #[test]
    fn welch_recovers_white_noise() {
        let mut rng = crate::shared::Rng::new(7);
        let mut uniform = || (rng.next_u64() >> 11) as f64 / (1u64 << 53) as f64 + 1e-300;
        let sigma = 1e-9;
        let x: Vec<f64> = (0..200_000).map(|_| sigma * (-2.0 * uniform().ln()).sqrt() * (2.0 * PI * uniform()).cos()).collect();
        let asd = welch(&[&x], 2048);
        let mean = asd.iter().map(|p| p.1).sum::<f64>() / asd.len() as f64;
        // White noise of variance σ² sampled every DT has one-sided ASD σ√(2·DT).
        assert!((mean / (sigma * (2.0 * DT).sqrt()) - 1.0).abs() < 0.03, "{mean:e}");
    }

    #[test]
    fn integrator_keeps_the_jacobi_constant() {
        // A point-mass orbit at 490 km in the rotating frame conserves v²/2 − Ω²(x²+y²)/2 − GM/r.
        let r0 = 6.868e6;
        let speed = (GM_E / r0).sqrt();
        let (x0, v0) = ([r0, 0.0, 0.0], [0.0, speed * 0.3 - OMEGA * r0, speed * 0.95]);
        let jacobi = |(x, v): (V3, V3)| dot(v, v) / 2.0 - OMEGA * OMEGA * (x[0] * x[0] + x[1] * x[1]) / 2.0 - GM_E / norm(x);
        let states = integrate(0.0, 2900, x0, v0, &mut |_, x, v| add(scale(x, -GM_E / norm(x).powi(3)), fictitious(x, v)));
        let (first, last) = (jacobi(states[0]), jacobi(states[2900]));
        assert!(((last - first) / first).abs() < 1e-11, "{first} {last}");
    }

    #[test]
    fn third_body_and_ephemerides() {
        // At the eruption the Sun is ~1 AU away, the Moon ~360–410 thousand km.
        let (sun, moon) = sun_moon(ERUPTION);
        assert!((norm(sun) / 1.496e11 - 1.0).abs() < 0.02 && (3.5e8..4.1e8).contains(&norm(moon)));
        // The Moon's tidal pull at the sub-lunar point at 490 km altitude: 2GM r/d³ ≈ 1.1e-6 m/s².
        let x = scale(moon, 6.868e6 / norm(moon));
        let a = norm(third_body(x, moon, GM_MOON));
        assert!((a / (2.0 * GM_MOON * 6.868e6 / norm(moon).powi(3)) - 1.0).abs() < 0.05, "{a:e}");
    }
}
