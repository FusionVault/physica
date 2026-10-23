//! ICGEM `.gfc` gravity-field models (format spec of February 2023), static or time-variable,
//! icgem1.0 or icgem2.0.
//!
//! Real files bend the spec, so the reader keys on `end_of_head`, ignores header lines it does not
//! know (COST-G `CMMNT`, Swarm `WARNING:`, free text), accepts any `errors` value (COST-G writes
//! `empirical`), and reads time columns from the end of each record so the number of error columns
//! never matters.

use std::collections::HashMap;
use std::fmt;

use geopotential::{Coefficients, Term, TermKind, TideSystem, TimeVariable};

/// A parsed `.gfc` file.
#[derive(Clone, Debug, PartialEq)]
pub struct Gfc {
    /// `modelname` or `product_name` from the header.
    pub name: Option<String>,
    /// The field. A static file has no `terms`, so `model.base` is the whole field; evaluate a
    /// time-variable one with [`TimeVariable::at`] at a Julian epoch (see [`julian_epoch`]).
    pub model: TimeVariable,
    /// Standard deviations of the `gfc`/`gfct` coefficients, when `errors` is not `no` (the first
    /// pair of error columns: calibrated if the file has both).
    pub sigma: Option<Coefficients>,
}

/// Why a `.gfc` file could not be read. Line numbers are 1-based.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum GfcError {
    /// A mandatory header key (`earth_gravity_constant`, `radius`, `max_degree`, `end_of_head`)
    /// is missing or unparsable.
    MissingHeader(&'static str),
    /// The file uses something this reader does not support.
    Unsupported(String),
    /// A data record is malformed.
    BadRecord {
        /// 1-based line number.
        line: usize,
        /// What is wrong with it.
        reason: &'static str,
    },
}

impl fmt::Display for GfcError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            GfcError::MissingHeader(key) => write!(f, "missing or invalid .gfc header key `{key}`"),
            GfcError::Unsupported(what) => write!(f, "unsupported .gfc content: {what}"),
            GfcError::BadRecord { line, reason } => write!(f, ".gfc line {line}: {reason}"),
        }
    }
}

impl std::error::Error for GfcError {}

/// The Julian epoch (years of 365.25 days from J2000.0, TT ignored) of a calendar date and time.
/// ICGEM time-variable terms use `(t − t0)` in these years: GOCO06s defines its annual period as
/// 365.25 d.
pub fn julian_epoch(year: i32, month: u32, day: u32, hour: u32, minute: u32) -> f64 {
    // Days from 1970-01-01 in the proleptic Gregorian calendar (H. Hinnant's days_from_civil).
    let y = i64::from(year) - i64::from(month <= 2);
    let era = y.div_euclid(400);
    let yoe = y - era * 400;
    let mp = (i64::from(month) + 9) % 12;
    let doy = (153 * mp + 2) / 5 + i64::from(day) - 1;
    let days = era * 146_097 + yoe * 365 + yoe / 4 - yoe / 100 + doy - 719_468;
    let mjd = days as f64 + 40_587.0 + (f64::from(hour) + f64::from(minute) / 60.0) / 24.0;
    2000.0 + (mjd - 51_544.5) / 365.25
}

/// Parse `yyyymmdd` or `yyyymmdd.hhmm` into a Julian epoch.
fn parse_date(s: &str) -> Option<f64> {
    let (date, time) = s.split_once('.').unwrap_or((s, "0000"));
    if date.len() != 8 || !date.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    let time = format!("{time:0<4}");
    let (year, month, day) = (date[..4].parse().ok()?, date[4..6].parse().ok()?, date[6..].parse().ok()?);
    let (hour, minute) = (time.get(..2)?.parse().ok()?, time.get(2..4)?.parse().ok()?);
    let valid = (1..=12).contains(&month) && (1..=31).contains(&day) && hour < 24 && minute < 60;
    valid.then(|| julian_epoch(year, month, day, hour, minute))
}

/// Parse a number, accepting Fortran `D` exponents.
fn number(s: &str) -> Option<f64> {
    s.replace(['D', 'd'], "e").parse().ok()
}

/// Read an ICGEM `.gfc` file.
///
/// ```
/// let text = "earth_gravity_constant 3.986004415e14\nradius 6378136.3\nmax_degree 2\n\
///             end_of_head\ngfc 0 0 1.0 0.0\ngfc 2 0 -4.84e-4 0.0\n";
/// let gfc = gravity_field::read_gfc(text).unwrap();
/// assert_eq!(gfc.model.base.c(2, 0), -4.84e-4);
/// ```
pub fn read_gfc(text: &str) -> Result<Gfc, GfcError> {
    let mut lines = text.lines().enumerate();
    let (mut gm, mut radius, mut nmax) = (None, None, None);
    let (mut tide, mut name, mut has_errors, mut icgem2) = (TideSystem::Unknown, None, false, false);
    let mut header_done = false;
    for (_, line) in lines.by_ref() {
        let mut fields = line.split_whitespace();
        let (Some(key), value) = (fields.next(), fields.next()) else { continue };
        match key {
            "end_of_head" => {
                header_done = true;
                break;
            }
            "earth_gravity_constant" => gm = value.and_then(number),
            "radius" => radius = value.and_then(number),
            "max_degree" => nmax = value.and_then(|v| v.parse::<usize>().ok()),
            "modelname" | "product_name" => name = value.map(str::to_owned),
            "errors" => has_errors = value != Some("no"),
            "format" => icgem2 = value == Some("icgem2.0"),
            "tide_system" => {
                tide = match value {
                    Some("zero_tide") => TideSystem::ZeroTide,
                    Some("tide_free") => TideSystem::TideFree,
                    Some("mean_tide") => TideSystem::MeanTide,
                    _ => TideSystem::Unknown,
                }
            }
            "norm" if value == Some("unnormalized") => {
                // shortcut: unnormalized models are rare; convert them when one is needed.
                return Err(GfcError::Unsupported("norm unnormalized".into()));
            }
            _ => {}
        }
    }
    if !header_done {
        return Err(GfcError::MissingHeader("end_of_head"));
    }
    let gm = gm.ok_or(GfcError::MissingHeader("earth_gravity_constant"))?;
    let radius = radius.ok_or(GfcError::MissingHeader("radius"))?;
    let nmax = nmax.ok_or(GfcError::MissingHeader("max_degree"))?;
    if nmax > geopotential::MAX_DEGREE {
        // shortcut: synthesis stops at degree 2700; truncate on read when such models are needed.
        return Err(GfcError::Unsupported(format!("max_degree {nmax} > {}", geopotential::MAX_DEGREE)));
    }

    let mut base = Coefficients::zeros(gm, radius, nmax);
    base.tide_system = tide;
    let mut sigma = has_errors.then(|| base.clone());
    let mut terms = Vec::new();
    let mut epochs: HashMap<(usize, usize), f64> = HashMap::new(); // icgem1.0 gfct reference epochs

    for (i, line) in lines {
        let f: Vec<&str> = line.split_whitespace().collect();
        let Some(&key) = f.first() else { continue };
        let bad = |reason| GfcError::BadRecord { line: i + 1, reason };
        let (n, m) = match (f.get(1).and_then(|s| s.parse().ok()), f.get(2).and_then(|s| s.parse().ok())) {
            (Some(n), Some(m)) if m <= n && n <= nmax => (n, m),
            (Some(_), Some(_)) => return Err(bad("degree/order outside max_degree")),
            _ => return Err(bad("bad degree/order")),
        };
        let (c, s) = match (f.get(3).copied().and_then(number), f.get(4).copied().and_then(number)) {
            (Some(c), Some(s)) => (c, s),
            _ => return Err(bad("bad coefficient")),
        };
        // Time columns sit at the end of the record: [t0 [t1]] [period].
        let from_end = |k: usize| f.len().checked_sub(k).filter(|&j| j >= 5).map(|j| f[j]);
        let date = |k| from_end(k).and_then(parse_date).ok_or_else(|| bad("bad or missing date"));
        let period = || from_end(1).and_then(number).filter(|p| *p > 0.0).ok_or_else(|| bad("bad period"));
        let mut term = |kind, epoch: f64, valid| terms.push(Term { n, m, kind, c, s, epoch, valid });
        let unbounded = f64::NEG_INFINITY..f64::INFINITY;
        let gfct_epoch = || epochs.get(&(n, m)).copied().ok_or_else(|| bad("record before its gfct"));

        match (key, icgem2) {
            ("gfc", _) => base.set(n, m, c, s),
            ("gfct", false) => {
                epochs.insert((n, m), date(1)?);
                base.set(n, m, c, s);
            }
            ("gfct", true) => {
                let (t0, t1) = (date(2)?, date(1)?);
                term(TermKind::Value, t0, t0..t1);
            }
            ("trnd" | "dot", false) => term(TermKind::Trend, gfct_epoch()?, unbounded),
            ("trnd" | "dot", true) => {
                let (t0, t1) = (date(2)?, date(1)?);
                term(TermKind::Trend, t0, t0..t1);
            }
            ("acos" | "asin", false) => {
                let period = period()?;
                let kind = if key == "acos" { TermKind::Cos { period } } else { TermKind::Sin { period } };
                term(kind, gfct_epoch()?, unbounded);
            }
            ("acos" | "asin", true) => {
                let (period, t0, t1) = (period()?, date(3)?, date(2)?);
                let kind = if key == "acos" { TermKind::Cos { period } } else { TermKind::Sin { period } };
                term(kind, t0, t0..t1);
            }
            _ => return Err(bad("unknown record key")),
        }
        // Error columns follow C, S, and must leave room for the record's own date columns.
        let date_cols = match (key, icgem2) {
            ("gfc", _) => Some(0),
            ("gfct", false) => Some(1),
            ("gfct", true) => Some(2),
            _ => None,
        };
        if let (Some(sigma), Some(date_cols)) = (sigma.as_mut(), date_cols) {
            if f.len() >= 7 + date_cols {
                match (number(f[5]), number(f[6])) {
                    (Some(sc), Some(ss)) => sigma.set(n, m, sc, ss),
                    _ => return Err(bad("bad error column")),
                }
            }
        }
    }
    Ok(Gfc { name, model: TimeVariable { base, terms }, sigma })
}

#[cfg(test)]
mod tests {
    use super::*;

    const COSTG: &str = include_str!("../tests/data/costg_2026-07_n10.gfc");
    const SWARM: &str = include_str!("../tests/data/swarm_2025-12_n10.gfc");
    const GOCO: &str = include_str!("../tests/data/goco06s_n4.gfc");

    #[test]
    fn julian_epochs() {
        assert_eq!(julian_epoch(2000, 1, 1, 12, 0), 2000.0);
        // GOCO06s states its reference epoch 2010-01-01 is MJD 55197.
        assert_eq!(julian_epoch(2010, 1, 1, 0, 0), 2000.0 + (55_197.0 - 51_544.5) / 365.25);
        assert_eq!(parse_date("20100101"), Some(julian_epoch(2010, 1, 1, 0, 0)));
        assert_eq!(parse_date("20020815.0817"), Some(julian_epoch(2002, 8, 15, 8, 17)));
        assert_eq!(parse_date("6.658e-13"), None);
        assert_eq!(parse_date("20101301"), None);
    }

    /// COST-G: no begin_of_head, CMMNT lines, `errors empirical`, `+…E+01` numbers.
    #[test]
    fn costg_monthly() {
        let g = read_gfc(COSTG).unwrap();
        assert_eq!(g.name.as_deref(), Some("GSM-2_2026182-2026212_GRFO_COSTG_BF01_0201"));
        let b = &g.model.base;
        assert_eq!((b.gm, b.radius, b.max_degree()), (3.986004415e14, 6378136.3, 90));
        assert_eq!(b.tide_system, TideSystem::TideFree);
        assert_eq!((b.c(0, 0), b.c(2, 0)), (1.0, -4.84165307066e-4));
        assert_eq!(b.s(2, 2), -1.40038200535e-6);
        assert!(g.model.terms.is_empty());
        let sigma = g.sigma.unwrap();
        assert_eq!((sigma.c(2, 1), sigma.s(2, 1)), (2.45322715392e-11, 1.62793405520e-11));
    }

    /// Swarm: `comment` lines, a stray `WARNING:` line, `errors no`.
    #[test]
    fn swarm_monthly() {
        let g = read_gfc(SWARM).unwrap();
        let b = &g.model.base;
        assert_eq!(b.max_degree(), 40);
        assert_eq!((b.c(2, 1), b.s(2, 1)), (-8.08481517811e-10, 1.48254384135e-9));
        assert!(g.sigma.is_none());
    }

    /// GOCO06s (icgem1.0): free text before the header, gfct with t0, then trnd/acos/asin.
    #[test]
    fn goco06s_time_variable() {
        let g = read_gfc(GOCO).unwrap();
        assert_eq!(g.name.as_deref(), Some("GOCO06s"));
        assert_eq!(g.model.base.tide_system, TideSystem::ZeroTide);
        assert_eq!(g.model.terms.len(), 3 * (15 - 3)); // trnd/acos/asin for every (n,m), 2 ≤ n ≤ 4
        assert_eq!(g.sigma.as_ref().unwrap().c(2, 0), 6.658069597318e-13);
        let t0 = julian_epoch(2010, 1, 1, 0, 0);
        let t = julian_epoch(2026, 7, 16, 12, 0);
        let w = core::f64::consts::TAU * (t - t0);
        let c20 = -4.841694947139e-4 - 2.154662425317e-11 * (t - t0)
            + 8.331958055607e-11 * w.cos()
            + 7.523339889891e-11 * w.sin();
        assert!((g.model.at(t).c(2, 0) - c20).abs() < 1e-18);
        // At t0 only the gfct value and the acos amplitude remain.
        let acos44 = g.model.terms.iter().find(|x| (x.n, x.m, x.kind) == (4, 4, TermKind::Cos { period: 1.0 }));
        let expect = -1.884945995629e-7 + acos44.unwrap().c;
        assert!((g.model.at(t0).c(4, 4) - expect).abs() < 1e-22);
    }

    #[test]
    fn icgem2_intervals() {
        let text = "earth_gravity_constant 1\nradius 1\nmax_degree 2\nformat icgem2.0\nerrors no\nend_of_head\n\
            gfct 2 0 1.0 0.0 20000101.0000 20100101.0000\n\
            gfct 2 0 2.0 0.0 20100101.0000 20200101.0000\n\
            trnd 2 0 0.5 0.0 20100101.0000 20200101.0000\n\
            acos 2 0 0.25 0.0 20100101.0000 20200101.0000 1.0\n";
        let m = read_gfc(text).unwrap().model;
        assert_eq!(m.at(julian_epoch(2005, 1, 1, 0, 0)).c(2, 0), 1.0);
        let (t0, t) = (julian_epoch(2010, 1, 1, 0, 0), julian_epoch(2012, 1, 1, 0, 0));
        let expect = 2.0 + 0.5 * (t - t0) + 0.25 * (core::f64::consts::TAU * (t - t0)).cos();
        assert!((m.at(t).c(2, 0) - expect).abs() < 1e-15);
        assert_eq!(m.at(julian_epoch(2021, 1, 1, 0, 0)).c(2, 0), 0.0);
    }

    #[test]
    fn missing_error_columns_are_not_read_as_dates() {
        let text = "earth_gravity_constant 1\nradius 1\nmax_degree 2\nerrors formal\nend_of_head\n\
            gfct 2 0 1.0 0.0 20100101\ngfc 2 1 1.0 2.0 3.0 4.0\n";
        let sigma = read_gfc(text).unwrap().sigma.unwrap();
        assert_eq!((sigma.c(2, 0), sigma.c(2, 1), sigma.s(2, 1)), (0.0, 3.0, 4.0));
    }

    #[test]
    fn rejects_bad_input() {
        let head = "earth_gravity_constant 1\nradius 1\nmax_degree 2\nend_of_head\n";
        let err = |body: &str| read_gfc(&format!("{head}{body}")).unwrap_err();
        assert_eq!(err("gfc 3 0 1 0\n"), GfcError::BadRecord { line: 5, reason: "degree/order outside max_degree" });
        assert_eq!(err("gfc 1 2 1 0\n"), GfcError::BadRecord { line: 5, reason: "degree/order outside max_degree" });
        assert_eq!(err("gfc 2 0 x 0\n"), GfcError::BadRecord { line: 5, reason: "bad coefficient" });
        assert_eq!(err("trnd 2 0 1 0\n"), GfcError::BadRecord { line: 5, reason: "record before its gfct" });
        assert_eq!(err("gfct 2 0 1 0 6.6e-13\n"), GfcError::BadRecord { line: 5, reason: "bad or missing date" });
        assert_eq!(err("xyz 2 0 1 0\n"), GfcError::BadRecord { line: 5, reason: "unknown record key" });
        assert_eq!(read_gfc("radius 1\nmax_degree 2\nend_of_head\n"), Err(GfcError::MissingHeader("earth_gravity_constant")));
        assert_eq!(read_gfc("gfc 0 0 1 0\n"), Err(GfcError::MissingHeader("end_of_head")));
        assert!(matches!(read_gfc("norm unnormalized\nend_of_head\n"), Err(GfcError::Unsupported(_))));
    }
}
