//! GRACE / GRACE-FO technical-note replacement series: TN-14 (`C̄20`, `C̄30` from satellite laser
//! ranging) and TN-13 (degree-1 / geocentre coefficients). Both are anonymous downloads from the
//! PO.DAAC docs bucket. Epochs are Julian years, as for [`read_gfc`](crate::read_gfc).

use std::fmt;

use crate::gfc::parse_date;

/// Why a technical-note file could not be read. Line numbers are 1-based.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TechNoteError {
    /// 1-based line number, 0 if the problem is the file as a whole.
    pub line: usize,
    /// What is wrong.
    pub reason: &'static str,
}

impl fmt::Display for TechNoteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "technical note line {}: {}", self.line, self.reason)
    }
}

impl std::error::Error for TechNoteError {}

/// One TN-14 month: SLR replacement values for the GRACE month starting at `start`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SlrC20C30 {
    /// Start of the solution span, Julian years.
    pub start: f64,
    /// End of the solution span, Julian years.
    pub end: f64,
    /// Replacement `C̄20`.
    pub c20: f64,
    /// Replacement `C̄30`, where TN-14 gives one (it recommends using it from August 2016).
    pub c30: Option<f64>,
}

/// One TN-13 month: degree-1 coefficients for the month starting at `start`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Degree1 {
    /// Start of the month, Julian years.
    pub start: f64,
    /// End of the month, Julian years.
    pub end: f64,
    /// `C̄10`.
    pub c10: f64,
    /// `C̄11`.
    pub c11: f64,
    /// `S̄11`.
    pub s11: f64,
}

fn mjd_to_julian(mjd: f64) -> f64 {
    2000.0 + (mjd - 51_544.5) / 365.25
}

/// Read TN-14 (`TN-14_C30_C20_GSFC_SLR.txt`): the rows after the `Product:` line, columns
/// MJD start, year, C20, ΔC20, σ, C30, ΔC30, σ, MJD end, year.
pub fn read_tn14(text: &str) -> Result<Vec<SlrC20C30>, TechNoteError> {
    let mut lines = text.lines().enumerate().skip_while(|(_, l)| !l.starts_with("Product:"));
    if lines.next().is_none() {
        return Err(TechNoteError { line: 0, reason: "no `Product:` line" });
    }
    let mut out = Vec::new();
    for (i, line) in lines {
        let f: Vec<f64> = match line.split_whitespace().map(str::parse).collect() {
            Ok(f) => f,
            Err(_) => return Err(TechNoteError { line: i + 1, reason: "non-numeric field" }),
        };
        match f.as_slice() {
            [] => continue,
            &[start, _, c20, _, _, c30, _, _, end, _] => out.push(SlrC20C30 {
                start: mjd_to_julian(start),
                end: mjd_to_julian(end),
                c20,
                c30: (!c30.is_nan()).then_some(c30),
            }),
            _ => return Err(TechNoteError { line: i + 1, reason: "expected 10 columns" }),
        }
    }
    Ok(out)
}

/// Read TN-13 (`TN-13_GEOC_*.txt`): `GRCOF2 1 m C S σC σS start end` records after
/// `end of header`, one `(1,0)` and one `(1,1)` per month.
pub fn read_tn13(text: &str) -> Result<Vec<Degree1>, TechNoteError> {
    let mut lines = text.lines().enumerate().skip_while(|(_, l)| !l.starts_with("end of header"));
    if lines.next().is_none() {
        return Err(TechNoteError { line: 0, reason: "no `end of header` line" });
    }
    let mut out: Vec<Degree1> = Vec::new();
    for (i, line) in lines {
        let f: Vec<&str> = line.split_whitespace().collect();
        if f.is_empty() {
            continue;
        }
        let bad = |reason| TechNoteError { line: i + 1, reason };
        let [_, n, m, c, s, _, _, t0, t1] = f[..] else { return Err(bad("expected 9 columns")) };
        let num = |v: &str| v.parse::<f64>().map_err(|_| bad("bad number"));
        let date = |v: &str| parse_date(v).ok_or_else(|| bad("bad date"));
        let (c, s, start, end) = (num(c)?, num(s)?, date(t0)?, date(t1)?);
        match (n, m) {
            ("1", "0") => out.push(Degree1 { start, end, c10: c, c11: f64::NAN, s11: f64::NAN }),
            ("1", "1") => match out.last_mut() {
                Some(d) if d.start == start && d.c11.is_nan() => (d.c11, d.s11) = (c, s),
                _ => return Err(bad("(1,1) record without its (1,0)")),
            },
            _ => return Err(bad("not a degree-1 record")),
        }
    }
    if out.iter().any(|d| d.c11.is_nan()) {
        return Err(TechNoteError { line: 0, reason: "a month has no (1,1) record" });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::julian_epoch;

    #[test]
    fn tn14_rows() {
        let rows = read_tn14(include_str!("../tests/data/tn14_excerpt.txt")).unwrap();
        assert_eq!(rows.len(), 4);
        assert_eq!(rows[0].start, mjd_to_julian(52368.0));
        assert_eq!((rows[0].c20, rows[0].c30), (-4.8416934147454e-4, None));
        let last = rows[3];
        assert_eq!((last.c20, last.c30), (-4.8416986716357e-4, Some(9.5700367112704e-7)));
        assert_eq!(last.start, julian_epoch(2026, 8, 1, 0, 0)); // MJD 61253
        assert!(read_tn14("no product line\n").is_err());
    }

    #[test]
    fn tn13_months() {
        let months = read_tn13(include_str!("../tests/data/tn13_excerpt.txt")).unwrap();
        assert_eq!(months.len(), 2);
        assert_eq!((months[0].c10, months[0].c11, months[0].s11), (5.120437146e-10, 2.897447347e-10, -6.826708442e-11));
        assert_eq!(months[1].start, julian_epoch(2026, 8, 1, 0, 0));
        let orphan = "end of header\nGRCOF2 1 1 1e-10 1e-10 0 0 20260801.0000 20260901.0000\n";
        assert_eq!(read_tn13(orphan).unwrap_err().reason, "(1,1) record without its (1,0)");
    }
}
