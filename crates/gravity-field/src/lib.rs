//! Readers for satellite gravity-field files, producing [`geopotential`] coefficient sets.
//!
//! - [`read_gfc`]: ICGEM `.gfc` models, static or time-variable (icgem1.0 and icgem2.0), with
//!   their formal or calibrated errors. Tested on real COST-G, Swarm and GOCO06s files.
//! - [`read_tn14`], [`read_tn13`]: the GRACE technical-note series that replace `C̄20`/`C̄30`
//!   (satellite laser ranging) and supply degree 1 (geocentre).
//! - [`julian_epoch`]: the time scale all of these use.
//!
//! ```no_run
//! let text = std::fs::read_to_string("GOCO06s.gfc").unwrap();
//! let gfc = gravity_field::read_gfc(&text).unwrap();
//! let july_2026 = gfc.model.at(gravity_field::julian_epoch(2026, 7, 16, 12, 0));
//! println!("C20 = {:e}", july_2026.c(2, 0));
//! ```

mod gfc;
mod technotes;

pub use gfc::{Gfc, GfcError, julian_epoch, read_gfc};
pub use technotes::{Degree1, SlrC20C30, TechNoteError, read_tn13, read_tn14};
