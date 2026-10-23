//! Readers for satellite gravity-field files, producing [`geopotential`] coefficient sets.
//!
//! - [`read_gfc`]: ICGEM `.gfc` models, static or time-variable (icgem1.0 and icgem2.0), with
//!   their formal or calibrated errors. Tested on real COST-G, Swarm and GOCO06s files.
//! - [`julian_epoch`]: the time scale ICGEM time-variable terms are evaluated in.
//!
//! ```no_run
//! let text = std::fs::read_to_string("GOCO06s.gfc").unwrap();
//! let gfc = gravity_field::read_gfc(&text).unwrap();
//! let july_2026 = gfc.model.at(gravity_field::julian_epoch(2026, 7, 16, 12, 0));
//! println!("C20 = {:e}", july_2026.c(2, 0));
//! ```

mod gfc;

pub use gfc::{Gfc, GfcError, julian_epoch, read_gfc};
