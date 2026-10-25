#![doc = include_str!("../README.md")]

mod coefficients;
mod field;
mod functional;
mod gfc;
mod legendre;
pub mod special;
mod synthesis;

pub use coefficients::{Coefficients, Term, TermKind, TideSystem, TimeVariable};
pub use field::{Field, geocentric};
pub use functional::{Functional, RHO_EARTH, RHO_WATER};
pub use gfc::{Gfc, GfcError, julian_epoch, parse_date, read_gfc};
pub use legendre::{Legendre, MAX_DEGREE};
pub use synthesis::{Partials, SecondPartials};
