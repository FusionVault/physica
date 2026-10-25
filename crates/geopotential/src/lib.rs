#![doc = include_str!("../README.md")]

pub mod bodies;
mod coefficients;
pub mod fft;
mod field;
mod functional;
mod gfc;
mod legendre;
pub mod normal;
pub mod planar;
pub mod slepian;
pub mod linalg;
pub mod special;
mod synthesis;

pub use coefficients::{Coefficients, Term, TermKind, TideSystem, TimeVariable};
pub use field::{Field, geocentric};
pub use functional::{Functional, RHO_EARTH, RHO_WATER};
pub use gfc::{Gfc, GfcError, julian_epoch, parse_date, read_gfc};
pub use legendre::{Legendre, MAX_DEGREE, gauss_legendre};
pub use synthesis::{Partials, SecondPartials, gauss_grid};
