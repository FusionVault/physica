#![doc = include_str!("../README.md")]

mod coefficients;
mod functional;
mod legendre;
mod synthesis;

pub use coefficients::{Coefficients, Term, TermKind, TideSystem, TimeVariable};
pub use functional::{Functional, RHO_EARTH, RHO_WATER};
pub use legendre::{Legendre, MAX_DEGREE};
pub use synthesis::{Partials, SecondPartials};
