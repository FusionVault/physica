//! Spherical-harmonic gravity fields in pure Rust: the engine behind geoid, gravity and
//! GRACE-style water-storage maps, and behind error spectra for simulated (quantum) sensors.
//!
//! - [`Legendre`]: fully normalized associated Legendre functions and their latitude
//!   derivatives, stable to degree 2700 at every latitude (Holmes & Featherstone 2002).
//! - [`Coefficients`]: `C̄nm, S̄nm` with their `GM`, radius and [`TideSystem`]; degree
//!   variances and amplitudes; tide-system conversion.
//! - [`TimeVariable`]: ICGEM `gfct`/`trnd`/`acos`/`asin` models evaluated at any epoch.
//! - [`Coefficients::partials`] (point potential and gradient),
//!   [`Coefficients::second_partials`] (and its Hessian, for gradiometry) and [`Coefficients::grid`]
//!   (any [`Functional`] on a latitude × longitude grid).
//!
//! Conventions follow IERS 2010 ch. 6 and Barthelmes (ICGEM STR 09/02, rev. 2013): geocentric
//! latitude, 4π normalization, no Condon–Shortley phase, and
//! `V = GM/r · Σ_n (R/r)^n Σ_m P̄nm(sin φ)(C̄nm cos mλ + S̄nm sin mλ)`.
//!
//! ```
//! use geopotential::{Coefficients, Functional};
//!
//! // A J2-only field: the geoid it implies at the equator and pole (spherical approximation).
//! let mut j2 = Coefficients::zeros(3.986004415e14, 6378136.3, 2);
//! j2.set(2, 0, -4.84165e-4, 0.0);
//! let n = j2.grid(&Functional::GeoidHeight, j2.radius, &[0.0, 90f64.to_radians()], &[0.0]);
//! assert!(n[0] > 0.0 && n[1] < 0.0); // bulge at the equator, depression at the pole
//! ```

mod coefficients;
mod functional;
mod legendre;
mod synthesis;

pub use coefficients::{Coefficients, Term, TermKind, TideSystem, TimeVariable};
pub use functional::{Functional, RHO_EARTH, RHO_WATER};
pub use legendre::{Legendre, MAX_DEGREE};
pub use synthesis::{Partials, SecondPartials};
