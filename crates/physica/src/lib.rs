//! The quantum-gravimeter crate family in one dependency, and the equations its domains share.
//!
//! Each crate is re-exported under its own name, so `physica::magnetics::main_field::wmm2025()`
//! and `magnetics::main_field::wmm2025()` are the same function. [`analogies`] collects the
//! identities that join the domains, each one a test.

pub use geopotential;
pub use gravimetry;
pub use gravitational;
pub use gravitics;
pub use magnetics;
pub use metrology;
pub use plasmoid;

pub mod analogies;
