#![doc = include_str!("../README.md")]

pub mod hfgw;
pub mod propulsion;
pub mod spacetime;
pub mod strong_field;
pub mod superconductors;

pub use gwsens::{C, G, HBAR, K_B, MU_0};

/// Elementary charge, C.
pub const E_CHARGE: f64 = 1.602_176_634e-19;
/// Electron mass, kg.
pub const M_E: f64 = 9.109_383_701_5e-31;
/// Vacuum permittivity, F/m.
pub const EPS_0: f64 = 8.854_187_812_8e-12;
/// Atomic mass unit, kg.
pub const AMU: f64 = 1.660_539_066_6e-27;
/// Standard gravity, m/s².
pub const G0: f64 = 9.806_65;
/// Planck length, m.
pub const PLANCK_LENGTH: f64 = 1.616_255e-35;
