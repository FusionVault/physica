//! Forward models for a buried body: the gravity and gravity gradient of a sphere and of a
//! horizontal cylinder, the shapes a ground survey actually fits.
//!
//! A gravimeter reads the vertical attraction (µGal); a gradiometer reads the tidal tensor
//! (Eötvös), in which a vehicle's own acceleration cancels because it is common to both proof
//! masses. The gradient falls off one power of distance faster than the attraction, so it
//! localises a shallow body much more sharply.
//!
//! Signs follow [`geopotential::SecondPartials::tensor_neu`]: `gz` is positive when the body
//! pulls down harder (a mass excess), and `gradient` is `Γ_uu`, positive above a mass excess.
//!
//! ```
//! use gravimetry::local::{EOTVOS, MICROGAL, sphere};
//!
//! // A 5 m void (Δρ = −2700 kg/m³) with its centre 20 m down, measured right above it.
//! let (gz, gradient) = sphere(5.0, 20.0, -2700.0, 0.0);
//! assert!((gz / MICROGAL + 23.6).abs() < 0.1);   // −23.6 µGal
//! assert!((gradient / EOTVOS + 23.6).abs() < 0.1); // −23.6 E
//! ```

use core::f64::consts::PI;

use metrology::G;

/// One microgal in m/s² (10⁻⁸), the unit a field gravimeter reports.
pub const MICROGAL: f64 = 1e-8;
/// One Eötvös in s⁻² (10⁻⁹), the unit a gradiometer reports.
pub const EOTVOS: f64 = 1e-9;

/// Vertical attraction (m/s²) and vertical gradient `Γ_uu` (s⁻²) of a buried sphere of radius `a`
/// (m) whose centre is at depth `depth` (m), with density contrast `density_contrast` (kg/m³),
/// measured at horizontal offset `offset` (m):
/// `g_z = GΔM z/(x²+z²)^(3/2)`, `Γ_uu = GΔM(2z² − x²)/(x²+z²)^(5/2)`, with `ΔM = (4/3)πa³Δρ`.
pub fn sphere(a: f64, depth: f64, density_contrast: f64, offset: f64) -> (f64, f64) {
    let dm = 4.0 / 3.0 * PI * a.powi(3) * density_contrast;
    let (x, z) = (offset, depth);
    let s2 = x * x + z * z;
    (G * dm * z / s2.powf(1.5), G * dm * (2.0 * z * z - x * x) / s2.powf(2.5))
}

/// The same for an infinite horizontal cylinder (a tunnel or pipe) of radius `a` (m) whose axis
/// is at depth `depth` (m): `g_z = 2Gλz/(x²+z²)`, `Γ_uu = 2Gλ(z² − x²)/(x²+z²)²`, with line
/// density `λ = πa²Δρ`. In two dimensions `Γ_nn = −Γ_uu` exactly.
pub fn cylinder(a: f64, depth: f64, density_contrast: f64, offset: f64) -> (f64, f64) {
    let lambda = PI * a * a * density_contrast;
    let (x, z) = (offset, depth);
    let s2 = x * x + z * z;
    (2.0 * G * lambda * z / s2, 2.0 * G * lambda * (z * z - x * x) / (s2 * s2))
}

/// The averaging time (s) a white-noise sensor needs to reach precision `target` from an
/// amplitude spectral density `asd` (units per √Hz): `(asd/target)²`.
pub fn averaging_time(asd: f64, target: f64) -> f64 {
    (asd / target).powi(2)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn buried_bodies() {
        // An air-filled tunnel, 1 m radius with its axis 1.9 m down (the Birmingham survey's
        // geometry): −44 µGal and −232 E right above it.
        let (gz, gradient) = cylinder(1.0, 1.9, -2000.0, 0.0);
        assert!((gz / MICROGAL + 44.14).abs() < 0.01 && (gradient / EOTVOS + 232.3).abs() < 0.1);
        // The gradient localises it: at one depth's offset it has already changed sign.
        assert!(cylinder(1.0, 1.9, -2000.0, 2.5).1 > 0.0);
        // In 2-D the trace vanishes, so Γ_nn = −Γ_uu.
        let void = sphere(5.0, 20.0, -2700.0, 0.0);
        assert!((void.0 / MICROGAL + 23.59).abs() < 0.01 && (void.1 / EOTVOS + 23.59).abs() < 0.01);
        // Both scale with the contrast and vanish without one.
        assert!(sphere(5.0, 20.0, 0.0, 3.0) == (0.0, 0.0));
        assert!((sphere(5.0, 20.0, -5400.0, 3.0).0 / sphere(5.0, 20.0, -2700.0, 3.0).0 - 2.0).abs() < 1e-12);
        // Far from a sphere the gradient is the point-mass 2GΔM/z³.
        let dm = 4.0 / 3.0 * PI * 125.0 * 2700.0;
        assert!((sphere(5.0, 100.0, 2700.0, 0.0).1 / (2.0 * G * dm / 1e6) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn sensor_averaging() {
        // The Exail AQG: 500 nm/s²/√Hz reaches 20 nm/s² in 10 minutes, and 2 µGal in 625 s.
        assert!((averaging_time(500e-9, 20e-9) - 625.0).abs() < 1e-9);
        assert!((500e-9 / 600f64.sqrt() / 20.4e-9 - 1.0).abs() < 0.01);
    }
}
