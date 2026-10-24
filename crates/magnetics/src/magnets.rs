//! Building and using magnets: the field of loops, solenoids and Helmholtz pairs, the forces and
//! stored energy that limit them, levitation, induction, and the magnetic resonance and sensing
//! numbers of medicine.
//!
//! Two numbers organise magnet engineering: the magnetic pressure `B²/2μ₀` (398 kPa at 1 T,
//! 806 MPa at 45 T) against the strength of the structure, and the virial theorem, which caps the
//! energy a magnet can store per kilogram of structure at `σ/ρ`.
//!
//! ```
//! use magnetics::magnets::{larmor_frequency, lorentz_line_load};
//!
//! // The CFS/MIT 20 T model coil: 40.5 kA in 20.1 T loads the conductor with 814 kN/m
//! // (Hartwig et al. 2024 publish 815 kN/m).
//! assert!((lorentz_line_load(40.5e3, 20.1) / 1e3 - 814.05).abs() < 0.01);
//! // A 3 T MRI scanner's proton frequency: 127.7 MHz.
//! assert!((larmor_frequency(3.0) / 1e6 - 127.732).abs() < 1e-3);
//! ```

use core::f64::consts::TAU;

use metrology::{GAMMA_P, HBAR, MU_0};

/// On-axis field (T) of a single loop of radius `r` (m) carrying `i` (A), at distance `z` from
/// its plane: `μ₀IR²/(2(R² + z²)^(3/2))`.
pub fn loop_bz(i: f64, r: f64, z: f64) -> f64 {
    MU_0 * i * r * r / (2.0 * (r * r + z * z).powf(1.5))
}

/// On-axis field (T) of a finite solenoid of length `l` and radius `r` (m), `n` turns per metre
/// carrying `i` (A), at distance `z` from its centre:
/// `(μ₀nI/2)[(z + L/2)/√((z + L/2)² + R²) − (z − L/2)/√((z − L/2)² + R²)]`.
pub fn solenoid_bz(z: f64, r: f64, l: f64, n: f64, i: f64) -> f64 {
    let end = |s: f64| s / (s * s + r * r).sqrt();
    MU_0 * n * i / 2.0 * (end(z + l / 2.0) - end(z - l / 2.0))
}

/// Field (T) at the centre of a Helmholtz pair of coils of radius `r` (m), `n` turns each,
/// carrying `i` (A): `(4/5)^(3/2) μ₀NI/R`, the most uniform field two loops can make.
pub fn helmholtz_b0(r: f64, n: f64, i: f64) -> f64 {
    0.8f64.powf(1.5) * MU_0 * n * i / r
}

/// The Lorentz force per metre of conductor (N/m) carrying `i` (A) across a field `b` (T): `IB`.
pub fn lorentz_line_load(i: f64, b: f64) -> f64 {
    i * b
}

/// Magnetic energy (J) stored in a volume `volume` (m³) of uniform field `b` (T): `B²V/(2μ₀)`.
pub fn stored_energy(b: f64, volume: f64) -> f64 {
    crate::mhd::magnetic_pressure(b) * volume
}

/// The virial limit on stored energy per kilogram of structure (J/kg) for a structural material
/// of working stress `sigma` (Pa) and density `rho` (kg/m³): `E/M ≤ σ/ρ`. Steel at 100 MPa allows
/// 12.5 kJ/kg, which is why superconducting energy storage never competed with batteries.
pub fn virial_specific_energy(sigma: f64, rho: f64) -> f64 {
    sigma / rho
}

/// Skin depth (m) of a conductor of resistivity `rho` (Ω m) and relative permeability `mu_r` at
/// frequency `f` (Hz): `√(2ρ/(ωμ₀μ_r))`. The same formula sets induction heating depths, the
/// depth ground currents reach, and how fast a planet's core lets a field change.
pub fn skin_depth(rho: f64, f: f64, mu_r: f64) -> f64 {
    (2.0 * rho / (TAU * f * MU_0 * mu_r)).sqrt()
}

/// The field-gradient product `B·dB/dz` (T²/m) that levitates a diamagnet of density `rho`
/// (kg/m³) and volume susceptibility `chi` (negative) against gravity `g`: `μ₀ρg/|χ|`.
/// Earnshaw's theorem does not forbid it because the moment is induced
/// ([Berry & Geim 1997](https://obelix.physik.uni-bielefeld.de/~schnack/Herbstakademie/CD-2008/frog-ejp.pdf)).
pub fn levitation_product(rho: f64, chi: f64, g: f64) -> f64 {
    MU_0 * rho * g / chi.abs()
}

/// Lift-to-drag ratio of an Inductrack (passive induction) maglev at speed `v`, for a track whose
/// transition speed is `v_t` (same units): `v/v_t`, well above the transition.
pub fn inductrack_lift_drag(v: f64, v_t: f64) -> f64 {
    v / v_t
}

/// Proton Larmor (NMR and MRI) frequency in field `b` (T), Hz: `γ_p B/2π`.
pub fn larmor_frequency(b: f64) -> f64 {
    GAMMA_P * b / TAU
}

/// The partial-body specific absorption rate limit (W/kg) in IEC 60601-2-33's normal operating
/// mode, for the ratio `r` of exposed to total patient mass: `10 − 8R`.
pub fn sar_partial_body(r: f64) -> f64 {
    10.0 - 8.0 * r
}

/// Energy resolution of a magnetometer, in units of ħ: `δB²V/(2μ₀ħ)` for field noise `db`
/// (T/√Hz) over sensing volume `volume` (m³). Atomic magnetometers approach the ħ scale.
pub fn energy_resolution(db: f64, volume: f64) -> f64 {
    db * db * volume / (2.0 * MU_0 * HBAR)
}

/// The final field (T) of explosive or electromagnetic flux compression of a field `b0` from
/// radius `r0` to `r`, with flux-conservation efficiency `lambda` (1 = ideal): `λB₀(r₀/r)²`.
pub fn flux_compression(b0: f64, r0: f64, r: f64, lambda: f64) -> f64 {
    lambda * b0 * (r0 / r).powi(2)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, rel: f64) -> bool {
        (a / b - 1.0).abs() < rel
    }

    #[test]
    fn coil_fields() {
        // A long solenoid is μ₀nI inside and half that at its end.
        let (n, i) = (1e4, 2.0);
        assert!(close(solenoid_bz(0.0, 0.01, 100.0, n, i), MU_0 * n * i, 1e-7));
        assert!(close(solenoid_bz(50.0, 0.01, 100.0, n, i), MU_0 * n * i / 2.0, 1e-6));
        // Helmholtz: coefficient 0.715542, equal to two loops at ±R/2, and flat at the centre.
        assert!(close(helmholtz_b0(1.0, 1.0, 1.0) / MU_0, 0.715_542, 1e-6));
        assert!(close(helmholtz_b0(0.3, 1.0, 5.0), 2.0 * loop_bz(5.0, 0.3, 0.15), 1e-12));
        let pair = |z: f64| loop_bz(1.0, 1.0, z - 0.5) + loop_bz(1.0, 1.0, z + 0.5);
        let h = 1e-3;
        assert!(((pair(h) - 2.0 * pair(0.0) + pair(-h)) / (h * h)).abs() < 1e-5 * pair(0.0));
    }

    #[test]
    fn forces_and_energy() {
        assert!(close(lorentz_line_load(40.5e3, 20.1), 814_050.0, 1e-12)); // published 815 kN/m
        assert!(close(stored_energy(16.5, 1.0), 1.083e8, 1e-3));
        assert!(close(virial_specific_energy(100e6, 8000.0), 12_500.0, 1e-12)); // Tixador: about 10 kJ/kg
        assert!(close(flux_compression(3.0, 25.0, 1.0, 1.0), 1875.0, 1e-12));
    }

    #[test]
    fn induction_and_levitation() {
        // Copper (σ = 5.96×10⁷ S/m): 8.42 mm at 60 Hz, 9.22 mm at 50 Hz, 65.2 µm at 1 MHz.
        let rho = 1.0 / 5.96e7;
        assert!(close(skin_depth(rho, 60.0, 1.0), 8.42e-3, 1e-3) && close(skin_depth(rho, 50.0, 1.0), 9.22e-3, 1e-3));
        assert!(close(skin_depth(rho, 1e6, 1.0), 65.2e-6, 1e-3));
        // Water floats where B dB/dz = 1400.9 T²/m (χ = −8.8×10⁻⁶).
        assert!(close(levitation_product(1000.0, -8.8e-6, 9.81), 1400.9, 1e-4));
        // Both of Post's Inductrack examples imply a 1.67 km/h transition speed.
        assert!(close(500.0 / 300.0, 50.0 / 30.0, 1e-12) && close(inductrack_lift_drag(500.0, 500.0 / 300.0), 300.0, 1e-12));
    }

    #[test]
    fn medicine_and_sensing() {
        for (b, mhz) in [(1.5, 63.866), (3.0, 127.732), (7.0, 298.042), (11.7, 498.16)] {
            assert!(close(larmor_frequency(b) / 1e6, mhz, 2e-5), "{b} T");
        }
        assert!((sar_partial_body(0.5) - 6.0).abs() < 1e-12);
        // Dang et al. 2010: 160 aT/√Hz over 0.45 cm³ is 44ħ.
        assert!(close(energy_resolution(160e-18, 0.45e-6), 43.5, 3e-3));
    }
}
