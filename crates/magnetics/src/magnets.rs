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

/// Radial field (T) at `r` (m, from the centre of a spherical conductor) of a current dipole `q`
/// (A m) at `r0`: `μ₀/(4π) (Q × (r − r₀))·r̂/|r − r₀|³`. The volume currents add nothing to the
/// radial component in a spherically symmetric conductor (Hämäläinen et al.,
/// [Rev. Mod. Phys. 65, 413 (1993)](https://doi.org/10.1103/RevModPhys.65.413), eqs. 30–32), so
/// this Biot–Savart term is exact there. A typical 10 nAm source 30 mm deep in a 100 mm head
/// peaks at 355 fT, inside the published "typically 50–500 fT".
pub fn meg_radial_field(q: [f64; 3], r0: [f64; 3], r: [f64; 3]) -> f64 {
    use crate::toroidal::{cross, dot};
    let d = [r[0] - r0[0], r[1] - r0[1], r[2] - r0[2]];
    MU_0 / (4.0 * core::f64::consts::PI) * dot(cross(q, d), r) / (dot(r, r).sqrt() * dot(d, d).powf(1.5))
}

/// Joule heating per metre (W/m) of a wire of diameter `d` (m) and resistivity `rho` (Ω m)
/// delivering power `p` (W) at voltage `v` (V): `(P/V)²ρ/(πd²/4)`. 8 µm copper is 343 Ω/m, so a
/// kilowatt at 230 V would dissipate 6.5 kW per metre and at 10 kV only 3.4 W per metre: the
/// "kilowatts over an 8 µm wire" (Moray King 2012) is ordinary at high voltage, and the claim as
/// told (voltage and source not found) is **unverified**.
pub fn wire_joule_per_length(p: f64, v: f64, d: f64, rho: f64) -> f64 {
    let i = p / v;
    i * i * rho / (core::f64::consts::PI * d * d / 4.0)
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

/// A molar susceptibility (m³/mol) from its CGS value `chi_cgs` (cm³/mol, as the CRC Handbook
/// lists them): `4π×10⁻⁶ χ`. Oxygen's +3449×10⁻⁶ cm³/mol at 293 K is 4.334×10⁻⁸ m³/mol, the most
/// of any elemental gas; the heavy rare earths are 8–51 times more per mole.
pub fn molar_susceptibility_si(chi_cgs: f64) -> f64 {
    4.0 * core::f64::consts::PI * 1e-6 * chi_cgs
}

/// A paramagnet's susceptibility at temperature `t` (K) from its value `chi` at `t0` (Curie's
/// law, `χ ∝ 1/T`).
pub fn curie_scale(chi: f64, t0: f64, t: f64) -> f64 {
    chi * t0 / t
}

/// The field-gradient product `B·dB/dz` (T²/m) that levitates a diamagnet of density `rho`
/// (kg/m³) and volume susceptibility `chi` (negative) against gravity `g`: `μ₀ρg/|χ|`.
/// Earnshaw's theorem does not forbid it because the moment is induced
/// ([Berry & Geim 1997](https://obelix.physik.uni-bielefeld.de/~schnack/Herbstakademie/CD-2008/frog-ejp.pdf)).
pub fn levitation_product(rho: f64, chi: f64, g: f64) -> f64 {
    MU_0 * rho * g / chi.abs()
}

/// The force (N) on a body of volume `volume` (m³) and volume susceptibility `chi` in a field
/// gradient: `F = (χV/μ₀) B ∂B/∂x`, for the gradient product `b_dbdx` (T²/m).
///
/// It is the same physics as [`levitation_product`], and it is the leading magnetic systematic in
/// precision gravity: a weakly magnetised source mass or nearby steel acts on a susceptible test
/// mass in step with the signal. Earth's own field gradient is far too weak to matter, but a lab
/// gradient of order 10⁻⁹ T²/m is enough to bias a Cavendish measurement at the 10⁻⁵ level, the
/// size of the spread between the measurements of G.
pub fn susceptibility_force(volume: f64, chi: f64, b_dbdx: f64) -> f64 {
    chi * volume / MU_0 * b_dbdx
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
    fn meg_and_thin_wires() {
        // A 10 nAm tangential dipole 30 mm below the surface of a 100 mm sphere. In the plane
        // normal to Q, B_r = μ₀Qd₀ sin θ/(4πD³), which peaks where Rd₀c² + (R² + d₀²)c − 3Rd₀ = 0,
        // c = cos θ: 354.9596 fT. The script's 400-point scan undershoots the maximum by 1.6 aT
        // and reports 354.958 fT, "355 fT".
        let (big_r, d0) = (0.100_f64, 0.070);
        let s = big_r * big_r + d0 * d0;
        let c = ((s * s + 12.0 * big_r * big_r * d0 * d0).sqrt() - s) / (2.0 * big_r * d0);
        let at = |th: f64| meg_radial_field([10e-9, 0.0, 0.0], [0.0, 0.0, d0], [0.0, big_r * th.sin(), big_r * th.cos()]);
        let peak = at(c.acos());
        let ft = peak * 1e15;
        assert!((ft - 355.0).abs() < 0.05 && (0.0..2e-3).contains(&(ft - 354.958)), "{ft}");
        assert!(peak > at(c.acos() + 1e-3) && peak > at(c.acos() - 1e-3) && (50e-15..500e-15).contains(&peak));
        // Directly above the dipole and along Q the radial field is zero.
        assert!(at(0.0).abs() < 1e-30 && meg_radial_field([10e-9, 0.0, 0.0], [0.0, 0.0, d0], [big_r, 0.0, 0.0]).abs() < 1e-30);
        // IACS copper, 1.7241×10⁻⁸ Ω m: 8 µm is 343.0 Ω/m (1 A), 6.48 kW/m for 1 kW at 230 V.
        let rho = 1.7241e-8;
        assert!(close(wire_joule_per_length(1.0, 1.0, 8e-6, rho), 343.0, 1e-3));
        assert!(close(wire_joule_per_length(1000.0, 230.0, 8e-6, rho), 6.48e3, 1e-3));
        assert!(close(wire_joule_per_length(1000.0, 10e3, 8e-6, rho), 3.43, 1e-3));
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
    fn molar_susceptibility() {
        // O₂: +3449×10⁻⁶ cm³/mol (CRC, 293 K) is 4.334×10⁻⁸ m³/mol, 1% above Moldover et al.'s
        // (2008) 42.92×10⁻⁹; at 273.15 K Curie's law gives 3700×10⁻⁶.
        let si = molar_susceptibility_si(3449e-6);
        assert!(close(si, 4.334e-8, 1e-3) && close(si, 42.92e-9, 0.011));
        assert!(close(curie_scale(3449e-6, 293.15, 273.15), 3700e-6, 1e-3));
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

    #[test]
    fn magnetic_systematics_in_precision_gravity() {
        // A 0.2 kg gold sphere (chi = -3.4e-5) against the Cavendish pull of 10 kg at 10 cm.
        let (volume, chi) = (0.2 / 19_300.0, -3.4e-5);
        let signal = metrology::G * 10.0 * 0.2 / 0.1f64.powi(2);
        assert!(close(signal, 1.335e-8, 1e-3));
        // Earth's own field gradient (50 uT over an Earth radius) is 2.5e-11 of the signal.
        let earth = 3.0 * 5e-5 * 5e-5 / 6.371e6;
        assert!(close(susceptibility_force(volume, chi, earth).abs() / signal, 2.47e-11, 1e-2));
        // A lab gradient of 1.05e-9 T^2/m, for instance 100 uT falling off over 10 m, biases it
        // at 2.2e-5: the size of the spread between the published measurements of G.
        let lab = 1.047e-9;
        assert!(close(susceptibility_force(volume, chi, lab).abs() / signal, 2.2e-5, 1e-2));
        // It is linear in the gradient and reverses with the sign of the susceptibility.
        assert!(close(susceptibility_force(volume, chi, 2.0 * lab) / susceptibility_force(volume, chi, lab), 2.0, 1e-12));
        assert!(susceptibility_force(volume, -chi, lab) == -susceptibility_force(volume, chi, lab));
    }
}
