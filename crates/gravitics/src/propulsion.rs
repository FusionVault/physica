//! "Propellantless" thrust: asymmetric capacitors ("lifters"), the EmDrive and Mach-effect
//! thrusters, with ion wind, the photon rocket and the balance effects that a measurement must
//! separate from them.
//!
//! - **Lifters** fly on ion wind, `F = I·d/μ`, and stop in vacuum. The Army Research Laboratory
//!   ([Bahder & Fazi 2003](https://arxiv.org/abs/physics/0211001)) and NASA (Canning et al. 2004)
//!   used the same law.
//! - **The EmDrive.** Shawyer's model gives `F = 2PQD_f/c`; Eagleworks reported 1.2 mN/kW, 360
//!   times the photon rocket's `1/c`. Dresden's balances saw nothing above 10–30 nN and showed how
//!   heating produces signals of the reported size
//!   ([Tajmar et al. 2022](https://doi.org/10.1007/s12567-021-00385-1)). A 1 K warming of a 15 cm
//!   cavity in air gives 150 µN of buoyancy alone.
//! - **Mach-effect thrusters.** Woodward's formula, in SI units, predicts about 0.1 µN for the
//!   device of [Fearn & Woodward 2013](https://arxiv.org/abs/1301.6178), which reported "a couple
//!   of microNewtons"; Dresden's balances traced such signals to vibration
//!   ([Monette et al. 2021](https://doi.org/10.1016/j.actaastro.2021.02.037)).
//!
//! ```
//! use gravitics::propulsion::photon_thrust;
//!
//! // Anything that beats light per watt needs artefact checks: Eagleworks claimed 360× it.
//! assert!((1.2e-3 / photon_thrust(1e3) - 360.0).abs() < 1.0);
//! ```

use core::f64::consts::PI;

use metrology::{C, E_CHARGE, G, G0};

/// Ion-wind thrust (N) of a current `current` (A) of ions with mobility `mobility` (m² V⁻¹ s⁻¹)
/// drifting across a gap `gap` (m): `F = I·d/μ`. In vacuum μ → ∞ and the thrust vanishes.
pub fn ion_wind_thrust(current: f64, gap: f64, mobility: f64) -> f64 {
    current * gap / mobility
}

/// Ion mobility `e/(2mν)` (m² V⁻¹ s⁻¹) of an ion of `mass` (kg) colliding `collision_rate` times a
/// second: the collision model in NASA's lifter report, which makes it the same law as
/// [`ion_wind_thrust`].
pub fn ion_mobility(mass: f64, collision_rate: f64) -> f64 {
    E_CHARGE / (2.0 * mass * collision_rate)
}

/// The most thrust (N) a current `current` (A) of ions of `mass` (kg) and `charge` (C), accelerated
/// through `voltage` (V) with nothing to collide with, can give: `I·√(2mV/q)`.
pub fn ballistic_thrust(current: f64, voltage: f64, mass: f64, charge: f64) -> f64 {
    current * (2.0 * mass * voltage / charge).sqrt()
}

/// Thrust (N) of a perfect photon rocket radiating `power` (W): `P/c`, 3.34 µN per kW.
pub fn photon_thrust(power: f64) -> f64 {
    power / C
}

/// Shawyer's EmDrive formula, `F = 2·P·Q·D_f/c` (N), for cavity power `power` (W), quality
/// factor `q` and design factor `design` (0.56 for a NASA-type cavity, up to 0.91).
pub fn shawyer_thrust(power: f64, q: f64, design: f64) -> f64 {
    2.0 * power * q * design / C
}

/// Buoyancy change (N) of a `volume` (m³) of air of density `air_density` (kg/m³) warmed by
/// `warming` (K) from `temperature` (K): `ρVgΔT/T`. A thermal artefact on any balance in air.
pub fn buoyancy_force(air_density: f64, volume: f64, warming: f64, temperature: f64) -> f64 {
    air_density * volume * G0 * warming / temperature
}

/// Woodward's predicted Mach-effect thrust in SI units (N): Eq. 10 of
/// [Fearn & Woodward 2013](https://arxiv.org/abs/1301.6178), `ω⁶m₀K_p²K_e x₀³V₀⁴/(8πGρ₀c²)`, times
/// the 4π the authors say it needs for SI units. For drive frequency `omega` (rad/s), mass `m0`
/// (kg), piezo and electrostriction constants `kp` (m/V) and `ke` (m/V²), stack length `x0` (m),
/// drive voltage `v0` (V) and density `rho0` (kg/m³).
pub fn mach_effect_thrust(omega: f64, m0: f64, kp: f64, ke: f64, x0: f64, v0: f64, rho0: f64) -> f64 {
    4.0 * PI * omega.powi(6) * m0 * kp * kp * ke * x0.powi(3) * v0.powi(4) / (8.0 * PI * G * rho0 * C * C)
}

/// General relativity's own bound on a Machian thrust from a device of mass `m` (kg) and size `l`
/// (m) driven at `f` (Hz): the Machian term is suppressed by `(GM/(Lc²))(2πfL/c)²`
/// ([Rodal, arXiv:2512.22687](https://arxiv.org/abs/2512.22687)), 8×10⁻³⁷ for a 0.1 kg, 1 cm, 50 kHz
/// drive. Set beside [`mach_effect_thrust`].
pub fn mach_effect_gr_bound(m: f64, l: f64, f: f64) -> f64 {
    G * m / (l * C * C) * (2.0 * PI * f * l / C).powi(2)
}

#[cfg(test)]
mod tests {
    use super::*;
    use metrology::AMU;

    fn close(a: f64, b: f64, rel: f64) -> bool {
        (a / b - 1.0).abs() < rel
    }

    #[test]
    fn lifters_are_ion_wind() {
        // ARL: 1 mA of N₂ ions (2.5 cm²/(V s)) across 4 cm gives 0.16 N.
        assert!(close(ion_wind_thrust(1e-3, 0.04, 2.5e-4), 0.16, 1e-9));
        // NASA's collision model gives a mobility of 1.7×10⁻⁴ for N₂ at 10¹⁰ collisions/s.
        let n2 = 28.0134 * AMU;
        assert!(close(ion_mobility(n2, 1e10), 1.72e-4, 5e-3));
        // Without air, the same ions at 30 kV give 1000× less.
        assert!(ion_wind_thrust(1e-3, 0.04, 2.5e-4) / ballistic_thrust(1e-3, 3e4, n2, E_CHARGE) > 1e3);
    }

    #[test]
    fn emdrive_and_its_artefacts() {
        assert!(close(photon_thrust(1e3), 3.336e-6, 1e-3));
        // Shawyer's formula: ~30 µN at Dresden's 11 W and loaded Q = 741 (they measured ≤ 30 nN);
        // ~1.6 mN at White's 60 W and Q = 7123 (who reported 30–130 µN).
        assert!(close(shawyer_thrust(11.0, 741.0, 0.56), 3.04e-5, 5e-3));
        assert!(close(shawyer_thrust(60.0, 7123.0, 0.56), 1.60e-3, 5e-3));
        // 1 K of warming of a 15 cm cube in air: ~140 µN (the paper: 150 µN).
        assert!(close(buoyancy_force(1.2, 0.15f64.powi(3), 1.0, 293.0), 1.36e-4, 5e-3));
    }

    #[test]
    fn woodward_predicts_a_tenth_of_a_micronewton() {
        // 32.35 kHz, the 46 g stack (19 mm) at 100 V: Eq. 10's printed prefactor 1.87×10²¹ and
        // "~8 nN"; with the SI 4π, about 0.1 µN, against "a couple of microNewtons" observed.
        let omega = 2.0 * PI * 32_350.0;
        let printed = mach_effect_thrust(omega, 0.046, 1.0, 1.0, 0.019, 100.0, 7.9e3) / (4.0 * PI);
        assert!(close(printed, 1.87e21, 5e-3));
        let si = mach_effect_thrust(omega, 0.046, 3.2e-10, 4e-11, 0.019, 100.0, 7.9e3);
        assert!(close(si / (4.0 * PI), 7.66e-9, 5e-3) && close(si, 9.62e-8, 5e-3));
    }

    #[test]
    fn machian_term_in_general_relativity() {
        assert!(close(mach_effect_gr_bound(0.1, 0.01, 5e4), 8.15498e-37, 1e-5));
    }
}
