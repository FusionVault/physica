//! Z-pinches and the dense plasma focus: Bennett equilibrium, the Pease–Braginskii current, Lee
//! model scalings, and the pinch's gravitational twin, the self-gravitating filament.
//!
//! A Z-pinch holds its plasma with its own current's field. In equilibrium the magnetic pressure
//! balances the plasma pressure: `μ₀I²/8π = N k_B(T_e + T_i)` (Bennett 1934). The same equation,
//! with `Gm²` in place of `μ₀e²u²/4π`, describes an isothermal filament held by its own gravity
//! (Ostriker 1964): its line mass is `2c_s²/G`. In a laboratory pinch the magnetic binding is
//! about 10²⁷ times the gravitational.
//!
//! ```
//! use plasmoid::pinch::bennett_temperature;
//!
//! // 1 MA through 10²⁰ ions per metre (10¹⁸ per cm), hydrogen: T_e = T_i = 1.56 keV.
//! assert!((bennett_temperature(1e6, 1e20, 1.0) / 1560.0 - 1.0).abs() < 0.005);
//! ```

use core::f64::consts::PI;

use metrology::{E_CHARGE as E, G, HBAR, K_B, M_E, M_H, M_P, MU_0 as MU0};

/// Bennett temperature (eV) with `T_e = T_i`, for current `current` (A), ion line density
/// `line_density` (m⁻¹) and ion charge `z`: `μ₀I²/(8πN(1 + z)e)`.
pub fn bennett_temperature(current: f64, line_density: f64, z: f64) -> f64 {
    MU0 * current * current / (8.0 * PI * line_density * (1.0 + z) * E)
}

/// Ion line density (m⁻¹) a current (A) confines at `T_e + T_i = t_sum` (eV), singly charged.
pub fn bennett_line_density(current: f64, t_sum: f64) -> f64 {
    MU0 * current * current / (8.0 * PI * t_sum * E)
}

/// Which Spitzer resistivity a calculation uses.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Resistivity {
    /// Across the field (NRL: `η⊥ = 1.03×10⁻⁴ Z lnΛ T^(−3/2)` Ω m, T in eV).
    Perpendicular,
    /// Along the field, `η⊥/1.96`.
    Parallel,
}

/// Pease–Braginskii current (A): above it bremsstrahlung beats Ohmic heating and a hydrogen pinch
/// collapses radiatively, `(π/μ₀)√(48 lnΛ/(σ₀A))·(1 + z)/z` (Klír 2005), about 1.4 MA.
pub fn pease_braginskii(ln_lambda: f64, z: f64, resistivity: Resistivity) -> f64 {
    let c = match resistivity {
        Resistivity::Perpendicular => 9.71e3,
        Resistivity::Parallel => 1.96 * 9.71e3,
    };
    let sigma0 = c / E.powf(1.5); // S/m per J^(3/2)
    let a = 1.69e-38 / E.sqrt(); // bremsstrahlung, W m³ J^(−1/2)
    PI / MU0 * (48.0 * ln_lambda / (sigma0 * a)).sqrt() * (1.0 + z) / z
}

/// Peak current (A) of a plasma focus on a bank of voltage `v0` (V), inductance `l0` (H) and
/// capacitance `c0` (F): `V₀/(Z₀ + 7 mΩ + 0.1 Z₀)` with `Z₀ = √(L₀/C₀)` and the run-down's dynamic
/// resistance (Lee 2009). It saturates near `V₀/7 mΩ` however large the bank.
pub fn lee_peak_current(v0: f64, l0: f64, c0: f64) -> f64 {
    let z0 = (l0 / c0).sqrt();
    v0 / (z0 + 7e-3 + 0.1 * z0)
}

/// Lee's fitted D–D neutron yield per shot from the peak current (MA), `1.8×10¹⁰ I^3.8`.
pub fn lee_yield(peak_ma: f64) -> f64 {
    1.8e10 * peak_ma.powf(3.8)
}

/// Axial run-down ("snowplow") speed (m/s) of a plasma focus with inner and outer electrodes of
/// radii `a`, `b` (m), current `current` (A) and fill density `rho` (kg/m³), with mass and
/// current fractions `fm`, `fc`.
pub fn snowplow_speed(current: f64, a: f64, b: f64, rho: f64, fm: f64, fc: f64) -> f64 {
    let c = b / a;
    fc / fm.sqrt() * (MU0 * c.ln() / (4.0 * PI * PI * (c * c - 1.0) * rho)).sqrt() * current / a
}

/// Electron cyclotron energy `ħeB/m_e` (eV) in a field `b` (T): 11.6 eV at 1 GG (10⁵ T), the
/// scale of LPPFusion's proposed magnetic suppression of bremsstrahlung (contested).
pub fn cyclotron_energy(b: f64) -> f64 {
    HBAR * E * b / M_E / E
}

/// The gravitational twin: critical line mass (kg/m) of an isothermal filament at temperature
/// `t` (K) and mean molecular weight `mu`, `2c_s²/G` (Ostriker 1964).
pub fn ostriker_line_mass(t: f64, mu: f64) -> f64 {
    2.0 * K_B * t / (mu * M_H) / G
}

/// The drift speed (m/s) at which two protons attract magnetically as strongly as they do
/// gravitationally: `m_p√(4πG/μ₀)/e`.
pub fn equal_attraction_drift() -> f64 {
    M_P * (4.0 * PI * G / MU0).sqrt() / E
}

/// How many times a pinch's magnetic self-binding exceeds its gravitational self-binding:
/// `(μ₀I²/4π)/(Gλ²)` for current `current` (A) and line mass `line_mass` (kg/m).
pub fn magnetic_over_gravitational(current: f64, line_mass: f64) -> f64 {
    MU0 * current * current / (4.0 * PI) / (G * line_mass * line_mass)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, rel: f64) -> bool {
        (a / b - 1.0).abs() < rel
    }

    #[test]
    fn bennett_and_pease_braginskii() {
        // Klír's practical form: T_B = 3.12 keV I² /((1 + z) N) with I in MA, N in 10¹⁸ cm⁻¹.
        assert!(close(2.0 * bennett_temperature(1e6, 1e20, 1.0) / 1e3, 3.12, 0.005));
        // FuZE 2019: 200 kA, 10²³ m⁻³, 3 mm radius: 2.2 keV (measured 1–2 keV).
        let n = 1e23 * PI * 3e-3 * 3e-3;
        assert!(close(bennett_temperature(200e3, n, 1.0), 2200.0, 0.02));
        assert!(close(bennett_line_density(1e6, 2000.0) * 2000.0 * E, MU0 * 1e12 / (8.0 * PI), 1e-12));
        assert!(close(pease_braginskii(10.0, 1.0, Resistivity::Perpendicular), 1.37e6, 0.01));
        assert!(close(pease_braginskii(10.0, 1.0, Resistivity::Parallel), 0.98e6, 0.01));
    }

    #[test]
    fn lee_model() {
        // Lee 2009 Table I (30 kV, 30 nH): peak currents for 0.45 kJ to 4.5 MJ banks.
        for (c0, published) in [(1e-6, 152e3), (10e-6, 447e3), (100e-6, 1156e3), (1e-3, 2321e3), (10e-3, 3407e3)] {
            assert!(close(lee_peak_current(30e3, 30e-9, c0), published, 0.015), "{c0}");
        }
        assert!(close(lee_yield(1.1), 2.6e10, 0.01));
        // Run-down at 1 MA in FF-1's fill (1.6×10²⁴ D/m³): of order 10 cm/µs.
        let v = snowplow_speed(1e6, 0.028, 0.05, 1.6e24 * 3.343_583_772e-27, 0.1, 0.7);
        assert!((5e4..3e5).contains(&v), "{v}");
    }

    #[test]
    fn the_pinch_and_its_gravitational_twin() {
        assert!(close(cyclotron_energy(1e5), 11.58, 0.005));
        // A molecular-cloud filament at 10 K: 16.6 solar masses per parsec.
        assert!(close(ostriker_line_mass(10.0, 2.3) * 3.0857e16 / 1.989e30, 16.6, 0.01));
        assert!(close(equal_attraction_drift(), 2.70e-10, 0.005));
        // A 1 MA deuterium pinch at 2 keV: magnetic binding 5.5×10²⁷ times the gravitational.
        let line_mass = bennett_line_density(1e6, 2000.0) * 3.343_583_772e-27;
        assert!(close(magnetic_over_gravitational(1e6, line_mass), 5.5e27, 0.02));
    }
}
