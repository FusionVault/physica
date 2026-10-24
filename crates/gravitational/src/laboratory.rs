//! Laboratory gravity, sized honestly: the gravitational field of a light or particle beam, the
//! gravitational waves a spinning body radiates, and the energy a wave carries.
//!
//! A beam of light or matter really does carry gravity. Parallel beams do not attract each other
//! (Tolman, Ehrenfest & Podolsky 1931), but a test mass at rest near a beam of power P is pulled
//! by `a = 4GP/(c³ρ)` while the beam passes ([Rätzel, Wilkens & Menzel 2016](https://arxiv.org/abs/1511.01023)).
//! This is a near field confined to the beam's surroundings, not a projected ray. Radiated waves
//! from lab rotors are suppressed by G/c⁵ (the quadrupole formula).
//!
//! ```
//! use gravitational::laboratory as lab;
//!
//! // A petawatt laser, felt 1 mm from the beam:
//! let a = lab::beam_acceleration(1e15, 1e-3);
//! assert!((a / 9.91e-18 - 1.0).abs() < 0.01); // 10⁻¹⁷ m/s²: no instrument can see it
//! ```

use metrology::{C, G};
use core::f64::consts::PI;

/// Acceleration of a test mass at rest at distance `rho` (m) from a beam of power `power` (W),
/// while the beam passes: `a = 4GP/(c³ρ)` (twice the Newtonian line-mass value).
pub fn beam_acceleration(power: f64, rho: f64) -> f64 {
    4.0 * G * power / (C.powi(3) * rho)
}

/// Tidal acceleration gradient next to the beam, s⁻²: `−4GP/(c³ρ²)`. This is
/// [Spengler, Rätzel & Braun 2022](https://arxiv.org/abs/2104.09209)'s Eq. 1, `R₀ₓ₀ₓ ≃ −4GP/(c⁵ρ²)`
/// in m⁻² (with x⁰ = ct), times c².
pub fn beam_tidal(power: f64, rho: f64) -> f64 {
    -4.0 * G * power / (C.powi(3) * rho * rho)
}

/// Velocity kick from a pulse of energy `energy` (J) passing at distance `rho` (m): the
/// acceleration integrated over the pulse duration, `Δv = 4GE/(c³ρ)`.
pub fn beam_velocity_kick(energy: f64, rho: f64) -> f64 {
    4.0 * G * energy / (C.powi(3) * rho)
}

/// Gravitational-wave power from a rod of mass `m` (kg) and length `l` (m) spinning about its
/// middle at `omega` (rad/s): `P = (2/45) G m² l⁴ ω⁶/c⁵`, W.
pub fn spinning_rod_power(m: f64, l: f64, omega: f64) -> f64 {
    2.0 / 45.0 * G * m * m * l.powi(4) * omega.powi(6) / C.powi(5)
}

/// Gravitational-wave power from a rigid body spinning at `omega` about a principal axis with
/// the other two moments of inertia differing by `delta_i` (kg m²):
/// `P = (32/5)(G/c⁵)(ΔI)² ω⁶`, W.
pub fn rotor_power(delta_i: f64, omega: f64) -> f64 {
    32.0 / 5.0 * G / C.powi(5) * delta_i * delta_i * omega.powi(6)
}

/// Wavelength of the waves from a body spinning at `omega`: they come out at twice the spin
/// frequency, so `λ = πc/ω`, m.
pub fn rotor_wavelength(omega: f64) -> f64 {
    PI * C / omega
}

/// Amplitude (N) of the harmonic at `n` times the rotation frequency of the force along `x` on a
/// test mass `m_test` (kg) at the origin from a rotor of `n` equal point masses `m_point` (kg)
/// spaced evenly on a circle of `radius` (m) about an axis parallel to z through
/// `(d cos φ, d sin φ, z)`: a direct sum of Newton's law over one revolution, the brute-force check
/// of the LIGO Newtonian calibrator ([Ross et al. 2021](https://arxiv.org/abs/2107.00141)).
pub fn rotor_force(m_test: f64, m_point: f64, n: usize, radius: f64, d: f64, phi: f64, z: f64) -> f64 {
    const STEPS: usize = 720;
    let (cx, cy) = (d * phi.cos(), d * phi.sin());
    let (mut re, mut im) = (0.0, 0.0);
    for i in 0..STEPS {
        let theta = 2.0 * PI * i as f64 / STEPS as f64;
        let fx: f64 = (0..n).map(|j| {
            let a = theta + 2.0 * PI * j as f64 / n as f64;
            let (x, y) = (cx + radius * a.cos(), cy + radius * a.sin());
            G * m_test * m_point * x / (x * x + y * y + z * z).powf(1.5)
        }).sum();
        let (s, c) = (n as f64 * theta).sin_cos();
        re += fx * c;
        im -= fx * s;
    }
    2.0 * re.hypot(im) / STEPS as f64
}

/// Displacement amplitude (m) of a free mass `mass` (kg) driven by a force of amplitude `force` (N)
/// at frequency `f`, well above its suspension resonance: `F/(M(2πf)²)`.
pub fn free_mass_displacement(force: f64, mass: f64, f: f64) -> f64 {
    force / (mass * (2.0 * PI * f).powi(2))
}

/// The strength `R = Gbρ sin α (r_max⁴ − r_min⁴)/(32π²)` (m⁵ s⁻²) of a two-sector rotor (sectors of
/// thickness `b`, density `rho`, opening angle `alpha`, radii `r_min`…`r_max`), as used by Virgo's
/// Newtonian calibrator ([Estevez et al. 2018](https://arxiv.org/abs/1806.06572)). The mirror
/// moves by `x f² = R·P` with [`sector_rotor_geometry`]'s P.
pub fn sector_rotor_strength(b: f64, rho: f64, alpha: f64, r_max: f64, r_min: f64) -> f64 {
    G * b * rho * alpha.sin() * (r_max.powi(4) - r_min.powi(4)) / (32.0 * PI * PI)
}

/// The geometric factor `P = [(15/Z − 6)² cos²φ + 36 sin²φ]^(1/2) Z^(−5/2) d^(−4)` (m⁻⁴), with
/// `Z = 1 + (z/d)²`, for a rotor at distance `d`, angle `phi` and height `z` from the mirror.
pub fn sector_rotor_geometry(d: f64, phi: f64, z: f64) -> f64 {
    let zz = 1.0 + (z / d).powi(2);
    ((15.0 / zz - 6.0).powi(2) * phi.cos().powi(2) + 36.0 * phi.sin().powi(2)).sqrt() * zz.powf(-2.5) / d.powi(4)
}

/// Efficiency `64GMv_s⁴/(15L²c⁵Γ)` of a bar of mass `m`, sound speed `vs` and length `l` as a
/// generator of gravitational waves, against its heating rate `gamma_heat` (s⁻¹, ω/Q): of order
/// 10⁻³³, so a lab generator is hopeless.
pub fn weinberg_efficiency(m: f64, vs: f64, l: f64, gamma_heat: f64) -> f64 {
    64.0 * G * m * vs.powi(4) / (15.0 * l * l * C.powi(5) * gamma_heat)
}

/// Energy flux of a plane wave of amplitude `h` at frequency `f`: `F = c³ω²h²/(32πG)`, W/m².
pub fn wave_flux(h: f64, f: f64) -> f64 {
    let w = 2.0 * PI * f;
    C.powi(3) * w * w * h * h / (32.0 * PI * G)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn beam_fixtures() {
        assert!((beam_acceleration(1e15, 1e-3) / 9.91e-18 - 1.0).abs() < 2e-3);
        // The LHC beam, 3.8×10¹² W, at 100 µm.
        assert!((beam_acceleration(3.8e12, 1e-4) / 3.77e-19 - 1.0).abs() < 2e-3);
        // A 100 J pulse at 1 mm: about 10⁻³⁰ m/s.
        assert!((beam_velocity_kick(100.0, 1e-3) / 9.9e-31 - 1.0).abs() < 0.01);
        // 1 PW at 1 mm: −1.10×10⁻³¹ m⁻², i.e. −9.91×10⁻¹⁵ s⁻²; the gradient of the acceleration.
        assert!((beam_tidal(1e15, 1e-3) / -9.91e-15 - 1.0).abs() < 2e-3);
        assert!((beam_tidal(1e15, 1e-3) * 1e-3 / -beam_acceleration(1e15, 1e-3) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn misner_thorne_wheeler_rod() {
        // 490 t, 20 m, spun to break-up at 28 rad/s: 2.27×10⁻²⁹ W at λ ≈ 3.36×10⁷ m.
        let p = spinning_rod_power(4.9e5, 20.0, 28.0);
        assert!((p / 2.27e-29 - 1.0).abs() < 0.01, "{p:e}");
        assert!((rotor_wavelength(28.0) / 3.36e7 - 1.0).abs() < 0.01);
        // The rod is a rotor with ΔI = m l²/12.
        let i = 4.9e5 * 400.0 / 12.0;
        assert!((rotor_power(i, 28.0) / p - 1.0).abs() < 1e-12);
    }

    #[test]
    fn flux_scales_as_h_squared() {
        assert!((wave_flux(1e-21, 100.0) / 1.5853e-3 - 1.0).abs() < 1e-4);
        assert!((wave_flux(2e-21, 100.0) / wave_flux(1e-21, 100.0) - 4.0).abs() < 1e-12);
    }

    #[test]
    fn newtonian_calibrators() {
        // LIGO's rotor (Ross et al. Table I): 19.10 pN at 2f, 9.017 pN at 3f by brute force.
        let (mirror, m, d, phi, z) = (39.680, 1.0558, 1.1802, 52.24f64.to_radians(), 0.010);
        assert!((rotor_force(mirror, m, 2, 0.06033, d, phi, z) / 19.10e-12 - 1.0).abs() < 1e-3);
        assert!((rotor_force(mirror, m, 3, 0.10478, d, phi, z) / 9.017e-12 - 1.0).abs() < 1e-3);
        let x = free_mass_displacement(19.16e-12, 39.68, 15.60);
        assert!((x / 5.026e-17 - 1.0).abs() < 1e-3 && (x / 3995.0 / 1.258e-20 - 1.0).abs() < 1e-3);
        // Virgo's 2018 rotor: x f² = 3.3×10⁻¹⁴ m Hz².
        let rp = sector_rotor_strength(0.052, 2800.0, 45f64.to_radians(), 0.190, 0.055) * sector_rotor_geometry(1.32, 0.241, 0.643);
        assert!((rp / 3.33101e-14 - 1.0).abs() < 1e-5);
        assert!((weinberg_efficiency(1000.0, 5e3, 1.5, 2.0 * PI * 1e3 / 1e6) / 5.19887e-33 - 1.0).abs() < 1e-5);
    }

    #[test]
    fn rotor_and_flux_conventions() {
        // The MTW beam as a rotor with I₁ = M(L²/12 + r²/4), I₂ = Mr²/2: 2.234×10⁻²⁹ W.
        let (m, l, r) = (4.9e5, 20.0, 1.0);
        assert!((rotor_power(m * (l * l / 12.0 + r * r / 4.0) - m * r * r / 2.0, 28.0) / 2.23378e-29 - 1.0).abs() < 1e-5);
        assert!((wave_flux(1e-22, 1e3) / 1.58532e-3 - 1.0).abs() < 1e-5);
        // The GW150914 paper's c³ḣ²/(16πG) is twice this flux.
        let hdot = 2.0 * PI * 250.0 * 1e-21;
        assert!((C.powi(3) * hdot * hdot / (16.0 * PI * G) / wave_flux(1e-21, 250.0) - 2.0).abs() < 1e-12);
    }
}
