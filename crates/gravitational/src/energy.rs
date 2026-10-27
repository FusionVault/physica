//! The energy gravitational waves carry: the Planck luminosity `c⁵/G`, the strain a given
//! luminosity produces, the quadrupole flux of a circular binary and its multipole coefficients,
//! and the orbital decay of binary pulsars (Peters & Mathews 1963).
//!
//! `c⁵/G = 3.6×10⁵² W` is the natural unit of gravitational-wave power, not a bound: the
//! brightest events reach about a thousandth of it (GW150914 at 9.6×10⁻⁴, GW190521 at
//! 1.0×10⁻³). The pulsar formula is the oldest confirmed prediction here; B1913+16's decay
//! matches it to 0.2% ([Weisberg & Huang 2016](https://arxiv.org/abs/1606.02744)) and the
//! double pulsar's to 0.01% ([Kramer et al. 2021](https://arxiv.org/abs/2112.06795)).
//!
//! ```
//! use gravitational::energy::{PLANCK_LUMINOSITY, peters_pbdot};
//! use metrology::SOLAR_MASS;
//!
//! // GW150914 peaked at 3.5×10⁵⁶ erg/s: under a thousandth of c⁵/G.
//! assert!((3.5e49 / PLANCK_LUMINOSITY / 9.6465e-4 - 1.0).abs() < 1e-4);
//! // The Hulse–Taylor pulsar loses 2.40×10⁻¹² s of orbital period per second.
//! let pbdot = peters_pbdot(1.438 * SOLAR_MASS, 1.390 * SOLAR_MASS, 0.322997448918 * 86_400.0, 0.6171340);
//! assert!((pbdot / -2.40263e-12 - 1.0).abs() < 1e-3);
//! ```

use core::f64::consts::{PI, SQRT_2};

use metrology::{C, G};

/// The Planck luminosity `c⁵/G`, W. A gravitational-wave source radiating `ε c⁵/G` has
/// luminosity fraction `ε`.
pub const PLANCK_LUMINOSITY: f64 = C * C * C * C * C / G;

/// Amplitude `h₀` (sinusoidal, optimally oriented) at distance `r` (m) of a source radiating the
/// fraction `fraction` of [`PLANCK_LUMINOSITY`] at frequency `f` (Hz), from the flux
/// `c³ω²h₀²/(32πG)`: `h₀ = 2√2 √ε c/(2πf r)`. At one wavelength from a source at `c⁵/G` this is
/// 0.45: the strain cannot be of order one in the wave zone.
pub fn strain_from_luminosity(fraction: f64, f: f64, r: f64) -> f64 {
    2.0 * SQRT_2 * fraction.sqrt() * C / (2.0 * PI * f * r)
}

/// Quadrupole luminosity of a circular binary in units of [`PLANCK_LUMINOSITY`]:
/// `(32/5) η² x⁵`, with `eta` the symmetric mass ratio and `x = (GMω/c³)^(2/3)` the
/// post-Newtonian parameter (1/6 at the innermost stable orbit).
pub fn circular_flux(eta: f64, x: f64) -> f64 {
    32.0 / 5.0 * eta * eta * x.powi(5)
}

/// The energy-flux coefficients `(a_ℓ, b_ℓ)` of the mass and current multipoles of order `l ≥ 2`
/// ([Blanchet's Living Review](https://arxiv.org/abs/1310.1528), Eq. 68):
/// `a_ℓ = (ℓ+1)(ℓ+2)/((ℓ−1)ℓ ℓ!(2ℓ+1)!!)` and `b_ℓ = 4ℓ(ℓ+2)/((ℓ−1)(ℓ+1)!(2ℓ+1)!!)`.
/// `a₂ = 1/5` is the quadrupole formula's `G/(5c⁵)`.
pub fn energy_flux_coefficients(l: u32) -> (f64, f64) {
    let factorial = |n: u32| (1..=n).map(f64::from).product::<f64>();
    let double_factorial = |n: u32| (1..=n).rev().step_by(2).map(f64::from).product::<f64>();
    let l_f = f64::from(l);
    let a = (l_f + 1.0) * (l_f + 2.0) / ((l_f - 1.0) * l_f * factorial(l) * double_factorial(2 * l + 1));
    let b = 4.0 * l_f * (l_f + 2.0) / ((l_f - 1.0) * factorial(l + 1) * double_factorial(2 * l + 1));
    (a, b)
}

/// Peters' enhancement of the radiated power by eccentricity `e`:
/// `(1 + 73e²/24 + 37e⁴/96)/(1 − e²)^(7/2)`.
pub fn peters_enhancement(e: f64) -> f64 {
    let e2 = e * e;
    (1.0 + 73.0 / 24.0 * e2 + 37.0 / 96.0 * e2 * e2) / (1.0 - e2).powf(3.5)
}

/// Orbital-period derivative (dimensionless, s/s) of a binary of masses `mp` and `mc` (kg),
/// period `pb` (s) and eccentricity `e`, from gravitational radiation:
/// `−(192π/5)(G/c³)^(5/3)(P_b/2π)^(−5/3) f(e) m_p m_c (m_p+m_c)^(−1/3)`.
pub fn peters_pbdot(mp: f64, mc: f64, pb: f64, e: f64) -> f64 {
    -192.0 * PI / 5.0 * (G / C.powi(3)).powf(5.0 / 3.0) * (pb / (2.0 * PI)).powf(-5.0 / 3.0) * peters_enhancement(e) * mp * mc
        / (mp + mc).cbrt()
}

/// The Shapiro range `r = Gm_c/c³` (s) of a binary pulsar's companion of mass `mc` (kg), the
/// post-Keplerian parameter that sets the size of the delay (Damour & Deruelle 1986).
pub fn shapiro_range(mc: f64) -> f64 {
    G * mc / C.powi(3)
}

/// The Shapiro shape `s = sin i` that general relativity predicts for pulsar and companion masses
/// `mp`, `mc` (kg), orbital period `pb` (s) and projected semi-major axis `x` (light-seconds):
/// `s = x n_b^(2/3) (GM/c³)^(−1/3) M/m_c`, with `n_b = 2π/P_b` and `M = m_p + m_c`.
pub fn shapiro_shape(mp: f64, mc: f64, pb: f64, x: f64) -> f64 {
    let m = mp + mc;
    x * (2.0 * PI / pb).powf(2.0 / 3.0) * (G * m / C.powi(3)).powf(-1.0 / 3.0) * m / mc
}

/// The advance of periastron per orbit, `6πGM/(c²a(1 − e²))` (rad), for semi-major axis `a` (m)
/// and eccentricity `e` about total mass parameter `gm`: the first term of Will's (2014) Eq. 64
/// with `γ = β = 1`, and the post-Keplerian `ω̇` once divided by the period. 42.98″ per century for
/// Mercury; 16.90° per year for the double pulsar.
pub fn periastron_advance(gm: f64, a: f64, e: f64) -> f64 {
    6.0 * PI * gm / (C * C * a * (1.0 - e * e))
}

/// The geodetic (de Sitter) precession, `(3/2)(GM)^(3/2)/(c²a^(5/2))` (rad/s), of a gyroscope, a
/// pulsar's spin or the Earth–Moon system's orbit on a circular orbit of radius `a` (m) about `gm`:
/// 6.6″ per year for Gravity Probe B, 1.9″ per century for the Moon's orbit round the Sun.
pub fn geodetic_precession(gm: f64, a: f64) -> f64 {
    1.5 * gm.powf(1.5) / (C * C * a.powf(2.5))
}

/// The Shapiro delay (s) at eccentric anomaly `u` of an orbit with eccentricity `e` and longitude
/// of periastron `omega` (rad): `Δ_S = −2r ln Λ_u`, with
/// `Λ_u = 1 − e cos u − s[sin ω(cos u − e) + √(1 − e²) cos ω sin u]`
/// ([Kramer et al. 2021](https://arxiv.org/abs/2112.06795), after Damour & Deruelle 1986).
/// On a circular orbit it peaks at superior conjunction at `2r z_s`, with `z_s = −ln(1 − s)`.
pub fn shapiro_delay(r: f64, s: f64, e: f64, omega: f64, u: f64) -> f64 {
    -2.0 * r * lambda_u(s, e, omega, u).ln()
}

/// `Λ_u` of [`shapiro_delay`].
fn lambda_u(s: f64, e: f64, omega: f64, u: f64) -> f64 {
    let (sw, cw) = omega.sin_cos();
    let (su, cu) = u.sin_cos();
    1.0 - e * cu - s * (sw * (cu - e) + (1.0 - e * e).sqrt() * cw * su)
}

/// The retardation `δΛ_ret`: the velocity terms inside the logarithm of Kopeikin & Schäfer's Eq. 130
/// ([gr-qc/9902030](https://arxiv.org/abs/gr-qc/9902030)), with their `x m_p/m_c` written as the
/// companion's projected semi-major axis `x_b`.
fn retardation(s: f64, e: f64, omega: f64, u: f64, x_b: f64, pb: f64) -> f64 {
    let (sw, cw) = omega.sin_cos();
    let (su, cu) = u.sin_cos();
    let (k, q) = (2.0 * PI * x_b / pb, (1.0 - e * e).sqrt());
    k / s * e * su - k * s / q * (sw * (cu - e) + q * cw * su) * (e * cw + ((cu - e) * cw - q * sw * su) / (1.0 - e * cu))
}

/// The Shapiro delay (s) to next-to-leading order,
/// `Δ_S = −2r ln(Λ_u + q_NLO(δΛ_len + δΛ_ret))` ([Kramer et al. 2021](https://arxiv.org/abs/2112.06795),
/// Eqs. 26, 32–33), for the orbit of [`shapiro_delay`] with projected semi-major axes `x` (pulsar)
/// and `x_b` (companion), in light-seconds, and orbital period `pb` (s). `q_nlo` scales both
/// corrections: 1 in general relativity, 0 to recover [`shapiro_delay`]. Kramer fit
/// `q_NLO = 1.15(13)` (Eq. 53), scaling this and [`aberration_bending`] together; the retardation
/// alone gives 1.32(24) in their summary table, 1.32(28) in their Eq. 55.
///
/// - Lensing, the bent path: `δΛ_len = 2rs/(x + x_b)`, Kramer's adaptation of Klioner & Zschocke's
///   Eq. 73 ([arXiv:1001.2133](https://arxiv.org/abs/1001.2133)), good to a few per cent outside
///   the eclipse. It shortens the delay, by 0.757 µs at the Double Pulsar's conjunction.
/// - Retardation, the companion moving while the signal crosses the system: `δΛ_ret`, the velocity
///   terms of Kopeikin & Schäfer's Eq. 130, as Kramer use it. Their text gives the focal parameter as
///   `a_R(1 − e²)^{1/2}`, but their Eq. 129 (and so Eq. 130) uses `a_R(1 − e²)`. For the Double
///   Pulsar it moves the delay by −1.082 to +1.189 µs within 0.2 rad of conjunction.
///
/// Against a direct integration of the delay of a ray past the moving companion (ours), Eq. 130 is
/// right to first order in its velocity. Like Kramer, this leaves out the factor `(1 − k·v_c)`
/// before the logarithm (KS99's Eq. 118), 1.2% of the retardation (13 ns for the Double Pulsar).
/// The terms of second order in the velocity reach about 100 ns at its conjunction, where
/// `Λ_u ≈ 6.6×10⁻⁵` is only about a hundred times `v²/2`. Refitting `s` absorbs most of that:
/// Eq. 130 then misses by 14 ns, and 6 ns with the prefactor.
#[allow(clippy::too_many_arguments)]
pub fn shapiro_delay_nlo(r: f64, s: f64, e: f64, omega: f64, u: f64, x: f64, x_b: f64, pb: f64, q_nlo: f64) -> f64 {
    let lensing = 2.0 * r * s / (x + x_b);
    -2.0 * r * (lambda_u(s, e, omega, u) + q_nlo * (lensing + retardation(s, e, omega, u, x_b, pb))).ln()
}

/// The aberration amplitude `𝒜 = Px/(P_b√(1 − e²)s²)` (s) of a pulsar of spin period `period`
/// (s), with its spin along the orbital angular momentum ([Kramer et al. 2021](https://arxiv.org/abs/2112.06795),
/// Eq. 29): the leading aberration delay is `𝒜[sin ψ + e sin ω]`. 3.65 µs for the Double Pulsar.
pub fn aberration_amplitude(period: f64, pb: f64, x: f64, e: f64, s: f64) -> f64 {
    period * x / (pb * (1.0 - e * e).sqrt() * s * s)
}

/// The aberration delay's lensing term `q_NLO 𝒟 cos ψ/Λ_u` (s), with `𝒟 = (P/π)r/(x + x_b)`
/// (Kramer et al. 2021, Eqs. 30, 31 and 34) and `ψ = ω + θ`, θ the true anomaly, for the orbit of
/// [`shapiro_delay_nlo`] and spin period `period` (s): the signal's bending changes its angle of
/// emission. Kramer scale it with the same `q_NLO`, since it cannot be told from retardation; their
/// timing code also moves ψ and `Λ_u` to where the companion is when the signal passes it, which
/// this leaves out. ±1.30 µs for the Double Pulsar.
#[allow(clippy::too_many_arguments)]
pub fn aberration_bending(period: f64, r: f64, s: f64, e: f64, omega: f64, u: f64, x: f64, x_b: f64, q_nlo: f64) -> f64 {
    let theta = 2.0 * ((1.0 + e).sqrt() * (u / 2.0).sin()).atan2((1.0 - e).sqrt() * (u / 2.0).cos());
    q_nlo * period / PI * r / (x + x_b) * (omega + theta).cos() / lambda_u(s, e, omega, u)
}

#[cfg(test)]
mod tests {
    use super::*;
    use metrology::{MPC, SOLAR_MASS};

    fn close(a: f64, b: f64, rel: f64) -> bool {
        (a / b - 1.0).abs() < rel
    }

    #[test]
    fn luminosity_and_strain() {
        assert!(close(PLANCK_LUMINOSITY, 3.62825e52, 1e-5) && close(PLANCK_LUMINOSITY / 4.0, 9.07064e51, 1e-5));
        assert!(close(3.7e49 / PLANCK_LUMINOSITY, 1.01977e-3, 1e-5)); // GW190521
        assert!(close(strain_from_luminosity(1e-3, 250.0, 410.0 * MPC), 1.34931e-21, 1e-5));
        assert!(close(strain_from_luminosity(1.0, 1.0, C), 0.450158, 1e-5));
        assert!(close(circular_flux(0.25, 1.0 / 6.0), 5.14403e-5, 1e-5));
    }

    #[test]
    fn multipole_coefficients() {
        for (l, a, b) in [(2, 1.0 / 5.0, 16.0 / 45.0), (3, 1.0 / 189.0, 1.0 / 84.0), (4, 1.0 / 9072.0, 4.0 / 14175.0)] {
            let (al, bl) = energy_flux_coefficients(l);
            assert!(close(al, a, 1e-12) && close(bl, b, 1e-12), "{l}: {al} {bl}");
        }
    }

    #[test]
    fn binary_pulsars() {
        assert!(close(peters_enhancement(0.6171340), 11.8568, 1e-5));
        let day = 86_400.0;
        // Masses quoted to 10⁻³: Weisberg & Huang's −2.40263(5)×10⁻¹² to that precision.
        let b1913 = peters_pbdot(1.438 * SOLAR_MASS, 1.390 * SOLAR_MASS, 0.322997448918 * day, 0.6171340);
        assert!(close(b1913, -2.40210e-12, 1e-5) && close(b1913, -2.40263e-12, 1e-3));
        let j0737 = peters_pbdot(1.338185 * SOLAR_MASS, 1.248868 * SOLAR_MASS, 0.1022515592973 * day, 0.087777023);
        assert!(close(j0737, -1.24781e-12, 1e-5));
    }

    #[test]
    fn periastron_advances_and_geodetic_precession() {
        use metrology::{AU, GM_SUN, YEAR};
        let arcsec = PI / 648_000.0;
        // Mercury, with JPL's mean elements at J2000 (a = 0.38709927 AU, e = 0.20563593): 42.98″
        // per century, Will's (2014) Eq. 65.
        let (a, e) = (0.387_099_27 * AU, 0.205_635_93);
        let period = 2.0 * PI * (a * a * a / GM_SUN).sqrt();
        assert!(close(periastron_advance(GM_SUN, a, e) * 100.0 * YEAR / period / arcsec, 42.98, 5e-4));
        // The double pulsar (Kramer et al. 2021, Table IV): P_b = 0.1022515592973 d, e = 0.087777023,
        // total mass 2.587052 M☉. The 1PN advance is 16.899367°/yr; they measure 16.899323(13),
        // and the difference is their 2PN term (4.39×10⁻⁴) less the Lense–Thirring one (4.83×10⁻⁴).
        let (pb, e, m) = (0.1022515592973 * 86_400.0, 0.087777023, 2.587052 * GM_SUN);
        let a = (m * (pb / (2.0 * PI)).powi(2)).cbrt();
        let omega_dot = periastron_advance(m, a, e) / pb * YEAR * 180.0 / PI;
        assert!(close(omega_dot, 16.899_367, 1e-6) && close(omega_dot - 16.899_323, 4.4e-5, 0.05), "{omega_dot}");
        // S2 round Sgr A* (GRAVITY 2020, Table E.1: 4.261×10⁶ M☉, 8246.7 pc, a = 125.058 mas,
        // e = 0.884649): 12.16′ per orbit against their 12.1′ (Eq. 1), measured 1.10 ± 0.19 of it.
        let s2 = periastron_advance(4.261e6 * GM_SUN, 125.058e-3 * 8246.7 * AU, 0.884_649) / arcsec / 60.0;
        assert!(close(s2, 12.156, 1e-3) && close(s2, 12.1, 0.06));
        // The Earth–Moon system as a gyroscope round the Sun: "about 2 arcseconds per century"
        // (Will 2014); Gravity Probe B at 642 km (Everitt et al. 2011): 6606.1 mas/yr predicted, with
        // the Earth's oblateness and the solar term, to 0.5% from the bare formula.
        assert!(close(geodetic_precession(GM_SUN, AU) * 100.0 * YEAR / arcsec, 1.919, 1e-3));
        let mas_per_year = YEAR * 180.0 / PI * 3.6e6;
        assert!(close(geodetic_precession(metrology::GM_EARTH, 6371e3 + 642e3) * mas_per_year / 6606.1, 1.0, 6e-3));
        // A circular orbit advances 6πGM/(c²a), and the advance diverges as e → 1 at fixed a.
        assert!(close(periastron_advance(GM_SUN, AU, 0.0), 6.0 * PI * GM_SUN / (C * C * AU), 1e-15));
        assert!(periastron_advance(GM_SUN, AU, 0.99) > 50.0 * periastron_advance(GM_SUN, AU, 0.0));
    }

    #[test]
    fn the_double_pulsar_shapiro_delay() {
        // Kramer et al. (2021), Table 3: r = 6.162(21) µs against GR's Gm_B/c³ = 6.1513 µs for
        // m_B = 1.248868 M☉, 0.5σ apart.
        let (ma, mb) = (1.338185 * SOLAR_MASS, 1.248868 * SOLAR_MASS);
        let r = shapiro_range(mb);
        assert!((r / 6.1513e-6 - 1.0).abs() < 1e-5 && ((6.162e-6 - r) / 0.021e-6 - 0.51).abs() < 0.01);
        // s = 1 − e^(−z_s) = 0.999936 for z_s = 9.65. GR's s from the masses (m_A = 1.338185),
        // P_b and x gives s_obs/s_GR = 1.000013, inside the paper's own test 1.00009(18).
        let measured = 1.0 - (-9.65f64).exp();
        assert!((measured - 0.999936).abs() < 1e-6);
        let s = shapiro_shape(ma, mb, 0.1022515592973 * 86_400.0, 1.415028603);
        assert!((measured / s - 1.00009).abs() < 0.00018 && (measured / s - 1.000013).abs() < 1e-6, "{s}");
        // On a circular orbit the delay peaks at superior conjunction (ω + u = π/2) at 2r z_s.
        let peak = shapiro_delay(r, measured, 0.0, 0.3, PI / 2.0 - 0.3);
        assert!((peak / (2.0 * r * 9.65) - 1.0).abs() < 1e-12 && peak > shapiro_delay(r, measured, 0.0, 0.3, 1.0));
    }

    /// The Double Pulsar (Kramer et al. 2021, Table 3): `(P_b, x, x_B, e, ω, r, s)`, with the
    /// measured r = 6.162 µs and z_s = 9.65, and x_B = x m_A/m_B.
    fn double_pulsar() -> (f64, f64, f64, f64, f64, f64, f64) {
        let x = 1.415028603;
        (0.1022515592973 * 86_400.0, x, x * 1.338185 / 1.248868, 0.087777023, 204.753686f64.to_radians(), 6.162e-6, 1.0 - (-9.65f64).exp())
    }

    /// Eccentric anomalies within ±0.2 rad of superior conjunction (ω + θ = π/2), 40 001 of them.
    fn near_conjunction(e: f64, omega: f64) -> impl Iterator<Item = f64> {
        let f = PI / 2.0 - omega;
        let u = 2.0 * ((1.0 - e).sqrt() * (f / 2.0).sin()).atan2((1.0 + e).sqrt() * (f / 2.0).cos());
        (0..=40_000).map(move |i| u - 0.2 + f64::from(i) * 1e-5)
    }

    #[test]
    fn the_double_pulsar_at_next_to_leading_order() {
        let (pb, x, xb, e, om, r, s) = double_pulsar();
        let spin = 1.0 / 44.05406864196281;
        // Kramer's Eq. 29: 𝒜 ≃ 3.65 µs (3.6503).
        let a = aberration_amplitude(spin, pb, x, e, s);
        assert!((a / 3.65e-6 - 1.0).abs() < 1e-3 && (a - 3.6503e-6).abs() < 5e-11, "{a}");
        // q_NLO = 0 is the leading order; at e = 0 the retardation is −(2πx_B/P_b)s sin ψ cos ψ,
        // ψ = ω + u (from KS99's Eq. 130).
        for u in [0.0, 1.0, 2.5, 4.0] {
            assert_eq!(shapiro_delay_nlo(r, s, e, om, u, x, xb, pb, 0.0), shapiro_delay(r, s, e, om, u));
            let circular = -2.0 * PI * xb / pb * s * (om + u).sin() * (om + u).cos();
            assert!((retardation(s, 0.0, om, u, xb, pb) - circular).abs() < 1e-18);
        }
        // Near conjunction (the research's regression values, gap4_shapiro_nlo.out): lensing
        // −0.7569 µs where Λ_u is least; retardation −1.0820/+1.1885 µs; both −1.4925/+0.8357 µs;
        // the aberration bending −1.3002/+1.2977 µs, with 𝒟 = 15.1891 ns.
        let lo = |u| shapiro_delay(r, s, e, om, u);
        let u_min = near_conjunction(e, om).min_by(|&a, &b| lambda_u(s, e, om, a).total_cmp(&lambda_u(s, e, om, b))).unwrap();
        let lensing = -2.0 * r * (lambda_u(s, e, om, u_min) + 2.0 * r * s / (x + xb)).ln() - lo(u_min);
        assert!((lensing + 0.7569e-6).abs() < 5e-11, "{lensing}");
        let range = |f: &dyn Fn(f64) -> f64| near_conjunction(e, om).map(f).fold((f64::MAX, f64::MIN), |(a, b), v| (a.min(v), b.max(v)));
        let ret = range(&|u| -2.0 * r * (lambda_u(s, e, om, u) + retardation(s, e, om, u, xb, pb)).ln() - lo(u));
        let both = range(&|u| shapiro_delay_nlo(r, s, e, om, u, x, xb, pb, 1.0) - lo(u));
        let bend = range(&|u| aberration_bending(spin, r, s, e, om, u, x, xb, 1.0));
        for ((got_lo, got_hi), (want_lo, want_hi)) in [(ret, (-1.0820, 1.1885)), (both, (-1.4925, 0.8357)), (bend, (-1.3002, 1.2977))] {
            assert!((got_lo * 1e6 - want_lo).abs() < 5e-5 && (got_hi * 1e6 - want_hi).abs() < 5e-5, "{got_lo} {got_hi}");
        }
        let d = aberration_bending(spin, r, s, e, om, 1.0, x, xb, 1.0) * lambda_u(s, e, om, 1.0);
        let theta = 2.0 * ((1.0 + e).sqrt() * 0.5f64.sin()).atan2((1.0 - e).sqrt() * 0.5f64.cos());
        assert!((d / (om + theta).cos() - 15.1891e-9).abs() < 5e-14);
    }

    #[test]
    fn retardation_is_the_moving_companion_to_first_order() {
        // An independent check of KS99's Eq. 130. A uniformly moving mass m delays a straight ray by
        // 2mγ²(1 − k·v)² ∫dτ/R (the boosted linearised Schwarzschild field, h_μν k^μ k^ν =
        // 4mγ²(1 − k·v)²/R), with R² = d² + γ²(v·d)² the distance in the mass's rest frame,
        // d = r + (k − v)τ and r = x_p − x_c at emission (KS99, after Eq. 117). The orbit is KS99's
        // Eqs. 120–125 with Ω = 0 and k = (0, 0, −1). With the velocities scaled by λ, Eq. 130 and its
        // prefactor (1 − k·v) (KS99's Eq. 118) must match the integral, up to a constant and a multiple
        // of k·v that timing absorbs (KS99 after Eq. 116), to second order in λ: the misfit is 0.98%,
        // 0.29% and 0.098% of the effect at λ = 0.1, 0.03 and 0.01 (ours). Without retardation it
        // would be 100%.
        let (pb, x, xb, e, om, r, s) = double_pulsar();
        let (ci, ar, lam) = ((1.0 - s * s).sqrt(), (x + xb) / s, 0.03);
        let p = [om.cos(), om.sin() * ci, om.sin() * s];
        let q = [-om.sin(), om.cos() * ci, om.cos() * s];
        let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        let delay = |rv: [f64; 3], v: [f64; 3]| {
            let (g2, w) = (1.0 / (1.0 - dot(v, v)), [-v[0], -v[1], -1.0 - v[2]]);
            let inv_r = |t: f64| {
                let d: [f64; 3] = core::array::from_fn(|j| rv[j] + w[j] * t);
                1.0 / (dot(d, d) + g2 * dot(v, d).powi(2)).sqrt()
            };
            2.0 * r * g2 * (1.0 + v[2]).powi(2) * crate::integrate_log(1e-9, 1e7, 200_000, inv_r)
        };
        let f = PI / 2.0 - om;
        let u_sc = 2.0 * ((1.0 - e).sqrt() * (f / 2.0).sin()).atan2((1.0 + e).sqrt() * (f / 2.0).cos());
        let us = (-10..=10).map(|i| u_sc + 0.02 * f64::from(i)).chain((0..7).map(|i| u_sc + 0.5 + f64::from(i) * (2.0 * PI - 1.0) / 6.0));
        let (mut kv, mut effect, mut misfit) = (vec![], vec![], vec![]);
        for u in us {
            let f = 2.0 * ((1.0 + e).sqrt() * (u / 2.0).sin()).atan2((1.0 - e).sqrt() * (u / 2.0).cos());
            let rv: [f64; 3] = core::array::from_fn(|j| ar * (1.0 - e * u.cos()) * (p[j] * f.cos() + q[j] * f.sin()));
            let speed = -lam * xb / (x + xb) * 2.0 * PI / pb * ar / (1.0 - e * e).sqrt();
            let v: [f64; 3] = core::array::from_fn(|j| speed * (-p[j] * f.sin() + q[j] * (f.cos() + e)));
            let lu = lambda_u(s, e, om, u);
            assert!(((dot(rv, rv).sqrt() - rv[2]) / ar - lu).abs() < 1e-12);
            let moving = delay(rv, v) - delay(rv, [0.0; 3]);
            kv.push(-v[2]);
            effect.push(moving);
            misfit.push(moving + 2.0 * r * (1.0 + v[2]) * (lu + retardation(s, e, om, u, xb, pb / lam)).ln() - 2.0 * r * lu.ln());
        }
        // The largest residual once a + b k·v is fitted by least squares.
        let detrended = |y: &[f64]| {
            let n = y.len() as f64;
            let (sx, sy, sxx, sxy) = kv.iter().zip(y).fold((0.0, 0.0, 0.0, 0.0), |(a, b, c, d), (k, y)| (a + k, b + y, c + k * k, d + k * y));
            let slope = (n * sxy - sx * sy) / (n * sxx - sx * sx);
            let icpt = (sy - slope * sx) / n;
            kv.iter().zip(y).map(|(k, y)| (y - icpt - slope * k).abs()).fold(0.0, f64::max)
        };
        let ratio = detrended(&misfit) / detrended(&effect);
        assert!(ratio < 5e-3 && detrended(&effect) > 2.5e-8, "{ratio}");
    }
}
