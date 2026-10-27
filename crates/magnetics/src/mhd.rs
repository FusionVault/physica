//! The numbers that decide how a magnetised plasma or conducting fluid behaves: Alfvén and
//! magnetosonic speeds, plasma beta, magnetic pressure, diffusion and the Reynolds, Hartmann,
//! Prandtl and Elsasser numbers, the kinetic scales where MHD ends, and the jump across a
//! perpendicular shock.
//!
//! These are ideal and resistive MHD in SI units. The Lundquist number, Sweet–Parker and
//! plasmoid reconnection rates live in `plasmoid::relaxation`.
//!
//! ```
//! use magnetics::mhd::{alfven_speed, mass_density};
//! use metrology::M_P;
//!
//! // PlasmaPy's example: 0.014 T and 5×10¹⁹ protons per m³ give 43 173.870 m/s.
//! let v = alfven_speed(0.014, mass_density(5e19, M_P, 1.0));
//! assert!((v - 43_173.870).abs() < 1e-3);
//! ```

use core::f64::consts::PI;

use metrology::{C, E_CHARGE, EPS_0, EV, M_E, MU_0};

/// Mass density (kg/m³) of a plasma of `n_i` ions per m³ of mass `m_i` (kg) and charge state `z`,
/// counting the `z` electrons per ion: `n_i(m_i + z m_e)`, as PlasmaPy does.
pub fn mass_density(n_i: f64, m_i: f64, z: f64) -> f64 {
    n_i * (m_i + z * M_E)
}

/// Alfvén speed `B/√(μ₀ρ)`, m/s, for field `b` (T) and mass density `rho` (kg/m³).
pub fn alfven_speed(b: f64, rho: f64) -> f64 {
    b / (MU_0 * rho).sqrt()
}

/// The fast and slow magnetosonic speeds `(v_f, v_s)` at angle `theta` (radians) to the field, from
/// the Alfvén speed `va` and sound speed `cs`:
/// `v² = ½[(v_A² + c_s²) ± √((v_A² + c_s²)² − 4v_A²c_s² cos²θ)]`.
pub fn magnetosonic(va: f64, cs: f64, theta: f64) -> (f64, f64) {
    let sum = va * va + cs * cs;
    let root = (sum * sum - 4.0 * va * va * cs * cs * theta.cos().powi(2)).max(0.0).sqrt();
    (((sum + root) / 2.0).sqrt(), ((sum - root) / 2.0).max(0.0).sqrt())
}

/// Magnetic pressure `B²/2μ₀`, Pa (equally, magnetic energy density, J/m³).
pub fn magnetic_pressure(b: f64) -> f64 {
    b * b / (2.0 * MU_0)
}

/// The field (T) whose magnetic pressure, or energy density, is `p` (Pa or J/m³): `√(2μ₀p)`,
/// the inverse of [`magnetic_pressure`]. 1 eV/cm³ needs 6.35 µG; air at 1 atm balances 0.505 T.
pub fn field_for_pressure(p: f64) -> f64 {
    (2.0 * MU_0 * p).sqrt()
}

/// Plasma beta, thermal over magnetic pressure: `2μ₀nk_BT/B²` for density `n` (m⁻³) of particles at
/// temperature `t_ev` (eV) in field `b` (T). Sum species yourself if they differ.
pub fn beta(n: f64, t_ev: f64, b: f64) -> f64 {
    n * t_ev * EV / magnetic_pressure(b)
}

/// Magnetic diffusivity `η = 1/(μ₀σ)`, m²/s, of conductivity `sigma` (S/m).
pub fn magnetic_diffusivity(sigma: f64) -> f64 {
    1.0 / (MU_0 * sigma)
}

/// Magnetic Reynolds number `UL/η`: advection over diffusion of the field. Dynamos need it large.
pub fn magnetic_reynolds(sigma: f64, u: f64, l: f64) -> f64 {
    u * l / magnetic_diffusivity(sigma)
}

/// Hartmann number `BL√(σ/(ρν))`: magnetic over viscous forces in a conducting fluid of density
/// `rho` and kinematic viscosity `nu` (m²/s), across half-width `l` (m).
pub fn hartmann_number(b: f64, l: f64, sigma: f64, rho: f64, nu: f64) -> f64 {
    b * l * (sigma / (rho * nu)).sqrt()
}

/// Magnetic Prandtl number `νμ₀σ`: viscosity over magnetic diffusivity (10⁻⁵ for liquid metals).
pub fn magnetic_prandtl(nu: f64, sigma: f64) -> f64 {
    nu * MU_0 * sigma
}

/// Elsasser number `σB²/(2ρΩ)`: Lorentz over Coriolis force in a rotating conductor (order one
/// in Earth's core).
pub fn elsasser(sigma: f64, b: f64, rho: f64, omega: f64) -> f64 {
    sigma * b * b / (2.0 * rho * omega)
}

/// Free decay time `r²μ₀σ/π²` (s) of the slowest dipole mode of a conducting sphere of radius `r`:
/// how long a field lasts with no dynamo.
pub fn free_decay_time(r: f64, sigma: f64) -> f64 {
    r * r * MU_0 * sigma / (PI * PI)
}

/// Inertial length `c/ω_p` (m) of particles of density `n` (m⁻³), mass `m` (kg) and charge state
/// `z`: below the ion value, ions decouple from the field and Hall MHD takes over. For electrons,
/// `ω_p` is `2π` times `gravitational::hf::plasma_frequency`, which `physica` checks.
pub fn inertial_length(n: f64, m: f64, z: f64) -> f64 {
    C / (n * (z * E_CHARGE).powi(2) / (EPS_0 * m)).sqrt()
}

/// Thermal gyroradius `√(2k_BT/m)/Ω` (m) of particles of temperature `t_ev` (eV), mass `m` and
/// charge state `z` in field `b` (T).
pub fn gyroradius(t_ev: f64, b: f64, m: f64, z: f64) -> f64 {
    (2.0 * t_ev * EV / m).sqrt() * m / (z * E_CHARGE * b)
}

/// The frequency (Hz) at which a flow of speed `v` carries past a spacecraft the larger of the
/// proton inertial length and gyroradius: where the solar-wind spectrum breaks at ion scales
/// ([Chen et al. 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC4459186/)). Density `n` (m⁻³),
/// proton temperature `t_ev` (eV), field `b` (T).
pub fn ion_break_frequency(v: f64, n: f64, t_ev: f64, b: f64) -> f64 {
    let scale = inertial_length(n, metrology::M_P, 1.0).max(gyroradius(t_ev, b, metrology::M_P, 1.0));
    v / (2.0 * PI * scale)
}

/// Attraction per unit length (N/m) between two parallel currents `i1`, `i2` (A) a distance `r`
/// (m) apart: `μ₀I₁I₂/(2πr)`. It falls off as `1/r`, more slowly than gravity's `1/r²` between
/// point masses, which is the basis of the claim that filamentary currents could do the work
/// attributed to dark matter.
///
/// The comparison is worth making in numbers. Two parallel gas filaments of line mass `λ` attract
/// gravitationally at `2Gλ₁λ₂/r`, also as `1/r`, so neither wins by geometry: what decides it is
/// the coefficient. At galactic scale the gravitational term dominates by about twelve orders of
/// magnitude even for mega-ampere currents, as this module's test shows.
pub fn line_current_force(i1: f64, i2: f64, r: f64) -> f64 {
    MU_0 * i1 * i2 / (2.0 * PI * r)
}

/// Gravitational attraction per unit length (N/m) between two parallel filaments of line mass
/// `lambda1`, `lambda2` (kg/m) a distance `r` (m) apart: `2Gλ₁λ₂/r`, the gravitational twin of
/// [`line_current_force`].
pub fn line_mass_force(lambda1: f64, lambda2: f64, r: f64) -> f64 {
    2.0 * metrology::G * lambda1 * lambda2 / r
}

/// Density compression `X = ρ₂/ρ₁` across a perpendicular MHD shock of Alfvén Mach number `m_a`,
/// upstream beta `beta` and adiabatic index `gamma`: the positive root of
/// `(2−γ)X² + [γ(1+β) + (γ−1)M_A²]X − (γ+1)M_A² = 0`, which tends to `(γ+1)/(γ−1)` (4) for
/// strong shocks. `None` unless the flow is faster than the fast speed, `M_A > √(1 + γβ/2)`.
pub fn perpendicular_shock(m_a: f64, beta: f64, gamma: f64) -> Option<f64> {
    if m_a * m_a <= 1.0 + gamma * beta / 2.0 {
        return None;
    }
    let (a, b, c) = (2.0 - gamma, gamma * (1.0 + beta) + (gamma - 1.0) * m_a * m_a, -(gamma + 1.0) * m_a * m_a);
    Some(if a.abs() < 1e-12 { -c / b } else { (-b + (b * b - 4.0 * a * c).sqrt()) / (2.0 * a) })
}

#[cfg(test)]
mod tests {
    use super::*;
    use metrology::{AMU, M_P};

    fn close(a: f64, b: f64, rel: f64) -> bool {
        (a / b - 1.0).abs() < rel
    }

    #[test]
    fn field_for_pressure_inverts_magnetic_pressure() {
        // The magnetics report's equipartition row: 1 eV/cm³ needs 6.35 µG; one atmosphere
        // balances 0.505 T (the plasmoids report's confinement cap).
        assert!(close(field_for_pressure(metrology::MEV), 6.35e-10, 1e-3));
        assert!(close(field_for_pressure(101_325.0), 0.5046, 1e-4));
        for b in [1e-9, 1.0, 45.0] {
            assert!(close(field_for_pressure(magnetic_pressure(b)), b, 1e-15));
        }
    }

    #[test]
    fn alfven_speeds() {
        // PlasmaPy's docstring values, protons and He²⁺ (α-particle mass, CODATA 2018).
        assert!(close(alfven_speed(0.014, mass_density(5e19, M_P, 1.0)), 43_173.870, 1e-7));
        assert!(close(alfven_speed(0.014, mass_density(5e19, 6.644_657_335_7e-27, 2.0)), 21_664.18, 1e-6));
        // A coronal loop, 11 G and 1.67×10⁻¹² kg/m³: 759 km/s.
        assert!(close(alfven_speed(11e-4, 1.67e-12), 759.3e3, 1e-3));
        // ITER: 5.3 T, 1.01×10²⁰ m⁻³ of D–T (2.5 u).
        assert!(close(alfven_speed(5.3, 1.01e20 * 2.5 * AMU), 7.30156e6, 1e-5));
    }

    #[test]
    fn waves_pressure_and_beta() {
        let (f, s) = magnetosonic(3.0, 1.0, 0.0);
        assert!((f - 3.0).abs() < 1e-12 && (s - 1.0).abs() < 1e-12);
        let (f, s) = magnetosonic(3.0, 1.0, core::f64::consts::FRAC_PI_2);
        assert!((f - 10f64.sqrt()).abs() < 1e-12 && s.abs() < 1e-7);
        for (b, p) in [(1.0, 397_887.357_5), (12.2, 5.922_155e7), (45.0, 8.057_219e8)] {
            assert!(close(magnetic_pressure(b), p, 1e-7));
        }
        // Chen et al.'s low-β solar wind: 3.9 cm⁻³, 1.64 eV protons, 13.3 nT.
        assert!(close(beta(3.9e6, 1.64, 13.3e-9), 0.014_559_8, 1e-5));
    }

    #[test]
    fn dimensionless_numbers() {
        // Liquid sodium: μ₀σ ≈ 10 s/m², so Rm = 100 needs UL = 10 m²/s.
        let sigma = 10.0 / MU_0;
        assert!(close(magnetic_reynolds(sigma, 1.0, 10.0), 100.0, 1e-12));
        // S = Ha·√Pm when the speed is the Alfvén speed.
        let (b, l, rho, nu) = (0.5, 0.1, 930.0, 7e-7);
        let s = alfven_speed(b, rho) * l / magnetic_diffusivity(sigma);
        assert!(close(s, hartmann_number(b, l, sigma, rho, nu) * magnetic_prandtl(nu, sigma).sqrt(), 1e-12));
        // Earth's core, illustrative inputs: Λ ≈ 2.5, free decay 49 kyr.
        assert!(close(elsasser(1e6, 2e-3, 1.1e4, 7.292e-5), 2.493, 1e-3));
        assert!(close(free_decay_time(3.48e6, 1e6) / (365.25 * 86_400.0), 48.9e3, 2e-3));
    }

    #[test]
    fn kinetic_scales() {
        // PlasmaPy: 5 m⁻³ of He⁺ and of electrons; the proton practical form 2.28×10⁷ cm/√(n/cm⁻³).
        assert!(close(inertial_length(5.0, 4.002_602 * AMU - M_E, 1.0), 2.029_86e8, 1e-5));
        assert!(close(inertial_length(5.0, M_E, 1.0), 2.376_53e6, 1e-5));
        assert!(close(inertial_length(1e6, M_P, 1.0) * 100.0, 2.277_11e7, 1e-5));
        // Chen et al. 2014, Table 1 means: predicted breaks 0.635 Hz (low β; observed 0.58 ± 0.19)
        // and 0.223 Hz (high β; observed 0.223 ± 0.068).
        assert!(close(ion_break_frequency(460e3, 3.9e6, 1.64, 13.3e-9), 0.635, 2e-3));
        assert!(close(ion_break_frequency(389e3, 17.1e6, 4.8, 1.14e-9), 0.223, 5e-3));
    }

    #[test]
    fn perpendicular_shocks() {
        assert!(close(perpendicular_shock(8.0, 1.0, 5.0 / 3.0).unwrap(), 3.615_425, 1e-6));
        assert!(close(perpendicular_shock(8.0, 2.1, 5.0 / 3.0).unwrap(), 3.483_387, 1e-6));
        assert!(close(perpendicular_shock(1e5, 1.0, 5.0 / 3.0).unwrap(), 4.0, 1e-6));
        assert!(perpendicular_shock(1.0, 1.0, 5.0 / 3.0).is_none());
    }

    /// Both line forces fall as 1/r, so the comparison is one of coefficients. At galactic scale
    /// gravity wins by twelve orders of magnitude even for mega-ampere Birkeland currents.
    #[test]
    fn filamentary_currents_against_gravity() {
        let pc = metrology::PARSEC;
        let lambda = metrology::SOLAR_MASS / pc; // one solar mass per parsec
        let magnetic = line_current_force(1e6, 1e6, pc);
        let gravitational = line_mass_force(lambda, lambda, pc);
        assert!(close(magnetic, 6.4815e-12, 1e-4) && close(gravitational, 17.9637, 1e-4));
        assert!(close(magnetic / gravitational, 3.606e-13, 1e-3));
        // Both scale as 1/r, so the ratio is scale-free: no separation makes magnetism win.
        let far = line_current_force(1e6, 1e6, 1e3 * pc) / line_mass_force(lambda, lambda, 1e3 * pc);
        assert!(close(far, magnetic / gravitational, 1e-12));
        // Matching gravity at a parsec would need a current of 1.7e12 A.
        let needed = (gravitational * 2.0 * PI * pc / MU_0).sqrt();
        assert!(close(needed, 1.6648e12, 1e-3), "{needed:e}");
        // In a laboratory pinch the balance is the other way round, as `plasmoid` records.
        assert!(line_current_force(1e6, 1e6, 0.01) > line_mass_force(1e-6, 1e-6, 0.01));
    }
}
