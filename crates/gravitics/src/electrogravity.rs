//! Claims that a pulsed high-voltage or time-varying magnetic field makes a local gravitational
//! field, computed against the classical forces that act in the same apparatus.
//!
//! The claim is in Pais's US Navy patents and in Zhang Xiangqian's US 12,417,870 (2025) with the
//! experiments behind it: an electromagnetic field that changes should make a gravitational, or
//! "vortex gravitational", field. Three rigs recur, all reported by Zhang:
//! - a light flake, or a nylon ring on a balance, around a wire with a gap, one side at 40–55 kV:
//!   the flake moves toward the positive wire, the ring weighs about 1 g less or more;
//! - a polyethylene ball on a cotton thread in a vacuum tank between two open coils, the upper
//!   one live and the lower one grounded: the ball turns;
//! - thin objects hung in the 18 mm air gap of a 0.1–0.5 T core: every material moves.
//!
//! Each rig has classical forces of the right size, which this module computes:
//! - **Polarisation.** A live electrode pulls any dielectric toward it, the dielectrophoretic
//!   (Kelvin) force ([`dielectrophoretic_force`], and [`sphere_force`] for a ball). It follows
//!   `∇E²`, so the field's sign does not matter, but **which electrode is live** does. When one
//!   side is grounded, swapping the leads moves the live side and the body follows it, which reads
//!   as "always toward the positive pole". [`crate::electrodes`] gives the field of the real
//!   geometry; for the nylon ring it predicts lighter with the top wire live, heavier with the
//!   bottom one live, about 0.1–0.7 g against the reported 0.8–1.2 g.
//! - **Charge.** Any ions, from corona in air or a glow discharge in a poor vacuum, charge a body
//!   to the Pauthenier limit ([`field_charge`]), and the field then pulls on that charge.
//! - **Torque.** A polarisation that lags a turning field gives a torque ([`rotation_torque`])
//!   that would wind a cotton thread ([`torsion_constant`]) hundreds of turns.
//! - **Magnetism.** In a magnet's gap, the susceptibility force
//!   (`magnetics::magnets::susceptibility_force`) and eddy currents move every material.
//!
//! Against these sit the standard prediction for the gravity of the field's own energy
//! ([`field_energy_gravity`], about 10⁻²⁷ m/s²), and the theory's own field. Zhang's field
//! `E_θ/e_r = A×R/c²` is the transverse field of an accelerating charge
//! ([`radiation_field`]), with the charge's acceleration read as a gravitational field, so the
//! charges' own acceleration ([`drift_velocity`]) bounds it.
//!
//! Three tests separate a real gravitational field from the classical forces:
//! - [`mass_independence`]: the classical forces scale with a body's volume, permittivity and
//!   charge, while gravity accelerates every body alike;
//! - the live-end test: drive the electrodes symmetrically at ±V/2 and the polarisation force on a
//!   centred body vanishes, while a field made by charges accelerating between the electrodes
//!   would not change;
//! - the screening test: a grounded metal box screens every electric force and no gravitational
//!   one, so an accelerometer boxed at the test body's place must read the claimed field.
//!
//! ```
//! use gravitics::electrogravity::{dielectrophoretic_force, mass_independence};
//!
//! // A 2 cm square of 0.15 mm plastic (ε_r = 2.2) near a 40 kV electrode, where
//! // E ≈ 10⁷ V/m falls off over centimetres: ∇E² ≈ 2×10¹⁶ V²/m³.
//! let volume = 0.02 * 0.02 * 0.15e-3;
//! let f = dielectrophoretic_force(volume, 2.2, 2e16);
//! assert!((f * 1e3 - 6.4).abs() < 0.1); // 6.4 mN: the weight of 650 mg
//!
//! // Polyethylene against tungsten, same shape: a classical force gives 20× the acceleration
//! // on the lighter body, a gravitational field exactly the same on both.
//! let ratio = mass_independence(950.0, 19_300.0);
//! assert!((ratio - 20.3).abs() < 0.1);
//! ```

use core::f64::consts::PI;

use metrology::{C, E_CHARGE, EPS_0, G, H};


/// The dielectrophoretic (Kelvin) force (N) on a neutral dielectric of volume `volume` (m³) and
/// relative permittivity `eps_r` in a field whose squared magnitude has gradient `grad_e2`
/// (V²/m³): `F = ½Vε₀(ε_r − 1)∇(E²)`, the dilute limit.
///
/// It points toward the strongest field. Because it depends on `E²`, flipping the sign of every
/// potential leaves it unchanged; moving which electrode is live does not.
pub fn dielectrophoretic_force(volume: f64, eps_r: f64, grad_e2: f64) -> f64 {
    0.5 * volume * EPS_0 * (eps_r - 1.0) * grad_e2
}

/// The Clausius–Mossotti factor `K = (ε_r − 1)/(ε_r + 2)` of a sphere, whose dipole in a field
/// `E` is `p = 4πε₀a³K E`. It runs from 0 for vacuum to 1 for a conductor.
pub fn clausius_mossotti(eps_r: f64) -> f64 {
    (eps_r - 1.0) / (eps_r + 2.0)
}

/// The dielectrophoretic force (N) on a sphere of radius `radius` (m): `F = 2πε₀a³K∇(E²)`.
/// For `ε_r → 1` it is [`dielectrophoretic_force`] with `V = 4πa³/3`.
pub fn sphere_force(radius: f64, eps_r: f64, grad_e2: f64) -> f64 {
    2.0 * PI * EPS_0 * radius.powi(3) * clausius_mossotti(eps_r) * grad_e2
}

/// The charge (C) a sphere of radius `radius` (m) collects from ions drifting along a field `e`
/// (V/m), the Pauthenier limit `q = 12πε₀a²E ε_r/(ε_r + 2)`; a conductor takes `12πε₀a²E`.
///
/// It needs ions: corona in air, or a glow discharge in a poor vacuum (`plasmoid::discharge`).
pub fn field_charge(radius: f64, eps_r: f64, e: f64) -> f64 {
    12.0 * PI * EPS_0 * radius * radius * e * eps_r / (eps_r + 2.0)
}

/// The largest torque (N m) a field `e` (V/m) can put on a lossy dielectric sphere of radius
/// `radius` (m) through its polarisation: `4πε₀a³E²(1 − K)/2`.
///
/// A static field gives an isotropic sphere no torque, since its dipole lines up with the field.
/// A field that turns while the sphere's charge relaxes leaves the dipole lagging, and the
/// out-of-phase part of the polarisability is at most `(1 − K)/2`, halfway between the
/// insulating limit `K` and the conducting limit 1 (electrorotation; Jones, *Electromechanics of
/// Particles*, 1995, ch. 4). Anisotropy, a seam or surface charge give torques of the same scale
/// times their own small factor.
pub fn rotation_torque(radius: f64, eps_r: f64, e: f64) -> f64 {
    2.0 * PI * EPS_0 * radius.powi(3) * e * e * (1.0 - clausius_mossotti(eps_r))
}

/// The torsion constant (N m/rad) of a thread, from the period `period` (s) of the free twisting
/// oscillation of the body it carries, of moment of inertia `inertia` (kg m²): `κ = 4π²I/T²`.
/// A torque `τ` twists the thread by `τ/κ` radians.
pub fn torsion_constant(inertia: f64, period: f64) -> f64 {
    4.0 * PI * PI * inertia / (period * period)
}

/// The angular momentum (kg m²/s) stored in crossed static fields over a volume `volume` (m³):
/// `L = ε₀ r (E × B) V`, for field magnitudes `e` (V/m) and `b` (T) at lever arm `lever` (m),
/// taken perpendicular.
///
/// Switching the fields off hands it to the apparatus as an impulse (the Feynman disc paradox).
/// It needs a real magnetic field: open coils fed 10 µA give about 10⁻¹⁰ T, far too little.
pub fn field_angular_momentum(e: f64, b: f64, lever: f64, volume: f64) -> f64 {
    EPS_0 * lever * e * b * volume
}

/// The transverse field (V/m) of a charge `charge` (C) accelerating at `acceleration` (m/s²), at
/// distance `distance` (m) and angle `angle` (rad) from the acceleration:
/// `E = qa sin θ/(4πε₀c²r)`.
///
/// This is J. J. Thomson's kink construction (Purcell, *Electricity and Magnetism*, app. B).
/// Zhang's `E_θ/e_r = A×R/c²` is the same field, derived the same way, with the charge's
/// acceleration `A` named a gravitational field. In Maxwell's theory it is an electric field and
/// acts only on charge and polarisation.
pub fn radiation_field(charge: f64, acceleration: f64, distance: f64, angle: f64) -> f64 {
    charge * acceleration * angle.sin() / (4.0 * PI * EPS_0 * C * C * distance)
}

/// Copper's conduction-electron density (m⁻³), one per atom (Kittel, table 6.1).
pub const COPPER_ELECTRONS: f64 = 8.47e28;

/// The drift velocity (m/s) of carriers of density `density` (m⁻³) carrying a current `current`
/// (A) through a cross-section `area` (m²): `v = I/(neS)`.
///
/// Given `dI/dt` in place of `I`, it returns the carriers' acceleration. Zhang's theory takes
/// that acceleration as the gravitational field (`A = −dv/dt`, his Eq. 14), so this bounds the
/// theory's prediction in its own terms: a body can gain at most the carriers' own change of
/// velocity.
pub fn drift_velocity(current: f64, density: f64, area: f64) -> f64 {
    current / (density * E_CHARGE * area)
}

/// The Newtonian pull (m/s²) of field energy `energy` (J) at distance `distance` (m): `GU/(c²r²)`.
/// General relativity's prediction for an electromagnetic field's own gravity, up to a factor of
/// 2 from the field's stresses.
pub fn field_energy_gravity(energy: f64, distance: f64) -> f64 {
    G * energy / (C * C * distance * distance)
}

/// The pressure (Pa) a field `e` (V/m) puts on a conductor's surface: `ε₀E²/2`. At 10⁶ V/m it
/// is 4.4 Pa, enough to shift an interferometer's mirror mounts by micrometres.
pub fn electric_pressure(e: f64) -> f64 {
    0.5 * EPS_0 * e * e
}

/// Zhang's five-constant equation, `4π²r³c²/(T²Ghν)`, for an orbit of radius `radius` (m) and
/// period `period` (s) about a mass whose Compton frequency `Mc²/h` is `frequency` (Hz).
///
/// It equals 1 exactly when Kepler's third law `GM = 4π²r³/T²` holds, so it is that law with the
/// central mass written as a frequency: an identity, not a new constraint linking the constants.
pub fn five_constants(radius: f64, period: f64, frequency: f64) -> f64 {
    4.0 * PI * PI * radius.powi(3) * C * C / (period * period * G * H * frequency)
}

/// The depolarisation factor along the long axis of a prolate spheroid whose long axis is
/// `aspect` times its short ones: `L = (1 − e²)/e³ (atanh e − e)`, with `e² = 1 − 1/aspect²`. A
/// sphere gives 1/3; the two short axes share the rest, `(1 − L)/2` each (Osborn, Phys. Rev. 67,
/// 351, 1945).
pub fn depolarisation(aspect: f64) -> f64 {
    let e2 = 1.0 - 1.0 / (aspect * aspect);
    if e2 < 0.09 {
        // (1 − e²) Σ e²ᵏ/(2k + 3): the closed form loses digits as e → 0.
        (1.0 - e2) * (0..24).map(|k| e2.powi(k) / (2 * k + 3) as f64).sum::<f64>()
    } else {
        let e = e2.sqrt();
        (1.0 - e2) / (e2 * e) * (e.atanh() - e)
    }
}

/// The polarisability per unit volume, in units of ε₀, of a dielectric prolate spheroid of axis
/// ratio `aspect`: `(ε_r − 1)/(1 + L(ε_r − 1))` along its long axis and across it. A sphere
/// gives `3K` both ways ([`clausius_mossotti`]); the difference is what lets a field turn a body
/// that is not quite round.
pub fn spheroid_polarisability(eps_r: f64, aspect: f64) -> (f64, f64) {
    let long = depolarisation(aspect);
    let alpha = |l: f64| (eps_r - 1.0) / (1.0 + l * (eps_r - 1.0));
    (alpha(long), alpha((1.0 - long) / 2.0))
}

/// The torque (N m) a field `e` (V/m) puts on a body of volume `volume` (m³) whose polarisability
/// per unit volume differs by `d_alpha` (units of ε₀) between its long and short axes, when the
/// long axis is `angle` (rad) from the field: `½ε₀VΔα E² sin 2θ`, turning the long axis toward
/// the field.
pub fn alignment_torque(volume: f64, d_alpha: f64, e: f64, angle: f64) -> f64 {
    0.5 * EPS_0 * volume * d_alpha * e * e * (2.0 * angle).sin()
}

/// The attraction (N) between a charge `charge` (C) and a flat dielectric of permittivity
/// `eps_r` at distance `distance` (m), through its image: `q²(ε_r − 1)/(16πε₀d²(ε_r + 1))`
/// (Jackson, *Classical Electrodynamics*, §4.4). A conductor gives `q²/(16πε₀d²)`.
pub fn image_force(charge: f64, distance: f64, eps_r: f64) -> f64 {
    charge * charge * (eps_r - 1.0) / (16.0 * PI * EPS_0 * distance * distance * (eps_r + 1.0))
}

/// Zhang's mass–charge conversion constant `k_mq` (C s/kg) in `q = k_mq dM/dt`, as his papers give
/// it: "experimentally calibrated" (Paper J, §5.5), `eG/c³` (Paper E, Eq. E42), and the
/// "canonical" `eG/(4πc³)` (Paper J, §5.5). They span 45 orders of magnitude.
pub const ZHANG_KMQ: [(&str, f64); 3] = [
    ("calibrated (Paper J)", 7.08e-11),
    ("eG/c³ (Paper E)", E_CHARGE * G / (C * C * C)),
    ("eG/(4πc³) (Paper J)", E_CHARGE * G / (4.0 * PI * C * C * C)),
];

/// Zhang's space–mass coupling constant `μ_g` (kg/m) and the gravitational constant it gives,
/// `G = 4πc²/μ_g` (Papers E and J): 6.685×10⁻¹¹, which the papers call within 0.16% of CODATA.
pub const ZHANG_MU_G: f64 = 1.6895e28;

/// The constant `f` (kg s/C) of Zhang's vortex gravitational field, `∇×A = B/f`, and of
/// `E = −f ∂A/∂t` (*Unified Field Theory*, 7th ed., ch. 34 and 36): his charge and mass constants
/// combined with `4πε₀G`, which in the papers' notation is `k_mq/(4πε₀G)`. The 2025 air-gap paper
/// calls the same constant `k`.
pub fn vortex_constant(k_mq: f64) -> f64 {
    k_mq / (4.0 * PI * EPS_0 * G)
}

/// The angular acceleration (rad/s²) Zhang's vortex field gives any body centred on a field line
/// of a magnetic field `b` (T), for his constant `f` ([`vortex_constant`]). From `∇×A = B/f`, a
/// uniform field makes `A_φ = Br/(2f)`, a rigid spin-up at `B/(2f)` whatever the body is made of.
pub fn vortex_spin(b: f64, f: f64) -> f64 {
    b / (2.0 * f)
}

/// How much more a classical polarisation force accelerates a light test body than a dense one of
/// the same shape: the density ratio `rho_dense/rho_light`.
///
/// A force set by volume and permittivity is the same on both, so the acceleration goes as `1/ρ`.
/// A gravitational field gives `a = E_g` whatever the body is made of, so the ratio is 1. Running
/// the same experiment with two bodies of equal shape and different density therefore separates
/// the two explanations, and it needs no absolute calibration.
pub fn mass_independence(rho_light: f64, rho_dense: f64) -> f64 {
    rho_dense / rho_light
}

/// The acceleration (m/s²) a force `force` (N) gives a body of volume `volume` (m³) and density
/// `rho` (kg/m³): `F/(ρV)`. Pair it with [`dielectrophoretic_force`] to turn a claimed
/// displacement into the acceleration a gravimeter would have to resolve.
pub fn acceleration(force: f64, volume: f64, rho: f64) -> f64 {
    force / (rho * volume)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, rel: f64) -> bool {
        (a / b - 1.0).abs() < rel
    }

    // The patent's ball: 0.35 g of polyethylene (950 kg/m³, ε_r = 2.3), radius 4.447 mm. The
    // field and its gradient at the ball come from `electrodes` (50 kV on the upper coil).
    const BALL_MASS: f64 = 0.35e-3;
    const BALL_RADIUS: f64 = 4.447e-3;
    const FIELD: f64 = 2.4991e5;
    const GRAD_E2: f64 = 1.2753e12;

    #[test]
    fn the_sheet_setup_is_dielectrophoresis() {
        let volume = 0.02 * 0.02 * 0.15e-3;
        let f = dielectrophoretic_force(volume, 2.2, 2e16);
        assert!(close(f, 6.375e-3, 1e-3), "{f:e}");
        // Flipping every potential does not change it: E -> -E leaves grad(E^2) alone.
        assert!(dielectrophoretic_force(volume, 2.2, 2e16) == f);
        // It vanishes in a uniform field and for vacuum-like material.
        assert!(dielectrophoretic_force(volume, 2.2, 0.0) == 0.0);
        assert!(dielectrophoretic_force(volume, 1.0, 2e16) == 0.0);
        // 6.4 mN lifts 650 mg, far more than any claimed gravitational effect.
        assert!(close(f / metrology::G0 * 1e6, 650.0, 1e-2));
    }

    #[test]
    fn the_sphere_formulas() {
        assert!(clausius_mossotti(1.0) == 0.0 && close(clausius_mossotti(1e12), 1.0, 1e-11));
        // Near ε_r = 1 the sphere force is the dilute formula for the sphere's volume.
        let (a, eps): (f64, f64) = (1e-3, 1.0 + 1e-6);
        let dilute = dielectrophoretic_force(4.0 / 3.0 * PI * a.powi(3), eps, 1e12);
        assert!(close(sphere_force(a, eps, 1e12), dilute, 1e-5));
        // The patent's ball sits 2.5x10^5 V/m from the live coil: pulled up by 1.9 uN, 5.5x10^-4
        // of its weight.
        let f = sphere_force(BALL_RADIUS, 2.3, GRAD_E2);
        assert!(close(f, 1.887e-6, 1e-3), "{f:e}");
        assert!(close(f / (BALL_MASS * metrology::G0), 5.50e-4, 1e-2));
    }

    #[test]
    fn ions_charge_the_ball_and_the_field_pulls_on_it() {
        // With a glow in the tank the ball can hold 0.88 nC, and the field then pulls 0.22 mN,
        // 6% of its weight: a swing of 3.7 degrees.
        let q = field_charge(BALL_RADIUS, 2.3, FIELD);
        assert!(close(q, 882.5e-12, 1e-3), "{q:e}");
        let ratio = q * FIELD / (BALL_MASS * metrology::G0);
        assert!(close(ratio, 0.0643, 1e-2), "{ratio}");
        // A conductor takes 12 pi eps0 a^2 E.
        assert!(close(field_charge(1.0, 1e15, 1.0), 12.0 * PI * EPS_0, 1e-12));
    }

    #[test]
    fn the_torque_dwarfs_the_thread() {
        let tau = rotation_torque(BALL_RADIUS, 2.3, FIELD);
        assert!(close(tau, 2.132e-7, 1e-3), "{tau:e}");
        // A conductor cannot lag: K = 1 leaves nothing out of phase.
        assert!(rotation_torque(BALL_RADIUS, 1e15, FIELD) < 1e-20);
        // A cotton thread that lets the ball twist back and forth every 30 s resists 1.2x10^-10
        // N m per radian, so the ceiling would wind it 280 turns; 1% of it, 2.8 turns.
        let inertia = 0.4 * BALL_MASS * BALL_RADIUS * BALL_RADIUS;
        let kappa = torsion_constant(inertia, 30.0);
        assert!(close(kappa, 1.2146e-10, 1e-3), "{kappa:e}");
        assert!(close(tau / kappa / (2.0 * PI), 279.4, 1e-3));
    }

    #[test]
    fn field_angular_momentum_needs_a_real_magnetic_field() {
        // 10^6 V/m across 10^-3 T over 10^-4 m^3 at a 5 cm lever arm.
        let l = field_angular_momentum(1e6, 1e-3, 0.05, 1e-4);
        assert!(close(l, 4.427e-14, 1e-3), "{l:e}");
        assert!(field_angular_momentum(0.0, 1e-3, 0.05, 1e-4) == 0.0);
        assert!(close(field_angular_momentum(2e6, 1e-3, 0.05, 1e-4) / l, 2.0, 1e-12));
        // The patent's open coils carry 10 uA: 7.4x10^-11 T at the ball (magnetics::solenoid_bz),
        // a millionth of Earth's field. Over the whole 10 cm x 20 cm tank that stores 1.3x10^-20
        // kg m^2/s, which would spin the ball at 5x10^-12 rad/s, one turn in 40 000 years.
        let tank = PI * 0.05 * 0.05 * 0.2;
        let l = field_angular_momentum(FIELD, 7.404e-11, 0.05, tank);
        assert!(close(l, 1.2866e-20, 1e-3), "{l:e}");
        let inertia = 0.4 * BALL_MASS * BALL_RADIUS * BALL_RADIUS;
        assert!(l / inertia < 5e-12, "{:e}", l / inertia);
    }

    #[test]
    fn the_theory_in_its_own_terms() {
        // The kink field: 1 nC at 10^12 m/s^2, seen side-on from 1 m: 10^-4 V/m.
        assert!(close(radiation_field(1e-9, 1e12, 1.0, PI / 2.0), 1e-4, 1e-8));
        assert!(radiation_field(1e-9, 1e12, 1.0, 0.0) == 0.0);
        // 10 uA in the patent's 1.5 mm copper: electrons drift at 4.2x10^-10 m/s. That is all
        // the velocity the theory's field can hand on; at the 30 kHz ripple of the supply it moves
        // the ball by 2x10^-15 m.
        let area = PI * 0.75e-3 * 0.75e-3;
        let v = drift_velocity(10e-6, COPPER_ELECTRONS, area);
        assert!(close(v, 4.170e-10, 1e-3), "{v:e}");
        assert!(v / (2.0 * PI * 3e4) < 3e-15);
        // The 55 kV weighing rig's claimed field: 1 g of 0.38 kg, 2.6x10^-2 m/s^2 held for as long
        // as the voltage is. A steady voltage accelerates no charge, so in the theory's own terms
        // the steady part is zero; 30 uA ramped in 1 ms gives 1.3x10^-6 m/s^2 while it lasts.
        let claimed = 1e-3 * metrology::G0 / 0.376;
        assert!(close(claimed, 0.0261, 1e-2));
        assert!(drift_velocity(30e-6 / 1e-3, COPPER_ELECTRONS, area) < 1e-4 * claimed);
    }

    #[test]
    fn general_relativity_for_the_same_rig() {
        // The coils store 7.0 mJ: as mass, 7.7x10^-20 kg, pulling 1.4x10^-27 m/s^2 at 6 cm.
        let g = field_energy_gravity(6.963e-3, 0.06);
        assert!(close(g, 1.436e-27, 1e-3), "{g:e}");
        // Maxwell stress on an interferometer mirror at 10^6 V/m: 4.4 Pa.
        assert!(close(electric_pressure(1e6), 4.427, 1e-3));
    }

    #[test]
    fn five_constants_is_keplers_law() {
        // The Moon: mean distance 384 399 km, sidereal month 27.321 661 days, about the
        // Earth-Moon mass (GM_EARTH x (1 + 1/81.3006)).
        let mass = metrology::GM_EARTH * (1.0 + 1.0 / 81.300_57) / G;
        let nu = mass * C * C / H;
        let ratio = five_constants(384_399e3, 27.321_661 * 86_400.0, nu);
        // 1 to the 0.3% by which the Sun perturbs the lunar orbit.
        assert!((ratio - 1.0).abs() < 5e-3, "{ratio}");
        // An exact Keplerian orbit gives exactly 1.
        let (r, gm): (f64, f64) = (1e7, 4e14);
        let t = 2.0 * PI * (r * r * r / gm).sqrt();
        assert!(close(five_constants(r, t, gm / G * C * C / H), 1.0, 1e-12));
    }

    #[test]
    fn spheroids_turn_in_a_field() {
        // A sphere: L = 1/3 and 3K both ways.
        assert!(close(depolarisation(1.0), 1.0 / 3.0, 1e-15));
        let (long, short) = spheroid_polarisability(2.3, 1.0);
        assert!(close(long, 3.0 * clausius_mossotti(2.3), 1e-14) && close(short, long, 1e-14));
        // Osborn's table: a 2:1 prolate spheroid has L = 0.1736 along its long axis.
        assert!(close(depolarisation(2.0), 0.17356, 1e-4), "{}", depolarisation(2.0));
        // The series and the closed form agree where they meet.
        let edge = 1.0 / (1.0f64 - 0.09).sqrt();
        assert!(close(depolarisation(edge * (1.0 - 1e-9)), depolarisation(edge * (1.0 + 1e-9)), 1e-8));
        // A needle of conductor: L → 0, polarisability without bound.
        assert!(spheroid_polarisability(1e9, 50.0).0 > 300.0);
        // A polyethylene ball 1% out of round: Δα = 0.003 27 ε₀, so a horizontal 5.6×10⁴ V/m
        // at 45° turns it with 1.7×10⁻¹¹ N m.
        let (long, short) = spheroid_polarisability(2.3, 1.01);
        let d = long - short;
        assert!(close(d, 3.2733e-3, 1e-4), "{d}");
        let tau = alignment_torque(4.0 / 3.0 * PI * 4.447e-3f64.powi(3), d, 5.63e4, PI / 4.0);
        assert!(close(tau, 1.6921e-11, 1e-4), "{tau:e}");
    }

    #[test]
    fn the_image_force() {
        // 1 nC 5 cm from glass (ε_r = 6): 6.4×10⁻⁷ N, and a conductor's q²/(16πε₀d²) as ε_r → ∞.
        assert!(close(image_force(1e-9, 0.05, 6.0), 6.4197e-7, 1e-4), "{:e}", image_force(1e-9, 0.05, 6.0));
        assert!(close(image_force(1e-9, 0.05, 1e12), 1e-18 / (16.0 * PI * EPS_0 * 0.0025), 1e-9));
    }

    #[test]
    fn zhangs_constants_in_his_own_terms() {
        let f: Vec<f64> = ZHANG_KMQ.iter().map(|&(_, k)| vortex_constant(k)).collect();
        // The three values of f: 9.5×10⁹, 5.3×10⁻³⁵ and 4.3×10⁻³⁶ kg s/C. Paper E's is exactly
        // e/(4πε₀c³): G cancels.
        assert!(close(f[0], 9.534e9, 1e-3), "{:e}", f[0]);
        assert!(close(f[1], E_CHARGE / (4.0 * PI * EPS_0 * C * C * C), 1e-12));
        assert!(close(f[2], f[1] / (4.0 * PI), 1e-12));
        // The patent rig's 7.4×10⁻¹¹ T at the ball: 3.9×10⁻²¹ rad/s² with the calibrated value,
        // 7×10²³ rad/s² with Paper E's.
        assert!(close(vortex_spin(7.404e-11, f[0]), 3.883e-21, 1e-3));
        assert!(vortex_spin(7.404e-11, f[1]) > 1e23);
        // Earth's own field, 5×10⁻⁵ T and always on, is 7×10⁵ times the rig's: with Paper E's
        // constant everything on Earth would spin up at 5×10²⁹ rad/s².
        assert!(vortex_spin(5e-5, f[1]) > 4e29);
        // G = 4πc²/μ_g = 6.6849×10⁻¹¹: 1 390 ppm above the highest of CODATA's 16 measurements
        // and 34σ from it (BIPM-01), 70σ from CODATA 2022's 6.674 30(15).
        let g = 4.0 * PI * C * C / ZHANG_MU_G;
        assert!(close(g, 6.685e-11, 2e-4), "{g:e}");
        let (top, u) = crate::anomalies::CODATA_2022_G.iter().map(|&(_, v, u)| (v, u)).fold((0.0, 0.0), |a, b| if b.0 > a.0 { b } else { a });
        assert!((g * 1e11 / top - 1.0) * 1e6 > 1390.0 && (g * 1e11 - top) / u > 34.0);
        assert!((g * 1e11 - 6.674_30) / 0.000_15 > 70.0);
        // Paper J's closures with its own numbers: (e k_mq)² μ_g/(ħc) is 6.9×10⁻⁵, not
        // α = 7.3×10⁻³, and (4π k_mq μ_g)⁻¹ is 6.7×10⁻²⁰, not ε₀ = 8.9×10⁻¹².
        let k = ZHANG_KMQ[0].1;
        let alpha = (E_CHARGE * k).powi(2) * ZHANG_MU_G / (metrology::HBAR * C);
        assert!(close(alpha, 6.88e-5, 1e-2), "{alpha:e}");
        assert!(close(1.0 / (4.0 * PI * k * ZHANG_MU_G), 6.65e-20, 1e-2));
    }

    #[test]
    fn mass_independence_is_the_discriminator() {
        // Polyethylene (950) against tungsten (19 300) of the same shape.
        assert!(close(mass_independence(950.0, 19_300.0), 20.32, 1e-3));
        // The same force on both bodies: the acceleration falls by exactly that ratio.
        let (volume, force) = (1e-6, 1e-6);
        let light = acceleration(force, volume, 950.0);
        let dense = acceleration(force, volume, 19_300.0);
        assert!(close(light / dense, mass_independence(950.0, 19_300.0), 1e-12));
        // A gravitational field would give both the same acceleration, so the ratio would be 1.
        assert!(mass_independence(950.0, 950.0) == 1.0);
        // The polyethylene ball of the reported setup, if the claimed force were 1 mN on 1 g:
        assert!(close(acceleration(1e-3, 1e-6, 1000.0), 1.0, 1e-12));
    }
}
