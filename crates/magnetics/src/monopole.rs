//! Magnetic monopoles: Dirac's quantised magnetic charge, its field, vector potential and field
//! energy, the energy it gains crossing a field, how a search turns events into a flux, and the
//! emergent and synthetic monopoles that condensed matter and cold atoms make.
//!
//! Dirac's condition `eg = nh/μ₀` makes the smallest magnetic charge `g_D = h/(μ₀e) =
//! 3.291×10⁻⁹ A m` (68.5 e in Gaussian units), and its flux `μ₀g_D = h/e` two superconducting flux
//! quanta. No free monopole has been found: Cabrera's single 1982 event was never repeated, and
//! the searches since bound the flux far below it. Spin ice, skyrmion lattices and Bose
//! condensates make real monopole *textures*, but their charge lives in the magnetisation, a
//! Berry phase or a superfluid velocity, never in B itself, where `∇·B = 0` holds for any current.
//!
//! ```
//! use magnetics::monopole::{dirac_charge, dirac_field, monopole_energy_gain};
//!
//! assert!((dirac_charge(1) / 3.291e-9 - 1.0).abs() < 2e-4);
//! // At one femtometre, a nucleon's size, a Dirac pole's field is 3.3×10¹⁴ T.
//! assert!((dirac_field(1e-15, 1) / 3.291e14 - 1.0).abs() < 2e-4);
//! // Crossing 1 G over 10 cm it gains 205 keV, which is why searches use magnets to accelerate them.
//! let kev = monopole_energy_gain(1, 1e-4, 0.1) / (1e3 * metrology::EV);
//! assert!((kev - 205.4).abs() < 0.1);
//! ```
//!
//! Sources: the PDG review "Magnetic monopoles" (Particle Data Group 2024, §95); Milton,
//! [Rep. Prog. Phys. 69, 1637 (2006)](https://arxiv.org/abs/hep-ex/0602040); Castelnovo,
//! Moessner & Sondhi, [Nature 451, 42 (2008)](https://doi.org/10.1038/nature06433).

use core::f64::consts::PI;

use metrology::{DIRAC_FLUX, E_CHARGE, H, HBAR, MU_0};

/// Magnetic charge `n·h/(μ₀e)`, A m, of `n` Dirac units: 3.291×10⁻⁹ A m for `n = 1`.
pub fn dirac_charge(n: i32) -> f64 {
    f64::from(n) * DIRAC_FLUX / MU_0
}

/// Field (T) at distance `r` (m) from a pole of `n` Dirac units: `μ₀g/(4πr²) = nħ/(2er²)`.
/// It depends on h and e alone, and its total flux is `n·h/e`.
pub fn dirac_field(r: f64, n: i32) -> f64 {
    f64::from(n) * HBAR / (2.0 * E_CHARGE * r * r)
}

/// The distance (m) at which a pole of `n` Dirac units makes the field `b` (T): `√(nħ/(2eB))`.
pub fn dirac_radius_for_field(b: f64, n: i32) -> f64 {
    (f64::from(n) * HBAR / (2.0 * E_CHARGE * b)).sqrt()
}

/// Vector potential (T m) at `r` (m, from the pole) of a pole of `n` Dirac units whose Dirac
/// string runs from the pole along `−n_hat` (a unit vector):
/// `A = (μ₀g/4π)(n̂ × r)/(r(r + r·n̂))`. Its curl is the radial field [`dirac_field`] everywhere
/// off the string, where it is singular. The opposite order, `r × n̂`, gives the field of the
/// opposite charge.
pub fn dirac_vector_potential(r: [f64; 3], n_hat: [f64; 3], n: i32) -> [f64; 3] {
    let len = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).sqrt();
    let along = r[0] * n_hat[0] + r[1] * n_hat[1] + r[2] * n_hat[2];
    let k = MU_0 * dirac_charge(n) / (4.0 * PI * len * (len + along));
    [
        k * (n_hat[1] * r[2] - n_hat[2] * r[1]),
        k * (n_hat[2] * r[0] - n_hat[0] * r[2]),
        k * (n_hat[0] * r[1] - n_hat[1] * r[0]),
    ]
}

/// Field energy (J) outside a radius `r0` (m) of a pole of `n` Dirac units:
/// `μ₀g²/(8πr₀) = n²ħc/(8αr₀)`, 3380 MeV·fm/r₀. It grows without bound as the core shrinks, so a
/// monopole's mass is set by its core (a GUT core gives about 17 M_X, PDG's 10¹⁷–10¹⁸ GeV).
pub fn monopole_field_energy(r0: f64, n: i32) -> f64 {
    let g = dirac_charge(n);
    MU_0 * g * g / (8.0 * PI * r0)
}

/// Energy (J) a pole of `n` Dirac units gains moving a distance `l` (m) along a field `b` (T):
/// `gBL`. Per metre it is the force `gB` on the pole.
pub fn monopole_energy_gain(n: i32, b: f64, l: f64) -> f64 {
    dirac_charge(n) * b * l
}

/// The isotropic flux (m⁻² s⁻¹ sr⁻¹) implied by `events` seen in time `t` (s) through a flat
/// loop of area `area` (m²) open to both hemispheres: `n/(t·A·2π)`.
pub fn flux_from_events(events: f64, t: f64, area: f64) -> f64 {
    events / (t * area * 2.0 * PI)
}

/// Emergent magnetic charge (A m) of a monopole in spin ice, for spins of moment `moment` (J/T)
/// on a pyrochlore lattice of cubic constant `a` (m): `2μ/a_d` with the diamond-lattice spacing
/// `a_d = √(3/2)·a` (Castelnovo, Moessner & Sondhi 2008). It is 1/8000 of `g_D` and sits in the
/// magnetisation H, not in B.
pub fn spin_ice_charge(moment: f64, a: f64) -> f64 {
    2.0 * moment / (1.5f64.sqrt() * a)
}

/// Circulation (m²/s) of the Dirac string of a synthetic monopole in a spinor Bose condensate of
/// atoms of mass `m_atom` (kg): `2h/M` (Ray et al., Nature 505, 657 (2014)). The "field" is the
/// superfluid's vorticity.
pub fn synthetic_monopole_circulation(m_atom: f64) -> f64 {
    2.0 * H / m_atom
}

/// The emergent field (T) that conduction electrons feel in a skyrmion lattice with one skyrmion
/// per `area` (m²): one flux `h/e` per skyrmion (Neubauer et al., PRL 102, 186602 (2009)), a Berry
/// phase rather than a field in B.
pub fn skyrmion_emergent_field(area: f64) -> f64 {
    DIRAC_FLUX / area
}

#[cfg(test)]
mod tests {
    use super::*;
    use metrology::{ALPHA, AMU, C, EV, M_E, MU_B, PARSEC};

    fn close(a: f64, b: f64, rel: f64) -> bool {
        (a / b - 1.0).abs() < rel
    }

    const MEV: f64 = 1e6 * EV;
    const FM: f64 = 1e-15;

    #[test]
    fn the_dirac_charge_and_its_flux() {
        let g = dirac_charge(1);
        assert!(close(g, 3.291e-9, 2e-4) && close(dirac_charge(-2), -2.0 * g, 1e-15));
        // In Gaussian units g_D/e = 1/(2α) = 68.52, and the Dirac condition is eg = h/μ₀.
        assert!(close(g / (E_CHARGE * C), 1.0 / (2.0 * ALPHA), 1e-9) && close(g / (E_CHARGE * C), 68.52, 1e-4));
        assert!(close(E_CHARGE * g, H / MU_0, 1e-15));
        // Its field μ₀g/(4πr²) carries the flux h/e through any sphere.
        for r in [1e-15, 5.291_772_109e-11, 1.0] {
            assert!(close(4.0 * PI * r * r * dirac_field(r, 1), DIRAC_FLUX, 1e-12));
            assert!(close(dirac_field(r, 1), MU_0 * g / (4.0 * PI * r * r), 1e-12));
        }
        // 3.291×10¹⁴ T at 1 fm, 1.1752×10⁵ T at the Bohr radius.
        assert!(close(dirac_field(1e-15, 1), 3.291e14, 2e-4));
        assert!(close(dirac_field(5.291_772_109e-11, 1), 1.1752e5, 2e-4));
    }

    #[test]
    fn the_radius_where_a_pole_makes_a_given_field() {
        // 10¹⁵ T at 0.5737 fm; Ivoilov & Urutskoev's 280 Oe and 460 Oe at 108.5 and 84.65 nm.
        assert!(close(dirac_radius_for_field(1e15, 1) / FM, 0.5737, 1e-3));
        assert!(close(dirac_radius_for_field(280e-4, 1) / 1e-9, 108.5, 2e-3));
        assert!(close(dirac_radius_for_field(460e-4, 1) / 1e-9, 84.65, 2e-3));
        let r = dirac_radius_for_field(7.0, 3);
        assert!(close(dirac_field(r, 3), 7.0, 1e-14));
    }

    #[test]
    fn the_vector_potential_curls_to_the_radial_field() {
        let curl = |p: [f64; 3], n_hat: [f64; 3]| {
            let d = 1e-6;
            let partial = |i: usize, j: usize| {
                let (mut a, mut b) = (p, p);
                a[j] += d;
                b[j] -= d;
                (dirac_vector_potential(a, n_hat, 1)[i] - dirac_vector_potential(b, n_hat, 1)[i]) / (2.0 * d)
            };
            [partial(2, 1) - partial(1, 2), partial(0, 2) - partial(2, 0), partial(1, 0) - partial(0, 1)]
        };
        let n_hats = [[0.0, 0.0, 1.0], [0.6, 0.0, 0.8], [0.0, -1.0, 0.0]];
        let points: [[f64; 3]; 3] = [[0.3, -0.4, 0.5], [-0.2, 0.7, 0.1], [1.1, 0.2, -0.3]];
        for (p, n_hat) in points.into_iter().zip(n_hats) {
            let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
            let b = curl(p, n_hat);
            for i in 0..3 {
                // B = μ₀g r̂/(4πr²), the field of a positive pole whatever the string's direction.
                assert!(close(b[i], dirac_field(r, 1) * p[i] / r, 1e-5), "{p:?} {i}");
            }
        }
        // At (0.3, −0.4, 0.5) m the z-component is 4.654×10⁻¹⁶ T; it is divergence-free there.
        assert!(close(curl([0.3, -0.4, 0.5], [0.0, 0.0, 1.0])[2], 4.654e-16, 1e-3));
        // The string: the potential diverges on the half-line along −n̂.
        assert!(!dirac_vector_potential([0.0, 0.0, -2.0], [0.0, 0.0, 1.0], 1)[0].is_finite());
    }

    #[test]
    fn field_energy_and_the_monopole_mass() {
        // μ₀g²/(8πr₀) = ħc/(8αr₀) = 3380.1 MeV fm / r₀.
        let k = monopole_field_energy(FM, 1);
        assert!(close(k, HBAR * C / (8.0 * ALPHA * FM), 1e-9) && close(k / MEV, 3380.1, 1e-4));
        // It equals m_p c² at r₀ = 3.6025 fm, and 1 J at 5.416×10⁻²⁵ m.
        assert!(close(monopole_field_energy(1.0, 1) / (938.272_088_16 * MEV) / FM, 3.6025, 1e-3));
        assert!(close(monopole_field_energy(1.0, 1), 5.416e-25, 1e-3));
        // A core of the GUT Compton length ħc/M_X gives M_X/(8α) = 17.13 M_X.
        assert!(close(1.0 / (8.0 * ALPHA), 17.129, 1e-4));
        assert!(close(monopole_field_energy(FM, 2), 4.0 * k, 1e-14));
    }

    #[test]
    fn energy_gained_in_a_field() {
        // 1 G over 10 cm: 205 keV (de Lacheze-Murel 2022). PDG: 3 µG over 0.526 kpc gives
        // 10¹¹ GeV = 16.02 J, the scale above which galactic monopoles are relativistic.
        assert!(close(monopole_energy_gain(1, 1e-4, 0.1) / (1e3 * EV), 205.4, 1e-3));
        assert!(close(monopole_energy_gain(1, 3e-10, 0.5257 * 1e3 * PARSEC), 1e11 * 1e9 * EV, 2e-3));
        // The force in 120 Oe (12 mT): 3.95×10⁻¹¹ N.
        assert!(close(monopole_energy_gain(1, 120e-4, 1.0), 3.949e-11, 1e-3));
    }

    #[test]
    fn cabrera_1982() {
        // One event in 151 days through a 20 cm² loop: 6.10×10⁻¹⁰ cm⁻² s⁻¹ sr⁻¹ (PRL 48, 1378).
        let f = flux_from_events(1.0, 151.0 * 86_400.0, 20e-4);
        assert!(close(f * 1e-4, 6.10e-10, 2e-3));
    }

    #[test]
    fn emergent_and_synthetic_monopoles() {
        // Dysprosium titanate: 10 µ_B spins, a = 3.54 Å (a_d = 4.3356 Å): 4.6 µ_B/Å = 4.278×10⁻¹³ A m
        // (Bramwell et al. 2009 measured 5 µ_B/Å).
        let q = spin_ice_charge(10.0 * MU_B, 3.54e-10);
        assert!(close(q, 4.278e-13, 1e-3) && close(q / (MU_B / 1e-10), 4.6, 5e-3));
        // Castelnovo's eq. 4: q_m/g_D = 20 r_e/a_d; g_D/q_m is about 8000.
        let r_e = ALPHA * HBAR / (M_E * C);
        assert!(close(q / dirac_charge(1), 20.0 * r_e / (1.5f64.sqrt() * 3.54e-10), 1e-6));
        assert!(close(dirac_charge(1) / q, 7692.0, 1e-3));
        // ⁸⁷Rb (86.909 180 529 u): 4πħ/M = 9.183×10⁻⁹ m²/s.
        assert!(close(synthetic_monopole_circulation(86.909_180_529 * AMU), 9.183e-9, 1e-3));
        assert!(close(synthetic_monopole_circulation(1.0), 4.0 * PI * HBAR, 1e-15));
        // One flux h/e per skyrmion: 2.5 T for one skyrmion per 1.654×10⁻¹⁵ m².
        assert!(close(skyrmion_emergent_field(1.6543e-15), 2.5, 1e-4));
    }
}
