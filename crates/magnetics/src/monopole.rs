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
//! Sources: the PDG review "Magnetic monopoles" (Particle Data Group 2025, review 94; the 2024
//! edition numbers the same equations 95.x); Milton,
//! [Rep. Prog. Phys. 69, 1637 (2006)](https://arxiv.org/abs/hep-ex/0602040); Castelnovo,
//! Moessner & Sondhi, [Nature 451, 42 (2008)](https://doi.org/10.1038/nature06433).

use core::f64::consts::PI;

use metrology::{ALPHA, C, DIRAC_FLUX, E_CHARGE, H, HBAR, M_E, MU_0};

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
/// `gBL`, the work of the force `gB` (derived from the field; PDG prints no such equation).
/// Galactic fields give about 10¹¹ GeV, the mass above which PDG's eq. 94.8 makes monopoles
/// non-relativistic.
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

/// The isotropic flux bound (m⁻² s⁻¹ sr⁻¹) set by the survival of the galactic magnetic field,
/// the Parker bound, for monopoles of mass `m_gev` (GeV): 10⁻¹⁵ cm⁻² s⁻¹ sr⁻¹, rising as
/// M/10¹⁷ GeV above 10¹⁷ GeV (PDG 2025, eq. 94.9).
pub fn parker_bound(m_gev: f64) -> f64 {
    1e-15 * (m_gev / 1e17).max(1.0) * 1e4
}

/// The extended Parker bound (m⁻² s⁻¹ sr⁻¹) from the survival of an early seed field:
/// `[M/10¹⁷ GeV + 3×10⁻⁶] × 10⁻¹⁶ cm⁻² s⁻¹ sr⁻¹` (PDG 2025, eq. 94.10).
pub fn parker_seed_bound(m_gev: f64) -> f64 {
    (m_gev / 1e17 + 3e-6) * 1e-16 * 1e4
}

/// Cross-section (m²) for a monopole moving at `beta` to catalyse nucleon decay (Callan–Rubakov),
/// `σ₀/β`, times `(β/β₀)^γ` below `low_beta = Some((β₀, γ))`, where the nucleus's own physics
/// suppresses or enhances it (IceCube 2014, eq. 4). Catalysis is a prediction for GUT monopoles
/// only, never observed; σ₀ is of the strong-interaction size, about 1 mb.
pub fn catalysis_cross_section(sigma0: f64, beta: f64, low_beta: Option<(f64, f64)>) -> f64 {
    let sigma = sigma0 / beta;
    match low_beta {
        Some((beta0, gamma)) if beta < beta0 => sigma * (beta / beta0).powf(gamma),
        _ => sigma,
    }
}

/// Mean distance (m) between catalysed decays, `1/(σn)`, for cross-section `sigma` (m²) and
/// target density `n` (m⁻³) (IceCube 2014, eq. 5).
pub fn catalysis_length(sigma: f64, n: f64) -> f64 {
    1.0 / (sigma * n)
}

/// Super-Kamiokande's bound (m⁻² s⁻¹ sr⁻¹) on GUT monopoles catalysing proton decay in the Sun,
/// `6.3×10⁻²⁴ (β/10⁻³)² / [(σ₀/1 mb)(f_π⁺/0.5)]` cm⁻² s⁻¹ sr⁻¹ for `sigma0_mb` (mb) and the
/// branching fraction `f_pi` to π⁺ (Ueno et al., Astropart. Phys. 36, 131 (2012), eq. 9). `None`
/// outside β < 10⁻³, where the Sun's gravitational capture, which sets the β², holds; it also
/// needs M ≲ 10¹⁷ GeV.
pub fn superk_catalysis_limit(beta: f64, sigma0_mb: f64, f_pi: f64) -> Option<f64> {
    (beta <= 1e-3).then(|| 6.3e-24 * (beta / 1e-3).powi(2) / (sigma0_mb * f_pi / 0.5) * 1e4)
}

/// Energy (J) of the bound state `N` (radial quantum number), total angular momentum `j`, of an
/// electron and a dyon of charge `z`e and `m` = `n`/2 for `n` Dirac units, the monopole far
/// heavier than the electron: `−Z²Ry/[N + ½ + √((j + ½)² − m²)]²` (Milton 2006, eq. 9.1, spin
/// omitted). `m = 0` is hydrogen; a Dirac pole (m = ½) binds hydrogen's ground state at 7.29 eV,
/// not 13.61 eV. `None` unless `j ≥ |m|`.
pub fn dyonic_level(z: f64, n_radial: u32, j: f64, m: f64) -> Option<f64> {
    let rydberg = ALPHA * ALPHA * M_E * C * C / 2.0;
    (j >= m.abs()).then(|| {
        let d = f64::from(n_radial) + 0.5 + ((j + 0.5).powi(2) - m * m).sqrt();
        -z * z * rydberg / (d * d)
    })
}

/// The binding parameter `γ̂ = (A/Z)·µ/µ_N` of a nucleus of charge `z`, mass number `a` and
/// magnetic moment `mu` (nuclear magnetons) in a monopole's field (Milton 2006, eq. 9.6, with
/// Table 1's convention, which equals Olaussen & Sollie's `S·g_M`).
pub fn gamma_hat(z: f64, a: f64, mu: f64) -> f64 {
    a / z * mu
}

/// Whether a nucleus binds to a monopole, by the regulator-independent criteria of
/// Bracci & Fiorentini (1984) and Olaussen & Sollie (1985), as Milton 2006 §9 collects them.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Binding {
    /// Bound, in the state of total angular momentum `j`.
    Binds { j: f64 },
    /// Not bound in any state.
    No,
    /// Not bound in the states the published criteria cover (the two lowest for spin above ½),
    /// so undecided: negative moments can bind in higher states, as ³He and ¹¹³Cd do.
    Undetermined,
}

/// Whether a nucleus of charge `z`, mass number `a`, spin `spin` and moment `mu` (nuclear
/// magnetons) binds to a monopole of `n` Dirac units (`m′ = n/2`, `l = |m′|Z`):
/// - a neutron (`z = 0`) binds in J = ½ if `|µ| > 3/(4|m′|)` (Milton 9.5);
/// - the lowest state J = l − S binds if `γ̂ > 1 + 1/(4l)`, for any spin (9.7; O–S 4.16a);
/// - spin ½: J = l + ½ binds if `|γ̂ − 1| > (J² + J − l²)/l` (9.8), and higher J bind less;
/// - spin > ½: J = l − S + 1 binds if Olaussen & Sollie's `λ±` (9.9) exceeds ¼.
///
/// Binding *energies* depend on the short-distance regulator and are not computed. The criteria
/// need `l ≥ S`; below that the result is `Undetermined`.
pub fn monopole_binding(z: f64, a: f64, spin: f64, mu: f64, n: i32) -> Binding {
    let m = f64::from(n).abs() / 2.0;
    if z == 0.0 {
        return if mu.abs() > 3.0 / (4.0 * m) { Binding::Binds { j: 0.5 } } else { Binding::No };
    }
    let l = m * z;
    if spin == 0.0 || mu == 0.0 {
        return Binding::No; // no moment, nothing to bind with: O–S λ₀ = −l
    }
    if l < spin {
        return Binding::Undetermined;
    }
    let g = gamma_hat(z, a, mu);
    if g > 1.0 + 1.0 / (4.0 * l) {
        return Binding::Binds { j: l - spin };
    }
    if spin == 0.5 {
        let j = l + 0.5;
        return if (g - 1.0).abs() > (j * j + j - l * l) / l { Binding::Binds { j } } else { Binding::No };
    }
    let (plus, minus) = os_lambdas(g, spin, l);
    if plus > 0.25 || minus > 0.25 { Binding::Binds { j: l - spin + 1.0 } } else { Binding::Undetermined }
}

/// Olaussen & Sollie's eigenvalues for J = l − S + 1 (their 4.16b; Milton 9.9), with
/// X = γ̂l/S. NaN where they are complex, which binds nothing.
fn os_lambdas(g: f64, spin: f64, l: f64) -> (f64, f64) {
    let x = g * l / spin;
    let root = ((l + 1.0).powi(2) + (2.0 * spin - 1.0 - l) * x + x * x / 4.0).sqrt();
    let mid = (spin - 0.5) * x - 2.0 * l - 1.0;
    (mid + root, mid - root)
}

/// Kazama–Yang–Goldhaber and Bloch corrections (K, B) for `n` Dirac units (Milton 2006 §8).
fn fast_corrections(n: i32) -> Option<(f64, f64)> {
    let k = if n.abs() == 1 { 0.406 } else { 0.346 };
    let b = match n.abs() {
        1 => 0.248,
        2 => 0.672,
        3 => 1.022,
        6 => 1.685,
        _ => return None,
    };
    Some((k, b))
}

/// Free-electron Fermi speed `(ħ/m_e)(3π²n_e)^(1/3)`, m/s.
fn fermi_speed(n_e: f64) -> f64 {
    HBAR / M_E * (3.0 * PI * PI * n_e).cbrt()
}

/// Energy loss (J/m) of a fast monopole of `n` Dirac units at speed `beta` in matter of electron
/// density `n_e` (m⁻³), mean excitation energy `i` (J) and density correction `delta`: Ahlen's
/// (1978) Bethe formula with `Ze → gβ`, `πn_e n²(ħc)²/(m_ec²)·[ln(2m_ec²β²γ²/I) + K/2 − δ/2 − ½ − B]`
/// (Patrizii & Spurio 2015, eq. 11; K and B from Milton 2006 §8). It does not fall with β like
/// a charge's: a Dirac pole ionises like a nucleus of Z = 68.5β. `None` outside 0.05 < β with
/// γ ≤ 100, or for an `n` without a published Bloch correction (1, 2, 3 and 6 have one).
pub fn monopole_dedx_fast(n: i32, beta: f64, n_e: f64, i: f64, delta: f64) -> Option<f64> {
    let gamma = 1.0 / (1.0 - beta * beta).sqrt();
    let (k, b) = fast_corrections(n)?;
    let mec2 = M_E * C * C;
    let hbar_c = HBAR * C;
    (beta > 0.05 && gamma <= 100.0).then(|| {
        let log = (2.0 * mec2 * (beta * gamma).powi(2) / i).ln();
        PI * n_e * f64::from(n * n) * hbar_c * hbar_c / mec2 * (log + k / 2.0 - delta / 2.0 - 0.5 - b)
    })
}

/// Energy loss (J/m) of a slow monopole of `n` Dirac units at `beta` in a degenerate electron gas
/// of density `n_e` (m⁻³): Ahlen & Kinoshita (1982),
/// `1.37 × (π/2) n_e n²(ħc)² β/(m_e c v_F) · [ln(2v_F/(αc)) − ½]`, the factor 1.37 from the
/// electron's spin (Patrizii & Spurio 2015, eq. 12). Good for metals, reasonable for Z ≥ 10;
/// it does not give scintillator light. `None` outside 10⁻³ ≤ β ≤ 10⁻², and where the logarithm
/// turns negative (v_F < 0.82αc, n_e below 1.28×10²⁹ m⁻³), where the electron gas is too
/// dilute for the model.
pub fn monopole_dedx_slow(n: i32, beta: f64, n_e: f64) -> Option<f64> {
    let v_f = fermi_speed(n_e);
    let log = (2.0 * v_f / (ALPHA * C)).ln() - 0.5;
    let hbar_c = HBAR * C;
    ((1e-3..=1e-2).contains(&beta) && log > 0.0)
        .then(|| 1.37 * PI / 2.0 * n_e * f64::from(n * n) * hbar_c * hbar_c * beta / (M_E * C * v_f) * log)
}

#[cfg(test)]
mod tests {
    use super::*;
    use metrology::{AMU, EV, GEV, MEV, MU_B, N_A, PARSEC};

    fn close(a: f64, b: f64, rel: f64) -> bool {
        (a / b - 1.0).abs() < rel
    }

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
        // 1 G over 10 cm: 205 keV (de Lacheze-Murel 2022). 3 µG over 0.526 kpc gives 10¹¹ GeV =
        // 16.02 J, PDG eq. 94.8's boundary between relativistic and slow galactic monopoles.
        assert!(close(monopole_energy_gain(1, 1e-4, 0.1) / (1e3 * EV), 205.4, 1e-3));
        assert!(close(monopole_energy_gain(1, 3e-10, 0.5257 * 1e3 * PARSEC), 1e11 * GEV, 2e-3));
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

    #[test]
    fn flux_bounds() {
        // PDG 2025 eqs 94.9 and 94.10, in cm⁻² s⁻¹ sr⁻¹ (×10⁻⁴ from SI).
        for (m, parker, seed) in [(1e15, 1e-15, 1.0003e-18), (1e17, 1e-15, 1.000_003e-16), (1e19, 1e-13, 1.000_000_03e-14)] {
            assert!(close(parker_bound(m) * 1e-4, parker, 1e-12) && close(parker_seed_bound(m) * 1e-4, seed, 1e-12), "{m}");
        }
        assert!(close(parker_seed_bound(0.0) * 1e-4, 3e-22, 1e-12));
        // Cabrera's event implies 6.1×10⁵ times the Parker bound.
        assert!(close(flux_from_events(1.0, 151.0 * 86_400.0, 20e-4) / parker_bound(1e16), 6.1e5, 2e-3));
    }

    #[test]
    fn catalysis() {
        // IceCube 2014: σ₀ = 1.7×10⁻²⁵ cm² at β = 10⁻³ is 1.7×10⁻²² cm²; 1.7×10⁻²³ cm² in ice's
        // hydrogen (6.131×10²² cm⁻³ at 0.917 g/cm³) is a 0.959 cm mean path (IceCube: "1 cm").
        assert!(close(catalysis_cross_section(1.7e-29, 1e-3, None), 1.7e-26, 1e-12));
        assert!(close(catalysis_cross_section(1.7e-29, 1e-3, Some((1e-4, 2.0))), 1.7e-26, 1e-12));
        assert!(close(catalysis_cross_section(1.0, 5e-4, Some((1e-3, 1.0))), 1.0 / 5e-4 * 0.5, 1e-12));
        let n_h = 2.0 * 917.0 / 0.018_015 * N_A;
        assert!(close(n_h, 6.131e28, 1e-3) && close(catalysis_length(1.7e-27, n_h), 0.959e-2, 1e-3));
        // Super-K, eq. 9: 6.3×10⁻²⁴ cm⁻² s⁻¹ sr⁻¹ at β = 10⁻³, falling as β², more than eight orders
        // below Parker; σ₀ = 10 mb tightens it tenfold; it does not apply at β ≥ 10⁻³.
        let f = |b: f64, s: f64| superk_catalysis_limit(b, s, 0.5).unwrap() * 1e-4;
        assert!(close(f(1e-3, 1.0), 6.3e-24, 1e-12) && close(f(1e-4, 1.0), 6.3e-26, 1e-12) && close(f(1e-3, 10.0), 6.3e-25, 1e-12));
        assert!((1e8..1e9).contains(&(parker_bound(1e17) / superk_catalysis_limit(1e-3, 1.0, 0.5).unwrap())));
        assert!(superk_catalysis_limit(2e-3, 1.0, 0.5).is_none());
    }

    #[test]
    fn hydrogen_with_a_dirac_pole() {
        // Milton eq. 9.1: m = 0 is hydrogen's −13.6057 eV; a Dirac pole (m = ½) makes it −7.29127 eV.
        let ev = |e: Option<f64>| e.unwrap() / EV;
        assert!(close(ev(dyonic_level(1.0, 0, 0.0, 0.0)), -13.605_693, 1e-7));
        assert!(close(ev(dyonic_level(1.0, 0, 0.5, 0.5)), -7.291_27, 1e-5));
        assert!(close(dyonic_level(1.0, 0, 0.5, 0.5).unwrap() / dyonic_level(1.0, 0, 0.0, 0.0).unwrap(), 0.535_898, 1e-5));
        assert!(dyonic_level(1.0, 0, 0.0, 0.5).is_none());
    }

    #[test]
    fn which_nuclei_bind() {
        // Milton's Table 1 γ̂ (LiveChart moments).
        for (z, a, mu, g) in [(1.0, 1.0, 2.792_847, 2.79), (2.0, 3.0, -2.127_625, -3.20), (4.0, 9.0, -1.177_43, -2.66), (13.0, 27.0, 3.6407, 7.56), (48.0, 113.0, -0.622_26, -1.46), (1.0, 2.0, 0.857_438, 1.71)] {
            assert!(close(gamma_hat(z, a, mu), g, 6e-3), "{z} {a}");
        }
        // Olaussen & Sollie's table, recomputed: ²⁷Al λ₀, λ± = 42.659, 35.513, 15.141 (printed
        // 42.8, 35.6, 15.2); ²³Na 20.001, 12.417, −2.416 (20.0, 12.4, −2.4).
        for (z, a, s, mu, want) in [(13.0, 27.0, 2.5, 3.6414, [42.659, 35.513, 15.141]), (11.0, 23.0, 1.5, 2.2175, [20.001, 12.417, -2.416])] {
            let (g, l) = (gamma_hat(z, a, mu), z / 2.0);
            let (plus, minus) = os_lambdas(g, s, l);
            assert!(close(g * l - l, want[0], 1e-4) && close(plus, want[1], 1e-4) && close(minus, want[2], 1e-3), "{z}");
        }
        // Milton's closed forms are where λ± = ¼: 9.12 (S = 3/2, l = 2) 1.556683 and 8.943317;
        // 9.13 (S = 5/2, l = 6.5) 1.373638 and 4.216105; 9.10 (S = 1, l = 2) 2.183824.
        for (g, s, l, plus) in [(1.556_683, 1.5, 2.0, true), (8.943_317, 1.5, 2.0, false), (1.373_638, 2.5, 6.5, true), (4.216_105, 2.5, 6.5, false), (2.183_824, 1.0, 2.0, true)] {
            let (p, m) = os_lambdas(g, s, l);
            assert!(((if plus { p } else { m }) - 0.25).abs() < 5e-5, "{g}");
        }
        // The table's verdicts: ¹H in J = 0, ³He in 3/2, ¹¹³Cd in 49/2, ²⁷Al in 4, ²H (|m′| = 1) in
        // 0, the neutron in ½; ⁵⁷Fe and ¹⁶O do not bind; ⁹Be and ¹⁷O are undecided by the
        // published criteria.
        use Binding::*;
        let cases = [
            ((1.0, 1.0, 0.5, 2.792_847, 1), Binds { j: 0.0 }),
            ((2.0, 3.0, 0.5, -2.127_625, 1), Binds { j: 1.5 }),
            ((48.0, 113.0, 0.5, -0.622_26, 1), Binds { j: 24.5 }),
            ((13.0, 27.0, 2.5, 3.6407, 1), Binds { j: 4.0 }),
            ((1.0, 2.0, 1.0, 0.857_438, 2), Binds { j: 0.0 }),
            ((0.0, 1.0, 0.5, -1.913_043, 1), Binds { j: 0.5 }),
            ((26.0, 57.0, 0.5, 0.0902, 1), No),
            ((8.0, 16.0, 0.0, 0.0, 1), No),
            ((4.0, 9.0, 1.5, -1.177_43, 1), Undetermined),
            ((8.0, 17.0, 2.5, -1.893_543, 1), Undetermined),
        ];
        for ((z, a, s, mu, n), want) in cases {
            assert_eq!(monopole_binding(z, a, s, mu, n), want, "{z} {a}");
        }
    }

    #[test]
    fn stopping_power() {
        // Silicon at 2.33 g/cm³: n_e = 6.99444×10²⁹ m⁻³, v_F = 3.179×10⁶ m/s; at the mantle's
        // 4.3 g/cm³, 3.900×10⁸ cm/s as Patrizii & Spurio print.
        let n_e = |rho: f64| rho * N_A * 14.0 / 0.028_085_5;
        assert!(close(n_e(2330.0), 6.994_44e29, 1e-5) && close(fermi_speed(n_e(4300.0)) * 100.0, 3.9e8, 1e-3));
        // Slow, β = 10⁻³: 9.8239×10⁻¹⁰ J/m (61.3 MeV/cm); without the spin factor 44.76 MeV/cm,
        // Milton's eq. 8.11 "45 GeV/cm × β" to 0.6%.
        let slow = monopole_dedx_slow(1, 1e-3, n_e(2330.0)).unwrap();
        assert!(close(slow, 9.8239e-10, 1e-4) && close(slow / 1.37 / MEV * 1e-2, 45.0, 6e-3));
        assert!(monopole_dedx_slow(1, 0.1, n_e(2330.0)).is_none() && monopole_dedx_slow(1, 5e-3, 1e29).is_none());
        assert!(close(monopole_dedx_slow(2, 5e-3, n_e(2330.0)).unwrap(), 4.0 * 5.0 * slow, 1e-12));
        // Fast, β = 0.5, I = 173 eV, δ = 0: 11.788 GeV/cm; its prefactor is the Bethe constant
        // K = 0.307075 MeV cm²/mol times (Z/A)ρ(n/2α)², a nucleus of charge g_D β/e = 68.5β.
        let fast = monopole_dedx_fast(1, 0.5, n_e(2330.0), 173.0 * EV, 0.0).unwrap();
        assert!(close(fast / GEV * 1e-2, 11.788, 1e-4));
        let mec2 = M_E * C * C;
        let bracket = (2.0 * mec2 * 0.25 / 0.75 / (173.0 * EV)).ln() + 0.406 / 2.0 - 0.5 - 0.248;
        let bethe = 0.307_075 * 14.0 / 28.0855 * 2.33 / (2.0 * ALPHA).powi(2) * MEV * 1e2;
        assert!(close(fast / bracket, bethe, 1e-5));
        assert!(monopole_dedx_fast(4, 0.5, 1e29, 173.0 * EV, 0.0).is_none() && monopole_dedx_fast(1, 0.01, 1e29, 173.0 * EV, 0.0).is_none());
    }
}
