//! Aneutronic fuels and their power balance: Q-values, the neutrons D–³He still makes, and why
//! p–¹¹B fights bremsstrahlung and synchrotron radiation.
//!
//! ```
//! use plasmoid::aneutronic::pb11_fusion_over_brems;
//! use plasmoid::fusion::Pb11Fit;
//!
//! // With equal electron and ion temperatures, p–¹¹B fusion never outruns bremsstrahlung: at
//! // 300 keV the best boron fraction gives P_fus/P_br = 0.46.
//! let best = (20..500).map(|k| pb11_fusion_over_brems(300.0, 300.0, f64::from(k) / 1000.0, Pb11Fit::TentoriBelloni).unwrap()).fold(0.0, f64::max);
//! assert!((best - 0.457).abs() < 0.002);
//! ```

use crate::fusion::{KEV, Pb11Fit, Reaction, bosch_hale, pb11_reactivity};

/// Q-values (MeV) from NIST atomic masses.
pub const Q_PB11: f64 = 8.6821;
/// D + ³He → ⁴He + p.
pub const Q_DHE3: f64 = 18.3531;
/// ³He + ³He → ⁴He + 2p.
pub const Q_HE3HE3: f64 = 12.8596;
/// p + ⁶Li → ⁴He + ³He.
pub const Q_PLI6: f64 = 4.0197;
/// D + T → ⁴He + n.
pub const Q_DT: f64 = 17.5893;
/// D + D → ³He + n, and its neutron's energy.
pub const Q_DDN: f64 = 3.2689;
/// D + D → T + p.
pub const Q_DDP: f64 = 4.0327;
const E_DD_NEUTRON: f64 = 2.45;
const E_DT_NEUTRON: f64 = 14.07;

const MEC2_KEV: f64 = 510.998_95;
/// `W₀ = α r_e² m_e c³` (W m³), the bremsstrahlung scale.
const W0: f64 = 1.4219e-36;

/// Fraction of D–³He fusion power carried by neutrons from the D–D side reactions, at `t` (keV)
/// with deuterium fraction `x_d` = n_D/(n_D + n_He3); `burn_tritium` adds the bred tritium burning
/// at once. `None` outside the Bosch–Hale fits.
pub fn dhe3_neutron_fraction(t: f64, x_d: f64, burn_tritium: bool) -> Option<f64> {
    let (n_d, n_h) = (x_d, 1.0 - x_d);
    let (dhe, ddn, ddp) = (bosch_hale(Reaction::DHe3, t)?, bosch_hale(Reaction::DDn, t)?, bosch_hale(Reaction::DDp, t)?);
    let mut total = n_d * n_h * dhe * Q_DHE3 + 0.5 * n_d * n_d * (ddn * Q_DDN + ddp * Q_DDP);
    let mut neutrons = 0.5 * n_d * n_d * ddn * E_DD_NEUTRON;
    if burn_tritium {
        total += 0.5 * n_d * n_d * ddp * Q_DT;
        neutrons += 0.5 * n_d * n_d * ddp * E_DT_NEUTRON;
    }
    Some(neutrons / total)
}

/// Bremsstrahlung power density (W/m³) of electrons of density `ne` (m⁻³) at `te` (keV) with
/// effective charge `z_eff`, including the relativistic and electron–electron terms (Putvinski
/// et al. 2019, Eq. 7, a fit to QED results within 1%). At low temperature it is
/// `5.355×10⁻³⁷ Z_eff n_e² √T`.
pub fn bremsstrahlung(ne: f64, te: f64, z_eff: f64) -> f64 {
    let t = te / MEC2_KEV;
    8.511 * W0 * ne * ne * t.sqrt() * (z_eff * (1.0 + 1.78 * t.powf(1.34)) + 2.12 * t * (1.0 + 1.1 * t + t * t - 1.25 * t.powf(2.5)))
}

/// Synchrotron power density (W/m³) with no wall reflection, `6.21×10⁻¹⁷ B² n_e T_e` (NRL), for
/// field `b` (T), density `ne` (m⁻³) and `te` (keV).
pub fn synchrotron(b: f64, ne: f64, te: f64) -> f64 {
    6.21e-17 * b * b * ne * te
}

/// p–¹¹B fusion power over bremsstrahlung (independent of density) at ion temperature `ti`,
/// electron temperature `te` (keV) and boron fraction `f_b` = n_B/(n_p + n_B).
pub fn pb11_fusion_over_brems(ti: f64, te: f64, f_b: f64, fit: Pb11Fit) -> Option<f64> {
    let (n_p, n_b) = (1.0 - f_b, f_b);
    let ne = n_p + 5.0 * n_b;
    let z_eff = (n_p + 25.0 * n_b) / ne;
    Some(n_p * n_b * pb11_reactivity(fit, ti)? * Q_PB11 * 1e3 * KEV / bremsstrahlung(ne, te, z_eff))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, rel: f64) -> bool {
        (a / b - 1.0).abs() < rel
    }

    #[test]
    fn deuterium_helium3_is_not_neutron_free() {
        assert!(close(dhe3_neutron_fraction(50.0, 0.5, false).unwrap(), 0.0131, 0.01));
        assert!(close(dhe3_neutron_fraction(50.0, 0.5, true).unwrap(), 0.0726, 0.01));
        assert!(close(dhe3_neutron_fraction(100.0, 0.25, false).unwrap(), 0.0034, 0.02));
    }

    #[test]
    fn radiation_losses() {
        // Low-temperature limit: NRL's 1.69×10⁻³² W cm³ eV^(−½) is 5.344×10⁻³⁷ W m³ keV^(−½).
        let low = bremsstrahlung(1.0, 1e-3, 1.0) / 1e-3f64.sqrt();
        assert!(close(low, 5.355e-37, 2e-3));
        // Synchrotron without reflection, at β = 1, 0.5 and 0.1 for Ti = 300, Te = 150 keV, f_B = 0.15:
        // 41, 82 and 410 times bremsstrahlung.
        let (n_p, n_b) = (0.85, 0.15);
        let ne_rel = n_p + 5.0 * n_b;
        let z_eff = (n_p + 25.0 * n_b) / ne_rel;
        let ne = 1e20;
        for (beta, expect) in [(1.0, 41.03), (0.1, 410.3)] {
            let b2 = 2.0 * 1.256_637_062_12e-6 * (ne / ne_rel * 300.0 + ne * 150.0) * KEV / beta;
            assert!(close(synchrotron(b2.sqrt(), ne, 150.0) / bremsstrahlung(ne, 150.0, z_eff), expect, 0.005));
        }
    }

    #[test]
    fn proton_boron_balance() {
        let ratio = |ti: f64, te: f64, f: f64| pb11_fusion_over_brems(ti, te, f, Pb11Fit::TentoriBelloni).unwrap();
        let (best, at) = (20..500).map(|k| (ratio(300.0, 300.0, f64::from(k) / 1000.0), f64::from(k) / 1000.0)).fold((0.0, 0.0), |a, b| if b.0 > a.0 { b } else { a });
        assert!(close(best, 0.4572, 2e-3) && (at - 0.105).abs() < 0.003);
        // Never above 1 with Te = Ti, between 70 and 500 keV; above 1 once the electrons run cold.
        assert!((70..=500).all(|t| ratio(f64::from(t), f64::from(t), 0.15) < 1.0));
        assert!(ratio(300.0, 100.0, 0.15) > 1.0);
    }
}
