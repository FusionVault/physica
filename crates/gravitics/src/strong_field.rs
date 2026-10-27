//! The Schwinger limit and strong-field QED: the field at which the vacuum breaks down into
//! electron–positron pairs, how close lasers and accelerators get, and the vacuum's nonlinear
//! optics below it.
//!
//! The critical field `E_S = m_e²c³/(eħ) = 1.32×10¹⁸ V/m` has never been reached in the lab: the
//! record laser (CoReLS, 1.1×10²³ W/cm²) is 1,454 times short in field. The strong-field regime
//! has been reached only in the electron's rest frame (SLAC E144 at χ ≈ 0.3, crystal channelling
//! at CERN up to χ = 2.4). Below E_S the vacuum is a weakly nonlinear medium (Euler–Heisenberg),
//! birefringent in a magnetic field; PVLAS's best measurement is still 6.8 times the predicted signal.
//!
//! ```
//! use gravitics::strong_field::{E_SCHWINGER, schwinger_rate};
//!
//! assert!((E_SCHWINGER / 1.32329e18 - 1.0).abs() < 1e-5);
//! // At a tenth of the critical field the vacuum makes 2.5×10⁴⁰ pairs per m³ per second; at a
//! // hundredth, 4×10⁻⁸⁵: the exponential suppression e^(−πE_S/E).
//! assert!((schwinger_rate(0.1 * E_SCHWINGER) / 2.4688e40 - 1.0).abs() < 1e-3);
//! ```

use core::f64::consts::PI;

use metrology::{ALPHA, C, E_CHARGE, EPS_0, GEV, HBAR, KEV, M_E, MEV, MU_0};

/// Electron rest energy, keV.
const MEC2_KEV: f64 = M_E * C * C / KEV;

/// The Schwinger critical field `m_e²c³/(eħ)`, V/m.
pub const E_SCHWINGER: f64 = M_E * M_E * C * C * C / (E_CHARGE * HBAR);
/// The critical magnetic field `E_S/c`, T.
pub const B_SCHWINGER: f64 = E_SCHWINGER / C;
/// Nuclear charge at which the 1s level of an extended nucleus dives into the negative-energy
/// continuum (published value; reached only transiently in heavy-ion collisions, not observed).
pub const Z_CRITICAL: f64 = 173.0;

/// How a laser's intensity relates to its peak field `E₀`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Intensity {
    /// Linearly polarised, cycle-averaged: `I = ε₀cE₀²/2`.
    Linear,
    /// `I = ε₀cE₀²` (circular polarisation, or the Gaussian-units `E²/4π` convention).
    Peak,
}

impl Intensity {
    fn factor(self) -> f64 {
        match self {
            Intensity::Linear => 0.5,
            Intensity::Peak => 1.0,
        }
    }
}

/// The intensity (W/m²) whose peak field is the Schwinger field.
pub fn critical_intensity(conv: Intensity) -> f64 {
    conv.factor() * EPS_0 * C * E_SCHWINGER * E_SCHWINGER
}

/// Peak field (V/m) of a laser of `intensity` (W/m²).
pub fn peak_field(intensity: f64, conv: Intensity) -> f64 {
    (intensity / (conv.factor() * EPS_0 * C)).sqrt()
}

/// Pair-production rate (m⁻³ s⁻¹) in a constant field `e` (V/m): the leading term of Schwinger's
/// series, `(αε₀E²/π²ħ) e^(−πE_S/E)`. This is the mean rate of pair creation.
pub fn schwinger_rate(e: f64) -> f64 {
    ALPHA * EPS_0 * e * e / (PI * PI * HBAR) * (-PI * E_SCHWINGER / e).exp()
}

/// The full series `(αε₀E²/π²ħ) Σ n⁻² e^(−nπE_S/E)`: the vacuum's decay (persistence) rate,
/// which exceeds the pair rate only near E_S.
pub fn vacuum_decay_rate(e: f64) -> f64 {
    let x = PI * E_SCHWINGER / e;
    ALPHA * EPS_0 * e * e / (PI * PI * HBAR) * (1..200).map(|n| (-(n as f64) * x).exp() / (n * n) as f64).sum::<f64>()
}

/// The suppression exponent `πE_S/E`.
pub fn schwinger_exponent(e: f64) -> f64 {
    PI * E_SCHWINGER / e
}

/// Normalised vector potential `a₀ = eE₀/(m_ecω)` of a laser of `intensity` (W/m²) and wavelength
/// `wavelength` (m).
pub fn a0(intensity: f64, wavelength: f64, conv: Intensity) -> f64 {
    let omega = 2.0 * PI * C / wavelength;
    E_CHARGE * peak_field(intensity, conv) / (M_E * C * omega)
}

/// Quantum nonlinearity `χ_e = (1 + β cos θ) γ a₀ ħω/(m_ec²)` of an electron of energy
/// `energy_gev` crossing a laser of strength `a0` and photon energy `photon_ev` at `angle` (rad,
/// 0 = head-on). χ ≈ 1 is the strong-field regime: the field in the electron's frame is E_S.
pub fn chi_e(energy_gev: f64, a0: f64, photon_ev: f64, angle: f64) -> f64 {
    let gamma = energy_gev * 1e6 / MEC2_KEV;
    let beta = (1.0 - 1.0 / (gamma * gamma)).sqrt();
    (1.0 + beta * angle.cos()) * gamma * a0 * photon_ev / (MEC2_KEV * 1e3)
}

/// The Ritus–Narozhny expansion parameter `αχ^(2/3)`; perturbation theory in the field is
/// conjectured to fail where it reaches 1 (χ ≈ 1600). Contested.
pub fn ritus_narozhny(chi: f64) -> f64 {
    ALPHA * chi.powf(2.0 / 3.0)
}

/// Pair production by a photon of quantum parameter `chi` (nonlinear Breit–Wheeler) is
/// suppressed as `e^(−8/(3χ))` for χ ≪ 1.
pub fn breit_wheeler_suppression(chi: f64) -> f64 {
    (-8.0 / (3.0 * chi)).exp()
}

/// Vacuum birefringence `Δn = 3A_eB²`, `A_e = 2α²ħ³/(45μ₀m_e⁴c⁵)`, in a field `b` (T); `None` above
/// 0.1 B_S, where the weak-field (Euler–Heisenberg) expansion fails.
pub fn birefringence(b: f64) -> Option<f64> {
    let a_e = 2.0 * ALPHA * ALPHA * HBAR.powi(3) / (45.0 * MU_0 * M_E.powi(4) * C.powi(5));
    (b <= 0.1 * B_SCHWINGER).then_some(3.0 * a_e * b * b)
}

/// The Euler–Heisenberg correction to the field energy, as a fraction: `(α/45π)(E/E_S)²`.
pub fn euler_heisenberg_fraction(e: f64) -> f64 {
    ALPHA / (45.0 * PI) * (e / E_SCHWINGER).powi(2)
}

/// Binding energy (keV) of the 1s electron of a point nucleus of charge `z`,
/// `m_ec²[1 − √(1 − (Zα)²)]`; `None` above Z = 1/α ≈ 137, where the point-nucleus Dirac equation
/// fails (an extended nucleus reaches the continuum at [`Z_CRITICAL`]).
pub fn dirac_1s_binding(z: f64) -> Option<f64> {
    let x = z * ALPHA;
    (x < 1.0).then(|| MEC2_KEV * (1.0 - (1.0 - x * x).sqrt()))
}

/// The Schwinger field of gapped graphene, a condensed-matter analogue: `Δ²/(eħv_F)` (V/m) for
/// gap `gap_ev` (eV) and Fermi velocity `v_f` (m/s).
pub fn graphene_critical_field(gap_ev: f64, v_f: f64) -> f64 {
    (gap_ev * E_CHARGE).powi(2) / (E_CHARGE * HBAR * v_f)
}

/// The field (T) whose `eB`, in the natural units lattice QCD and heavy-ion papers use, is
/// `eb_gev2` (GeV²): `B = eB/(eħc²)`. 1 GeV² is 1.69×10¹⁶ T and `m_π²` 3.3×10¹⁴ T; lattice QCD
/// finds hadrons still confined at 9 GeV² (1.5×10¹⁷ T; D'Elia et al. 2022).
pub fn field_from_eb(eb_gev2: f64) -> f64 {
    eb_gev2 * GEV.powi(2) / (E_CHARGE * HBAR * C * C)
}

/// The field (T) quoted as `b_ev2` (eV²) in Heaviside–Lorentz natural units, where the energy
/// density is `B²/2`: `B = b·√(μ₀e⁴/(ħc)³)`, so 1 T is 195.35 eV². Quoting B this way, rather
/// than eB, differs by the charge `√(4πα)`: the source of heavy-ion papers' "eB ≈ m_π² ≈ 10¹⁸ G",
/// which is 3.3 times too large as a field.
pub fn field_from_natural(b_ev2: f64) -> f64 {
    b_ev2 * (MU_0 * E_CHARGE.powi(4) / (HBAR * C).powi(3)).sqrt()
}

/// The field (T) whose energy density `B²/2μ₀` equals `Λ⁴/(ħc)³` for a QCD scale `lambda_mev`
/// (MeV): `√(2μ₀Λ⁴/(ħc)³)`. It gives 10¹⁵ T only at Λ = 372 MeV; at PDG 2018's three-flavour
/// Λ⁽³⁾ = 332 ± 17 MeV (eq. 9.24d; later editions quote only α_s) it is 7.98×10¹⁴ T. Equal energy
/// density is not a deconfinement threshold: lattice QCD keeps hadrons confined to fields a
/// hundred times higher.
pub fn lambda4_field(lambda_mev: f64) -> f64 {
    (2.0 * MU_0 * (lambda_mev * MEV).powi(4) / (HBAR * C).powi(3)).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, rel: f64) -> bool {
        (a / b - 1.0).abs() < rel
    }

    #[test]
    fn qcd_field_scales() {
        // eB in GeV²: 1 → 1.6904×10¹⁶ T (Bali et al. 2012); m_π² = (140 MeV)² → 3.3132×10¹⁴ T
        // (Skokov et al. 2009); 9 → 1.520×10¹⁷ T (D'Elia et al. 2022); 10¹⁵ T is 0.0592 GeV².
        assert!(close(field_from_eb(1.0), 1.6904e16, 1e-4) && close(field_from_eb(0.14 * 0.14), 3.3132e14, 1e-4));
        assert!(close(field_from_eb(9.0), 1.520e17, 1e-3) && close(1e15 / field_from_eb(1.0), 0.059_157, 1e-4));
        // Heaviside–Lorentz: 1 T = 195.35 eV²; MoEDAL's Pb–Pb 7.6 and 4.0 GeV² are 3.89 and
        // 2.05×10¹⁶ T; and B quoted as eB differs by e_HL = √(4πα).
        assert!(close(field_from_natural(195.35), 1.0, 1e-4));
        assert!(close(field_from_natural(7.6e18), 3.890e16, 1e-3) && close(field_from_natural(4.0e18), 2.048e16, 1e-3));
        assert!(close(field_from_eb(1.0) / field_from_natural(1e18), 1.0 / (4.0 * PI * ALPHA).sqrt(), 1e-9));
        // The audit's Λ⁴ criterion: 2.90, 6.52 and 8.37×10¹⁴ T at Λ = 200, 300, 340 MeV;
        // 10¹⁵ T needs Λ = 371.67 MeV.
        for (l, b) in [(200.0, 2.8957e14), (300.0, 6.5154e14), (332.0, 7.9794e14), (340.0, 8.3686e14), (371.67, 1e15)] {
            assert!(close(lambda4_field(l), b, 1e-3), "{l} MeV");
        }
    }

    #[test]
    fn critical_fields() {
        assert!(close(E_SCHWINGER, 1.32329e18, 1e-5) && close(B_SCHWINGER, 4.41401e9, 1e-5));
        // 4.6×10²⁹ W/cm² (Di Piazza et al.) with I = ε₀cE²; half that cycle-averaged.
        assert!(close(critical_intensity(Intensity::Peak) * 1e-4, 4.64811e29, 1e-4));
        assert!(close(critical_intensity(Intensity::Linear) * 1e-4, 2.32406e29, 1e-4));
    }

    #[test]
    fn pair_production() {
        let prefactor = ALPHA * EPS_0 * E_SCHWINGER * E_SCHWINGER / (PI * PI * HBAR);
        assert!(close(prefactor, 1.08704e56, 1e-4));
        // π × prefactor in cm⁻³ s⁻¹ is Melissinos's 3.45×10⁵⁰ (ratio 0.990).
        assert!(close(PI * prefactor * 1e-6, 3.415e50, 1e-3));
        assert!(close(vacuum_decay_rate(E_SCHWINGER) / schwinger_rate(E_SCHWINGER), 1.011016, 1e-5));
        assert!(close(schwinger_rate(0.05 * E_SCHWINGER), 1.4017e26, 1e-3));
        assert!(close(schwinger_rate(0.01 * E_SCHWINGER), 3.9683e-85, 1e-3));
    }

    #[test]
    fn lasers_and_accelerators() {
        // CoReLS 2021: 1.1×10²³ W/cm² at 0.8 µm, peak field 9.10×10¹⁴ V/m = E_S/1454.
        let i = 1.1e27;
        assert!(close(a0(i, 0.8e-6, Intensity::Linear), 226.84, 1e-4) && close(a0(i, 0.8e-6, Intensity::Peak), 160.40, 1e-4));
        assert!(close(E_SCHWINGER / peak_field(i, Intensity::Linear), 1454.0, 1e-3));
        assert!(close(schwinger_exponent(9.104e14), 4566.0, 1e-3));
        // LUXE: 16.5 GeV electrons, a₀ = 7.9, 1.55 eV photons at 17.2° from head-on: χ = 1.5.
        assert!(close(chi_e(16.5, 7.9, 1.55, 17.2f64.to_radians()), 1.5127, 1e-3));
        assert!(close(ritus_narozhny(1604.18), 1.0, 1e-4));
        assert!(close(breit_wheeler_suppression(1.0), 6.9483e-2, 1e-4));
    }

    #[test]
    fn vacuum_optics_and_heavy_atoms() {
        // PVLAS's field, 2.5 T: Δn = 2.48×10⁻²³ (QED); PVLAS measured (12 ± 17)×10⁻²³.
        assert!(close(birefringence(2.5).unwrap(), 2.48375e-23, 1e-4));
        assert!(birefringence(B_SCHWINGER).is_none());
        assert!(close(euler_heisenberg_fraction(0.1 * E_SCHWINGER), 5.16e-7, 2e-3));
        assert!(close(dirac_1s_binding(82.0).unwrap(), 101.581, 1e-4) && close(dirac_1s_binding(92.0).unwrap(), 132.28, 1e-4));
        assert!(dirac_1s_binding(138.0).is_none());
        // Gapped graphene: 0.2 eV and 10⁶ m/s give 6×10⁷ V/m (Schmitt et al. 2023: ≲ 6×10⁷).
        assert!(close(graphene_critical_field(0.2, 1e6), 6.077e7, 1e-3));
    }
}
