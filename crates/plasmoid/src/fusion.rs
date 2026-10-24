//! Getting nuclei through the Coulomb barrier: tunnelling and Gamow peaks, reaction rates,
//! ignition, electron screening, muon catalysis and gravitational confinement in stars, with
//! the checks that test claims of low-energy fusion.
//!
//! Energies are in keV and reaction rates in m³/s unless a name says otherwise.
//!
//! ```
//! use plasmoid::fusion::{Reaction, bosch_hale, gamow_peak, gamow_energy, M_D, M_T, reduced_mass};
//!
//! // D–T at 10 keV: the reactions happen near 31 keV, three times the temperature.
//! let e_g = gamow_energy(1.0, 1.0, reduced_mass(M_D, M_T));
//! let (peak, _) = gamow_peak(e_g, 10.0);
//! assert!((peak - 30.9).abs() < 0.1);
//! // Bosch & Hale's reactivity there: 1.136×10⁻²² m³/s.
//! assert!((bosch_hale(Reaction::DT, 10.0).unwrap() / 1.1362e-22 - 1.0).abs() < 1e-3);
//! ```

use core::f64::consts::PI;

/// Fine-structure constant (CODATA 2018).
pub const ALPHA: f64 = 7.297_352_569_3e-3;
/// `e²/4πε₀`, keV fm.
pub const E2: f64 = 1.439_964_548;
/// keV in joules.
pub const KEV: f64 = 1.602_176_634e-16;
/// Atomic mass unit, keV.
pub const AMU_KEV: f64 = 931_494.102_42;
/// Nuclear masses (CODATA 2018), keV.
pub const M_P: f64 = 938_272.088_16;
/// Deuteron, keV.
pub const M_D: f64 = 1_875_612.942_4;
/// Triton, keV.
pub const M_T: f64 = 2_808_921.132_5;
/// Helion (³He nucleus), keV.
pub const M_HE3: f64 = 2_808_391.607_4;
/// ¹¹B nucleus: the atomic mass less five electrons and their 0.67 keV binding, keV.
pub const M_B11: f64 = 11.009_305_166 * AMU_KEV - 5.0 * 510.998_95 + 0.6704;

/// Reduced mass of two particles, keV.
pub fn reduced_mass(m1: f64, m2: f64) -> f64 {
    m1 * m2 / (m1 + m2)
}

/// Height of the Coulomb barrier (keV) between charges `z1`, `z2` touching at `r_fm` (fm):
/// `1.44 MeV fm · Z₁Z₂/R`. The radius is the caller's convention.
pub fn coulomb_barrier(z1: f64, z2: f64, r_fm: f64) -> f64 {
    E2 * 1e3 * z1 * z2 / r_fm
}

/// Sommerfeld parameter `η = Z₁Z₂α√(m_rc²/2E)` for reduced mass `mr` and energy `e` (keV).
pub fn sommerfeld(z1: f64, z2: f64, mr: f64, e: f64) -> f64 {
    z1 * z2 * ALPHA * (mr / (2.0 * e)).sqrt()
}

/// Tunnelling (Gamow) factor `e^(−2πη)`.
pub fn tunnelling(z1: f64, z2: f64, mr: f64, e: f64) -> f64 {
    (-2.0 * PI * sommerfeld(z1, z2, mr, e)).exp()
}

/// Gamow energy `E_G = 2m_rc²(παZ₁Z₂)²`, keV (`e^(−2πη) = e^(−√(E_G/E))`).
pub fn gamow_energy(z1: f64, z2: f64, mr: f64) -> f64 {
    2.0 * mr * (PI * ALPHA * z1 * z2).powi(2)
}

/// The Gamow peak and its width at temperature `kt` (keV): `E₀ = (E_G kT²/4)^(1/3)`,
/// `Δ = 4√(E₀kT/3)`.
pub fn gamow_peak(e_g: f64, kt: f64) -> (f64, f64) {
    let e0 = (e_g * kt * kt / 4.0).cbrt();
    (e0, 4.0 * (e0 * kt / 3.0).sqrt())
}

/// The reactions Bosch & Hale (1992) fit.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Reaction {
    /// D + T → ⁴He + n (17.59 MeV).
    DT,
    /// D + D → ³He + n (3.27 MeV).
    DDn,
    /// D + D → T + p (4.03 MeV).
    DDp,
    /// D + ³He → ⁴He + p (18.35 MeV).
    DHe3,
}

// (B_G keV^½, m_rc² keV, C1…C7 with C1 in cm³/s, validity keV), Bosch & Hale 1992 Table VII.
fn bh(r: Reaction) -> (f64, f64, [f64; 7], (f64, f64)) {
    match r {
        Reaction::DT => (34.3827, 1_124_656.0, [1.17302e-9, 1.51361e-2, 7.51886e-2, 4.60643e-3, 1.35000e-2, -1.06750e-4, 1.36600e-5], (0.2, 100.0)),
        Reaction::DDn => (31.3970, 937_814.0, [5.43360e-12, 5.85778e-3, 7.68222e-3, 0.0, -2.96400e-6, 0.0, 0.0], (0.2, 100.0)),
        Reaction::DDp => (31.3970, 937_814.0, [5.65718e-12, 3.41267e-3, 1.99167e-3, 0.0, 1.05060e-5, 0.0, 0.0], (0.2, 100.0)),
        Reaction::DHe3 => (68.7508, 1_124_572.0, [5.51036e-10, 6.41918e-3, -2.02896e-3, -1.91080e-5, 1.35776e-4, 0.0, 0.0], (0.5, 190.0)),
    }
}

/// The Bosch–Hale form `C₁θ√(ξ/(m_rc²T³))e^(−3ξ)` with θ and ξ built from C₂…C₇.
fn bosch_hale_form(e_g: f64, mr: f64, c: &[f64; 7], t: f64) -> f64 {
    let theta = t / (1.0 - t * (c[1] + t * (c[3] + t * c[5])) / (1.0 + t * (c[2] + t * (c[4] + t * c[6]))));
    let xi = (e_g / (4.0 * theta)).cbrt();
    c[0] * theta * (xi / (mr * t.powi(3))).sqrt() * (-3.0 * xi).exp()
}

/// Maxwellian reactivity ⟨σv⟩ (m³/s) at ion temperature `t` (keV), from Bosch & Hale (1992);
/// `None` outside the fit's range (0.2–100 keV; D–³He 0.5–190 keV).
pub fn bosch_hale(r: Reaction, t: f64) -> Option<f64> {
    let (b_g, mr, c, (lo, hi)) = bh(r);
    (lo..=hi).contains(&t).then(|| 1e-6 * bosch_hale_form(b_g * b_g, mr, &c, t))
}

/// Cross section (mb) at centre-of-mass energy `e` (keV): Bosch & Hale's S-factor Padé fit,
/// `σ = S(E)/(E e^(B_G/√E))`. D–T is the low-energy fit, valid to 550 keV.
pub fn bosch_hale_sigma(r: Reaction, e: f64) -> f64 {
    let (b_g, a, b): (f64, [f64; 5], [f64; 4]) = match r {
        Reaction::DT => (34.3827, [6.927e4, 7.454e8, 2.050e6, 5.2002e4, 0.0], [6.38e1, -9.95e-1, 6.981e-5, 1.728e-4]),
        Reaction::DDp => (31.3970, [5.5576e4, 2.1054e2, -3.2638e-2, 1.4987e-6, 1.8181e-10], [0.0; 4]),
        Reaction::DDn => (31.3970, [5.3701e4, 3.3027e2, -1.2706e-1, 2.9327e-5, -2.5151e-9], [0.0; 4]),
        Reaction::DHe3 => (68.7508, [5.7501e6, 2.5226e3, 4.5566e1, 0.0, 0.0], [-3.1995e-3, -8.5530e-6, 5.9014e-8, 0.0]),
    };
    let s = (a[0] + e * (a[1] + e * (a[2] + e * (a[3] + e * a[4])))) / (1.0 + e * (b[0] + e * (b[1] + e * (b[2] + e * b[3]))));
    s / (e * (b_g / e.sqrt()).exp())
}

/// The p–¹¹B reactivity fits.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Pb11Fit {
    /// Nevins & Swain (2000), without its 148 keV resonance term, as published tables use it.
    NevinsSwain,
    /// Tentori & Belloni (2023), from the higher Sikora–Weller cross section.
    TentoriBelloni,
}

/// p + ¹¹B → 3 ⁴He reactivity (m³/s) at `t` (keV), in the Bosch–Hale form; `None` outside
/// 50–500 keV (Nevins–Swain) or 70–500 keV (Tentori–Belloni).
pub fn pb11_reactivity(fit: Pb11Fit, t: f64) -> Option<f64> {
    let (c, lo) = match fit {
        Pb11Fit::NevinsSwain => ([4.4467e-14, -5.9357e-2, 2.0165e-1, 1.0404e-3, 2.7621e-3, -9.1653e-6, 9.8305e-7], 50.0),
        Pb11Fit::TentoriBelloni => ([9.7827e-14, -5.1610e-2, 1.3240e-1, 3.8446e-4, 1.2499e-3, -5.6715e-6, 1.1615e-6], 70.0),
    };
    let mr = reduced_mass(M_P, M_B11);
    (lo..=500.0).contains(&t).then(|| bosch_hale_form(gamow_energy(1.0, 5.0, mr), mr, &c, t))
}

/// Bremsstrahlung coefficient `C_B` (W m³ keV^−½) of `P = C_B Z_eff n_e² √T`, as Wurzel & Hsu use.
pub const C_BREMS: f64 = 5.34e-37;

/// The Lawson product nτ_E (s/m³) a D–T plasma at `t` (keV) needs for fuel gain `q` (use
/// `f64::INFINITY` for ignition): `3T/[(f_α + 1/Q)⟨σv⟩E_F/4 − C_B√T]` (Wurzel & Hsu 2022).
/// `None` where bremsstrahlung wins or outside the reactivity fit.
pub fn lawson(t: f64, q: f64) -> Option<f64> {
    let sv = bosch_hale(Reaction::DT, t)?;
    let den = (3.5 / 17.6 + 1.0 / q) * sv * 17.6e3 * KEV / 4.0 - C_BREMS * t.sqrt();
    (den > 0.0).then(|| 3.0 * t * KEV / den)
}

/// How electron screening `ue` (keV) enhances a rate at energy `e`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Screening {
    /// `e^(πηU_e/E)`, valid for U_e ≪ E.
    FirstOrder,
    /// The tunnelling factor at E + U_e over that at E.
    Exact,
}

/// Rate enhancement from electron screening of `ue` (keV) at energy `e` (keV).
pub fn screening(z1: f64, z2: f64, mr: f64, e: f64, ue: f64, form: Screening) -> f64 {
    match form {
        Screening::FirstOrder => (PI * sommerfeld(z1, z2, mr, e) * ue / e).exp(),
        Screening::Exact => tunnelling(z1, z2, mr, e + ue) / tunnelling(z1, z2, mr, e),
    }
}

/// Salpeter's weak-screening factor in a star: `e^(0.188 Z₁Z₂ ζ √ρ T₆^(−3/2))`, with density
/// `rho` in g/cm³ and temperature `t6` in 10⁶ K.
pub fn salpeter(z1z2: f64, zeta: f64, rho: f64, t6: f64) -> f64 {
    (0.188 * z1z2 * zeta * rho.sqrt() * t6.powf(-1.5)).exp()
}

/// Muon lifetime, s.
pub const MUON_LIFETIME: f64 = 2.196_981_1e-6;

/// Fusions per muon: `1/(ω_s + λ₀/λ_c)` for alpha-sticking probability `sticking` and cycling
/// rate `cycle_rate` (s⁻¹).
pub fn muon_yield(sticking: f64, cycle_rate: f64) -> f64 {
    1.0 / (sticking + 1.0 / (MUON_LIFETIME * cycle_rate))
}

/// Fusions per muon needed to repay the energy `muon_cost_mev` spent making each muon
/// (17.6 MeV per D–T fusion).
pub fn muon_breakeven(muon_cost_mev: f64) -> f64 {
    muon_cost_mev / 17.6
}

/// D–D reactions per second per watt, from the mean of the two branches (3.651 MeV); half of
/// them make a neutron.
pub fn dd_per_watt() -> f64 {
    1.0 / ((3.268 + 4.033) / 2.0 * 1e3 * KEV)
}

/// Neutrons per second that `watts` of D–D fusion must make, and how many times more that is
/// than `claimed` neutrons per second.
pub fn dd_neutron_deficit(watts: f64, claimed: f64) -> (f64, f64) {
    let expected = watts * dd_per_watt() / 2.0;
    (expected, expected / claimed)
}

/// Fusion power density (W/m³) of densities `n1`, `n2` (m⁻³) reacting at ⟨σv⟩ `sv` (m³/s) and
/// releasing `q_mev` each; `identical` halves it for like particles (D–D).
pub fn fusion_power(n1: f64, n2: f64, sv: f64, q_mev: f64, identical: bool) -> f64 {
    n1 * n2 * sv * q_mev * 1e3 * KEV / if identical { 2.0 } else { 1.0 }
}

/// Surface field (V/m), potential (V) and stored energy (J) of a charge `q` (C) on a sphere of
/// radius `r` (m).
pub fn charged_sphere(q: f64, r: f64) -> (f64, f64, f64) {
    const K: f64 = 8.987_551_792_3e9;
    (K * q / (r * r), K * q / r, K * q * q / (2.0 * r))
}

/// Magnetic field (T) at the centre of a spherical shell of charge `q` (C) and radius `r` (m)
/// spinning at `omega` (rad/s): `(2/3)μ₀Qω/(4πR)`.
pub fn spinning_shell_field(q: f64, r: f64, omega: f64) -> f64 {
    const MU0: f64 = 1.256_637_062_12e-6;
    2.0 / 3.0 * MU0 * q * omega / (4.0 * PI * r)
}

/// The Widom–Larsen electron-capture rate as its authors give it, `1.86×10⁻³(β − β₀)²` Hz, for
/// mass-enhancement `beta` above threshold `beta0`. Ciuchi et al. (2012) get about 300× less.
pub fn widom_larsen_rate(beta: f64, beta0: f64) -> f64 {
    1.86e-3 * (beta - beta0).powi(2)
}

/// Gravitational confinement: the Lane–Emden polytrope of index `n`, integrated to its first
/// zero. Returns `(ξ₁, −ξ₁²θ′(ξ₁), ρ_c/ρ̄)`.
pub fn lane_emden(n: f64) -> (f64, f64, f64) {
    // Series start, then RK4 on (θ, θ′).
    let h = 1e-4;
    let mut xi = 1e-3;
    let (mut th, mut dth) = (1.0 - xi * xi / 6.0, -xi / 3.0);
    let f = |xi: f64, th: f64, dth: f64| -> (f64, f64) { (dth, -th.max(0.0).powf(n) - 2.0 * dth / xi) };
    loop {
        let (k1a, k1b) = f(xi, th, dth);
        let (k2a, k2b) = f(xi + h / 2.0, th + h / 2.0 * k1a, dth + h / 2.0 * k1b);
        let (k3a, k3b) = f(xi + h / 2.0, th + h / 2.0 * k2a, dth + h / 2.0 * k2b);
        let (k4a, k4b) = f(xi + h, th + h * k3a, dth + h * k3b);
        let (nth, ndth) = (th + h / 6.0 * (k1a + 2.0 * k2a + 2.0 * k3a + k4a), dth + h / 6.0 * (k1b + 2.0 * k2b + 2.0 * k3b + k4b));
        if nth <= 0.0 {
            let frac = th / (th - nth);
            let (x1, d1) = (xi + frac * h, dth + frac * (ndth - dth));
            return (x1, -x1 * x1 * d1, -x1 / (3.0 * d1));
        }
        (xi, th, dth) = (xi + h, nth, ndth);
    }
}

/// Kelvin–Helmholtz time (s): how long gravity alone could power a star of mass `m` (kg),
/// radius `r` (m) and luminosity `l` (W), `GM²/(2RL)`.
pub fn kelvin_helmholtz(m: f64, r: f64, l: f64) -> f64 {
    6.674_30e-11 * m * m / (2.0 * r * l)
}

/// Virial temperature (K) of a star of mass `m` (kg) and radius `r` (m) with mean molecular
/// weight `mu`: `GMμm_u/(k_B R)`.
pub fn virial_temperature(m: f64, r: f64, mu: f64) -> f64 {
    6.674_30e-11 * m * mu * 1.660_539_066_60e-27 / (1.380_649e-23 * r)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, rel: f64) -> bool {
        (a / b - 1.0).abs() < rel
    }

    #[test]
    fn barriers_and_tunnelling() {
        assert!(close(coulomb_barrier(1.0, 1.0, 4.0), 360.0, 1e-3));
        assert!(close(coulomb_barrier(1.0, 1.0, 3.243), 444.0, 1e-3));
        assert!(close(coulomb_barrier(1.0, 5.0, 4.514), 1595.0, 1e-3));
        let dt = reduced_mass(M_D, M_T);
        assert!(close(sommerfeld(1.0, 1.0, dt, 10.0), 1.7304, 1e-4));
        assert!(close(tunnelling(1.0, 1.0, dt, 10.0), 1.8969e-5, 1e-3));
        let dd = reduced_mass(M_D, M_D);
        assert!(close(tunnelling(1.0, 1.0, dd, 1.0), 2.3147e-14, 1e-3));
        let pb = reduced_mass(M_P, M_B11);
        assert!(close(tunnelling(1.0, 5.0, pb, 10.0), 2.2848e-21, 2e-3));
        // Gamow energies and peaks; √E_G(D–T) is Bosch & Hale's B_G = 34.3827.
        assert!(close(gamow_energy(1.0, 1.0, dt), 1182.2, 1e-4));
        assert!(close(gamow_energy(1.0, 1.0, dt).sqrt(), 34.3827, 1e-5));
        assert!(close(gamow_energy(1.0, 5.0, pb), 22589.1, 1e-4));
        let (e0, w) = gamow_peak(gamow_energy(1.0, 1.0, dt), 10.0);
        assert!(close(e0, 30.918, 1e-4) && close(w, 40.607, 1e-4));
        let (e0, w) = gamow_peak(gamow_energy(1.0, 1.0, reduced_mass(M_P, M_P)), 1.3357);
        assert!(close(e0, 6.036, 1e-3) && close(w, 6.557, 1e-3));
    }

    #[test]
    fn bosch_hale_reproduces_its_tables() {
        // UWFDM-1268 Table 6 (cm³/s), within 0.3%.
        let table = [
            (1.0, [6.86e-21, 9.93e-23, 1.02e-22, 3.05e-26]),
            (10.0, [1.14e-16, 6.02e-19, 5.78e-19, 2.12e-19]),
            (100.0, [8.45e-16, 2.68e-17, 2.24e-17, 1.72e-16]),
        ];
        for (t, row) in table {
            for (r, v) in [Reaction::DT, Reaction::DDn, Reaction::DDp, Reaction::DHe3].iter().zip(row) {
                let sv = bosch_hale(*r, t).unwrap() * 1e6;
                assert!(close(sv, v, 4e-3), "{r:?} {t}: {sv:e}");
            }
        }
        assert!(close(bosch_hale(Reaction::DT, 10.0).unwrap(), 1.1362e-22, 1e-4));
        assert!(close(bosch_hale(Reaction::DDn, 0.1 * 2.0).unwrap() / bosch_hale(Reaction::DDn, 0.2).unwrap(), 1.0, 1e-12));
        assert!(bosch_hale(Reaction::DT, 150.0).is_none() && bosch_hale(Reaction::DHe3, 150.0).is_some());
        // Cross sections: D–T peaks at 5.07 b near 65 keV; 27 mb at 10 keV.
        assert!(close(bosch_hale_sigma(Reaction::DT, 10.0), 27.02, 1e-3));
        let peak = (100..1300).map(|k| bosch_hale_sigma(Reaction::DT, f64::from(k) * 0.1)).fold(0.0, f64::max);
        assert!(close(peak, 5065.7, 1e-3));
        assert!(close(bosch_hale_sigma(Reaction::DDn, 50.0), 16.487, 1e-3));
        assert!(close(bosch_hale_sigma(Reaction::DDp, 50.0), 15.573, 1e-3));
    }

    #[test]
    fn cross_sections_average_to_the_reactivity() {
        // ∫σv over a Maxwellian matches the reactivity fit to about 1%.
        let mr = 1_124_656.0 * KEV / (299_792_458.0f64).powi(2);
        let t = 10.0;
        let n = 20_000;
        let e_max = 40.0 * t;
        let mut sum = 0.0;
        for k in 1..n {
            let e = e_max * k as f64 / n as f64;
            sum += bosch_hale_sigma(Reaction::DT, e) * 1e-31 * e * (-e / t).exp();
        }
        let sv = (8.0 / (PI * mr)).sqrt() * (t * KEV).powf(-1.5) * sum * (e_max / n as f64) * KEV * KEV;
        assert!(close(sv, bosch_hale(Reaction::DT, t).unwrap(), 0.012), "{sv:e}");
    }

    #[test]
    fn proton_boron() {
        let tb = |t| pb11_reactivity(Pb11Fit::TentoriBelloni, t).unwrap();
        let ns = |t| pb11_reactivity(Pb11Fit::NevinsSwain, t).unwrap();
        assert!(close(tb(100.0), 7.0650e-23, 2e-3) && close(tb(300.0), 4.3870e-22, 2e-3) && close(tb(500.0), 5.6190e-22, 2e-3));
        // Published NS tables at 75 and 500 keV; the new fit is 30% and 50% higher at 300 and 500.
        assert!(close(ns(75.0), 2.455e-23, 2e-3) && close(ns(500.0), 3.741e-22, 2e-3));
        assert!((tb(300.0) / ns(300.0) - 1.296).abs() < 0.005 && (tb(500.0) / ns(500.0) - 1.502).abs() < 0.005);
        assert!(pb11_reactivity(Pb11Fit::TentoriBelloni, 60.0).is_none());
    }

    #[test]
    fn lawson_criterion() {
        // Ignition: minimum nTτ 2.97×10²¹ keV s/m³ near 14 keV; ignition temperature 4.3 keV.
        let ts: Vec<f64> = (0..1981).map(|k| 1.0 + 0.05 * f64::from(k)).filter(|t| *t <= 100.0).collect();
        let (t_min, ntt) = ts.iter().filter_map(|&t| lawson(t, f64::INFINITY).map(|n| (t, n * t))).min_by(|a, b| a.1.total_cmp(&b.1)).unwrap();
        assert!(close(ntt, 2.968e21, 2e-3) && (t_min - 14.25).abs() < 0.06, "{t_min} {ntt:e}");
        let (t1, nt1) = ts.iter().filter_map(|&t| lawson(t, 1.0).map(|n| (t, n * t))).min_by(|a, b| a.1.total_cmp(&b.1)).unwrap();
        assert!(close(nt1, 4.677e20, 2e-3) && (t1 - 13.65).abs() < 0.06);
        assert!(lawson(4.25, f64::INFINITY).is_none() && lawson(4.35, f64::INFINITY).is_some());
    }

    #[test]
    fn screening_and_stars() {
        let dd = reduced_mass(M_D, M_D);
        assert!(close(screening(1.0, 1.0, dd, 10.0, 0.025, Screening::FirstOrder), 1.0125, 1e-3));
        assert!(close(screening(1.0, 1.0, dd, 1.0, 0.3, Screening::FirstOrder), 111.0, 2e-3));
        assert!(close(screening(1.0, 1.0, dd, 1.0, 0.3, Screening::Exact), 47.47, 2e-3));
        let dhe = reduced_mass(M_D, M_HE3);
        assert!(close(screening(1.0, 2.0, dhe, 1.0, 0.219, Screening::FirstOrder), 1862.0, 3e-3));
        assert!(close(salpeter(1.0, 1.291, 150.0, 15.7), 1.0489, 1e-4));
        let (xi, mass, conc) = lane_emden(3.0);
        assert!(close(xi, 6.897, 1e-3) && close(mass, 2.018, 1e-3) && close(conc, 54.18, 2e-3), "{xi} {mass} {conc}");
        assert!(close(kelvin_helmholtz(1.989e30, 6.957e8, 3.828e26) / 3.156e13, 15.7, 0.01));
        assert!(close(virial_temperature(1.989e30, 6.957e8, 0.61), 1.40e7, 0.01));
    }

    #[test]
    fn muons_and_claims() {
        assert!(close(muon_yield(0.0047, 2.5287e8), 153.85, 1e-3));
        assert!(close(muon_breakeven(4700.0), 267.0, 2e-3) && close(muon_breakeven(6000.0), 341.0, 2e-3));
        assert!(close(dd_per_watt(), 1.710e12, 2e-3));
        // Fleischmann & Pons: 1.751 W of excess heat as D–D needs 1.5×10¹² n/s; they saw 4×10⁴.
        let (expected, deficit) = dd_neutron_deficit(1.751, 4e4);
        assert!(close(expected, 1.497e12, 2e-3) && close(deficit, 3.74e7, 2e-3));
        // Pais's fusion device: D–D at 20 keV, 10²⁰ m⁻³, in a 0.15 m sphere gives ~200 W, not GW.
        let p = fusion_power(1e20, 1e20, 5.002e-24, 3.651, true);
        assert!(close(p, 1.463e4, 2e-3) && close(p * 4.0 / 3.0 * PI * 0.15f64.powi(3), 207.0, 0.01));
        let (field, volts, energy) = charged_sphere(1.0, 0.15);
        assert!(close(field, 3.99e11, 2e-3) && close(volts, 5.99e10, 2e-3) && close(energy, 3.00e10, 2e-3));
        assert!(close(spinning_shell_field(1.0, 0.15, 1e4), 4.44e-3, 2e-3));
        assert!(close(widom_larsen_rate(20.0, 2.59), 0.564, 2e-3));
    }
}
