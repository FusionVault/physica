//! Distances in a Friedmann universe: what a lensed quasar's time delays, a standard siren, a
//! supernova or the CMB's acoustic scale measure.
//!
//! [`Cosmology`] carries curvature, CPL dark energy and radiation, with massive neutrinos through
//! the exact Fermi–Dirac integral, and [`PLANCK18`] is astropy's realization of Planck 2018. Its
//! `E(z) = H(z)/H₀` is
//! `√(Ω_m(1+z)³ + Ω_k(1+z)² + Ω_de f_de(z) + Ω_γ(1+z)⁴[1 + ρ_ν/ρ_γ(z)])`, the line-of-sight
//! distance `D_C = (c/H₀)∫dz/E` and the transverse distance `D_M = (c/H₀)S_k(D_C H₀/c)`, with
//! `S_k(x) = sinh(√Ω_k x)/√Ω_k` (open), `x` (flat) or `sin(√−Ω_k x)/√−Ω_k` (closed)
//! ([Hogg 1999](https://arxiv.org/abs/astro-ph/9905116), Eqs. 14–16). Between redshifts the
//! angular diameter distance is `S_k(χ₂ − χ₁)/(1 + z₂)`; Hogg's Eq. 19 gives the same only for
//! `Ω_k ≥ 0` and picks the wrong branch in a closed universe once `√−Ω_k χ > π/2`.
//!
//! The free functions are the flat, radiation-free ΛCDM case, which is all a time-delay lens
//! needs: a lens at `z_d` in front of a source at `z_s` has the time-delay distance
//! `(1 + z_d)D_d D_s/D_ds`, which scales as `1/H₀` ([Wong et al. 2020](https://arxiv.org/abs/1907.04869)).
//!
//! ```
//! use gravitational::H100;
//! use gravitational::cosmology::{PLANCK18, angular_diameter_distance, time_delay_distance};
//! use metrology::MPC;
//!
//! // RXJ1131−1231 (z_d = 0.295, z_s = 0.654) at H₀ = 70 km/s/Mpc and Ω_m = 0.3: 2334 Mpc,
//! // against the measured 2096 (+98/−83), so this lens alone says H₀ = 70 × 2334/2096 = 78.
//! let ddt = time_delay_distance(0.295, 0.654, 0.7 * H100, 0.3) / MPC;
//! assert!((ddt - 2334.06).abs() < 0.01);
//! assert!(angular_diameter_distance(0.0, 0.295, 0.7 * H100, 0.3) / MPC > 800.0);
//!
//! // Radiation sets the CMB's angle: Planck's sound horizon r* = 144.57 Mpc at z* = 1089.80
//! // subtends 100θ* = 1.04124, against the measured 1.04119 ± 0.00029.
//! let theta = 100.0 * 144.57 * MPC / PLANCK18.transverse_distance(0.0, 1089.80);
//! assert!((theta - 1.04119).abs() < 0.00029);
//! ```

use core::f64::consts::PI;

use metrology::{C, EV, G, HBAR, K_B};

use crate::{H100, integrate_log};

/// A Friedmann–Lemaître universe. Curvature is whatever the densities leave:
/// `Ω_k = 1 − Ω_m − Ω_de − Ω_r`, with `Ω_r` set by `t_cmb`, `n_eff` and `m_nu`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Cosmology {
    /// The Hubble constant, s⁻¹ (e.g. `0.7 * H100`).
    pub h0: f64,
    /// Cold matter (baryons and dark matter, without the massive neutrinos) today.
    pub omega_m: f64,
    /// Dark energy today.
    pub omega_de: f64,
    /// The dark energy's equation of state today, `w(a) = w₀ + w_a(1 − a)` (Chevallier–Polarski–Linder).
    pub w0: f64,
    /// Its rate of change, `w_a`.
    pub wa: f64,
    /// The CMB temperature today, K; 0 switches all radiation off, as in astropy.
    pub t_cmb: f64,
    /// The effective number of neutrino species, shared equally by the three in `m_nu`.
    pub n_eff: f64,
    /// The three neutrino masses as rest energies, J (e.g. `0.06 * EV`).
    pub m_nu: [f64; 3],
}

/// Planck 2018 (TT,TE,EE+lowE+lensing+BAO, Table 2) as astropy's `Planck18`: H₀ = 67.66,
/// Ω_m = 0.30966 without the 0.06 eV neutrino (0.3111 with it), T_CMB = 2.7255 K, N_eff = 3.046,
/// and flat, so `omega_de` is `1 − Ω_m − Ω_γ − Ω_ν` ([Planck 2018 VI](https://arxiv.org/abs/1807.06209);
/// [astropy](https://github.com/astropy/astropy/blob/main/astropy/cosmology/data/Planck18.ecsv)).
pub const PLANCK18: Cosmology = Cosmology { h0: 0.6766 * H100, omega_m: 0.30966, omega_de: 0.688_846_270_597_353_6, w0: -1.0, wa: 0.0, t_cmb: 2.7255, n_eff: 3.046, m_nu: [0.0, 0.0, 0.06 * EV] };

impl Cosmology {
    /// Flat ΛCDM without radiation: the universe of the free functions.
    const fn flat(h0: f64, omega_m: f64) -> Self {
        Self { h0, omega_m, omega_de: 1.0 - omega_m, w0: -1.0, wa: 0.0, t_cmb: 0.0, n_eff: 0.0, m_nu: [0.0; 3] }
    }

    /// The photon density today, `Ω_γ = (π²/15)(k_BT)⁴/(ħc)³ ÷ 3H₀²c²/(8πG)`.
    fn omega_photons(&self) -> f64 {
        PI * PI / 15.0 * (K_B * self.t_cmb).powi(4) / (HBAR * C).powi(3) * 8.0 * PI * G / (3.0 * self.h0 * self.h0 * C * C)
    }

    /// `ρ_ν/ρ_γ` at redshift `z`: `(7/8)(4/11)^{4/3}(N_eff/3)` per species, times its exact
    /// Fermi–Dirac ratio `F(y)/F(0)`, `F(y) = ∫x²√(x² + y²)/(eˣ + 1)dx`, `F(0) = 7π⁴/120`, with
    /// `y = m_νc²/(k_BT_ν)` and `T_ν = (4/11)^{1/3}T_CMB(1 + z)` (Komatsu et al. 2011, Eqs. 24–25).
    fn neutrino_ratio(&self, z: f64) -> f64 {
        let t_nu = (4.0f64 / 11.0).cbrt() * self.t_cmb * (1.0 + z);
        let fermi = |m: f64| {
            let y = m / (K_B * t_nu);
            // e^−60 is past round-off; below 10⁻⁶ the integrand is x²y, worth 10⁻¹⁹y.
            if y == 0.0 { 1.0 } else { integrate_log(1e-6, 60.0, 600, |x| x * x * (x * x + y * y).sqrt() / (x.exp() + 1.0)) / (7.0 * PI.powi(4) / 120.0) }
        };
        7.0 / 8.0 * (4.0f64 / 11.0).powf(4.0 / 3.0) * self.n_eff / 3.0 * self.m_nu.iter().map(|&m| fermi(m)).sum::<f64>()
    }

    /// The radiation density `Ω_γ(1 + ρ_ν/ρ_γ)` at redshift `z`, in units of today's critical density
    /// and without the `(1 + z)⁴`.
    fn radiation(&self, z: f64) -> f64 {
        if self.t_cmb == 0.0 { 0.0 } else { self.omega_photons() * (1.0 + self.neutrino_ratio(z)) }
    }

    /// The curvature density `Ω_k = 1 − Ω_m − Ω_de − Ω_r`; positive is open.
    pub fn curvature(&self) -> f64 {
        1.0 - self.omega_m - self.omega_de - self.radiation(0.0)
    }

    /// `E(z) = H(z)/H₀`.
    pub fn e(&self, z: f64) -> f64 {
        self.e_with(z, self.curvature())
    }

    fn e_with(&self, z: f64, omega_k: f64) -> f64 {
        let x = 1.0 + z;
        let de = x.powf(3.0 * (1.0 + self.w0 + self.wa)) * (-3.0 * self.wa * z / x).exp();
        (x * x * (self.omega_m * x + omega_k + self.radiation(z) * x * x) + self.omega_de * de).sqrt()
    }

    /// The line-of-sight comoving distance (m) from redshift `z1` to `z2`, `(c/H₀)∫dz/E(z)`.
    pub fn comoving_distance(&self, z1: f64, z2: f64) -> f64 {
        let omega_k = self.curvature();
        // Over x = 1 + z, so the log-grid Simpson rule starts above zero and reaches z* ≈ 1090.
        C / self.h0 * integrate_log(1.0 + z1, 1.0 + z2, 2000, |x| 1.0 / self.e_with(x - 1.0, omega_k))
    }

    /// The transverse comoving distance (m) of `z2` seen from `z1`: `(c/H₀)S_k(χ₂ − χ₁)`.
    pub fn transverse_distance(&self, z1: f64, z2: f64) -> f64 {
        let (omega_k, dh) = (self.curvature(), C / self.h0);
        let chi = self.comoving_distance(z1, z2) / dh;
        let s = omega_k.abs().sqrt();
        dh * if omega_k > 0.0 { (s * chi).sinh() / s } else if omega_k < 0.0 { (s * chi).sin() / s } else { chi }
    }

    /// The angular diameter distance (m) of an object at `z2` seen from `z1 < z2` (`z1 = 0` for
    /// the observer): `D_M(z₁, z₂)/(1 + z₂)`.
    pub fn angular_diameter_distance(&self, z1: f64, z2: f64) -> f64 {
        self.transverse_distance(z1, z2) / (1.0 + z2)
    }

    /// The luminosity distance (m) to redshift `z`, `(1 + z)D_M`.
    pub fn luminosity_distance(&self, z: f64) -> f64 {
        (1.0 + z) * self.transverse_distance(0.0, z)
    }

    /// The time-delay distance (m) of a lens at `zd` and a source at `zs`: `(1 + z_d)D_d D_s/D_ds`.
    pub fn time_delay_distance(&self, zd: f64, zs: f64) -> f64 {
        let (dd, ds, dds) = (self.angular_diameter_distance(0.0, zd), self.angular_diameter_distance(0.0, zs), self.angular_diameter_distance(zd, zs));
        (1.0 + zd) * dd * ds / dds
    }
}

/// The transverse comoving distance (m) to redshift `z` in a flat ΛCDM universe with Hubble
/// constant `h0` (s⁻¹, e.g. `0.7 * H100`) and matter density `omega_m`, radiation neglected
/// (its `Ω_r ≈ 9×10⁻⁵` moves `D_M` by about 10⁻⁴ at `z = 2`).
pub fn comoving_distance(z: f64, h0: f64, omega_m: f64) -> f64 {
    Cosmology::flat(h0, omega_m).comoving_distance(0.0, z)
}

/// The angular diameter distance (m) of an object at redshift `z2` seen from redshift `z1 < z2`
/// (`z1 = 0` for the observer), in the same universe.
pub fn angular_diameter_distance(z1: f64, z2: f64, h0: f64, omega_m: f64) -> f64 {
    Cosmology::flat(h0, omega_m).angular_diameter_distance(z1, z2)
}

/// The time-delay distance (m) of a lens at redshift `zd` and a source at `zs`:
/// `(1 + z_d)D_d D_s/D_ds`.
pub fn time_delay_distance(zd: f64, zs: f64, h0: f64, omega_m: f64) -> f64 {
    Cosmology::flat(h0, omega_m).time_delay_distance(zd, zs)
}

#[cfg(test)]
mod tests {
    use super::*;
    use metrology::MPC;

    #[test]
    fn the_two_closed_forms() {
        // Einstein–de Sitter (Ω_m = 1): D_M = (2c/H₀)(1 − 1/√(1 + z)); de Sitter (Ω_m = 0): cz/H₀.
        let h0 = 0.7 * H100;
        for z in [0.1, 1.0, 3.0] {
            let eds = 2.0 * C / h0 * (1.0 - 1.0 / (1.0f64 + z).sqrt());
            assert!((comoving_distance(z, h0, 1.0) / eds - 1.0).abs() < 1e-13, "z = {z}");
            assert!((comoving_distance(z, h0, 0.0) / (C * z / h0) - 1.0).abs() < 1e-14, "z = {z}");
        }
        assert_eq!(comoving_distance(0.0, h0, 0.3), 0.0);
    }

    #[test]
    fn h0licow() {
        // Wong et al. (2020), Tables 1–2: the six lenses' (z_d, z_s) and median D_Δt (Mpc) with
        // its 16th–84th percentiles. D_Δt ∝ 1/H₀, so each gives H₀ = 70 D_Δt(70)/D_Δt; their
        // inverse-variance mean at Ω_m = 0.3 is the paper's flat-ΛCDM 73.3 +1.7/−1.8 to 0.1.
        // The paper's own value is the product of the six posteriors, four of them joint in
        // (D_Δt, D_d) through kernel density estimates of unpublished chains; Table 2 prints only
        // the marginals, so that likelihood cannot be rebuilt from the paper.
        let lenses = [
            (0.6304, 1.394, 5156.0, 296.0, 236.0),
            (0.295, 0.654, 2096.0, 98.0, 83.0),
            (0.4546, 1.693, 2707.0, 183.0, 168.0),
            (0.745, 1.789, 5769.0, 589.0, 471.0),
            (0.6575, 1.662, 4784.0, 399.0, 248.0),
            (0.311, 1.722, 1470.0, 137.0, 127.0),
        ];
        let (mut sum, mut weights) = (0.0, 0.0);
        for (zd, zs, ddt, up, down) in lenses {
            let h = 70.0 * time_delay_distance(zd, zs, 0.7 * H100, 0.3) / MPC / ddt;
            let sigma = h * 0.5 * (up + down) / ddt;
            sum += h / (sigma * sigma);
            weights += 1.0 / (sigma * sigma);
        }
        let (h0, sigma) = (sum / weights, weights.powf(-0.5));
        assert!((h0 - 73.3).abs() < 0.1 && (sigma - 1.75).abs() < 0.1, "{h0} ± {sigma}");
    }

    /// `Cosmology { h0: h·H100, omega_m, omega_de }` with no radiation and Λ.
    fn lcdm(h: f64, omega_m: f64, omega_de: f64) -> Cosmology {
        Cosmology { omega_de, ..Cosmology::flat(h * H100, omega_m) }
    }

    #[test]
    fn external_codes_in_astropys_tests() {
        // icosmo at z = 1, h = 0.7 (astropy's cosmo_{flat,open,closed}.ecsv): D_M (Mpc) for
        // (Ω_m, Ω_Λ) = (0.3, 0.7), (0.3, 0.1) and (2, 0.1), printed to 8, 7 and 8 digits.
        for (om, ol, dm) in [(0.3, 0.7, 3303.8288), (0.3, 0.1, 2979.346), (2.0, 0.1, 2152.7954)] {
            let got = lcdm(0.7, om, ol).transverse_distance(0.0, 1.0) / MPC;
            assert!((got / dm - 1.0).abs() < 5e-8, "{om} {ol}: {got}");
        }
        // Kayser, Helbig & Schramm's code: flat (70.4, 0.272), D_A(1, 2) = 646.22968662822018 Mpc.
        let da = lcdm(0.704, 0.272, 0.728).angular_diameter_distance(1.0, 2.0) / MPC;
        assert!((da / 646.229_686_628_220_2 - 1.0).abs() < 1e-14, "{da}");
        // Mathematica's w₀w_a luminosity distances (Ω_m = 0.2, Ω_de = 0.8, w₀ = −1.1, w_a = 0.2),
        // printed to six digits; the CPL density scale at z = 0.5 to every digit.
        let cpl = Cosmology { w0: -1.1, wa: 0.2, ..lcdm(0.7, 0.2, 0.8) };
        for (z, dl) in [(0.2, 1004.0), (0.4, 2268.62), (0.9, 6265.76), (1.2, 9061.84)] {
            assert!((cpl.luminosity_distance(z) / MPC / dl - 1.0).abs() < 3e-6, "z = {z}");
        }
        let scale = (cpl.e(0.5).powi(2) - 0.2 * 1.5f64.powi(3)) / 0.8;
        assert!((scale - 0.924_631_066_952_902_1).abs() < 1e-15, "{scale}");
    }

    #[test]
    fn exact_neutrinos() {
        // ρ_ν/ρ_γ ÷ (0.22710731766 N_eff) for three 0.25 eV neutrinos at T_CMB = 3 K, from the
        // exact integral in Mathematica (astropy's test_massivenu_density), at z = 0, 1, 2, 10, 1000.
        let nu = Cosmology { t_cmb: 3.0, n_eff: 3.0, m_nu: [0.25 * EV; 3], ..lcdm(0.75, 0.25, 0.75) };
        let per = 7.0 / 8.0 * (4.0f64 / 11.0).powf(4.0 / 3.0) * 3.0;
        for (z, ratio) in [(0.0, 429.924), (1.0, 214.964), (2.0, 143.312), (10.0, 39.1005), (1000.0, 1.11086)] {
            assert!((nu.neutrino_ratio(z) / per / ratio - 1.0).abs() < 2e-6, "z = {z}");
        }
        // Massless neutrinos are radiation: (7/8)(4/11)^{4/3} each.
        assert!((Cosmology { m_nu: [0.0; 3], ..nu }.neutrino_ratio(5.0) / per - 1.0).abs() < 1e-15);
    }

    #[test]
    fn planck18() {
        // Flat by construction, Ω_m + Ω_ν = 0.3111 as Planck quotes, and the CMB angle:
        // 100θ* = 100r*/D_M(z*) = 1.04124 against Table 2's 1.04119 ± 0.00029, but 1.03639
        // (16σ off) when the same universe has no radiation.
        let p = PLANCK18;
        assert!(p.curvature().abs() < 1e-15, "{}", p.curvature());
        let omega_nu = p.omega_photons() * p.neutrino_ratio(0.0);
        assert!((p.omega_photons() / 5.402015e-5 - 1.0).abs() < 1e-6 && (omega_nu / 1.439709e-3 - 1.0).abs() < 1e-6);
        assert!((p.omega_m + omega_nu - 0.3111).abs() < 5e-5);
        let theta = |c: Cosmology| 100.0 * 144.57 * MPC / c.transverse_distance(0.0, 1089.80);
        assert!((theta(p) - 1.04124).abs() < 5e-6, "{}", theta(p));
        assert!((theta(Cosmology::flat(p.h0, 0.30966 + 0.0014)) - 1.03639).abs() < 5e-6);
        // D_C(1) = 3395.6346 Mpc, from the research script.
        assert!((p.comoving_distance(0.0, 1.0) / MPC - 3395.6346).abs() < 1e-4);
    }

    #[test]
    fn curvature_and_the_free_functions() {
        // Mattig's closed form (Hogg Eq. 17, Ω_Λ = 0) for open, flat and closed universes.
        for om in [0.3, 1.0, 2.0] {
            let c = lcdm(0.7, om, 0.0);
            for z in [0.5, 2.0, 5.0] {
                let mattig = C / c.h0 * 2.0 * (2.0 - om * (1.0 - z) - (2.0 - om) * (1.0 + om * z).sqrt()) / (om * om * (1.0 + z));
                assert!((c.transverse_distance(0.0, z) / mattig - 1.0).abs() < 1e-13, "{om} {z}");
            }
        }
        // Closed, Ω_m = 3: D_A(2, 50) = 40.5883 Mpc from S_k(χ₂ − χ₁); Hogg's Eq. 19 gives 30.0828.
        let closed = lcdm(0.7, 3.0, 0.0);
        assert!((closed.angular_diameter_distance(2.0, 50.0) / MPC - 40.5883).abs() < 1e-4);
        let (dh, ok) = (C / closed.h0, closed.curvature());
        let (m1, m2) = (closed.transverse_distance(0.0, 2.0), closed.transverse_distance(0.0, 50.0));
        let hogg = (m2 * (1.0 + ok * m1 * m1 / (dh * dh)).sqrt() - m1 * (1.0 + ok * m2 * m2 / (dh * dh)).sqrt()) / 51.0;
        assert!((hogg / MPC - 30.0828).abs() < 1e-4);
        // The flat functions are the radiation-free limit of the struct.
        let flat = lcdm(0.7, 0.3, 0.7);
        assert_eq!(flat.curvature(), 0.0);
        assert_eq!(time_delay_distance(0.295, 0.654, flat.h0, 0.3), flat.time_delay_distance(0.295, 0.654));
        assert_eq!(comoving_distance(1.0, flat.h0, 0.3), flat.transverse_distance(0.0, 1.0));
    }
}
