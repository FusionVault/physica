//! Tests of the weak equivalence principle: the Eötvös parameter η, what the best experiments
//! measured, and how strongly a new long-range force may couple to baryon or lepton number
//! before MICROSCOPE would have seen it.
//!
//! A new force from a boson coupled to a charge `Q` (baryon number `B`, lepton number `L`, or
//! `B − L`) pulls on each body in proportion to its charge per unit mass, which differs between
//! elements, so two bodies of different composition fall differently. MICROSCOPE compared
//! titanium and platinum in orbit for two years and found no difference at 10⁻¹⁵
//! ([Touboul et al. 2022](https://arxiv.org/abs/2209.15487)); in
//! [Fayet's (2025)](https://arxiv.org/abs/2507.02723) reading this bounds a `B − L` coupling at
//! 1.1×10⁻²⁵, 10²⁵ times weaker than electromagnetism.
//!
//! The identifiers say `eta`, since `gravimetry::local::EOTVOS` is the unit of gravity gradient.
//!
//! ```
//! use gravitics::equivalence::{MICROSCOPE_LIMITS, Mediator, TESTS, coupling_bound};
//!
//! // MICROSCOPE's η(Ti, Pt) = (−1.5 ± 2.3 ± 1.5)×10⁻¹⁵: 2.7×10⁻¹⁵ in all.
//! assert!((TESTS[0].combine() - 2.746e-15).abs() < 1e-18);
//! // The Earth's (B − L) per atomic mass unit, and Ti's less Pt's (Fayet 2025, Eqs. 6–7).
//! let g = coupling_bound(Mediator::Vector, 0.5138, -0.05625, MICROSCOPE_LIMITS);
//! assert!((g / 1.067e-25 - 1.0).abs() < 1e-3); // "|g_B−L| < 1.1×10⁻²⁵"
//! ```

use metrology::{ALPHA, AMU, C, G, HBAR};

/// The Eötvös parameter `η = 2(a_A − a_B)/(a_A + a_B)` of two bodies falling with accelerations
/// `a_a` and `a_b` (Touboul et al. 2022, Eq. 1).
pub fn eta(a_a: f64, a_b: f64) -> f64 {
    2.0 * (a_a - a_b) / (a_a + a_b)
}

/// One measurement of η with its statistical and systematic errors (1σ). Where a paper quotes
/// one combined error it is in `stat` and `syst` is zero.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Test {
    /// The experiment and the two materials.
    pub name: &'static str,
    /// η of the first material against the second.
    pub eta: f64,
    /// The statistical error.
    pub stat: f64,
    /// The systematic error.
    pub syst: f64,
    /// Where it comes from.
    pub source: &'static str,
}

impl Test {
    /// The two errors in quadrature.
    pub fn combine(&self) -> f64 {
        self.stat.hypot(self.syst)
    }
}

/// Published tests of the equivalence principle, the tightest first: in orbit (MICROSCOPE's
/// Ti–Pt result and its Pt–Pt control), a torsion balance, and atoms in free fall.
pub const TESTS: [Test; 4] = [
    Test { name: "MICROSCOPE Ti–Pt", eta: -1.5e-15, stat: 2.3e-15, syst: 1.5e-15, source: "https://arxiv.org/abs/2209.15487" },
    Test { name: "MICROSCOPE Pt–Pt (control)", eta: 0.0, stat: 1.1e-15, syst: 2.3e-15, source: "https://arxiv.org/abs/2209.15487" },
    Test { name: "Eöt-Wash Be–Ti, towards the Earth", eta: 0.3e-13, stat: 1.8e-13, syst: 0.0, source: "https://arxiv.org/abs/0712.0607" },
    Test { name: "Stanford ⁸⁵Rb–⁸⁷Rb", eta: 1.6e-12, stat: 1.8e-12, syst: 3.4e-12, source: "https://arxiv.org/abs/2005.11624" },
];

/// MICROSCOPE's one-sided 95% limits on η(Ti, Pt): below 4.5×10⁻¹⁵ if positive, below 6.5×10⁻¹⁵
/// in size if negative. The result leans negative, so a positive η is the more constrained
/// (Fayet 2025, Eq. 16).
pub const MICROSCOPE_LIMITS: (f64, f64) = (4.5e-15, 6.5e-15);

/// `α/(G u²) = e²/(4πε₀ G u²)`, 1.2536×10³⁶: the electric repulsion of two unit charges over the
/// gravitational attraction of two atomic mass units, the scale of Fayet's Eq. 7.
pub const ALPHA_OVER_G_U2: f64 = ALPHA * HBAR * C / (G * AMU * AMU);

/// The spin of the boson that carries a new force.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Mediator {
    /// Spin 1, like the photon: like charges repel.
    Vector,
    /// Spin 0: like charges attract.
    Scalar,
}

/// The η a new long-range force puts between two bodies falling towards the Earth (Fayet 2025,
/// Eq. 7): `η = ∓(α/(G u²))ε²(Q/A_r)_⊕Δ(Q/A_r)`, the upper sign for a vector mediator. `epsilon`
/// is the coupling in units of `e`, `earth` the Earth's charge per atomic mass unit and
/// `difference` the first body's less the second's. For a range much longer than the Earth's
/// radius; a finite range `λ` multiplies η by `(1 + r/λ)e^(−r/λ)` at distance `r`
/// ([`yukawa_factor`](crate::modified::yukawa_factor) less one) and by the Earth's form factor
/// (see [`sphere_form_factor`]) (Fayet's Eq. 94).
pub fn eta_from_coupling(mediator: Mediator, epsilon: f64, earth: f64, difference: f64) -> f64 {
    let sign = match mediator {
        Mediator::Vector => -1.0,
        Mediator::Scalar => 1.0,
    };
    sign * ALPHA_OVER_G_U2 * epsilon * epsilon * earth * difference
}

/// The largest coupling `|g| = |ε|√(4πα)` that one-sided limits `(positive, negative)` on η
/// allow, for the same inputs as [`eta_from_coupling`]: the mediator and the sign of the charge
/// difference fix the sign of η, and so which limit applies. With [`MICROSCOPE_LIMITS`] and a
/// `B − L` charge it is 1.07×10⁻²⁵ for a vector and 1.28×10⁻²⁵ for a scalar, which Fayet rounds
/// to 1.1 and 1.3 (his Eqs. 24–25). A finite range raises it by `F^(−1/2)`, where `F` is the
/// factor that range puts on η.
pub fn coupling_bound(mediator: Mediator, earth: f64, difference: f64, (positive, negative): (f64, f64)) -> f64 {
    let per_epsilon2 = eta_from_coupling(mediator, 1.0, earth, difference);
    let limit = if per_epsilon2 > 0.0 { positive } else { negative };
    (limit / per_epsilon2.abs()).sqrt() * (4.0 * core::f64::consts::PI * ALPHA).sqrt()
}

/// How much more strongly a homogeneous sphere of radius `R` sources a Yukawa force of range `λ`
/// than a point of the same charge, at `x = R/λ`: `φ(x) = 3(x cosh x − sinh x)/x³` (Fayet 2025,
/// Eq. 37). It is 1 for a long range and grows as `(3/2)eˣ/x²`, since the charge near the surface
/// is nearer. A layered body is the mass-weighted sum of its shells' terms (his Eq. 80).
pub fn sphere_form_factor(x: f64) -> f64 {
    if x.abs() < 1.0 {
        // Σ 3x²ᵏ/((2k + 1)!(2k + 3)): the closed form loses every digit to cancellation as x → 0.
        let (mut term, mut sum) = (1.0, 0.0);
        for k in 0..12 {
            sum += 3.0 * term / f64::from(2 * k + 3);
            term *= x * x / f64::from((2 * k + 2) * (2 * k + 3));
        }
        sum
    } else {
        3.0 * (x * x.cosh() - x.sinh()) / x.powi(3)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::modified::{boson_mass_ev, yukawa_factor};

    #[test]
    fn eotvos_parameters() {
        // η is antisymmetric, and Δa/g to first order in Δa.
        let g = 7.95;
        assert!(eta(3.0, 1.0) == 1.0 && eta(1.0, 2.0) == -eta(2.0, 1.0) && eta(g, g) == 0.0);
        assert!((eta(g + 1e-6, g) / (1e-6 / g) - 1.0).abs() < 1e-7);
        // Touboul et al. (2022), Eq. 5 and the control: η(Ti, Pt) = [−1.5 ± 2.3 ± 1.5]×10⁻¹⁵,
        // 2.746×10⁻¹⁵ in quadrature, a differential pull of 2.2×10⁻¹⁴ m/s² at 710 km.
        let [ti_pt, pt_pt, eot_wash, stanford] = TESTS;
        assert!((ti_pt.combine() / 2.745_906e-15 - 1.0).abs() < 1e-6);
        assert!((ti_pt.combine() * metrology::GM_EARTH / (6371e3 + 710e3f64).powi(2) / 2.18e-14 - 1.0).abs() < 0.01);
        assert!(ti_pt.eta.abs() < ti_pt.combine() && pt_pt.eta.abs() < pt_pt.combine());
        // Schlamminger et al. (2008), Eq. 2, and Asenbaum et al. (2020): 3.85×10⁻¹² in all.
        assert!(eot_wash.eta.abs() < eot_wash.combine() && (stanford.combine() / 3.847e-12 - 1.0).abs() < 1e-3);
    }

    #[test]
    fn fayets_coupling_bounds() {
        assert!((ALPHA_OVER_G_U2 / 1.2536e36 - 1.0).abs() < 5e-5);
        // Eq. 6: Δ(Q/A_r) of Ti less Pt. Eq. 7 prints the combined coefficients ∓1.00×10³³ (B),
        // ±3.623×10³⁴ (B − L) and ∓3.482×10³⁴ (L); they give the Earth's charges per mass. Its L
        // and B − L add up to 1.0008 nucleons per u, between ¹⁶O's 1.0003 and ⁵⁶Fe's 1.0012; B's
        // own coefficient and Δ(B/A_r) carry two or three digits, so its 1.010 is good to 1%.
        let (db, dl): (f64, f64) = (0.00079, 0.05704);
        let dbl = db - dl;
        assert!((dbl + 0.05625).abs() < 1e-15);
        let earth = |coefficient: f64, d: f64| coefficient / (ALPHA_OVER_G_U2 * d.abs());
        let (eb, ebl, el) = (earth(1.00e33, db), earth(3.623e34, dbl), earth(3.482e34, dl));
        assert!((ebl + el - 1.0008).abs() < 1e-4 && (eb - (ebl + el)).abs() < 0.01);
        assert!((eta_from_coupling(Mediator::Vector, 1.0, ebl, dbl) / 3.623e34 - 1.0).abs() < 1e-12);
        // Eqs. 24–25 (95% CL), and note 40's unrounded 1.07 and 1.28 for B − L.
        for (mediator, q, d, want) in [
            (Mediator::Vector, ebl, dbl, 1.1e-25), (Mediator::Vector, eb, db, 7.7e-25), (Mediator::Vector, el, dl, 1.3e-25),
            (Mediator::Scalar, ebl, dbl, 1.3e-25), (Mediator::Scalar, eb, db, 6.4e-25), (Mediator::Scalar, el, dl, 1.1e-25),
        ] {
            let got = coupling_bound(mediator, q, d, MICROSCOPE_LIMITS);
            assert!((got / want - 1.0).abs() < 0.05, "{mediator:?}: {got:e} against {want:e}");
        }
        let vector = coupling_bound(Mediator::Vector, ebl, dbl, MICROSCOPE_LIMITS);
        let scalar = coupling_bound(Mediator::Scalar, ebl, dbl, MICROSCOPE_LIMITS);
        assert!((vector / 1.067e-25 - 1.0).abs() < 1e-3 && (scalar / 1.283e-25 - 1.0).abs() < 1e-3);
        // The two differ by √(6.5/4.5), the dissymmetry of the result.
        assert!((scalar / vector - (6.5f64 / 4.5).sqrt()).abs() < 1e-12);
    }

    #[test]
    fn the_form_factor_of_a_sphere() {
        // Fayet's Table IV, homogeneous sphere, to its printed digits.
        for (x, want) in [(0.5, 1.0252), (1.0, 1.1036), (2.0, 1.4616), (3.0, 2.2428), (5.0, 7.1243), (10.0, 297.36), (20.0, 1.728e6), (30.0, 1.722e10), (50.0, 3.049e18), (100.0, 3.992e39), (200.0, 2.696e82)] {
            assert!((sphere_form_factor(x) / want - 1.0).abs() < 2.5e-4, "x = {x}: {}", sphere_form_factor(x));
        }
        // 1 + x²/10 + x⁴/280 + x⁶/15120 near zero, where the closed form has no digits left; the
        // series and the closed form meet at x = 1.
        for x in [0.0, 1e-8, 1e-3, 0.1] {
            let x2: f64 = x * x;
            let want = 1.0 + x2 / 10.0 + x2 * x2 / 280.0 + x2 * x2 * x2 / 15120.0;
            assert!((sphere_form_factor(x) - want).abs() < 1e-8 * x2 * x2 * x2 + 2e-16, "{x}");
        }
        let (below, above) = (sphere_form_factor(1.0 - 1e-15), sphere_form_factor(1.0));
        assert!((below - above).abs() < 4e-15);
        assert!((sphere_form_factor(-2.0) - sphere_form_factor(2.0)).abs() < 1e-15);
    }

    #[test]
    fn a_layered_earth_at_satellite_height() {
        // Fayet's five-shell Earth (Eq. 80): Φ(x) = Σ yᵢ φ(Rᵢx/R), Table IV's Φ(3) = 2.0021.
        let shells = [(0.5, 1.0), (0.1333, 0.9953), (0.1715, 0.8948), (0.1929, 0.5462), (0.0023, 0.1916)];
        let earth = |x: f64| shells.iter().map(|&(y, s)| y * sphere_form_factor(s * x)).sum::<f64>();
        assert!((earth(3.0) / 2.0021 - 1.0).abs() < 1e-5 && (earth(0.0) - 1.0).abs() < 1e-15);
        // Table V: at MICROSCOPE's 710 km, a boson of mass m weakens η by F = (1 + r/λ)e^(−r/λ)Φ,
        // so the coupling limits grow by G = F^(−1/2).
        let (radius, r) = (6371e3, 6371e3 + 710e3);
        let hbar_c = HBAR * C / metrology::EV; // eV m
        for (m, f, g) in [(1e-14, 0.957, 1.022), (2e-14, 0.867, 1.074), (5e-14, 0.574, 1.320), (1e-13, 0.281, 1.886), (2e-13, 0.0927, 3.28)] {
            let lambda = hbar_c / m;
            assert!((boson_mass_ev(lambda) / m - 1.0).abs() < 1e-12);
            let big_f = (yukawa_factor(r, 1.0, lambda) - 1.0) * earth(radius / lambda);
            // To the table's three digits.
            assert!((big_f / f - 1.0).abs() < 2e-3 && (big_f.powf(-0.5) / g - 1.0).abs() < 2e-3, "m = {m:e}: F = {big_f}");
        }
    }
}
