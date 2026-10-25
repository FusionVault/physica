# Gravity experiments, gravimetry and gravitational-wave physics: a formula catalogue for the physica crates (2022–2026)

Every number marked *(computed)* was recomputed by `gravity_frontiers_checks.py` in this folder (standard-library Python, run as `python3 -I gravity_frontiers_checks.py`; all assertions pass). Item IDs (T1, E1, …) are used in the closing inferences. Status labels: **established**, **recent and replicated**, **recent and unreplicated**, **contested**.

Coverage check: I read the "What is inside" tables of `gravimetry`, `gravitational`, `gravitics` and `geopotential`, grepped all crate sources, and skimmed `research/topics/gravity-anomalies/report.md`, `gravitational-waves/report.md` and `gravimetry-frontiers/report.md`. None of the crates model Earth tides, Kerr ringdown, remnant mass or spin, horizon area, tidal deformability, normal gravity or W0. The only tide-related code is the permanent-tide `C̄20` conversion.

## Earth tides for gravimeters (Longman, the tidal potential, Love numbers, δ, the permanent tide, ocean loading)

### Takeaway
None of the crates computes a tide correction for a gravimeter. Already covered: `geopotential::Coefficients::to_tide_system`, which shifts `C̄20` only, with k = 0.3 (Rapp 1989). Longman's 1959 closed form now has a public test vector reproduced to 10⁻¹⁰ mGal. The IERS 2010 Love numbers and Mäkinen's 2021 permanent-tide formulas give every other constant needed. A `gravimetry::tides` module is the highest-value gap in this area: an absolute gravimeter cannot reach µGal without it.

### Cited Findings

**T1. Longman lunisolar gravity tide.** New `gravimetry::tides::longman(lat, lon, h, t) -> (moon, sun, total)`. Links gravimetry, geopotential and metrology. Status: established (foundational, 1959).
- Formula ([Longman 1959, JGR 64, 2351](https://doi.org/10.1029/JZ064i012p02351), as implemented in [LongmanTide](https://github.com/jrleeman/LongmanTide)). These are the vertical tidal accelerations, positive upward, multiplied by the gravimetric factor 1 + h₂ − (3/2)k₂ with Longman's h₂ = 0.612 and k₂ = 0.303:
  - g_m = (μMr/d³)(3cos²θ − 1) + (3/2)(μMr²/d⁴)(5cos³θ − 3cosθ)
  - g_s = (μSr/D³)(3cos²φ − 1)
  - The zenith angles come from the station and mean orbital elements: cosθ = sinλ sinI sin l + cosλ[cos²(I/2)cos(l − χ) + sin²(I/2)cos(l + χ)]. The Sun's cosφ uses the obliquity ω, its longitude l₁ = h + 2e₁sin(h − p₁) and χ₁ = t + h.
  - The distances: 1/d = 1/c + a′e cos(s − p) + a′e² cos 2(s − p) + (15/8)a′me cos(s − 2h + p) + a′m² cos 2(s − h), with a′ = 1/(c(1 − e²)). The station radius is r = a/√(1 + 0.006738 sin²λ) + H.
  - The polynomials in Julian centuries T from 1899-12-31 12:00 give the mean longitudes s, p, h, N and p₁. The constants are cgs: μ = 6.673×10⁻⁸, M = 7.3537×10²⁵ g, S = 1.993×10³³ g, e = 0.05490, m = 0.074804, c = 3.84402×10¹⁰ cm, c₁ = 1.495×10¹³ cm, a = 6.378270×10⁸ cm, i = 0.08979719 rad and ω = 23.452°. Longman counts longitude positive west.
- Test vector ([LongmanTide `tests/test_basics.py`](https://github.com/jrleeman/LongmanTide/blob/master/tests/test_basics.py)): latitude 40.7914°, longitude 282.1414° E, 370 m, 2015-04-23 00:00 UT gives moon 0.0324029651226, sun −0.0288682178454 and total 0.00353474727722 mGal. My independent re-implementation reproduces all three to better than 10⁻⁹ mGal *(computed)*.
- Use it as the fast, dependency-free correction (accuracy of order a few µGal). T2 and T3 are the precise route.

**T2. Body-tide gravimetric factors from Love numbers.** New `gravimetry::tides::{delta, gamma}`, or `geopotential::love`. Links gravimetry and geopotential. Status: established.
- δₙ = 1 + (2/n)hₙ − ((n+1)/n)kₙ (the gravity factor), and γₙ = 1 + kₙ − hₙ (tilt and levelling, the "diminishing factor").
- IERS 2010 nominal values: k₂₀ = 0.30190, k₂₁ = 0.29830 − 0.00144i and k₂₂ = 0.30102 − 0.00130i (anelastic), with k⁺ = −0.00089, −0.00080 and −0.00057; k₃ₘ = 0.093 (0.094 for m = 3) ([IERS 2010 Table 6.3 as transcribed in Orekit](https://github.com/CS-SI/Orekit/blob/develop/src/main/resources/assets/org/orekit/IERS-conventions/2010/tab6.3.txt)). The displacement numbers are h(0) = 0.6078, h(2) = −0.0006, h₃ = 0.292, l(0) = 0.0847, l(2) = 0.0002 and l₃ = 0.015 ([IERS 2010 ch. 7 constants in Orekit `IERSConventions.java`](https://github.com/CS-SI/Orekit/blob/develop/src/main/java/org/orekit/utils/IERSConventions.java); [IERS TN36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html)).
- Tests: δ₂(0.6078, 0.30190) = 1.15495; δ₂(Longman's 0.612, 0.303) = 1.1575; δ₃(0.292, 0.093) = 1.0707; γ₂ = 0.6941 *(computed)*. Mäkinen gives 1 + k − h ≈ 0.69 ([Mäkinen 2021, Eq. 52](https://arxiv.org/abs/2006.08440)) and uses δ = 1.16 as the conventional factor ([same, App. A](https://doi.org/10.1007/s00190-021-01541-5)).
- Frequency dependence: the IERS tables of δk_f, for example O₁ δk_fR = −83×10⁻⁵ and M₂ +0.00004, are in [Orekit's tab6.5a–c transcriptions](https://github.com/CS-SI/Orekit/tree/develop/src/main/resources/assets/org/orekit/IERS-conventions/2010). Orekit also documents two corrections to the printed IERS tables: the P₁ out-of-phase sign, and the K₁ value of −0.80 instead of −0.78 mm ([Orekit tab7.3a header](https://github.com/CS-SI/Orekit/blob/develop/src/main/resources/assets/org/orekit/IERS-conventions/2010/tab7.3a.txt)). Any Rust port should carry those corrections.

**T3. Rigid-Earth scale of the tide, and its link to wave detectors.** Doc test in `physica::analogies`. Links gravimetry and gravitational. Status: established.
- The peak vertical lunar tide is 2GM_☾R/d³ = 109.99 µGal, the solar one 50.51 µGal, and their ratio 0.459 *(computed, mean distances)*.
- The lunar tidal gradient at Earth is 2GM_☾/d³ = 1.73×10⁻⁴ E *(computed)*. That is the same tensor as `gravitational::detectors::tidal_acceleration` (Δa = ½ω²hL). The analogy is one line: a tide is a Newtonian "strain" field.

**T4. The permanent tide in gravity, geoid and heights.** New `geopotential::tide` (station quantities), complementing `to_tide_system` (coefficients). Links geopotential and gravimetry. Status: established (formalised 2021).
- Conventional permanent tidal potential: W_T(r, φ) = A(r/a)²(sin²φ − 1/3), with A = −2.9166 m²s⁻². Equivalently A′ = 2A/3 = −1.9444 m²s⁻² on P₂, and A″ = A′/√5 = −0.86956 on the normalised P̄₂. The ellipsoidal form is W_T(ϕ, h) = (1 + 2h/a)(0.9722 − 2.8841 sin²ϕ − 0.0195 sin⁴ϕ) m²s⁻² ([Mäkinen 2021, J. Geod. 95:106, Eqs. 18–19](https://doi.org/10.1007/s00190-021-01541-5); [arXiv:2006.08440](https://arxiv.org/abs/2006.08440)).
- Its contribution to free-fall acceleration is g_T(ϕ) = −30.49 + 90.95 sin²ϕ + 0.31 sin⁴ϕ µGal (Eq. 20). In the spherical approximation, −2A/a·(sin²φ − 1/3) gives −30.485 µGal at the equator and +60.97 µGal at the pole, against 60.77 published; the ellipsoidal terms account for the difference *(computed)*.
- The mean-tide geoid lies above the zero-tide geoid by −0.296(sin²φ − 1/3) m (Heikkinen 1978, Eq. A-4), which is 9.9 cm at the equator; A/γ_e gives 0.0994 m *(computed)*. Ekman's tide-free to zero-tide correction to a geopotential-number height is −(−h₂ + k₂)(99 − 296 sin²ϕ) mm, with h₂ = 0.605 and k₂ = 0.30190 (Eq. A-8). That is 30.0 mm at the equator *(computed)*.
- The Honkasalo term c_ls = 0.037(1 − 3 sin²ϕ) mGal with δ = 1.20 (Eq. A-1) is the legacy correction that makes IGSN71 gravity mean-tide ([Mäkinen 2021, App. A](https://arxiv.org/abs/2006.08440)).
- The IHRF recommendation is to model in zero-tide and add mean-tide afterwards ([Mäkinen 2021 abstract](https://arxiv.org/abs/2006.08440)).

### Inferences
- A `gravimetry::tides` module of T1, T2 and T4 joins gravimetry and geopotential directly: absolute-gravity values need the zero-tide or mean-tide label that `geopotential::TideSystem` already carries for coefficients. Reusing that enum for station gravity and heights gives one tide-system vocabulary across both crates.
- `metrology` lacks GM_MOON. T1 and T3 both need it; add it beside `GM_SUN`. Longman's own cgs constants should stay inside `longman()` so the test vector holds exactly.
- Implementation order: Longman first (exact test vector), then δ and γ, then the permanent tide. A harmonic tidal-potential catalogue (HW95, KSM03) is a larger job; see Gaps.

### Gaps
- **Tidal-potential catalogues.** Unverified this session: HW95 (Hartmann & Wenzel 1995, GRL, doi:10.1029/95GL03324; about 12 935 waves) and KSM03 (Kudryavtsev 2004, J. Geod.; about 28 806 terms). I could not fetch their wave tables or a station amplitude table to test against, and IERS TN36 ch. 6–7 was unreachable (obspm.fr refused the connection and iers.org returned 404). The values above come from Orekit's transcription of the IERS tables, not from the PDF.
- **Ocean-tide loading.** Not fetched: Farrell 1972 (Rev. Geophys., doi:10.1029/RG010i003p00761) Green's-function tables and the Bos–Scherneck loading service. The Newtonian part of the gravity loading Green's function has the closed form G·m/(4R² sin(ψ/2)), from Σ n Pₙ(cos ψ) = −1/(4 sin(ψ/2)) (my derivation). The elastic part needs load Love numbers h′ₙ, which `gravimetry::love_numbers` lacks (it holds only k′ₙ from Wahr et al. 1998). No published station value was found to test against.
- **Body-tide models.** No 2022–2026 paper changing the δ factors was found. DDW99 latitude-dependent δ values are not verified here.

## Small-mass and tabletop gravity, including G (2021–2026)

### Takeaway
Much of this area is already in the repository. Already covered: `gravitics::anomalies::{CODATA_2022_G, gls_mean, birge_ratio}`; `gravimetry::atom::aharonov_bohm_rate` (Overstreet 2022); `gravitics::quantum::bmv_phases`; and the doc-test values for Westphal 2021, Rosi 2014 and Fuchs 2024 in `gravity-anomalies/report.md`. I found no peer-reviewed source smaller than Westphal's 92 mg sphere, and no peer-reviewed 2024–2025 G value, so this question adds little new code.

### Cited Findings
- **Levitated test masses (Fuchs et al. 2024, Sci. Adv.).** A 0.43 mg magnetically levitated particle sensed a 2.4 kg source about half a metre away, at 26.7 Hz, with Q above 10, force noise of 0.5 fN/√Hz, and a quoted coupling of 1030 aN ([PMC10889343](https://pmc.ncbi.nlm.nih.gov/articles/PMC10889343/); [arXiv:2303.03545](https://arxiv.org/abs/2303.03545)). Status: recent and unreplicated.
  - Point masses: Gm·M/r² at 0.5 m gives about 2.8×10⁻¹⁶ N (280 aN). Reaching 1030 aN needs r ≈ 0.26 m *(computed)*. The repository report already flags a 1030 vs 30 aN inconsistency inside the paper ([gravity-anomalies report](../../gravity-anomalies/report.md)). Do not use this value as a test.
- **G.** CODATA 2022 is already covered. A torsion-pendulum group reports two experiments agreeing at 12 ppm and aims for 2 ppm; these are targets, not results ([APS April 2024 abstract](https://meetings-archive.aps.org/apr/2024/d19/4)). A 2026 NIST result is reported only second-hand ([Starts With A Bang](https://startswithabang.substack.com/p/a-new-experiment-deepens-the-physics)). The repo report records 6.67366(20) vs 6.67387, unresolved. Status: contested (the G scatter).
- **Atom-interferometric G.** Rosi et al. 2014 (G = 6.67191(99)×10⁻¹¹) is already a doc test. `gravimetry::atom::phase_sensitivity(k_eff, T)` returns the single-diffraction k_eff·T², which gives 4.12×10⁵ rad per m/s² for Rosi's T = 0.16 s. This is already covered and documented in `atom.rs`.

### Inferences
- The only new code worth adding here is a test of `phase_sensitivity` against Rosi's 0.547870 rad phase. It is already listed as a doc test in the gravity-anomalies report, so nothing new is proposed.
- The Planck-mass source (21.76 µg) remains the next milestone. No result at or below 1 mg as a *source* was found up to October 2026.

### Gaps
- I found no 2023–2026 paper measuring gravity from a source smaller than 92 mg. One extended search returned only Fuchs 2024.
- The NIST 2026 G value and uncertainty were not obtained from the primary source.

## Equivalence-principle tests (MICROSCOPE, atom interferometers, torsion balances) and an Eötvös calculator

### Takeaway
No crate has an Eötvös-parameter calculator; MICROSCOPE is only a proposed doc test. The final MICROSCOPE η, together with Fayet's 2025 charge-per-mass formula, turns into a short, fully testable `gravitics::equivalence` module. It reproduces the published B−L coupling bounds of 1.1×10⁻²⁵ and 1.3×10⁻²⁵ to 3%, and it reuses the existing `gravitics::modified::yukawa_factor`.

### Cited Findings

**P1. Eötvös parameter and its error budget.** New `gravitics::equivalence::{eotvos, combine}`. Links gravitics and metrology. Status: established.
- η = 2(a_A − a_B)/(a_A + a_B) ≈ Δa/g.
- Results to test against:
  - MICROSCOPE final: η(Ti, Pt) = [−1.5 ± 2.3 (stat) ± 1.5 (syst)]×10⁻¹⁵ ([Touboul et al. 2022, PRL 129, 121102, arXiv:2209.15487](https://arxiv.org/abs/2209.15487); [CQG 39, 204009, doi:10.1088/1361-6382/ac84be](https://arxiv.org/abs/2209.15487v1)). Its quadrature is 2.746×10⁻¹⁵, which at 710 km altitude (g = 7.95 m/s²) is a differential acceleration of 2.18×10⁻¹⁴ m/s² *(computed)*.
  - Asenbaum et al. 2020, ⁸⁵Rb–⁸⁷Rb: η = [1.6 ± 1.8 (stat) ± 3.4 (sys)]×10⁻¹², with 2 s of free fall and a sensitivity of 5.4×10⁻¹¹/√Hz ([PRL 125, 191101, arXiv:2005.11624](https://arxiv.org/abs/2005.11624)). Its quadrature is 3.85×10⁻¹² *(computed)*.
  - Albers et al. 2020, Rb–K: η = −1.9×10⁻⁷ with σ = 3.2×10⁻⁷ ([arXiv:2003.00939](https://arxiv.org/abs/2003.00939)).
- Status: MICROSCOPE established (recent, unreplicated at 10⁻¹⁵; space follow-ups are proposals). The atom results are recent and replicated at the 10⁻⁷–10⁻¹² levels.

**P2. From η to a new-force coupling (MICROSCOPE final).** New `gravitics::equivalence::coupling_bound`. Links gravitics (`modified::yukawa_factor`, `boson_mass_ev`) and metrology. Status: recent and unreplicated (2025 reanalysis of 2022 data).
- δ_Ti-Pt = ∓(α/(G_N u²)) ε² (Q/A_r)_⊕ Δ(Q/A_r)_Ti-Pt, with g = εe = ε√(4πα); the upper sign is for a spin-1 mediator and the lower for spin-0 ([Fayet 2025, PRD, Eq. 7, arXiv:2507.02723](https://arxiv.org/abs/2507.02723); [PRD 10.1103/17wc-xxjm](https://journals.aps.org/prd/pdf/10.1103/17wc-xxjm)).
- Inputs: α/(G_N u²) ≃ 1.2536×10³⁶; Δ(B/A_r) = 0.00079, Δ(L/A_r) = 0.05704, Δ((B−L)/A_r) = −0.05625 (Eq. 6); combined coefficients ∓1.00×10³³ε_B², ±3.623×10³⁴ε²_{B−L} and ∓3.482×10³⁴ε_L².
- Bounds at 95% CL: δ < 4.5×10⁻¹⁵ for δ > 0 and |δ| < 6.5×10⁻¹⁵ for δ < 0 (Eq. 16). Spin-1: |g_{B−L}| < 1.1×10⁻²⁵, |g_B| < 7.7×10⁻²⁵, |g_L| < 1.3×10⁻²⁵. Spin-0: |g_{B−L}| < 1.3×10⁻²⁵, |g_B| < 6.4×10⁻²⁵, |g_L| < 1.1×10⁻²⁵. These hold for m below about 10⁻¹⁴ eV/c² ([Fayet 2025](https://arxiv.org/html/2507.02723v2)).
- Tests *(computed)*:
  - e²/(4πε₀Gu²) = 1.2536×10³⁶.
  - 0.00079 − 0.05704 = −0.05625.
  - Spin-1: √(4.5×10⁻¹⁵/3.623×10³⁴)·√(4πα) = 1.067×10⁻²⁵, against 1.1 published.
  - Spin-0: with 6.5×10⁻¹⁵ the bound is 1.283×10⁻²⁵, against 1.3.
  - The implied Earth (B−L)/A_r is 0.514.
- Finite range: the homogeneous-sphere form factor is φ(x) = 3(x cosh x − sinh x)/x³ with x = R/λ (Eq. 37). The general Φ(x) = ⟨sinh(mr)/(mr)⟩ is weighted by ρ(r) (Eq. 32), and at the satellite there is an extra e^(−mz) with z ≃ 710 km (Eq. 55). This links to `geopotential` if an Earth density profile is ever added.

### Inferences
- P1 and P2 make a compact, fully testable module that joins gravitics (Yukawa range and boson mass already exist) to metrology (u, α, G). It also answers the "Eötvös calculator" request with real published limits rather than toy numbers.
- `gravimetry::local::EOTVOS` (the gradient unit, 10⁻⁹ s⁻²) and the Eötvös *parameter* η share a name. The new module should use `eta` in identifiers to avoid confusion.

### Gaps
- Torsion-balance EP results from 2022–2026 were not found. The latest in the repo is Eöt-Wash 2008, η(Be−Ti) = (0.3 ± 1.8)×10⁻¹³.
- No 2023–2026 atom-interferometer EP result tighter than Asenbaum 2020 was found in one search.
- The 2022 CQG companion paper's own B−L numbers were not read. A search summary also gave a garbled "2.2×10⁻¹³" figure for the 2018 bound, so do not use it.

## Atom-interferometer gravimetry and gradiometry 2022–2026

### Takeaway
Already covered: the space-sensor noise model `gravimetry::atom::AtomInterferometer` (CARIOQA, double diffraction), `gravitational::detectors::AtomGwDetector` (AION-km), `gravitational::ultralight` (MAGIS mid-band) and `gravimetry::local::{sphere, cylinder, averaging_time}`. The 2022 field gradiometer of Stray et al. supplies a clean ground test: 150 E gives 17.5 mrad, and the projection noise is about 44 E/√Hz. It also shows that `AtomInterferometer` hard-codes double diffraction, so it cannot represent ground sensors without a factor-2 switch.

### Cited Findings

**A1. Ground gradiometer phase and projection noise (single diffraction).** Extend `gravimetry::atom` with a diffraction-order field or a ground constructor. Links gravimetry and gravitational (detectors). Status: recent and replicated (outdoor 2022; the method dates from the 1990s).
- Δφ = k_eff·ΔΓ·L·T²; the per-shot gradient noise is √2/(C√N·k_eff·L·T²), and the ASD is that times √(T_cycle).
- Stray et al. 2022 ([Nature 602, 590, doi:10.1038/s41586-021-04315-3](https://pmc.ncbi.nlm.nih.gov/articles/8866129/)):
  - Instrument: 1 m baseline, T = 85 ms, about 10⁵ atoms, 0.7 Hz.
  - Noise: quantum-projection limit about 44 E/√Hz; achieved 466 ± 8 E/√Hz outdoors, reaching 20 E in 10 min.
  - Target: a 2 m × 2 m tunnel, centre (1.89 −0.59/+2.3) m deep, detected at SNR 8 with the centre located to ±0.19 m. The modelled peak was about 150 E, a 17.5 mrad phase.
  - Calibration: a test mass gave 205 ± 13.1 E against 202 E modelled.
- Tests *(computed)*:
  - k_eff·150 E·1 m·(0.085 s)² = 17.45 mrad, against 17.5 published.
  - With C = 1, √2/√10⁵/(k_eff·0.085²)/√0.7 = 45.9 E/√Hz, against about 44.
  - `averaging_time(466, 20)` = 543 s, consistent with "20 E within 10 min".
  - A line-mass (2D) model of the void (Δρ = −1800 kg/m³, 4 m²) gives a top-minus-bottom difference over 1 m of 176 E with the lower sensor at ground level, or 149 E with it 0.2 m up. The paper's sensor height was not extracted, so treat this as a consistency check only.
- Already-covered mismatch: `AtomInterferometer::phase` and `gradient_noise` use 2k_eff (microgravity double diffraction). With Stray's inputs they would give half the noise (about 23 E/√Hz).

**A2. Commercial absolute quantum gravimeters.** Doc data only. Status: recent and replicated.
- The AQG B01 field version repeated better than 50 nm s⁻² against an FG5, with no significant drift over weeks ([Cooke et al. 2021, Geosci. Instrum. Methods Data Syst.](https://doaj.org/article/4b862c019a904fb4a5c66eb53a408450)).
- A 2026 systems review reports an AQG-A01 campaign white noise of about 750 nm s⁻² Hz⁻½ and a one-day standard deviation of 9.4 nm s⁻² ([arXiv:2608.21436](https://arxiv.org/pdf/2608.21436)).
- Test: `averaging_time(750e-9, 9.4e-9)` gives 6.4×10³ s, well inside a day, so a one-day σ of 9.4 nm s⁻² implies non-white noise or a duty cycle below 100%. This is my inference, not the paper's.

**A3. Long-baseline projects.**
- MAGIS-100: Fermilab completed the laser laboratory for the 100 m strontium interferometer in a former access shaft; construction began in 2023 ([Newswise/Fermilab](https://www.newswise.com/articles/fermilab-completes-laser-lab-construction-for-world-s-largest-vertical-atom-interferometer); [design, arXiv:2104.02835](https://arxiv.org/abs/2104.02835)). A secondary source says commissioning is expected by 2028.
- Already covered: `AtomGwDetector` (AION-km, 2.08×10⁻²¹ peak strain) and the `ultralight::MIDBAND` mass window. Status: in construction.

**A4. Gravity-gradient phase of an atom gravimeter, and its compensation.** Candidate `gravimetry::atom::gradient_bias`. Status: established; foundational (Peters et al. 2001, Metrologia 38, 25; Roura 2017 PRL 118, 160401), not re-fetched.
- Δφ = k_eff T²[g − Γ(z₀ + v₀T − (7/12)gT²)]. The bias vanishes when the π-pulse k_eff is shifted by Δk/k = ΓT²/2.
- Internal test only: for Γ = 3.086×10⁻⁶ s⁻² and T = 0.16 s, the (7/12)gT²Γ term is 4.5×10⁻⁷ m/s² (45 µGal) *(my arithmetic)*.

### Inferences
- A1 is a small change to existing code (add `diffraction: f64` with default 2.0, or a `ground()` constructor), and it gives `gravimetry::atom` a real field test. The tunnel case also exercises `gravimetry::local::cylinder` against a published survey: the first gravimetry-crate test against a real underground target.
- Gradiometer, Earth tide and GW detector are the same tidal tensor (T3, A1 and `detectors::tidal_acceleration`). The lunar gradient of 1.7×10⁻⁴ E is five orders below Stray's 20 E resolution, while the 110 µGal tide dominates absolute gravimeters. A one-paragraph `physica::analogies` doc test can show both.

### Gaps
- No 2023–2026 Exail-published field accuracy figure for the AQG-B series was found; the IUGG 2023 AQG-B10 abstract had no numbers in the excerpt.
- AION-10 (Oxford) status was not found.
- The CARIOQA-PMP Phase A results were not re-checked; CARIOQA is already covered in `gravimetry::atom`.
- A4 needs its sources fetched before coding.

## Gravitational waves 2022–2026: ringdown, remnants, the area law, tidal deformability, O4

### Takeaway
The `gravitational` crate stops at the inspiral. Already covered: chirp mass, Schwarzschild ISCO, SNR, peak luminosity, memory, propagation, and PTA backgrounds in `gravitational::sgwb` (Hellings–Downs, NANOGrav 2023). Four verified, testable additions complete the merger and ringdown:
- **Kerr QNM fits** (Berti, Cardoso & Will 2006): they reproduce GW150914's 251 Hz and 4.0 ms.
- **Remnant spin** (Hofmann, Barausse & Rezzolla 2016): it gives exactly 0.68646 for q = 1 and 0.686 for GW250114.
- **Remnant mass** (Barausse, Morozova & Rezzolla 2012): 62.6 M☉ for GW250114, against 62.7 published.
- **Kerr horizon areas**: the GW250114 area test, A_f/A_i = 1.57.

### Cited Findings

**G1. Kerr quasinormal-mode fits.** New `gravitational::ringdown::{qnm, frequency, damping_time}`. Links gravitational and metrology (`SOLAR_MASS`, G, C). Status: established (2006 fits to Leaver-method data).
- Fits: M·ω_R = F = f₁ + f₂(1 − j)^f₃ and Q = q₁ + q₂(1 − j)^q₃; then ω_I = F/(2Q), f = F·c³/(2πGM) and τ = Q/(πf) ([Berti, Cardoso & Will 2006, PRD 73, 064030, Table VIII, arXiv:gr-qc/0512160](https://arxiv.org/abs/gr-qc/0512160)).
- Coefficients (f₁, f₂, f₃, q₁, q₂, q₃; maximum error over j ∈ [0, 0.99]):
  - (2,2,0): 1.5251, −1.1568, 0.1292, 0.7000, 1.4187, −0.4990 (1.85%, 0.88%)
  - (2,2,1): 1.3673, −1.0260, 0.1628, 0.1000, 0.5436, −0.4731 (1.56%, 1.69%)
  - (2,2,2): 1.3223, −1.0257, 0.1860, −0.1000, 0.4206, −0.4256
  - (2,1,0): 0.6000, −0.2339, 0.4175, −0.3000, 2.3561, −0.2277
  - (3,3,0): 1.8956, −1.3043, 0.1818, 0.9000, 2.3430, −0.4810 (1.36%, 0.42%)
  - (4,4,0): 2.3000, −1.5056, 0.2244, 1.1929, 3.1191, −0.4825
  - Table IX (l = 3) and Table X (l = 4) have all m. They were extracted from the PDF with pdftotext and match the paper's typeset tables.
- Exact values to test the fit against (Table II, Mω as real, imaginary):
  - j = 0: (2,2,0) is 0.3737, 0.0890 and (2,2,1) is 0.3467, 0.2739.
  - j = 0.7: (2,2,0) is 0.5326, 0.0808 and (2,2,1) is 0.5212, 0.2442.
  - j = 0.9: (2,2,0) is 0.6716, 0.0649.
- The fits agree within their stated maxima; for example, at j = 0 the fit gives 0.3683 against the exact 0.3737, an error of 1.45% *(computed)*. Tests should use the fit's tolerance, or ship Table II itself and interpolate.
- Event test: GW150914's detector-frame remnant of M_f = 68 M☉ and a_f = 0.67 predicts f₂₂₀ = 251 ± 8 Hz and τ₂₂₀ = 4.0 ± 0.3 ms ([LVC, PRL 116, 221101, arXiv:1602.03841](https://arxiv.org/abs/1602.03841)). The fit gives 248.4 Hz and 4.06 ms, inside both intervals *(computed)*. For GW250114 (62.7 M☉ × 1.09, χ = 0.68) the fit gives 249.0 Hz and 4.10 ms, with no published number to compare *(computed)*.

**G2. Kerr ISCO (needed by G3 and G4, and to generalise `inspiral::f_isco`).** Extend `gravitational::inspiral`. Status: established (Bardeen, Press & Teukolsky 1972, as printed in HBR16 Eqs. 2–6).
- r_ISCO = 3 + Z₂ − sign(a)√((3 − Z₁)(3 + Z₁ + 2Z₂)), with Z₁ = 1 + (1 − a²)^(1/3)[(1 + a)^(1/3) + (1 − a)^(1/3)] and Z₂ = √(3a² + Z₁²).
- E_ISCO = √(1 − 2/(3r_ISCO)) and L_ISCO = (2/(3√3))[1 + 2√(3r_ISCO − 2)] ([HBR16](https://arxiv.org/abs/1605.01938)).
- Tests: a = 0 gives r = 6, E = √(8/9) and L = 2√3. Already covered: `inspiral::f_isco` handles the Schwarzschild case only.

**G3. Remnant spin.** New `gravitational::remnant::final_spin`. Status: established; NR-calibrated, 619 simulations, rms 0.002 for aligned spins ([Hofmann, Barausse & Rezzolla 2016, ApJL 825, L19, arXiv:1605.01938](https://arxiv.org/abs/1605.01938)).
- Aligned spins (Eq. 7): a_fin = a_tot + ν[L_ISCO(a_eff) − 2a_tot(E_ISCO(a_eff) − 1)] + Σ_{i=0}^{3} Σ_{j=0}^{4} k_ij ν^(2+i) a_eff^j.
- Definitions: q = M₂/M₁ ≤ 1, ν = q/(1+q)², a_tot = (a₁ + a₂q²)/(1+q)² and a_eff = a_tot + ξν(a₁ + a₂).
- Coefficients for n_M = 3, n_J = 4 (Table 1):
  - k01 3.39221, k02 4.48865, k03 −5.77101, k04 −13.0459
  - k10 35.1278, k11 −72.9336, k12 −86.0036, k13 93.7371, k14 200.975
  - k20 −146.822, k21 387.184, k22 447.009, k23 −467.383, k24 −884.339
  - k30 223.911, k31 −648.502, k32 −697.177, k33 753.738, k34 1166.89
  - ξ = 0.474046
- k00 is fixed by Eq. 11, √3/2 + Σᵢ k_i0/4^(2+i) = 0.68646 (Scheel et al. 2009's NR value). That gives k00 = −5.97723 *(computed)*, matching the paper's footnote "k00 ≈ −5.9". The paper prefers the 6-coefficient n_M = 1, n_J = 2 set (k01 −1.2019, k02 −1.20764, k10 3.79245, k11 1.18385, k12 4.90494, ξ 0.41616) to avoid overfitting.
- Precessing spins (Eqs. 13–17): |a_fin| = (1+q)⁻²[|a₁|² + |a₂|²q⁴ + 2|a₁||a₂|q² cos α + 2(|a₁| cos β + |a₂| q² cos γ)|ℓ|q + |ℓ|²q²]^(1/2), with ℓ built from Eq. 13 using ν^(1+i).
- Tests *(computed)*: q = 1 without spin gives 0.68646 (by construction). GW250114 (33.6 + 32.2 M☉, no spin) gives 0.6862, against the published 0.68 ± 0.01. The test-particle limit (q = 10⁻⁴, a₁ = 0.9) gives 0.90006.

**G4. Remnant mass (radiated energy).** New `gravitational::remnant::radiated_energy`. Links `energy`. Status: established ([Barausse, Morozova & Rezzolla 2012, ApJ 758, 63, Eqs. 12, 16–18, arXiv:1206.3803](https://arxiv.org/abs/1206.3803)).
- E_rad/M = [1 − Ẽ_ISCO(ã)]ν + 4ν²[4p₀ + 16p₁ã(ã + 1) + Ẽ_ISCO(ã) − 1], with ã = (|a₁| cos β + q²|a₂| cos γ)/(1 + q)², p₀ = 0.04827 ± 0.00039 and p₁ = 0.01707 ± 0.00032. M_f = M(1 − E_rad/M).
- Tests *(computed)*: q = 1 without spin gives E_rad = p₀ exactly, so M_f/M = 0.95173. GW250114 gives M_f = 62.63 M☉, against 62.7 (+1.0 −1.1) published ([Wikipedia: GW250114, citing the PRL supplement](https://en.wikipedia.org/wiki/GW250114)).

**G5. Kerr horizon area and the area-law test.** New `gravitational::remnant::horizon_area`; links gravitics (`spacetime` horizon temperatures and lifetimes via S = A/(4ℓ_P²)). Status: GW250114 recent and unreplicated as a single event; the theorem itself is established.
- A = 8π(GM/c²)²(1 + √(1 − χ²)).
- Tests *(computed)*: GW250114 has A_i = 2.374×10⁵ km² (non-spinning progenitors) and A_f = 3.734×10⁵ km² (source frame), so A_f/A_i = 1.573. A press report gives "about 240 000" and "about 400 000 km²" ([Ars Technica via Wikipedia](https://en.wikipedia.org/wiki/GW250114)); the first matches and the second is loose rounding. For equal masses without spin, A_f/A_i = 1.564.
- Published: network SNR 80; m₁ = 33.6 (+1.2 −0.8) and m₂ = 32.2 (+0.8 −1.3) M☉; χ ≤ 0.26; e ≤ 0.03. The remnant area exceeds the initial area even with up to five merger cycles excised; ringdown mode frequencies are within ±30% of Kerr ([LVK, PRL 135, 111403, arXiv:2509.08054, doi:10.1103/kw5g-d732](https://arxiv.org/abs/2509.08054)).
- The area increase holds at 4.4σ for a start time of −40 t_M, at 3.4σ or more for all start times after −250 t_M, and above 5σ for t ≥ −10 t_M ([arXiv:2509.08054 PDF, via search excerpt](https://arxiv.org/pdf/2509.08054)).
- The companion spectroscopy paper reports the fundamental ℓ = m = 4 mode constrained to tens of percent and the quadrupole frequency to a few percent ([arXiv:2509.08099](https://arxiv.org/abs/2509.08099)). The first-overtone significance is 4.1σ ([Wikipedia: GW250114](https://en.wikipedia.org/wiki/GW250114)).

**G6. Tidal deformability.** New `gravitational::tides::{lambda, lambda_tilde, tidal_phase}`. Links gravitational and gravimetry (the same Love-number concept as T2: k₂ of a neutron star against k₂ of the Earth). Status: established (formula); GW170817 values recent and unreplicated (one event).
- Λ = (2/3)k₂(c²R/(GM))⁵ and Λ̃ = (16/13)[(m₁ + 12m₂)m₁⁴Λ₁ + (m₂ + 12m₁)m₂⁴Λ₂]/(m₁ + m₂)⁵. Equal masses give Λ̃ = Λ exactly.
- GW170817: Λ̃ = 300 (+420 −230) as a 90% highest-posterior-density interval with low-spin priors, giving the upper edge Λ̃ ≤ 720; with high-spin priors Λ̃ ∈ (0, 630) ([LVC, PRX 9, 011001](https://doi.org/10.1103/PhysRevX.9.011001); [Strathprints copy](https://strathprints.strath.ac.uk/id/eprint/68545)). Λ₁.₄ = 190 (+390 −120) is cited by later work as LVC's EOS-insensitive value ([arXiv:2505.14822](https://arxiv.org/pdf/2505.14822)).
- Tests *(computed)*:
  - The Newtonian n = 1 polytrope has k₂ = (15 − π²)/(2π²) = 0.25991, a textbook value.
  - 1.4 M☉ at 11.9 km has compactness 0.1737, so Λ = 190 needs k₂ = 0.0451.
  - Λ̃(1.46, 1.27; 200, 450) = 301.7.

**G7. O4 catalogue (data for tests and docs, not formulas).** Status: recent.
- GWTC-4.0 adds 128 O4a candidates (p_astro ≥ 0.5) for 218 in total, with source properties for 86 at a FAR below 1/yr. GW230518 and GW230529 are NSBH candidates. GW231123_135430 is the heaviest, at a median 137 M☉ component; the lightest is GW230627 at 5.79 M☉ ([arXiv:2508.18082](https://arxiv.org/abs/2508.18082); ApJL 1004, L22).

### Inferences
- G1–G5 form one `gravitational::ringdown` and `remnant` layer with a closed chain of tests. Progenitor masses and spins (G3, G4) give the remnant; the remnant gives the QNMs (G1), which are checked against GW150914; the remnant and progenitors give the area test (G5), checked against GW250114. Every link reuses `metrology` constants.
- G5 joins gravitational and gravitics: `gravitics::spacetime` already has Hawking temperatures and lifetimes, which are functions of the same M and χ.
- G6 joins gravitational and gravimetry or geopotential: the Love number k₂ appears as the Earth's 0.30190 (T2), the Moon's k₃ = 0.0163 (in the gravimetry-frontiers report) and a neutron star's k₂ of about 0.05–0.1. A `physica::analogies` doc test can state Λ for the Earth: k₂ = 0.3 and C = GM/(Rc²) = 6.96×10⁻¹⁰ give Λ ≈ 10⁴⁵, which shows why only compact stars have measurable Λ (my arithmetic, not run in the script).
- Pulsar timing is already covered in `gravitational::sgwb`; nothing new is proposed.

### Gaps
- GW250114's measured f₂₂₀ and τ₂₂₀ in Hz and ms were not found, nor its numerical area-law significance in the PRL main text. The 4.4σ figure comes from a search excerpt of the PDF.
- The GW230529 component masses were not verified.
- The Λ₁.₄ = 190 attribution (LVC PRL 121, 161101, arXiv:1805.11581) is cited by later work. I did not open the primary.
- Scheel et al. 2009 is cited through HBR16 for a_f = 0.68646 ± 0.00004. The corresponding NR mass (about 0.9516) is my recollection and is not used as a test.
- UIB2016 (Jiménez-Forteza et al. 2017) final-mass and spin fits, which LALSuite uses, were not fetched: git.ligo.org blocked the automated fetch.

## Relativistic geodesy from the gravity side (geopotential differences and height systems)

### Takeaway
`geopotential` computes the gravitational potential V but has no centrifugal term, no normal gravity and no W0, so it cannot yet produce the geopotential numbers or heights that clocks measure. Four small, verified pieces close that gap:
- the IHRS W0;
- GRS80 Somigliana normal gravity;
- the geopotential number C = W0 − W;
- the conversion Δf/f = ΔW/c².

The permanent-tide corrections (T4) go with them. Clock physics itself is left to the light researcher.

### Cited Findings

**R1. W0, the gravity potential and geopotential numbers.** New `geopotential::height::{gravity_potential, geopotential_number, normal_height}`. Links geopotential, metrology and gravimetry. Status: established (IAG 2015 resolution; realisation strategy 2021).
- W = V + ½ω²(x² + y²); C_P = W0 − W_P; normal height H* = C/γ̄ and orthometric height H = C/ḡ.
- W0(IHRS) = 62 636 853.4 m²s⁻² (IAG Resolution 1, Prague 2015) ([IAG resolutions](https://geodesy.science/wp-content/uploads/Prague-Czech-Republic-2015.pdf); [Springer 2024 chapter](https://link.springer.com/10.1007/1345_2024_244); [Sánchez et al. 2021, IHRS realisation strategy](https://e-docs.geo-leo.de/entities/publication/bcb2325a-06e4-4260-96e8-256a6d48bf94)).
- A chronometric-geodesy review notes that the IHRS W0 is inconsistent with the value underlying Terrestrial Time ([Delva, Denker & Lion, "Chronometric geodesy: methods and applications", arXiv:1804.09506](https://arxiv.org/abs/1804.09506)).
- Test *(computed)*: L_G·c² with L_G = 6.969290134×10⁻¹⁰ gives 62 636 856.0 m²s⁻². That is 2.60 m²s⁻² above the IHRS value, or 0.265 m in height. The L_G value is IAU 2000 B1.9, from memory, and the review confirms the inconsistency.
- The `gravimetry::sim::EARTH_ROTATION` constant (7.292115×10⁻⁵ rad/s) already exists and should move to `metrology`.

**R2. Normal gravity (GRS80 Somigliana).** New `geopotential::normal::gamma(lat)`. Status: established (Moritz 2000; constants from memory but self-consistent to 10⁻⁹).
- γ(φ) = γ_e(1 + k sin²φ)/√(1 − e² sin²φ), with γ_e = 9.7803267715 m/s², k = 0.001931851353 and e² = 0.00669438002290.
- Test *(computed)*: γ(90°) = 9.8321863685 m/s², the GRS80 polar value, to better than 10⁻⁹. U0(GRS80) = 62 636 860.850 m²s⁻² sits 0.76 m above the IHRS W0 in height *(computed; U0 from memory)*.

**R3. Clock-measured potential differences (gravity side only).** Doc test. Status: recent and unreplicated for this link; chronometric levelling in general is replicated.
- Δf/f = ΔW/c². The 457 km strontium fibre link measured ΔW = 3918.1(2.6) m²s⁻², against a geodetic 3915.88(0.30) ([PR Applied 21, L061001 (2024)](https://journals.aps.org/prapplied/abstract/10.1103/PhysRevApplied.21.L061001), already in the gravimetry-frontiers report).
- Tests *(computed)*: 3918.1 m²s⁻² is 4.359×10⁻¹⁴ in fractional frequency, about 399.5 m of height. A 1 mm height change is about 1.1×10⁻¹⁹, as already noted in the repo.

### Inferences
- R1, R2 and T4 together let `geopotential` output exactly what a clock comparison measures: a zero-tide or mean-tide geopotential number relative to a named W0. They also link to `gravimetry`, because normal and orthometric heights need γ̄ and ḡ, and to `metrology` for c, ω and W0.
- The 2.6 m²s⁻² gap between the IHRS W0 and L_G·c² must be explicit in the API, as two named constants. A clock-based height compared against an IHRS height would otherwise carry a silent 26.5 cm offset.

### Gaps
- Moritz 2000 (GRS80) and IAU 2000 B1.9 were not fetched. Their constants pass internal consistency checks (γ_p to 10⁻⁹; L_G·c² gives the IERS W0 to 0.1 m²s⁻²), but cite the primaries before shipping.
- The Helmert orthometric correction (ḡ ≈ g + 0.0424×10⁻⁵ H) and the 2026 IAG paper on commission and omission errors at clock sites were not checked this session.
