# Laboratory gravity generation, manipulation and measurement that physically works (or that GR allows): formulas, published test values, and status labels

Scope: fields, forces and strains people actually create or measure with controlled lab masses; what GR allows in principle (gravitomagnetism, frame dragging, gravitational waves); and the tabletop quantum-gravity frontier. Every item is labelled **[mainstream]**, **[contested]** or **[refuted]**.

Conventions:
- SI units throughout. G = 6.67430×10⁻¹¹ m³ kg⁻¹ s⁻², c = 299 792 458 m/s, ħ = 1.054571817×10⁻³⁴ J s.
- "computed" means I recomputed the number with `python3 -I`. The script is saved next to these notes as `lab_gravity_checks.py`; its sections A–M match the labels used below.
- Published values are quoted exactly as given.
- Most numbers come from arXiv full texts (PDF converted locally with pdftotext). Journal references are given where the text or the search result stated them. Where I did not verify a journal reference, I say so.

Sister notes, so this file does not repeat them in full:
- `../Gravitational waves and gravity beams/gravity_beams_and_lab_generation.md`: lab GW generation, Gertsenshtein effect, gravity of light beams, GP-B and LARES.
- `../Gravity beams fringe physics and DIRDs/superconductor_gravity_claims_math.md`: Podkletnov, Tajmar, Ning Li, DeWitt, London moment.

---

## Q1. Newtonian calibrators (NCal) at LIGO and Virgo: formula for the injected force, displacement and strain, and the published amplitudes

### Takeaway
A rotating rotor with two-fold (quadrupole) and three-fold (hexapole) mass patterns is the most precisely controlled man-made time-varying gravitational field in routine use **[mainstream]**.
- It moves a 40 kg LIGO test mass by about 5×10⁻¹⁷ m at 15.6 Hz (a strain of about 1.3×10⁻²⁰, computed). The force is about 19 pN, predicted to under 1%.
- Virgo's O4 NCal is now the absolute calibration reference, with a 1-σ uncertainty of 0.17% on the injected amplitude.
- The closed-form point-mass formulas below reproduce the published model values to 0.2% (computed), so they are ready-made unit tests for a Rust calculator.

### Cited Findings

**LIGO Hanford NCal: Ross et al. 2021, "Initial Results from the LIGO Newtonian Calibrator", arXiv:2107.00141 [mainstream]** — [arXiv:2107.00141](https://arxiv.org/abs/2107.00141)
- **Hardware:** an aluminium disk 25.4 cm in diameter and 5.08 cm high, with five tungsten cylinders (r_c = 1.968 cm, l_c = 5.08 cm).
  - Two cylinders sit at r_q = 6.033 cm (quadrupole, force at 2f). Three sit at r_h = 10.478 cm (hexapole, force at 3f).
  - The rotor is mounted in air outside the end-X test-mass vacuum chamber.
  - The force from a two-fold pattern falls as ~1/d⁴; the force from a three-fold pattern falls as ~1/d⁵.
- **Table I parameters (mean ± uncertainty):**

  | Parameter | Value |
  |---|---|
  | Tungsten cylinder mass m | 1.0558 kg ± 0.3 g |
  | Test mass M | 39.680 kg ± 10 g |
  | Rotor radial position ρ | 1180.2 mm ± 2.4 mm |
  | Rotor angle Φ | 52.24° ± 0.08° |
  | Rotor height z | 10.0 mm ± 3.0 mm |
  | Test-mass radius | 170.11 mm |
  | Test-mass length | 199.8 mm |

- **Point-mass formulas (eqs. 2–3)**, with ζ = 1 + (z/d)², θ = 2πft, and z/d ≪ 1 (ρ ≈ d):
  - `F_x(2f) = (9 G M m r_q² / (2 d⁴)) · ζ^(−5/2) · [ (5/(6ζ) − 2/3) cos(2θ − Φ) + (5/(6ζ)) cos(2θ − 3Φ) ]`
  - `F_x(3f) = (15 G M m r_h³ / (2 d⁵)) · ζ^(−7/2) · [ (7/(8ζ) − 3/4) cos(3θ − 2Φ) + (7/(8ζ)) cos(3θ − 4Φ) ]`
  - Units: N, with m in kg and r, d in m.
  - The fraction layout was garbled by text extraction. I reconstructed it, and the reconstruction reproduces the paper's Table II to 0.2% (see Inferences).
- **Predicted amplitudes (Table II):**

  | Model | F_x(2f) (pN) | F_x(3f) (pN) |
  |---|---|---|
  | Finite element | 19.18 ± 0.14 | 9.07 ± 0.09 |
  | Multipole (l ≤ 11) | 19.16 ± 0.14 | 9.06 ± 0.09 |
  | Point mass | 19.04 ± 0.15 | 8.97 ± 0.09 |

  Model uncertainty is about 0.75% (2f) and about 0.95% (3f).
- **Measured amplitudes (Table IV), corrected for torque coupling:**
  - Rotor frequencies 7.80, 8.55 and 9.56 Hz (lines at 15.60/23.40, 17.11/25.66 and 19.11/28.67 Hz). The 4.16, 5.17 and 6.36 Hz runs did not give useful 2f values.
  - F_x(2f) values include 19.28 ± 0.30 pN (±1.56%) and 19.10 ± 0.47 pN (±2.44%). The combined 2f value is 19.28 ± 1.12%.
  - F_x(3f) values range from 8.79 ± 0.32 to 9.15 ± 0.24 pN.
  - There is one excursion at 10.34 Hz: 16.59 ± 2.33 pN, where the detector noise is worse.
- **Corrections and conversion:**
  - The apparent force from torque coupling is 0.339 ± 0.026 pN (2f) and 0.199 ± 0.015 pN (3f). Magnetic coupling is < 0.04 pN.
  - Force is recovered from strain as `F_x = L h / S` (eq. 6). Above 10 Hz the test mass is treated as free, `S = −1/(M (2πf)²)`, with an error of no more than 0.4%.
- **G from LIGO:** the same data give **G = (6.56 ± 0.41)×10⁻¹¹ m³ kg⁻¹ s⁻²**, "the first gravitational constant measurement utilizing a gravitational wave detector".

**Virgo NCal, first tests: Estevez et al. 2018, arXiv:1806.06572 [mainstream]** — [arXiv:1806.06572](https://arxiv.org/abs/1806.06572)
- **Two point masses m at ±r, rotor on the beam axis:**
  - `F = (9 G M m r² / (2 d⁴)) cos(2θ)` (eq. 4).
  - Free-mass displacement amplitude at f_h = 2f_rotor: `a(f_h) = R·P / f_h²`, with `R = G m r² / (8π²)` and `P = 9/d⁴` (eqs. 5–7).
  - "This motion is independent of the mirror mass."
- **Off-axis rotor** (angle Φ, height z), eq. 8:
  - `P = [ (15/Z − 6)² cos²Φ + 6² sin²Φ ]^(1/2) · Z^(−5/2) · d^(−4)`, with `Z = 1 + (z/d)²`.
  - At Φ = π/2 the signal is still 2/3 of the maximum. At z = d/2, P no longer depends on Φ.
- **Extended rotor** (two sectors of opening angle α, thickness b, density ρ, radii r_min to r_max), eq. 10:
  - `R = (G b ρ / (8π²)) · sin α · (r_max⁴ − r_min⁴)/4`. This is maximal at α = π/2.
- **Test configuration and result:**
  - Rotor: two 45° Al-7075 sectors (2.80 kg/litre), r_min = 55.0 mm, r_max = 190.0 mm, thickness 52.0 mm, "each sector weighs about 2.06 kg".
  - Position: d = 1.32 m, Φ = 0.241 rad, z = 0.643 m below the beam plane.
  - Expected displacement: **3.3×10⁻¹⁴ (1 Hz/f)² m**.
  - The line SNR was about 50 in a 5 s FFT. Two models agreed to better than 2.2%.

**Virgo O3 and O4 NCal: Aubin et al. 2024, "The Virgo Newtonian calibration system for the O4 observing run", arXiv:2406.10028 [mainstream]** — [arXiv:2406.10028](https://arxiv.org/abs/2406.10028)
- **Earlier results quoted:** the O3 test NCal accuracy was 1.4%, with a 3% difference from the photon calibrator (PCal). The LIGO O3 NCal reached 0.8% total uncertainty.
- **O4 layout:**
  - Two triplets of rotors around the north-end mirror: "Far" aluminium rotors at 2.1 m and "Near" PVC rotors at 1.7 m.
  - Nominal angles ϕ = 34.5° and twist ψ = 12°.
  - Each rotor has two 90° sectors, diameter 208 mm, thickness 104.4 mm. Rotor mass is 5 kg (Al) or 2.5 kg (PVC).
- **Operation:** lines at about 36 Hz (rotor at about 18 Hz). PVC rotors are limited to 50 Hz.
- **Uncertainties:**
  - The signal goes as distance⁻⁴, so a 0.5 mm survey error gives 4δD/D = 0.12% (Near) and 0.1% (Far).
  - The uncertainty on G contributes 0.002%.
  - **Total 1-σ uncertainty: 0.17%** on the injected mirror displacement for the Near pair at about 36 Hz. This is better than PCal.

**Other NCal context [mainstream]:**
- KAGRA's 2f+3f rotor (Inoue et al.) and Matone et al. are the design precursors. Using two frequencies lets d be inferred from the ratio of force amplitudes. — [Ross et al. 2021](https://arxiv.org/abs/2107.00141)

### Inferences
- **Formula check (computed):**
  - My reconstruction of Ross eq. 2 gives 19.07 pN and eq. 3 gives 8.978 pN. The paper's point-mass values are 19.04 and 8.97 pN.
  - A brute-force point-mass sum over the actual 2- and 3-cylinder layout gives 19.10 and 9.02 pN.
  - So the closed forms match a direct sum to about 0.2–0.5%.
  - The Estevez form `[(15/Z−6)² cos²Φ + 36 sin²Φ]^(1/2)` is algebraically identical to Ross's two-cosine form. Use `|A e^(−iΦ) + B e^(−3iΦ)|` with A + B = 5/(3ζ) − 2/3 and A − B = −2/3, then multiply by 9.
- **Displacement and strain injected at LIGO (computed):** using F = 19.16 pN, F₃ = 9.06 pN, M = 39.68 kg and an assumed arm length L = 3995 m (standard LIGO value, not taken from the paper):

  | Line | Displacement x (m) | Strain h |
  |---|---|---|
  | 15.60 Hz | 5.03×10⁻¹⁷ | 1.26×10⁻²⁰ |
  | 23.40 Hz | 1.06×10⁻¹⁷ | 2.64×10⁻²¹ |
  | 19.12 Hz | 3.35×10⁻¹⁷ | 8.4×10⁻²¹ |

- **Estevez 2018 check (computed):** R·P = 3.331×10⁻¹⁴ m Hz², matching the published 3.3×10⁻¹⁴. However, the sector volume from the stated dimensions gives 1.89 kg per sector, not "about 2.06 kg". The geometry may include extra material not described; flag this in tests.
- **Virgo O4 rough estimate (computed):** assumptions are r_min ≈ 0 (consistent with the 4.98 kg computed Al rotor mass), z = 0, no twist, f_h = 36 Hz, and L = 3000 m.
  - Far Al rotor: x ≈ 2.3×10⁻¹⁸ m (h ≈ 7.8×10⁻²²).
  - Near PVC rotor (density about 1400 kg/m³ assumed): x ≈ 2.7×10⁻¹⁸ m (h ≈ 9×10⁻²²).
  - These are order-of-magnitude only; the paper does not print the absolute amplitude.
- **Calculator recipe:**
  1. `F(2f)` and `F(3f)` from the formulas above.
  2. `x = F/(M(2πf)²)`, valid well above the 0.6 Hz suspension resonance (Virgo; Estevez 2018).
  3. `h = x/L`.
  - Test vectors: Ross Table II/IV, Estevez 3.3×10⁻¹⁴ (1 Hz/f)² m, and G = (6.56 ± 0.41)×10⁻¹¹.
- **Physical meaning:** NCals are near-field Newtonian gravity, not gravitational waves (the GW wavelength at 15 Hz is 2×10⁷ m). They are nonetheless the "gravity generators that work", producing about 10⁻²⁰ effective strain from 5 kg of rotating metal.

### Gaps
- I did not fetch Estevez et al. 2021 (Virgo O3 NCal, arXiv:2011.03728) or the Datrier et al. derivation note. The O3 numbers come only via Aubin et al. 2024.
- Journal references for Ross 2021 (believed to be PRD 104, 082006) and Estevez 2018 (CQG) were not verified on publisher pages.
- No published absolute Virgo O4 displacement amplitude was found in the text; the O4 estimate above is my own.

---

## Q2. Measuring the gravity of small and controlled masses: Westphal 2021 (90 mg), atom-interferometer G (Rosi 2014), gravitational Aharonov–Bohm (Overstreet 2022), and 2022–2024 successors

### Takeaway
The smallest single source mass whose gravity has been resolved is a 92.1 mg gold sphere (Westphal et al., 2021) **[mainstream]**.
- The signal is about 3×10⁻¹⁰ m/s² of modulated acceleration, giving G to 1% statistics but 9% low.
- Atom interferometry gives G = 6.67191(99)×10⁻¹¹ (150 ppm) from a 516 kg tungsten source (Rosi et al., 2014).
- Overstreet et al. (2022) showed that a 1.25 kg tungsten mass imprints a gravitational Aharonov–Bohm phase on ⁸⁷Rb wave packets 25 cm apart.
- A 0.43 mg levitated magnet detected the gravity of 2.45 kg brass masses (Fuchs et al., 2024).
- The next milestone is a Planck-mass (22 µg) source.

### Cited Findings

**Westphal, Hepach, Pfaff, Aspelmeyer, "Measurement of gravitational coupling between millimetre-sized masses", arXiv:2009.09546v2 (Nature 2021) [mainstream]** — [arXiv:2009.09546](https://arxiv.org/abs/2009.09546)
- **Masses:**
  - Source: gold sphere, r = (1.07 ± 0.04) mm, m_s = (92.1 ± 0.1) mg.
  - Test mass: m_t = (90.7 ± 0.1) mg. Counterbalance: m_a = 91.5 mg, on a 37.9 mm glass capillary beam with a 4 µm silica fibre.
  - The pendulum is in vacuum at 6×10⁻⁷ mbar.
- **Pendulum:** f₀ = 3.59 mHz, Q = 4.9, effective mass 183.8 mg. Readout noise is about 2 nm/√Hz.
- **Drive:** the source is moved by about 3 mm at f_mod = 12.7 mHz. Both the linear signal and the 2f_mod = 25.4 mHz signal (from nonlinearity of the potential) are observed.
  - Centre separation varied from 2.5 mm to 5.8 mm, with a minimum surface gap of about 0.4 mm.
  - The source is discharged to below 8×10⁴ e. Electrostatic forces are below 3% of gravity.
- **Scale:** "at a center separation of 2.5 mm the static gravitational force between our masses is expected to be 9 × 10⁻¹⁴ N". A human 2.5 m away exerts the same force.
- **Results:**
  - Single fit: G_fit = (5.89 ± 0.20)×10⁻¹¹ N m²/kg².
  - Combined: **G_comb = (6.04 ± 0.06)×10⁻¹¹** (1% statistical). This is about 9% below CODATA, "fully covered by the known systematic" uncertainties below 10%.
  - The modulated acceleration resolved is ≈ 3×10⁻¹⁰ m/s², with systematic uncertainty ≈ 3×10⁻¹¹ m/s² and statistical uncertainty ≈ 3×10⁻¹² m/s² over 350 h.
- **Outlook:** a source at the Planck mass, m_P ≈ 22 µg = 2.4×10⁻⁴ m_s, "should become possible for Q ≳ 20,000". Also, SNR ≈ 100 would allow about 1 mg sources at SNR ≈ 1.
- **Previous small sources** (cited in the same paper): a hole-patterned 5 cm Pt attractor with a few hundred mg of mass modulation, and a single 700 mg tungsten sphere.

**Rosi, Sorrentino, Cacciapuoti, Prevedelli, Tino, "Precision measurement of the Newtonian gravitational constant using cold atoms", arXiv:1412.7954 (Nature 510, 518, 2014) [mainstream]** — [arXiv:1412.7954](https://arxiv.org/abs/1412.7954)
- **Result:** **G = 6.67191(99)×10⁻¹¹ m³ kg⁻¹ s⁻²** (150 ppm), written as 6.67191(77)(65), statistical then systematic. It is 1.5 combined σ from the then-CODATA value 6.67384(80)×10⁻¹¹.
- **Source mass:** 24 tungsten-alloy cylinders, total about 516 kg. Each has diameter 99.90 mm and height 150.11 mm, in a hexagonal arrangement on two platforms.
  - Configuration C: cylinder centres 40 mm and 261 mm from the lower-cloud apogee. Configuration F: −74 mm and 377 mm.
- **Interferometer:** ⁸⁷Rb Raman gradiometer with T = 160 ms. Apogees are about 60 cm and 90 cm above the MOT, and phase goes as φ = k·gT².
- **Signal:** sensitivity 3×10⁻⁹ g per √(1 s). The measured differential phase shift is **0.547870(63) rad**, of which 97.0% is from the cylinders, 2.8% from the supports and 0.2% from other moving masses. Data: 100 h in July 2013.

**Overstreet, Asenbaum, Curti, Kim, Kasevich, "Observation of a gravitational Aharonov–Bohm effect", Science 375, 226 (2022) [mainstream]** — [Science DOI 10.1126/science.abl7152](https://www.science.org/doi/10.1126/science.abl7152); [Science News summary](https://www.sciencenews.org/article/quantum-particles-gravity-spacetime-aharonov-bohm-effect)
- **Setup:** ⁸⁷Rb wave packets about 25 cm apart. One packet passes near a kg-scale tungsten source mass. The phase shift deviates from the deflection-induced contribution, as quantum mechanics predicts. Search snippets give closest approach as about 7.5 cm.
- **Companion paper (same group plus Giacomini), "Inference of gravitational field superposition from quantum measurements", arXiv:2209.02214v3 (2023):** — [arXiv:2209.02214](https://arxiv.org/abs/2209.02214)
  - Source mass: tungsten, **1.25 kg semicircular ring**. The upper arm is launched to 4 cm above it; closest approach 7 cm.
  - Interferometer: 52ħk, 25 cm wave-packet separation. The gradiometer reference is 24 cm below.
  - Measured phase: **−0.23, −0.24, −0.27 rad** for arm probability 0.25, 0.5, 0.75 (slope −0.07 ± 0.05, p = 0.39; reduced χ² = 1.1 against QM).
  - A semiclassical "field-energy" model predicts −0.198, −0.374, −0.394 rad and is excluded.
  - The authors' claim that this is "indirect evidence that gravity has quantum features" is an interpretation **[contested]**.
- Wikipedia (secondary) also states a 1.25 kg tungsten source. — [Wikipedia](https://en.wikipedia.org/wiki/Gravitational_Aharonov-Bohm_effect)

**Fuchs et al., "Measuring gravity with milligram levitated masses", Sci. Adv. 10, eadk2949 (2024) [mainstream]** — [PMC full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC10889343)
- **Test mass:** 0.43 mg (three 0.25 mm Nd₂Fe₁₄B cubes plus a glass bead), levitated in a type-I tantalum superconducting trap (T_c = 4.48 K).
- **Mode:** z-mode at 26.7 Hz with Q = 9.13×10⁶ (lower bound; τ = 1.09×10⁵ s), k = 12×10⁻³ N/m, γ/2π = 2.9 µHz.
- **Source:** three 2.45 kg brass masses on a wheel, at a vertical separation of 48 ± 4 cm. The drive is at the mode resonance.
- **Signal:**
  - The introduction states "a coupling of 1030 aN with a force noise of 0.5 fN/√Hz". The discussion states "a 30-aN gravitational signal at 27-Hz". The paper does not reconcile these two values.
  - Agreement with simulation is given as a factor of 0.35 ± 0.02.
  - Operation was below 100 mK, with a mode temperature of 3 K.

**Dynamic near-field gravity at Hz frequencies [mainstream]** (relayed from the sister note, not re-fetched):
- Brack et al. 2022, Nature Physics 18, 952: a beam vibrating near 42 Hz drives a detector beam through gravity alone, with about 1 pm displacement and agreement with Newtonian modelling to about 2%.
- A rotating-bar follow-up appeared in Communications Physics (2023). — [APS Physics](https://physics.aps.org/articles/v15/121)

**Reference value:** CODATA G = 6.67430(15)×10⁻¹¹ m³ kg⁻¹ s⁻² (quoted in Westphal et al.). — [arXiv:2009.09546](https://arxiv.org/abs/2009.09546)

### Inferences
- **Westphal (computed):**
  - Static force at 2.5 mm: G m_s m_t/d² = 8.92×10⁻¹⁴ N (paper: 9×10⁻¹⁴ N).
  - An 80 kg person at 2.5 m gives 7.7×10⁻¹⁴ N on m_t, consistent with the paper's comparison.
  - Linear modulation with d₀ = 4.2 mm and d_m = 1.6 mm: `2 G m_s d_m / d₀³` = 2.66×10⁻¹⁰ m/s² (paper: ≈ 3×10⁻¹⁰).
  - Planck mass √(ħc/G) = 21.76 µg.
- **Rosi (computed):** with k_eff = 2·(2π/780.24 nm), k_eff T² = 4.12×10⁵ rad per m/s². So 0.547870 rad corresponds to a source-induced change in differential acceleration of about 1.33×10⁻⁶ m/s², roughly 1.4×10⁻⁷ g, between configurations C and F. This is an inference: the actual extraction uses the full trajectory model.
- **Overstreet order of magnitude (computed):** m_Rb·GM/(ħr) for M = 1.25 kg is 1.63 rad/s at r = 7 cm and 0.46 rad/s at 25 cm. An arm-to-arm potential difference of about 1 rad/s, applied for a fraction of a second near the mass, gives phases of about 0.1–1 rad. This is consistent with the measured −0.23 to −0.27 rad.
- **Calculator set:**
  - `F = G m₁ m₂ / r²` (with finite-size corrections for spheres: exact for uniform spheres outside contact).
  - Modulation harmonics: `G m_s d_m / d₀³ · 2` at f_mod; the 2f term is `3 G m_s d_m² / (2 d₀⁴)` (my Taylor expansion, not from the paper).
  - Atom-interferometer phase: `φ = k_eff a T²`.

### Gaps
- I could not read the Science 2022 paper itself (paywall). Its quoted phase versus pulse timing, and the 12ħk or other configurations, are not confirmed here; the 52ħk numbers come from the arXiv companion.
- The Fuchs et al. 1030 aN versus 30 aN inconsistency is in the paper itself (as extracted). I did not resolve which number is the measured gravitational force amplitude.
- No source-mass result below 90 mg was found up to Oct 2026.

---

## Q3. Short-range tests of the inverse-square law: the Yukawa parameterisation, Eöt-Wash bounds, Casimir-limited searches, and α–λ limits at several λ

### Takeaway
Gravity obeys 1/r² to the extent that any gravitational-strength (|α| = 1) Yukawa term must have λ < 38.6 µm (95% CL; Lee et al. 2020, Eöt-Wash) **[mainstream]**.
- Below about 10 µm the limits weaken quickly: |α| < 1.4×10⁴ at 10 µm.
- Casimir and patch-potential backgrounds take over below a few µm, and neutron scattering below about 10 nm.
- The dark-energy length λ_d = (ħc/ρ_d)^(1/4) ≈ 85 µm (computed 84.9 µm) is now probed down to below half its value.

### Cited Findings

**Parameterisation [mainstream]:**
- Potential: `V(r) = V_N(r)[1 + α exp(−r/λ)]`. — [Lee et al. 2020, arXiv:2002.11761](https://arxiv.org/abs/2002.11761)
- Point-mass force (derivative of the above; computed): `F(r) = G m₁ m₂/r² · [1 + α (1 + r/λ) e^(−r/λ)]`.
- Sphere–plate force difference used by Decca's group: `ΔF_h(z) = −4π² G α λ³ e^(−z/λ) R K_t K_s`, where K_t and K_s are density and thickness factors of test and source. — [Chen et al. 2016, arXiv:1410.7267](https://arxiv.org/abs/1410.7267)

**Eöt-Wash, Lee, Adelberger, Cook, Fleischer, Heckel 2020, "New test of the gravitational 1/r² law at separations down to 52 µm", arXiv:2002.11761 (PRL 2020) [mainstream]** — [arXiv:2002.11761](https://arxiv.org/abs/2002.11761)
- **Apparatus:** a stationary torsion pendulum and a rotating attractor, with 18-fold and 120-fold hole patterns (torques at 18ω and 120ω).
  - Tungsten test bodies cut from 50 µm foil; detector and attractor thicknesses 54 and 99 µm.
  - A 10 µm gold-coated BeCu electrostatic shield foil sits between them.
  - Separations 52 µm to 3.0 mm.
- **Calibration test value:** measured torque (2.137 ± 0.009) fN·m versus expected (2.112 ± 0.005) fN·m, giving γ = 0.988 ± 0.005.
- **Results:**
  - Newtonian fit P = 0.654. The best Yukawa fit is Δχ² = 3.3 at λ = 7.1 µm (not significant).
  - **Any |α| = 1 Yukawa must have λ < 38.6 µm (95% CL).**
  - Dilaton or heavy-graviton mass > 5.1 meV; radion unification mass > 7.1 TeV; largest extra dimension (toroidal radius) < 30 µm.
  - It also quotes λ_d = (ħc/ρ_d)^(1/4) ≈ 85 µm with ρ_d ≈ 3.8 keV/cm³, and "percent-level measurements of G_N at separations down to about 50 µm".

**Earlier Eöt-Wash: Kapner et al. 2007, PRL 98, 021101 [mainstream, superseded]** — [arXiv:hep-ph/0611184](https://arxiv.org/abs/hep-ph/0611184)
- Separations 9.53 mm to 55 µm. |α| ≤ 1 down to λ = 56 µm (95% CL); extra dimension R ≤ 44 µm.

**About 10 µm: Geraci, Smullin, Weld, Chiaverini, Kapitulnik 2008, PRD 78, 022002 [mainstream]** — [arXiv:0802.2350](https://arxiv.org/abs/0802.2350)
- "|α| > 14,000 are excluded at λ = 10 microns" (95% CL); best bound in 5–15 µm at that time.

**0.4–4 µm: Sushkov, Kim, Dalvit, Lamoreaux 2011, PRL 107, 171101 [mainstream]** — [arXiv:1108.2547](https://arxiv.org/abs/1108.2547)
- Casimir torsion pendulum, separations 0.7–7 µm. New best bound for λ from 0.4 to 4 µm; (4+n)-dimensional Planck scale M* > 70 TeV. The abstract gives no α value.

**40–8000 nm, Casimir-less: Chen, Tham, Krause, López, Fischbach, Decca 2016, PRL 116, 221102 [mainstream]** — [arXiv:1410.7267](https://arxiv.org/abs/1410.7267)
- **Method:** a rotating Au/Si source mass under a uniform Au overlayer (λ_p = 135 nm structure) suppresses the Casimir difference by more than 10⁶.
- **Measured Casimir force:** F_C(200 nm) = 34 pN, scaling as z^(−2.78).
- **Result:** up to 10³ improvement in α, largest at λ = 300 nm. The text says |α| < 1 is reached only for λ > 50 µm.

**Neutron quantum states and neutron scattering [mainstream]:**
- Analytic limits from neutron quantum states in Earth's gravity: for λ > 10 µm, α < 2×10¹¹ (90% CL). — [A. Westphal et al. 2007, arXiv:hep-ph/0703108](https://arxiv.org/abs/hep-ph/0703108)
- Pendellösung interferometry in silicon improved Yukawa strength limits by an order of magnitude over 20 pm–10 nm. — [Heacock et al. 2021, Science, arXiv:2103.05428](https://arxiv.org/abs/2103.05428)
- Reviews say neutron scattering gives the strongest limits for λ < 9 nm. Lateral and normal Casimir measurements cover 9–31 nm, and Casimir-less measurements cover up to about 8 µm. — [Klimchitskaya & Mostepanenko review search summary, arXiv:1704.02920](https://arxiv.org/pdf/1704.02920); [arXiv:2103.06576](https://arxiv.org/pdf/2103.06576)

**Astrophysical complement [mainstream]:**
- Neutron-star cooling bounds scalar–nucleon couplings g_N ≲ 5×10⁻¹⁴ for ranges 10⁻⁶–10⁻¹² m (eV ≲ m_φ ≲ MeV). — [Fiorillo et al., PRL 135, 211003 (2025), arXiv:2506.19906](https://arxiv.org/abs/2506.19906)

### Inferences
- **α–λ test table:**

  | λ | Limit on \|α\| | Source |
  |---|---|---|
  | 38.6 µm | 1 | Lee 2020 |
  | 56 µm | 1 | Kapner 2007, superseded |
  | 10 µm | 1.4×10⁴ | Geraci 2008 |
  | > 10 µm (neutron quantum states) | 2×10¹¹ | Westphal 2007, a much weaker method |

  Below 4 µm, only relative improvements are quoted in the abstracts I read (see Gaps).
- **Lee 2020 checks (computed):**
  - λ_d = (ħc/ρ_d)^(1/4) = 84.9 µm (paper ≈ 85 µm; the printed "4√" is a fourth root, not a factor of 4).
  - Boson mass: ħc/λ for λ = 38.6 µm is 5.11 meV (paper 5.1 meV).
  - Force factor at r = λ with α = 1: 1 + 2/e = 1.736.
- **Coupling conversion (computed):** a scalar Yukawa between nucleons with coupling g maps to `α = g² ħc / (4π G m_N²)`. So the neutron-star bound g_N ≈ 5×10⁻¹⁴ corresponds to α ≈ 3.4×10¹⁰ at λ ~ 1 µm–1 pm. The convention (factor 4π, and nucleon versus atomic mass) should be checked before using it as a test value.
- **Why the curve is Casimir-limited:** at a 200 nm gap, F_C ≈ 34 pN. Gravity between the same µm-scale bodies is about 10⁻¹⁹–10⁻¹⁷ N, so any Yukawa search must subtract the Casimir force differentially (Decca) or shield it (Eöt-Wash foil). Electrostatic patches and vibration set the 52 µm floor at Eöt-Wash.

### Gaps
- I did not obtain tabulated |α| values at λ = 1 mm, 1 µm, 100 nm, 10 nm or 1 nm from a primary source. The best route is the Lee 2020 Supplemental Material (+α and −α curves) or a recent PDG review.
- The HUST (Tan et al. 2016/2020, PRL) bounds, believed to be λ ≈ 48–59 µm at |α| = 1, were not retrieved.
- I found no 2021–2026 Eöt-Wash result below 52 µm; the authors stated they were adding active vertical vibration isolation.

---

## Q4. Robert Forward's GR-based "gravity control" papers, gravitomagnetic fields from mass currents, his dipole and torus designs, and frame dragging from rotating lab masses

### Takeaway
Forward's 1961 and 1963 papers are legitimate linearised GR **[mainstream physics]**. Mass currents produce gravitomagnetic fields; accelerated mass flow in a torus produces a gravitational "dipole" field.
- The coupling coefficient is G/c² = 7.4×10⁻²⁸ m/kg (computed). This matches the "order 10⁻²⁸" quoted for his dipole term.
- Any lab realisation is therefore about 20–25 orders of magnitude below Earth's 1 g. Forward himself framed the designs as guidelines, not devices.
- Frame dragging from the strongest feasible lab rotor (a 10 t tungsten sphere at its burst speed) is about 5×10⁻²¹ rad/s (computed). That is about 10⁶ times weaker than the Earth's frame dragging GP-B measured (6.0×10⁻¹⁵ rad/s).
- No lab frame-dragging experiment from Braginsky, Caves and Thorne (1977) has been performed.

### Cited Findings

**Bibliography [mainstream]:**
- Forward, R. L., "General Relativity for the Experimentalist", Proc. IRE 49, 892 (1961). A 2023 review says Forward used the electromagnetic analogy to refer to a formalism useful to deal with GR "in experimentally realisable conditions". — [Ruggiero 2023 review, arXiv:2304.02167](https://arxiv.org/abs/2304.02167)
- Forward, R. L., "Guidelines to Antigravity", Am. J. Phys. 31(3), 166–170 (1963), DOI 10.1119/1.1969340. — [AIP page](https://pubs.aip.org/aapt/ajp/article/31/3/166/1046748/Guidelines-to-Antigravity); [DOI as cited on Phys-L](https://www.phys-l.org/archives/2016/7_2016/msg00018.html)

**Forward's torus dipole term [mainstream as GR; the numbers are not practical]:**
- As reproduced by a 2013 paper, Forward's calculated gravitational dipole is `G_D = −(µ_G/4π) d/dt[ u (r/R)² ]` (eq. 28 there).
  - µ_G = 4πG/c² is the "gravitomagnetic equivalent of the permeability".
  - r and R are the internal and external radii of a confined toroidal flow, and u is the flow velocity.
  - The coefficient is described as "of the order of 10⁻²⁸", which makes measurement "the main difficulty". — [arXiv:1306.4893](https://arxiv.org/abs/1306.4893)
- A secondary source describes Forward's case of a fluid toroidal mass undergoing minor-axis rotational acceleration, which by gravitomagnetic arguments pulls matter through the throat. — [search summary of a gravitoelectromagnetism PDF, users.ntua.gr](http://users.ntua.gr/kyriakos/Gravitoelectromagnetism.pdf)

**GEM field definitions** (one common convention) **[mainstream]:**
- Potentials: `Φ = GM/r`, `A = (G/c)(J × x)/r³`. Fields: `B = ∇ × A`, `E = −∇Φ − (2/c)∂A/∂t`.
- Equation of motion (stationary case): `dv/dt = −E − 2(v/c) × B`.
- There are no gravitoelectromagnetic waves; GWs are the TT part. — [Ruggiero 2023, arXiv:2304.02167](https://arxiv.org/abs/2304.02167)

**Laboratory gravitomagnetism proposals:**
- Braginsky, Caves and Thorne (1977), "Laboratory experiments to test relativistic gravity", PRD 15, 2047, proposed seven lab experiments. One searches for frame dragging by a rotating body using a high-Q toroidal microwave cavity. — [Caltech record](https://authors.library.caltech.edu/records/re5v4-65528)
- "None of the proposed experiments … were performed up today" (2023). GINGER ring-laser gyroscopes aim at the Earth's Lense–Thirring and de Sitter terms. — [Ruggiero 2023](https://arxiv.org/abs/2304.02167)

**Measured frame dragging, Earth as the source [mainstream]** (relayed from the sister note):
- GP-B frame-dragging drift: −37.2 ± 7.2 mas/yr against GR's −39.2. — [Everitt et al. 2011, arXiv:1105.3456](https://arxiv.org/abs/1105.3456)
- LARES-2 with LAGEOS (Nature 2026): about 10⁻³ relative uncertainty. — [Nature 2026](https://www.nature.com/articles/s41586-026-10715-0)

**Lab "gravitomagnetic" claims far above GR [refuted or contested]:**
- Rotating-superconductor anomalies (Tajmar) were attributed to fields far larger than GR predicts. Such fields would conflict with neutron-star properties, and later Tajmar et al. work narrowed the constraints. — [Ruggiero 2023](https://arxiv.org/abs/2304.02167). Details are in the sister note `superconductor_gravity_claims_math.md`.

### Inferences
- **Dipole coefficient (computed):** G/c² = 7.426×10⁻²⁸ m/kg. This is exactly the "order 10⁻²⁸" coefficient, and it is the universal suppression factor of every Forward-type device.
- **Induction estimate (my linearised-GR scaling, not from Forward's text):** a mass m accelerating at a_f at distance r induces a gravitoelectric acceleration of order `a_ind ≈ 4 G m a_f / (c² r)`. This follows from h₀ᵢ ≈ −(4G/c³)∫ρvᵢ/r d³x′ and a = −c ∂ₜh₀ᵢ.
  - For m = 10⁴ kg, a_f = 10⁴ m/s² and r = 1 m: a_ind ≈ 3×10⁻¹⁹ m/s² (computed).
  - For 1 g = 9.8 m/s² with a_f = 10⁶ m/s², you need m/r ≈ 3.3×10²¹ kg/m (computed). Earth's m/r is about 9.4×10¹⁷ kg/m, so this is a compactness between Earth and a white dwarf.
  - So Forward's "antigravity" generators are physically allowed but need stellar-density, high-acceleration matter flows **[mainstream conclusion; inference]**.
- **Lense–Thirring from a lab sphere (computed):**
  - Use `Ω_LT(pole) = 2GJ/(c² r³)` with J = (2/5)MR²ω. Take a tungsten sphere with R = 0.5 m and M = 1.01×10⁴ kg, spun until the surface speed equals √(σ/ρ) for σ = 1 GPa (ω = 455 rad/s).
  - Result: Ω_LT ≈ 5.5×10⁻²¹ rad/s at the pole surface, or Ω_LT/ω ≈ 1.2×10⁻²³.
  - Inside a thin shell (M = 10⁴ kg, R = 1 m, ω = 300 rad/s): `Ω = 4GMω/(3c²R)` ≈ 3.0×10⁻²¹ rad/s.
  - GP-B's 39.2 mas/yr is 6.02×10⁻¹⁵ rad/s, so lab sources are about 10⁶ times weaker than Earth's field in low orbit. Earth-based gyroscopes also see Earth's 7.3×10⁻⁵ rad/s Coriolis background.
- **Rust calculators:**
  - `omega_lt(J, r_vec)`: (G/(c² r³))[3(J·r̂)r̂ − J], the standard gyroscope precession formula. It is not in the fetched text; it comes from GP-B practice.
  - `thirring_shell(M, R, ω)`.
  - `forward_coeff()` = G/c².
  - Test values: 5.47×10⁻²¹ rad/s and 2.97×10⁻²¹ rad/s for the configurations above (computed).

### Gaps
- I could not obtain the full text of Forward 1961 (IEEE) or Forward 1963 (AIP). Forward's own numbers are not verified, including any "1 g" design or the neutron-star-density torus often attributed to him, and the exact definition and units of u in the dipole formula.
- The Braginsky–Caves–Thorne cavity sensitivity numbers were not read.

---

## Q5. Laboratory generation of gravitational waves: the best realistic sources, power and strain, and 2015–2026 proposals

### Takeaway
No lab-generated GW has ever been detected **[mainstream]**.
- Rotating bars give powers of 10⁻⁴⁰ W (1 kg, 100 Hz) to 2×10⁻²⁹ W (the 490 t MTW beam) (computed).
- Electromagnetic schemes reach strains of about 10⁻³² (Gertsenshtein, 10 T, 100 m) or about 10⁻³⁶ (petawatt twisted light at µm). These are physically sound but 10 or more orders of magnitude from detection.
- Any lab detector below about 1 GHz sits in the near zone, where Newtonian coupling (the NCal and Brack et al. 2022 signals) dominates the radiative part.
- Claims of "laser-like" GW generation in superconductors are **[contested]**.

### Cited Findings
Relayed from the sister note `gravity_beams_and_lab_generation.md` §2, with sources as cited there:
- **Quadrupole formula** for a rigid body spinning at Ω: `P = (32/5)(G/c⁵)(I₁ − I₂)² Ω⁶` **[mainstream]**. — [Nordita/Stockholm lecture exercise](https://indico.fysik.su.se/event/8276/contributions/13166/attachments/5761/7341/stockholm24ex.pdf)
- **MTW steel beam:** about 20 m long, 1 m radius, 4.9×10⁵ kg, spun at the tensile limit **[mainstream]**. — [Chiao et al. arXiv:1712.08680](https://arxiv.org/pdf/1712.08680)
- **Weber's lab Hertz analysis:** up to 10⁻¹³ erg/s (10⁻²⁰ W), "much less than the receiver sensitivity" **[mainstream]**. — [Kolosnitsyn & Rudenko arXiv:1504.06548](https://arxiv.org/abs/1504.06548)
- **Gertsenshtein Hertz experiment [mainstream physics, unrealised]:** — [arXiv:1504.06548](https://arxiv.org/abs/1504.06548)
  - Amplitude: `h_max = 2G H₀ H̃₀ L/(c⁴ k)` (Gaussian units), enhanced by finesse.
  - Parameters H₀ = 10 T, L = 100 m, λ = 1 cm, 1 − R = 10⁻⁸ and a 10 kW pump give h ~ 10⁻³². Superconducting resonators with finesse ~ 10¹⁰ give h ~ 10⁻³⁰.
  - The receiver would see about 10⁻³ photons/s at h = 10⁻³², needing about 4×10⁵ s of integration at 4 K.
- **High-power lasers [mainstream physics]:**
  - Laser-driven shock waves: analytic metric perturbations (Kadlecová, Klimo, Weber, Korn, EPJ D 71, 89, 2017). — [arXiv:1602.08904](https://arxiv.org/abs/1602.08904)
  - Twisted light: h ∼ 10⁻³⁶ at µm distances for 1 PW. — [Atonga et al. arXiv:2309.04191](https://arxiv.org/html/2309.04191v2)
- **Superconducting "GW laser" [contested]:**
  - Chiao et al. 2017 claim coherent generation above a threshold, using a "relative gravitational permeativity" and Meissner-like GW reflection. — [arXiv:1712.08680](https://arxiv.org/abs/1712.08680)
  - Chiao, Fitelson and Speliotopoulos (2003) found no YBCO transducer signal, with an upper limit of 1.6×10⁻⁵. — [arXiv:gr-qc/0304026](https://arxiv.org/abs/gr-qc/0304026)
- **Near-field dynamic gravity (measured; not radiation) [mainstream]:**
  - Hirakawa group rotating-bar quadrupole–quadrupole coupling up to 4.2 m (Nature 283, 184, 1980). — [Nature 283, 184](https://www.nature.com/articles/283184a0.pdf)
  - Brack et al. 2022 (42 Hz, ~1 pm). — [APS Physics](https://physics.aps.org/articles/v15/121)
- **UHF-GW review:** no known astrophysical sources at MHz–GHz. — [Aggarwal et al. 2021, Living Rev. Relativ. 24, 4, arXiv:2011.12414](https://arxiv.org/abs/2011.12414)

### Inferences
- **Recomputed (computed):**
  - MTW beam (I₁ = M(L²/12 + r²/4) = 1.645×10⁷ kg m², I₂ = Mr²/2, Ω = 28 rad/s): **P = 2.23×10⁻²⁹ W**. The far-zone strain at one GW wavelength (πc/Ω = 3.4×10⁷ m) is h = 4G(I₁−I₂)Ω²/(c⁴r) ≈ **1.2×10⁻⁴¹**.
  - The "strain at 1 m" from the same formula is 4.2×10⁻³⁴. It is formal only, because 1 m is deep in the near zone.
  - A 1 kg, 0.2 m bar rotating at 100 Hz: P = 1.2×10⁻⁴⁰ W.
- **Comparison:** the LIGO NCal (Q1) produces an effective test-mass strain of about 1.3×10⁻²⁰ from 5 kg at 1.2 m. That is about 14 orders of magnitude above the formal radiative h at 1 m of the MTW beam, because it is near-field Newtonian.
- **Best radiative source on paper:** Gertsenshtein in a high-field, high-finesse cavity (h ~ 10⁻³⁰–10⁻³², about 10⁻¹² W). Its problems are EM isolation and the lack of a GHz GW detector with the required sensitivity.
- **Power per mass:** P ∝ M²Ω⁶ L⁴. Material strength caps ΩL (the surface speed) at about √(σ/ρ) ≈ 200–1000 m/s, so mechanical sources cannot close the gap. This is a scaling argument (inference).

### Gaps
- I did not re-fetch the 2015–2026 proposal papers for this note; numbers are relayed from the sister note.
- I found no lab GW-generation experiment reporting a positive or limit result after 2003, apart from Chiao's null.

---

## Q6. Tabletop quantum gravity: the BMV entanglement proposal, experimental progress to 2026, gravitational decoherence (Diósi–Penrose, Donadi 2021), and single-graviton detection (Tobar 2024; Dyson 2013)

### Takeaway
The BMV proposal needs m ≈ 10⁻¹⁴ kg masses in 250 µm superpositions, 200 µm apart, for about 2.5 s at 10⁻¹⁵ Pa and 0.15 K **[mainstream proposal]**.
- That gives entangling phases of order 0.1–0.8 rad.
- No gravitationally induced entanglement has been observed as of Oct 2026. Enabling steps are levitated-particle ground-state control (2025) and gravity from 0.4–90 mg sources.
- Donadi et al. (2021) excluded the parameter-free Diósi–Penrose model (R₀ > 0.54×10⁻¹⁰ m against the "natural" 0.05×10⁻¹⁰ m) **[refuted for that version]**.
- Tobar et al. (2024) showed that a ground-state-cooled acoustic resonator of about 15 kg could register the stimulated absorption of single gravitons from GW170817-like signals **[mainstream proposal]**. Whether a click proves quantisation is **[contested]**.
- Dyson's 2013 conjecture is that single gravitons are undetectable in LIGO-type interferometers.

### Cited Findings

**Bose, Mazumdar, Morley, Ulbricht, Toroš, Paternostro, Geraci, Barker, Kim, Milburn, "Spin entanglement witness for quantum gravity", arXiv:1707.06050 (PRL 119, 240401, 2017) [mainstream proposal]** — [arXiv:1707.06050](https://arxiv.org/abs/1707.06050)
- **Phases:** `φ_RL ~ G m₁ m₂ τ/(ħ(d − Δx))`, `φ_LR ~ G m₁ m₂ τ/(ħ(d + Δx))`, `φ ~ G m₁ m₂ τ/(ħ d)`, with Δφ_RL = φ_RL − φ and Δφ_LR = φ_LR − φ.
  - Entanglement exists unless Δφ_LR + Δφ_RL = 2nπ; it is maximal for π.
- **Explicit scheme:**
  - Micro-diamonds with an NV spin, m ~ 10⁻¹⁴ kg, R ~ 1 µm, in traps d ~ 450 µm apart.
  - Stern–Gerlach splitting with ∂ₓB ~ 10⁶ T/m for τ_acc ~ 500 ms gives Δx ~ 250 µm, so d − Δx ~ 200 µm.
  - The Casimir–Polder potential is ~0.1 of the gravitational potential (diamond ε ≈ 5.7). The masses then fall freely for τ ~ 2.5 s.
  - The mutual acceleration Gm/(d − Δx)² ~ 10⁻¹⁶ m/s² is "truly negligible".
- **Reported phases:** "Δφ_LR ∼ −0.2 and Δφ_RL ∼ 0.7", giving "W ∼ 1.16".
  - The witness is printed (arXiv v1) as `W = |⟨σx ⊗ σz⟩ − ⟨σy ⊗ σz⟩|`, with W > 1 taken to mean entanglement.
- **Decoherence budget:** P = 10⁻¹⁵ Pa and T = 0.15 K. Electronic spin coherence of 1 s below 77 K. Collisional decoherence time is comparable to the total fall time τ + 2τ_acc ~ 3.5 s. Centre-of-mass temperatures up to 100 K are tolerable.
- **Companion proposal:** Marletto & Vedral (2017). The Overstreet et al. companion paper calls both "BMV experiments" and says their theoretical implication "is still an open question" **[contested interpretation]**. — [arXiv:2209.02214](https://arxiv.org/abs/2209.02214)

**Progress to 2026 [mainstream; no GIE result yet]:**
- A 2026 search found no completed gravitationally-induced-entanglement experiment. Enabling results cited in recent preprints:
  - M. Rossi, A. Militaru et al., "Quantum delocalization of a levitated nanoparticle", PRL 135, 083601 (2025).
  - M. Kamba, N. Hara, K. Aikawa, "Quantum squeezing of a levitated nanomechanical oscillator", Science 389, 1225 (2025).
  - — [reference list as seen via arXiv:2602.19306](https://arxiv.org/pdf/2602.19306)
- Diamagnetic microchip traps for Stern–Gerlach µm superpositions of two nanodiamonds, with EM screening. — [arXiv:2411.02325](https://arxiv.org/abs/2411.02325); [PRA](https://journals.aps.org/pra/abstract/10.1103/qdyh-rbyq)
- 2026 theory preprints on GIE with mechanically constrained systems (e.g. a carbon-nanotube pendulum) and on shielded setups. — [arXiv:2605.00967](https://arxiv.org/pdf/2605.00967); [arXiv:2604.22593](https://arxiv.org/pdf/2604.22593)
- **Present source-mass frontier:** 92.1 mg (Westphal 2021); 0.43 mg test mass (Fuchs 2024). See Q2.
- **State of the art for superpositions** (cited in Donadi 2021):
  - 0.5 m spatial superposition, but of Rb atoms (1.42×10⁻²⁵ kg).
  - Macromolecules above 25 kDa (~10⁻²³ kg).
  - Phonon superpositions of about 10¹⁶ C atoms over 10⁻¹¹ m, lasting ~10⁻¹² s.
  - — [Donadi et al. arXiv:2111.13490](https://arxiv.org/abs/2111.13490)

**Diósi–Penrose test: Donadi, Piscicchia, Curceanu, Diósi, Laubenstein, Bassi, "Underground test of gravity-related wave function collapse", arXiv:2111.13490 (Nature Physics 2021) [refuted: parameter-free DP]** — [arXiv:2111.13490](https://arxiv.org/abs/2111.13490)
- **Collapse time:** `τ_DP = ħ/ΔE_DP`, with `ΔE_DP(d) = −8πG ∫dr ∫dr′ µ(r)[µ(r′ + d) − µ(r′)]/|r − r′|`, normalisation as printed (eqs. 1–2).
  - Point masses diverge, so the mass density is smeared with size R₀.
- **Radiation from DP-induced diffusion of protons:** `dΓ/dω = 2 G e² N² N_a / (3 π^(3/2) ε₀ c³ R₀³ ω)` (eq. 4). N is the atomic number and N_a the number of atoms; electrons are conservatively neglected.
- **Experiment:**
  - HPGe coaxial p-type detector at LNGS: 8.0 cm × 8.0 cm crystal, 375 cm³ active volume, shielded with 62 kg of electrolytic copper and other layers.
  - About 62 days of data (August 2014 and August 2015). The analysis region is 1000–3800 keV; a ≈ 1.8×10⁻²⁹ m³ appears in the likelihood.
- **Results:**
  - **R₀ > 0.54×10⁻¹⁰ m (95%)**, about 3 orders of magnitude above previous bounds.
  - Previous bounds: R₀ ≥ (40.1 ± 0.5)×10⁻¹⁵ m from GW detectors and ≳ 10⁻¹³ m from neutron stars.
  - Penrose's "natural" R₀ for Ge (from the Schrödinger–Newton solution) is 0.05×10⁻¹⁰ m. The paper concludes the natural parameter-free DP model is ruled out.
- **Penrose's own mirror proposal:** a ~10⁻¹² kg mirror would have τ_DP ≃ 0.002–0.013 s.

**Single gravitons: Tobar, Manikandan, Beitel, Pikovski, "Detecting single gravitons with quantum sensing", arXiv:2308.15440v2 (Nature Communications 2024) [mainstream proposal]** — [arXiv:2308.15440](https://arxiv.org/abs/2308.15440)
- **Spontaneous emission** (n = 1 → 0) for a bar resonator, eq. 4:
  - `Γ_spon = 8 G M L² ω_l⁴/(l⁴ π⁴ c⁵) = 8π G ρ v_s⁴ R²/(L c⁵)`, with v_s = Lω_l/(lπ).
  - Nb bar (ρ = 8570 kg/m³, v_s ≈ 0.5×10⁴ m/s, L = 1 m, R = 0.5 m): **Γ_spon ≈ 10⁻³³ Hz**. The atomic 3d → 1s rate is about 10⁻⁴⁰ Hz.
- **Stimulated absorption or emission**, eq. 5:
  - `Γ_stim = M L² ω_l² h²/(4 l⁴ π⁵ ħ) = M v_s² h²/(4 l² π³ ħ)`.
  - "Using an aluminum bar of mass 1800 kg, a gravitational wave with amplitude h = 5 × 10⁻²² will result in **Γ_stim ≈ 1 Hz**."
- **Optimal mass for a chirping binary** (eq. A45): `M ≈ (24π²/5) · ħ/(h₀² v_s²) · (G M_c/(2c³))^(5/3) · ω^(8/3)`.
  - Beryllium at ω = 2π × 100 Hz with h₀ ~ 2×10⁻²² (GW170817): **M ~ 15 kg**.
- **Table I and protocol:**
  - Other entries: about 250 kg (Al), about 9 t and about 6 t (Nb, CuAl6%), more than 52 t, more than 20 t, and about 100 kg (quartz, 5.5 MHz, h₀ ~ 10⁻¹⁶, hypothetical UHF source).
  - Temperatures from 1 mK down to 0.1 µK. Q factors 10⁷–10¹³.
  - Protocol: continuous quantum-jump sensing of the first excited level, correlated with LIGO events. A simulated example uses M = 21.73 kg, ω/2π = 100 Hz, h₀ = 2×10⁻²² and a 4 s duration.
- "Currently detected GWs consist of ≳ 10³⁶ gravitons" (Appendix A3).
- **Interpretation [contested]:**
  - Carney, Domcke and Rodd (2023) agree such a detector is possible in Dyson's sense. However, "the same signal could just as well be explained via classical gravitational waves", so it would not prove quantisation. — [arXiv:2308.12988](https://arxiv.org/abs/2308.12988)
  - The Tobar authors frame it as the analogue of the photo-electric effect. — [arXiv:2308.15440](https://arxiv.org/abs/2308.15440)

**Dyson 2013, "Is a graviton detectable?" (Poincaré Prize Lecture, ICMP 2012; IAS 2013 essay) [mainstream conjecture]:**
- To detect a single graviton with LIGO-type apparatus, the mirrors would have to be so heavy that they collapse into a black hole. Dyson conjectures that single gravitons may be unobservable by any conceivable apparatus. — [IAS essay](https://www.ias.edu/ideas/2013/dyson-graviton); [abstract PDF (fetch refused; content from search)](https://icmp12.aau.dk/media/15234/abstract-dyson.pdf)
- A 2026 lecture restates the step: a single graviton's arm-length change is about the Planck length, and the uncertainty principle then forces black-hole collapse. — [Soda 2026 slides](https://indico.ipmu.jp/event/519/contributions/9084/attachments/5337/7334/QuantumCrossRoad.pdf)

### Inferences
- **BMV phases (computed):** m = 10⁻¹⁴ kg, d = 450 µm, Δx = 250 µm, τ = 2.5 s.
  - φ = 0.352, φ_RL = 0.791, φ_LR = 0.226 rad, so Δφ_RL = 0.440 and Δφ_LR = −0.126. The paper's quoted values are about 0.7 and −0.2, about 1.6× larger.
  - The mutual acceleration is 1.67×10⁻¹⁷ m/s² (paper: "~10⁻¹⁶").
- **Witness check (computed; important for implementers):**
  - The printed operator |⟨σx⊗σz⟩ − ⟨σy⊗σz⟩| = |(x₁ − y₁)z₂| can reach √2 on product states. So W > 1 does not certify entanglement for that operator.
  - The operator |⟨σx⊗σz⟩ − ⟨σy⊗σy⟩| is bounded by 1 on separable states (Cauchy–Schwarz). Its optimum over local bases for this state is **1 + |sin((Δφ_LR + Δφ_RL)/2)|** (numerically verified).
  - That gives 1.156 with my computed phases, matching the published "W ∼ 1.16", and 1.247 with the paper's rounded phases.
  - Likely the intended witness is the σy⊗σy form with locally readjusted bases, as the paper notes. Check the published PRL version before coding.
- **Diósi–Penrose test numbers (computed):**
  - Eq. 4 for the Ge crystal (N_a = 1.66×10²⁵, N = 32), with no efficiency or background: ≈ 11 photons/day in 1–3.8 MeV at R₀ = 0.54×10⁻¹⁰ m. At the "natural" R₀ = 0.05×10⁻¹⁰ m it is ≈ 1.3×10⁴/day. The R₀⁻³ scaling explains why a low-background crystal is decisive.
  - Uniform-sphere Penrose estimate, using E_Δ = GM²(6/(5R) − 1/d) for d ≥ 2R (normalisation differs from Donadi's printed −8πG form; check it): a 10⁻¹⁴ kg diamond-scale sphere (R = 1 µm, d = 10R) gives τ ≈ 14 ms. A 1 mg tungsten sphere gives about 3×10⁻²¹ s.
- **Tobar checks (computed):**
  - Γ_stim(1800 kg Al, v_s = 5.4 km/s, h = 5×10⁻²²) = 1.00 Hz.
  - Γ_spon(Nb bar) = 9.3×10⁻³⁴ Hz.
  - Be optimal mass (v_s = 13 km/s, M_c = 1.188 M_⊙, 100 Hz, h₀ = 2×10⁻²²) = 12.8 kg (paper ~15 kg).
  - Graviton energy at 100 Hz: ħω = 6.6×10⁻³² J.
- **Dyson-scale numbers (computed):** a 150 Hz GW with h = 10⁻²¹ carries `F = c³ω²h²/(32πG)` = 3.6×10⁻³ W/m². That is about 3.6×10²⁷ gravitons per m² in 0.1 s, consistent with Tobar's ≳ 10³⁶ per detection once the effective area (~λ²) is included.

### Gaps
- I did not fetch the Marletto–Vedral PRL, the published BMV PRL (to check the witness typo), or the Nature Communications version of Tobar.
- I found no source for Penrose's response to Donadi (whether his own proposal implies the radiation). Treat "Penrose's model refuted" as applying to the parameter-free DP master-equation version only.
- No dedicated 2025–2026 GIE experimental result was found; the Rossi and Kamba references come from a preprint reference list seen in search, not read directly.

---

## Q7. Superconductors as gravitational-wave detectors or transducers: Chiao's claims (contested), the experiment, and the standard rebuttals, with numbers

### Takeaway
Chiao's claims that superconductors convert EM waves into GWs with efficiency of order unity, or act as GW mirrors, rest on heuristic arguments **[contested]**: (c/v_s)⁴ ≈ 10²⁰ enhancement, "impedance matching" to Z_g = 4πG/c ≈ 2.8×10⁻¹⁸ m²/(kg s), and a "Heisenberg–Coulomb" effect.
- The only experiment (YBCO at 12 GHz, 77 K) saw nothing: η < 1.6×10⁻⁵.
- Kiefer & Weber (2005) judged the arguments heuristic and inapplicable to the classical Weinberg efficiency formula (η ≈ 10⁻³⁴ for bars). Licht's classical calculation found no strong enhancement.
- No positive evidence exists to 2026. Later supporters (gravitational Casimir effect, neutron-star mirrors) are also theory only.

### Cited Findings

**Chiao's claims [contested]:**
- "Superconductors as quantum transducers and antennas for gravitational and electromagnetic radiation". Ginzburg–Landau arguments give conversion efficiency "on the order of unity", with EM and GW fields impedance-matched about a hundred coherence lengths beneath the surface. — [arXiv:gr-qc/0204012](https://arxiv.org/abs/gr-qc/0204012)
- Same claim plus a Meissner-like expulsion of the Lense–Thirring field. — [arXiv:gr-qc/0208024](https://arxiv.org/abs/gr-qc/0208024)
- Thin type-I superconducting films reflect microwave GWs via a "Heisenberg–Coulomb" effect. — [Minter, Wegter-McNelly, Chiao, arXiv:0903.0661](https://arxiv.org/abs/0903.0661); [Chiao, Minter, Wegter-McNelly, arXiv:0903.3280](https://arxiv.org/abs/0903.3280)
- Electron-coated superfluid-helium "Millikan oil drops" as Planck-mass transducers. — [arXiv:gr-qc/0702100](https://arxiv.org/abs/gr-qc/0702100); [arXiv:0904.3956](https://arxiv.org/abs/0904.3956); [IJMPD 16, 2309 (2008), arXiv:gr-qc/0606118](https://arxiv.org/abs/gr-qc/0606118)
- Laser-like GW generation from dynamical Casimir in SRF cavities. — [arXiv:1712.08680](https://arxiv.org/abs/1712.08680)
- Press coverage (2009) reported the claimed effect as about 42 orders of magnitude larger than in an ordinary conducting sheet, and asked why sensitive experiments had not seen it. — [MIT Technology Review](https://www.technologyreview.com/2009/03/23/214718/if-superconducting-sheets-reflected-gravitational-waves/amp/)

**The experiment [mainstream measurement, null]:**
- Chiao, Fitelson and Speliotopoulos (2003) ran a first Hertz-type experiment with YBCO as generator and detector. Upper limit on conversion efficiency: **1.6×10⁻⁵** at liquid-nitrogen temperature. — [arXiv:gr-qc/0304026](https://arxiv.org/abs/gr-qc/0304026)
- Per Kiefer & Weber: 12 GHz EM quadrupole radiation, two Faraday cages; the result assumes equal conversion in both directions. Chiao blamed residual microwave and far-infrared losses or ohmic dissipation. — [Kiefer & Weber, arXiv:gr-qc/0408010](https://arxiv.org/abs/gr-qc/0408010)

**Standard rebuttal: Kiefer & Weber, "On the interaction of mesoscopic quantum systems with gravity", Ann. Phys. (Leipzig) 14, 253 (2005) [mainstream critique]** — [arXiv:gr-qc/0408010](https://arxiv.org/abs/gr-qc/0408010)
- **Weinberg bar emission efficiency:** `η_rad = Γ_grav/Γ_heat = 64 G M v_s⁴/(15 L² c⁵ Γ_heat) ≈ 10⁻³⁴` (eq. 37).
- **Chiao's heuristic:** replace v_s by c, giving (c/v_s)⁴ ≈ 10²⁰ for v_s ≈ 10⁻⁵c. Relevant frequencies are about 10 MHz and above, where "no detectable astrophysical sources are available".
- **Impedances:** EM, Z₀ = √(µ₀/ε₀) ≈ 377 Ω. GW, `Z_g = √(µ_g/ε_g) = 4πG/c ≈ 2.79×10⁻¹⁸ m²/(kg s)` (eqs. 38–39).
- **Verdict:** the impedance-matching "expectations are, however, based on heuristic considerations and not on detailed calculations". Weinberg's formula rests on classical sound-wave properties, and "it is not clear whether these arguments can be applied to quantum fluids". There is no non-relativistic derivation that includes both EM and GW fields.
- **Licht's classical Ginzburg–Landau calculation:** EM generation by GWs in superconductors is "not strongly enhanced compared to that in normal conductors and should not lead to detectable results".
- **Recommendation:** "the interaction first be further explored on a theoretical basis before new experiments are carried out."

**Later supporting theory (no experiment) [contested]:**
- Gravitational Casimir effect from partial GW reflection by superconductors. — [Bahamonde, Faizal, Quach, Norte, IJMPD 29 (2020), arXiv:2005.06640](https://arxiv.org/abs/2005.06640)
- Neutron stars as GW mirrors. — [Wei et al., ApSS 365, 148 (2020), arXiv:1911.04201](https://arxiv.org/abs/1911.04201)
- The mainstream UHF-GW review mentions the superconducting mirror and Anandan–Chiao detectors only as proposals. — [Aggarwal et al. arXiv:2011.12414](https://arxiv.org/abs/2011.12414)

### Inferences
- **Size of the claim (computed):**
  - The Coulomb-to-gravity force ratio between two electrons is e²/(4πε₀ G m_e²) = **4.17×10⁴²**. This matches the "42 orders of magnitude" enhancement in the press account. The Heisenberg–Coulomb argument effectively proposes that the electrical stiffness of the Cooper-pair fluid, rather than its mass, sets the GW response.
  - The equivalence principle makes electrons and ions fall together under a GW. Kiefer & Weber and Licht find no mechanism that converts the 10⁴² ratio into coupling.
- **Weinberg check (computed):** η for M = 1000 kg, v_s = 5 km/s, L = 1.5 m and Γ_heat = ω/Q (1 kHz, Q = 10⁶) is 5×10⁻³³, the same order as Kiefer & Weber's "≈ 10⁻³⁴". Multiplying by 10²⁰ (Chiao's heuristic) still gives η ~ 10⁻¹³. That is far below unity and below the 1.6×10⁻⁵ experimental limit, which is therefore non-constraining for the conservative estimate.
- **Z_g check (computed):** 4πG/c = 2.798×10⁻¹⁸ m²/(kg s), matching Kiefer & Weber.
- **Status:** order-unity transduction and GW mirrors remain **[contested]**, with no positive measurement in 23 years. The specific YBCO 77 K test set only η < 1.6×10⁻⁵. Related superconductor "gravity" anomalies (Podkletnov shielding; Tajmar's 10⁻⁴ω gravitomagnetic London moment) are **[refuted / not reproduced]**; see the sister note.

### Gaps
- I found no peer-reviewed comment written specifically against Minter–Wegter-McNelly–Chiao (2009/2010) beyond Kiefer & Weber (2005, which predates it) and Licht. A targeted citation search (INSPIRE) is needed.
- The 12 GHz figure is from Kiefer & Weber's description; I did not read the Chiao–Fitelson–Speliotopoulos experimental details (geometry, power, noise) directly.
