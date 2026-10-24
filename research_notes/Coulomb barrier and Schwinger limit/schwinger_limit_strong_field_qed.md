# The Schwinger limit and strong-field QED: formulas, experiments and facilities

Conventions in these notes:
- **Units.** SI unless stated. Di Piazza et al. and Melissinos use Gaussian or natural units; conversions are given where they matter. Constants: CODATA 2022, with e, h and c exact.
- **Status labels.** **[mainstream]** is standard, peer-reviewed QED. **[contested]** is published but disputed in its interpretation, or preliminary. **[refuted]** has been contradicted by later experiments.
- **"Computed".** These values were recomputed with `python3 -I sfqed_checks.py`, a stdlib-only script in the session scratch directory (`/tmp/claude-1000/-home-dan-code-FusionVault-quantum-gravimeter/7beb35b3-c723-440c-9450-66a991641392/scratchpad/sfqed/sfqed_checks.py`). "Ratio" means computed ÷ published. Every "computed" line is the note author's own check, not a published claim.
- **Reading.** Primary PDFs were read in full text (pdftotext) for Di Piazza 2012, the 2025 ESPPU strong-field-QED input, the LUXE 2023 proceedings, Schützhold–Gies–Dunne 2008, Berdyugin 2022, Melissinos 1998, Rafelski 2016 and Maltsev 2019. The other sources are arXiv abstract pages or search-result snippets, marked "(snippet)". Paywalled journal pages were not read.

---

## Q1. What are the critical fields E_S, B_S and I_S, and the Schwinger pair-production rate? Worked values at E = 0.1, 0.05 and 0.01 E_S

### Takeaway
E_S = m²c³/(eħ) = 1.3233×10¹⁸ V/m, B_S = 4.4140×10⁹ T, and I_S is 4.65×10²⁹ W/cm² or 2.32×10²⁹ W/cm² depending on the intensity convention. The Schwinger series gives the vacuum-decay rate; its n = 1 term alone is the mean pair-production rate. Both carry the exp(−πE_S/E) factor. That factor makes the rate fall from 2.5×10⁴⁰ m⁻³s⁻¹ at 0.1 E_S, to 1.4×10²⁶ at 0.05 E_S, to 4×10⁻⁸⁵ at 0.01 E_S (all computed).

### Cited Findings
**Critical field [mainstream]**
- Di Piazza et al. (RMP 84, 1177, 2012) define F_cr = m²/|e| (natural units) as "1.3 × 10¹⁶ V/cm". Their symbol table defines I0 = E0²/4π and I_cr = F_cr²/4π, "≈ 4.6 × 10²⁹ W/cm²". The arXiv text, as extracted, prints the magnetic equivalent as "4.4 × 10¹⁴ G". — [Di Piazza et al., arXiv:1111.3886](https://arxiv.org/abs/1111.3886)
  - Computed: B_S = 4.414×10¹³ G, so the printed 10¹⁴ G looks like a typo, or an extraction artefact, for 10¹³ G.
- Physical meaning: a field of order F_cr gives an e⁺e⁻ pair energy of about 2m over one Compton wavelength λ_C (Heisenberg–Euler 1936; Sauter 1931; Schwinger 1951). — [Di Piazza et al.](https://arxiv.org/abs/1111.3886)
- Schützhold–Gies–Dunne write the pair rate as suppressed like exp{−πE_S/E}, with E_S = m²c³/(ħq). — [arXiv:0807.0754](https://arxiv.org/abs/0807.0754) (PRL 101, 130404, 2008)
- Melissinos gives E_c = m²c³/(eħ) = 1.3×10¹⁶ V/cm. He states that reaching Υ = 1 with no incident particle "would require an intensity of the em flux I = 5 × 10²⁹ W/cm²". — [Melissinos, hep-ph/9805507](https://arxiv.org/abs/hep-ph/9805507)

**Formulas and test vectors (SI)**

| Quantity | Exact form | Computed (CODATA 2022) | Published test value | Ratio |
|---|---|---|---|---|
| E_S | m_e²c³/(eħ) | 1.32329×10¹⁸ V/m | 1.3×10¹⁶ V/cm (Di Piazza; Melissinos) | 1.018 (rounding) |
| B_S | m_e²c²/(eħ) = E_S/c | 4.41401×10⁹ T = 4.414×10¹³ G | 4.4×10¹³ G (correct order) | 1.003 |
| I_S, "E²/4π" convention | ε0 c E_S² | 4.648×10²⁹ W/cm² | 4.6×10²⁹ W/cm² (Di Piazza I_cr = F_cr²/4π) | 1.010 |
| I_S, cycle-averaged linear | ε0 c E_S²/2 | 2.324×10²⁹ W/cm² | (no published value read) | — |
| ƛ_C | ħ/(m_e c) | 3.86159×10⁻¹³ m | — | — |
| check | eE_S ƛ_C / m_ec² | 1 exactly | — | — |

**Schwinger rate [mainstream]**
- Cohen & McGady say Schwinger's result for the rate w is "an infinite series", while the pair-production rate Γ "is just the first term of that series". — [Cohen & McGady, arXiv:0807.1117](https://arxiv.org/abs/0807.1117)
  - w is the vacuum-decay rate, 2 Im L_eff/ħ. Γ is Nikishov's mean number of pairs per volume per time.
- Melissinos (after Brezin–Itzykson) gives the pair probability per unit time and volume for an oscillating field: w ≃ (αE²/πħ) exp(−πm²c³/(eħE)) for η ≫ 1, and w ∝ (eE/2mω0c)^(4mc²/ħω0) for η ≪ 1 (the multiphoton regime). His interpolating form is w = (αE²/πħ)·[g + ½g′/η]⁻¹·exp(−πg(η)E_c/E), with g(η) = (4/π)∫₀¹ dy √[(1−y²)/(1+y²/η²)].
  - Published numerical prefactor (his eq. 8′): "αE_c²/πħ = 3.45 × 10⁵⁰ cm⁻³ s⁻¹". — [hep-ph/9805507](https://arxiv.org/abs/hep-ph/9805507)

Constant-field formulas in SI, with x = E/E_S:
- **Spinor QED, vacuum decay:** w(E) = (e²E²/(4π³ħ²c)) Σₙ₌₁^∞ n⁻² exp(−nπ/x) = (αε0E²/(π²ħ)) Σₙ n⁻² exp(−nπ/x).
- **Spinor QED, pair rate (n = 1 term):** Γ(E) = (αε0E²/(π²ħ)) exp(−π/x) = (c/(4π³ƛ_C⁴)) x² exp(−π/x).
- **Scalar QED:** half the prefactor and alternating signs, (e²E²/(8π³ħ²c)) Σ (−1)ⁿ⁺¹ n⁻² exp(−nπ/x). Computed scalar/spinor ratio at x = 0.1: 0.5.
- **Prefactor at x = 1:** c/(4π³ƛ_C⁴) = 1.08704×10⁵⁶ m⁻³s⁻¹ = 1.087×10⁵⁰ cm⁻³s⁻¹ (computed).
- **Cross-check against Melissinos:** his αE_c²/(πħ), with E in Heaviside–Lorentz normalisation (αε0E_S²/(πħ)), computes to 3.415×10⁵⁰ cm⁻³s⁻¹ against the published 3.45×10⁵⁰ (ratio 0.990). The Schwinger constant-field prefactor is exactly this number divided by π.

**Worked rates (computed)**

| E/E_S | exp(−π/x) | Γ_pair (m⁻³ s⁻¹) | Γ_pair (cm⁻³ s⁻¹) | w_sum/Γ_pair | pairs in (1 µm)³ during 1 µm/c |
|---|---|---|---|---|---|
| 1 | 4.321×10⁻² | 4.698×10⁵⁴ | 4.698×10⁴⁸ | 1.0110 | 1.6×10²² |
| 0.5 | 1.867×10⁻³ | 5.075×10⁵² | 5.075×10⁴⁶ | 1.00047 | 1.7×10²⁰ |
| 0.2 | 1.507×10⁻⁷ | 6.553×10⁴⁷ | 6.553×10⁴¹ | 1.0000 | 2.2×10¹⁵ |
| **0.1** | 2.271×10⁻¹⁴ | **2.469×10⁴⁰** | 2.469×10³⁴ | 1.0000 | 8.2×10⁷ |
| **0.05** | 5.158×10⁻²⁸ | **1.402×10²⁶** | 1.402×10²⁰ | 1.0000 | 4.7×10⁻⁷ |
| 0.02 | 6.042×10⁻⁶⁹ | 2.627×10⁻¹⁶ | 2.6×10⁻²² | 1.0000 | 8.8×10⁻⁴⁹ |
| **0.01** | 3.651×10⁻¹³⁷ | **3.968×10⁻⁸⁵** | 4.0×10⁻⁹¹ | 1.0000 | 1.3×10⁻¹¹⁷ |

- Computed: one pair in a (1 µm)³ × (1 µm/c) volume needs E = 0.0644 E_S, which corresponds to I = ε0cE²/2 = 9.6×10²⁶ W/cm² for a static field of that size.
- Rust test vectors: exp(−π) = 4.3214×10⁻² at x = 1. w_sum/Γ_pair = 1.01102 at x = 1 (the n ≥ 2 terms add 1.1%) and is 1 to better than 10⁻⁶ for x ≤ 0.2.

### Inferences
- Below about 0.3 E_S the full sum and the first term agree to better than 10⁻⁴. A calculator only needs the series near x ≈ 1. Label the two outputs clearly: "vacuum decay" for the sum, "mean pairs" for the n = 1 term.
- The intensity convention changes the reported I_S by a factor of 2. Always store E and derive I with an explicit convention flag: linear, cycle-averaged ε0cE²/2, or "E²/4π" (ε0cE²).
- The table describes a uniform static E field. A plane wave alone has F = G = 0 and produces no pairs, as Melissinos notes. Real laser pairs at sub-critical fields come from colliding beams, multi-pulse foci, or an incident high-energy particle (Q2–Q3).
- Standard results not re-derived from a fetched source in this session:
  - For general fields the rate uses the invariants ε, η (with ε² − η² = E² − c²B², εη = cE·B). The Nikishov/Schwinger general form is Γ = (e²εη/(4π²ħ²c)) coth(πη/ε) exp(−πE_S/ε), with the full series having 1/n coth(nπη/ε) terms. Check the factors against Fedotov et al. 2023 §2 before coding it.

### Gaps
- No published table of the Schwinger rate at 0.1, 0.05 and 0.01 E_S was found. These rows are computed only. The prefactor is tied to a published number only through Melissinos's 3.45×10⁵⁰ cm⁻³s⁻¹.
- Schwinger 1951 (Phys. Rev. 82, 664) and Fedotov et al. 2023 (too large for the fetch tool, >10 MB) were not read in full text.

---

## Q2. Laser parameters a₀, χ_e and χ_γ, the Ritus–Narozhny conjecture, record intensity, and planned facilities with their a₀ and χ

### Takeaway
a₀ = eE₀/(mcω) and χ = E_rest/E_S. For a head-on or angled collision, χ_e = (1 + β cos θ)·γ·a₀·ħω/(mc²). This formula reproduces the LUXE design values to within 1.5%.
- **Record:** the intensity record is still CoReLS, (1.1 ± 0.2)×10²³ W/cm² in 2021, which is a₀ ≈ 230 (linear convention, 0.8 µm assumed). That is still 2×10⁶ below I_S.
- **Facilities:** PW-class lasers give a₀ ≈ 100–300 (ESPPU 2025 table). χ ≫ 1 needs a GeV-class electron beam: χ ≈ 0.12 per unit a₀ at 10 GeV.
- **Ritus–Narozhny:** the Ritus–Narozhny regime αχ^(2/3) ≈ 1 means χ ≈ 1600, far beyond any planned facility. The conjecture's meaning is **[contested]**.

### Cited Findings
**Definitions [mainstream]**
- **Di Piazza's ξ0:** "ξ0 = 6.0 √(I0[10²⁰ W/cm²]) λ0[µm] = 7.5 √(I0[10²⁰ W/cm²])/ω0[eV], where I0 = E0²/4π is the wave's peak intensity".
- **Di Piazza's χ0:** χ0 = (|e|/m³)√(−(F0,μν p0^ν)²) = (p0,−/m)(E0/F_cr).
  - For an ultrarelativistic electron counter-propagating with the wave: "χ0 = 5.9 × 10⁻² ε0[GeV] √(I0[10²⁰ W/cm²])".
  - The photon analogue is κ0 = 5.9×10⁻² ω[GeV] √(I0[10²⁰ W/cm²]).
  - Source: [Di Piazza et al.](https://arxiv.org/abs/1111.3886)
- The ESPPU 2025 input defines χ = F_rf/F_qed, the field in the rest frame of the charge in units of the Schwinger field. In a plane wave "χ = ξη", with η = k·p/m². In the quasi-static regime (ξ ≫ 1, χ ≪ 1) the nonlinear Breit–Wheeler probability "scales as ∼ exp(−8/3χ)". — [Sarri et al., ESPPU input, arXiv:2504.02608](https://arxiv.org/abs/2504.02608)
- Classical radiation reaction (Landau–Lifshitz energy loss): Δγ/γ = √(π/2) τ0 τl ωl² γ ξ² / (1 + √(π/2) τ0 τl ωl² γ ξ²), with τ0 = 2e²/3m "≃ 6.4 × 10⁻²⁴ s".
  - The extracted text drops a superscript in the numerator. ωl² is required for the expression to be dimensionless.
  - Source: [ESPPU input](https://arxiv.org/abs/2504.02608)
- Gaunt-factor approximations for quantum radiation reaction: g(χ) ≈ (1 + 12χ + 31χ² + 3.7χ³)^(−4/9), and g(χ) ≈ (1 + 4.8(1+χ)ln(1+1.7χ) + 2.44χ²)^(−2/3). — [ESPPU input](https://arxiv.org/abs/2504.02608)

**Computed checks of the definitions**

| Check | Computed | Published | Ratio |
|---|---|---|---|
| a₀ coefficient, linear, I = ε0cE₀²/2 | 0.85493 λ[µm]√(I/10¹⁸) | 0.855 (common rule) | 1.000 |
| a₀ coefficient, I = ε0cE₀² ("E²/4π") | 6.045 λ[µm]√(I/10²⁰) | 6.0 (Di Piazza) | 1.0075 |
| ħω at 800 nm | 1.5498 eV | 1.55 eV (LUXE) | 1.000 |
| Di Piazza χ0 coefficient (1 GeV, 10²⁰ W/cm², 800 nm, head-on) | 0.0574 | 0.059 | 0.973 |
| τ0 = 2r_e/3c | 6.266×10⁻²⁴ s | 6.4×10⁻²⁴ s (ESPPU) | 0.979 |
| g(χ=1), approximation 1 / approximation 2 | 0.1795 / 0.1811 | — | — |
| exp(−8/3χ) at χ = 0.1 / 0.3 / 1 | 2.6×10⁻¹² / 1.4×10⁻⁴ / 6.9×10⁻² | — | — |

**Ritus–Narozhny conjecture [contested]**
- **Statement:** in fields that can be treated as locally constant and crossed, the contribution of a class of dressed loop diagrams "scales as αχ^(2/3)". When this reaches order one, "the perturbative expansion in the quantum photon field breaks down".
- **Evidence for:** bubble-chain corrections to electron forward scattering have been confirmed to scale as M_{n+1}/M_n ∼ αχ^(2/3). The resummed bubble chain "remains unbounded as χ grows", so it cannot be the only relevant correction at large χ.
- **Evidence against:** other fully quantum examples show the Furry expansion converging. In a constant magnetic field the all-loop Heisenberg–Euler limit has a different field dependence.
- **Status:** the ESPPU input concludes "it is simply not known how QED behaves there".
- Source: [ESPPU input](https://arxiv.org/abs/2504.02608)
- Mironov, Meuren & Fedotov (PRD 102, 053005, 2020) resummed bubble-type polarisation corrections to the electron self-energy in a constant crossed field. They found that perturbation theory breaks down already at g = αχ^(2/3) ≳ 1, while the leading-order result stays dominant until g ≫ 1 (snippet). — [arXiv:2003.06909](https://arxiv.org/abs/2003.06909)
- Original statement: Fedotov, "Conjecture of perturbative QED breakdown at αχ^(2/3) ≳ 1", J. Phys. Conf. Ser. 826, 012027 (2017); Ritus, Sov. Phys. JETP (1970). Both are listed as references in the [ESPPU input](https://arxiv.org/abs/2504.02608); neither was read.
- Computed values:
  - αχ^(2/3) = 1 at χ = α^(−3/2) = 1604.2.
  - Reaching χ = 1600 head-on with a₀ = 100 at 800 nm needs a 1.35 TeV electron.
  - At LUXE phase-1 (χ = 4.45), αχ^(2/3) = 0.020.

**Record intensity [mainstream]**
- CoReLS (IBS, Korea), Yoon et al., "Realization of laser intensity over 10²³ W/cm²", Optica 8, 630 (2021), DOI 10.1364/OPTICA.420520.
  - Peak intensity of (1.1 ± 0.2)×10²³ W/cm² over 80 consecutive shots at 0.1 Hz.
  - f/1.1 off-axis parabola with two-stage adaptive optics; spot 1.1 µm FWHM.
  - A CoReLS presentation lists 2.7 PW and 55.6 J for this shot (snippet).
  - Sources: [IBS repository](https://pr.ibs.re.kr/handle/8788114/10375) (snippet); [Optica news release](https://www.optica.org/en-us/about/newsroom/news_releases/2021/researchers_produce_laser_pulses_with_record-break/) (snippet)
- A search for a newer intensity record (2025–2026) found none. ELI's L4 ATON news is a peak-power record, not an intensity record (snippet). — [ELI news](https://www.eli-laser.eu/news/eli-s-l4-aton-laser-achieves-5-petawatt/)
- Computed for 1.1×10²³ W/cm² at an assumed 0.8 µm: a₀ = 227 (linear convention) or 160 (E²/4π convention). E₀/E_S = 6.9×10⁻⁴. I_S(ε0cE²/2)/I = 2.1×10⁶.

**Facilities (status as of 2025–2026)**
- **ESPPU 2025 Table 1** ("indicative" ξ, for a 2 µm FWHM spot holding 60% of the energy). Columns: E (J), τ (fs), ξ, repetition rate. Designs are marked in italics in the original table. — [ESPPU input](https://arxiv.org/abs/2504.02608)
  - FACET-II: 0.3 J, 30 fs, ξ = 6, 10 Hz
  - GEMINI: 15 J, 45 fs, ξ = 32, 1/20 Hz
  - BELLA: 40 J, 40 fs, ξ = 55, 1 Hz
  - EPAC: 30 J, 30 fs, ξ = 55, 10 Hz
  - CALA: 70 J, 27 fs, ξ = 90, 1 Hz
  - ZEUS: 75 J, 25 fs, ξ = 100, 1/60 Hz
  - APOLLON: 70 J, 22 fs, ξ = 100, 1/60 Hz
  - CoReLS: 80 J, 20 fs, ξ = 110, 1/10 Hz
  - ELI-NP: 200 J, 20 fs, ξ = 180, 1/60 Hz
  - SULF: 300 J, 23 fs, ξ = 200, 1/60 Hz
  - ELI-BL: 1500 J, 150 fs, ξ = 230, 1/60 Hz
  - NSF OPAL: 2×500 J, 20 fs, ξ = 280, 1/300 Hz
  - VULCAN2020: 400 J, 20 fs, ξ = 330, 1/300 Hz
  - SEL: 25,000 J, 25 fs, ξ = 1,800
- The same input says that a 5 GeV electron beam meeting 10²³ W/cm² head-on "would result in a quantum parameter χ ≈ 10". It also notes that LWFA electron beams above 10 GeV have been reported with a 500 TW driver. — [ESPPU input](https://arxiv.org/abs/2504.02608)
- **ZEUS** (University of Michigan): first official experiment at 2 PW, announced May 2025, with plans to reach 3 PW (snippet). — [Michigan Engineering news](https://news.engin.umich.edu/2025/05/the-us-has-a-new-most-powerful-laser)
- **ELI-Beamlines L4 ATON:** 5.1 PW (786 J compressed to 154 fs) in a September 2025 test campaign, with plans "to increase power toward 10 PW in 2026". — [ELI news, 16 Oct 2025](https://www.eli-laser.eu/news/eli-s-l4-aton-laser-achieves-5-petawatt/)
- **SEL** (Station of Extreme Light, SHINE, Shanghai):
  - A July 2025 circular describes the "world first 100 PW-class laser facility" as under construction, with the laser–XFEL collision chamber installed (snippet). Earlier completion dates of 2023 and then 2025 slipped (snippet).
  - Sources: [SHINE Forum circular](https://sari.cas.cn/kjdt/xshd/202507/t20250703_7879010.html); [Wikipedia](https://en.wikipedia.org/wiki/Station_of_Extreme_Light)
- **LUXE** (DESY / European XFEL):
  - Beam: 16.5 GeV, 1.5×10⁹ e⁻ per bunch, 17° crossing angle, 1 Hz laser. Laser: 800 nm (1.55 eV), 30 fs, w0 ≥ 3 µm.
  - Phase-0: 40 TW, peak intensity < 1.33×10²⁰ W/cm², ξ < 7.9, χ_e < 1.5.
  - Phase-1: 350 TW, < 12×10²⁰ W/cm², ξ < 23.6, χ_e < 4.45.
  - Sources: [Yap for LUXE, arXiv:2302.04699](https://arxiv.org/abs/2302.04699) (Table 1); TDR: [arXiv:2308.00515](https://arxiv.org/abs/2308.00515), Eur. Phys. J. ST 233(10), 1709–1974; CDR: [arXiv:2102.02032](https://arxiv.org/abs/2102.02032), Eur. Phys. J. ST 230(11), 2445–2560 (2021)
  - Funding status (ESPPU 2025): the beamline to the LUXE area "have been funded by the European Union, but significant investment is still required" for the laser and detectors. — [ESPPU input](https://arxiv.org/abs/2504.02608)
- **E-320** (SLAC FACET-II): about 10 TW laser, 10 GeV beam (ESPPU gives the SLAC linac maximum as 13 GeV), "ξ ≈ 1 − 5". ESPPU parameter points: E320 (ξ, η) = (2, 0.15); proposed E320II (16, 0.15). — [ESPPU input](https://arxiv.org/abs/2504.02608)
- **Apollon, EP-OPAL, ELI-NP and ZEUS:** they appear only as ξ values in the ESPPU table. No facility-published χ values were found.

**Computed facility checks**
- LUXE ξ from intensity (linear convention, 0.8 µm): 7.888 (published 7.9, ratio 0.998) and 23.69 (published 23.6, ratio 1.004).
- LUXE χ_e with 17.2° crossing: 1.513 (published 1.5, ratio 1.009) and 4.519 (published 4.45, ratio 1.016).
- ESPPU Table 1 ξ, reproduced with peak I = 0.6·E·0.9394/τ · 2/(πw0²), w0 = FWHM/√(2 ln 2), and the "E²/4π" convention at 0.8 µm: ratios of 0.90–0.98 for FACET-II, GEMINI, ZEUS, APOLLON, CoReLS, ELI-NP and OPAL. ELI-BL and VULCAN reproduce to 0.96–0.98 only with λ ≈ 1.05 µm.
- **SEL:** the tabulated ξ = 1,800 reproduces (ratio 0.95) only from the tabulated 25 kJ in 25 fs, which is about 1 EW. SEL is publicly a 100 PW (50–100 PW) design. For 2.5 kJ in 25 fs, ξ ≈ 540.
- ESPPU's "χ ≈ 10" (5 GeV, 10²³ W/cm², head-on, 800 nm): 12.8 with the linear a₀ and 9.1 with the E²/4π a₀.
- χ per unit a₀ for a 10 GeV electron head-on at 800 nm: 0.1187. ZEUS-class ξ ≈ 100 therefore gives χ ≈ 12, and ELI-NP ξ ≈ 180 gives χ ≈ 21. These are idealised peak values and ignore radiation losses before the electron reaches the peak.

### Inferences
- There are at least three a₀ conventions in common use, and they differ by up to √2:
  1. LUXE: linear, I = ε0cE₀²/2.
  2. Di Piazza and the ESPPU table: I = E₀²/4π = ε0cE₀².
  3. "a₀ rms" (E144 used E_rms; see Q3).
- A Rust calculator should take the field amplitude E₀ as primary and carry an explicit polarisation and convention enum. Most of the 1.3×–1.4× disagreements in Q2 and Q3 are convention mismatches, not physics.
- The 2.7% gap on Di Piazza's χ0 coefficient (0.0574 against 0.059) is unexplained. It may be rounding in E_cr or ω. Treat 0.059 as approximate.
- The SEL row in the ESPPU table looks internally consistent with an energy that is not SEL's public 100 PW design. Quote SEL a₀ as "≈ 500–800 for 100 PW, depending on convention" (computed), not 1,800.

### Gaps
- No facility-published peak χ was found for ELI-NP, ELI-BL, ZEUS, SEL, EP-OPAL or Apollon. They depend on the electron source, so only the computed "χ per unit a₀" scaling is given.
- Apollon's and ELI-NP's current operating power in 2026 was not checked against facility pages. The EP-OPAL design status (NSF design study) was not verified beyond the ESPPU table.
- The ZEUS news page returned HTTP 403, so its numbers come from search snippets. One snippet's "25 quintillionths of a second" pulse length is implausible (25 fs is 25 quadrillionths) and was not used.

---

## Q3. Which experiments have reached the strong-field regime: SLAC E144, E-320, CERN NA63 crystals, and STAR 2021?

### Takeaway
- **E144 (1997):** the first laboratory nonlinear Breit–Wheeler pairs, 106 ± 14 positrons at χ ≈ 0.3 with about 5 laser photons absorbed **[mainstream]**.
- **NA63:** aligned crystals at CERN reached χ up to 1.4 (radiation emission) and 2.4 (trident), in agreement with LCFA **[mainstream]**.
- **E-320:** preliminary conference reports of nonlinear Compton scattering and quantum radiation reaction at χ ≳ 0.5 (2025). No peer-reviewed paper was found **[contested/preliminary]**.
- **STAR (2021):** 6085 exclusive e⁺e⁻ pairs from photon collisions with a 16.8 ± 2.5% cos 4Δφ modulation. The "linear Breit–Wheeler with real photons" and "vacuum birefringence" interpretations are **[contested]**.

### Cited Findings
**SLAC E144 [mainstream]**
- "A signal of 106±14 positrons above background" from 46.6 GeV electrons colliding with terawatt pulses at 527 nm (Nd:glass, frequency doubled) at the SLAC Final Focus Test Beam. Published in PRL 79, 1626 (1997). A control run at 49.1 GeV gave 22 ± 10 positrons (snippet). — [OSTI record](https://www.osti.gov/biblio/491608) (snippet)
- **Two-step interpretation (Melissinos):**
  - First, a photon backscatters "to become a high energy γ-ray (ωγ ∼ 29 GeV)", which then "scatters from at least four laser photons to produce the pair".
  - The positron yield "varies as η^2n with n = 5.1 ± 0.2".
  - The lab field is "E ∼ 3 × 10¹⁰ V/cm", which "reaches near-critical value in the rest frame" (E* = 2γE ≃ 1.8×10⁵ E).
  - Fit to exp(−πg/Υ): πg(η) = 2.01 ± 0.12 ± 0.4. After converting to peak fields and the γ-photon frame: "πg(η) = 1.77 ± 0.35 observed", against "2.12 predicted".
  - An estimate from Schwinger-type formulas predicts about 4 pairs per laser shot, against an observed 0.1 per shot.
  - E144 defined η with E_rms.
  - Source: [Melissinos, hep-ph/9805507](https://arxiv.org/abs/hep-ph/9805507)
- Parameter points for E144: ESPPU plots (ξ, η) = (0.36, 0.83). The LUXE proceedings say E144 "reached χ = 0.25 and ξ = 0.4". — [ESPPU input](https://arxiv.org/abs/2504.02608); [arXiv:2302.04699](https://arxiv.org/abs/2302.04699)
- Computed:
  - ħω(527 nm) = 2.3526 eV.
  - η = (1 + cos 17°)γħω/mc² = 0.821 (ESPPU 0.83, ratio 0.990).
  - χ = 0.36 × 0.821 = 0.296 (ESPPU 0.30, ratio 0.986).
  - LUXE's (ξ = 0.4, χ = 0.25) pair is not self-consistent with η = 0.82, which gives 0.33. Probably an rms/peak convention mismatch.
  - Threshold for a 29 GeV photon at 17°: n ≥ 2m²/(ω_γ ħω(1+cos θ)) = 3.91, so at least 4 laser photons, matching the published "at least four".
  - Significance of 106/14: 7.6σ.

**E-320 at FACET-II [contested/preliminary]**
- ExHILP 2025 abstract ("first quantitative measurements of quantum radiation reaction"):
  - The 10 GeV compressed FACET-II beam collides with "10 TW-scale NIR laser pulses".
  - May 2025 runs reached a peak a₀² ≳ 10, giving χ ≳ 0.5 at a collision angle of about 30°.
  - The electron spectrum is exponential below the linear Compton edge, with energy losses down to at least 5 GeV. The analysis is marked preliminary.
  - Source: [ExHILP 2025 contribution 935](https://indico.eli-laser.eu/event/198/contributions/935/)
- APS-DPP 2024 talk: a strongly nonlinear Compton signal with an intensity-dependent redshift of the Compton edges, and an average χ ≥ 0.12 (snippet). A 2025 abstract reports higher-order nonlinear Compton scattering and progress towards observing nonlinear Breit–Wheeler positrons (snippet). — [APS DPP 2024 TM10.7](https://meetings-archive.aps.org/dpp/2024/tm10/7/) (snippet); search summary
- A preprint estimates about 0.01–0.1 nonlinear Breit–Wheeler pairs per collision, with large backgrounds (snippet). — [arXiv:2506.04992](https://arxiv.org/abs/2506.04992) (snippet only)
- The LPHYS 2026 abstract repeats the rest-frame enhancement "~4γ²" and the recoil thresholds of χ ≳ 0.1 and ≳ 0.3. It uses electro-optic sampling with real-time Bayesian analysis for shot diagnostics. — [LPHYS 2026 talk 3978](https://www.lasphys.com/workshops/book-of-abstracts/2026/talks/3978.html)
- Computed: a₀ = √10 at 30° and 10 GeV gives χ = 0.350, which is 0.70 of the published "≳0.5". a₀ = √20 gives χ = 0.495. So the published pair (a₀² ≈ 10, χ ≈ 0.5) is consistent only if a₀² is a cycle-averaged (a₀,peak²/2) quantity.

**CERN NA63 aligned crystals [mainstream]**
- Wistisen, Di Piazza, Knudsen & Uggerhøj, Nat. Commun. 9, 795 (2018):
  - Emission spectra from ultrarelativistic positrons in silicon. "Only a full quantum theory of radiation reaction is capable of explaining the experimental results." — [arXiv:1704.01080](https://arxiv.org/abs/1704.01080)
  - Parameters (snippet): 178.2 GeV positrons; Si crystals 3.8 mm and 10 mm thick along ⟨111⟩; χ ≤ 1.4; 0.7 ≲ ξ ≲ 7.
- Nielsen et al., PRL 130, 071601 (2023): trident e⁻ → e⁻e⁺e⁻ "probes values for the strong field parameter χ up to 2.4". Data and LCFA "show remarkable agreement" over "almost 3 orders of magnitude in yield". — [arXiv:2211.02390](https://arxiv.org/abs/2211.02390)
  - Setup (snippet): 200 GeV electrons on a 400 µm Ge crystal along ⟨110⟩. Differential follow-up in PRD 108, 052013 (2023), [arXiv:2307.03464](https://arxiv.org/abs/2307.03464).
- ESPPU places NA63 at (ξ, η) = (20, 0.1), so χ = 2 (computed). It notes that questions remain "about the contribution when the photon propagator is off-shell". — [ESPPU input](https://arxiv.org/abs/2504.02608)
- Computed effective lab-frame crystal fields from E = χE_S/γ: 5.3×10¹² V/m (Si, 178.2 GeV, χ = 1.4) and 8.1×10¹² V/m (Ge, 200 GeV, χ = 2.4).

**All-optical experiments (Gemini, CoReLS) [mainstream for the observation, contested for the quantitative quantum-radiation-reaction model]**
- "Laser-based experiments have provided strong initial evidence" of quantum radiation reaction, but "precision measurements and validation of theoretical models and approximations are still outstanding". ESPPU parameter points: Gemini1 (ξ, η) = (1…2, 0.006); Gemini2 (24.7, 0.003…0.02). — [ESPPU input](https://arxiv.org/abs/2504.02608)

**STAR 2021 [observation mainstream; interpretation contested]**
- Adam et al. (STAR), PRL 127, 052302 (2021):
  - "observation of 6085 exclusive electron-positron pairs in ultra-peripheral Au+Au collisions" at √s_NN = 200 GeV.
  - "a large fourth-order angular modulation of cos 4Δφ = (16.8 ± 2.5)%".
  - "The experimental results have implications for vacuum birefringence."
  - Source: [arXiv:1910.12400](https://arxiv.org/abs/1910.12400)
- Brandenburg, Seger, Xu & Zha, "Observation of the Breit-Wheeler Process and Vacuum Birefringence in Heavy-Ion Collisions", Rep. Prog. Phys. (DOI 10.1088/1361-6633/acdae4), review the decades-long effort and "the signatures of vacuum birefringence". — [arXiv:2208.14943](https://arxiv.org/abs/2208.14943)
- Computed: the significance of the modulation from zero is 16.8/2.5 = 6.7σ.

### Inferences
- **χ reached so far:**
  - Laser–beam experiments: E144 about 0.3, E-320 about 0.5 (preliminary), Gemini below 0.3.
  - Crystals: NA63 up to 2.4.
  - Only crystals have χ > 1 with published, peer-reviewed agreement with theory. LUXE phase-1 (4.45) would be the first laser experiment clearly above χ = 1.
- **E144:** it is fairly called the first nonlinear Breit–Wheeler and trident observation. Calling it "Schwinger pair production" (Melissinos's framing) is a minority interpretation; it was perturbative multiphoton (ξ < 1), not tunnelling. Label that framing [contested].
- **STAR:** the photons are quasi-real Weizsäcker–Williams photons, and "vacuum birefringence" is inferred from the polarisation dependence of the Breit–Wheeler cross-section, not measured as a refractive-index difference. Present it as [contested] interpretation on top of a [mainstream] measurement.

### Gaps
- A search summary attributed "~1.3×10¹⁸ W/cm², Υ ≈ 0.3" for E144 to an unnamed review. It was not traced to a specific source, so it is left out of the findings.
- The E144 PRL full text (Burke et al., 1997) and Bamber et al. PRD 60, 092004 (1999) were not read. The SLAC-PUB copy failed on its TLS certificate and OSTI refused connections. Peak laser intensity, η_max, and the nonlinear-Compton harmonic orders come only from secondary sources.
- No peer-reviewed E-320 result paper was found as of October 2026.
- No published critique of STAR's "real photon" or "vacuum birefringence" claims was read. The contested label rests on the review's own framing and on general knowledge.

---

## Q4. Vacuum birefringence: Euler–Heisenberg Δn, the PVLAS final limit, and magnetar X-ray polarisation as evidence

### Takeaway
- **Euler–Heisenberg:** Δn = n∥ − n⊥ = 3A_eB², with A_e = 2α²ħ³/(45μ0m_e⁴c⁵) = 1.3247×10⁻²⁴ T⁻². That gives Δn = 3.97×10⁻²⁴ at 1 T and 2.48×10⁻²³ at 2.5 T (computed) **[mainstream theory]**.
- **PVLAS (2020, final):** Δn = (12 ± 17)×10⁻²³ at 2.5 T, an uncertainty about 7× the QED value. There is no laboratory detection.
- **Magnetars:** IXPE's 4U 0142+61 result (polarisation degree 14% → 41% with a 90° swing) is described by its authors as "further support" for vacuum birefringence. Optical RX J1856.5−3754 (16.43 ± 5.26%) was the earlier claim. Both are model-dependent **[contested]**.

### Cited Findings
- **Low-field Euler–Heisenberg Lagrangian [mainstream]:** L_EH = −F/4π + (α/360π²)(4F² + 7G²)/F_cr² (Gaussian natural units). The vacuum refractive-index change from a plane wave of intensity I0 is "of the order of (α/45π)(I0/I_cr)". — [Di Piazza et al.](https://arxiv.org/abs/1111.3886)
- **PVLAS final [mainstream; null result]:** Ejlli et al., "The PVLAS experiment: a 25 year effort to measure vacuum magnetic birefringence", Physics Reports 871, 1–74 (2020) (journal per snippet).
  - At B = 2.5 T: "Δn(PVLAS) = (12 ± 17) × 10⁻²³" and "|Δκ|(PVLAS) = (10 ± 28) × 10⁻²³".
  - QED predicts Δn(QED) = 2.5×10⁻²³ at 2.5 T. The uncertainty on Δn is about a factor of 7 above that prediction.
  - Source: [arXiv:2005.12913](https://arxiv.org/abs/2005.12913)
  - The brief's "Eur. Phys. J. C 2020" citation does not match; the final result is the 2020 Physics Reports paper. The earlier EPJ C paper is Della Valle et al. 2016, [arXiv:1510.08052](https://arxiv.org/abs/1510.08052).
- **Earlier PVLAS (2014, new apparatus):** unitary-field Δn_u = (4 ± 20)×10⁻²³ T⁻², with QED Δn_u = 4×10⁻²⁴ T⁻² (snippet). — [arXiv:1406.6518](https://arxiv.org/abs/1406.6518)
- **IXPE 4U 0142+61** (Taverna et al., Science 378, 646, 2022, DOI 10.1126/science.add0080):
  - Linear polarisation "(12±1)%" over 2–8 keV; "(14±1)% at 2–4 keV and (41±7)% at 5.5–8 keV".
  - The degree drops below sensitivity at 4–5 keV, where "the polarization angle swings by ~90°".
  - The authors say the observations "lend further support to the presence of the quantum mechanical effect of vacuum birefringence".
  - Source: [arXiv:2205.08898](https://arxiv.org/abs/2205.08898)
- **IXPE 1RXS J170849.0−400910** (Zane et al. 2023), from snippets:
  - About 35% polarisation integrated over 2–8 keV, rising from about 20% (2–3 keV) to about 80% (6–8 keV) at constant angle.
  - The small emitting regions "do not allow a robust conclusion" about vacuum birefringence.
  - Source: [arXiv:2301.12919](https://arxiv.org/abs/2301.12919)
- Lai interprets the 4U 0142+61 swing as photon mode conversion at the QED "vacuum resonance" in a magnetar atmosphere, needing B ≳ 5×10¹⁴ G for 1RXS J1708 (snippet). — [arXiv:2209.13640](https://arxiv.org/abs/2209.13640)
- **RX J1856.5−3754** (Mignani et al., MNRAS 2017, DOI 10.1093/mnras/stw2798): VLT optical polarimetry gives "P.D. = 16.43% ± 5.26%" and position angle 145.39° ± 9.44°. This "is large enough to support the presence of vacuum birefringence". — [arXiv:1610.08323](https://arxiv.org/abs/1610.08323)

**Computed test vectors**

| Quantity | Computed | Published | Ratio |
|---|---|---|---|
| A_e = 2α²ħ³/(45μ0m_e⁴c⁵) | 1.32467×10⁻²⁴ T⁻² | 1.32×10⁻²⁴ (commonly quoted; not in a fetched text) | 1.0035 |
| Δn_u = 3A_e | 3.974×10⁻²⁴ T⁻² | 4×10⁻²⁴ T⁻² (PVLAS 2014) | 0.994 |
| Δn at 1 T | 3.974×10⁻²⁴ | — | — |
| Δn at 2.5 T | 2.484×10⁻²³ | 2.5×10⁻²³ (PVLAS 2020) | 0.994 |
| n∥ − 1 = (7/2)A_eB² at 2.5 T | 2.898×10⁻²³ | — | — |
| n⊥ − 1 = 2A_eB² at 2.5 T | 1.656×10⁻²³ | — | — |
| PVLAS σ / Δn_QED | 6.84 | "about a factor 7" | 0.98 |
| PVLAS (meas − QED)/σ | 0.56σ | — | — |
| PVLAS per T²: 12/6.25, 17/6.25 | (1.92 ± 2.72)×10⁻²³ T⁻² | — | — |
| SI EH coefficient 2α²ε0²ħ³/(45m⁴c⁵) | 1.305×10⁻⁵² m³/J | — | — |
| B/B_S at 10¹⁴ G | 2.27 | — | — |

SI form used: L_EH = (2α²ε0²ħ³/(45m⁴c⁵))[(E² − c²B²)² + 7c²(E·B)²]. The split n∥ = 1 + (7/2)A_eB² and n⊥ = 1 + 2A_eB² is the standard PVLAS convention. It is consistent with the published 2.5×10⁻²³, but was not read in a fetched text.

### Inferences
- Laboratory vacuum birefringence is still undetected (2026). The best limit sits about 7× above QED. Astrophysical "evidence" requires atmosphere or condensed-surface emission models plus field geometry. Label it [contested], not [mainstream] detection.
- Magnetar dipole fields of order 10¹⁴–10¹⁵ G exceed B_S. The weak-field Δn = 3A_eB² formula does not apply there. A calculator must guard B ≪ B_S, or use the strong-field Euler–Heisenberg expressions, which are not covered here.
- PVLAS's |Δκ| (dichroism) limit also constrains axion-like particles and millicharged particles. That is out of scope here.

### Gaps
- The PVLAS Physics Reports full text and its exact unitary-field final number were not read; the per-T² values above are computed.
- No independent critical assessment of the 4U 0142+61 vacuum-birefringence inference was read beyond Lai (snippet) and Zane et al. (snippet).

---

## Q5. Supercritical fields in heavy-ion collisions (Z₁ + Z₂ > 173): the GSI positron-line episode and current status

### Takeaway
For extended nuclei the 1s level reaches the negative-energy continuum at Z ≈ 173; the point-nucleus Dirac limit is Z = 1/α = 137 (computed) **[mainstream]**.
- **Claims:** GSI's EPOS and ORANGE groups reported narrow positron lines (about 250–400 keV) and later e⁺e⁻ sum-energy lines in U+Cm, U+U and U+Th collisions (1983 onwards).
- **Outcome:** improved experiments (APEX at Argonne, EPOS II, ORANGE) failed to reproduce them. A new-particle explanation was excluded by e⁺e⁻ (Bhabha) scattering. The consensus attributes the earlier signals to nuclear pair conversion. Spontaneous vacuum decay has never been observed **[refuted as a vacuum-decay or new-particle signal]**.
- **Today:** theory (Maltsev et al. 2019) proposes impact-parameter-resolved pair measurements at CRYRING@GSI/FAIR. No such experiment has been reported.

### Cited Findings
- **Critical charge:** the point at which the 1s wavefunction joins the negative-frequency continuum "has been determined to be about Z = 173". The vacuum becomes charged "(charge 2 for Z > 172) in supercritical fields". — [Rafelski, Kirsch, Müller, Reinhardt & Greiner, arXiv:1604.08690](https://arxiv.org/abs/1604.08690)
- The Shabaev group: Z ≈ 173 for extended nuclei, against about 137 for a point nucleus (snippet). — [arXiv:2008.05005](https://arxiv.org/abs/2008.05005)
- **Timescales:**
  - Supercritical K-vacancy lifetime against positron emission: τ_e+ ≃ 10⁻¹⁸–10⁻¹⁹ s.
  - Non-sticking collision time just below the Coulomb barrier: τ_coll ≃ 2R_cr/v ≃ 0.25×10⁻²⁰ s, with R_cr ≃ 35 fm.
  - "The emission time for positrons is typically 100 times longer", so the expected yield is about 1%. The Sommerfeld parameter is Z₁Z₂α/v > 500.
  - Source: [Rafelski et al.](https://arxiv.org/abs/1604.08690)
- **Maltsev et al.** (PRL 123, 113401, 2019): U–U critical distance R_cr ≈ 32.6 fm; calculations at closest approach R_min = 16.5 fm, with the nuclei "about 1–2 fm away from touching". Supercritical signature: at fixed R_min, a pair probability that increases as the collision energy decreases, and a "pronounced decrease of dP/dη at η ≈ 1" as a function of Z. The authors write that "such experiments seem feasible with the CRYRING facility at GSI/FAIR, where storing of bare uranium nuclei at low energies is anticipated". — [arXiv:1903.08546](https://arxiv.org/abs/1903.08546)
- In slow collisions (about 0.1c) the supercritical field lasts only about 10⁻²¹ s, much shorter than the spontaneous resonance lifetime (snippet). — [Chinese Physics C article](https://cpc.ihep.ac.cn/article/id/74b39231-07ce-466b-b9b5-bad577ff5111)

**GSI positron lines [refuted]**
- Original claim: Schweppe et al., "Observation of a Peak Structure in Positron Spectra from U+Cm Collisions", PRL 51, 2261 (1983). — [GSI repository record](https://repository.gsi.de/record/50234)
- Narrow peaks at about 250, 330 and 400 keV were reported in systems with Z₁ + Z₂ = 180–188, on a broad background. Back-to-back correlated emission was absent, which excludes a two-body decay of a neutral 1.5–2 MeV state (snippet). — [arXiv:0908.1502](https://arxiv.org/abs/0908.1502)
- An early reading as the e⁺e⁻ decay of a neutral particle near 1.8 MeV/c² was ruled out by Bhabha-scattering experiments (snippet). — [arXiv:nucl-ex/9706009](https://arxiv.org/abs/nucl-ex/9706009)
- **Rafelski et al.'s history:**
  - Questions arose "when the experimental groups found positrons produced in subcritical systems", and electrons were then seen with the positrons.
  - "All efforts were ended when improved experiments failed to find peaked positron lines where earlier experiments had seen them."
  - "The consensus view today is that the earlier intriguing observations were due to highly system dependent nuclear excitations converting into pairs."
  - The new-particle reading "has been ruled out in the study of e⁺e⁻ resonant scattering".
  - Source: [arXiv:1604.08690](https://arxiv.org/abs/1604.08690)
- APEX (Argonne): "No evidence is found for narrow structures similar to those previously reported". Its upper limits on a neutral particle of 1.4–2.1 MeV/c² are "significantly less than previously reported" (snippet). A comment argued that APEX was less sensitive to extremely narrow lines (snippet). — [OSTI 196563](https://www.osti.gov/biblio/196563) (snippet)
- The later double-ORANGE ²³⁸U + ¹⁸¹Ta measurement shows no line structures within its statistics. It notes that internal pair conversion "can be considered a possible mechanism" for some earlier weak lines (snippet). — [arXiv:nucl-ex/0004009](https://arxiv.org/abs/nucl-ex/0004009)
- Maltsev et al.: "no sign of the spontaneous pair production or the diving phenomenon had been found" at GSI. — [arXiv:1903.08546](https://arxiv.org/abs/1903.08546)

**Computed**
- Point-nucleus collapse at Z = 137.036.
- Point-nucleus 1s binding: 101.6 keV for Z = 82 and 132.3 keV for Z = 92. Finite nuclear size reduces these; these are point-nucleus test values only.
- The velocity implied by τ = 2R/v with R = 35 fm and τ = 0.25×10⁻²⁰ s is 0.093c, consistent with "~0.1c".
- e²/(4πε0) = 1.43996 MeV·fm. The U+U Coulomb energy at R = 16.5 fm is 738.7 MeV. Head-on U+U reaching 16.5 fm needs E_lab ≈ 6.21 MeV/u.

### Inferences
- The "supercritical vacuum decay" prediction itself is still [mainstream] theory. What is [refuted] is the 1980s claim that the GSI lines were its signal, or a new 1.8 MeV boson.
- A clean experiment needs either nuclear sticking times ≳ 10⁻²⁰ s or the differential dP/db signature of Maltsev et al.
- No FAIR or CRYRING supercritical-collision experiment was found to be scheduled or reported by October 2026. "FAIR plans" means theory proposals and the bare-U storage capability, not an approved run.

### Gaps
- No FAIR, SPARC or CRYRING experiment proposal or schedule for supercritical U+U collisions was found. Shabaev-group talks (2021–2023) are slides only.
- The APEX PRL (Ahmad et al. 1995/1997) and the EPOS II and ORANGE refutation papers were not read in full text; OSTI refused connections.

---

## Q6. The dynamically assisted Schwinger mechanism and condensed-matter analogues (graphene, Berdyugin et al. 2022)

### Takeaway
Adding a weak, fast field to a strong, slow one lowers the tunnelling exponent once the combined Keldysh parameter γ = mω/(qE) exceeds π/2. At γ = π the exponent falls by about 40%, which turns an exp(−πE_S/E) suppression of 2×10⁻¹⁴ at 0.1 E_S into about 5×10⁻⁹ (computed) **[mainstream theory, not yet observed in QED]**.
- **Graphene:** Berdyugin et al. (Science 2022) report a Schwinger-like interband e–h pair-production regime in graphene and its superlattices. It has an E^(3/2) production rate, j_c = n e v_F critical currents, and dV/dI ∝ j^(−1/3) or j^(+1/3) scaling **[mainstream as an analogue]**.

### Cited Findings
**Dynamically assisted Schwinger mechanism** (Schützhold, Gies & Dunne, PRL 101, 130404, 2008)
- Field model: E(t) = E/cosh²(Ωt) + ε/cosh²(ωt), with E_S ≫ E ≫ ε and m ≫ ω ≫ Ω.
- Combined Keldysh parameter: γ = mω/(qE), built from the fast pulse's frequency and the slow field's strength.
- For γ < π/2 the result reduces to Schwinger's, A_inst = πm²/(qE).
- For γ ≥ π/2: A_inst ≈ (m²/qE)[2 arcsin(π/(2γ)) + (π/(2γ²))√(4γ² − π²)].
- At γ = π the action "is reduced by about 40%". For γ → ∞ it decays as 1/γ.
- Rate: ℑΓ ∼ exp{−A_inst}. "a reduction of 40% in the exponent could make the difference between a suppression of 10⁻¹⁰ and 10⁻⁶".
- Source: [arXiv:0807.0754](https://arxiv.org/abs/0807.0754)
- Computed checks of eq. (12):
  - A(π/2)/π = 1.0000.
  - 1 − A(π)/π = 0.391 (published "about 40%").
  - A(100)·100/(2π) = 0.99996, matching the 1/γ tail (asymptote 2π/γ).
  - At E = 0.1 E_S, γ = π needs an assisting photon of ħω = π·0.1·mc² = 160.5 keV.
  - The exponent then goes from exp(−10π) = 2.27×10⁻¹⁴ to exp(−10·A(π)) = 4.91×10⁻⁹.
- Computed: the ordinary Keldysh parameter for an 800 nm field at 0.1 E_S is γ = (ħω/mc²)(E_S/E) = 3.0×10⁻⁵. An optical laser alone is deep in the tunnelling regime.

**Graphene Schwinger analogue**
- Allor, Cohen & McGady (PRD 78, 096009, 2008): "The rate of pair production per unit area for four species of massless fermions" is π⁻² ħ^(−3/2) c̃^(−1/2) (qE)^(3/2), where c̃ is the 2D "speed of light", the Fermi velocity v_F. — [arXiv:0708.1471](https://arxiv.org/abs/0708.1471)
- **Berdyugin et al.**, "Out-of-equilibrium criticalities in graphene superlattices", Science 375, 430 (2022), DOI 10.1126/science.abi8627 **[mainstream analogue]**:
  - Critical currents set in at j_c ≈ n e v_F, when the drift velocity reaches v_F.
  - Above j_c, interband carriers are generated by the Schwinger (Zener–Klein) mechanism "at a rate ∝ E^(3/2)". This gives extra current ∝ V^(3/2) and dV/dI ∝ j^(−1/3).
  - At the neutrality point the e–h plasma gives n ∝ E^(3/4) and dV/dI ∝ j^(+1/3), "in quantitative agreement with our experiment".
  - Example numbers:
    - Monolayer graphene: v_F = 1×10⁶ m/s, n = 0.4×10¹⁶ m⁻², constriction 0.4 µm wide, T = 2 K. Peaks persist up to about 50 K.
    - dV/dI minima: "j_m ≈ 0.05 mA µm⁻¹", which "translates into n = j_m/ev_F ≈ 3×10¹⁰ cm⁻²".
    - Superlattices: characteristic v_F of "several 10⁴ m s⁻¹".
  - Graphene mimics a QED effect "predicted to occur at enormous fields of 10¹⁸ V m⁻¹".
  - Source: [arXiv:2106.12609](https://arxiv.org/abs/2106.12609)
- Computed:
  - n = j_m/(e v_F) = 3.12×10¹⁰ cm⁻² (published ≈ 3×10¹⁰, ratio 1.04).
  - j_c = n e v_F at n = 0.4×10¹⁶ m⁻² is 0.641 mA/µm.
  - Allor–Cohen–McGady rate (4 species): 6.0×10²⁴, 1.90×10²⁶ and 6.0×10²⁷ m⁻²s⁻¹ at E = 10⁴, 10⁵ and 10⁶ V/m for v_F = 10⁶ m/s. With v_F = 3×10⁴ m/s, multiply by √(10⁶/3×10⁴) = 5.77.

### Inferences
- Dynamically assisted pair production in real QED would need roughly 100 keV photons overlapping a field of order 0.1 E_S. That is still far beyond the 10⁻³ E_S of the CoReLS record. No QED observation exists; it is a theoretical enhancement.
- In graphene the "mass gap" is zero or small, so there is no exp(−πE_S/E) threshold. The analogue tests the E^(3/2) prefactor scaling and the pair-plasma dynamics, not the non-perturbative exponential.
- Gapped analogues, such as Landau–Zener in semiconductors, trapped-ion or cold-atom simulators, and superlattice minibands, would test the exponent. Berdyugin is the high-profile 2022 example.

### Gaps
- Berdyugin's supplementary derivation (ref. 29 in the paper), with the exact prefactor and the annihilation model, was not read.
- No 2023–2026 follow-ups (graphene, bilayer, cold atoms, trapped ions) were searched.
- No experimental test of dynamically assisted QED pair production was found. The LUXE γ-laser mode and XFEL plus optical-laser proposals (e.g. SEL collision chamber, HIBEF) were not checked for explicit dynamically assisted Schwinger programmes.
