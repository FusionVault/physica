# Testing Puthoff's polarizable-vacuum (PV) / exponential metric against GR

Scope and conventions.
- PV metric (Puthoff 2002, eqs 17 and 20): ds² = −(c²/K)dt² + K(dx² + dy² + dz²), K = exp(2GM/rc²), with r the isotropic radius. Unless stated, units are G = c = M = 1.
- "PV metric" means the static exponential geometry, which is also Yilmaz's and Papapetrou's. "PV theory" means Puthoff's (Dicke's) scalar Lagrangian dynamics, which is a different object. Several tests separate the two.
- Every number marked **[computed]** comes from my scripts in this folder:
  - `pv_metric_tests.py` → `pv_metric_tests_output.txt`
  - `pv_qnm_check.py` → `pv_qnm_check_output.txt`
  - `pv_extra_values.py` → `pv_extra_values_output.txt`
  - All three run with `python3 -P <script>`.
- Downloads are in `pv_sources/<id>/` and were treated as data only.
- A parallel worker's `puthoff_metric_checks.out.txt` in this folder independently gives the same expansions, 4π versus 15π/4, a PV third-order coefficient of 47.97 (my fit: 48.0000), and δ_shadow = +0.0463.

## 1. Second post-Newtonian order: coefficients and which experiments constrain them

### Takeaway
- PV and GR agree at 1PN (β = γ = 1). They first differ in the U² term of g_ij, where the 2PN parameter is δ = 4/3 for PV against 1 for GR, and in the U³ term of g₀₀.
- Every solar-system and S2 2PN observable differs from GR by 10⁻³ to 10⁻⁷ of current precision, so no existing weak-field experiment discriminates.
- The only proposed weak-field experiment that would is LATOR. Its quoted 2×10⁻³ accuracy on the second-order deflection gives about 33σ, because PV's 2PN deflection is 1/15 = 6.7% larger.

### Cited Findings
- **Puthoff's own expansions** ([Puthoff 2002, arXiv:gr-qc/9909037](https://arxiv.org/abs/gr-qc/9909037), eqs 21 and 45–46):
  - g₀₀ = e^{−2GM/rc²} = 1 − 2(GM/rc²) + 2(GM/rc²)² − …; g_ii = −(1 + 2GM/rc² + …).
  - Comparing with the isotropic Schwarzschild eqs 42–43, the paper says "to the order of expansion shown, the two metric tensors are identical. Since the three classical tests of GR do not require terms beyond these…".
- **Both expansions, from Ibison 2006 eq 3** ([arXiv:0705.0080](https://arxiv.org/abs/0705.0080), CQG 23, 577):
  - Yilmaz/exponential: g₀₀ = 1 − 2φ + 2φ² − (4/3)φ³, g₁₁ = −1 − 2φ − 2φ².
  - GR: g₀₀ = 1 − 2φ + 2φ² − (3/2)φ³, g₁₁ = −1 − 2φ − (3/2)φ².
  - Ibison: "To date these give observationally indistinguishable predictions for red-shift, light bending and perihelion advance."
- **Boonserm et al.** ([arXiv:1805.03781](https://arxiv.org/abs/1805.03781), PRD 98, 084048, §1): in weak fields g_ab = η_ab + (2m/r)δ_ab, which "exactly matches the lowest-order weak-field expansion of general relativity". "Strong-field behaviour… and even medium-field behaviour (2m/r ~ 1), is markedly different."
- **Psaltis et al. 2020** ([arXiv:2010.01055](https://arxiv.org/abs/2010.01055), PRL 125, 141104):
  - In areal coordinates, −g_tt = 1 − 2/r + 2(β̄−γ̄)/r² − 2ζ/r³ + ….
  - For a metric that terminates at 2PN, r_ph = 3 + 5ζ/9 and r_sh = 3√3(1 + ζ/9). The shadow therefore "probes the behavior of the spacetime at least at the 2PN order".
- **S2 orbit, as quoted in EHT Sgr A* Paper VI, subsection "Tests with the S2 Orbit"** ([arXiv:2311.09484](https://arxiv.org/abs/2311.09484)):
  - GRAVITY gives f_SP = 1.1 ± 0.19, so 2γ − β = 1.3 ± 0.57.
  - The 2PN term is Δφ₂ = −6π[GM/(2a(1−e²)c²)]²(10 − e²), with Δφ₂/Δφ₁ ≃ 7×10⁻⁴.
  - "the second post-Newtonian term would have to be ∼285 times larger than the Schwarzschild prediction in order for it to cause deviations detectable with current instruments."
- **Cassini** ([Will 2014 Living Review, arXiv:1403.7377](https://arxiv.org/abs/1403.7377), Shapiro time-delay section):
  - γ − 1 = (2.1 ± 2.3)×10⁻⁵, with the spacecraft at 8.43 AU and closest approach 1.6 R_⊙.
  - The table "Current limits on the PPN parameters" gives γ − 1: 2.3×10⁻⁵ (Cassini) and 2×10⁻⁴ (VLBI deflection); β − 1: 8×10⁻⁵ (perihelion).
- **LATOR (proposed)** ([Turyshev et al., arXiv:gr-qc/0410044](https://arxiv.org/abs/gr-qc/0410044)):
  - "The second order light deflection is approximately 1700 pm and with 10 pm accuracy and the adopted measurement strategy it could be measured with accuracy of ∼2×10⁻³, including first ever measurement of the PPN parameter δ."
  - Its stated target for γ is 1 part in 10⁸.
- **[computed] Metric series.** Fitted Taylor coefficients reproduce Ibison's eq 3 exactly. In the 2PN convention g_ij = 1 + 2γU + (3/2)δU², PV has δ = 4/3 and GR has δ = 1.
  - In areal radius, with the EHT/Psaltis convention g_tt = −1 + 2/R − 2κ₁/R² + 2κ₂/R³ − 2κ₃/R⁴, PV has κ₁ = 0, κ₂ = 1/6, κ₃ = −1/3. These follow from −g_tt = (U_A/u)², with u the Lambert-W inverse of u e^{−u} = U_A.
- **[computed] Light bending to 4th order**, with exact null-geodesic quadrature and a polynomial fit in x = GM/(c²b):
  - GR: α = 4x + 11.780972x² + 42.6667x³ + 170.09x⁴, which is 4, 15π/4, 128/3 and 3465π/64.
  - PV: α = 4x + 12.566371x² + 48.0000x³ + 201.06x⁴, which is **4, 4π, 48, 64π**.
  - The PV 2PN coefficient matches the general form π(2 + 2γ − β + 3δ/4): 15π/4 for GR and 4π for PV.
  - Optical-medium view: PV's index is n = 1 + 2u + 2u² and GR's isotropic index is n = 1 + 2u + (7/4)u². Then α₂ = π(n₂ + 2).
- **[computed] Sun-grazing ray (b = R_⊙):**
  - 1PN: 1.751190″.
  - 2PN: GR 10.947 µas, PV 11.677 µas, difference **0.730 µas**.
  - Equivalent to a shift in γ of (π/8)GM/(bc²) = 8.3×10⁻⁷, or 0.004σ of VLBI's 2×10⁻⁴.
  - LATOR at 2×10⁻³ of the 2PN term would give (1/15)/0.002 ≈ **33σ**.
- **[computed] Shapiro delay at 2PN.**
  - Method: the exact travel time T = ∫G(w) w dw/√(w² − b²), with optical radius w = n r and G = n dr/dw.
  - At fixed b, a full pass gives ΔT = (π/2)(GM/c²)²/(bc) = 16.42 ps at R_⊙. That comparison is not physical, because the outgoing direction also differs.
  - With emitter and receiver fixed (ΔT|_Θ = ΔT|_b − bΔα|_b), **ΔT = (π/4)(GM/c²)²/(bc)**: 8.21 ps at R_⊙ and 5.13 ps at Cassini's 1.6 R_⊙.
  - In Cassini's Doppler observable this mimics δγ = (π/8)GM/(bc²) = 5.2×10⁻⁷, against Cassini's σ = 2.3×10⁻⁵ (0.02σ).
- **[computed] Periastron advance, gauge-invariant form.** For near-circular orbits, K = Ω_φ/Ω_r as a function of x = (GMΩ/c³)^{2/3}:
  - GR: K = 1 + 3x + 13.5x² + 67.5x³, matching (1 − 6x)^{−1/2}.
  - PV: **K = 1 + 3x + 14.5x² + 76.5x³**.
  - The extra 2PN advance is therefore 2πx² per orbit. For Mercury that is 4.1×10⁻¹⁵ rad per orbit, or **0.35 µas per century**, against the β − 1 limit of 8×10⁻⁵ (about 3 mas per century).
- **[computed] S2 at identical E and L in both metrics** (a = 125 mas at 8246.7 pc, M = 4.261×10⁶ M_⊙, e = 0.88466; p = 5328 M):
  - Δφ_GR = 12.182′ per orbit.
  - PV − GR = **2.7×10⁻⁷ rad per orbit**, which is 7.6×10⁻⁵ of the precession and 4×10⁻⁴σ of the GRAVITY measurement.
- **[computed] Static redshift.** It is identical through 2PN, since β = 1 for both.
  - At the same areal radius, z_PV − z_GR ≈ U_A³/6. That is 6×10⁻¹¹ at S2's pericentre (U_A = 1/1400).
  - At a neutron-star surface with U_A = 0.172: z_GR = 0.23466 and z_PV = 0.23712, a 1.0% difference.

### Inferences
- **2PN status: untested.** PV passes every current weak-field test only because those tests stop at 1PN. That is a statement of indistinguishability, not a confirmation.
- The ranked 2PN discriminators are:
  1. Second-order solar light deflection: 6.7% of a 10.9 µas effect. It is the only one with a credible proposed mission (LATOR, ~33σ), but LATOR was never flown and its 2×10⁻³ is a design claim.
  2. The 2PN Shapiro term: ~5–8 ps one-way. It needs γ-equivalent precision of a few ×10⁻⁸.
  3. Periastron advance: negligible everywhere in the solar system and for S2.
  - The redshift difference appears only at 3PN.
- **Neutron-star test.** In principle, NICER-type pulse-profile fits of mass and radius could see the 1% surface-redshift difference and the PV light bending near the surface. But PV has no published neutron-star interior, so the mass–radius relation is unconstrained and the test is not yet well posed.

### Gaps
- I did not open Epstein & Shapiro (1980) or Richter & Matzner for the general 2PN deflection formula. The formula π(2 + 2γ − β + 3δ/4) is checked numerically here (both 4π and 15π/4 reproduced) but not quoted from them.
- LATOR does not define its PPN δ in the text I read, so identifying it with the coefficient of (3/2)U² in g_ij is my assumption.
- I did not open GRAVITY 2020 for S2's eccentricity; I used 0.88466 from memory. The PV − GR result is insensitive to it at the stated precision. I found no ASTROD or Gaia document giving a 2PN forecast.

## 2. Black-hole shadows: PV against the EHT

### Takeaway
- PV's static shadow is 2e GM/c², which is **4.627% larger** than Schwarzschild's 3√3 GM/c² (δ_PV = +0.0463).
- Against the EHT Sgr A* fiducial results, PV sits at 1.40σ (VLTI prior) and 0.96σ (Keck prior). Schwarzschild sits at 0.89σ and 0.44σ.
- Against the Keck+VLTI combination of Vagnozzi et al., PV is at 1.6σ: outside 1σ, inside 2σ. Against M87*, PV is at 0.33σ.
- **Current data mildly prefer GR, with Gaussian likelihood ratios of 1.4–2.5, but do not exclude PV.**
- No paper has published a dedicated EHT fit of the exponential metric. Makukov & Mychelkin (2018) and Mychelkin et al. (2024) say only that the shadow is "5%" or "over 5%" larger.
- The wormhole that the EHT and Vagnozzi exclude is **not** the PV metric.

### Cited Findings
- **Boonserm et al. §7 (ISCO and photon sphere)** ([arXiv:1805.03781](https://arxiv.org/abs/1805.03781)): the photon potential V₀ = e^{−4m/r}L²/r² "has a single peak at r = 2m", with V₀,max = L²/(2me)². This "unstable photon sphere at r = 2m" corresponds to areal r_s = 2m e^{1/2} ≈ 3.297m.
- **Makukov & Mychelkin 2018** ([arXiv:1809.05290](https://arxiv.org/abs/1809.05290), PRD 98, 064050, subsection "Shadows of compact objects"):
  - "for the antiscalar case… R_sh = 2eM = 5.437M, which… proves to be 5% larger than in the vacuum case."
  - Abstract: "the antiscalar solution predicts 5% larger shadows". They "restrict… to the static case, leaving the effect of rotation on the shadow image for later studies".
- **Mychelkin, Suliyeva & Makukov 2025** ([arXiv:2510.15391](https://arxiv.org/abs/2510.15391), IJMPD 35, 2550085, introduction) on their 2024 GRG paper: "an evaluation… of the strong-field shadow effect yielded a notable increase in shadow size (over 5%) compared to the vacuum solution", which "might serve as a potentially observable signature, provided the central mass estimate is refined." I did not open the 2024 GRG paper (doi 10.1007/s10714-024-03227-7).
- **Turimov, Turaev, Ahmedov & Stuchlík 2022** (Phys. Dark Univ. 35, 100946, [doi:10.1016/j.dark.2021.100946](https://doi.org/10.1016/j.dark.2021.100946)). Known only from search-result abstract text; not opened: "the size of the shadow of the wormhole is greater than the shadow of the Schwarzschild black hole."
- **EHT Sgr A* Paper VI** ([arXiv:2311.09484](https://arxiv.org/abs/2311.09484)):
  - Table 2 gives the Schwarzschild deviation δ (columns GRMHD / analytic Kerr / analytic non-Kerr):

    | Pipeline | Prior | GRMHD | Analytic Kerr | Analytic non-Kerr |
    |---|---|---|---|---|
    | eht-imaging | VLT(I) | −0.08 +0.09/−0.09 | −0.05 +0.09/−0.11 | −0.07 +0.10/−0.09 |
    | eht-imaging | Keck | −0.04 +0.09/−0.10 | −0.02 +0.11/−0.11 | −0.02 +0.10/−0.09 |
    | SMILI | VLT(I) | −0.10 +0.12/−0.10 | | |
    | SMILI | Keck | −0.06 +0.13/−0.10 | | |
    | DIFMAP | VLT(I) | −0.12 +0.10/−0.08 | | |
    | DIFMAP | Keck | −0.08 +0.09/−0.09 | | |
    | mG-ring | VLT(I) | −0.17 +0.11/−0.10 | | |
    | mG-ring | Keck | −0.13 +0.11/−0.11 | | |

  - Fiducial bounds (section "Metric tests from shadow size"): Keck 3√3(1 − 0.14) ≲ r̃_sh ≲ 3√3(1 + 0.05), i.e. 4.5–5.5 M; VLTI 3√3(1 − 0.17) ≲ r̃_sh ≲ 3√3(1 + 0.01), i.e. 4.3–5.3 M.
  - The measured shadow diameters span "46.9–50.0 µas" for the most likely values, with a 68% envelope of "41.7–55.6 µas" (subsection "The Diameter of the Black-Hole Shadow").
  - The Kerr range is "∼8%" (δ-posterior figure caption).
  - Priors (subsection "Priors on θ_g"): Keck M = 3.951×10⁶ M_⊙ at R₀ = 7935 pc; VLTI M = 4.297×10⁶ M_⊙ at 8277 pc.
  - The paper bounds "of order ∼1 and ∼5" the 1PN and 2PN coefficients κ₁ and κ₂, e.g. JP: −3.1 ≲ κ₂ ≲ 1.5.
  - "we do not use these post-Newtonian expansions in order to calculate black-hole shadows, which would have been inappropriate".
  - Its wormhole example is a Morris–Thorne metric with g_tt = −exp(2Φ), Φ = −r_t/r_A in **areal** radius and r_t = M. It is marked "×" (entirely ruled out) in Table 4, with listed κ₁ = −1, κ₂ = −2/3.
  - The very presence of a central brightness depression rules out compact objects without photon spheres (subsection "Surfaceless Horizonless Compact Objects"). PV has one, at r = 2M.
- **Vagnozzi et al.** ([arXiv:2205.07787](https://arxiv.org/abs/2205.07787); I did not check the journal reference), "Horizon-scale tests of gravity theories and fundamental physics from the EHT image of Sagittarius A*":
  - Averaging Keck and VLTI gives −0.125 ≲ δ ≲ 0.005 (1σ) and −0.19 ≲ δ ≲ 0.07 (2σ), i.e. 4.55 ≲ r_sh/M ≲ 5.22 (1σ) and 4.21 ≲ r_sh/M ≲ 5.56 (2σ).
  - Their Morris–Thorne wormhole with A(r) = e^{−2M/r} in areal r has "r_sh/M = e ≈ 2.718", excluded at >5σ (Fig. caption).
- **EHT M87* Paper VI** ([arXiv:1906.11243](https://arxiv.org/abs/1906.11243), shadow-size test section): "δ = −0.01 ± 0.17 (68% credible intervals) for the stellar and δ = 0.78 ± 0.3 for the gas dynamics priors."
- **Psaltis 2020** ([arXiv:2010.01055](https://arxiv.org/abs/2010.01055)): the M87 shadow is consistent with Kerr to "within ∼17%"; gives the r_sh = 3√3(1 + ζ/9) formula.
- **BHEX photon-ring paper** ([arXiv:2406.09498](https://arxiv.org/abs/2406.09498), "Lensing bands" subsection): an n = 1 ring diameter "would provide an estimate of the black hole mass-to-distance ratio accurate to within ∼10% based on purely geometric considerations derived from strong lensing alone", improving with mild source assumptions.
- **[computed] PV shadow.**
  - b_c = 2e = 5.436564 against √27 = 5.196152, so δ_PV = 0.046267.
  - Sgr A*: angular diameter 55.72 µas (VLTI prior) against GR 53.25 µas, and 53.44 µas (Keck) against 51.08 µas.
  - Psaltis's truncated-2PN formula with ζ = κ₂ = 1/6 gives only +1.85%. Most of PV's +4.63% comes from orders beyond 2PN.
- **[computed] Tension**, as z = (δ_PV − δ̂)/σ₊:

  | Data | PV | Schwarzschild | Likelihood ratio GR:PV |
  |---|---|---|---|
  | eht-imaging VLTI GRMHD (fiducial) | 1.40σ | 0.89σ | 1.80 |
  | eht-imaging Keck GRMHD (fiducial) | 0.96σ | 0.44σ | 1.43 |
  | Most PV-favourable row (eht-imaging Keck analytic Kerr) | 0.60σ | 0.18σ | 1.18 |
  | Least PV-favourable row (mG-ring VLTI) | 1.97σ | 1.55σ | 2.09 |
  | Vagnozzi et al. Keck+VLTI combination | 1.63σ | 0.92σ | 2.49 |
  | M87* stellar-dynamics prior | 0.33σ | 0.06σ | 1.05 |

  - The PV–Schwarzschild separation is 0.51 fiducial σ. Discriminating them needs σ_δ ≈ 0.015 for 3σ, or 0.009 for 5σ.

### Inferences
- **Shadow status: mildly disfavoured but alive** (about 1–2σ for Sgr A*, 0.3σ for M87*).
- Two caveats tilt this further against PV:
  - All the δ posteriors centre on negative values.
  - Spin can only shrink a Kerr shadow (by up to about 7–8%), so GR spans δ ∈ [−0.08, 0], while static PV sits at +0.046.
- A spinning PV object might also shrink. Makukov & Mychelkin (2018) built an exact rotating Papapetrou generalisation but have not computed its shadow, so how PV's prediction depends on spin is unknown.
- The literature's "5% larger" is a rounding of 4.63%, and "over 5%" is not supported by 2e/√27.
- **Resolving the sign discrepancy.** For Φ = −M/r_A, −g_tt = e^{−2/r} gives κ₁ = +1 and κ₂ = +2/3 in the EHT's own convention (that convention returns q̄²/2 for RN and −ᾱ²/4 for KS, both as listed). EHT Table 4's "−1, −2/3" for the MT wormhole therefore looks like a sign slip.
- **Do not confuse that wormhole with PV.** The excluded MT metric has the exponential redshift in *areal* radius (r_sh = e), whereas PV has it in *isotropic* radius (r_sh = 2e). The EHT's exclusion of the MT wormhole does not apply to PV.
- **What would decide it.** A ~1.5% (3σ) shadow-size measurement with spin marginalised. That needs a reduction of the α₁ (ring-to-shadow) calibration that GRMHD-calibrated ring sizes cannot deliver. A direct n = 1 or n = 2 photon-ring measurement (BHEX, lunar baselines) is the clean route. BHEX's geometric-only ~10% is not enough, but its source-modelled accuracy might be.

### Gaps
- I found no paper that fits the exponential/PV/Yilmaz metric to EHT visibilities or posteriors. All published statements are size comparisons.
- I could not open Turimov et al. 2022 or Mychelkin et al. 2024 (GRG) to check their shadow numbers.
- I have not checked an EHT-style analysis of the 2018 or 2021 M87* and Sgr A* campaigns for updated δ.

## 3. Strong deflection: Bozza coefficients for the exponential metric

### Takeaway
- **[computed] PV:** ā = 1 exactly, b̄ = −0.4658176, u_m = 2e.
- **GR, reproduced:** ā = 1, b̄ = ln[216(7 − 4√3)] − π = −0.4002300, u_m = 3√3.
- Because ā = 1 in both, the photon-ring Lyapunov exponent is γ = π in both, and successive subrings shrink by e^π = 23.14 in both. **Subring ratios cannot tell PV from GR; only the absolute size (u_m) and b̄ (subring offsets) can.**
- A 2023 paper by Bozza's method exists, but its numbers were not accessible.

### Cited Findings
- **Bozza 2002** ([arXiv:gr-qc/0208075](https://arxiv.org/abs/gr-qc/0208075), §2):
  - Metric ds² = A dt² − B dx² − C dΩ².
  - Photon sphere C′/C = A′/A; u_m = √(C_m/y_m) with y = A.
  - β_m = C_m(1 − y_m)²(C″_m y_m − C_m A″_m)/(2y_m²C′_m²).
  - R(z, x₀) = 2√(By)(1 − y₀)√C₀/(CA′).
  - ā = R(0, x_m)/(2√β_m) and b̄ = −π + b_R + ā log(2β_m/y_m), with b_R the subtracted regular integral.
  - The strong-field limit is α(θ) = −ā log(θD_OL/u_m − 1) + b̄.
  - Schwarzschild (§3.1): ā = 1, b_R = 2 log[6(2 − √3)] = 0.9496, b̄ = −0.4002, giving α = −log(2θD_OL/(3√3) − 1) + log[216(7 − 4√3)] − π.
- **The 2023 paper.** "Strong lensing in the exponential wormhole spacetimes" ([ScienceDirect S1384107623000416](https://www.sciencedirect.com/science/article/abs/pii/S1384107623000416)) applies Bozza's method to the exponential metric "for M = 0.6, 0.8, 0.9, and 1". Known only from the search snippet; the page returned HTTP 403 and no coefficient values were visible.
- **Boonserm et al. §6 (effective refractive index)** ([arXiv:1805.03781](https://arxiv.org/abs/1805.03781)):
  - Null geodesics are exactly those of the refractive index n(r) = e^{2m/r}, "well defined all the way down to r = 0".
  - The isotropic Schwarzschild index is n = (1 + m/2r)³/|1 − m/2r|.
  - "Strong field lensing phenomena are markedly different from Schwarzschild."
- **[computed] PV analytic pieces** (isotropic x = r): x_m = 2, y_m = e⁻¹, C_m = 4e, C′_m = 2e, C″_m = e, A′_m = 1/(2e), A″_m = −1/(4e).
  - So β_m = e(1 − 1/e)² = 1.086161 and R(0, x_m) = 2√e(1 − 1/e), giving **ā = 1**.
  - b_R = 0.8999782 and ln[2(e−1)²] = 1.7757969, so **b̄ = −π + b_R + ln[2(e − 1)²] = −0.46581756**.
- **[computed] Independent check.** α(b) + ln(b/b_c − 1), extrapolated to b → b_c from the exact integral, gives b̄ = −0.46581758 for PV and −0.40023006 for GR, agreeing with the formula to 2×10⁻⁸.
- **[computed] Image and subring structure:**
  - (b₄π − b_c)/(b₅π − b_c) = 23.140960 for PV and 23.140939 for GR, against e^π = 23.140693.
  - The photon-orbit Lyapunov exponent is λ = Ω_c exactly in both, so γ = λπ/Ω_c = π.
  - b(α=2π)/b_c − 1: GR 1.2563×10⁻³, PV 1.1768×10⁻³.
  - b(α=3π)/b_c − 1: GR 5.409×10⁻⁵, PV 5.066×10⁻⁵.
- **[computed] Bozza observables for Sgr A*** (VLTI θ_g = 5.124 µas):

  | | θ∞ | s | r |
  |---|---|---|---|
  | GR | 26.627 µas | 33.32 nas | 6.822 mag |
  | PV | 27.859 µas | 32.65 nas | 6.822 mag |

- **[computed] Test values for `exact_deflection`:**
  - b = 6M: GR 1.719388310, PV 1.999746113.
  - b = 10M: GR 0.590395788, PV 0.610437437.
  - b = b_c(1 + 10⁻³): GR 6.5106448, PV 6.4454619.
  - b = b_c(1 + 10⁻⁶): GR 13.4152856, PV 13.3496988.

### Inferences
- ā = 1 for PV is a closed-form result: C″_m A_m − C_m A″_m = 2 and A_m B_m = 1, the same as Schwarzschild. It means the logarithmic divergence rate and the e^{π} subring demagnification are identical.
- **Strong-lensing discriminators** are therefore the critical curve (+4.6%) and b̄ (−0.4658 against −0.4002). b̄ moves the n ≥ 1 subrings' offsets from the critical curve by about 6% (e^{Δb̄} = 0.936).
- **Library check.** These values make a sharp doc test for the repo's PV `exact_deflection` near b_c, which is currently tested only for "> 10".

### Gaps
- The 2023 New Astronomy paper's ā and b̄ could not be read, so I have no published PV b̄ to compare with −0.46582.
- The photon-ring (Johnson et al. 2020) statement that γ = π for Schwarzschild is taken from the repo's docstring, not from a source I opened. My computation reproduces it.

## 4. Horizon and wormhole nature: throat, energy conditions, ringdown, echoes

### Takeaway
- The exponential metric has no horizon. It is a traversable wormhole:
  - throat at isotropic r = M, areal radius eM;
  - photon sphere at isotropic r = 2M, areal 3.297M;
  - ISCO at isotropic (3 + √5)M = 5.236M, areal 6.338M.
- **Energy conditions.** It violates the null Ricci convergence condition. In GR language it is sourced by a phantom (negative-kinetic-energy) scalar.
- **No published QNMs.** No paper I found gives quasinormal modes of the exponential metric.
- **[computed] Test-field QNMs** (a fixed geometry with GR-like perturbation equations): the spin-2 axial ℓ = 2 mode is Mω = 0.3448 − 0.0851i, against Schwarzschild's 0.3737 − 0.0890i, i.e. −7.7% in frequency and −4.4% in damping. This is larger than the "few percent" GW250114 bound on the quadrupolar frequency.
- **No echoes.** Every potential I checked is a single barrier, and waves pass through the throat into the other universe. PV therefore predicts no echoes, so echo non-detections do not test it.

### Cited Findings
- **Boonserm et al.** ([arXiv:1805.03781](https://arxiv.org/abs/1805.03781)):
  - **Throat (§2).** A(r) = 4πr²e^{2m/r}, dA/dr = 8π(r − m)e^{2m/r}; "has a minimum at r = m", A″ = 8πe² > 0. The exponential metric "has no horizons, g_tt ≠ 0, and so is not a black hole", with the region r < m "an infinite-volume 'other universe' that exhibits the 'underhill effect'; time runs slower on the other side."
  - **Areal radius (§3.2).** The minimum is r_s = m e, double-covered. The explicit form is ds² = −e^{2W(−m/r_s)}dt² + dr_s²/(1 + W(−m/r_s))² + r_s²dΩ².
  - **Ricci convergence (§5).** R_ab = −(2m²/r⁴)diag{0,1,0,0}, "This clearly will not be satisfied for all timelike, null, or spacelike vectors"; "the violation of the null Ricci convergence condition is crucial for understanding the flare out".
  - **GR interpretation (§9).** G_ab = −½{∇_aΦ∇_bΦ − ½g_ab(∇Φ)²} with Φ = 2m/r, i.e. "a negative kinetic energy massless scalar field, a 'ghost' or 'phantom' field". Also "a special sub-case of the Brans class IV solution of Brans–Dicke theory".
  - **ISCO (§7).** r_ISCO = (3 + √5)m ≈ 5.236067977m; areal 6.337940263m; L_min ≈ 3.523216438m. They note "Schwarzschild ISCOs for massive particles at r_s = 6m have already been seen by astronomers; this might place interesting bounds".
  - **Regge–Wheeler potentials (§8, "inverse Cowling approximation… we do not a priori know the spacetime dynamics").**
    - Tortoise coordinate dr* = e^{2m/r}dr.
    - 𝒱_S = e^{−4m/r}[ℓ(ℓ+1) + (S² − 1)((1 − m/r)² − 1)]/r².
    - Spin 0, ℓ = 0: a peak at (1 + 1/√3)m and a trough at (1 − 1/√3)m.
    - Spin 1: a maximum at r = 2m of ℓ(ℓ+1)/(2me)².
  - **Their conclusion (§10):** "If one wishes to replace all the black hole candidates astronomers have identified with traversable wormholes, then certainly a careful phenomenological analysis of this quite radical proposal… should be carried out."
- **Mychelkin, Suliyeva & Makukov** ([arXiv:2510.15391](https://arxiv.org/abs/2510.15391), abstract):
  - At the throat "the topological Gauss–Bonnet invariant swaps sign, and the value of the Keplerian frequency for circular geodesics becomes singular".
  - The curvature invariants are extremal at r = M/2.
- **Makukov & Mychelkin 2018** ([arXiv:1809.05290](https://arxiv.org/abs/1809.05290)):
  - They "obtain the exact rotational generalization of the antiscalar Papapetrou spacetime as a viable alternative to the Kerr black hole".
  - ω_isco = 0.0633326/M against 0.0680414/M, "the 6.92% difference" (first noted by Watt & Misner), "characterizing the upper cut-off in the gravitational wave spectrum".
- **Turimov et al. 2022** (search abstract only, not opened):
  - The "ISCO… MBO and photon sphere… are less than that obtained in the Schwarzschild spacetime".
  - "energy efficiency… is relatively large than that in the Schwarzschild spacetime."
  - The exponential wormhole "is a special case of the well-known JNW spacetime".
- **GW250114** ([LVK, arXiv:2509.08099](https://arxiv.org/abs/2509.08099), abstract):
  - The ℓ = m = 2 fundamental and first overtone "match the Kerr prediction to tens of percent".
  - Their parameterised IMR fit bounds "the quadrupolar frequency to within a few percent of the GR prediction".
- **Schwarzschild reference values** ([Berti, Cardoso & Starinets 2009, arXiv:0905.2975](https://arxiv.org/abs/0905.2975), Pöschl–Teller discussion): Mω = 0.1105 − 0.1049i (ℓ = s = 0) and 0.3737 − 0.0890i (ℓ = s = 2), from the "numerical result" (Leaver; Berti).
- **[computed] Circular orbits:**
  - **ISCO.** PV MΩ = 0.0633326 (Makukov's value reproduced), E = 0.9452132, L = 3.523216 (Boonserm's value reproduced). GR MΩ = 0.0680414, E = 0.9428090, so PV is **−6.920%** in Ω_ISCO.
  - **Efficiency.** 1 − E_ISCO is **5.479% (PV) against 5.719% (GR)**, i.e. *lower* in PV, contrary to Turimov's abstract.
  - **f_GW(ISCO) for 10 M_⊙:** GR 439.7 Hz, PV 409.3 Hz.
  - **Static-emitter redshift:** at the photon sphere 1 + z = 1.6487 (PV) against 1.7321 (GR); at the ISCO 1.2104 against 1.2247.
- **[computed] Test-field ringing.** Time-domain evolution on the r* line with Boonserm's potentials, fitted over windows starting at 60–100 M; output in `pv_qnm_check_output.txt`.
  - **Validation.** GR spin-2 ℓ = 2 gives 0.3737 − 0.0890i, matching Berti et al.
  - **PV results:**
    - spin-2 ℓ = 2: **0.3448 − 0.0851i** (stable to 0.0004 across windows);
    - spin-0 ℓ = 2: 0.4638 − 0.0924i, against GR 0.4834 − 0.0968i;
    - spin-1 ℓ = 1: about 0.236 − 0.086i, against GR about 0.248 − 0.092i (noisier).
  - **Potential shape.** All six potentials have one maximum; PV's spin-2 ℓ = 2 peak is at r = 2.44 (isotropic).
  - **Eikonal limit.** Since λ = Ω_c in both metrics, PV/GR = 3√3/(2e) = 0.9558 for both the real and imaginary parts.

### Inferences
- **Horizon status: decisive in principle, but no observation yet probes the horizon itself in a way that excludes PV.**
- **EHT.** The images need a photon sphere and a dark interior. PV has the photon sphere, and infalling matter crosses the throat rather than hitting a surface, so the arguments in the EHT section "Are there Viable Alternatives to an Event Horizon?" about thermalising and reflecting surfaces do not obviously apply.
- **Ringdown proxy.** If PV's remnant rang like a GR spin-2 field on the static exponential metric, a 7.7% frequency shift at fixed mass would exceed GW250114's few-percent bound. But:
  1. GW remnants are spinning, and no QNMs of the rotating Papapetrou metric exist.
  2. In PV *theory* the radiation is a single scalar, longitudinal mode (§5). So this proxy tests a hybrid "PV geometry + GR waves", not PV.
- **Echoes.** With a single barrier and the potential vanishing as r* → −∞ (r → 0), there is no cavity, so no echoes and no long-lived trapped modes in these channels. Echo searches cannot test PV.
- **Accretion.** The ISCO frequency (−6.9%) and efficiency (5.48% against 5.72%) are degenerate with a slightly retrograde Kerr spin in X-ray continuum fitting or QPO work. This is not a clean test without an independent spin.

### Gaps
- No published QNM spectrum (test-field or gravitational) of the exponential or rotating Papapetrou metric was found. Mine are test-field only.
- Gravitational perturbations of PV theory (the nonlinear K field, Puthoff eq 59) have never been computed, as far as I found.
- I did not open LVK echo-search or GWTC-3 ringdown papers. Turimov et al.'s definition of "efficiency" is unknown, so the discrepancy with my 5.48% is unresolved.

## 5. Dynamics and radiation: binary pulsars, speed and polarisation of gravity, critiques

### Takeaway
- **Radiation.** PV theory radiates a single longitudinal scalar mode, with monopole and dipole terms absent. Its binary power is **2/3 of GR** (Ibison eq 65), at most 1.5% eccentricity-dependent:
  - PSR B1913+16 predicts 0.6732 against the measured 0.9983 ± 0.0016, about **200σ off**.
  - The double pulsar (validated to 1.3×10⁻⁴ at 95%) puts it about **5000σ off**.
- **Speed of gravity.** K-waves travel at c/K, exactly the speed of light, so GW170817's speed bound is passed trivially.
- **Polarisation.** The pure-scalar polarisation is disfavoured by Bayes factor >1000 (GW170814).
- **Frame dragging [inference].** PV has no g₀ᵢ, which forces PPN α₁ = −8 and **zero Lense–Thirring precession**, against GP-B (20%) and LAGEOS (10%), and |α₁| < 10⁻⁴.
- **Misner's critique.** "Yilmaz cancels Newton" targets Yilmaz's field equations, not Puthoff's Lagrangian.

### Cited Findings
- **Puthoff 2002** ([arXiv:gr-qc/9909037](https://arxiv.org/abs/gr-qc/9909037)):
  - **Lagrangian (eqs 54–55).** The K-field Lagrangian must take "the standard Lorentz-invariant form for propagational disturbances of a scalar": L_K = −(λ/K²)[(∇K)² − (1/(c/K)²)(∂K/∂t)²], with λ = c⁴/32πG and f(K) = 1/K².
  - **Equations.** Eq 58 adds a "dielectric force ∝ ∇K"; eq 59 is the field equation (∇² − (K/c)²∂²_t)√K = −(K/4λ)[mass, EM and K-field energy densities]. The static solution (eqs 61–62) is K = e^{2GM/rc²}.
  - **His own caveat:** "important topics such as gravitational radiation and frame-dragging effects were not addressed."
- **Ibison 2003, "An Investigation of the Polarizable Vacuum Cosmology"** ([arXiv:astro-ph/0302273](https://arxiv.org/abs/astro-ph/0302273), sections "Gravitational radiation" and "Developments"; no journal reference):
  - The PV metric "always has the form ds² = c²dt²/K − K dx²", with the K-field action "of a scalar wave equation with speed c/K".
  - **Multipoles.** The monopole and dipole terms are "proportional respectively to the total energy and momentum… constant in time and cannot contribute to radiation".
  - **Polarisation.** "the radiation is longitudinally polarized", with "just one polarization degree of freedom".
  - **Eq 65.** P_PV/P_GR = ½(51e⁴ + 396e² + 128)/(37e⁴ + 292e² + 96). "exactly 2/3" at e = 0, changing "only by 1.5% across the full range of eccentricity".
  - **Conclusion.** Since PSR 1913+16 agrees with GR "to better than 1%, … it follows that PV is incompatible with observation". He suggests vacuum anisotropy or the Yilmaz theory might restore two transverse polarisations.
- **Hulse–Taylor** ([Weisberg & Huang 2016, arXiv:1606.02744](https://arxiv.org/abs/1606.02744)): observed/GR orbital decay is "0.9983 ± 0.0016"; e = 0.6171340(4) (table "Orbital Parameters").
- **Double pulsar** ([Kramer et al. 2021, arXiv:2112.06795](https://arxiv.org/abs/2112.06795), PRX 11, 041050): e = 0.088; GR's quadrupole formula validated "at a level of 1.3×10⁻⁴ (95% conf.)".
- **GW170817** ([LVK+Fermi+INTEGRAL, arXiv:1710.05834](https://arxiv.org/abs/1710.05834), ApJL 848, L13): the speed of gravity minus the speed of light lies "between −3×10⁻¹⁵ and +7×10⁻¹⁶ times the speed of light".
- **GW170814** ([arXiv:1709.09660](https://arxiv.org/abs/1709.09660), PRL 119, 141101): "Bayes' factors of more than 200 and 1000 in favor of the purely tensor polarization against purely vector and purely scalar, respectively."
- **Misner, "Yilmaz Cancels Newton"** ([arXiv:gr-qc/9504050](https://arxiv.org/abs/gr-qc/9504050), Nuovo Cim. B114, 1079, 1999): with Yilmaz's G^{μν} = 8πG(T^{μν} − t^{μν}) and the Bianchi identities, "the Euler hydrodynamic equation in the Newtonian limit is modified to remove all gravitational forces", allowing "a solar system in which the Sun and the planets are permanently at rest".
- **Reply** ([Alley, Aschan & Yilmaz, arXiv:gr-qc/9506082](https://arxiv.org/abs/gr-qc/9506082)): Misner's claim "is based on a false premise. With the correct premise the conclusion… applies to general relativity and not to the Yilmaz theory."
- **Watt & Misner 1999** ([arXiv:gr-qc/9910032](https://arxiv.org/abs/gr-qc/9910032), abstract): a theory whose "spacetime metric is derived from a single scalar field Φ" with a "non-linear flat-space four-dimensional wave equation" is used as a numerical-relativity test bed. It "does not generate results which are exactly identical quantitatively to those of general relativity". Makukov & Mychelkin cite it for the exponential metric's ISCO.
- **PPN framework** ([Will 2014, arXiv:1403.7377](https://arxiv.org/abs/1403.7377)):
  - PPN metric: g₀ᵢ = −½(4γ + 3 + α₁ − α₂ + ζ₁ − 2ξ)Vᵢ − ½(1 + α₂ − ζ₁ + 2ξ)Wᵢ − ½(α₁ − 2α₂)wⁱU − α₂wʲU_ij.
  - Eq 42 (Lense–Thirring): Ω_LT = −½(1 + γ + α₁/4)[J − 3n(n·J)]/r³.
  - Gravity Probe B agreed with GR to "20 percent for the frame-dragging effect"; LAGEOS with GRACE/CHAMP gave "a 10 percent confirmation".
  - Table "Current limits on the PPN parameters": α₁ < 10⁻⁴ (LLR) and 7×10⁻⁵ (PSR J1738+0333); α₂ < 2×10⁻⁹; ξ < 4×10⁻⁹.
- **[computed] Ibison's ratio.** 0.66667 (e = 0), 0.67321 (B1913+16), 0.66693 (J0737−3039). The eq 64 normalisation reproduces Peters–Mathews (11.853275 at e = 0.6171). That is 203σ from 0.9983 ± 0.0016 and about 5000σ from the double-pulsar bound.

### Inferences
- **Radiation status: PV theory is excluded.** This is the cleanest published test. The deficit is a factor 2/3 in a quantity measured to 10⁻³–10⁻⁴, and it does not depend on any strong-field assumption. It applies to Puthoff's or Dicke's scalar theory, not to the exponential metric as a geometry (the GR + phantom-scalar reading radiates like GR).
- **Speed of gravity: passes.** PV's K-wave speed equals the light speed by construction (eq 55 against eq 24, v_L = c/K), so GW170817's bound is satisfied identically. The polarisation content (one scalar mode) is what fails: GW170814, BF > 1000.
- **Frame dragging: my new constructed test.**
  - **Argument.** If PV's metric keeps the diagonal form c²dt²/K − K dx² for moving or rotating sources, as Ibison states and Puthoff's eqs 17 and 59 imply (there is no vector potential), then g₀ᵢ ≡ 0 in PV's frame. Matching to Will's g₀ᵢ, the sum of the Vᵢ and Wᵢ coefficients is gauge-invariant under t → t + λχ̇. It gives 4γ + 4 + α₁ = 0, i.e. **α₁ = −8**, and α₂ = −1 + ζ₁ − 2ξ.
  - **Consequences.** Ω_LT ∝ (1 + γ + α₁/4) = 0, so PV predicts **no Lense–Thirring precession**. GP-B (20%) and LAGEOS (10%) then exclude PV at about 5σ and 10σ. |α₁| < 10⁻⁴ excludes α₁ = −8 by about 10⁵. If ζ₁ and ξ are small, α₂ ≈ −1 against 2×10⁻⁹.
  - **Caveat.** This needs PV to have a preferred frame, which it does as formulated in a background Euclidean space with speed c/K. It is my derivation, not a published one.
- **Misner's critique** concerns Yilmaz's modified Einstein equations. Puthoff's PV has a well-defined action, so the Bianchi-identity argument does not carry over, and the dispute (Misner against Alley, Aschan & Yilmaz) is unresolved in the literature I read. PV's decisive failures are radiation and gravitomagnetism instead.

### Gaps
- **Dipole radiation from strongly self-gravitating bodies.** PV's analogue of scalar "sensitivities" (body-dependent K charges) has not been computed, so I cannot quantify it for neutron-star–white-dwarf systems.
- **Full PPN parameter set.** The ζᵢ, ξ and the Nordtvedt parameter η_N of PV theory are not published. My α₁ result uses only g₀ᵢ.
- **Yilmaz radiation.** I did not find a computation of binary-pulsar radiation in the full Yilmaz tensor theory. Ibison only conjectures it could fix the 2/3.

## 6. Cosmology and laboratory analogues

### Takeaway
- **PV cosmology:**
  - It has no big bang, starting from a non-zero scale factor.
  - It cannot produce BBN light elements or the CMB as thermal relics.
  - It is otherwise tunable to the age and q₀ (Ibison 2003).
- **Yilmaz cosmology:** consistent only with an empty universe or a decaying negative-energy vacuum (Ibison 2006). This is contested by Robertson's 2015 SN Ia fit.
- **Lab analogues.** Because PV's null geodesics are *exactly* rays in a medium with n = e^{2GM/rc²}, a graded-index or waveguide analogue (as in Sheng et al. 2013 for GR-like profiles) can verify PV's light-bending mathematics, b_c = 2e·a and b̄ = −0.4658. It cannot tell whether spacetime is PV.

### Cited Findings
- **Ibison 2003** ([arXiv:astro-ph/0302273](https://arxiv.org/abs/astro-ph/0302273), abstract and "Summary of cosmological findings"):
  - "a big bang is not a solution: the initial scale factor may be small but not zero, and must have zero first derivative".
  - The Friedmann-like equation has "the ȧ² term… multiplied by −1/3".
  - "PV can explain neither present proportions of light elements nor the cosmic microwave background as the result of high temperatures in the very early history".
  - It "predicts a maximum initial temperature that is too small for BBN".
  - Predictions are "compatible with current estimates for the age and deceleration parameter, and the maximum observable redshift, provided the matter density Ω_m is close to the luminous matter density".
- **Ibison 2006** ([arXiv:0705.0080](https://arxiv.org/abs/0705.0080), CQG 23, 577):
  - Flat Yilmaz cosmology "is consistent only with either a completely empty universe, or with a negative energy vacuum that decays to produce a constant density of matter".
  - "The theory does not support the Cosmological Principle for curved space K ≠ 0".
- **Robertson 2015** ([arXiv:1507.07809](https://arxiv.org/abs/1507.07809), abstract; no journal reference):
  - "A redshift-luminosity distance relation in excellent agreement with observations is calculated here for SNe 1a using the Yilmaz gravitational theory". It needs no Λ but "still requires the existence of non-baryonic dark matter".
  - This conflicts with Ibison 2006's flat-space result, and I did not resolve the conflict.
- **Boonserm et al. §6** ([arXiv:1805.03781](https://arxiv.org/abs/1805.03781)): after dropping the conformal factor, dŝ² = −e^{−4m/r}dt² + dx², so n(r) = e^{2m/r}. Via Fermat's principle it "completely characterizes the focussing/defocussing of null geodesics", and "falls naturally into the 'analogue spacetime' programme".
- **Sheng et al. 2013**, "Trapping light by mimicking gravitational lensing" ([arXiv:1309.7706](https://arxiv.org/abs/1309.7706); the arXiv record has no journal reference, and I did not verify the published version):
  - A "microstructured optical waveguide around a microsphere" experimentally demonstrates "far-field gravitational lensing effects and the critical phenomenon in close proximity to the photon sphere".
  - Their method matches "the effective refractive index of the metric n = √(g_rr/g₀₀)" to the measured profile.
  - Table 1 gives a photon-sphere radius of 28.5 µm (AFM) and 27.96 µm (interference); the critical impact parameter is √2 times the photon-sphere scale for their fitted profile.
- **Genov, Zhang & Zhang 2009**, "Mimicking celestial mechanics in metamaterials", Nature Physics 5, 687–692, [doi:10.1038/nphys1338](https://doi.org/10.1038/nphys1338). Bibliographic details only, from search results. The Nature, HKU and LBL pages all returned 403 or a login redirect, so I could not check its index profile or whether a PV-like exp(2a/r) profile was considered.

### Inferences
- **Cosmology status: PV fails** by giving no hot early universe (no BBN, no thermal CMB). This is a strong, published, PV-specific failure. Yilmaz cosmology is contested.
- **Analogue status: a verification of mathematics, not of nature.** Any medium with n(r) = n∞·exp(2a/r) reproduces the full PV null-geodesic problem exactly. A laboratory analogue could check:
  - photon sphere at r = 2a;
  - critical impact parameter 2e·a;
  - ā = 1 and b̄ = −0.4658;
  - deflection series 4x + 4πx² + 48x³ with x = a/b.
- **Ease against GR.** It is easier than a Schwarzschild analogue, because n stays finite outside r = 0, whereas the isotropic Schwarzschild index diverges at r = m/2.
- **Index contrast needed.** e ≈ 2.72 at the photon sphere and e² ≈ 7.4 at the throat, relative to n∞. Waveguide effective-index engineering (Sheng-type) or transformation optics would be needed to reach the photon sphere.
- **Design-scale test.** The weak-field 2PN term (π/4)(a/b)², the PV − GR difference, is 6.7% of the second-order deflection, which is large enough for an analogue.

### Gaps
- I did not find any published analogue (optical, acoustic or waveguide) built specifically for the exponential or PV index.
- I did not check Puthoff's later PV papers (for example the 2010 metric-engineering papers) for a PV cosmology different from Dicke's and Ibison's. No PV CMB or BAO analysis exists that I found.

## 7. Implementation: what the Rust library has and should compute

### Takeaway
- `gravitics::spacetime` already covers 2PN bending (`light_bending`, 4π vs 15π/4) and the exact PV deflection (`exact_deflection`, b_c = 2e).
- It lacks:
  - a metric-aware shadow;
  - strong-deflection coefficients;
  - ISCO and periastron functions;
  - the radiation ratio;
  - the third-order bending check for PV (48x³).
- Each item below gives the test value computed here, so it can be a doc test, per the repo's tested-examples rule.

### Cited Findings
- **Existing code** (`crates/gravitics/src/spacetime.rs`):
  - `Metric::{GeneralRelativity, PolarizableVacuum}`.
  - `light_bending` (2PN).
  - `exact_deflection` (Darwin integral; PV via F = e^{4u}/b² − u²; returns `None` inside b_c).
  - `shadow_diameter(mass, distance)`, GR-only 6√3.
  - `shadow_deviation`.
  - `kerr_pole_on_shadow`.
  - The existing tests check the PV limb difference, PV capture at 2e and PV > GR at 10 rg.
- **Bozza's formulas for any static spherical metric** ([arXiv:gr-qc/0208075](https://arxiv.org/abs/gr-qc/0208075)) and **Ibison's eq 65** ([arXiv:astro-ph/0302273](https://arxiv.org/abs/astro-ph/0302273)) are the closed forms to implement.
- **[computed] Test values** (units GM/c² unless stated):

  | # | Function (proposed unless noted) | GR | PV |
  |---|---|---|---|
  | 1 | `critical_impact_parameter(metric)`, used by a metric-aware `shadow_diameter(mass, distance, metric)` | √27 = 5.196152 | 2e = 5.436564, δ = +0.046267 |
  | 1a | Sgr A* (4.297×10⁶ M_⊙, 8277 pc) | 53.254 µas | 55.718 µas |
  | 2 | `strong_deflection(metric) -> (ā, b̄)` | (1, ln[216(7−4√3)]−π = −0.40023004) | (1, −0.46581756) |
  | 2a | Doc test: `exact_deflection` at b_c(1+10⁻⁶) | 13.4152856 | 13.3496988 |
  | 2b | Doc test: residual α + ln ε at ε = 10⁻⁶ (O(ε ln ε) offset ≈ 6×10⁻⁶ from b̄) | −0.400225 | −0.465812 |
  | 3 | Photon-ring subring ratio (b₄π−b_c)/(b₅π−b_c), reusing the existing doc test with `Metric::PolarizableVacuum` | 23.140939 | 23.140960 |
  | 4 | Far-field series as a doc test on `exact_deflection` (3rd term in x³) | 128/3 | 48 |
  | 5 | `isco(metric) -> (areal radius, MΩ, E)`, with efficiency 1 − E | (6, 6^{−3/2} = 0.0680414, √(8/9) = 0.9428090), 5.719% | (6.337940, 0.0633326, 0.9452132), 5.479% |
  | 6 | `periastron_ratio(metric, r)` = Ω_φ/Ω_r for circular orbits | (1 − 6x)^{−1/2}: 1.58114 at x = 0.1 | series 1 + 3x + 14.5x² + 76.5x³; isotropic r = 9 → x = 0.09965, K = 1.61645; r = 19 → x = 0.04996, K = 1.19927 |
  | 7 | `surface_redshift(metric, GM/(c²R_areal))`; PV needs u e^{−u} = U_A on the W₀ branch, z = e^u − 1 | U_A = 0.172 → 0.2346620 | 0.2371177 |
  | 8 | `shapiro_2pn_difference(mass, b)`, fixed endpoints: (π/4)(GM/c²)²/(bc) | | 8.211 ps at the solar limb |
  | 9 | `pv_radiated_power_ratio(e)`, Ibison eq 65 | | 2/3 (e = 0), 0.67321 (e = 0.617134), 0.66693 (e = 0.0878) |
  | 10 | `regge_wheeler(metric, s, l, r)`, Boonserm §8 potentials | ℓ = 2, s = 2 peak at r = 3.28 | ℓ = 2, s = 2 peak 0.13188 at r = 2.4365 |
  | 10a | Eikonal QNM ω ≈ [(l + ½) − i(n + ½)]/b_c (valid for both, since λ = Ω_c) | | |
  | 10b | Time-domain QNMs, kept as an example, not library | 0.3737 − 0.0890i | 0.3448 − 0.0851i |

  - Row 1: to keep the change small, add `critical_impact_parameter` and leave the GR `shadow_diameter` in place.
  - Row 3: the existing doc-test pattern works unchanged with `Metric::PolarizableVacuum`.
  - Row 4: `light_bending` plus 48x³ should match `exact_deflection` far out, just as GR does with 128/3.
  - Row 5 requires PV's (3 + √5) ISCO.
  - Row 6 (GR) matches `(1−6x)^(−1/2)`, e.g. 1.58114 at x = 0.1.
- **[computed] Exact deflection values** for further doc tests:
  - b = 7M: GR 1.127639105, PV 1.216647657.
  - b = 20M: GR 0.236135995, PV 0.239033624.
  - b = 100M: GR 0.041222540, PV 0.041306741.

### Inferences
- **Scoreboard** (both ways computed; status labels per the repo's fringe-physics rule):

  | Test | GR | PV | Data | Status |
  |---|---|---|---|---|
  | 1PN: redshift, bending 1.7512″, Mercury 42.98″/cy, Shapiro, Cassini | ✓ | ✓, identical | Will 2014 | Both pass |
  | 2PN solar bending (10.95 vs 11.68 µas) and Shapiro (8 ps) | ✓ | not yet tested | Cassini, VLBI below by 10²–10³ | Undecided; LATOR-class mission → ~33σ |
  | S2 2PN precession | ✓ | Δ = 2.7×10⁻⁷ rad per orbit | | Undecided (10⁻³–10⁻⁴ of σ) |
  | EHT Sgr A* shadow (+4.6%) | | 1.0–1.6σ off (2.0σ worst row) | GR within 0.4–0.9σ | Mildly disfavoured |
  | M87* shadow | | 0.33σ | | Not constraining |
  | Strong lensing ā | 1 | 1 | | Indistinguishable |
  | Strong lensing b̄ | −0.400 | −0.466 | | Untested |
  | ISCO Ω (−6.9%), efficiency (5.48 vs 5.72%) | | | | Degenerate with spin; not a clean test yet |
  | Ringdown | | test-field −7.7% | | Proxy only; PV theory has no tensor ringdown |
  | Binary-pulsar decay | 1 | 2/3 | 0.9983 ± 0.0016; 1.3×10⁻⁴ | **PV theory excluded** (~200σ; ~5000σ) |
  | Polarisation | | pure scalar | GW170814, BF > 1000 | **Disfavoured** |
  | GW speed | | = c | GW170817 | Pass |
  | Frame dragging (my inference) | | zero; α₁ = −8 | GP-B 20%, LAGEOS 10%, α₁ < 10⁻⁴ | **PV theory excluded if the preferred-frame argument holds** |
  | Cosmology | | no hot big bang, no BBN or CMB | Ibison 2003 | **PV cosmology fails** |

- **Where PV stands, honestly.**
  - The *exponential metric* as a static geometry, sourced in GR by a phantom scalar, is alive. It is mildly disfavoured (about 1–2σ) by the Sgr A* shadow, and first-class discriminators exist: shadow precision near 1.5%, b̄, and 2PN solar bending.
  - Puthoff's *PV theory* (scalar K-field dynamics) is excluded by binary pulsars and, by my derivation, by frame dragging, independently of any strong-field question.
  - The way to "verify" PV is to compute these numbers, which this file does, and to fund the two experiments that can still discriminate the geometry: percent-level horizon-scale imaging, and second-order solar light deflection.
- **Per the "rigs are examples" rule.** The time-domain QNM fitter and the analogue-waveguide design are examples. The closed forms (b_c, Bozza coefficients, ISCO, Ibison's ratio, surface redshift) belong once in the library.

### Gaps
- A rotating-PV shadow, needed to compare with spin-marginalised EHT posteriors, requires Makukov & Mychelkin's rotating Papapetrou metric, which I did not implement or check.
- I have not verified whether the repo's `Metric` enum should distinguish "PV geometry" from "PV dynamics". The radiation and frame-dragging results apply only to the latter.
