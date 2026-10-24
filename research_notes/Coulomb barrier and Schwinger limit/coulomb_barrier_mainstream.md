# Getting nuclei through the Coulomb barrier: mainstream formulas, published test values and experiments

Conventions used throughout: E is centre-of-mass (CM) energy unless marked "lab"; T is temperature in keV (kT); ⟨σv⟩ is in cm³/s unless marked m³/s.
"(computed)" means recomputed for these notes with `python3 -I` using CODATA-2018 constants (α = 7.2973525693e-3, ħc = 197.3269804 MeV fm, so e²/4πε₀ = 1.439965 MeV fm; m_p = 938272.088, m_d = 1875612.942, m_t = 2808921.133, m_³He = 2808391.607 keV).
Scripts: `/tmp/claude-1000/-home-dan-code-FusionVault-quantum-gravimeter/7beb35b3-c723-440c-9450-66a991641392/scratchpad/calc/coulomb.py` and `.../calc/pp.py` (scratch, session-local).
Labels: [mainstream], [contested], [refuted].

## Q1. Coulomb barrier height, classical turning point, Sommerfeld parameter, Gamow factor, Gamow energy, Gamow peak

### Takeaway
The barrier (~0.36–0.48 MeV for D–T, ~1.6–1.9 MeV for p–¹¹B, depending on the radius convention) is far above plasma energies of 10–15 keV. Fusion happens only by tunnelling. The tunnelling probability is P = exp(−2πη) = exp(−√(E_G/E)), and the Maxwellian-weighted rate peaks at E₀ = (E_G (kT)²/4)^{1/3}. For D–T at 10 keV that gives E₀ ≈ 30.9 keV and a 1/e width ≈ 40.6 keV. For p+p at the solar centre (1.55×10⁷ K) it gives E₀ ≈ 6.04 keV and width ≈ 6.56 keV (all computed). Every element of this is [mainstream].

### Cited Findings
- [mainstream] Cross section factorisation: σ(E) = S(E)/E · exp[−2πη(E)], with Sommerfeld parameter η(E) = Z₁Z₂α/v and v = √(2E/μ) (ħ = c = 1). The S-factor is slowly varying for non-resonant reactions. — [Adelberger et al., Solar Fusion II, Rev. Mod. Phys. 83, 195 (2011), Eq. 7](https://arxiv.org/abs/1004.2318); same form with η(E) = Z₁Z₂α√(μ/2E) and P(E) ≡ exp(−2πη(E)) in [Acharya et al., Solar Fusion III (arXiv:2405.06470v3, 2025), Eq. 24](https://arxiv.org/abs/2405.06470).
- [mainstream] Gamow constant as used by Bosch & Hale: B_G = πα Z₁Z₂ √(2 m_r c²) [keV^{1/2}], so 2πη = B_G/√E and the Gamow energy E_G = B_G² = 2 m_r c² (παZ₁Z₂)². — [Heltemes, Moses & Santarius, UWFDM-1268 (2005), Eq. 2](https://fti.neep.wisc.edu/fti.neep.wisc.edu/pdf/fdm1268.pdf)
- [mainstream] Published B_G and m_r c² (Bosch & Hale 1992). D–T (³H(d,n)⁴He): B_G = 34.3827 keV^{1/2}, m_r c² = 1124656 keV. D–D (both branches): 31.3970, 937814. D–³He: 68.7508, 1124572. — [MFE Formulary notebook (Bosch–Hale implementation)](https://notebook.community/MFEFormulary/MFEFormulary/ipython/BH_ReactionRates); [UWFDM-1268 Table 5](https://fti.neep.wisc.edu/fti.neep.wisc.edu/pdf/fdm1268.pdf)
- [mainstream] Gamow peak (SF-II Eqs. 8–9):
  ```
  E0/kT  = (π Z1 Z2 α / √2)^(2/3) · (μ/kT)^(1/3)          (natural units)
  ΔE0/kT = 4 · sqrt(E0 / (3 kT))                          (full width at 1/e, Gaussian approx.)
  <σv>   = sqrt(2/(μ kT)) · (ΔE0/kT) · f0 · Seff · exp(−3E0/kT)
         = 1.301e-14 cm^3/s · (Z1 Z2 / A)^(1/3) · f0 · (Seff / MeV b) · T9^(−2/3) · exp(−3E0/kT),   A = A1A2/(A1+A2)
  Seff   = S(0)[1 + 5kT/(36E0)] + S'(0) E0 [1 + 35kT/(36E0)] + ½ S''(0) E0² [1 + 89kT/(36E0)]
  ```
  — [Solar Fusion II, Eq. 8](https://arxiv.org/abs/1004.2318)
- [mainstream] Engineering approximations: E₀ ≈ 0.1220 (Z₁²Z₂² μ T₉²)^{1/3} MeV and ΔE₀ ≈ 0.2368 (Z₁²Z₂² μ T₉⁵)^{1/6} MeV, with μ in u and T₉ in GK. These are valid only for non-resonant σ or broad resonances. — [Rauscher, Phys. Rev. C 81, 045807 (2010)](https://elog.gsi.de/esr/E127/190925_135951/PhysRevC.81.045807.pdf) (seen via search summary); also used in [arXiv:2512.01846](https://arxiv.org/pdf/2512.01846)
- [mainstream] The Solar Fusion II reference temperature for the Sun's centre is 1.55×10⁷ K. pp S-factor: S(0) = (4.01 ± 0.04)×10⁻²² keV b and S′(0) = (4.49 ± 0.05)×10⁻²⁴ b. — [Solar Fusion II, Table I](https://arxiv.org/abs/1004.2318)
- [mainstream] Textbook barrier test value: treating D and T as touching spheres of radius 2.0 fm (separation 4.0 fm) gives a barrier of 5.76×10⁻¹⁴ J ≈ 360 keV. — [NCERT-style worked problem (esaral)](https://www.esaral.com/q/consider-the-d-t-reaction-deuterium-tritium-fusion-56156). A university lecture also gives ≈0.36 MeV for D–T and notes the D–T cross section peaks near 110 keV (lab), well below the barrier. — [Lehigh ME362 lecture 3](https://www.lehigh.edu/~eus204/teaching/ME362/lectures/lecture03.pdf) (via search summary)
- Conflicting loose statements exist. Wikipedia's fusion article says "the energy necessary to overcome the Coulomb barrier is 0.1 MeV" for D–T, and its colliding-beam article gives ~100 keV. Neither is a barrier height; both are roughly the energy where σ_DT becomes large. — [Wikipedia: Nuclear fusion](https://en.wikipedia.org/wiki/Nuclear_fusion); [Wikipedia: Colliding beam fusion](https://en.wikipedia.org/wiki/Colliding_beam_fusion)

### Inferences
- Formulas a Rust implementation needs (units explicit):
  ```
  V_C [MeV]    = 1.439965 · Z1 Z2 / R[fm],  R = r0 (A1^(1/3) + A2^(1/3))  (convention!)
  r_tp [fm]    = 1.439965 · Z1 Z2 / E[MeV]                 classical turning point
  η            = Z1 Z2 α sqrt(μc² / 2E)                     (μc², E in the same unit)
  2πη          = B_G / sqrt(E[keV]),   B_G = π α Z1 Z2 sqrt(2 μc²[keV])
  E_G          = B_G² = 2 μc² (π α Z1 Z2)²
  E0           = (E_G (kT)² / 4)^(1/3);   Δ = 4 sqrt(E0 kT / 3);   depth of Gamow suppression = exp(−3E0/kT)
  ```
  The E₀ form is algebraically identical to SF-II Eq. 8 (computed check: same numbers to ≤10⁻⁴).
- Reduced masses and Gamow constants from CODATA-2018 masses (computed) vs Bosch–Hale published values:

  | pair | m_r c² computed (keV) | BH m_r c² | B_G computed (keV^½) | BH B_G | E_G computed (keV) |
  |---|---|---|---|---|---|
  | D–D | 937806.5 | 937814 | 31.3969 | 31.3970 | 985.8 |
  | D–T | 1124647.3 | 1124656 | 34.3826 | 34.3827 | 1182.2 |
  | D–³He | 1124562.5 | 1124572 | 68.7626 | 68.7508 | 4728.3 |
  | p–¹¹B | 859604.5 | — | 150.2967 | — | 22589.1 |
  | p–p | 469136.0 | — | 22.2065 | — | 493.1 |

  For D–³He the computed B_G differs from BH by 0.017%, presumably from older constants or masses. For bit-for-bit reproduction of BH tables, use BH's tabulated B_G and m_r c², not CODATA-derived ones.
- Barrier heights (computed; convention-dependent):

  | pair | r0 = 1.2 fm: R, V_C | r0 = 1.4 fm: R, V_C | NCERT 2 fm spheres |
  |---|---|---|---|
  | D–D | 3.024 fm, 476 keV | 3.528 fm, 408 keV | 360 keV |
  | D–T | 3.243 fm, 444 keV | 3.783 fm, 381 keV | 360 keV (published test value) |
  | D–³He | 3.243 fm, 888 keV | 3.783 fm, 761 keV | 720 keV |
  | p–¹¹B | 3.869 fm, 1861 keV | 4.514 fm, 1595 keV | — |
- Turning point and tunnelling factor P = exp(−2πη) (computed). r_tp is 1440·Z₁Z₂ fm at 1 keV, 144·Z₁Z₂ fm at 10 keV and 14.4·Z₁Z₂ fm at 100 keV.

  | pair | η, P at 1 keV | η, P at 10 keV | η, P at 100 keV |
  |---|---|---|---|
  | D–D | 4.9970, 2.315e-14 | 1.5802, 4.876e-5 | 0.4997, 4.330e-2 |
  | D–T | 5.4722, 1.169e-15 | 1.7304, 1.897e-5 | 0.5472, 3.212e-2 |
  | D–³He | 10.944, 1.370e-30 | 3.4608, 3.601e-10 | 1.0944, 1.032e-3 |
  | p–¹¹B | 23.921, 5.333e-66 | 7.5643, 2.285e-21 | 2.3920, 2.970e-7 |
- Gamow peaks (computed, exact μ; values from the 0.1220/0.2368 approximation with integer-A μ in brackets):

  | case | kT | E₀ (keV) | Δ (keV) | exp(−3E₀/kT) |
  |---|---|---|---|---|
  | p+p, Sun centre 1.55×10⁷ K | 1.336 keV | 6.036 [6.020] | 6.557 [6.549] | 1.29e-6 |
  | p+p, 1.57×10⁷ K | 1.353 keV | 6.088 [6.071] | 6.628 [6.619] | 1.37e-6 |
  | D–T | 10 keV | 30.92 [30.84] | 40.61 [40.56] | 9.37e-5 |
  | D–T | 15 keV | 40.51 [40.42] | 56.93 [56.87] | 3.03e-4 |
  | D–D | 10 keV | 29.10 | 39.40 | 1.62e-4 |
  | D–³He | 10 keV | 49.08 | 51.16 | 4.04e-7 |
  | p–¹¹B | 15 keV | 108.3 | 93.1 | 3.91e-10 |
  | p–¹¹B | 300 keV | 798 | 1130 | 3.42e-4 |

  At T₆ = 15 and μ = 0.5 the approximation gives the commonly quoted ≈5.9 keV for p+p (computed: 0.122·(0.5·0.015²)^{1/3} MeV = 5.89 keV).
- Worked solar rate (computed, inputs from SF-II Table I, f₀ = 1.049 from Q3, T = 1.55×10⁷ K): S_eff = 4.463×10⁻²² keV b and ⟨σv⟩_pp = 1.597×10⁻⁴³ cm³/s from the 1.301×10⁻¹⁴ form, or 1.593×10⁻⁴³ cm³/s from the exact prefactor. With n_p = 3.07×10²⁵ cm⁻³ (assumed ρ = 150 g/cm³, X = 0.34), the proton lifetime against pp burning is ≈6.5×10⁹ yr.
- At 10 keV the D–T barrier is 30–45× the plasma temperature, and the classical turning point (144 fm) is ~40× the nuclear radius. Ignition-relevant reactions therefore come from the Gamow-peak tail at 3–4 kT.

### Gaps
- I could not access Atzeni & Meyer-ter-Vehn, Clayton or Iliadis directly, so the 0.1220/0.2368 constants come from search summaries of Rauscher (2010) and a 2025/26 arXiv paper. There is no authoritative single "published" barrier number because the radius convention differs between sources. The only explicit published barrier test value found is the 360 keV textbook problem.
- No published Gamow-peak width for the Sun was found in a primary source fetched here; the 5.9/6.4 keV pair is folklore that I could only recompute.

## Q2. S-factors and reactivities (Bosch–Hale 1992; p–¹¹B), Lawson criterion and triple product, NIF ignition and records

### Takeaway
The Bosch–Hale fits reproduce published tables to ≤0.3% (computed). D–T ⟨σv⟩ at 10 keV is 1.136×10⁻¹⁶ cm³/s (computed; published table 1.14×10⁻¹⁶). With f_c = 0.2 and C_B = 5.34×10⁻³⁷ W m³ keV^{−½}, the Wurzel–Hsu ignition triple product is reproduced: minimum n T τ_E ≈ 2.97×10²¹ keV s m⁻³ at 14.25 keV (published: 2.9×10²¹ at 14 keV), and the ideal ignition temperature is 4.30 keV (published: 4.3 keV). For p–¹¹B, Tentori & Belloni (2023) raise the Nevins–Swain reactivity by 12% at 100 keV and 50% at 500 keV, and both are reproduced (computed). NIF: 2.05 MJ in gave 3.15 MJ out on 5 Dec 2022, and the record is 8.6 ± 0.45 MJ from 2.08 MJ (gain 4.13) on 7 Apr 2025. By June 2026 there had been 11 ignition shots.

### Cited Findings
- [mainstream] Bosch–Hale reactivity parameterisation, with T in keV and the result in cm³/s:
  ```
  θ  = T / ( 1 − T(C2 + T(C4 + T C6)) / (1 + T(C3 + T(C5 + T C7))) )
  ξ  = ( B_G² / (4θ) )^(1/3)
  <σv> = C1 · θ · sqrt( ξ / (m_r c² · T³) ) · exp(−3ξ)
  ```
  — [UWFDM-1268 Eqs. 2–5](https://fti.neep.wisc.edu/fti.neep.wisc.edu/pdf/fdm1268.pdf); [MFE Formulary notebook](https://notebook.community/MFEFormulary/MFEFormulary/ipython/BH_ReactionRates). Original: H.-S. Bosch & G.M. Hale, Nucl. Fusion 32, 611–631 (1992) — [OSTI record](https://www.osti.gov/etdeweb/biblio/5161054). An erratum exists — [Semantic Scholar](https://www.semanticscholar.org/paper/ERRATUM:-Improved-formulas-for-fusion-and-thermal-Bosch-Hale/ba3cd4b3e38e46e7157d67f932039c16a540def6).
- [mainstream] Coefficients (exact digits from the MFE notebook):

  | reaction | C1 | C2 | C3 | C4 | C5 | C6 | C7 |
  |---|---|---|---|---|---|---|---|
  | ³H(d,n)⁴He | 1.17302e-9 | 1.51361e-2 | 7.51886e-2 | 4.60643e-3 | 1.35000e-2 | −1.06750e-4 | 1.36600e-5 |
  | ²H(d,p)³H | 5.65718e-12 | 3.41267e-3 | 1.99167e-3 | 0 | 1.05060e-5 | 0 | 0 |
  | ²H(d,n)³He | 5.43360e-12 | 5.85778e-3 | 7.68222e-3 | 0 | −2.96400e-6 | 0 | 0 |
  | ³He(d,p)⁴He | 5.51036e-10 | 6.41918e-3 | −2.02896e-3 | −1.91080e-5 | 1.35776e-4 | 0 | 0 |

  — [MFE Formulary notebook](https://notebook.community/MFEFormulary/MFEFormulary/ipython/BH_ReactionRates). The same values at 3 s.f. are in [UWFDM-1268 Table 5](https://fti.neep.wisc.edu/fti.neep.wisc.edu/pdf/fdm1268.pdf).
- [mainstream] Reactivity validity ranges: 0.2–100 keV for T(d,n), D(d,n) and D(d,p), and 0.5–190 keV for ³He(d,p). — [UWFDM-1268 Table 7](https://fti.neep.wisc.edu/fti.neep.wisc.edu/pdf/fdm1268.pdf). The MFE notebook's "T range" column (0.5–550, 0.5–4900, 0.5–5000, 0.3–900 keV) is actually the cross-section (S-fit) energy range.
- [mainstream] Published check values of the Bosch–Hale model in cm³/s ([UWFDM-1268 Table 6](https://fti.neep.wisc.edu/fti.neep.wisc.edu/pdf/fdm1268.pdf)):

  | T (keV) | T(d,n)⁴He | D(d,n)³He | D(d,p)T | DD total | ³He(d,p)⁴He |
  |---|---|---|---|---|---|
  | 1 | 6.86e-21 | 9.93e-23 | 1.02e-22 | 2.01e-22 | 3.05e-26 |
  | 2 | 2.98e-19 | 3.11e-21 | 3.15e-21 | 6.26e-21 | 1.40e-23 |
  | 5 | 1.37e-17 | 9.13e-20 | 9.02e-20 | 1.82e-19 | 6.36e-21 |
  | 10 | 1.14e-16 | 6.02e-19 | 5.78e-19 | 1.18e-18 | 2.12e-19 |
  | 20 | 4.33e-16 | 2.60e-18 | 2.40e-18 | 5.00e-18 | 3.48e-18 |
  | 50 | 8.65e-16 | 1.13e-17 | 9.84e-18 | 2.12e-17 | 5.55e-17 |
  | 100 | 8.45e-16 | 2.68e-17 | 2.24e-17 | 4.93e-17 | 1.72e-16 |
- [mainstream] Older NRL Plasma Formulary values (based on Duane) in cm³/s, at T = 1, 2, 5, 10, 20, 50, 100 keV. DT: 5.5e-21, 2.6e-19, 1.3e-17, 1.1e-16, 4.2e-16, 8.7e-16, 8.5e-16. DD total: 1.5e-22, 5.4e-21, 1.8e-19, 1.2e-18, 5.2e-18, 2.1e-17, 4.5e-17. D³He: 1.0e-26, 1.4e-23, 6.7e-21, 2.3e-19, 3.8e-18, 5.4e-17, 1.6e-16. — [UWFDM-1268 Table 1, reproducing the NRL Formulary](https://fti.neep.wisc.edu/fti.neep.wisc.edu/pdf/fdm1268.pdf). (The NRL PDF itself returned HTTP 403.)
- [mainstream] Machine-precision check values from a Bosch–Hale implementation: D(d,n)³He at T = 0.1 keV gives 1.65721857e-37 m³/s (below the fit's validity range), and at T = 100 keV gives 2.68167420e-23 m³/s. — [MFE Formulary notebook printed output](https://notebook.community/MFEFormulary/MFEFormulary/ipython/BH_ReactionRates)
- [mainstream] Bosch–Hale cross-section S-factor fit (E in keV CM, S in keV·mb, σ in mb):
  ```
  S(E) = (A1 + E(A2 + E(A3 + E(A4 + E A5)))) / (1 + E(B1 + E(B2 + E(B3 + E B4))))
  σ(E) = S(E) / (E · exp(B_G / sqrt(E)))
  ```
  Coefficients, identified by B_G and by energy range, since the column headers on the extracted page are misaligned:

  | | D(d,p)T | D(d,n)³He | D–³He | D–T (low range) |
  |---|---|---|---|---|
  | B_G | 31.3970 | 31.3970 | 68.7508 | 34.3827 |
  | A1 | 5.5576e4 | 5.3701e4 | 5.7501e6 | 6.927e4 |
  | A2 | 2.1054e2 | 3.3027e2 | 2.5226e3 | 7.454e8 |
  | A3 | −3.2638e-2 | −1.2706e-1 | 4.5566e1 | 2.050e6 |
  | A4 | 1.4987e-6 | 2.9327e-5 | 0 | 5.2002e4 |
  | A5 | 1.8181e-10 | −2.5151e-9 | 0 | 0 |
  | B1 | 0 | 0 | −3.1995e-3 | 6.38e1 |
  | B2 | 0 | 0 | −8.5530e-6 | −9.95e-1 |
  | B3 | 0 | 0 | 5.9014e-8 | 6.981e-5 |
  | B4 | 0 | 0 | 0 | 1.728e-4 |
  | range (keV) | 0.5–5000 | 0.5–4900 | 0.3–900 | 0.5–550 |
  | (ΔS)max | 2.0% | 2.5% | 2.2% | 1.9% |

  — [Wikipedia: Nuclear fusion (Bosch–Hale table)](https://en.wikipedia.org/wiki/Nuclear_fusion). The D–T column matches [LLNL pyPICfusion `_BoschHale_library.py`](https://github.com/llnl/pyPICfusion/blob/master/_BoschHale_library.py) (via search summary), and the ranges match the [MFE notebook](https://notebook.community/MFEFormulary/MFEFormulary/ipython/BH_ReactionRates).
- [mainstream] Published cross-section check value: "the cross section for D-D fusion is effectively zero below 20 keV, and reaches 16.3 mb at 100 keV" (deuteron lab energy). — [Geuther, Danon & Saglime, PRL 96, 054803 (2006)](https://pages.github.rpi.edu/danony/html/Papers/PhysRevLett_96_054803%20Nuclear%20Reactions%20Induced%20by%20a%20Pyroelectric%20Accelerator.pdf)
- [mainstream] p–¹¹B: [Tentori & Belloni, Nucl. Fusion 63, 086001 (2023)](https://iopscience.iop.org/article/10.1088/1741-4326/acda4b).
  - High-temperature fit (70–500 keV, error <1%): a Bosch–Hale-type form with ξ = (E_G/4θ)^{1/3} and fit parameters P₁…P₇:

    | | P₁ [keV m³ s⁻¹] | P₂ [keV⁻¹] | P₃ [keV⁻¹] | P₄ [keV⁻²] | P₅ [keV⁻²] | P₆ [keV⁻³] | P₇ [keV⁻³] |
    |---|---|---|---|---|---|---|---|
    | Nevins & Swain 2000 | 4.4467e-14 | −5.9357e-2 | 2.0165e-1 | 1.0404e-3 | 2.7621e-3 | −9.1653e-6 | 9.8305e-7 |
    | Tentori & Belloni | 9.7827e-14 | −5.1610e-2 | 1.3240e-1 | 3.8446e-4 | 1.2499e-3 | −5.6715e-6 | 1.1615e-6 |

  - Low-temperature fit (10–70 keV, error <2%): a non-resonant term plus the 148 keV resonance term.
  - Their reactivity is higher than Nevins–Swain above 30 keV: "12% at 100 keV and 50% at 500 keV", more than 20% at 200 keV, and 30% at 300 keV.
  - New S-factor fit for E ≤ 0.400 MeV: C₀ = 197, C₁ = 0.269 and C₂ = 2.54×10⁻⁴ MeV b.
- [mainstream] Lawson / triple product, with n in m⁻³, T in keV, ε_F the fusion energy in J and f_c the charged fraction:
  ```
  P_B = C_B n_i n_e Z² T_e^(1/2) V,    C_B = 5.34e-37 W m³ keV^(-1/2)
  ignition (steady state, no conduction):  (1/4) f_c n² <σv>_DT ε_F ≥ C_B n² T^(1/2)
  n T τ_E = 3T² / [ (f_c + 1/Q_fuel) <σv> ε_F / 4 − C_B T^(1/2) ]      (Eq. 22)
  ```
  - "When T < 4.3 keV, S_B > S_c and ignition is not possible."
  - Table II, minimum n_i τ_E for D+T: 2.5×10¹⁹ m⁻³ s at 26 keV for Q_fuel = 1, and 1.6×10²⁰ m⁻³ s at 26 keV for Q_fuel = ∞.
  - Table III, minimum n_i T τ_E for D+T: 4.6×10²⁰ keV s m⁻³ at 14 keV for Q_fuel = 1, and 2.9×10²¹ at 14 keV for Q_fuel = ∞.
  - Catalysed D–D: 1.1×10²³ at 52 keV. D+³He: 5.2×10²² at 68 keV. p+¹¹B: no entry ("–").
  - Reactivities are integrated from Bosch–Hale cross sections. Q_sci = P_F/P_ext. The 1997 National Academies definition of ICF ignition is fusion energy exceeding laser energy (Q_sci > 1).
  — [Wurzel & Hsu, Phys. Plasmas 29, 062103 (2022)](https://arxiv.org/abs/2105.10954)
- [mainstream] NIF shot record as listed by LLNL:

  | date | laser energy on target | fusion yield | target gain |
  |---|---|---|---|
  | 5 Dec 2022 | 2.05 MJ | 3.15 MJ | — |
  | 30 Jul 2023 | 2.05 MJ | 3.88 MJ | — |
  | 8 Oct 2023 | 1.9 MJ | 2.4 MJ | — |
  | 30 Oct 2023 | 2.2 MJ | 3.4 MJ | — |
  | 12 Feb 2024 | 2.2 MJ | est. 5.2 MJ | — |
  | 18 Nov 2024 (6th) | 2.2 MJ | 4.1 MJ | — |
  | 23 Feb 2025 (7th) | 2.05 MJ | 5.0 MJ | 2.44 |
  | 7 Apr 2025 (8th, record) | 2.08 MJ | 8.6 ± 0.45 MJ | 4.13 |
  | 22 Jun 2025 (LANL-led) | — | 2.4 ± 0.09 MJ | — |
  | 1 Oct 2025 (10th) | 2.065 MJ | 3.5 ± 0.17 MJ | 1.74 |
  | 20 Jun 2026 (11th) | — | 7.9 ± 0.4 MJ | ≈3.8 |

  — [LLNL/LMF "Achieving Fusion Ignition"](https://lmf.llnl.gov/science/achieving-fusion-ignition)
  - LLNL's FY2025 annual report gives the record shot's peak power as 456 TW and its gain as 4.13. It reports the 1 Oct 2025 yield as 3.6 MJ, which conflicts with the 3.5 MJ above. — [LLNL Annual FY2025](https://annual.llnl.gov/fy-2025/national-ignition-facility-2025)
  - LLNL experts estimate that a target gain near 15 is needed for a test facility to reach zero net electric demand. — [LLNL annual report PDF](https://annual.llnl.gov/sites/annual/files/2026-02/llnlannual25-3-national-ignition-facility.pdf) (via search summary)

### Inferences
- Recomputed Bosch–Hale values (computed, cm³/s) agree with UWFDM-1268 Table 6 within the 3-s.f. rounding (max |Δ| = 0.3%), e.g.:
  - T(d,n) at 10 keV: 1.1362e-16; at 1 keV: 6.8569e-21.
  - D(d,n) at 10 keV: 6.0227e-19. D(d,p) at 10 keV: 5.7813e-19. ³He(d,p) at 10 keV: 2.1261e-19.
  - The MFE endpoints are reproduced exactly: 1.6572186e-37 and 2.6816742e-23 m³/s.
  - Inside the fit range (≤100 keV), D–T ⟨σv⟩ peaks at T = 66.7 keV with 8.9465e-16 cm³/s (computed). Wikipedia's "about 70 keV" agrees.
  - NRL (Duane-based) values differ from Bosch–Hale by −20% for D–T at 1 keV, −67% for D–³He at 1 keV and +8% for D–³He at 10 keV. Use Bosch–Hale as the reference.
- S-fit checks (computed):
  - D–T σ peaks at E_cm = 64.8 keV (deuteron lab ≈ 108 keV) with σ = 5.066 b, matching the "peak near 110 keV" lecture value.
  - σ_D(d,n) at E_cm = 50 keV (lab 100 keV) is 16.49 mb and σ_D(d,p) is 15.57 mb. Geuther's "16.3 mb at 100 keV" therefore matches the n-branch within 1%.
  - D–T at 10 keV CM: S = 1.425×10⁴ keV·b and σ = 27.0 mb.
- Independent consistency test for a Rust crate (computed): a numerical Maxwellian average of the S-fit cross sections, ⟨σv⟩ = √(8/πμ)(kT)^{−3/2}∫σE e^{−E/kT}dE, matches the reactivity fit to 0.2–1.2% for D–T and D–D at 2, 10 and 50 keV. For D–³He it matches to 0.6–3.1%, with the largest gap at 10 keV.
- Lawson reproduction (computed) with f_c = 3.5/17.6, ε_F = 17.6 MeV, C_B as above and BH reactivity:
  - Q_fuel = ∞: min nTτ = 2.968×10²¹ keV s m⁻³ at 14.25 keV, and min nτ = 1.567×10²⁰ s m⁻³ at 26.0 keV.
  - Q_fuel = 1: min nTτ = 4.677×10²⁰ at 13.65 keV, and min nτ = 2.527×10¹⁹ at 25.7 keV.
  - Ideal ignition temperature: 4.305 keV.
  - All match Wurzel–Hsu Tables II–III and the 4.3 keV statement to the published 2 s.f., which makes them good end-to-end integration tests.
- p–¹¹B reproduction (computed). The formula is the Bosch–Hale form with E_G = 22589.1 keV (B_G² from CODATA) and m_r c² = 859604.5 keV. For Nevins–Swain I added the commonly used resonance term 5.41×10⁻²¹ T^{−3/2} exp(−148/T) m³/s; that term is my assumption, from memory of Nevins & Swain 2000.
  - TB/(NS + res) − 1 = +12.6% at 100 keV, +21.4% at 200 keV, +29.4% at 300 keV and +50.1% at 500 keV. This reproduces all four published statements, so both the form and the coefficients are confirmed.
  - Values (m³/s): TB(100 keV) = 7.07e-23, TB(300 keV) = 4.39e-22, TB(500 keV) = 5.62e-22. NS(+res) peaks at ≈3.75e-22 near 476 keV, while TB is still rising at 500 keV.
  - Without the resonance term, the ratio at 100 keV would be +14.8%.
- NIF gains (computed) are 1.54 (Dec 2022), 1.89 (Jul 2023), 2.44 (Feb 2025) and 4.13 (Apr 2025). "Target gain" counts only laser energy delivered to the target, not wall-plug energy.
- The June 2026 shot's implied laser energy is 7.9/3.8 ≈ 2.08 MJ (computed).

### Gaps
- Bosch & Hale (1992) itself is paywalled (IOP), so Table VIII was not seen directly. The 3-s.f. UWFDM table and the MFE notebook are secondary reproductions. Also not obtained: the erratum's content, the D–T high-energy S-fit (550–4700 keV), and the stated reactivity-fit accuracy.
- The NRL Formulary PDF returned 403. Its table is used via the UWFDM reproduction.
- The Tentori–Belloni equations are images on the IOP page, so the exact LT formula (NR plus resonance terms) and E_G/m_r c² values were not transcribed. The Nevins–Swain resonance-term coefficient is from memory, but its use is validated by the reproduction above.
- LLNL's "about 300 MJ wall-plug per shot" figure was not found on the pages fetched.

## Q3. Electron screening: U_e, enhancement f = exp(πηU_e/E), Salpeter weak screening, measured gas-target values

### Takeaway
Bound electrons lower the barrier by a constant energy U_e. The lab enhancement is f = exp(πη(E)U_e/E) [mainstream]. For d+d with gas targets, U_e ≈ 20–25 eV, consistent with the adiabatic limit, gives only ~1% at 10 keV but ~40% at 1 keV (computed). Gas-target values for some reactions exceed the adiabatic limit by ~2×: ³He(d,p) has U_e = 219 ± 15 eV vs U_ad = 119 eV [contested; unexplained]. Metal hosts give "several 100 eV", with values varying between labs [contested]. In stars, Salpeter weak screening gives f₀ ≈ 1.05 for solar p+p (computed).

### Cited Findings
- [mainstream] Lab enhancement: f_lab(E) ≡ σ_lab(E)/σ_b(E) ∼ exp(πη(E)U_e/E) ≥ 1, for U_e ≪ E (Assenbaum, Langanke & Rolfs 1987). Here U_e is estimated from the difference in atomic binding energies between the compound atom and the entrance-channel atoms. — [Solar Fusion II, Eq. 15](https://arxiv.org/abs/1004.2318)
- [mainstream] The equivalent form S_lab(E) = S_b(E) exp(πη(E)U_e/E) comes from expanding the Gamow exponent to first order after E → E + U_e. It neglects a small normalisation correction. — [Solar Fusion III, Eq. 25](https://arxiv.org/abs/2405.06470)
- [mainstream] Salpeter (1954) weak screening:
  - Debye potential V(r) = (αZ₁Z₂/r) exp(−r/R_D).
  - Enhancement f₀ ∼ exp(Z₁Z₂α/(R_D kT)) = exp(0.188 Z₁Z₂ ζ ρ₀^{1/2} T₆^{−3/2}), where ζ = [Σᵢ XᵢZᵢ²/Aᵢ + (f′/f)Σᵢ XᵢZᵢ/Aᵢ]^{1/2}, f′/f ≈ 0.92 for electron degeneracy in the solar core, ρ₀ is in g/cm³ and T₆ in 10⁶ K.
  - Validity requires Z₁Z₂α/(R_D kT) ≪ 1, satisfied in the solar core for Z₁Z₂ ≲ 10.
  - Field-theory analysis expects only ~4% corrections to the standard formula for solar conditions, though "controversies have not entirely died out".
  — [Solar Fusion II, Eqs. 13–14](https://arxiv.org/abs/1004.2318); [Solar Fusion III, Eqs. 27–28](https://arxiv.org/abs/2405.06470)
- [mainstream] d+d with a gaseous target (Greife et al. 1995) "were compatible with the adiabatic limit (U_e = 20 eV, Bracci et al., 1991, 1990)". — [Solar Fusion III §XII.A](https://arxiv.org/abs/2405.06470). Conference slides quote the Greife gas-target value as U_e = 25 ± 5 eV (Z. Phys. A351, 1995). — [Czerski slides (COMEX5)](https://comex5.ifj.edu.pl/slides/czerski.pdf) (via search summary)
- [contested] ³He(d,p)⁴He: data represented by U_e = 219 ± 15 eV (Aliotta et al. 2001), "significantly larger" than the adiabatic U_ad = 119 eV. Extraction depends on the assumed bare σ_b(E). — [Solar Fusion II](https://arxiv.org/abs/1004.2318)
- [contested] ³He(³He,2p)⁴He: the fit gave U_e = 305 ± 90 eV, and screening raises S_lab by ~40% at the lowest data point. — [Solar Fusion III](https://arxiv.org/abs/2405.06470)
- [contested] Experimental U_e tend to exceed theory by a factor ~2. Proposed explanations include low-energy stopping-power errors (smaller stopping powers were verified experimentally), Trojan-Horse σ_b constraints that lower U_e (Tumino et al.), and clusterisation or polarisation. — [Solar Fusion II](https://arxiv.org/abs/1004.2318); [Solar Fusion III](https://arxiv.org/abs/2405.06470)
- [contested] For d+d in metals, screening energies "varied among themselves and could reach values of several 100 eV". The gas adiabatic limit does not apply to implanted targets. Huke et al. (2008) identified O and C contamination as the source of some very large values. Metal experiments are said to confirm the Debye-model trend U_e(T) ∝ T^{−1/2}. — [Solar Fusion III](https://arxiv.org/abs/2405.06470); [Solar Fusion II](https://arxiv.org/abs/1004.2318)
- [mainstream] The simple atomic estimate for D(d,p)t is U_ad ≈ 27.2 eV (2 Ry). Other model values are 39 eV for neutral atoms and 52 eV for an ionised projectile. — search summary of [ar5iv 0804.3886](https://ar5iv.arxiv.org/html/0804.3886) / [ar5iv 0805.4538](https://ar5iv.arxiv.org/html/0805.4538) (secondary; not verified in full text)

### Inferences
- Enhancement table (computed). f is the first-order formula; f_exact = P(E+U_e)/P(E) applies the shift inside the Gamow factor only. The two differ significantly once U_e/E ≳ 0.05, so a calculator should expose both.

  | reaction, U_e | E = 1 keV | 2 keV | 5 keV | 10 keV | 20 keV |
  |---|---|---|---|---|---|
  | d+d, 20 eV | 1.369 (1.363) | 1.117 | 1.0285 | 1.0100 | 1.0035 |
  | d+d, 25 eV | 1.481 (1.470) | 1.149 | 1.0357 | 1.0125 | 1.0044 |
  | d+d, 300 eV (metal-like) | 111 (47.5) | 5.29 (4.47) | 1.524 (1.497) | 1.161 | 1.054 |
  | d+³He, 119 eV | 59.8 (42.9) | 4.25 (3.99) | 1.442 | 1.138 | 1.047 |
  | d+³He, 219 eV | 1862 (653) | 14.3 (11.7) | 1.961 (1.920) | 1.269 | 1.088 |
- Salpeter f₀ for solar p+p (computed, assumptions stated): ρ = 150 g/cm³, X = 0.34, Y = 0.64, Z = 0.02 with an assumed Z²/A ≈ 3.5 for metals, and f′/f = 0.92. These give ζ = 1.291 and f₀ = 1.049 at T₆ = 15.7, or 1.050 at T₆ = 15.5. Screening matters at the 5% level in the Sun.
- In a 10–15 keV fusion plasma, screening corrections to ⟨σv⟩ are negligible. The Gamow peak sits at 30–40 keV, so even U_e = 300 eV changes σ by only a few percent there (computed). Screening is a big lever only for cold-target, low-E beam physics.

### Gaps
- No primary-text value of Greife (1995) was read. The 25 ± 5 eV figure comes only from conference slides (search summary), while SF-III says the data are "compatible" with 20 eV. Treat the gas d+d value as 20–25 eV.
- No published Salpeter f₀ number for solar p+p was found in the fetched text, so my 1.049 depends on assumed composition and the metal Z²/A.

## Q4. Muon-catalysed fusion: cycle rate, α-sticking, fusions per muon, energy breakeven arithmetic

### Takeaway
The best measured yield is ≈150 d–t fusions per muon (LAMPF, Jones et al. 1986). That is ≈2.6 GeV of fusion energy against an estimated ≈4.7–6 GeV cost to make each muon, while scientific breakeven needs ≈270–340 fusions per muon (computed) [mainstream]. The ceiling is set by effective α-sticking ω_s ≈ 0.4–0.6% (initial ω_s⁰ = 0.912%) and the cycle time relative to the 2.197 μs muon lifetime: X_μ = 1/(ω_s + λ₀/λ_c).

### Cited Findings
- [mainstream] Up to 150 ± 20 fusions per muon were observed (Jones et al., PRL 1986) in high-density cold D–T. — [MuFusE paper, arXiv:2606.05333 (2026)](https://arxiv.org/abs/2606.05333)
- [mainstream] Nature review abstract: catalysis yields of ~150 fusions per muon. — [Jones, Nature 321, 127 (1986)](https://www.nature.com/articles/321127a0) (abstract via search summary)
- [mainstream] The yield measured experimentally at ~1000 °C and liquid-hydrogen density (LHD ≡ 4.25×10²² atoms/cm³) is "less than about 170 per muon".
  - Muon production costs "about 5 GeV", so about 285 fusions would be needed for "scientific break-even (1/3 of the commercial break-even)".
  - Initial sticking ω_s⁰ = 0.912% (Hu, Hale & Cohen 1994), with μα⁺ formed at 3.47 MeV. τ_μ = 2.197×10⁻⁶ s.
  - Effective sticking is ω_s^eff = ω_s⁰(1 − R), where R is the reactivation (stripping) fraction.
  - Their computed ω_s^eff falls from 0.580% at 0.1 eV to 0.347% at 1000 eV plasma temperature [theory].
  — [Pahlavani & Motevalli, Acta Phys. Pol. B 40, 319 (2009)](https://www.actaphys.uj.edu.pl/R/40/2/319/pdf)
- [mainstream] Cost per muon "about 5–6 GeV" in the literature (Jändel).
  - Their own estimate is 4.70 GeV per muon: a 3.61 GeV deuteron beam on tungsten giving 0.77 π⁻/μ⁻ per beam particle (simulation).
  - Effective sticking p ≅ 0.47% and t_c/τ ≅ 0.18% at solid density give X = 1/(p + t_c/τ) ≈ 154, consistent with 150.
  - [contested: single-group design study] Q = 1.87 for an optimised configuration, and Q_elec = 14% at 150 fusions per muon, assuming 18% accelerator efficiency (PSI), 50% muon delivery and >60% heat-to-electric conversion.
  — [Kelly, Hart & Rose, J. Phys. Energy 3, 035003 (2021)](https://iopscience.iop.org/article/10.1088/2515-7655/abfb4b)
- [mainstream] Jackson (1957) estimated ~1% sticking, implying ~100 fusions per muon. Later measurements give ~0.3–0.5%, implying up to ~200 (even 350). One commonly quoted electrical cost is ≈6 GeV per muon, uncited on that page. — [Wikipedia: Muon-catalyzed fusion](https://en.wikipedia.org/wiki/Muon-catalyzed_fusion)
- [mainstream] Measured high-density ω_s is "typically 10–50% lower than those predicted by standard theoretical models". The cycling rate λ_c generally increases with temperature, and the highest rates are in cryogenic solid D–T. — [MuFusE, arXiv:2606.05333](https://arxiv.org/abs/2606.05333)
- [mainstream] Fusion inside the muonic molecule proceeds at rates "of order 10¹² s⁻¹". — [Koonin & Nauenberg preprint (1989)](https://newenergytimes.com/v2/library/1989/1989KooninNauenberg.shtml)

### Inferences
- Calculator (computed):
  - λ₀ = 1/τ_μ = 4.5517×10⁵ s⁻¹.
  - Kelly's t_c/τ = 0.18% implies λ_c = λ₀/0.0018 ≈ 2.5×10⁸ s⁻¹.
  - Then X_μ = 1/(0.0047 + 0.0018) = 153.8.
  - Sticking-only caps 1/ω_s: 110 at 0.912%, 213 at 0.47%, 250 at 0.40%, 333 at 0.30%.
- Breakeven arithmetic (computed): fusions per muon = E_μ/17.6 MeV.
  - Energy breakeven: 267 at 4.70 GeV, 284 at 5 GeV, 341 at 6 GeV. Pahlavani gives "about 285" at 5 GeV.
  - Including an electricity-to-heat round trip at η = 0.4 (0.33): 668 (809), 710 (861) and 852 (1033) respectively.
  - 150 × 17.6 MeV = 2.64 GeV ≈ 25 m_μc². Muon-catalysed fusion therefore beats the muon rest-mass energy but not its production cost.
- Even at the theoretical sticking floor, μCF stays below 1/ω_s ≈ 200–330 unless reactivation improves or muon cost falls by ~2×. That is why it remains [mainstream] science without a power case.

### Gaps
- No primary measured value of λ_c (s⁻¹) at specific density and temperature was read; the ≈2.5×10⁸ s⁻¹ is derived from Kelly et al.'s ratio. The Breunlich et al. (1989) and Petitjean (1992) reviews and the PSI/RIKEN-RAL sticking measurements were not fetched.
- The 5 GeV per muon estimate's provenance ("1980 estimate", Jändel) was not traced to its primary source.

## Q5. Pyroelectric fusion, lattice confinement fusion (NASA Glenn), fusors

### Takeaway
All three devices produce genuine D–D beam-target or near-beam-target neutrons at ~10²–10⁸ fusions/s. Their fusion power is ~10⁻¹⁴ to ~10⁻⁶ of the beam or electrical input (computed; ~10⁻⁶ is the pyroelectric beam alone, ignoring heater energy). Pyroelectric fusion is [mainstream], having been replicated independently. NASA's bremsstrahlung "lattice confinement" neutrons are [contested]: one group, with enhanced-screening and Oppenheimer–Phillips interpretations. Fusors are [mainstream], and losses on the grid and to collisions keep them 5+ orders of magnitude below breakeven.

### Cited Findings
- [mainstream] Naranjo, Gimzewski & Putterman, "Observation of nuclear fusion driven by a pyroelectric crystal", Nature 434, 1115–1117 (2005). Abstract: a deuteron beam "(> 100 keV and >4 nA)" striking a deuterated target produces "a neutron flux over 400 times the background level". — [RePEc abstract page](https://ideas.repec.org/a/nat/nature/v434y2005i7037d10.1038_nature03575.html) (numbers seen in the search snippet; the fetched page text was qualitative only)
- [mainstream] Wikipedia's details of the 2005 experiment:
  - LiTaO₃ crystal "heated from −34 to 7 °C" over a few minutes.
  - A tungsten needle gives "about 25 gigavolts per meter".
  - About 100 keV ion energy (a modelling estimate).
  - "Around 1000 fusion reactions per second" on an ErD₂ target.
  - The apparatus requires far more energy than it produces.
  — [Wikipedia: Pyroelectric fusion](https://en.wikipedia.org/wiki/Pyroelectric_fusion)
- [mainstream] RPI replication used paired LiTaO₃ crystals (20 mm diameter × 10 mm), a tungsten tip and a deuterated polystyrene target, with no cryogenics.
  - X-ray endpoint 200 keV, i.e. a 200 kV potential.
  - Net 138 ± 7 neutron counts per cycle over an average 120 s emission, against background 0.14 ± 0.01 cps.
  - Detector intrinsic efficiency ≈0.4 and solid-angle factor ≈0.41.
  - Source strength ≈10,500 neutrons per cooling phase, with an inferred ion current of 0.55 nA. Theory gives ≈18,950 neutrons per heating cycle per nA.
  - Naranjo et al. reported 15,300 neutron counts per cycle, i.e. 85,000 source neutrons at detection efficiency 0.18.
  — [Geuther, Danon & Saglime, PRL 96, 054803 (2006)](https://pages.github.rpi.edu/danony/html/Papers/PhysRevLett_96_054803%20Nuclear%20Reactions%20Induced%20by%20a%20Pyroelectric%20Accelerator.pdf)
- [mainstream] Tornow, Shafroth & Brownridge reported neutrons without a tip, with ions up to 200 keV at 3 nA. — [APS April 2007 abstract](https://meetings-archive.aps.org/apr/2007/k15/4) (via search summary)
- [contested] NASA Glenn bremsstrahlung experiment ("lattice confinement fusion"): [Steinetz et al., Phys. Rev. C 101, 044610 (2020); NASA/TP-20205001616](https://ntrs.nasa.gov/api/citations/20205001616/downloads/TP-20205001616.pdf). The theory companion is [Pines et al., Phys. Rev. C 101, 044609 (2020), "Nuclear Fusion Reactions in Deuterated Metals"](https://journals.aps.org/prc/abstract/10.1103/PhysRevC.101.044609) (citation from the TP's reference 12).
  - Method: a 2.9 MeV, 15 mA DC electron beam (range 2.5–3.0 MeV, 5–30 mA) on a bremsstrahlung converter irradiates ErD₃ and TiD₂.
  - Photons above 2.226 MeV photodissociate D, giving photoneutrons below 400 keV (average 145 keV). These elastically heat deuterons (average 64 keV), which fuse with cold lattice D.
  - Fuel: 5×10²⁴ D atoms in each case, as 480 g ErD₃ in 16 vials or 216 g TiD₂. Number densities: ErD₃ 0.8×10²³ D/cm³, TiD₂ 1×10²³ D/cm³.
  - Case-control: unfueled Er and Ti exposures under the same conditions.
  - Measured: 2.45 MeV neutrons at 1.5 ± 0.3×10³ n/s (TS1575) and 1.6 ± 0.3×10³ n/s (TS1576) with EJ-309, scaled to the sample location. The prediction using shell and plasma screening (λ_sc = 4.16×10⁻¹⁰ cm) is 1.2×10³ n/s.
  - Detection: absolute efficiency ≈7 counts (EJ-309) or ≈2 (stilbene) per 10⁶ fusion neutrons. Events above an 8σ pulse-shape cut were accepted.
  - TiD₂ produced 1.31× the neutrons of ErD₃.
  - They also report ~4 and ~5 MeV neutrons, attributed to "boosted" (heated-deuteron) fusion or Oppenheimer–Phillips stripping.
- [mainstream] A student-built Farnsworth–Hirsch fusor (TU Eindhoven) produced "just below 10⁶ neutrons per second" at −43 kV and 30 mA in deuterium. — [FuseNet](https://fusenet.eu/node/361)
- [mainstream] Hirsch's late-1960s device (six ion sources, up to 150 kV) reached nearly 10⁸ fusions/s in deuterium, under 0.1 W even for D–T, against kilowatt-scale consumption. — [Chemie.de Lexikon (German Wikipedia mirror)](https://www.chemie.de/lexikon/Farnsworth-Hirsch-Fusor.html). Original: Hirsch, J. Appl. Phys. 38, 4522 (1967).
- [mainstream] Breakeven arguments:
  - "no fusor has ever come close to break-even energy output".
  - Grid losses "can be at least five orders of magnitude higher than the energy released".
  - Rider (MIT 1995; Phys. Plasmas 4, 1039, 1997) showed that non-equilibrium (non-Maxwellian) plasmas lose energy to thermalisation and bremsstrahlung faster than fusion can supply it.
  - Highest "fusor-like" flux: 3×10¹¹ D–D n/s (Phoenix Nuclear Labs, 2013 press release). That device is a beam-target generator, not a gridded fusor.
  — [Wikipedia: Fusor](https://en.wikipedia.org/wiki/Fusor)
- Hobbyist reports put simple d–d fusors rarely above ~10⁵ n/s. — [fusor.net forum](https://www.fusor.net/old-boards/songs.com/msgs/msg-907.html) (low-quality source)

### Inferences
- Energy balances (computed; D–D mean Q ≈ 3.65 MeV per fusion, and ≈2 fusions per neutron because the branches are ~50/50):
  - **Pyroelectric**: beam ≥ 100 keV × 4 nA = 0.4 mW. ~1000 fusions/s ≈ 5.8×10⁻¹⁰ W, so P_fus/P_beam ≈ 1.5×10⁻⁶. 85,000 neutrons per cycle ≈ 1.0×10⁻⁷ J of fusion energy per cycle. RPI's 10,500 neutrons over ~120 s is ≈90 n/s.
  - **NASA LCF**: beam 2.9 MV × 15 mA = 43.5 kW. 1.5×10³ n/s gives ≈1.8×10⁻⁹ W of fusion power, so P_fus/P_beam ≈ 4×10⁻¹⁴. The theory calculation used 16 × 450 µA = 7.2 mA, which still gives ~10⁻¹³.
  - **Fusor (Eindhoven)**: 1.29 kW in, ≈1.2×10⁻⁶ W out, so Q ≈ 9×10⁻¹⁰. Hirsch's ~10⁸ D–D fusions/s is ≈6×10⁻⁵ W, or ≈2.8×10⁻⁴ W at the same rate in D–T.
- Physically, all three are beam-target fusion. The fast ion (pyro and fusor) or the photoneutron-heated deuteron (LCF) carries 10–200 keV, and the thick-target yield is set by ∫σ(E)/(dE/dx) dE. The Bosch–Hale σ above therefore suffices to estimate them; the RPI paper does exactly this and matches its yield to first order.
- LCF: the 2.45 MeV signal is plausible from conventional physics, since photoneutron-driven d–d is expected. The claimed need for enhanced screening (with the predicted 1.2×10³ n/s relying on it) and the 4–5 MeV neutron interpretation are what remain [contested].

### Gaps
- I did not find an independent replication of the NASA Glenn results or a published critique in the fetched sources.
- I could not read the Naranjo Nature full text (paywalled), so the 25 GV/m and ~1000/s figures rest on Wikipedia, and the beam numbers on an abstract snippet.
- Hirsch's 1967 original neutron rates were not read; the ~10⁸ figure is from a lexicon mirror.

## Q6. Baseline tunnelling rate for molecular D₂ at room conditions (Koonin & Nauenberg 1989) and related baselines

### Takeaway
Ground-state molecular fusion rates are p+d ≈ 10⁻⁵⁵ s⁻¹, d+d ≈ 3×10⁻⁶⁴ s⁻¹ and d+t ≈ 10⁻⁶⁹ s⁻¹ per molecule [mainstream]. A mole of D₂ would wait ~10³² years for one fusion (computed). Reaching cold-fusion claim rates requires an effective electron mass of 5–10 m_e. Muonic molecules (m_μ ≈ 207 m_e) fuse at ~10¹² s⁻¹, about 76 orders of magnitude faster (computed).

### Cited Findings
- [mainstream] Koonin & Nauenberg, "Calculated fusion rates in isotopic hydrogen molecules", Nature 339, 690 (1989).
  - They used the Kolos–Wolniewicz Born–Oppenheimer potential, represented by a Morse form for 1.1 < r < 3 a.u. plus a seven-term interpolation at smaller r.
  - Table 2, log₁₀ rate (s⁻¹) at m*/m_e = 1: p+d −55.0, p+t −57.8, d+d −63.5, d+t −68.9. The text gives d+d as "3×10⁻⁶⁴ s⁻¹".
  - p+d rises to −36.0, −19.0 and −10.4 for m*/m_e = 2, 5 and 10.
  - The equilibrium separation scales as m_e/m*.
  - "Enhancements of the electron mass by factors of 5–10 would be required": ≈5 m_e for the Jones et al. claims and ≈10 m_e for Fleischmann–Pons. They know of no plausible mechanism.
  - Muon-catalysed rates are "of order 10¹² s⁻¹".
  — [Koonin & Nauenberg preprint (archived)](https://newenergytimes.com/v2/library/1989/1989KooninNauenberg.shtml)
- [mainstream] The new d+d value was described as about ten orders of magnitude above an earlier estimate. The earlier Brigham Young estimate was ~10⁻⁷⁰ per D₂ molecule per second. — [New Energy Times archive / APS 1989 press materials](https://newenergytimes.com/v2/sr/ColdFusion1989/APS-Baltimore/19890515-AIP-Press-Release.shtml) (via search summary; it describes the increase as "hundred million times", which conflicts with the preprint's "ten orders")
- [mainstream] A 2025 review cites γ_DD = 3×10⁻⁶⁴ s⁻¹ and notes that Koonin–Nauenberg neglect state-dependent parts of the potential. — [arXiv:2501.08338](https://arxiv.org/pdf/2501.08338) (via search summary)
- [mainstream for the baseline; contested for the proposed cascade] Observable fusion at ambient conditions would need "an enhancement of >40 orders of magnitude". The authors argue individual known mechanisms give up to 30 orders each and could be cascaded. — [Metzler, Hunt & Galvanetto, arXiv:2208.07245 (2022)](https://arxiv.org/abs/2208.07245)

### Inferences
- Baselines (computed):
  - d+d at 3×10⁻⁶⁴ s⁻¹: 1.8×10⁻⁴⁰ s⁻¹ per mole of D₂, i.e. one fusion per ~1.7×10³² yr per mole.
  - p+d at 10⁻⁵⁵: one per 5×10²³ yr per mole.
  - d+t at 10^−68.9: one per 4×10³⁷ yr per mole.
  - The ratio of μ-molecular to e-molecular d+d rates is 10¹²/3×10⁻⁶⁴ ≈ 10^75.5.
- Compare the Sun: a core proton lasts ≈6.5×10⁹ yr against p+p (computed, Q1), even though the weak interaction gives p+p an S-factor ~23 orders of magnitude smaller than d+d. S_pp(0) = 4.0×10⁻²² keV b, against S ≈ 55 keV b for D(d,p)T as E→0 from the Bosch–Hale A1 (computed). Temperature at the Gamow peak, not molecular confinement, is what makes stellar fusion go.
- Claims of d–d fusion in metals at rates requiring m* ≈ 5–10 m_e (Fleischmann–Pons-type excess heat attributed to d–d fusion) are generally regarded as [refuted]. This label is my inference from the K–N calculation and the absence of matching neutron or tritium yields; no consensus review (e.g. DOE 1989/2004) was fetched here. Enhanced-screening beam experiments (Q3) are a separate, [contested] topic.

### Gaps
- The Jones et al. (1989) claimed rate (per d–d pair per second) and Fleischmann–Pons yields were not retrieved from primary sources, so those numbers are omitted.
- The Nature version of Koonin–Nauenberg was not accessible; the values are from the archived preprint, which may differ slightly from the published table.
