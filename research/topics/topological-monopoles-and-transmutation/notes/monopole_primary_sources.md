# Primary sources for `magnetics::monopole` (issue #24)

Read 2026-10-27, before `magnetics::monopole`'s remaining functions were coded (#24, #8). Every number is recomputed by [`monopole_sources_checks.py`](monopole_sources_checks.py) (stdlib, `python3 -I`; all checks pass). The papers themselves are not stored here.
Scanned 1980s preprints were read from the KEK preprint-scan server (IDs from INSPIRE `KEKSCAN`); equations were read from the page images, not OCR.

| Source | Status |
|---|---|
| Milton, Rep. Prog. Phys. 69, 1637 (2006), arXiv:hep-ex/0602040 | full text opened |
| Bracci & Fiorentini, NPB 232, 236 (1984), preprint IFUP TH-83/2 | full scan opened (KEK 8302263) |
| Olaussen & Sollie, NPB 255, 465 (1985), preprint CERN-TH-4037/84 | full scan opened (KEK 8412217) |
| Goebel, Monopole '83 p. 333, MAD/TH/146 | full scan opened (KEK 8401109) |
| Olaussen, Olsen, Overbo, Osland, PRL 52, 325 (1984), DESY-83-092 | full scan opened (KEK 8311038), only checked for 263 keV |
| Sivers, PRD 2, 2048 (1970) | **abstract only** (APS paywall; INSPIRE abstract) |
| Patrizii & Spurio, ARNPS 65, 279 (2015), arXiv:1510.07125 | full text opened |
| Ahlen, PRD 17, 229 (1978); Ahlen & Kinoshita, PRD 26, 2347 (1982) | **abstracts only** (APS paywall; OSTI pages unreachable). Formulas taken as reprinted by Milton §8 and Patrizii–Spurio |
| PDG 2025 rev. 94 "Magnetic Monopoles" (Milstead & Weinberg, rev. Aug 2025), https://pdg.lbl.gov/2025/reviews/rpp2025-rev-mag-monopole-searches.pdf | opened; also PDG 2024 rev. 95 (rev. Mar 2024) for comparison |
| Super-K, Ueno et al., arXiv:1203.0940 | opened. **Journal is Astropart. Phys. 36, 131 (2012)**, doi:10.1016/j.astropartphys.2012.05.008 — not PRD 85 |
| IceCube, EPJC 74, 2938 (2014), arXiv:1402.3460 | opened |
| PDG QCD review 2025 (rev. Aug 2025), 2024, 2022, 2020, 2018 | opened. Only the 2018 edition prints Λ values |

---

## 1. Monopole–nucleus binding (Milton 2006 §9)

### Units and definitions (Milton)
- Gaussian units, α = e²/ħc; "natural units ħ = c = 1" in field-theory contexts (§1, p. 2).
- (3.1) e g = m′ħc, m′ integer or half-integer; (3.2) dyons e₁g₂ − e₂g₁ = −m′ħc.
- (9.2) spin ½: H_S = −(eħ/2µc) γ σ·B, γ = 1 + κ = g/2.
- (9.3) µ = (e/m_p c) γ S, nucleus charge Ze, mass M = A m_p (monopole mass ≫ M).
- (9.4) **l = |m′| Z**.
- (9.6) **γ̂ = (A/Z) γ, κ̂ = γ̂ − 1**.
- J = total angular momentum of the nucleus–monopole system (orbital + field + spin); lowest J = l − S.

**Convention trap:** (9.3) taken literally gives γ = µ/(2S µ_N), but Milton's Table 1 lists γ = µ/µ_N (²H 0.857, ²⁷Al 3.63), and the criteria come from Olaussen–Sollie with g_M = Aµ/(SZ) (their 4.17, µ in µ_N). The γ̂ that makes 9.7–9.13 correct is
**γ̂ = S·g_M = (A/Z)·µ/µ_N** (µ is the tabulated moment in nuclear magnetons). For S = ½ the two readings agree. For S > ½, (9.3) is off by a factor 2S. Our check script uses (A/Z)µ/µ_N, which is **correct**.

### Equations as printed (Milton)
- (9.1) Dyon–dyon Coulomb binding, q = e₁e₂ + g₁g₂ < 0, m′ = −(e₁g₂ − e₂g₁)/ħc, from H_NR (3.126) (ħ = 1, µ = reduced mass):
  E_Nj = −(µq²/2) [N + ½ + ((j + ½)² − m′²)^½]⁻²; N principal quantum number; j ≥ |m′|.
- (9.5) neutron (Z = 0): binds in J = ½ if |γ| > 3/(4|m′|).
- (9.7) lowest state J = l − ½ (and J = l − S for any S, §9.2): **γ̂ > 1 + 1/(4l)**.
- (9.8) higher J (spin ½): **|κ̂| > κ_c = (J² + J − l²)/l**; e.g. J = l + ½ → κ_c = 2 + 3/(4l).
- (9.9) next state J = l − S + 1: λ± > ¼ with λ± = (S − ½)(γ̂/S) l − 2l − 1 ± [(1 + l)² + (2S − 1 − l)(γ̂/S) l + ¼((γ̂/S) l)²]^½.
- (9.10) S = 1: γ_c = (3/4l)(3 + 16l + 16l²)/(9 + 4l).
- (9.11a,b) λ+ > ¼ for γ̂ > γ_c⁻; λ− > ¼ for γ̂ > γ_c⁺.
- (9.12) S = 3/2: **(γ_c)∓ = (3/4l)(6 + 4l ∓ √(33 + 32l))**. Text: ⁹Be (γ̂ = −2.66) cannot bind; 3 > γ_c⁻ > 1.557, 3 < γ_c⁺ < 8.943 over 2|m′| = 1…∞.
- (9.13) S = 5/2: **(γ_c)∓ = (36 + 28l ∓ √(1161 + 1296l + 64l²))/(12l)**. Text: ²⁷Al (γ̂ = 7.56) binds; 1.67 > γ_c⁻ > 1.374, 1.67 < γ_c⁺ < 4.216.
- Validity: non-relativistic, point dipole, l ≥ S (exceptions ²H, ⁸Li, ¹⁰B with |n| = 1); binding *energies* need a short-distance regulator (hard core / form factor) and are "very sensitive" to it. The criteria themselves are regulator-independent (O–S §1).

### Table 1 (Milton) — γ, γ̂, J, E_b, notes, ref
n ½ −1.91 — ½, 350 keV NR,hc [146 Sivers] · ¹H ½ 2.79 2.79 l−½=0: 15.1 keV NR,hc [142 B–F]; 320 keV NR,hc [146]; 50–1000 keV NR,FF [147 O–S]; 263 keV R [148,149] · ²H 1 0.857 1.71 l−1=0 (|m′|=1) 130/λ keV R,IM [151,152] · ³He ½ −2.13 −3.20 l+½=3/2 13.4 keV NR,hc [142] · ²⁷Al 5/2 3.63 7.56 l−5/2=4: 2.6 MeV NR,FF [148,149]; 560 keV NR,hc [153 Goebel] · ¹¹³Cd ½ −0.62 −1.46 l+½=49/2 6.3 keV NR,hc [142].

### Where each piece comes from (originals opened)
- **9.7 and 9.8 (spin ½)**: Bracci & Fiorentini [142] eqs 2.19 and 2.27–2.28 (Gaussian-like, ħ = c = 1, n = ±1 only). µ = eK S/m (2.3), m = nuclear mass, so **K = A·µ/µ_N**; q = Ze/(2e)·n = Z/2.
  - (2.19) λ_I = |q| − K/2 for J = |q| − ½; bound iff λ_I < −¼ ⇔ K > 2|q| + ½ ⇔ γ̂ > 1 + 1/(4l). ✔ = Milton 9.7.
  - (2.24) 2×2 matrix; (2.27) λ_II = (j² − q²) − (j² − |q|K + K²/4)^½, j = J + ½; bound iff λ_II < −¼, with (2.28a) j² ≥ (|q| + 1)². Algebra (and numerical root-finding in the script) gives exactly |κ̂| > (J² + J − l²)/l. ✔ = Milton 9.8.
  - **Sign convention differs from O–S**: B–F bind for λ < −¼; O–S bind for λ > ¼.
  - B–F Table I recomputed (−λ_I, −λ_II): ¹H 0.895 (0.90), ¹³C 1.55 (1.57), ¹⁹F 20.49 (20.47), ³¹P 10.02 (10.04), ³He 1.539 (1.54), ¹¹³Cd 10.44 (10.48). B–F binding energies match Milton: ¹H 15.1 keV, ³He 13.4 keV, ¹¹³Cd 6.3 keV ✔. **Printed typo found**: ³H λ_I printed 2.97; formula gives 3.97 (the same row's λ_II = 2.21 confirms K).
- **9.9–9.13 (general S)**: Olaussen & Sollie [147] §4, ħ = c = 1.
  - (1.2) µ = g_M (Ze/2M) S; (4.17) **g_M = Aµ/(SZ)**, µ in nuclear magnetons; q = Zeg, |q| = |nZ/2| (= Milton l), X = g_M|q| = γ̂ l/S.
  - (4.16a) j = |q| − S: **λ₀ = S X − |q|** (binding iff λ₀ > ¼ ⇔ Milton 9.7 for every S).
  - (4.16b) j = |q| + 1 − S: **λ± = (S − ½)X − 2|q| − 1 ± √((|q| + 1)² + (2S − 1 − |q|)X + X²/4)** = Milton 9.9.
  - (4.14) E_B ≈ −λ_n/(2M r₀²) ∼ A^(−2/3), r₀ ≈ A^(1/3) fm (order of magnitude).
  - Solving λ± = ¼ numerically reproduces 9.10, 9.12, 9.13 to 1e−9 at l = 2, 4, 6.5, 20. ✔
  - O–S Table (|eg| = ½) recomputed: ²⁷Al S 5/2, µ 3.6414, |q| 6.5 → X 19.66 (19.7), λ₀ 42.66 (42.8), λ+ 35.51 (35.6), λ− 15.14 (15.2), E_B 2.6 MeV. ²³Na 17.00/20.00/12.42/−2.42 (17.0/20.0/12.4/−2.4). ⁵⁵Mn, ⁸⁵Rb, ¹³⁵Ba, ¹³⁸La all within rounding.
- **Goebel [153]**: proton λ 0.90; ²⁷Al |q| 6.5, S 2.5, X 19.5, ground-state binding 0.56 MeV (hard core) ✔ Milton's 560 keV.
- **Olaussen et al. PRL 52 [149]**: monopole–proton binding 263 keV (plus a 938 MeV E = 0 state) ✔.
- **Sivers [146]**: abstract only; Milton says it "seems to contain errors". Its 350 keV / 320 keV not checked.

### Milton errors / caveats found
1. (9.3) vs Table 1: γ convention differs by 2S for S > ½ (see above). Use γ̂ = (A/Z)µ/µ_N.
2. §9.2 text "1.67 < γ_c⁺ < 4.216" for Al: (9.13) gives γ_c⁺ → 3 as l → ∞ (γ_c⁻ → 5/3). The tight range is 3 < γ_c⁺ ≤ 4.216. Conclusion (Al binds, γ̂ = 7.56) unchanged.
3. (9.12) limits: both γ_c∓ → 3 as l → ∞ ✔ as printed.

### Comparison with `monopoles_checks.py` §6–7
`ghat`, `lowest_binds`, `kappa_c`, `gc32`, `gc52`, `E_dyonic` all **match** the originals. Recomputed: γ_c(l=2, S=3/2) = 1.55668/8.94332; γ_c(l=6.5, S=5/2) = 1.37364/4.21610; E(Z=1, N=0, j=½, m′=½) = −7.29127 eV (CODATA Ry 13.605693123 eV). The check's comment "eq 9.13 large-l limit γ_c⁻ → 5/3 (published 1.67)" is right; Milton's 1.67 for γ_c⁺ is the slip.
Note for Rust: 9.12/9.13 only cover the two lowest J; a negative γ̂ (⁹Be, ¹⁷O) fails them, but the S = ½ rule 9.8 shows negative moments can bind in higher J (³He, ¹¹³Cd). The general-S higher-J case needs the full O–S matrix (4.15), which no one tabulated, so return "not determined" there rather than "no binding".

### Test vectors (Rust)
| fn | input | expect | tol |
|---|---|---|---|
| γ̂ | ²⁷Al Z 13, A 27, µ 3.6414 | 7.5629 | 1e-4 |
| λ₀, λ± (O–S 4.16) | ²⁷Al S 2.5, Z 13, A 27, µ 3.6414, n 1 | 42.659, 35.513, 15.141 | 1e-4 (printed 42.8/35.6/15.2, 0.5%) |
| λ₀, λ± | ²³Na S 1.5, Z 11, A 23, µ 2.2175 | 20.001, 12.417, −2.416 | 1e-4 |
| γ_c (9.12) | l = 2 | 1.556683, 8.943317 | 1e-6 |
| γ_c (9.13) | l = 6.5 | 1.373638, 4.216105 | 1e-6 |
| γ_c (9.10) | l = 2 | 2.183824 | 1e-6 |
| κ_c (9.8) | J 1.5, l 1 / J 24.5, l 24 | 2.75 / 2.03125 | exact |
| E (9.1) | Z 1, N 0, j ½, m′ ½ | −7.29127 eV (0.535898 Ry) | 1e-5 |
| B–F λ_II | ³He Z 2, A 3, µ −2.13, j 2 | −1.5385 | 1e-4 |

---

## 2. Monopole stopping power

### Patrizii & Spurio, arXiv:1510.07125, §III A–B (Gaussian units; g = magnetic charge in Gaussian, g_D = e/2α = 68.5e)
- **Fast, β > 0.05 and γ ≤ 100 (Eq. 11)**, from Ahlen [28] by Ze → gβ:
  dE/dx = (4πN_e g²e²/m_ec²) [ln(2m_ec²β²γ²/I) + k/2 − ½ − δ/2 − B_m]
  N_e electron density, I mean ionisation potential, δ density effect, k QED (KYG) correction, B_m Bloch correction; "(Ze)² = g_D²β² ∼ (68.5β)²". Above βγ ≈ 10⁴ radiative processes dominate.
- **Slow, 10⁻³ ≲ β ≲ 10⁻² (Eq. 12)**, Ahlen & Kinoshita [32], free degenerate electron gas (good for metals; "reasonable" for Z ≥ 10):
  dE/dx = (2πN_e g²e²β/(m_e c v_F)) [ln(2m_e v_F a₀/ħ) − ½], v_F = (ħ/m_e)(3π²N_e)^(1/3) ≃ 3.9×10⁸ cm/s for **Si at the Earth-mantle density 4.3 g/cm³** (not lab Si, 2.33).
  **"×1.37" source**: P&S state that, per [32], coupling to the electron magnetic moment "gives a multiplicative factor of 1.37 to Eq. 12". For conductors an extra conduction-electron term is to be added. Visible-photon transitions are largely suppressed, so Eq. 12 does not give scintillator light.
- Fig. 2 (Si 4.3 g/cm³, Fe 11.5 g/cm³) curves for n = 1, 2, 3; not digitised.

### Milton §8 (same constants, the only open source for them)
- (8.10) −dE/dx = (4πg²e²/c²m_e) N_e [ln(2m_ec²β²γ²/I) + K(|n|)/2 − δ/2 − ½ − B(|n|)], valid β > 0.1.
  **K = 0.406 (n = 1), 0.346 (n ≥ 2)** (Kazama–Yang–Goldhaber); **B = 0.248, 0.672, 1.022, 1.685 for n = 1, 2, 3, 6** (Bloch). Scale by n² apart from K, B.
- (8.11) for β < 0.01, from eq. (60) of Ahlen–Kinoshita "as an approximation for all materials": **−dE/dx = (45 GeV/cm) n² β**. An empirical polynomial joins 0.01 < β < 0.1.

### SI form for Rust (exact identity, no α left)
With g = n e/(2α) in Gaussian, g²e² = n²(ħc)²/4:
- fast: dE/dx = π N_e n² (ħc)²/(m_ec²) · [ln(2m_ec²β²γ²/I) + K/2 − δ/2 − ½ − B]
- slow: dE/dx = 1.37 · (π/2) N_e n² (ħc)² β/(m_e c v_F) · [ln(2v_F/(αc)) − ½]  (2m_e v_F a₀/ħ = 2v_F/(αc))
- Equivalently prefactor = K_Bethe (Z/A) ρ (n/2α)², K_Bethe = 0.307075 MeV cm²/mol.

### Comparison with `strange_radiation_checks.py` §4
`dEdx_slow_GeV_cm` **matches** Eq. 12 × 1.37 (it takes g²e² = (n/2α)² e⁴ with e² = 1.43996448×10⁻⁷ eV·cm, correct Gaussian). `fermi_velocity_cm_s` reproduces 3.900×10⁸ cm/s ✔. The fast formula is only described, not coded; K and B were not in the notes and should be added.
**Cross-check**: Eq. 12 alone in Si (2.33 g/cm³) at β = 10⁻³ gives 44.8 MeV/cm, the Ahlen–Kinoshita "45 GeV/cm·β" value; ×1.37 gives 61.3 MeV/cm. At 4.3 g/cm³: 91.6 / 125.5 MeV/cm. So the repo's 61 MeV/cm is Eq. 12 × 1.37, about 1.4× Milton's all-materials approximation.

### Test vectors (Rust)
| fn | input | expect |
|---|---|---|
| v_F | Si, ρ 4.3 g/cm³ (Z 14, A 28.0855) | 3.8996×10⁸ cm/s (P&S 3.9×10⁸) |
| v_F | Si, ρ 2.33 | 3.17916×10⁶ m/s; N_e 6.99444×10²⁹ m⁻³ |
| slow ×1.37 | Si 2.33, β 10⁻³, n 1 | 61.316 MeV/cm = 9.8239×10⁻¹⁰ J/m |
| slow ×1 | Si 2.33, β 10⁻³ | 44.76 MeV/cm (≈ Milton 8.11: 45 MeV/cm, within 1%) |
| fast | Si 2.33, I 173 eV, δ 0, n 1, β 0.5 | 11.788 GeV/cm = 1.88869×10⁻⁷ J/m |
| fast | same, β 0.1 / 0.2 / 0.8 | 5.934 / 8.307 / 14.591 GeV/cm |
| (g_D/e)² | — | 4694.7 |

No printed point value of a full fast-monopole dE/dx curve was found in an open source (Ahlen 1978 and the MACRO/Derkaoui tables are paywalled); the 45 GeV/cm·β line is the only printed number.

---

## 3. PDG "Magnetic Monopoles"

Used: **PDG 2025 update, review 94** (Milstead & Weinberg, revised August 2025; in S. Navas et al., PRD 110, 030001 (2024) and 2025 update). PDG 2024 numbers the same equations **95.1–95.12** (revised March 2024); formulas identical. 2025 adds the colour-magnetic (Z₃) discussion and a ref. Natural Heaviside–Lorentz units (e = √(4πα)).

- (94.1) **B = Q_M r̂/(4πr²)**.
- (94.2) **Q_E^min Q_M^min = 2π**; Dirac charge Q_M^D ≡ 2π/e; 6π/e if no colour-magnetic charge.
- (94.3) L = m r × v − (1/4π) Q_E Q_M r̂ (PDG 2024's text layer drops the 1/4π fraction bar; the 2025 PDF shows it).
- (94.4) M_mon ∼ 4πM_V/e²; (94.6) M_mon ∼ 4πM_X/g², "Mmon ∼ 10¹⁷–10¹⁸ GeV" for 10¹⁶ GeV unification.
- (94.8) **v_mag ∼ c for M ≲ 10¹¹ GeV; 10⁻³c (10¹⁷ GeV/M)^½ for M ≳ 10¹¹ GeV**, B ∼ 3 µG.
  **There is no "energy gained in a field" equation in PDG**; 94.8 is a velocity. Its second branch is ½Mv² = 5×10¹⁰ GeV for every M (independent of M), i.e. it encodes a fixed energy gain ≈ g_D·B·L. g_D × 3 µG × 1 kpc = 1.90×10¹¹ GeV; 10¹¹ GeV needs L = 0.526 kpc (derived, as in the notes). So `monopole_energy_gain = n g_D B L` is derived from 94.1 + work, not printed by PDG; cite it as derived.
- (94.9) **Parker: F < 10⁻¹⁵ cm⁻² sr⁻¹ s⁻¹ (M ≲ 10¹⁷ GeV); F < 10⁻¹⁵ (M/10¹⁷ GeV) cm⁻² sr⁻¹ s⁻¹ (M ≳ 10¹⁷ GeV)**.
- (94.10) **Seed (extended Parker): F < [M/(10¹⁷ GeV) + 3×10⁻⁶] × 10⁻¹⁶ cm⁻² sr⁻¹ s⁻¹** (confirmed on the rendered page; bracket encloses both terms). Floor 3×10⁻²² as M → 0; exactly Parker/10 at 10¹⁷ GeV; linear in M at all masses.
- (94.11) F_uniform = 1.3×10⁻¹⁶ Ω_M (10¹⁷ GeV/M)(v/10⁻³c).
- (94.12) catalysis bounds F (σ_ΔB β/10⁻²⁷ cm²) ∼ 10⁻¹⁸–10⁻²⁹; σ_ΔB β ∼ 10⁻²⁷ cm².
- §94.5.3 Super-K quoted as 6.3×10⁻²⁴ (β/10⁻³)² — PDG drops the σ₀ and f_π factors.

### Super-Kamiokande 2012 (arXiv:1203.0940; Astropart. Phys. 36, 131)
- (6) f_p = ∫ n_M v_rel σ ρ_p N_A d³x, **σ = σ₀/β_rel²** (σ₀ ∼ 1 mb), β_rel = 1.7×10⁻³ (thermal protons, T_c = 1.5×10⁷ K; the footnote's "m/s" is a typo for units of c), ρ_p = 50 g/cm³.
- (7) N_M (σ₀/1 mb)(f_π⁺/0.5) < 7.1×10¹⁷ (90% C.L.).
- (8) N_M = πR⊙² (1 + (β_esc/β_M)²) ε(M_M, β_M) 4πF_M t⊙, β_esc = 2×10⁻³, R⊙ = 7.0×10¹⁰ cm, t⊙ = 4.6×10⁹ yr.
- (9) **F_M (σ₀/1 mb)(f_π⁺/0.5) < 6.3×10⁻²⁴ (β_M/10⁻³)² cm⁻² s⁻¹ sr⁻¹ (90% C.L.)**, valid for **β_M < 10⁻³** and **M_M ≲ 10¹⁷ GeV** (capture efficiency ε drops above).
- Recomputed: (7)+(8) with ε = 1 and the small-β form (β_esc/β)² gives 6.32×10⁻²⁴ at β = 10⁻³ ✔. With the full 1 + (β_esc/β)² it is 5.06×10⁻²⁴, so Eq. 9 is the β ≪ β_esc limit (25% loose, conservative, at β = 10⁻³). β² scaling comes from gravitational focusing, not from σ.
- **Mismatch with `superk_limit`** in monopoles_checks.py: it omits the (f_π⁺/0.5) factor and the validity limits. Formula otherwise matches.

### IceCube 2014 (arXiv:1402.3460; EPJC 74, 2938)
- (4) **σ_cat = σ₀/β for β ≥ β₀; σ₀/β · F(β) for β < β₀, F(β) = (β/β₀)^γ** (γ, β₀ nucleus-dependent; sign of γ enhances or suppresses).
- (5) **λ_cat = 1/(σ_cat n)**, n = "particle density of the medium" (n not given numerically).
- σ_cat = 1.7×10⁻²² / 1.7×10⁻²³ cm² ↔ λ_cat = 1 mm / 1 cm; estimates 10⁻²⁷–10⁻²¹ cm².
- Our `catalysis_length = β/(σ₀ n)` **matches** for β ≥ β₀; the β < β₀ branch is missing (add F(β) with β₀, γ inputs). Note the β-scaling differs from Super-K's σ₀/β².
- Implied n = 5.88×10²² cm⁻³; ice H nuclei (0.917 g/cm³) 6.131×10²² give λ = 0.959 mm / 0.959 cm (4% off; inference, not stated by IceCube).

### Test vectors (Rust)
| fn | input | expect |
|---|---|---|
| parker_bound(M GeV) | 10¹⁵ / 10¹⁷ / 10¹⁹ | 1e-15 / 1e-15 / 1e-13 cm⁻² sr⁻¹ s⁻¹ |
| seed_bound(M) | 10¹⁵ / 10¹⁷ / 10¹⁹ / 0 | 1.0003e-18 / 1.000003e-16 / 1.00000003e-14 / 3e-22 (the 10¹⁹ entry first read 1.0000003e-14; corrected when the Rust test failed) |
| v_mag/c (94.8) | 10¹⁰ / 10¹³ / 10¹⁷ GeV | 1 / 0.1 / 1e-3 |
| superk_limit(β, σ₀ mb, f_π) | (1e-3, 1, 0.5) / (1e-4, 1, 0.5) / (1e-3, 10, 0.5) | 6.3e-24 / 6.3e-26 / 6.3e-25 |
| sigma_cat(σ₀, β) | (1.7e-25 cm², 1e-3), β ≥ β₀ | 1.7e-22 cm² |
| catalysis_length | σ_cat 1.7e-23 cm², n 6.131e22 cm⁻³ | 0.959 cm |
| Q_M^D/e = 2π/e² (HL) | — | 1/(2α) = 68.518 = g_D/(ec) |

---

## 4. Λ_QCD

- **PDG 2025 QCD review (rev. Aug 2025) prints no Λ value.** It quotes only α_s(m_Z²) = **0.1180 ± 0.0009** (eq. 9.25, "PDG 2025 average"), and says it has "become standard practice" to quote α_s rather than Λ (§9.1.1). 2024, 2022 and 2020 editions likewise print no Λ.
- Last printed values: **PDG 2018 QCD review (S. Bethke et al., rev. Sept 2017), eqs 9.24a–d**, for α_s(M_Z²) = 0.1181 ± 0.0011, MS-bar, 4-loop running (their eq. 9.5) with 3-loop matching at m_c, m_b, m_t = 1.3, 4.2, 173 GeV:
  Λ^(6) = 89 ± 6 MeV; **Λ^(5) = 210 ± 14 MeV**; **Λ^(4) = 292 ± 16 MeV**; **Λ^(3) = 332 ± 17 MeV**.
  URL: https://pdg.lbl.gov/2018/reviews/rpp2018-rev-qcd.pdf
- Recomputed with the 4-loop formula (b₀…b₃, MS-bar): α_s = 0.1181 → Λ^(5) = 209.5 MeV ✔ (vs 210). With PDG 2025's 0.1180 at m_Z = 91.188 GeV → **Λ^(5) ≈ 208.4 MeV (derived, not printed)**.
- Suggested default: `LAMBDA_QCD_NF3 = 332 MeV` cited as "PDG 2018 eq. 9.24d" (n_f = 3 is the regime relevant to deconfinement / strong-field QCD). The notes' 200–340 MeV range is consistent with these.

---

## Summary of mismatches
1. Milton (9.3) γ definition vs Table 1/criteria (factor 2S for S > ½); our code uses the right one (A/Z)µ/µ_N.
2. Milton §9.2 "1.67 < γ_c⁺" should be 3 < γ_c⁺ (limit of 9.13). Conclusion unchanged.
3. Bracci–Fiorentini Table I ³H λ_I misprint (2.97 for 3.97). Not used by us.
4. `superk_limit` lacks the (f_π⁺/0.5) factor and its validity range (β < 10⁻³, M ≲ 10¹⁷ GeV). Journal is Astropart. Phys. 36, 131, not PRD 85.
5. `catalysis_length` lacks IceCube's β < β₀ branch.
6. PDG has no energy-gain equation; 94.8 is v_mag. `monopole_energy_gain` should be cited as derived.
7. PDG 2025 gives no Λ; cite PDG 2018 eq. 9.24.
8. Ahlen 1978 / Ahlen–Kinoshita 1982 not opened (paywall). The 1.37, K and B values rest on Patrizii–Spurio and Milton.
