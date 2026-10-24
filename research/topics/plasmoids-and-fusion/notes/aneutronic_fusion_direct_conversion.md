# Aneutronic fusion and direct energy conversion: formulas, test values and physical limits

Conventions. T in keV unless stated; ⟨σv⟩ in m³/s or cm³/s as stated; powers in W m⁻³. Labels: **[mainstream]** = standard physics or a peer-reviewed measurement that nobody disputes; **[contested]** = published but disputed, model-dependent, company-reported, or a projection that has not been measured; **[refuted]** = contradicted by later analysis or data. "computed" = recomputed in `aneutronic_checks.py` (in this folder; standard library only; `python3 -I aneutronic_checks.py` runs in about 1 s, and every published anchor is an `assert`).

## 1. Aneutronic and low-neutron reactions: Q-values, products, neutron fractions

### Takeaway
The Q-values all recompute from NIST atomic masses to better than 0.1%. "Aneutronic" fuels are only *low*-neutron. In p–¹¹B, side reactions put under about 0.2% of the energy into neutrons. In D–³He, the D–D side reactions put about 0.3–1.3% into neutrons, rising to about 2–7% if the bred tritium also burns (computed).

### Cited Findings
- [mainstream] The reactions and their energy release: D+³He → ⁴He + p + 18.3 MeV; D+⁶Li → 2⁴He + 22.4 MeV; p+⁶Li → ⁴He + ³He + 4.0 MeV; ³He+⁶Li → 2⁴He + p + 16.9 MeV; ³He+³He → ⁴He + 2p + 12.86 MeV; p+⁷Li → 2⁴He + 17.2 MeV; p+¹¹B → 3⁴He + 8.7 MeV; p+¹⁵N → ¹²C + ⁴He + 5.0 MeV. Wikipedia attributes this table to Harms & Schoepf (2000) — [Wikipedia: Aneutronic fusion](https://en.wikipedia.org/wiki/Aneutronic_fusion)
- [mainstream] Putvinski et al. take a yield of Y = 8.68 MeV for p–¹¹B — [Putvinski, Ryutov & Yushmanov 2019, Nucl. Fusion 59 076018 (TAE copy, p.2)](https://tae.com/wp-content/uploads/fliphtml5/139/files/basic-html/page2.html); Liu et al. also use 8.68 MeV — [Liu et al. 2024, arXiv:2405.13260](https://arxiv.org/html/2405.13260v1)
- [mainstream] The atomic masses used in the check, in u: ¹H 1.00782503223(9); ²H 2.01410177812(12); ³H 3.0160492779(24); ³He 3.0160293201(25); ⁴He 4.00260325413(6); ⁶Li 6.0151228874(16); ⁷Li 7.0160034366(45); ¹¹B 11.00930536(45); ¹⁴N 14.00307400443(20) — [NIST Atomic Weights and Isotopic Compositions](https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&all=all&ascii=ascii2&isotype=some)
- [mainstream] The p–¹¹B neutron sources:
  - ¹¹B + α → ¹⁴N + n + 157 keV is the main source. At least 0.1% of reactions produce a neutron, and these neutrons carry under 0.2% of the total energy (Heindler & Kernbichler 1989).
  - ¹¹B(p,n)¹¹C has Q = −2.8 MeV, so it needs energetic protons.
  - ¹¹B(p,γ) has a branching ratio of about 10⁻⁴, giving 4, 12 and 16 MeV gammas.
  - The definition of "aneutronic" is ≤1% of the energy in neutrons (New Jersey A2731, 2006).
  - Source: [Wikipedia: Aneutronic fusion](https://en.wikipedia.org/wiki/Aneutronic_fusion)
- [mainstream] Wang et al. give the ¹¹B(p,n)¹¹C channel a branching of about 10⁻⁵, with the threshold quoted as about 2.765 MeV — [Wang, Li, Wu & Cui 2026, arXiv:2601.00241](https://arxiv.org/html/2601.00241v1)
- [contested — no source on the page] Wikipedia says the neutron energy fraction of D–³He is "probably several percent" — [Wikipedia: Aneutronic fusion](https://en.wikipedia.org/wiki/Aneutronic_fusion)
- [mainstream, but no per-value source on the page] Ideal ignition temperatures: D–T 13.6 keV, D–³He 58 keV, p–⁶Li 66 keV, p–¹¹B 123 keV. Reactivity peaks: p–¹¹B at about 600 keV, D–T at about 66 keV (Lerner & Terry 2007) — [Wikipedia: Aneutronic fusion](https://en.wikipedia.org/wiki/Aneutronic_fusion)

### Inferences
Formulas, all computed and asserted in the script.
- **Q-value.** `Q = (Σ M_atomic(reactants) − Σ M_atomic(products)) × 931.49410242 MeV/u`. Electron counts balance in every reaction listed, so atomic masses can be used directly.
- **Two-body product split** (zero initial kinetic energy, non-relativistic): `E_a = Q · m_b/(m_a + m_b)`.

Computed Q-values, MeV (published value in brackets):

| Reaction | Q, computed |
|---|---|
| p+¹¹B → 3α | 8.6821 (8.68) |
| D+³He → α+p | 18.3531 (18.35) |
| ³He+³He → α+2p | 12.8596 (12.86) |
| p+⁶Li → α+³He | 4.0197 (4.0) |
| D+D → T+p | 4.0327 |
| D+D → ³He+n | 3.2689 |
| D+T → α+n | 17.5893 |
| p+⁷Li → 2α | 17.3462 |
| D+⁶Li → 2α | 22.3728 |
| α+¹¹B → ¹⁴N+n | 0.1581 (0.157) |
| p+¹¹B → ¹¹C+n | −2.7648 (−2.8). This uses an ¹¹C mass recalled from AME2020 (mass excess 10650.3 keV), not fetched. |

Product energies, computed:
- D–³He: p 14.66 MeV, α 3.69 MeV.
- p–⁶Li: α 1.73 MeV, ³He 2.29 MeV.
- D–D: T 1.01 + p 3.02 MeV; ³He 0.82 + n 2.45 MeV.
- p–¹¹B: mean α energy Q/3 = 2.89 MeV. This is a three-body continuum, not one line.

Gamow energy, `E_G = 2 m_r c² (π α Z₁ Z₂)²`, computed:
- p–¹¹B: 22601 keV. Nevins–Swain use 22589 keV, because nuclear rather than atomic masses give a slightly smaller m_r.
- D–³He: √E_G = 68.77 keV^½ (Bosch–Hale: 68.7508).
- ³He–³He: E_G = 23625 keV, about the same barrier as p–¹¹B but without its 0.6 MeV resonance.
- p–⁶Li: E_G = 7607 keV.

D–³He neutron power fraction, computed from the Bosch–Hale fits. 2.45 MeV D–D neutrons only, and in brackets the value if every bred T burns at once (adding 14.07 MeV neutrons):

| T | 50:50 D:³He | 1:3 D:³He |
|---|---|---|
| 50 keV | 1.31% (7.3%) | 0.45% (2.6%) |
| 70 keV | 1.05% (5.8%) | 0.36% (2.1%) |
| 100 keV | 1.01% (5.5%) | 0.34% (2.0%) |

So "several percent" holds only when the tritium burns. Pumping the tritium out, or running ³He-rich, brings the fraction below 1%.

### Gaps
- I found no accessible reactivity fit or table for ³He–³He or p–⁶Li. NACRE II ([arXiv:1310.7099](https://arxiv.org/pdf/1310.7099)) covers ⁶Li(p,α)³He as an astrophysical rate, but I did not open it.
- I did not obtain a sourced value for the p–⁶Li secondary-reaction neutron fraction.
- The ¹¹C mass is recalled, not fetched.

## 2. Cross sections and reactivities (Nevins–Swain, Sikora–Weller, Tentori–Belloni, Bosch–Hale)

### Takeaway
Nevins–Swain and Tentori–Belloni both use the same Bosch–Hale/Peres functional form, and both coefficient sets are published, so a Rust calculator can share one function. The newer p–¹¹B data (Sikora–Weller 2016 and the measurements after it) raise the reactivity by about 12% at 100 keV and about 50% at 500 keV. The computed values reproduce published numbers to 4 significant figures.

### Cited Findings
- [mainstream] **Common form** (Peres / Bosch–Hale), with T in keV:
  ```
  θ  = T / [1 − T(C2 + T(C4 + T·C6)) / (1 + T(C3 + T(C5 + T·C7)))]
  ξ  = (B_G² / (4θ))^(1/3)
  ⟨σv⟩ = C1 · θ · sqrt(ξ / (m_r c² · T³)) · exp(−3ξ)
  ```
  Tentori–Belloni write ξ ≡ (E_G/4θ)^{1/3} with E_G = 22.589 MeV. Their high-temperature (HT) fit covers 70–500 keV, and a separate low-temperature (LT) form covers 10–70 keV — [Tentori & Belloni 2023, Nucl. Fusion 63 086001](https://iopscience.iop.org/article/10.1088/1741-4326/acda4b)
- [mainstream] **p–¹¹B coefficients** from Tentori–Belloni Table 2. P1 is in keV m³ s⁻¹ and the output is in m³/s. The pairs are Nevins & Swain 2000 → Tentori & Belloni 2023:
  - P1: 4.4467e-14 → 9.7827e-14
  - P2: −5.9357e-2 → −5.1610e-2
  - P3: 2.0165e-1 → 1.3240e-1
  - P4: 1.0404e-3 → 3.8446e-4
  - P5: 2.7621e-3 → 1.2499e-3
  - P6: −9.1653e-6 → −5.6715e-6
  - P7: 9.8305e-7 → 1.1615e-6
  - m_r c² = 859526 keV and E_G = 22589 keV — [Tentori & Belloni 2023](https://iopscience.iop.org/article/10.1088/1741-4326/acda4b); the same values appear in [Ahmad et al. 2026, Front. Nucl. Eng.](https://www.frontiersin.org/journals/nuclear-engineering/articles/10.3389/fnuen.2026.1714531/full)
- [mainstream] **Tentori–Belloni vs Nevins–Swain.** TB is "decidedly larger above approximately 30 keV". The gap is 12% at 100 keV, "higher than 20%" at 200 keV, 30% at 300 keV and 50% at 500 keV. With suprathermal effects included it could reach about 50% at about 300 keV — [Tentori & Belloni 2023](https://iopscience.iop.org/article/10.1088/1741-4326/acda4b)
- [mainstream] **Tentori–Belloni S-factor breakpoints** (piecewise, CM energy): 0.400, 0.668 and 9.76 MeV. Resonance parameters in TB vs NS: E0 = 640.9 vs 581.3 keV, δE0 = 85.5 vs 85.7 keV. The narrow 148 keV resonance (A_L = 1.82e4 MeV b, δE_L = 2.35 keV) is carried over from NS — [Tentori & Belloni 2023](https://iopscience.iop.org/article/10.1088/1741-4326/acda4b)
- [mainstream] **Published test values** (Ahmad et al. 2026, Fig. 1 range 75–500 keV, m³/s): NS 2.455e-23 and 3.741e-22; TB 2.777e-23 and 5.619e-22 — [Ahmad et al. 2026](https://www.frontiersin.org/journals/nuclear-engineering/articles/10.3389/fnuen.2026.1714531/full)
- [mainstream] **Sikora & Weller 2016.** The new ¹¹B(p,α)αα rate is about 10–15% above Angulo et al. (NACRE 1999) between 200 and 600 keV — [J. Fusion Energy, doi:10.1007/s10894-016-0069-y](https://link.springer.com/article/10.1007/s10894-016-0069-y). Wang et al. 2026 prefer the newer Mazzucconi et al. data for 0.34–4.73 MeV, and add a resonance at about 4.7 MeV that is negligible below T ≈ 0.4 MeV. A ±10% change in the strength of the 0.6 MeV resonance moves the reactivity by about as much as the spread between existing evaluations — [Wang et al. 2026](https://arxiv.org/html/2601.00241v1)
- [mainstream] **Putvinski et al. 2019.** New cross sections (about 3% error bars) give about +20% reactivity, and kinetic effects in the proton tail add about +10%, for about +30% in total. Burnout of the tail costs under 1% — [Putvinski et al. 2019, pp.2–4](https://tae.com/wp-content/uploads/fliphtml5/139/files/basic-html/page4.html)
- [mainstream] **Bosch & Hale 1992**, "Improved formulas for fusion cross-sections and thermal reactivities", Nucl. Fusion 32, 611–631. It gives R-matrix-based fits for D(d,n)³He, D(d,p)T, T(d,n)⁴He and ³He(d,p)⁴He — [INIS record](https://inis.iaea.org/records/4crxt-2ed10)
- [mainstream] **Bosch–Hale cross section**, as implemented in WarpX: `S(E) = (A1 + E(A2 + E(A3 + E(A4 + E·A5)))) / (1 + E(B1 + E(B2 + E(B3 + E·B4))))` and `σ = S/(E·exp(B_G/√E))`, with σ in mb and E in keV (CM).
  - D–³He: A1–A3 = 5.7501e6, 2.5226e3, 4.5566e1; B1–B3 = −3.1995e-3, −8.5530e-6, 5.9014e-8.
  - D–T: A = 6.927e4, 7.454e8, 2.050e6, 5.2002e4, 0; B = 6.38e1, −9.95e-1, 6.981e-5, 1.728e-4.
  - D–D (n): A = 5.3701e4, 3.3027e2, −1.2706e-1, 2.9327e-5, −2.5151e-9.
  - D–D (p): A = 5.5576e4, 2.1054e2, −3.2638e-2, 1.4987e-6, 1.8181e-10.
  - Source: [WarpX BoschHaleFusionCrossSection.H](https://warpx.readthedocs.io/en/22.11/_static/doxyhtml/_bosch_hale_fusion_cross_section_8_h_source.html)

### Inferences
**Computed p–¹¹B ⟨σv⟩, m³/s.** NS is shown without its 148 keV resonance term; TB is the HT fit.

| T, keV | NS | TB | TB/NS − 1 |
|---|---|---|---|
| 75 | 2.4554e-23 | 2.7771e-23 | — |
| 100 | 6.1526e-23 | 7.0654e-23 | 0.126 |
| 150 | 1.5761e-22 | 1.8664e-22 | 0.176 |
| 200 | 2.4274e-22 | 2.9575e-22 | 0.214 |
| 300 | 3.3852e-22 | 4.3872e-22 | 0.294 |
| 400 | 3.6982e-22 | 5.1469e-22 | 0.390 |
| 500 | 3.7412e-22 | 5.6194e-22 | 0.501 |

- These match Ahmad et al. at 75 and 500 keV to 4 digits, and match TB's stated 12/>20/30/50% gaps. Ahmad et al.'s values therefore omit the NS resonance term.
- The NS narrow-resonance term `5.41e-15 T^-1.5 exp(-148/T)` cm³/s is recalled from the NS paper, not seen in a fetched source. It adds 5.8% at 70 keV, 2.0% at 100 keV and 0.1% at 500 keV.
- The TB HT fit is still rising at 500 keV, its validity edge.

**Computed D-ion reactivities, cm³/s.** These use the Bosch–Hale C1–C7 reactivity coefficients, which I **recalled** from BH 1992 and did not fetch:
- D–³He: B_G = 68.7508, m_r c² = 1124572, C1–C5 = 5.51036e-10, 6.41918e-3, −2.02896e-3, −1.91080e-5, 1.35776e-4, and C6 = C7 = 0.
- D–T: B_G = 34.3827, m_r c² = 1124656, C1–C7 = 1.17302e-9, 1.51361e-2, 7.51886e-2, 4.60643e-3, 1.35000e-2, −1.06750e-4, 1.36600e-5.
- D–D (n): B_G = 31.3970, m_r c² = 937814, C1, C2, C3, C5 = 5.43360e-12, 5.85778e-3, 7.68222e-3, −2.964e-6.
- D–D (p): B_G = 31.3970, m_r c² = 937814, C1, C2, C3, C5 = 5.65718e-12, 3.41267e-3, 1.99167e-3, 1.05060e-5.

The check is independent: Simpson integration of the WarpX/BH *cross-section* fits over a Maxwellian, `⟨σv⟩ = sqrt(8/(π m_r)) (kT)^-3/2 ∫ σ(E) E e^(−E/kT) dE`. The two agree within 1.3% everywhere, except D–³He at 10 keV, where they differ by 3.3%.

| Reaction | T, keV | Fit | Integral |
|---|---|---|---|
| D–³He | 10 | 2.1261e-19 | 2.0592e-19 |
| D–³He | 20 | 3.4821e-18 | 3.4677e-18 |
| D–³He | 50 | 5.5539e-17 | 5.4914e-17 |
| D–³He | 100 | 1.7185e-16 | 1.7227e-16 |
| D–T | 10 | 1.1362e-16 | 1.1418e-16 |
| D–D (n) | 50 | 1.1330e-17 | 1.1376e-17 |
| D–D (p) | 50 | 9.8383e-18 | 9.9695e-18 |

For Rust: D–T at 10 keV = 1.136e-16 cm³/s is a good unit test. It matches the BH table value as I recall it.

### Gaps
- I could not open the Tentori–Belloni low-temperature (10–70 keV) formula; its equations are images. The NS original paper (Nucl. Fusion 40, 865) and its low-temperature branch were not fetched either.
- The Bosch–Hale 1992 Table VII (reactivity coefficients) and Table VIII (tabulated ⟨σv⟩) were not reachable. The UW-FTI copy failed DNS and the NRL Plasma Formulary returned 403. The coefficients above are recalled and validated only by the integral cross-check.
- Validity ranges as I recall them, unverified: BH D–³He reactivity 0.5–190 keV, cross section 0.3–900 keV.
- I did not obtain the Sikora–Weller peak cross-section values.
- I found no paper by "Ostrovsky et al. 2025". The 2025–26 updates I did find are Wang et al. 2026, Ochs et al. 2026 and Ahmad et al. 2026.

## 3. Power balance: bremsstrahlung, synchrotron, P_fus/P_br vs T_e/T_i; Rider; Putvinski; 2025–26 updates

### Takeaway
For a thermal p–¹¹B plasma with T_e = T_i, bremsstrahlung beats fusion at every temperature. The computed best P_fus/P_br is about 0.49 with the Tentori–Belloni reactivity. Fusion can exceed bremsstrahlung only if T_e ≲ 0.33–0.5 T_i. Putvinski calls that self-consistent regime "marginal" ignition; Rider showed that holding it out of equilibrium costs more power than it saves. Unreflected synchrotron emission would be 40–400× bremsstrahlung, so the scheme also needs near-perfect wall reflection or opacity.

### Cited Findings
- [mainstream] **Putvinski et al. 2019, Eq. 7** (fit to QED results, better than 1%), with t = T_e/m_ec²:
  `P_Brem/W0 = 8.511 t^{1/2} [Z_eff(1 + 1.78 t^{1.34}) + 2.12 t (1 + 1.1t + t² − 1.25 t^{2.5})]`
  The first bracket term is electron–ion and the 2.12 t term is electron–electron. The authors say Rider's formula "does not give accurate dependence on Zeff" (it is good at Z_eff ≈ 3). The page OCR is garbled, and the bracket structure is a reconstruction — [Putvinski et al. 2019, p.3](https://tae.com/wp-content/uploads/fliphtml5/139/files/basic-html/page3.html)
- [mainstream] **Putvinski, other results:** ion–electron exchange `P_ei = (3/2) ν_ei (T_i − T_e)` with a relativistic correction (about 10% above the classical rate); T_e typically 100–200 keV; about 10% of the alpha energy goes to electrons; optimum n_B/n_i = 15% and T ≈ 300 keV.
  - Conclusion: "Bremstrahlung radiation does not preclude ignition in pB11 plasmas". Ignition is "marginal", but high Q is reachable over a broad range of T, under near-perfect confinement.
  - [contested] This is the TAE position, not shared by the ICF and Rider-type analyses below.
  - Source: [Putvinski et al. 2019, p.4](https://tae.com/wp-content/uploads/fliphtml5/139/files/basic-html/page4.html)
- [mainstream] **Rider form, as quoted by Wang et al. (Eq. 10):**
  `P_brem/n_e² = 5.172×10⁻⁴³ √T_e {Z_eff[1 + 0.7936 t + 1.874 t²] + (3/√2) t}` MW m³, with Z_eff = Σn_iZ_i²/Σn_iZ_i.
  - With the reference composition f₁ = 1/2, f₂ = 1/10 and T_e = T, bremsstrahlung exceeds fusion over the whole range.
  - With a self-consistent T_e, fusion exceeds bremsstrahlung between T ≈ 0.17 and 0.38 MeV, peaking near 0.25 MeV.
  - Source: [Wang et al. 2026](https://arxiv.org/html/2601.00241v1)
- [mainstream] **Ochs, Kolmes, Glasser & Fisch 2026** (Phys. Plasmas 33, 012703), the inertial-fusion case:
  - Model: total bremsstrahlung from the Svensson 1982 fit (≲2%); spectrum from Bethe–Heitler with the Elwert correction; Sikora–Weller ⟨σv⟩ (≲3.5%) times a kinetic enhancement.
  - With T_e = T_i: fusion exceeds bremsstrahlung above 150 keV for f_b = 0.1, around 200 keV for f_b = 0.3, and never for f_b = 0.5.
  - Without reabsorption the optimum is T_i ≈ 160 keV, f_b = 1/4, Q_sci ≈ 2.
  - Breakeven needs n_eR ≳ 3×10²⁵ cm⁻² (about 100 g/cm²) and n_e²R ≳ 10⁵³ cm⁻⁵. That is about 100× today's areal density at about 1000× the pressure.
  - Source: [Ochs et al. 2025/2026, arXiv:2511.10885](https://arxiv.org/html/2511.10885v1)
- [mainstream] **Ahmad et al. 2026.** At T_e = T_i there is no net energy. At T_e = 0.5 T_i, TB gives net energy for T_i = 190–330 keV and NS gives none. At T_e = 0.25 T_i, TB gives 125–500 keV and NS 140–500 keV.
  - Minimum n_eτT: 1.289e22 m⁻³ s keV at 245 keV (no radiation, TB); 1.171e23 at 250 keV (T_e = 0.5 T_i, TB); 1.503e22 at 270 keV (T_e = 0.25 T_i, TB).
  - Some of the paper's numbers are internally inconsistent.
  - Source: [Ahmad et al. 2026](https://www.frontiersin.org/journals/nuclear-engineering/articles/10.3389/fnuen.2026.1714531/full)
- [mainstream] **Rider 1995 thesis**, "Fundamental limitations on plasma fusion systems not in thermodynamic equilibrium", MIT EECS, advisor L. M. Lidsky — [DSpace@MIT](https://dspace.mit.edu/handle/1721.1/11412). The Phys. Plasmas companion is 4, 1039 (1997), doi:10.1063/1.872556.
- [mainstream] **Liu et al. 2024, re-running Rider's case** (T_i = 300 keV, ⟨E_e⟩ = 35 keV, T_e = 23.33 keV):
  - With Rider's ⟨σv⟩: P_ie/P_fus = 36 and P_brem/P_fus = 0.5.
  - With the 2016 ⟨σv⟩: P_ie/P_fus = 19.5 and P_brem/P_fus = 0.27.
  - Thermal crossing at T_e = 136.1 keV, where P_fus − P_brem = 8.31% of P_fus.
  - ENN scenario (T_i = 150 keV, T_i/T_e = 4): P_brem/P_fus = 0.80 and P_ie/P_fus = 11.
  - The recirculating power is "substantially greater than the fusion power", which confirms Rider.
  - Source: [Liu et al. 2024](https://arxiv.org/html/2405.13260v1)
- [contested — no source on the page] Wikipedia says "some calculations" give bremsstrahlung at least 1.74× fusion for p–¹¹B and 1.39× for ³He–³He. It also says p–¹¹B needs 500× the nTτ of D–T and has about 2100× lower power density at equal pressure — [Wikipedia: Aneutronic fusion](https://en.wikipedia.org/wiki/Aneutronic_fusion)

### Inferences
**Prefactor checks** (computed):
- The Putvinski coefficient 8.511 = (16/3)√(2π/3)·(2√3/π), the non-relativistic result with the Gaunt factor.
- W0 = α r_e² m_e c³ n_e² = 1.422e-36 n_e² W m³ in SI, which is e⁶n_e²/(ħ m_e c²) in Gaussian units.
- The low-t limit is 5.355e-37 Z_eff n_e² T_keV^{1/2} W m⁻³. This matches the NRL Formulary 1.69e-32 W cm³ eV^-1/2 (= 5.344e-37) to 0.2%; the NRL value is recalled, because the PDF returned 403.
- Rider's 5.172e-37 is 3.4% lower than this. I could not trace why.

**Relativistic enhancement** of Putvinski over the non-relativistic formula at Z_eff = 3 (computed): ×1.16 at T_e = 50 keV, ×1.63 at 150 keV, ×2.56 at 300 keV. The Putvinski and Rider forms agree within 1.4% from 50 to 300 keV at Z_eff = 3.

**Synchrotron**, computed from Larmor with ⟨v⊥²⟩ = 2kT/m:
`P_sync = e⁴B²n_e kT_e/(3π ε₀ m_e³c³) = 6.21e-17 B²[T] n_e T_e[keV]` W m⁻³, which matches NRL's 6.21e-28 in cgs-eV units. This is the optically thin upper bound.
- For p–¹¹B at T_i = 300 keV, T_e = 150 keV, f_B = 0.15: P_sync/P_br = 41 at β = 1, 82 at β = 0.5 and 410 at β = 0.1.
- Even at β = 1, n_e = 1e20 needs B = 3.7 T.
- So synchrotron loss must be cut by more than 100× by opacity and wall reflection (Trubnikov-type treatment), or killed by very high β (FRC). This is a hard requirement for magnetic p–¹¹B.

**P_fus/P_br** (computed; Putvinski Eq. 7 bremsstrahlung, density cancels, f_B = n_B/n_i = 0.15, n_e = n_p + 5n_B), shown as NS / TB:

| T_i, keV | T_e/T_i = 1 | 0.5 | 0.33 | 0.25 |
|---|---|---|---|---|
| 150 | 0.36 / 0.42 | 0.65 / 0.77 | 0.88 / 1.03 | 1.05 / 1.23 |
| 200 | 0.40 / 0.49 | 0.80 / 0.97 | 1.10 / 1.34 | 1.34 / 1.62 |
| 300 | 0.34 / 0.44 | 0.76 / 0.99 | 1.12 / 1.45 | 1.40 / 1.81 |
| 500 | 0.19 / 0.28 | 0.47 / 0.71 | 0.76 / 1.14 | 1.00 / 1.51 |

- At T_e = T_i = 300 keV the best f_B is 0.105, giving P_fus/P_br = 0.457 (TB).
- With TB and Putvinski's +10% kinetic enhancement, T_i = 300 keV and T_e = 150 keV give a ratio of about 1.08. That fits "marginal ignition".
- **Conflict:** Ochs et al. report fusion above bremsstrahlung for T_e = T_i above 150 keV at f_b = 0.1. My computation gives at most about 0.5, a factor of about 2 apart. Possible causes: a different definition of f_b, the kinetic enhancement factor, or the Svensson fit. I could not resolve it.
- Liu et al.'s Rider case was not reproduced: I get P_br/P_fus = 0.30 (NS) and 0.23 (TB) at f_B = 0.1, against their 0.5 and 0.27. Their composition is not stated on the page.

**Physical limit.** At T_e ≈ T_i, thermal p–¹¹B cannot ignite on bremsstrahlung alone. Every positive balance depends on keeping T_e well below T_i. Rider's point is that the power P_ie flowing from ions to electrons is then 10–36× the fusion power, so it must be recirculated at about 90%+ efficiency. That is why bremsstrahlung recovery (Section 5) carries the whole argument.

### Gaps
- The Svensson 1982 total-bremsstrahlung fit coefficients were not obtained.
- I did not obtain a Trubnikov synchrotron formula with wall reflectivity, or Putvinski's Fig. 4 numbers.
- Rider's thesis full text and exact tables were not read. The "1.74×" figure is unsourced on Wikipedia.
- No "Ostrovsky et al. 2025" paper was found.

## 4. Laser-driven p–¹¹B (HB11, Margarone 2020, Bonvalet), avalanche claims and critiques

### Takeaway
The best laser shots give about 10⁹–10¹⁰ α/sr. That is a fusion energy out per laser energy in of about 10⁻⁵–10⁻⁴ (computed, assuming isotropic emission), roughly 10⁶–10⁷ short of the target gain of about 100 needed for a power plant. Analyses, including one supported by HB11 itself, find the yields consistent with beam fusion and higher cross sections, not an "avalanche".

### Cited Findings
- [mainstream] **Margarone et al. 2020** (Front. Phys. 8:343), LFEX at ILE Osaka:
  - Laser: 1.4 kJ, 2.6 ps, ~50 µm spot, ~2–3×10¹⁹ W/cm².
  - Pitcher–catcher target: 25 µm CH, then 2 mm BN.
  - Protons above 10¹³ sr⁻¹ (1–10 MeV), with about 7% laser-to-proton conversion.
  - Alphas: **5×10⁹ sr⁻¹ (±20%)** at 8–10 MeV. Results are described as preliminary.
  - Source: [Margarone et al. 2020](https://www.frontiersin.org/articles/10.3389/fphy.2020.00343/full)
- [mainstream] **Yields by experiment:**
  - Belyaev 2005: about 10³ per pulse, later corrected to 10⁵ α/sr.
  - Labaune 2013 (LULI, pitcher–catcher): 10⁷ α/sr.
  - Giuffrida 2020 (PALS, BN in-target, PRE 101 013204): 10¹⁰ α/sr/shot.
  - PALS 2024, 550 J, 350 ps, ~10¹⁶ W/cm², B₁₈H₂₂ at 0.9 g/cm³: 1.7×10⁹ α/sr/shot.
  - Bonvalet 2021 (PRE 103, 053202) is listed as a pure-B pitcher–catcher experiment, but its yield is not given in the text.
  - Source: [arXiv:2407.06729](https://ar5iv.labs.arxiv.org/html/2407.06729)
- [mainstream] Picciotto 2014 reached about 10⁹/sr. Giuffrida 2020 estimated about 10¹¹ α per pulse in total — [search summary of PRX/PRE records](https://doaj.org/article/053da7de879a42b082e2832a5291f0fd)
- [contested — preprint] **Wei et al. 2023**, boron-doped foam heated by hohlraum X-rays plus a TNSA proton beam: up to 10¹⁰ α/sr per shot, "the highest yield value normalized to the laser energy". The yield per proton "exceeds the classical expectation of beam-target reactions by up to four orders of magnitude" — [arXiv:2308.10878](https://arxiv.org/abs/2308.10878)
- [contested→refuted for present experiments] **Avalanche claims** (Hora et al. 2017, Laser Part. Beams 35, 730, "Road map to clean energy using laser beam ignition of boron-hydrogen fusion").
  - Belloni et al. 2018 (Phys. Plasmas 25, 020701) find an avalanche time of almost 1 µs on the 675 keV resonance in PALS plasmas, ten times longer than claimed. They conclude no avalanche can develop there, even with idealised confinement. This is an abstract-level summary — [search summary](https://www.frontiersin.org/researchtopic/11751)
  - Belloni 2021 (PPCF 63 055020) models the α–p elastic chain at n_e = 10²⁴–10²⁸ cm⁻³, T_e ≲ 100 keV. He finds k∞ "of the order of 10⁻²" even at the most favourable parameters, where a value near 1 would be needed — [arXiv:2012.14533](https://arxiv.org/abs/2012.14533)
- [mainstream] **Mehlhorn et al.** (poster supported by HB11 Energy; co-authors include Margarone, Batani and Belloni):
  - The yields are "consistent with their generation via beam fusion reactions". "Beam fusion reactions do not scale to net gain".
  - "Higher cross section appears to explains higher yields, not avalanche".
  - Thick-target yield for 1 MeV protons in BN at STP is 6.56e-5. "Gains >= 100 are required for economic power production".
  - Plant example: 280 MJ per pulse at 5 Hz from a 2 MJ laser; gain 140; 1400 MW fusion → 560 MWe at 40%; 140 MW recirculating (7% laser); 420 MW net.
  - Source: [Mehlhorn et al., HB11 poster PDF](https://hb11.energy/wp-content/uploads/2022/01/Study-of-Proton-Boron-Fusion-Burn-Driven-by-Short-Pulse-Lasers-by-T-Mehlhorn.pdf)
- [contested] Wikipedia records that HB11's Margarone et al. 2022 was claimed as the first commercial demonstration of fusion, with no figures — [Wikipedia: Aneutronic fusion](https://en.wikipedia.org/wiki/Aneutronic_fusion)

### Inferences
Fusion energy out per laser energy in, computed as `E_fus/E_L = (4π N_sr/3) × 8.68 MeV / E_L`. This assumes isotropic emission, which overstates beam-target yields that are forward-peaked:
- Margarone 2020: 2.1×10⁻⁵.
- PALS B₁₈H₂₂: 1.8×10⁻⁵.
- Giuffrida 2020: 7.7×10⁻⁵, taking 10¹¹ α total and an assumed ~600 J.

Beam fusion ceiling: 6.56e-5 reactions per 1 MeV proton gives E_fus/E_p = 5.7×10⁻⁴. Even a 100%-efficient proton driver cannot break even by beam–target fusion, so the gain must come from thermonuclear burn of a compressed, hot fuel, where Ochs et al.'s ρR ≈ 100 g/cm² applies.

Mehlhorn plant check, computed: Q_eng = 0.4 × 140 × 0.07 = 3.92 (they publish 560/140 = 4). Q_eng = 1 needs a gain of 35.7. Their "3.525 mg ¹¹B burned per 280 MJ pulse" recomputes to 3.68 mg, a 4% difference.

### Gaps
- I did not obtain Bonvalet et al. 2021's yield numbers, the Hora et al. avalanche numbers, or HB11's current (2024–26) claims.
- The LFEX and PALS laser energies for Giuffrida 2020 are not in the fetched text.

## 5. Direct energy conversion: venetian blind, TWDEC, inductive, photon; thermodynamics

### Takeaway
Measured electrostatic direct conversion reached 48–65%. The 70–90% figures for venetian-blind, periodic-focusing, TWDEC and ICC converters are design projections. Helion's 95% round-trip figure is a company-reported subscale circuit result. Direct conversion is not Carnot-limited because the ions are a directed, low-entropy flux. What limits it is the spread of ion energies against a finite number of collector stages, plus secondary electrons, space charge and losses in the power conditioning.

### Cited Findings
- [mainstream] **Venetian blind:**
  - About 60% "economically compatible", technical maximum up to 70%, for 100–150 keV D–T plasma (Barr, Burleigh, Dexter, Moir & Smith 1974).
  - Beam experiments: 65% peak, 50% minimum (Moir et al. 1977; Barr, Moir & Hamilton 1982).
  - TMX plasma converter: 48% (Barr & Moir 1983).
  - Post 1969: efficiency depends on the number of collectors.
  - Rosenbluth & Hinton 1994 questioned charge neutrality over a Debye length of "a few centimeters".
  - Source: [Wikipedia: Direct energy conversion](https://en.wikipedia.org/wiki/Direct_energy_conversion)
- [mainstream] **Barr & Moir 1983** tested single-stage and two-stage venetian-blind converters up to 100 kV and 6 kW, for about 80 h in total — [Fusion Sci. Technol., ANS](https://wx1.ans.org/pubs/journals/fst/a_20820)
- [contested — projection] **Periodic electrostatic focusing:** about 60% economic, up to 90% technical, at about 600 keV (Smith et al. 1972) — [Wikipedia: Direct energy conversion](https://en.wikipedia.org/wiki/Direct_energy_conversion)
- [contested — projection] **TWDEC** for 14.7 MeV D–³He protons:
  - ARTEMIS design: 75% (Momota et al. 1992). Parameters: 10 m cavity, 10 T, 155 MHz.
  - Maximum projected 90% (Rostoker, Binderbauer & Monkhorst 1997).
  - Source: [Wikipedia: Direct energy conversion](https://en.wikipedia.org/wiki/Direct_energy_conversion)
  - Momota, Miley & Nadler 2000 estimate about 60% overall for a TWDEC on a D–³He IEC core — [search summary, NIFS/ANS records](https://www.nifs.ac.jp/report/nifs641.html)
  - Kobe University experiments (Shibata et al. 2019) report only that "Higher efficiency was obtained with dual-frequency modulation", with no numbers — [Shibata et al., PFR 14 2405027](https://www.jspf.or.jp/PFR/PFR_articles/pfr2019/pfr2019_14-2405027.html)
  - NIFS-250 (Iiyoshi, Momota et al. 1993) gives no efficiency figures in its abstract — [NIFS-250](https://www.nifs.ac.jp/report/nifs250.html)
- [contested — projection] **Inverse cyclotron converter (TAE):** 90% maximum projected, at 5 MHz and 0.6 T — [Wikipedia: Direct energy conversion](https://en.wikipedia.org/wiki/Direct_energy_conversion)
- [mainstream] **Magnetic direct conversion cycle:** 62% from alpha energy to work (Oliphant, Ribe & Coultas 1973, Nucl. Fusion 13, 529) — [Wikipedia: Direct energy conversion](https://en.wikipedia.org/wiki/Direct_energy_conversion)
- [contested — company-reported, subscale] **Helion:**
  - "recover energy at over 95% round-trip efficiency for over 1 million pulses" (2015) — [Helion: About](https://www.helionenergy.com/about)
  - Net electricity is possible "If 95% of the input energy after the pulse is recovered". With high-β FRCs, "for every unit of energy applied by the magnetic field, we get nearly one unit of energy in the plasma itself". Pulses last under 1 ms. No formula or component efficiencies are given — [Helion blog](https://www.helionenergy.com/blog/how-to-make-fusion-electricity-without-ignition)
- [contested] **X-ray conversion:** Lawrenceville Plasma Physics claims 81% photoelectric conversion (US patent 7482607, 2009). Quimby 1977 says "hundreds or thousands of layers are needed" to absorb the X-rays — [Wikipedia: Direct energy conversion](https://en.wikipedia.org/wiki/Direct_energy_conversion); [Wikipedia: Aneutronic fusion](https://en.wikipedia.org/wiki/Aneutronic_fusion)

### Inferences
**Ideal N-stage electrostatic collector**, computed. Ion flux `f(E) ∝ (E−s)e^{−(E−s)}` (E in units of T) for E > s, where s is the sheath or ambipolar drop. Each ion lands on the highest collector at or below its energy and returns qV. The collector voltages are optimised by search.

| s | N=1 | N=2 | N=4 | N=8 |
|---|---|---|---|---|
| 0 | 0.419 | 0.593 | 0.741 | 0.824 |
| 3 | 0.642 | 0.778 | 0.874 | 0.928 |
| 5 | 0.731 | 0.840 | 0.911 | 0.951 |

The analytic check for N = 1, s = 0: V* = (1+√5)/2 and η = 0.420. The measured 48% (one stage) and 65% are what this ideal predicts once secondary electrons, imperfect transmission and grid interception are added. The **limit is the beam's energy spread relative to N, not temperature**: a monoenergetic beam with one matched collector reaches η → 1 in principle.

**Why direct conversion is not Carnot-limited.** Carnot bounds heat engines working between two reservoirs. A directed ion beam with a narrow energy spread carries little entropy per joule, so it can be decelerated against an electric field, ideally with no rejected heat. What limits it instead:
1. The energy spread across N stages, as in the table.
2. Space charge and neutralisation.
3. Secondary electron emission and sputtering.
4. Interception on the electrode grids.
5. Losses in the DC→AC inverter and power conditioning.
6. For TWDEC, bunching efficiency (the phase spread).

**Bremsstrahlung is the exception.** It is broadband X-rays, so converting it either goes through heat, which is Carnot-limited, or through layered photoelectric converters, which the 81% claim relies on.

**Inductive recovery** is circuit physics: the round-trip efficiency is set by coil resistance, switches and capacitor loss. A 95% round trip on the *driver* energy is plausible for a circuit, but it is not the same as converting the *fusion* energy at 95%.

### Gaps
- Moir & Barr 1973 (Nucl. Fusion 13, 35) and Post 1969 were not read directly. Their numbers come from Wikipedia's summary.
- I found no measured TWDEC efficiency, and no measured conversion efficiency for an X-ray converter.
- Helion has published no peer-reviewed data on recovering fusion-product energy.

## 6. Overall plant efficiency arithmetic: recirculating power, Q_eng vs Q_sci

### Takeaway
For magnetic confinement, `Q_eng = η_c η_d (Q_sci + 1)`. For inertial confinement, `Q_eng = η_c η_d G`. With a steam cycle (η_c ≈ 0.4), magnetic fusion needs Q_sci ≈ 2.6–4 just to break even. NIF's record Q_sci = 4.13 is still only about 1% of engineering breakeven at its wall-plug efficiency. High-efficiency driver-energy recovery (Helion's argument) lowers the Q_sci needed into the 0.06–0.5 range.

### Cited Findings
- [mainstream] **Engineering breakeven:** `Q_E = P_fus/P_heat = 1/(η_heat · f_recirc · η_elec · (1 − f_ch))`, with f_ch = 0.2 for D–T.
  - Typical values: η_heat = 0.6–0.7, η_elec = 0.35–0.4, conversion loop 0.20–0.25.
  - Practical example: η_heat = 0.7, η_elec = 0.4, f_recirc = 0.2 gives Q ≈ 22 (Entler 2015; McCracken & Stott 2005).
  - Magnetic confinement typically needs Q of 5–8; inertial confinement needs 50–100.
  - Source: [Wikipedia: Fusion energy gain factor](https://en.wikipedia.org/wiki/Fusion_energy_gain_factor)
- [mainstream] **NIF:**
  - 5 Dec 2022: 3.15 MJ from 2.05 MJ on target, Q = 1.54, with "over 400 MJ" of electrical energy drawn.
  - 2025: 8.6 MJ from 2.08 MJ, Q = 4.13. NIF has reached ignition seven times.
  - Source: [Wikipedia: Fusion energy gain factor](https://en.wikipedia.org/wiki/Fusion_energy_gain_factor)
- [mainstream] **Definitions:** Q_eng = power delivered to the grid / recirculating electric power used to heat, sustain, control and assemble the plasma. Q_sci = fusion power / applied heating power — [Helion blog](https://www.helionenergy.com/blog/how-to-make-fusion-electricity-without-ignition)
- [mainstream] **HB11-type laser plant arithmetic:** gain 140, 7% laser, 40% conversion, 4:1 gross-to-recirculating, 420 MW net — [Mehlhorn et al.](https://hb11.energy/wp-content/uploads/2022/01/Study-of-Proton-Boron-Fusion-Burn-Driven-by-Short-Pulse-Lasers-by-T-Mehlhorn.pdf)

### Inferences
Formulas, all computed:
- **Magnetic, steady state.** Heating power P_h is absorbed. Output = η_c(P_fus + P_h); input = P_h/η_d + P_aux. With P_aux = 0, Q_eng = η_c η_d (Q_sci + 1), so breakeven needs Q_sci = 1/(η_c η_d) − 1.
  - η_c = 0.4, η_d = 0.7 → Q_sci = 2.57.
  - η_c = 0.4, η_d = 0.5 → Q_sci = 4.0.
  - η_c = 0.6, η_d = 0.7 → Q_sci = 1.38.
  - η_c = 0.8, η_d = 0.9 → Q_sci = 0.39.
  - The Entler formula with the Wikipedia inputs gives 22.3 (published about 22).
- **Inertial** (driver energy not recovered): Q_eng = η_c η_d G.
  - NIF 2022: yield/wall-plug = 3.15/400 = 0.79%. At η_c = 0.4, Q_eng = 0.0032.
  - NIF 2025: 0.0086, assuming the same ~400 MJ wall-plug draw (an assumption).
  - HB11 example: 3.92.
- **Pulsed with inductive recovery** (Helion-type). E_out = η_rec E_in + η_dc f_ch Q_sci E_in, so breakeven needs Q_sci = (1 − η_rec)/(η_dc f_ch). With η_dc = 0.9 and f_ch = 0.95, η_rec = 0.95 / 0.90 / 0.80 / 0.60 gives Q_sci = 0.058 / 0.117 / 0.234 / 0.468. The claim depends on η_rec over the *whole* cycle, including plasma formation and losses, which has not been shown publicly [contested].
- **p–¹¹B with bremsstrahlung.** If a fraction f_x of the output is X-rays converted thermally at 0.4, and the charged particles are converted directly at 0.8, then η_c = 0.8 − 0.4 f_x. That is 0.60 at f_x = 0.5 and 0.48 at f_x = 0.8. Since Section 3 gives P_br ≳ P_fus, f_x ≥ 0.5 is realistic, which removes most of direct conversion's advantage over a steam cycle.

### Gaps
- I did not obtain an authoritative wall-plug energy for NIF 2025.
- I found no peer-reviewed full-cycle Q_eng analysis for Helion's Polaris.
- The ARTEMIS overall plant efficiency number (Momota 1992) was not read in the original.
