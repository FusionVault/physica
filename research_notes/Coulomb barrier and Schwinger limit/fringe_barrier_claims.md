# Fringe, contested and refuted claims of beating the Coulomb barrier or reaching the Schwinger regime, with checks against mainstream bounds

How to read these notes:
- **Status labels.** **[mainstream]** agrees with standard nuclear physics and QED and has been reproduced. **[contested]** is disputed in the literature, with no settled refutation. **[refuted]** has a named refutation, a failed replication or a retraction, or our own calculation shows it is internally inconsistent by many orders of magnitude (stated each time).
- **"(computed)".** A number recomputed in `python3 -I` from the authors' own figures, with CODATA 2018 constants. The script is `/tmp/claude-1000/-home-dan-code-FusionVault-quantum-gravimeter/7beb35b3-c723-440c-9450-66a991641392/scratchpad/calc/barrier_checks.py`, and its asserts pass. Computed values are the note author's inferences, not published claims.
- **Recalled inputs.** Three inputs were not re-fetched this session and are marked RECALLED where used:
  - the Bosch–Hale S(0) values 53.7 and 55.6 keV b for the n and p branches of D–D;
  - the neutron fluence-to-dose coefficient, about 4×10⁻¹⁰ Sv cm² at 2.5 MeV;
  - the Pd density (12.02 g/cm³) and molar mass (106.42 g/mol).

  Q-values used: 3.269 MeV (d(d,n)³He), 4.033 MeV (d(d,p)t) and 23.847 MeV (d(d,γ)⁴He).
- **"(snippet)".** The fact comes from a search-engine summary of the page, not from reading the page itself. Treat these as weaker.
- **Builds on** `research_notes/Gravity beams fringe physics and DIRDs/baker_li_hfgw_and_pais_patents.md` (Q4 there) and does not repeat it. That note already covers:
  - the Pais S_max formula and its 10²⁴–10³³ W m⁻² reproductions;
  - the $508k NAWCAD test and its "Pais Effect could not be proven" result;
  - the craft patent's eq. 4 (ν ≈ 4.2×10²¹ Hz).

## Q1. Pais's "Plasma compression fusion device" (US 2019/0295733) and his "vacuum polarization" equations

### Takeaway
The fusion application was published on 26 Sep 2019 and never granted. The USPTO mailed two non-final actions in 2020, and the application was abandoned on 14 Jul 2021 for failure to respond. No test of the device is on record.

Its numbers do not survive a check (computed):
- **Charge.** One coulomb on a 0.3 m fusor gives a 4×10¹¹ V/m surface field (400× Pais's own vacuum-breakdown figure) and 3×10¹⁰ J of stored energy.
- **Field.** Spinning that charge at 10⁴ rad/s makes about 4 mT, not the claimed 10⁶ T.
- **Field energy.** 10⁶ T in a 0.3 m sphere would hold 5.6×10¹⁵ J (1.3 Mt TNT).
- **Fusion power.** At the patent's own density (10²⁰ m⁻³) and D–D temperature (20 keV), D–D fusion yields about 160 W, not GW–TW.

The "vacuum polarization" side fails by a similar margin. Even at the 5.6×10¹⁵ V/m implied by Pais's own craft-patent example, the Euler–Heisenberg correction is 9×10⁻¹⁰ of the field energy, and the Schwinger pair rate is 10⁻²⁶⁸ m⁻³ s⁻¹ (computed). Status: **[refuted]** by calculation and by the absence of any test. The Navy test of the related field generator was null (prior notes).

### Cited Findings
- **Bibliographic record.** Application US 15/928,703, "Plasma Compression Fusion Device". Inventor Salvatore Pais; assignee US Department of the Navy; filed and priority date 22 Mar 2018; published 26 Sep 2019. Google Patents gives the legal status "Abandoned", with these events:
  - non-final action mailed 16 Aug 2020;
  - non-final action mailed 11 Dec 2020;
  - abandoned 14 Jul 2021 "for failure to respond to an office action".

  Source: [US20190295733A1](https://patents.google.com/patent/US20190295733A1/en)
- **Conflicting press.** Atomic Insights (Rod Adams, 13 Oct 2019) headlined that the USPTO "has issued a patent" ([Atomic Insights](https://atomicinsights.com/?p=17049)). A commenter on that page notes it is a pre-grant publication, not an issued patent, and that publication does not require a working device ([Atomic Insights comments](https://atomicinsights.com/?p=17049)). Futurism, citing The Drive, said in 2019 that it "sailed through" without rejection ([Futurism](https://futurism.com/the-byte/navy-scientist-ufo-tech-fusion), snippet). Both statements predate the 2020 non-final actions shown in the Google Patents file history.
- **Equations in the description** (as extracted by Google Patents; the OCR garbles layout). Source: [US20190295733A1](https://patents.google.com/patent/US20190295733A1/en)
  - Eq. 1 (Lawson): nTτ_E ≥ 3×10²¹ keV s m⁻³.
  - Eq. 2: energy gain ~ B³.
  - Eq. 3: fusion power density ~ P²·B⁴.
  - Eqs. 4–5: the two D–D branches, ²H+²H → ³H (1.01 MeV) + p (3.02 MeV) and ³He (0.82 MeV) + n (2.45 MeV), each 50%.
  - Eqs. 6–7: the S_max = f_G(σ²/ε0)[…] forms copied from the HFGW patent.
  - Eq. 8: E_T ≈ m R_v² Ω² exp(2Ωt).
  - Eq. 9: S_max ≈ (Q²/ε0)(R_v²/R_s⁵) Ω exp(2Ωt).
  - Eq. 10: B/R_v ~ curl B.
  - Eq. 11: B_MAX ≈ μ0 σ R_ω ω³ t_op².
- **Stated numbers.** Source: [US20190295733A1](https://patents.google.com/patent/US20190295733A1/en)
  - plasma density about 10²⁰ m⁻³;
  - D–T above 175 million °C and D–D about 232 million °C; p–¹¹B about 2 billion °C;
  - "possibly exceeding 30 Tesla" needed, with B_MAX of order 10⁶ T claimed achievable;
  - spin about 10⁴ rad/s, based on a lab test of 10 cm discs at 100,000 RPM;
  - energy gain of order 10¹⁸ "derived from Equation 2";
  - output of gigawatts to terawatts from kW–MW input;
  - vacuum about 10⁻⁵ Torr; device diameter 0.3–2 m; 10–15 cm B₄C or W shielding;
  - "at least one Coulomb" per dynamic fusor.

  The claims recite no output figure, and the document presents no experimental data.
- **Terms not used.** "Schwinger", "vacuum polarization" and "Pais effect" do not appear in the fusion application. — [US20190295733A1](https://patents.google.com/patent/US20190295733A1/en)
- **Peer-reviewed paper.** Pais published the same concept as "The Plasma Compression Fusion Device—Enabling Nuclear Fusion Ignition", IEEE Transactions on Plasma Science, DOI 10.1109/TPS.2019.2942997 (received 11 Apr 2018; accepted 11 Jul 2019). Affiliation: NAVAIR/NAWCAD. It is a design paper; no experimental results were found. — [IEEE Xplore](https://ieeexplore.ieee.org/abstract/document/8871349/authors) (snippet)
- **Critique.** Rod Adams reached no verdict of his own. The quantitative objections come from commenters on his page ([Atomic Insights](https://atomicinsights.com/?p=17049)):
  - one coulomb on a picofarad-scale object would need gigavolts and would arc over;
  - toroidal coils have almost no external field;
  - roughly 75 MW of neutron power in a 1 GW(th) D–D machine would activate and melt the components;
  - GW–TW output needs a far larger heat-transfer area.
- **Craft-patent formulas** (already in the prior notes, Q4):
  - ε0E_S² = ħν⁴/(8π²c³) gives ν ≈ 4.2×10²¹ Hz;
  - the S_max surface field 5.6×10¹⁵ V/m sits 240× below the Schwinger field 1.32×10¹⁸ V/m;
  - NAWCAD spent $508k with a null result.

  Source: [prior notes](../Gravity%20beams%20fringe%20physics%20and%20DIRDs/baker_li_hfgw_and_pais_patents.md)
- **Mainstream benchmark: vacuum birefringence.** PVLAS (2020 final report) did not reach the QED vacuum birefringence. It measured Δn = (12±17)×10⁻²³ at 2.5 T against Δn(QED) = 2.5×10⁻²³; the uncertainty is about 7× the prediction. — [arXiv:2005.12913](https://arxiv.org/abs/2005.12913v1) (snippet) **[mainstream]**
- **Mainstream benchmark: laser intensity.** The record laser intensity is 1.1(±0.2)×10²³ W/cm² (CoReLS 4 PW, f/1.1, 1.1 µm FWHM spot; Yoon et al., Optica 8, 630–635, 2021). — [IBS record](https://pr.ibs.re.kr/handle/8788114/10375) (snippet) **[mainstream]**
- **Mainstream benchmark: SLAC E144** (Burke et al., PRL 79, 1626, 1997). It collided 46.6 GeV electrons with terawatt 527 nm pulses and saw 106 ± 14 positrons above background via nonlinear Breit–Wheeler pair production. The laser parameter was ξ ≈ 0.4. — [SLAC-PUB-7564](https://www.slac.stanford.edu/pubs/slacpubs/7500/slac-pub-7564.pdf) (snippet); [LUXE LoI](https://arxiv.org/pdf/1909.00860) (snippet) **[mainstream]**

### Inferences
**Fusion application (computed unless noted)**
- **Charge.** 1 C on a sphere of radius 0.15 m gives E_surf = 4.0×10¹¹ V/m, V = 6.0×10¹⁰ V and U = 3.0×10¹⁰ J (about 7 t TNT). At r = 1 m it gives 9.0×10⁹ V/m and 4.5×10⁹ J. Even the largest fusor is 9× above the 10⁹ V/m vacuum limit Pais states in US 10,135,366 (prior notes). **[refuted]**
- **Field from spinning charge.** A uniformly charged shell (1 C, R = 0.15 m, ω = 10⁴ rad/s) has an interior field B = (2/3)μ0Qω/(4πR) = 4.4×10⁻³ T. The claimed 10⁶ T is 2.2×10⁸ times larger.
  - Pais's eq. 11 is dimensionally μ0σv·(ωt)². It reaches 10⁶ T only through the unbounded (ωt)² factor: about 1.2 s of operation at 10⁴ rad/s (our arithmetic).
  - Nothing supplies the energy. A 10⁶ T field has a pressure of 4.0×10¹⁷ Pa (3.9×10¹² atm) and holds 5.6×10¹⁵ J (1.3 Mt TNT) in a 0.3 m sphere, against kW–MW input. **[refuted]**
- **Gain and power density.** "Gain ~ B³" and "power density ~ β²B⁴" look like tokamak scalings at fixed β and size (our inference; not sourced). Turning a proportionality into the number 10¹⁸ requires G = (B/1 T)³, which has no basis. At the 30 T the patent itself cites, the same formula gives 2.7×10⁴.
- **Lawson criterion.** 232 million °C is 20.0 keV. At n = 10²⁰ m⁻³ the patent's own eq. 1 needs τ_E ≥ 1.5 s, and 3×10²¹ keV s m⁻³ is a D–T figure; D–D needs more.
  - Confinement pressure is not the obstacle. 2nkT = 6.4×10⁵ Pa needs only 1.3 T at β = 1.
- **Power.** Treat the D–D reactivity at 20 keV with constant RECALLED S-factors: ⟨σv⟩ ≈ 3.9×10⁻²⁴ m³ s⁻¹. At 10²⁰ m⁻³ that gives 1.1×10⁴ W m⁻³, or 160 W in a 0.3 m sphere: 6×10⁶ to 6×10⁹ below "GW–TW".
  - Producing 1 GW in that sphere would need n ≈ 2.5×10²³ m⁻³ (1.6 GPa, 63 T at β = 1), held for seconds.
  - A 1 GW D–D source would also emit about 8.5×10²⁰ neutrons/s.
- **Neutron and gamma counts:** see Q2 (8.5×10¹¹ n s⁻¹ per watt).

**"Vacuum polarization" and Schwinger regime (computed)**
- **Benchmark scales.**
  - E_S = m²c³/(eħ) = 1.3233×10¹⁸ V/m; B_S = 4.414×10⁹ T.
  - I_S = 2.32×10²⁹ W/cm² (cε0E²/2) or 4.65×10²⁹ (cε0E²).
  - The record laser (1.1×10²³ W/cm²) reaches E0 = 9.1×10¹⁴ V/m = 6.9×10⁻⁴ E_S, short by 2.1×10⁶ in intensity.
- **Euler–Heisenberg.** The fractional correction to the field energy is ΔL/L0 = (α/45π)(E/E_S)². Our A_e = 2α²ħ³/(45μ0m⁴c⁵) = 1.325×10⁻²⁴ T⁻² gives Δn = 3A_eB² = 2.48×10⁻²³ at 2.5 T, which matches PVLAS's 2.5×10⁻²³. The correction is:
  - 3×10⁻²³ at the 10⁹ V/m breakdown limit;
  - 2.4×10⁻¹¹ at the record laser field;
  - 9.4×10⁻¹⁰ at Pais's 5.6×10¹⁵ V/m (σ = 5×10⁴ C/m²);
  - 2.7×10⁻¹² for 10⁶ T.

  Any "vacuum engineering" effect is therefore a part in 10⁹ or less, even granting fields 10⁶ times above breakdown.
- **Schwinger pair rate.** Leading term w = (eE)²/(4π³ħ²c)·exp(−πE_S/E). log10 w (m⁻³ s⁻¹) is:
  - −1.8×10⁹ at 10⁹ V/m;
  - −1933 at the record laser field (ignoring the fact that a plane wave makes no pairs);
  - −268 at 5.6×10¹⁵ V/m;
  - +40 at 0.1 E_S;
  - +55 at E_S.

  The process switches on only within a factor of about 10 of E_S. Pais's devices sit 10⁹ (breakdown-limited) to 240 (his own σ) below it. **[refuted]**
- **The ħν⁴/(8π²c³) relation.** It is the standard zero-point energy density up to an angular cutoff ν **[mainstream identity]**.
  - At the Compton cutoff ω_C = mc²/ħ it equals 1.80×10²² J/m³.
  - ε0E_S²/2 = 7.75×10²⁴ J/m³ is exactly (π/α) = 430.5 times that (computed; an identity, since E_S and ω_C share m, c and ħ).
  - Equating the two therefore just returns ν = 5.4 ω_C (4.2×10²¹ s⁻¹). It is dimensional numerology, not a dynamical condition for pair creation or "polarization".

### Gaps
- I did not read the contents of the 2020 USPTO non-final actions, so the examiner's grounds are unknown. Patent Center was not queried.
- I did not read the full text of the IEEE TPS paper (paywalled); I relied on the abstract only.
- I found no NAWCAD or Navy statement specific to the fusion device, and no record that it was ever built.
- The constant-S D–D reactivity is an approximation. Bosch–Hale fits were not fetched, and the true ⟨σv⟩ at 20 keV may differ by a factor of about 2.

## Q2. Cold fusion / LENR: Fleischmann–Pons, the DOE reviews, Google 2019, ARPA-E, EU programmes, SPAWAR, Widom–Larsen, Rossi, Mills

### Takeaway
Fleischmann–Pons (F–P) claimed 1.75 W of excess heat from a 0.4×10 cm Pd rod together with about 4×10⁴ neutrons/s. If the heat were D–D fusion, it would come with 1.5×10¹² n/s (computed). The neutrons are therefore missing by 3.7×10⁷, and F–P said so themselves ("10¹¹–10¹⁴" reactions needed versus "1–2×10⁴" measured).

Both DOE reviews (1989 and 2004) rejected a dedicated programme. In 2004 the 18 reviewers split about evenly on excess heat, and two-thirds found the nuclear evidence inconclusive.

Recent tests:
- Google's 2015–2018 programme (about 30 researchers, $10M) found no anomalies.
- For ARPA-E's $9.96M, 8-project topic (2023) I found no published confirmation of excess heat. Its visible output is beam-driven fusion-rate work.
- DARPA's 2026 MARRS programme sets its state of the art at about 0.1 reactions s⁻¹ g⁻¹, which is 6×10⁻¹⁴ W/g (computed).

Widom–Larsen fails its own rate arithmetic by about 300×. The Rossi Lugano COP was inflated by an emissivity error (corrected COP 1.13–1.32). Hydrino states are not square-integrable.

### Cited Findings
**Fleischmann & Pons 1989** (J. Electroanal. Chem. 261, 301–308; errata in vol. 263, pp. 187–188). Read from the New Energy Times copy, whose margin annotations are Krivit's, not the authors'. Source for all items below: [F–P 1989 PDF](https://www.newenergytimes.com/v2/library/1989/1989Fleischmann-PrelimNote-SK-MARKUP.pdf)
- **Table 1, excess heating, 0.4×10 cm rod:**

  | Current density | Excess heat | Specific rate |
  |---|---|---|
  | 8 mA cm⁻² | 0.153 W | 0.122 W cm⁻³ |
  | 64 mA cm⁻² | 1.751 W | 1.39 W cm⁻³ |
  | 512 mA cm⁻² | 26.8 W | 21.4 W cm⁻³ |

  The 512 mA cm⁻² row was "measured on electrodes of length 1.25 cm and rescaled to 10 cm".
- **Heat claims.** Enthalpy "can exceed 10 W cm⁻³" for more than 120 h. More than 4 MJ cm⁻³ was liberated, which the authors call "inconceivable" except by nuclear processes. A 1×1×1 cm cube partly melted (melting point 1554 °C) and destroyed part of the fume cupboard.
- **Neutrons.** The flux "is of the order 4 × 10⁴ s⁻¹ for a 0.4 cm × 10 cm rod electrode polarised at 64 mA cm⁻²". The raw count was "3 times above" background, from a Bonner-sphere monitor with an estimated efficiency of about 2.4×10⁻⁴, further reduced by about 100 for geometry, counted over 50 h.
- **Tritium.** Reaction (v) d(d,p)t runs "to the extent of 1–2 × 10⁴ atoms s⁻¹", whereas "the data on enthalpy generation would require rates … in the range 10¹¹–10¹⁴ atoms s⁻¹". The authors conclude "other nuclear processes must be involved".
- **Rate per pair.** They ask whether a fusion rate of 10⁻¹⁹ s⁻¹ per D–D pair is possible at about 1 eV.
- **Gamma line.** The erratum changes the capture-gamma line in reaction (vii) from "2.5 MeV" to "2.225 MeV". It relabels the Fig. 1A y-axis from 0–25,000 counts to 0–1000 counts, and adds Marvin Hawkins as co-author.
- **Wikipedia's figure conflicts.** Wikipedia gives the F–P neutron claim as "about 4,000 neutrons per second". That conflicts with the primary text's 4×10⁴. — [Wikipedia: Cold fusion](https://en.wikipedia.org/wiki/Cold_fusion)

**The "missing neutrons" and "missing gamma" arguments**
- **Missing neutrons.** With the known branching ratios, 1 W would mean about 10¹² neutrons/s, a level that "would have been fatal to the researchers". — [Wikipedia: Cold fusion](https://en.wikipedia.org/wiki/Cold_fusion)
- **Missing gamma.** The DOE 2004 report says D–D fusion produces p, n, t, ³He, ⁴He and high-energy gammas "in well known proportions". The ⁴He + γ channel occurs "approximately only once for every 10⁷ D + D fusion reactions". Proponents' claim that the 24 MeV appears as lattice heat without gammas drew "serious concerns". — [DOE 2004 report](https://lenr-canr.org/acrobat/DOEreportofth.pdf)
- **Petrasso et al. (MIT), Nature 339, 183–185 (1989).** The Utah gamma peak was narrower than the detector resolution and lacked the Compton edge. The neutron rate F–P derived from it was too high by a factor of about 40. — [AIP press release, May 1989](https://newenergytimes.com/v2/sr/ColdFusion1989/APS-Baltimore/19890515-AIP-Press-Release.shtml) (snippet); citation via [Nature 344, 401 (1990)](https://www.nature.com/articles/344401a0) (snippet)
- **Koonin & Nauenberg, Nature 339, 690 (1989).** The D₂ molecular d+d fusion rate is 3×10⁻⁶⁴ s⁻¹; electron-mass enhancements of 5–10 would be needed to reach the claimed rates. — [K–N preprint copy](https://newenergytimes.com/v2/library/1989/1989KooninNauenberg.shtml) (snippet)
- **Metzler, Hunt & Galvanetto (arXiv:2208.07245; New J. Phys. 26, 101202, 2024).** Observable D–D fusion in solids would need ">40 orders of magnitude" of enhancement. — [arXiv:2208.07245](https://arxiv.org/abs/2208.07245)

**DOE reviews**
- **1989 ERAB panel.** Co-chaired by John Huizenga and Norman Ramsey. DOE's response (James Decker to ERAB chairman John Landis, 22 Mar 1990) accepted the report's two thrusts:
  - "skepticism regarding the scientific validity of cold fusion", with "unresolved scientific issues";
  - "no need to establish special cold fusion programs", with research "at a modest level of effort" through normal competitive funding.

  Source: [DOE/ERAB transmittal, OSTI via Black Vault](https://documents.theblackvault.com/documents/doe/osti-coldfusion/ERAB.pdf)
- **Panel size and date.** The panel was described as 22 members, with its final report delivered in November 1989. — [Morrison newsletter / search summary](https://newenergytimes.com/v2/archives/DROM/17.shtml) (snippet)
- **2004 review ("Report of the Review of Low Energy Nuclear Reactions", 1 Dec 2004).** Source for all items below: [DOE 2004 report](https://lenr-canr.org/acrobat/DOEreportofth.pdf)
  - **Process.** Requested by Hagelstein, McKubre, Nagel, Chubb and Hekman. Nine mail reviewers plus nine panel reviewers met on 23 Aug 2004, for 18 comments in total.
  - **Excess heat.** Reviewers were "split approximately evenly".
  - **Repeatability.** Most said the effects "are not repeatable", that the magnitude "has not increased in over a decade", and that many experiments were poorly documented.
  - **Helium.** ⁴He was reported in 5 of 16 heat-producing cells, "very close to … background". Air contamination was raised.
  - **Nuclear evidence.** "Two-thirds of the reviewers commenting on Charge Element 1 did not feel the evidence was conclusive for low energy nuclear reactions, one found the evidence convincing, and the remainder … somewhat convinced."
  - **Recommendation.** No reviewer recommended a focused federal programme; they supported individual peer-reviewed proposals.
  - **Conclusions.** They were "similar to those found in the 1989 review".
  - **Beam experiments.** Low-energy beam experiments on deuterated foils (the screening work, Q4) were "generally viewed … favorably" but called "peripheral".

**Google-funded programme, Berlinguette et al., "Revisiting the cold case of cold fusion", Nature 570, 45–52 (2019)**
- **Scale.** About 30 researchers, led by Matthew Trevithick (Google) and Curtis Berlinguette (UBC), with $10M from Google. — [Physics World, 28 May 2019](https://physicsworld.com/a/cold-case-of-cold-fusion-reopened-by-google-funded-scientists)
  - About 12 projects ran from 2015 to 2018 and produced 28 peer-reviewed papers. — [IEEE Spectrum](https://spectrum.ieee.org/google-nuclear-energy)
- **Loading.** Normal Pd takes about 0.7 H per metal atom, and the claimed threshold is 0.875. The team exceeded the threshold "only twice" in "several dozen tests". — [Physics World](https://physicsworld.com/a/cold-case-of-cold-fusion-reopened-by-google-funded-scientists)
- **Heat.** None of 420 nickel-powder samples "provided evidence of excess heat". — [Physics World](https://physicsworld.com/a/cold-case-of-cold-fusion-reopened-by-google-funded-scientists)
- **Tritium.** Berkeley pulsed-plasma runs (Schenkel), lasting weeks, saw no enhanced tritium. — [Physics World](https://physicsworld.com/a/cold-case-of-cold-fusion-reopened-by-google-funded-scientists)
- **Conclusion.** They "found no evidence of anomalous effects claimed by proponents of cold fusion". — [Physics World](https://physicsworld.com/a/cold-case-of-cold-fusion-reopened-by-google-funded-scientists). They "did not find credible evidence of anomalies". — [IEEE Spectrum](https://spectrum.ieee.org/google-nuclear-energy)
- **Volume and pages** as cited by a later paper. — [scientific.net citation](https://www.scientific.net/KEM.881.51) (snippet)

**ARPA-E LENR Exploratory Topic (announced Feb 2023)**
- **Awards** ([New Energy Times listing, 1 Mar 2023](https://news.newenergytimes.net/2023/03/01/u-s-department-of-energys-arpa-e-funds-low-energy-nuclear-reaction-research/)):

  | Awardee | Amount | Focus |
  |---|---|---|
  | Stanford | $1,500,000 | LENR-active nanoparticles in D₂ gas |
  | MIT | $2,000,000 | Testing nuclear-anomaly claims in gas-loaded metal–H |
  | Amphionic | $295,924 | — |
  | Energetics Technology Center | $1,500,000 | Pd–D co-deposition |
  | LBNL | $1,500,000 | Excitation energies below 500 eV |
  | Texas Tech | $1,150,000 | — |
  | University of Michigan | $1,108,412 | Pd–Ni nanocomposite heat claims |
  | University of Michigan | $902,213 | n, γ and ion emission measurements |

  The total is $9,956,549. Director Evelyn Wang framed the question: "does this area show promise … Or can we conclusively show that it does not?"
- **Annual report.** The FY2023 Annual Report (Aug 2025) lists 8 projects and $10.0M. — [ARPA-E FY2023 report](https://arpa-e.energy.gov/sites/default/files/2025-09/ARPA-E%20FY%202023%20Annual%20Report.pdf) (snippet)
- **Visible outputs.**
  - The LBNL/UC Davis project "Quantifying Nuclear Reactions in Metal Hydrides at Low Energies" ($1.5M over 2.5 years) — [LBNL ATAP](https://atap.lbl.gov/?p=1493) (snippet) — led to Karahadian, Colborne, Persaud, Schenkel & Munday, arXiv:2512.06212 (Dec 2025). D–D fusion below 2.5 keV in Pd and Ti hydrides is "enhanced by over 10¹⁸ compared to theoretical bare-nucleus fusion rates", but "the reaction rates reported here are low". — [arXiv:2512.06212](https://arxiv.org/abs/2512.06212v1) **[contested]** (not peer-reviewed)
  - Texas Tech (Gillespie, Lin, Jones, Duncan) showed that untreated stainless steel 316 holds 15–240 pmol of ⁴He (0.07–1.20 pmol/cm²), so ⁴He as a nuclear-ash marker needs pretreatment. Its ARPA-E funding is not stated on the page. — [arXiv:2405.02307](https://arxiv.org/abs/2405.02307)
- **No final results.** I found no ARPA-E final programme report and no published confirmation of excess heat from any of the eight teams (see Gaps).
- **DARPA MARRS (follow-on).** "Mechanisms for Amplification of fusion Reaction Rates in Solids". Programme manager Thomas Schenkel (DSO); proposers' day 20 Jan 2026. Source for all items below: [DARPA MARRS proposers' day](https://www.darpa.mil/sites/default/files/attachment/2026-01/marrs-program-proposers-day-presentation.pdf)
  - **Objective.** Raise rates "from ~0.1 (state-of-the-art) to 10⁶ fusion reactions/second/gram": 10³ s⁻¹ g⁻¹ at the end of Phase 1 (18 months) and 10⁶ at the end of Phase 2 (18 months).
  - **Detection limits on the slides.** About 10⁹ reactions s⁻¹ g⁻¹ by heat, about 0.1 s⁻¹ by neutrons and about 10⁻³ s⁻¹ by protons. The slides also state "1 Watt~10¹² reactions/second".
  - **Screening design targets.** "Ue≥600 eV", and a theory slide predicting "≥1 Watt/gram for 10⁻¹³ barn in 'properly loaded metal'". The measured Pd yield "near room temperature corresponds to a cross section of ~2×10⁻¹³ barn".
  - It cites ARPA-E LENR, Chen et al. (Nature 2025), Czerski (PRC 106, L011601, 2022) and Metzler et al. (2024).

**EU programmes**
- **CleanHME** ("Clean Energy from Hydrogen-Metal Systems", grant 951974, H2020 FET Proactive). Coordinator University of Szczecin; 1 Aug 2020 to 31 Jan 2025; total €5,623,722.50, EU contribution €5,509,447.50; 17 participants in 9 countries. CORDIS lists no results. — [CORDIS 951974](https://cordis.europa.eu/project/id/951974)
- **CleanHME claims.** The project says "first experiments have shown that the total heat energy produced exceeds by many orders of magnitude the chemical energy". — [Polish Science](https://www.polishscience.pl/en/?p=141620) (snippet)
- **HERMES.** Reported as a €4M award from the same FET round. — [Science|Business](https://sciencebusiness.net/news/projects-touting-potential-new-hydrogen-energy-source-among-winners-latest-eu-research-awards) (snippet)

**Navy SPAWAR CR-39 tracks**
- **Claim.** Mosier-Boss, Szpak, Gordon & Forsley, "Triple tracks in CR-39 as the result of Pd–D co-deposition: evidence of energetic neutrons", Naturwissenschaften 96, 135–142 (2009). Triple tracks attributed to ¹²C(n,n′)3α would imply neutrons of at least 9.6 MeV. — [paper copy](https://newenergytimes.com/v2/library/2008/2008BossTripleTracks.pdf) (snippet); [Wikipedia: Pamela Mosier-Boss](https://en.wikipedia.org/wiki/Pamela_Mosier-Boss)
- **Kowalski's challenge.** Kowalski challenged the 2007 CR-39 interpretation, and the authors replied in EPJ Applied Physics. A 2010 paper reported Pd/D tracks "indistinguishable" from those of a DT neutron source. — [EPJAP reply](https://resolve.cambridge.org/core/journals/the-european-physical-journal-applied-physics/article/reply-to-comment-on-the-use-of-cr39-in-pdd-codeposition-experiments-a-response-to-kowalski/5D8D5610FB89C14454595ABEDD10712E) (snippet); [EPJAP 2010 abstract](https://www.cambridge.org/core/journals/the-european-physical-journal-applied-physics/article/abs/comparison-of-pdd-codeposition-and-dt-neutron-generated-triple-tracks-observed-in-cr39-detectors/D7A063731E77E4D92ADC997A5A2AA2AA) (snippet)
- **Validation.** Wikipedia says such claims cannot be validated without quantitative neutron analysis. — [Wikipedia: Cold fusion](https://en.wikipedia.org/wiki/Cold_fusion) **[contested]**

**Widom–Larsen theory**
- **Claim.** Widom & Larsen, "Ultra low momentum neutron catalyzed nuclear reactions on metallic hydride surfaces", EPJC 46, 107 (2006), arXiv:cond-mat/0505026. "Heavy" surface electrons (mass shifted up by local fields) are captured by protons, e + p → n + ν; "No Coulomb barriers exist for the weak interaction neutron production." — [arXiv:cond-mat/0505026](https://arxiv.org/abs/cond-mat/0505026)
- **Critique.** Ciuchi, Maiani, Polosa, Riquer, Ruocco & Vignati, EPJC 72, 2193 (2012), arXiv:1209.6501. Source for all items below: [arXiv:1209.6501](https://arxiv.org/abs/1209.6501)
  - W–L's rate is Γ ≈ 1.86×10⁻³(β−β0)² Hz ≈ 0.56 Hz for β = m*/m_e ≈ 20, with β0 ≈ 2.59.
  - The correct calculation gives 1.8×10⁻³ Hz at β = 20 and 6.9×10⁻⁷ Hz at β = 5.2. W–L's value would need β = 61 (m* ≈ 30 MeV).
  - W–L's own estimator gives β = 1.01 for Ωp = 0.1 eV and 1.8 for their Ωp = 0.8 eV, "below threshold for nuclear transmutation".
  - Conclusion: a rate "about 300 times smaller", and it is "questionable" that such β can be realised.

**Rossi E-Cat, Lugano report (Levi et al., Oct 2014)** Source for the items below: [Lugano report](https://www.mech.kth.se/~hanno/Lugano_report.pdf)
- **Run.** 32 days in March 2014, at about 1260 °C then about 1400 °C.
- **COP.** "About 3.2 and 3.6"; per file, e.g. COP = (2128.32 + 307.98)/(815.86 − 37.77) = 3.13 ± 8%.
- **Energy.** Net about 1.5 MWh (abstract) or 1,618,194 Wh = 5825 MJ (eq. 28). With a 1 g charge that is 5.8×10⁶ MJ/kg.
- **Isotopes.**
  - ⁷Li fell from a natural 93% to 7.9% (SIMS) or 42.5% (ICP-MS).
  - ⁶²Ni rose from 3.6% to 98.7% (SIMS) or 99.3% (ICP-MS).
- **Radiation.** No radiation signs were found.
- **Temperature method.** Temperature came from an IR camera with alumina emissivity taken from the literature.
- **Custody.** Rossi intervened during fuel insertion, startup, shutdown and spent-fuel extraction. Dorigo called this "a major flaw, which totally invalidates any conclusions". — [Wikipedia: Energy Catalyzer](https://en.wikipedia.org/wiki/Energy_Catalyzer)
- **Replica critique.** Valat, Goldwater, Greenyer, Higgins & Hunt (Martin Fleischmann Memorial Project, itself a pro-LENR group), J. Condensed Matter Nucl. Sci. 21, 81–96 (2016), used an IR-camera-calibrated replica. — [MFMP paper](https://jcmns.org/article/72418-investigations-of-the-lugano-hotcat-reactor/attachment/150447.pdf)
  - File 4: excess 102.9 W, COP 1.13, actual surface temperature about 944 °C (reported 1257 °C).
  - File 8: excess 289.8 W, COP 1.32, about 1045 °C (reported 1404 °C).
  - Conclusion: "at least an order of magnitude lower excess heat than reported". The emissivity function used was wrong, and the dummy-run power appears to have been in error by 251%.
- **2013 test critique.** The earlier 2013 test was criticised as "far from being an independent one"; resistive heating alone reproduces its temperature curve. — [arXiv:1306.6364](https://arxiv.org/abs/1306.6364) (snippet)
- **Industrial Heat lawsuit.** Industrial Heat said it could not reproduce the results after three years. Rossi sued in April 2016 over an $89M fee, and the case settled on 5 Jul 2017 with undisclosed terms. No independent or peer-reviewed tests exist. — [Wikipedia: Energy Catalyzer](https://en.wikipedia.org/wiki/Energy_Catalyzer)

**Mills hydrino**
- **Claimed states.** Binding energies W_q = W1/q² with W1 = 13.6 eV and q a fraction.
- **Rathke's critique,** New J. Phys. 7, 127 (2005):
  - the n < 1 Schrödinger solutions "are not square integrable";
  - the stability bound for the non-relativistic hydrogen atom is about 20 eV at most;
  - Mills's wave equation is not Lorentz invariant;
  - his radial solution f(r) = δ(r − r_n)/r is not a solution;
  - there is "no theoretical support of the hydrino hypothesis".

  Source: [Rathke, ar5iv](https://ar5iv.labs.arxiv.org/html/quant-ph/0505150) **[refuted]**
- **SunCell power claims.** A 2020 report associated with the company (Booker) claims 300 kW and 200 kW. — [peak-oil.org summary](https://news.peak-oil.org/2020/01/validation-report-of-300000w-of-power-produced-by-the-suncell-brilliant-light-power/) (snippet). I found no independent peer-reviewed validation. **[contested]**

### Inferences
**Missing products (computed)**
- **Per watt.** Mean D–D Q = 3.651 MeV, so 1 W = 1.710×10¹² reactions/s = 8.55×10¹¹ n/s.
- **F–P's 1.751 W.** That implies 1.50×10¹² n/s against the 4×10⁴ claimed: a deficit of 3.7×10⁷. For 26.8 W: 2.3×10¹³ n/s.
- **Dose.** The 1.751 W case gives a flux of 1.2×10⁷ n cm⁻² s⁻¹ at 1 m. With the RECALLED dose coefficient that is about 17 Sv/h, a lethal dose within roughly 15–30 minutes. **[refuted]** for D–D as the heat source.
- **Gamma budget.** If the heat came from d+d → ⁴He, it needs 2.6×10¹¹ ⁴He/s per W. In the standard channel each would carry a 23.8 MeV gamma. Proponents must posit a branch with no gamma and no other channels, which every measured d+d energy contradicts (DOE 2004). **[contested]** at best for ⁴He-as-ash.
- **Energy per atom.** F–P's 4 MJ cm⁻³ is 367 eV per Pd atom: about 100× chemistry, but still "only" 10⁻⁴ of a 23.8 MeV reaction per Pd atom.
- **Rate per pair.** Assuming D/Pd ≈ 1 in the rod, the heat needs 3.5×10⁻¹¹ s⁻¹ per D pair, which is 1.2×10⁵³ times Koonin–Nauenberg's 3×10⁻⁶⁴ s⁻¹. F–P's neutrons alone need 9×10⁻¹⁹ s⁻¹ (they quote 10⁻¹⁹), which is 10⁴⁵ above K–N. This is consistent with Metzler's ">40 orders".

**Programme bounds (computed)**
- **MARRS rates in watts.** The SOA 0.1 s⁻¹ g⁻¹ is 5.8×10⁻¹⁴ W/g. The Phase 2 goal of 10⁶ s⁻¹ g⁻¹ is 5.8×10⁻⁷ W/g. The heat-detection floor of 10⁹ s⁻¹ g⁻¹ is 5.8×10⁻⁴ W/g.
- **F–P against the SOA.** F–P's 1.39 W cm⁻³ (0.116 W/g, 2.0×10¹¹ s⁻¹ g⁻¹) is 2×10¹² times the rate the mainstream LENR-adjacent programme now treats as state of the art. Even MARRS's Phase 2 goal would be 10³ below calorimetric visibility.
- **Programme trajectory.** The ARPA-E and DARPA work has moved from testing heat claims to beam-assisted fusion-rate measurement. That trajectory is itself evidence about the heat claims (inference).

**Theory and claim checks**
- **Widom–Larsen** (computed). Their eq. 5 at β = 20 gives 0.564 Hz, which reproduces; Ciuchi's 1.8×10⁻³ Hz is 313× lower. The threshold (m_n − m_p − m_e)c² = 0.782 MeV gives β0 = (m_n − m_p)/m_e = 2.531; Ciuchi quote 2.59, a small definitional difference. **[refuted]** as a quantitative mechanism.
- **Lugano energy** (computed). 5825 MJ from 1 g is 1.3×10⁵ times petrol (about 46 MJ/kg; RECALLED). p + ⁷Li → 2α (17.35 MeV) on the 0.011 g of ⁷Li would release only 2.6 GJ, 45% of the claimed total, even at complete burn-up.
- **Lugano gammas.** p + ⁵⁸Ni → ⁵⁹Cu, the route to ⁶²Ni enrichment, would leave positron emitters and 511 keV gammas, yet no radiation was detected (our inference; ⁵⁹Cu decay data not fetched).
- **Lugano correction.** The MFMP correction lowers the File-4 net output from 1658 W to 103 W (16×). Status: COP 3.2–3.6 **[refuted]**; residual excess **[contested]**, and claimed only by a pro-LENR group.
- **Hydrino energies** (computed from Mills's W_q formula). q = 1/2 binds at 54.4 eV and q = 1/4 at 217.6 eV, releasing 40.8 and 204 eV. Both exceed Rathke's quoted ~20 eV stability bound.

### Gaps
- The Nature 2019 paper itself was not read (the publisher redirected to a login), so the maximum D/Pd, calorimeter sensitivity and neutron limits come from Physics World, a secondary source.
- The ERAB 1989 report body (DOE/S-0073) was not text-extractable: the copy read is the DOE transmittal and response letter. I have no ERAB numbers for neutron or tritium upper limits.
- **ARPA-E outcomes:**
  - no ARPA-E final report or project-closeout page was found;
  - the ARPA-E topic page returned 404;
  - MIT, Stanford, Michigan, ETC and Amphionic outcomes were not found;
  - DARPA cites "Karahadian et al., https://doi.org/10.2172/2569797 (2025)", but OSTI refused the connection.
- HERMES's grant number, dates and results were not verified, and no CleanHME results or deliverables were found.
- No independent SPAWAR CR-39 replication, numbers (track counts) or programme end date were found.
- Petrasso's Nature paper and the K–N paper were seen only through secondary copies and snippets.
- The neutron dose coefficient, petrol energy density and ⁵⁹Cu decay data were RECALLED, not fetched.

## Q3. Bubble fusion / sonofusion (Taleyarkhan)

### Takeaway
Taleyarkhan et al. (Science 2002) claimed 4–8×10⁴ n/s of 2.5 MeV neutrons and tritium from cavitated chilled deuterated acetone. Their own tritium implied 7×10⁵ n/s, 9–17× more than they measured, whereas D–D predicts about 1:1. The claimed power is about 10⁻⁷ W (computed).

ORNL's Shapira and Saltmarsh found no correlated neutrons, and an upper limit at least 10⁴ below what the tritium needs. Naranjo showed the 2006 spectrum matched ²⁵²Cf. Purdue (2008) found misconduct on two counts: false claims of independent confirmation, and an added co-author. ONR debarred Taleyarkhan for 28 months. Status: **[refuted]** (no independent replication; misconduct findings did not rule on validity).

### Cited Findings
- **Science 295, 1868–1873 (8 Mar 2002).** Source for all items below: [NET copy of Taleyarkhan et al. 2002](https://www.newenergytimes.com/v2/library/2002/2002-TAL-EvidenceForNuclearEmissions.shtml)
  - **Neutrons.** About 4–8×10⁴ n/s of D–D neutrons, from a ~4% rise in 2.5 MeV counts: ~4 SD for one run, >10 SD aggregated.
  - **Tritium.** +7.1 cpm (7 h runs, ~4 SD collectively) and +14.6 cpm (12 h runs, >4.5 SD), with 1 SD ≈ 3.5 cpm.
  - **Coincidence.** A coincidence peak within 0–2 µs of the sonoluminescence flash, seen only in chilled (~0 °C) C₃D₆O.
  - **Conditions.** Drive at 19.3 kHz; onset about −7 bar, amplitude about ±15 bar. Nucleation was seeded by 14 MeV neutrons from a pulsed generator at ~10⁶ n/s. Simulated temperatures were 10⁶–10⁷ K.
  - **Tritium–neutron mismatch.** The tritium-implied D–D neutron rate is about 7×10⁵ n/s, roughly 9–17× the direct measurement.
  - **Caveat in the paper itself.** Independent attempts "yielded smaller neutron emission and no evidence for real coincidences".
- **Shapira & Saltmarsh, PRL 89, 104302 (2002).** No 2.5 MeV neutron emission correlated with SL; any emission is "at least 4 orders of magnitude too small" to explain the tritium. — [arXiv:1307.3217 (a pro-Taleyarkhan comment summarising it)](https://arxiv.org/pdf/1307.3217) (snippet); [Wikipedia: Bubble fusion](https://en.wikipedia.org/wiki/Bubble_fusion)
- **Later history.** Source for all items below: [Wikipedia: Bubble fusion](https://en.wikipedia.org/wiki/Bubble_fusion)
  - Galonsky (Science 2002) argued the neutrons were too energetic for D–D.
  - The 2004 PRE paper claimed replication, with neutron emission "many acoustic cycles later".
  - Xu & Butt (Nucl. Eng. Des., 2005) were former students and did not attempt coincidence measurements.
  - The 2006 PRL claimed about 5×10³–10⁴ n/s.
  - Naranjo found that spectrum statistically inconsistent with D–D and consistent with ²⁵²Cf.
  - Putterman and Suslick found no evidence of fusion (BBC Horizon, 2005).
  - The 2002 patent application was rejected in 2005.
  - Purdue first cleared Taleyarkhan (7 Feb 2007), then reopened the case after Rep. Brad Miller's report. On 18 Jul 2008 a five-institution committee substantiated 2 of 12 allegations, without ruling on "the validity of the experimental results".
  - He lost his named professorship on 27 Aug 2008, and ONR debarred him for 28 months (to Sept 2011).

### Inferences
- **Power** (computed). 4×10⁴, 8×10⁴ and 7×10⁵ n/s correspond to D–D powers of 4.7×10⁻⁸, 9.4×10⁻⁸ and 8.2×10⁻⁷ W. Even if true, this is a neutron-source claim, not an energy source.
- **Ratio** (computed). The tritium/neutron ratio is 8.8–17.5 where D–D gives about 1. The two signatures cannot both be D–D products at the stated efficiencies. That internal mismatch matches the Shapira–Saltmarsh ≥10⁴ bound.
- **Background.** The seeding source (14 MeV neutrons at 10⁶ n/s) is about 10–25× brighter than the claimed signal. A 4% excess in a pulse-height window therefore depends on how the detector separates 14 MeV recoils, which Galonsky's and Naranjo's critiques target (inference).

### Gaps
- The Shapira–Saltmarsh PRL and Naranjo comment were not read in full, and their numeric upper limits (n/s) were not extracted.
- The Purdue committee report and ONR debarment notice were not read; I relied on Wikipedia's summary.
- I found no post-2011 bubble-fusion replication attempts.

## Q4. Enhanced electron screening in metals: measured U_e against theory, and its status

### Takeaway
Beam–target d(d,p)t measurements at keV energies in deuterated metals consistently give screening energies U_e of about 200–700 eV. Examples: Ta 322, Zr 297, Pd 313, Cu 470 and Pt 670 eV. A D₂ gas target gives 15 ± 5 eV, close to the 20–30 eV adiabatic limit.

The data are real and reproduced by several groups **[mainstream]**. The interpretation is open **[contested]**:
- the Debye model's "~300 eV" reproduces numerically but was excluded by Huke et al. (2008), who blamed deuteron-density dynamics and surface contamination;
- the Szczecin group now invokes a ⁴He 0⁺ threshold resonance, defects (340 eV in ZrD₂), and "thermal fusion in ion tracks";
- the 2024 review says no suitable theory exists.

Screening is a multiplicative boost that matters only near and below a few keV: it is 1.5× at 5 keV for U_e = 300 eV (computed). It does not make room-temperature fusion visible. With U_e = 300 eV it gives σ ≈ 10⁻¹⁸ b at 300 K in a naive model, and that would already exceed the measured cold rates.

### Cited Findings
- **Review slides on electron screening** (Edinburgh, ERC talk, no date on the slides; cites 2024 papers). Source for all items below: [Edinburgh slides](https://indico.ph.ed.ac.uk/event/373/contributions/3866/attachments/1856/2896/electron-screening.pdf)
  - **Gas target.** D₂ gas: "Ue = 15 ± 5 eV in fair agreement with adiabatic limit (Ue = 20–30 eV)" (Greife et al., Z. Phys. A 351, 107, 1995).
  - **Metals** (Raiola et al., EPJA 19, 283, 2004; 55 samples): Cu 470 eV and Pt 670 eV ("factor ~20–25 higher"); Nd and Hf < 30 eV.
  - **Pattern.** Elements in the same group show similar U_e: a "large effect ~300 eV (ca. 10x higher than in D₂ gas)" and a "small effect ~30 eV".
  - **Debye model.** R_D = 69√(T/(n_eff ρ_a)) m, and "R_D ~ 1/10 R_a ⇒ U_e,D ~ 300 eV for d+d in metals".
  - **Lithium.** ⁶,⁷Li(p,α) results: U_e = 3790 ± 330, 1280 ± 60 and 185 ± 150 eV (Cruz et al., PLB 624, 181, 2005). The extraction does not show which value goes with which target. The slides call it "additional (Debye) enhancement" in Li metal and PdLi.
  - **Puzzle.** "Despite theoretical developments over past 30 years, discrepancies remain … ('electron screening puzzle')". It is called a puzzle for many reactions (Spitaleri et al., PLB 755, 275, 2016). Candidate explanations include stopping powers extrapolated from high energy, channelling, target type, and dynamic or clustering effects.
- **Sample count conflict.** A search summary gives 58 deuterated samples for Raiola 2004 and says insulators and semiconductors show a small effect. — [search summary](https://www.research.ed.ac.uk/en/publications/enhanced-electron-screening-in-ddpt-for-deuterated-metals/) (snippet). This conflicts with the slides' "55 samples".
- **Huke, Czerski, Heide, Ruprecht, Targosz & Żebrowski, PRC 78, 015803 (2008), arXiv:0805.4538.** Source for all items below: [arXiv:0805.4538](https://arxiv.org/abs/0805.4538)
  - **Table I (eV):**

    | Host | U_e (eV) |
    |---|---|
    | Ta | 322 ± 15 |
    | Zr | 297 ± 8 |
    | Al | 190 ± 15 |
    | Sr | 350–800 |
    | Na | ≲150 |
    | Pd | 313 ± 2 |

  - **Conclusion.** The previously reported "unusual changes in the reaction yields are mainly due to deuteron density dynamics". Surface contamination from residual vacuum gas also matters. The Debye–Hückel model "could be clearly excluded".
  - **Overview figure.** Fig. 11 compiles Garching 1990, Tohoku 1998 and 2002 (including 30 nm PdO) and Bochum 2002–2004, with values up to about 800 eV.
- **Czerski, Huke, Heide & Ruprecht, EPL 68, 363 (2004).** Experiment against dielectric-function theory:

  | Host | Experiment (eV) | Theory (eV) |
  |---|---|---|
  | Zr | 295 ± 7 | ~112 |
  | Pd (PdD₀.₂) | 296 ± 15 | ~134 |
  | Ta | 302 ± 13 | ~136 |

  Source: [search summary of the EPL table](https://www.newenergytimes.com/v2/library/2004/2004Czerski-The2Hdp3HReaction.pdf) (snippet)
- **Czerski et al., EPL 113, 22001 (2016).** A ZrD₂ measurement below 25 keV found the rise in enhancement "much weaker than that determined before". It can be explained by screening plus a suggested 0⁺ threshold resonance in ⁴He, with surface impurities raising U_e. — [search summary](https://avesis.ogu.edu.tr/yayin/31f1baf4-bb50-48ec-9cd0-be4b167f564d/screening-and-resonance-enhancements-of-the-h-2d-ph-3-reaction-yield-in-metallic-environments) (snippet)
- **Czerski, Dubey, Kowalska, Haridas Das, Kaczmarski, Targosz-Ślęczka & Valat, arXiv:2409.02112** (v1 Aug 2024; v2 Dec 2025; no journal reference). Source for all items below: [arXiv:2409.02112](https://arxiv.org/abs/2409.02112v2)
  - ZrD₂ measured down to E_cm = 675 eV in ultra-high vacuum.
  - The thick-target yield falls about 7 orders of magnitude, then plateaus below about 2.5 keV.
  - U_e = 340 eV, against "about 100 eV" for defect-free material.
  - They claim "first direct observation" of thermal D–D fusion in ion tracks. **[contested]**
- **Review status, Vesić, "Electron Screening in Laboratory Nuclear Reactions", Particles 7, 834–855 (2024).** Values in several cases are "more than an order of magnitude above" the adiabatic prediction, and "a suitable theory has not yet been found" (wording from the group's conference abstract). — [preprints.org](https://preprints.org/manuscript/202409.0669/v1) (snippet); [IJS news](https://f2.ijs.si/en/news/2024101410394302/two-interesting-review-articles-recently-published-in-mdpi-journals) (snippet)
- **DOE 2004 view.** Reviewers "generally viewed … favorably" these beam/foil screening studies but called them peripheral to cold fusion. — [DOE 2004 report](https://lenr-canr.org/acrobat/DOEreportofth.pdf)

### Inferences
**Enhancement factors** (computed; d+d Gamow energy E_G = 985.8 keV, 2πη = 31.40/√E[keV]). Using f = σ(E+U_e)/σ(E):

| E_cm | U_e = 25 eV | U_e = 300 eV | U_e = 600 eV |
|---|---|---|---|
| 10 keV | 1.012 | 1.16 | 1.33 |
| 5 keV | 1.036 | 1.50 | 2.17 |
| 2.5 keV | 1.10 | 2.99 | 7.6 |
| 1 keV | 1.47 | 47 | 717 |

So the "anomaly" is a factor of 1.2–3 at the 2.5–10 keV energies where it is measured. The large factors appear only where the bare yield is itself tiny.

**Debye check** (computed). T = 300 K, ρ_a = 6×10²⁸ m⁻³ and n_eff = 1 give R_D = 4.9×10⁻¹² m and U_D = 295 eV, reproducing the slides' "~300 eV". R_D = 0.05 Å is far below any electron screening length in a degenerate metal. Classical Debye screening of conduction electrons is not physically applicable (our inference), consistent with Huke's exclusion. **[refuted]** as an explanation.

**Room-temperature extrapolation** (naive thermal model: σ = (S/E)·exp(−√(E_G/(E+U_e))), pair energy kT = 25.9 meV, D/Pd = 1, relative speed 2.2×10³ m/s, RECALLED S = 109 keV b):

| U_e | σ at 300 K |
|---|---|
| 0 eV | 10⁻²⁶⁷⁵ b |
| 27 eV | 10⁻⁷⁶ b |
| 100 eV | 10⁻³⁶·⁵ b |
| 300 eV | 10⁻¹⁸·³ b |
| 600 eV | 10⁻¹¹ b |

- **MARRS check.** The model reproduces MARRS's "≥1 W/g for 10⁻¹³ barn" (we get 2.5 W/g).
- **U_e needed in this model:**
  - 168 eV for MARRS's measured SOA (0.1 s⁻¹ g⁻¹);
  - 270 eV for the MARRS Phase 2 goal (10⁶ s⁻¹ g⁻¹);
  - 423 eV for F–P's heat;
  - 464 eV for 1 W/g;
  - 228 eV for F–P's neutrons.
- **Inference.** If the 300–700 eV measured at keV energies applied unchanged to thermal deuterons, Pd–D would already produce about 2×10⁷ (300 eV) to 9×10¹⁵ (700 eV) reactions s⁻¹ g⁻¹ (computed), and neutron counters would see it. They do not (SOA 0.1 s⁻¹ g⁻¹). So keV-energy U_e values cannot be extrapolated to thermal energies as a constant shift, or the naive model is wrong (localized deuterons at about 2.9 Å spacing are not a 300 K gas). Either way, claims that the screening anomaly "explains" cold-fusion heat are **[refuted]** by the neutron-rate bound.

**"10¹⁸" enhancement check** (computed). How large an enhancement over bare is depends on E. 10¹⁸ needs U_e of only 51 eV at E_cm = 50 eV, 195 eV at 100 eV and 990 eV at 200 eV. The headline figure in arXiv:2512.06212 is mainly a statement about how small bare sub-keV cross-sections are. Whether it needs physics beyond the measured U_e depends on the energies of the plateau events, which I could not extract.

### Gaps
- The Raiola 2002 and 2004 tables (all 55–58 hosts) and the Kasagi (Tohoku) PdO and Pd values were not extracted beyond the slide and figure overview.
- The Vesić 2024 full text was blocked (MDPI returned 403). Its compiled table and conclusions were seen only via snippets.
- The Greife 1995 gas value (15 ± 5 eV per the slides) was not checked against the original. Some secondary sources, and the task framing, quote about 25 eV for the adiabatic limit; the slides give 20–30 eV.
- I did not read the Czerski PRC 106, L011601 (2022) threshold-resonance parameters.
- The naive thermal model is the note author's construction. A proper lattice calculation (Koonin–Nauenberg style, with localized deuterons) was not done.

## Q5. Claims from 2020–2026: "breaking the Coulomb barrier" at low energy, or Schwinger-limit effects on a tabletop

### Takeaway
The 2020–2026 record moved toward mainstream beam–target physics:
- UBC and LBNL (Nature 2025) showed that electrochemical loading raises beam-driven D–D rates in Pd by 15(2)%, with energy out below 10⁻¹⁰ of energy in **[mainstream]**;
- LBNL and UC Davis (arXiv, Dec 2025) report sub-keV yields "over 10¹⁸" above bare-nucleus theory but low absolute rates **[contested]**;
- Szczecin claims thermal D–D fusion in ion tracks **[contested]**;
- NASA Glenn's 2020 "lattice confinement fusion" is photon-beam-driven **[contested]**;
- DARPA's 2026 MARRS programme targets 10⁶ reactions s⁻¹ g⁻¹, about 6×10⁻⁷ W/g (computed).

The exotic claims fared badly. Holmlid's laser-driven ultra-dense-hydrogen meson paper was retracted (PLOS ONE). Rossi's 2019 E-Cat SK demonstration is unvalidated.

I found no 2020–2026 tabletop claim of reaching the QED Schwinger field other than the Pais patents (Q1). The real Schwinger-type results are condensed-matter analogues (graphene, Nature Physics 2023) **[mainstream]**. The strongest lasers remain 2×10⁶ below the Schwinger intensity.

### Cited Findings
- **Chen, Maiwald, Schauer, …, Schenkel, Stolar, Brown & Berlinguette, Nature 644, 640–645 (2025), DOI 10.1038/s41586-025-09042-7.** Fifteen authors, from UBC, TRIUMF and LBNL.
  - Electrochemical loading of a Pd target bombarded by deuterium ions increased D–D fusion rates by "15(2)%". — [RePEc abstract](https://ideas.repec.org/a/nat/nature/v644y2025i8077d10.1038_s41586-025-09042-7.html)
  - The "Thunderbird Reactor" combines a microwave plasma thruster with plasma-immersion ion implantation into a Pd foil that is also the cathode of a D₂O electrolysis cell. — [Physics World](https://physicsworld.com/?p=123517)
  - The reactor produces less than about 10⁻¹⁰ of the energy it consumes. Physics World's text shows "10-10", most likely a lost superscript. — [Physics World](https://physicsworld.com/?p=123517)
  - Physics World calls it a "megaelectronvolt" accelerator. DARPA's slide labels the same demonstration a "beam driven fusion experiment (15 keV)". — [DARPA MARRS](https://www.darpa.mil/sites/default/files/attachment/2026-01/marrs-program-proposers-day-presentation.pdf). The two sources conflict on the beam energy; MeV refers to the products. **[mainstream]**
- **Karahadian et al., arXiv:2512.06212 (5 Dec 2025).** Electrochemical loading combined with ion bombardment produced D–D fusion below 2.5 keV in Pd and Ti hydrides, "enhanced by over 10¹⁸" relative to bare-nucleus theory, with low rates and no energy claim. — [arXiv:2512.06212](https://arxiv.org/abs/2512.06212v1) **[contested]**
- **Czerski et al., arXiv:2409.02112.** A plateau below 2.5 keV, interpreted as thermal D–D fusion in ion tracks, with U_e = 340 eV (Q4). — [arXiv:2409.02112](https://arxiv.org/abs/2409.02112v2) **[contested]**
- **NASA Glenn "lattice confinement fusion".** Two Phys. Rev. C 101 papers (April 2020), including Steinetz et al., "Novel nuclear reactions observed in bremsstrahlung-irradiated deuterated metals", PRC 101, 044610. Source for the items below: [Wikipedia: Lattice confinement fusion](https://en.wikipedia.org/wiki/Lattice_confinement_fusion); [NASA Glenn LCF page](https://www.nasa.gov/glenn/glenn-expertise-space-exploration/lattice-confinement-fusion/) (snippet)
  - ErD and TiD samples were irradiated by a dynamitron electron beam on a Ta converter (Wikipedia: about 2.2 MeV).
  - The proposed mechanism: deuteron photodissociation, then neutron-heated deuterons undergoing screened fusion or stripping.
  - The reactions were not self-sustaining, and the electron source was "energetically expensive". Wikipedia tags one claim as dubious. **[contested]**
- **Holmlid, ultra-dense hydrogen H(0).**
  - PLOS ONE retracted "Mesons from laser-induced processes in ultra-dense hydrogen H(0)" (2017). Reviewers cited conflict with baryon-number conservation and decay constants matching electronics pickup near intense lasers. — [PMC retraction notice](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6386439/) (snippet)
  - Hansen & Engelen argued that the 2022 muon source "does not produce any muons". — [Energy Sustain. Soc. comment](https://link.springer.com/doi/10.1186/s13705-023-00403-6) (snippet)
  - According to a critic's page, the University of Gothenburg removed its H(0) pages as of 1 Jan 2024. — [Nikhef page (critic)](https://www.nikhef.nl/~h02/udh3.html) (snippet). **[refuted]**
- **Rossi E-Cat SK** (31 Jan 2019, live stream). Critics called it "giant claims … with no validation". — [Wikipedia: Energy Catalyzer](https://en.wikipedia.org/wiki/Energy_Catalyzer) **[refuted/unvalidated]**
- **DARPA MARRS (2026)** cites Chen et al. 2025, Karahadian et al. 2025, Czerski 2022 and Metzler 2024 as the "recent advances" it builds on. — [DARPA MARRS](https://www.darpa.mil/sites/default/files/attachment/2026-01/marrs-program-proposers-day-presentation.pdf)
- **Press framing.** Coverage described the work as "pursuit of cold fusion". — [ANS Nuclear News, Aug 2025](https://www.ans.org/news/2025-08/article-7308/university-adds-electrochemical-boost-to-pursuit-of-cold-fusion/) (snippet). Berlinguette called it "the first clear experimental link between electrochemistry and fusion science". — [Physics World](https://physicsworld.com/?p=123517)
- **Schwinger analogue in graphene.** Schmitt, Plaçais et al., "Mesoscopic Klein-Schwinger effect in graphene", Nature Physics 19, 830–835 (2023). They report a universal 1D Schwinger conductance at pinch-off in ballistic graphene transistors, an analogue of pair creation for Dirac quasiparticles. — [ar5iv 2207.13400](https://ar5iv.labs.arxiv.org/html/2207.13400) (snippet) **[mainstream]** (analogue, not QED vacuum)
- **Laser and accelerator benchmarks.** Record intensity 1.1×10²³ W/cm² (2021); E144 positrons (1997); PVLAS at about 7× the QED prediction (2020) (see Q1). — [IBS](https://pr.ibs.re.kr/handle/8788114/10375) (snippet); [SLAC-PUB-7564](https://www.slac.stanford.edu/pubs/slacpubs/7500/slac-pub-7564.pdf) (snippet); [arXiv:2005.12913](https://arxiv.org/abs/2005.12913v1) (snippet) **[mainstream]**
- **History framing.** Messinger, Metzler & Price, "Gatekeeping: a Partial History of Cold Fusion" (arXiv:2601.09996, Jan 2026), describes a "modest renaissance" with US and EU funding. It is an opinion and history piece, not data. — [arXiv:2601.09996](https://arxiv.org/abs/2601.09996)

### Inferences
- **The new results are not cold fusion.** The 2025 Nature and arXiv results measure beam-driven D–D fusion, where the deuteron carries keV of kinetic energy. The condensed matter changes the rate by a factor (1.15×, or "10¹⁸" at sub-keV). They support the long-standing screening and loading physics of Q4, not the 1989 heat claims.
- **MARRS against heat claims** (computed). Even MARRS's Phase 2 target (10⁶ s⁻¹ g⁻¹ = 5.8×10⁻⁷ W/g) is about 2×10⁵ below F–P's 0.116 W/g. Any future validated result in this programme would still leave the historic heat claims unsupported by 5–12 orders of magnitude.
- **Schwinger regime** (computed; see Q1). To reach it on a tabletop, a static device would need fields about 10⁹ above vacuum breakdown, and a laser would need about 2×10⁶ more intensity than the 2021 record. Without a boosted electron beam (E144-type) or a strong-field facility (LUXE or ELI class), no tabletop "Schwinger-limit" claim is physically credible.
- **Search coverage.** I found no such 2020–2026 claim to test. Absence from a web search is weak evidence.
- **NASA LCF energy budget.** The photodisintegration threshold of the deuteron is 2.224 MeV (implied by F–P's corrected gamma line of 2.225 MeV). A beam endpoint near 2.2–3 MeV therefore lies just above the energy needed to make photoneutrons. Separating D–D fusion neutrons (2.45 MeV) from photoneutrons is the crux, and I could not check it without the PRC data (inference).

### Gaps
- The NASA PRC 2020 neutron yields, beam endpoint (Wikipedia says about 2.2 MeV) and any published critique were not retrieved.
- The beam energy in Chen et al. (2025) is unresolved (DARPA says 15 keV, Physics World says MeV), as are the absolute rates.
- The Karahadian preprint's absolute yields, energies and implied U_e were not extracted (abstract only).
- I did not search Japanese (Clean Planet/Tohoku) or Italian (ENEA) 2020–2026 excess-heat claims, or Brillouin Energy.
- I did not find the Holmlid retraction date or the full text of the PLOS notice (snippet only).
