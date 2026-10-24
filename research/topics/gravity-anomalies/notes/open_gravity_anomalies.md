# Open gravity anomalies and unexplained measurements: calculable models and the data that test them

Labels: **[mainstream]** = accepted by the field; **[contested]** = published and under live dispute; **[refuted]** = tested and failed. Every item also says **open** or **resolved**.
"Computed" means I recalculated it with `python3 -I` in `/tmp/claude-1000/-home-dan-code-FusionVault-quantum-gravimeter/7beb35b3-c723-440c-9450-66a991641392/scratchpad/calc/anomalies.py` (session scratch, so it may be deleted). Published values are quoted as printed. Constants used in computed lines: c = 299 792 458 m/s, ω⊕ = 7.292115×10⁻⁵ rad/s, R⊕ = 6371 km, GM☉ = 1.32712440018×10²⁰ m³ s⁻², 1 au = 1.495978707×10¹¹ m, H0 = 70 km/s/Mpc. These are standard values; I added them myself and they are not quoted from the sources below.

## 1. The Newtonian constant G: spread of precision measurements, CODATA values, and proposed explanations

### Takeaway
G is the worst-known major constant. The 16 CODATA input values span about 550 ppm (computed), while single experiments claim 12 to 20 ppm. CODATA 2018 and 2022 both give G = 6.674 30(15)×10⁻¹¹ m³ kg⁻¹ s⁻² (22 ppm) only after multiplying every uncertainty by 3.9. The spread is **[mainstream]** and **open**, and is put down to unidentified systematics. The 2026 NIST repeat of the BIPM apparatus came out about 250 to 280 ppm *lower* than BIPM-14, which points to a hidden systematic rather than new physics. The claim that G varies in step with length of day (Anderson et al. 2015) is **[refuted]**: statistical reanalyses strongly favour unknown noise over a periodic signal.

### Cited Findings
- CODATA 2022: no new G input datum since 2018, "the 2022 and 2018 recommended values of G are therefore identical." The 3.9 expansion factor from 2018 is reused. — [CODATA 2022, Mohr et al., arXiv:2409.03787 (Rev. Mod. Phys. 97, 025002 (2025))](https://arxiv.org/abs/2409.03787)
- Recommended value: G = 6.674 30(15)×10⁻¹¹ m³ kg⁻¹ s⁻², u_r = 2.2×10⁻⁵. After expansion, χ² = 12.9 for ν = 15, p = 0.61, Birge ratio R_B = 0.93. Before expansion, the largest normalized residuals were 7.7 (BIPM-14), 6.8 (JILA-18) and 4.8 (BIPM-01). — [CODATA 2022 PDF, Sec. XIV–XV](https://arxiv.org/pdf/2409.03787)
- CODATA 2022 Table XXX, the 16 inputs (G in 10⁻¹¹ m³ kg⁻¹ s⁻²; relative u). Correlations: r(NIST-82, LANL-97) = 0.351, r(HUST-05, HUST-09) = 0.134, r(HUST-09, HUSTT-18) = 0.068. — [CODATA 2022 PDF](https://arxiv.org/pdf/2409.03787)

  | ID | Method | G | u_r |
  |---|---|---|---|
  | NIST-82 (Luther & Towler) | fiber torsion, dynamic | 6.672 48(43) | 6.4e-5 |
  | TR&D-96 (Karagioz & Izmailov) | fiber torsion, dynamic | 6.672 9(5) | 7.5e-5 |
  | LANL-97 (Bagley & Luther) | fiber torsion, dynamic | 6.673 98(70) | 1.0e-4 |
  | UWash-00 (Gundlach & Merkowitz) | fiber torsion, dynamic compensation | 6.674 255(92) | 1.4e-5 |
  | BIPM-01 (Quinn et al.) | strip torsion, compensation + static | 6.675 59(27) | 4.0e-5 |
  | UWup-02 (Kleinevoß) | suspended body | 6.674 22(98) | 1.5e-4 |
  | MSL-03 (Armstrong & Fitzgerald) | strip torsion, compensation | 6.673 87(27) | 4.0e-5 |
  | HUST-05 (Hu et al.) | fiber torsion, dynamic | 6.672 22(87) | 1.3e-4 |
  | UZur-06 (Schlamminger et al.) | stationary body, weight change | 6.674 25(12) | 1.9e-5 |
  | HUST-09 (Luo; Tu) | fiber torsion, dynamic | 6.673 49(18) | 2.7e-5 |
  | BIPM-14 (Quinn et al. 2013, 2014) | strip torsion, compensation + static | 6.675 54(16) | 2.4e-5 |
  | LENS-14 (Prevedelli; Rosi et al.) | double atom interferometer gradiometer | 6.671 91(99) | 1.5e-4 |
  | UCI-14 (Newman et al.) | cryogenic torsion, dynamic | 6.674 35(13) | 1.9e-5 |
  | HUSTT-18 (Li et al., time-of-swing) | fiber torsion, dynamic | 6.674 184(78) | 1.2e-5 |
  | HUSTA-18 (Li et al., angular-acceleration feedback) | fiber torsion, dynamic compensation | 6.674 484(77) | 1.2e-5 |
  | JILA-18 (Parks & Faller) | suspended body, displacement | 6.672 60(25) | 3.7e-5 |
- BIPM (Quinn, Speake, Parks, Davis) 2013 published G = 6.67545(18)×10⁻¹¹ — [cited in arXiv:1504.07622](https://arxiv.org/pdf/1504.07622). An erratum found three errors (one was a density-gradient correction applied twice) that together raise G by 13 ppm — [PRL 113, 039901 (2014) erratum](https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevLett.113.039901/fulltext). CODATA uses the corrected 6.675 54(16) — [CODATA 2022](https://arxiv.org/pdf/2409.03787).
- Rosi et al. 2014, cold-atom interferometry: G = 6.67191(99)×10⁻¹¹ m³ kg⁻¹ s⁻², 150 ppm. — [Rosi et al., Nature 510, 518 (2014), arXiv:1412.7954](https://arxiv.org/abs/1412.7954)
- Li et al. 2018 (HUST), two methods: 6.674184(78) and 6.674484(78). The two differ by 2.7σ, which "suggests there could be sources of error unaccounted for." — [Wikipedia: Gravitational constant](https://en.wikipedia.org/wiki/Gravitational_constant). CODATA lists the second as 6.674 484(77). — [CODATA 2022](https://arxiv.org/pdf/2409.03787)
- Brack et al. (2022) measured G dynamically with resonating beams, but at u_r = 1.7×10⁻² it is "not competitive." Speake et al. (2023) tested whether stray AC magnetic fields caused the BIPM discrepancy and found that "unlikely." Atom-interferometry and displacement-sensor approaches are encouraged. — [CODATA 2022 PDF](https://arxiv.org/pdf/2409.03787)
- **NIST 2026 replication of the BIPM apparatus** (Schlamminger, Chao, Lee, Shakarji, Newell, Stirling, Cochrane, Speake; Metrologia 63, 025012 (2026), DOI 10.1088/1681-7575/ae570f, published 16 April 2026). The abstract gives G = (6.67366 ± 0.00020)×10⁻¹¹ m³ kg⁻¹ s⁻², "lower by 2.81×10⁻⁴ (relative) than the BIPM determination." The page states u_r = 3.0×10⁻⁶, which does not match ±0.00020 (that is 3.0×10⁻⁵; computed), so the page appears to contain a typo. — [NIST publication page](https://www.nist.gov/node/1911461)
  - **Conflict:** APS *Physics* news (1 May 2026) gives the NIST result as 6.67387×10⁻¹¹ with 57 ppm uncertainty, "250 ppm lower than the BIPM value," and within 1σ of the world average. Schlamminger is quoted: "We were not able to identify a single 'smoking gun,' so it's probably a multitude of causes." The article describes BIPM-14 as about 200 ppm above the world average. — [Physics 19, 64 (2026)](https://physics.aps.org/articles/v19/64). I could not resolve which number is final (the Metrologia full text was not accessible). Both are lower than BIPM-14.
  - A 2024 CPEM conference version identified a gas-pressure effect that had affected earlier NIST runs. — [NIST record via search](https://www.nist.gov/publications/redetermination-gravitational-constant-bipm-torsion-balance-nist)
- A Bayesian reanalysis with explicit systematic-error modelling recommends G = 6.6740 ± 0.0015 ×10⁻¹¹ m³ kg⁻¹ s⁻², which is "significantly larger" in uncertainty than CODATA 2018. — [Rinaldi, Middleton, Del Pozzo, Gair, arXiv:2209.07416](https://arxiv.org/pdf/2209.07416)
- **Time-variation claim:** Anderson et al. (EPL 110, 10002, 2015) reported a 5.9-year periodicity in G measurements correlated with the length of day (LOD). — [Pitkin comment, arXiv:1505.06725](https://arxiv.org/abs/1505.06725v1)
  - Schlamminger, Gundlach & Newman (2015) recompiled 35 years of G data. A 5.9-yr sine still fits better than a straight line, but their corrections "significantly weaken the correlation." — [arXiv:1505.01774](https://arxiv.org/abs/1505.01774v1)
  - Pitkin (2015), Bayesian model comparison: a constant plus an extra unknown Gaussian noise term is favoured over sinusoidal models by a Bayes factor of at least about e³⁰. — [arXiv:1505.06725](https://arxiv.org/abs/1505.06725v1)
  - A frequentist reanalysis (Desai, EPL 2016) supports Pitkin. — [arXiv:1607.03845](https://arxiv.org/pdf/1607.03845)
- An independent bound on secular variation from lunar laser ranging: Ġ/G = (7.1 ± 7.6)×10⁻¹⁴ yr⁻¹, cited as the most recent LLR limit (about 10× tighter than the earlier (4 ± 9)×10⁻¹³ yr⁻¹). — [Belgacem et al., arXiv:1812.11181](https://arxiv.org/pdf/1812.11181)

#### Test vectors (computed, reproduce CODATA exactly)
- Generalized-least-squares mean of the 16 Table XXX values with the 3 correlations: **G = 6.674300**. Unexpanded u = 0.000038 (5.7 ppm). **χ² = 195.7** for ν = 15 (Birge ratio 3.61). ×3.9 gives **u = 0.00015 (22.1 ppm)** and χ²/3.9² = **12.9**, matching CODATA's χ² = 12.9.
- Normalized residuals (y − Ḡ)/u_i: BIPM-14 **+7.75**, JILA-18 **−6.80**, BIPM-01 **+4.78**. These match CODATA's 7.75, −6.80 and 4.8.
- Spread max − min = 0.00368 (551 ppm). HUSTA-18 − HUSTT-18 = 45 ppm = 2.7σ. Quinn 2013 → 2014 = +13.5 ppm.
- NIST-2026 vs BIPM-14: −282 ppm = −7.3σ (Metrologia numbers), or −250 ppm = −4.0σ (APS numbers). NIST-2026 vs CODATA: −96 ppm, or −64 ppm on the APS number.
- Rosi 2014 vs CODATA 2022: −358 ppm.

### Inferences
- The NIST replication used the same apparatus but new mass metrology and environment, and it did not reproduce the BIPM value. That makes an apparatus- or site-specific systematic the most economical explanation for the high BIPM-01/BIPM-14 results. A gravitational or time-varying cause is less likely, because the same physics would have to act differently on one instrument in two labs.
- A periodic G with LOD period would need ΔG/G of order 10⁻⁴ over about 6 years, if the scatter in the table were real signal. LLR limits secular drift to about 10⁻¹³ yr⁻¹. Purely periodic variation is not directly bounded by that number, so the decisive rebuttal is statistical (Pitkin, Desai), not LLR.
- For a Rust calculator: implement GLS with a full covariance, and a Birge-ratio expansion as a separate step. The test vectors above are a strong end-to-end check.

### Gaps
- I could not open the Metrologia 2026 full text to settle 6.67366(20) vs 6.67387 (57 ppm). They may be different stages (preprint vs final, or one method vs combined); this is unverified.
- Rothleitner & Schlamminger 2017 (Rev. Sci. Instrum. 88, 111101) is cited by CODATA as a key review. I did not fetch its numbers.
- I did not find the fitted amplitude of the Anderson et al. 2015 sinusoid in an accessible source.

## 2. The Earth flyby anomaly

### Takeaway
Six Earth flybys from 1990 to 2005 showed mm/s-level changes in asymptotic speed. Anderson et al. (2008) fitted them with ΔV∞/V∞ = K(cos δᵢ − cos δₒ), K = 2ω⊕R⊕/c = 3.099×10⁻⁶. The data are real tracking residuals **[mainstream]**. The formula is **[contested]** and failed its out-of-sample tests: Rosetta 2007/2009, Juno 2013 and BepiColombo 2020 all came out null. ESA's 2026 reanalysis reproduces the 2005 Rosetta anomaly, rules out ephemeris, gravity-field and relativistic causes, and proposes albedo-triggered outgassing. Status: **open**, now usually treated as spacecraft-specific rather than new physics.

### Cited Findings
- Empirical formula: ΔV∞/V∞ = K(cos δᵢ − cos δ_f), K = 2ω⊕R⊕/c = 3.099×10⁻⁶. δᵢ and δ_f are the incoming and outgoing declinations of the asymptotic velocity. The formula fits flybys below 2000 km altitude, "perhaps fortuitously"; "a phenomenological fit is not an explanation." — [Nieto & Anderson, "Earth Flyby Anomalies", arXiv:0910.1321](https://arxiv.org/pdf/0910.1321)
- Constants as stated in the literature: K = 3.099×10⁻⁶ with ω_E = 7.292115×10⁻⁵ rad/s and R_E = 6371 km. — [Mbelek, arXiv:0809.1888](https://arxiv.org/pdf/0809.1888v3)
- Anderson et al. 2008 PRL 100, 091102, Table I (reproduced). Columns: altitude (km), V∞ (km/s), δᵢ, δₒ (deg), observed ΔV∞ ± σ (mm/s), predicted (mm/s). — [reproduced in arXiv:2406.06604, Table 1](https://arxiv.org/pdf/2406.06604)

  | Date | Craft | Alt | V∞ | δᵢ | δₒ | Obs ΔV∞ | σ | Pred |
  |---|---|---|---|---|---|---|---|---|
  | 1990-12-08 | Galileo I | 960 | 8.949 | −12.52 | −34.15 | 3.92 | 0.3 | 4.12 |
  | 1992-12-08 | Galileo II | 303 | 8.877 | −34.26 | −4.87 | −4.6 | 1 | −4.67 |
  | 1998-01-23 | NEAR | 539 | 6.851 | −20.76 | −71.96 | 13.46 | 0.01 | 13.28 |
  | 1999-08-18 | Cassini | 1175 | 16.01 | −12.92 | −4.99 | −2 | 1 | −1.07 |
  | 2005-03-04 | Rosetta I | 1956 | 3.863 | −2.81 | −34.29 | 1.8 | 0.03 | 2.07 |
  | 2005-08-02 | MESSENGER | 2347 | 4.056 | 31.44 | −31.92 | 0.02 | 0.01 | 0.06 |
- **Uncertainty conflicts between secondary sources:** Wikipedia and Páramos & Hechenblaikner give Galileo I 3.92 ± 0.08, NEAR 13.46 ± 0.13 and Rosetta 1.80 ± 0.05 (Wikipedia: 1.82 ± 0.05). — [Wikipedia: Flyby anomaly](https://en.wikipedia.org/wiki/Flyby_anomaly); [Páramos & Hechenblaikner, arXiv:1210.7333](https://ar5iv.arxiv.org/html/1210.7333). Nieto & Anderson's own text gives Galileo I as (4.13 ± 0.03) mm/s. — [arXiv:0910.1321](https://arxiv.org/pdf/0910.1321)
- Later flybys: Rosetta 2007 (V∞ 9.36 km/s, 5322 km) and Rosetta 2009 (9.38 km/s, 2483 km) gave ΔV∞ ≈ 0. — [arXiv:1210.7333](https://ar5iv.arxiv.org/html/1210.7333)
- Juno 2013-10-09 was predicted to show an anomaly of order 6 mm/s by the empirical formula (Iorio, range-rate). Standard J₂, post-Newtonian and Lense–Thirring signatures are ≤ about 25 µm/s. — [Iorio, Adv. Space Res. 54, 2441 (2014), arXiv:1311.4218](https://arxiv.org/abs/1311.4218). The Anderson-type value is about +6 mm/s; Busack's alternative gives −7 mm/s. — [arXiv:1312.1139](https://arxiv.org/pdf/1312.1139)
- Juno result: "The 2013 flyby of Juno on the way to Jupiter yielded no anomalous acceleration" (table entry 0 ± 0.8 mm/s). Source: Thompson et al., *Reconstruction of Earth flyby by the Juno spacecraft*, AAS 14-435 (2014). — [Wikipedia: Flyby anomaly](https://en.wikipedia.org/wiki/Flyby_anomaly)
- BepiColombo, April 2020: "no unexpected accelerations were found" (Iess et al., Space Sci. Rev. 217, 21 (2021)). Wikipedia's table lists 0 ± 0.003 mm/s, which I could not confirm in the primary source. — [Wikipedia: Flyby anomaly](https://en.wikipedia.org/wiki/Flyby_anomaly)
- Acedo (2017) also says the Anderson formula does not predict the null results for Rosetta II/III and Juno. He fits a tentative empirical radial-velocity/c force that decays with altitude (scale about 1000 km) and says further speculation should wait for more data. — [Acedo, arXiv:1711.02875](https://ar5iv.labs.arxiv.org/html/1711.02875)
- Páramos & Hechenblaikner: an Earth-oblateness error is implausible, because the needed J₂ change is outside its constrained range. Poorly modelled spacecraft behaviour should not be dismissed. — [arXiv:1210.7333](https://ar5iv.arxiv.org/html/1210.7333)
- **2026, ESA/ESOC (Hoffmann & Budnik, ISSFD 2026):** orbit determination with ESOC's new deep-space navigation system reproduces the Rosetta 2005 anomaly. It rules out ephemeris errors, gravity-field modelling, relativistic effects and propagator settings. Solar radiation pressure and albedo pressure are insufficient. Albedo illumination of thermally sensitive surfaces "could have triggered outgassing" matching cruise-phase outgassing events. A data package is planned. No numbers are given in the abstract. — [ISSFD 2026 abstract](https://www.issfd2026.fr/data/onglet12/module0/modalPreview.php?langue=fr&paramProjet=20840)
- Other proposed explanations: Earth-bound dark matter, special or general relativity modifications, light-speed anisotropy, frame-dragging (cannot account for it), range-proportional telemetry delay, gravity-model truncation (Juno needed at least 50×50). — [Wikipedia: Flyby anomaly](https://en.wikipedia.org/wiki/Flyby_anomaly); [Nieto & Anderson](https://arxiv.org/pdf/0910.1321)

#### Test vectors (computed)
- K = 2·7.292115e-5·6.371e6 / 299792458 = **3.0993×10⁻⁶**.
- ΔV∞ = K·V∞·(cos δᵢ − cos δₒ), with V∞ in m/s giving m/s. Galileo I **+4.123**, Galileo II **−4.674**, NEAR **+13.279**, Cassini **−1.068**, Rosetta I **+2.067**, MESSENGER **+0.055** mm/s. All match the paper's predicted column to its 2 decimals.
- Observed ΔV/V: NEAR **1.965×10⁻⁶**, Galileo I 4.380×10⁻⁷, MESSENGER 4.93×10⁻⁹.

### Inferences
- A formula with one universal constant and no spacecraft parameters should have predicted Juno at about +6 mm/s; it got 0. The near-null for BepiColombo is stronger still. These are the out-of-sample tests, and the formula fails them. That points to a spacecraft- or tracking-specific cause, consistent with ESA's outgassing hypothesis.
- ESA's 2026 work says conventional *gravitational* modelling is ruled out for Rosetta 2005. The remaining candidate is non-gravitational (thermal or outgassing), the same class of answer that solved Pioneer (Section 4).

### Gaps
- Juno's V∞ and declinations, and the exact Thompson et al. residual, were not retrieved from the primary AAS paper.
- No numbers from the ESA 2026 reanalysis are public in the abstract.
- Hayabusa2 (2015), OSIRIS-REx (2017) and other later Earth flybys: I found no published anomaly values.

## 3. Gravitational shielding / absorption: Majorana, Allais, eclipse gravimetry, LLR bounds

### Takeaway
Majorana's absorption law F = F_N·exp(−h∫ρ dl) needs h about 10⁻¹³ m²/kg to match his 1920 claim. Lunar laser ranging bounds h < 10⁻²² m²/kg (Eckhardt 1990), and later LLR gives h = (3 ± 5)×10⁻²² m²/kg. That is about 9 orders of magnitude below Majorana's value. Majorana's shielding is **[refuted]** and **resolved**. The Allais pendulum effect and the 1997 eclipse gravimeter dips are **[contested]** and unreplicated: a 2009 superconducting-gravimeter campaign saw nothing, and atmospheric or local causes match the size. Read as shielding, the 7 µGal 1997 signal would need h about 10⁷ times the LLR bound (computed).

### Cited Findings
- Majorana's law: F = (GMm/r²)·e^(−h∫ρ(r)dr), with h an extinction coefficient in m²/kg (1920). — [Wikipedia: Gravitational shielding](https://en.wikipedia.org/wiki/Gravitational_shielding)
- Majorana's first experiments gave h ≈ 6.7×10⁻¹² cm²/g, later revised to about half that. This comes from a search snippet of the PDF; the abstract has no numbers. — [von Borzeszkowski & Chrobok, arXiv:gr-qc/0310081](https://arxiv.org/pdf/gr-qc/0310081)
- Bounds on h (m²/kg): — [Wikipedia: Gravitational shielding](https://en.wikipedia.org/wiki/Gravitational_shielding)
  - Poincaré 1908 (lunar motion): h ≤ 10⁻¹⁸
  - Unnikrishnan & Gillies 2000 (Zürich G experiment): h < 4.3×10⁻¹⁵, the best laboratory bound
  - Yang & Wang 2002 (1997 eclipse data): constraint 6×10⁻¹⁹
  - **Eckhardt 1990 (lunar ranging): h < 10⁻²²**
  - **Williams, Turyshev & Boggs 2009 (LLR): h = (3 ± 5)×10⁻²²**
- In cgs, LLR gives h ≤ 1.0×10⁻²¹ cm²/g, and the Zürich bound is h ≤ 4.3×10⁻¹⁴ cm²/g. Both are consistent with the m²/kg values (1 cm²/g = 0.1 m²/kg). — [arXiv:hep-ph/0204284](https://arxiv.org/pdf/hep-ph/0204284)
- Allais reported anomalous Foucault-pendulum precession in the 30 June 1954 eclipse and paraconical-pendulum anomalies in the 2 October 1959 eclipse. Later results are mixed:
  - Saxl & Allen 1970 reported anomalies.
  - Kuusela 1991 and Jun et al. 1991 (1990 eclipse) found none.
  - Salva 2011 (2010 eclipse, automated Foucault pendulum) found no change above 0.3°/h.
  - The NASA 1999 campaign was never published.
  - Van Flandern & Yang estimate the gravity anomaly is "about a factor of 100,000 too small" to explain Allais's precession.
  - Status per Wikipedia: controversial, inconsistent results for over five decades. — [Wikipedia: Allais effect](https://en.wikipedia.org/wiki/Allais_effect)
- Wang et al. 2000 (PRD 62, 041101): a LaCoste-Romberg D gravimeter at Mohe during the 9 March 1997 eclipse saw near-symmetric dips around first and last contact. The anomaly is (7.0 ± 2.7)×10⁻⁸ m/s², and "there may be a shielding property of gravitation." — [arXiv:1003.4947](https://arxiv.org/abs/1003.4947)
- Unnikrishnan, Mohapatra & Gillies 2001 (PRD 63, 062002): "Anomalous gravity data during the 1997 total solar eclipse do not support the hypothesis of gravitational shielding." — [Wikipedia: Gravitational shielding](https://en.wikipedia.org/wiki/Gravitational_shielding). A search summary attributes to this paper the suggestion that seismic disturbance from crowds moving before and after the eclipse could explain the dips. I did not verify this in the paper. — [allais.wiki copy of PRD 63, 062002](https://allais.wiki/itscom/blackprior/misc/unnikrishnan.pdf)
- Van Flandern & Yang 2003 (PRD 67, 022002): cooling of the upper atmosphere during the eclipse produces a moving air-mass excess giving about 7×10⁻⁸ m/s², the right magnitude. — [arXiv:physics/0610197 summary](https://arxiv.org/pdf/physics/0610197v1)
- 22 July 2009 eclipse (Hubei): a superconducting gravimeter (Jiufeng) plus LCR-G and LCR-ET gravimeters. "No 'gravity anomaly valleys' of abnormal amplitude were observed." — [Cui, Fang, Wang, Geomatics & Inf. Sci. Wuhan Univ. 36(11), 1332 (2011)](https://ch.whu.edu.cn/en/article/id/700)

#### Test vectors (computed)
- Model: fractional loss = 1 − exp(−h·Σ), with Σ the diametral column density 2Rρ.
  - Earth: Σ⊕ = 2·6371 km·5514 kg/m³ = **7.026×10¹⁰ kg/m²**.
  - Moon: Σ☾ = 2·1737.4 km·3344 kg/m³ = **1.162×10¹⁰ kg/m²**. The radii and densities are standard values I assumed.
- Majorana h = 6.7×10⁻¹³ m²/kg: loss through Earth **4.60%**, through Moon **0.78%**.
- Zürich bound 4.3×10⁻¹⁵: Earth 3.0×10⁻⁴, Moon 5.0×10⁻⁵.
- Eckhardt bound 10⁻²²: Earth **7.0×10⁻¹²**, Moon **1.16×10⁻¹²**.
- Sun's pull at 1 au: GM☉/au² = **5.930×10⁻³ m/s²**.
- The 1997 dip of 7.0×10⁻⁸ m/s², read as lunar shielding of the Sun, implies h = 7.0e-8 / 5.930e-3 / 1.162e10 = **1.0×10⁻¹⁵ m²/kg**, about 10⁷ × the LLR bound.
- Maximum eclipse signal allowed by LLR: 10⁻²²·Σ☾·g☉ = **6.9×10⁻¹⁵ m/s² (6.9×10⁻⁷ µGal)**, far below any gravimeter.

### Inferences
- The LLR bound makes any eclipse gravimeter "shielding" signal unobservable by about 10 orders of magnitude. Reported µGal eclipse dips therefore cannot be absorption and must be environmental: atmospheric mass, thermal tilt, or human activity. The 2009 superconducting-gravimeter null supports this.
- Majorana's h would make the Moon's orbit feel a few-percent modulation in the Sun's pull during eclipses, and would change Earth's effective mass seen through the planet. That is why lunar dynamics (Poincaré, then LLR) kill it.

### Gaps
- Eckhardt 1990 (PRD 42, 2144 is my recollection, not verified) was not fetched directly. The 10⁻²² m²/kg value is via Wikipedia and is consistent with hep-ph/0204284.
- I found no published superconducting-gravimeter results for the 2017 or 2024 North American eclipses.
- The Yang & Wang 2002 eclipse constraint (6×10⁻¹⁹) is quoted only from Wikipedia, and its derivation is unclear. It conflicts with my computed 1.0×10⁻¹⁵ needed to *explain* the 1997 dip, because a bound and an explanation use different assumptions.

## 4. The Pioneer anomaly and its thermal resolution (worked example of a solved anomaly)

### Takeaway
Pioneer 10/11 showed a constant-looking sunward acceleration a_P = (8.74 ± 1.33)×10⁻¹⁰ m/s². A finite-element thermal model fed with flight telemetry (Turyshev et al. 2012) reproduces its size, direction and slow decay. Once thermal recoil is included, "no anomalous acceleration remains." New-physics readings are **[refuted]**; the thermal explanation is **[mainstream]**; status **resolved**. The lesson: the decisive evidence was the *time dependence* (a decaying heat source), which a new constant force could not mimic.

### Cited Findings
- a_P = (8.74 ± 1.33)×10⁻⁸ cm/s² (= 8.74×10⁻¹⁰ m/s²) directed toward the Sun, from Pioneer 10 data 1987–1998 plus Pioneer 11. Extensive searches for systematic errors had not found a cause. — [Anderson, Laing, Lau, Liu, Nieto, Turyshev, PRD 65, 082004 (2002), arXiv:gr-qc/0104064](https://arxiv.org/abs/gr-qc/0104064)
- Thermal resolution: the thermal recoil force's magnitude, time behaviour and direction resemble the anomaly. A parameterized thermal-recoil fit to the Doppler data shows no statistically significant difference from the thermal model, and "once the thermal recoil force is properly accounted for, no anomalous acceleration remains." — [Turyshev, Toth, Kinsella, Lee, Lok, Ellis, PRL 108, 241101 (2012), arXiv:1204.2507](https://arxiv.org/abs/1204.2507)
- Model and numbers from the paper: — [ar5iv of arXiv:1204.2507](https://ar5iv.labs.arxiv.org/html/1204.2507)
  - a(t) = 2^(−(t−t₀)/τ)·a₀, with t₀ = 1 Jan 1980.
  - Thermal model: τ = 36.9 ± 6.7 yr, fit a₀ = (7.4 ± 2.5)×10⁻¹⁰ m/s².
  - Doppler: τ = 28.8 ± 2.0 yr, a₀ = (10.1 ± 1.0)×10⁻¹⁰ m/s².
  - ²³⁸Pu half-life 87.74 yr.
  - Efficiency factors: thermal η_rtg = 0.0104, η_elec = 0.406; Doppler best fit η_rtg = 0.0144, η_elec = 0.480.
  - The thermal error budget is 26.3%, dominated by RTG surface degradation (25%).
  - The 1σ error ellipses overlap, so "no statistically significant acceleration anomaly exists."

#### Test vectors (computed)
- cH₀ (H₀ = 70) = **6.80×10⁻¹⁰ m/s²**; a_P / cH₀ = 1.29. This is the numerical coincidence that fuelled cosmological readings.
- Thermal a₀ / a_P = 0.85.
- Doppler model at t₀ + 20 yr: 10.1·2^(−20/28.8) = **6.24×10⁻¹⁰ m/s²**.

### Inferences
- The decay time τ ≈ 29–37 yr is much shorter than the ²³⁸Pu half-life (87.74 yr). That fits electrical power and RTG output falling faster than the isotope decays (thermocouple degradation), which a thermal origin explains and a cosmological constant acceleration cannot.
- For the flyby anomaly (Section 2), the same method applies: a time-resolved, spacecraft-specific thermal or outgassing model tested against the Doppler residuals. ESA's 2026 Rosetta reanalysis is this step.

### Gaps
- Independent confirmations (Rievers & Lämmerzahl 2011; Francisco, Bertolami, Gil, Páramos 2012) were not fetched in this session.

## 5. Low-acceleration gravity: MOND a₀, the radial acceleration relation, and the wide-binary dispute

### Takeaway
The radial acceleration relation (RAR) is an empirical law: across 2693 points in 153 galaxies, the observed acceleration is a single function of the baryonic acceleration, g_obs = g_bar/(1 − e^(−√(g_bar/g†))), with g† = 1.20 ± 0.02 (random) ± 0.24 (syst) ×10⁻¹⁰ m/s². The RAR is **[mainstream]** as an observation; MOND as its explanation is **[contested]**. The wide-binary test in Gaia DR3 is **[contested]** and **open**:
- Chae finds a gravity boost of about 1.3–1.6 below about 10⁻⁹ m/s².
- Banik et al. find Newton preferred at 19σ.
- Chae's 2026 high-quality 36-binary result is driven by two systems.

### Cited Findings
- RAR: "the same relation is followed by 2693 points in 153 galaxies," and it persists where dark matter dominates. — [McGaugh, Lelli, Schombert, PRL 117, 201101 (2016), arXiv:1609.05917](https://arxiv.org/abs/1609.05917)
- Fit: g_obs = F(g_bar) = g_bar / (1 − e^(−√(g_bar/g†))), with g† = 1.20 ± 0.02 (random) ± 0.24 (systematic) ×10⁻¹⁰ m s⁻². The systematic term comes from a 20% normalization uncertainty in Υ⋆. Residuals are Gaussian with width 0.11 dex (rms 0.13 dex including outliers), and the data are "consistent with negligible intrinsic scatter." Υ⋆(3.6 µm) is 0.50 for disks and 0.7 for bulges. — [ar5iv of arXiv:1609.05917](https://ar5iv.labs.arxiv.org/html/1609.05917)
- Chae 2023 (ApJ 952, 128): 26,615 Gaia DR3 wide binaries within 200 pc. δ_obs-newt = 0.034 ± 0.007 at g_N ≈ 10^−8.91 and 0.109 ± 0.013 at g_N ≈ 10^−10.15 m/s². At 10^−10.15, g_obs/g_pred = 10^(√2·δ) = 1.43 ± 0.06. Combined significance 10σ. — [arXiv:2305.04613](https://arxiv.org/abs/2305.04613)
- Banik, Pittordis, Sutherland, Famaey, Ibata, Mieske, Zhao 2024 (MNRAS 527, 4573): 8611 wide binaries within 250 pc, separations 2–30 kAU. α_grav = −0.021 (+0.065/−0.045), consistent with Newton. MOND is excluded at 16σ (MCMC), and Newton is preferred over the best MOND model at 19σ. MOND predicted about 20% higher orbital velocities at large separation (external field effect included). — [arXiv:2311.03436](https://arxiv.org/abs/2311.03436)
- Chae 2025 (ApJ, Bayesian 3D orbit modelling, 312 binaries with radial velocities): — [arXiv:2502.09373](https://arxiv.org/abs/2502.09373)
  - High-acceleration regime (10^−7.9 to 10^−6.7): Γ = 0.000 ± 0.011.
  - MOND regime (10^−11 to 10^−9.5, N = 35): Γ = 0.085 ± 0.040, γ_g = 1.48 (+0.33/−0.23).
  - MOND + transition regime (N = 111): Γ = 0.063 ± 0.015, γ_g = 1.34 (+0.10/−0.08).
  - The MOND-regime Γ agrees with the MOND prediction of about 0.07.
- Chae, Lee, Hernandez, Orlov, Lim, Turnshek, Lee (arXiv:2601.21728, v3 Sept 2026): 36 binaries within 150 pc, RV error < 100 m/s. Γ = 0.102 (+0.028/−0.021) for g < 10^−9.1, and Γ = 0.132 (+0.035/−0.027) for g < 10^−9.5. γ = 1.60 (+0.24/−0.14), which "rules out Newton and is higher than relevant MOND predictions." The signal is dominated by two systems (HD 189739/HD 189760; TYC 259-236-1/TYC 259-906-1) whose 3D velocities exceed Newtonian escape velocity. — [arXiv:2601.21728](https://arxiv.org/abs/2601.21728)
- Commentary: the dispute centres on sample selection (very wide pairs vs contamination) and velocity-uncertainty cuts. Banik: "I disagree very much with the results of Chae." No consensus. — [Physics World](https://physicsworld.com/?p=109496); [Big Think](https://develop.bigthink.com/?p=452954)

#### Test vectors (computed)
- RAR with g† = 1.20×10⁻¹⁰. Columns: g_bar → g_obs (ratio):

  | g_bar | g_obs | ratio |
  |---|---|---|
  | 1e-12 | **1.1462e-11** | 11.462 |
  | 1e-11 | **3.9881e-11** | 3.9881 |
  | 1e-10 | **1.6705e-10** | 1.6705 |
  | 1e-9 | **1.0590e-9** | 1.0590 |
  | 1e-8 | **1.0001e-8** | 1.0001 |
- Small-g_bar limit of the same formula: g_obs → √(g_bar·g†). Example: √(1e-12·1.2e-10) = 1.0954e-11, which approaches the exact 1.1462e-11 as g_bar → 0.
- Flat rotation speed in that limit: v⁴ = G·M·g†. For M = 10¹¹ M☉: v = **199.8 km/s**. Radius where g_N = g†: √(GM/g†) = 10.8 kpc.
- Wide binary with 2 M☉ total: g_N = g† at r = √(GM/g†) = **9,942 au** (about 10 kAU), inside Banik's 2–30 kAU window.
- Conversions:
  - Chae 2023: 10^(√2·0.109) = **1.426**; 10^(√2·0.034) = 1.117.
  - Chae 2025/2026 use γ_g = 10^(2Γ): Γ = 0.085 → **1.479**, 0.063 → 1.337, 0.102 → 1.600, 0.132 → 1.837, MOND 0.07 → 1.380. These match the published γ values.
  - Banik's 20% velocity boost corresponds to a gravity boost of 1.2² = **1.44**.
- cH₀/2π (H₀ = 70) = **1.082×10⁻¹⁰ m/s²**, close to g†.

### Inferences
- Both camps agree that MOND predicts a gravity boost of about 1.4 at about 10 kAU separations. The disagreement is purely about data quality: contamination by hidden triples and chance alignments versus genuine signal. That makes it a decidable question with better radial velocities and astrometry.
- The 2026 result rests on two systems moving faster than escape velocity. That pattern is the expected fingerprint of unresolved companions or projection effects, so it needs follow-up spectroscopy before it counts.
- For a calculator: implement F(g_bar) and its inverse numerically. Check that F → g_bar for g_bar ≫ g† and F → √(g_bar g†) for g_bar ≪ g†.

### Gaps
- Gaia DR4 timing and its expected effect on the wide-binary test were not verified in this session.
- Hernandez's and Saad & Ting's independent wide-binary results, and Banik's formal rebuttal of Chae 2025/2026, were not fetched.
- I did not fetch cluster- or cosmology-scale tests of MOND, which are outside this scope.

## 6. Antimatter gravity and the equivalence principle: ALPHA-g, MICROSCOPE, Eöt-Wash

### Takeaway
- ALPHA-g (2023) found antihydrogen falls down: a_g = [0.75 ± 0.13 (stat+syst) ± 0.16 (sim)] g, which is 1.2σ from +1 g and 8.5σ from −1 g (computed). Repulsive antigravity is **[refuted]**; the sign is **resolved**, while the magnitude at the percent level is **open**.
- MICROSCOPE's final result η(Ti,Pt) = [−1.5 ± 2.3 (stat) ± 1.5 (syst)]×10⁻¹⁵ shows no equivalence-principle violation **[mainstream]**.
- Eöt-Wash torsion balances give η(Be,Ti) = (0.3 ± 1.8)×10⁻¹³ and bound gravitational-strength Yukawa forces to λ < 38.6 µm.

### Cited Findings
- ALPHA-g: "the local gravitational acceleration of antihydrogen is directed towards Earth and has magnitude a_g = [0.75 ± 0.13 (stat. + syst.) ± 0.16 (sim.)]g," "consistent with the attractive gravitational force between matter and Earth." The next step targets about 1% using laser-cooled antihydrogen. — [CERN Courier](https://cerncourier.com/?p=109604); paper: E. K. Anderson et al., Nature 621, 716–722 (2023), doi:10.1038/s41586-023-06527-1 — [sci.news summary](https://www.sci.news/physics/antimatter-gravity-12309.html)
- Method: about 100 antihydrogen atoms per trap fill, released over 20 s, with up/down escape fractions compared to simulation across bias-field settings. — [sci.news](https://www.sci.news/physics/antimatter-gravity-12309.html)
- Earlier ALPHA bounds were ±7500% (2012) and ±100% (2013); the 2023 margin is about 20%. — [Wikipedia: Gravitational interaction of antimatter](https://en.wikipedia.org/wiki/Gravitational_interaction_of_antimatter)
- MICROSCOPE final: η(Ti,Pt) = [−1.5 ± 2.3 (stat) ± 1.5 (syst)]×10⁻¹⁵, no violation of the weak equivalence principle. — [Touboul et al., PRL 129, 121102 (2022), arXiv:2209.15487](https://arxiv.org/abs/2209.15487)
- Eöt-Wash rotating torsion balance (Be–Ti): — [Schlamminger, Choi, Wagner, Gundlach, Adelberger, PRL 100, 041101 (2008), arXiv:0712.0607](https://arxiv.org/abs/0712.0607)
  - η = (0.3 ± 1.8)×10⁻¹³ toward Earth
  - Δa = (0.6 ± 3.1)×10⁻¹⁵ m/s²
  - galactic dark-matter η = (−4 ± 7)×10⁻⁵
  - space-fixed differential acceleration < 8.8×10⁻¹⁵ m/s² (95%)
- Eöt-Wash 1/r² test: separations 52 µm to 3.0 mm. Gravitational-strength Yukawa interactions are limited to ranges < 38.6 µm (95% CL). — [Lee, Adelberger, Cook, Fleischer, Heckel, PRL 124, 101101 (2020), arXiv:2002.11761](https://arxiv.org/abs/2002.11761)

#### Test vectors (computed)
- ALPHA-g with the two errors in quadrature: σ = √(0.13² + 0.16²) = **0.206 g**. (1 − 0.75)/σ = **1.21σ**; (0.75 + 1)/σ = **8.5σ**.
- MICROSCOPE: σ = √(2.3² + 1.5²) = **2.75×10⁻¹⁵**.

### Inferences
- ALPHA-g's 0.75 g is not a hint of reduced antimatter gravity at this precision (1.2σ). Simulation systematics dominate the error, so the 1% upgrade has to improve modelling, not just counting statistics.
- Theories in which antimatter or Cooper pairs couple anomalously to gravity at the 10⁻⁴ level (see Section 7) must be checked for consistency with these bounds. Those bounds constrain bulk ordinary matter, though, not exotic condensates.

### Gaps
- The full Nature paper (paywalled redirect) was not read. The uncertainty breakdown is from CERN Courier.
- GBAR and AEgIS 2024–2026 results were not found in this session.
- The LLR Nordtvedt-parameter bound and atom-interferometer equivalence-principle tests (e.g. Asenbaum 2020) were not fetched.

## 7. The Tate anomaly in the Cooper-pair mass

### Takeaway
Tate et al. (1989) measured the London moment of a rotating niobium ring and inferred a Cooper-pair mass m*/2m_e = 1.000084(21). Theory predicts about 8 ppm *below* 2m_e, so the discrepancy is about 92 ppm, or 4.4σ (computed). It has never been repeated with similar precision. Status: **open** and **[contested]**. The Tajmar–de Matos gravitomagnetic explanation needs a field of B_g = 1.84×10⁻⁴·ω, about 10 orders of magnitude above anything known, and is **[contested]**, with no independent confirmation.

### Cited Findings
- Tate, Cabrera, Felch & Anderson, "Precise Determination of the Cooper-Pair Mass," PRL 62, 845–848 (1989): m₀/2m_e = 1.000084(21). — [Gallerati, Modanese, Ummarino review, arXiv:2206.07574](https://arxiv.org/pdf/2206.07574)
- Theory predicts a value about 8 ppm *less* than 2m_e; the measurement is 84 ppm above. — [Tajmar & de Matos, arXiv:gr-qc/0607086 (search snippet)](https://arxiv.org/abs/gr-qc/0607086)
- Tajmar & de Matos's Ginzburg–Landau analysis with a gravitomagnetic term gives B_g = 2ω(m* − m)/m. For thick Nb this is **B_g = 1.84×10⁻⁴ ω**, "very large" compared with Earth's gravitomagnetic field of about 10⁻¹⁴ rad/s. — [arXiv:2206.07574](https://arxiv.org/pdf/2206.07574)
- A US Air Force-funded test of high-temperature superconductors (BSCCO, YBCO) found no gravitational anomalies down to the facility's resolution. This does not test Tate's niobium setup. — [DTIC ADA461570](https://apps.dtic.mil/sti/html/tr/ADA461570/index.html)

#### Test vectors (computed)
- (1.000084 − 0.999992)/0.999992 = **92 ppm**; ÷ 21 ppm = **4.4σ**.
- B_g / ω = 2·9.2×10⁻⁵ = **1.84×10⁻⁴**, which reproduces the published value.
- London moment B = 2m_eω/e = **1.1371×10⁻¹¹ T per rad/s**. With m* in place of 2m_e it is B = m*ω/e, the quantity Tate measured.

### Inferences
- A 92 ppm error in a single 1989 experiment that was never repeated looks like an undiagnosed systematic or a theory gap in the work-function and contact-potential corrections. This is the same pattern as the G spread: one precision apparatus, no replication.
- If the gravitomagnetic reading were right, rotating superconductors would produce frame-dragging fields 10¹⁰ × Earth's. That should be detectable with modern gyroscopes and ring lasers. I found no independently confirmed detection in this session.

### Gaps
- **"Hoang et al. 2020" could not be found.** No 2020 paper on the Cooper-pair mass or London moment by that name turned up. The assignment may refer to a different author or year; I did not invent a result.
- Tajmar's own later null results (ca. 2008–2011) and the Canterbury ring-laser null were not fetched.
- The exact theoretical value (0.999992) is from a secondary snippet, not Tate's paper.

## 8. Anomalies reported 2023–2026 in gravity measurements (gravimeters, clocks, satellite geodesy) and how they fared

### Takeaway
Nothing reported from 2023 to 2026 overturns standard gravity. The main developments:
1. **NIST 2026** failed to reproduce BIPM's high G. Systematics are now favoured; the spread is still open.
2. **Wide binaries** remain contested and open (Section 5).
3. **ESA 2026** reanalysed the Rosetta flyby. The anomaly is real in the data and not gravitational; the outgassing hypothesis is open.
4. A **GRACE/SLR 2007 transient gravity anomaly** at the Atlantic–Africa margin (EGU 2025) is a geophysical hypothesis, not new physics. It is contested and open.
5. **Lattice atom-interferometer** (2024) and other fifth-force searches came out null.

The older **LLR anomalous lunar eccentricity rate** is still listed as open but has shrunk with better tidal models.

### Cited Findings
- NIST 2026 G redetermination: see Section 1. Result 6.67366(20) (Metrologia abstract) or 6.67387 at 57 ppm (APS), 250–281 ppm below BIPM-14; "no single 'smoking gun'." — [NIST](https://www.nist.gov/node/1911461); [APS Physics 19, 64](https://physics.aps.org/articles/v19/64)
- Speake et al. (2023) tested stray AC magnetic fields as the cause of the BIPM G offsets and concluded it was unlikely. — [CODATA 2022](https://arxiv.org/pdf/2409.03787)
- GRACE + SLR, 2003–2015: a significant transient north–south gravity anomaly at the Atlantic/Africa boundary, peaking January 2007, timescale 2–3 years. Surface hydrology and ocean mass do not explain it, so the origin is suggested to be solid-Earth. The proposed cause is vertical displacement of the perovskite → post-perovskite transition near the African LLSVP base. It coincides with the 2007 geomagnetic jerk. Presented at EGU25, followed up at EGU26. — [EGU25-8808](https://meetingorganizer.copernicus.org/EGU25/EGU25-8808.html); [EGU26-12420](https://meetingorganizer.copernicus.org/EGU26/EGU26-12420.html)
- Rosetta flyby ESA reanalysis 2026: see Section 2. — [ISSFD 2026](https://www.issfd2026.fr/data/onglet12/module0/modalPreview.php?langue=fr&paramProjet=20840)
- Lattice atom interferometry (Berkeley-led, Nature 2024): the attraction of a small source mass agrees with Newtonian gravity, ruling out screened fifth-force theories over their natural parameter space. — [Physics Today](https://physicstoday.aip.org/news/atom-interferometer-measures-centimeter-scale-gravity)
- A 2025 proposal for a 10 m-scale interferometer would extend chameleon and symmetron bounds by 1 to 1.5 orders of magnitude. This is a proposal, not a measurement. — [arXiv:2511.09750](https://arxiv.org/abs/2511.09750v1)
- Chae et al. 2026 wide-binary high-quality sample: see Section 5. — [arXiv:2601.21728](https://arxiv.org/abs/2601.21728)
- LLR lunar eccentricity: Williams & Boggs (2009, 16,941 ranges 1970–2008, DE421) found an anomalous de/dt of about (9 ± 3)×10⁻¹² yr⁻¹. Williams & Dickey (2003) had about (1.6 ± 0.5)×10⁻¹¹ yr⁻¹. With DE430 and better tides, Iorio reports about (5 ± 2)×10⁻¹² yr⁻¹. Neither standard tidal physics nor tested modified-gravity or cosmological models explain it. — [Iorio, arXiv:1102.0212](https://ar5iv.labs.arxiv.org/html/1102.0212); [Iorio, arXiv:1404.6537](https://arxiv.org/pdf/1404.6537); [Anderson & Nieto, arXiv:0907.2469](https://arxiv.org/pdf/0907.2469). Note: the DE430-era number comes via Iorio, not the JPL team directly.

### Inferences
- The pattern across 2023–2026 is the Pioneer pattern: anomalies tied to a single instrument or analysis (BIPM G, Rosetta 2005, Tate) shrink or move to systematics when someone independently re-does them.
- The live exceptions are statistical-astronomy tests (wide binaries) and slowly accumulating geodetic series (LLR eccentricity rate), where the "instrument" is a population or a decades-long model.
- For this project (satellite gravimetry, GRACE-FO): the GRACE 2007 transient is the one recent satellite-geodesy "anomaly." It is a mass-redistribution signal inside Earth, so it is a target for the Rust tooling as a geophysics validation case, not a gravity-physics test.

### Gaps
- No 2023–2026 optical-clock or gravimeter result claiming a gravity anomaly turned up. A dedicated clock-network search was not run.
- I found no peer-reviewed journal version of the GRACE 2007 anomaly, only EGU abstracts.
- The 2024 lattice-interferometer result's exact citation and bounds were not fetched; it is from Physics Today coverage.
- JPL's own current statement on the LLR eccentricity rate was not retrieved.
