# Gravity and the Schwinger effect: Unruh, Hawking, Gibbons–Hawking, their analogues, and vacuum-inertia claims

Conventions. SI units with CODATA 2018 constants. M☉ = 1.9884×10³⁰ kg, R☉ = 6.957×10⁸ m, L☉ = 3.828×10²⁶ W, 1 yr = 3.15576×10⁷ s.
"computed" means I recomputed the value in `gravity_unruh_hawking_checks.py` (same folder; run it with `python3 -I`; it uses only the standard library and contains asserts).
"Published" values are quoted exactly as the source gives them. Labels: **[mainstream]**, **[contested]**, **[refuted]**.
Related sibling notes: `../Gravity beams fringe physics and DIRDs/baker_li_hfgw_and_pais_patents.md` covers the Pais patents in detail, and `dirds_gravity_and_metric_engineering.md` covers Puthoff's PV metric and DIRD 15. These notes do not repeat them.

## Q1. Unruh, Schwinger-as-thermal, Hawking, charged black holes and de Sitter: formulas and test numbers

### Takeaway
All four effects share one structure: a horizon or an acceleration a fixes a temperature ħa/(2πck_B), and pair creation is Boltzmann-suppressed at that temperature. The Schwinger exponent πE_S/E equals mc²/(2k_BT_U) exactly, where T_U is the Unruh temperature for a = eE/m. So the flat-space Schwinger rate is a Boltzmann factor for one rest mass at **2T_U**, not at T_U. Whether this is "the same effect" or a formal correspondence is debated (the "factor-of-2 problem"). Every gravitational version is unmeasurably cold for astrophysical objects: a solar-mass black hole is at 61.7 nK, the cosmological horizon at about 2×10⁻³⁰ K, and Earth's surface gravity corresponds to 4×10⁻²⁰ K.

### Cited Findings

**Formula sheet (each row has a test value)**

| # | Quantity | Exact form (SI) | Source | Published test value | Computed |
|---|---|---|---|---|---|
| 1 | Unruh temperature [K] | T_U = ħa/(2πck_B), a = proper acceleration [m s⁻²] | [Wikipedia, Unruh effect](https://en.wikipedia.org/wiki/Unruh_effect); [Crispino, Higuchi & Matsas, RMP 80, 787 (2008), arXiv:0710.5373](https://arxiv.org/abs/0710.5373) | "≈ 4.06 × 10⁻²¹ K·s²·m⁻¹ × a"; "2.47×10²⁰ m⋅s⁻² corresponds approximately to a temperature of 1 K" (Wikipedia) | 4.055×10⁻²¹ K s² m⁻¹; a(1 K) = 2.466×10²⁰ m s⁻² = 2.515×10¹⁹ g₀ |
| 2 | Schwinger critical field [V m⁻¹] | E_S = m_e²c³/(eħ) | [Schmitt et al., Nat. Phys. 19, 830 (2023), arXiv:2207.13400](https://arxiv.org/abs/2207.13400) | "1.32 × 10¹⁸ V/m" with Δ_S = mc² = 511 keV | 1.323×10¹⁸ V m⁻¹ |
| 3 | Acceleration and Unruh T at E_S | a_S = eE_S/m = mc³/ħ; T_U(E_S) = m_ec²/(2πk_B) | from rows 1 and 2 | Crispino Sec. IV.D: near the Schwinger limit "the Unruh temperature becomes comparable to the electron mass"; such fields give accelerations "as high as 10²⁸ g⊕" | a_S = 2.327×10²⁹ m s⁻²; T_U = 9.44×10⁸ K = 81.3 keV |
| 4 | Schwinger pair rate (spin ½, 3+1 D) [pairs m⁻³ s⁻¹] | Γ = (eE)²/(4π³ħ²c) Σₙ n⁻² exp(−nπE_S/E) | [Schwinger, Phys. Rev. 82, 664 (1951)](https://doi.org/10.1103/PhysRev.82.664). Schmitt et al. give the d-dimensional form w ∝ Σₙ (E/n)^((d+1)/2) e^(−πnE_S/E), which for d = 3 is E²/n² | — | Γ(E_S) = 4.75×10⁵⁴; Γ(0.1E_S) = 2.47×10⁴⁰ m⁻³ s⁻¹ |
| 5 | Schwinger ↔ Unruh identity | πE_S/E ≡ m c²/(2k_B T_U(a = eE/m)) | Cai & Kim, arXiv:1407.4569, Eq. (15): flat 2D "N_S = e^(−m/(2T_U))"; their Sec. 5: the dS/AdS results "reduce to the Schwinger formula with T=2T_U" | — | identity checked to 10⁻¹² for three fields |
| 6 | Hawking temperature [K] | T_H = ħc³/(8πGMk_B) | [Wikipedia, Hawking radiation](https://en.wikipedia.org/wiki/Hawking_radiation) (cites Hawking 1974, *Nature* 248, 30) | "about 3 solar masses, this temperature is 10⁻⁷ K"; Earth mass "10⁻¹² K" | 1 M☉: 6.17×10⁻⁸ K = 61.7 nK; 3 M☉: 2.06×10⁻⁸ K (Wikipedia's 10⁻⁷ K is about 5× high, order-of-magnitude only) |
| 7 | Photon-only evaporation time [s] | t_ev = 5120πG²M³/(ħc⁴); P = ħc⁶/(15360πG²M²) | Wikipedia, Hawking radiation | "≈ 2.140×10⁶⁷ years × (M/M☉)³" | 2.095×10⁶⁷ yr (2% below the published figure; the cause is not identified); P(1 M☉) = 9.0×10⁻²⁹ W |
| 8 | Page 1976 power, M ≫ 10¹⁷ g | P = 2×10⁻⁴ ħc⁶G⁻²M⁻² (81% ν, 17% γ, 2% graviton). Integrating gives τ = G²M³/(3·2×10⁻⁴ ħc⁴) | [Page, PRD 13, 198 (1976), abstract](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.13.198) | "a hole of mass M≫10¹⁷ g should emit a total power output of 2×10⁻⁴ℏc⁶G⁻²M⁻²"; a PBH has decayed within the present age "if and only if its initial mass was M<(5±1)×10¹⁴ g" | τ = 8.72×10⁻²⁷ (M/1 g)³ s; 1 M☉: 2.17×10⁶⁶ yr, 9.65× shorter than row 7 because neutrinos dominate and greybody factors are included |
| 9 | PBH evaporating now | M* | MacGibbon, Carr & Page, PRD 78, 064043 (2008), as cited on Wikipedia, Hawking radiation | "(5.00±0.04)×10¹¹ kg" | T_H(5×10¹¹ kg) = 2.45×10¹¹ K = 21.1 MeV. The massless-only Page τ gives 34.5 Gyr, so the extra species emitted at about 20 MeV (μ, π, …) are what bring τ down to the age of the universe |
| 10 | Black hole in equilibrium with the CMB | M = ħc³/(8πGk_B T_CMB) | row 6 with T_CMB = 2.7255 K ([Fixsen 2009](https://arxiv.org/abs/0911.1955)) | — | 4.50×10²² kg (about 0.6 lunar masses). Heavier holes currently gain mass from the CMB |
| 11 | Gibbons–Hawking temperature [K] | T_GH = ħH/(2πk_B) | [Gibbons & Hawking, PRD 15, 2738 (1977)](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.15.2738): an observer's detector "would see thermal radiation apparently coming from the cosmological event horizon" (paraphrased abstract); Cai & Kim Eq. (14), T_GH = H/(2π) | Planck 2018: H₀ = 67.4 ± 0.5 km s⁻¹ Mpc⁻¹, Ω_m = 0.315 ± 0.007 ([arXiv:1807.06209](https://arxiv.org/abs/1807.06209)) | T_GH(H₀) = 2.66×10⁻³⁰ K; asymptotic de Sitter (H₀√Ω_Λ, Ω_Λ = 0.685) = 2.20×10⁻³⁰ K |
| 12 | Schwinger effect in dS₂ | N_dS = (e^(−(S_μ−S_λ)) + e^(−2S_μ))/(1 − e^(−2S_μ)), with S_μ = 2πγ_dS/H, S_λ = 2πqE/H², γ_dS² = (qE/H)² + m² − H²/4 (natural units). Effective T_dS = √(T_GH² + T_U²) + T_U, T_U = qE/(2πm) | [Cai & Kim, arXiv:1407.4569](https://ar5iv.labs.arxiv.org/html/1407.4569), Eqs. (6), (14), (15), (17), (30), (61)–(62) | strong-field limit T_dS ≈ 2T_U | H = 10⁻³m: N_dS/e^(−πm²/qE) = 1.000 for qE = 0.5 and 2 (m = 1). At E = 0, N_dS → e^(−2πm/H) = e^(−m/T_GH) (ratio 1.17 at H = 0.2m, from the −H²/4 term) |
| 13 | Extremal Reissner–Nordström horizon field | E_h = c⁴/(√(4πε₀) G^(3/2) M), with r₊ = GM/c² and Q = M√(4πε₀G) | derived from the RN metric; discharge physics from [Gibbons, CMP 44, 245 (1975)](https://link.springer.com/doi/10.1007/BF01609829) (paywalled; abstract seen via search) | Gibbons: for holes "more massive than 10¹⁷ g" the charge loss "is governed by a Schwinger type formula" | E_h = E_S at M = 1.06×10³⁶ kg = 5.3×10⁵ M☉. For a 10 M☉ extremal hole E_h = 7.1×10²² V m⁻¹ (5×10⁴ E_S) |

**Status statements and the factor-of-2 debate**
- Crispino et al. 2008 (Sec. IV): "the Unruh effect does not need experimental confirmation any more than free Quantum Field Theory does". **[mainstream]** — [arXiv:0710.5373](https://arxiv.org/abs/0710.5373)
- Nikishov & Ritus (1988), as reviewed in Crispino Sec. IV.D: photon emission by an electron in a constant field E obeys the detailed balance P(p⊥→p′⊥)/P(p′⊥→p⊥) = exp(−βΔE), with β⁻¹ = a/2π and a = eE/m_e. This is the cleanest published bridge between QED in a constant field and Unruh thermality. The review warns that the link between Rindler energy and transverse energy "is not entirely clear". **[mainstream as a result; the interpretation is contested]** — [Crispino et al.](https://arxiv.org/abs/0710.5373)
- Labun & Rafelski (PRD 86, 041701, 2012) assign a vacuum-fluctuation temperature T_M to a constant field. T_M equals T_U for gyromagnetic factor g = 1, and they conjecture that pair production is captured by the g = 1 case. **[contested]** — [arXiv:1203.6148](https://arxiv.org/abs/1203.6148)
- Volovik describes the naive match exp(−πm²/qE) = exp(−πm/a) as a "factor of 2 problem": the Schwinger exponent implies twice the Unruh temperature. In 2024 he argues that the Unruh, de Sitter and Hawking temperatures all acquire a factor of 2 from coherence between processes. **[contested]** — [arXiv:2206.02799](https://arxiv.org/abs/2206.02799) (JETP Lett., DOI 10.1134/S0021364022601968); [arXiv:2409.03316](https://arxiv.org/abs/2409.03316)
- Parentani & Massar (1997), "The Schwinger mechanism, the Unruh effect and the production of accelerated black holes", is the standard reference connecting the three. **[mainstream]** — [INSPIRE 416584](https://inspirehep.net/literature/416584)
- For an accelerated detector in de Sitter space the literature uses T = √(T_U² + T_GH²), that is √(a² + H²)/2π (Narnhofer–Peter–Thirring 1996; Deser–Levin for AdS). This is summarised in Kim's review. **[mainstream]** — [Kim, arXiv:1512.03142](https://arxiv.org/abs/1512.03142); Narnhofer et al. are listed in [Crispino et al.](https://arxiv.org/abs/0710.5373) references. (The search-engine summary I saw called this a "geometric mean". The formula is a root-sum-square.)
- Kim (2016), "Schwinger Effect, Hawking Radiation, and Unruh Effect", extends the thermal reading to charged black holes. **[mainstream]** — [arXiv:1602.05336](https://arxiv.org/abs/1602.05336) (not read in full)
- Hawking radiation from astrophysical or primordial black holes has not been detected. The Fermi telescope has seen no primordial-black-hole bursts ("none have been detected"), and no micro black holes have been seen at CERN. **[mainstream status]** — [Wikipedia, Hawking radiation](https://en.wikipedia.org/wiki/Hawking_radiation)

**Gravity-scale Unruh numbers (computed)**
- T_U(a = g₀ = 9.80665 m s⁻²) = 3.98×10⁻²⁰ K. T_U at the Sun's surface gravity (274.0 m s⁻², [NASA](https://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html)) = 1.11×10⁻¹⁸ K.

### Inferences
- For calculators, implement the identity in row 5 as a unit test: `schwinger_exponent(E) == m c² / (2 k_B T_U(eE/m))`. Expose two separate functions, `t_unruh(a)` and `t_schwinger_eff(E) = 2*t_unruh(eE/m)`. This keeps the factor of 2 explicit rather than buried.
- "Gravity as Schwinger" is a mathematical correspondence (Euclidean periodicity, an exponent with a horizon), not a mechanism that lets gravitational fields produce pairs at useful rates. Pair production at a horizon needs κ ~ mc³/ħ, that is T_H ~ m_ec²/(2πk_B) = 9.4×10⁸ K. For electrons this is a black hole of about 1.3×10¹⁴ kg (computed from row 6). No laboratory gravitational field comes within 10²⁸ of the required acceleration (2.3×10²⁹ m s⁻²).
- Charged astrophysical black holes cannot stay near-extremal below about 5×10⁵ M☉ (row 13). Above that mass they are discharged by accreting plasma rather than by Schwinger emission. This is consistent with the standard assumption Q ≈ 0 used for the Kerr metric.

### Gaps
- I could not access the full text of Gibbons 1975 (paywalled), so the exact discharge criterion and mass threshold as he states them are unverified. Row 13 is my own derivation.
- I did not obtain the original Hawking 1975 (CMP 43, 199) or Page 1976 coefficient tables. Page's lifetime coefficient (8.72×10⁻²⁷ s g⁻³) was computed from the abstract's rounded "2×10⁻⁴". Page's own, more precise value should be checked before using this as a test vector at better than 1%.
- I could not identify the source of the 2% gap between my 2.095×10⁶⁷ yr and Wikipedia's 2.140×10⁶⁷ yr (constants or year definition; not resolved).
- The journal reference for Cai & Kim 1407.4569 was not shown on the page I read (probably JHEP 2014; unverified).

## Q2. Laboratory Unruh proposals and claimed detections (Bell–Leinaas, Chen–Tajima, 2010s–2020s)

### Takeaway
No direct, undisputed laboratory detection of the Unruh effect exists as of 2026 (Wikipedia: "It is currently not clear whether the Unruh effect has actually been observed, since the claimed observations are disputed"). Storage-ring spin polarisation (Sokolov–Ternov, 92.38%) can be described in the electron's frame with an effective temperature, but circular motion does not give a true Unruh thermal bath. The two polarisation formulas differ "substantially". Positron-channelling data (Lynch et al. 2021) is the strongest "real-system" claim and remains contested. Everything since then is a proposal: cavities, mirrors, Josephson junctions, geometric phases, superradiance.

### Cited Findings
- **Sokolov–Ternov (1963) baseline:** P(t) = P₀(1 − e^(−t/τ)), with P₀ = 8/(5√3) "≈ 92.38%". The build-up time in natural units is τ = 8m²ρ³/(5√3 e²γ⁵), which in SI is τ = 8m_e²cρ³/(5√3 αħ²γ⁵). **[mainstream]** — [Crispino et al. Sec. IV.A](https://arxiv.org/abs/0710.5373)
  - Published: polarisation equilibrium "could be achieved in a couple of hours, half an hour and 10 minutes" at LEP, HERA and SPEAR respectively.
  - Computed (bend-only, with my assumed bending radii, which are not from a source): SPEAR 3.7 GeV with ρ = 12.7 m gives 4.9 min; LEP 45.6 GeV with ρ = 3096 m gives 4.1 h. The ring-average factor R/ρ > 1 lengthens both.
- **Bell & Leinaas (1983, Nucl. Phys. B 212, 131, "Electrons as accelerated thermometers"; 1987, NPB 284, 488):** they proposed using the spin as a thermometer for comoving observers. For linear acceleration they give Γ± = ∓(8μ²/3) ΔE(ΔE² + a²)/(1 − exp(±2πΔE/a)), with μ = ge/(4m_e). — [Crispino et al.](https://arxiv.org/abs/0710.5373)
  - Linac case: the review states "Stanford linac with an accelerating field of 10 MV/m", acceleration "2 × 10¹⁶ g⊕", and T "= 0.7 × 10⁻³ K". The polarisation build-up time is "much larger than the lifetime of the Universe".
  - Computed: 10 MV m⁻¹ gives a = eE/m_e = 1.76×10¹⁸ m s⁻² = 1.79×10¹⁷ g₀, and T_U = 7.1 mK. The review's a and T are self-consistent (2×10¹⁶ g₀ gives 0.80 mK) but match 1 MV m⁻¹, not 10 MV m⁻¹. Treat this as an inconsistency in the review.
  - Circular motion: the review gives an effective temperature β⁻¹/a = (ΔE/a)/ln[1 + 4√3(ΔE/a)exp(2√3ΔE/a)]. Its limits are β⁻¹/a → 1/(2√12) (ΔE ≪ a) and 1/(2√3) (ΔE ≫ a), compared with Unruh's 1/(2π). The "thermal" polarisation P₁ = (1 − e^(−πg))/(1 + e^(−πg)) (Eq. 4.4) and the exact Derbenev–Kondratenko P₂ (Eq. 4.5) "show substantial differences". The effective temperature "cannot be considered as the temperature of a legitimate thermal bath". **[mainstream that the analogy is imperfect]**
  - Computed for LEP (45.6 GeV, ρ ≈ 3.1 km): a = γ²c²/ρ = 2.3×10²³ m s⁻², giving T_U ≈ 940 K, or about 1700 K for the ΔE ≫ a circular effective temperature.
- **Akhmedov & Singleton (2007, JETP Lett.)** argue that the Unruh effect "physically coincides" with the experimentally verified Sokolov–Ternov effect, so that it may already have been observed. Their companion paper says the exponential Unruh term is "hidden", and its v3 says it is unobservable. **[contested]** — [arXiv:0705.2525](https://arxiv.org/abs/0705.2525); [arXiv:hep-ph/0610391](https://arxiv.org/abs/hep-ph/0610391); [arXiv:hep-th/0607026](https://arxiv.org/abs/hep-th/0607026) ("No Unruh-like Radiation for circular motion")
- **Rogers (1988) Penning trap:** v = 0.6c, B = 1.5×10⁵ G, ω_z ≈ 5×10¹⁰ s⁻¹, ω₀ ≈ 2×10¹² s⁻¹, a = γ²vω₀ ≈ 6×10¹⁹ g⊕, "a/2π ≈ 2 K". Computed: ω₀ = 2.11×10¹² s⁻¹, a = 6.05×10¹⁹ g₀, T_U = 2.4 K (agrees). Never realised. **[proposal]** — [Crispino et al. Sec. IV.B](https://arxiv.org/abs/0710.5373)
- **Scully et al. (2003) / Belyanin et al. (2006), atoms through a cavity:** α ≈ 10⁸ Hz (3×10¹⁵ g⊕) and ΔE ≈ 10¹⁰ Hz give an excitation ratio e^(−2πΔE/α) ~ 10⁻²⁰⁰. Sharp cavity boundaries raise the ratio to α/(2πΔE) ~ 10⁻³, but that enhancement is a boundary effect, not the Unruh effect. I did not reproduce the 10⁻²⁰⁰. **[proposal]** — [Crispino et al. Sec. IV.C](https://arxiv.org/abs/0710.5373)
- **Chen & Tajima (PRL 83, 256, 1999), "Testing Unruh Radiation with Ultraintense Lasers":** Unruh radiation "is detectable, in principle, beyond the conventional radiation (most notably the Larmor radiation) background noise". — [APS abstract](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.83.256)
  - Model (via the review): a standing wave, a = 2a₀ω₀cos ω₀t, a₀ = eE₀/(m_eω₀), Larmor energy per half cycle ΔI_L = (4π/3)e²a₀²ω₀. The estimate is ΔI_U/ΔI_L ≈ (9/π²)(ħω₀/m_ec²) a₀ log(a₀/π) "≈ 3 × 10⁻⁴" for ω₀ ≈ 2×10¹⁵ s⁻¹ and a₀ ≈ 100. — [Crispino Sec. IV.D](https://arxiv.org/abs/0710.5373)
  - Computed: with ln the ratio is 8.1×10⁻⁴; with log₁₀ it is 3.5×10⁻⁴. The published 3×10⁻⁴ matches log₁₀, so check the original before coding it. Peak a = 2a₀ω₀c = 1.2×10²⁶ m s⁻², T_U = 4.9×10⁵ K (42 eV). **[proposal; never done]**
  - Schützhold et al. (2006) proposed correlated two-photon emission inside forward and backward cones (where Larmor radiation has a blind spot) as the signature. **[proposal]** — [Crispino et al.](https://arxiv.org/abs/0710.5373)
- **Present laser frontier:** Yoon et al., Optica 8, 630 (2021) reached "(1.1 ± 0.2) × 10²³ W/cm²" with a 1.1 μm spot. Computed: E_peak = √(2I/(cε₀)) = 9.1×10¹⁴ V m⁻¹ = E_S/1454, so the Schwinger exponent is πE_S/E ≈ 4.6×10³ and spontaneous pair production is nil. — [IBS repository](https://pr.ibs.re.kr/handle/8788114/10375); [Optica OPN](https://www.optica-opn.org/home/newsroom/2021/may/adaptive_optics_spurs_record_laser_intensity)
- **Lynch, Cohen, Hadad & Kaminer, PRD 104, 025015 (2021):** analyses positrons channelled in silicon. They find "a Larmor formula and power spectrum that are both thermalized by the acceleration" and claim "the first observation of acceleration-induced thermality in a non-analogue system". Wikipedia attributes the data to CERN NA63 and says the effect "may have been observed". **[contested]** — [arXiv:1903.00043](https://arxiv.org/abs/1903.00043); [Wikipedia, Unruh effect](https://en.wikipedia.org/wiki/Unruh_effect)
- **Other disputes:** Smolyaninov claims Unruh radiation has been observed; Ford & O'Connell argue that a uniformly accelerated charge emits none. **[contested]** — [Wikipedia, Unruh effect](https://en.wikipedia.org/wiki/Unruh_effect)
- **Quantum simulation:** Hu, Feng, Zhang & Chin, Nat. Phys. (2019), s41567-019-0537-1, "Quantum Simulation of Coherent Hawking-Unruh Radiation". A parametrically driven BEC of "≈10⁵ atoms" emits pairs, and the extracted temperature and entropy "agree with Unruh's predictions". **[mainstream as a simulation; not the gravitational or vacuum effect]** — [arXiv:1807.07504](https://arxiv.org/abs/1807.07504)
- **Proposals 2011–2025 (none realised):**
  - Martín-Martínez, Fuentes & Mann, PRL (2011): Berry-phase detection with current technology. — [Wikipedia](https://en.wikipedia.org/wiki/Unruh_effect)
  - Gregori et al., EPJC (2024): "Measuring Unruh radiation from accelerated electrons". — [Strathprints 89183](https://strathprints.strath.ac.uk/89183)
  - Barman, Ghosh & Majhi (2024/25): geometric phase enhanced by mirrors at lower acceleration. — [arXiv:2405.07711](https://arxiv.org/abs/2405.07711)
  - Zheng et al. (Dec 2024): a circular Unruh–DeWitt detector in a cylindrical cavity with superradiance. They note that "experimental confirmation remains elusive". — [arXiv:2412.17353](https://arxiv.org/abs/2412.17353)
  - Deswal et al., PRL (2025), DOI 10.1103/6z1l-kkmk: a time-resolved, superradiantly amplified burst between mirrors that needs "significantly lower" acceleration. — [phys.org coverage](https://phys.org/news/2025-11-faint-quantum-space.html)
  - Hiroshima University, PRL (July 2025): "Circular-Motion Fulling-Davies-Unruh Effect in Coupled Annular Josephson Junctions", in which fluxon–antifluxon decay shows up as a voltage jump. — [Hiroshima U. news](https://www.hiroshima-u.ac.jp/node/93889) (the search also returned [node 92582](https://www.hiroshima-u.ac.jp/node/92582); I did not check which node is the PRL item)
  - **[proposals]**

### Inferences
- Test vectors for an `unruh` crate: (a) 1 K ↔ 2.466×10²⁰ m s⁻²; (b) Rogers' 2.4 K; (c) the 10 MV m⁻¹ linac at 7.1 mK (flag the review's 0.7 mK); (d) Sokolov–Ternov P₀ = 0.92376 and τ_ST ∝ ρ³γ⁻⁵ (SPEAR at 4.9 min bend-only); (e) the circular β⁻¹/a curve with limits 0.1443 and 0.2887.
- Every proposal faces the same scaling: T_U = 1 K needs about 2.5×10¹⁹ g. Only charged particles in strong EM fields (rings, lasers, crystal channels) get there, and their ordinary Larmor or synchrotron background is described completely by inertial-frame QED. So any "detection" is a reinterpretation of QED data, not a new effect. This is why Crispino et al. treat the Unruh effect as needing no confirmation.

### Gaps
- I did not read the Bell–Leinaas 1983/1987 originals or Unruh's 1998 "Acceleration radiation for orbiting electrons". The numbers come from the 2008 review.
- I found no published comment or rebuttal specifically refuting Lynch et al. 2021, only Wikipedia's "disputed" status. The NA63 dataset details (positron energy, crystal) are not in the abstract.
- I did not verify the 2025 PRLs (Deswal; Hiroshima) beyond their abstracts and press coverage, and I have no numerical thresholds for them.
- I did not examine the AnaBHEL laser "flying mirror" analogue project (Chen, Mourou et al.).

## Q3. Analogue gravity: BEC, water-wave, optical and graphene "Schwinger" experiments

### Takeaway
Analogue experiments test the kinematics of horizons (QFT on an effective metric), not gravity itself. Steinhauer's Rb-87 BEC has measured spontaneous Hawking correlations (2016, entanglement contested), a thermal spectrum at the predicted T_H = 0.35 nK (2019), and stationarity (2021). Water-wave (2011) and fibre-optic (2019) analogues observed stimulated, classical Hawking-type mode conversion. Graphene superlattices (Berdyugin 2022) and ballistic graphene transistors (Schmitt 2023) show "Schwinger-like" electron–hole pair creation with E_S reduced from 1.3×10¹⁸ to about 6×10⁷ V m⁻¹, because v_F ≈ c/300 and the effective gap is about 0.2 eV.

### Cited Findings
- **Analogue Hawking temperature (Steinhauer formula):** T_H = ħg/(2πck_B) with surface gravity g = c(dv/dx − dc/dx) evaluated at the sonic horizon. For stationary 1D flow (nv constant): k_BT_H = −(ħ/2π)(c/n · dn/dx + dc/dx)|ₓ₌₀ (Eq. 1). **[mainstream]** — [de Nova, Golubkov, Kolobov & Steinhauer, Nature 569, 688 (2019), arXiv:1809.00913](https://arxiv.org/abs/1809.00913)
- **Steinhauer 2019 numbers:** "The experiment is repeated 7400 times". The "predicted Hawking temperature of 0.35 nK, or 0.12 mc²_out". The correlation band lies at "5.5º from the hydrodynamic angle". The spectrum is found "to be thermal at the Hawking temperature implied by the analogue surface gravity", and the correlations "imply that the analogue black hole has no analogue firewall". **[mainstream]** — [arXiv:1809.00913](https://arxiv.org/abs/1809.00913)
  - Computed (Rb-87, m = 86.909 u): k_BT_H = 0.12mc_out² gives c_out = 0.53 mm s⁻¹ and κ = 2πk_BT_H/ħ = 288 s⁻¹. A real 1 M☉ black hole (61.7 nK) is 176× hotter than this analogue.
- **Steinhauer 2016, Nat. Phys. 12, 959:** "We observe spontaneous Hawking radiation, stimulated by quantum vacuum fluctuations". "high energy pairs are entangled, while the low energy pairs are not". — [arXiv:1510.00621](https://arxiv.org/abs/1510.00621)
  - Disputed by Leonhardt, "Questioning the recent observation of quantum Hawking radiation", Annalen der Physik (2018). Leonhardt raised "severe doubts" on statistical significance and on consistency with theoretical (Heisenberg) bounds. — [arXiv:1609.03803](https://arxiv.org/abs/1609.03803)
  - Steinhauer replied that the note "suffers from technical difficulties which invalidate its main claims". — [arXiv:1609.09017](https://arxiv.org/abs/1609.09017)
  - In 2019 Leonhardt called the thermal-spectrum work "an important step for the community". — [Physics World](https://physicsworld.com/a/thermal-spectrum-of-analogue-black-hole-puts-hawking-radiation-in-a-new-light/)
  - **[2016 entanglement: contested; 2019 thermality: mainstream-accepted as an analogue result]**
- **Kolobov, Golubkov, de Nova & Steinhauer, Nat. Phys. 17, 362 (2021):** spontaneous Hawking radiation is stationary, "observing such a system at six different times". The end is "marked by the formation of an inner horizon". Simulations favour Bogoliubov–Cherenkov–Landau stimulation "rather than black-hole lasing". **[mainstream]** — [arXiv:1910.09363](https://arxiv.org/abs/1910.09363)
- **Water waves, Weinfurtner, Tedford, Penrice, Unruh & Lawrence, PRL 106, 021302 (2011):** a streamlined obstacle in open-channel flow blocks long waves and converts them to short waves. The measured amplitudes of the converted waves support "the thermal nature" of the conversion. This is stimulated and classical, not spontaneous. **[mainstream]** — [arXiv:1008.1911](https://arxiv.org/abs/1008.1911)
- **Optics, Drori, Rosenberg, Bermudez, Silberberg & Leonhardt, PRL 122, 010404 (2019):** Kerr-induced moving horizons in fibre: "We have observed in our experiment that the probe stimulates Hawking radiation". **[mainstream, stimulated]** — [arXiv:1808.09244](https://arxiv.org/abs/1808.09244)
  - The earlier optical "white hole" claim (Belgiorno et al. 2010) is contested: "the results remain debated, and the experiment's status as a genuine confirmation remains in doubt". **[contested]** — [Wikipedia, Hawking radiation](https://en.wikipedia.org/wiki/Hawking_radiation)
- **Graphene superlattices, Berdyugin et al., Science 375, 430 (2022):** "The criticalities develop upon the velocity of electron flow reaching the Fermi velocity", with "a marked anomaly caused by the Schwinger-like production of hot electron-hole plasma". **[mainstream as a condensed-matter result; "Schwinger" is an analogy]** — [arXiv:2106.12609](https://arxiv.org/abs/2106.12609)
  - Press coverage puts the electron velocity at about c/300. Computed: c/v_F = 300 for v_F = 10⁶ m s⁻¹. — [StatNano/phys.org](https://statnano.com/world-news/90592/Cosmic-physics-mimicked-on-table-top-as-graphene-enables-Schwinger-effect)
- **Ballistic graphene transistors, Schmitt et al., Nat. Phys. 19, 830 (2023), "Mesoscopic Klein-Schwinger effect":**
  - 1D rate w₁d = (2e/h) E ln[1/(1 − e^(−πE_S/E))] (Eq. 1).
  - Pair current I_S = 4V ln[1/(1 − e^(−πV_S/V))] g_sg_v e²/h, with V_S = ΛE_S = ΛΔ_S²/(eħv_F).
  - Differential conductance G_S = 4{(πV_S/V)/(e^(πV_S/V) − 1) + ln[1/(1 − e^(−πV_S/V))]} g_sg_v e²/h (Eq. 3), with affine form G_S ≈ 1.20(V/(0.62V_S) − 1) g_sg_v e²/h.
  - Graphene E_S ≲ 6×10⁷ V m⁻¹ for an effective gap < 0.2 eV.
  - Measured G₀ ≃ 0.18 ± 0.02 mS against predicted 1.20 × 4e²/h = 0.186 mS.
  - Caveats: the gap Δ_S ≃ μ_s is an ansatz with "no full theory", and Zener breakdown at V_Z ≈ 1.8–1.9 V_g bounds the window.
  - The published abstract gives the pinch-off length as Λ ≳ 0.1 μm (the arXiv version says about 1 μm).
  - **[mainstream as an analogue; contested in detail (gap ansatz)]**
  - Computed: E_S = (0.2 eV)²/(eħ·10⁶ m s⁻¹) = 6.08×10⁷ V m⁻¹; G₀ = 0.186 mS.
  - Sources: [arXiv:2207.13400](https://arxiv.org/abs/2207.13400); [NIMS record](https://mdr.nims.go.jp/datasets/e22fb4d3-2924-4ba9-8c5f-729bb274aa8c)

### Inferences
- Rust test vectors: `graphene_es(0.2 eV, 1e6 m/s) = 6.08e7 V/m`; `schwinger_1d_conductance` at V → ∞ approaching the affine slope; `G0 = 1.20*4e²/h = 0.186 mS`. For the BEC: `t_hawking_analogue(kappa=288 s⁻¹) = 0.35 nK`.
- Analogues confirm that horizons plus vacuum fluctuations give thermal pair emission, independently of the microphysics. They say nothing about whether gravity can be engineered: the analogue "metric" is a property of the fluid or medium, and Einstein's equations play no part.

### Gaps
- I did not read Barceló, Liberati & Visser (Living Rev. Relativ. 2011) or Euvé et al. (PRL 2016, water-tank correlations), so their numbers are not included.
- The Weinfurtner 2011 abstract gives no fitted temperature. I did not obtain one.
- Berdyugin et al.: I did not obtain critical current densities or field values (they are not in the abstract).
- I did not locate a published critique of Steinhauer 2019 or Kolobov 2021.

## Q4. Fringe claims: zero-point-field inertia (HRP), Puthoff's vacuum engineering, strong-field "gravity control" (Pais), and the bounds they violate

### Takeaway
Haisch–Rueda–Puthoff (1994) put the observed mass in by hand: m = Γ_zħω_c²/(2πc²) has two free parameters, and with a Planck cutoff the electron needs Γ_z = 1.4×10⁻⁶⁵ s, about 10⁻²² Planck times. The 1998 version replaces these with a free spectral function η(ω). Most hadron mass sits in electrically neutral gluon fields, so an electromagnetic-ZPF origin of inertia fails (Woodward 2001). Puthoff's ZPF-gravity (1989) did not survive analysis by its own co-workers (Cole, Rueda & Danley 2001). Any scheme that lets a Planck-cutoff ZPF gravitate is off by 10¹²¹ from the observed dark-energy density. Strong-field "gravity control" fails on GR: even a field at the Schwinger limit filling a 1 m sphere produces GM/(rc²) ≈ 3×10⁻¹⁹, and QED destroys such a field in about 10⁻¹⁷ s.

### Cited Findings
- **HRP 1994 (PRA 49, 678), "Inertia as a zero-point-field Lorentz force":** "inertia is an electromagnetic resistance arising from the known spectral distortion of the ZPF in accelerated frames". The abstract also offers a "physically rigorous version of Mach's principle". **[refuted, EM-only version]** — [APS abstract](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.49.678)
- **The mass formula as the authors restate it (AIAA 98-3143, 1998):** m_zp = Γ_zħω_c²/(2πc²) (Eq. 5). Γ_z is "a free parameter and Γ_z ≠ Γ_e" (Γ_e = 6.25×10⁻²⁴ s), and ω_c "is also a free parameter". — [Haisch, Rueda & Puthoff, arXiv:physics/9807023](https://arxiv.org/abs/physics/9807023)
  - The authors acknowledge that the 1994 approach "had two weaknesses": complex mathematics and dependence "upon a simplified model" (Planck-oscillator partons).
  - The newer Rueda–Haisch form is m_zp = (V₀/c²)∫η(ω)ρ_zp dω, with ρ_zp = ħω³/(2π²c³) (Eq. 6) and η(ω) a free scattering fraction.
  - The same paper speculates that "The energy released in fusion would be coming from the ZPF". This is a fusion-relevant fringe claim. **[contested; no distinct prediction]**
- Computed: with ω_c = ω_P = √(c⁵/ħG) = 1.855×10⁴³ s⁻¹, Eq. 5 needs Γ_z = 2πmc²/(ħω_c²) = 1.42×10⁻⁶⁵ s for the electron and 2.60×10⁻⁶² s for the proton. The Planck time is 5.39×10⁻⁴⁴ s. The fitted Γ_z simply encodes each mass.
- **Woodward (2001), Found. Phys. 31, 819–835, "Gravity, Inertia, and Quantum Vacuum Zero Point Fields":** hadron inertial masses reside in "electrically chargeless" gluons, so the claim fails for hadrons. Non-gravitational ZPFs cannot cause inertia. **[mainstream critique]** — [DOI 10.1023/A:1017500513005](https://dx.doi.org/10.1023/A%3A1017500513005) (abstract via search results)
  - Reply: Dobyns, Rueda & Haisch, "The Case for Inertia as a Vacuum Effect: A Reply to Woodward and Mahood" (Found. Phys. 2000), points to gluon and W/Z zero-point fields as well. — [arXiv:gr-qc/0002069](https://arxiv.org/abs/gr-qc/0002069)
  - Rueda & Haisch 2005 (Annalen der Physik) still treat only the EM field and extend the idea to passive gravitational mass and the weak equivalence principle. — [arXiv:gr-qc/0504061](https://arxiv.org/abs/gr-qc/0504061)
  - **[contested, generalised version: no quantitative prediction published]**
- **The equivalence-principle bound any mass-from-vacuum model must meet:** MICROSCOPE final result η(Ti,Pt) = [−1.5 ± 2.3(stat) ± 1.5(syst)]×10⁻¹⁵ (PRL 129, 121102, 2022). A material-dependent η(ω) would generically violate this unless it is tuned to be universal. **[mainstream]** — [arXiv:2209.15487](https://arxiv.org/abs/2209.15487)
- **Mass–energy bookkeeping in nuclear reactions is confirmed directly:** Rainville et al., Nature (2005) found 1 − Δmc²/E = −1.4(4.4)×10⁻⁷ ("four-tenths of 1 part in 1 million"). A ZPF reinterpretation of the fusion mass defect therefore has no room to differ numerically. **[mainstream]** — [NIST](https://www.nist.gov/node/399321)
- **Puthoff 1989 (PRA 39, 2333), gravity as a ZPF-induced van der Waals-type force from Zitterbewegung:**
  - Original paper. — [APS](https://link.aps.org/doi/10.1103/PhysRevA.39.2333)
  - Carlip's comment (PRA 47, 3452, 1993) argued that the derivation was flawed.
  - Puthoff's reply (PRA 47, 3454) says Carlip's reanalysis neglects physical constraints, producing "an incorrect 1/r-type force law", and that inverse-square gravity emerges once cutoffs are included. — [APS 47, 3454](https://link.aps.org/doi/10.1103/PhysRevA.47.3454)
  - Cole, Rueda & Danley (PRA 63, 054101, 2001) analysed the van der Waals route and concluded it does not yield Newtonian gravity, "at least not with this particular starting point". — [BU PDF](https://www.bu.edu/simulation/publications/lab/pdf/PRA2001.pdf)
  - **[refuted]**
- **Puthoff "vacuum engineering":**
  - Puthoff, Little & Ibison, "Engineering the Zero-Point Field and Polarizable Vacuum For Interstellar Flight" (JBIS 2002). — [arXiv:1012.5264](https://arxiv.org/abs/1012.5264)
  - DIRD 15 (2010) and the PV model K = exp(2GM/rc²) are covered in the sibling note `dirds_gravity_and_metric_engineering.md`. There the PV metric matches GR at 1PN but differs at 2PN: light bending of 4π(M/b)² versus GR's 15π/4. **[contested → refuted at 2PN if tested; no engineering demonstration]**
- **Cosmological-constant bound (computed):**
  - ZPF energy density to cutoff ω_c is u = ħω_c⁴/(8π²c³). With ω_c = ω_P this is 5.87×10¹¹¹ J m⁻³.
  - The observed dark-energy density is Ω_Λ·3H₀²c²/(8πG) = 5.25×10⁻¹⁰ J m⁻³ (Planck 2018 inputs).
  - The ratio is 10¹²¹. The cutoff that would match observation is ħω ≈ 6.7 meV.
  - Any proposal that treats Planck-cutoff ZPF energy as a gravitating, usable reservoir contradicts the observed expansion rate. **[mainstream: the problem is real; reference Weinberg, RMP 61, 1 (1989), [DOI](https://doi.org/10.1103/RevModPhys.61.1), not re-read]**
- **Pais patents (Navy, 2016–2019):**
  - S_max = f_G(σ²/ε₀)×(velocity term), with numbers reproducible only for surface fields 10²–10⁷ above vacuum breakdown.
  - NAWCAD spent $508k, and the "Pais Effect" "could not be proven".
  - Pais's Eq. 4 sets ε₀E² at the Schwinger field equal to a vacuum-energy flux. This has no basis in Schwinger physics, whose rate goes as exp(−πE_S/E).
  - **[refuted]** — sibling note [`baker_li_hfgw_and_pais_patents.md`](../Gravity%20beams%20fringe%20physics%20and%20DIRDs/baker_li_hfgw_and_pais_patents.md) and its sources (US 10,322,827 B2 etc.)
- **Bounds on strong-field "gravity control" (computed):**
  - A field at E_S stores u = ε₀E_S²/2 = 7.75×10²⁴ J m⁻³, equivalent to 8.6×10⁷ kg m⁻³.
  - A sphere filled with it gives GM/(rc²) = 2.7×10⁻¹⁹ at r = 1 m and 2.7×10⁻¹³ at r = 1 km. GR couples only to stress-energy, and EM field energy is a weak source.
  - QED limit: at E = E_S the pair rate is 4.75×10⁵⁴ m⁻³ s⁻¹. Its rest-mass drain alone is 7.8×10⁴¹ W m⁻³, which would empty the field energy in 1.0×10⁻¹⁷ s, before counting the cascade from acceleration.
  - The best laboratory field (9.1×10¹⁴ V m⁻¹, see Q2) is 1454× below E_S.

### Inferences
- The labels correspond to bounds as follows.
  - (i) HRP inertia, EM-ZPF only: [refuted] by the QCD origin of hadron mass (Woodward's gluon argument) and fine-tuned by construction (free Γ_z, ω_c, η(ω)). A version tuned to universality evades MICROSCOPE only by giving up predictive content.
  - (ii) ZPF gravity (Puthoff 1989): [refuted] by Cole–Rueda–Danley 2001.
  - (iii) "Tap the vacuum for propulsion or energy": contradicts the 10¹²¹ bound, unless only Casimir-type differences are used. Those are tiny, and the sibling DIRD note covers them.
  - (iv) Schwinger-regime gravity control (Pais): [refuted]. The required fields cannot exist for longer than about 10⁻¹⁷ s, and even if they did, their gravity is about 10⁻¹⁹ of c².
- A calculator suite only needs `zpf_energy_density(omega_c)`, `hrp_gamma_z(mass, omega_c)` and `field_curvature(E, r)` to reproduce every bound above.
- These ideas connect to Unruh physics. HRP's "spectral distortion of the ZPF in accelerated frames" is exactly the classical (stochastic electrodynamics) shadow of the Unruh effect; Boyer's quasi-Planckian spectrum is cited in their 1998 paper. Mainstream QFT says that for a uniformly accelerated charge the Unruh bath produces only zero-energy Rindler-photon exchange that reproduces Larmor radiation (Crispino Sec. III.C), not an extra inertial force. So the same physics that underlies Q1 rules out the HRP reading.

### Gaps
- I did not read the full HRP 1994 paper (paywalled) or Carlip's 1993 comment. The Carlip critique is described only through Puthoff's reply abstract and a citing paper.
- I found no quantitative published refutation of the 1998–2005 Rueda–Haisch formulation beyond Woodward's argument. The weak-equivalence-principle argument above is my own inference.
- I did not read Woodward & Mahood (Found. Phys. 29, 899, 1999), only the Dobyns et al. reply.
- I found no other credible claim (beyond Pais) of using Schwinger-regime fields for gravity control.

## Q5. Gravity's role in fusion: stellar confinement as the one working "barrier breaker"

### Takeaway
Gravity does not lower the Coulomb barrier. It confines a hot plasma for millions to billions of years, so tunnelling at the Gamow peak works even though the pp barrier (about 720 keV at 2 fm) is 530 k_BT_c. Simple hydrostatics gives the central conditions to within a factor of about 2. The pp rate at T_c gives about 180 W m⁻³ (computed, unscreened), compared with the model's 276.5 W m⁻³. The Sun's "Lawson product" n_ekT τ ≈ 4×10⁴⁶ keV s m⁻³ is 10²⁵ above D–T ignition, because the first pp step is a weak interaction.

### Cited Findings
- **Solar central values (NASA Sun fact sheet):** mass 1,988,400×10²⁴ kg; radius 695,700 km; luminosity 382.8×10²⁴ J/s; central pressure "2.477 x 10¹¹ bar" (= 2.477×10¹⁶ Pa); central temperature "1.571 x 10⁷ K"; central density "1.622 x 10⁵ kg/m³"; mean density 1408 kg m⁻³. **[mainstream]** — [NASA NSSDC](https://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html)
- **Solar core (Wikipedia):**
  - Central fusion power "is estimated by models to be about 276.5 watts/m³", which is "similar to an active compost heap".
  - Central hydrogen fraction is about 33–34% by mass (sources differ).
  - CNO "generates less than 10% of the total solar energy".
  - The first pp step has a characteristic time of "about one billion years".
  - **[mainstream]** — [Wikipedia, Solar core](https://en.wikipedia.org/wiki/Solar_core)
- **pp S-factor (Solar Fusion II, Adelberger et al., RMP 83, 195, 2011):**
  - S₁₁(0) = (4.01 ± 0.04)×10⁻²² keV b; S′₁₁(0) = (4.49 ± 0.05)×10⁻²⁴ b; S′/S = 11.2 ± 0.1 MeV⁻¹.
  - "Gamow peak energy is ∼6 keV for temperatures characteristic of the Sun's center"; T_c ≈ 15.5×10⁶ K.
  - Solar Fusion III (Acharya et al., RMP 97, 035002, 2025) gives updated recommendations; I did not read its values.
  - **[mainstream]** — [arXiv:1004.2318](https://arxiv.org/abs/1004.2318); [arXiv:2405.06470](https://arxiv.org/abs/2405.06470)
- **D–T ignition for comparison:** the minimum triple product is about 3×10²¹ keV s m⁻³ at T = 14 keV (Wikipedia, after Wesson 2004). — [Wikipedia, Lawson criterion](https://en.wikipedia.org/wiki/Lawson_criterion)
- **Hydrostatic and virial estimates (computed, in the checks script):**
  - Order of magnitude: GM²/R⁴ = 1.13×10¹⁵ Pa. A uniform-density sphere gives P_c = 3GM²/(8πR⁴) = 1.35×10¹⁴ Pa, 184× low.
  - n = 3 polytrope (Eddington standard model), solved by RK4: ξ₁ = 6.897, −ξ₁²θ′(ξ₁) = 2.018, ρ_c/ρ̄ = 54.18. This gives ρ_c = 7.64×10⁴ kg m⁻³ (NASA: 1.622×10⁵) and P_c = GM²/[4π(n+1)θ′(ξ₁)²R⁴] = 1.245×10¹⁶ Pa. The SSM value is 1.99× higher. My guess at the cause is extra central condensation from helium enrichment, which I have not checked.
  - Temperature: T ~ GMμm_u/(k_BR) with μ = 0.61 gives 1.40×10⁷ K. The ideal-gas T from NASA's P_c and ρ_c with μ_c = 0.85 (my assumption, from X_c ≈ 0.34 and Y ≈ 0.64) gives 1.56×10⁷ K, compared with 1.571×10⁷ K published.
  - Kelvin–Helmholtz time GM²/(2RL) = 15.7 Myr (the energy-confinement time scale).
- **Barrier and Gamow window for p + p at T_c = 1.571×10⁷ K (computed):**
  - k_BT_c = 1.354 keV.
  - Coulomb barrier e²/(4πε₀·2 fm) = 720 keV = 532 k_BT_c.
  - Gamow energy E_G = 2m_rc²(παZ₁Z₂)² = 493.1 keV.
  - Gamow peak E₀ = (E_G(k_BT)²/4)^(1/3) = 6.09 keV (matches "∼6 keV"). Width Δ = 4√(E₀k_BT/3) = 6.63 keV; peak suppression exp(−3E₀/k_BT) = 1.4×10⁻⁶.
- **Burning rate (computed; no screening, S(E) = S₀(1 + (S′/S)E), numerical integral):**
  - ⟨σv⟩_pp = √(8/(πμ)) (k_BT)^(−3/2) ∫S(E)exp(−E/k_BT − √(E_G/E))dE = 1.6×10⁻⁴⁹ m³ s⁻¹.
  - With n_p = X_cρ_c/m_p = 3.30×10³¹ m⁻³ (X_c = 0.34): reaction rate n_p²⟨σv⟩/2 = 8.7×10¹³ m⁻³ s⁻¹; proton lifetime 1/(n_p⟨σv⟩) = 6.0 Gyr; power density at 13.1 MeV per pp reaction (ppI: 26.2 MeV deposited per 2 pp) = 183 W m⁻³.
  - This compares with the model's 276.5 W m⁻³ (ratio 0.66). The plausible causes are electron screening (a few %), the X_c and T_c choices (the rate goes roughly as T⁴), the ppII/ppIII branches and CNO.
  - The 6.0 Gyr lifetime conflicts with Wikipedia's "about one billion years". The definitions may differ (per proton versus per pair); unresolved.
- **"Lawson product" of the solar centre (computed):** n_e k_BT_c τ_KH = 4.4×10⁴⁶ keV s m⁻³, about 1.5×10²⁵ times D–T ignition.

### Inferences
- In the Sun, gravity's role in fusion is confinement and compression to n ~ 10³¹–10³² m⁻³ at about 1.35 keV, not barrier lowering. The barrier is "broken" by quantum tunnelling at the Gamow peak (about 10⁻⁶ suppression at 6 keV). The pp rate is then throttled further by the weak interaction. S₁₁ ~ 4×10⁻²² keV b is about 10²⁵–10²⁶ below a strong-interaction S-factor such as d+t (about 2.5×10⁴ keV b; that value is from memory and not sourced here). Gravity can afford this only because τ ~ 10⁷ yr (thermal) to 10¹⁰ yr (nuclear). **[mainstream]**
- Unruh and Hawking temperatures associated with stellar gravity are negligible. The Sun's surface gravity corresponds to T_U = 1.1×10⁻¹⁸ K. Stellar heating is classical virial heating, T ~ GMμm_u/(k_BR) ≈ 1.4×10⁷ K, so there is no "quantum-gravity assist" to fusion anywhere short of black-hole horizons.
- Suggested calculator set with test vectors:
  - `gamow_energy(Z1, Z2, mr)`: 493.1 keV for pp.
  - `gamow_peak(EG, T)`: 6.09 keV at 15.71 MK.
  - `sigma_v(S0, S1, T)`: 1.6×10⁻⁴⁹ m³ s⁻¹.
  - `lane_emden(n)`: n = 3 gives ξ₁ = 6.89685, ρ_c/ρ̄ = 54.18.
  - `kh_time(M, R, L)`: 15.7 Myr.
  - Cross-check against NASA P_c, T_c and ρ_c.

### Gaps
- I did not obtain the S₁₁ value recommended in Solar Fusion III (2025). The numbers above use Solar Fusion II.
- I did not obtain the original SSM (Bahcall/Serenelli) central composition. X_c = 0.34 comes from Wikipedia, and μ_c = 0.85 is my own derivation from it.
- I did not reconcile the 183 W m⁻³ versus 276.5 W m⁻³ gap quantitatively (no screening or branch-ratio model was implemented).
- The source of the "about one billion years" pp time scale (Byrne 1995, via Wikipedia) was not checked.
