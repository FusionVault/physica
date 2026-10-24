# Alternative theories of gravity and inertia, their bounds, and vacuum/"free energy" claims

Labels: **[mainstream]** is established physics or a standard, unrefuted framework. **[contested]** is published and live but disputed, or untested. **[refuted]** means its tested predictions failed, or it conflicts with conservation laws or data.

"Computed" means I recomputed the number in `alternative_gravity_checks.py`, which sits in this folder; run it with `python3 -I alternative_gravity_checks.py`. Its asserts pass. Constants used: CODATA G, ħ and c; GM☉ = 1.32712440018×10²⁰ m³ s⁻²; R☉ = 6.957×10⁸ m; and H0 = 67.4 km/s/Mpc unless stated.

---

## Q1. The PPN framework: γ and β, their best bounds, and how each theory maps onto them

### Takeaway
GR has γ = β = 1. Cassini gives γ − 1 = (2.1 ± 2.3)×10⁻⁵. MESSENGER gives β − 1 = (−1.6 ± 1.8)×10⁻⁵. LLR and MICROSCOPE pin the equivalence-principle and Nordtvedt combinations at 10⁻⁴ to 10⁻¹⁵. Among the theories here, scalar–tensor gravity and massive gravity predict PPN deviations; Brans–Dicke needs ω > 40 000. TeVeS, Einstein–Cartan and emergent gravity have γ = β = 1 or no PPN limit at all, so they must be tested elsewhere.

### Cited Findings
**Defining formulas [mainstream]**
- Perihelion shift, Will's Eq. (65): ω̇ = 42.″98 [ (2 + 2γ − β)/3 + 3×10⁻⁴ (J2/10⁻⁷) ] arcsec per century. The per-orbit form is Δω = (6πm/p)[(2 + 2γ − β)/3 + …]. — [Will 2014, Living Rev. Rel. 17, 4](https://arxiv.org/abs/1403.7377)
- Nordtvedt parameter: ηN = 4β − γ − 3. Light deflection: δθ = ½(1 + γ) × 1.7505″ at the solar limb. — [Will 2014](https://arxiv.org/abs/1403.7377)
- Brans–Dicke: γ = (1 + ω)/(2 + ω) and β = 1. In the Damour–Esposito-Farèse (DEF) parametrisation, ln A(φ) = α0(φ − φ0) + ½β0(φ − φ0)² + …. Brans–Dicke corresponds to α0² = 1/(2ωBD + 3), equivalently 1/(2 + ωBD) = 2α0²/(1 + α0²). "All versions of scalar-tensor gravity predict that γ ≤ 1." — [Will 2014](https://arxiv.org/abs/1403.7377)

**Best bounds [mainstream]**
- Cassini 2002 Doppler (Bertotti, Iess & Tortora 2003): γ − 1 = (2.1 ± 2.3)×10⁻⁵. This used X- and Ka-band links, at closest approach 1.6 R☉, with the spacecraft 8.43 AU from the Sun. "Scalar-tensor theories must have ω > 40,000." — [Will 2014](https://arxiv.org/abs/1403.7377)
- Table 4 of Will 2014 lists the current PPN limits:

  | Parameter | Effect | Limit | Source |
  |---|---|---|---|
  | γ−1 | time delay | 2.3×10⁻⁵ | Cassini |
  | γ−1 | light deflection | 2×10⁻⁴ | VLBI |
  | β−1 | perihelion | 8×10⁻⁵ | J2☉ = (2.2 ± 0.1)×10⁻⁷ |
  | β−1 | Nordtvedt | 2.3×10⁻⁴ | assumes ηN = 4β − γ − 3 |
  | ξ | spin precession | 4×10⁻⁹ | millisecond pulsars |
  | α1 | orbital polarisation | 10⁻⁴ | LLR |
  | α1 | orbital polarisation | 7×10⁻⁵ | PSR J1738+0333 |
  | α2 | spin precession | 2×10⁻⁹ | millisecond pulsars |
  | α3 | pulsar acceleration | 4×10⁻²⁰ | pulsar Ṗ statistics |
  | ζ1 | — | 2×10⁻² | combined |
  | ζ2 | binary acceleration | 4×10⁻⁵ | P̈ of PSR 1913+16 |
  | ζ3 | Newton's third law | 10⁻⁸ | lunar acceleration |

  — [Will 2014](https://arxiv.org/abs/1403.7377)
- Perihelion fits that include MESSENGER, adopting the Cassini γ: β − 1 = (−4.1 ± 7.8)×10⁻⁵. Mars Reconnaissance Orbiter gives β − 1 = (0.4 ± 2.4)×10⁻⁴. LLR gives |ηN| = (4.4 ± 4.5)×10⁻⁴. — [Will 2014](https://arxiv.org/abs/1403.7377)
- MESSENGER ranging (Genova et al. 2018, Nat. Commun. 9, 289), with γ fixed to 1:
  - β − 1 = (−1.6 ± 1.8)×10⁻⁵ (Table 2 gives −1.625×10⁻⁵).
  - η = (−6.6 ± 7.2)×10⁻⁵.
  - J2☉ = (2.246 ± 0.022)×10⁻⁷.
  - Ġ/G is below 4×10⁻¹⁴ per year, with a 2σ uncertainty of 5×10⁻¹⁴ per year.
  - The GR perihelion contribution used is about 42.98″ per Julian century.

  — [Genova et al. 2018 (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC5773540)
- LLR, new APOLLO and OCA infrared data (Biskupek, Müller & Torre):
  - β − 1 = (6.2 ± 7.2)×10⁻⁵ and γ − 1 = (1.7 ± 1.6)×10⁻⁴.
  - Δ(mg/mi)EM = (−2.1 ± 2.4)×10⁻¹⁴.
  - Ġ/G0 = (−5.0 ± 9.6)×10⁻¹⁵ per year.

  — [arXiv:2012.12032](https://arxiv.org/abs/2012.12032)
- The INPOP17a lunar ephemeris, fitted to 48 years of LLR data, finds "no violation of the principle of equivalence at the 10⁻¹⁴ level". — [Viswanathan et al. 2018, MNRAS 476, 1877](https://arxiv.org/abs/1710.09167)
- MICROSCOPE final result for the WEP: η(Ti, Pt) = [−1.5 ± 2.3(stat) ± 1.5(syst)]×10⁻¹⁵ at 1σ. — [Touboul et al. 2022, PRL 129, 121102](https://arxiv.org/abs/2209.15487)
- Hulse–Taylor gravitational-wave damping agrees with GR "to better than half a percent". — [Will 2014](https://arxiv.org/abs/1403.7377)

**How each theory maps onto PPN**
- GR: γ = β = 1, and all other parameters are 0. [mainstream] — [Will 2014](https://arxiv.org/abs/1403.7377)
- Brans–Dicke and scalar–tensor: γ < 1 as above. Binary pulsars additionally bound the DEF β0, which governs spontaneous scalarisation (Will's Fig. 9, after Freire et al. 2012). [mainstream as a framework; pure BD is squeezed to ω > 4×10⁴] — [Will 2014](https://arxiv.org/abs/1403.7377)
- TeVeS "naturally predicts γ = 1 in the post-Newtonian limit, so the Cassini measurements are irrelevant here". Its tightest bounds come from neutron-star–white-dwarf binaries. — [Will 2014](https://arxiv.org/abs/1403.7377). Bekenstein computed β and γ "to agree with solar system measurements". — [Bekenstein 2004, PRD 70, 083509](https://arxiv.org/abs/astro-ph/0403694). In generalised (Æther-type) TeVeS, β, γ and ξ are as in GR. α1 and α2 can satisfy solar-system bounds only "within a modest range of small values" of the couplings. — [Sagi 2009, PRD 80, 044032](https://arxiv.org/abs/0905.4001). [refuted for original TeVeS by GW170817; see Q2]
- Massive gravity (linear Fierz–Pauli): this suffers the van Dam–Veltman–Zakharov (vDVZ) discontinuity, because "the predicted perihelion advance, for example, violates experiment" as m_g → 0. Viable versions need non-linear (Vainshtein) screening. — [Will 2014 §3.3.7](https://arxiv.org/abs/1403.7377); [de Rham 2014, Living Rev. Rel. 17, 7](https://arxiv.org/abs/1401.4173)
- Einstein–Cartan: torsion is algebraically tied to spin density and vanishes outside matter, so the vacuum geometry and the PPN parameters are those of GR. — [Wikipedia: Einstein–Cartan theory](https://en.wikipedia.org/wiki/Einstein%E2%80%93Cartan_theory)

### Inferences
- Computed:
  - Brans–Dicke ω bound. The Cassini 2σ lower edge is γ − 1 = −2.5×10⁻⁵, so ω > 1/(2.5×10⁻⁵) − 2 ≈ 4.0×10⁴. This reproduces Will's "ω > 40,000".
  - In Brans–Dicke at ω = 40 000, ηN = 4β − γ − 3 = 1/(2 + ω) = 2.5×10⁻⁵.
- Computed: 4GM☉/(c²R☉) = 1.751″ with R☉ = 6.957×10⁸ m. Will's 1.7505″ uses a slightly larger R☉. A Rust calculator should reproduce 1.7505″ to 4 significant figures only when it uses Will's radius.
- Computed:
  - The GR perihelion of Mercury comes out at 42.98″ per century, using a = 5.7909050×10¹⁰ m, e = 0.205630 and P = 87.9691 d. This matches Will and Genova.
  - A theory with γ = ½ and β = 1, the effective weak-field value of linear massive gravity without Vainshtein, gives 28.65″ per century. That is about 14″ short.
  - The γ = ½ value comes from the standard vDVZ result, where light bending is ¾ of GR. I did not verify it in a fetched source.
- Computed: the GR round-trip Shapiro delay for the Cassini 2002 geometry is 2(1 + γ)(GM/c³) ln(4 r_E r_C / b²) = 262 µs, with r_E = 1 AU, r_C = 8.43 AU and b = 1.6 R☉.
- From background knowledge, not verified in a fetched text: the DEF formulas are γ − 1 = −2α0²/(1 + α0²) and β − 1 = β0 α0²/[2(1 + α0²)²]. The first is consistent with Will's 1/(2 + ω) = 2α0²/(1 + α0²).

**Test values for a Rust PPN calculator.** "Published" values are quoted from sources; "computed" values come from the script.

| Input | Output | Kind |
|---|---|---|
| γ = 1, β = 1 | Mercury 42.98″/cy | published (Will Eq. 65) and computed |
| ω = 40 000 | γ − 1 = −2.4998×10⁻⁵; ηN = 2.4998×10⁻⁵ | computed |
| γ = 1 | solar-limb deflection 1.7505″ (Will) / 1.751″ (R☉ = 6.957e8) | published / computed |
| Cassini geometry | 262.4 µs round-trip | computed |

### Gaps
- I did not find a post-2014 Living Reviews update by Will. Table 4 above is the 2014 version, and newer β numbers come from Genova 2018 and Biskupek 2020.
- I did not open the INPOP21a or INPOP23 and the latest EPM ephemeris PPN fits. The BepiColombo radio-science PPN results, expected to reach about 10⁻⁶ in γ and β, were not found as of 2026.
- I did not verify that Biskupek et al. was published in a journal. I cite the arXiv version only.

---

## Q2. Scalar–tensor, massive gravity, Einstein–Cartan, MOND/AQUAL/TeVeS and Verlinde: equations, added parameters, distinctive predictions and bounds

### Takeaway
- Brans–Dicke survives only as ω > 4×10⁴ [mainstream, constrained].
- The graviton mass is below about 2×10⁻²³ eV/c² from gravitational waves (GWTC-3/4), and below about 3×10⁻²³ eV/c² from planetary ephemerides [mainstream bounds].
- Einstein–Cartan is untestable with present data: torsion matters only above about 10⁴⁷ g/cm³ [mainstream theory, no signal].
- MOND fits galaxy rotation curves with one constant, a0 ≈ 1.2×10⁻¹⁰ m s⁻². However, original TeVeS and other "dark-matter emulators" were killed by GW170817. The wide-binary test is disputed, and the Solar-System quadrupole is in tension with the radial acceleration relation (RAR) [contested].
- Verlinde's weak-field formula fails Solar-System perihelia by 7 orders of magnitude [refuted in its current form].

### Cited Findings
**Brans–Dicke and scalar–tensor [mainstream; constrained]**
- The equations and the ω > 40 000 bound are given in Q1. In the limit ωBD → ∞ (α0 → 0) the theory becomes indistinguishable from GR. Expected relic deviations in attractor models are |γ − 1| ~ 10⁻⁵ to a few ×10⁻⁷. — [Will 2014](https://arxiv.org/abs/1403.7377)

**Massive gravity [mainstream framework; graviton mass bounded]**
- The dispersion relation is v_g/c ≈ 1 − ½(c/(λ_g f))², with λ_g = h/(m_g c). — [Will 2014](https://arxiv.org/abs/1403.7377)
- LVK GWTC-3, as originally posted (v1, Dec 2021): m_g ≤ 1.27×10⁻²³ eV/c² at 90% credibility. — [arXiv:2112.06861v1](https://arxiv.org/abs/2112.06861v1). **Conflict:** the revised v3, published as Phys. Rev. D 112, 084080 (2025), states m_g ≤ 2.42×10⁻²³ eV/c² at 90%. — [arXiv:2112.06861v3](https://arxiv.org/abs/2112.06861). I did not check the reason for the revision, which loosened the bound. Cite the version explicitly.
- LVK GWTC-4.0 Paper II (March 2026 preprint, 91 events including 42 from O4a): m_g ≤ 1.92×10⁻²³ eV/c² at 90%. — [arXiv:2603.19020](https://arxiv.org/abs/2603.19020)
- Planetary ephemerides, Yukawa potential V = −(GM/r)e^(−r/λ_g):
  - INPOP17b (Bernus et al. 2019, PRL 123, 161103): λ_g > 1.83×10¹³ km, i.e. m_g < 6.76×10⁻²³ eV/c² at 90% CL. The paper also warns that bounds from post-fit residuals "are generically overestimated". — [arXiv:1901.04307](https://arxiv.org/abs/1901.04307)
  - INPOP19a: λ_g ≤ 3.43×10¹³ km is excluded at 99.7% (m_g > 3.62×10⁻²³ eV/c²). The 90% limit is λ_g > 3.93×10¹³ km (m_g < 3.16×10⁻²³ eV/c²). — [Bernus et al. 2020, arXiv:2006.12304](https://arxiv.org/abs/2006.12304)
- Will 2018: from the perihelion advance of Mars (MRO), λ_g > 1.2–2.2×10¹⁴ km, "surpassing the gravitational-wave bound by an order of magnitude". Contrast this with Bernus's caution above. — [Will 2018, CQG 35, 17LT01](https://arxiv.org/abs/1805.10523)
- Review of all the bound types, including Yukawa, dispersion and fifth force, and their theory dependence. — [de Rham et al. 2017, RMP 89, 025004](https://arxiv.org/abs/1606.08462)
- GW170817 gives an independent bound on the graviton mass in bimetric theories. — [Baker et al. 2017, PRL 119, 251301](https://arxiv.org/abs/1710.06394)

**Einstein–Cartan (torsion) [mainstream theory; no measurable deviation]**
- The Cartan equation is T_ab^c + g_a^c T_bd^d − g_b^c T_ad^d = κ σ_ab^c, with κ = 8πG/c⁴. Torsion does not propagate and vanishes outside matter. Big-bounce claims (Popławski) are conjectural. — [Wikipedia: Einstein–Cartan theory](https://en.wikipedia.org/wiki/Einstein%E2%80%93Cartan_theory)
- Spin–torsion effects are comparable to mass effects when n ≃ m/(κħ²). That is about 10⁴⁷ g cm⁻³ for electron-like matter and 10⁵⁴ g cm⁻³ for nucleon-like matter, "never reached in normal situations". — [Ruggiero & Tartaglia 2003, gr-qc/0306029](https://arxiv.org/abs/gr-qc/0306029)
- Lorentz-violation searches constrain 19 of the 24 independent torsion components, down to about 10⁻³¹ GeV. — [Kostelecký, Russell & Tasson 2008, PRL 100, 111102](https://arxiv.org/abs/0712.4393)

**MOND, AQUAL and TeVeS [contested; original TeVeS refuted]**
- Milgrom's law involves an acceleration constant "of the order of the square-root of the cosmological constant in natural units". Relativistic attempts include TeVeS, GEA and BIMOND. — [Famaey & McGaugh 2012, Living Rev. Rel. 15, 10](https://arxiv.org/abs/1112.3960)
- The radial acceleration relation (RAR) uses 2693 points in 153 SPARC galaxies. The fit is g_obs = g_bar / (1 − e^(−√(g_bar/g†))), with g† = 1.20 ± 0.02 (random) ± 0.24 (systematic) ×10⁻¹⁰ m s⁻². — [McGaugh, Lelli & Schombert 2016, PRL 117, 201101](https://arxiv.org/abs/1609.05917)
- Wide binaries, a disputed result:
  - Chae (2023) finds an anomaly below g_N ≈ 10⁻⁹ m s⁻². He reports δ = 0.034 ± 0.007 at 10^−8.91 and δ = 0.109 ± 0.013 at 10^−10.15 m s⁻², and g_obs/g_pred = 1.43 ± 0.06, at 10σ in total. — [arXiv:2305.04613](https://arxiv.org/abs/2305.04613)
  - Contradicted by Banik et al. (2024, MNRAS 527, 4573), using 8611 wide binaries. They find "Newtonian dynamics is preferred at 19σ", with α_grav = −0.021 (+0.065/−0.045), which "excludes MOND at 16σ". — [arXiv:2311.03436](https://arxiv.org/abs/2311.03436)
- Solar-System quadrupole, from Cassini: AQUAL and QUMOND need a much sharper transition than the RAR allows. That is an 8.7σ tension, reduced to 1.9σ if galaxies with bulges are dropped. The quadrupole bound implies no deviation in nearby wide binaries. — [Desmond, Hees & Famaey 2024, MNRAS 530, 1781](https://arxiv.org/abs/2401.04796)
- TeVeS PPN: γ = β = 1 (Q1).
- GW170817 versus TeVeS:
  - "Dark matter emulators", where matter couples to a different metric than gravitational waves, are ruled out. The dark-matter Shapiro delay along the line of sight is estimated at about 400 days, against the 1.7 s observed. The paper gives |γ_GW − γ_EM| < 9.8×10⁻⁸. — [Boran et al. 2018, PRD 97, 041501](https://arxiv.org/abs/1710.06168)
  - The conclusions "extend to … TeVeS and other MOND-like gravities". — [Ezquiaga & Zumalacárregui 2017, PRL 119, 251304](https://arxiv.org/abs/1710.05901)
- AeST (Skordis & Złośnik 2021) is a newer relativistic MOND built to fit the CMB and the linear matter power spectrum, and it is ghost-free at second order. [contested] — [PRL 127, 161302](https://arxiv.org/abs/2007.00082)

**Verlinde's emergent gravity [weak-field formula refuted; overall contested]**
- The apparent dark-matter formula is M_D²(r) = (cH0 r²/(6G)) d(M_b(r) r)/dr, with a0 = cH0. For a point mass this equals MOND with a0 = cH0/6. — [Brouwer et al. 2017, MNRAS 466, 2547](https://arxiv.org/abs/1612.03034)
- Galaxy–galaxy lensing (KiDS × GAMA, 33 613 isolated galaxies) agrees with the emergent-gravity prediction "despite requiring no free parameters". — [Brouwer et al. 2017](https://arxiv.org/abs/1612.03034)
- KiDS-1000 lensing RAR agrees with MOND and emergent-gravity predictions overall. However, there is a ≥6σ difference between early- and late-type galaxies at equal stellar mass, which current modified-gravity theories "would be unable to explain". — [Brouwer et al. 2021, A&A 650, A113](https://arxiv.org/abs/2106.11677)
- RAR test: emergent gravity fits only with substantially lowered stellar M/L ratios. It predicts residuals correlated with radius, which are not observed. — [Lelli, McGaugh & Schombert 2017](https://arxiv.org/abs/1702.04355)
- Solar System: the predicted secular perihelion advances "are discrepant with the data by seven orders of magnitude, ruling out the present emergent gravity weak-field formula with high confidence". — [Hees et al. 2017, PRD 95, 064019](https://arxiv.org/abs/1702.04358)

### Inferences
- Computed λ_g = h/(m_g c), with hc = 1.23984×10⁻⁶ eV m:

  | Bound | λ_g |
  |---|---|
  | 1.27×10⁻²³ eV | 9.76×10¹³ km |
  | 2.42×10⁻²³ eV | 5.12×10¹³ km |
  | 1.92×10⁻²³ eV | 6.46×10¹³ km |
  | 3.16×10⁻²³ eV | 3.92×10¹³ km (matches Bernus's 3.93×10¹³ km; asserted in the script) |
  | 6.76×10⁻²³ eV | 1.83×10¹³ km (matches Bernus 2019) |

  Will's Mars bound of 1.2–2.2×10¹⁴ km corresponds to m_g < 1.03×10⁻²³ to 5.6×10⁻²⁴ eV.
- Computed: for m_g = 1.92×10⁻²³ eV at f = 100 Hz, 1 − v_g/c = 1.08×10⁻²¹. The Yukawa factor at 1 AU is 1 − 2.3×10⁻⁹, which is why ephemerides need about 10 AU baselines.
- Computed MOND numbers:
  - cH0 = 6.55×10⁻¹⁰ m s⁻², so cH0/2π = 1.04×10⁻¹⁰ and cH0/6 = 1.09×10⁻¹⁰. This is the "a0 ~ √Λ" coincidence.
  - The MOND radius of the Sun, √(GM☉/a0), is 7.0×10³ AU. This is why wide binaries at 2–30 kAU are the critical test.
  - The baryonic Tully–Fisher relation (BTFR), V⁴ = G M_b a0, gives 200 km/s for M_b = 10¹¹ M☉.
  - RAR test values, with g† = 1.20×10⁻¹⁰:

    | g_bar | g_obs |
    |---|---|
    | 10⁻¹² | 1.146×10⁻¹¹ |
    | 1.2×10⁻¹⁰ | 1.898×10⁻¹⁰ |
    | 10⁻⁸ | 1.000×10⁻⁸ |

  - The deep-MOND limit √(g_bar a0) is asserted in the script.
- Computed Einstein–Cartan threshold: ρ = m²c⁴/(8πGħ²) gives 3.6×10⁴⁷ g cm⁻³ for electrons and 1.2×10⁵⁴ g cm⁻³ for neutrons. This reproduces Ruggiero & Tartaglia's 10⁴⁷ and 10⁵⁴ when κ = 8πG/c⁴ is used. Neutron-star cores are about 10¹⁵ g cm⁻³, roughly 39 orders lower.
- Theories that touch gravity in a way a gravimeter could see come out as follows:
  - Brans–Dicke and scalar–tensor theories modify G_eff and add a scalar ("breathing") polarisation.
  - Massive gravity adds a Yukawa suppression.
  - MOND matters only below 10⁻¹⁰ m s⁻², i.e. 11 orders below g at Earth's surface, so a terrestrial or GRACE-type gravimeter cannot test MOND directly. The external field of the Galaxy also screens it.

### Gaps
- I did not fetch Verlinde 2016 (arXiv:1611.02269), Bekenstein & Milgrom 1984 (AQUAL: ∇·[μ(|∇Φ|/a0)∇Φ] = 4πGρ, from background knowledge), or the Chae 2024–2025 and Hernandez rebuttals to Banik. The wide-binary dispute was still open in the sources found.
- I did not determine why GWTC-3 v3 loosened the graviton bound from 1.27 to 2.42×10⁻²³ eV.
- Gravity Probe B numbers were not fetched.

---

## Q3. Mach's principle and inertia: Sciama, Woodward, Haisch–Rueda–Puthoff, and the GR analyses

### Takeaway
- Frame-dragging in GR is mainstream and "partly Machian".
- Sciama's inertial induction is a heuristic. It amounts to the cosmological coincidence GM_H/(R_H c²) ~ 1, not a tested theory [contested].
- Woodward's Mach effect, analysed covariantly by Rodal (published in GRG 2019, posted to arXiv December 2025), appears in GR only as a term suppressed by (U_N/c²)(ωL/c)². That is about 10⁻³⁶ for a lab device [GR analysis mainstream; propulsion claim refuted].
- Haisch–Rueda–Puthoff zero-point-field inertia remains unconfirmed and is criticised on basic grounds: hadron mass resides in the electrically neutral gluon sector [contested to refuted].

### Cited Findings
- **Sciama 1953** (MNRAS 113, 34), "inertial induction", adds an acceleration term F = G m_A m_B **a** /(r c²). The Bondi–Samuel list states Mach's principle as Ω = 4πρGT² "is a definite number, of order unity", where T is the Hubble time. [contested] — [Wikipedia: Mach's principle](https://en.wikipedia.org/wiki/Mach%27s_principle)
- **GR and Mach [mainstream]:**
  - Lense–Thirring "certainly satisfies … 'matter there influences inertia here'".
  - At a 1993 Tübingen conference, experts answered "Is general relativity perfectly Machian?" with 3 yes and 22 no.

  — [Wikipedia: Mach's principle](https://en.wikipedia.org/wiki/Mach%27s_principle)
- **Woodward's Mach effect:**
  - Woodward predicted a response ∝ ∂²_t M_loc(t), "amplified by a Sciama-scale cosmic potential Φ/c² ∼ −1".
  - Rodal's analysis:
    - In GR the relevant Landau–Lifshitz terms H^00 ∂²_t H^μν are not independent sources. For lab devices they are suppressed by ~(U_N/c²)(ωL/c)² ≪ 1, "with no enhancement by any cosmological potential".
    - In Hoyle–Narlikar conformal gravity the scalar response is Poisson-like, suppressed by (ωL/c)². Baryon-number conservation forbids an oscillating rest-mass monopole.
    - Conclusion: "any Mach-effect thrust is far too small for practical propulsion".
  - Woodward's reply in GRG disputes Rodal's earlier claim. — [Woodward, GRG, DOI 10.1007/s10714-019-2654-7](https://link.springer.com/article/10.1007/s10714-019-2654-7)
  - Woodward died on 9 August 2025. — [Wikipedia: James Woodward](https://en.wikipedia.org/wiki/James_Woodward_(physicist))

  [Rodal's GR analysis: mainstream. Woodward's propulsion claim: refuted in the sense of no viable mechanism; the experiments are covered in the neighbouring notes.] — [Rodal, arXiv:2512.22687; GRG 51, 64 (2019)](https://arxiv.org/abs/2512.22687)
- **Haisch–Rueda–Puthoff [contested]:**
  - Rueda & Haisch (1998) propose that inertia "may be at least in part" a reaction force from the electromagnetic zero-point field (ZPF). They relate this to Sakharov's induced-gravity conjecture. — [Found. Phys. 28, 1057](https://arxiv.org/abs/physics/9802030)
  - A Foundations of Physics critique argues that hadron inertial mass resides in "chargeless, photon-like gluons", so electromagnetic ZPF cannot cause inertia. — [DOI 10.1023/A:1017500513005](https://dx.doi.org/10.1023/A%3A1017500513005), as summarised by search. I did not read the full text.
  - Carlip (PRA 47, 3452, 1993) reported "a serious computational error" in Puthoff's 1989 zero-point-fluctuation gravity paper. Puthoff replied that other parameter values restore it. — [Wikipedia: Harold E. Puthoff](https://en.wikipedia.org/wiki/Harold_E._Puthoff), a secondary source.

### Inferences
- Computed Sciama's condition with today's numbers: G ρ_c/H0² = 3/(8π) = 0.119, 4πGρ_c T² = 3/2 exactly, and GM_H/(R_H c²) = ½ for the Hubble sphere at critical density. So "Φ/c² ~ −1" holds in order of magnitude in any flat FRW universe. It is a coincidence built into the Friedmann equation, not independent evidence for inertial induction.
- Computed Rodal's suppression factor for a typical Mach-effect thruster stack (M = 0.1 kg, L = 1 cm, f = 50 kHz): (ωL/c)² = 1.1×10⁻¹⁰, and (U_N/c²)(ωL/c)² = 8×10⁻³⁷.
- HRP: the mainstream account is that about 99% of visible-matter mass is QCD binding energy, computed ab initio by lattice QCD. That leaves little room for an electromagnetic-ZPF origin of inertia. This is from my background knowledge; no lattice-QCD source was fetched.

### Gaps
- I did not obtain Sciama's 1953 paper text, Haisch–Rueda–Puthoff 1994 (PRA 49, 678), Carlip 1993, or Levin's critique directly.
- The Woodward thrust experiments, including Tajmar's Dresden tests, are left to the other notes, as assigned.

---

## Q4. Fringe unified theories with gravity-control claims: Heim (Dröscher–Häuser), Kozyrev and Russian "torsion fields", Weber-type electrodynamics (Assis)

### Takeaway
- **Heim theory:** one peer-reviewed paper (Z. Naturforsch. 1977). Its predicted extra neutrinos were excluded at LEP, and its propulsion claims rested on Tajmar's 2006 signals, which Tajmar later attributed to artefacts [refuted].
- **Torsion fields:** state-funded in the USSR from 1986 to 1991, denounced by Academician Aleksandrov in 1991, and opposed by the Academy's pseudoscience commission from 1998. The related gyroscope "antigravity" claims failed three independent 1990 replications [refuted].
- **Weber-type gravity:** any conservative Weber-like gravity force gets either the perihelion or light bending right, but not both (Bunchaft & Carneiro). Weber electrodynamics' predicted electron-mass shift inside a charged shell was not seen [refuted].

### Cited Findings
**Heim / Dröscher–Häuser [refuted]**
- Heim's sole peer-reviewed summary is "Vorschlag eines Weges einer einheitlichen Beschreibung der Elementarteilchen", Z. Naturforsch. 32a (1977). The mass formula was programmed at DESY in 1982, but its derivation was not widely published. The theory predicts two extra neutrinos, which "experiments at LEP have shown … do not exist" (L3 1992).
- Dröscher & Häuser published AIAA papers in 2002 and 2004-3700 ("Guidelines for a space propulsion device based on Heim's quantum theory"). In 2008 the AIAA Nuclear and Future Flight Propulsion Technical Committee published a statement on Extended Heim Theory in Aerospace America. New Scientist (2006) called the physics "obscure and largely unrecognized", "incomplete at best", and not having passed normal peer review.

— [Wikipedia: Heim theory](https://en.wikipedia.org/wiki/Heim_theory)
- Häuser received the AIAA Nuclear and Future Flight Propulsion award for AIAA 2004-3700. — [idw press release](https://idw-online.de/en/news111911)
- Claimed device: a ring rotating above a superconducting coil. About 25 T would be needed to cancel Earth's pull on a 150-tonne craft, and a faster ring in a stronger field would produce repulsive "gravitophotons". This needs "a magnetic coil several metres in diameter capable of sustaining an enormous current density". — [The Register, 6 Jan 2006](https://www.theregister.com/2006/01/06/hyperdrive/)
- Dröscher & Häuser cite "experiments at AIT Seibersdorf since 2006" (Tajmar) as support. — [Bentham chapter](https://benthambooks.com/book/9781608053995/chapter/100225/)
- The empirical anchor:
  - Tajmar's 2006–2007 frame-dragging-like signals reached k ≈ 2.6×10⁻⁵, about 10²¹ above GR.
  - By 2011 an improved setup gave a limit of 4×10⁻¹¹, about 10⁵ below the prediction.
  - Tajmar attributed the earlier signals to acoustic and vibration artefacts.

  — see the neighbouring note [superconductor_gravity_claims_math.md](../Gravity%20beams%20fringe%20physics%20and%20DIRDs/superconductor_gravity_claims_math.md), with sources there.

**Kozyrev and Russian "torsion fields" [refuted]**
- Kozyrev's first "unconventional" work dates from 1958. His detection approach used changes in the conductivity and mechanical properties of materials. — [Kernbach 2013, arXiv:1312.1148](https://arxiv.org/abs/1312.1148). Kernbach is sympathetic to the field; treat this as documentary history.
- The Akimov–Shipov programme claims:
  - spin–spin "torsion" signals at up to about 10⁹ c, carrying information but no energy, "in explicit violation of special relativity";
  - a 1987 Ministry of Defence request, with 500 million roubles reportedly allocated;
  - Ministry of Science funding in 1992–95 and Ministry of Defence funding in 1996–97;
  - a 1994 claim that copper resistivity was cut to 1/80. Independent tests found exposed and unexposed samples identical.
  - A 2008 "torsion engine" on the Yubileiny satellite did not change the orbit "by even a micron" (Kruglyakov).

  — [Wikipedia: Torsion field (pseudoscience)](https://en.wikipedia.org/wiki/Torsion_field_(pseudoscience))
- Chronology from Kernbach:
  - mid-1986: Ryzhkov's order "Take steps to organize the works";
  - 22 Dec 1989: Decree N724 created a Centre of unconventional technologies headed by A. E. Akimov;
  - May 1991: Academician E. B. Aleksandrov's memo to the Supreme Soviet science committee: "defense agencies and KGB funded disjointed pseudoscientific closed developments";
  - 26 Jun 1991: the State Committee for Science and Technology (SCST USSR) created the centre "Vent";
  - 4 Jul 1991: the same SCST issued a resolution on the "vicious practice of financing …".

  — [Kernbach, arXiv:1312.1148](https://arxiv.org/abs/1312.1148)
- The RAS Commission on Pseudoscience was set up in 1998 under the RAS Presidium, at Vitaly Ginzburg's initiative, with Eduard Kruglyakov as first chair. — [Wikipedia: Commission on Pseudoscience](https://en.wikipedia.org/wiki/Commission_on_Pseudoscience)
- In May 1998 Kruglyakov published a critique of torsion "research" in Rossiyskaya Gazeta, as reported by search of [UFN 1999](https://ufn.ru/ufn99/ufn99_3/ufn993m.pdf).
- Kruglyakov's 2001 book criticises Akimov and the funding by the Science and Technology Commission, the Ministry of Defence and the FSB. — [Kowalski, review of Kruglyakov](https://msuweb.montclair.edu/~kowalskil/cf/93kruglyakov.html)
- **Spinning-gyroscope weight claims**, the closest quantitative test of Kozyrev-type "rotation alters weight" ideas:
  - Hayasaka & Takeuchi (PRL 63, 2701, 1989) reported up to about 0.005% weight loss for right-spinning rotors.
  - Faller et al. (PRL 64, 825, 1990) found no anomalous change.
  - Quinn & Picard (Nature 343, 732, 1990) saw about 5% of the required effect, reduced more than tenfold by friction and temperature corrections.
  - Nitschke & Wilmarth (PRL 64, 2115, 1990) set a limit more than 2 orders below the claim.

  — [Nitschke & Wilmarth](https://link.aps.org/doi/10.1103/PhysRevLett.64.2115); [Quinn & Picard](https://link.springer.com/article/10.1038/343732a0); [JILA record for Faller et al.](https://jila.colorado.edu/bibcite/reference/4429); [arXiv:1506.02689](https://arxiv.org/abs/1506.02689)

**Weber-type electrodynamics and gravity (Assis) [refuted]**
- Weber force: F_W = (q1 q2 / 4πε0 r²) r̂ (1 − ṙ²/(2c²) + r r̈/c²).
- A Weber-like force F = k(r̂/r²)(1 − μṙ² + γ r r̈) is conservative if and only if γ = 2μ.
- Gravity versions:

  | Force | μ | γ | Result |
  |---|---|---|---|
  | Tisserand | 1/c² | 2/c² | 3/8 of the then-known anomaly |
  | Assis | 3/c² | 6/c² | correct perihelion, "twice" the light deflection |
  | Raguza | 9/c² | 6/c² | both perihelion and deflection, but non-conservative |

- The perihelion depends only on γ: γ = 6/c² gives it "no matter the value of μ". Light deflection needs 2γ − μ = 3/c². "No conservative gravitational Weber-like force can yield simultaneously the precession of the perihelion of Mercury and the gravitational deflection of light."

— [Bunchaft & Carneiro 1997, gr-qc/9708047](https://arxiv.org/abs/gr-qc/9708047)
- Assis & Graneau (Found. Phys. 1996) derive inertia from Weber gravity with distant matter, Machian "relational mechanics". This "requires a large amount of uniformly distributed matter". — via [PhilPapers](https://philpapers.org/rec/ASSNFO); [Assis, Relational Mechanics](https://www.ifi.unicamp.br/~assis/Relational-Mechanics.pdf)
- Weber's prediction that a charge q inside a shell at potential U changes its inertial mass by qU/(3c²):
  - Mikhailov (1999) claimed a positive result.
  - Junginger & Popovic (Can. J. Phys. 82, 731, 2004) repeated it with optical coupling, found no voltage-dependent frequency shift, and concluded that the force law "is not valid".
  - Lőrincz and Tajmar later showed Mikhailov's experimental concept was flawed. A small group still pursues the idea.

  — [search summary incl. McMaster record](https://experts.mcmaster.ca/scholarly-works/1436968); [Weikert & Tajmar 2019, arXiv:1902.05419](https://arxiv.org/abs/1902.05419)

### Inferences
- Computed: I integrated the Weber-like orbit equation numerically with RK4, using r̈(1 + γGM/r) = rφ̇² − (GM/r²)(1 − μṙ²) and GM/(ac²) = 10⁻⁴. The perihelion advance as a fraction of GR is:

  | μ, γ (in units of 1/c²) | Advance / GR |
  |---|---|
  | Tisserand (1, 2) | 0.3333 |
  | Weber original (½, 1) | 0.1667 |
  | Assis (3, 6) | 0.9998 |
  | Raguza (9, 6) | 0.9998 |
  | (0, 6) | 0.9998 |

  So Δω = (γc²/6) × 6πGM/(c²a(1 − e²)), independent of μ, which confirms Bunchaft & Carneiro. Tisserand's 1/3 of 43″ is 14.3″, i.e. 3/8 of the 38″ known in his day. These are good Rust test values.
- Computed: Weber's electron-mass shift inside a 10 kV shell is eU/(3m_e c²) = 0.65%. That is a large, easily measured effect, so the Junginger–Popovic null is decisive at that level.
- Heim gives no published, independently derivable numerical prediction that has been confirmed. Its one checkable particle-physics prediction, extra neutrinos, is excluded, and its propulsion numbers (25 T, 150 t) have no derivation in a peer-reviewed physics journal.
- The torsion-field claims (10⁹ c, signals carrying information but no energy) contradict special relativity and the Einstein–Cartan result that torsion is confined to matter. Einstein–Cartan torsion is non-propagating and needs about 10⁴⁷ g/cm³ to matter (Q2).

### Gaps
- I could not access the AIAA 2004-3700 full text (paywalled), so the exact Dröscher–Häuser field equations and gravitophoton coupling are unquoted.
- No primary Kozyrev numbers were found: neither his "course of time" c₂ nor his gyroscope weight-change magnitudes. Secondary sources only.
- I found no published text of a 1991 RAS resolution. The 1991 actions documented are the SCST resolutions and Aleksandrov's memo. The 1998 commission is documented.
- I did not access Mikhailov 1999 or the Lőrincz–Tajmar paper directly.

---

## Q5. Vacuum energy and "free energy": Forward 1984, Casimir-cavity claims (Moddel 2021), quantum inequalities, and the Noether and second-law limits

### Takeaway
- Casimir forces are real [mainstream].
- Forward's 1984 "vacuum-fluctuation battery" is a legitimate one-shot device. Re-expanding the plates costs at least the energy gained, so there is no net cycle [mainstream].
- Moddel's group reports pW currents from MIM diodes beside Casimir cavities. No independent replication was found, and a 4-hour output exceeds the cavity's entire Casimir energy by about 10¹¹. Any genuine output must therefore be a continuous energy flow, which in thermal equilibrium the second law forbids [contested].
- Quantum inequalities cap sampled negative energy at about 10⁻²⁶ J m⁻³ for a 1 ns sample.
- The dark-energy density is 5×10⁻¹⁰ J m⁻³ and uniform. Even if it could be "harvested", 1 km³ holds 0.5 J.

### Cited Findings
- **Forward 1984, PRB 30, 1700**, "Extracting electrical energy from the vacuum by cohesion of charged foliated conductors" [mainstream, one-shot]:
  - It proposes a "vacuum-fluctuation battery": a stack of like-charged plates whose electrostatic repulsion nearly balances the Casimir attraction, so the plates move together and the Casimir force adds energy to the electric field.
  - It is recharged by briefly making the electrical repulsion slightly stronger, which re-expands the stack.

  — [APS abstract](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.30.1700)
- **Moddel et al. 2021**, "Optical-cavity-induced current", Symmetry 13, 517 (not PRB or PRR) [contested]:
  - A submicron optical cavity on one side of a metal–insulator–metal (MIM) junction (Ni/NiOx/Al2O3) gives current at zero applied voltage, at the nA level.
  - Current rises as the cavity thins.
  - Maximum power is 1.4 pW in 0.02 µm² for a 33 nm cavity, with P_max = |I_SC V_OC|/4.
  - Eight types of artefact test were negative.
  - A 4-hour continuous run showed no decay. The authors used this to exclude trapped-charge depletion, which could give "20 nA for 3.2 µs".
  - Their interpretation is suppression of vacuum modes upsetting the balance of ZPF-excited electron injection.

  — [arXiv:2101.03085](https://arxiv.org/abs/2101.03085)
- No independent 2023–2026 replication was found. Later descriptions come from Moddel's group and its patents, e.g. US11258379B2. — [Google Patents](https://patents.google.com/patent/US11258379B2/en)
- **Moddel & Dmitriyeva 2019** (Atoms 7, 51):
  - "None of the proposals have been reliably demonstrated."
  - Nonlinear processing of the ZPF and mechanical extraction with Casimir cavities "are shown to violate thermodynamics principles".
  - The SED atom-pumping route does not, "but may face other obstacles". Initial results are "inconclusive".

  — [arXiv:0910.5893](https://arxiv.org/abs/0910.5893)
- **Quantum inequalities** (Ford & Roman) [mainstream]:
  - For a massless scalar field in 4D Minkowski space, with Lorentzian sampling time t0, ρ̂ ≥ −3/(32π² t0⁴) in units ħ = c = 1. The same bound holds for the flux.
  - These are "of the form required to prevent large-scale violations of the second law of thermodynamics".
  - Any inertial observer "cannot see an arbitrarily large negative energy density which lasts for an arbitrarily long period of time".

  — [Ford & Roman 1995, PRD 51, 4277](https://arxiv.org/abs/gr-qc/9410043); [Ford & Roman 1997, PRD 55, 2082](https://arxiv.org/abs/gr-qc/9607003)
- **Work-extraction limit** [mainstream]: the maximum work extractable by cyclic unitary operations, the "ergotropy", is fixed by the density matrix and the Hamiltonian. — [Allahverdyan, Balian & Nieuwenhuizen 2004](https://arxiv.org/abs/cond-mat/0401574)

### Inferences
- **Casimir arithmetic** (computed, ideal plates): |E|/A = π²ħc/(720 d³) and P = π²ħc/(240 d⁴).

  | d | Energy per area | Pressure |
  |---|---|---|
  | 1 µm | 4.3×10⁻¹⁰ J/m² | 1.3×10⁻³ Pa |
  | 100 nm | 4.3×10⁻⁷ J/m² | 13 Pa |
  | 10 nm | 4.3×10⁻⁴ J/m² | 1.3×10⁵ Pa (about 1.3 atm) |

  Collapsing a 1 m² plate pair from 1 µm to 10 nm yields at most 4.3×10⁻⁴ J, once. Real metals give less below the plasma wavelength (Lifshitz theory). This is a known limit of the ideal formula; I did not compute it.
- **Moddel 2021 energy budget** (computed):
  - Power density: 1.4 pW / 0.02 µm² = 70 W/m².
  - 4 hours at 1.4 pW = 2.0×10⁻⁸ J. The ideal Casimir energy of a 33 nm × 0.02 µm² cavity is 2.4×10⁻¹⁹ J, a ratio of 8×10¹⁰.
  - So the claimed output cannot be released Casimir energy. It must be a steady power flow. If the device sits in thermal equilibrium, a DC output from a single reservoir would violate the second law and detailed balance, which is Moddel & Dmitriyeva's own criterion.
  - The physically allowed explanations are an unrecognised non-equilibrium source (thermal gradients, light, electrochemistry, or instrument offsets) or a genuinely new non-equilibrium effect. Either way, it is not "free" vacuum energy.
- **Quantum-inequality magnitudes** (computed, SI form ρ ≥ −3ħ/(32π² c³ t0⁴)):

  | t0 | Bound |
  |---|---|
  | 1 ns | −3.7×10⁻²⁶ J m⁻³ |
  | 1 µs | −3.7×10⁻³⁸ J m⁻³ |
  | 1 s | −3.7×10⁻⁶² J m⁻³ |

- **Dark energy** (computed from Planck 2018: Ω_Λ = 0.6847, H0 = 67.36): ρ_Λc² = 5.2×10⁻¹⁰ J m⁻³, which is 0.52 J per km³ and 5.7×10¹¹ J in Earth's volume.
  - It is the vacuum's ground state and spatially uniform, with no gradient to run an engine.
  - By the ergotropy result, a passive state, the ground state included, yields zero work under any cyclic unitary process.
  - Planck parameters are from background knowledge of [Planck 2018, arXiv:1807.06209](https://arxiv.org/abs/1807.06209); I did not refetch them.
- **Gravity "free energy"** (computed): lifting 1 kg by 1 m costs 9.81 J, and the round trip in a static field returns exactly 0.
  - The most any mass can give by falling from infinity to Earth's surface is GM/R = 62.6 MJ/kg, once.
  - Energy conservation follows from time-translation symmetry (Noether). A static spacetime has a timelike Killing vector and a conserved energy, so a closed cycle in a static gravitational field does zero net work. Textbook; no source fetched.
  - Devices claiming net output from gravity must therefore be drawing on a time-dependent source (tides, Earth's rotation) or be in error.

### Gaps
- I could not read the Forward 1984 full text (APS paywall), so his numerical example of energy per stack is not quoted.
- The Moddel & Dmitriyeva 2019 full text returned 403 from MDPI, so I used only the arXiv abstract.
- I found no published critique specifically of the 2021 MIM result, and no independent replication.

---

## Q6. The speed of gravity and GW polarisations as tests

### Takeaway
GW170817 with GRB 170817A fixes −3×10⁻¹⁵ ≤ (v_gw − c)/c ≤ +7×10⁻¹⁶ [mainstream]. That eliminated most late-universe Horndeski, Galileon, Einstein-Æther, Hořava and TeVeS-type theories. Network polarisation tests favour pure tensor modes, by log₁₀ Bayes factors of about 21 against vector and 23 against scalar for GW170817 [mainstream]. GWTC-3 and GWTC-4.0 find no non-GR polarisation content.

### Cited Findings
- **Speed bound:** GW170817 constrains the difference between the speed of gravity and the speed of light to "between −3×10⁻¹⁵ and +7×10⁻¹⁶ times the speed of light". It also gives new Lorentz-invariance bounds and a Shapiro-delay equivalence-principle test. — [LVK, Fermi & INTEGRAL 2017, ApJL 848, L13](https://arxiv.org/abs/1710.05834)
- Other statements of the bound and what it eliminates:
  - |c_g/c − 1| ≤ 5×10⁻¹⁶. This eliminates the covariant Galileon and most Horndeski and beyond-Horndeski dark-energy models, and extends to Einstein-Æther, Hořava, generalised Proca, TeVeS and MOND-like theories. — [Ezquiaga & Zumalacárregui 2017](https://arxiv.org/abs/1710.05901)
  - Deviations are below a few ×10⁻¹⁵. Of the quartic and quintic beyond-Horndeski functions, only a conformal coupling survives. — [Creminelli & Vernizzi 2017, PRL 119, 251302](https://arxiv.org/abs/1710.05877)
  - Stringent constraints on general scalar–tensor and vector–tensor theories. — [Baker et al. 2017](https://arxiv.org/abs/1710.06394)
  - |γ_GW − γ_EM| < 9.8×10⁻⁸ from the dark-matter Shapiro delay. — [Boran et al. 2018](https://arxiv.org/abs/1710.06168)
- **Polarisations:** generic metric theories allow up to six polarisations.
  - GW170814 was the first three-detector event. Bayes factors were "more than 200 and 1000 in favor of the purely tensor polarization against purely vector and purely scalar, respectively". — [LVK 2017, PRL 119, 141101](https://arxiv.org/abs/1709.09660)
  - GW170817, with the sky position fixed to NGC 4993, gave log₁₀ Bayes factors of +20.81 ± 0.08 (tensor vs vector) and +23.09 ± 0.08 (tensor vs scalar). — [LVK, "Tests of GR with GW170817", arXiv:1811.00364](https://arxiv.org/abs/1811.00364)
  - GWTC-3: "no evidence for dispersion of GWs, non-GR modes of polarization, or post-merger echoes". — [arXiv:2112.06861](https://arxiv.org/abs/2112.06861)
  - GWTC-4.0 (91 events): "no evidence … for non-GR modes of polarization". — [arXiv:2603.19019](https://arxiv.org/abs/2603.19019)
- A massive spin-2 graviton has "3 additional helicity states", which is the origin of vDVZ. — [Will 2014](https://arxiv.org/abs/1403.7377)

### Inferences
- Computed:
  - Upper bound: Δv/v ≈ cΔt/D with Δt = 1.74 s and D = 26 Mpc (the conservative distance) gives +6.5×10⁻¹⁶, which rounds to the published +7×10⁻¹⁶.
  - Lower bound: allowing the gravitational wave to be emitted up to 10 s before the gamma rays gives (1.74 − 10) s × c/D = −3.1×10⁻¹⁵, the published −3×10⁻¹⁵.
  - These are reproducible Rust test values.
- How the theories here map onto the polarisation and speed tests:
  - Brans–Dicke and scalar–tensor theories predict an extra scalar breathing mode. They are not excluded by pure-mode tests while the scalar amplitude is small, consistent with ω > 4×10⁴.
  - TeVeS-type "dark-matter emulators" fail on speed and Shapiro delay.
  - Einstein–Cartan, with non-propagating torsion, predicts GR waves.
  - Heim, torsion-field and Weber-type theories make no quantitative gravitational-wave predictions that can be confronted with the data.

### Gaps
- I did not extract the GW170817 Shapiro-delay equivalence-principle numbers from ApJL 848 L13 (−2.6×10⁻⁷ ≤ γ_GW − γ_EM ≤ 1.2×10⁻⁶ is my recollection, unverified). The Boran et al. figure above is the verified one.
- I did not fetch the mixed-polarisation, as opposed to pure-mode, bounds from GWTC-3 or GWTC-4.0.
