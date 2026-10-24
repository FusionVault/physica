# Replication erases gravity's anomalies, except the statistical ones

As of October 2026, no measurement contradicts general relativity (GR) or Newtonian gravity in any regime that has been tested. The anomalies that once looked like openings all follow one pattern. Those tied to a single apparatus or spacecraft shrink to systematics once someone independently redoes the work. Those that stay open are statistical tests on populations, or decades-long geodetic series:

- NIST rebuilt the BIPM torsion balance and got a G **250–282 ppm lower** than BIPM had.
- The Pioneer acceleration turned out to be thermal recoil.
- ESA reproduced the 2005 Rosetta flyby anomaly and ruled out every gravitational cause.
- Lunar laser ranging (LLR) puts gravitational shielding about **10⁹ times below** Majorana's 1920 claim.

What is genuinely unsettled is a short list. It has three parts.

- **Open questions in data:** the 551 ppm spread in measured G (a metrology problem), the Gaia wide-binary test of MOND (Chae against Banik, open but decidable), the lunar eccentricity rate, the size of antimatter's gravity at the percent level, and the 1989 Tate Cooper-pair mass, which nobody has repeated.
- **GR predictions not yet observed:** gravitational-wave memory, primordial waves, scalar or vector admixtures in wave polarisations, entanglement induced by gravity, and a single graviton.
- **Alternative theories confined to tight boxes:** the post-Newtonian (PPN) parameters γ and β lie within about 2×10⁻⁵ of 1; waves travel at c to 10⁻¹⁵; the graviton is lighter than 1.92×10⁻²³ eV/c² (from waves) or 1.0×10⁻²⁴ eV/c² (from ephemerides); and gravitational-strength Yukawa forces must have ranges below 38.6 µm.

The gravity people actually manipulate is Newtonian and small. LIGO's rotor pushes a 40 kg mirror with a **19 pN** force predicted to better than 1%, and a 92 mg sphere is the smallest source whose gravity has been measured. GR caps everything beyond that: the strongest feasible laboratory rotor drags frames at about 5×10⁻²¹ rad/s and radiates 10⁻²⁹ W. The fringe theories (Heim, Russian torsion fields, Weber gravity) and the vacuum-energy devices each fail at one step that a calculator can encode. About sixty calculators from the notes, each with a test vector, are ready for `gwsens`, `geo5d` and `gravitics`. The strongest end-to-end tests are the CODATA G fit (χ² = 195.7 reproduced exactly), the Anderson flyby table, the LIGO calibrator force and the decay rates of two binary pulsars.

*Labels and numbers.* Every item carries a label:

- **[mainstream]**: accepted physics, or a collaboration measurement.
- **[contested]**: published and under live dispute, preprint-only, or model-dependent.
- **[refuted]**: tested and failed, or in conflict with conservation laws or data.

Every item also says whether it is **open** or **resolved**. Values marked *(computed)* were recomputed with the four stdlib-only scripts beside the notes, all run with `python3 -I`; every one runs cleanly:

- [anomalies_checks.py](notes/anomalies_checks.py)
- [alternative_gravity_checks.py](notes/alternative_gravity_checks.py)
- [lab_gravity_checks.py](notes/lab_gravity_checks.py)
- [deep_gravitational_wave_checks.py](notes/deep_gravitational_wave_checks.py)

Computed values are our inferences, not published claims. Three further numbers were checked while writing this report and are marked *(checked here)*. One of them corrects the notes.

## Single-apparatus anomalies shrink when a second lab rebuilds the apparatus

Eight anomaly families were examined. The ones that dissolved had something in common: one instrument, one analysis or one spacecraft, and then a test that a constant new force could not pass, whether a replication, an out-of-sample prediction or a time dependence.

### G spans 551 ppm, and NIST's rebuild of the BIPM balance points at systematics

G is the worst-known of the major constants. CODATA 2022 received no new input after 2018, so its value is identical to 2018's: **G = 6.674 30(15)×10⁻¹¹ m³ kg⁻¹ s⁻² (22 ppm)**. CODATA reached that value only after multiplying every input uncertainty by 3.9 ([CODATA 2022](https://arxiv.org/abs/2409.03787)). The 16 inputs span **551 ppm** *(computed)*, although individual experiments claim 12–20 ppm. The methods include fibre and strip torsion balances, suspended bodies, a stationary-body weight change and an atom-interferometer gradiometer.

A generalised least-squares fit with the three published correlations reproduces CODATA exactly *(computed)*:

- mean 6.674300, with an unexpanded uncertainty of 5.7 ppm;
- **χ² = 195.7 for 15 degrees of freedom**, a Birge ratio of 3.61, which falls to CODATA's own 12.9 after the ×3.9 expansion;
- normalised residuals of +7.75 for BIPM-14, −6.80 for JILA-18 and +4.78 for BIPM-01, CODATA's three worst ([CODATA 2022 PDF](https://arxiv.org/pdf/2409.03787)).

The scatter shows up even inside single groups. HUST's two 2018 methods differ by 45 ppm, or 2.7σ ([Wikipedia: Gravitational constant](https://en.wikipedia.org/wiki/Gravitational_constant)). BIPM's 2013 value moved by 13.5 ppm in an erratum that found, among other errors, a density-gradient correction applied twice ([PRL 113, 039901](https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevLett.113.039901/fulltext)). Cold-atom interferometry gives 6.67191(99), which is 358 ppm below CODATA ([Rosi et al. 2014](https://arxiv.org/abs/1412.7954)). A Bayesian reanalysis that models unknown systematics recommends **6.6740 ± 0.0015**, an uncertainty ten times CODATA's ([Rinaldi et al.](https://arxiv.org/pdf/2209.07416)).

The decisive event of 2026 was NIST operating the BIPM apparatus itself, with new mass metrology and a new environment. The two published accounts give different numbers:

| Source | G (10⁻¹¹ m³ kg⁻¹ s⁻²) | Uncertainty | vs BIPM-14 *(computed)* | vs CODATA *(computed)* |
|---|---|---|---|---|
| Metrologia abstract ([NIST](https://www.nist.gov/node/1911461)) | 6.67366(20) | ±0.00020 | **−282 ppm, −7.3σ** | −96 ppm |
| APS news ([Physics 19, 64](https://physics.aps.org/articles/v19/64)) | 6.67387 | 57 ppm | **−250 ppm, −4.0σ** | −64 ppm |

The APS account adds that the NIST value sits within 1σ of the world average. The NIST page also states u_r = 3.0×10⁻⁶, which cannot match ±0.00020 (that is 3.0×10⁻⁵ *(computed)*), so the page has a typo. Which number is final is unresolved, because the notes could not open the Metrologia full text. No single cause was identified. An earlier conference version had found a gas-pressure effect in NIST's own earlier runs ([NIST record](https://www.nist.gov/publications/redetermination-gravitational-constant-bipm-torsion-balance-nist)), and a 2023 test ruled out stray AC magnetic fields ([CODATA 2022 PDF](https://arxiv.org/pdf/2409.03787)).

The inference is direct. The same apparatus in a different lab gave a different answer, so a systematic tied to the apparatus or its site is the economical explanation. New physics would have to act differently on one instrument in two places. **The spread is [mainstream] and open; its cause is unidentified.**

The claimed 5.9-year oscillation of G in step with the length of day is **[refuted] and resolved**:

- a recompilation of 35 years of data significantly weakened the correlation ([Schlamminger, Gundlach & Newman](https://arxiv.org/abs/1505.01774v1));
- a Bayesian model comparison favours a constant plus unknown noise by a Bayes factor of **at least e³⁰** ([Pitkin](https://arxiv.org/abs/1505.06725v1));
- a frequentist reanalysis agrees ([Desai](https://arxiv.org/pdf/1607.03845)).

Secular drift is bounded independently:

- LLR gives Ġ/G = (7.1 ± 7.6)×10⁻¹⁴ yr⁻¹ ([Belgacem et al.](https://arxiv.org/pdf/1812.11181));
- newer APOLLO and OCA laser ranging gives (−5.0 ± 9.6)×10⁻¹⁵ yr⁻¹ ([Biskupek et al.](https://arxiv.org/abs/2012.12032));
- MESSENGER gives below 4×10⁻¹⁴ yr⁻¹ ([Genova et al. 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC5773540)).

None of these drift limits constrains a periodic signal, so the rebuttal of the length-of-day claim rests on statistics, not on LLR.

### The flyby formula predicted Juno at +6 mm/s and got zero

Six Earth flybys between 1990 and 2005 showed millimetre-per-second changes in asymptotic speed. Anderson et al. fitted them with **ΔV∞/V∞ = K(cos δᵢ − cos δₒ)**, where δᵢ and δₒ are the incoming and outgoing declinations and **K = 2ω⊕R⊕/c = 3.099×10⁻⁶**. The authors themselves call this a phenomenological fit, not an explanation ([Nieto & Anderson](https://arxiv.org/pdf/0910.1321)).

The formula reproduces the paper's predicted column to two decimals *(computed)*, against the observed values in Table I ([via arXiv:2406.06604](https://arxiv.org/pdf/2406.06604)):

| Flyby | Predicted ΔV∞ *(computed)* | Observed ΔV∞ |
|---|---|---|
| Galileo I | +4.123 mm/s | |
| Galileo II | −4.674 mm/s | |
| NEAR | +13.279 mm/s | 13.46 ± 0.01 mm/s |
| Cassini | −1.068 mm/s | |
| Rosetta I | +2.067 mm/s | |
| MESSENGER | +0.055 mm/s | |

Secondary sources disagree on the uncertainties. Galileo I is listed as ±0.3 or ±0.08 mm/s, and Nieto and Anderson's own text gives (4.13 ± 0.03) mm/s ([Wikipedia: Flyby anomaly](https://en.wikipedia.org/wiki/Flyby_anomaly); [Páramos & Hechenblaikner](https://ar5iv.arxiv.org/html/1210.7333)).

Every out-of-sample test failed:

- Rosetta's 2007 and 2009 flybys gave ΔV∞ ≈ 0 ([Páramos & Hechenblaikner](https://ar5iv.arxiv.org/html/1210.7333)).
- Juno in 2013 was predicted to show **about +6 mm/s** ([Iorio](https://arxiv.org/abs/1311.4218)). It showed **0 ± 0.8 mm/s**.
- BepiColombo in 2020 showed no unexpected acceleration.

The Juno and BepiColombo results are reported via [Wikipedia: Flyby anomaly](https://en.wikipedia.org/wiki/Flyby_anomaly), citing the primary papers. An independent fit also finds that the formula cannot accommodate these nulls ([Acedo 2017](https://ar5iv.labs.arxiv.org/html/1711.02875)).

In 2026, ESA's new deep-space navigation system reproduced the 2005 Rosetta anomaly. It ruled out ephemeris errors, gravity-field modelling, relativistic effects and propagator settings, and found solar and albedo radiation pressure insufficient. It proposes that albedo illumination of thermally sensitive surfaces triggered outgassing, and a data package is planned ([ISSFD 2026](https://www.issfd2026.fr/data/onglet12/module0/modalPreview.php?langue=fr&paramProjet=20840)). The status labels are:

- the tracking residuals: **[mainstream]**;
- the universal formula: **[contested]**, and it failed its predictions;
- the anomaly itself: **open**, now treated as specific to each spacecraft.

Once ESA releases the package, anyone with an orbit-determination code can reanalyse it.

### Shielding is excluded ten orders below Majorana, so eclipse dips are environmental

Majorana's absorption law is F = F_N·exp(−h∫ρ dl), where h is an extinction coefficient in m²/kg ([Wikipedia: Gravitational shielding](https://en.wikipedia.org/wiki/Gravitational_shielding)). His experiments implied **h ≈ 6.7×10⁻¹³ m²/kg**, taken from a search snippet of the source ([von Borzeszkowski & Chrobok](https://arxiv.org/pdf/gr-qc/0310081)). With assumed standard radii and densities, the diametral column densities are 7.026×10¹⁰ kg/m² through Earth and 1.162×10¹⁰ kg/m² through the Moon. Majorana's h would then remove **4.6% of gravity through Earth and 0.78% through the Moon** *(computed)*.

The bounds on h (m²/kg) are ([Wikipedia: Gravitational shielding](https://en.wikipedia.org/wiki/Gravitational_shielding); consistent in cgs with [hep-ph/0204284](https://arxiv.org/pdf/hep-ph/0204284)):

| Source | Bound |
|---|---|
| Poincaré, lunar motion | 10⁻¹⁸ |
| Zürich G experiment (best laboratory bound) | 4.3×10⁻¹⁵ |
| Eckhardt, LLR | **< 10⁻²²** |
| Williams, Turyshev & Boggs, LLR | **(3 ± 5)×10⁻²²** |

Majorana shielding is therefore **[refuted] and resolved**, by a factor of 10⁹–10¹⁰.

The eclipse anomalies follow from that. Allais reported anomalous pendulum precession during the 1954 and 1959 eclipses. Later attempts were mixed, and an automated Foucault pendulum in 2010 saw no change above 0.3°/h ([Wikipedia: Allais effect](https://en.wikipedia.org/wiki/Allais_effect)). A LaCoste-Romberg gravimeter at Mohe saw near-symmetric dips of **(7.0 ± 2.7)×10⁻⁸ m/s²** around the 1997 eclipse ([Wang et al.](https://arxiv.org/abs/1003.4947)). Three later studies undercut that result:

- the same eclipse data were judged not to support shielding ([Wikipedia: Gravitational shielding](https://en.wikipedia.org/wiki/Gravitational_shielding));
- cooling of the upper atmosphere produces a moving air-mass excess of about 7×10⁻⁸ m/s², the right size ([Van Flandern & Yang, summary](https://arxiv.org/pdf/physics/0610197v1));
- a superconducting gravimeter during the 2009 eclipse saw no anomaly ([Cui, Fang & Wang](https://ch.whu.edu.cn/en/article/id/700)).

Two numbers settle it *(computed)*. Read as the Moon shielding the Sun's 5.930×10⁻³ m/s², the 1997 dip needs **h = 1.0×10⁻¹⁵ m²/kg, about 10⁷ times the LLR bound**. The largest eclipse signal LLR allows is **6.9×10⁻¹⁵ m/s², or 6.9×10⁻⁷ µGal**.

Wikipedia's 6×10⁻¹⁹ "constraint" from the same eclipse data uses different assumptions from this explanation-style estimate, and the notes did not resolve its derivation. The Allais effect and the eclipse dips are **[contested] and open but unreplicated**. For this project the result is a practical rule: a µGal eclipse signal in any gravimeter series is atmospheric, thermal or human, never absorption.

### Pioneer shows how anomalies die: through their time dependence

Pioneer 10 and 11 showed a constant-looking sunward acceleration of **a_P = (8.74 ± 1.33)×10⁻¹⁰ m/s²** that years of searching could not attribute ([Anderson et al. 2002](https://arxiv.org/abs/gr-qc/0104064)). A finite-element thermal model driven by flight telemetry reproduced its size, direction and slow decay. Once thermal recoil is included, no anomalous acceleration remains ([Turyshev et al. 2012](https://arxiv.org/abs/1204.2507)). The decaying model a(t) = a₀·2^(−(t−t₀)/τ) gives two fits ([ar5iv](https://ar5iv.labs.arxiv.org/html/1204.2507)):

| Fit | τ | a₀ |
|---|---|---|
| Thermal model | 36.9 ± 6.7 yr | (7.4 ± 2.5)×10⁻¹⁰ m/s² |
| Doppler data | 28.8 ± 2.0 yr | (10.1 ± 1.0)×10⁻¹⁰ m/s² |

The 1σ error ellipses overlap, and the thermal error budget of 26.3% is dominated by degradation of the radioisotope generators' surfaces.

The decay is the fingerprint. Both values of τ are far shorter than plutonium-238's 87.74-year half-life, consistent with output falling faster than the isotope decays, as degrading thermocouples would cause. A cosmological constant acceleration cannot decay at all. The cH₀ coincidence that fuelled cosmological readings is a_P/cH₀ = 1.29 for H₀ = 70. Only **63 W** of directed thermal emission would supply a_P on an assumed 241 kg craft *(computed)*.

New-physics readings of Pioneer are **[refuted]**, the thermal explanation is **[mainstream]**, and the case is **resolved**. The same method, a time-resolved model specific to the spacecraft, is what ESA has now applied to Rosetta.

### Wide binaries are MOND's one live test, and contamination will decide it

The radial acceleration relation (RAR) is a genuine empirical law. Across **2693 points in 153 galaxies**, the observed acceleration is a single function of the baryonic one: g_obs = g_bar/(1 − e^(−√(g_bar/g†))), with **g† = 1.20 ± 0.02 (random) ± 0.24 (systematic) ×10⁻¹⁰ m/s²**. The residuals are 0.11 dex and consistent with negligible intrinsic scatter ([McGaugh, Lelli & Schombert 2016](https://arxiv.org/abs/1609.05917); [ar5iv](https://ar5iv.labs.arxiv.org/html/1609.05917)).

Test values *(computed)*:

- g_bar = 10⁻¹² m/s² gives 1.1462×10⁻¹¹ m/s², a boost of 11.46; 10⁻¹⁰ gives 1.6705×10⁻¹⁰; 10⁻⁸ gives 1.0001×10⁻⁸. The low-acceleration limit tends to √(g_bar·g†).
- The baryonic Tully–Fisher relation v⁴ = GMg† gives **199.8 km/s** for 10¹¹ M☉.
- Newtonian acceleration equals g† at 10.8 kpc for 10¹¹ M☉, and at **9,942 au** for a binary of 2 M☉ in total.
- cH₀/2π is 1.08×10⁻¹⁰ m/s² for H₀ = 70, close to g†.

The RAR as an observation is **[mainstream]**. MOND as its explanation is **[contested]**. Neither can be reached by a gravimeter: g† is 11 orders below g at Earth's surface, and the Galaxy's external field screens it anyway.

The wide-binary test in Gaia DR3 is **[contested] and open**:

- **Chae (2023)**, 26,615 binaries within 200 pc: a gravity boost of **1.43 ± 0.06** at g_N ≈ 10^−10.15 m/s², with 10σ combined significance ([arXiv:2305.04613](https://arxiv.org/abs/2305.04613)).
- **Banik et al. (2024)**, 8611 binaries at 2–30 kAU: α_grav = −0.021 (+0.065/−0.045). They exclude MOND at 16σ and prefer Newton at **19σ**. MOND predicted about a 20% velocity boost, which is a 1.44 gravity boost *(computed)* ([arXiv:2311.03436](https://arxiv.org/abs/2311.03436)).
- **Chae (2025)**, Bayesian 3D orbit modelling of 312 binaries with radial velocities: Γ = 0.000 ± 0.011 at high acceleration and **Γ = 0.085 ± 0.040 (γ_g = 1.48)** in the MOND regime ([arXiv:2502.09373](https://arxiv.org/abs/2502.09373)).
- **Chae et al. (2026)**, 36 binaries with radial-velocity errors below 100 m/s: **γ = 1.60 (+0.24/−0.14)**, above MOND's own prediction. The result is dominated by two systems whose 3D velocities exceed Newtonian escape ([arXiv:2601.21728](https://arxiv.org/abs/2601.21728)).
- **Solar-System quadrupole**, measured by Cassini: AQUAL and QUMOND need a far sharper transition than the RAR allows. The tension is 8.7σ, or 1.9σ without bulge galaxies, and the bound implies no deviation in nearby wide binaries ([Desmond, Hees & Famaey 2024](https://arxiv.org/abs/2401.04796)).

Commentary frames the dispute as sample selection and velocity cuts ([Physics World](https://physicsworld.com/?p=109496)).

Both camps agree on the prediction, a boost of about 1.4 at around 10 kAU, so the argument is purely about data quality. Velocities above escape speed are the expected signature of unresolved companions or projection effects. Spectroscopic follow-up of those two systems, and Gaia DR4, can settle it.

### Antimatter falls down, and composition never mattered

ALPHA-g found that antihydrogen accelerates toward Earth at **a_g = [0.75 ± 0.13 (stat+syst) ± 0.16 (simulation)] g** ([CERN Courier](https://cerncourier.com/?p=109604)). The experiment used about 100 atoms per trap fill, released over 20 s ([sci.news](https://www.sci.news/physics/antimatter-gravity-12309.html)). Earlier bounds were ±7500% (2012) and ±100% (2013) ([Wikipedia](https://en.wikipedia.org/wiki/Gravitational_interaction_of_antimatter)). Combined in quadrature the uncertainty is 0.206 g, which puts the result **1.2σ from +1 g and 8.5σ from −1 g** *(computed)*.

Repulsive antigravity is therefore **[refuted] and resolved**. The magnitude at the percent level is **open**: the simulation systematic dominates, so the planned 1% laser-cooled run must improve modelling, not just counting statistics.

Ordinary matter shows no composition dependence at all:

- MICROSCOPE: η(Ti, Pt) = [−1.5 ± 2.3 (stat) ± 1.5 (syst)]×10⁻¹⁵, about ±2.75×10⁻¹⁵ in quadrature *(computed)* ([Touboul et al. 2022](https://arxiv.org/abs/2209.15487)).
- Eöt-Wash, Be–Ti: η = (0.3 ± 1.8)×10⁻¹³ toward Earth, and (−4 ± 7)×10⁻⁵ toward galactic dark matter ([Schlamminger et al. 2008](https://arxiv.org/abs/0712.0607)).
- 48 years of LLR: no violation at 10⁻¹⁴ ([Viswanathan et al. 2018](https://arxiv.org/abs/1710.09167)).

The weak equivalence principle is **[mainstream] and resolved at 10⁻¹⁵**. These bounds constrain bulk matter, though, not exotic condensates, which leaves the Tate result standing alone.

### The Tate anomaly is a 92 ppm result from 1989 that nobody repeated

Tate and colleagues measured the London moment of a rotating niobium ring and inferred a Cooper-pair mass of **m*/2m_e = 1.000084(21)** ([review, arXiv:2206.07574](https://arxiv.org/pdf/2206.07574)). Theory predicts about 8 ppm *below* 2m_e; that value comes from a secondary snippet ([Tajmar & de Matos](https://arxiv.org/abs/gr-qc/0607086)). The discrepancy is therefore **92 ppm, or 4.4σ** *(computed)*. The London moment itself is 2m_eω/e = 1.1371×10⁻¹¹ T per rad/s.

Tajmar and de Matos explain the excess with a gravitomagnetic term, B_g = 2ω(m* − m)/m = **1.84×10⁻⁴ ω**. That is about 10¹⁰ times Earth's gravitomagnetic field of roughly 10⁻¹⁴ rad/s ([arXiv:2206.07574](https://arxiv.org/pdf/2206.07574)). The evidence against that reading:

- An Air Force-funded search in high-temperature superconductors found nothing down to its resolution, though it did not test niobium ([DTIC ADA461570](https://apps.dtic.mil/sti/html/tr/ADA461570/index.html)).
- The rotating-superconductor anomalies the gravitomagnetic reading leaned on implied fields far above GR, and later work from the same group narrowed them ([Ruggiero 2023](https://arxiv.org/abs/2304.02167)).

The Tate result is **[contested] and open**, and the gravitomagnetic explanation is **[contested]** with no independent confirmation. The pattern matches the BIPM G result: one precision apparatus and no replication. The economical reading is an undiagnosed systematic, or a gap in the work-function and contact-potential corrections. A field 10¹⁰ times Earth's would be easy to see with a modern ring-laser gyroscope, so a repeat is a clean, cheap target for a parallel track.

### 2023–2026 brought one geophysical puzzle and no new physics

GRACE and satellite laser ranging found a transient north–south gravity anomaly at the Atlantic–Africa margin, peaking in January 2007 with a timescale of 2–3 years. Surface hydrology and ocean mass do not explain it. The proposed cause is vertical displacement of the perovskite to post-perovskite transition near the base of the African large low-shear-velocity province, and the anomaly coincides with the 2007 geomagnetic jerk ([EGU25-8808](https://meetingorganizer.copernicus.org/EGU25/EGU25-8808.html); [EGU26-12420](https://meetingorganizer.copernicus.org/EGU26/EGU26-12420.html)). It is **[contested] and open**, a mass-redistribution signal inside Earth rather than a gravity-physics test. That makes it a natural validation case for this project's satellite tooling.

A 2024 lattice atom interferometer found that a small source mass attracts exactly as Newton predicts. That excludes screened fifth forces over their natural parameter space ([Physics Today](https://physicstoday.aip.org/news/atom-interferometer-measures-centimeter-scale-gravity)), and a 10 m design proposes 1–1.5 more orders ([arXiv:2511.09750](https://arxiv.org/abs/2511.09750v1)).

The anomalous lunar eccentricity rate has shrunk as tidal models improved: (1.6 ± 0.5)×10⁻¹¹ yr⁻¹ in 2003, (9 ± 3)×10⁻¹² yr⁻¹ with DE421, and about (5 ± 2)×10⁻¹² yr⁻¹ with DE430 ([Iorio, arXiv:1102.0212](https://ar5iv.labs.arxiv.org/html/1102.0212); [Iorio, arXiv:1404.6537](https://arxiv.org/pdf/1404.6537)). Neither standard tides nor the modified-gravity models tested explain it. It remains **open**; the notes give it no label, and this report treats it as **[contested]**.

## Alternative theories survive only inside a 10⁻⁵ PPN box and a 10⁻¹⁵ speed bound

Every alternative theory has to pass three filters before its novel predictions matter: the Solar-System PPN parameters, the gravitational-wave speed, and the graviton mass. Most fail at least one.

### Solar-System numbers every metric theory must reproduce

The PPN formulas are ([Will 2014](https://arxiv.org/abs/1403.7377)):

- perihelion shift ω̇ = 42.98″ × [(2 + 2γ − β)/3 + 3×10⁻⁴(J₂/10⁻⁷)] per century;
- limb deflection ½(1 + γ) × 1.7505″;
- Nordtvedt parameter η_N = 4β − γ − 3.

The best current limits:

| Parameter | Limit | Effect | Source |
|---|---|---|---|
| γ − 1 | (2.1 ± 2.3)×10⁻⁵ | Shapiro delay, Cassini | [Will 2014](https://arxiv.org/abs/1403.7377) |
| β − 1 | (−1.6 ± 1.8)×10⁻⁵, with γ = 1 and J₂☉ = (2.246 ± 0.022)×10⁻⁷ | MESSENGER ranging | [Genova et al. 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC5773540) |
| β − 1, γ − 1 | (6.2 ± 7.2)×10⁻⁵, (1.7 ± 1.6)×10⁻⁴ | LLR | [Biskupek et al.](https://arxiv.org/abs/2012.12032) |
| ξ | 4×10⁻⁹ | spin precession, millisecond pulsars | [Will 2014](https://arxiv.org/abs/1403.7377) |
| α₁ | 7×10⁻⁵ | PSR J1738+0333 | [Will 2014](https://arxiv.org/abs/1403.7377) |
| α₂ | 2×10⁻⁹ | spin precession | [Will 2014](https://arxiv.org/abs/1403.7377) |
| α₃ | 4×10⁻²⁰ | pulsar Ṗ statistics | [Will 2014](https://arxiv.org/abs/1403.7377) |
| ζ₃ | 10⁻⁸ | lunar acceleration | [Will 2014](https://arxiv.org/abs/1403.7377) |

Hulse–Taylor orbital damping agrees with GR to better than half a percent.

The scripts reproduce the textbook numbers *(computed)*:

- **Mercury:** 42.98″ per century from a = 5.7909050×10¹⁰ m, e = 0.205630 and P = 87.9691 d.
- **Limb deflection:** 1.751″ with R☉ = 6.957×10⁸ m; Will's 1.7505″ uses a slightly larger radius.
- **Cassini Shapiro delay:** 262.4 µs round trip for the 2002 conjunction geometry.
- **Brans–Dicke:** with γ = (1 + ω)/(2 + ω), Cassini's 2σ edge gives **ω > 4.0×10⁴**, reproducing Will's 40,000. At that ω, γ − 1 = −2.4998×10⁻⁵ and η_N = +2.4998×10⁻⁵.
- **Linear massive gravity:** without Vainshtein screening, its effective γ = ½ gives 28.65″ per century, about 14″ short of Mercury's observed shift. The γ = ½ value is standard van Dam–Veltman–Zakharov (vDVZ) lore that the notes did not verify in a fetched source.

All scalar–tensor theories predict γ ≤ 1. Brans–Dicke is **[mainstream] as a framework and open, but squeezed** to ω > 4×10⁴. The Damour–Esposito-Farèse expressions for γ and β were written from background knowledge and need checking before coding. TeVeS predicts γ = β = 1, which makes Cassini irrelevant to it ([Will 2014](https://arxiv.org/abs/1403.7377); [Bekenstein 2004](https://arxiv.org/abs/astro-ph/0403694)). Einstein–Cartan torsion vanishes outside matter, so its PPN parameters are those of GR.

### The graviton is lighter than 2×10⁻²³ eV, but the number depends on the version

Massive gravity predicts dispersion, v_g/c ≈ 1 − ½(c/(λ_g f))², with λ_g = h/(m_g c) ([Will 2014](https://arxiv.org/abs/1403.7377)). Linear Fierz–Pauli theory suffers the vDVZ discontinuity, so viable versions need non-linear Vainshtein screening ([de Rham 2014](https://arxiv.org/abs/1401.4173)). The quoted gravitational-wave bound has moved between catalogue versions:

| Bound (90%) | Source | Status |
|---|---|---|
| m_g ≤ 1.27×10⁻²³ eV/c² | GWTC-3, v1 (Dec 2021) ([arXiv:2112.06861v1](https://arxiv.org/abs/2112.06861v1)) | superseded |
| m_g ≤ 2.42×10⁻²³ eV/c² | GWTC-3, revised v3, Phys. Rev. D 112, 084080 (2025) ([arXiv:2112.06861](https://arxiv.org/abs/2112.06861)) | current for GWTC-3 |
| 2.23×10⁻²³ eV/c² | GWTC-3 bound as restated in the GWTC-4.0 paper ([arXiv:2603.19020](https://arxiv.org/abs/2603.19020)) | discrepancy with v3 unexplained |
| **m_g ≤ 1.92×10⁻²³ eV/c²** | GWTC-4.0 (91 events) ([arXiv:2603.19020](https://arxiv.org/abs/2603.19020)) | newest |

The revision that loosened the GWTC-3 bound is not explained in the texts read.

Planetary ephemerides model the same physics as a Yukawa potential −(GM/r)e^(−r/λ_g):

| Bound | Source |
|---|---|
| m_g < 6.76×10⁻²³ eV/c² (90%) | INPOP17b ([Bernus et al. 2019](https://arxiv.org/abs/1901.04307)) |
| **m_g < 3.16×10⁻²³ eV/c²** (90%), λ_g > 3.93×10¹³ km | INPOP19a ([Bernus et al. 2020](https://arxiv.org/abs/2006.12304)) |
| λ_g > 1.2–2.2×10¹⁴ km, i.e. m_g < 1.03×10⁻²³ to 5.6×10⁻²⁴ eV/c² *(computed)* | Mars perihelion ([Will 2018](https://arxiv.org/abs/1805.10523)) |
| **m_g < 1.01×10⁻²⁴ eV/c²** (99.7%) | Mariani et al. 2023, cited by LVK as about an order of magnitude stronger than the GW bound ([arXiv:2603.19020](https://arxiv.org/abs/2603.19020)) |

Bernus warns that bounds derived from post-fit residuals are generically overestimated. The graviton mass is **[mainstream] as a bound and open**.

Reference values *(computed)*:

- λ_g = 6.46×10¹⁶ m (6.46×10¹³ km) at 1.92×10⁻²³ eV, and 3.92×10¹³ km at 3.16×10⁻²³ eV, matching Bernus's 3.93×10¹³ km.
- 1 − v_g/c = 1.08×10⁻²¹ at 100 Hz.
- The extra delay over 1 Gpc at 100 Hz is 1.1×10⁻⁴ s.

The notes state that the Yukawa factor at 1 au is 1 − 2.3×10⁻⁹. Recomputed here, r/λ_g = 1.496×10¹¹ m / 6.458×10¹⁶ m, so the factor is **1 − 2.3×10⁻⁶** *(checked here)*. The notes' value is off by 10³.

### GW170817 removed most relativistic dark-energy and MOND theories in one event

The 1.74 s delay between GW170817 and GRB 170817A fixes **−3×10⁻¹⁵ ≤ (v_GW − c)/c ≤ +7×10⁻¹⁶**. The bound assumes the conservative distance D = 26 Mpc and allows the gamma rays to have been emitted up to 10 s after merger ([LVK, Fermi & INTEGRAL 2017](https://arxiv.org/abs/1710.05834)). The script reproduces both edges: 1.74 s × c/D = +6.50×10⁻¹⁶ and (1.74 − 10) s × c/D = −3.09×10⁻¹⁵ *(computed)*. The same paper bounds the equivalence-principle difference in Shapiro delay to −2.6×10⁻⁷ ≤ γ_GW − γ_EM ≤ 1.2×10⁻⁶.

The speed bound eliminated a whole class of theories. It removes the covariant Galileon and most Horndeski and beyond-Horndeski dark-energy models, and extends to Einstein-Æther, Hořava, generalised Proca, TeVeS and other MOND-like gravities. Only narrow loopholes survive ([Ezquiaga & Zumalacárregui 2017](https://arxiv.org/abs/1710.05901); [Creminelli & Vernizzi 2017](https://arxiv.org/abs/1710.05877)).

"Dark-matter emulators" fail separately. In these theories matter and gravitational waves follow different metrics, which predicts a Shapiro delay of about 400 days against the 1.7 s observed, and gives |γ_GW − γ_EM| < 9.8×10⁻⁸ ([Boran et al. 2018](https://arxiv.org/abs/1710.06168)). Original TeVeS is **[refuted] and resolved**. AeST, a newer relativistic MOND built to fit the CMB and ghost-free at second order, is **[contested] and open** ([Skordis & Złośnik 2021](https://arxiv.org/abs/2007.00082)).

Verlinde's emergent gravity gives a MOND-like point-mass law with a₀ = cH₀/6 (1.09×10⁻¹⁰ m/s² for H₀ = 67.4 *(computed)*). Its record is mixed:

- Galaxy–galaxy lensing of 33,613 isolated galaxies agreed with it with no free parameters ([Brouwer et al. 2017](https://arxiv.org/abs/1612.03034)).
- KiDS-1000 found a difference of at least 6σ between early- and late-type galaxies at equal stellar mass, which no current modified-gravity theory explains ([Brouwer et al. 2021](https://arxiv.org/abs/2106.11677)).
- The RAR fits only with lowered mass-to-light ratios, and the theory predicts radius-correlated residuals that are not observed ([Lelli et al. 2017](https://arxiv.org/abs/1702.04355)).
- Its Solar-System perihelion advances miss by **seven orders of magnitude** ([Hees et al. 2017](https://arxiv.org/abs/1702.04358)).

The weak-field formula is **[refuted]**, and the wider programme is **[contested]**.

Einstein–Cartan torsion is **[mainstream] theory with no measurable deviation**. Torsion does not propagate, and spin–torsion effects rival mass effects only near **3.6×10⁴⁷ g/cm³ for electron-like matter and 1.2×10⁵⁴ g/cm³ for nucleon-like matter** *(computed; matching [Ruggiero & Tartaglia](https://arxiv.org/abs/gr-qc/0306029))*. That is about 39 orders of magnitude above neutron-star cores. Lorentz-violation searches constrain 19 of the 24 torsion components to about 10⁻³¹ GeV ([Kostelecký, Russell & Tasson 2008](https://arxiv.org/abs/0712.4393)).

### Mach's principle supplies coincidences, not a mechanism

Frame dragging in GR is **[mainstream]** and partly Machian. Sciama's 1953 inertial induction adds a term G m_A m_B **a**/(rc²) and relies on Ω = 4πρGT² being of order unity ([Wikipedia: Mach's principle](https://en.wikipedia.org/wiki/Mach%27s_principle)). For any flat Friedmann universe at critical density, that quantity is **exactly 3/2**, and GM_H/(R_H c²) is **exactly ½** *(computed)*. The "Φ/c² ≈ −1" of inertial induction is therefore built into the Friedmann equation, not independent evidence. It is **[contested]**: a heuristic.

Woodward's Mach-effect thrust received a covariant analysis from Rodal. In GR the relevant terms are suppressed by (U_N/c²)(ωL/c)², with no cosmological enhancement, and baryon conservation forbids an oscillating rest-mass monopole ([Rodal, arXiv:2512.22687](https://arxiv.org/abs/2512.22687)). For a typical 0.1 kg, 1 cm stack driven at 50 kHz, the suppression is **8×10⁻³⁷** *(computed)*. Woodward disputed Rodal's earlier claim in the same journal ([GRG reply](https://link.springer.com/article/10.1007/s10714-019-2654-7)). The GR analysis is **[mainstream]**, and the propulsion claim is **[refuted]**: there is no viable mechanism.

Haisch–Rueda–Puthoff zero-point-field inertia ([Rueda & Haisch 1998](https://arxiv.org/abs/physics/9802030)) faces two criticisms. Hadron inertia resides in chargeless gluons, so an electromagnetic zero-point field cannot be its source ([critique](https://dx.doi.org/10.1023/A%3A1017500513005)). Separately, a computational error was reported in Puthoff's 1989 zero-point gravity paper ([Wikipedia: Harold E. Puthoff](https://en.wikipedia.org/wiki/Harold_E._Puthoff)). It remains **[contested to refuted]** and unconfirmed.

### Heim, torsion fields and Weber gravity each fail a checkable test

**Heim theory** has one peer-reviewed summary (Z. Naturforsch. 1977). Its one checkable particle-physics prediction, two extra neutrinos, was excluded at LEP ([Wikipedia: Heim theory](https://en.wikipedia.org/wiki/Heim_theory)). The Dröscher–Häuser propulsion papers claimed about 25 T could cancel Earth's pull on a 150 t craft, with no derivation in a physics journal ([The Register](https://www.theregister.com/2006/01/06/hyperdrive/)). They also cited the rotating-superconductor signals that their own experimenter later attributed to artefacts ([Bentham chapter](https://benthambooks.com/book/9781608053995/chapter/100225/)). **[refuted]**, resolved.

**Russian "torsion fields"** were state-funded from 1986 to 1991. A 1989 decree created a centre for unconventional technologies, an academician's memo of May 1991 denounced closed pseudoscientific funding, and the State Committee for Science and Technology issued a resolution against the practice in July 1991 ([Kernbach, arXiv:1312.1148](https://arxiv.org/abs/1312.1148)). The Academy's Commission on Pseudoscience followed in 1998 ([Wikipedia](https://en.wikipedia.org/wiki/Commission_on_Pseudoscience)). The physical claims fail on their own terms:

- Signals at up to 10⁹ c that carry information but no energy violate special relativity.
- A claimed 80-fold cut in copper resistivity was not seen by independent tests.
- A 2008 "torsion engine" on the Yubileiny satellite did not change the orbit ([Wikipedia: Torsion field](https://en.wikipedia.org/wiki/Torsion_field_(pseudoscience))).

The closest quantitative cousin, a weight loss of up to 0.005% in right-spinning gyroscopes, failed three 1990 replications:

- a null result ([Faller et al.](https://jila.colorado.edu/bibcite/reference/4429));
- about 5% of the required effect, reduced more than tenfold by friction and temperature corrections ([Quinn & Picard](https://link.springer.com/article/10.1038/343732a0));
- a limit more than two orders below the claim ([Nitschke & Wilmarth](https://link.aps.org/doi/10.1103/PhysRevLett.64.2115)).

**[refuted]**, resolved.

**Weber-type gravity** has a clean no-go result. A force F = k(r̂/r²)(1 − μṙ² + γrr̈) is conservative only if γ = 2μ. The perihelion depends only on γ, while light deflection needs 2γ − μ = 3/c². So no conservative Weber-like force gets both Mercury and light bending right ([Bunchaft & Carneiro 1997](https://arxiv.org/abs/gr-qc/9708047)). A numerical RK4 orbit integration confirms Δω/Δω_GR = γc²/6, independent of μ *(computed)*:

| Force | μ, γ (×1/c²) | Advance / GR |
|---|---|---|
| Tisserand | 1, 2 | 0.3333 |
| Weber's original | ½, 1 | 0.1667 |
| Assis | 3, 6 | 0.9998 |

Weber electrodynamics also predicts a 0.65% inertial-mass shift for an electron inside a 10 kV shell *(computed)*, a large effect, and a repeat with optical coupling found none ([McMaster record](https://experts.mcmaster.ca/scholarly-works/1436968); [Weikert & Tajmar 2019](https://arxiv.org/abs/1902.05419)). **[refuted]**, resolved.

### Vacuum "free energy" fails on bookkeeping before it reaches physics

Casimir forces are **[mainstream]**. Forward's 1984 "vacuum-fluctuation battery" is a legitimate one-shot device: re-expanding the plates costs at least the energy gained ([PRB 30, 1700](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.30.1700)). Collapsing a 1 m² ideal plate pair from 1 µm to 10 nm yields at most **4.3×10⁻⁴ J, once** *(computed)*.

Moddel's group reported zero-bias current from metal–insulator–metal diodes beside submicron optical cavities: 1.4 pW in 0.02 µm² for a 33 nm cavity, steady over 4 hours ([Symmetry 13, 517, arXiv:2101.03085](https://arxiv.org/abs/2101.03085)). The budget rules out the cavity as the source *(computed)*:

- The power density is 70 W/m².
- Four hours at 1.4 pW is 2.0×10⁻⁸ J, about **8×10¹⁰ times** the ideal Casimir energy of the cavity (2.4×10⁻¹⁹ J).
- So the output cannot be released Casimir energy. It would have to be a continuous flow, which the second law forbids from a single reservoir in equilibrium.

The same group's 2019 review concluded that no proposal had been reliably demonstrated, and that mechanical Casimir extraction violates thermodynamics ([Moddel & Dmitriyeva](https://arxiv.org/abs/0910.5893)). No independent replication was found ([US11258379B2](https://patents.google.com/patent/US11258379B2/en)). The result is **[contested] and open**, but whatever it is, it is not free vacuum energy.

The general limits are **[mainstream]**:

- Quantum inequalities cap sampled negative energy density at −3ħ/(32π²c³t₀⁴), which is **−3.7×10⁻²⁶ J/m³ for a 1 ns sample** *(computed)* ([Ford & Roman 1995](https://arxiv.org/abs/gr-qc/9410043)).
- Ergotropy theory shows that a passive state, the ground state included, yields zero work under any cyclic unitary process ([Allahverdyan et al. 2004](https://arxiv.org/abs/cond-mat/0401574)). Dark energy (5.2×10⁻¹⁰ J/m³, so 0.52 J per km³ *(computed)*) is such a ground state, and is uniform as well.
- A closed cycle in a static gravitational field does zero net work. Falling from infinity to Earth's surface yields 62.6 MJ/kg, once. These are textbook results from Noether's theorem and were not separately sourced.

## Laboratory gravity works at piconewtons, and general relativity caps the rest

The gravity that people create and control in the laboratory is near-field Newtonian. It is measured to better than a percent, and it is the right benchmark for any claim of "gravity control".

### Newtonian calibrators are the gravity generators that work

LIGO Hanford's calibrator is an aluminium disk carrying tungsten cylinders: two at 6.033 cm, which force the mirror at twice the rotation frequency, and three at 10.478 cm, which force it at three times. It is mounted 1.18 m from a 39.68 kg test mass ([Ross et al. 2021](https://arxiv.org/abs/2107.00141)). The predicted 2f force is **19.18 ± 0.14 pN** from finite elements and 19.04 pN for point masses. The measurement gave 19.28 pN ± 1.12%, and the same data yield **G = (6.56 ± 0.41)×10⁻¹¹**, the first G measured with a gravitational-wave detector.

The closed-form point-mass formulas, reconstructed from garbled text extraction, give 19.07 pN and 8.978 pN; a brute-force sum over the cylinders gives 19.10 pN and 9.02 pN *(computed)*. The free-mass response x = F/(M(2πf)²) gives **5.0×10⁻¹⁷ m at 15.6 Hz, a strain of 1.26×10⁻²⁰** on an assumed 3995 m arm *(computed)*.

Virgo's calibrator is more precise still. A 2018 test predicted 3.3×10⁻¹⁴ (1 Hz/f)² m of mirror motion, which the formula reproduces as 3.331×10⁻¹⁴. The stated dimensions, however, give 1.89 kg per sector rather than the "about 2.06 kg" printed ([Estevez et al. 2018](https://arxiv.org/abs/1806.06572)). For the O4 run the calibrator became the absolute reference, with **0.17% total 1σ uncertainty** on the injected displacement, better than the photon calibrator, against 1.4% in O3 ([Aubin et al. 2024](https://arxiv.org/abs/2406.10028)).

At 15 Hz the gravitational wavelength is 2×10⁷ m, so these signals are near-field gravity, not waves. They are nonetheless the only gravity generators in routine use, and they produce about 10⁻²⁰ effective strain from a few kilograms of rotating metal. **[mainstream]**, working.

### The smallest source ever weighed is 92 mg; the Planck mass is next

A **92.1 mg gold sphere** is the smallest source whose gravity has been resolved ([Westphal et al. 2021](https://arxiv.org/abs/2009.09546)):

- The static force at 2.5 mm is 9×10⁻¹⁴ N (8.92×10⁻¹⁴ N *(computed)*), the same as a person 2.5 m away.
- The modulated acceleration was about 3×10⁻¹⁰ m/s²; the linear term computes to 2.66×10⁻¹⁰.
- G came out at **(6.04 ± 0.06)×10⁻¹¹**, 9% low but inside the known systematics.

The next milestone is a Planck-mass source (21.76 µg *(computed)*), which needs Q ≳ 20,000.

Other small-mass measurements:

- **Atom interferometry:** a 516 kg tungsten source produced a measured phase of 0.547870(63) rad and G = 6.67191(99)×10⁻¹¹ ([Rosi et al. 2014](https://arxiv.org/abs/1412.7954)). With k_eff·T² = 4.12×10⁵ rad per m/s², that phase is about 1.3×10⁻⁶ m/s² of differential acceleration *(computed; an inference, since the real extraction uses full trajectories)*.
- **Gravitational Aharonov–Bohm effect:** a 1.25 kg tungsten ring imprinted phases of −0.23 to −0.27 rad on ⁸⁷Rb wave packets 25 cm apart, and excluded a semiclassical model that predicted −0.198 to −0.394 rad ([Overstreet et al. 2022](https://www.science.org/doi/10.1126/science.abl7152); [arXiv:2209.02214](https://arxiv.org/abs/2209.02214)). The authors' reading that this is indirect evidence for quantum features of gravity is **[contested]**.
- **Levitated test mass:** a 0.43 mg magnet in a superconducting trap sensed 2.45 kg brass masses ([Fuchs et al. 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC10889343)). The paper quotes both 1030 aN and 30 aN for the signal and does not reconcile them.
- **Dynamic near-field gravity:** a beam vibrating near 42 Hz drove a detector beam through gravity alone, about 1 pm, matching Newtonian modelling to about 2% ([APS Physics](https://physics.aps.org/articles/v15/121)).

All are **[mainstream]**.

### The inverse-square law holds to 38.6 µm; below that, Casimir forces dominate

Short-range tests use the parameterisation V = V_N[1 + αe^(−r/λ)], whose point-mass force is F = (Gm₁m₂/r²)[1 + α(1 + r/λ)e^(−r/λ)]. Eöt-Wash measured from 52 µm to 3.0 mm with 18- and 120-fold attractors and found that **any |α| = 1 Yukawa must have λ < 38.6 µm (95%)** ([Lee et al. 2020](https://arxiv.org/abs/2002.11761)). The same data give:

- a dilaton or heavy-graviton mass above 5.1 meV;
- a largest extra dimension below 30 µm;
- a calibration ratio of 0.988 ± 0.005.

The dark-energy length (ħc/ρ_d)^(1/4) is 84.9 µm and the boson mass ħc/λ at 38.6 µm is 5.11 meV *(computed)*. The paper's printed "4√" is a fourth root, not a factor of 4. The tested range now reaches below half the dark-energy length.

The limits weaken fast at shorter range:

| Range | Limit on \|α\| | Source |
|---|---|---|
| 56 µm (superseded) | 1 | [Kapner et al. 2007](https://arxiv.org/abs/hep-ph/0611184) |
| 10 µm | 1.4×10⁴ | [Geraci et al. 2008](https://arxiv.org/abs/0802.2350) |
| 0.4–4 µm | best bound at the time; no α value in the abstract | Casimir torsion pendulum ([Sushkov et al. 2011](https://arxiv.org/abs/1108.2547)) |
| peak at 300 nm | up to 10³ better, but \|α\| < 1 only above 50 µm | Casimir-suppressed source; 34 pN of Casimir force had to be cancelled at 200 nm ([Chen et al. 2016](https://arxiv.org/abs/1410.7267)) |
| 20 pm–10 nm | 10× improvement | neutron Pendellösung ([Heacock et al. 2021](https://arxiv.org/abs/2103.05428)) |

Neutron-star cooling bounds scalar–nucleon couplings to g_N ≲ 5×10⁻¹⁴ at 10⁻⁶–10⁻¹² m ([Fiorillo et al. 2025](https://arxiv.org/abs/2506.19906)). That maps to α ≈ 3.4×10¹⁰ *(computed; the 4π and nucleon-mass convention needs checking)*.

The law is **[mainstream] and open below 38.6 µm**. The obstacle is Casimir and patch-potential backgrounds, not gravity.

### Forward's "antigravity" is correct GR with a 7×10⁻²⁸ coupling

Forward's "General Relativity for the Experimentalist" (1961) and "Guidelines to Antigravity" (1963) are legitimate linearised GR **[mainstream]** ([Ruggiero 2023](https://arxiv.org/abs/2304.02167); [AJP 31, 166](https://pubs.aip.org/aapt/ajp/article/31/3/166/1046748/Guidelines-to-Antigravity)). His toroidal dipole is G_D = −(µ_G/4π)d/dt[u(r/R)²], with µ_G = 4πG/c². Its coefficient is "of the order of 10⁻²⁸" ([arXiv:1306.4893](https://arxiv.org/abs/1306.4893)), which is exactly **G/c² = 7.43×10⁻²⁸ m/kg** *(computed)*.

A linearised-GR scaling written for this purpose, a ≈ 4Gma_f/(c²r), shows the size of the problem *(computed; an inference)*. Ten tonnes accelerating at 10⁴ m/s² induce 3×10⁻¹⁹ m/s² at 1 m. Inducing 1 g with a_f = 10⁶ m/s² needs m/r ≈ 3.3×10²¹ kg/m, a compactness between Earth (9.4×10¹⁷ kg/m) and a white dwarf.

Frame dragging gives the same verdict *(computed)*:

- A 10 t tungsten sphere spun to its burst speed drags frames at **5.5×10⁻²¹ rad/s** at its pole, 1.2×10⁻²³ of its own spin rate.
- A thin shell of 10⁴ kg spun at 300 rad/s gives 3.0×10⁻²¹ rad/s inside.
- Gravity Probe B's measured −37.2 ± 7.2 mas/yr, against GR's −39.2 mas/yr ([Everitt et al. 2011](https://arxiv.org/abs/1105.3456)), is 6.0×10⁻¹⁵ rad/s. That is a million times the laboratory value.
- LARES-2 and LAGEOS have since measured Earth's frame dragging to about 10⁻³ ([Nature 2026](https://www.nature.com/articles/s41586-026-10715-0)).

None of the 1977 Braginsky–Caves–Thorne laboratory proposals has been performed ([Caltech record](https://authors.library.caltech.edu/records/re5v4-65528); [Ruggiero 2023](https://arxiv.org/abs/2304.02167)). These are **[mainstream] proposals and open**.

### Laboratory gravitational waves sit ten or more orders from detection

No laboratory-generated gravitational wave has ever been detected. The quadrupole formula for a rigid rotor, P = (32/5)(G/c⁵)(I₁ − I₂)²Ω⁶ ([lecture exercise](https://indico.fysik.su.se/event/8276/contributions/13166/attachments/5761/7341/stockholm24ex.pdf)), sets the scale *(computed)*:

| Source | Power | Strain |
|---|---|---|
| MTW 490 t steel beam spun at its tensile limit ([Chiao et al.](https://arxiv.org/pdf/1712.08680)) | **2.2×10⁻²⁹ W** | 1.2×10⁻⁴¹ at one wavelength (3.4×10⁷ m) |
| 1 kg, 0.2 m bar at 100 Hz | 1.2×10⁻⁴⁰ W | |

Weber's own laboratory analysis gave about 10⁻²⁰ W ([Kolosnitsyn & Rudenko](https://arxiv.org/abs/1504.06548)). The best concepts on paper are electromagnetic:

- **Gertsenshtein "Hertz experiment"** (10 T, 100 m, a 1 − R = 10⁻⁸ cavity and a 10 kW pump): h ~ 10⁻³². That would deliver about 10⁻³ photons/s to a receiver and need 4×10⁵ s of integration. A finesse of 10¹⁰ reaches 10⁻³⁰ ([Kolosnitsyn & Rudenko](https://arxiv.org/abs/1504.06548)).
- **Petawatt twisted light:** about 10⁻³⁶ at µm distances ([Atonga et al.](https://arxiv.org/html/2309.04191v2)).

These schemes are physically sound but unrealised. Material strength caps the rim speed of any rotor at roughly √(σ/ρ), so mechanical sources cannot close the gap; this is a scaling inference. The LIGO calibrator's 1.3×10⁻²⁰ effective strain is about 14 orders above the formal radiative strain of the MTW beam at 1 m, because it is a near field. Chiao's superconducting "GW laser" is **[contested]**, and his only experiment saw nothing ([Chiao, Fitelson & Speliotopoulos](https://arxiv.org/abs/gr-qc/0304026)).

### Tabletop quantum gravity has one proposal, one refutation and one detector design

The BMV entanglement proposal would place 10⁻¹⁴ kg diamonds in 250 µm superpositions, 450 µm apart, for about 2.5 s at 10⁻¹⁵ Pa and 0.15 K ([Bose et al. 2017](https://arxiv.org/abs/1707.06050)). The phases φ ~ Gm₁m₂τ/(ħd) compute to 0.352, 0.791 and 0.226 rad, so Δφ_RL = 0.440 and Δφ_LR = −0.126 *(computed)*. The paper quotes about 0.7 and −0.2, roughly 1.6 times larger. The mutual acceleration is 1.7×10⁻¹⁷ m/s².

The paper's witness also needs care. The printed operator |⟨σx⊗σz⟩ − ⟨σy⊗σz⟩| can reach √2 on product states, so W > 1 certifies nothing. The σy⊗σy form is bounded by 1 on separable states, and its optimum over local bases is **1 + |sin((Δφ_LR + Δφ_RL)/2)|**. That gives 1.156 with the computed phases, matching the published "W ∼ 1.16" *(computed)*. Check the published version before coding.

No gravitationally induced entanglement had been observed as of October 2026. The enabling steps are 2025 levitated-particle delocalisation and squeezing ([reference list via arXiv:2602.19306](https://arxiv.org/pdf/2602.19306)), together with diamagnetic chip traps ([arXiv:2411.02325](https://arxiv.org/abs/2411.02325)). The proposal is **[mainstream] and open**.

The parameter-free Diósi–Penrose collapse model is **[refuted]**. An underground germanium detector set **R₀ > 0.54×10⁻¹⁰ m** against Penrose's natural 0.05×10⁻¹⁰ m ([Donadi et al. 2021](https://arxiv.org/abs/2111.13490)). At ideal efficiency, the radiation formula predicts about 11 photons per day at the bound and 1.3×10⁴ at the natural value, a scaling of R₀⁻³ *(computed)*. Generalised versions remain open.

Tobar et al. showed that a ground-state-cooled bar of about 15 kg could absorb single gravitons from a GW170817-like signal ([arXiv:2308.15440](https://arxiv.org/abs/2308.15440)). The design is **[mainstream] as a proposal**, though whether a click proves quantisation is **[contested]** (see the graviton section below). Dyson's conjecture that single gravitons are undetectable in LIGO-type interferometers is a **[mainstream] conjecture and open** ([IAS essay](https://www.ias.edu/ideas/2013/dyson-graviton)).

### Superconducting transducers rest on heuristics and one null experiment

Chiao claimed that superconductors convert electromagnetic waves into gravitational waves with efficiency of order unity, or reflect gravitational waves outright ([gr-qc/0204012](https://arxiv.org/abs/gr-qc/0204012); [Minter, Wegter-McNelly & Chiao](https://arxiv.org/abs/0903.0661)). The arguments rest on three ingredients:

- impedance matching to Z_g = 4πG/c = **2.80×10⁻¹⁸ m²/(kg s)**;
- a (c/v_s)⁴ ≈ 10²⁰ enhancement;
- a "Heisenberg–Coulomb" effect. The electron Coulomb-to-gravity ratio is e²/(4πε₀Gm_e²) = 4.17×10⁴², matching the "42 orders" in press coverage ([MIT Technology Review](https://www.technologyreview.com/2009/03/23/214718/if-superconducting-sheets-reflected-gravitational-waves/amp/)).

The only experiment, using YBCO at 12 GHz and 77 K, set **η < 1.6×10⁻⁵** ([gr-qc/0304026](https://arxiv.org/abs/gr-qc/0304026)). Kiefer and Weber judged the arguments heuristic and not based on detailed calculations, and Licht's classical Ginzburg–Landau calculation found no strong enhancement ([Kiefer & Weber 2005](https://arxiv.org/abs/gr-qc/0408010)).

Weinberg's classical efficiency for a 1000 kg bar is **5×10⁻³³** *(computed)*. Even multiplied by 10²⁰ it is 10⁻¹³, so the YBCO limit does not constrain the conservative estimate. Later supporting work is theory only ([Bahamonde et al. 2020](https://arxiv.org/abs/2005.06640)). The transducer claims are **[contested] and open**, with no positive result in 23 years.

## Gravitational waves carry a thousandth of c⁵/G yet pass through matter untouched

The deep physics of gravitational waves is settled at leading order and verified at 10⁻⁴. The open items are effects GR predicts but nobody has yet seen.

### Energy, luminosity and a misattributed bound

The Isaacson effective stress tensor gives a plane-wave energy flux of **F = c³h₀²ω²/(32πG)**. That is 1.5 mW/m² for h₀ = 10⁻²² at 1 kHz, comparable to full moonlight; the formula gives 1.585 *(computed)* ([Flanagan & Hughes 2005](https://arxiv.org/abs/gr-qc/0501041)). The GW150914 properties paper used c³|ḣ|²/(16πG) and found about 10⁵ erg s⁻¹ m⁻² at 250 Hz ([PRL 116, 241102](https://arxiv.org/abs/1602.03840)). That is 0.0198 W/m² *(computed)*, twice the cycle-averaged form, against 0.0174 W/m² from the peak luminosity spread over 4π(410 Mpc)².

The multipole flux coefficients are exact rationals *(computed from [Blanchet's Living Review](https://arxiv.org/abs/1310.1528))*:

| Flux | Coefficients |
|---|---|
| Energy | a₂ = 1/5 (the quadrupole formula), b₂ = 16/45, a₃ = 1/189, b₃ = 1/84, a₄ = 1/9072, b₄ = 4/14175 |
| Angular momentum | a′₂ = 2/5 |
| Recoil | a″₂ = 2/63, c″₂ = 16/45 |

The Peters–Mathews decay formula is verified by two binary pulsars:

- **Hulse–Taylor:** intrinsic Ṗ_b/GR = **0.9983 ± 0.0016**, with a GR value of −2.40263(5)×10⁻¹² ([Weisberg & Huang 2016](https://arxiv.org/abs/1606.02744)). The script gives f(e) = 11.857 and −2.4021×10⁻¹², limited by the masses being quoted to 10⁻³.
- **Double pulsar:** **0.999963(63)**, confirming wave emission to 1.3×10⁻⁴ ([Kramer et al. 2021](https://arxiv.org/abs/2112.06795)). The script gives −1.24781×10⁻¹².

These formulas are **[mainstream] and resolved**.

The events themselves:

- GW150914 radiated **3.0 M☉c² (5.4×10⁴⁷ J)** and peaked at **3.5×10⁵⁶ erg/s (200 M☉c²/s)** ([PRL 116, 241102](https://arxiv.org/abs/1602.03840)). That is 9.65×10⁻⁴ c⁵/G *(computed)*.
- GW190521, at 150 M☉, peaked at 3.7×10⁵⁶ erg/s, or 1.02×10⁻³ c⁵/G ([ApJL 900, L13](https://arxiv.org/abs/2009.01190)).
- GW250114 reached network SNR 80 and confirmed Hawking's area law ([arXiv:2509.08054](https://arxiv.org/abs/2509.08054)). That paper quotes SNR 26 for GW150914 against 24 in the discovery paper, a convention difference that tests must not trip over.

Peak luminosity is therefore scale-free: about 10⁻³ c⁵/G whatever the mass. The leading-order flux at the innermost stable circular orbit is 5.1×10⁻⁵ c⁵/G, about 20 times below the merger peak *(computed)*. Any tool that extrapolates the quadrupole formula through merger is badly wrong.

The Planck luminosity c⁵/G = **3.63×10⁵² W** is a unit, not a theorem:

- Suitable initial data can produce arbitrarily large luminosity; requiring no past horizon appears to keep it below c⁵/G; critical collapse reaches about 0.2 c⁵/G ([Cardoso et al. 2018](https://arxiv.org/abs/1803.03271)).
- Any bound must depend on the source's internal structure ([Jowsey & Visser 2021](https://arxiv.org/abs/2105.06650)).

The unit is **[mainstream]**; the bound is **[contested] and open**. The same paper reproduces a note in which Dyson denies ever proposing c⁵/G as a limit, so the attribution is **[refuted]**. It traces to his 1962 "Gravitational Machines" essay. Its Eq. 11, W = 128V¹⁰/(5Gc⁵), reproduces his own 2×10³⁷ erg/s ([Dyson 1963 scan](https://www.gwern.net/doc/science/1963-dyson.pdf)). But the Peters formula for the same binary gives **3.3×10³⁸ erg/s, 16 times more** *(computed)*. This discrepancy was found in the notes, is **[contested]** until checked against the printed book, and makes a useful documented counter-example.

Gravitational-wave recoil is real: GW200129's remnant was kicked at about **1542 km/s**, or 0.5% of c ([Varma et al. 2022](https://arxiv.org/abs/2201.01302)). That event's evidence is **[contested]**, because a data-quality reanalysis exists that the notes did not read.

### Memory is predicted, hereditary, and still undetected

The nonlinear (Christodoulou) memory is a permanent offset in the + polarisation. It is sourced by the waves' own energy flux and depends on the source's whole past. Its leading post-Newtonian form is **h₊^(mem) = (2ηMx/R)·(1/96)sin²Θ(17 + cos²Θ)**, which vanishes face-on ([Favata 2010](https://arxiv.org/abs/1003.3486)). Cut at the innermost stable orbit and edge-on, the formula gives **1.1×10⁻²² for GW150914** and 1.4×10⁻¹⁵ for a 10⁹ + 10⁹ M☉ binary at 1 Gpc *(computed)*. Both are lower bounds, because the merger adds more.

No memory of any kind had been detected as of October 2026, so it is **[mainstream] theory and open**:

- **Forecast:** at design sensitivity, LIGO needs about 35 (or 90) GW150914-like events for SNR 3 (or 5) ([Lasky et al. 2016](https://arxiv.org/abs/1605.01415)).
- **GWTC-3:** stacking 88 events gave ln BF = 0.01 ([Cheung, Lasky & Thrane 2024](https://arxiv.org/abs/2404.11919)).
- **GWTC-5.0:** a hierarchical analysis measures the memory coefficient at **0.26 (+4.09/−4.08)** against GR's 1, and needs about 2000 events for 1σ ([Mitman, Isi & Farr 2026](https://arxiv.org/abs/2605.27500)).
- **(2,0) mode:** a July 2026 preprint finds log₁₀B = 1.38 ± 0.79 for the full (ℓ = 2, m = 0) mode in GWTC-4.0, and estimates 166 events for decisive evidence ([Rosselló-Sastre et al.](https://arxiv.org/abs/2607.04909)). It is **[contested]**, and it is not a memory detection: it tests a harmonic that includes an oscillatory part.
- **Pulsar timing:** NANOGrav's best memory Bayes factor is 3.1, and no sky region reaches below h₀ = 10⁻¹⁴ ([Agazie et al. 2025](https://arxiv.org/abs/2502.18599)). The timing-residual ramp for h₀ = 10⁻¹⁴ and B = ½ is 0.79 µs after 5 years *(computed)*.
- **Spin memory** is a prediction with no measurement, not even a limit ([Nichols 2017](https://arxiv.org/abs/1702.03300)).

### Propagation: c to 10⁻¹⁵, and absorption negligible by twenty orders

The speed and graviton-mass bounds are covered above. On lensing, wave optics replaces geometrical optics when the wavelength exceeds the lens's Schwarzschild radius. The point-lens amplification uses w = 8πGM_Lz f/c³ ([Takahashi & Nakamura 2003](https://arxiv.org/abs/astro-ph/0305055)). Reference values *(computed)*:

- w = 1.24 for 100 M☉ at 100 Hz, and w = 1 at 80.8 M☉.
- |F(y = 0)|² = πw/(1 − e^(−πw)) = 3.2835 at w = 1.
- At y = 1, the magnifications are μ± = 1.17082 and −0.17082, and the time delay is 4099 s for 10⁸ M☉.

The LVK O4a searches found no strong lensing. One point-mass outlier, **GW231123**, has waveform uncertainties that complicate its interpretation, so it is **[contested]** ([arXiv:2512.16347](https://arxiv.org/abs/2512.16347)).

Matter barely absorbs gravitational waves. A viscous fluid damps them at γ = 16πGη/c² ([Baym, Patil & Pethick 2017](https://arxiv.org/abs/1707.05192)). That is 3.7×10⁻²⁹ s⁻¹ in water, an absorption length of **2.6×10¹¹ Gpc** *(computed)*. The formula holds only for ωτ ≪ 1. Feeding it a mantle viscosity of 10²¹ Pa s gives a spurious 54 au, which is why code must refuse that regime. Other propagation effects:

- Free-streaming neutrinos cut primordial tensor power by 35.6% ([Weinberg 2004](https://arxiv.org/abs/astro-ph/0306304)).
- Dark-matter dispersion is negligible ([Flauger & Weinberg 2018](https://arxiv.org/abs/1801.00386)).

Gertsenshtein conversion in transit is tiny. The rate is P ≈ 8.3×10⁻⁵⁰ (B/1 G)²(L/1 cm)²; CMB distortions then limit backgrounds near 10¹⁰ Hz to Ω_GW ≲ 1 ([Pshirkov & Baskaran 2009](https://arxiv.org/abs/0903.4160)). In the LIGO band, conversion to a propagating photon is forbidden because 100 Hz lies below the interstellar plasma frequency (1.56 kHz at 0.03 cm⁻³ *(computed)*).

One conflict the notes did not catch. Rothman and Boughn's mixing length, L = 2/(G^(1/2)B) (2.25 Mpc at 1 G), and its probability P = sin²(G^(1/2)BD/2) ([Rothman & Boughn 2006](https://arxiv.org/abs/gr-qc/0601043)) have a small-angle limit **exactly four times smaller** than Pshirkov's prefactor and the SI form 4πGB²L²/(μ₀c⁴) used in `gwsens` *(checked here)*. For 5 µG over 8 kpc, the Rothman–Boughn form gives 3.1×10⁻¹⁶ and the SI form 1.26×10⁻¹⁵. The SI form implies a mixing length of 1.13 Mpc at 1 G.

### Polarisations: pure scalar and pure vector are dead, admixtures unconstrained

A generic metric theory allows six polarisations, while GR allows two tensor modes. For L-shaped interferometers the breathing and longitudinal scalar modes respond identically, so at most five modes can be separated. Separating them needs more detectors than modes, and only two were operating in O1 and O4a ([GWTC-3 TGR](https://arxiv.org/abs/2112.06861); [GWTC-4.0 TGR I](https://arxiv.org/abs/2603.19019)). With GW170817's sky position fixed, tensor beats vector by log₁₀B = +20.81 ± 0.08 and scalar by +23.09 ± 0.08 ([LVC 2019](https://arxiv.org/abs/1811.00364)). NANOGrav's correlations follow the tensor Hellings–Downs curve, with Bayes factors of 200–1000 ([Agazie et al. 2023](https://arxiv.org/abs/2306.16213)).

The GWTC-4.0 Bayes factors against pure tensor are ([Table 9](https://arxiv.org/abs/2603.19019)):

| Hypothesis | GWTC-3.0 | GWTC-4.0 |
|---|---|---|
| Scalar | −5.76 ± 0.42 | **−14.72 ± 0.59** |
| Vector | −2.39 ± 0.41 | **−5.33 ± 0.58** |
| Tensor + scalar | −0.07 ± 0.41 | −0.20 ± 0.57 |
| Tensor + vector | −0.04 ± 0.41 | +0.10 ± 0.57 |
| Tensor + vector + scalar | −0.19 ± 0.41 | −0.31 ± 0.57 |

Pure scalar and pure vector are **excluded (resolved)**. Mixtures sit at zero within ±0.6, so a small scalar or vector admixture is **[mainstream] question, open**. Brans–Dicke's breathing mode, at ω > 4×10⁴, lives exactly there. Settling it needs LIGO-India, or Virgo and KAGRA operating fully.

### Gravitons: 10³⁶ per wave, about one quantum per bar

Every detected wave is classical in the counting sense:

- N = h₀²c⁵/(32πGħω²) gives **3.9×10³⁶ gravitons** for h₀ = 10⁻²¹ at 150 Hz *(computed)* ([Tobar et al. 2024](https://arxiv.org/abs/2308.15440)).
- The occupation per de Broglie cell is 2×10³⁵ for h = 10⁻²² at 1 kHz ([Carney, Domcke & Rodd 2024](https://arxiv.org/abs/2308.12988)).
- By the same formula, a transient at CAST's sensitivity near 10¹⁸ Hz would have 10⁻⁵ quanta per cell, a single-graviton regime.

Graviton–matter cross-sections sit near ℓ_P² = **2.6×10⁻⁶⁶ cm², 3.9×10⁻⁴² of Thomson's** *(computed)*. A detector is capped at about a Jupiter mass (2.3×10²⁷ kg) before it collapses ([Rothman & Boughn 2006](https://arxiv.org/abs/gr-qc/0601043)). Thermal gravitons from the Sun would register about four times over the Sun's lifetime in an Earth-sized detector ([Quanta 2024](https://www.quantamagazine.org/it-might-be-possible-to-detect-gravitons-after-all-20241030/)).

Tobar's rates *(computed)*:

- **Stimulated:** Γ_stim = v_s²Mh²/(4l²π³ħ) is **1.00 Hz** for an 1800 kg aluminium bar at h = 5×10⁻²².
- **Spontaneous:** Γ_spon = 8πGρv_s⁴R²/(Lc⁵) is 9.3×10⁻³⁴ Hz for niobium.
- **Optimal mass:** 12.8 kg of beryllium at 100 Hz and 220 kg of aluminium at 150 Hz, against the paper's ~15 kg and ~250 kg.

An independent classical check agrees: a 2300 kg bar absorbs **1.15×10⁻³¹ J, or 0.19ħω**, from GW170817 *(computed)*. The formulas are **[mainstream]**. Carney et al. show that a classical wave acting on a quantised detector produces the same clicks, so proving quantisation would need sub-Poissonian counting statistics. The interpretation is **[contested]**.

### Primordial waves: r < 0.034 caps inflation at 1.4×10¹⁶ GeV

The bound on the tensor-to-scalar ratio has tightened from 0.036 (BICEP/Keck 2021 ([PRL 127, 151301](https://arxiv.org/abs/2110.00483))) to 0.032 (Planck PR4 + BK18 ([Tristram et al.](https://arxiv.org/abs/2112.07961))) and **r < 0.034** in the end-of-2025 combination of all CMB data. That combination also gives n_s = 0.9682 ± 0.0032 ([Balkenhol et al.](https://arxiv.org/abs/2512.10613)).

With V = (3π²/2)A_s r M_pl⁴ and H = πM_pl√(A_s r/2) ([Baumann TASI](https://arxiv.org/abs/0907.5424)), and ln(10¹⁰A_s) = 3.044 ([Planck 2018](https://arxiv.org/abs/1807.06209)), the scale of inflation is capped *(computed)*:

| r | V^(1/4) (GeV) | H_inf (GeV) |
|---|---|---|
| 0.036 | **1.408×10¹⁶** | 4.70×10¹³ |
| 0.032 | 1.368×10¹⁶ | 4.43×10¹³ |
| 0.01 | 1.022×10¹⁶ | 2.48×10¹³ |

The present-day background is Ω_GW h² ≤ **1.3×10⁻¹⁶** *(computed from [Caprini & Figueroa](https://arxiv.org/abs/1801.04268))*, neglecting the drop in g* and the neutrino damping, both of which lower it. Primordial waves are **[mainstream] bounds and open**: not yet detected.

The nanohertz background is about seven orders above that, at **A = 2.4 (+0.7/−0.6)×10⁻¹⁵**. Its Hellings–Downs evidence is about 3–4σ ([NANOGrav 15-yr](https://arxiv.org/abs/2306.16213)). How to read it:

- **Supermassive black-hole binaries:** fit only with parameters at the edges of expectation, **[contested]** ([NANOGrav](https://arxiv.org/abs/2306.16220)).
- **New physics:** fits with Bayes factors of 10–100, which the collaboration says should not be regarded as evidence, **[contested]** ([NANOGrav](https://arxiv.org/abs/2306.16219)).
- **MeerKAT:** 3–3.4σ under one noise model and insignificant under another ([Miles et al. 2025](https://arxiv.org/abs/2412.01153)).

### Standard sirens cannot yet arbitrate the Hubble tension

GW170817 alone gives **H₀ = 70.0 (+12.0/−8.0) km s⁻¹ Mpc⁻¹**; the naive v/d = 3017/43.8 is 68.9 ([Nature 551, 85](https://arxiv.org/abs/1710.05835)). GWTC-4.0 dark sirens combined with GW170817 give **75.4 (+12.8/−9.1)**, with Ξ₀ = 1.4 (+1.0/−0.4) against GR's 1 ([ApJL 1007, L17](https://arxiv.org/abs/2509.04348)). Independent reanalyses that use a mass feature, or cross-correlate with DES, sit between 72 and 83. They were read from summaries only and are **[contested]** ([arXiv:2601.03257](https://arxiv.org/abs/2601.03257v1)). The 68% intervals cover both 67 and 73, so siren cosmology is **[mainstream] and open**.

### What a wave does to matter is set by G/c³, not by its energy flux

The tidal acceleration ½ω²hL is **1.2×10⁻¹⁵ m/s² across 1 m** for h = 10⁻²¹ at 250 Hz, and a 4 km arm changes by 2×10⁻¹⁸ m *(computed from [Flanagan & Hughes](https://arxiv.org/abs/gr-qc/0501041))*. A GW150914 analogue at 8 kpc would deliver **4.6×10⁷ W/m², 3.4×10⁴ times the solar constant**, yet stretch a metre at only 6.3×10⁻¹¹ m/s² *(computed)*. The bar that absorbs 0.19 quanta from GW170817 gains Δv = 1.7×10⁻⁴³ m/s, because its integrated cross-section is 4×10⁻²⁵ m² Hz ([Aguiar 2011](https://arxiv.org/abs/1009.1138)). Recomputed, that is 4.32×10⁻²⁵ m² Hz, though the paper's quoted area ratio implies a different physical area.

The luminosity scale maps onto a strain ceiling, **h₀ = 2√2·√℘·c/(ωr)** for L = ℘c⁵/G *(computed)*. With ℘ = 10⁻³ at 250 Hz and 410 Mpc it gives 1.35×10⁻²¹, reproducing GW150914's 1.0×10⁻²¹. Strain reaches order unity only within about a wavelength of a Planck-luminosity source.

Two further results:

- A collimated wave beam would give thrust P/c, 3.3×10⁻⁹ N/W, the same as a photon rocket.
- Ten years of superconducting-gravimeter data, read through Earth's normal modes at 0.3–5 mHz, bound the background to **Ω_GW between 0.035 and 0.15** ([Coughlin & Harms 2014](https://arxiv.org/abs/1406.1147)). That is the one published gravitational-wave bound produced by gravimeters, the instruments this project already processes.

## The ledger: every item, its label and status, and where a parallel track pays

| Item | Label | Status | Deciding number |
|---|---|---|---|
| Spread of G measurements | [mainstream] | open | 551 ppm; Birge ratio 3.61 before ×3.9 |
| NIST 2026 rebuild of BIPM balance | [mainstream] | open (cause unidentified) | −250 to −282 ppm vs BIPM-14 |
| G tracks length of day | [refuted] | resolved | Bayes factor ≥ e³⁰ for noise |
| Earth flyby residuals | [mainstream] | open, spacecraft-specific | NEAR 13.46 mm/s |
| Anderson flyby formula | [contested] | failed out of sample | Juno: +6 predicted, 0 ± 0.8 observed |
| Majorana gravitational shielding | [refuted] | resolved | h < 10⁻²² vs 6.7×10⁻¹³ m²/kg |
| Allais pendulum effect | [contested] | open, unreplicated | 2010 null < 0.3°/h |
| 1997 eclipse gravimeter dips | [contested] | open, unreplicated | needs 10⁷ × the LLR bound |
| Pioneer as new physics | [refuted] | resolved | thermal τ ≈ 29–37 yr |
| Pioneer thermal recoil | [mainstream] | resolved | 26.3% thermal error budget |
| Radial acceleration relation | [mainstream] | resolved as data | g† = 1.20×10⁻¹⁰ m/s² |
| MOND as its explanation | [contested] | open | quadrupole tension 8.7σ |
| Wide-binary gravity boost | [contested] | open, decidable | 1.4–1.6 vs Newton at 19σ |
| Repulsive antimatter gravity | [refuted] | resolved | 8.5σ from −1 g |
| Antimatter gravity magnitude | [mainstream] | open at % level | 0.75 ± 0.21 g |
| Weak equivalence principle | [mainstream] | resolved at 10⁻¹⁵ | η = (−1.5 ± 2.75)×10⁻¹⁵ |
| Tate Cooper-pair mass | [contested] | open | 92 ppm, 4.4σ |
| Tajmar–de Matos gravitomagnetism | [contested] | unconfirmed | B_g = 1.84×10⁻⁴ ω |
| GRACE 2007 Atlantic transient | [contested] | open (geophysics) | peak January 2007 |
| LLR eccentricity rate | [contested]* | open, shrinking | (5 ± 2)×10⁻¹² yr⁻¹ |
| Screened fifth forces | [mainstream] null | natural space excluded | lattice interferometer 2024 |
| PPN γ, β | [mainstream] | GR holds | ~2×10⁻⁵ |
| Brans–Dicke / scalar–tensor | [mainstream], constrained | open | ω > 4×10⁴ |
| Graviton mass | [mainstream] bound | open | < 1.92×10⁻²³ eV (GW); < 1.01×10⁻²⁴ (ephemerides) |
| Einstein–Cartan torsion | [mainstream] theory | untestable | 3.6×10⁴⁷ g/cm³ |
| Original TeVeS | [refuted] | resolved | 400-day Shapiro delay vs 1.7 s |
| AeST | [contested] | open | survives c_g = c |
| Verlinde weak-field formula | [refuted] | resolved | perihelia off by 10⁷ |
| Sciama inertial induction | [contested] | heuristic | 4πGρ_cT² = 3/2 exactly |
| Woodward Mach-effect thrust | [refuted] (GR analysis [mainstream]) | resolved | suppression 8×10⁻³⁷ |
| Zero-point-field inertia | [contested to refuted] | unconfirmed | gluon mass |
| Heim theory | [refuted] | resolved | extra neutrinos excluded at LEP |
| Torsion fields | [refuted] | resolved | 1990 gyroscope nulls |
| Weber-type gravity | [refuted] | resolved | perihelion or light bending, not both |
| Forward vacuum battery | [mainstream] | one-shot only | 4.3×10⁻⁴ J per m² |
| Moddel cavity current | [contested] | unreplicated | 8×10¹⁰ × the cavity's Casimir energy |
| Quantum inequalities | [mainstream] | resolved | −3.7×10⁻²⁶ J/m³ at 1 ns |
| Speed of gravity | [mainstream] | resolved | −3×10⁻¹⁵ to +7×10⁻¹⁶ |
| Newtonian calibrators | [mainstream] | working | 0.17% (Virgo O4) |
| 92 mg source mass | [mainstream] | resolved | G 9% low, within systematics |
| Atom-interferometer G | [mainstream] | resolved | 150 ppm |
| Gravitational Aharonov–Bohm | [mainstream]; quantum reading [contested] | resolved / open | −0.23 to −0.27 rad |
| Levitated mg gravimetry | [mainstream] | force value unreconciled | 1030 vs 30 aN in paper |
| Inverse-square law | [mainstream] | open below 38.6 µm | \|α\| < 1.4×10⁴ at 10 µm |
| Forward GR designs | [mainstream] physics | impractical | G/c² = 7.4×10⁻²⁸ m/kg |
| Laboratory frame dragging | [mainstream] proposal | never performed | 5.5×10⁻²¹ rad/s |
| Laboratory GW generation | [mainstream] | none detected | 2.2×10⁻²⁹ W |
| Superconducting GW laser | [contested] | open | no positive result |
| BMV entanglement test | [mainstream] proposal | not done | W ≈ 1.16 |
| Parameter-free Diósi–Penrose | [refuted] | resolved | R₀ > 0.54×10⁻¹⁰ m |
| Single-graviton bar (Tobar) | [mainstream] proposal; proof [contested] | open | ~15 kg Be |
| Dyson undetectability | [mainstream] conjecture | open | Jupiter-mass cap |
| Superconducting transducers | [contested] | open, null | η < 1.6×10⁻⁵ |
| Wave energy and decay formulas | [mainstream] | resolved | 0.999963(63) |
| c⁵/G as a luminosity bound | [contested] | open, no theorem | observed peaks ~10⁻³ |
| "Dyson's bound" attribution | [refuted] | resolved | Dyson's own denial |
| Dyson Eq. 11 prefactor | [contested] | check print | 16× below Peters |
| GW200129 recoil | [mainstream] physics; event [contested] | open | 1542 km/s |
| Displacement memory | [mainstream] theory | undetected | 0.26 ± 4.1 vs 1 |
| (2,0)-mode evidence | [contested] | preprint | log₁₀B = 1.38 ± 0.79 |
| Spin memory | [mainstream] theory | undetected | no limit exists |
| Wave lensing | [mainstream] | none confirmed | GW231123 [contested] |
| Absorption by matter | [mainstream] | resolved | 2.6×10¹¹ Gpc in water |
| Gertsenshtein in transit | [mainstream] | resolved; 4× convention conflict | 8.3×10⁻⁵⁰ G⁻² cm⁻² |
| Pure scalar / vector polarisation | [mainstream] | excluded | log₁₀B = −14.7, −5.3 |
| Polarisation admixtures | [mainstream] question | open | log₁₀B ≈ 0 ± 0.6 |
| Primordial B-modes | [mainstream] bound | open | r < 0.034 |
| Nanohertz background | [mainstream] 3–4σ; origin [contested] | open | A = 2.4×10⁻¹⁵ |
| Siren H₀ | [mainstream] | open | 75.4 (+12.8/−9.1) |
| Earth normal-mode GW bound | [mainstream] | resolved | Ω_GW < 0.035–0.15 |

\* The notes give the eccentricity rate no explicit label; [contested] is this report's assignment.

The open items fall into four classes, and only some offer a real parallel track:

1. **Metrology:** the G spread and the Tate mass. Here the decisive move is an independent apparatus, not a new theory. NIST has just shown what that looks like. A clean G calculator, with full covariance, Birge expansion and one switch per input, is the right tool for evaluating each new datum.
2. **Statistics on populations:** wide binaries, the nanohertz background's origin, siren H₀, memory stacking and polarisation admixtures. The decisive move there is better data plus honest treatment of contamination. Both sides of the wide-binary dispute agree on the prediction, which makes it the most decidable open question in low-acceleration gravity.
3. **Undetected GR predictions:** memory, spin memory, primordial waves, entanglement induced by gravity, and single gravitons. These are not anomalies; they are queued experiments.
4. **Bounded parameters:** PPN, graviton mass, the inverse-square law and antimatter magnitude. These only tighten.

The fringe theories offer no parallel track. Each fails at a computable step: a neutrino count, a perihelion, a conservation law, or an energy budget eight to eleven orders short. Several items connect directly to this project's instruments:

- the Coughlin–Harms normal-mode bound uses superconducting-gravimeter series;
- the GRACE 2007 transient is a satellite-geodesy validation case;
- the LLR shielding bound turns any eclipse signal in a gravimeter record into an environmental flag;
- the calibrator formulas are Newtonian forward modelling with sub-percent published answers.

## Conclusion

The pattern worth keeping is the test, not any single anomaly. Every anomaly that died was killed by something a constant new force cannot fake: a second apparatus (G), an out-of-sample prediction (Juno), or a time dependence (Pioneer). Every anomaly still alive has either one apparatus and no repeat (Tate), or rests on statistics where contamination and signal look alike (wide binaries). That gives the crates a design rule: ship every claimed effect beside the replication or bound that tests it, and let the test assert the ratio. The rule also gives a fast triage for future claims. Ask whether anyone else has rebuilt the apparatus, whether the effect decays like its hardware, and whether it predicted an event it had not seen.

The second takeaway is about scale. Laboratory gravity that works is piconewton-scale Newtonian force, calibrated to 0.17%. Laboratory gravity that GR allows beyond that is 10⁻²¹ rad/s of frame dragging and 10⁻²⁹ W of radiation. Waves themselves carry kilowatts per square metre near a Galactic merger, but matter's cross-section, about 10⁻²⁵ m² Hz, wastes almost all of it. "Gravity manipulation" therefore fails on coupling, not energy. The honest frontier is measurement: smaller sources toward the Planck mass, shorter ranges below 38.6 µm, first detections of memory and gravitationally induced entanglement, and gravimeters reading Earth's normal modes as a wave detector. Two errors surfaced while writing this report, a 10³ slip in the notes' Yukawa factor at 1 au and a factor-of-4 conflict between two published Gertsenshtein formulas, on top of the 16× gap the notes found in Dyson's 1963 formula. That is a reminder that the test vectors below need checking against typeset sources before they become regression tests.

## What goes into the crates

All calculators are written clean-room from the cited formulas in SI units, so the MIT licence holds. Tolerances:

- formulas to 10⁻⁶ relative error against the script values;
- published comparisons to their stated rounding;
- each claim beside the bound it violates, following the rule from the earlier claims work.

Several entries reuse existing functions:

- `gwsens`: `lab::wave_flux`, `lab::rotor_power`, `hf::gertsenshtein_probability`, `inspiral::chirp_mass`;
- `gravitics`: `spacetime::quantum_inequality`, `casimir_energy_density`, `casimir_pressure`, `light_bending`, `superconductors::london_moment`, `ring_frame_dragging`, `propulsion::mach_effect_thrust`;
- `satgrav`: `atom::AtomInterferometer::phase`.

`geo5d` is not in this workspace, so check it for an existing Lense–Thirring function before adding one. Values marked *(checked here)* are this report's own recomputations.

### gwsens: gravitational waves

| Function (module) | Formula | Test vector | Label |
|---|---|---|---|
| `planck_luminosity()` (`energy`) | c⁵/G | 3.62825×10⁵² W; c⁵/4G = 9.07064×10⁵¹ W | mainstream unit, not a bound |
| `strain_from_luminosity(frac, f, r)` (`energy`) | h₀ = 2√2·√℘·c/(2πf·r) | (10⁻³, 250 Hz, 410 Mpc) → 1.34931×10⁻²¹; (1, 1 Hz, r = c·1 s) → 0.450158 | mainstream, derived |
| `circular_flux(eta, x)` (`energy`) | (32/5)η²x⁵ in units of c⁵/G | (¼, 1/6) → 5.14403×10⁻⁵ | mainstream |
| `flux_coefficients(l)` (`energy`) | a_ℓ, b_ℓ, a′_ℓ, b′_ℓ, a″_ℓ, b″_ℓ, c″_ℓ (Blanchet Table 1) | a₂ = 1/5, b₂ = 16/45, a₃ = 1/189, b₃ = 1/84, a₄ = 1/9072, b₄ = 4/14175, a′₂ = 2/5, a″₂ = 2/63, c″₂ = 16/45 | mainstream |
| `peters_enhancement(e)` (`energy`) | (1 + 73e²/24 + 37e⁴/96)/(1−e²)^(7/2) | e = 0.6171340 → 11.8568 | mainstream |
| `peters_pbdot(mp, mc, pb_days, e)` (`energy`) | −(192π/5)T☉^(5/3)(P_b/2π)^(−5/3) f(e) m_p m_c (m_p+m_c)^(−1/3) | B1913+16 (1.438, 1.390, 0.322997448918, 0.6171340) → −2.4021×10⁻¹² (published −2.40263(5)×10⁻¹², tolerance 10⁻³); J0737−3039 (1.338185, 1.248868, 0.1022515592973, 0.087777023) → −1.24781×10⁻¹² (tolerance 10⁻⁴) | mainstream |
| `energy_spectrum(chirp, f)` (`inspiral`) | dE/df = (π^(2/3)/3)G^(2/3)M_c^(5/3)f^(−1/3) | 1.188 M☉ at 900 Hz, ÷ 4π(40 Mpc)² → 2.66691×10⁻⁷ J m⁻² Hz⁻¹ | mainstream |
| `memory_0pn(m1, m2, d, x, theta)` (`memory`) | (2ηMx/R)·sin²Θ(17 + cos²Θ)/96, with M = G(m₁+m₂)/c² | (36, 29 M☉, 410 Mpc, 1/6, π/2) → 1.10657×10⁻²²; (10⁹, 10⁹ M☉, 1 Gpc) → 1.41236×10⁻¹⁵; Θ = 0 → 0 exactly | mainstream (lower bound) |
| `pta_memory_residual(t, t0, h0, psi, alpha, l)` (`memory`) | B·h₀[(t−t₀)Θ(t−t₀) − (t−tᵢ)Θ(t−tᵢ)], B = ½cos2ψ(1−cos α), tᵢ = t₀ + (l/c)(1+cos α) | B = ½, h₀ = 10⁻¹⁴, 5 yr before the pulsar term → 0.78894 µs | mainstream |
| `speed_bounds(delay, d, max_emission)` (`propagation`) | (delay·c/d, (delay − max_emission)·c/d) | (1.74 s, 26 Mpc, 10 s) → (+6.50199×10⁻¹⁶, −3.08658×10⁻¹⁵) | mainstream |
| `compton_wavelength(m_ev)` (`propagation`) | hc/(m_g c²) | 1.92×10⁻²³ eV → 6.45751×10¹⁶ m; 2.42×10⁻²³ → 5.12331×10¹⁶ m; 1.01×10⁻²⁴ → 1.22757×10¹⁸ m; 3.16×10⁻²³ → 3.924×10¹⁶ m (Bernus: 3.93×10¹³ km, to 1%) | mainstream |
| `massive_graviton_delay(m_ev, f, d)` (`propagation`) | (d/2c)(m_g c²/hf)²; with 1 − v/c = ½(c/(λ_g f))² | (1.92×10⁻²³, 100 Hz, 1 Gpc) → 1.1092×10⁻⁴ s; 1 − v/c = 1.078×10⁻²¹; doc test e^(−1 au/λ_g) = 1 − 2.3166×10⁻⁶ *(checked here)* | mainstream |
| `lens_parameter(m_lens, f)` (`propagation`) | w = 8πGM_Lz f/c³ | (100 M☉, 100 Hz) → 1.23791; w = 1 at 80.7813 M☉ | mainstream |
| `point_lens_geometric(y, m_lens)` (`propagation`) | μ± = ½ ± (y²+2)/(2y√(y²+4)); Δt = (4GM/c³)[y√(y²+4)/2 + ln((√(y²+4)+y)/(√(y²+4)−y))] | y = 1 → μ± = 1.17082, −0.17082; Δt = 4098.91 s at 10⁸ M☉ (4GM/c³ = 1970.2 s) | mainstream |
| `wave_optics_peak(w)` (`propagation`) | \|F(y=0)\|² = πw/(1 − e^(−πw)) | w = 1 → 3.28348 | mainstream |
| `viscous_damping(eta, omega, tau)` (`propagation`) | γ = 16πGη/c²; returns `None` if ωτ > 1 | water 10⁻³ Pa s → 3.7328×10⁻²⁹ s⁻¹ (c/γ = 2.60277×10¹¹ Gpc); air 1.8×10⁻⁵ → 6.71903×10⁻³¹; a mantle viscosity of 10²¹ Pa s at Hz frequencies must be refused | mainstream |
| `NEUTRINO_DAMPING` const (`propagation`) | power loss 0.356 for modes entering in radiation era | amplitude factor √0.644 = 0.802496 | mainstream |
| `gertsenshtein_coherent(b, l)` (`hf`) | sin²(√(4πG/μ₀)·B·L/c²), whose small-angle limit is the existing `gertsenshtein_probability` | 1 G, 1 cm → 8.26272×10⁻⁵⁰; mixing length at 1 G = 1.127 Mpc; 5 µG over 8 kpc → 1.2588×10⁻¹⁵ (Rothman–Boughn's form gives 3.1469×10⁻¹⁶, exactly 4× less) *(checked here)* | mainstream; resolve the factor of 4 first |
| `plasma_frequency(n_e)` (`hf`) | 8.98 kHz × √(n_e/cm⁻³) | 0.03 cm⁻³ → 1555.38 Hz (guard: no propagating photon below it) | mainstream |
| `graviton_count(h0, f)` (`quanta`) | h₀²c⁵/(32πGħω²) | (10⁻²¹, 150 Hz) → 3.85283×10³⁶ (published ≈ 4×10³⁶) | mainstream |
| `occupation_per_cell(h, f)` (`quanta`) | [c²h²ω²/(32πG)]/(ħω) × (c/f)³ | (10⁻²², 1 kHz) → 2.15031×10³⁵ (published 2×10³⁵); (7×10⁻²⁸, 10¹⁸ Hz) → 1.05365×10⁻⁵; equals 1 at 4.63715×10¹² Hz for h = 10⁻³⁰ | mainstream |
| `stimulated_rate(vs, m, h, l)` (`quanta`) | v_s²Mh²/(4l²π³ħ) | (5.4 km/s, 1800 kg, 5×10⁻²², 1) → 1.00326 Hz (published ≈ 1 Hz) | mainstream |
| `spontaneous_rate(rho, vs, r, l)` (`quanta`) | 8πGρv_s⁴R²/(Lc⁵) | (8570, 5 km/s, 0.5 m, 1 m) → 9.27562×10⁻³⁴ Hz (published ≈ 10⁻³³) | mainstream |
| `optimal_bar_mass(vs, h0, chirp, f)` (`quanta`) | (24π²/5)ħ/(h₀²v_s²)·(GM_c/2c³)^(5/3)·ω^(8/3) | Be (13 km/s, 2×10⁻²², 1.19 M☉, 100 Hz) → 12.8458 kg (paper ~15); Al (5.4 km/s, 150 Hz) → 219.501 kg (paper ~250) | mainstream; interpretation contested |
| `bar_cross_section(m, vs, l, r, poisson, n, theta)` (`quanta`) | (8/π)(M/n²)(G/c)(v/c)² sin⁴θ[1 + σ_P(1−2σ_P)(nπR/L)²/2] | LSU (2300 kg, 5.4 km/s, 3 m, 0.6 m, 0.345) → 4.31988×10⁻²⁵ m² Hz (published 4×10⁻²⁵); × `energy_spectrum` → 1.15207×10⁻³¹ J = 0.193189 ħω, Δv = 1.67083×10⁻⁴³ m/s | mainstream |
| `tidal_acceleration(h, f, l)` (`lab`) | ½(2πf)²hL | (10⁻²¹, 250 Hz, 1 m) → 1.2337×10⁻¹⁵ m/s²; Galactic GW150914 analogue → 6.32272×10⁻¹¹ per m | mainstream |
| `ncal_force_2f(m_test, m_cyl, r_q, d, phi, z)` (`lab`) | (9GMm r_q²/2d⁴)ζ^(−5/2)\|A e^(−iΦ) + B e^(−3iΦ)\|, A = 5/(6ζ) − 2/3, B = 5/(6ζ), ζ = 1 + (z/d)² | Ross Table I (39.680, 1.0558, 0.06033, 1.1802, 52.24°, 0.010) → 19.07 pN (paper's point mass 19.04, brute force 19.10) | mainstream; formula reconstructed |
| `ncal_force_3f(m_test, m_cyl, r_h, d, phi, z)` (`lab`) | (15GMm r_h³/2d⁵)ζ^(−7/2)\|A e^(−2iΦ) + B e^(−4iΦ)\|, A = 7/(8ζ) − 3/4, B = 7/(8ζ) | same table, r_h = 0.10478 → 8.978 pN (paper 8.97, brute force 9.017) | mainstream; formula reconstructed |
| `free_mass_displacement(force, mass, f)` (`lab`) | F/(M(2πf)²), valid well above the suspension resonance | (19.16 pN, 39.68 kg, 15.60 Hz) → 5.026×10⁻¹⁷ m; ÷ 3995 m → 1.258×10⁻²⁰ | mainstream |
| `virgo_ncal_rp(b, rho, alpha, rmax, rmin, d, phi, z)` (`lab`) | R·P, with R = Gbρ sinα(r_max⁴ − r_min⁴)/(32π²) and P = [(15/Z−6)²cos²Φ + 36 sin²Φ]^(1/2) Z^(−5/2) d^(−4) | (0.052, 2800, 45°, 0.190, 0.055, 1.32, 0.241, 0.643) → 3.331×10⁻¹⁴ m Hz² (published 3.3×10⁻¹⁴) | mainstream |
| `weinberg_efficiency(m, vs, l, gamma_heat)` (`lab`) | 64GMv_s⁴/(15L²c⁵Γ) | (1000 kg, 5 km/s, 1.5 m, 2π·10³/10⁶) → 5.199×10⁻³³ (Kiefer–Weber ≈ 10⁻³⁴) | mainstream |
| `inflation_scale(r, a_s)` (`sgwb`) | V^(1/4) = (3π²A_s r/2)^(1/4)M_pl; H = πM_pl√(A_s r/2) | A_s = 2.0989×10⁻⁹, M_pl = 2.43532×10¹⁸ GeV: r = 0.036 → (1.40841×10¹⁶, 4.70261×10¹³) GeV; 0.034 → (1.38842×10¹⁶, 4.57012×10¹³); 0.01 → (1.02247×10¹⁶, 2.47849×10¹³) | mainstream |
| `inflation_omega(r, a_s, omega_rad_h2)` (`sgwb`) | Ω_rad h²·rA_s/24 (k ≫ k_eq) | (0.036, 2.0989×10⁻⁹, 4.17867×10⁻⁵) → 1.31559×10⁻¹⁶ (upper bound) | mainstream |

Tests only, using existing functions:

- `rotor_power` reproduces the MTW beam at 2.234×10⁻²⁹ W (I₁ = M(L²/12 + r²/4), I₂ = Mr²/2).
- `wave_flux(1e-22, 1e3)` returns 1.58532 mW/m².
- The GW150914 PRL convention c³ḣ²/(16πG) is exactly twice `wave_flux`, giving 0.0198165 W/m².
- Dyson's Eq. 11, 128V¹⁰/(5Gc⁵), gives 2.07316×10³⁷ erg/s against Peters's 3.31706×10³⁸ erg/s, a ratio of 16. This is a documented counter-example.
- Luminosity ratios: GW150914 at 9.64651×10⁻⁴ c⁵/G and GW190521 at 1.01977×10⁻³ c⁵/G.
- Cross-sections: ℓ_P² = 2.61228×10⁻⁶⁶ cm², ℓ_P²/σ_T = 3.92679×10⁻⁴², and a Jupiter-mass cap of 2.29715×10²⁷ kg.

### geo5d: relativity and PPN

| Function (module) | Formula | Test vector | Label |
|---|---|---|---|
| `perihelion_advance(a, e, period, gm, gamma, beta)` (`ppn`) | 6πGM/(c²a(1−e²)) · (2 + 2γ − β)/3 per orbit; the J₂ term is added separately as 3×10⁻⁴(J₂/10⁻⁷) | Mercury (5.7909050×10¹⁰ m, 0.205630, 87.9691 d, GM☉, 1, 1) → 42.98″/century; (γ = ½, β = 1) → 28.65″ | mainstream |
| `limb_deflection(gamma, gm, radius)` (`ppn`) | ½(1+γ)·4GM/(c²R) | γ = 1, R = 6.957×10⁸ m → 1.751″ (Will's 1.7505″ uses a larger R☉) | mainstream |
| `shapiro_round_trip(gamma, gm, r1, r2, b)` (`ppn`) | 2(1+γ)(GM/c³) ln(4r₁r₂/b²) | (1, GM☉, 1 au, 8.43 au, 1.6 R☉) → 262.4 µs | mainstream |
| `brans_dicke_gamma(omega)` (`ppn`) | (1+ω)/(2+ω) | ω = 4×10⁴ → γ − 1 = −2.4998×10⁻⁵ | mainstream, constrained |
| `brans_dicke_omega_min(gamma_minus_one_floor)` (`ppn`) | −1/(γ−1) − 2 | Cassini 2σ edge −2.5×10⁻⁵ → 4.0×10⁴ | mainstream |
| `nordtvedt(beta, gamma)` (`ppn`) | 4β − γ − 3 | (1, Brans–Dicke at 4×10⁴) → 2.4998×10⁻⁵ | mainstream |
| `PPN_LIMITS` const (`ppn`) | Cassini, MESSENGER and LLR values, plus Will's Table 4 | γ−1 (2.1 ± 2.3)×10⁻⁵; β−1 (−1.6 ± 1.8)×10⁻⁵; ξ 4×10⁻⁹; α₁ 7×10⁻⁵; α₂ 2×10⁻⁹; α₃ 4×10⁻²⁰; ζ₃ 10⁻⁸ | mainstream data |
| `gyro_precession(j, r)` (`frame_dragging`) | (G/(c²r³))[3(J·r̂)r̂ − J] | tungsten sphere (R = 0.5 m, M = 1.011×10⁴ kg, ω = 455.3 rad/s) at its pole → 5.466×10⁻²¹ rad/s; Ω/ω = 1.201×10⁻²³ | mainstream |
| `shell_interior(m, r, omega)` (`frame_dragging`) | 4GMω/(3c²R) | (10⁴ kg, 1 m, 300 rad/s) → 2.97×10⁻²¹ rad/s | mainstream |

Doc test: Gravity Probe B's 39.2 mas/yr is 6.022×10⁻¹⁵ rad/s. The Damour–Esposito-Farèse formulas γ − 1 = −2α₀²/(1 + α₀²) and β − 1 = β₀α₀²/[2(1 + α₀²)²] belong here only after checking them against a fetched source.

### gravitics: frontier and exotic physics

| Function (module) | Formula | Test vector | Label |
|---|---|---|---|
| `gls_mean(values, cov)` → (mean, u, χ²) (`anomalies`) | Ḡ = (1ᵀC⁻¹y)/(1ᵀC⁻¹1), u = (1ᵀC⁻¹1)^(−½), χ² = rᵀC⁻¹r | CODATA 2022's 16 inputs with r = 0.351, 0.134, 0.068 → 6.674300, 0.000038 (5.7 ppm), χ² = 195.7 (ν = 15) | mainstream |
| `birge_expand(u, chi2, dof, factor)` (`anomalies`) | R_B = √(χ²/ν); u·k; χ²/k² | Birge ratio 3.61; k = 3.9 → 0.00015 (22.1 ppm), χ² 12.9; residuals +7.75 (BIPM-14), −6.80 (JILA-18), +4.78 (BIPM-01) | mainstream |
| `CODATA_2022_G` const (`anomalies`) | the 16-row table of identifier, value and uncertainty | spread 0.00368, i.e. 551 ppm | mainstream data |
| `flyby_anderson(v_inf, dec_in, dec_out)` (`anomalies`) | K·V∞(cos δᵢ − cos δₒ), K = 2ω⊕R⊕/c = 3.0993×10⁻⁶ | Galileo I +4.123, Galileo II −4.674, NEAR +13.279, Cassini −1.068, Rosetta +2.067, MESSENGER +0.055 mm/s | contested; doc must state that Juno was predicted at +6 and observed 0 ± 0.8 |
| `column_density(radius, density)` (`anomalies`) | 2Rρ | Earth (6371 km, 5514) → 7.026×10¹⁰ kg/m²; Moon (1737.4 km, 3344) → 1.162×10¹⁰ | mainstream geometry |
| `shielding_loss(h, column)` (`anomalies`) | 1 − e^(−hΣ) | Majorana 6.7×10⁻¹³: Earth 4.598×10⁻², Moon 7.755×10⁻³; LLR 10⁻²²: 7.026×10⁻¹², 1.162×10⁻¹²; eclipse bound 10⁻²²·Σ☾·5.9301×10⁻³ = 6.9×10⁻¹⁵ m/s² | refuted (Majorana) |
| `pioneer_thermal(a0, tau, t)` (`anomalies`) | a₀·2^(−t/τ) | (10.1×10⁻¹⁰, 28.8 yr, 20 yr) → 6.24×10⁻¹⁰ m/s² | mainstream |
| `tate_excess(measured, theory, sigma)` (`anomalies`) | (m − m_th)/m_th and its significance; B_g/ω = 2·excess | (1.000084, 0.999992, 21 ppm) → 92 ppm, 4.4σ, B_g/ω = 1.84×10⁻⁴; pair with `london_moment(1.0)` = 1.1371×10⁻¹¹ T | contested |
| `rar(g_bar, g_dagger)` (`modified`) | g_bar/(1 − e^(−√(g_bar/g†))) | g† = 1.20×10⁻¹⁰: 10⁻¹² → 1.1462×10⁻¹¹; 10⁻¹¹ → 3.9881×10⁻¹¹; 10⁻¹⁰ → 1.6705×10⁻¹⁰; 1.2×10⁻¹⁰ → 1.898×10⁻¹⁰; 10⁻⁹ → 1.0590×10⁻⁹; 10⁻⁸ → 1.0001×10⁻⁸; 10⁻¹⁴ within 1% of √(g_bar g†) | mainstream data, contested theory |
| `btfr_velocity(m, a0)` (`modified`) | (GMa₀)^(1/4) | 10¹¹ M☉, 1.2×10⁻¹⁰ → 199.8 km/s | contested |
| `mond_radius(m, a0)` (`modified`) | √(GM/a₀) | 2 M☉ → 9942 au; 1 M☉ → 7030 au; 10¹¹ M☉ → 10.8 kpc | contested |
| `einstein_cartan_density(m)` (`modified`) | m²c⁴/(8πGħ²) | electron → 3.593×10⁴⁷ g/cm³; neutron → 1.215×10⁵⁴ g/cm³ | mainstream |
| `weber_perihelion_ratio(gamma_c2)` (`modified`) | Δω/Δω_GR = γc²/6, independent of μ | Tisserand (γ = 2) → 1/3; Weber (1) → 1/6; Assis (6) → 1 (the RK4 script gives 0.3333, 0.1667, 0.9998) | refuted (cannot also give 2γ − μ = 3/c²) |
| `weber_mass_shift(q, u, m)` (`modified`) | qU/(3mc²) | electron at 10 kV → 6.523×10⁻³ | refuted by experiment |
| `force_factor(r, alpha, lambda)` (`yukawa`) | 1 + α(1 + r/λ)e^(−r/λ) | (λ, 1, λ) → 1.736 | mainstream |
| `dark_energy_length(rho)` (`yukawa`) | (ħc/ρ)^(1/4) | 3.8 keV/cm³ → 84.89 µm (published ≈ 85) | mainstream |
| `boson_mass_ev(lambda)` (`yukawa`) | ħc/λ | 38.6 µm → 5.112 meV (published 5.1) | mainstream |
| `alpha_from_coupling(g)` (`yukawa`) | g²ħc/(4πGm_N²) | 5×10⁻¹⁴ → 3.368×10¹⁰ | mainstream; convention to check |
| `ISL_LIMITS` const (`yukawa`) | (λ, \|α\|) pairs | (38.6 µm, 1); (56 µm, 1, superseded); (10 µm, 1.4×10⁴) | mainstream data |
| `bmv_phases(m, d, dx, tau)` (`quantum`) | Gm²τ/(ħr) at r = d, d − Δx, d + Δx | (10⁻¹⁴ kg, 450 µm, 250 µm, 2.5 s) → 0.352, 0.791, 0.226 rad; Δφ_RL = 0.440, Δφ_LR = −0.126 | mainstream proposal |
| `bmv_witness_optimum(d_lr, d_rl)` (`quantum`) | 1 + \|sin((Δφ_LR + Δφ_RL)/2)\| | (−0.126, 0.440) → 1.156 (published ≈ 1.16); (−0.2, 0.7) → 1.247 | mainstream; never ship the printed σx⊗σz − σy⊗σz form |
| `dp_photon_rate(r0, z, n_atoms, e_lo, e_hi)` (`quantum`) | ∫ 2Ge²Z²N_a/(3π^(3/2)ε₀c³R₀³ω) dω = prefactor·ln(e_hi/e_lo) | Ge 375 cm³ (Z = 32, N_a = 1.655×10²⁵), 1–3.8 MeV, unit efficiency: R₀ = 0.54×10⁻¹⁰ → 10.67/day; 0.05×10⁻¹⁰ → 1.345×10⁴/day | refuted (parameter-free model) |
| `penrose_time_sphere(m, r, d)` (`quantum`) | ħ/[GM²(6/(5R) − 1/d)], d ≥ 2R | (10⁻¹⁴ kg, 1 µm, 10 µm) → 14.36 ms | contested; normalisation differs from Donadi's |
| `forward_coupling()` (`spacetime`) | G/c² | 7.426×10⁻²⁸ m/kg | mainstream |
| `induced_acceleration(m, a_flow, r)` (`spacetime`) | 4Gma_f/(c²r) | (10⁴ kg, 10⁴ m/s², 1 m) → 2.97×10⁻¹⁹ m/s²; 9.8 m/s² at a_f = 10⁶ needs m/r = 3.299×10²¹ kg/m | mainstream; scaling is an inference |
| `mach_effect_gr_bound(m, l, f)` (`propulsion`) | (GM/(Lc²))·(2πfL/c)² | (0.1 kg, 1 cm, 50 kHz) → 8.155×10⁻³⁷; ship beside `mach_effect_thrust` | mainstream (Rodal) |

Doc tests only (no public function):

- **Antimatter:** ALPHA-g σ = 0.206 g, giving 1.21σ and 8.5σ.
- **Equivalence principle:** MICROSCOPE quadrature 2.75×10⁻¹⁵.
- **NIST against BIPM-14:** −282 ppm (−7.3σ) and −250 ppm (−4.0σ).
- **Wide binaries:** 10^(2Γ) = 1.479, 1.337, 1.600, 1.837 and 1.380; Chae 2023's 10^(√2δ) = 1.426 and 1.117.
- **Verlinde:** cH₀/6 = 1.091×10⁻¹⁰ at H₀ = 67.4.
- **Sciama:** 3/(8π) = 0.1194, 3/2 and ½.
- **Vacuum energy:** `quantum_inequality(1e-9)` = −3.718×10⁻²⁶; Casimir energy per area π²ħc/(720d³) = 4.334×10⁻¹⁰ J/m² at 1 µm; Moddel's 2.016×10⁻⁸ J against 2.412×10⁻¹⁹ J (ratio 8.359×10¹⁰); dark energy 5.245×10⁻¹⁰ J/m³.
- **Superconductors:** Z_g = 4πG/c = 2.798×10⁻¹⁸; e²/(4πε₀Gm_e²) = 4.166×10⁴².
- **Small masses:** Westphal's static force 8.921×10⁻¹⁴ N and modulation 2Gm_s d_m/d₀³ = 2.655×10⁻¹⁰ m/s²; Planck mass 21.76 µg.
- **Aharonov–Bohm scale:** m_Rb·GM/(ħr) = 1.631 rad/s at 7 cm.
- **Rosi's single-diffraction phase:** k_eff T² = 4.123×10⁵ rad per m/s², which is half of what `satgrav::atom::AtomInterferometer::phase` returns, because that function assumes double diffraction.

### Verify before coding

| Item | Problem | Action |
|---|---|---|
| Ross Eqs. 2–3 | fractions reconstructed from garbled extraction | check the typeset PDF; the 0.2% match to Table II is reassuring but not proof |
| BMV witness | printed operator certifies nothing | check the published PRL |
| Penrose sphere time | normalisation differs from Donadi's −8πG form | rederive before shipping |
| α from coupling g | 4π and nucleon-mass convention | fix the convention against the neutron-star paper |
| Gertsenshtein sin² form | Rothman–Boughn is 4× below Pshirkov and the 2025 Review | reconcile before adding `gertsenshtein_coherent` |
| Dyson Eq. 11 | 16× below Peters | check the printed book |
| Aguiar bar ratio | 4.5×10²⁴ published vs 8.3×10²⁴ computed with area 2RL | use the cross-section, not the ratio, as the test |
| Tobar optimal masses | 12.8 kg vs ~15 kg; 220 kg vs ~250 kg | test the formula to 10⁻⁶ and the paper to 25% |
| NIST 2026 G | 6.67366(20) vs 6.67387 (57 ppm) | do not hard-code either until the Metrologia text is read |
| Fuchs 2024 force | 1030 aN vs 30 aN | exclude from tests |
| Estevez sector mass | 1.89 kg computed vs ~2.06 kg printed | keep the R·P test, skip mass |
| Yukawa factor at 1 au | notes say 1 − 2.3×10⁻⁹ | use 1 − 2.3166×10⁻⁶ *(checked here)* |
| GWTC-3 graviton mass | 1.27, 2.42 or 2.23×10⁻²³ eV | cite the version in the constant's doc |
| DEF and vDVZ γ = ½ | background knowledge, not fetched | source before coding |
| Tate theory value | 0.999992 from a snippet | source Tate 1989 before fixing the test |

### What stays out

Leave out the following:

- models of Heim gravitophotons, torsion-field generators, Weber electrodynamics beyond the two refutation checks, zero-point-field inertia, Moddel's cavity current, Chiao's order-unity transducer, and the superconducting gravitational-wave laser, none of which yields a derivable number to test;
- the Anderson flyby formula as physics: keep it only as a labelled reproduction beside the Juno null;
- Verlinde beyond the a₀ doc test;
- the full wave-optics lensing factor F(f), which needs a confluent hypergeometric function of complex argument and should wait until a tested ₁F₁ exists;
- the six polarisation antenna patterns and the Coughlin–Harms normal-mode estimator, deferred because the notes did not extract their formulas, though the latter is the most valuable gravimetry-to-waves link once they are sourced.
