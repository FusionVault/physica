# Deep physics of gravitational waves: what they carry, memory, propagation, gravitons, primordial waves, sirens and limits

Scope: everything except detector noise curves, inspiral SNR, stochastic power-law curves and HFGW limits, which are in `reports/Gravitational waves and gravity beams.md`. SI units unless stated. Constants: G = 6.67430×10⁻¹¹, c = 299 792 458, ħ = 1.054571817×10⁻³⁴, GM☉ = 1.32712440018×10²⁰ m³ s⁻², 1 Mpc = 3.0857×10²² m. Every value marked **computed** was recomputed with the stdlib-only script `deep_gravitational_wave_checks.py` in this folder (`python3 -I deep_gravitational_wave_checks.py`; it ends with assert-based checks against the published numbers). Labels: **[mainstream]** = standard GR or a collaboration measurement; **[contested]** = disputed, preprint-only, or model-dependent; **[refuted]** = wrong or retracted.

## 1. Energy, momentum and angular momentum carried by GWs

### Takeaway
The energy, momentum and angular momentum a GW carries are fixed by the Isaacson effective stress tensor and its multipole expansion. Those formulas are verified to 10⁻⁴ by binary pulsars, and they give ≈3 M☉c² for GW150914. Comparable-mass black-hole mergers peak at about 10⁻³ c⁵/G whatever their mass. The "Dyson c⁵/G luminosity bound" is a misattribution: c⁵/G is a unit, not a proven limit. Dyson's own 1963 formula turns out to be 16× low against the quadrupole formula (computed here).

### Cited Findings
- **[mainstream] Isaacson effective stress tensor.** In TT gauge, T^GW_ab = (1/32π)⟨∇_a h_cd ∇_b h^cd⟩ with G = c = 1, Eq. 5.38 of [Flanagan & Hughes 2005](https://arxiv.org/abs/gr-qc/0501041); in SI, t_μν = (c⁴/32πG)⟨∂_μh^TT_ij ∂_νh^TT_ij⟩. For the plane wave h_xx = −h_yy = h₀cos(ωt−ωz) the averaged energy density equals the energy flux, T^tt = T^tz = h₀²ω²/32π (Eq. 5.40), i.e. **F = c³h₀²ω²/(32πG)** in SI. Published test: F = 1.5 mW m⁻² (h₀/10⁻²²)²(f/1 kHz)², "comparable to the flux of reflected sunlight from a full moon" (Eq. 5.41). **Computed: 1.585 mW m⁻²**.
- **[mainstream] Flux in terms of the two polarisations.** dE/(dA dt) = (c³/16πG)⟨ḣ₊² + ḣ×²⟩. The GW150914 properties paper uses ≈ c³|ḣ|²/(16πG) ∼ 10⁵ erg s⁻¹ m⁻² for |h| ≈ 10⁻²¹ at 250 Hz, which gives a luminosity of ∼10⁵⁶ erg/s at the inferred distance ([LVC, PRL 116, 241102](https://arxiv.org/abs/1602.03840)). **Computed: 0.0198 W m⁻²** (= 1.98×10⁵ erg s⁻¹ m⁻²). The peak luminosity spread over 4π(410 Mpc)² gives **0.0174 W m⁻²** (computed).
- **[mainstream] Multipole expansion of the fluxes** ([Blanchet, Living Rev. Relativ., Eqs. 179–181 and Table 1](https://arxiv.org/abs/1310.1528), after Thorne 1980). The energy flux is F = Σ_{ℓ≥2} G/c^{2ℓ+1} [a_ℓ U⁽¹⁾_L U⁽¹⁾_L + (b_ℓ/c²) V⁽¹⁾_L V⁽¹⁾_L]. Here U_L and V_L are the radiative mass- and current-type moments; at leading order U_L = d^ℓI_L/dt^ℓ and V_L = d^ℓJ_L/dt^ℓ.
  - Energy flux: a_ℓ = (ℓ+1)(ℓ+2)/[(ℓ−1)ℓ ℓ!(2ℓ+1)!!] and b_ℓ = 4ℓ(ℓ+2)/[(ℓ−1)(ℓ+1)!(2ℓ+1)!!].
  - Angular-momentum flux: G_i = ε_iab Σ G/c^{2ℓ+1}[a′_ℓ U_aL−1 U⁽¹⁾_bL−1 + (b′_ℓ/c²) V_aL−1 V⁽¹⁾_bL−1], with a′_ℓ = (ℓ+1)(ℓ+2)/[(ℓ−1)ℓ!(2ℓ+1)!!] and b′_ℓ = 4ℓ²(ℓ+2)/[(ℓ−1)(ℓ+1)!(2ℓ+1)!!].
  - Linear-momentum (recoil) flux: a″_ℓ = 2(ℓ+2)(ℓ+3)/[ℓ(ℓ+1)!(2ℓ+3)!!], b″_ℓ = 8(ℓ+3)/[(ℓ+1)!(2ℓ+3)!!] and c″_ℓ = 8(ℓ+2)/[(ℓ−1)(ℓ+1)!(2ℓ+1)!!] (the c″ term couples mass and current moments).
  - **Computed test values:** a₂ = 1/5, which is the quadrupole formula P = (G/5c⁵)⟨I⃛_ij I⃛_ij⟩. Also b₂ = 16/45, a₃ = 1/189, b₃ = 1/84, a₄ = 1/9072, b₄ = 1/3543.75 and a′₂ = 2/5, the latter giving dJ_i/dt = (2G/5c⁵)ε_iab⟨Ï_ac I⃛_bc⟩. For recoil, a″₂ = 2/63 and c″₂ = 16/45.
- **[mainstream] Circular-binary flux at leading order:** P = (32/5)(c⁵/G)η²x⁵, with x = (GMω/c³)^{2/3} and η = m₁m₂/M². **Computed:** an equal-mass binary at the ISCO (x = 1/6) gives 5.14×10⁻⁵ c⁵/G (the general formula is in [Blanchet LRR](https://arxiv.org/abs/1310.1528)).
- **[mainstream] Peters–Mathews orbital decay, tested by the Hulse–Taylor pulsar.** Ṗ_b = −(192π/5)T☉^{5/3}(P_b/2π)^{−5/3} f(e) m_p m_c (m_p+m_c)^{−1/3}, with T☉ = GM☉/c³ and f(e) = (1 + 73e²/24 + 37e⁴/96)/(1−e²)^{7/2}.
  - PSR B1913+16 ([Weisberg & Huang 2016](https://arxiv.org/abs/1606.02744)): m_p = 1.438 ± 0.001, m_c = 1.390 ± 0.001 M☉, e = 0.6171340(4), P_b = 0.322997448918(3) d. The observed Ṗ_b is −2.423(1)×10⁻¹² and the Galactic term −0.025(4)×10⁻¹², so the intrinsic value is −2.398(4)×10⁻¹². The GR value is **−2.40263(5)×10⁻¹²**, and the ratio observed/GR is **0.9983 ± 0.0016**.
  - **Computed:** f(e) = 11.857 and Ṗ_b(GR) = −2.4021×10⁻¹². This agrees to 2×10⁻⁴, which is the limit set by the masses being quoted to only 10⁻³.
- **[mainstream] Double pulsar J0737−3039A/B** ([Kramer et al., PRX 11, 041050 (2021)](https://arxiv.org/abs/2112.06795)): P_b = 0.1022515592973(10) d, e_T = 0.087777023(61), Ṗ_b(obs) = −1.247920(78)×10⁻¹², m_A = 1.338185(+12/−14) M☉, m_B = 1.248868(+13/−11) M☉ and Ṗ_b(ext) = −1.68(+11/−10)×10⁻¹⁶. The result is **Ṗ_b^GW/Ṗ_b^GW,GR = 0.999963(63)**, which confirms GW emission to 1.3×10⁻⁴ at 95% confidence. **Computed** Peters value: −1.24781×10⁻¹² (use as a test to about 10⁻⁴).
- **[mainstream] GW150914 energetics** ([PRL 116, 241102](https://arxiv.org/abs/1602.03840)): E_rad = **3.0 (+0.5/−0.4) M☉c² = 5.3 (+0.9/−0.8)×10⁴⁷ J**. The peak luminosity was **3.5 (+0.5/−0.4)×10⁵⁶ erg s⁻¹ = 200 (+30/−20) M☉c² s⁻¹**, and about 5% of the total mass-energy was radiated, roughly half of it in the merger. The discovery abstract gives 3.0 ± 0.5 M☉c², a peak strain of 1.0×10⁻²¹, matched-filter SNR 24, a false-alarm rate below 1 per 203 000 yr and significance above 5.1σ ([PRL 116, 061102](https://arxiv.org/abs/1602.03837)). **Computed:** 3.0 M☉c² = 5.361×10⁴⁷ J, 200 M☉c² s⁻¹ = 3.574×10⁴⁹ W, and L_peak/(c⁵/G) = **9.65×10⁻⁴**.
- **[mainstream] GW190521:** E_rad = **7.6 (+2.2/−1.9) M☉c²** and ℓ_peak = **3.7 (+0.7/−0.9)×10⁵⁶ erg s⁻¹** ([LVC, ApJL 900, L13](https://arxiv.org/abs/2009.01190)). The component masses are 85 and 66 M☉, the remnant 142 M☉ and the distance 5.3 Gpc ([PRL 125, 101102](https://arxiv.org/abs/2009.01075)). **Computed:** ℓ_peak/(c⁵/G) = 1.02×10⁻³.
- **[mainstream] GW170817:** E_rad **> 0.025 M☉c²** (this is the inspiral energy up to about 600 Hz, where tidal effects become strong) at D_L = 40 (+8/−14) Mpc ([PRL 119, 161101](https://arxiv.org/abs/1710.05832)). **Computed:** > 4.47×10⁴⁵ J.
- **[mainstream] GW250114, the loudest event to date:** network SNR 80 (77–80 across pipelines), with m₁ = 33.6 and m₂ = 32.2 M☉, spins ≤ 0.26 and e ≤ 0.03. Its strain amplitude of ∼10⁻²¹ is similar to GW150914's. The remnant's horizon area exceeds the initial total, as Hawking's area law requires, even with up to five of the strongest merger cycles excised, and the ringdown frequencies are constrained to ±30% of Kerr ([LVK, PRL 2025, arXiv:2509.08054](https://arxiv.org/abs/2509.08054)). That paper quotes SNR 26 for GW150914, against 24 in the discovery PRL, a convention difference that should be flagged in any test.
- **[mainstream as a unit; the bound is contested] Planck luminosity.** L_P = c⁵/G = **3.628×10⁵² W** and c⁵/4G = 9.07×10⁵¹ W (computed). Cardoso, Ikeda, Moore and Yoo find in full nonlinear simulations that suitable initial data can give arbitrarily large luminosity. Requiring no past horizon, however, seems to limit luminosity to below c⁵/G (a thought experiment gives L ≲ L_P/2). Critical collapse gives the largest value simulated so far, ≈0.2 L_P. Their table lists inspiralling black holes at ∼7×10⁴⁹ W (2×10⁻³ L_P) and a gamma-ray burst at ∼3×10⁴⁷ W (10⁻⁵ L_P) ([PRD 97, 084013 (2018)](https://arxiv.org/abs/1803.03271)). Jowsey and Visser find that Vaidya spacetimes and an evaporating constant-density star also give arbitrarily large luminosity. They conclude that any bound must depend on the source's internal structure, and they list strong, medium and weak conjectures of ℘ ≤ 1/4, ℘ ≤ 1/2 and ℘ = O(1) in L = ℘c⁵/G ([Jowsey & Visser 2021](https://arxiv.org/abs/2105.06650)).
- **[refuted as attribution] "Dyson's bound c⁵/G".** Jowsey and Visser reproduce a note from Dyson in which he denies ever proposing c⁵/G as a luminosity limit. He traces the idea to his 1962 essay "Gravitational Machines" (Chapter 12 of *Interstellar Communication*, 1963), whose Eq. 11 gives W = 128V¹⁰/(5Gc⁵) for two equal masses in circular orbit. That expression tends to 128c⁵/5G as V → c, independent of mass ([Jowsey & Visser 2021](https://arxiv.org/abs/2105.06650)).
- **[contested: a prefactor discrepancy found here] Dyson's Eq. 11 does not match the quadrupole formula.** Dyson's chapter defines V = (GM/4R)^{1/2} with R the orbital radius of each star, and takes M = 1 M☉, R = 2×10⁹ cm, V = 1.3×10⁸ cm/s, a period of 100 s and W = 2×10³⁷ erg/s ([Dyson 1963 scan](https://www.gwern.net/doc/science/1963-dyson.pdf)). **Computed:**
  - His formula with his own numbers reproduces 2.07×10³⁷ erg/s.
  - The Peters quadrupole formula for the same binary (separation 4×10⁹ cm) gives **3.32×10³⁸ erg/s = (2048/5)V¹⁰/(Gc⁵)**, so Dyson's prefactor is 16× low.
  - Both forms are invalid as V → c, so neither "128c⁵/5G" nor "2048c⁵/5G" is a physical limit.
- **[mainstream; event-specific evidence contested] GW recoil.** For GW200129, Varma et al. constrain the remnant kick to v_f ∼ **1542 (+747/−1098) km/s**, or v_f ≳ 698 km/s as a one-sided limit at 90%. This is the first large kick measured for an individual event, and kicks up to about 5000 km/s are possible in general ([PRL 128, 191102 (2022)](https://arxiv.org/abs/2201.01302)). **Computed:** v/c = 5.1×10⁻³, so the radiated linear momentum is about 0.5% of M_f c.

### Inferences
- The peak luminosity of a comparable-mass BBH is scale-free. GW150914 (65 M☉) reached 9.65×10⁻⁴ L_P and GW190521 (150 M☉) 1.02×10⁻³ L_P (computed), because L ∝ (c⁵/G)·F(η, spins) with no mass scale. A Rust `peak_luminosity()` should therefore return a multiple of c⁵/G, not a value in watts per solar mass.
- The leading-order flux at the ISCO (5.1×10⁻⁵ L_P) underestimates the merger peak (≈10⁻³ L_P) by about 20× (computed). Any tool that extrapolates the quadrupole formula through merger will be badly wrong. Use NR-fit numbers for the peak.
- Good regression tests for a `gw_energy` crate: the Blanchet coefficients (exact rationals), the Flanagan–Hughes 1.5 mW m⁻² flux, and the two pulsar Ṗ_b values. Add a Dyson-Eq.-11 test as a documented counter-example (a factor of 16).

### Gaps
- Isaacson 1968 and Thorne 1980 (Rev. Mod. Phys.) are paywalled. Their formulas were taken from Flanagan & Hughes and Blanchet's open reviews, and Maggiore's textbook was not accessible.
- The task's "3.6×10⁵⁶ erg/s" for GW150914 matches common quotations of the discovery paper, but the abstract fetched here did not state it. The properties PRL gives 3.5 (+0.5/−0.4)×10⁵⁶.
- No published value for the angular momentum radiated by GW150914 was found.
- A later data-quality reanalysis of GW200129's precession and kick evidence is known to exist but was not re-read, so treat the kick as [contested].
- Dyson's original Eq. 11 was read from a scanned copy, so a typographical exponent cannot be excluded. The factor of 16 should be checked against the printed book.

## 2. Nonlinear (Christodoulou) memory, displacement memory and spin memory

### Takeaway
GR predicts a permanent, + polarised offset (memory) of about 0.05–0.1 ηM/R, which vanishes face-on. As of October 2026 no memory of any type has been detected. LVK stacking of GWTC-3 gives ln BF = 0.01, and a hierarchical analysis of GWTC-5.0 measures the memory amplitude as 0.26 ± 4.1 (GR = 1). Estimates of the catalogue needed range from about 170 events (full (2,0) mode, optimistic) to about 2000 (memory term alone). Pulsar-timing limits on memory bursts are no better than about 10⁻¹⁴.

### Cited Findings
- **[mainstream] The nonlinear memory integral** (Blanchet–Damour, Christodoulou 1991; Thorne): δh^TT_jk = (4/R)∫_{−∞}^{T_R}dt′ ∫ (dE^gw/dt′dΩ′) [n′_j n′_k/(1 − n′·N)]^TT dΩ′ with G = c = 1 (multiply by G/c⁴ in SI). The memory is hereditary, depending on the whole past of the source, and Thorne showed that it is the linear memory of the radiated gravitons treated as unbound particles ([Favata 2010, CQG 27, 084036, Eq. 12](https://arxiv.org/abs/1003.3486)).
- **[mainstream] Post-Newtonian inspiral memory** ([Favata 2010, Eqs. 17–18](https://arxiv.org/abs/1003.3486)): h₊^(mem) = (2ηMx/R) Σ_n x^{n/2} H₊^{(n/2,mem)}, with the leading term **H₊^(0) = (1/96) sin²Θ(17 + cos²Θ)**. The 0.5PN and 1.5PN terms vanish, and memory appears only in the + polarisation. The merger memory follows from the (2,2) mode alone, **h₊^(mem) ≈ (R/192π) sin²Θ(17 + cos²Θ)∫|ḣ₂₂|² dt** (Eq. 21). Favata estimated that Advanced LIGO could see memory from single events only out to ≲ 20 Mpc, while LISA could see it to z ≲ 2.
- **Computed test values (0PN formula truncated at x = 1/6, edge-on), which are lower bounds because merger adds more:**
  - GW150914 (36 + 29 M☉, 410 Mpc): ηM/R = 1.87×10⁻²¹ and h_mem ≈ **1.1×10⁻²²**.
  - An equal-mass SMBHB of 10⁹ + 10⁹ M☉ at 1 Gpc: h_mem ≈ **1.4×10⁻¹⁵**.
  - Θ = 0 must return exactly zero.
- **[mainstream] LIGO forecasts:** at design sensitivity an SNR 3 (5) memory detection needs about **35 (90) GW150914-like events** ([Lasky et al., PRL 117, 061102 (2016)](https://arxiv.org/abs/1605.01415)). Grant and Nichols forecast that LVK at O4 sensitivity for 1.5 yr plus O5 for 1 yr can detect displacement memory in the population. Cosmic Explorer could see it in individual loud events and detect **spin memory** in a population within 2 yr ([Grant & Nichols, PRD (v4 2024)](https://arxiv.org/abs/2210.16266)).
- **[mainstream] GWTC-3 stacking:** 88 BBH events give **ln BF = 0.01**, which the authors describe as favouring the no-memory hypothesis, so memory cannot yet be measured. They expect an O(2000)-event catalogue to be needed, and low-frequency non-Gaussian noise complicates beyond-GR amplitude tests ([Cheung, Lasky & Thrane 2024](https://arxiv.org/abs/2404.11919)).
- **[mainstream method, null result] GWTC-5.0 hierarchical inference:** the "memory enhancement factor", the coefficient of the supermomentum-flux term, is **0.26 (+4.09/−4.08)** at 68%, consistent with GR's value of 1. About 2000 detections are forecast to constrain it away from zero at 1σ. Multiplying Bayes factors naively across events can favour wrong conclusions ([Mitman, Isi & Farr, arXiv:2605.27500 v2, June 2026](https://arxiv.org/abs/2605.27500)).
- **[contested: preprint, waveform-systematics dependent] The full (ℓ=2, m=0) mode in GWTC-4.0** (IMRPhenomTHM_20, oscillatory plus memory parts): cumulative **log₁₀B = 1.38 ± 0.79** in favour of the (2,0) mode and stacked **SNR_memory = 0.89 (+0.29/−0.11)** after about 7.5 months of O4a. Decisive evidence is optimistically estimated to need **N_events = 166 (+82/−55)** ([Rosselló-Sastre et al., arXiv:2607.04909, July 2026](https://arxiv.org/abs/2607.04909)).
- **[mainstream] Pulsar timing.** The NANOGrav 15-yr memory search found at most **BF = 3.1** for a memory event, and only in the last three years of data, where too few TOAs follow a burst for it to be credible. Strain limits depend on sky position and epoch. They are not always better than the 12.5-yr limits, mostly because the common red noise rose, and **even the best sky regions do not reach below h₀ = 10⁻¹⁴** ([Agazie et al., ApJ 987 (2025), arXiv:2502.18599](https://arxiv.org/abs/2502.18599)). Signal model: δt(t) = B h_mem(t), with B = ½cos(2ψ)(1 − cos α) and h_mem(t) = h₀[(t−t₀)Θ(t−t₀) − (t−t_i)Θ(t−t_i)], where t_i = t₀ + (|p|/c)(1 + cos α). **Computed:** with B = 1/2 and h₀ = 10⁻¹⁴, the residual ramps by **0.79 µs** after 5 yr.
- **[mainstream theory, undetected] Spin memory** is a change in the magnetic-parity part of the time integral of the strain, sourced by the angular-momentum flux per solid angle. For quasi-circular inspirals it grows secularly through the inspiral, unlike displacement memory, which grows little before merger ([Nichols, PRD 95, 084048 (2017)](https://arxiv.org/abs/1702.03300)).

### Inferences
- The two 2026 estimates of the catalogue needed (about 166 versus about 2000 events) measure different things. One tests the presence of the whole (2,0) harmonic, which includes an oscillatory ringdown part. The other tests the memory coefficient alone. Report both, and never present the first as a memory detection.
- For a Rust `memory` module, implement Favata Eq. 17–18 (0PN) and Eq. 21 (integral over |ḣ₂₂|²) and use the PTA ramp model as a second backend. Unit tests: zero face-on, the 17/96 edge-on factor, and the 0.79 µs ramp.
- Memory from GW150914 (10⁻²²) is a few per cent of the oscillatory strain (10⁻²¹), and its power sits below 10 Hz. That explains why single-event detection is hopeless in LVK but natural for LISA and PTAs.

### Gaps
- Favata's 2009 ApJL saturation value for the final memory (the NR-calibrated coefficient in front of ηM/R) was not extracted. The 1.1×10⁻²² above is the 0PN value cut at the ISCO, not a fit.
- The Hübner et al. GWTC-1/GWTC-2 memory searches and the EPTA/PPTA memory limits were not fetched.
- No measurement of spin memory exists, even as a limit.

## 3. Propagation: speed, dispersion, graviton mass, lensing, matter, Gertsenshtein

### Takeaway
Gravity travels at c to within −3×10⁻¹⁵ to +7×10⁻¹⁶ (GW170817). The graviton mass is below 1.92×10⁻²³ eV/c² from GWTC-4.0, and below 1.0×10⁻²⁴ from planetary ephemerides. No lensed GW is confirmed, though GW231123 is an outlier. Absorption by matter is negligible by more than 20 orders of magnitude: viscous damping lengths are about 10¹¹ Gpc even in water. Gertsenshtein conversion during propagation reaches only about 10⁻¹⁶ in galactic fields, and in the LIGO band it is forbidden in interstellar plasma.

### Cited Findings
- **[mainstream] Speed of gravity:** the observed delay between GW170817 and GRB 170817A is **(+1.74 ± 0.05) s**, which bounds **−3×10⁻¹⁵ ≤ (v_GW − v_EM)/v_EM ≤ +7×10⁻¹⁶**. The analysis conservatively uses D = 26 Mpc, the lower distance bound, and allows the GRB to have been emitted up to 10 s after merger. A Shapiro-delay test of the equivalence principle gives **−2.6×10⁻⁷ ≤ γ_GW − γ_EM ≤ 1.2×10⁻⁶**, using only the Milky Way potential beyond 100 kpc (2.5×10¹¹ M☉) ([LVC + Fermi + INTEGRAL, ApJL 848, L13 (2017)](https://arxiv.org/abs/1710.05834)). **Computed:** 1.74 s × c/26 Mpc = 6.50×10⁻¹⁶ and (10 − 1.74) s × c/26 Mpc = 3.09×10⁻¹⁵, which reproduce both published bounds.
- **[mainstream] What c_T = c ruled out:** with |c_g/c − 1| ≤ 5×10⁻¹⁶, a large class of scalar-tensor dark-energy models is strongly disfavoured, including the covariant Galileon, Horndeski theory and most beyond-Horndeski generalisations. The same applies to other theories with a varying c_g, among them Einstein-Aether, Hořava gravity, generalised Proca, TeVeS and MOND-like gravities. Three loopholes survive: the simplest Horndeski terms, a conformal transformation, or a background-independent compensation ([Ezquiaga & Zumalacárregui, PRL 119, 251304 (2017)](https://arxiv.org/abs/1710.05901)).
- **[mainstream] Graviton mass and dispersion.** The LVK modified dispersion relation is E² = p²c² + A_α p^α c^α with α ∈ {0, 0.5, 1, 1.5, 2.5, 3, 3.5, 4}, where α = 0 is a massive graviton with A₀ = m_g²c⁴.
  - **GWTC-4.0** (42 O4a plus 49 earlier events): **m_g ≤ 1.92×10⁻²³ eV/c²** at 90%. This is a 1.16× improvement on the GWTC-3.0 bound of **2.23×10⁻²³** as recomputed in 2025–26. The planetary-ephemeris bound is **1.01×10⁻²⁴ eV/c² at 99.7%** (Mariani et al. 2023), an order of magnitude stronger ([LVK, GWTC-4.0 TGR II, arXiv:2603.19020 v2, July 2026](https://arxiv.org/abs/2603.19020)).
  - **Conflict:** the GWTC-3 TGR paper, revised as v3 in November 2025, states **m_g ≤ 2.42×10⁻²³ eV/c²**. It calls this 1.31× better than the Solar-System bound of 3.16×10⁻²³ and 1.28× better than a corrected GWTC-2 bound of 3.09×10⁻²³ ([LVK, arXiv:2112.06861 v3](https://arxiv.org/abs/2112.06861)). The often-quoted 1.27×10⁻²³ is therefore superseded.
  - **Computed:** Compton wavelengths h/(m_g c) are 6.46×10¹⁶ m (1.92×10⁻²³ eV), 5.12×10¹⁶ m (2.42×10⁻²³ eV) and 1.23×10¹⁸ m (1.01×10⁻²⁴ eV). The extra delay is Δt ≈ (D/2c)(m_gc²/hf)², or **1.1×10⁻⁴ s** for 100 Hz over 1 Gpc in flat space at the GWTC-4 bound.
- **[mainstream] Lensing of GWs (wave optics)** ([Takahashi & Nakamura, ApJ 595, 1039 (2003)](https://arxiv.org/abs/astro-ph/0305055)). Wave optics replaces geometrical optics when λ_GW exceeds the lens's Schwarzschild radius, which matters for lenses below about 10⁸ M☉ for LISA sources.
  - Point-mass amplification factor: F(f) = exp[πw/4 + i(w/2)(ln(w/2) − 2φ_m(y))] Γ(1 − iw/2) ₁F₁(iw/2, 1; iwy²/2), with **w = 8πGM_Lz f/c³**, φ_m = (x_m − y)²/2 − ln x_m and x_m = (y + √(y²+4))/2.
  - Geometrical limit: F = |μ₊|^{1/2} − i|μ₋|^{1/2}e^{2πifΔt}, with μ± = ½ ± (y²+2)/(2y√(y²+4)) and Δt = (4GM_Lz/c³)[y√(y²+4)/2 + ln((√(y²+4)+y)/(√(y²+4)−y))]. The typical delay is 2×10³ s (M_Lz/10⁸ M☉).
  - **Computed:** w = 1.24 for 100 M☉ at 100 Hz (w = 1 at 80.8 M☉). |F(y=0)|² = πw/(1 − e^{−πw}) = **3.2835** at w = 1 (from |Γ(1−iw/2)|² = (πw/2)/sinh(πw/2)). At y = 1, μ± = 1.17082 and −0.17082, and Δt = 4099 s for 10⁸ M☉, against 4GM/c³ = 1970 s.
- **[mainstream; GW231123 contested] Lensing searches.** The GWTC-4.0 searches of O4a used a saddle-point phase shift, pairs of events with matching frequency evolution, and sub-threshold counterparts. They found no evidence of strong lensing and used the null result to constrain lensed-event rates and the high-redshift BBH merger rate. The point-mass distortion search found one outlier, **GW231123_135430**, whose waveform uncertainties complicate interpretation ([LVK, arXiv:2512.16347 v3, Feb 2026](https://arxiv.org/abs/2512.16347)). An independent O1–O4a reanalysis finds a lensing probability below 0.6% for every pair and a lensing fraction below 1.4% at 90% ([Barsode, Maity & Ajith, arXiv:2607.08466](https://arxiv.org/abs/2607.08466); search-result summary only).
- **[mainstream] Absorption by matter is negligible.** Hawking (1966) showed that a viscous fluid damps GWs at **γ = 16πGη** (c = 1), i.e. 16πGη/c² in SI. Baym, Patil and Pethick rederive this from the Boltzmann equation, identify collisional and Landau damping, and show that collisions suppress the anisotropic stress. Damping within one expansion time is ≲ 1/(ωτ_U), about 10⁻²¹ for ω ∼ 10³ s⁻¹ and τ_U ∼ 10¹⁸ s, so interstellar and intergalactic damping is negligible except possibly for primordial waves ([PRD 96, 084033 (2017)](https://arxiv.org/abs/1707.05192)). **Computed:** γ = 3.73×10⁻²⁹ s⁻¹ for water (η = 10⁻³ Pa s), an absorption length c/γ of **2.6×10¹¹ Gpc**. For air the length is 1.4×10¹³ Gpc.
- **[mainstream] Dispersion by dark matter:** cold dark matter produces a small frequency-dependent change in GW speed, too small to detect in the near future ([Flauger & Weinberg, PRD 97, 123506 (2018)](https://arxiv.org/abs/1801.00386)).
- **[mainstream] Neutrino damping of primordial waves:** free-streaming neutrinos reduce the squared amplitude by **35.6%** for modes entering during radiation domination, independent of cosmological parameters, and lower the CMB tensor correlations by 9.0–10.7% ([Weinberg, PRD 69, 023503 (2004)](https://arxiv.org/abs/astro-ph/0306304)). **Computed:** the amplitude factor is √0.644 = 0.802.
- **[mainstream] Resonant absorption cross-section of a bar:** ∫σ_n dν = (8/π)(M/n²)(G/c)(V/c)² sin⁴θ_i [1 + σ_P(1−2σ_P)(nπR/L)²/2]. Here V is the sound speed and σ_P the Poisson ratio. For the LSU bar (2300 kg, 3 m, R = 0.6 m, 5.4 km/s, σ_P = 0.345) this is **4×10⁻²⁵ m² (Hz), 4.5×10²⁴ times smaller than the physical cross-section** ([Aguiar 2011, Res. Astron. Astrophys. 11, 1, Eq. 4.1](https://arxiv.org/abs/1009.1138)). **Computed:** 4.32×10⁻²⁵ m² Hz. Taking 2RL as the physical area gives a ratio of 8.3×10²⁴, so the paper must use a different area.
- **[mainstream] Gertsenshtein conversion during propagation.**
  - Conversion probability: P(g→γ) ≃ **8.3×10⁻⁵⁰ (B/1 G)²(L_coh/1 cm)²**. In plasma the coherence length is L_coh ≃ 3×10¹⁴ (f/10¹⁰ Hz)(n_e/1 cm⁻³)⁻¹ cm.
  - Total over a structure: P ≃ 7.2×10⁻¹¹ (B/1 G)²(f/10¹⁰ Hz)(n_e/1 cm⁻³)⁻¹(L_Σ/1 Mpc).
  - CMB spectral distortions then bound HFGWs to Ω_GW ≲ 1 near 10¹⁰ Hz, about 25 orders of magnitude better than the direct limits available in 2009 ([Pshirkov & Baskaran, PRD 80, 042002 (2009)](https://arxiv.org/abs/0903.4160)).
  - Rothman and Boughn give the vacuum mixing length L = 2/(G^{1/2}B) in Gaussian units with c = 1, roughly a megaparsec for 1 G, and P = sin²(G^{1/2}BD/2). For a photon crossing a 5 µG galactic field to the Galactic Centre, P ∼ 10⁻¹⁵. Strong fields add vacuum-polarisation refraction, which quenches conversion when DB²ω ≥ 10⁴³ (Dyson; Zel'dovich 1974) ([Rothman & Boughn, Found. Phys. 36, 1801 (2006)](https://arxiv.org/abs/gr-qc/0601043)).
  - **Computed:** the SI form 4πGB²L²/(μ₀c⁴) gives 8.263×10⁻⁵⁰ for 1 G and 1 cm, matching Pshirkov's prefactor and the report's 8.263×10⁻³⁸ for 1 T and 1 m. The vacuum mixing length is **2.25 Mpc** at 1 G. P = 3.1×10⁻¹⁶ for 5 µG over 8 kpc, and P = 7.9×10⁻¹⁹ for 1 µG over 1 kpc with full coherence. In plasma at f = 10¹⁰ Hz and n_e = 10⁻⁶ cm⁻³, L_coh = 9.7×10⁻⁵ Mpc.

### Inferences
- In the LIGO band, Gertsenshtein conversion cannot produce a propagating photon in the interstellar medium. f ∼ 100 Hz lies below the ISM plasma frequency (f_p ≈ 8.98 kHz × √(n_e/cm⁻³) = 1.56 kHz for n_e = 0.03 cm⁻³, computed). The process matters only for f far above MHz.
- The viscous-damping formula assumes the hydrodynamic regime, ωτ ≪ 1. Plugging in a geophysical mantle viscosity (∼10²¹ Pa s, a property of Maxwell times of order 1000 yr) would wrongly predict γ = 3.7×10⁻⁵ s⁻¹, an absorption length of only 8×10¹² m, about 54 AU (computed). At Hz frequencies rock is elastic, and absorption goes through the resonant cross-section above. A Rust function should take the relaxation time and refuse ωτ > 1.
- The speed bound uses an asymmetric emission-delay prior. Code should store both bounds and the assumptions (D = 26 Mpc, 0–10 s intrinsic delay), not a single "5×10⁻¹⁶".

### Gaps
- The v1 (2021) GWTC-3 graviton-mass value could not be retrieved, because the fetch was truncated by the author list. That 1.27×10⁻²³ was the original value is from memory and unverified. The difference between 2.42 and 2.23×10⁻²³ in the two LVK papers is not explained in the text read.
- The Barsode et al. lensing numbers come from a search summary only. The LVK O4a lensing rate constraints were not extracted.
- Dyson's 1969 work on the Earth's seismic response to GWs was not fetched.

## 4. Polarisations beyond GR

### Takeaway
A general metric theory allows six polarisations; GR allows two, both tensor. Ground interferometers cannot separate the breathing from the longitudinal scalar mode, so at most five modes are testable, and mode separation needs more detectors than modes. Pure-scalar and pure-vector hypotheses are excluded decisively (log₁₀B = −14.7 and −5.3 over GWTC-4.0, and +23/+21 for tensor over scalar/vector with GW170817's known sky position). Tensor-plus-extra mixtures remain unconstrained, except through Occam penalties.

### Cited Findings
- **[mainstream]** GR allows only two of the six polarisation states of a generic metric theory: two tensor (+, ×), two vector and two scalar (breathing, longitudinal). For L-shaped interferometers the breathing and longitudinal responses are degenerate, so the maximum number of separately testable modes is five. Earlier analyses therefore tested only pure hypotheses. The GWTC-3 analysis used null-stream projection with one or two basis modes ([LVK GWTC-3 TGR §VII, arXiv:2112.06861 v3](https://arxiv.org/abs/2112.06861)).
- **[mainstream] Network requirement:** the number of detectors must exceed the number of columns of the beam-pattern matrix F = [f₊ f× f_b]. In O1 and O4a only two detectors operated, so the analysis was limited to a single effective mode with complex scaling factors relative to the basis mode ([LVK GWTC-4.0 TGR I §5, arXiv:2603.19019](https://arxiv.org/abs/2603.19019)).
- **[mainstream] GW170817, sky position fixed to NGC 4993:** log₁₀B of pure tensor over pure vector is **+20.81 ± 0.08** and over pure scalar **+23.09 ± 0.08**. Without the sky position, the evidence against scalar is +5.84 ± 0.09 and the vector test is inconclusive at +0.72 ± 0.09 ([LVC, PRL 123, 011102 (2019)](https://arxiv.org/abs/1811.00364)).
- **[mainstream] GWTC-4.0 polarisation Bayes factors** (log₁₀B^X_T; negative favours tensor), from [GWTC-4.0 TGR I, Table 9](https://arxiv.org/abs/2603.19019):

| Hypothesis | GWTC-3.0 | O4a | GWTC-4.0 |
|---|---|---|---|
| S | −5.76 ± 0.42 | −8.97 ± 0.41 | **−14.72 ± 0.59** |
| V | −2.39 ± 0.41 | −2.94 ± 0.40 | **−5.33 ± 0.58** |
| TS | −0.07 ± 0.41 | −0.13 ± 0.40 | −0.20 ± 0.57 |
| TV | −0.04 ± 0.41 | +0.14 ± 0.40 | +0.10 ± 0.57 |
| VS | −2.75 ± 0.41 | −2.46 ± 0.40 | −5.21 ± 0.58 |
| TVS | −0.19 ± 0.41 | −0.13 ± 0.40 | −0.31 ± 0.57 |

- **[mainstream]** GWTC-3 found that mixed hypotheses with a tensor component cannot be ruled out conclusively, but the population is consistent with pure tensor polarisation ([GWTC-3 TGR](https://arxiv.org/abs/2112.06861)).
- **[mainstream] Pulsar timing:** the NANOGrav 15-yr correlations follow the Hellings–Downs (tensor) pattern, favoured over an uncorrelated common process by Bayes factors of 200–1000 ([Agazie et al., ApJL 951, L8 (2023)](https://arxiv.org/abs/2306.16213)).

### Inferences
- Mixed hypotheses (TS, TV, TVS) sit at log₁₀B ≈ 0 with uncertainty ±0.6. The current network therefore cannot see a small scalar or vector admixture, and LIGO-India or a fully operating Virgo and KAGRA is required. A Rust polarisation module should expose the six antenna patterns and a null-stream projector, and its tests should check that the breathing and longitudinal patterns are degenerate for 90° arms.
- The S and V Bayes factors grow roughly linearly with the number of events (GWTC-3 to GWTC-4 roughly tripled S). This is the expected stacking behaviour for a consistently wrong hypothesis.

### Gaps
- The explicit antenna-pattern formulas for vector and scalar modes were not extracted from a fetched source; Nishizawa et al. 2009 or Isi et al. are the usual references.
- PTA tests of non-Einsteinian correlations (scalar-transverse and others) in the NANOGrav, EPTA and MPTA data were not fetched.
- The GW170814 three-detector Bayes factors were not re-fetched.

## 5. Gravitons: Dyson's question, Tobar et al. 2024, cross-sections

### Takeaway
Every detected GW contains about 10³⁶–10³⁸ gravitons per cubic wavelength, so interferometers see classical fields. Graviton-matter cross-sections are about ℓ_P² ≈ 2.6×10⁻⁶⁶ cm², roughly 4×10⁻⁴² of the Thomson cross-section. Dyson and Rothman–Boughn concluded that single-graviton detection from natural thermal sources is effectively impossible, with about 4 events over the Sun's life in an Earth-sized detector. Tobar et al. (2024) show that stimulated absorption of one quantum in a ground-state-cooled bar of about 15–250 kg during a GW170817-like event is within reach. Carney et al. show that such a click would not prove the field is quantised.

### Cited Findings
- **[mainstream] Dyson (Poincaré Prize Lecture, ICMP Aalborg, 6 August 2012; Int. J. Mod. Phys. A 28, 1330041, 2013):** he argued that LIGO-type interferometers cannot detect single gravitons and left the question open for two other detector types ([search summary of the ICMP abstract and IAS repository record](https://albert.ias.edu/20.500.12111/9703)). His example was thermal gravitons from the Sun: a detector as large as the Earth running for the Sun's 5×10⁹ yr lifetime would see such an event only about four times ([Quanta Magazine, 30 Oct 2024](https://www.quantamagazine.org/it-might-be-possible-to-detect-gravitons-after-all-20241030/)).
- **[mainstream] Rothman and Boughn's detection arithmetic** ([Found. Phys. 36, 1801 (2006)](https://arxiv.org/abs/gr-qc/0601043)):
  - Detection criterion: nσλ ≥ 1, which leads to (σ/4π m_p)(M_d L_γ τ_s)/(R² ε_γ) ≥ 1 (Eq. 2.4).
  - The detector mass is capped at about a Jupiter mass, M_d ∼ (α/α_g)^{3/2} m_p, before gravitational collapse.
  - Gravito-Compton scattering has σ ∼ G²m² ∼ 10⁻¹¹⁰ cm² (∼10⁻⁴⁵ ℓ_P²). Gravito-atomic processes have σ ∼ ħG ∼ 10⁻⁶⁶ cm².
  - Gravito-ionisation of hydrogen has σ = 3072πħG(ak)⁵/[5(1 + a²k²)⁵] (Eq. 5.5).
  - The best natural case, a Jupiter-mass detector in close orbit around a neutron star, gives about 10⁻² detections per year.
  - A LIGO-class wave with h ∼ 10⁻²⁰ at 10³ Hz carries about 10³⁸ gravitons per cubic wavelength.
  - They judge that Dyson's conjecture is very likely true once realistic detectors are considered.
- **Computed (the same arithmetic):**
  - ℓ_P² = Għ/c³ = **2.612×10⁻⁶⁶ cm²**, and (Gm_e/c²)² = 4.58×10⁻¹¹¹ cm².
  - (α/α_g)^{3/2} m_p = 2.30×10²⁷ kg, about 1.2 Jupiter masses.
  - The maximum ionisation cross-section, at ak = 1, is 60.3 ℓ_P² = 1.58×10⁻⁶⁴ cm².
  - **Graviton-to-photon cross-section ratio:** ℓ_P²/σ_Thomson = 2.6×10⁻⁶⁶/6.65×10⁻²⁵ = **3.9×10⁻⁴²**.
- **[mainstream] Carney, Domcke and Rodd ([PRD 109, 044009 (2024)](https://arxiv.org/abs/2308.12988)):**
  - Occupation per de Broglie volume: nλ_dB³ = πh²M_Pl²/(2f²) ≃ **2×10³⁵ (h/10⁻²²)²(1 kHz/f)²**, with M_Pl the reduced Planck mass and ħ = c = 1.
  - CAST's stochastic sensitivity h_c ≃ 5×10⁻²⁸ at 10¹⁸ Hz corresponds to ρ ≃ 2×10¹⁸ keV/cm³ and n ≃ 4×10¹⁷ gravitons/cm³, but nλ_dB³ ∼ 10⁻⁵. A CAST click from a transient with h ≃ 7×10⁻²⁸(T_m/T_gw)^{1/4} would therefore be a single-graviton detection.
  - Because a classical wave acting on a quantised detector produces the same discrete clicks, proving quantisation needs sub-Poissonian counting statistics, which they consider likely impossible in practice.
  - **Computed:** nλ³ = 2.15×10³⁵ (h = 10⁻²², 1 kHz) and 1.05×10⁻⁵ (h = 7×10⁻²⁸, 10¹⁸ Hz). For h = 10⁻³⁰, nλ³ = 1 at f = 4.6×10¹² Hz.
- **[contested in interpretation; formulas mainstream] Tobar, Manikandan, Beitel and Pikovski, Nat. Commun. 15, 7229 (2024)** ([arXiv:2308.15440](https://arxiv.org/abs/2308.15440)):
  - **Spontaneous emission rate:** Γ_spon = 8GML²ω_l⁴/(l⁴π⁴c⁵) = 8πGρv_s⁴R²/(Lc⁵) (Eq. 4). It is ≈ 10⁻⁴⁰ Hz for an atom's 3d→1s transition (Weinberg) and ≈ 10⁻³³ Hz for a niobium bar (ρ = 8570 kg/m³, v_s = 5 km/s, L = 1 m, R = 0.5 m). **Computed: 9.28×10⁻³⁴ Hz.**
  - **Stimulated rate:** Γ_stim = ML²ω_l²h²/(4l⁴π⁵ħ) = v_s²Mh²/(4l²π³ħ) (Eq. 5). For a 1800 kg aluminium bar with h = 5×10⁻²², Γ_stim ≈ 1 Hz. **Computed: 1.003 Hz** at v_s = 5.4 km/s.
  - **Optimal mass:** M = π²ħω³/[v_s²χ(h,ω,t)²] with χ = ∫₀ᵗ ḧ(s)e^{iωs}ds (Eq. 6). For a slow chirp, M ≈ (24π²/5)(ħ/(h₀²v_s²))(GM_c/2c³)^{5/3}ω^{8/3}, which gives ∼15 kg of beryllium at 100 Hz with h₀ ∼ 2×10⁻²². **Computed: 12.8 kg** (Be, v_s = 13 km/s, M_c = 1.19 M☉), and **220 kg** for aluminium at 150 Hz, against the paper's ∼250 kg.
  - **Table I:**

| Source | f | Material | Mass |
|---|---|---|---|
| GW170817 | 100 Hz | Be | ∼15 kg |
| GW170817 | 150 Hz | Al | ∼250 kg |
| GW170608 | 175 Hz | Nb | ∼9 t |
| GW150914 | 200 Hz | CuAl6% | ∼6 t |
| BH superradiance (A0620−00) | 33 kHz | sapphire | ∼100 kg |

    These assume T ≈ 1 mK and Q ≈ 10¹⁰ for transients.
  - The single-phonon probability peaks at P_max = 1/e ≈ 0.37 for |β| = 1.
  - **Graviton count:** N = h₀²c⁵/(32πGħν²) (Eq. A22, with ν the angular frequency) gives N ≈ 4×10³⁶ for h₀ = 10⁻²¹ at 150 Hz. **Computed: 3.85×10³⁶.** Detected GWs contain ≳ 10³⁶ gravitons.
  - Quanta reports that roughly one in three suitable waves would excite a 15 kg beryllium bar by one quantum, confirmed by coincidence with LIGO ([Quanta 2024](https://www.quantamagazine.org/it-might-be-possible-to-detect-gravitons-after-all-20241030/)).

### Inferences
- An independent cross-check computed here supports Tobar et al.'s regime. The classical energy a 2300 kg LSU-type bar absorbs from GW170817 at 900 Hz is 1.2×10⁻³¹ J, or **0.19 ħω** (Section 8), so a resonant bar absorbs about one quantum from a nearby BNS. The physics is classical absorption cross-section times fluence; the "graviton" reading needs the quantum-jump readout and the interpretation caveat.
- "Graviton–photon cross-section" has two distinct meanings that code must not conflate. Absorption and scattering by matter go as ℓ_P² (gravito-atomic) or G²m² (gravito-Compton). Graviton–photon conversion in a field B goes as GB²L²/c⁴ (Section 3).

### Gaps
- Dyson's 2013 paper itself could not be read: the IAS repository returned 403 and the ICMP server refused connections. His own numbers (graviton luminosity of the Sun, the cross-section he used) are therefore taken from Quanta's paraphrase. "About 4 events" is a secondary-source figure.

## 6. Primordial gravitational waves

### Takeaway
B-mode bounds give r < 0.036 (BICEP/Keck 2021), r < 0.032 (Planck PR4 + BK18) and r < 0.034 (all CMB, end of 2025). This caps the inflation energy scale at V^{1/4} ≲ 1.4×10¹⁶ GeV, H_inf ≲ 4.7×10¹³ GeV and Ω_GW h² ≲ 1.3×10⁻¹⁶ today. The nanohertz PTA background (A ≈ 2.4×10⁻¹⁵) is about 3–4σ evidence of HD correlations. It fits supermassive black-hole binaries with somewhat stretched parameters. Cosmological fits score Bayes factors of 10–100 over SMBHBs, but these are explicitly not evidence for new physics.

### Cited Findings
- **[mainstream] BICEP/Keck 2021:** r₀.₀₅ < **0.036** at 95% with σ(r) = 0.009, from all BICEP2, Keck Array and BICEP3 data through 2018 combined with WMAP and Planck ([PRL 127, 151301](https://arxiv.org/abs/2110.00483)).
- **[mainstream] Planck PR4 + BK18 (Tristram et al., PRD 105, 083524 (2022)):** r < **0.032** at 95%. The history of the bound runs PR2+BK 0.09, BK15 0.07, PR3+BK15 0.065, PR4 alone 0.056 and BK18 0.036 ([arXiv:2112.07961](https://arxiv.org/abs/2112.07961)).
- **[mainstream] End-of-2025 combination (Planck + SPT + ACT + BICEP/Keck):** **r < 0.034** at 95% and n_s = 0.9682 ± 0.0032. Adding DESI BAO shifts n_s to 0.9728 ± 0.0029, which favours monomial potentials with N* ∼ 50 over Starobinsky R² and Higgs inflation; the n_s shift itself is marginal ([Balkenhol et al., OJAp, arXiv:2512.10613 v2, June 2026](https://arxiv.org/abs/2512.10613)).
- **[mainstream] Energy-scale formulas** ([Baumann TASI lectures, Eqs. 216–218](https://arxiv.org/abs/0907.5424)): Δ_t² = 2H²/(π²M_pl²) with the reduced Planck mass, r = Δ_t²/Δ_s², and **V^{1/4} ∼ (r/0.01)^{1/4} × 10¹⁶ GeV**. The exact slow-roll forms are V = (3π²/2)A_s r M_pl⁴ and H = πM_pl√(A_s r/2). Planck 2018 gives **ln(10¹⁰A_s) = 3.044 ± 0.014** ([Planck 2018 VI](https://arxiv.org/abs/1807.06209)).
- **Computed** (A_s = 2.0989×10⁻⁹, M_pl = 2.4353×10¹⁸ GeV):

| r | V^{1/4} (GeV) | H_inf (GeV) |
|---|---|---|
| 0.036 | **1.408×10¹⁶** | 4.70×10¹³ |
| 0.034 | 1.388×10¹⁶ | 4.57×10¹³ |
| 0.032 | 1.368×10¹⁶ | 4.43×10¹³ |
| 0.01 | 1.022×10¹⁶ | 2.48×10¹³ |
| 0.001 | 5.75×10¹⁵ | 7.84×10¹² |

- **[mainstream] Present-day energy density** ([Caprini & Figueroa, CQG 35, 163001, Eq. 215](https://arxiv.org/abs/1801.04268)): Ω_GW(k) = (3/128)Ω_rad P_h(k)[(1/2)(k_eq/k)² + 16/9], which tends to Ω_rad P_h/24 for k ≫ k_eq. **Computed:** with P_h = rA_s (n_T ≈ 0) and Ω_rad h² = 4.18×10⁻⁵ (Ω_γh² = 2.47×10⁻⁵ from the same review, plus N_eff = 3.046 neutrinos), Ω_GW h² ≤ **1.3×10⁻¹⁶** for r < 0.036. This neglects the g* drop and the 35.6% neutrino damping, both of which lower it.
- **[mainstream] NANOGrav 15-yr** ([ApJL 951, L8](https://arxiv.org/abs/2306.16213)), from 67 pulsars:
  - A GWB power law is favoured over pulsar-intrinsic noise by BF > 10¹⁴, and over an uncorrelated common process by BF = 200–1000.
  - The false-alarm probability for the observed BF is p = 10⁻³ (≈3σ); a frequentist correlation statistic gives p = 5×10⁻⁵ to 1.9×10⁻⁴ (≈3.5–4σ).
  - Amplitude: **A = 2.4 (+0.7/−0.6)×10⁻¹⁵** at 1 yr⁻¹ for an f^{−2/3} h_c spectrum.
- **[contested] Interpretation as SMBHBs:** astrophysical models reproduce the amplitude and shape. The amplitude, however, needs many parameters at the edges of expected values, or a few notably different from standard expectations ([NANOGrav, ApJL 952, L37](https://arxiv.org/abs/2306.16220)).
- **[contested] Interpretation as new physics:** inflation, scalar-induced GWs, first-order phase transitions, cosmic strings and domain walls can all reproduce the signal, except stable field-theory cosmic strings. Many fit better than SMBHB-only, with BF = 10–100, but these depend on SMBHB modelling and "should not be regarded as evidence for new physics" ([NANOGrav, ApJL 951, L11](https://arxiv.org/abs/2306.16219)).
- **[contested] MeerKAT PTA, 4.5 yr** (83 pulsars, about 250 000 TOAs): one noise model gives HD significance of 3–3.4σ, while slightly different assumptions give an insignificant result. The cause is that a few nearby high-precision pulsars dominate. Amplitudes: h_c,yr = 7.5×10⁻¹⁵ at α = −0.26, and 4.8×10⁻¹⁵ at α = −2/3 ([Miles et al., MNRAS 536, 1489 (2025), arXiv:2412.01153](https://arxiv.org/abs/2412.01153)). In the June 2023 coordinated release, each IPTA member reported 2–4σ evidence ([MPTA maps paper](https://arxiv.org/abs/2412.01214)).

### Inferences
- Vacuum inflationary waves (Ω h² ≲ 10⁻¹⁶) are about 7 orders of magnitude below the PTA signal (Ω h² ≈ 3.6×10⁻⁹ at 1/yr, from the existing report). An inflationary PTA explanation therefore needs a strongly blue tensor tilt or a non-standard post-inflation history, which is why such models appear only in the new-physics fit set.
- A Rust `inflation` helper needs only A_s, r and the reduced Planck mass. Test against the table above to 10⁻⁴.

### Gaps
- No IPTA DR3 combined result (2025–26) was found. The SPT-3G-only r < 0.25 and the forecast σ(r) ∼ 0.004 (BK + SPT delensing) come from search snippets of talks and preprints that were not read in full.

## 7. Standard sirens: H₀ from GW170817 and dark sirens

### Takeaway
GW170817 alone gives H₀ = 70.0 (+12.0/−8.0) km s⁻¹ Mpc⁻¹. The GWTC-4.0 dark sirens (GLADE+), combined with GW170817, give 75.4 (+12.8/−9.1), with the dark sirens alone giving 79.6 (+20.8/−16.3). GW cosmology is still at the ±12–15% level and cannot yet arbitrate the 67-versus-73 Hubble tension.

### Cited Findings
- **[mainstream] Bright siren GW170817** ([LVC et al., Nature 551, 85 (2017)](https://arxiv.org/abs/1710.05835)): **H₀ = 70.0 (+12.0/−8.0) km s⁻¹ Mpc⁻¹** (68.3%), from d = 43.8 (+2.9/−6.9) Mpc and the NGC 4993 Hubble-flow velocity v_H = 3017 ± 166 km s⁻¹. The distance–inclination degeneracy dominates the error. **Computed:** v_H/d = 3017/43.8 = 68.9, which is the naive point estimate; the published posterior peaks higher because of the inclination prior.
- **[mainstream] GWTC-4.0 cosmology** ([LVK, ApJL 1007, L17, 10 Aug 2026, arXiv:2509.04348](https://arxiv.org/abs/2509.04348)), using candidates through O4a with the GLADE+ catalogue:
  - Combined with GW170817: **H₀ = 75.4 (+12.8/−9.1) km s⁻¹ Mpc⁻¹** at 68.3%, (+24.5/−13.3) at 90%.
  - Dark sirens alone: 79.6 (+20.8/−16.3).
  - GW170817 reanalysed with the same priors: 77.8 (+25.7/−11.6).
  - Modified GW propagation: **Ξ₀ = 1.4 (+1.0/−0.4)** (GR is Ξ₀ = 1).
  - Spectral-siren bounds tightened by about 45% relative to GWTC-3.0.
- **[contested: independent analyses, read from search summaries only]** Using GWTC-4.0 with a mass feature near 63 M☉: spectral sirens 78.8 (+19.0/−15.3) and dark sirens 82.5 (+16.8/−14.3) at 68% ([arXiv:2601.03257](https://arxiv.org/abs/2601.03257v1)). Cross-correlation with DES Y6: 71.7 (+22.7/−18.3), or 73.1 (+11.9/−8.5) with GW170817 ([Glasgow eprint 389795](https://eprints.gla.ac.uk/389795/)).

### Inferences
- The GWTC-4 central value (75) sits nearer SH0ES (∼73) than Planck (∼67), but the 68% interval covers both. No tension claim is justified.
- Formula for a calculator: v_H = H₀ d_L at low z, and d_L = (1+z)(c/H₀)∫dz/E(z) in general. Modified propagation scales d_L^GW/d_L^EM through Ξ₀; the full parameterisation was not extracted.

### Gaps
- The number of bright sirens needed for a 2% H₀ (the "∼50 events" forecast) was not re-sourced. The GWTC-3 dark-siren value was not re-fetched, and the exact event count in GWTC-4 cosmology (141 events with FAR < 0.25 yr⁻¹) comes from a search summary.

## 8. GWs and gravity manipulation: what a wave does to matter, and the hard limits

### Takeaway
The strongest waves ever to reach Earth carry about 0.02 W m⁻². That is comparable to moonlight, yet they stretch a metre by about 10⁻²¹ m and accelerate it by about 10⁻¹⁵ m s⁻². A resonant bar absorbs about 10⁻³¹ J (≈0.2 quanta) from a nearby neutron-star merger, a momentum kick of about 10⁻⁴³ m/s. Even a Galactic black-hole merger, delivering about 3×10⁴ times the solar constant at peak, would give only 6×10⁻¹¹ m s⁻² per metre. GWs carry large energy fluxes that matter cannot absorb, because the coupling is set by G/c³ through a cross-section near 10⁻²⁵ m² Hz. The luminosity scale c⁵/G maps directly onto a strain ceiling, h ≲ 2√(2℘)·c/(ωr).

### Cited Findings
- **[mainstream] Tidal action (geodesic deviation):** ẍ_j = ½ḧ^TT_jk x_k in the proper detector frame ([Flanagan & Hughes 2005](https://arxiv.org/abs/gr-qc/0501041)). For h = h₀cos ωt the tidal acceleration amplitude is a = ½ω²h₀L. **Computed:** h = 10⁻²¹ at 250 Hz gives a = 1.23×10⁻¹⁵ m s⁻² across 1 m, and a 4 km LIGO arm changes by ΔL = hL/2 = 2×10⁻¹⁸ m.
- **[mainstream] Energy flux at Earth from the strongest signals:** GW150914 delivered ∼10⁵ erg s⁻¹ m⁻² (= 0.01 W m⁻²) at peak ([PRL 116, 241102](https://arxiv.org/abs/1602.03840)), and typical bursts carry ∼1.5 mW m⁻², comparable to full moonlight ([Flanagan & Hughes](https://arxiv.org/abs/gr-qc/0501041)). GW250114 has a similar strain, ∼10⁻²¹ ([arXiv:2509.08054](https://arxiv.org/abs/2509.08054)). **Computed** (GW150914):
  - Peak flux 0.0174 W m⁻² from L/(4πD²), or 0.0198 W m⁻² from c³ḣ²/(16πG).
  - Total fluence 2.7×10⁻⁴ J m⁻², an isotropic average.
- **[mainstream] Absorbed energy and momentum.** **Computed** from the bar cross-section (Aguiar Eq. 4.1, Section 3) and the Newtonian inspiral spectrum dE/df = (π^{2/3}/3)G^{2/3}M_c^{5/3}f^{−1/3}, which follows from the quadrupole flux and the ḟ formula. For GW170817 (M_c = 1.188 M☉, 40 Mpc), the spectral fluence at 900 Hz is 2.67×10⁻⁷ J m⁻² Hz⁻¹. An optimally oriented 2300 kg aluminium bar would absorb:
  - **E_abs = 1.15×10⁻³¹ J**,
  - that is **0.19 ħω** and 8×10⁻⁸ k_B × (0.1 K),
  - momentum E/c = 3.8×10⁻⁴⁰ kg m s⁻¹,
  - and a centre-of-mass Δv of 1.7×10⁻⁴³ m/s.
- **[mainstream] The whole Earth as a GW detector:** ten years of superconducting-gravimeter data from the Global Geodynamics Project, read through the Earth's normal modes at 0.3–5 mHz, constrain the background to Ω_GW between 0.035 and 0.15 ([Coughlin & Harms, PRD 90, 042005 (2014)](https://arxiv.org/abs/1406.1147)). This is directly relevant to the project's gravimetry tooling.
- **[mainstream] Momentum carried away by the source** (recoil) is the only macroscopic GW "propulsion" ever observed: about 1500 km/s for GW200129, i.e. ≈0.5% of M_f c (Section 1, [Varma et al. 2022](https://arxiv.org/abs/2201.01302)). For any collimated GW beam, the thrust per watt is P/c = 3.34×10⁻⁹ N W⁻¹, the same as a photon rocket (computed).
- **[mainstream as a unit; the bound is contested] Maximum luminosity and strain.** c⁵/G = 3.63×10⁵² W. Observed BBH peaks reach about 10⁻³ of it, and numerical critical collapse reaches 0.2 ([Cardoso et al. 2018](https://arxiv.org/abs/1803.03271)). No theorem bounds L in GR ([Jowsey & Visser 2021](https://arxiv.org/abs/2105.06650)).
- **Computed (derived here from F = c³h₀²ω²/32πG = ℘c⁵/(G·4πr²)):** **h₀ = 2√2·√℘·c/(ωr)**.
  - With ℘ = 10⁻³ at 250 Hz and r = 410 Mpc, this gives **h₀ = 1.35×10⁻²¹**, which reproduces GW150914's observed 1.0×10⁻²¹ and is a good closed-form sanity test.
  - At r = λ with ℘ = 1 it gives h₀ = 0.45: the strain becomes order unity only within about a wavelength of a Planck-luminosity source.
- **Computed: the loudest plausible wave.** A GW150914 analogue in the Milky Way at 8 kpc would give h ≈ 5.1×10⁻¹⁷ and a peak flux of 4.6×10⁷ W m⁻² (3.4×10⁴ solar constants). The tidal acceleration across 1 m would still be only 6.3×10⁻¹¹ m s⁻².

### Inferences
- "Gravity manipulation" with GWs fails on absorption, not on available energy. A strong wave carries kW m⁻²-class power within a few kpc of a merger, but the absorption cross-section of a 2.3 t bar is about 10⁻²⁵ m² Hz, roughly 10²⁵ times smaller than its geometric area. A GW170817-like merger at 8 kpc would deposit only about 3×10⁻²⁴ J (≈5×10⁶ quanta) in such a bar. This scales the 40 Mpc result by (40 Mpc/8 kpc)² = 2.5×10⁷.
- In a crate, keep three distinct quantities and test each one: flux (W m⁻²), tidal acceleration (m s⁻² m⁻¹) and absorbed energy (J, via Σ × spectral fluence). Many fringe claims of GW "beams" or "propulsion" confuse the first with the third.
- Earth normal modes are the only geophysical GW detector with a published Ω bound. Superconducting-gravimeter time series, already within the project's domain, can reproduce the Coughlin–Harms analysis.

### Gaps
- The Newtonian spectrum dE/df was used from the standard derivation, not from a fetched source; Sathyaprakash & Schutz 2009 in Living Reviews is the usual reference. The absorbed-energy figure uses an isotropic-average fluence and the optimal bar orientation, so it is an order-of-magnitude estimate.
- No published "strongest GW ever at Earth" figure in W m⁻² was found beyond the GW150914 PRL estimate. GW250114's flux was not published and would be similar.
- Dyson's 1969 calculation of the Earth's seismic response to GWs was not fetched.
