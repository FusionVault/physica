# Gravity beams: the gravitational fields of light and matter beams, beam-like exact solutions, and generating gravitational waves in the lab

How the labels are used throughout:
- **[A] ESTABLISHED**: standard GR, peer-reviewed, consistent with experiment. Formulas we can compute in Rust.
- **[B] SPECULATIVE-SERIOUS**: mainstream authors and journals, but untested or far beyond current sensitivity.
- **[C] FRINGE/REFUTED**: failed replication, refuted, or not accepted by the community.

Every number marked "computed" was checked with a short Python calculation (G = 6.674e-11, c = 2.998e8, SI units). The checked values are given in the Inferences sections.

---

## 1. The gravitational field of light beams (Tolman–Ehrenfest–Podolsky, Bonnor, pp-waves and gyratons, Rätzel 2016, spin and OAM light)

### Takeaway
GR gives light beams a well-defined and fully established gravitational field [A]:
- Co-propagating beams do not attract each other. Counter-propagating beams do attract.
- A test mass at rest next to a beam of power P feels an acceleration a = 4GP/(c³ρ). For a petawatt pulse at 1 mm this is about 1e-17 m/s², and it lasts only while the pulse passes.
- Spin and orbital-angular-momentum (OAM) corrections are smaller still, around 1e-29 m/s².

No laser-gravity effect is measurable with any foreseeable technology. The most promising "beam" source is the LHC proton beam. Its signal-to-noise ratio (S/N) reaches about 16 in one week only after a high-luminosity upgrade, an optimised tungsten pendulum, and neglecting all technical noise.

### Cited Findings
**Classical foundations [A]**
- Tolman, Ehrenfest and Podolsky (1931), "On the Gravitational Field Produced by Light", Phys. Rev. 37, 602–615, DOI 10.1103/PhysRev.37.602. In the weak-field approximation, test rays moving parallel to a light pencil or pulse "do so with uniform unit velocity", so co-propagating light does not attract. Counter-propagating pulses do attract. — [INSPIRE record](https://inspirehep.net/literature/46619); [Caltech record](https://authors.library.caltech.edu/records/a22aa-d4x62)
- Bonnor (1969), "The gravitational field of light", Commun. Math. Phys. 13, 163–174, DOI 10.1007/BF01645484:
  - It gives an exact Einstein solution for a steady light beam (now called the "Bonnor beam").
  - A second solution has two parallel co-moving beams, which do not interact.
  - A uniform beam is gravitationally stable.
  - The solutions are plane-fronted gravitational waves (pp-waves). — [search summary of the Springer record](https://doi.org/10.1007/bf01645484); [Wikipedia: Bonnor beam](https://en.wikipedia.org/wiki/Bonnor_beam)
- Faraoni and Dumse (1999), Gen. Rel. Grav. 31, 91–105 (arXiv:gr-qc/9811052), explain why parallel pp-waves superpose linearly while antiparallel ones scatter and interact nonlinearly. For parallel beams the gravitoelectric force is "exactly cancelled by a gravitomagnetic force". For equal energy densities, light-to-light attraction is twice matter-to-light attraction and four times matter-to-matter attraction. — [arXiv:gr-qc/9811052](https://arxiv.org/abs/gr-qc/9811052)
- Aichelburg and Sexl (1971), "On the gravitational field of a massless particle", Gen. Rel. Grav. 2, 303–312 (bibcode 1971GReGr...2..303A):
  - The field of a single photon-like particle is obtained by ultra-boosting Schwarzschild (v→1, m→0 at fixed energy).
  - The result is an impulsive pp-wave.
  - In vacuum the pp-wave field equation reduces to a 2-D Laplace equation for the profile function. — [Wikipedia: Aichelburg–Sexl ultraboost](https://en.wikipedia.org/wiki/Aichelburg%E2%80%93Sexl_ultraboost); [arXiv:1606.04259 review](https://arxiv.org/pdf/1606.04259); [arXiv:1204.0744](https://arxiv.org/pdf/1204.0744)
- Gyratons: Frolov and Fursaev (2005), "Gravitational field of a spinning radiation beam-pulse in higher dimensions", Phys. Rev. D 71, 104034 (hep-th/0504027).
  - The linearised field of a spinning radiation beam can be obtained by boosting the Lense–Thirring metric to near light speed.
  - They give an exact D-dimensional solution containing floor(D/2) arbitrary functions of retarded time u.
  - In 4-D the metric has a du·dx cross term, which signals dragging into rotation by the beam's spin. — [arXiv:hep-th/0504027](https://arxiv.org/abs/hep-th/0504027); [search summary](https://arxiv.org/pdf/hep-th/0504027)

**Rätzel, Wilkens and Menzel (2016), NJP 18, 023009, arXiv:1511.01023 [A: linearised GR]**
- Metric: h₀₀ = h_zz = −h₀z = h^p, with h^p = κ ln(ζ(b)/ζ(a)) and κ = 4GP/c⁵ (their Eqs. 7, 14, 30). — [arXiv HTML](https://arxiv.org/html/1511.01023v1)
- Acceleration of a test particle at rest is (c²/2)∂ₓh^p, so |ẍ| ≈ 4GP/(c³x) near the pulse (Eq. 31). A counter-propagating massless particle feels 2c²∂ₓh^p, four times as much (Eq. 33). Co-propagating massless particles feel nothing. — [arXiv HTML](https://arxiv.org/html/1511.01023v1); [NJP abstract](https://content.cld.iop.org/journals/1367-2630/18/2/023009/revision1/abstract-transcript-NJPaa1155va.txt)
- All physical effects are confined to spherical shells expanding at c from the emission and absorption events:
  - The emission shell attracts test masses at rest towards the beam axis.
  - The absorption shell repels them.
  - Linear polarisation modulates the field at the norm of E. Circular polarisation gives no modulation. — [arXiv abstract](https://arxiv.org/abs/1511.01023)
- Numbers:
  - The paper gives |ẍ| ≈ 10⁻¹⁸ m/s² for P ≈ 10¹⁵ W at x = 2.5 mm.
  - At ρ = 1 cm, the curvature of a 10¹⁵ W beam is comparable to that of an astrophysical GW with ω = 10³ and h₊ = 10⁻²². — [arXiv HTML](https://arxiv.org/html/1511.01023v1)

**Spengler, Rätzel and Braun (2022), "Perspectives of measuring gravitational effects of laser light and particle beams", NJP 24 (5), arXiv:2104.09209 [A for the physics, B for the measurement]**
- Formulas:
  - Curvature next to a cw beam of power P: R₀ₓ₀ₓ ≃ −4GP/(c⁵ρ²) (Eq. 1).
  - A non-relativistic sensor feels the equivalent Newtonian potential Φ = (4GP/c³) ln ρ while the pulse is passing, and zero at other times (Eq. A29).
  - Quote: "there is no gravitational wave generated as this type of source is not quadrupolar in nature. Rather one can detect the gravitational near-field." — [arXiv:2104.09209 PDF](https://arxiv.org/pdf/2104.09209v5)
- Sources (Table I):
  - Laser pulses in a cavity: 3·10¹⁴ W pulse power, average 2·10¹⁰ W.
  - cw laser plus cavity: 2·10¹¹ W.
  - LHC: 10¹⁴ W pulse power, average 3.8·10¹² W, waist 16 µm. That is 2808 bunches of 1.15·10¹¹ protons at 6.5 TeV, each bunch about 10⁵ J. — [arXiv:2104.09209 PDF](https://arxiv.org/pdf/2104.09209v5)
- Detector amplitudes (Table II, 5 mK):
  - Rod: lasers 2·10⁻²⁵ m, LHC 8·10⁻²⁴ m, against a sensitivity of 1·10⁻¹⁷ m/√Hz.
  - Pendulum: LHC 1·10⁻²⁴ m, against a sensitivity of 2·10⁻¹⁴ m/√Hz.
  - Liquid-helium detector: the LHC force is 3.5 orders of magnitude below the detection limit after one week.
  - Monolithic pendulum: the LHC signal is 5 orders of magnitude below the sensitivity. — [arXiv:2104.09209 PDF](https://arxiv.org/pdf/2104.09209v5)
- Best case: an optimised oscillator gives S/N ≈ 0.6 with the LHC in one week. With the high-luminosity LHC (×10) and tungsten test masses (×2.9) it gives S/N ≈ 16, "neglecting all sources of technical noise". The authors call S/N > 1 "a necessary condition" for what "would still make for a very difficult experiment". — [arXiv:2104.09209 PDF](https://arxiv.org/pdf/2104.09209v5); [edoc HU Berlin](https://edoc.hu-berlin.de/items/3203340f-6c8c-4595-b59c-c5e275a21a40/full)

**Spin and OAM of light [A for the theory, effects unmeasurable]**
- Schneiter, Rätzel and Braun (2018), "The gravitational field of a laser beam beyond the short wavelength approximation", DOI 10.1088/1361-6382/aadc81 (Class. Quantum Grav.), arXiv:1804.08706.
  - Treating the beam as an exact Maxwell solution with diffraction, they find frame dragging from the beam's spin angular momentum, with a sign set by helicity.
  - Example: P₀ ∼ 10¹⁵ W, w₀ ∼ 1 mm, divergence ∼ 10⁻³, test particle at x = w₀ moving at 10 m/s gives a transverse acceleration of ±10⁻²⁹ m/s².
  - Overlapping co-propagating rays now feel a nonzero relative acceleration, ∼ −10⁻³¹ m/s². The authors say this is "too small to be detected with current technology". — [arXiv abstract](https://arxiv.org/abs/1804.08706); [arXiv HTML](https://arxiv.org/html/1804.08706)
- The same group showed that the polarisation of a test beam rotates in a laser's gravitational field (gravitational Faraday effect and optical activity). The effect appears only beyond geometric optics. — [arXiv:1812.04505](https://arxiv.org/html/1812.04505v2)
- Strohaber:
  - "Frame dragging with optical vortices" (arXiv:1112.3414) gives OAM-induced frame dragging from Laguerre–Gaussian (LG) beams.
  - arXiv:1807.00933 (2018) gives weak-field GR of LG beams with spin, OAM and spin-orbit frame-dragging terms, and proposes a gravitational Aharonov–Bohm analogue.
  - No lab-scale numbers appear in the abstract. — [arXiv:1112.3414](https://arxiv.org/pdf/1112.3414); [arXiv:1807.00933](https://arxiv.org/abs/1807.00933)
- Atonga, Martineau, Aboushelbaya, Barrau, … Norreys (2023, v2 2024), "Gravitational waves from high-power twisted light", arXiv:2309.04191. This is the main 2023–2024 OAM paper.
  - It uses Bessel pulses with ℓ = 0 and ±1 at 351 nm, with NIF-ARC (0.4–1.7 kJ, 1.3–38 ps) and SEL (100 PW: 1.5 kJ in 15 fs) parameters.
  - Result: h ∼ 1.90×10⁻³⁶ (P/1 PW)(r/1.36 µm)⁻¹ at ω_g = 2ω ≈ 1.07×10¹⁶ s⁻¹.
  - GW intensity: 1.44×10⁻⁵ W/m² near the source, 1.01×10⁻¹⁹ W/m² at 10 m.
  - The radiated power scales as ω_g² L² P₀² G/c⁷ (their Eq. 86).
  - The authors judge the chances of detecting light's gravity as a Newtonian or tidal force to be "extremely low".
  - The effective mass of a 1 MJ laser pulse is about 10⁻¹¹ kg.
  - No journal reference was found. — [arXiv:2309.04191 HTML](https://arxiv.org/html/2309.04191v2)
- The "gravitational spin Hall effect" covers the opposite direction: gravity acting on light's spin, not the field light produces. Oancea et al., Phys. Rev. D (2020), arXiv:2003.04553, derive polarisation-dependent ray equations from covariant WKB (Wentzel–Kramers–Brillouin) theory. — [arXiv:2003.04553](https://ar5iv.labs.arxiv.org/html/2003.04553); [overview arXiv:1904.09963](https://www.arxiv.org/abs/1904.09963)
- Other 2025 papers on structured light in gravitational fields, all about gravity acting on light:
  - Gravitational Faraday effect and OAM, arXiv:2501.15846.
  - GWs changing photon OAM states, arXiv:2504.16452.
  - Orbital Hall effect of twisted particles, arXiv:2510.03450. — [2501.15846](https://arxiv.org/pdf/2501.15846); [2504.16452](https://arxiv.org/html/2504.16452v2); [2510.03450](https://arxiv.org/pdf/2510.03450)

### Inferences
- **Formulas for Rust [A]**, all in SI units:
  - Beam acceleration: `a(P, ρ) = 4 G P / (c³ ρ)`.
  - Tidal curvature: `R0x0x(P, ρ) = −4 G P / (c⁵ ρ²)`.
  - Equivalent potential: `Φ = (4GP/c³) ln ρ`.
  - Metric: `h00 = hzz = −h0z = (8GP/c⁵) ln ρ`, up to an additive constant. This follows from Φ = (c²/2)h₀₀ and matches Rätzel's κ = 4GP/c⁵.
  - Velocity kick from one pulse of energy E passing at distance ρ (a·τ with E = Pτ): `Δv = 4 G E / (c³ ρ)`. This is our derivation.
  - The factor 4G (instead of 2G for a Newtonian line mass λ = P/c³) is the "light gravitates twice as strongly" factor from (ρ + p) with p = ρc² for radiation.
  - Aichelburg–Sexl in the same normalisation: ds² = −du dv + dx² + dy² − (8Gp/c³) δ(u) ln(ρ/ρ₀) du². Check the sign and factor conventions against the original before coding.
- **Computed numbers:**
  - a(1 PW, 1 mm) = 9.9e-18 m/s². a(1 PW, 2.5 mm) = 4.0e-18 m/s², consistent in order of magnitude with Rätzel's "~1e-18".
  - a(100 kW cw, 1 mm) = 9.9e-28 m/s².
  - a(cavity 2e11 W, 100 µm) = 2.0e-20 m/s².
  - a(LHC 3.8e12 W, 100 µm) = 3.8e-19 m/s².
  - Equivalent point mass for 1 PW at 2.5 mm: 3.7e-13 kg. One extraction of Rätzel's text gave "1e-13 g"; that may be a units slip by the extractor, so check the paper itself.
  - Curvature of 1 PW at 1 cm: 1.1e-33 m⁻². A GW with h = 1e-22 at 1 kHz has hω²/(2c²) = 2.2e-32 m⁻², so Rätzel's comparison holds to an order of magnitude.
  - LHC stored energy: 3.36e8 J, mass-equivalent 3.7e-9 kg, average power E·c/L_ring = 3.78e12 W. This matches Spengler's Table I.
  - Velocity kick from a 100 J (1 PW × 100 fs) pulse at 1 mm: Δv ≈ 1e-30 m/s. Hopeless.
- **Can any of it be measured?** No laser effect can be measured with any foreseeable technology. The effects are about 10 orders of magnitude below the best detectors.
  - The LHC near-field effect is the only "beam gravity" proposal with S/N > 1 even on paper, and only after neglecting technical noise.
  - A lab measurement of dynamic, time-varying Newtonian gravity at 42 Hz already exists (Brack et al. 2022, Section 2), so the near-field detection method is real. The LHC signal is still about 4–5 orders of magnitude weaker per unit test mass than that demonstration.
- **For the user's "gravity beam" idea:** a beam of light or relativistic matter does carry a gravitational field. That field is a tube-shaped near field (a logarithmic 2-D potential), not a radiated beam of gravity. Co-moving objects feel nothing, and nothing is projected ahead of the pulse.

### Gaps
- I did not read the Atonga et al. conclusion section (truncated fetch), so their final detectability statement and journal status are unverified.
- I found no 2023–2025 paper that computes a measurable gravitational field from OAM light. The 2025 OAM papers found are about gravity acting on light.
- I did not see Bonnor's 1970 spinning-null-fluid paper or the Frolov–Israel–Zelnikov follow-up gyraton papers, so their exact citations are not verified.
- I did not see the explicit gyraton metric functions in the abstract.

---

## 2. Generating gravitational waves in the laboratory (quadrupole formula, Weber-era rotating bars, gravitational Hertz experiments, 2015–2026 proposals)

### Takeaway
The quadrupole formula makes lab GW generation hopeless for mechanical sources:
- A 490-tonne, 20 m steel beam spun to breakup radiates about 2×10⁻²⁹ W [A]. Its wavelength is about 34,000 km, so any lab detector sits deep in the near field.
- Only GHz-range electromagnetic schemes, mainly the Gertsenshtein effect (conversion between EM waves and GWs in a strong magnetic field), allow a true "Hertz" geometry with generator and receiver in each other's far field.
- The most optimistic peer-reviewed estimate is Kolosnitsyn and Rudenko (2015): h ~ 10⁻³²–10⁻³⁰ and about 14 detected photons/s, assuming ideal conditions and ignoring all technical noise [B].

No lab-generated GW has ever been detected. What has been measured is dynamic near-field Newtonian gravity: Hirakawa et al. 1980, and Brack et al. 2022 at 42 Hz.

### Cited Findings
**Quadrupole formula and classic numbers [A]**
- The GW luminosity from the quadrupole formula has prefactor G/(5c⁵) times the squared third time-derivative of the traceless mass quadrupole. A standard exercise applies it to a spinning steel rod. — [Nordita/Stockholm 2024 lecture exercise](https://indico.fysik.su.se/event/8276/contributions/13166/attachments/5761/7341/stockholm24ex.pdf)
- Misner, Thorne and Wheeler's discouraging example is the Einstein-style steel beam: about 20 m long, 1 m radius, 4.9×10⁵ kg, spun end-over-end at the speed its tensile strength allows. — [Chiao et al. arXiv:1712.08680](https://arxiv.org/pdf/1712.08680)
  - Another paper gives the same beam a radiated power "of order 10⁻²³ erg/s" alongside an internal power flow of about 2×10⁸ erg/s. — [arXiv:0707.2318](https://arxiv.org/pdf/0707.2318). This conflicts with the commonly quoted 2.2×10⁻²² erg/s; see Inferences.
- Weber's first lab "Hertz experiment" analysis found a generated power of up to 10⁻¹³ erg/s (10⁻²⁰ W), "much less than the receiver sensitivity". The general conclusion: under lab "slow motion and weak field" conditions, an effective generator is "extremely difficult", and the Hertz experiment is more Gedanken than practical. — [Kolosnitsyn & Rudenko arXiv:1504.06548](https://arxiv.org/pdf/1504.06548)

**Gravitational Hertz experiment via the Gertsenshtein effect [B]**
- Kolosnitsyn and Rudenko (2015), "Gravitational Hertz experiment with electromagnetic radiation in a strong magnetic field", arXiv:1504.06548 (Physica Scripta). — [arXiv abstract](https://arxiv.org/abs/1504.06548); [PDF](https://arxiv.org/pdf/1504.06548)
  - Field equation: □h = −16πG τ/c⁴, sourced by the Maxwell stress tensor.
  - Generated amplitude, Gaussian units (Eq. 10): h_max = 2G H₀ H̃₀ L/(c⁴k). A Fabry–Pérot (FP) cavity enhances it by the finesse F (Eq. 12).
  - Parameters: H₀ = 10 T, L = 100 m, λ = 1 cm, 1 − R = 10⁻⁸, 10 kW pump in 1 cm². These give h ~ 10⁻³². Superconducting resonators with F ~ 10¹⁰ give h ~ 10⁻³⁰.
  - Inverse effect: the receiver photon flux is N_S ≈ 10⁻³ s⁻¹ at h = 10⁻³². With F = 10⁹ (h ≈ 2×10⁻³¹) it rises to about 14 s⁻¹.
  - Against a 4 K black-body background of about 6×10⁶ photons in the cavity, they estimate 4×10⁵ s (about 4 days) of integration.
  - Their positive conclusion holds "only at the level of 'potential opportunity'": all of the GW hits the receiver, there are no losses, and instrumental noise is ignored.
  - The GW is emitted along the EM beam axis, in both directions, as two anti-phase components. This phase-matched geometry makes the emission beam-like (see Section 3).
- Gertsenshtein (1962), "Wave resonance of light and gravitational waves", is the pioneering paper. A 2001 analysis found the efficiency of a cavity-resonator generator "many orders below the observability threshold". — [search summary of a historical review](https://producciocientifica.uv.es/documentos/648b52a6b367423619da04c7?lang=fr)

**Superconductor-based generation proposals [B, with a null experiment]**
- Chiao, Fitelson and Speliotopoulos (2003) tested a YBCO "transducer" of microwaves to GWs and back. They found no signal and set an upper limit of 1.6×10⁻⁵ on the conversion efficiency at 77 K. The preprint was never published in a journal. Kiefer and Weber called the underlying coupling rule "questionable". — [Gallerati, Modanese & Ummarino 2022 review, arXiv:2206.07574](https://arxiv.org/html/2206.07574v3)
- Chiao et al. (2017), "Dynamical Casimir effect and the possibility of laser-like generation of gravitational radiation" (arXiv:1712.08680), claim that coherent GW generation is possible "above a certain threshold" if superconducting RF cavities have high Q for both EM and gravitational microwaves. This relies on an unproven "relative gravitational permeativity" and a Meissner-like GW reflection. No journal reference was found. — [arXiv:1712.08680](https://arxiv.org/abs/1712.08680)
- A 2007 proposal (arXiv:0710.1378) uses two superconducting sphere-ring systems in separate dewars as a microwave Hertz experiment. — [arXiv:0710.1378](https://arxiv.org/pdf/0710.1378)

**High-power-laser GW generation [B]**
- Kadlecová, Klimo, Weber and Korn (2017), "Gravitational wave generation by interaction of high power lasers with matter using shock waves", Eur. Phys. J. D 71, 89, DOI 10.1140/epjd/e2017-70586-y (arXiv:1602.08904). They derive analytic metric perturbations and radiated power for laser-driven shock waves. Part II covers ablation and piston models (arXiv:1608.04925). — [arXiv:1602.08904](https://arxiv.org/abs/1602.08904v4); [Part II](https://arxiv.org/pdf/1608.04925)
- Atonga et al. (2023/24), twisted light: h ∼ 10⁻³⁶ at µm distances for 1 PW (Section 1). The radiated power loss is "too tiny to have a substantial impact on the dynamics of the system". — [arXiv:2309.04191 HTML](https://arxiv.org/html/2309.04191v2)

**Near-field lab gravity at Hz to kHz [A, measured]: this is not GW radiation**
- Nature 283, 184 (1980) tested the inverse-square law to R = 4.2 m with a "low frequency dynamical quadrupole–quadrupole interaction" between two masses. The result agreed with the inverse-square law. This is the Hirakawa group's rotating-bar method; a follow-up is Ogawa, Tsubono and Hirakawa, PRD 1982. — [Nature 283, 184](https://www.nature.com/articles/283184a0.pdf); [PRD 26 issue](https://journals.aps.org/prd/issues/26/4)
- Brack et al. (2022), "Dynamic measurement of gravitational coupling between resonating beams in the hertz regime", Nature Physics 18, 952–957, DOI 10.1038/s41567-022-01642-8 (ETH Zürich).
  - A transmitter beam vibrating near 42 Hz drives a high-Q detector beam through gravity alone. The displacement is about 1 pm, and the result agrees with Newtonian modelling to about 2%.
  - A rotating-bar follow-up appeared in Communications Physics (2023), DOI 10.1038/s42005-023-01389-5. — [APS Physics news](https://physics.aps.org/articles/v15/121); [APS 2023 abstract](https://meetings-archive.aps.org/apr/2023/b14/3)

**Detector-side context: ultra-high-frequency GWs (UHF-GW) [A/B]**
- Aggarwal et al. (2021), "Challenges and opportunities of gravitational-wave searches at MHz to GHz frequencies", Living Rev. Relativ. 24, 4, DOI 10.1007/s41114-021-00032-5. This is the mainstream UHF-GW review. It notes there are no known astrophysical sources in the band, which makes it a window for physics beyond the Standard Model. — [UWA repository](https://research-repository.uwa.edu.au/en/publications/challenges-and-opportunities-of-gravitational-wave-searches-at-mh/)
- A 2025 follow-up review, "Challenges and Opportunities of Gravitational Wave Searches above 10 kHz" (arXiv:2501.11723), was cited in a reference list; I did not read it. — [search result list](https://arxiv.org/pdf/2506.20198)
- Kolosnitsyn and Rudenko quote the relic GW background above about 10¹¹ Hz as h ~ 10⁻³² Hz⁻¹/², with BBN limits of h ~ 10⁻³⁴–10⁻³⁷ Hz⁻¹/² at 10¹²–10¹⁵ Hz. That is far below any detector. — [arXiv:1504.06548](https://arxiv.org/pdf/1504.06548)

### Inferences
- **Quadrupole formulas for Rust [A]:**
  - General: `P = (G/5c⁵) ⟨d³Q_ij/dt³ · d³Q_ij/dt³⟩`.
  - Rigid body spinning at Ω about a principal axis: `P = (32/5)(G/c⁵)(I₁−I₂)² Ω⁶`.
  - Thin rod spun about its centre (I = ML²/12): `P = (2/45) G M² L⁴ Ω⁶ / c⁵`.
  - Strain amplitude at distance d (the pulsar-style formula): `h₀ = 4 G (I₁−I₂) Ω² / (c⁴ d)`.
  - GW energy flux: `F = c³ ω² (h₊² + h×²) / (32π G)`.
  - Natural luminosity scale: c⁵/G = 3.6e52 W.
- **Computed numbers:**
  - MTW beam (4.9e5 kg, 20 m, Ω = 28 rad/s): P = 2.3e-29 W (2.3e-22 erg/s), matching the textbook value usually quoted from MTW §36. The arXiv:0707.2318 figure of "~1e-23 erg/s" appears to assume a different Ω or prefactor.
  - The same beam: GW wavelength πc/Ω = 3.4e7 m. Strain at one wavelength: h ≈ 1e-41.
  - 1 t, 2 m rod at 50 Hz: P = 1.9e-32 W. 10 kg, 1 m rod at 100 Hz: P = 7.5e-36 W.
  - Kolosnitsyn h = 1e-32 at 30 GHz over 1 cm²: F ≈ 1.4e-8 W/m², about 1.4e-12 W of GW power.
  - Single-pass Gertsenshtein conversion probability (SI): `P_conv = 4πG B² L² / (μ₀ c⁴)`. This equals G B²L²/c⁴ in Gaussian units and is consistent with Kolosnitsyn's Eq. 10 combined with the flux formula. B = 10 T, L = 100 m gives 8.3e-32; L = 1 m gives 8.3e-36.
- **Pump leakage in Kolosnitsyn's scheme:** the generator pump (10 kW at 30 GHz, about 5e26 photons/s) is at the same frequency as the expected signal (about 14 photons/s). A real experiment would need roughly 10⁻²⁶–10⁻²⁷ (about 260–270 dB) EM isolation between generator and receiver. Their "potential" analysis does not address this, and it is probably the real obstacle. This is our inference.
- **Lab detectors sit in the near field:** a lab detector is in the GW far zone only if r ≫ λ_GW = c/f. At 1 kHz, λ = 300 km. So any table-top "GW" experiment below about 1 GHz measures near-field (Newtonian or tidal) coupling, which dominates the radiative part by powers of (λ/r). A genuine lab Hertz experiment needs f ≳ 1–10 GHz, which leaves only EM or microwave schemes.
- **Overall [B]:** these proposals are physically sound in principle. In practice they sit 10 or more orders of magnitude (mechanical) or several orders plus unsolved isolation problems (Gertsenshtein) from a demonstration. A "lab GW Hertz experiment" remains unrealised.

### Gaps
- I could not open the Kadlecová et al. numerical h and power values (only the abstract was available).
- The MTW textbook page was not checked directly. The 2.2e-29 W value is reproduced by our calculation, not quoted from the book.
- I found no 2024–2026 paper claiming a feasible, noise-budgeted lab generator-plus-detector. The 2024–2025 papers found are detector-only: superconducting levitated detector arXiv:2408.01583; EM cavities as bars, Berlin et al. 2023. — [arXiv:2408.01583](https://arxiv.org/pdf/2408.01583)

---

## 3. GW "beams", focusing and lensing, directional emission, and "gravitational wave lasers"

### Takeaway
- GW lensing by masses is established theory. LIGO–Virgo–KAGRA searched all O3 binary-black-hole events and found no confident lensing detection [A].
- Directional ("collimated") GW emission requires a coherent source much larger than λ_GW: phase-matched EM–GW conversion along a beam, or a phased array of quadrupoles [B].
- A "GW laser" (gaser) has a serious version and a fringe version:
  - Serious [B]: stimulated GW emission in astrophysical boson clouds (Liu 2024), and Schützhold's 2025 PRL on graviton stimulated emission or absorption by light, which is a detection scheme.
  - Fringe [C]: Fontana's high-Tc "gaser" and Chiao's dynamical-Casimir generator.
- Nothing has been demonstrated.

### Cited Findings
- LVK (2023), "Search for gravitational-lensing signatures in the full third observing run of the LIGO-Virgo network", arXiv:2304.08393. It covers all O3 binary-black-hole events; none of the searches gave significant evidence of lensing, and the non-detection was used to constrain lensing rates. A follow-up of O3 candidates likewise reports "no confident detections". — [arXiv:2304.08393](https://arxiv.org/abs/2304.08393v1); [UWA listing](https://research-repository.uwa.edu.au/en/publications/search-for-lensing-signatures-in-the-gravitational-wave-observati/); [follow-up arXiv:2306.03827](https://ar5iv.labs.arxiv.org/html/2306.03827)
- Coherent multi-quadrupole emitters (atoms, molecules, superconductor vortices at 10²²–10²⁴ cm⁻³) were proposed for 10¹⁰–10¹⁵ Hz. The tiny G/c⁵ factor would be offset by the number of emitters, by ω⁶, and by "a very sharp beaming of the output signal". This is the HFGW programme associated with R. M. L. Baker's 2003 and 2007 conferences. — [Kolosnitsyn & Rudenko arXiv:1504.06548](https://arxiv.org/pdf/1504.06548)
- In Gertsenshtein conversion the GW is generated synchronously with, and collinear to, the EM wave. It leaves the FP cavity along ±x, which makes the source intrinsically beam-like. — [arXiv:1504.06548](https://arxiv.org/pdf/1504.06548)
- Atonga et al. claim that the frequency, polarisation and emission direction of GWs from twisted light are controllable through pulse parameters and optics. Their single-source angular pattern is g_D(θ) = (1/16)[3 − 4cos2θ + cos4θ]. — [arXiv:2309.04191](https://arxiv.org/html/2309.04191v2)
- The "gaser":
  - It is hypothetical. Gravitons, as spin-2, have quadrupole rather than dipole transitions.
  - Fontana (cond-mat/0208276) argues stimulated emission is the only viable mechanism for quantum-system GW production and proposes high-Tc superconductors as the active medium [C/B].
  - Liu (arXiv:2401.16096) argues stimulated GW emission can occur in "gravitational atoms" (black holes surrounded by ultralight-boson clouds) [B, astrophysical].
  - A July 2026 preprint (arXiv:2607.02068) proposes "graviton lasing" with a BEC or superfluid helium in an optical cavity [B, unreviewed]. — [Wikipedia: Gravity laser](https://en.wikipedia.org/wiki/Gravity_laser); [cond-mat/0208276](https://arxiv.org/pdf/cond-mat/0208276); [arXiv:2401.16096](https://arxiv.org/pdf/2401.16096); [arXiv:2607.02068](https://arxiv.org/pdf/2607.02068)
- Schützhold (2025), "Stimulated Emission or Absorption of Gravitons by Light", Phys. Rev. Lett. 135, 171501, DOI 10.1103/xd97-c6d7 [B].
  - Setup: a two-stage Mach–Zehnder interferometer in which a passing astrophysical GW shifts photon frequencies, followed by about 10⁶ km of phase accumulation.
  - At LIGO-like strain this gives a phase of about 10⁻⁷. That would show energy exchange in units larger than ħω, a signature of stimulated graviton processes.
  - It would not be conclusive proof of gravitons, because a semiclassical interaction gives the same phase. — [APS Physics Viewpoint v18/175](https://physics.aps.org/articles/v18/175); [HZDR press release](https://www.hzdr.de/db/Cms?pOid=76137)

### Inferences
- **Beaming limit:** a single compact quadrupole radiates a broad pattern (sin⁴-type or quadrupole lobes) and cannot be collimated. Collimation requires a coherent source of size D ≫ λ_GW, with angular width θ ~ λ_GW/D (ordinary diffraction). At 1 kHz that means D ≫ 300 km. Only GHz-and-above schemes can be "beams" at lab scale, which is why HFGW proposals cluster there.
- **Lensing in Rust:** wave-optics GW lensing matters when λ_GW is comparable to the lens Schwarzschild radius, 2GM_lens/c². In the geometric-optics limit the lensing is the same as for light, and these formulas are standard.
- **What "GW lasers" would actually take:** population inversion on quadrupole transitions with lab-relevant gain. Graviton emission rates per atom are so small (G/c⁵ suppression) that no proposed lab medium has a demonstrated or credible threshold. The astrophysical boson-cloud version is serious because black-hole superradiance supplies enormous occupation numbers.

### Gaps
- I did not retrieve wave-optics lensing primary sources (for example Takahashi & Nakamura 2003) or the exact O3 lensing-rate constraints.
- I did not read the 2026 graviton-lasing preprint beyond its search summary.

---

## 4. Fringe claims: Baker / Li "HFGW" generators and detectors, Podkletnov and Modanese gravity shielding and "impulse gravity generator", Tajmar's gravitomagnetism, Ning Li

### Takeaway
None of these claims has been independently reproduced.
- **Baker / Li:** the US government's JASON panel (2008) judged the HFGW "threat" not credible for communication, imaging or propulsion.
- **Podkletnov:** his weight-loss claims failed replication at NASA Marshall, at Hathaway et al. (to 0.001%), and elsewhere. His "impulse gravity generator" and "64c" speed claims come only from the inventor and his co-author.
- **Tajmar:** he withdrew his own early anomalies after better vibration and thermal isolation, attributing them to acoustic noise and helium effects. An independent ring-laser test (Graham et al. 2008) was null.
- **Ning Li:** her superconductor theory was criticised as wrong "by many orders of magnitude", and NASA tests were null.

Scientific consensus: these are [C].

### Cited Findings
**R. M. L. Baker, F.-Y. Li et al., HFGW and the "Li–Baker detector" [C for the generation and propulsion claims; the detector's sensitivity claims are disputed]**
- Li et al. (2008), Eur. Phys. J. C 56, 407–423. This is the Li–Baker coupled EM-system HFGW detector paper. — [Li–Baker 2008 PDF hosted by Baker](https://www.drrobertbaker.com/docs/Li-Baker%20et%20al.%20HFGW%20Detector%202008.pdf); [search summary](https://katalog.lib.cas.cz/KNAV/EdsRecord/edb,28723840)
- The JASON report "High Frequency Gravitational Waves", JSR-08-506 (October 2008), was commissioned by the Office of the Director of National Intelligence (ODNI). It concluded "No foreign threat in HFGW is credible", covering communication, object detection or imaging, vehicle propulsion, or any other practical use. It also said "There is no substitute for seeking expert scientific and technical opinion in such matters." — [FAS Secrecy News, Dec 15 2008](https://fas.org/publication/jason_study); [report PDF link](https://irp.fas.org/agency/dod/jason/gravwaves.pdf)
- Baker's rebuttal, from an interested party (GravWave LLC):
  - JASON based its analysis mainly on the Gertsenshtein effect and found GravWave's analysis of the Li–Baker detector wrong "by an order of magnitude".
  - Baker says the detector is "not based on the Gertsenshtein effect" but on "a different effect found by Dr. Fangyu Li".
  - He disputes JASON's dismissal of HFGW propulsion.
  - Baker reports about 35 papers and six HFGW patents. — [Baker Q&A on JASON](https://www.drrobertbaker.com/docs/Q%20&%20A%20JASON.pdf); [KurzweilAI](https://www.kurzweilai.net/researchers-challenge-gravitational-waves-report-for-director-of-national-intelligence)
- Baker's patent filing alleges a missing 1/μ₀ factor in JASON's Eq. 4, the EM-to-GW power conversion. That is an advocate's assertion, and a patent is not peer review. — [USPTO 10322827](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/10322827)

**Podkletnov (gravity shielding) and Podkletnov & Modanese ("impulse gravity generator") [C]**
- Podkletnov and Nieminen (1992, Physica C) reported:
  - A 145 mm YBCO disk, levitated and spun to several thousand rpm, with test masses above it losing up to 0.05% of their weight.
  - Critics attributed this to systematic errors.
  - A 1997 version with a 27 cm toroidal disk claimed 0.3–0.5%. It was accepted by J. Phys. D, then withdrawn, and posted to arXiv instead. — [Gallerati, Modanese & Ummarino, Front. Phys. 10, 941858 (2022), arXiv:2206.07574](https://arxiv.org/html/2206.07574v3). Note that Modanese is Podkletnov's co-author, so this review is not neutral, yet it still reports the null replications.
- Replications:
  - Woods et al. (2001): disks up to 3 inches, no effect.
  - Hathaway et al. (Physica C): disk rotated mechanically only, at lower speed and temperature. No weight change or gravity-like force down to 0.001%.
  - NASA Marshall and the University of Alabama in Huntsville: only static tests on a 4-inch single-layer disk were published, with "limited results".
  - Junker (from 2004, half scale): inconclusive. — [arXiv:2206.07574](https://arxiv.org/html/2206.07574v3); [search summary](https://arxiv.org/pdf/2206.07574)
- NASA Marshall (Noever & Koczor 1998, RF-illuminated levitated YBCO, gravimeter):
  - Any apparent change was 3–5×10⁻⁵ cm/s², "significantly lower" than the claimed 2.1% modification.
  - The peak coincided with the superconductor's AC-resistance peak, which suggests magnetic-field leakage into the instrument.
  - A related 1999 paper covers fabricating 30 cm annular disks. — [NTRS 19990064220](https://ntrs.nasa.gov/api/citations/19990064220); [NTRS 19990068007](https://ntrs.nasa.gov/api/citations/19990068007)
- Impulse gravity generator:
  - Podkletnov & Modanese (2001), arXiv:physics/0108005, used high-voltage discharges through a YBCO electrode. A 2002 "second experiment" claimed a beam carrying about 1 mJ with momentum far above E/c.
  - A 2012 book chapter by the same two authors claims a 63 ± 1 ns delay over 1211 m between two piezo sensors with Rb clocks, that is "64c". Some versions of the abstract print "631 ns", an inconsistency.
  - The authors themselves say this "does not necessarily imply superluminal information transmission".
  - There is no independent replication. Tajmar & Lorincz saw accelerometer effects, which a 2016 follow-up attributed to acoustic noise once the EM shielding was improved. — [arXiv:physics/0108005](https://arxiv.org/abs/physics/0108005v2); [Bentham chapter](https://www.benthamdirect.com/content/books/9781608053995.chapter-8); [Eurekaselect](https://www.eurekaselect.com/public/chapter/3862); [arXiv:2206.07574](https://arxiv.org/html/2206.07574v3)

**Tajmar and de Matos: gravitomagnetic London moment [C, refuted largely by Tajmar's own later work]**
- The claim: a rotating superconductor produces a gravitomagnetic field B_g ≈ 1.84×10⁻⁴ ω, invoked to explain Tate's Cooper-pair mass excess (m′/2mₑ = 1.000084(21)). The 2022 review notes this is enormous compared with Earth's frame dragging (about 10⁻¹⁴ rad/s). — [arXiv:2206.07574](https://arxiv.org/html/2206.07574v3); [Tajmar & de Matos gr-qc/0603033](https://arxiv.org/abs/gr-qc/0603033); [gr-qc/0607086](https://arxiv.org/pdf/gr-qc/0607086)
- Claimed signals:
  - 2007: accelerometer signals of about 10⁻⁴–10⁻³ m/s² near a rotating 69 mm Nb ring.
  - 2008: gyroscope signals from Al, Nb and YBCO, only for clockwise rotation.
  - The group reported "all systematic effects… at least 3 orders of magnitude below" the signal. — [arXiv:2206.07574](https://arxiv.org/html/2206.07574v3); [gr-qc/0610015](https://arxiv.org/pdf/gr-qc/0610015); [0707.3806](https://arxiv.org/abs/0707.3806v2)
- Refutation:
  - Graham et al. (2008, Canterbury) used a lead disk and a large ring-laser gyroscope. The result was within 1σ of null, against about 11σ for the Tajmar group.
  - 2009: liquid helium co-rotating with the sample affected the results.
  - Tajmar et al. (2011) reinterpreted the anomalies as acoustic noise, after new equipment measured signals about two orders of magnitude lower. — [arXiv:2206.07574](https://arxiv.org/html/2206.07574v3)
- Tajmar, "Evaluation of Enhanced Frame-Dragging in the Vicinity of a Rotating Niobium Superconductor, Liquid Helium and a Helium Superfluid", Supercond. Sci. Technol. 24, 125011 (2011). His TU Dresden group states that after proper vibration and thermal isolation, the gyroscope signal fell below the noise level. — [TU Dresden Breakthrough Propulsion Physics page](https://tu-dresden.de/ing/maschinenwesen/ilr/rfs/forschung/forschungsfelder/raumfahrtantriebe-und-neue-konzepte/breakthrough-propulsion-physics)
- A dissenting framing: Hauser and Dröscher (a heterodox theory) attribute the differing results to thermal environments. — [search summary](https://arxiv.org/pdf/0908.4495)

**Ning Li ("AC gravity") [C]**
- Li & Torr (1991), Phys. Rev. D 43, 457: a theoretical paper on gravitomagnetic fields in pure superconductors. — [DOI 10.1103/PhysRevD.43.457](https://dx.doi.org/10.1103/PhysRevD.43.457)
- E. G. Harris (1999, Found. Phys. Lett.) commented on the Torr & Li work: because of unrealistic assumptions, the calculated fields are "too large by many orders of magnitude". — [DOI 10.1023/A:1021621425670](https://dx.doi.org/10.1023%2FA%3A1021621425670)
- Li's claim of an "11 kilowatt gravity-like beam" ("AC gravity") is known only from low-quality secondary sources. — [mirror page (low quality)](https://ufologie.patrickgross.org/htm/sciegra.htm)
- A Tajmar & de Matos review states that these experiments "could not be reproduced by others up to now", and that the theories "were either shown to be wrong or are often based on difficult to prove assumptions". — [search summary, gr-qc/0607086 and related](https://arxiv.org/pdf/gr-qc/0607086)

### Inferences
- **Size of the Tajmar claim:** GR frame dragging inside a rotating shell of mass M and radius R is Ω_drag ≈ (4GM/(3c²R)) ω. For M = 1 kg and R = 0.1 m this is about 1e-26 ω (computed). The claimed B_g ≈ 1.8e-4 ω is therefore about 22 orders of magnitude above GR. Extraordinary claims of this size need extraordinary replication, and the replication went the other way.
- **Baker's propulsion and communication claims:** even the most optimistic peer-reviewed lab generator (Kolosnitsyn: about 1e-12 W of GWs from a 10 kW, 10 T, 100 m superconducting setup) gives a thrust of P/c ≈ 3e-21 N. This is an inference, but it matches JASON's conclusion that the claims are not credible.
- **The "64c" claim** conflicts with the multi-messenger bound from GW170817 and GRB 170817A, which limits the difference between GW and light speeds to about 1e-15. That bound is well known (LVK/Fermi/INTEGRAL, ApJL 848, L13, 2017) but was not re-fetched here.
- **Summary for the user:** the [C] claims share a pattern. Effects at the 1e-4–1e-2 level, with no GR mechanism, vanished as controls improved or were never reproduced. The real 2020s frontier (Section 1: Spengler; Section 2: Brack) operates at 1e-19–1e-25 m-scale signals and is honest about it.

### Gaps
- I could not extract the JASON report text: the PDF returned empty, possibly because it is a scan. Its specific technical calculations (detector noise, generator efficiency) are known here only through the FAS summary and Baker's rebuttal.
- I found no neutral peer-reviewed paper that quantifies the error in the Li–Baker detector's sensitivity.
- I could not identify the NASA NTRS record for the Ning Li co-authored null result (reported in a search snippet as "<2 parts in 10⁸ of g") or the venue of the Hathaway paper, so both are unverified.
- I found no primary source for Ning Li's later career or the AC Gravity LLC claims.

---

## 5. Gravitomagnetism and frame dragging as real "beam-like" effects (Gravity Probe B, LAGEOS/LARES, LARES-2)

### Takeaway
Gravitomagnetism is real, measured physics [A]:
- **Gravity Probe B:** measured frame dragging of −37.2 ± 7.2 mas/yr against a GR prediction of −39.2 (about 19% precision), and the geodetic effect at −6601.8 ± 18.3 against −6606.1 mas/yr.
- **LARES-2 with LAGEOS** (Ciufolini et al., Nature 2026): measured Earth's Lense–Thirring node drag to about one part in a thousand. A critic (Iorio) disputes the error budget.

These are the same h₀ᵢ terms that make a light beam's gravity cancel for co-moving test light.

### Cited Findings
- Everitt et al. (2011), Phys. Rev. Lett. 106, 221101 (arXiv:1105.3456):
  - Geodetic drift: −6601.8 ± 18.3 mas/yr, against GR's −6606.1.
  - Frame dragging: −37.2 ± 7.2 mas/yr, against GR's −39.2. — [arXiv:1105.3456](https://arxiv.org/abs/1105.3456)
- Ciufolini, Paris, Pavlis, Ries, Matzner, Paolozzi, Ortore, Bianco, Kuzmicz-Cieslak, Gurzadyan and Penrose (2023), "First results of the LARES 2 space experiment to test the general theory of relativity", Eur. Phys. J. Plus 138 (11), DOI 10.1140/epjp/s13360-023-04696-6 (arXiv:2311.13268).
  - LARES 2 was launched in 2022 on VEGA C. The results are "in complete agreement" with GR.
  - The abstract gives no numeric ratio. It describes the result as preliminary first evidence. — [arXiv:2311.13268](https://arxiv.org/abs/2311.13268); [MDSOAR record](https://mdsoar.org/items/3a5ebd8b-d10d-491b-b7d4-f7fd31480736)
- Ciufolini et al. (2026), "LARES-2 satellite measures frame-dragging effect around the Earth", Nature, DOI 10.1038/s41586-026-10715-0.
  - It combines LARES-2, LAGEOS and GRACE over about 3 years.
  - Relative uncertainty is about one part in a thousand, an order of magnitude better than before.
  - Extended data reportedly gives about 2×10⁻³ for the secular trend.
  - Dominant concerns are Earth's gravity field (the combination cancels even-zonal errors), lunisolar tides (K1) and non-gravitational forces. — [phys.org, July 2026](https://phys.org/news/2026-07-validate-einstein.html); [Nature](https://www.nature.com/articles/s41586-026-10715-0); [search summary incl. Extended Data](https://www.scientificamerican.com/article/einsteins-greatest-theory-triumphs-again-in-landmark-frame-dragging-measurement/)
- Critique and rebuttal:
  - Iorio (2025), "Will LAGEOS and LARES 2 succeed in accurately measuring frame-dragging?", Eur. Phys. J. C, DOI 10.1140/epjc/s10052-025-13964-x. Iorio argues that J₂ biases the result indirectly through errors in all the other parameters.
  - Ciufolini et al. (2024, EPJ C), "On the high accuracy to test dragging of inertial frames with the LARES 2 space experiment", reply that the critic's claims are "not reproducible in any independent analyses". They also note that the critic's earlier papers made opposite statements.
  - A 2019 Monte Carlo study gives a 1.3×10⁻³ systematic budget. — [EPJC 2025](https://link.springer.com/article/10.1140/epjc/s10052-025-13964-x); [Ciufolini et al. 2024 rebuttal](https://mdsoar.org/items/968436ce-e6dc-43d5-bf3e-0801e9fb4072/full); [arXiv:1911.06183](https://www.arxiv.org/abs/1911.06183); [arXiv:1910.10224](https://arxiv.org/pdf/1910.10224)
- Link to beams: for parallel pp-waves, the gravitoelectric force is exactly cancelled by the gravitomagnetic force (Faraoni & Dumse). Gyratons are boosted Lense–Thirring metrics (Frolov & Fursaev). — [gr-qc/9811052](https://arxiv.org/abs/gr-qc/9811052); [hep-th/0504027](https://arxiv.org/abs/hep-th/0504027)

### Inferences
- **Frame-dragging formulas for Rust [A]**, standard textbook GR checked numerically here:
  - Lense–Thirring node precession: `Ω̇_LT = 2 G J / (c² a³ (1−e²)^{3/2})`. With J⊕ = 5.86e33 kg m²/s, a = 12,270 km, e = 0.0045 (LAGEOS), this gives 30.7 mas/yr (computed). The LAGEOS + LARES-2 node sum doubles it while cancelling the classical J₂ precession, because the inclinations are supplementary.
  - Schiff gyroscope precession: `Ω = (3GM/(2c²r³))(r × v) + (G/(c²r³))[3(J·r̂)r̂ − J]`.
  - Simplified orbit averages for a circular polar orbit at 642 km: geodetic `(3/2)(GM)^{3/2}/(c² a^{5/2})` gives 6638 mas/yr; frame dragging `GJ/(2c²a³)` gives 41 mas/yr (computed). GP-B's full predictions of 6606 and 39.2 include Earth's oblateness and the guide-star geometry.
  - Gravitomagnetic field in linearised GR: `h₀ᵢ = −(4G/c³) ∫ ρvᵢ/|x−x′| d³x′` (sign convention dependent). This gives B_g ~ G J/(c² r³) for a spinning body.
- **Why this is the honest "gravity beam" story:** moving mass-energy really does produce velocity-dependent, magnetic-like gravity. It has been measured at the 0.1–20% level using the whole Earth as the source. A lab source (kg, m, kHz) is about 25 orders of magnitude weaker. Lab "gravitomagnetic" claims at 1e-4 ω (Section 4) are therefore physically implausible.

### Gaps
- I could not read the Nature 2026 paper itself (paywall). The exact measured ratio μ = Ω_obs/Ω_GR with its 1σ is not confirmed beyond "about 1e-3 relative uncertainty".
- The LARES-2 and LAGEOS orbital elements (inclinations about 70.16° and 109.84°) were not fetched. The figures above are recalled from memory and are approximate.
- I did not fetch earlier LAGEOS/LARES results (for example the 2016 and 2019 LARES at about 5% and 2%).
