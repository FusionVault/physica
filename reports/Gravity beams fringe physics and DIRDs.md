# Fringe gravity's sums reproduce, its physics never does

Almost every headline number in the fringe gravity literature can be reproduced from its authors' own formulas, usually to within a few percent. Almost every claim then fails at one identifiable physics step, by a margin that can be computed. Li and Baker's high-frequency gravitational-wave (HFGW) detector arithmetic matches its published values to 0–4%. But its claimed 3×10³ s integration time ignores the shot noise of the reference beam. The homodyne quantum limit needs **4×10¹⁷–3×10¹⁸ s, a factor of 10¹⁴–10¹⁵**, which the 2008 JASON panel called "more than 20 orders of magnitude" too optimistic in its own framing **[refuted]**. The other headline claims fail the same way:

- Baker's generator fluxes need a beam whose solid angle is **3×10¹¹ times** smaller than diffraction allows.
- Pais's patents need surface fields **5.6×10⁶ times** above the vacuum-breakdown limit his own patent states.
- Podkletnov's 2% weight loss is **7.7×10⁶ times** larger than general relativity (GR) predicts, and of the opposite sign.
- Tajmar's frame-dragging signal was **3.6×10²¹ times** GR before Tajmar's own group attributed it to vibration.

The nine gravity papers in the Defense Intelligence Agency's Defense Intelligence Reference Documents (DIRDs) quote textbook GR correctly, but get their own new numbers wrong, once by **37 orders of magnitude**. The documented classified gravity work is geodesy: the Navy's Trident gravity gradiometer and the Geosat marine-gravity data, both declassified in the 1990s. "Born secret" covers only nuclear Restricted Data. For the crates, the 2023–2026 mainstream work already has homes in `gwsens` and `satgrav::gw`: the high-frequency ceiling and new laboratory limits, atom-interferometer detectors, LARES-2 frame dragging and atmospheric waves in GOCE and GRACE-FO data. This pass adds about a dozen mainstream calculators and a small, labelled `claims` module, where every test pins the published number next to the GR number it must be compared with. Generators, shields, warp energies and zero-point devices stay out as physics.

*Labels.* **[mainstream]** means agrees with standard GR, electrodynamics or QFT. **[contested]** means published but disputed, with no settled refutation. **[refuted]** means a named null result, published refutation, or an arithmetic or dimensional error. The DIRD notes' [GR], [SPEC] and [FRINGE] map onto these labels. Numbers marked *(check)* were recomputed for this report with the three scripts in the notes folder: [baker_li_checks.py](../research_notes/Gravity%20beams%20fringe%20physics%20and%20DIRDs/baker_li_checks.py), [superconductor_gravity_checks.py](../research_notes/Gravity%20beams%20fringe%20physics%20and%20DIRDs/superconductor_gravity_checks.py) and [dirds_checks.py](../research_notes/Gravity%20beams%20fringe%20physics%20and%20DIRDs/dirds_checks.py), all run with `python3 -I`. They are our inferences, not published claims.

## The 2023–2026 mainstream work already has a home in two crates

Any fringe claim has to be read against the mainstream baseline, so that comes first. The reference for gravitational waves above 10 kHz is the 2025 Living Review by Aggarwal et al. (v2, February 2026). Its key number is a bound, not a sensitivity. Gravitational waves present before Big Bang nucleosynthesis count as extra radiation, so **Ω_GW h² < 1.1×10⁻⁶**. In the one-sided convention this means a strain of about **1.3×10⁻³⁰ at 1 GHz**. The review states that no proposal above the LIGO band reaches below this bound ([Review 2025](https://arxiv.org/abs/2501.11723)). A 2026 analysis finds that backgrounds above the MHz band lie below the quantum limit itself ([Guo et al.](https://arxiv.org/abs/2501.18146)) **[mainstream]**.

The new hardware limits are real, but they sit 18–31 orders of magnitude above that ceiling:

- QUEST set the first broadband 13–80 MHz limits, near 3×10⁻²⁰ Hz⁻¹ᐟ² ([Patra et al. 2025](https://arxiv.org/abs/2410.09175)).
- An ABRACADABRA-10cm transient search excludes primordial-black-hole sources only within 46 km ([Pappas et al.](https://arxiv.org/abs/2505.02821)).
- A CAPP-12T reanalysis reached h₀ ≈ 3.9×10⁻²¹ at 5.3 GHz ([Kim et al. 2025](https://arxiv.org/abs/2511.17817)).
- Laser-cavity spectroscopy constrains 100 kHz–100 MHz ([Budker et al. 2026](https://arxiv.org/abs/2607.12617)).
- The revived MAGO cavity found mechanical Q factors "well below commonly assumed theoretical values" ([Dokuyucu et al. 2026](https://arxiv.org/abs/2607.02476)).

The 2021 bulk-acoustic-wave "events" are not gravitational-wave candidates. A gravitational-wave reading would need Ω_GW ≈ 10⁸ ([Lasky & Thrane](https://arxiv.org/abs/2110.13319)) **[refuted as GWs]**.

**No atom-interferometer gravitational-wave detector operates in 2026.** MAGIS-100 expects commissioning by 2028 ([Fermilab](https://news.fnal.gov/tag/magis-100/)). AION-10 has a 2025 technical design report and "is not expected to detect GWs" ([AION-10 TDR](https://arxiv.org/abs/2508.03491)). AICE is a 2025 letter of intent for a 100 m shaft at CERN ([Baynham et al.](https://arxiv.org/abs/2509.11867)). The km-scale roadmap starts in 2035 ([TVLBAI](https://arxiv.org/abs/2503.21366)). The gating hardware is large momentum transfer: about 400ħk has been demonstrated ([Wilkason et al.](https://arxiv.org/abs/2205.06965)), against the 5,000–40,000ħk the designs assume.

The same tidal algebra converts any gradiometer to strain, h = 2√S_Γ/ω². That puts a GOCE-class instrument at 5×10⁻¹¹ at 0.1 Hz, 8–12 orders short of LISA ([Junca et al.](https://arxiv.org/abs/1902.05337)) **[mainstream]**.

The serious, mainstream version of a "gravity beam" is the near field of a light or particle beam. A test mass at distance ρ feels a = 4GP/(c³ρ), which is **9.9×10⁻¹⁸ m/s² for a petawatt beam at 1 mm** ([Rätzel et al.](https://arxiv.org/abs/1511.01023)). The best-case LHC experiment reaches S/N ≈ 0.6 in a week, before any technical noise ([Spengler et al.](https://arxiv.org/abs/2104.09209)) **[mainstream]**.

The measured relative of "moving mass makes a field" is frame dragging. LARES-2, combined with LAGEOS and GRACE gravity fields, measured it to about one part in a thousand in 2026 ([Ciufolini et al., Nature 2026](https://www.nature.com/articles/s41586-026-10715-0)) **[mainstream]**. Iorio's challenge to that error budget is **[contested]** ([Iorio 2025](https://link.springer.com/article/10.1140/epjc/s10052-025-13964-x)).

The "gravity waves" these satellites have actually recorded are atmospheric buoyancy waves:

- GOCE recorded Tohoku's infrasound as air-density changes of about 11% at 270 km ([ESA](https://www.esa.int/Applications/Observing_the_Earth/GOCE/GOCE_the_first_seismometer_in_orbit)).
- The 2022 Tonga Lamb wave circled the globe at 318 m/s ([Wright et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9492537/)).
- GRACE-FO and Swarm accelerometer densities show 200–450 m/s thermospheric waves ([Li et al. 2023](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2022GL102265)).

What causes those thermospheric waves is a live 2025 dispute. One model attributes them to secondary gravity waves ([Li et al. 2025](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2024AV001470)); another gives the Lamb wave about 25% ([Liu et al. 2023](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2023GL103682)) **[contested]**.

On the software side, Rust waveforms arrived on 22 October 2026 as the GPL-3.0 `waveformrs` crate ([crates.io](https://crates.io/crates/waveformrs)). The open gap is permissive sensitivity, signal-to-noise and HFGW tooling, which no language has.

The first pass of that plan is already in the working tree. `gwsens` has `strain`, `curves`, `sgwb`, `inspiral`, `hf` and `lab` modules. The `hf` module carries the ceiling constants, `gertsenshtein_probability`, coherent and stochastic photon rates, cavity noise and a `Status` registry. The `lab` module carries the beam near field, rotor and rod power, and `wave_flux`. `satgrav::gw` adds the Sr-87 wavenumber, the atom-interferometer response, and gradiometer, satellite-pair and clock strain. Two things are still missing. No Lense–Thirring code turned up in this workspace; the earlier plan placed it in the user's `geo5d`. The Tonga test against GRACE-FO ranging residuals is also outstanding, because it needs Level-1B readers.

This pass found a smaller set of 2023–2026 items on the fringe edge:

- a 2023 Li–Baker-type simulation that never mentions the JASON critique;
- the Living Review's §5.5.4 paragraph recording concerns about "heterodyne enhancement";
- the Dresden group's 2021–2022 null results on superconductor, Mach-effect and EmDrive claims;
- a December 2025 GR analysis of Mach-effect thrust;
- the 2024 subluminal warp shell and the MIT-licensed Warp Factory code;
- a September 2026 re-release of the DIRDs.

The sections below take each one in turn.

## Li and Baker's detector arithmetic holds to 4%, but its noise floor is off by 10¹⁴

### Every Li-group paper computes the same cross term

From 2003 to 2023, every paper by Fangyu Li's group computes one quantity, which they call the "perturbative photon flux" (PPF). It is the first-order cross term between two fields. One is the weak field a GW induces in a static magnet (the inverse Gertsenshtein effect). The other is a "synchro-resonant" 5–10 GHz Gaussian microwave beam. The standard parameter set is a 10 W beam with a 5 cm waist, 3 T over 6 m, 5 GHz and h = 10⁻³⁰. It yields **N_x^(1) ≈ 8.2×10² photons/s**, with a claimed required accumulation time of only **10³–10⁵ s** ([Li et al., EPJC 56, 2008](https://arxiv.org/abs/0806.1989); [Li et al., PRD 80, 064013, 2009](https://arxiv.org/abs/0909.4118)). Earlier work had predicted about 10³ s⁻¹ through 10⁻² m² at 10¹⁰ Hz ([Li, Tang & Shi, PRD 67, 2003](https://arxiv.org/abs/gr-qc/0306092)). A 2016 paper scans tensor indices that push the flux from 1.4×10⁻⁴ to 2.85×10⁷ s⁻¹ ([Li, Wang & Wen, CPC 2016](https://arxiv.org/abs/1509.02443)). The 2020 "3DSR" paper claims all six possible GW polarisations are separable ([Li et al., EPJC 80, 2020](https://arxiv.org/abs/1712.00766)).

The arithmetic holds up. From the authors' own formulas *(check)*:

- The PRD 2009 chain reproduces to within 0–4%. The second-order power is 2.39×10⁻⁴⁰ W against the published 2.3×10⁻⁴⁰. The photon rate is 7.2×10⁻¹⁷ s⁻¹ against 7.0×10⁻¹⁷. The plane-wave field is 868 V/m against 8.7×10², and the cross-term power is 6.91×10⁻²⁰ W against 6.9×10⁻²⁰.
- The graviton flux through the beam waist is 3.758×10¹⁶ s⁻¹ against the published 3.77×10¹⁶.
- The one outlier is the beam's peak field. The standard Gaussian-beam relation gives 1.385×10³ V/m, 10% above the published 1.26×10³, and the papers do not state their convention.

These reproductions need the papers' own rounded inputs, k = 100 m⁻¹ and hν = 3.3×10⁻²⁴ J.

The pure inverse-Gertsenshtein part is standard physics, and the authors say so. For 10 T over 10 m at h = 10⁻³⁰, the conversion ratio is 1.4×10⁻³³ and the rate 7×10⁻¹⁷ photons/s, one photon per 1.4×10¹⁶ s. In their own words, such a lab detector "would not be useful" ([PRD 2009](https://arxiv.org/abs/0909.4118)) **[mainstream]**. Their ratio is 1.8× the first-order SI result 4πGB²L²/(μ₀c⁴) = 8.3×10⁻³⁴ *(check)*.

### JASON's verdict: an engineering error, not a gravitation error

JASON's report JSR-08-506 (October 2008, lead author D. Eardley, for the ODNI and DIA) heard Baker, Stephenson and others and reviewed about 1000 pages. Its conclusions are four:

- the Li–Baker analysis is "incorrect by many orders of magnitude";
- detecting either relic or artificial HFGW is infeasible for the foreseeable future;
- no foreign HFGW threat is credible, whether for communication, radar, propulsion or any other use;
- weak-field GR is complete and reliable.

Its §6 says the homodyne sensitivity estimates are **"too optimistic by more than 20 orders of magnitude"**, and that this is an error of electrical engineering, not of gravitation physics ([JASON JSR-08-506](https://web.archive.org/web/2010/http://www.fas.org/irp/agency/dod/jason/gravwaves.pdf)) **[refuted]**.

The point of disagreement is narrow, and both sides write the same algebra. The detected count is N_d = Ṅ₀t + 2(Ṅ₀Ṅ_GW)^½ t cos δ + Ṅ_GW t. That is JASON eq. 3-32 and Li et al. eqs. 32–35. JASON counts the shot noise √(Ṅ₀t) of every reference photon that carries the interference term. With that noise the reference-beam rate cancels out of the result, leaving t > 1/Ṅ_GW: the standard quantum limit of homodyne detection. Li et al. instead put the receiver at x = 0, where the transverse Poynting flux of the beam vanishes by symmetry. They then count only "other" background (scattering, diffraction, thermal), giving Δt_min = N_other/N₁² ([PRD 2009 eqs. 50–58](https://arxiv.org/abs/0909.4118)). JASON calls neglecting the reference-beam noise "physically untenable": the interference term cannot exist unless E₀ is present at the same time.

Our check supports JASON's reading. The published "signal" at x = 3.5 cm (3.54×10² s⁻¹) equals 2√(N_BPF·Ṅ_GW) = **3.565×10² s⁻¹, a match to 0.7%**, where N_BPF is the local background flux of 1.24×10²² s⁻¹ *(check)*. So the PPF is the homodyne cross term. The 2025 Living Review states that linear-signal sensitivity is independent of the background field. It records the Woods (2012) and Eardley (2008) concerns about the Li–Baker scheme, notes that it relies on technology that does not exist, and says no experimental results suggest feasibility ([Review 2025 §5.3.1, §5.5.4](https://arxiv.org/abs/2501.11723)).

### How big the miss is

For the EPJC and PRD parameters (3 T, 6 m, 0.01 m², 5 GHz, h = 10⁻³⁰), the pure conversion rate is Ṅ_GW = 2.56×10⁻¹⁸ s⁻¹ with Li's field formula, or 3.20×10⁻¹⁹ s⁻¹ with the standard ½hBωL field. JASON's floor is therefore **t ≥ 3.9×10¹⁷ to 3.1×10¹⁸ s**, against Li's claimed 3×10³ s: **a ratio of 1.3×10¹⁴ to 1.0×10¹⁵**. Held at 3×10³ s, the homodyne limit allows only **h ≥ 1.1×10⁻²³, seven orders of magnitude worse** than the claimed 10⁻³⁰ *(check)*. JASON's own "10²⁵" is the factor 2E₀/E_GW by which the cross term exceeds |E_GW|² for its parameter set, which recomputes to 2.3×10²⁵ *(check)*.

O(1) convention factors explain every near-miss in the reproduction, and none matters against 10¹⁴:

- The Li and Baker papers use a GW energy flux c³ω²A²/(8πG), **4× the standard single-polarisation c³ω²h²/(32πG)** **[mainstream]**.
- The induced field appears as hBωL in PRD 2009 but carries a ½ in the 3DSR paper.
- JASON's Gaussian-units conversion F = 4πGB²L²/c⁴ is 4π higher than the first-order result GB²L²/c⁴.

JASON's own arithmetic slips all err in the proponents' favour *(check)*:

- its E_GW is 9×10⁻²² V/m, not 10⁻²² V/m;
- its copper cavity gives 1.8×10⁻³ graviton/s, not 2×10⁻²;
- its stored energy is P·Q/ω = 0.32 J, not about 2 J.

### The other side of the record

Baker's public rebuttal says JASON found the work wrong "by an order of magnitude". JASON's text says "many" and "more than 20" ([Baker Q&A](https://www.drrobertbaker.com/docs/Q%20&%20A%20JASON.pdf)). The rebuttal's characterisation is therefore **[refuted]** by the report itself. A full-text search of the 2008, 2009, 2016 and 2020 Li-group papers finds no occurrence of "JASON" ([PRD 2009](https://arxiv.org/abs/0909.4118); [EPJC 2020](https://arxiv.org/abs/1712.00766)).

The line continues into this period. Wu, Li and Jiang (Chinese Physics C, 2023) simulate a modulated-beam Li–Baker detector at h = 10⁻³¹. They report an SNR of 259–320 at chosen transverse positions and argue that the background "can be eliminated" there, without citing JASON ([Wu et al., CPC 2023](https://cpc.ihep.ac.cn/article/doi/10.1088/1674-1137/ace9c3)). The JASON argument applies unchanged, because an SNR that ignores the reference beam's shot noise is not a quantum-limited SNR. So the paper is **[refuted]** by an argument it does not address. That reading is our inference.

## Baker's and Pais's generators need beams 10¹¹ too narrow and fields 10⁶ too strong

### Baker's formula is right; his multipliers are not

Robert Baker's "jerk" formula, P = 1.76×10⁻⁵²(2rΔf/Δt)² W, uses the correct constant: 32G/(5c⁵) = 1.764×10⁻⁵² *(check)*. For a circular binary it equals the Peters result exactly ([Baker, Li & Li 2006](https://www.drrobertbaker.com/docs/AIP;%20HFGW%20Laser%20Generator.pdf)) **[mainstream]**. For a single piezoelectric film bulk acoustic resonator (FBAR) pair it gives 1.10×10⁻³³ W. The standard quadrupole result for a linearly driven pair, (16G/15c⁵)r²F₀²ω², gives 1.81×10⁻³³ W, a ratio of 0.61 *(check)*.

The headline numbers come from three multipliers, each departing from mainstream physics:

**Coherent N² summation.** The 2012 double-helix array has N = 1.55×10¹³ FBAR pairs over 20 m. Summed coherently it gives 2.64×10⁻⁷ W ([Baker & Baker, Physics Procedia 38, 2012](https://inspirehep.net/files/6bf245c473964d406575a245362a0ede)); the standard quadrupole sum gives 4.4×10⁻⁷ W. But the full N² gain holds only when all sources sit within about λ/2, and this array is 0.7 m × 20 m against a 6 cm wavelength. Incoherently the same array radiates **2.8×10⁻²⁰ W** *(check)* **[contested]**.

**A sub-diffraction "needle" beam.** Baker's cap area of 1.71/N m² at 1 m implies a half-angle of 1.9×10⁻⁷ rad. Any wave obeying the free wave equation from a 0.7 m aperture at 6 cm spreads by at least 1.22λ/D ≈ 0.10 rad. **The solid angle, and so the flux, is overstated by 3×10¹¹** *(check)*. The diffraction-limited flux at 1 m is about 1.3×10⁻⁵ W/m², not 6×10⁵ W/m². That is a strain of about 1.8×10⁻³⁰ rather than the published 1.98×10⁻²⁵. The label is **[contested]**: our diffraction check refutes the beam, but no published rebuttal exists.

**A strain relation off by 2.** Baker's A = (8πGF/c³ω²)^½ = 1.28×10⁻¹⁸√F/ν reproduces his own numbers: 4.80×10⁻²⁸ against the published 4.9×10⁻²⁸ for the FBAR array ([Baker 2009](https://www.drrobertbaker.com/docs/HFGW%27s%20Power%20Requirements.pdf)). The standard single-polarisation relation gives 9.60×10⁻²⁸ for the same flux *(check)*.

Two internal checks also fail. The 2006 laser paper's eq. 1 equates a rotating-rod expression with the dumbbell jerk form, and the two sides differ by **a factor of 144** *(check)*. Baker's PSR 1913+16 "calibration", 1.01×10²⁵ W, does not reproduce: Peters–Mathews gives 7.75×10²⁴ W, and his own jerk form gives 2.06×10²⁵ W *(check)*.

JASON's ceilings for electromagnetically driven generators are:

| Generator | Power | Rate |
|---|---|---|
| Copper cavity | ~2×10⁻²⁷ W | ~10⁻³–10⁻² graviton/s |
| Superconducting cavity | ~10⁻²³ W | |
| Mechanical limit | ~10⁻²⁰ W | |
| Container-explosion limit | ~10⁻¹⁸ W | |

For propulsion, 1 g from a GW focus 100 m ahead needs **10³⁸ W, more than 10²⁵ times world power use** ([JASON](https://web.archive.org/web/2010/http://www.fas.org/irp/agency/dod/jason/gravwaves.pdf)) **[mainstream]**. Baker replies that his generators are mechanical, not Gertsenshtein-type ([Baker Q&A](https://www.drrobertbaker.com/docs/Q%20&%20A%20JASON.pdf)). The quadrupole numbers above answer that on his own terms.

### DIRD 21 restates Baker eighteen months after JASON

Baker's DIRD 21, "High-Frequency Gravitational Wave Communications" (April 2010), repeats the jerk formula: 2.82×10⁻¹³ W for r = 10 m, Δf = 4×10⁸ N and Δt = 2×10⁻¹⁰ s. It makes no rebuttal of JASON. It also adds two arithmetic errors *(check)*:

- Its eq. 2 prints A ∝ √(F/ν). The worked examples use the dimensionally correct √F/ν, and only that form reproduces 2×10⁻³³ at 7×10⁻¹¹ W/m².
- Its Shannon capacity of "1.9×10⁶ bps" is B·ln(1+S/N). With log₂, **B·log₂(1+S/N) = 2.82×10⁶ bps** ([DIRD 21](https://documents2.theblackvault.com/documents/dia/AAWSAP-DIRDs/DIRD_21-DIRD_High-Frequency_Gravitational_Wave_Communications.pdf)) **[refuted, arithmetic]**.

The link budget's physical inputs have no empirical support: the 10¹⁴ W/m² generator, the N² gain and the 10⁻³² detector. JASON had already judged that no practical use of HFGW is credible **[refuted]**.

### Pais's patents need fields beyond his own breakdown limit

Salvatore Pais's Navy patents rest on one expression, S_max = f_G(σ²/ε₀)·v, where σ is a surface charge density and v a velocity-like term. It reproduces his headline numbers. With σ = 5×10⁴ C/m² on a 2 m disc at 30,000 rpm, the craft patent's example gives **1.77×10²⁴ W/m²** (published "order 10²⁴") ([US 10,144,532](https://patents.google.com/patent/US10144532B2/en)). The HFGW patent's eq. 2 gives 1.1×10³² against the claimed 10³³ ([US 10,322,827](https://patents.google.com/patent/US10322827B2/en)) *(check)*.

But the implied surface field σ/ε₀ = **5.6×10¹⁵ V/m is 5.6×10⁶ times the 10⁹ V/m vacuum limit** stated in his own US 10,135,366 ([patent](https://patents.google.com/patent/US10135366B2/en)). Capped at that limit, the same formula gives 5.6×10¹⁰ W/m², **3×10¹³ times lower** *(check)*. The expression has no c in it, grows without bound with operating time, and has no link to mechanical input power. A uniformly spinning, axisymmetric charge is a static field configuration and radiates nothing **[refuted]**. That last argument is our inference; the press critiques are qualitative.

The HFGW patent's eq. 4 uses JASON's Gaussian-units formula with SI inputs. Pais is right that the SI version needs 1/μ₀. But the SI-correct 8.3×10¹⁵ W needs B = 10¹⁰ T, an energy density B²/2μ₀ = 4.0×10²⁵ J/m³ *(check)*. The craft patent's eq. 4 equates ε₀E² at the Schwinger field to a cut-off zero-point density, giving ν = 4.2×10²¹ Hz. That equation has no basis in Schwinger pair production, whose rate scales as exp(−πE_S/E).

The Navy's record closes the question empirically. The USPTO examiner first rejected the craft patent, noting it would need about 10⁹ T and 10¹⁸ V/m. It was allowed after the Naval Aviation Enterprise's chief technology officer warned of Chinese work ([TWZ 2019](https://www.twz.com/28729/docs-show-navy-got-ufo-patent-granted-by-warning-of-similar-chinese-tech-advances)). NAWCAD then tested the "High Energy Electromagnetic Field Generator" from October 2016 to September 2019, for **$508,000**, and the "Pais Effect" "could not be proven" ([TWZ, NAWCAD statement](https://www.twz.com/39012/the-navy-finally-speaks-up-about-its-bizarre-ufo-patent-experiments)). A September 2018 spin test of a capacitor at 100,000 rpm "did not measure any anomalous electromagnetic effects" ([TWZ FOIA](https://www.twz.com/38937/navy-ufo-patent-documents-talk-of-spacetime-modification-weapon-detail-experimental-testing)). No test of the HFGW generator itself was reported **[refuted as tested]**.

## Superconductor gravity claims sit 10⁶ to 10²⁹ above general relativity

### Podkletnov's shield has the wrong sign and the wrong geometry

Podkletnov reported weight losses of 0.05–0.3% in 1992. In 1997 he reported up to **2.1%**, for 25–30 s while braking a levitated 275 mm two-layer YBCO toroid from 5000 rpm under 3–4 MHz fields. He gave the same percentage for every material, and the same loss 3 m above the disc ([Podkletnov 1997](https://arxiv.org/abs/cond-mat/9701074)). The paper gives no formula.

GR makes a definite prediction. At an assumed 6.38 g/cm³ the disc's surface density is 63.8 kg/m², so its own pull *adds* to the weight by at most 2πGσ/g = 2.7×10⁻⁹. That is **7.7×10⁶ times smaller than the claim, and of the opposite sign**. The rotation-dependent part is about 4×10⁻²³, so the claimed rotation effect exceeds GR by **10²⁰–10²¹** *(check)*.

Geometry adds a second contradiction. For a uniform Earth, the fraction of g arriving through a downward cone of half-angle θ is 1 − cos³θ. That reproduces Modanese's own "shadow" table to within about 20% ([Modanese, supr-con/9601001](https://arxiv.org/abs/supr-con/9601001)). At 3 m, 22 disc radii up, the fraction is 0.31%. So a 2% shield near the disc should read about 0.007% there, not 2% *(check)*. Modanese himself called the missing height dependence "a severe challenge".

Every replication was null:

| Null test | Bound on Δg/g | Below the 2.1% claim by |
|---|---|---|
| NASA MSFC static, Li et al. 1997 | < 2×10⁻⁸ | ≈10⁶ |
| Noever et al. 1998, rotating field | ≲ 10⁻⁶ | ≈2×10⁴ |
| Tajmar et al. 2006, accelerometers | < 5×10⁻⁶ | ≈4×10³ |
| Hathaway, Cleveland & Bao 2003 | < 10⁻⁵ | ≈2×10³ |
| Woods et al. 2001 (BAE/Sheffield) | ±4×10⁻⁵ | ≈500 |

Sources: [NTRS 19990039542](https://ntrs.nasa.gov/citations/19990039542); [NTRS 19990019627](https://ntrs.nasa.gov/citations/19990019627); [Tajmar et al.](https://arxiv.org/abs/gr-qc/0603033); [Gallerati, Modanese & Ummarino 2022](https://arxiv.org/abs/2206.07574); [Woods](https://arxiv.org/abs/physics/0204031).

NASA's RF-illuminated test saw a 3–5×10⁻⁸ g *increase* that tracked the AC-resistance peak, an electromagnetic artefact ([NTRS 19990064220](https://ntrs.nasa.gov/citations/19990064220)). NASA's Breakthrough Propulsion Physics (BPP) review judged the approach "non-viable" ([Millis, NASA/TM-2004-213082](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20040070788.pdf)). GR's +3×10⁻⁹ is consistent with every null.

One caveat is honest. No replication combined every prescribed condition: the two-layer toroid, magnetic levitation and rotation, the 10⁵ Hz drive plus MHz fields, and temperatures below 70 K. So the shield is "not reproduced" rather than falsified under identical conditions **[refuted for every configuration tested]**. The burden stays with the claimant, who never released data. Checking would be cheap: a 0.05% change is about 490 mGal, enormous next to µGal gravimeter noise.

### The impulse generator's tables are pendulum geometry

The 2001 "impulse gravity generator" reported 18.5 g pendulums deflected by up to 142 mm after 2 MV discharges, with the same result at 150 m in another building ([Podkletnov & Modanese](https://arxiv.org/abs/physics/0108005)). Every row of both tables follows from geometry alone: Δh = l − √(l² − Δl²) and ΔE = mgΔh, with l = 0.8 m. At 2 MV this gives Δh = 12.70 mm, ΔE = 2.305×10⁻³ J and a kick of 0.499 m/s, or "about 10³ g" over 10⁻⁴ s (509 g exactly) *(check)*. No physics model enters.

The physics then fails four ways *(check)*:

- The pendulum's momentum is **1.2×10⁹ times E/c**, but any radiation, electromagnetic or gravitational, carries p = E/c.
- The Marx bank stores 2.5 kJ and can deliver 2.5 mC, not the paper's 0.1 C and 10⁵ J, a factor of 40.
- Treating the whole bank energy as mass-energy, GR gives a velocity change of 5×10⁻³⁰ m/s, **10²⁹ below the claim**.
- The 2012 "64c" speed (63 ns over 1211 m) means Δv/c = 63. The GW170817 bound, |Δv/c| < 7×10⁻¹⁶, is about **9×10¹⁶ times** tighter ([LIGO/Virgo/Fermi/INTEGRAL](https://arxiv.org/abs/1710.05834); the notes quote this bound without re-fetching it). Clock synchronisation over 1.2 km is the obvious systematic.

Scaled-down repeats by Tajmar and Lőrincz traced their accelerometer signals to acoustic noise ([Gallerati et al. 2022](https://arxiv.org/abs/2206.07574)) **[contested, unsubstantiated]**.

Modanese's theory reproduces only in order of magnitude. It needs a local cosmological constant **10²⁰ times** the observed value, and a vacuum energy density 10^15.8 times the observed one. It makes no prediction of the weight-loss percentage that could be checked *(check)*; ([Modanese, EPL 1996](https://arxiv.org/abs/hep-th/9505094)) **[contested]**.

### Tajmar's frame dragging: 10²¹ above GR, then withdrawn by Tajmar

In 2006 de Matos and Tajmar reported tangential accelerations of about 100 µg next to a niobium ring accelerated at 1500 rad/s², a coupling of −6.6 to −9.6×10⁻⁸ s². They proposed a "gravitomagnetic London moment" B_g = 2ωρ*/ρ, from a Proca-type graviphoton mass ([Tajmar et al.](https://arxiv.org/abs/gr-qc/0603033); [Tajmar & de Matos](https://arxiv.org/abs/gr-qc/0603032)).

The coefficient reproduces: m_e²/(μ₀e²g₀) = 2.623×10⁻¹⁸ against the published 2.624×10⁻¹⁸ *(check)*. But the two papers silently use Nb London depths of about 28 nm and 42 nm. The measured frame-dragging ratio is k = Ω_drag/ω ≈ 2.6×10⁻⁵. GR for the same 0.349 kg ring gives k = 4GM/[c²(R_o + R_i)] = **7.2×10⁻²⁷**, so the claim is **3.65×10²¹ (10^21.6) times GR**, or 10^21.4–10^21.9 depending on convention *(check)*. The paper's "30 orders" and ESA's "10²⁰" ([ESA 2006](https://www.esa.int/esaCP/SEM0L6OVGJE_index_0.html)) do not reproduce from the stated setup. A GR source with that k would need white-dwarf compactness.

Tajmar's own group closed the case:

- By 2007 the signal did not track the superconducting transition and showed "parity violation". The data excluded the Cooper-pair hypothesis by 5 orders ([Tajmar et al. 2007](https://arxiv.org/abs/0707.3806)).
- A Canterbury ring-laser test was null to within 1σ; Tajmar's signal would have been 11σ there ([Gallerati et al. 2022](https://arxiv.org/abs/2206.07574)).
- In 2011 an improved gyroscope set a limit of **4×10⁻¹¹**, 800× below the 2007 signal and 10⁵ below the prediction ([Tajmar, SST 24, 125011](https://iopscience.iop.org/article/10.1088/0953-2048/24/12/125011)).
- The earlier signals were attributed to "acoustic noise and vibration artefacts" ([Tajmar et al. 2022](https://www.frontiersin.org/articles/10.3389/fphy.2022.892215/full)).

**[Refuted, including by the original author.]** Even the 2011 null limit sits 5.6×10¹⁵ above GR. A laboratory measurement of frame dragging remains about 16 orders out of reach *(check)*.

### Ning Li's equations are mainstream and negligible; the enhancement is not

Li and Torr's 1991 PRD equations are standard linearised gravito-electromagnetism (GEM) bookkeeping: inside a superconductor B + (m/q)B_g → 0 ([Li & Torr, PRD 43, 457](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.43.457)). Taken at face value they give an induced B_g of **4.2×10⁻³² s⁻¹ per tesla**, and a residual field from Earth's gravitomagnetism of about **5.7×10⁻²⁶ T** *(check)* **[mainstream, unmeasurable]**.

The 1993 claim was that lattice ions, not Cooper pairs, carry the angular momentum and yield detectable fields. Harris showed those fields were overestimated "by many orders of magnitude" through "unrealistic assumptions" ([Harris 1999](https://inspirehep.net/literature/1359113)) **[refuted, theory]**. The exact factor in Harris's paper was not retrieved. The NASA static test with Li as lead author found < 2×10⁻⁸. Her company AC Gravity LLC took a **$448,970** DoD grant in 2001–2002 and published no results ([Wikipedia](https://en.wikipedia.org/wiki/Ning_Li_(physicist))).

### The serious theory predicts nothing a gravimeter can see

Ummarino and Gallerati couple time-dependent Ginzburg–Landau theory to gravito-Maxwell equations ([EPJC 77, 549](https://arxiv.org/abs/1710.01267); [Gallerati & Ummarino 2022](https://arxiv.org/abs/2203.09417)). Their numbers reproduce *(check)*: g* = 1.92×10⁻¹¹ for YBCO at 77 K, and the Josephson period of 123.6 s for a 300 nm junction. They also expose a hidden **4π change** in η between the 2017 and 2021–22 papers.

Their gravity "reduction" is transient. It decays within ητ ≈ **4.3 ps** for YBCO, so a static weighing experiment of the Podkletnov kind predicts *zero* steady weight change. The theory also conflicts with mainstream physics. In equilibrium, gravity on electrons in a metal is balanced by the tiny Schiff–Barnhill field m_e g/e = 5.6×10⁻¹¹ V/m ([Schiff & Barnhill 1966](https://doi.org/10.1103/PhysRev.151.1067)). Merging (m/e)E_g into a "generalised" electric field makes gravity screenable like charge, which violates the universality of free fall, because m/e differs between species **[contested]**.

The mainstream baselines are clean. The London moment B = 2m_eω/e = 1.137×10⁻¹¹ T per rad/s is measured; it is 7.39 nT at 650 rad/s. GR shifts it by a fraction of 7.2×10⁻²⁷. DeWitt's gravitational-drag current is about 8×10⁻²⁹ A, far below the femtoamp floor of a SQUID. The gravitational penetration depth for niobium is 8×10¹³ m ([Kiefer & Weber](https://arxiv.org/abs/gr-qc/0408010)) *(check)* **[mainstream]**.

The Tate Cooper-pair mass anomaly is a real but shrinking condensed-matter discrepancy, not evidence for gravity **[mainstream, open]**. It was 92 ppm at 4.4σ in 1989 ([Tate et al.](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.62.845)) and is now 18 ± 2.1 ppm for Nb (Hoang et al. 2020, via [Tajmar et al. 2022](https://www.frontiersin.org/articles/10.3389/fphy.2022.892215/full)).

## The public record shows unclassified relativity, null tests and declassified geodesy

### Seventy years of "gravity control" end in ion wind and thermal drift

The 1956 "Electrogravitics Systems" report (GRG-013/56) is a compilation of trade-press claims about Convair, Lear, Sperry, Glenn Martin and others "studying" gravity. Its "CONFIDENTIAL" marking rests on an editor's note in a fringe reprint ([Valone scan](https://archive.org/details/thomas-valone-electrogravitics-systems)) **[contested, unverified]**.

What the record documents is ordinary science. Martin's research institute RIAS hired relativists whose output was standard GR ([TWZ](https://www.twz.com/30499/the-truth-is-the-military-has-been-researching-anti-gravity-for-nearly-70-years)). The USAF funded GR from 1956 until the Mansfield Amendment of 1969–70, and, as far as the historian Goenner knows, "the supported projects were not classified" ([Goenner](https://arxiv.org/abs/1607.03319)) **[mainstream]**. The USAF's 1972 "Project Outgrowth" found that the advanced propulsion ideas "lack theoretical and technical merit".

T. T. Brown's asymmetric capacitors ("lifters") are fully explained by ion drift, **F = I·d/μ**. With 1 mA over 4 cm and a nitrogen-ion mobility of 2.5×10⁻⁴ m²/(V·s), that gives 0.16 N, the right order to lift a lifter ([Bahder & Fazi, ARL-TR-3005](https://arxiv.org/abs/physics/0211001)). NASA's collision model is the same law in other clothes. It implies a mobility of about 1.7×10⁻⁴ ([Canning et al., NASA/CR-2004-213312](https://ntrs.nasa.gov/api/citations/20040171929/downloads/20040171929.pdf)) *(check)*.

The vacuum tests are all null:

- Talley's torsion pendulum (USAF, 1991) saw no force in vacuum.
- Canning et al. saw no rotation below 300 Torr.
- Tajmar's sealed Faraday box put any anomaly **at least five orders below corona wind** ([Tajmar 2004](https://tu-dresden.de/die_tu_dresden/fakultaeten/fakultaet_maschinenwesen/ilr/rfs/forschung/folder.2007-08-21.5231434330/ag_raumfahrtantriebe/Biefeld-Brown%20Effect%20AIAA%20Journal%20Revised.pdf)).

**[Refuted.]**

NASA's BPP programme (1996–2002) ran 16 tasks. Six were not viable, six were unresolved and four merited sequels ([Millis 2004](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20040070788.pdf)). BAE's Project Greenglow had been public since the 2002 Jane's report and that 2004 NASA memo, long before its 2015 book and 2016 BBC programme. By the accounts available it produced no breakthrough. Boeing's GRASP is known from one leaked-briefing report, and Boeing denied company funding ([Wikipedia: Podkletnov](https://en.wikipedia.org/wiki/Eugene_Podkletnov)).

The two "unresolved" BPP threads were later closed in Dresden:

- **Mach-effect thrusters.** Woodward's own theory, ⟨F⟩ = ω⁶m₀K_p²K_e x₀³V₀⁴/(8πGρ₀c²), predicts **7.0 nN** for his stack *(check)*, against about 100 µN reported ([Fearn & Woodward](https://arxiv.org/abs/1301.6178)). Dresden balances with < 5 nN noise found only vibration artefacts ([Monette, Kößling & Tajmar 2021](https://ui.adsabs.harvard.edu/abs/2021AcAau.182..531M/abstract)). A December 2025 preprint finds the effect suppressed in covariant GR by (U/c²)(ωL/c)² ([Rodal](https://arxiv.org/abs/2512.22687)). The thrust is **[refuted as tested]** and the theory **[contested]**.
- **EmDrive.** Eagleworks reported 1.2 mN/kW, **360× the photon-rocket limit of 1/c = 3.34 µN/kW** *(check)*. Dresden's battery-powered test saw nothing above a **10–30 nN** floor at 11 W, at or below the photon limit. It reproduced the White-like signal as thermal expansion stressing the balance pivots ([Tajmar, Neunzig & Weikert, CEAS 2022](https://d-nb.info/1244138614/34)) **[refuted]**.

### Secrecy on the record is real, but none of it is about gravity control

Invention-secrecy orders are counted publicly. There were **6,543 in effect at the end of FY2025**, and 356 new orders in FY2024, of which 94 were "Type 2" classified DoD orders. No subject breakdown is published ([FAS](https://sgp.fas.org/othergov/invention/stats.html)). Linking them to propulsion or gravity would be speculation.

"Born secret" is the Atomic Energy Act's Restricted Data doctrine, which covers nuclear weapons and special nuclear material and was tested in *United States v. Progressive* ([Wikipedia](https://en.wikipedia.org/wiki/Born_secret)). No source in these notes links it to gravity research. Senator Reid's June 2009 request for Special Access Program status for AAWSAP was declined ([AARO Vol. I](https://documents2.theblackvault.com/documents/dod/AARO_Historical_Record_Report_Volume_1_2024.pdf)). JASON's HFGW report does carry a classified appendix. The Pais episode ran the other way: a service lab used its CTO to push non-working concepts into print, and tested them only afterwards.

The classified-then-declassified gravity technology that *is* documented is the same technology this repository uses:

- **The Trident gradiometer.** Bell Aerospace's rotating-accelerometer gravity gradiometer was built for the Navy's Trident submarines, to compensate their inertial navigation for gravity. It was declassified in the early 1990s and became BHP's Falcon airborne gradiometer (first survey October 1999) and Bell Geospace's FTG ([Bell, Scientific American 1998](https://www.scientificamerican.com/article/gravity-gradiometry); [Dransfield 2007](https://www.911metallurgist.com/wp-content/uploads/2015/10/Airborne-Gravity-Gradiometry-in-The-Search-for-Mineral-Deposits.pdf)).
- **Geosat.** The Navy classified Geosat's geodetic altimetry. It was released in stages and fully declassified on **22 June 1995**, after ERS-1's geodetic phase. That release enabled the Sandwell–Smith marine gravity grids ([Sandwell et al.](https://www.springerprofessional.de/marine-gravity-from-satellite-altimetry-over-ocean-and-sea-ice/13386424); [Wikipedia: Geosat](https://en.wikipedia.org/wiki/Geosat)) **[mainstream]**.

The gradiometer numbers are implementable and testable. Lee's requirement of 1 ppb accelerometer matching follows from Γ·L/a = (10 E)(0.1 m)/(1 m/s²) = 10⁻⁹, where 1 E = 10⁻⁹ s⁻². Falcon's 6 E/√Hz at 55 m/s, filtered to 300 m, gives 2.6 E RMS, against the 2.5 E published ([Dransfield 2007](https://www.911metallurgist.com/wp-content/uploads/2015/10/Airborne-Gravity-Gradiometry-in-The-Search-for-Mineral-Deposits.pdf)) *(check)*.

## The gravity DIRDs quote textbook relativity correctly and compute new numbers wrongly

### Nine documents, mostly from one firm

Of the 37 DIRDs the DIA released in March 2022 (the 38th is a classified laser-weapons report), nine concern gravity or spacetime. They are DIRDs 13, 14, 15, 18, 19, 21, 24, 29 and 36 ([Black Vault](https://www.theblackvault.com/documentarchive/the-advanced-aerospace-weapon-system-applications-program-aawsap-documentation/)). Author names are redacted. The attributions come from a 2018 list published by George Knapp: five to Eric Davis, and one each to Puthoff, Obousy, Hathaway, Baker and Winterberg ([parabunk transcription](http://parabunk.blogspot.com/2018/07/aawsap-dirds.html)). A September 2026 Pentagon release reportedly included the same 37 documents; the official page was not found ([Gray News](https://www.wdbj7.com/2026/09/19/pentagon-releases-another-batch-ufo-files-including-documents-secretive-government-program/)).

The DIRDs reproduce wherever they quote standard results. Obousy's Lobo–Visser warp energies, E = −(c⁴/G)v²R²σ, match on every row: −3.03×10⁵⁰ J at v = 1 for R = 50 m ([DIRD 13](https://documents2.theblackvault.com/documents/dia/AAWSAP-DIRDs/DIRD_13-DIRD_Warp_Drive_Dark_energy_and_the_Manipulation_of_Extra_Dimensions.pdf)). Davis's wormhole results also match ([DIRD 18](https://documents2.theblackvault.com/documents/dia/AAWSAP-DIRDs/DIRD_18-DIRD_Traversible_Wormholes_Stargates_and_Negative_Energy.pdf)):

- the Morris–Thorne–Yurtsever Casimir plate spacing, 4.05×10⁻¹⁸√r m;
- the wormhole mass, −1.35×10²⁷ kg per metre of throat;
- Ford's objection to gravitational squeezing, with λ_c = 2.42×10¹¹ m.

Puthoff's eq. 4 is exactly the Schwarzschild metric in harmonic coordinates ([DIRD 15](https://documents2.theblackvault.com/documents/dia/AAWSAP-DIRDs/DIRD_15-DIRD_Advanced_Space_Propulsion_Based_on_Vacuum_Spacetime_Metric_Engineer.pdf)). Davis's levitation energy of 62.5 MJ/kg checks ([DIRD 19](https://documents2.theblackvault.com/documents/dia/AAWSAP-DIRDs/DIRD_19-DIRD_Antigravity_for_Aerospace_Applications.pdf)). All of these are *(check)* **[mainstream]**.

### Where they compute something new, it goes wrong

The original calculations contain errors peer review would have caught *(check)*:

- **DIRD 13, eq. 2.4.** Inverting the speed bound for v = 1 gives R = 8.6×10⁻¹³ m. The document says 3.37×10²⁴ m, which is R²/(GM/c²), a dimensionally wrong expression and **an error of about 37 orders of magnitude** **[refuted, arithmetic]**.
- **DIRD 13, Table 3.** The "dimensional warp" table quotes 10⁴² J for a bubble said to hold 100 m³ at 10⁴² J/m³, which is 10⁴⁴ J. Its v = 100 row breaks the table's own v² scaling.
- **DIRD 18.** The gravitational-squeezing density recomputes as −2.53×10⁻²², a 13% mismatch. Davis quotes Ford's objection, which shifts the density by **48 orders of magnitude**, and calls it an "open question" **[contested]**. His flat-faced "stargate" also omits Visser's result: the exotic matter does not vanish but moves to the cube's edges, as negative-tension strings of −c⁴/(8G) = −1.51×10⁴³ J/m. A 1 m cube needs −2.0×10²⁷ kg of them ([Visser 1989](https://arxiv.org/abs/0809.0907)).
- **DIRD 24.** The Planck-cutoff zero-point energy of "10¹¹³ J/m³" lies between the angular (5.9×10¹¹¹) and cyclic (9.2×10¹¹⁴) conventions. Either is **10¹²¹–10¹²⁴ above the observed 5.7×10⁻¹⁰ J/m³**. That gap is the textbook cosmological-constant problem **[mainstream]**, not an energy source. The extraction devices are **[refuted/unsupported]** ([DIRD 24](https://documents2.theblackvault.com/documents/dia/AAWSAP-DIRDs/DIRD_24-DIRD_Concepts_for_Extracting_Energy_from_the_Quantum_Vacuum.pdf)).

Two DIRDs hold up better than their reputation. Hathaway's DIRD 14 is a candid negative review: interactions between gravity and electromagnetism are "many orders of magnitude smaller than those required" ([DIRD 14](https://documents2.theblackvault.com/documents/dia/AAWSAP-DIRDs/DIRD_14-DIRD_The_Role_of_Superconductors_in_Gravity_Research.pdf)). DIRD 36's homodyne tomography of negative-energy states has since been realised as sub-cycle quantum electrodynamics, with squeezed vacuum sampled below the vacuum noise level ([Riek et al., Nature 2017](https://ideas.repec.org/a/nat/nature/v541y2017i7637d10.1038_nature21024.html)) **[mainstream]**.

The official verdict is AARO's 2024 history. It says the papers "were never thoroughly peer reviewed", and that DIA cancelled the programme in 2012 "due to lack of merit and the utility of the deliverables" ([AARO Vol. I](https://documents2.theblackvault.com/documents/dod/AARO_Historical_Record_Report_Volume_1_2024.pdf)). Sean Carroll's comment on DIRD 13 runs the same way: "real physics", but not something that will "connect with engineering anytime soon, probably anytime ever" ([ScienceAlert](https://www.sciencealert.com/here-s-what-theoretical-physicists-think-of-the-us-military-s-warp-drive-research)).

### Puthoff's metric is a clean test of any geodesic code

Puthoff's polarisable-vacuum metric, K = exp(2GM/rc²), is not in DIRD 15. It is in his 2002 Foundations of Physics paper ([Puthoff 2002](https://arxiv.org/abs/gr-qc/9909037)), and it is exactly the exponential metric. It matches GR at first post-Newtonian order: 1.7512″ of solar light bending, 42.98″/century for Mercury, and Pound–Rebka's 2.456×10⁻¹⁵ *(check)*.

Integrating the null geodesics numerically separates the two at second order. The exponential metric gives a second-order coefficient of **4π (12.5658)**; Schwarzschild gives **15π/4 (11.7804)** *(check)*. For a Sun-grazing ray the difference is only 0.73 µas, below current VLBI and Gaia precision. In the strong field the exponential metric is a horizonless traversable wormhole, with its throat at r = GM/c² and areal radius e·GM/c² (4014 m for one solar mass) ([Boonserm et al. 2018](https://arxiv.org/abs/1805.03781)). Black-hole observations are therefore the decisive test; that inference is ours **[contested]**.

### Warp-drive physics since 2019 tightened the DIRDs' case against them

**Superluminal warp drives violate the null energy condition for some observer** ([Santiago, Schuster & Visser 2022](https://arxiv.org/abs/2105.03079)) **[mainstream]**. That makes the 2021 "positive-energy" solitons valid only for Eulerian observers **[contested]** ([Lentz](https://arxiv.org/abs/2006.07125); [Fell & Heisenberg](https://arxiv.org/abs/2104.06488)).

The legitimate new territory is subluminal. Bobrick and Martire showed that "any warp drive requires propulsion" ([CQG 2021](https://arxiv.org/abs/2102.06824)). Fuchs et al. built a **2.365 Jupiter-mass shell** with compactness 0.333 that carries a 0.02c shift with no energy-condition violation ([Fuchs et al. 2024](https://arxiv.org/abs/2405.02709)). They analysed it with Warp Factory, which is **MIT-licensed MATLAB** ([GitHub](https://github.com/NerdsWithAttitudes/WarpFactory); [Helmerich et al. 2024](https://arxiv.org/abs/2404.03095)). A collapsing warp bubble would emit a short, high-frequency GW burst ([Clough, Dietrich & Khan 2024](https://astro.theoj.org/article/121868-what-no-one-has-seen-before-gravitational-waveforms-from-warp-drive-collapse)) **[mainstream]**.

Pfenning and Ford's quantum-inequality bound reproduces to about 10%: 6.9×10⁶² v kg for a 100 m bubble against the published 6.2×10⁶². Their "quarter of a solar mass" for a 1 m wall recomputes as 0.56 M_⊙ ([Pfenning & Ford](https://arxiv.org/abs/gr-qc/9702026)) *(check)*. None of this touches measurable gravity except through planet-scale masses.

## One table ranks every claim against general relativity

| Claim (source) | Published | Reproduces from own formula? | Mainstream value | Gap | Label |
|---|---|---|---|---|---|
| Li–Baker integration time (PRD 2009) | Δt_min = 3×10³ s at h = 10⁻³⁰ | yes, 3.12×10³ s (4%) | homodyne limit 3.9×10¹⁷–3.1×10¹⁸ s | **1.3×10¹⁴–10¹⁵** in time; 7 orders in strain | refuted (JASON 2008) |
| Li–Baker "signal" at 3.5 cm | 3.54×10² s⁻¹ | yes, as 2√(N₀Ṅ_GW) = 3.565×10² | it is the homodyne cross term | — | refuted |
| Li/Baker GW flux | c³ω²A²/8πG | — | c³ω²h²/32πG | 4× in flux, 2× in strain | refuted (convention) |
| Pure inverse-Gertsenshtein (PRD 2009) | 1.4×10⁻³³ | yes, 1.5×10⁻³³ | 8.3×10⁻³⁴ | 1.8× | mainstream |
| Baker jerk constant | 1.76×10⁻⁵² | yes, 1.764×10⁻⁵² | 32G/5c⁵ | none | mainstream |
| Baker coherent array (2012) | 2.64×10⁻⁷ W | yes | incoherent 2.8×10⁻²⁰ W | up to N = 1.6×10¹³ | contested |
| Baker needle beam | 6×10⁵ W/m² at 1 m | yes | diffraction-limited 1.3×10⁻⁵ W/m² | **3×10¹¹** in solid angle | contested (our check refutes) |
| Baker laser eq. 1 | — | sides differ by 144 | — | 144× | refuted (internal) |
| Baker PSR 1913+16 | 1.01×10²⁵ W | no (jerk form 2.06×10²⁵) | 7.75×10²⁴ W | 1.3× | refuted (calibration) |
| Pais S_max (craft) | ~10²⁴ W/m² | yes, 1.77×10²⁴ | field limit 10⁹ V/m | **5.6×10⁶** in field, 3×10¹³ in S | refuted (Navy null) |
| Pais HFGW power | 10¹⁰ W | yes (Gaussian formula, SI inputs) | needs B = 10¹⁰ T, 4×10²⁵ J/m³ | unphysical input | refuted |
| Podkletnov 1997 disc | −2.1% | no formula | GR +2.7×10⁻⁹ | **7.7×10⁶**, opposite sign; rotation term 10²⁰ | refuted for tested configs |
| Podkletnov height | same 2% at 3 m | — | shadow geometry: 0.007% | ~300× | refuted (geometry) |
| Impulse generator | 0.5 m/s kick | yes, pure geometry | GR 5×10⁻³⁰ m/s | **10²⁹**; p/(E/c) = 1.2×10⁹ | contested, unsubstantiated |
| "64c" impulse speed | 63 ns over 1211 m | yes, 63.1 ns | Δv/c < 7×10⁻¹⁶ | 9×10¹⁶ | contested, unsubstantiated |
| Tajmar 2006 ring | k ≈ 2.6×10⁻⁵ | yes, coefficient to 0.04% | k_GR = 7.2×10⁻²⁷ | **3.65×10²¹** | refuted (Tajmar 2011) |
| Li–Torr PRD 1991 | — | 4.2×10⁻³² s⁻¹/T | same | — | mainstream, negligible |
| Torr–Li enhancement | "detectable" | — | Harris: "many orders" too large | not retrieved | refuted (theory) |
| Ummarino–Gallerati | ps–ns transients | yes (with a 4π η change) | static weight change: 0 | — | contested |
| Mach-effect thrust | ~100 µN | own theory gives 7.0 nN | Dresden < 5 nN | 1.4×10⁴ | refuted as tested |
| EmDrive (Eagleworks) | 1.2 mN/kW | — | photon limit 3.34 µN/kW; Dresden ≤ 30 nN | 360× photon limit | refuted |
| DIRD 13 eq. 2.4 | R = 3.37×10²⁴ m | no: 8.6×10⁻¹³ m | — | **37 orders** (dimensional) | refuted (arithmetic) |
| DIRD 21 Shannon | 1.9×10⁶ bps | no: 2.82×10⁶ | — | ln vs log₂ | refuted (arithmetic) |
| DIRD 18 squeezing | −2.23×10⁻²² J/m³ | 13% off | Ford: −1.8×10⁻⁷⁰ | 48 orders | contested |
| DIRD 24 zero-point density | 10¹¹³ J/m³ | between conventions | observed 5.7×10⁻¹⁰ | 10¹²¹–10¹²⁴ | mainstream (Λ problem); devices refuted |
| Puthoff polarisable vacuum | 1.75″, 43″/cy | yes | GR equal at 1PN; 4π vs 15π/4 at 2PN | 0.73 µas | contested (alternative theory) |
| Pfenning–Ford warp bound | 6.2×10⁶² v kg | 6.9×10⁶² (10%) | — | — | mainstream |

## Conclusion

The arithmetic is never where these claims fail. Each one fails at a single step that a calculator can encode as a guard:

| Claim | Failing step | Guard |
|---|---|---|
| Li–Baker detector | noise model | homodyne quantum limit |
| Baker generators | coherence and beam width | diffraction |
| Pais patents | material limits | vacuum breakdown |
| Podkletnov shield | sign, magnitude and geometry | GR pull and shadow fraction |
| Tajmar ring | ad hoc graviphoton mass | GR frame dragging |
| Ummarino–Gallerati | screening of gravity | equivalence principle |
| DIRD 13 | dimensional slip | unit check |

That suggests a design rule for the crates. Never ship a function that only reproduces a claim. Ship it beside the bound it violates, and make the test assert the ratio. The O(1) convention factors (a 4× flux, a ½ in the field, a 4π in Gaussian units, a 4π change in η, two London depths) explain every case where a reproduction misses by less than a factor of 10. They must be explicit parameters for the reproductions to pass, but they are irrelevant beside gaps of 10¹⁴–10²⁹.

The second insight is historical, and it points back at this project. Seventy years of "gravity control" left nothing behind but ion wind, thermal drift and vibration artefacts. Meanwhile, the gravity technology governments actually classified was gravity mapping itself: the Trident gradiometer and Geosat's marine gravity, declassified in the 1990s, and the ancestors of the Falcon instrument and the satellite fields these crates process. The honest frontier is the same one the mainstream notes found. Frame dragging is measured by LARES-2 and sits 16 orders out of laboratory reach. The polarisable-vacuum model differs from GR by 0.73 µas. A gradiometer needs 1 ppb accelerometer matching. Each is real, small and checkable. The events that could change any verdict here are cheap to recognise: a full-condition Podkletnov replication (0.05% is 490 mGal), a Li-group response to JASON, or a peer-reviewed number from Harris or Woods. None has appeared since JASON reported in 2008.

## What goes into the crates

> **As built (2026-10-24).** The labelled calculators went into a new crate, first published as
> `gravity-claims` and renamed [`gravitics`](../crates/gravitics), not a `gwsens::claims` module, because
> propulsion, superconductors and warp drives are not gravitational-wave sensitivity. It has four
> modules: `hfgw`, `superconductors`, `propulsion` and `spacetime`. Each test asserts the
> published value and its gap to general relativity or the replications. The Gertsenshtein
> conversion and the standard wave flux are reused from `gwsens`, and Lense–Thirring stays in
> `geo5d`. The tables below are the full proposal. Rows not yet built remain candidates:
> Peters–Mathews power, gradiometer matching, the Baker diffraction check, Ummarino–Gallerati and
> the DIRD document checks.

All calculators are written from the cited formulas, not translated from GPL code, so the MIT licence holds. Test the formulas to 10⁻⁶ relative error against the scripts' values. Test published comparisons to their stated rounding (4% for Li–Baker, which needs the papers' own k = 100 m⁻¹ and hν = 3.3×10⁻²⁴ J as inputs). DIRD reproductions use the DIRDs' own constants where they state them: G = 6.673×10⁻¹¹, c = 3×10⁸, ħ = 1.055×10⁻³⁴, as `dirds_checks.py` does. Several entries need no new function, only a test that combines existing ones.

### Mainstream additions to existing modules

| Function (module) | Formula | Test vector (from the notes) | Label |
|---|---|---|---|
| `homodyne_cross_rate(n_ref, n_sig)` (`gwsens::hf`) | 2√(N₀·Nₛ) | N₀ = 1.24×10²², Nₛ = 2.562×10⁻¹⁸ → 3.565×10² s⁻¹ (published 3.54×10²) | mainstream |
| `homodyne_min_time(n_sig)` (`gwsens::hf`) | t = 1/Nₛ (JASON eq. 3-34; cross term over shot noise at SNR 2) | Nₛ = 2.562×10⁻¹⁸ → 3.903×10¹⁷ s; Nₛ = 3.203×10⁻¹⁹ → 3.122×10¹⁸ s | mainstream |
| test only (`gwsens::hf`) | `coherent_photon_rate` is the standard ½hBωL, E²/2Z₀ rate; Li's convention is exactly 8× | (1, 0.01 m², 3 T, 6 m, 5 GHz, 10⁻³⁰) → 3.203×10⁻¹⁹ s⁻¹ standard, 2.562×10⁻¹⁸ Li; strain limit at 3×10³ s: 10⁻³⁰·√(3.903×10¹⁷/3×10³) = 1.14×10⁻²³ | mainstream |
| `em_generator_power(u, b, l)` (`gwsens::hf`) | `gertsenshtein_probability(b, l)` · U·c/L = 4πGB²UL/(μ₀c³) | 0.2 J, 10 T, 0.3 m → 1.486×10⁻²⁸ W; ×4π = 1.867×10⁻²⁷ W, JASON's Gaussian-units value (published 2×10⁻²⁷) | mainstream |
| test only (`gwsens::hf`) | Pais eq. 4 in SI is `gertsenshtein_probability(1e10, 1.0)` · 10³³ W | 8.263×10¹⁵ W; as-written Gaussian misuse 1.038×10¹⁰ W; B²/2μ₀ = 3.979×10²⁵ J/m³ | refuted input |
| `strain_from_flux(flux, f)` (`gwsens::lab`) | h = √(32πG·F/(c³ω²)), the inverse of `wave_flux` | 3.53 W/m², f = c/0.061 m → 9.601×10⁻²⁸ | mainstream |
| `linear_pair_power(r, f0, omega)` (`gwsens::lab`) | (16G/15c⁵)r²F₀²ω² | 0.25 m, 2 N, 2π·2.5 GHz → 1.813×10⁻³³ W; ×N² (N = 1.55×10¹³) = 4.357×10⁻⁷ W; ×N = 2.811×10⁻²⁰ W | mainstream |
| `diffraction_half_angle(lambda, d)` (`gwsens::lab`) | 1.22λ/D | 0.06 m, 0.7 m → 0.1046 rad | mainstream |
| `ring_frame_drag(m, r_out, r_in)` (`gwsens::lab`) | k = Ω_drag/ω = 4GM/[c²(R_o+R_i)] | 0.3489 kg, 0.075, 0.069 m → 7.198×10⁻²⁷; × r·ω̇ at 0.036 m, 1500 rad/s² = 3.887×10⁻²⁵ m/s² | mainstream |
| `london_moment(omega)` (`gwsens::lab`) | B = 2m_eω/e | 1.137×10⁻¹¹ T per rad/s; 650 rad/s → 7.391×10⁻⁹ T; GR shift × `ring_frame_drag` = 5.32×10⁻³⁵ T | mainstream |
| `annulus_axis_gravity(sigma, r_in, r_out, z)` (`gwsens::lab`) | 2πGσ[z/√(z²+R_i²) − z/√(z²+R_o²)] | σ = 63.8 kg/m², 0.040, 0.1375 m: Δg/g₀ = 9.579×10⁻¹⁰ at 0.025 m, 2.619×10⁻¹² at 3 m; sheet bound 2πGσ/g₀ = 2.728×10⁻⁹ | mainstream |
| `peters_mathews_power(m1, m2, a, e)` (`gwsens::inspiral`) | (32/5)G⁴m₁²m₂²(m₁+m₂)/(c⁵a⁵)·f(e), f(e) = (1 + 73e²/24 + 37e⁴/96)/(1−e²)^(7/2) | 1.4398, 1.3886 M_⊙, 1.95×10⁹ m, e = 0.6171 → 7.749×10²⁴ W | mainstream |
| `gradiometer_matching(gamma_asd, baseline, accel_asd)` (`satgrav::noise`) | Γ·L/a | 10 E/√Hz, 0.1 m, 1 m s⁻²/√Hz → 1.0×10⁻⁹ | mainstream |
| `asd_to_rms(asd, speed, wavelength)` (`satgrav::noise`) | ASD·√(v/λ) | 6 E/√Hz, 55 m/s, 300 m → 2.569 E (Falcon published 2.5) | mainstream |
| registry entry (`gwsens::hf::REGISTRY`) | one `Status::Disputed` row for the Li–Baker claim (h ≈ 10⁻³⁰ at 5 GHz, 1 Hz band) | appears on the landscape plot beside the ceiling | refuted |

The Lense–Thirring node rate stays where the earlier plan put it, in `geo5d`'s relativity module: 2GJ/(c²a³(1−e²)^(3/2)) gives 30.67 mas/yr for LAGEOS. `ring_frame_drag` and Earth's equatorial Ω_LT = GJ/(c²R³) = 1.683×10⁻¹⁴ rad/s give two more fixtures for the same physics.

### A labelled `claims` module

Put these in a new `gwsens::claims` module. Every doc comment opens with its label and names the bound it violates, and every test asserts the published value *and* the ratio to the mainstream function above.

| Function | Formula | Test vector | Label, and the bound it is tested against |
|---|---|---|---|
| `li_baker_dt_min(n_other, n1)` | N_other/N₁² | 2.1×10⁹, 8.2×10² → 3.123×10³ s (published 3.0×10³); 2.1×10¹¹ → 3.123×10⁵ s | refuted; `homodyne_min_time` ratio ≥ 1.30×10¹⁴ |
| `baker_jerk_power(r, df, dt)` | (32G/5c⁵)(2rΔf/Δt)² | constant 1.764×10⁻⁵²; (0.25, 2, 4×10⁻¹⁰) → 1.100×10⁻³³ W; (10, 4×10⁸, 2×10⁻¹⁰) → 2.82×10⁻¹³ W (DIRD 21); 2r = 0.5, Δf = 3.1×10¹³ → 2.643×10⁻⁷ W | mainstream formula, refuted array use; ratio to `linear_pair_power` = 0.607 |
| `baker_strain(flux, nu)` | √(8πG·F/(c³ω²)) = 1.256×10⁻¹⁸√F/ν | (3.53, c/0.061) → 4.801×10⁻²⁸ (published 4.9×10⁻²⁸); (7×10⁻¹¹, 5×10⁹) → 2.14×10⁻³³; (3, 1.2×10¹⁴) → 1.85×10⁻³²; (3.55×10⁻⁹, 5×10⁹) → 1.496×10⁻³² | refuted convention; `strain_from_flux` is exactly 2× |
| `baker_needle_half_angle(n)` | √(1.71/(πN)) | N = 1.55×10¹³ → 1.874×10⁻⁷ rad | contested; (`diffraction_half_angle(0.06, 0.7)`/this)² = 3.11×10¹¹ |
| `pais_smax(f_g, sigma, v)` with `surface_field(sigma)` = σ/ε₀ | f_G·σ²·v/ε₀ | (1, 5×10⁴, 2·3141.6) → 1.774×10²⁴ W/m²; field 5.647×10¹⁵ V/m; σ = ε₀·10⁹ → 5.563×10¹⁰ W/m²; (1, 1, (10⁷)³), i.e. eq. 2 → 1.129×10³² | refuted; field ≥ 5.6×10⁶ × the patent's own 10⁹ V/m |
| `pendulum_impulse(l, dl, m, t)` | Δh = l − √(l²−Δl²); ΔE = mgΔh; v = √(2gΔh); F = v/t | (0.8, 0.142, 0.0185, 10⁻⁴) → 12.70 mm, 2.305×10⁻³ J, 0.499 m/s, 4.99×10³ m/s²; all 14 table rows to rounding | geometry mainstream, claim contested; GR kick 5.157×10⁻³⁰ m/s (ratio 9.7×10²⁸); p/(E/c) = 1.201×10⁹ |
| `shadow_fraction(h)` | 1 − cos³(atan(1/h)) | 0.6464 (h = 1), 0.2845 (2), 0.01481 (10), 0.003143 (21.8) | mainstream geometry; refutes a 2% reading at 3 m |
| `cooper_density_ratio(lambda_l, rho)` | m_e²/(μ₀e²λ_L²ρ) | (27.6 nm, 8570) → 3.94×10⁻⁶; (39 nm) → 1.974×10⁻⁶; m_e²/(μ₀e²g₀) = 2.623×10⁻¹⁸ | contested input to Tajmar–de Matos |
| `coupling_to_k(coupling, r)` | \|coupling\|·g₀/r | (9.64×10⁻⁸ s², 0.036 m) → 2.626×10⁻⁵ | refuted; ÷ `ring_frame_drag` = 3.648×10²¹ |

Four propulsion one-liners are optional. They are not gravity, but they complete the public-record ledger. `ion_drift_thrust(i, d, mu)` = I·d/μ gives 0.16 N for (10⁻³, 0.04, 2.5×10⁻⁴) **[mainstream]**. `photon_thrust(p)` = P/c gives 3.336×10⁻⁶ N per kW and 3.67×10⁻⁸ N at 11 W **[mainstream]**. `shawyer_thrust(p, q, df)` = 2PQD_f/c gives 3.045×10⁻⁵ N at (11 W, 741, 0.56) and 1.597×10⁻³ N at (60 W, 7123, 0.56) **[refuted]**. `mach_effect_thrust(omega, m0, kp, ke, x0, v0, rho0)` = ω⁶m₀K_p²K_e x₀³V₀⁴/(8πGρ₀c²) gives 6.99×10⁻⁹ N at (2π·32.35 kHz, 0.042 kg, 3.2×10⁻¹⁰, 4×10⁻¹¹, 0.019 m, 100 V, 7900 kg/m³) **[contested theory, refuted as tested]**.

### Doc tests only, without a public function

These checks are one or two multiplications each. Keep them as documented assertions, so the claims ledger carries the numbers without growing the API.

| Check | Values to assert | Label |
|---|---|---|
| Marx bank | 20 × 25 nF at 100 kV: 2500 J, 1.25×10⁻⁹ F, 2.5×10⁻³ C, against the paper's 0.1 C | refuted (40×) |
| 64c delay | 1211/(64c) = 6.312×10⁻⁸ s; light time 4.039×10⁻⁶ s | contested |
| Li–Torr | μ_g(m_e/e)/μ₀ = 4.222×10⁻³² s⁻¹ per T; (m_e/e)·10⁻¹⁴ = 5.686×10⁻²⁶ T | mainstream, negligible |
| Ummarino–Gallerati | g* = 1.919×10⁻¹¹ (YBCO, 77 K); η = 4πμ₀σD = 1.263×10⁻², while Bi2212 needs μ₀σD = 3.491×10⁻⁴; recovery ητ = 4.318×10⁻¹² s; Josephson 2m_e gℓ/ħ = 0.05083 s⁻¹ at 300 nm | contested |
| Equilibrium field | Schiff–Barnhill m_e g/e = 5.576×10⁻¹¹ V/m | mainstream |
| DeWitt and penetration depth | 8πε₀Gm_eMv/(ed) = 7.96×10⁻²⁹ A (1 kg, v/d = 300π); λ_g = 7.96×10¹³ m (Nb, 39 nm) | mainstream |
| Modanese | required Λ is 10^19.96 × observed | contested |

### Optional, outside the gravimetry core

Implement these only if the workspace gains a relativity module. `casimir_pressure(d)` = −π²ħc/(240d⁴) gives −1.300×10⁻³ Pa at 1 µm, with energy density −4.334×10⁻⁴ J/m³. The sphere–plate force is −2.723×10⁻¹⁰ N at R = 100 µm and d = 100 nm. The weight of the Casimir energy of 1 cm² plates at 10 nm is 4.730×10⁻²⁴ N **[mainstream]**. `qi_min_energy_density(t0)` = −3ħc/(32π²(ct₀)⁴) gives −3.718×10⁻²⁶ J/m³ at 1 ns; use it as a guard on any negative-energy input **[mainstream]**.

A null-geodesic deflection integrator for isotropic metrics, n(u) = √(B/A), should reproduce four checks, as `dirds_checks.py` does with composite Gauss–Legendre:

- the second-order coefficients 12.566 (exponential metric) and 11.780 (Schwarzschild) to 2×10⁻²;
- 1.7512″ of solar light bending;
- 42.98″/century for Mercury;
- the exponential metric's throat at areal radius 4014 m per solar mass.

The polarisable-vacuum model itself is **[contested]**. The DIRD warp and wormhole one-liners belong in regression tests at most, never in an API:

- Lobo–Visser −3.035×10⁵⁰ J at v = 1;
- v_max = 1.722×10⁻¹⁴;
- the corrected R = 8.611×10⁻¹³ m;
- the Morris–Thorne–Yurtsever spacing 1.281×10⁻¹⁶ m at 1 km;
- Visser's edge tension −1.513×10⁴³ J/m;
- the Morris–Thorne throat tension 4.815×10⁴⁰ Pa at 10 m;
- the corrected Shannon capacity 2.816×10⁶ bps.

### What stays out

Leave out the following as physics:

- any model of a gravity-wave generator, shield, thruster or communicator;
- Obousy's extra-dimension warp (no metric);
- the gravitational-squeezing density (Ford's objection);
- Modanese's local Λ (it predicts no checkable number);
- the Ummarino–Gallerati "screening" (it violates the equivalence principle);
- Puthoff's "craft signature" table (metric identities with no source);
- Lentz's soliton (unspecified form factor and the observer dispute);
- every zero-point-energy device.

A metric → stress-energy → energy-condition evaluator would be the most interesting new Rust code in this area, validated on Schwarzschild, Alcubierre, Morris–Thorne and the Fuchs shell. But it has no data to test against beyond published figures and no link to satellite gravimetry. If it is wanted, it belongs in a separate crate, written clean-room from the papers, with MIT-licensed Warp Factory as the reference.
