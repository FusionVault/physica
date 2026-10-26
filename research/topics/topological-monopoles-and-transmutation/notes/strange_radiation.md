# Strange radiation, Urutskoev's exploding foils, Lochak's leptonic monopole, and the related Russian anomaly literature: claims, numbers, conventional physics, calculators

Scope: the signatures listed in B. W. Greenyer's talk "The Power and Signature of the Topological Monopole" (TeslaTech, 15 Aug 2026; transcript lines 410–446, 676–735, 893–1031, 1453–1528). General LENR material (DOE reviews, Google 2019, electron screening) is covered in `research/topics/coulomb-and-schwinger/notes/fringe_barrier_claims.md` and not repeated here.

Every computed number below is printed and, where a published value exists, asserted by `strange_radiation_checks.py` in this folder (`python3 -I strange_radiation_checks.py`, standard library, runs in well under a second, non-zero exit on failure). Downloaded primary sources are in `dl_strange/` (git-ignored). Status labels: **mainstream**, **contested**, **refuted**, **unverified**.

## Q1. Urutskoev, Liksonov & Tsinoev (2002) and later work: what was reported, with what numbers, and who replicated it

### Takeaway
Urutskoev's group (RECOM, Kurchatov Institute) reported "alien" elements and a depleted ⁴⁸Ti fraction in the residue of titanium foils exploded in water, film tracks they called "strange radiation", and a magnetic-field effect on the tracks. They detected no gamma rays, neutrons or residual radioactivity. The reports come from the group itself and its collaborators; the published element list matches the setup's own materials. One independent-laboratory check (Nantes, with Lochak's foundation) used stainless-steel electrodes, and its "alien" Fe/Cr/Ni follows the steel's composition. Status: **unverified** (transmutation, monopole origin); the absence of radiation is the group's own measured result.

### Cited Findings
- Set-up: a capacitor bank at U ~ 4.8 kV storing W ~ 50 kJ discharged into Ti foils welded to Ti electrodes, inside a polyethylene chamber filled with distilled water; current pulse ~0.15 ms; glow above the cover ~5 ms; light energy of the "ball-like plasma formation" estimated at ~1 kJ; "more than 100 tests" for the glow dynamics — [Urutskoev, Liksonov & Tsinoev, AFLB 27(4) 701 (2002)](https://fondationlouisdebroglie.org/AFLB-274/aflb274p701.pdf)
- Spectra: >1000 lines per shot, no O or N lines (only weak traces), most abundant Ti, Fe, Cu, Zn, Cr, Ni, Ca, Na; Cu and Zn attributed by the authors to sliding discharges on cables — [Urutskoev 2002](https://fondationlouisdebroglie.org/AFLB-274/aflb274p701.pdf)
- Starting foil 99.71643 % Ti with listed impurities (Fe 0.10613 %, Ni 0.04193 %, Ca 0.03399 %, Cr 0.00844 % …), natural Ti isotope ratio; residue mass ~0.5 g, visibly inhomogeneous; test 226: Ti 92 % of atoms in the products; "alien" elements in "more than 200 experiments"; mean Ti transformation 4 % over 24 tests (N 169–240); Fig. 10 plots up to ~40 % decrease of ⁴⁸Ti weight fraction for 21 experiments, drawn on the assumption that all lost Ti is ⁴⁸Ti — [Urutskoev 2002](https://fondationlouisdebroglie.org/AFLB-274/aflb274p701.pdf)
- "Regularities" claimed: transformation mainly of even-even isotopes; binding-energy difference between original and product mix ΔE < 0.1 MeV/atom "clearly due to mass-spectrometric measurements errors"; products in the ground state, "no appreciable radioactivity"; 10¹⁹–10²⁰ transformations per shot; no significant X-ray or gamma flux — [Urutskoev 2002](https://fondationlouisdebroglie.org/AFLB-274/aflb274p701.pdf)
- Particle detection: plastic scintillators at 0.4 m and 0.8 m gave ~100 ns pulses; the delay between them gave a "radiation propagation velocity" V ~ 20–40 m/s — [Urutskoev 2002](https://fondationlouisdebroglie.org/AFLB-274/aflb274p701.pdf)
- Film/emulsion: RF-ZMP fluorographic film, RM-1MD medical film, type R nuclear plates (100 µm), wrapped in black paper, 20 cm–4 m from the explosion; "caterpillar" tracks 1–3 mm long, ~20 µm wide, energy "E ~ 700 MeV assuming Coulomb braking" from the darkened area; tracks lie in the emulsion plane (entry/exit depths differ by <10–15 µm); a source velocity β = 10⁻³ inferred from track length and pulse duration; with H ~ 20 G at the explosion site, "comet" tracks 300–1300 µm, six in 4 cm², "E ~ 1 GeV"; three stacked films near a 1.2 kG Sm-Co magnet: coincident darkenings, ~700 MeV; some marks "very similar to scratches or ink spots", which the authors could not unconditionally exclude as artefacts — [Urutskoev 2002](https://fondationlouisdebroglie.org/AFLB-274/aflb274p701.pdf)
- Urutskoev's 2022 retrospective: all nuclear detectors at the institute were "silent"; the emulsion tracks came first; a magnetic field changed "caterpillar" to "comet" tracks; the energy needed to make one Dirac monopole pair per transformation exceeded the capacitor bank, which they called a "catastrophic contradiction with the conservation law"; >400 experiments before publication; "such a device [a monopole counter] has not been developed to this day" — [Urutskoev, Rybakov & Samsonenko, AFLB 47(1) 153 (2022)](https://fondationlouisdebroglie.org/AFLB-471/aflb471m969.pdf)
- A search summary of Urutskoev's later work says that some shots with the same stored energy were quieter and showed no plasma "ball" and no isotopic distortion (I did not open that paper; treat as unverified) — [search result, not opened](https://www.researchgate.net/scientific-contributions/LI-Urutskoev-34693546)
- Dubna check (Kuznetsov, Mishinsky, Penkov, Arbuzov, Zhemenik, JINR): foil explosions in water, H₂O₂ and glycerin with Ti, Ta, Ta+Pb; reports increased Na, Al, K, Ca, Mn, Fe, Ni and others after Ta-foil explosions; "no gamma or beta activity was detected in any sample"; argues the reaction domain must be atomic-sized. (The PDF uses a shifted font; numbers in it did not extract, so none are quoted here.) — [Kuznetsov et al., AFLB 28, 173 (2003)](https://fondationlouisdebroglie.org/AFLB-282/aflb282p173.pdf). Lochak cites it as vol. 24 — [Lochak 2007, ref. 52](https://arxiv.org/abs/0801.2752); the PDF path and Priem et al. give vol. 28 — [Priem 2008, ref. 6](https://fondationlouisdebroglie.org/AFLB-331/aflb331m634.pdf)
- Nantes check (Priem, Racineux with Lochak, Daviau, Fargue, Karatchentcheff, Lehn): Magneform 360 µF bank, 3 kJ at 4.2 kV; 22 ml chamber; 0.3 g T40 Ti fuse (99.81 at.% Ti, 0.13 % Fe); **stainless-steel electrodes (18-8: Fe 72 %, Cr 20 %, Ni 8 %)**; EDX of residues: Ti fell and Fe rose in every shot (e.g. shot 2: Fe 5.16, Cr 1.1, Ni 0.52, Cu 0.41, Zn 0.55 %; shot 3: Fe 6.6, Cr 1.3, Ni 0.67 %); authors: "results in good agreement" with Urutskoev but "it is indispensable to repeat the tests using another material… to determine the origin of the iron" — [Priem et al., AFLB 33, 129 (2008)](https://fondationlouisdebroglie.org/AFLB-331/aflb331m634.pdf)
- Later Nantes/Paris work looked for deuterium instead: gas D/H 150 ppm reference vs 168, 217, 226 ppm in three shots (CEA Saclay, ±1 ppm); 13 water samples did not deviate from the evaporation ("meteoric") line; an earlier claimed water enrichment 157 → 194 ppm is now called "doubtful", the same samples later read 390 ppm (and the reference 200 ppm), "most probably a pollution" — [de Lacheze-Murel et al., AFLB 47, 123 (2022)](https://fondationlouisdebroglie.org/AFLB-471/aflb471m990.pdf)
- Ivoilov (Kazan) low-energy version: carbon electrodes in water, ≤40 A at 80 V, ~150 Oe Helmholtz field; X-ray film 10–15 cm away; 0–10 tracks per 10 min over 35 cm², "sometimes tracks are absent during some days"; Fe and Ni films "totally absorb", glass, Si, Ge "reflect"; ~5 % of tracks come as centre-symmetric pairs on opposite film sides; **control films left 10 min in a 20 kOe magnet, never near the arc, showed "the same tracks"** (called cosmic "control background"); a Sr-90 (1.2×10⁶ Bq) or Cs-137 (2.0×10⁸ Bq) source in the field doubled the track count over that background; >20 such experiments — [Ivoilov, AFLB 31, 115 (2006)](https://fondationlouisdebroglie.org/AFLB-311/aflb311m484.pdf)
- Lochak's summary adds: tracks are "hairless" (no delta electrons); track width ~30 µm at ~0.5 m and ~5 µm at 2 m; ammonium nitrate decomposed in Al containers but not steel ones; Chernobyl ²³⁵U "enriched up to 27 %" offered as a monopole effect; biological effects on several hundred mice (Pryakhin) — [Lochak, Z. Naturforsch. 62a, 231 (2007), DOI 10.1515/zna-2007-5-602](https://arxiv.org/abs/0801.2752)
- Savvatimova & Rodionov reported "strange" nuclear tracks from glow-discharge work at RCCNT&BL-13 (2005); the conference summary says the tracks' nature "might probably disclose the mechanism of Cold Nuclear Transmutation" — [RCCNT&BL-13 summary, Bazhutov](https://newenergytimes.com/v2/conferences/2005/RCCNTBL13/RCCNTBL-13-Summary.shtml)

### Inferences
- Energy budget (computed): 10¹⁹–10²⁰ transformations × the authors' own 0.1 MeV/atom error bound = 160–1602 kJ, versus a 50 kJ bank. Their bound does not exclude a net nuclear energy change several times the stored energy, so "ΔE ≈ 0" is not established by their data. 10¹⁹–10²⁰ ⁴⁸Ti atoms is 0.8–8 mg.
- Isotope bookkeeping (computed): if 4 % of Ti is "burned" and all of it is ⁴⁸Ti (natural 73.72 %, IAEA), the remaining ⁴⁸Ti fraction is 72.62 %. A ⁴⁸Ti fraction of 65 % (quoted secondhand) would need 24.9 % burnout, far above the 4 % mean.
- Contamination fingerprint (computed from Priem Table 3): shot 2's Fe implies 7.2 % electrode steel; the observed Cr and Ni are 0.77 and 0.91 of what that steel fraction would bring; for shot 3, 9.2 % steel and 0.71, 0.91. Within EDX accuracy, the "alien" Fe, Cr and Ni follow 18-8 steel. In Urutskoev's 2002 list, Fe, Cr, Ni, Cu, Zn, Ca and Na are also common materials in the setup, its cables, water and handling.
- The two velocity estimates in the 2002 paper differ by 7495–14990× (computed: β = 10⁻³ is 2.998×10⁵ m/s, against the 20–40 m/s time-of-flight figure). A Dirac-charged particle in Earth's field also gains 1.03 MeV per metre (computed), so a 20–40 m/s magnetic particle would need very strong drag or a very large mass.
- Ivoilov's controls are central: identical tracks appeared on films that were only stored in a 20 kOe magnet. A film-handling or environmental mechanism (static discharge, pressure or chemical marks, or cosmic rays) therefore makes the tracks without any discharge. A blind protocol (randomised exposed/sham films, scored by someone who does not know which is which) is the obvious missing test.

### Gaps
- No independent group with non-steel, non-Ti-contaminating electrodes and blind isotope analysis was found; no Western replication of the Ti isotope shift was found.
- I did not find a published critique by Volkovich, Stepanov or others. Volkovich appears only in the 2002 acknowledgements. I did not find Skvortsov's work.
- The mass-spectrometer type, isobar corrections (⁵⁰Cr 4.345 %, ⁵⁰V 0.250 %, ⁴⁸Ca 0.187 %, IAEA) and per-shot uncertainties for the Ti isotope data are not given in the 2002 paper.

## Q2. The Fe-57 Mössbauer ("NMR") claim, and what a monopole near ⁵⁷Fe nuclei would do

### Takeaway
The claim is a change of about ±0.25 kOe (±8×10⁻⁴) in the mean hyperfine field of ⁵⁷Fe foils placed near the N or S pole of a magnet. That is the size of the shift a temperature difference of 6–7 K would cause. The 2002 report does not mention temperature control. A Dirac pole's own field is enormous close by (5.4 kT at the Fe nearest-neighbour distance, computed). But it is radial, so it averages to zero along the in-plane magnetisation and cannot give a sign-selective shift of the mean field in a simple geometry. Status: **unverified**; the size of the effect is at the level of known systematics.

### Cited Findings
- 2002: three ⁵⁷Fe foils, two near opposite poles of a ~1 kG magnet 70 cm from the explosion, one standard; H = 330.42 kG (standard), 330.66 kG (N, +0.24), 330.13 kG (S, −0.29), error 0.012 kG; "relative change ~8·10⁻⁴"; added line broadening, e.g. Γ₁ 0.334 → 0.363/0.366 mm/s (±0.003); no quadrupole shift — [Urutskoev 2002, pp. 722–723](https://fondationlouisdebroglie.org/AFLB-274/aflb274p701.pdf)
- 2004 (Kazan): 20 µm foils enriched to 94 % ⁵⁷Fe, 15×15 mm; conversion-electron detection (probe depth ~0.2 µm, vs ~10 µm for γ/X-ray detection); 20 mCi ⁵⁷Co(Cr) source; detector held at 35 ± 0.5 °C, vibrator at 32 ± 0.5 °C; 3 h spectra, σ(H_eff) = ±50 Oe; H_eff ≈ 330 kOe, line ratio ≈ 3:4:1:1:4:3 (in-plane magnetisation); 14 experiments; first effects vanished when the samples were re-measured 2 days later; December 2002 series: sample No. 6(N) decayed with τ = 16 ± 5 h from +280 ± 70 Oe, No. 5(S) only its first point deviated, by 460 ± 70 Oe, sample 7(L) and the control were flat for 80 h; the authors quote a hyperfine-field temperature coefficient of −40 Oe/K and an isomer-shift coefficient of +0.65 µm/s/K, citing Preston et al. 1962 — [Ivoilov & Urutskoev, AFLB 29 HS3, 1177 (2004)](https://fondationlouisdebroglie.org/AFLB-297/aflb297m331.pdf)
- Lochak describes the result as a shift "in the opposite direction" for the other pole, and says the spectrum returns to normal if the sample is left for three days — [Lochak 2007](https://arxiv.org/abs/0801.2752)
- α-Fe at room temperature: H_n = 330 ± 3 kOe; μ₁/μ₀ = −1.715 ± 0.004; Debye temperature 400 ± 30 K — [Preston, Hanna & Heberle, Phys. Rev. 128, 2207 (1962), DOI 10.1103/PhysRev.128.2207](https://journals.aps.org/pr/abstract/10.1103/PhysRev.128.2207)
- ⁵⁷Fe moments (ENSDF): ground 1/2⁻, μ = +0.0902 μN; 14.4129 keV 3/2⁻ state, μ = −0.15531 μN, T½ = 98.3 ns — [IAEA LiveChart API](https://nds.iaea.org/relnsd/v1/data?fields=levels&nuclides=57fe)
- Dirac field: B = Q_M r̂/(4πr²) with the quantisation Q_E Q_M = 2π (natural units) — [PDG 2024, Magnetic Monopoles, Eqs. 95.1–95.2](https://pdg.lbl.gov/2024/reviews/rpp2024-rev-mag-monopole-searches.pdf)
- The idea of using ferromagnetic domains as monopole traps comes from Martemyanov & Khakimov, ZhETF 62, 35 (1972) (as cited by both Kazan and Kurchatov papers; I did not open it) — [Ivoilov & Urutskoev 2004, ref. 2](https://fondationlouisdebroglie.org/AFLB-297/aflb297m331.pdf)

### Inferences
- Sextet (computed from the ENSDF moments at 33.0 T): lines at ±0.832, ±3.072, ±5.313 mm/s, outer splitting 10.625 mm/s; μ₁/μ₀ from ENSDF is −1.722, within 0.4 % of Preston's −1.715. The 2002 shifts move the outer lines by +3.86 µm/s and −4.66 µm/s. Velocity-scale drift of that size has to be controlled by recalibrating against a standard foil in the same run.
- Temperature equivalent (computed at −40 Oe/K): 2002 N +0.24 kG is 6.0 K, S −0.29 kG is 7.2 K; 2004 No. 6(N) 280 Oe is 7.0 K, No. 5(S) 460 Oe is 11.5 K. The 2004 detector was held to ±0.5 K, but the sample's history between irradiation and measurement (transport, an hours-long relaxation, the first point only) is where an uncontrolled variable would act.
- Second-order Doppler (computed, Debye θ = 400 K at 308 K): |dδ/dT| = 0.673 µm/s/K, which matches the 0.65 µm/s/K the Kazan authors used. Their isomer-shift argument is quantitatively sound. But they also report that δ(t) changed "very close [to the] methodical error", so the temperature cross-check is weak.
- Monopole field (computed): 3.29×10¹⁴ T at 1 fm, 3.29×10¹² T at 10 fm, 5.35×10³ T at the Fe nearest-neighbour distance (2.48 Å), 1.18×10⁵ T at the Bohr radius; it falls to the α-Fe hyperfine field (33 T) at 3.16 nm. A nucleus within a few nm of a trapped pole would have its lines moved far out of the measured sextet. It would not take part in a small shift of the main lines.
- Geometry: the pole field is radial. Averaged over the nuclei around one pole, its component along the magnetisation is zero to first order. A uniform sheet of same-sign poles gives a field normal to the foil, perpendicular to the in-plane magnetisation, so |B| rises at second order whatever the sign. A +0.024 T rise needs a 1.26 T normal field, a sheet of 6.09×10¹⁰ Dirac poles/cm² (computed). A sign-dependent ±0.025 T mean shift is therefore not what simple trapped-monopole geometry predicts; a sign-independent rise plus broadening would be.
- Interaction energy of a pole with the ⁵⁷Fe ground-state moment at 10 fm: |μ·B| = 9.4 keV (computed), the scale of the "monopole bound to a paramagnetic nucleus" picture used in the talk.

### Gaps
- The Russian original (Ivoilov & Urutskoev, Prikladnaya Fizika 2004 no. 5) and any spectra with error bars from other groups were not found. I found no independent Mössbauer replication.
- The −40 Oe/K slope was taken from the Kazan paper's citation of Preston 1962; the Preston abstract I opened does not state it.

## Q3. Lochak's light leptonic monopole: equations, predicted properties, relation to Dirac quantisation

### Takeaway
Lochak's monopole is the massless Dirac equation with a second, chiral (γ₅) gauge coupling to an axial pseudo-potential. It describes a massless, spin-½, weakly interacting particle with magnetic charge, a "magnetically excited neutrino". Its charge obeys a Dirac-type relation eg/ħc = m′, a half-integer projection. With n = 1 the pole strength is the usual g_D = 68.5 e (Gaussian). The theory does not predict a mass or a lifetime. It is not part of the Standard Model; Stumpf proposed an extension. No detector has measured a light monopole. Status: **contested** (a self-consistent proposal, outside mainstream particle physics, with no confirmed experimental signature).

### Cited Findings
- Equation: γ^μ(∂_μ − (g/ħc) γ₅ B_μ)Ψ = 0 (Lochak's Eq. 6.23), where B_μ is a pseudo-potential with E = curl **B** and H = ∇W + (1/c)∂**B**/∂t (Eqs. 4.7–4.8, after de Broglie and Cabibbo–Ferrari); chiral gauge invariance holds only for m₀ = 0, so "conservation of magnetism is weaker than conservation of electricity"; the conserved magnetic current is K_μ = i g Ψ̄γ_μγ₅Ψ (Eq. 7.1) — [Lochak, Z. Naturforsch. 62a, 231–246 (2007)](https://arxiv.org/abs/0801.2752) (DOI 10.1515/zna-2007-5-602, Crossref-checked)
- Weyl split into two uncoupled two-component equations (Eq. 8.4) for left and right monopoles of opposite chirality, which "reduces to a pair of two-component neutrino equations if g = 0"; charge conjugation does not reverse g, so "we cannot create or annihilate pairs of monopoles with charges g and −g", which Lochak uses to argue against vacuum-polarisation catastrophes and for monopoles being "abundant in nature" — [Lochak 2007](https://arxiv.org/abs/0801.2752)
- Dirac relation: eg/ħc = n/2 (Eq. 5.2); in his derivation D = eg/ħc = m′ ∈ {−j, …, j} (Eq. 11.10), with m′ the projection of total angular momentum on the charge–pole axis — [Lochak 2007](https://arxiv.org/abs/0801.2752)
- A nonlinear generalisation gives a mass term F(ρ) and both bradyon and tachyon dispersion relations (Section 12) — [Lochak 2007](https://arxiv.org/abs/0801.2752)
- Predictions stated: monopoles produced by electromagnetic pulses or arcs, nuclear transmutations, shorter beta lifetimes, emission of monopoles instead of neutrinos in a magnetic field — [Lochak 2007, abstract](https://arxiv.org/abs/0801.2752)
- In Gaussian units the elementary pole is 137/2 times the electron charge, 3.29×10⁻⁹ in SI (A·m); one pole moving 0.1 m in 1 G gains 205 keV — [de Lacheze-Murel et al. 2022, footnote 1](https://fondationlouisdebroglie.org/AFLB-471/aflb471m990.pdf)
- Stumpf placed the monopole in an extended Standard Model (Z. Naturforsch. 66a–67a, 2010–2012); Lochak & Stumpf, *The Leptonic Magnetic Monopole* (Elsevier 2015); Urutskoev: the theory "is in no way connected with the Standard Model"; transformations in ~1 eV plasma are collective and "cannot be described by the single-particle equation" — [Urutskoev et al. 2022](https://fondationlouisdebroglie.org/AFLB-471/aflb471m969.pdf)
- Mainstream constraints: no confirmed monopole observation; electron magnetic charge < 10⁻²⁴ Q_D (induction); muon g−2 gives a model-dependent mass limit of 120 GeV; LHC limits up to ~4 TeV for pair production; bulk matter < ~10⁻²⁹ monopoles per nucleon; electroweak theory has no topologically stable monopoles — [PDG 2024](https://pdg.lbl.gov/2024/reviews/rpp2024-rev-mag-monopole-searches.pdf)

### Inferences
- g_D = h/(μ₀e) = 3.291×10⁻⁹ A·m, g_D/(ec) = 1/(2α) = 68.52, (g/e)² = 4695 (all computed, matching the published 3.29×10⁻⁹, 68.5 and ~4,700).
- Collider pair-production limits assume production of g/−g pairs. Lochak's C-symmetry argument says that channel does not exist, so those limits do not directly test his particle. Any light particle with Dirac magnetic charge still couples strongly to photons (α_m = 1/(4α) ≈ 34), so loop effects (muon g−2, PDG) and SQUID induction loops remain direct tests. Inference; not computed here.
- A massless monopole travels at c. That contradicts β ≈ 10⁻³ and 20–40 m/s in the experiments that are offered as evidence for it.

### Gaps
- Lochak's 1985 IJTP paper (24, 1019) and the 2015 book were not opened; their DOIs were not verified.
- No calculation of a light monopole's ionisation, its trapping in Fe, or its effect on beta decay rates was found in the theory papers.

## Q4. Klimov's vortex plasma reactor and Filimonenko: excess energy, transmutation, claimed decay suppression

### Takeaway
Klimov's (JIHT RAS) published JCMNS papers report COP 2–10, new elements (Li, Ca, Fe, Cu, Zn …), a changed Li isotope ratio and soft X-rays. I found **no** published measurement in which a gamma source's decay was suppressed for ~2 days; that claim from the talk remains unsourced. In the 2024 paper the measured COP is 1.62–2.82, and it reaches the quoted range only after an unmeasured loss factor is applied. Filimonenko's 1950s claims: no primary source found. Status: Klimov excess heat and transmutation **unverified**; gamma-suppression claim **unverified** (no source); Filimonenko **unverified** (no source).

### Cited Findings
- 2016 (JCMNS 19, 155): RF plasmoid over water; Li and Ca concentrations up 10²–10³×; Li-7/Li-6 13.192 ± 0.004 → 18.996 ± 0.012; COP 4–10 behind shock waves, 2–4 in swirl flow — [Klimov, JCMNS 19, 155–163 (2016)](https://jcmns.org/article/72385.pdf)
- 2016 (JCMNS 19, 145): COP 2–10 at 1–10 kW; soft X-rays 1–10 keV; "dusty plasma-gas flow behind reactor's nozzle is radioactive at long distance"; ~1 keV/atom specific energy release; transmutation "proved" by optical, ion-mass, X-ray and α spectroscopy — [Klimov et al., JCMNS 19, 145–154 (2016)](https://coldfusioncommunity.net/wp-content/uploads/2018/08/145_JCMNS-Vol19.pdf)
- 2024 (JCMNS 38, 287): COP 0.73–0.80 without relaxation oscillations (no sediment), 1.62–2.82 with them (sediment ~17 mg/s); "taking into account these factors [unmeasured losses] increases the indicated values of the COP by 1.5–2 times"; calls for isotopic analysis as further work — [Klimov, Altunin, Kulikovskii, JCMNS 38, 287–294 (2024)](https://jcmns.org/article/124964.pdf)
- 2020 Zvenigorod abstract (Belov, Evstigneev, Klimov, Ryabkov): cold neutron flux ~10⁷ s⁻¹, decaying 5–30 min after switch-off (from a search summary; not opened) — [Zvenigorod 2020 abstract](https://www.fpl.gpi.ru/Zvenigorod/XLVII/Cm/en/KH-Klimov_e.docx)

### Inferences
- The 2024 numbers (computed): measured COP 1.62–2.82 becomes 2.43–5.64 only after the authors' estimated ×1.5–2 for unmeasured losses. The excess-heat claim therefore depends on an estimate, not on a closed calorimetric balance.
- Li-7/Li-6 rose ×1.440 (computed). Lithium is a common trace contaminant and fractionates in electrochemical and plasma processes, so isotope changes alone do not discriminate between transmutation and fractionation.

### Gaps
- No primary source for "gamma generator… suppressed its observable decay rate for ~2 days" (talk, lines 692–700) was found in Klimov's JCMNS papers or in searches (English and Russian).
- No primary source for Ivan Filimonenko's 1950s "warm fusion" or radiation suppression was found (English and Russian searches).

## Q5. Bogdanovich et al. (2019), Nesterovich, the "Moscow Nuclear Physics Institute", and the 2008 conference ("superconductivity of the Dirac monopole", Vishnevsky)

### Takeaway
**Confirmed:** Bogdanovich, Volkov, Len' and Nesterovich, "Video Recording of Long-Lived Plasmoids near Objects Exposed to Remote and Direct Effects of High-Current Pinch Discharges", Technical Physics 64(4), 465–469 (2019). Search results and MEPhI's repository place the authors at **MEPhI** (National Research Nuclear University, "Moscow Engineering Physics Institute"). So the talk's "Moscow Nuclear Physics Institute" is most likely MEPhI; the affiliation was not seen in the paper itself. **Not confirmed:** I found no "Vishnevsky (2008)" or "superconductivity of the Dirac monopole" contribution. The program and abstracts book of the 2008 conference (RCCNT&BL-15, Dagomys, 1–8 Oct 2008) contains no "monopole", "Dirac" or "Vishnevsky" (Вишнев…) in its text layer. The 2008 ball-lightning theory there is A. N. Vlasov's magnetohydrodynamic current-ring model. The audit's attribution should be marked unverified until the user supplies the citation from the slide. Status: Bogdanovich plasmoid lifetimes **unverified** (full text not read).

### Cited Findings
- Crossref: Technical Physics 64(4), 465–469, April 2019, authors Bogdanovich, Volkov, Len', Nesterovich, DOI 10.1134/S1063784219040066; Russian original Zh. Tekh. Fiz. 89(4) 507–511, DOI 10.21883/JTF.2019.04.47304.257-18, Ioffe Institute. Abstract: plasmoids formed during periodic discharge in a flow of liquid (water), during and after the discharge; "long-lived luminous formations not only in the air, but also on the surface and above the surface of nearby objects"; their properties have "features not only of macro-, but also micro-objects, such as Dirac's monopole" — [Crossref record](https://api.crossref.org/works/10.21883/JTF.2019.04.47304.257-18)
- MEPhI affiliation and funding (Russian Science Foundation 17-72-20228), surname also spelled "Nestorovich" — [Math-Net.Ru entry (search summary)](https://www.mathnet.ru/eng/jtf5634); [MEPhI open repository (search summary)](https://openrepository.mephi.ru/entities/publication/8129ceef-e7d5-4d4c-a233-07733cd4707e)
- A related 2017 paper by Bogdanovich, Nestorovich et al. models pinch plasmoids as unclosed wormholes using Einstein–Maxwell solutions (Russian Physics Journal; search summary only) — [bakhtiniada.ru](https://bakhtiniada.ru/1064-8887/article/view/238927)
- RCCNT&BL-15: organised by the Russian Physical Society, Nuclear Society of Russia, Mendeleev Chemical Society, MSU Physics, RUDN and the RAS Committee on Ball Lightning; the program committee includes B. U. Rodionov (MEPhI); ball-lightning theory talks include A. N. Vlasov (Ryazan) on an MHD current-ring model with claimed plasma conductivity ~5×10¹⁸ Ω⁻¹m⁻¹ — [RCCNT&BL-15 Program & Abstracts](https://newenergytimes.com/v2/conferences/2008/RCCNTBL15/Program_Abstracts.pdf)
- The talk's "Afanasiev and Dubovik" toroidal-moment / Aharonov–Bohm references and the "synchrotron shadow" image were not traced.

### Inferences
- "Two-day-living plasmoid on metal" (talk line ~1490) may come from the Bogdanovich video paper. The lifetimes in that paper could not be read: the Ioffe and MEPhI servers refused the connection.
- Rodionov (MEPhI) co-reported "strange tracks" with Savvatimova (Q1). The MEPhI link in the talk may run through him rather than through a "Vishnevsky".

### Gaps
- Full text of Bogdanovich 2019 (plasmoid lifetimes, sizes, the "orange-shell" spheres with axial holes attributed to Nesterovich): not obtained.
- The RCCNT&BL full proceedings (as opposed to the abstracts book) were not found; a Vishnevsky paper could be there or at another 2008 meeting.

## Q6. Beta-decay and EC rate changes: the claims, the refutations, the real effects, and what a 2-day suppression would mean

### Takeaway
Real, measured decay-rate changes come from changing the atomic electrons: bound-state beta decay of bare ¹⁸⁷Re (×1.28×10⁹ faster) and ¹⁶³Dy, and electron-capture changes with chemistry (⁷Be 0.83 %). Claims of seasonal, solar-linked modulation at the 0.1 % level are refuted by multi-laboratory data, which limit annual amplitudes to 6×10⁻⁶–8×10⁻⁵. A total 2-day suppression of a ⁶⁰Co source would leave a permanent +7.2×10⁻⁴ activity step, 9–120× above those limits, so routine metrology would see it. Status: ionisation and chemistry effects **mainstream**; solar modulation **refuted**; monopole-induced changes **unverified**.

### Cited Findings
- Bare ¹⁸⁷Re: T½ = 32.9 ± 2.0 yr versus 42 Gyr for neutral atoms; log ft 7.87 ± 0.03 — [Bosch et al., PRL 77, 5190 (1996), DOI 10.1103/PhysRevLett.77.5190](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.77.5190)
- First bound-state beta decay (bare ¹⁶³Dy⁶⁶⁺) — [Jung et al., PRL 69, 2164 (1992), DOI 10.1103/PhysRevLett.69.2164](https://api.crossref.org/works/10.1103/PhysRevLett.69.2164) (title and DOI verified; numbers not read this session)
- ⁷Be: T½ = 52.68 ± 0.05 d inside C₆₀ versus 53.12 ± 0.05 d in Be metal, a 0.83 % difference — [Ohtsuki et al., PRL 93, 112501 (2004), DOI 10.1103/PhysRevLett.93.112501](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.93.112501)
- Claim: periodic fluctuations in ³²Si (BNL) and ²²⁶Ra (PTB) decay rates correlated with the Earth–Sun distance — [Jenkins et al., Astropart. Phys. 32, 42 (2009), DOI 10.1016/j.astropartphys.2009.05.004](https://arxiv.org/abs/0808.3283)
- Refutation: 14 laboratories, sources followed for 200 days to four decades; annual oscillation amplitude < 0.0006 %–0.008 %; deviations "attributable to instabilities in the instrumentation"; the Jenkins reading of the PTB ²²⁶Ra data was ~0.15 %, Pommé's own fit 0.083(2) %, attributed to instrumentation — [Pommé et al., Phys. Lett. B 761, 281 (2016), DOI 10.1016/j.physletb.2016.08.038](https://pmc.ncbi.nlm.nih.gov/articles/PMC5207040/)
- Lochak cites the Re result and chemical effects "generally less than 1 %" and links them to Kadomtsev's prediction that superstrong fields change beta lifetimes; Urutskoev's uranium-solution work claims explosions accelerate β⁻ processes — [Lochak 2007](https://arxiv.org/abs/0801.2752); [Urutskoev et al. 2022](https://fondationlouisdebroglie.org/AFLB-471/aflb471m969.pdf)
- Half-lives used: ⁶⁰Co 1925.28 d, ¹³⁷Cs 30.08 y — [IAEA LiveChart](https://nds.iaea.org/relnsd/v1/data?fields=ground_states&nuclides=60co)

### Inferences
- If decays stop for Δt and the nuclei are kept, the activity afterwards is higher by exp(λΔt) − 1 (computed): 7.20×10⁻⁴ for ⁶⁰Co and 1.26×10⁻⁴ for ¹³⁷Cs after 2 days. A gamma line that is merely absorbed or redirected, without any change in decays, leaves no such step. Measuring before, during and after separates the two.
- A pole's field reaches the critical scale m_e²c²/(eħ) = 4.41×10⁹ T only within 273 fm of the nucleus (computed). A Kadomtsev-type field effect would therefore need the pole bound well inside the K shell, not merely nearby.

### Gaps
- No primary measurement of monopole- or discharge-induced decay-rate change was found with stated counting statistics.

## Q7. Track physics: what a moving Dirac monopole would leave, and what else makes such marks

### Takeaway
A fast Dirac monopole ionises like a nucleus of Z_eff = 68.5β, about 4700 times a unit charge (Ahlen). A slow one (β ~ 10⁻³) loses ~60–170 MeV/cm by the Ahlen–Kinoshita formula (computed). The claimed 0.7–1 GeV per 1–3 mm track is 14–163× more than that formula gives at the claimed β = 10⁻³. The real tracks lie in the emulsion plane, are "hairless", and appear on control films stored in a magnet. Those properties are not those of a straight ionising particle. Status: monopole interpretation **unverified**; Ahlen formulas **mainstream**.

### Cited Findings
- Fast monopoles (β > 0.05): dE/dx = (4πN_e g²e²/m_ec²)[ln(2m_ec²β²γ²/I) + K/2 − 1/2 − δ/2 − B_m]; "a monopole with magnetic charge g_D behaves as an equivalent electric charge (Ze)² = g_D²β² ∼ (68.5β)²" (Eq. 11). Slow monopoles (10⁻³ ≲ β ≲ 10⁻²), Ahlen–Kinoshita: dE/dx = (2πN_e g²e²β/m_ecv_F)[ln(2m_ev_Fa₀/ħ) − 1/2] (Eq. 12), times ~1.37 for the electron-spin coupling; v_F ≃ 3.9×10⁸ cm/s for Si at 4.3 g/cm³; at this low energy transfer "transitions yielding photons in the visible range are largely suppressed" — [Patrizii & Spurio, arXiv:1510.07125 (2015)](https://arxiv.org/abs/1510.07125)
- Original formulas: [Ahlen, PRD 17, 229 (1978), DOI 10.1103/PhysRevD.17.229](https://api.crossref.org/works/10.1103/PhysRevD.17.229); [Ahlen & Kinoshita, PRD 26, 2347 (1982), DOI 10.1103/PhysRevD.26.2347](https://api.crossref.org/works/10.1103/PhysRevD.26.2347) (both Crossref-checked; not opened)
- Ionisation ~4,700× an electron's; NTDs (CR-39 etc.) record the restricted energy loss as etch-pit cones — [Mavromatos & Mitsou, arXiv:2005.05100](https://arxiv.org/abs/2005.05100)
- Monopole energy loss is "several thousand times" a unit charge's; early NTD and induction candidates were never confirmed — [PDG 2024](https://pdg.lbl.gov/2024/reviews/rpp2024-rev-mag-monopole-searches.pdf)
- Urutskoev cites Matsumoto's "mesh-like traces on nuclear emulsions during cold fusion" (Fusion Tech. 23, 103, 1993) as "exactly the same tracks" — [Urutskoev 2002](https://fondationlouisdebroglie.org/AFLB-274/aflb274p701.pdf)

### Inferences
- Computed (Eq. 12 × 1.37, β = 10⁻³): Si at 2.33 g/cm³, 61 MeV/cm, giving 6–18 MeV over 1–3 mm; AgBr at 6.47 g/cm³, 169 MeV/cm, giving 17–51 MeV (free-electron-gas model; order of magnitude only for AgBr). The claimed 700–1000 MeV is 14–163× higher, and the paper's "Coulomb braking" calibration of darkening is not described.
- v_F for Si at 4.3 g/cm³ = 3.900×10⁸ cm/s (computed, matches Patrizii & Spurio).
- A straight-moving particle cannot follow the emulsion plane for 1–10 mm at a depth constant to 10–15 µm. Both Lochak and Ivoilov invoke trapping between the layers. Surface marks from static discharge, pressure or chemistry are confined to the surface by construction and are the standard alternative. Ivoilov's magnet-only controls show that the effect occurs without any discharge.
- A blind test: an etched CR-39 stack (calibrated on known ions) next to the discharge and a matched one far away. A Dirac monopole makes collinear etch cones through every sheet; surface artefacts do not.

### Gaps
- An authoritative radiographic-film artefact reference (static "tree/crown/smudge" marks) was not opened; search returned only teaching material.
- Matsumoto 1993 not opened.

## Q8. The "Pripyat aurora": air fluorescence from ionising radiation

### Takeaway
Ionising radiation makes air glow in the N₂ bands. The 337 nm band alone yields 5.61 photons per MeV deposited (AIRFLY). In a 1–100 Gy/h field over a 100 m path this gives 0.12–11.7 kR in the 337 nm band (computed), from below to well above the faintest visible aurora class (IBC I = 1 kR). A glow near the damaged reactor is therefore expected from known physics without monopoles. The specific claim, light seen by "thousands" in Pripyat and spectrally unusual, has no primary source. Status: ionised-air glow **mainstream**; the monopole interpretation **unverified**.

### Cited Findings
- AIRFLY: absolute yield of the 337 nm band 5.61 ± 0.06 (stat) ± 0.21 (syst) photons/MeV at 1013 hPa and 293 K — [Ave et al., Astropart. Phys. 42, 90 (2013), DOI 10.1016/j.astropartphys.2012.12.006](https://arxiv.org/abs/1210.6734)
- Auroral brightness classes: IBC I 1 kR (Milky Way), II 10 kR, III 100 kR, IV 1000 kR (full moon) — [Australian Space Weather Services, "The Aurora"](https://www.sws.bom.gov.au/Category/Educational/Space%20Weather/The%20Aurora/Aurora.pdf)
- Contamination dose rates at the Chernobyl plant on 20 May 1986 were 10–100 mR/h at Units 1–2 and 20–600 mR/h in the turbine hall (plant gallery, secondary) — [ChNPP gallery](https://gallery.chnpp.gov.ua/en/79-2011-06-10-13-28-26/-q-q)

### Inferences
- Column brightness (computed, 337 nm band only, isotropic, ρ_air = 1.204 kg/m³): 1 Gy/h gives 0.12 kR, 10 Gy/h 1.17 kR, 100 Gy/h 11.7 kR over a 100 m path. The full 300–430 nm emission is several bands, so these are lower bounds for blue-UV light. The calculator takes the dose rate as an input; no verified 26 April 1986 near-core dose rate was found.
- Lochak's paper offers a different Chernobyl monopole claim (²³⁵U "enriched up to 27 %"), not an aurora.

### Gaps
- No primary eyewitness or instrumental record of a Pripyat "aurora" or its spectrum was found; no verified near-core dose rate for the first night (INSAG-7 / UNSCEAR not opened).

## Q9. Calculators for the crates

### Takeaway
Fifteen calculators can be written with published test values. The monopole and Mössbauer ones belong in `magnetics`, decay, isotope and composition bookkeeping in `metrology`, air fluorescence in `plasmoid`, and the claim registry with status labels in `physica`. Each row below is exercised by `strange_radiation_checks.py`.

### Cited Findings

| Function (Rust signature) | Formula | Primary source | Published test value | Target crate::module | Status |
|---|---|---|---|---|---|
| `dirac_charge(n: i32) -> f64` (A·m) | g = n h/(μ₀e); g/(ec) = n/(2α) | PDG 2024 Eq. 95.2; Lochak 2007 Eq. 5.2 | 3.29×10⁻⁹ A·m (de Lacheze-Murel 2022); 68.5 (Mavromatos-Mitsou Eq. 11) | magnetics::monopole | mainstream (formula) |
| `monopole_field(r_m: f64, n: i32) -> f64` (T) | B = μ₀g/(4πr²) | PDG 2024 Eq. 95.1 | flux μ₀g_D = h/e (exact); no published point value | magnetics::monopole | mainstream |
| `monopole_energy_gain(b_t: f64, l_m: f64, n: i32) -> f64` (J) | W = gBL | de Lacheze-Murel 2022 fn. 1 | 205 keV for 1 G, 0.1 m | magnetics::monopole | mainstream |
| `pole_moment_energy(mu_nuclear: f64, r_m: f64) -> f64` (eV) | U = μ μ_N μ₀g/(4πr²) | built from PDG Eq. 95.1 + ENSDF moments | none (9.4 keV at 10 fm, illustrative) | magnetics::monopole | mainstream formula; bound-state claim unverified |
| `monopole_dedx_fast(beta, n_e, i_ev, n) -> f64` | Patrizii-Spurio Eq. 11 (Ahlen 1978) | Ahlen PRD 17, 229 | (g/e)² ≈ 4,700 (Mavromatos-Mitsou) | magnetics::monopole | mainstream |
| `monopole_dedx_slow(beta, rho, z, a, n) -> f64` | Patrizii-Spurio Eq. 12 ×1.37 (Ahlen-Kinoshita 1982) | Ahlen & Kinoshita PRD 26, 2347 | v_F(Si, 4.3 g/cm³) = 3.9×10⁸ cm/s | magnetics::monopole | mainstream |
| `fe57_sextet(b_t: f64) -> [f64; 6]` (mm/s) | v = c[−(μ₁/I₁)m₁ + (μ₀/I₀)m₀]μ_N B/E_γ | ENSDF moments (IAEA); Preston 1962 | μ₁/μ₀ = −1.715 ± 0.004 (Preston); H = 330 ± 3 kOe | magnetics::mossbauer | mainstream |
| `hyperfine_temperature_equivalent(dh_oe: f64, slope_oe_per_k: f64) -> f64` | ΔT = ΔH/(dH/dT) | Ivoilov & Urutskoev 2004 (slope from Preston 1962) | −40 Oe/K; claimed ΔH/H ≈ 8×10⁻⁴ (Urutskoev 2002) | magnetics::mossbauer | mainstream; claim unverified |
| `effective_field(b_hf, b_par, b_perp) -> f64` | \|B\| = √((B_hf + b∥)² + b⊥²) | vector addition; geometry argument (Q2) | none (sign-symmetry self-check) | magnetics::mossbauer | mainstream |
| `second_order_doppler_coeff(t_k, theta_d, a_amu) -> f64` (µm/s/K) | C_v(Debye)/(2Mc) | Preston 1962 (θ_D = 400 K) | 0.65 µm/s/K (Ivoilov & Urutskoev 2004) | magnetics::mossbauer | mainstream |
| `activity_step_after_suppression(t_half_d, dt_d) -> f64` | exp(ln2 Δt/T½) − 1 | decay law; IAEA half-lives | compare with Pommé 2016 limit 6×10⁻⁶–8×10⁻⁵ | metrology::decay | mainstream; claim unverified |
| `ec_half_life_shift(t1, t2) -> f64` | (T₁ − T₂)/T₁ | Ohtsuki 2004 | 0.83 % (⁷Be in C₆₀ vs Be) | metrology::decay | mainstream |
| `burnout_isotope_fraction(f, a48) -> f64` and `transmutation_energy_bound(n, de_mev) -> f64` | (a − f)/(1 − f); E = NΔE | Urutskoev 2002 Fig. 10, p. 723–724 | IAEA Ti abundances (sum 100 %); 50 kJ bank | metrology::composition | unverified (the claim); arithmetic mainstream |
| `steel_fingerprint(fe, cr, ni, steel) -> (f64, f64, f64)` | k = Fe/Fe_steel; Cr/(k Cr_steel), Ni/(k Ni_steel) | Priem et al. 2008 Table 3 | 18-8 steel 72/20/8 %; obs/expected 0.71–0.91 | metrology::composition | mainstream (method) |
| `air_fluorescence_column(dose_gy_h, path_m, yield_per_mev) -> f64` (R) | Ẋ = (Ḋ ρ/MeV) Y L / 10¹⁰ | AIRFLY 2013; Rayleigh/IBC scale | 5.61 photons/MeV (337 nm); IBC I = 1 kR | plasmoid::fluorescence | mainstream |

### Inferences
- Claim records for `physica` (status labels): Ti-foil transmutation (unverified); strange-radiation tracks as monopoles (unverified); ⁵⁷Fe hyperfine shift (unverified; size equals a 6–7 K systematic); Lochak monopole theory (contested); Klimov COP > 1 and transmutation (unverified); gamma suppression for 2 days (unverified, no source); solar modulation of decay (refuted); bound-state beta decay and EC chemistry (mainstream); ionised-air glow (mainstream).
- Rigs-as-examples rule: the Urutskoev bank (50 kJ, 4.8 kV), the Kazan arc (40 A, 80 V) and the Nantes Magneform shot (3 kJ, 4.2 kV) belong in examples that call the library functions, not in the library.

### Gaps
- No published numerical test value for a light-monopole cross-section, trapping energy or decay-rate effect exists to calibrate a "Lochak" calculator. Only the Dirac charge and the classical field relations are testable.
