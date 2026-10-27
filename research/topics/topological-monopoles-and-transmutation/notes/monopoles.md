# Magnetic monopoles: Dirac, GUT, emergent and synthetic, and the monopole claims in Greenyer's TeslaTech 2026 talk

Scope: real (Dirac, 't Hooft–Polyakov/GUT), emergent (spin ice, chiral-magnet hedgehogs) and synthetic (spinor BEC) monopoles; searches and bounds; Callan–Rubakov catalysis; strong-field QCD; monopole–nucleus binding; and each monopole claim in Bob W. Greenyer, "The Power and Signature of the Topological Monopole" (TeslaTech, 15 Aug 2026; transcript /home/dan/Desktop/transcript_1793004282390.txt, lines 159–440, 680–735, 1453–1495, 2510–2525). Extends research/topics/magnetics-and-mhd/report.md (`magnetics::frontier`: `dirac_charge`, `DIRAC_FLUX`, `flux_from_events`, `parker_bound`, `spin_ice_charge`), which is not redone here.

Every computed number below is recomputed and asserted in `monopoles_checks.py` (same folder; 160 checks, standard library, `python3 -I`, exits 1 on failure). Downloads are in `notes/dl_monopoles/` (git-ignored). Constants are CODATA 2022 from [NIST allascii.txt](https://physics.nist.gov/cuu/Constants/Table/allascii.txt). Nuclear data (moments, masses, mass excesses, Q values, radii, abundances) are from the [IAEA LiveChart API](https://www-nds.iaea.org/relnsd/v1/data?fields=ground_states&nuclides=57fe) (ENSDF/AME), fetched 2026-10-26. All DOIs below were resolved on api.crossref.org. Status labels: mainstream / contested / refuted / unverified, and "derived" for numbers I computed that no source prints.

## 1. Dirac 1931: quantization, string, field, flux, energy; 't Hooft–Polyakov and the GUT mass scale

### Takeaway
A Dirac monopole carries g_D = h/(μ₀e) = 3.291×10⁻⁹ A·m (flux μ₀g_D = h/e = 2Φ₀ = 4.136×10⁻¹⁵ Wb; 68.52 e in Gaussian units). Its field B = nh/(4πe r²) is 3.29×10¹⁴ T at 1 fm and 1.18×10⁵ T at the Bohr radius. A smooth, finite-core version exists only as the 't Hooft–Polyakov soliton of a broken non-abelian gauge theory, and in GUTs it has mass ~10¹⁷–10¹⁸ GeV. Mainstream.

### Cited Findings
- Dirac, "Quantised singularities in the electromagnetic field", Proc. R. Soc. Lond. A 133, 60–72 (1931); DOI checked on Crossref — [doi:10.1098/rspa.1931.0130](https://doi.org/10.1098/rspa.1931.0130). The full text was not opened. The equations below are quoted from the PDG and Castelnovo et al., not from Dirac.
- PDG review "94. Magnetic Monopoles" (Milstead & Weinberg, revised August 2025). Field eq. 94.1, B = Q_M r̂/(4πr²). A vector potential whose curl is B "must be singular along some line running from the origin to spatial infinity" (the Dirac string). Requiring the string to be unobservable for every pair of charges gives eq. 94.2, Q_E^min·Q_M^min = 2π (natural units), so Q_M^D ≡ 2π/e. The quantization also follows from the extra angular momentum −Q_EQ_M r̂/4π (eq. 94.3) being quantized in half-integer units. That one condition covers every electric charge, which is why one monopole quantizes all charge. Without colour-magnetic charge the minimum becomes 6π/e, because quark charges then enter — [PDG 2025](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-mag-monopole-searches.pdf).
- In SI form, Castelnovo, Moessner & Sondhi write the condition as e·q = nh/μ₀, "whence any monopoles in our universe must be quantised in units of q_D = h/μ₀e" — [arXiv:0710.5515](https://arxiv.org/abs/0710.5515) (Nature 451, 42 (2008), [doi:10.1038/nature06433](https://doi.org/10.1038/nature06433)).
- 't Hooft–Polyakov: a triplet Higgs breaks SU(2) to U(1). The Higgs field "twisting" at infinity "cannot be smoothly deformed to a spatially uniform vacuum configuration" (topological stability), which gives Q_M = 4π/e and M_mon ~ 4πM_V/e² (eq. 94.4). The standard electroweak theory has no topologically stable monopole. SU(5) gives Q_M = 2π/e with M_mon ~ 4πM_X/g² (eq. 94.6), "for a unification scale of 10¹⁶ GeV … M_mon ∼ 10¹⁷–10¹⁸ GeV" — [PDG 2025](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-mag-monopole-searches.pdf).
- Precision limits: the electron's magnetic charge is < 10⁻²⁴ Q_M^D (induction method), and muon g−2 gives a model-dependent mass floor of 120 GeV — [PDG 2025](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-mag-monopole-searches.pdf).
- Computed (checks §1): g_D = 3.29106×10⁻⁹ A·m; 1/(2α) = 68.518; μ₀g_D = h/e = 2Φ₀ = 4.13567×10⁻¹⁵ Wb; B(r) = μ₀g/(4πr²) ≡ nħ/(2er²): B(1 fm) = 3.291×10¹⁴ T, B(a₀) = 1.1753×10⁵ T. B reaches 10¹⁵ T at r = 0.574 fm *(derived)*.

### Inferences
- "Topological" means two different things here. For Dirac's point monopole it is the integer first Chern number of the U(1) bundle (the flux quantum h/e), with a point singularity. For the 't Hooft–Polyakov/GUT monopole it is the winding of the Higgs field (π₂), which gives a smooth core of size ~1/M_X. The talk's line "the Dirac monopole … is actually a topological monopole that cannot be smoothed away and at large distances appears like a point" fits the GUT soliton. A pure Dirac monopole is singular rather than smooth. Status: mainstream, with this terminology caveat.
- Every monopole bound inside or next to a nucleus puts that nucleus in a field above 10¹⁴ T (3.3×10¹⁴ T at 1 fm). This is the same order as the talk's 10¹⁵ T (see §5 and §7).

### Gaps
- I did not open Dirac's paper text, so no equation number from it is given.

## 2. Monopole-catalysed baryon decay (Rubakov, Callan), its searches, and the "John G. Cramer 1983" item

### Takeaway
Callan–Rubakov catalysis (p + M → e⁺ + π⁰ + M with σβ ~ 10⁻²⁷ cm², a strong-interaction size) is mainstream theory for GUT monopoles. The PDG calls it model-dependent and "not even a universal property of all GUT monopoles". No search has seen it. The Cramer item is his Analog column AV-01, "When Proton Meets Monopole" (written 16 Dec 1983, published July 1984). It says the monopole is "unaffected", not that it gets stronger.

### Cited Findings
- Callan, "Dyon-fermion dynamics", Phys. Rev. D 26, 2058–2068 (1982) — [doi:10.1103/PhysRevD.26.2058](https://doi.org/10.1103/PhysRevD.26.2058). Rubakov, "Adler-Bell-Jackiw anomaly and fermion-number breaking in the presence of a magnetic monopole", Nucl. Phys. B 203, 311–348 (1982) — [doi:10.1016/0550-3213(82)90034-7](https://doi.org/10.1016/0550-3213(82)90034-7). (Neither text was opened. Rubakov's 1981 JETP Lett. note has no Crossref DOI that I found.)
- PDG: a GUT monopole core contains the superheavy gauge fields that mediate baryon-number violation. The surprise found by Callan and Rubakov is that the rate "is not suppressed by powers of the gauge boson mass". Cross-sections for "p + monopole → e⁺ + π⁰ + monopole are essentially geometric; i.e., σ_ΔB β ∼ 10⁻²⁷ cm²". The PDG adds that catalysis "is model-dependent and is not even a universal property of all GUT monopoles". Astrophysical catalysis bounds are F·(σβ/10⁻²⁷ cm²) ∼ 10⁻¹⁸–10⁻²⁹ cm⁻² sr⁻¹ s⁻¹ (eq. 94.12) — [PDG 2025](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-mag-monopole-searches.pdf).
- IceCube: σ_cat = σ₀/β for β ≥ β₀, with a further F(β) factor below β₀ (eq. 4). Current estimates run "10⁻²⁷ cm² to 10⁻²¹ cm²", with λ_cat = 1/(σ_cat·n) (eq. 5). σ_cat = 1.7×10⁻²² and 1.7×10⁻²³ cm² "correspond to λ_cat = 1 mm and λ_cat = 1 cm". Limits: Φ₉₀ ≤ 10⁻¹⁸ (10⁻¹⁷) cm⁻² s⁻¹ sr⁻¹ for σ_cat = 10⁻²² (10⁻²⁴) cm², "three orders of magnitude below the Parker bound" — [arXiv:1402.3460](https://arxiv.org/abs/1402.3460) (EPJC 74, 2938 (2014)).
- Super-Kamiokande (neutrinos from catalysed decay in the Sun): F_M(σ₀/1 mb) < 6.3×10⁻²⁴ (β_M/10⁻³)² cm⁻² s⁻¹ sr⁻¹ at 90% C.L., "more than eight orders of magnitude more stringent" than 1×10⁻¹⁵ — [arXiv:1203.0940](https://arxiv.org/abs/1203.0940) (Astropart. Phys. 36, 131 (2012)).
- Other catalysis searches (Soudan, MACRO, IMB, Baikal, IceCube) give β-dependent limits of ~10⁻¹⁸–10⁻¹⁴ cm⁻² sr⁻¹ s⁻¹ — [PDG 2025](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-mag-monopole-searches.pdf).
- Cramer, "When Proton Meets Monopole", Alternate View AV-01: "Published in the July-1984 issue of Analog … written and submitted 12/16/83". Cramer's monopoles are Planck-mass (~10⁻⁸ g) black-hole remnants "stuck" with magnetic flux, not GUT solitons. Two short quotes: "But the monopole, the cause of it all, is unaffected." and "the monopole is the analog of a chemical catalyst". He also sketches monopole-catalysed power and a magnetic "hemi-bottle" rocket that collects the monopoles "for re-use in the next engine cycle". He says proton life expectancy is cut "by a factor of 10¹²". His only catalysis reference is "F. Wilczek, Phys. Review Letters 48, 1144 (1982)" — [npl.washington.edu/AV/altvw01.html](https://npl.washington.edu/AV/altvw01.html).
- Crossref gives PRL 48, 1144 as Wilczek, "Magnetic Flux, Angular Momentum, and Statistics" (the anyon paper), not a catalysis paper — [doi:10.1103/PhysRevLett.48.1144](https://doi.org/10.1103/PhysRevLett.48.1144).
- Computed (checks §5): the IceCube pair (1.7×10⁻²³ cm², 1 cm) implies n = 5.88×10²² cm⁻³. This matches the density of hydrogen nuclei in ice, 6.13×10²² cm⁻³ at an assumed 0.917 g/cm³, to 4%. 10⁻¹⁵/6.3×10⁻²⁴ = 1.59×10⁸. Cramer's 10³²/10¹² = 10²⁰ yr = 10¹⁰ × the 10¹⁰-yr age he quotes, which matches his "10 billion times the present age". His printed "about 10¹⁰ years" therefore looks like a typo for 10²⁰.

### Inferences
- Talk vs source. The talk says Cramer showed a monopole "will catalyze the decay of baryons and it won't affect the monopole. In fact, it will actually make it stronger." The first half matches Cramer, and Cramer did propose catalysis as an energy and propulsion source. "Stronger" is not in the column, and in Callan–Rubakov the monopole's charge is conserved. Status: catalysis mainstream (theory, unobserved); "makes it stronger" unsupported. The citation should read Cramer, Analog, July 1984 (written 1983).
- From IceCube's numbers, the published σ↔λ conversion uses free protons (hydrogen nuclei) as targets *(inferred; IceCube does not state n)*.
- Catalysis needs a GUT core, not a field threshold. Nothing in the Callan–Rubakov mechanism depends on reaching a given |B|.

### Gaps
- Rubakov 1981 (JETP Lett. 33, 644) and Callan's 1983 NPB paper were not opened. Neither the IceCube nor the Super-K paper prints its target density n.

## 3. Searches and bounds (Cabrera, Parker, MACRO, IceCube, ANTARES, MoEDAL, PDG)

### Takeaway
No confirmed monopole exists. Flux limits sit far below the Parker bound: MACRO 1.4×10⁻¹⁶, ANTARES ~7×10⁻¹⁸, IceCube 2.0×10⁻¹⁹ cm⁻² s⁻¹ sr⁻¹. In matter the limit is ≲10⁻²⁹ monopoles per nucleon. At colliders MoEDAL excludes Schwinger-produced monopoles up to 75 GeV (1–3 g_D) and 80 GeV (2–45 g_D), and Drell–Yan/photon-fusion searches reach ~4 TeV.

### Cited Findings
- Parker bound F < 10⁻¹⁵ cm⁻² sr⁻¹ s⁻¹ for M ≲ 10¹⁷ GeV (eq. 94.9); seed-field bound (eq. 94.10); galactic acceleration v_mag ~ c for M ≲ 10¹¹ GeV and 10⁻³c(10¹⁷ GeV/M)^½ above that (eq. 94.8) — [PDG 2025](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-mag-monopole-searches.pdf); Parker, ApJ 160, 383 (1970), [doi:10.1086/150442](https://doi.org/10.1086/150442).
- Bound in matter: "upper limit on the monopoles per nucleon ratio of ∼10⁻²⁹" (Jeon & Longo). Early induction (Cabrera) and track-detector candidates were never confirmed. NOvA: 2×10⁻¹⁴ for 6×10⁻⁴ < β < 10⁻³. Auger/RICE/ANITA reach 2.5×10⁻²¹ for ultra-relativistic monopoles. LHC mass limits reach "around 4 TeV" — [PDG 2025](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-mag-monopole-searches.pdf).
- Cabrera, PRL 48, 1378 (1982), [doi:10.1103/PhysRevLett.48.1378](https://doi.org/10.1103/PhysRevLett.48.1378). Its implied flux of 6.10×10⁻¹⁰ cm⁻² s⁻¹ sr⁻¹ is already in the repo (`flux_from_events`).
- MACRO final: "1.4×10⁻¹⁶ cm⁻² s⁻¹ sr⁻¹" for 4×10⁻⁵ < β < 1 — [arXiv:hep-ex/0207020](https://arxiv.org/abs/hep-ex/0207020) (EPJC 25, 511 (2002)).
- IceCube (2886 days): "below 2.0×10⁻¹⁹ cm⁻² s⁻¹ sr⁻¹" for 0.750c–0.995c — [arXiv:2109.13719](https://arxiv.org/abs/2109.13719) (PRL 128, 051101 (2022)).
- ANTARES (10 yr): "∼ 7×10⁻¹⁸ cm⁻² s⁻¹ sr⁻¹" for ultra-relativistic monopoles, β ≥ 0.55 — [arXiv:2202.13786](https://arxiv.org/abs/2202.13786) (JHEAp 34, 1 (2022)).
- MoEDAL Schwinger search, Pb–Pb at 5.02 TeV/nucleon, 0.235 nb⁻¹: "masses up to 75 GeV/c² were excluded" for 1g_D ≤ g ≤ 3g_D, "the first lower mass limit for finite-size MMs from a collider search". Peak field "B ≈ 10¹⁶ T, with an inverse decay time ω ≈ 10²⁶ s⁻¹". The fit gives B = 7.6 ± 0.3 GeV², ω = 73 ± 3 GeV (natural units c = ħ = ε₀ = 1) — [arXiv:2106.11933](https://arxiv.org/abs/2106.11933) (Nature 602, 63 (2022), [doi:10.1038/s41586-021-04298-1](https://doi.org/10.1038/s41586-021-04298-1)).
- MoEDAL scan of the CMS beam pipe, Pb–Pb at 2.76 TeV, 184.07 µb⁻¹: excludes masses up to 80 GeV for 2–45 Dirac units, with B_max = 4.0 ± 0.2 GeV² and ω = 40 ± 1.7 GeV. Point-like (GUT/BSM) production in such collisions "is exponentially suppressed by a factor of e^−O(500)" — [arXiv:2402.15682](https://arxiv.org/abs/2402.15682) (PRL 133, 071803 (2024)).
- Computed (checks §3, §5): 7.6 GeV² (Heaviside–Lorentz, 1 T = 195.35 eV²) = 3.89×10¹⁶ T. 4.0 GeV² = 2.05×10¹⁶ T. 73 GeV/ħ = 1.11×10²⁶ s⁻¹, a field lifetime of 9.0×10⁻²⁷ s. MACRO/Parker = 0.14; IceCube/Parker = 2×10⁻⁴. In a 3 µG field, a Dirac monopole gains the PDG's 10¹¹ GeV (16.02 J) over 0.526 kpc *(derived; the PDG does not state a length)*.

### Inferences
- The PDG eq. 94.8 boundary (M ≈ 10¹¹ GeV) matches g_D·B·L ≈ Mc² for L ≈ 0.5 kpc, a reasonable galactic coherence length. This gives a test vector for a `monopole_energy_gain` function.

### Gaps
- I did not extract the PDG flux-limit curves (Fig. 94.1) numerically.

## 4. Emergent and synthetic monopoles: spin ice, spinor BECs, knots, chiral-magnet hedgehogs

### Takeaway
All four are real, peer-reviewed observations. None is a free magnetic charge of the vacuum. The spin-ice charge is a real H-field source of 2μ/a_d ≈ 4.6 µ_B/Å = 4.28×10⁻¹³ A·m (measured 5 µ_B/Å), ~1/7700 of g_D. The BEC "monopoles" live in a synthetic field (superfluid vorticity), and the chiral-magnet ones in an emergent field acting on conduction electrons. In each case "topological" refers to an order-parameter winding number.

### Cited Findings
- Spin ice: an inverted dumbbell makes "two adjacent sites with net magnetic charge Q_α = ±q_m = ±2µ/a_d" (eq. 3), with µ ≈ 10 µ_B. Eq. 4 gives q_m = (2µ/µ_b)(αλ_C/2πa_d)q_D ≈ q_D/8000. The authors stress that spin-ice strings "are observable and thus not quantised" (no Dirac condition) — [arXiv:0710.5515](https://arxiv.org/abs/0710.5515).
- Bramwell et al. (muon spin rotation, Dy₂Ti₂O₇): "determine the elementary unit of magnetic charge to be 5 µB Å⁻¹", against the predicted 2µ/a = 4.6 µ_B/Å with a_d = 4.3356 Å — [arXiv:0907.0956](https://arxiv.org/abs/0907.0956) (Nature 461, 956 (2009), [doi:10.1038/nature08500](https://doi.org/10.1038/nature08500)).
- Ray et al. 2014 (spin-1 ⁸⁷Rb BEC): A* and B* = ħ∇×A* "are related to the superfluid velocity, v_s, and vorticity, Ω = ∇ × v_s". Eq. 2 gives Ω = −(ħ/M)r̂′/r′² + (4πħ/M)δ(x′)δ(y′)Θ(z′)r̂′, "a monopole attached to a semi-infinite vortex line singularity, of phase winding 4π" (the Dirac string, seen directly in images) — [arXiv:1408.3133](https://arxiv.org/abs/1408.3133) (Nature 505, 657 (2014), [doi:10.1038/nature12954](https://doi.org/10.1038/nature12954)).
- Ray et al. 2015: "Observation of isolated monopoles in a quantum field", Science 348, 544–547 — a monopole-like spin texture in ⁸⁷Rb ([Crossref summary](https://doi.org/10.1126/science.1258289)). The 2017 PRX follow-up explains that in the polar phase π₂(O_P) ≅ ℤ allows singular point defects ('t Hooft–Polyakov analogues, no string). In the ferromagnetic phase π₂ = 0, so the Dirac monopoles there "are not point defects … but in the associated synthetic magnetic field". The paper observes an isolated monopole decaying into a Dirac-monopole configuration — [arXiv:1611.07766](https://arxiv.org/abs/1611.07766) (PRX 7, 021023 (2017)); theory in Tiurev et al., [doi:10.1103/PhysRevA.93.033638](https://doi.org/10.1103/PhysRevA.93.033638).
- Hall et al. 2016, "Tying quantum knots": knot solitons in the polar phase of a ⁸⁷Rb spinor BEC, a "topologically non-trivial element of the third homotopy group" showing the Hopf fibration; ⁸⁷Rb mass m = 1.443×10⁻²⁵ kg — [arXiv:1512.08981](https://arxiv.org/abs/1512.08981) (Nature Phys. 12, 478 (2016), [doi:10.1038/nphys3624](https://doi.org/10.1038/nphys3624)).
- Chiral magnets: in MnSi the skyrmion lattice gives an emergent field B_eff = Φ₀Φ with "Φ₀ = h/e … the flux quantum for a single electron". The integrated skyrmion density per cell is −1, giving "B_eff ≈ 2.5 T", which was measured through the topological Hall effect (Δρ_xy ≈ 4.5 nΩcm) — [arXiv:0902.1933](https://arxiv.org/abs/0902.1933) (PRL 102, 186602 (2009)). Milde et al. imaged skyrmion lines merging. In the bulk this needs "a hedgehog-like spin structure with the properties of the elusive magnetic monopole" (an emergent monopole, or Bloch point) — [Science 340, 1076 (2013)](https://doi.org/10.1126/science.1234657).
- Computed (checks §9): q_m = 2(10µ_B)/4.3356 Å = 4.613 µ_B/Å = 4.278×10⁻¹³ A·m, equal to the repo's 2µ/(√(3/2)·3.54 Å). The measured 5 µ_B/Å = 4.637×10⁻¹³ A·m. Eq. 4 is exactly q_m/q_D = 20r_e/a_d, giving q_D/q_m = 7693 ("≈8000"). ⁸⁷Rb mass from AME = 1.443×10⁻²⁵ kg. The Dirac-string circulation is 4πħ/M = 9.183×10⁻⁹ m²/s *(derived)*. B_eff = 2.5 T with h/e per skyrmion means one skyrmion per 1.654×10⁻¹⁵ m² *(derived)*.

### Inferences
- The talk says these were "predicted and simulated in Bose-Einstein condensate". Status: mainstream, and understated, because they were created and imaged (2014, 2015, 2017). "At the extreme such structures with requisite flux intensity can cause the collapse of baryons": the BEC monopole's "flux" is vorticity (units m²/s), and the emergent-magnet flux acts only on spin-polarised carriers. Neither carries real magnetic charge, so neither can drive Callan–Rubakov catalysis, which needs a GUT core. Status: refuted for analogue monopoles; for real GUT monopoles see §2.
- A spin-ice monopole is a real source of H inside the crystal, ~7700× weaker than g_D, and it cannot leave the lattice.

### Gaps
- No primary source gives a BEC synthetic-monopole "charge" in SI. 4πħ/M is my reading of Ray et al. eq. 2. The Neubauer extract does not state the skyrmion lattice period, so I give the implied area rather than a forward test.

## 5. "With sufficient magnetic strength, approximately 10¹⁵ T, baryonic matter is disassembled into quark-gluon plasma"

### Takeaway
10¹⁵ T is indeed of order the QCD field scale: eB = 0.059 GeV² ≈ 3 m_π², and B²/2µ₀ equals Λ⁴ for Λ ≈ 372 MeV. But lattice QCD shows a static field of that size does not deconfine cold matter. At eB = 9 GeV² (1.5×10¹⁷ T, 150× larger) the theory is still confined below T_c = 63 MeV. Heavy-ion collisions briefly reach 10¹⁴–4×10¹⁶ T. Status: refuted as a zero-temperature deconfinement threshold. The scale coincidence is real.

### Cited Findings
- Skokov, Illarionov & Toneev (UrQMD): "eB_y ∼ 10⁻¹·m_π² for the SPS energy range and eB_y ∼ m_π² for the RHIC energies … LHC … eB_y ∼ 15 m_π²". Their conversion is "m_π² = 140² × 0.512·10¹⁴ Gauss ≈ 10¹⁸ Gauss" — [arXiv:0907.1396](https://arxiv.org/abs/0907.1396) (IJMPA 24, 5925 (2009), [doi:10.1142/S0217751X09047570](https://doi.org/10.1142/S0217751X09047570)).
- Deng & Huang (HIJING): Au+Au at RHIC "can reach eB ∼ m_π² ∼ 10¹⁸ G", with sizeable event-by-event fields in all directions — [arXiv:1201.5108](https://arxiv.org/abs/1201.5108) (PRC 85, 044907 (2012)).
- Bali et al. (N_f = 2+1, physical masses, continuum-extrapolated): "the transition temperature T_c significantly decreases with increasing magnetic field" (inverse catalysis), and "the analytic crossover … persists up to our largest magnetic fields eB ≈ 1 GeV²" — [arXiv:1111.4956](https://arxiv.org/abs/1111.4956) (JHEP 02 (2012) 044, [doi:10.1007/JHEP02(2012)044](https://doi.org/10.1007/JHEP02(2012)044)).
- D'Elia et al.: the crossover "turns into a first order transition in the presence of large enough magnetic background fields". The critical endpoint lies between eB = 4 GeV² (T_c = 98 ± 3 MeV) and eB = 9 GeV² (T_c = 63 ± 5 MeV) — [arXiv:2111.11237](https://arxiv.org/abs/2111.11237) (PRD 105, 034511 (2022)).
- MoEDAL ultraperipheral Pb–Pb fields of ≈10¹⁶ T lasting ~1/ω (§3) — [arXiv:2106.11933](https://arxiv.org/abs/2106.11933).
- Magnetars: the repo already has SGR 1806−20 at 2.06×10¹⁵ G (2.1×10¹¹ T) (`pulsar_surface_field`, research/topics/magnetics-and-mhd/report.md).
- Computed (checks §2–3): eB = 1 GeV² ⇔ B = 1.690×10¹⁶ T. eB = m_π² (140 MeV) ⇔ 3.313×10¹⁴ T. Skokov's identity treats B (Heaviside–Lorentz) = m_π², not eB = m_π². The two differ by 1/e = 3.302, so "eB ∼ m_π² ∼ 10¹⁸ G" in both heavy-ion papers is loose by that factor *(derived)*. 10¹⁵ T ⇔ eB = 0.0592 GeV² = 3.02 m_π². Bali's 1 GeV² = 1.69×10¹⁶ T; D'Elia's 9 GeV² = 1.52×10¹⁷ T. B²/2µ₀ at 10¹⁵ T = 3.98×10³⁵ J/m³ = 2483 MeV/fm³. That is 22.4× the mass-energy density of a ⁵⁶Fe nucleus (110.7 MeV/fm³, uniform sphere of radius √(5/3)·3.7377 fm). The field whose energy density equals ⁵⁶Fe's mass-energy density is 2.11×10¹⁴ T, and the one equal to its binding-energy density is 2.05×10¹³ T.

### Inferences
- Computed from the standard side, the claim fails on lattice QCD. At T ≈ 0 the vacuum stays in the confined, chirally broken phase for every field yet simulated, up to 1.5×10¹⁷ T. The field lowers T_c (inverse catalysis) but does not remove it. A field above 10¹⁵ T would at most lower the temperature needed to make QGP, still well above 60 MeV per particle.
- Computed from the claim's side: the energy-density argument gives a real ordering. At 10¹⁵ T the field energy density is ~20× nuclear mass-energy density, so any static structure holding such a field would itself weigh far more than the matter it "disassembles". The heavy-ion fields that exceed 10¹⁵ T last only ~10⁻²⁶ s.
- "Electron-quark-gluon plasma" condensing back into "standard stable chemical elements": no primary source found. Status: unverified.

### Gaps
- I opened no source giving Λ_QCD in MeV (the PDG 2025 QCD review defines Λ only as an integration constant), so §10 gives the Λ⁴ condition as a function of Λ. Lattice results at T → 0 for eB > 9 GeV² were not found.

## 6. The bond-strength ladder made quantitative

### Takeaway
The ordering in the talk (van der Waals < hydrogen < metallic ≈ covalent ≪ β < α < cluster < nucleon rest energy) is correct, spanning 0.012 eV to 938 MeV. The field at which magnetism matches each rung depends on the criterion. μ_B·B = E gives 212 T (vdW) to 7.8×10⁴ T (covalent). B²/2µ₀ = E/V (V = (3.761 Å)³) gives 9.6–185 T for chemistry. For nuclear rungs μ_N·B = E gives 5.9×10¹¹–3.0×10¹⁶ T. The standard result in strong fields is the opposite of disassembly: atoms and molecules become more tightly bound.

### Cited Findings
- vdW: Ar₂ well depth −99.2(10) cm⁻¹ at R = 3.761(3) Å (spectroscopic, Herman et al.), against ab initio −99.19 to −99.60 cm⁻¹ — [arXiv:2304.14719](https://arxiv.org/abs/2304.14719) (Lang, Przybytek & Lesiuk 2023).
- Hydrogen bond: water dimer "D0 = 1105 ± 10 cm−1 (13.2 ± 0.12 kJ/mol)" — Rocher-Casterline et al., J. Chem. Phys. 134, 211101 (2011), [doi:10.1063/1.3598339](https://doi.org/10.1063/1.3598339) ([PDF](https://dornsife.usc.edu/reisler-lab/wp-content/uploads/sites/443/2024/07/101.pdf)).
- Metallic: Δ_fH°gas(Al) = 330.0 ± 4.0 kJ/mol (CODATA, Cox, Wagman et al.) — [NIST WebBook Al](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7429905&Mask=1); Δ_fH°gas(Fe) = 415.47 kJ/mol (Chase 1998) — [NIST WebBook Fe](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7439896&Mask=1). For an element, Δ_fH° of the gas is its 298 K cohesive (sublimation) enthalpy.
- Covalent: Δ_fH°gas(H) = 217.998 ± 0.006 kJ/mol — [NIST WebBook H](https://webbook.nist.gov/cgi/cbook.cgi?ID=C12385136&Mask=1); Δ_fH°gas(C) = 716.68 ± 0.45 kJ/mol — [NIST WebBook C](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7440440&Mask=1).
- Nuclear (IAEA LiveChart/AME): Q_β(³H) = 18.59202 keV, Q_β(n) = 782.347 keV, Q_α(²¹⁰Po) = 5407.5307 keV, B/A(⁵⁶Fe) = 8790.3563 keV; mass excesses ²²³Ra 17233.245, ²⁰⁹Pb −17614.574, ¹⁴C 3019.89328 keV — [LiveChart](https://www-nds.iaea.org/relnsd/v1/data?fields=ground_states&nuclides=223ra). Proton m_pc² = 938.27208943 MeV (CODATA).
- Strong fields bind matter more tightly. At B = 10¹² G (10⁸ T) the H atom is bound by 161 eV (13.6 eV at B = 0), and H₂ dissociation is Q₂ ≃ 46 eV against 4.75 eV at zero field (Lai, Table I and text). B₀ = 2.3505×10⁹ G — [arXiv:astro-ph/0009333](https://arxiv.org/abs/astro-ph/0009333) (Rev. Mod. Phys. 73, 629 (2001), [doi:10.1103/RevModPhys.73.629](https://doi.org/10.1103/RevModPhys.73.629)). Near B₀ a new "paramagnetic bonding" appears, drawing "spin-parallel hydrogen atoms or ground-state helium atoms … together into pairs" — Lange et al., Science 337, 327 (2012), [doi:10.1126/science.1219703](https://doi.org/10.1126/science.1219703) (editor's summary).
- Computed ladder (checks §4; chemistry V = (3.761 Å)³, an assumption):

| Rung | Energy E | μB = E field | B²/2µ₀ = E/V field |
|---|---|---|---|
| vdW, Ar₂ D_e | 0.012299 eV | 212.5 T (µ_B) | 9.65 T |
| H-bond, (H₂O)₂ D₀ | 0.13700 eV | 2367 T | 32.2 T |
| metallic, Al cohesive | 3.4202 eV | 5.909×10⁴ T | 160.9 T |
| metallic, Fe cohesive | 4.3060 eV | 7.439×10⁴ T | 180.5 T |
| covalent, C–C (716.68/2 kJ/mol) | 3.7139 eV | 6.416×10⁴ T | 167.7 T |
| covalent, H–H (2×217.998 kJ/mol) | 4.5188 eV | 7.807×10⁴ T | 184.9 T |
| β, ³H Q | 18.592 keV | 5.898×10¹¹ T (µ_N) | — |
| β, n Q | 782.35 keV | 2.482×10¹³ T | — |
| α, ²¹⁰Po Q | 5.4075 MeV | 1.715×10¹⁴ T | — |
| binding, ⁵⁶Fe B/A | 8.7904 MeV | 2.788×10¹⁴ T | 2.05×10¹³ T (nuclear density) |
| cluster, ²²³Ra → ²⁰⁹Pb + ¹⁴C Q | 31.828 MeV (from mass excesses) | 1.010×10¹⁵ T | — |
| nucleon rest energy | 938.272 MeV | 2.976×10¹⁶ T (µ_N); 1.066×10¹⁶ T (µ_p = 2.7928 µ_N) | 2.11×10¹⁴ T (⁵⁶Fe mass density) |

- The atomic unit of field, ħ/(ea₀²) = 2.3505×10⁵ T (CODATA), equals Lai's B₀ — [CODATA](https://physics.nist.gov/cuu/Constants/Table/allascii.txt).

### Inferences
- On the claim's side, the talk's 10¹⁵ T sits where µ_N·B equals cluster-decay Q values (1.0×10¹⁵ T for 31.8 MeV) and somewhat below the nucleon rest-energy rung. As an order-of-magnitude ladder of "magnetic energy per moment against binding energy" it is consistent.
- On the standard side, matching an energy scale does not make bonds fail. Lai and Lange show fields of 10⁵–10⁸ T increase atomic and molecular binding by an order of magnitude or more and create new bonds. "Bonds will fail based on their strength … broken by excitement and structure-based vacuum energy gradients" has no primary source. Status: energies mainstream; the field-driven disassembly sequence unverified, and contradicted for atoms and molecules by strong-field atomic physics.
- The anecdotes (gecko, gypsum in tornadoes, Shoulders, Usherenko) give no field measurements, so they cannot be placed on the ladder.

### Gaps
- The energy-density column needs a volume per bond. I used one common volume (Ar₂ R³) for every chemistry rung, which is an assumption. Bond-specific volumes would raise the covalent fields by up to ~(3.761/0.74)^{3/2} (not computed).

## 7. Monopole–nucleus "paramagnetic bound states", spectral lines, NMR/Mössbauer, and oxygen

### Takeaway
Bound states of monopoles with nuclei having a suitable magnetic moment are in the mainstream literature: p 15–320 keV, ²⁷Al 0.56–2.6 MeV. The binding criteria, however, exclude ⁵⁷Fe, ¹⁶O (no moment) and probably ¹⁷O. If an electron orbits a nucleus carrying a monopole, its levels shift by eV: the ground state goes from −13.6 to −7.29 eV, and spin Zeeman terms are ~6.8 eV. So "different spectral lines" would follow if monopoles were bound. Natural matter holds ≲10⁻²⁹ monopoles per nucleon, and the cited Fe-57 data are a single unreplicated Mössbauer study. Status: binding theory mainstream but "primitive" (Milton's word); the talk's applications unverified.

### Cited Findings
- Milton's review §9. Spin-½ binding in the lowest state J = l − ½ requires γ̂ > 1 + 1/4l (eq. 9.7), with γ̂ = (A/Z)γ and l = |m′|Z. Higher states require |κ̂| > (J² + J − l²)/l (eq. 9.8). "The proton, with κ = 1.79, will bind." For S = 3/2, eq. 9.12 says ⁹Be (γ̂ = −2.66) "cannot have binding". For S = 5/2, eq. 9.13 says ²⁷Al (γ̂ = 7.56) binds, with 1.67 > γ_c⁻ > 1.374 and γ_c⁺ < 4.216. Table 1 (binding energies, hard core or form factor): n 350 keV; ¹H 15.1 keV (Bracci & Fiorentini), 320 keV, 50–1000 keV, 263 keV (relativistic); ³He 13.4 keV; ²⁷Al 2.6 MeV (form factor) or 560 keV (hard core); ¹¹³Cd 6.3 keV. Milton's caveats: the theory "is rather primitive". Results are "very sensitive" to the regulator. Even so, "monopoles will bind to all nuclei, even, for example, Be, because the magnetic field in the vicinity of the monopole is so strong that the monopole will disrupt the nucleus". Kazama–Yang deeply bound states imply "electrons can bind to monopoles". The Coulomb–monopole spectrum (eq. 9.1) is E_Nj = −(µq²/2)[N + ½ + ((j+½)² − m′²)^½]⁻² — [arXiv:hep-ex/0602040](https://arxiv.org/abs/hep-ex/0602040) (Rep. Prog. Phys. 69, 1637 (2006), [doi:10.1088/0034-4885/69/6/R02](https://doi.org/10.1088/0034-4885/69/6/R02)).
- MoEDAL: aluminium, "a nucleus with large and positive anomalous magnetic moment", was chosen as the trapping material for this reason. A state bound at the keV level "will be stable for around 10,000 years" in the lattice — [arXiv:1405.7662](https://arxiv.org/abs/1405.7662) (IJMPA 29, 1430050 (2014)).
- Primary binding papers (verified on Crossref, not opened): Sivers, PRD 2, 2048 (1970), [doi:10.1103/PhysRevD.2.2048](https://doi.org/10.1103/PhysRevD.2.2048); Kazama, Yang & Goldhaber, PRD 15, 2287 (1977), [doi:10.1103/PhysRevD.15.2287](https://doi.org/10.1103/PhysRevD.15.2287); Bracci & Fiorentini, PLB 124, 493 (1983), [doi:10.1016/0370-2693(83)91559-9](https://doi.org/10.1016/0370-2693(83)91559-9); Goebel, "Binding of Nuclei to Monopoles", NATO ASI (1984), [doi:10.1007/978-1-4757-0375-7_30](https://doi.org/10.1007/978-1-4757-0375-7_30); Olaussen & Sollie, NPB 255, 465 (1985), [doi:10.1016/0550-3213(85)90147-6](https://doi.org/10.1016/0550-3213(85)90147-6).
- Nuclear moments (LiveChart/ENSDF): ⁵⁷Fe 1/2⁻, µ = +0.0902 µ_N, abundance 2.119%; ¹⁶O 0⁺, no moment, 99.757%; ¹⁷O 5/2⁺, µ = −1.893543 µ_N, 0.03835%; ²⁷Al 5/2⁺, +3.6407 µ_N; ⁹Be −1.17743; ³He −2.12762531; ¹¹³Cd −0.62226; ²H 0.857438231 — [LiveChart](https://www-nds.iaea.org/relnsd/v1/data?fields=ground_states&nuclides=17o).
- Urutskoev's "NMR" data are in fact Mössbauer data. Ivoilov & Urutskoev used conversion-electron Mössbauer spectroscopy (CEMS, depth ~0.2 µm) on 20 µm iron foils enriched to 94% ⁵⁷Fe, 15×15 mm, exposed to wire explosions in liquid (150 kA, pulsed 120–150 Oe on the samples). H_eff ≈ 330 kOe, and the 3-hour standard deviation is ±50 Oe. Sample 6(N) relaxed with τ = 16 ± 5 h from +280 ± 70 Oe. Sample 5(S) showed a single early point displaced by 460 ± 70 Oe. The control and sample 7(L) showed nothing. "Second measurement … 2 days later did not show any difference". dH_eff/dT = −40 Oe/K, with temperature stabilised to ±0.5 °C. They conclude the field "separates N- and S-monopoles" — [Ann. Fond. L. de Broglie 29 (2004) 1177](https://fondationlouisdebroglie.org/AFLB-297/aflb297m331.pdf). The 2008 overview interprets these as Lochak's light leptonic monopole — [arXiv:0801.2752](https://arxiv.org/abs/0801.2752) (not opened beyond the search abstract).
- Bound-in-matter limit ∼10⁻²⁹ monopoles per nucleon — [PDG 2025](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-mag-monopole-searches.pdf).
- Computed (checks §6–8). γ̂ reproduces Milton's Table 1 for ¹H, ³He, ⁹Be, ²⁷Al, ¹¹³Cd and ²H to 0.6%, using γ̂ = (A/Z)·µ/µ_N; Table 1's γ column is µ/µ_N. Eq. 9.12 at l = 2 gives 1.557/8.943 and eq. 9.13 at l = 6.5 gives 1.374/4.216, as printed. ⁵⁷Fe: γ̂ = 0.198 < 1.019 (lowest state fails), and |κ̂| = 0.80 < κ_c(J = l + ½) = 2.058, so no binding by these criteria. ¹⁷O: γ̂ = −4.024, negative like ⁹Be, so no binding by the same reading (γ_c⁻(l = 4) = 1.295). ¹⁶O has no dipole moment. Eq. 9.1 with m′ = ½ (electron with a Dirac monopole on a Z = 1 nucleus, spinless) gives −7.291 eV, 0.536 of the no-monopole value. The spin Zeeman scale µ_B·B_Dirac(a₀) is 6.80 eV. A triplet O₂ molecule (2µ_B) at 3 Å from a Dirac monopole has |µB| = 0.423 eV *(derived)*. A Dirac monopole's field equals the reported 280 Oe at 108.4 nm and 460 Oe at 84.6 nm. Filling the 0.2 µm × 15 × 15 mm CEMS layer at that spacing needs 8.43×10⁹ monopoles: 1.50 kg if they were GUT monopoles (10¹⁷ GeV), or 202.6 J of pair-creation energy at MoEDAL's 75 GeV floor. 280 Oe equals a 7 K temperature step. The shifts are 4.0σ and 6.6σ of the quoted ±70 Oe. A Dirac monopole in 120 Oe feels 3.95×10⁻¹¹ N.

### Inferences
- On the claim's side: a nucleus–monopole pair does change atomic spectra at the eV scale, the same scale as optical and soft-X-ray lines. "Different spectral lines" would therefore be a genuine signature, and the talk's physics here is qualitatively sound. The direction is uncertain. The spinless orbit is less bound (pushed out, not "pulled in"), while the spin–monopole term (~6.8 eV) can be attractive for one spin state (Kazama–Yang).
- On the standard side: for an XRF spectrum to show "none of the normal lines", a sizeable fraction of atoms would need a bound monopole. Even a 10⁻⁹ fraction is 10²⁰ times the natural bound-in-matter limit of 10⁻²⁹ per nucleon. No independent measurement of the Hutchison sample exists.
- Fe-57: the binding criteria in the review the talk relies on say ⁵⁷Fe does not bind a monopole through its dipole moment. A uniformly scattered population of radial fields averages to zero along the in-plane hyperfine axis, which broadens lines more than it shifts the mean H_eff. A 280 Oe step is also within 7 K of temperature drift. Status: unverified (one group, one venue, effect seen only in the first hours, no replication found).
- Chernobyl "monopoles bound to the oxygen" causing aurora: ¹⁶O (99.76%) has no nuclear moment, and ¹⁷O fails the criterion. Electronic (paramagnetic) attraction of O₂ to a monopole is about 0.4 eV at 3 Å, which is above kT, so a weak molecular "paramagnetic bound state" is possible if monopoles existed. No primary source for the aurora observation was found. Status: unverified.
- Mainstream "paramagnetic bonding" (Lange 2012) needs ~10⁵ T and involves no monopoles. It should not be confused with the talk's "paramagnetic bound states".

### Gaps
- None of the Bracci–Fiorentini, Goebel, Olaussen–Sollie or Sivers papers were opened, so Table 1 numbers are as reported by Milton. The S > ½ binding formula (eq. 9.9) did not extract cleanly, so I did not reproduce the treatment of negative γ̂ for S = 5/2 directly. The ¹⁷O conclusion follows the review's Be example and is labelled as such.

## 8. The "magnetic flux capacitor": field energy of a monopole outside a core r₀

### Takeaway
U(r₀) = µ₀g²/(8πr₀) = n²ħc/(8αr₀) = 3380 MeV·fm/r₀. For a GUT core r₀ = ħc/M_X this gives U = M_X/(8α) = 17.1 M_X, which is 1.7×10¹⁷ GeV for M_X = 10¹⁶ GeV and inside the PDG's 10¹⁷–10¹⁸ GeV. For a fixed charge the stored energy grows only as the core shrinks. Mainstream (derived).

### Cited Findings
- Field and charge as in §1 — [PDG 2025](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-mag-monopole-searches.pdf); M_mon ~ 4πM_X/g² (eq. 94.6), "M_mon ∼ 10¹⁷–10¹⁸ GeV" — same source.
- The talk: "Phase singularities produce a core with increasing magnetic energy. This is the flux capacitor. A topological magnetic flux, not electric, capacitance." (transcript ~l. 248–256).
- Computed (checks §1): µ₀g_D²/(8π) = ħc/(8α) = 5.4155×10⁻²⁵ J·m = 3380.1 MeV·fm. U = m_pc² at r₀ = 3.602 fm. U = 1 J at r₀ = 5.42×10⁻²⁵ m. 1/(8α) = 17.13.

### Inferences
- The classical self-energy reproduces the PDG GUT mass to within the order of magnitude. This independently checks why GUT monopoles are so heavy.
- A "capacitor" that keeps increasing its magnetic energy would need either more magnetic charge (quantized in g_D, and conserved) or a shrinking core. Neither process is described in the talk. Status: as a description of a real monopole, the energy formula is mainstream. "Increasing magnetic energy" from vortex focusing is unverified.

### Gaps
- None.

## 9. The talk's remaining monopole claims, both ways

### Takeaway
Monopoles made by "fractal electrodynamic structures", "experimentally established": no primary evidence was found, and classical electrodynamics forbids it for any current distribution. The "1.36 kT" paper is a particle-in-cell simulation reporting "up to 1000 Tesla". The "166× D–D fusion" figure does not come out of any standard comparison, which gives 157×, 257× or 1028×.

### Cited Findings
- Yue et al., "Dynamics of moving electron vortices and magnetic ring in laser plasma interaction", Phys. Plasmas 28, 042303 (2021). The work is "multi-dimensional Particle-in-Cell simulations". "In three dimensions, the structure appears as an expanding magnetic ring with an internal magnetic field up to 1000 Tesla", carrying >5% of the laser energy up to 50 µm — [doi:10.1063/5.0034098](https://doi.org/10.1063/5.0034098) (Crossref abstract).
- The talk calls this "a mainstream paper replicating 1976 work in Japan … 1.36 thousand Tesla on that ring … coming to a point" (transcript ~l. 382–396).
- Talk: "Over the past 10 years, it has been experimentally established that fractal electrodynamic structures can lead to the formation of magnetic charges in the form of topological monopoles" (~l. 181–185). No citation is given.
- Talk: "938 million electron volts per proton, or 16 times if it's an oxygen. Just one proton per nucleon gives you 166 times more energy than deuterium deuterium fusion" (~l. 2515–2520).
- Computed (checks §5) from AME mass excesses: D+D → T+p Q = 4032.66 keV; D+D → ³He+n 3268.91 keV; D+D → ⁴He 23846.5 keV. m_pc² ÷ (mean D–D Q per reaction) = 257.0; ÷ (mean per nucleon) = 1028; ÷ (D+D → ⁴He per nucleon) = 157.4. The ¹⁶O atom's mass-energy is 14899.2 MeV, not 16 × 938.272 = 15012.35 MeV (0.76% less).

### Inferences
- Magnetic charge from currents: B = ∇×A for any current density, so ∇·B = 0 and no arrangement of currents, toroidal or fractal, can make net magnetic charge. Toroidal (anapole) moments are higher multipoles with zero monopole term. The field of a current ring is finite at its centre for finite current. Status: refuted within Maxwell electrodynamics; unverified as an experimental claim (no primary source).
- Yue 2021: the field is ~10³ T in a simulation, 12 orders of magnitude below the 10¹⁵ T the talk needs. Status: the paper is mainstream; the talk's description of it as an experimental replication is not supported by its abstract.
- 938 MeV per catalysed nucleon is the correct order for Callan–Rubakov (§2), and the comparison with fusion is right in direction. The specific 166 is unexplained (status: unverified arithmetic). The closest standard value is 157 (per nucleon, D+D → ⁴He).

### Gaps
- I did not find which figure in Yue 2021 shows 1.36 kT, or which 1976 Japanese work is meant. I found no primary source for Pripyat aurora observations, the Klimov/Filimonenko decay-suppression claims, or Moscow (MEPhI) monopole observations. These were out of the time budget.

## 10. The user's audit: its formulas checked

### Takeaway
The audit (research/external/2026-10-26-greenyer-audit-and-crate.md) has five problems. (1) Its g_D = 2πħ/e is the flux h/e in Wb, and putting that into B = µ₀g/(4πr²) makes every field µ₀ ≈ 1.26×10⁻⁶ too small. (2) Its vector form of the Dirac potential has the opposite sign to its own spherical form. (3) The "Λ⁴" criterion gives 10¹⁵ T only for Λ ≈ 372 MeV. (4) Its DisassemblyRegime thresholds follow neither the Zeeman nor the energy-density criterion. (5) "Kramer 1983" should be Cramer (Analog, July 1984), and p + M → e⁺ + π⁰ + M is correct as written but has no field threshold.

### Cited Findings
- Audit text (data, unverified): "g_n = n 2πħ/e = n g_D". A(r) = (µ₀g/4π)(1 − cosθ)/(r sinθ) φ̂ = (µ₀g/4πr)(r × n̂)/(r + r·n̂), "Dirac string … along the negative z-axis", B = µ₀g r/(4πr³). The code's docstring says "g_D = 2πħ/e (A*m or Weber/mu_0)". B_crit ≈ 10¹⁵ T "where the magnetic energy density equals the QCD confinement scale Λ_QCD⁴". DisassemblyRegime: ≥10 T vdW, ≥100 T hydrogen bonds, ≥10³ T metallic, ≥10⁴ T covalent, ≥10⁷ T isotope stability, ≥10¹⁰ T nuclear cluster, ≥10¹⁵ T EQGP. "Monopole-Catalyzed Baryon Decay (Callan-Rubakov / Kramer 1983) … p + M → e⁺ + π⁰ + M" — /home/dan/code/FusionVault/quantum-gravimeter/research/external/2026-10-26-greenyer-audit-and-crate.md (lines 58–70, 140–145, 337–360, 536–575).
- SI Dirac condition e·q = nh/µ₀, q_D = h/(µ₀e) — [arXiv:0710.5515](https://arxiv.org/abs/0710.5515); PDG field eq. 94.1 — [PDG 2025](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-mag-monopole-searches.pdf).
- Computed (checks §2):
  1. 2πħ/e = 4.13567×10⁻¹⁵ Wb = h/e. The audit's B divided by the correct B is exactly µ₀ = 1.2566×10⁻⁶. The two consistent SI forms are: (a) g in A·m, e·g = nh/µ₀, g_D = h/(µ₀e) = 3.291×10⁻⁹ A·m, B = µ₀g/(4πr²); or (b) q_m in Wb (= µ₀g), e·q_m = nh, q_D = h/e, B = q_m/(4πr²) with no µ₀. The audit's condition eg/(2πħ) = n is form (b), but its field formula is form (a).
  2. A finite-difference curl at (0.3, −0.4, 0.5) gives, for the spherical form (µ₀g/4π)(1 − cosθ)/(r sinθ) φ̂, curl_z = +µ₀g z/(4πr³) (correct, string on −z). For the vector form (µ₀g/4πr)(r × n̂)/(r + r·n̂), curl_z = −µ₀g z/(4πr³). With n̂ = ẑ, r × ẑ = −ρφ̂, so the vector form equals minus the spherical form. Its denominator r + r·n̂ vanishes on the half-line along −n̂, so the string is on −n̂, as the audit says. The correct vector form is A = (µ₀g/4π)(n̂ × r)/(r(r + r·n̂)), giving B = +µ₀g r/(4πr³) with the string on −n̂. The audit's Rust `vector_potential` uses r × n̂ and has the same sign error.
  3. B²/(2µ₀) = Λ⁴/(ħc)³ gives B = 2.90×10¹⁴ T (Λ = 200 MeV), 6.52×10¹⁴ T (300 MeV) and 8.37×10¹⁴ T (340 MeV). B = 10¹⁵ T needs Λ = 371.7 MeV. The alternative criterion eB = Λ² gives 10¹⁵ T at Λ = 243.2 MeV. Either way the equality is an order-of-magnitude coincidence, not a threshold (§5: lattice QCD keeps cold matter confined to ≥1.5×10¹⁷ T).
  4. Against the §6 ladder: vdW 10 T matches the energy-density field for Ar₂ in (3.761 Å)³ (9.65 T), but the Zeeman field is 212 T. H-bond 100 T lies between 32.2 T (energy density) and 2367 T (Zeeman). Metallic 10³ T lies between 161–181 T and 5.9–7.4×10⁴ T. Covalent 10⁴ T is 54× the energy-density field (184.9 T for H–H) and below the Zeeman field (7.8×10⁴ T). Isotope 10⁷ T is 4.8 orders below the smallest nuclear Zeeman field (³H Q, 5.9×10¹¹ T). Nuclear 10¹⁰ T is 3.3 orders below the field matching ⁵⁶Fe binding-energy density (2.05×10¹³ T). Only the first threshold, and only under the energy-density criterion, can be derived. The rest are neither the energy-density nor the µB criterion applied consistently. Strong-field atomic physics (Lai 2001) also says the regimes would strengthen, not break, atomic and molecular bonds.
  5. The reaction p + M → e⁺ + π⁰ + M is quoted correctly from the Callan–Rubakov literature ([PDG](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-mag-monopole-searches.pdf)). The physics is due to Callan (PRD 26, 2058, 1982) and Rubakov (NPB 203, 311, 1982). The talk's "John G. Cramer 1983" is Cramer's Analog column AV-01, written December 1983 and published July 1984 ([npl.washington.edu](https://npl.washington.edu/AV/altvw01.html)). It is a popular account, not a source of the mechanism, and the audit's "Kramer" is a misspelling. Catalysis has no field threshold "B_core ≥ B_crit". The audit's `MonopoleCatalyzedBaryonDecay` yield of 16 × 938.272 = 15012.35 MeV per ¹⁶O overstates the atom's 14899.2 MeV mass-energy by 0.76%.

### Inferences
- If the audit crate is kept as an example, change three things: (i) `fundamental_dirac_charge` returns h/(µ₀e) in A·m, or rename it a flux and drop µ₀ from the field; (ii) use n̂ × r in `vector_potential`; (iii) `B_CRITICAL_EQGP` and the regime thresholds should be labelled "speaker's claim; refuted (lattice QCD)" rather than presented as physics.

### Gaps
- I did not check the audit's other sections (toroidal moments, ¹⁷O reaction Q values, O₂ susceptibility force, ARV phased array). Those are outside this monopole brief.

## 11. Calculators for the crates

### Takeaway
The following functions implement the verified physics. Each has a published or derived test value asserted in monopoles_checks.py. The target is mostly the planned `magnetics::frontier`, which is already specified in the magnetics report but not yet in `crates/magnetics/src`. Field-unit helpers go in the existing `gravitics::strong_field`.

### Cited Findings
Sources as in §1–10.

| Function signature | Formula | Primary source | Published / derived test value | Target crate::module | Status |
|---|---|---|---|---|---|
| `dirac_field(r: f64, n: i32) -> f64` | nh/(4πe r²) = µ₀ n g_D/(4πr²) | PDG eq. 94.1; Castelnovo 2008 (q_D = h/µ₀e) | 1 fm → 3.291×10¹⁴ T; a₀ → 1.175×10⁵ T (derived) | magnetics::frontier | mainstream |
| `dirac_radius_for_field(b: f64, n: i32) -> f64` | √(nħ/(2eB)) | as above | 10¹⁵ T → 0.574 fm; 280 Oe → 108.4 nm; 460 Oe → 84.6 nm (derived) | magnetics::frontier | mainstream |
| `dirac_vector_potential(r: Vec3, n_hat: Vec3, n: i32) -> Vec3` | (µ₀g/4π)(n̂ × r)/(r(r + r·n̂)); string on −n̂ | PDG §94.1 (string); audit corrected | numerical curl at (0.3, −0.4, 0.5) = +µ₀g z/(4πr³) = 4.654×10⁻¹⁶ T (r in m, g = g_D) | magnetics::frontier | mainstream |
| `monopole_field_energy(r0: f64, n: i32) -> f64` | n²ħc/(8α r₀) = µ₀g²/(8πr₀) | derived from PDG eq. 94.1 | 3380.1 MeV·fm; r₀ = 3.602 fm for m_pc²; M_X/(8α) = 17.13 M_X → 1.71×10¹⁷ GeV ∈ PDG 10¹⁷–10¹⁸ | magnetics::frontier | mainstream (derived) |
| `monopole_energy_gain(n: i32, b: f64, l: f64) -> f64` | n g_D B L | PDG eq. 94.8 | 3 µG × 0.526 kpc → 10¹¹ GeV = 16.02 J; 120 Oe force 3.95×10⁻¹¹ N | magnetics::frontier | mainstream |
| `catalysis_length(sigma0: f64, beta: f64, n: f64) -> f64` | β/(σ₀ n) for β ≥ β₀ (IceCube eqs. 4–5) | IceCube 2014 | σ_cat 1.7×10⁻²³ cm², λ 1 cm ⇒ n = 5.88×10²² cm⁻³ (≈ ice H-nucleus density 6.13×10²², 4%) | magnetics::frontier | mainstream (catalysis unobserved) |
| `superk_catalysis_limit(beta: f64, sigma0_mb: f64) -> f64` | 6.3×10⁻²⁴(β/10⁻³)²/(σ₀/1 mb) | Super-K 2012 | (10⁻³, 1) → 6.3×10⁻²⁴; Parker/limit = 1.59×10⁸ | magnetics::frontier | mainstream |
| `field_from_eb(eb_gev2: f64) -> f64` | eB/(eħc²) | Bali 2012; D'Elia 2022; Skokov 2009 | 1 GeV² → 1.690×10¹⁶ T; m_π² (140 MeV) → 3.313×10¹⁴ T; 9 GeV² → 1.52×10¹⁷ T; 10¹⁵ T → 0.0592 GeV² | gravitics::strong_field | mainstream |
| `field_from_natural(b_ev2: f64) -> f64` | B_HL/195.353 eV² per T | MoEDAL 2022/2024; Skokov identity | 7.6 GeV² → 3.89×10¹⁶ T (paper "≈10¹⁶ T"); 4.0 GeV² → 2.05×10¹⁶ T; 1 MeV² → 0.512×10¹⁴ G | gravitics::strong_field | mainstream |
| `lambda4_field(lambda_mev: f64) -> f64` | √(2µ₀) Λ²/(ħc)^{3/2} | audit's criterion; derived | 200 → 2.90×10¹⁴ T; 340 → 8.37×10¹⁴ T; 10¹⁵ T ⇔ Λ = 371.7 MeV | gravitics::strong_field | derived; as a deconfinement threshold refuted (lattice) |
| `zeeman_field(energy_ev: f64, moment: f64) -> f64` and existing `equipartition_field(u)` = √(2µ₀u) | E/µ; √(2µ₀E/V) | CODATA µ_B, µ_N; ladder sources §6 | Ar₂ 212.5 T / 9.65 T; H–H 7.807×10⁴ T / 184.9 T; ⁵⁶Fe B/A 2.788×10¹⁴ T; ²²³Ra→¹⁴C 1.010×10¹⁵ T; m_p 2.976×10¹⁶ T | magnetics::frontier (ladder table as a tested example in physica::analogies) | energies mainstream; "field breaks bond" unverified/contradicted (Lai 2001) |
| `monopole_binds(z: u32, a: u32, spin: f64, mu_n: f64, m: f64) -> Option<f64>` | γ̂ = (A/Z)µ/µ_N, l = mZ; γ̂ > 1 + 1/4l (9.7); \|κ̂\| > (J² + J − l²)/l (9.8); γ_c±(l) for S = 3/2 (9.12), 5/2 (9.13) | Milton 2006 §9; LiveChart moments | γ̂: ¹H 2.79, ²⁷Al 7.56, ⁹Be −2.66, ¹¹³Cd −1.46; γ_c(Be, l = 2) = 1.557/8.943; γ_c(Al, l = 6.5) = 1.374/4.216; ⁵⁷Fe 0.198 → no binding; ¹⁷O −4.02 | magnetics::frontier | contested (theory "primitive", regulator-dependent) |
| `dyonic_level(z: f64, n: u32, j: f64, m: f64) -> f64` | −Z²Ry/(N + ½ + √((j+½)² − m²))² (spinless) | Milton eq. 9.1 | m = 0 → −13.6057 eV; H + Dirac monopole → −7.291 eV (×0.536) | magnetics::frontier | mainstream (spin term omitted) |
| `spin_ice_charge(mu, a)` (existing): add measured vector | 2µ/a_d | Castelnovo 2008; Bramwell 2009 | 4.613 µ_B/Å = 4.278×10⁻¹³ A·m (theory 4.6); measured 5 µ_B/Å = 4.637×10⁻¹³; q_D/q_m = 20r_e/a_d → 7693 | magnetics::frontier | mainstream |
| `synthetic_monopole_circulation(m_atom: f64) -> f64` | 4πħ/M (Dirac-string vorticity flux) | Ray 2014 eq. 2; Hall 2016 (m_Rb) | ⁸⁷Rb: 9.183×10⁻⁹ m²/s | magnetics::frontier | mainstream (synthetic field) |
| `skyrmion_emergent_field(area: f64) -> f64` | (h/e)/A per skyrmion | Neubauer 2009 | 2.5 T ⇔ 1.654×10⁻¹⁵ m² | magnetics::frontier | mainstream (emergent field) |
| `catalysis_yield(atomic_mass_u: f64) -> f64` | M c² of the nucleus/atom (not A·m_p) | AME via LiveChart; PDG catalysis | ¹⁶O 14899.2 MeV (audit's 15012.35 too high by 0.76%); m_p/D–D: 257 (per reaction), 1028 (per nucleon), 157.4 (D+D→⁴He per nucleon) — talk's 166 unreproduced | plasmoid::fusion | mainstream for GUT monopoles; talk's 166 unverified |
| `urutskoev_monopole_count(dh_oe: f64, area: f64, depth: f64) -> f64` (example, not library) | V/(4π r(ΔH)³/3) | Ivoilov & Urutskoev 2004 | 280 Oe, 15×15 mm, 0.2 µm → 8.43×10⁹; 1.50 kg if GUT; 202.6 J at 75 GeV pairs; 280 Oe ≡ 7 K | magnetics::frontier example | unverified (claim) |

### Inferences
- Under the user's "rigs are examples" rule, the Urutskoev count and the bond ladder belong in examples or tests, not library functions. Formulas live once in the libraries.
- `dirac_charge`, `DIRAC_FLUX`, `flux_from_events`, `parker_bound` and `spin_ice_charge` already exist in the specification. `dirac_field`, `dirac_radius_for_field` and `monopole_field_energy` should call `dirac_charge` rather than restate h/(µ₀e).

### Gaps
- `crates/magnetics/src` has no `frontier.rs` yet (checked 2026-10-26), so the module path is the planned one. No Rust was written (out of scope).

## Primary-source check (2026-10-27)

Before the remaining calculators were coded (#24), each formula above that came from a review was read in its original where one could be opened; see [monopole_primary_sources.md](monopole_primary_sources.md) and its check script. Everything used here matches. Corrections:
- **Milton 2006, eq. 9.3** implies γ = µ/(2Sµ_N), but his Table 1 and the criteria (Olaussen & Sollie's g_M = Aµ/(SZ)) need **γ̂ = (A/Z)µ/µ_N**, as used here. They differ by 2S for S > ½.
- **Milton §9.2's "1.67 < γ_c⁺" for ²⁷Al** should read 3 < γ_c⁺ ≤ 4.216, the large-l limit of his eq. 9.13. Al still binds.
- **Bracci & Fiorentini, Table I**: the ³H λ_I is printed 2.97; the formula gives 3.97. Not used here.
- **Super-Kamiokande 2012** is Astropart. Phys. 36, 131 (not PRD). Eq. 9 carries a factor (f_π⁺/0.5) and holds for β < 10⁻³ and M ≲ 10¹⁷ GeV. `superk_limit` above now has the factor.
- **IceCube 2014, eq. 4** has a second branch below β₀, σ₀/β·(β/β₀)^γ; the catalysis-length row above covers only β ≥ β₀.
- **PDG eq. 94.8** is the galactic speed, not an energy gain. `monopole_energy_gain` is derived from the field, and the 0.526 kpc row is our inference.
- **Λ_QCD**: PDG editions from 2020 to 2025 print only α_s. The last printed Λ is PDG 2018 eq. 9.24, with Λ⁽³⁾ = 332 ± 17 MeV.
- **Ahlen 1978 and Ahlen & Kinoshita 1982** remain unopened (paywalled). The slow factor 1.37 and the K and B constants rest on Patrizii & Spurio and Milton §8.
- **Ahlen & Kinoshita's free-gas formula** turns negative below n_e = 1.28×10²⁹ m⁻³ (v_F < 0.82αc), so the crate returns `None` there.
