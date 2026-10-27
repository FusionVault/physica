# Transmutation claims in Greenyer's "The Power and Signature of the Topological Monopole" (TeslaTech, 15 Aug 2026): primary sources, and the nuclear and energy accounting for each

Scope: the nine slide items, the talk's other transmutation claims (transcript lines ~205–235, 310–330, 730–745, 1290–1330, 1580–1760, 2500–2660 of `/home/dan/Desktop/transcript_1793004282390.txt`), and the user's own audit of the talk (`research/external/2026-10-26-greenyer-audit-and-crate.md`). It builds on `research/topics/coulomb-and-schwinger/notes/fringe_barrier_claims.md` (DOE 2004, Google 2019, screening, Widom–Larsen, Rossi, MFMP's critique) and does not repeat it.

Method. Masses are AME2020 atomic mass excesses ([AME2020 mass_1.mas20.txt, IAEA AMDC](https://www-nds.iaea.org/amdc/ame2020/mass_1.mas20.txt); Wang, Huang, Kondev, Audi, Naimi, *Chinese Physics C* 45, 030003 (2021), [doi:10.1088/1674-1137/abddaf](https://doi.org/10.1088/1674-1137/abddaf), Crossref-checked), 42 nuclides embedded in `transmutation_claims_checks.py`. Q = ΣΔ(in) − ΣΔ(out) on atomic masses, with A and Z conserved (the script refuses an unbalanced reaction); electron capture (EC) and β⁻ are carried by the atomic masses, each emitted positron costs 2mₑc² of kinetic Q (it comes back when the positron annihilates). Isotopic abundances are NIST's ([NIST Atomic Weights and Isotopic Compositions](https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=ascii2&isotype=some)). Tunnelling uses the formulas of `crates/plasmoid/src/fusion.rs` (`gamow_energy` E_G = 2m_rc²(παZ₁Z₂)², `gamow_peak` E₀ = (E_G kT²/4)^⅓, `tunnelling` e^(−√(E_G/E))); the Maxwell-averaged factor is e^(−3E₀/kT). Every number marked *(computed here)* is asserted in the check script (98 checks, 0.02 s, `python3 -I`). Labels: **mainstream** (textbook physics), **contested** (published, open dispute), **unverified** (claimed, no independent check found), **refuted** (independently tested and not confirmed, or contradicts a measured bound by a large factor).

## 1. Joseph Papp (1960s): the Papp noble-gas engine

### Takeaway
Papp's two patents describe a noble-gas "fuel" fired by sparks and radioactive electrodes, attribute the energy to nuclear reactions, and give no measured output; no working engine has been shown to anyone outside his circle. Label: **unverified** (no measurement exists to test), and any nuclear release would carry a radiation signature that the patents do not report.

### Cited Findings
- US 3,680,431 (filed 1 Nov 1968, granted 1 Aug 1972), "Method and means for generating explosive forces": a fill of roughly 40–60% argon, 30–40% xenon, 6–8% neon at no more than about 3–4 atm; the abstract says the reaction involves "fission, spallation, transfer reactions and cascade effects" and the mechanism "is not completely understood"; the patent gives no energy per charge, no power and no yield measurement — [US3680431A, Google Patents](https://patents.google.com/patent/US3680431A/en)
- US 4,428,193 (filed 4 Sep 1980, granted 31 Jan 1984), "Inert gas fuel, fuel preparation apparatus and system for extracting useful work from the fuel": preferred fuel 35.6% He, 26.3% Ne, 16.9% Ar, 12.7% Kr, 8.5% Xe at one atmosphere per cylinder; electrodes holding about 4 g thorium-232 and about 2 g "rubidium-37" (the patent labels elements by atomic number); presented as a thermonuclear process with an ignition "about 100,000,000 °F" lasting about 10⁻⁶ s; it states no measured output power, efficiency or horsepower, and no run time — [US4428193A, Google Patents](https://patents.google.com/patent/US4428193A/en)
- The 1966 (or 1968) Caltech demonstration that Feynman attended, the explosion and the settled lawsuit are documented only in secondary accounts and Feynman's own later write-up; Wikipedia notes no contemporary confirmation of the demonstration — search summaries of [Josef Papp, Wikipedia](https://en.wikipedia.org/wiki/Josef_Papp) and [Museum of Hoaxes reprint of the 1989 San Jose Mercury News article](https://hoaxes.org/comments/papparticle.html) (not opened in full)

### Inferences
- Noble gases have no chemistry to release, so any work beyond the electrical input must be nuclear. A 100 cm³ cylinder at 1 atm and 293 K holds 4.2×10⁻³ mol (2.5×10²¹ atoms) *(arithmetic, not scripted)*; the energy content is not the limit if the process were nuclear — the radiation is. Every nuclear channel at the MeV scale emits gammas or neutrons in proportion to power (for D–D, `plasmoid::fusion::dd_per_watt` gives 1.7×10¹² reactions/s per watt, half of them with a neutron); the patents report radioactive sources as inputs, not emissions as outputs.
- The 10⁸ °F ignition (≈ 5.6×10⁷ K, kT ≈ 4.8 keV) for 1 µs in 1 atm of gas is the patent's assertion; at that kT noble-gas nuclei (Z ≥ 2) have negligible fusion reactivity, so the "fusion" route would need the D–T-class tunnelling that only hydrogen isotopes have.
- Both ways: the claim is tested by a calorimetric run with neutron and gamma monitors; none is published.

### Gaps
- No measured horsepower, run time or radiation figure exists in either patent; I found no independent test report.

## 2. John Hutchison (1979+): disruption, transmutation, levitation

### Takeaway
The repository already refutes the levitation reading: Hutchison's 27 kg ball needs 9.65 kPa against at most 40 Pa of electrostatic pressure at breakdown (`research/topics/magnetics-and-mhd/report.md`, `frontier` table, label refuted). George Hathaway's sample analyses exist only as secondhand descriptions. Label: **unverified** (transmutation), **refuted** (levitation by the electrostatic forces available).

### Cited Findings
- The repository's calculation: "Hutchison's 27 kg ball needs 9.65 kPa against 40 Pa electrostatic at breakdown" — `research/topics/magnetics-and-mhd/report.md` (row `drain_time`, `reactionless_power`)
- Hathaway's spectrographic results on Hutchison's metal samples are described only in a podcast blurb ("considerable anomalies in the metallic crystalline structure and distribution of naturally occurring impurities") and a 1996 Planetary Association for Clean Energy compilation that lists Hathaway among assessments — search summary of [Tim Ventura Interviews episode](https://podcast.app/hutchison-effect-metallurgy-spectrographic-analysis-e371145643) (not opened)

### Inferences
- "Anomalies in the distribution of impurities" are not by themselves transmutation: segregation of impurities to grain boundaries and fracture surfaces is ordinary metallurgy. A transmutation claim needs a mass balance or a non-natural isotope ratio (see section 18).

### Gaps
- No primary Hathaway report, its method, its controls or its data were found; I could not open the 1996 compilation.

## 3. Ken Shoulders (1980s+): EVs boring through glass and alumina

### Takeaway
Shoulders' patent describes micrometre "EVs" of ~10¹¹ electrons with little ion neutralisation; an unneutralised cluster of that charge and size would sit at ~288 MV, five orders above the ~2 kV used to launch it, so the EV must be almost fully neutralised or the stated charge is an overestimate. Label: **unverified** (the boring is documented by Shoulders' micrographs; the charge figure is internally inconsistent with the launch voltage unless neutralised).

### Cited Findings
- US 5,018,180 (priority 3 May 1988, granted 21 May 1991), "Energy conversion using high charge density": EVs typically about 1 µm across (0.1 µm to larger), "on the order of 10¹¹ electron charges" with little positive-ion compensation, rings of ten beads ~10¹² charges up to ~20 µm across moving at roughly c/10, launched by cathode pulses around 2 kV (some designs under 100 V); a passing EV "may cut into the guide material" — [US5018180A, Google Patents](https://patents.google.com/patent/US5018180A/en)

### Inferences
- Surface potential kQ/R for 10¹¹ e in a 1 µm sphere: 288 MV *(computed here)*; self-energy (3/5)QV ≈ 2.8 J, i.e. ~170 MeV per electron. A 2 kV pulse can supply at most 2 keV per electron. So either the cluster is neutralised to about 1 part in 10⁵ (then it is a dense neutral plasma, not a bare charge) or the 10¹¹ figure is an overestimate. The existing `plasmoid::fusion::charged_sphere(q, r)` already returns these quantities.
- Boring through glass needs only local energy deposition of order the vaporisation energy of the channel, which any dense arc spot also delivers; it is not a nuclear signature.

### Gaps
- Shoulders' other patents (5,054,046; 5,123,039; 5,054,047; 5,148,461) were not opened.

## 4. M. Solin (1992): "superconducting condensate" / "quantum coherent nuclear reactor" patent

### Takeaway
I could not find the 1992 Solin patent. A secondary history says Mikhail Solin's electron-beam experiments on molten zirconium (1976) were "patented in the post soviet era panic in 1992". Label: **unverified** (primary source not located).

### Cited Findings
- Mikhail Solin's electron-beam work on molten zirconium dates to 1976 and was patented in 1992; no number, title or claims are given — [Remote View, "A Century of Cold Fusion in Russia: Ekaterina Belousova and Alexander Parkhomov"](https://remoteview.substack.com/p/a-century-of-cold-fusion-in-russia)
- A Russian review lists Solin as leading a group at JVA Corporation, Ekaterinburg, and a bibliography lists a Solin paper with V. G. Khlopin Radium Institute co-authors on "D-D Fusion Stimulation by a High-Current Arc Discharge in Gas-Filled Metal" — search summaries of [Goryachev, LENR-CANR](https://lenr-canr.org/acrobat/GoryachevIorganizati.pdf) and [New Energy Times cold fusion bibliography](https://www.newenergytimes.com/v2/archives/fic/BIB_R11/FBIB_3N4.PDF) (not opened)

### Gaps
- Patent number, office (likely RU), claims, and any measurement. A search of FIPS/Espacenet for inventor "Solin" (Солин) 1992 is the next step; Google Patents' inventor search needs a browser.

## 5. Mahadevan Srinivasan (ICCF-20, Sendai, 2016): "four tonnes of excess metal per day for 11 weeks" in an 11 MW carbon-arc furnace

### Takeaway
I found no primary record (paper, abstract or slides) of the furnace claim. The only trace is a forum tribute, and it is ambiguous between "4 tonnes over 11 weeks" and "4 tons/day" — a factor of 77. Taken as stated in the talk (4 t/day of iron via Si+Si or C+O fusion), the fusion power is 0.87–4.1 TW, 8×10⁴–4×10⁵ times the furnace's 11 MW, with a steady ⁵⁶Ni/⁵⁶Co activity of ~10¹³ Ci; an 11 MW ferrosilicon furnace makes only ~29–31 t of alloy a day. Label: **unverified** (no primary source); as a nuclear claim at the stated rate, **refuted** by the energy and radiation budget; as an accounting anomaly of 52 kg/day, it is within ordinary mass-balance uncertainty.

### Cited Findings
- A LENR Forum tribute says Srinivasan revealed at ICCF-20 in Sendai that an industrial arc furnace in India had produced 4 tonnes of excess metals (Si, Fe) over 11 weeks with no ill effects on workers, likened it to Ohsawa's carbon-arc reactions and MFMP's NOVA experiments, and gives his death as 31 August 2020; a separate post in the same thread says "4 tons/day" — search-result summary of [RIP Dr M. 'Chino' Srinavasan, LENR Forum](https://www.lenr-forum.com/forum/thread/6373-rip-dr-m-chino-srinavasan/) (the page returned HTTP 403; not read directly)
- ICCF-20 was held in Sendai, 2–7 October 2016, sponsored by Clean Planet and Tohoku University — search summary of [LENR-CANR, ICCF20](https://lenr-canr.org/wordpress/?p=1752) (not opened)
- The talk itself: "four tons of excess metal synthesized per day for 11 weeks continuously in a carbon arc furnace in Coimbatore", "11 MW three-phase", "fusion of carbon and oxygen to silicon, fusion of two silicons to iron" — transcript lines 215–221, 1675–1684
- FeSi75 in a submerged arc furnace takes 8.46–9.18 MWh per tonne of alloy (seven reductant mixes; one plant makes about 75 t/day at 55.2 MVA) — [Etemadi, Koohestani & Tajally, *Heliyon* 9, e13956 (2023), doi:10.1016/j.heliyon.2023.e13956](https://pmc.ncbi.nlm.nih.gov/articles/PMC9984837/)
- Silicon burning in stars proceeds by photodisintegration and α-capture up to ⁵⁶Ni at 2.7–3.5 GK, and ⁵⁶Ni decays to ⁵⁶Co and then ⁵⁶Fe with half-lives of about 6 and 77 days — [Silicon-burning process, Wikipedia](https://en.wikipedia.org/wiki/Silicon-burning_process) (citing Clayton 1983 and Woosley, Heger & Weaver 2002)

### Inferences
- Q-values *(computed here)*: ¹²C+¹⁶O→²⁸Si +16.756 MeV; ²⁸Si+²⁸Si→⁵⁶Ni +10.922 MeV; ⁵⁶Ni→⁵⁶Fe by two EC +6.700 MeV (⁵⁶Ni→⁵⁶Co 2.133, ⁵⁶Co→⁵⁶Fe 4.567). Z balance: 2 Si (Z = 28) give ⁵⁶Ni, not ⁵⁶Fe (Z = 26); the iron needs two weak decays. Per ⁵⁶Fe from 2(C+O): 51.133 MeV.
- 4 t/day of ⁵⁶Fe = 4.31×10²⁸ atoms/day. Via Si+Si only: 8.72×10¹¹ W = 7.9×10⁴ × 11 MW; from C+O all the way: 4.08×10¹² W. If the iron passed through ⁵⁶Ni, the steady activity equals the production rate, 5.0×10²³ Bq = 1.35×10¹³ Ci *(computed here)*. At "4 t over 11 weeks" the Si+Si power is still 1.13×10¹⁰ W, 10³ × the furnace.
- Tunnelling at arc temperatures: Maxwell-averaged factor for ¹²C+¹⁶O is 10^−2708 at 5000 K and 10^−2149 at 10⁴ K; for ²⁸Si+²⁸Si 10^−8773 at 5000 K *(computed here)*. For comparison the measured ¹²C+¹⁶O barrier is ~7.8 MeV ([Goasduff et al., arXiv:1107.1124](https://ar5iv.arxiv.org/html/1107.1124)); arc kT is ~0.4–0.9 eV.
- Both ways: an 11 MW furnace at 8.46–9.18 MWh/t makes 28.8–31.2 t/day of FeSi75 *(computed here)*. "4 t/day excess" is 13–14% of that output; "4 t in 11 weeks" (52 kg/day) is 0.17%. A 0.2% discrepancy is within the assay and weighing uncertainty of raw materials (quartzite, coke, scrap, electrode paste); a 13% one would be a gross accounting error. Without the mass balance, assays and isotope data, neither can be tested, and the "no ill effects" observation is itself the strongest evidence against a nuclear origin at the stated rate.
- A nuclear origin via α-conjugate fusion (C+O, Si+Si) would make essentially pure ⁵⁶Fe and ²⁸Si; natural iron is 91.75% ⁵⁶Fe, 5.85% ⁵⁴Fe, 2.12% ⁵⁷Fe, natural silicon 92.2% ²⁸Si ([NIST](https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=ascii2&isotype=some)). An ⁵⁴Fe/⁵⁶Fe ratio diluted below 0.064 in the "excess" metal would be a signature; natural ratios would point to unaccounted feed.

### Gaps
- No ICCF-20 paper, abstract or slide by Srinivasan on this furnace was found; the furnace's product (FeSi? SiC? CaC₂?), the year, the mass-balance method and any isotope data are unknown. The forum page could not be read directly (403) and is itself secondhand.

## 6. Takaaki Matsumoto (1989–2001): the "Nattoh" model, "itons", "electro-nuclear collapse", carbon from lead

### Takeaway
Matsumoto's Fusion Technology notes (1989–1993) propose deuteron clusters ("Nattoh"), new particles ("itons"), and traces on nuclear emulsions and copper plates that he attributes to neutron "gravity decay" and tiny black and white holes; none gives a nuclear mass balance. Lead → carbon is endothermic by 24 MeV per ²⁰⁸Pb even with 22 β⁻ decays. Label: **unverified** (no independent replication found); the black/white hole interpretation is not testable from the published abstracts.

### Cited Findings
- "'Nattoh' Model for Cold Fusion", *Fusion Technology* 16(4), 532–534 (Dec 1989), doi:10.13182/FST89-A29119: small deuteron clusters and many-body fusion — [ANS](https://wx1.ans.org/pubs/journals/fst/a_29119)
- "Cold Fusion Experiments with Ordinary Water and Thin Nickel Foil", *Fusion Technology* 24(3), 296–306 (Nov 1993), doi:10.13182/FST93-A30205: excess heat from ordinary water, emulsion tracks of electrons, protons and deuterons, "micro-explosions" from neutron-nucleus decay, traces attributed to tiny black and white holes — [ANS](https://wx1.ans.org/pubs/journals/fst/a_30205)
- "Experiments of One-Point Cold Fusion", *Fusion Technology* 24(3), 332–339 (Nov 1993), doi:10.13182/FST93-A30209: discharges in H₂O/D₂O with 0.6 mol/L K₂CO₃, Pt pin anode on a Cu plate cathode, ring-shaped spots attributed to "gravity decay" of single and di-neutrons and "itonic hydrogen clusters" — [ANS](https://wx1.ans.org/pubs/journals/fst/a_30209)
- A December 1990 paper predicts the iton with a rest mass 2 to 26 electron masses; a 1990s letter to the editor questions whether the one-point traces are charge clusters rather than fusion — search summary of [ANS Fusion Technology listings](https://wx1.ans.org/pubs/journals/fst/a_29308) (not opened; the letter's PDF is a scan that could not be read)
- The talk: "carbon synthesis from the lead … ring products consisted of conventional elements mainly carbon not dependent on the plasma material" (transcript lines 1531–1537); Crossref also lists "Progresses in NATTOH model and new particles emitted during cold fusion", AIP Conf. Proc. 228, 827–839 (1991), doi:10.1063/1.40720

### Inferences
- ²⁰⁸Pb → 17 ¹²C + ⁴He (A 208, Z 82 → 104 needs 22 β⁻): Q = −24.17 MeV *(computed here)*; the energy must come from the discharge. Carbon is also the most common surface contaminant in SEM/EDX work (pump oil, adhesive tape, adventitious hydrocarbons), so "mainly carbon whatever the plasma material" is the result contamination predicts.
- A charge-cluster reading (Shoulders' EVs) of the ring spots is the critique raised in the journal itself; it needs no nuclear process.

### Gaps
- Matsumoto's 2001 book was not found; the "electro-nuclear collapse" claim could not be read in the original.

## 7. A. F. Kladov (2003): cavitation transmutation / destruction of matter

### Takeaway
Kladov's Russian patent RU 2054604 ("Method for producing energy", filed 2 Jul 1993, published 20 Feb 1996) claims bubble-collapse energy "several orders of magnitude higher than the binding energy" of nucleons; I found only secondhand quotes of the 2003 "Cavitation destruction of matter" text. Label: **unverified**; the patent's mechanism contradicts the measured bound on cavitation temperatures.

### Cited Findings
- RU 2054604 by A. F. Kladov, with the stated mechanism that the energy released in bubble collapse is "several orders of magnitude higher than the binding energy of elementary particles (nucleons) in the nucleus"; the 2003 work lists three modes (adding nucleons, splitting into fragments, breaking into the smallest fragments or radiation) — search summaries of [Remote View, "Method for obtaining energy, Anatoly Kladov"](https://remoteview.substack.com/p/method-for-obtaining-energy-anatoly) and [Remote View, "O-Day – Blowing Bubbles"](https://remoteview.substack.com/p/o-day-blowing-bubbles) (blogs; not opened)

### Inferences
- Nucleon binding is ~8 MeV per nucleon; single-bubble sonoluminescence temperatures are of order 10⁴ K (kT ≈ 1 eV), seven orders of magnitude lower. At 10⁴ K the ¹⁶O+¹⁶O Maxwell-averaged tunnelling factor is 10^−2740 *(computed here)*. The sonofusion record (Taleyarkhan) is covered in `fringe_barrier_claims.md` Q3.

### Gaps
- The patent text (FIPS) and the 2003 report were not opened; no isotope or mass data were found.

## 8. S. V. Adamenko (Proton-21, Kyiv, 2000–2006; book 2007): copper → lead, gold and "every element"

### Takeaway
Adamenko's group reports >4000 (later 10,000–30,000) relativistic-electron-beam shots of up to 1 kJ (100–300 J typical) on metal anodes, with products from across the periodic table, non-natural isotope ratios (S, Ti, Fe; ⁵⁷Fe > ⁵⁶Fe) and no residual radioactivity. Making lead from copper is endothermic by ≥ 267 MeV per ²⁰⁸Pb, so 300 J could make at most 7×10¹² Pb atoms (2.4 ng) even at 100% efficiency, against the lab's claim that new nuclei are "comparable with the total number of initial copper nuclei". Label: **unverified** (no independent peer-reviewed replication found; the review that calls the evidence "impressive" is a book review).

### Cited Findings
- Adamenko & Adamenko, "Analysis of laboratory nucleosynthesis products", arXiv:nucl-ex/0307011 (2003): more than 4000 shots on light, medium and heavy targets, total energy up to 1 kJ; "Elements with nuclei three and more times heavier than the nucleus of the target main element are detected in the products"; isotope ratios differ significantly from natural — [arXiv](https://arxiv.org/abs/nucl-ex/0307011)
- Book: *Controlled Nucleosynthesis: Breakthroughs in Experiment and Theory*, eds. Adamenko, Selleri, van der Merwe, Springer 2007, Fundamental Theories of Physics 156, doi:10.1007/978-1-4020-5874-5 (Crossref-checked). Dolan's review: experiments 2000–2006; a fast REB on a hemispherical anode tip (Cu, Fe, Pb…) explodes it; ion energy from optical spectra ~2 kJ "much greater than the REB energy of 300 J"; implanted ⁶⁰Co or ^110m^Ag activity typically 5–10% lower after the shot; "S, Ti, and Fe, showed significant shifts from their natural isotopic abundances"; superheavy (A > 300) claims; proposed "Coulomb collapse" at n_e > n_cr ~ 2×10³² cm⁻³ for Z = 92; 773 pages; "the experimental evidence is impressive" — [Dolan, New Energy Times #32 (2009)](https://www.newenergytimes.com/v2/books/Reviews/AdamenkoByDolan.pdf)
- Kowalski's compilation: input "100-200 joules in a beam of electrons"; "The total number of new nuclei is comparable with the total number of initial copper nuclei"; Fe-57 exceeding Fe-56 in products; product radioactivity "not different from the background activity"; confirmation in US and Russian labs asserted but not cited; a Ukrainian patent and a PCT application — [Kowalski, Montclair State University page 217](https://msuweb.montclair.edu/~kowalskil/cf/217kiev.html)
- US patent application 2005/0200256 (filed 19 May 2003, published 15 Sep 2005): compressing a substance by impact in axisymmetric relativistic vacuum diodes with a plasma cathode and an anode-enhancer — search summary of [Justia, inventor Stanislav Adamenko](https://patents.justia.com/inventor/stanislav-adamenko) (page returned 403)
- The talk: "a copper pin, hits it with 300 J … nearly 100% lead in a blob … 83% gold"; "the awarded 2023 patent … an iron hemisphere comprising of 93% by weight iron and admixtures of silicon" (transcript lines 317–323, 735–738)

### Inferences
- Nucleon budget *(computed here)*: ²⁰⁸Pb (Z 82, N 126) is more neutron-rich than copper (N/Z 1.17–1.24), so 4 ⁶³Cu (A 252, Z 116) → ²⁰⁸Pb leaves 34 p + 10 n. With strong processes only (spares as 5 α + 24 p): Q = −427.6 MeV; the least endothermic route (spares as 11 α after 12 EC): −267.2 MeV. Gold: 3 ⁶³Cu + ⁶⁵Cu → ¹⁹⁷Au + 10 α + 17 p: −381.0 MeV; via 9 EC to ¹⁹⁷Au + 14 α + n: −274.9 MeV. Binding per nucleon peaks near Fe/Ni, so every Cu → Pb/Au path costs energy.
- At 300 J and −267 MeV per Pb: ≤ 7.0×10¹² atoms = 2.4 ng per shot, a lead sphere of at most 7.4 µm diameter *(computed here)*, with 100% of the beam energy going into the transmutation. A µm-scale "blob of nearly 100% lead" seen in SEM is therefore not excluded by energy alone; "new nuclei comparable to the copper nuclei" for any target heavier than a few ng is. Both ways: microgram-scale lead and gold are routinely present as contamination (solders, electrodes, sputter coatings); the iron hemisphere of 93% Fe with Si is the composition of steel/iron debris.
- Lower ⁶⁰Co activity by 5–10% after the shot is what loss of implanted material from an exploding anode predicts; the claim needs total-activity accounting of all debris.

### Gaps
- No 2023 Adamenko patent was found (the talk's "awarded 2023 patent"); the book's product masses per shot and its SIMS isotope tables were not opened.

## 9. Bogdanovich et al. (2019): video of long-lived plasmoids

### Takeaway
Covered in depth by another researcher; bibliographic check only. Label: **unverified** (interpretation).

### Cited Findings
- B. Yu. Bogdanovich, N. V. Volkov, N. A. Len', A. V. Nesterovich (MEPhI), "Video Recording of Long-Lived Plasmoids near Objects Exposed to Remote and Direct Effects of High-Current Pinch Discharges", *Technical Physics* 64(4), 465–469 (2019), doi:10.1134/S1063784219040066 (Crossref-checked); luminous formations persisting after periodic discharges in a water flow, discussed with reference to Dirac's monopole — search summary of [MEPhI repository](https://openrepository.mephi.ru/entities/publication/8129ceef-e7d5-4d4c-a233-07733cd4707e)

## 10. Ohsawa (1964) and Kervran: "biological transmutation" and carbon-arc iron

### Takeaway
Kervran's reactions (Na+O→K, K+H→Ca, Mg+O→Ca) are each exothermic by 8–20 MeV per event, so they would emit MeV gammas or particles; at body temperature their tunnelling factors are ~10^−11 400. Ohsawa and Kushi's carbon-arc "iron" was identified by magnetism; even a macrobiotic publication notes magnetic carbon is possible and that the Na→K experiment was never reproduced. Label: **refuted** (biological transmutation at the stated rates), **unverified** (Ohsawa arc iron; see section 11 for the 1994 tests).

### Cited Findings
- "Unfortunately, the first experiment [Na → K] was never reproduced by others, despite numerous attempts"; the conclusion of Ohsawa and Kushi that iron had formed rested on magnetism, which "alone is therefore no proof at all" because magnetic forms of carbon exist — search summary of [Ohsawa Macrobiotics, transmutation article (2014)](https://ohsawamacrobiotics.com/wp-content/uploads/2022/05/transmutation-mg.pdf) (not opened)
- Kervran's ideas shaped Ohsawa's macrobiotic diet; his claims contradict known physics — search summary of [Corentin Louis Kervran, Wikipedia](https://en.wikipedia.org/wiki/Corentin_Louis_Kervran) (not opened)

### Inferences
- *(computed here)*: ²³Na+¹⁶O→³⁹K +19.540 MeV; ³⁹K+¹H→⁴⁰Ca +8.328 MeV; ²⁴Mg+¹⁶O→⁴⁰Ca +16.176 MeV. ²³Na+¹⁶O at 310 K: Maxwell-averaged tunnelling 10^−11 397. A body transmuting even 1 mg/day of Na to K would release 1 mg/23 g mol⁻¹ × N_A × 19.5 MeV ≈ 82 MJ/day ≈ 950 W *(arithmetic)* as penetrating radiation.
- 2 ¹²C + 2 ¹⁶O balances to A 56 but Z 28: the product is ⁵⁶Ni, radioactive (6 d), not ⁵⁶Fe; see section 11.

### Gaps
- Ohsawa's 1964 primary report was not found.

## 11. Sundaresan & Bockris (Texas A&M, 1994) and Singh et al. (BARC, 1994): iron from carbon arcs in water

### Takeaway
Sundaresan & Bockris saw iron rise from ~2 ppm in the rods to up to 286 ppm in the detritus, none when N₂ replaced dissolved O₂; Singh et al. saw 50–2000 ppm, dropping to 20–100 ppm when the trough was covered, with natural iron isotope ratios and Si, Ni, Al, Cr following the same pattern, and called it unclear whether iron was synthesised. Natural isotope ratios and co-varying Ni/Cr/Al/Si are the contamination signature; α-conjugate synthesis (C+O, 4 ¹⁷O) would give pure ⁵⁶Fe. Label: **contested** (published, unresolved), leaning to contamination on the isotope evidence.

### Cited Findings
- R. Sundaresan, J. O'M. Bockris, "Anomalous Reactions During Arcing Between Carbon Rods in Water", *Fusion Technology* 26(3P1), 261–265 (Nov 1994), doi:10.13182/FST94-A30330 (Crossref-checked): spectroscopically pure carbon rods, highly purified water, 20–25 A intermittently for several hours; iron from about 2 ppm to as much as 286 ppm in the detritus; rods cool beyond 2 cm from the tips; adsorption from water or air ruled out; weak correlation with arcing time; no iron when nitrogen replaced dissolved oxygen; temperature rise said to match the proposed reaction — [ANS abstract](https://wx1.ans.org/pubs/journals/fst/a_30330)
- M. Singh, M. D. Saksena, V. S. Dixit, V. B. Kartha, "Verification of the George Oshawa Experiment for Anomalous Production of Iron from Carbon Arc in Water", *Fusion Technology* 26(3P1), 266–270 (Nov 1994), doi:10.13182/FST94-A30331 (Crossref-checked): DC arc 1–20 h, iron 50–2000 ppm in the first series, 20–100 ppm with the trough fully covered; isotopic abundances matched natural iron; Si, Ni, Al, Cr followed similar patterns; whether iron "is really being synthesized through transmutation" left open — [ANS abstract](https://wx1.ans.org/pubs/journals/fst/a_30331)
- The talk: "if they boiled the water and put nitrogen there, they saw zero synthesis of iron. If they put oxygen in there, they saw the iron synthesize"; "They never looked for helium" (transcript lines 1665–1671)

### Inferences
- Energy per µg of ⁵⁶Fe from 2(C+O): 88.2 kJ *(computed here)*; a 25 A × ~25 V arc (~600 W) for 5 h delivers ~11 MJ, the nuclear heat of ~120 µg of iron. Whether the detritus held more or less than that is not in the abstract (gap). If it passed through ⁵⁶Ni, 1 µg is 0.39 Ci of ⁵⁶Ni *(computed here)* plus a 77-day ⁵⁶Co tail — trivially detectable with a survey meter; no radiation is reported.
- Both ways on the O₂/N₂ result: dissolved oxygen also oxidises and dissolves iron-bearing parts (clamps, electrode holders, vessel fittings) and changes arc chemistry; covering the trough cut the iron 10–20× in Singh et al., which points to airborne or splash contamination, not oxygen as a nuclear fuel. The covaried Ni, Cr, Al, Si is the fingerprint of stainless steel and dust.
- Natural isotope ratios (Singh) are the decisive test: 2(¹²C+¹⁶O), 2 ²⁸Si and 4 ¹⁷O all make only A = 56, so synthesised iron would lack ⁵⁴Fe and ⁵⁷Fe (5.85% and 2.12% of natural iron).

### Gaps
- Full texts (detritus mass, absolute µg of Fe, calorimetry) were not opened (paywalled).

## 12. The speaker's reaction 4 ¹⁷O → ¹²C + ⁵⁶Fe

### Takeaway
The reaction balances (A 68, Z 32) and is exothermic by 57.37 MeV *(computed here)*; it is a five-body rearrangement with no known mechanism, ¹⁷O is 0.038% of oxygen, and it would make isotopically pure ⁵⁶Fe. Label: **unverified** (proposed by the speaker; no measurement of ¹⁷O depletion or pure ⁵⁶Fe exists).

### Cited Findings
- The talk: "You get oxygen 17 formed first and then the four oxygen 17s go to carbon 12 and iron 56" (transcript lines 1663–1665)
- Natural oxygen: ¹⁶O 99.757%, ¹⁷O 0.038%, ¹⁸O 0.205% — [NIST](https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=ascii2&isotype=some)

### Inferences
- Even the pairwise first step ¹⁷O+¹⁷O has E_G = 34.1 GeV; at 10⁴ K and with 4-body coincidence required, no rate exists. If ¹⁷O must first be made from ¹⁶O by p + e⁻ + ¹⁶O → ¹⁷O + ν (+3.361 MeV, a weak process), the rate is that of electron capture in flight, not a strong-interaction rate.
- At Srinivasan's 4 t Fe/day, natural-abundance ¹⁷O would require 1.2×10⁴ t of oxygen per day *(computed here)*.

## 13. ¹⁶O + ¹⁶O fusion, "millions of eV per fusion of two oxygens"

### Takeaway
Correct as nuclear data: the ¹⁶O+¹⁶O channels release up to 16.54 MeV (³²S+γ) and mostly 7.7–9.6 MeV (³¹P+p, ²⁸Si+α), but only at stellar oxygen-burning temperatures of 1.5–2.6 GK; 5% of events give neutrons. Label: **mainstream** (reaction energies), **refuted** (as an arc or plasmoid process at ≤10⁴ K).

### Cited Findings
- Branches and Q: ²⁸Si+α +9.593 MeV (34%), ³¹P+p +7.676 MeV (56%), ³¹S+n +1.459 MeV (5%), ³⁰Si+2p +0.381 MeV, ³⁰P+d −2.409 MeV (5%); ³²S+γ +16.539 MeV; ²⁴Mg+2α −0.393 MeV; ignition at (1.5–2.6)×10⁹ K — [Oxygen-burning process, Wikipedia](https://en.wikipedia.org/wiki/Oxygen-burning_process) (citing Woosley, Heger & Weaver, *Rev. Mod. Phys.* 74 (2002), and Clayton 1983)

### Inferences
- AME2020 reproduces every listed Q within 3 keV except ³¹S+n: 1.497 MeV against the listed 1.459 (38 keV) *(computed here)*; the listed figure probably comes from an older ³¹S mass; the check uses a stated 40 keV tolerance for it.
- At 10⁴ K the Maxwell-averaged tunnelling factor is 10^−2740 *(computed here)*. The 5% neutron branch means oxygen burning makes 4.2×10¹⁰ n/s per watt *(computed here, branch-weighted mean Q 7.51 MeV)*, which no plasmoid experiment has reported.

## 14. "938 MeV per proton … 16 times for oxygen … 166 times more than D–D" (monopole-catalysed baryon collapse)

### Takeaway
Monopole catalysis (Rubakov–Callan) of nucleon decay is a mainstream prediction for GUT monopoles, with strong-interaction-scale cross sections, but no GUT monopole has been found and the energy ceiling is the atom's rest energy: 14,899 MeV for ¹⁶O, not 16 × 938.272 = 15,012 MeV. The "166×" could not be reproduced with any standard D–D convention (ratios are 64, 157, 257 or 1028). Label: **mainstream** (catalysis as theory), **refuted** (16 m_p as the yield; it exceeds E = mc² by 0.73%), **unverified** (any laboratory monopole).

### Cited Findings
- The talk: "that's 938 million electron volts per proton, or 16 times if it's an oxygen. Just one proton per nucleon gives you 166 times more energy than deuterium deuterium fusion" (transcript lines 2515–2520); "Typically, 1 kg will give you 30 33 days" (line 2522)
- C. G. Callan, "Dyon-fermion dynamics", *Phys. Rev. D* 26, 2058–2068 (1982), Crossref-checked; V. A. Kuzmin & V. A. Rubakov, "On the fate of superheavy magnetic monopoles in a neutron star", *Phys. Lett. B* 125, 372–376 (1983), Crossref-checked
- Catalysis cross sections are set by the strong-interaction scale (~m_p⁻²) rather than the monopole's geometric size, ~32 orders of magnitude larger, with current estimates 10⁻²⁷–10⁻²¹ cm² and IceCube flux limits of 10⁻¹⁷–10⁻¹⁸ cm⁻² s⁻¹ sr⁻¹ depending on the cross section — search summaries of [Brandenberger & Perivolaropoulos (1988)](https://olympias.lib.uoi.gr/jspui/bitstream/123456789/17319/1/perivolaropoulos-1988-Proton%20decay%20catalyzed%20by%20superconducting%20cosmic%20strings.pdf) and [IceCube, arXiv:1402.3460](https://arxiv.org/pdf/1402.3460) (not opened)
- The repository's `magnetics::frontier` already holds Dirac charge, Parker bounds and Cabrera's flux (`research/topics/magnetics-and-mhd/report.md`)

### Inferences
- *(computed here)*: 16 m_pc² = 15,012.35 MeV; M(¹⁶O atom)c² = 14,899.17 MeV; per kilogram the talk's figure gives 9.053×10¹⁶ J/kg = 1.00728 c². The excess comes from counting 8 neutrons as protons and ignoring the 127.6 MeV binding.
- Ratios against D–D *(computed here)*: m_p / mean D–D Q (3.651 MeV) = 257; per nucleon (938.272 vs 0.913 MeV) = 1028; one proton per ¹⁶O nucleus (58.6 MeV per nucleon) = 64.3; m_p vs D+D→⁴He per nucleon (23.85/4) = 157. None gives 166.
- Monopole catalysis would not "collapse" all 16 nucleons into energy in one event unless the monopole stays in the nucleus for 16 catalyses; the products (e⁺, π, γ, ν) carry 100s of MeV each and would be seen by any detector. "1 kg in 33 days" at the c² ceiling is 31.5 GW of gamma and pion power.

### Gaps
- Rubakov's original 1981 JETP Letters paper was not located on Crossref here; the cross-section numbers above are from search summaries.

## 15. Huang et al. (2024, *Scientific Reports*): water cavitation, ²²Ne without ²⁰Ne/²¹Ne, ¹⁷O, carbon

### Takeaway
The paper (Greenyer a co-author) reports COPx up to 4.26 and quadrupole-MS signals at m/z 22, 19, 33, 45 interpreted as ²²Ne, H₂¹⁷O, ¹⁶O¹⁷O and ¹²C¹⁶O¹⁷O, with no m/z 21 and no excess m/z 20; it reports peak ratios, not concentrations, gives no Q-values, and carries an editor's note that its conclusions are subject to criticisms under consideration. ¹⁷O+¹⁷O→¹²C+²²Ne is +6.41 MeV, not the audit's 17.13. Each m/z used has a common molecular isobar. Label: **contested** (published with an editor's expression of concern; no independent replication found).

### Cited Findings
- B.-J. Huang, …, R. W. Greenyer, "Water can trigger nuclear reaction to produce energy and isotope gases", *Sci. Rep.* 14, 214 (2 Jan 2024), doi:10.1038/s41598-023-50824-8 (Crossref-checked; preprint doi:10.21203/rs.3.rs-3466805/v1, 2023). Editor's note: "Readers are alerted that the conclusions presented in this article are subject to criticisms … that are being considered by the Editors." Reactors: VCS with ~3 kW compressor input, ~21 bar, ~2 mm gap, pulsed flow; new reactors up to 10 kW, 190 °C. Max COPx 4.26 (VCS), 2.55 (DHX); gas-sampled excess-heat runs COPx 1.10–1.61. Two quadrupole MS (Extorr XT200M, XT300M), resolution better than 0.5 amu; 14 samples from 8 reactors; results as intensity ratios, no concentrations. m/z 22: K22 up to 56.0; CO₂⁺⁺ said to be < 2% of m/z 44 and removed with a Ca(OH)₂ absorber; no m/z 21; m/z 20 below the air reference (calibrated with ⁴⁰Ar⁺⁺). ¹⁷O inferred from m/z 19 (H₂¹⁷O, with an HDO correction), 33 (¹⁶O¹⁷O) and 45 (¹²C¹⁶O¹⁷O, K45 2.1–75). Proposed reactions (eqs. 6–9) without Q-values; GM counter max 0.4 µSv/h, about twice background — [PMC10762109](https://pmc.ncbi.nlm.nih.gov/articles/PMC10762109/)
- A 2026 preprint by overlapping authors reports ¹⁷O production with excess heat across 11 reactors and 30 mol% ¹⁷O-enriched water from a 32 L system (not peer-reviewed) — search summary of [preprints.org 202607.0090](https://www.preprints.org/manuscript/202607.0090) (not opened)
- Natural neon: ²⁰Ne 90.48%, ²¹Ne 0.27%, ²²Ne 9.25% — [NIST](https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=ascii2&isotype=some)

### Inferences
- *(computed here)*: p + e⁻ + ¹⁶O → ¹⁷O + ν: +3.361 MeV (eq. 7 is charge- and lepton-balanced if "¹H" is a bare proton; the ν carries the energy, which is consistent with "no radiation" but makes it a weak-interaction rate); ¹⁶O(p,γ)¹⁷F: +0.600 MeV (the strong route, which gives β⁺-unstable ¹⁷F); ¹⁷O+¹⁷O→¹²C+²²Ne: +6.407 MeV; chain per ²²Ne: 13.13 MeV. The ¹⁷O+¹⁷O Gamow energy is 34.1 GeV; at the reactors' ≤ 190 °C (or 10⁴ K in a collapsing bubble) the rate is nil.
- Isobars at the m/z used: 22 — CO₂⁺⁺ (they absorb CO₂); 19 — H₃O⁺, ¹⁸OH⁺, F⁺; 33 — HO₂⁺, ¹⁶O¹⁶OH⁺, ³³S⁺; 45 — ¹³C¹⁶O₂⁺ (1.07% of CO₂), CHO₂⁺. Because K45 rises with K44 (CO₂ itself up to 71×), m/z 45 tracks ¹³CO₂ unless the 45/44 ratio departs from 0.012; the paper reports K-ratios to air, not 45/44.
- Both ways: an isotope-ratio mass spectrometer (or a high-resolution MS separating ²²Ne⁺ at 21.9914 u from ¹²C¹⁶O₂⁺⁺ at 21.9949 u, Δm/m ≈ 1.6×10⁻⁴) would settle m/z 22; ¹⁷O enrichment in the water by IRMS or laser spectroscopy is a direct test of the 2026 preprint's 30 mol%.

### Gaps
- No published critique text (the editors' concerns) or replication was found; the preprint was not opened.

## 16. Tungsten → calcium + xenon (Ohmasa gas, 10 s), Parkhomov's long run, LION, Langmuir 1912, "carbon the first alpha-conjugate product", "Ohmasa gas fissions heavy elements to safe elements"

### Takeaway
W → Ca + Xe is exothermic (46–85 MeV) but forbidden in practice by a ~150 MeV Coulomb barrier (spontaneous fission of W has never been seen); with stable products only ⁴⁸Ca/⁴⁶Ca/⁴⁴Ca pair with Xe, so W-derived calcium would be ⁴⁸Ca-rich, while natural calcium is 96.9% ⁴⁰Ca. EDX cannot see isotopes, and Ca is the most common contaminant (the talk itself notes calcium oxide). Label: **unverified** (no published isotope data found); the "fission to safe elements" proposal needs a barrier-crossing mechanism no measurement supports.

### Cited Findings
- The talk: "as it goes from tungsten to say calcium and xenon, just the calcium occupies twice … the atomic volume of the original nuclei of tungsten" (lines 1302–1306); "classic stellar synthesis products from the affected tungsten in 10 seconds" (line 1584); "Since a [Ohmasa] gas appears to fission heavy elements to safe elements, it should be suitable for treating dangerous products of the fission industry" (lines ~1620–1625)
- Parkhomov's M7 reactor ran 225 days (Oct 2017–May 2018; the talk says 221) with about 4,100 MJ excess from 1.2 g of nickel powder (2.1 MeV per Ni atom) and calcium in the SiC tube rising from about 1% to 23%, attributed to C + Si fusion via the weak interaction; presented at ICCF-22, Assisi, 2019 — [Remote View, "A Century of Cold Fusion in Russia"](https://remoteview.substack.com/p/a-century-of-cold-fusion-in-russia) (secondary)
- Natural W: ¹⁸⁰W 0.12%, ¹⁸²W 26.50%, ¹⁸³W 14.31%, ¹⁸⁴W 30.64%, ¹⁸⁶W 28.43%; natural Ca: ⁴⁰Ca 96.941%, ⁴²Ca 0.647%, ⁴³Ca 0.135%, ⁴⁴Ca 2.086%, ⁴⁶Ca 0.004%, ⁴⁸Ca 0.187%; natural Xe ¹³⁶Xe 8.86%, ¹³⁴Xe 10.44%, ¹³²Xe 26.91% — [NIST](https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=ascii2&isotype=some)

### Inferences
- Balanced splits into stable isotopes *(computed here)*: ¹⁸⁴W→⁴⁸Ca+¹³⁶Xe +84.95 MeV; ¹⁸²W→⁴⁸Ca+¹³⁴Xe +84.10; ¹⁸²W→⁴⁶Ca+¹³⁶Xe +81.32; ¹⁸⁰W→⁴⁸Ca+¹³²Xe +83.87. ¹⁸³W and ¹⁸⁶W have no split into stable Ca + stable Xe at all (max A = 48 + 136 = 184). ⁴⁰Ca needs ¹⁴⁰Xe–¹⁴⁶Xe partners (all β⁻-radioactive): ¹⁸⁰W→⁴⁰Ca+¹⁴⁰Xe +58.20, ¹⁸⁴W→⁴⁰Ca+¹⁴⁴Xe +46.01. The point-Coulomb barrier, Z₁Z₂ = 1080, at 1.2(48^⅓+136^⅓) = 10.5 fm to 11.5 fm, is ~135–148 MeV *(arithmetic, not scripted)*, so Q < V_B: tungsten is stable against this split, which is why it is not observed.
- Decisive test: ⁴⁸Ca/⁴⁰Ca in the "new" calcium (natural 0.0019) and xenon in the gas phase (a noble gas, invisible to EDX). Natural Ca ratios mean contamination; ⁴⁰Ca-dominated "W-derived" calcium would have to come with radioactive ¹⁴⁰⁻¹⁴⁶Xe → Cs → Ba → La chains, i.e. measurable activity.
- "Fission to safe elements" of actinide waste: fission of U/Pu fragments is what makes the waste; a process that split them into stable light products would have to beat the same barrier and conserve A and Z, and most two-body splits of A ≈ 90–140 into stable pairs do not exist for the same reason as ¹⁸³W/¹⁸⁶W above.
- Parkhomov's Ca 1% → 23% in a SiC tube: ²⁸Si+¹²C→⁴⁰Ca is +13.354 MeV *(computed here)*, so 1 g of Ca made this way would release 32 GJ *(arithmetic)*, against the 4.1 GJ of reported excess over the whole run; the point-Coulomb barrier (Z₁Z₂ = 84, R = 1.2(28^⅓+12^⅓) = 6.4 fm) is ~19 MeV against furnace kT ≈ 0.1 eV. Ceramic tubes and cements commonly carry Ca, which migrates to hot surfaces. Isotopes do not discriminate here (C+Si gives only ⁴⁰Ca, as natural Ca mostly is); the mass balance does.
- "Carbon is the first alpha-conjugate product" (⁸Be is unbound by 91.8 keV against 2α *(computed here)*, so ¹²C is the first stable α-conjugate): **mainstream** as nuclear structure; it says nothing about the carbon's origin. Carbon is the dominant EDX contaminant.

### Gaps
- No MFMP/Ohmasa tungsten paper with isotope data, no LION reactor report, and no Langmuir 1912 source for "7,000 times the volume" of hydrogen were found (hydrogen dissolved in and adsorbed on tungsten from its manufacture is the standard explanation; Wendt & Irion outgassed their wires for this reason, section 17).

## 17. Wendt & Irion (1922): exploding tungsten wires and helium, and the refutation

### Takeaway
Wendt & Irion saw the helium yellow line and 0.18–3.62 cc of unabsorbed gas from exploding 0.5–0.8 mg tungsten wires with a ~0.1 µF condenser at 22–37 kV, and called the work preliminary. Splitting 0.713 mg of ¹⁸⁴W into helium costs 58.8 MJ (it is also charge-unbalanced without 18 β⁻ per atom); each shot stored ~50 J. Allison & Harkins (1924) reported "the absence of helium". Label: **refuted**.

### Cited Findings
- G. L. Wendt & C. E. Irion, "Experimental Attempts to Decompose Tungsten at High Temperatures", *J. Am. Chem. Soc.* 44, 1887–1894 (1922), doi:10.1021/ja01430a006 (Crossref-checked): wires 0.035 mm diameter, ~4 cm, 0.5–0.7 mg; 300 cc Pyrex bulbs; condenser of 100 glass plates, capacity "about 0.1 microfarad", charged to ~30 kV (holding 45 kV); wires outgassed 15 h to remove hydrogen and helium "which may arise from an impurity of thorium in the metal"; the helium yellow line "always brilliant"; table of 21 explosions, 0.18–3.62 cc of gas; "If the entire mean weight of 0.713 mg. of tungsten had been decomposed into helium this should have occupied almost exactly 4 cc"; in CO₂, 0.713 mg gives 1.01 cc unabsorbed in KOH; "nothing is proved beyond the importance of the problem" — [Zenodo scan of the article](https://zenodo.org/record/1558743)
- S. K. Allison & W. D. Harkins, "The absence of helium from the gases left after the passage of electrical discharges …", *J. Am. Chem. Soc.* 46, 814–824 (1924), doi:10.1021/ja01669a004 (Crossref title; full text not opened)

### Inferences
- *(computed here)*: their 4 cc is right (3.999 cc). ¹⁸⁴W (Z 74) → 46 ⁴He (Z 92) needs 18 β⁻: Q = −157.25 MeV per atom; 0.713 mg → 58.8 MJ against ½CV² = 49.6 J at 31.5 kV, a factor of 1.2×10⁶. The helium line is explained by the thorium impurity they themselves name, by air leaks, or by helium from the glass.

## 18. Molar magnetic susceptibility: "nothing comes close to oxygen"

### Takeaway
O₂ is the only strongly paramagnetic elemental gas, at +3449×10⁻⁶ cm³/mol (CRC, 293 K) = 4.33×10⁻⁸ m³/mol, within 1% of Moldover et al.'s 42.92×10⁻⁹ m³/mol. Among elements at room temperature, Tb, Ho, Dy, Er and Tm are 8–51 times larger per mole (Nd, Pr also exceed O₂), and Fe, Co, Ni and Gd are ferromagnetic at room temperature (Fe, Co, Ni far above it), so "nothing comes close" holds for gases only. Label: **mainstream** for "O₂ is the most paramagnetic elemental gas"; **refuted** for "among all elements".

### Cited Findings
- O₂ gas, zero-pressure limit, 293.15 K: χ_m = (42.92 ± 0.06)×10⁻⁹ m³ mol⁻¹, second magnetic virial coefficient (−1.8 ± 0.5) cm³ mol⁻¹ — [Moldover, May & Schmidt, *Phys. Rev. A* 78 (2008), NIST page](https://www.nist.gov/node/643416)
- Molar susceptibilities, SI m³/mol (Wolfram ElementData): Tb 2.161385×10⁻⁶, Ho 9.05467×10⁻⁷, Dy 8.85625×10⁻⁷, Er 6.30566×10⁻⁷, Tm 3.36179×10⁻⁷, Nd 6.9235×10⁻⁸, Pr 5.9604×10⁻⁸, O 4.27184×10⁻⁸, Eu 4.1942×10⁻⁸, Ce 3.0826×10⁻⁸; Gd, Fe, Co, Ni listed N/A (ferromagnetic); N₂ −1.5×10⁻¹⁰, Xe −5.65×10⁻¹⁰, W +8.44×10⁻¹⁰, C −7.45×10⁻¹¹ — [periodictable.com, Molar Magnetic Susceptibility](https://periodictable.com/Properties/A/MolarMagneticSusceptibility.al.html)
- The CRC table "Magnetic Susceptibility of the Elements and Inorganic Compounds" lists χ_m in 10⁻⁶ cm³/mol (CGS) — [Taylor & Francis chapter listing](https://www.taylorfrancis.com/books/9781439855126/chapters/10.1201%2Fb17379-51) (not opened; the 3449 value is the user's audit's and the coordinator's)

### Inferences
- *(computed here)*: CGS→SI is ×4π×10⁻⁶ (the audit's conversion is correct): 3449×10⁻⁶ → 4.334×10⁻⁸ m³/mol, 1.0% above Moldover. Per mole Tb/O₂ = 50.6, Dy/O₂ = 20.7; per atom (O₂ has two O atoms) double those. By Curie's law the 273.15 K ("STP") value is 7% higher, ≈ 3700×10⁻⁶; the audit labels the 293 K CRC value "at STP".
- Gd's Curie point is near room temperature, Dy's and the other heavy rare earths' far below it; the speaker's "only ferromagnetic below room temperature" fits them but not Fe, Co, Ni.

### Gaps
- The CRC table itself could not be opened (the FNAL and e-magnetica copies returned 403); 3449×10⁻⁶ is consistent with Moldover to 1%.

## 19. Detection and contamination: what distinguishes transmutation

### Takeaway
Elemental appearance (EDX, SEM) cannot separate transmutation from contamination; isotope ratios and mass balance can. Each pathway in the talk predicts a specific, testable isotopic signature; only Singh et al. (natural Fe) and Adamenko (non-natural S, Ti, Fe) report isotope data, and Huang et al. report molecular-ion ratios, not isotope ratios.

### Cited Findings
- Natural abundances for the comparison: Fe 5.845/91.754/2.119/0.282% (A 54/56/57/58); Ca 96.941% ⁴⁰Ca, 0.187% ⁴⁸Ca; Ne 90.48/0.27/9.25% (20/21/22); O 99.757/0.038/0.205%; Cu 69.15/30.85% (63/65); Pb 1.4/24.1/22.1/52.4% (204/206/207/208) — [NIST](https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=ascii2&isotype=some)
- Isotope data reported: natural Fe (Singh et al. 1994, [ANS](https://wx1.ans.org/pubs/journals/fst/a_30331)); shifted S, Ti, Fe and ⁵⁷Fe > ⁵⁶Fe (Adamenko, [Dolan review](https://www.newenergytimes.com/v2/books/Reviews/AdamenkoByDolan.pdf), [Kowalski](https://msuweb.montclair.edu/~kowalskil/cf/217kiev.html)); none for Sundaresan & Bockris (abstract), Srinivasan, Matsumoto, Kladov, Ohsawa, the W→Ca+Xe work

### Inferences
- Signatures per pathway: α-conjugate fusion (C+O, Si+Si, 4 ¹⁷O) → pure ⁵⁶Fe, ²⁸Si; W→Ca+Xe → ⁴⁸Ca-rich Ca and Xe gas; Cu→Pb (any route) → Pb with no ²⁰⁴Pb unless made specifically (²⁰⁴Pb is primordial-only); ¹⁷O+¹⁷O→²²Ne → ²²Ne with no ²⁰Ne (as Huang claim) — which a resolving MS must separate from CO₂⁺⁺. Natural ratios in the product favour contamination; non-natural ratios need SIMS/ICP-MS with blanks, as Adamenko's lab says it did.
- Every exothermic route above emits MeV radiation per event; a mass of new element m (g) of mass number A and Q per event releases m/A × N_A × Q. For µg-to-tonne claims this is kJ to TJ. Absence of gammas, neutrons and activation is a strong constraint on any nuclear reading.

## 20. The user's audit: its numbers checked

### Takeaway
Of the audit's five nuclear numbers, 12C+16O→28Si (16.76 MeV) is right, 2 ²⁸Si→56Fe (10.92 MeV) is the value for ⁵⁶Ni, 4 ¹⁷O→⁵⁶Fe+¹²C (49.47) and 2 ¹⁷O→¹²C+²²Ne (17.13) are wrong (57.37 and 6.41), and 16 × 938.272 MeV per ¹⁶O exceeds the atom's rest energy. Its O₂ susceptibility and CGS→SI factor are right (the "at STP" label is not). The specific energies it derives inherit the Q errors, and the monopole one exceeds c².

### Cited Findings
- Audit values: 4 ¹⁷O → ⁵⁶Fe + ¹²C Q ≈ +49.47 MeV; ¹²C+¹⁶O→²⁸Si +16.76 MeV; 2 ²⁸Si → ⁵⁶Fe + 2e⁺ + 2ν +10.92 MeV; 2 ¹⁷O → ¹²C + ²²Ne +17.13 MeV; monopole decay of ¹⁶O 16 × 938.272 MeV; O₂ χ_m ≈ +3449×10⁻⁶ cm³/mol, "×4π×10⁻⁶ for SI"; `specific_energy_joules_per_kg_o2` = Q/(n_O · A · u) — `research/external/2026-10-26-greenyer-audit-and-crate.md` lines 46, 127–143, 297–299, 594–626

### Inferences
| Audit item | Audit value | Correct (AME2020) | Likely source of the error |
|---|---|---|---|
| 4 ¹⁷O → ¹²C + ⁵⁶Fe | +49.47 MeV | **+57.372 MeV** | Not reproduced: ¹⁶O in place of ¹⁷O does not balance A, ⁵⁶Ni in place of ⁵⁶Fe does not balance Z; source unidentified |
| ¹²C + ¹⁶O → ²⁸Si | +16.76 MeV | +16.756 MeV | correct |
| 2 ²⁸Si → ⁵⁶Fe + 2e⁺ + 2ν | +10.92 MeV | **+15.578 MeV** kinetic (+17.622 MeV incl. annihilation); 10.922 MeV is 2 ²⁸Si → **⁵⁶Ni** | Q for the Ni compound nucleus quoted with the Fe label; Z: 2 × 14 = 28 = Ni |
| 2(C+O) → ⁵⁶Fe, "44.44 MeV per ⁵⁶Fe" | 44.44 | 44.434 to ⁵⁶Ni; **51.133** to ⁵⁶Fe (atomic) | same Ni/Fe slip |
| 2 ¹⁷O → ¹²C + ²²Ne | +17.13 MeV | **+6.407 MeV** (chain from ¹⁶O + 2p + 2e⁻: 13.13) | Possibly 13.13 (the full chain) mistyped as 17.13; unconfirmed |
| Monopole decay of ¹⁶O | 16 × 938.272 = 15,012.35 MeV | ceiling **14,899.17 MeV** = M(¹⁶O)c² | counts 8 neutrons as protons, drops 127.6 MeV binding; gives 1.0073 c² per kg |
| O₂ χ_m | +3449×10⁻⁶ cm³/mol "at STP" | 3449×10⁻⁶ at 293 K (≡ 4.334×10⁻⁸ m³/mol; Moldover 4.292×10⁻⁸, −1%); ≈ 3700×10⁻⁶ at 273 K | label only |
| J/kg O (audit formula) | Ohsawa 7.02×10¹³; Srinivasan 1.34×10¹⁴; Huang 4.86×10¹³; monopole 9.053×10¹⁶ | **8.14×10¹³; 1.542×10¹⁴ (to Fe); 1.818×10¹³; ≤ 8.988×10¹⁶ (= c²)** | inherits the Q errors; the monopole case breaks E = mc² |

All corrected values are asserted in `transmutation_claims_checks.py` (section 2, 8, 9).

## 21. Calculators for the crates

### Takeaway
One new module, `plasmoid::nuclear`, covers every accounting above: a small embedded AME2020 table, a reaction-balance and Q-value function that refuses unbalanced reactions, and helpers for power, activity and specific energy; plus one new function in `plasmoid::fusion` (thermal tunnelling) and one unit conversion in `magnetics::magnets`. Masses needed (42): n, ¹H, ²H, ³H, ³He, ⁴He, ¹²C, ¹⁶O, ¹⁷O, ¹⁷F, ²⁰Ne, ²²Ne, ²³Na, ²⁴Mg, ²⁷Al, ²⁷Si, ²⁸Si, ³⁰Si, ³⁰P, ³¹P, ³¹S, ³²S, ³⁹K, ⁴⁰Ca, ⁴⁴Ca, ⁴⁶Ca, ⁴⁸Ca, ⁵⁶Fe, ⁵⁶Co, ⁵⁶Ni, ⁶³Cu, ⁶⁵Cu, ¹⁸⁰W, ¹⁸²W, ¹⁸⁴W, ¹³²Xe, ¹³⁴Xe, ¹³⁶Xe, ¹⁴⁰Xe, ¹⁴⁴Xe, ¹⁹⁷Au, ²⁰⁸Pb (values and atomic masses in the check script's `AME` table; the two columns cross-check to 0.0012 keV).

### Cited Findings
- AME2020 — [IAEA AMDC](https://www-nds.iaea.org/amdc/ame2020/mass_1.mas20.txt); NIST isotopic compositions — [NIST](https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=ascii2&isotype=some); oxygen-burning Q-values — [Wikipedia (Woosley et al. 2002; Clayton 1983)](https://en.wikipedia.org/wiki/Oxygen-burning_process); existing `plasmoid::fusion` (`gamow_energy`, `gamow_peak`, `tunnelling`, `coulomb_barrier`, `dd_per_watt`, `charged_sphere`)

### Inferences
| Function (crate::module) | Formula | Primary source | Test value | Status |
|---|---|---|---|---|
| `Nuclide { z, a, mass_excess_kev }`, `AME2020: &[Nuclide]`, `nuclide(z, a)` (`plasmoid::nuclear`, new) | table lookup; `(M − A)u = Δ` cross-check | AME2020 | 42 entries; `(M−A)u − Δ` < 0.01 keV | mainstream |
| `q_value(lhs, rhs, positrons, weak) -> Result<f64, Imbalance>` (`plasmoid::nuclear`) | ΣΔ_in − ΣΔ_out − 2mₑc²·positrons; requires ΔA = 0, ΔZ = positrons + weak | AME2020 | D+D→³He+n 3.2689 (3.268); D+D→T+p 4.0327 (4.033); ¹⁶O+¹⁶O→²⁸Si+α 9.594 (9.593), →³¹P+p 7.678 (7.676), →³²S 16.542 (16.539); ¹²C+¹⁶O→²⁸Si 16.756 (16.76) | mainstream |
| same, slide/talk reactions | as above | AME2020 | 2 ²⁸Si→⁵⁶Ni 10.922; ⁵⁶Ni→⁵⁶Fe 6.700 (weak 2); 4 ¹⁷O→¹²C+⁵⁶Fe 57.372; ¹⁷O+¹⁷O→¹²C+²²Ne 6.407; ²³Na+¹⁶O→³⁹K 19.540; ¹⁸⁴W→⁴⁸Ca+¹³⁶Xe 84.949; 4 ⁶³Cu→²⁰⁸Pb+11α (weak 12) −267.2; ¹⁸⁴W→46α (weak −18) −157.25 | mainstream data; claims labelled per section |
| `stable_splits(parent, z1)` (`plasmoid::nuclear`) | pairs (A₁, A − A₁) with both stable | AME2020 + NIST (stability = nonzero abundance) | ¹⁸⁴W, Z₁ = 20 → only (48, 136); ¹⁸³W, ¹⁸⁶W → none | mainstream |
| `production_power(kg_per_s, a, q_mev)` (`plasmoid::nuclear`) | ṁ/(A u)·Q | — | 4 t/day ⁵⁶Fe via Si+Si → 8.72×10¹¹ W; via 2(C+O) → 4.08×10¹² W | mainstream arithmetic |
| `max_atoms(energy_j, q_mev)` (`plasmoid::nuclear`) | E/|Q| | — | 300 J at −267.2 MeV → 7.01×10¹² Pb atoms, 2.42 ng | mainstream arithmetic |
| `steady_activity(rate_per_s)`, `activity(n, half_life)` (`plasmoid::nuclear`) | A = R; A = N ln2/t½ | half-lives: Wikipedia silicon burning (6 d, 77 d) | 4 t/day → 1.35×10¹³ Ci; 1 µg ⁵⁶Ni → 0.389 Ci | mainstream |
| `specific_energy(q_mev, fuel_mass_u)` (`plasmoid::nuclear`), with `assert ≤ c²` | Q/(m u) | AME2020 | 4 ¹⁷O→Fe+C 8.14×10¹³ J/kg; whole ¹⁶O = c² = 8.988×10¹⁶; 16 m_p/16 u = 1.0073 c² (rejected) | mainstream |
| `log10_thermal_tunnelling(e_g, kt)` (`plasmoid::fusion`, add) | −3E₀/(kT ln 10), E₀ from `gamow_peak` | fusion.rs; D–T E₀ = 30.9 keV at 10 keV | ¹²C+¹⁶O: −2707.5 (5000 K), −2149.0 (10⁴ K); ¹⁶O+¹⁶O −2740.4 (10⁴ K); ²⁸Si+²⁸Si −8773 (5000 K); ²³Na+¹⁶O −11397 (310 K) | mainstream |
| `chi_molar_cgs_to_si(chi)` (`magnetics::magnets`, add) | 4π×10⁻⁶ χ | Moldover et al. 2008 | 3449×10⁻⁶ → 4.334×10⁻⁸ m³/mol (42.92×10⁻⁹, 1%) | mainstream |
| `curie_scale(chi, t0, t)` (`magnetics::magnets`, add) | χ·t₀/t | — | 3449×10⁻⁶ at 293 K → 3700×10⁻⁶ at 273.15 K | mainstream (paramagnets) |
| example `examples/transmutation_audit.rs` (`physica`) | calls the above for each talk claim | this note | sections 5, 8, 11–17, 20 numbers | example, not library (per "rigs are examples") |
| reuse `charged_sphere(q, r)` (`plasmoid::fusion`, existing) | kQ/R, (3/5)Q²/(4πε₀R) | US 5,018,180 | 10¹¹ e, 1 µm → 288 MV | existing |

### Gaps
- Not checked: whether `metrology` exposes `M_E` in keV and `AMU` in keV with the CODATA 2018 values the script uses (fusion.rs implies yes); the half-lives of ⁵⁶Ni/⁵⁶Co are from a secondary page (6 d, 77 d), and an ENSDF value should replace them in the crate.

## Built into `plasmoid::nuclear` (2026-10-27, #22)

- **AME2020 table:** every one of this note's 43 hand-entered masses matches the evaluation's `mass_1.mas20` exactly. The crate's table adds ⁶Li, ⁷Li, ¹¹B, ¹¹C and ¹⁴N.
- **Iron mass:** the script took 4 t of iron as 55.845 u per atom, natural iron's atomic weight. An α-conjugate route makes pure ⁵⁶Fe at 55.935 u, so the figures are now 4.31×10²⁸ atoms/day, 8.722×10¹¹ W (Si+Si) and 4.083×10¹² W (C+O). With 55.845 u the earlier 8.74×10¹¹ W was right, and the script's 8.72 passed only on a loose 2×10⁻³ tolerance.
- **p + ¹¹B → ¹¹C + n** is −2.7640 MeV from AME2020 (the research sweep quoted −2.7648).
- **`aneutronic::Q_PB11`** was 8.6821 MeV, from NIST's AME2012-based masses. AME2020 moves ¹¹B by 0.2 µu, giving 8.6819 MeV. The crate now uses 8.6819, and a test pins all seven aneutronic Q constants to the table.
- **Natural abundances:** NIST's full table, 288 isotopes in 84 elements, each element summing to 1. `stable_splits` takes "occurs in nature" from it, so ⁴⁸Ca and ¹³⁶Xe count.
