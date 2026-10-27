# Ball lightning, plasmoids, Shoulders' charge clusters, ectons and iron microspheres

Scope: the physics behind the ball-lightning, plasmoid, EVO and microsphere claims in Greenyer's TeslaTech talk (15 Aug 2026; transcript lines cited as "talk L…"). Every recomputed number is asserted in `ball_lightning_and_clusters_checks.py` (39 checks, standard library, about 0.01 s, `python3 -I`, all pass). Downloads are in `dl_ball/` (read only). The prior topic `research/topics/plasmoids-and-fusion/` already covers the ball-lightning virial cap (`plasmoid::relaxation::virial_energy_cap`), Finkelstein–Rubinstein, Paiva and the silicon chemistry. Those findings are cross-referenced here, not repeated. Labels: **[mainstream] [contested] [refuted] [unverified]**.

## Q1. Ball-lightning observations, laboratory analogues and energies

### Takeaway
The only instrumented natural ball lightning (Qinghai, 2012) lasted 1.64 s. It had an apparent diameter of about 5 m and showed neutral Si, Fe and Ca lines from soil throughout its life. Its 99.4 Hz flicker matches the second harmonic of a 50 Hz power line 20 m away. The data fit the Abrahamson–Dinniss model, where the ball is a network of oxidising Si nanoparticles. Abrahamson–Dinniss's own numbers recompute exactly: 100 g m⁻³ of Si stores 3.24 MJ m⁻³ (paper: 3.3), and a 300 mm ball releases 1,117 W (paper: about 1,100). Laboratory analogues (Si arcs, microwave fireballs) last milliseconds to seconds and are all chemically or externally powered. **[mainstream; one natural event]**

### Cited Findings
- Cen, Yuan & Xue, PRL 112, 035001 (2014), DOI 10.1103/PhysRevLett.112.035001 (Crossref-checked). The event occurred at 21:54:59 on 23 July 2012, during a thunderstorm on the Qinghai Plateau, at 0.902 km (±0.0077 km). The "entire luminous duration of the BL is 1.64 s". The ball formed at the bottom of a cloud-to-ground channel. The camera's focal length was 3.3 mm and its pixel pitch 4 µm, giving "about 1.1 m" per pixel. The centre moved nine pixels between 20 ms and 1160 ms, an "average 2D speed … about 8.6 m/s". The apparent diameter was about 5 m. The authors call this "the luminance range", not the real size. The colour went purple, then orange at 80 ms, white from 160 to 1100 ms, and red after 1120 ms. "Si I, Fe I, and Ca I exist during most lifetime of the BL". The intensity oscillated at "99.4 Hz", which the authors attribute to 35 kV lines 20 m away. Sound ranging (a 2.64 s delay at 340.51–341.15 m/s) gave 0.898–0.9 km — [PRL PDF](https://fisherp.mit.edu/wp-content/uploads/2020/04/PhysRevLett.112.035001.pdf). Recomputed: 1.093 m per pixel, 8.63 m/s, and 898.9–900.6 m *(check)*.
- Cen et al.'s introduction gives the textbook ranges: diameter 1–100 cm, speed of a few m/s, lifetime 1–10 s. It summarises the analogues as follows. Paiva et al. ran 20–25 V, 100–140 A, 60 Hz arcs in pure silicon and got 1–4 cm balls lasting 2–5 s. Stephan & Massey made 0.1–1 mm Si balls at about 3140 K. Dikhtyar & Jerby made fireballs of about 3 cm that lasted 30–40 ms after the microwave power was turned off — [Cen 2014](https://fisherp.mit.edu/wp-content/uploads/2020/04/PhysRevLett.112.035001.pdf). Paiva PRL 98, 048501 (2007) and Dikhtyar & Jerby PRL 96, 045002 (2006) were Crossref-checked; only their metadata was opened.
- Abrahamson & Dinniss, Nature 403, 519 (2000), DOI 10.1038/35000525 (Crossref-checked). Their reference case is a ball 300 mm across, lasting about 10 s, with the luminosity of a 100 W lamp. They found soil-strike deposits of chains of 5–70 nm spheres. The modelled ball has 100 g m⁻³ of Si loading in 25 nm spheres. Lifetimes before melting are about 30, 10 and 3 s at 1,200, 1,300 and 1,400 K. The heat release is 79.0 kW m⁻³ at 1,473 K, with an activation energy of 189 kJ mol⁻¹. "The total chemical energy density is 3.3 MJ m⁻³", against 15 MJ m⁻³ quoted for average ball lightning. Total heat release is "here 1,100 W" — [Nature PDF](https://scienceinschool.web.cern.ch/wp-content/uploads/2009/08/issue12_nature_abrahamson2000.pdf). Recomputed: 3.24 MJ m⁻³, 1,117 W, and an Arrhenius ratio k(1400)/k(1300) = 3.49 against the lifetime ratio of 3.33. The 1,200/1,300 K ratio is looser (4.29 against 3) *(check)*.
- Ito, Tamura, Cappelli & Hamaguchi, PRE 80, 067401 (2009), DOI 10.1103/PhysRevE.80.067401. Silicon-arc balls had temperatures of 2800–3400 K. Most lasted about 0.1 s, and more buoyant ones lasted several seconds. Their model is a 100 µm molten Si core inside a 1.5 mm, 90%-porous nanowire network — [OUKA PDF](https://ir.library.osaka-u.ac.jp/repo/ouka/all/78475/PhysRevE_80_06_067401.pdf).
- Witness statistics as summarised by Donoso & Rañada (2003): typical diameter 10–50 cm; lifetime distribution peaking at 2–5 s, with a mean of about 10 s or more; visible output 10–150 W; released energy "estimated between 10 kJ and more than 1 MJ". Submarine battery short circuits produced balls at about 150 kA, with an estimated 200–400 kJ — [arXiv:physics/0312122](https://arxiv.org/pdf/physics/0312122).
- Ohtsuki & Ofuruton, "Plasma fireballs formed by microwave interference in air", Nature 350, 139 (1991), DOI 10.1038/350139a0. A correction appeared in Nature 353, 868 (Crossref-checked). The full text was not opened. A 2.45 GHz, 5 kW source is reported second-hand in the plasmoid notes.
- Kapitsa's hypothesis is cited by Cen et al. as Sov. Phys. JETP 30, 973 (1970). Delayen & Mammosser (1999) recount that Kapitza (1958) made stable filament rings, lasting "a few seconds", in He-filled spherical resonators at about 1.5 GHz — [PAC 1999 MOP113](https://jacow.org/p99/papers/MOP113.pdf).
- Cross-reference: the virial cap 3p₀V gives 1.273 kJ for a 20 cm ball and 19.9 MJ for a 5 m ball. Silicon oxidation fits the energy budget easily — `research/topics/plasmoids-and-fusion/report.md` §"Ball lightning cannot be a self-confined plasmoid" *(check)*.

### Inferences
- The Qinghai ball's 5 m is a luminance size. Cen et al. themselves note that silicon balls with a 1–4 cm apparent size are lit by a core of about 1 mm. So the Qinghai event does not need a 5 m self-confined plasmoid.
- Fe lines in the only natural spectrum come from soil, which is about 5% Fe by mass (see Q7). Iron in lightning-related debris is therefore expected without any exotic source.
- The Abrahamson–Dinniss store is 10.7 times the virial cap per unit volume. This is allowed, because the cap limits field and thermal energy, not chemical energy *(check)*.

### Gaps
- Stenhoff's book and the Bychkov and Shmatov papers were not opened. No energy-density table is quoted from them.
- The full texts of Paiva 2007, Ohtsuki 1991 and Dikhtyar & Jerby 2006 were not opened. Their numbers here come via Cen 2014.

## Q2. The virial theorem and "metastable structures in virial equilibrium"

### Takeaway
Shafranov's virial theorem says that an ideal-MHD equilibrium of finite extent must be held by external stresses: externally supplied currents, wall forces, gravity or ambient pressure. In air the only available external stress is atmospheric pressure, which caps field-plus-thermal energy at 3p₀V. Faddeev, Freyhult, Niemi and Rajan extended this result to a two-fluid field theory. Their extension keeps the ban on isolated self-confined plasma "in an otherwise empty space". Stable knotted solitons are allowed only inside a uniform background plasma, and only with a nontrivial electric field. Air at room temperature is not such a background. The talk's phrase "metastable … in virial equilibrium" is coherent only if it means held by air pressure or material strength, and then the 3p₀V cap applies. **[mainstream theorem; the self-confined reading is refuted; the knotted-soliton loophole is contested]**

### Cited Findings
- Faddeev et al., arXiv:physics/0009061 (2000). They write that Shafranov's theorem "implies that nontrivial magnetohydrodynamical equilibrium configurations must be supported by externally supplied currents". Scaling the energy gives ∫d³x (3p + B²/2) = 0, a sum of positive terms, "and we conclude that … non-trivial localized equilibrium configurations do not exist in ideal magnetohydrodynamics". In the two-component field theory they find that "in an otherwise empty space an initially localized plasma configuration can not be confined solely by its internal electromagnetic interactions". With a background of constant asymptotic density, however, "a virial argument can not exclude the existence of stable finite energy solitons". These solitons "are accompanied by a nontrivial electric field" and are suggested for coronal loops, "maybe even ball lightning" — [arXiv PDF](https://arxiv.org/pdf/physics/0009061).
- Donoso & Rañada (2003) write out the magnetic virial theorem with the moment of inertia I, the kinetic, internal and magnetic energies, and a surface pressure term −3p_ext V. They state that a BL's magnetic energy "must be of tens of kJ, even reaching several MJ". Without balance it would explode "in a small fraction of a second, so that the virial theorem would exclude indeed any electromagnetic model". Their way out is force-free, linked "magnetic knots", in which the plasma sits in streamers filling about one ppm of the volume — [arXiv:physics/0312122](https://arxiv.org/pdf/physics/0312122). Primary papers: Rañada & Trueba, Nature 383, 32 (1996), DOI 10.1038/383032a0; Rañada, Soler & Trueba, PRE 62, 7181 (2000), DOI 10.1103/PhysRevE.62.7181 (both Crossref-checked, not opened).
- Faddeev & Niemi, "Stable knot-like structures in classical field theory", Nature 387, 58 (1997), DOI 10.1038/387058a0 (Crossref-checked, not opened).
- Finkelstein & Rubinstein, Phys. Rev. 135, A390 (1964): air pressure bounds stored energy below some reported values — cross-referenced from `plasmoids-and-fusion/notes/plasmoids_frc_reconnection.md`.
- A minority view (Soshnikov, arXiv:1007.4377) holds that "the virial theorem limiting BL power does not apply" to a moving, field-fed, open ball. It gives a maximum energy density of 10⁸–10⁹ J m⁻³ — [arXiv PDF](https://arxiv.org/pdf/1007.4377). **[contested minority]**

### Inferences
- A force-free field makes the volume magnetic term equal the boundary Maxwell-stress term, so the stress still has to be taken up at the boundary. The knot models move the problem onto the boundary; they do not remove it. Inside a finite volume, the theorem requires a boundary stress of at least about B²/2μ₀. In air the most that pressure can confine is 0.505 T, from √(2μ₀p₀) in the plasmoid report.
- The talk's resonance or "feeding off oxygen" picture (talk L2365–2400) is an open, chemically fed system. It is consistent with the Abrahamson–Dinniss class of models, not with a self-confined "virial equilibrium".

### Gaps
- Shafranov's original statement (Reviews of Plasma Physics, vol. 2, 1966) was not opened. The statement here comes via Faddeev et al.

## Q3. Shoulders' EVs and charge clusters against Mesyats's ectons

### Takeaway
Shoulders' patent gives EVs of about 1 µm carrying "typically on the order of 10¹¹ electron charges", with at most one ion per 100,000 electrons, moving at about 0.1c. The mainstream explosive-emission literature (Mesyats) describes the same order of charge packet. Mesyats's ectons emit "individual portions of electrons (10¹¹–10¹² particles)" in 10⁻⁹–10⁻⁸ s from cathode micro-explosions, and they leave micron craters. A 1.6 A ecton running for 10 ns carries exactly 1.0×10¹¹ electrons *(check)*. The packet size is therefore mainstream. The claimed compactness and near-zero neutralisation are not physically possible. A bare 10¹¹-electron, 1 µm sphere would sit at 288 MV with a surface field of 5.8×10¹⁴ V m⁻¹ and hold 2.77 J. It would Coulomb-explode in about 7×10⁻¹⁷ s, and binding it magnetically at 0.1c needs ions on 99% of electrons (Budker), not 0.001%. **[packet charge: mainstream (ecton); stable, compact, nearly unneutralised EV: refuted by standard electrostatics; "makes atoms disappear": unverified]**

### Cited Findings
- US 5,018,180, "Energy conversion using high charge density" (priority 3 May 1988; filed 3 May 1989; assignee Jupiter Toy Co.; inventor K. R. Shoulders). Size: "on the order of one micrometer in lateral dimension, but can be larger or as small as 0.1 micrometer". Charge: "typically on the order of 10¹¹ electron charges"; "the minimum charge observed for a one micrometer EV is 10⁸". A ring of ten beads "may include 10¹² electron charges". Ions: "upper limit of one ion per 100,000 electron charges". Speed: "approximately one-tenth the speed of light". On dielectric guides, a channel of about 20 µm "will accommodate EV passage without boring", while the EV "may effectively cut into the guide material". The witness plate is a foil that "will sustain visible damage" — [Google Patents](https://patents.google.com/patent/US5018180A/en) (only the first 100 kB of the page was read).
- Mesyats, "Ecton or electron avalanche from metal", Phys. Usp. 38, 567 (1995), DOI 10.1070/PU1995v038n06ABEH000089 (Crossref-checked). An ecton is driven by energy density "(10⁴ J g⁻¹)". Electrons are emitted for 10⁻⁹–10⁻⁸ s at "a total energy of 10⁻⁸ J". Cathode plasma density exceeds 10²⁰ cm⁻³. Craters have r_cr ≈ 10⁻⁴ cm (Fig. 1 shows a Mo crater about 1 µm across, made in 10 ns). The specific action is h̄ = ∫j²dt, with Table 1 giving Cu = 2×10⁹ A² s cm⁻⁴. The current density is j = (h̄/t_d)^½, "about 10⁹ A cm⁻² for t_d ~ 10⁻⁹ s". The paper gives t_e (Eq. 2.12), the electrons per ecton n_e = i³(256π²a²e h̄ sin⁴(θ/4))⁻¹ (Eq. 2.14) and the specific erosion γ_m = (2/3)ρ(a/h̄)^½ (Eq. 2.15). The model gives γ_m ≈ 0.6×10⁻⁴ g C⁻¹ against measured (0.4–1.15)×10⁻⁴ g C⁻¹, with the threshold i_th = 1.6 A and t_c ≈ 22 ns for silver. The conclusion: "individual portions of electrons (10¹¹–10¹² particles) … called ectons". Plyutto's collective ion acceleration may be due to "electron clusters in an ecton" — [UFN PDF](https://ufn.ru/ufn95/ufn95_6/ufn956a.pdf).
- Recomputed against Mesyats: j = 10⁹ A cm⁻² exactly. A packet of 10¹¹ e (1.6×10⁻⁸ C) at γ_m = 0.6×10⁻⁴ g C⁻¹ erodes 0.96 pg, which needs 9.6×10⁻⁹ J at 10⁴ J g⁻¹, matching Mesyats's 10⁻⁸ J per ecton *(check)*.
- Budker condition for a self-pinched beam: the electron-neutralisation fraction must exceed 1/γ² — Lampe review, [arXiv:physics/0409157](https://arxiv.org/pdf/physics/0409157) (search excerpt only, not opened).

### Inferences
- Derived for a uniform sphere of 10¹¹ e with d = 1 µm *(check)*:
  - Coulomb self-energy 3kQ²/5R = 2.77 J; surface potential 288 MV; surface field 5.8×10¹⁴ V m⁻¹.
  - Electron density 1.9×10²⁹ m⁻³, which is 2.25 times the conduction-electron density of copper.
  - Explosion time √3/ω_p = 7.0×10⁻¹⁷ s.
  - Kinetic energy at 0.1c: 2.57 keV per electron, 41 µJ in total.
  - Budker needs f_e > 0.99, about 10⁵ times the patent's upper ion limit.
  - Holding the cloud with a static field (the Brillouin limit) would take 2.0×10⁵ T.
- The self-energy is 2.9×10⁸ times the energy of one ecton. A 1 µm "EV" carrying the patent's charge therefore cannot be a bare electron cluster. Either it is nearly neutral (a dense cathode-plasma blob, as in Mesyats), or it is much larger than 1 µm, or the charge is overstated. All three are testable, by Faraday-cup charge, by size from witness-plate craters, and by ion spectra.
- Boring into ceramics and glass (talk L345) is what explosive-emission plasma does on dielectric surfaces. Mesyats §4.2 and §6.4 cover ectons and dielectric breakdown, and surface flashover leaves channels. The 41 µJ of kinetic energy per packet is enough to evaporate about 4 ng of oxide at about 10⁴ J g⁻¹, roughly a 1 µm × 1 mm channel. No anomaly is needed for micron channels.

### Gaps
- Shoulders' book "EV — A Tale of Discovery" (1987) and the 5,054,046, 5,123,039 and 5,153,901 patents were not opened, so no measured crater volumes or witness-plate data are quoted. The talk's claims about alumina and quartz boring and about atoms disappearing (talk L1104–1118) have no primary data here. **[unverified]**
- The mercury-wetted tungsten tip (talk L2390): no primary source was opened. Liquid-metal Taylor-cone ion sources are mainstream, but their link to EV stability is unverified.

## Q4. Usherenko super-deep penetration (SDP)

### Takeaway
SDP is a real materials-processing literature, mainly from Belarus (Usherenko, from 1970s onwards; J. Eng. Phys. Thermophys. 2002–2004). Powder particles of 10–100 µm fired at 1–3 km s⁻¹ leave traces far deeper than classical cratering predicts. The mechanism is debated (microjets, channel collapse, electromagnetic effects). The talk's "sand grains 30 cm into steel … with transmuted matter" was not found in any source opened. **[effect: contested but published; transmutation: unverified]**

### Cited Findings
- Crossref (metadata only): Usherenko, "Modern notions of the effect of superdeep penetration", J. Eng. Phys. Thermophys. 75, 753 (2002), DOI 10.1023/A:1016886300510. Usherenko, Koval' & Usherenko, "Formation of high-velocity microjets in different variants of superdeep penetration", 76, 109 (2003), DOI 10.1023/A:1022975510554. Same authors, "Estimation of the energy expended for superdeep penetration", 77, 641 (2004), DOI 10.1023/B:JOEP.0000036512.09918.b0.
- Nikitin, Nikitin, Velichko & Nikitina, Crystallography 69(5), 811–820 (2024). They propose an explanation for effects "at the collision of a stream of metal particles with a size of 10–100 microns, moving at a speed of 1–3 km/s". A particle charged to "about 10⁻⁹ C" presses on the channel wall "with a force of about 500 MPa", and the resulting combination "leads to a significant reduction in the particle-wall friction force" — [journal page](https://rjmseer.com/0023-4761/article/view/673716).

### Inferences
- The SDP effect is about anomalous tracks of micro-particles in metals. Explaining it needs no nuclear physics, and every proposed mechanism is mechanical or electromagnetic.

### Gaps
- No opened source gave penetration depths in calibres. The 10²–10⁴ figure and the 30 cm figure stay unverified. The ITMO and BNTU full texts would not download (one had a TLS certificate error, one PDF was truncated). Independent replications outside the Usherenko group were not found.

## Q5. "Cavity lights" and the "Freiberger" ball-lightning model

### Takeaway
"Freiberger" in the talk is David Fryberger of SLAC (died May 2026). In SLAC-PUB-6473 (1994) he proposed that ball lightning is a vorton–antivorton plasma powered by catalysed nucleon decay. The "cavity lights" were found at Jefferson Lab and first reported in 1999 (Delayen & Mammosser). In 1.5 GHz superconducting niobium cavities at 2 K, operated in field emission, glowing filaments and later small luminous objects (≤1.5 mm) were seen orbiting for seconds. The conventional explanations considered were trapped particulates, field-emitted charge or residual-gas ionisation, and Kapitza-type RF plasma filaments. Fryberger's own conventional small-particle model (NIM A 2001) produces the elliptical orbits but needs a particle "highly unlikely" to be "made of any known material". The joint 2009 NIM A paper calls the objects unexplained by established physics. No independent replication was found. **[observation: published, single group; vorton/nucleon-decay model: unverified, no mainstream support]**

### Cited Findings
- Fryberger, "A model for ball lightning", SLAC-PUB-6473 (October 1994), invited talk at the First International Workshop on the Unidentified Atmospheric Light Phenomena in Hessdalen: "The core … of BL in this model is comprised of a vorton-antivorton plasma. The energy of BL … derives from nucleon decay catalyzed by this plasma" — [INSPIRE PDF](https://inspirehep.net/files/6460d66897b84a2be962137f01570058).
- Delayen & Mammosser, "Light emission phenomena in superconducting niobium cavities", PAC 1999. The setup was a single-cell 1500 MHz cavity with a CCD camera in the beam pipe. "When operating the cavity in field emission, glowing filaments of light were observed trapped by RF fields in closed-orbit trajectories". Closed filaments were 3 cm to more than 7 cm across and oscillated at "a few tens of hertz"; open filaments moved at "a few meters per second". The effect was repeated on a CEBAF five-cell cavity. Conventional options discussed: particles heated on the surface and trapped by the RF, with the glow from "interaction with field emission, or by ionization of the residual helium gas". However, "It is not apparent … how the fields of a TM01 mode can produce the closed trajectories". Luminous spots had also been seen in room-temperature cavities deliberately contaminated with particles. Kapitza's 1958 RF filament rings are compared — [JACoW PDF](https://jacow.org/p99/papers/MOP113.pdf).
- Fryberger, "A small particle model as a possible explanation of recently reported cavity lights", NIM A 459, 29–43 (2001), DOI 10.1016/S0168-9002(00)00987-6, SLAC-PUB-8405 (Crossref-checked). At an excitation "equivalent to an average acceleration of ~4 MV/m", a polarisable small particle sits in a cylindrical harmonic-oscillator well that has elliptical orbits. However, "it is highly unlikely that a small object orbiting in such a potential well can be made of any known material" — [INSPIRE record](https://inspirehep.net/files/4626dc3f537a8ff9180f5d320804c4e1).
- Anthony, Delayen, Fryberger, Goree, Mammosser, Szalata & Weisend, "Experimental studies of light emission phenomena in superconducting RF cavities", NIM A 612, 1–45 (2009), DOI 10.1016/j.nima.2009.08.084 (Crossref-checked). There were nine runs, and the objects were "≤1.5 mm in size, freely moving about in the vacuum space". Up to seven objects moved together in one run. "these mobile luminous objects are without explanation within the realm of established physics" — [ORNL record](https://impact.ornl.gov/en/publications/experimental-studies-of-light-emission-phenomena-in-superconducti/).
- Fryberger, SLAC-PUB-13583 (April 2009), "A ball lightning model as a possible explanation of recently reported cavity lights". The objects were at 2 K with no beam and a very good vacuum; reflections in the walls show they are not touching the walls. The standard small-particle model "can explain some features" but is "incomplete", which motivated the vorton ball-lightning model — [INSPIRE PDF](https://inspirehep.net/files/5792f16537cf52514fffba8608775a5c).

### Inferences
- The talk's details are mixed up. "1998" should be 1999 (first report) to 2009, "superconducting chamber" is a superconducting RF cavity, and "DOE-funded" fits SLAC and JLab as DOE laboratories. "Baryon decay" matches Fryberger's 1994 model. The talk's "Shverblis and … Nevsky" line could not be traced to a source.
- Particulate contamination plus field emission is the default mainstream hypothesis for light in SRF cavities. Fryberger's analysis says that hypothesis fails quantitatively for these orbits, and no later group published a test. This is a genuinely open, unreplicated anomaly. It is not evidence for nucleon decay, because no nucleon-decay signature (pions, gammas, missing mass) was reported.

### Gaps
- The SLAC obituary appeared only in a search snippet, not opened. The full text of NIM A 612 (45 pp.) was not read, so its run-by-run speeds and sizes are not quoted. No replication or critique paper was found.

## Q6. Long-lived glowing objects on metals; Egely; Bogdanovich

### Takeaway
No primary source was found for these: plasmoids "living 2 days on metal" (talk L1494), free-electron-laser-hardened landing gear with glowing objects lasting two months, the 1988 Moscow-basement structure (Zverev, talk L2365), and plasmoid triads orbiting a common centre (talk L913–1000). They stay anecdotal. **[unverified]** Bogdanovich et al. (MEPhI) 2019 is covered in detail by another researcher. Cross-reference that note rather than this one.

### Cited Findings
- None located in primary literature within this search budget.

### Inferences
- The standard yardsticks for any such claim: unpowered electronic or plasma structures on metal decay in ns–µs (ectons stop in 10⁻⁹–10⁻⁸ s; Mesyats 1995). An object glowing for days needs a stored-energy source, and its power can be bounded from the light output. Its persistence then becomes a calorimetry question, which a calculator can frame.

### Gaps
- Egely's reactor, Bogdanovich 2019 (deferred to the other note), the FEL anecdote and the Zverev 1988 event all lack primary sources here.

## Q7. Iron-rich "crenellated" microspheres and "FeO2"

### Takeaway
Iron-rich microspheres are a known product of many high-temperature processes:
- lightning (iron fulgurites made of metallic Fe spheres; lightning spherules with Fe–Ti-oxide compositions inherited from the target rock),
- electrical discharge machining (α-Fe with some wüstite, hollow spheres at high energy),
- cosmic dust (I-type spherules of wüstite and magnetite with dendritic textures).

Dendritic, cellular or "orange-peel" surfaces are standard rapid-solidification textures. A 50 µm Fe droplet loses its latent heat in about 4–24 ms *(check)*. "FeO2" is a real compound, but only as a pyrite-type phase synthesised at 76 GPa and 1,800 K. A spherule with 93 wt% Fe (the talk's Adamenko example) can hold at most FeO₀.₂₆ even if all the remainder were oxygen. It is mostly metallic iron, not FeO₂ (63.6 wt% Fe) *(check)*. **[Fe microspheres: mainstream, ubiquitous; "FeO2 ratio" as a signature: refuted by composition; microspheres as a monopole signature: unverified]**

### Cited Findings
- Iron: melting point 1811 K, heat of fusion 13.81 kJ mol⁻¹, liquid density 6.98 g cm⁻³ at the melting point, atomic weight 55.845. Iron is about 5% of the crust's mass. FeO, Fe₃O₄ and Fe₂O₃ are the common oxides, all non-stoichiometric — [Wikipedia: Iron](https://en.wikipedia.org/wiki/Iron).
- Hu et al., "FeO2 and FeOOH under deep lower-mantle conditions and Earth's oxygen–hydrogen cycles", Nature 534, 241 (2016), DOI 10.1038/nature18018 (Crossref-checked). The press release describes "a highly stable, pyrite-structured iron oxide (FeO2) at 76 gigapascals and 1,800 kelvin" — [Carnegie](https://carnegiescience.edu/feo2-and-feooh-under-deep-lower-mantle-conditions-and-earths-oxygen-hydrogen-cycles).
- Boyd, Genge & Tomkins, "Microspherules formed by lightning strikes", Meteoritics & Planetary Science 61, 743–767 (2026), DOI 10.1111/maps.70123 (Crossref-checked). The spherules on Western Australian rock are "non-chondritic and are either dominated by silicates or Fe-Ti-bearing oxides". Their textures range "from glassy to coarse-grained", and their compositions match the underlying gneiss. The authors judge them "not compatible with either microtektites … or micrometeorites" — [abstract page](https://karmaka.de/?p=40532).
- The cited Cen 2014 spectrum shows Fe I lines in natural ball lightning, attributed to soil (Q1).
- Search excerpts only, not opened:
  - A University of Nebraska State Museum bulletin (1959) describes an "iron fulgurite … composed mainly of tiny spheres and filaments of metallic iron" — [UNL](https://digitalcommons.unl.edu/museumbulletin/13).
  - In EDM of iron, "besides the main fraction of α-Fe, a small percentage of wustite" forms in water, and none in kerosene; in mild-steel EDM debris, high energy input "showed formation of hollow spheres" — [IAEA INIS](https://inis.iaea.org/records/rg06v-ta916); [NML eprint](https://eprints.nmlindia.org/820).
  - The criteria for I-type cosmic spherules include "dendritic crystal morphology, the presence of wüstite and magnetite", and Ni-bearing phases. Terrestrial Fe spherules can show similar dendrites and voids, so Ni content is what distinguishes them — [cosmochemistry-papers summary](https://cosmochemistry-papers.com/2017/09/13/the-mineralogy-and-petrology-of-i-type-cosmic-spherules-implications-for-their-sources-origins-and-identification-in-sedimentary-rocks/).

### Inferences
- Fe mass fractions *(check)*: FeO 77.7%, Fe₃O₄ 72.4%, Fe₂O₃ 69.9%, FeO₂ 63.6%. An EDS spectrum with "Fe:O ≈ 1:2" by atom on a sub-micron oxide skin over metal would be an artefact of the interaction volume and the light-element quantification. It would not show a bulk phase that is unstable below tens of GPa.
- Lumped solidification time *(check)*:
  - Radiation alone (ε = 1, 1811 K): t = ρLd/(6εσT⁴), giving 24 ms for d = 50 µm.
  - Still-gas conduction (Nu = 2, assumed k_gas = 0.067 W m⁻¹ K⁻¹, ΔT = 1500 K): t = ρLd²/(12kΔT), giving 3.6 ms for 50 µm, and scaling as d², so 0.04 ms at 5 µm.
  - Millisecond solidification with large undercooling produces dendritic and cellular ("crenellated") surfaces and shrinkage holes. This is the expected texture for arc, lightning, cavitation-erosion and EDM droplets.
- The Adamenko hemisphere (93% Fe with Si) on a Cu target, Fe spheres from an Egely reactor, Fe spheres in a lightning slug, and SiO₂/CaO-coated Fe spheres (talk L700–790) all match the known products of arcs and melts on ordinary materials. Those materials carry Fe as a common impurity, and soil and dust are about 5% Fe. To count as a signature, a sample would need Fe that is absent from all inputs (shown by blank runs), or anomalous Fe isotope ratios. Neither is reported in the talk.

### Gaps
- The Adamenko patent image, the MFMP samples and the lightning slug were not examined. No published Fe-isotope data on such spheres were found. Cavitation-erosion spherule compositions were not searched.

## Q8. Tesla's "electrical particles of matter" and Wendt & Irion (1922)

### Takeaway
Wendt & Irion (JACS 1922) reported helium from tungsten wires exploded by a capacitor discharge. Allison & Harkins (JACS 1924) published a paper whose title states the null result: "The absence of helium from the gases left after the passage of electrical discharges: I, between fine wires in a vacuum". Rutherford commented in Nature in April 1922. The claim is regarded as refuted. **[refuted]** No primary Tesla text on "electrical particles of matter" was located. **[unverified]**

### Cited Findings
- Wendt & Irion, "Experimental attempts to decompose tungsten at high temperatures", JACS 44, 1887–1894 (1922), DOI 10.1021/ja01430a006. Wendt also published "The decomposition of tungsten", Science 55, 567 (1922). Allison & Harkins, JACS 46, 814–824 (1924), DOI 10.1021/ja01669a004. Rutherford, "Disintegration of elements", Nature 109, 418 (1922), DOI 10.1038/109418a0. Metadata for all four is Crossref-checked — [Crossref](https://api.crossref.org/works/10.1021/ja01669a004). Full texts are paywalled and were not opened.
- Briscoe, Robinson & Stephenson, "The electrical explosion of tungsten wires", J. Chem. Soc. 127, 240 (1925): metadata only — [RSC](https://pubs.rsc.org/en/content/articlelanding/1925/ct/ct9252700240). A DTIC exploding-wire review says the repeats found no helium (search excerpt only).

### Inferences
- In modern terms, an exploding wire deposits roughly 10⁴ J g⁻¹ (Mesyats's figure for electrical explosion), which is about 19 eV per tungsten atom *(check)*. Nuclear rearrangement needs MeV per nucleus. That is a shortfall of about 10⁵, which is why the null replications were expected.

### Gaps
- The quantitative helium yields (cm³ per mg of W) appear only on a secondary alchemy-history site, which was not opened, so no number is quoted. No modern repeat with mass spectrometry was found.

## Calculators for the crates

All test values are asserted in `ball_lightning_and_clusters_checks.py`. "Derived" means the value has no published counterpart and is checked against an independent algebraic form.

| Function signature | Formula | Primary source | Test value | Target crate::module | Status |
|---|---|---|---|---|---|
| `virial_energy_cap(volume, p0) -> f64` (exists) | 3p₀V | Finkelstein & Rubinstein 1964 | 20 cm → 1.273 kJ; 5 m → 19.9 MJ | plasmoid::relaxation | mainstream (exists) |
| `pixel_footprint(pitch, distance, focal) -> f64` | pitch·distance/focal | Cen et al. 2014 | 4 µm, 902 m, 3.3 mm → 1.093 m ("about 1.1 m") | metrology (or physica example) | mainstream |
| `sound_range(delay, c_sound) -> f64` | Δt·c | Cen et al. 2014 | 2.64 s × 340.51–341.15 → 898.9–900.6 m (0.898–0.9 km) | metrology | mainstream |
| `silicon_store(loading_kg_m3) -> f64` | loading·ΔH_f(SiO₂)/M_Si | Abrahamson & Dinniss 2000 | 0.1 kg m⁻³ → 3.24 MJ m⁻³ (paper 3.3) | plasmoid::discharge | mainstream |
| `arrhenius_ratio(ea, t1, t2) -> f64` | exp[(E_a/R)(1/T₁ − 1/T₂)] | Abrahamson & Dinniss 2000 | 189 kJ mol⁻¹, 1300→1400 K → 3.49 (lifetimes 10/3 s) | physica (shared) | mainstream |
| `ecton_current_density(h_bar, t_d) -> f64` | (h̄/t_d)^½ | Mesyats 1995 | 10⁹ A²s cm⁻⁴, 1 ns → 10⁹ A cm⁻² | plasmoid::discharge | mainstream |
| `ecton_electrons(i, a, h_bar, theta) -> f64` | i³/(256π²a²e h̄ sin⁴(θ/4)) | Mesyats 1995 Eq. 2.14 | no worked value in paper; use `packet_electrons` | plasmoid::discharge | mainstream (needs a test vector) |
| `ecton_erosion(rho, a, h_bar) -> f64` | (2/3)ρ(a/h̄)^½ | Mesyats 1995 Eq. 2.15 | model 0.6×10⁻⁴ g C⁻¹; measured (0.4–1.15)×10⁻⁴ | plasmoid::discharge | mainstream |
| `packet_electrons(current, duration) -> f64` | I·t/e | Mesyats 1995 | 1.6 A × 10 ns → 1.0×10¹¹ (portion 10¹¹–10¹²) | plasmoid::discharge | mainstream |
| `charge_cluster(n, diameter) -> (energy, potential, field, density)` | 3kQ²/5R, kQ/R, kQ/R², N/V | Shoulders US 5,018,180 inputs | 10¹¹ e, 1 µm → 2.77 J, 288 MV, 5.8×10¹⁴ V m⁻¹, 1.9×10²⁹ m⁻³ | plasmoid::fusion (next to `charged_sphere`) | derived; tests the EV claim |
| `coulomb_explosion_time(n, radius) -> f64` | √(mR³/kQe) = √3/ω_p | standard electrostatics | 10¹¹ e, 0.5 µm → 7.0×10⁻¹⁷ s | plasmoid::fusion | derived |
| `budker_neutralisation(beta) -> f64` | 1/γ² = 1 − β² | Budker condition (via Lampe review) | β = 0.1 → 0.99 (claimed ≤10⁻⁵; shortfall 9.9×10⁴) | plasmoid::pinch | mainstream; refutes the unneutralised EV |
| `brillouin_field(n) -> f64` | √(2mₑn/ε₀) | standard non-neutral plasma | 1.9×10²⁹ m⁻³ → 2.0×10⁵ T | plasmoid::pinch | derived |
| `energy_per_atom(specific_energy, molar_mass) -> f64` | w·M/F (eV) | Mesyats 1995 (10⁴ J g⁻¹) | W → 19.05 eV per atom | physica (examples) | derived; frames Wendt–Irion |
| `iron_mass_fraction(n_fe, n_o) -> f64` | nM_Fe/(nM_Fe + oM_O) | stoichiometry | FeO 0.777; Fe₃O₄ 0.724; Fe₂O₃ 0.699; FeO₂ 0.636; 93 wt% → ≤FeO₀.₂₆ | physica (examples) | mainstream; refutes "FeO2" for 93% Fe |
| `droplet_solidify_radiative(d, rho, latent, t, eps) -> f64` | ρLd/(6εσT⁴) | Fe data (Wikipedia) | 50 µm Fe, ε = 1 → 23.6 ms | physica (examples) | derived |
| `droplet_solidify_conductive(d, rho, latent, k_gas, dt) -> f64` | ρLd²/(12kΔT) | Nu = 2 sphere | 50 µm, k = 0.067, ΔT = 1500 K → 3.58 ms | physica (examples) | derived; k_gas is an assumption |
