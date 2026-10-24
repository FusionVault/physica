# Gravity-control and "breakthrough propulsion" programmes: the public record, the tests, and the calculable physics

Legend for every claim: **[mainstream]** means accepted physics or uncontested documentation; **[contested]** means published but disputed or not replicated; **[refuted]** means tested and shown false or explained as an artefact. **(doc)** marks a documented fact from a primary or near-primary source. **(rumour)** marks a secondhand or unverified claim. Numbers marked "computed" are my own arithmetic from the formulas the sources give, checked with a short Python script. They are not quoted from the sources.

## 1. 1950s electrogravitics: the 1956 report, Martin/RIAS, the industry programmes, T. T. Brown, and what modern tests showed

### Takeaway
The 1950s "gravitics" boom is documented, but almost all of it comes from trade-press and newsletter claims: Aviation Report, Aviation Studies' GRG-013/56, and the Herald Tribune. What is documented as real is unclassified general-relativity research funded by RIAS and the USAF. Every modern test of Brown's asymmetric capacitors (Talley 1991, Bahder & Fazi 2003, Canning et al. 2004, Tajmar 2004) finds the thrust is ion drift / corona wind in air, with no force in vacuum. The working formula is F = I·d/μ.

### Cited Findings
**The "Electrogravitics Systems" report (1956)**
- (doc) The title page reads: "Electrogravitics Systems: An examination of electrostatic motion, dynamic counterbary and barycentric control", prepared by the "Gravity Research Group, Special Weapons Study Unit, Aviation Studies (International) Limited, London, 29-31 Cheval Place, Knightsbridge", Report GRG-013/56, February 1956. Its contents are Introductory Notes, Discussion, Conclusions, App. I "Aviation Report Extracts" and App. II "Electrostatic Patents". The copy carries a "Technical Library, Wright-Patterson Air Force Base" stamp. — [Valone (ed.), *Electrogravitics Systems: Reports on a New Propulsion Methodology*, Integrity Research Institute, 1st ed. 1994 (archive.org scan)](https://archive.org/details/thomas-valone-electrogravitics-systems)
- (doc, describing claims) The report defines electrogravitics as "a synthesis of electrostatic energy use for propulsion ... and gravitics, or dynamic counterbary". Its company claims are compiled from the trade newsletter *Aviation Report*: "Convair is extensively committed to the work with several rigs"; Lear Inc. and the Sperry division of Sperry-Rand "have a division ... working on gravity research"; Curtiss-Wright, Lockheed, Boeing and North American "have not yet declared themselves". It also says "Glenn Martin say gravity control could be achieved in six years, but ... it would entail a Manhattan District type of effort". The *Aviation Report* of 15 Nov 1955 named Glenn Martin, Convair, Sperry-Rand, Sikorsky, Bell, Lear Inc. and Clark Electronics. — [same scan](https://archive.org/details/thomas-valone-electrogravitics-systems). **[contested/unverified]**: these are journalists' statements that companies were "studying" gravity. They are not evidence of hardware.
- (rumour) The editor Valone says the word "CONFIDENTIAL" was written under a black streak on the original title page, and that Paul LaViolette retrieved the copy from the Wright-Patterson technical library, although it is "not actually listed on their library computer". — [Valone scan](https://archive.org/details/thomas-valone-electrogravitics-systems). Fringe sites give a reference number "13-1-00034-5879" and "declassified 1990", which I could not check against any library catalogue. — [search summary of hobbyist pages, e.g. geocities mirror](https://www.geocities.ws/marksrealm/project080.html)
- (rumour) LaViolette's claim that 1950s electrogravitics led to a B-2 auxiliary propulsion system appears in the same compilation ("The U.S. Antigravity Squadron", 1993). It has no independent confirmation. **[refuted/unsupported]** — [Valone scan](https://archive.org/details/thomas-valone-electrogravitics-systems)

**Glenn L. Martin / RIAS and the industry programmes (historians' view)**
- (doc) RIAS was founded in Baltimore in 1955 by George Bunker of the Glenn L. Martin Company. It was later renamed Martin Marietta Laboratories and dropped basic research in 1973. — [Wikipedia: Research Institute for Advanced Studies](https://en.wikipedia.org/wiki/Research_Institute_for_Advanced_Studies)
- (doc, oral history) Louis Witten said at the 27th Texas Symposium roundtable (2013) that he was recruited by George S. Trimble, Martin's VP for Aviation and Advanced Propulsion, for work Trimble framed as "the discovery of anti-gravity". Witten's RIAS-era output was ordinary general relativity, e.g. Phys. Rev. 113, 357; Phys. Rev. 115, 206; Rev. Mod. Phys. 34, 442. TWZ found no record confirming Witten's claimed Wright Field contract. — [The War Zone (TWZ), "The Truth Is The Military Has Been Researching 'Anti-Gravity' For Nearly 70 Years"](https://www.twz.com/30499/the-truth-is-the-military-has-been-researching-anti-gravity-for-nearly-70-years)
- (doc) Ansel Talbert's three-part New York Herald Tribune series ran 20–22 November 1956: "Conquest of Gravity Aim of Top Scientists in U.S.", "Space-Ship Marvel Seen If Gravity is Outwitted" and "New Air Dream – Planes Flying Outside Gravity". It quoted William Lear, Lawrence Bell, Igor Sikorsky, George Trimble and others. TWZ's own caption dates the front page to 1955, which conflicts with its text. — [TWZ](https://www.twz.com/30499/the-truth-is-the-military-has-been-researching-anti-gravity-for-nearly-70-years)
- (doc, historian) From September 1956 the USAF funded gravitation research through the General Physics Laboratory of the Aeronautical Research Laboratories at Wright-Patterson AFB. The funding ended with the Mansfield Amendment (1969/70). Goenner writes: "To my knowledge, the supported projects were not classified." He draws on J. Goldberg, "United States Air Force support of general relativity: 1956–1972", *Einstein Studies* 3, 89–102 (1988). — [Goenner, arXiv:1607.03319](https://arxiv.org/abs/1607.03319). Joshua Goldberg ran that laboratory from 1956 to 1962. — [TWZ](https://www.twz.com/30499/the-truth-is-the-military-has-been-researching-anti-gravity-for-nearly-70-years) **[mainstream]**
- (doc) The 1972 USAF study "Advanced Propulsion Concepts – Project Outgrowth" (ed. F. Mead, for AFRPL Edwards) concluded that "all of the ideas discussed lack theoretical and technical merit". — [TWZ quoting DTIC AD-750554](https://apps.dtic.mil/dtic/tr/fulltext/u2/750554.pdf)

**T. Townsend Brown and the Biefeld–Brown effect**
- (doc) Brown's patents: GB 300,311 (1928) "A method of and an apparatus or machine for producing force or motion"; US 1,974,483 (1934); US 2,949,550 (1960) "Electrokinetic apparatus"; US 3,018,394 (1962); US 3,022,430 (1962); US 3,187,206 (1965); US 3,196,296 (1965). — [Wikipedia: Thomas Townsend Brown](https://en.wikipedia.org/wiki/Thomas_Townsend_Brown). Bahder & Fazi note that the 1960 patent "makes no reference to gravitational effects", while the 1965 patent claims a force persists in vacuum. — [Bahder & Fazi, arXiv:physics/0211001](https://arxiv.org/abs/physics/0211001)
- (rumour) A physicist (not named) who watched Brown's disc demonstration in the early 1950s put it down to "electric wind". Brown's claimed vacuum demonstration at SNCASO in France (1955–56) has no published data. — [Wikipedia: T. T. Brown, citing Schatzkin, *The Man Who Mastered Gravity*](https://en.wikipedia.org/wiki/Thomas_Townsend_Brown)
- (doc) De Seversky's "Ionocraft" patent (US 3,130,945, 1964) reports 50 mN from a 0.1 m² electrode array with a 5 cm gap at 20 kV and 0.5 mA. It notes the thrust varied with humidity and pressure, "exactly what one would expect from a corona wind". — [Tajmar 2004 preprint (TU Dresden)](https://tu-dresden.de/die_tu_dresden/fakultaeten/fakultaet_maschinenwesen/ilr/rfs/forschung/folder.2007-08-21.5231434330/ag_raumfahrtantriebe/Biefeld-Brown%20Effect%20AIAA%20Journal%20Revised.pdf)

**Modern tests: Army Research Laboratory, Bahder & Fazi (ARL-TR-3005, June 2003)**
- (doc) The authors built four "lifters". Lifter 1 (~5 g): ~37 kV, ~1.5 mA, failed to lift. Lifter 2: 6 cm gap, 30 kV, ~1.5 mA, "hovered vigorously". The wire electrode radius was r1 = 0.0025 in and the foil-edge radius r2 = 0.125 in. — [arXiv:physics/0211001](https://arxiv.org/abs/physics/0211001); [official ARL-TR-3005 at GovInfo](https://www.govinfo.gov/content/pkg/GOVPUB-D101-PURL-LPS125022/pdf/GOVPUB-D101-PURL-LPS125022.pdf)
- (doc) **Ballistic ion-wind upper bound (their eq. 3):** F = I·√(2mV/q). With electrons, V = 40 kV and I = 1 mA it lifts M = 6.8×10⁻⁵ g (eq. 5). With Cu⁺ ions they report 0.002 g (eq. 6). Either way it is "at least three orders of magnitude too small". — [arXiv:physics/0211001](https://arxiv.org/abs/physics/0211001)
- (doc) **Ion-drift formula (their eqs. 10–14):** F = Q·V/l and I = Q·v/l, so F = I·V/v. With v = μE and E = V/l this gives **F = I·l/μ**. The mobilities are μ_e = 620 cm²/(V·s) at E = 10⁴ V/cm and μ(N₂ ions) = 2.5 cm²/(V·s). With I = 10⁻³ A and l = 0.04 m: F = 6.4×10⁻⁴ N for electrons (lifts 0.064 g, eqs. 15–16) and **F = 0.16 N for N₂ ions** (eq. 17, ~16 g). This is the right order of magnitude for a lifter. Yet the authors concluded that "the physical basis for the Biefeld-Brown effect is not understood", and they left vacuum testing as future work. — [arXiv:physics/0211001](https://arxiv.org/abs/physics/0211001)
- (doc) In vacuum, their force is F = (ε₀/2)[∮_{S+} E² n dS + ∮_{S−} E² n dS] (eq. 7). — [arXiv:physics/0211001](https://arxiv.org/abs/physics/0211001)

**NASA: Canning, Melcher & Winet, NASA/CR-2004-213312, "Asymmetrical Capacitors for Propulsion" (ISR Fairmont WV, Oct 2004, contract NAS3-00124)**
- (doc) **Model (eqs. 5–9):** F = qE = qV/d = tIV/d. They assume 10¹⁰ collisions/s at 1 atm (each one stops the ion), which gives **F = 2·d·(10¹⁰ s⁻¹)·m·I/e**. With d = 7 cm, m(N₂) = 4.7×10⁻²⁶ kg and I = 3.5 mA, F = 1.44 N (a stated maximum). Averaging over field-line paths adds about 38%. The impulse per ion is ∫dh/α, i.e. proportional to the vertical extent of its path (eqs. 10–12). — [NTRS 20040171929 (PDF)](https://ntrs.nasa.gov/api/citations/20040171929/downloads/20040171929.pdf)
- (doc) Data in Table 1, at ~1 atm with the cylinder grounded: Device 3/A ran at 40 kV, 220 µA across the ACT (240 µA total), 110 RPM. Device 3/D: 39 kV, 89 µA, 90 RPM. Device 4/A: 42 kV, 130 µA, 128 RPM. Device 4/D: 46 kV, 307 µA, 138 RPM. The thrust direction did not depend on polarity, and similar forces appeared in N₂ and argon. — [NTRS PDF](https://ntrs.nasa.gov/api/citations/20040171929/downloads/20040171929.pdf)
- (doc) **Vacuum null:** they pumped below 10⁻⁴ Torr and raised the pressure in steps. "No device showed signs of rotation at a pressure less than 300 Torr", with one exception: a single 1/8-turn jerk during an arc at 44 kV and 5.5×10⁻⁵ Torr, which they attributed to ejected material or outgassing. Their conclusion: "we find no evidence to support [new physics] ... their operation is fully explained by a very simple theory that uses only electrostatic forces and the transfer of momentum by multiple collisions." **[refuted: electrogravitic interpretation]** — [NTRS PDF](https://ntrs.nasa.gov/api/citations/20040171929/downloads/20040171929.pdf)

**Other in-vacuum or isolated nulls**
- (doc) Talley, "Twenty First Century Propulsion Concept", PL-TR-91-3009 (Veritay Technology for USAF Phillips Lab, contract F04611-89-C-0023, Feb 1989–Jul 1990): torsion-pendulum tests in vacuum. "Direct experimental results show that under high vacuum conditions… no detectable propulsive force was electrostatically induced by applying a static potential difference… between test device electrodes". Also: "If such a force still exists and lies below the threshold of measurements in this program, then the force may be too small to be attractive for many, if not most, space propulsion applications" (p. 91). — [quoted in Canning et al., NTRS PDF](https://ntrs.nasa.gov/api/citations/20040171929/downloads/20040171929.pdf); [WorldCat record](https://search.worldcat.org/oclc/227770672). The tests used voltages up to ~19 kV at ~1 µTorr, and Talley saw a transient force only during electrical breakdown. — [search summary of report abstract](https://search.worldcat.org/oclc/227770672); [Tajmar 2004](https://tu-dresden.de/die_tu_dresden/fakultaeten/fakultaet_maschinenwesen/ilr/rfs/forschung/folder.2007-08-21.5231434330/ag_raumfahrtantriebe/Biefeld-Brown%20Effect%20AIAA%20Journal%20Revised.pdf)
- (doc) Tajmar, "The Biefeld-Brown Effect: Misinterpretation of Corona Wind Phenomena", *AIAA Journal* 42(2), 315–318 (2004). He ran a corona discharge inside a sealed, grounded Faraday box on a pendulum, with laser position readout of ±0.1 µm. The force accuracy was **10.5 µN**, at 0–20 kV. "No thrust was detected." Any anomalous force "must be at least five orders of magnitude below corona wind phenomena". — [TU Dresden preprint](https://tu-dresden.de/die_tu_dresden/fakultaeten/fakultaet_maschinenwesen/ilr/rfs/forschung/folder.2007-08-21.5231434330/ag_raumfahrtantriebe/Biefeld-Brown%20Effect%20AIAA%20Journal%20Revised.pdf) (Millis's TM cites it with a slightly different journal label, "AIAA Journal of Propulsion & Power, 42".)
- (doc) Two earlier AFAL studies: the 1988 Veritay "21st Century Propulsion Concept" (DTIC ADA197537) found ion effects "negligible". The 1989 "Electric Propulsion Study" (DTIC ADA227121) reviewed ionic wind and the Mach effect. — [TWZ](https://www.twz.com/30499/the-truth-is-the-military-has-been-researching-anti-gravity-for-nearly-70-years)

### Inferences
- **Rust-ready ion-wind formulas.**
  - Drift thrust is `F = I*d/mu`. Thrust per power is `F/P = d/(mu*V)`, because P = IV.
  - Canning's collision model implies a mobility mu = e/(2·m·ν) ≈ 1.7×10⁻⁴ m²/(V·s) for N₂ with ν = 10¹⁰ s⁻¹ (computed). That is consistent with Bahder's 2.5×10⁻⁴. The two "different" NASA and ARL formulas are therefore the same law.
  - Computed checks:
    - Bahder eq. 17 gives 0.160 N.
    - Canning eq. 9 gives 1.44 N.
    - De Seversky's 0.5 mA over 5 cm predicts ~0.10–0.125 N, against the 50 mN measured: right order, with inefficiency.
    - Bahder lifter 2 at 30 kV, 1.5 mA and 6 cm: F ≈ 0.36 N at P = 45 W, so F/P = d/(μV) ≈ 8 mN/W (with μ = 2.5×10⁻⁴).
  - Good unit tests: `ballistic_upper_bound(I, V, m, q) = I*sqrt(2*m*V/q)` must be ≪ `drift(I, d, mu)`. In vacuum, μ → ∞ means F → 0.
- **Arithmetic slip in Bahder & Fazi.** Their Cu⁺ ballistic number in eq. 6 (0.002 g) recomputes to ≈0.023 g (2.3×10⁻⁴ N; computed). That is about 2 orders below a 3–5 g lifter, not 3. The electron value in eq. 5 (6.8×10⁻⁵ g) checks out. Their conclusion is unchanged.
- The 1950s record supports the claim that aerospace executives publicly expressed interest, and that the USAF and industry funded unclassified GR. It does not support any classified hardware programme. The only "classified" marking claimed for GRG-013/56 rests on the editor's note in a fringe compilation.

### Gaps
- I could not open the full text of PL-TR-91-3009, so Talley's torsion-pendulum force resolution in N is unknown.
- I found no primary document for "Project Winterhaven" (Brown's 1952 proposal) or for a dated ONR/Navy evaluation of Brown. Only the biography-based account of an "electric wind" verdict surfaced.
- Goldberg's 1988 *Einstein Studies* chapter, and any historian's account specifically of Convair, Bell or Lear gravity "rigs", were not retrieved. Nick Cook's *The Hunt for Zero Point* (2001) was not consulted directly, so treat its claims as unverified.

## 2. NASA Breakthrough Propulsion Physics (BPP, 1996–2002) and *Frontiers of Propulsion Science* (2009)

### Takeaway
BPP ran three "grand challenges": propellantless propulsion, hyper-fast travel and breakthrough power. Of 16 tasks (BPP's and others'), 6 were found not viable, 6 were unresolved and 4 merited sequels. Biefeld–Brown and the rotating-superconductor (Podkletnov) approach were explicitly judged non-viable. Woodward's transient-inertia work was left "unresolved". The 2009 AIAA book compiles this work.

### Cited Findings
- (doc) Millis, "Prospects for Breakthrough Propulsion From Physics", NASA/TM-2004-213082 (May 2004, NASA Glenn). BPP sponsored 8 tasks (5 competitive, 2 in-house, 1 minor grant), which produced 13 peer-reviewed articles. "Of the 16 specific tasks ... 6 were found not to be viable, 6 remain unresolved or have debatable findings, and 4 have clear opportunities for sequels." — [NTRS 20040070788](https://ntrs.nasa.gov/citations/20040070788) ([PDF](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20040070788.pdf))
- (doc) Task verdicts in the TM — [NTRS PDF](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20040070788.pdf):
  - **Schlicher thruster** (a specially terminated coax): "Tests observed no such thrust". **[refuted]**
  - **Cavendish test of superconductor claims** (Robertson, AIAA-2001-3364): the RF "coupled too strongly to supporting instrumentation and prevented any discernable results". The result was inconclusive and no sequel was planned.
  - **Woodward transient inertia:** Woodward's later papers "indicated that the effect was much smaller than originally reported". The independent test (Cramer, AIAA-2001-3908) was "unable to resolve any discernable effect". "This transient inertia approach is considered unresolved." **[contested]**
  - **EM torsion theory:** the test was null, but "missed a critical characteristic", so it stays unresolved.
  - **Superluminal tunnelling:** "the entire signal is still light-speed limited". Not viable. **[mainstream]**
  - **Vacuum energy (Casimir, MEMS cavities):** net propulsive forces are possible in principle "even thought [sic] the forces are impractically small". This line has sequels.
  - **Podkletnov:** NASA MSFC tried a full replication from 1995 to 2002 but could not finish the hardware. The privately funded Hathaway, Cleveland & Bao replication (*Physica C* 385 (2003) 488–500) "found no evidence of a gravity-like force to the limits of the apparatus sensitivity", a sensitivity "50 times better than that available to Podkletnov". The TM concludes: "this rotating, RF-pumped superconductor approach is considered non-viable." **[refuted at that sensitivity]**
  - **Podkletnov "force-beam":** the claim is 4×10⁻⁴ to 23×10⁻⁴ J imparted to a distant 18.5 g pendulum. "Unsubstantiated by reliable independent sources." **[contested/unsubstantiated]**
  - **Slepian drive / Abraham–Minkowski** (ISR, Congressional earmark): an Air Force Academy assessment "concluded that no net propulsive forces are expected". Brito's device data were "not sufficiently above the noise". Unresolved.
  - **Biefeld–Brown:** "all rigorous experimental tests indicate that the observed thrust is attributable to ion wind" (citing Talley PL-TR-91-3009 and Tajmar 2004). Lifters and ACTs "are not viable candidates for breakthrough physics propulsion." **[refuted]**
  - **ESA gravity-control study:** it recommended equivalence-principle tests, a probe for the Pioneer anomaly, and gravitomagnetic experiments in quantum materials.
- (doc) The TM quotes Aviation Week: "At least one large aerospace company is embarking on ZPE (quantum vacuum) research in response to a Defense Dept. request". It adds that "Given the private and protected nature of such sponsorship, it is not known to what extent these results will be disseminated." — [NTRS PDF](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20040070788.pdf)
- (doc) *Frontiers of Propulsion Science* is Progress in Astronautics and Aeronautics vol. 227 (AIAA, Reston VA, 2009; 737 pp.; ISBN 978-1-56347-956-4), edited by M. G. Millis and E. W. Davis, with 18 contributing authors. — [Aeronautical Journal review (Cambridge)](https://www.cambridge.org/core/journals/aeronautical-journal/article/frontiers-of-propulsion-science-progress-in-astronautics-and-aeronautics-series-vol-227-edited-by-m-g-millis-and-e-w-davis-american-institute-of-aeronautics-and-astronautics-1801-alexander-bell-drive-suite-500-reston-va-201914344-usa-2009-737pp-illustrated-9995-aiaa-members-12995-nonmembers-isbn-9781563479564/10715D1687BAC6184386F0E40904EB06); [NSS review](https://nss.org/frontiers-of-propulsion-science-by-marc-g-millis-and-eric-w-davis/). The table of contents includes Gilster's history, Frisbee on the limits of interstellar flight, Millis's "Prerequisites for space drive science", Davis's review of gravity control, chapters on superconductor experiments and "antigravity" devices, and a "null-results review of an electrogravity patent". — [UTB library TOC](https://vufind.katalog.k.utb.cz/Record/73788/TOC)
- (doc) Millis's companion IAC-10 paper "Progress in Revolutionary Propulsion Physics" summarises the book's assertions. It sits behind a paywall. — [IAF directory](https://iafastro.directory/iac/archive/browse/IAC-10/C4/8/8526/)

### Inferences
- BPP's verdicts line up with every later test. The concepts BPP left "unresolved" (Woodward and EmDrive-like microwave cavities) were later closed as artefacts by the Dresden SpaceDrive group (sections 5 and 6).
- For code: the photon-rocket bound F/P = 1/c is the benchmark BPP and Dresden both use to judge "propellantless" claims (section 6).

### Gaps
- BPP's total budget is not in the TM. Figures circulate (around $1.2M), but I did not verify them.
- I did not obtain the chapter-by-chapter conclusions or chapter authors of *Frontiers of Propulsion Science*.

## 3. BAE Systems "Project Greenglow" (mid-1990s to ~2005)

### Takeaway
Greenglow was a small, company-funded BAE study led by Ron Evans. It was publicly described as early as 2002–2004, so it was not first disclosed in 2017. It funded assessments of Podkletnov's claims (null), microwave-thruster work and gravito-electromagnetic theory. By the accounts available it produced no breakthrough. Public profile rose with Evans's 2015 book and BBC Horizon in March 2016.

### Cited Findings
- (doc) NASA's 2004 TM describes it: "Project Greenglow, British Aerospace Systems... sponsored a modest project to look at a variety of breakthrough propulsion approaches. Headed by Dr. Ron Evans, incremental research tasks were supported". The tasks were:
  - assessments of Podkletnov's gravity-shield claims ("null findings"; ref. 70 is R. C. Woods, "Review of Claims of Interaction Between Gravitation and High-Temperature Superconductors", AIP Conf. Proc. 699, STAIF 2004, pp. 1085–1092);
  - microwave thrusters (P. D. Smith & Vinogradov, "Microwave Thrust Mechanism", AIAA-2001-3802);
  - theoretical gravitation (Tucker & Clark, "Gauge Symmetry and Gravito-Electromagnetism", *Class. Quantum Grav.* 17 (2000) 4125; Burton et al., "Towards the Control of Matter with Gravity", AIAA-2001-3912);
  - vacuum forces and "what-if" assessments.

  "It is not known if, or at what level, this project will continue." — [NASA/TM-2004-213082](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20040070788.pdf)
- (doc) A Sheffield replication team (Woods et al., "Gravity Modification by High Temperature Superconductors", AIAA 2001-3363) never reached the conditions Podkletnov specified when he visited in 2000. — [Wikipedia: Eugene Podkletnov](https://en.wikipedia.org/wiki/Eugene_Podkletnov). Note: the link between that Sheffield work and BAE funding is my inference from Millis's ref. 70. It is not stated in the sources I read.
- (doc) The book is Ronald Evans, *Greenglow: & the Search for Gravity Control*, Matador, April 2015, 320 pp., ISBN 978-1-78462-023-3. Evans "persuaded BAE Systems to sponsor a small programme of research to investigate some ideas linked with gravity control." — [Schweitzer Fachinformationen listing](https://www.schweitzer-online.de/buch/Evans/Greenglow/9781784620233/A24255982/); [NetGalley](https://www.netgalley.co.uk/catalog/book/69589)
- (doc) BBC Horizon, "Project Greenglow – The Quest for Gravity Control", aired 23 March 2016, narrated by Peter Capaldi. The BBC blurb says BAE began the project "in the mid 1990s". The programme covered Podkletnov and the EmDrive. — [Centauri Dreams (Millis), 2016](https://centauri-dreams.org/?p=35243)
- (rumour/secondary) Evans reportedly wrote on retiring in 2005 that no major breakthroughs had been made. In 2016 BAE declined to discuss the project. — [search summary of Inverse article](https://inverse.com/article/13238-if-project-greenglow-worked-engineers-would-have-built-a-flying-saucer)
- (doc) Jane's (July 2002) reported that a Boeing GRASP briefing said Boeing, BAE Systems and Lockheed Martin had tried to approach Podkletnov directly. — [Wikipedia: Eugene Podkletnov](https://en.wikipedia.org/wiki/Eugene_Podkletnov)

### Inferences
- The premise that Greenglow was "made public in 2017 or later" does not hold up. It was in Jane's (2002) and a NASA TM (2004), and Evans's book (2015) and BBC Horizon (2016) followed. No BAE or government release of Greenglow records was found.
- Greenglow's measurable outputs are the same Podkletnov nulls and microwave-thruster work covered in sections 4 and 6. Nothing Greenglow-specific has published numbers.

### Gaps
- I found no BAE primary document, budget, or list of funded universities and contracts. The Sheffield, Lancaster and Smith & Vinogradov links to BAE money are inferred from Millis's reference list.

## 4. Boeing GRASP (2002) and the Podkletnov replication interest

### Takeaway
GRASP ("Gravity Research for Advanced Space Propulsion") is known from a single Jane's Defence Weekly report by Nick Cook (July 2002), based on a leaked internal document. Boeing denied funding it with company money but would not comment on "black projects". Podkletnov's claims (a 0.05–0.3% weight loss in 1992, 2% in a withdrawn 1996 paper) failed a full replication at 50× better sensitivity. **[refuted]**

### Cited Findings
- (doc) Podkletnov & Nieminen, "A possibility of gravitational force shielding by bulk YBa₂Cu₃O₇₋ₓ superconductor", *Physica C* 203(3–4), 441–444 (10 Dec 1992), doi:10.1016/0921-4534(92)90055-H. Samples were reported to lose "from 0.05 to 0.3% of [their] weight, depending on the rotation speed". — [Wikipedia: Eugene Podkletnov](https://en.wikipedia.org/wiki/Eugene_Podkletnov); [search summary](https://arxiv.org/pdf/cond-mat/0303647v1)
- (doc) A 1996 *J. Phys. D* submission claiming 2% was withdrawn after acceptance. The co-author Vuorinen disavowed it. A 1997 preprint (arXiv:cond-mat/9701074) is believed to be substantially the same paper. — [Wikipedia: Eugene Podkletnov](https://en.wikipedia.org/wiki/Eugene_Podkletnov)
- (doc) Hathaway, Cleveland & Bao, "Gravity modification experiment using a rotating superconducting disk and radio frequency fields", *Physica C* 385 (2003) 488–500: no gravity-like force at a sensitivity 50× better than Podkletnov's. — [NASA/TM-2004-213082](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20040070788.pdf)
- (doc, journalism) Nick Cook in Jane's Defence Weekly (July 2002) quoted the GRASP briefing: "If gravity modification is real, it will alter the entire aerospace business." The briefing focused on Podkletnov's "impulse gravity generator" and its possible use as a weapon. The project was attributed to Phantom Works, Seattle, and Jane's reported that Phantom Works head George Muellner found the Russian work convincing. — [Wikipedia: Eugene Podkletnov](https://en.wikipedia.org/wiki/Eugene_Podkletnov); [2002 mailing-list repost of Jane's/BBC text](https://pages.uoregon.edu/stevev/sd-archive/sd-2002-3/msg00011.html); [wissenschaft.de](https://wissenschaft.de/artikel/project-grasp-boeing-will-schwerkraft-teilweise-aufheben)
- (doc) The BBC News story "Boeing tries to defy gravity" (July 2002) said the work was run by "the top-secret Phantom Works". — [TechSpot repost](https://www.techspot.com/community/topics/boeing-tries-to-defy-gravity.2242/). The sources disagree on its date: Wikipedia gives 7 July 2002, while another repost dates the Jane's piece 29 July 2002.
- (doc) Space.com (Jim Banke, 31 July 2002) reported that Boeing denied funding GRASP with company money and said it could not comment on "black projects". — [Wikipedia: Eugene Podkletnov](https://en.wikipedia.org/wiki/Eugene_Podkletnov)
- (doc) The Podkletnov & Modanese "impulse gravity generator" (arXiv physics/0108005) claims 4×10⁻⁴ to 23×10⁻⁴ J delivered to an 18.5 g pendulum. — [NASA/TM-2004-213082](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20040070788.pdf); [ar5iv physics/0108005](https://ar5iv.arxiv.org/html/physics/0108005) **[contested/unsubstantiated]**

### Inferences
- Computed: an energy E given to the 18.5 g pendulum implies an initial speed v = √(2E/m) of ≈0.21 to 0.50 m/s. That is large and easy to measure, so an independent replication would be cheap to judge. None has appeared in the sources reviewed.
- For code: a "shielding fraction" claim is directly testable with a gravimeter. A 0.05% weight change is 4.9×10⁻³ m/s² ≈ 490 mGal, enormous next to a modern gravimeter's µGal noise. Computed: 0.0005 × 9.81 m/s².

### Gaps
- The original Jane's article and the GRASP document are not publicly available. Only reposts exist.
- The NASA MSFC static test (Li et al., *Physica C* 281, 1997) did not surface with numbers.

## 5. Mach-effect thrusters (Woodward): predicted-thrust formula and the SpaceDrive nulls

### Takeaway
Woodward's theory predicts a mass fluctuation δm₀ ∝ (1/4πGρ₀c²)·∂P/∂t. Combined with electrostriction it gives a time-averaged thrust ⟨F⟩ = ω⁶m₀K_p²K_e x₀³V₀⁴/(8πGρ₀c²), which is ~8 nN for Woodward's own PZT stack. The Dresden SpaceDrive group, with balances reaching <5 nN noise after averaging, found the claimed µN thrusts are vibration, thermal and liquid-metal-contact artefacts. **[refuted as tested; theory contested]**

### Cited Findings
- (doc) **Woodward mass-fluctuation formula** (Fearn & Woodward, "Experimental Null test of a Mach Effect Thruster", arXiv:1301.6178, 2013, eq. 1):
  δm₀ ≈ (1/4πG)·[(1/ρ₀c²)·∂P/∂t − (1/ρ₀c²)²·P²/V], where P is the power into an accelerated object of volume V and density ρ₀. — [arXiv:1301.6178](https://arxiv.org/abs/1301.6178)
- (doc) Their eq. 2 simplifies this to δm₀ ≈ (1/4πGρ₀c²)·F·a ≈ m₀a²/(4πGρ₀c²): "the simplest Mach effect depends on the square of the acceleration". — [arXiv:1301.6178](https://arxiv.org/abs/1301.6178)
- (doc) **Thrust prediction** (eqs. 3–10):
  - Drive: i = i₀cos ωt and V = V₀cos ωt at resonance.
  - Piezo length: x = x₀(1 + K_pV), giving δm₀ = (m₀K_p²x₀²V₀²ω⁴/8πGρ₀c²)(1 − cos 2ωt) (eq. 6).
  - Electrostriction: x = x₀(1 + K_eV²), with acceleration ẍ = −2ω²K_e x₀V₀² cos 2ωt (eq. 8).
  - Time-averaged thrust: **⟨F⟩ = ω⁶m₀K_p²K_e x₀³V₀⁴/(8πGρ₀c²)** (eq. 10).

  Example: f = 32.35 kHz (ω⁶ = 7.05×10³¹), m₀ = 0.042 kg, x₀ = 19 mm, V₀ = 100 V, ρ₀ = 7.9×10³ kg/m³. This gives F ≈ 1.87×10²¹·K_p²K_e. With K_p = 3.2×10⁻¹⁰ m/V and K_e ≈ 4×10⁻¹¹ m/V², "a thrust on the order of 8 nN is predicted", while "the observed thrust ... is a couple of microNewtons". The authors attribute the difference to unmodelled mechanical resonance. — [arXiv:1301.6178](https://arxiv.org/abs/1301.6178)
- (doc) Woodward's torsion-balance results show "a particular force trace amounting up to 100 µN" at a resonance of 36 kHz and about 30 W. — [Monette, Kößling, Neunzig & Tajmar, "The SpaceDrive Project – Mach-Effect-Thruster Experiments on High-Precision Balances in Vacuum", Space Propulsion 2020+1, SP2020_00278 (Mar 2021)](https://tu-dresden.de/ing/maschinenwesen/ilr/rfs/ressourcen/dateien/forschung/folder-2007-08-21-5231434330/ag_raumfahrtantriebe/SPC-Mach-Effect-Thruster-Experiments-on-High-Precision-Balances-in-Vacuum.pdf?lang=en) **[contested]**
- (doc) Dresden instruments in the same SP2020+1 paper:
  - a torsion balance and an inverted counterbalanced double pendulum, both with "background noise lower than 5 nN after averaging";
  - calibration by voice coil (1–100 µN pulses) and laser interferometry (attocube IDS3010: 1 pm resolution, 2 nm noise);
  - devices in mu-metal shielding, swept from 20 to 50 kHz;
  - beam vibration recorded at up to 25 MHz sampling.
- (doc) Artefacts the same paper identified:
  - a Galinstan liquid-metal-contact repulsion of ≈1.5 µN drift plus ≈1 µN steady force;
  - a thermal drift of ≈30 nN;
  - switching transients just under 100 nN;
  - "significant vibration when exciting the actuator around its longitudinal resonance, regardless of the mounting and isolation".

  Conclusion: "the claimed thrusts using an MET or MEGA thruster device consist in vibrational artefacts ... Our tests reveal the presence of mechanical artefacts but no thrust." — [same paper](https://tu-dresden.de/ing/maschinenwesen/ilr/rfs/ressourcen/dateien/forschung/folder-2007-08-21-5231434330/ag_raumfahrtantriebe/SPC-Mach-Effect-Thruster-Experiments-on-High-Precision-Balances-in-Vacuum.pdf?lang=en) **[refuted]**
- (doc) The peer-reviewed version is Monette, Kößling & Tajmar, "Experimental investigation of Mach-Effect thrusters on torsion balances", *Acta Astronautica* 182 (May 2021) 531–546, doi:10.1016/j.actaastro.2021.02.037. Its highlights say the claim is "identified as a vibration artefact" and that "varying driving conditions do not significantly affect the forces observed, contrary to the theory". — [ADS](https://ui.adsabs.harvard.edu/abs/2021AcAau.182..531M/abstract); [Crossref metadata via search](https://www.sciencedirect.com/science/article/abs/pii/S0094576521001119)
- (doc) An earlier SpaceDrive balance (IAC-18, Tajmar, Kößling, Monette, Weikert) was designed to resolve ≥100 nN for thrusters up to 10 kg. It found "magnetic interaction from not sufficiently shielded cables or thrusters" to be a major source of false µN signals. — [IAF IAC-18 record](https://iafastro.directory/iac/archive/browse/IAC-18/C4/7-C3.5/46021/); [statnano/NBF summary](https://statnano.com/world-news/425/The-SpaceDrive-Project-First-Results-on-EMDrive-and-Mach-Effect-Thrusters-|-NextBigFuturecom)
- (doc, theory, preprint) Rodal (Space Studies Institute), arXiv:2512.22687 (Dec 2025): in covariant GR the candidate ∂²M/∂t² term is suppressed by ~(U/c²)(ωL/c)² ≪ 1. Both GR and Hoyle–Narlikar scalar–tensor gravity predict any Mach-effect thrust "is far too small for practical propulsion". — [arXiv:2512.22687](https://arxiv.org/abs/2512.22687) **[not yet peer-reviewed]**
- (doc) NASA BPP (2004) had earlier rated the Woodward approach "unresolved" (see section 2). — [NASA/TM-2004-213082](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20040070788.pdf)

### Inferences
- Computed check of eq. 10: the prefactor ω⁶m₀x₀³V₀⁴/(8πGρ₀c²) is 1.71×10²¹, against the paper's 1.87×10²¹. Using K_p and K_e as quoted gives F ≈ 7.0 nN. So Woodward's own theory predicts ≈7–8 nN. That sits right at the Dresden noise floor (<5 nN averaged) and ~10⁴× below the reported ~100 µN.
- The steep ω⁶V₀⁴ scaling makes vibration artefacts (which also scale steeply with drive amplitude at resonance) hard to separate. This is the Dresden diagnosis.
- Rust: `mach_thrust(omega, m0, kp, ke, x0, v0, rho0) = omega.powi(6)*m0*kp*kp*ke*x0.powi(3)*v0.powi(4)/(8.0*PI*G*rho0*C*C)`. Reproduce ≈7 nN as a unit test.

### Gaps
- I did not read Tajmar's model paper "Mach-Effect thruster model" (*Acta Astronautica* 141 (2017) 8–16) or Woodward's 2019–2021 "MEGA drive" rebuttals (>100 µN reported in the press), so I cannot give their formulas or numbers.
- The full text of the 2021 Acta paper (exact upper bound in nN per device) was not accessed. The <5 nN figure comes from the companion conference paper.

## 6. EmDrive: Shawyer's formula, Eagleworks' 1.2 mN/kW, and the Dresden bound below the photon limit

### Takeaway
Shawyer predicts F = 2PQD_f/c. Eagleworks (White et al. 2017) reported 1.2 ± 0.1 mN/kW in vacuum, ≈360× the photon-rocket limit. Tajmar, Neunzig & Weikert (2021) used a battery-powered EmDrive on an inverted counterbalanced double pendulum and saw **no thrust above a 10–30 nN noise floor at ≤11 W in the cavity**. They reproduced the White-like signal (≈1–20 µN) as a thermal-expansion artefact. Any anomalous thrust is at or below 1/c ≈ 3.3 nN/W. **[refuted]**

### Cited Findings
- (doc) **Shawyer's formula**, as given in the Dresden paper's eq. 1: F = 2·P·Q_u·D_f/c. Here P is the cavity input power and Q_u is the unloaded quality factor. D_f is a design factor that runs from 0.56 for a NASA-type geometry to 0.91 for optimised designs. Shawyer quotes Q_u of "several thousands" for first-generation designs and "several ten thousands" for second-generation ones. The paper notes the loaded Q is more realistic and that many reports do not say which Q they used. — [Tajmar, Neunzig & Weikert, "High-accuracy thrust measurements of the EMDrive and elimination of false-positive effects", *CEAS Space Journal* 14(1), 31–44 (2022; online July 2021), doi:10.1007/s12567-021-00385-1, open access (DNB copy)](https://d-nb.info/1244138614/34); Shawyer 2015, *Acta Astronautica* 116, 166–174, doi:10.1016/j.actaastro.2015.07.002 — [Wikipedia: EmDrive](https://en.wikipedia.org/wiki/EmDrive) **[refuted]**
- (doc) **Photon-rocket limit:** "F/P = 1/c = 3.3 nN/W". — [CEAS paper](https://d-nb.info/1244138614/34) **[mainstream]**
- (doc) **Eagleworks:** White, March, Lawrence, Vera, Sylvester, Brady & Bailey, "Measurement of Impulsive Thrust from a Closed Radio-Frequency Cavity in Vacuum", *J. Propulsion & Power* 33(4), 830–841 (2017), doi:10.2514/1.B36120. Torsion pendulum in vacuum, 40–80 W, thrust-to-power ratio "1.2±0.1 mN/kW". — [Semantic Scholar record](https://api.semanticscholar.org/graph/v1/paper/DOI:10.2514/1.B36120); [Wikipedia: EmDrive](https://en.wikipedia.org/wiki/EmDrive). The Dresden table summarises White's run as 1937 MHz, loaded Q 7123, 40–80 W, thrusts of 30–130 µN, "around an order of magnitude below Eq. (1)". — [CEAS paper](https://d-nb.info/1244138614/34). An earlier Eagleworks conference paper reported 91.2 µN at 17 W (Brady et al., AIAA 2014-4029). — [Wikipedia: EmDrive](https://en.wikipedia.org/wiki/EmDrive) **[refuted as artefact]**
- (doc) Other claims tabulated or reviewed — [CEAS paper](https://d-nb.info/1244138614/34); [Wikipedia: EmDrive](https://en.wikipedia.org/wiki/EmDrive):
  - Shawyer: 20 mN at 850 W (2002), never published in detail.
  - Yang (NWPU): 70–720 mN at 80–2500 W. A later pendulum test found nothing above a 0.7 mN uncertainty at 230 W.
  - Fetta: ~10 mN at 10.5 W with a superconducting cavity at Q ~10⁷, "poorly" documented.
  - Tajmar 2015: ~20 µN in vacuum, with an anomalous signal in the null orientation. A later torsion-balance test stayed within ±20 µN.
  - 2018: a 10,000× attenuator did not remove the "thrust"; the fault was traced to cable and magnetic shielding.
- (doc) **Dresden balance and noise:**
  - inverted counterbalanced double pendulum with voice-coil damping, ~15 s response, for thrusters up to 10 kg;
  - "sub-nN resolution could be obtained" with averaging;
  - final cavity tests: "All thrust values were below our noise level of 10–30 nN";
  - vacuum 10⁻² mbar (enough to suppress buoyancy) in a chamber 1.2 × 1.5 × 2.5 m;
  - battery-powered EmDrive on liquid-metal contacts with twisted-pair cabling;
  - photon benchmark at 11 W cavity power: 36 nN; 200 nN when counting the amplifier's 60 W electrical input.

  — [CEAS paper](https://d-nb.info/1244138614/34) **[mainstream measurement]**
- (doc) **The thermal artefact that reproduces White.**
  - Position A (cavity on a cage aligned with the platform): "a plateau and a drift similar to the ones reported in White et al. with 1.5 µN for an input power of 10 W". Screwing the cavity to the cage rods raised this to 20 µN. "This was not a real thrust but originated from mechanical stress due to the thermally expanding cavity, which acted on the balance's pivots."
  - Position B (cavity on dampers): ≈200 nN, from centre-of-gravity shifts.
  - Position C (hung from a bolt bearing): "the thrust artefacts disappeared and only noise was left."

  — [CEAS paper](https://d-nb.info/1244138614/34)
- (doc) **Buoyancy artefact in air:** "a typical cavity volume of 15 × 15 × 15 cm³, a temperature change of only one degree Celsius creates a lifting force of 150 µN!" — [CEAS paper](https://d-nb.info/1244138614/34)
- (doc) **Bounds.** "From Eq. (1), we would have expected around 1 mN ... around 5 orders of magnitude above our measurements. Following the data from White et al., we would have expected at least 10 µN, which is still three orders of magnitude above our data." The abstract says the data "limit any anomalous thrust to below the force equivalent from classical radiation for a given amount of power" and rule out previous results "by at least two orders of magnitude". The conference version and Wikipedia's quote say "at least 3 orders of magnitude". — [CEAS paper](https://d-nb.info/1244138614/34); [Semantic Scholar abstract](https://api.semanticscholar.org/graph/v1/paper/DOI:10.1007/s12567-021-00385-1); [Wikipedia: EmDrive](https://en.wikipedia.org/wiki/EmDrive)
- (doc) The Dresden group published three papers in 2021 claiming that every published positive result had been a false positive. — [Wikipedia: EmDrive](https://en.wikipedia.org/wiki/EmDrive)

### Inferences
- **Computed numbers.**
  - Photon limit: 1/c = 3.336 µN/kW. White's 1.2 mN/kW is ≈360× that. The Dresden text says "a factor of 500", which looks like rounding or an error.
  - Shawyer's eq. 1 at the Dresden operating point (P = 11 W, loaded Q = 741, D_f = 0.56) gives ≈30 µN, not the "around 1 mN" the paper states. Q ≈ 2×10⁴ would be needed for 1 mN. Their "~1 mN" probably assumes a Shawyer-class unloaded Q. Either way the measured ≤30 nN is ≥10³ below.
  - Shawyer's eq. 1 at White's conditions (60 W, Q = 7123) gives ≈1.6 mN, about 10× the 30–130 µN White reported.
- **Buoyancy check:** ΔF = ρ_air·V·g·ΔT/T. With V = 3.375×10⁻³ m³ and ΔT = 1 K at 293 K this gives ≈1.4×10⁻⁴ N, matching the paper's 150 µN.
- **Published thrust-measurement noise floors, for implementers:**

  | Experiment | Noise floor or accuracy | Power | Result |
  |---|---|---|---|
  | Tajmar 2004 (corona box) | 10.5 µN | 0–20 kV | null |
  | Dresden EmDrive 2021 | 10–30 nN (sub-nN averaged) | ≤11 W cavity | null; artefacts of 200 nN–20 µN identified |
  | Dresden MET 2021 | <5 nN averaged | ~tens of W, 20–50 kHz | artefacts only |
  | IAC-18 balance | ≥100 nN (design) | – | – |
  | Eagleworks 2017 | claimed 1.2 ± 0.1 mN/kW | 40–80 W | disputed (thermal) |
  | Yang | 0.7 mN uncertainty | 230 W | null |
  | Canning 2004 | rotation (qualitative) | 40–50 kV | null below 300 Torr |

- Rust: `shawyer(p, q, df) = 2.0*p*q*df/C`, `photon(p) = p/C`, `buoyancy(rho, vol, dT, T) = rho*vol*G0*dT/T`. Flag any measured F/P > 1/c as needing artefact checks: thermal centre-of-mass shift, pivot stress, cable-field magnetic coupling, buoyancy, outgassing.

### Gaps
- I did not obtain White et al.'s own balance noise figure or error budget. The JPP paper is paywalled.
- I did not read the third 2021 Dresden paper (on laser and superconducting EmDrives, CEAS 2022) for its bounds.

## 7. Secrecy and scrutiny on the public record: Invention Secrecy Act, Pais patents, classified gravity research, "born secret"

### Takeaway
Secrecy-order counts are public: 6,543 in force at the end of FY2025. They are not broken down by subject, and I found no public evidence of gravity- or propulsion-specific orders. The Pais "UFO patents" went the other way: the Navy pushed them through to publication, then a $508k NAWCAD test (2016–2019) could not show the "Pais Effect". The documented classified gravity work is geodetic and navigational: the Trident gravity gradiometer and the Geosat altimetry-derived marine gravity (declassified June 1995). "Born secret" applies only to the Atomic Energy Act's Restricted Data, and I found no documented link to gravity research.

### Cited Findings
**Invention Secrecy Act statistics (FAS)**
- (doc) The USPTO data compiled by Steven Aftergood — [FAS, Invention Secrecy Statistics FY05–FY25](https://sgp.fas.org/othergov/invention/stats.html):

  | | FY21 | FY22 | FY23 | FY24 | FY25 |
  |---|---|---|---|---|---|
  | Total orders in effect | 5,976 | 6,057 | 6,155 | 6,471 | 6,543 |
  | New orders | 61 | 87 | 125 | 356 | 102 |
  | Rescinded | 0 | 6 | 27 | 40 | 30 |

  Sponsors of new orders in FY24 were Navy 168, AF 113, DTSA 46, Army 12, Foreign 9, NASA/DARPA 5, NSA 2 and DOE 1. New DoD orders in FY24 split into Type 1 (export control) 145, Type 2 (classified) 94 and Type 3 102. "John Doe" orders on private inventors numbered 29, 1, 25, 0 and 18 for FY21–25. **[mainstream]**
- (doc) The screening guidance is the Patent Security Category Review List (PSCRL), supplied by the defence agencies, DOE, NASA and others and treated as confidential. FAS obtained a 1971 edition (declassified; 21 technology groups) and a 1994 edition under FOIA. One cited 1971 item is solar photovoltaics above 20% efficiency. — [FAS invention secrecy programme page](https://sgp.fas.org/othergov/invention/program.html); [FAS admin page](https://sgp.fas.org/othergov/invention/admin.html) (the PSCRL PDF at fas.org/sgp/othergov/invention/pscrl.pdf returned empty when I tried it)

**Pais Navy patents (2018–2019) and their review**
- (doc) The patents and applications — [Wikipedia: Salvatore Pais](https://en.wikipedia.org/wiki/Salvatore_Pais):
  - US 10,144,532 B2 "Craft using an inertial mass reduction device" (filed 2016, application 15/141,270; issued 4 Dec 2018; expired 9 Jan 2023 for non-payment);
  - US 10,135,366 B2 "Electromagnetic field generator" (issued 20 Nov 2018);
  - US 10,322,827 B2 "High frequency gravitational wave generator" (issued 18 Jun 2019);
  - US 2019/0058105 A1 "Piezoelectricity-induced room temperature superconductor" (published 21 Feb 2019);
  - US 2019/0295733 A1 "Plasma compression fusion device" (published 26 Sep 2019; abandoned).
- (doc) The USPTO first rejected the craft patent. It was granted after James Sheehy, Chief Technology Officer of the Naval Aviation Enterprise, wrote that it was operable and warned of similar Chinese work. — [Wikipedia: Salvatore Pais](https://en.wikipedia.org/wiki/Salvatore_Pais). Emails obtained by The War Zone show the inventions passed "technical and marketing reviews" in 2015 and the NAWCAD Invention Evaluation Board in 2016–2017. — [TWZ: "Emails Show Navy's 'UFO' Patents Went Through Significant Internal Review"](https://www.twz.com/37134/emails-show-navys-ufo-patents-went-through-significant-internal-review-resulted-in-a-demo)
- (doc) NAWCAD tested the "High Energy Electromagnetic Field Generator" from October 2016 to September 2019 for $508,000, of which $462,000 was salaries. The "Pais Effect" could not be proven, and no further research followed. Pais says the tests never reached the 1 coulomb of charge he considers necessary. — [TWZ](https://www.twz.com/37134/emails-show-navys-ufo-patents-went-through-significant-internal-review-resulted-in-a-demo); [Wikipedia: Salvatore Pais](https://en.wikipedia.org/wiki/Salvatore_Pais) **[refuted as tested]**

**Documented classified gravity research and data**
- (doc) Gravity gradiometry was developed by Bell Aerospace for the US Navy's Trident programme as a "stealth" technique for navigating ballistic-missile submarines. Bell, Anderson & Pratson describe it so in *The Leading Edge* (1997). — [Duke Scholars record](https://scholars.duke.edu/publication/794010). R. E. Bell calls it "formerly classified" in "Gravity Gradiometry", *Scientific American* 278(6):74 (June 1998). — [Scientific American](https://www.scientificamerican.com/article/gravity-gradiometry). It was declassified "in the early 1990s" and moved to oil and gas exploration. — [Bell Geospace company profile](https://www.rigzone.com/directory/company/3514/Bell-Geospace-Inc/) **[mainstream]**
- (doc) The accuracy of submarine-launched ballistic missiles depended on precise gravity knowledge along the trajectory. A gradiometer (Metzger, Bell) was made a requirement of the Trident inertial navigation system. — [Brazilian Journal of Geophysics review](https://www.sbgf.org.br/revista/index.php/rbgf/article/download/327/125)
- (doc) Geosat Geodetic Mission altimetry, the main source of high-resolution marine gravity before GRACE and GOCE, was classified by the US Navy. It was released in stages: 1990 (Antarctic waters), 1992 (south of 30°S), then globally in July 1995. — [Wikipedia: Geosat](https://en.wikipedia.org/wiki/Geosat); [GlobalSecurity](https://www.globalsecurity.org/space/systems/geosat.htm). Sandwell, Yale, McAdoo & Smith: "The high density coverage obtained by ERS-1 during its geodetic mapping phase (4/94–3/95) prompted the US Navy to declassify all of the Geosat altimeter data (6/22/95)." — [Sandwell et al. 1996 chapter (Springer)](https://www.springerprofessional.de/marine-gravity-from-satellite-altimetry-over-ocean-and-sea-ice/13386424). The release enabled the Scripps/NOAA v7.2 2-minute marine gravity grid. — [NOAA NCEI accession 0054150](https://www.ncei.noaa.gov/archive/accession/0054150) **[mainstream]**
- (doc) By contrast, USAF-funded GR research from 1956 to 1969/70 was, as far as the historian Goenner knows, unclassified. — [arXiv:1607.03319](https://arxiv.org/abs/1607.03319)

**"Born secret"**
- (doc) Restricted Data under the Atomic Energy Act covers the design, manufacture or use of atomic weapons and the production or use of special nuclear material. Because the definition covers "all data" not specifically declassified, it is read as making such information "born secret" even when privately created, "unlike any other forms of controlled national defense information". The test case was *United States v. Progressive* (1979), which was dismissed as moot. — [Wikipedia: Born secret](https://en.wikipedia.org/wiki/Born_secret); [Wikipedia: Restricted Data](https://en.wikipedia.org/wiki/Restricted_Data); [Wikipedia: U.S. v. Progressive](https://en.wikipedia.org/wiki/United_States_v._Progressive,_Inc.)

### Inferences
- **Born secret:** none of the sources I read connect "born secret" or Restricted Data to gravity research. The doctrine's statutory scope (nuclear weapons and special nuclear material) excludes it unless the work involves nuclear material. Gravity-related secrecy on the record runs through ordinary national-security classification (Navy geodesy and navigation) and Invention Secrecy Act orders, not the AEA.
- The FY24 spike to 356 new orders, with 94 Type 2 "classified" DoD orders, is real but has no subject breakdown. Linking it to propulsion or gravity would be speculation.
- The Pais episode is the opposite of secrecy: a service lab used its CTO to force publication of non-working concepts, and tested them only afterwards.

### Gaps
- I could not retrieve the PSCRL 1971 or 1994 text to check whether any category mentions gravity, gravitation, inertia or field propulsion.
- I found no FOIA-released list of individual secrecy orders with gravity or propulsion titles.
- I found no primary document showing DMA or NIMA classification of geopotential models (e.g. pre-EGM96 DoD models) or their declassification dates. This is plausible but unsourced here.

## 8. Gravity gradiometry's classified history: Bell's rotating-accelerometer GGI, its declassification, and BHP's Falcon (Eötvös performance)

### Takeaway
Bell Aerospace (Niagara Falls; later Textron, then Lockheed Martin) built the rotating-accelerometer gravity gradiometer instrument (GGI): four accelerometers on a spinning disc, with the gradient signal at twice the spin frequency. It was built for Navy SSBN inertial-navigation gravity compensation, declassified in the early 1990s, and became BHP's Falcon AGG (first survey October 1999) and Bell Geospace's FTG and Air-FTG. Falcon's noise fell from ~30–40 E/√Hz in early trials to ~2.5–3.5 E RMS in survey data by 2005–2007 (1 E = 10⁻⁹ s⁻²).

### Cited Findings
- (doc) History, from Dransfield (BHP Billiton), "Airborne Gravity Gradiometry in the Search for Mineral Deposits", *Proceedings of Exploration 07*, Paper 20 — [911metallurgist PDF](https://www.911metallurgist.com/wp-content/uploads/2015/10/Airborne-Gravity-Gradiometry-in-The-Search-for-Mineral-Deposits.pdf):
  - From 1970s prototypes (Forward; Trageser; Metzger) "the US Navy selected the Bell Aerospace Gravity Sensors System (GSS) for gravity compensation of its inertial navigation systems."
  - In 1983 the Air Force Geophysics Laboratory chose the same GGI for a Defense Mapping Agency airborne Gravity Gradiometer Survey System. Its first airborne test (1988) had GGI performance problems, and DMA never adopted it.
  - In the early 1990s Lockheed Martin built a lower-noise GGI for arms-control verification.
  - BHP and Lockheed Martin signed an agreement in 1996. The first airborne gravity-gradient survey (Falcon) flew over the Bathurst Camp, New Brunswick, in October 1999.
  - Bell Geospace took a ship-borne FTG in 1998 and flew its first Air-FTG in early 2003.
- (doc) Bell Aerospace's 1976 NASA study treated the rotating-accelerometer gradiometer as the basis for an orbital instrument. — [NTRS 19770007712](https://ntrs.nasa.gov/citations/19770007712). Falcon "grew out of" technology "developed for the US Navy in the 1970s by Bell Aerospace (now Lockheed Martin)" that "became available for commercial application at the end of the 'cold war'". — [BHP 2001 Falcon abstract](https://www.bhp.com/-/media/bhp/documents/investors/reports/2001/mining2001abstractfalcon8november0.pdf). Early Gulf of Mexico surveys used the USNS Vanguard's installed FTG. — [SPE / JPT summary](https://jpt.spe.org/gravity-oil-and-gas-exploration)
- (doc) **How the GGI works.**
  - Each "complement" has four accelerometers equi-spaced on a circle, sensitive axes tangential.
  - It "rejects both common mode acceleration and rotations about the axis perpendicular to the plane".
  - "Rotation of the complement about the perpendicular axis moves the gradient signal to twice the rotation frequency".
  - The FTG has three GGIs on mutually orthogonal spin axes, each with one complement on a ~15 cm diameter circle. The AGG (Falcon) has one GGI with a larger complement, giving "intrinsic noise power ... eight times lower" per GGI and "two and two-thirds smaller" overall than the FTG.
  - The AGG's near-horizontal disc cuts coupling to vertical turbulence. Typical vertical aircraft acceleration has twice the power of horizontal.

  — [Dransfield 2007](https://www.911metallurgist.com/wp-content/uploads/2015/10/Airborne-Gravity-Gradiometry-in-The-Search-for-Mineral-Deposits.pdf)
- (doc) **Requirements.**
  - "For a gravity gradiometer with a 10 cm accelerometer separation to achieve a useful noise limit of 10 Eö/√Hz in a typical low-level survey acceleration environment of 1 m s⁻²/√Hz requires a matching of one part per billion (Lee, 2001)."
  - Rotational-acceleration cancellation needs accelerometer pairs matched to ~10 ppb.
  - Pseudo-gradients from rotation rates need three-axis rotational control at the ~10 µrad level.
  - A navigation and terrain error below 1 m keeps the gradient error to ~1.8 E at 80 m height.

  — [Dransfield 2007](https://www.911metallurgist.com/wp-content/uploads/2015/10/Airborne-Gravity-Gradiometry-in-The-Search-for-Mineral-Deposits.pdf)
- (doc) **Performance in Eötvös.**
  - Early Falcon trials: ~30–40 E/√Hz, about an order of magnitude above what low-clearance surveys need. — [Lee, "FALCON gravity gradiometer technology", ASEG Extended Abstracts 2001](https://www.tandfonline.com/doi/abs/10.1071/ASEG2001ab068) (from a search summary of the abstract)
  - Noise densities quoted by Dransfield: Falcon 6 Eö/√Hz, Air-FTG 13 Eö/√Hz.
  - Falcon survey noise: 2.5 Eö RMS filtered to 300 m wavelength at 55 m/s (Boggs et al. 2007). The 2005 average was 3.5 Eö RMS at 300 m.
  - Air-FTG: from 15+ to 5.4 to 3.5 Eö RMS at 800 m wavelength over 3 years (Murphy et al. 2007).
  - Falcon's bandwidth gives "spatial resolution nearly 3 times better than the Air-FTG".

  — [Dransfield 2007](https://www.911metallurgist.com/wp-content/uploads/2015/10/Airborne-Gravity-Gradiometry-in-The-Search-for-Mineral-Deposits.pdf). Lee's ASEG 2001 paper also defines the Eötvös as 10⁻⁹ s⁻² and gives noise in E²/Hz. — [CSIRO ASEG2001ab019](https://publish.csiro.au/ex/ASEG2001ab019) **[mainstream]**
- (doc) A Dransfield & Lee paper gives 17 Eö/√Hz for the GDD noise density. — [search summary citing ResearchGate "The FALCON airborne gravity gradiometer survey systems"](https://www.researchgate.net/publication/327069609_The_FALCON_airborne_gravity_gradiometer_survey_systems)

### Inferences
- Computed: Γ_noise·L/a_env = (10×10⁻⁹ s⁻²)(0.1 m)/(1 m/s²) = 10⁻⁹. This reproduces Lee's 1 ppb matching requirement. Write it as `matching = gamma_noise * baseline / accel_env`.
- Rust model of a rotating GGI. Four tangential accelerometers at radius R, spin rate Ω. The summed output is modulated at 2Ω, and its demodulated amplitude is proportional to (Γ_xx − Γ_yy) and Γ_xy. Converting a noise density (E/√Hz) to RMS needs the survey filter: σ ≈ ASD·√(bandwidth), where the bandwidth ≈ v/λ for ground speed v and filter wavelength λ.
  - Computed: 6 E/√Hz × √(55/300 Hz) ≈ 2.6 E RMS, consistent with the 2.5 E RMS Falcon figure.
- Link to this repo: the Navy-to-mining lineage of the GGI is the one well-documented "classified gravity technology" story. Satellite gradiometry (GOCE-class) differs from airborne AGG in noise, altitude and bandwidth, but I did not source satellite noise figures in this note.

### Gaps
- I found no declassified Navy document giving the GSS/GGI's Trident-era noise in Eötvös, or the exact declassification date and authority (only "early 1990s" and "end of the cold war").
- The full texts of Lee (2001) and Boggs et al. (2007) were not retrieved. The 30–40 E/√Hz early figure comes from a search summary of the abstract.
- I did not find Air-FTG's own noise density in E/√Hz from a Bell Geospace primary source.
