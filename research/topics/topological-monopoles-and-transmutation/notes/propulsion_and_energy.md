# Propulsion and energy claims in Greenyer, "The Power and Signature of the Topological Monopole" (TeslaTech, 15 Aug 2026)

Scope: the ARV craft and its simulator, Holt 1979, Forward 1963, the MHD and British Rail patents,
"oxygen as the fuel", relic neutrinos, vacuum/Casimir energy, high-frequency gravitational waves,
and the user's own audit simulator. Every number marked *computed* is reproduced and asserted in
`propulsion_and_energy_checks.py` (89 checks, standard library, < 0.1 s, `python3 -I`, exit 1 on
failure). Downloads are in `notes/dl_propulsion/` (git-ignored). Transcript line numbers refer to
`/home/dan/Desktop/transcript_1793004282390.txt` (auto-captions, so wording is approximate).

What the talk claims (transcript): the "alien reproduction vehicle" "could do up to eight times the
speed of light" with large oxygen tanks outside the cabin (l. 48–87); "The oxygen is the fuel …
There's about 1 ton of it" (l. 1820–1831); the simulator has 48 capacitor plates, 2.5 MV, 30–150 kHz
"as it was specified on Mark McCandlish's specifications", 250 ns current pulses, "7.3 teraamps"
electron feed, 950 kg of oxygen, "28 days of fuel of permanent acceleration", and "seems to get in
front of its own light cone" (l. 1832–2240); "1 kg will give you 30 33 days … like 900 days" with
coasting (l. 2525–2535); "938 million electron volts per proton … 166 times more energy than
deuterium deuterium fusion" (l. 2514–2522); relic neutrinos as a "superfluid condensate … four to
five times the mass of all ordinary matter" (l. 1199–1217) and waterspouts "running with the flux of
relic neutrinos" (l. 1230–1240); energy emitted as "light, electrons, coherent matter waves, and
high-frequency gravitational waves" (l. 297–305); vacuum energy/Casimir with "transmission of large
power over thin wires, shielding and propulsion" (l. 1131–1153).

## Holt 1979 and Forward 1963: what the cited sources actually say

### Takeaway
Holt's NASA TM-80961 is a 13-page speculative concept paper. It contains no equations and no design
numbers: no plate count, voltages, frequencies or oxygen. Its mechanism is a "resonance" between
pulsed electromagnetic fields and "space-time metrics", which teleports the craft to a "distant
space-time point". This is *unverified* and has no quantitative content to test. Forward 1963 is
mainstream general relativity (gravitomagnetic and induced accelerations, suppressed by G/c²), and
the repository already implements it.

### Cited Findings
- NTRS record 19800010907: A. C. Holt, "Field resonance propulsion concept", NASA-TM-80961 / JSC-16073, August 1979 (NTIS N80-19184). The abstract proposes "a resonance between coherent, pulsed electromagnetic wave forms, and gravitational wave forms (or space-time metrics)", using "magnetic field line merging, hydromagnetic wave effects, free-electron lasers, laser generation of megagauss fields" — [NTRS](https://ntrs.nasa.gov/citations/19800010907)
- From the full text (PDF opened): the foreword says the concept was presented at a special session of the 15th Joint AIAA/SAE/ASME Propulsion Conference (June 1979). It says the concept "was developed as the result of private, unofficial research. NASA is not involved in UFO research." The paper rests on "two assumptions" (spacetime is a "projection" of a higher-dimensional space, and an EM–gravity relationship exists). It cites UFO reports (e.g. "calculated speeds of 22,000 mph or more", "absence of sonic booms"). In the proposed implementation, high-powered lasers generate oppositely directed megagauss fields whose merging point is made to oscillate by "alternately pulsing adjacent laser sets". It recommends UFO radiation measurements and "Regressive hypnosis". It ends: "the concept utilizes two unproven assumptions and currently unverified magnetohydrodynamic processes" — [NTRS PDF](https://ntrs.nasa.gov/api/citations/19800010907/downloads/19800010907.pdf)
- The text contains no equations, and its only numbers are 2–3 kG sunspot fields, "10% of the sun's total output" for some flares, and 22,000 mph. The figures (Fig. 3, "Field resonance system"; Fig. 4, "winged disc") cannot be read from the text layer — [NTRS PDF](https://ntrs.nasa.gov/api/citations/19800010907/downloads/19800010907.pdf)
- Later papers by Holt: "Hydromagnetics and future propulsion systems", AIAA Student Journal (NTRS 19800052277). It mentions "superconducting magnets and a configuration of tunable free-electron lasers" — [NTRS](https://ntrs.nasa.gov/citations/19800052277). Also "Prospects for a breakthrough in field-dependent 'propulsion'", AIAA Paper 80-1233 (NTRS 19800054809), which proposes "gravimagnetic systems and field resonance systems" — [NTRS](https://ntrs.nasa.gov/citations/19800054809)
- Forward, "Guidelines to Antigravity", Am. J. Phys. 31, 166–170 (1963), DOI 10.1119/1.1969340 (Crossref-checked). The abstract says the paper treats "little known aspects of Einstein's general theory of relativity", and "Three distinctly different non-Newtonian gravitational forces are described" — [AIP](https://pubs.aip.org/aapt/ajp/article/31/3/166/1046748/Guidelines-to-Antigravity)
- The repository already has `gravitics::spacetime::induced_acceleration(m, a_flow, r) = 4Gm a_f/(c² r)`, citing Forward 1963. Its doc comment notes that one g needs m/r ≈ 3×10²¹ kg/m (crates/gravitics/src/spacetime.rs, lines 39–45).

### Inferences
- Holt's paper does not describe 48 plates, 2.5 MV, 30–150 kHz, oxygen fuel or "six yin-yang pairs". The talk attributes the six/twelve "phasing regions" to a later Holt paper the speaker "does have access to" (l. 1866–1872), which I could not identify. The 1979 TM only says the concept "is currently undergoing revision". Holt's actual hardware is lasers, megagauss fields and field-line merging. Reconnection is mainstream plasma physics (`plasmoid::relaxation` and `magnetics::mhd` cover it), but nothing in it gives propulsion without reaction mass.
- Label: Holt field resonance: *unverified* (no testable quantity); its reconnection ingredient: *mainstream*. Forward 1963: *mainstream* (already in the crate). The talk's "Robert Ford 1962" is Robert L. Forward, Am. J. Phys. 1963.

### Gaps
- I did not open the AJP full text of Forward 1963 (PDF only, behind a paywall). Its numerical examples are not quoted here; the repository's formula is taken from its existing citation.
- I did not find the full text of AIAA 80-1233, so the "12 phase sections" attribution is unverified.

## The 1967 vectored MHD patent and the British Rail "space vehicle" patent

### Takeaway
US 3,322,374 (King, 1967) is an air-breathing induction-MHD craft: it ionizes air and pushes it with a travelling magnetic field. That is physically sound in principle (reaction mass = air), but it gives no numbers. GB 1310990 (British Railways Board, 1973) is a laser-fusion craft that directs charged fusion products with magnets. It is a rocket whose exhaust is the fusion products, also without numbers. Neither patent relies on oxygen fuel or on propellantless thrust.

### Cited Findings
- US3322374A, "Magnetohydrodynamic propulsion apparatus", James F. King Jr. Filed 30 Sep 1964, granted 30 May 1967. An ionizer at the nose makes the air conductive, and stacked superconducting driving rings excited by polyphase AC create a travelling field that induces eddy currents in the ionized air, which is pushed opposite to the travel direction. Vectoring comes from sectored ionizer arcs (sector control), field shaping, and shorted skirt coils. There are 9 claims, and "no numbers for voltages, field strengths, frequencies, currents, or thrust" — [Google Patents](https://patents.google.com/patent/US3322374A/en)
- GB1310990A, "Space vehicle", applicant British Railways Board. Priority/filing 11 Dec 1970, published 21 Mar 1973, lapsed 1976 for non-payment of renewal fees (legal events). In the abstract, a fusion zone below the platform is fed liquid fuel ignited by lasers; electromagnets (possibly superconducting) deflect the charged fusion products; some products land on insulated electrodes to generate power; a homopolar generator powers the lasers; electrode voltages and magnet fields steer the thrust. The Google Patents page names no inventor and gives no numbers — [Google Patents](https://patents.google.com/patent/GB1310990A/en)
- FlightGlobal (2 Sep 1996) quotes a BR Board spokesman calling it "a patent for a flying saucer" and says the inventor is "now retired" (unnamed), filed 1973. Conflict: it implies the patent was still live in 1996, while the Google Patents legal events show a 1976 lapse — [FlightGlobal](https://www.flightglobal.com/spacecraft-on-the-right-lines-for-br/1517.article)

### Inferences
- The talk says the BR craft "has a nuclear reactor on a board to provide the … ions" (l. 1812–1815). The patent abstract says laser fusion, and the charged products are both the power source and the exhaust. That makes it a fusion rocket, bounded by the exhaust-velocity accounting below (fusion products have u ≈ 0.03–0.1 c).
- The King patent is air-breathing: its thrust per power is bounded by actuator-disk momentum theory (computed below: ideal ~1.1 MW to hover 5 t on a 3.5 m-radius disk), and it gives no thrust in vacuum. Status of both mechanisms: *mainstream physics, unbuilt devices*.

### Gaps
- The inventor's name "Charles Osmond Frederick" for GB1310990 appears only in search-result snippets of press articles. Neither the Google Patents page nor the FlightGlobal article I opened names him. I did not open the Espacenet full-text specification or claims.
  - **Resolved (2026-10-26, coordinator):** the printed specification, which the user supplied, names him. Its front page reads "(72) Inventor CHARLES OSMOND FREDERICK", "(54) SPACE VEHICLE", application No. 59083/70 filed Dec. 1970, complete specification published 21 March 1973, classification G21B 1/00; B64G 1/00, applicant British Railways Board. The text confirms the design read above: pulsed laser fusion above 1000 Hz, liquid fuel pumped through a nozzle, and superconducting magnets that deflect the charged fusion products onto insulated electrodes.

## "Oxygen as the fuel": energy per kilogram, both ways

### Takeaway
Per kilogram, oxygen with hydrogen gives 1.5×10⁷ J/kg O (chemistry, mainstream). ¹⁶O+¹⁶O fusion gives up to 16.54 MeV per reaction, 5.0×10¹³ J/kg (mainstream nuclear physics, but it needs stellar carbon/oxygen-burning temperatures). Complete conversion gives c² = 9.0×10¹⁶ J/kg (no known mechanism; proton decay is unobserved). The talk's "millions of electron volts per fusion of two oxygens" is correct. The "938 MeV per proton … 16 times if it's an oxygen" overstates the full mass–energy of ¹⁶O by 0.76%. The "166 times more energy than D–D" cannot be reproduced: full conversion is 257× D–D per reaction and ~1028× per unit mass.

### Cited Findings
- AME2020 atomic mass excesses (keV): n 8071.31806; ¹H 7288.971064; ²H 13135.722895; ³H 14949.81090; ³He 14931.21888; ⁴He 2424.91587; ¹²C 0; ¹⁶O −4737.00217; ¹⁷O −808.76421; ²²Ne −8024.716; ²⁸Si −21492.79711; ³¹P −24440.54442; ³¹S −19042.531; ³²S −26015.53714; ⁵⁶Fe −60607.163 — [AME2020 mass_1.mas20, IAEA AMDC](https://www-nds.iaea.org/amdc/ame2020/mass_1.mas20.txt)
- ΔfH°(H₂O, gas) = −241.826 ± 0.040 kJ/mol (CODATA review, Cox, Wagman et al.) — [NIST WebBook](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Mask=1)
- *Computed* Q values (atomic masses; electrons balance): ¹⁶O+¹⁶O → ³²S 16.5415 MeV; → ²⁸Si+α 9.5939; → ³¹P+p 7.6776; → ³¹S+n 1.4972; D+D → ³He+n 3.2689; D+D → T+p 4.0327. The user's audit values: 4 ¹⁷O → ⁵⁶Fe+¹²C is 57.372 MeV (audit says 49.47); 2 ¹⁷O → ¹²C+²²Ne is 6.407 MeV (audit says 17.13); ¹²C+¹⁶O → ²⁸Si 16.756 MeV (audit 16.76, correct) — from the AME2020 values above.
- *Computed* specific energies: H₂/O₂ → H₂O(g) 1.342×10⁷ J/kg of propellant (1.511×10⁷ J/kg of oxygen), mass fraction ε = 1.49×10⁻¹⁰. ¹⁶O+¹⁶O → ³²S 4.99×10¹³ J/kg, ε = 5.55×10⁻⁴. 4 ¹⁷O → Fe+C 8.14×10¹³ J/kg of ¹⁷O. Complete conversion c² = 8.988×10¹⁶ J/kg. The audit's "16 × m_p c² per ¹⁶O" is 9.053×10¹⁶ J/kg, 1.0076 × E = mc² of a ¹⁶O atom.
- The talk ties oxygen's role to its magnetic susceptibility (l. 1720–1735). The audit quotes χ_m ≈ +3449×10⁻⁶ cm³/mol for O₂ — [user audit, unverified](../../../external/2026-10-26-greenyer-audit-and-crate.md)

### Inferences
- Paramagnetism does not change nuclear Q values. Q is fixed by masses, which AME2020 gives to better than 1 keV for every species here. The branching between the ¹⁶O+¹⁶O channels and the cross sections at accessible energies are the other researcher's scope.
- The audit's Q values for the ¹⁷O pathways are wrong in opposite directions (−14% and +167%). The corrected values are used in the audit section below.
- Status: O₂ chemistry and ¹⁶O+¹⁶O Q values *mainstream*. O+O fusion as a craft fuel at room conditions *unverified* (no measured rate). Baryon collapse / complete conversion of ordinary oxygen: *unverified*. Monopole catalysis of proton decay (Callan–Rubakov) is a mainstream theoretical prediction, but no GUT monopole has been observed, and the talk's "topological monopole" is not a GUT monopole (other researchers' scope).

### Gaps
- I did not reconstruct "166×" from any combination of inputs. Candidates (m_p c²/⟨Q_DD⟩ = 257; per-mass ratio ~1028; one proton per ¹⁶O nucleus vs D–D per unit mass ≈ 64) all differ.

## Rocket accounting for 950 kg of oxygen

### Takeaway
Take a 5000 kg craft (the audit's value) carrying 950 kg of oxygen (the simulator's value). Even complete conversion into a perfectly collimated photon exhaust gives only Δv = 0.208 c. At 1 g that lasts 74.6 days of proper time. O+O fusion exhaust gives Δv = 2.1×10⁶ m/s (0.7% c), and H₂/O₂ chemistry gives 1.2 km/s. "Eight times light speed" is impossible under any of these. "28 days of permanent acceleration" from 950 kg would need a photon drive of 35 TW (2.7 g). "1 kg gives 30–33 days" is 33 GW, worth only 110 N of photon thrust (2 milli-g on 5 t).

### Cited Findings
- Photon-rocket thrust per power is 1/c = 3.336 µN/W, already in `gravitics::propulsion::photon_thrust` (crates/gravitics/src/propulsion.rs). The repository has no rocket equation (grep for "rocket", "Tsiolkovsky", "exhaust" across crates/ found none).
- Relativistic rocket (Ackeret 1946; the standard form is in the Acta Astronautica discussion of the photon-rocket limit) — [Haug 2017, DOI 10.1016/j.actaastro.2017.03.011](https://doi.org/10.1016/j.actaastro.2017.03.011) (Crossref-checked title only; formula derived below, not quoted)
- *Computed* (formulas derived in the check script from energy–momentum conservation in the instantaneous rest frame): Δv = c·tanh((u/c)·ln R), with R = m₀/m₁ = 5000/4050 and u = c·√(ε(2−ε)) the exhaust speed when a fraction ε of the propellant rest mass becomes exhaust kinetic energy. For a photon exhaust (ε = 1) this reduces exactly to (R²−1)/(R²+1).
  - Photon rocket: Δv = 0.2077 c. Proper time at 1 g: τ = (c/g) ln R = 74.6 days.
  - ¹⁶O+¹⁶O → ³²S exhaust: u = 9.99×10⁶ m/s, Δv = 2.10×10⁶ m/s.
  - H₂/O₂ (950 kg O + 120 kg H₂): u = 5.18 km/s ideal, Δv = 1.21 km/s.
  - 950 kg over 28 days as a photon rocket: P = 3.53×10¹³ W, constant proper acceleration 2.66 g.
  - 1 kg fully converted over 31.5 days: P = 3.30×10¹⁰ W, photon thrust 110 N.

### Inferences
- The talk's two durations are mutually inconsistent at a fixed burn: 950 kg × 30–33 days/kg is 78–86 years of burn, not 28 days. They can be reconciled only if "28 days" is a higher burn setting (the talk mentions a "burn rate" control, l. 1913).
- Even with energy unlimited, a craft that carries its reaction mass cannot exceed c. In the air it can borrow reaction mass (the King patent; ion wind), but then its thrust stops in vacuum.
- Status: rocket bounds *mainstream*. ARV ">c": *refuted* as a propulsion claim within special relativity. A warp metric is a different question, covered by `gravitics::spacetime::alcubierre_energy` (negative energy of order solar masses).

### Gaps
- The craft's real mass is not given in the talk. The 5000 kg is the audit's assumption; the results scale through R.

## Thrust-to-power limits: photon, air-breathing EHD/MHD, ion wind in air versus vacuum

### Takeaway
In vacuum, nothing that only radiates can beat F/P = 1/c = 3.3 mN/kW. In air, ion wind (EHD) reaches tens of N/kW because air is the reaction mass. Masuyama's MIT thesis measured 68.43 mN/W, and MIT flew a 5 m-span EHD aeroplane at 40 kV. The repository's `ion_wind_thrust` (F = I·d/μ) already shows the thrust vanishing as mobility → ∞ (vacuum).

### Cited Findings
- Masuyama, MIT S.M. thesis (2012), "Performance characterization of electrohydrodynamic propulsion devices": "A thrust-to-power ratio as high as 68.43 mN/W was obtained for a d = 21 cm thruster", and a longer gap "generates more thrust per input power" — [DSpace@MIT 1721.1/77114](https://dspace.mit.edu/handle/1721.1/77114)
- Masuyama & Barrett, "On the performance of electrohydrodynamic propulsion", Proc. R. Soc. A 469, 20120623 (2013), DOI 10.1098/rspa.2012.0623 (Crossref-checked; publisher page returned 403, so no number is quoted from it). Press accounts cite 110 N/kW; I did not verify this against the paper.
- Xu, He, Strobel, Gilmore et al., "Flight of an aeroplane with solid-state propulsion", Nature 563, 532–535 (2018), DOI 10.1038/s41586-018-0707-9 (Crossref-checked). The abstract reports a fixed-wing aeroplane with "a five-metre wingspan" flown "ten times" in "steady-level flight", carrying all batteries and an "ultralight high-voltage (40-kilovolt) power converter" — [Europe PMC abstract](https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:%2210.1038/s41586-018-0707-9%22&resultType=core&format=json)
- Repository: `ion_wind_thrust(I, d, μ) = I·d/μ`, `ballistic_thrust`, `photon_thrust`. The tests note that "Without air, the same ions at 30 kV give 1000× less" (crates/gravitics/src/propulsion.rs, tests).
- *Computed*: one-dimensional EHD theory gives F/P = d/(μV) from F = Id/μ and P = IV. The ideal actuator disk gives P = F^{3/2}/√(2ρA): hovering 5000 kg on a 3.5 m-radius disk in sea-level air needs 1.12 MW.

### Inferences
- The talk's "unlimited" fuel in air (l. 1920) is consistent only with air-breathing operation, where oxygen (or air) is reaction mass. The thrust-per-power limits then are the EHD/actuator-disk ones, not anything specific to oxygen.
- Status: EHD in air *mainstream (measured)*. EHD/"lifter" thrust in vacuum *refuted* (it is ion wind). Biefeld–Brown as a gravitational effect is *refuted* (the repository already covers this).

### Gaps
- The thrust-to-power figure from the 2018 Nature paper is not in its abstract, and I could not open the full text.

## The simulator's electrical numbers: 48 plates at 2.5 MV, 30–150 kHz, 250 ns, 7.3 TA

### Takeaway
2.5 MV across plates big enough to fit the craft stores only about a kilojoule and gives about a kilonewton of electrostatic attraction between the plates (an internal force, so zero net thrust). A 7.3 TA current is not physical at this scale. At 3.5 m its self-field is 4.2×10⁵ T (348× the indoor record of 1200 T), with a magnetic pressure of 6.9×10¹⁶ Pa. A ring carrying it stores 2.6×10²⁰ J, 3.1 times the complete mass–energy of all 950 kg of oxygen. An electron beam carrying it would need βγ ≈ 4×10⁸ to stay under the Alfvén limit.

### Cited Findings
- Talk: 2.5 MV, 30–150 kHz, 250 ns pulses, 7.3 TA "even at normal feed rate … So, it's just absurd" (transcript l. 1960–1995).
- Record indoor field of 1200 T by electromagnetic flux compression: Nakamura, Ikeda, Sawabe et al., Rev. Sci. Instrum. 89, 095106 (2018), DOI 10.1063/1.5044557 (Crossref-checked; the number is from the title).
- Alfvén, "On the motion of cosmic rays in interstellar space", Phys. Rev. 55, 425 (1939), DOI 10.1103/PhysRev.55.425 (Crossref-checked). I_A = 4πε₀mc³βγ/e = 17.0 kA·βγ is *computed* in the script.
- *Computed*, with an assumed geometry: an annulus 0.5–3.5 m (37.7 m², i.e. 48 wedges of 0.785 m²) and a 1 m gap, chosen so that 2.5 MV/m stays below air breakdown (~3 MV/m, not sourced here). Results: C = 334 pF, U = ½CV² = 1.04 kJ, pressure ε₀E²/2 = 27.7 Pa, force 1.04 kN, displacement current CVω = 472 A at 90 kHz (reactive).
- *Computed*: B = μ₀I/(2πr) = 4.17×10⁵ T at 3.5 m. Pressure B²/2μ₀ = 6.92×10¹⁶ Pa. Ring inductance μ₀R(ln(8R/a) − 2) = 9.89 µH (R = 3.5 m, a = 0.4 m). Energy ½LI² = 2.64×10²⁰ J = 3.09 × (950 kg)c². The electron mass flow I·m_e/e is 41.5 kg/s.

### Inferences
- Forces between capacitor plates are internal stresses and cancel in the craft's total momentum. Only momentum carried away by fields (≤ P/c) or matter can propel it.
- The simulator's "7.3 TA" and "2.5 MV" come from a display, not a circuit model. Status: as a physical current, *refuted* by energy accounting (it exceeds the craft's whole mass–energy budget).

### Gaps
- McCandlish's ARV specifications (diameter, plate area, gap) were not located in a primary document, so the geometry above is an assumed test case.

## Relic neutrinos: "superfluid condensate … four to five times all ordinary matter"

### Takeaway
The cosmic neutrino background is real (indirectly detected). It has T_ν = 1.95 K and n ≈ 336–340 cm⁻³ in total. Its mass density is Ω_ν h² = Σm_ν/93.14 eV < 6.9×10⁻⁴ (DESI DR2, Σm_ν < 0.064 eV), which is under 3.1% of the baryons, not 400–500%. The "four to five times ordinary matter" figure matches cold dark matter (Ω_c/Ω_b = 5.36), not neutrinos. Relic neutrinos are fermions and cannot Bose-condense. Kapusta's Higgs-exchange BCS pairing gives T_c "less than 1 K by many orders of magnitude" at cosmological densities. Their largest calculated mechanical effect is ~10⁻²⁶ cm/s², 6×10²⁷ times smaller than g. Status: the neutrino background is *mainstream*; "four to five times" is *refuted*; a "condensate" is *refuted* for the known neutrinos; driving waterspouts is *refuted*.

### Cited Findings
- PDG 2024, "Neutrinos in Cosmology" (Lesgourgues & Verde): "Tν0 ≃ 1.7 × 10⁻⁴ eV ≃ 1.9 K"; Ων = Σmν/(93.12 h² eV) (Eq. 26.2); "total neutrino average number density today: n0ν = 339.5 cm⁻³". Neutrinos "cannot be all the dark matter". Direct detection (PTOLEMY) depends on local clustering, "an effect however expected to be small" — [PDG 2024 rev. 26](https://pdg.lbl.gov/2024/reviews/rpp2024-rev-neutrinos-in-cosmology.pdf)
- DESI DR2 Results II (arXiv:2503.14738): Ων h² = Σmν/93.14 eV (Eq. 7); "95% upper limits on the sum of neutrino masses, finding Σmν < 0.064 eV assuming ΛCDM and Σmν < 0.16 eV in the w0wa model"; the oscillation lower limit is "Σmν > 0.059 eV" — [arXiv:2503.14738](https://arxiv.org/abs/2503.14738)
- Planck 2018 VI (A&A 641, A6, 2020, DOI 10.1051/0004-6361/201833910, Crossref-checked): Ωc h² = 0.120 ± 0.001, Ωb h² = 0.0224 ± 0.0001, Σmν < 0.12 eV; T₀ = 2.7255 K (Fixsen 2009) — [arXiv:1807.06209](https://arxiv.org/abs/1807.06209). Fixsen, ApJ 707, 916 (2009), DOI 10.1088/0004-637X/707/2/916 (Crossref-checked).
- KATRIN: "m_ν < 0.45 eV at 90 % confidence level" from 259 days of data — [arXiv:2406.13516](https://arxiv.org/abs/2406.13516). The local relic overdensity limit is "9.7e10 (1.1e11) at a 90% (95%) confidence level" — [arXiv:2202.04587](https://arxiv.org/abs/2202.04587)
- Kapusta, "Neutrino Superfluidity", PRL 93, 251801 (2004), DOI 10.1103/PhysRevLett.93.251801 (Crossref-checked). Dirac neutrinos "display BCS superfluidity for any nonzero mass" via Higgs exchange. For "a neutrino of mass 1 eV" the exponent is "about −10⁴⁶", and at cosmological densities "the critical temperature Tc is less than 1 K by many orders of magnitude no matter what the neutrino mass is" — [arXiv:hep-th/0407164](https://arxiv.org/abs/hep-th/0407164)
- Duda, Gelmini & Nussinov, "Expected signals in relic neutrino detectors", PRD 64, 122001 (2001), DOI 10.1103/PhysRevD.64.122001 (Crossref-checked). The G_F effects give accelerations "of the order of f · 10⁻²⁷ cm/sec²" with f "at most about 17", at most "f · 10⁻²⁶ cm/sec²" for non-clustered neutrinos, "at most thirteen orders of magnitude smaller" than the smallest measurable 10⁻¹² cm/s². The G_F² effects give f·10⁻³⁴ cm/s² (relativistic) — [arXiv:hep-ph/0107027](https://arxiv.org/abs/hep-ph/0107027)
- Tremaine & Gunn, "Dynamical role of light neutral leptons in cosmology", PRL 42, 407 (1979), DOI 10.1103/PhysRevLett.42.407 (Crossref-checked; not opened).
- *Computed*: T_ν = (4/11)^{1/3}·2.7255 K = 1.9454 K. n_γ = 410.7 cm⁻³. n_ν (3 flavours, ν+ν̄, instantaneous decoupling) = 336.0 cm⁻³, 1% below the PDG 339.5. Ω_ν/Ω_b = 0.031 (DESI) and 0.058 (Planck). "4.5× baryons" would need Σm_ν = 9.4 eV, 147× the DESI limit and ~7× the most KATRIN allows (Σ ≤ 3 × 0.45 = 1.35 eV for quasi-degenerate masses). Energy density at Σm = 0.064 eV: 1.16×10⁻¹² J/m³. Ram pressure on a perfect absorber moving at ~370 km/s: 1.8×10⁻¹⁸ Pa. Pauli ceiling for 0.05 eV neutrinos bound below 550 km/s: 3.4 cm⁻³ per species (g = 2), below the cosmic mean of ~56 per state.

### Inferences
- A waterspout needs upward accelerations of order g. Neutrino forces fall short by ≥ 5.8×10²⁷ (DGN's upper estimate), and even total absorption of their ram pressure gives 10⁻¹⁸ Pa. The neutrinos cannot drive water.
- "Consciousness imparts information into the condensate" (l. 1206–1212) has no measurable content and is not assessed.

### Gaps
- I did not open Tremaine–Gunn's paper. The Pauli ceiling above is the generic phase-space bound, not their numerical coefficient.

## Vacuum energy capture and Casimir forces

### Takeaway
The Casimir force is real and measured. It is 1.3 mPa at 1 µm and 1.3 bar at 10 nm. The energy released by letting plates collapse is a one-time π²ħc/(720 d³) per unit area: 4.3×10⁻⁴ J/m² at 10 nm, so matching 1 kg of H₂/O₂ would need about 3×10⁴ km² of plates. Because the force is conservative (it comes from an energy E(d)), any closed cycle returns zero net work. Forward's 1984 "vacuum-fluctuation battery" is explicitly recharged by putting the work back. Moddel's Casimir-cavity measurements report conductance changes attributed to zero-point fluctuations, not net power from a cycle. Status: Casimir force *mainstream*; cyclic vacuum-energy extraction *unverified* (and contradicted by the conservative-force argument); the talk's "power over thin wires, shielding and propulsion" is *unverified*.

### Cited Findings
- Forward, "Extracting electrical energy from the vacuum by cohesion of charged foliated conductors", Phys. Rev. B 30, 1700–1702 (1984), DOI 10.1103/PhysRevB.30.1700 (Crossref-checked). The abstract: "Any pair of conducting plates at close distances (< μm) experience an attractive Casimir force". "The battery can be recharged by making the electrical forces slightly stronger than the Casimir force" — [APS abstract](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.30.1700)
- Moddel, Weerakkody, Doroski et al., "Casimir-cavity-induced conductance changes", Phys. Rev. Research 3, L022007 (2021), DOI 10.1103/PhysRevResearch.3.L022007 (Crossref-checked). They "interpret the conductance changes in terms of a ΔEΔt uncertainty-principle-like limit to the injection of hot carriers from zero-point fluctuations"; the abstract has no energy-balance discussion — [APS abstract](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.3.L022007)
- Related, not opened: Moddel et al., "Optical-Cavity-Induced Current", Symmetry 13, 517 (2021), DOI 10.3390/sym13030517. Dmitriyeva & Moddel, "Test of zero-point energy emission from gases flowing through Casimir cavities", Physics Procedia 38, 8 (2012), DOI 10.1016/j.phpro.2012.08.007.
- Lamoreaux, "Demonstration of the Casimir Force in the 0.6 to 6 µm Range", PRL 78, 5 (1997), DOI 10.1103/PhysRevLett.78.5 (Crossref-checked; not opened).
- Repository: `gravitics::spacetime::casimir_pressure`, `casimir_energy_density`, `casimir_sphere_plate`, `zero_point_energy_density`, `quantum_inequality` already exist (crates/gravitics/src/spacetime.rs).
- *Computed*: P(1 µm) = 1.300×10⁻³ Pa; P(10 nm) = 1.300×10⁵ Pa; E/A(10 nm) = 4.334×10⁻⁴ J/m²; E/A(1 µm) = 4.334×10⁻¹⁰ J/m².

### Inferences
- The only thing missing from the crate is the one-shot collapse energy per area, `casimir_energy_per_area(d) = π²ħc/(720 d³)`. It equals `casimir_energy_density(d)·d`, so it is a one-line addition and possibly unnecessary.

### Gaps
- I did not open a 2022–2026 Moddel paper claiming net power output, so I cannot say whether such a claim exists or whether it was replicated.

## High-frequency gravitational waves from the craft

### Takeaway
The quadrupole formula (the repository's `gravitational::laboratory`) caps GW emission from anything craft-sized at absurdly small powers. A 5000 kg, 7 m rod spinning at 10 rev/s gives 4.5×10⁻³³ W. The plates' ~1 kJ field energy oscillating across 3.5 m at 90 kHz gives of order 10⁻⁴⁴ W. GW emission from the craft is real in principle (*mainstream*) but undetectable by 20+ orders of magnitude. As momentum it is bounded by P/c anyway, so it cannot propel.

### Cited Findings
- Repository: `gravitational::laboratory::spinning_rod_power = (2/45) G m² l⁴ ω⁶/c⁵`, `rotor_power`, `quadrupole_strain`, `burst_energy`. `gravitational::hf` holds the cosmological ceiling Ω h² ≈ 10⁻⁶ and the published HFGW limits (Aggarwal et al. 2025 Living Review, arXiv:2501.11723) (crates/gravitational/src/laboratory.rs, hf.rs).
- *Computed*: G/c⁵ = 2.756×10⁻⁵³ W⁻¹. Rod: 4.52×10⁻³³ W. Field-energy quadrupole: (G/c⁵)(M L² ω³)² with M = 1043 J/c² gives 1.8×10⁻⁴⁴ W (order of magnitude only).

### Inferences
- No new calculator is needed. The existing ones cover this.

### Gaps
- None for scope.

## "Ahead of its own light cone": what a simulator showing that means

### Takeaway
A simulator that integrates v ← v + a·dt (Newtonian) with constant thrust passes c after t = c/a. The user's audit code does exactly this (no Lorentz factor). With its own numbers (a ≈ 2.06 g) it passes c after 172 days. Seeing the craft "in front of its own light cone" is therefore a numerical artefact of a non-relativistic integrator, not a physical prediction. With relativistic kinematics, rapidity grows linearly and v < c always.

### Cited Findings
- Talk: "The weird thing is the simulation actually it seems to get in front of its own light cone" (transcript l. 2215–2219).
- User audit `ArvFlightSimulator::step`: `self.velocity = self.velocity.add(acceleration.scale(dt_sec))` with no relativistic correction — [user audit, unverified](../../../external/2026-10-26-greenyer-audit-and-crate.md)

### Inferences
- If superluminal travel is wanted, the repository's warp calculators (`gravitics::spacetime::alcubierre_energy`) give the negative-energy cost.

### Gaps
- The talk's own simulator code (O-Animator) was not available.

## The user's audit: its simulator checked

### Takeaway
With its own defaults the audit's simulator gives F = P/v_eff = 3.18×10¹⁰ W / 3.15×10⁵ m/s = 1.0×10⁵ N, i.e. 2.06 g on 5000 kg. Inputs: 1 kg of oxygen per 33 days at the "baryon collapse" yield of 9.05×10¹⁶ J/kg, and v_eff = 90 kHz × 3.5 m = 315 km/s. In vacuum this is 952 times the photon-rocket ceiling P/c. It also exceeds by 959× the most the burned mass could give if all of it left at c (ṁc = 105 N). With no matter or radiation leaving the craft, F = P/v_eff creates momentum from nothing, which violates momentum conservation. With the corrected 4 ¹⁷O → Fe+C yield the same formula gives only 91 N. The audit's "1 µF at 2.5 MV over 250 ns" is a 10 MA displacement current storing 3.1 MJ (12.5 TW if dumped in 250 ns), and 1 µF at a 1 m gap needs 0.11 km² of plate.

### Cited Findings
- Audit code (unverified data): `craft_mass_kg: 5_000.0`, `oxygen_fuel_kg: 950.0`, `oxygen_burn_rate_kg_s: 1.0 / (33.0 * 86_400.0)`, `frequency_hz: 90_000.0`, `major_radius: 3.5`. `let v_eff = (self.frequency_hz * self.core_monopole.major_radius).max(1.0); let thrust_newtons = power_watts / v_eff;`. Hover sets `mode_efficiency` to 0 and velocity is held. In atmosphere and water modes the burn is counted but tanks are not depleted. `effective_capacitance_f = 1.0e-6`, `dv_dt = potential_volts / pulse_rise_time_s`. Q values used: 49.47 (4 ¹⁷O → Fe+C), 17.13 (2 ¹⁷O → C+²²Ne), 16 × 938.272 MeV per ¹⁶O. Specific energy uses 17 u and 16 u per oxygen atom — [user audit](../../../external/2026-10-26-greenyer-audit-and-crate.md)
- *Computed* with the audit's numbers: ṁ = 3.507×10⁻⁷ kg/s, P = 3.175×10¹⁰ W, F = 1.008×10⁵ N, a = 2.056 g, F/(P/c) = c/v_eff = 951.7, ṁc = 105.1 N. Tank life at that rate is 85.8 years (not 28 days). The Newtonian velocity passes c after 172 days. Corrected 4 ¹⁷O pathway (57.37 MeV): P = 2.86×10⁷ W, F = 90.6 N. Displacement current: 1.0×10⁷ A; ½CV² = 3.125 MJ; 1.25×10¹³ W over 250 ns; ε₀A/d = 1 µF needs A/d = 1.13×10⁵ m.
- *Computed* air-breathing comparison: ideal actuator-disk thrust from 31.8 GW on a 38.5 m² disk in sea-level air is 4.6×10⁷ N, and hovering 5 t needs 1.1 MW. In air, the audit's 10⁵ N is allowed by momentum conservation (air is the reaction mass); in vacuum it is not.

### Inferences
- F = P/v is the power delivered to a body already moving at speed v (P = F·v). Used as a thrust law at fixed P with v_eff a constant, it ignores where the momentum goes. A consistent vacuum model must expel momentum: F = ṁ_exhaust·u_e with energy ½ṁu_e² ≤ P (non-relativistic), giving F ≤ 2P/u_e. Pure radiation gives F ≤ P/c. A fuel-mass-limited exhaust is bounded by Δv = c tanh((u/c) ln R).
- To fix the simulator with the least code: in vacuum use the relativistic rocket (thrust ≤ P/c for radiation, or ṁ_e u_e with exhaust mass taken from the tank). Integrate rapidity instead of velocity. In air, cap thrust with actuator-disk momentum theory. "Hover holds velocity" is fine in free space (no thrust needed to coast) but wrong in a gravity field: holding altitude needs F = mg continuously, at the actuator-disk power above.
- Status of the audit's model: *refuted* in vacuum (it violates momentum conservation by ~950×). In air it is *unverified but allowed in magnitude* (the mechanism is unspecified).

### Gaps
- The audit's "Srinivasan" step 2 ²⁸Si → ⁵⁶Fe + 2e⁺ + 2ν (+10.92 MeV) is left to the masses researcher. A quick check shows that its electron and positron bookkeeping needs care (atomic 2M(²⁸Si) − M(⁵⁶Fe) = 17.62 MeV, before subtracting 4 m_e c² for the positrons and spare electrons).

## Calculators for the crates

Only items not already in the repository. Existing calculators reused: `gravitics::propulsion::{photon_thrust, ion_wind_thrust, ballistic_thrust}`, `gravitics::spacetime::{induced_acceleration, casimir_pressure, casimir_energy_density, alcubierre_energy}`, `gravitational::laboratory::{spinning_rod_power, rotor_power}`, `magnetics::mhd::magnetic_pressure`, `magnetics::magnets::stored_energy`.

| Function signature | Formula | Primary source | Published / test value | Target crate::module | Status |
|---|---|---|---|---|---|
| `relativistic_rocket_dv(u: f64, mass_ratio: f64) -> f64` | Δv = c·tanh((u/c) ln R) | Ackeret 1946 (standard; derived in the script) | u = c, R = 5000/4050 → 0.2077 c, equal to (R²−1)/(R²+1) | gravitics::propulsion | mainstream |
| `exhaust_speed(eps: f64) -> f64` | u = c√(ε(2−ε)), ε = Q/(m c²) | energy–momentum conservation | ε = 5.551×10⁻⁴ (¹⁶O+¹⁶O → ³²S) → 9.988×10⁶ m/s; ε = 1.494×10⁻¹⁰ (H₂/O₂) → 5181 m/s | gravitics::propulsion | mainstream |
| `photon_rocket_proper_time(a: f64, mass_ratio: f64) -> f64` | τ = (c/a) ln R | special relativity | 1 g, R = 1.2346 → 74.56 days | gravitics::propulsion | mainstream |
| `actuator_disk_power(thrust: f64, rho: f64, area: f64) -> f64` | P = F^{3/2}/√(2ρA) | momentum theory (Rankine–Froude) | 49 033 N, 1.225 kg/m³, 38.48 m² → 1.118 MW | gravitics::propulsion | mainstream |
| `reaction_q(lhs: &[f64], rhs: &[f64]) -> f64` (mass excesses, keV) | Q = ΣΔ_in − ΣΔ_out | AME2020 (Wang et al. 2021) | ¹⁶O+¹⁶O → ³²S 16.5415 MeV; 4 ¹⁷O → ⁵⁶Fe+¹²C 57.372 MeV; 2 ¹⁷O → ¹²C+²²Ne 6.407 MeV | plasmoid::fusion (shared with the masses researcher) | mainstream |
| `specific_energy(q_kev: f64, reactant_mass_u: f64) -> f64` | e = Q/(m_reactants) (J/kg) | AME2020; NIST ΔfH°(H₂O,g) = −241.826 kJ/mol | O+O → ³²S 4.989×10¹³ J/kg; H₂/O₂ 1.342×10⁷ J/kg | plasmoid::fusion | mainstream |
| `relic_neutrino_temperature(t_cmb: f64) -> f64` | T_ν = (4/11)^{1/3} T_CMB | PDG 2024 ch. 26 | 2.7255 K → 1.9454 K (PDG "≃ 1.9 K") | gravitational::cosmology | mainstream |
| `relic_neutrino_density(t_cmb: f64) -> f64` | n = 3·(3/11)·2ζ(3)/π²·(kT/ħc)³ | PDG 2024 ch. 26 | 336.0 cm⁻³ (PDG 339.5 cm⁻³ with decoupling corrections, within 1.1%) | gravitational::cosmology | mainstream |
| `omega_nu_h2(sum_m_ev: f64) -> f64` | Σm_ν/93.14 eV | DESI DR2 Eq. 7 (PDG: 93.12 eV) | 0.064 eV → 6.871×10⁻⁴ (0.031 of Ω_b h² = 0.0224) | gravitational::cosmology | mainstream |
| `pauli_bound_density(m_ev: f64, v_max: f64, g: f64) -> f64` | n ≤ g(4π/3)(m v/h)³ | Tremaine & Gunn 1979 (phase-space argument) | 0.05 eV, 550 km/s, g = 2 → 3.39 cm⁻³ | gravitational::cosmology | mainstream |
| `relic_neutrino_acceleration_max() -> f64` (constant) | ≲ 17 × 10⁻²⁶ cm/s² | Duda, Gelmini & Nussinov 2001 | 1.7×10⁻²⁷ m/s², g/a = 5.8×10²⁷ | gravitational::cosmology | mainstream (bound) |
| `casimir_energy_per_area(d: f64) -> f64` | π²ħc/(720 d³) | Casimir 1948; Forward 1984 (one-shot) | 10 nm → 4.334×10⁻⁴ J/m²; 1 µm → 4.334×10⁻¹⁰ J/m² | gravitics::spacetime | mainstream; cyclic extraction unverified |
| `parallel_plate(area: f64, gap: f64, v: f64) -> (f64, f64, f64)` | C = ε₀A/d, U = ½CV², F = ε₀AV²/(2d²) | Maxwell (textbook) | 37.70 m², 1 m, 2.5 MV → 334 pF, 1043 J, 1043 N | gravitics::electrodes | mainstream |
| `displacement_current(c: f64, v: f64, t_rise: f64) -> f64` | I = C·V/t | Maxwell | 1 µF, 2.5 MV, 250 ns → 1.0×10⁷ A | gravitics::electrodes | mainstream |
| `ring_inductance(r: f64, a: f64) -> f64` | L = μ₀R(ln(8R/a) − 2) | standard thin-ring formula | 3.5 m, 0.4 m → 9.889 µH; ½LI² at 7.3 TA = 2.635×10²⁰ J = 3.09 × (950 kg)c² | magnetics::magnets | mainstream |
| `alfven_current(beta_gamma: f64) -> f64` | I_A = 4πε₀m_ec³βγ/e | Alfvén 1939 | βγ = 1 → 17.045 kA; 7.3 TA needs βγ = 4.28×10⁸ | plasmoid::pinch | mainstream |

Status summary: Holt field resonance *unverified* (no quantities); Forward 1963 *mainstream*; King MHD and BR fusion patents *mainstream physics, unbuilt*; oxygen chemistry and O+O Q values *mainstream*; O+O fusion as an ambient-condition fuel and baryon collapse *unverified*; ARV > c and "28 days" with 950 kg at the implied thrusts *refuted* by rocket accounting; 7.3 TA *refuted* by energy accounting; relic-neutrino "4–5× ordinary matter", "condensate" and waterspout forcing *refuted*; Casimir force *mainstream*, cyclic vacuum-energy capture *unverified*; HFGW emission *mainstream but negligible*; the audit simulator's vacuum thrust law *refuted* (momentum conservation).
