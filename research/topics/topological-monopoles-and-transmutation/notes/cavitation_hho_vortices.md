# Cavitation, oxyhydrogen (HHO / Brown's gas / Ohmasa gas) and vortex phenomena in Greenyer's "topological monopole" talk

Scope: the established physics, with formulas and published test values, behind the cavitation,
oxyhydrogen and whirlwind items in B. W. Greenyer, "The Power and Signature of the Topological
Monopole" (TeslaTech, 15 Aug 2026; transcript `/home/dan/Desktop/transcript_1793004282390.txt`),
plus the primary source and a status label for each claim. Labels: **mainstream**, **contested**,
**refuted**, **unverified**. Every computed number below is asserted in
`cavitation_hho_vortices_checks.py` (standard library, `python3 -I`, about 0.3 s, all checks pass).
Downloads (untrusted, read only) are in `dl_cavitation/`. "Opened" means I read the text (full text,
abstract via OpenAlex/Europe PMC, or patent text); DOIs were checked against Crossref.

## 1. Bubble dynamics, sonoluminescence temperatures, erosion jets, and whether collapse reaches nuclear energies

### Takeaway
Rayleigh collapse, the Rayleigh–Plesset equation, jet impact and water hammer explain cavitation
erosion with jets of about 100 m/s and pressures from 0.1 GPa up to several GPa (mainstream). The
highest measured single-bubble temperatures are 4,000–15,000 K (about 1.3 eV). At 15,000 K the D–D
Gamow exponent is about 173, so reactivity is about 10⁻⁹³ m³/s. That is 64 orders of magnitude
below the 10⁶–10⁷ K that Taleyarkhan's simulations assumed. Taleyarkhan's bubble-fusion claim
failed independent replication, and a 2008 Purdue committee found research misconduct in follow-up
claims: **refuted** as a replicated result.

### Cited Findings
- **Rayleigh collapse time** `t_c = τ R₀√(ρ/Δp)` with τ = 0.91468 for an empty spherical cavity in an inviscid, incompressible liquid. The original is Lord Rayleigh, Phil. Mag. 34, 94 (1917) ([doi](https://doi.org/10.1080/14786440808635681); Crossref-checked, not opened). The constant and its Gamma-function form are from Klotz, "Bubble dynamics in N dimensions", Phys. Fluids 25, 082109 (2013), which says "For N = 3, τ reduces to the familiar result of 0.91468" ([arXiv:1206.5753](https://arxiv.org/abs/1206.5753), opened; [doi](https://doi.org/10.1063/1.4817803)). Script: τ = √(3/2)·B(5/6,½)/3 = 0.914681. Integrating the Rayleigh–Plesset energy form `Ṙ² = (2Δp/3ρ)((R₀/R)³−1)` gives the same value to 10⁻³.
- **The SBSL review.** Brenner, Hilgenfeldt & Lohse, Rev. Mod. Phys. 74, 425 (2002), abstract (opened via OpenAlex): the evidence "favors a description of sonoluminescence caused by adiabatic heating of the bubble at collapse, leading to partial ionization of the gas … and to thermal emission such as bremsstrahlung", and "Shock waves inside the bubble do not seem to play a prominent role" ([doi](https://doi.org/10.1103/RevModPhys.74.425)). Full text not opened (bot-blocked), so no typical R₀, P_a or peak values are quoted from it.
- **Measured temperatures.** Flannigan & Suslick, Nature 434, 52 (2005), abstract (opened via [Europe PMC](https://europepmc.org/article/MED/15744295)). The study used SBSL in concentrated aqueous H₂SO₄ with Ar, and measured Ar emission temperatures of "4,000–15,000 K". The Ar excited states are ">13 eV" and "cannot be thermally populated" at those temperatures, and O₂⁺ cannot be made thermally, "We therefore conclude that these emitting species must originate from collisions with high-energy electrons, ions or particles from a hot plasma core" ([doi](https://doi.org/10.1038/nature03361)). Script: kT(15,000 K) = 1.2926 eV, and exp(−13 eV/kT) = 4.3×10⁻⁵.
- **Vapour limits the temperature.** Didenko & Suslick, Nature 418, 394 (2002), abstract: "temperatures attained in single-bubble cavitation in liquids with significant vapour pressures will be substantially limited by the endothermic chemical reactions of the polyatomic species inside the collapsing bubble" ([doi](https://doi.org/10.1038/nature00895)).
- **Jet speeds.** Plesset & Chapman, J. Fluid Mech. 47, 283 (1971), abstract: velocities scale as (Δp/ρ)^½. For Δp/ρ = 10⁶ cm²/s² (≈1 atm over the density of water) the jet reached "about 130m/sec" for a bubble touching the wall and "170m/sec" for one half a radius from it, "of a magnitude which can explain cavitation damage" ([doi](https://doi.org/10.1017/S0022112071001058)). Script: jet coefficients 13 and 17 times √(Δp/ρ).
- **Erosion pressures.** Philipp & Lauterborn, J. Fluid Mech. 361, 75 (1998), abstract: damage occurs at γ = s/R_max < 2. Jet impact matters only at γ ≲ 0.7, where "the impact velocity rises to 83 m s⁻¹, corresponding to a water hammer pressure of about 0.1 GPa". The largest erosive force comes from collapse in contact with the wall, "where pressures of up to several GPa act". At γ ≲ 1.7 the bubble becomes toroidal and the damage is a ring of pits ([doi](https://doi.org/10.1017/S0022112098008738)). Script: Joukowsky `p = ρcv` = 998 × 1483 × 83 = 0.12 GPa.
- **Toroidal collapse is ordinary hydrodynamics.** Lauterborn & Kurz, Rep. Prog. Phys. 73, 106501 (2010), abstract: aspherical collapse near boundaries "introduces a new effect: jet formation, the self-piercing of a bubble with a high velocity liquid jet" ([doi](https://doi.org/10.1088/0034-4885/73/10/106501)).
- **Taleyarkhan et al., Science 295, 1868 (2002).** The abstract reports tritium and "evidence for neutron emission near 2.5 million electron volts" in cavitated deuterated acetone, with "Hydrodynamic shock code simulations … (10⁶ to 10⁷ kelvin)" ([doi](https://doi.org/10.1126/science.1067589)).
- **Shapira & Saltmarsh, PRL 89, 104302 (2002).** Using "the same cavitation apparatus" and a larger detector, they found "no evidence for 2.5-MeV neutron emission correlated with sonoluminescence". Any emission is "at least 4 orders of magnitude too small" to explain the reported tritium ([doi](https://doi.org/10.1103/PhysRevLett.89.104302)). A Technical Comment by the same authors appeared in Science 297, 1603 (2002) ([doi](https://doi.org/10.1126/science.297.5587.1603a)).
- **Camara, Hopkins, Suslick & Putterman, PRL 98, 064301 (2007).** With a 15 atm drive and neutron-generator seeding, the upper bound is "a factor of 10 000 less than the signal reported by Taleyarkhan et al." ([doi](https://doi.org/10.1103/PhysRevLett.98.064301)). Naranjo's Comment is PRL 97, 149403 (2006) ([doi](https://doi.org/10.1103/PhysRevLett.97.149403)).
- **Misconduct finding** (secondary, Wikipedia's chronology): on 18 July 2008 a Purdue committee found two of twelve allegations substantiated. These were overstating independent confirmation and adding an author to answer a reviewer. The committee did not investigate the 2002 paper itself. Sanctions followed in August 2008, and the Office of Naval Research debarred him until September 2011, citing Reich, Nature news, 23 Nov 2009 ([Bubble fusion](https://en.wikipedia.org/wiki/Bubble_fusion); not checked against Purdue primary documents).
- **Gamow suppression,** computed with the repository's `plasmoid::fusion` formulas: the Bosch–Hale D(d,n) S-factor and reactivity fits, ported to Python and matched to the repository test value σ(50 keV) = 16.487 mb.
  - E_G = (31.397 keV^½)² = 985.8 keV. The Maxwellian ⟨σv⟩ integrated numerically from the S-factor agrees with the Bosch–Hale reactivity fit at 1 keV to 1%.
  - At 15,000 K the Gamow peak is at 74.4 eV with exponent 3E₀/kT = 172.7, and ⟨σv⟩_DDn = 1.2×10⁻⁹³ m³/s.
  - Assume a generous hot spot: n_D = 10²⁹ m⁻³ (above liquid density), radius 1 µm, lasting 1 ns. That gives 2.5×10⁻⁶² fusions per collapse.
  - At Taleyarkhan's assumed 10⁷ K (0.862 keV), ⟨σv⟩ = 4.2×10⁻²⁹ m³/s and the same hot spot gives about 880 per collapse.

### Inferences
- **Fusion is a temperature question.** The measured temperatures (≤15,000 K, 1.3 eV) sit 3 orders of magnitude below the 10⁶–10⁷ K that bubble-fusion needed. Reactivity falls by more than 10⁶⁴ between the two. Whether D–D fusion occurs in a collapse therefore turns entirely on whether 10⁶–10⁷ K is reached. Spectroscopy has not shown it, and three replication attempts saw no neutrons at 10⁻⁴ of the claimed level.
- **Ideal compression.** Adiabatic compression of argon (γ = 5/3, `T = T₀(R₀/R)^{3(γ−1)}`) reaches 15,000 K from 300 K at a radius ratio of 7.07. That is a modest compression, consistent with the review's adiabatic-heating picture.
- **Toroids and pit rings.** The "toroidal structures" and "heart shapes" Greenyer shows on cavitated metal (transcript lines 1296–1323, 1542–1560) match in form what Philipp & Lauterborn report for γ ≲ 1.7: toroidal bubbles and circular damage rings. The standard explanation is jet-driven toroidal collapse, not a magnetic monopole. Label: the morphology is **mainstream**; the monopole interpretation is **unverified**.
- **Kladov's energy claim fails by six orders.** Kladov claims bubble collapse exceeds nucleon binding energies (section 2). The measured hot-spot energies are kT ≈ 1.3 eV, about 1.6×10⁻⁷ of 8 MeV per nucleon.

### Gaps
- I could not open the Brenner–Hilgenfeldt–Lohse full text, so typical R₀, P_a and the hard-core radius are not quoted.
- Flannigan & Suslick, Nature Physics 6, 598 (2010) ([doi](https://doi.org/10.1038/nphys1701)) is reported in press coverage as >16,000 K. I did not open it.
- The Purdue misconduct report and the ONR debarment notice were not read in primary form.

## 2. Kladov (2003) "cavitation transmutation/destruction of matter", Claytor (LANL), and the Huang et al. water-cavitation experiment

### Takeaway
I found no primary "Kladov 2003" text. The citable primary is Kladov's Russian patent RU 2054604 C1
(filed 1993, published 1996). It claims cavitation in an "ultrasonic activator" releases "the
energy of strong interactions", with heat output 2–3 times power input. Its one fully stated
example has N₂/N₁ = 6.5/6 = 1.08, within ordinary calorimetric error of a pump that turns all its
shaft work into heat: **unverified**. I found no cavitation or transmutation publication by
T. Claytor. Huang et al. (Sci. Rep. 14, 214, 2024, with Greenyer as co-author) runs pulsed hot-water
flow in heat exchangers with pressure spikes above 720 bar. Collapse of vapour cavities
(condensation-induced water hammer) explains spikes of that size by conventional hydrodynamics. The
paper's nuclear interpretation carries an editor's note that it is under criticism: **contested**.

### Cited Findings
- **Kladov's patent.** RU 2054604 C1, "Energy generation method", inventor А. Ф. Кладов, filed 1993-07-02, published 1996-02-20 (Google Patents machine translation, opened; [patent](https://patents.google.com/patent/RU2054604C1/en)).
  - **Pressure rule:** cavitation under constant P₁ and alternating P₂ with `P₁ = (0.3–0.7)(P₂+P₃)` and `P₂+P₃−P₁ = (1–10)σ`. P₃ is the vapour pressure and σ the liquid's tensile strength.
  - **Claimed mechanism:** the collapse energy is "several orders of magnitude higher than the binding energy of elementary particles (nucleons) in the nucleus".
  - **Example 1:** tap water at 20 °C, σ ≈ 0.35 MPa, P₁ = 0.15 MPa, P₂ = 0.5 MPa, electrical input N₁ = 6 kW, heat N₂ = 6.5 kW, "a coefficient (k) of 1.08". The heat is measured as `N = Δt·G·c`.
  - **General claims:** the text says power consumption is "on average two to three times less than the heat output", and that ionizing radiation "does not exceed the natural background". It also reports adding stable isotopes and measuring γ, β and neutron counts (figures 6–9, not reproduced in the text).
  - **Device:** a multi-stage centrifugal pump with perforated rotor and stator rings, a rotary cavitation heater design.
  - Script: P₂ = 0.5 MPa and P₁ = 0.15 MPa follow from the rule, and k = 1.083.
- **Context in the talk.** Greenyer lists "Anatoly Klyosov's … cavitation, transmutation, and destruction of matter from 2003" (transcript ~line 224; the slide says "A. Kladov (2003)").
- **Claytor.** Searches found only that T. N. Claytor (LANL 1986–2012) is known for tritium production in a modified low-energy experiment, and that he received the 2012 Preparata Medal ([Los Alamos Daily Post](https://ladailypost.com/los-alamos-scientist-awarded-giuliano-preparata-medal/), search snippet only, not opened). I found nothing linking him to cavitation.
- **Huang et al. bibliographic details** (opened, [PMC10762109](https://pmc.ncbi.nlm.nih.gov/articles/PMC10762109/); [doi](https://doi.org/10.1038/s41598-023-50824-8), Crossref: Sci. Rep. 14, 214, 2024). The 19 authors include Robert William Greenyer (Martin Fleischmann Memorial Project), who "supervised physical reaction mechanism". The editor's note (5 Jan 2024): "Readers are alerted that the conclusions presented in this article are subject to criticisms". Crossref shows no retraction or update record as of 2026-10-26.
- **Huang et al. apparatus.**
  - Reactor VCS: a triple-pipe heat exchanger about 30 m long, with water at about 21 bar in a gap of about 2 mm, heated by R22 vapour from a 3 kW compressor (150–160 °C outlet). Flow was pulsed 2–10 cycles/min at an average 1.2 L/min.
  - Later reactors used a once-through electric boiler (≤10 kW, water up to 190 °C) and "resonator" sections.
  - The paper credits cavitation to "implosion of water vapor bubbles" but gives no pump or nozzle specification.
- **Huang et al. reported results.**
  - Peak COPx was 4.26 (VCS) with "pressures above 720 bar that ruptured gauges and pipes", and 2.55 (DHX) above 240 bar.
  - The gas-sampled runs had COPx 1.0–1.61, with excess heat classed as COPx > 1.05.
  - Gas analysis used quadrupole MS (Extorr XT200M/XT300M, better than 0.5 amu resolution).
  - The paper itself describes the mechanism as speculative.

### Inferences
- **Huang's pressure spikes.**
  - Stopping a flow suddenly makes a Joukowsky pressure `Δp = ρcΔv`. Reaching 720 bar needs only Δv = 48.6 m/s of arrested liquid.
  - A vapour void collapsing under the stated 21 bar reaches U = 99 m/s by R₀/R = 2 (Rayleigh energy form), and its water hammer is about 1,465 bar.
  - Pulsed flow of water near or above saturation in a long heat exchanger is the classic setting for condensation-induced water hammer.
  - So the reported >720 bar spikes and ruptured gauges have a standard hydrodynamic explanation and are not, by themselves, evidence of nuclear energy release. The same transients load the calorimetry, so the COPx values need careful error analysis. The nuclear accounting is covered by another researcher.
- **Kladov's numbers.** The pressure work available, (P₁+P₂) times the volume per molecule, is about 1.2×10⁻⁴ eV per water molecule. Even the hottest measured collapse (≈1.3 eV, section 1) is about 10⁻⁷ of nuclear binding per nucleon. The patent's mechanism is inconsistent with the energy densities its own pressures set. Its quantified example (k = 1.08) is consistent with COP = 1 inside ordinary measurement error for an electrically driven pump, since all shaft work ends up as heat. Russian rotary "cavitation heat generators" with COP > 1 claims are a large patent family (the patent page lists about 50 later citing patents), but I found no independent calorimetry.

### Gaps
- The primary "Kladov 2003" text (perhaps a conference paper titled "Cavitation destruction of matter") was not found. Secondary blog posts quote it but are not citable.
- The patent's radiation figures (6–9) and its table of examples 2–12 were not readable in the machine translation.
- I found no Claytor cavitation paper, and could not identify the "heart shape" replication Greenyer cites.
- Huang et al. give no pressure-transducer traces or a pump/nozzle description, so the water-hammer explanation cannot be tested against their data.

## 3. Oxyhydrogen: properties, electrolysis, Brown's-gas claims, Santilli's "magnecules", HHO cutting

### Takeaway
A stoichiometric 2H₂:O₂ mixture has an equilibrium adiabatic flame temperature of about 3,075 K at
1 atm. That value is published (CEARUN) and reproduced here from GRI-Mech 3.0 polynomials. It is
far below the 4,890–4,931 K that ignoring dissociation would give. Electrolysis needs at least
1.229 V (reversible) and 1.481 V (thermoneutral, HHV). A cell at 2.0 V stores at most 74% of its
electrical input as gas heating value, so any "COP > 1" from burning electrolysis gas contradicts
the first law. Brown's patents describe ordinary oxyhydrogen and atomic-hydrogen welding. The
"implosion" is the 1,866-fold volume collapse when the burnt gas condenses. Radiation remediation,
Santilli's "magnecules", and "0.6 m of steel cut with HHO in China" are **unverified** (Santilli
**contested**).

### Cited Findings
- **Thermochemistry** (NIST WebBook, CODATA values of Cox, Wagman & Medvedev, opened): ΔfH°(H₂O, l) = −285.830 ± 0.040 kJ/mol; ΔfH°(H₂O, g) = −241.826 ± 0.040 kJ/mol; S°(H₂O, l) = 69.95, S°(H₂) = 130.680, S°(O₂) = 205.152 J/(mol K) ([water](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Units=SI&Mask=3), [H₂](https://webbook.nist.gov/cgi/cbook.cgi?ID=C1333740&Units=SI&Mask=3), [O₂](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7782447&Units=SI&Mask=3)). Script:
  - E_rev = −ΔG/2F = 1.2289 V and E_tn = −ΔH/2F = 1.4812 V.
  - HHV 141.79 MJ/kg and LHV 119.96 MJ/kg; 12.75 MJ per m³ of H₂ at 0 °C, 1 atm.
  - The maximum fuel-cell work is 237.14 kJ per mol H₂.
- **Adiabatic flame temperature.** Marzouk, Eng. Technol. Appl. Sci. Res. 13(4), 11437 (2023) ([doi](https://doi.org/10.48084/etasr.6132), opened via [arXiv:2503.11826](https://arxiv.org/abs/2503.11826)). For oxy-hydrogen from 298.15 K at 1 atm:
  - Equilibrium: 3,074.51 K (NASA CEARUN) and 3,076.92 K (GRI-Mech 3.0 in Cantera).
  - Complete combustion: 4,930.56 K (CEARUN) and 4,890.21 K (GRI-Mech, Excel).
  - Air–hydrogen equilibrium: 2,378.62 K (CEARUN).
  - Script (six species H₂O, H₂, O₂, OH, H, O; GRI-Mech 3.0 NASA-7 data from [thermo30.dat](http://combustion.berkeley.edu/gri-mech/version30/files30/thermo30.dat)): equilibrium 3,078.7 K (within 0.06% of Cantera), complete combustion 4,890.2 K. About 42% of the product moles are dissociated at the flame temperature. The polynomials reproduce NIST's ΔfH(H₂O, g) to 10⁻⁵.
- **Secondary oxyhydrogen properties** (Wikipedia, citing Calvert 2008, the NASA Glenn safety manual, and Moyle et al., AIChE J. 1960; [Oxyhydrogen](https://en.wikipedia.org/wiki/Oxyhydrogen), opened):
  - about 2,800 °C maximum at stoichiometry;
  - autoignition about 570 °C;
  - spark ignition energy about 0.007 mJ;
  - flammable at 4–95% H₂.
  - The article also records that "Brown's gas"/HHO claims are associated with pseudoscience (Ball, Nature news 2006, 2007).
- **Yull Brown's patents.** US 4,014,777 (1977) and US 4,081,656 (1978), priority 1973 (opened). They describe "welding, brazing or the like utilizing a mixture of hydrogen and oxygen generated in substantially stoichiometric proportions in an electrolytic cell", plus "atomic welding": the mixture passes through an arc between tungsten electrodes, and atomic H and O "on recombination generate an intensely hot flame". A later Brown patent, AU 590309 B2, is titled "Producing a vacuum by implosive reaction of hydrogen and oxygen" (title from the citation list on the [US4014777A page](https://patents.google.com/patent/US4014777A/en); not opened).
- **Santilli's HHO.** R. M. Santilli, "A new gaseous and combustible form of water", Int. J. Hydrogen Energy 31, 1113 (2006) ([doi](https://doi.org/10.1016/j.ijhydene.2005.11.006)). Published responses:
  - J. M. Calo, Comments, IJHE 32, 1309 (2007) ([doi](https://doi.org/10.1016/j.ijhydene.2006.11.004));
  - J. Kadeisvili, Rebuttal, IJHE 33, 918 (2008) ([doi](https://doi.org/10.1016/j.ijhydene.2007.10.030));
  - Cloonan, "A chemist's view of J.M. Calo's comments", IJHE 33, 922 (2008) ([doi](https://doi.org/10.1016/j.ijhydene.2007.11.009)).
  - All four were Crossref-checked; their contents were not opened.
- **Greenyer's claims.**
  - "HHO cutting through in China 0.6 of a meter of steel", attributed to "a coherent matter wave beam" (transcript lines 1415–1422).
  - "Yull Brown did radiation remediation" (line 1542).
  - Ohmasa gas "appears to fission heavy elements to safe elements … as was … claimed to have been shown in the 1980s by Yull Brown" (lines 1627–1636).

### Inferences
- **Energy bound.** At 2.0 V and 100% Faradaic efficiency, the gas carries 1.481/2.0 = 74% of the electrical input as HHV. Burning it returns at most that, so a closed electrolysis–combustion loop cannot have COP > 1 (**mainstream**). Fuel-cell use is bounded by 237 kJ/mol H₂.
- **"Variable flame temperature" on different targets.** Without dissociation the flame would be at 4,890–4,931 K; equilibrium dissociation (H, O, OH) brings it to 3,075 K. Radicals recombine on a cold or catalytic surface and release heat there, which is the principle of Brown's own atomic-welding claim and of Langmuir's atomic-hydrogen torch. Heat transfer then depends strongly on the target's conductivity and catalytic activity. That explains target-dependent behaviour without new physics: an inference, as I did not open a measurement.
- **"Implosion".** 2H₂ + O₂ (3 mol, 67.2 L at STP) burns to 2 H₂O, and once condensed the products occupy 36 cm³, a 1,866-fold shrink. A burnt stoichiometric charge therefore leaves a near-vacuum. That is the conventional reading of Brown's "implosive reaction" title (**mainstream** physics; the patent itself not opened).
- **HHO cutting 0.6 m of steel.** Oxy-fuel cutting is sustained by iron burning in an oxygen jet; the preheat flame only starts it, so the fuel gas matters little to the maximum thickness once the cut runs. Hydrogen is an established preheat gas. A 0.6 m oxygen cut is within the range industrial oxy-fuel cutting is generally described as covering (a search snippet from TWI says "1mm to 1000mm"; not opened). So the cut, if real, needs no "matter wave beam". No primary source for the China demonstration was found: **unverified**, with a conventional explanation available.

### Gaps
- No primary source found for the China 0.6 m demonstration, for Brown's radiation-remediation claims, or for a measured laminar flame speed of undiluted stoichiometric H₂/O₂ (heat-flux studies found use diluted mixtures).
- Oxy-fuel thickness limits are not quoted from an opened primary source. TWI's page refused connection and thefabricator.com refused (403).
- Calo's critique and the Santilli replies were not read.

## 4. Ryushin Ohmasa: patents, "Ohmasa gas", COP claims, and what the vibrator does

### Takeaway
Ohmasa's patents (as Omasa, assignee Japan Techno) describe alkaline electrolysis with
low-frequency mechanical vibration stirring: a 10–200 Hz motor drives blades at 0.1–15 mm amplitude
and 200–1,000 strokes per minute. This is a mixing and bubble-detachment technique (**mainstream**).
The one embodiment with stated voltage, current density and gas rate is consistent with ordinary
Faradaic electrolysis at 74% HHV efficiency (COP < 1). Two later patents make energy claims:
- a fuel cell giving "2 to 3.5 times the electrical power";
- a CO₂ + Ohmasa-gas product with heating value 44.8 MJ/m³.

The CO₂ claim exceeds the chemical bound of its hydrogen content by 6.3–8.6 times: **unverified**,
and inconsistent with energy conservation unless the gas analysis is wrong. I did not find the
patents with "COP 2–8.7" or "11.7".

### Cited Findings
- **US 7,459,071 B2** (priority 2001-05-02, granted 2008-12-02; [patent](https://patents.google.com/patent/US7459071B2/en), opened): "Hydrogen-oxygen gas generator…". Its operating conditions:
  - motor 10 to 200 Hz;
  - blade tip amplitude "0.1 to 15.0 millimeters and a frequency of 200 to 1,000 times per minute";
  - electrolyte 5–10% alkali hydroxide, pH 7–10, 20–70 °C, current density 5–20 A/dm²;
  - electrode spacing 1–20 mm (search summary of the abstract).
- **US 7,459,071 B2, first embodiment:**
  - stainless cell 270 × 1660 × 390 mm (H), 250 W three-phase vibration motor;
  - 50 Pt-plated Ti anodes and 50 Ti cathodes, KOH 8% at 55 °C, pH 10;
  - "2.0 volts (direct current)", "5 A/dm²", collection "1,000 liters per hour".
  - Later embodiments give 1,200–3,000 L/h. The text contains no cavitation and no over-unity claim.
- **CA 2,522,593 A1** ("Fuel for fuel cell…", [patent](https://patents.google.com/patent/CA2522593A1/en), opened):
  - Composition of the vibration-electrolysis gas: H₂ 55–70 mol%, H 0.12–0.45, H₃ and HD 0.03–0.14, OH 0.3–1.2, ¹⁶O 1.0–4.2, O₂ 5–27 mol%. These are from mass-spectrometer peaks, whose plotted heights for two masses were scaled up "100 times" and "10 times".
  - Claim: the gas is "capable of generating 2 to 3.5 times the electrical power" as fuel-cell fuel compared with conventional hydrogen. One example gives 0.5 W at 0.6 V against 0.15–0.2 W.
- **US 2016/0145521 A1** (priority 2013-06-21; [patent](https://patents.google.com/patent/US20160145521A1/en), opened):
  - Method: CO₂ plus "OHMASA-GAS" at 0.1–10 MPa and 5–50 °C, without catalyst, gives a "new combustible gas". Its heat is "44.8 MJ/m3 (10,700 Kcal/m3)" at CO₂ 20% / Ohmasa gas 80%, and "43.9 MJ/m3 (10,500 Kcal/m3)" at 40/60.
  - Evidence offered: a "bluish sharp flame", combustion temperatures of 300–500 °C, and "No liquid such as waterdrops and steam" after compression.
  - Background claims: the gas "contains H2O clusters" and "does not explode". It cites Japanese patents 3975467, 4076953 and 4599387.
- **Greenyer's claims:** "a matter transformation cloud off observed COP of 2 to 8.7 in his original patents and up to 11.7 in the patent I've translated" (transcript lines 1273–1276). Ohmasa vibrator plates "transmute both ions in solution and metals" (lines 1612–1622). Ohmasa gas on tungsten for "10 seconds" leaves "a strontium ball" (lines 488–502).

### Inferences
- **Faraday check of the first embodiment** (assuming litres at 0 °C, 1 atm; the patent does not say):
  - 1,000 L/h of 2:1 gas needs 1,594 A, which is 319 dm² of anode at 5 A/dm². That fits within the 526 dm² that 50 plates filling the cell cross-section could offer, so the numbers are self-consistent with ordinary electrolysis.
  - Electrical input is 2.0 V × 1,594 A ≈ 3.2 kW, plus the 250 W motor, giving HHV efficiency of 74% or less (COP < 1).
- **The CO₂ claim.** CO₂ has no heating value, so the feed's heating value is at most that of its H₂:
  - 20/80: 0.8 × (2/3) × 12.75 = 6.80 MJ/m³ (7.14 if the gas were 70% H₂). The claimed 44.8 MJ/m³ is 6.3–6.6 times that.
  - 40/60: 5.10–5.36 MJ/m³ against 43.9 claimed, 8.2–8.6 times.
  - Normalizing per m³ of product after a methanation-type reaction cannot help, because converting CO₂ + 4H₂ to CH₄ + 2H₂O makes water. The patent says no liquid appeared, and Sabatier chemistry needs a catalyst and heat that the patent excludes.
  - These ratios fall inside Greenyer's quoted "COP 2 to 8.7". That may be coincidence; I could not confirm which patents he meant.
- **Radicals in the gas.** H, OH and O at mole-percent levels in a stored room-temperature gas would recombine within microseconds to milliseconds. Peaks at m/z 1, 16 and 17 are standard electron-impact fragments of H₂O, and H₃⁺ at m/z 3 forms in ion sources. The claimed radical composition is therefore most simply a mass-spectrometry artefact (inference; no independent analysis found).
- **Does the vibrator cavitate?** Peak blade-tip speed is 2π × (1,000/60) × 15 mm = 1.57 m/s. The cavitation number σ = (p − p_v)/(½ρv²) is about 66, taking p_v up to 20 kPa as an assumed value above water's vapour pressure at 55 °C. Inception needs σ of order 1, so the blades do not cause hydrodynamic cavitation by speed alone. Pits on vibrator plates that Greenyer calls "cavitation spots" (lines 1545–1552) more likely come from electrolytic corrosion, gas-bubble attack or fatigue. Because the patent says blades "flutter", local edge cavitation in their vortices cannot be excluded without examining a plate.

### Gaps
- No patent or paper with COP 2–8.7 or 11.7 was found. Justia's inventor listing was blocked and searches returned only the patents above. Japanese originals (J-PlatPat) were not searched.
- I found no independent analysis of Ohmasa gas, of the 2019 tungsten test (Greenyer, Egely), or of vibrator-plate transmutation.

## 5. Vortex physics: tornado pressure and hydrogen bonds, waterspouts, dust-devil fields, "black water"

### Takeaway
Measured tornado winds reach about 132 m/s. A Rankine vortex then has a core deficit ρv² of about
200 hPa, twice the largest measured 100 hPa (Manchester, SD, 2003). These pressures and flows carry
about 1% of one water–water hydrogen bond of kinetic energy per molecule, and 10⁻⁵ of it as pressure
work on liquid water. Hydrogen bonds therefore do not "give up" in a tornado; airborne gypsum board
fails mechanically and by wetting (**mainstream**). The "hydrogen bonds giving up" mechanism is
**refuted** on energy grounds. Dust devils carry real charge: fields above 4.3 kV/m, ULF magnetic
signals about 0.1 nT, about 10⁶ charges per cm³ (**mainstream**). The dark waterspout patch is
Golden's "dark spot" stage, a sea-surface pattern (**mainstream**). Water "turned to water gas"
running "with the flux of relic neutrinos" is **unverified**, with no supporting measurement.

### Cited Findings
- **Water-dimer bond energy.** Rocher-Casterline et al., J. Chem. Phys. 134, 211101 (2011), abstract: water dimer "D(0) = 1105 ± 10 cm(-1) (13.2 ± 0.12 kJ/mol)" ([doi](https://doi.org/10.1063/1.3598339)). Script: 1105 cm⁻¹ = 13.22 kJ/mol.
- **Tornado winds.** Wurman et al., BAMS 88, 31 (2007), abstract: Doppler-on-Wheels-constrained winds "up to 132 m s⁻¹" ([doi](https://doi.org/10.1175/BAMS-88-1-31)).
- **Tornado pressure measurements.** Karstens et al., Mon. Wea. Rev. 138, 2570 (2010): probe pressure intercepts in nine tornadoes since 2002, with "a comparative range of observed pressure deficits" ([doi](https://doi.org/10.1175/2010MWR3201.1); abstract only, full text blocked).
- **The Manchester record** (secondary): the 24 June 2003 Manchester, SD tornado gave "the largest drop in atmospheric pressure ever recorded, 100 hPa (mb) in less than one minute", with Guinness listing it as the record ([Tim Samaras](https://en.wikipedia.org/wiki/Tim_Samaras)).
- **Rankine vortex** (cyclostrophic balance `dp/dr = ρv²/r`, solid-body core inside r_c, 1/r outside): the total centre deficit is ρv_max². The script integrates the profile and recovers ρv_max² to 10⁻³. At 132 m/s and ρ = 1.15 kg/m³ that is 200 hPa, and 100 hPa corresponds to 93 m/s.
- **Waterspout life cycle.** Golden, J. Appl. Meteor. 13, 676 (1974), abstract: five stages. The first is "the dark spot, a prominent light-colored disc on the sea surface surrounded by a dark patch … which represents a complete vortex column from cloud base to sea surface", followed by the spiral pattern, spray ring, mature spray vortex, and decay ([doi](https://doi.org/10.1175/1520-0450(1974)013<0676:TLCOFK>2.0.CO;2)).
- **Dust-devil electric fields.** Farrell et al., J. Geophys. Res. 109, E03004 (2004) ([doi](https://doi.org/10.1029/2003JE002088); [author copy](https://faculty.washington.edu/dcatling/papers_mine/Farrell-matador2004.pdf), opened):
  - In a direct pass the vertical field "exceeded 4.3 kV/m" and saturated the field mill. These are values "larger than any previously reported from a dust devil", against the earlier ≈0.6 kV/m of Freier (1960) and Crozier (1964).
  - Charge density "n ≈ 10⁶ el/cc" from the dipole estimate. "If this measured vertical E-field persists … the total vertical potential would be nearly 0.8 MV".
  - The fields are "0.1% of terrestrial atmospheric breakdown" and "10–20% of the Martian breakdown potential (of 20 kV/m)".
  - ULF magnetic field just outside the dust devil: B ≈ 0.1 nT. Diameter inferred: 10 m.
  - Script: Farrell's estimate `n ≈ 4πε₀r²E/(e·V/2)` with r = 50 m and a 10 m × 100–200 m column gives 0.95–1.9 × 10⁶ cm⁻³; 4.3/20 = 21.5%; and 0.8 MV/4.3 kV/m = 186 m of column.
- **Greenyer's claims:** gypsum "disassembling in the air" because "the hydrogen bonds are just giving up" (transcript lines 338–352). Waterspout water "partially … turned to water gas and is running with the flux of relic neutrinos", and "black water" in "that 1969 NOAA … photo" with smoke flares (lines 1236–1248).

### Inferences
- **Energy per molecule.** Water moving at 132 m/s has bulk kinetic energy of 0.157 kJ/mol, which is 1.19% of D₀. Thermal RT at 298 K is 16 times larger, and it does not break hydrogen bonds in liquid water. A 100 hPa pressure change does work of 1.4×10⁻⁵ D₀ on liquid water. Neither can dehydrate gypsum or break its bonds.
- **Gypsum board.** Its mid-air disintegration is most simply aerodynamic and impact loading of a weak, brittle, often wetted plaster panel. Dynamic pressure ½ρv² is about 100 hPa at 132 m/s, the same order as the core deficit.
- **"Black water".** The dark patch around a light disc is Golden's "dark spot" stage, a pattern of surface roughness and spray seen from above. The 1969 photo Greenyer shows probably comes from the NOAA Florida Keys programme that Golden analysed. I did not identify the photo, so the match is unconfirmed.
- **Vortex electrodynamics.** Dust devils are a documented natural case of a vortex charging and producing kV/m fields and nT-scale ULF signals, via triboelectric grain charging and mass separation. That is the measured electrodynamics of a whirlwind and needs no monopole.

### Gaps
- I found no published in-situ electric-field measurement inside a waterspout. The searches returned only Russian storm studies and over-ocean buoy field mills.
- The Karstens et al. table of pressure deficits and the Lee et al. 2004 Manchester probe paper were not opened.
- The 1969 NOAA photograph was not identified.

## 6. The "Polifka air grinder" (Windhexe) and comminution

### Takeaway
The device is the Windhexe of Francis D. "Frank" Polifka: US 6,971,594 B1 (priority 1999, granted
2005, expired). Compressed air at 10–600 psi is injected tangentially into a cylindrical-over-
conical chamber. The vortex makes particles collide with each other and the walls, which grinds and
dries materials including glass and aluminium. The patent says nothing about magnetic output. A
magnetic grey powder from glass and aluminium cans would most simply contain iron worn from the
chamber, which the patent makes of "heavy metal material" for "wear life". The grinding is
**mainstream**; "massive magnetic properties" is **unverified**.

### Cited Findings
- **The patent.** US 6,971,594 B1, "Apparatus and method for circular vortex air flow material grinding", inventor Francis D. Polifka, original assignee Vortex Dehydration Technology LLC, priority 1999-03-23, filed 2000-03-16, published 2005-12-06, "Expired – Fee Related" ([patent](https://patents.google.com/patent/US6971594B1/en), opened).
- **Apparatus:**
  - Compressed air enters through angled slots "relatively circumferentially … so as to generate a circular vortex flow of air and material … for grinding and drying".
  - Air at "about 10 to about 600 pounds per square inch (psi)" and "about 5 to about 12,000 cubic feet per minute".
  - Uses listed include glass, grain, paper, plastic, aluminum and granite. "Glass can be processed into a fine powder or silica."
  - The sidewalls and cover "are made of the same heavy metal material to ensure a relatively long wear life".
- **Greenyer's claims:** "particles coming out that have got a massive magnetic … properties", glass bottles "2 seconds later it comes out as a fine gray free flowing powder. Same with aluminum cans", "We're replicating that" (transcript lines 1225–1235).

### Inferences
- **Grinding mechanism.** Particle–particle and particle–wall impact at near-sonic jet speeds, as in conventional fluid-energy (jet) mills, is the standard explanation.
- **Comminution energy.** Rittinger's law (E ∝ new surface area, E = K_R(1/x₂ − 1/x₁)) is the textbook relation. I did not open a primary source for K_R or for glass's fracture surface energy, so no number is given.
- **Grey magnetic powder.** Glass (soda-lime silica) and aluminium are not ferromagnetic. A magnetic fraction points to iron from the steel walls, can coatings or contaminants. A simple test: magnetically separate the powder and measure Fe by XRF or ICP before and after, plus chamber-wall mass loss.

### Gaps
- No measurement of the powder's magnetism, composition or particle size was found. Secondary pages (Wikipedia "Windhexe", a 2002 Seattle Times article) were seen only as search snippets.
- No primary Rittinger constant or glass surface energy was opened.

## 7. Chladni plates and Helmholtz resonators (the standing-wave tool, transcript lines 138–148)

### Takeaway
Both are classical acoustics with published test data (**mainstream**). A Helmholtz resonator rings
at `f = (c/2π)√(A/(V L′))` with end correction `L′ = L + 1.5a` for an unflanged neck; a smartphone
tea-bottle experiment recovers c = 345(5) m/s. Plate modes follow Kirchhoff theory,
`ω = (λ²/a²)√(D/ρh)`, with λ² tabulated by Leissa (NASA SP-160). The script reproduces the
clamped-circular values from the Bessel frequency equation.

### Cited Findings
- **Helmholtz bottle.** Monteiro, Stari, Cabeza & Martí, "A bottle of tea as a universal Helmholtz resonator" ([arXiv:1805.04014](https://arxiv.org/abs/1805.04014), opened):
  - formula `f = (c/2π)√(A/(V L′))` with "L'=L+ 1.5 a" for an unflanged outer end, valid for dimensions much smaller than the wavelength;
  - neck d = 29.0(2) mm, L = 70(1) mm, cavity V = 435 cm³, wavelengths "between 1 and 1.5 m";
  - fitted c = 345(5) m/s against a "reference value, 344 m/s" at 22.3 °C.
  - Script: c = √(γRT/M) = 344.6 m/s at 22.3 °C; the full bottle rings at 223 Hz (λ = 1.54 m); λ = 1.5 m and 1.0 m need air volumes of 410 and 182 cm³.
- **Plate modes.** Leissa, *Vibration of Plates*, NASA SP-160 (1969) ([NTRS 19700009156](https://ntrs.nasa.gov/citations/19700009156), opened):
  - Table 2.1, clamped circular plate, λ² = ωa²√(ρh/D) (n nodal diameters, s nodal circles): (0,0) 10.2158 (the scan reads "10.2168"), (1,0) 21.26, (2,0) 34.88, (0,1) 39.771. These are roots of `J_n(λ)I_{n+1}(λ) + I_n(λ)J_{n+1}(λ) = 0`.
  - Table 2.5, completely free circular plate, ν = 0.33: (2,0) 5.253, (0,1) 9.084, (3,0) 12.23.
  - Script reproduces the four clamped values to ≤3×10⁻⁴.

### Inferences
- **Fit to the crates.** These formulas carry Greenyer's standing-wave and phonon-focusing language into testable numbers: resonance frequency, mode shapes and nodal radii. They belong with `gravimetry::atmosphere::sound_speed`, which already exists. Phase-conjugate "standing waves" that pick out where structures form (transcript lines 236–262) can be checked against computed nodal patterns of the actual cell or plate.

### Gaps
- Free-square-plate (classic Chladni) frequency parameters were not extracted; Leissa covers them in later chapters.
- The free-circular-plate frequency equation was not coded. It needs Poisson's ratio, and its Leissa values are quoted only.

## Calculators for the crates

Status describes the physics, not the fringe claim it tests. "Test value" is what the crate's doc
test should assert; each is asserted in `cavitation_hho_vortices_checks.py`.

| Function (proposed) | Formula | Primary source | Published test value | Target crate::module | Status |
|---|---|---|---|---|---|
| `rayleigh_collapse_time(r0, rho, dp) -> f64` | `t = 0.91468 R₀√(ρ/Δp)`, τ = √(3/2)·B(5/6,½)/3 | Rayleigh 1917; Klotz, Phys. Fluids 25, 082109 (2013) | τ = 0.91468 | plasmoid::cavitation (new) | mainstream |
| `rayleigh_wall_speed(r0, r, rho, dp) -> f64` | `Ṙ² = (2Δp/3ρ)((R₀/R)³−1)` | Rayleigh 1917; Brenner et al., RMP 74, 425 (2002) | integrates to τ = 0.91468 | plasmoid::cavitation | mainstream |
| `rayleigh_plesset(r, rdot, t, params) -> f64` (R̈) | `RR̈ + 1.5Ṙ² = (p_B − p_∞(t))/ρ − 4νṘ/R − 2S/(ρR)` | Brenner et al. 2002, sec. II (not opened; use Lauterborn & Kurz 2010) | empty-cavity limit gives τ = 0.91468 | plasmoid::cavitation | mainstream |
| `jet_speed(dp, rho, coeff) -> f64` | `v = k√(Δp/ρ)`, k ≈ 13 (touching), 17 (γ = 1.5) | Plesset & Chapman, JFM 47, 283 (1971) | 130 and 170 m/s at Δp/ρ = 10⁶ cm²/s² | plasmoid::cavitation | mainstream |
| `water_hammer(rho, c, dv) -> f64` | `p = ρcΔv` (Joukowsky) | Philipp & Lauterborn, JFM 361, 75 (1998) | 83 m/s gives about 0.1 GPa | plasmoid::cavitation | mainstream |
| `adiabatic_collapse_temperature(t0, ratio, gamma) -> f64` | `T = T₀(R₀/R)^{3(γ−1)}` | Brenner et al. 2002 (adiabatic picture) | 15,000 K from 300 K at R₀/R = 7.07 (Ar); compare 4,000–15,000 K, Flannigan & Suslick, Nature 434, 52 (2005) | plasmoid::cavitation | mainstream |
| `fusion_per_collapse(reaction, kt, n, volume, duration) -> f64` | `½n²⟨σv⟩Vτ` via `fusion::bosch_hale` or Gamow-peak integration below 0.2 keV | Bosch & Hale 1992 (in repo); Taleyarkhan 2002 vs Shapira & Saltmarsh 2002, Camara 2007 | 2.5×10⁻⁶² per collapse at 15,000 K for the generous hot spot; ⟨σv⟩ = 1.2×10⁻⁹³ m³/s | plasmoid::fusion (reuse `gamow_peak`, `bosch_hale_sigma`) | mainstream calc; bubble fusion refuted |
| `maxwell_reactivity_from_s(s_factor, b_g, mr, kt) -> f64` | `⟨σv⟩ = √(8/πm)(kT)^{-3/2}∫S e^{−E/kT−B_G/√E}dE` | Bosch & Hale 1992 | equals `bosch_hale(DDn, 1 keV)` to 1% | plasmoid::fusion | mainstream |
| `electrolysis_voltages() -> (f64, f64)` | `E_rev = −ΔG/2F`, `E_tn = −ΔH/2F` | NIST WebBook (CODATA, Cox et al. 1984) | 1.229 V, 1.481 V | plasmoid::conversion | mainstream |
| `electrolyser_hhv_efficiency(v_cell, faradaic) -> f64` | `η = 1.481 V · η_F / V_cell` | NIST; Ohmasa US 7,459,071 B2 (2.0 V, 1,000 L/h) | 0.741 at 2.0 V; 1,594 A for 1,000 L/h | plasmoid::conversion (Ohmasa as example) | mainstream; COP > 1 claims refuted by first law |
| `hhv_per_m3(x_h2, co2_fraction) -> f64` | `(1−f_CO₂)·x_H₂·ΔH/V_m` | NIST; Ohmasa US 2016/0145521 A1 | 12.75 MJ/m³ (H₂, STP); claim/bound = 6.6 (20% CO₂), 8.6 (40% CO₂) | plasmoid::conversion (example) | bound mainstream; claim unverified |
| `adiabatic_flame_temperature(fuel, oxidiser, p) -> f64` | Equilibrium of H₂O, H₂, O₂, OH, H, O with NASA-7 polynomials; enthalpy balance | Marzouk, ETASR 13, 11437 (2023); GRI-Mech 3.0 thermo | 3,076.92 K (Cantera), 3,074.51 K (CEARUN); 4,890.21 K without dissociation | plasmoid::conversion (or new combustion module) | mainstream |
| `condensation_shrink(t, p) -> f64` | `3V_m(T,p)/(2M_w/ρ_w)` | CODATA molar volume 22.413969 L/mol | 1,866 times at STP | plasmoid::conversion | mainstream ("implosion") |
| `cavitation_number(p, p_v, rho, v) -> f64` | `σ = (p − p_v)/(½ρv²)` | standard definition (no primary opened) | Ohmasa blade 1.57 m/s gives σ = 66 | plasmoid::cavitation | mainstream |
| `rankine_pressure_deficit(vmax, rho) -> f64` | `Δp = ρv_max²` (cyclostrophic, Rankine) | Wurman et al., BAMS 88, 31 (2007); Karstens et al., MWR 138, 2570 (2010) | 200 hPa at 132 m/s; 100 hPa ↔ 93 m/s | gravimetry::atmosphere | mainstream |
| `molecule_kinetic_energy(m_molar, v) -> f64` vs `WATER_DIMER_D0` | `½Mv²`; D₀ = 1105 cm⁻¹ | Rocher-Casterline et al., JCP 134, 211101 (2011) | D₀ = 13.2 kJ/mol; 132 m/s gives 1.19% of D₀ | metrology (constant) + gravimetry::atmosphere | mainstream; "H-bonds giving up" refuted |
| `dust_devil_charge_density(e, r, volume) -> f64` | `n ≈ 4πε₀r²E/(e·V/2)` | Farrell et al., JGR 109, E03004 (2004) | ≈10⁶ el/cm³ for E = 4.3 kV/m, r = 50 m | magnetics::space_weather or plasmoid::discharge | mainstream |
| `helmholtz_frequency(c, area, volume, neck, radius) -> f64` | `f = (c/2π)√(A/(V(L+1.5a)))` | Monteiro et al., arXiv:1805.04014 | c = 344 m/s at 22.3 °C; 223 Hz for 435 cm³ | gravimetry::atmosphere (beside `sound_speed`) | mainstream |
| `clamped_circular_plate_lambda2(n, s) -> f64` | roots of `J_n I_{n+1} + I_n J_{n+1} = 0`; `ω = λ²/a²·√(D/ρh)` | Leissa, NASA SP-160 (1969), Table 2.1 | 10.2158, 21.26, 34.88, 39.771 | physica::analogies (or new acoustics module) | mainstream |
| `kladov_pressure_window(sigma, p3) -> (f64, f64)` | `P₁ = (0.3–0.7)(P₂+P₃)`, `P₂+P₃−P₁ = (1–10)σ` | RU 2054604 C1 | P₂ = 0.5 MPa, P₁ = 0.15 MPa, k = 1.08 | example only (rigs are examples), plasmoid::cavitation | claim unverified |
