# Applied magnetics and MHD engineering: formulas, published test cases, records, and claim checks

Units are SI throughout. Conversions used: 1 G = 10⁻⁴ T, 1 MG = 100 T, 1 psi = 6894.757 Pa, 1 kn = 0.514444 m/s, 1 lb = 0.45359237 kg, 1 ft = 0.3048 m, 1 Wh = 3600 J, 1 mS/cm = 0.1 S/m, μ₀ = 4π×10⁻⁷ H/m (the 2019 SI value differs by 5×10⁻¹⁰ relative). Every test value tagged **[check §N]** is recomputed in `applied_magnetics_checks.py` in this folder. Run it from another directory: `cd /tmp && python3 -I <path>/applied_magnetics_checks.py`. Run on 2026-10-24: 49 published-value checks pass. Lines marked ILLUSTRATION in that script use assumed inputs and are not tests.

## 1. Magnet design: solenoid/Helmholtz/Bitter/hybrid formulas, magnetic pressure and hoop stress, stored energy, record fields, CFS/MIT TFMC, REBCO Jc scaling

### Takeaway
The design formulas are textbook and exact: on-axis solenoid and Helmholtz fields, the Fabry factor for resistive coils, magnetic pressure B²/2μ₀ and the BJR hoop stress. Published devices give several clean test cases. The best is the TFMC Lorentz line load, where I·B = 40.5 kA × 20.1 T = 814 kN/m against a published 815 kN/m. Records as of 2026-10: 45.22 T steady hybrid (Hefei, 2022); 42.02 T resistive (Hefei, 2024); 35.6 T all-superconducting user magnet (SECUF, Jan 2026); 48.7 T from an HTS test coil inside a resistive magnet (MagLab, Sept 2025); 100.75 T non-destructive pulsed (LANL, 2012). The "45.5 T" figure is a test coil in a resistive background magnet. It is not an all-superconducting record.

### Cited Findings
**Records and design numbers [mainstream]**
- Hahn et al., Nature 570, 496 (2019): 45.5 T DC from a REBCO test coil giving 14.4 T inside a 31.1 T resistive magnet. The tape had a 30 µm substrate, the winding ran at about 1260 A/mm² (1.26×10⁹ A/m²), and the coil used no-insulation winding — [IDEAS/Nature record](https://ideas.repec.org/a/nat/nature/v570y2019i7762d10.1038_s41586-019-1293-1.html). The coil quenched after the brief record run, and teardown showed edge damage to the REBCO from high stress — [Physics World](https://physicsworld.com/a/little-big-coil-creates-record-breaking-continuous-magnetic-field/). **[check §1]**
- MagLab "Little Big Coil" LBC9, announced 11 Sept 2025: 17.6 T inside the 31 T resistive magnet, 48.7 T total. The coil is about 2 in (51 mm) tall and 1.5 in (38 mm) across, wound from about 720 ft (219 m) of REBCO with no insulation — [MagLab news](https://nationalmaglab.org/news-events/news/a-prototype-miniature-superconducting-magnet/). It survived stress above 150,000 psi (1.03 GPa) and is described as a pilot for superconducting magnets above 50 T — [MagLab ASC highlight](https://nationalmaglab.org/magnet-development/applied-superconductivity-center/research/science-highlights/path-to-superconducting-coils/). **[check §1]**
- All-superconducting magnets. MagLab 32 T reached field on 8 Dec 2017, using REBCO inner coils (SuperPower) inside Nb₃Sn/NbTi outer coils (Oxford Instruments). It opened to users in 2021 — [MagLab](https://nationalmaglab.org/news-events/news/new-world-record-magnet-fulfills-superconducting-promise). A 32.35 T record followed in 2019, then 35.1 T from the Hefei Institute of Plasma Physics in Sept 2025 — [CAS Hefei](https://english.hf.cas.cn/nr/bth/202509/t20250929_1059018.html); [Cold Facts](https://cold-facts.org/2026/01/06/china-sets-world-record-with-35-1-tesla-superconducting-magnet-advancing-fusion-and-energy-research/). SECUF (Beijing) reported 35.6 T with a 35 mm usable aperture on 24 Jan 2026 — [IOP CAS](https://english.iop.cas.cn/news/202604/t20260408_1155395.html). Some coverage gives 35.5 T for the Hefei magnet, so the sources conflict.
- Steady hybrid. Hefei SHMFF reached 45.22 T on 12 Aug 2022 with 26.9 MW input: 34.22 T from the resistive insert plus about 11 T from the superconducting outsert, in a 32 mm bore. It was first built in 2016 as a 40 T magnet — [CAS](https://english.cas.cn/head/202208/t20220815_311432.shtml); [phys.org](https://phys.org/news/2022-08-china-world-strongest-steady-magnetic.html). It beat the MagLab 45 T hybrid of 1999 — [Silicon Republic](https://siliconrepublic.com/machines/china-world-record-strong-stable-magnetic-field). **[check §1]**
- MagLab 45 T hybrid: the magnet page lists an 11.5 T superconducting outsert, a 33.5 T resistive insert and 30 MW — [MagLab 45 T page](https://nationalmaglab.org/about-the-maglab/around-the-lab/meet-the-magnets/meet-the-45-tesla-hybrid-magnet/). The lab's Q&A says 33 MW — [MagLab Q&A](https://nationalmaglab.org/about-the-maglab/facts-figures/maglab-q-a/). The original design split was 31 T resistive plus at least 14 T outsert, and the outsert was degraded by an "unprotected" quench in July 2000 — [IEEE TAS 2003](https://ieeexplore.ieee.org/document/1211857).
- Resistive (Bitter-type): Hefei reached 42.02 T at 32.3 MW on 22 Sept 2024, beating MagLab's 41.4 T of 2017 — [CAS Hefei](https://english.hf.cas.cn/nr/headline/202409/t20240922_690291.html).

**CFS/MIT SPARC Toroidal Field Model Coil (TFMC), test 5 Sept 2021 [mainstream]**
- The coil reached 20.1 T peak on conductor at 40.5 kA, with 20 K joints. It carried 815 kN/m of Lorentz load on the REBCO stacks, with almost 1 GPa of stress taken by the case. It has 15 demountable joints at 0.5–2.0 nΩ in fields up to 12 T, uses 600 W of cooling at 20 K with up to 70 g/s flow at 20 bar-a, and its leads are rated to 50 kA. An intentional open-circuit quench was run at 31.5 kA. The coil is about 3 m across and the programme ran 2018–2021 — [Hartwig et al., arXiv 2308.12301 / IEEE TAS 2024](https://arxiv.org/abs/2308.12301). **[check §3]**
- Mass 9,270 kg, more than 270 km of REBCO, designed for up to 110 MJ stored energy — [CERN Indico abstract](https://indico.cern.ch/event/975584/contributions/4430919). CFS reports 20 T and 110 MJ in the 2021 test — [CFS](https://cfs.energy/news-and-media/cfs-mit-high-field-magnet-technology/).

**REBCO critical-current scaling [mainstream]**
- Senatore et al., SuST 29, 014002 (2016), covering six manufacturers at 4.2–77 K, up to 19 T and 0°/45°/90°: Jc(T,B) = Jc(0,0)·exp(−T/T*)·B^(−α). The power law holds from 0.5 to 19 T at T ≤ 30 K, but only from 0.1 to 1 T at 77 K. The overall relation is valid from about 0.1 T up to about 60% of B_irr. At 0.5 T, T* lies between 20 and 30 K. α(0°) is nearly constant from 5 to 40 K, at about 0.55 for SuperOx and BHTS and about 0.75 for SuperPower — [arXiv 1512.01930](https://arxiv.org/abs/1512.01930). A change ΔT multiplies Jc by exp(−ΔT/T*) whatever the starting temperature (same paper). **[check §5]**
- CERN (Succi, Ballarino et al., IEEE TAS 2025) validated a field-and-temperature scaling law on commercial tapes. They combined 4.2 K transport with magnetisation data up to Tc to find a single-vortex to collective-pinning crossover field B₀ — [Southampton ePrints](https://eprints.soton.ac.uk/500774).

**Formulas and their sources**
- Thin-solenoid on-axis field: B(z) = (μ₀nI/2)[(z+L/2)/√((z+L/2)²+R²) − (z−L/2)/√((z−L/2)²+R²)]. This comes from integrating the loop field B = μ₀IR²/2(R²+z²)^{3/2} (standard). The script checks the long-solenoid limit μ₀nI and the end field μ₀nI/2. Helmholtz pair: B₀ = (4/5)^{3/2}μ₀NI/R, with d²B/dz² = 0 at the centre **[check §4]**.
- Bitter/Fabry: Mulhall (1972) uses the "well known Fabry formula", with geometry factor G and a space factor λ — [Fermilab conf. paper](https://lss.fnal.gov/conf/C720919/p360.pdf). Derived here for uniform J (inner radius a₁, α = a₂/a₁, β = half-length/a₁): B₀ = μ₀·G·√(Pλ/(ρa₁)), with G = √(β/(2π(α²−1)))·ln[(α+√(α²+β²))/(1+√(1+β²))]. A numerical maximum is G_SI = 0.1426 at α ≈ 3.1, β ≈ 1.9, which is 0.179 in Montgomery's practical units (G·4π/10). That matches the literature value I remember **[check §4]**.

### Inferences
- **Rust calculators (inputs → published output):** `lorentz_line_load(I, B)` with 40.5 kA and 20.1 T gives 814 kN/m (published 815). `magnetic_pressure(B)` gives 806 MPa at 45 T and 3.98 GPa at 100 T. `rebco_jc(T, B, Jc00, T_star, alpha)`. `fabry_field(P, lambda, rho, a1, alpha, beta)`. `helmholtz_b0`. `solenoid_bz`. A Senatore fit cannot give absolute Jc without one transport point, so the crate should take Jc at a reference (T, B) and scale from it.
- Taking 110 MJ at 40.5 kA gives an inductance L = 2E/I² ≈ 0.134 H. That is my number, not a published one. The TFMC stores 11.9 kJ/kg, about 9% of the virial limit for steel at 1 GPa (127 kJ/kg) **[check §3]**.
- Magnetic pressure at 48.7 T is 944 MPa, close to the >1.03 GPa that LBC9 survived. That is why every gain beyond about 45 T now depends on tape strength and screening-current stress, not on Jc.
- An ideal uniform-J resistive coil with a₁ = 16 mm (32 mm bore), ρ = 2×10⁻⁸ Ω·m, λ = 0.8 and 26.9 MW gives an upper bound of about 46 T. Hefei's insert gives 34.22 T, which is about 74% of that ideal. That is plausible given non-uniform J and cooling channels. This is an illustration, not a test.
- Hybrid power is falling: 45 T at 30–33 MW (MagLab, 1999) against 45.22 T at 26.9 MW (Hefei, 2022).

### Gaps
- Published Fabry-factor tables (Montgomery 1969) were not fetched, so the 0.179 value rests on my derivation and memory.
- No published inductance or turn count for the TFMC turned up (the Whyte et al. IEEE TAS 2024 assessment was not read). Whether 110 MJ was measured or is the design value is unclear.
- Absolute REBCO Jc at given (T, B) for specific tapes is only in the paper's figures. I found no tabulated test point. The Fleiter–Ballarino CERN parameterisation is internal (EDMS) and was not found.
- No primary MagLab specification sheet gives the as-built resistive-insert power of the 45 T hybrid (30 MW and 33 MW are both quoted).

## 2. Superconducting magnetic energy storage: virial limit and published systems

### Takeaway
The virial theorem limits any magnet's stored energy by the stress its structure carries: E ≤ σ·V_structure, so M_structure ≥ ρE/σ [mainstream]. At realistic stresses that caps SMES at roughly 10–50 kJ/kg (3–14 Wh/kg). Lithium-ion cells reach about 150–250 Wh/kg (my comparison, not sourced). This is why SMES survived only as short, high-power grid stabilisers, such as BPA's 30 MJ / 10 MW unit of 1983.

### Cited Findings
- Tixador (IEEE CSC, 2008): the virial theorem ties the minimum structural mass to the stored energy. At 100 MPa working stress, SMES energy density is limited mechanically to about 10 kJ/kg. Force-free configurations exist only for infinite structures — [Tixador review slides](https://snf.ieeecsc.org/files/ieeecsc/slides/CR5_Final3_012008.pdf). **[check §6]**
- At 400 MPa and 8000 kg/m³, conventional solenoid or toroidal SMES reaches about 25–50 kJ/kg, and force-balanced coils about 50–100 kJ/kg in theory — [CERN Indico 2025 abstract](https://indico.cern.ch/event/1431972/contributions/6420406).
- Bullard and Haugan estimate about 30 Wh/kg (108 kJ/kg) for HTS SMES without extra structure, and perhaps about 150 Wh/kg with light, strong reinforcement, still under the virial limit — [CERN Indico 2017](https://indico.cern.ch/event/578092/contributions/2537694).
- Force-reduced magnets are argued to save structural mass, contrary to earlier claims — [Virginia Tech thesis](https://vtechworks.lib.vt.edu/items/620cc411-aace-496c-a73d-36743ecec6b7).
- BPA Tacoma SMES: a 30 MJ (8.33 kWh) coil with a 10 MW converter, installed 1982–83 to damp the 0.35 Hz Pacific AC Intertie oscillation. It exchanged up to 11 MW at 0.1–1.0 Hz on a 13.8 kV bus — [Rogers et al. abstracts](https://typeset.io/authors/j-rogers-gyvx79zmlv); [Springer chapter](https://link.springer.com/doi/10.1007/978-1-4613-0639-9_65). **[check §6]**

### Inferences
- **Rust calculator:** `virial_min_mass(E, sigma, rho, q=1.0)` = q·ρE/σ. Test: σ/ρ = 100 MPa / 8000 = 12.5 kJ/kg, against Tixador's "about 10 kJ/kg", which implies q ≈ 1.25 for a practical solenoid. q = 1 is the ideal bound for structure purely in tension, and solenoids need extra compression structure, so q is greater than 1.
- 30 Wh/kg needs about 864 MPa of average structural stress at steel density. That is plausible only if the REBCO's Hastelloy substrate (σ above 1 GPa) is counted as structure **[check §6]**.
- At full converter power the BPA coil holds 30 MJ / 10 MW = 3 s. SMES competes on power (MW in milliseconds), not energy.

### Gaps
- Primary data for later SMES units (Japanese Chubu/ISTEC, the American Superconductor D-SMES, Korean and Chinese HTS SMES) were not retrieved. Nor was the formal virial derivation (Moon 1982 or Eyssa) with the exact solenoid and toroid q factors.

## 3. Magnetic levitation: Earnshaw and its loopholes, the levitated frog, Inductrack, EDS and EMS maglev

### Takeaway
Earnshaw's theorem forbids stable static levitation of permanent magnets and paramagnets. It does not apply to induced diamagnetism, so a frog floats where B·dB/dz = μ₀ρg/|χ| = 1400.9 T²/m [mainstream]. Berry and Geim's paper gives a complete numeric test case that reproduces to 0.3%. Moving-conductor EDS (Inductrack, SCMaglev) gets stability from induced currents. Lift/drag rises with speed as kvL/R. Both of Post's published Inductrack numbers imply the same transition speed of 1.67 km/h.

### Cited Findings
- Berry and Geim, Eur. J. Phys. 18, 307 (1997): Earnshaw's theorem does not apply to induced magnetism, so a total-energy minimum (gravity plus magnetism) can exist for diamagnets. Paramagnets (χ ≈ 10⁻³) can be balanced but are laterally unstable. The equilibrium condition gives B·B′ = −1400.9 T²/m for water (χ = −8.8×10⁻⁶, ρ = 1000 kg/m³). In the Nijmegen Bitter magnet (4 MW dissipated), the inflection point is at z_i = 78 mm with B = 0.63 B₀ and gradient −8.15 B₀ T/m, so B₀ is predicted at 16.5 T. The predicted stable zone is 67.5–87.5 mm, and a 3 mm plastic sphere was stable from 69±1 to 86±1 mm. The induced current in a 10 mm frog is about 1.5 A, with an internal field of about 10⁻⁵ B₀ ≈ 1.5 G. Frogs were held for up to 30 min with no noticeable biological effect — [Berry & Geim 1997 PDF](https://obelix.physik.uni-bielefeld.de/~schnack/Herbstakademie/CD-2008/frog-ejp.pdf). **[check §7]**
- Inductrack (Post, LLNL; patent US6629503): Lift/Drag = kv(L/R) with k = 2π/λ. NdFeB Halbach arrays give about 1.0 T peak, and lift can exceed 40 t/m² of magnet area with magnet weight a few percent of the load. L/D is about 300 at 500 km/h. L/D = 30 at 50 km/h means about 6500 N drag and more than 90 kW drag power for a 20,000 kg urban car. The transition speed is "a few meters/second" — [Google Patents US6629503](https://patents.google.com/patent/US6629503). OSTI says "a few km/h", with auxiliary wheels needed at rest — [OSTI 791522](https://www.osti.gov/biblio/791522). Wikipedia gives a test model levitating above 35 km/h (22 mph), L/D 200:1 at 500 km/h and lift of 50× magnet weight — [Wikipedia Inductrack](https://en.wikipedia.org/wiki/Inductrack). **[check §8]**
- SCMaglev (superconducting EDS): levitates about 100 mm once above 150 km/h. The crewed record is 603 km/h (L0, 21 Apr 2015). Earlier records were 590 km/h (16 Apr 2015), 581 km/h (2 Dec 2003), 552 km/h (1999) and 531 km/h (1997). The relative passing speed record is 1026 km/h (2004) — [Wikipedia SCMaglev](https://en.wikipedia.org/wiki/SCMaglev).
- Shanghai Transrapid (EMS, feedback-stabilised): 431 km/h in service, the highest operating passenger speed from 2002 to 2021. It covers 30.5 km in 7 min 20 s — [Wikipedia Maglev](https://en.wikipedia.org/wiki/Maglev). A 501 km/h test run on 12 Nov 2003 is reported only by weak sources.

### Inferences
- **Rust calculators:** `diamag_levitation_bgradb(rho, chi, g)` gives 1400.9 T²/m. `levitation_center_field(bgradb, b_ratio, grad_coeff)` gives 16.5 T. `inductrack_lift_drag(v, v_t)` = v/v_t.
- **Inductrack consistency [check §8]:** both of Post's examples (300 at 500 km/h and 30 at 50 km/h) give v_t = R/(kL) = 1.67 km/h (0.46 m/s). That agrees with Post's "1–2 mph" estimate on Wikipedia, but not with the patent's "a few m/s". The patent text likely means the half-lift point loosely. A crate should take v_t as an input. Its examples also reproduce the 20 t lift (6500 × 30 / g = 19.9 t) and 90.3 kW drag power.
- Loopholes, by mechanism:
  - Diamagnetism: χ < 0, with superconductors as the extreme case χ = −1.
  - Spin stabilisation: the Levitron, referenced by Berry and Geim.
  - Active feedback: EMS, as in Transrapid.
  - Induced eddy currents in a moving conductor: EDS, as in SCMaglev and Inductrack. These need a minimum speed or wheels.

### Gaps
- Primary Transrapid EMS gap (about 10 mm) and SCMaglev field and drag data were not retrieved. The Halbach-array formula B₀ = B_r(1−e^{−kd})·sin(π/M)/(π/M) is standard (Halbach 1980, as used by Post) but was not in any fetched page. General Atomics test-track numbers were not found.

## 4. MHD generators: power density, efficiency, test facilities, why they stalled

### Takeaway
A segmented Faraday generator delivers P = σu²B²k(1−k) per unit volume at local electrical efficiency k. The power peaks at k = ½ [mainstream; derived from J = σ(E + u×B) with E = −kuB]. Large plants did deliver tens of MW: Avco Mk V 32 MW (1965) and U-25 20.4 MW (1970s). But enthalpy extraction stayed around 6–8%, and the work stalled on electrode corrosion, electrothermal instability, seed recovery and cost. Closed-cycle disk generators reached 15.7–18% in Fuji-1 blow-down runs and about 30% only in shock tubes.

### Cited Findings
- Avco Mk V, 1965: 32 MW in a one-minute run on alcohol with oxygen at 45 kg/s — [encyclopedia.com](https://www.encyclopedia.com/environment/encyclopedias-almanacs-transcripts-and-maps/magnetohydrodynamics). Another source gives 550 MW thermal for 32 MWe — [AGIR](https://www.agir.ro/univers-ingineresc/numar-5-2007/conversia-mhd-a-energiei-si-cercetari-rom%C3%82nesti-%C3%8En-domeniu_1702.html). Wikipedia says 35 MW gross with about 8 MW used to drive its magnet — [Wikipedia MHD generator](https://en.wikipedia.org/wiki/MHD_generator). The primary Avco/USAF report is [DTIC AD0817747](https://apps.dtic.mil/sti/pdfs/AD0817747.pdf). **[check §9]**
- U-25 (Moscow): a 250 MW gas combustor and 20 MW in 1972 (encyclopedia.com, above). Full rated power of 20.4 MW is reported in a 1979 review — [Natural Resources Forum 1979](https://ideas.repec.org/a/wly/natres/v3y1979i2p133-145.html). In 1974 RAND reported 5 MW of a planned 25 MW — [RAND R1404](https://www.rand.org/pubs/reports/R1404.html). Wikipedia: 25 MW design, 6 MW delivered by 1974. U-25B used a superconducting magnet and gave 1.4 MW. The US DOE CDIF ran 4000 h at a 50 MW demonstration (1992) — [Wikipedia](https://en.wikipedia.org/wiki/MHD_generator). **[check §9]**
- Japanese closed-cycle disk generators (Tokyo Tech):
  - Fuji-1 blow-down: 18% enthalpy extraction with Cs-seeded Ar — [OSTI 381435](https://www.osti.gov/biblio/381435). 15.7% at 0.46 MPa, and a maximum output of 516.7 kW at 0.6 MPa — [OSTI 442748](https://www.osti.gov/biblio/442748).
  - Shock tube: 30.2% in a new shock-tube-driven disk, and 38% with Cs-seeded He and inlet swirl — [OSTI 381436](https://www.osti.gov/biblio/381436); [OSTI 207708](https://www.osti.gov/biblio/207708).
  - Wikipedia credits Fuji-1 itself with 30.2% and about 100 MW/m³, and says Tokyo Tech held a 22% efficiency record in 1994 — [Wikipedia](https://en.wikipedia.org/wiki/MHD_generator). The OSTI sources conflict with this.
  - JSME (2015): extraction reached the 30% target, but isentropic efficiency stayed low — [JSME newsletter](https://jsme-fed.org/newsletters-e/2015_12/no7.html).
- Why it stalled — [Wikipedia](https://en.wikipedia.org/wiki/MHD_generator):
  - Electrothermal instability "crippled" non-equilibrium (closed-cycle) programmes.
  - Oxygen migrating through slag corrodes metallic anodes.
  - Seed (potassium) recovery adds cost.
  - The plants needed exotic ceramics and never had long-duration tests.
  - Other technologies matched the efficiency at lower lifecycle cost.
  - Open-cycle coal channels ran near 17%, against about 40% for Rankine plants. "Up to 60%" was a combined-cycle projection.

### Inferences
- **Rust calculators:** `faraday_mhd_power_density(sigma, u, B, k)`, `mhd_load_factor_efficiency(k) = k` and `mhd_enthalpy_extraction(P_e, P_th)`. Test cases: Mk V 32/550 gives 5.8%, U-25 20.4/250 gives 8.2% **[check §9]**. With assumed σ = 10 S/m, u = 1000 m/s, B = 5 T and k = ½, the power density is 62.5 MW/m³, the same order as the "about 100 MW/m³" claimed for Fuji-1 (illustration).
- Continuous-electrode and Hall generators need the Hall parameter β = ω_ce·τ_e. In standard texts (Rosa) the continuous-electrode Faraday power falls by about 1/(1+β²), and the Hall generator recovers about β²/(1+β²) of the segmented case. These forms are not verified against a fetched source, so treat them as a to-do.
- Honest scorecard: MW-scale open-cycle plants extracted about 6–8% (on nameplate). Mk V's own magnet took about 25% of gross output, and copper magnets set the economics. Superconducting magnets (U-25B) arrived only at the 1.4 MW scale.

### Gaps
- No fetched primary source states P = σu²B²k(1−k) or the Hall forms. The Rosa text was not online. The formula here is derived from Ohm's law, which is exact for an ideal segmented channel.
- The Mk V thermal input (400 vs 550 MW) and run length (one vs three minutes) conflict between sources. Fuji-1 power density was not confirmed from a primary source.

## 5. MHD propulsion (Yamato-1), liquid-metal pumps, electromagnetic flowmeters

### Takeaway
Seawater conducts only 4.29 S/m, about 10⁷ times less than copper, so a seawater MHD thruster wastes most of its electrical input as Joule heat. Its electrical efficiency is η_e = σuB²/(σuB² + f), where f = JB is the thrust density, and only fields of 10–20 T or more help. Yamato-1 (1992, 4 T, 2 × 8 kN design) reached about 6–8 kn at about 15% system efficiency [mainstream, low-performing]. Electromagnetic flowmeters use the same physics in reverse: U = B·D·v.

### Cited Findings
- Yamato-1: 185 t displacement, 30 m long, 10.39 m beam, 2.69 m draft. It has two 4 T MHD thrusters of 8 kN each. Speed is 8 kn (15 km/h, 4.1 m/s) in the infobox, with no citation. The overall system efficiency of 15% comes from Sasakawa (1997). The ship was completed in 1991 and first ran in Kobe harbour in June 1992 — [Wikipedia Yamato-1](https://en.wikipedia.org/wiki/Yamato-1). One sea trial logged 6.2 kn (3.2 m/s) — [Seattle Times, 20 Jul 1992](https://archive.seattletimes.com/archive/19920720/1503104/engineering----may-the-force-be-with-you). Each unit gave about 8000 N of Lorentz force at 4 T and 4000 A — [AMSM review (Arak Univ.)](https://masm.araku.ac.ir/article_252207_f2dac30dc6d60c7986558e2f976a8e57.pdf?lang=en). Sasakawa's account says the magnets were the main weight problem (more than 100 t at first, later about 30 t) — [Ship & Ocean Foundation PDF](https://www.spf.org/en/_opri_media/publication/docs/yamato-1.pdf) (from a search summary; the PDF was too large to fetch). **[check §10]**
- Parametric studies:
  - Speed performance rises in proportion to seawater conductivity and with B², and efficiency favours larger systems with longer channels — [DTIC ADA218318](https://apps.dtic.mil/sti/tr/pdf/ADA218318.pdf) (attribution from a search summary).
  - Thruster efficiency rises with field strength and thruster diameter, and falls with wall roughness and flow velocity — [TRID 443209](https://trid.trb.org/View/443209).
  - A 2023 design baseline uses a 6.5 T two-pole superconducting magnet with a 260 mm warm bore — [Magnetohydrodynamics 59(2) 2023](https://mhd.sal.lv/contents/2023/2/MG.59.2.7.R.html).
- Seawater conductivity reference: C(35, 15 °C, 0 dbar) = 42.914 mS/cm = 4.2914 S/m (TEOS-10 convention; PSS-78 itself does not fix the value) — [TEOS-10 gsw_C_from_SP](https://mail.teos-10.org/pubs/gsw/html/gsw_C_from_SP.html). Poisson used 42.933 mS/cm, and EURAMET quotes 42.896 mS/cm — [Sea-Bird AN14](https://pscfiles.apl.uw.edu/woodgate/BeringStraitArchive/BeringStraitMooringData/BeringStraitMoorings2007to2009IPY_versionMar10/SeasoftForWavesProcessingSoftware/website/application_notes/AN14.htm); [EURAMET](https://www.euramet.org/technical-committees/tc-projects/details/project/study-on-the-traceability-of-salinity-measurements-in-seawater). **[check §10]**
- Electromagnetic flowmeter: U = K·B·D·v, where D is the electrode spacing (the pipe inner diameter), so Q = πUD/(4KB) — [Huddersfield thesis](https://eprints.hud.ac.uk/6884/1/25_T_Leeungculsatien.pdf). The field is constantly reversed to cancel electrode electrochemistry, so permanent magnets cannot be used. The fluid must be ionically conducting and the pipe insulated (rubber-lined) — [Wikipedia Magnetic flow meter](https://en.wikipedia.org/wiki/Magnetic_flow_meter).

### Inferences
- **Rust calculators:** `mhd_thrust(I, B, d) = I·B·d`. For Yamato-1, F/(IB) = 8000/(4000·4) = 0.5 m effective electrode length. `mhd_thruster_electrical_efficiency(sigma, u, B, f)` and `em_flowmeter_voltage(B, D, v)`, which gives 0.5 mV at 5 mT, 0.1 m and 1 m/s (illustration).
- Illustration at σ = 4.29 S/m, u = 5 m/s and f = 2×10⁴ N/m³: η_e = 1.7% at 4 T, 9.7% at 10 T and 30% at 20 T **[check §10]**. Joule loss scales as J²/σ while useful work scales as JuB. The fix is B², which explains why designers moved to 6.5 T or more. Seawater MHD is [mainstream physics, impractical engineering]. The "silent MHD submarine" remains fiction-grade at present fields.
- Liquid-metal pumps work far better because σ is about 10⁶–10⁷ S/m (sodium, NaK, lithium), about 10⁶ times seawater. This is my inference, not sourced here.

### Gaps
- Yamato-1 electrode voltage, input power, duct geometry and the measured split between electrical and jet efficiency are in Motora/Takezawa et al. and Sasakawa (1997). The large PDF could not be fetched.
- No published EM liquid-metal pump test case (flow, head, efficiency) was retrieved.

## 6. Electromagnetic launchers, induction heating, eddy-current braking

### Takeaway
Railgun force is F = ½L′I², with L′ ≈ 0.4–0.7 µH/m for single-turn rails. The US Navy fired a world-record 33 MJ shot on 10 Dec 2010 at about 2.5 km/s. The programme was halted in 2021 after about US$500 M, on barrel life (12–24 rounds) and integration [mainstream physics; programme failed on engineering]. Induction heating and eddy braking follow skin depth δ = √(2ρ/(ωμ)), which gives 8.42 mm for copper at 60 Hz.

### Cited Findings
- ONR shot of 33 MJ at NSWC Dahlgren, 10 Dec 2010, a world record — [ONR release](https://www.onr.navy.mil/media-center/press-releases/2010/electromagnetic-railgun-new-record) (page refused connection when fetched; the record is confirmed by [Guinness](https://www.guinnessworldrecords.com/world-records/most-powerful-electromagnetic-rail-gun) and [Military.com](https://www.military.com/defensetech/2010/12/13/bang)). Military.com calls it nearly triple the 10.64 MJ of Jan 2008, and ONR claims at least 110 nmi (204 km) range.
- NavWeaps figures — [NavWeaps](https://navweaps.com/Weapons/WNUS_Rail_Gun.php) **[check §11]**:
  - Round: 15 kg (33 lb) projectile, 20 kg (44 lb) saboted round.
  - Speed and energy: about 2,500 m/s (8,200 fps) nominal; saboted muzzle energy about 63 MJ against about 150 MJ breech energy.
  - Rate and power: 6–12 rounds/min; about 30 MW needed at 12 rounds/min.
  - Barrel life: 12–24 rounds (prototype).
  - Test history: lab firing began 31 Jan 2008 at 10.64 MJ with a 3.2 kg (7 lb) slug at 2,520 m/s (8,268 fps).
  - Programme end: cancelled July 2021.
- BAE's 32 MJ prototype was delivered to Dahlgren in January 2012. Its range is quoted as up to 220 mi (354 km) for the BAE prototype and 100 nmi for a 32 MJ system. Forces wear out the barrel after one or two dozen shots — [Wikipedia EMRG](https://en.wikipedia.org/wiki/Electro-Magnetic_Laboratory_Rail_Gun).
- The FY2021 request was US$9.5 M. The FY2022 request zeroed the line items, and the Navy "paused" the work at the end of 2021. Total spending was about US$500 M — [Task & Purpose](https://taskandpurpose.com/news/navy-electromagnetic-railgun-dead/); [Military & Aerospace Electronics](https://militaryaerospace.com/sensors/article/14204856/electromagnetic-railgun-budget-request-hypersonic).
- Single-turn railguns have dL/dz between 0.4 and 0.7 µH/m. N-turn designs scale as about N², and six-turn devices reach about 40 µH/m — [US patent 4833965](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/4833965). L′ increases as the armature moves toward the muzzle — [IEEE-linked study](https://dspace.tbzmed.ac.ir/items/b498a851-d734-4439-bed9-ef47acdd2b08).
- Skin depth δ = √(2ρ/(ωμ)) ≈ 503√(ρ/(μ_r f)). Copper: 9.22 mm at 50 Hz, 8.42 mm at 60 Hz, 652 µm at 10 kHz, 65.2 µm at 1 MHz. The page text also says "about 8.5 mm at 60 Hz" — [Wikipedia Skin effect](https://en.wikipedia.org/wiki/Skin_effect). **[check §11]** The table matches σ = 5.96×10⁷ S/m (pure Cu at 20 °C). The 8.5 mm figure matches the article's σ = 58.5 MS/m.

### Inferences
- **Rust calculators:**
  - `railgun_force(Lp, I)`.
  - `muzzle_energy(m, v)`. Tests: 44 lb at 2500 m/s gives 62.4 MJ against "about 63" (1%). 7 lb at 8268 fps gives 10.08 MJ against a published 10.64 MJ. The 5% gap suggests the published figure includes armature or sabot mass or a higher velocity.
  - `skin_depth(rho, f, mu_r)`, tested at 8.42 mm, 9.22 mm and 65.2 µm.
- 33 MJ at 2.5 km/s implies about 10.6 kg launched. Spread over about 10 m of rail, the average force is 3.3 MN, which needs I ≈ 3.6 MA at L′ = 0.5 µH/m (illustration). Energy efficiency is about 63/150 = 42% breech to muzzle on NavWeaps figures.
- Coilguns (F = I₁I₂ dM/dx) have no comparable published record. The crate can implement the mutual-inductance force on top of the loop formulas from §1.

### Gaps
- No primary US Navy figure for peak current or rail length. No published coilgun test case with energy, mass and velocity was retrieved. No eddy-current brake test case (force against speed for a given conductor and field) was retrieved. The formula family exists (low-speed drag ∝ σv; above the magnetic Reynolds crossover, drag ∝ 1/v), but it is not sourced here.

## 7. Pulsed and explosive fields: flux compression, single-turn coils, the 100 T non-destructive record, Z/MagLIF

### Takeaway
Conserving flux gives B = λB₀(r₀/r)², where λ ≤ 1 is the fraction of flux retained. The records span three categories:
- Explosive flux compression: about 2800 T (VNIIEF, Sarov; 170 kg TNT). Sakharov's group reported "hundreds of tesla" by 1952 and about 2500 T for MK-1 in a secondary source.
- Indoor electromagnetic flux compression: 1200 T (ISSP Tokyo, 2018).
- Non-destructive multi-shot: 100.75 T (LANL, 22 Mar 2012). Single-turn coils reach 100–300 T in a semi-destructive mode.

Magnetic pressure rises from 4 GPa at 100 T to 3.1 TPa at 2800 T, which is why only destructive methods pass about 100 T.

### Cited Findings
- Sakharov's MK-1 (axial field compression) experiments started in spring 1952. They reached "millions of gauss (hundreds of teslas)" from a 3 T (30 kG) seed. MK-2 reached 100 MA (10⁸ A) and 200 T, converting up to 20% of the explosive energy to magnetic energy. A three-module disk generator 1 m across delivered 100 MJ at 256 MA, and VNIIEF systems use up to 25 modules — [Wikipedia EPFCG](https://en.wikipedia.org/wiki/Explosively_pumped_flux_compression_generator). The 1965 Doklady paper "Magnetic cumulation" (Sakharov, Lyudaev, Pavlovskii et al.) describes MK-1 and MK-2 — [mathnet.ru](https://mathnet.ru/eng/dan31766). Sakharov's 1966 review gives 15–25 MG (1500–2500 T) in individual US and USSR experiments — [UFN 1966](https://ufn.ru/ufn66/ufn66_2/ufn662i.pdf).
- 28 MG (2800 T), the highest explosively driven MFC field, reported by Boyko et al. (VNIIEF) — [Velikovich, APS DPP 2012](https://archive.aps.org/dpp/2012/pt3/1). Cascades reach 1000–1500 T routinely, and the 2800 T record used 170 kg TNT (stated as 680 MJ), destroying everything within several metres — [arXiv 2605.11379 review (2026)](https://arxiv.org/abs/2605.11379). **[check §12]**
- ISSP 1200 T: Nakamura, Ikeda, Sawabe, Matsuda and Takeyama, Rev. Sci. Instrum. (2018), doi:10.1063/1.5044557. Electromagnetic flux compression with a liner imploding at about 5 km/s, measured by Faraday rotation — [Sci.News](https://www.sci.news/physics/strongest-magnetic-field-achieved-indoors-06420.html); [ISSP](https://www.issp.u-tokyo.ac.jp/maincontents/highlight/2018/28_en.html); [Physics World](https://physicsworld.com/a/indoor-record-for-magnetic-field-strength-is-smashed-by-physicists-in-japan/).
- Single-turn coils: 100–300 T for microseconds. They are semi-destructive because the coil explodes outward and the sample survives — [LNCMI](https://lncmi.cnrs.fr/en/megagauss-fields/); [Gen et al. arXiv 2102.10254](https://arxiv.org/abs/2102.10254). ISSP horizontal STC: 300 T in a 5 mm bore, 200 T in 10 mm, about 6 µs, 200 kJ banks — [ISSP Activity Report 2024](https://www.issp.u-tokyo.ac.jp/maincontents/docs/ar2024_Part6.pdf).
- LANL 100.75 T on 22 Mar 2012, multi-shot (up to one pulse per hour). It uses seven coil sets of about 18,000 lb (8.2 t) and a 1200 MJ motor generator, and is described as about 2,000,000 times Earth's field — [DOE](https://www.energy.gov/articles/teslas-lab-los-alamos-powerful-magnets-come-full-circle); [MagLab 100 T page](https://nationalmaglab.org/about-the-maglab/around-the-lab/meet-the-magnets/meet-the-100-tesla-pulsed-magnet/). **[check §1]**
- Sandia Z / MagLIF (details in `../../plasmoids-and-fusion/notes/z_pinch_and_plasma_focus.md`):
  - Z drives up to 27 MA in about 100 ns.
  - Gomez 2014: 10 T preseed, 19 MA, up to 2×10¹² D–D neutrons — [Wikipedia MagLIF](https://en.wikipedia.org/wiki/Magnetized_liner_inertial_fusion).
  - Gomez 2020: 15.9 T, 20 MA, 1.1×10¹³ D–D neutrons.
  - Fuel magnetisation from secondary DT neutrons: more than 10¹⁰ secondary DT neutrons at a ρR of only 2 mg/cm² (Gomez et al. PRL 113, 155003) — [APS accepted ms](https://link.aps.org/accepted/10.1103/PhysRevLett.113.155003); method in [Schmit et al. PRL 113, 155004](https://link.aps.org/accepted/10.1103/PhysRevLett.113.155004).
  - BR ≈ 0.4 MG·cm (0.4 T·m) — [Sandia APS DPP abstract](https://meetings-archive.aps.org/dpp/2017/bp11/117/) (attribution from a search summary). It falls from about 0.5 to about 0.2 MG·cm as preheat rises (Lewis et al. 2021) — [OSTI 1820420](https://www.osti.gov/pages/servlets/purl/1820420).

### Inferences
- **Rust calculators:** `flux_compression_field(B0, r0, r, lambda)`. `magnetic_pressure(B)`: 4.04 GPa at 100.75 T, 35.8 GPa at 300 T, 573 GPa at 1200 T, 3.12 TPa at 2800 T. `tnt_energy(kg)`: 170 kg gives 711 MJ at 4.184 MJ/kg, against the published 680 MJ, which implies about 4.0 MJ/kg **[check §12]**.
- A 3 T seed compressed 25× in radius gives 1875 T ideally, so 2800 T needs either more convergence or a cascade (staged) design. Cascades are what VNIIEF used. This is illustration only.
- MagLIF: a 10 T seed in a liner of about 2.3 mm radius reaches BR = 0.4 T·m at r ≈ 130 µm if flux is conserved. That is B ≈ 3000 T at a convergence of about 17. This is illustration only, and liner radius, Nernst losses and λ are not sourced here.

### Gaps
- The exact date and venue of the 2800 T record (Boyko et al., Megagauss conference, around 1998–2001) and its sample bore were not confirmed. The MK-1 "25 MG in 1964" figure comes only from an encyclopedia mirror, and the date is unverified.
- No published λ (flux-loss) values for a specific generator were retrieved (e.g., Los Alamos Ranchero or Fowler's results). The script uses λ only as an input.

## 8. Magnetic space propulsion: Hall, MPD, VASIMR, magnetic nozzles, magnetic sails, M2P2, E-sails

### Takeaway
Every electric thruster obeys η = T·g₀·I_sp/(2P). Published numbers reproduce well:
- SPT-100B: 83 mN, 1.35 kW, 1600 s → 48%.
- VX-200: 5.7 N, 5000 s, 200 kW RF → 70% against a claimed 72%.
- AEPS: more than 580 mN at 12.5 kW and more than 2600 s → about 59%.

Hall thrusters are mature, flown hardware [mainstream]. The X3 holds the Hall thrust record of 5.4 N (2017). VASIMR's thruster numbers held up in short runs [mainstream], but its 39-day Mars claim needs reactors at about 1 kg/kW that do not exist [contested]. Magsail physics is sound, but the original drag was overestimated 3.1× [mainstream, corrected]. The small M2P2 "bubble" thrust claim did not survive kinetic analysis [refuted for the small device]. E-sail thrust is about 500 nN/m at 1 AU [mainstream theory; no flight demonstration found].

### Cited Findings
- **SPT-100:** 1,350 W, 1,600 s. The base model gives 80 mN at 47%, the SPT-100B 83 mN at 48% — [Wikipedia SPT-100](https://en.wikipedia.org/wiki/SPT-100). NASA gives about 1600 s at η = 0.50 — [NTRS 19940020707](https://ntrs.nasa.gov/citations/19940020707). **[check §13]**
- **NASA/UMich X3** nested Hall thruster: 5–102 kW on xenon; 5.4 ± 0.1 N at 99 kW and 400 V. Efficiency ran 0.54–0.67 and I_sp 1800–2650 s over the map — [NTRS 20190001137](https://ntrs.nasa.gov/citations/20190001137). The previous record was 3.3 N (university press claim) — [Space.com](https://www.space.com/38444-mars-thruster-design-breaks-records.html).
- **AEPS/HERMeS** 12.5 kW: more than 580 mN, more than 2600 s at 600 V — [NTRS 20190029263](https://ntrs.nasa.gov/api/citations/20190029263/downloads/20190029263.pdf). 13.3 kW total including the PPU, and a maximum of 600 mN — [Wikipedia AEPS](https://en.wikipedia.org/wiki/Advanced_Electric_Propulsion_System).
- **Self-field MPD (Maecker):** T = (μ₀I²/4π)[ln(r_a/r_c) + ¾] — [MDPI Energies 12, 2428 (2019)](https://www.mdpi.com/1996-1073/12/12/2428/htm). Choueiri (JPP 1998) found it too simple, especially at low current, where electrode current distribution matters — [Princeton](https://collaborate.princeton.edu/en/publications/scaling-of-thrust-in-self-field-magnetoplasmadynamic-thrusters/). Princeton argon data span 8–50 kA and 2–36 g/s — [NTRS 19720019032](https://ntrs.nasa.gov/api/citations/19720019032/downloads/19720019032.pdf).
- **VASIMR VX-200:**
  - 72 ± 6% thruster efficiency at 4900 ± 300 s and 200 kW of coupled RF on argon, with ion-cyclotron heating efficiency 85 ± 7% — [Longmier et al., JPP 2014](https://www.adastrarocket.com/technical-papers-archives/Ben-JPP-2014.pdf).
  - The 2011 papers give 72% at 5000 s and 5.7 N — [Berkeley SPRG seminar abstract](https://cdaweb.gsfc.nasa.gov/pub/data/fast/documents/archived_websites/sprg.ssl.berkeley.edu/seminars/seminars_042611.html).
  - Ad Astra (2015) quotes 200 kW electrical for 5 N (40 kW/N) — [Wikipedia VASIMR](https://en.wikipedia.org/wiki/VASIMR).
  - The VX-200SS ran 88 h continuously at 80 kW in July 2021, stopping 12 h short of the 100 h NASA NextSTEP goal because of a facility sensor fault. Thrust stayed within ±10% — [SpaceNews](https://spacenews.com/?p=84168); [AIAA ASCEND 2021](https://www.adastrarocket.com/wp-content/uploads/2021/12/2021_AIAA_ASCEND_VASIMR.pdf).
  - **Claim:** 39-day Earth–Mars transit at about 200 MW — [The Space Review](https://thespacereview.com/article/1690/1).
  - **Critique (Zubrin, 2011):**
    - Reactor mass: at α = 20 kg/kW, 200 MW of reactor weighs 4000 t against about 600 t assumed for the whole mission. That mass gives 6–8 months at best.
    - Flown reactor: Topaz gave 10 kW at about 100 kg/kW.
    - Efficiency: about 50% in short test-stand burns, against about 70% for ion engines with thousands of flight hours.
    — [Universe Today](https://www.universetoday.com/articles/zubrin-claims-vasimr-is-a-hoax). **[check §13]**
- **Magnetic sail:** Andrews and Zubrin (proposed 1988; papers 1990–91) used a 50–100 km superconducting loop of about 100 t. Drag F = C_d(ρu²/2)·πR_mp². Magnetopause pressure balance gives B_mp = u√(ρμ₀C_so), with C_so = 2, 1 or ½ depending on the reference. Freeland (2015, Project Icarus) found thrust optimistic by a factor of 3.1 from a numerical integration error — [Wikipedia Magnetic sail](https://en.wikipedia.org/wiki/Magnetic_sail); [NTRS 19910012840](https://ntrs.nasa.gov/citations/19910012840). Order of magnitude: F ~ πR_m²ρv² — [arXiv 2510.21743](https://arxiv.org/abs/2510.21743).
- **M2P2** (Winglee and Slough, 2000): a 15 km standoff (20 km cross-section) bubble was to reach 50–80 km/s after about 3 months. That was extrapolated from much smaller simulations — [JAXA/JSTS 2004 review](https://www.jstage.jst.go.jp/article/jsts/20/2/20_2_1/_article). Khazanov et al. (2003) found MHD invalid at that scale (below the ion Larmor radius and skin depth), with kinetic simulation giving negligible thrust — [NTRS 20030067378](https://ntrs.nasa.gov/citations/20030067378). The JAXA magnetoplasma sail (MPS) moved to 10–100 m coils for sub-newton thrust, and the action/reaction and finite-Larmor issues remain open — [JAXA repository](https://jaxa.repo.nii.ac.jp/records/5915).
- **E-sail (Janhunen):** about 0.5 mN/km (500 nN/m) at 1 AU for a 10 g/km tether — [arXiv 1404.7430](https://arxiv.org/abs/1404.7430). About 580 nN/m — [Copernicus preprint 2012](https://gi.copernicus.org/preprints/2/429/2012/gid-2-429-2012.pdf). 50–100 nN/m at 20 kV in the 2007 formulation — [Janhunen 2007](https://electric-sailing.fi/Firenze2007.pdf). Thrust scales as about 1/r (or r^(−7/6) in the 2007 paper).

### Inferences
- **Rust calculators:**
  - `thruster_efficiency(T, Isp, P)`. Tests: SPT-100B 48.2%, SPT-100 46.5%, VX-200 69.9% against 72%, AEPS 59%.
  - `thrust_from_power(eta, P, Isp)`. The 2014 72% at 4900 s implies 6.0 N, against 5.7 N measured in 2011. That is consistent within the ±6% and ±300 s quoted.
  - `maecker_thrust(I, ra_rc)`: 94 N at 20 kA with r_a/r_c = 5 (illustration).
  - `magsail_drag(C_d, rho, u, R)`: 26.6 N for R = 100 km with C_d = 1 at n = 5 cm⁻³ and 450 km/s, and 8.6 N after Freeland's correction (illustration).
  - `esail_thrust(L, r_AU)`: 1 N needs about 2000 km of tether (20 kg) at 1 AU.
  - `nuclear_electric_mass(P, alpha)`: 4000 t.
  - All **[check §13–14]**.
- X3 consistency: 5.4 N at 99 kW is 54.5 mN/kW. Efficiencies of 0.54–0.67 imply I_sp of 2019–2505 s, inside the published 1800–2650 s range. The record holds up internally.
- The VASIMR thruster claims held up. The Mars architecture is a statement about an imagined reactor (α ≈ 1 kg/kW), not about the thruster. The "72%" is RF-coupled thruster efficiency and leaves out the DC-to-RF conversion loss.

### Gaps
- No primary measured MPD thrust test point (thrust at a given current and electrode radii) was retrieved for a Maecker check. Magnetic-nozzle detachment efficiencies (VX-200 or Helicon) were not retrieved. No flight E-sail result was found (ESTCube-1's tether did not deploy, which is not confirmed here).

## 9. Medical and sensing: NMR/MRI Larmor frequency, SAR limits, SQUID and optically pumped magnetometers, magnetocaloric cooling

### Takeaway
Proton Larmor frequency is γ_p/2π = 42.577 478 461 MHz/T (CODATA 2022). The water-shielded value is 42.576 385 43 MHz/T. SAR limits in normal mode are 2 W/kg whole-body and 3.2 W/kg head (IEC 60601-2-33) [mainstream]. The Romalis SERF result (160 aT/√Hz in 0.45 cm³) reproduces its published 44ħ energy resolution to 1.2%. Gd gives about 3.6 K/T near room temperature.

### Cited Findings
- CODATA 2022: γ_p/2π = 42.577 478 461(18) MHz/T. Shielded proton (H₂O, sphere, 25 °C): 42.576 385 43(17) MHz/T, or γ′_p = 2.675 153 194(11)×10⁸ s⁻¹T⁻¹ — [NIST wall chart 2022](https://www.physics.nist.gov/cuu/pdf/wall_2022.pdf); [arXiv 2409.03787](https://arxiv.org/abs/2409.03787). **[check §15]**
- IEC 60601-2-33 (2nd and 3rd editions):
  - Normal mode: 2 W/kg whole-body, 3.2 W/kg head. Partial-body limit (10 − 8R) W/kg, where R is the exposed-to-total mass ratio.
  - First-level mode: 4 W/kg whole-body, derated 0.25 W/kg per °C at higher temperature.
  - Any 10 s window: up to 2× the stated values.
  — [CMRR IEC SAR guidelines](https://www.cmrr.umn.edu/media/forms/CMRR_IEC-SAR_Guidelines.pdf); [US 7860551](https://patents.google.com/patent/US7860551). The 4th edition (2022) and its 2023 amendment exist — [VDE preview](https://www.vde-verlag.de/iec-normen/preview-pdf/info_iec60601-2-33{ed4.0}b.pdf). **[check §15]**
- SERF magnetometers: Dang, Maloof and Romalis, APL 97, 151110 (2010): 160 aT/√Hz as a gradiometer, 0.45 cm³ volume, 44ħ energy resolution — [arXiv 0910.2206](https://arxiv.org/abs/0910.2206). Kominis et al., Nature 422, 596 (2003): 0.54 fT/√Hz in 0.3 cm³, with a fundamental limit below 0.01 fT/√Hz — [doi 10.1038/nature01484](https://www.doi.org/10.1038/NATURE01484). **[check §15]**
- SQUID and OPM in practice:
  - A SQUID MEG system reports about 5 fT/√Hz white noise — [arXiv 1310.4619](https://arxiv.org/abs/1310.4619).
  - QuSpin Gen-2 zero-field OPM: under 15 fT/√Hz (3–100 Hz), typically 7–10 fT/√Hz — [QuSpin](https://quspin.com/products-qzfm-gen2-arxiv/).
  - "1 fT/√Hz" for low-Tc SQUIDs is a best case — [Chinese Physics B 28, 040703 (2019)](https://cpb.iphy.ac.cn/article/2019/1983%5Ccpb_28_4_040703.html) (attribution from a search summary).
- Magnetocaloric: Gd gives about 3.6 K adiabatic ΔT at 1 T near room temperature, scaling about as B^0.7. Brown's regenerative (AMR) experiment reached a 46 K span with the hot end at 319 K and 7 T — [Risø DTU, arXiv 0902.0812](https://arxiv.org/abs/0902.0812).

### Inferences
- **Rust calculators:**
  - `larmor_frequency(B, shielded)`: 63.866 MHz at 1.5 T, 127.732 MHz at 3 T, 298.042 MHz at 7 T.
  - `sar_partial_body_limit(R, mode)`: 6 W/kg at R = 0.5 in normal mode.
  - `magnetometer_energy_resolution(dB, V) = dB²V/(2μ₀)`: 43.5ħ against the published 44ħ.
  - `gd_adiabatic_dT(B) = 3.6·B^0.7`: 14 K at 7 T (illustration; no source states this single-pass value).
- MRI's 2 W/kg whole-body limit is roughly the human basal metabolic rate (about 1.2–1.5 W/kg; my comparison). That is why the derating is tied to core temperature.

### Gaps
- The 4th-edition (2022) IEC limit values were not confirmed as unchanged. No primary SQUID white-noise benchmark (best low-Tc device) was retrieved. No primary magnetocaloric ΔT table (Pecharsky/Gschneidner Gd data) was retrieved.

## 10. Fringe and over-claimed devices: published claim against measured value

### Takeaway
Most of the over-claims in this field are scale-up extrapolations (VASIMR's 39 days, M2P2, open-cycle MHD "60%", railgun range), not physics violations. Steorn Orbo-type "magnet motors" are the outright refuted class. Several record headlines also mix categories, such as test coil against user magnet or shock tube against blow-down rig.

### Cited Findings
- **Steorn Orbo** (permanent-magnet "free energy") [refuted]: Steorn's own 22-member jury, chaired by Ian MacDonald, ruled unanimously in June 2009 that its "attempts to demonstrate the claim have not shown the production of energy." The 2007 London demonstration was cancelled for technical difficulties — [Wikipedia Orbo](https://www.wikipedia.org/wiki/Orbo); [Silicon Republic](https://siliconrepublic.com/?p=459296); [Museum of Hoaxes](https://hoaxes.org/weblog/comments/6498).
- **VASIMR 39-day Mars** [contested]: the claim needs 200 MW at about 1 kg/kW. What exists is 200 kW in pulses (72 ± 6% thruster efficiency) and 88 h at 80 kW (2021), with no flight — [The Space Review](https://thespacereview.com/article/1690/1); [Universe Today](https://www.universetoday.com/articles/zubrin-claims-vasimr-is-a-hoax); [SpaceNews](https://spacenews.com/?p=84168). A rebuttal argues Zubrin omits ongoing space-reactor work — [The Space Review](https://thespacereview.com/article/1896/1).
- **M2P2** [refuted at the proposed scale]: the claim was a 15 km bubble from a compact device. Kinetic analysis found negligible thrust, and Wikipedia notes the field falloff was optimistic and thrust "dramatically overestimated" — [NTRS 20030067378](https://ntrs.nasa.gov/citations/20030067378); [Wikipedia](https://en.wikipedia.org/wiki/Magnetic_sail).
- **Magsail (Andrews–Zubrin)** [mainstream, corrected]: drag was too high by 3.1× (Freeland 2015) — [Wikipedia](https://en.wikipedia.org/wiki/Magnetic_sail).
- **Yamato-1** [mainstream, underperformed]: design 8 kn with 8 kN per thruster. A logged trial gave 6.2 kn, and the system efficiency was about 15% — [Wikipedia](https://en.wikipedia.org/wiki/Yamato-1); [Seattle Times](https://archive.seattletimes.com/archive/19920720/1503104/engineering----may-the-force-be-with-you).
- **Open-cycle coal MHD "up to 60%"** [projection, never achieved]: channels ran near 17%, and MW plants extracted about 6–8% — [Wikipedia](https://en.wikipedia.org/wiki/MHD_generator); see §4.
- **Fuji-1 "30.2%"** [conflated]: Wikipedia credits Fuji-1, but OSTI attributes 30.2% to a new shock-tube-driven disk. Fuji-1 blow-down peaked at 18% — [OSTI 381436](https://www.osti.gov/biblio/381436); [OSTI 381435](https://www.osti.gov/biblio/381435).
- **Navy railgun "110 nmi"** [physics fine, programme failed]: 33 MJ was achieved, but barrel life was 12–24 rounds and the programme ended in 2021 — [NavWeaps](https://navweaps.com/Weapons/WNUS_Rail_Gun.php); [Task & Purpose](https://taskandpurpose.com/news/navy-electromagnetic-railgun-dead/).
- **Record categories** [mainstream, but headlines blur them]:
  - The 45.5 T (2019) and 48.7 T (2025) results are HTS test coils inside 31 T resistive magnets. The 2019 coil was damaged on quench — [Physics World](https://physicsworld.com/a/little-big-coil-creates-record-breaking-continuous-magnetic-field/).
  - Hefei's "world record steady field" (45.22 T) is for a working user hybrid — [Silicon Republic](https://siliconrepublic.com/machines/china-world-record-strong-stable-magnetic-field).
  - The 2025 Chinese all-superconducting figure is reported as both 35.1 T and 35.5 T, then 35.6 T at SECUF in 2026 (§1).
  - X3's "record" comes from press releases, though it is internally consistent (§8).
- **Inductrack L/D**: the patent gives about 300 at 500 km/h and Wikipedia about 200. The patent's "few m/s" transition speed conflicts with the 1.67 km/h implied by its own L/D examples (§3).

### Inferences
- A crate-level `claim_check` pattern helps here: store `{claim, measured, source}` triples and let tests assert the measured value from published inputs. VASIMR, X3, SPT-100, TFMC, Berry–Geim and Inductrack all have complete input sets, given in §§1–8.
- The common thread is magnetic pressure B²/2μ₀ and conductor σ. Devices fail where the claim needs stress beyond about 1 GPa (magnets above 50 T, SMES above 50 kJ/kg) or σB² beyond what seawater or plasma supplies (MHD ships, generators).

### Gaps
- The Searl effect generator, Howard Johnson motor and similar devices were not researched. No peer-reviewed measurement exists to set against them, and the Steorn verdict is the only documented test found. EmDrive, Mach-effect and other propellantless claims are excluded because the project already covers them.
