# Patents (2020–2026) with physics that standard calculations can check

Scope note: ten patents were read from their own text (USPTO scanned PDFs read page by page, or Google Patents' description) and checked. Their numbers are recomputed in `patents_checks.py` beside this file; run it with `python3 -I patents_checks.py`. Google Patents and Justia blocked automated access partway through, so several planned patents (CFS, Tokamak Energy, Helion, General Fusion, AOSense, Exail/iXblue, optical clocks) were found but not read. They are listed under Gaps. Inputs that a patent does not state are marked ASSUMED in the script. The verdict labels are neutral: *consistent*, *optimistic*, *internal error*, or *contradicts a bound*.

## Gravity and inertial sensing: do stated sensitivities sit inside the shot-noise and k_eff·T² limits?

### Takeaway
The one atom-gradiometer patent read with numbers, Birmingham/Delta G's US 11,269,111, states geometry rather than a sensitivity. Its 15 cm drop and 300 ms cycle put a 1 E target about 160 s of averaging away at the shot-noise limit. That is consistent with standard physics, and with the 20 E in 10 minutes the same group measured outdoors. The existing crate functions already cover every check needed.

### Cited Findings
- **US 11,269,111 B2, "Gravity gradiometer".** Inventors Lamb, De Villiers, Holynski and Bongs. Original assignee University of Birmingham, now listed as DELTA G Ltd. Priority 2017-12-15, granted 2022-03-08. Stated numbers: gradiometer baseline 0.5–2 m (preferred 0.8–1.2 m); Raman beam 1.8–2.2 cm; cooling beam 5.8–6.2 cm; measurement period 300 ms; collision-free time ≥3 s (~10⁻⁹ mbar); MOT gradient up to 15 G/cm; Rb natural linewidth 6 MHz; 6.8 GHz laser frequency; 15 cm drop; maximum tilt sin⁻¹(1/15) = 3.8°; target accuracy 1 E. The patent claims "10–100 fold" more sensitivity than current technology without giving a number. — [Google Patents US11269111B2](https://patents.google.com/patent/US11269111B2/en)
- The related field instrument (Stray et al. 2022) detected a 2 m tunnel with SNR 8 and a statistical uncertainty of 20 E within 10 min at a 0.7 Hz measurement rate. — [Nature/PMC 8866129](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8866129/); [Laser Focus World](https://www.laserfocusworld.com/test-measurement/article/14270785/quantum-gravity-gradient-sensor-detects-outdoor-tunnel)
- Recomputed (`patents_checks.py`):
  - The tilt is asin(1/15) = 3.82° (patent 3.8°, consistent).
  - The Rb D2 natural linewidth is 1/(2π·26.24 ns) = 6.07 MHz (patent 6 MHz, consistent); see [Steck, Rb-87 D line data](https://steck.us/alkalidata/rubidium87numbers.pdf).
  - A 15 cm drop lasts 175 ms, so with dropped atoms T ≈ 87 ms and k_eff T² = 1.23×10⁵ rad/(m s⁻²).
  - With ASSUMED N = 10⁶ atoms and contrast 0.5 on a 1 m baseline, shot noise gives 23 E per shot and 12.6 E/√Hz at the patent's 300 ms cycle, so 1 E takes about 158 s.
- Other granted 2020–2023 atom-interferometer patents surfaced but were not read: iXblue US 11,175,139 B2 (hybrid cold-atom/classical accelerometer, inventor Barrett, 2021-11-16); iXblue EP 3893003 B1 (2D-grating interferometer, 2022-12-14); AOSense US 10,535,981 B1 (2020-01-14) and US 11,852,485 B2 (Kasevich, 2023-12-26). — [Google Patents query results](https://patents.google.com/?q=(%22atom+interferometer%22+gravimeter)&status=GRANT&after=publication:20200101)

### Inferences
- **Verdict for US 11,269,111: consistent.** Its stated geometry supports its 1 E goal at the shot-noise limit after a few minutes of averaging. It states no sensitivity that could contradict a bound.
- **Calculator: already covered.** `gravimetry::atom::phase_sensitivity` (k_eff T²), `gravimetry::atom::gradient_noise`/`gradient_asd`, `gravimetry::local::averaging_time` and `gravimetry::local::cylinder` (the tunnel signal). A doc-test could reproduce the 3.82° tilt and the 158 s to 1 E from the patent's own geometry.

### Gaps
- No Exail/Muquans, AOSense, Infleqtion or AtomQuantum patent that states a µGal/√Hz figure was read. The AOSense and iXblue PDFs are scanned; US 11,175,139 was downloaded but its numeric pages were not reached before the stop.
- The University of Glasgow's Wee-g MEMS gravimeter is described as using a "patented spring", but no patent number was found. Published performance: about 10 µGal/√Hz noise floor and 3 µGal bias instability after 300 s. — [AQT](https://www.aqt.ac.uk/?p=924); [Glasgow eprint](https://eprints.gla.ac.uk/378022)
- No superconducting-gravimeter patent (GWR) was found.

## Magnetics: HTS magnets, optically pumped and NV magnetometers

### Takeaway
Both magnetometer patents read (Sandia/QuSpin and a fibre-coupled NV patent) agree with atomic physics in their key numbers. The Sandia/QuSpin zero-field offset of 7 kHz is 15–20% below what the patent's own buffer-gas table predicts (8.2–8.8 kHz). The NV patent's "sub-nT/√Hz" sits about 10–20 times above its photon shot-noise floor. The CFS/MIT magnet patent was identified, but its numbers were not extracted.

### Cited Findings
- **US 11,543,474 B1, "Optically pumped gradient magnetometer".** Assignees National Technology & Engineering Solutions of Sandia, LLC and QuSpin, Inc. Inventors Schwindt, Jau, Campbell and Shah. Filed 2021-07-20, granted 2023-01-03. Stated numbers:
  - two 1 cm³ ⁸⁷Rb cells, separated by the microwave wavelength, 4.4 cm;
  - ⁸⁷Rb gyromagnetic ratio "7 Hz/nT"; the stretched-state |1,1⟩–|2,2⟩ coherence shifts 3×7 = 21 Hz/nT, so the beat note is f_b = 7 kHz + 21 Hz/nT × ΔB;
  - N₂ buffer at 30 torr and 15 torr gives a measured zero-field offset of 7 kHz at about 100 °C;
  - Table 1 (from Vanier et al. 1982), for N₂: β = 547.9 Hz/torr, δ = 0.52 Hz/°C/torr, γ = −0.0013 Hz/°C²/torr;
  - coherence time about 1 ms, 500 Hz repetition, 795 nm pump square-wave modulated at 200 kHz, 12 µW probe at 780 nm, cells at 80–120 °C.
  — [USPTO PDF 11543474](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/11543474)
- Recomputed for US 11,543,474:
  - The Breit–Rabi low-field coefficient of |1,1⟩→|2,2⟩ with g_J = 2.00233 and g_I = −0.000995 is 21.015 Hz/nT, and g_F(F=2) gives 6.996 Hz/nT. Both are consistent with the patent.
  - c/ν_hf = 4.39 cm (consistent).
  - The patent's own Table 1 gives 15 torr × 547.9 Hz/torr = 8.22 kHz at 0 °C and 8.80 kHz at 100 °C. If the fill pressures were measured at 20 °C, the prediction is 8.20 kHz. The measured 7 kHz is about 15–20% lower.
- **US 11,852,701 B2, fibre-coupled NV-diamond sensor** (front page not read; the title in USPTO search is "Defect centre-based sensor"). Stated numbers:
  - lenses f₁ = 20 mm and f₂ = 8 mm; excitation NA ≈ 0.17, collection NA 0.23, fibre NA 0.27, and a 20% loss;
  - bias field 0.5–2 mT;
  - 0.3–1 mW of fluorescence on the photodiode gives 2.5 V; the patent estimates "on the order of 10⁸" NV⁻ centres from a 12 ns radiative lifetime;
  - spatial resolution 4 µm instead of 1 µm; sample 0.1–0.3 mm thick, volume 10 µm³–2 mm³;
  - claims "sub-nT/√Hz sensitivity", a dynamic range of about 1 T, Earth's field as 0.065 mT, MIT pulses of 1–10 mT above 100 Hz, and >1 W of 532 nm light and >1 W of microwaves in pulsed mode.
  — [USPTO PDF 11852701](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/11852701)
- Recomputed for US 11,852,701:
  - 0.3–1 mW at 700 nm (ASSUMED band centre) is 1.1–3.5×10¹⁵ photons/s. R × 12 ns gives 1.3–4.2×10⁷ NV, so the patent's 10⁸ is consistent as a lower-bound estimate that assumes unit collection.
  - With ASSUMED 1 MHz linewidth and 1% contrast, the CW-ODMR shot-noise limit, (4/3√3)·Δν/(γ C √R), is 46–84 pT/√Hz. The formula is from [Barry et al., Rev. Mod. Phys. 92, 015004 (2020)](https://doi.org/10.1103/RevModPhys.92.015004).
  - The energy resolution at 1 nT/√Hz in 2 mm³ is 7.5×10¹² ħ, far above the ħ-scale bound.
- **US 12,293,871 B2, "Grooved, stacked-plate superconducting magnets and electrically conductive terminal blocks…".** Assignees MIT and Commonwealth Fusion Systems LLC. Inventors Labombard, Granetz, Irby, Vieira, Beck, Brunner, Doody, Greenwald, Hartwig, Michael, Mumgaard, Radovinsky, Shiraiwa, Sorbom, Wright and Zhou. Filed 2023-09-28, granted 2025-05-06. It is a continuation of US 11,810,712. — [USPTO PDF 12293871](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/12293871)
- The matching published coil test is 20 T at 40.5 kA, 815 kN/m (Hartwig et al. 2024). The crate already reproduces this: `magnetics::magnets::lorentz_line_load(40.5e3, 20.1)` = 814 kN/m. — [MIT Technology Review](https://www.technologyreview.com/2024/04/23/1090425/mits-superconducting-magnets-are-ready-for-fusion/)

### Inferences
- **US 11,543,474: consistent on Zeeman physics.** The 21 Hz/nT and 4.4 cm are exact. There is mild internal tension on the 7 kHz offset against the patent's own Table 1 (about 17% low). Likely causes are cell-to-cell temperature or fill-pressure differences; it is not a physics violation.
- **Calculator (new).** `magnetics::magnets::alkali_zeeman_coefficient(g_j, g_i, I, F, m_F, F', m_F')` (low-field Breit–Rabi) and `buffer_gas_shift(p, T, β, δ, γ)`. Only the proton Larmor frequency exists today (`magnets::larmor_frequency`).
- **US 11,852,701: consistent.** The sub-nT/√Hz claim lies above the photon shot-noise limit for any plausible contrast and linewidth.
- **Calculator (new).** `magnetics::magnets::odmr_shot_noise(linewidth, contrast, photon_rate)`. The comparison with the ħ bound is already covered by `magnetics::magnets::energy_resolution`.
- **CFS 20 T coil: already covered** by `lorentz_line_load`, but only for the coil test, not for the patent's own numbers (see Gaps).

### Gaps
- The numeric pages of US 12,293,871 (field, current, turn-to-turn resistance, quench numbers) were downloaded but not read before the stop.
- Tokamak Energy's US 11,101,060 B2 ("Partially-insulated HTS coils") was downloaded (23 pages, scanned) but not read.
- No magnetic-navigation patent was checked.
- The assignee and grant date of US 11,852,701 were not captured. The grant date is probably 2023-12-26, since the number sits next to AOSense's US 11,852,485 B2, issued that day; this is an inference.

## Plasma and fusion: densities, temperatures and gains against Lawson, Bennett, beta and Greenwald

### Takeaway
Zap Energy's sheared-flow patent mostly states engineering and simulation numbers. One of them is a checkable unit conversion that is wrong: it gives 1 Torr·L as 2.5×10¹⁹ molecules at 273 K, against 3.54×10¹⁹ by the ideal-gas law. TAE's tunable-beam patent was read only through a search summary. No 2020–2026 Helion or General Fusion patent with numbers was read.

### Cited Findings
- **US 11,758,640 B2, "Apparatus and method for extended plasma confinement".** Assignee Zap Energy, Inc. Inventors Shumlak, Nelson and Meier. Filed 2022-06-15, granted 2023-09-12 (terminal disclaimer). Stated numbers:
  - simulated shear factors vsfa = 0.25, 0.5 and 0.75 keep the pinch to t = 14, 26.4 and 37 τ_A. Without shear, confinement degrades after 5 τ_A, with significant loss by 8 τ_A;
  - gas valves at z = −50, −75 and −100 cm (eight per ring, 45° apart); orifice ≥0.075 in; plenum ≥1 cm³;
  - pre-ionized plasma 1–10 eV at up to 50×10³ m/s, with injector supplies of 100–500 V and 1–500 A;
  - "10 Torr-L of neutral gas … in the same time as 500 Torr-L of pre-ionized gas (1/50)"; "1 Torr-L is proportional to 2.5×10¹⁹ molecules at 273 K";
  - RF pre-ionization at 400 kHz, 13.56 MHz or 2.45 GHz.
  — [USPTO PDF 11758640](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/11758640)
- Standard criterion: sheared axial flow suppresses the m = 1 kink when dv_z/dr ≥ 0.1 k V_A. — [Shumlak & Hartman, PRL 75, 3285 (1995)](https://doi.org/10.1103/PhysRevLett.75.3285)
- FuZE measurements (secondary summary of the PRL): about 200 kA pinch current, 1–2 keV, about 10¹⁷ cm⁻³, and neutrons for about 5 µs. — [OSTI 1961753](https://www.osti.gov/pages/biblio/1961753-sustained-neutron-production-from-sheared-flow-stabilized-pinch)
- Recomputed for Zap and FuZE:
  - 1 Torr·L at 273.15 K is 3.535×10¹⁹ molecules. The patent's 2.5×10¹⁹ is 0.71 of that, a 29% error.
  - The Bennett line density at 200 kA with Te + Ti = 2 keV (ASSUMED sum) is 6.24×10¹⁸ m⁻¹, which at 10²³ m⁻³ means a radius of 4.5 mm.
- **US 11,195,627 B2 (TAE Technologies, granted 2021-12-07).** Four of eight neutral-beam injectors tune from about 15 to about 40 keV during a shot at constant beam current. The patent claims up to 2× heating power and about 5× lower charge-exchange losses. These numbers come from a search-engine summary of the USPTO PDF, not a direct read. — [USPTO PDF 11195627](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/11195627); [Justia listing](https://patents.justia.com/patent/11195627)
- Recomputed for TAE: at constant current the injected power rises ×2.67 from 15 to 40 keV, so "up to 2×" absorbed heating is consistent with some shine-through.
- Context, not a 2020–2026 patent: Helion's FRC patent US 9,524,802 B2 dates from 2016. Its published compression target is a 10 keV FRC at 10²³ m⁻³ in a field of up to 20 T. Pressure balance at β = 1 checks: B²/2µ₀ = 1.59×10⁸ Pa equals n k T = 1.6×10⁸ Pa. — [arXiv 2501.03425](https://arxiv.org/html/2501.03425v1); [Google Patents US9524802B2](https://patents.google.com/patent/US9524802B2/en)

### Inferences
- **US 11,758,640: consistent in physics.** The simulated stabilisation grows monotonically with shear, as Shumlak & Hartman predict. There is an **internal error** in the Torr·L conversion (29% low).
- **Calculator: already covered.**
  - `plasmoid::pinch::bennett_line_density` / `bennett_temperature` (the Bennett check);
  - `plasmoid::discharge::gas_density` (the Torr·L check);
  - `magnetics::mhd::alfven_speed` (τ_A).
- **Calculator (new).** A one-line `plasmoid::pinch::shear_threshold(k, v_a)` = 0.1 k V_A would turn "vsfa" into a testable statement.
- **TAE: consistent** (power scaling only).

### Gaps
- No Helion, General Fusion or Zap patent stating an energy gain, Lawson triple product or Greenwald fraction was read. Google Patents and Justia blocked automated access (503, 403 and Cloudflare challenges) before the Helion and General Fusion searches could run.
- The TAE PDF (55 scanned pages) was downloaded, but its plasma parameters (Te, n, B, τ_E) were not read directly. One search summary gave Te ≈ 30–50 eV in one regime and "above about 1.0 keV" in the claims family; neither is verified.

## Light: laser–plasma accelerators, optical clocks and optical gravimetry

### Takeaway
LBNL/UC's frozen-gas plasma-channel patent gets its laser optics right: a Rayleigh range of 2.45 mm for a 25 µm spot at 800 nm. Its cryogenic recipe has a phase-diagram error, saying CO₂ freezes onto the channel wall at "−50 °C or lower". Solid CO₂ cannot form above its −56.6 °C triple point at these pressures, and at the stated 10 psi it deposits only below about −83 °C. No optical-clock patent with numbers was read.

### Cited Findings
- **US 11,438,999 B2, "Devices and methods for creating plasma channels for laser plasma acceleration".** Assignee The Regents of the University of California. Inventors Leemans, Swanson, Mao, Syversrud, Gonsalves and Sipla. Filed 2020-11-09, granted 2022-09-06. Stated numbers:
  - the Rayleigh range is "only a few millimeters for an 800 nanometer laser with a beam spot size of 25 microns";
  - in a hollow channel with wall density 10¹⁷ cm⁻³, TeV particles reach about 0.01 mm·mrad normalized emittance; in a near-hollow channel at about 10¹⁰ cm⁻³, emittance growth approaches 10⁻⁶ mm·mrad;
  - channel 1–50 cm long and 100 µm–5 mm across (about 1 mm preferred);
  - "when nitrous oxide is used … cooled to about −155 °C or lower; when carbon dioxide is used … about −50 °C or lower";
  - gas supply 0.1–150 psi (about 10 psi); frozen layer 10–500 µm; second gas H₂, He, N₂ or Ar at 5–30 torr;
  - laser about 80 µm across, at a wavelength printed as "about 800 microns" (a typo for nanometres), 1 nJ–10 J;
  - 10–50 kV discharge timed 100–600 ns before the laser; 500 A discharge pulse; ice shells 74 µm and 50 µm.
  — [USPTO PDF 11438999](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/11438999)
- Recomputed for US 11,438,999:
  - π w²/λ = 2.45 mm (consistent).
  - The cold wave-breaking field at 10¹⁷ cm⁻³ is 96√n₀ V/m = 30.4 GV/m; see [Esarey, Schroeder & Leemans, Rev. Mod. Phys. 81, 1229 (2009)](https://doi.org/10.1103/RevModPhys.81.1229).
  - Clausius–Clapeyron from the triple points (ΔH_sub about 26 kJ/mol for CO₂ and 23 kJ/mol for N₂O, both approximate): at 10 psi, CO₂ deposits below −83 °C and N₂O below −94 °C.
  - 5–30 torr of H₂ at 293 K gives n_e = 3.3×10¹⁷–2.0×10¹⁸ cm⁻³ if fully ionized.

### Inferences
- **Verdict for US 11,438,999.** The optics are consistent and the N₂O temperature is consistent (conservative). The **CO₂ "−50 °C" contradicts the CO₂ phase diagram**: −50 °C lies above the triple point, so no solid forms at any pressure the device uses. The "800 microns" is a typo.
- **Calculator (new).** A small `plasmoid::wakefield` module: plasma frequency, cold wave-breaking field E₀ = m_e c ω_p/e, Rayleigh range, dephasing and depletion lengths, and matched spot size, from Esarey et al. 2009. The nearest existing code is `gravitics::strong_field` (a₀, with explicit intensity conventions), which this module would reuse.

### Gaps
- No optical-clock patent with stated stability numbers was located. Vector Atomic's iodine-clock work (below 10⁻¹⁴ for days, Roslund et al., Nature 2024) appears in papers, not in a patent found here. — [Physics Today](https://physicstoday.aip.org/news/keeping-accurate-time-while-on-the-ocean)
- No optical-gravimetry patent was read. A 2019 Maryland LWFA patent (US 10,524,344, critical density 2×10²⁰–5×10²¹ cm⁻³) falls just outside the date range and was not read.

## Fringe gravity and propulsion: granted patents beyond Pais and US 12,417,870

### Takeaway
Four granted US patents from 2021–2025 were computed in their own terms against the standard bound.
- **Purvis** claims 7.72 N from 144 W, 1.6×10⁷ times the photon-rocket ceiling P/c.
- **Unspace**'s rotating-mass gravitational-wave "beam" is short of its own 6×10⁻¹¹ J by about 10²⁵ even for its most asymmetric rotor, and a 1–3 mm waist cannot form a beam at a 2×10⁵ m wavelength.
- **Aurigema & Buhler** measured 237 µN, 7×10⁻⁴ of the one-sided electrostatic force in their own rig, where conservation of the Maxwell stress predicts zero net force.
- **EdDrive**'s arithmetic is self-consistent, but a centre of mass moved by internal forces stays fixed in free space.

Each claim is decided by momentum conservation or by the quadrupole formula plus diffraction.

### Cited Findings
- **US 11,511,891 B2, "System and method for generating forces using asymmetrical electrostatic pressure".** Inventors and applicants Andrew Neil Aurigema and Charles Raymond Buhler IV, no corporate assignee. Filed 2019-11-19 (provisional 2018-11-19), granted 2022-11-29. Stated numbers and claims:
  - momentum P(t) = Mv = ε₀ t [E₂²A₂ − E₁²A₁], so F = ε₀[E₂²A₂ − E₁²A₁] (its eqs. 10–11). This rests on its own postulate d(Ut) = 0 (eq. 4);
  - "maximum allowed electric field in air ~10⁶ V/m … electrostatic pressure on the order of ~1 Pa";
  - test article: nine blades 0.25 in high, 0.25 in apart and 4 in long, a 0.25 in gap, epoxy-filled, a ground electrode of 4"×2" = 8 in², and a supply up to +40 kVDC. At +25 kV the average force was 237 µN at 0.1 mN resolution;
  - enlarging the ground from 8 to 15 in² was expected to raise the force about 1.8×; it measured 421 µN (×1.77);
  - force ∝ V²; optimal ground of about 23 triangles, with a 6× gain claimed;
  - for AC drive: "the frequency component … plays no part".
  — [USPTO PDF 11511891](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/11511891); [Google Patents](https://patents.google.com/patent/US11511891)
- Recomputed for US 11,511,891:
  - At 1 MV/m the Maxwell pressure ε₀E²/2 is 4.43 Pa (ε₀E² = 8.85 Pa), so the patent's "~1 Pa" is a factor 4–9 low.
  - In the test gap, 25 kV/6.35 mm = 3.94 MV/m gives 68.6 Pa, or 354 mN on 8 in². The measured 237 µN is 6.7×10⁻⁴ of that one-sided force.
  - The area ratio 15/8 = 1.875 against the measured 1.776.
  - Standard physics: for a static, isolated, charge-neutral system, the Maxwell stress integrated over a closed surface at infinity vanishes (the field falls as 1/r³), so the net force is zero (Jackson, *Classical Electrodynamics*, §6.7). A 7×10⁻⁴ imbalance is the size of systematic forces such as image and Coulomb forces on grounded surroundings, ion wind, and dielectric relaxation.
- Secondary source on the same team (Exodus Propulsion): a progression to "one full Earth gravity" on 30–40 g articles, against "about 10 mN" as the highest on a stacked system. 1 g on 30–40 g needs 294–392 mN, so the two statements disagree by 30–40×. — [The Debrief](https://thedebrief.org/nasa-veterans-propellantless-propulsion-drive-that-physics-says-shouldnt-work-just-produced-enough-thrust-to-defeat-earths-gravity/)
- **US 11,961,666 B2, "Pulsed e-field propulsion system".** Inventor James W. Purvis, no corporate assignee. Filed 2020-08-06, granted 2024-04-16. Stated numbers:
  - final design 7.72 N from 144 W (18.65 W/N), about 6 in × 4 in; both optimum designs are about 6–7 N at about 20 W/N;
  - the mechanism is Jefimenko's causal fields, with an "impulse" field of (µ₀/4π)·İ·geometry. İ ≈ 10⁷ C/s², about 1 C of isolated charge, a coil of 1 µH, 1 Ω and a 2 µs decay, and 4L/(R²C) ≈ 2;
  - cited prior art: Biefeld–Brown 50 mN from 35–45.5 W; Cannae 8–10 mN from 10.5 W; NASA Eagleworks 0.091 mN from 17 W; Woodward MET 6–7 µN from 100 W;
  - the thrust figures are predictions; no measured thrust of the final design is reported.
  — [Google Patents US11961666B2](https://patents.google.com/patent/US11961666B2/en)
- Recomputed for US 11,961,666:
  - P/c = 4.80×10⁻⁷ N, so the claim is 1.61×10⁷ times the photon-rocket ceiling.
  - Its own circuit checks out: the envelope 2L/R = 2 µs, C = 2 µF follows from 4L/(R²C) = 2, and (µ₀/4π)·10⁷ A/s = 1 V/m per unit geometry factor.
  - The cited prior art runs from 1.6×10³ (Eagleworks) to 3.7×10⁵ (Biefeld–Brown in air) times P/c.
- **US 11,149,719 B2, "EdDrive propellantless propulsion system".** Inventor Edward Von Bargen, no corporate assignee. Filed 2019-09-09, granted 2021-10-19. Stated numbers: 1 kg prototype on a 12 V battery; 36 in in the first second; Δv about 2.045 mph; 914.4 mN "per cycle" and 76.2 mN average over a 12 s cycle; projectile about 50% of host mass; gears engage 1–2 ms after rebound. — [Google Patents US11149719B2](https://patents.google.com/patent/US11149719B2/en)
- Recomputed for US 11,149,719:
  - 0.9144 m in 1 s and 0.9142 m/s give 0.914 N for 1 kg, and 914.4/12 = 76.2 mN, so the arithmetic is self-consistent.
  - Standard physics: internal forces leave the centre of mass fixed, and the hull oscillates by m_p/(M+m_p) × stroke (0.10 m for an ASSUMED 0.30 m stroke). Net drift on a surface comes from friction asymmetry.
- **US 12,494,696 B2, "Semi-coherent gravitational wave beam generation…".** Assignee Unspace LLC. Inventors Lentz, Crombe, Stephens and Laske. Priority 2022-09-15, granted 2025-12-09. Stated numbers:
  - tungsten/aluminium rotatable masses with R 50–200 mm and T 10–100 mm (up to metres in some embodiments);
  - 0.1–100,000 RPM (32,000 RPM in narrower embodiments) and tip speed 11–90.88 m/s;
  - frequency 0–1700 Hz (500–1700 Hz preferred) and wavelength 2×10⁵–6×10⁵ m;
  - beam waist 0.0010–0.0030 m and aggregate beam energy >6×10⁻¹¹ J;
  - communication through the Earth, across the US and the Atlantic, and to Mars.
  - No strain, radiated power or quadrupole formula is given.
  — [Google Patents US12494696B2](https://patents.google.com/patent/US12494696B2/en)
- Recomputed for US 12,494,696:
  - c/f = 1.76×10⁵–6.0×10⁵ m, so the wavelength statement is internally consistent.
  - For the most favourable shape, an ASSUMED 77 kg tungsten bar of 0.4 m × 0.1 m × 0.1 m at 90.88 m/s tip speed (ω = 454 rad/s, f_GW = 145 Hz), the quadrupole power (2/45)GM²L⁴ω⁶/c⁵ is 1.65×10⁻³⁶ W. The claimed 6×10⁻¹¹ J is 3.6×10²⁵ times the energy radiated in 1 s.
  - At 32,000 RPM the tip would reach 670 m/s, and ρv² = 8.7 GPa exceeds tungsten's strength of about 1–1.5 GPa.
  - Axisymmetric discs radiate exactly zero.
  - A 1.5 mm waist at λ = 1.76×10⁵ m has a Rayleigh range of 4×10⁻¹¹ m, so the radiation cannot be collimated into a beam.

### Inferences
- **US 11,511,891: contradicts a bound** (closed-system momentum conservation through the Maxwell stress tensor). The measured force is 7×10⁻⁴ of the local electrostatic force, which sets the systematic-control level a decisive test needs. The decisive experiment is a torsion balance in vacuum with the supply on board.
- **Calculator.** `gravitics::propulsion::electrostatic_pressure(E)` = ε₀E²/2, plus a doc-test that the closed-surface Maxwell stress of a dipole vanishes. The existing `gravitics::electrogravity` image force and `gravitics::propulsion::ion_wind_thrust` give the systematics.
- **US 11,961,666: contradicts a bound** (P/c), by 1.6×10⁷. **Calculator: already covered** by `gravitics::propulsion::photon_thrust`, plus the patent's numbers as a test row.
- **US 11,149,719: contradicts momentum conservation** for free-space use; its table-top numbers are self-consistent. **Calculator (new):** `gravitics::propulsion::internal_mass_excursion(m_host, m_proj, stroke)`, which is zero net and bounded oscillation.
- **US 12,494,696: contradicts two bounds**, the quadrupole formula (by about 10²⁵) and diffraction (a waist far below the wavelength). Its wavelength and frequency are self-consistent. **Calculator: already covered** by `gravitational::laboratory::spinning_rod_power`, `rotor_power` and `rotor_wavelength`. New: a one-line `gravitational::laboratory::rayleigh_range(w0, λ)`, or a cross-domain analogy in `physica::analogies` with the Gaussian beam.

### Gaps
- US 12,226,246 B2 (Weng-Dah Ken, "Ultra-precision timing clock method", 2025-02-18; "virtual space-time lattice") was found in the gravitational-wave search but not read.
- US 12,512,772 B2 (application US 2023/0170827 A1, "Magnetic levitation based low-gravity system", priority 2021-11-30) was downloaded but not read. It is a diamagnetic-levitation simulator, not gravity shielding, and `magnetics::magnets::levitation_product` already covers the physics. — [Google Patents US20230170827A1](https://patents.google.com/patent/US20230170827A1/en)
- No granted 2020–2026 patent claiming gravity shielding or inertial-mass reduction, other than the Pais patents, turned up in the searches run. Additional Navy/Pais patents were left out under the no-military rule.
- The number of US 2020 "Exodus" patent that The Debrief cites (WO 2020/159603 A2) was not reconciled with US 11,511,891.
