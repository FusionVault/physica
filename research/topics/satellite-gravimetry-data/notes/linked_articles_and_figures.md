# Three linked quantum satellite gravimetry articles: content, figures and how to reproduce them

Scope: (1) Open Access Government "Why quantum satellite gravimetry?"; (2) Zingerle et al. 2024, AGU Earth and Space Science, DOI 10.1029/2024EA003630; (3) phys.org Dec 2018 news item, plus the Hannover cold-atom gradiometer paper Douch et al. 2018 (Adv. Space Res.) suggested by the caller. Full texts were read where reachable. Wiley full text was read in a real browser, and all 12 Zingerle figures were viewed at full resolution. Some pages had bot protection (HAL Anubis, mediaTUM Anubis, ScienceDirect) and were NOT bypassed. GFZ servers refused connections.

**Important correction to the brief:** the phys.org 2018 item is NOT about Leibniz University Hannover / QUEST. It is a NASA Goddard + AOSense press release, and it cites no peer-reviewed paper. Douch et al. 2018 is a separate Hannover/LNE-SYRTE work that is still relevant, so it is analysed in its own section below.

---

## Q1. Article 1 — "Why quantum satellite gravimetry?" (Open Access Government): what it says, who wrote it, and what figures it has

### Takeaway
This is a 2-page opinion and policy piece by Bart Root (TU Delft). It has **no data figures, no equations and no reference list**. Nothing in it can be reproduced as a visualization. Its value is context: the mission lineage, why quantum sensors help (drift and bias stability), and pointers to CARIOQA, QSG4EMT, NASA's quantum gravity gradiometer pathfinder and MAGIC.

### Cited Findings
- Author: B.C. (Bart) Root, Assistant Professor, Delft University of Technology (Aerospace Engineering / Space Engineering, Planetary Exploration), ORCID 0000-0001-7742-1434. Published in the Open Access Government April 2026 issue, pp. 196–197. First published 27 Feb 2026, modified 16 Apr 2026. DOI 10.56367/OAG-050-12256. Licence CC BY-NC-ND 4.0. — [OAG article](https://www.openaccessgovernment.org/article/why-quantum-satellite-gravimetry/205621/)
- Argument: space gravimetry gives global coverage (a full Earth map in about one month). Accelerometers separate non-gravitational forces, and orbit determination or inter-satellite ranging captures gravity. Mission history: CHAMP 2000, GRACE 2002, GOCE 2009, GRACE-FO 2018, and MAGIC (ESA–NASA) as the successor. — [OAG article](https://www.openaccessgovernment.org/article/why-quantum-satellite-gravimetry/205621/)
- Quantum opportunities: laser-frequency/ranging improvements could push distance uncertainty below 0.1 nm (projected), versus nanometre-level for GRACE-FO. Cold-atom interferometry (BECs, typically rubidium) could replace or supplement electrostatic accelerometers, with better stability (less drift and bias). The main barrier is low TRL, so the author argues for continued investment. — [OAG article](https://www.openaccessgovernment.org/article/why-quantum-satellite-gravimetry/205621/)
- Projects named: NASA Quantum Gravity Gradiometer Pathfinder; ESA/Horizon Europe-backed CARIOQA (cold-atom accelerometer pathfinder); QSG4EMT (quantum gravimetry mission architecture studies for mass-transport monitoring); MAGIC. No grant numbers are given. — [OAG article](https://www.openaccessgovernment.org/article/why-quantum-satellite-gravimetry/205621/)
- Only image: a decorative iStock header ("Space probe is flying towards the sun…", © Elen11/iStock, elements furnished by NASA). It has no axes, units or data. — [OAG article](https://www.openaccessgovernment.org/article/why-quantum-satellite-gravimetry/205621/)
- QSG4EMT is an **ESA** study contract (ESA Contract No. 4000138395/22/NL/SD), 2022–2024, prime TU Munich (lead R. Pail), with University of Bonn (J. Kusche) and Trieste (C. Braitenberg) participating. Final report dated Oct 2024 (~385 pp.) plus an executive summary. It covers ll-SST constellations, quantum gradiometry and combined concepts. — [TUM QSG4EMT page](https://asg.ed.tum.de/iapg/qsg4emt); [QSG4EMT final report PDF](https://www.asg.ed.tum.de/fileadmin/w00cip/iapg/QSG4EMT/QSG4EMT_FR_V1.0.pdf); [Exec summary](https://asg.ed.tum.de/fileadmin/w00cip/iapg/QSG4EMT/QSG4EMT_ExecutiveSummary.pdf); [Bonn ESA-QSG4EMT](https://www.igg.uni-bonn.de/apmg/de/research/research-projects/esa-qsg4emt); [Trieste DMG](https://dmg.units.it/en/node/43399)
- NASA QGGPf (Quantum Gravity Gradiometer Pathfinder) is a NASA ESTO Instrument Incubator Program project. It uses a single sensor with laser-cooled atoms as proof masses. Secondary press reports cite rubidium atoms, ~0.25 m³, ~125 kg, launch near the end of the decade, and Infleqtion/AOSense/Vector Atomic as partners (unverified secondary numbers). — [NASA ESTO Chiow PDF](https://esto.nasa.gov/files/quantum/Chiow.pdf); [APS 2025 abstract](https://meetings-archive.aps.org/smt/2025/mar-f07/9/); [Quantum Insider (secondary)](https://thequantuminsider.com/2025/04/16/nasa-to-launch-first-space-based-quantum-gravity-sensor-to-map-earths-hidden-shifts/)

### Inferences
- The QSG4EMT final report (TUM, lead Pail, who is also a co-author on Zingerle et al. 2024) is the most likely place to find reproducible figures for the topics Root's article raises. Zingerle et al. 2024 is effectively a peer-reviewed sibling of that work: same TUM group and the same simulator lineage.
- The OAG article calls QSG4EMT/CARIOQA "ESA and Horizon Europe". The evidence suggests CARIOQA is the Horizon Europe part and QSG4EMT is ESA-funded.

### Gaps
- I did not open the QSG4EMT final report or the OAG PDF issue. Their figures are not documented here.
- No CARIOQA paper was opened. The only citation found is Lévèque et al. 2022, ICSO, DOI 10.1117/12.2690536, taken from Zingerle's reference list.

---

## Q2. Article 2 — Zingerle et al. 2024 (ESS, DOI 10.1029/2024EA003630): authors, funding, methods and simulation setup

### Takeaway
This is a closed-loop simulation study by TUM, LUH and GFZ. It covers CAI/hybrid versus electrostatic accelerometers for (a) GRACE-like single-pair and MAGIC-like double-pair SST and (b) GOCE-like SGG.
- **SST:** temporal aliasing dominates every realistic scenario. Instrument gains only appear when the time-variable signal is attenuated to under 10% (GRACE-FO noise) or under 0.1% (CAI12).
- **SGG:** a future CAI gradiometer with a CAI gyroscope beats GOCE by more than 10× and raises the resolvable degree from 223 to 263. It is still essentially blind to time-variable gravity once gyro noise is included.

### Cited Findings
**Bibliographic and funding**
- Title "The Benefits of Future Quantum Accelerometers for Satellite Gravimetry". Authors: P. Zingerle (TUM IAPG), M. Romeshkani (LUH IfE), J. Haas (GFZ Hydrology), T. Gruber (TUM), A. Güntner (GFZ and Univ. Potsdam), J. Müller (LUH IfE), R. Pail (TUM). ESS 11(9), e2024EA003630. Received 2024-03-14, accepted 2024-08-03, published 2024-09-01. CC BY 4.0. — [Crossref](https://api.crossref.org/works/10.1029/2024EA003630)
- Funders: BMWi/BMWK projects 50EE2220 (A/B/C); DFG TerraQ (Project-ID 434617780, SFB 1464); DLR Q-BAGS (50WM2181). Computations ran at LRZ and the LUH cluster. Open access via Projekt DEAL. — [Wiley full text](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024EA003630); [Crossref](https://api.crossref.org/works/10.1029/2024EA003630)
- Data availability: supplementary data are at TUM mediaTUM, DOI 10.14459/2024MP1737747. **Simulation software is not shared (IP rights)**; the methodology is published. — [Wiley full text](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024EA003630)
- Open PDF copies: Wiley pdfdirect (GOLD OA), TUM mediaTUM doc 1761409, and a GFZ pubman item. — [Semantic Scholar API](https://api.semanticscholar.org/graph/v1/paper/DOI:10.1029/2024EA003630?fields=title,externalIds,openAccessPdf); [mediaTUM PDF](https://mediatum.ub.tum.de/doc/1761409/1761409.pdf); [GFZ pubman](https://gfzpublic.gfz-potsdam.de/pubman/item/item_5027867)

**Instrument scenarios (Section 2)**
- Accelerometers: SuperSTAR (GRACE/-FO), MicroSTAR (NGGM-class), Alvarez@350 (drag-free, 350 km) and Alvarez@500 (no drag-free, 500 km), which are scaled-down LISA-Pathfinder-type sensors; ONERA "Dalin"; CAI11 (white 1e-11 m/s²/√Hz, "realistic"); CAI12 (white 1e-12 m/s²/√Hz, "future"). All three axes are assumed equal. Hybrid noise is "optimally combining" an EA with a CAI. — [Wiley full text](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024EA003630)
- LRI variants: LRI_GraceFO (after Hauk et al. 2023), LRI_NG30 (projected to 2030), LRI_NG33 (projected to 2033). — [Wiley full text](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024EA003630)
- SST classes (Table 1): GRACE-FO = SuperSTAR + LRI_GraceFO; NGGM = MicroSTAR + LRI_NG30; CAI11 = Alvarez@500 + LRI_NG33 (similar: CAI11 + LRI_NG33); CAI12 = (Dalin+) CAI12 + LRI_NG33 (similar: Alvarez@350 + LRI_NG33 for SP, Dalin + LRI_NG33 for DP). Above ~0.01 Hz the LRI limits; below that, the accelerometer dominates. For CAI12 the LRI dominates the whole spectrum. "NGGM" here is an instrument class, not the ESA NGGM mission. — [Wiley full text](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024EA003630)
- SGG classes (Table 2):

  | Class | Accelerometer | Gyroscope / attitude |
  |---|---|---|
  | GOCE | 3D MicroSTAR × 0.2 (about 5× better) | reconstructed from 3D gradiometer + R Gyro |
  | Future EA | 3D Alvarez@350 or Dalin | – |
  | Realistic CAI | EA + CAI11 | R Gyro, 1e-8 rad/s/√Hz (after Douch et al. 2018) |
  | Future CAI | EA + CAI12 | CAI Gyro, 1e-9 rad/s/√Hz (after Savoie et al. 2018) |

  Only diagonal tensor components are used (six single-axis CAI accelerometers). — [Wiley full text](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024EA003630)

**SST simulation setup (Section 3)**
- Simulator: an adapted TUM full-scale SST simulator (Daras 2016; Daras et al. 2015). Range-rate observations use the short-arc integral-equation approach (Mayer-Gürr 2006). The stochastic model is arc-wise covariance matrices built from the instrument ASDs. Estimation is unconstrained least squares with a static SH parametrization. — [Wiley full text](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024EA003630)
- Orbits:
  - Single pair: polar, 89° inclination, 463 km.
  - Double pair: the same polar pair plus an inclined pair at 70°, 432 km.
  - Both have a ~7-day sub-cycle and correspond to scenario "3d_H" of the ESA NGGM/MAGIC science support study (Heller-Kaikov et al. 2023).
  - Inter-satellite distance is not stated in the paper; the intro mentions ~220 km for GRACE/-FO.
  - Source: [Wiley full text](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024EA003630)
- Forward ("true world") models: static GOCO05s; non-tidal ESA ESM AOHIS (Dobslaw et al. 2015); ocean tides EOT11a. All are taken to d/o 120. — [Wiley full text](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024EA003630)
- "Full-noise" scenario: all time-variable signals in the forward model, with realistic de-aliasing in the backward model (AO de-aliasing product from the ESA ESM plus GOT4.7 tides removed). "Product-only" scenario: static field only, so the instrument is the only error source. — [Wiley full text](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024EA003630)
- Retrieval periods and degrees: 31 days by default (Section 3 results), estimated to d/o 120. Weekly (7-day) solutions are used for hydrology: d/o 120 for the double pair, d/o 100 for the single pair (ground-track limit). Long-term product-only weekly runs cover the whole year 2002. — [Wiley full text](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024EA003630)
- Temporal-aliasing test: all time-variable components are multiplied by an attenuation factor (100%, 10%, 5%, 1%, 0%). With GRACE-FO noise the signal must drop below 10% before sensor noise matters; with CAI12 it must drop below 0.1%. — [Wiley full text](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024EA003630)
- Key SST results:
  - Product-only SP: NGGM to CAI11 gives about 1 order of magnitude improvement, driven by low-frequency (<0.01 Hz) accelerometer noise; the gain is largest in sectorial coefficients, meaning less striping.
  - DP beats SP by at least 1 order of magnitude, and DP is less sensitive to low-frequency noise.
  - In full-noise runs, all scenarios except GRACE-FO SP are aliasing-dominated, and all DP full-noise runs are practically identical.
  - Source: [Wiley full text](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024EA003630)

**SGG simulation setup (Section 4)**
- Setup:
  - 2 months of data at 5 s sampling.
  - GOCE-like orbit: 96.7° inclination, 250 km.
  - Forward model: GOCO05s plus ESA ESM **HIS** (non-tidal, without AO), both to d/o 180.
  - Least squares on Vxx, Vyy, Vzz with a static SH parametrization; formal errors agree with empirical errors.
  - Degree errors exclude near-zonal (polar-gap) coefficients (Sneeuw 2000).
  - Source: [Wiley full text](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024EA003630)
- Results:
  - Realistic CAI is gyro-limited and worse than GOCE.
  - Future CAI is ~2 orders better than Realistic CAI and more than 1 order better than GOCE.
  - Alvarez@350 matches future CAI at low degree but is worse than GOCE at high degree.
  - Time-variable gravity is observable up to ~d/o 35 only if gyro noise is excluded.
  - Max resolvable d/o is 263 (~76 km) for future CAI versus 223 (~90 km) for GOCE (2-month data), using d_min ≈ 6371 km·π/l_max.
  - The GOCE comparison curve comes from 2 months of real GOCE data (Pail et al. 2011).
  - Source: [Wiley full text](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024EA003630)

**Hydrology (Section 5)**
- Inputs:
  - Weekly product-only solutions for 2002 are synthesized to 1° grids of 7-day mean TWS anomalies (EWH).
  - Each grid is truncated at the d/o where the mean weekly error exceeds the mean HIS signal (Table 3): GRACE 42 SP / 92 DP; NGGM 73 / 120; CAI11 100 / 120; CAI12 100 / 120.
  - Reference: ESA ESM HIS, 1995–2006, d/o 120. A 12-year linear trend is removed per grid cell.
  - Perfect AO/HIS separation is assumed, and aliasing and omission are neglected.
  - Source: [Wiley full text](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024EA003630)
- Basins: 474 of the 520 GRDC "Major River Basins of the World" (46 dropped because they are too small for 1°). Latitude-dependent area averages are computed per basin, then RMSD against HIS over 2002. — [Wiley full text](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024EA003630)
- Extremes:
  - Thresholds are percentiles of the 1995–2006 HIS reference: 99/98/95 for wet and 1/2/5 for dry.
  - Metrics: I (basins with a detected extreme), II (correct-week peaks CP), III (incorrect peaks IP), IV (mean TWS deviation of correct peaks), V (false-positive basins), VI (peaks in false-positive basins).
  - A "peak" is any time step crossing the threshold.
  - Source: [Wiley full text](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024EA003630)
- Selected numbers (Table 4, double pair, dry 5% / wet 95%):

  | Scenario | Score (dry / wet) | IV mean deviation (mm, dry / wet) | V false-positive basins (dry / wet) |
  |---|---|---|---|
  | GRACE | 0.14 / 0.48 | 8.4 / 11.6 | 55 / 31 |
  | NGGM | 0.49 / 0.72 | 3.0 / 3.6 | 30 / 11 |
  | CAI11 | 0.7 / 0.77 | 2.0 / 3.1 | 15 / 7 |
  | CAI12 | 0.7 / 0.8 | 1.7 / 2.8 | 6 / 6 |

  The ESM reference has 205 dry and 173 wet basins with extremes. — [Wiley full text](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024EA003630)
- Against the MAGIC MRD (threshold 60 mm, target 6 mm EWH weekly at 166–217 km, i.e. d/o 120/92): NGGM, CAI11 and CAI12 double pairs beat the 6 mm target for extremes; the double pair is the biggest single gain. — [Wiley full text](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024EA003630)

### Inferences
- The paper text cross-references "Figure 1b" and "Figure 1a" where it clearly means Figures 2b and 2a (Section 3.2, full-noise paragraphs). This is an editorial slip; treat Figure 2 as the degree-error figure.
- The "optimal combination" for hybrid accelerometers is not given as a formula. A reasonable reproduction is inverse-variance PSD fusion, S_hyb(f) = [S_EA(f)⁻¹ + S_CAI(f)⁻¹]⁻¹, with the CAI branch cut off above its Nyquist (the figure shows the CAI lines ending near 3·10⁻² Hz). This is my assumption, not the paper's.

### Gaps
- Inter-satellite distances for the 3d_H pairs, the arc length of the short-arc approach, and the exact SST sampling rate are not stated in the paper. They are in Heller-Kaikov et al. 2023 and Daras 2016, which I did not open.
- The contents of the supplementary dataset (mediaTUM 10.14459/2024MP1737747) could not be inspected: Anubis bot-check and fetch timeouts. Whether it holds degree-error curves, SH coefficients or TWS grids is unknown.

---

## Q3. Article 2 — every formula in Zingerle et al. 2024 (exact, from the article's MathJax TeX)

### Takeaway
The paper has 9 numbered equations: 5 accelerometer ASDs, 3 LRI ASDs and a hydrological detection score. Gyro and CAI levels are given inline. Degree-error and EWH conversion formulas are not written out; the paper defers to Heller-Kaikov et al. 2023.

### Cited Findings
All from the [Wiley full text](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024EA003630). These are one-sided ASDs; units are m/(s²·√Hz) for accelerometers and m/√Hz for the LRI.
- (1) Alvarez@350: `acc_A350(f) = 3e-8·f² + 4e-13·sqrt(1 + 7e-4/f + (3e-4/f)²)`
- (2) Alvarez@500: `acc_A500(f) = 3e-8·f² + 5e-13·sqrt(1 + (1/f)^(2/3))`
- (3) Dalin: `acc_D(f) = 2.4e-13·f^(-1/4) + 1e-18·f^(-2) + 8e-11·f²`
- (4) SuperSTAR: `acc_S(f) = 1e-10·sqrt(1 + 0.005 Hz / f)`
- (5) MicroSTAR: `acc_M(f) = (1/√2)·1e-11·sqrt( (1e-3/f)⁴ / ((1e-5/f)⁴ + 1) + 1 + (f/1e-1)⁴ )`
- CAI11: `1e-11` (white); CAI12: `1e-12` (white).
- (6) LRI_GraceFO: `lri_G(f) = 2.2e-10/sqrt(f) + 1e-12/f²`
- (7) LRI_NG30: `lri_NG30(f) = 2.2e-10/sqrt(f) + 1e-13/f²`
- (8) LRI_NG33: `lri_NG33(f) = 1.1e-10/sqrt(f) + 5e-14/f²`
- Gyros: R Gyro `1e-8 rad/s/√Hz`; CAI Gyro `1e-9 rad/s/√Hz`. GOCE accelerometer = MicroSTAR × 0.2.
- (9) `Score = (CP − IP) / N_ref` (CP = correct-week peaks, IP = wrong-week peaks, N_ref = reference peaks).
- Resolution rule: `d_min ≈ 6371 km · π / l_max`.
- Figure 1 caption: accelerometer ASDs are multiplied by √2 to express them as line-of-sight (inter-satellite) acceleration.

### Inferences
These are standard relations, needed for reproduction but not printed in the paper. I checked them numerically against Figure 1.
- **LRI to acceleration:** `lri_acc(f) = (2πf)² · lri(f)`. Checks: LRI_GraceFO at 0.1 Hz gives ≈3.1e-10, and at 1e-3 Hz gives ≈4e-11 m/s²/√Hz, matching Fig. 1a.
- **Product (combined) noise (Fig. 1b):** `prod(f) = sqrt( (√2·acc(f))² + ((2πf)²·lri(f))² )`. Checks:
  - GRACE-FO at 1e-3 Hz ≈ 3.5e-10, as in Fig. 1b.
  - CAI12 + LRI_NG33 at 1e-3 Hz ≈ 2.5e-12, flat at low frequency because the LRI dominates, as in the paper.
- **Top axis of Fig. 1 ("Corresponding SH degree"):** `l ≈ f / f_rev` with `f_rev = 1/T_orbit`. For 463 km, T ≈ 5620 s, so f_rev ≈ 1.78e-4 Hz, which matches the dotted vertical line at degree 1.
- **EWH SH coefficients** (Wahr et al. 1998 convention): `ΔC^EWH_lm = (a·ρ_E / (3·ρ_w)) · (2l+1)/(1+k'_l) · ΔC_lm`, with a ≈ 6378 km, ρ_E ≈ 5517 kg/m³, ρ_w = 1000 kg/m³, and k'_l the load Love numbers.
- **Degree error:** `σ_l = sqrt( Σ_m (ΔC_lm² + ΔS_lm²) )` in the chosen functional (geoid: multiply ΔC by R; EWH: as above). The paper does not say whether it uses degree amplitude or degree RMS (÷(2l+1)). Use Heller-Kaikov et al. 2023's convention.
- **Basin-size thresholds in Fig. 8 equal (πR/l_max)²:** d/o 120 gives 166.8 km → ~2.8·10⁴ km², and d/o 92 gives 217.6 km → ~4.7·10⁴ km², exactly the dashed-line labels.
- **Polar-gap exclusion for SGG degree errors** (GOCE-like i = 96.7°): the usual rule of thumb (Sneeuw & van Gelderen) drops orders m < |90° − i|·l (in radians ≈ 0.117·l). The paper only says "near-zonal coefficients excluded (Sneeuw 2000)", so the exact cut is an assumption.

### Gaps
- No explicit hybrid-combination formula, no observation-equation formulas (short-arc range-rate), and no gradient observation equations are printed in the paper.

---

## Q4. Article 2 — every figure: what is plotted, axes, units, colour scales, projections, and how to reproduce it

### Takeaway
All 12 figures were viewed at full resolution.
- **Figs 1, 2, 4, 5, 6, 7** are MATLAB-style log-y line plots (ASD vs frequency; degree error vs degree).
- **Fig 3** pairs SH coefficient triangles (log colour, jet-like colormap) with global EWH error maps (diverging red–white–blue, −2..2 mm, pseudo-cylindrical Robinson/Equal-Earth-like projection).
- **Figs 8–12** are Python/Excel-style hydrology statistics: log–log scatter, box plot, bars, and categorical basin maps.

Figures 1, 3 (as an analogue), 5, 7 and the hydrology statistics can be reproduced from public data plus a least-squares simulator. Exact reproduction of Figs 2, 4 and 6 needs a full SST/SGG closed-loop simulator, and the authors' code is not public.

### Cited Findings
Captions are paraphrased from the [Wiley full text](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024EA003630). Axis details come from my reading of the full-resolution figure images on the same page.

- **Fig 1 – instrument ASDs.**
  - What it shows: (a) accelerometer ASDs (solid: electrostatic; bold dotted: CAI11 orange, CAI12 red) and LRI ASDs expressed as range acceleration (thin dashed), with accelerometers ×√2. (b) Product noise for the four SST classes: GRACE-FO red, NGGM blue, CAI11 green, CAI12 purple.
  - X axis: frequency [Hz], log, ~1.4e-4 to 1e-1. Top axis: "Corresponding SH degree [1]" (1…400+), with a dotted vertical line at degree 1.
  - Y axis: acc/sqrt(f) [m/s²/sqrt(Hz)], log, ~5e-13 to 1e-9.
  - CAI dotted lines stop near 3e-2 Hz.
  - Data: Eqs. 1–8 only.
- **Fig 2 – degree errors in EWH, monthly (31-day), d/o 120.**
  - Black: monthly reference signal (ESM HIS + AO de-aliasing errors). Coloured: scenario minus truth. Solid = SP, dashed = DP.
  - (a) Product-only. X: Degree 2–120. Y: [mm] EWH, log 1e-4 to 10.
  - (b) Full-noise. X: Degree 2–60. Y: log 0.3 to 10. Orange "Error free (SP/DP)" curves are runs with negligible white range-rate noise.
- **Fig 3 – coefficient triangles and spatial errors**, product-only SP, NGGM (a, b) vs CAI11 (c, d), in EWH.
  - Triangles: x = Order −120..120 (negative = sine S_lm, positive = cosine C_lm); y = Degree 0..120, increasing downward. Log colour bar ~1e-5 to 1e2 [mm] EWH with a jet-like colormap.
  - Maps: global error field, diverging red (negative) / white / blue (positive), −2..2 [mm] EWH, pseudo-cylindrical projection with coastlines and dashed graticule.
  - NGGM shows strong N–S stripes; CAI11 is much fainter.
- **Fig 4 – aliasing attenuation**, double-pair, (a) GRACE-FO and (b) CAI12.
  - Curves: 100% (full noise), 10%, 5%, 1%, 0% (product-only), plus black monthly HIS reference.
  - X: Degree 2–120. Y: [mm] EWH, log 1e-4 to 1e2.
- **Fig 5 – SGG static field in geoid height.**
  - (a) Future CAI, Realistic CAI, GOCE, Future EA (Alvarez@350, Dalin, dashed) and black reference GOCO05s.
  - (b) Gyro impact: Dalin + CAI11 + CAI Gyro vs Dalin + CAI11 + R Gyro.
  - X: Degree 0–180. Y: "Geoid (m)", log ~1e-5 to 1e2.
- **Fig 6 – SGG time-variable gravity in EWH.**
  - Curves: Future CAI (without gyro noise) blue, Future CAI (with gyro noise) dark green, Only CAI12 + CAI Gyro light green, and the black thick "Ref monthly HIS".
  - X: Degree 0–180. Y: [mm] EWH, log 1–100.
  - The blue curve crosses the HIS signal near degree ~35.
- **Fig 7 – max recoverable d/o.**
  - Curves: Future CAI (green) vs GOCE (magenta) vs GOCO05s reference (black) in Geoid (m), log ~1e-5 to 1e2. X: Degree 0–280.
  - Dashed vertical lines at 223 (GOCE) and 263 (Future CAI), where each error curve crosses the signal.
- **Fig 8 – RMSD vs basin size** (scatter, log–log).
  - X: Basin size [km²], ~2e3 to 1e7. Y: RMSD [mm EWH]: about 1e-1 to 1e2 for single pairs (left), 1e-2 to 1e2 for double pairs (right).
  - Markers (matplotlib default colours): CAI12 blue dots, CAI11 orange, NGGM green ×, GRACE red +.
  - Legend means: SP 7.88 / 7.89 / 16.12 / 29.29 mm; DP 0.52 / 1.29 / 3.02 / 11.23 mm.
  - Right panel: dashed lines at ~2.8·10⁴ km² ("d/o 120", purple) and ~4.7·10⁴ km² ("d/o 92", red).
- **Fig 9 – box plot of DP RMSD** for the 237 basins larger than ~47,000 km² (d/o 92).
  - Y: RMSD [mm EWH], log 1e-2 to 1e2. Matplotlib boxplot: orange median line, green triangle mean, circles for outliers.
- **Fig 10 – bars, DP.**
  - Left: number of basins with extremes in 2002 (0–210) for HIS ref., GRACE, NGGM, CAI11, CAI12; dry 5% orange, wet 95% blue.
  - Right: mean deviation from the ESM reference [mm EWH], 0–12.
- **Fig 11 – basin maps of Score (Eq. 9)**, dry extremes 5%, DP, four panels (CAI12, CAI11, NGGM, GRACE).
  - Categorical classes: −inf (grey; extremes in the scenario but none in the reference), < −10.0 (magenta), 0 to −10.0 (red), 0–0.2 (yellow), 0.2–0.4 (light green), 0.4–0.6 (green), 0.6–0.8 (light blue), 0.8–1.0 (blue), no extremes (white).
  - Base map: GRDC basin polygons on a plain lon/lat (equirectangular-looking) world map.
- **Fig 12 – false positives, DP**, dual-axis bars.
  - Left axis: number of false-positive basins (0–80 dry; 0–60 wet). Right axis: mean number of false-positive peaks (0–8).
  - Scenarios: GRACE, NGGM, CAI11, CAI12. Dry is orange/red; wet is blue/dark blue.

### Inferences — reproduction recipes (my proposals, not from the paper)
- **Fig 1:** fully reproducible from Eqs. 1–8.
  - Evaluate on `f = logspace(log10(1.4e-4), -1, N)`. Multiply accelerometer ASDs by √2 and LRI by (2πf)².
  - Product = root-sum-square (Q3).
  - Top axis: l = f·T_orbit, with T_orbit ≈ 5620 s at 463 km.
  - Truncate CAI lines at ~3e-2 Hz to match the figure.
- **Figs 2 and 4 (SST closed loop):**
  1. Propagate two GRACE-like orbits (89°/463 km) and, for DP, an extra 70°/432 km pair, with the ~220 km separation assumed.
  2. Forward-model range-rates from GOCO05s (+ ESM AOHIS 6-hourly + EOT11a for full-noise) to d/o 120.
  3. Add coloured noise generated from the product ASD.
  4. Estimate a static SH set to d/o 120 by least squares with noise-whitened (arc-wise covariance) observations; for full-noise, remove ESM AO de-aliasing + GOT4.7 first.
  5. Compute per-degree EWH error vs truth.
  6. Reference curve: monthly HIS (+ AO de-aliasing error) signal degree amplitude in EWH.
  - Cheaper approximation: a semi-analytical / "lumped coefficient" or formal-error (normal-matrix) propagation (Sneeuw 2000). This gives product-only curves but not aliasing.
- **Fig 3:** from any simulated or formal error set ΔC/ΔS to d/o 120:
  - Triangle: plot `log10|ΔC^EWH_lm|` at (m, l), with S coefficients at −m.
  - Map: synthesize Σ ΔC^EWH·P̄_lm(sinφ)·cos/sin(mλ) on a 1° or 0.5° grid; diverging colormap −2..2 mm; Robinson or Equal Earth projection.
- **Figs 5 and 7 (SGG):**
  - Simulate 2 months at 5 s on a 96.7°/250 km orbit; compute Vxx, Vyy, Vzz in the LNOF from GOCO05s (+ HIS for Fig 6) to d/o 180 (to ~280 for Fig 7).
  - Add accelerometer-pair noise (gradient noise ASD ≈ √2·acc_ASD / baseline; GOCE baseline ≈ 0.5 m) plus gyro-driven centrifugal error (2ω·n_ω terms, see the Douch poster in Q6).
  - Solve by least squares, compute the geoid degree error `R·sqrt(Σ_m ΔC²+ΔS²)` excluding polar-gap orders, and plot against the GOCO05s signal degree amplitude.
  - Fig 7: find the crossing degree, then d_min = πR/l.
- **Figs 8–12 (hydrology):**
  1. Weekly TWS grids (1°, EWH) for 2002, truncated at the Table 3 cut-offs; detrend each cell by the 1995–2006 linear trend.
  2. Area-average over GRDC basin polygons with cos(φ) weights.
  3. RMSD over 2002 vs HIS, then scatter against basin area. Box plot for basins > 47,000 km².
  4. Percentile thresholds from 1995–2006 HIS; CP/IP counts per basin; Score per Eq. 9; choropleth with the classes above.
  - The paper's error grids are not public, so a reproducer must replace them with own simulated errors added to HIS (e.g. errors from a formal-error covariance).

### Gaps
- Exact colour maps (MATLAB `jet` vs `turbo`) and the exact map projection of Fig 3 are inferred from the image, not stated. The projection looks like Robinson or Equal Earth.
- The noise realization (random seed) and the exact orbit initial elements for 3d_H are unavailable.

---

## Q5. Article 3 — phys.org (20 Dec 2018) "Team creates and demonstrates first quantum sensor for satellite gravimetry": what it says and what lies behind it

### Takeaway
This is a NASA Goddard press release (by Lori Keesey) about a NASA Goddard + AOSense atom-interferometer gravity gradiometer prototype. It has no data figures (two photos) and no peer-reviewed paper. The closest technical sources are:
- a 2017 APS DAMOP presentation (Rakholia et al., NTRS 20170005272), giving the terrestrial prototype at 8·10⁻² E/√Hz, a 7·10⁻⁵ E/√Hz microgravity extrapolation, and a claimed ~1 cm EWH monthly at 200 km;
- a 2014 AGU talk (Luthcke et al., G23C-07).

It is **not** related to the Hannover/QUEST group.

### Cited Findings
- Date 20 Dec 2018. Source: NASA Goddard Space Flight Center (Lori Keesey). People named: Babak Saif (Goddard optical physicist), Lee Feinberg, John Mather. Partner: AOSense, Inc. (Sunnyvale, CA; quantum sensors and atomic clocks since 2004). Steven Chu is credited with inventing atom interferometry in the late 1980s. — [phys.org](https://phys.org/news/2018-12-team-quantum-sensor-satellite-gravimetry.html)
- Funding: NASA SBIR, Instrument Incubator Program, and Goddard IRAD. — [phys.org](https://phys.org/news/2018-12-team-quantum-sensor-satellite-gravimetry.html)
- Claims:
  - The prototype is much smaller than room-sized predecessors.
  - Its target is time-varying gravity, mainly water-mass change.
  - It could replace the two-satellite GRACE-FO system or add accuracy on a second satellite in a complementary orbit.
  - Mather suggests planetary interiors and new gravitational-wave frequency bands.
  - No sensitivity, orbit or resolution numbers are given.
  - Source: [phys.org](https://phys.org/news/2018-12-team-quantum-sensor-satellite-gravimetry.html)
- Images: (1) an atom-path-control demonstration spelling "NASA" (credit AOSense); (2) the terrestrial proof-of-concept gravity gradiometer (credit AOSense). Neither has data axes. — [phys.org](https://phys.org/news/2018-12-team-quantum-sensor-satellite-gravimetry.html)
- Original NASA release: [nasa.gov feature](https://www.nasa.gov/feature/goddard/2018/nasa-industry-team-creates-and-demonstrates-first-quantum-sensor-for-satellite-gravimetry)
- Technical backing: "Development of an Atom Interferometer Gravity Gradiometer for Earth Sciences", APS DAMOP, June 2017, report GSFC-E-DAA-TN43191.
  - Authors: A. Rakholia, A. Sugarbaker, A. Black, M. Kasevich (AOSense); B. Saif, S. Luthcke, L. Callahan, B. Seery, L. Feinberg, J. Mather, R. Keski-Kuha (GSFC).
  - Terrestrial prototype target 8·10⁻² E/√Hz; two atom sources with a 2 m baseline; T = 300 ms; 12ħk; N = 10⁶ atoms at 3 nK (Raman sideband cooling + magnetic lensing).
  - Microgravity extrapolation 7·10⁻⁵ E/√Hz.
  - Simulations claim monthly TVG of 1 cm EWH at 200 km resolution, 1–2 orders better than GRACE.
  - No orbit parameters or figures were visible.
  - Source: [NTRS 20170005272](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20170005272.pdf)
- Earlier: Luthcke, Saif, Black, Rowlands, "Earth Time Variable Gravity from a Spaceborne Cold Atom Gravity Gradiometer", AGU FM 2014, G23C-07.
  - Describes the Cold Atom Gravity Gradiometer (CAGG), selected by NASA ESTO IIP, measuring a single gradient tensor component.
  - Target noise floor 7·10⁻⁵ E/√Hz, compared with ~3·10⁻³ E/√Hz for GOCE short-term (≈40× higher).
  - A single satellite instead of GRACE's two; built-in rotation compensation.
  - Simulated TVG recovery for various orbits and multi-satellite configurations, with no numbers in the abstract.
  - Source: [AGU FM14 abstract](https://studio.m-anage.com/agu/fm14/webprogram/Paper29042.html)
- Successor programme: NASA QGGPf (JPL-led, ESTO IIP), with Infleqtion/AOSense/Vector Atomic named in press. — [NASA ESTO Chiow PDF](https://esto.nasa.gov/files/quantum/Chiow.pdf); [Quantum Computing Report](https://quantumcomputingreport.com/nasa-awards-infleqtion-20m-follow-on-contract-for-quantum-gravity-gradiometer-pathfinder-mission/)

### Inferences
- There is nothing to reproduce visually from the phys.org item. To visualize its claims, plot the NASA/AOSense noise levels against GOCE and against Zingerle's/Douch's gradiometer levels on a gradient-noise ASD chart. Use the gradient-noise ↔ accelerometer-noise relation n_grad = √2·n_acc / L.
- Unit conversion: 1 E = 1e-9 s⁻²; 1 mE = 1e-12 s⁻².
- The NASA numbers (7e-5 E/√Hz = 0.07 mE/√Hz) are about 70× better than the Hannover/SYRTE 5 mE/√Hz concept and below the Douch-poster requirement of 0.1 mE/√Hz for TVG. They are extrapolations, not demonstrated.

### Gaps
- No peer-reviewed journal paper from the Goddard/AOSense team was found. The NTRS PDF's slides (and any simulated EWH maps) could not be read; only metadata and abstract were.

---

## Q6. The Hannover paper — Douch et al. 2018 (Adv. Space Res. 61(5):1307–1323, DOI 10.1016/j.asr.2017.12.005) and its December-2018 sibling poster

### Takeaway
The full text was not directly readable: the HAL copy is behind an Anubis bot-check, which I did not bypass, ScienceDirect blocked access, and Semantic Scholar lists it as closed. I documented it from a detailed JRC technical report section that reproduces its main figure, and from the open follow-up Trimeche et al. 2019 (CQG, arXiv:1903.09828), which uses the same simulation pipeline.
- Key setup: CAI gradiometer (Carraz et al. 2014 design), white noise 5 mE/√Hz, 239 km orbit, 8-month (extrapolated) mission, EIGEN-6C4 truth, least squares.
- Key result: the quasi-inertial 3-axis mode beats the full GOCE mission; the nadir 1-axis (Vyy) mode does not.

A Hannover/Stuttgart poster (Douch, Knabe, Wu, Müller, Heinzel, created Dec 2018) is the closest real "Dec 2018 Hannover" item. It has reproducible TVG figures: gradient PSDs, geoid degree variances, EWH maps with a 3° Gaussian filter, and SNR triangles.

### Cited Findings
**Douch et al. 2018 (via secondary sources)**
- Authors: Karim Douch, Hu Wu, Christian Schubert, Jürgen Müller (Leibniz Univ. Hannover), Franck Pereira dos Santos (LNE-SYRTE, Paris). Adv. Space Res. 61(5), 1307–1323. HAL copy hal-01826234 (file article_CAIgrad_final.pdf, ~2.85 MB). — [Web search summary of HAL/RG records](https://hal.sorbonne-universite.fr/hal-01826234); [Semantic Scholar](https://api.semanticscholar.org/graph/v1/paper/DOI:10.1016/j.asr.2017.12.005?fields=title,externalIds,openAccessPdf)
- Concept: the Carraz et al. 2014 CAI gradiometer, sensitivity 3.5 mE/√Hz. This corresponds to a white one-sided ASD of 3.5·√2 ≈ 5 mE/√Hz, flat at low frequency (unlike GOCE). — [JRC report EUR 30371 EN (Travagnin 2020), §5.3.1–5.3.2](https://publications.jrc.ec.europa.eu/repository/bitstream/JRC121833/cold_atom_interferometry_satellite-based_gravimetry.pdf)
- Monthly TVG: detecting time-variable Vxx, Vyy, Vzz (ESA ESM signal) with SNR > 1 needs better than 0.1 mE/√Hz, so the paper focuses on the static field with GOCE as benchmark. A search-engine summary instead said "below about 5 mE/√Hz". **Conflict:** the JRC text (0.1 mE/√Hz) is the more detailed source and agrees with the Douch poster. — [JRC report](https://publications.jrc.ec.europa.eu/repository/bitstream/JRC121833/cold_atom_interferometry_satellite-based_gravimetry.pdf)
- Rotation constraint: the phase is `ΔΦ = k_eff·a·T² + 2·k_eff·(Ω × v)·T²`. Rotation perpendicular to the measuring axis must stay below ~1 µrad/s, against the GOCE-like orbital rate of ~1 mrad/s. So only two modes are possible:
  - nadir pointing: 1-axis instrument measuring only Vyy (cross-track);
  - quasi-inertial pointing: 3-axis instrument measuring Vxx, Vyy, Vzz.
  - Source: [JRC report](https://publications.jrc.ec.europa.eu/repository/bitstream/JRC121833/cold_atom_interferometry_satellite-based_gravimetry.pdf)
- Simulation: 8-month mission at 239 km; orbit error 2.6 cm (GOCE-like); star-tracker and gyro noise included; gradiometer noise 5 mE/√Hz. Reference GOCE solution: diagonal gradients from the whole mission (Nov 2009–Oct 2013, ~47 months). Compared against EIGEN-6C4 and GOCO05S. — [JRC report](https://publications.jrc.ec.europa.eu/repository/bitstream/JRC121833/cold_atom_interferometry_satellite-based_gravimetry.pdf)
- Main figure (reproduced as JRC Fig. 30):
  - Geoid-height degree variance, CAI vs GOCE. Solid = true errors (recovered − input), dashed = formal errors. "Comb" = combined-gradient solution.
  - The JRC notes the EIGEN-6C4 and GOCO05S black curves are identical and suspects an error.
  - Result: quasi-inertial beats nadir, especially above degree ~50. The 8-month quasi-inertial solution beats full-mission GOCE and also GOCE's 8 months at ≤239 km.
  - Source: [JRC report](https://publications.jrc.ec.europa.eu/repository/bitstream/JRC121833/cold_atom_interferometry_satellite-based_gravimetry.pdf)
- Douch et al. 2018 is the source of the "realistic gyroscope" level (~1e-8 rad/s/√Hz) used by Zingerle et al. 2024. — [Wiley full text](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024EA003630)

**Trimeche et al. 2019 (open follow-up; same pipeline, Douch, Wu and Müller co-authors)**
- Citation: CQG 36, 215004 (2019), DOI 10.1088/1361-6382/ab4548, arXiv:1903.09828. — [arXiv abs](https://arxiv.org/abs/1903.09828)
- Gravity recovery setup (Section VIII):
  - 71 days of real GOCE orbit (2 Mar–10 May 2013), mean altitude 239 km, positions spline-interpolated from 0.1 Hz to 1 Hz.
  - Truth: EIGEN-6C4 to d/o 360. Orbit noise σ = 2.6 cm.
  - Attitude: 3 star trackers + 1 gyro fused with a Wiener filter in the frequency domain (GOCE-style).
  - Gradiometer noise: 5 mE/√Hz white, plus centrifugal-correction error from angular-rate noise.
  - Gradients rotated GRF → Earth-fixed. Least squares to d/o 240 (58,077 unknowns).
  - Results extrapolated to 8 months by √time scaling.
  - Source: [arXiv HTML](https://arxiv.org/html/1903.09828)
- Figures: 15 (workflow), 16 (ASD of angular-velocity noise), 17a–d (formal-error triangles for Vxx, Vyy, Vzz and combined, log scale), 18 (degree medians of true errors in geoid height, 3-axis nadir vs 1-axis Vyy), 19 (8-month CAI vs GOCE full-mission and low-altitude subset). Result: about 2× improvement above degree 50, more at low degrees; Vzz is ~2× better than Vxx/Vyy, and the combined solution is dominated by Vzz. — [arXiv HTML](https://arxiv.org/html/1903.09828)
- Funding: ESA contract 4000112677/14/NL/MP ("Study of a Cold Atom Interferometer Gravity Gradiometer Sensor and Mission Concepts"); Lower Saxony "Vorab"/QUANOMET (QT3); DLR/BMWi 50 WP 1431 and 1700. — [arXiv HTML](https://arxiv.org/html/1903.09828)

**Douch, Knabe, Wu, Müller, Heinzel poster "What is Required to Recover the Time-variable Gravitational Field Using Satellite Gradiometry?"** (A0, PDF created 5 Dec 2018; affiliations GIS Stuttgart, IfE Hannover, AEI Hannover)
- Setup:
  - Only Vxx, Vyy, Vzz at 1 Hz over 29 days.
  - Static truth EIGEN-6C4 to d/o 360; time-variable truth = ESA updated ESM to d/o 180 (April 2006).
  - Gradiometer noise PSD rises as f⁴ above 0.01 Hz and as f⁻² below 1e-4 Hz (flat between).
  - Two polar circular orbits, 303 km and 361 km, with a 29-solar-day repeat.
  - Three total-noise levels: nominal, 5× and 10×.
  - Satellite dynamics from the HPS simulator (Pelivan et al. 2012); common-mode accelerations derived from GOCE data.
  - Source: [Douch_2018 poster PDF](https://www.gis.uni-stuttgart.de/forschung/doc/Douch_2018.pdf)
- Formulas on the poster:
  - `Vxx = Γxx + ωy² + ωz²` (likewise for yy, zz).
  - Error model: `V̂xx = Vxx + n_xx − 2ωy·n_ωy − 2ωz·n_ωz − 2dθz·Vxy + 2dθy·Vxz + ∇Vxx·dr + deal_err + (CMRR term × common-mode a_x / L)`.
  - `Γxx = (a_x^A1 − a_x^A4)/L`; measured `ã_x = s_x·a_x + n`.
  - `CMRR = (ŝ_j·s_i − s_j·ŝ_i)/(ŝ_j·s_i + s_j·ŝ_i)` for pairs 14, 25, 36.
  - `SNR(l,m) = log10( median_k |C^{t-var}_{l,m}(k)| / err(|C_{l,m}|) )`.
  - Source: [Douch_2018 poster PDF](https://www.gis.uni-stuttgart.de/forschung/doc/Douch_2018.pdf)
- Figures on the poster:
  - Fig 1: ASD of common-mode accelerations vs required CMRR.
  - Fig 2: star configuration of 3 accelerometer pairs.
  - Fig 3: √PSD of static and time-variable gradients along orbit vs requirement.
  - Fig 4: geoid degree variance of errors for 3 noise levels vs average hydrology signal.
  - Fig 5: recovered TVG fields (EWH, 3° Gaussian filter) and differences to the error-free 29-day mean.
  - Fig 6: SNR(l,m) triangle.
  - Fig 7: √PSDs of error contributors.
  - Tables: required position error (cm), CMRR, and angular-rate noise (rad/s/√Hz).
  - Source: [Douch_2018 poster PDF](https://www.gis.uni-stuttgart.de/forschung/doc/Douch_2018.pdf)
- Conclusions:
  - At 303 km, 0.1 mE/√Hz in the LNOF gives positive SNR to ~d/o 40.
  - 1 mE/√Hz is no better than GRACE.
  - Requirements: 0.1 arcsec/√Hz on the cross-track rotation angle; ~1 cm radial position knowledge for Vzz; de-aliasing errors must shrink below 2·10⁻³ Hz.
  - Source: [Douch_2018 poster PDF](https://www.gis.uni-stuttgart.de/forschung/doc/Douch_2018.pdf)

### Inferences
- To reproduce the Douch 2018 headline figure:
  1. Simulate 1 Hz diagonal gradients on a 239 km near-polar (GOCE-like 96.7°) orbit from EIGEN-6C4 or GOCO05s.
  2. Add 5 mE/√Hz white noise per component plus a centrifugal-term error from a ~1e-8 rad/s/√Hz gyro.
  3. Solve to d/o 240 by least squares, either full or with a block-diagonal/semi-analytical normal matrix.
  4. Plot the geoid degree median or RMS of true errors (solid) and formal errors (dashed).
  5. Scale from 71 days to 8 months by √(t1/t2).
  6. GOCE comparison curves: use the published GOCE time-wise (TIM) model formal errors, or EIGEN-6C4 − GOCE-only model differences.
- For a TVG visualization like the poster's Fig. 5, apply a 3° Gaussian smoother (Jekeli/Wahr weights W_l) to monthly ESM-based coefficient errors before EWH synthesis.

### Gaps
- Douch et al. 2018's own figure list, captions, equations, funding acknowledgements (QUEST / SFB 1128 geo-Q / CNES are plausible but **unverified**) and exact orbit inclination were not accessible.
- The original ScienceDirect/HAL text should be read by a human with access before claiming exact figure contents.

---

## Q7. Reference signals and public data: what is used and is it downloadable?

### Takeaway
All forward and reference models used by Zingerle (and Douch/Trimeche) are public:
- ESA ESM (updated, AOHIS + de-aliasing error products, 6-hourly, d/o 180, 1995–2006) from GFZ ISDC, DOI 10.5880/GFZ.1.3.2014.001;
- GOCO05s and EIGEN-6C4 static models (ICGEM);
- EOT11a and GOT4.7 tides;
- GRDC Major River Basins.

The ESM README was only seen via search-result snippets; GFZ servers refused direct connections during this session.

### Cited Findings
- Zingerle's forward models: GOCO05s (Mayer-Gürr et al. 2015), ESA ESM AOHIS (Dobslaw et al. 2015), EOT11a (Savcenko & Bosch 2012). De-aliasing: ESM AO product + GOT4.7 (Ray 2008). Hydrology reference: ESM HIS 1995–2006. Basins: GRDC (2020) "Major River Basins of the World". — [Wiley full text](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024EA003630)
- Updated ESA ESM:
  - ICGEM-format ASCII SH coefficients, d/o 180, 6-hour step, covering 1995–2006.
  - Components: A (IB-corrected atmosphere), O, H, I, S, plus their sum AOHIS.
  - Files reachable via DOI 10.5880/GFZ.1.3.2014.001.
  - Reference: Dobslaw et al. 2015, J. Geod. 89:505–513, DOI 10.1007/s00190-014-0787-8.
  - Source: [GFZ ISDC ESAESM2 README (via search snippet)](https://isdc-data.gfz.de/future_missions/esaesm/ESAESM2/README.txt); [GFZ ISDC ESM page](https://isdc.gfz.de/esmdata/esaesm/); [GFZ pubman Dobslaw 2015](https://gfzpublic.gfz.de/pubman/item/item_1010904)
- Douch poster: static EIGEN-6C4 d/o 360 and updated ESM d/o 180 (April 2006). Trimeche 2019: EIGEN-6C4 d/o 360 and the real GOCE orbit Mar–May 2013. — [Douch poster](https://www.gis.uni-stuttgart.de/forschung/doc/Douch_2018.pdf); [arXiv HTML](https://arxiv.org/html/1903.09828)

### Inferences
- GOCO05s, EIGEN-6C4, EOT11a and GOT4.7 are distributed in ICGEM `.gfc` format (ICGEM portal, GFZ). A single `.gfc` parser plus a fully-normalized Legendre synthesizer covers every forward and reference model needed.
- GRDC basin polygons come from the GRDC/BfG portal (registration/terms may apply).
- The "monthly HIS reference" curve in Figs 2/4/6: take 6-hourly HIS for a month, average to a monthly mean, subtract the long-term mean (or use deviations), convert to EWH, and compute the per-degree RMS over months. This is my reconstruction; the paper does not spell out its averaging.
- The "AO de-aliasing errors" component of the reference curve likely comes from the ESM "DEAL"/"AOerr" error products that accompany the updated ESM. I could not confirm the file names.

### Gaps
- I could not open the GFZ ISDC directory listing (ECONNREFUSED) to confirm exact AOHIS / DEAL / AOerr file names, or whether "ESAESM2" is a newer version than Dobslaw 2015.
- The GRDC basin dataset licence and download URL were not checked.

---

## Q8. Most relevant follow-up papers: what each article cites and what cites it

### Takeaway
Most of the network goes through Zingerle et al. 2024 and Douch et al. 2018. Article 1 cites nothing formally, and Article 3 cites no paper.
- **Most useful for reproduction:** Heller-Kaikov et al. 2023 (MAGIC 3d_H orbits, EWH degree-error convention), Daras & Pail 2017 (temporal aliasing), Dobslaw et al. 2015 (ESM), Abrykosov et al. 2019 (hybrid accelerometers), HosseiniArani et al. 2024/2025 (hybrid CAI + EA, open on arXiv), Romeshkani et al. 2024 (quantum gradiometry), Trimeche et al. 2019 (open, same pipeline as Douch).

### Cited Findings
- **Zingerle 2024 is cited by** (8 per Semantic Scholar, Oct 2026):
  - Chen et al. 2026, quantum gravity satellites in VLEO, static field (10.1029/2026ea005368)
  - Middendorf et al. 2026, extreme wet events in TWS from future missions (10.1016/j.asr.2026.07.075)
  - Doğu & Erkeç 2026 (10.1109/RAST69551.2026.11672430)
  - de Souza, Trombettoni, Braitenberg 2026, ultracold quantum gravimeters for geophysicists (arXiv 2601.00676)
  - Romeshkani et al. 2025, accelerometer data transplant (10.1029/2025EA004417)
  - Romeshkani et al. 2025, ML for transplanting accelerometer data (arXiv 2508.03522; 10.1007/s00190-026-02108-y)
  - HosseiniArani et al. 2025, combined classical and quantum accelerometers (10.1029/2024EA004187)
  - Romeshkani et al. 2024, quantum gravimetry for future satellite gradiometry (10.1016/j.asr.2024.11.072)
  - Source: [Semantic Scholar API](https://api.semanticscholar.org/graph/v1/paper/DOI:10.1029/2024EA003630?fields=title,externalIds,openAccessPdf,citationCount,citations.title,citations.year,citations.externalIds,citations.authors)
- HosseiniArani et al. "Combined Classical and Quantum Accelerometers…" preprint, arXiv 2405.11259. The hybrid instrument improves gravity solutions by 1–2 orders of magnitude at low and high degrees. — [arXiv HTML](https://arxiv.org/html/2405.11259v1)
- **Zingerle 2024 cites** 44 references (Crossref). Key DOIs: 10.1016/j.asr.2019.01.034; 10.1007/s00190-022-01659-0; 10.1103/PhysRevLett.120.061101; 10.1016/j.actaastro.2015.06.021; 10.1093/gji/ggad472; 10.1002/2017JB014250; 10.1007/s00190-014-0764-2; 10.1007/s00190-014-0787-8; 10.1016/j.asr.2017.12.005; 10.1093/gji/ggad266; 10.1029/2023EA003098; 10.1016/j.asr.2023.12.067; 10.1016/j.asr.2011.06.005; 10.1117/12.2690536 (CARIOQA, Lévèque et al. 2022); 10.1007/s00190-020-01401-8; 10.1007/s00190-011-0467-x; 10.1126/sciadv.aau7948; 10.1007/s00190-017-1042-x; 10.1007/s00190-011-0461-3; 10.12762/2016.AL12-11; 10.5880/G3P.2023.001; supplementary data 10.14459/2024MP1737747. Non-DOI items include Daras 2016 (DGK C 770), Mayer-Gürr 2006, Sneeuw 2000 (DGK C 527), Pail et al. 2015 (DGK 320), Knabe 2023 thesis. — [Crossref](https://api.crossref.org/works/10.1029/2024EA003630)
- **Douch 2018 is cited by** 42 (Semantic Scholar). Most relevant:
  - Trimeche et al. 2019, CQG (arXiv 1903.09828)
  - MOCASS (10.1007/s10712-019-09566-4; 10.1007/s00024-021-02756-5)
  - Lévèque et al. 2021 "Gravity field mapping using laser-coupled quantum accelerometers in space" (arXiv 2011.03382; 10.1007/s00190-020-01462-9)
  - Müller & Wu 2020 (10.1007/s00190-020-01401-8)
  - Hybrid electrostatic–atomic accelerometer 2022 (arXiv 2206.00634; 10.3390/rs14143273)
  - Quantum satellite mission with atomic clocks and CAI gradiometers 2022 (10.1088/2058-9565/aca8cc)
  - MOCAST+ 2023 (10.1007/s10712-022-09760-x)
  - New-generation accelerometer gradiometry 2023 (10.1016/j.asr.2023.08.023)
  - MaQuIs Mars quantum gravity mission 2023 (10.1016/j.pss.2023.105800)
  - "A path towards effective exploitation of quantum sensors in future satellite gravity missions" (arXiv 2404.07835; 10.1016/j.asr.2025.12.053)
  - Tiangong + CAI mission 2024 (10.1016/j.geog.2024.06.004)
  - Zingerle 2024; Romeshkani 2024
  - Review "Observation of the Earth gravity field from space… quantum physics" 2025 (10.1007/s40766-025-00076-1)
  - ESA's next-generation gravity mission concepts 2020 (10.1007/s12210-020-00875-0)
  - Source: [Semantic Scholar API](https://api.semanticscholar.org/graph/v1/paper/DOI:10.1016/j.asr.2017.12.005?fields=title,externalIds,openAccessPdf,citationCount,citations.title,citations.year,citations.externalIds)
- JRC overview with figure reproductions of Carraz 2014, Douch 2018 and Trimeche 2019: Travagnin 2020, EUR 30371 EN, CC BY 4.0. — [JRC report](https://publications.jrc.ec.europa.eu/repository/bitstream/JRC121833/cold_atom_interferometry_satellite-based_gravimetry.pdf)
- Article 1 context: ESA "ESA activities and perspectives on Quantum Space Gravimetry" (SPIE 12777, 127773K, DOI 10.1117/12.2690535); QSG4EMT final report. — [SPIE](https://proceedings.spiedigitallibrary.org/conference-proceedings-of-spie/12777/127773K/ESA-activities-and-perspectives-on-Quantum-Space-Gravimetry/10.1117/12.2690535.full); [QSG4EMT FR](https://www.asg.ed.tum.de/fileadmin/w00cip/iapg/QSG4EMT/QSG4EMT_FR_V1.0.pdf)

### Inferences
- The Crossref list for Zingerle is alphabetical, which lets the DOIs be matched to the in-text authors. **This mapping is my inference from the ordering, not checked against each DOI:**
  - Abrykosov 2019 → 10.1016/j.asr.2019.01.034
  - Alvarez 2022 → 10.1007/s00190-022-01659-0
  - Armano 2018 (LISA Pathfinder) → 10.1103/PhysRevLett.120.061101
  - Christophe 2015 → 10.1016/j.actaastro.2015.06.021
  - Daras 2024 (MAGIC) → 10.1093/gji/ggad472
  - Daras & Pail 2017 → 10.1002/2017JB014250
  - Daras 2015 → 10.1007/s00190-014-0764-2
  - Dobslaw 2015 → 10.1007/s00190-014-0787-8 (confirmed by GFZ)
  - Heller-Kaikov 2023 → 10.1093/gji/ggad266
  - Hauk 2023 → 10.1029/2023EA003098
  - Kornfeld 2019 → 10.2514/1.A34326
  - Kupriyanov 2024 → 10.1016/j.asr.2023.12.067
  - Lenoir 2011 → 10.1016/j.asr.2011.06.005
  - Müller & Wu 2020 → 10.1007/s00190-020-01401-8
  - Pail 2011 → 10.1007/s00190-011-0467-x
  - Savoie 2018 → 10.1126/sciadv.aau7948
  - Siemes 2018 → 10.1007/s00190-017-1042-x
  - Stummer 2011 → 10.1007/s00190-011-0461-3
  - Touboul 2016 → 10.12762/2016.AL12-11
- For a developer, the most reproduction-relevant open texts are Zingerle 2024 (CC BY), Trimeche 2019 (arXiv), HosseiniArani 2024 (arXiv 2405.11259), Lévèque 2021 (arXiv 2011.03382) and the JRC report (CC BY).

### Gaps
- Citation counts change over time (Semantic Scholar snapshot, Oct 2026).
- No citing-paper list exists for the OAG opinion piece or the phys.org news item.
