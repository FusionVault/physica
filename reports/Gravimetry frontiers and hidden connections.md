# Turn 259 gravity months into discoveries

The project's best next step is not finer maps. It is a small set of spherical-harmonic, uncertainty-aware detectors, each paired with a test fixed in advance on the 259 COST-G monthly fields. Three engines carry most of the new work:

- a matched-filter scan that needs no destriping (no GRACE precedent found in our searches);
- a minimum-norm "where water moves" transport field (no GRACE precedent found);
- a per-coefficient Kalman smoother that is exactly a spherical Matérn Gaussian process (a partial novelty: the method has relatives, the framing is new).

Between them they support about twenty hypotheses that could each be proven wrong. The hypotheses link mass change to El Niño, atmospheric rivers, heat, fire, earthquakes, Earth rotation, the core and dams.

Two of the tests have deadlines:

- **El Niño forecast.** A strong El Niño is under way (Oceanic Niño Index (ONI) +2.16 in July–September 2026). A dated forecast of the 2026–27 drop in global land water can be committed now, before any processing centre publishes the target months.
- **Kamchatka earthquake.** The Mw 8.8 event of July 2025 has 12 post-event months in our record. The literature has only a conference abstract on it.

User evidence points the same way. Operational users rank continuity, latency and calibrated uncertainty above resolution. No current tool combines several solutions, user-drawn regions, uncertainty and an up-to-date archive. The project's real advantage is statistical discipline, because published climate-link studies routinely assume four to five times more independent samples than autocorrelated monthly series contain.

Animations should pair playback that marks the data gaps with static small multiples and visible uncertainty. A staged plan below links each crate feature to a figure, a viewer feature and a test with a pass/fail gate.

**Every novelty verdict in this report comes from limited general web searches made on 24 October 2026, not from Google Scholar, Web of Science or Scopus. Treat "no prior found" as "not found by us", and re-verify in Google Scholar before claiming priority.**

## Users want basin numbers they can trust, and research now rewards error models

### Operational users rank continuity, latency and uncertainty above resolution

Gravity-derived water storage now feeds three kinds of operational product:

- **Weekly US drought indicators.** NASA Goddard publishes groundwater and soil-moisture drought indicators every week by assimilating GRACE-FO storage into a land model, at 0.125° over the continental US and 0.25° globally ([NASA/NDMC](https://nasagrace.unl.edu/)).
- **European groundwater.** The G3P groundwater product moved into operational use in the Copernicus Climate Change Service (C3S) in January 2025 ([GFZ G3P](https://www.gfz.de/en/section/earth-system-modelling/projects/g3p)).
- **A new C3S climate-variable service.** A C3S "Terrestrial Water Storage and Groundwater" service now delivers climate data records. Its authors say users **"need products with defined uncertainty and provenance"** ([GSTM2026 B.6](https://meetingorganizer.copernicus.org/GSTM2026/session/59417)).

ESA's requirements document for the future MAGIC mission (which pairs the US GRACE-C with ESA's Next Generation Gravity Mission, NGGM) sets out user priorities. In order: a continuous record of **at least 30 years**, then higher resolution, then daily-to-weekly sampling, then shorter latency. It puts GRACE-FO's average latency at **about 50 days** ([MAGIC MRD](https://www.asi.it/wp-content/uploads/2021/03/MAGIC_NGGM_MCDO_MRD_v1_0-signed2.pdf)).

Resolution has hard limits that users have learned the hard way. The US Geological Survey (USGS) found large mismatches in groundwater trends between wells and the GRACE assimilation product, and concluded that GRACE can support, but not replace, well monitoring ([USGS](https://www.usgs.gov/publications/groundwater-drought-united-states-spatial-and-temporal-variability)). A survey of scientists and practitioners found a gap between what researchers supply and what practitioners want, and called for transparency to build trust ([Xu et al. 2023](https://www.sciencedirect.com/science/article/abs/pii/S0022169423007874)). The public responds to simple framing: GFZ's World Water Day release used annual national balances in gigatonnes (Gt), deficit animations and a "water storage spiral" ([GFZ](https://www.gfz.de/en/press/news/details/jahresbilanz-2023-des-gesamtwasserspeichers-in-deutschland)).

Continuity is a live risk, and it affects our data directly:

- **Mission end date.** Sources disagree on when GRACE-FO ends. The CEOS database says September 2026 ([CEOS](https://database.eohandbook.com/database/missionsummary.aspx?missionID=686)); WMO OSCAR says 2029 or later ([OSCAR](https://space.oscar.wmo.int/satellites/view/grace_c_grace_fo)).
- **Operating mode.** GRACE-FO has flown in continuous "wide-pointing" mode since July 2023, with its laser interferometer in diagnostic mode ([Landerer 2025](https://cpaess.ucar.edu/sites/default/files/grace-fo-2025-presentation/landerer-felix.pdf)).
- **Successor launch.** A NASA Inspector General audit says GRACE-C will launch **before July 2029** within a $658M baseline, with 61% of reserves already spent ([NASA OIG](https://oig.nasa.gov/news/audit-of-nasas-upcoming-twin-satellite-mission-identifies-cost-and-schedule-risks/)). This supersedes the December 2028 date in our earlier report.
- **Quicklook data.** GRACE-C plans quicklook data with about **3-day latency** ([GSTM2026 A.3](https://meetingorganizer.copernicus.org/GSTM2026/session/59410)).

So software that assumes "monthly" is a liability. The ingestion layer should accept any time interval, and the processing chain should run as a stream.

### Existing tools fail on access, installation and side-by-side comparison

The recurring complaints in the trackers are mundane, which makes them actionable.

- **Fragile data access.** gravity-toolkit's scripted downloads broke repeatedly as NASA's PO.DAAC archive moved to the cloud, most recently in July 2025 ([#161](https://github.com/tsutterley/gravity-toolkit/issues/161)). Users also need both ICGEM `.gfc` and PO.DAAC `shm` formats ([#152](https://github.com/tsutterley/gravity-toolkit/issues/152)).
- **Unclear processing choices.** Users keep asking which low-degree coefficients to replace with satellite laser ranging (SLR) values, and where to find glacial isostatic adjustment (GIA) files ([#128](https://github.com/tsutterley/gravity-toolkit/issues/128); [#170](https://github.com/tsutterley/gravity-toolkit/issues/170)).
- **Installation and threading in SHTOOLS.** Open issues report that pyshtools forces a downgrade to numpy 1.x on Windows ([#505](https://github.com/SHTOOLS/SHTOOLS/issues/505)) and segfaults when called from several threads ([#385](https://github.com/SHTOOLS/SHTOOLS/issues/385)). A request for the Hessian of a spherical-harmonic (SH) expansion has gone unanswered since 2022 ([#356](https://github.com/SHTOOLS/SHTOOLS/issues/356)). Our `geopotential` crate already computes second partials, so that feature only needs documenting.
- **An overloaded calculator.** ICGEM's calculation service was refusing new jobs "due to high load" on the day of research ([ICGEM](https://icgem.gfz.de/calcgrid)).

Web portals cover only part of what users need:

- **GravIS** limits its reliable results to areas of **at least 100,000 km²**, leaves GIA errors out of its uncertainties, keeps most of its processing code proprietary, and shows only GFZ and COST-G solutions ([Dahle et al. 2025](https://essd.copernicus.org/articles/17/611/2025/essd-17-611-2025.html)).
- **The CU Boulder mascon tool** still serves outdated mascon versions with no uncertainty ([CCAR](https://ccar.colorado.edu/grace/about.html)).
- **Google Earth Engine's** JPL mascons stopped at **2024-09-30** ([Earth Engine](https://developers.google.com/earth-engine/datasets/catalog/NASA_GRACE_MASS_GRIDS_V04_MASCON?authuser=0)).

A 2026 review of GRACE data assimilation adds the decisive methodological point. It found **larger differences between different mascon solutions than between SH solutions**, and asked providers for "standardized Level-3 products and reliable error estimates" ([Springer et al. 2026](https://hess.copernicus.org/articles/30/985/2026/)).

The Rust ecosystem has almost no gravity tooling. The awesome-georust list has no gravity, geoid or SH entries ([awesome-georust](https://github.com/pka/awesome-georust)). Our crates have tiny user bases (`geopotential` had 40 downloads, `gravity-field` 8), so breaking changes cost nothing yet ([crates.io](https://crates.io/crates/geopotential)). Pure-Rust NetCDF and HDF5 readers now exist ([netcdf-reader](https://crates.io/crates/netcdf-reader)), which makes reading standard files in the browser plausible.

| Priority | Feature | Who needs it | Evidence | What the crates lack today |
|---|---|---|---|---|
| 1 | Polygon and basin time series with propagated uncertainty, several solutions side by side | Water managers, C3S users, journalists | Mascon solutions differ more than SH ones; GravIS lacks drawn polygons and other centres | Basin kernels (needs SH analysis), uncertainty propagation, reading CSR/GFZ/JPL/ITSG solutions |
| 2 | Standard Level-3 outputs with provenance: CF NetCDF, Cloud-Optimised GeoTIFF (COG), CSV/JSON | Assimilation modellers, C3S, GIS users | Calls for standardised Level-3 products and for "defined uncertainty and provenance" | Writers; a provenance record (filter, GIA, TN-13/14 choices) |
| 3 | Python wheels with no Fortran (PyO3, abi3 stable ABI) | The Python GRACE community | SHTOOLS numpy-2 and thread-safety failures | Bindings |
| 4 | Global drought layer: percentiles against 2002–present, onset and breakpoints | Drought monitors outside the US, NGOs | NASA's drought product is model-assimilated and US-centred | The state-space smoother and change-point detection (Section 2) |
| 5 | Ingestion at any cadence: monthly, 10-day, daily, GRACE-C quicklook, ESA-SING simulated fields | Flood and drought early warning; MAGIC preparation | MAGIC latency targets; GRACE-C 3-day quicklook | Remove the "monthly" assumptions from the time axis |
| 6 | Robust access plus an open mirror of the archive | Everyone | PO.DAAC breakage; ICGEM load refusals | The 259 fields are already held; CC BY 4.0 allows mirroring with attribution |
| 7 | Recompute in the browser (filter radius, scales) via WASM | Educators, analysts without Python | Today's portals serve fixed, server-made products | `geopotential` compiled to WASM |
| 8 | Groundwater separation (storage minus soil moisture, snow and surface water) | Groundwater managers | The top practitioner interest; G3P | Ancillary land data; defer |

### 2023–2026 research rewards error modelling over super-resolution

Machine learning on GRACE has mostly turned to water rather than gravity, and its strongest claims are contested.

- **Downscaling adds little real detail.** A 2024 critique shows that resampled 0.5° grids carry no physical information below GRACE's native ~3° scale ([Sci. Rep. 2024](https://www.nature.com/articles/s41598-024-60366-2)). At GSTM 2026, a downscaled product "added little" for estimating aquifer specific yield across India ([GSTM2026 B.4](https://meetingorganizer.copernicus.org/GSTM2026/session/59415)).
- **Simple baselines are hard to beat.** A NASA Goddard study found that **linear regression outperforms LSTM and Temporal Fusion Transformer** networks for predicting storage ([Nie et al. 2025](https://arxiv.org/abs/2510.10799)).

The defensible ML results are of two kinds. The first is reconstructions checked against independent gravity. TWSTORE combines climate-predictor models with low-degree SLR and DORIS gravity back to 1984, and warns that "GRACE-derived results should not be extrapolated to the past" ([Hacker et al. 2026](https://essd.copernicus.org/articles/18/1747/2026/)). The second is short-lead nowcasts that bridge the latency: a Bayesian convolutional network reaches **r = 0.92–0.95** over the 1–3 months before the data arrive ([Mo et al. 2025](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2024GL112677?af=R)). The lesson for us is to build evaluation tooling with mandatory linear baselines, not more networks.

The high-frequency frontier has moved below the monthly Level-2 coefficients:

- **5-day mascons.** CSR's 5-day mascons show elevated antecedent storage before **80–85% of 3,272 catalogued floods**. The fields are available only on request ([Rateb et al. 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11180122/)).
- **Range-rate residuals as event detectors.** Post-fit inter-satellite range-rate residuals, used as line-of-sight gravity differences, pick out cyclones and monsoon floods ([Cuadrat-Grzybowski et al. 2025](https://arxiv.org/abs/2503.04227)).
- **Multi-satellite combinations.** GOCO2025s combines **kinematic orbits of 18 low-Earth-orbit satellites** with GRACE-FO laser ranging, and publishes its full normal equations ([GOCO2025s](https://essd.copernicus.org/preprints/essd-2026-639/essd-2026-639.pdf)).

The clearest advance in estimation is treating background-model errors as random processes rather than subtracting the models:

- AOe07 supplies a 6-hourly error series for the atmosphere and ocean de-aliasing model, AOD1B ([Shihora et al. 2024](https://link.springer.com/article/10.1007/s00190-024-01832-7)).
- Non-stationary error covariances matter in closed-loop simulation ([Abrykosov & Pail 2025](https://link.springer.com/article/10.1007/s00190-025-01980-4)).
- GFZ now co-estimates the de-aliasing coefficients in its operational Level-2 processing ([EGU25-16562](https://meetingorganizer.copernicus.org/EGU25/EGU25-16562.html)).

Regional methods and GIA also moved:

- **Slepian functions.** An altimetry-regularised Slepian method improved Greenland accuracy by **25–58%** over the classical version ([Gong et al. 2025](https://link.springer.com/article/10.1007/s00190-025-02022-9)).
- **3-D GIA.** Using 3-D mantle viscosity in the GIA correction shifts the Antarctic loss estimate by about **12 Gt/yr** ([JPL poster](https://postdocs.jpl.nasa.gov/media/2025%20RPD/2025%20PRD-E-025.pdf); [Hightower et al. 2026](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2026JB033791)). A 2024 study found limited sensitivity, so the size of the effect is contested ([J. Geodyn. 2024](https://www.sciencedirect.com/science/article/abs/pii/S0264370724000309)).

Future missions and quantum pathfinders are all funded:

- **NGGM** had its final design (Phase B2/C) secured at ESA's November 2025 Ministerial, for a launch around 2032 ([ESA FutureEO](https://climate.esa.int/documents/3235/1_FutureEO_post-CM25_-_Floberghagen.pdf)).
- **CARIOQA**, Europe's quantum accelerometer pathfinder, entered Phase B in October 2025, and its Phase A orbit and accelerometer data are public ([EGU26-9366](https://meetingorganizer.copernicus.org/EGU26/EGU26-9366.html); [Zenodo](https://zenodo.org/records/17580678)).
- **NASA's quantum gradiometer pathfinder** received a further $20M in August 2026 ([Infleqtion](https://infleqtion.com/nasa-awards-infleqtion-20-million-as-worlds-first-quantum-gravity-mission-advances-toward-flight/)).
- **FEATS-G**, a €500k GFZ-led study, will run end-to-end quantum-gravimetry simulations from 2026 to 2029 ([GFZ](https://www.gfz.de/en/section/earth-system-modelling/projects/feats-g-quantifying-the-potential-of-quantum-space-gravimetry)).

The only full open recovery suite remains GROOPS, which is monolithic C++ under the GPL ([GROOPS](https://github.com/groops-devs/groops)). Mission-impact studies still use crude error proxies. A 2026 hydrology study added plain Gaussian noise of 1.9–42 mm to a proxy field ([Liaqat et al. 2026](https://hess.copernicus.org/articles/30/4969/2026/)).

This opens three roles for `satgrav`:

- **A temporal-aliasing layer.** AOe07-style perturbations injected into the closed loop, which today has no aliasing at all.
- **A light mission-error emulator.** Degree-dependent, anisotropic, correlated noise for GRACE-FO, GRACE-C, NGGM, MAGIC and hybrid quantum concepts, published in ESA-SING formats.
- **Validation against CARIOQA's public Phase A data.**

Two smaller extensions fit geo5d and the planetary case.

- **Clock-based height measurement (chronometric levelling).** A 457 km fibre link between two strontium clocks measured a geopotential difference of **3918.1(2.6) m²/s²**, against a geodetic value of 3915.88(0.30) ([PR Applied 2024](https://journals.aps.org/prapplied/abstract/10.1103/PhysRevApplied.21.L061001)). ESA's ACES clock ensemble launched to the ISS on 21 April 2025 ([ESA](https://www.esa.int/Science_Exploration/Human_and_Robotic_Exploration/ACES_Atomic_Clock_Ensemble_in_Space)).
  - A 1 mm geoid change is a fractional frequency of about 1.1×10⁻¹⁹ (g·1 mm/c²; our arithmetic). Hydrological corrections from COST-G therefore matter only once clocks pass the low-10⁻¹⁸ level that transportable strontium clocks now approach ([TerraQ](https://arxiv.org/pdf/2503.09272)).
  - geo5d should still offer clock-site differences computed from a gravity model, with commission and omission errors, as a 2026 IAG paper studies ([IAG 2026](https://link.springer.com/chapter/10.1007/1345_2026_308)).
- **Planetary fields.** GRAIL's measured degree-3 lunar Love number, **k₃ = 0.0163 ± 0.0007**, is 72% above the spherically symmetric value ([Park et al. 2025](https://www.nature.com/articles/s41586-025-08949-5)). That is a reason to make synthesis body-agnostic, so `satgrav` doubles as a GRAIL- or Mars-pair simulator.

## Seven spherical-harmonic algorithms have no GRACE precedent we could find

Our data shape what is optimal. Each COST-G file is a weighted combination of several analysis-centre solutions (nine in July 2026). It carries `errors empirical`: one calibrated standard deviation (sigma) per coefficient, and **no covariance matrix** ([COST-G July 2026 header](https://icgem.gfz.de/getseries/02_COST-G_/Grace-Grace-FO_RL02.1/GSM-2_2026182-2026212_GRFO_COSTG_BF01_0201.gfc)).

Fully normalised SH functions are orthonormal on the sphere. So decompositions such as empirical orthogonal functions (EOF), dynamic mode decomposition (DMD) and linear inverse models (LIM) can run directly on the 259 × 8,281 coefficient matrix, with no gridding or area weights. Dividing each coefficient by its sigma first turns plain PCA into noise-weighted PCA.

Two infrastructure pieces are missing and block several methods: SH analysis (grid to coefficients) and SH rotation. Both are cheap at degree 90. Rotation can be ported from MIT-licensed FastTransforms but not from GPL ducc0 ([FastTransforms](https://github.com/MikaelSlevinsky/FastTransforms); [ducc0](https://pypi.org/project/ducc0/)).

### A matched-filter scan detects abrupt mass changes without destriping

GRACE earthquake studies fit steps at times and places known in advance. The literature also documents that standard destriping "tends to distort the true coseismic signals in both amplitude and spatial pattern" ([EPS 2015](https://earth-planets-space.springeropen.com/articles/10.1186/s40623-015-0188-0)).

A generalised likelihood-ratio test (GLRT) in the coefficient domain avoids both problems (method M1). For a step template $g$ (a disc or cap load converted to Stokes coefficients with our Love numbers), the statistic

$$T=\frac{\big(\sum_{nm} g_{nm}\Delta d_{nm}/\sigma^2_{nm}\big)^2}{\sum_{nm} g_{nm}^2/\sigma^2_{nm}}$$

is χ²₁ under Gaussian diagonal noise. Here $\Delta d$ is the after-minus-before difference of the seasonally adjusted coefficients.

The key trick is that, for an axisymmetric template $K(\psi)=\sum_n k_nP_n(\cos\psi)$, the addition theorem gives $g_{nm}=k_n\bar Y_{nm}(\theta_0,\lambda_0)/(2n+1)$. So the numerator at every candidate centre at once is a **single SH synthesis** of $k_n\Delta d_{nm}/((2n+1)\sigma^2_{nm})$. A full scan over epochs and radii costs seconds.

Real stripes are correlated, so the diagonal-noise p-value is optimistic. The thresholds should therefore come from scan-statistic Monte Carlo: time-shuffled or sign-flipped residual months, taking the family-wise maximum over space, time and radius. Kulldorff-style scan statistics exist in other fields ([JOSS](https://joss.theoj.org/papers/10.21105/joss.00515)).

The only matched filters we found in GRACE work are for accelerometer calibration. Nearest neighbours are Slepian localisation without destriping ([Wang et al. 2012](https://geoweb.princeton.edu/people/simons/PDF/reprints/GRL-2012a.pdf)) and a GFZ order-dependent filter for Mw 8.0–8.2 events ([GFZ](https://gfzpublic.gfz-potsdam.de/pubman/item/item_5020737)). M1 is the engine for five of the hypotheses below (E1, E2, E3, E6 and R3).

### A minimum-norm transport field shows where water moves

Gravity measures mass change, not velocity. But the continuity equation $\partial\sigma/\partial t+\nabla\!\cdot\!F=0$ has a unique least-motion solution, the curl-free flux $F=\nabla\psi$ with $\nabla^2\psi=-\dot\sigma$. In spherical harmonics this is diagonal:

$$\psi_{nm}=\frac{R^2\,\dot\sigma_{nm}}{n(n+1)}\quad(n\ge1).$$

$F$ then follows from the first partials that `geopotential` already provides (method M2). Among all fluxes with the observed divergence, $\nabla\psi$ minimises $\int|F|^2$: any divergence-free addition only adds energy, because its cross term integrates to zero on a closed sphere. High-degree stripes are damped by $1/n$, so flux maps are smoother than the anomaly maps behind them.

The same algebra gives a closed-form comparison metric, the H⁻¹ (linearised Wasserstein) distance $\sum R^2|\Delta\sigma_{nm}|^2/(n(n+1))$. Unlike RMS, it does not penalise a slightly shifted blob twice. That makes it a fairer score for COST-G against CSR, GFZ, JPL and ITSG, and for gap fillers.

Related prior work exists but differs. Lateral fluxes have been inferred only through process models with adjoints ([EGU26-4049](https://meetingorganizer.copernicus.org/EGU26/EGU26-4049.html)) or proposed as a DMD network between sub-basins ([Uni Stuttgart](https://www.gis.uni-stuttgart.de/forschung/projekte/data-driven-modelling)). Kusche et al. warn that short-term lateral transport is hard to infer from noisy monthly data ([Kusche et al. 2016](https://hal-insu.archives-ouvertes.fr/insu-01355127/file/Kusche_et_al-2016-Geophysical_Research_Letters.pdf)).

The honest label is therefore "minimum transport consistent with GRACE". It is physically meaningful over the ocean. Over land, where precipitation minus evaporation feeds mass vertically, it is a kinematic diagnostic only. Adding ERA5's P − E as a source term turns it into a mapped water-budget residual (lateral river plus groundwater transport), which we did not find published either.

### A per-coefficient Kalman smoother is an exact spherical Gaussian process

Gap filling is crowded. The yardstick is the 2024 J. Geodesy benchmark, in which a plain autoregressive (AR) model won with errors of **0.9–1.7 cm equivalent water height (EWH)**, while Swarm-only fields scored 7.4–7.6 cm ([Lecomte et al. 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11585520/)). No method in that benchmark returns calibrated uncertainties.

Method M3 runs, for each of the 8,281 coefficients, a structural time-series model:

- **State:** level, slope, and annual and semi-annual terms.
- **Process noise:** a random walk, pooled per degree and fitted by maximum likelihood.
- **Observation noise:** $(s_n\sigma_{nm}(t))^2$, where $s_n$ is a per-degree calibration factor that is also estimated.

A Rauch–Tung–Striebel (RTS) smoother then gives every month, including missing ones, a smoothed value and a posterior variance. The cost is linear in coefficients × months, and each new month is one filter step.

The framing is the publishable part. An isotropic Gaussian process on the sphere is diagonal in SH. A temporal Matérn kernel has an exact state-space form. So M3, with degree-tied hyperparameters, is an exact separable space–time Matérn Gaussian process on the sphere, at linear cost. As a side result, the estimated $s_n$ answers a question the 2026 assimilation review cares about: are the COST-G sigmas calibrated?

Kalman smoothing of SH coefficients exists for ITSG's daily fields, with process noise taken from models ([ITSG README](https://ftp.tugraz.at/pub/ITSG/GRACE/ITSG-Grace2018/daily_kalman/README.txt)). State-space models exist for regional series ([Didova et al. 2016](https://www-app2.gfz.de/grace-references/b2hd-2016JGeod..90.1207D.html)). The monthly, per-coefficient, sigma-calibrated variant is a partial novelty.

Adding Swarm (to about degree 12–15) and SLR (degree ≤ 5–10) monthly fields as extra observations, each with its own sigma, gives a fused estimate across the 2017–18 gap (M4). Swarm is "the backbone" of that gap in existing combinations ([J. Geod. 2024](https://link.springer.com/article/10.1007/s00190-024-01888-5)).

### Five cheaper methods complete the set, and several popular ones should be skipped

| ID | Method | What it adds | Novelty (web-limited) | Validation anchor | Effort |
|---|---|---|---|---|---|
| M1 | SH matched-filter (GLRT) scan with scan-statistic nulls | Blind detection of abrupt mass changes, no destriping | No prior found | Blind recovery of Sumatra 2004, Maule 2010, Tohoku 2011; injection curves | Medium |
| M2 | Minimum-norm transport field plus H⁻¹ metric | "Where water moves" maps; a fairer comparison metric | No prior found | Synthetic moving blob; Greenland outflow points seaward | Low |
| M3 | Per-coefficient state-space smoother = spherical Matérn Gaussian process | Gap filling with posterior variance; sigma calibration | Partial (framing new) | Lecomte remove-restore test; 68/95% coverage | Low–medium |
| M4 | Kalman fusion of COST-G, Swarm and SLR | An honest bridge across 2017–18 | Partial | Withhold 11 months in 2015–16 | Medium |
| M5 | Bayesian online change-point detection (BOCPD) with known noise variance | Streaming alerts as each month arrives | No prior found ([BOCPD](https://arxiv.org/pdf/0710.3742)) | Synthetic steps; mean run length to a false alarm | Low |
| M6 | Diagonal plus low-rank error covariance from residual EOFs | A cheap substitute for VADER-style anisotropic filtering without full covariances | Partial | Compare with ITSG-Grace2018 full covariances | Low–medium |
| M7 | Linear inverse model / principal oscillation patterns, and optimised or bagged (BOP) DMD on whitened EOFs | Forecasts, oscillatory modes, gap-aware decomposition with error bars | LIM: no prior found; DMD on GRACE exists ([Libero et al.](https://meetingorganizer.copernicus.org/EGU24/EGU24-18454.html)) | 1–6-month skill against damped persistence | Low–medium |
| M8 | Needlets with sigma-aware thresholding; region-adaptive Wiener filter from Slepian multitaper spectra | Scale-dependent destriping; filters that are weaker over ice, stronger over deserts | Needlets: no prior found; Wiener: partial ([Sasgen 2006](https://www-app2.gfz-potsdam.de/grace-references/b2hd-2006StGG...50.499S.html)) | Closed-loop leakage against Gaussian, DDK and Wiener filters | Low / medium |
| M9 | Slow feature analysis | Slow drivers (depletion, decadal ocean) without a trend model | No prior found | Synthetic slow driver | Low |
| M10 | PCMCI+ with partial correlation; Liang–Kleeman information flow | Separates direct from mediated links between climate drivers and storage | Application not found ([review](https://arxiv.org/pdf/2404.05746)) | Synthetic networks with known causal graphs | Medium |
| M11 | Sparse (L1) source localisation, FISTA (fast iterative shrinkage-thresholding) | Compact sources (glaciers, reservoirs) that L2 methods blur | Partial ([sparse DDK](https://research.tudelft.nl/files/125307388/remotesensing_14_02810_v2.pdf)) | Injected loads | Medium |

M6 deserves emphasis because it repairs M1's main weakness. The leading EOFs of month-to-month residuals capture the coherent stripe pattern. The Woodbury identity then makes weighted least-squares and matched-filter weights cost O(N·k²) with k ≈ 20–50, instead of handling 549 MB covariance matrices per month. The time-variable decorrelation filter VADER cut cumulative geoid errors by about 15% using formal full covariances ([Horvath et al. 2018](https://gfzpublic.gfz-potsdam.de/pubman/item/item_3766971_3/component/file_3806902/3766971.pdf)). M6 tests whether residuals alone can recover most of that.

Several popular directions should be skipped:

- **Deep-learning downscaling.** The critiques above apply.
- **Neural fields.** Hybrid SH-plus-SIREN models beat linear SH for the static EGM2008 field ([EGU26-14971](https://meetingorganizer.copernicus.org/EGU26/EGU26-14971.html)), but none model time variation, and they add little at degree 90.
- **Nonlinear dependence methods.** Convergent cross mapping, nearest-neighbour transfer entropy and recurrence quantification are fragile at 259 strongly seasonal samples.
- **Koopman forecasting.** It was less robust than linear models with short training records ([EGU24-2554](https://meetingorganizer.copernicus.org/EGU24/EGU24-2554.html)).
- **Deep-learning destriping of SH coefficients.** We found no paper on it, so it is a genuine niche. It belongs after M6, benchmarked against DDK, VADER and recent physics-informed filters.

## A shared protocol stops 24 years of data from fooling us

The climate-link literature for gravity data is strongest on sign and location and weakest on significance.

- **Too few independent samples assumed.** Phillips et al. used a 95% threshold of r = 0.21, which implies about 90 independent months out of 96. The authors themselves called the degrees of freedom "unclear" and warned of "false correlations" ([Phillips et al. 2012](https://escholarship.org/content/qt9p53g2p8/qt9p53g2p8.pdf?t=n0u7rd)). A recent COST-G Greenland study reports **p = 10⁻²⁷⁴** for correlations between cumulative indices and a trending mass series ([Cambier et al. 2025](https://dx.doi.org/10.3390/rs17213552)). Correlating two integrated, trending series is the textbook setting for spurious regression.
- **Weak results reported as significant.** A 2025 claim of a solar-cycle imprint on land water rests on **r = 0.12**, described as "significant at over 80%" ([Sci. Rep. 2025](https://www.nature.com/articles/s41598-025-99880-2)).
- **Methods that cannot show causation.** Wavelet coherence is the dominant method, and even its users note it "cannot establish causality" ([Rateb & Scanlon 2025](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2025AV001684)).

Our own arithmetic shows the scale of the problem. Use the standard AR(1) effective-sample-size correction, N_eff ≈ N(1 − r₁r₂)/(1 + r₁r₂). A deseasonalised storage series with lag-1 autocorrelation 0.9 and an ENSO index at 0.95 give **N_eff ≈ 20 from N ≈ 250**. The 95% critical correlation is then **about 0.44, not 0.12**. These autocorrelations are illustrative; compute them per cell.

Spatial multiplicity compounds this. At 300 km smoothing the land surface holds only about 530 independent footprints. A sweep over ten indices, thirteen lags and every 2° cell runs hundreds of thousands of tests, so maps stippled at a pointwise 5% level will be stippled by chance. False-discovery-rate (FDR) control is recommended because it is only mildly sensitive to spatial correlation ([Wilks 2016](https://ams.confex.com/ams/96Annual/webprogram/Paper276937.html); [Ventura et al. 2004](https://www.stat.berkeley.edu/~paciorek/files/jclim/vent.etal.2004.pdf)).

The protocol for every test follows from these numbers.

**Signal handling**

- **Deseasonalise.** Fit annual, semi-annual, 161-day and 322-day terms, not a monthly climatology. The 161-day alias of the S2 tide persists in the RL06 release over West Antarctica ([Remote Sens. 2022](https://doi.org/10.3390/rs14215403)).
- **Difference or detrend.** Analyse first differences dS/dt or detrended anomalies. This removes linear GIA exactly, which matters because our grids do not remove GIA.
- **Respect the gap.** Never interpolate across 2017-07 to 2018-05.
- **Use true epochs.** Use each file's actual mid-epoch, not its calendar label. The 259 fields cover only 257 calendar labels: 2011-10 and 2015-04 each carry two fields, and 35 calendar months, in 19 gaps, have none (`viewer/data/tws_index.json`).

**Validation and significance**

- **Split by mission.** Fit on GRACE (2002–2017) and test on GRACE-FO (2018–2026). They are independent missions.
- **Null distributions.** Draw at least 1,000 phase-randomised or IAAFT (iterative amplitude-adjusted Fourier transform) surrogates of the index while keeping the gravity field intact, and record the maximum statistic over cells and lags. Apply Benjamini–Hochberg FDR at q = 0.1 to maps. Report effect sizes with 12-month block-bootstrap intervals.
- **Pre-register forecasts.** Commit forecasts as dated git tags before the target data exist.

Three controls run before any discovery is claimed.

- **Q0: noise covariates.** This control matters more than it looks, because several possible confounders overlap the key windows:
  - GRACE-FO's accelerometer "transplant" data affect the C20 and C30 coefficients ([EGU2020-3017](https://meetingorganizer.copernicus.org/EGU2020/EGU2020-3017.html)), and transplants get harder near solar maximum ([EGU24-10401](https://meetingorganizer.copernicus.org/EGU24/EGU24-10401.html)).
  - Thermospheric density rose about **two orders of magnitude** toward the solar-cycle-25 maximum ([EGU25-20153](https://meetingorganizer.copernicus.org/EGU25/EGU25-20153.html)).
  - Wide-pointing mode has run since July 2023.
  - COST-G's contributing centres and weights change from month to month. In July 2026 they were GFZ 0.04, CSR 0.20, JPL 0.06, ITSG 0.12, LUH 0.10, AIUB 0.06, CNES 0.06, HUST 0.19 and APM-SYSU 0.17 ([COST-G series](https://icgem.gfz.de/sp/02_COST-G_/Grace-Grace-FO_RL02.1)).

  All four confounders overlap the 2023–24 El Niño step and the Kamchatka window. The test regresses each month's ocean residual RMS, its high-degree (n > 40) degree variance and the C20/C30 scatter on F10.7 solar flux, the Kp geomagnetic index, transplant and wide-pointing flags, and the centre weights parsed from file headers. It uses autocorrelation-robust (HAC) errors. The output is a per-month noise weight. A null means months can be treated as equal beyond their sigmas. We found no study that quantifies monthly Level-2 noise against F10.7.
- **Q1: negative controls.** Solar-cycle and quasi-biennial oscillation (QBO) indices, plus a deliberately scrambled index, run through every pipeline. Any "detection" flags trend leakage or gap handling, not a discovery.
- **Q2: sigma calibration** (from M3). If s_n ≈ 1 at every degree and withheld months fall inside the 68% and 95% bands at the nominal rates, COST-G's sigmas are trustworthy. Either answer is publishable.

## Twenty-four hypotheses, ranked by payoff, risk and urgency

Each hypothesis below states the data, the test, what a null result would look like, and its novelty status. These searches were limited: one general web index, a shared budget of about 200 queries across six threads (59 for algorithms, about 45 for climate), run on 24 October 2026. "No prior found" means not found by those queries. Run targeted Google Scholar and Web of Science searches for each before investing, and again before claiming priority.

### Weather and climate: the El Niño ratchet leads a dozen tests

**W1, the El Niño ratchet (time-critical).**

- **Background.** Global land storage outside the ice sheets dropped abruptly between May 2014 and March 2016 and never recovered ([Rodell et al. 2024](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11671563/)). It set record lows of **−6,697 km³ in 2023 and −7,404 km³ in 2024** ([NREE 2023](https://www.nature.com/articles/s43017-024-00545-x); [NREE 2024](https://www.nature.com/articles/s43017-025-00659-w)). The barystatic sea-level rise (from added ocean mass) was about **6.2 mm in both the 2014–16 and 2023–24 El Niños**, faster in the second ([Jin et al. 2026](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2025JC023533)). La Niña works the other way: in 2010–11, water stored on land lowered sea level by about 7 mm ([Fasullo et al. 2013](https://climate.fas.harvard.edu/sites/projects.iq.harvard.edu/files/climate/files/fasullo.etal_2013.pdf)).
- **Hypothesis.** The response is asymmetric. Strong El Niños remove more water, faster, than La Niñas restore, which produces a staircase.
- **Data.** Continental land storage with ice sheets and glacier regions masked, and the monthly Niño 3.4 anomaly from NOAA's index list ([NOAA PSL](https://psl.noaa.gov/data/climateindices/list/)).
- **Test.** A distributed-lag regression of dS/dt on positive and negative Niño 3.4 separately (lags 0–9 months, plus a trend term). The null is equal coefficients (β⁺ = β⁻), tested with block bootstrap and surrogates.
- **Pre-registration.** ONI already stands at **+1.39, +1.80 and +2.16** for the three-month seasons May–July, June–August and July–September 2026 ([CPC ONI](https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt)). Fit on 2002–2025, then commit continental anomaly forecasts for December 2026, March 2027 and June 2027 with 90% intervals now. August 2026 fields from CSR and ITSG are already public, so they cannot count as out-of-sample.
- **What a null looks like.** The positive and negative coefficients agree within their intervals, and the 2026–27 drop falls short of the forecast's lower bound or recovers within 12 months. That would make the staircase a trend plus red noise.
- **Novelty.** The steps are documented; we found no test of the asymmetry. Power is limited to two or three strong events, so pre-registration is the protection against fooling ourselves.

**W2, are the steps real change-points?**

- **Background.** Two published break dates disagree: 2011–12 for extremes ([Rateb & Scanlon 2025](https://www.eurekalert.org/news-releases/1112360)) and 2014–16 for mean storage. Neither, as far as we found, was tested against a red-noise or long-memory null.
- **Data.** Global non-ice land storage.
- **Test.** Pruned exact linear time (PELT) or likelihood-ratio change-point detection with trend and seasonal terms. Calibrate critical values against AR(1) and ARFIMA (long-memory) surrogates fitted to 2002–2014. Repeat without north-eastern South America to see whether the step is global, and repeat on GRACE-FO alone for 2023–24, with the Q0 centre-weight covariate included.
- **What a null looks like.** The observed step statistic sits below the 95th percentile of the surrogate maxima.
- **Novelty.** We found no formal test against red noise.

**W3, are atmospheric rivers storage-efficient?**

- **Background.** Atmospheric-river (AR) storms drove California's water-year-2023 gains ([J. Hydrol. 2025](https://www.sciencedirect.com/science/article/abs/pii/S0022169425001428)). East Antarctica gained about 695 Gt in 2021–23 through AR-fed snowfall ([ScienceDaily 2026](https://www.sciencedaily.com/releases/2026/09/260906170128.htm)).
- **Hypothesis.** At equal precipitation anomaly, months with more AR days produce larger dS/dt.
- **Data.** Monthly AR-day counts from the Guan–Waliser ERA5 catalogue ([PANGAEA](https://doi.pangaea.de/10.1594/PANGAEA.957161), extended by the [UCLA v4 release](https://dataverse.ucla.edu/dataverse/ar)); GPCP precipitation; dS/dt.
- **Test.** Cold-season regressions dS/dt = a·P′ + b·AR′ for basins larger than 2×10⁵ km² (Sierra and Central Valley, Pacific Northwest, central Chile, Norway, East Antarctic drainage). Permute years within seasons; cross-validate by leaving one water year out; report partial R² and variance inflation.
- **What a null looks like.** b ≈ 0 with negligible partial R².
- **Novelty.** We found no global AR–storage regression.

**W4, deep memory and heat.**

- **Background.** Storage–temperature coupling is controlled by root-zone moisture ([Chen et al. 2019](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2019JD030324)).
- **Hypothesis.** April–May storage predicts summer hot-day counts beyond what surface soil moisture and precipitation predict.
- **Test.** The cross-validated gain in R² per cell, FDR across cells, with irrigated cells masked.
- **What a null looks like.** No cells survive FDR, even in the deep-rooted transition zones (eastern South America, southern Africa, the Great Plains).
- **Novelty.** We found no test of GRACE's added skill for heat.

**W5, fire beyond the Amazon.**

- **Background.** Storage deficits preceded southern Amazon fire seasons by 1–5 months ([Chen et al. 2013](https://agupubs.onlinelibrary.wiley.com/doi/10.1002/jgrg.20046)).
- **Test.** Spearman regression of log burned area on pre-season storage, controlling for precipitation and ENSO, in boreal Canada and Siberia, south-eastern Australia, the Mediterranean and African savannas. Fit on 2002–2016, test on 2018–2025, and pre-register the 2026–27 Amazon and Australian seasons alongside W1.
- **What a null looks like.** Leave-one-season-out skill no better than precipitation alone.
- **Novelty.** An extension of a single regional study.

**W6–W11, lower-priority climate tests.**

- **W6** replicates the finding that the CO₂ growth rate is sensitive to land water storage ([Humphrey et al. 2018](https://ideas.repec.org/a/nat/nature/v560y2018i7720d10.1038_s41586-018-0424-4.html)). Fit on GRACE, predict GRACE-FO including 2023–24, and run a Chow test. The null is a slope indistinguishable from zero after 2018.
- **W7** asks whether central-Pacific and eastern-Pacific El Niños leave distinct global storage fingerprints. Only regional evidence exists ([Nicaragua](https://ui.adsabs.harvard.edu/abs/2022RemS...14.6012J); [Yangtze](https://doaj.org/article/e69a08b36ab74379b4974be61e21a023)). The test compares the two regression maps by pattern correlation against surrogates, and a null means the maps differ no more than chance. It is exploratory given the few events.
- **W8** maps each cell's "memory" τ as a leaky integrator of its dominant index, TWS(t) ≈ Σₖ β e^(−k/τ) I(t − k). The null is no cross-validated gain over instantaneous regression. Hydraulic memory is known for precipitation forcing, not index forcing, so this is a partial novelty, and τ sets the maximum useful forecast lead for W1, W4 and W5.
- **W9** replicates the warm-pool and atmospheric-river explanation for East Antarctica. Fit 2002–2020 and predict 2021–23 from basin dM/dt. The null is a failed prediction.
- **W10** extends a US finding that pre-season storage changes crop-yield sensitivity to rain ([PMC6850578](https://pmc.ncbi.nlm.nih.gov/articles/PMC6850578/)) to global panels.
- **W11** tests whether tropical land storage helps predict ENSO 3–9 months ahead. A null is expected, and would usefully bound the strength of land–atmosphere feedback.

### Solid Earth: Kamchatka 2025 is the most time-critical target

**E1, the Kamchatka Mw 8.8 earthquake (time-critical).**

- **Background.** GRACE-FO monthly fields show a coseismic signal of "nearly 20 µGal" peak-to-peak from the 29 July 2025 earthquake. The only published work is a 2026 conference abstract that proposes reading any excess over the seismic moment as aseismic slip ([EGU26-15981](https://meetingorganizer.copernicus.org/EGU26/EGU26-15981.html)).
- **Expected pattern.** Earlier great earthquakes followed a two-stage pattern: gravity kept falling for months after the event, then recovered ([Heki group](https://eprints.lib.hokudai.ac.jp/dspace/handle/2115/67855)).
- **Data.** COST-G to 2026-07 (12 full post-event months); the USGS finite-fault model; aftershocks from the USGS catalogue ([FDSN](https://earthquake.usgs.gov/fdsnws/event/1/)).
- **Test.** Treat July 2025 as a mixed month (only about 8% of it follows the event). Fit seasonal and alias terms on 2018 to June 2025. Fit a step plus an exponential or logarithmic recovery within about 2,500 km, at n ≤ 60 without destriping (M1). Pass a dislocation forward model through the same filter, and take null amplitudes from the same months at 100 random ocean locations, weighted by the Q0 noise weights.
- **What a null looks like.** A step below the 95th percentile of the random-location null. That would mean a Mw 8.8 event is undetectable in COST-G months during solar maximum and wide-pointing mode, which is itself a useful threshold.
- **Novelty.** Only an abstract precedes us, so this is the most likely to be published by someone else soon.

**E2, a blind global detector with a stated false-alarm rate.**

- **Background.** About ten megathrust events since 2002 have known GRACE signals ([AGU FM15](https://agu.confex.com/agu/fm15/webprogram/Paper85405.html)). Chignik 2021 (Mw 8.2) sits "near the coseismic detection threshold" ([ESSOAr](https://essopenarchive.org/doi/full/10.1002/essoar.10510424.1)).
- **Hypothesis.** The M1 scan recovers every shallow Mw ≥ 8.5 event blind, and detection probability falls below 50% near Mw 8.2–8.3.
- **Data.** COST-G, the USGS catalogue at Mw ≥ 7.5, and Global CMT focal mechanisms.
- **Test.** A ROC curve, a logistic fit of detection probability against Mw, and a null from 1,000 time-shuffles. Hydrological extremes count as "true but non-seismic" detections, to be classified rather than scored as false alarms.
- **What a null looks like.** No operating point recovers the Mw ≥ 8.5 events at an acceptable false-alarm rate.
- **Novelty.** We found no published global change-point detector for GRACE with a false-alarm rate.

**E3, falsifying the pre-seismic signal.**

- **Background.** A gravity-gradient anomaly before Tohoku was reported by Panet et al. ([2018](https://www-app2.gfz.de/grace-references/b2hd-2018NatGe..11..367P.html)). It was judged not unique in time or space by an independent permutation test ([Wang & Bürgmann 2019](https://www-app2.gfz-potsdam.de/grace-references/b2hd-2019GeoRL..46.7323W.html)), then defended with cross-model consistency ([Panet et al. 2022](https://www.documentation.ird.fr/hor/fdi:010085978)). A 2026 preprint on Turkey 2023 reports no significance tests at all ([arXiv 2606.17938](https://arxiv.org/abs/2606.17938)).
- **Test.** Pre-register a 6-month pre-event window and a threshold for Tohoku 2011, Maule 2010, Illapel 2015, Chignik 2021 and Kamchatka 2025. Compare against random times and places along the same subduction zones, and require consistency between COST-G and individual-centre solutions.
- **What a null looks like (the expected outcome).** Pre-event anomalies are no more frequent than at random times and places.
- **Novelty.** A test of a contested claim. We found none run on a combined product.

**E4, hydrological loading and seismicity worldwide.**

- **Background.** Seasonal water loading modulates small-earthquake rates at the kilopascal level in California ([Johnson et al. 2017](https://central.scec.org/node/5455)), Nepal (winter rates about twice summer rates; [Bettinelli](https://pubs.usgs.gov/of/2010/1099/chanard/of2010-1099_chanard.pdf)), the Himalayan thrust at about 2 kPa ([Panda et al. 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6175945/)) and New Madrid ([Craig et al. 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5732252/)). In Taiwan, shallow and deep events are out of phase ([Hsu et al. 2021](https://escholarship.org/content/qt0h27b1fx/qt0h27b1fx.pdf)), so a single global statistic could cancel real signals.
- **Test.** Compute GRACE-derived load Coulomb stress on Global CMT fault planes for M ≥ 4.5 events, declustered. Run Schuster periodicity tests per region, stratified by depth and faulting type, with FDR. California is the positive control.
- **What a null looks like.** Uniform p-values after FDR while the California control is still recovered.
- **Novelty.** We found no uniform global GRACE-based test.

**E5 and E6, human-made loads.**

- **E5** asks whether groundwater depletion changes long-term earthquake rates. It tests single-case claims statistically: Lorca ([González et al.](https://www.earth-prints.org/handle/2122/8270)), the Himalaya ([Kundu et al.](https://www-app2.gfz.de/grace-references/b2hd-2015GeoRL..4210607K.html)) and Morocco ([SSA](https://www.seismosoc.org/news/groundwater-depletion-may-have-been-significant-trigger-for-2023-morocco-earthquake/)). The design is a difference-in-differences comparison of declustered rates between depletion regions and paired controls. The null is no relative rate change.
- **E6** uses the Grand Ethiopian Renaissance Dam (GERD), a **74 km³** reservoir filled in yearly steps from 2020 to 2024 ([GSTM2022-68](https://meetingorganizer.copernicus.org/GSTM2022/GSTM2022-68.html)), as a ground truth of known timing and location. Point-mass forward models pass through the same filter. The output is an empirical minimum detectable mass for COST-G, and also an anthropogenic positive for E2. The null is no detection even cumulatively, which would set a floor of several tens of km³.

### Earth rotation and the core: close the budget, then look at the residual

**R1, closing the polar-motion and length-of-day budget.**

- **Background.** Surface mass transport explains **83 ± 23%** of 2003–2015 polar motion ([Adhikari & Ivins 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4846461)). Climate-driven mass change lengthened the day by **1.33 ms per century** over 2000–2018 ([JPL](https://www.jpl.nasa.gov/news/nasa-funded-studies-explain-how-climate-is-changing-earths-rotation/)). The Chandler wobble nearly vanished after 2015 and then re-excited ([Shi et al. 2025](https://link.springer.com/article/10.1007/s00190-025-02021-w)).
- **Data.** COST-G ΔC21 and ΔS21 (polar motion) and the SLR-derived C20 from TN-14 (length of day); degree 1 from TN-13; the atmosphere and ocean mass removed by de-aliasing, restored from the GAC product; model motion terms from GFZ's effective angular momentum (EAM) series; IERS EOP C04 observations.
- **Test.** Compare excitation functions in complex band-pass, with coherence and residual spectra at the annual and Chandler periods.
- **What a null looks like.** A residual that stays coherent in the Chandler band after hydrology and cryosphere are included. That would point to an unmodelled excitation.
- **Novelty.** An extension to 2026 with an independent combined product.

**R2, GRACE-calibrated hydrology for rotation forecasts.** GRACE's latency makes it useless for direct operational forecasting. It could instead bias-correct the model hydrological angular momentum inside forecasts that already reach **3.02 mas** for 10-day polar motion ([GFZ](https://gfzpublic.gfz-potsdam.de/pubman/item/item_3247891_5/component/file_4001888/3247891.pdf)). The test is a hindcast skill comparison, and the null is no gain. We found no study of this use.

**R3, replicating the January 2007 core anomaly (high risk).**

- **Background.** A conference abstract reports a rapid north–south gravity anomaly over the eastern Atlantic in January 2007, coincident with a geomagnetic jerk ([EGU26-12420](https://meetingorganizer.copernicus.org/EGU26/EGU26-12420.html)). It follows earlier claims of core signals ([Mandea et al. 2012](https://www-app2.gfz-potsdam.de/grace-references/b2hd-2012PNAS..109..19129M.html)). An independent expert notes that the observed signal is at least 10× what the core should produce ([Eos](https://eos.org/articles/can-nasas-gravity-satellites-detect-motions-in-earths-core)). The same group found that solutions disperse by about **0.34 cm EWH**, or roughly 20% of the low-degree signal ([Lecomte et al. 2023](https://www-app2.gfz-potsdam.de/grace-references/b2hd-2023JGRB..12826503L.html)).
- **Test.** Run M1 on COST-G, a processing chain independent of the GRGS solution used in the claim, after removing the spread of a hydrology-and-ocean model ensemble. Cross-correlate any step with geomagnetic secular acceleration and with the R1 residual.
- **What a null looks like (the likely outcome).** No step outside the model spread. That result is still informative because it comes from an independent product.

**R4, the inner-core rotation bound.** Updating the inner-core rotation bounds from the S22 coefficient (0.09°, 0.3° and 0.4° at about 4, 6 and 12-year periods; [Lecomte et al. 2023](https://doaj.org/article/ba63cd986f8941a898d65fc646885c7d)) with 2023–26 data is cheap but should tighten them only marginally.

### Ranked register

| Rank | ID | Hypothesis | Key data beyond COST-G | Decisive test | Null looks like | Novelty (web-limited) |
|---|---|---|---|---|---|---|
| 1 | Q0 | Noise depends on solar activity, transplant periods, wide-pointing mode and centre weights | F10.7, Kp, file-header weights | HAC regression of ocean residual RMS | No covariate matters | No prior found |
| 2 | W1 | El Niño ratchet, pre-registered for 2026–27 | Niño 3.4, ONI | Asymmetric distributed-lag regression; dated forecast | β⁺ = β⁻; drop recovers within 12 months | Asymmetry untested |
| 3 | E1 | Kamchatka 2025 co- and postseismic signal, aseismic excess | USGS finite fault | Step plus recovery fit against a 100-location null | Below the 95th percentile of the null | Abstract-only prior |
| 4 | E2 | Blind detector with Mw threshold and false-alarm rate | USGS catalogue, Global CMT | ROC; logistic P(detect) against Mw | No useful operating point | No prior found |
| 5 | W2 | 2014–16 and 2023–24 steps are true change-points | — | PELT against AR(1) and ARFIMA surrogates | Inside the surrogate range | No formal test found |
| 6 | R1 | Rotation budget closes through 2026, including the Chandler anomaly | IERS C04, EAM, GAC | Band-pass excitation comparison | Coherent Chandler residual | Extension |
| 7 | W3 | Atmospheric-river precipitation is storage-efficient | AR catalogue, GPCP | Partial regression, leave-one-year-out | b ≈ 0 | No prior found |
| 8 | W4 | Pre-summer storage predicts hot days | ESA CCI soil moisture, Tmax | Cross-validated R² gain with FDR | No cells survive | No prior found |
| 9 | E3 | Pre-seismic anomalies are not special | Catalogue | Pre-registered permutation test | No difference from random (expected) | Contested-claim test |
| 10 | W8 | Storage memory map τ | Indices | Leaky integrator against instantaneous model | No cross-validated gain | Partial |
| 11 | W5 | Pre-season storage predicts burned area beyond the Amazon | GFED5 or MODIS burned area | Out-of-sample rank regression | No skill beyond precipitation | Extension |
| 12 | E4 | Seasonal loading phase-locks seismicity, stratified | Global CMT | Schuster test per region, FDR | Uniform p-values; control passes | No global test found |
| 13 | Q2 | COST-G sigmas are calibrated | — | Estimate s_n; coverage on withheld months | s_n ≈ 1 (also useful) | Not checked |
| 14 | E6 | GERD sets the minimum detectable mass | Altimetry levels | Matched forward model | No detection | Partial |
| 15 | W7 | ENSO diversity fingerprints | Niño 3/Niño 4 | Map difference against surrogates | Maps indistinguishable | Regional only |
| 16 | W6 | CO₂ growth sensitivity holds after 2018 | NOAA GML | Chow test | Slope ≈ 0 | Replication |
| 17 | W9 | East Antarctic warm-pool and AR link | ERSST, AR counts | 2002–20 fit predicts 2021–23 | Prediction fails | Replication |
| 18 | R2 | GRACE-calibrated hydrology improves rotation forecasts | EAM forecasts | Hindcast skill | No gain | No study found |
| 19 | E5 | Depletion changes long-term earthquake rates | Regional catalogues | Difference-in-differences | No relative change | No statistical test found |
| 20 | R3 | January 2007 core-linked step | CHAOS-7, model ensemble | M1 step against model spread | No step (likely) | Conference-only claim |
| 21 | W10 | Pre-season storage modulates yield response | GDHY or FAOSTAT yields | Panel regression | No interaction | US-only prior |
| 22 | R4 | Tighter inner-core S22 bound | — | Lecomte method, multi-model hydrology | Marginal change (expected) | Update |
| 23 | W11 | Land storage helps predict ENSO | Warm water volume | PCMCI+ or conditional Granger | Null (expected) | No prior found |
| 24 | Q1 | Solar and QBO imprints | Sunspots, QBO | Same pipeline as all tests | Null required | Negative control |

## Animations should show change, gaps and doubt together

The evidence on animation is clear and slightly uncomfortable:

- **Animation is for presenting, not analysing.** Animation is the fastest way to present and the most enjoyable, but it causes many errors; for analysis, static small multiples were most accurate ([Robertson et al. 2008](https://faculty.cc.gatech.edu/~stasko/papers/infovis08-anim.pdf)). On phones, small multiples beat animation in **7 of 9 tasks** ([arXiv:1907.03919](https://arxiv.org/pdf/1907.03919)).
- **Animated uncertainty works.** Hypothetical outcome plots (HOPs), which animate individual draws from the uncertainty distribution, let viewers compare quantities more accurately than error bars ([Hullman lab](https://www.hullmanlab.northwestern.edu/paper%20discrete-outcome/2015/01/01/hops.html)). They also help untrained viewers judge trends when evidence is weak ([Kale et al. 2018](https://mucollective.northwestern.edu/files/2018-HOPsTrends-InfoVis.pdf)).
- **Colour maps distort.** Non-uniform colour maps such as jet can distort what viewers see by **more than 7%**. Colour-vision deficiency affects 8% of men ([Crameri et al. 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7595127)).
- **The reference set has a gap.** NASA's reference storage animations use one variable, a two-hue dry/wet colour key, a fixed projection and 30 fps ([SVS 4415](https://svs.gsfc.nasa.gov/4415)). None of the pages we fetched documents how they interpolate between months or shows uncertainty. That is the gap our viewer can fill honestly.

| Rank | Animation | What it shows | Inputs from our crates | Honesty guard |
|---|---|---|---|---|
| 1 | Gap-aware monthly playback with a small-multiples calendar | 24 years of anomalies; tap a thumbnail to jump | Existing grids as binary float16 | Interpolate data, not colours; no tween across missing months; time-proportional playhead |
| 2 | Seasonal "breathing" map | Annual phase (cyclic romaO colours) × amplitude (lightness); a synthetic year loop shows seasonal waves sweeping | Per-cell harmonic fit | Caption the smoothing bias: amplitudes drop 25–40% at 800 km smoothing ([CU](https://sealevel.colorado.edu/index.php/bibcite/reference/302)) |
| 3 | Noise-floor mask plus HOP mode | Which features exceed noise; noise "twinkles", robust features persist | Sigmas, M3 posterior, M6 draws | Hatch cells where the anomaly is below σ; 2–4 draws per second |
| 4 | Pre-registered forecast reveal | W1's dated cone, then each new month landing inside or outside it | Committed forecast, new fields | The cone's git date shown on screen |
| 5 | Apparent transport particles, ocean first | M2 flux as fading trails; land shown as pulsing sources and sinks | Poisson-solved ψ and gradients on a 2° grid | Label "minimum transport consistent with GRACE"; state the speed scale |
| 6 | Detection-map playback | M1 z-scores through time, with earthquakes and GERD fills overlaid | Scan output, catalogues | Show family-wise thresholds, not raw p-values |
| 7 | Spherical-harmonic build-up | Partial sums to degree L with a cursor on the degree spectrum and a πR/L km readout | Per-degree grids | Explain ringing; then animate the Gaussian filter weights |
| 8 | Coefficient-triangle playback | Stripes as fixed-order structure; destriping on/off | Monthly ΔC/ΔS to degree 60 | Log-scaled diverging colours |
| 9 | Ground track, aliasing and error-spectrum race | Coverage filling in; a fast sine aliasing to slow; sensor versus aliasing error curves | `Pair::positions`, closed-loop runs, noise CSVs | Carries the message that temporal aliasing, not sensor noise, limits missions ([Zingerle et al. 2024](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024EA003630)) |
| 10 | Rotating exaggerated geoid | The "potato", with an exaggeration slider | Static geoid grid | "×k" always on screen; ESA's video used ×7000 ([Sci. Am.](https://www.scientificamerican.com/blog/observations/the-geoid-why-a-map-of-earths-gravity-yields-a-potato-shaped-planet/)) |

The first five reshape the existing viewer. Idea 5 may be new: we found no particle or flux animation derived from GRACE, though our search ran out before a targeted query. Idea 4 turns scientific honesty into the narrative itself.

**Honesty rules**

- Interpolate in millimetres inside the shader, then apply the colour map.
- Never ease between months, because easing invents sub-monthly dynamics.
- Fix colour limits during playback.
- Show the 2017–18 gap as a hatched "GRACE → GRACE-FO" frame and mark missing months on the slider.
- Render nearest-cell at full zoom so the 2° grid and 300 km smoothing are not dressed up as detail.
- Burn the date, units, colour bar and "COST-G RL02.1, CC BY 4.0" into every exported clip.

**Colour**

- Use Crameri vik for anomalies and romaO for phase ([Crameri cyclic maps](https://www.fabiocrameri.ch/cycliccolourmaps/)).
- Smooth the lightness reversal at the centre of the diverging map, as Kovesi recommends ([Kovesi 2015](https://arxiv.org/abs/1509.03700)).

**Accessibility**

- Under `prefers-reduced-motion`, start paused and replace rotation and scaling with fades ([MDN](https://developer.mozilla.org/en-US/docs/Web/CSS/@media/prefers-reduced-motion)).
- Every auto-playing animation needs a pause control to meet WCAG 2.2.2 ([W3C](https://www.w3.org/WAI/WCAG22/Understanding/pause-stop-hide.html)).
- Keep pulses below three per second.

**Technology**

Hand-written WebGL2 beats every library for a 180×90 field:

- WebGL2 reaches 97.13% of users, against 87.82% for WebGPU ([caniuse](https://caniuse.com/webgl2); [caniuse](https://caniuse.com/webgpu)).
- Library costs are high for one bitmap: deck.gl is 575 KB gzipped (and its globe view is experimental), CesiumJS is 1.77 MB, and MapLibre 6 ships as modules only ([GlobeView](https://deck.gl/docs/api-reference/core/globe-view); [MapLibre v6](https://github.com/maplibre/maplibre-gl-js/releases/tag/v6.0.0)).

Upload all months as one R16F texture array, about 8.4 MB on the GPU. Ship binary float16 instead of today's 5.7 MB base64 JSON. Port d3's Equal Earth inverse to the shader ([d3-geo](https://github.com/d3/d3-geo/blob/main/src/projection/equalEarth.js)). Bilinear-filter by hand, because linear filtering of 32-bit float textures works on only 46–48% of iOS devices ([web3dsurvey](https://web3dsurvey.com/webgl2/extensions/OES_texture_float_linear)).

For export, render frames off-screen and encode them with Mediabunny's `CanvasSource` or CCapture.js v2 ([Mediabunny](https://mediabunny.dev/guide/media-sources); [CCapture](https://cdn.jsdelivr.net/npm/ccapture.js@2.0.0/README.md)), with a gifenc fallback. mp4-muxer is deprecated ([README](https://cdn.jsdelivr.net/npm/mp4-muxer@5.2.2/README.md)). The Rust side can write GIF and APNG files for reports through the `gif` and `png` crates.

## A six-stage build plan pairs every algorithm with a test

The plan follows three rules:

- **Order by deadline.** Do the time-critical commitments first.
- **Build the missing harmonic infrastructure next,** because four methods depend on it.
- **Deliver four things per stage:** one crate feature, one figure, one viewer improvement and one research result, each closed by a gate that checks against independent references rather than our own earlier output.

Generic methods (smoothers, detectors, transport, basin kernels) belong in the published crates. Hypothesis pipelines, readers for external data, pre-registration files and surrogate tests belong in an unpublished `apps/lab` binary next to `gravviz`. That keeps the 1.85-MSRV libraries free of climate-index plumbing.

| Stage | Crates | gravviz figures | Viewer | Research output | Gate |
|---|---|---|---|---|---|
| 0. Commit and control (now, before the forecast target months exist anywhere) | `gravimetry`: seasonal fit with annual, semi-annual, 161-day and 322-day terms; trend with σ. `apps/lab`: N_eff, block bootstrap, phase-randomised and IAAFT surrogates, Benjamini–Hochberg FDR, NOAA PSL/CPC index readers that record download dates; COST-G header weight parser | Global land storage staircase with the dated W1 cone; noise against F10.7 | Real mid-epoch time axis | W1 and W5 forecasts as a dated git tag; Q0 noise weights; Q1 controls | Negative controls null; injected synthetic steps recovered at the nominal rate |
| 1. Harmonic infrastructure and honest playback | `geopotential`: Gauss–Legendre analysis, rotation (FastTransforms port), batched ring synthesis. `gravimetry`: cap and disc loads, basin kernels, transport ψ and H⁻¹ (M2) | Line-integral-convolution (LIC) transport map; H⁻¹ matrix of COST-G against CSR, GFZ, JPL and ITSG | Binary float16 data; WebGL2 renderer; gap-aware clock; small multiples; breathing map; polygon mean with CSV download and a warning below 10⁵ km² | M2 descriptive figure and short paper | Analysis round trip ≈ 10⁻¹⁴; synthetic blob flux recovered; Greenland flux points seaward; viewer pixels match current output |
| 2. Smoother, streaming and errors | `gravimetry`: per-coefficient smoother (M3), BOCPD (M5), diagonal-plus-low-rank covariance (M6), then Swarm and SLR fusion (M4); CF NetCDF writer with σ variables and provenance attributes | Remove-restore bars against the AR benchmark; s_n by degree | Noise-floor hatching; HOP mode; ±1σ ribbons | Q2; gap-filling benchmark | Coverage near 68/95%; RMS reported honestly against 0.9–1.7 cm |
| 3. Detection | `gravimetry::detect`: GLRT scan (M1) with scan-statistic nulls, using M6 noise; dislocation templates via rotation | Detection probability against Mw; Kamchatka map and series | Detection-map playback with catalogue overlay | E1, E2, E3, E6 | Blind recovery of Sumatra, Maule and Tohoku; false-alarm rate measured on GRACE-FO-only years |
| 4. Climate and rotation | `apps/lab` readers for GPCP, AR catalogue, burned area, IERS C04 and EAM; LIM and DMD (M7); PCMCI+ and Liang–Kleeman screening (M10) | Ratchet regression; change-point surrogate histogram; τ map; excitation budget | Anomaly pulses; trend and time-of-emergence map; scrollytelling | W2–W8, R1, E4 | Fit on GRACE, pass on GRACE-FO; FDR q = 0.1; Q1 still null |
| 5. Products and reach | Reading CSR, GFZ, JPL and ITSG monthly, ITSG daily and CNES 10-day fields; GIA ensemble with spread; COG/CSV output; PyO3 abi3 wheels; `geopotential` compiled to WASM. `satgrav`: AOe07 aliasing layer, mission error emulator in ESA-SING format, body parameter. geo5d: clock-site differences with commission, omission and time-variable terms | Error-spectrum race; ground-track coverage; aliasing explainer | In-browser filter radius; SH build-up; spread across solutions; clip export | R2, R3, R4, W9; emulator release for hydrologists | Wheels install on Windows, macOS ARM and Linux; xarray round trip; Zingerle and CARIOQA checkpoints still pass |

Stage 0 is small but has a hard deadline. Its value disappears once the target months exist anywhere. The remaining stages follow from dependencies:

- **Stage 1 before Stage 3.** Basin kernels, the H⁻¹ metric and dislocation templates all need SH analysis and rotation.
- **Stage 2 before Stage 3.** The detector needs both the M3 residual stream and M6 correlated noise to set honest thresholds.
- **User-facing wins ride along early.** Polygon averages with CSV download (Stage 1) and CF NetCDF with uncertainty and provenance (Stage 2) address the top two user requests long before Stage 5's multi-solution reading and Python wheels.

Engineering cautions:

- **Builds.** Run the first WASM and maturin builds, and anything pulling C-backed NetCDF, under a resource cap (`systemd-run --user --scope -p TasksMax=300 -p MemoryMax=8G -- timeout 600 …`), as the project's machine rules require.
- **Licences.** Port only from MIT, BSD or Apache sources. Use GROOPS and ducc0 only as black-box reference implementations. Check the ITSG covariance licence before adding any of it to test fixtures.
- **Index versions.** Record each climate index's version and download date, because CPC revises the Dipole Mode Index for up to two months and ONI uses sliding base periods ([NOAA PSL](https://psl.noaa.gov/data/climateindices/list/)).

## Conclusion

The novel opportunity here is not a new satellite product. It is turning the long gravity record into a disciplined experimental instrument. Three cheap, harmonic-native engines carry most of the novel work. The matched-filter scan serves the Kamchatka, detector, precursor, dam and core tests. The transport field serves the flux maps, the comparison metric and the particle animation. The calibrated smoother serves gap filling, uncertainty, streaming alerts and the answer to whether COST-G's sigmas tell the truth.

Several hidden confounders turned out to sit on top of the most exciting windows: solar maximum, wide-pointing mode, accelerometer transplants and shifting centre weights. That makes the noise control a scientific result in its own right rather than housekeeping. Several of the best tests are expected to return nulls: the precursor, the core step, solar cycles and land-to-ENSO feedback. Published carefully, those nulls would correct a literature that has overstated significance for a decade.

The windows are also closing. GRACE-FO's end date is uncertain, GRACE-C launches by mid-2029, and the 2026–27 El Niño is happening now. So the project's lasting product is a pipeline that ingests data at any cadence, updates its forecasts and detectors with each new field, and states uncertainty and novelty status plainly. That pipeline will carry over unchanged to GRACE-C quicklook fields and, eventually, to the first quantum-sensor solutions.
