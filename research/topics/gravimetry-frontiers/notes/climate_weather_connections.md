# Satellite-gravimetry mass change and weather/climate: established links, contested claims, gaps, and testable hypotheses

Scope note for the report writer: the research was done on 2026-10-24 using about 45 web searches and fetches. The session's web-search budget ran out before every lead could be followed, so "not found" below means "not found in this search", not "does not exist". Labels used throughout:
**[Established]** = several independent studies agree; **[Single-study]** = one paper or group; **[Contested]** = conflicting results or a published critique or retraction; **[Gap]** = nothing found.
Our own data, for reference: 259 COST-G monthly fields (2002-04 to 2026-07, degree 90), processed to 2° TWS grids, 300 km smoothing, GIA **not** removed.

---

## 1. Teleconnections: TWS, ocean bottom pressure (OBP) and ice mass vs ENSO, IOD, PDO, AMO, NAO, SAM, MJO, QBO, solar cycle

### Takeaway
ENSO is by far the best-documented driver of interannual GRACE TWS. It is strongest in the tropics and subtropics (Amazon, Orinoco, La Plata, Borneo, Australia, southern Africa), with lags of 0 to about 6 months, and it synchronises wet and dry extremes across continents. IOD (East Africa, Australia), PDO/IPO (decadal modulation) and AMO are documented regionally. NAO and Greenland blocking are documented for Greenland ice mass, and SAM for Southern Ocean OBP. MJO, QBO and the solar cycle are essentially unstudied with GRACE, or rest on weak single claims. Most studies use pointwise correlation, EOF/ICA or wavelet coherence with weak significance testing. Causal-discovery methods have hardly been applied.

### Cited Findings
**ENSO and land TWS**
- [Established] Phillips et al. 2012 (GRL, doi:10.1029/2012GL052495) used CSR RL04, 2003–2010, against the Multivariate ENSO Index. The peak correlation gives R² = 0.65 for the Amazon and Borneo, rising to 0.76 when a Hilbert-transform term is added to capture lag. Tropical/ITCZ regions correlate negatively with MEI (Venezuela/Colombia 0.62, Costa Rica 0.55, Indian subcontinent 0.45, Borneo 0.65) and arid or temperate regions positively (SE United States 0.44). The largest lags are along the Amazon, Orange, Niger and Mekong rivers. — [Phillips et al. 2012, GRL](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2012GL052495); full text read via [eScholarship PDF](https://escholarship.org/content/qt9p53g2p8/qt9p53g2p8.pdf?t=n0u7rd)
- The same paper's significance threshold was r = 0.21 at 95% (0.26 with lags up to 3 months), which implies roughly 90 independent months out of 96. The authors themselves say "the significant level is questionable due to the shortness of the timeseries and the unclear degree of freedom", and that "there are likely to be false correlations in some regions". Its suggestion that some Greenland mass loss is ENSO-related rests on 8 years containing a strong trend and three La Niñas. — [Phillips et al. 2012 PDF](https://escholarship.org/content/qt9p53g2p8/qt9p53g2p8.pdf?t=n0u7rd)
- [Established] Ni et al. 2018 (Surveys in Geophysics 39:1–22, doi:10.1007/s10712-017-9421-7) found interannual TWS strongly correlated with ENSO over much of the globe, with maximum cross-correlation up to about 0.70. I did not find the specific lags. — [GFZ GRACE bibliography](https://www-app2.gfz.de/grace-references/b2hd-2018SGeo...39....1N.html); [Sci. Rep. 2021 summary of Ni et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8245421/)
- A global GRACE + land-surface-model analysis (J. Hydrology 2021, "Links between global terrestrial water storage and large-scale modes of climatic variability") found the strongest ENSO correlations in the Amazon, Orinoco and La Plata basins, plus significant mid- and high-latitude correlations. — [ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0022169421004662)
- A review in Surveys in Geophysics (2022) states that most interannual variability in GRACE TWS can be related to climate indices such as ENSO. — [Springer review](https://link.springer.com/article/10.1007/s10712-022-09754-9)
- [Single-study, recent] Rateb & Scanlon 2025 (AGU Advances, doi:10.1029/2025AV001684) used GRACE/GRACE-FO 2002–2024, defining extremes as above the 90th or below the 10th percentile, and reconstructed the 11-month 2017–2018 gap probabilistically. They found ENSO the dominant driver of global TWS extremes, synchronising droughts and pluvials across continents through a 2–3-year oscillatory cycle. The wet-tail extremal pattern index is coherent with ENSO indices in the 2–4-year band, and the authors note that "wavelet coherence alone cannot establish causality". Examples: El Niño dry extremes in South Africa (mid-2000s) and the Amazon (2015–16), and La Niña wet extremes in Australia, SE Brazil and South Africa (2010–11). — [AGU Advances](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2025AV001684); [EurekAlert release](https://www.eurekalert.org/news-releases/1112360)
- [Single-study] Rateb & Scanlon 2025 also report a regime change around 2011–2012: wet extremes dominated before 2011 and dry extremes after 2012. They attribute it to a decadal Pacific pattern modulating ENSO, but the press text does not name it. — [EurekAlert](https://www.eurekalert.org/news-releases/1112360)
- [Single-study] Yangtze basin (Earth and Space Science, 2023): eastern-Pacific-type ENSO correlates with TWS at r = 0.56–0.65 with 5–6-month lags. — [DOAJ record](https://doaj.org/article/e69a08b36ab74379b4974be61e21a023)
- [Single-study] Nicaragua (Remote Sensing 2022): the TWS–ENSO resonance periods of about 2 and 4 years are tied to the central-Pacific (quasi-2-year) and eastern-Pacific (quasi-4-year) modes. — [ADS](https://ui.adsabs.harvard.edu/abs/2022RemS...14.6012J)

**IOD, PDO, AMO, NAO, AO and multi-index studies**
- [Established regionally] Africa, using independent component analysis of GRACE: the first independent component of the climate indices is linked to NAO and to significant TWS reductions over southern Africa, while IOD is linked to TWS increases over equatorial East Africa. — [STOTEN 2018 / PubMed 29710593](https://pubmed.ncbi.nlm.nih.gov/29710593/)
- China: ENSO, IOD and PDO significantly affected TWS over more than 80% of regions, and NAO over about 40% (Sustainability 2019, doi:10.3390/su11236646). — [MDPI](https://doi.org/10.3390/su11236646)
- Yunnan: AO and ENSO influence TWS more strongly than PDO and sunspots. The authors note that sunspot effects on TWS are "under-examined". — [J. Hydrology 2019](https://www.sciencedirect.com/science/article/abs/pii/S0022169419304330)
- Global, six ocean-related teleconnections (Water Resources Management 2022, Li, Zha, Shi, Zhong, Tso, Wu): ENSO and AMO are the two most dominant. — [NERC open archive](https://nora.nerc.ac.uk/id/eprint/531665/)
- An EGU 2020 abstract (Tiewei Li, Beijing Normal University) analysed GRACE TWS against 14 teleconnections for 2003–2016. I did not confirm a journal version. — [EGU2020-6613](https://meetingorganizer.copernicus.org/EGU2020/EGU2020-6613.html)
- Sudd wetland (century-scale reconstruction): wavelet coherence with ENSO and IOD at 8–16 months and 4–7 years, with IOD apparently amplifying ENSO. — [ScienceDirect 2025](https://www.sciencedirect.com/science/article/pii/S2666017225000380)
- African aquifers: GRACE storage, hydrologic extremes and teleconnections (Environmental Research Letters 2021). — [IOPscience](https://iopscience.iop.org/article/10.1088/1748-9326/ac3bfc)
- A Climate Dynamics paper ("Analysis of the interannual variability in satellite gravity solutions: detection of climate modes fingerprints…") found that PDO, ENSO, NPGO and AMO explain 30–50% of interannual ice-mass variability in Greenland and West Antarctica, while AO and SAM have a negligible effect. — [Climate Dynamics, doi:10.1007/s00382-021-05953-z](https://link.springer.com/article/10.1007/s00382-021-05953-z) (abstract-level only; the paywall blocked the full text)

**East Africa / Lake Victoria and the Caspian**
- [Established] Lake Victoria rose about 1.4 m in 2019–2020, attributed mainly to successive positive-IOD rainfall, as monitored by GRACE-FO and ERA5 (Sensors 2021, doi:10.3390/s21134304). The same summary quotes "approximately 1750 gigatonne". **Flag:** 1.4 m over the lake's roughly 68,800 km² is only about 96 km³, so that figure probably refers to a different quantity (basin TWS?) and needs checking in the paper. — [Sensors 2021](https://doi.org/10.3390/s21134304); the lake peaked at 13.42 m in May 2020, above 13.41 m in 1964 — [ArcGIS story map](https://storymaps.arcgis.com/stories/bd820937c06845faa86f7f8944d56f47)
- The record 2020 Lake Victoria levels have a possible anthropogenic-climate-change contribution (Earth System Dynamics 2024). — [ESD 15, 225](https://esd.copernicus.org/articles/15/225/2024/esd-15-225-2024.html). The flood wave's timing through the Sudd is followed after the 2019/20 positive IOD — [HESS 2026](https://hess.copernicus.org/articles/30/5297/2026/)
- [Established, figures vary] The Caspian is declining (about 6.1 cm/yr over 1995–2023 in one study, nearly 7 cm/yr in 2020–2025 per a press summary), driven mainly by reduced Volga inflow and evaporation. One GRACE study's "TWSA ≈ −1,200 mm by 2020" looks implausibly large for 300 km-smoothed GRACE and should be checked. — [Duku et al. 2026, Earth's Future](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2025EF008028); [ResearchGate preprint on 2003–2024 trend shifts](https://www.researchgate.net/publication/405543374_Trend_shifts_and_hydroclimatic_drivers_of_Caspian_Sea_Level_2003-2024_revealed_by_GRACE_and_multi-source_observations)

**Ice sheets**
- [Single-study, uses COST-G] Cambier et al. 2025 (Remote Sensing 17:3552) applied EOF and wavelet analysis to COST-G Greenland mass. Mode 1 holds annual and interannual frequencies (4–7 and 11 years). At the annual scale, mass correlates with temperature (r = −0.88), NAO (r = 0.74) and the Greenland Blocking Index (r = −0.85). The authors report r = −0.91 with "p = 10^−274" for combined modes against cumulative indices, and conclude that blocking precedes and intensifies mass loss. — [Remote Sensing 2025](https://dx.doi.org/10.3390/rs17213552); [EGUsphere preprint](https://egusphere.copernicus.org/preprints/2025/egusphere-2025-558/)
- [Single-study] Interannual Antarctic mass variations over 2003–2017 were driven by precipitation and "highly correlated with ENSO" (EPSL 2021). — [ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0012821X21000558)
- [Single-study, 2026] East Antarctica gained about 695 Gt over 2021–2023, the largest gain in the GRACE record, against a long-term loss of about 140.5 Gt/yr. Model experiments attribute it to warming of the western tropical Pacific and eastern Indian Ocean warm pool, which drove a wave train and more atmospheric rivers into Queen Mary and Wilkes Land (Institute of Oceanology, CAS, published in Nature per the press). — [ScienceDaily, Sep 2026](https://www.sciencedaily.com/releases/2026/09/260906170128.htm)
- [Preprint] Atmospheric rivers and winter sea ice drive the recent reversal in Antarctic mass loss, and accumulation accelerated further through December 2024. — [arXiv 2510.03590](https://arxiv.org/pdf/2510.03590)

**Ocean bottom pressure**
- [Established] SAM drives Southern Ocean OBP and Antarctic Circumpolar Current transport. GRACE-based ACC transports correlate more strongly with SAM than any tide-gauge record tested ([Bergmann et al. 2012, JGR Oceans](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2012JC007872)). Southern Indian Ocean transport correlates with SAM at r = 0.61, and at 0.88 with 45–65°S zonal winds ([Makowski et al. 2015](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2014jc010575)). A roughly 1.5-year OBP variation in all basins is phase-locked with SAM, and GRACE interannual OBP exceeds model predictions ([GRACE/ECCO modes](https://www.researchgate.net/publication/253587740_Global_modes_of_variability_and_propagation_in_ocean_bottom_pressure_from_GRACE_and_ECCO)).
- SAM, SOI and IOD each affect southern Indian Ocean OBP through Ekman convergence. — [Frontiers in Marine Science 2022](https://frontiersin.org/articles/10.3389/fmars.2022.916592/full). Interannual South Pacific OBP mechanisms are covered in [Climate Dynamics 2022](https://link.springer.com/article/10.1007/s00382-022-06198-0)
- Near Drake Passage, OBP is contaminated by West Antarctic mass loss and GRACE noise. — [Makowski et al. 2015](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2014jc010575)
- AMOC: synthetic tests show interannual AMOC (20–60°N) is recoverable from about 20 years of GRACE-like OBP with errors under 1 Sv, and real GRACE estimates at 26°N agree with RAPID. Leakage across steep slopes and land hydrology near coasts are the main problems. — [Landerer et al. 2015, Ocean Science 11:953](https://os.copernicus.org/articles/11/953/2015/); [US CLIVAR abstract](https://usclivar.org/sites/default/files/meetings/2015/abstracts/Landerer_Felix_2_AMOC2017_86317.pdf)

**MJO, QBO and the solar cycle**
- [Gap] I found no GRACE study linking TWS to the MJO. MJO–extreme-rainfall work, for example over Indonesia, uses gauges and APHRODITE. — [arXiv 2007.10574](https://arxiv.org/pdf/2007.10574)
- [Contested/weak, single study] Scientific Reports 2025 argues that land water volume change has an 11-year solar-cycle imprint through a QBO → ENSO → rainfall chain. The reported correlation with sunspot number is r = 0.12 with 115 degrees of freedom, called "significant at over 80%", which is below any conventional threshold. — [Sci. Rep. 2025](https://www.nature.com/articles/s41598-025-99880-2)

### Inferences
- The ENSO–TWS link is robust in sign and location, but almost every published significance level assumes far more independent samples than exist. Phillips' implied ~90 degrees of freedom contrasts with an effective sample of order 10–30 for monthly TWS against ONI once autocorrelation is included (see Section 5). Many of the weaker mid- and high-latitude "teleconnections" in the literature, including ENSO–Greenland, are therefore probably not significant.
- Correlating cumulative indices with a trending ice-mass series, as in the Cambier et al. p = 10^−274, is the textbook setting for spurious regression between integrated series. Those p-values should not be taken at face value. The physical blocking → melt link is well supported by other work.
- Rodell & Li 2023 (intensity tracks global temperature more than ENSO; Section 4) and Rateb & Scanlon 2025 (ENSO dominates extremes) use different metrics: summed intensity vs extremes occurrence/synchrony. They are complementary rather than strictly contradictory, but the tension is worth testing directly.
- ENSO diversity (central- vs eastern-Pacific events) has regional GRACE evidence (Nicaragua, Yangtze) but no global analysis. This is a cheap, open target.
- OBP–climate links are established mainly for SAM. OBP against ENSO, IOD and AO beyond the Southern Ocean is thin; a 2026 Nature Index overview mentions EOF modes tied to SAM, ENSO and PDO, but I did not check the primary sources.

### Gaps
- No GRACE study of MJO aggregate effects. A monthly-mean RMM amplitude by phase could be regressed on Maritime Continent and East African TWS, but it would be confounded with ENSO.
- No rigorous GRACE–QBO or GRACE–solar test. With about two solar cycles in the record (2002–2026), the solar hypothesis cannot reach meaningful power.
- No global CP-vs-EP ENSO TWS fingerprint study.
- Arctic Ocean OBP against AO with GRACE was not found in these searches.
- I could not read the full text of Rateb & Scanlon 2025 (HTTP 403) or of the Climate Dynamics fingerprint paper (paywall), so their exact significance methods are unknown.

---

## 2. Predictability: does TWS lead drought, heatwaves, monsoon failure, wildfire, crop yield, floods, lake levels or dust? Does GRACE add skill?

### Takeaway
Antecedent TWS has documented lead skill for floods (Missouri 2011: 5–11 months; much shorter in a global analysis), Amazon fire-season severity (1–5 months) and hydrological drought (1–6 months in data-assimilation systems), and for the interannual CO2 growth rate (a contemporaneous sensitivity). The added value over soil moisture and precipitation is real but "modest" where it has been quantified, and comes from the deep (groundwater and root-zone) memory. Heatwaves, crop yield, monsoon failure and dust have almost no GRACE predictive studies. These are the clearest openings.

### Cited Findings
- [Single-study, often cited] Reager et al. 2014 (Nature Geoscience): GRACE basin storage indicated flood predisposition 5–11 months ahead for the 2011 Missouri floods. An unusually high discharge (about a 1-in-20-month event) was flagged 11 months ahead, and a 1-in-100-month event 5 months ahead. Model snow and soil moisture were useful only about 2 months ahead. — [Nature Geoscience](http://www.nature.com/ngeo/journal/v7/n8/full/ngeo2203.html); [NASA JPL](https://www.jpl.nasa.gov/news/nasa-satellite-data-give-early-clues-to-flood-danger/)
- [Contested in magnitude] A preliminary global analysis of daily and monthly GRACE found lead times of about 2 months for annual flooding, much shorter than the Missouri headline. — [AGU FM 2020 abstract](https://ui.adsabs.harvard.edu/abs/2020AGUFMH064.0008R/abstract)
- [Single-study] Five-day GRACE/GRACE-FO solutions compared with about 3,300 Dartmouth Flood Observatory events (2002–2021) show saturation-excess floods have detectable antecedent storage increases 15–50 days before and during the flood. Monthly products have 40–60-day latency. — [Sci. Rep. 2024](https://www.nature.com/articles/s41598-024-64491-w)
- [Single-study] Chen et al. 2013 (JGR Biogeosciences): GRACE TWS deficits preceded severe fire seasons in the southern Amazon, with the strongest links from April to August, i.e. 1–5 months before the September peak (2002–2011). — [JGR-B](https://agupubs.onlinelibrary.wiley.com/doi/10.1002/jgrg.20046). The UC Irvine fire-season forecasts used GRACE TWS as a supporting signal next to SST: [2014](https://www.ess.uci.edu/%7Eamazonfirerisk/ForecastWeb/SAMFSS2014.html), [2015](https://www.ess.uci.edu/~amazonfirerisk/ForecastWeb/SAMFSS2015.html)
- [Established, nature of skill] GRACE assimilation into the Catchment model gave "modest, but statistically significant" improvements in hydrological skill and drought indicators for the US Drought Monitor (Houborg et al. 2012). The added value comes from deep soil moisture and groundwater, which the monitors otherwise lacked. — [NASA NTRS](https://ntrs.nasa.gov/citations/20100026436); [NASA JPL drought page](https://grace.jpl.nasa.gov/applications/drought-monitoring/); seasonal groundwater forecast initialisation over the US — [J. Hydrometeorology 2020](https://journals.ametsoc.org/view/journals/hydr/21/1/jhm-d-19-0096.1.xml)
- Operational products: weekly GRACE-DA groundwater and root-zone indicators exist for the US (V4.0, 0.125°) and globally (V3.0, 0.25°). They are meant to show climate-driven wetness, not pumping-driven aquifer depletion. — [data.nasa.gov global V3.0](https://data.nasa.gov/dataset/groundwater-and-soil-moisture-conditions-from-grace-and-grace-fo-data-assimilation-l4-7-da-31ce9); [NASA GSFC](https://earth.gsfc.nasa.gov/hydro/data/groundwater-and-soil-moisture-conditions-grace-data-assimilation)
- A NASA drought-prediction project found skill globally at 1-month lead and up to 3 months regionally, with skill depending on active teleconnections (ENSO, IOD). — [NASA NTRS](https://ntrs.nasa.gov/search.jsp?R=20160002258)
- Sub-seasonal to seasonal TWS forecasting (HESS 2026, FEWS NET LDAS): ROC scores above 0.6 for tercile forecasts over more than 50% of the domain at 1–6-month leads, with skill very sensitive to precipitation variability. — [HESS 30, 1097](https://hess.copernicus.org/articles/30/1097/2026/)
- Li et al. 2024 (GRL) forecast global land water storage up to 1 year ahead with "acceptable skill", performing well up to 6 months over most land. — [GRL 2024](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2024GL109101?af=R)
- [Established, single high-profile paper with follow-ups] Humphrey et al. 2018 (Nature 560:628): the interannual atmospheric CO2 growth rate is strongly sensitive to GRACE TWS, with drier years showing faster growth, independent of tropical temperature effects. Carbon-cycle models underestimate the link. — [IDEAS/Nature](https://ideas.repec.org/a/nat/nature/v560y2018i7720d10.1038_s41586-018-0424-4.html); a Reading follow-up (Petch et al.) covers tropical forests — [Reading blog 2024](https://blogs.reading.ac.uk/weather-and-climate-at-reading/2024/land-surface-water-controls-on-atmospheric-co2-growth)
- Heatwaves: Chen et al. 2019 (JGR Atmospheres) used GRACE TWS for global soil moisture–air temperature coupling and found it controlled mainly by root-zone moisture. GRACE helped most where vegetation is deep-rooted (eastern South America, southern Africa). They note that prior hot-extreme coupling studies used precipitation indices or modelled soil moisture rather than GRACE. — [JGR-A 2019](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2019JD030324). A review reports TWS falling 1–3 months after rises in temperature and PET — [Geocarto 2025 review](https://www.tandfonline.com/doi/full/10.1080/10106049.2025.2482706)
- Central European drought: GRACE-FO deficits of 112 ± 10.5 Gt (2018) and 145 ± 12 Gt (2019), equal to 73% and 94% of the mean seasonal amplitude. 2019 was the largest deficit in the record, and no GRACE/GRACE-FO offset was found. — [Boergens et al. 2020, GFZ PDF](https://gfzpublic.gfz-potsdam.de/pubman/item/item_5002572_4/component/file_5002739/5002572.pdf)
- Crop yield (US maize, counties): higher pre-season GRACE storage and soil moisture worsened yield loss from excessive rain but eased loss from extreme drought. Only 53% of county-years with excessive rain had negative yield anomalies. — [PMC6850578](https://pmc.ncbi.nlm.nih.gov/articles/PMC6850578/)
- Monsoon: a single-paper finding is that La Niña monsoon rains filled SE Asian soils to capacity while El Niño caused modest decreases, giving a nonlinear precipitation–storage relation. — [Rateb & Scanlon 2025 (ResearchGate)](https://www.researchgate.net/publication/397797513_Dynamics_and_Couplings_of_Terrestrial_Water_Storage_Extremes_From_GRACE_and_GRACE-FO_Missions_During_2002-2024). Groundwater "hydraulic memory" is long (accumulated annual precipitation) in drylands and short (monthly precipitation) in humid aquifers — [ResearchGate](https://www.researchgate.net/publication/343918275_Climate-groundwater_dynamics_inferred_from_GRACE_and_the_role_of_hydraulic_memory)
- Dust: GRACE documents large Middle East groundwater loss (Tigris–Euphrates and western Iran: 17.3 ± 2.1 mm/yr, about 91.3 ± 10.9 km³, 2003–2009) — [Voss et al. 2013, WRR](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1002/wrcr.20078). I found no study linking GRACE storage to dust emission.
- Current context: CPC's ONI is +1.39 (MJJ 2026), +1.80 (JJA) and +2.16 (JAS 2026), so a strong El Niño is under way. — [CPC oni.ascii.txt](https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt). Press reports intense Amazon burning before the peak — [Gizmodo](https://gizmodo.com/the-super-el-nino-hasnt-even-peaked-and-its-already-burning-up-the-amazon-2000826295) (news source, low weight)

### Inferences
- The robust predictive mechanism is memory. TWS integrates P − ET − R, so it carries information that monthly precipitation and surface soil moisture do not, at lead times set by the storage e-folding time. Predictability should therefore be highest where groundwater or deep root-zone storage dominates and where a large teleconnection (ENSO, IOD) is active.
- "Does GRACE add skill?" is mostly answered for drought monitoring and initialisation (modest but significant) and barely tested as a predictor of downstream impacts (heat, fire outside the Amazon, yields, dust). Those tests are cheap with monthly TWS plus public impact datasets. The right design is incremental skill (cross-validated R² or ROC gain) over a baseline of antecedent precipitation (GPCP) plus surface soil moisture (ESA CCI).
- The ongoing 2026 El Niño (ONI +2.16 in JAS) offers a rare chance for **pre-registered** out-of-sample prediction. Our record ends in 2026-07, so predictions for Aug 2026 to 2027 TWS, fire and CO2 growth can be written down now and checked when the COST-G fields arrive.

### Gaps
- No GRACE-based test of heatwave predictability beyond soil moisture was found.
- No global (non-Amazon) test of pre-season TWS against burned area was found. Canada 2023, Australia 2019–20 and Siberia were not checked.
- No GRACE predictor of monsoon onset or failure was found, i.e. pre-monsoon TWS against subsequent all-India or West African rainfall.
- No GRACE–dust link was found.
- No quantitative "with vs without GRACE" drought-skill comparison beyond Houborg's "modest" result was found.

---

## 3. Atmospheric rivers (ARs) and extreme precipitation in GRACE

### Takeaway
AR-driven mass gains are clearly detected in monthly GRACE/GRACE-FO at the basin and season scale: California water year 2023 (cross-checked with GNSS) and East Antarctic snowfall 2009, 2011 and 2021–2023. Event-scale detection needs daily or 5-day solutions, which see storage build-up 15–50 days before saturation-excess floods. I found no global study relating monthly AR frequency or intensity (for example from the Guan–Waliser catalogue) to TWS change. That is an open niche.

### Cited Findings
- [Established] California water year 2023: GNSS loading shows TWS gains in key watersheds over the first 6 months that were more than twice average and among the highest in about two decades, relieving meteorological and then (with a lag) hydrological drought. — [CW3E summary of Martens et al.](https://cw3e.ucsd.edu/gnss-geodesy-quantifies-water%E2%80%90storage-gains-and-drought-improvements-in-california-spurred-by-atmospheric-rivers/)
- GRACE/GRACE-FO over the Sacramento, San Joaquin and Tulare basins showed the largest winter net gain in about 22 years (October 2022 to March 2023), though groundwater still reflected drought. — [NASA Earth Observatory](https://earthobservatory.nasa.gov/images/151438/satellites-show-gains-in-california-water)
- J. Hydrology 2025 attributes California's WY2023 storage gains to high-category AR landfalls and record snowfall, using GNSS, GRACE and NLDAS. — [ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0022169425001428); [AGU24 abstract](https://studio.m-anage.com/agu/agu24/meetingapp.cgi/Paper/1622328)
- Antarctica: nine ARs produced 80% of the exceptional 2009/2011 snow accumulation at Princess Elisabeth station. — [ScienceDaily 2015](https://www.sciencedaily.com/releases/2015/01/150120112206.htm). ARs account for up to 70% of extreme snowfall events in East Antarctica since the 1980s, and AR-fed snowfall drove the 2021–2023 GRACE mass gain of about 695 Gt — [ScienceDaily 2026](https://www.sciencedaily.com/releases/2026/09/260906170128.htm); [arXiv 2510.03590](https://arxiv.org/pdf/2510.03590)
- Sub-monthly GRACE: 5-day solutions detect flood precursors 15–50 days ahead ([Sci. Rep. 2024](https://www.nature.com/articles/s41598-024-64491-w)). Reconstructed daily TWS captures storage changes from extreme rainfall in the Haihe basin ([Remote Sensing 2023](https://dx.doi.org/10.3390/rs15061564)). Daily downscaling over CONUS was tested on Hurricane Harvey ([J. Hydrology 2024](https://www.sciencedirect.com/science/article/abs/pii/S0022169424015907)). High-frequency GRACE variability can be reconstructed from filtered precipitation anomalies ([Humphrey et al. 2016, Surveys in Geophysics](https://link.springer.com/article/10.1007/s10712-016-9367-1)).
- AR catalogues: Guan–Waliser ERA5 global hourly ARs for 1979–2021 are on PANGAEA, doi:10.1594/PANGAEA.957161 ([PANGAEA](https://doi.pangaea.de/10.1594/PANGAEA.957161)). The tARget algorithm is at version 4 (Guan & Waliser 2024, Scientific Data), with the database on [UCLA Dataverse](https://dataverse.ucla.edu/dataverse/ar). ARTMIP Tier 2 holds multi-algorithm ERA5 AR tags for 2000–2019 ([NCAR GDEX d651018](https://gdex.k8s.ucar.edu/datasets/d651018/)).

### Inferences
- Monthly COST-G at degree 90 with 300 km smoothing can resolve seasonal AR-integrated gains in large basins (above about 200,000 km²: Central Valley plus the Sierra, Pacific Northwest, central Chile, Norway, East Antarctic drainage basins). Single AR events (1–3 days) cannot be resolved; only their monthly aggregate can.
- A physically motivated, untested question is whether AR precipitation is more "storage-efficient" than non-AR precipitation. Snowpack and saturation could raise dS/P, or runoff losses could lower it. This is testable with monthly dS/dt, GPCP and AR-frequency counts (Hypothesis H3).

### Gaps
- No global or multi-region study regresses monthly GRACE dS/dt on AR frequency or IVT outside California and Antarctica.
- The PANGAEA ERA5 catalogue ends in 2021. Coverage for 2022–2026, including California 2023, depends on the v4 Dataverse release, whose ERA5 end date I did not verify.

---

## 4. Global hydrological-cycle intensification, global land-storage shifts (2014–2016 drop, 2023–2024 records) and regime shifts

### Takeaway
There is strong agreement that global non-ice land TWS dropped abruptly between May 2014 and March 2016, coincident with the 2014–16 El Niño and NE South American drought, and has not recovered. 2023 and 2024 set successive record lows, and the 2023–24 El Niño produced a barystatic sea-level jump as large as 2015–16 but faster. Whether this is a "permanent regime shift" is a hypothesis, not established. Claims of a soil-moisture regime change (Seo et al. 2025) and of a 10% rise in evapotranspiration (Pascolini-Campbell 2021, apparently retracted) are contested. Extreme-event intensity tracks global temperature in one study and ENSO in another.

### Cited Findings
- [Established] Rodell et al. 2024 (Surveys in Geophysics, doi:10.1007/s10712-024-09860-w): global land TWS excluding ice sheets and glaciers declined rapidly from May 2014 to March 2016. By 2023 its upper range was still about 1 cm EWH below the earlier range. The decline began with a record NE South America drought, and later droughts on five continents prevented recovery. A divergence between two barystatic sea-level estimates after about 2015 is not clearly linked to it. — [PMC full text](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11671563/); [NASA GRACE-FO news](https://gracefo.jpl.nasa.gov/news/287/nasa-satellites-reveal-abrupt-drop-in-global-freshwater-levels/)
- An EGU25 abstract (ETH Zürich) poses the regime-shift hypothesis for the 2015–16 decline coinciding with El Niño. — [EGU25-17534](https://meetingorganizer.copernicus.org/EGU25/EGU25-17534.html)
- [Established] Global TWS hit record lows of −6,697 km³ in 2023 (−759 km³ vs 2022) and −7,404 km³ in 2024 (−796 km³ vs 2023). In 2024, gains in Africa, Australia, Europe and parts of Antarctica were outweighed by losses in NW Canada, South America, southern Africa and Greenland. — [Nature Reviews Earth & Environment 2023](https://www.nature.com/articles/s43017-024-00545-x); [2024](https://www.nature.com/articles/s43017-025-00659-w). A NASA version gives the 2023 anomaly as −9.94 cm EWH — [NTRS PDF](https://ntrs.nasa.gov/api/citations/20240004868/downloads/Li_waterstorage_final.pdf)
- [Single-study] Jin et al. 2026 (JGR Oceans): barystatic sea-level rise was about 6.2 mm in both the 2014–16 and 2023–24 El Niños, but much faster in 2023–24, driven mainly by land-water loss in the Americas and Africa, with South America fastest. — [JGR Oceans](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2025JC023533)
- [Established] La Niña does the reverse. In 2010–11, anomalous low-latitude land storage, dominated by Australia's internally draining basins, dropped global mean sea level by about 7 mm over about 18 months. IOD and SAM modulated the La Niña rainfall. — [Fasullo et al. 2013 PDF](https://climate.fas.harvard.edu/sites/projects.iq.harvard.edu/files/climate/files/fasullo.etal_2013.pdf); [AGU release](https://news.agu.org/press-release/global-sea-level-rise-dampened-by-australia-floods); Boening et al. 2012, GRL ("The 2011 La Niña: so strong, the oceans fell"), as summarised by [Skeptical Science](https://skepticalscience.com/news.php?n=2164)
- [Single-study, high-profile] Chandanpurkar et al. 2025 (Science Advances 11:eadx0298), April 2002 to April 2024: drying land area grows by about twice California's area each year. "Mega-drying" regions have formed by the linking of earlier hotspots since the 2014(–16) El Niño, mostly in the Northern Hemisphere. Groundwater depletion is 68% of TWS loss over non-glaciated land, and continents now contribute more to sea-level rise than the ice sheets. — [Science Advances](https://www.science.org/doi/10.1126/sciadv.adx0298); [ASU news](https://news.asu.edu/20250725-environment-and-sustainability-new-global-study-shows-freshwater-disappearing-alarming)
- [Contested] Seo et al. 2025 (Science 387:1408, doi:10.1126/science.adq6529) use ERA5 soil moisture: a loss of about 1,614 Gt over 2000–2002 and a further 1,009 Gt over 2003–2016, consistent with about 4.4 mm of GMSL rise and about 45 cm of pole drift, described as a "permanent" regime change. A non-peer-reviewed critique says supplementary datasets other than ERA5 do not match. — [PolyU record](https://research.polyu.edu.hk/en/publications/abrupt-sea-level-rise-and-earths-gradual-pole-shift-reveal-perman/); [critique (Substack)](https://bioticregulation.substack.com/p/we-are-losing-soil-moisture-why)
- [Established-ish] Rodell & Li 2023 (Nature Water 1:241, doi:10.1038/s44221-023-00040-5) found 1,056 TWS extreme events over 2002–2021. Total intensity tracks global mean temperature more closely than ENSO. The wet-to-dry ratio fell from southern Europe to SW China and rose in the Sahel (5–20°N). The most intense event was the central African pluvial from 2019. — [NASA NTRS](https://ntrs.nasa.gov/citations/20230003509); [Climate Signals](https://climatesignals.org/scientific-reports/changing-intensity-hydroclimatic-extreme-events-revealed-grace-and-grace-fo)
- [Contested/retracted] Pascolini-Campbell et al. 2021 (Nature 593:543, doi:10.1038/s41586-021-03503-5) reported a 10 ± 2% rise in land ET over 2003–2019 from a GRACE water balance. An IDEAS/RePEc listing marks it as **"RETRACTED ARTICLE"**. I could not retrieve the reason before the search budget ran out. — [IDEAS listing](https://ideas.repec.org/a/nat/nature/v593y2021i7860d10.1038_s41586-021-03503-5.html)
- [Single-study] Rateb & Scanlon 2025 report a shift in global extremes from wet-dominated to dry-dominated around 2011–2012. — [EurekAlert](https://www.eurekalert.org/news-releases/1112360)

### Inferences
- Two independent "step" analyses (Rodell 2014–16; Rateb 2011–12) put the break at different dates. With 24 years of strongly autocorrelated data and an 11-month gap in 2017–18 right after the 2014–16 drop, changepoint dates are poorly constrained. Neither study (as far as I could see) tested the step against a red-noise or long-memory null. That is a concrete, testable gap.
- Step-like declines during strong El Niños, with little La Niña recovery, would amount to an ENSO-asymmetric "ratchet" layered on groundwater depletion. I found no study that tests ENSO-phase asymmetry in global TWS response directly. The 2026 El Niño gives a third event to test it.
- GIA matters for global and regional trends (Canada, Fennoscandia, Antarctica, Patagonia, Alaska), but a linear GIA term cannot create a step change. Step and changepoint analyses on GIA-uncorrected data are valid if the null model includes a linear trend.

### Gaps
- No formal significance test of the 2014–16 step against red-noise or long-memory surrogates was found.
- The reason for the Pascolini-Campbell 2021 retraction was not retrieved.
- Whether the Seo et al. soil-moisture regime change is consistent with GRACE TWS (rather than only ERA5) was not resolved.
- Global hydrological intensification is shown in GRACE by Rodell & Li 2023, but whether TWS variance or the size of the seasonal cycle has changed over 2002–2026 was not found as a tested result. A 2020 Remote Sensing paper proposes emerging variability changes as a target for future missions ([doi:10.3390/rs12233898](https://doi.org/10.3390/rs12233898)), but I did not read it.

---

## 5. Methods for robust detection: field significance, effective degrees of freedom, detrending, surrogates, causal discovery; pitfalls

### Takeaway
GRACE–climate studies mostly use pointwise correlation, EOF/ICA and wavelet coherence. Field-significance control (FDR) and autocorrelation-aware degrees of freedom are rarely reported. Causal discovery (PCMCI, CCM, transfer entropy) has been benchmarked on hydrometeorological data but applied to GRACE only regionally or for predictor ranking. The main risks of spurious discovery are: trends (groundwater depletion, ice loss, uncorrected GIA); integrated or cumulative indices; few independent ENSO cycles; tidal or 161-day aliasing; the 2017–18 gap; and the multiplicity of grid cells × indices × lags.

### Cited Findings
- Wilks 2016 (BAMS 97:2263, doi:10.1175/BAMS-D-15-00267.1), "The stippling shows statistically significant grid points": testing each grid point without a multiplicity adjustment routinely overstates results. FDR control is recommended because it is only modestly sensitive to spatial autocorrelation. Earlier work: Wilks 2006 (JAMC) on field significance and FDR, and Ventura et al. 2004 (J. Climate), where Benjamini–Hochberg still controlled false rejections under realistic spatial correlation. — [AMS 2016 abstract](https://ams.confex.com/ams/96Annual/webprogram/Paper276937.html); [Ventura et al. 2004 PDF](https://www.stat.berkeley.edu/~paciorek/files/jclim/vent.etal.2004.pdf)
- [Contested detail] A 2026 preprint argues that global FDR can mislead for spatially coherent signals and that domain restriction changes the threshold. A reviewer notes that α_FDR = 0.1 is a choice, not a property of FDR. — [EGUsphere 2026-2203 AC1](https://egusphere.copernicus.org/preprints/2026/egusphere-2026-2203/egusphere-2026-2203-AC1-supplement.pdf); [AC2](https://egusphere.copernicus.org/preprints/2026/egusphere-2026-2203/egusphere-2026-2203-AC2-supplement.pdf)
- An example of the degrees-of-freedom problem in the GRACE literature: Phillips 2012 used r_crit = 0.21 and acknowledged "unclear degree of freedom" (Section 1). — [Phillips 2012 PDF](https://escholarship.org/content/qt9p53g2p8/qt9p53g2p8.pdf?t=n0u7rd)
- Causal-discovery benchmarks: Ombadi et al. 2020 (WRR) compared Granger causality, transfer entropy, the PC algorithm and CCM against a bucket model with known causality. — [WRR](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020WR027251). A 2025 comparison of CCF, CCM, TE and PCMCI+ on hydrological systems — [SERRA 2025](https://link.springer.com/article/10.1007/s00477-025-02977-3). PCMCI+ was most stable across climates, and TWS emerged as a causal driver — [EGUsphere 2025-4650](https://egusphere.copernicus.org/preprints/2025/egusphere-2025-4650/egusphere-2025-4650.pdf). PCMCI cannot orient contemporaneous links, which matters for monthly data where the response lag is under a month — same comparison sources.
- GRACE + causality so far: an ESSD 2026 pre-GRACE TWS reconstruction ranks predictors with Granger causality — [ESSD 18, 1747](https://essd.copernicus.org/articles/18/1747/2026/). Pixelwise Granger tests between vegetation and TWS — [Remote Sensing 2025](https://doi.org/10.3390/rs17101701). A Heilongjiang groundwater study using PCMCI and TE was mentioned in search results, but I did not obtain its URL or results. I found no global PCMCI or CCM study of TWS against multiple climate indices.
- Wavelet coherence is the dominant method in recent GRACE–teleconnection papers (Rateb 2025, Sudd, Yunnan, Greenland/COST-G), and Rateb explicitly notes it cannot establish causality. — [AGU Advances](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2025AV001684)
- GRACE-specific artefacts: the S2 alias (161 days) is still visible in RL06 over West Antarctica. Mascon regularisation may amplify aliasing on land. Part of the 161-day signal may come from accelerometer errors linked to the 161-day β′ cycle. Tide errors can contaminate Antarctic trends. — [Remote Sensing 2022, doi:10.3390/rs14215403](https://doi.org/10.3390/rs14215403); [Seo et al. 2008 (alias error)](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2006JB004747); [S2 tide aliasing in GRACE solutions, J. Geodesy 2008 (CSR)](https://filedrop.csr.utexas.edu/pub/ggfc/papers/jog_2008.pdf)
- The GRACE to GRACE-FO gap is 11 months (2017–2018), and Rateb filled it with probabilistic reconstruction. — [EurekAlert](https://www.eurekalert.org/news-releases/1112360). No GRACE/GRACE-FO offset was found globally — [Boergens 2020](https://gfzpublic.gfz-potsdam.de/pubman/item/item_5002572_4/component/file_5002739/5002572.pdf)
- Leakage: coastal land–ocean leakage and steep-slope leakage limit OBP and AMOC estimates, and land hydrology contaminates coastal OBP. — [Landerer et al. 2015](https://os.copernicus.org/articles/11/953/2015/)

### Inferences
(The standard references in this subsection come from background knowledge and were not re-fetched this session: Bretherton et al. 1999, J. Climate 12:1990, doi:10.1175/1520-0442(1999)012<1990:TENOSD>2.0.CO;2; Livezey & Chen 1983, MWR 111:46; Ebisuzaki 1997, J. Climate 10:2147; Runge et al. 2019, Science Advances, doi:10.1126/sciadv.aau4996; Sugihara et al. 2012, Science, doi:10.1126/science.1227079.)
- **Effective sample size.** For two AR(1) series with lag-1 autocorrelations r₁ and r₂, N_eff ≈ N(1 − r₁r₂)/(1 + r₁r₂) (Bretherton 1999). For deseasonalised, undetrended TWS (r₁ ≈ 0.9) against ONI (r₂ ≈ 0.95) with N ≈ 250, N_eff ≈ 250 × 0.145/1.855 ≈ 20. The 95% critical |r| is then about 0.44, not about 0.12. These r values are typical and illustrative; compute them for real per grid cell.
- **Spatial multiplicity.** A 2° grid has about 4,000–5,000 land cells, but 300 km Gaussian smoothing plus degree-90 truncation gives only of order 500–1,000 independent footprints (land area ≈ 1.5×10⁸ km² ÷ π(300 km)² ≈ 530). Testing about 10 indices × 13 lags (0–12 months) per cell means about 10⁵ tests. Pointwise α = 0.05 maps will be stippled by chance. Use FDR (α_FDR = 0.1, as Wilks suggests) plus a global surrogate test: phase-randomise the index (Ebisuzaki) or apply IAAFT, keep the TWS field intact so spatial correlation is preserved, and record the maximum statistic over lags and cells for each surrogate.
- **Trends and GIA.** Uncorrected GIA is a near-linear trend. Detrending or first-differencing (dS/dt) removes it exactly for interannual work, but **any trend, step-size or long-term result in Canada, Fennoscandia, Greenland, Antarctica, Alaska or Patagonia needs a GIA model** (for example ICE-6G_D or Caron et al. 2018). Groundwater depletion and ice loss are also trends: correlating two trending series gives spurious r ≈ ±1. Cumulative indices, as in the Cambier p = 10^−274, are integrated (random-walk-like) series, which makes the spurious-regression risk extreme.
- **Seasonality.** Remove a fitted annual, semiannual and 161-day term rather than monthly climatology alone, so that the S2 alias does not leak into "teleconnection" bands. Teleconnections are season-locked (e.g. IOD in SON), so season-stratified tests leave only about 24 samples per season.
- **Too few events.** 2002–2026 contains fewer than about 8 El Niño events and 2–3 strong ones (2015–16, 2023–24, 2026–27); count them from the ONI file. Any ENSO-asymmetry or CP-vs-EP result is power-limited. Report effect sizes with bootstrap confidence intervals, not just p-values.
- **A built-in out-of-sample split.** GRACE (2002–2017) and GRACE-FO (2018–2026) are independent missions. Fit on one and test on the other. Also pre-register forecasts for the 2026–27 El Niño before the fields are released.
- **Causal discovery.** With monthly data, most physical responses are under a month, so links are contemporaneous and PCMCI+ cannot orient them. CCM assumes deterministic nonlinear dynamics and is fragile with strong seasonality and short records. Use PCMCI+ on deseasonalised, detrended or differenced series with a small, physically motivated variable set: indices, regional TWS, GPCP precipitation. Treat the results as hypothesis generation, then confirm with an out-of-sample split.

### Gaps
- I did not find any GRACE paper that reports Livezey–Chen or FDR field significance for TWS–climate index maps. That does not mean none exist; the search budget ran out.
- I did not find a GRACE study using IAAFT or phase-randomised surrogates for teleconnection significance.

---

## 6. Public data needed (sources, URLs, formats)

### Takeaway
All the climate indices needed are free, plain-text and monthly from NOAA PSL and CPC, with no account. The Guan–Waliser AR catalogue is free from PANGAEA or UCLA Dataverse. ERA5 needs a CDS account. NASA GES DISC products, such as GRACE-DA drought indicators, normally need a free Earthdata login (not verified here).

### Cited Findings
- NOAA PSL climate-index list: https://psl.noaa.gov/data/climateindices/list/. Paths (prefix https://psl.noaa.gov): Niño 3.4 anomaly `/data/correlation/nina34.anom.data` (CPC, ERSSTv5); ONI `/data/correlation/oni.data`; MEI.v2 `/enso/mei/data/meiv2.data` (bimonthly; "Jan" = Dec–Jan); PDO `/data/correlation/pdo.data`; AMO unsmoothed `/data/correlation/amon.us.data` (long: `amon.us.long.data`; Kaplan SST); NAO `/data/correlation/nao.data` (CPC) or `jonesnao.data` (CRU); AO `/data/correlation/ao.data`; SAM/AAO `/data/correlation/aao.data`; QBO `/data/correlation/qbo.data` (30 hPa equatorial zonal wind, NCEP/NCAR R1); solar 10.7 cm flux `/data/correlation/solar.data`; TNA `/data/correlation/tna.data`; WHWP `/data/correlation/whwp.data`. Format: a header line with start and end years, then one row per year with the year and 12 monthly values, then a missing-value line at the end. DMI, sunspots, Greenland blocking and NPGO are not on this list. — [PSL list](https://psl.noaa.gov/data/climateindices/list/)
- CPC ONI: https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt. Columns are SEAS (3-month season), YR, TOTAL (°C), ANOM. The latest rows are MJJ 2026 +1.39, JJA +1.80, JAS +2.16. — [CPC](https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt)
- DMI: the PSL GCOS-WGSP version is recomputed from HadISST1.1, with 1981–2010 climatology, monthly, 1870 to present — [PSL DMI page](https://psl.noaa.gov/gcos_wgsp/Timeseries/DMI). The east pole file is `psl.noaa.gov/gcos_wgsp/Timeseries/Data/dmieast.had.long.data` — [IRI Data Library](https://iridl.ldeo.columbia.edu/SOURCES/.NOAA/.PSL/.DMI/.DMI-EAST). The CPC DMI is ERSSTv5-based with 1991–2020 climatology, and values can be revised for up to 2 months — [CPC DMI monthly](https://www.cpc.ncep.noaa.gov/products/international/ocean_monitoring/IODMI/DMI_month.html)
- AR catalogues: [PANGAEA ERA5 Guan–Waliser, 1979–2021, hourly](https://doi.pangaea.de/10.1594/PANGAEA.957161); [UCLA Dataverse tARget v4 (2024)](https://dataverse.ucla.edu/dataverse/ar); [ARTMIP Tier 2 (NCAR GDEX d651018)](https://gdex.k8s.ucar.edu/datasets/d651018/)
- GRACE-DA drought indicators (GES DISC, global V3.0 at 0.25°, weekly): [data.nasa.gov](https://data.nasa.gov/dataset/groundwater-and-soil-moisture-conditions-from-grace-and-grace-fo-data-assimilation-l4-7-da-31ce9)

### Inferences
- The PSL PDO description says "Using data from 1948 to 2002" (the EOF basis period). Use the same index version throughout and record the download date, because CPC and PSL revise values: the CPC DMI by up to 2 months, and ONI uses sliding 30-year base periods.
- The ONI is a 3-month running mean. Using it at monthly resolution adds autocorrelation, so lower N_eff accordingly, or use the monthly Niño 3.4 anomaly instead.

### Gaps (sources known from background knowledge; URLs **not verified** this session)
- GPCP monthly precipitation v2.3 (2.5°): NOAA PSL gridded page `https://psl.noaa.gov/data/gridded/data.gpcp.html`, or NCEI CDR.
- SPEI drought index: SPEIbase at `https://spei.csic.es/database.html`.
- Burned area: GFED5 (`https://www.globalfiredata.org`) or MODIS MCD64A1 (Earthdata login).
- CO2 growth rate: NOAA GML Mauna Loa and global monthly, `https://gml.noaa.gov/ccgg/trends/`.
- Soil moisture: ESA CCI SM (`https://climate.esa.int/en/projects/soil-moisture/`).
- Daily Tmax: CPC Global Unified Temperature on PSL, or Berkeley Earth.
- Sunspots: SILSO (`https://www.sidc.be/SILSO/datafiles`).
- QBO alternative: FU Berlin station-based winds.
- SAM station index: Marshall (`https://legacy.bas.ac.uk/met/gjma/sam.html`).
- Greenland Blocking Index: NOAA PSL `/data/timeseries/daily/GBI/` (path unverified).
- RMM MJO index: BoM.
- Crop yields: GDHY (Iizumi) or FAOSTAT.
- GIA models: ICE-6G_D, Caron 2018 (JPL).
- ENSO E/C indices (Takahashi et al. 2011): IGP Peru.
- Verify all of these before scripting downloads.

---

## 7. Ranked testable hypotheses (monthly TWS + public indices)

### Takeaway
Ranking weighs novelty (not found in the literature searched), feasibility with our COST-G 2° grids, and robustness to the pitfalls in Section 5. The top three can be done with data already in hand plus free text files. H1 includes a live, pre-registered test of the 2026–27 El Niño, which is the strongest protection against fooling ourselves. **Common protocol for all:** deseasonalise by fitting annual, semiannual and 161-day terms; analyse first differences dS/dt (which removes linear GIA) or detrended anomalies; use gap-aware statistics with no interpolation across 2017-07 to 2018-05; fit on GRACE (2002–2017) and test on GRACE-FO (2018–2026); compute significance with index surrogates (phase-randomised or IAAFT, at least 1,000), taking the maximum over lags and cells; use FDR at q = 0.1 for maps; report effect sizes with block-bootstrap confidence intervals.

### Cited Findings
(The evidence behind each hypothesis is cited in Sections 1–5; key anchors are repeated here.)
- Steps in 2014–16 and 2023–24, and the 2010–11 La Niña rebound: [Rodell 2024](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11671563/), [Jin 2026](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2025JC023533), [Fasullo 2013](https://climate.fas.harvard.edu/sites/projects.iq.harvard.edu/files/climate/files/fasullo.etal_2013.pdf). ONI is +2.16 in JAS 2026 — [CPC](https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt)
- AR-driven gains: [California 2023](https://www.sciencedirect.com/science/article/abs/pii/S0022169425001428); [Antarctica 2021–23](https://www.sciencedaily.com/releases/2026/09/260906170128.htm)
- Fire: [Chen 2013 Amazon](https://agupubs.onlinelibrary.wiley.com/doi/10.1002/jgrg.20046). Heat coupling: [Chen 2019](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2019JD030324). CO2: [Humphrey 2018](https://ideas.repec.org/a/nat/nature/v560y2018i7720d10.1038_s41586-018-0424-4.html). Yield: [PMC6850578](https://pmc.ncbi.nlm.nih.gov/articles/PMC6850578/)

### Inferences (the hypotheses)

**H1. The "El Niño ratchet": global land TWS responds asymmetrically to ENSO. Strong El Niños remove water faster and more than La Niñas restore it, producing step declines (2015–16, 2023–24) that the 2026–27 El Niño will repeat. [Novel: asymmetry test not found; steps documented]**
- Data: global and continental land TWS (mask Greenland, Antarctica and glacier cells), Niño 3.4 monthly anomaly, ONI.
- Test: a distributed-lag regression of dS/dt on Niño3.4⁺ and Niño3.4⁻ separately, with lags of 0–9 months and a linear trend term. H0: β⁺ = β⁻, tested by block bootstrap (12-month blocks) plus surrogates.
- Pre-registration: before the Aug 2026 to 2027 COST-G fields appear, publish (e.g. as a dated git commit) the predicted continental TWS anomalies for Dec 2026, Mar 2027 and Jun 2027, with 90% intervals from the 2002–2025 model. Falsified if global land TWS falls by less than the lower bound or recovers within 12 months.
- Pitfalls: only 2–3 strong events, and the groundwater-depletion trend; include the trend in the null.

**H2. The 2014–16 drop is a genuine changepoint against a red-noise or long-memory null, and 2023–24 is a second one. [Gap: no formal test found]**
- Data: global non-ice land TWS (monthly).
- Test: likelihood-ratio or PELT changepoint on deseasonalised series with a trend. Calibrate the critical value with AR(1) and ARFIMA surrogates fitted to the 2002–2014 segment. Repeat excluding NE South America to see whether the step is global or regional, and test the Rateb 2011–12 date against the Rodell 2014–16 date.
- Pitfalls: the gap begins about 15 months after the step ends, and GRACE's final year had degraded accelerometer data (background knowledge, not verified here). Repeat on GRACE-FO-only data for 2023–24.

**H3. AR precipitation is more storage-efficient than non-AR precipitation: at equal precipitation anomaly, months with more AR days produce larger dS/dt in mid-latitude west-coast and snow-fed basins. [Novel: no global AR–TWS regression found]**
- Data: Guan–Waliser AR catalogue (ERA5; PANGAEA to 2021, v4 for later), aggregated to monthly AR-day counts per 2° cell; GPCP monthly P; COST-G dS/dt.
- Test: per region (California/Sierra, Pacific NW, central Chile, Norway/W Europe, UK, East Antarctic basins), fit dS/dt = a·P′ + b·AR′ + ε, cold-season only. H0: b = 0. Significance by permuting years within season, cross-validated by leaving one water year out.
- Pitfalls: AR count and P are collinear, so report variance inflation and partial R². Leakage dominates in small basins, so use basins larger than about 2×10⁵ km².

**H4. Pre-summer TWS deficit predicts summer hot-day counts beyond surface soil moisture and precipitation (deep-memory heat predictability). [Novel: no GRACE added-skill test for heat found]**
- Data: April–May TWS anomaly; ESA CCI surface soil moisture; GPCP; a daily Tmax product (CPC Global, Berkeley Earth or ERA5) giving JJA (DJF in the Southern Hemisphere) days above the 90th percentile.
- Test: per cell, the cross-validated R² gain of adding TWS to a model with soil moisture and precipitation, with FDR across cells. Expect positive gains in deep-rooted transitional zones (eastern South America, southern Africa, US Great Plains, Eastern Europe), per Chen 2019.
- Pitfalls: warming trends in both series, so detrend both. Irrigation regions give TWS declines with local cooling, so mask irrigated cells or test them separately.

**H5. Pre-fire-season TWS predicts burned area beyond the Amazon (boreal Canada and Siberia, SE Australia, Mediterranean, African savannas). [Novel extension of Chen 2013]**
- Data: TWS anomaly 1–5 months before each region's fire peak; GFED5 or MCD64A1 monthly burned area; GPCP; Niño 3.4.
- Test: a rank (Spearman) regression of log burned area on pre-season TWS, controlling for pre-season P and ENSO. Fit on 2002–2016 and test on 2018–2025. Pre-register 2026–27 Amazon and Australian fire seasons under the current El Niño.
- Pitfalls: about 24 fire seasons per region; human ignition and land-use trends; Canada 2023 alone could dominate, so use leave-one-out.

**H6. Out-of-sample replication of Humphrey 2018: the TWS–CO2 growth-rate sensitivity found in GRACE also holds in GRACE-FO (2018–2026), including the 2023–24 El Niño. [Replication with independent processing]**
- Data: global and tropical land TWS (detrended); NOAA GML monthly global CO2 growth rate; tropical land temperature.
- Test: regress the detrended CO2 growth rate on TWS and on temperature, fitted on 2002–2016 and predicted for 2018–2026, then compare the slope across periods with a Chow test.
- Pitfalls: CO2 growth-rate noise from smoothing, and volcanic or anthropogenic emission changes (COVID-19 2020).

**H7. ENSO diversity: central-Pacific and eastern-Pacific El Niños leave distinct, separable global TWS fingerprints, for example NE South America more sensitive to CP and the Yangtze to EP. [Gap: only regional studies]**
- Data: Niño 4 and Niño 3 (or E/C indices); TWS grids.
- Test: multiple regression maps on orthogonalised E and C indices at lags of 0–6 months. Use a global surrogate test for whether the two maps differ (pattern correlation of the maps against a surrogate distribution).
- Pitfalls: very few CP and EP events; report this as exploratory with power estimates.

**H8. A TWS "memory map": each region behaves as a leaky integrator of its dominant climate index, TWS(t) ≈ Σₖ β e^(−k/τ) I(t−k). τ is longest in groundwater-dominated drylands and East Africa (cumulative IOD, e.g. Lake Victoria 2019–20) and short in humid tropics. [Partly novel: hydraulic memory known for precipitation forcing, not for index forcing]**
- Test: a per-cell nonlinear fit of β and τ for ENSO and IOD (and NAO in Europe/Africa), compared with an instantaneous-regression model by cross-validated error. Map τ with bootstrap confidence intervals.
- Use: τ sets the maximum useful forecast lead (it links to H1, H4 and H5).

**H9. Replicate the 2026 Antarctic result: interannual East Antarctic (Queen Mary and Wilkes Land) mass flux lags Indo-Pacific warm-pool SST and AR landfall counts. [Single-study replication]**
- Data: COST-G basin mass flux (dM/dt, which removes linear GIA); warm-pool SST index (from ERSST; Niño 4 and DMI-east as proxies); Guan–Waliser Antarctic AR counts.
- Test: lagged regression with surrogates. Out of sample: does the relation fitted to 2002–2020 predict 2021–2023?
- Pitfalls: Antarctic GIA uncertainty (avoid trends), S2/161-day aliasing over West Antarctica (fit and remove it), and leakage from 300 km smoothing near the coast.

**H10. Pre-season TWS modulates crop-yield response to growing-season rain globally: wet antecedent storage amplifies excess-rain losses and buffers drought losses. [Extension of a US-only study]**
- Data: GDHY or FAOSTAT yields (maize, wheat, rice, soy); pre-season TWS; growing-season GPCP.
- Test: panel regression of yield anomaly on P, TWS and their interaction, with country and year effects and errors clustered by year.
- Pitfalls: technology trends in yield (detrend), and national aggregation blurring 2° signals.

**H11. (Expected null, worth a short test) Tropical land TWS (Amazon, Congo, Maritime Continent) contains information about ENSO's evolution 3–9 months later beyond ENSO's own persistence. [Novel; land→ocean feedback is physically weak]**
- Test: PCMCI+ or conditional Granger on Niño 3.4, regional TWS and WWV (warm water volume), cross-validated. A clean null would usefully bound land–atmosphere feedback.

**H12. (Control or "do-not-chase") Solar-cycle and QBO imprints on global land TWS. [Single weak claim, r = 0.12]**
- Run it as a methods control. With about two solar cycles, surrogate tests should show no detectable signal. A "significant" result would most likely signal a flaw in the pipeline (trend leakage, gaps), so it is useful as a negative control alongside a deliberately scrambled index.

### Gaps
- I could not verify whether H3, H4, H5, H7 or H8 have been done in literature published after my searches or outside the indexed results. Run a targeted literature search (Google Scholar, Web of Science) for each before investing effort.
- Power analysis is needed before H1, H7 and H9, because few strong ENSO events and about 24 seasons mean modest effects (r ≈ 0.3) may be undetectable.
- GIA is not removed in our data. H2 (steps) and anything quoting trend magnitudes in glaciated or rebounding regions need a GIA model before results are reported.
