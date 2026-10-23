# Satellite-gravimetry mass change vs solid-Earth, Earth-rotation, geomagnetic and human-activity phenomena (GRACE/GRACE-FO 2002–2026)

Scope note: web research only, done 2026-10-24. The WebSearch budget ran out before the last planned searches, so several well-known papers that I could not check this session appear under **Gaps** as "unverified leads" and not as findings. "Established", "single-study" and "contested" labels are given in the Inferences sections. The data context is the project's 259 COST-G monthly fields (2002-04 to 2026-07, degree 90), the TN-13/TN-14 low-degree replacements, and our own Rust SH synthesis.

## 1. Seismicity: does hydrological/seasonal (and GRACE-measured) loading modulate earthquake rates?

### Takeaway
Seasonal hydrological loading modulates small-earthquake rates at the few-kPa level in several well-instrumented regions: California, the Himalaya, Taiwan, New Madrid and SE Tibet. The effect on large earthquakes is not established. Most studies use GPS-constrained or regional load models. GRACE is usually a supporting input, and I found no uniform global GRACE-based test.

### Cited Findings
- **California (established, GPS-based):** Johnson, Fu & Bürgmann (2017) used 9 years of GPS vertical deformation to model monthly hydrospheric loading and its stress on small-earthquake fault planes. Earthquakes occur more often when the load stress favours rupture, so seismicity rates are "modestly modulated" by hydrological loading. Science 356:1161–1164, doi:10.1126/science.aak9547 — [SCEC record](https://central.scec.org/node/5455)
- Annual peak-to-peak hydrological stresses on northern California faults can exceed 1 kPa. The largest occur on dipping reverse faults in the Coast Ranges and along the eastern Sierra front — [SCEC poster 2016](https://central.scec.org/node/5456)
- A NASA-hosted abstract comparing 7 seasonal loading sources (2006–2014) found hydrology to be the largest source of seasonal stress in California. It found a modest increase in microseismicity timing with an ~8 kPa seasonal mean-normal-stress perturbation — [NASA ESD abstract](https://esdpubs.nasa.gov/node/69024)
- Johnson says there is no increase in large earthquakes from these low-amplitude stresses, only a correlation with small ones — [UC Berkeley news](https://news.berkeley.edu/?p=63090)
- **Himalaya (established for micro-seismicity):** Bettinelli et al. (2008) found Nepal seismicity in winter to be about twice the summer rate. They attributed the geodetic seasonal signal mostly to monsoon water storage, using GRACE-derived land-water mass. EPSL 266:332–344, doi:10.1016/j.epsl.2007.11.021 — [USGS OF 2010-1099 abstract](https://pubs.usgs.gov/of/2010/1099/chanard/of2010-1099_chanard.pdf); [CNIG presentation](https://cnig.gouv.fr/IMG/documents_wordpress/2018/GEOPOS_2018_KChanard.pdf) (the presentation itself lists "does GRACE-derived seasonal stress correlate with seismicity?" as an open question)
- Panda et al. (2018) estimated hydrological-load stress perturbations on the Main Himalayan Thrust of ~2 kPa, ~0.02% of a typical stress drop. Shear-stress perturbations are ~3× larger on the mid-crustal ramp than on the flat, and Coulomb stress peaks in May–June. They propose "fault resonance" in which seasonal stress modulates slow slip on the ramp, which then modulates micro-seismicity. Nature Communications 9:4140, doi:10.1038/s41467-018-06371-2 — [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC6175945/)
- Ader & Avouac (2013) developed Schuster-spectrum periodicity and declustering methods for Himalayan catalogues. EPSL 377:97–105, doi:10.1016/j.epsl.2013.06.032 — [cited in Panda et al. 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6175945/)
- **Mixed evidence:** A 2004–2015 regional study found that Mw 3–4 earthquakes in the Eastern Himalaya increase significantly as terrestrial water storage (TWS) declines. It found insufficient evidence for long-term rate changes in other magnitude bins and zones — [IIGM repository](https://repository.iigm.res.in/bitstreams/8cdb8e13-7d5e-4e8f-acb6-cab8b16577aa/download)
- Geodesy and Geodynamics (2022) used GRACE and GNSS to assess seasonal loading and seismicity in Nepal — [DOAJ](https://doaj.org/article/0f40782c17fa4ba7a3a01680aa1c241c)
- **SE Tibet (GRACE-constrained):** She, Fu & Xu (2022), "Seasonal terrestrial water load modulation of seismicity at the southeastern margin of the Tibetan Plateau constrained by GNSS and GRACE data", GJI 230(3):1966–1979, doi:10.1093/gji/ggac168 (results not retrieved) — [GFZ GRACE bibliography](https://www-app2.gfz-potsdam.de/grace-references/b2hd-2022GeoJI.230.1966S.html)
- **Taiwan:** Hsu et al. (2021) report that western-Taiwan seismicity peaks in Feb–Apr and is lowest in Jul–Sep, in phase with annual water unloading. Deep eastern events behave the same way, but shallow (<18 km) eastern events are anti-phased and not explained by loading or pore pressure. The authors infer modestly higher hazard in low-load seasons. Science Advances 7:eabf7282, doi:10.1126/sciadv.abf7282 — [eScholarship PDF](https://escholarship.org/content/qt0h27b1fx/qt0h27b1fx.pdf); [Academia Sinica summary](https://www.earth.sinica.edu.tw/en/research/fields/detail/28)
- **New Madrid (stable continental interior):** Craig, Chanard & Calais (2017) show annual and multi-year micro-earthquake rates tracking water loading of the upper Mississippi embayment, through small elastic crustal stresses. Nature Communications 8, doi:10.1038/s41467-017-01696-w — [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5732252/)
- GRACE month-to-month variability is "highly contaminated with outliers, measurement errors and uncertainties" (Humphrey et al. 2016, citing Bonin et al. 2012). This matters for seasonal-stress studies — [PMC4944666](https://pmc.ncbi.nlm.nih.gov/articles/PMC4944666/)

### Inferences
- **Established:** seasonal modulation of *small* earthquakes by hydrological load at the ~1–10 kPa level in specific regions. **Not established:** any modulation of M≥6 timing, or any global relation.
- Most positive studies drive their stress models with GPS loading inversions or regional models. GRACE enters as a cross-check or for the long-wavelength part. A globally uniform test is the obvious gap: GRACE-derived load → elastic stress on GCMT focal planes → Schuster/phase test with one declustering method and a per-region null. The scale mismatch is real, since GRACE sees ~300 km while stress is evaluated on km-scale faults. For a smooth surface load, though, the elastic stress at seismogenic depth is dominated by long wavelengths, which partly justifies using GRACE.
- The Taiwan result, with opposite phase at different depths, warns that a single global "load-phase" statistic can cancel real signals. Any global test should stratify by depth and faulting type.

### Gaps
- Japan snow-load seasonality (Heki 2003, EPSL) and Japan Sea-side seasonal seismicity: not verified this session (search budget exhausted).
- No global GRACE-based seasonal-seismicity study found. Absence is not proof, because only a few searches covered this.
- I could not check whether Bettinelli's 2:1 winter/summer ratio survives the Ader & Avouac declustering; the abstract was not retrieved.

## 2. Coseismic, postseismic, pre-seismic, slow-slip and volcanic gravity signals

### Takeaway
Coseismic and postseismic gravity from Mw≳8.3 megathrust events is established, about ten events in the GRACE era plus the 2025 Mw 8.8 Kamchatka event in GRACE-FO. The Mw 7.8 Turkey 2023 sequence is below the usual detection threshold. The pre-Tohoku gravity "precursor" is contested: a single group reports it, and an independent statistical test (Wang & Bürgmann 2019) found it not unique. I found no published GRACE detection of a slow-slip event or of an eruption's own mass.

### Cited Findings
- **First 2-D coseismic detection:** Han et al. (2006) measured a ±15 µGal change from Sumatra-Andaman 2004. Vertical displacement of density interfaces and crust/mantle dilatation contributed roughly equally at GRACE wavelengths. Science 313:658–662, doi:10.1126/science.1128661 — [GFZ record](https://www-app2.gfz-potsdam.de/grace-references/b2hd-2006Sci...313..658H.html); [AGRIS abstract](https://agris.fao.org/search/en/records/65de59f70f3e94b9e5ce1272)
- **Event list:** about 20 Mw≥8.0 earthquakes occurred 2002–2017. GRACE detected regional co- and postseismic changes for Sumatra-Andaman 2004 (Mw 9.2), Nias 2005 (8.5), Kuril 2006/07 (8.3), Bengkulu 2007 (8.5), Samoa-Tonga 2009 (8.1), Maule 2010 (8.8), Tohoku-Oki 2011 (9.0), Wharton Basin 2012 (8.6) and Okhotsk 2013 (8.3, deep ~600 km). These include megathrust, strike-slip and normal faulting — [AGU FM15 abstract](https://agu.confex.com/agu/fm15/webprogram/Paper85405.html)
- **Postseismic pattern (Heki group):** for the three M9-class events, gravity drops coseismically, keeps dropping for a few months, then rises over longer times. These are two components with different time constants and polarities, tentatively afterslip and Maxwell viscoelastic relaxation — [Hokkaido Univ. repository](https://eprints.lib.hokudai.ac.jp/dspace/handle/2115/67855). Cumulative postseismic change can exceed the coseismic change — [AGU FM15](https://agu.confex.com/agu/fm15/webprogram/Paper85405.html)
- Tohoku-Oki coseismic gravity is detectable with only one month of post-event data, but standard destriping distorts coseismic signals — [Earth Planets Space 2015](https://earth-planets-space.springeropen.com/articles/10.1186/s40623-015-0188-0)
- Spherical-harmonic solutions recover mostly long-wavelength coseismic signal and are hit harder by striping. Mascons are more stable but attenuate coseismic signals — [Applied Sciences 2026](https://www.mdpi.com/2076-3417/16/9/4066)
- **Detection threshold:** GRACE range data should allow source-parameter estimation for Mw≥8.0 — [NASA GSFC geodesy page](https://science.gsfc.nasa.gov/earth/geodesy/researchareas/121). The Mw 8.2 Chignik (2021) event is "near the coseismic detection threshold" for GRACE-era missions. The Kuril 2006/07 and Samoa 2009 events had large, long-wavelength postseismic signals — [ESSOAR preprint](https://essopenarchive.org/doi/full/10.1002/essoar.10510424.1)
- **Kamchatka, 29 July 2025, Mw 8.8:** GRACE-FO monthly fields show a coseismic signal of "nearly 20 µGal" peak-to-peak at mission resolution. Braitenberg, Maurizio, Javed & Panet invert it for dislocation on the seismological fault plane, and read any excess over the seismic estimate as aseismic slip. Detection is "challenging due to the relatively high noise level". This is an EGU26 abstract only, with no filtering, solution or slip values given — [EGU26-15981](https://meetingorganizer.copernicus.org/EGU26/EGU26-15981.html). The magnitude is listed as Mw 8.7–8.8, and other sources give up to 8.9 — [Wikipedia](https://en.wikipedia.org/wiki/2025_Kamchatka_earthquake); [RJES](https://rjes.ru/ru/nauka/article/114798/view)
- **Turkey 2023 (Mw 7.8 + 7.7):** I found no peer-reviewed GRACE-FO coseismic detection. An arXiv preprint (Mousavi, Mohammadigheymasi & Moradi, June/Aug 2026) claims pre-rupture "systematic increase in entropy measures alongside a gradual decrease in EWH". Its abstract gives no anomaly magnitudes and no significance tests — [arXiv 2606.17938](https://arxiv.org/abs/2606.17938). The rupture was ~300 km on the East Anatolian Fault plus ~180 km on secondary faults — [He et al. 2023 GRL via DOAJ](https://doaj.org/article/016c2f618bc3477a9e1dfa7e5b3b4217)
- **Earthquake signals as nuisance:** Bonin et al. (2026) built an empirical earthquake model (iterative PCA on CSR mascons) for Tohoku 2011 and the Andaman-Sumatra events of 2004, 2007 and 2012. Removing it cut near-epicentre ocean-mascon RMS from >50 cm to <5 cm. Earth and Space Science 13(3), e2025EA004219 — [GFZ record](https://www-app2.gfz.de/grace-references/b2hd-2026ES...1304219B.html). A companion ocean mascon product with earthquakes and barystatic/GRD removed is in ESSD 18:3481–3505 (2026) — [ESSD](https://essd.copernicus.org/articles/18/3481/2026/essd-18-3481-2026-metrics.html)
- **Pre-seismic claim (contested):** Panet, Bonvalot, Narteau, Remy & Lemoine (2018) reported gravity-gradient changes around Japan a few months before Tohoku-Oki and interpreted them as large-scale deep slab motion. Nature Geoscience 11:367–373, doi:10.1038/s41561-018-0099-3 — [GFZ record](https://www-app2.gfz.de/grace-references/b2hd-2018NatGe..11..367P.html); [Eos](https://eos.org/research-spotlights/months-of-gravity-changes-preceded-the-tohoku-earthquake)
- Wang & Bürgmann (2019) built a hypothesis test for short-term gradient anomalies at any time and place. They concluded the proposed precursor is not unique in time or space, so it does not support accelerated subduction, and they suggest temporally correlated GRACE errors or other processes. GRL 46:7323–7332, doi:10.1029/2019GL082682 — [GFZ record](https://www-app2.gfz-potsdam.de/grace-references/b2hd-2019GeoRL..46.7323W.html); abstract via [Semantic Scholar API](https://api.semanticscholar.org/graph/v1/paper/DOI:10.1029/2019GL082682?fields=title,abstract,year,venue)
- Panet et al. (2022) reply with a cross-model consistency method. They find the pre-Tohoku dipolar anomaly "nearly unique in space and time" in GRGS03, and the pre-seismic increase "highly consistent" between GRGS03 and CSR06 regardless of ocean dealiasing. JGR Solid Earth 127:e2022JB024542 — [IRD record](https://www.documentation.ird.fr/hor/fdi:010085978); [GFZ record](https://www-app2.gfz-potsdam.de/grace-references/b2hd-2022JGRB..12724542P.html). A GSTM2022 abstract compares Tohoku with Maule 2010 and proposes systematic monitoring of the Pacific belt — [GSTM2022-76](https://meetingorganizer.copernicus.org/GSTM2022/GSTM2022-76.html)
- A JpGU 2019 attempt reproduced pre-seismic gradients "fairly similar" to Panet et al. 2018, "although we still see certain differences due to unknown reasons" — [JpGU 2019](https://confit.atlas.jp/guide/event/jpgu2019/subject/E_SGD01-P03/date). Commentators asked whether the signal is unique to Tohoku or needs an M9 to be seen — [Earth Magazine](https://www.earthmagazine.org/article/gravity-changes-portended-2011-tohoku-earthquake/)
- **Slow slip:** gravity is suggested as "a promising tool to reveal silent earthquakes", but I found no demonstrated GRACE detection — [IPGP seminar](https://www.ipgp.fr/actus-et-agenda/agenda/seminaires/can-geodynamics-benefit-from-satellite-gravity-the-case-study-of-andaman-sumatra-2004-earthquake)
- **Volcanoes:** a 2024 Surveys in Geophysics study finds a future quantum-gravity mission could detect deep-ocean volcanic mass changes of ≳6.5 km³. Hunga Tonga–Hunga Ha'apai and Fani Maoré meet that threshold, and the authors find GRACE insufficient — [Surv. Geophys. 2024](https://link.springer.com/article/10.1007/s10712-024-09839-7); [GFZ record](https://www-app2.gfz.de/grace-references/b2hd-2024SGeo...45.1331B.html)

### Inferences
- **Established:** co- and postseismic GRACE signals for Mw≳8.3 shallow megathrust events, and the two-time-constant postseismic pattern. **Single-group/contested:** the pre-Tohoku precursor. Its originators find it robust across solutions; an independent test finds it not unique. The community has not converged, and the claim should be treated as doubtful until tested on other events and on independent combined products such as COST-G.
- The Kamchatka 2025 event (Mw 8.8, oceanic setting with low hydrological noise) has about 12 months of post-event GRACE-FO data in our COST-G set (to 2026-07). Only an EGU abstract exists so far, which makes it the most time-sensitive novel target.
- The Turkey 2023 "precursor" preprint has no significance testing, and the event is below the ~Mw 8 threshold. Flag it as doubtful.
- Slow-slip detection with monthly GRACE is implausible for typical Mw 6–7 slow-slip events, which sit below the coseismic threshold. The only realistic GRACE "slow slip" signal is aseismic excess slip after great earthquakes, as the Kamchatka abstract proposes.

### Gaps
- Matsuo & Heki (2011, Tohoku coseismic, GRL) and Heki & Matsuo (2010, Maule, GRL): unverified leads; DOIs not checked this session.
- Illapel 2015 (Mw 8.3) and South Sandwich 2021 (Mw 8.1) GRACE/GRACE-FO detections: not found.
- The full Wang & Bürgmann statistics (false-alarm rate) and the Panet 2022 Maule result: not retrieved (paywall).
- A 2020 study reportedly found GRACE signals at only the Icelandic volcanoes (Eyjafjallajökull, Grímsvötn) among 13 VEI≥4 eruptions, attributed to glacier response. I could not identify the paper, so it is not cited.

## 3. Earth rotation: hydrological/cryospheric excitation of polar motion and LOD

### Takeaway
That surface mass transport, much of it GRACE-measured, explains most of 2003–2015 polar motion and the post-2000 eastward drift of the pole is established. Groundwater pumping and glacier loss explain the 1990s drift change, but rely on models before 2002. Climate-driven LOD lengthening (~1.33 ms/century for 2000–2018) is now a published single-group result. Whether GRACE improves operational LOD prediction remains largely untested in the literature I found. GRACE-FO-era work links the Chandler wobble's post-2015 amplitude anomaly to hydrology and cryosphere.

### Cited Findings
- Adhikari & Ivins (2016): TWS plus cryosphere changes explain ~83 ± 23% of the observed polar motion amplitude (2003–2015) and match the mean direction to 5.9° ± 7.6°. Water storage alone accounts for the decadal-like variability. Since ~2000 the spin axis has drifted toward Greenwich, ~75° east of its long-term direction. Science Advances 2:e1501693, doi:10.1126/sciadv.1501693 — [PMC4846461](https://pmc.ncbi.nlm.nih.gov/articles/PMC4846461)
- Deng et al. (2021): the polar drift changed from southward to eastward in 1995, and drift speed for 1995–2020 was ~17× that of 1981–1995. Glacier loss is the main driver and groundwater a contributor; pre-2002 TWS comes from glacier and pumping estimates, not GRACE. GRL 48:e2020GL092114, doi:10.1029/2020GL092114 — [GFZ record](https://www-app2.gfz.de/grace-references/b2hd-2021GeoRL..4892114D.html); [AGU release](https://news.agu.org/press-release/climate-has-shifted-the-axis-of-the-earth/)
- Seo et al. (2023): a model with ~2,150 Gt of groundwater redistribution (1993–2010; a prior model estimate) matches observed polar drift. Without it the model misses by 78.5 cm, or 4.3 cm/yr. The water corresponds to >6 mm sea-level equivalent, and mid-latitude withdrawals (western North America, NW India) matter most. GRL, doi:10.1029/2023GL103509 — [AGU release](https://news.agu.org/press-release/weve-pumped-so-much-groundwater-that-weve-nudged-the-earths-spin)
- Kiani Shahvandi et al. (2024, PNAS, 15 July 2024): ice and groundwater mass change lengthened the day by 0.3–1.0 ms/century over 1900–2000 and by 1.33 ms/century over 2000–2018, using GRACE/GRACE-FO plus earlier mass-balance studies. Under high emissions this could reach 2.62 ms/century by 2100, above the Moon's tidal ~2.4 ms/century — [JPL news](https://www.jpl.nasa.gov/news/nasa-funded-studies-explain-how-climate-is-changing-earths-rotation/)
- Kiani Shahvandi et al. (2024, Nature Geoscience, 12 July 2024): with machine learning on the 120-year record, ~90% of recurring polar motion fluctuations (1900–2018) are attributed to groundwater, ice and sea level, and the rest to interior dynamics — [JPL news](https://www.jpl.nasa.gov/news/nasa-funded-studies-explain-how-climate-is-changing-earths-rotation/); [ETH press](https://ethrat.mironet.ch/research/how-climate-change-is-altering-the-earths-rotation.html). A 2026 JGR Solid Earth follow-up extends this to the Late Pliocene, finding the current lengthening rate exceptional over 3.6 Myr (press summary only) — [earth.com](https://earth.com/news/earth-days-getting-longer-rotation-slows-climate-change)
- **Chandler wobble:** Shi, Zhou, Chen & Xu (2025, J. Geodesy) report anomalous behaviour after 2015, with near-disappearance and re-excitation. Using GRACE (2002-04 to 2015-12) and GRACE-FO (2018-06 to 2024-12), they find GRACE/GFO captures hydrological and cryospheric Chandler-band forcing better than models, though discrepancies remain at times — [Springer](https://link.springer.com/article/10.1007/s00190-025-02021-w); [GFZ record](https://www-app2.gfz.de/grace-references/b2hd-2025JGeod..99...97S.html)
- Jin, Chambers & Tapley (2010, JGR): TWS and ocean bottom pressure are the major polar-motion excitations after atmosphere. GRACE fits best at annual period (Px); ECCO/ECMWF models fit better intraseasonally — [USF digital commons](https://digitalcommons.usf.edu/msc_facpub/190)
- Nastula & Gross (2015): Chandler period 430.9 ± 0.7 solar days and Q = 127 (56–255), from GRACE/SLR and models — [NASA ESD record](https://esdpubs.nasa.gov/sgpubs/node/69386)
- **LOD/EOP prediction:** GFZ combined atmosphere, ocean and hydrology angular-momentum forecasts reach 10-day errors of 3.02 mas (PM) and 5.39 mas (ΔUT1), 34.5% and 44.7% better than Bulletin A — [GFZ publication](https://gfzpublic.gfz-potsdam.de/pubman/item/item_3247891_5/component/file_4001888/3247891.pdf). ETH ML extension of the EAM forecasts improves ultra-short-term LOD prediction by up to 19% — [J. Geod. 2023](https://link.springer.com/article/10.1007/s00190-023-01745-x). Hydrological model excitations "are generally deemed uncertain". An ESD supplement compares GRACE+SLR+DORIS hydrological χ with the COST-G Level-3 product and GLDAS — [ESD 2024 supplement](https://esd.copernicus.org/preprints/esd-2024-21/esd-2024-21-AC1-supplement.pdf)
- **Low-degree caveat:** GRACE-FO C20 and C30 are sensitive to the GRACE-D accelerometer "transplant" approach — [EGU2020-3017](https://meetingorganizer.copernicus.org/EGU2020/EGU2020-3017.html); [Behzadpour 2021 JGR via TU Graz](https://repository.tugraz.at/records/b8kfr-fp535)

### Inferences
- **Established:** GRACE-era polar motion is dominated by surface mass transport (multiple groups). **Single-group but plausible:** quantitative climate LOD rates and the ML attribution over 120 years. **Model-dependent:** the pre-2002 groundwater attribution.
- Our COST-G set extends Adhikari & Ivins' 2003–2015 budget to 2026-07, covering the post-2015 Chandler anomaly and the 2019–2024 hydrological extremes. Polar-motion mass excitation comes from ΔC21/ΔS21 and LOD mass excitation from ΔC20 (textbook relations). For a closure test, ΔC20 must come from SLR (TN-14) and degree-1 from TN-13. The atmosphere/ocean mass removed by AOD1B dealiasing must be added back (GAC), and motion terms come from model EAM (e.g., GFZ ESMGFZ).
- GRACE latency (months) makes it useless for operational LOD forecasting directly. It is more promising as calibration or bias correction of model hydrological angular momentum (HAM), and that use appears untested.

### Gaps
- DOIs for both Kiani Shahvandi 2024 papers, Agnew (2024, Nature, global warming delaying a negative leap second), and Chen/Wilson's C20–LOD comparisons: unverified leads.
- IERS EOP 20 C04 documentation could not be fetched; the IERS page returned 404 and hpiers.obspm.fr refused the connection. The series is standard, but its citation and URL must be confirmed separately.
- I found no study that isolates the effect of GRACE-derived HAM on LOD prediction skill.

## 4. Core and geomagnetism: can GRACE see core flow, jerks, decadal LOD or inner-core rotation?

### Takeaway
Claims exist and come mostly from one IPGP/CNES/Strasbourg network: Mandea et al. 2012, the Lecomte et al. 2021 poster, and Gaugne Gouranton et al. 2025–2026 on a January 2007 Eastern-Atlantic anomaly coinciding with the 2007 jerk. The same network's 2023 uncertainty study concludes that extracting a core signal is "very difficult", and an independent expert (Dumberry) notes the observed signal is ≥10× the expected core signal. Status: contested, not established.

### Cited Findings
- Mandea, Panet, Lesur, de Viron, Diament & Le Mouël (2012) report shared interannual variability between geomagnetic secular acceleration (CHAMP) and GRACE (CNES/GRGS) gravity over a region centred on Africa, from the Atlantic to the Indian Ocean. PNAS 109:19129–19133, doi:10.1073/pnas.1207346109 — [GFZ record](https://www-app2.gfz-potsdam.de/grace-references/b2hd-2012PNAS..109..19129M.html); [ScienceDaily](https://sciencedaily.com/releases/2012/10/121022145248.htm)
- Lecomte et al. presented a time-varying pattern linked to fluid-core motion and gravity–magnetic correlations at AGU FM 2021 (poster). Dumberry (not involved) said ≥90% of the measured signal is surface processes, that the gravity change is ≥10× what the core is expected to produce, and that the correlation could be coincidence — [Eos, Dec 2021](https://eos.org/articles/can-nasas-gravity-satellites-detect-motions-in-earths-core)
- Lecomte, Rosat, Mandea, Boy & Pfeffer (2023): GRACE solutions disperse by ~0.34 cm EWH (~20% of total low-degree signal), hydrological models by 0.89–2.10 cm. Residual hydrology in small regions can contaminate low degrees, so extracting core information is very difficult. JGR Solid Earth 128:e2023JB026503 — [GFZ record](https://www-app2.gfz-potsdam.de/grace-references/b2hd-2023JGRB..12826503L.html); [ESSOAR preprint](https://essopenarchive.org/doi/full/10.22541/essoar.167768139.90331914/v1)
- Lecomte, Rosat, Mandea & Dumberry (2023), inner-core differential rotation from SLR+GRACE/GRACE-FO S22: upper bounds of 0.09°, 0.3° and 0.4° at ~4, 6 and 12-year periods (for 90 m inner-core topography). These are likely overestimates because S22 is still hydrology-contaminated. GRL 50:e2023GL104790 — [DOAJ](https://doaj.org/article/ba63cd986f8941a898d65fc646885c7d); [GFZ record](https://www-app2.gfz-potsdam.de/grace-references/b2hd-2023GeoRL..5004790L.html)
- Gaugne Gouranton, Panet, Mandea, Greff-Lefftz & Rosat (EGU26-12420): a rapid north–south gravity anomaly over the Eastern Atlantic in January 2007, found by pattern recognition, evolving over months to years, "at least partly" solid-Earth in origin and coincident with the 2007 geomagnetic jerk. Hypothesised cause: decimetric CMB topography from perovskite→post-perovskite transition shifts near the African LLSVP. No amplitudes or statistics given, and no peer-reviewed paper cited — [EGU26-12420](https://meetingorganizer.copernicus.org/EGU26/EGU26-12420.html). An EGU25 abstract uses GRACE + SLR 2003–2015 to constrain CMB mass redistribution — [EGU25-8808](https://meetingorganizer.copernicus.org/EGU25/EGU25-8808.html). A poster states that the inferred CMB topography changes are too large to be explained by core-flow dynamic pressure — [EGU 2025 OSPP poster](https://www.egu.eu/media/awards/posters/2025/09/15/52c03157-5cf7-490a-831d-693b7aee7442.pdf)
- **Decadal LOD and jerks:** Holme & de Viron (2013, Nature) identified a ~5.9-year LOD oscillation (~0.1 ms) whose extrema or features coincide with geomagnetic jerks within ±6 months. They suggest fluid-core motions or inner-core coupling as origins — [authors' manuscript](https://www.liverpool.ac.uk/~holme/nature_sub.pdf); [Science News](https://www.sciencenews.org/article/every-six-years-earth-spins-slightly-faster-and-then-slower). A 2020 analysis gives P = 5.85 ± 0.06 yr with no damping over 1962–2019 and says the oscillation's stability is disputed — [arXiv 2007.08252](https://arxiv.org/pdf/2007.08252)

### Inferences
- **Contested/doubtful:** any GRACE core-flow detection. The strongest new claim (Jan 2007 anomaly) exists only as conference abstracts. **Testable with our data:** whether a January-2007 north–south Eastern-Atlantic step appears in COST-G (a different processing chain from GRGS) after removing GLDAS/WGHM hydrology and ocean models, and whether similar steps follow later jerks.
- The S22 inner-core bound can be updated with three extra years (2023–2026). Hydrology contamination will probably dominate, so expect only a marginal tightening unless hydrology is removed with multiple models and the spread is carried as uncertainty.

### Gaps
- Dates of post-2007 geomagnetic jerks (e.g. 2011, 2014, 2017, 2020) and the seismological inner-core rotation results (Yang & Song 2023 Nat. Geosci.; Wang & Vidale 2024 Nature): not verified this session.
- The Holme & de Viron 2013 DOI was not checked (believed to be 10.1038/nature12282; unverified).

## 5. Human activity: groundwater depletion, dams, irrigation, lakes, mining and induced seismicity

### Takeaway
GRACE detection of large aquifer depletion (Middle East, India, California, North China Plain) and the attribution framework (Rodell et al. 2018, 34 trends) are established. Single-dam detection is marginal: for GERD I found only simulation and feasibility studies. Links from pumping-induced unloading to earthquakes are single-event or model-based (Lorca 2011; Gorkha 2015; Central Valley/Parkfield; Morocco 2023) and remain debated.

### Cited Findings
- Rodell et al. (2018) quantified 34 GRACE TWS trends (2002–2016) and attributed them to natural variability, unsustainable groundwater use, climate change or combinations. Some, such as NW China and the Okavango, had lacked attribution. Antarctica lost −127.6 ± 39.9 Gt/yr and Greenland −279.0 ± 23.2 Gt/yr. Nature 557:651–659, doi:10.1038/s41586-018-0123-1 — [NASA NTRS PDF](https://ntrs.nasa.gov/api/citations/20180003130/downloads/20180003130.pdf); [IDEAS](https://ideas.repec.org/a/nat/nature/v557y2018i7707d10.1038_s41586-018-0123-1.html)
- Voss et al. (2013): Tigris–Euphrates–western Iran storage fell −27.2 mm/yr (143.6 km³, Jan 2003–Dec 2009), mainly groundwater. About 60% came from pumping and was linked to the 2007 drought. This was the second-fastest loss after India. WRR 49:904–914, doi:10.1002/wrcr.20078 — [UCAR](https://impacts.ucar.edu/en/publications/groundwater-depletion-in-the-middle-east-from-grace-with-implicat/); [AGU release](https://news.agu.org/press-release/groundwater-levels-drop-at-alarming-rate-in-large-swath-of-middle-east)
- Chao, Wu & Li (2008): ~10,800 km³ is impounded in reservoirs since ~1900, offsetting ~30 mm of sea-level rise (~0.55 mm/yr over 50 years). Correcting for it gives a near-constant ~2.46 mm/yr rise over 80 years. Science 320:212–214, doi:10.1126/science.1154580 — [preprint PDF](https://ftp.soest.hawaii.edu/coastal/Climate%20Articles/Chao_2008%20Water%20impoundment.pdf); [CU sea level](https://sealevel.colorado.edu/index.php/bibcite/reference/277)
- **GERD:** Chen et al. (2021) simulated GERD-induced deformation for 5/10/15-year fills of ~70 km³ and examined the limits of GRACE-FO detection (conclusions not retrieved). Remote Sensing 13:874, doi:10.3390/rs13050874 — [Crossref](https://api.crossref.org/works/10.3390%2FRS13050874). Abdelmohsen et al. used GRACE/GRACE-FO, Sentinel-1/2, MODIS and altimetry to quantify GERD seepage; the planned reservoir is 1,904 km² and 74 km³ — [GSTM2022-68](https://meetingorganizer.copernicus.org/GSTM2022/GSTM2022-68.html)
- **Lorca 2011 (Mw 5.1):** shallow (2–4 km) slip on the Alhama de Murcia Fault coincides with positive Coulomb stress from groundwater extraction (water table down ~250 m since 1960). González et al., Nature Geoscience 5:821–825, doi:10.1038/NGEO1610 — [Earth-prints](https://www.earth-prints.org/handle/2122/8270); [LiveScience](https://www.livescience.com/24161-earthquake-triggered-by-humans.html). González: "Even without the groundwater extraction, the earthquake was overdue" — [ScienceNews](https://sciencenews.org/?p=9077)
- **Indo-Gangetic depletion and Himalaya:** Kundu et al. (2015) argue that GRACE-measured unloading probably influenced Gorkha 2015 and other Main Himalayan Thrust events. GRL 42:10,607–10,613 — [GFZ record](https://www-app2.gfz.de/grace-references/b2hd-2015GeoRL..4210607K.html). Kundu's later wording is "may have modestly advanced the earthquake clock". The same group estimates groundwater unloading at ~10–14% of tectonic stress accumulation for the 2023 Morocco earthquake, and says more work is needed — [SSA news](https://www.seismosoc.org/news/groundwater-depletion-may-have-been-significant-trigger-for-2023-morocco-earthquake/)
- **Central Valley:** Amos et al. (2014, Nature) find 1–3 mm/yr uplift around the southern San Joaquin Valley from groundwater loss (~160 km³ over >150 yr). The resulting reduction in San Andreas normal stress is proposed as a mechanism for Parkfield microseismicity seasonality. Amos called the stress changes "exceedingly small" — [UC Berkeley news](https://news.berkeley.edu/?p=41787); [Earth Magazine](https://www.earthmagazine.org/article/humans-causing-californias-mountains-grow/)
- Groundwater pumping also measurably moved the pole (Seo et al. 2023, Section 3) — [AGU release](https://news.agu.org/press-release/weve-pumped-so-much-groundwater-that-weve-nudged-the-earths-spin)

### Inferences
- **Established:** GRACE detection of basin-scale depletion. **Single-event or contested:** depletion-triggered earthquakes. The press often overstates them, for example "probably all earthquakes" — [Gulf Times](https://gulf-times.com/story/481838/Depletion-of-groundwater-linked-to-Nepal-earthquak).
- GERD's planned 74 km³, filled in annual monsoon steps in 2020–2024, is a known-timing, known-location artificial mass load. That makes it a natural ground-truth test for an abrupt-change detector. The detectability question (signal vs Blue Nile seasonal variability) is open in the literature I found.
- The rotational effect of reservoir impoundment (Chao 2008 scale ~10,800 km³) is resolvable in principle in C20/C21 trends. Single dams like GERD are not.

### Gaps
- Rodell et al. 2009 (NW India, Nature), Famiglietti et al. 2011 (Central Valley, GRL), Feng et al. 2013 (North China Plain, WRR) and the Three Gorges GRACE studies: unverified leads.
- Actual GERD fill volumes per year (altimetry/Sentinel) were not retrieved.
- Mining and lake drainage/filling (e.g., Caspian Sea decline): no sources gathered.
- Whether any study links *GRACE-observed* depletion trends to a *measured* change in seismicity rate (as opposed to a modelled stress change) was not found.

## 6. Ocean and solid-Earth coupling: GIA, ocean bottom pressure, tidal aliasing, ice-dynamics events

### Takeaway
GIA is the dominant systematic uncertainty in Antarctic GRACE mass balance. Low-viscosity regions (Amundsen Sea) need ~10% upward revision and a correction that is itself time-variable. S2 tidal aliasing at ~161 days remains the most persistent tide alias, strongest around Antarctica. Large earthquakes are a major contaminant of GRACE ocean-bottom-pressure products.

### Cited Findings
- Barletta et al. (2018, Science): ~41 mm/yr uplift in the Amundsen Sea Embayment implies mantle viscosity ~4×10¹⁸ Pa·s, a GIA response time of decades to a century, and a ~10% upward revision of gravimetric ice loss — [Penn State record](https://pure.psu.edu/en/publications/observed-rapid-bedrock-uplift-in-amundsen-sea-embayment-promotes-/). The low-viscosity GIA rate in the Amundsen sector exceeds +15 Gt/yr and is accelerating, which could mask ice-loss acceleration — [EGU2018-10349](https://meetingorganizer.copernicus.org/EGU2018/EGU2018-10349.pdf)
- Forward and inverse GIA models systematically underestimate uplift in low-viscosity, thin-lithosphere regions — [PMC5111427](https://pmc.ncbi.nlm.nih.gov/articles/PMC5111427)
- A geologically constrained GIA model gives −69 ± 18 Gt/yr for Antarctica (2002–2010); allowing for GIA and destriping errors widens this to −126 to −29 Gt/yr. Nature, doi:10.1038/nature11621 — [Springer](https://link.springer.com/article/10.1038/nature11621)
- **Tidal aliasing:** S2 errors alias to ~161 days and appear near Antarctic ice shelves, Greenland and Hudson Bay. RL04 reduced but did not remove them, and C20 carried enlarged S2 errors (Chen, Wilson & Seo 2009, J. Geod.) — [KOPRI](https://repository.kopri.re.kr/handle/201206/6084). An RL06 evaluation (2002–2017) of the 161.0-day S2, 171.2-day P1 and 322.1-day S1 aliases finds S2 mostly in West Antarctica, with an eastward-circling pattern of about 2 years. Mascons show unusually large alias error, attributed to regularisation. Remote Sensing 14:5403, doi:10.3390/rs14215403 — [MDPI](https://www.mdpi.com/2072-4292/14/21/5403)
- Simulated S2 errors are ~3× measurement noise at n<15, while K1, O1 and M2 average down — [CSR GSTM abstract](https://www2.csr.utexas.edu/grace/GSTM/abstracts/abstract_submissions/ray_1127425201.html). AOD1B dealiasing is documented in TN-06 — [PO.DAAC TN-06](https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-docs/grace/open/docs/TN-06_AOD1B.pdf)
- Earthquake signals dominate ocean mascons near epicentres (>50 cm RMS before correction) — [Bonin et al. 2026 record](https://www-app2.gfz.de/grace-references/b2hd-2026ES...1304219B.html)

### Inferences
- For any solid-Earth test near Antarctica, Greenland or Hudson Bay, 161-day and 322-day aliases are a confounder that must be fitted alongside annual and semi-annual terms. GRACE-FO alias periods differ because of the orbit, so they should be recomputed rather than copied from GRACE.
- A time-variable (accelerating) GIA correction in West Antarctica is an under-used idea. Most processing still subtracts a constant GIA rate.

### Gaps
- Ice-shelf collapse events (Larsen B 2002, Conger 2022) and GRACE-seen glacier-acceleration responses: no sources gathered.
- Sea-level fingerprint validation with GRACE OBP and validation of ocean-bottom-pressure recorders: not researched this session.

## 7. Space weather and atmosphere: any credible geophysical link, and nuisance effects on GRACE

### Takeaway
I found no credible evidence that space weather causes geophysical mass change seen by GRACE. Solar activity matters as a nuisance: thermospheric density rises about two orders of magnitude toward solar maximum, and storms such as May 2024 (G5) stress orbit and drag modelling. GRACE-FO's GRACE-D accelerometer failed in 2018, and the transplant data used in its place affect C20 and C30 in particular, with harder force modelling near the solar cycle 25 maximum.

### Cited Findings
- GRACE-D accelerometer noise rose by orders of magnitude in June 2018, with common-mode errors on all axes. The SDS "transplant" derives GRACE-D accelerations from GRACE-C with time and attitude corrections plus modelled thruster responses — [EGU2020-3017](https://meetingorganizer.copernicus.org/EGU2020/EGU2020-3017.html); [Behzadpour et al. 2021 JGR via TU Graz](https://repository.tugraz.at/records/b8kfr-fp535)
- Low-degree zonals C20 and C30 are sensitive to the transplant. TU Graz's remove-restore transplant improves C20/C30 and reduces the SLR–GRACE-FO C30 offset — [TU Graz alternative L1B](https://tugraz.at/institute/ifg/downloads/gravity-field-models/alternative-grace-fo-l1b). Rising solar activity near the cycle 25 peak makes force-model transplants harder, and thruster leakage affects both satellites — [EGU24-10401](https://meetingorganizer.copernicus.org/EGU24/EGU24-10401.html)
- GRACE-FO accelerometer-derived densities (2018–2024) rise by ~2 orders of magnitude toward solar maximum, with storm disturbances mainly at high latitudes — [EGU25-20153](https://meetingorganizer.copernicus.org/EGU25/EGU25-20153.html)
- The May 2024 (Gannon) storm reached G5, the most severe since October 2003. Atmospheric density is the main uncertainty in storm-time orbit work (GRACE-FO, TerraSAR-X) — [UCL Discovery](https://discovery.ucl.ac.uk/id/eprint/10204526). GRACE-FO and CHAMP accelerometer densities serve as truth for 45 storms in a density-method comparison — [arXiv 2408.16805](https://arxiv.org/abs/2408.16805v1)

### Inferences
- **Not credible / no evidence found:** any space-weather → geophysics mass link. **Credible nuisance:** solar-activity-dependent noise in GRACE-FO monthly fields, especially C20, C30 and possibly short-wavelength noise. This should be a control variable in every detector and precursor test, because months near solar maximum (2023–2025) overlap the Kamchatka window.

### Gaps
- I found no study quantifying monthly Level-2 noise or degree-variance vs F10.7/Kp, or degradation in May 2024 specifically.

## 8. Detection methods for abrupt events in gravity time series, and validation catalogues

### Takeaway
GRACE earthquake work mostly uses event-specific step fitting at known times and places, or empirical PCA models (Bonin 2026). The pre-seismic debate centres on hypothesis tests for anomaly uniqueness (Wang & Bürgmann 2019; Panet 2022). I found no published blind, global change-point detector for GRACE validated against a catalogue with a reported false-alarm rate. Methods exist in GPS offset detection and GRACE hydrology breakpoints.

### Cited Findings
- Change-point detection (CPD) works offline and online for automatic GPS offset detection — [SciTePress 2019](https://scitepress.org/Papers/2019/79515/)
- A GRACE TWS breakpoint method uses EMD to find extreme-return periods and then a statistical test for abrupt-change times (hydrology, not seismic) — [Annals of Geophysics 2014 via DOAJ](https://doaj.org/article/b20cf7413dd84e2bb0e7888cf3bab7a1)
- Bonin et al. (2026) use iterative PCA to build an empirical earthquake model for removal — [GFZ record](https://www-app2.gfz.de/grace-references/b2hd-2026ES...1304219B.html)
- Wang & Bürgmann (2019) built a hypothesis test for anomaly significance "at any time and place" — [GFZ record](https://www-app2.gfz-potsdam.de/grace-references/b2hd-2019GeoRL..46.7323W.html). Panet et al. (2022) use cross-solution spatial-overlap consistency — [IRD](https://www.documentation.ird.fr/hor/fdi:010085978)
- Destriping distorts coseismic signals — [EPS 2015](https://earth-planets-space.springeropen.com/articles/10.1186/s40623-015-0188-0)
- **Validation catalogue, USGS ComCat FDSN service:** `https://earthquake.usgs.gov/fdsnws/event/1/query` with starttime/endtime, min/maxmagnitude, rectangle or circle, depth, and orderby, in quakeml/geojson/csv/kml/text formats. The hard cap is 20,000 events per query (HTTP 400 above it); there is also a `count` method — [USGS FDSN event docs](https://earthquake.usgs.gov/fdsnws/event/1/)
- Known GRACE-detected Mw≥8 event list for positives — [AGU FM15](https://agu.confex.com/agu/fm15/webprogram/Paper85405.html)

### Inferences
- A detector can be validated against three catalogues. Seismic positives: ComCat Mw≥8 2002–2026, including Kamchatka 2025 and Chignik 2021. Anthropogenic positives: GERD annual fills 2020–2024, with timing from Sentinel/altimetry. Hydrological extremes (floods, droughts) count as "true but non-seismic" detections that must be classified, not counted as false alarms. The null comes from random time/location permutation, as in Wang & Bürgmann, and from cross-centre consistency, as in Panet 2022. COST-G is itself a multi-centre combination, so cross-centre consistency is partly built in.
- Candidate statistics: a per-pixel or per-Slepian-region step+exp/log template matched filter (the postseismic shape is known from Heki's two-time-constant model), CUSUM/PELT change-point with seasonal + 161-day + 322-day nuisance terms, and a solar-activity covariate.

### Gaps
- No published global false-alarm rate for GRACE abrupt events was found. That gap is the opportunity.

## 9. Ranked, testable hypotheses using the 259 COST-G monthly fields plus public data

### Takeaway
The most novel and feasible tests are, in order: (1) the 2025 Kamchatka Mw 8.8 co- and postseismic signal, (2) a blind, catalogue-validated global abrupt-change detector with an empirical Mw threshold and false-alarm rate, and (3) a falsification test of the pre-seismic-gravity claim across all great earthquakes in an independent combined product. Rotation-budget closure through 2026, a uniform global hydro-stress/seismicity test, and core-signal reproducibility follow. Run the solar-activity noise control first.

### Cited Findings
- Only an EGU abstract exists for the Kamchatka GRACE-FO signal (~20 µGal p-p) — [EGU26-15981](https://meetingorganizer.copernicus.org/EGU26/EGU26-15981.html)
- The precursor is contested (Panet 2018/2022 vs Wang & Bürgmann 2019) — [GFZ](https://www-app2.gfz.de/grace-references/b2hd-2018NatGe..11..367P.html); [GFZ](https://www-app2.gfz-potsdam.de/grace-references/b2hd-2019GeoRL..46.7323W.html)
- The polar-motion budget is published to 2015 (Adhikari & Ivins) and the Chandler anomaly to 2024 (Shi et al. 2025) — [PMC4846461](https://pmc.ncbi.nlm.nih.gov/articles/PMC4846461); [GFZ](https://www-app2.gfz.de/grace-references/b2hd-2025JGeod..99...97S.html)
- The January 2007 core-related anomaly is conference-only — [EGU26-12420](https://meetingorganizer.copernicus.org/EGU26/EGU26-12420.html)

### Inferences (the ranked list; rank = novelty × feasibility × falsifiability)

**0. Control (run first): solar-activity and transplant noise in COST-G.**
H0: High-degree noise (degree-variance n>40 over open ocean) and C20/C30 scatter in COST-G GRACE-FO months correlate with F10.7/Kp and with transplant-only periods. Data: COST-G, TN-14, NOAA F10.7, GFZ Kp. Test: regress ocean-pixel RMS residual (after an annual + semi-annual + 161/322-day fit) on monthly F10.7, with an HAC (heteroskedasticity- and autocorrelation-consistent) standard error. The output is a per-month noise weight used by H1–H3.

**1. Kamchatka 2025 (Mw 8.8) co- and postseismic gravity, and the aseismic excess.** (novel: abstract only; feasible: ~12 post-event months)
H1: COST-G shows a coseismic step consistent with the USGS finite-fault model. The step is followed by the Heki-type pattern: further decrease over a few months, then recovery. Gravity-inferred moment exceeds the seismic moment, implying afterslip.
Data: COST-G 2002-04..2026-07, USGS finite-fault and ComCat aftershocks, ocean model (AOD1B GAD), TN-13/14.
Test: regional (≈2,500 km radius) gravity disturbance at n≤60 without destriping (or with a matched, documented filter). Fit seasonal + alias terms on 2018–2025-06. Fit a step + exponential/log term at 2025-07/08. Do a matched-filter amplitude against the forward model (spherical, layered-Earth dislocation) through the same filter. Use bootstrapped null amplitudes from the same months at 100 random oceanic locations. Falsified if the step is not above the 95th percentile of the null.

**2. Blind global abrupt-change detector with empirical Mw threshold and false-alarm rate.** (novel: no published global CPD with FAR found)
H2: A matched-filter/CUSUM scan of COST-G recovers all Mw≥8.5 shallow events (Sumatra 2004, Nias 2005, Bengkulu 2007, Maule 2010, Tohoku 2011, Wharton 2012, Kamchatka 2025) blind, at a false-alarm rate below X per year. Detection probability drops below 50% at Mw≈8.2–8.3, testing the Chignik "near threshold" claim. GERD fills 2020–2024 are not detected (or are detected only cumulatively).
Data: COST-G, ComCat (FDSN query with minmagnitude=7.5), GCMT mechanisms, GERD fill timing from altimetry/Sentinel.
Test: ROC curve; detection probability vs Mw (logistic fit); null by time-shuffling 1,000 times.

**3. Pre-seismic gravity falsification test across great earthquakes.** (contested claim; new data product)
H3: In COST-G, pre-event gradient anomalies within 6 months before Tohoku 2011, Maule 2010, Illapel 2015, Chignik 2021 and Kamchatka 2025 are not more frequent or stronger than at random times and places along the same subduction zones.
Data: COST-G and ComCat. Use gravity gradients via our SH synthesis.
Test: replicate the Wang & Bürgmann-style permutation test. Add the Panet 2022 consistency criterion by splitting COST-G vs individual-centre solutions (CSR, GFZ, JPL, ITSG from ICGEM). Pre-register the window and threshold before looking. A positive result in all five events would be a major finding; a null result strengthens the "doubtful" label.

**4. Rotation budget closure 2002–2026, including the post-2015 Chandler anomaly.**
H4a: COST-G ΔC21/ΔS21 (+GAC + model motion terms) closes geodetic polar-motion excitation (IERS EOP 20 C04) at the annual period and for 2016–2024 Chandler-band amplitude changes. Hydrology and cryosphere account for the re-excitation, extending Shi et al. 2025 to 2026-07 with an independent product.
H4b: SLR/TN-14 ΔC20 trends 2002–2026 give a climate LOD contribution consistent with ~1.33 ms/century (2000–2018) and show acceleration after 2018.
Data: COST-G, TN-13/14, AOD1B GAC, GFZ ESMGFZ EAM (motion terms), IERS C04.
Test: complex-band-pass χ comparison, coherence and residual spectrum. The geodetic-minus-surface residual is the input to H8.

**5. Uniform global hydro-load Coulomb stress vs seismicity phase.**
H5: Across all GCMT M≥4.5 events (2002–2026) with mechanisms, the phase of GRACE-derived load Coulomb stress shows significant clustering (Schuster p<0.01) in regions of large seasonal TWS amplitude (monsoon Asia, Amazon margins, Taiwan, Himalaya). Results are stratified by depth (Taiwan anti-phase lesson) and faulting type. Regions with low TWS amplitude show none.
Data: COST-G TWS (degree 60, Love numbers already in our gravimetry crate), GCMT/ComCat, declustered by Zaliapin–Ben-Zion or Reasenberg.
Test: elastic half-space or spherical Love-number stress at the hypocentre for the GCMT plane, then a per-region Schuster test with FDR (false-discovery-rate) control. Compare against a GPS-based published result (California, Johnson 2017) as a positive control.

**6. Human-made mass steps as ground truth (GERD, Three Gorges) and their rotational fingerprint.**
H6: GERD cumulative storage (2020–2024) is detectable in COST-G only after removing Blue Nile hydrology with a model ensemble. The detection SNR scales with the cumulative fill, which gives an empirical minimum-detectable point mass for COST-G at degree 90.
Data: COST-G, GLDAS/WGHM, Hydroweb/DAHITI altimetry for reservoir level, Sentinel-1/2 area.
Test: compare a point-mass forward model through the same filter with regional averages; use a seasonal-climatology null.

**7. Groundwater depletion trend vs secular seismicity-rate change.**
H7: Where GRACE-derived depletion produces positive long-term Coulomb stress rates on mapped faults (NW India/Himalayan front, Central Valley margins, North China Plain, Tigris–Euphrates/Zagros), declustered M≥Mc rates trend upward relative to paired control regions. This tests Kundu-, Amos- and González-type claims statistically rather than case by case.
Data: COST-G trends, regional catalogues (ComCat/ISC), fault databases.
Test: difference-in-differences of rates, using completeness-magnitude checks per period.

**8. Core/jerk gravity reproducibility (high novelty, high risk).**
H8: The January 2007 north–south Eastern-Atlantic gravity step (EGU26-12420) appears in COST-G at comparable amplitude and is not explained by the hydrology/ocean model ensemble. Analogous steps follow later jerks. The geodetic-minus-surface polar-motion residual from H4 shows a feature at the same time.
Data: COST-G (independent of GRGS), GLDAS/WGHM/ocean models, CHAOS-7 geomagnetic secular acceleration.
Test: regional pattern-matched step detection with the H2 detector, an uncertainty envelope from the model spread (Lecomte 2023 numbers), and cross-correlation with secular acceleration maps. Expect a likely null given Dumberry's ≥10× objection, but a null in an independent product is itself informative.

**9. Inner-core S22 bound update (low expected gain).**
H9: Adding 2023–2026 COST-G S22 tightens Lecomte et al.'s (2023) upper bounds (0.09°/0.3°/0.4°) only marginally, because hydrology dominates.
Test: their method plus multi-model hydrology removal.

### Gaps
- The ranking is my judgement, not a published one. Any of H1–H3 could already be in press. Re-search before committing, especially for Kamchatka 2025 and COST-G-based precursor tests.
- The H0 control's actual effect size is unknown. If solar-max months are much noisier, H1 sensitivity drops.
- IERS EOP 20 C04 and GFZ ESMGFZ EAM data documentation URLs were not verified this session.
