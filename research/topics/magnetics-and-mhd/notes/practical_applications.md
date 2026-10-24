# Practical applications of the crate family, with published numbers to reproduce

Scope: who uses each kind of calculation, for what decision and at what scale, plus published numeric cases (inputs and the published answer) that a short worked example could reproduce. Findings are grouped by crate. Within each crate, applications are ranked by how directly the crate's own functions reproduce the published number (1 = closed-form or table check; higher = needs external data or models). Company figures are labelled **[claim]**. Dates are given where the source has them. Research done 2026-10-24.

## gravimetry: how is satellite gravimetry used, and what published numbers can a worked example reproduce?

### Takeaway
GRACE/GRACE-FO mass change feeds real decisions: weekly US drought maps, California groundwater (SGMA) reports, India's water planning, WWF's finance risk tool, and the ice-sheet numbers behind sea-level assessments. The numbers easiest to reproduce are basin mass trends in Gt/yr, mm/yr to km³ conversions, and coseismic gravity of about ±15 µGal. Drought percentiles and flood lead times need a land model as well.

### Cited Findings

**Who uses it (operational users, 2024–2026)**
- US Drought Monitor (NDMC, USDA, NOAA, NASA): GRACE-FO data feed a NASA hydrologic model that makes weekly groundwater and soil-moisture indicators, "a key input" to the weekly USDM maps. California DWR has put GRACE-FO Central Valley groundwater estimates in its semi-annual public reports since 2024, under the Sustainable Groundwater Management Act. Texas runs a STRAWS dashboard with the Texas Water Development Board. India's Ministry of Water Resources used GRACE-informed assessments (2016 and 2025 releases). NGS uses GRACE for the geodetic datums (DGEOID2022). The World Bank Climate Knowledge Portal and the WWF Water Risk Filter ("used by industry and finance") also use it. A 2022 GRACE-FO user survey names farmers, commodity interests, energy producers, insurers and supply-chain consultants as known users — [NASA JPL GRACE applications](https://grace.jpl.nasa.gov/applications/overview/)
- The GRACE-DA groundwater drought indicator is a percentile against a 1948–2012 climatology (other product versions say 1948–2014), on a 0.125° grid over North America, issued weekly (Mondays) from April 2002 — [Drought.gov](https://www.drought.gov/data-maps-tools/groundwater-and-soil-moisture-conditions-grace-data-assimilation); [nasagrace.unl.edu](https://nasagrace.unl.edu)

**Ranked applications with reproducible numbers**

1. **Ice-sheet mass loss and its sea-level contribution (most direct: monthly fields to mascon or basin sums to Gt/yr to mm of sea level)**
   - IMBIE (Otosaka et al., ESSD 15:1597–1616, April 2023), combining 50 satellite estimates for 1992–2020: the two ice sheets raised global mean sea level by 21.0 ± 1.9 mm. The loss rate rose from 105 Gt/yr (1992–1996) to 372 Gt/yr (2016–2020). Greenland averaged 169 ± 9 Gt/yr, ranging from 86 Gt/yr (2017) to 444 Gt/yr (2019). West Antarctica lost 82 ± 9 Gt/yr, the Antarctic Peninsula 13 ± 5 Gt/yr, and East Antarctica gained 3 ± 15 Gt/yr — [NOAA repository record of Otosaka et al. 2023](https://repository.library.noaa.gov/view/noaa/69781)
   - From GRACE/GRACE-FO, Greenland lost about 264 Gt/yr over 2002–2023, about 0.8 mm/yr of sea-level rise — [NASA SVS](https://svs.gsfc.nasa.gov/31156/). Vital Signs gives 265 Gt/yr since 2002 — [earth.gov Vital Signs](https://earth.gov/sealevel/vital-signs/greenland-ice-sheet/). The JPL mascon trend for 04/2002–04/2024 is −267.49 ± 21.00 Gt/yr — [NASA SVS data file](https://svs.gsfc.nasa.gov/webapps/sealevelgraphs/Updated/data/greenland_mass_200204_202404.txt). This file is the best single target for a worked example, because inputs and output come from one product.
   - Sea-level budget: an ensemble of 90 GRACE/GRACE-FO processing variants gives a global barystatic (ocean-mass) rate of about 2.00–2.45 mm/yr ± 0.27 mm/yr (95%) for Jan 2005–Dec 2019. The choice of GIA model moved budget closure by 0.20–0.30 mm/yr. This is from a search summary of [Ocean Science 19:321 (2023)](https://os.copernicus.org/articles/19/321/2023/); I did not confirm which paper reports the 90-variant figure.

2. **Groundwater depletion in named basins (direct: equivalent water height in mm/yr times area gives km³/yr; 1 km³ of water = 1 Gt)**
   - Middle East (Tigris–Euphrates–western Iran; Voss et al. 2013, *WRR* 49:904): total water storage fell 143.6 km³ from Jan 2003 to Dec 2009, −27.2 mm/yr equivalent water height. Groundwater fell 17.3 ± 2.1 mm/yr, or 91.3 ± 10.9 km³, about 60% of the loss and the second-fastest groundwater loss on Earth after India — [AGU press release, 12 Feb 2013](https://news.agu.org/press-release/groundwater-levels-drop-at-alarming-rate-in-large-swath-of-middle-east); [UCAR dataset](https://data.ucar.edu/dataset/groundwater-depletion-in-the-middle-east-from-grace-with-implications-for-transboundary-water-m)
   - California Central Valley (Liu, Famiglietti et al., *Nat. Commun.* Dec 2022): groundwater loss of 1.86 km³/yr over 1961–2021, 2.41 km³/yr over 2003–2021 and 8.58 km³/yr over 2019–2021, nearly five times the 60-year average — [PMC9763392](https://pmc.ncbi.nlm.nih.gov/articles/PMC9763392). Earlier GRACE total-water result (2011): the Sacramento and San Joaquin basins lost more than 30 km³ after October 2003, and the San Joaquin lost 3.5 km³/yr. That figure is total water, not groundwater alone — [NASA JPL news](https://www.jpl.nasa.gov/news/nasa-data-reveal-major-groundwater-loss-in-california/)
   - Northern India: Rodell et al. 2018 (*Nature*, "Emerging trends in global freshwater availability") quantified 34 GRACE trends for 2002–2016. A secondary summary gives northern India at 19.2 ± 1.1 Gt/yr, attributed to irrigation pumping — [Springer Nature "Behind the paper"](https://communities.springernature.com/posts/measuring-changes-in-freshwater-availability-from-space-a-brief-and-narrow-history); [TABLE debates](https://www.tabledebates.org/node/9450). An earlier GRACE study found a net loss of 109 km³, about double India's largest surface reservoir — [UN-SPIDER](https://www.un-spider.org/news-and-events/news/grace-satellites-show-global-scarcity-groundwater); [NASA JPL](https://www.jpl.nasa.gov/news/satellites-unlock-secret-to-northern-indias-vanishing-water/). A figure of "54 km³/yr" comes from a different study and should not be credited to Rodell.

3. **Earthquake gravity (direct with the Okada/Okubo functions plus degree truncation and filtering)**
   - 2004 Sumatra–Andaman (Han et al., *Science* 313:658–662, 2006, doi:10.1126/science.1128661): GRACE saw a ±15 µGal gravity change. Because GRACE sees only long wavelengths, vertical displacement and crustal dilatation contributed about equally — [Ohio State news](https://news.osu.edu/satellite-data-reveals-gravity-change-from-sumatran-earthquake/); [AGRIS abstract](https://agris.fao.org/search/en/records/65de59f70f3e94b9e5ce1272)
   - 2011 Tohoku-Oki: Matsuo & Heki 2011 (*GRL* 38, L00G12) report back-arc gravity decreases from crustal dilatation, scaling linearly with seismic moment — [Hokkaido repository](https://eprints.lib.hokudai.ac.jp/dspace/handle/2115/48676). Wang et al. 2012 (*GRL*) used Slepian localisation to estimate a total moment of about 4.59 × 10²² N·m (Mw 9.07) through July 2011 — [Princeton PDF](https://geoweb.princeton.edu/people/simons/PDF/reprints/GRL-2012a.pdf). A 2015 study found that standard destriping can distort coseismic amplitude and pattern, which is a test case for the crate's stripe filter — [Earth, Planets and Space](https://earth-planets-space.springeropen.com/articles/10.1186/s40623-015-0188-0)

4. **Gap filling (direct test for the Kalman filler)**
   - There is an 11-month GRACE/GRACE-FO gap, July 2017 to May 2018 (GRACE ended Oct 2017; GRACE-FO launched May 2018) — [J. Geodetic Science 2022](https://www.degruyterbrill.com/document/doi/10.1515/jogs-2022-0129/html). Forward–backward autoregressive bridging has been published — [EGU22-111](https://meetingorganizer.copernicus.org/EGU22/EGU22-111.html). A worked example: withhold 11 known months from GRACE, fill them, and score against truth.

5. **Reservoir storage**
   - Three Gorges, Apr 2002–May 2010: GRACE minus the WGHM model gave the total volume gain to within 14% (3.2 km³) of in-situ data. Monthly changes explained 76% of in-situ variance (r² = 0.76) with an uncertainty of 4.62 km³ — [eScholarship](https://escholarship.org/uc/item/2n63b635). The planned volume was about 40 km³ by 2009 — [EGU2008 abstract](https://www.meetings.copernicus.org/www.cosis.net/abstracts/EGU2008/06276/EGU2008-A-06276.pdf)

6. **Flood potential (indirect: needs a runoff model)**
   - Reager, Thomas & Famiglietti 2014 (*Nat. Geosci.* 7:588, doi:10.1038/ngeo2203): for the 2011 Missouri River floods, GRACE storage could have stretched regional warning lead time from 2 to up to 5 months. For the 2011 Columbia floods, warning was possible about 3 months ahead. Flash floods are out of reach — [UCI press release](https://news.uci.edu/press-releases/satellites-reveal-possible-catastrophic-flooding-months-in-advance-uci-finds); [JPL flood-potential page](https://grace.jpl.nasa.gov/applications/flood-potential). Context: Rocky Mountain snowpack was about 212% of normal — [2011 Missouri River flood](https://en.wikipedia.org/wiki/2011_Missouri_River_flood)

7. **Atmospheric pressure waves and sensor noise (direct but simple)**
   - Hunga Tonga 2022 Lamb wave: about 314 ± 3 m/s. The peak was 780 Pa at the nearest station (757 km), and more than 100 Pa was recorded worldwide (preprint) — [ESSOAr preprint](https://essopenarchive.org/doi/full/10.1002/essoar.10511401.1). Other estimates give 318 m/s — [Wright et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9492537/). The wave made four passages around Earth (plus three antipodal) over 6 days — [Newswise on the *Science* paper](https://www.newswise.com/articles/as-big-as-it-gets-hunga-volcano-comparable-to-krakatoa)
   - Cold-atom accelerometer noise floors in m s⁻² Hz⁻¹ᐟ²: SYRTE lab 5e-8; commercial iXblue/Muquans 5e-7; CARIOQA pathfinder target 1e-10; future quantum space gravimetry 1e-12 — [Lévèque et al. 2022](https://arxiv.org/pdf/2211.01215)

### Inferences
- Voss et al. check: 143.6 km³ / 7 yr / 0.0272 m/yr and 91.3 km³ / 7 yr / 0.0173 m/yr both give a region of about 754,000 km². The two published numbers agree, so a worked example can recover the basin area from them (computed).
- Greenland at 264 Gt/yr ÷ 361.8 Gt per mm (ocean area about 3.618 × 10⁸ km²) gives 0.73 mm/yr. NASA's "≈0.8 mm/yr" therefore uses rounding or a different ocean area. IMBIE's 21 mm corresponds to about 7,600 Gt (computed; the ocean area is a standard value, not taken from these sources).
- The ±15 µGal Sumatra signal is the obvious first test for the Okada/Okubo code at GRACE resolution (degree 60 or less with a Gaussian filter), because the published value is quoted at that resolution.
- Insurance and finance use GRACE through intermediaries (WWF Water Risk Filter, consultants). No source showed GRACE used directly in pricing.

### Gaps
- No primary source found for an insurance or parametric product triggered by GRACE data. A Financial Protection Forum piece discusses research only ([financialprotectionforum.org](https://financialprotectionforum.org/node/809)).
- The search summary gave the Matsuo & Heki 2011 Tohoku back-arc decrease as "~7 mGal or less", which is implausibly large for GRACE. The abstract probably says µGal. Check the PDF before quoting.
- A 2019 Greenland record of "532 ± 58 Gt/yr" appeared only in a search summary without an identifiable primary source.
- No hydropower or agriculture case with numbers was found, beyond the user categories listed on the JPL page.
- The Frederikse et al. 2020 barystatic value was not retrieved.

## geopotential: how are the geoid and gravity models used, and which named-place numbers can be reproduced?

### Takeaway
Geoid models turn GNSS ellipsoidal heights into usable heights: NAPGD2022/GEOID2022 in the US, AUSGeoid2020 with GDA2020 in Australia, and EGM96/EGM2008 worldwide. The best worked-example targets are GeographicLib's published EGM96 test value (Timbuktu, 28.7068 m) and the −106 m Indian Ocean geoid low. Applied gravity (INS deflection compensation, quantum gravimeters, gradiometers, tunnel detection, volcano and reservoir monitoring) gives µGal and Eötvös-level targets.

### Cited Findings

1. **Geoid height at a named place (direct: spherical-harmonic synthesis from ICGEM coefficients)**
   - GeographicLib GeoidEval with the default EGM96 5′ grid: `echo 16:46:33N 3:00:34W | GeoidEval` (Timbuktu) returns 28.7068 m geoid height, with north and east slopes of −0.02e-6 and −1.73e-6 — [GeoidEval man page](https://manpages.ubuntu.com/manpages/focal/man1/GeoidEval.1.html)
   - The world's largest geoid low, south of India, reaches about −106 m — [Srinivas et al. 2012, J. Earth Syst. Sci. (repository copy)](https://repository.iigm.res.in/bitstreams/9a52a220-4176-4087-8e70-5973e2dd045a/download). EGM2008 is complete to degree and order 2159, with extra coefficients to degree 2190 — [GMT EGM2008 dataset page](https://www.generic-mapping-tools.org/remote-datasets/earth-geoid.html). In Thailand, EGM2008 against GPS/levelling showed 1.012 m RMS with a mean offset of 0.997 m (local datum offset) — [ACRS 2012](https://acrs-aars.org/proceeding/ACRS2012/Proceeding%20ACRS%202012/Technical%20Sessions/G1%20Photogrammetry%20(4)/G1-2.pdf)

2. **National height systems (partly direct: the geoid part is computable, the old-datum part needs NGS data)**
   - US: NAPGD2022 replaces NAVD 88 with four products: GM2022, GEOID2022 (gravimetric and time-dependent), DEFLEC2022 and GRAV2022. Using the experimental xGeoid18b, the NAPGD2022 − NAVD 88 differences are −0.722 m in California, −0.274 m in Louisiana, −0.646 m in Michigan and −0.307 m in Maine. At the Monroe, NC CORS (NCMR) the difference is −0.26 m (174.40 m vs 174.66 m) — [GPS World](https://www.gpsworld.com/napgd2022-is-coming-soon-how-will-the-new-geopotential-datum-impact-products-and-services/). The 2023 NGS plan put the beta in mid-2024 to mid-2025 and the official announcement in mid-2025 — [FGCS Feb 2023](https://www.fgdc.gov/organization/working-groups-subcommittees/fgcs/meetings/feb-2023/2-nsrs-modernization-efforts-fgcs-20230207.pdf). The release has slipped — [Oregon Surveyor 2023](https://associationpublications.com/flipbook/plso/2023/MarApr/17)
   - Australia: AUSGeoid2020 gives ellipsoid-to-AHD separation onshore with uncertainty. Its geometric component absorbs an AHD-versus-quasigeoid offset of about −0.5 m in the southwest to +0.5 m in the northeast. Against independent GNSS-levelling the standard deviation is 27 mm. It must be paired with GDA2020, while AUSGeoid09 goes with GDA94 — [Geoscience Australia](https://www.ga.gov.au/geodesy/ausgeoid/files.jsp); [Brown et al. 2018, J. Geod.](https://link.springer.com/article/10.1007/s00190-018-1202-7); [Janssen & Watson 2019](https://www.spatial.nsw.gov.au/__data/assets/pdf_file/0003/222474/2019_Janssen_and_Watson_Azimuth586_AUSGeoid2020.pdf)

3. **Inertial-navigation gravity compensation (direct: deflection of the vertical from EGM2008 coefficients)**
   - Gravity-induced INS errors sit almost entirely in the horizontal channels and come mostly from deflection of the vertical (DOV). Compensating from EGM2008 cut rugged-terrain horizontal error by up to 48.85% of the regional maximum. The paper recommends computing the gravity vector within 1 s and updating it at least every 100 s for marine navigation — [Sensors 2016, PMC5191156](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5191156/)
   - DOV is typically a few arcseconds except in rugged terrain — [ISPRS](https://WWW.ISPRS.org/proceedings/XXXV/congress/comm5/papers/678.pdf). DOV-induced INS position error reached about 3 km on one route, and gravity-aided navigation degraded beyond about 9.6 nmi of error (after 44 h in simulation) — [MATEC 2018](https://www.matec-conferences.org/articles/matecconf/pdf/2018/19/matecconf_eecr2018_07006.pdf). Stored DOV grids for navigators are patented — [US6249745](https://patents.google.com/patent/US6249745)

4. **Quantum and absolute gravimeters (direct: noise and averaging-time arithmetic)**
   - Exail AQG-A and AQG-B **[datasheet]**: sensitivity 500 nm/s²/√Hz (≈50 µGal/√Hz). Precision at a quiet site is 50 nm/s² in 1.5 min, 20 nm/s² in 10 min and 10 nm/s² in 40 min. Long-term stability ≤20 nm/s², trueness ≤150 nm/s², repeatability ≤50 nm/s², cycle rate >1 Hz, power 350 W (A) or 500 W (B), operating range 18–24 °C (A) or 0–40 °C (B). Named deployments are Etna (since August 2020), BKG (2025 *J. Geod.*), and the WET-CAG2021 comparison at Wettzell — [Exail](https://www.exail.com/product/quantum-gravimeters). The first unit (AQG01) reached 2 µGal after 1000 s — [EGU2016-11824](https://opac.geologie.ac.at/ais312/dokumente/EGU/EGU2016-11824.pdf). There is "µGal-level" reproducibility across 16 units — [arXiv:2405.10844](https://www.arxiv.org/pdf/2405.10844)
   - Etna: an AQG-B about 2.5 km from the active craters gave the first absolute atom-gravimeter time series on an active volcano, tracking bulk mass changes from magma and gas — [Antoni-Micollier et al., GRL 2022](https://www.earth-prints.org/bitstream/2122/15317/5/Geophysical%20Research%20Letters%20-%202022.pdf). The NEWTON-g project added a MEMS gravimeter array anchored to the AQG — [CORDIS](https://cordis.europa.eu/article/id/430166-using-gravity-to-explore-volcanoes)

5. **Volcano, reservoir and CO₂ monitoring (direct forward model: point mass or prism to µGal)**
   - Kīlauea, 2009–2015: 1.9 × 10¹¹ kg accumulated about 1.3 km below sea level beneath the caldera, attributed to magma densification — [USGS, Koymans et al. 2022](https://www.usgs.gov/publications/microgravity-change-during-2008-2018-kilauea-summit-eruption-nearly-a-decade). After the 2018 collapse (Nov 2018 to Mar 2019), gravity rose by up to 44 µGal, modelled as mass about 1 km deep filling void space — [KNMI listing of Poland et al.](https://www.knmi.nl/research/publications/post-collapse-gravity-increase-at-the-summit-of-kilauea-volcano-hawai-i)
   - Sleipner CO₂ storage, North Sea: seafloor gravity repeatability was 4.3 µGal (2003) and 3.5 µGal (2005). The predicted time-lapse signal was 3–16 µGal/yr — [OSTI](https://www.osti.gov/biblio/841869); [OSTI 902109](https://www.osti.gov/biblio/902109)

6. **Civil engineering and exploration gradiometry (direct: gradient of a buried cylinder or prism, in Eötvös)**
   - Stray et al. 2022 (*Nature* 602:590): a cold-atom gravity gradiometer surveyed an 8.5 m line at 0.5 m spacing over a multi-utility tunnel at Birmingham. It detected the tunnel at SNR 8, with the tunnel centre estimated at about 1.9 m depth. The team spun out Delta-g — [Physics World](https://physicsworld.com/a/quantum-gravity-gradient-sensor-used-outdoors-to-find-tunnel/); [Physics Today](https://physicstoday.aip.org/news/quantum-sensor-identifies-underground-tunnel)
   - Bell Geospace **[company claim]**: more than 2 million line-km of airborne full-tensor gradiometry acquired and more than 150 interpreted projects — [Bell Geospace](https://info.bellgeo.com/minex). 1 Eö = 0.1 mGal/km — [patent 7065449](https://patents.justia.com/patent/7065449). Early FALCON trials ran at about 30–40 E/√Hz — [Lee, ASEG 2001](https://www.tandfonline.com/doi/abs/10.1071/ASEG2001ab068)

### Inferences
- Timbuktu is the cleanest reproduction target, because the published input (lat/lon) and output (28.7068 m, EGM96) are explicit. Expect differences at the mm to cm level between direct synthesis and GeographicLib's interpolated 5′ grid.
- The NAVD 88 differences (−0.27 to −0.72 m) mostly reflect NAVD 88's own tilt and bias, not geoid error, so `geopotential` reproduces only the geoid side.
- Averaging check on the AQG: 500 nm/s² / √600 s = 20 nm/s², which matches the datasheet's "20 nm/s² in 10 min" (computed). The same arithmetic gives 2 µGal (20 nm/s²) at 625 s, which is the datasheet's own 10-minute figure (corrected 2026-10-24; an earlier version said 6,250 s).

### Gaps
- No verified global EGM2008 maximum (often said to be near New Guinea); it would have to be computed from the grid.
- I did not find the measured µGal amplitudes for Koymans et al. 2009–2015, the Etna AQG signal sizes, or the Troll/Alnes values.
- No satellite orbit-determination numbers (for example J2 nodal precession) were researched.
- No independent noise specification was found for Bell Geospace FTG.
- No AOSense datasheet was retrieved.
- I did not confirm whether NAPGD2022 had been formally released by October 2026.

## gravitational: how are gravitational-wave calculators used, and which published numbers can be reproduced?

### Takeaway
Pulsar timing gives the most exact target: the Hulse–Taylor orbital decay follows from published orbital elements through Peters–Mathews to 0.02% with rounded masses. Compact-binary parameters (GW150914, GW170817) and per-detector BNS ranges (O4: LIGO about 150–175 Mpc, Virgo about 50–60 Mpc) are the operational numbers. Cosmic Explorer, Einstein Telescope and LISA trade studies are what design calculators like the crate are used for.

### Cited Findings

1. **Binary pulsars (most direct: closed form)**
   - PSR B1913+16 (Weisberg & Huang 2016, ApJ 829:55), 9,257 TOAs over 35 years. Pb = 0.322997448918(3) d, e = 0.6171340(4), ⟨ω̇⟩ = 4.226585(4) deg/yr, γ = 0.004307(4) s, m_p = 1.438 ± 0.001 M☉, m_c = 1.390 ± 0.001 M☉. Observed Ṗb = −2.423(1) × 10⁻¹², intrinsic (after galactic correction) −2.398(4) × 10⁻¹², GR −2.40263(5) × 10⁻¹², ratio 0.9983 ± 0.0016 — [arXiv:1606.02744](https://arxiv.org/abs/1606.02744)
   - Double pulsar (Kramer et al. 2021), 16 years of timing: the orbit shrinks about 7 mm/day, in agreement with GR to 0.013% — [WVU GWAC](https://gwac.wvu.edu/blog/news-events/2021/12/13/einstein-wins-again); [APS Physics viewpoint](https://physics.aps.org/articles/pdf/10.1103/Physics.14.173)

2. **Compact-binary events (direct: chirp mass, frequency evolution, energy radiated, horizon distance)**
   - GW150914: 36 (+5/−4) and 29 ± 4 M☉, remnant 62 ± 4 M☉ with spin 0.67 (+0.05/−0.07), D_L 410 (+160/−180) Mpc, z 0.09, sky area 610 deg² — [arXiv:1602.03840](https://arxiv.org/pdf/1602.03840)
   - GW170817 (PRL 119:161101): network SNR 32.4, D_L 40 (+8/−14) Mpc, rest-frame chirp mass 1.188 (+0.004/−0.002) M☉, 28 deg² — [arXiv:1710.05832](https://arxiv.org/pdf/1710.05832.pdf)
   - GW250114 (2025) was the loudest event so far, at SNR about 77–80 (secondary sources). It gave the tightest test so far of Hawking's area theorem — [APS press release, Sept 2025](https://www.aps.org/about/news/2025/09/clearest-signal-hawkings-area-theorem); [arXiv:2507.08789](https://arxiv.org/pdf/2507.08789)

3. **Observing-run BNS range (direct if published ASDs are loaded; the range is an integral over the noise curve)**
   - O4 targets: LIGO 160–190 Mpc, Virgo 80–115 Mpc, KAGRA >1 Mpc at start. Achieved by April 2024: LIGO 155–175 Mpc, Virgo 55–60 Mpc. KAGRA aimed for about 10 Mpc after the earthquake recovery — [LVK observing plans](https://observing.docs.ligo.org/plan); [ICRR](https://gwcenter.icrr.u-tokyo.ac.jp/en/archives/1650).   - O4 ran 24 May 2023 to 18 Nov 2025 with about 250 real-time candidates. KAGRA counts about 300 significant signals. O5 is planned for 2027–2028, and a six-month intermediate run is expected from Nov 2026 — [LIGO](https://ligo.org/ligo-virgo-kagra-collaboration-successfully-wraps-up-its-fourth-observing-run/); [ICRR](https://gwcenter.icrr.u-tokyo.ac.jp/en/archives/2069)

4. **Pulsar-timing background (direct: h_c(f) = A (f/yr⁻¹)^(−2/3) to Ω_gw)**
   - NANOGrav 15-yr (June 2023): A = 2.4 (+0.7/−0.6) × 10⁻¹⁵ at f = 1/yr with γ = 13/3. The Hellings–Downs evidence is about 3σ (Bayesian, p = 10⁻³) and 3.5–4σ (frequentist, p = 5 × 10⁻⁵ to 1.9 × 10⁻⁴). The Bayes factor is 200–1000 — [arXiv:2306.16213](https://arxiv.org/pdf/2306.16213). A 2026 reanalysis gives A = 2.1 (+0.6/−0.5) × 10⁻¹⁵ at 3.16–3.32σ (from an aggregator; check against the primary) — [emergentmind listing](https://www.emergentmind.com/papers/2606.28554)

5. **Next-generation detector trade studies (indirect: needs design noise curves)**
   - Cosmic Explorer: BNS horizon z = 3.1 for CE Stage 1 (40 km) against 0.17 for LIGO A+ (BBH 26 against 2.6) — [Reitze et al. 2019, arXiv:1907.04833](https://arxiv.org/pdf/1907.04833). The reference design is one 40 km and one 20 km detector; the 20 km one keeps sensitivity near 3.7 kHz for post-merger signals — [CE design document](https://dcc.cosmicexplorer.org/CE-P2300001-v1/public). A 2025 network study gives BNS reach z ≈ 2.5 for CE 40 km at SNR 12 — [arXiv:2507.20258](https://arxiv.org/pdf/2507.20258)
   - Einstein Telescope: about €2 bn. The designs under debate are a 10 km triangle and a double-L with 15 km arms. Bid books are due end of 2026 and the site choice is in 2027. Lusatia and Sardinia agreed a joint double-L bid in July 2026 — [Physics World](https://physicsworld.com/a/physicists-set-to-decide-location-for-next-generation-einstein-telescope/); [TU Dresden, 17 Jul 2026](https://tu-dresden.de/kooperation/news/einstein-telescope-lausitz-und-sardinien-besiegeln-gemeinsame-bewerbung?set_language=en)
   - LISA: ESA adopted it on 25 Jan 2024, with launch about 2035 on Ariane 6 and 2.5 million km arms — [AEI](https://www.aei.mpg.de/1114814/esa-gives-go-ahead-for-flagship-gravitational-wave-observatory-in-space); [Wikipedia](https://en.wikipedia.org/wiki/Laser_Interferometer_Space_Antenna)
   - Atom-interferometer detectors: no atom-interferometer GW detector operates in 2026. MAGIS-100 at Fermilab (100 m shaft) expects commissioning by 2028 — [Fermilab](https://news.fnal.gov/tag/magis-100/). The AION-10 design report says it "is not expected to detect GWs" — [AION-10 TDR](https://arxiv.org/abs/2508.03491)

### Inferences
- Hulse–Taylor check (computed): Peters–Mathews with Pb, e, m_p = 1.438 and m_c = 1.390 gives Ṗb = −2.4021 × 10⁻¹², which is 0.02% from the published −2.40263 × 10⁻¹². The difference comes from the 3-digit masses, so this is an ideal doc-tested example.
- GW150914 radiated about 36 + 29 − 62 ≈ 3 M☉c² (computed from the cited masses).
- O4 BNS range is the best "decision" number. Range sets the expected detection rate, roughly ∝ range³, which is how observing-run trade-offs are argued.

### Gaps
- Peak strain 1.0 × 10⁻²¹ for GW150914 was not confirmed in these searches.
- No source was found on how GW calculators are used in education and outreach.
- The exact Cosmic Explorer Horizon Study numbers (arXiv:2109.09882) were not retrieved; the BNS reach varies from z ≈ 2 to 3 depending on assumptions.
- Final O4 per-detector median ranges were not retrieved. A search summary said O4b peaked at about 165 Mpc for LIGO and about 50 Mpc for Virgo, but I could not tie that to a specific LVK page.

## plasmoid: how are fusion and plasma calculators used, and which published design points can be reproduced?

### Takeaway
The cleanest reproducible fusion numbers are the peer-reviewed Wurzel & Hsu triple-product tables (n·T·τ_E multiplies back exactly) and energy-yield arithmetic (JET's 69 MJ needed 0.205 mg of burned D-T; NIF's 8.6 MJ is 3.05 × 10¹⁸ reactions). The Lawson criterion now appears in regulation: the NRC's 2026 fusion rule sits under the byproduct-material framework, not reactor licensing. Company design points (SPARC Q ≈ 11, Infinity Two 350 MWe, Helion 50 MWe) are **[claims]** that a 0-D power balance can check.

### Cited Findings

1. **Lawson and triple-product tables (most direct: arithmetic plus reactivity)**
   - Wurzel & Hsu 2025 (Phys. Plasmas 32:112106, v6 Nov 2025), "achieved" values: TFTR 83546 (1995): T_i 43 keV, n_i 6.6 × 10¹⁹ m⁻³, τ_E* 0.28 s, giving 7.9 × 10²⁰ keV m⁻³ s. JT-60U E31872 (1998): 16.8 keV, 4.8 × 10¹⁹, 0.69 s, giving 5.6 × 10²⁰. JET 42976 (1997): 28 keV, 0.51 s, giving 4.7 × 10²⁰ with Q_sci 0.63. JET 99972 (21 Dec 2021): 11 keV, 8.5 × 10¹⁹, 0.23 s, giving 2.2 × 10²⁰. JET 104600 (9 Oct 2023): 7 keV, 9 × 10¹⁹, 0.23 s, giving 1.4 × 10²⁰. W7-X (7 Dec 2017): 3.5 keV, 8 × 10¹⁹, 0.22 s, giving 6.2 × 10¹⁹. KSTAR 2014 and EAST 2015: 9.6 × 10¹⁸. MagLIF z3289 (2018): 2.6 keV, nτ 4.3 × 10²⁰, giving 1.1 × 10²¹. FuZE 2022: 2.4 keV, 4.9 × 10²³ m⁻³, 1 µs, giving 1.2 × 10¹⁸. C-2W 122588 (2023): about 1 keV, 1.69 × 10¹⁹, 2 ms, giving 3.4 × 10¹⁶. NIF N210808: 10.9 keV, nτ 1.4 × 10²¹, giving 1.6 × 10²² keV m⁻³ s — [arXiv:2505.03834v6](https://arxiv.org/html/2505.03834v6). The 2022 original notes that the triple product alone does not map to a unique Q, while (T, nτ) pairs do — [arXiv:2105.10954](https://arxiv.org/pdf/2105.10954)

2. **NIF ignition shots (direct: gain and reaction count)**
   - Wurzel & Hsu Table 4: N221204 (5 Dec 2022) 2.0 MJ in, 3.1 MJ out, Q_sci 1.51. N230729 2.0 MJ in, 3.9 MJ out, Q 1.89. N250406 (7 Apr 2025) 2.1 MJ in, 8.6 MJ out, Q 4.13. N250622 2.1 MJ in, 2.4 MJ out, Q 1.17 — [arXiv:2505.03834v6](https://arxiv.org/html/2505.03834v6)
   - LLNL: the 7 Apr 2025 shot delivered 2.08 MJ at 456 TW peak, yielded 8.6 MJ, and reached a record target gain of 4.13. The previous record was 5.2 MJ (Feb 2024) — [LLNL Annual Report FY2025](https://annual.llnl.gov/fy-2025/national-ignition-facility-2025). The yield uncertainty is ±0.45 MJ, and LLNL credits a better diamond capsule — [LLNL, Apr 2026](https://lasers.llnl.gov/news/target-breakthrough-enabled-fusion-record-nif). LLNL experts forecast that a gain approaching 15 would be needed for a zero-net-electric test facility (search summary) — [LLNL IFE-STAR](https://www.llnl.gov/article/52906/llnl-experts-foster-national-fusion-energy-ecosystem-ife-star-conference)

3. **JET final D-T record (direct: energy to fuel burned)**
   - 69.26 MJ over 5.2 s from 0.21 mg of fuel (DTE3, Aug–Oct 2023, announced Feb 2024). The previous record was 59 MJ (2021) — [DIFFER](https://www.differ.nl/news/JET-energy-record); [UK Gov](https://www.gov.uk/government/news/jets-final-tritium-experiments-yield-new-fusion-energy-record). Wurzel & Hsu list shot 104522 at P_in 35 MW, P_F 13 MW, Q_sci 0.37 — [arXiv:2505.03834v6](https://arxiv.org/html/2505.03834v6)

4. **Tokamak design points (0-D power balance) [claims and projections]**
   - SPARC (Creely et al., JPP 86(5), 2020): B₀ 12.2 T, R 1.85 m, a 0.57 m, goal Q > 2. With H98 = 1 it projects Q ≈ 11 and P_fus ≈ 140 MW; Q > 2 still holds at H98 = 0.7 — [JPP](https://resolve.cambridge.org/core/journals/journal-of-plasma-physics/article/overview-of-the-sparc-tokamak/DD3C44ECD26F5EACC554811764EF9FF0); [VTT](https://cris.vtt.fi/en/publications/overview-of-the-sparc-tokamak/). Physics-based modelling gave Q ≈ 9 — [APS DPP 2020](https://meetings-archive.aps.org/dpp/2020/jo08/4/)
   - SPARC status: the first of 18 TF magnets (about 24 t, 20 T) was installed in Jan 2026. First plasma and the Q > 1 target are both set for 2027 — [TechCrunch, 6 Jan 2026](https://techcrunch.com/2026/01/06/commonwealth-fusion-systems-installs-reactor-magnet-lands-deal-with-nvidia)

5. **Other company design points [claims]**
   - Type One Energy Infinity Two: D-T stellarator with 800 MW fusion power and a nominal 350 MWe, set out in 7 peer-reviewed JPP papers (March 2025), aimed at TVA — [Type One](https://typeoneenergy.com/insights-announcements/type-one-energy-first-realistic-unified-fusion-power-plant-design-basis/); [NEI Magazine](https://www.neimagazine.com/news/infinity-two-clears-first-major-design-review/)
   - Proxima Fusion Stellaris: quasi-isodynamic HTS stellarator in Fusion Eng. Des. (26 Feb 2025). The Alpha demo is planned for 2031 — [Proxima](https://proximafusion.com/press-news/proxima-fusion-and-partners-publish-stellaris-fusion-power-plant-concept-to-bring-limitless-safe-clean-energy-to-the-grid). IEEE Spectrum gives a 1 GW plant target — [IEEE Spectrum](https://spectrum.ieee.org/amp/stellarator-2671694687)
   - Helion: the Microsoft contract requires Orion to be grid-connected and operating by end-2028. In Sept 2026 a Helion official said the plant "will not, in no universe, have 50 megawatts hit full operation in 2028". The earlier claim was at least 50 MWe after ramp-up — [Axios, 8 Sep 2026](https://www.axios.com/pro/climate-deals/2026/09/08/helion-energy-fusion-net-electricity-date-year)
   - Zap Energy (FuZE, sheared-flow Z-pinch): Thomson scattering measured T_e = 1–3 keV on axis, peaking together with the neutron output (Levitt et al., PRL 132:155101, 2024) — [Zap Energy](https://www.zapenergy.com/updates/2024/04/elevated-electron-temp). The FuZE-Q diagnostics aim to infer about 1 keV — [APS DPP 2024](https://archive.aps.org/dpp/2024/tp12/109)

6. **Investor and regulatory context (not reproducible, but shows where the numbers are used)**
   - IAEA World Fusion Outlook 2025: more than 160 devices operating, under construction or planned; about 40 countries; private investment above USD 10 bn — [WNN](https://www.world-nuclear-news.org/articles/fusion-entering-implementation-phase-grossi). The FIA 2025 survey counts USD 9.766 bn across 53 companies, and F4E counts €13 bn across 77 companies (search summary) — [ANS investment tag](https://www.ans.org/news/tag-investment)
   - The US NRC proposed rule "Regulatory Framework for Fusion Machines" (Federal Register, 26 Feb 2026) amends 10 CFR Part 30 (byproduct material) under ADVANCE Act §205. It defines fusion machines as particle accelerators, with "low offsite dose" meaning <1 rem (10 mSv), and comes with draft NUREG-1556 Vol. 22 — [Federal Register](https://www.govinfo.gov/content/pkg/FR-2026-02-26/html/2026-03865.htm); [Orrick](https://www.orrick.com/en/Insights/2026/03/NRC-Proposed-Fusion-Rule-Further-Clarifies-Path-for-Commercial-Deployment)
   - Industrial plasma, briefly: SpaceX's Starlink V2-mini argon Hall thruster gives 2.4× the thrust and 1.5× the specific impulse of the earlier krypton thruster (SpaceX figures via TechCrunch, Mar 2023) — [TechCrunch](https://techcrunch.com/2023/03/03/spacexs-acquisition-of-swarm-is-paying-off-with-new-starlink-thrusters)

### Inferences
- The triple products recompute exactly from the table (computed): JET 99972 gives 2.15 × 10²⁰, TFTR 7.95 × 10²⁰, JT-60U 5.56 × 10²⁰, FuZE 1.18 × 10¹⁸, C-2W 3.38 × 10¹⁶ and W7-X 6.16 × 10¹⁹ keV m⁻³ s. That makes them ideal doc-tests, with the crate's ⟨σv⟩ supplying the ignition curve.
- JET check (computed): 69.26 MJ ÷ 17.589 MeV = 2.46 × 10¹⁹ D-T reactions, times 5.030 u = 0.205 mg. That matches the "0.21 mg of fuel", so the published fuel mass is the mass burned.
- NIF check (computed): 8.6 MJ ÷ 17.589 MeV = 3.05 × 10¹⁸ reactions (neutrons). For N221204, 3.15 MJ gives 1.12 × 10¹⁸.
- Infinity Two's 350 MWe from 800 MW fusion is a 44% fusion-to-net-electric ratio, which a power-balance function can sanity-check **[claim-derived]**.

### Gaps
- I did not find an investor document that shows a triple-product or Lawson calculation explicitly. The regulatory documents use dose criteria, not Lawson.
- TAE (Norm, Copernicus) and Helion Polaris have no peer-reviewed 2025–2026 performance numbers. Polaris D-T claims appeared only in secondary sources.
- The Lawson ignition threshold as a single number (about 3 × 10²¹ keV m⁻³ s near 14 keV) was not confirmed from the fetched text.
- Plasma-etching numbers were not researched.

## magnetics (planned): how are geomagnetic models, MHD, magnets and space weather used, and which numbers can be reproduced?

### Takeaway
The WMM2025 test-value table is a ready-made doc-test: 12 field rows plus secular-variation rows at given date, height, latitude and longitude. The WMM is the navigation standard for DoD, UK MoD, NATO and IHO, and it sits inside phone compasses. Directional drilling needs declination to 0.1°. Grid operators plan for an 8 V/km 1-in-100-year geoelectric field, and the May 2024 storm drove up to 113 A in New Zealand transformers. The existing Rust geomagnetic crates are outdated or alpha.

### Cited Findings

1. **WMM/IGRF field synthesis (most direct: official test values)**
   - WMM2025 test values **[NCEI, Jan 2025]**, computed in double precision (single precision can differ by up to 0.1 nT), heights above the WGS84 ellipsoid. At 2025.0, h = 0 km, 80°N 0°E: X 6521.6, Y 145.9, Z 54791.5, H 6523.2, F 55178.5 nT, I 83.21°, D 1.28°. At 0°N 120°E: X 39677.8, Y −109.6, Z −10580.2, F 41064.3 nT, I −14.93°, D −0.16°. At 80°S 240°E: X 6117.5, Y 15751.9, Z −52022.5, F 54698.2 nT, I −72.00°, D 68.78°, GV −51.22°. At h = 100 km, 80°N 0°E: F 52964.9 nT, D 0.85°. At 2027.5, h = 0, 80°N 0°E: F 55253.9 nT, D 2.59°. Secular variation at 2025.0, 80°N 0°E: Ẋ −8.3, Ẏ 59.5, Ż 31.1 nT/yr, Ḋ 0.52°/yr — [NCEI WMM2025 test values PDF](https://www.ncei.noaa.gov/sites/default/files/2025-02/WMM2025testvalues.pdf)
   - WMM2025 was released 17 Dec 2024 by NGA, NOAA, BGS and the UK DGC, covering 2025–2030. It comes as standard (3,300 km resolution at the equator) and WMMHR2025 (about 300 km). It is the standard for navigation and attitude/heading systems for the US DoD, UK MoD, UK Hydrographic Office, NATO and IHO — [NGA](https://www.nga.mil/news/NGA_NOAA_BGS_Publish_World_Magnetic_Model_2025.html). NCEI says its accuracy affects phone and GPS-device performance — [NCEI](https://www.ncei.noaa.gov/news/world-magnetic-model-2025-released); [NESDIS technical report release, Aug 2025](https://www.nesdis.noaa.gov/news/noaabritish-geological-survey-release-2025-world-magnetic-model-report)
   - The north magnetic pole moved about 36 km/yr over the past year and the south pole about 9 km/yr, down from 50–60 km/yr in the early 2000s — [NCEI](https://www.ncei.noaa.gov/node/8576)
   - Existing Rust crates: `wmm` 0.2.3 (no_std, declination only, model valid 2020–2025, last release Nov 2022) — [docs.rs](https://docs.rs/crate/wmm/latest); `igrf` 0.2.0 ("not well tested", possible index panics) — [docs.rs](https://docs.rs/crate/igrf/latest); `ferromagnetic` 0.0.1-alpha.5 — [docs.rs](https://docs.rs/ferromagnetic)

2. **Directional drilling (direct: declination and dip accuracy to wellbore position)**
   - BGS: MWD magnetic surveys need real-time field estimates. The desired accuracy is 0.1° in declination, 0.05° in dip and 50 nT in total field. BGGM is re-released every year and covers the core plus large-scale crust. IFR1 adds local crustal anomalies, which can be "several tenths of a degree" and cause systematic errors that grow with depth. IIFR uses nearby observatory data to remove the external field — [BGS IFR](https://geomag.bgs.ac.uk/data_service/directionaldrilling/ifr.html); [BGS directional drilling](https://geomag.bgs.ac.uk/data_service/directionaldrilling/home.html); [BGS IFR 2025 brochure](https://esc.bgs.ac.uk/documents/Geomag_Services_IFR1_IFR2_2025-OnLineVersion.pdf)

3. **Space weather and grid GIC (indirect: needs a ground-conductivity model; the benchmark scaling is direct)**
   - NERC TPL-007 benchmark: a peak geoelectric field of 8 V/km at 60° geomagnetic latitude over resistive ground, a 1-in-100-year event based on the 1989 storm. It scales as E_peak = 8 × α × β V/km. A supplemental 12 V/km, not spatially averaged, was added at FERC's direction — [NERC FAQ](https://www.nerc.com/pa/Stand/Project201303GeomagneticDisturbanceMitigation/common_questions_June12.pdf). Transformer thermal assessments apply above 75 A/phase (benchmark) and 85 A/phase (supplemental) (2020 slides) — [ERCOT](https://www.ercot.com/files/docs/2020/05/15/GMD_Intro_and_Preliminary_Benchmark_GMDVA_Voltage_Findings.pdf). An advocacy critique says spatial averaging gives only 2.4 V/km for Washington DC — [STG Coalition comments to DOE](https://www.energy.gov/sites/default/files/2022-07/STG-Coalition-Comments-SEAB-13June2022-Final.pdf)
   - March 1989: the Hydro-Québec grid collapsed at 2:44 a.m. on 13 March and was dark for more than 9 hours — [Hydro-Québec](https://hydroquebec.com/learning/notions-de-base/tempete-mars-1989.html). About 6 million people lost power — [Wikipedia](https://en.wikipedia.org/wiki/March_1989_geomagnetic_storm)
   - Sweden: on 30 Oct 2003 Malmö lost power for about 50,000 customers for 20–50 min — [Wik et al. 2009, Ann. Geophys.](https://angeo.copernicus.org/articles/27/1775/2009/angeo-27-1775-2009.pdf). The largest Swedish GIC was almost 300 A in a 400 kV transformer neutral on 6 Apr 2000 — [Lund University](https://portal.research.lu.se/en/publications/calculation-of-geomagnetically-induced-currents-in-the-400-kv-pow/). Another abstract says 270 A — [EGU2007](https://www.meetings.copernicus.org/www.cosis.net/abstracts/EGU2007/03121/EGU2007-J-03121.pdf)
   - May 2024 "Gannon" storm: provisional Dst minimum −412 nT around 02:00 UT on 11 May — [flash data report, arXiv:2407.07665](https://arxiv.org/pdf/2407.07665). At New Zealand's Halfway Bush substation, 220 kV transformers saw GIC above 50 A, peaking at 113 A, with about 200 A modelled without mitigation — [Mac Manus et al. 2025, Space Weather](https://nora.nerc.ac.uk/id/eprint/539520). In the US, GIC was measured at 47 sites, and TVA's modelled values correlated at r > 0.8 — [arXiv:2507.07009](https://arxiv.org/html/2507.07009v4). NERC found the bulk power system stable, with equipment impacts confined to the US Northeast, Mid-Atlantic and Western Canada — [NERC after-action plan](https://nerc.com/pa/RAPA/GMD/Documents/Event_After_Action_Plan_May_2024.pdf)

4. **MRI field strengths (direct: B²/2μ₀, Larmor frequency, fringe-field falloff)**
   - FDA 2003 guidance: static fields up to 8 T are non-significant risk for adults, children and infants over 1 month. Studies above 8 T need FDA (IDE) and IRB approval — [FDA 2003 guidance](https://users.fmrib.ox.ac.uk/~peterj/safety_docs/fda_guidelines_2003.pdf); [PMC8170917](https://pmc.ncbi.nlm.nih.gov/articles/PMC8170917)
   - Iseult (CEA NeuroSpin) at 11.7 T produced its first human brain images in April 2024, from about 20 volunteers. Hospital MRIs rarely exceed 3 T — [AuntMinnie Europe](https://www.auntminnieeurope.com/clinical-news/mri/article/15667839/frenchled-team-unveils-117-tesla-mri-scans-of-human-brain); [AFP via Malay Mail](https://www.malaymail.com/news/life/2024/04/03/worlds-most-powerful-mri-scans-first-images-of-human-brain/127009)

5. **Maglev (weak link: speeds only)**
   - The Shanghai Transrapid's regular-service record is 431 km/h — [La Sicilia](https://www.lasicilia.it/news/attualita/3014968/giappone-lora-della-supervelocita-il-treno-che-tocca-603-5-km-h-e-lenigma-europeo.html). The JR Central L0 reached 603 km/h in 2015 and levitates about 10 cm above about 150 km/h. The Chuo Shinkansen is planned for 500 km/h, with opening dates of 2034–2037 depending on the source — [AutoPapo](https://autopapo.com.br/curta/trem-mais-rapido-do-mundo-faz-603-km-h-mas-esta-prestes-a-ser-ultrapassado/); [Bloomberg Línea](https://www.bloomberglinea.com/estilo-de-vida/el-fin-de-una-era-china-prepara-el-tren-bala-mas-rapido-del-mundo-que-desafia-a-japon/)

### Inferences
- The test table is self-consistent (computed): atan2(145.9, 6521.6) = 1.28° and atan2(54791.5, 6523.2) = 83.21°. A planned `magnetics` crate can doc-test all 12 rows plus 12 SV rows to 0.1 nT.
- Drilling consequence (computed): a 0.1° azimuth error over a 5 km horizontal step-out moves the bottom hole about 8.7 m sideways. This is why 0.1° is the target and why crustal anomalies of tenths of a degree matter.
- MRI (computed): at 11.7 T the magnetic pressure B²/2μ₀ is about 54.5 MPa, and the proton Larmor frequency is about 498 MHz (using γ/2π = 42.577 MHz/T, a standard constant).
- Because `wmm` stops at 2025 and `igrf` warns of panics, a tested WMM2025/IGRF-14 Rust crate fills a real gap.

### Gaps
- No IGRF-14 test values were fetched (they would mirror the WMM table).
- I did not get a declination for a named city from an official source, only the test grid points.
- No ISCWSA MWD error-model numbers were found.
- No measured Hydro-Québec dB/dt, and no 2003 Swedish GIC amperes.
- NERC's final May 2024 report, with US transformer amperes, was not located.
- The polar "blackout zone" threshold (2000 nT horizontal field) was not confirmed.
- No MHD application numbers (for example liquid-metal pumps or solar wind) were researched.

## gravitics: how have frontier-physics claims been computed against known physics, and which numbers support research planning?

### Takeaway
Agencies have used "compute the claim against the limit": NASA's Breakthrough Propulsion Physics project (1996–2002), the NASA Marshall Podkletnov replication, the Dresden EmDrive tests and the JASON hafnium review. Each turned an extraordinary claim into an effect size and compared it with a measured bound. The planning-grade numbers are the Schwinger field (1.32 × 10¹⁸ V/m) against the 1.1 × 10²³ W/cm² laser record, the Eöt-Wash 38.6 µm Yukawa bound, and G's 2.2 × 10⁻⁵ scatter.

### Cited Findings

1. **Scatter of G (direct: weighted mean and Birge ratio)**
   - CODATA 2022: σ_G/G = 2.2 × 10⁻⁵, from 16 G measurements with no new inputs since 2018. A large expansion factor of 3.9 is applied to the input uncertainties — [CODATA 2022, arXiv:2409.03787](https://arxiv.org/pdf/2409.03787); [space-based G paper](https://arxiv.org/pdf/2609.09485). CODATA G = 6.67430(15) × 10⁻¹¹ m³ kg⁻¹ s⁻², as quoted in [Westphal et al.](https://arxiv.org/abs/2009.09546)

2. **Short-range gravity, Eöt-Wash style (direct: Yukawa α–λ exclusion from a torque budget)**
   - Lee, Adelberger, Cook, Fleischer & Heckel (2020): a torsion balance over a rotating attractor with 18-fold and 120-fold symmetries, at separations of 52 µm to 3.0 mm. Gravitational-strength Yukawa interactions are excluded for λ > 38.6 µm at 95% confidence — [arXiv:2002.11761](https://arxiv.org/pdf/2002.11761). The thesis quotes 39 µm and found a new systematic from test-mass magnetic susceptibility — [UW thesis](https://digital.lib.washington.edu/researchworks/handle/1773/45542)
   - Tabletop quantum gravity: the smallest source mass with resolved gravity is a 92.1 mg gold sphere (Westphal et al., *Nature* 2021) — [arXiv:2009.09546](https://arxiv.org/abs/2009.09546)

3. **Schwinger limit and laser-facility planning (direct)**
   - E_S = m²c³/(eħ) = 1.3233 × 10¹⁸ V/m and B_S = 4.414 × 10⁹ T. I_S is 4.65 × 10²⁹ or 2.32 × 10²⁹ W/cm² depending on convention; Di Piazza et al. use I_cr ≈ 4.6 × 10²⁹ W/cm² — [Di Piazza et al., RMP 84:1177](https://arxiv.org/abs/1111.3886)
   - Laser record: (1.1 ± 0.2) × 10²³ W/cm² from the 4 PW CoReLS laser (2021, *Optica*), focused f/1.1 to a spot about 1.1 µm across — [Guinness](https://www.guinnessworldrecords.com/world-records/692568-highest-intensity-laser-beam); [IBS](https://pr.ibs.re.kr/handle/8788114/10375)
   - SLAC E-144 (1997): 46.6 GeV electrons hit a 527 nm laser at about 1.3 × 10¹⁸ W/cm², producing 106 ± 14 positrons above background by multiphoton Breit–Wheeler (n ≈ 5 photons) — [SLAC-PUB-7444](https://www.slac.stanford.edu/pubs/slacpubs/7250/slac-pub-7444.pdf); [PRL 79:1626](https://www.lle.rochester.edu/media/timeline/documents/1997/PhysRevLett.79.1626.pdf)

4. **Claim-versus-bound evaluations by agencies and reviewers (direct: effect size against null bound)**
   - NASA Breakthrough Propulsion Physics (1996–2002): about $1.2 M total per Wikipedia. Of 14 topics, 6 were found nonviable, 4 merited more work and 4 were unresolved. It ended in 2002 when work below TRL 3 was cut — [Wikipedia](https://en.wikipedia.org/wiki/Breakthrough_Propulsion_Physics_Project). A June 2002 NASA Glenn notice lists about $1.7 M of in-house work over 7 years (8 tasks) and about $1.4 M for an MSFC gravity-shield test, figures that conflict with the Wikipedia total — [SpaceRef](https://spaceref.com/?p=92981). Task selection multiplied scoring criteria rather than adding them — [Millis 2004, NASA NTRS](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20040070788.pdf)
   - Podkletnov claimed a 0.05–0.3% weight loss over a rotating superconductor (1992). Hathaway, Cleveland & Bao (*Physica C* 2003), close to the NASA Marshall effort, saw no effect down to the 0.001% level — [review, arXiv:2206.07574](https://arxiv.org/pdf/2206.07574)
   - EmDrive: Tajmar's TU Dresden group (three CEAS Space Journal papers, 2021) traced the apparent thrust to thermal effects. They said their measurements "refute all EmDrive claims by at least 3 orders of magnitude" — [Wikipedia](https://en.wikipedia.org/wiki/EmDrive); [Big Think](https://bigthink.com/the-future/emdrive-debunked)
   - Hafnium-178m2 isomer "bomb": JASON called Collins's 1999 PRL result "a priori implausible" and "inconclusive, at best", and IDA called it flawed — [Fourmilog summary](https://www.fourmilab.ch/fourmilog/archives/2006-06/000715.html). Argonne found no triggering, and a 2008 DARPA/LLNL TRIP review found the results "statistically marginal and inconsistent" — [Physics Today](https://physicstoday.aip.org/news/conflicting-results-on-a-long-lived-nuclear-isomer-of-hafnium-have-wider-implications)

5. **MOND and wide binaries (indirect: needs Gaia data)**
   - Banik et al. 2024 (MNRAS 527:4573), 8,611 wide binaries: Newtonian dynamics preferred at 19σ, with MOND excluded at 16σ — [arXiv:2311.03436](https://arxiv.org/abs/2311.03436). The dispute with Chae and others was still open in the sources found.

### Inferences
- Laser gap (computed): √(1.1 × 10²³ / 4.65 × 10²⁹) ≈ 4.9 × 10⁻⁴, or √(1.1 × 10²³ / 2.32 × 10²⁹) ≈ 6.9 × 10⁻⁴. So the best lasers reach about 0.05–0.07% of E_S in the lab frame. Pair-production experiments (E-144, LUXE) therefore rely on the electron-frame boost, not on lasers alone, which is the key planning number.
- Podkletnov and Hathaway compare a 0.05–0.3% claim with a 0.001% null, a factor of 50–300. Tajmar's "3 orders of magnitude" is the same kind of number. A `gravitics` "claim ratio" function could make these one-line doc-tests.
- G's 2.2 × 10⁻⁵ uncertainty after a ×3.9 expansion suggests a Birge-ratio example built on the 16 CODATA inputs.

### Gaps
- No JASON report text was found that covers gravity-control or propulsion claims; only the hafnium review was located. The DIRDs material is in the repo's fringe-physics notes, not re-checked here.
- The current status of ELI-NP's 10 PW operation was not confirmed.
- The LUXE start date was not researched.
- No source shows a laser facility citing Schwinger arithmetic in a procurement or planning document.
- The MOND a₀ value was not re-sourced here.

## Why Rust: what published uses of Rust in space, geodesy or scientific software support Rust crates over Python?

### Takeaway
Rust has flight heritage: Nyx/ANISE/Hifitime were the flight-dynamics backbone for Firefly's Blue Ghost 1 Moon landing on 2 Mar 2025, and ANISE is now described as TRL 9. There is a Rust onboard-software framework (sat-rs) for Stuttgart's ROMEO satellite, and the Ferrocene compiler is safety-qualified (ISO 26262 and IEC 61508, 2023). The geomagnetic niche has only stale or alpha crates.

### Cited Findings
- Nyx, ANISE and Hifitime (open-source Rust) were "the backbone" of the flight-dynamics toolset for Firefly's Blue Ghost Mission 1, which landed in Mare Crisium on 2 Mar 2025 — [Nyx Space blog](https://nyxspace.com/blog/category/blue-ghost/). ANISE rewrites NAIF SPICE's core functions. It was used from launch to landing and is "now considered NASA TRL 9". The vendor **[claim]** says it predicted the first eclipse to within 1 s. Nyx is AGPL-3.0 and ANISE/Hifitime are MPL-2.0, with Python bindings via PyO3 — [ANISE GitHub](https://github.com/nyx-space/anise)
- sat-rs is a Rust flight-software framework from the University of Stuttgart IRS, with an STM32 bare-metal example and a `spacepackets` crate for ECSS/CCSDS packets. The ROMEO satellite's onboard software is written in Rust, and sat-rs is meant to replace the institute's C++ framework — [sat-rs README](https://egit.irs.uni-stuttgart.de/rust/sat-rs/src/commit/1061ab52870cfe065d35f5e2e746824c1a76576a/README.md)
- Ferrocene: TÜV SÜD qualified it under ISO 26262 and IEC 61508 (Oct 2023, version 23.06 on Rust 1.68). Sources disagree on SIL 3 vs SIL 4 — [Ferrous Systems](https://ferrous-systems.com/blog/officially-qualified-ferrocene/). A core-library subset was certified to IEC 61508 SIL 2 in Dec 2025 — [Business Wire](https://www.businesswire.com/news/home/20251203818102/en), and ISO 26262 ASIL B followed in 26.02.0 — [Ferrous Systems](https://ferrous-systems.com/blog/ferrocene-26-02-0/)
- Perkel, "Why scientists are turning to Rust", *Nature* 588:185–186 (Dec 2020), exists, but its content was paywalled and not read — [doi:10.1038/d41586-020-03382-2](https://doi.org/10.1038/d41586-020-03382-2)
- Existing Rust geomagnetic crates are limited: `wmm` is WMM2020-only, `igrf` is "not well tested" and `ferromagnetic` is alpha — [docs.rs wmm](https://docs.rs/crate/wmm/latest); [docs.rs igrf](https://docs.rs/crate/igrf/latest); [docs.rs ferromagnetic](https://docs.rs/ferromagnetic)

### Inferences
- The strongest argument is reuse in flight and embedded contexts (no_std, certified toolchain) together with PyO3 bindings for Python users, which is the pattern ANISE followed to TRL 9.
- These crates' test-value doc-tests (WMM table, GeoidEval Timbuktu, Hulse–Taylor Ṗb, Wurzel–Hsu triple products) give the kind of verification evidence that a TRL or safety-case argument needs.

### Gaps
- I found no published case of a WebAssembly scientific tool in geodesy or GW written in Rust.
- No peer-reviewed benchmark of Rust against Python for geodesy or gravity-field synthesis was found.
- ROMEO's launch and on-orbit status were not confirmed.
