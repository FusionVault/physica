# Pair open GRACE data with simulated quantum sensors

No satellite has released quantum gravity data yet, and none will before about 2030. Every classical mission named in the Open Access Government article does have open data that a program can fetch today. GRACE (2002–2017), GRACE-FO (2018 to now, with August 2026 the latest monthly field at NASA PO.DAAC), GOCE, CHAMP, Swarm and satellite laser ranging all publish spherical-harmonic or gridded gravity products. Most can be downloaded anonymously over HTTPS from GFZ's ICGEM, TU Graz, CSR, GSFC and ESA's Swarm server; the rest need a free NASA Earthdata or ESA login. The quantum side so far consists of hardware pathfinders (Europe's CARIOQA and NASA's QGGPf, both aiming at 2030), a 75 ms atom interferometer on the China Space Station, and simulation studies. A useful Rust application must therefore do two jobs with one engine. It must turn real GRACE-era coefficients into geoid, gravity and water-storage maps, and it must push published quantum-sensor noise models through to the same error spectra and maps. Of the three linked articles, only Zingerle et al. 2024 has data figures. Its instrument-noise figure can be rebuilt exactly from equations printed in the paper. Its map and static-field figures can be rebuilt in kind from public models plus a least-squares error propagator. Its aliasing figures need a closed-loop orbit simulator that the authors did not release. The phys.org story is a NASA Goddard/AOSense press release, not Hannover work, and neither it nor the opinion piece has anything to plot. geo5d v0.2.0 already has a working point-wise Holmes–Featherstone engine and an ICGEM parser. It lacks time-variable fields, grids, filters and water conversion. The strongest plan keeps geo5d as the geodesy façade and adds four companion crates, `geopotential`, `gravity-field`, `gravimetry` and `satgrav`, all unclaimed on crates.io as of 22 October 2026, plus an unpublished egui/plotters viewer. The work is staged in five phases, each closed by tests that check the code against independent references.

## All open gravity data is classical, and GRACE-FO still delivers monthly fields

The GRACE-FO data stream is healthy in October 2026, so the classical record that any quantum mission will be judged against is still growing. PO.DAAC's operational GRACE-FO release remains **RL06.3**, with open-ended monthly collections from CSR, GFZ and JPL. The latest month published on 22 October 2026 was **August 2026** (`GSM-2_2026213-2026243_GRFO_UTCSR_BA01_0603`, posted 1 October). Level-1B daily bundles reached 24 September 2026 ([CMR granules](https://cmr.earthdata.nasa.gov/search/granules.json?short_name=GRACEFO_L2_CSR_MONTHLY_0063&provider=POCLOUD&sort_key=-start_date&page_size=2); [CMR L1B](https://cmr.earthdata.nasa.gov/search/granules.json?short_name=GRACEFO_L1B_ASCII_GRAV_JPL_RL04&provider=POCLOUD&sort_key=-start_date&page_size=2)). Reprocessing is arriving unevenly. **GRACE RL07** collections from CSR and JPL appeared at PO.DAAC in 2026 for the 2002–2017 period ([CMR UMM](https://cmr.earthdata.nasa.gov/search/collections.umm_json?short_name=GRACE_L2_CSR_MONTHLY_0070&provider=POCLOUD)). GRACE-FO RL07 was planned for "Spring '26" but is not yet in the catalogue ([Landerer, GRACE-FO STM 2025](https://cpaess.ucar.edu/sites/default/files/grace-fo-2025-presentation/landerer-felix.pdf)). That is a reason for software to discover releases through the catalogue rather than hard-code them. The mission's longer-term status is partly unresolved. JPL called GRACE-FO "healthy" in its final extended year through September 2026, and GFZ funds operations through **December 2029** to overlap with GRACE-C ([Landerer 2025](https://cpaess.ucar.edu/sites/default/files/grace-fo-2025-presentation/landerer-felix.pdf)). No published outcome of NASA's 2026 Senior Review was found, although data were still flowing in October.

| Mission (agency) | Span | Open products | Best programmatic source | Login |
|---|---|---|---|---|
| GRACE (NASA/DLR) | 2002-04 → 2017-06 | L1B; monthly SH RL06 and new RL07; L3 grids; mascons | PO.DAAC; ICGEM; TU Graz ITSG-Grace2018 | Earthdata for PO.DAAC payloads; others anonymous |
| GRACE-FO (NASA/GFZ) | 2018-05 → Aug 2026 | L1B daily; monthly SH RL06.3 (CSR/JPL/GFZ); 1° grids; JPL/CSR/GSFC mascons | PO.DAAC; CSR and GSFC mascon sites; ICGEM COST-G (to Jul 2026); TU Graz ITSG operational (to Aug 2026) | Mixed |
| GOCE (ESA) | 2009–2013 | Gradients `EGG_NOM_2_`, `EGG_TRF_2_`; orbits `SST_PSO_2_`; models | ESA goce-ds; models also on ICGEM as `GO_CONS_GCF_2_*` | ESA EO Sign In for native files |
| Swarm (ESA) | 2013-12 → 2025-12 | Monthly GPS-derived SH to degree 40, released quarterly | swarm-diss JSON API; ICGEM COST-G Swarm | Anonymous |
| CHAMP (GFZ) | 2000–2010 | CHAMP2013s monthly series | ICGEM | Anonymous |
| SLR (LAGEOS, LARES …) | 2002 → Aug 2026 | TN-14 C20/C30; CSR 5×5 monthly | PO.DAAC docs bucket; CSR filedrop | Anonymous |
| Static combinations | — | 184 models, e.g. EGM2008, XGM2019e, GOCO06s, GOCO2025s | ICGEM | Anonymous |
| Chinese GRACE-like / quantum | — | No public release found | — | — |

Sources for the table: [CMR collection search](https://cmr.earthdata.nasa.gov/search/collections.json?keyword=GRACE&provider=POCLOUD&page_size=200), [GOCE dissemination](https://goce-ds.eo.esa.int/), [Swarm listing](https://swarm-diss.eo.esa.int/?do=list&maxfiles=500&pos=0&file=swarm%2FLevel2longterm%2FEGF), [ICGEM temporal list](https://icgem.gfz.de/sl/temporal), [ICGEM static table](https://icgem.gfz.de/tom_longtime), [TN-14](https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-docs/gracefo/open/docs/TN-14_C30_C20_GSFC_SLR.txt).

Swarm deserves attention because it is the only gravity source spanning the 2017–2018 gap between GRACE and GRACE-FO. Its **145 monthly fields** are anonymous downloads of about 19 kB each ([swarm-diss](https://swarm-diss.eo.esa.int/?do=list&maxfiles=500&pos=0&file=swarm%2FLevel2longterm%2FEGF)). They agree with GRACE to about **4 cm equivalent water height (EWH)** over land at 750 km smoothing ([Teixeira da Encarnação et al. 2020](https://essd.copernicus.org/articles/12/1385/2020/)). On the quantum side, **no in-orbit cold-atom gravimetry data has been released**. China's space-station results are accelerometer and gyroscope demonstrations, not gravity-field products ([Li et al. 2025, NSR](https://academic.oup.com/nsr/article/12/4/nwaf012/7951977)). No portal exists for China's "Gravimetry augment and Mass change exploring" missions either.

### Three access tiers decide how the fetcher must authenticate

Programmatic access splits into three tiers, and a Rust client needs a different code path for each. **Discovery needs no account.** NASA's CMR REST API returns collection and granule JSON over plain HTTPS. The total hit count arrives in the `CMR-Hits` header, and each granule lists HTTPS, S3 and OPeNDAP links ([live CMR query](https://cmr.earthdata.nasa.gov/search/granules.json?short_name=GRACEFO_L2_CSR_MONTHLY_0063&provider=POCLOUD&sort_key=-start_date&page_size=2)). **Anonymous downloads** cover a lot of ground:

- ICGEM, for every `.gfc` static model and temporal series. Static models sit at `https://icgem.gfz.de/getmodel/gfc/<sha256>/<MODEL>.gfc`. The hash is opaque, so it must be scraped from the `/tom_longtime` table. Series follow `/getseries/<group>/<series>/[subdir/]<file>.gfc` ([ICGEM static table](https://icgem.gfz.de/tom_longtime); [ICGEM temporal](https://icgem.gfz.de/sl/temporal)).
- TU Graz's ITSG directory listings ([TU Graz](https://ftp.tugraz.at/pub/ITSG/GRACE/ITSG-Grace_operational/monthly/)).
- CSR's single-file mascon NetCDF, about 113.5 MB ([CSR mascons](https://www2.csr.utexas.edu/grace/RL06_mascons.html)).
- GSFC mascons ([GSFC](https://earth.gsfc.nasa.gov/geo/data/grace-mascons)).
- Swarm's JSON listing plus `?do=download&file=…` endpoint.
- PO.DAAC's docs bucket and public sidecars, which hold TN-13, TN-14, release notes and the `.txt`/`.tif` siblings of the Level-3 grids.

**A free NASA Earthdata Login** is required for every PO.DAAC data payload. An unauthenticated GET returns HTTP 302 to `urs.earthdata.nasa.gov/oauth/authorize`, so the client must follow a cookie-preserving redirect with basic auth on the URS hop, or send a bearer token. Raw bearer-token behaviour on archive URLs was not tested. The Python `earthaccess` library reads credentials from environment variables, then `~/.netrc`, then a prompt, which is a sensible precedent for a Rust client ([earthaccess docs](https://earthaccess.readthedocs.io/en/latest/user/howto/authenticate/)). Direct S3 access works only from AWS us-west-2 ([s3credentials README](https://archive.podaac.earthdata.nasa.gov/s3credentialsREADME)). **ESA EO Sign In** is needed only for native GOCE files ([goce-ds](https://goce-ds.eo.esa.int/)). GFZ's ISDC moved to HTTPS-only in June 2025 and switched off FTP at the end of 2025 ([Dahle 2025](https://cpaess.ucar.edu/sites/default/files/grace-fo-2025-presentation/dahle-christoph-sds.pdf)). It refused connections during this research, so its directory paths are unverified.

The lowest-friction default is therefore ICGEM, TU Graz, CSR and Swarm with no credentials at all. Earthdata Login becomes an opt-in for PO.DAAC Level-1/2/3 files. Mascons are the easiest Level-3 input because they need no filtering. The JPL file `GRCTellus.JPL.200204_202607.GLO.RL06.3M.MSCNv04CRI.nc` holds the whole record on a 0.5° grid (time 258 × lat 360 × lon 720). It stores `lwe_thickness` in cm, with time in days since 2002-01-01 and the 2004.000–2009.999 mean removed ([PO.DAAC virtual-zarr refs](https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-public/virtual_collections/TELLUS_GRAC-GRFO_MASCON_CRI_GRID_RL06.3_V4/TELLUS_GRAC-GRFO_MASCON_CRI_GRID_RL06.3_V4_virtual_https.json)).

### Parsers must tolerate real files, not just the specification

Three text formats and one NetCDF family cover almost everything.

**ICGEM `.gfc`** (spec dated February 2023). The header ends at `end_of_head`, with mandatory `earth_gravity_constant`, `radius`, `max_degree` and `errors`, and optional `tide_system` and `norm`. Data keys are `gfc`, `gfct`, `trnd`/`dot`, `asin` and `acos`. In `icgem2.0`, every time-variable record carries a `t0 t1` validity interval in `yyyymmdd.hhmm`. The time evaluation is $G(t)=G(t_0)+\mathrm{trnd}\,(t-t_0)+\sum_k[\mathrm{asin}_k\sin(2\pi(t-t_0)/p_k)+\mathrm{acos}_k\cos(2\pi(t-t_0)/p_k)]$, with $t$ in years ([ICGEM format PDF](https://icgem.gfz.de/files/6defd9f87a64b16ec96bd657a74d0292ae3447cc7418da5104441102d8cb2b58.pdf)). Real files break the spec:

- COST-G monthly files have no `begin_of_head`, use `CMMNT` lines, use the non-standard `errors empirical`, and write numbers like `+0.10000000000E+01` ([COST-G July 2026](https://icgem.gfz.de/getseries/02_COST-G_/Grace-Grace-FO_RL02.1/GSM-2_2026182-2026212_GRFO_COSTG_BF01_0201.gfc)).
- Swarm ZIPs carry a stray `WARNING:` line in the header.
- GOCO06s reorders its header keys ([GOCO06s.gfc](https://icgem.gfz.de/getmodel/gfc/32ec2884630a02670476f752d2a2bf1c395d8c8d6d768090ed95b4f04b0d5863/GOCO06s.gfc)).

A robust parser should key on `end_of_head`, treat unknown header lines as comments, and infer the column count from the record key plus the `errors` and `format` values.

**GRACE Level-2 files** have a YAML header followed by records of the form `GRCOF2 L M C S σC σS t_start t_end`. The file name encodes the product, mission, centre, observation type, truncation degree and release. In `BA01`, for example, `B` means MWI range-rate, `A` means 60×60, and `01` means unconstrained ([CSR RL06.3 notes](https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-docs/gracefo/open/docs/GRACE-FO_L2-ReleaseNotes_CSR_RL06.3.txt)).

**TN-14** is fixed-width Fortran with `NaN` for missing C30.

GRCOF2 and `gfc` are near-isomorphic, so one coefficient-set type can serve both. The exact RL06.3 YAML key names were not inspected, because the files sit behind Earthdata Login.

Licensing is permissive across the board. NASA collections carry `UseConstraints: None` and ask users to cite the dataset DOI, for example **10.5067/TEMSC-3JC634** for JPL mascons ([CMR UMM](https://cmr.earthdata.nasa.gov/search/collections.umm_json?short_name=TELLUS_GRAC-GRFO_MASCON_CRI_GRID_RL06.3_V4&provider=POCLOUD)). COST-G and AOD1B RL07 are **CC BY 4.0** ([COST-G header](https://icgem.gfz.de/getseries/02_COST-G_/Grace-Grace-FO_RL02.1/GSM-2_2026182-2026212_GRFO_COSTG_BF01_0201.gfc)). ESA data are free but require the credit line "Data provided by the European Space Agency" ([ESA T&C](https://earth.esa.int/eogateway/documents/20142/1564626/Terms-and-Conditions-for-the-use-of-ESA-Data.pdf)). The licences for ITSG, CSR, GSFC and Swarm EGF products were not stated on the pages fetched. An app that caches them should show their citations and avoid claiming a licence.

## Quantum gravimetry will not reach orbit before 2030

Every space gravity mission so far has used electrostatic accelerometers. The quantum record in orbit is short.

The **China Space Station Atom Interferometer** (launched November 2022) ran with an interrogation time of **T = 75 ms** and reached a single-shot acceleration resolution of **1 µm/s²** ([Li et al. 2025](https://academic.oup.com/nsr/article/12/4/nwaf012/7951977)). NASA's Cold Atom Lab on the ISS ran pathfinder Mach–Zehnder and Bragg-recoil experiments ([Williams et al. 2024, arXiv:2402.14685](https://www.arxiv.org/abs/2402.14685)). Its successor BECCAL has slipped to **2027** ([DLR, April 2025](https://www.dlr.de/en/latest/news/2025/dlr-at-the-space-symposium-2025)). Every mission study assumes T of **1–5 s**, roughly 13–67 times longer than anything demonstrated in orbit. Closing that gap is what the pathfinders are for.

**CARIOQA** is Europe's single-axis BEC accelerometer pathfinder. It entered **Phase B in October 2025** for 24 months, aiming at TRL 6, under CNES and DLR leadership with 15 partner organisations ([Biskupek et al., EGU26-9366](https://meetingorganizer.copernicus.org/EGU26/EGU26-9366.html)). A 2026 Springer chapter places a dedicated satellite in the **early 2030s**, with a three-year life and simulations at 500 km ([Schilling et al. 2026](https://link.springer.com/chapter/10.1007/1345_2026_314)). That is lower than the 600 km sun-synchronous orbit of the 2022 definition. The instrument target is a quantum-projection-noise floor below **10⁻¹⁰ m s⁻² Hz⁻¹ᐟ²**, with T = 1 s, contrast about 0.8, 10⁵ atoms and a 3.5 s cycle ([Lévèque et al. 2022](https://arxiv.org/pdf/2211.01215)).

NASA's **Quantum Gravity Gradiometer Pathfinder** (QGGPf, JPL-led) is a single-axis, cross-track, two-source rubidium BEC gradiometer. Its target is **below 100 mE/√Hz from 0.1 mHz to 0.1 Hz**, in a package of 150 kg and 350 W ([NASA QGGPf overview 2025](https://assets.science.nasa.gov/content/dam/science/cds/science-enabling-technology/events/2025/quantum-day/4_SMD_Quantum_Day_Quantum_Gravity_Gradiometer.pdf)). Infleqtion won a **$20 M follow-on on 27 August 2026**, bringing NASA's investment to $40 M and targeting a **2030 LEO launch** ([Infleqtion](https://ir.infleqtion.com/news-events/press-releases/detail/207/nasa-awards-infleqtion-20-million-as-worlds-first-quantum-gravity-mission-advances-toward-flight)).

Classical successors fly first. GRACE-C is on schedule for **December 2028**, using laser ranging only ([GSTM2026 A.3](https://meetingorganizer.copernicus.org/GSTM2026/session/59410)). ESA's NGGM targets **2032**, and together the two form the MAGIC constellation ([ESA LPS25](https://lps25.esa.int/lps25-presentations/presentations/2516/_2516.pdf)). ESA's roadmap places a classical-plus-quantum hybrid after MAGIC, with a full quantum mission after that ([Carraz et al., EGU23-14485](https://meetingorganizer.copernicus.org/EGU23/EGU23-14485.html)).

The simulation literature contains the result that most shapes what a Rust tool should compute: **temporal aliasing, not sensor noise, limits today's gravimetry**. Zingerle et al. 2024 ran closed-loop simulations with realistic background models. For a satellite pair, sensor noise only becomes visible once the time-variable signal is scaled down to under **10% with GRACE-FO-class sensors**, and under **0.1% with a 10⁻¹² m s⁻² Hz⁻¹ᐟ² atom accelerometer** ([Zingerle et al. 2024](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024EA003630)). Zahzam et al. likewise found that the benefit of a hybrid accelerometer "largely vanishes" once aliasing is included, unless the constellation has multiple pairs ([Zahzam et al. 2022](https://arxiv.org/pdf/2206.00634)).

The quantum advantage shows up instead in three places:

- **Gradiometry of the static field.** A future atom gradiometer raises the resolvable degree from **223 to 263** (about 90 km to 76 km) relative to GOCE ([Zingerle et al. 2024](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024EA003630)).
- **Multi-pair constellations.** In ESA's QSG4EMT study, weekly recoverable degree rises from about 15 for one pair to about 70 for six pairs ([Zingerle, Gruber, Pail & Daras 2024, DOAJ](https://doaj.org/article/eca7c3c5bd254f5696e933d48d968dc8)).
- **Basin-scale water storage.** For 70% of basins, the TWS error threshold falls from 19.5 cm (GRACE-FO) to 2.3 cm (MAGIC), 1.5 cm (CAI3) and 1.0 cm (CAI6) ([Kusche et al. 2025](https://academic.oup.com/gji/article/242/2/ggaf195/8151505)).

The EU–Japan FEATS-G project, coordinated by GFZ from 2026, will extend these end-to-end simulations ([GFZ FEATS-G](https://www.gfz.de/en/section/earth-system-modelling/projects/feats-g-quantifying-the-potential-of-quantum-space-gravimetry)).

For a simulator, these studies reduce to a small set of presets:

| Preset | Sensor | White-noise floor | Key parameters | Source |
|---|---|---|---|---|
| CARIOQA-like | 1-axis along-track CAI accelerometer | 3×10⁻¹⁰ (conservative) / 8×10⁻¹¹ m s⁻² Hz⁻¹ᐟ² | T = 1 s, C ≈ 0.8, N = 10⁵–10⁶, T_c = 3.5 s, ~500 km | [Lévèque 2022](https://arxiv.org/pdf/2211.01215) |
| QGGPf-like | 1-axis cross-track CAI gradiometer | < 100 mE/√Hz, 0.1 mHz–0.1 Hz | BEC, ~1 m sensor head | [NASA 2025](https://assets.science.nasa.gov/content/dam/science/cds/science-enabling-technology/events/2025/quantum-day/4_SMD_Quantum_Day_Quantum_Gravity_Gradiometer.pdf) |
| GOCE-class CAI-GG | 3-axis CAI gradiometer | 5 mE/√Hz | T = 5 s, D = 0.5 m, 239 km | [Trimeche 2019](https://arxiv.org/pdf/1903.09828) |
| Zingerle CAI11 / CAI12 | Hybrid EA + CAI accelerometer | 10⁻¹¹ / 10⁻¹² m s⁻² Hz⁻¹ᐟ² | with LRI_NG33 ranging | [Zingerle 2024](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024EA003630) |
| GRICE | Laser-linked CAI pair | 6×10⁻¹⁰ m s⁻² Hz⁻¹ᐟ² | T = 0.5 s, 373 km | [Lévèque 2021](https://arxiv.org/pdf/2011.03382) |

## One harmonic engine and four conversion layers carry every formula

The geodesy has been settled for decades: IERS 2010, Barthelmes' ICGEM report, Wahr et al. 1998 and Holmes & Featherstone 2002. What changed in 2024–2026 is the sensor noise models. The architecture follows from that. Build one spherical-harmonic (SH) synthesis that returns the potential and its derivatives, then layer degree-dependent factors on top for gravity functionals, water storage, error spectra and sensors. Formulas below marked *derived* come from the research notes rather than a fetched source, and need a numerical check before use.

### The geopotential sum and its stable recursion

Everything starts from the IERS expansion in **geocentric** spherical coordinates, using 4π-normalised associated Legendre functions (ALFs) without the Condon–Shortley phase ([IERS 2010, ch. 6](https://iers-conventions.obspm.fr/content/chapter6/icc6.pdf)):

$$V(r,\varphi,\lambda)=\frac{GM}{r}\sum_{n=0}^{N}\Big(\frac{R}{r}\Big)^{n}\sum_{m=0}^{n}\big[\bar C_{nm}\cos m\lambda+\bar S_{nm}\sin m\lambda\big]\,\bar P_{nm}(\sin\varphi),\qquad \bar P_{nm}=\sqrt{(2-\delta_{m0})(2n+1)\tfrac{(n-m)!}{(n+m)!}}\;P_{nm}$$

Geodetic input $(\phi,\lambda,h)$ must first go through ECEF coordinates, with $X=(N_\phi+h)\cos\phi\cos\lambda$ and $Z=((1-e^2)N_\phi+h)\sin\phi$. Take the geocentric latitude from `atan2(Z, √(X²+Y²))` rather than `asin`, to keep precision near the poles ([Wikipedia: geographic coordinate conversion](https://en.wikipedia.org/wiki/Geographic_coordinate_conversion)). With $t=\sin\varphi$ and $u=\cos\varphi$, the ALFs follow a column recursion ([gravity-toolkit](https://raw.githubusercontent.com/tsutterley/gravity-toolkit/main/gravity_toolkit/associated_legendre.py)):

$$\bar P_{00}=1,\quad \bar P_{11}=\sqrt3\,u,\quad \bar P_{mm}=u\sqrt{\tfrac{2m+1}{2m}}\,\bar P_{m-1,m-1},\quad \bar P_{m+1,m}=\sqrt{2m+3}\;t\,\bar P_{mm},\quad \bar P_{nm}=a_{nm}\,t\,\bar P_{n-1,m}-b_{nm}\,\bar P_{n-2,m}$$

$$a_{nm}=\sqrt{\frac{(2n-1)(2n+1)}{(n-m)(n+m)}},\qquad b_{nm}=\sqrt{\frac{(2n+1)(n+m-1)(n-m-1)}{(n-m)(n+m)(2n-3)}}$$

Near the poles, $u^m$ underflows at high order. Holmes & Featherstone avoid this by recursing on $\bar P_{nm}/u^m$, seeding sectorals with a scale factor of about $10^{-280}$, and applying $u^m$ by Horner's scheme over order. This stays complete to **degree 2700 at all latitudes**, except for first derivatives exactly at a pole ([Holmes & Featherstone 2002](https://www.doi.org/10.1007/S00190-002-0216-2)). Beyond about 2700, extended-exponent "X-numbers" are needed ([Fukushima 2012 summary](https://science.ecosyste.ms/projects/28565)). The per-order sum $\sum_n c_n\bar P_{nm}$, with $c_n=(R/r)^n\bar C_{nm}$, collapses to a Clenshaw recurrence (*derived*, hand-checked):

$$y_{N+1}=y_{N+2}=0,\qquad y_n=c_n+a_{n+1,m}\,t\,y_{n+1}-b_{n+2,m}\,y_{n+2},\qquad \sum_{n=m}^{N}c_n\bar P_{nm}=y_m\,\bar P_{mm}$$

For polar orbits, the derivative with respect to colatitude should avoid the $1/u$ factor (*derived*). Use $\tfrac{d\bar P_{nm}}{d\theta}=\tfrac12\big[k_m\sqrt{(n+m)(n-m+1)}\,\bar P_{n,m-1}-\sqrt{(n-m)(n+m+1)}\,\bar P_{n,m+1}\big]$, with $k_1=\sqrt2$, $k_{m\ge2}=1$, and $d\bar P_{n0}/d\theta=-\sqrt{n(n+1)/2}\,\bar P_{n1}$.

Grids reuse this work. For each latitude ring, compute $a_m=\sum_n c_n\bar P_{nm}$ (and $b_m$ from the sine coefficients), form $F_m=a_m-ib_m$, and run one real inverse FFT along longitude. That costs $O(L^3)$ for the whole grid, against $O(L^4)$ point by point. On a Gauss–Legendre grid of $(L+1)\times(2L+1)$ nodes, analysis is exact ([SHTOOLS grid formats](https://shtools.github.io/SHTOOLS/grid-formats.html)):

$$\bar C_{nm}=\frac{1}{4\pi}\sum_i w_i\,\bar P_{nm}(t_i)\Big[\frac{2\pi}{N_\lambda}\sum_j f_{ij}\cos m\lambda_j\Big]$$

Two benchmarks justify skipping exotic transforms. SHTns showed that vectorised $O(L^3)$ Gauss–Legendre transforms beat asymptotically fast ones **up to L = 1023** ([Schaeffer 2013](https://arxiv.org/abs/1202.6522)). GeographicLib's ring reuse cuts EGM2008 evaluation from **78 ms to 92 µs per point** ([GeographicLib gravity](https://geographiclib.sourceforge.io/C++/doc/gravity.html)).

### From potential to geoid, gravity and gradients

ICGEM computes every functional from the disturbing potential $T=W-U$ ([Barthelmes STR09/02 rev. 2013](https://icgem.gfz.de/docs/str-0902-revised.pdf)). The normal field has only even zonal terms, $\bar C^U_{n0}=-J_n/\sqrt{2n+1}$. It must be rescaled to the model's constants before subtraction:

$$C^U_{\ell m}=\hat C^U_{\ell m}\,\frac{GM^U}{GM}\Big(\frac{R^U}{R}\Big)^{\ell}$$

The height anomaly on the ellipsoid (Bruns, eq. 118), and the spherical-approximation gravity disturbance and anomaly (eqs. 125–126), are:

$$\zeta=\frac{GM}{r_e\,\gamma}\sum_{\ell}\Big(\frac{R}{r_e}\Big)^{\ell}\sum_m\bar P_{\ell m}\big(C^T_{\ell m}\cos m\lambda+S^T_{\ell m}\sin m\lambda\big),\qquad \delta g_{sa},\ \Delta g_{sa}=\frac{GM}{r^2}\sum_\ell\Big(\frac Rr\Big)^{\ell}\{(\ell+1),\ (\ell-1)\}\sum_m\bar P_{\ell m}\big(C^T_{\ell m}\cos m\lambda+S^T_{\ell m}\sin m\lambda\big)$$

The 2009 edition of the report prints these with a **wrong leading minus sign**; the 2013 revision is correct ([2009 edition](https://icgem.gfz.de/docs/str-0902.pdf)). Barthelmes warns that the spherical approximation "is not accurate enough for most practical applications", with errors of 0.69 mGal wrms. The exact form, $|\nabla W|-|\nabla U|$, uses the radial, latitude and longitude partials of the same sum plus the centrifugal terms $\Phi_r=\omega^2r\cos^2\varphi$ and $\Phi_\varphi=-\omega^2r^2\cos\varphi\sin\varphi$.

Normal gravity on the ellipsoid follows Somigliana's closed form, $\gamma=\gamma_e(1+k\sin^2\phi)/\sqrt{1-e^2\sin^2\phi}$. For GRS80, $\gamma_e=9.7803267715$, $\gamma_p=9.8321863685$ m s⁻², $k=0.001931851353$ and $U_0=62636860.850$ m² s⁻² ([Moritz, GRS80](https://ciencias.ulisboa.pt/sites/default/files/fcul/dep/dqb/doc/GRS80_Moritz.pdf)). Two convention constants decide whether results match other tools:

- **The zero-degree term**, about **−0.41 m** for EGM2008 relative to WGS84 ([GeographicLib](https://geographiclib.sourceforge.io/C++/doc/gravity.html)).
- **The tide-system offset** of $\bar C_{20}$: zero-tide minus tide-free is **−4.1736×10⁻⁹** for EGM2008 ([IERS 2010](https://iers-conventions.obspm.fr/content/chapter6/icc6.pdf)).

For gradiometry, the local north–east–up tensor (*derived*) is:

$$V_{zz}=V_{rr},\quad V_{xx}=\tfrac{V_r}{r}+\tfrac{V_{\varphi\varphi}}{r^2},\quad V_{yy}=\tfrac{V_r}{r}-\tfrac{\tan\varphi}{r^2}V_\varphi+\tfrac{V_{\lambda\lambda}}{r^2\cos^2\varphi}$$

Its trace vanishes outside the masses, which gives a free unit test. Its $1/\cos\varphi$ singularities are why polar-orbit work uses the Petrovskaya–Vershkov non-singular forms or a Cartesian Hessian ([Petrovskaya & Vershkov 2007](https://link.springer.com/article/10.1007/s00190-006-0087-z)).

### From monthly coefficients to water

Wahr et al. (1998, eq. 14) convert coefficient changes to surface mass. Dividing by $\rho_w=1000$ kg/m³ gives equivalent water height ([Wahr et al. 1998](https://geodynamicsprogram.whoi.edu/wp-content/uploads/sites/14/2018/10/wahr_et_al.1998_117046.pdf)):

$$\Delta\sigma(\theta,\phi)=\frac{a\rho_{ave}}{3}\sum_{l}\sum_{m}\bar P_{lm}(\cos\theta)\,\frac{2l+1}{1+k_l}\,\big(\Delta C_{lm}\cos m\phi+\Delta S_{lm}\sin m\phi\big),\qquad \rho_{ave}=5517\ \mathrm{kg/m^3}$$

The PREM load Love numbers are $k_2=-0.303$, $k_3=-0.194$, $k_4=-0.132$, $k_{10}=-0.069$, $k_{20}=-0.051$, $k_{50}=-0.027$, $k_{100}=-0.014$ and $k_{200}=-0.007$. Linear interpolation stays within 0.05%. In the centre-of-figure frame $k_1=+0.027$; in the centre-of-mass frame it is −1.

Smoothing uses Jekeli's Gaussian, inserted as $2\pi W_l$:

$$W_0=\tfrac{1}{2\pi},\quad W_1=\tfrac{1}{2\pi}\Big[\tfrac{1+e^{-2b}}{1-e^{-2b}}-\tfrac1b\Big],\quad W_{l+1}=-\tfrac{2l+1}{b}W_l+W_{l-1},\quad b=\tfrac{\ln2}{1-\cos(r/a)}$$

The recursion is unstable at high $l$, which is why gravity-toolkit clamps it at $10^{-10}$ ([gauss_weights.py](https://raw.githubusercontent.com/tsutterley/gravity-toolkit/main/gravity_toolkit/gauss_weights.py)). Stripes are removed in one of two ways:

- **Swenson–Wahr destriping:** a moving quadratic fit across degree for each order, even and odd degrees separately, with window half-width $\mathrm{round}(15e^{-m/10})$ and applied for $m\ge5$ ([destripe_harmonics.py](https://raw.githubusercontent.com/tsutterley/gravity-toolkit/main/gravity_toolkit/destripe_harmonics.py)).
- **DDK decorrelation:** $W_a=(E^{-1}+aS^{-1})^{-1}E^{-1}$ with signal power $l^{-4}$, where $a$ runs from 10¹⁴ (DDK1) to 5×10⁹ (DDK8) ([ICGEM DDK note](https://icgem.gfz.de/docs/Note_on_the_DDK_filters.pdf)).

The monthly pipeline then runs in this order:

1. Rescale to common $GM$ and $R$.
2. Replace $C_{20}$ always, and $C_{30}$ after August 2016, with TN-14 values. Insert degree 1 from TN-13 ([TN-14](https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-docs/gracefo/open/docs/TN-14_C30_C20_GSFC_SLR.txt); [JPL geocenter](https://grace.jpl.nasa.gov/data/get_data/geocenter)).
3. Subtract a reference-period mean.
4. Subtract a glacial isostatic adjustment rate such as ICE-6G_D. The choice of GIA model is the main driver of trend differences ([Peltier et al. 2012](https://www.atmosp.physics.utoronto.ca/%7epeltier/pubs_recent/Peltier,%20W.R.%20%20et%20al.,%20Comment%20on%20Chambers%20et%20al%20(2010),%20JGR-Solid%20Earth%20117,%20B11403,%202012.pdf)).
5. Filter.
6. Convert with Wahr's factor.
7. Synthesise.
8. Restore leakage with a gain factor $k=\arg\min\sum(S_T-kS_F)^2$ fitted on a model filtered the same way ([Landerer & Swenson 2012](https://grace.jpl.nasa.gov/system/downloadable_items/3_Landerer.2011WR011453.pdf)).

A basin mean can be computed directly from coefficients (*derived*). With $\vartheta^{c,s}_{lm}$ the 4π coefficients of the basin mask:

$$\overline{\mathrm{EWH}}=\frac{a\rho_{ave}}{3\rho_w\vartheta_{00}}\sum_{l,m}\frac{2l+1}{1+k_l}\tilde W_l\big(\vartheta^c_{lm}\Delta C_{lm}+\vartheta^s_{lm}\Delta S_{lm}\big)$$

Its variance is $\mathbf g^\top\Sigma\mathbf g$. Trends and seasonal amplitudes come from a least-squares fit of $c_0+c_1t+a_1\cos\omega t+b_1\sin\omega t+\dots$

### From coefficients to error spectra

Every comparison figure in the quantum literature is a degree spectrum. The degree variance is $c_n=\sum_m(\bar C_{nm}^2+\bar S_{nm}^2)$ and the per-coefficient degree RMS is $\sqrt{c_n/(2n+1)}$ ([SHTOOLS](https://shtools.github.io/SHTOOLS/real-spherical-harmonics.html); [Wahr 1998 eq. 24](https://geodynamicsprogram.whoi.edu/wp-content/uploads/sites/14/2018/10/wahr_et_al.1998_117046.pdf)). The degree amplitude in a chosen functional is $\kappa_n\sqrt{c_n}$, and the cumulative error is $[\sum_{n\le N}\kappa_n^2c_n]^{1/2}$, with these transfer factors at $r=R$ ([Barthelmes 2013](https://icgem.gfz.de/docs/str-0902-revised.pdf); [gravity-toolkit units](https://raw.githubusercontent.com/tsutterley/gravity-toolkit/main/gravity_toolkit/units.py)):

| Functional | $\kappa_n$ |
|---|---|
| Geoid height | $R$ |
| Gravity disturbance | $\frac{GM}{R^2}(n+1)$ |
| Gravity anomaly | $\frac{GM}{R^2}(n-1)$ |
| Radial gradient $T_{rr}$ | $\frac{GM}{R^3}(n+1)(n+2)$ |
| Equivalent water height | $\frac{R\rho_{ave}}{3\rho_w}\frac{2n+1}{1+k_n}$ |

At satellite altitude, multiply by $(R/r)^{n+1}$. Kaula's rule of about $10^{-5}/n^2$ per coefficient supplies a signal model for omission error ([AGU FM15](https://agu.confex.com/agu/fm15/webprogram/Paper80878.html)). The omission error comes out at roughly $63.7\,\mathrm{m}/N$ (*derived*), which is about 1 m at $N=60$.

Two mappings link sensors to degrees. The first is the observation frequency, $n\approx f\,T_{orbit}$: at 463 km, $f_{rev}\approx1.78\times10^{-4}$ Hz, so 10 mHz maps to about degree 57. This matches Encarnação et al. ([arXiv:2404.07835](https://arxiv.org/abs/2404.07835)). The second is resolution, $d_{min}\approx\pi R/\ell_{max}$ ([Zingerle 2024](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024EA003630)).

A fast formal-error path (*derived*) builds normal equations $N=A^\top WA$, with $W$ the inverse noise covariance from the sensor ASD (or a whitening filter), then reads off $\sigma(\bar C_{nm})=\sqrt{\mathrm{diag}(N^{-1})}$. White atom-interferometer noise makes $W$ diagonal, which is a real computational advantage over electrostatic sensors with $1/f$ noise.

### From atoms to acceleration noise

A space atom interferometer is a π/2–π–π/2 light-pulse sequence. Its phase is $\Delta\phi=n_k\,k_{eff}\,a\,T^2$, and it is read out as $P=P_0-\tfrac{C}{2}\cos\Delta\phi$. For ⁸⁷Rb at 780 nm, $k_{eff}=4\pi/\lambda\approx1.61\times10^7$ m⁻¹. In microgravity, $n_k=2$, because double diffraction is required there ([Lévèque 2021](https://arxiv.org/pdf/2011.03382); [Zahzam 2022](https://arxiv.org/pdf/2206.00634)). The per-shot and spectral noise follow from quantum projection noise, generalised by Encarnação et al. with an entanglement exponent $\alpha$ (0 for the standard quantum limit, 1 for the Heisenberg limit):

$$\sigma_\phi=\frac{1}{C\,N^{(1+\alpha)/2}},\qquad \sigma_{a}^{shot}=\frac{\sigma_\phi}{n_k\,k_{eff}\,T^2},\qquad \mathrm{ASD}_a=\sigma_a^{shot}\sqrt{T_c},\qquad T_c=\begin{cases}2T+T_{prep}&\text{sequential}\\ T_{prep}&\text{interleaved}\end{cases}$$

A gradiometer differences two clouds that share the same Raman beams, separated by a baseline $D$. Its phase is $\Delta\phi=n_k k_{eff}\Gamma DT^2$ and its noise is $\sigma_\Gamma=\sqrt2\,\sigma_\phi/(n_kk_{eff}DT^2)$ ([Trimeche 2019](https://arxiv.org/pdf/1903.09828)).

In frequency, the acceleration transfer function is a low-pass, $|H_a(\omega)|=(4n_kk_{eff}/\omega^2)\sin^2(\omega T/2)$ (*derived*), which tends to $n_kk_{eff}T^2$ as $\omega\to0$. Dead time aliases white acceleration noise into the phase as $\sigma_\Phi^2=(k^2T^4/2)(2T_c/3T-1)S_a/\tau$ ([Cheinet et al.](https://ar5iv.labs.arxiv.org/html/physics/0510197)).

Rotation is the hard constraint. Without a counter-rotating mirror, keeping contrast at nadir pointing (1.2 mrad/s) would need a cloud temperature below $T_0=m/(16k_Bk^2\Omega^2T^4)\approx$ **3 fK**, which is unachievable. The gravity-gradient bound $T_1=m/(4k_Bk^2T_{zz}^2T^6)$ is about **100 pK** ([Trimeche 2019](https://arxiv.org/pdf/1903.09828)). With active compensation, the phase reduces to $2k_{eff}T^2[a_x+x_0\Omega_y^2]$ ([HosseiniArani et al. 2024](https://arxiv.org/html/2405.11259v1)).

Hybrid accelerometers fuse the two sensor types in one of two ways:

- **Bias loop:** $b_{i+1}=b_i+G_b(a_{EA,i}-a_{CAI,i})$ ([Zahzam 2022](https://arxiv.org/pdf/2206.00634)).
- **Extended Kalman filter:** gain $K=P^-H/(R+H^2P^-)$ with $H=2k_{eff}T^2$ ([HosseiniArani 2024](https://arxiv.org/html/2405.11259v1)).

Neither form is what Zingerle used. Their "optimal combination" is unpublished, so inverse-variance PSD fusion, $S_{hyb}=(S_{EA}^{-1}+S_{CAI}^{-1})^{-1}$, is a stand-in assumption.

The mission-level noise curves come straight from Zingerle et al.'s printed equations. Units are m s⁻² Hz⁻¹ᐟ² for accelerometers and m Hz⁻¹ᐟ² for ranging ([Zingerle 2024](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024EA003630)):

| Eq. | Instrument | ASD |
|---|---|---|
| 1 | Alvarez@350 (drag-free) | $3\!\times\!10^{-8}f^2+4\!\times\!10^{-13}\sqrt{1+7\!\times\!10^{-4}/f+(3\!\times\!10^{-4}/f)^2}$ |
| 2 | Alvarez@500 | $3\!\times\!10^{-8}f^2+5\!\times\!10^{-13}\sqrt{1+(1/f)^{2/3}}$ |
| 3 | Dalin (ONERA) | $2.4\!\times\!10^{-13}f^{-1/4}+10^{-18}f^{-2}+8\!\times\!10^{-11}f^2$ |
| 4 | SuperSTAR (GRACE-FO) | $10^{-10}\sqrt{1+0.005/f}$ |
| 5 | MicroSTAR (NGGM) | $\tfrac{1}{\sqrt2}10^{-11}\sqrt{\frac{(10^{-3}/f)^4}{(10^{-5}/f)^4+1}+1+(f/0.1)^4}$ |
| — | CAI11 / CAI12 | $10^{-11}$ / $10^{-12}$ (white) |
| 6 | LRI GRACE-FO | $2.2\!\times\!10^{-10}/\sqrt f+10^{-12}/f^2$ |
| 7 | LRI 2030 | $2.2\!\times\!10^{-10}/\sqrt f+10^{-13}/f^2$ |
| 8 | LRI 2033 | $1.1\!\times\!10^{-10}/\sqrt f+5\!\times\!10^{-14}/f^2$ |

Combined "product" noise along the line of sight (*derived*, checked against Fig. 1) is $\sqrt{(\sqrt2\,\mathrm{acc})^2+((2\pi f)^2\,\mathrm{lri})^2}$. For gradiometers, the gradient noise is $\sqrt2\,n_{acc}/L$, plus centrifugal errors $-2\omega_y n_{\omega y}-2\omega_z n_{\omega z}$ from gyro noise, since $V_{xx}=\Gamma_{xx}+\omega_y^2+\omega_z^2$ ([Douch et al. 2018 poster](https://www.gis.uni-stuttgart.de/forschung/doc/Douch_2018.pdf)).

## Only one of the linked articles has figures worth reproducing

**The Open Access Government piece** is a two-page opinion article by Bart Root of TU Delft, published February 2026 under CC BY-NC-ND. It has **no data figures, equations or references** ([OAG](https://www.openaccessgovernment.org/article/why-quantum-satellite-gravimetry/205621/)). Its value is as context: the lineage from CHAMP to GRACE-FO to MAGIC, the argument that cold atoms cut drift and bias, and pointers to CARIOQA, NASA's gradiometer pathfinder and QSG4EMT. QSG4EMT is an ESA contract led by TU Munich, not a Horizon Europe project. Its 385-page final report is the deepest public source behind the article ([TUM QSG4EMT](https://asg.ed.tum.de/iapg/qsg4emt)).

**The phys.org item** (20 December 2018) is a NASA Goddard press release about an AOSense atom-interferometer gradiometer prototype. It is not Leibniz University Hannover work and cites no paper. Its only images are two photographs ([phys.org](https://phys.org/news/2018-12-team-quantum-sensor-satellite-gravimetry.html)). The numbers behind it come from a 2017 APS presentation: a terrestrial target of **8×10⁻² E/√Hz** with a 2 m baseline and T = 300 ms, a microgravity extrapolation to **7×10⁻⁵ E/√Hz**, and a claimed monthly 1 cm EWH at 200 km resolution ([NTRS 20170005272](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20170005272.pdf)). None of these have been demonstrated in orbit. The faithful way to visualise the article is a gradient-noise chart that places the AOSense levels next to GOCE (10–20 mE/√Hz), the Carraz/Trimeche concept (3.5–5 mE/√Hz), QGGPf (100 mE/√Hz), and the roughly **0.1 mE/√Hz** needed to see monthly time-variable gravity ([JRC EUR 30371](https://publications.jrc.ec.europa.eu/repository/bitstream/JRC121833/cold_atom_interferometry_satellite-based_gravimetry.pdf)).

**Zingerle et al. 2024** ("The Benefits of Future Quantum Accelerometers for Satellite Gravimetry", *Earth and Space Science* 11(9), CC BY 4.0) is the substantive source. It simulates:

- GRACE-like single pairs (89°, 463 km) and MAGIC-like double pairs (adding 70°, 432 km), recovered to degree 120 from 31-day or 7-day data.
- A GOCE-like gradiometer (96.7°, 250 km, 2 months at 5 s sampling) recovered to degree 180–280.

Forward models are GOCO05s, the ESA Earth System Model (ESM) AOHIS series, and EOT11a tides. The authors **did not share their simulator (IP rights)**. Supplementary data sit at mediaTUM, DOI 10.14459/2024MP1737747, but their contents could not be inspected ([Wiley full text](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024EA003630)). Its twelve figures fall into four tiers:

| Fig. | Content | Axes | How to reproduce | Tier |
|---|---|---|---|---|
| 1 | Accelerometer, LRI and combined ASDs | log f 1.4×10⁻⁴–0.1 Hz (top axis: SH degree); log ASD 5×10⁻¹³–10⁻⁹ | Evaluate Eqs. 1–8; ×√2 for accelerometers, ×(2πf)² for LRI; RSS; top axis n = f·T_orbit (~5620 s); CAI lines stop near 3×10⁻² Hz | **Exact** |
| 3 | SH error triangles + global EWH error maps (NGGM vs CAI11) | Order −120..120 vs degree 0..120, log colour 10⁻⁵–10² mm; map −2..2 mm, Robinson/Equal-Earth-like | Formal or simulated ΔC, ΔS → EWH factor; triangle with S at −m; synthesise to a 1° grid; diverging colour map | In kind |
| 5, 7 | Static-field geoid degree errors for the SGG classes; maximum degree 223 vs 263 | Degree 0–180/280; geoid m, log 10⁻⁵–10² | Normal-equation propagation of 5 s gradient noise on a 96.7°/250 km orbit; polar-gap orders excluded; crossing degree → πR/ℓ | In kind |
| 2, 4, 6 | SST degree errors (product-only vs full noise); aliasing attenuation 100%→0%; SGG time-variable | Degree 2–120; EWH mm, log | Closed loop: orbits, range-rates from GOCO05s + 6-hourly ESM + tides, coloured noise, de-aliasing, least squares | Simulator only |
| 8–12 | Basin RMSD vs area, box plot, extremes bars, score maps (Eq. 9), false positives | log area 2×10³–10⁷ km²; RMSD mm | 1° weekly TWS from HIS plus own error realisations; GRDC basins with cos φ weights; percentile thresholds; Score = (CP − IP)/N_ref | In kind (errors are ours) |

Several numeric checks tie a reproduction to the paper. All of them were computed from the printed equations:

- LRI_GraceFO as acceleration is **3.1×10⁻¹⁰** at 0.1 Hz and **4×10⁻¹¹** at 1 mHz.
- GRACE-FO combined noise at 1 mHz is **3.5×10⁻¹⁰**.
- CAI12 with LRI_NG33 is **2.5×10⁻¹²** at 1 mHz, flat because ranging dominates.
- The basin-size markers in Fig. 8, 2.8×10⁴ and 4.7×10⁴ km², equal $(\pi R/\ell)^2$ for ℓ = 120 and 92.

The paper's text refers to "Figure 1a/b" where it means Figure 2, and it never says whether its degree errors are amplitudes or RMS. Use Heller-Kaikov et al. 2023's convention.

The **Hannover work** the brief had in mind is Douch et al. 2018 (*Adv. Space Res.* 61(5), DOI 10.1016/j.asr.2017.12.005). It is closed-access and was read only through the JRC report, which reproduces its main figure. That figure shows geoid degree errors for a **5 mE/√Hz** atom gradiometer on a 239 km orbit over an 8-month extrapolated mission. A quasi-inertial 3-axis instrument beats the full 47-month GOCE mission, and nadir-pointing 1-axis $V_{yy}$ does not ([JRC report](https://publications.jrc.ec.europa.eu/repository/bitstream/JRC121833/cold_atom_interferometry_satellite-based_gravimetry.pdf)).

Its open sibling, Trimeche et al. 2019, describes the full pipeline: 71 days of real GOCE orbit, EIGEN-6C4 truth to degree 360, recovery to degree 240 with 58,077 unknowns, and √time scaling to 8 months. It reports about a **2× gain over GOCE above degree 50** ([arXiv:1903.09828](https://arxiv.org/pdf/1903.09828)). The Douch–Knabe December 2018 poster adds two more reproducible targets: the gradiometer error model and an SNR triangle, $\log_{10}(\mathrm{median}|C^{tvar}_{lm}|/\mathrm{err})$ ([poster](https://www.gis.uni-stuttgart.de/forschung/doc/Douch_2018.pdf)).

Every forward model these papers use is public. GOCO05s, EIGEN-6C4, EOT11a and GOT4.7 are on ICGEM. The updated ESA ESM (6-hourly, degree 180, 1995–2006) is at DOI 10.5880/GFZ.1.3.2014.001. GRDC river basins come from BfG under terms that were not checked. The ESM file names could not be confirmed because GFZ's servers refused connections ([GFZ ESM page](https://isdc.gfz.de/esmdata/esaesm/)).

## Companion crates around geo5d beat a bigger geo5d

geo5d v0.2.0 was published on 14 September 2026 under Apache-2.0 by crates.io user `dan-garden` (GitHub FusionVault). It is about 850 lines across 7 files, with **0 reverse dependencies and 29 downloads**, so breaking changes cost nothing now ([crates.io API](https://crates.io/api/v1/crates/geo5d)). Its `SphericalHarmonicModel` already does a lot:

- Holds 4π-normalised coefficients in triangular `n(n+1)/2+m` order.
- Synthesises point by point with a Holmes–Featherstone column recursion that returns first and second derivatives.
- Provides the geopotential, gravity, disturbing potential, quasigeoid, disturbance, spherical-approximation anomaly, deflections and an ENU gradient tensor.
- Parses `.gfc` files ([gravity.rs](https://github.com/FusionVault/geo5d/blob/main/src/gravity.rs)).

The changelog claims |∇V| agrees with pyshtools to **0.05 µGal** on EGM96. It also claims the geoid agrees with GeographicLib to "a uniform ~9 cm" ([CHANGELOG](https://github.com/FusionVault/geo5d/blob/main/CHANGELOG.md)). GeographicLib itself reproduces NGA's EGM2008 synthesis to **78 pm** ([GeographicLib](https://geographiclib.sourceforge.io/C++/doc/gravity.html)). A 9 cm offset is therefore a convention to explain, not agreement.

Four defects matter for satellite gravimetry:

- **Time-variable models are read wrongly, with no warning.** `from_gfc` reads `gfct` lines as static values and drops `trnd`, `asin` and `acos`. A model such as GOCO06s, which uses all four, is loaded incorrectly. The parser also ignores the `norm`, `tide_system` and `errors` headers.
- **Possible underflow.** The sectoral seed multiplies $u^m$ directly into a $10^{280}$-scaled value. That may underflow at high order near the poles, and it is validated only to degree 360.
- **Panics at the trust boundary.** Constructors use `assert!` where they should return `Result`.
- **Two copies of geo3d.** geo5d pulls geo3d 0.2 through geo4d 0.1.1, while geo3d is at 0.3.4. A crate using both gets two incompatible `Ecef` types ([crates.io geo3d](https://crates.io/api/v1/crates/geo3d)).

The underflow is a hypothesis; no code was run to confirm it.

No mature Rust crate covers geodetic SH gravity. satkit stops at degree 70 for orbit work. astrodyn_gravity is a JEOD port tied to Bevy. sphrs and apollo-sht use graphics or complex orthonormal conventions, and apollo-sht needs Rust 1.95. The `egm2008` and `egm96` crates are interpolation grids. A crates.io search for "mascon" returns nothing ([crates.io search](https://crates.io/api/v1/crates?q=mascon)). The niche is open, and geo5d's engine is the right seed. Folding grids, I/O, hydrology, sensors and GUI into geo5d would bury a clean geodesy library under C-backed and GPU dependencies. The better split is a small workspace in which geo5d stays the façade for normal gravity, heights, relativity and time transfer, and delegates synthesis to a new core:

| Crate (name free on 2026-10-22) | Scope | Key dependencies | MSRV |
|---|---|---|---|
| `geopotential` | no_std + alloc SH core: ALFs (Holmes–Featherstone with Horner; X-numbers later), normalisation and phase conversions, static plus time-variable coefficient sets, point/ring/grid synthesis, GLQ and Driscoll–Healy analysis, degree spectra and transfer factors | `libm`; features `fft` (realfft 3.5), `rayon` (1.12) | 1.85 |
| `gravity-field` | ICGEM gfc v1/v2 (all record keys), GRCOF2 + YAML, TN-13/TN-14, Love numbers, mascon NetCDF; optional fetchers for CMR, ICGEM, TU Graz, swarm-diss with Earthdata auth from env/`.netrc` | `netcdf-reader` 0.9 (pure Rust, default); features `netcdf-c` (netcdf 0.12), `hdf5` (hdf5-pure 0.47 for GSFC mascons, MSRV 1.89) | 1.85 |
| `gravimetry` | GRACE L2 chain: C20/C30/degree-1/GIA, Gaussian/DDK/destriping, EWH, basin means, trend/seasonal fits, gain factors, spherical-cap mascons | `geopotential`, `ndarray` 0.17 | 1.85 |
| `satgrav` | Mission and sensor simulation: Zingerle/EA/LRI ASDs, atom-interferometer models, hybrid fusion, noise synthesis from ASDs, normal-equation error propagation for SST and SGG | `faer` 0.24, `rustfft` | 1.85 |
| `geo5d` (existing) | Normal gravity, heights, relativity, transfer; depends on `geopotential` | geo4d (bumped to geo3d 0.3) | 1.85 |
| `apps/gravviz` (unpublished) | CLI figures + interactive viewer | plotters 0.3.7, proj4rs 0.2.1, eframe/egui 0.37, egui_plot 0.38, wgpu 30 | 1.95+ |

Sources for the table: crate metadata from the [crates.io API](https://crates.io/api/v1/crates/geopotential); pure-Rust NetCDF from [netcdf-rust README](https://github.com/roteiro-gis/netcdf-rust); viz MSRVs from [crates.io egui](https://crates.io/crates/egui).

The other free names are `gravity-field`, `gravimetry`, `satgrav`, `geopotential`, `shgrav`, `mascon`, `grace-l2`, `gravfield` and `geo6d`. The taken names are `gravity`, `grace`, `gfc` and `egm`; `egm` holds an empty template ([crates.io](https://crates.io/crates/egm)). Avoid names that look like official ports (`shtools`, `shtns`, `gravity-toolkit`, `icgem`). Publish real 0.1.0 crates rather than placeholders, and re-check availability on the day you publish.

Licensing shapes what can be ported. gravity-toolkit (MIT, v1.2.8, October 2026) is the closest functional match for the GRACE chain ([polargeodesy/gravity-toolkit](https://github.com/polargeodesy/gravity-toolkit)). GeographicLib's GravityModel/GravityCircle (MIT) gives fast ring synthesis with well-documented conventions. SHTOOLS (BSD-3) covers grids and localisation ([SHTOOLS](https://github.com/SHTOOLS/SHTOOLS)). Port from these, with NOTICE attribution. Use GROOPS (GPL-3.0), ducc0 (GPL-2.0), SHTns (CeCILL) and the ICGEM calculation service only as black-box oracles. The ICGEM service currently says it accepts no new calculations because of load ([ICGEM calcgrid](https://icgem.gfz.de/calcgrid)). Keep nyx-space (AGPL) and the GPL colour-map crate `prismatica` out of the dependency tree. Dual-license the new crates **MIT OR Apache-2.0** and relicense geo5d the same way while it has no users. Pin the pure-Rust NetCDF readers, which are about seven months old, and round-trip test them against real mascon files.

The visualisation stack should split along the line between publication and exploration. **Static figures** follow one route. Synthesise a grid in Rust. For each output pixel, inverse-project with proj4rs (`robin`, `eqearth` or `moll`, which its docs mark as spherical). Colour with embedded Crameri scientific maps (`vik` for diverging errors, `batlow` for sequential; MIT per secondary sources). Composite axes, colour bars and log-scale spectra in plotters. Also write CF NetCDF and GeoTIFF (`netcdf-writer`, `geotiff-writer`) so that PyGMT or Cartopy can produce journal-grade coastlines and graticules without rebuilding them in Rust ([proj4rs projections](https://github.com/3liz/proj4rs/blob/main/PROJECTIONS.md); [Crameri via colorlab](https://github.com/beddalumia/colorlab/blob/master/crameri/README.md)). That means deliberately not copying Zingerle's jet-like colour map, which is perceptually misleading.

**Interactive exploration** uses eframe/egui with egui_plot for spectra and time series, plus a wgpu textured sphere for the globe. It runs natively and as WASM. Rerun gives an instant developer view of fields and orbit tracks. If the goal becomes a browser globe with basemaps, keep Rust as a WASM compute core feeding CesiumJS (Apache-2.0) or MapLibre GL JS (globe since v5.0) ([MapLibre roadmap](https://maplibre.org/roadmap/maplibre-gl-js/globe-view)). Those GUI crates need Rust **1.95–1.97**, which is the main reason the viewer must live outside the 1.85-MSRV libraries ([crates.io bevy](https://crates.io/crates/bevy)).

## Five gated phases take the toolkit from parser to globe

Each phase ends with tests that compare against closed forms, published constants or independent implementations, never against our own earlier output. The first two phases need no network data at all.

| Phase | Deliverable | Validation gate (must pass before the next phase) |
|---|---|---|
| 0. Groundwork | Workspace with `resolver = "3"`; bump geo4d/geo5d to geo3d 0.3; `Result` instead of `assert!`; dual licence; CI on MSRV 1.85 and stable 1.99.1; `cargo-semver-checks`; crates.io Trusted Publishing; real 0.1.0 releases to claim names | geo5d's existing tests pass on geo3d 0.3; `cargo tree -d` shows one geo3d |
| 1. `geopotential` | ALFs, Clenshaw, pole-safe derivatives, time-variable coefficient sets, point/ring/grid synthesis, GLQ/DH analysis, spectra; geo5d delegates to it | Closed forms $\bar P_{20}=\tfrac{\sqrt5}{2}(3t^2-1)$, $\bar P_{22}=\tfrac{\sqrt{15}}{2}u^2$; sum rule $\sum_m\bar P_{nm}^2=2n+1$ to n = 2700 at 0°, 45°, 89.9°, 89.9999°; Laplace trace ≈ 0; GRS80 $U_0$, $\gamma_e$, $\gamma_p$ to 10 digits; GLQ round trip ≈ 10⁻¹³; EGM2008 geoid vs GeographicLib at sub-mm after the −0.41 m and tide conventions (explaining the 9 cm); IERS $\bar C_{20}$ offset −4.1736×10⁻⁹ |
| 2. `gravity-field` | Tolerant gfc/GRCOF2/TN parsers, Love-number tables, mascon NetCDF reader, catalogue fetchers behind a `fetch` feature | Real fixtures parse: GOCO06s, COST-G Jul 2026 (`errors empirical`, `+…E+01`), Swarm Dec 2025 ZIP (stray `WARNING`), TN-14 with `NaN` C30, the TN-13 GRCOF2 line; gfct/trnd/asin/acos time evaluation matches hand computation; JPL mascon dims 258×360×720; live CMR/swarm-diss tests `#[ignore]`d by default |
| 3. `gravimetry` | Corrections, filters, EWH, basin means, fits, gain factors, cap mascons | $W_0=1/2\pi$; Love interpolation < 0.05%; same GSM file through gravity-toolkit matches; reproduce JPL's 1° GRD-3 land grid for Aug 2026 (300 km Gaussian, destriped, ICE-6G_D, TN-14, Sun 2016 geocentre); basin mean from coefficients equals cos φ-weighted grid integration; cap coefficients equal numerical quadrature |
| 4. `satgrav` | ASD library, atom-interferometer models, hybrid fusion, ASD-shaped noise, normal-equation error propagation for gradiometry and simplified SST | Zingerle Fig. 1 checkpoints (3.1×10⁻¹⁰, 4×10⁻¹¹, 3.5×10⁻¹⁰, 2.5×10⁻¹²); CARIOQA parameters give 2.3×10⁻¹⁰ (paper: 3×10⁻¹⁰ conservative) and 7.3×10⁻¹¹ at N = 10⁶ (paper: 8×10⁻¹¹); GRICE gives 5.8×10⁻¹⁰ (paper: 6×10⁻¹⁰); Trimeche gives 3.5 mE/shot and a 1087 rad gradient phase; d_min = 90/76 km at degrees 223/263; 5 mE/√Hz at 239 km lands near Trimeche's ~2× over GOCE above degree 50 |
| 5. `gravviz` | Fig. 1 replica, degree-error and triangle plots, projected EWH maps, GRACE-FO mascon time slider, basin time series, egui globe | Numeric snapshot tests on the rendered grids; same NetCDF rendered by PyGMT matches visually; Fig. 3 analogue shows NGGM striping fading under CAI11 |

The CARIOQA, GRICE and Trimeche checks were computed from the papers' own parameters ([Lévèque 2022](https://arxiv.org/pdf/2211.01215); [Lévèque 2021](https://arxiv.org/pdf/2011.03382); [Trimeche 2019](https://arxiv.org/pdf/1903.09828)). GRICE's "2.3 mrad per shot" is internally inconsistent; use 3.3 mrad.

An optional sixth phase, a closed-loop SST simulator, is the only route to Zingerle's Figs. 2, 4 and 6. It needs orbit propagation, range-rate forward modelling from 6-hourly ESM fields, de-aliasing, and arc-wise noise covariances, which is months of work. GROOPS can serve as a GPL black-box reference but not as a source to translate. The semi-analytic path in Phase 4 reproduces the product-only curves, though not aliasing, which is the main effect those figures demonstrate. Two engineering cautions apply throughout. First, ship small fixtures (a truncated EGM96, one GRACE month) under `tests/data` rather than bundling EGM2008's 4.7 M coefficients. Second, run the first builds of the GPU stack and any C-backed `netcdf-src` feature under a resource cap, since those are the heaviest compiles in the tree.

## Conclusion

The question "where is the quantum satellite data?" has a clear answer: it does not exist yet, and the first pathfinder data, around 2030–2031, will come from single-axis instruments. CARIOQA's 10⁻¹⁰ m s⁻² Hz⁻¹ᐟ² target equals the GRACE-FO accelerometer's floor ([Zahzam 2022](https://arxiv.org/pdf/2206.00634)); its gain is drift-free stability, not a lower noise floor. Building for it is still not premature. The binding constraint the papers reveal is temporal aliasing, not sensor noise, so the most valuable thing a Rust toolkit can offer is not a better atom-interferometer formula. It is honest error spectra that put sensor noise, ranging noise and aliased signal on one degree axis. A tool built that way can show why a 10⁻¹² m s⁻² Hz⁻¹ᐟ² accelerometer barely helps a single pair, and why multi-pair constellations and gradiometry are where quantum sensors earn their cost. Few public tools make that argument visible, and Zingerle's own simulator is closed.

The investment compounds because the real-data and simulated branches share one engine. The coefficient container that loads an August 2026 GRACE-FO month is the same one that holds simulated CAI12 errors. The EWH factor is the same in both. CARIOQA or QGGPf results will most likely arrive as Level-2 coefficients or accelerometer series in the existing SDS and ICGEM formats, as GRACE-C's will. A crate family that parses those formats correctly, including the time-variable `gfct`/`trnd` records geo5d currently drops, will be ready the day the first quantum field is published. The Rust ecosystem has no geodetic spherical-harmonic gravity crate and no mascon crate today, so finishing that crate family carefully would fill a real gap.
