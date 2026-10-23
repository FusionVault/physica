# Satellite gravimetry methods and missions: the 2023–2026 research frontier (status 24 October 2026)

Conventions used in every bullet:
- **Status** tags: `[established]` = several independent groups, operational or used in products; `[emerging]` = 1–3 groups, published or preprint, not yet adopted operationally; `[open problem]` = gap with no convincing solution yet. `[foundational, pre-2023]` marks older work cited only for context.
- **Open?** = whether code and/or data are public, as stated on the source page ("not stated" means the page I read gave no code/data statement).
- The previous report (GRACE/GOCE data access, CARIOQA/QGGPf basics, Zingerle et al. 2024) is not repeated here.
- Research limit: the shared web-search budget ran out partway through. Later items were checked by fetching known primary URLs, not by new searches. The Gaps sections list what is still unverified.

---

## 1. Machine learning in satellite gravimetry (downscaling, gap filling, pre-2002 reconstruction, forecasting, denoising, neural fields, foundation models)

### Takeaway
ML work on GRACE is now mostly about water storage, not gravity. The most defensible results are: (a) reconstructions trained on GRACE that are validated against independent data such as streamflow, sea-level budgets or SLR/DORIS gravity, and (b) short-lead forecasts (1–3 months) that bridge the 40–60-day latency. Spatial downscaling below GRACE's ~300 km native resolution is still contested, because the fine-scale detail comes from the auxiliary hydrological model rather than from gravity. In 2025 a linear-regression baseline matched or beat LSTM and Temporal Fusion Transformer models for TWS prediction. I found no peer-reviewed deep-learning destriping of Level-2 spherical-harmonic (SH) fields. Neural-field representations of gravity exist only for small bodies and regional geoids, and no Earth-observation foundation model treats gravity as a first-class input.

### Cited Findings

**Spatial and temporal downscaling (status: emerging; spatial downscaling is contested)**
- Uz, Akyilmaz, Shum et al. (J. Hydrology, Dec 2024) downscaled monthly GRACE/GRACE-FO mascon TWSA to **daily** sampling over CONUS with deep learning. They validated against an independent observation and a land hydrology model, and say the daily output is also implicitly enhanced spatially. [emerging; Open? not stated] — [ScienceDirect S0022169424015907](https://www.sciencedirect.com/science/article/abs/pii/S0022169424015907)
- Uz, Atman, Akyılmaz, Shum (GRACE/GRACE-FO Science Team Meeting 2024, abstract GSTM2024-54) downscaled JPL mascon TWSA from ~300 km to **50 km** (Apr 2002–Dec 2022). The deep-learning model uses a "soft constrained" loss meant to keep the GRACE signal intact, and learns spatiotemporal correlations from WaterGAP (WGHM). It was checked against glacier loss, ENSO variability, InSAR subsidence and river altimetry. The abstract gives no numbers. [emerging; Open? not stated] — [GSTM2024-54](https://meetingorganizer.copernicus.org/GSTM2024/GSTM2024-54.html)
- Dombos, Gou, Soja (ETH), **D-SHIFT**, arXiv 2605.00652 (May 2026). The model is trained on monthly data to map low-resolution SH solutions to mascon-like fields, then applied to **daily SH solutions**. Monthly global mean RMSE is ~2.3 cm against mascons. The Greenland basin-mean trend is −10.5 cm/yr, against −12.0 cm/yr for CSR monthly, and coastal mass-loss patterns are reproduced better. [emerging; Open? no code/data statement on abstract page] — [arXiv:2605.00652](https://arxiv.org/abs/2605.00652)
- Regional CNN-LSTM-attention and spatiotemporal NN downscaling to 0.1° (Xinjiang, Tengger Desert) report correlation with well levels of only ~0.57. [emerging] — [J. Hydrology Regional Studies 2025](https://www.sciencedirect.com/science/article/pii/S2214581825008110); [Water 18(7):781](https://doi.org/10.3390/w18070781)
- **Critique:** a 2024 *Scientific Reports* study (Australian Water Outlook) argues that resampled 0.5°/0.25° GRACE grids are spatially correlated from the native ~3° resolution, so they "do not contain physical information at a finer scale". Finer estimates need assimilated local hydrological information. [open problem] — [Sci. Rep. 2024, s41598-024-60366-2](https://www.nature.com/articles/s41598-024-60366-2)
- A spectral-combination downscaling paper states that the output spatial resolution "is completely dependent upon the hydrological model". [open problem] — [J. Hydrology 2022, S0022169422012057](https://www.sciencedirect.com/science/article/pii/S0022169422012057)
- A 2026 systematic review of GRACE downscaling and ML for groundwater exists. [review] — [Hydrology 13(5):135](https://doi.org/10.3390/hydrology13050135)

**Filling the 2017–2018 gap (status: established as practice, but every method relies on something other than gravity)**
- The gap is about 11 months, Jul 2017–May 2018. [established] — [Remote Sens. 14(7):1565 (2022)](https://doi.org/10.3390/rs14071565)
- A transformer with causal-convolution self-attention fills the gap (Remote Sensing of Environment 2024). [emerging; Open? not stated] — [RSE S0034425724004917](https://www.sciencedirect.com/science/article/abs/pii/S0034425724004917)
- A two-step reconstruction (2025, Int. J. Digital Earth). [emerging] — [T&F 17538947.2025.2468418](https://www.tandfonline.com/doi/full/10.1080/17538947.2025.2468418)
- Singular-spectrum-analysis gap filling for downscaled TWSA (CIDR). [emerging] — [arXiv:2504.08454](https://arxiv.org/pdf/2504.08454)
- **Physics-based alternative:** combining Swarm high-low SST with SLR gives a continuous gravity series across the gap (J. Geod. 2024). Swarm is "the backbone" during the gap, and ocean RMS error rises ~20–30% without GRACE. A TU Graz processing change roughly doubled Swarm solution quality, but COST-G adopted the improved orbits only from 2020 onward. This bears directly on the gap months in our 259 COST-G fields. [established/emerging; Open? open-access paper] — [J. Geod. 2024, 10.1007/s00190-024-01888-5](https://link.springer.com/article/10.1007/s00190-024-01888-5)
- A 2026 bias-aware ML downscaling study (Bug River basin) finds the largest disagreement between ML imputation and ARIMA occurs during the 2017–2018 gap. [emerging] — [Remote Sens. 18(17):2909](https://doi.org/10.3390/rs18172909)

**Reconstruction before 2002 (status: emerging; trends and accelerations disagree between products)**
- **GTWS-MLrec** (Yin, Slater, Khouakhi, Yu, Liu, Li, Pokhrel, Gentine; ESSD 15, 5597–5615, 2023; doi:10.5194/essd-15-5597-2023) gives monthly 0.25° TWSA for **1940–2022**. Of MARS, GLM, ANN, GAM and RF, the RF and ANN models did best, using 11 ERA5 variables plus runoff, soil moisture, impervious area, LAI and SIF with a 3-month lag. Independent checks: annual streamflow at 10,168 gauges PCC ≈ 0.58–0.60; water balance in 341 basins PCC ≈ 0.58; sea-level-based ocean mass PCC > 0.7. [emerging; **Open: data** doi:10.5281/zenodo.10040927] — [ESSD](https://essd.copernicus.org/articles/15/5597/2023/)
- **TWSTORE** (Hacker, Gutknecht, Löcher, Kusche; ESSD 18, 1747–1781, 9 Mar 2026). Monthly fields from **1984** at 0.5°, effective resolution ~330 km. Climate-predictor models (ARX, multiple linear regression, random forest) are trained in a **GRACE-derived EOF basis**, then combined by **variance component estimation** with low-degree (≈ 5–10) **SLR + DORIS** gravity fields (Löcher et al. 2025). The authors warn that "GRACE-derived results should not be extrapolated to the past". Their trends and accelerations differ notably from other reconstructions (HG19, Li21, CHR22, PA25, MA25). [emerging, first use of pre-GRACE geodetic gravity in a reconstruction; **Open: data** doi:10.5281/zenodo.15827789 and 10.5281/zenodo.16643628] — [ESSD 18, 1747](https://essd.copernicus.org/articles/18/1747/2026/)
- **ML-TWiX** (Scientific Data 2026) extends TWSA back to 1980 with a **deep ensemble** that also produces uncertainty estimates. Details not verified (paywall redirect). [emerging] — [Sci. Data s41597-026-06604-w](https://www.nature.com/articles/s41597-026-06604-w)
- Chanda et al. (ESSD 17, 2575, 2025): CNN, SVR, extra-trees and stacking reconstruction of long-term global TWSA. [emerging] — [ESSD 17, 2575](https://essd.copernicus.org/articles/17/2575/2025/)
- Arzoumanidis, Johannsen, Middendorf, Eicker, Dehbi, arXiv 2606.23833 (Jun 2026): a spatiotemporal graph neural network (MTGNN) driven by daily ERA5 P, ET and runoff, reconstructing South America back to **1940**. Grid-cell r = 0.69, basin r = 0.94, within 0.025 correlation of the best baseline (GTWS-MLrec, RM-REC, GRAiCE) while using one-half to one-tenth as many predictors. All models are weak in arid regions. [emerging; **Open: code** github.com/hcu-cml/MTGNN-TWS-Reconstruction-GRACE] — [arXiv:2606.23833](https://arxiv.org/abs/2606.23833)

**Forecasting TWS (status: emerging; gains over simple baselines are contested)**
- Li et al. (GRL 2024, doi:10.1029/2024GL109101) used ML on prior hydrometeorological drivers to forecast non-seasonal GRACE TWS **up to 12 months ahead**, including exceptional droughts. [emerging] — [GRL](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2024GL109101)
- Li et al. (WRR 2026, doi:10.1029/2025WR041710): an observation-only forecast out to 12 months, evaluated 2010–2024 and benchmarked against ECMWF SEAS5. [emerging] — [WRR](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2025WR041710)
- Mo et al. (GRL 2025, doi:10.1029/2024GL112677): a **Bayesian CNN** nowcasts the 1–3 months of GRACE latency with r = 0.92–0.95, flagging hydrological drought up to 3 months before GRACE data arrive. [emerging] — [GRL](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2024GL112677?af=R)
- Li et al. (WRR 2025, doi:10.1029/2024WR037926): LSTM combined with CLM simulations and GRACE data assimilation for reanalysis and forecasting. [emerging] — [WRR](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2024WR037926?af=R)
- **Contested:** Nie, Kumar et al. (NASA GSFC), arXiv 2510.10799 (Oct 2025), "Rethinking deep learning: linear regression remains a key benchmark in predicting terrestrial water storage". On the open HydroGlobe dataset, linear regression **outperforms LSTM and TFT**. The authors say DL gains for land-surface states are "not guaranteed". [contested; data open (HydroGlobe), code not stated] — [arXiv:2510.10799](https://arxiv.org/abs/2510.10799)
- The GB-ClimaX reproducibility release (vision transformer for Great Basin TWS, 3-month lead, against LSTM, ConvLSTM and linear baselines) concludes that antecedent storage and hydrologic regime limit DL gains. [emerging; **Open: code/data** on Zenodo] — [Zenodo 22803285](https://zenodo.org/records/22803285)

**Denoising and destriping (status: no deep-learning SH destriping found; conventional filters are still advancing)**
- I found no paper applying CNN or U-Net deep learning to destripe GRACE Level-2 SH coefficients. Searches returned only general image/seismic U-Net work. [open problem] — [search result set incl. arXiv:2209.14973](https://arxiv.org/pdf/2209.14973)
- Stepwise iterative enhanced destriping (GJI 2025, ggaf336) combines Gaussian smoothing with **variational mode decomposition** and suppresses noise better than other filters, especially over oceans. [emerging; **Open: code** Zenodo 16625922] — [GJI 243(1) ggaf336](https://academic.oup.com/gji/article/243/1/ggaf336/8242174); [Zenodo 16625922](https://zenodo.org/records/16625922)
- A **physically informed spatial filter (PISP)** for destriping (GJI 2026, ggag048) models the structure of the striping noise instead of only suppressing it. [emerging] — [GJI 245(1) ggag048](https://academic.oup.com/gji/article/245/1/ggag048/8445404)
- Dual low-pass filtering of Level-2 stripes (Acta Geod. Geophys. 2024). [emerging] — [s40328-024-00448-9](https://link.springer.com/article/10.1007/s40328-024-00448-9)
- A latitude/longitude correlation spatial filter for N–S stripes (J. Geod. 2022). [foundational, pre-2023] — [s00190-022-01614-z](https://link.springer.com/article/10.1007/s00190-022-01614-z)
- In the ESA MAGIC studies, the time-variable decorrelation filter (**VADER/VDK**) usually beats DDK, and double-pair solutions need much weaker filtering than single-pair ones. [emerging] — [TUM IAPG MAGIC results](https://www.asg.ed.tum.de/iapg/magic/results-4)
- ML has been used at the instrument level: **Noise2Noise** unsupervised denoising of GRACE-FO accelerometer data beat traditional smoothing (Sensors 2023) [emerging] — [PMC10347034](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10347034/); an LSTM forecasting model for GRACE accelerometer data (2023) [emerging] — [arXiv:2308.08621](https://arxiv.org/pdf/2308.08621)

**Neural fields / implicit representations of gravity (status: emerging for small bodies; open problem for Earth time-variable gravity)**
- **PINN-GM-III**, Martin & Schaub, arXiv 2312.10257 (v3 Feb 2025, submitted to J. Astronautical Sci.), addresses feature divergence, low-altitude sampling bias, numerical instability and extrapolation error. It proposes six evaluation metrics and is tested on a heterogeneous-density asteroid. [emerging; **Open: code** GravNN repo] — [arXiv:2312.10257](https://arxiv.org/abs/2312.10257); [GitHub MartinAstro/GravNN](https://github.com/MartinAstro/GravNN)
- **MasconCube**, Fanti & Izzo (ESA ACT), arXiv 2509.08607 (v3 Apr 2026), is an *explicit* 3-D grid of point masses optimized self-supervised. It trains ~40× faster than GeodesyNets and is tested on Bennu, Eros and Itokawa. [emerging; code not linked on abstract page] — [arXiv:2509.08607](https://arxiv.org/abs/2509.08607)
- "Physics-informed neural networks for geoid modeling", J. Geod. (online 17 Jan 2026, doi:10.1007/s00190-025-02017-6), from the Chinese Academy of Surveying and Mapping. A likely precursor (EGU24-4068: Jiang, Tu, Dang) embeds the Laplace equation of the disturbing potential and the gravity-anomaly equation in a CNN loss (Colorado test). Authorship of the journal paper is not verified. [emerging] — [J. Geod.](https://link.springer.com/article/10.1007/s00190-025-02017-6); [CASM news](https://casm.ac.cn/kjcx/kjjz/202602/12/243878.html); [EGU24-4068](https://meetingorganizer.copernicus.org/EGU24/EGU24-4068.html)
- Implicit neural representations for 3-D gravity inversion (Computers & Geosciences 2025). [emerging] — [S0098300425002328](https://www.sciencedirect.com/science/article/abs/pii/S0098300425002328)
- PINN gravity models need "extensive computational resources and training time" (per the MasconCube paper and the PINN geoid literature). [established caveat] — [arXiv:2509.08607](https://arxiv.org/html/2509.08607)

**Foundation models (status: open problem)**
- **AlphaEarth Foundations** (Google DeepMind, arXiv 2507.22291) includes GRACE among its inputs. Its supplement reports GRACE had "no significant negative effect" on loss or reconstruction, which is not evidence of benefit. [emerging] — [arXiv:2507.22291](https://arxiv.org/pdf/2507.22291)
- EarthPT (700M-parameter autoregressive EO time-series transformer) works on optical reflectance, not gravity. [related] — [arXiv:2309.07207](https://arxiv.org/pdf/2309.07207)
- Yu, Peng, Lin, Keenan, Bishop, arXiv 2608.15282 (Aug 2026), review EO foundation models for ecohydrology. Pretraining is dominated by optical and active-microwave data. The abstract does not mention GRACE; a search snippet suggests the PDF lists GRACE as an input, which is unverified. [review] — [arXiv:2608.15282](https://arxiv.org/abs/2608.15282)

### Inferences
- The robust ML contribution is **ML as a predictor and gap-filler checked against independent gravity**: TWSTORE's SLR/DORIS + VCE fusion, and validation against streamflow and sea-level budgets. "Super-resolution" checked only against GRACE itself is not robust. A Rust toolkit can contribute evaluation tooling, not more models:
  - EOF/PCA of SH or gridded series (TWSTORE works in a GRACE-EOF basis);
  - variance component estimation for combining fields;
  - a standard benchmark set built from the 259 COST-G fields, with mandatory linear-regression/ARX baselines, as Nie et al. 2025 argue.
- Deep-learning SH destriping is an open niche. A defensible version would learn the stripe pattern from the formal VCMs or from ocean residuals, using the COST-G series as training data. It must be benchmarked against DDK, VADER/VDK, VMD-stepwise and PISP, which a Rust `gravimetry` crate could implement as baselines.
- Neural fields compete with SH only where SH converge poorly (inside the Brillouin sphere, small bodies). For Earth time-variable gravity, SH or mascons stay the natural basis. A fast SH synthesis crate is useful to ML users mainly as a fast generator of training data and as a physics layer.

### Gaps
- I could not verify any 2023–2026 paper on DL denoising of Level-2 SH coefficients or of daily Kalman fields. This may exist in Chinese-language journals or recent conference abstracts.
- I did not check the claimed "abrupt decline in global TWS since 2014–2016" (Rodell and co-authors) or related drying claims, which are contested. The search budget was exhausted.
- Quantitative skill and full methods for ML-TWiX (Sci. Data 2026) and the 2024 RSE transformer gap filler are unverified (paywalls).
- I found no gravity-specific foundation model, and no neural-field model of Earth's *time-variable* field.

---

## 2. Sub-monthly and daily gravity, high-low SST from Swarm and other LEOs, combination with altimetry and GNSS loading

### Takeaway
Sub-monthly gravity is still research-grade, not operational. ITSG-Grace2018 daily Kalman fields (d/o 40) remain the reference. New ideas since 2023 include:
- CSR 5-day mascons (latency projected at 3–5 days, available only on request);
- line-of-sight gravity differences (LGD) and K-band range-rate post-fit residuals for event detection (hurricanes, cyclones, monsoon floods);
- ML transfer of mascon detail into daily fields (D-SHIFT).

I found no operational daily product from GFZ or COST-G. On the high-low side, GOCO2025s (TU Graz, 2026) combines kinematic orbits of 18 LEOs with GRACE-FO laser interferometer (LRI) data. Joint GNSS-loading + GRACE inversions have become routine regionally (Brazil, California, China).

### Cited Findings

**Daily and 5-day solutions (status: emerging)**
- ITSG-Grace2018 daily solutions are Kalman-smoothed and estimated to d/o 40 (~500 km half-wavelength) with no missing days. Days with too few observations lean toward the background and annual signal. [established, foundational 2018–2019; **Open: data**] — [TU Graz ITSG-Grace2018](https://www.tugraz.at/institute/ifg/downloads/gravity-field-models/itsg-grace2018)
- A flood-detection study using the Kalman-smoothed daily fields (Apr 2002–Aug 2016) with extreme precipitation reports GRACE detected ~81% of Dartmouth Flood Observatory (DFO) floods. I attribute this to ESSD 15, 521 (2023) from a search summary; verify the exact figure. [emerging] — [ESSD 15, 521 (2023)](https://essd.copernicus.org/articles/15/521/2023/)
- **CSR 5-day (5D) mascons:** Rateb, Save, Sun, Scanlon, *Sci. Rep.* 14:13841 (2024), doi:10.1038/s41598-024-64491-w.
  - Built from GRACE/GRACE-FO tracking data only, on a geodesic grid of 42,107 mascons (~120 km spacing).
  - SLR C20/C30, ICE6G-D GIA; 283 missing 5D epochs filled by a Bayesian (Stan) framework.
  - Against 3,272 DFO floods (2002–2021), 80–85% of floods were preceded by elevated antecedent storage. Response coincidence was 0.35–0.71, and Katrina is visible in 5D but not monthly data.
  - Latency of 3–5 days (one day with "quick look") is projected, not demonstrated.
  - [emerging; **Open:** derived data doi:10.18738/T8/6HKCGW; 5D fields **on request only**; code **not released**] — [PMC11180122](https://pmc.ncbi.nlm.nih.gov/articles/PMC11180122/)
- **K/Ka-band range-rate post-fit residuals as line-of-sight gravity differences (LGD):** Cuadrat-Grzybowski, Teixeira da Encarnação, Visser (TU Delft), arXiv 2503.04227 (Mar 2025); JGR Solid Earth 2025, doi:10.1029/2025JB031615. CSR residuals show the Zapiola Rise dipole modes, Cyclone Oswald (2013) and the Bangladesh 2007 monsoon, and correlate highly with ITSG daily fields. The authors propose them as input to a high-frequency mass-change framework finer than Kalman methods. [emerging; residuals provided by CSR, code not stated] — [arXiv:2503.04227](https://arxiv.org/abs/2503.04227); [JGR](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2025JB031615)
- A follow-up "data-driven high-frequency mass change models from GRACE orbit residuals" exists (title only seen). [emerging] — [ResearchGate 406969352](https://www.researchgate.net/publication/406969352_Novel_Data-driven_High-Frequency_Mass_Change_Models_from_GRACE_orbit_residuals)
- Li et al., GRL 2025 (doi:10.1029/2025GL116973): an LGD-based method gives **daily** TWS changes during Hurricanes Katrina and Laura. GRACE KBR and GRACE-FO LRI perform comparably for hurricane floods. [emerging] — [GRL](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2025GL116973)
- Storm Daniel (Thessaly, Sep 2023) shows up even in **monthly** JPL GRACE-FO data as a ~10 cm residual peak in a small catchment. [emerging] — [IAG Symp. chapter 1345_2026_307](https://link.springer.com/chapter/10.1007/1345_2026_307)
- **GFZ/COST-G daily:** I found no operational daily product. GravIS (ESSD 17, 611, 2025) is monthly. GFZ operational GRACE-FO Level-2 moved from RL06.1 to **RL06.3 in June 2024**. Near-real-time LRI+SLR work (NELOS) is research-stage. [status: open] — [GravIS ESSD 2025](https://essd.copernicus.org/articles/17/611/2025/); [GSTM2024-89 GFZ L2/L3 overview](https://presentations.copernicus.org/GSTM2024/GSTM2024-89_presentation.pdf); [GFZ projects](https://www.gfz.de/en/section/global-geomonitoring-and-gravity-field/projects)
- The GRACE-FO Science Team Meeting 2025 (planned Boulder, 7–9 Oct 2025) switched to remote-only. I found no published summary. — [UCAR CPAESS](https://cpaess.ucar.edu/meetings/grace-fo-2025-science-team-meeting)

**High-low SST, multi-LEO and combined static/temporal models**
- **GOCO2025s** (Öhlinger, Mayer-Gürr, Krauss, Dumitraschkewitz, Strasser, Süsser-Rechberger, Brockmann, Niedermaier; ESSD preprint essd-2026-639, discussion since 29 Sep 2026):
  - static field to d/o 300, with constrained secular, annual and semi-annual terms to d/o 200;
  - 23 years of GOCE, GRACE and GRACE-FO data plus **kinematic orbits of 18 LEOs** and SLR to 9 geodetic satellites, combined at normal-equation level;
  - **GRACE-FO LRI dominates d/o 70–200**; processed with GROOPS.
  - [emerging; **Open: data** doi:10.3217/f48p8-8h651 including full normal equations and uncertainties] — [ESSD preprint](https://essd.copernicus.org/preprints/essd-2026-639/essd-2026-639.pdf)
- Swarm-only monthly fields are a COST-G product used to bridge the gap. SLR constrains the very low degrees and high-low SST gives higher resolution at lower precision. [established] — [J. Geod. 2024, 10.1007/s00190-024-01888-5](https://link.springer.com/article/10.1007/s00190-024-01888-5); [COST-G](https://geodesy.science/item/cost-g/)
- Sentinel-1/-2/-3 GPS kinematic positions as pseudo-observations give "mostly uninterrupted" large-scale time-variable gravity (Bern). [emerging] — [EGU22 abstract](https://ui.adsabs.harvard.edu/abs/2022EGUGA..24.5202G/abstract)
- A network precise-orbit solution of 8 LEOs (GRACE-FO, Sentinel-3A/B, Sentinel-6A, Swarm) improves absolute position from ~4 cm (3 satellites) to ~2–2.5 cm 1D RMS. [emerging, 2022] — [PMC8838239](https://pmc.ncbi.nlm.nih.gov/articles/PMC8838239/)

**Combination with GNSS loading (status: established regionally)**
- Zhong et al., WRR 2024 (doi:10.1029/2024WR037538), Brazil: GRACE-based **geopotential differences** jointly inverted with GNSS vertical displacements via Green's functions. Uncertainty is lower than GRACE-only mascons, and GRACE trends replace GNSS trends at 103 stations to suppress tectonic contamination. [established/emerging] — [WRR](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024WR037538)
- California (Remote Sensing of Environment 311:114303, 2024): a hybrid **GNSS + GRACE + InSAR** inversion that separates groundwater loss in the Central Valley. [emerging] — [VTechWorks PDF](https://vtechworks.lib.vt.edu/server/api/core/bitstreams/ca880660-8e8e-4dda-a3d6-40495e48dbb3/content)
- Guangdong (GJI 2026, ggag401): 49 GNSS stations plus GRACE/GFO mascons (2017–2024) with objective-space weighting to reduce boundary artefacts. [emerging] — [GJI ggag401](https://academic.oup.com/gji/advance-article/doi/10.1093/gji/ggag401/8853772)
- Yellow River (Chin. J. Geophys. 2025): GRACE-derived "virtual stations" fill sparse GNSS networks. [emerging] — [CJG 2025](http://www.geophy.cn/en/article/doi/10.6038/cjg2025T0202)
- Carlson et al. 2022: the GNSS+GRACE framework with wavelet separation. [foundational, pre-2023] — [PMC9287077](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9287077/)

### Inferences
- The high-frequency frontier is moving **below Level-2** (KBRR/LRI residuals, LGD). A toolkit that works only on SH Level-2 fields can still help by:
  - reading ITSG daily fields and CSR 5D grids;
  - computing LGD-like along-track gravity differences *from* SH fields along simulated orbits (the `satgrav` closed-loop simulator already has the orbit machinery);
  - reproducing the event-detection statistics (flood coincidence rates).
- GNSS loading is a natural next crate feature. `gravimetry` already has Love numbers, so predicted vertical/horizontal elastic loading displacements (load Love numbers h′, l′) from the COST-G fields at GNSS station coordinates are a small step. They would enable the GRACE+GNSS joint inversions above.
- GOCO2025s ships full normal equations, which allows combination and error-propagation experiments without reprocessing. A Rust reader for its time-variable ICGEM format (secular/annual/semi-annual terms) fits the `gravity-field` crate.

### Gaps
- I could not confirm whether a newer ITSG release (e.g., an "ITSG-Grace2024" or operational daily stream) exists; the search budget was exhausted.
- I could not confirm whether COST-G has issued a new release (e.g., RL02) or a daily or near-real-time combination.
- Combination with **altimetry** for ocean mass/sea-level budgets after 2023 was not covered beyond the altimetry-regularized Slepian method (section 3).

---

## 3. Regularization and estimation advances (mascons, Slepian/regional, Bayesian/GP, uncertainty quantification, leakage, GIA)

### Takeaway
The clearest 2024–2026 advance in estimation is **stochastic modelling of background-model errors**. The AOe07 error series for AOD1B RL07, time-variable error variance-covariance matrices (VCMs) in closed-loop simulation, and co-estimation of AOD coefficients with VCM constraints in GFZ Level-2 processing move aliasing treatment from "remove the model" to "weight by the model's uncertainty". Slepian methods gained altimetry-informed regularization: Greenland accuracy improved by 25–58% over classical Slepian. Bayesian methods appear for gap filling and reconstruction, but I found **no Gaussian-process GRACE paper**. 3-D (laterally heterogeneous) GIA now has quantified GRACE impacts: the Antarctic loss rate changes by roughly 12 Gt/yr. I found no ICE-7G result.

### Cited Findings

**Background-model uncertainty and aliasing (status: emerging, entering operations at GFZ)**
- **AOe07** (Shihora et al., J. Geod. 2024, doi:10.1007/s00190-024-01832-7) is a 6-hourly error series for **AOD1B RL07**, 1995–2020, built as the difference between atmosphere-ocean models. It is intended both for monthly processing and for simulations of future missions. [emerging/established for simulation; Open: AOe07 distributed by GFZ] — [J. Geod. 2024](https://link.springer.com/article/10.1007/s00190-024-01832-7)
- **Abrykosov & Pail, J. Geod. 99:54 (2025), doi:10.1007/s00190-025-01980-4**:
  - derive stationary and *non-stationary* AO error VCMs from AOe07 and propagate them in GRACE-type closed-loop simulations;
  - AO uncertainty varies in time as well as space, so stationary weighting alone is insufficient;
  - in a Bender double pair, the inclined pair's intrinsic de-aliasing keeps long-wavelength errors in their own band, so weighting can target them.
  - Numbers not verified (paywall). [emerging] — [J. Geod.](https://link.springer.com/article/10.1007/s00190-025-01980-4); [GFZ GRACE references](https://www-app2.gfz.de/grace-references/b2hd-2025JGeod..99...54A.html)
- Murböck, Dahle, Panafidina, Hauk, Wilms, Neumayer, Flechtner (EGU25-16562): AOD stochastic modelling is **implemented in GFZ GRACE/GRACE-FO Level-2 processing**. AOD coefficients are co-estimated using error VCMs as constraints, and different temporal-correlation lengths are tested. No numbers in the abstract. [emerging→operational] — [EGU25-16562](https://meetingorganizer.copernicus.org/EGU25/EGU25-16562.html)
- Kvas & Mayer-Gürr (2019) used stationary AO uncertainty in real processing (ITSG). [foundational, pre-2023] — cited in [GOCO2025s preprint](https://essd.copernicus.org/preprints/essd-2026-639/essd-2026-639.pdf) and [EGU summary of Abrykosov & Pail](https://link.springer.com/article/10.1007/s00190-025-01980-4)

**Slepian, regional and leakage methods**
- **Improved Slepian method** (Gong, Ran, Han et al., J. Geod. 100:3, Dec 2025, doi:10.1007/s00190-025-02022-9):
  - satellite-altitude pseudo-observations separate the Greenland ice sheet from ocean and island signals, with **altimetry as the regularization matrix**;
  - accuracy improves by **25–58%** over classical Slepian (validated against the input–output method);
  - classical Slepian underestimates the trend by ~80 Gt/yr; the improved method is comparable to mascons.
  - [emerging] — [J. Geod.](https://link.springer.com/article/10.1007/s00190-025-02022-9); [ADS](https://ui.adsabs.harvard.edu/abs/2025JGeod.100....3G)
- A Slepian approach for mass-term sea level in the South China Sea with buffer-zone leakage correction (Int. J. Appl. Earth Obs. Geoinf. 2024). [emerging] — [S1569843224004199](https://www.sciencedirect.com/science/article/pii/S1569843224004199)
- Slepian inversion directly from GRACE inter-satellite geopotential differences (Chin. J. Geophys. 2023). [emerging] — [CJG 2023](http://en.dzkx.org/article/doi/10.6038/cjg2023Q0970)
- Mascon regularization and error characterization (GSFC). [foundational, pre-2023] — [NTRS 20190004972](https://ntrs.nasa.gov/citations/20190004972)

**Bayesian, Gaussian-process and uncertainty quantification**
- A Bayesian (MCMC) additive reconstruction of GRACE series: trend, interannual, annual and semi-annual terms. It captures 99% of variability at basin scale and 78% at 1° grid, and is proposed for nowcasting. [foundational 2022; emerging] — [Rateb et al. 2022, ESS, doi:10.1029/2021EA002162](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2021EA002162)
- The CSR 5D mascons fill gaps with a Bayesian framework (Stan). [emerging] — [Sci. Rep. 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11180122/)
- The ML-TWiX deep ensemble outputs TWSA uncertainty; TWSTORE distributes uncertainties from variance component estimation. [emerging] — [Sci. Data 2026](https://www.nature.com/articles/s41597-026-06604-w); [ESSD 2026](https://essd.copernicus.org/articles/18/1747/2026/)
- I found no 2023–2026 paper using a **Gaussian process** prior on GRACE mass change. The closest hierarchical-GP analogue is climate-model bias correction (GMD 2024). [open problem] — [GMD 17, 5733](https://gmd.copernicus.org/articles/17/5733/2024/)
- GOCO2025s notes that observation-level stochastic modelling plus background-model uncertainties yields "realistic formal errors". [established practice at TUG] — [GOCO2025s](https://essd.copernicus.org/preprints/essd-2026-639/essd-2026-639.pdf)

**GIA advances (status: emerging)**
- Hightower et al., JGR Solid Earth 2026 (doi:10.1029/2026JB033791) use **3-D mantle viscosity from seismic tomography** for the GRACE(-FO) GIA correction. Global-mean ocean and Greenland trends barely change, but regional patterns change and the Antarctic loss estimate falls. A JPL postdoc poster quantifies the Antarctic change at **~12 Gt/yr**. [emerging] — [JGR](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2026JB033791); [JPL poster](https://postdocs.jpl.nasa.gov/media/2025%20RPD/2025%20PRD-E-025.pdf)
- 2024 (J. Geodynamics): Antarctic GIA mass-change estimates show **limited sensitivity** to lateral viscosity variations; only two earlier studies had used 3-D models for Antarctic GIA corrections. **This partly contradicts the 2026 result**. [contested] — [S0264370724000309](https://www.sciencedirect.com/science/article/abs/pii/S0264370724000309)
- A coupled ice-sheet–solid-Earth model (The Cryosphere 18, 4233, 2024): weak West Antarctic structure is hard to reconcile with Northern Hemisphere viscosities. [emerging] — [TC 2024](https://tc.copernicus.org/articles/18/4233/2024/)
- van der Wal et al. 2015: 3-D composite viscosity gives smaller Antarctic loss than 1-D. [foundational, pre-2023] — [EPSL 2015](https://www.sciencedirect.com/science/article/pii/S0012821X15000059)
- The CSR 5D product still uses **ICE6G-D**. [established practice] — [Sci. Rep. 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11180122/)

### Inferences
- For `satgrav`, the most valuable addition is a **background-model error layer**: AOe07-style stochastic AO errors injected into the closed-loop simulation, plus optional VCM-weighted estimation. The user's simulator currently has no temporal aliasing, so it cannot reproduce today's dominant error source or the single- vs double-pair difference. A first step would read AOe07 SH time series and add them as "truth minus model" perturbations.
- Slepian functions (concentration problem on a spherical cap or region) and regional leakage correction are a clean, reusable Rust contribution. I found no maintained Rust or C implementation; most code is MATLAB (Simons) or Python. Altimetry-regularized or pseudo-observation variants could build on that.
- GIA is moving from single models to ensembles and 3-D Earth structure. A toolkit should let users choose or ensemble GIA SH rate fields (ICE6G-D, 3-D variants) and report the spread, not hard-code one.

### Gaps
- I could not check 2023–2026 updates to JPL, CSR or GSFC mascon releases (e.g., RL06.3 mascons, RL07 plans) or new mascon regularization variants; the search budget was exhausted.
- I found no ICE-7G GRACE-correction result; ICE-7G status is unknown.
- The quantitative gains from Abrykosov & Pail (2025) and Murböck et al. (2025) are unverified.
- An "FES2022 + AOD1B RL07 reduce noise over oceans/deserts" result appeared in a search summary without a clear primary source; not included as a finding.

---

## 4. Next-generation missions, temporal-aliasing mitigation and simulation software

### Takeaway
- **GRACE-C** is in development: a NASA OIG audit (Sep 2026) says it is on track to launch before **July 2029**, within its $658M baseline, but has burned 61% of its reserves. A GFZ page still says Dec 2028.
- **NGGM** (ESA inclined pair, 397 km / 70°): final design (Phase B2/C) was funded at ESA's Nov 2025 Ministerial (CM25), launch ~2032, forming MAGIC with GRACE-C.
- **Quantum:** CARIOQA Phase B started Oct 2025 (24 months, target TRL 6). NASA's QGGPf got a further $20M in Aug 2026 (launch ~2030, technology demo). A new Horizon Europe study, **FEATS-G** (GFZ-led, Sep 2026–Aug 2029), will run end-to-end quantum-gravimetry simulations.
- **Aliasing mitigation** has converged on Bender double pairs, Wiese-type daily low-degree co-estimation, ocean-tide error stochastic modelling, AOe07 VCM weighting and VADER filtering.
- **Open software:** GROOPS (GPL v3, C++) is the only full open gravity-recovery and simulation suite I found.

### Cited Findings

**GRACE-C (status: in development)**
- NASA OIG audit (Sep 2026):
  - on track to launch **before July 2029** within a **$658M** baseline (excluding DLR);
  - $50.5M of $82.5M project reserves spent as of May 2026 (61.2%), with $13.6M unallocated;
  - a Feb 2026 JPL independent estimate projects a $9.3M overrun against $628M, and euro exposure adds ~$8M;
  - Tesat (lasers), ONERA (accelerometers), DLR mission control.
  - [established] — [NASA OIG](https://oig.nasa.gov/news/audit-of-nasas-upcoming-twin-satellite-mission-identifies-cost-and-schedule-risks/); [oversight.gov report](https://www.oversight.gov/reports/audit/nasas-management-gravity-recovery-and-climate-experiment-continuity-grace-c-mission); [Exterra summary](https://www.exterrajsc.com/p/audit-finds-grace-c-gravity-mission)
- GFZ lists GRACE-C launch as **Dec 2028**, which conflicts with OIG's "before July 2029" commitment. Treat Dec 2028 as the project target and July 2029 as the agency commitment. — [GFZ S-NGGM](https://www.gfz.de/en/section/global-geomonitoring-and-gravity-field/projects/s-nggm-studies-for-next-generation-gravity-missions)

**NGGM / MAGIC (status: in development)**
- At LPS25 (Jun 2025), NGGM was at the end of Phase B1, with the system requirements review (ISRR) planned for Sep 2025. Orbit ~397 km, 70° inclination; MAGIC targets ≥4 years of combined operations from 2032. [established] — [LPS25 status (Daras)](https://lps25.esa.int/lps25-presentations/presentations/806/848_Daras_update.pdf); [LPS25 industry status](https://lps25.esa.int/lps25-presentations/presentations/2516/_2516.pdf)
- **CM25 (Bremen, Nov 2025, €22.3bn total)**: ESA's post-CM25 FutureEO deck lists "Final design (phase B2/C) of the New Generation Gravity Mission 'NGGM'" as secured. It also notes "NGGM — Kick-off Bridging phase held. ADS GmbH and TAS-I working jointly". [established] — [ESA FutureEO post-CM25 (Floberghagen)](https://climate.esa.int/documents/3235/1_FutureEO_post-CM25_-_Floberghagen.pdf); [ESA CM25 news](https://www.esa.int/About_Us/Corporate_news/ESA_Member_States_commit_to_largest_contributions_at_Ministerial)
- eoHandbook lists NGGM launch Dec 2032 and end of life Dec 2040 (provisional). — [eoHandbook](https://database.eohandbook.com/database/missionsummary.aspx?missionID=1386)
- An ESA story frames NGGM's science case partly around detecting **AMOC weakening**. — [ESA](https://www.esa.int/Applications/Observing_the_Earth/FutureEO/New_ESA_gravity_mission_to_detect_weakening_ocean_conveyor)
- **Science-impact simulation:** Liaqat, Camici, Leopardi, Gaona, Brocca, HESS 30, 4969 (2026), estimate **precipitation from TWS** (SM2RAIN inversion).
  - With NGGM/MAGIC-like errors (<10 mm, 5-day sampling) the mean correlation is R > 0.78, against R = 0.40 for GRACE-FO-like (25 mm, 30-day).
  - Satisfactory performance (R ≥ 0.7) needs 5-day sampling with errors ≲ 15 mm.
  - The study adds Gaussian noise of 1.9, 4.2, 19 and 42 mm to a proxy; no real GRACE data are used.
  - [emerging; data on request] — [HESS 2026](https://hess.copernicus.org/articles/30/4969/2026/)

**Aliasing mitigation (status: established in simulation; partly operational)**
- From the ESA MAGIC/Science results page (TUM):
  - altitude is the main performance driver; an inclined pair below ~400 km reduces retrieval errors;
  - ocean-tide error **stochastic modelling** (Abrykosov et al. 2021) helps both 3- and 7-day double-pair retrievals;
  - **classic Wiese co-parameterization** (daily low-degree fields, Wiese et al. 2011) slightly beats the data-driven multistep self-de-aliasing method (DMD, Abrykosov et al. 2022a), and both beat no co-parameterization;
  - co-estimating trends over 10 years gives no benefit;
  - VDK beats DDK filtering.
  - [established in simulation] — [TUM IAPG MAGIC](https://www.asg.ed.tum.de/iapg/magic/results-4)
- Iterating ocean-tide co-estimation (feeding the estimated tide model back as de-aliasing) cuts ocean-tide aliasing by >50% for a Bender-type NGGM. [foundational, 2018] — [Hauk & Pail 2018 (ResearchGate)](https://www.researchgate.net/publication/325699815_Treatment_of_ocean_tide_aliasing_in_the_context_of_a_next_generation_gravity_field_mission)
- Bender "self-de-aliasing" lets short-period AO signals be estimated rather than only removed. [foundational, 2017] — [Daras & Pail 2017, JGR, doi:10.1002/2017JB014250](https://agupubs.onlinelibrary.wiley.com/doi/10.1002/2017JB014250)
- AOe07 VCM weighting and GFZ operational AOD co-estimation are covered in section 3. — [J. Geod. 2025](https://link.springer.com/article/10.1007/s00190-025-01980-4); [EGU25-16562](https://meetingorganizer.copernicus.org/EGU25/EGU25-16562.html)
- MAGIC retrieval optimized with background-model uncertainty information (2025). [emerging] — [ResearchGate 389736469](https://www.researchgate.net/publication/389736469_Optimized_gravity_field_retrieval_for_the_MAGIC_mission_concept_using_background_model_uncertainty_information)
- Enhanced electrostatic and optical accelerometry for future missions. [emerging, 2023] — [arXiv:2310.14875](https://arxiv.org/html/2310.14875)

**Quantum pathfinders and studies (new since the previous report)**
- **CARIOQA:**
  - Phase B began **Oct 2025** for 24 months, aiming for TRL 6 of payload and platform critical technologies (CNES + DLR coordination). Phase A (early 2024–mid 2025) found a pathfinder "technically feasible within the next decade". No performance numbers or new simulations were reported (Biskupek et al., EGU26-9366). [established status] — [EGU26-9366](https://meetingorganizer.copernicus.org/EGU26/EGU26-9366.html); [GMV](https://www.gmv.com/en-es/communication/news/carioqa-phase-b-launches-europe-strengthens-leadership-quantum-sensing-space)
  - Exail delivered a space-grade laser source to CNES. — [Exail](https://www.exail.com/news/exail-delivers-space-grade-laser-source-to-cnes-for-europe-carioqa-quantum-pathfinder-mission)
  - **Phase A orbit and accelerometer data are public on Zenodo.** [Open: data] — [Zenodo 17580678](https://zenodo.org/records/17580678)
  - A CARIOQA space-weather/thermosphere application paper (IAG Symposia 2026). — [Springer 1345_2026_314](https://link.springer.com/chapter/10.1007/1345_2026_314)
- **NASA QGGPf:** Infleqtion got a **$20M follow-on** (press release 25 Aug 2026), bringing NASA funding to $40M, for hardware integration, subsystem fabrication and microgravity ground testing. Launch target ~2030. Monarch Quantum supplies photonic light engines (Feb 2026). It is a technology demo, "not a finished science instrument". These figures are company press releases, not NASA statements. [status] — [Infleqtion](https://infleqtion.com/nasa-awards-infleqtion-20-million-as-worlds-first-quantum-gravity-mission-advances-toward-flight/); [PR Newswire, Monarch](https://www.prnewswire.com/news-releases/monarch-quantum-selected-to-deliver-integrated-photonics-quantum-light-engines-for-nasas-first-planned-space-based-quantum-gravity-gradiometer-302697215.html); [Quantum Insider 2026](https://thequantuminsider.com/2026/02/10/infleqtion-and-nasa-to-fly-the-worlds-first-quantum-gravity-sensor-to-space/)
- **FEATS-G** ("Quantifying the potential of quantum space gravimetry"):
  - Horizon Europe grant 101295973, **€500k**, **1 Sep 2026–31 Aug 2029**, GFZ (contacts H. Dobslaw, R. Schachtschneider);
  - end-to-end simulations of quantum gravimetry *and* gradiometry concepts using CARIOQA and industry sensor specifications, benchmarked against GRACE/GRACE-FO, GRACE-C and NGGM/MAGIC;
  - partner count conflicts: GFZ says 13 European + 5 Japanese institutions (17 listed), the search snippet of CORDIS says 13 + 7;
  - software and data release not stated.
  - [emerging] — [GFZ FEATS-G](https://www.gfz.de/en/section/earth-system-modelling/projects/feats-g-quantifying-the-potential-of-quantum-space-gravimetry); [CORDIS 101295973](https://cordis.europa.eu/project/id/101295973)
- "Quantum gravity gradiometry for future mass change science", *EPJ Quantum Technology* (2025, doi:10.1140/epjqt/s40507-025-00338-1), describes a pathfinder to validate atom-interferometric gravity remote sensing. Details not read (paywall redirect). [emerging] — [EPJ QT](https://link.springer.com/article/10.1140/epjqt/s40507-025-00338-1)

**Simulation and processing software**
- **GROOPS** (TU Graz) is C++ under **GPL v3** on GitHub. It does gravity recovery from satellite and terrestrial data, GNSS orbit determination, SLR and GNSS constellations, with MPI parallelism. It was used for ITSG-Grace2018 and GOCO2025s. [established; **Open: code**] — [GitHub groops-devs/groops](https://github.com/groops-devs/groops); [GSTM2024-55](https://meetingorganizer.copernicus.org/GSTM2024/GSTM2024-55.html)
- End-to-end simulators for GRACE-C/NGGM exist (Daras et al. 2024, cited by GOCO2025s), but I found no open-source release. — [GOCO2025s intro](https://essd.copernicus.org/preprints/essd-2026-639/essd-2026-639.pdf)

### Inferences
- MAGIC aliasing performance depends on an established toolbox: Wiese daily low-degree co-estimation, ocean-tide/AO error stochastic models and VADER. Only GROOPS implements it openly, as a monolithic C++/GPL code. A permissively licensed (MIT/Apache) Rust library offering composable pieces (Bender geometry, Wiese parameterization, AOe07 injection, VCM weighting, DDK/VADER) would fill a real gap for FEATS-G-style studies and for students.
- The HESS 2026 study shows mission-impact work is often done with simple Gaussian-noise proxies. `satgrav` could provide a more realistic but still lightweight error emulator: degree-dependent, anisotropic and correlated noise per mission (GRACE-FO, GRACE-C, NGGM, MAGIC, CAI-hybrid). It would plug straight into such hydrological studies.
- The release of CARIOQA Phase A orbit and accelerometer data makes a public validation case for a Rust hybrid CAI + electrostatic accelerometer noise model.

### Gaps
- I could not open the CM25 subscription tables for MAGIC-specific funding, or the September 2025 NGGM ISRR outcome.
- I found no 2025–2026 quantitative CAI-mission simulation results beyond those in the previous report; FEATS-G has only just started.
- I could not access the EPJ QT 2025 quantum-gradiometry paper's numbers.

---

## 5. Chronometric and relativistic geodesy with optical clocks, and the link to gravity-field models (geo5d)

### Takeaway
Chronometric levelling has moved from laboratory demos to long-baseline, multi-country campaigns:
- A **457 km fibre-linked Sr-clock comparison** measured ΔW = 3918.1(2.6) m²/s², about 27 cm in height (PR Applied 2024).
- A **10-clock, 6-country, 45-day** coordinated comparison over fibre and satellite links (Optica 2025).
- Transportable Sr clocks below 3×10⁻¹⁸ (TerraQ), where 10⁻¹⁸ corresponds to ~1 cm.

**ACES** launched to the ISS on 21 Apr 2025 and was switched on 28 Apr 2025. A Sep 2026 simulation finds its optical link could deliver decimetre-level height differences in days and centimetre level over longer periods; no flight geopotential result is published yet. A 2026 IAG paper explicitly studies combining clock networks with global geopotential models and their commission errors. That is exactly geo5d's chronometric-levelling use case.

### Cited Findings
- **Long-distance chronometric levelling** (PR Applied 21, L061001, 3 Jun 2024; arXiv 2309.14953): two Sr optical lattice clocks over a 457 km interferometric fibre give ΔW = **3918.1(2.6) m²s⁻²**, against geodetic **3915.88(0.30)**, i.e. ~27 cm height uncertainty. [established demonstration] — [PR Applied](https://journals.aps.org/prapplied/abstract/10.1103/PhysRevApplied.21.L061001); [arXiv:2309.14953](https://arxiv.org/abs/2309.14953)
- Chronometric levelling with transportable optical clocks beyond the laser coherence limit (2023). [emerging] — [arXiv:2310.08835](https://arxiv.org/pdf/2310.08835)
- **Lindvall et al., Optica 12, 843 (2025), doi:10.1364/OPTICA.561754**: **ten optical clocks in six countries** (FI, FR, DE, IT, UK, JP) measured simultaneously over fibre and satellite links for 45 days. It reports 38 frequency ratios, four measured directly for the first time. Satellite links carry uncertainties "up to 100 times larger" than the clocks, which motivates fibre. The redefinition of the second is targeted for 2030. [established] — [arXiv:2505.06763](https://arxiv.org/abs/2505.06763); [Optica OPN news](https://www.optica-opn.org/home/newsroom/2025/june/optical_clock_network_spans_six_countries/); [Observatoire de Paris](https://lte.observatoiredeparis.psl.eu/Global-Optical-Clock-Network-Redefines-Precision-Timing)
- An international optical clock comparison over the European fibre network (Hg-199 at LNE-OP against Yb⁺, Yb and Sr) reports improved ratio uncertainties. [emerging] — [arXiv:2604.27963](https://arxiv.org/pdf/2604.27963); [PR Research](https://journals.aps.org/prresearch/abstract/10.1103/l4bh-ryxs)
- **TerraQ**: transportable Sr lattice clocks below **3×10⁻¹⁸**. Clocks in the 10⁻¹⁸ regime are needed for cm-level heights; PTB lattice clocks are being upgraded to low 10⁻¹⁸. [emerging] — [arXiv:2503.09272](https://arxiv.org/pdf/2503.09272)
- **ACES**: launched **21 Apr 2025** (SpaceX CRS-32), mounted on Columbus 25 Apr, first switched on **28 Apr 2025**. It is designed to measure continental geopotential differences between ground clocks. [established] — [ESA ACES](https://www.esa.int/Science_Exploration/Human_and_Robotic_Exploration/ACES_Atomic_Clock_Ensemble_in_Space); [arXiv:2411.02912](https://arxiv.org/pdf/2411.02912)
- **ACES height-difference simulation**: Lachmann, Müller, Vollmair, Schlicht (LUH/TUM), arXiv 2609.21786 (18 Sep 2026, submitted to *Remote Sensing*).
  - Full simulation of the microwave (MWL) and optical (ELT) links with clock and link noise models;
  - a **slope-based estimator** extracts the redshift;
  - quasi-common-view, non-common-view and split non-common-view schemes;
  - ELT reaches **dm in a few days and cm over longer periods**; MWL is worse;
  - operational optical data are not yet available.
  - [emerging; code not stated] — [arXiv:2609.21786](https://arxiv.org/abs/2609.21786)
- "How accurately can an optical clock network estimate absolute gravity potential when combined with global geopotential models?" (IAG Symposia 2026, Springer 1345_2026_308) studies propagation of GGM commission error into clock-derived absolute potential. Authors and numbers not verified. [emerging] — [Springer](https://link.springer.com/chapter/10.1007/1345_2026_308)
- Acta Physica Sinica review (2026) of optical-clock precision measurement, including chronometric levelling against classical geodesy. [review] — [doi:10.7498/aps.75.20251594](https://wulixb.iphy.ac.cn/en/article/doi/10.7498/aps.75.20251594)
- The concept goes back to Bjerhammar (1985). [foundational] — [arXiv:2503.09272](https://arxiv.org/pdf/2503.09272); see also [Atomic Clocks for Geodesy, arXiv:1803.01585](https://arxiv.org/pdf/1803.01585)

### Inferences
- geo5d's chronometric levelling could add three things:
  - **forward-modelling ΔW between clock sites from a GGM** (`geopotential` synthesis plus the normal potential), with commission error from the coefficient sigmas and omission error from a Kaula or degree-variance extrapolation, as in the 2026 IAG chapter;
  - **time-variable potential corrections** from the COST-G monthly fields: hydrology changes W by mm-level equivalent heights, comparable to 10⁻¹⁹ fractional frequency, so it matters below 10⁻¹⁸;
  - tidal and permanent-tide conventions, using the existing Love numbers.
  This would make the crates a reference implementation for comparing clock and geodetic ΔW, as in the 457 km test (agreement of 2.2 m²/s² within 2.6 m²/s² uncertainty).
- The ACES simulation's slope estimator and its common-view/non-common-view schemes are simple enough to port as a small module. It would connect `satgrav`-style orbit simulation (ISS orbit) with geo5d redshift.

### Gaps
- I found no ACES flight geopotential results as of Oct 2026, nor I-SOC status updates (search budget exhausted).
- I found no 2023–2026 result on clock-based **height-system unification** (IHRS/IHRF) specifically.
- I did not verify the authors and numbers of the IAG 2026 clock+GGM chapter.

---

## 6. Planetary and lunar gravimetry with reusable methods

### Takeaway
The most significant new planetary result is GRAIL's **degree-3 tidal Love number** k₃ = 0.0163 ± 0.0007, 72% above a spherically symmetric Moon. It implies a 2–3% nearside/farside difference in mantle shear modulus (Park et al., *Nature* 2025). Quantum gravimetry for planets is concept-stage: MaQuIs for Mars (GRACE-like pair, cold-atom interferometer for drift correction). I found no confirmed lunar quantum-gravimeter mission. Small-body gravity is where neural and mascon representations (PINN-GM-III, GeodesyNets, MasconCube) are maturing. These are the reusable methods.

### Cited Findings
- **Park et al., Nature 2025, doi:10.1038/s41586-025-08949-5**, "Thermal asymmetry in the Moon's mantle inferred from monthly tidal response":
  - **k₃ = 0.0163 ± 0.0007**, ~72% above the value for a spherically symmetric Moon;
  - explained by a 2–3% nearside/farside variation in mantle shear modulus, i.e. a **100–200 K** nearside thermal anomaly.
  - [established result; GRAIL data public via PDS] — [Nature](https://www.nature.com/articles/s41586-025-08949-5); [PMC12119328](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12119328/); [Physics Today](https://physicstoday.aip.org/news/the-moons-interior-is-surprisingly-irregular)
- GRAIL-era k₂ = 0.02416 ± 0.00022 (R = 1738 km), Williams et al. 2014. [foundational] — [JGR Planets 2014](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2013JE004559)
- **MaQuIs** (Mars Quantum Gravity Mission; *Planetary and Space Science* 239, 105800, Oct 2023): a GRACE/GRAIL-like pair ~200 km apart with a laser interferometer (1064–1550 nm). Cold-atom interferometers give drift-free correction of the inertial measurement unit and are the least mature item. Simulations suggest unprecedented resolution of Mars' static and time-variable field. [concept] — [DLR elib PDF](https://elib.dlr.de/198790/1/MaQuIs___Mars_Quantum_Gravity_Mission_in_review.pdf); [TU Delft portal](https://research.tudelft.nl/en/publications/maquisconcept-for-a-mars-quantum-gravity-mission/); [Open Access Government 2025](https://www.openaccessgovernment.org/article/maquis-pioneering-quantum-space-exploration-to-unlock-mars-interior-and-atmosphere/201817/)
- Lunar quantum sensing: a Q-CTRL/Fleet Space proposal for quantum gravimetry and magnetometry from a lunar-orbit nanosatellite (planned for 2023) has unclear status. [unverified] — [AZoQuantum](https://www.azoquantum.com/News.aspx?newsID=7842)
- Lunar gravitational-wave detector concepts (LGWA, LSGA) reuse lunar seismometry and gravimetry. [concept, 2024 review] — [Phil. Trans. A 382:20230066](https://royalsocietypublishing.org/rsta/article/382/2271/20230066/112635/Opportunities-and-limits-of-lunar-gravitational)
- Small-body methods: PINN-GM-III [emerging; code open] — [arXiv:2312.10257](https://arxiv.org/abs/2312.10257); MasconCube (~40× faster than GeodesyNets) [emerging] — [arXiv:2509.08607](https://arxiv.org/abs/2509.08607); PINNs for asteroid Eros [emerging] — [ResearchGate 383436930](https://www.researchgate.net/publication/383436930_Physics-Informed_Neural_Networks_PINNs_for_gravity_field_modelling_and_representation_With_application_to_asteroid_EROS)
- Cold-atom gravimetry for planetary missions (2020 assessment). [foundational, pre-2023] — [PSS 2020](https://www.sciencedirect.com/science/article/abs/pii/S0032063320303238)

### Inferences
- What carries over from planetary work to the Rust crates:
  - **body-agnostic SH synthesis** (GM, reference radius, normalization), so lunar (GRGM1200-class) and Mars fields read and synthesize the same way as Earth fields;
  - **degree-n tidal Love-number forward models** (k₂, k₃) for synchronous satellites;
  - **inter-satellite LGD/range-rate simulation** shared by GRACE, GRAIL and MaQuIs geometries; `satgrav` already simulates GRACE pairs, so a "body" parameter makes it a GRAIL/MaQuIs simulator almost for free;
  - a mascon/point-mass representation (MasconCube-style) as an alternative basis for regional or small-body fields.

### Gaps
- I found no 2024–2026 lunar quantum-gravimetry mission study (ESA or NASA) or confirmed flight plan.
- I did not check BepiColombo, Juno or VERITAS gravity-method papers from 2023–2026 (budget exhausted).
- The k₂ value in Park et al. 2025 was not visible in the snippets I read.
