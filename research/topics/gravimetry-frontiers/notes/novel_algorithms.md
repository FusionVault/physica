# Novel algorithms for satellite-gravimetry time series and spherical-harmonic fields

Scope: algorithms from signal processing, dynamical systems, statistics, ML and numerical analysis that could run on our 259 COST-G monthly fields (2002–2026, degree/order 90, 8,281 real coefficients per month) inside the Rust crates `geopotential`, `gravity-field`, `gravimetry` and `satgrav`. For each method the notes give what it does, prior use in gravimetry (with the search terms used), why it could help, data and compute needs, and a validation sketch.

Facts about our own data that shape every proposal:
- The COST-G monthly files are weighted combinations of about 11 centre solutions "combined by VCE on solution level". They carry `errors empirical`, meaning one calibrated sigma per coefficient, with no full covariance matrix. Licence: CC BY 4.0 (DOI 10.5880/COST-G.ICGEM_02_L2). Sources: the local file header `GSM-2_2005274-2005304_GRAC_COSTG_BF01_0201.gfc` and the [COST-G Jul 2026 .gfc](https://icgem.gfz.de/getseries/02_COST-G_/Grace-Grace-FO_RL02.1/GSM-2_2026182-2026212_GRFO_COSTG_BF01_0201.gfc) cited in our earlier report.
- So any "optimal filter with calibrated covariances" built on COST-G alone is **diagonal**. Full monthly covariances do exist for ITSG-Grace2018: the VADER filter study used the full monthly covariance matrices of ITSG-Grace2014 ([Horvath et al. 2018, GFZ PDF](https://gfzpublic.gfz-potsdam.de/pubman/item/item_3766971_3/component/file_3806902/3766971.pdf)).
- The crates already provide SH synthesis with first and second partials, degree variances, Love numbers, Gaussian smoothing, destriping, anomaly computation, TN-13/TN-14 readers, and FFT colouring and whitening. They do **not** yet provide SH analysis (grid to coefficients, quadrature) or SH rotation. Several proposals below need one or both.

Method note: 59 web searches (one general web index, not Google Scholar or Web of Science) plus about 15 page fetches. The shared search budget for this session (200 searches, shared by all agents) ran out before the numerics and some secondary checks were finished. "None found" therefore means "none surfaced by the stated queries". It is not a literature-review guarantee.

---

## 1. Space–time decompositions: EOF/PCA, rotated EOF, complex/Hilbert EOF, MSSA, POP, DMD/Koopman, SFA, ICA. Which have been applied to GRACE?

### Takeaway
EOF/PCA, rotated and extended EOF, ICA (including complex ICA), MSSA (including rotated MSSA), complex EOF (for ocean bottom pressure) and DMD have all been applied to GRACE, so none of them is novel as such. No GRACE application was found for POP/LIM, Koopman/kernel analog forecasting, or slow feature analysis. The novel openings are those three, plus DMD variants that natively handle irregular sampling and give uncertainties (optimized DMD and BOP-DMD), run directly on sigma-whitened SH coefficients.

### Cited Findings
- **EOF/PCA is established.** Wouters & Schrama (2007, GRL) applied EOF filtering to GRACE SH solutions ([UU portal](https://research-portal.uu.nl/en/publications/improved-accuracy-of-grace-gravity-solutions-through-empirical-or/)). Seo et al. (2017, JGR Solid Earth) used extended EOFs (EEOF) to separate Greenland ice-mass signal from meridional stripe noise, with less attenuation than conventional filters ([KOPRI repository](https://repository.kopri.re.kr/handle/201206/6054)). Rangelova & Sideris compared PCA/EOF with MSSA over North America ([PDF](https://topo.auth.gr/wp-content/uploads/2024/10/1-15-Rangelova.pdf)).
- **Cyclostationary EOF** has been applied to GRACE/GRACE-FO since 2003 to separate trends from internal variability (Buzzanga et al. 2025, *Commun. Earth Environ.*) ([GFZ GRACE bibliography](https://www-app2.gfz.de/grace-references/b2hd-2025ComEE...6..246B.html)). An EOF basis from GRACE is used to reconstruct TWS back to 1984 (TWSTORE, ESSD 2026) ([ESSD](https://essd.copernicus.org/articles/18/1747/2026/essd-18-1747-2026.html)).
- **Caution on EOFs:** only a few modes carry clear physical meaning, and one process can spread over several modes. This is the usual motivation for rotation ([TC preprint supplement](https://tc.copernicus.org/preprints/tc-2019-211/tc-2019-211-supplement.pdf)).
- **ICA is established.** Forootan & Kusche (J. Geod. 2012) compared PCA, rotated PCA and ICA on GRACE TWS and found that ICA extracts the original signals better from a GRACE-type superposition ([Cardiff ORCA](https://orca.cardiff.ac.uk/id/eprint/94862/); [EGU2011 abstract](https://meetingorganizer.copernicus.org/EGU2011/EGU2011-7246.pdf)). A complex ICA (CICA) for non-stationary patterns followed in *Surv. Geophys.* 2018 ([GFZ bibliography](https://www-app2.gfz-potsdam.de/grace-references/b2hd-2018SGeo...39..435F.html)).
- **MSSA is established.** Zotov et al. applied MSSA to GRACE Stokes coefficients as early as EGU 2010 ([abstract](https://meetingorganizer.copernicus.org/EGU2010/EGU2010-6288.pdf)). Wang et al. (GJI 2020; J. Hydrol. 2021) used improved MSSA for post-processing and for bridging the gap ([GFZ bib 2021](https://www-app2.gfz.de/grace-references/b2hd-2021JHyd..59425972W.html)). Gauer, Chanard & Fleitout used MSSA across multiple Level-2 solutions to bridge the gap and filter stripes ([ESSOAr](https://essopenarchive.org/doi/full/10.1002/essoar.10512455.1)). Rotated MSSA over the Tibetan Plateau reduced mode mixing and found interannual periods of about 6.1, 3.4 and 2.5 years ([J. Appl. Geod.](https://www.degruyterbrill.com/document/doi/10.1515/jag-2022-0002/html)).
- **Complex (Hilbert) EOF has been applied to GRACE, but for the ocean.** Yu, Chao, García-García & Luo (2018, JGR Oceans) found a leading CEOF mode in the Argentine Gyre: a counter-clockwise dipole with a period of about 25 days, missed by the de-aliasing ocean model ([HUST summary](https://ggg.hust.edu.cn/info/1238/3572.htm)). An EGU 2020 abstract extended GRACE CEOF to the other gyres ([EGU2020-4022](https://meetingorganizer.copernicus.org/EGU2020/EGU2020-4022.html)). The physical interpretation of the 25-day signal is disputed ([EGU2007 abstract](https://www.meetings.copernicus.org/www.cosis.net/abstracts/EGU2007/02443/EGU2007-J-02443.pdf)). Phillips et al. (2012, GRL) used a Hilbert transform to handle lags when correlating GRACE TWS with ENSO ([CU Experts](https://experts.colorado.edu/display/pubid_17410)).
- **DMD has been applied to GRACE TWS.**
  - Libero, Tartakovsky & Ciriello (EGU 2024) tested "a suite of different DMD algorithms" on 20 years of monthly GRACE TWS anomalies, identifying trend and seasonal modes. The abstract does not say whether the data were gridded or SH, and reports no forecast skill ([EGU24-18454](https://meetingorganizer.copernicus.org/EGU24/EGU24-18454.html)).
  - The journal paper is Libero & Ciriello (2025), *Adv. Water Resour.* 203:105015, "Dominant spatiotemporal structures in total water storage anomalies" ([GFZ bibliography](https://www-app2.gfz.de/grace-references/b2hd-2025AdWR..20305015L.html)).
  - A Stuttgart Geodetic Institute project proposes DMD to estimate a linear "network of water flux between the different sub-regions" from sub-basin TWS ([GIS Stuttgart project page](https://www.gis.uni-stuttgart.de/forschung/projekte/data-driven-modelling)).
- **DMD variants relevant to GRACE.**
  - Optimized DMD by variable projection (Askham & Kutz) handles "data which may be collected at unevenly spaced sample times" and shows "less bias in the presence of noise than standard DMD" ([arXiv:1704.02343](https://arxiv.org/abs/1704.02343)).
  - BOP-DMD bags many optimized-DMD fits from one snapshot set to give spatial and temporal uncertainty quantification and probabilistic forecasts ([arXiv:2107.10878](https://arxiv.org/abs/2107.10878)).
- **Koopman methods in climate, none in GRACE.**
  - Giannakis & Slawinska used data-driven Koopman eigenfunctions to recover ENSO as a coherent oscillator ([EGU2018-11002](https://meetingorganizer.copernicus.org/EGU2018/EGU2018-11002.pdf)).
  - Kernel analog forecasting predicted Niño 3.4 out to 10 months, reported as 3–8 months longer than a linear inverse model ([NSF PAR](https://par.nsf.gov/search/award_ids:1854383)).
  - Koopman ensembles beat linear models on 2000-year CESM runs but were less robust with short training records ([EGU24-2554](https://meetingorganizer.copernicus.org/EGU24/EGU24-2554.html)). This matters for our 259-month record.
- **Slow feature analysis (SFA) in geophysics, none in GRACE.** SFA has been used to extract slowly varying driving forces from non-stationary series, such as Arosa ozone ([Wang & Chen, NPG 2015](https://npg.copernicus.org/articles/22/377/2015/)), and to study intrinsic climate oscillations in proxy records ([Entropy 2021](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8070030/)). The results cluster in one research group.
- **No POP or LIM application to GRACE TWS was found.** The closest linear-predictive work is a 4-parameter linear recursive model driven by precipitation and temperature, which reconstructs daily TWS (Xie & Yi, HESS 2026) ([HESS](https://hess.copernicus.org/articles/30/3221/2026/hess-30-3221-2026.pdf)).
- **Tensor-product (not tensor-factorization) space–time models exist:** spherical wavelets × Legendre wavelets or B-splines (Fengler et al. 2007; Nutz & Wolf 2008; Schmidt et al. 2007) ([Nutz & Wolf](https://kluedo.ub.rptu.de/frontdoor/index/index/docId/1924); [Schmidt EGU2007](https://www.meetings.copernicus.org/www.cosis.net/abstracts/EGU2007/04079/EGU2007-J-04079.pdf)). No CP or Tucker decomposition of GRACE was found.

### Inferences
Implementation note that applies to every method in this section:
- Fully normalised real SH functions are orthonormal on the sphere (with the 4π convention). The Euclidean inner product of coefficient vectors therefore equals the area-weighted spatial inner product (Parseval).
- So EOF, DMD, LIM and SFA can all run **directly on the 259 × 8,281 coefficient matrix**, with no gridding and no cos-latitude weights. The SVD of that matrix is trivial on a laptop.
- To use the calibrated errors, first **whiten**: divide each coefficient by its sigma, or by a monthly-median sigma per (n, m). This turns plain PCA into a noise-weighted (generalised) PCA. A single linear change of variables reuses all downstream code.

Catalogue:

**D1. Optimized/BOP-DMD on sigma-whitened SH coefficients across gaps.** Novelty: *partial*. DMD on GRACE exists; the irregular-time, uncertainty-quantified, coefficient-domain variant was not found.
- *Does:* fits x(t) ≈ Σ_k b_k φ_k e^{ω_k t} with continuous-time eigenvalues ω_k. Variable projection lets the sample times be the true mid-epochs of each month, so the 2017–18 gap and the battery-era missing months need no interpolation. Bagging gives confidence intervals on the frequencies, damping rates and modes.
- *Why:* it gives annual, semi-annual and alias-period modes together with their damping, with error bars, from one fit. The same fit then interpolates across gaps (a deterministic gap filler) and forecasts.
- *Data/compute:* 259 snapshots; truncate to the top 20–50 whitened PCs first. Variable projection is a small nonlinear least-squares problem. Bagging multiplies the cost by about 100, which is still seconds to minutes.
- *Validate:*
  1. Use the Lecomte et al. (2024) remove-restore protocol (section 5). The target is to beat the AR benchmark of about 0.9–1.7 cm EWH.
  2. Check that the recovered annual eigenvalue period is 365.25 d within its CI.
  3. Synthetic test: known damped oscillators plus noise drawn from the COST-G sigmas.
- *Rust fit:* new module in `gravimetry` (or a new crate `gravimetry-ts`). Needs a dense SVD and an eigen-solver, such as `nalgebra` or `faer`.

**D2. POP / Linear Inverse Model (LIM) in EOF space.** Novelty: *none found* for GRACE.
- *Does:* fits x_{t+1} = A x_t + ξ_t on the leading k PCs of deseasonalised anomalies, with A = C(1) C(0)⁻¹. POPs are the eigenvectors of A. The noise covariance Q = C(0) − A C(0) Aᵀ.
- *Why:* Lecomte et al. found a univariate AR model was the best gap filler (section 5). A LIM is its multivariate, physically interpretable generalisation. It gives a forecast, a gap filler (via a Kalman smoother on the same state-space) and oscillatory modes such as ENSO-linked TWS dipoles.
- *Data/compute:* about 250 transitions, so keep k ≤ ~20. Missing months need EM or Kalman-based estimation of A. The compute is trivial.
- *Validate:*
  1. Forecast skill at 1–6 months versus persistence and damped persistence.
  2. Remove-restore gap test.
  3. Tau-test: check that A estimated at lag 2 is consistent with A at lag 1 squared, the standard LIM consistency check (background knowledge, not sourced here).

**D3. Slow feature analysis (SFA) for slow drivers.** Novelty: *none found* for GRACE.
- *Does:* finds projections w that minimise the variance of the month-to-month difference relative to total variance. This is a generalised eigenproblem: cov(Δx) w = λ cov(x) w.
- *Why:* it isolates multi-year drivers, such as aquifer depletion and decadal ocean-mass modes, without a parametric trend model. That makes it a model-free complement to BEAST trends.
- *Data/compute:* the same reduced PC matrix as D2. One small eigenproblem.
- *Validate:* synthetic test with an injected slow driver plus seasonal signal and sigma-noise. Compare with linear trends and with the leading EOF.
- *Risk:* the published uses come mostly from one group; scientific value is uncertain.

**D4. Koopman / kernel analog forecasting.** Novelty: *none found* for GRACE.
- Promising in climate, but the EGU24-2554 finding (less robust with short records) suggests 259 months is marginal.
- Treat it as a stretch goal behind D1/D2. The baseline is the LIM, the same baseline the ENSO Koopman work beat.

**D5. Complex/Hilbert EOF on land TWS.** Novelty: *partial*. CEOF on GRACE ocean signals exists; a land application to propagating flood waves or drought fronts was not found.
- Cheap to add once EOF exists: analytic signal per coefficient series via FFT Hilbert transform (`satgrav::noise::fft` already exists), then a complex SVD.
- *Validate:* Amazon seasonal flood wave propagation (west to east phase progression) against published phase maps. The value is visualisation (phase-propagation maps), not a new estimator.

**Skip as not novel:** ICA, rotated EOF, MSSA (implement only as baselines). **Skip as poor fit:** CP/Tucker tensor factorisation. The (C/S, n, m, t) layout is not a natural low-rank tensor.

### Gaps
- The Libero & Ciriello (2025) full text was not read. It is unknown whether they used optimized DMD, irregular sampling, or SH-domain inputs. Check before claiming novelty for D1.
- Search terms used: "dynamic mode decomposition GRACE terrestrial water storage"; "Koopman operator GRACE gravity time series"; "Koopman operator hydrology climate ENSO prediction data-driven"; "complex EOF Hilbert GRACE propagating mass signal"; "\"complex empirical orthogonal function\" GRACE"; "complex EOF ocean bottom pressure GRACE Argentine Gyre"; "principal oscillation pattern GRACE gravity"; "linear inverse model LIM GRACE terrestrial water storage prediction"; "slow feature analysis geophysical time series climate"; "rotated EOF GRACE terrestrial water storage modes"; "independent component analysis GRACE Forootan Kusche"; "multichannel singular spectrum analysis GRACE gap filling"; "tensor decomposition GRACE gravity spatiotemporal"; "functional data analysis GRACE terrestrial water storage curves". Functional data analysis (functional PCA) on GRACE: none found.

---

## 2. Spectral-domain methods unique to spherical harmonics

### Takeaway
Slepian localisation, spherical wavelets, Wiener filtering, DDK and VADER anisotropic decorrelation, sparse (L1) DDK, daily Kalman smoothing and per-coefficient Fourier analysis are all established. Three openings remain:
- needlets, for which no GRACE use was found;
- per-coefficient state-space (Kalman/GP) smoothing that uses each month's calibrated sigma as heteroscedastic observation noise;
- region-adaptive Wiener filters built from Slepian-multitaper local spectra.

With COST-G's diagonal sigmas, any covariance-optimal filter is diagonal unless correlated error is modelled separately (section 7).

### Cited Findings
- **Slepian functions (established):** Harig & Simons for Greenland (2012), Antarctica (2015) and GRL 2016 (fewer parameters by keeping only the most concentrated functions) ([GRL 2016 PDF](https://geoweb.princeton.edu/people/simons/PDF/reprints/GRL-2016.pdf); [AGU 2014](https://agu.confex.com/agu/fm14/webprogram/Paper23038.html)). Also Longuevergne et al. 2010 for the High Plains aquifer, and Iceland at −9.3 ± 1.0 Gt/yr for 2002–2016 (von Hippel & Harig 2019) ([Frontiers](https://www.frontiersin.org/journals/earth-science/articles/10.3389/feart.2019.00171/full)).
- **Spherical multitaper (established in general geophysics):** Wieczorek & Simons 2005 (localised admittance/coherence), 2007 (minimum-variance multitaper with Slepian tapers) ([GJI 2005](https://geoweb.princeton.edu/people/simons/BIB/Wieczorek+2005-GJI.html); [JFAA 2007](https://link.springer.com/doi/10.1007/s00041-006-6904-1)). Simons, Dahlen & Wieczorek 2006 (SIAM Review) laid the Slepian foundations ([PDF](https://geoweb.princeton.edu/people/simons/PDF/reprints/SIAM-2006.pdf)). Dahlen & Simons 2008 compared periodogram, maximum-likelihood and multitaper estimators ([arXiv:0705.3083](https://arxiv.org/pdf/0705.3083)). No time-variable-gravity application of local multitaper spectra surfaced.
- **Spherical wavelets (established; needlets none found):** Schmidt et al. 2006 built regional wavelet multiresolution models of the Amazon from GRACE range-rates ([EarthRef](https://earthref.org/ERR/94011/)). Fengler et al. 2007 (J. Geod.) combined spherical wavelets in space with Legendre wavelets in time ([GFZ bib](https://www-app2.gfz.de/pb1/op/grace/references/b2hd-2007JGeod..81....5F.html)). A search for "needlet GRACE gravity localization" returned no needlet papers ([closest hit](https://jgit.kntu.ac.ir/article-1-167-en.html)).
- **Wiener filtering (established):**
  - Sasgen, Martinec & Fleming (2006) designed a Wiener filter from signal and noise degree spectra inferred from the GRACE degree-power spectrum. For April 2004 it closely resembles a Gaussian of about 440 km half-width ([GFZ bib](https://www-app2.gfz-potsdam.de/grace-references/b2hd-2006StGG...50.499S.html)).
  - An anisotropic Wiener-type filter following Klees et al. (2008) and Siemes et al. (2013) has been applied to GRACE trends ([TU Delft repository](https://repository.tudelft.nl/file/File_940d6e90-0d2a-409a-ad43-663571c5188c)).
  - A regularisation-type filter using the SH variance–covariance matrices appears in an EGU 2009 abstract ([EGU2009-7148](https://meetingorganizer.copernicus.org/EGU2009/EGU2009-7148.pdf)).
- **DDK and VADER (established):**
  - DDK (Kusche et al. 2009) uses error covariance plus a signal model, giving a spatially anisotropic filter ([Ju et al. 2023](https://www.irsm.cas.cz/materialy/acta_content/2023_doi/Ju_AGG_2023_0005.pdf); [ICGEM note](https://icgem.gfz-potsdam.de/docs/Note_on_the_DDK_filters.pdf)).
  - The DDK1–8 matrices come from a single month's (August 2003) error covariance, so they are static ([GSTM session](https://www2.csr.utexas.edu/grace/GSTM/a1.html)).
  - VADER, the time-variable decorrelation filter, reduced cumulative geoid errors by about 15% versus static decorrelation. For comparable strength it gave better destriping and smaller formal errors than DDK ([Horvath et al. 2018](https://gfzpublic.gfz-potsdam.de/pubman/item/item_3766971_3/component/file_3806902/3766971.pdf)).
- **Sparse filtering (partially explored):** Qian et al. (2022) replaced DDK's Tikhonov term with a weighted L1 norm (a "sparse DDK") ([TU Delft PDF](https://research.tudelft.nl/files/125307388/remotesensing_14_02810_v2.pdf)). Fischer & Michel (2012) developed sparse matching-pursuit regularisation for inverse gravimetry, with South American mass variations as a case study ([PDF](https://www.uni-siegen.de/fb6/geomathe/publications/fischer_michel_invprobl_2012.pdf)). Compressive sensing has been applied to GRACE-based groundwater monitoring in Turkey ([ITU](https://research.itu.edu.tr/en/publications/efficient-monitoring-of-groundwater-level-changes-using-compressi/)). LASSO has been used for sparse spherical-radial-basis regional gravity models ([Crossref record](https://api.crossref.org/works/10.3390%2FRS15194870)).
- **Kalman filtering of SH coefficients (established for daily fields):**
  - ITSG's daily solutions (degrees 2–40) use the Kurtenbach et al. (2012) Kalman smoother. Its process-noise correlations come from WGHM, ECMWF and OMCT model output for 1976–2000 ([ITSG-Grace2014 README](https://ftp.tugraz.at/pub/ITSG/GRACE/ITSG-Grace2014/README.txt)).
  - ITSG-Grace2018 stabilises the daily solutions with an order-3 AR model ([ITSG-Grace2018 daily README](https://ftp.tugraz.at/pub/ITSG/GRACE/ITSG-Grace2018/daily_kalman/README.txt)).
  - State-space models with stochastic trend and seasonal terms have been applied to GRACE regional series (Didova et al. 2016, J. Geod.) ([GFZ bib](https://www-app2.gfz.de/grace-references/b2hd-2016JGeod..90.1207D.html); [GSTM 2015](https://www2.csr.utexas.edu/grace/GSTM/2015/a1.html)). A Kalman framework for regional mass-change series also exists (Wöhnke et al. 2024, GEM) ([HCU](https://forschungsportal.hcu-hamburg.de/en/publications/kalman-filter-framework-for-a-regional-mass-change-model-from-gra/)).
- **Per-coefficient spectral analysis (established, limited):** Fourier analysis of selected GFZ coefficient series selects frequencies for a Gauss–Markov model with polynomial and periodic terms ([KIT](https://publikationen.bibliothek.kit.edu/1000171404/153024675); [Acta Geod. Geophys. 2024](https://link.springer.com/article/10.1007/s40328-024-00446-x)). A GRACE book chapter links a variability drop in zonal degrees 17–20 to pointing noise near 3.2–3.7 mHz, using the orbit-period degree-to-frequency mapping ([Austria-Forum web book](https://austria-forum.org/web-books/en/gracegravity00en2018isds/000168)).
- **SH-domain SSA for gaps (established):** Yi & Sneeuw (JGR Solid Earth) filled GRACE gaps with SSA in the SH domain ([Uni Stuttgart news](https://www.gis.uni-stuttgart.de/en/institute/news/news/Filling-the-Data-Gaps-Within-GRACE-Missions-Using-Singular-Spectrum-Analysis/)).
- **Errors are strongly correlated:** "errors in the SH coefficients are highly correlated" because of anisotropic K-band sensitivity ([J. Geod. 2017](https://link.springer.com/article/10.1007/s00190-017-1063-5)).

### Inferences
**S1. Per-coefficient state-space smoother with calibrated heteroscedastic noise.** Novelty: *partial*. Kalman on SH exists for daily fields with model-derived process noise, and state-space models exist for regional series. A monthly, per-coefficient smoother with process noise learned from the 259-month record itself and sigma_nm(t) as observation noise was not found.
- *Does:* for each of the 8,281 coefficients, runs a structural time-series model:
  - state = level, slope, and annual and semi-annual cos/sin terms;
  - random-walk process noise q, pooled per degree (or degree band) and fitted by maximum likelihood or EM;
  - observation noise = (s_n · sigma_nm(t))², with s_n a per-degree sigma-calibration factor that is also estimated.
  The RTS smoother then gives a smoothed field and a posterior variance for every month, including the missing ones.
- *Why:* it is temporally optimal under its model and handles gaps natively. It yields gap-fill uncertainties for free, which the Lecomte benchmark methods do not provide. It runs in O(8,281 × 259) per pass and is streaming-capable: a new GRACE-FO month is one filter step.
- *Connection:* a temporal Matérn-ν GP (ν = 1/2, 3/2) has an exact state-space form (Hartikainen & Särkkä; background knowledge, not sourced here). An *isotropic* spatial GP on the sphere is diagonal in SH (degree-variance spectrum). So S1, with degree-pooled hyperparameters, *is* an exact separable spatio-temporal Matérn GP on the sphere, at linear cost. This is a crisp, publishable framing (see section 5).
- *Limitation:* temporal-only. Stripes that persist between months are not removed. Combine with the existing `gravimetry::destripe` or with S3.
- *Validate:*
  1. Lecomte remove-restore (single months and an 11-month gap), RMS in cm EWH on a 1° grid.
  2. Coverage test: the fraction of withheld months inside the 1σ/2σ posterior bands should be about 68%/95%.
  3. Check the estimated s_n against 1 (are COST-G sigmas well calibrated?).

**S2. Needlet (and Slepian) multiscale analysis with sigma-aware thresholding.** Novelty: *none found* for needlets on GRACE.
- *Does:* multiplies coefficients by smooth band-pass windows b(n/B^j), each covering about a factor-B band of degrees. The windows form a tight frame, so the scales sum back to the field. Each scale is then synthesised on a grid. Noise variance per scale and location is closed-form from sigma: Σ_n b²(n/B^j) Σ_m σ²_nm Ȳ²_nm(θ, λ). That allows per-scale z-score thresholding (wavelet shrinkage).
- *Why:* scale-dependent destriping (stripes live at high degree and order), multiscale anomaly maps, and good spatial localisation without a regional basis. The needlet construction itself comes from CMB analysis (background knowledge, not verified here).
- *Data/compute:* J ≈ 6–8 scales for n ≤ 90. One synthesis per scale per month. Uses existing synthesis plus `degree_variances`.
- *Validate:* synthetic basin signals with sigma-noise and correlated stripe noise. Compare leakage and RMS against Gaussian 300 km, DDK3–5 and Sasgen-Wiener. Real-data check: Greenland and Amazon mass series against published values.

**S3. Region-adaptive Wiener filter from Slepian-multitaper local spectra.** Novelty: *partial*. Global Wiener and anisotropic Wiener filters exist; a filter whose signal spectrum is estimated locally (Greenland ≠ Sahara) from multitaper spectra was not found.
- *Does:* estimates the local signal degree spectrum S_n(region) from the 259-month record with Slepian tapers, and the local noise spectrum from sigma. It builds w_n = S_n/(S_n + N_n) per region and blends the regions with a partition of unity.
- *Why:* a single global Gaussian is too strong over ice sheets and too weak over deserts.
- *Needs:* a Slepian basis. For polar caps the localisation kernel is block-diagonal by order m, so the eigenproblems are tiny. Arbitrary regions need SH analysis/quadrature of products of Ylm, a missing capability.
- *Validate:* closed-loop with a synthetic truth (hydrology model or ESM), comparing leakage and RMS against a global Wiener filter.

**S4. Diagonal VADER (time-variable Wiener from COST-G sigmas).** Novelty: *low* (it is diagonal VADER), but useful and very cheap.
- *Does:* w_nm(t) = S_nm / (S_nm + σ²_nm(t)), with S_nm the empirical signal variance per coefficient from the record, after removing the sigma-variance.
- *Why:* noisier months (early mission, the 2016–17 battery period, early GRACE-FO) are automatically smoothed more.
- *Caution:* diagonal weights cannot remove correlated stripe structure; run after destriping.

**Not novel, implement as infrastructure only:** Slepian basis for regional mass (high value for users), Gaussian/DDK, SH-domain SSA.

### Gaps
- No source found for needlets in any gravimetry setting. The needlet-in-CMB background is from prior knowledge and was not verified this session.
- The search did not establish whether any group runs per-coefficient monthly Kalman smoothing with calibrated sigma as R(t). Check the Graz and Bonn literature, such as ITSG monthly processing and Kurtenbach's thesis, before claiming S1.
- Search terms used: "Slepian functions GRACE mass change Harig Simons"; "needlets spherical wavelets GRACE time-variable gravity"; "needlet GRACE gravity localization"; "spherical multitaper spectral analysis gravity Wieczorek Simons"; "compressed sensing sparse recovery GRACE gravity field"; "DDK filter full covariance Wiener regularization GRACE Kusche"; "Wiener filter GRACE Sasgen optimal filtering signal noise degree variance"; "VADER filter time-variable decorrelation GRACE Horvath anisotropic"; "Kalman filter spherical harmonic coefficients GRACE daily ITSG Kurtenbach"; "state space model Kalman filter time-variable trend GRACE mass time series Didova"; "degree order coefficient time series spectral analysis GRACE Stokes coefficients periodogram".

---

## 3. Change-point and anomaly detection in mass fields

### Takeaway
GRACE break detection so far uses simple tests: Pettitt, EMD-based breakpoints, BEAST as preprocessing, and extreme-value return levels. No GRACE use was found for Bayesian online change-point detection (BOCPD), PELT, matched-filter or GLRT detectors in the coefficient domain, or space–time scan statistics.

The strongest publishable gap is an **SH-domain matched-filter/GLRT scan**. The literature documents that destriping distorts coseismic signals. A matched filter weighted by the calibrated sigmas needs no destriping, and for isotropic templates the scan over all locations reduces to one SH synthesis per epoch.

### Cited Findings
- **BOCPD:** Adams & MacKay (2007) is an exact online method for the run length since the last change ([arXiv:0710.3742](https://arxiv.org/pdf/0710.3742)). An extension detects mean and slope (baseline) shifts in drifting series ([arXiv:2201.02325](https://arxiv.org/pdf/2201.02325)). No GRACE application surfaced.
- **BEAST:** a Bayesian model-averaging decomposition into trend, seasonality and abrupt change, with probabilities of change (Zhao et al. 2019) ([Rbeast CRAN](https://ftp.fau.de/cran/web/packages/Rbeast/index.html)). It has been used on GRACE only as a piecewise-linear detrending step before gap filling (Qian et al., *Remote Sensing* 2024) ([MDPI](https://www.mdpi.com/2072-4292/16/19/3693)).
- **Other break tests on GRACE:**
  - Pettitt's test on 285 Indian grid points for 2002–2020 ([Sustainability 2023](https://www.mdpi.com/2071-1050/15/2/1572)). Pettitt detects only a single change and is weak in small samples ([bootstrap Pettitt, arXiv:2411.05233](https://arxiv.org/pdf/2411.05233)).
  - An EMD-plus-breakpoint method on the 35 largest basins matched SPI anomalies and documented extremes ([Remote Sensing 2017](https://www.mdpi.com/2072-4292/9/12/1287)).
  - A claimed global TWS regime shift in 2015–2016 linked to El Niño ([EGU25-17534](https://meetingorganizer.copernicus.org/EGU25/EGU25-17534.html)).
- **Extreme value theory (established):** Middendorf, Dobslaw, Jensen & Eicker (JGR Solid Earth 2025) fitted GEV distributions to GRACE TWS for 1-in-10-year dry return levels and compared them with CMIP6 ([GFZ](https://gfzpublic.gfz.de/pubman/item/item_5037428)).
- **Coseismic detection is filter-sensitive:**
  - Commonly used destriping "tends to distort the true coseismic signals in both amplitude and spatial pattern" (Tohoku, *Earth Planets Space* 2015) ([EPS](https://earth-planets-space.springeropen.com/articles/10.1186/s40623-015-0188-0)).
  - Wang, Shum, Simons et al. (GRL 2012) avoided destriping and used Slepian localisation with a 350 km Gaussian ([PDF](https://geoweb.princeton.edu/people/simons/PDF/reprints/GRL-2012a.pdf)).
  - Han & Simons (2008) made the Sumatra–Andaman signal detectable with spatiospectral windows ([PDF](https://igpppublic.ucsd.edu/~shearer/Files/Sumatra_Papers/han_jgr08.pdf)).
  - A GFZ V-shape SH-domain filter improved detection of Mw 8.0–8.2 events in the GRACE-FO period ([GFZ](https://gfzpublic.gfz-potsdam.de/pubman/item/item_5020737)).
  - GRACE detected Sumatra-Andaman, Maule and Tohoku-Oki. Detection depends on event type, area and record length ([Stuttgart thesis](https://elib.uni-stuttgart.de/items/64132fac-f96b-4e98-95f6-1cc8882bad83/full); [Matsuo & Heki 2011](https://eprints.lib.hokudai.ac.jp/repo/huscap/all/48676/2011GL049018.pdf)).
- **Matched filters in GRACE: only for accelerometer calibration.** Pagiatakis & Tzamali used penumbra transitions as known pulses ([EGU23-10423](https://meetingorganizer.copernicus.org/EGU23/EGU23-10423.html)). No gravity-field matched filter surfaced.
- **Scan statistics:** Kulldorff space–time and expectation-based scan statistics are implemented in R `scanstatistics` ([JOSS](https://joss.theoj.org/papers/10.21105/joss.00515)). They are compared with multivariate SPC and tensor methods for satellite-image anomalies ([arXiv:2309.09878](https://arxiv.org/pdf/2309.09878)). No gravity application surfaced.
- **Flood detection from daily GRACE:** Kalman-filtered daily solutions tracked the 2004 and 2007 Ganges–Brahmaputra floods, the first hydrological evaluation of daily GRACE solutions ([Gouweleeuw et al., HESS 2018](https://hess.copernicus.org/articles/22/2867/2018/)).

### Inferences
**C1. SH-domain GLRT / matched-filter scan for abrupt localised mass steps.** Novelty: *none found*. Highest-priority publishable candidate.
- *Model:* Δd_nm = a · g_nm(θ₀, λ₀) + ε_nm, with ε_nm ~ N(0, σ²_nm), where:
  - Δd = (mean of k months after t₀) − (mean of k months before t₀), after removing a seasonal/trend fit;
  - g = template (disc or cap load of radius r, converted to Stokes coefficients with the existing Love numbers);
  - σ²_nm = propagated variance of the difference.
- *Statistic:* T = (Σ g Δd/σ²)² / Σ g²/σ², which is χ²₁ under the null.
- *Key trick:* for an axisymmetric template K(ψ) = Σ_n k_n P_n(cos ψ), centring at (θ₀, λ₀) gives g_nm = k_n Ȳ_nm(θ₀, λ₀)/(2n+1) (addition theorem). So the numerator over *all* centres is a **single SH synthesis** of the field k_n Δd_nm/((2n+1) σ²_nm). The denominator is another synthesis of order-dependent squared-Legendre weights. A full global scan over t₀ × radii costs about 259 × R × 2 syntheses, which is seconds.
- *Why:* it is optimal (Neyman–Pearson) under Gaussian diagonal noise, needs **no destriping** (avoiding the documented distortion), uses the calibrated sigmas directly, and gives p-values.
- *Correlated stripes break the diagonal assumption.* Calibrate thresholds empirically, Kulldorff-style:
  - Monte Carlo nulls from time-shuffled or sign-flipped residual months;
  - report family-wise p-values over space × time × radius.
  This combines matched filtering with scan statistics (C3).
- *Extensions:*
  - double-couple coseismic templates (needs dislocation Green's functions and SH rotation; FastTransforms or ducc0 have rotations, see section 7);
  - ramp templates for drought onset;
  - a sequential (CUSUM) version for operational monitoring.
- *Validate:*
  1. Blind recovery of Sumatra (Dec 2004), Maule (Feb 2010) and Tohoku (Mar 2011) locations and epochs.
  2. Detection-threshold curve (minimum detectable Mw or Gt versus radius) from injection tests into real residuals.
  3. Hydrology events: reservoir impoundments and lake drainages with known volumes.
  4. False-alarm rate on GRACE-FO-only years.

**C2. BOCPD on basin, Slepian or PC series with known noise.** Novelty: *none found*.
- *Does:* run-length posterior with a Gaussian likelihood. Observation variance is known from propagated sigma: for a basin average with SH weights w_nm, the variance is Σ w²_nm σ²_nm. The level (and slope, via the baseline-shift extension) is unknown with a conjugate prior. Hazard 1/λ, with λ around 60 months.
- *Why:* online, probabilistic alerts as each GRACE-FO month arrives. It fits a streaming pipeline (S1 provides the deseasonalised residual).
- *Compute:* O(run-length) per month, about 200 lines of Rust.
- *Validate:* synthetic steps and ramps with sigma noise; average run length to false alarm. Real events: the coseismic list above and the 2015–16 global regime shift claim (EGU25-17534).

**C3. Gaussian space–time scan statistic.** Novelty: *none found*.
- Kulldorff scan with a Gaussian (not Poisson) likelihood over cylinders (cap × time window). It is essentially C1 with box templates in time, so implement it as a C1 option.

**C4. PELT / multiple change points on coefficient or PC series.** Novelty: *not specifically searched*. No GRACE PELT hit appeared in the breakpoint searches.
- A cheap offline companion to BOCPD with an exact penalised segmentation. Low priority.

### Gaps
- Not verified: whether Okubo-type coseismic gravity templates in SH form are already used in a likelihood-ratio framework under another name, such as "template fitting" in the Chao/Han literature. A targeted Scholar search on "GRACE coseismic likelihood ratio" and "template fitting" is needed before publication.
- No dedicated search for "PELT GRACE" or "CUSUM GRACE" (search budget).
- Search terms used: "Bayesian online change point detection GRACE water storage"; "BEAST Bayesian change point GRACE terrestrial water storage breakpoint"; "breakpoint detection GRACE terrestrial water storage Pettitt abrupt change"; "coseismic gravity change detection GRACE matched filter earthquake"; "\"matched filter\" GRACE gravity"; "spatio-temporal scan statistic satellite gravity anomaly detection"; "extreme value theory GRACE drought flood return period terrestrial water storage"; "GRACE flood detection daily gravity solutions Gouweleeuw".

---

## 4. Causal and nonlinear dependence: PCMCI, CCM, transfer entropy, recurrence analysis, mutual information

### Takeaway
GRACE dependence studies use correlation, Hilbert-lagged correlation and Granger causality. No GRACE application was found for PCMCI, convergent cross mapping, transfer entropy, Liang–Kleeman information flow or recurrence quantification. All are novel *applications*, but the gain is scientific (needs external climate data), not algorithmic. With 259 months, linear or Gaussian variants (PCMCI with ParCorr, Liang–Kleeman) are statistically safer than kernel or nearest-neighbour estimators (CCM, KSG transfer entropy, RQA).

### Cited Findings
- **PCMCI:** constraint-based causal discovery for multivariate time series (PC-stable parent selection, then momentary conditional independence tests). PCMCI+ handles contemporaneous links, and LPCMCI handles hidden confounders ([review arXiv:2404.05746](https://arxiv.org/pdf/2404.05746)). No PCMCI–GRACE paper surfaced. A causal-inference TWS reconstruction for the Ganges has an unspecified method ([EGU24-971](https://meetingorganizer.copernicus.org/EGU24/EGU24-971.html)).
- **Granger causality is used with GRACE:** vegetation–TWS in China ([Xiao et al. 2023](https://www-app2.gfz-potsdam.de/grace-references/b2hd-2023JHyd..62630313X.html)) and in Australia ([Chen, Guan & Batelaan 2022](https://www-app2.gfz-potsdam.de/grace-references/b2hd-2022JHyd..61328336C.html)).
- **CCM in hydrology, not with GRACE:** soil moisture → precipitation, with the strongest lag around one month ([Sci. Rep. 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6093863/)). Groundwater–streamflow in Victoria, with the authors warning that link variability "can be due to some limitations of the method" ([Monash](https://research.monash.edu/en/publications/identifying-causal-interactions-between-groundwater-and-streamflo/)).
- **Transfer entropy in hydrology, not with GRACE:** directional information flow among precipitation, runoff, ET, soil moisture and SWE in the Ebinur Basin, from VIC model output ([Frontiers in Water 2025](https://www.frontiersin.org/journals/water/articles/10.3389/frwa.2025.1622980/pdf)).
- **ENSO–TWS coupling is correlation-based:** Phillips et al. 2012 ([CU](https://experts.colorado.edu/display/pubid_17410)); correlations of about −0.4 to −0.6 in a 2026 HESS paper ([HESS](https://hess.copernicus.org/articles/30/4771/2026/)).
- **RQA:** a nonlinear method quantifying recurrences in phase space, with no stationarity assumption ([Wikipedia](https://en.wikipedia.org/wiki/Recurrence_quantification_analysis); [crqa, arXiv:2006.01954](https://arxiv.org/pdf/2006.01954)). It has been used on gravitational-wave signals but not GRACE.
- **Liang–Kleeman information flow:** a search pairing it with GRACE returned only general GRACE papers and Granger-causality studies. The closest general reference returned was a trend/seasonal/extremes decomposition of GRACE TWS ([Humphrey et al. 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4944666/)), not an information-flow study.

### Inferences
**K1. PCMCI on basin TWS × climate drivers.** Novelty: *none found* (application).
- Variables: deseasonalised basin or PC series (from section 1 or 2), ENSO/IOD indices, GPCP or ERA5 precipitation, temperature. Lags 0–6 months; about 10 variables is feasible with ParCorr at n ≈ 250.
- *Value:* separates direct from mediated ENSO → precipitation → TWS links, and detects teleconnection paths. Rust has no PCMCI; ParCorr-based PCMCI is a few hundred lines (regressions plus partial-correlation tests).
- *Validate:* synthetic VAR systems with a known graph and GRACE-like sigma noise; physical sanity checks (P → TWS at lag 0–2).

**K2. Liang–Kleeman information flow.** Novelty: *none found*.
- Closed-form from sample covariances of series and their differences, so it is cheap and well-suited to short records (background knowledge, not sourced here). A good first screen before K1.

**K3. CCM, transfer entropy, RQA, mutual information.** Novelty: *none found*, but low priority.
- Seasonality creates spurious CCM links through shared forcing, and nearest-neighbour estimators are noisy at n ≈ 259.
- If attempted: deseasonalise, use surrogate-data significance (phase-randomised or twin surrogates), and report only robust links.

### Gaps
- No search specifically for "mutual information GRACE" (search budget). Treat it as unverified.
- External drivers (precipitation, ENSO indices) are needed and are not in our current data set.
- Search terms used: "PCMCI causal discovery GRACE terrestrial water storage"; "convergent cross mapping GRACE groundwater precipitation causality"; "transfer entropy GRACE water storage ENSO"; "transfer entropy terrestrial water storage soil moisture information theory"; "Liang-Kleeman information flow GRACE causality water storage"; "recurrence quantification analysis GRACE gravity".

---

## 5. Gap filling and super-resolution: GPs on the sphere, neural fields, PINNs, ensemble Kalman, the 2017–2018 gap

### Takeaway
Gap filling is a crowded field: SSA/MSSA, AR, PCA with priors, BEAST + ML, Bayesian CNNs, GAIN, spatio-temporal GNNs, GP–random-forest, and EnKF with hydrology models. A 2024 J. Geodesy benchmark ranked a plain AR model best (about 0.9–1.7 cm EWH) and Swarm-only worst (7.4–7.6 cm). That benchmark is the yardstick to beat.

What is missing:
- an exact, linear-cost spatio-temporal GP on the sphere that uses calibrated heteroscedastic errors and returns calibrated uncertainties (S1);
- multi-source Kalman fusion of COST-G, Swarm and SLR in the SH domain, weighting each by its own sigma.

Neural fields for gravity exist only for *static* fields so far. PINNs with an explicit continuity residual were not found.

### Cited Findings
- **The 2017–18 benchmark** (Lecomte, Rosat & Mandea, J. Geod. 2024) ([PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC11585520/)):
  - *Methods compared:* Swarm replacement; CTAS fits; PCA with Swarm or GLDAS priors; iterative ICA; SSA; improved MSSA; AR with GLDAS prior.
  - *Protocol:* remove-restore; latitude-weighted RMS on a 1° grid; bootstrapped.
  - *Results:* AR best at about 0.9–1.7 cm EWH. SSA second for single months; MSSA second for the 11-month gap (1.9 cm at degree 40, versus 2.1 cm for Wang et al. 2021). CTAS about 2.0 cm. Swarm replacement 7.4/7.6 cm. SLR solutions had about 2.6–3.6 cm error, larger than the methods' accuracy.
  - They cite Löcher & Kusche (2020): GRACE EOFs plus SLR for low degrees.
- **Other gap-filling approaches:**
  - BEAST detrending + drivers (Qian 2024, correlations above 0.9 in most basins) ([MDPI](https://www.mdpi.com/2072-4292/16/19/3693)).
  - Bayesian CNNs ([arXiv:2101.09361](https://arxiv.org/pdf/2101.09361)).
  - GAIN imputation for the Pearl River ([MDPI](https://www.mdpi.com/2072-4292/16/24/4727/xml)).
  - Spatio-temporal GNNs for South America ([arXiv:2606.23833](https://arxiv.org/pdf/2606.23833)).
  - GP + random forest for groundwater gaps and 1 km downscaling, R = 0.88 against wells ([Remote Sensing 2026](https://mdpi.com/2072-4292/18/16/2702)).
  - A CNN mapping SH coefficients directly to TWS ([ESS 2023](https://www-app2.gfz.de/grace-references/b2hd-2023ESS...1003023Z.html)).
- **GPs on the sphere (theory):**
  - Spherical Matérn as the solution of a fractional SPDE with the Laplace–Beltrami operator ([gstlearn SPDE on the sphere](https://soft.mines-paristech.fr/gstlearn/1.11.1/references/SPDE_Sphere.md)); surface-FEM error analysis ([arXiv:2102.08822](https://ar5iv.labs.arxiv.org/html/2102.08822)).
  - Matérn on spheres has smoothness restrictions, which the 𝓕-family addresses ([arXiv:2101.05394](https://ar5iv.arxiv.org/html/2101.05394)).
  - Kriging is "also known as Gaussian process regression" ([Wikipedia](https://en.wikipedia.org/wiki/Kriging)).
  - Spatio-temporal GPR can be run with a Kalman smoother ([arXiv:2209.12565](https://arxiv.org/pdf/2209.12565)).
- **Hierarchical Bayesian combination (established):**
  - Zammit-Mangion et al. combine GRACE, altimetry and GPS in a hierarchical framework for Antarctica; Finn Lindgren, an SPDE developer, is co-author on the 2015 paper ([Bristol](https://research-information.bris.ac.uk/en/publications/resolving-the-antarctic-contribution-to-sea-level-rise-a-hierarch/); [EGU2016-13180](https://meetingorganizer.copernicus.org/EGU2016/EGU2016-13180.pdf?pdf), East Antarctica 57 ± 20 Gt/yr).
  - The ERC GlobalMass project (2016–2022) used a Bayesian hierarchical model for the sea-level budget ([CORDIS](https://cordis.europa.eu/project/id/694188/reporting)).
- **Neural fields:**
  - Hybrid SH + SIREN models of the *static* EGM2008 field (Restrepo Berrío, Kusche, Springer, Rußwurm) beat linear SH up to about 200,000 parameters (≈ degree 446). Inference took 0.028 s versus 205.9 s for 250,000 points. The work is 2D only, with no time variation and no GRACE ([EGU26-14971](https://meetingorganizer.copernicus.org/EGU26/EGU26-14971.html)).
  - Herglotz-NET is an INR with harmonic positional encoding for spherical data ([arXiv:2502.13777](https://arxiv.org/html/2502.13777v1)).
  - INRs have been used for 3D gravity inversion ([arXiv:2510.17876](https://arxiv.org/html/2510.17876v1)).
- **PINN-like downscaling:**
  - Self-supervised GRACE ocean-bottom-pressure downscaling to 0.25°, with GRACE as the mass-conservation constraint ([arXiv:2404.05818](https://arxiv.org/pdf/2404.05818)).
  - Custom-loss TWS downscaling from 300 km to 50 km with WaterGAP correlations ([GSTM2024-54](https://meetingorganizer.copernicus.org/GSTM2024/GSTM2024-54.html)).
  - No explicit continuity-residual PINN surfaced.
- **Spherical neural operators:** SFNO generalises FNOs to the sphere ([arXiv:2306.03838](https://ar5iv.labs.arxiv.org/html/2306.03838)), with a GPU SHT in NVIDIA torch-harmonics ([tutorial](https://nvidia.github.io/torch-harmonics/tutorials/train_spherical_neural_operator.html)). No GRACE use found.
- **EnKF data assimilation (established):**
  - WaterGAP with gridded GRACE and error correlations (Eicker 2014) ([HCU](https://forschungsportal.hcu-hamburg.de/en/publications/calibrationdata-assimilation-approach-for-integrating-grace-data-/)).
  - Rhine: groundwater correlation 0.31 → 0.53, RMSE 8.4 → 5.4 cm ([HESS 2015](https://hess.copernicus.org/articles/19/2079/2015/)).
  - Full error covariances in PCR-GLOBWB ([HESS 2017](https://hess.copernicus.org/articles/21/2053/2017/)).
- **Daily Kalman fields with predictions in gaps "tend toward the mean trend and annual signal"** ([ITSG-Grace2014 README](https://ftp.tugraz.at/pub/ITSG/GRACE/ITSG-Grace2014/README.txt)).

### Inferences
**G1. Exact spherical spatio-temporal Matérn GP = S1.** Novelty: *partial / framing new*.
- An isotropic GP on the sphere has covariance diagonal in (n, m), with variance depending only on n. A separable temporal Matérn has a state-space form. So the full GP posterior factorises into 8,281 scalar Kalman smoothers with degree-tied hyperparameters, at cost O(N_coeff · T) rather than O((N_grid · T)³).
- Heteroscedastic, calibrated observation noise σ_nm(t) enters exactly.
- *Deliverable:* gap-filled fields with posterior variances, which no method in the Lecomte benchmark provides.
- *Non-isotropy* (ice sheets versus deserts) is the main weakness. A fix is a mixture over S3 regions, or a "warped" variance field (non-stationary amplitude via a spatial multiplier). That breaks diagonality, but CG with SH-domain operators remains cheap.

**G2. Multi-source Kalman fusion in the SH domain.** Novelty: *partial*. Löcher & Kusche use EOFs plus SLR; Lecomte et al. use Swarm as a PCA prior.
- In the S1 state-space, add Swarm (degree ≤ about 40, but useful to about 12–15) and SLR (degree ≤ about 5–10) monthly fields as extra observations with *their own* sigmas and per-source calibration factors estimated by ML.
- Expected to beat CTAS within the gap and to give honest uncertainties.
- *Data:* Swarm EGF (`gravity-field` already parses it) and an SLR series (needs a reader).
- *Validate:* withhold an 11-month GRACE block in 2015–16 (Swarm overlap) and compare against AR, MSSA and Swarm-only using the Lecomte metric.

**G3. Neural fields / PINNs.** Novelty: time-variable SH neural fields and continuity-residual PINNs not found.
- Heavy for a Rust stack and not clearly better than G1 at degree 90. A defensible path is a small Rust MLP (e.g. `candle` or `burn`) on SH positional encodings with a sigma-weighted loss, as a research toy. Defer.

**Not recommended as new work:** EnKF with hydrology models (established, needs a land model); CNN/GNN gap fillers (crowded field).

### Gaps
- No source was obtained for the geodetic equivalence "least-squares collocation ≡ kriging ≡ GP regression". It is well known in geodesy but was not verified here (search budget).
- It is unknown whether Lecomte's AR used per-coefficient or per-grid AR. Read the full text before designing the comparison.
- Search terms used: "Gaussian process regression GRACE gap filling spatiotemporal"; "Matérn Gaussian process sphere SPDE stochastic partial differential equation GRACE"; "implicit neural representation spherical harmonics neural field gravity"; "physics-informed neural network GRACE downscaling mass conservation"; "spherical Fourier neural operator GRACE terrestrial water storage"; "GRACE GRACE-FO gap bridging Swarm SLR machine learning comparison 2017 2018"; "ensemble Kalman filter GRACE data assimilation hydrology model"; "Zammit-Mangion hierarchical Bayesian Antarctica GRACE altimetry SPDE INLA"; "GlobalMass Bayesian hierarchical model sea level budget GRACE".

---

## 6. Inverse problems: mass-flux and transport inference (Helmholtz on the sphere), sea-level fingerprints, sparse source localisation

### Takeaway
Three things already exist: inferring lateral water fluxes from GRACE with process models (diffusion–advection with an adjoint), DMD network fluxes between sub-basins (proposed), and the atmospheric moisture-flux water budget.

A **model-free, minimum-norm transport field** computed directly from the SH mass tendency was not found. It is the Helmholtz/Poisson solution ψ_nm = R² σ̇_nm / (n(n+1)), F = ∇ψ. It is nearly free to compute with our crates, naturally suppresses high-degree noise, and gives "where water moves" maps. It extends to linearised and unbalanced optimal transport, also not found for GRACE.

Sea-level fingerprint inversion is established, so treat it as infrastructure. Sparse source localisation is partially explored (sparse DDK, matching pursuit).

### Cited Findings
- **Process-model lateral fluxes:** Ramillien, Darrozes & Seoane (EGU 2026) couple a surface diffusion–advection model with WGHM assimilation of gridded GRACE. They estimate subsurface diffusivity and P−ET biases by adjoint-based optimisation to "diagnose effective surface diffusivity and lateral water fluxes" ([EGU26-4049](https://meetingorganizer.copernicus.org/EGU26/EGU26-4049.html)).
- **Data-driven sub-basin flux network:** linear storage–drainage model; DMD "seems to be a promising algorithm to determine this matrix" ([Uni Stuttgart GIS](https://www.gis.uni-stuttgart.de/forschung/projekte/data-driven-modelling)).
- **Difficulty:** Kusche et al. (GRL 2016) note GRACE captures lateral redistribution at grid scale. They caution that short-term lateral transports are hard to infer because month-to-month variability is noisy and resolution limited ([HAL PDF](https://hal-insu.archives-ouvertes.fr/insu-01355127/file/Kusche_et_al-2016-Geophysical_Research_Letters.pdf)).
- **Water-budget framing:** vertically integrated moisture convergence (P − E) has been compared with GRACE storage change (R − dS/dt) ([EGU2009-10219](https://meetingorganizer.copernicus.org/EGU2009/EGU2009-10219.pdf)). Helmholtz decomposition of *atmospheric* moisture-flux divergence exists ([Frontiers 2022](https://www.frontiersin.org/journals/earth-science/articles/10.3389/feart.2022.847317/epub)). No Helmholtz split applied to GRACE-derived mass tendency surfaced.
- **Ocean transport from GRACE (established, a related concept):** GRACE ocean bottom pressure gives lower-NADW (AMOC) transport anomalies agreeing with RAPID/MOCHA at 26°N ([Landerer, US CLIVAR](https://usclivar.org/sites/default/files/amoc/2015/LandererAMOC2014.pdf)). Synthetic tests show interannual AMOC anomalies recoverable with errors under 1 Sv in many locations ([Bentel et al., Ocean Sci. 2015](https://os.copernicus.org/articles/11/953/2015/)).
- **Laplace–Beltrami eigenfunctions:** spherical harmonics are eigenfunctions of the Laplace–Beltrami operator, with r²∇²Y = −ℓ(ℓ+1)Y ([Wikipedia](https://en.wikipedia.org/wiki/Spherical_harmonics)).
- **Optimal transport:** Wasserstein-2 equals the minimum kinetic energy of a density flow (Benamou–Brenier) ([arXiv:1611.09946](https://arxiv.org/pdf/1611.09946)). Dynamic *unbalanced* W1-type transport allows mass creation and destruction ([ESAIM: COCV 2019](https://www.esaim-cocv.org/10.1051/cocv/2018017)). No OT–GRACE work surfaced.
- **Sea-level fingerprints (established):**
  - Hsu & Velicogna (GRL 2017) derived fingerprints from GRACE and validated them with tropical ocean-bottom-pressure records ([UCAR](https://impacts.ucar.edu/en/publications/detection-of-sea-level-fingerprints-derived-from-grace-gravity-da/)). Of about 1.8 mm/yr, roughly 43% came from Greenland, 16% from Antarctica and 30% from glaciers ([NASA](https://sealevel.nasa.gov/news/95/nasauci-find-evidence-of-sea-level-fingerprints)).
  - Adhikari et al. (ESSD 2019) give monthly fingerprints from RL06 Stokes coefficients, with and without rotational feedback, in CM and CF frames ([ESSD](https://essd.copernicus.org/articles/11/629/2019)).
  - Joint GRACE + altimetry fingerprint inversions exist; about 100 fingerprints explained 1.56 of 1.94 mm/yr for 2002–2009 ([EGU2013-10331](https://meetingorganizer.copernicus.org/EGU2013/EGU2013-10331.pdf); [CU ref](https://sealevel.colorado.edu/bibcite/reference/113)).
- **Sparse sources:** sparse DDK ([Qian 2022](https://research.tudelft.nl/files/125307388/remotesensing_14_02810_v2.pdf)); matching pursuit for gravimetry ([Fischer & Michel 2012](https://www.uni-siegen.de/fb6/geomathe/publications/fischer_michel_invprobl_2012.pdf)). An EEG analogy notes that L2 "overestimate[s] the activated spatial areas for focal sources" ([HAL](https://hal-univ-tlse3.archives-ouvertes.fr/hal-01261357)). No L1 point-mass source localisation from GRACE coefficients surfaced.

### Inferences
**F1. Minimum-norm mass-transport field ("where water moves").** Novelty: *none found*. High visual and publication value, trivial compute.
- *Physics:* the continuity equation on the sphere is ∂σ/∂t + ∇·F = s, where σ is surface mass density (kg m⁻²) and s is the vertical source (P − E for land, net freshwater flux for the ocean).
- *Transport-only explanation (s = 0):* take the irrotational F = ∇ψ with ∇²ψ = −σ̇. In SH, ψ_nm = R² σ̇_nm / (n(n+1)) for n ≥ 1. σ̇ comes from finite differences or the S1 smoother's slope state, and σ_nm from the existing `gravimetry` Love-number anomaly code. Then F = ∇ψ is evaluated with `geopotential`'s existing first partials (θ and λ derivatives).
- *n = 0 must vanish (global mass conservation),* which holds for GSM monthly fields. Degree 1 comes from TN-13 and C20 from TN-14; the readers exist.
- *Why it is the right estimator:* among all fields F with the observed divergence, ∇ψ has minimum ∫|F|². Adding any divergence-free G adds ∫|G|², because the cross term ∫∇ψ·G = −∫ψ ∇·G = 0 on a closed sphere. So it is the least-motion explanation of the observed mass change. Units: F in kg m⁻¹ s⁻¹ (mass flux per unit width).
- *Noise robustness:* |F| scales like σ̇_n / n. High-degree stripes are damped by 1/n relative to the mass map, so flux maps are inherently smoother than the anomaly maps they come from.
- *Extensions:*
  - (a) *Water-budget residual flux:* use s = P − E from ERA5 (external data). Solving ∇²ψ = s − σ̇ gives a lateral (river + groundwater + ocean) transport field. This mapped version of the basin-integrated budget literature was not found.
  - (b) *Masked transport:* constrain F to land or to ocean with a masked Poisson problem. Solve with CG using SH synthesis/analysis per iteration, which needs SH analysis (missing).
  - (c) *Density-weighted (linearised W2):* ∇·(ρ ∇ψ) = −σ̇, with ρ a reference mass field such as the mean ice or water load. This makes transport prefer paths through existing mass.
  - (d) *Full unbalanced optimal transport* between monthly fields, using Sinkhorn iterations whose Gibbs kernel is a heat kernel on the sphere. That kernel is diagonal in SH, so each iteration reuses `gravimetry::gaussian_weights`-style smoothing. Convolutional Sinkhorn via heat kernels is from Solomon et al. (2015), background knowledge not verified here.
- *Byproduct metric:* the flux energy Σ R² |Δσ_nm|² / (n(n+1)) is the H⁻¹ (negative-Sobolev, linearised-Wasserstein) distance between two fields. It is a closed-form, displacement-tolerant metric for comparing solutions (COST-G vs CSR/GFZ/JPL) or validating gap fillers. Unlike RMS, it does not doubly penalise a slightly shifted blob. Not found used in gravimetry.
- *Validate:*
  1. *Synthetic:* transport a Gaussian blob along a known path on the sphere and check that recovered F matches the true divergent part.
  2. *Model truth:* take a hydrology or ocean model with known lateral fluxes (e.g. ESM or OMCT components, AOD1B), synthesise its σ, and compare F with the model's divergent flux component.
  3. *Real-data sanity:* Greenland outflow should point seaward; Amazon seasonal flux should reverse sign between rising and falling seasons; in the global land–ocean exchange, the n=0 ocean/land split should match barystatic sea-level estimates.
- *Caveat:* F is a kinematic diagnostic, not the true water path. Rivers and groundwater flow are not minimum-norm, and s = 0 attributes P − E to transport. Present it as "minimum transport consistent with GRACE".

**F2. Sea-level equation solver (infrastructure).** Novelty: *none*. Needed for fingerprints, barystatic budgets and GIA-consistent ocean masks.
- Pseudo-spectral SLE iteration: ocean function × SH analysis, Love numbers (existing), optional rotational feedback.
- *Validate:* against the Adhikari et al. (2019) published fingerprints.
- Needs SH analysis on a grid and an ocean-function grid.

**F3. Sparse (L1 / group-LASSO) localisation of mass sources from coefficients.** Novelty: *partial*.
- *Does:* dictionary of cap/disc templates at grid centres and multiple radii. Weighted L1 with sigma-whitened residuals, solved by FISTA. Dictionary products use one synthesis and one analysis per iteration.
- *Why:* compact sources (glaciers, reservoirs, earthquakes) are blurred by L2 methods.
- *Pairing:* C1's matched-filter map is the first iteration of matching pursuit, so F3 is the natural "detect then localise" second stage.
- *Validate:* injected point and disc loads into real residual months; localisation error versus magnitude.

### Gaps
- The full text of EGU26-4049 is unavailable; its flux formulation and resolution are unknown. If it also outputs flux vectors, F1's novelty is "model-free and minimum-norm" rather than "first flux map".
- No literature check on whether the H⁻¹ distance has been used to compare GRACE solutions, or whether OT is used in sea-level or ice-sheet comparisons (search budget exhausted before the "Wasserstein forecast verification geoscience" query ran).
- Search terms used: "horizontal water mass flux divergence from GRACE moisture transport Helmholtz"; "lateral water transport inferred from GRACE storage change divergence flux potential"; "GRACE ocean bottom pressure meridional transport AMOC Landerer"; "optimal transport Wasserstein distance GRACE mass change"; "sea level fingerprint inversion GRACE sources Hsu Velicogna"; "sea level equation inversion GRACE Adhikari Ivins fingerprint self-consistent"; "sparse L1 regularization point mass GRACE source localization"; "topological data analysis persistent homology satellite gravity" (no gravity use of TDA found).

---

## 7. Efficient numerics: fast SH transforms (FFT/GPU/WASM), streaming algorithms, full-covariance uncertainty propagation

### Takeaway
At degree 90, transforms are not the bottleneck. The modern fast transforms (SHTns, ducc0, FastTransforms) target degrees of thousands. The real bottlenecks are:
- covariance storage: 8,281² ≈ 68.6 M entries ≈ 549 MB per month in float64;
- Monte Carlo nulls for the scan statistics (about 10⁵–10⁶ syntheses).

Recommended engineering, in order:
1. SH analysis (quadrature) and rotation, both missing and needed by S3, F1b, F2 and F3.
2. A batched, cache-friendly synthesis for Monte Carlo work.
3. A "diagonal + low-rank" error-covariance model estimated from residuals, as a cheap substitute for full covariances.
4. Coefficients-to-browser WASM synthesis for the viewer.

### Cited Findings
- **SHTns** (Schaeffer): Gauss–Legendre fast SHT with explicit vectorisation and on-the-fly associated Legendre functions. O(N³), tested to N = 1023, and reported to perform "much better than any third party implementation" ([arXiv:1202.6522](https://arxiv.org/abs/1202.6522)).
- **ducc0 `sht`:** derived from libsharp and "significantly enhanced". It adds arbitrary grids via non-uniform FFTs, exact analysis on Clenshaw–Curtis/Fejér/McEwen–Wiaux grids, and iterative analysis on arbitrary grids. It is up to 2× faster at high band limits with Ishioka (2018) recurrences, multithreaded and SIMD, and its rotation code comes from FastTransforms. C++17, GPLv2+; the FFT parts are also BSD ([PyPI ducc0](https://pypi.org/project/ducc0/)).
- **FastTransforms** (Slevinsky; MIT licence): `sph2fourier` converts SH to bivariate Fourier series (degree 1023 in about 0.027 s; degree 8191 in about 4.9 s on 18 threads), with relative errors of about 5e-16 to 5e-15. It also does `sph_rotation` ([GitHub](https://github.com/MikaelSlevinsky/FastTransforms); Slevinsky 2019, ACHA 47:585–606).
- **GPU SHT:**
  - Hupca et al. (2010) reported up to 18× over S2HAT on one core and 5.5× over four cores, with FFTs, not Legendre, limiting performance ([arXiv:1010.1260](https://ar5iv.arxiv.org/html/1010.1260)).
  - Szydlarski et al. (2011) covered GPUs and multicore CPUs ([arXiv:1106.0159](https://arxiv.org/pdf/1106.0159)).
  - The HPEC 2018 MonarchSHT abstract identifies the Legendre transform as the dominant cost ([HPEC](https://abstractbook.ieee-hpec.org/s14.htm)).
  - NVIDIA torch-harmonics provides differentiable GPU SHTs ([tutorial](https://nvidia.github.io/torch-harmonics/tutorials/train_spherical_neural_operator.html)).
- **Covariance at scale:**
  - Propagating full monthly covariances to Greenland mascons with optimal weighting cut random noise by 35–69% per drainage system. The propagated matrices are ill-conditioned and were inverted approximately via eigen-decomposition ([PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6952056/); [J. Geod. 2017](https://link.springer.com/article/10.1007/s00190-017-1063-5)).
  - GFZ calls monthly covariance propagation to TWS "computationally expensive" and uses a non-stationary anisotropic spatial covariance model instead ([GFZ](https://gfzpublic.gfz.de/rest/items/item_5009734_5/component/file_5010009/content)).
  - A Stuttgart thesis tested full versus block-diagonal versus diagonal covariances in Bayesian regularisation ([Uni Stuttgart](https://elib.uni-stuttgart.de/items/aeb375de-dde5-4833-97c9-41e6f1125b04/full)).
- **Full covariances are available for ITSG** (used by VADER); COST-G gives only empirical sigmas (see the header notes at the top of this file).

### Inferences
**N1. SH analysis + rotation in `geopotential`.** Novelty: *none* (infrastructure), but a hard prerequisite.
- Gauss–Legendre or Driscoll–Healy/Clenshaw–Curtis quadrature analysis, which is exact for band-limited fields. The existing Legendre recursion plus an FFT over longitude gives O(N³) at N = 90, which is trivial.
- Rotation via Wigner-d (Risbo or Trapani–Navaza recursions; background knowledge) or by porting FastTransforms' algorithm, which is MIT-licensed. ducc0 is GPL, so it is not portable into MIT/Apache crates.
- *Validate:* round trip synthesis → analysis to about 1e-14; rotation composition identities.

**N2. Diagonal + low-rank error covariance from residuals.** Novelty: *partial* (VADER uses formal full covariances; residual-based low-rank models were not found).
- *Model:* Σ_t ≈ D_t + U Λ Uᵀ, where D_t = diag(σ²_nm(t)) and U holds the leading EOFs of month-to-month *residuals* after the S1 signal model is removed. Stripes are coherent across months, so they show up in U.
- *Payoffs:*
  - Woodbury makes Wiener, GLS and matched-filter weights O(N k²) with k ≈ 20–50.
  - It fixes C1's diagonal-noise assumption and gives an anisotropic, VADER-like filter *without* full covariances.
- *Validate:* compare against ITSG-Grace2018 full monthly covariances (correlation of eigenvectors; filter output differences); variance calibration on withheld months.

**N3. Streaming pipeline.** Each new GRACE-FO month should trigger:
- one S1 Kalman step: O(N_coeff);
- one BOCPD update per tracked series;
- a rank-one incremental SVD update for EOFs and N2 (Brand-style incremental SVD; background knowledge, no source found because of the search budget);
- one C1 scan at the newest t₀.
Everything stays O(N k) per month, which suits an always-on monitoring service.

**N4. Batched synthesis for Monte Carlo.** Scan-statistic nulls need about 1,000 permutations × 259 epochs × a few radii, i.e. about 10⁶ degree-90 syntheses.
- Batch many fields per Legendre evaluation: compute P̄_nm(cos θ) once per latitude ring and apply it to a matrix of coefficient vectors (a BLAS-3 style GEMM), then batched FFTs. This uses the cache far better than repeated single syntheses.
- A GPU is optional: published GPU gains are about 5–18×, and the bottleneck shifts to FFTs.

**N5. WASM synthesis in the viewer.**
- The viewer currently ships precomputed zlib+base64 grids.
- Shipping coefficients instead costs 259 × 8,281 × 4 bytes ≈ 8.6 MB as float32, before compression.
- Synthesising in-browser from Rust compiled to WASM makes the filter radius, Wiener/needlet scale, flux arrows (F1) and detection maps (C1) interactive.
- At degree 90 a 1° grid is about 65 k points × 8,281 terms, so it needs a proper Legendre-ring + FFT synthesis, not point-by-point evaluation.

### Gaps
- No benchmarks of Rust SH libraries were gathered, and crates.io was not searched (search budget). Check before writing N1 from scratch.
- No sources were gathered for streaming/incremental SVD in geoscience or for randomized low-rank covariance in gravimetry (search budget exhausted at "ducc0", "FastTransforms", "Wasserstein verification", "LSC≡GP" and "streaming incremental SVD" queries; ducc0 and FastTransforms were then read directly).
- Search terms used: "SHTns GPU spherical harmonic transform performance"; "GRACE full variance-covariance propagation mass estimates efficient monthly covariance matrices".

---

## 8. Which candidates are most novel and most worth implementing first?

### Takeaway
Rank by novelty (no prior GRACE use found) × usefulness × fit to our data (259 COST-G months with diagonal calibrated sigmas, SH-native crates):
1. **C1 SH-domain GLRT/matched-filter scan** for abrupt mass events;
2. **F1 minimum-norm transport field** ("where water moves"), plus the H⁻¹ comparison metric;
3. **S1/G1 per-coefficient state-space ≡ exact spherical Matérn GP** for gap filling with uncertainties, plus G2 Swarm/SLR fusion;
4. **C2 BOCPD** streaming alerts;
5. **D1 optimized/BOP-DMD** across gaps;
6. **D2 LIM/POP**.

The infrastructure prerequisite is **N1 (SH analysis + rotation)**, and the cross-cutting enabler is **N2 (diagonal + low-rank error covariance)**.

### Cited Findings
- AR is the method to beat for gap filling (about 0.9–1.7 cm EWH), with the remove-restore protocol as the agreed test ([Lecomte et al. 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11585520/)).
- Destriping distorts coseismic signals, which motivates a filter-free detector ([EPS 2015](https://earth-planets-space.springeropen.com/articles/10.1186/s40623-015-0188-0); [Wang et al. 2012](https://geoweb.princeton.edu/people/simons/PDF/reprints/GRL-2012a.pdf)).
- Lateral transport inference so far is model-based (adjoint diffusion–advection) or proposed as a DMD sub-basin network ([EGU26-4049](https://meetingorganizer.copernicus.org/EGU26/EGU26-4049.html); [Uni Stuttgart](https://www.gis.uni-stuttgart.de/forschung/projekte/data-driven-modelling)).
- DMD on GRACE exists (Libero et al.), so only its irregular-sampling and uncertainty-quantified variants are new ([EGU24-18454](https://meetingorganizer.copernicus.org/EGU24/EGU24-18454.html); [arXiv:1704.02343](https://arxiv.org/abs/1704.02343); [arXiv:2107.10878](https://arxiv.org/abs/2107.10878)).

### Inferences
Shortlist table:

| # | Algorithm | Novelty for GRACE | Needs beyond current crates | Effort | Validation anchor |
|---|---|---|---|---|---|
| C1 | SH-domain GLRT/matched-filter + scan-statistic nulls | None found | Cap-load templates (Love numbers exist); N2 for correlated noise; N1 rotation for double-couple templates | Medium | Blind recovery of Sumatra 2004, Maule 2010, Tohoku 2011; injection detection curves |
| F1 | Minimum-norm transport field ψ = R²σ̇/(n(n+1)), F = ∇ψ; H⁻¹ metric | None found | None for the basic version (partials exist); N1 for masked or OT versions; ERA5 for the budget version | Low | Synthetic moving blob; model truth with known fluxes; Greenland/Amazon sanity |
| S1/G1 | Per-coefficient state-space smoother ≡ spherical Matérn GP with calibrated σ(t) | Partial (framing new) | ML/EM for hyperparameters | Low–medium | Lecomte remove-restore; posterior coverage |
| G2 | Kalman fusion of COST-G + Swarm + SLR in SH domain | Partial | SLR reader | Medium | 2015–16 withheld block vs AR/MSSA/Swarm |
| C2 | BOCPD with known observation variance | None found | None | Low | Synthetic steps; known events; false-alarm run length |
| D1 | Optimized/BOP-DMD on σ-whitened coefficients | Partial | Dense linear algebra (faer/nalgebra) | Medium | Remove-restore; annual-period CI |
| D2 | LIM/POP in EOF space | None found | Same as D1 | Low | 1–6-month forecast skill vs persistence; gap test |
| S2 | Needlets with σ-aware thresholding | None found | None | Low | Synthetic leakage vs Gaussian/DDK/Wiener |
| S3 | Region-adaptive Wiener from Slepian multitaper | Partial | N1 for arbitrary regions | Medium | Closed-loop synthetic |
| N2 | Diagonal + low-rank residual error covariance | Partial | None | Low–medium | vs ITSG-Grace2018 full covariances |
| K1/K2 | PCMCI, Liang–Kleeman on basins × climate indices | None found (application) | External climate data | Medium | Synthetic known graphs |
| D3 | Slow feature analysis | None found | None | Low | Synthetic slow driver |
| F3 | L1 sparse source localisation | Partial | N1 | Medium | Injected loads |
| — | Koopman, CCM, TE, RQA, CP/Tucker, neural fields, PINNs | None found / weak fit | Various | High or low value | Defer |

Suggested sequence:
1. N1 (analysis/rotation), then F1 (fast, visual, likely a short paper or viewer feature).
2. Then S1, building the residual stream that C2 and N2 consume.
3. Then N2, then C1 (the main paper candidate).
4. D1/D2 as an "operator-learning on GRACE coefficients" study benchmarked under the same Lecomte protocol as S1/G2.

Shared validation harness to build once:
- a remove-restore runner (Lecomte metric, 1° latitude-weighted RMS, bootstrap);
- a synthetic-truth generator (model field + noise drawn from COST-G sigmas + correlated stripe noise from N2's U);
- an event catalogue (coseismic epochs and locations from the cited detections; Ganges–Brahmaputra 2007 flood; the 2015–16 ENSO regime shift);
- a cross-centre comparison (COST-G vs CSR/GFZ/JPL) using both RMS and the F1 H⁻¹ metric.

### Gaps
- Every "none found" verdict rests on general-web searches (59 queries; terms listed per section) and should be re-checked in Google Scholar / Web of Science before any publication claims priority. The ones that matter most are C1 (GRACE + likelihood ratio, template fitting, Chao, Han), F1 (GRACE + flux potential, Poisson, Helmholtz, H⁻¹, Wasserstein) and S1 (per-coefficient Kalman + calibrated errors).
- Not checked: whether COST-G sigmas are well calibrated month by month (S1 estimates the scale factor s_n, which answers this as a side result), and whether ITSG full covariances can be legally redistributed in test fixtures (check the licence before committing any; see the project's rule on third-party data).
