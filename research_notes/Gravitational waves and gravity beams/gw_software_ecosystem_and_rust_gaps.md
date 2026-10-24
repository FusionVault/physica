# Gravitational-wave data-analysis software ecosystem (2026) and gaps a Rust crate could fill

Research date: 2026-10-24. Versions and dates come from live queries of the PyPI JSON API, the crates.io API, the GitHub/GitLab APIs and the GWOSC event API, made on that day. Formulas were checked against the primary papers (PDF text extracted locally) or against the LALSuite source code.

## 1. What the main tools do, how they are licensed, and how well they are maintained in 2026

### Takeaway
The Python stack (LALSuite, PyCBC, GWpy, Bilby, the PTA and LISA tools) is mature, active and mostly copyleft (GPL) at the core, with permissive (MIT) tools on top. **Rust is no longer empty.** Since mid-2026 several crates from people inside IGWN have appeared:
- `waveformrs` (Colm Talbot, GPL-3.0, published 2026-10-22): TaylorF2, IMRPhenomD/X*, NRTidal and NRSur7dq4, all validated against LALSimulation.
- `bilby_rust`: antenna response and GMST.
- `distmarg`: distance-marginalised likelihood.
- `igwn-ligolw`: LIGO_LW XML.
- `libgrav`: binary parameter conversions plus a WASM binding; not on crates.io.

A new crate therefore should **not** re-implement waveforms first. The parts still uncovered are sensitivity curves, the SNR/horizon toolkit, SGWB/PTA tooling, file I/O (GWF, GWOSC HDF5) and the high-frequency band.

### Cited findings

**Python: core LIGO-Virgo-KAGRA (LVK) stack**
- **LALSuite** 7.26.16 (PyPI, uploaded 2026-10-02), GPL-2.0-or-later, hosted at git.ligo.org/lscsoft/lalsuite. It is the C reference library (LALSimulation waveforms, detector geometry) — [PyPI](https://pypi.org/project/lalsuite/), [repo](https://git.ligo.org/lscsoft/lalsuite)
- **PyCBC** 2.11.0 (2026-04-15), GPL-3.0, needs Python ≥ 3.11. GitHub gwastro/pycbc has 392 stars and was last pushed 2026-10-23. It describes itself as the "core package to analyze gravitational-wave data, find signals, and study their parameters… used in the first direct detection (GW150914)" — [PyPI](https://pypi.org/project/PyCBC/), [GitHub](https://github.com/gwastro/pycbc)
- **GWpy** 4.1.0 (2026-10-13), GPL-3.0-or-later. Development is on GitLab; the GitHub mirror has 383 stars and was pushed 2026-08-28 — [PyPI](https://pypi.org/project/gwpy/), [GitHub mirror](https://github.com/gwpy/gwpy)
- **Bilby** 2.8.2 (2026-08-10), MIT. bilby-dev/bilby was pushed 2026-10-20 (129 stars). It is "a unified framework for stochastic sampling packages and gravitational-wave inference" — [PyPI](https://pypi.org/project/bilby/), [GitHub](https://github.com/bilby-dev/bilby). Its compiled helper **bilby.cython** 0.5.4 (2026-04-02) is MIT — [PyPI](https://pypi.org/project/bilby.cython/)
- **gwsurrogate** 1.2.0 (2026-09-10). PyPI says MIT; GitHub reports the licence as "NOASSERTION". It gives access to NR surrogate waveforms — [PyPI](https://pypi.org/project/gwsurrogate/), [GitHub](https://github.com/sxs-collaboration/gwsurrogate). **sxs** 2026.0.7 (2026-09-30) is MIT — [PyPI](https://pypi.org/project/sxs/)
- **pyseobnr** 0.3.7 (2026-09-21), GPL-3.0-or-later — [PyPI](https://pypi.org/project/pyseobnr/)
- JAX-based tools: **ripple** (rippleGW 0.4.1, 2026-10-14, MIT, "differentiable gravitational waveforms in JAX") and **JimGW** 0.5.2 (2026-10-14, MIT) — [ripple GitHub](https://github.com/GW-JAX-Team/ripple), [PyPI JimGW](https://pypi.org/project/JimGW/)
- Other active packages:
  - **gwosc** client 0.8.3 (2026-06-25, MIT) — [PyPI](https://pypi.org/project/gwosc/)
  - **gwpopulation** 1.4.0 (2026-10-23, MIT) — [PyPI](https://pypi.org/project/gwpopulation/)
  - **dingo-gw** 0.10.3 (2026-10-12, MIT) — [PyPI](https://pypi.org/project/dingo-gw/)
  - **pesummary** 1.8.1 (2026-09-21, MIT) — [PyPI](https://pypi.org/project/pesummary/)
  - **ligo.skymap** 2.5.4 (2026-06-06, GPL-3.0-or-later) — [PyPI](https://pypi.org/project/ligo.skymap/)
  - **pygwinc** (PyPI name `gwinc`) 0.6.2 (2024-10-15, Unlicense), the detector noise-budget code — [PyPI](https://pypi.org/project/gwinc/)

**Python: detector-network Fisher and forecast tools**
- **GWFish**: the PyPI release is 0.0.1 from 2022-01-19 and is stale. GitHub janosch314/GWFish is active (pushed 2026-04-20, 44 stars), described as "Simulation of detector networks with Fisher-matrix PE". GitHub detects no licence and there is no LICENSE file in the repository root; only the old PyPI classifier says MIT — [PyPI](https://pypi.org/project/GWFish/), [GitHub](https://github.com/janosch314/GWFish)
- **gwbench** 0.8.9 (2025-10-17). GitLab sborhanian/gwbench is GPL-2.0+, last active 2026-07-08 — [PyPI](https://pypi.org/project/gwbench/), [GitLab](https://gitlab.com/sborhanian/gwbench)
- **gwfast** 1.1.2 (2024-11-04), GPL-3.0. GitHub was last pushed 2025-07-09, so it is slowing. It ships copies of the CE and ET ASD files — [PyPI](https://pypi.org/project/gwfast/), [GitHub](https://github.com/CosmoStatGW/gwfast), [psds readme](https://gwfast.readthedocs.io/en/latest/psds_readme_link.html)

**Python: continuous waves and stochastic background**
- **PyFstat** 2.4.0 (2026-08-19), MIT, pushed 2026-10-23: "a python package for gravitational wave analysis with the F-statistic" — [PyPI](https://pypi.org/project/PyFstat/), [GitHub](https://github.com/PyFstat/PyFstat)
- **pygwb** 1.5.1 (2024-07-10), MIT, git.ligo.org/pygwb/pygwb. No PyPI release in more than 2 years — [PyPI](https://pypi.org/project/pygwb/)

**Python: pulsar timing arrays (PTA)**
- **enterprise-pulsar** 3.5.0 (2026-07-14), MIT, Python 3.10 to below 3.15, pushed 2026-07-15. It is a PTA "noise analysis, gravitational-wave searches, and timing model analysis" code — [PyPI](https://pypi.org/project/enterprise-pulsar/), [GitHub](https://github.com/nanograv/enterprise)
- **enterprise_extensions** 3.0.3 (2025-06-05), MIT — [PyPI](https://pypi.org/project/enterprise-extensions/)
- **PTArcade** 1.1.5 (2025-06-27), MIT, Python ≥ 3.9 and below 3.13. GitHub pushed 2026-10-20. It is "a simple wrapper of ENTERPRISE and ceffyl" for new-physics searches — [PyPI](https://pypi.org/project/ptarcade/), [GitHub](https://github.com/andrea-mitridate/PTArcade)
- **discovery** (NANOGrav), MIT, pushed 2026-10-13: "Next-generation pulsar-timing-array data analysis". It is not on PyPI as `discovery-pulsar` — [GitHub](https://github.com/nanograv/discovery)
- Other PTA tools:
  - **PINT** (pint-pulsar 1.1.8, 2026-10-08, BSD) — [PyPI](https://pypi.org/project/pint-pulsar/)
  - **holodeck-gw** 1.6 (2026-09-29, MIT) — [PyPI](https://pypi.org/project/holodeck-gw/)
  - **ceffyl** 1.41.2 (2025-12-04) — [PyPI](https://pypi.org/project/ceffyl/)
  - **hasasia** 1.2.3 (2021-09-15, MIT; PTA sensitivity curves, stale) — [PyPI](https://pypi.org/project/hasasia/)

**Python: LISA**
- **lisaanalysistools** 1.2.8 (2026-02-10), Apache-2.0. Its repository lisa-analysis-tools/lisa-analysis-tools was pushed 2026-10-24. Note that the PyPI name `lisatools` is a separate, older package (0.5.1, 2023-01-24, MIT) — [PyPI lisaanalysistools](https://pypi.org/project/lisaanalysistools/), [GitHub](https://github.com/lisa-analysis-tools/lisa-analysis-tools), [PyPI lisatools](https://pypi.org/project/lisatools/)
- **BBHx** 1.2.3 (2025-10-03). The licence sources conflict: the PyPI classifier says MIT, while the GitHub licence file is GPL-3.0 — [PyPI](https://pypi.org/project/bbhx/), [GitHub](https://github.com/lisa-analysis-tools/BBHx)
- **fastemriwaveforms** 2.2.0 (2026-10-08), MIT — [PyPI](https://pypi.org/project/fastemriwaveforms/)
- LISA simulation chain from the LISA Consortium (gitlab.in2p3.fr), all BSD-3-Clause and active in 2026:
  - **lisaorbits** 3.1.0 (2026-07-02)
  - **lisagwresponse** 3.0.2 (2026-01-19)
  - **pytdi** 2.2.1 (2026-08-05)
  - **lisainstrument** 2.3.0 (2026-06-28)
  - Sources: [lisaorbits](https://pypi.org/project/lisaorbits/), [pytdi](https://pypi.org/project/pytdi/), [lisainstrument](https://pypi.org/project/lisainstrument/), [lisagwresponse](https://pypi.org/project/lisagwresponse/)
- The LDC software has no PyPI package named `ldc`. The PyPI API returned "not found" — [PyPI query](https://pypi.org/pypi/ldc/json)

**Python: sensitivity-curve tools**
- **LEGWORK** 1.0.0 (2026-09-16), MIT, 45 stars, pushed 2026-09-16: "predictions about stellar-origin sources in space-based gravitational wave detectors" — [PyPI](https://pypi.org/project/legwork/), [GitHub](https://github.com/TeamLEGWORK/LEGWORK)
- **gwent** 0.5.0 (2023-07-31), MIT. GitHub was last pushed 2026-01-22. It is the "Gravitational Wave dEtector desigN Toolkit" — [PyPI](https://pypi.org/project/gwent/), [GitHub](https://github.com/ark0015/gwent)
- **GWplotter** (Moore, Cole & Berry), a web tool at http://rhcole.com/apps/GWplotter — [Moore, Cole & Berry 2015, arXiv:1408.0740](https://arxiv.org/abs/1408.0740)
- Thrane and Romano's power-law integrated (PI) curve code is at LIGO-P1300115 — [DCC](https://dcc.ligo.org/LIGO-P1300115/public) (link from the paper and a search snippet)

**Julia and C/C++**
- Julia:
  - **PostNewtonian.jl** (moble), MIT, pushed 2026-10-08, 19 stars: "orbital dynamics and waveforms for binary black-hole systems in the post-Newtonian approximation" — [GitHub](https://github.com/moble/PostNewtonian.jl)
  - **GWecc.jl** (PTA signals from eccentric supermassive black-hole binaries), GPL-3.0, last pushed 2024-06 — [GitHub](https://github.com/abhisrkckl/GWecc.jl)
  - **GWSC.jl**, sensitivity curves, mainly for PTAs — [GitHub](https://github.com/bingining/GWSC.jl)
  - **binado/GravitationalWaves.jl** and **AstroSGWB.jl** are new and have 0 stars — [GitHub](https://github.com/binado/GravitationalWaves.jl)
  - No mature, general Julia GW data-analysis stack turned up.
- C/C++:
  - LALSuite is C.
  - The **cWB** (coherent WaveBurst) pipeline library is C++, GPL-3.0, on GitLab gwburst/public/library, last active 2026-10-19 — [GitLab](https://gitlab.com/gwburst/public/library)
  - **PhaseTracer** (C++, GPL-3.0) computes GW spectra from first-order phase transitions — [GitHub](https://github.com/PhaseTracer/PhaseTracer)
  - **synthlisa** (C++ LISA simulator) has been stale since 2017 — [GitHub](https://github.com/vallis/synthlisa)

**Rust: every GW-relevant crate and repository found (crates.io searched for "gravitational wave", "ligo", "pulsar timing", "lisa", "post-newtonian", "inspiral", "bilby", "gwf", "healpix" and others)**
- **waveformrs** 0.2.0 (crates.io and PyPI, 2026-10-22), GPL-3.0-or-later, by ColmTalbot.
  - Models: TaylorF2, IMRPhenomD, IMRPhenomXAS/XHM/XP/XPHM/XPNR, NRTidalv2 and NRTidalv3 variants, and the time-domain surrogate NRSur7dq4.
  - "Every model is checked against LALSimulation in the test suite."
  - Features: batches across threads, fixed frequency grids and Python bindings. Tested on macOS and Linux.
  - Licence: GPL because "much of it is translated from LALSuite and the GNU Scientific Library".
  - Dependencies are only num-complex, once_cell and an optional rustfft.
  - Sources: [crates.io](https://crates.io/crates/waveformrs), [GitHub](https://github.com/ColmTalbot/waveformrs)
- **bilby_rust** 0.1.2 (2026-10-21; first release 2024-08-21). crates.io says MIT; there is no licence file on GitHub. It implements Bilby geometry in Rust via PyO3:
  - `antenna_response`, `antenna_response_all_model` (all polarisation modes) and `time_delay_from_geocenter_vectorized`
  - `frequency_dependent_detector_tensor` and `greenwich_mean_sidereal_time`
  - Sources: [crates.io](https://crates.io/crates/bilby_rust), [GitHub](https://github.com/ColmTalbot/bilby_rust)
- **distmarg** 0.1.0 (2026-10-21), MIT, by ColmTalbot. It computes a distance-marginalised GW likelihood in the (ρ_mf, v*) basis. A pure-Rust implementation and a separate array-API Python implementation share no code and are cross-checked through mpmath "oracle" fixtures — [crates.io](https://crates.io/crates/distmarg), [GitHub](https://github.com/ColmTalbot/distance_marginalization)
- **igwn-ligolw** 0.1.0 (2026-05-23), MIT, by M. Coughlin: a Rust-native reader and writer for LIGO_LW XML — [crates.io](https://crates.io/crates/igwn-ligolw), [docs.rs](https://docs.rs/igwn-ligolw)
- **arrakis-schema** 0.5.0 (2026-07-27), LGPL-3.0-or-later: schemas for the IGWN "Arrakis" API — [crates.io](https://crates.io/crates/arrakis-schema)
- **ligo_hires_gps_time** 0.7.0 (2025-10-28), GPL-3.0-or-later, from LIGO CDS: a high-resolution GPS-epoch time library — [crates.io](https://crates.io/crates/ligo_hires_gps_time)
- **gwrs** 0.1.3 (2025-07-06), MIT OR Apache-2.0, 1714 downloads, 0 stars. The repository holds only basic TimeSeries, Segments, Channel and array types, with no commits since 2025-07, so it is effectively dormant — [crates.io](https://crates.io/crates/gwrs), [GitHub](https://github.com/open-physics/gwrs)
- **psrutils** 0.2.7 (2025-09-05), MIT OR Apache-2.0: parsers for tempo2 `.par` and `.tim` files plus an MJD type — [crates.io](https://crates.io/crates/psrutils), [GitHub](https://github.com/SGullin/psrutils)
- **libgrav** (Daniel Williams, GitHub transientlunatic/libgrav), Apache-2.0 on GitHub and MIT OR Apache-2.0 in Cargo.toml, last pushed 2026-10-07.
  - It calls itself "the spiritual successor to the LALSuite library".
  - Current content: only binary parameter conversions (Mtot, q, η, Mc, χ_eff, χ_p), constants and `uom` unit safety.
  - It claims bindings for Python, WASM (`grav-wasm` on npm), Julia, C, Fortran, Go, R and MATLAB.
  - The `grav` crate is **not on crates.io**, and a query of the npm registry for `grav-wasm` returned no published package.
  - Sources: [GitHub](https://github.com/transientlunatic/libgrav), [crates.io lookup](https://crates.io/api/v1/crates/grav)
- Hobby and small repositories:
  - **rustylisa** (LISA waveform generator, MIT, stale since 2024-11) — [GitHub](https://github.com/William-Metz/rustylisa)
  - **gw_fopt** (GW spectrum from first-order phase transitions, Apache-2.0, 2026-08) — [GitHub](https://github.com/hoquoctrung98/gw_fopt)
- Related general crates exist: **healpix** / **cdshealpix** (HEALPix tessellation) — [crates.io cdshealpix](https://crates.io/crates/cdshealpix)

### Inferences
- The Rust GW space is moving very fast in late 2026, and largely through one LVK developer (Colm Talbot: waveformrs, bilby_rust, distmarg) plus Coughlin (igwn-ligolw) and Williams (libgrav). Any new crate should coordinate with these or interoperate with them rather than compete.
- Licensing matters. A crate translated from LALSuite must be GPL; waveformrs states this explicitly. A crate written cleanly from the published papers (formulas below) can be MIT/Apache, which suits WASM, embedded use and mixed-licence projects. This is the main remaining reason to write a permissive TaylorF2 rather than depend on waveformrs.
- The weak spots in Python are in maintenance, not capability:
  - pygwb has no release since 2024-07.
  - gwent's last release was 2023 and hasasia's was 2021.
  - GWFish's PyPI package is a 2022 placeholder with no licence file.
  - gwfast has slowed.
  - PTArcade caps Python below 3.13.

  Sensitivity-curve and forecast tooling is therefore the most fragmented and least maintained layer.

### Gaps
- I did not check maintainers' statements on long-term support (for example LALSuite's migration plans), or whether waveformrs builds for `wasm32-unknown-unknown`. Its dependency list suggests it could, but it is only tested on macOS and Linux.
- I did not verify GstLAL, SPIIR or MBTA (other C/C++ search pipelines) in this session.
- I found no Rust GWF (frame-file) reader. The standard implementations, FrameL and frameCPP, are C/C++, but I did not verify their 2026 versions.

## 2. Open data: GWOSC, PTA and LISA releases and their formats

### Takeaway
GWOSC now hosts O1 to O4b. O4a (GWTC-4.0) was released on 2025-08-26 and O4b (GWTC-5.0, with Virgo) on 2026-05-26. All strain data is CC BY 4.0, in GWF and HDF5, at 16384 Hz and 4096 Hz, in fixed 4096-s files with 1 Hz DQ and injection bitmasks. The cumulative GWTC list in the event API holds 391 entries. On the PTA side, NANOGrav 15-yr and PPTA DR3 are public, and IPTA DR3 is reported as released in a 2026 proceeding but without a data paper I could find. For LISA, Sangria-HM (2024-11) is the latest confirmed public challenge. Mojito Light (4 Hz, 2 yr) was planned for the end of 2025, and its public release is not confirmed.

### Cited findings

**GWOSC: O4a and O4b**
- O4a ran from 2023-05-24 to 2024-01-16. The data release and GWTC-4.0 were published 2025-08-26 — [GWOSC search summary / O4a page](https://gwosc.org/O4/O4a), [EGO news](https://www.ego-gw.it/blog/2025/08/26/gwtc-4-0-updated-gravitational-wave-catalog-released/)
- O4b dates and detectors:
  - The run lasted from 10 Apr 2024 15:00 UTC (GPS 1396796418) to 28 Jan 2025 17:00 UTC (GPS 1422118818).
  - The release starts earlier, at 6 Apr 2024 05:37:12 UTC (GPS 1396417050), to include data used for SN2024ggi.
  - Detectors: H1, L1 and V1.
  - Sample rates: "the original 16 KHz and the downsampled 4KHz".
  - Detections are in GWTC-5.0.
  - Data DOI: 10.7935/8emv-ag54. First release May 26, 2026. "Data released under a CC BY 4.0 License".
  - Source channels: H1/L1:DCS-CALIB_STRAIN_CLEAN_AR01 and V1:Hrec_hoftRepro1AR_16384Hz.
  - O4 is split into O4a, O4b and O4c.
  - Source: [GWOSC O4b](https://gwosc.org/O4/O4b)
- The O4 technical details:
  - Data come as GWF and HDF5, at both 16384 Hz and 4096 Hz.
  - Downsampling uses `scipy.signal.decimate`.
  - The 4 kHz data are corrupted near Nyquist by the anti-aliasing filter. Use the 16 kHz data at about 1700 Hz and above.
  - Calibration is valid from 10 Hz to 5000 Hz. DQ is less well studied above 2 kHz.
  - Files have a fixed duration of 4096 s with fixed boundaries.
  - DQ and hardware-injection flags are "summarized in 1 Hz vectors" (bitmasks) in both formats.
  - GWF channels are `{ifo}:GWOSC-4KHZ_R1_STRAIN`, `…_DQMASK` and `…_INJMASK` (16KHZ equivalents exist).
  - "The HDF5 group names for O4 are the same as used on O3", plus a meta attribute holding the GWF channel name.
  - Per-file metadata includes a BNS range computed with gwpy `sensemon_range` from 10 Hz.
  - Source: [GWOSC O4 technical details](https://gwosc.org/O4/o4_details)
- GWTC-4.0 contents:
  - It includes O4a transients with p_astro > 0.5, and PE for events with FAR < 1/yr.
  - Event strain is available as GWF, HDF or TXT from the Event Portal.
  - PE samples are on Zenodo, using PESummary parameter names.
  - The initial PE has "a known issue in the likelihood used for inference" (Methods paper §5.10). "(reweighted)" and "(rerun)" results are provided, and the default PE is the rerun.
  - GW230630_070659 was later marked as having evidence of instrumental origin.
  - The page lists the Open Data paper as "ApJS 1004 2329 (2026)", quoted as shown.
  - Source: [GWOSC GWTC-4.0](https://gwosc.org/GWTC-4.0/)
- Event-API counts on 2026-10-24: GWTC-1-confident 11, GWTC-2.1-confident 54, GWTC-3-confident 35, GWTC-4.0 129, GWTC-5.0 161, cumulative "GWTC" 391 — [GWOSC event API](https://gwosc.org/eventapi/json/GWTC/). A seminar listing gives "161 new events… bringing the total to 390" — [CFT seminar](https://www.cft.edu.pl/nauka/seminaria/390-and-counting-highlights-from-the-ligo-virgo-kagra-recent-o4b-run-data-release/). The difference of 1 is probably the GW230630 instrumental flag or another counting convention; I did not check.
- O4c: a LIGO memo (revised June 2025) set the planned end of O4 as 18 Nov 2025 and said an extension would not be part of the O4c release. I found no confirmed O4c release date — [LIGO memo M2300059](https://dcc-backup.ligo.org/public/0186/M2300059/005/20240911%20Memo_%20Public%20Data%20Releases%20for%20O4%20and%20O5.pdf)
- Per-event strain files for GW150914 (GWTC-1) are reached from the event JSON — [GWOSC event JSON](https://gwosc.org/eventapi/json/GWTC-1-confident/GW150914/v3)

**PTAs**
- NANOGrav 15-yr: 67 pulsars, published in ApJL as part of the 15-year focus issue, DOI 10.3847/2041-8213/acdac6 — [arXiv:2306.16213](https://arxiv.org/abs/2306.16213)
- PPTA DR3: 32 pulsars from the 64-m Parkes "Murriyang" telescope, spanning up to 18 yr — [arXiv:2306.16230](https://arxiv.org/abs/2306.16230)
- IPTA DR3: an APS April 2025 abstract describes it as "in preparation" with 121 pulsars and CHIME, LOFAR and NenuFAR data — [APS 2025 abstract](https://archive.aps.org/smt/2025/apr-b22/1). A 2026 *J. Phys. Conf. Ser.* proceeding (via search snippet) refers to its release, with 130 pulsars, baselines up to 27 yr and coverage of 30 to 5107 MHz. **The pulsar counts conflict and I found no IPTA DR3 data paper.**
- Formats: tempo2 `.par` (timing model) and `.tim` (TOAs) files. These can be parsed in Rust with psrutils — [psrutils README](https://github.com/SGullin/psrutils)

**LISA**
- "Challenge 2a HM (Sangria-HM)" was released 15 Nov 2024 to "introduce new researchers to LISA data analysis", with entries due 1 Mar 2025. The LDC site does not list data formats, cadence or a licence — [LISA Data Challenges](https://lisa-ldc.in2p3.fr/)
- Mojito, LISA's first Common Dataset:
  - It was planned "by the end of 2025", with a "Mojito Light release in dec 2025".
  - A May 2026 DDPC talk describes it as 4 Hz data for 2 years with realistic TDI and ESA orbits, as a "LEGO-type dataset with 37 bricks".
  - The search did not confirm a public release.
  - Sources: [DDPC LISA France, May 2026](https://indico.in2p3.fr/event/38944/contributions/172717/attachments/101391/157416/DDPC_LISAFrance_260504.pdf), [LISA DA update 2025](https://indico.in2p3.fr/event/36513/contributions/162592/attachments/97006/149166/GroupeAnalyseDeDonnees_AG-GdR-OG_251014.pdf)

**Detector noise-curve data files (for validation)**
- Cosmic Explorer: CE-T2000017-v8 "Cosmic Explorer Strain Sensitivity".
  - Files: `cosmic_explorer_strain.txt` (40 km baseline), `cosmic_explorer_20km_strain.txt` (20 km, compact-binary tuned), `cosmic_explorer_20km_pm_strain.txt` (20 km, post-merger tuned) and `cosmic_explorer_40km_lf_strain.txt` (40 km, low-frequency tuned).
  - The strain curves are for a source 15° off normal incidence.
  - Source: [CE DCC T2000017](https://dcc.cosmicexplorer.org/T2000017/public)
- Einstein Telescope:
  - ET-D (Hild et al. 2011) — [ET TDS r=14065](https://apps.et-gw.eu/tds/?content=3&r=14065)
  - The 2023 design-comparison ASDs (10, 15 and 20 km; HF-only and HF+LF cryogenic) are listed by gwfast, which points to ET TDS r=18213 — [gwfast psds readme](https://gwfast.readthedocs.io/en/latest/psds_readme_link.html), [Branchesi et al. 2023, arXiv:2303.15923](https://arxiv.org/abs/2303.15923)
- aLIGO O4 "aligo_O4high.txt" is in LIGO-T2000012 (used by Chen et al.) — [Chen et al., arXiv:1709.08079](https://arxiv.org/abs/1709.08079), [LIGO-T2000012](https://dcc.ligo.org/LIGO-T2000012/public)

### Inferences
- A Rust crate that reads GWOSC HDF5 needs only the generic `hdf5` crate plus the documented GWOSC layout. It also has to decode the 1 Hz DQ and INJ bitmasks, which is a small, well-specified, testable task.
- Reading GWF natively in Rust would mean writing a frame-format parser (with its compression) from scratch, a larger job. HDF5 should come first.

### Gaps
- I did not fetch the exact HDF5 dataset paths, such as `strain/Strain`, `meta/GPSstart` and `quality/simple/DQmask`. GWOSC says the O4 names "are the same as used on O3", but the O3 page was not fetched. Check it on the O3 technical-details page before coding.
- I could not confirm the LDC data formats (HDF5 is likely) or the Mojito licence and download location.
- I did not retrieve the NANOGrav 15-yr data-set timespan (commonly quoted as about 16 years) or the free-γ fit from a fetched source.

## 3. Algorithms and exact formulas for a validated core

### Takeaway
Every formula below comes from a primary source (paper PDF text or LALSuite source) and can be implemented in roughly 1–2k lines of Rust with no external numerical dependencies beyond an FFT and quadrature. Each item has a reference implementation to compare against (LAL/PyCBC/Bilby/hasasia/legwork/Thrane–Romano code).

### Cited findings

**(a) TaylorF2: 3.5PN non-spinning phase and Newtonian amplitude (stationary-phase approximation)**
- Sathyaprakash and Schutz 2009 (LRR 12:2), Eqs. 127–129, with ν = η the symmetric mass ratio and units G = c = 1:
  - H(f) = A f^(−7/6) exp[i Ψ(f) + iπ/4]
  - A = (C/(D π^(2/3))) √(5ν/24) M^(5/6)
  - Ψ(f) = 2π f t_C + Φ_C + (3/(128 ν)) Σ_k α_k (π M f)^((k−5)/3)

  Source: [Sathyaprakash & Schutz 2009, arXiv:0903.0338](https://arxiv.org/abs/0903.0338)
- The coefficients α_k (same source, Eq. 129):
  - α0 = 1
  - α1 = 0
  - α2 = 3715/756 + 55ν/9
  - α3 = −16π
  - α4 = 15293365/508032 + 27145ν/504 + 3085ν²/72
  - α5 = π(38645/756 − 65ν/9)[1 + ln(6^(3/2) π M f)]
  - α6 = 11583231236531/4694215680 − 640π²/3 − (6848/21)γ_E + (−15737765635/3048192 + 2255π²/12)ν + 76055ν²/1728 − 127825ν³/1296 − (6848/63) ln(64 π M f)
  - α7 = π(77096675/254016 + 378515ν/1512 − 74045ν²/756)

  Source: [arXiv:0903.0338](https://arxiv.org/abs/0903.0338)
- The same coefficients in **LALSimulation's code**, `XLALSimInspiralTaylorF2Phasing_*` in `LALSimInspiralPNCoefficients.c`:
  - 2PN: `5.*(74.3/8.4 + 11.*eta)/9.`
  - 3PN: `-16π`
  - 4PN: `5.*(3058.673/7.056 + 5429./7.*eta + 617.*eta*eta)/72.`
  - 5PN: `5./9.*(772.9/8.4-13.*eta)*π`, with a log coefficient `5./3.*(772.9/8.4-13.*eta)*π` multiplying ln v
  - 6PN: `11583.231236531/4.694215680 - 640/3 π² - 684.8/2.1 γ + η(-15737.765635/3.048192 + 225.5/1.2 π²) + 76.055/1.728 η² - 127.825/1.296 η³ + (-684.8/2.1)·ln 4`, with a log coefficient `-684.8/2.1` multiplying ln v
  - 7PN: `π(770.96675/2.54016 + 378.515/1.512 η - 740.45/7.56 η²)`
  - Prefactor `pfaN = 3/(128 η)`, with the phase series divided by v⁵.
  - Spin-orbit, spin-spin, quadrupole-monopole and tidal (5PN to 7.5PN) terms are added on top.
  - Source: [LALSimInspiralPNCoefficients.c](https://git.ligo.org/lscsoft/lalsuite/-/raw/master/lalsimulation/lib/LALSimInspiralPNCoefficients.c)
- LAL's TaylorF2 conventions (`LALSimInspiralTaylorF2.c`):
  - `phasing += shft*f - 2*phi_ref - ref_phasing`
  - `amp0 = -4 m1 m2 / r · MRSUN · MTSUN · sqrt(π/12)`
  - `amp = amp0·sqrt(-dEnergy/flux)·v`
  - `h̃ = amp·[cos(phasing − π/4) − i sin(phasing − π/4)]`, which is a **negative-exponent** convention, opposite to S&S Eq. 127's `+iΨ`.
  - The default end frequency is the Schwarzschild ISCO, `vISCO = 1/√6`.
  - Source: [LALSimInspiralTaylorF2.c](https://git.ligo.org/lscsoft/lalsuite/-/raw/master/lalsimulation/lib/LALSimInspiralTaylorF2.c)

**(b) Matched filter, inner product and SNR**
- The inner product is ⟨Ã|B̃⟩ = 4 Re ∫₀^∞ Ã*(f) B̃(f)/S_n(f) df (Eq. 14, after Finn 1992). The optimal (Wiener) filter is K̃ = h̃/S_n (Eq. 15). The optimal SNR is ρ² = ∫₀^∞ 4|h̃(f)|²/S_n(f) df (Eq. 16). S_n is the **one-sided** PSD; the two-sided PSD is S_n/2 (footnote) — [Moore, Cole & Berry 2015, arXiv:1408.0740](https://arxiv.org/abs/1408.0740)
- Characteristic strain is h_c(f)² = 4 f² |h̃(f)|² (Eq. 17). Noise amplitude is h_n(f)² = f S_n(f) (Eq. 18). Then ρ² = ∫ d(ln f) [h_c/h_n]² (Eq. 19) — [arXiv:1408.0740](https://arxiv.org/abs/1408.0740)

**(c) Antenna patterns**
- For an L-shaped (90°) interferometer, with sky angles (θ, φ) in the detector frame and polarisation angle ψ (S&S Eq. 57):
  - F+ = ½(1 + cos²θ) cos 2φ cos 2ψ − cos θ sin 2φ sin 2ψ
  - F× = ½(1 + cos²θ) cos 2φ sin 2ψ + cos θ sin 2φ cos 2ψ
  - The tensor form is F+ = d:e+, F× = d:e× (Eq. 56), and the maximum of either is 1.
  - For non-perpendicular arms, the detector-frame x and y are chosen so the arm bisector lies along the bisector of the coordinate axes.
  - The bar detector has F+ = sin²θ cos 2ψ (Eq. 58).
  - Source: [arXiv:0903.0338](https://arxiv.org/abs/0903.0338)
- The overlap reduction function for interferometers with opening angle δ is normalised as γ_IJ(f) = (5/sin²δ) Γ_IJ(f) (Eq. 16), which shows that the response scales with sin δ. The ET and BBO-type geometry uses δ = 60° — [Thrane & Romano 2013, arXiv:1310.5300](https://arxiv.org/abs/1310.5300)
- The ET 2023 design paper says the triangle geometry, the arm angle and the number of interferometers must be accounted for separately when its ASDs are used — [arXiv:2303.15923](https://arxiv.org/abs/2303.15923) (via search summary)

**(d) Horizon distance and range**
- The horizon d_h is the maximum distance at which a binary has SNR ≥ 8, placed overhead and face-on. The range R is the radius of the Euclidean sphere with the same volume as the redshifted sensitive volume. At low z, R ≈ d_h/2.264 (Finn & Chernoff 1993), historically called the "sensemon distance". f_p = 2.264 is independent of the noise curve and the masses. In FLRW cosmology, solve for the z at which a face-on, overhead binary of redshifted mass (1+z)M has SNR 8 — [Chen et al. 2021, arXiv:1709.08079](https://arxiv.org/abs/1709.08079)

**(e) Chirp mass from frequency evolution**
- From the leading-order (quadrupole) frequency evolution: M_c = (m1 m2)^(3/5)/(m1+m2)^(1/5) = (c³/G)[(5/96) π^(−8/3) f^(−11/3) ḟ]^(3/5). This is the GW150914 discovery paper's estimator, giving M_c ≈ 30 M☉ — [Abbott et al. 2016, PRL 116 061102](https://doi.org/10.1103/PhysRevLett.116.061102). *The formula is recalled from the PRL body. Only the abstract was re-fetched in this session.*

**(f) Stochastic background**
- Equations from Thrane & Romano 2013:
  - Ω_gw(f) = (2π²/(3H0²)) f³ S_h(f) (Eqs. 3–4)
  - h_c(f) = √(f S_h(f)) (Eq. 5), so Ω_gw = (2π²/(3H0²)) f² h_c² (Eq. 6)
  - Power laws: Ω_gw = Ω_β (f/f_ref)^β (Eq. 7) and h_c = A_α (f/f_ref)^α (Eq. 8), with Ω_β = (2π²/(3H0²)) f_ref² A_α² and β = 2α + 2 (Eq. 9)
  - Flat (inflationary) spectrum: β = 0, α = −1. Binary coalescences: β = 2/3, α = −2/3 (Eqs. 10–11).
  - f_ref is typically 1/yr for PTAs and 100 Hz for ground-based detectors.
  - Source: [arXiv:1310.5300](https://arxiv.org/abs/1310.5300)
- Cross-correlation SNR for an isotropic, unpolarised background (same paper):
  - ρ = √(2T) [∫ df Σ_{I<J} Γ_IJ²(f) S_h²(f)/(P_nI(f) P_nJ(f))]^(1/2) (Eqs. 21–22)
  - Effective network noise S_eff(f) = [Σ_{I<J} Γ_IJ²/(P_nI P_nJ)]^(−1/2) (Eq. 23)
  - For M identical, co-located detectors, S_eff = √(2/(M(M−1))) S_n with S_n = P_n/R (Eqs. 26–27)
  - Source: [arXiv:1310.5300](https://arxiv.org/abs/1310.5300)
- **Hellings–Downs** (Thrane & Romano Eq. 20, normalised so ζ_II = 1), with x = (1 − cos ψ_IJ)/2:
  - ζ_IJ = (3/2) x ln x − x/4 + 1/2 + δ_IJ/2
  - The PTA ORF in residual units is Γ_IJ(f) = ζ_IJ / (3(2πf)²) (Eq. 19)
  - Source: [arXiv:1310.5300](https://arxiv.org/abs/1310.5300)
- **Power-law integrated (PI) curve** (Thrane & Romano Eqs. 29–30):
  - Ω_β = (ρ/√(2T)) [∫_{fmin}^{fmax} df (f/f_ref)^(2β)/Ω_eff²(f)]^(−1/2), where Ω_eff is S_eff converted with Eq. 3
  - Ω_PI(f) = max_β [Ω_β (f/f_ref)^β]
  - The paper uses β ∈ {−8, −7, …, 8}, ρ = 1 and T = 1 yr; the choice of f_ref does not change the curve
  - A predicted power law tangent to the curve has ρ = 1, and generally ρ = Ω_β^pred/Ω_β
  - Source: [arXiv:1310.5300](https://arxiv.org/abs/1310.5300)

**(g) Analytic noise curves**
- **aLIGO** (S&S Table 1): S_h = S0 [x^(−4.14) − 5x^(−2) + 111(1 − x² + 0.5x⁴)/(1 + 0.5x²)], with x = f/215 Hz, S0 = 1.0×10⁻⁴⁹ Hz⁻¹ and lower cutoff f_s = 20 Hz — [arXiv:0903.0338](https://arxiv.org/abs/0903.0338)
- **Einstein Telescope**, ET-B-era fit with S&S's rounded coefficients:
  - S_h = S0 [x^α + a0 x^β + b0 (1 + b1x + b2x² + b3x³ + b4x⁴ + b5x⁵ + b6x⁶)/(1 + c1x + c2x² + c3x³ + c4x⁴)]
  - x = f/200 Hz, S0 = 1.5×10⁻⁵² Hz⁻¹, f_s = 10 Hz
  - α = −4.1, β = −0.69, a0 = 186, b0 = 233, b1 = 31, b2 = −65, b3 = 52, b4 = −42, b5 = 10, b6 = 12, c1 = 14, c2 = −37, c3 = 19, c4 = 27
  - Source: [arXiv:0903.0338](https://arxiv.org/abs/0903.0338)
- The same table has the old LISA instrumental fit, S0 = 9.2×10⁻⁴⁴, f0 = 1 mHz, (x/10)^(−4) + 173 + x². It is superseded by Robson et al. — [arXiv:0903.0338](https://arxiv.org/abs/0903.0338)
- **LISA**, Robson, Cornish & Liu 2019 (CQG 36 105011):
  - Full model: S_n(f) = (10/(3L²)) [P_OMS + 2(1 + cos²(f/f*)) P_acc/(2πf)⁴] (1 + (6/10)(f/f*)²) (Eq. 13)
  - Eq. 1 is the simpler form with 4P_acc/(2πf)⁴, plus S_c
  - L = 2.5 Gm and f* = c/(2πL) = 19.09 mHz
  - P_OMS = (1.5×10⁻¹¹ m)² [1 + (2 mHz/f)⁴] Hz⁻¹ (Eq. 10)
  - P_acc = (3×10⁻¹⁵ m s⁻²)² [1 + (0.4 mHz/f)²][1 + (f/8 mHz)⁴] Hz⁻¹ (Eq. 11)
  - Galactic confusion: S_c(f) = A f^(−7/3) e^(−f^α + β f sin(κf)) [1 + tanh(γ(f_k − f))] Hz⁻¹ (Eq. 14), with A = 9×10⁻⁴⁵. This is half the value in the paper's Ref. [8] because the curves are two-channel.
  - Fit parameters by mission duration:

    | Duration | α | β | κ | γ | f_k (Hz) |
    |---|---|---|---|---|---|
    | 6 mo | 0.133 | 243 | 482 | 917 | 0.00258 |
    | 1 yr | 0.171 | 292 | 1020 | 1680 | 0.00215 |
    | 2 yr | 0.165 | 299 | 611 | 1340 | 0.00173 |
    | 4 yr | 0.138 | −221 | 521 | 1680 | 0.00113 |

  - Source: [Robson, Cornish & Liu, arXiv:1803.01944](https://arxiv.org/abs/1803.01944)
- **Cosmic Explorer** has no analytic fit in the primary source. The official product is tabulated ASDs (CE-T2000017-v8, see §2) — [CE DCC](https://dcc.cosmicexplorer.org/T2000017/public)

### Inferences
- **The horizon formula follows from (a), (b) and (d).** Here ρ_th = 8 and the integral runs to f_ISCO = c³/(6^(3/2) π G M):

  d_h = (1/ρ_th) √(5/6) (c/π^(2/3)) (G M_c/c³)^(5/6) [∫_{fmin}^{f_ISCO} f^(−7/3)/S_n(f) df]^(1/2)

  It comes from the Newtonian TaylorF2 amplitude |h̃| = (1/D)√(5/24) π^(−2/3) c (G M_c/c³)^(5/6) f^(−7/6), for a face-on, overhead source. At cosmological distances replace M_c with (1+z)M_c and solve for z, as Chen et al. do.
- **The Newtonian chirp time** is τ0 = (5/256)(πf)^(−8/3)(G M_c/c³)^(−5/3). This is the integral of ḟ from (e). It is a cheap self-consistency test for the TaylorF2 t(f).
- **Time-maximised SNR** is computed with an inverse FFT: z(t) = 4 ∫ s̃(f) h̃*(f)/S_n(f) e^(2πift) df, and ρ(t) = |z(t)|/σ with σ² = ⟨h|h⟩. This is the standard FINDCHIRP construction (Allen et al., PRD 85 122006, arXiv:gr-qc/0509116; *citation from prior knowledge, not re-fetched*). PyCBC's `matched_filter` and `sigma` are the reference implementations.
- **ET triangle.** Each of the three 60° interferometers should have F_A^(i) = sin 60° × F_A^L evaluated in its own arm-bisector frame, with the frames rotated by 120° from each other. This follows from the bisector convention (S&S) and the sin δ scaling in Thrane & Romano Eq. 16. Validate against `lal.ComputeDetAMResponse` with LAL's ET detector definitions, or against bilby_rust.
- **PTA conversions.** The PTA timing-residual power-law index γ relates to α by γ = 3 − 2α, so α = −2/3 gives γ = 13/3. This matches NANOGrav's use of γ = 13/3 for an f^(−2/3) h_c spectrum.
- **Sign-convention traps for validation:**
  - S&S use e^(+iΨ+iπ/4); LAL uses e^(−i(Ψ−π/4)).
  - S_n is one-sided, and two-sided is half.
  - Robson's 10/3 factor already includes sky and polarisation averaging over two channels.
  - CE strain curves assume 15° off normal incidence.
  - S&S's ET coefficients are rounded, so results differ slightly from the unrounded fit.

### Gaps
- I did not fetch the H1–L1 interferometer overlap reduction function (Flanagan 1993; Allen & Romano 1999), which is needed for ground-based PI curves. It is not in the sourced material here.
- I did not get exact equation numbers or the paper's own spin and tidal coefficient forms from Buonanno et al. 2009 (arXiv:0907.0700). The LAL source is the operative reference.
- The LIGO "aLIGO design" curves (T1800044, aLIGOZeroDetHighPower / T0900288) were not fetched. Only their existence is known, and they are used through PyCBC/LAL.

## 4. Reference numbers for validation tests

### Takeaway
Published anchors exist for GW150914, horizon distances (Chen et al., aLIGO O4high curve) and NANOGrav 15-yr. Hellings–Downs values, LISA curve values, Ω_gw h² and chirp-time values can be derived exactly from the cited formulas. These derived numbers, computed here in double-precision Python, should be regenerated with the reference code (LAL/PyCBC/hasasia/legwork) and used as fixtures.

### Cited findings
- **GW150914 (PRL abstract):**
  - Matched-filter SNR 24.
  - FAR below 1 per 203,000 yr, significance above 5.1σ.
  - Source-frame masses 36 (+5/−4) and 29 (+4/−4) M☉, final mass 62 ± 4 M☉, radiated energy 3.0 ± 0.5 M☉c².
  - Distance 410 (+160/−180) Mpc, z = 0.09 (+0.03/−0.04).
  - Peak strain 1.0×10⁻²¹, frequency sweep 35 → 250 Hz.
  - Source: [arXiv:1602.03837 / PRL 116 061102](https://arxiv.org/abs/1602.03837)
- **GW150914 (GWTC-1, GWOSC event API v3):**
  - GPS 1126259462.4.
  - m1 = 35.6 (+4.7/−3.1) and m2 = 30.6 (+3.0/−4.4) M☉.
  - Source-frame chirp mass 28.6 (+1.7/−1.5) M☉.
  - Network matched-filter SNR 25.2.
  - D_L = 440 (+150/−170) Mpc, z = 0.09.
  - Final mass 63.1, χ_eff = −0.01.
  - Source: [GWOSC event API](https://gwosc.org/eventapi/json/GWTC-1-confident/GW150914/v3)
- **Horizon and range (Chen et al. Table 1, aLIGO O4high curve from T2000012, IMRPhenomD, SNR 8):**

  | Binary (source frame) | Horizon d_h | Range R |
  |---|---|---|
  | 1.4 + 1.4 M☉ | z = 0.10, D_L = 449.4 Mpc (comoving 410.4) | 179.9 Mpc |
  | 10 + 10 M☉ | D_L = 2617.7 Mpc (z = 0.46) | 772.6 Mpc |
  | 30 + 30 M☉ | D_L = 8014.8 Mpc (z = 1.14) | 1572.7 Mpc |

  - A GW150914-like source has horizon z = 1.21 (D_L = 8553 Mpc).
  - Median "reach" distances d^r_50 are 188 Mpc for BNS and 2276 Mpc for 30+30.
  - Cosmology: Ω_m = 0.3065, h = 0.679.
  - Web calculator: https://users.rcc.uchicago.edu/~dholz/gwc/
  - Source: [arXiv:1709.08079](https://arxiv.org/abs/1709.08079)
- **NANOGrav 15-yr:**
  - A = 2.4 (+0.7/−0.6)×10⁻¹⁵ at f = 1/yr for γ = 13/3 (median, 90% credible interval), from 67 pulsars.
  - Hellings–Downs versus a common uncorrelated process: Bayes factor 200–1000.
  - Significance: p ≈ 10⁻³ (about 3σ) from the Bayes factors, and p = 5×10⁻⁵ to 1.9×10⁻⁴ (about 3.5–4σ) from the frequentist statistic.
  - Source: [arXiv:2306.16213](https://arxiv.org/abs/2306.16213)
- Thrane & Romano's figure assumptions (fixtures for PI curves): T = 1 yr, ρ = 1, β from −8 to 8; a 0.25 Hz bin width is quoted as typical for LIGO — [arXiv:1310.5300](https://arxiv.org/abs/1310.5300)

### Inferences (numbers computed here from the formulas in §3; regenerate with the reference code before using them as fixtures)
- **Hellings–Downs**, from Thrane & Romano Eq. 20 for distinct pulsars:

  | Angle | ζ |
  |---|---|
  | 0⁺ | 0.5 |
  | 10° | 0.44250 |
  | 30° | 0.21163 |
  | 60° | −0.08236 |
  | 90° | −0.14486 |
  | 120° | −0.01114 |
  | 150° | 0.16971 |
  | 180° | 0.25 |

  The zeros are at 49.32° and 121.79°, and the minimum is −0.15190 at 82.48°.
- **NANOGrav in energy-density units.** With Ω_gw(f) = (2π²/(3H0²)) f² h_c², A = 2.4×10⁻¹⁵ and f = 1/yr (Julian year), Ω_gw(f_yr) h² ≈ 3.6×10⁻⁹. The 90% band A = 1.8–3.1×10⁻¹⁵ gives 2.0–6.0×10⁻⁹.
- **GW150914 detector-frame chirp mass:** (1+z)·M_c,src = 1.09 × 28.6 ≈ 31.2 M☉. M_c computed from the median m1 and m2 is 28.72 M☉, not 28.6, because medians do not commute with the nonlinear chirp-mass function. A test must not expect exact agreement here.
- **BNS 1.4+1.4 M☉:**
  - M_c = 1.21877 M☉
  - ḟ(100 Hz) = 17.364 Hz/s, and inverting via (e) recovers M_c exactly
  - τ0(20 Hz) = 157.9 s and τ0(10 Hz) = 1002.4 s
  - f_ISCO(2.8 M☉) = 1570.4 Hz
- **Horizon of the S&S analytic aLIGO fit** (1.4+1.4, SNR 8, 20 Hz to f_ISCO, Euclidean): about 496.5 Mpc, giving a range of about 219 Mpc via /2.264. This is about 10% above the 449 Mpc that Chen et al. report for the O4high curve. The fit is approximate, so use tabulated curves through PyCBC/LAL for tight tests, and keep the fit only as a smoke test.
- **LISA, Robson Eq. 13 (instrument only) and Eq. 14 (4-yr confusion added):**

  | f | √S_n instrument (Hz⁻¹ᐟ²) | with 4-yr confusion (Hz⁻¹ᐟ²) |
  |---|---|---|
  | 0.1 mHz | 4.60×10⁻¹⁷ | 4.63×10⁻¹⁷ |
  | 1 mHz | 1.28×10⁻¹⁹ | 2.88×10⁻¹⁹ |
  | 3 mHz | 1.74×10⁻²⁰ | 1.77×10⁻²⁰ |
  | 10 mHz | 1.20×10⁻²⁰ | 1.20×10⁻²⁰ |
  | 100 mHz | 4.61×10⁻²⁰ | 4.61×10⁻²⁰ |

  The instrumental minimum is about 1.18×10⁻²⁰ at about 7.3 mHz. Compare with legwork or lisatools, which implement Robson's curve.
- **ET fit vs aLIGO fit (√S):** at 10, 100 and 1000 Hz the ET fit gives 5.7×10⁻²⁴, 3.2×10⁻²⁵ and 6.2×10⁻²⁵ Hz⁻¹ᐟ², against 1.8×10⁻²², 2.9×10⁻²⁴ and 1.4×10⁻²³ for the aLIGO fit. The ET fit's Euclidean BNS horizon (no sin 60°, no redshift) of about 6 Gpc is unphysical, and shows that 3G horizons **require** the cosmological treatment of Chen et al.

### Gaps
- I found no published horizon or range for the aLIGO "design" curve (T1800044) in a fetched source. Search snippets mentioned "428 Mpc" and "430 Mpc" design BNS horizons from unidentified papers, and they are not used here. Compute this with PyCBC against the T1800044 file, or use the Chen et al. O4high numbers above.
- I did not re-fetch the GW150914 "M_c ≈ 30 M☉" figure from the f–ḟ estimate in the PRL body.

## 5. What the community lacks that Rust would suit, and a scope we can validate

### Takeaway
Waveform generation in Rust is now taken: waveformrs covers TaylorF2 through IMRPhenomXPHM and NRSur7dq4, validated against LAL, but under GPL-3. The best open niche is a **permissively licensed, dependency-light "sensitivity and SNR" crate**. Its parts:
- analytic and tabulated PSDs (aLIGO, ET, CE, LISA with confusion)
- inner product, SNR, horizon, range and reach
- antenna patterns for L-shaped and triangular detectors
- SGWB and PTA tools (Ω ↔ S_h ↔ h_c, ORFs, Hellings–Downs, PI curves, NANOGrav power-law conversions)
- optionally a clean-room TaylorF2 3.5PN for self-contained use

It should compile to `no_std`/WASM for browser plotters and offer an extensible high-frequency (MHz–GHz) curve registry. Nothing comparable exists in Rust, and the Python equivalents (gwent, hasasia, pygwb, GWFish) are the least maintained part of the ecosystem.

### Cited findings
- High-frequency GWs:
  - The 2025 Living Review "Challenges and Opportunities of Gravitational Wave Searches above 10 kHz" (Aggarwal et al., LRR 28, 10 (2025); v2 Feb 2026) compares detector concepts with signal models. No code repository is listed on its arXiv page — [arXiv:2501.11723](https://arxiv.org/abs/2501.11723)
  - The only tooling found is a **web plotter** by Muia, Ringwald and Tamarit, "High frequency gravitational wave plotter: Noise-equivalent strain" (Zenodo 10.5281/zenodo.15720342). It is based on that review's data and accepts user CSV curves — [Shplot (Mainz)](https://incandenza-01.zdv.uni-mainz.de/Shplot/) (via search summary)
- Existing sensitivity-curve tools are Python or web, band-specific, and partly stale:
  - gwent: last release 2023
  - hasasia: last release 2021, PTA only
  - LEGWORK: active, but LISA and stellar sources only
  - GWplotter: a web app from 2015
  - GWSC.jl: Julia, PTA-oriented
  - Sources: [gwent](https://pypi.org/project/gwent/), [hasasia](https://pypi.org/project/hasasia/), [LEGWORK](https://github.com/TeamLEGWORK/LEGWORK), [GWSC.jl](https://github.com/bingining/GWSC.jl), [arXiv:1408.0740](https://arxiv.org/abs/1408.0740)
- Validation culture in the new Rust crates is high:
  - waveformrs checks every model against LALSimulation and keeps a `lal-quirks.md` file.
  - distmarg cross-validates independent Rust and Python implementations against mpmath fixtures.
  - Sources: [waveformrs](https://github.com/ColmTalbot/waveformrs), [distmarg](https://github.com/ColmTalbot/distance_marginalization)
- WASM:
  - libgrav advertises a `grav-wasm` npm binding with a vectorised Float64Array API, but only for parameter conversions, and the npm package was not found in the registry — [libgrav](https://github.com/transientlunatic/libgrav)
  - No Rust or WASM SNR or sensitivity-curve tool was found on crates.io or GitHub — [crates.io search](https://crates.io/search?q=gravitational%20wave)
- Existing Rust building blocks to reuse rather than rewrite:
  - bilby_rust: antenna response and GMST, MIT
  - igwn-ligolw: LIGO_LW XML, MIT
  - psrutils: par/tim parsing, MIT/Apache
  - ligo_hires_gps_time: GPS time, GPL-3
  - cdshealpix: sky maps
  - Sources: [bilby_rust](https://github.com/ColmTalbot/bilby_rust), [igwn-ligolw](https://crates.io/crates/igwn-ligolw), [psrutils](https://crates.io/crates/psrutils), [ligo_hires_gps_time](https://crates.io/crates/ligo_hires_gps_time)

### Inferences
- **Recommended crate 1, `gwsens` (working name), MIT/Apache-2.0, `no_std` + `alloc`, optional `wasm-bindgen`:**
  1. **psd**:
     - S&S aLIGO and ET fits
     - Robson LISA with the confusion table
     - a log-log interpolator for tabulated ASDs (CE T2000017, ET-D/2023, aLIGO O4high)
     - conversions S_n ↔ h_n ↔ Ω_eff
  2. **snr**:
     - the inner product
     - optimal SNR for a supplied h̃
     - a TaylorF2-Newtonian fast path ρ ∝ √∫f^(−7/3)/S_n
     - horizon (Euclidean and FLRW, which needs a small z ↔ D_L integrator), range (2.264) and reach
  3. **antenna**: F+ and F× for the L-shape and the 60° triangle (three rotated detectors), plus polarisation-angle averaging.
  4. **sgwb**:
     - power-law conversions (Thrane–Romano Eqs. 3–9)
     - S_eff for networks (Eq. 23)
     - PI curves (Eqs. 29–30) for any detector, including PTAs (Hellings–Downs Eq. 19–20) and user-supplied HFGW curves
  5. **pta**: γ ↔ α, A ↔ Ω_gw h², and NANOGrav 15-yr presets.

  **Validation targets:**
  - LAL/PyCBC `sigma` and analytic PSDs
  - Chen et al. Table 1 (using the O4high file)
  - Thrane–Romano code (LIGO-P1300115)
  - legwork or lisatools for LISA
  - hasasia or enterprise for Hellings–Downs
  - the GW150914 and NANOGrav numbers in §4

  Every formula is closed-form or a 1-D quadrature, so tests can be tight: relative 1e-10 against the formula, and about 1e-6 against reference code with matched integration grids.
- **Optional crate 2, `taylorf2` (permissive)**: non-spinning 3.5PN phase and Newtonian amplitude, written from S&S Eq. 127–129 rather than translated from LAL code, so it can stay MIT. Validate against LAL TaylorF2 with spin orders off and amplitude order 0, and against waveformrs's TaylorF2. The sign-convention notes in §3 apply. Writing from the paper avoids GPL and keeps the browser/WASM path permissive. If GPL is acceptable, depend on waveformrs instead.
- **Optional crate 3, GWOSC HDF5 reader**: strain plus the DQ and INJ bitmask decoder, through the `hdf5` crate. It has a narrow, well-documented spec (§2), and its test fixtures can be the per-event 32-s files. **A GWF reader is a real gap but a larger project**, best deferred.
- **Rust's advantages here:**
  - **WASM** for interactive plotters. GWplotter-style pages and the HFGW plotter are currently JavaScript or Python web apps.
  - **Deterministic, bit-reproducible curves** for papers, with no SciPy-version drift in quadrature.
  - **Embedded or edge use**, where a TaylorF2-Newtonian template bank or a horizon calculator fits in `no_std`.
  - **A single source of truth** that can be bound to Python with PyO3 and to JavaScript.

  Real-time detection in production remains dominated by GstLAL, PyCBC Live, MBTA and SPIIR, so "edge detection" is a demonstration rather than a gap the community is asking to have filled.
- **HFGW.** No standard tool covers MHz–GHz. A Rust registry of noise-equivalent strain curves would fill that gap if it shipped data digitised from LRR 28:10 or the Zenodo record with provenance, used a common (f, S_n or h_c) schema, and converted directly to Ω_gw and PI curves. It is the most novel contribution in this list, but it also carries the highest risk: validation is limited to reproducing the review's figures, because no reference code exists.

### Gaps
- I did not determine the HFGW plotter's data licence or exact Zenodo contents.
- I did not check whether waveformrs or bilby_rust build for `wasm32`, or whether their authors plan permissive relicensing or `no_std`.
- I did not survey existing browser GW visualisations (for example GWOSC web apps or educational "chirp" apps) to measure WASM demand. That demand is inferred, not sourced.
