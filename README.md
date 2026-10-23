# quantum-gravimeter

Rust crates for satellite gravimetry. They compute geoid, gravity and water-storage fields from open
GRACE / GRACE-FO / GOCE / Swarm data, and the same error spectra and maps for simulated
quantum (cold-atom) sensors.

| Crate | Status | Scope |
|---|---|---|
| [`geopotential`](crates/geopotential) | 0.1.0, unpublished | Legendre functions to degree 2700, coefficients with tide systems, GM/R rescaling and ICGEM time-variable terms, point and grid synthesis, functionals, degree spectra |
| [`gravity-field`](crates/gravity-field) | 0.1.0, unpublished | ICGEM `.gfc` reader (static and time-variable icgem1.0/2.0, with errors), tested on real COST-G, Swarm and GOCO06s files; TN-14 (C20/C30) and TN-13 (degree 1) readers. Planned: GRACE Level-2 and mascon NetCDF readers, data fetchers |
| [`gravimetry`](crates/gravimetry) | 0.1.0, unpublished | Load Love numbers, Gaussian smoothing, Swenson–Wahr destriping, baseline means and anomalies. Planned: GIA, DDK, basins |
| [`satgrav`](crates/satgrav) | 0.1.0, unpublished | Noise models of electrostatic and cold-atom accelerometers and laser ranging (Zingerle et al. 2024), pair line-of-sight noise, hybrid stand-in; atom-interferometer physics (phase, shot noise, accelerometer and gradiometer noise, checked against CARIOQA, GRICE and Trimeche et al.); closed-form gradiometry error propagation; noise colouring and whitening (FFT); a closed-loop GRACE-type pair simulation (line-of-sight gravity difference, whitened least squares). Planned: time-variable gravity and aliasing, orbit dynamics |
| [`gravviz`](apps/gravviz) | app, not published | The figures below, and the data for the viewer |
| [`viewer`](viewer) | web page ([published](https://claude.ai/artifact/J7pKCUXkF2R25Ck9ZyBVB5), private) | Month-by-month water-storage maps (flat or globe) and the 2002–2026 series at any place |

## Figures

[`figures/noise_spectra.svg`](figures/noise_spectra.svg) reproduces Fig. 1 of Zingerle et al.
(2024) from the paper's equations. [`figures/water_2026-07.svg`](figures/water_2026-07.svg) is the
July 2026 GRACE-FO water-storage anomaly. [`figures/gradiometry_errors.svg`](figures/gradiometry_errors.svg)
compares GOCE-like and cold-atom gradiometer missions against the GOCO06s signal.
[`figures/error_triangle_2026-07.svg`](figures/error_triangle_2026-07.svg) shows the July 2026
field's calibrated error per coefficient.
[`figures/sst_closed_loop.svg`](figures/sst_closed_loop.svg) is the closed-loop simulation of four
GRACE-type mission classes against a real month's errors (Zingerle et al. 2024, Fig. 2a in kind). The line figures and the map have a CSV of the plotted
values beside them.

```bash
cargo run --release -p gravviz -- noise figures
cargo run --release -p gravviz -- water-map MONTH.gfc BASELINE_DIR TN-14.txt TN-13.txt figures
cargo run --release -p gravviz -- gradiometry GOCO06s.gfc figures
cargo run --release -p gravviz -- triangle MONTH.gfc figures
cargo run --release -p gravviz -- series MONTHS_DIR BASELINE_DIR TN-14.txt TN-13.txt viewer/data
cargo run --release -p gravviz -- sst figures REFERENCE_MONTH.gfc   # ~5 min, ~1 GB per class
```

`series` writes the viewer's data (`viewer/data/`, not committed) from every month in MONTHS_DIR.
ICGEM serves the whole COST-G series as one zip (the "zip" link on the series page), which avoids
its per-file rate limit.

The water map needs ICGEM COST-G monthly files (`https://icgem.gfz.de/series/02_COST-G_/Grace-Grace-FO_RL02.1`) and the TN-14/TN-13 files from the
[PO.DAAC docs bucket](https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-docs/gracefo/open/docs/).

The research behind this, covering open data access, quantum missions, formulas and the plan, is in
[`reports/`](reports/Satellite%20gravimetry%20data%20and%20Rust%20tooling.md), with source notes in
[`research_notes/`](research_notes/).

```bash
cargo test
```

Licensed under either of [MIT](LICENSE-MIT) or [Apache-2.0](LICENSE-APACHE) at your option.
