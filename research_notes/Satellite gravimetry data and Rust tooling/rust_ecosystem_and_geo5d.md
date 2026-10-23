# Rust ecosystem for satellite-gravimetry computation and visualization, and an evaluation of geo5d

All crates.io figures (version, last update, total and 90-day downloads, license, MSRV) come from the crates.io API, queried on 2026-10-22. GitHub figures come from the GitHub REST API on the same day. "Last" means the crate's `updated_at`, which is normally the date of the last publish.

## 1. What is geo5d, who owns it, and should it be extended or flanked by companion crates?

### Takeaway
geo5d (v0.2.0, 2026-09-14, Apache-2.0) is a small, two-day-old work-in-progress physical-geodesy crate by Daniel Garden (crates.io `dan-garden`, GitHub org FusionVault). It already contains a working point-wise spherical-harmonic gravity engine: Holmes–Featherstone recursion with first and second derivatives, geopotential, geoid/height anomaly, anomalies, deflections, the gradient tensor, and an ICGEM `.gfc` parser. It has none of the satellite-gravimetry layer: no grid synthesis or analysis, no time-variable or GRACE Level-2 handling, no filtering, no EWH or mascons, no NetCDF. The best structure keeps geo5d as the "physics and relativity" façade and adds companion crates (an SH core, data I/O, GRACE processing, viz), with geo5d's synthesis moved into, or delegating to, the new SH core.

### Cited Findings
**Identity and metadata**
- Description: "The geopotential & relativistic dimension of geodesy in pure f64: normal gravity, height systems, chronometric leveling (clock frequency ↔ geopotential ↔ height), and relativistic clock rates. Built on geo4d/geo3d. (Work in progress.)" Keywords: clocks, geodesy, geopotential, gravity, relativity. Categories: mathematics, aerospace, science::geo. — [crates.io API geo5d](https://crates.io/api/v1/crates/geo5d)
- Versions: 0.1.0 (2026-09-14T12:52Z) and 0.2.0 (2026-09-14T13:25Z), both published the same day. 29 total downloads, 0 reverse dependencies, edition 2021, `rust-version = 1.85`, license Apache-2.0 (single license, not dual MIT/Apache). ~850 Rust code lines in 7 files per the crates.io linecount. Crate size is 29.6 kB. — [crates.io API geo5d](https://crates.io/api/v1/crates/geo5d)
- Dependencies (0.2.0): `geo4d ^0.1.1` (normal) and `serde ^1` (optional). The only feature is `serde = ["dep:serde", "geo4d/serde"]`. — [crates.io API dependencies](https://crates.io/api/v1/crates/geo5d/0.2.0/dependencies); [Cargo.toml](https://github.com/FusionVault/geo5d/blob/main/Cargo.toml)
- docs.rs built 0.2.0 successfully. Public items: re-exported geo4d/geo3d types (`Duration, Ecef, Ellipsoid, Epoch, Geodetic, StateVector, Trajectory, Vec3, GRS80, WGS84`), `GeopotentialNumber` and `NormalGravity`, plus these modules:
  - `gravity::{SphericalHarmonicModel, GravityModel, GfcError}`
  - `relativity::{height_difference, geopotential_difference, fractional_frequency_shift, fractional_rate_tt, fractional_rate_tcg(_of), height_uncertainty, sync_fractional_offset, accumulated_offset_seconds}`
  - `transfer::{sagnac_correction, shapiro_delay, earth_shapiro_delay, gnss_eccentricity_correction, one_way_transfer, two_way_sagnac_residual, proper_time_offset}`
  - `constants` (GRS80/WGS84 parameters, `W0_IAU`, `W0_IHRS`, `L_G`, `L_B`, `GM_EARTH`, `GM_SUN`, …)

  — [docs.rs geo5d 0.2.0 all items](https://docs.rs/geo5d/0.2.0/geo5d/all.html)
- Repository: [github.com/FusionVault/geo5d](https://github.com/FusionVault/geo5d). Public, Apache-2.0, created 2026-09-14, last push 2026-09-14, 0 stars, 2 commits. Files: `src/{constants,gravity,height,lib,normal,relativity,transfer}.rs`, `examples/chronometric.rs`, a CI workflow, CHANGELOG, NOTICE. `gravity.rs` is 27 kB, the largest file. — [GitHub API repo](https://api.github.com/repos/FusionVault/geo5d); [GitHub tree](https://github.com/FusionVault/geo5d)

**The gravity engine as implemented (source read from `src/gravity.rs`)**
- The `GravityModel` trait has 4 methods: `gm()`, `reference_radius()`, `potential(Ecef)` and `acceleration(Ecef) -> Vec3`. Both return the gravitational part only, without the centrifugal term. — [gravity.rs](https://github.com/FusionVault/geo5d/blob/main/src/gravity.rs)
- `SphericalHarmonicModel { gm, a, omega, nmax, cnm: Vec<f64>, snm: Vec<f64> }` holds fully normalized (4π, no Condon–Shortley phase) coefficients in triangular index `n(n+1)/2+m`. Methods:
  - `new` (panics via `assert!` if arrays are too short), `grs80_normal_field()`, `wgs84_normal_field()`, `from_gfc(&str)`, `with_omega`, `max_degree`
  - `geopotential` (W = V + centrifugal), `gravity` (vector), `gravity_magnitude`, `disturbing_potential`
  - `geoid_undulation(position, normal, gravity, w0)`, `gravity_disturbance`, `gravity_anomaly` (spherical approximation −∂T/∂r − 2T/r), `deflections` (ξ, η), `gradient_tensor_enu` (Marussi tensor in ENU)

  — [gravity.rs](https://github.com/FusionVault/geo5d/blob/main/src/gravity.rs)
- Synthesis is point-by-point. A double loop over m and n with a Holmes–Featherstone column recursion uses a global scale factor `HF_SCALE = 1e280`. The sectorial is computed as `√3·u·1e280·∏ u·√((2k+1)/2k)`, and the loop accumulates 12 lumped sums per order, so cost is O(nmax²) per point. There is no grid or FFT path and no parallelism. — [gravity.rs](https://github.com/FusionVault/geo5d/blob/main/src/gravity.rs)
- `from_gfc` reads `earth_gravity_constant`, `radius` and `max_degree`, then lines beginning `gfc` or `gfct` (`n m C S`). It ignores every other key and line type. The rotation rate defaults to GRS80's. — [gravity.rs](https://github.com/FusionVault/geo5d/blob/main/src/gravity.rs)
- `geoid_undulation` returns the height anomaly / quasigeoid ζ = T/γ − (W₀−U₀)/γ. The docs say the classical geoid differs by the Bouguer term, up to ~1 m in the Himalaya, which is out of scope. Geodetic→ECEF is done with `geo4d::WGS84` whichever normal field is passed. — [gravity.rs](https://github.com/FusionVault/geo5d/blob/main/src/gravity.rs)

**Validation claims (from the author's docs)**
- The 0.2.0 changelog says |∇V| matches pyshtools to 0.05 µGal (5×10⁻¹¹) across a global grid on EGM96. After adding the zero-degree term, the geoid "agrees with GeographicLib on the EGM96 geoid to a uniform ~9 cm (a W₀/tide-system convention constant) away from high terrain". 0.2.0 had two breaking changes: `geoid_undulation` now takes `w0`, and `NormalGravity` gained a `u0` field. — [CHANGELOG.md](https://github.com/FusionVault/geo5d/blob/main/CHANGELOG.md)
- Internal checks: GRS80 even-zonal synthesis reproduces Somigliana normal gravity, the gradient tensor satisfies Laplace's trace identity, and an analytic (n−1) anomaly test. Sign conventions are unit-tested against the Tokyo Skytree clock-redshift experiment (Takamoto et al. 2020). — [README via crates.io](https://crates.io/crates/geo5d)
- Stated roadmap: a classical-geoid option via the Bouguer term, higher-degree synthesis benchmarks, and the ns-level periodic TDB−TT series. Out of scope: X-ray pulsar navigation and Conformal Geometric Algebra. — [README via crates.io](https://crates.io/crates/geo5d)
- The README is internally inconsistent. It says "Not yet published; the remaining steps are an external cross-check against GeographicLib on a real EGM model", yet the crate is published and the same README and CHANGELOG report that cross-check as done. — [README via crates.io](https://crates.io/crates/geo5d); [CHANGELOG.md](https://github.com/FusionVault/geo5d/blob/main/CHANGELOG.md)

**The family it builds on**
- geo4d 0.1.1 (2026-09-14, Apache-2.0, MSRV 1.85, 39 downloads) covers 4D geodesy:
  - TAI-based ns `Epoch` with TAI/GPS/UTC/TT scales
  - `Geodetic4`/`Ecef4`/`StateVector`/`Trajectory` with linear, Hermite and Catmull-Rom interpolation (CCSDS OEM-style)
  - `Helmert14` (ITRF and GDA2020), the ITRF2020 `PlateMotion` model
  - closest approach and conflict, collision probability (Chan)

  It depends on `geo3d ^0.2.2`. — [crates.io API geo4d](https://crates.io/api/v1/crates/geo4d); [geo4d README](https://crates.io/crates/geo4d)
- geo3d 0.3.4 (2026-10-20, Apache-2.0, MSRV 1.85, 532 downloads, ~7.4k code lines) is the base layer. It has had 10 releases since 2026-09-13. Features: `std` (default), `libm` (no_std), `footprint`, and interop with `glam ^0.33`, `nalgebra ^0.35`, `mint`, `geo-types`. Categories include `no-std::no-alloc` and `wasm`. Its dev-dependencies (`geographiclib-rs =0.2.7`, `map_3d =0.1.7`, `sgp4 =2.4.0`, `geo ^0.31`) are used for cross-checks. — [crates.io API geo3d](https://crates.io/api/v1/crates/geo3d)

**Owner and their other crates**
- Sole owner: user `dan-garden` ("Daniel Garden", GitHub dan-garden). The crates.io account was created 2026-09-08. — [crates.io owners API](https://crates.io/api/v1/crates/geo5d/owners)
- The same user publishes 28 crates. The FusionVault GitHub org has 4 public repos: datafusion-opensearch, geo3d, geo4d, geo5d. — [crates.io user crates API](https://crates.io/api/v1/crates?user_id=440573&per_page=100); [GitHub org repos API](https://api.github.com/orgs/FusionVault/repos)
  - Geodesy: geo3d, geo4d, geo5d
  - FusionVault data platform: fv-compute, fv-compute-container, fv-compute-wasm, fv-connect, fv-contract, fv-mac, fv-plan, fv-streams (+ -connectors/-engine/-exchange/-formats/-ops/-runtime/-sql/-state/-types at 0.8.0), fv-value, fv-value-datafusion, datafusion-opensearch
  - Utilities: downsampling (LTTB/Visvalingam)
  - pageindex, pageindex-cli, pageindex-flash, pageindex-pdf

### Inferences
- geo5d is very probably the requester's own crate: owner "Daniel Garden", local user `dan`. If so, refactoring, renaming and re-licensing are within the developer's control. Nobody depends on it (0 reverse deps, 29 downloads), so breaking changes now cost almost nothing.
- **Fit for satellite gravimetry.** geo5d's point synthesis is fine for single-point or orbit evaluation and for GRACE-scale degrees (nmax ≈ 60–120). It is missing everything that turns static geodesy into satellite gravimetry:
  - Time-variable fields are silently mishandled. `gfct` lines are read as static, and the ICGEM `trnd`/`acos`/`asin` terms are dropped. The `norm`, `tide_system` and `errors` header keys, and the σ columns, are ignored.
  - There is no GRACE Level-2 (GRCOF2/YAML-header) reader, no Love numbers, and no degree-1/C20/C30 replacement or GIA.
  - There are no filters (Gaussian, DDK, destriping), no EWH conversion, and no mascon support.
  - There is no gridded synthesis (per-latitude Legendre reuse plus a longitude FFT), no analysis (grid → coefficients), and no degree or error-degree variances.
  - There is no NetCDF/GeoTIFF output and no rayon.
- **Correctness risks to test before reuse.**
  - (a) The `1e280` scaling with u^m multiplied directly into the sectorial may underflow for high m near the poles. It is only validated on EGM96 (n=360), not at EGM2008's n=2190. Holmes & Featherstone's own scheme factors u^m out and applies it via Horner.
  - (b) The "uniform ~9 cm" geoid offset against GeographicLib should be resolved: W₀, zero-degree term, or tide system. GeographicLib reproduces NGA to 78 pm on EGM2008 (§4), so 9 cm is a convention difference, not agreement.
  - (c) `assert!` panics in `new()` and `from_gfc` silently drop malformed `nmax` data. Use `Result` at the trust boundary.
- **Dependency drift.** geo5d → geo4d 0.1.1 → geo3d ^0.2.2, while geo3d is now at 0.3.4. A downstream crate that uses geo3d 0.3 types directly alongside geo5d gets two semver-incompatible copies of geo3d with non-interchangeable `Ecef`/`Geodetic` types. Bump geo4d and geo5d before building on them.
- **Recommendation: companion crates, not one big geo5d.**
  - Keep geo5d as the façade for normal gravity, heights, relativity and time transfer, plus a `GravityModel` trait boundary.
  - Put the SH math (normalizations, ALF recursions, coefficient containers with time variability, grid synthesis and analysis, spectra) in a new no_std-friendly core crate that geo5d depends on. This replaces geo5d's private `synthesize`, so there is one engine to validate.
  - Put formats (ICGEM gfc v2 incl. time-variable terms, GRACE L2, Love numbers, TN-13/TN-14 tables, mascon NetCDF) in an I/O crate with optional C-backed features.
  - Put GRACE processing (filters, EWH, mascons, time series) in its own crate.
  - Keep visualization in a separate app crate. Its MSRV is much higher (§5) and it should not burden the library MSRV of 1.85.
- The family's geo3d→geo4d→geo5d naming suggests continuing the pattern, but satellite gravimetry is not a "6th dimension". Descriptive names (§3) will be easier to find on crates.io.

### Gaps
- I did not download the `.crate` tarballs or the `tests/` directory, so the full test suite and exact validation tolerances were not inspected. The repo tree lists no `tests/` directory; tests are probably inline in the modules.
- I could not confirm whether the README's validation scripts (pyshtools and GeographicLib comparisons) are committed anywhere. None appear in the repo tree.
- I did not verify the underflow hypothesis numerically, since no code was run per the constraints.

## 2. Existing Rust crates to depend on (SH/Legendre, geodesy, NetCDF/HDF5, numerics, orbits, gravity)

### Takeaway
No mature Rust crate covers geodetic spherical-harmonic gravity fields. The closest are satkit (orbit-oriented, degree ≤70), astrodyn_gravity (JEOD port), geo5d itself, and the brand-new apollo-sht (generic complex SHT). Rust-native "EGM" crates are only interpolation grids or empty placeholders. The general substrate is mature: ndarray/faer/nalgebra, rustfft/realfft, rayon, proj4rs, hifitime/anise/sgp4, the georust netcdf/hdf5 bindings, and the new pure-Rust NetCDF-4/HDF5 readers. The geodesy-normalized SH core has to be written, and geo5d is the starting point.

### Cited Findings
**Spherical harmonics / Legendre / gravity in Rust**

| Crate | Version (last) | Downloads (90 d) | License / MSRV | Notes |
|---|---|---|---|---|
| sphrs | 0.2.2 (2023-05-14) | 18,825 (541) | MIT OR Apache-2.0 | Graphics/physics SH and solid harmonics, inspired by Google's spherical-harmonics and Robin Green's "Gritty Details". Not geodesy 4π normalization, not built for high degree. [crates.io](https://crates.io/crates/sphrs) |
| spherical-harmonics | 0.1.0 (2019-06-02) | 1,850 (1) | Apache-2.0 | Abandoned. [crates.io](https://crates.io/crates/spherical-harmonics) |
| apollo-sht | 0.9.1 (2026-10-19, created 2026-08-03) | 157 | MIT, MSRV 1.95 | Complex orthonormal SH with Condon–Shortley phase, Gauss–Legendre quadrature SHT plans, optional wgpu backend through its own "Hephaestus/Leto/Mnemosyne" ecosystem. Very new. [crates.io](https://crates.io/crates/apollo-sht) |
| astrodyn_gravity | 0.2.0 (2026-06-09) | 424 | MIT OR Apache-2.0, MSRV 1.89 | Port of NASA JEOD v5.4.0 gravity: Gottlieb normalized recursion, used to 360×360 (GGM05C), tides, third body, PN corrections. Bevy-oriented "astrodyn" workspace. [crates.io](https://crates.io/crates/astrodyn_gravity) |
| satkit | 0.24.1 (2026-10-01) | 170,294 (56,021) | MIT OR Apache-2.0 | Orbit toolkit. EGM96, EGM2008, JGM2, JGM3 and ITU GRACE16 compiled in, SH to degree/order 70, solid tides (IERS 2010), "ICGEM reference values" in tests, regression-tested against GMAT. [crates.io](https://crates.io/crates/satkit); [README](https://github.com/ssmichael1/satkit) |
| egm | 0.1.0 (2026-08-21) | 20 | MIT | Placeholder: `src/lib.rs` is the `cargo new` template `add()` function. [crates.io](https://crates.io/crates/egm); [lib.rs](https://github.com/jbangelo/egm-rs/blob/HEAD/src/lib.rs) |
| egm2008 | 0.1.0 (2023-11-01) | 55,310 (13,253) | BSD-3-Clause | FlightAware. Embedded geoid grid generated with NGA's Fortran, plus interpolation. Not SH. [crates.io](https://crates.io/crates/egm2008) |
| egm96 | 0.3.2 (2026-10-17) | 15,005 (4,618) | Zlib | Geoid height above WGS84 from EGM96. [crates.io](https://crates.io/crates/egm96) |
| kinavis-geoid | 0.1.0 (2026-10-17) | 10 | MIT OR Apache-2.0 | EGM96 grid embedded. [crates.io](https://crates.io/crates/kinavis-geoid) |
| japan-geoid | 0.6.0 (2025-12-08) | 21,758 | — | GSI Japan geoid grid. [crates.io search](https://crates.io/search?q=geoid) |

- Searching crates.io for "mascon" returns 0 crates and "icgem" returns 2 (geo5d, satkit). — [crates.io API search](https://crates.io/api/v1/crates?q=mascon)

**Geodesy and projections**

| Crate | Version (last) | Downloads (90 d) | License / MSRV | Notes |
|---|---|---|---|---|
| geo | 0.33.1 (2026-04-20) | 24.6 M (4.9 M) | MIT OR Apache-2.0, 1.88 | georust 2D primitives/algorithms. [crates.io](https://crates.io/crates/geo) |
| proj | 0.31.0 (2025-08-29) | 881,563 | MIT OR Apache-2.0, 1.85 | Bindings to C PROJ. [crates.io](https://crates.io/crates/proj) |
| proj4rs | 0.2.1 (2026-10-09) | 1.23 M (410,535) | MIT OR Apache-2.0 | Pure Rust, WASM JS API (proj4js replacement). Projections include `moll`, `robin`, `eqearth`, `ortho`, `laea`, `sinu`, `eqc`, `geos`, `tpers` and more. NTv2/GeoTIFF grid support. [crates.io](https://crates.io/crates/proj4rs); [PROJECTIONS.md](https://github.com/3liz/proj4rs/blob/main/PROJECTIONS.md); [README](https://github.com/3liz/proj4rs) |
| proj-core | 0.12.0 (2026-10-18, created 2026-03-21) | 7,228 | MIT OR Apache-2.0, 1.85 | Pure Rust with embedded EPSG. [crates.io](https://crates.io/crates/proj-core) |
| geodesy (Knudsen) | 0.15.0 (2026-03-03) | 192,308 (81,830) | MIT OR Apache-2.0 on crates.io (GitHub reports Apache-2.0), 1.87 | Operators include `gridshift`, `deformation`, `deflection` (coarse deflection of the vertical from a geoid grid) and geoid tide-system conversion. [crates.io](https://crates.io/crates/geodesy); [Rumination 002](https://github.com/busstoptaktik/geodesy/blob/main/ruminations/002-rumination.md) |
| geographiclib-rs | 0.2.7 (2026-02-17) | 21.7 M | MIT, 1.70 | Ports only the geodesic part (`Geodesic`, `GeodesicLine`, `PolygonArea`). No GravityModel. [docs.rs items](https://docs.rs/geographiclib-rs/0.2.7/geographiclib_rs/all.html) |
| nav-types | 0.5.2 (2023-07-10) | 728,198 | MIT | Stale. [crates.io](https://crates.io/crates/nav-types) |
| map_3d | 0.1.7 (2026-04-19) | 1.12 M | Apache-2.0 | [crates.io](https://crates.io/crates/map_3d) |
| wbprojection | 0.3.3 (2026-08-04) | 18,374 | MIT OR Apache-2.0 | Whitebox, PROJ-inspired, includes Mollweide. [crates.io](https://crates.io/crates/wbprojection) |
| d3_geo_rs | 3.2.4 (2026-09-04) | 69,575 | MIT | Port of d3-geo. [crates.io](https://crates.io/crates/d3_geo_rs) |

**NetCDF / HDF5**

| Crate | Version (last) | Downloads (90 d) | License / MSRV | Notes |
|---|---|---|---|---|
| netcdf | 0.12.1 (2026-07-26) | 356,197 (53,871) | MIT OR Apache-2.0, 1.77 | Bindings to C netCDF. `netcdf-src` builds from source. [crates.io](https://crates.io/crates/netcdf) |
| hdf5-metno | 0.15.0 (2026-09-18) | 2.74 M (712,021) | MIT OR Apache-2.0, 1.85.1 | Maintained fork of the C-binding `hdf5` crate. The original `hdf5` 0.8.1 has not been updated since 2021-11-21. [crates.io](https://crates.io/crates/hdf5-metno); [crates.io hdf5](https://crates.io/crates/hdf5) |
| netcdf-reader / netcdf-writer | 0.9.3 (2026-10-21, created 2026-03/07) | 28,116 / 3,400 | MIT OR Apache-2.0, 1.81 | Pure Rust, no C. Reads and writes CDF-1/2/5 and NetCDF-4 (HDF5-backed), chunked layouts, unlimited dims, deflate/shuffle/fletcher32 filters (lz4/nbit/scaleoffset read). Validated against libhdf5/netcdf-c. Companion `hdf5-reader`/`hdf5-writer`. [crates.io](https://crates.io/crates/netcdf-reader); [README](https://github.com/roteiro-gis/netcdf-rust) |
| hdf5-pure | 0.47.0 (2026-09-26, created 2026-03-14) | 134,183 | MIT OR Apache-2.0, 1.89 | Pure Rust read/write/edit-in-place, WASM and no_std+alloc, deflate/shuffle/LZF/scale-offset, SWMR. Optional ndarray. [crates.io](https://crates.io/crates/hdf5-pure) |
| rust-hdf5 | 0.7.8 (2026-10-20) | 12,500 | MIT, 1.89 | Pure Rust read/write/SWMR. [crates.io](https://crates.io/crates/rust-hdf5) |
| netcdf3 | 0.6.1 (2026-05-03) | 60,088 | MIT OR Apache-2.0 | Classic NetCDF-3 only. [crates.io](https://crates.io/crates/netcdf3) |

**Numerics**

| Crate | Version (last) | Downloads (90 d) | License / MSRV | Notes |
|---|---|---|---|---|
| ndarray | 0.17.2 (2026-01-10) | 135.8 M | MIT OR Apache-2.0, 1.64 | [crates.io](https://crates.io/crates/ndarray) |
| nalgebra | 0.35.0 (2026-05-24) | 97.0 M | Apache-2.0 only, 1.89 | [crates.io](https://crates.io/crates/nalgebra) |
| faer | 0.24.4 (2026-06-24) | 5.23 M (1.54 M) | MIT, 1.84 | Repository on Codeberg. [crates.io](https://crates.io/crates/faer) |
| rustfft | 6.4.1 (2025-09-18) | 33.1 M | MIT OR Apache-2.0, 1.61 | [crates.io](https://crates.io/crates/rustfft) |
| realfft | 3.5.0 (2025-06-12) | 18.9 M | MIT | [crates.io](https://crates.io/crates/realfft) |
| ndrustfft | 0.6.2 (2026-01-06) | 257,582 | MIT | [crates.io](https://crates.io/crates/ndrustfft) |
| rayon | 1.12.0 (2026-04-14) | 614 M | MIT OR Apache-2.0, 1.80 | [crates.io](https://crates.io/crates/rayon) |
| num-complex | 0.4.6 (2025-01-25) | 428.6 M | MIT OR Apache-2.0 | [crates.io](https://crates.io/crates/num-complex) |

**Orbits and time**

| Crate | Version (last) | Downloads (90 d) | License / MSRV | Notes |
|---|---|---|---|---|
| sgp4 | 2.4.0 (2026-02-23) | 410,433 | MIT | [crates.io](https://crates.io/crates/sgp4) |
| nyx-space | 2.6.0 (2026-09-12) | 79,483 (3,348) | AGPL-3.0-or-later | [crates.io](https://crates.io/crates/nyx-space); [GitHub](https://github.com/nyx-space/nyx) |
| anise | 0.10.6 (2026-10-19) | 227,379 | MPL-2.0 | [crates.io](https://crates.io/crates/anise) |
| hifitime | 4.3.1 (2026-08-07) | 1.34 M | MPL-2.0 | [crates.io](https://crates.io/crates/hifitime) |
| brahe | 1.8.0 (2026-10-16) | 4,054 | MIT, 1.95 | [crates.io](https://crates.io/crates/brahe) |
| lox-space | 0.1.0-alpha.51 (2026-10-14) | 19,842 | MPL-2.0, 1.90 | [crates.io](https://crates.io/crates/lox-space) |
| satkit | 0.24.1 | — | — | See the SH table above. |

### Inferences
- **SH core.** Write it yourself, starting from geo5d's Holmes–Featherstone code. Use realfft/rustfft for longitude FFTs in gridded synthesis and analysis, and rayon behind a feature for per-latitude parallelism. Leave sphrs and apollo-sht out: different normalization (orthonormal complex with Condon–Shortley phase), and apollo-sht is young with heavy ecosystem deps and MSRV 1.95.
- **Cross-checks.** satkit (MIT/Apache, compiled-in EGM2008/GRACE16 to degree 70) is a useful Rust-native check for low-degree acceleration on orbits.
- **Licences.** Avoid nyx-space (AGPL-3.0) as a dependency of an Apache/MIT library. anise and hifitime are MPL-2.0, which is file-level copyleft and acceptable as unmodified dependencies. geo4d already provides TAI/GPS/TT epochs, so hifitime is optional.
- **I/O.** Default to the pure-Rust `netcdf-reader` (and `hdf5-pure`) for WASM friendliness and zero C builds. Offer `netcdf`/`hdf5-metno` C bindings behind a feature for edge-case compression and performance. The pure-Rust readers are ~7 months old, so pin versions and add round-trip tests against real GRACE mascon files.
- **Linear algebra.** Use faer for heavier work (least squares and normal equations if gravity recovery is ever done), ndarray for gridded fields, and nalgebra only where geo3d interop requires it (geo3d already offers a `nalgebra` feature).

### Gaps
- I did not inspect every candidate's API docs (for example whether satkit exposes a reusable SH evaluator outside its propagator, or whether proj4rs `moll`/`robin` are spherical-only on ellipsoidal input; PROJECTIONS.md marks them "spherical").
- Maintenance quality (open-issue counts, bus factor) was only checked via last-release dates and downloads.

## 3. crates.io name availability (checked 2026-10-22)

### Takeaway
Most candidate names are free: gravity-field, geopotential, shgrav, grace-rs, gravimetry, geoid, icgem, mascon, sphharm, legendre-norm and others. Taken: `egm` (an empty placeholder), `gravity`, `grace`, `gfc`, `egm2008`, `egm96`.

### Cited Findings
- Method: `GET https://crates.io/api/v1/crates/{name}`, where HTTP 404 means available. The API treats `-` and `_` as the same name: `proj_core` resolves to `proj-core` and `hdf5_metno` to `hdf5-metno`. A capitalized lookup (`Geoid`, `ICGEM`) also returned 404, matching the lowercase result. — [crates.io API](https://crates.io/api/v1/crates/proj_core)
- **AVAILABLE (404):** gravity-field (and gravity_field), geopotential, shgrav, grace-rs, gravimetry, geoid, icgem, mascon, sphharm (and sph_harm), legendre-norm, shtools, grace-fo, geo6d, gravfield, gravity-model, geoid-rs, sht, shtns, gravity-toolkit, geogravity, satgrav, grace-l2, sh-synthesis. — [crates.io API](https://crates.io/api/v1/crates/geopotential)
- **TAKEN:**
  - `egm`: 0.1.0, 2026-08-21, jbangelo/egm-rs, 20 downloads, empty `add()` template — [crates.io](https://crates.io/crates/egm); [lib.rs](https://github.com/jbangelo/egm-rs/blob/HEAD/src/lib.rs)
  - `gravity`: 0.1.3, 2017-05-01, "Rust bindings for the Gravity API" — [crates.io](https://crates.io/crates/gravity)
  - `grace`: 0.1.0, 2021-01-01, Diggsey, "Intercept shutdown requests" — [crates.io](https://crates.io/crates/grace)
  - `gfc`: 1.0.0, 2023-03-03, a random-file CLI — [crates.io](https://crates.io/crates/gfc)
  - `egm2008`: FlightAware — [crates.io](https://crates.io/crates/egm2008)
  - `egm96`: micahcc — [crates.io](https://crates.io/crates/egm96)

### Inferences
- Candidate names by role:
  - SH core: `shgrav`, `sphharm` or `geopotential`. `geopotential` is the most descriptive and matches geo5d's keyword.
  - Data I/O: `gravity-field`, or `icgem` if the crate is only an ICGEM-format parser.
  - Processing: `gravimetry`, or `grace-l2` for a GRACE-specific crate.
  - Mascons: `mascon`.
- Avoid names of existing non-Rust projects or services (`shtools`, `shtns`, `gravity-toolkit`, arguably `icgem`, a GFZ service name). They could be read as official ports.
- crates.io names cannot be reserved without publishing. If names matter, publish minimal real 0.0.x/0.1.0 crates soon. Do not publish empty placeholders, which crates.io policy discourages, as the `egm` example shows.
- `geo6d` is free if the developer wants to keep the family numbering, but it says nothing about the content.

### Gaps
- Availability can change at any moment. Re-check just before publishing.
- I did not check crates.io's reserved-name list or its trademark and squatting policy text in this session.

## 4. Reference implementations to port or validate against, with licenses

### Takeaway
Port from permissive sources: SHTOOLS (BSD-3), gravity-toolkit (MIT) and GeographicLib (MIT), keeping attribution. Use GPL/CeCILL code (GROOPS GPL-3.0, ducc0 GPL-2.0, SHTns CeCILL-2.1, nyx AGPL) and the ICGEM calculation service only as external validation oracles. Reading their outputs is fine; translating their code into an Apache/MIT crate is not.

### Cited Findings
- **SHTOOLS/pyshtools.** v4.14.1, released 2026-02-22 (PyPI upload 2026-02-23), BSD-3-Clause, Fortran-95 and Python. It supports all standard normalizations and phase conventions and both regular (Driscoll–Healy) and Gauss–Legendre grids. Transforms are "proven to be accurate up to about degree 2800" natively, and beyond that with the DUCC0 backend. It also offers multitaper/Slepian localization, xarray/netcdf I/O, and gravity and magnetic field classes. 431 stars, last push 2026-05-18. — [SHTOOLS README](https://github.com/SHTOOLS/SHTOOLS); [GitHub API](https://api.github.com/repos/SHTOOLS/SHTOOLS); [PyPI](https://pypi.org/project/pyshtools/)
- **gravity-toolkit (Tyler Sutterley).** Now at `polargeodesy/gravity-toolkit` (the old `tsutterley/...` URL redirects). MIT, PyPI v1.2.8 (uploaded 2026-10-02), last push 2026-10-13. Scope: "Python tools for obtaining and working with Level-2 spherical harmonic coefficients" from GRACE and GRACE-FO. Package modules:
  - Readers: `read_GRACE_harmonics`, `read_gfc_harmonics`, `read_SLR_harmonics` (+ `SLR/`), `read_love_numbers`, `read_GIA_model`
  - Filtering: `destripe_harmonics`, `gauss_weights`, `gen_averaging_kernel`
  - Synthesis and analysis: `harmonic_summation`, `clenshaw_summation`, `gen_stokes`, `gen_harmonics`, `harmonic_gradients`, `fourier_legendre`, `associated_legendre`, `legendre`
  - Loads: `gen_spherical_cap`, `gen_disc_load`, `gen_point_load`, `ocean_stokes`
  - Other: `geocenter`, `mascons`, `sea_level_equation`, `degree_amplitude`, `grace_date`, `grace_find_months`, `time_series/`, `units`

  The sister project `polargeodesy/model-harmonics` is also MIT. — [GitHub API](https://api.github.com/repos/polargeodesy/gravity-toolkit); [module listing](https://github.com/polargeodesy/gravity-toolkit/tree/main/gravity_toolkit); [PyPI](https://pypi.org/project/gravity-toolkit/); [docs](https://gravity-toolkit.readthedocs.io/en/latest/)
- **GROOPS (TU Graz ITSG).** GPL-3.0, C++, "A software toolkit for gravity field recovery and GNSS processing". Latest release tag 2025-11-15, last push 2026-09-12, 248 stars. It has cookbooks for GRACE gravity-field recovery, POD-based gravity, kinematic LEO orbits and regional geoid, and is driven by config files plus a GUI. — [GitHub API](https://api.github.com/repos/groops-devs/groops); [GROOPS docs](https://groops-devs.github.io/groops/html/index.html)
- **GeographicLib GravityModel (C++, docs v2.7).** Models egm84 (n=180), egm96 (360), egm2008 (2190), wgs84 and grs80 (20). Computes `GeoidHeight`, `SphericalAnomaly`, `Disturbance` and `T`. `GravityCircle` is about 800× faster per point for egm2008 when evaluating along a parallel. Data comes as a text `.egm` file plus a little-endian binary `.egm.cof`. Accuracy against NGA's synthesis: 1.1 mm (egm84), 23 nm (egm96), 78 pm for egm2008 (23 pm after fixing an NGA issue). The docs list NGA code discrepancies, for example a missing 1/r term in the egm2008 disturbing potential. The geographiclib GitHub repo is MIT (last push 2026-08-24), and PyPI `geographiclib` 2.1 is MIT. — [GeographicLib gravity docs](https://geographiclib.sourceforge.io/C++/doc/gravity.html); [GitHub API](https://api.github.com/repos/geographiclib/geographiclib); [PyPI](https://pypi.org/project/geographiclib/)
- **SHTns (N. Schaeffer, CNRS).** CeCILL 2.1, described as GPL-compatible, written in C. PyPI `shtns` 3.7.5 was uploaded 2026-02-13, and the home page is bitbucket.org/nschaeff/shtns. Its authors ask for citation of the 2013 G³ paper (doi 10.1002/ggge.20071). — [SHTns README mirror](https://vsl.cis.udel.edu/trac/civl/browser/CIVL/examples/omp/shtns/README.md); [anaconda listing](https://anaconda.org/hcc/shtns); [PyPI](https://pypi.org/project/shtns/)
- **ducc0 (M. Reinecke, MPCDF).** GPL-2.0, PyPI 0.42.0 uploaded 2026-10-20. GitHub mirror `mreineck/ducc` last pushed 2026-10-21. It is also the high-degree backend of pyshtools. — [GitHub API](https://api.github.com/repos/mreineck/ducc); [PyPI](https://pypi.org/project/ducc0/); [SHTOOLS README](https://github.com/SHTOOLS/SHTOOLS)
- **ICGEM calculation service (GFZ).** An online service that computes gravity-field functionals from static and temporal models, on user-selected grids or points, with optional coefficient-range selection or Gaussian filtering. Functionals include height anomaly, geoid, gravity disturbance, gravity anomalies (classical and modern), simple Bouguer anomaly and surface gravity. Formulas are documented in GFZ STR09/02 (Barthelmes 2013). Newer tools: the G3 Browser for temporal series (2023) and model-difference computation (2024). — [GFZ ICGEM page](https://www.gfz.de/en/section/global-geomonitoring-and-gravity-field/data-products-services/icgem); [IAG Geodesist's Handbook 2024 ICGEM entry](https://media.geodesy.science/iag/geodesists-handbook/2024/iag-gh2024_icgem.pdf); [ESSD 2019 ICGEM paper](https://essd.copernicus.org/articles/11/647/2019/essd-11-647-2019.html)
- **Fukushima X-numbers and ultra-high-degree ALFs.**
  - A Polish-journal comparison of underflow-avoidance methods found the Fukushima method fastest at equal accuracy for EGM2008's ~2.4 M Legendre functions — [journals.pan.pl](https://journals.pan.pl//Content/107428/PDF/292.pdf)
  - Goli, Foroughi & Novák (2022) propose "successive ratios" and "Midway" alternatives benchmarked against X-numbers — [Stud. Geophys. Geod. doi 10.1007/s11200-022-0830-9](https://springerlink.bibliotecabuap.elogim.com/10.1007/s11200-022-0830-9)
  - Enomoto (2015) found a Fourier-based method's errors much smaller than X-numbers for truncations 39–10239 — [J-STAGE SOLA](https://www.jstage.jst.go.jp/article/sola/11/0/11_2015-033/_article)
- **NGA EGM2008 synthesis.** NGA's Fortran synthesis routine is distributed from the official EGM2008 site. One source calls it HARMONIC_SYNTH.f, and the coefficient file is `EGM2008_to2190_TideFree`. — [UTM EGM2008 evaluation note](https://people.utm.my/zainalabidin/2016/09/03/egm2008-evaluation-procedure); [Wikipedia EGM](https://en.wikipedia.org/wiki/Earth_Gravitational_Model)
- **NASA JEOD.** The Gottlieb algorithm is already ported to Rust in `astrodyn_gravity` (MIT/Apache). — [crates.io astrodyn_gravity](https://crates.io/crates/astrodyn_gravity)

### Inferences
- Validation ladder for the Rust SH core:
  1. Closed-form normal field (already in geo5d).
  2. GeographicLib egm96/egm2008 geoid and anomaly, at pm–nm agreement, so ~9 cm offsets must be explained.
  3. pyshtools grids, the Driscoll–Healy/GLQ round trip, degree variances, and `SHGravCoeffs` functionals.
  4. The ICGEM calculation service for grid functionals with a Gaussian filter and the zero-degree/tide conventions.
  5. gravity-toolkit for the GRACE L2 chain: C20/C30 replacement, degree-1, GIA, destriping, Gaussian smoothing, Love-number EWH conversion, mascon handling.
  6. GROOPS or ducc0/SHTns outputs only as black-box numerical references.
- The best port candidates are gravity-toolkit (MIT, Python, the closest functional match to "GRACE processing") and GeographicLib's GravityModel/GravityCircle (MIT, C++, for fast parallel-circle synthesis and its well-documented conventions). SHTOOLS (BSD-3) is good for grid conventions and Slepian/localization. Keep a NOTICE file with upstream copyright notices.
- For degrees above ~2700 (beyond EGM2008, e.g. topographic models to 5400 and up), implement Fukushima X-numbers or a comparable extended-range method from the papers. Algorithms from papers carry no code licence.
- NGA's code is a US-government product, likely public domain. GeographicLib documents its bugs, so use NGA's code only as a historical reference.

### Gaps
- I could not fetch Fukushima's 2012 J. Geodesy paper or locate official Fortran code and its licence terms.
- I could not confirm the exact current file name (`hsynth_WGS84.f` vs `HARMONIC_SYNTH.f`) or the licence of NGA's EGM2008 synthesis code from NGA directly.
- The GeographicLib gravity page itself shows no licence statement. MIT is inferred from the GitHub repo and PyPI metadata.
- The SHTns licence was confirmed from mirrors and third-party listings, not the upstream LICENSE file.
- ducc0's exact SPDX ("GPL-2.0" vs "GPL-2.0-or-later") was not checked against the source headers.
- "GrapeGrav" is not a recognizable project, and I found nothing under that name.

## 5. Visualization stack in Rust and recommended pragmatic stacks

### Takeaway
(a) **Static, publication-quality maps:**
1. Compute grids in Rust.
2. Project them with proj4rs (pure Rust; `moll`, `robin`, `eqearth`).
3. Colour with embedded Crameri scientific colour maps (MIT, secondhand) or colorous/colorgrad.
4. Rasterize, then draw axes, colorbars and log-scale degree-variance plots with plotters (SVG/PNG).
5. Write NetCDF (netcdf-writer) and GeoTIFF (geotiff-writer) so the same grids can go to GMT/PyGMT for journal figures.

(b) **Interactive 3D globe:** an eframe/egui app with egui_plot for the 2D panels and a custom wgpu (or three-d) textured-sphere globe. It runs natively and in WASM. For rapid developer visualization, log meshes to Rerun. If a browser globe with basemaps is required, use a Rust→WASM compute core feeding CesiumJS (Apache-2.0) or MapLibre GL JS (globe projection since v5.0) through wasm-bindgen.

### Cited Findings
**2D plotting and GUI**

| Crate | Version (last) | Downloads (90 d) | License / MSRV | Notes |
|---|---|---|---|---|
| plotters | 0.3.7 (2024-09-08) | 253.5 M (59.1 M) | MIT | No release in over two years. GitHub last push 2026-04-13. [crates.io](https://crates.io/crates/plotters); [GitHub API](https://api.github.com/repos/plotters-rs/plotters) |
| egui / eframe | 0.37.0 (2026-10-15) | 27.6 M / 21.5 M | MIT OR Apache-2.0, 1.95 | [crates.io](https://crates.io/crates/egui) |
| egui_plot | 0.38.0 (2026-10-19) | 8.85 M | MIT OR Apache-2.0, 1.95 | [crates.io](https://crates.io/crates/egui_plot) |
| walkers | 0.63.0 (2026-10-16) | 1.24 M | MIT, 1.95 | Slippy-map widget for egui. [crates.io](https://crates.io/crates/walkers) |
| rsplot | 0.6.1 (2026-10-17) | 409 | MIT OR Apache-2.0 | silx-style scientific plotting for egui on wgpu. [crates.io](https://crates.io/crates/rsplot) |
| plotly | 0.14.1 (2026-02-15) | 4.59 M | MIT | [crates.io](https://crates.io/crates/plotly) |
| charming | 0.6.0 (2025-06-17) | 1.58 M | MIT OR Apache-2.0 | ECharts. [crates.io](https://crates.io/crates/charming) |
| map2fig | 0.7.8 (2026-06-27) | 294 | MIT | Publication-quality HEALPix maps with Mollweide. [crates.io](https://crates.io/crates/map2fig) |
| galileo | 0.2.1 (2025-07-11) | 8,849 | MIT OR Apache-2.0 | Map rendering engine, slow release cadence. [crates.io](https://crates.io/crates/galileo) |

**3D and GPU**

| Crate | Version (last) | Downloads (90 d) | License / MSRV | Notes |
|---|---|---|---|---|
| wgpu | 30.0.1 (2026-08-22) | 43.0 M | MIT OR Apache-2.0, 1.87 | [crates.io](https://crates.io/crates/wgpu) |
| bevy | 0.20.0 (2026-10-08) | 8.22 M | MIT OR Apache-2.0, 1.97.1 | [crates.io](https://crates.io/crates/bevy) |
| bevy_egui | 0.43.1 (2026-10-10) | — | MIT | [crates.io](https://crates.io/crates/bevy_egui) |
| kiss3d | 0.47.0 (2026-10-04) | 594,677 | BSD-3-Clause | [crates.io](https://crates.io/crates/kiss3d) |
| three-d | 0.19.0 (2026-04-17) | 473,179 | MIT | "including web". [crates.io](https://crates.io/crates/three-d) |
| rerun | 0.40.0 (2026-10-21) | 2.22 M (553,344) | MIT OR Apache-2.0, 1.96 | GitHub Apache-2.0, 11.6k stars. [crates.io](https://crates.io/crates/rerun); [GitHub API](https://api.github.com/repos/rerun-io/rerun) |

**Colormaps**

| Crate | Version (last) | Downloads | License / MSRV | Notes |
|---|---|---|---|---|
| colorous | 1.0.16 (2025-03-03) | 1.52 M | Apache-2.0 | d3-scale-chromatic schemes. [crates.io](https://crates.io/crates/colorous) |
| colorgrad | 0.9.0 (2026-09-28) | 1.64 M | MIT OR Apache-2.0, 1.85 | [crates.io](https://crates.io/crates/colorgrad) |
| prismatica | 1.0.0 (2026-03-31) | — | GPL-3.0-only | "308 scientific colormaps as compile-time Rust constants". [crates.io](https://crates.io/crates/prismatica) |
| colorcet | 0.2.2 (2025-12-13) | — | CC-BY-4.0 AND PSF-2.0 AND MIT | [crates.io](https://crates.io/crates/colorcet) |

- No crate is named `scientific-colourmaps` (404). — [crates.io API](https://crates.io/api/v1/crates/scientific-colourmaps)
- Crameri's Scientific colour maps:
  - Licensed MIT according to the pycmap catalogue and the colorlab README. This is secondhand: the Zenodo licence field was not seen. — [pycmap catalog](https://pycmap.readthedocs.io/en/latest/catalog/diverging/crameri%3Abroc/); [colorlab README](https://github.com/beddalumia/colorlab/blob/master/crameri/README.md)
  - Concept DOI 10.5281/zenodo.1243862 — [Zenodo v7.0.0](https://zenodo.org/records/4491293)
  - Suggested citation: Crameri, Shephard & Heron (2020) *Nature Communications* — [colorlab README](https://github.com/beddalumia/colorlab/blob/master/crameri/README.md)

**Raster and NetCDF output**

| Crate | Version (last) | Downloads | License / MSRV | Notes |
|---|---|---|---|---|
| png | 0.18.1 (2026-02-14) | — | MIT OR Apache-2.0 | [crates.io](https://crates.io/crates/png) |
| image | 0.25.10 (2026-03-10) | — | MIT OR Apache-2.0 | [crates.io](https://crates.io/crates/image) |
| tiff | 0.11.3 (2026-02-10) | — | MIT | [crates.io](https://crates.io/crates/tiff) |
| geotiff-writer | 0.9.0 (2026-10-18) | — | MIT OR Apache-2.0 | Pure-Rust GeoTIFF and COG writer. [crates.io](https://crates.io/crates/geotiff-writer) |
| gdal | 0.19.0 (2025-12-23) | — | MIT | Bindings to C GDAL. [crates.io](https://crates.io/crates/gdal) |
| netcdf-writer | 0.9.3 | — | — | Pure Rust, see §2. |

**Web and globe**

| Item | Version (last) | Downloads | License | Notes |
|---|---|---|---|---|
| wasm-bindgen | 0.2.129 (2026-09-25) | — | MIT OR Apache-2.0 | [crates.io](https://crates.io/crates/wasm-bindgen) |
| CesiumJS | 1.146 (2026-10-01) | — | Apache-2.0 | [GitHub API](https://api.github.com/repos/CesiumGS/cesium) |
| deck.gl | v9.4.0 (2026-09-05) | — | MIT | [GitHub API](https://api.github.com/repos/visgl/deck.gl) |
| globe.gl | — | — | MIT | [GitHub API](https://api.github.com/repos/vasturiano/globe.gl) |
| MapLibre GL JS | v6.15.0 (2026-10-18) | — | GitHub licence detector: NOASSERTION | Globe projection shipped in GL JS 5.0.0 (January 2025). MapLibre Native does not have globe yet (August 2026 newsletter). [GitHub API](https://api.github.com/repos/maplibre/maplibre-gl-js); [MapLibre globe roadmap](https://maplibre.org/roadmap/maplibre-gl-js/globe-view); [MapLibre newsletter Aug 2026](https://maplibre.org/news/2026-09-02-maplibre-newsletter-august-2026/) |

- pyshtools makes its publication maps with Cartopy and PyGMT. This is the de-facto geoscience standard the Rust output should be able to feed. — [SHTOOLS README](https://github.com/SHTOOLS/SHTOOLS)

### Inferences
- **(a) Static maps and spectra.** Use plotters despite its stale release: it is still the most-downloaded Rust plotting crate and supports log axes plus SVG and bitmap backends.
  - Render the field itself as a projected raster yourself: inverse-project each output pixel with proj4rs (`moll`/`robin`/`eqearth`), sample the grid, and apply a colour LUT. Then composite in plotters.
  - Embed Crameri LUTs (vik, roma, batlow, …, each a 256×3 table) with attribution rather than depending on GPL `prismatica`.
  - Also write CF-compliant NetCDF and GeoTIFF so GMT/PyGMT/Cartopy can produce final journal figures. This avoids rebuilding coastlines, graticules and labelled projections in Rust.
- **(b) Interactive globe.** eframe/egui is the lowest-friction single codebase for desktop and web. egui_plot covers degree-variance and time-series panels, and walkers gives a 2D slippy map. Draw the globe as a UV/icosphere mesh with the field as a texture or vertex colours, via egui's wgpu paint callback or three-d. Bevy is heavier (MSRV 1.97.1, ECS) and only worth it if the app becomes a 3D scene with orbits and satellites.
- **Developer tooling.** Rerun gives an immediate 3D globe or point-cloud view of synthesized fields and orbit tracks with no UI code.
- **Browser product.** If basemaps, terrain and imagery matter, keep Rust as a WASM compute module and let CesiumJS (Apache-2.0) or MapLibre (globe since 5.0) render, passing image tiles or textures across wasm-bindgen.
- **MSRV.** The viz crates need Rust 1.95–1.97 (egui, egui_plot, walkers 1.95; rerun 1.96; bevy 1.97.1). Keep them in a separate app crate so the libraries can stay at MSRV 1.85 like the geo3d/4d/5d family.

### Gaps
- I did not verify whether Rerun has a native geospatial or map view suitable for gridded fields, or whether deck.gl's GlobeView is still flagged experimental in v9.4.
- The MapLibre GL JS licence was not confirmed, because GitHub reports NOASSERTION. It is commonly BSD-3-Clause; check LICENSE.txt.
- I did not evaluate coastline and graticule data sources (for example Natural Earth) or Rust shapefile/GeoJSON readers for map overlays.
- I did not confirm Crameri's licence on Zenodo directly.

## 6. Crate publishing best practices (workspace, features, no_std, docs.rs, licensing, CI, semver, MSRV)

### Takeaway
Use a Cargo workspace with small, focused library crates and a separate app crate. Dual-license MIT OR Apache-2.0, the ecosystem norm; geo5d is currently Apache-2.0 only. Set `rust-version` and use the MSRV-aware resolver. Gate heavy and C-backed dependencies behind features, and publish from CI with crates.io Trusted Publishing (OIDC). The current stable toolchain is Rust 1.99.1, released 2026-10-13.

### Cited Findings
- Current stable Rust is 1.99.1 (58cf8f9c3, 2026-10-13), and the stable channel manifest is dated 2026-10-15. — [static.rust-lang.org stable channel manifest](https://static.rust-lang.org/dist/channel-rust-stable.toml)
- Rust 1.84.0 stabilized an MSRV-aware resolver. Opt in with `incompatible-rust-versions = "fallback"` in `.cargo/config.toml` or with `resolver = "3"` (requires MSRV ≥ 1.84). Edition 2024 implies resolver 3. A virtual workspace must set `resolver` explicitly. — [Rust 1.84.0 release post](https://blog.rust-lang.org/2025/01/09/Rust-1.84.0/); [Edition guide: Rust-version aware resolver](https://doc.rust-lang.org/edition-guide/rust-2024/cargo-resolver.html)
- crates.io Trusted Publishing (RFC 3691) lets GitHub Actions exchange an OIDC token for a short-lived publish token. It needs `permissions: id-token: write` and a one-time per-crate configuration after the first manual publish. Later updates added GitLab.com support, a "Trusted Publishing Only" mode, and blocked `pull_request_target`/`workflow_run` triggers. — [RFC 3691](https://rust-lang.github.io/rfcs/3691-trusted-publishing-cratesio.html); [crates.io development update, Jan 2026](https://blog.rust-lang.org/2026/01/21/crates-io-development-update); [Alpha-Omega blog](https://alpha-omega.dev/blog/trusted-publishing-secure-rust-package-deployment-without-secrets/)
- The family's current practice:
  - geo5d: edition 2021, `rust-version = "1.85"`, Apache-2.0, `[lints.clippy] all = warn`, CI workflow, CHANGELOG with "plain semver; tags vX.Y.Z", `exclude` for CI files. Package metadata is deliberately not inherited from a workspace so that each repo publishes standalone. — [geo5d Cargo.toml](https://github.com/FusionVault/geo5d/blob/main/Cargo.toml); [CHANGELOG](https://github.com/FusionVault/geo5d/blob/main/CHANGELOG.md)
  - geo3d: `std` default, `libm` for no_std, and optional interop features — [crates.io API geo3d](https://crates.io/api/v1/crates/geo3d)
- Ecosystem licence norm: most foundational crates surveyed are "MIT OR Apache-2.0" (rayon, ndarray, egui, wgpu, rustfft, satkit, sphrs, geo, proj4rs). nalgebra is Apache-2.0 only and faer is MIT only. — [crates.io API data, §2/§5 tables](https://crates.io/crates/rayon)

### Inferences
- **Suggested workspace** (each crate publishable on its own, with versions bumped independently):

  ```
  gravity-workspace/
    crates/sh-core   (e.g. `geopotential` or `shgrav`): no_std + alloc (libm feature, like geo3d), f64 ALFs (Holmes–Featherstone with u^m Horner; X-number feature later), normalization/phase conversions, coefficient containers (static + time-variable), point synthesis, grid synthesis/analysis (features: `fft` = realfft, `rayon`), degree variances.
    crates/io        (e.g. `gravity-field` / `icgem`): ICGEM gfc 2.0 incl. gfct/trnd/acos/asin and tide/norm headers, GRACE L2 GRCOF2+YAML, Love numbers, TN-13/TN-14; features `netcdf-pure` (default, netcdf-reader) / `netcdf-c` (netcdf) / `hdf5`.
    crates/grace     (e.g. `gravimetry` / `grace-l2`): Gaussian/DDK/destriping filters, EWH via load Love numbers, GIA/degree-1/C20 handling, mascon grids, time series.
    geo5d            (existing): depends on sh-core; keeps normal gravity, heights, relativity, transfer; `GravityModel` trait as the boundary.
    apps/viewer      (not necessarily published): eframe/egui + egui_plot + wgpu globe; plotters CLI for static figures. Higher MSRV allowed.
  ```

- **Library hygiene.** Dual-license new crates MIT OR Apache-2.0, and consider relicensing geo5d the same way while it has no users. Add `[package.metadata.docs.rs] all-features = true` so feature-gated APIs are documented. Run `cargo-semver-checks` in CI before each release. Test the MSRV job at 1.85 for the libraries and the stable/beta jobs at 1.99. Keep `#![forbid(unsafe_code)]` as the family does. Return `Result` from parsers and constructors instead of panicking.
- **Data.** Do not embed large coefficient files (EGM2008 is ~4.7 M coefficients, per geo5d's docs) in crates. crates.io packages should stay small. Fetch them at runtime, or ship small test fixtures (EGM96 truncated, a single GRACE month) under `tests/data`, and exclude big files with `exclude`.
- **Versioning.** Stay on 0.x until the SH-core API (coefficient container, normalization enum, grid types) settles. Breaking changes in 0.x are minor bumps and should be recorded in the CHANGELOG, as geo5d already does.

### Gaps
- I did not fetch current Cargo-book or docs.rs pages for the exact `[package.metadata.docs.rs]` keys, the 2026 status of `cargo-semver-checks`, or crates.io's naming and squatting policy text. Those recommendations rest on general knowledge, not on sources retrieved in this session.
- I did not inspect geo5d's CI workflow contents (the `.github/workflows/ci.yml` fetch failed), so I could not confirm whether it tests MSRV or no_std.
