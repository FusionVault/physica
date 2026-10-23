# quantum-gravimeter

Rust crates for satellite gravimetry. They compute geoid, gravity and water-storage fields from open
GRACE / GRACE-FO / GOCE / Swarm data, and the same error spectra and maps for simulated
quantum (cold-atom) sensors.

| Crate | Status | Scope |
|---|---|---|
| [`geopotential`](crates/geopotential) | 0.1.0, unpublished | Legendre functions to degree 2700, coefficients with tide systems and ICGEM time-variable terms, point and grid synthesis, functionals, degree spectra |
| `gravity-field` | planned | ICGEM `.gfc`, GRACE Level-2, TN-13/14 and mascon NetCDF readers; data fetchers |
| `gravimetry` | planned | GRACE processing chain: C20/degree-1/GIA, filters, water height, basins |
| `satgrav` | planned | Sensor noise models (electrostatic, laser ranging, atom interferometry), error propagation |

The research behind this, covering open data access, quantum missions, formulas and the plan, is in
[`reports/`](reports/Satellite%20gravimetry%20data%20and%20Rust%20tooling.md), with source notes in
[`research_notes/`](research_notes/).

```bash
cargo test
```

Licensed under either of [MIT](LICENSE-MIT) or [Apache-2.0](LICENSE-APACHE) at your option.
