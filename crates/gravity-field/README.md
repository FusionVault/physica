# gravity-field

Readers for satellite gravity-field files, producing
[`geopotential`](https://crates.io/crates/geopotential) coefficient sets.

- **`read_gfc`:** ICGEM `.gfc` models, static or time-variable (icgem1.0 and icgem2.0), with
  their formal or calibrated errors. It is tested on real COST-G, Swarm and GOCO06s files, which
  bend the format spec: missing `begin_of_head`, `CMMNT` and `WARNING:` lines, `errors empirical`,
  and free text before the header. Malformed records fail with their line number.
- **`read_tn14`, `read_tn13`:** the GRACE / GRACE-FO technical-note series that replace C̄20
  and C̄30 (satellite laser ranging) and supply degree 1 (geocentre).
- **`julian_epoch`:** the time scale (years of 365.25 days) these use.

```rust,no_run
let text = std::fs::read_to_string("GOCO06s.gfc").unwrap();
let gfc = gravity_field::read_gfc(&text).unwrap();
let july_2026 = gfc.model.at(gravity_field::julian_epoch(2026, 7, 16, 12, 0));
println!("C20 = {:e}", july_2026.c(2, 0));
```

Models are free from [ICGEM](https://icgem.gfz.de) (each COST-G series also comes as one zip),
and the technical notes from the
[PO.DAAC docs bucket](https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-docs/gracefo/open/docs/).

Licensed under either of MIT or Apache-2.0 at your option. The test fixtures are real data
excerpts; their sources and licences are in `tests/data/README.md`.
