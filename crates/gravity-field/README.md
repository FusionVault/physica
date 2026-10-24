# gravity-field

Read the files that satellite gravity fields come in:
- ICGEM `.gfc` models, from static models like GOCO06s to every monthly GRACE field;
- the GRACE technical notes that correct the lowest degrees.

You get [`geopotential`](https://crates.io/crates/geopotential) coefficients, ready to evaluate.

![The error of every coefficient of the July 2026 GRACE-FO field](https://raw.githubusercontent.com/FusionVault/quantum-gravimeter/main/figures/error_triangle_2026-07.svg)

*The calibrated error of each coefficient of the July 2026 COST-G field, read with this crate.*

## Install

```bash
cargo add gravity-field
```

## Example

Download a model first, for example `GOCO06s.gfc` from the
[ICGEM model table](https://icgem.gfz.de/tom_longtime).

```rust,no_run
use gravity_field::{julian_epoch, read_gfc};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let text = std::fs::read_to_string("GOCO06s.gfc")?;
    let gfc = read_gfc(&text)?;

    // Time-variable models apply their trends and annual terms; static ones ignore the date.
    let field = gfc.model.at(julian_epoch(2026, 7, 16, 12, 0));

    println!("{} to degree {}", gfc.name.unwrap_or_default(), field.max_degree());
    println!("C20 = {:e}", field.c(2, 0));
    Ok(())
}
```

## What is inside

| Function | Reads |
|---|---|
| `read_gfc` | ICGEM `.gfc` models, static or time-variable (icgem1.0 and 2.0), with their errors ([format](https://icgem.gfz.de/docs/ICGEM-Format-2023.pdf)) |
| `read_tn14` | NASA's TN-14: C̄20 and C̄30 from satellite laser ranging ([Loomis et al. 2020](https://doi.org/10.1029/2019GL085488)) |
| `read_tn13` | NASA's TN-13: degree 1, the motion of Earth's centre of mass ([Sun et al. 2016](https://doi.org/10.1002/2016JB013073)) |
| `julian_epoch` | Dates in the decimal years these files use |

The files are free:
- **Models:** from [ICGEM](https://icgem.gfz.de). Each COST-G monthly series also comes as one
  zip.
- **Technical notes:** from NASA's
  [PO.DAAC](https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-docs/gracefo/open/docs/).

## How it is checked

- **Real files** from COST-G, GOCO06s and NASA. These bend the format in practice: missing header
  markers, comment and warning lines, `errors empirical`, and free text before the header.
- **Errors:** a malformed record fails with its line number.
- **Fixtures:** their sources are listed in `tests/data/README.md`.

Part of [quantum-gravimeter](https://github.com/FusionVault/quantum-gravimeter).
