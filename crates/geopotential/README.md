# geopotential

Earth's gravity field from spherical-harmonic coefficients. Read a model from an ICGEM `.gfc`
file, from static models such as GOCO06s to every monthly GRACE field. Then evaluate the
potential, the gravity vector or the gravity-gradient tensor at any point, or draw maps of the
geoid, gravity anomalies or water height. Pure Rust, no required dependencies.

![The geoid from the GOCO06s satellite model](https://raw.githubusercontent.com/FusionVault/physica/main/figures/geoid_goco06s.svg)

*The geoid from GOCO06s, drawn with this crate.*

## Install

```bash
cargo add geopotential
```

## Example: a field and what it does

```rust
use geopotential::{Coefficients, Functional};

// A field with a single term: J2, Earth's flattening.
let mut field = Coefficients::zeros(3.986004415e14, 6_378_136.3, 2);
field.set(2, 0, -4.84165e-4, 0.0);

// Geoid height at the equator and at the North Pole.
let lats = [0.0, 90f64.to_radians()];
let geoid = field.grid(&Functional::GeoidHeight, field.radius, &lats, &[0.0]);
assert!(geoid[0] > 0.0); // the equatorial bulge
assert!(geoid[1] < 0.0); // flatter at the poles

// Add the central term (C̄00 = 1) and find gravity on a satellite 400 km above northern Italy.
field.set(0, 0, 1.0, 0.0);
let r = 6_378_137.0 + 400e3;
let (lat, lon) = (45f64.to_radians(), 7f64.to_radians());
let [north, east, up] = field.partials(r, lat, lon).acceleration_neu(r, lat);
assert!(up < -8.0); // pointing down, about 8.7 m/s²
assert!(north.abs() > 0.0 && east == 0.0);
```

## Example: read a real model

Download a model first, for example `GOCO06s.gfc` from the
[ICGEM model table](https://icgem.gfz.de/tom_longtime).

```rust,no_run
use geopotential::{julian_epoch, read_gfc};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let gfc = read_gfc(&std::fs::read_to_string("GOCO06s.gfc")?)?;

    // Time-variable models apply their trends and annual terms; static ones ignore the date.
    let field = gfc.model.at(julian_epoch(2026, 7, 16, 12, 0));
    println!("{} to degree {}, C20 = {:e}", gfc.name.unwrap_or_default(), field.max_degree(), field.c(2, 0));
    Ok(())
}
```

![The error of every coefficient of the July 2026 GRACE-FO field](https://raw.githubusercontent.com/FusionVault/physica/main/figures/error_triangle_2026-07.svg)

*The calibrated error of each coefficient of the July 2026 COST-G field, read with `read_gfc`.*

## What is inside

| Item | What it does |
|---|---|
| `Legendre` | Normalized Legendre functions and their derivatives, stable to degree 2700 at every latitude ([Holmes & Featherstone 2002](https://doi.org/10.1007/s00190-002-0216-2)) |
| `Coefficients` | The coefficients with their `GM`, radius and tide system; degree variances, rescaling, truncation and tide-system conversion ([IERS Conventions 2010](https://www.iers.org/TN36), ch. 6) |
| `TimeVariable` | Models with trends and periodic terms, evaluated at any date ([ICGEM](https://icgem.gfz.de)) |
| `partials`, `second_partials` | Potential, gravity vector and gravity-gradient tensor at a point |
| `grid` with `Functional` | Maps of geoid height, gravity anomaly or disturbance, radial gradient or equivalent water height ([Barthelmes 2013](https://doi.org/10.2312/GFZ.b103-0902-26); [Wahr et al. 1998](https://doi.org/10.1029/98JB02844)) |
| `read_gfc`, `Gfc` | ICGEM `.gfc` models, static or time-variable (icgem1.0 and 2.0), with their errors ([format](https://icgem.gfz.de/docs/ICGEM-Format-2023.pdf)) |
| `julian_epoch`, `parse_date` | The decimal-year dates these files use |
| `special` | Complete elliptic integrals `K` and `E` by the arithmetic–geometric mean, and Bulirsch's general complete integral `cel` (Bulirsch 1969), which the fields of rings of mass, charge and current share |
| `linalg` | Pivoted Gaussian elimination for general systems and Cholesky for the normal equations of least squares |
| `fft` | The radix-2 fast Fourier transform, for time series and gridded fields |

**Conventions:**
- angles in radians, geocentric latitude;
- fully normalized coefficients with no Condon–Shortley phase, as in ICGEM files;
- the optional `serde` feature serializes `Coefficients`.

Models are free from [ICGEM](https://icgem.gfz.de); each COST-G monthly series also comes as one
zip.

## Applications

- **Height systems.** Turning a GNSS ellipsoidal height into a usable one needs the geoid; this
  is the synthesis behind national models such as GEOID2022 and AUSGeoid2020.
- **Orbit design and determination.** `J2` from a real model sets the nodal precession that makes
  an orbit sun-synchronous; `examples/sun_synchronous_orbit.rs` reproduces Landsat's 98.2°.
- **Inertial navigation.** The deflection of the vertical is what an accelerometer mistakes for
  horizontal acceleration, and compensating it from a model is standard practice.
- **Gradiometry.** `SecondPartials::tensor_neu` gives the tidal tensor a gradiometer measures, in
  which a vehicle's own acceleration cancels.
- **Reading published models.** `read_gfc` takes any ICGEM `.gfc` file, static or time-variable.

## How it is checked

- **Legendre functions:** against closed forms, the addition theorem to degree 2700 close to the
  pole, and the Legendre equation.
- **Derivatives:** against finite differences.
- **The field:** against Laplace's equation, and the IERS tide offset of C̄20.
- **An independent implementation:** the degree-300 GOCO06s model agrees with
  [geo5d](https://crates.io/crates/geo5d) to about 1e-15.
- **Real files** from COST-G and GOCO06s, which bend the format in practice: missing header
  markers, comment and warning lines, `errors empirical`, and free text before the header. A
  malformed record fails with its line number. The fixtures' sources are listed in
  `tests/data/README.md`.

One of eight crates in [physica](https://github.com/FusionVault/physica):
[`physica`](https://crates.io/crates/physica) (all of them, plus cross-domain analogies) ·
[`metrology`](https://crates.io/crates/metrology) (constants) ·
[`geopotential`](https://crates.io/crates/geopotential) (the gravity field) ·
[`gravimetry`](https://crates.io/crates/gravimetry) (missions and GRACE) ·
[`gravitational`](https://crates.io/crates/gravitational) (waves) ·
[`gravitics`](https://crates.io/crates/gravitics) (frontier gravity) ·
[`magnetics`](https://crates.io/crates/magnetics) (magnetism and MHD) ·
[`plasmoid`](https://crates.io/crates/plasmoid) (fusion plasma).
