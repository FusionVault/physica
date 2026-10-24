# geopotential

Earth's gravity field from spherical-harmonic coefficients. Evaluate the potential, gravity
vector or gravity-gradient tensor at any point, or draw maps of the geoid, gravity anomalies or
water height. Pure Rust, no required dependencies.

![The geoid from the GOCO06s satellite model](https://raw.githubusercontent.com/FusionVault/quantum-gravimeter/main/figures/geoid_goco06s.svg)

*The geoid from GOCO06s, drawn with this crate.*

## Install

```bash
cargo add geopotential
```

## Example

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

To load a real model such as GOCO06s or a monthly GRACE field, use
[`gravity-field`](https://crates.io/crates/gravity-field).

## What is inside

| Item | What it does |
|---|---|
| `Legendre` | Normalized Legendre functions and their derivatives, stable to degree 2700 at every latitude ([Holmes & Featherstone 2002](https://doi.org/10.1007/s00190-002-0216-2)) |
| `Coefficients` | The coefficients with their `GM`, radius and tide system. Also degree variances, rescaling, truncation and tide-system conversion ([IERS Conventions 2010](https://www.iers.org/TN36), ch. 6) |
| `TimeVariable` | Models with trends and periodic terms, evaluated at any date ([ICGEM](https://icgem.gfz.de)) |
| `partials`, `second_partials` | Potential, gravity vector and gravity-gradient tensor at a point |
| `grid` with `Functional` | Maps of geoid height, gravity anomaly or disturbance, radial gradient or equivalent water height ([Barthelmes 2013](https://doi.org/10.2312/GFZ.b103-0902-26); [Wahr et al. 1998](https://doi.org/10.1029/98JB02844)) |

**Conventions:**
- angles in radians, geocentric latitude;
- fully normalized coefficients with no Condon–Shortley phase, as in ICGEM files;
- the optional `serde` feature serializes `Coefficients`.

## How it is checked

- Legendre functions against closed forms, the addition theorem to degree 2700 close to the
  pole, and the Legendre equation.
- Derivatives against finite differences.
- The field against Laplace's equation.
- The IERS tide offset of C̄20.
- The real degree-300 GOCO06s model agrees with an independent implementation
  ([geo5d](https://crates.io/crates/geo5d)) to about 1e-15.

Part of [quantum-gravimeter](https://github.com/FusionVault/quantum-gravimeter), with
[`gravity-field`](https://crates.io/crates/gravity-field) (file readers),
[`gravimetry`](https://crates.io/crates/gravimetry) (GRACE processing) and
[`satgrav`](https://crates.io/crates/satgrav) (missions and quantum sensors).
