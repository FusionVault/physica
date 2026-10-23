# geopotential

Spherical-harmonic gravity fields in pure Rust, with no required dependencies (optional `serde`). It is the engine behind
geoid, gravity and GRACE water-storage maps, and behind error spectra for simulated satellite
missions.

- **Legendre functions:** fully normalized (4π, no Condon–Shortley phase), with first and
  second latitude derivatives, stable to degree 2700 at every latitude (Holmes & Featherstone 2002)
  and exact at the poles.
- **Coefficients:** with their `GM`, radius and tide system, plus degree variances and amplitudes,
  tide-system conversion (IERS 2010), `GM`/radius rescaling and truncation.
- **Time-variable models:** ICGEM `gfct`/`trnd`/`acos`/`asin` terms evaluated at any epoch.
- **Synthesis:** potential with its gradient and Hessian at a point, and grids of geoid height,
  gravity disturbance and anomaly, radial gradient or equivalent water height.

```rust
use geopotential::{Coefficients, Functional};

// A J2-only field: the geoid it implies at the equator and the pole.
let mut j2 = Coefficients::zeros(3.986004415e14, 6378136.3, 2);
j2.set(2, 0, -4.84165e-4, 0.0);
let n = j2.grid(&Functional::GeoidHeight, j2.radius, &[0.0, 90f64.to_radians()], &[0.0]);
assert!(n[0] > 0.0 && n[1] < 0.0);
```

Read ICGEM `.gfc` files with [`gravity-field`](https://crates.io/crates/gravity-field).

**Validation.** The tests cover closed forms, the addition theorem to degree 2700 near the pole,
finite differences, Laplace's equation, and the IERS C̄20 tide offset. Separately, synthesis of
the real degree-300 GOCO06s model agreed with an independent implementation to about 1e-15 in
potential and 1e-6 µGal in gravity.

Licensed under either of MIT or Apache-2.0 at your option.
