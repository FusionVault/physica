# metrology

The physical constants the quantum-gravimeter crates share, defined once: the exact SI defining
constants, the CODATA 2018 values of G, α, μ₀, ε₀ and the particle masses, standard gravity, and
the astronomical units of mass, distance and time.

## Install

```bash
cargo add metrology
```

## Example

```rust
use metrology::{C, G, HBAR, SOLAR_MASS};

// The Schwarzschild radius of the Sun, 2GM/c²: about 2.95 km.
let r = 2.0 * G * SOLAR_MASS / (C * C);
assert!((r / 2953.3 - 1.0).abs() < 1e-3);

// The Planck length, √(ħG/c³).
assert!(((HBAR * G / C.powi(3)).sqrt() / metrology::PLANCK_LENGTH - 1.0).abs() < 1e-5);
```

## What is inside

| Group | Constants | Source |
|---|---|---|
| SI defining constants (exact) | `C`, `H`, `HBAR`, `E_CHARGE`, `K_B`, `N_A`, `R_GAS` | [SI Brochure 2019](https://www.bipm.org/en/publications/si-brochure) |
| Measured | `G`, `ALPHA`, `MU_0`, `EPS_0`, `M_E`, `M_P`, `AMU`, `M_H`, `PLANCK_LENGTH` | [CODATA 2018](https://physics.nist.gov/cuu/Constants/) |
| Energy | `EV`, `KEV` | |
| Earth and astronomy | `G0`, `GM_EARTH`, `GM_SUN`, `SOLAR_MASS`, `PARSEC`, `MPC`, `YEAR` | [IERS Conventions 2010](https://www.iers.org/TN36); [IAU 2015 B3](https://www.iau.org/static/resolutions/IAU2015_English.pdf) |

## How it is checked

A test checks that the derived constants agree with each other: `ħ = h/2π`, `R = N_A k_B`,
`μ₀ε₀c² = 1`, `α = e²/4πε₀ħc`, the Planck length, and the solar mass.

One of six crates in [quantum-gravimeter](https://github.com/FusionVault/quantum-gravimeter):
[`metrology`](https://crates.io/crates/metrology) (constants) ·
[`geopotential`](https://crates.io/crates/geopotential) (the gravity field) ·
[`gravimetry`](https://crates.io/crates/gravimetry) (missions and GRACE) ·
[`gravitational`](https://crates.io/crates/gravitational) (waves) ·
[`gravitics`](https://crates.io/crates/gravitics) (frontier gravity) ·
[`plasmoid`](https://crates.io/crates/plasmoid) (fusion plasma).
