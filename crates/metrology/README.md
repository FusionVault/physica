# metrology

The physical constants the physica crates share, defined once: the exact SI defining
constants, the CODATA 2018 values of G, α, μ₀, ε₀ and the particle masses, the flux quanta,
standard gravity, and the astronomical units of mass, distance and time; plus the one-line
laws of radioactive decay.

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
| Measured | `G`, `ALPHA`, `MU_0`, `EPS_0`, `M_E`, `M_P`, `AMU`, `M_H`, `MU_B`, `MU_N`, `GAMMA_P`, `PLANCK_LENGTH` | [CODATA 2018](https://physics.nist.gov/cuu/Constants/) |
| Flux quanta (exact) | `FLUX_QUANTUM` (h/2e), `JOSEPHSON_CONSTANT` (2e/h), `DIRAC_FLUX` (h/e, one Dirac monopole) | [SI Brochure 2019](https://www.bipm.org/en/publications/si-brochure); [PDG monopole review](https://pdg.lbl.gov/) |
| Energy | `EV`, `KEV` | |
| Molecular binding | `WATER_DIMER_D0` (one hydrogen bond, 13.22 kJ/mol) | [Rocher-Casterline et al. 2011](https://doi.org/10.1063/1.3598339) |
| `decay` | `activity`, `activity_step_after_suppression`, `ec_half_life_shift`, `CURIE` | the decay law; [Pommé et al. 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC5207040/) and [Ohtsuki et al. 2004](https://doi.org/10.1103/PhysRevLett.93.112501) for the tests |
| Derived couplings | `EINSTEIN_KAPPA` (8πG/c⁴) | [Einstein field equations](https://en.wikipedia.org/wiki/Einstein_field_equations) |
| Earth and astronomy | `G0`, `GM_EARTH`, `EARTH_ROTATION`, `MOON_EARTH_MASS_RATIO`, `GM_MOON`, `WGS84_A`, `WGS84_F`, `GM_SUN`, `SOLAR_MASS`, `AU`, `PARSEC`, `MPC`, `YEAR` | [IERS Conventions 2010](https://www.iers.org/TN36); [IAU 2015 B3](https://www.iau.org/static/resolutions/IAU2015_English.pdf) |
| Reference ellipsoids and potentials | `GRS80_F`, `GRS80_GM`, `W0_IAU`, `W0_IHRS`, `L_G`, `L_B`, `EARTH_ANGULAR_MOMENTUM` | [Moritz 2000](https://doi.org/10.1007/s001900050278); IAU 2000 B1.9 and 2006 B3; [IAG 2015 Resolution 1](https://iag.dgfi.tum.de/fileadmin/IAG-docs/IAG_Resolutions_2015.pdf) |

## Applications

Anywhere a calculation must agree with another: the crates in this family all take their
constants from here, so a result computed through two of them matches to the last digit. Pin the
version to pin the constants, which is what reproducing a published number needs.

## How it is checked

A test checks that the derived constants agree with each other: `ħ = h/2π`, `R = N_A k_B`,
`μ₀ε₀c² = 1`, `α = e²/4πε₀ħc`, the Planck length, the solar mass, both magnetons, and the flux
quanta (`h/e = 2Φ₀`, and μ₀ times the Dirac charge 3.291×10⁻⁹ A m). The decay tests reproduce
radium's 0.988 Ci per gram, ⁶⁰Co's 7.2×10⁻⁴ step after a two-day stop, and ⁷Be's 0.83%
half-life shift in C₆₀.

One of eight crates in [physica](https://github.com/FusionVault/physica):
[`physica`](https://crates.io/crates/physica) (all of them, plus cross-domain analogies) ·
[`metrology`](https://crates.io/crates/metrology) (constants) ·
[`geopotential`](https://crates.io/crates/geopotential) (the gravity field) ·
[`gravimetry`](https://crates.io/crates/gravimetry) (missions and GRACE) ·
[`gravitational`](https://crates.io/crates/gravitational) (waves) ·
[`gravitics`](https://crates.io/crates/gravitics) (frontier gravity) ·
[`magnetics`](https://crates.io/crates/magnetics) (magnetism and MHD) ·
[`plasmoid`](https://crates.io/crates/plasmoid) (fusion plasma).
