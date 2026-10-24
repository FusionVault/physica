# gravitics

*Gravitics* is what the 1950s called the search for gravity control. This crate computes that
whole parallel track of physics both ways:
- **The proposal's own model,** from its authors' formulas and inputs, reproducing their published
  numbers: superconductor gravity, high-frequency gravitational-wave generators and detectors,
  propellantless thrusters, warp drives and the polarizable vacuum.
- **The standard prediction,** from general relativity or standard physics, for the same setup.
- **What experiments measured.**

Use it to explore a proposal, to design an experiment that could test one, or to see how large
an effect must be before an instrument can see it.

![Each proposal's own numbers against general relativity and experiment](https://raw.githubusercontent.com/FusionVault/quantum-gravimeter/main/figures/claims_gap.svg)

## Install

```bash
cargo add gravitics
```

## Example: one proposal, both tracks

```rust
use gravitics::superconductors::{TAJMAR_REPORTED, ring_frame_dragging};

// Tajmar & de Matos's spinning niobium ring (0.349 kg, radii 75 and 69 mm).
// Their measurement: the frame turned with 2.6×10⁻⁵ of the ring's rotation.
// General relativity, for the same ring: 7.2×10⁻²⁷.
let gr = ring_frame_dragging(0.349, 0.075, 0.069);
println!("reported {TAJMAR_REPORTED:.1e}, general relativity {gr:.1e}"); // 2.6e-5, 7.2e-27
assert!(TAJMAR_REPORTED / gr > 1e21);
```

```rust
use gravitics::propulsion::{photon_thrust, shawyer_thrust};

// An EmDrive cavity at Dresden's settings: Shawyer's model against a photon rocket.
let model = shawyer_thrust(11.0, 741.0, 0.56); // 30 µN
let light = photon_thrust(11.0);               // 0.037 µN
assert!(model / light > 800.0); // Dresden measured less than 0.03 µN
```

## The proposals

| Proposal | Its own model or report | General relativity or standard physics | What experiments found |
|---|---|---|---|
| Tajmar: frame dragging near a spinning superconductor | k = 2.6×10⁻⁵ | 7.2×10⁻²⁷ | Tajmar's own later limit is 800× below the first signal; traced to vibration |
| Podkletnov: weight change above a spinning superconductor | up to −2.1% | a weight increase of at most 3×10⁻⁹ | no change above 2×10⁻⁸ (NASA Marshall) |
| Li–Baker: a GHz gravitational-wave detector | a signal in 3×10³ s | one signal photon in 3×10¹⁸ s, against the reference beam's shot noise (JASON's argument) | none built |
| Pais: 10²⁴ W/m² from a charged disc | σ = 5×10⁴ C/m² | a surface field 5.6×10⁶ times the breakdown limit his patent states | Navy test (2016–19): effect not shown |
| EmDrive | 1.2 mN/kW (Eagleworks) | a photon rocket gives 3.3 µN/kW | nothing above 10–30 nN (Dresden) |
| Mach-effect thruster | "a couple of µN" | the authors' own formula in SI: about 0.1 µN | vibration artefacts (Dresden) |
| Lifters ("electrogravitics") | lift | ion wind, `F = I·d/μ` | no force in vacuum |
| Warp drive (Alcubierre) | travel faster than light | negative energy of a galaxy or more, in walls quantum field theory allows | — |
| Polarizable vacuum (Puthoff) | gravity as a vacuum index | agrees with general relativity at first order, differs at second | 0.73 µas apart for light grazing the Sun |

## What is inside

| Module | What it computes |
|---|---|
| `hfgw` | High-frequency gravitational waves: the Li–Baker detector (flux conventions, graviton rate, the interference signal, the shot-noise time) and the Pais patents (radiated flux, surface field, magnetic conversion as he uses it and in SI) |
| `superconductors` | The London moment; frame dragging of a ring in general relativity against Tajmar's report; the Cooper-pair density and gravitational London depth; a perfect shield's shadow; Podkletnov's pendulum and Marx bank; the measured bounds |
| `propulsion` | Ion wind and ion mobility, the vacuum (ballistic) limit, the photon rocket, Shawyer's EmDrive formula, buoyancy on a balance, Woodward's Mach-effect thrust |
| `spacetime` | Alcubierre warp-bubble energy ([Pfenning & Ford 1997](https://arxiv.org/abs/gr-qc/9702026)), the quantum inequality, Casimir plates and spheres, wormhole tensions, light bending in Puthoff's polarizable vacuum against general relativity, and horizon temperatures (Unruh, Hawking, Gibbons–Hawking) with black-hole lifetimes and zero-point energy |
| `strong_field` | The Schwinger limit: the critical field, the pair-production and vacuum-decay rates, laser a₀ and the quantum parameter χ (with explicit intensity conventions), Ritus–Narozhny, Breit–Wheeler suppression, vacuum birefringence, Euler–Heisenberg, supercritical atoms, and the graphene analogue. The Schwinger exponent is exactly a Boltzmann factor at twice the Unruh temperature, and a test pins that identity |

![The energy a warp bubble needs](https://raw.githubusercontent.com/FusionVault/quantum-gravimeter/main/figures/warp_energy.svg)

*A warp bubble's negative energy grows as its wall thins. Quantum field theory allows walls only
about 100 Planck lengths thick, where a 100 m bubble needs billions of times the ordinary matter
in the observable universe.*

![Thrust per watt](https://raw.githubusercontent.com/FusionVault/quantum-gravimeter/main/figures/thrust_per_watt.svg)

*Thrust per watt: the propellantless proposals, a photon rocket, ion wind, and the best
measurement.*

The research behind every number, with its sources and the status of each proposal, is in the
[research report](https://github.com/FusionVault/quantum-gravimeter/blob/main/reports/Gravity%20beams%20fringe%20physics%20and%20DIRDs.md).

## How it is checked

Every test reproduces a published number from the paper's own inputs, mostly to better than 1%:
- the Li–Baker graviton flux and Table I signal;
- the Pais disc and conversion numbers;
- Tajmar's predicted field for niobium;
- Podkletnov's pendulum tables;
- Woodward's printed 1.87×10²¹ prefactor;
- Shawyer's thrust at the Dresden and Eagleworks settings;
- the Army lifter thrust;
- Pfenning & Ford's warp energy;
- the quantum inequality;
- the Casimir pressure at 1 µm;
- the Sun's 1.7512″ light bending;
- the Schwinger field and pair rates, CoReLS's a₀, LUXE's χ, PVLAS's QED birefringence, and
  Dirac binding for lead and uranium;
- Unruh, Hawking and Gibbons–Hawking temperatures.

The formulas were independently re-derived from their sources and every test value recomputed;
the Casimir formulas are also checked against each other.

Formerly published as `gravity-claims`. Part of
[quantum-gravimeter](https://github.com/FusionVault/quantum-gravimeter).
