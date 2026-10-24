# gravity-claims

Reality checks for extraordinary gravity claims. Each function reproduces a published number from
its authors' own formula, and sits next to what general relativity, standard physics or the
replications say. The tests pin both, so you can see exactly how far off each claim is.

Nothing here supports the claims. Most were made by people with real publications, which is why
their arithmetic can be checked at all.

## Install

```bash
cargo add gravity-claims
```

## Example: how far off are they?

```rust
use gravity_claims::{hfgw, propulsion, superconductors as sc};

// Tajmar's spinning niobium ring: claimed frame dragging against general relativity.
let gr = sc::ring_frame_dragging(0.349, 0.075, 0.069);
assert!(sc::TAJMAR_CLAIM / gr > 1e21);

// The EmDrive: Eagleworks' 1.2 mN/kW against a perfect photon rocket.
assert!(1.2e-3 / propulsion::photon_thrust(1e3) > 300.0);

// The Li–Baker detector: one signal photon in 3.9×10¹⁷ s, not the 3×10³ s claimed.
assert!(hfgw::shot_noise_time(2.56e-18, 2.0) / 3e3 > 1e14);
```

## The claims

| Claim | Published | Standard physics or replication | Gap |
|---|---|---|---|
| Tajmar: frame dragging near a spinning superconductor | k = 2.6×10⁻⁵ | General relativity: 7.2×10⁻²⁷ | 3.6×10²¹ |
| Podkletnov: weight loss above a spinning superconductor | up to 2.1% | No change above 2×10⁻⁸ (NASA Marshall) | 10⁶ |
| Li–Baker: GHz gravitational-wave detector | signal in 3×10³ s | One photon in 3.9×10¹⁷ s (JASON) | 10¹⁴ in time |
| Pais: 10²⁴ W/m² from a charged disc | σ = 5×10⁴ C/m² | 10⁹ V/m breakdown, his own figure | 5.6×10⁶ in field |
| EmDrive (Eagleworks) | 1.2 mN/kW | No thrust above 10–30 nN (Dresden); photon rocket 3.3 µN/kW | 360× light |
| Mach-effect thruster (Woodward) | ~100 µN | His own formula: ~7 nN; vibration artefacts (Dresden) | 10⁴ |
| Lifters ("electrogravitics") | lift | Ion wind, `F = I·d/μ`; nothing in vacuum | none needed |

## What is inside

| Module | What it checks |
|---|---|
| `hfgw` | The Li–Baker detector (wave flux conventions, graviton rate, the interference "signal", JASON's shot-noise time) and the Pais patents (radiated flux, surface field, the Gertsenshtein conversion as he uses it and in SI) |
| `superconductors` | The London moment; frame dragging of a ring in general relativity against Tajmar's claim; the Cooper-pair density and gravitational London depth; a perfect shield's shadow; Podkletnov's pendulum and Marx bank; the replication bounds |
| `propulsion` | Ion wind and ion mobility, the ballistic (vacuum) limit, the photon rocket, Shawyer's EmDrive formula, buoyancy artefacts, Woodward's Mach-effect thrust |
| `spacetime` | Alcubierre warp-bubble energy ([Pfenning & Ford 1997](https://arxiv.org/abs/gr-qc/9702026)), the quantum inequality, Casimir plates and spheres, wormhole tensions, and light bending in Puthoff's polarizable vacuum against general relativity |

The research behind every number, with its sources and status (mainstream, contested or
refuted), is in the
[research report](https://github.com/FusionVault/quantum-gravimeter/blob/main/reports/Gravity%20beams%20fringe%20physics%20and%20DIRDs.md).

## How it is checked

Every test reproduces a published number from the paper's own inputs, mostly to better than 1%:
- the Li–Baker graviton flux and Table I signal;
- JASON's integration time;
- the Pais disc and conversion numbers;
- Tajmar's ρ*/ρ;
- Podkletnov's pendulum tables;
- Woodward's 7 nN;
- Shawyer's thrust at the Dresden and Eagleworks settings;
- the Army lifter thrust;
- Pfenning & Ford's warp energy;
- the quantum inequality;
- the Casimir pressure at 1 µm;
- the Sun's 1.7512″ light bending.

The Casimir formulas are also checked against each other (pressure as the derivative of energy,
and the proximity-force sphere).

Part of [quantum-gravimeter](https://github.com/FusionVault/quantum-gravimeter).
