# Roadmap

Where the work stands and what comes next, in the order agreed on 2026-10-24. Work through the
phases in this order; do not start a later phase early.

## Phase A: the structure (roadmap steps 1–3)

The plan is the "What goes into the crates" section of the
[magnetics and MHD report](topics/magnetics-and-mhd/report.md). Done so far:

- `geopotential::Field` (the common field trait) and `geopotential::geocentric` (WGS84);
  `metrology` gained `WGS84_A`, `WGS84_F`, `MU_B`, `GAMMA_P` and `AU`; `gravimetry` gained
  `sea_level_equivalent` and `OCEAN_AREA`.
- The `magnetics` crate: `main_field` (WMM2025 through geopotential's engine; all 100 NOAA test
  rows), `mhd`, `stability`, `flow`, `dynamo`, `space_weather` and `magnets`, all tested.
- Runnable examples, run by `cargo test`, in geopotential, gravimetry, gravitational, gravitics,
  plasmoid and magnetics.

Also done in Phase A: `gravitics::gem` (the gravitomagnetic dipole as a `Field`, three
conventions, Blackett, the extremal charge-to-mass ratio), `plasmoid::tokamak` (q95, Greenwald,
Troyon, IPB98 with ITER and SPARC), and the `physica` umbrella crate (re-exports, seven tested
analogies, the three-fields example). The last duplicated constants now come from `metrology`.

Left out of Phase A: `magnetics::{devices, astro, frontier}`. Writing them was stopped twice
by an automated safety filter, so they wait; their items can be revisited in Phase C at the
level of general published physics.

## Phase B: review the external research (in this session)

Go through every item in both batches in [`external/`](external/README.md) and write a
`<batch>-review.md` beside each, following the method in that README.

**Done 2026-10-24.** Both batches reviewed in full, with a `<batch>-review.md` and a
`<batch>-checks.py` beside each; every item has a verdict. The accepted items, for Phase C:

- `gravitics::electrogravity` (new): dielectrophoretic force, field angular momentum, and the
  mass-independence discriminator — the Pais/Zhang claim as a claim-vs-bound calculator.
- `geopotential::SecondPartials::tensor_neu` (local gradient tensor) and `gravimetry::local`
  (sphere, cylinder, prism forward models).
- `gravitational::laboratory::{quadrupole_strain, burst_energy}` and `metrology::EINSTEIN_KAPPA`
  (8πG/c⁴); a lab-HFGW example.
- `gravitics::anomalies::tension` and the same-site G tests; `magnetics::magnets::susceptibility_force`.
- Ampère line-current force beside the Ostriker filament (galactic dark-matter check).
- `magnetics::space_weather` dynamo EMF; capped AGW growth in `gravimetry::atmosphere`; PEGS
  timing in `gravimetry::coseismic`.
- `gravimetry::atom` phase and Aharonov–Bohm helpers; ultralight-DM mid-band; plasma index and
  the Gordon-metric note.
- To source before coding: chameleon/symmetron screening, the PEGS template amplitudes, the
  Pais patent formula, Peratt 1986, DeWitt 1966.

Two errors in our own notes were found and fixed: the AQG averaging time (625 s, not 6,250 s)
and the "3 km navigation error" (flagged to-source).

The approach to fringe material, as the user set it: **nothing fringe is brushed off.** Many
of these ideas come from established, high-profile scientists, and some may already be verified.
For every claim:

- compute it in the proposal's own terms, from its authors' formulas and inputs;
- compute the standard prediction for the same setup;
- set both against what experiments measured;
- try to verify it: find the primary source, recompute its numbers, and look for any
  replication or measurement that supports it;
- record it in the crates as a calculator, whatever the verdict, labelled with its status.

That is how `gravitics` already works ("each proposal's own model next to general relativity and
the measurements"). Items with dual-use framing (detection, cloaking) are implemented at the
level of the public physics of the effect.

## Phase C: implement what the review produces — done 2026-10-24

Every accepted item is implemented and tested (293 workspace tests pass; clippy and rustdoc
clean):

| Item | Where |
|---|---|
| Einstein coupling `8πG/c⁴` | `metrology::EINSTEIN_KAPPA` |
| Local gravity-gradient tensor (north, east, up), checked against finite differences in a true local frame | `geopotential::SecondPartials::tensor_neu` |
| Buried-body forward models (sphere, horizontal cylinder) and sensor averaging | `gravimetry::local` |
| Wave growth with a dissipation ceiling | `gravimetry::atmosphere::wave_growth` |
| Prompt elastogravity: lead time over the P-wave, and the dipole amplitude | `gravimetry::coseismic::{prompt_lead, prompt_gravity}` |
| Atom-interferometer phase sensitivity and the gravitational Aharonov–Bohm rate | `gravimetry::atom::{phase_sensitivity, aharonov_bohm_rate}` |
| Laboratory quadrupole strain and burst energy | `gravitational::laboratory::{quadrupole_strain, burst_energy}` |
| Plasma critical density and refractive index (the Gordon-metric link) | `gravitational::hf::{critical_density, plasma_index}` |
| Ultralight dark matter: Compton frequency, the mid-band in mass, coherence length | `gravitational::ultralight` |
| Dielectrophoretic force, field angular momentum, the mass-independence discriminator | `gravitics::electrogravity` |
| Pairwise tension, and the test that the G scatter follows the apparatus not the site | `gravitics::anomalies::tension` |
| The DeWitt gravitomagnetic suppression `Gm_e²/(ħc)` | `gravitics::gem::ELECTRON_GRAVITATIONAL_COUPLING` |
| Susceptibility force in a field gradient, and the bias it puts on G | `magnetics::magnets::susceptibility_force` |
| Ionospheric dynamo field | `magnetics::space_weather::dynamo_field` |
| Ampère line-current force beside its gravitational twin (the galactic comparison) | `magnetics::mhd::{line_current_force, line_mass_force}` |

Still to source before coding (unchanged): chameleon/symmetron screening, the PEGS template
amplitudes, the Pais patent's own formula, Peratt 1986, DeWitt 1966, the rectangular-prism
forward model, magnon HFGW limits.

Also still open: `magnetics::{devices, astro, frontier}`, specified in the magnetics report at
lines 181–211 with test values and sources. Four attempts to write them were stopped by an
automated filter, so they are left for the user to add from that specification.

## Phase D: documentation and release — done 2026-10-24

- Five new examples (`gravitics/electrogravity_claim`, `gravimetry/buried_tunnel`,
  `gravitational/lab_hfgw`, `plasmoid/sparc`, `magnetics/storm`), seventeen in all, each run by
  `cargo test`.
- READMEs for `magnetics` and `physica`, an Applications section in every crate README, and the
  root README updated for eight crates.
- Published 2026-10-24: `metrology` 0.2.0, `geopotential` 0.3.0, `gravimetry` 0.6.0,
  `gravitational` 0.2.0, `gravitics` 0.4.0, `plasmoid` 0.3.0, and the new `magnetics` 0.1.0 and
  `physica` 0.1.0, each tagged `<crate>-vX.Y.Z`.

303 workspace tests pass; clippy and rustdoc clean.

## Phase E: Zhang's unified field theory and US 12,417,870 — done 2026-10-25

Asked for by the user: a 3D rebuild of the patent's rig from real physics, the physics the patent
needs, and Zhang's own constant `k`, with every formula in the crates. Reviewed in
[`external/2026-10-24-zhang-unified-field-theory-review.md`](external/2026-10-24-zhang-unified-field-theory-review.md),
with the theory's equations and four rigs computed both ways.

- `gravitics::electrodes` (new): electrodes of rings and straight wire segments, the field and
  charges of any coil-and-lead rig.
- `gravitics::electrogravity`: the sphere's dielectrophoresis, the Pauthenier charge, the
  polarisation-torque ceiling, a thread's torsion, spheroid depolarisation and the torque on a
  body not quite round, the image force, the Thomson–Purcell kink field (Zhang's `E_θ/e_r`), the
  carriers' drift (the theory's own bound), the gravity of field energy, the electric pressure,
  the five-constant identity, and Zhang's `k_mq` (three values), `μ_g` and vortex constant
  `f = k_mq/(4πε₀G)`.
- `plasmoid::discharge` (new): Paschen breakdown, the mean free path, Peek's corona onset, gas
  density and speed, slip-corrected drag on a moving or turning sphere, ion mobility, the
  floating potential.
- `physica/examples/gravity_patent.rs`: the rig simulated at 1 ms steps under every force and
  torque, plus the weighing rig. It is a calibration case, so it is an example, not library code.
- [`web/gravity-patent.html`](../web/gravity-patent.html): the same simulation in 3D, checked
  against the example's pinned values on load.
- Three statements in the batch-2 review and this review's first version corrected from the
  primary sources and the 3D solution (Setup 1's polarity, Setup 2's magnetic field, the
  symmetric drive with leads).
- Found: Zhang's G is 34σ above every CODATA measurement; Paper J's closures do not hold with its
  own numbers; the two closed-form `k_mq` would spin everything in Earth's field.
- Still open: the weighing rig's seven-cable ratio (×7.1 against ×4.4 from a cage model; it
  needs the real layout in 3D), the paper behind the "calibrated" `k_mq`, and the Michelson and
  80 V rigs, which have no published details.

Published 2026-10-25: `gravitics` 0.5.0, `plasmoid` 0.4.0 and `physica` 0.2.0, then patch
releases of the other five for the repository's new name (`metrology` 0.2.1, `geopotential`
0.3.1, `gravimetry` 0.6.1, `gravitational` 0.2.1, `magnetics` 0.1.1), each tagged
`<crate>-vX.Y.Z`.

## Phase F: cross-domain frontiers — planned 2026-10-25

Research into 2022–2026 papers and patents in plasma physics, magnetics, gravity and light:
[`topics/cross-domain-frontiers/report.md`](topics/cross-domain-frontiers/report.md), with seven
notes and eight passing check scripts. It proposes 42 additions, each pinned to a published number,
in six phases after a phase of ten fixes to the published crates. Tracked as GitHub issues in the
[Cross-domain frontiers milestone](https://github.com/FusionVault/physica/milestone/1), with the
build order in [#10](https://github.com/FusionVault/physica/issues/10):

| Issue | Phase |
|---|---|
| [#1](https://github.com/FusionVault/physica/issues/1) | 0: fix the ten defects in the published crates |
| [#2](https://github.com/FusionVault/physica/issues/2) | 1: the shared potential-field core in `geopotential` |
| [#3](https://github.com/FusionVault/physica/issues/3) | 2: measuring gravity on the ground; the equivalence principle |
| [#4](https://github.com/FusionVault/physica/issues/4) | 3: black holes, lensing and light in curved spacetime |
| [#5](https://github.com/FusionVault/physica/issues/5) | 4: plasma kinetics, collisions, ignition, dense matter and lasers |
| [#6](https://github.com/FusionVault/physica/issues/6) | 5: magnetometers, magnets, the main field, optics and aeromagnetics |
| [#7](https://github.com/FusionVault/physica/issues/7) | 6: `physica` analogies and examples |
| [#8](https://github.com/FusionVault/physica/issues/8) | Verify before coding |
| [#9](https://github.com/FusionVault/physica/issues/9) | Deferred items and what would bring them in |

## Phase G: the topological-monopole talk — researched 2026-10-26

Bob Greenyer's TeslaTech 2026 talk, the user's slide of nine names and the user's own audit:
[`topics/topological-monopoles-and-transmutation/report.md`](topics/topological-monopoles-and-transmutation/report.md),
with seven notes and seven passing check scripts. Every claim is computed both ways and labelled.
Tracked in the [Topological monopole talk milestone](https://github.com/FusionVault/physica/milestone/2):
[#27](https://github.com/FusionVault/physica/issues/27) `metrology` first, then [#24](https://github.com/FusionVault/physica/issues/24) `magnetics::monopole` (a new module, the user's
decision) and [#22](https://github.com/FusionVault/physica/issues/22) `plasmoid::nuclear`; [#20](https://github.com/FusionVault/physica/issues/20), [#21](https://github.com/FusionVault/physica/issues/21), [#23](https://github.com/FusionVault/physica/issues/23),
[#25](https://github.com/FusionVault/physica/issues/25) and [#26](https://github.com/FusionVault/physica/issues/26) by area; [#28](https://github.com/FusionVault/physica/issues/28), the one example and claim registry, last.

## One build order (2026-10-26)

Every report, review and PR was swept against the crates and the issues on 2026-10-26. What no
issue tracked is now in the [Research backlog milestone](https://github.com/FusionVault/physica/milestone/3)
([#29](https://github.com/FusionVault/physica/issues/29)–[#41](https://github.com/FusionVault/physica/issues/41)), and the build order for all open work, across the three milestones,
is [#10](https://github.com/FusionVault/physica/issues/10). Issues carry `P1`–`P3`, `crate: …`, `claims`, `blocked`, `on hold` and `tracking` labels.

## Open

- **Retired crates:** `satgrav`, `gravity-field` and `gravity-claims` are deleted; `gwsens` was
  yanked by the user on 2026-10-24 (deletion was blocked because `gravitics` 0.1.0 and 0.2.0
  depend on it).
- **Decisions for the user** ([#30](https://github.com/FusionVault/physica/issues/30)): whether to download IGRF-14's coefficient file and
  EGM96 for the IGRF and geoid examples; whether the retired `gravity-field` readers come back
  ([#35](https://github.com/FusionVault/physica/issues/35)).
- **On hold at the user's request:** `magnetics::{devices, astro, frontier}`, every specified
  calculator listed in [#42](https://github.com/FusionVault/physica/issues/42). The specification is in
  [`topics/magnetics-and-mhd/report.md`](topics/magnetics-and-mhd/report.md). Its monopole rows
  moved to `magnetics::monopole` ([#24](https://github.com/FusionVault/physica/issues/24)).

## Later

- W5's pre-registration is time-critical ([#29](https://github.com/FusionVault/physica/issues/29)).
- The PPN calculators, once planned for `geo5d`, are now `gravitational::ppn` ([#32](https://github.com/FusionVault/physica/issues/32)).
- The "verify before coding" items of every report are on [#8](https://github.com/FusionVault/physica/issues/8).
- W1 scoring when COST-G's December 2026 fields appear, about March 2027 ([#41](https://github.com/FusionVault/physica/issues/41)).
