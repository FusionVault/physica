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

## Phase E: Zhang's unified field theory and US 12,417,870 — done 2026-10-24

Asked for by the user: a 3D rebuild of the patent's rig from real physics, and the physics the
patent needs. Reviewed in
[`external/2026-10-24-zhang-unified-field-theory-review.md`](external/2026-10-24-zhang-unified-field-theory-review.md),
with the theory's 18 equations and four rigs computed both ways.

- `gravitics::electrodes` (new): coaxial ring electrodes, the field and charges of coil and wire
  rigs.
- `gravitics::electrogravity`: the sphere's dielectrophoresis, the Pauthenier charge, the
  polarisation-torque ceiling, a thread's torsion, the Thomson–Purcell kink field (Zhang's
  `E_θ/e_r`), the carriers' drift (the theory's own bound), the gravity of field energy, the
  electric pressure, and the five-constant identity.
- `plasmoid::discharge` (new): Paschen breakdown, the mean free path, Peek's corona onset.
- `physica/examples/gravity_patent.rs`, and the 3D model [`web/gravity-patent.html`](../web/gravity-patent.html).
- Two statements in the batch-2 review corrected from the primary sources (Setup 1's polarity,
  Setup 2's magnetic field).
- Still open: the constant `k` in `∇×A = B/k` (perhaps in the "mass–charge conversion constant"
  paper), the seven-cable ratio of the weighing rig (×7.1 against ×4.4 from a cage model; it
  needs the real layout in 3D), and the Michelson and 80 V rigs, which have no published details.

Not yet published: the new modules wait for the next release of `gravitics`, `plasmoid` and
`physica`.

## Open

- **Retired crates:** `satgrav`, `gravity-field` and `gravity-claims` are deleted; `gwsens` was
  yanked by the user on 2026-10-24 (deletion was blocked because `gravitics` 0.1.0 and 0.2.0
  depend on it).
- **Decision for the user:** whether to download IGRF-14's coefficient file and EGM96 for the
  IGRF and geoid examples.
- **On hold at the user's request:** `magnetics::{devices, astro, frontier}`. The specification
  is in [`topics/magnetics-and-mhd/report.md`](topics/magnetics-and-mhd/report.md), lines 181–211
  (formulas, test vectors, labels) and 275–278 (worked examples), with the research behind it in
  that topic's `notes/` and its check scripts.

## Later

- `geo5d`: the relativity and PPN calculators (a proposed session).
- The "Verify before coding" items in the gravity-anomalies and magnetics reports.
- W1 scoring when COST-G's December 2026 fields appear (about March 2027).
