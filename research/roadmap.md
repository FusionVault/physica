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

## Phase C: implement what the review produces

Write each reviewed item into its crate and module, with a published test value or the claim's
own numbers as the test, plus examples where one reproduces a published result.

## Phase D: documentation and release (roadmap steps 4–5), once

4. **Examples and documentation:** the remaining examples from the report's table, an
   Applications section in every README, the `magnetics` and `physica` READMEs, and the root
   README.
5. **Release:** workspace tests and clippy, commit, push, publish `magnetics` and `physica`, and
   new releases of `metrology`, `geopotential`, `gravimetry`, `gravitics` and `plasmoid`.

Open decisions for the user: retire `gwsens` by yanking it (or delete `gravitics` and
republish), and whether to download IGRF-14's coefficient file and EGM96 for the IGRF and geoid
examples.

## Later

- `geo5d`: the relativity and PPN calculators (a proposed session).
- The "Verify before coding" items in the gravity-anomalies and magnetics reports.
- W1 scoring when COST-G's December 2026 fields appear (about March 2027).
