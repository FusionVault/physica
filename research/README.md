# Research

Everything behind the crates, in one structure:

- **`topics/<topic>/`**: a research report (`report.md`) and the notes it was written from
  (`notes/`). The notes include Python check scripts (`*_checks.py`, standard library only) that
  recompute every number the crates are tested against.
- **`tests/<id>/`**: one folder per test of the data. A registered test has its pre-registration
  (`prereg.md`), any freeze of its pipeline (`freeze.md`) and input checksums; every test that has
  been run has its `result.md`, reported whether it hit or missed.
- **`audits/`**: independent re-derivations of the crates' formulas, and what they changed.
- **`external/`**: research the user did elsewhere, kept verbatim until it is reviewed
  against the crates ([index](external/README.md)).
- **[`roadmap.md`](roadmap.md)**: what is in progress and what comes next.

The figures and analysis outputs are in [`figures/`](../figures), made by the
[`lab`](../lab) crate.

## Topics

| Topic | Report | What it covers |
|---|---|---|
| Satellite gravimetry data | [Pair open GRACE data with simulated quantum sensors](topics/satellite-gravimetry-data/report.md) | The data, the missions, the formulas and the Rust ecosystem |
| Gravimetry frontiers | [Turn 259 gravity months into discoveries](topics/gravimetry-frontiers/report.md) | New tests linking gravity data to climate, earthquakes and the solid Earth |
| Gravitational waves | [Chart gravitational waves honestly, from nanohertz to terahertz](topics/gravitational-waves/report.md) | Detectors, high-frequency waves, atom interferometers, gravity beams, and atmospheric waves seen from orbit |
| Fringe physics and the DIRDs | [Fringe gravity's sums reproduce, its physics never does](topics/fringe-physics-and-dirds/report.md) | Li–Baker, Pais, Podkletnov, Tajmar, the DIA "DIRD" papers and the public record of gravity programmes |
| The Coulomb barrier and the Schwinger limit | [Tunnelling crosses the barrier; the vacuum holds firm](topics/coulomb-and-schwinger/report.md) | Fusion tunnelling, cold-fusion claims, strong-field QED, and the Unruh–Schwinger link to gravity |
| Gravity anomalies | [Replication erases gravity's anomalies, except the statistical ones](topics/gravity-anomalies/report.md) | G, flybys, shielding, Pioneer, MOND, alternative theories, lab gravity and deep gravitational-wave physics |
| Plasmoids and fusion | [Pinches confine; bremsstrahlung still beats boron](topics/plasmoids-and-fusion/report.md) | Z-pinches, the plasma focus, plasmoids and reconnection, aneutronic fusion and direct conversion |
| Magnetics and MHD | [One harmonic engine joins gravity and magnetism](topics/magnetics-and-mhd/report.md) | Geomagnetic main-field models on the gravity engine, MHD, magnets and devices, space weather, gravity–magnetism links and claims, and the practical uses of every crate |
| Cross-domain frontiers | [Shared kernels, not a new crate, unify physica](topics/cross-domain-frontiers/report.md) | 2022–2026 papers and patents in plasma, magnetics, gravity and light: 42 testable additions in six phases, the shared potential-field core, ten defects in the published crates, and patent verdicts |
| Frontier gaps and metric engineering | [Pulsars exclude Puthoff's theory, not his metric](topics/frontier-gaps-and-metric-engineering/report.md) | The Phase 1 and Phase 3 gaps with methods and test values; Puthoff's 2010 metric-engineering paper audited claim by claim; every test of the polarizable-vacuum metric, strong field included; the warp-drive and wormhole literature since 1994; eleven defects in the shipped crates |
| Topological monopoles and transmutation | [Monopole mathematics holds; the talk's laboratory signatures are unverified](topics/topological-monopoles-and-transmutation/report.md) | Greenyer's TeslaTech 2026 talk computed both ways: Dirac and GUT monopoles, catalysis and search bounds, strange radiation, the slide's nine transmutation names, ball lightning and charge clusters, toroidal moments, propulsion, cavitation and oxyhydrogen; the user's audit checked item by item; about 120 calculators across seven issues |
| Toroidal novelty | [Toroidal moments hold one publishable result, about pulses, and two errata](topics/toroidal-novelty/report.md) | Whether toroidal moments, anapoles, monopole duality and their gravitational analogues hold anything new: four lines of literature search, proof and numerical check |

## Tests

| Test | Question | Files | Status |
|---|---|---|---|
| W1 | How much water does land lose in an El Niño, and does it come back? | [prereg](tests/W1/prereg.md), [inputs](tests/W1/inputs.sha256) | Registered 2026-10-24; scored when COST-G's fields from December 2026 on are out (about March 2027) |
| E1 | Can monthly gravity see the 2025 Kamchatka earthquake, and is there more than the slip explains? | [result](tests/E1/result.md) | Detected; matches the USGS fault model, with no excess |
| Q0 | What makes a monthly field noisy? | [result](tests/Q0/result.md) | Solar activity and satellite modes |
| Q2 | Are the stated errors of the monthly fields right? | [result](tests/Q2/result.md) | Too small, by about 2× above degree 10 |
| M1 | Can a blind search of every month and place find earthquakes? | [result](tests/M1/result.md) | Sumatra and Tohoku, with no false alarms |
| T1 | Does GRACE-FO's laser ranging feel the 2022 Tonga eruption's pressure wave? | [prereg](tests/T1/prereg.md), [freeze](tests/T1/freeze.md), [result](tests/T1/result.md) | Not detected: 1.75 ± 1.52 hPa (Z = 1.15, p = 0.12; 95% limit 4.3 hPa, against the expected 2.1–3.9) |

**How registration works.** A forecast is committed and tagged before the data that decides it
exists (W1), or before the analysis is run on it (T1), and the tag is pushed so its date is
public. The registrations and the freeze are kept exactly as tagged, so the paths inside them
refer to the repository at their tag:
[`prereg/W1-2026-10-24`](https://github.com/FusionVault/physica/tree/prereg/W1-2026-10-24),
[`prereg/T1-2026-10-24`](https://github.com/FusionVault/physica/tree/prereg/T1-2026-10-24) and
[`prereg/T1-frozen`](https://github.com/FusionVault/physica/tree/prereg/T1-frozen).

## Audits

- [Formula audit, 2026-10-24](audits/formula-audit-2026-10-24.md): an independent re-derivation of
  every formula in the crates; four errors and several convention problems found and fixed.
