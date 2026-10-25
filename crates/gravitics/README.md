# gravitics

*Gravitics* is what the 1950s called the search for gravity control. This crate computes that
whole parallel track of physics, and gravity's open anomalies, three ways:
- **The proposal's own model,** from its authors' formulas and inputs, reproducing their published
  numbers: superconductor gravity, high-frequency gravitational-wave generators and detectors,
  propellantless thrusters, warp drives, the polarizable vacuum, MOND and Weber-like gravity.
- **The standard prediction,** from general relativity or standard physics, for the same setup.
- **What experiments measured.**

It also carries the general-relativity side that light near strong gravity needs: exact light
deflection at any field strength, deflection by a body's whole multipole field, and black-hole
shadows.

Use it to explore a proposal, to design an experiment that could test one, or to see how large
an effect must be before an instrument can see it.

![Each proposal's own numbers against general relativity and experiment](https://raw.githubusercontent.com/FusionVault/physica/main/figures/claims_gap.svg)

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

## Example: an anomaly against the test that decides it

```rust
use gravitics::anomalies::{LLR_ABSORPTION_BOUND, MAJORANA_ABSORPTION, column_density, flyby_anderson, shielding_loss};

// Majorana's gravitational shielding would hide 4.6% of the Sun behind the Earth at eclipse;
// lunar laser ranging allows 7×10⁻¹².
let earth = column_density(6_371e3, 5514.0);
assert!(shielding_loss(MAJORANA_ABSORPTION, earth) > 0.04);
assert!(shielding_loss(LLR_ABSORPTION_BOUND, earth) < 1e-11);

// Anderson's flyby formula reproduces NEAR's 13.3 mm/s, and predicted +6 mm/s for Juno,
// which showed 0 ± 0.8 mm/s.
let near = flyby_anderson(6851.0, (-20.76f64).to_radians(), (-71.96f64).to_radians());
assert!((near * 1e3 - 13.28).abs() < 0.01);
```

## The proposals

| Proposal | Its own model or report | General relativity or standard physics | What experiments found |
|---|---|---|---|
| Tajmar: frame dragging near a spinning superconductor | k = 2.6×10⁻⁵ | 7.2×10⁻²⁷ | Tajmar's own later limit is 800× below the first signal; traced to vibration |
| Podkletnov: weight change above a spinning superconductor | up to −2.1% | a weight increase of at most 3×10⁻⁹ | no change above 2×10⁻⁸ (NASA Marshall) |
| Li–Baker: a GHz gravitational-wave detector | a signal in 3×10³ s | one signal photon in 3×10¹⁸ s, against the reference beam's shot noise (JASON's argument) | none built |
| Pais: 10²⁴ W/m² from a charged disc | σ = 5×10⁴ C/m² | a surface field 5.6×10⁶ times the breakdown limit his patent states | Navy test (2016–19): effect not shown |
| Zhang: a varying field makes gravity (US 12,417,870) | a ring weighs ±1 g at 55 kV: 2.6×10⁻² m/s²; his vortex `∇×A = B/f` with his three constants spins things up at 10⁻²¹ rad/s², or at 10²⁹ rad/s² in Earth's own field; his G is 6.685×10⁻¹¹ | the live wire's polarisation pull, 0.1–0.7 g toward whichever side is live; a charged ball on a thread turns and swings by tens of degrees and a quarter millimetre; general relativity: 10⁻²⁷ m/s² | his G is 34σ above every measurement; the rigs not independently replicated |
| EmDrive | 1.2 mN/kW (Eagleworks) | a photon rocket gives 3.3 µN/kW | nothing above 10–30 nN (Dresden) |
| Mach-effect thruster | "a couple of µN" | the authors' own formula in SI: about 0.1 µN; general relativity suppresses the Machian term by 10⁻³⁶ | vibration artefacts (Dresden) |
| Lifters ("electrogravitics") | lift | ion wind, `F = I·d/μ` | no force in vacuum |
| Warp drive (Alcubierre) | travel faster than light | negative energy of a galaxy or more, in walls quantum field theory allows | — |
| Polarizable vacuum (Puthoff) | gravity as a vacuum index | agrees with general relativity at first order, differs at second | 0.73 µas apart for light grazing the Sun |
| Forward's antigravity | an induced pull from accelerated mass | real in general relativity, with coupling G/c² = 7.4×10⁻²⁸ m/kg | — |

## The anomalies

| Anomaly | The claim | The test | Status |
|---|---|---|---|
| The value of G | 16 measurements span 551 ppm | Birge ratio 3.6 against their errors | open: a measurement problem |
| Earth flybys | Anderson's formula fits six | Juno: +6 mm/s predicted, 0 ± 0.8 observed | failed out of sample |
| Gravitational shielding | Majorana's 6.7×10⁻¹² cm²/g | lunar laser ranging: below 10⁻²² m²/kg | refuted |
| Pioneer | a constant 8.7×10⁻¹⁰ m/s² | it decayed with the plutonium's heat | resolved: thermal recoil |
| MOND | one acceleration scale, g† = 1.2×10⁻¹⁰ m/s² | wide binaries | contested |
| Short-range forces | new physics below a millimetre | Eöt-Wash: no gravitational-strength force beyond 38.6 µm | open below 38.6 µm |
| Diósi–Penrose collapse | spontaneous X-rays | none seen underground | parameter-free version refuted |

## What is inside

| Module | What it computes |
|---|---|
| `electrogravity` | The pulsed-field gravity claims against the forces in the same rig: dielectrophoresis (dilute and on a sphere), the Pauthenier charge, the polarisation-torque ceiling, spheroid depolarisation and the torque on a body not quite round, the image force, a thread's torsion, field angular momentum, the Thomson–Purcell kink field that Zhang's theory relabels, the carriers' drift (the theory's own bound), the gravity of field energy, Zhang's constants (`k_mq`, `μ_g`, his vortex constant and spin-up), and the mass-independence discriminator |
| `electrodes` | The field of high-voltage rigs: coil turns as charged rings and wires as straight segments, with potential and field in closed form (elliptic integrals, logarithms) and the charge on each piece for given potentials |
| `gem` | Gravitoelectromagnetism in three conventions with the spin-2 factor, the gravitomagnetic dipole as a `Field`, the extremal charge-to-mass ratio, and Blackett's magnetism-from-spin |
| `anomalies` | CODATA 2022's sixteen G measurements with their correlations, their generalised-least-squares mean and Birge ratio; Anderson's flyby formula; column densities and shielding against Majorana and lunar laser ranging; Pioneer's thermal decay |
| `modified` | The radial acceleration relation ([McGaugh et al. 2016](https://arxiv.org/abs/1609.05917)), Tully–Fisher velocity and MOND radius; Yukawa forces with the current limits ([Lee et al. 2020](https://arxiv.org/abs/2002.11761)), the dark-energy length and boson mass; Einstein–Cartan density; Weber-like perihelion and mass shift |
| `quantum` | Bose–Marletto–Vedral entanglement phases ([Bose et al. 2017](https://arxiv.org/abs/1707.06050)); the Diósi–Penrose X-ray rate ([Donadi et al. 2021](https://arxiv.org/abs/2111.13490)) |
| `hfgw` | High-frequency gravitational waves: the Li–Baker detector (flux conventions, graviton rate, the interference signal, the shot-noise time) and the Pais patents (radiated flux, surface field, magnetic conversion as he uses it and in SI) |
| `superconductors` | The London moment; frame dragging of a ring in general relativity against Tajmar's report; the Cooper-pair density and gravitational London depth; a perfect shield's shadow; Podkletnov's pendulum and Marx bank; the measured bounds |
| `propulsion` | Ion wind and ion mobility, the vacuum (ballistic) limit, the photon rocket, Shawyer's EmDrive formula, buoyancy on a balance, Woodward's Mach-effect thrust and general relativity's bound on it ([Rodal](https://arxiv.org/abs/2512.22687)) |
| `spacetime` | Alcubierre warp-bubble energy ([Pfenning & Ford 1997](https://arxiv.org/abs/gr-qc/9702026)), the quantum inequality, Casimir plates and spheres, wormhole tensions, light bending in Puthoff's polarizable vacuum against general relativity, to second order and exactly at any field strength (Darwin's integral, checked against [Bozza 2002](https://arxiv.org/abs/gr-qc/0208075) and [Iyer & Petters 2007](https://arxiv.org/abs/gr-qc/0611086)), Forward's induced acceleration, and horizon temperatures (Unruh, Hawking, Gibbons–Hawking) with black-hole lifetimes and zero-point energy; deflection by a body's whole gravity field from `geopotential` coefficients (the Born integral, with J₂ adding 239 µas at Jupiter's limb); black-hole shadows: the Schwarzschild diameter, the EHT's deviation δ, and the Kerr shadow seen from the pole ([EHT Sgr A* Paper VI](https://arxiv.org/abs/2311.09484); Teo 2003) |
| `strong_field` | The Schwinger limit: the critical field, the pair-production and vacuum-decay rates, laser a₀ and the quantum parameter χ (with explicit intensity conventions), Ritus–Narozhny, Breit–Wheeler suppression, vacuum birefringence, Euler–Heisenberg, supercritical atoms, and the graphene analogue. The Schwinger exponent is exactly a Boltzmann factor at twice the Unruh temperature, and a test pins that identity |

![The energy a warp bubble needs](https://raw.githubusercontent.com/FusionVault/physica/main/figures/warp_energy.svg)

*A warp bubble's negative energy grows as its wall thins. Quantum field theory allows walls only
about 100 Planck lengths thick, where a 100 m bubble needs billions of times the ordinary matter
in the observable universe.*

![Thrust per watt](https://raw.githubusercontent.com/FusionVault/physica/main/figures/thrust_per_watt.svg)

*Thrust per watt: the propellantless proposals, a photon rocket, ion wind, and the best
measurement.*

The research behind every number, with its sources and the status of each proposal, is in the
reports on [fringe physics and the DIRDs](https://github.com/FusionVault/physica/blob/main/research/topics/fringe-physics-and-dirds/report.md),
[the Coulomb barrier and the Schwinger limit](https://github.com/FusionVault/physica/blob/main/research/topics/coulomb-and-schwinger/report.md)
and [where gravity physics is still open](https://github.com/FusionVault/physica/blob/main/research/topics/gravity-anomalies/report.md).

## Applications

- **Evaluating a proposal.** Compute a claim in its authors' own terms, the standard prediction
  for the same apparatus, and the measurement, as agency reviews have done; the gap is the
  answer. `examples/electrogravity_claim.rs` does this for the pulsed-field gravity claims and
  names the experiment that would settle them; `physica`'s `gravity_patent` example rebuilds the
  rig of US 12,417,870 with `electrodes` and simulates the ball under every force and torque.
- **Designing an experiment.** `anomalies::tension` and `magnets`-style systematics show what a
  precision gravity measurement must control; `modified` carries the current short-range bounds.
- **Planning a strong-field campaign.** `strong_field` gives the Schwinger field, the quantum
  parameter a beam reaches, and how far the best lasers remain from vacuum breakdown.
- **Teaching the analogy honestly.** `gem` gives gravitoelectromagnetism in three conventions
  that agree, with the spin-2 factor where it belongs.

## How it is checked

Every test reproduces a published number from the paper's own inputs, mostly to better than 1%:
- **Proposals:** the Li–Baker graviton flux and Table I signal; the Pais disc and conversion
  numbers; Tajmar's predicted field for niobium; Podkletnov's pendulum tables; Woodward's printed
  1.87×10²¹ prefactor; Shawyer's thrust at the Dresden and Eagleworks settings; the Army lifter
  thrust.
- **Spacetime:** Pfenning & Ford's warp energy; the quantum inequality; the Casimir pressure at
  1 µm; the Sun's 1.7512″ light bending; Unruh, Hawking and Gibbons–Hawking temperatures.
- **Strong fields:** the Schwinger field and pair rates, CoReLS's a₀, LUXE's χ, PVLAS's QED
  birefringence, and Dirac binding for lead and uranium.
- **Anomalies:** CODATA 2022's G (6.674 30, Birge ratio 3.61, χ² = 195.7 for 15 degrees of
  freedom) rebuilt from its inputs; Anderson's six flyby predictions; Majorana and LLR shielding.
- **Alternatives:** the RAR's published shape and deep-MOND limit; Eöt-Wash's 5.1 meV boson;
  Ruggiero & Tartaglia's Einstein–Cartan densities; the Weber perihelion ratios against a direct
  orbit integration.
- **Quantum gravity:** Bose et al.'s phases; Donadi et al.'s germanium rate.

The formulas were independently re-derived from their sources and every test value recomputed;
the Casimir formulas are also checked against each other.

One of eight crates in [physica](https://github.com/FusionVault/physica):
[`physica`](https://crates.io/crates/physica) (all of them, plus cross-domain analogies) ·
[`metrology`](https://crates.io/crates/metrology) (constants) ·
[`geopotential`](https://crates.io/crates/geopotential) (the gravity field) ·
[`gravimetry`](https://crates.io/crates/gravimetry) (missions and GRACE) ·
[`gravitational`](https://crates.io/crates/gravitational) (waves) ·
[`gravitics`](https://crates.io/crates/gravitics) (frontier gravity) ·
[`magnetics`](https://crates.io/crates/magnetics) (magnetism and MHD) ·
[`plasmoid`](https://crates.io/crates/plasmoid) (fusion plasma).
