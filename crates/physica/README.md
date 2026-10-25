# physica

The whole physica family in one dependency, plus the equations its domains share.

Seven crates cover Earth's gravity field, satellite gravimetry, gravitational waves, frontier
gravity proposals, fusion plasma, magnetism and the constants behind them. They are published
separately so a project can take only what it needs; this crate re-exports all of them and adds
the one thing none of them can hold alone — the identities that join them, each written as a
test.

## Install

```bash
cargo add physica
```

Every crate is re-exported under its own name, so `physica::magnetics::main_field::wmm2025()` and
`magnetics::main_field::wmm2025()` are the same function.

## Example: three fields of the Earth through one interface

Gravity, the main magnetic field and the gravitomagnetic field of Earth's spin are all potential
fields, so one trait evaluates all three.

```rust
use physica::geopotential::{Coefficients, Field};
use physica::gravitics::gem::GravitomagneticDipole;
use physica::magnetics::main_field::wmm2025;

let mut gravity = Coefficients::zeros(3.986_004_415e14, 6_378_136.3, 2);
gravity.set(0, 0, 1.0, 0.0);
gravity.set(2, 0, -4.841_694_95e-4, 0.0);
let magnetism = wmm2025().at(2025.0);
let spin = GravitomagneticDipole { j: 5.86e33 };

let fields: [&dyn Field; 3] = [&gravity, &magnetism, &spin];
let r = 6_371_200.0 + 450e3;
let values: Vec<f64> = fields.iter().map(|f| f.magnitude(r, 0.5, 1.0)).collect();
assert!((8.0..8.8).contains(&values[0]));        // m/s²
assert!((20_000.0..60_000.0).contains(&values[1])); // nT
assert!(values[2] < 1e-5);                        // m/s², and far smaller
```

## Example: the same equation in two domains

```rust
use physica::metrology::{C, G, MU_0};
use physica::plasmoid::fusion::spinning_shell_field;

// Inside a spinning charged shell the magnetic field is (2/3)(μ₀/4π)Qω/R; inside a spinning
// massive shell, frames are dragged at 4GMω/(3c²R). Swap μ₀/4π for G/c² and charge for mass, and
// gravity is exactly twice as strong — because the graviton has spin 2.
let (m, r, omega) = (1e4, 1.0, 300.0);
let analogue = spinning_shell_field(m, r, omega) * (G / (C * C)) / (MU_0 / (4.0 * core::f64::consts::PI));
let frame_dragging = 4.0 * G * m * omega / (3.0 * C * C * r);
assert!((frame_dragging / analogue - 2.0).abs() < 1e-12);
```

## What is inside

| Re-export | What it covers |
|---|---|
| [`metrology`](https://crates.io/crates/metrology) | The physical constants every crate shares |
| [`geopotential`](https://crates.io/crates/geopotential) | Earth's gravity field from spherical harmonics, the `.gfc` reader, the shared `Field` trait |
| [`gravimetry`](https://crates.io/crates/gravimetry) | Gravity missions and quantum sensors, GRACE processing, earthquake gravity, local forward models, atmospheric waves |
| [`gravitational`](https://crates.io/crates/gravitational) | Gravitational waves from nanohertz to terahertz: detectors, inspirals, energy, propagation, gravitons, backgrounds |
| [`gravitics`](https://crates.io/crates/gravitics) | Frontier gravity computed against general relativity and experiment, the Schwinger limit, gravity's anomalies, gravitoelectromagnetism |
| [`magnetics`](https://crates.io/crates/magnetics) | Earth's magnetic field, magnetohydrodynamics, magnets, space weather |
| [`plasmoid`](https://crates.io/crates/plasmoid) | Fusion reactivities and Lawson, pinches, plasmoids, aneutronic power, tokamak limits |
| `analogies` | The cross-domain identities, as doc tests |

## The analogies

Each is an exact identity a test asserts, not a loose resemblance.

| Analogy | The shared equation |
|---|---|
| One spherical-harmonic engine | gravity `V` with `GM, R, C̄` ↔ the magnetic field with `a², a, g/√(2n+1)` |
| Spinning dipoles | both are twice as strong at the pole as at the equator |
| The spin-2 factor | frame dragging inside a spinning shell is exactly 2× its electromagnetic analogue |
| One charge-to-mass constant | `√(4πε₀G)` behind the extremal black hole, the weak-gravity scale, and the drift at which two protons' magnetism matches their gravity — whose product with the proton ratio is exactly `c` |
| Pinch and filament | Bennett's `μ₀I²/4π` ↔ Ostriker's `Gλ²` |
| Schwinger and Unruh | pair-creation suppression `e^(−πE_S/E)` is a Boltzmann factor at the Unruh temperature |
| Gravity's vacuum constants | `1/(4πG)` and `4πG/c²` give `1/√(εμ) = c` |
| One plasma, three crates | the plasma frequency and Alfvén speed written in `gravitational`, `magnetics` and `plasmoid` agree: exactly, and to the 2.7×10⁻⁴ the electrons' mass makes |

## Applications

Taking the whole family at once suits work that crosses domains: evaluating Earth's gravity and
magnetic fields along one orbit, checking a frontier proposal against both the relativity and the
plasma physics it invokes, or teaching where one equation reappears under another name. For a
single domain, depend on that crate alone.

## How it is checked

Every analogy above is a doc test, so `cargo test` asserts each identity to 10⁻¹² or its stated
tolerance. The `examples/` directory evaluates three fields of the Earth through one trait and
reproduces NOAA's published WMM2025 value at 100 km altitude. `gravity_patent` rebuilds the rig
of US 12,417,870 (a ball on a thread between high-voltage coils, claimed to move by a
gravitational field) and simulates it: the field of its coils and leads, every force and torque
on the ball from five crates (charge, polarisation, gas, magnetism, the field's own gravity, and
the claim's own model with its constants), integrated at 1 ms. It is a calibration case for the
physics, so it lives in the examples, not the library. A 3D model of the same rig, computed
with the same equations, is in the repository's
[`web/`](https://github.com/FusionVault/physica/tree/main/web) directory.

One of eight crates in [physica](https://github.com/FusionVault/physica):
[`physica`](https://crates.io/crates/physica) (all of them, plus cross-domain analogies) ·
[`metrology`](https://crates.io/crates/metrology) (constants) ·
[`geopotential`](https://crates.io/crates/geopotential) (the gravity field) ·
[`gravimetry`](https://crates.io/crates/gravimetry) (missions and GRACE) ·
[`gravitational`](https://crates.io/crates/gravitational) (waves) ·
[`gravitics`](https://crates.io/crates/gravitics) (frontier gravity) ·
[`magnetics`](https://crates.io/crates/magnetics) (magnetism and MHD) ·
[`plasmoid`](https://crates.io/crates/plasmoid) (fusion plasma).
