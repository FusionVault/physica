# satgrav

Model gravity-mapping satellites, from today's GRACE-FO to missions with quantum (cold-atom)
sensors. It covers how noisy the instruments are, what errors the gravity field ends up with, and
a full simulation of a satellite pair.

![Instrument noise of electrostatic and cold-atom accelerometers and laser ranging](https://raw.githubusercontent.com/FusionVault/quantum-gravimeter/main/figures/noise_spectra.svg)

*Instrument noise, reproducing Figure 1 of
[Zingerle et al. (2024)](https://doi.org/10.1029/2024EA003630) from its equations.*

## Install

```bash
cargo add satgrav
```

## Example: a quantum accelerometer against GRACE-FO

```rust
use satgrav::{Accelerometer, AtomInterferometer, Ranging, pair_noise};

// A cold-atom accelerometer like CARIOQA's: 1 s between pulses, 80% contrast,
// 100,000 atoms per shot, one shot every 3.5 s.
let atoms = AtomInterferometer { t: 1.0, contrast: 0.8, atoms: 1e5, cycle: 3.5 };
let white = atoms.acceleration_asd(); // about 2.3e-10 m/s² per √Hz

// Noise of a satellite pair at 0.1 mHz, the longest wavelengths, in m/s² per √Hz.
let quantum = pair_noise(Accelerometer::AtomInterferometer { white }, Ranging::Lri2033, 1e-4);
let today = pair_noise(Accelerometer::SuperStar, Ranging::GraceFo, 1e-4);
println!("quantum {quantum:.1e}, GRACE-FO {today:.1e}"); // 3.2e-10 against 1.0e-9
assert!(quantum < today / 3.0);
```

## What is inside

| Item | What it does |
|---|---|
| `Accelerometer`, `Ranging`, `pair_noise` | Noise spectra of electrostatic and cold-atom accelerometers and laser ranging, and of a whole satellite pair ([Zingerle et al. 2024](https://doi.org/10.1029/2024EA003630), eqs. 1–8) |
| `AtomInterferometer` | Phase, shot noise, and accelerometer and gradiometer noise of a cold-atom sensor in orbit ([Lévèque et al. 2021](https://arxiv.org/abs/2011.03382); [CARIOQA](https://arxiv.org/abs/2211.01215); [Trimeche et al. 2019](https://doi.org/10.1088/1361-6382/ab4548)) |
| `gradiometry_errors` | The error of each coefficient a gradiometer mission like GOCE achieves |
| `noise` | Noise with any spectrum, and its whitening |
| `sim` | A closed-loop simulation of a GRACE-type pair: simulate, estimate the field back, and compare its errors with the predicted ones |

![Closed-loop simulation of four mission classes](https://raw.githubusercontent.com/FusionVault/quantum-gravimeter/main/figures/sst_closed_loop.svg)

*The simulated error of four mission classes, from today's to a quantum pair, against the
errors of a real month.*

## How it is checked

- **Noise models:** they reproduce Figure 1 of Zingerle et al. (2024).
- **Cold-atom sensors:** they reproduce the published CARIOQA, GRICE and Trimeche et al. (2019)
  sensitivities.
- **Gradiometry errors:** the closed form matches brute-force least squares.
- **Simulator:** it recovers a known field exactly without noise. With noise, a χ² test confirms
  that its actual errors match its predicted ones.

**Not modelled yet:**
- time-variable gravity, and so temporal aliasing;
- orbit dynamics and attitude errors;
- multi-axis gradiometers.

Part of [quantum-gravimeter](https://github.com/FusionVault/quantum-gravimeter).
