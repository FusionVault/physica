# gwsens

Gravitational-wave sensitivity in one place, from nanohertz to terahertz:
- strain conventions;
- detector noise curves;
- signal-to-noise and how far a detector can see;
- stochastic backgrounds and pulsar timing;
- high-frequency limits;
- lab-gravity reality checks.

Pure Rust, SI units, no dependencies.

## Install

```bash
cargo add gwsens
```

## Example: how far can Advanced LIGO see two neutron stars?

```rust
use gwsens::{MPC, SOLAR_MASS, curves, inspiral};

let (m1, m2) = (1.4 * SOLAR_MASS, 1.4 * SOLAR_MASS);
let horizon = inspiral::horizon(m1, m2, 8.0, curves::aligo, 20.0);
let range = horizon / inspiral::HORIZON_TO_RANGE;
println!("horizon {:.0} Mpc, average range {:.0} Mpc", horizon / MPC, range / MPC);
```

## Example: the ceiling on high-frequency waves

Any gravitational-wave background from the early Universe counts as extra radiation, so Big
Bang nucleosynthesis and the CMB cap it. No detector built or proposed above 10 kHz reaches the
cap.

```rust
use gwsens::{hf, strain::{Convention, hc_from_omega}};

let ceiling = hc_from_omega(hf::CEILING_BBN_CMB, 1e9, 1.0, Convention::OneSided);
println!("at 1 GHz no background can exceed h_c = {ceiling:.1e}"); // 1.3e-30

for limit in hf::REGISTRY {
    println!("{:<45} {:.1e} {} ({:?})", limit.name, limit.value, limit.quantity, limit.status);
}
```

## Example: a light beam's gravity

```rust
use gwsens::lab;

// A petawatt laser pulls a mass 1 mm away with 10⁻¹⁷ m/s², only while it passes.
println!("{:.1e} m/s²", lab::beam_acceleration(1e15, 1e-3));
```

## What is inside

| Module | What it does |
|---|---|
| `strain` | Ω_GW ↔ S_h ↔ h_c, with the one-sided and two-sided conventions made explicit ([Moore, Cole & Berry 2015](https://arxiv.org/abs/1408.0740)) |
| `curves` | Advanced LIGO and Einstein Telescope fits ([Sathyaprakash & Schutz 2009](https://arxiv.org/abs/0903.0338)); LISA with galactic confusion ([Robson, Cornish & Liu 2019](https://arxiv.org/abs/1803.01944)); tabulated curves |
| `inspiral` | Chirp mass, time to merger, ISCO, matched-filter SNR, horizon and range ([Chen et al. 2021](https://arxiv.org/abs/1709.08079)) |
| `sgwb` | Power-law-integrated sensitivity curves ([Thrane & Romano 2013](https://arxiv.org/abs/1310.5300)), Hellings–Downs, pulsar-timing backgrounds ([NANOGrav 2023](https://arxiv.org/abs/2306.16213)) |
| `hf` | The cosmological ceiling, magnetic conversion of waves to photons ([Ejlli et al. 2019](https://arxiv.org/abs/1908.00232)), cavity noise, source lines, and a registry of limits with their status ([Aggarwal et al. 2025](https://arxiv.org/abs/2501.11723)) |
| `lab` | A beam's gravitational near field ([Rätzel et al. 2016](https://arxiv.org/abs/1511.01023)), waves from spinning bodies, wave energy flux |

## How it is checked

Every formula is tested against numbers published in, or recomputed from, the papers above:
- the LISA curve at five frequencies;
- the binary-neutron-star chirp mass, time to merger and ISCO;
- NANOGrav's energy density;
- Hellings–Downs values;
- Ejlli et al.'s photon-rate prefactor;
- an ADMX-like cavity;
- the petawatt beam;
- the Misner–Thorne–Wheeler spinning rod.

All formulas are written from the papers, not translated from other codes, so the crate stays
MIT-licensed.

Part of [quantum-gravimeter](https://github.com/FusionVault/quantum-gravimeter).
