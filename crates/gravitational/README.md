# gravitational

Gravitational waves from nanohertz to terahertz, in one place: what they carry, how they travel,
how sensitive every kind of detector is, and how far a detector can see. Strain conventions,
noise curves, inspirals and binary pulsars, memory, the graviton's mass and number, stochastic
backgrounds, high-frequency limits, and the gravity of lab sources, each tested against published
numbers. Pure Rust, SI units.

## Install

```bash
cargo add gravitational metrology
```

## Example: how far can Advanced LIGO see two neutron stars?

```rust
use gravitational::{curves, inspiral};
use metrology::{MPC, SOLAR_MASS};

let (m1, m2) = (1.4 * SOLAR_MASS, 1.4 * SOLAR_MASS);
let horizon = inspiral::horizon(m1, m2, 8.0, curves::aligo, 20.0);
let range = horizon / inspiral::HORIZON_TO_RANGE;
println!("horizon {:.0} Mpc, average range {:.0} Mpc", horizon / MPC, range / MPC);
```

## Example: the strongest and the weakest

The brightest events radiate about a thousandth of the Planck luminosity `c⁵/G`; the oldest
confirmed prediction is the slow decay of a binary pulsar's orbit. A detectable wave holds about
10³⁶ gravitons per wavelength-sized cell, yet matter is so transparent that water absorbs it
over 10¹¹ Gpc.

```rust
use gravitational::{energy, propagation, quanta};
use metrology::{C, MPC, SOLAR_MASS};

assert!((3.5e49 / energy::PLANCK_LUMINOSITY - 9.6e-4).abs() < 1e-5); // GW150914 at its peak

let pbdot = energy::peters_pbdot(1.438 * SOLAR_MASS, 1.390 * SOLAR_MASS, 0.322997448918 * 86_400.0, 0.6171340);
println!("Hulse–Taylor: {pbdot:.3e} s/s"); // −2.402e-12; measured −2.40263(5)e-12

println!("{:.1e} gravitons", quanta::graviton_count(1e-21, 150.0)); // 3.9e36
let water = propagation::viscous_damping(1e-3, 1e3, 1e-12).unwrap();
println!("absorbed over {:.1e} Gpc of water", C / water / (1e3 * MPC)); // 2.6e11
```

## Example: the ceiling on high-frequency waves

Any gravitational-wave background from the early Universe counts as extra radiation, so Big
Bang nucleosynthesis and the CMB cap it. No detector built or proposed above 10 kHz reaches the
cap.

```rust
use gravitational::{hf, strain::{Convention, hc_from_omega}};

let ceiling = hc_from_omega(hf::CEILING_BBN_CMB, 1e9, 1.0, Convention::OneSided);
println!("at 1 GHz no background can exceed h_c = {ceiling:.1e}"); // 1.3e-30

for limit in hf::REGISTRY {
    println!("{:<45} {:.1e} {} ({:?})", limit.name, limit.value, limit.quantity, limit.status);
}
```

## Example: a gravity sensor as a wave detector

A wave of strain h stretches a baseline L, so two free masses feel a tidal acceleration ½ω²hL.
Gradiometers, atom interferometers and satellite pairs all measure that; the baseline decides
how well.

```rust
use gravitational::detectors::{AtomGwDetector, K_SR87, gradiometer_strain};

// AION-km: two strontium clouds 2 km apart (Badurina et al. 2020).
let aion = AtomGwDetector { baseline: 2000.0, t: 5.0, lmt: 40_000.0, phase_noise: 3e-6, k: K_SR87, diamonds: 1.0 };
println!("{:.1e} per √Hz at 0.1 Hz", aion.peak_strain()); // 2.1e-21

// A GOCE-class gradiometer (10 mE/√Hz) read the same way: 5e-11 per √Hz at 0.1 Hz.
assert!(gradiometer_strain(1e-11, 0.1) > 1e10 * aion.peak_strain());
```

## What is inside

| Module | What it computes |
|---|---|
| `strain` | Ω_GW ↔ S_h ↔ h_c, with the one-sided and two-sided conventions made explicit ([Moore, Cole & Berry 2015](https://arxiv.org/abs/1408.0740)) |
| `cosmology` | Flat ΛCDM comoving, angular-diameter and time-delay distances; the six H0LiCOW lenses give back H₀ = 73.3 ([Wong et al. 2020](https://arxiv.org/abs/1907.04869)) |
| `curves` | Advanced LIGO and Einstein Telescope fits ([Sathyaprakash & Schutz 2009](https://arxiv.org/abs/0903.0338)); LISA with galactic confusion ([Robson, Cornish & Liu 2019](https://arxiv.org/abs/1803.01944)); tabulated curves |
| `inspiral` | Chirp mass, time to merger, ISCO, the energy spectrum, matched-filter SNR, horizon and range ([Chen et al. 2021](https://arxiv.org/abs/1709.08079)); the Kerr ISCO with its energy and angular momentum, and the equatorial photon orbit (Bardeen, Press & Teukolsky 1972) |
| `energy` | The Planck luminosity, strain from luminosity, the circular-binary flux and the multipole flux coefficients ([Blanchet](https://arxiv.org/abs/1310.1528)), Peters' eccentricity enhancement and binary-pulsar decay, and the Shapiro range, shape and delay ([Kramer et al. 2021](https://arxiv.org/abs/2112.06795)) |
| `ringdown` | Kerr quasinormal modes for `l = 2, 3, 4`, every `m`, three overtones, from the fits of [Berti, Cardoso & Will 2006](https://arxiv.org/abs/gr-qc/0512160), checked against their exact tables; ringing frequency and damping time |
| `remnant` | The final spin ([Hofmann, Barausse & Rezzolla 2016](https://arxiv.org/abs/1605.01938)), the radiated energy ([Barausse, Morozova & Rezzolla 2012](https://arxiv.org/abs/1206.3803)) and the horizon area of a merger's black hole, for aligned spins |
| `tides` | Neutron-star tidal deformability `Λ` and the binary `Λ̃` the waveform measures, with GW170817's bounds ([LVC 2019](https://doi.org/10.1103/PhysRevX.9.011001)) |
| `memory` | The leading-order memory of a binary ([Favata 2010](https://arxiv.org/abs/1003.3486)) and its ramp in pulsar timing |
| `propagation` | Speed bounds from GW170817, graviton Compton wavelength, dispersion and delay with the current mass bounds, point-lens lensing in geometrical and wave optics ([Takahashi & Nakamura 2003](https://arxiv.org/abs/astro-ph/0305055)), viscous absorption with its validity check ([Baym et al. 2017](https://arxiv.org/abs/1707.05192)) |
| `quanta` | Graviton counts and occupation numbers, stimulated and spontaneous single-graviton rates and the optimal bar mass ([Tobar et al. 2024](https://arxiv.org/abs/2308.15440)), and a bar's absorption cross-section ([Aguiar 2011](https://arxiv.org/abs/1009.1138)) |
| `sgwb` | Power-law-integrated sensitivity curves ([Thrane & Romano 2013](https://arxiv.org/abs/1310.5300)), Hellings–Downs, pulsar-timing backgrounds ([NANOGrav 2023](https://arxiv.org/abs/2306.16213)), and the energy scale and background of inflation |
| `hf` | The cosmological ceiling, magnetic conversion of waves to photons ([Ejlli et al. 2019](https://arxiv.org/abs/1908.00232)), the plasma cut-off, cavity noise, source lines, and a registry of limits with their status ([Aggarwal et al. 2025](https://arxiv.org/abs/2501.11723)) |
| `detectors` | Gravity sensors as wave detectors: atom-interferometer strain response ([Graham et al. 2016](https://arxiv.org/abs/1606.01860), Eq. 9), and the strain noise of gradiometers, optical clocks ([Kolkowitz et al. 2016](https://arxiv.org/abs/1606.01859)) and satellite pairs |
| `ultralight` | Ultralight dark matter as an oscillating signal: the Compton frequency `mc²/h`, the mass window of the mid-band between LISA and the ground detectors (MAGIS-100, AION), and the coherence length |
| `laboratory` | A beam's gravitational near field ([Rätzel et al. 2016](https://arxiv.org/abs/1511.01023)), waves from spinning bodies, wave energy flux, and the Newtonian calibrators of LIGO ([Ross et al. 2021](https://arxiv.org/abs/2107.00141)) and Virgo ([Estevez et al. 2018](https://arxiv.org/abs/1806.06572)) |

## Applications

- **Detector design.** Horizon and range decide a detector's event rate, which is how
  next-generation designs are argued; `examples/detector_range.rs` works it through.
- **Testing general relativity.** `examples/binary_pulsar.rs` reproduces the Hulse–Taylor orbital
  decay from published orbital elements.
- **Bounding a background.** `sgwb` and `hf` convert between the strain, spectral density and
  energy-density conventions papers quote, and `hf::REGISTRY` carries the published limits with
  their status.
- **Sizing a laboratory source.** `laboratory::quadrupole_strain` says what a lab source radiates,
  and `examples/lab_hfgw.rs` sets it against the best detectors.
- **Multi-messenger timing.** `propagation` gives the speed bounds, dispersion and lensing delays.

## How it is checked

Every formula is tested against numbers published in, or recomputed from, the papers above:
- **Detectors:** the LISA curve at five frequencies; the published peak strains of AION-10,
  AION-100, AION-km and AICE; Kolkowitz et al.'s clock strain.
- **Binaries:** the binary-neutron-star chirp mass, time to merger and ISCO; the decay of
  B1913+16 (to its quoted mass precision) and of the double pulsar; GW170817's spectral fluence.
- **Energy and quanta:** the Planck luminosity; Blanchet's flux coefficients exactly; Tobar et
  al.'s rates and bar masses; Carney et al.'s occupation number; Aguiar's LSU bar cross-section.
- **Propagation:** GW170817's speed bounds; the graviton Compton wavelength against Bernus et al.;
  lensing magnifications and delays; viscous damping in water and air, and its refusal for a
  mantle.
- **Backgrounds:** NANOGrav's energy density; Hellings–Downs values; inflation's energy scale for
  r = 0.036, 0.034 and 0.01.
- **The lab:** the petawatt beam; the Misner–Thorne–Wheeler spinning rod; LIGO's rotor force
  (19.10 pN by direct summation) and Virgo's x f² = 3.3×10⁻¹⁴ m Hz²; Ejlli et al.'s photon-rate
  prefactor; an ADMX-like cavity.

All formulas are written from the papers, not translated from other codes, so the crate stays
MIT-licensed. The research behind the newer modules is in the report on
[where gravity physics is still open](https://github.com/FusionVault/physica/blob/main/research/topics/gravity-anomalies/report.md).

One of eight crates in [physica](https://github.com/FusionVault/physica):
[`physica`](https://crates.io/crates/physica) (all of them, plus cross-domain analogies) ·
[`metrology`](https://crates.io/crates/metrology) (constants) ·
[`geopotential`](https://crates.io/crates/geopotential) (the gravity field) ·
[`gravimetry`](https://crates.io/crates/gravimetry) (missions and GRACE) ·
[`gravitational`](https://crates.io/crates/gravitational) (waves) ·
[`gravitics`](https://crates.io/crates/gravitics) (frontier gravity) ·
[`magnetics`](https://crates.io/crates/magnetics) (magnetism and MHD) ·
[`plasmoid`](https://crates.io/crates/plasmoid) (fusion plasma).
