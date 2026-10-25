# gravimetry

Satellite gravimetry from instrument to map. Model the satellites, from today's GRACE-FO to
missions with quantum (cold-atom) sensors, and simulate what gravity field they recover. Then
turn the monthly fields of GRACE and GRACE-FO into maps of water storage, and model the gravity
change an earthquake causes. Built on [`geopotential`](https://crates.io/crates/geopotential).

![Water storage anomaly in July 2026](https://raw.githubusercontent.com/FusionVault/physica/main/figures/water_2026-07.svg)

*Water stored on land and in the ocean in July 2026, against the 2004–2009 average.*

## Install

```bash
cargo add gravimetry geopotential
```

## Example: from a monthly field to water height

The standard chain of the GRACE Tellus maps
([Landerer & Swenson 2012](https://doi.org/10.1029/2011WR011453)):
1. subtract a baseline;
2. remove the stripes;
3. smooth over 300 km;
4. convert to centimetres of water.

```rust
use geopotential::{Coefficients, Functional};

let (gm, radius) = (3.986004415e14, 6_378_136.3);

// A month and the baseline it is compared with (here: one changed coefficient).
let mut month = Coefficients::zeros(gm, radius, 60);
month.set(2, 0, 1e-10, 0.0);
let baseline = Coefficients::zeros(gm, radius, 60);

let change = gravimetry::anomaly(&month, &baseline);
let change = gravimetry::destripe(&change);
let change = gravimetry::smooth(&change, &gravimetry::gaussian_weights(300e3, radius, 60));

let love = gravimetry::love_numbers(60).unwrap();
let water = change.grid(&Functional::WaterHeight { love: &love }, radius, &[90f64.to_radians()], &[0.0]);
println!("{:.0} mm of water at the North Pole", water[0] * 1000.0); // 19 mm
```

With real data, read each month with `geopotential::read_gfc` and the low-degree corrections
with `technotes::read_tn14` and `read_tn13`. The
[full chain](https://github.com/FusionVault/physica/blob/main/lab/src/main.rs)
replaces C̄20, C̄30 and degree 1 before subtracting the baseline.

## Example: a quantum accelerometer against GRACE-FO

![Instrument noise of electrostatic and cold-atom accelerometers and laser ranging](https://raw.githubusercontent.com/FusionVault/physica/main/figures/noise_spectra.svg)

*Instrument noise, reproducing Figure 1 of
[Zingerle et al. (2024)](https://doi.org/10.1029/2024EA003630) from its equations.*

```rust
use gravimetry::atom::AtomInterferometer;
use gravimetry::missions::{Accelerometer, Ranging, pair_noise};

// A cold-atom accelerometer like CARIOQA's: 1 s between pulses, 80% contrast,
// 100,000 atoms per shot, one shot every 3.5 s, double diffraction in orbit.
let atoms = AtomInterferometer { t: 1.0, contrast: 0.8, atoms: 1e5, cycle: 3.5, diffraction: 2 };
let white = atoms.acceleration_asd(); // about 2.3e-10 m/s² per √Hz

// Noise of a satellite pair at 0.1 mHz, the longest wavelengths, in m/s² per √Hz.
let quantum = pair_noise(Accelerometer::AtomInterferometer { white }, Ranging::Lri2033, 1e-4);
let today = pair_noise(Accelerometer::SuperStar, Ranging::GraceFo, 1e-4);
println!("quantum {quantum:.1e}, GRACE-FO {today:.1e}"); // 3.2e-10 against 1.0e-9
assert!(quantum < today / 3.0);
```

## Example: the gravity change of an earthquake

![The 2025 Kamchatka earthquake in GRACE-FO against the USGS fault model](https://raw.githubusercontent.com/FusionVault/physica/main/figures/e1_kamchatka.svg)

*Left: the gravity step the July 2025 Kamchatka earthquake left in GRACE-FO. Right: the USGS
fault model through the same filter, computed with this crate.*

```rust
use gravimetry::coseismic::{Fault, Medium};

// A thrust fault: 100 × 50 km, 8 m of slip, its centre 20 km deep.
let fault = Fault {
    lat: 52.5, lon: 160.0, depth: 20e3,
    strike: 217.0, dip: 15.0, rake: 90.0,
    length: 100e3, width: 50e3, slip: 8.0, opening: 0.0,
};

// Uplift and gravity change at the surface straight above it.
let above = fault.surface_change(0.0, 0.0, &Medium::default());
println!("uplift {:.2} m, gravity {:+.0} µGal", above.uplift, above.gravity * 1e8); // 1.20 m, +115 µGal

// What a satellite sees: surface masses on a grid, as coefficients, smoothed like the data.
let lats: Vec<f64> = (0..50).map(|i| (47.0 + 0.2 * i as f64).to_radians()).collect();
let lons: Vec<f64> = (0..70).map(|j| (153.0 + 0.2 * j as f64).to_radians()).collect();
let at_sea = |_lat: f64, _lon: f64| true;
let masses = gravimetry::coseismic::surface_masses(&[fault], &lats, &lons, &Medium::default(), at_sea);
let field = gravimetry::mass_coefficients(&masses, 3.986004415e14, 6_378_136.3, 60);
let field = gravimetry::smooth(&field, &gravimetry::gaussian_weights(300e3, field.radius, 60));
```

For a real earthquake, `coseismic::read_coulomb` reads the `coulomb.inp` file that the
[USGS](https://earthquake.usgs.gov/) publishes with every finite-fault model.

## Example: fill the gaps in a monthly record

The GRACE record has gaps: missing months, and almost a year between the two missions. The
Kalman smoother treats each coefficient as a slowly wandering signal plus noise. It gives the best
estimate and its uncertainty for every month, gaps included.

```rust
use gravimetry::kalman::{Process, smooth};

// Five months, the middle one missing, each observed with an error of 1.
let times = [0.0, 1.0, 2.0, 3.0, 4.0].map(|month| month / 12.0);
let observed = [Some((1.0, 1.0)), Some((2.0, 1.0)), None, Some((2.5, 1.0)), Some((2.0, 1.0))];

// A signal that varies by about 2 over about 6 months.
let result = smooth(&times, &observed, Process { sd: 2.0, tau: 0.5 });
println!("month 3: {:.2} ± {:.2}", result.mean[2], result.sd[2]);
```

`kalman::fit` picks the signal's size and correlation time, and checks the stated errors, by
maximum likelihood over many series at once.

## What is inside

| Item | What it does |
|---|---|
| `missions` | Noise spectra of electrostatic and cold-atom accelerometers and laser ranging, and of a whole satellite pair ([Zingerle et al. 2024](https://doi.org/10.1029/2024EA003630), eqs. 1–8) |
| `atom` | Phase, shot noise, and accelerometer and gradiometer noise of a cold-atom sensor, in orbit (double diffraction) or on the ground (single), its phase sensitivity and the gravitational Aharonov–Bohm rate ([Lévèque et al. 2021](https://arxiv.org/abs/2011.03382); [Stray et al. 2022](https://doi.org/10.1038/s41586-021-04315-3); [CARIOQA](https://arxiv.org/abs/2211.01215); [Trimeche et al. 2019](https://doi.org/10.1088/1361-6382/ab4548)) |
| `gradiometry` | The error of each coefficient a gradiometer mission like GOCE achieves |
| `local` | Gravity and gravity gradient of a buried sphere and horizontal cylinder, the shapes a ground survey fits, and the averaging time an instrument needs to see them |
| `sim` | A closed-loop simulation of a GRACE-type pair: simulate, estimate the field back, and compare its errors with the predicted ones |
| `noise` | Noise with any spectrum, its whitening, and the FFT |
| `technotes` | NASA's TN-14 (C̄20 and C̄30 from laser ranging, [Loomis et al. 2020](https://doi.org/10.1029/2019GL085488)) and TN-13 (degree 1, [Sun et al. 2016](https://doi.org/10.1002/2016JB013073)), from [PO.DAAC](https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-docs/gracefo/open/docs/) |
| `love_numbers` | How much the solid Earth sags under a load, per degree ([Wahr et al. 1998](https://doi.org/10.1029/98JB02844), Table 1) |
| `gaussian_weights`, `smooth` | Averaging over a chosen radius ([Wahr et al. 1998](https://doi.org/10.1029/98JB02844), after Jekeli 1981) |
| `destripe` | Removes the north–south stripes of GRACE noise ([Swenson & Wahr 2006](https://doi.org/10.1029/2005GL025285)), ported from [gravity-toolkit](https://github.com/polargeodesy/gravity-toolkit) |
| `mean`, `anomaly` | Baseline means and differences, across `GM` and radius conventions |
| `mass_coefficients` | Spherical-harmonic coefficients of point masses |
| `transport_weights` | The smallest mass flow that carries a change between months |
| `coseismic` | Uplift and gravity change of faults in an elastic half-space ([Okada 1985](https://doi.org/10.1785/BSSA0750041135); [Okubo 1992](https://doi.org/10.1029/92JB00178)); reads USGS `coulomb.inp` files and corrects for seawater pushed aside by seafloor uplift |
| `kalman` | A Kalman smoother for series with gaps and per-month errors, and a maximum-likelihood fit of its settings ([Rauch et al. 1965](https://doi.org/10.2514/3.3166)) |
| `atmosphere` | Sound, buoyancy and Lamb waves in the air, the mass of a pressure change, and the ocean's forced response to a fast pressure pulse ([Wright et al. 2022](https://doi.org/10.1038/s41586-022-05012-5); [Proudman 1929](https://doi.org/10.1111/j.2044-8312.1929.tb05286.x)) |

![Closed-loop simulation of four mission classes](https://raw.githubusercontent.com/FusionVault/physica/main/figures/sst_closed_loop.svg)

*The simulated error of four mission classes, from today's to a quantum pair, against the
errors of a real month.*

## Applications

- **Drought and groundwater.** GRACE and GRACE-FO mass change feeds weekly drought indicators and
  statutory groundwater reporting; the processing chain here is the one those products use.
- **Ice sheets and sea level.** `examples/ice_and_sea_level.rs` turns Greenland's published mass
  trend into its sea-level contribution.
- **Earthquakes.** `coseismic` models the gravity change of a fault and reads USGS finite-fault
  files; `prompt_lead` gives how far ahead of the P-wave a light-speed gravity signal arrives.
- **Site surveys and monitoring.** `local` gives the gravity and gradient of a buried void or
  dense body, and the averaging time a quantum gravimeter needs to resolve it.
- **Mission design.** `missions`, `atom` and `sim` size the sensors and simulate what field a
  satellite pair recovers, including cold-atom instruments.

## How it is checked

- **Noise models:** they reproduce Figure 1 of Zingerle et al. (2024).
- **Cold-atom sensors:** they reproduce the published CARIOQA, GRICE and Trimeche et al. (2019)
  sensitivities.
- **Gradiometry errors:** the closed form matches brute-force least squares.
- **Simulator:** it recovers a known field exactly without noise; with noise, a χ² test confirms
  that its actual errors match its predicted ones.
- **Technical notes:** real NASA files, excerpted in `tests/data`.
- **Love numbers:** against the published table. **Smoothing weights:** against their
  closed-form limit. **Stripe removal:** on synthetic patterns.
- **Fault model:** Okada's (1985) published check values; mass conservation (a buried fault's
  equivalent mass sums to zero); on the 2025 Kamchatka fault model, its uplift matches the
  USGS's own surface displacement with r = 0.89.
- **Point masses:** their potential matches Newton's law.
- **Kalman smoother:** it matches the exact Gaussian-process solution to 1e-10, its fit
  recovers the settings of simulated data, and on the real record withheld months are predicted
  with calibrated errors.
- **Atmosphere:** the isothermal values at 250 K (sound 317 m/s, cut-off 0.0217 rad/s, buoyancy
  0.0196 rad/s); the Lamb wave is trapped at every wavelength; the Tonga Lamb speed of 318.2 m/s
  is sound at 252 K; a slow pressure pulse gives the inverse barometer.

**Not included yet:** glacial isostatic adjustment, DDK filters, basin averages and leakage
corrections; time-variable gravity in the simulator, so no temporal aliasing.

One of eight crates in [physica](https://github.com/FusionVault/physica):
[`physica`](https://crates.io/crates/physica) (all of them, plus cross-domain analogies) ·
[`metrology`](https://crates.io/crates/metrology) (constants) ·
[`geopotential`](https://crates.io/crates/geopotential) (the gravity field) ·
[`gravimetry`](https://crates.io/crates/gravimetry) (missions and GRACE) ·
[`gravitational`](https://crates.io/crates/gravitational) (waves) ·
[`gravitics`](https://crates.io/crates/gravitics) (frontier gravity) ·
[`magnetics`](https://crates.io/crates/magnetics) (magnetism and MHD) ·
[`plasmoid`](https://crates.io/crates/plasmoid) (fusion plasma).
