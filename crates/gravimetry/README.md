# gravimetry

Turn monthly gravity fields from GRACE and GRACE-FO into maps of water storage, and model the
gravity change an earthquake causes. Built on
[`geopotential`](https://crates.io/crates/geopotential).

![Water storage anomaly in July 2026](https://raw.githubusercontent.com/FusionVault/quantum-gravimeter/main/figures/water_2026-07.svg)

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

With real data, read the month and the technical notes with
[`gravity-field`](https://crates.io/crates/gravity-field). The
[full chain](https://github.com/FusionVault/quantum-gravimeter/blob/main/apps/gravviz/src/main.rs)
also replaces C̄20, C̄30 and degree 1.

## Example: the gravity change of an earthquake

![The 2025 Kamchatka earthquake in GRACE-FO against the USGS fault model](https://raw.githubusercontent.com/FusionVault/quantum-gravimeter/main/figures/e1_kamchatka.svg)

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

![Are COST-G's stated errors right?](https://raw.githubusercontent.com/FusionVault/quantum-gravimeter/main/figures/kalman.svg)

*The smoother fitted to 24 years of monthly fields: how far COST-G's stated errors need scaling
per degree, and how well withheld months are predicted.*

## What is inside

| Item | What it does |
|---|---|
| `love_numbers` | How much the solid Earth sags under a load, per degree ([Wahr et al. 1998](https://doi.org/10.1029/98JB02844), Table 1) |
| `gaussian_weights`, `smooth` | Averaging over a chosen radius ([Wahr et al. 1998](https://doi.org/10.1029/98JB02844), after Jekeli 1981) |
| `destripe` | Removes the north–south stripes of GRACE noise ([Swenson & Wahr 2006](https://doi.org/10.1029/2005GL025285)), ported from [gravity-toolkit](https://github.com/polargeodesy/gravity-toolkit) |
| `mean`, `anomaly` | Baseline means and differences, across `GM` and radius conventions |
| `coseismic` | Uplift and gravity change of faults in an elastic half-space ([Okada 1985](https://doi.org/10.1785/BSSA0750041135); [Okubo 1992](https://doi.org/10.1029/92JB00178)). It also reads USGS `coulomb.inp` files and corrects for seawater pushed aside by seafloor uplift |
| `mass_coefficients` | Spherical-harmonic coefficients of point masses |
| `kalman` | A Kalman smoother for series with gaps and per-month errors, and a maximum-likelihood fit of its settings ([Rauch et al. 1965](https://doi.org/10.2514/3.3166)) |
| `transport_weights` | The smallest mass flow that carries a change, the flow field in the [GRACE Mass Atlas](https://github.com/FusionVault/quantum-gravimeter#run-the-viewer) |

## How it is checked

- **Love numbers:** against the published table.
- **Smoothing weights:** against their closed-form limit.
- **Stripe removal:** on synthetic patterns.
- **Fault model:**
  - Okada's (1985) published check values;
  - mass conservation: a buried fault's equivalent mass sums to zero;
  - on the 2025 Kamchatka fault model, its uplift matches the USGS's own surface displacement
    with r = 0.89.
- **Point masses:** their potential matches Newton's law.
- **Kalman smoother:**
  - it matches the exact Gaussian-process solution to 1e-10;
  - its fit recovers the settings of simulated data;
  - on the real record, withheld months are predicted with calibrated errors.

**Not included yet:** glacial isostatic adjustment, DDK filters, basin averages and leakage
corrections.

Part of [quantum-gravimeter](https://github.com/FusionVault/quantum-gravimeter).
