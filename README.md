# quantum-gravimeter

Rust crates for satellite gravimetry: turn open GRACE / GRACE-FO / GOCE / Swarm gravity data into
geoid, gravity and water-storage maps, and simulate how future quantum (cold-atom) sensors would
change the errors.

![July 2026 water storage anomaly from GRACE-FO](figures/water_2026-07.svg)

## The crates

| Crate | What it does |
|---|---|
| [`geopotential`](https://crates.io/crates/geopotential) ([docs](https://docs.rs/geopotential)) | The spherical-harmonic engine:<br>• Legendre functions to degree 2700<br>• gravity potential, gradient and Hessian at a point<br>• grids of geoid, gravity anomaly or water height<br>• degree spectra, tide systems |
| [`gravity-field`](https://crates.io/crates/gravity-field) ([docs](https://docs.rs/gravity-field)) | Reads gravity-field files:<br>• ICGEM `.gfc` models, static or time-variable<br>• GRACE technical notes TN-13 and TN-14 |
| [`gravimetry`](https://crates.io/crates/gravimetry) ([docs](https://docs.rs/gravimetry)) | GRACE processing: Love numbers, baseline means, stripe removal, Gaussian smoothing |
| [`satgrav`](https://crates.io/crates/satgrav) ([docs](https://docs.rs/satgrav)) | Missions and sensors:<br>• accelerometer and laser-ranging noise models<br>• atom-interferometer physics<br>• error propagation<br>• a closed-loop simulation of a GRACE-type satellite pair |
| [`gravviz`](apps/gravviz) | Command-line tool that makes the figures below and the viewer's data (not published) |
| [`viewer`](viewer) | A web page for exploring every month from 2002 to 2026 on a map or globe |

[`geo5d`](https://crates.io/crates/geo5d) (normal gravity, height systems, relativistic geodesy) uses
`geopotential` for its gravity synthesis.

## Quick start

Add the crates you need:

```bash
cargo add geopotential gravity-field gravimetry satgrav
```

**Gravity at a point from a published model.** Download a model such as GOCO06s from the
[ICGEM static models table](https://icgem.gfz.de/tom_longtime), then:

```rust
use gravity_field::{julian_epoch, read_gfc};

let model = read_gfc(&std::fs::read_to_string("GOCO06s.gfc")?)?.model;
let field = model.at(julian_epoch(2026, 7, 1, 0, 0)); // applies the trend and annual terms

let (r, lat, lon) = (6_378_137.0 + 400e3, 45f64.to_radians(), 7f64.to_radians());
let p = field.partials(r, lat, lon); // geocentric radius (m), latitude, longitude (rad)
let [north, east, up] = p.acceleration_neu(r, lat);
println!("V = {:.3} m²/s², g = ({north:.3e}, {east:.3e}, {up:.6}) m/s²", p.v);
```

**One GRACE month as water height.** This is the same chain NASA's GRACE Tellus maps use:
1. Swap in the better C₂₀ from satellite laser ranging.
2. Subtract a baseline mean.
3. Cut to degree 60.
4. Remove the stripes.
5. Smooth over 300 km.
6. Convert to water height.

```rust
use geopotential::Functional;
use gravity_field::{read_gfc, read_tn14};

let read = |path: &str| -> Result<geopotential::Coefficients, Box<dyn std::error::Error>> {
    Ok(read_gfc(&std::fs::read_to_string(path)?)?.model.base)
};
let mut month = read("costg/GSM-2_2026182-2026212_GRFO_COSTG_BF01_0201.gfc")?;
let slr = read_tn14(&std::fs::read_to_string("TN-14_C30_C20_GSFC_SLR.txt")?)?;
let july = slr.iter().find(|row| (2026.49..2026.51).contains(&row.start)).ok_or("no TN-14 month")?;
month.set(2, 0, july.c20, 0.0);

let baseline = vec![read("costg/GSM-2_2004001-2004013_GRAC_COSTG_BF01_0201.gfc")? /* , … every 2004–2009 month */];
let mean = gravimetry::mean(&baseline).ok_or("empty baseline")?;
let anomaly = gravimetry::destripe(&gravimetry::anomaly(&month, &mean).truncated(60));
let anomaly = gravimetry::smooth(&anomaly, &gravimetry::gaussian_weights(300e3, month.radius, 60));

let love = gravimetry::love_numbers(60).ok_or("degree above the table")?;
let lats: Vec<f64> = (-89..=89).map(|d| f64::from(d).to_radians()).collect();
let lons: Vec<f64> = (-180..180).map(|d| f64::from(d).to_radians()).collect();
let metres = anomaly.grid(&Functional::WaterHeight { love: &love }, month.radius, &lats, &lons);
```

[`apps/gravviz/src/main.rs`](apps/gravviz/src/main.rs) has the full version, which also replaces
C₃₀ and adds degree 1 from TN-13.

**A cold-atom accelerometer against GRACE-FO.**

```rust
use satgrav::{Accelerometer, AtomInterferometer, Ranging, pair_noise};

// CARIOQA-like: T = 1 s, contrast 0.8, 1e5 atoms, 3.5 s cycle.
let cai = AtomInterferometer { t: 1.0, contrast: 0.8, atoms: 1e5, cycle: 3.5 };
let white = cai.acceleration_asd(); // ≈ 2.3e-10 m s⁻² Hz⁻¹ᐟ²
let quantum = pair_noise(Accelerometer::AtomInterferometer { white }, Ranging::Lri2033, 1e-4);
let today = pair_noise(Accelerometer::SuperStar, Ranging::GraceFo, 1e-4);
// At 0.1 mHz (the longest wavelengths): 3.2e-10 against 1.0e-9.
```

These three examples are compiled against the published crates.

## Reproduce the figures

You need Rust 1.85 or later. All data below is free; none of it needs an account.

**1. Get the data** (about 140 MB unpacked):

```bash
mkdir -p data && cd data
# Every COST-G monthly GRACE/GRACE-FO field, 2002–2026, as one zip (ICGEM rate-limits single files).
curl -L -o costg.zip https://icgem.gfz.de/getseries/02_COST-G_/Grace-Grace-FO_RL02.1
unzip -q -j costg.zip -d costg
mkdir -p baseline && cp costg/GSM-2_200[4-9]*.gfc baseline/
# GRACE technical notes: C20/C30 from laser ranging (TN-14) and degree 1 (TN-13).
curl -LO https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-docs/gracefo/open/docs/TN-14_C30_C20_GSFC_SLR.txt
curl -LO https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-docs/gracefo/open/docs/TN-13_GEOC_CSR_RL0603.txt
# GOCO06s, a static model (for the gradiometry figure).
curl -L -o GOCO06s.gfc https://icgem.gfz.de/getmodel/gfc/32ec2884630a02670476f752d2a2bf1c395d8c8d6d768090ed95b4f04b0d5863/GOCO06s.gfc
cd ..
```

**2. Make the figures** (each writes an SVG, plus a CSV of the plotted values for the line figures
and the map):

| Command | Output |
|---|---|
| `cargo run --release -p gravviz -- noise figures` | [`noise_spectra.svg`](figures/noise_spectra.svg): instrument noise, Zingerle et al. (2024) Fig. 1, from its equations |
| `cargo run --release -p gravviz -- water-map data/costg/GSM-2_2026182-2026212_GRFO_COSTG_BF01_0201.gfc data/baseline data/TN-14_C30_C20_GSFC_SLR.txt data/TN-13_GEOC_CSR_RL0603.txt figures` | [`water_2026-07.svg`](figures/water_2026-07.svg): the July 2026 water-storage map above |
| `cargo run --release -p gravviz -- triangle data/costg/GSM-2_2026182-2026212_GRFO_COSTG_BF01_0201.gfc figures` | [`error_triangle_2026-07.svg`](figures/error_triangle_2026-07.svg): that month's error per coefficient |
| `cargo run --release -p gravviz -- gradiometry data/GOCO06s.gfc figures` | [`gradiometry_errors.svg`](figures/gradiometry_errors.svg): GOCE-like and cold-atom gradiometers against GOCO06s |
| `cargo run --release -p gravviz -- sst figures data/costg/GSM-2_2026182-2026212_GRFO_COSTG_BF01_0201.gfc` | [`sst_closed_loop.svg`](figures/sst_closed_loop.svg): closed-loop simulation of four mission classes (about 5 minutes and 1 GB of memory per class) |

## Run the viewer

The viewer is a single web page that reads precomputed monthly grids:

```bash
cargo run --release -p gravviz -- series data/costg data/baseline data/TN-14_C30_C20_GSFC_SLR.txt data/TN-13_GEOC_CSR_RL0603.txt viewer/data
python3 -m http.server --directory viewer 8000
```

Then open <http://localhost:8000>. A local server is needed because browsers block a page opened
straight from disk from reading files next to it.

## What is checked

- **`geopotential`** is tested against:
  - closed-form Legendre functions;
  - the addition theorem to degree 2700 at 89.9999°;
  - the Legendre equation;
  - finite differences;
  - Laplace's equation;
  - the IERS 2010 tide offset.

  On the real degree-300 GOCO06s model it agrees with an independent implementation to about
  1e-15.
- **`gravity-field`** is tested on real COST-G, GOCO06s and NASA files, and on a Swarm-layout file.
- **`satgrav`** reproduces:
  - Zingerle et al. (2024), Fig. 1;
  - the published CARIOQA, GRICE and Trimeche et al. (2019) instrument numbers.

  Its closed-loop simulator recovers a known field without noise. With noise, a χ² test confirms
  that the actual errors match the predicted ones.

Run everything with:

```bash
cargo test --workspace
```

## Limits

- **Glacial isostatic adjustment is not removed**, so Hudson Bay and Scandinavia show crustal
  rebound as water gain.
- **No leakage correction yet.**
- **The simulator has no time-variable gravity**, so it cannot show temporal aliasing.

The [research report](reports/Satellite%20gravimetry%20data%20and%20Rust%20tooling.md) covers the
data sources, the quantum missions (CARIOQA, NASA's gradiometer pathfinder), the formulas and the
plan; its sources are in [`research_notes/`](research_notes/).

## Licence and data credits

The code is licensed under either of [MIT](LICENSE-MIT) or [Apache-2.0](LICENSE-APACHE) at your
option.

Data credits:
- **COST-G monthly fields:** CC BY 4.0, doi:10.5880/COST-G.ICGEM_02_L2.
- **GOCO06s:** CC BY 4.0, doi:10.5880/ICGEM.2019.002.
- **GRACE technical notes:** NASA, no use constraints.
- **Coastlines:** Natural Earth, public domain.
- **Stripe removal:** ported from [gravity-toolkit](https://github.com/polargeodesy/gravity-toolkit)
  (MIT).
