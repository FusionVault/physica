# quantum-gravimeter

Rust tools for satellite gravimetry. Read the gravity fields measured by GRACE, GRACE-FO and
GOCE, turn them into maps of water, ice and earthquakes, and simulate the quantum (cold-atom)
sensors of future missions.

![Water storage anomaly in July 2026, from GRACE-FO](figures/water_2026-07.svg)

## The crates

| Crate | What it does |
|---|---|
| [`geopotential`](crates/geopotential) · [docs](https://docs.rs/geopotential) | Earth's gravity field from spherical harmonics: gravity, geoid and water height at a point or on a map |
| [`gravity-field`](crates/gravity-field) · [docs](https://docs.rs/gravity-field) | Reads gravity-field files: ICGEM `.gfc` models and the GRACE technical notes |
| [`gravimetry`](crates/gravimetry) · [docs](https://docs.rs/gravimetry) | GRACE processing (smoothing, stripe removal, water height), earthquake gravity models and mass transport |
| [`satgrav`](crates/satgrav) · [docs](https://docs.rs/satgrav) | Gravity missions and quantum sensors: noise models, error predictions, a satellite-pair simulator, and the sensors as gravitational-wave detectors |
| [`gwsens`](crates/gwsens) · [docs](https://docs.rs/gwsens) | Gravitational-wave sensitivity: strain conventions, LIGO, Einstein Telescope and LISA noise curves, horizons, stochastic backgrounds, high-frequency limits and lab-gravity reality checks |
| [`gravviz`](apps/gravviz) | The command-line tool that makes every figure here (not published) |
| [`viewer`](viewer) | The GRACE Mass Atlas, an interactive web page for all of it |

[`geo5d`](https://crates.io/crates/geo5d) (heights, normal gravity, relativistic geodesy) uses
`geopotential` for its gravity.

## Quick start

```bash
cargo add geopotential gravity-field
```

Download a gravity model, for example `GOCO06s.gfc` from the
[ICGEM model table](https://icgem.gfz.de/tom_longtime), then:

```rust
use gravity_field::read_gfc;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let model = read_gfc(&std::fs::read_to_string("GOCO06s.gfc")?)?.model.base;

    // Gravity 400 km above northern Italy (radius in metres, angles in radians).
    let (r, lat, lon) = (6_378_137.0 + 400e3, 45f64.to_radians(), 7f64.to_radians());
    let [north, east, up] = model.partials(r, lat, lon).acceleration_neu(r, lat);

    println!("gravity: {north:.2e} north, {east:.2e} east, {up:.4} up (m/s²)");
    Ok(())
}
```

Each crate's README has more examples.

## Gallery

| | |
|---|---|
| ![Geoid](figures/geoid_goco06s.svg) **The geoid** from GOCO06s | ![Kamchatka](figures/e1_kamchatka.svg) **An earthquake from space:** Kamchatka 2025 against the USGS fault model |
| ![El Niño](figures/w1_land_storage.svg) **El Niño drains the continents:** global land water and a registered forecast | ![Noise](figures/q0_noise.svg) **What makes a month noisy:** solar activity and satellite modes |
| ![Instrument noise](figures/noise_spectra.svg) **Quantum against classical sensors** ([Zingerle et al. 2024](https://doi.org/10.1029/2024EA003630), Fig. 1) | ![Closed loop](figures/sst_closed_loop.svg) **Simulated missions** against a real month |

## Research

The [research report](reports/Gravimetry%20frontiers%20and%20hidden%20connections.md) proposes
new tests of how gravity data connects to climate and earthquakes; a
[second report](reports/Gravitational%20waves%20and%20gravity%20beams.md) covers gravitational
waves, high-frequency waves and gravity beams. Each forecast is committed as a
git tag before the data that decides it exists, and results are reported whether they hit or
miss.

| Test | Question | Status |
|---|---|---|
| [W1](prereg/W1-2026-10-24.md) | How much water does land lose in an El Niño, and does it come back? | Registered 2026-10-24; scored from early 2027 |
| [E1](results/E1-2026-10-24.md) | Can monthly gravity see the 2025 Kamchatka earthquake, and is there more than the slip explains? | Detected. Matches the USGS fault model, with no excess |
| [Q0](results/Q0-2026-10-24.md) | What makes a monthly field noisy? | Solar activity and satellite modes |
| [Screen](figures/screen_indices.svg) | Which climate patterns move global land water? | Only El Niño, once 27 tests are accounted for |
| [Q2](results/Q2-2026-10-24.md) | Are the stated errors of the monthly fields right? | Too small, by about 2× above degree 10 |
| [M1](results/M1-2026-10-24.md) | Can a blind search of every month and place find earthquakes? | Sumatra and Tohoku, with no false alarms; smaller ones need a targeted test |
| [T1](prereg/T1-2026-10-24.md) | Does GRACE-FO's laser ranging feel the 2022 Tonga eruption's pressure wave? | Registered 2026-10-24 before any data was fetched; feasible on paper (SNR 25–55 per good crossing) |

**Findings so far:**
- **El Niño drains the land.** Global land loses about 620 km³ of water per °C of El Niño
  (p < 0.001). La Niña's effect is smaller and not significant. The forecast for the 2026 El
  Niño is a loss of about 1,300 km³ by December, about 3.7 mm of sea level.
- **The Kamchatka earthquake is seen from orbit.** The Mw 8.8 earthquake left a
  −3.5 / +3.5 µGal dipole, larger than at any of 100 random places.
  - It has 0.91 ± 0.16 times the amplitude of the USGS fault model
    ([Okubo 1992](https://doi.org/10.1029/92JB00178) half-space), with correlation 0.59, so the
    slip explains it.
  - After-slip is not significant yet (p = 0.08).
- **The sun adds noise.** Solar activity raises the short-wavelength noise of a month by 30% per
  standard deviation of solar flux. At the 300 km scale used above, no factor matters, so solar
  maximum does not explain those results.
- **A blind search finds the two biggest earthquakes and nothing false.** Searching every month
  and place, with no catalogue, found Sumatra 2004 and Tohoku 2011 and no false alarms in 254
  months. No flood or drought made a step that sudden.
- **The stated errors are too small.** Above degree 10, the monthly fields scatter 2 to 2.5 times
  more than their stated errors. Scaled up, a Kalman smoother predicts withheld months with
  calibrated errors (rms 1.07), so it can fill the gaps in the record and map its uncertainty.

## Reproduce everything

You need Rust 1.85 or later. All the data is free, and none of it needs an account.

**1. Download the data** (about 150 MB):

```bash
mkdir -p data && cd data

# Every COST-G monthly field, 2002–2026, as one zip
curl -L -o costg.zip https://icgem.gfz.de/getseries/02_COST-G_/Grace-Grace-FO_RL02.1
unzip -q -j costg.zip -d costg
mkdir -p baseline && cp costg/GSM-2_200[4-9]*.gfc baseline/

# GRACE technical notes: C20 and C30 from laser ranging (TN-14), degree 1 (TN-13)
curl -LO https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-docs/gracefo/open/docs/TN-14_C30_C20_GSFC_SLR.txt
curl -LO https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-docs/gracefo/open/docs/TN-13_GEOC_CSR_RL0603.txt

# A static model
curl -L -o GOCO06s.gfc https://icgem.gfz.de/getmodel/gfc/32ec2884630a02670476f752d2a2bf1c395d8c8d6d768090ed95b4f04b0d5863/GOCO06s.gfc

# Climate indices (NOAA PSL), and solar and geomagnetic activity (GFZ Potsdam)
mkdir -p indices
for f in nina34.anom pdo amon.us nao ao aao qbo solar; do
  curl -L -o indices/$f.data https://psl.noaa.gov/data/correlation/$f.data
done
curl -L -o indices/dmi.had.long.data https://psl.noaa.gov/gcos_wgsp/Timeseries/Data/dmi.had.long.data
curl -L -o indices/Kp_ap_Ap_SN_F107_since_1932.txt https://kp.gfz-potsdam.de/app/files/Kp_ap_Ap_SN_F107_since_1932.txt

# Every earthquake of magnitude 8 or more since 2002 (USGS)
curl -L -o usgs_m8_2002_2026.csv "https://earthquake.usgs.gov/fdsnws/event/1/query?format=csv&starttime=2002-04-01&endtime=2026-08-01&minmagnitude=8.0&orderby=time-asc"

# The USGS fault model of the 2025 Kamchatka earthquake
mkdir -p usgs_us6000qw60
for f in coulomb.inp surface_deformation.disp; do
  curl -L -o usgs_us6000qw60/$f https://earthquake.usgs.gov/product/finite-fault/us6000qw60_1/us/1763589364802/$f
done
cd ..
```

**2. Make the figures.** Each command writes an SVG and, for most, a CSV of the plotted values.
Set three shortcuts first:

```bash
MONTH=data/costg/GSM-2_2026182-2026212_GRFO_COSTG_BF01_0201.gfc
TN14=data/TN-14_C30_C20_GSFC_SLR.txt
TN13=data/TN-13_GEOC_CSR_RL0603.txt
```

| Command (`cargo run --release -p gravviz -- …`) | Figure |
|---|---|
| `geoid data/GOCO06s.gfc figures` | [The geoid](figures/geoid_goco06s.svg) |
| `water-map $MONTH data/baseline $TN14 $TN13 figures` | [Water in July 2026](figures/water_2026-07.svg) |
| `triangle $MONTH figures` | [Error of each coefficient](figures/error_triangle_2026-07.svg) |
| `noise figures` | [Instrument noise](figures/noise_spectra.svg) |
| `gradiometry data/GOCO06s.gfc figures` | [Gradiometer missions](figures/gradiometry_errors.svg) |
| `sst figures $MONTH` | [Simulated missions](figures/sst_closed_loop.svg) (about 5 minutes and 1 GB per mission) |
| `w1 data/costg data/baseline $TN14 $TN13 data/indices/nina34.anom.data figures` | [El Niño and land water](figures/w1_land_storage.svg) |
| `screen figures/w1_land_storage.csv data/indices figures` | [Climate-index screen](figures/screen_indices.svg) |
| `q0 data/costg data/baseline $TN14 $TN13 data/indices/Kp_ap_Ap_SN_F107_since_1932.txt figures` | [Noise](figures/q0_noise.svg) |
| `e1 data/costg data/baseline $TN14 $TN13 data/usgs_us6000qw60 figures/q0_noise.csv figures` | [Kamchatka earthquake](figures/e1_kamchatka.svg) |
| `kalman data/costg data/baseline $TN14 $TN13 figures viewer/data` | [Are COST-G's errors right?](figures/kalman.svg), plus the gap-free record for the viewer |
| `m1 data/costg data/baseline $TN14 $TN13 figures/kalman_degrees.csv data/usgs_m8_2002_2026.csv figures viewer/data` | [Blind search for sudden mass changes](figures/m1_detector.svg) |
| `tonga figures` | [The Tonga Lamb wave as GRACE-FO would see it](figures/tonga_lamb_wave.svg) |

**3. Run the viewer.** The GRACE Mass Atlas has five views:
- **Water:** every month on a map or globe, with animated mass flow and an average year. The
  Kalman-smoothed record fills the gaps and shows its uncertainty.
- **El Niño:** the forecast, and the climate-index screen.
- **Earthquake:** an animated map of the Kamchatka earthquake.
- **Noise:** what makes a month noisy.
- **Missions:** the instrument figures.

```bash
cargo run --release -p gravviz -- series data/costg data/baseline $TN14 $TN13 viewer/data
cp figures/{w1_land_storage.csv,w1_monthly.csv,w1_forecast.csv,screen_indices.csv,e1_frames.json,q0_noise.csv,q0_effects.csv,noise_spectra.svg,sst_closed_loop.svg,gradiometry_errors.svg,error_triangle_2026-07.svg,kalman.svg} viewer/data/
python3 -m http.server --bind 127.0.0.1 --directory viewer 8000
```

Open <http://localhost:8000>, and stop the server with Ctrl-C when you are done.

## How it is checked

```bash
cargo test --workspace
```

The tests compare against:
- closed-form results;
- published tables, such as Okada's (1985) check values and the Wahr et al. (1998) Love numbers;
- an independent implementation (geo5d, to about 1e-15);
- the published figures of Zingerle et al. (2024), CARIOQA and GRICE;
- simulations that recover a known field.

Every example in the crate READMEs is compiled and run as a test.

**Limits:**
- **Glacial isostatic adjustment is not removed**, so Hudson Bay and Scandinavia show crustal
  rebound as water gain.
- **No leakage correction yet.**
- **The simulator has no time-variable gravity**, so it cannot show temporal aliasing.

## Data and sources

**Data:**
- **COST-G monthly fields:** [doi:10.5880/COST-G.ICGEM_02_L2](https://doi.org/10.5880/COST-G.ICGEM_02_L2), CC BY 4.0.
- **GOCO06s:** [doi:10.5880/ICGEM.2019.002](https://doi.org/10.5880/ICGEM.2019.002), CC BY 4.0.
- **GRACE technical notes:** NASA, via
  [PO.DAAC](https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-docs/gracefo/open/docs/).
- **Climate indices:** [NOAA PSL](https://psl.noaa.gov/data/climateindices/list/).
- **Solar and geomagnetic activity:** GFZ Potsdam,
  [doi:10.5880/Kp.0001](https://doi.org/10.5880/Kp.0001), CC BY 4.0.
- **Kamchatka fault model:** [USGS](https://earthquake.usgs.gov/earthquakes/eventpage/us6000qw60/finite-fault).
- **Coastlines and land:** [Natural Earth](https://www.naturalearthdata.com/).

**Key papers:**
- Wahr et al. (1998), [water from gravity](https://doi.org/10.1029/98JB02844).
- Swenson & Wahr (2006), [stripe removal](https://doi.org/10.1029/2005GL025285).
- Landerer & Swenson (2012), [GRACE Tellus maps](https://doi.org/10.1029/2011WR011453).
- Holmes & Featherstone (2002), [high-degree Legendre functions](https://doi.org/10.1007/s00190-002-0216-2).
- Okada (1985) [uplift](https://doi.org/10.1785/BSSA0750041135) and Okubo (1992)
  [gravity](https://doi.org/10.1029/92JB00178) from faults.
- Zingerle et al. (2024), [quantum accelerometers for gravimetry](https://doi.org/10.1029/2024EA003630).

The first [research report](reports/Satellite%20gravimetry%20data%20and%20Rust%20tooling.md)
covers the data, the missions and the formulas. Its sources are in
[`research_notes/`](research_notes/).
