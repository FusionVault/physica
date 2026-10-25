# magnetics

Magnetism and magnetohydrodynamics, computed: Earth's magnetic field from the World Magnetic
Model, the numbers that decide how a magnetised plasma behaves, magnet engineering up to its
stress limits, space weather from the solar wind to the power grid, and the machines that turn
conducting flow into power or thrust.

Earth's main field is a potential field in spherical harmonics, exactly like its gravity field, so
this crate evaluates it with [`geopotential`](https://crates.io/crates/geopotential)'s engine
rather than a second one: give it `GM = a²` and `C̄ₙₘ = gₙᵐ/√(2n+1)` and the same code serves both.
All 100 rows of NOAA's official WMM2025 test file reproduce to 0.001 nT.

## Install

```bash
cargo add magnetics geopotential metrology
```

## Example: a compass correction anywhere on Earth

```rust
use magnetics::main_field::wmm2025;

// Declination, inclination and total field in London, mid-2026.
let e = wmm2025().elements(2026.5, 51.5074f64.to_radians(), (-0.1278f64).to_radians(), 0.0);
println!("D {:+.2}°, I {:+.2}°, F {:.0} nT", e.declination.to_degrees(), e.inclination.to_degrees(), e.f);
assert!((0.0..3.0).contains(&e.declination.to_degrees()));

// The model also gives the dipole and where its axis meets the sphere.
let (lat, lon) = wmm2025().geomagnetic_pole(2025.0);
assert!((lat.to_degrees() - 80.79).abs() < 0.01);
```

## Example: the numbers that set a plasma's behaviour

```rust
use magnetics::mhd::{alfven_speed, beta, mass_density, perpendicular_shock};
use metrology::M_P;

// PlasmaPy's reference case: 0.014 T through 5×10¹⁹ protons per m³.
let v = alfven_speed(0.014, mass_density(5e19, M_P, 1.0));
assert!((v - 43_173.870).abs() < 1e-3);

// A strong perpendicular shock compresses the plasma, approaching 4 for a monatomic gas.
let x = perpendicular_shock(8.0, 1.0, 5.0 / 3.0).unwrap();
assert!((x - 3.6154).abs() < 1e-3);
assert!(perpendicular_shock(1.0, 1.0, 5.0 / 3.0).is_none()); // too slow to shock
assert!(beta(3.9e6, 1.64, 13.3e-9) < 0.02); // the solar wind at low beta
```

## Example: what a magnet can hold, and what a storm induces

```rust
use magnetics::magnets::{larmor_frequency, lorentz_line_load, virial_specific_energy};
use magnetics::space_weather::{dynamo_field, geoelectric_halfspace, standoff, dynamic_pressure};

// The CFS/MIT 20 T model coil: 40.5 kA in 20.1 T loads the conductor at 814 kN/m.
assert!((lorentz_line_load(40.5e3, 20.1) / 1e3 - 814.05).abs() < 0.01);
// The virial theorem caps stored energy per kilogram of structure: 12.5 kJ/kg in steel at 100 MPa.
assert!((virial_specific_energy(100e6, 8000.0) - 12_500.0).abs() < 1.0);
// A 3 T MRI scanner's proton frequency.
assert!((larmor_frequency(3.0) / 1e6 - 127.73).abs() < 0.01);

// A typical solar wind puts the magnetopause 10 Earth radii upstream.
assert!((standoff(29_733.4e-9, dynamic_pressure(5e6, 400e3), true) - 10.08).abs() < 0.01);
// A 1000 nT variation over resistive ground induces 13 V/km, above NERC's 8 V/km benchmark.
let (e, _) = geoelectric_halfspace(1000e-9, 300.0, 1e-4);
assert!(e * 1e3 > 8.0);
assert!(dynamo_field(2.04, 3e-5) > 0.0); // a neutral wind across the field drives a current
```

## What is inside

| Module | What it computes |
|---|---|
| `main_field` | Earth's main field from WMM2025 (bundled) or any NOAA `.COF` model, through `geopotential`'s engine: the field elements D, I, H, F, X, Y, Z and their rates, the dipole moment, the geomagnetic pole, and `MainField` as a [`geopotential::Field`] |
| `mhd` | Alfvén and magnetosonic speeds, plasma beta, magnetic pressure, magnetic diffusivity and Reynolds number, Hartmann, Prandtl and Elsasser numbers, free decay time, inertial length and gyroradius, the ion spectral break ([Chen et al. 2014](https://doi.org/10.1002/2014GL062009)), the perpendicular shock jump, and parallel line forces beside their gravitational twins |
| `stability` | Kadomtsev's sausage and kink criteria, magnetic Rayleigh–Taylor cut-off, Furth–Killeen–Rosenbluth tearing growth |
| `flow` | The Hartmann velocity profile and pressure gradient, and the thin-walled pipe limit that sizes liquid-metal blankets |
| `dynamo` | The αΩ dynamo number and growth rate, with the published laboratory self-excitation thresholds (Riga, Karlsruhe, VKS) |
| `anomaly` | Magnetic anomalies of uniformly magnetised bodies by Poisson's relation, `B = (μ₀/4π)(T − I tr T)·M` from any body's gravity-gradient tensor in `geopotential::bodies`, inside or out; the sphere and the infinite cylinder in closed form; and the total-field reading `\|F + B\| − \|F\|` without cancellation; the demagnetising factors of any ellipsoid in one Carlson integral and the magnetisation `(I + χN)⁻¹(χH₀ + M_r)` with the Königsberger ratio ([Chen, Brug & Goldfarb 1991](https://doi.org/10.1109/20.102932); [Takahashi & Oliveira 2017](https://doi.org/10.5194/gmd-10-3591-2017)) |
| `magnets` | Loop, solenoid and Helmholtz fields; Lorentz line load, stored energy and the virial limit; skin depth; diamagnetic levitation and Inductrack lift-to-drag; the force on a susceptible body in a field gradient; MRI Larmor frequency, SAR limits and magnetometer energy resolution; flux compression |
| `space_weather` | Solar-wind dynamic pressure, magnetopause standoff (Chapman–Ferraro and [Shue et al. 1998](https://doi.org/10.1029/98JA01103)), the Parker spiral, ring-current energy from Dst, the plasmapause, dipole L-shells, the induced geoelectric field and apparent resistivity, the ionospheric dynamo, and the NERC benchmarks |

**Conventions:** SI throughout, with nanotesla for geomagnetic fields; angles in radians;
geodetic latitude and height above the WGS84 ellipsoid for the field elements.

## Applications

- **Navigation and heading.** The World Magnetic Model is the standard for the US and UK
  defence ministries, NATO and the International Hydrographic Organization, and it is what a
  phone compass corrects by. `main_field::elements` gives declination and its drift anywhere.
- **Directional drilling.** Magnetic surveying needs declination to about 0.1°, because that
  error alone moves the bottom of a well ~8.7 m over a 5 km reach.
- **Grid resilience.** Operators plan for NERC TPL-007's once-in-a-century geoelectric field;
  `space_weather` gives the induced field for a given ground conductivity and storm.
- **Fusion and liquid-metal engineering.** `mhd`, `stability` and `flow` give the stability
  limits and the pressure drop that size a blanket channel.
- **Magnet design.** `magnets` gives the field, the force, the stored energy and the structural
  limit that decides what a coil can survive.
- **Medicine and sensing.** Larmor frequencies, SAR limits and magnetometer energy resolution.

## How it is checked

- **The main field:** every one of the 100 rows of NOAA's `WMM2025_TestValues.txt` (bundled) to
  0.001 nT and 0.005°, plus the report's Table 6 and the published geomagnetic pole and dipole
  moment. The Schmidt and fully-normalised evaluations agree to 6×10⁻¹¹ nT.
- **MHD:** PlasmaPy's published Alfvén speeds exactly; a coronal loop; ITER; Chen et al.'s
  observed solar-wind spectral breaks (0.58 ± 0.19 Hz and 0.223 ± 0.068 Hz); the shock jump
  against the conservation laws and its strong-shock limit of 4.
- **Stability and flow:** Kadomtsev's marginal kink radius on the Bennett profile (`r* = √3 r_p`,
  `β* = 1/3`); Khomenko et al.'s 38 km and 155 km prominence wavelengths; the FKR coefficient
  0.547; KIT's published duct pressure gradient.
- **Dynamos:** the αΩ onset at dynamo number 2 recovered by bisection on the dispersion relation.
- **Magnets:** the CFS/MIT coil's 815 kN/m; Berry and Geim's levitation condition and 16.5 T
  field; copper skin depths; CODATA's proton gyromagnetic ratio; Dang et al.'s 44ħ.
- **Space weather:** the factor `4^(1/6)` between the compressed and uncompressed standoff; Shue
  et al.'s values; the Dessler–Parker–Sckopke energy for the May 2024 and March 1989 storms; the
  magnetotelluric identity recovering the ground resistivity exactly.

Every example above is compiled and run as a test, and so are the programs in `examples/`.

One of eight crates in [physica](https://github.com/FusionVault/physica):
[`physica`](https://crates.io/crates/physica) (all of them, plus cross-domain analogies) ·
[`metrology`](https://crates.io/crates/metrology) (constants) ·
[`geopotential`](https://crates.io/crates/geopotential) (the gravity field) ·
[`gravimetry`](https://crates.io/crates/gravimetry) (missions and GRACE) ·
[`gravitational`](https://crates.io/crates/gravitational) (waves) ·
[`gravitics`](https://crates.io/crates/gravitics) (frontier gravity) ·
[`magnetics`](https://crates.io/crates/magnetics) (magnetism and MHD) ·
[`plasmoid`](https://crates.io/crates/plasmoid) (fusion plasma).
