# geopotential

Earth's gravity field from spherical-harmonic coefficients. Read a model from an ICGEM `.gfc`
file, from static models such as GOCO06s to every monthly GRACE field. Then evaluate the
potential, the gravity vector or the gravity-gradient tensor at any point, or draw maps of the
geoid, gravity anomalies or water height. Pure Rust, no required dependencies.

![The geoid from the GOCO06s satellite model](https://raw.githubusercontent.com/FusionVault/physica/main/figures/geoid_goco06s.svg)

*The geoid from GOCO06s, drawn with this crate.*

## Install

```bash
cargo add geopotential
```

## Example: a field and what it does

```rust
use geopotential::{Coefficients, Functional};

// A field with a single term: J2, Earth's flattening.
let mut field = Coefficients::zeros(3.986004415e14, 6_378_136.3, 2);
field.set(2, 0, -4.84165e-4, 0.0);

// Geoid height at the equator and at the North Pole.
let lats = [0.0, 90f64.to_radians()];
let geoid = field.grid(&Functional::GeoidHeight, field.radius, &lats, &[0.0]);
assert!(geoid[0] > 0.0); // the equatorial bulge
assert!(geoid[1] < 0.0); // flatter at the poles

// Add the central term (C̄00 = 1) and find gravity on a satellite 400 km above northern Italy.
field.set(0, 0, 1.0, 0.0);
let r = 6_378_137.0 + 400e3;
let (lat, lon) = (45f64.to_radians(), 7f64.to_radians());
let [north, east, up] = field.partials(r, lat, lon).acceleration_neu(r, lat);
assert!(up < -8.0); // pointing down, about 8.7 m/s²
assert!(north.abs() > 0.0 && east == 0.0);
```

## Example: read a real model

Download a model first, for example `GOCO06s.gfc` from the
[ICGEM model table](https://icgem.gfz.de/tom_longtime).

```rust,no_run
use geopotential::{julian_epoch, read_gfc};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let gfc = read_gfc(&std::fs::read_to_string("GOCO06s.gfc")?)?;

    // Time-variable models apply their trends and annual terms; static ones ignore the date.
    let field = gfc.model.at(julian_epoch(2026, 7, 16, 12, 0));
    println!("{} to degree {}, C20 = {:e}", gfc.name.unwrap_or_default(), field.max_degree(), field.c(2, 0));
    Ok(())
}
```

![The error of every coefficient of the July 2026 GRACE-FO field](https://raw.githubusercontent.com/FusionVault/physica/main/figures/error_triangle_2026-07.svg)

*The calibrated error of each coefficient of the July 2026 COST-G field, read with `read_gfc`.*

## What is inside

| Item | What it does |
|---|---|
| `Legendre` | Normalized Legendre functions and their derivatives, stable to degree 2700 at every latitude ([Holmes & Featherstone 2002](https://doi.org/10.1007/s00190-002-0216-2)) |
| `Coefficients` | The coefficients with their `GM`, radius and tide system; degree variances, rescaling, truncation and tide-system conversion ([IERS Conventions 2010](https://www.iers.org/TN36), ch. 6) |
| `TimeVariable` | Models with trends and periodic terms, evaluated at any date ([ICGEM](https://icgem.gfz.de)) |
| `partials`, `second_partials` | Potential, gravity vector and gravity-gradient tensor at a point; the geopotential and gravity with the centrifugal term, and the deflections of the vertical |
| `normal` | The normal gravity of a level ellipsoid from its four defining constants: GRS80 and WGS84 as Somigliana's formula on and above the ellipsoid, `U₀`, the zonal harmonics and the normal field as coefficients, Bruns' formula with its zero-degree term, and dynamic and normal heights from geopotential numbers ([Moritz 2000](https://doi.org/10.1007/s001900050278); Heiskanen & Moritz 1967) |
| `grid` with `Functional` | Maps of geoid height, gravity anomaly or disturbance, radial gradient or equivalent water height ([Barthelmes 2013](https://doi.org/10.2312/GFZ.b103-0902-26); [Wahr et al. 1998](https://doi.org/10.1029/98JB02844)) |
| `Coefficients::analysis`, `gauss_grid`, `gauss_legendre` | The exact inverse of `grid`: any functional on a Gauss–Legendre grid back to coefficients, by Gauss–Legendre quadrature |
| `slepian` | Slepian functions of a polar cap from Grünbaum's commuting tridiagonal matrix, and of any region (a cap anywhere, a polygon) from its concentration matrix, integrated exactly in longitude and on clustered Gauss–Legendre panels in colatitude; the basin kernel that averages a field over a region, with the least leakage ([Simons, Dahlen & Wieczorek 2006](https://doi.org/10.1137/S0036144504445765); [Swenson & Wahr 2002](https://doi.org/10.1029/2001JB000576)) |
| `read_gfc`, `Gfc` | ICGEM `.gfc` models, static or time-variable (icgem1.0 and 2.0), with their errors ([format](https://icgem.gfz.de/docs/ICGEM-Format-2023.pdf)) |
| `julian_epoch`, `parse_date` | The decimal-year dates these files use |
| `bodies` | Potential, attraction and gradient tensor of a uniform body in closed form: a rectangular prism ([Nagy, Papp & Benedek 2000](https://doi.org/10.1007/s001900000116)), finite on its faces and edges, with a Gauss–Legendre cubature of Bernstein-ellipse order outside it that holds 10⁻¹⁴ on needles and plates where the corner sums lose 10⁻⁸ ([Fukushima 2019](https://doi.org/10.1093/gji/ggz449) diagnosed the loss); and a triangulated polyhedron ([Werner & Scheeres 1997](https://doi.org/10.1007/BF00053511)), whose trace is an inside test, with a mesh check and a cubature over its tetrahedra far away |
| `planar` | Level grids in the wavenumber domain: upward continuation, reduction to the pole (plain or Wiener-damped) and to the equator, any change of field direction, Baranov pseudo-gravity (Blakely 1995, ch. 12), or any kernel through `Grid::filter`, with edge padding and a cosine taper against the periodic wrap; and Euler deconvolution of a profile or a grid for a source's position and the background, with its covariance, in moving windows, and as Euler inversion, which fits the data too and picks the structural index ([Reid & Thurston 2014](https://doi.org/10.1190/geo2013-0235.1); [Uieda et al. 2025](https://doi.org/10.1093/gji/ggaf114)) |
| `special` | Complete elliptic integrals `K` and `E` by the arithmetic–geometric mean; Bulirsch's general complete integral `cel` (Bulirsch 1969); Carlson's symmetric `R_F`, `R_D`, `R_J`, `R_C` with their principal values ([Carlson 1995](https://arxiv.org/abs/math/9409227)); Legendre's incomplete `F`, `E` and `Π` for any `m`, with `Π` beyond its pole as a principal value ([DLMF §19.25](https://dlmf.nist.gov/19.25)); and `ring`, the average of `1/r` over a ring and its gradient, which rings of mass, charge and current share |
| `linalg` | Pivoted Gaussian elimination for general systems, Cholesky for the normal equations of least squares, and symmetric eigenproblems by Jacobi's method or, for large matrices, Householder reduction and the implicit QL iteration |
| `fft` | The radix-2 fast Fourier transform, for time series, and `fftn` in any number of dimensions with `wavenumber`, for gridded fields |

**Conventions:**
- angles in radians, geocentric latitude;
- fully normalized coefficients with no Condon–Shortley phase, as in ICGEM files;
- the optional `serde` feature serializes `Coefficients`.

Models are free from [ICGEM](https://icgem.gfz.de); each COST-G monthly series also comes as one
zip.

## Applications

- **Height systems.** Turning a GNSS ellipsoidal height into a usable one needs the geoid; this
  is the synthesis behind national models such as GEOID2022 and AUSGeoid2020, and `normal` turns
  a levelled geopotential number into a dynamic or normal height.
- **Orbit design and determination.** `J2` from a real model sets the nodal precession that makes
  an orbit sun-synchronous; `examples/sun_synchronous_orbit.rs` reproduces Landsat's 98.2°.
- **Inertial navigation.** The deflection of the vertical is what an accelerometer mistakes for
  horizontal acceleration, and compensating it from a model is standard practice.
- **Gradiometry.** `SecondPartials::tensor_neu` gives the tidal tensor a gradiometer measures, in
  which a vehicle's own acceleration cancels.
- **Reading published models.** `read_gfc` takes any ICGEM `.gfc` file, static or time-variable.

## How it is checked

- **Legendre functions:** against closed forms, the addition theorem to degree 2700 close to the
  pole, and the Legendre equation.
- **Derivatives:** against finite differences.
- **The field:** against Laplace's equation, and the IERS tide offset of C̄20.
- **An independent implementation:** the degree-300 GOCO06s model agreed with the former geo5d
  crate's synthesis to about 1e-15.
- **Normal gravity:** every derived constant that Moritz (2000) prints for GRS80 and the WGS84
  document prints for WGS84 (γ_e, γ_p, k, U₀, J₂–J₈, the even zonals), to their digits; the
  synthesized normal field is Somigliana's gravity to 10⁻¹⁰ m/s²; the free-air gradient 0.3086
  mGal/m.
- **Bodies:** Fukushima's (2019) Table 4 prism to 10⁻¹⁴; a 1×1×100 needle and a 100×100×1 plate
  to 3×10⁻¹⁴ from 0.3 off their sides to 100 Brillouin radii; a cube mesh against the prism
  across the polyhedron's hand-over and out to 1000 sizes; a dented, non-convex mesh.
- **Special functions:** Carlson's check values and principal values; Abramowitz & Stegun's
  Tables 17.5, 17.6 and 17.9 within their rounding, DLMF's special values and periods, Π beyond
  its pole against DLMF 19.25.16.
- **Planar:** padding against the periodic wrap; reduction to the equator at 10° inclination to
  1%; Euler deconvolution exact for a dipole, a point mass and a line mass, and Euler inversion's
  misfit picking the structural index.
- **Slepian functions:** a cap away from the pole has the polar cap's spectrum to 10⁻¹²; a box
  and a triangle have their exact areas and Shannon numbers; the basin mean of P̄₂₀ over a cap in
  closed form.
- **Real files** from COST-G and GOCO06s, which bend the format in practice: missing header
  markers, comment and warning lines, `errors empirical`, and free text before the header. A
  malformed record fails with its line number. The fixtures' sources are listed in
  `tests/data/README.md`.

One of eight crates in [physica](https://github.com/FusionVault/physica):
[`physica`](https://crates.io/crates/physica) (all of them, plus cross-domain analogies) ·
[`metrology`](https://crates.io/crates/metrology) (constants) ·
[`geopotential`](https://crates.io/crates/geopotential) (the gravity field) ·
[`gravimetry`](https://crates.io/crates/gravimetry) (missions and GRACE) ·
[`gravitational`](https://crates.io/crates/gravitational) (waves) ·
[`gravitics`](https://crates.io/crates/gravitics) (frontier gravity) ·
[`magnetics`](https://crates.io/crates/magnetics) (magnetism and MHD) ·
[`plasmoid`](https://crates.io/crates/plasmoid) (fusion plasma).
