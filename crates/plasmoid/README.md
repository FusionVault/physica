# plasmoid

Fusion and self-organising plasma physics, computed. Every function reproduces a published number,
and the proposals that promise more (cold fusion, Pais's fusion device, magnetically suppressed
bremsstrahlung, overunity plasmoids) are set against the limit they have to beat.

- **Getting through the Coulomb barrier:** tunnelling and Gamow peaks, Bosch–Hale reaction rates,
  the p–¹¹B fits, the Lawson criterion, electron screening, muon catalysis, and gravity's way
  of doing it in stars.
- **Pinches and the plasma focus:** Bennett equilibrium, the Pease–Braginskii current, Lee's
  plasma-focus scalings, and the pinch's gravitational twin, the self-gravitating filament.
- **Plasmoids:** Taylor relaxation (`∇×B = λB`, what a plasma settles into by itself),
  spheromaks, field-reversed configurations, and reconnection, slow (Sweet–Parker) and fast
  (plasmoid-mediated).
- **Aneutronic power:** Q-values, the neutrons D–³He still makes, and p–¹¹B's race against
  bremsstrahlung and synchrotron radiation.
- **Electricity without a steam cycle:** multi-stage direct converters, inductive recovery in
  pulsed machines, and the plant arithmetic.

## Install

```bash
cargo add plasmoid metrology
```

## Example: the Coulomb barrier, and how fusion gets past it

```rust
use plasmoid::fusion::{Reaction, M_D, M_T, bosch_hale, coulomb_barrier, gamow_energy, gamow_peak, reduced_mass};

// D–T touching at 3.24 fm: a 444 keV barrier.
let barrier = coulomb_barrier(1.0, 1.0, 3.243);
// At 10 keV the reactions happen at the Gamow peak, 31 keV: tunnelling, 14× below the top.
let (peak, _) = gamow_peak(gamow_energy(1.0, 1.0, reduced_mass(M_D, M_T)), 10.0);
println!("barrier {barrier:.0} keV, Gamow peak {peak:.1} keV");
// Bosch & Hale's reaction rate there.
assert!((bosch_hale(Reaction::DT, 10.0).unwrap() / 1.1362e-22 - 1.0).abs() < 1e-3);
```

## Example: the pinch and its gravitational twin

```rust
use metrology::{PARSEC, SOLAR_MASS};
use plasmoid::pinch::{bennett_line_density, magnetic_over_gravitational, ostriker_line_mass};

// The equation that holds a 1 MA Z-pinch together is the one that holds a molecular-cloud
// filament together, with gravity in place of magnetism (Bennett 1934, Ostriker 1964).
let cloud = ostriker_line_mass(10.0, 2.3) * PARSEC / SOLAR_MASS; // 16.6 solar masses per parsec
let pinch = bennett_line_density(1e6, 2000.0) * 3.3436e-27;     // deuterium at 2 keV, kg/m
assert!((cloud - 16.6).abs() < 0.2);
// In the pinch, magnetism beats the plasma's own gravity by 5.5×10²⁷.
assert!(magnetic_over_gravitational(1e6, pinch) > 1e27);
```

## What is inside

| Module | What it computes | Key sources |
|---|---|---|
| `fusion` | Coulomb barrier, Sommerfeld and Gamow factors, Gamow peak, Bosch–Hale reactivities and cross sections (D–T, D–D, D–³He), p–¹¹B (Nevins–Swain, Tentori–Belloni), Lawson, screening (first-order and exact), Salpeter, muon catalysis, Lane–Emden, and the low-energy fusion checks | Bosch & Hale 1992 ([implementation notes](https://notebook.community/MFEFormulary/MFEFormulary/ipython/BH_ReactionRates)); [Wurzel & Hsu 2022](https://arxiv.org/abs/2105.10954); [Tentori & Belloni 2023](https://doi.org/10.1088/1741-4326/acda4b) |
| `pinch` | Bennett temperature and line density, Pease–Braginskii, Lee's peak current and yield, snowplow speed, cyclotron energy, the Ostriker filament | [Klír 2005](https://arxiv.org/abs/physics/0703207); [Lee 2009](https://www.plasmafocus.net/IPFS/2009/09PP8APL.pdf); Ostriker 1964 (ApJ 140, 1056) |
| `relaxation` | Bessel-function model, spheromak eigenvalues, Taylor-state energy, FRC beta and s-parameter, Lundquist number, Sweet–Parker and plasmoid reconnection, the virial cap on a free plasmoid | [Taylor 1986](https://doi.org/10.1103/RevModPhys.58.741); [Tuszewski 1988](https://zenodo.org/records/1235740); [Loureiro et al. 2007](https://arxiv.org/abs/astro-ph/0703631) |
| `aneutronic` | Q-values, D–³He neutron fraction, relativistic bremsstrahlung, synchrotron, the p–¹¹B fusion-to-bremsstrahlung ratio | Putvinski, Ryutov & Yushmanov 2019 (Nucl. Fusion 59, 076018) |
| `discharge` | Gas discharges and the thin gas around them: Paschen breakdown and its minimum with air's Townsend constants, the mean free path, Peek's corona onset, the floating potential in a glow, ion mobility, gas density and molecular speed, and slip-corrected drag on a sphere moving or turning at any Knudsen number: whether a high-voltage rig "in vacuum" is free of ions, and how long a body in it keeps swinging | Raizer, *Gas Discharge Physics* (1991), table 4.1; Peek 1929; Hirschfelder, Curtiss & Bird 1954; Davies 1945; Lieberman & Lichtenberg |
| `conversion` | N-stage direct-collector efficiency, engineering gain of magnetic and inertial plants, pulsed inductive breakeven, laser p–¹¹B gain | Moir & Barr 1973 (Nucl. Fusion 13, 35); [Margarone et al. 2020](https://doi.org/10.3389/fphy.2020.00343) |

The research behind it, with every claim's status, is in the reports on the
[Coulomb barrier and the Schwinger limit](https://github.com/FusionVault/physica/blob/main/research/topics/coulomb-and-schwinger/report.md)
and on [Z-pinches, plasmoids and aneutronic fusion](https://github.com/FusionVault/physica/blob/main/research/topics/plasmoids-and-fusion/report.md).

## Applications

- **Judging a design point.** `tokamak` checks a machine against the safety-factor, density,
  beta and confinement limits it must respect; `examples/sparc.rs` runs ITER and SPARC.
- **Reading a shot.** `examples/triple_products.rs` places published results against the Lawson
  requirement at their own temperature, and `examples/fusion_yields.rs` turns a yield into the
  reactions and fuel behind it.
- **Choosing a fuel.** `aneutronic` and `conversion` weigh p–¹¹B and D–³He against bremsstrahlung
  and set out what direct conversion returns without a steam cycle.
- **Pinches and self-organisation.** `pinch` and `relaxation` give the equilibria, scalings and
  reconnection rates, with the gravitational twin of the Bennett relation alongside.

## How it is checked

- **Bosch–Hale:** all of UWFDM-1268's Table 6 (D–T, both D–D branches, D–³He at 1, 10 and
  100 keV) within 0.4%; the D–T cross-section peak (5.07 b); and the Maxwellian average of the
  cross section against the reactivity fit.
- **p–¹¹B:** the published Nevins–Swain and Tentori–Belloni tables, and their 30% and 50% gaps.
- **Lawson:** Wurzel & Hsu's 2.97×10²¹ keV s m⁻³ at 14 keV and the 4.3 keV ignition temperature.
- **Pinches:** Klír's 3.12 keV Bennett coefficient, the 1.4 MA Pease–Braginskii current, and Lee's
  Table I peak currents within 1.5%.
- **Plasmoids:** Ji & Daughton's Lundquist numbers and Bhattacharjee et al.'s rates.
- **Aneutronic power:** NRL's bremsstrahlung limit and Putvinski's synchrotron ratios.
- **Laser shots:** Margarone's and the PALS shots' energy returns.

One of eight crates in [physica](https://github.com/FusionVault/physica):
[`physica`](https://crates.io/crates/physica) (all of them, plus cross-domain analogies) ·
[`metrology`](https://crates.io/crates/metrology) (constants) ·
[`geopotential`](https://crates.io/crates/geopotential) (the gravity field) ·
[`gravimetry`](https://crates.io/crates/gravimetry) (missions and GRACE) ·
[`gravitational`](https://crates.io/crates/gravitational) (waves) ·
[`gravitics`](https://crates.io/crates/gravitics) (frontier gravity) ·
[`magnetics`](https://crates.io/crates/magnetics) (magnetism and MHD) ·
[`plasmoid`](https://crates.io/crates/plasmoid) (fusion plasma).
