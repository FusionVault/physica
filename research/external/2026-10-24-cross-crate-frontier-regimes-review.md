# Review: cross-crate frontier regimes (batch of 2026-10-24)

The batch is [`2026-10-24-cross-crate-frontier-regimes.md`](2026-10-24-cross-crate-frontier-regimes.md).
Each item is reviewed as the roadmap sets out: computed in its own terms, beside the standard
prediction and the measurements, with an attempt to verify it, and a verdict on what goes into
which crate. Numbers marked *(computed)* come from the check script
[`2026-10-24-cross-crate-frontier-regimes-checks.py`](2026-10-24-cross-crate-frontier-regimes-checks.py).

Verdicts: **covered** (the crates already do it), **new** (computable, goes into a crate),
**to source** (needs a primary source before coding), **out** (with the reason).

---

## Item 1. The G scatter and the magnetic-susceptibility blind spot

### The claim

CODATA's 16 measurements of G disagree beyond their error bars (σ_G/G = 2.2×10⁻⁵ after a 3.9×
expansion), and Lee et al. 2020 found a new systematic from test-mass magnetic susceptibility.
The batch proposes that environment-dependent physics biases each lab differently: a chameleon
or symmetron scalar whose range depends on local density and field energy, or a coupling between
the local geomagnetic gradient and gravitational mass (including gravitomagnetic or Einstein–Cartan
spin–torsion terms). It proposes regressing the 16 residuals on each lab's coordinates, local
magnetic field and local geopotential.

### What the crates already do

| Piece | Where | Status |
|---|---|---|
| The 16 CODATA 2022 inputs, their correlations, the GLS mean 6.674 300 4, χ² = 195.7, Birge ratio 3.61, the ×3.9 expansion to 22.1 ppm | `gravitics::anomalies::{CODATA_2022_G, CODATA_2022_G_CORRELATIONS, gls_mean, birge_ratio}`; example `gravitics/examples/codata_g.rs` | covered |
| Lee et al. 2020, λ < 38.6 µm, and the thesis's new susceptibility systematic | `gravitics::modified::YUKAWA_LIMITS`; the systematic is recorded in `research/topics/magnetics-and-mhd/notes/practical_applications.md` | covered (the systematic as a note) |
| The magnetic field vector at any lab and date (2025–2030) | `magnetics::main_field::wmm2025().elements(...)` | covered for current epochs; past epochs need IGRF-14 (download pending the user's decision) |
| Gravity, geoid and gravity gradient at a lab | `geopotential` with a real model (`read_gfc`) | covered |
| Einstein–Cartan spin–torsion density | `gravitics::modified::einstein_cartan_density`: spin–torsion rivals mass energy only near 3.6×10⁴⁷ g/cm³ | covered |
| Earth's gravitomagnetic field at a lab | `gravitics::gem::GravitomagneticDipole` | covered |

### Testing the claim with the table itself *(computed)*

A bias set by the site (geomagnetic latitude, local field, local geopotential) is the same for
every apparatus at that site. The CODATA table has two sites with repeated measurements:

- **HUST, Wuhan:** four inputs, HUST-05, HUST-09, HUSTT-18 and HUSTA-18, spanning **339 ppm**.
  HUSTT-18 and HUSTA-18 come from the same lab in the same year with two methods, and differ by
  **45 ppm, 2.7σ**.
- **BIPM, Sèvres:** BIPM-01 and BIPM-14 agree within **7.5 ppm, 0.16σ**, and both sit about
  +190 ppm above the mean.

So the scatter follows the **apparatus, not the place**: a site-level field cannot produce a
339 ppm spread within one site. This does not rule out magnetic effects. It says they must enter
through the apparatus (test-mass susceptibility, magnetised source masses, nearby steel), which is
exactly the class of systematic Lee et al. found.

### The magnetic force on a test mass *(computed)*

The force on a test mass of volume susceptibility χ and volume V in a field B with gradient
∂B/∂x is

  **F = (χV/μ₀) B ∂B/∂x**,

the same physics as diamagnetic levitation (`magnetics::magnets::levitation_product`). For a
0.2 kg gold test mass (χ = −3.4×10⁻⁵) against a Cavendish signal of 10 kg at 10 cm
(1.33×10⁻⁸ N):

| Field | B ∂B/∂x (T²/m) | Force (N) | Fraction of the signal |
|---|---|---|---|
| Earth's own dipole gradient (B = 50 µT, ∂B/∂x ≈ 3B/R⊕) | 1.2×10⁻¹⁵ | 3.3×10⁻¹⁹ | 2.5×10⁻¹¹ |
| What a 2.2×10⁻⁵ bias needs | 1.05×10⁻⁹ | 2.9×10⁻¹³ | 2.2×10⁻⁵ |
| For example, B = 100 µT with ∂B/∂x = 10 µT/m | 1.0×10⁻⁹ | — | ≈ 2×10⁻⁵ |

Earth's field is six orders of magnitude too weak to matter. A lab field of twice Earth's
with a 10 µT/m gradient (a weakly magnetised source mass or nearby steel) is enough. In a
torsion balance the relevant quantity is the part of this force that moves with the source
masses, so the bias depends on how each apparatus modulates its sources; that is why it is an
apparatus systematic.

### The other proposed couplings *(computed or covered)*

- **Gravitomagnetism at the lab.** Earth's equatorial field is 5.04×10⁻⁶ m/s² (Mashhoon
  convention). A test mass moving at 1 mm/s feels 2(v/c)B = 3.4×10⁻¹⁷ m/s² against the 6.7×10⁻⁸
  m/s² Cavendish acceleration of a 0.2 kg mass: 5×10⁻¹⁰ of the signal. Too small by 10⁴–10⁵.
- **Einstein–Cartan spin–torsion.** It matters near 10⁴⁷ g/cm³; lab matter is 10⁴⁶ times less
  dense. Negligible in standard Einstein–Cartan theory.
- **Chameleon and symmetron scalars.** These are the live possibility: their range depends on
  the ambient density, so the vacuum chamber's size and walls change the force. The quantitative
  test is the thin-shell factor of the test and source masses. **To source:** Khoury & Weltman
  2004 (Phys. Rev. D 69, 044026) for the chameleon thin-shell condition, Hinterbichler & Khoury
  2010 for the symmetron, and Burrage & Sakstein 2018 (Living Reviews in Relativity) for the
  current lab bounds, including atom-interferometer bounds.

### Verdicts and what goes into the crates

| Piece | Verdict | Target |
|---|---|---|
| Pairwise tension `(a − b)/√(u_a² + u_b²)` and the same-site tests (HUST 339 ppm, HUSTT–HUSTA 2.7σ, BIPM 0.16σ) | **new** | `gravitics::anomalies::tension`, with the site tests |
| Susceptibility force `F = (χV/μ₀) B ∂B/∂x` and the bias it causes in G | **new** | `magnetics::magnets::susceptibility_force`; a `gravitics::anomalies` doc example of the G bias |
| Gravitomagnetic and Einstein–Cartan couplings at the lab | **covered**, shown negligible | example tying `gem` to the CODATA signal |
| Chameleon and symmetron thin-shell screening and lab bounds | **to source**, then **new** | `gravitics::modified` (screened fifth forces) |
| The 16-lab regression on site coordinates, local field and geopotential | **new**, as a formal null test | `lab` subcommand: needs each lab's published coordinates and epoch, and IGRF-14 for pre-2025 fields; with 16 points and three regressors, report p-values corrected for the number of regressors tried |
| The claim that a correlation would be "a publishable frontier result" | the same-site test already shows a site-level effect cannot be the main cause; an apparatus-level magnetic systematic is consistent with the data and with Lee et al. | — |

---

## Item 2. High-frequency gravitational waves from HED plasmas and the Gertsenshtein effect

### The claim

High-energy-density experiments (NIF's 8.6 MJ shot, MagLIF and FuZE Z-pinches) and strong static
fields (SPARC's 12–20 T magnets) should radiate or convert high-frequency gravitational waves
(MHz–THz): through the quadrupole of an asymmetric implosion (the third time derivative of the
mass quadrupole), and through Gertsenshtein–Zel'dovich conversion of an intense laser in a strong
magnetic field, `P(γ ↔ g) ∝ (G/c⁴)B₀²L²`. The vacuum's Euler–Heisenberg nonlinearity adds a
small effective stress-energy. The batch proposes computing the lab strain h(f) and comparing it
with resonant-cavity and magnon detectors. It cites JASON's 2008 HFGW review (JSR-08-506) and
the DIRDs.

### What the crates already do

| Piece | Where | Status |
|---|---|---|
| Gertsenshtein conversion `P = 4πGB²L²/(μ₀c⁴)` | `gravitational::hf::gertsenshtein_probability` | covered |
| The same conversion written as an axion coupling `√(16πG)` (identical to 10⁻⁹) | recomputed in the magnetics notes; an analogy for `physica` | covered (notes) |
| Graviton counts in a wave | `gravitational::quanta::graviton_count`, `occupation_per_cell` | covered |
| Wave energy flux, rotor power, the beam's near-field gravity | `gravitational::laboratory` | covered |
| Euler–Heisenberg energy correction `(α/45π)(E/E_S)²` and vacuum birefringence | `gravitics::strong_field::{euler_heisenberg_fraction, birefringence}` | covered |
| The Schwinger field, CoReLS's record intensity and the field ratio | `gravitics::strong_field`; example `laser_vs_schwinger` | covered |
| High-frequency detector limits (QUEST, Holometer, CAPP-12T at 5.3 GHz: h₀ < 3.9×10⁻²¹, MAGE, light-shining-through-a-wall, CAST) | `gravitational::hf::REGISTRY` | covered |
| Li–Baker HFGW detectors and JASON's shot-noise argument; the Pais generator patent | `gravitics::hfgw`; the fringe-physics report | covered |

### The numbers in the claim's own terms *(computed)*

The quadrupole formula gives the strain at distance r as **h ≈ 2GÏ/(c⁴r)** and the radiated
power as **P = (G/5c⁵)⟨I⃛ᵢⱼI⃛ᵢⱼ⟩**. For a burst of duration τ, I⃛ ≈ Ï/τ, so the energy is about
(G/5c⁵)(Ï/τ)²τ. The coupling is G/c⁴ = 8.26×10⁻⁴⁵ s²/(kg m) (8πG/c⁴ = 2.08×10⁻⁴³).

For NIF, take the capsule's ~20 kJ of implosion kinetic energy stagnating in ~100 ps, with a P₂
asymmetry ε, so Ï ≈ 2εE_kin:

| Source | Strain at 1 m | Energy radiated | Gravitons per shot (at 10 GHz) |
|---|---|---|---|
| NIF implosion, ε = 0.01 | 6.6×10⁻⁴² | 8.8×10⁻³⁹ J | 1.3×10⁻¹⁵ |
| NIF implosion, ε = 0.1 | 6.6×10⁻⁴¹ | 8.8×10⁻³⁷ J | 1.3×10⁻¹³ |
| The 8.6 MJ neutron burst, 1% anisotropic (memory-like step) | 2.8×10⁻³⁹ | — | — |
| Gertsenshtein: a 4 J, 0.8 µm laser pulse through 2 m of 20 T | — | — | 2.1×10⁻¹⁵ |

The best detector near these frequencies, the CAPP-12T haloscope recast at 5.3 GHz, reaches
h₀ ≈ 3.9×10⁻²¹, which is **about 20 orders of magnitude** above the lab strain at 1 m. A
Gertsenshtein probability of 1.3×10⁻³⁴ in a SPARC-like field means one graviton per ~10¹⁵ laser
pulses. The physics is standard and the effects are real; their size is what the calculation
shows.

### Verdicts and what goes into the crates

| Piece | Verdict | Target |
|---|---|---|
| Quadrupole strain `h = 2GÏ/(c⁴r)` and burst energy `(G/5c⁵)(Ï/τ)²τ` for any lab source | **new** | `gravitational::laboratory::{quadrupole_strain, burst_energy}` |
| The Einstein coupling 8πG/c⁴ = 2.077×10⁻⁴³ | **new** (a derived constant) | `metrology::EINSTEIN_KAPPA` |
| NIF, Z-pinch and SPARC estimates as worked numbers | **new** | example `gravitational/examples/lab_hfgw.rs`, using `plasmoid`'s shot data and `hf::REGISTRY` for the detector gap |
| Laser-pulse conversion in a strong magnet | **covered** (`gertsenshtein_probability`, `graviton_count`) | the same example |
| Euler–Heisenberg vacuum stress-energy | **covered** | — |
| Magnon-based HFGW limits (Ito et al. 2020, and later) | **to source**, then **new** | `gravitational::hf::REGISTRY` |
| The implosion's asymmetry ε and stagnation time for a named NIF shot | **to source** (published P₂ and bang-time widths) | the example's inputs |

## Item 3. Quantum gradiometry and mass-anomaly detection

Reviewed at the level of the general published physics, as the roadmap sets for items with
dual-use framing: the gravity and gradient of buried bodies, why gradiometers beat the
equivalence principle, sensor averaging, and how the deflection of the vertical enters inertial
navigation.

### The claim

An accelerometer cannot tell vehicle acceleration from a horizontal gravity anomaly (the
equivalence principle), so an uncompensated deflection of the vertical (DOV) of a few
arcseconds causes inertial-navigation errors (the batch says 3 km). A gravity gradiometer
measures the tidal tensor Γᵢⱼ = ∂²V/∂xᵢ∂xⱼ (1 E = 10⁻⁹ s⁻²), in which a rigid vehicle's own
acceleration cancels, leaving the curvature `R₀ᵢ₀ⱼ`. Gravity cannot be shielded, and cold-atom
sensors (10⁻¹⁰–10⁻¹² m s⁻²/√Hz, gradiometry below 0.01 E/√Hz) cross a threshold for detecting
buried voids and dense masses at a distance. The batch then proposes specific military
detection scenarios and a gravity-plus-magnetic fusion pipeline.

### What the crates already do

| Piece | Where | Status |
|---|---|---|
| Potential and all second partials of a field model, stable to degree 2700; the radial gradient on a map | `geopotential::Coefficients::{second_partials, grid}` with `Functional::RadialGradient` | covered (the full local tensor is new; see below) |
| Gradiometer mission errors per coefficient (GOCE-class) | `gravimetry::gradiometry` | covered |
| Cold-atom accelerometer and gradiometer noise | `gravimetry::atom::AtomInterferometer` | covered |
| A gradiometer read as a strain detector (`½ω²hL`) | `gravitational::detectors::{tidal_acceleration, gradiometer_strain}` | covered |
| Gravitational shielding: Majorana's coefficient against lunar laser ranging (< 10⁻²² m²/kg) | `gravitics::anomalies::{shielding_loss, LLR_ABSORPTION_BOUND}` | covered: "gravity cannot be shielded" holds to 10⁻¹² for the whole Earth |
| The Exail AQG, Stray et al.'s tunnel (SNR 8), the INS DOV papers | `research/topics/magnetics-and-mhd/notes/practical_applications.md` | covered (notes) |

### The numbers *(computed)*

**Buried bodies.** A horizontal cylinder of radius a with its axis at depth z and density
contrast Δρ (a tunnel or pipe; line mass λ = Δρπa²) gives, at horizontal offset x,

  **Δg_z = 2Gλz/(x² + z²)**, **Γ_zz = 2Gλ(z² − x²)/(x² + z²)²**, with Γ_xx = −Γ_zz (Laplace in 2-D).

A sphere of mass contrast ΔM gives Δg_z = GΔMz/(x² + z²)^(3/2) and
Γ_zz = GΔM(2z² − x²)/(x² + z²)^(5/2).

| Body | Δg_z above it | Γ_zz above it |
|---|---|---|
| Air-filled tunnel, a = 1 m, axis at 1.9 m, Δρ = −2000 kg/m³ | −44 µGal | −232 E |
| Void sphere, a = 5 m, centre at 20 m, Δρ = −2700 kg/m³ | −23.6 µGal | −23.6 E |

These are what civil engineering surveys (Stray et al. 2022) and volcano and reservoir monitoring
measure.

**Sensor averaging.** White noise of amplitude spectral density S averages down as S/√t. The
AQG's 500 nm s⁻²/√Hz gives 53 nm/s² in 90 s, 20 nm/s² in 10 min and 10 nm/s² in 40 min,
matching its datasheet. Reaching 2 µGal (20 nm/s²) takes **625 s**.
**Correction to our notes:** `practical_applications.md` says "2 µGal at about 6,250 s"; the
right figure is 625 s.

**Why gradiometers sidestep the equivalence principle.** Two proof masses a baseline L apart
share the vehicle's acceleration and differ only by ΓL: 1 E over 0.5 m is a differential
acceleration of 5×10⁻¹⁰ m/s².

**The deflection of the vertical in inertial navigation.** A deflection θ makes the navigator
read a horizontal acceleration gθ: 1″ is 47.5 µm/s² (4.75 mGal). Schuler tuning (period
`2π√(R/g)` = 84.4 min) bounds the resulting position error at about Rθ, peaking at 2Rθ:

| θ | Horizontal gravity error | Schuler-bounded peak position error |
|---|---|---|
| 1″ | 47.5 µm/s² | 62 m |
| 5″ | 238 µm/s² | 309 m |
| 10″ | 475 µm/s² | 618 m |

**Check on the claim:** "a few arcseconds causes 3 km of position error" is not what a
Schuler-bounded error gives (hundreds of metres). The 3 km in our notes comes from one simulated
route over 44 hours (MATEC 2018), where errors also accumulate through other channels.
**To source:** read that paper before quoting 3 km.

### Verdicts and what goes into the crates

| Piece | Verdict | Target |
|---|---|---|
| The full gravity-gradient tensor in the local frame (north, east, up), in s⁻² and Eötvös | **new** | `geopotential::SecondPartials::tensor_neu` |
| The deflection of the vertical from a model against a normal field | **new** | `geopotential` (needs a normal-gravity reference; GRS80 zonals are already used in `lab`) |
| Local forward models: sphere and horizontal cylinder (Δg_z, Γ_zz, Γ_xx); the rectangular prism (Nagy 1966) | **new** (prism **to source**) | `gravimetry::local` |
| White-noise averaging S/√t and time-to-precision | **new** (one line) | `gravimetry::atom` or an example |
| Schuler period and the DOV-to-position-error bound | **new** | `geopotential` example (navigation) |
| Correction: AQG time to 2 µGal is 625 s, not 6,250 s | **fix** | `practical_applications.md` |
| "A few arcseconds causes 3 km" | **to source** (MATEC 2018) | — |
| Specific detection scenarios (submarines, hardened sites, sensor fusion) | **out** as applications; the general forward models and noise formulas above contain all of the physics | — |

## Item 4. Prompt elastogravity signals and ionospheric MHD coupling

### The claim

A large earthquake or explosion redistributes mass, and the gravity change propagates at the
speed of light — a prompt elastogravity signal (PEGS) arriving minutes before the P-wave. Its
size is `Δg ∼ GΔM d/r³` (quadrupolar), at the nGal–µGal level, right at cold-atom sensor noise
floors, so it could give speed-of-light early warning. Separately, the Hunga Tonga Lamb wave,
amplifying as it climbs into thinning air, drives an ionospheric dynamo `J = σ(v × B)` and
launches Alfvén waves to the conjugate hemisphere — a geomagnetic precursor.

### What the crates already do

| Piece | Where | Status |
|---|---|---|
| Coseismic gravity change of a fault (Okada/Okubo half-space) | `gravimetry::coseismic` | covered (the static change, not the prompt transient) |
| The atmosphere's acoustic-gravity dispersion, Lamb wave, the mass of a pressure change, the ocean's forced response | `gravimetry::atmosphere` | covered |
| Cold-atom sensor noise floors | `gravimetry::atom` | covered |
| The geomagnetic field vector at 100–300 km (for v × B) | `magnetics::main_field` | covered |
| F-region Alfvén speed, the ionospheric conductivity scales | `magnetics::mhd::alfven_speed` | covered |
| Tohoku's moment, Hunga Tonga's 780 Pa / 318 m/s Lamb wave | `research/topics/.../practical_applications.md`, and the T1 pipeline in `lab` | covered (data) |

### The numbers *(computed)*

**PEGS timing and size.** At 1000 km, light takes 3.3 ms, the P-wave ~143 s and the Lamb wave
~3145 s, so a light-speed signal **leads the P-wave by about 143 s**. The monopole scale
`Δg = GΔM/r²` gives 0.67 nm/s² for ΔM = 10¹³ kg at 1000 km, 6.7 nm/s² for 10¹⁴ kg. The real
coseismic signal is quadrupolar and smaller; published PEGS (Vallée et al. 2017, *Science*;
Juhel et al. 2018) are ~0.1–1 nm/s² — at the edge of the best cold-atom floors. **To source:**
the Vallée/Juhel amplitudes and the Zhang et al. 2020 detectability study for the exact
template, rather than the monopole toy here.

**Acoustic-gravity amplification.** Energy-flux conservation `½ρv² c = const` with scale height
H ≈ 7 km gives `v ∝ exp(z/2H)`: a factor 1.3×10³ by 100 km. Extrapolated to 300 km it gives
2×10⁹, which is **unphysical** — molecular viscosity and thermal conduction dissipate the wave
first. A calculator must cap this, which is the same lesson as the mantle-viscosity guard in
`gravitational::propagation::viscous_damping`.

**Ionospheric dynamo.** The Lamb wave's surface particle velocity is `v₀ = p/(ρ₀c)` = 2.0 m/s
for 780 Pa. Amplified to the dynamo region and crossed with a 30,000 nT field, the dynamo field
`uB` is of order tens of mV/m (77 mV/m for an uncapped 100 km growth — an upper bound, since the
same dissipation applies). The F-region Alfvén speed `B/√(μ₀ρ)` for 10¹² O⁺/m³ is **164 km/s**,
so Alfvén waves reach the conjugate point in tens of seconds.

### Verdicts and what goes into the crates

| Piece | Verdict | Target |
|---|---|---|
| PEGS arrival-time gap `r/v_seismic − r/c` and the monopole/quadrupole `Δg(r)` scale | **new** | `gravimetry::coseismic` (a prompt-signal helper) |
| The light-speed lead over P and Lamb waves | **new** (one line) | the same |
| AGW velocity amplification `exp(z/2H)` **with a dissipation cap** | **new**, capped | `gravimetry::atmosphere` |
| Ionospheric dynamo field `uB` and the F-region Alfvén transit | **new** | `magnetics::space_weather` (dynamo EMF) + existing `alfven_speed` |
| The PEGS template amplitude for Tohoku | **to source** (Vallée 2017, Juhel 2018, Zhang 2020) | the example's inputs |
| Early-warning operational framing | kept to the physics (timing, amplitude, noise floor); the warning-system design is **out** | — |

## Item 5. Condensates, gravitomagnetism and atom-interferometer phases

*Not reviewed yet.*

## Item 5. Condensates, gravitomagnetism and atom-interferometer phases

### The claim

The 1990s superconductor-antigravity claims (Podkletnov, Tajmar) failed because the
gravitomagnetic coupling is suppressed by `Gm_e²/(ħc) ∼ 10⁻⁴⁵`. But the real frontier is that a
macroscopic quantum condensate (Cooper pairs, a BEC) is an ultra-coherent phase amplifier: the
DeWitt–London Cooper-pair momentum includes a gravitomagnetic term `2m_e A_g`. Atom
interferometers like MAGIS-100 measure the phase `Δφ = (1/ħ)∮p_μ dx^μ`, sensitive to ultralight
scalar/vector dark matter in the 0.01–10 Hz LISA–LIGO mid-band and to gravitational
Aharonov–Bohm and frame-dragging phases.

### What the crates already do

| Piece | Where | Status |
|---|---|---|
| The gravitomagnetic suppression `Gm_e²/(ħc)` and permeability `4πG/c²` | `gravitics::gem::GM_PERMEABILITY` (and the Podkletnov/Tajmar claim-vs-bound calculators) | covered |
| Frame dragging (Lense–Thirring) and its `Field` | `gravitics::gem::GravitomagneticDipole` | covered |
| Atom-interferometer GW detectors (MAGIS, AION) strain response | `gravitational::detectors::AtomGwDetector` | covered |
| The atom's Aharonov–Bohm / single-diffraction phase, Rosi's k_eff T² | the gravity-anomalies report and its checks (Rosi `k_eff T²` = 4.12×10⁵ rad per m/s²; Overstreet `mGM/(ħr)` = 1.63 rad/s) | covered (report), **new** as functions |
| London moment `B = 2mω/e` | `gravitics::superconductors::london_moment` | covered |
| Ultralight-dark-matter oscillation of constants | — | **new** |

### The numbers *(computed)*

- `Gm_e²/(ħc)` = 1.75×10⁻⁴⁵: this is why a spun superconductor's weight change is unmeasurable
  (as the Hathaway null showed), confirming the batch's own explanation.
- **Mid-band mapping.** An ultralight boson of mass m oscillates fundamental "constants" at its
  Compton frequency `f = mc²/h`. The LISA–LIGO gap maps to boson masses:
  f = 0.01 Hz ↔ 4.1×10⁻¹⁷ eV, 1 Hz ↔ 4.1×10⁻¹⁵ eV, 10 Hz ↔ 4.1×10⁻¹⁴ eV.
- **Atom-gravimeter phase.** `Δφ = k_eff g T²` (Kasevich–Chu). Rosi's k_eff T² = 4.12×10⁵ rad
  per m/s², so g₀ gives 4.0×10⁶ rad. (Matches the gravity-anomalies report.)
- **Gravitational Aharonov–Bohm.** `mGM/(ħr)` for ⁸⁷Rb near a 1.25 kg mass at 7 cm is 1.63 rad/s
  (Overstreet et al. 2022). (Matches the report.)

### Verdicts and what goes into the crates

| Piece | Verdict | Target |
|---|---|---|
| `Gm_e²/(ħc)` as the gravitomagnetic suppression | **new** (one constant) | `gravitics::gem` |
| Compton frequency `f = mc²/h` of an ultralight boson, and the mass↔frequency map | **new** | `gravitational` (a dark-matter mid-band helper) |
| Atom-gravimeter phase `k_eff g T²` and the Aharonov–Bohm phase rate `mGM/(ħr)` | **new** as functions (now only in the report) | `gravimetry::atom` |
| DeWitt–London Cooper-pair gravitomagnetic momentum | **to source** (DeWitt 1966, PRL 16, 1092) then a doc note; the measurable effect is the suppressed London-moment shift | `gravitics::gem` doc |
| "Condensate as a phase amplifier" as an experiment design | the measurable quantities (phase sensitivity, mid-band masses) are covered above; the device concept is **out** | — |

---

## Summary: what batch 1 adds

Ranked by how directly each becomes a tested calculator.

1. **`gravitics::anomalies::tension`** and the same-site G tests (HUST 339 ppm, BIPM 0.16σ): the
   data already shows the scatter is per-apparatus, not per-site.
2. **`geopotential::SecondPartials::tensor_neu`** (the local gravity-gradient tensor) and
   **`gravimetry::local`** (sphere, cylinder, prism forward models): foundational and widely
   useful.
3. **`gravitational::laboratory::{quadrupole_strain, burst_energy}`** and
   **`metrology::EINSTEIN_KAPPA`** (8πG/c⁴): the lab-HFGW machinery, with the honest ~20-order
   gap to detectors.
4. **`magnetics::space_weather`** dynamo EMF and **`gravimetry::atmosphere`** capped AGW growth;
   **`gravimetry::coseismic`** PEGS timing.
5. **`gravimetry::atom`** phase and Aharonov–Bohm helpers; **`gravitational`** ultralight-DM
   mid-band; the `magnetics::magnets::susceptibility_force` tie to the G bias.
6. **Scalar-field screening** (chameleon/symmetron) in `gravitics::modified`: the one genuinely
   new physics item, pending primary sources.

Two errors in our own notes were found and fixed while reviewing: the AQG averaging time
(2 µGal at 625 s, not 6,250 s) and the "3 km navigation error" (a single simulated 44 h route,
not a general consequence; flagged to-source).

All "new" items are standard, published physics. Nothing here supports the stronger
frontier readings: the G scatter is an apparatus systematic, lab HFGW strains are ~20 orders
below detection, and the gravitomagnetic couplings are suppressed by 10⁴⁵ exactly as the batch
itself notes.
