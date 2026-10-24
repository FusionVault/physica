# Review: BAASS, Peratt and Pais (batch of 2026-10-24)

The batch is [`2026-10-24-baass-peratt-pais-video.md`](2026-10-24-baass-peratt-pais-video.md),
a summary of a livestream mapped onto the crates. Reviewed as the roadmap sets out: each claim
computed in its own terms, beside the standard prediction and the measurements, with a verdict on
what goes into which crate, and the fringe pieces paired with the conventional explanation and
the bound they must beat. Numbers marked *(computed)* come from
[`2026-10-24-baass-peratt-pais-video-checks.py`](2026-10-24-baass-peratt-pais-video-checks.py).

Verdicts: **covered**, **new**, **to source**, **out** (with the reason).

---

## Pillar 1. Peratt's Z-pinch plasmoids, scale invariance and galactic rotation

### The claim

Pulsed megajoule currents through plasma undergo m = 0 and m = 1 MHD instabilities and
self-organise into a stack of ~9 toroidal plasmoids firing axial jets. Maxwell + MHD is
scale-invariant, verified over 14 orders of magnitude and extrapolated by Alfvén over 28, to the
Hubble distance. Peratt argues "dark matter" is unnecessary: intergalactic Birkeland currents
attract filamentarily (∝ 1/r, against gravity's 1/r²) and flatten galaxy rotation curves.

### What the crates already do

| Piece | Where | Status |
|---|---|---|
| Bennett equilibrium, Pease–Braginskii, Lee plasma-focus scalings, the Ostriker self-gravitating filament | `plasmoid::pinch` | covered |
| Taylor relaxation (`∇×B = λB`), spheromak eigenvalues, FRC β and s, reconnection | `plasmoid::relaxation` | covered |
| m = 0 and m = 1 stability (Kadomtsev) on the Bennett profile | `magnetics::stability` | covered |
| Alfvén speed, Lundquist number, magnetic Reynolds number | `magnetics::mhd`, `plasmoid::relaxation::lundquist` | covered |
| FuZE, MagLIF, C-2W parameters; the MOND exclusion at the wide-binary scale | `plasmoid` (notes); `gravitics` | covered |
| Galactic field strength, equipartition field, synchrotron | `magnetics::astro` (planned; see batch-1 note) | partly |

### The numbers in the claim's own terms *(computed)*

- **Scale invariance.** The Alfvén speed `B/√(μ₀ρ)` is 93 km/s for FuZE (3 T, 4.9×10²³ m⁻³) and
  65 km/s for a galactic filament (3 µG, 10⁴ m⁻³) — similar across 23 orders of magnitude in
  density. Lab (10⁻³ m) to Hubble (10²⁶ m) is **29 orders**; Alfvén's "28" uses 10⁻² m to 10²⁶ m.
  The scale-invariance of ideal MHD is standard and real.
- **The dark-matter argument, computed.** Two 1 MA Birkeland currents a parsec apart attract at
  `μ₀I₁I₂/(2πr)` = 6.5×10⁻¹² N/m. Two gas filaments of 1 M☉/pc line mass at the same separation
  attract gravitationally at `2Gλ₁λ₂/r` = 18 N/m. The magnetic force is **3.6×10⁻¹³ of gravity**.
  Even at a mega-amp, intergalactic currents are ~12 orders too weak to replace the dark-matter
  mass (9.3×10¹⁰ M☉ enclosed for a flat 200 km/s curve at 10 kpc *(computed)*). This is the
  honest result of the claim in its own terms, and it is why the mainstream does not use it. The
  crates already carry the counter-evidence: `gravitics` records the 16σ wide-binary exclusion of
  MOND, and the Pease–Braginskii current caps how much a pinch can carry.

### Verdicts

| Piece | Verdict | Target |
|---|---|---|
| MHD scale-invariance demonstrated across lab and astrophysical parameters (Alfvén speed, Lundquist, Reynolds at both ends) | **new** as an example | `magnetics` or `physica` example: "one MHD across 23 orders" |
| Ampère line-current force `μ₀I₁I₂/(2πr)` beside the Ostriker gravitational line force | **new** | `magnetics::mhd` (or beside `plasmoid::pinch::magnetic_over_gravitational`), as a labelled galactic-scale comparison |
| The 9-plasmoid / toroid stack as the saturated Taylor state | **to source** (Peratt 1986, *IEEE Trans. Plasma Sci.*); the stability and relaxation pieces exist, the specific 9-bead count needs the paper | `plasmoid` doc note |
| "Birkeland currents replace dark matter" | computed as **3.6×10⁻¹³ of gravity** at galactic scale — set in the crate as a claim-vs-bound calculator, labelled refuted at that scale, beside the wide-binary exclusion | `gravitics` or `plasmoid` |
| Peratt's underground-test diagnostics as the data source | **out** (no public numbers; not needed for the physics) | — |

## Pillar 2. BAASS, the DIRDs, the Gordon metric and cloaking

Reviewed at the level of the general physics (the Gordon optical metric, plasma refraction); the
document provenance and the radar-evasion application are out.

### The claim

A DIA document set (BAASS 10-month report, 38 DIRDs by Puthoff, Davis et al.) pairs
spacetime-metric engineering (wormholes, warp metrics, vacuum energy) with a simulation of a 10 m
sphere in a 1 m graded dielectric shell reducing its radar cross-section. The physics link is the
Gordon optical metric: curved spacetime `g_μν` is isomorphic to flat space filled with a medium
of `n(r) = √(ε_r μ_r) ≡ 1/g₀₀`, so a plasma or metamaterial shell with a graded index refracts
radar around the hull (transformation optics).

### What the crates already do

| Piece | Where | Status |
|---|---|---|
| Plasma frequency `(1/2π)√(n_e e²/ε₀m_e)` | `gravitational::hf::plasma_frequency` | covered |
| Puthoff's polarizable-vacuum metric and its light bending `K = e^(2GM/rc²)` against GR | `gravitics::spacetime::{Metric::PolarizableVacuum, light_bending}` | covered |
| Warp-bubble energy, the quantum inequality, Casimir, wormhole tensions (the DIRD topics) | `gravitics::spacetime` | covered |

### The numbers *(computed)*

- **Plasma cutoff.** `f_p = 8.98 GHz × √(n_e/10¹⁸ m⁻³)`. The critical density for 10 GHz (X-band)
  is 1.24×10¹⁸ m⁻³: below it the plasma refracts the wave, above it reflects.
- **Underdense refraction.** The index `n = √(1 − (f_p/f)²)` at 10 GHz is 0.96 at 10¹⁷ m⁻³, 0.77
  at 5×10¹⁷, 0.44 at 10¹⁸ — a usable graded index just below the critical density.
- **Gordon metric.** A ray turns through ~90° over a shell of thickness L when `dn/dr ∼ 1/L`, an
  order-unity gradient for a 1 m shell. This is the transformation-optics statement; the
  polarizable-vacuum `light_bending` is the same `n ↔ metric` identity for gravity.

### Verdicts

| Piece | Verdict | Target |
|---|---|---|
| Plasma frequency, critical density, underdense index `n(f)` | **covered** (`plasma_frequency`); add the index and critical density | `gravitational::hf` or a new `plasma` helper |
| The Gordon metric `n = 1/√g₀₀` beside the polarizable vacuum | **new** (a doc note tying them together) | `gravitics::spacetime` doc |
| The DIRD metric-engineering topics | **covered** (warp, wormhole, Casimir already in `spacetime`) | — |
| Radar-cross-section reduction / cloaking as an application | **out**; the plasma index and ray-bending physics above are the general content | — |
| The BAASS/DIRD documents themselves | **out** (provenance, not physics) | — |

## Pillar 3. The Chinese "Pais effect" patent and the Zhang high-voltage setups

This is the batch's strongest candidate for a crate calculator: a concrete claim (pulsed
high-voltage fields make a "vortex gravitational field") with two experimental setups, exactly
the claim-against-bound form `gravitics` already uses for Podkletnov, Tajmar and the EmDrive.

### The claim

A 2025 Chinese patent and Zhang Jikuan's replications, echoing Pais's US Navy patents, say that
rapidly varying high-voltage fields (∂B/∂t, high dV/dt) make a local gravito-electric or "vortex"
field. Setup 1: a 0.15 mm plastic sheet between HV wires at 40 kV deflects toward the positive
terminal and, reportedly, keeps deflecting the same way when the polarity is reversed. Setup 2: a
polyethylene ball on a thread in a vacuum Faraday cage between offset coils at 20 kV feels a
rotational torque. Gary Stephenson is cited for a 20 kV nonlinearity threshold.

### The conventional forces, computed *(computed)*

**Setup 1 is dielectrophoresis.** A neutral dielectric of volume V and permittivity ε_r in a
non-uniform field feels the Kelvin force

  **F_DEP = ½ V ε₀(ε_r − 1) ∇(E²)**.

For a 2 cm × 2 cm × 0.15 mm sheet (ε_r = 2.2) near a wire at E ≈ 10⁷ V/m with ∇E² ≈ 2×10¹⁶ V²/m³,
this is **6.4 mN**, the weight of ~650 mg — far larger than any claimed gravitational effect.
Decisively, F_DEP depends on **E²**, so reversing the polarity (E → −E) leaves it unchanged: the
sheet deflects toward the high-gradient (positive-tip) region **regardless of polarity**, which
is exactly the reported behaviour. The claim's own signature is the classical giveaway.

**Setup 2 is field angular momentum.** Crossed static E and B fields store angular momentum
`L_EM = ε₀∫ r × (E × B) dV` (the Feynman disk paradox). Switching the field off transfers that
torque to the apparatus. For E ≈ 10⁶ V/m, B ≈ 10⁻³ T, lever arm 0.05 m, volume 10⁻⁴ m³, that is
~4×10⁻¹⁴ kg m²/s of stored angular momentum delivered as an impulse — a real torque on a
suspended body with no new physics. Residual electrostatics on a "neutral" ball (patch charge,
induced dipole in the offset field) add more.

> **Correction (2026-10-24), after reading the primary sources.** See
> [the review of Zhang's sources](2026-10-24-zhang-unified-field-theory-review.md#7-corrections-to-the-batch-2-review).
> Zhang's own accounts say Setup 1's sheet moves toward the positive wire, and does so again after
> the leads are swapped, so it changes direction. Dielectrophoresis still explains it, through
> which side is live rather than through polarity independence. Setup 2's open coils pass about
> 10 µA, about 7×10⁻¹¹ T, so field angular momentum is negligible there. The rotation budget is
> the polarisation torque, `gravitics::electrogravity::rotation_torque`.

### The discriminator the crate should encode

A genuine gravito-electric field `E_g` gives `a = E_g`, **independent of material**. The
classical forces do not:

- F_DEP scales with `V(ε_r − 1)` — volume and permittivity, independent of mass density ρ.
- The Feynman torque scales with the stored field momentum, independent of the ball's mass.

So the test is **mass-independence**: swap the polyethylene ball (ρ = 950, ε_r = 2.2) for a
tungsten one (ρ = 19,300) of identical outer geometry and surface conductance. A real g-field
gives the same acceleration; a DEP/Feynman artefact gives a force that barely changes while the
acceleration drops by the density ratio (~20×). Equivalently, a shielded cold-atom gradiometer
(1×10⁻¹⁰ m/s²) outside the cage would see a true metric perturbation as an Eötvös gradient that
penetrates mu-metal, while an electrostatic force leaves none.

### Verdicts

| Piece | Verdict | Target |
|---|---|---|
| Dielectrophoretic force `F_DEP = ½Vε₀(ε_r−1)∇E²`, with its polarity-independence | **new** | `gravitics` (a new `electrogravity` module for the claim-vs-bound calculators) |
| Field angular momentum `L_EM = ε₀∫ r×(E×B)dV` and the switch-off torque | **new** | same module |
| The mass-independence discriminator (DEP/Feynman ∝ volume; gravity ∝ ρV) | **new** — this is the calculator's point | same module |
| The claim as a labelled entry (its own model vs the classical force vs what would prove it) | **new**, status contested/unverified, in the `gravitics` style ("a parallel track, not wished away") | the proposals table |
| Pais's and the Chinese patent's own field formula | **to source** (US 10,144,532 and the 2025 patent) before coding the patent's own numbers | — |
| Stephenson's 20 kV threshold | **out** (a stated voltage, no formula to verify) | — |
| Dual-use ("moving target object") framing | **out**; the physics above is the content | — |

## Pillar 4. Bearden's 4D conservation and Helion's inductive recovery

### The claim

Bearden argues standard EE forces energy conservation in 3D (killing the dipole in a closed
loop), while GR conserves energy-momentum in 4D (`∇_μ T^μν = 0`), so energy can be drawn from the
vacuum. Helion's grant for pulsed-power semiconductors is cited as real pulsed-power engineering:
colliding FRC plasmoids compressed by a pulsed field, with the expanding plasma inducing back-EMF
into the coils (`V = −N dΦ/dt`) and recharging the capacitors through solid-state switches, no
steam turbine.

### What the crates already do

| Piece | Where | Status |
|---|---|---|
| Direct conversion of charged fusion products (N-stage collector) | `plasmoid::conversion::collector_efficiency` | covered |
| Engineering gain of magnetic and inertial plants | `plasmoid::conversion::{q_eng_magnetic, q_eng_inertial}` | covered |
| Pulsed inductive breakeven | `plasmoid::conversion` | covered |
| Bosch–Hale reactivities for D–T vs D–³He | `plasmoid::fusion::bosch_hale` | covered |
| Covariant conservation and the Einstein coupling `8πG/c⁴` | — (batch 1, item 2 proposes `metrology::EINSTEIN_KAPPA`) | partly |

### The numbers *(computed)*

- **The covariant conservation law is real and is not a free-energy source.** `∇_μ T^μν = 0`
  does let the metric exchange energy with matter (cosmological redshift, Hawking radiation, the
  dynamical Casimir effect are the textbook cases). But the coupling is `8πG/c⁴ = 2.08×10⁻⁴³`
  s²/(kg m): the "stiffness" is `1/κ = 4.8×10⁴² N`. To extract lab-scale energy from a
  time-varying metric needs energy densities ~43 orders beyond anything pulsed power reaches.
  This is the honest version of Bearden's "4D inflow", and it refutes the free-energy reading.
- **Helion's inductive recovery is ordinary engineering.** With round-trip efficiency η_rec, the
  electrical return fraction is η_rec × (charged-particle fraction): 0.95 for an aneutronic,
  all-charged reaction, but only ~0.65 for D–T (80% of its energy is in neutrons). This is why
  direct recovery suits D–³He, and it is already `plasmoid::conversion`'s domain.

### Verdicts

| Piece | Verdict | Target |
|---|---|---|
| `8πG/c⁴` and the metric "stiffness" `1/κ`, as the honest form of "4D energy" | **new** (a constant) | `metrology::EINSTEIN_KAPPA` (shared with batch 1 item 2) |
| Covariant conservation as a labelled claim-vs-bound entry (energy exchange is real, lab extraction is 10⁴³ off) | **new** | `gravitics` doc, beside the free-energy device checks |
| Inductive recovery `η_rec × charged-fraction`, D–T vs D–³He | **covered** / **new** (the explicit charged-fraction weighting) | `plasmoid::conversion` |
| The Helion grant, dollar figures, semiconductor claims | **out** (business, not physics) | — |

---

## Summary: what batch 2 adds

Ranked by how directly each becomes a tested calculator.

1. **`gravitics::electrogravity`** (new module): dielectrophoretic force `½Vε₀(ε_r−1)∇E²`, field
   angular momentum `ε₀∫r×(E×B)dV`, and the mass-independence discriminator. This turns the
   Pais/Zhang claim into the same claim-vs-bound calculator the crate already uses for
   Podkletnov, Tajmar and the EmDrive — exactly the "parallel track, not wished away" the user
   asked for.
2. **Ampère line-current force** beside the Ostriker gravitational one, as the galactic-scale
   dark-matter check (magnetic force 3.6×10⁻¹³ of gravity): `magnetics::mhd` / `plasmoid::pinch`.
3. **`metrology::EINSTEIN_KAPPA`** (8πG/c⁴) and the covariant-conservation claim-vs-bound note.
4. **Plasma index and the Gordon-metric note**: mostly covered by `plasma_frequency` and
   `light_bending`; add the underdense index `n(f)` and critical density.
5. **The MHD scale-invariance example** (one MHD across 23 orders) and the 9-plasmoid Taylor-state
   note (to source: Peratt 1986).

Everything "new" is standard physics or a claim computed honestly against its conventional
rival. The three strong frontier readings all fail their own numbers: Birkeland currents are
3.6×10⁻¹³ of gravity at galactic scale, the "vortex gravity" force is dielectrophoresis to the
sheet's own polarity signature, and "4D vacuum energy" is a real but 10⁴³-suppressed
metric–matter coupling. Each goes into the crates as a calculator so the margin is reproducible,
labelled with its status.
