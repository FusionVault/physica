# Review: Zhang Xiangqian's unified field theory and US 12,417,870

The sources are listed in [`2026-10-24-zhang-unified-field-theory.md`](2026-10-24-zhang-unified-field-theory.md).
The method is the roadmap's: compute each claim in its own terms, compute the standard prediction
for the same apparatus, set both against what was measured, and put every piece in the crates as
a labelled calculator. Numbers marked *(computed)* come from
[`2026-10-24-zhang-unified-field-theory-checks.py`](2026-10-24-zhang-unified-field-theory-checks.py)
and are pinned by tests in `gravitics::electrodes`, `gravitics::electrogravity`,
`plasmoid::discharge` and `physica/examples/gravity_patent.rs`.

Verdicts: **covered**, **new**, **to source**, **out** (with the reason).

## The short answer

The patent gives no mechanism, no number for the force, and no vacuum level. Rebuilt from its
own embodiment, leads included, and simulated with every force and torque (section 2b), the rig
turns and swings the ball by classical physics alone, if the ball carries charge. With
0.2 nC and a 0.1 nC patch, which is what handling leaves on a polymer ball, switching on rocks
it through tens of degrees and swings it a quarter of a millimetre. A neutral, round ball does
not turn at all and swings 14 µm. The "near vacuum" the patent leaves unstated decides the rest:
from about 4 to 150 Pa the gas breaks down, clamps the field in the tank and drains the ball, and
nothing moves.

The theory behind it now has its constant. Zhang's book links the gravitational and magnetic
fields by `∇×A = B/f` with `f = k_mq/(4πε₀G)`, and his papers give `k_mq` three values 45 orders
apart. With the "calibrated" one, his vortex spins the ball up at 10⁻²¹ rad/s². With either
closed form, Earth's own field would spin everything on Earth at 10²⁹–10³⁰ rad/s², so those are
ruled out by everyday experience. His G, 6.685×10⁻¹¹, sits 34σ above the highest of the 16 CODATA
measurements. Read as the acceleration of the conduction electrons, his field moves the ball by
10⁻¹⁵ m. General relativity gives the field's energy a pull of 10⁻²⁷ m/s².

Reading the primary sources also corrected two statements in the
[batch-2 review](2026-10-24-baass-peratt-pais-video-review.md) (see the last section).

## 1. The theory, equation by equation

| Equation | In standard terms | Verdict |
|---|---|---|
| `B = V×E/c²` | The magnetic field of a charge moving at `v ≪ c`. Correct. | covered (textbook) |
| Faraday's law from `B = V×E/c²` and a Lorentz transformation | Reproduces `∇×E = −∂B/∂t`, a standard result. | covered |
| `dB/dt = −A×E/c² − V×(dE/dt)/c²` | The time derivative of the first equation for one charge. `A = −dv/dt` is that charge's acceleration (his Eq. 14); the "gravitational field" is a relabelling of it. | **new**: `electrogravity::drift_velocity` gives the carriers' acceleration that bounds it |
| `E_θ/e_r = A×R/c²` (accelerating positive charge makes a reverse gravitational field) | Exactly J. J. Thomson's kink construction of the radiation field, `E = qa sin θ/(4πε₀c²r)` (Purcell, app. B): the ratio of radiation to Coulomb field is `ar sin θ/c²` *(computed: equal to 10⁻⁴ V/m in both forms for 1 nC at 10¹² m/s², 1 m)*. Correct physics, relabelled. In Maxwell's theory it is an electric field, acting only on charge and polarisation. | **new**: `electrogravity::radiation_field` |
| `∇×A = B/k` (vortex gravitational field) | The book derives `∇×A = B/f` and `E = −f ∂A/∂t` (ch. 34, 36) with `f` built from his charge and mass constants and `4πε₀G`: in the papers' notation `f = k_mq/(4πε₀G)`, the constant the 2025 air-gap paper calls `k`. A uniform `B` then spins any body up at `B/(2f)` about the field line. With the three published `k_mq`, `f` is 9.5×10⁹, 5.3×10⁻³⁵ and 4.3×10⁻³⁶ kg s/C *(computed)*. | **new**: `electrogravity::{vortex_constant, vortex_spin, ZHANG_KMQ}`; the two closed forms refuted (Earth's field would spin everything at 10²⁹–10³⁰ rad/s²), the calibrated one invisible |
| `P = (C − V)M`, `F = (C − V)dm/dt`, `q = k_mq dM/dt`, `g = dC/dt` | Definitions. `k_mq` is 7.08×10⁻¹¹ C s/kg "experimentally calibrated" (Paper J, citing a paper we have not seen), `eG/c³` = 3.97×10⁻⁵⁵ (Paper E, Eq. E42) and `eG/(4πc³)` = 3.16×10⁻⁵⁶ (Paper J's "canonical" value). | **new**: `ZHANG_KMQ` carries all three, labelled |
| `G = 4πc²/μ_g` with `μ_g` = 1.6895×10²⁸ kg/m (Papers E, J) | 6.6849×10⁻¹¹: 1 390 ppm above the highest of CODATA 2022's 16 measurements (BIPM-01, 34σ) and 70σ from the adopted value *(computed)*. The papers call it "within 0.16%"; the measurements agree with each other to a few hundred ppm. | **new**: `ZHANG_MU_G`, with the test against `anomalies::CODATA_2022_G`; refuted |
| Paper J's closures `α = (e k_mq)² μ_g/(ħc)` and `ε₀ = (4π k_mq μ_g)⁻¹` | With the paper's own numbers they give 6.9×10⁻⁵ (not 7.3×10⁻³) and 6.7×10⁻²⁰ (not 8.9×10⁻¹²) *(computed)*: the paper's own numbers do not close. | recorded in the `zhangs_constants_in_his_own_terms` test |
| `m = ω²r³/G`; `4π²r³c²/(T²Ghν) = 1` | Kepler's third law, and that law with the central mass written as its Compton frequency `Mc²/h`. The Moon's orbit gives 0.99728 *(computed)*, off by the Sun's perturbation. An identity, not a new link between the constants. | **new**: `electrogravity::five_constants`, with the test |
| Fields propagate faster than light | GW170817 bounds the speed of gravity to `−3×10⁻¹⁵ < (v − c)/c < 7×10⁻¹⁶`. | covered: `gravitational::propagation::speed_bounds`; refuted for gravitational fields |
| The "quasi-excited" craft (5000 kg at 0.45 g) | No formula. | **out**: nothing to compute |

## 2. The patent's rig, rebuilt *(computed)*

Geometry from the patent and the book: coils 12 turns of 1.5 mm copper, radius 3.5 cm, from 6.5
to 17.5 cm above and below the centre of a 10 cm tank; upper coil at 50 kV, lower grounded; ball
0.35 g of polyethylene (950 kg/m³, radius 4.447 mm, ε_r = 2.3) at the centre. The coils are
solved as 24 charged rings (`gravitics::electrodes`; its ring field agrees with direct summation
to 10⁻¹²). This first pass treats the coils alone and the ball at rest; section 2b adds the leads
and the motion.

| Quantity | Value | What it means |
|---|---|---|
| Charge on the live coil / induced on the grounded one | 278.5 nC / −57.5 nC | 5.6 pF; 7.0 mJ stored |
| Field at the ball | 2.50×10⁵ V/m, downward | single-ended: strongest toward the live coil |
| `∇(E²)` at the ball | 1.28×10¹² V²/m³ | |
| Pull toward the live coil (`sphere_force`) | 1.9 µN = 5.5×10⁻⁴ of the weight | lightens the thread; no swing, by symmetry |
| Charge if ions drift one way, as in a corona (`field_charge`) | up to 0.88 nC (Pauthenier) | the field then pulls 0.22 mN = 6.4% of the weight, a 3.7° swing; a glow instead drains the ball to a few pC (2b) |
| Polarisation-torque ceiling (`rotation_torque`) | 2.1×10⁻⁷ N m | a thread with a 30 s free twist (`torsion_constant`, 1.2×10⁻¹⁰ N m/rad) would be wound 279 turns; 1% of the ceiling gives 2.8 turns |
| B at the ball from 10 µA (`magnetics::magnets::solenoid_bz`) | 7.4×10⁻¹¹ T | 1.5×10⁻⁶ of Earth's field |
| Field angular momentum over the tank | 1.3×10⁻²⁰ kg m²/s | would spin the ball at 5×10⁻¹² rad/s, one turn in 40 000 years |
| The theory in its own terms (`drift_velocity`) | electrons drift at 4.2×10⁻¹⁰ m/s at 10 µA | at the 30 kHz ripple, a displacement of 2×10⁻¹⁵ m |
| General relativity (`field_energy_gravity`) | 7.0 mJ is 7.7×10⁻²⁰ kg | a pull of 1.4×10⁻²⁷ m/s² at 6 cm |

A static axially symmetric field puts no torque about the thread on an isotropic ball. Turning
needs symmetry broken: the staggered open ends (FIG. 3), the helical winding, a seam in the ball,
patch charge, or a charged, twisted thread. Any of these at 1% of the ceiling turns the ball by
turns. The rig has no measured torque to set against the ceiling, so the verdict is that the
classical budget exceeds what is needed by orders of magnitude; which piece acts is a question
for the measurements in section 5.

## 2b. The rig with its leads, simulated *(computed)*

The axisymmetric picture above leaves out the leads, and they matter. Solved with the leads as
wire segments (`gravitics::electrodes`, 93 pieces; an independent standard-library solver in the
checks script agrees to six digits), the field at the ball is 2.67×10⁵ V/m, tilted 12.2° off
vertical: a sideways part of 5.6×10⁴ V/m. `physica/examples/gravity_patent.rs` then integrates
the ball's swing and twist at 1 ms steps, with every force and torque below.
`web/gravity-patent.html` runs the same model and checks itself against the example's pinned
values on load.

At 50 kV and 1 Pa, 3 s after switching on, with 0.2 nC on the ball, a 0.1 nC patch and 1% out of
round:

| Force on the ball | Fraction of its weight |
|---|---|
| Coulomb, on its 0.2 nC | 1.6×10⁻² |
| Polarisation pull | 3.8×10⁻⁴ |
| The patch dipole in the field's gradient | 8.8×10⁻⁵ |
| Gas drag (at Kn ≈ 1.5) | 1.5×10⁻⁷ |
| Image in the glass wall | 5×10⁻⁸ |
| Buoyancy | 1.3×10⁻⁸ |
| Lorentz, in the coil's and Earth's fields | 2×10⁻¹⁵ |
| Diamagnetism | 7×10⁻¹⁷ |
| Gravity of the field's 12 mJ (general relativity) | 6×10⁻²⁹ |

The torques about the thread are 9×10⁻⁹ N m from the patch dipole, 10⁻¹¹ from the shape, and
10⁻¹⁰ per radian from the thread. So the patch turns the ball to face the leads' sideways field,
and it rocks there with a period of about 3 s.

| Case (4 s on) | Largest turn | Largest swing |
|---|---|---|
| The defaults above | 68° | 0.25 mm |
| A neutral, round ball | 0° | 14 µm |
| The PTFE ball, same charge | 55° | 0.13 mm |
| The symmetric drive | 13° | 0.21 mm |
| At 10 Pa (the gas breaks down past 574 V) | 0° | 0.1 µm |

So the reported rotation has a classical explanation that needs only charge on the ball, and it
makes a prediction: a neutralised ball should stop turning, and the turn should follow the leads'
layout, not the coils' current.

## 3. What the patent leaves out *(computed)*

This is the physics the patent needs and does not state, and what the crates now compute:

- **The vacuum level.** In a 10 cm tank, Townsend breakdown needs 857 V at 100 Pa and 307 V at
  10 Pa, against a minimum for air of 305 V at 1.1 Pa m. Below 4.1 Pa no Townsend breakdown fits
  in 10 cm (`plasmoid::discharge::paschen`). The mean free path is 6.6 mm at 1 Pa and reaches the
  tank's size below 0.07 Pa (`mean_free_path`). A rough vacuum of 5–100 Pa therefore sits near
  the Paschen minimum. There a 30 kHz ripple can light an electrodeless glow through the walls
  and charge the ball. Excluding ions with margin needs well under 1 Pa, where the mean free path
  approaches the tank's size, and that takes a turbo pump.
- **The open ends in air.** Peek's corona onset for the 0.75 mm copper is 6.3 MV/m against a
  surface field of about 12 MV/m at 50 kV, so the bare break points are in corona, as Zhang notes
  for his flake rig (`peek_onset`).
- **Which coil is live.** The lower coil is grounded, so the drive is single-ended. The
  polarisation pull points to the live coil whatever the sign. Driving the coils alone at ±25 kV
  makes `∇E²` vanish at the centre (`electrodes` test). With the leads it no longer vanishes, but
  it falls and changes sign (the example's `the_symmetric_drive` test).
- **The ball's charge.** The patent says nothing about it, and it is the biggest force in the
  rig: 0.2 nC gives 1.6% of the weight against 0.04% for the polarisation pull. A glow drains it
  to the floating potential, a few picocoulombs (`floating_potential`).
- **The thread.** Its torsion constant, from the free twist period (`torsion_constant`), sets
  how little torque is visible.

## 4. The other rigs, computed both ways

### 4.1 The flake between gapped wires (2023)

**The claim.** The flake moves toward the positive wire, and toward the positive wire again when
the leads are swapped. In the theory, the positive charges accelerating toward the negative wire
make a gravitational field pointing back at the positive one.

**The standard prediction.** A high-voltage supply ("high-voltage pack" here) usually has one
output near ground. Swapping the leads then moves the live side, and the polarisation pull
follows the live side. A net charge on the flake (triboelectric, or from the corona Zhang reports
at the joints) does the same, through `qE`. Either reproduces "toward the positive wire" without
new physics. The batch-2 review's reading, that the flake moves the same way whatever the
polarity, does not match Zhang's own account; see the corrections below.

**Verdict.** **covered** (`dielectrophoretic_force`, `field_charge`, `electrodes`), with the
live-end test as the decider.

### 4.2 The nylon ring on a balance (Zhang & Xu)

**The claim, in its own terms.** About 1 g of a 0.376 kg ring is a field of 2.6×10⁻² m/s²
(2 700 ppm of g) *(computed)*, held for as long as the voltage is on. In the theory the field is
the acceleration of the charges. A steady voltage accelerates none, so the steady part of the
prediction is zero. The 1 ms ramp of 30 µA gives 1.3×10⁻⁶ m/s² while it lasts, 2×10⁴ below the
claim.

**The standard prediction.** The two wires solved as 400 rings, the upper live at 55 kV, the
lower grounded *(computed)*. The ring is pulled up, toward the live wire, by 0.10–0.37 g
(between the normal-field factor `1 − 1/ε_r` and the dilute `ε_r − 1` for nylon, ε_r = 3.6),
or 0.20–0.72 g if the ring's inner face comes within 5 mm of the cable. With the lower wire live
it is pulled down by the same amount. With ±27.5 kV it is not pulled at all. The model's line
charge is 537 nC/m, against 511 nC/m for an isolated wire of the same length, so the model holds
together.

| | Reported | Standard (polarisation) | Zhang's theory |
|---|---|---|---|
| Upper wire positive | −0.83 to −1.20 g | −0.10 to −0.72 g | lighter, by 1.3×10⁻⁶ m/s² during the 1 ms ramp (its field is the charges' acceleration), then nothing |
| Upper wire negative | +0.76 to +1.08 g | +0.10 to +0.72 g | heavier |
| Held while the voltage is steady | yes | yes (a static field) | no (needs acceleration) |
| Seven cables against one | ×7.1 | ×4.4 (seven live cables as a 13 mm cage) | ×7 (linear superposition) |
| Same drive at ±27.5 kV | not done | zero | unchanged |

The odd part of the reported change is −0.97 g and the even part −0.05 g. A net charge on the
ring would also give an odd part: 47 nC at 2×10⁵ V/m, 5% of the surface charge that would break
down the air, so it is easily reached by handling nylon.

Standard electrostatics gets the sign pattern exactly and the size to within about a factor of
3 from a geometry the paper only partly gives. The paper does not give the length of the wires,
the distance to ground, or where the cable sits inside the ring. The seven-cable ratio is the one
number the crude cage model does not reach (×4.4 against ×7.1). It is recorded here as open, not
explained away: it needs the real layout of the seven cables, solved in 3D.

**Verdict.** **new** (`electrodes`, the example's `ring_pull`), with three tests:
1. **The live-end test.** Drive ±27.5 kV. Polarisation gives zero; a field from accelerating
   charges does not care where the ground is.
2. **The screening test.** A grounded metal box screens every electric force and no
   gravitational one. An accelerometer, or a second balance, inside a grounded box fixed at the
   ring's place must read 2.6×10⁻² m/s² if the field is gravitational. A phone's accelerometer
   resolves 10⁻³ m/s².
3. **Mass independence** (`mass_independence`). Swap the nylon ring for one of the same shape
   and different density.

### 4.3 The air gap of a 0.5 T core (Zhang & Xu, 2025)

**The standard prediction** *(computed)*. Taking `B ∂B/∂x ≈ B²/gap = 14 T²/m`, the
susceptibility force (`magnetics::magnets::susceptibility_force`) on the water-based samples
(leaf, pork skin, paper; χ ≈ −9×10⁻⁶) is 1.0–1.3% of their weight. That is a visible deflection
of about a millimetre on a 10 cm thread, and it acts with DC too. The aluminium plate at 50 Hz
(skin depth 11.6 mm, `skin_depth`) carries eddy currents whose time-averaged repulsion is about
0.08 N, four times its weight. With DC it gets the same impulse at switch-on and switch-off,
which is when the paper sees motion. Ceramics often hold paramagnetic iron at χ ~ 10⁻⁴, ten times
water. Every listed material moves, as reported.

**The claim in its own terms.** With `f = k_mq/(4πε₀G)`, 0.5 T in the gap spins a sample up at
2.6×10⁻¹¹ rad/s² with the calibrated constant, and at 5×10³³ rad/s² with Paper E's
*(computed)*. The first is invisible; the second would shred the core and everything near any
magnet. Neither matches a visible, material-dependent motion.

**Verdict.** **covered** (`susceptibility_force`, `skin_depth`).

### 4.4 The Michelson fringe shift (2025, log entry only)

At 10⁶ V/m the field presses on a conductor at `ε₀E²/2` = 4.4 Pa (`electrogravity::electric_pressure`),
4.4×10⁻⁴ N on 1 cm². That is enough to move a mirror mount by micrometres. The gravity of the
field energy changes the path by about 10⁻⁴⁶ m. **To source**: no paper, voltages or geometry.

### 4.5 Claims without data

The April 2024 mass change on a balance is probably the rig in 4.2. The 80 V/30 A version has no
published rig: **to source**. The site's three converging lines of evidence are as follows.
Pais's patents are covered in `gravitics::hfgw`, beside the Navy's 2016–19 test that did not show
the effect. The 2021 Brookhaven result (STAR, Breit–Wheeler pair production in heavy-ion
collisions) is QED, covered in `gravitics::strong_field`, and says nothing about gravity. And
this rig. None of the three is evidence for the others.

## 5. The experiment that would settle it

Everything below fits the patent's own budget (about $300–500, per the site):
1. Pump the tank below 0.1 Pa and record the pressure.
2. Swap which coil is grounded, or drive them at ±25 kV. The live-side pull falls and changes
   sign; a gravitational field made by the charges should not care where the ground is.
3. Measure the ball's charge with a Faraday cup before and after a run, then neutralise it (an
   ionising blower in air, before pumping down). The simulation says a neutral, round ball turns
   by nothing and swings 14 µm.
4. Time the free twist of the ball on its thread (`torsion_constant`) and record the twist angle
   and its timing at switch-on and switch-off.
5. Repeat with a PTFE ball of the same size: ε_r 2.1 against 2.3, density 2200 against 950.
   Polarisation forces give the polyethylene ball 2.6× the acceleration (density ratio 2.32 ×
   Clausius–Mossotti ratio 1.13); a gravitational field gives both the same.
6. For the weighing rig, put an accelerometer in a grounded box at the ring's place.

## 6. What went into the crates

| Item | Where |
|---|---|
| Electrodes of rings and straight wire segments: potential, field (elliptic integrals by the AGM, logarithms for segments) and the charges for given potentials | `gravitics::electrodes` (new) |
| Clausius–Mossotti factor, DEP on a sphere, Pauthenier charge, the polarisation-torque ceiling, a thread's torsion constant, spheroid depolarisation and the torque on a body not quite round, the image force | `gravitics::electrogravity` |
| Thomson–Purcell kink field (Zhang's `E_θ/e_r`), carrier drift (the theory's own bound), gravity of the field energy, electric pressure, the five-constant identity, Zhang's `k_mq`, `μ_g` and vortex constant | `gravitics::electrogravity` |
| Paschen's law and its minimum, air's Townsend constants, mean free path, Peek's corona onset, gas density and molecular speed, slip-corrected drag on a moving and a turning sphere, ion mobility, the floating potential | `plasmoid::discharge` (new) |
| The patent's rig simulated (geometry, every force and torque, 1 ms steps) and the weighing rig; a calibration case, so an example rather than library code | `physica/examples/gravity_patent.rs` |
| The 3D model of the same simulation, checked against the example on load | [`web/gravity-patent.html`](../../web/gravity-patent.html) |

## 7. Corrections to the batch-2 review

1. **Setup 1's polarity.** The batch said the sheet "keeps deflecting the same way when the
   polarity is reversed", and the review took that as the dielectrophoretic signature. Zhang's
   own accounts (7th edition, ch. 37; the repository's experiment log) say it moves toward the
   positive wire, and does so again after the leads are swapped. That means it changes direction.
   Dielectrophoresis still explains it, but through the live side, not through polarity
   independence. `gravitics::electrogravity`'s docs and example are corrected.
2. **Setup 2's mechanism.** The review computed field angular momentum with an assumed 1 mT.
   The patent's open coils pass about 10 µA. That gives 7×10⁻¹¹ T, and an angular momentum that
   would turn the ball at 5×10⁻¹² rad/s. The rotation budget is the polarisation torque, together
   with the charge from any glow in the tank.
3. **The symmetric drive.** The first version of this review said ±25 kV makes the pull at the
   centre vanish. That holds for the coils alone. The leads break the symmetry, so the pull falls
   and changes sign without vanishing.
