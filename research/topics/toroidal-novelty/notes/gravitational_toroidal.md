# Gravitational analogues of toroidal moments and anapoles

Every number marked "recomputed" comes from `gravitational_toroidal_checks.py` in this folder
(numpy from the user site; `python3 -P gravitational_toroidal_checks.py`; it exits 1 on any
failure, and it passes). Constants are the ones in `crates/metrology`. Downloads are in `dl_grav/`
and are untrusted. Labels: **known** (found in a source we opened, or textbook); **new to our
search** (not found, after the searches listed at the end, which is not the same as new);
**uncertain** (seen only as an abstract or snippet, or not checked).

Conventions: linearised gravity in harmonic gauge, signature (−+++),
`□h̄^μν = −16πG T^μν/c⁴`; time dependence `e^(−iωt)`, `k = ω/c`; `ρ = T⁰⁰/c²`, mass current
`Σ_i = T^0i/c`, stresses `T_ij`; conservation `ρ̇ + ∂_kΣ_k = 0`, `Σ̇_i + ∂_jT_ij = 0`.

What the crates already have: `magnetics::toroidal` (Dubovik's `T`, `p_eff = p + (ik/c)T`, ω⁶
power, the anapole, `aharonov_bohm_phase`), `gravitics::gem` (GEM conventions, `GM_PERMEABILITY`
= 4πG/c², the dipole field), `gravitational::energy::energy_flux_coefficients` (Blanchet's
`a_ℓ`, `b_ℓ`; `a₃ = 1/189`, `b₂ = 16/45` are used below), `gravimetry::atom::aharonov_bohm_rate`
(scalar gravitational AB only) and `AtomInterferometer`, and the "spin-2 factor" analogy in
`physica::analogies`.

## Verdict in one paragraph

Nothing here is a research-paper result. Three things are provable and checked, and none was
found in the literature in this form: (1) **the gravitational toroidal dipole cannot radiate**. In
the l = 1 mass moment the toroidal term is cancelled exactly by the stress term, through momentum
conservation. In EM the same term is Dubovik's `−Ṫ/c²` and radiates as ω⁶. (2) **At l = 2 the
"toroidal" correction is not a property of the mass current.** The stresses that drive the current
enter at the same order. The same closed current can be made to radiate or not by changing its
stresses alone, and a lone closed current with a nonzero toroidal dipole can radiate nothing at any
frequency. (3) **The compact gravitomagnetic Aharonov–Bohm (AB) phase** round a toroidal mass
current is `16πGmṀ(R − √(R² − a²))/(ħc²)`, 4× the naive analogue. For one tonne at 1 km/s it is
2×10⁻¹² rad, about 10¹⁸ shots of a CARIOQA-class interferometer.

Each is either known in substance (no dipole gravitational waves; only the TT stresses radiate;
Stachel–Aharonov–Bohm (Stachel–AB) holonomy round a gravitational solenoid; the Devaney–Wolf
construction) or a short consequence of it. The honest ceiling is a pedagogical note (Eur. J. Phys.
or Am. J. Phys. level), and only after reading Damour & Iyer 1991 and Dubovik & Tugushev 1990, which
we did not open. A referee is likely to find (1) and (3) there or in the textbooks.

## 1. The GEM analogue: toroidal mass-current moments and a gravitational anapole

### Takeaway
Gravitomagnetism (GEM) obeys magnetostatics with `J → ρv`. A closed poloidal mass current
therefore has a toroidal moment `T_g = (1/10)∫[(x·Σ)x − 2r²Σ] = ½∫(x·Σ)x` and no gravitomagnetic
dipole. Its exterior gravitomagnetic field is zero (ideal winding) or a contact term (point
moment). The exterior vector potential is a pure gauge with nonzero holonomy. Lynden-Bell & Katz
built the exact general-relativistic (GR) version in 2012. A non-radiating oscillating
("anapole") gravitational source exists, and the Devaney–Wolf construction makes its whole
exterior field vanish. Unlike EM, gravity needs no electric-dipole partner: a lone closed current
with a toroidal dipole can be exactly non-radiating once its stresses are chosen.

### Cited Findings
- **[known] Gravitational solenoid, exact GR.** "a current of matter moves up on the inside of a
  toroidal shell and down on the outside"; the metric is "static outside the torus but
  stationary inside", with the gravitomagnetic field confined to the torus. Their field equation
  is `curl H = −2κJ`, with κ = 8πG/c⁴ (so −16πJ for G = c = 1), against EM's `curl H = 4πj`. Their
  `H = −8I grad φ` replaces EM's `2I`: a factor of 4. They give energy conditions for the shell.
  [Lynden-Bell & Katz 2012, CQG 29, 115010, arXiv:1205.1616](https://arxiv.org/abs/1205.1616),
  eqs. 1.5, 1.7, 2.14; opened (`dl_grav/1205.1616.txt`).
- **[known] Its exterior is detectable only through the holonomy.** "the gravitational
  Stachel-Aharanov-Bohm effect will still affect the interference of light some of which passes
  through the torus and some outside" —
  [Lynden-Bell & Katz 2014, MNRAS, arXiv:1312.6043](https://arxiv.org/abs/1312.6043), §3 end;
  opened. Original idea: Stachel, PRD 26, 1281 (1982), "Globally stationary but locally static
  space-times: a gravitational analog of the Aharonov-Bohm effect" (cited there; not opened).
- **[known, EM side] A static toroidal dipole has a contact field only, and a non-radiating source
  can have localised potentials.** Already in `magnetics::toroidal` (Dubovik & Cheshkov 1974;
  Nemkov, Basharin & Fedotov). For EM the Devaney–Wolf criterion is that the transverse Fourier
  components vanish on the shell `|q| = k`:
  [Nemkov et al., arXiv:1605.09033](https://arxiv.org/abs/1605.09033) (snippet; not opened).
- **[uncertain] A GR-inspired optical toroidal AB effect** (light in toroidal metamaterials, by
  analogy with stationary metrics; no gravitational toroidal moment) —
  [Besharat, Miri & Nouri-Zonoz 2019, J. Phys. Commun. 3, 115019, arXiv:1904.00270](https://arxiv.org/abs/1904.00270),
  abstract only.
- **[uncertain] Fringe-adjacent.** Stephenson, Rieken & Bhargava (MG15 proceedings, 2022,
  pp. 1934–1939; INSPIRE abstract only) propose measuring "acceleration fields along the axis of
  symmetry of a toroidal coil" from an SMES toroid's gravitomagnetic field with a quantum
  gravimeter. In GEM the exterior field of an ideal toroid is zero, and the force on a test mass at
  rest is zero in any field (it scales with the mass's velocity). The electron mass current is
  `NI·m_e/e`, so the 24 A-turn toroid of `magnetics::toroidal` carries 1.36×10⁻¹⁰ kg/s
  (recomputed). Computed both ways: GR gives no axial acceleration. The claim conflicts with GEM as
  summarised in the abstract, so its status is **contested (abstract only)**.

### Inferences
- **[new to our search] Gravitational anapole, Devaney–Wolf form (Theorem 3 below).** Take any
  smooth compact symmetric `F_ij`. The source `T_ij = (∇² + k²)F_ij`, with `Σ` and `ρ` fixed by
  conservation, radiates nothing, and its harmonic-gauge field is `h̄^μν ∝ F^μν`, zero outside the
  source. Its Newtonian quadrupole `−2∫F_<ij>/c²` is not zero: the 1/c² corrections cancel it.
  Recomputed: P/P_ref = 8.6×10⁻³⁰, exterior field 4×10⁻¹⁷ of its scale. It is the routine
  carry-over of Devaney–Wolf to a tensor Helmholtz equation. We found no source that states it for
  gravity.
- **[new to our search] Gravity's anapole needs no dipole partner.** The EM anapole is
  `p = −(ik/c)T`. In gravity the l = 1 sector never radiates (Q2), so no partner is needed there. A
  single centred Gaussian vortex cell (a closed poloidal mass current with toroidal dipole
  `T_z = 0.108i` against a scale of 0.26 in script units) has stress `T_ij = x_(iU_j)`. That stress is
  longitudinal in Fourier space (`q_iW_j + δ_ij f` terms only), so the TT projection kills it at
  every frequency. Recomputed: power 9×10⁻²⁶ of a comparable pair. "Longitudinal stresses don't
  radiate" is textbook. The point that a lone closed toroidal current can be such a source is the
  part we did not find stated.
- **[known in substance] Static exterior.** For a stationary closed current, `∇×b_g` vanishes
  outside the sources, so a toroidal mass-current moment couples only to mass currents passing
  through it. Physically, a test body outside feels the mass monopole and the stresses' active mass
  (`ρ + T_kk/c²`), never the toroidal moment.

### Gaps
- Dubovik & Tugushev (Phys. Rep. 187, 145, 1990) may discuss gravity; not opened.
- Energy conditions for the anapole source: its oscillating part must ride on a static background
  heavy enough to keep `ρ > 0`. We did not work out how heavy.
- Papapetrou (Ann. IHP A 11, 57, 1969), "On the existence of gravitational fields which are
  stationary initially and finally", may contain exact non-radiating time-dependent solutions. The
  page had no abstract; not opened.

## 2. Gravitational-wave content: where a toroidal moment enters

### Takeaway
The toroidal-type term is the `−(4(2ℓ+1)/(c²(ℓ+1)(2ℓ+3)))∂_t∫x̂_iL σ_i` term of the 1PN mass
moments (Blanchet & Damour 1989): **known**. It is a 1/c² correction to the mass moments, the
analogue of EM's `E_lm = Q̇_lm + k²T_lm`, but with three differences, all proved and checked here.
(a) At l = 1 it is cancelled exactly by the stress part of σ, so a gravitational toroidal dipole
never radiates. In EM it radiates `μ₀ω⁴|(ik/c)T|²/(12πc)`.
(b) At l = 2, for a closed current (ρ static), the radiative moment is
`I_ij = c⁻²[∫x̂_ij T_kk − (4/3)Ȧ_ij]`, with toroidal moment `A_ij = ∫(x·Σ)x̂_ij`. The stresses
`T_kk` enter at the same order, and with the same weight in our example. The toroidal term alone
gives 16× the power, the stress term alone 9×, and together the exact value. At fixed mass current
the power rises as ω⁸; the current quadrupole's rises as ω⁶.
(c) For a circular binary the toroidal-type term contributes `−(32/21)(1 − 3ν)` to the 1PN flux
coefficient `−1247/336 − 35ν/12`. That is 9% of it for equal masses and 41% in the test-mass limit.
The split is convention-dependent.

### Cited Findings
- **[known] The 1PN mass moments.**
  `I_L = ∫x̂_L σ + (1/(2c²(2ℓ+3)))∂_t²∫r²x̂_L σ − (4(2ℓ+1)/(c²(ℓ+1)(2ℓ+3)))∂_t∫x̂_iL σ_i + O(c⁻⁴)`,
  with `σ = (T⁰⁰ + T^ii)/c²` and `σ_i = T^0i/c`. The all-order form, with a `δ_ℓ+2 x̂_ijL Σ_ij` term,
  is Theorem 7. The current moments are `J_L = ∫ε_ab<i_ℓ x̂_L−1>a σ_b + O(c⁻²)`.
  [Blanchet, Living Rev. Rel., arXiv:1310.1528](https://arxiv.org/abs/1310.1528), eqs. 135, 138,
  140; opened (`dl_grav/1310.1528.txt`). Blanchet & Damour 1989 is the original (not opened).
  Thorne 1980 (RMP 52, 299) has the same slow-motion formula (not opened).
- **[known] The 1PN circular flux** `F = (32c⁵/5G)ν²x⁵[1 + (−1247/336 − 35ν/12)x + …]` and Kepler
  `Ω² = (Gm/r³)[1 + (−3 + ν)γ + …]` — same review, opened.
- **[known] Only two families of radiative moments in linearised gravity.** Mass and current; the
  other source moments are gauge (`W_L, X_L, Y_L, Z_L`, eq. 141 of the review). Zschocke gives a
  detailed proof of the 10-moment STF theorem:
  [arXiv:1401.1056](https://arxiv.org/abs/1401.1056), opened. Damour & Iyer, PRD 43, 3259 (1991),
  "Multipole analysis for electromagnetism and linearized gravity with irreducible Cartesian
  tensors", is the primary source and the most likely place for prior art on the toroidal question.
  **Not opened.**
- **[known, EM] Toroidal moments are not an independent radiating family.** They are the
  `O(k²)` part of the exact electric moments (Fernandez-Corbaton, Nanz & Rockstuhl,
  arXiv:1507.00755, already cited in `magnetics::toroidal`).
- **[uncertain] Prior student work.** A 2018 University of Puget Sound summer project "uses a
  linearized form of Einstein's General Relativity to find the quadrupole moment from an
  oscillating toroidal mass and charge current" (search snippet; the site failed DNS, not opened).

### Inferences
- **Theorem 1 [the proof is new to our search; the content is known: no dipole GWs].** For any
  conserved compact source, the O(c⁻²) part of the l = 1 mass moment vanishes at every instant.
  *Proof.* Let `N_i = ∫ρr²x_i`, `A_i = ∫(x·Σ)x_i`, `B_i = ∫r²Σ_i`. Conservation gives
  `Ṅ_i = 2A_i + B_i`, `Ȧ_i = ∫(x_jT_ij + x_iT_jj)` and `Ḃ_i = 2∫x_jT_ij`. With
  `∫x̂_ijΣ_j = A_i − B_i/3`, eq. 138 at ℓ = 1 gives
  `c²I_i^(1) = ∫x_iT_jj + (1/10)N̈_i − (6/5)∂_t(A_i − B_i/3) = ∫x_iT_jj − Ȧ_i + ½Ḃ_i = 0.` ∎
  The same computation in EM (no stress, weight 3/10 in place of 6/5) gives
  `(1/10)N̈ − (3/10)∂_t(A − B/3) = −∂_t(A − 2B)/10 = −Ṫ_i`, Dubovik's toroidal dipole. That
  reproduces `p_eff = p + (ik/c)T` of `magnetics::toroidal::effective_dipole`. The EM weight 3/10
  is inferred: it reproduces both coefficients of Dubovik's T from one parameter. We did not take
  it from Damour & Iyer. Recomputed on a random conserved source: gravity residual 8.5×10⁻¹⁵ of the
  term sizes. Without `T_jj` in σ the residual is 0.37. EM matches `iωT_i` to 10⁻⁸. So **the
  gravitational toroidal dipole radiates nothing**: momentum conservation ties its rate of change to
  stress moments that cancel it. EM's charge current has no such constraint.
- **Theorem 2 [new to our search as a statement; it follows directly from eq. 138].** If ρ is
  static and Σ closed, then `∫x̂_ijkΣ_k = (7/5)A_ij`, and by momentum conservation
  `Ȧ_ij = ∫x̂_ijT_kk + 2S_ij` with `S_ij = ∫x_kT_k<i x_j>`. So
  `I_ij^(1) = c⁻²[∫x̂_ijT_kk − (4/3)Ȧ_ij] = −c⁻²[(1/3)∫x̂_ijT_kk + (8/3)S_ij]`.
  The radiation is a functional of the stresses, and the toroidal moment is one of their moments.
  Recomputed against the exact linearised power `P = (Gω²/4πc⁵)∮Λ_ij,kl T̃_ij T̃*_kl dΩ` for two
  opposed vortex cells: 0.9994 at kd ≈ 0.05 and 0.9977 at 0.1. Toroidal term alone: 16.0×.
  Stress alone: 9.0×. The amplitudes are 4 and −3 times the true one. P(2ω)/P(ω) = 63.96 at fixed
  stress, so 2⁸ at fixed mass current.
  Adding a divergence-free stress `S_ij = ε_ikz ε_jmz ∂_k∂_m g` changes neither Σ nor ρ, yet turns a
  silent cell into a radiator with `I_ij = ∫x̂_ijS_kk` (recomputed: 0.9996). **The same mass current
  radiates differently with different stresses.** There is no gravitational "toroidal family" that
  the mass current fixes on its own.
- **Binary share [new to our search; arithmetic on a known formula; convention-dependent].** We
  took the Blanchet–Damour point-particle quadrupole
  `Σm_A[1 + (3/2)v_A²/c² − Gm_B/(rc²)]x̂_A + (1/14c²)∂_t²Σm_Ar_A²x̂_A − (20/21c²)∂_tΣm_Av_A^k x̂_A^ijk`,
  plus the Newtonian octupole and current quadrupole, with the 1PN Kepler law and spectral time
  derivatives. That reproduces `−1247/336 − 35ν/12` to 10⁻⁵ at ν = 0.25, 0.1 and 0.01
  (recomputed). By hand, its CM form is Blanchet's `(29/42)(1−3ν)v²`, `−(1/7)(5−8ν)Gm/r`,
  `(11/21)(1−3ν)r²v^<ij>` and `−(4/7)(1−3ν)(r·v)x^<i v^j>`. Dropping only the −20/21 term changes
  the coefficient by exactly `−(32/21)(1 − 3ν)`: −0.381 of −4.440 for equal masses, −1.524 of −3.711
  as ν → 0. Splitting off Dubovik-style instead, absorbing the r² term, gives `−(44/21)(1 − 3ν)`
  (by hand, not in the script). Different splits give different numbers, so it is bookkeeping, not a
  measurable "toroidal luminosity".

### Gaps
- Damour & Iyer 1991 not opened (see above). That is the main novelty risk for Theorems 1–2.
- Theorem 2 covers ρ_osc = 0 at O(c⁻²) only. The general toroidal/stress split at ℓ ≥ 3 and at 2PN
  (non-compact terms) is not done.
- Any physical system where the stress-dependence is large is unidentified. In a self-gravitating
  binary the "stresses" are gravitational and come out as the known 1PN terms.

## 3. Gravitomagnetic AB phase outside a toroidal mass current

### Takeaway
Outside an ideal toroidal mass current the gravitomagnetic field is zero, but
`∮h_0i dx^i ≠ 0` on loops that link the tube. A slow particle of mass m gains
**`Δφ = 16πGmṀ(R − √(R² − a²))/(ħc²)`**, where Ṁ is the poloidal mass flux through the hole. In the
thin limit this is `16GmT_g/(ħc²R²)` in terms of the toroidal moment. It is the EM toroid formula
with `μ₀ → 4 × GM_PERMEABILITY` and `q → m`. The effect is **known qualitatively**
(Stachel 1982; Lynden-Bell & Katz 2014 for light; Chiao et al. 2013 and Mo & Modesto 2026 for
rotating cylinders). The closed-form torus phase and its size against atom interferometers were not
found. The size is hopeless: 2.0×10⁻¹² rad for a tonne at 1 km/s, 10⁻¹⁷ of Earth's Sagnac phase
over 1 m², and 4×10¹⁸ shots, or about 4×10¹¹ years, for a CARIOQA-class sensor to reach SNR 1.

### Cited Findings
- **[known] Rotating-cylinder gravitational AB phase.** `Δφ = 8πGmμωρ₀²/(ħc²)` for a thin
  shell of linear density μ, angular velocity ω and radius ρ₀, with the gravitomagnetic field zero
  outside: [Mo & Modesto, arXiv:2608.15587 (Aug 2026)](https://arxiv.org/abs/2608.15587),
  eqs. 4.17, 4.27–4.28; opened. Our torus formula reproduces it in the limit R → ∞ (recomputed to
  10⁻⁶).
- **[known] Thought experiment with a mass-current "solenoid".** `B_g = 4μ_g I_g0` inside, zero
  outside; the phase is `(m c/ħ)∮h_0i dx^i`, with `μ_g = 9.31×10⁻²⁷` SI, which is
  `gravitics::gem::GM_PERMEABILITY` —
  [Chiao et al. 2013, arXiv:1301.4270](https://arxiv.org/abs/1301.4270), eqs. 6, 10–11, 28;
  opened. No laboratory number given.
- **[known] Measured scalar gravitational AB.** Overstreet et al. 2022 (Science 375, 226), as in
  `gravimetry::atom::aharonov_bohm_rate`: 1.631 rad/s for ⁸⁷Rb at 7 cm from 1.25 kg (recomputed).
  Ho & Morgan, Aust. J. Phys. 47, 245 (1994), propose a gravito-electric version only (abstract).
- **[known] Sensitivity.** CARIOQA per-shot phase noise `1/(C√N)` = 3.95×10⁻³ rad (C = 0.8,
  N = 10⁵, 3.5 s cycle), from the `gravimetry::atom` test (recomputed).

### Inferences
- **[new to our search] The torus formula.** Lagrangian `L ⊃ mc h_0i v^i` and
  `h_0i = −(4G/c³)∫ρv/|x − x'|` (by hand, this reproduces the Lense–Thirring `g_0φ = −2GJ sin²θ/(c³r)`) give
  `|Δφ| = (4Gm/ħc²)∮a·dl`, with `∇×∇×a = 4πρv`. Ampère gives `b_φ = 2Ṁ/s` inside the tube, and the
  flux through the cross-section is `4πṀ(R − √(R² − a²))`. Checked by a Neumann double sum over 96
  discrete turns: 0.05% for two different linking loops, and 9×10⁻¹⁹ for a loop that does not link
  the tube. The ratio to the naive analogue is exactly 4, matching Lynden-Bell & Katz's 16π against
  4π. That makes it the house "spin-2 factor" in a potential: 4 here, where frame-dragging
  precession has 2 in `physica::analogies`.
- **Sizes (recomputed; ring masses and speeds are assumed design inputs, not sourced).** N
  spinning hoops round a torus with total mass M and rim speed v give, in the thin limit,
  `Δφ = 4GmMv(a/R)/(ħc²)`.
  - Overstreet's 1.25 kg at 100 m/s, a/R = ½: 2.5×10⁻¹⁶ rad, 2×10²⁶ shots.
  - 1 t at 1 km/s, a/R = ½: 2.03×10⁻¹² rad, 3.8×10¹⁸ shots (4×10¹¹ yr).
  - The electron current of `magnetics::toroidal`'s 24 A-turn toroid gives ⁸⁷Rb 3.2×10⁻²⁸ rad. For
    an electron the gravitomagnetic phase is `4α_G/α` = 9.6×10⁻⁴³ of the EM AB phase.
  - Earth's Sagnac phase for the same atom over 1 m² is 2.0×10⁵ rad, so reversing the rings' spin
    would have to resolve 10⁻¹⁷ of it.
- **[new to our search] No astrophysical case found.** A natural toroidal mass flow needs a closed
  probe path through it. The Sun's meridional circulation has two mirror cells, so its toroidal
  dipole vanishes by reflection symmetry (by hand; no data used).

### Gaps
- The atom loop must link the tube, which means a ring-shaped interferometer threaded through the
  hole. Its geometry and systematics (Newtonian phase of the rings' mass, Sagnac drift, vibration)
  are not modelled.
- Light: the Sagnac-type time difference `Δt = 2∮h_0i dx^i/c` would follow from the same Γ. We did
  not evaluate it.

## 4. Novelty check (frank verdicts)

| Item | Verdict | Evidence |
|---|---|---|
| GEM toroidal mass-current moment, zero exterior field | **Known** | Lynden-Bell & Katz 2012 (exact GR); GEM = magnetostatics with J → ρv |
| Gravitational "anapole" term | Not found by name | the queries below; nothing on arXiv or INSPIRE |
| Devaney–Wolf non-radiating gravitational source, compact h̄ | **New to our search**, routine | no gravity source found; the EM criterion is known |
| Lone closed toroidal current exactly non-radiating via longitudinal stress | **New to our search**, routine | TT-projection argument is textbook |
| Toroidal-type term in GW moments | **Known** | Blanchet LRR eq. 138 (Blanchet & Damour 1989) |
| No toroidal-dipole GW radiation; cancellation needs T^ii | **Content known** (no dipole GWs); explicit toroidal proof new to our search | Theorem 1; Damour & Iyer 1991 not opened |
| l = 2 toroidal correction fixed by stresses, not current | **New to our search** as a statement | Theorem 2; follows directly from eq. 138 |
| Toroidal share of 1PN binary flux −(32/21)(1−3ν) | **New to our search**, convention-dependent, trivial | checked against −1247/336 − 35ν/12 |
| Gravitomagnetic AB round a toroidal mass current | **Known** qualitatively | Stachel 1982; LB&K 2014; Chiao 2013; Mo & Modesto 2026 (cylinder) |
| Closed-form torus phase, ×4 rule, atom-interferometer size | **New to our search**, a transcription | none found with numbers |
| SMES toroid "axial acceleration" (Stephenson et al.) | **Contested** (abstract only) | GEM exterior field zero; force needs velocity |

**Strongest honest claim:** "In linearised gravity a toroidal mass-current moment is not an
independent radiator. Its dipole is cancelled by momentum conservation. Its quadrupole effect is
inseparable from the stresses that drive it, to the point that one closed current can radiate or
be silent depending on its stresses. Its only exterior signature is a gravitomagnetic AB phase,
16πGmṀ(R − √(R² − a²))/(ħc²), which is 10¹⁸ shots beyond present atom interferometers for any
laboratory rotor." Every part is checked numerically. None of it would survive review as a
discovery; it could be a short pedagogical paper once Damour & Iyer 1991 is read.

## 5. Implementation

1. **`gravimetry::atom::gravitomagnetic_ab_phase(m, mass_flux, major, minor) -> f64`** — returns
   `16πGmṀ(R − √(R² − a²))/(ħc²)`. Status: **mainstream** (Stachel; Lynden-Bell & Katz;
   Mo & Modesto). Put it next to `aharonov_bohm_rate`. Tests, all from the script:
   - The cylinder limit equals Mo & Modesto eq. 4.28.
   - The thin limit is `16GmT_g/(ħc²R²)`, with `T_g` in the form of
     `magnetics::toroidal::toroidal_solenoid_moment` (`ṀπRa²/2`).
   - 1 t at 1 km/s with a/R = ½ gives 2.0325×10⁻¹² rad, and `AtomInterferometer` CARIOQA needs
     `(3.95×10⁻³/Δφ)²` ≈ 3.8×10¹⁸ shots.
   - The Neumann-sum check can reuse the `solenoid` helper in `toroidal.rs`'s tests, scaled to mass
     current.
2. **`physica::analogies`, a new block "Toroidal solenoids in both fields".** Compute
   `magnetics::toroidal::aharonov_bohm_phase(q, MU_0·N·I·(R − √(R² − a²)))` and
   `gravitomagnetic_ab_phase(m, Ṁ, R, a)` with `q → m`, `NI → Ṁ` and `MU_0 → gravitics::gem::GM_PERMEABILITY`.
   Assert that the ratio is exactly 4, next to the existing spinning-shell factor of 2. The
   cross-crate test belongs in `physica`, which depends on `gravimetry`, `gravitics` and
   `magnetics`. `gravimetry` depends only on `geopotential` and `metrology`, so the new
   function's own tests use bare constants.
3. **Theorems 1–3: no crate function.** They are statements about moments, not calculators. At most,
   add one sentence to the `magnetics::toroidal` module doc: "gravity has no counterpart of the ω⁶
   toroidal-dipole radiation: in the ℓ = 1 mass moment the toroidal term is cancelled by momentum
   conservation". Link this note. The binary split is convention-dependent and should not become an
   API.

## Search log (2026-10-27)

WebSearch:
- "gravitational anapole toroidal moment linearized gravity"
- "\"toroidal\" mass current gravitomagnetic Aharonov-Bohm"
- "non-radiating sources linearized gravity gravitational waves anapole"
- "Lynden-Bell Katz gravitational toroidal solenoid …"
- "\"toroidal\" moment gravitational radiation multipole Thorne Blanchet"
- "Dubovik toroid moments gravitation general relativity"
- "gravitomagnetic Aharonov-Bohm effect atom interferometer proposal rotating mass phase estimate"
- "Ho Morgan … gravitational Aharonov-Bohm …"
- "Damour Iyer 1991 multipole analysis … toroidal moments"
- "\"nonradiating\" sources \"linearized gravity\" … Devaney Wolf"
- "\"gravitational anapole\""
- "\"toroidal dipole\" gravity … mass current"
- "toroidal multipole moments of stress-energy tensor gravitational radiation … anapole nonradiating" (extended mode)
- "Gbur \"nonradiating sources\" … gravitational"

arXiv full-text search (abstracts):
- "anapole gravitational": 2 hits, both EM or dark matter
- "toroidal gravitomagnetic": 3 hits, none relevant except 2306.04832 (rotating dust)
- the arXiv API was rate-limited

INSPIRE:
- anapole and gravitational wave: 2
- title anapole and title gravit*: 1 (axion)
- "toroidal moment" and gravitomagnetic: 0
- nonradiating and "linearized gravity": 0
- "toroidal dipole" and gravitational: 0
- toroidal and gravitomagnetic: 4
- "gravitational solenoid": 1 (1612.00905, not opened)
- toroidal and Aharonov-Bohm and gravitational: 2
- citations of 1205.1616: 4 (1402.6171, 1412.7742, 1312.6043, 1904.00270)

Opened: 1205.1616, 1312.6043, 1301.4270, 2608.15587, 1310.1528, 1401.1056; 0802.3346 (Biemond)
was opened and is not relevant (an azimuthally rotating torus, not a toroidal current).

Not opened, and worth reading before any write-up: Damour & Iyer PRD 43, 3259 (1991); Dubovik &
Tugushev Phys. Rep. 187, 145 (1990); Stachel PRD 26, 1281 (1982); Thorne RMP 52, 299 (1980);
Papapetrou Ann. IHP 11, 57 (1969); the Puget Sound 2018 project.
