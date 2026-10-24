# Superconductor–gravity claims: published mathematics, checkable numbers and GR baselines

Conventions in these notes. SI units unless stated. Constants used in my checks: G = 6.67430e-11 m³ kg⁻¹ s⁻², c = 2.99792458e8 m/s, m_e = 9.1093837e-31 kg, e = 1.602176634e-19 C, μ0 = 1.25663706e-6 H/m, ħ = 1.054571817e-34 J s, g0 = 9.80665 m/s².

"Computed" means I recomputed the number myself with `superconductor_gravity_checks.py` (saved next to this file; run with `python3 -I superconductor_gravity_checks.py`). Its printed values are the test vectors quoted below.

Status labels: **mainstream** (standard physics, agrees with experiment), **contested** (published, but not accepted and not independently confirmed), **refuted** (null replications or the authors' own retraction or reinterpretation).

Gravitomagnetic (GEM) field conventions differ by factors of 2 between authors. In these notes, k ≡ Ω_drag/ω is the dimensionless frame-dragging ratio: the local inertial-frame rotation rate divided by the source's angular velocity. Tajmar's 2006 "B_g" equals 2Ω (from his Larmor convention B_g = −2ω). His 2007 "B_g" equals Ω, the gyroscope precession rate. Both are given where relevant.

## Q1. Podkletnov & Nieminen (1992, 1997), Podkletnov & Modanese (2001 impulse generator, 2003 JLTP, 2012 "64c"), and Modanese's theory

### Takeaway
Podkletnov reported 0.05–0.3% (1992) and up to 2.1% (1997) weight loss above a levitated, rotating, RF-driven two-layer YBCO disc. He later reported a "gravity impulse" from 2 MV discharges that swung 18.5 g pendulums by up to 142 mm, with a claimed propagation speed of 64c. General relativity predicts a weight *increase* of at most about 3×10⁻⁹ from the disc's own mass, so the claim exceeds GR by about 10⁷ and has the opposite sign. The rotation-dependent GR part is smaller again, by about 10²⁰. Modanese's "theory" is qualitative: a condensate locally modifies an effective cosmological term. It needs a local Λ about 20 orders of magnitude larger than the observed value and produces no number that can be checked against the measured weight loss. Status: **refuted/contested**. No independent replication; nulls are listed in Q2.

### Cited Findings
**1992 Physica C 203, 441 (Podkletnov & Nieminen)**
- The disc was 145 mm in diameter and 6 mm thick, cooled with liquid helium and levitated over a solenoid (Meissner effect). A SiO₂ sample of about 5 g, 15 mm above the disc, lost about 0.05% of its weight. When the disc was rotated by lateral AC fields the loss rose to 0.3%. No shielding was seen when the disc rested on a fixed support. — [Modanese, EPL 35 (1996) 413, arXiv hep-th/9505094](https://arxiv.org/abs/hep-th/9505094)
- The sample mass was 5.48 g on an electro-optical balance, with a loss "from 0.05 to 0.3% … depending on the rotation speed". This is the abstract as relayed by a search-engine summary; the page itself returned 403 to me and I did not see the primary text. — [CORE record](https://core.ac.uk/works/1028378/)
- A 1996 paper submitted to J. Phys. D (with Vuorinen) claimed 2% and was withdrawn after acceptance. The cond-mat/9701074 preprint is believed to be essentially that paper. — [Wikipedia: Eugene Podkletnov](https://en.wikipedia.org/wiki/Eugene_Podkletnov)

**1997 report (arXiv cond-mat/9701074v3, Podkletnov alone)** — [arXiv](https://arxiv.org/abs/cond-mat/9701074)
- **Disc.** Toroidal disc: outer diameter 275 mm, inner diameter 80 mm, thickness 10 mm. Made of dense YBa₂Cu₃O₇₋ₓ, with the top layer melt-textured by induction heating to 1200 °C.
  - Upper layer: 6–7 mm, orthorhombic, T_c = 94 K, 75% c-axis oriented, J_c > 15,000 A/cm² at 75 K.
  - Lower layer: T_c = 60 K, porous, about 40% tetragonal. It is not superconducting in operation.
  - Lattice constants (nm): upper a = 0.381, b = 0.387, c = 1.165.
  - Estimated operating current density: 5,000–7,000 A/cm².
  - A 1 mm magnetic foil was glued on top to help rotation.
- **Fields, temperature and drive.** Fields up to 2 T; frequencies 10³–10⁸ Hz; disc temperature 40–70 K, "lower than 70 K in stationary conditions". Rotation was driven by 10⁵ Hz current (8–10 A, 1.2 mm wire) in two solenoids, up to 5000 rpm. Levitation height was ≥15 mm, with a maximum of 30–35 mm, which is reachable only below 70 K.
- **Weighing.** Dupont DTA/TGA balance with arms up to 220 cm and a sensitivity of about 10⁻⁶ g for 10–50 g samples. Samples sat 25–3000 mm above the cryostat inside a closed glass tube. Checks were done with an analytical balance (0.01 g on 100–250 g samples), a torsion balance and a spring balance.
- **Results.**
  - Levitating, not rotating: 0.05–0.07%.
  - At 5000 rpm: typically 0.3–0.5%, best 0.5%.
  - Above the inner edge: 0.1–0.25%.
  - While braking from 5000 to 3500 rpm: 1.9–2.1% for 25–30 s.
  - Maximum effect at 3.2–3.8 MHz.
  - The same percentage loss for all materials (metal, glass, plastic, wood) and masses.
  - The same loss at 3 m height inside the disc's vertical projection; none below the cryostat.
  - Air-pressure drop of 5 mmHg (8 mmHg while braking).
  - No effect with static fields only, with a metal disc, or with a single-phase superconductor disc.
- **Tables.** At 2 MHz, rpm 4000/4200/4400/4600/4800/5000 → 0.17/0.19/0.20/0.21/0.22/0.23%. At 4300 rpm, MHz 3.1/3.2/3.3/3.4/3.5/3.6 → 0.22/0.23/0.24/0.26/0.29/0.32%.
- **Systematics the paper itself discusses.**
  - Diamagnetism: a gradient of 0.17 T/cm gives 16 dyn on 1 g of NaCl, i.e. 1.6%; 0.1% needs about 0.04 T/cm.
  - Archimedes force: about 0.1% of weight at 1 kg/dm³.
  - RF levitation: thick Cu, Al and steel screens (300 mm diameter, 50 mm thick) did not change the effect.
  - No formulas for the effect are given.

**2001 "Impulse Gravity Generator" (Podkletnov & Modanese, arXiv physics/0108005v2)** — [arXiv](https://arxiv.org/abs/physics/0108005)
- **Emitter.** Two-layer YBCO disc, 80–120 mm in diameter and 7–15 mm thick, with Ce and Ag additions; the SC layer was 4 mm (emitter 1) or 8 mm (emitter 2). T_c = 87–90 K, trapped field 0.5 T at 77 K, J_c > 5×10⁴ A/cm². Operated at 50–70 K.
- **Chamber.** Quartz cylinder about 1 m in diameter and 1.5 m long, with an electrode gap of 0.15–0.40 m, a 1.05 m solenoid giving 0.9 T, and a small coil to freeze field into the emitter.
- **Marx generator.** Twenty 25 nF capacitors charged to 50–100 kV, giving up to 2 MV, with a peak current of about 10⁴ A and a discharge lasting 10⁻⁵–10⁻⁴ s. Above 500 kV the discharge became a "flat glowing discharge".
- **Detectors.** Pendulums made of spheres 10–25 mm in diameter on 800 mm (or 500 mm) cotton strings in evacuated glass cylinders, at 6 m and at 150 m (in another building). Shielding: a 0.3 m brick wall plus a 1 × 1.2 × 0.025 m steel sheet, and a further 0.8 m of brick for the 150 m station. Test masses were 10–50 g; the tables use an 18.5 g sphere. Results at 3–6 m and at 150 m were "identical".
- **Table 1 (emitter 2; 12-shot average, σ = 5–7%).** Voltage (kV) → Δl (mm) / Δh (mm) / ΔE (10⁻⁴ J):

  | Voltage (kV) | Δl (mm) | Δh (mm) | ΔE (10⁻⁴ J) |
  |---|---|---|---|
  | 500 | 56.5 | 2.0 | 3.6 |
  | 750 | 91.3 | 5.2 | 9.5 |
  | 1000 | 110.4 | 7.7 | 13.9 |
  | 1250 | 123.0 | 9.5 | 17.3 |
  | 1500 | 131.6 | 10.9 | 19.8 |
  | 1750 | 137.6 | 11.9 | 21.7 |
  | 2000 | 142.0 | 12.7 | 23.1 |

- **Table 2 (emitter 1).** Δl = 40.0, 70.9, 85.3, 94.6, 100.8, 104.7, 107.1 mm at the same voltages, with ΔE = 1.8 … 13.1 ×10⁻⁴ J. The beam varied by up to 12–15% across its cross-section. A microphone with about 16 kHz bandwidth also responded.
- **The paper's only quantitative formula, eq. (1).** Ft = √(2gl[1 − √(1 − (d/l)²)]), where F is the impulse "strength" (m/s²), t is the duration, l is the pendulum length and d is the half-amplitude. With t = 10⁻⁴ s the paper gets "F ∼ 10³ g".
- **Charge and energy argument.** "Q up to 0.1 C", U = QV/2 ∼ 10⁵ J, implying an effective capacitance of about 10⁻⁷ F against at most 10⁻¹⁰ F for a normal conductor.
- **Modanese's model in the same paper.**
  - Fluctuation action: ΔS = Λ τ M r² Q (eq. 2, with ħ = c = 1). τ is the fluctuation duration, M the virtual ± masses, r their separation and Q a shape factor.
  - The condensate subtracts 8πG·L(x) from Λ (eq. 8): S_int + S_Λ = ½∫d⁴x h_μν ∂^μφ*∂^νφ + ½∫d⁴x Tr h (Λ/8πG − L).
  - GL Lagrangian: L = −(1/2m)|−i∇ψ + 2eAψ|² − αψ*ψ − ½β(ψ*ψ)² (eq. 11). For solutions this becomes L = −(1/2m)[ħ²(∇ρ)² + ħ²ρ∇²ρ − mβρ²] (eq. 12).
  - Estimates: L₁ = ½βn_p² ∼ 10⁴ J/m³ and L₂ ∼ ρ|ρ″|/2m ∼ 10⁷ J/m³, using ρ ∼ 10²⁷ m⁻³, ξ ∼ 10⁻⁹ m and κ ∼ 10².
  - Fluctuation amplitude: A ∝ 1/|Λ/8πG − L| (eq. 13), so the model requires Λ/8πG > L₂ at atomic scales.

**2003 J. Low Temp. Phys. 132, 239–259 (Podkletnov & Modanese; arXiv physics/0209051)** — [arXiv](https://arxiv.org/abs/physics/0209051)
- The cathode is 10 cm in diameter, the capacitor array is charged up to 2000 kV, and peak currents are of order 10⁴ A at 50–70 K.
- Above 500 kV there is a flat discharge and a collimated radiation pulse carrying at least 10⁻³ J, described as not electromagnetic.
- Wikipedia notes that this version appears to drop the explicit "gravity" label. I did not verify this against the full text.

**2012 Bentham book chapter, "Study of light interaction with gravity impulses and measurements of the speed of gravity impulses"**
- Laser test: a laser crossed the beam at 0.1°, giving an interaction region of about 57 m, and the laser spot intensity fell by 7–10% during discharges.
- Speed test: two piezoelectric sensors 1211 m apart, read with synchronised rubidium atomic clocks.
- Sources: [Gallerati, Modanese & Ummarino review, arXiv 2206.07574](https://arxiv.org/abs/2206.07574). The chapter listing ([Bentham Direct](https://www.benthamdirect.com/content/books/9781608053995.chapter-8), [EurekaSelect](https://www.eurekaselect.com/public/chapter/3862)) reports 64c over 1211 m, with a delay of 63 ± 1 ns (one listing prints "631 ns", which is a typo). I saw these listings only through a search summary because the publisher's page returned 403.

**Modanese, EPL 35 (1996) 413, "Theoretical analysis of a reported weak gravitational shielding effect"** — [arXiv hep-th/9505094](https://arxiv.org/abs/hep-th/9505094)
- **Classical bound.** Media shield gravity by less than "one part in 10¹⁰ or less" (Majorana-type tests). Einstein's equations R_μν − ½g_μνR = −8πG T_μν (eq. 1) cannot produce about 0.01 g from an apparatus of that mass-energy.
- **Euclidean QG interaction energy (eqs. 2–3).** E = lim_{T→∞} −(ħ/T) log ⟨exp(−ħ⁻¹ Σᵢ mᵢ ∫ dsᵢ)⟩_{S_g}.
- **Gravitational action (eq. 4).** S_g = ∫d⁴x √g [Λ/8πG − R/8πG + ¼ a R_μνρσR^μνρσ].
- **Condensate vertex (eq. 7).** L_hφ0 = κ ∂_μφ₀* ∂_νφ₀ h^μν, with κ ∼ 10⁻³³ cm. It is "too small (by several magnitude orders)".
- **Key coupling (eqs. 8–9).** S = ∫d⁴x √g [Λ/8πG + ½μ²(x) − R/8πG] + …, with μ²(x) = [∂_μφ₀]*[∂^μφ₀] + m²|φ₀|². The condensate therefore adds positively to the cosmological term.
- **Numbers.**
  - Observed |Λ|G < 10⁻¹²⁰, so |Λ|/8πG < 10¹² cm⁻⁴.
  - Pair density 10²⁰ cm⁻³ and m_pair ∼ 10¹⁰ cm⁻¹ give φ₀ ∼ 10⁵ cm⁻¹, m²|φ₀|² ∼ 10³⁰ cm⁻⁴ and ∂φ₀ ∼ 10¹² cm⁻².
  - Scale-dependent Λ: |Λ|G ∼ (ℓ₀p)^γ (eq. 11). The required value is |Λ| ∼ 10⁻³⁶ cm⁻² at scales of 10⁻⁸–10⁻⁶ cm, corresponding to a graviton mass of about 10⁻¹⁸ cm⁻¹ and γ > 2.
- Modanese states that Torr & Li's claimed anomalously strong gravitoelectric field "contrasts with the full equations". He gives the gravitational London depth in a neutron star (10¹⁷ kg/m³) as about 12 km, citing Lano, hep-th/9603077.
- **Update (supr-con/9601001, Jan 1996).** Modanese's geometric "shadow" model gives shielding/max-shielding against height h (in disc radii) of 0.62, 0.34, 0.18, 0.102, 0.050, 0.029, 0.018 for h = 1, 2, 3, 4, 6, 8, 10, with tan θ = 1/h for the shadow cone. He calls the absence of height dependence a "severe challenge". — [arXiv supr-con/9601001](https://arxiv.org/abs/supr-con/9601001)

### Inferences
- **Computed: Table 1 reproduces exactly from pendulum geometry.** Δh = l − √(l² − Δl²) with l = 0.8 m, and ΔE = m g Δh with m = 18.5 g. At 2 MV: Δh = 12.70 mm, ΔE = 23.05×10⁻⁴ J (published 23.1), v = Ft = 0.499 m/s. With t = 10⁻⁴ s this gives F = 509 g, which the paper rounds to "∼10³ g". Every row of both tables matches to rounding. So Tables 1 and 2 are only a geometric conversion of the observed deflection Δl into Δh and ΔE; no physics model enters.
- The paper's small-angle form "Ft ∼ (g/l)d" is dimensionally wrong. It should be √(g/l)·d = 2πd/T.
- **Computed: momentum and energy do not match any wave.** The pendulum's momentum is p = 9.2×10⁻³ kg m/s, while E/c = 7.7×10⁻¹² kg m/s, so p/(E/c) = 1.2×10⁹. Any radiation (electromagnetic or gravitational) obeys p = E/c, as the authors themselves note.
- **Computed: the charge claim is inconsistent with the Marx bank.** The bank (20 × 25 nF erected in series, so 1.25 nF at 2 MV) stores ½CV² = 2.5 kJ and can deliver Q = CV = 2.5 mC. The paper's Q = 0.1 C and U ∼ 10⁵ J are 40× more than the bank holds. The paper attributes this to an anomalous effective capacitance of 10⁻⁷ F.
- **Computed GR baseline, impulse generator.** Take the whole 2.5 kJ bank energy as mass-energy, Δm = 2.8×10⁻¹⁴ kg, at 6 m. Then a = GΔm/r² = 5.2×10⁻²⁶ m/s², and over 10⁻⁴ s Δv = 5×10⁻³⁰ m/s. The claim (0.5 m/s) exceeds GR by about 10²⁹. Even the Newtonian pull of a 1-tonne apparatus at 6 m (1.9×10⁻⁹ m/s²) is 2.7×10¹² times smaller than the claimed F.
- **Computed: 64c.** 1211 m / 64c = 63.1 ns, compared with 4.04 µs for light; this matches the "63 ± 1 ns" listing. Synchronising two rubidium clocks 1.2 km apart to well under 63 ns, and identifying trigger latency, is the obvious systematic. GR predicts gravitational disturbances travel at c; GW170817/GRB 170817A bound −3×10⁻¹⁵ < Δv/c < +7×10⁻¹⁶ ([LIGO/Virgo/Fermi/INTEGRAL, arXiv 1710.05834](https://arxiv.org/abs/1710.05834); I quote this from memory and did not re-fetch it in this session). A 64c signal means Δv/c = 63, about 10¹⁶ times beyond that bound, although the claimed "impulse" is not a gravitational wave in GR's sense.
- **Computed: the laser crossing geometry fits.** 57 m = D/tan(0.1°) with D = 0.10 m, so the 57 m interaction length implies a beam about 10 cm wide.
- **Computed GR baseline for the 1997 disc.** Assume a YBCO density of 6.38 g/cm³; then the disc mass is 3.47 kg and the surface density σ = 63.8 kg/m².
  - The disc's own Newtonian pull adds to the weight, by at most 2πGσ/g = 2.7×10⁻⁹. On axis above the hole it is 9.6×10⁻¹⁰ at 25 mm and 2.6×10⁻¹² at 3 m.
  - Claimed: −2.1×10⁻². Ratio about 7.7×10⁶, with the opposite sign.
  - Rotation at 5000 rpm: ω = 523.6 rad/s, rim speed 72 m/s, rim acceleration about 3,800 g, kinetic energy 4.9 kJ ≡ 5.4×10⁻¹⁴ kg. This changes the GR pull by a fraction of about 10⁻¹⁴, so Δg/g ∼ 4×10⁻²³, and the claimed rotation-dependent effect exceeds GR by about 10²⁰–10²¹.
  - A stationary test mass feels no gravitomagnetic force (F = m v × B_g with v = 0).
- **Computed: height dependence contradicts any shadow-type shield.** For a uniform Earth, the fraction of g that arrives through a downward cone of half-angle θ is exactly 1 − cos³θ (∫cos φ · 2R_E cos φ dΩ normalised). With tan θ = 1/h this gives 0.646, 0.285, 0.146, 0.087, 0.040, 0.023, 0.0148 for h = 1–10. These are within about 20% of Modanese's table; his values are higher for h ≥ 2, as expected from a centrally condensed Earth.
- At 3 m above the 275 mm disc (h ≈ 22 radii) the fraction is 0.31%. A 2.1% shield near the disc should therefore read about 0.007% at 3 m, yet Podkletnov reported the same 2%. A 2% effect at 3 m would require the "shield" to act on directions that do not pass through the disc.
- **Computed: Modanese's numbers reproduce only in order of magnitude.**
  - m_pair = 5.2×10¹⁰ cm⁻¹, φ₀ = √(n/m) = 4.4×10⁴ cm⁻¹, m n = 5.2×10³⁰ cm⁻⁴.
  - The required local Λ ∼ 10⁻³⁶ cm⁻² is 10^20.0 times the observed Λ ≈ 1.1×10⁻⁵⁶ cm⁻².
  - L₂ ∼ ħ²n/(2mξ²) = 3×10⁶ J/m³ reproduces "∼10⁷".
  - L₁ does **not** reproduce from the paper's stated κ–β relation: I get 4×10⁶ J/m³. It reproduces as the GL condensation energy B_c²/2μ0 = 1.6×10⁴ J/m³ for B_c = 0.2 T.
  - The requirement Λ/8πG > L₂ needs a vacuum energy density about 10^15.8 times the observed 5.3×10⁻¹⁰ J/m³.
  - The model has no free-parameter-free prediction of the weight-loss percentage, so it cannot be checked against data.
- **Calculator spec (Podkletnov module).**
  - `pendulum_impulse(l, dl, m, t) -> (dh, dE, v, F)`; test vector: l = 0.8, dl = 0.142, m = 0.0185, t = 1e-4 → (12.70 mm, 2.305e-3 J, 0.499 m/s, 4.99e3 m/s²).
  - `annulus_pull(sigma, Ri, Ro, z)` for the GR baseline.
  - `shadow_fraction(h) = 1 − cos³(atan(1/h))`.
  - `marx(C_stage, n, V_stage) -> (E, C_series, Q)`; test vector: (25e-9, 20, 1e5) → (2500 J, 1.25e-9 F, 2.5e-3 C).
  - `delay(L, v_over_c)`; test vector: (1211, 64) → 63.1 ns.

### Gaps
- I could not see the 1992 Physica C full text. The rotation speed, temperature and exact balance specification for 1992 are unverified; the disc size comes from Modanese's EPL summary.
- The 2012 chapter (64c) full text was inaccessible (403). The clock-synchronisation method, error budget and sensor latency are unknown to me.
- No photographs, videos or independent data exist for either Podkletnov experiment, per [arXiv 2206.07574](https://arxiv.org/abs/2206.07574) and [Tajmar et al. 2022](https://www.frontiersin.org/articles/10.3389/fphy.2022.892215/full).
- Unnikrishnan, "Does a superconductor shield gravity?" (Physica C 266, 133, 1996) and de Podesta & Bull (Physica C 253, 199, 1995): I could not retrieve their abstracts or arguments.

## Q2. Replication attempts and the null bounds they set

### Takeaway
No replication has ever reproduced the effect. The null bounds on fractional weight change are:
- NASA MSFC static test: < 2×10⁻⁸ (about 20 nano-g).
- Woods et al. (BAE/Sheffield): ±4×10⁻⁵.
- Hathaway et al.: < 1×10⁻⁵ (0.001%).
- Tajmar's accelerometers: < 5×10⁻⁶ (0.0005%).

These are below Podkletnov's 2×10⁻² by factors of 500 to 10⁶. However, none of the replications reproduced the full set of his conditions (two-layer 275 mm disc, magnetic levitation and rotation at 5000 rpm, MHz RF), and Podkletnov's defenders use this to keep the claim nominally "untested". Tajmar's 2022 thrust balance and the Tajmar–Lőrincz tests found no force-beam effects; the Poher-type accelerometer signals were electromagnetic or acoustic artefacts. Boeing's GRASP was a press report about an internal evaluation, not an experiment. Status: **refuted for every configuration actually tested**.

### Cited Findings
- **NASA MSFC/UAH: Li, Noever, Robertson, Koczor & Brantley, Physica C 281, 260–267 (1997), "Static test for a gravitational force coupled to type II YBCO superconductors".** Bulk YBCO was levitated in a DC field and a gravimeter measured g. Changes were "less than 2 parts in 10⁸" of normal gravity. — [NTRS record](https://ntrs.nasa.gov/citations/19990039542)
- **Noever, Koczor & Roberson (NASA MSFC, 1998).** A 15 cm YBCO disc sat in a three-phase AC motor stator with the field rotating at up to 12,000 rpm (< 60 G at the disc edge, < 10 G at the centre), with and without DC levitation. They report an "apparent gravity increase of less than 1×10⁻⁶" and no detectable effect of field, rotation or temperature. — [NTRS 19990019627](https://ntrs.nasa.gov/citations/19990019627)
- **Noever & Koczor (NASA MSFC, 1998).** RF illumination at 1–15 MHz and ≥ 90 W from spiral antennas of 15 cm radius produced an apparent *increase* of 3–5×10⁻⁵ cm/s² (≈ 3–5×10⁻⁸ g), "significantly lower" than the 2.1% claim. They correlate it with the 2–7 MHz AC-resistance peak, i.e. an RF heating/EM artefact. — [NTRS 19990064220](https://ntrs.nasa.gov/citations/19990064220)
- **NASA/UAH disc programme.** Single-layer discs were made in 1999; Podkletnov was hired as a consultant and a bi-layer disc was made. Only static tests on a purchased 4-inch single-layer disc were published. — [arXiv 2206.07574](https://arxiv.org/abs/2206.07574)
- **NASA BPP programme.** It ran 1996–2002 with a total of $1.2 M, and found "no evidence of the effect" (citing Millis 2005, Ann. NY Acad. Sci. 1065, 441; Robertson et al. 2001 AIAA; Hathaway et al. 2003). — [Wikipedia: Breakthrough Propulsion Physics Project](https://en.wikipedia.org/wiki/Breakthrough_Propulsion_Physics_Project)
- **Woods, Cooke, Helme & Caldwell, AIAA 2001-3363 (37th JPC, Salt Lake City; BAE Systems/Sheffield).**
  - Discs up to 3 inches; many of the stated conditions were not met. — [arXiv 2206.07574](https://arxiv.org/abs/2206.07574)
  - Equipment resolution was "of the order of ±0.004%", and the "subset of … conditions examined did not produce gravity modification". — [Woods, arXiv physics/0204031](https://arxiv.org/abs/physics/0204031)
- **Hathaway, Cleveland & Bao, Physica C 385 (2003) 488–500, "Gravity modification experiment using a rotating superconducting disk and radio frequency fields".** The disc was supported and rotated only mechanically, at lower speed and lower temperature than Podkletnov's. "No weight modifications or gravity-like forces were detected down to the 0.001% level." — [arXiv 2206.07574](https://arxiv.org/abs/2206.07574)
  - A secondary blog says the disc ran at about 400–800 rpm and that rapid induction heating cracked a disc within 3 s. This is unverified. — [Lifeboat blog](https://lifeboat.com/blog/2013/01/gravity-modification-what-went-wrong)
  - Hathaway later wrote the history chapter "Gravitational experiments with superconductors: history and lessons" in *Frontiers of Propulsion Science* (AIAA 2009, pp. 229–246). — [arXiv 2206.07574 ref. 22](https://arxiv.org/abs/2206.07574)
- **Tajmar et al. (2006).** Vertical accelerometers (±5 µg) saw no change above Nb, Pb, YBCO or BSCCO rings during the phase transition or rotation, so any shielding "must be < 0.0005% compared to claims of up to 2%". — [arXiv gr-qc/0603033](https://arxiv.org/abs/gr-qc/0603033)
- **Tajmar, Plesescu & Seifert, Meas. Sci. Technol. 21, 015111 (2010).** Weighing during the phase transition of bulk Nb and high-T_c superconductors on a magnetic suspension balance showed no anomaly. — [Tajmar et al. 2022, Front. Phys. 10:892215](https://www.frontiersin.org/articles/10.3389/fphy.2022.892215/full)
- **Impulse-generator replications.**
  - Junker (Inst. für Gravitationsforschung, from 2004, 50% scale): inconclusive.
  - Tajmar & Lőrincz scaled-down tests (2.2 kV, not in vacuum): with good EM shielding, the accelerometer effects were attributed to acoustic noise (2016).
  - Source: [arXiv 2206.07574](https://arxiv.org/abs/2206.07574).
- **Poher test.** Lőrincz & Tajmar (52nd AIAA JPC) found that with double Faraday cages and isolated data acquisition, "the Poher effect was only an electromagnetic artefact". — [Tajmar et al. 2022](https://www.frontiersin.org/articles/10.3389/fphy.2022.892215/full)
- **Tajmar, Neunzig & Kößling, Front. Phys. 10:892215 (2022).** A cryogenic thrust balance with coaxial current leads tested BSCCO and YBCO up to 15 A. "No force was seen within our resolution of around 100 nN", i.e. about 5×10⁻⁹ N/A against Poher's claimed about 1 N/A. The paper's reviewers were Unnikrishnan and Hathaway. — [Front. Phys.](https://www.frontiersin.org/articles/10.3389/fphy.2022.892215/full)
- **Boeing GRASP ("Gravity Research for Advanced Space Propulsion").**
  - BBC, 7 July 2002, and Nick Cook, Jane's Defence Weekly, 29 July 2002: the project was based on an internal briefing document. Boeing, BAE and Lockheed reportedly tried to approach Podkletnov.
  - Boeing told Space.com (31 July 2002) that GRASP was not funded with company money.
  - Source: [Wikipedia: Eugene Podkletnov](https://en.wikipedia.org/wiki/Eugene_Podkletnov).
  - The GRASP paper focused on the impulse gravity generator, according to secondary reposts. — [phys-l archive repost](https://phys-l.org/archives/2002/08_2002/msg00059.html)
- **ESA gravity-control review.** Tajmar & Bertolami (physics/0412176, 2005) concluded "there are currently no good grounds" for taking such effects seriously. — [Nature News, 26 Jan 2005](https://www.nature.com/news/2005/050124/full/news050124-8.html)

### Inferences
- Bound ladder for the calculator, as fractional |Δg/g| upper limits:

  | Test | Bound | Factor below the 2.1% claim |
  |---|---|---|
  | Li et al. 1997 | 2e-8 | ≈ 1×10⁶ |
  | Noever 1998 (AC motor) | ≲ 1e-6 | |
  | Tajmar 2006 | 5e-6 | ≈ 4×10³ |
  | Hathaway 2003 | 1e-5 | ≈ 2×10³ |
  | Woods 2001 | 4e-5 | ≈ 500 |
  | Podkletnov 1997 claim | 2.1e-2 (sign −) | — |

- The GR prediction (Q1) is +3×10⁻⁹, which is below every one of these instruments except roughly Li 1997. GR is therefore consistent with every null.
- No replication combined all of Podkletnov's prescribed conditions (two-layer 275 mm toroid, magnetic levitation and magnetic rotation at about 5000 rpm, 10⁵ Hz drive plus 3–4 MHz field, < 70 K). The claim is therefore "not reproduced" rather than "falsified under identical conditions". The burden stays on the claimant: no data, photos or independent apparatus were ever provided.

### Gaps
- The Hathaway 2003 primary text (rotation speed, RF power, balance type) was not accessible; the 0.001% figure is from the 2022 review.
- I did not retrieve the primary text of the Woods AIAA paper (disc size, temperature, rpm). There is a title discrepancy ("superconductors" vs "semiconductors") in a secondary blog.
- The Robertson et al. 2001 AIAA paper and Koczor & Noever's 1999 AIAA disc-fabrication paper were not retrieved. Their sensitivities are unknown.
- I found no primary Jane's text.

## Q3. Tajmar & de Matos: gravitomagnetic London moment, Cooper-pair mass anomaly, measurements, retraction, and the GR comparison

### Takeaway
de Matos & Tajmar proposed B_g = 2ω ρ*/ρ, about 10⁻⁶–10⁻⁴ ω, to explain Tate's 92 ppm Cooper-pair mass excess. In 2006 they reported about 100 µg tangential accelerations next to an angularly accelerated (1500 rad/s²) Nb ring at 4 K, giving a coupling of −6.6 to −9.6 ×10⁻⁸ s², i.e. k ≈ 2.6×10⁻⁵.

GR gives k = 2GM/(c²R̄) = 7.2×10⁻²⁷ for their 0.35 kg Nb ring, so the claim exceeds GR by 3.6×10²¹ (10^21.6). This reproduces "about 22 orders". Tajmar's own figures were "30 orders" in the paper and 10²⁰ in ESA's press release.

By 2007 the signals did not track T_c and showed "parity violation". The data excluded their own Cooper-pair hypothesis by 5 orders of magnitude. By 2011 an improved setup set a limit of 4×10⁻¹¹, about 800× below the 2007 gyro signal and 10⁵ below the prediction, and Tajmar attributed the earlier signals to acoustic and vibration artefacts. Graham et al.'s ring laser saw nothing significant. Status: **refuted, including by the original author**.

### Cited Findings
- **London-moment basis.** B = −(2m*/e*)ω (eq. 1). Tate et al. measured m*/2m_e = 1.000084(21) for Nb, against a theoretical 0.999992. — [Tajmar et al., arXiv gr-qc/0603033](https://arxiv.org/abs/gr-qc/0603033)
- **Primary Tate result.** h/m′ = 2SΔν; Δν and S measured to 5 ppm precision and 21 ppm accuracy; (m′/2m_e)_Nb = 1.000084(21) against a theoretical 0.999992. — [Tate, Cabrera, Felch & Anderson, PRL 62, 845 (1989)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.62.845)
- **Modified London moment.** B = −(2m*/e*)ω − (m*/e*)B_g (eq. 2). To explain Tate: B_g = 2ω(Δm/m*) = 1.84×10⁻⁴ ω (eq. 3). Earth's B_g is about 10⁻¹⁴ rad/s. — [arXiv gr-qc/0603033](https://arxiv.org/abs/gr-qc/0603033)
- **"Gravitomagnetic London moment" (Proca/graviphoton argument).**
  - Proca-gravity curl gives B_g = B₀,g e^(−x/λ_g) − 2ωρ_m μ₀g λ_g² (eq. 7).
  - Equivalent graviphoton mass: 1/λ_g² = (m_g c/ħ)² = μ₀g ρ_m (eq. 8). For normal matter this yields g = −a and B_g = −2ω (eq. 9, equivalence principle plus gravitomagnetic Larmor theorem).
  - For the coherent part: B_g = 2ω ρ*_m/ρ_m (eq. 11) and g = −(ρ*_m/ρ_m) a (eq. 13).
  - Fluxoid: ∮(m*v + e*A + m*A_g)·dl = nh/2 (eq. 14), giving B = −(2m/e)ω − (m/e)B_g (eq. 15).
  - Identification Δm*/m* = ρ*_m/ρ_m (eq. 17). Tate's 9.2×10⁻⁵ against Nb ρ*/ρ = 3.95×10⁻⁶ is "a factor of 23 away".
  - Source: [Tajmar & de Matos, arXiv gr-qc/0603032](https://arxiv.org/abs/gr-qc/0603032)
- **Induction law and predicted signal (gr-qc/0603033).**
  - rot g = −∂B_g/∂t (eq. 5).
  - g = −(r/2)Ḃ_g φ̂ = −r(ρ*_m/ρ_m)ω̇ φ̂ = −[2.624×10⁻¹⁸/(λ_L²ρ_m)] r ω̇ φ̂, with g in units of g0 (eq. 6).
  - Predicted Nb coupling (g/ω̇) is −0.7×10⁻⁸ s² from the density ratio and −37×10⁻⁸ s² from Tate.
  - Source: [arXiv gr-qc/0603033](https://arxiv.org/abs/gr-qc/0603033)
- **Setup (gr-qc/0603033).**
  - Rings: Nb 150×138×15 mm; Pb of the same size; YBCO 160×130×15 mm; BSCCO-2212 150×138×15 mm.
  - Drive: up to 6500 rpm with the electric motor; up to about 4500 rpm and faster acceleration with the air motor.
  - Sensors: AppliedMEMS 1500S accelerometers in an evacuated Faraday cage. Noise σ was 5 µg vertical and 15 µg radial/tangential, giving a coupling noise of ±2×10⁻⁸ s².
- **Results (gr-qc/0603033).**
  - Nb at about 1500 rad/s² produced about 100 µg, i.e. −6.6×10⁻⁸ (signal-to-noise 3.3).
  - Mean coupling −9.64×10⁻⁸ s²: 13× the eq. (4) prediction and 3.8× below Tate's.
  - Pb gave about 84% of Nb, as predicted.
  - At the T_c transition while rotating: a 277 µg peak at about 550 rad/s, giving ≥ 5×10⁻⁸ s².
  - About 250 runs. The paper says the signals are "30 orders of magnitude higher than what general relativity predicts".
- **ESA press release (March 2006).** The field is "100 millionths" of g and "one hundred million trillion times" (10²⁰) larger than GR. — [ESA: "Towards a new test of general relativity?"](https://www.esa.int/esaCP/SEM0L6OVGJE_index_0.html)
- **STAIF 2007 / gr-qc/0610015.**
  - Prediction: B_g = 2ω ρ*/ρ = 3.9×10⁻⁶ ω for Nb at T = 0.
  - Ring: OD 150 mm, 6 mm wall, 15 mm high; up to 650 rad/s and about 1500 rad/s²; in-ring accelerometer at r = 3.6 cm; requirements of 2×10⁻⁴ rad/s for gyros and 1 µg for accelerometers.
  - Measured couplings:
    - −3.9 ± 0.52 ×10⁻⁸ g rad⁻¹ s² (35–75 µg at 1200 rad/s²);
    - superconducting −2.26 ± 0.3 ×10⁻⁸ against normal −1.24 ± 1 ×10⁻⁹;
    - T → 0 extrapolation −14.4 ± 2.8 ×10⁻⁹.
  - Above-ring laser gyro: about 1/4 of the eq. (1) prediction, signal-to-noise better than 16:1; the reference gyro noise was 3×10⁻⁵ rad/s.
  - Source: [arXiv gr-qc/0610015](https://arxiv.org/abs/gr-qc/0610015)
- **2007 reanalysis (arXiv 0707.3806, "Search for frame-dragging-like signals close to spinning superconductors").**
  - The signals appear below a critical temperature that "does not coincide with the material's superconducting temperature" and show parity violation (mainly one rotation sense).
  - Phenomenology: B_g = γ B_g0 with B_g0 = 4Gmω/[c²(R_o + R_i)] at the ring centre (eq. 4), giving γ ≈ 1.2–1.7×10¹⁸.
  - Table 1 (B_g/ω ×10⁸, CW/CCW): YBCO+Al 5.3/−1.2; Nb+Al 3.2/−0.4; Al 3.8/−0.7; Graham et al. Pb −5.3 ± 8.5 / 37.7 ± 13.2; GP-B upper limit < 0.1.
  - GP-B "rules out all present models by up to 4 orders of magnitude" and "rules out our initial Cooper-pair mass anomaly hypothesis … by 5 orders of magnitude". Theory values (×10⁸): Tajmar–de Matos 395, Dröscher–Hauser 130, de Matos–Beck 1.6.
  - The accelerometer and gyro results disagree by a factor of 100 under the induction law.
  - Source: [arXiv 0707.3806](https://arxiv.org/abs/0707.3806)
- **2011 refutation.** M. Tajmar, "Evaluation of enhanced frame-dragging in the vicinity of a rotating niobium superconductor, liquid helium and a helium superfluid", Supercond. Sci. Technol. 24, 125011 (2011).
  - SRS-1000 gyro; no anomalies within 3× the ±5×10⁻⁸ rad/s noise; electric motor up to 5000 rpm.
  - Limits: 4×10⁻¹¹ for Nb or liquid He and 3×10⁻¹⁰ for superfluid He, "up to two orders of magnitude below any previous result".
  - Source: [IOP](https://iopscience.iop.org/article/10.1088/0953-2048/24/12/125011)
- The refined measurements "ruled out a gravity-like behavior and attributed the obtained results most likely to acoustic noise and vibration artefacts". — [Tajmar et al. 2022](https://www.frontiersin.org/articles/10.3389/fphy.2022.892215/full)
- In 2011 Tajmar re-interpreted the results because new equipment measured gyroscope signals about 2 orders of magnitude lower. — [arXiv 2206.07574](https://arxiv.org/abs/2206.07574)
- **Intermediate result (2010).** The rotating-liquid-helium fibre-optic gyroscope test gave an anomalous signal of similar magnitude but opposite sign, at resolution about 3× the gyro resolution. — [arXiv 0911.1033](https://arxiv.org/abs/0911.1033)
- **Graham et al. 2008 (Canterbury Ring Laser Group, Pb disc, UG2 ring laser).** "Results declared within 1σ from a null result; for the results of the Tajmar group the figure would be 11σ." They also reported a "parity violation". — [arXiv 2206.07574](https://arxiv.org/abs/2206.07574)
  - Tajmar's conversion of their CCW data used a ×1.7×10⁶ dipolar extrapolation and gives 3.8 ± 1.3 ×10⁻⁷, with γ ≈ 1.9×10¹⁸. — [arXiv 0707.3806](https://arxiv.org/abs/0707.3806)
- The 2022 review says that by GR the Tajmar-group field "would be equivalent to one produced by a white dwarf star". — [arXiv 2206.07574](https://arxiv.org/abs/2206.07574)
- **Kiefer & Weber's assessment.** It is "strongly questionable" whether such effects can be separated from systematics. Tate's analysis is very sensitive to the effective Cooper-pair mass, and the classical coupling gives B_g ≈ G(2m/e)²ω ≪ B (Gaussian units). — [Kiefer & Weber, arXiv gr-qc/0408010](https://arxiv.org/abs/gr-qc/0408010)

### Inferences
- **Computed: Tajmar's eq. (6) coefficient reproduces.** ρ*λ_L² = m_e²/(μ0e²) = 2.572×10⁻¹⁷ kg/m, and divided by g0 this is 2.623×10⁻¹⁸ (paper 2.624×10⁻¹⁸).
  - Nb ρ*/ρ = 3.94×10⁻⁶ requires λ_L = 27.6 nm. With λ_L = 39 nm or 41.6 nm it is 1.97×10⁻⁶ or 1.74×10⁻⁶.
  - The −0.7×10⁻⁸ s² prediction is reproduced (6.4×10⁻⁹) with ρ*/ρ = 1.74×10⁻⁶ at r = 3.6 cm.
  - The Tate value −37×10⁻⁸ is reproduced as 3.4×10⁻⁷.
  - So the two Tajmar papers silently use different λ_L for Nb: about 28 nm and about 42 nm.
- **Computed: measured k.** k = |coupling|·g0/r = 9.64×10⁻⁸ × 9.81/0.036 = 2.63×10⁻⁵, which matches their "Δm*/m* ≈ 2.6×10⁻⁵". Also 100 µg/1500 rad/s² = 6.67×10⁻⁸ s², as in the paper.
- **Computed GR prediction for the same ring.**
  - Nb ring OD 150, ID 138, height 15 mm, density 8570 kg/m³: M = 0.349 kg, R̄ = 72 mm.
  - Lense–Thirring Biot–Savart at the centre of a mass ring: Ω_drag = 2GMω/(c²R̄), so k_GR = 7.20×10⁻²⁷. This equals Tajmar's own B_g0/ω = 4Gm/[c²(R_o + R_i)].
  - Ratio k_meas/k_GR = 3.65×10²¹ (10^21.56). With the 2006 B_g = 2Ω convention it is 7.3×10²¹ (10^21.86).
  - Direct acceleration comparison at r = 3.6 cm and ω̇ = 1500 rad/s²: a_GR = r k_GR ω̇ = 3.9×10⁻²⁵ m/s², against 9.8×10⁻⁴ m/s² claimed, a ratio of 2.5×10²¹.
  - **The "≈ 22 orders" prior estimate is reproduced: 21.4–21.9 depending on convention.** The paper's "30 orders" and ESA's "10²⁰" are not reproducible from the stated setup.
- **Computed: what source would GR need?** k = 2GM/(c²R) = 2.6×10⁻⁵ requires compactness GM/(c²R) = 1.3×10⁻⁵. That is white-dwarf-like (Sirius B ≈ 2.5×10⁻⁴), consistent with the review's "white dwarf" phrase. Packed into the ring's 4.07×10⁻⁵ m³, the mass would be 1.3×10²¹ kg at a density of 3×10²⁵ kg/m³.
- **Computed: the 2007 gyro claims.** B_g/ω = 3.2×10⁻⁸ for Nb gives γ = 4.4×10¹⁸ against the bare Nb ring; the paper's 1.2×10¹⁸ includes the Al holder mass. That is about 18.6 orders above GR.
- **Computed: the 2011 limit.**
  - 4×10⁻¹¹ is 800× below the 2007 Nb gyro signal.
  - It is 9.9×10⁴ below the de Matos ρ*/ρ prediction (3.95×10⁻⁶) and 4.6×10⁶ below the Tate-derived 1.84×10⁻⁴.
  - It is still 5.6×10¹⁵ above GR. A laboratory frame-dragging measurement remains about 15.7 orders out of reach.
- **Physics objection (mainstream).** The Proca step is the weak point. Linearised GR has no graviton mass, and the observed Λ corresponds to λ_g ∼ Hubble scale. Choosing 1/λ_g² = μ₀g ρ_m is an ad hoc ansatz that makes the "mass" depend on local density. In GR the gravitomagnetic Larmor theorem (B_g = −2ω inside a rotating frame) is a coordinate effect, not a source-generated field, so it cannot leak outside a ring.
- **Calculator spec (Tajmar module).**
  - `gr_ring_k(M, Ro, Ri) = 4GM/(c²(Ro+Ri))`.
  - `gr_tangential_accel(r, k, omega_dot) = r·k·omega_dot`.
  - `tajmar_k(rho_star_over_rho)` returns k. The 2006 B_g = 2ω·k.
  - `rho_star(lambda_L) = m_e²/(μ0 e² λ_L²)`.
  - `coupling_s2(k, r) = k·r/g0`.
  - `ratio = k_claim/k_GR`.
  - Test vectors: (0.349 kg, 0.075, 0.069) → 7.20e-27; claim 2.63e-5 → ratio 3.65e21; λ_L = 27.6 nm, ρ = 8570 → ρ*/ρ = 3.94e-6.

### Gaps
- I did not retrieve Graham et al. 2008 (Physica C 468, 383?) directly. Its title, pages, disc mass and rpm, and the native null bound are unverified. A search snippet quoted "1.7 ± 3.8 ×10⁻⁷ assuming a dipolar distribution", which I could not confirm. The ring laser is UG2, 21 × 39.7 m (834 m²; Hurst et al., J. Appl. Phys. 2009) per search summaries only.
- The full text of the 2011 SST paper (vibration mechanism, exact definition of "coupling" for the 4×10⁻¹¹ limit) was not read.
- The STAIF journal citation for gr-qc/0610015 (AIP Conf. Proc. volume and pages) is unverified.

## Q4. Ning Li & Douglas Torr: gravitomagnetic fields in superconductors, Harris's critique, and the history

### Takeaway
Li & Torr's 1991 PRD equations are standard linearised-GEM bookkeeping: inside a pure superconductor B + (m/q)B_g → 0. Taken at face value they give tiny numbers: an induced B_g of 4×10⁻³² s⁻¹ per tesla, and a residual B of about 6×10⁻²⁶ T from Earth's B_g. Their 1992–1993 papers then claimed that lattice-ion spin alignment yields detectable gravitomagnetic and gravitoelectric fields. Harris (Found. Phys. Lett. 12, 201, 1999) showed the fields were "too large by many orders of magnitude because of unrealistic assumptions". Li's NASA collaborations produced only a null static test (< 2×10⁻⁸ g). Her company AC Gravity LLC took a $448,970 DoD grant (2001–2002) and published no results. Status: PRD 1991 formalism **mainstream but negligible**; enhanced-field claims **refuted (theoretically) and never demonstrated**.

### Cited Findings
- **PRD 43, 457 (15 Jan 1991), "Effects of a gravitomagnetic field on pure superconductors" (Li & Torr, UAH CSPAR).**
  - B(z) ≈ −(m²μ_g/q²μ) B₀ − (m/q) B_g,0.
  - B_g(z) ≈ B_g,0 + (μ_g m/qμ) B₀.
  - B + (m/q)B_g decays exponentially to zero over λ.
  - A small residual uniform B pervades the superconductor.
  - Source: [APS abstract](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.43.457). The equations were relayed by a fetch summary and should be confirmed against the PDF.
- The conclusion was that a pure superconductor always carries a residual B produced by Earth's gravitomagnetic field. — [Ruggiero & Tartaglia, arXiv gr-qc/0207065](https://arxiv.org/abs/gr-qc/0207065)
- **PRB 46, 5489 (1 Sep 1992), "Gravitational effects on the magnetic attenuation of superconductors".** A "σ–μ model" in which perfect conductivity and zero permeability become related once gravity is included. The abstract has no numbers. — [APS abstract](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.46.5489)
- **Torr & Li, Found. Phys. Lett. 6, 371–383 (1993), "Gravitoelectric-electric coupling via superconductivity".** They claim that "the carriers of quantized angular momentum in superconductors are not the Cooper pairs but the lattice ions". — [Springer](https://link.springer.com/article/10.1007/BF00665654) (abstract via search summary)
- **Harris, Found. Phys. Lett. 12(2), 201–208 (1999).** The fields Torr & Li calculate are overestimated "by many orders of magnitude" because of "unrealistic assumptions". — [INSPIRE](https://inspirehep.net/literature/1359113)
- Harris also wrote the standard GEM analogy paper: Am. J. Phys. 59, 421 (1991). — [Gallerati & Ummarino, arXiv 2203.09417 ref. 66](https://arxiv.org/abs/2203.09417)
- Modanese (1996) argued that Torr & Li's anomalously strong gravitoelectric field "contrasts with the full equations". — [arXiv hep-th/9505094](https://arxiv.org/abs/hep-th/9505094)
- **History.**
  - Li proposed that aligned rotating lattice ions in a high-T_c superconducting disc could produce a strong gravitomagnetic field and repulsive force, and sought NASA and Army Aviation and Missile Command funding.
  - She left UAH in 1999 to found AC Gravity LLC, which received a DoD grant of $448,970 in 2001. The grant ended in 2002 with no public results.
  - In 2003 she gave a talk, "Measurability of AC gravity fields", at a MITRE conference.
  - She was injured in 2014 and died on 27 July 2021.
  - Source: [Wikipedia: Ning Li (physicist)](https://en.wikipedia.org/wiki/Ning_Li_(physicist))
  - The grant was from the Army Aviation and Missile Command for a "Gravito-Electro Magnetic Superconductivity Experiment" that ended on 25 Sep 2002. — secondary Spanish-language article relayed by search ([Muy Interesante](https://muyinteresante.okdiario.com/?p=188419)); unverified
- NASA MSFC's static test with Li as lead author found < 2×10⁻⁸ (Physica C 281, 260, 1997). — [NTRS](https://ntrs.nasa.gov/citations/19990039542)

### Inferences
- **Computed: the PRD 1991 terms with μ_g = 4πG/c² = 9.33×10⁻²⁷ m/kg and m/q = m_e/e = 5.69×10⁻¹² kg/C.**
  - Induced B_g per tesla of applied field: μ_g(m/q)/μ0 = 4.2×10⁻³² s⁻¹ T⁻¹.
  - The −(m²μ_g/q²μ)B₀ term is about 2.4×10⁻⁴³ B₀.
  - Residual B = (m/q)B_g,Earth with B_g,Earth ≈ 1.7×10⁻¹⁴ rad/s (computed GJ/(c²R³) = 1.68×10⁻¹⁴) gives about 10⁻²⁵ T.
  - All of these are unmeasurable, and they are what the standard theory gives.
- The "enhanced" Li–Torr fields require replacing the Cooper-pair mass current by coherently aligned nuclear spins acting as gravitomagnetic dipoles. Even full nuclear-spin alignment gives an angular-momentum density of order nħ ≈ 10²⁹ m⁻³ × 10⁻³⁴ J s = 10⁻⁵ kg m⁻¹ s⁻¹. The GR gravitomagnetic field of that is about (G/c²)×10⁻⁵ ≈ 10⁻³² s⁻¹. This is my order-of-magnitude estimate, consistent with Harris's "many orders of magnitude".
- **Calculator spec (Li–Torr module).**
  - `li_torr_internal(B0, Bg0) -> (B, Bg)` using the PRD formulas with μ = μ0 and μ_g = 4πG/c².
  - Test vector: B0 = 1 T, Bg0 = 1e-14 s⁻¹ → B ≈ −5.7e-26 T (magnitude), Bg ≈ 1e-14 + 4.2e-32 s⁻¹.

### Gaps
- I could not read the full PRD/PRB/FPL texts, so Li–Torr's specific "predicted field strengths" (for example any lattice-ion enhancement factor) are unverified. The often-quoted large enhancement number was not found in any retrieved source.
- I could not read Harris 1999 beyond the abstract, so the exact number of orders of magnitude is unknown to me.
- I have no primary record of the NASA "Delta-G" programme beyond the NTRS abstracts above. The DoD contract number and final report status (reportedly never delivered) are unverified.

## Q5. Legitimate theory: DeWitt, the London moment and Cooper-pair mass, Kiefer & Weber, Ummarino & Gallerati

### Takeaway
- **Mainstream.** The London moment B = 2m_eω/e (1.137×10⁻¹¹ T per rad/s) is measured. DeWitt's minimal coupling puts the GR gravitomagnetic potential into the fluxoid, so the London moment measures rotation relative to local inertial frames. The GR correction for lab rings is a fraction of about 10⁻²⁶, and DeWitt's induced current is about 10⁻²⁸ A, far below the fA SQUID floor. The gravitational London penetration depth is about 10¹⁴ m, so there is no gravitomagnetic Meissner expulsion.
- **Tate anomaly.** Tate's 92 ppm excess (4.4σ) has been re-measured: Hoang et al. 2020 find an 18 ± 2.1 ppm gap for Nb. It remains a small open discrepancy in condensed-matter corrections, not evidence for gravity. Tajmar's own data exclude a gravitomagnetic explanation by 5 orders of magnitude.
- **Contested.** Ummarino & Gallerati's TDGL-plus-gravito-Maxwell papers predict only *transient* changes lasting picoseconds to nanoseconds: Eg/g* → 1 within about 4 ps for YBCO, and Δg ≈ 1×10⁻⁸ to 2×10⁻⁴ m/s² outside a 30 cm disc for 1.6 ns or 7.5 ps. Their framework treats the gravitoelectric field like an electric field that the condensate can screen, which conflicts with the equivalence principle and Schiff–Barnhill equilibrium.

### Cited Findings
- **DeWitt, "Superconductors and gravitational drag", PRL 16, 1092–1093 (1966).** — [INSPIRE](https://inspirehep.net/literature/50381)
  - The gravitomagnetic potential A_g produces a Berry-type phase ΔS ∼ ∮A_g·dl (eq. 43); the coupling is H_int ∼ A·A_g (eq. 44).
  - DeWitt's induced current in a ring near a rotating mass is I ∼ ε0c⁴κ m M v/(e d), with κ = 8πG/c⁴ (eq. 71). For M = 1 kg and v = 300×2π×d/2 this gives I ≈ 10⁻²⁹ A, "out of experimental reach" (SQUIDs reach about fA).
  - Tajmar–de Matos-type quantisation including A_g (Gaussian units): ∮(m v_s + (e/c)A − (2m/c)A_g)·dl = nh/2 (eq. 72). Thick ring: B = −(2mc/e)ω + (2m/e)B_g (eq. 74). Classical coupling B_g ≈ G(2m/e)²ω ≪ B (eq. 76).
  - Tajmar–de Matos's implied fields are B_g ≈ 7000 c/2 Hz (Hildebrandt's London-moment data) and 1.65×10⁻⁵ c/2 Hz (Tate). These are "strongly questionable".
  - Source: [Kiefer & Weber, Ann. Phys. (Leipzig) 14, 253 (2005), arXiv gr-qc/0408010](https://arxiv.org/abs/gr-qc/0408010)
- DeWitt predicted Meissner-like expulsion of the Lense–Thirring field. A 2022 paper notes technical shortcomings in his flux-quantisation interpretation and his current estimate. — search summary of [Front. Phys. 2022, 10.3389/fphy.2022.823592 (Inan)](https://www.frontiersin.org/articles/10.3389/fphy.2022.823592/full); not read in full
- **London moment.** B = −(2mc/e)ω (Gaussian) deep inside, using quantisation with n = 0. Thin ring: h/2m = 2SΔω gives the Cooper-pair mass m′ = 2m (eqs. 64–66). — [Kiefer & Weber](https://arxiv.org/abs/gr-qc/0408010)
- **Tate anomaly.** PRL 62, 845 (1989): (m′/2m_e) = 1.000084(21) against a theoretical 0.999992 (which includes relativistic corrections). — [APS](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.62.845)
  - Corrections were attempted by Liu (1998) and Jiang & Liu (2001), but "the discrepancy … remained unsolved". — [arXiv 0707.3806](https://arxiv.org/abs/0707.3806)
- **Hoang et al., Mater. Lett. 262, 127176 (2020), "High-precision measurement of Cooper-pair mass using rotating spherical-shell superconductor".** Up to six superconductors were measured; for Nb "the gap to predictions is now 18 ppm with an accuracy of only 2.1 ppm". A "simple explanation" is that the theory needs revision. — [Tajmar et al. 2022](https://www.frontiersin.org/articles/10.3389/fphy.2022.892215/full)
- Bulyzhenkov argues that flux quantisation cannot measure the inertial Cooper-pair mass in principle. — [arXiv 0907.0438](https://arxiv.org/abs/0907.0438) (via search summary)
- **Ummarino & Gallerati, EPJC 77, 549 (2017), "Superconductor in a weak static gravitational field" (arXiv 1710.01267v4).** — [arXiv](https://arxiv.org/abs/1710.01267)
  - **Gravito-Maxwell equations (eq. 34).** ∇·E_g = 4πG ρ_g; ∇·B_g = 0; ∇×E_g = −∂B_g/∂t; ∇×B_g = (4πG/c²) j_g + (1/c²)∂E_g/∂t. Force: dv/dt = E_g + v×B_g (eq. 37). Definitions: E_g = −½G₀₀ᵢ, A_g = ¼h̄₀ᵢ, B_g = ∇×A_g.
  - **Generalised fields (eq. 38).** E = E_e + (m/e)E_g, B = B_e + (m/e)B_g, φ = φ_e + (m/e)φ_g, A = A_e + (m/e)A_g. Vacuum constants: ε_g = (1/4πG)(e²/m²) and μ_g = (4πG/c²)(m²/e²) (eq. 41).
  - **Gravitational London depth.** λ_g = √(μ0λ_e²/μ_g) = √(c²/(4πG m n_s)) (eq. 59). Combined depth λ = λ_gλ_e/√(λ_g² + λ_e²) ≈ λ_e, with λ_g/λ_e ≈ 10²¹ (eqs. 61–62).
  - **TDGL (eq. 65).** (ħ²/2mD)(∂_t + 2ieφ/ħ)ψ − aψ + b|ψ|²ψ + (1/2m)(iħ∇ + 2eA/c)²ψ = 0 for ψ, plus the Maxwell equation for A. Here a = a₀(T − T_c) and b = b(T_c) (eq. 66).
  - **Dimensionless scales (eq. 68).** ξ = h/√(2m|a|); λ = √(bmc²/(4π|a|e²)); H_c; κ = λ/ξ; τ = λ²/D; η = 4πσD/(ε0c²).
  - **Gravity term.** φ = −g*x with g* = λκ m g/(√2 e H_c D) ≪ 1 (eqs. 72–73). Static solution ψ₀ = tanh(κx/√2)·tanh(κ(x − L)/√2) (eq. 77). General field E_g/g* = (1 − e^(−Pt) − …) with P = |ψ₀|²/η (eqs. 86–88).
  - **Simplified case J = 0 (eq. 93).** E_g/g* = 1 − exp(−|ψ₀|² t/η), with t in units of τ. The effect is larger for small η (high normal-state resistivity, small D = v_F ℓ/3), i.e. in "bad", impure samples.
  - **Table 1.**

    | Quantity | YBCO | Pb |
    |---|---|---|
    | T_c | 89 K | 7.2 K |
    | T* | 77 K | 6.3 K |
    | ξ | 3.6e-9 m | 1.7e-7 m |
    | λ | 3.3e-7 m | 7.8e-8 m |
    | σ⁻¹ | 4e-7 Ω m | 2.5e-9 Ω m |
    | H_c | 0.2 T | 0.018 T |
    | κ | 94.4 | 0.48 |
    | τ | 3.4e-10 s | 6.1e-15 s |
    | η | 1.27e-2 | 6.6e3 |
    | D | 3.2e-4 m²/s | 1 m²/s |
    | ℓ | 6e-9 m | 1.7e-6 m |
    | v_F | 1.6e5 m/s | 1.83e6 m/s |

  - **Table 2 (g*).** YBCO: 2.6e-12 (0 K), 9.8e-12 (70 K), 2e-11 (77 K), 2.8e-7 (87 K). Pb: 1e-17 (0 K), 1.4e-17 (4.2 K), 8.2e-17 (6.26 K), 2.2e-15 (7.10 K). For Pb with L = 6λ the system relaxes after about 10⁵ τ.
  - Conclusion: "the reduction of the gravitational field in a superconductor, if it exists, is a transient phenomenon". They suggest pulsed magnetic fields to switch superconductivity on the t_int timescale.
- **Follow-ups**, summarised in [Gallerati & Ummarino, Symmetry 14, 554 (2022), arXiv 2203.09417](https://arxiv.org/abs/2203.09417) and [Gallerati, Modanese & Ummarino review, arXiv 2206.07574](https://arxiv.org/abs/2206.07574).
  - **Josephson AC from gravity (CQG 37, 217001, 2020; arXiv 2009.04967).** ΔV_g = (m/e)gℓ and ω = 2eΔV/ħ = 2mgℓ/ħ (eq. 25 of the 2206 review). For ℓ = 1 nm, ω ≈ 1.7×10⁻⁴ s⁻¹ and T ≈ 3.7×10⁴ s. For an LTSC with ℓ = 300 nm, V_g ≈ 1.67×10⁻¹⁷ V, ω ≈ 0.05 s⁻¹ and T ≈ 123 s.
  - **Outside the sample, fluctuation regime (Symmetry 11, 1341, 2019; arXiv 1910.13897).** Disc with R = 15 cm and h = 3 cm, measured at d = 0.1 cm.
    - In (ξ₀ = 360 nm, T_c = 3.410 K, ΔT = 10⁻³ K): peak Δ = 1.24×10⁻⁸ m/s² at τ₀ = 1.64 ns.
    - Ba₀.₄K₀.₆Fe₂As₂ (ξ₀ = 1.20 nm, T_c = 37.0 K, ΔT = 0.1 K): Δ = 1.74×10⁻⁴ m/s² at τ₀ = 7.5×10⁻³ ns.
    - Scalings: Δ ∝ ξ(T)⁻¹ (eq. 96) and τ₀ ∝ (T − T_c)⁻¹ (eq. 97).
  - **Inside the sample (Entropy 23, 193, 2021; arXiv 2102.01489).** E_g/g* ≈ 1 − t/η − (1/η)∂γ₂/∂x at t → 0⁺ (eq. 134), with τ = λ²/D and λ(T) ≈ λ₀/√((T_c − T)/T_c). For Bi₂Sr₂CaCu₂O₈ (T_c ≈ 109 K, λ₀ ≈ 500 nm, σ⁻¹ ≈ 3.6×10⁻⁶ Ω m, D ≈ 10⁻³ m²/s) at 105 K: τ ≈ 6.8×10⁻⁹ and η ≈ 3.5×10⁻⁴, giving "a reduction of the local gravitational field of the order of 2×10⁻⁵".
  - **Vortex lattice with EM fields (arXiv 2110.07335).** E_z = 4πB_c(T)D/(ηL) ≈ 77 V/m for Bi₂Sr₂Ca₃Cu₃O₁₀ (T_c ≈ 107 K, L = 15 cm, T = 102 K, B_c ≈ 0.32 T), plus a transient E_g/g dip over less than 1 τ (eq. 165, Fig. 9).
- **Critique of the Ummarino–Gallerati GEM conventions.** Behera, EPJC 77, 822 (2017). — [arXiv 1709.04352](https://arxiv.org/abs/1709.04352)
- **Gravity-induced fields in metals (mainstream).** Experiments include Witteborn–Fairbank, Lockhart et al. PRL 38, 1220 (1977), and Rieger PRB 2, 825 (1970) on gravitationally induced electric fields in metals. — [arXiv 2203.09417 refs. 59–61](https://arxiv.org/abs/2203.09417)

### Inferences
- **Computed: Ummarino–Gallerati numbers reproduce.**
  - g*(YBCO, 77 K) = 1.92×10⁻¹¹ and g*(Pb, 6.26 K) = 8.20×10⁻¹⁷, using m = m_e and e (not 2m_e, 2e). τ_YBCO = 3.40×10⁻¹⁰ s.
  - η_YBCO = 4πμ0σD = 1.263×10⁻² and η_Pb = 6.3×10³ (paper 6.6×10³; the difference is resistivity rounding).
  - **Unit inconsistency.** The 2017 paper uses η = 4πμ0σD, while the 2021–22 papers use η = μ0σD: the Bi2212 value 3.49×10⁻⁴ reproduces only without the 4π. This is a factor-4π convention change between papers.
  - Bi2212 τ = 6.81×10⁻⁹ s reproduces. The vortex-lattice E_z = 76.8 V/m reproduces, and it equals 4πB_c/(μ0σL), independent of D.
  - Josephson: ω(1 nm) = 1.69×10⁻⁴ s⁻¹ with T = 3.71×10⁴ s; ω(300 nm) = 0.0508 s⁻¹ with T = 123.6 s; V_g = 1.673×10⁻¹⁷ V. All reproduce, using m = m_e and Cooper-pair charge 2e in ω = 2eΔV/ħ.
- **Computed: the transient lifetime in eq. (93).** For YBCO at 77 K, the e-folding time is ητ = 4.3×10⁻¹² s, so any "reduction" decays within about 4 ps of the transition. The claimed outside-sample Δ are 1.2×10⁻⁸ m/s² (≈ 1.3×10⁻⁹ g) for 1.6 ns and 1.7×10⁻⁴ m/s² (≈ 1.8×10⁻⁵ g) for 7.5 ps. A static-weighing experiment (Podkletnov-style, seconds to minutes) predicts **zero** steady weight change, so these papers do not support the Podkletnov claims. The predicted signals are too brief for any gravimeter: superconducting gravimeters average over seconds.
- **Mainstream objection to the U–G Josephson and "screening" results (my inference).** In equilibrium the electrochemical potential is uniform. Gravity on electrons in a metal is balanced by a tiny electrostatic field E = m_e g/e = 5.6×10⁻¹¹ V/m (computed; Schiff & Barnhill, Phys. Rev. 151, 1067 (1966), [doi](https://doi.org/10.1103/PhysRev.151.1067)). A static gravitational potential difference across a junction therefore does not produce a voltage-driven AC Josephson current. Combining (m/e)E_g with E_e into one "generalised field" also makes gravity screenable like electricity, which violates universality of free fall, because the factor m/e differs between species.
- **Computed: legitimate magnitudes for calculators.**
  - London moment: B/ω = 2m_e/e = 1.137×10⁻¹¹ T s; at 650 rad/s, B = 7.39 nT.
  - GR frame-drag correction to B for the Nb ring at 650 rad/s: (2m_e/e) k_GR ω = 5.3×10⁻³⁵ T, a fraction of 7.2×10⁻²⁷. This is the GR "gravitomagnetic London moment" and is 22 orders below Tate's 9.2×10⁻⁵.
  - Tate excess: 92 ppm at 21 ppm, i.e. 4.4σ. Hoang 2020 Nb: 18 ppm at 2.1 ppm, i.e. 8.6σ.
  - DeWitt current with my constants: I = 8πε0Gm_eMv/(ed) = 8.0×10⁻²⁹ A for M = 1 kg and v/d = 300π s⁻¹ (Kiefer & Weber: "≈ 10⁻²⁹ A"; same order).
  - Nb λ_L = 39 nm gives n_s = 9.3×10²⁷ m⁻³, ρ* = 0.0169 kg/m³ and λ_g = √(c²/(4πGρ*)) = 8.0×10¹³ m. Then λ_g/λ_L = 2.0×10²¹, which matches Ummarino's "≈ 10²¹".
  - For neutron-star density (10¹⁷ kg/m³) this formula gives 33 km, against "≈ 12 km" quoted by Modanese from Lano. The definitions differ by a numerical factor (√(16πG) vs √(4πG) conventions).
- **Calculator spec (legitimate-theory module).**
  - `london_moment(omega) = 2 m_e ω/e`.
  - `london_moment_gr(omega, Omega_drag) = 2 m_e (ω − Ω_drag)/e`. This is my reading of DeWitt minimal coupling: rotation relative to local inertial frames. Sign conventions vary by author.
  - `cooper_mass_from_tate(S, dnu) = h/(2 S dnu)`.
  - `lambda_g(rho_star) = c/√(4πGρ*)`.
  - `dewitt_current(M, v, d)`.
  - `ug_gstar(lambda, kappa, Hc, D)`.
  - `ug_eta(sigma, D, convention)`.
  - `ug_tau(lambda, D)`.
  - `ug_eg_ratio_simplified(psi0, t, eta) = 1 − exp(−psi0² t/eta)`.
  - `ug_josephson(ell) -> (omega, period, Vg)`.
  - `schiff_barnhill() = m_e g/e`.

### Gaps
- I did not read the full DeWitt 1966 letter. Its exact formula, the sign of the drag term and the Meissner-expulsion claim are relayed through Kiefer & Weber and a 2022 Frontiers paper.
- No retrieved source gives a definitive *resolution* of the Tate anomaly. The 2020 Hoang et al. result shrinks it but leaves an 18 ppm gap. Theoretical corrections (Cabrera & Peskin 1989; Liu 1998; Jiang & Liu 2001; Capelle & Gross 1999) were not retrieved.
- I did not open the original U–G 2019/2020/2021 papers; their numbers are taken from the authors' 2022 reviews. The exact meaning of "reduction of the order of 2×10⁻⁵" (the fraction of E_g/g* and the time window) is unclear from the review text.
- Primary material on Graham et al. and Unnikrishnan is missing (see Q1 and Q3).
