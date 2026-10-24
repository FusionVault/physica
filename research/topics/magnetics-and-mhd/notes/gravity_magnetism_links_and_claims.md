# Gravity and magnetism: GEM, frontier couplings and fringe magnetic claims, each as formula + test value + bound

Every number marked "script" is recomputed by `gravity_magnetism_checks.py` (in this folder). It uses the standard library only, and its asserts pass. It was run with `python3 -I` from the session scratchpad, outside this folder. Labels: **[mainstream]** is established physics, **[contested]** is published but disputed or unreplicated, **[refuted]** contradicts a measurement or a conservation law. SI is used throughout. Gaussian and Heaviside–Lorentz (HL) forms are given wherever a source uses them.

This file does not redo work covered elsewhere. That includes the London moment and Tajmar, Ning Li, Podkletnov and the Li–Baker HFGW claims (`../../fringe-physics-and-dirds/notes/superconductor_gravity_claims_math.md`, `baker_li_hfgw_and_pais_patents.md`). It also includes Lense–Thirring and Gravity Probe B (`../../gravity-anomalies/notes/lab_gravity_manipulation.md`), Gertsenshtein conversion and HFGW (`../../gravitational-waves/notes/high_frequency_gravitational_waves.md`), Schwinger and birefringence (`../../coulomb-and-schwinger/notes/schwinger_limit_strong_field_qed.md`), and EmDrive, Mach-effect, Biefeld–Brown and the Bennett pinch (other notes in the tree).

## Q1. Gravitoelectromagnetism (GEM): conventions, the factors of 2 and 4, a rotating body's field, the charges and the moment, and a Maxwell-to-GEM table with test values

### Takeaway
[mainstream] All GEM conventions are the same physics: linearised GR with h̄₀₀ = 4Φ/c² and h̄₀ᵢ = −2Aᵢ/c². The factor 2 comes from spin-2 gravity, which makes the gravitomagnetic charge twice the gravitoelectric charge. Different authors put that 2 in the field definition, in the Ampère law, or in the force law. If you rescale B so that the field equations look exactly like Maxwell's, a 4 appears in the force law. The quantities that do not depend on convention are the gyroscope precession Ω_P and the force 2m(v/c)×B (Mashhoon). These are what LAGEOS/LARES (61.3 mas/yr, 2026) and GP-B (39.2 mas/yr) test. In the script, all three conventions give the same force to 10⁻¹².

### Cited Findings
- **Mashhoon convention** [mainstream] — [Mashhoon, gr-qc/0311030](https://arxiv.org/abs/gr-qc/0311030):
  - Potentials: h̄₀₀ = 4Φ/c² and h̄₀ᵢ = −2Aᵢ/c². The metric is ds² = −c²(1 − 2Φ/c²)dt² − (4/c)(A·dx)dt + (1 + 2Φ/c²)δᵢⱼdxⁱdxʲ.
  - Far field of a rotating body: Φ ≈ GM/r and A ≈ (G/c) J×x/r³.
  - Field definitions: E = −∇Φ − (1/c)∂(A/2)/∂t and B = ∇×A. Both "have dimensions of acceleration".
  - Field equations: ∇·E = 4πGρ, ∇×(B/2) = (1/c)∂E/∂t + (4πG/c)j, ∇×E = −(1/c)∂(B/2)/∂t, ∇·(B/2) = 0.
  - Force law: F = −mE − 2m(v/c)×B.
  - Charges: a test particle has q_E = −m and q_B = −2m. A source has Q_E = M and Q_B = 2M. The ratio "is always 2, since linearized gravity is a spin-2 field"; for Maxwell's spin-1 field it is unity.
- **Mashhoon: moments, precession and energy** [mainstream] — [Mashhoon, gr-qc/0311030](https://arxiv.org/abs/gr-qc/0311030):
  - A test gyroscope has gravitomagnetic dipole moment μ_g = −S/c. A source has gravitomagnetic dipole moment J/c.
  - Precession: Ω_P = (GJ/c²r³)[3(Ĵ·r̂)r̂ − Ĵ] = B/c. The interaction energy is S·Ω_P.
  - Gravitational Larmor theorem: Ω_L = q_B B/(2mc).
  - GEM Poynting vector: S = −(c/2πG) E×B.
  - Spin energy scales: ħΩ_⊕ ≈ 10⁻¹⁹ eV and ħΩ_P ≈ 10⁻²⁹ eV at Earth, ≈ 10⁻²⁷ eV near Jupiter.
- **Ruggiero & Tartaglia convention** [mainstream] — [Ruggiero & Tartaglia, gr-qc/0207065](https://arxiv.org/abs/gr-qc/0207065):
  - Potentials: h⁰⁰ = 4Φ/c², h⁰ˡ = −2Aˡ/c², with **Φ = −GM/r** (the opposite sign to Mashhoon) and Aₗ = (G/c)ε_nkl Sⁿxᵏ/r³.
  - Fields: E_G = −∇Φ − (1/2c)∂A/∂t and B_G = ∇∧A (eq. 15).
  - Field equations: ∇·E_G = −4πGρ and ∇∧(B_G/2) = (1/c)∂E_G/∂t − (4πG/c)j.
  - Dipole field: B_G = −4(G/c)[3r(r·S) − Sr²]/(2r⁵) (eq. 24).
  - Geodesic: dv/dt = E_G + (v/c)∧B_G (eq. 26).
  - Their footnote: "the effective gravitomagnetic charge is twice the gravitoelectric one; this is a remnant of the fact that the linear approximation of GR involves a spin-2 field".
- **"Maxwell-like SI" convention** (pointer only; it cites Mashhoon et al. 2001 *LNP* 562 and Clark & Tucker 2000 *CQG* 17) — [Wikipedia: Gravitoelectromagnetism](https://en.wikipedia.org/wiki/Gravitoelectromagnetism):
  - Field equations: ∇·E_g = −4πGρ_g and ∇×B_g − (1/c²)∂E_g/∂t = −(4πG/c²)J_g.
  - Force law: F = m(E_g + v×4B_g).
  - Rotating body: B_g = (G/2c²)[L − 3(L·r̂)r̂]/r³.
  - Its Earth formula, B_g = 2πR g₀/(5c²T), assumes a uniform sphere.
  - The page states: "There is no scaling choice that allows all the GEM and EM equations to be perfectly analogous".
- [mainstream] LARES 2 first results (2023, *EPJ Plus*) are "in complete agreement" with GR — [MDSOAR record](https://mdsoar.org/items/3a5ebd8b-d10d-491b-b7d4-f7fd31480736).
- [mainstream] *Nature* 655, 332–335 (2026) is the combined LARES 2 + LAGEOS measurement — [Sapienza IRIS](https://iris.uniroma1.it/handle/11573/1771060).
  - The measured node drift is 61.3 mas/yr, matching GR — [Physics World](https://physicsworld.com/a/frame-dragging-measurement-around-the-earth-sets-new-precision-record/).
  - The data are about 200,000 laser ranges from July 2022 to June 2025 — [The Debrief](https://thedebrief.org/scientists-achieve-most-precise-measurement-yet-of-einsteins-frame-dragging-effect/).
  - The stated uncertainty conflicts between sources: about 1 part in 10³ per [Scientific American](https://www.scientificamerican.com/article/einsteins-greatest-theory-triumphs-again-in-landmark-frame-dragging-measurement/) and Physics World, and 0.2% per The Debrief. Earlier tests were "several percent ... as much as 10 percent".
- [mainstream] GP-B frame dragging measured −37.2 ± 7.2 mas/yr, against GR's −39.2 — [Everitt et al. 2011, arXiv:1105.3456](https://arxiv.org/abs/1105.3456).

### Inferences
- **Earth test values (script).**
  - J_⊕ = 5.86×10³³ kg m²/s, assuming C = 0.3307 MR². This matches 5.87×10³³ in [Biemond's table](https://vixra.org/pdf/2006.0221v1.pdf).
  - Mashhoon B: 5.03×10⁻⁶ m/s² at the equator and 1.01×10⁻⁵ m/s² at the pole.
  - Ω_P = B/c = 1.68×10⁻¹⁴ rad/s at the equator.
  - Maxwell-like B_g: 8.38×10⁻¹⁵ s⁻¹ at the equator.
  - The force per mass on a body moving at 7.9 km/s is 2.65×10⁻¹⁰ m/s² in all three conventions, equal to 10⁻¹².
- **Factor-of-2 trap in Ruggiero & Tartaglia.** Their eq. (24) is −2 times the curl of their own eq. (10) potential. It is consistent with their force law (26) but not with their definition (15). Use (24) and (26) together.
- **Uniform-sphere trap.** The Wikipedia Earth formula overstates the exact B_g by 0.4/0.3307 = 1.21 (script).
- **Permeability: the factor-of-4 trap.**
  - μ_g = 4πG/c² = 9.33×10⁻²⁷ m/kg goes with "force ×4". This is the value Li & Torr use.
  - If B is rescaled so the force has factor 1, μ_g becomes 16πG/c² = 3.73×10⁻²⁶ m/kg.
  - ε_g = 1/(4πG) = 1.19×10⁹ kg s²/m³, and 1/√(ε_g μ_g) = c exactly (script).
- **Lense–Thirring node rates, Ω̇ = 2GJ/[c²a³(1−e²)^{3/2}] (script).**
  - LAGEOS 30.7 mas/yr, LAGEOS 2 31.5, LARES 118, LARES 2 30.7 mas/yr.
  - LAGEOS + LARES 2 = 61.35 mas/yr, against the published 61.3 (orbital elements assumed; see Gaps).
- **GP-B check (script).** The polar-orbit average GJ/(2c²a³) at a = R_⊕ + 642 km, projected onto the guide-star declination, gives 39.17 mas/yr. That matches GR's 39.2. The measurement sits −0.28σ from GR.
- **Spin coupling.** ħΩ_P at the pole is 2.2×10⁻²⁹ eV (script), matching Mashhoon's ~10⁻²⁹ eV.
- **Recommended crate design.** Store one internal convention (Mashhoon) and convert explicitly at the API boundary, e.g. `gem::Convention::{Mashhoon, RuggieroTartaglia, MaxwellSI}`. Test every function against Ω_P, which is the same in all conventions. This extends `gravitics::superconductors::ring_frame_dragging`.

**Maxwell → GEM map** (test values from the script)

| Maxwell (SI) | GEM, Mashhoon | GEM, Maxwell-like SI | Test value |
|---|---|---|---|
| charge q | q_E = −m, q_B = −2m (source +M, +2M) | m | ratio q_B/q_E = 2 (spin 2) |
| ε₀ | 1/(4πG) (attractive sign handled by q_E = −m) | −1/(4πG) | 1.19×10⁹ kg s² m⁻³ |
| μ₀ | 4πG/c² with ∇×(B/2) | 4πG/c² (force ×4), or 16πG/c² if rescaled | 9.33×10⁻²⁷ / 3.73×10⁻²⁶ m/kg |
| E | E = −∇Φ − (1/c)∂(A/2)/∂t, Φ = +GM/r | E_g = g | 9.81 m/s² at Earth's surface |
| B | B = ∇×A (m/s²) | B_g (s⁻¹) | Earth equator 5.03×10⁻⁶ m/s² ↔ 8.38×10⁻¹⁵ s⁻¹ |
| qv×B | −2m(v/c)×B | 4m v×B_g | 2.65×10⁻¹⁰ m/s² at 7.9 km/s |
| dipole moment μ | source J/c; gyroscope μ_g = −S/c | L, with dipole prefactor G/2c² | J_⊕ = 5.86×10³³ kg m²/s |
| Larmor ω = qB/2m | Ω_L = q_B B/(2mc) = −B/c | — | 1.68×10⁻¹⁴ rad/s at the equator |
| spin precession | Ω_P = B/c = (GJ/c²r³)[3(Ĵ·r̂)r̂ − Ĵ] | — | GP-B 39.2 mas/yr; LAGEOS node 30.7 mas/yr |
| Poynting E×H | S = −(c/2πG)E×B | S_g = −(c²/4πG)E_g×4B_g | — |
| spin energy −μ·B | S·Ω_P | — | ħΩ_P ≈ 2×10⁻²⁹ eV at the pole |

### Gaps
- The Ruggiero & Tartaglia journal version was not checked. The internal factor −2 between their eqs. (10) and (24) is my reading of the arXiv text.
- The orbital elements are from memory, not re-fetched: LAGEOS 12270 km and e = 0.0045; LAGEOS 2 12163 km and 0.014; LARES 7821 km and 0.0008; LARES 2 12266 km and about 10⁻⁴. The IM Peg declination of 16.84° is also from memory. The 61.3 mas/yr and 39.2 mas/yr matches support them.
- The *Nature* 2026 abstract was not read, so the 0.1% versus 0.2% uncertainty is unresolved. Clark & Tucker 2000 was not read.

## Q2. Blackett's 1947 hypothesis and the Schuster–Wilson hypothesis: β for the Earth, Sun, 78 Vir, pulsars and planets, and the 1952 null experiment

### Takeaway
[refuted] Blackett's hypothesis, μ = β√G J/(2c) (Gaussian), equivalent to **μ = (β/2)√(4πε₀G) J in SI**, gives β ≈ 0.3 for the Earth and ≈ 1.2 for 78 Vir. Across the bodies now measured, however, β ranges over about 5×10⁶. Mars's missing dipole forces β ≲ 2×10⁻⁵, against 0.31 for the Earth. Blackett's own 1952 magnetometer experiment was negative, and he abandoned the idea in the 1950s. The surviving form is Barrow & Gibbons's conjecture that β ≤ O(1). It holds for charged rotating black holes.

### Cited Findings
- [refuted] Schuster (*Proc. Phys. Soc.* 24, 121, 1912) proposed that rotating masses are magnets. Wilson (*Proc. R. Soc. A* 104, 451, 1923) "found no such effect". Blackett revived the idea in *Nature* 159, 658 (1947). "By the 1950s even Blackett felt it had been refuted" (Lovell 1975) — [Wikipedia: Blackett effect](https://en.wikipedia.org/wiki/Blackett_effect).
- [refuted] Wilson's residual charge was Q* = β G^½ m. He tested a swinging (translating) bar and saw no field. Blackett's β was 0.3, 1.14 and 1.16 for the Earth, Sun and 78 Vir. His negative experiment used a 10×10 cm gold cylinder at rest, co-rotating with the Earth — [Biemond, viXra 2006.0221 (not peer-reviewed; compilation)](https://vixra.org/pdf/2006.0221v1.pdf).
- [refuted] Blackett 1952, "A negative experiment relating to magnetism and the Earth's rotation", *Phil. Trans. R. Soc. A* 245, 309–370. Its "highly sensitive magnetometer ... gave a negative result and thus disproved the general validity of the earlier theory" — [Royal Society Making Science (referee reports)](https://makingscience.royalsociety.org/items/rr_77_79/referees-report-by-edward-crisp-bullard-on-a-paper-a-negative-experiment-relating-to-magnetism-and-the-earths-rotation-by-patrick-maynard-stuart-blackett).
- [refuted] Campanelli's convention is μ = β√G J/(2c), Gaussian.
  - Sirag (*Nature* 278, 535, 1979) extended the fit to Mercury, Venus, Jupiter, Saturn, the Moon and Her X-1.
  - Mars's dipole upper limit μ ≲ 2×10²⁰ G cm³ gives **β ≲ 2×10⁻⁵**, "the strongest constraint". Inputs: M = 6.4×10²³ kg, R = 3390 km, k = 0.925, P = 1.03 d.
  - The Barnett effect always dominates in the laboratory.
  - B_max ~ 10⁻¹³ G (β/10⁻⁵)(σ_max/1 MPa)^½(ρ/1 g cm⁻³)^½(R/1 m), so a laboratory test is not feasible.
  - — [Campanelli, arXiv:2010.02734](https://arxiv.org/abs/2010.02734)
- [contested] Barrow & Gibbons write μ/J = β√G/c, with no factor 2. For electrons β = N^½, where N = e²/(Gm²).
  - Kerr–Newman black holes have g = 2 and β < 1. Kerr–Kaluza–Klein holes have β = v/c ≤ 1, and Kerr–Sen holes β ≤ 2.
  - They conjecture that the Schuster–Wilson–Blackett number is bounded above by a number of order unity.
  - — [Barrow & Gibbons, PRD 95, 064040 / arXiv:1701.06343](https://arxiv.org/abs/1701.06343)
- [contested] β* values in Biemond's Table 1:
  - Moon 0.22, Mercury 0.86, Venus 0.038, Earth 0.31, Mars 0.0029 (crustal), Jupiter 0.073, Saturn 0.014.
  - 96 isolated pulsars 0.061 (mean), 3 binary millisecond pulsars 6.4×10⁻⁶, 10 isolated white dwarfs 32.
  - Sun 0.2 with B_p = 8 G; 78 Vir 1.16 with 1.5 kG.
  - Surdin's rotating brass and tungsten cylinders gave mean-square fields "in reasonable agreement" with the relation, but with no sign and with fluctuations.
  - — [Biemond](https://vixra.org/pdf/2006.0221v1.pdf)

### Inferences
- **SI form.** μ = (β/2)√(4πε₀G) J, with √(4πε₀G) = 8.62×10⁻¹¹ C/kg (script). This is the Reissner–Nordström extremal charge-to-mass ratio (Q8). Blackett's body therefore acts like mass carrying a "gravitational charge" β√(4πε₀G)M.
- **β recomputed from Biemond's (μ, J) (script).**
  - Earth 0.313, Sun 0.200, 78 Vir 1.16, Jupiter 0.073, Saturn 0.014, Venus 0.038, Moon 0.22.
  - Pulsars 0.060, millisecond pulsars 6.5×10⁻⁶, white dwarfs 31.
  - The spread max/min is 4.7×10⁶. The Mars dipole limit gives β < 2.4×10⁻⁵, against the Earth's 0.31, a factor of 1.3×10⁴. A universal β is refuted.
- **Convention trap.** Blackett/Campanelli β = 2 × Barrow–Gibbons β. For the electron β = 4.1×10²¹ in the Blackett convention and 2.0×10²¹ in Barrow–Gibbons (= N^½) (script).
- **Laboratory ceiling.** A tungsten sphere of R = 1 m at its burst speed, with β = 10⁻⁵, gives 1.7×10⁻¹⁵ T (script). Barnett magnetisation swamps this.
- **Calculator:** `blackett_moment(J, beta, Convention)` and `blackett_beta(mu, J)`. Use the Mars bound as the "mainstream bound to beat".

### Gaps
- The 1947 and 1952 full texts were not read. I did not find the 1952 predicted field and magnetometer sensitivity.
- The Sun's β differs between Blackett (1.14, from the 1940s "general field") and Biemond (0.2 at 8 G). This was not resolved.
- Babcock's 78 Vir field reached me only through secondary compilations.

## Q3. Magnetic monopoles: Dirac quantisation, the mass scale, the Parker bound, MoEDAL/ATLAS/IceCube, Cabrera 1982, and spin ice

### Takeaway
[mainstream theory, no detection]
- The Dirac charge is g_D = e/(2α) ≈ 68.5e in Gaussian units, or **h/(μ₀e) = 3.29×10⁻⁹ A·m in SI**. Its flux is μ₀g_D = h/e = 2Φ₀ = 4.14×10⁻¹⁵ Wb.
- GUT monopoles weigh about 10¹⁷–10¹⁸ GeV.
- The Parker bound is F < 10⁻¹⁵ cm⁻² s⁻¹ sr⁻¹. Neutrino telescopes reach 2×10⁻¹⁹ for relativistic monopoles.
- Collider searches exclude Drell–Yan monopoles up to about 4 TeV. MoEDAL's Schwinger searches reach about 80 GeV with no perturbative assumption.
- Cabrera's 1982 single event implies 6.1×10⁻¹⁰, about 6×10⁵ times the Parker bound [refuted as a flux].
- Spin-ice "monopoles" are real quasiparticles [mainstream]. Their charge is about 1/7700 of g_D.

### Cited Findings
- [mainstream] "magnetic charge g to be an integer multiple of the Dirac charge g_D ≃ 68.5e" — [Zhang et al. 2026, arXiv:2605.21281](https://arxiv.org/abs/2605.21281).
  - That paper builds a self-consistent Parker bound from the galactic dynamo's lowest eigenmode. It "substantially modif[ies] the standard extended Parker bound at low and intermediate masses".
- [mainstream] For a unification scale of 10¹⁶ GeV, GUT monopoles have M ~ 10¹⁷–10¹⁸ GeV — [PDG 2025 review: Magnetic Monopoles](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-mag-monopole-searches.pdf).
  - Parker bound (eq. 94.9): F < 10⁻¹⁵ cm⁻² sr⁻¹ s⁻¹ for M ≲ 10¹⁷ GeV, and 10⁻¹⁵(M/10¹⁷ GeV) above that.
  - Extended seed-field bound (eq. 94.10): F < (M/10¹⁷ GeV + 3×10⁻⁶)×10⁻¹⁶.
  - Galaxy-cluster arguments give a bound about 3 orders lower, "although less secure".
  - Collider "charge-dependent mass limits extend up to around 4 TeV".
- [mainstream] PDG 2025 listings — [PDG Listings](https://pdg.lbl.gov/2025/listings/rpp2025-list-magnetic-monopole.pdf):
  - Flux limits (cm⁻² sr⁻¹ s⁻¹):
    - IceCube (Abbasi 2022): < 2×10⁻¹⁹ for 0.86 < β < 0.995.
    - ANTARES (2017): < 1.5×10⁻¹⁸ for β > 0.6.
    - MACRO (2002): < 1.4×10⁻¹⁶ for 1.1×10⁻⁴ < β < 1.
    - Kobayashi 2023 Parker-type bound: < 3×10⁻²¹ for M < 10¹⁸ GeV.
    - Adams 1993 seed bound: 10⁻¹⁶(M/10¹⁷ GeV).
    - Parker 1970: < 10⁻¹⁶ as listed.
  - The **Cabrera 1982** candidate gives a flux of 6×10⁻¹⁰ and has "single Dirac charge within ±5%". Cabrera 1983 gives < 4×10⁻¹¹.
  - Collider cross-sections:
    - ATLAS (Aad 2023): σ < 4×10⁻⁴¹ cm² for 200–4000 GeV and 1–2 g_D.
    - MoEDAL Schwinger searches: Pb–Pb at 5.02 TeV covered masses 0–70 GeV and charges 1–3 g_D (2022). The 2024 search covered 0–80 GeV and charges 1–50 g_D.
    - MoEDAL dyons: mass limits 870–2040 GeV.
- [mainstream] Cabrera 1982 used "a 20-cm²-area superconducting loop". It recorded "a single candidate event" in "five runs totaling 151 days" and set an "upper limit of 6.1×10⁻¹⁰ cm⁻² sec⁻¹ sr⁻¹" — [Cabrera, PRL 48, 1378](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.48.1378).
- [mainstream] Newer flux searches:
  - NOvA (2025): 2×10⁻¹⁶ cm⁻² s⁻¹ sr⁻¹ for M > 10¹³ GeV and 0.005 < β < 0.8 — [arXiv:2512.20294](https://arxiv.org/abs/2512.20294).
  - ANTARES complete dataset: about 10⁻¹⁸ for relativistic monopoles — [arXiv:2505.23929](https://arxiv.org/pdf/2505.23929).
- [mainstream] Spin ice (Castelnovo, Moessner & Sondhi, *Nature* 451, 42, 2008):
  - The charges are Q = ±q_m = ±2μ/a_d, with μ ≈ 10 μ_B and a_d = √(3/2) a, where a ≃ 3.54 Å.
  - The ratio to the Dirac charge is about 8000.
  - Their proposed detection is a Stanford-type (Cabrera) loop experiment.
  - — [arXiv:0710.5515](https://arxiv.org/abs/0710.5515)

### Inferences
- **Unit traps (script).**
  - q_D = h/(μ₀e) = 3.291×10⁻⁹ A·m.
  - e/(2α) = 68.52 e in Gaussian units. In Gaussian units g carries charge units, so g_D = 68.5 × 4.8×10⁻¹⁰ esu.
  - μ₀q_D = h/e = 4.14×10⁻¹⁵ Wb = 2Φ₀. Each loop turn sees a 2Φ₀ jump.
- **Cabrera's flux (script).** 1/(t·A·2π sr) = 6.10×10⁻¹⁰ reproduces the PRL value exactly. It is 6.1×10⁵ times the Parker bound, 4.4×10⁶ times MACRO's limit, and 3×10⁹ times IceCube's (relativistic only). The event stays unexplained, but no flux at that level exists [refuted].
- **Parker bounds (script, PDG eqs. 94.9/94.10).** At M = 10¹⁵ GeV the bounds are 10⁻¹⁵ and 10⁻¹⁸. At 10¹⁷ GeV they are 10⁻¹⁵ and 10⁻¹⁶.
- **GUT mass.** M_X/α_GUT with M_X = 10¹⁶ GeV and α = 1/40 gives 4×10¹⁷ GeV, about 7×10⁻¹⁰ kg. Monopoles heavier than about 10¹⁷ GeV are bound by the Galaxy's gravity, not ejected by its field (PDG).
- **Spin ice (script).** q_m = 4.28×10⁻¹³ A·m, and q_D/q_m = 7.7×10³, consistent with "~8000".
- **Calculators:** `dirac_charge(n)`, `flux_from_events(n, t, area)` (acceptance 2πA for a flat loop), `parker_bound(M)`, `parker_seed_bound(M)`.

### Gaps
- The measured spin-ice charge (Bramwell et al. 2009) was not fetched.
- The 2026 self-consistent Parker numbers were not extracted, only the abstract.
- No dedicated ATLAS 2024–25 mass-limit paper was read beyond the PDG entries.

## Q4. Axions and magnetic fields: haloscope power, the CAST limit, and Gertsenshtein conversion in the same magnets

### Takeaway
[mainstream] The haloscope signal is P = g²_aγγ (ρ_a/m_a) B² V C Q (natural HL units). With ADMX's 2018 reference parameters it gives **1.9×10⁻²² W**; the script reproduces 1.94×10⁻²² W. CAST's 9 T, 9.26 m magnet limits g_aγγ < 0.66×10⁻¹⁰ GeV⁻¹ (2017) and < 5.8×10⁻¹¹ GeV⁻¹ (2024). In the same magnet, a graviton converts like an axion with effective coupling √(16πG) = 5.8×10⁻¹⁹ GeV⁻¹. That is 10⁸ times weaker in coupling and 10¹⁶ times weaker in probability: 5.7×10⁻³⁴ against 5.7×10⁻¹⁸.

### Cited Findings
- [mainstream] ADMX, PRL 120, 151301 (2018):
  - P_axion = 1.9×10⁻²² W × (V/136 L)(B/6.8 T)²(C/0.4)(g_γ/0.97)²(ρ_a/0.45 GeV cm⁻³)(f/650 MHz)(Q/50,000).
  - g_γ = −0.97 (KSVZ) and 0.36 (DFSZ).
  - Masses covered: 2.66–2.81 μeV (645–680 MHz). The search excludes DFSZ axions at 0.45 GeV/cm³.
  - — [arXiv:1804.05750](https://arxiv.org/abs/1804.05750)
- [mainstream] The same formula is eq. (5) of the ADMX design paper, which also gives the radiometer equation S/N = (P/k_B T_sys)√(t/b) — [Khatiwada et al., arXiv:2010.00169](https://arxiv.org/abs/2010.00169).
- [mainstream] CAST 2017 used a 9 T LHC test magnet with L = 9.26 m and found g_aγ < 0.66×10⁻¹⁰ GeV⁻¹ at 95% CL — [CAST, Nature Physics 13, 584 / arXiv:1705.02290](https://arxiv.org/abs/1705.02290).
- [mainstream] CAST 2024 (*PRL* 133, 221005; arXiv:2406.16840) set 5.8×10⁻¹¹ GeV⁻¹ at 95% CL for m_a ≲ 0.02 eV. It used an IAXO-pathfinder Micromegas detector with Xe gas and doubled the exposure — [INSPIRE](https://inspirehep.net/literature/2801541).
- [mainstream] The project already uses the Gertsenshtein probability P = 4πG B²L²/(μ₀c⁴) — `../../fringe-physics-and-dirds/notes/baker_li_checks.py` line 26.

### Inferences
- **Unit conversion (script).** 1 T = 195.4 eV² in HL natural units. In Gaussian natural units multiply by √(4π). Both (gBL/2)² and 4πGB²L² assume HL.
- **ADMX reproduction (script).** The formula reproduces 1.94×10⁻²² W, within 2% of the paper. That match also checks the KSVZ relation g = g_γ(α/π)m_a/Λ², with Λ = 78 MeV assumed. At 650 MHz the KSVZ coupling is g = 1.0×10⁻¹⁵ GeV⁻¹, and the signal is 450 photons/s. DFSZ gives 2.7×10⁻²³ W.
- **Conversion probabilities in the CAST magnet (script).** P_aγ = (gBL/2)² is 7.4×10⁻¹⁸ at the 2017 limit and 5.7×10⁻¹⁸ at the 2024 limit. P_gγ = 5.74×10⁻³⁴. The graviton's equivalent coupling √(16π)/M_Pl = 5.81×10⁻¹⁹ GeV⁻¹ reproduces P_gγ to 10⁻³, so the two calculators can share one function.
- **Calculator:** `haloscope_power(V, B, C, g, rho, m, Q)` and `conversion_probability(coupling, B, L)`. Use `coupling = sqrt(16 pi G)` for gravitons. This links to `gravitics::hfgw::pais_conversion` and the HFGW notes, which already recast ADMX-SLIC and ABRA limits as GW strain limits.

### Gaps
- The 2020 and 2021 ADMX PRLs were not fetched, so later reference parameters (for example 7.6 T) are unverified.
- Λ = 78 MeV is inferred from the reproduction, not read from the paper.

## Q5. Magnetism in strong gravity: magnetars, Blandford–Znajek jet power for M87*, and neutron-star mergers

### Takeaway
[mainstream]
- **Blandford–Znajek power.** In SI, P = kΦ²Ω_H²/(μ₀c), with k = 1/(6π) ≈ 0.054 for a monopole field. This equals (k/4π)Φ²Ω²/c in Gaussian units and kΦ²Ω² in Tchekhovskoy's HL units.
- **M87\*.** EHT's 1–30 G horizon-scale fields give at most about 2×10⁴² erg/s. The EHT's jet-power floor is 10⁴² erg/s, and estimates span 10⁴²–10⁴⁵ erg/s. Meeting the floor needs B_H ≳ 23 G, which favours MAD flows, as the EHT found.
- **Magnetars.** Fields reach 10¹¹ T, about 36 times the Schwinger field. Their external field energy is about 2×10⁴⁷ erg.
- **Mergers.** Kelvin–Helmholtz amplification takes 10¹³ G to ≳10¹⁶ G.

### Cited Findings
- [mainstream] Tchekhovskoy, Narayan & McKinney 2010 (*ApJ* 711, 50) work in "Heaviside-Lorentzian units" with G = c = 1:
  - The original BZ77 scaling is P = kΦ²_tot a²/(16M²), with k = (6π)⁻¹ ≈ 0.054 (monopolar) and ≈ 0.044 (paraboloidal).
  - The more accurate BZ2 form is P = kΦ²Ω_H².
  - The BZ6 form is kΦ²(Ω_H² + αΩ_H⁴ + βΩ_H⁶), with α ≈ 1.38 and β ≈ −9.2.
  - — [arXiv:0911.2228](https://arxiv.org/abs/0911.2228)
- [mainstream] EHT M87* Paper VIII (*ApJL* 910, L13, 2021) finds near-horizon B "of order 1–30 G", Ṁ = (3–20)×10⁻⁴ M_⊙ yr⁻¹, and consistent models that are MAD — [arXiv:2105.01173](https://arxiv.org/abs/2105.01173).
- [mainstream] EHT M87* Paper V (*ApJL* 875, L5, 2019):
  - Jet-power estimates "range from 10⁴² to 10⁴⁵ erg s⁻¹".
  - The adopted constraint is P_jet,min = 10⁴² erg s⁻¹. Non-spinning models "do not produce sufficiently powerful jets".
  - Successful jets are powered "through mechanisms akin to the Blandford–Znajek process".
  - The models use M = 6.2×10⁹ M_⊙.
  - — [arXiv:1906.11242](https://arxiv.org/abs/1906.11242)
- [mainstream] Magnetars have fields up to ~10¹⁵ G at the surface.
  - SGR 1806−20's spin-down gives B_p ≈ 1.6×10¹⁵ G and an external field energy of 2×10⁴⁷ erg — [Magnetar giant flares, astro-ph/0604386](https://arxiv.org/pdf/astro-ph/0604386).
  - Its 27 December 2004 giant flare released about 2×10⁴⁶ erg (isotropic) — [astro-ph/0609595](https://arxiv.org/pdf/astro-ph/0609595).
- [mainstream] Kiuchi et al. 2015 used 17.5 m grid spacing. A 10¹³ G seed field was amplified by "at least about 10³", with saturation energy ≳ 4×10⁵⁰ erg — [arXiv:1509.09205](https://arxiv.org/abs/1509.09205).
  - A 2017 run reached about 1% of the thermal energy — [arXiv:1710.01311](https://arxiv.org/abs/1710.01311).

### Inferences
- **4π trap.** The script checks that the SI and Gaussian BZ forms agree to 10⁻⁹. Gaussian Φ (G cm²) = 10⁸ × Φ (Wb).
- **M87\* (script; M = 6.5×10⁹ M_⊙ assumed, Φ ≈ 2πr_H²B_H).**
  - Ω_H = ac/(2r_H).
  - At a = 0.9: 1.9×10³⁹ erg/s at 1 G and 1.7×10⁴² erg/s at 30 G.
  - 10⁴² erg/s needs B_H ≈ 23 G at a = 0.9–0.94, or 32 G at a = 0.5.
  - Ṁc² = 1.7×10⁴³–1.1×10⁴⁴ erg/s, so the observed jets need an efficiency of 1% to above 100%. That is the MAD signature.
- **Magnetar (script).** The external dipole energy is πB_p²R³/(3μ₀) = 2.1×10⁴⁷ erg (Gaussian form B²R³/12). The 2004 flare was 9% of it. B_p/B_Schwinger = 36, which links to the Schwinger notes.
- **Merger (script).** The virial ceiling on the field is about 1.4×10¹⁴ T (1.4×10¹⁸ G). Putting 4×10⁵⁰ erg into a sphere of r = 15 km gives an rms B of 2.7×10¹⁶ G.
- **The field's own gravity (script).** A field's energy curves spacetime with radius R_c = c²√(μ₀/4πG)/B. That is 3.5×10¹⁸ m at 1 T and 2.2×10⁷ m (about 2000 R_NS) for a magnetar. The gravity of magnetic energy is a ~10⁻⁷ correction even there.
- **Calculators:** `bz_power(Phi, a, M, k, order)`, `horizon_omega(a, M)`, `dipole_energy(B_p, R)`, `field_curvature_radius(B)`.

### Gaps
- Sgr A* EHT polarisation numbers (2024) were not fetched.
- The MAD flux normalisation φ ≈ 50 (Gaussian) was not fetched.
- M = 6.5×10⁹ M_⊙ (EHT Paper VI) is from memory; Paper V uses 6.2×10⁹.
- The Kaspi & Beloborodov review was not read.

## Q6. Diamagnetic levitation and Earnshaw: the one real "antigravity"

### Takeaway
[mainstream] Earnshaw's theorem does not forbid stable levitation of induced (diamagnetic) moments. Levitation needs B dB/dz = μ₀ρg/|χ|.
- Water needs 1401 T²/m (Berry–Geim, with χ = −8.8×10⁻⁶), or 1362 T²/m with χ = −9.05×10⁻⁶. That was achieved at B₀ = 16.5 T in a 4 MW Bitter magnet, which levitated a live frog.
- Pyrolytic graphite (χ⊥ ≈ −4.5×10⁻⁴, ρ = 2200 kg/m³) needs only about 60 T²/m. NdFeB arrays reach this, so graphite floats at room temperature.
- This is a magnetic force on magnetisation, not a change in gravity or mass.

### Cited Findings
- [mainstream] Berry & Geim, *Eur. J. Phys.* 18, 307 (1997):
  - "Earnshaw's theorem does not apply to induced magnetism".
  - The levitation condition is B(z)B′(z) = −1400.9 T² m⁻¹ for water, with χ = −8.8×10⁻⁶ (Kaye & Laby) and ρ = 1000 kg/m³.
  - The Bitter magnet has its inflection point at z_i = 78 mm, with B(z_i) = 0.63B₀ and gradient −8.15B₀ T/m. This predicts B₀ = 16.5 T. The stable zone runs 67.5–87.5 mm.
  - The magnet "consumed 4 MW, but ... this is power dissipated in the coils, not power required for levitation". A persistent superconducting magnet needs no power at all.
  - A 10 mm frog carries an equivalent current of about 1.5 A, giving B ≈ 10⁻⁵B₀ ≈ 1.5 G inside it.
  - — [PDF](https://obelix.physik.uni-bielefeld.de/~schnack/Herbstakademie/CD-2008/frog-ejp.pdf)
- [mainstream] Pyrolytic graphite levitates stably at room temperature about 1 mm above crossed NdFeB arrays. Thickness beyond about 0.5 mm adds little lift — [The Wonders of Levitation, arXiv:1010.5761](https://arxiv.org/pdf/1010.5761).
- [mainstream] A 2023 study cites Simon, Heflinger & Geim (*Am. J. Phys.* 69, 702, 2001) for pyrolytic graphite's χ. The rendered digits are corrupted but read as (−85, −85, −450)×10⁻⁶, i.e. χ⊥ = −4.5×10⁻⁴. The study gives ρ = 2.2 g/cm³ — [PMC10304365](https://pmc.ncbi.nlm.nih.gov/articles/PMC10304365/).
- The diamagnetic micro-traps used for gravity-entanglement proposals are already in `../../gravity-anomalies/notes/lab_gravity_manipulation.md` ([arXiv:2411.02325](https://arxiv.org/abs/2411.02325)).

### Inferences
- **Water and the frog (script).**
  - μ₀ρg/|χ| = 1401 T²/m, reproducing eq. (27). With χ = −9.05×10⁻⁶ it is 1362 T²/m.
  - Unit trap: χ_SI = 4π χ_cgs, so water's χ_cgs is −7.2×10⁻⁷.
  - B₀ = √(1400.9/(0.63×8.15)) = 16.52 T. The frog's internal field is 1.45 G.
- **Pyrolytic graphite (script).** It needs 60 T²/m. A 0.5 T surface field that falls to zero over 4 mm gives about 62 T²/m, which matches millimetre levitation heights.
- **Why it mimics weightlessness.** The force per unit mass is (|χ|/μ₀ρ) B∇B. It is the same for every water-like tissue, so levitation imitates free fall (microgravity simulation). The force acts on magnetisation, not on mass, and no gravitational effect is involved. A 16.5 T field stores 1.08×10⁸ J/m³.
- **Calculator:** `levitation_product(rho, chi)`, `required_central_field(profile)`, `is_stable(D1, D2)`.

### Gaps
- Simon–Heflinger–Geim 2001 was not read directly; the graphite χ comes through a corrupted rendering.
- The NdFeB gradient is illustrative, not taken from a measured magnet.

## Q7. Fringe-but-calculable magnetic claims tested against conservation laws and measurements

### Takeaway
[refuted] Every claim with numbers fails conservation by orders of magnitude.
- **MEG.** The patent's "COP 3.44" device has a 39 cm³ NdFeB magnet holding about 10–20 J. At the claimed excess (82 W for two coils) the magnet would empty in about 0.3 s, and the patent itself admits demagnetisation.
- **Roshchin–Godin.** Their "Searl-type" converter claims 7 kW self-running and a 35% weight loss, i.e. 1.2 kN.
  - The magnets store about 10⁴ J and the rotor about 5×10⁴ J, a few seconds of output.
  - A 1.2 kN reactionless lift breaks momentum conservation. An air jet would need about 30 kW, and photons 3.6×10¹¹ W.
  - The claimed 0.05 T "magnetic walls" at 15 m are 4×10⁴ times the magnets' dipole field there. They would store about 2×10⁵ J.
- **Hutchison, Searl and Townsend Brown's magnetic claims.** There are no quantitative primary data to test.

### Cited Findings
- [refuted] **Bearden et al., US 6,362,718 B1**, "Motionless electromagnetic generator" — [Google Patents](https://patents.google.com/patent/US6362718B1/en):
  - Issued 26 March 2002. Inventors Patrick, Bearden, Hayes, Moore and Kenny.
  - The only numbers are extrapolated to 100 V: 140 mA in, 14 W in, and 48 W out per coil at 4000 V. That is "COP of 3.44 for each output coil".
  - Test conditions: 87.5 kHz, 40-turn input and 450-turn output windings, measured data only up to 75 V.
  - The magnet is NdFeB, 40 × 25.4 × 38.1 mm.
  - The device "receives static energy from the magnetic flux of the permanent magnet" and is described as an "open dissipative system", "not ... a perpetual motion machine".
  - "Continued operation of the electromagnetic generator 10 causes demagnetization of the permanent magnet".
  - Status: "Expired – Lifetime" (2020).
- [refuted] Search found no independent peer-reviewed test of the MEG. Hobbyist advice holds that the effect needs a "non linear load" and that a resistor will not show it — [CR4 forum thread](https://cr4.globalspec.com/thread/64399/Information-about-the-MEG) (low-quality source; quoted only as the stated condition).
- [refuted] **Roshchin & Godin**, *Tech. Phys. Lett.* 26(12), 1105–1107 (2000), DOI 10.1134/1.1337268 — [full text (IRI mirror)](https://www.integrityresearchinstitute.org/Propulsion/Roshchin_Godin.PDF); [ADS](https://ui.adsabs.harvard.edu/abs/2000TePhL..26.1105R/abstract):
  - Setup: rare-earth magnets with B_r = 0.85 T, H_c ≈ 600 kA/m and "specific magnetic energy ~150 J/m³" (sic). The stator holds 110 kg of magnets and the rollers 115 kg; the working body is about 1 m across; the platform weighs 350 kg.
  - Running: weight changes began at 200 rpm. At about 550 rpm motor current reached zero and reversed (the "critical regime"). The load was ten 1 kW heaters. The maximum was 7 kW, with ΔG = 35% of 350 kg. The effect reversed sign with rotation direction, with the counter-clockwise onset at about 600 rpm.
  - Side effects: "magnetic walls" of about 0.05 T, 5–8 cm thick, 50–60 cm apart and out to 15 m. The surface ran 6–8 °C below 22 ± 2 °C, as measured with an alcohol thermometer. A corona discharge and ozone were seen.
  - Explanation offered: Dyatlov's "vacuum domain" model. Work was said to be continuing at Glushko NPO Energomash. The setup is credited to J.A. Thomas's 1994 Searl book.
- [refuted] A 2026 preprint uses the Roshchin–Godin report "only as benchmarks and not as evidence". It "reports no anomalous-weight measurement" — [Iadicicco & Verolino, arXiv:2608.21403](https://arxiv.org/abs/2608.21403).
- [refuted] Hutchison claims "the levitation of heavy objects and the fusion of metal and wood". He "is the only one who can produce the effects, but not even he can replicate them—at least not in the presence of unbiased observers" — [Skeptic's Dictionary](https://skepdic.com/hutchisonhoax.html).
- [contested] Ning Li: Li & Torr's GEM equations are standard, and the predicted fields are negligible. Already covered in `../../fringe-physics-and-dirds/notes/superconductor_gravity_claims_math.md` (Q4).

### Inferences
- **MEG (script).**
  - The magnet volume is 3.87×10⁻⁵ m³. Its stored energy is B_r²/2μ₀·V = 22 J at an assumed B_r = 1.2 T, or 12 J at an assumed (BH)max of 300 kJ/m³.
  - The excess claimed for two coils is 82 W, which drains the magnet in 0.27 s. One coil gives about 0.6 s.
  - A static field does no net work around a closed cycle (hysteresis only removes energy). A sustained COP > 1 therefore means a measurement error, typically pulsed, non-sinusoidal power measured into a non-linear load.
- **Roshchin–Godin (script; ρ = 7500 kg/m³ assumed).**
  - **Stored magnet energy.** The magnet volume is 0.030 m³, holding 8.6 kJ (B_r²/2μ₀) or 4.5 kJ if "150" means 150 kJ/m³. At 7 kW that lasts 0.6–1.2 s.
  - **Rotor kinetic energy.** The 115 kg roller ring at r = 0.5 m and 550 rpm holds 47.7 kJ, or 6.8 s at 7 kW. "Energetically autonomous" operation therefore needs a hidden input.
  - **The 1.2 kN lift.** Momentum must leave the device somehow. An ideal air jet over a 1 m disc would need 30 kW (more than the claimed output) and a 25 m/s downdraft. A photon jet would need 3.6×10¹¹ W.
  - **"Magnetic walls".** The magnets' dipole field at 15 m is 1.2 μT, against 0.05 T claimed. Shells of 0.05 T through a 2.5 m-high, 15 m-radius volume would hold 2.1×10⁵ J, which is 25–50 times the magnets' own stored energy.
  - **Cooling as a power source.** Drawing work from room air at a single temperature breaks the Kelvin–Planck second law.
- **Hutchison's "60-lb cannonball" (claim from secondary search snippets; script).**
  - A 27 kg iron sphere needs 9.65 kPa of support under it. Electrostatic pressure at air breakdown (3 MV/m) is 40 Pa, 240 times too little.
  - A magnetic pressure of 0.16 T would do, but a ferromagnet cannot be held stably by static fields (Earnshaw applies; Berry–Geim).
- **Permanent-magnet "over-unity" motors.** The same argument as the MEG: fixed magnets give a conservative force field, so the work around any closed path is zero.
- **Calculator:** an `overunity` module of conservation checks:
  - `stored_magnet_energy(V, Br)`
  - `drain_time(E, P_excess)`
  - `reactionless_power(F, area)`, using actuator-disk and photon limits
  - `dipole_field(m, r)`

### Gaps
- **Searl SEG:** no primary measurement found; the only source in the chain is the 1994 Thomas book.
- **Hutchison:**
  - The 60-lb cannonball figure, Marc Millis's (NASA BPP) "no independent verification" remark and the 2005 admission (no reproduction since about 1991) come only from secondary search snippets I did not verify.
  - So do Hutchison's own timeline claims of a 1983 Army Intelligence/Los Alamos demonstration and a 1986 Canadian intelligence investigation. I found no primary government report.
- **Townsend Brown:** I found no reliable record of later *magnetic* claims; his patents in the record are electrostatic and are covered elsewhere.
- **Official evaluations:** none found for the MEG or Searl. No later Energomash publication on the Roshchin–Godin converter was found.

## Q8. Unified-field attempts joining gravity and electromagnetism, where they yield a checkable number

### Takeaway
- **Kaluza–Klein** [mainstream as a framework] fixes the extra-dimension radius from e: R = √(16πG)ħ/(ec) (Gaussian) = 1.34×10⁻³³ m, which reproduces Klein's 0.8×10⁻³⁰ cm circumference. Its testable remnant is a drift of α; optical clocks bound |α̇/α| to about 10⁻¹⁹/yr.
- **Einstein–Maxwell** gives the extremal charge-to-mass ratio √(4πε₀G) = 8.62×10⁻¹¹ C/kg. The electron exceeds it by 2×10²¹ (weak-gravity conjecture). It is the same constant as in Blackett's SI form.
- **Weyl's** real gauge factor [refuted] became the imaginary Aharonov–Bohm phase eΦ/ħ, which is measured.
- **Weinberg–Witten** is a no-go theorem and gives no number.

### Cited Findings
- [mainstream] In Kaluza–Klein theory, Kerr–KK black holes have g = 2 − v² and β = v/c ≤ 1, and Kerr–Sen holes have β ≤ 2 — [Barrow & Gibbons, arXiv:1701.06343](https://arxiv.org/abs/1701.06343).
- [mainstream] The ¹⁷¹Yb⁺ E3/E2 clock comparison gives (1/α)(dα/dt) = 1.8(2.5)×10⁻¹⁹/yr, and coupling to the solar potential (c²/α)(dα/dΦ) = −2.4(3.0)×10⁻⁹ — [Filzinger dissertation, PTB 2026 (via search summary)](https://oar.ptb.de/resources/show/10.7795/110.20260227).

### Inferences
- **Kaluza–Klein radius (script).** R = √(16πG·4πε₀)ħ/(ce) = √(16π/α)ℓ_P = 1.341×10⁻³³ m, and 2πR = 8.4×10⁻³¹ cm, matching Klein's 1926 value of 0.8×10⁻³⁰ cm (*Nature* 118, 516; not re-fetched). If α ∝ 1/R², the clock bound gives |Ṙ/R| < 3.4×10⁻¹⁹/yr at 2σ. A static KK circle is untestable directly, but any slow drift of the radion is tightly bound.
- **Einstein–Maxwell and the weak-gravity conjecture (script).** e/(m√(4πε₀G)) is 2.04×10²¹ for the electron and 1.11×10¹⁸ for the proton. This is the weak-gravity-conjecture margin (Arkani-Hamed, Motl, Nicolis & Vafa, hep-th/0601001; not re-fetched). It equals Barrow–Gibbons's electron β = N^½, so Q2, Q8 and the GEM "charge" √G·m are one constant.
- **Weyl to Aharonov–Bohm (script).** Weyl's 1918 scale factor exp(∮φ·dx) was dropped after Einstein objected that spectral lines would depend on an atom's history. London's 1927 imaginary version is the AB phase 2πΦ/(h/e), with h/e = 4.14×10⁻¹⁵ Wb and h/2e = 2.07×10⁻¹⁵ Wb. Tonomura et al. 1986 (*PRL* 56, 792; not re-fetched) measured the superconducting half-quantum phase. This is the one Weyl-derived number that has been verified.
- **Weinberg–Witten** (*Phys. Lett. B* 96, 59, 1980; not re-fetched) forbids massless spin > 1 particles with a Lorentz-covariant conserved stress tensor. It rules out an emergent graviton built from ordinary fields, but supplies no number for a calculator.
- **Calculator:** `kk_radius()`, `extremal_charge_to_mass()`, `wgc_ratio(q, m)` and `ab_phase(flux)`. These belong in a `unification` module that shares √(4πε₀G) with `blackett`.

### Gaps
- Klein 1926, Weinberg–Witten, Tonomura 1986 and the WGC paper were cited from known bibliographic data and not re-fetched in this session.
- The α̇ bound comes from a PTB dissertation summary, not the PRL.
- The Einstein–Weyl objection was not given a quantitative bound from a primary source.
