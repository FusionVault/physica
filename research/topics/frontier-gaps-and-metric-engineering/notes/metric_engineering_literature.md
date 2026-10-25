# Spacetime metric engineering, 1994–2026: warp drives, traversable wormholes, quantum inequalities and claimed experiments

Scope: closed forms and published numbers that `gravitics::spacetime` can compute and test. This note builds on `research/topics/fringe-physics-and-dirds/notes/dirds_gravity_and_metric_engineering.md` (the "DIRD review"). It does not repeat what that review already established: the PF Eq. 28 discrepancies, the Puthoff PV checks, the Fuchs shell parameters, the Visser cube and the VKD integrals. It adds the papers and numbers that review listed as gaps.

**Method.** Every source was opened from its arXiv TeX (`me_sources/<id>/`), from NASA NTRS PDFs, or through the Springer HTML for White et al. 2021. Downloads were treated as data, and none of them was executed. Every published number was recomputed by scripts in `me_scripts/`:
- `metric_engineering_checks.py` (+`.out`): published numbers versus computed values.
- `vdb_alpha_check.py` (+`.out`): the Van Den Broeck parameter puzzle.
- `metric_to_energy_conditions.py` (+`.out`): a general metric → Einstein tensor → energy-condition evaluator, validated on Schwarzschild, Alcubierre, Natário, Morris–Thorne and the Santiago–Schuster–Visser closed forms.
- `hawking_ellis_check.py` (+`.out`): exact energy conditions from eigenvalues, compared with observer sampling.

**Labels.**
- **[GR]**: standard GR/QFT, accepted.
- **[SPEC]**: valid mathematics, physically speculative or disputed.
- **[FRINGE]**: unsupported or contradicted.
- **[checked]**: I reproduced the number.
- **[mismatch]**: I did not reproduce it; the ratio is given.

Units: v is in units of c; c²/G = 1.3466×10²⁷ kg/m; c⁴/G = 1.2104×10⁴⁴ J/m; l_P = 1.616×10⁻³⁵ m.

---

## Q1. Warp-drive energy: Alcubierre, Pfenning–Ford, Van Den Broeck, Natário, Krasnikov — what do they give and do the published numbers check out?

### Takeaway
Every superluminal metric in this family needs negative Eulerian energy. Four things are new here.
- **Alcubierre's wall.** For Alcubierre's own tanh wall the exact total is **E = −(v²/36)R²σ·c²/G**. That is 2/3 of Pfenning–Ford Eq. 28, which assumes a piecewise-linear wall, and 1/36 of the Lobo–Visser estimate "−v²R²σ".
- **A floor on the energy.** For a flat interior of radius r₀ there is a hard lower bound, **|E| ≥ (v²/12) r₀ c²/G**. This is 5.6×10⁻³ M_⊙ for r₀ = 100 m at v = 1.
- **Natário's drive is far worse.** It scales as **R⁴σ³/180**: about 3,000 M_⊙ for R = 100 m and a 1 m wall, roughly 8,000 times Alcubierre.
- **Van Den Broeck's quantum-inequality check fails with the parameters it lists.** It used α = 10³⁴ and Δ̃ = 10⁻³² m. With the printed α = 10¹⁷ and Δ̃ = 10⁻¹⁵ m, the QI is violated by a factor of 2×10³³.

### Cited Findings
**Alcubierre (1994), CQG 11, L73** [GR] — [arXiv:gr-qc/0009013](https://arxiv.org/abs/gr-qc/0009013)
- Metric:
  - ds² = −dt² + (dx − v_s f(r_s)dt)² + dy² + dz².
  - f = [tanh(σ(r_s+R)) − tanh(σ(r_s−R))]/[2 tanh(σR)].
- Expansion θ = v_s (x_s/r_s) df/dr_s.
- Eulerian energy density, as printed: "T^{αβ}n_αn_β = α²T⁰⁰ = (1/8π)G⁰⁰ = −(1/8π) v_s²ρ²/(4r²)(df/dr_s)²".
  - The paper: "everywhere negative implies that the weak and dominant energy conditions are violated… similarly… the strong energy condition is also violated".
- **SI form:** ρ_E = −(c⁴/32πG) v² f′² (y²+z²)/r_s².
- **[checked]** A finite-difference Einstein tensor reproduces this to 10⁻⁶ at five points (`metric_to_energy_conditions.out`). At (0, R, 0) with σR = 8 it is exactly −1/(2π) in G = c = 1, since f′(R) = −σ/2.

**Total Eulerian energy (exact, any shape function)** [GR]
- E = −(v²/12)(c²/G)∫₀^∞ r² f′² dr.
- The same integral appears as Pfenning–Ford's unnumbered equation before Eq. 28 ([gr-qc/9702026](https://arxiv.org/abs/gr-qc/9702026)), as Lobo–Visser's "M_warp" ([gr-qc/0406083](https://arxiv.org/abs/gr-qc/0406083) §2), and in Bobrick–Martire's appendix ([2102.06824](https://arxiv.org/abs/2102.06824)).

**Three versions of the wall [checked numerically, σR = 2 to 20,000]** (`metric_engineering_checks.out` §1):
- **PF Eq. 28 (piecewise-linear wall Δ):** E = −(v²/12)(R²/Δ + Δ/12) c²/G. This is the library's `alcubierre_energy`.
- **Alcubierre's tanh wall, σR ≫ 1:** ∫r²f′²dr → R²σ/3 exactly, so **E = −(v²/36)R²σ c²/G**. PF match the two walls by slope, Δ = [1+tanh²σR]²/(2σ tanh σR) ≈ 2/σ (PF §2). With that matching the tanh energy is **0.6667 × PF Eq. 28** at every σ tested.
- **Lobo–Visser:** "we can estimate M_warp ≈ −v²R²σ" (LV §2.1). That estimate is **36× the exact value**. DIRD 13's "E_warp = −v²c⁴R²σ/G" inherits the 36× (see the DIRD review §Q2).

**Pfenning & Ford (1997)** [GR] — [gr-qc/9702026](https://arxiv.org/abs/gr-qc/9702026)
- The quantum-inequality wall bound is "Δ ≤ (3/4)√(3/π) v_b/α²". With α = 1/10 the paper writes "Δ ≤ 10² v_b L_Planck".
  - **[checked]** The coefficient is 73.3, so "10²" is rounded up.
- The other text discrepancies (6.2×10⁶⁵ g against 6.9×10⁶⁵ g; "6.2×10⁷⁰ L_P" against 3.2×10⁷⁰; "a quarter" against 0.56 M_⊙) are already in the library docs and the DIRD review. They are reconfirmed in `.out` §1.

**Minimum energy for a given flat interior** [GR; closed form from Bobrick–Martire's variational result]
- BM, Appendix "Optimizing the energy requirements": the Euler–Lagrange equation for ℒ = r²f′² gives "f̄(r_s) = min(r₀/r, 1)". Using it "reduces the energy requirement… by about a factor of three".
- **[checked]** ∫_{r₀}^∞ r²(r₀/r²)²dr = r₀, so **E_min = −(v²/12) r₀ c²/G**. This is the electrostatic capacitance of a sphere, so it is a true lower bound for any f with f = 1 on r ≤ r₀.
  - For r₀ = 100 m and v = 1: 1.122×10²⁸ kg = 5.6×10⁻³ M_⊙ = 5.9 M_Jupiter.
- **[checked]** BM's "factor of three" is the ratio to Alcubierre's own example, σR = 8: I get 2.68.
- A tanh wall with no flat interior (σR ≈ 1.23) reaches 0.667 R. That beats BM's "optimum" only because it gives up the flat region.

**Van Den Broeck (1999), CQG 16, 3973** [GR/SPEC] — [gr-qc/9905084](https://arxiv.org/abs/gr-qc/9905084) (arXiv v3, Sept 1999; the CQG text was not checked)
- Metric: ds² = −dt² + B²(r_s)[(dx − v_s f dt)² + dy² + dz²].
  - B = α(−(n−1)wⁿ + nw^{n−1}) + 1, with w = (R̃+Δ̃−r)/Δ̃.
  - The Eulerian density is "T⁰⁰ = (1/8π)(B′²/B⁴ − 2B″/B³ − 4B′/(B³r))".
- Listed numbers (§3): "α = 10¹⁷, Δ̃ = 10⁻¹⁵ m, R̃ = 10⁻¹⁵ m, R = 3×10⁻¹⁵ m", which gives an inner "pocket" 200 m across.
- **Energies reproduced [checked]:**
  - E_IV = −(1/12)v_s²[(R+Δ/2)²/Δ + Δ/12] with Δ = 100 l_P gives 6.25×10²⁹ kg (paper "6.3×10²⁹ v_s kg").
  - E_II,− = −1.38×10³⁰ kg (paper −1.4×10³⁰).
  - E_II,+ = 4.87×10³⁰ kg (paper 4.9×10³⁰).
  - T⁰⁰ changes sign at w = 0.9812 (paper 0.981).
  - E_II scales as αΔ̃, as the paper says: "energies do not vary much with α if R̃ = Δ̃ and αR̃ = 100 m".
- **[mismatch] The quantum-inequality check.** The paper gives:
  - a T⁰⁰ peak of "−4.9×10²/Δ̃²" at "w = 0.349";
  - "r_c,min = Δ̃/72.5 = 1.4×10⁻³⁴ m";
  - "−6.6×10⁹³ kg/m³" against the QI bound "−9.2×10⁹⁴ kg/m³".

  What I find (`vdb_alpha_check.out`):
  - With the printed α = 10¹⁷ the peak is −176/Δ̃² at w = 0.576, and r_c = Δ̃/43.4 = 2.3×10⁻¹⁷ m.
  - With **α = 10³⁴ and Δ̃ = 10⁻³² m** (still αR̃ = 100 m) I get exactly −492/Δ̃² at w = 0.3490 and r_c = Δ̃/72.46. With τ₀ = 0.1 r_c the comparison is then 6.6×10⁹³ against 9.2×10⁹⁴ kg/m³, matching the paper.
  - With the printed parameters, |T⁰⁰| = 2.4×10⁵⁹ kg/m³ against a QI limit of 1.2×10²⁶. That is **violated by 2×10³³**.

  So "the QI is amply satisfied" holds only for a Planck-scale wall, Δ̃ ~ 600 l_P. The energies are unaffected.
- Status: the energy reduction is real arithmetic [GR]. The Alcubierre-sector wall is still Planck-thin. Warp Factory notes that VdB "does not address the key issue which creates the energy condition violations" ([2404.03095](https://arxiv.org/abs/2404.03095) §4.2).

**Natário (2002), "Warp drive with zero expansion", CQG 19, 1157** [GR] — [gr-qc/0110086](https://arxiv.org/abs/gr-qc/0110086)
- Shift field (polar axis x): **X = −2v_s f cosθ e_r + v_s(2f + r f′) sinθ e_θ**, with "f(r) = ½ for large r and f(r) = 0 for small r". It is divergence-free, so θ = 0.
- Energy density: "ρ = −(1/16π)K_ijK^ij = −(v_s²/8π)[3(f′)²cos²θ + (f′ + (r/2)f″)² sin²θ]".
- General theorem (Natário §1): ρ = (1/16π)(θ² − K_ijK^ij), so a zero-expansion drive has ρ ≤ 0 everywhere.
- **[checked]** The finite-difference Einstein tensor matches the formula to 6 digits at four points.
- **Total energy, my integration [no published number found]:**
  - E = −(c²/G)[(v²/2)∫f′²r²dr + (v²/3)∫(f′ + rf″/2)²r²dr].
  - Equivalently −(v²/12)(c²/G)∫(4f′ + rf″)²r²dr, which follows from SSV's vorticity identity (below). The two agree to 10⁻⁹.
  - For a thin tanh wall, with f = (1 − f_Alc)/2: **E ≈ −(v²/180)R⁴σ³ c²/G**.
  - R = 100 m, σ = 2 m⁻¹, v = 1 gives **−3.0×10³ M_⊙**, against −0.38 M_⊙ for Alcubierre with the same R and σ (ratio 8,000).
  - Reason: the incompressible flow must squeeze the whole outside stream through the wall, so the tangential shift reaches ~vRσ/4, i.e. 50v for R = 100 m and σ = 2 m⁻¹.

**Krasnikov (1998) and Everett & Roman (1997)**
- **Krasnikov** [SPEC] — [gr-qc/9511068](https://arxiv.org/abs/gr-qc/9511068).
  - 2D metric ds² = −(dt − dx)(dt + k dx), with k = 1 − (2−δ)θ_ε(t−x)[θ_ε(x) − θ_ε(x+ε−D)].
  - The paper concedes that creating it "requires violation of the weak energy condition". It cites Everett–Roman for "huge amounts (e.g. 10³² M_galaxy) of 'negative energy'".
- **Everett & Roman** [GR] — [gr-qc/9702049](https://arxiv.org/abs/gr-qc/9702049)
  - 4D tube: ds² = −dt² + (1−k)dx dt + k dx² + dρ² + ρ²dφ².
  - Energy density seen by a static observer (ER Eq. T): T_tt = [32π(1+k)²]⁻¹[−4(1+k)k_ρ/ρ + 3k_ρ² − 4(1+k)k_ρρ].
  - With θ_ε(ξ) = ½{tanh[2(2ξ/ε − 1)] + 1}, their fitted end-cap value is "T = (1/8π)·0.271(3.504 − 4.034n …)/(ε²(n−1)(1.018+0.491δ)²) ≈ −1/(8πε²)".
  - **[checked]** Recomputing Eq. T directly agrees with their fit to 0.1% at four (n, δ) pairs. T_tt·8πε² = −1.056 at large n.
  - QI wall bound "ε ≲ l_P/σ² … 10⁴ l_P ≈ 10⁻³¹ m". **[checked]** 1.6×10⁻³¹ m.
  - Energy estimate "E ≈ −αρ_max D/ε" gives "−α10⁶³ g" for a 1 m × 1 m tube with ε = 100 l_P, and "−10³² M_galaxy" to the nearest star.
    - **[checked to order of magnitude]** T·V = (1/8πε²)(2παερD) is actually αρD/(4ε), so ER drop a factor of 4. They also round 1 m to 10³⁵ l_P; it is 6.2×10³⁴.
    - My values are 9.6×10⁶⁶ α m_P, and 1.7×10³² M_galaxy for the stellar tube (without the 1/4).
  - Loophole they compute: a nearly-flat tube with η = 2−δ ≈ l_P²ε²/τ₀⁴ ≈ 10⁻⁶⁶ satisfies the QI. Light inside then beats c by only 1 part in 10⁶⁶, "completely unobservable". **[checked]** 2.6×10⁻⁶⁶.

### Inferences
- The library's `alcubierre_energy` (PF Eq. 28) is correct for PF's linear wall. For Alcubierre's tanh wall it overstates the energy by 1.5×, given PF's own Δ–σ matching. The docs should say so, or a σ-based variant should be added.
- The capacitance bound E_min = −(v²/12) r₀ c²/G is the cleanest number to compute. Any shape function with a flat interior of radius r₀ needs at least this much Eulerian negative energy. Thicker walls cannot go below it, which is why White's "thick wall" savings (Q4) have a floor.
- Natário-type drives, sometimes presented as "gentler" because nothing is compressed, need orders of magnitude more negative energy for thin walls.

### Gaps
- I found no published total-energy number for Natário's drive. The R⁴σ³/180 law is mine, and it is cross-checked two ways.
- Van Den Broeck's published CQG text was not opened. The α = 10³⁴ hypothesis explains every QI-section number exactly, but I could not confirm which α the journal version prints.
- The journal references for Krasnikov (PRD 57, 4760) and Everett–Roman (PRD 56, 2100) were not confirmed on a fetched page.

---

## Q2. Positive-energy and "physical" warp drives (Lentz; Bobrick–Martire; Fell–Heisenberg; Santiago–Schuster–Visser; Fuchs et al.) — claims, formulas, consensus

### Takeaway
The consensus is settled. In the Natário class (unit lapse, flat slices, any shift), the Eulerian energy density is a 3-divergence minus the vorticity squared. So ∫ρ d³x = −(1/32π)∫|∇×v|² ≤ 0, and every superluminal drive of this kind violates the NEC.
- **Lentz's "positive energy":** contradicted by that identity, and refuted by direct computation in a 2025 preprint (sign error).
- **Fell–Heisenberg:** their own paper concedes WEC, SEC and DEC violations for non-Eulerian observers.
- **What survives** is subluminal: Bobrick–Martire's positive-energy shells, and Fuchs et al.'s 2.4 M_Jupiter shell at 0.02c with no violations found.
- **Newest work.** Garattini–Zatrimaylov (2025) evades the identity with a shift that does not vanish at infinity. A 2026 toolkit (warpax) shows that Eulerian or sampled checks underestimate violations.

### Cited Findings
**Santiago, Schuster & Visser (2022), PRD 105, 064038** [GR] — [arXiv:2105.03079](https://arxiv.org/abs/2105.03079). These closed forms are what a library should implement for shift-only metrics.
- Conventions: metric −dt² + δ_ij(dx^i − v^i dt)(dx^j − v^j dt), Eulerian n^a = (1, v^i), K_ij = v_(i,j).
- "ρ = (1/16π)(K² − tr(K²))". Equivalently ρ = (1/16π){∂_i(v_i v_{j,j} − v_j v_{i,j}) − ½ω·ω}, with ω = ∇×v (Eq. E:divergence).
- Flux: "f_i = (1/16π)(∇×(∇×v))_i". It is zero for gradient flows.
- Stresses: "T_ij = (1/8π)(ℒ_nK_ij + KK_ij − 2(K²)_ij − (ℒ_nK + ½K² + ½tr K²)δ_ij)".
- NEC-relevant combination: "ρ + p̄ = (1/24π)(−2ℒ_nK + K² − 3 tr(K²))".
- With fall-off, "∫ρ d³x = −(1/32π)∫(ω·ω)d³x ≤ 0". For zero-vorticity (Lentz/Fell–Heisenberg) drives "∫ρ d³x = 0", so "if the warp field has positive energy density anywhere, it must have negative energy density somewhere else".
  - The fall-off needed is Φ_i = o(r^{−1/2}).
- NEC theorem: the NEC implies dK/dτ ≤ −(3/2)tr([K^tf]²) ≤ 0 along Eulerian worldlines. A localized bubble forces K to return to zero after passing, so "any warp drive in the generic Natário form will inevitably violate the NEC". "Because the NEC is violated then so are the WEC, SEC, and DEC."
- SSV note that their theorem does not cover non-unit-lapse or Van Den Broeck metrics, which "cannot simply be dismissed out of hand".
- **[checked]**
  - ρ + p̄ from the finite-difference Einstein tensor matches SSV's formula to 4–5 digits at six Alcubierre points, with K = +v_(i,j) and ℒ_n = ∂_t + v^i∂_i. The opposite sign of K fails.
  - The vorticity identity reproduces Alcubierre's total exactly: ω² = v²f′²(y²+z²)/r², giving −(v²/12)∫r²f′².
  - It also reproduces Natário's total to 10⁻⁹ (Q1).
  - On the axis of the front wall Eulerian ρ ≈ 0, yet the NEC minimum is −0.68 (G = c = 1, v = 1, σR = 8). This is SSV's point: Eulerian positivity says nothing about the NEC.

**Lentz (2021), "Breaking the warp barrier: hyper-fast solitons in Einstein–Maxwell-plasma theory", CQG 38, 075015, doi:10.1088/1361-6382/abe692** [FRINGE as to positivity] — [arXiv:2006.07125](https://arxiv.org/abs/2006.07125). The journal reference is confirmed in the Celmaster–Rubin bibliography.
- Claim: a shift potential N_i = ∂_iφ obeying a wave equation, ∂_x²φ + ∂_y²φ − (2/v_h²)∂_z²φ = ρ. Under an l₁-norm ansatz this gives E = (1/16π)(2∂_z²φ(ρ + (2/v_h²)∂_z²φ) − 4(∂_z∂_xφ)²) (Eq. Edens). The abstract claims solitons "sourced by purely positive energy densities".
- Energy: "E_tot ~ C v_s² R²/w", C ~ O(1). R = 100 m and w = 1 m give "E_tot ~ (few)×10⁻¹ M_⊙ v_s".
  - **[checked/inconsistent]** R²/w·c²/G = 6.77 M_⊙, so 0.3 M_⊙ implies C ≈ 0.044.
  - The text's result is linear in v_s while its formula is quadratic.
- **Refutation** [preprint]: Celmaster & Rubin (Nov 2025), "Violations of the Weak Energy Condition for Lentz Warp Drives" — [arXiv:2511.18251](https://arxiv.org/abs/2511.18251).
  - "We demonstrate that Lentz's claim is incorrect… there are spacetime regions where the energy density is negative."
  - They trace "the illusion of a positive energy density" to "an overall sign error in the expression for φ_L". They add that Lentz's Fig. 3 positivity is an artefact of its colour display (γ = 1). They also say a modified, corrected geometry "still violates the WEC even in the Eulerian reference frame".
  - I found no journal version.
- Warp Factory's "Lentz-inspired" metric shows energy-condition "violations… consistently at the interfaces of the rhomboid structure". The cause is "non-zero pressure terms at all boundaries" ([2404.03095](https://arxiv.org/abs/2404.03095) §4.4).

**Fell & Heisenberg (2021), "Positive energy warp drive from hidden geometric structures", CQG 2021 (doi 10.1088/1361-6382/ac0e47; volume and article number not confirmed)** [SPEC; Eulerian-only] — [arXiv:2104.06488](https://arxiv.org/abs/2104.06488)
- Numbers (§Numerics):
  - "ρ_max ≈ 3.2×10²⁶ kg/m³";
  - "E_total ≈ 9.25×10⁴³ J";
  - "E_⊙ ≈ 1.78×10⁴⁷ J";
  - hence "four orders of magnitude smaller than the rest mass energy of the sun".
- **[checked]**
  - M_⊙c² = 1.787×10⁴⁷ J.
  - The ratio is 5.2×10⁻⁴, so 3.3 orders rather than four.
  - E_total = 1.03×10²⁷ kg ≈ 0.54 M_Jupiter.
  - ρ_max is 1.4×10⁹ times nuclear density. FH themselves: "any attempt… will form a black hole".
- FH's own concessions: "the WEC is violated in compact regions within the energy distribution" for non-Eulerian frames. They also write that the configurations "violate the SEC" and that "the full dominant energy condition is also violated".
- SSV on FH: "the first explicit example… the velocity potential Φ is certainly C⁰ but is only piecewise differentiable… the metric is discontinuous… not viable".

**Bobrick & Martire (2021), "Introducing physical warp drives", CQG 38, 105009** [GR] — [arXiv:2102.06824](https://arxiv.org/abs/2102.06824)
- Main claims: "any warp drive… is a shell of regular or exotic material moving inertially", and "any warp drive requires propulsion".
- Flattening the drive by α_X gives "E → E/α_X"; choosing α_X = 1+v² removes the v² dependence (§3.3).
- Optimal Alcubierre shape f = min(r₀/r, 1), "about a factor of 3" (Q1). The abstract's "two orders of magnitude" combines flattening with shape.
- BM say extreme flattening allows superluminal solutions that satisfy the quantum inequalities. They then add that such solutions "more likely probe the limits of applicability of the quantum inequalities", since superluminal motion violates the ANEC.

**Fuchs, Helmerich, Bobrick, Sellers, Melcher & Martire (2024), "Constant velocity physical warp drive solution", CQG 41, 095013** [GR, numerical] — [arXiv:2405.02709](https://arxiv.org/abs/2405.02709). The parameters are in the DIRD review: R₁ = 10 m, R₂ = 20 m, M = 4.49×10²⁷ kg = 2.365 M_J, β_warp = 0.02, no violations.
- **[checked]** M/M_J = 2.365 and 2GM/(c²R₂) = 0.333.
- The method is a TOV shell with smoothed density and pressure, plus a compact shift inside. It is subluminal and constant-velocity only; acceleration is not modelled.

**Warp Factory (Helmerich et al. 2024), CQG, doi 10.1088/1361-6382/ad2e42** — [arXiv:2404.03095](https://arxiv.org/abs/2404.03095); [repo](https://github.com/NerdsWithAttitudes/WarpFactory) (MIT). Source read as text only.
- Alcubierre test case: "v_s = 0.1 c, R = 300 m, σ = 0.015 m⁻¹", 1 m grid, 4th-order finite differences.
- The code builds the comoving shift as `beta{1} = v*(1-fs)` with Alcubierre's shape function, and converts with `c^4/(8*pi*G)` (files `metricGet_AlcubierreComoving.m`, `shapeFunction_Alcubierre.m`, `einE.m`).
- **[mismatch, unresolved]**
  - Analytically, and by my finite-difference evaluation of that exact comoving metric, the peak is ρ_E = **−6.775×10³⁵ J/m³** on the ring at r = 300 m.
  - The paper's Fig. "Energy density for the Alcubierre metric" colour bar bottoms out near **−2.7×10³⁶ J/m³**, 4.0× larger.
  - The ring position and width match σ = 0.015. A factor of 4 = 2² suggests the figure was made with v = 0.2 or σ = 0.03, not the captioned values.
  - Total Eulerian energy for the captioned case: −5.1×10²⁶ kg.

**2025–2026 developments**
- **Garattini & Zatrimaylov (2025), "Warp drive in a de Sitter universe"** [SPEC; preprint] — [arXiv:2502.13153](https://arxiv.org/abs/2502.13153)
  - "if one takes a shift vector that is both irrotational… and asymptotically non-vanishing (or falling off slower than r^{−1/2}), it's possible to produce ρ that is non-negative everywhere".
  - With the bubble moving radially at the Hubble flow, the WEC and NEC are "satisfied up to a total divergence term that averages to zero". So they are locally violated.
- **Le (2026), "Observer-robust energy condition verification for warp drive spacetimes"** — [arXiv:2602.18023](https://arxiv.org/abs/2602.18023), via the arXiv abstract page; v6 dated Sept 2026, submitted to CQG.
  - The warpax toolkit (Python, JAX autodiff) uses Hawking–Ellis classification plus optimisation over observers. It reports that more than 96% of grid points are type I, where an algebraic eigenvalue check is exact.
  - "For the Rodal metric… Eulerian analysis misses DEC violations at over 28% of grid points and WEC violations at over 15%."
  - Alcubierre WEC violations come out "~90,000 times larger" at rapidity cap ζ = 5, scaling as e^{2ζ}.
  - **[consistent]** With the NEC violated, my sampled WEC minimum grows as γ²: −39.4 at β = 0.99, against an NEC minimum of −0.796 (×γ² = 50).
- **Barzegar (2026), SISSA talk** — [indico](https://indico.sissa.it/event/175/sessions/722/). The abstract claims no-go results for "Restricted Warp" models, including NEC violation and no globally hyperbolic superluminal model. Only the search snippet was seen.

### Inferences
- Current consensus [GR]:
  - Superluminal warp drives of the Natário class violate the NEC, and with it every pointwise condition.
  - "Positive energy" claims rest on the Eulerian observer alone (Lentz, FH), on non-decaying shifts (Garattini–Zatrimaylov), or on a sign error (Lentz, per Celmaster–Rubin).
  - Subluminal, positive-ADM-mass shells (BM, Fuchs) are legitimate GR solutions. They are better described as massive moving shells with a flat interior than as "drives": no propulsion mechanism exists, and the shell is 2.4 M_J.
- For the library, energy-condition output should come from the Hawking–Ellis eigenvalues (exact at type I points), with the NEC on normalized null vectors as the primary number. "WEC minima over observers" are unbounded wherever the NEC fails, so they should not be reported as a magnitude.

### Gaps
- The Celmaster–Rubin preprint has no journal status yet. I did not reproduce their Lentz computation; their source PNGs were not used.
- I did not obtain warpax's licence or its result for the Fuchs-type shell. The "Rodal metric" paper itself was not opened.
- The 4× Warp Factory figure discrepancy was not resolved: I did not run the MATLAB code, deliberately, as downloaded code is not executed.

---

## Q3. Traversable wormholes: Morris–Thorne, thin shells, Maldacena et al., and quantum-inequality bounds (with checks of the library's tension functions)

### Takeaway
- **The library's tension functions are correct.** `wormhole_throat_tension` = c⁴/(8πG b₀²) is the general Morris–Thorne radial tension at any throat, whatever the redshift function. `wormhole_edge_tension` = −c⁴/(8G) matches Visser's cube. A 1 m thin-shell wormhole needs −2.7×10²⁷ kg, the same order as the 1 m cube's 2.0×10²⁷ kg.
- **Quantum inequalities confine exotic matter.** They limit a QFT-supported throat to about 10⁻³¹ m. Alternatively, the exotic matter must sit in a band of about 10⁻²¹ m for a 1 m throat.
- **Maldacena–Milekhin's "humanly traversable" wormhole** reproduces to ≤6%. Its mouths are r_e ≈ 1.5×10⁷ m and have the mass of about 10⁴ Suns.
- **The Fewster–Eveson bound is tighter.** For Lorentzian sampling it is 9/64 of the Ford–Roman bound now in the library.

### Cited Findings
**Morris–Thorne energy conditions** [GR] — Ford & Roman (1996), PRD 53, 5496 — [gr-qc/9510071](https://arxiv.org/abs/gr-qc/9510071)
- For Φ = 0 and b = r₀²/r: "ρ = p_r = −P = −r₀²/(8πr⁴)", with curvature radius r_c = r₀ at the throat.
- **[checked]** The finite-difference evaluator in Cartesian coordinates gives ρ, p_r and ρ+p_r = −2r₀²/(8πr⁴) to 6 digits at r = 1.2 and 2.0. The NEC minimum equals the radial ρ + p_r.
- General throat: p_r = (1/8π)[−b/r³ + 2(1−b/r)Φ′/r] → −1/(8πr₀²) at b = r₀ for any Φ. So the library's tension τ₀ = c⁴/(8πG b₀²) holds for every Morris–Thorne throat. **[checked]** 4.815×10⁴⁰ Pa at 10 m.
- NEC at the throat: Fewster & Roman (2005) write the null-contracted density as "(b₀′ − 1)/(8πr₀²l_p²) < 0" — [gr-qc/0507013](https://arxiv.org/abs/gr-qc/0507013).

**Thin-shell wormholes: Visser (1989), Nucl. Phys. B 328, 203** [GR] — [arXiv:0809.0927](https://arxiv.org/abs/0809.0927); **Poisson & Visser (1995)** — [gr-qc/9506083](https://arxiv.org/abs/gr-qc/9506083)
- Static shell at radius a in Schwarzschild mass M:
  - "σ = −(1/2πa)√(1−2M/a)";
  - surface tension "ϑ = −(1/4πa)(1−M/a)/√(1−2M/a)", which is a pressure.
- Dynamic case: σ = −(1/2πa)√(1−2M/a+ȧ²).
- Poisson–Visser stability: stable iff 2M/a₀ + (M²/a₀²)/(1−2M/a₀) + (1+2β₀²)(1−3M/a₀) < 0. So "large traversable wormhole (a₀ > 3M) are stable only for β₀² < 0".
- **SI [checked]:** σ(a = 1 m, M = 0) = −(c⁴/G)/(2π·1 m) = −1.926×10⁴³ J/m². The total 4πa²σ/c² = −2a·c²/G·√(1−2GM/ac²) = −2.69×10²⁷ kg.
  - For comparison, Visser's 1 m cube has 12 edges at −1.513×10⁴³ J/m, i.e. −2.02×10²⁷ kg (DIRD review). The two constructions agree in order of magnitude.

**Quantum inequalities: Ford–Roman, Fewster–Eveson, Fewster–Roman** [GR/QFT]
- FR Lorentzian QI (FR96 Eq. QI): "(τ₀/π)∫⟨T_μνu^μu^ν⟩dτ/(τ²+τ₀²) ≥ −3/(32π²τ₀⁴)". This is the library's `quantum_inequality`.
- **Fewster & Eveson (1998), PRD 58, 084010** — [gr-qc/9805024](https://arxiv.org/abs/gr-qc/9805024)
  - General massless 4D bound: "ρ ≥ −(1/16π²)∫dt (f^{1/2}″(t))²", for any smooth sampling function f.
  - "For the Lorentzian function… = −27/(2048π²τ⁴), which is 9/64 of Ford and Roman's result".
  - **[checked]** Quadrature gives 27/(2048π²) exactly. A normalized Gaussian weight e^{−t²/τ²}/(√π τ) gives −3/(64π²τ⁴), i.e. half the FR constant.
  - In SI the general bound is ⟨ρ⟩_f ≥ −(ħ/(16π²c³))∫(f^{1/2}″)²dt.
  - Kontou & Sanders (2020) restate this as the standard minimally-coupled scalar QEI — [arXiv:2003.01815](https://arxiv.org/abs/2003.01815) §3.3.
- **Throat-size bound** (FR96 §IV, b = r₀²/r, τ₀ = f r₀): "r₀ ≲ l_p/(2f²)… for f ≈ 0.01, r₀ ≲ 10⁴ l_p = 10⁻³¹ m". **[checked]** The exact form is r₀ ≤ √(3/4π) l_P/f², giving 7.9×10⁻³² m.
- **"Absurdly benign" band** (FR96, Eq. AB;Ba0): "a₀ ≲ (r₀/(8f⁴l_p))^{1/3} l_p". With r₀ = 1 m this is "10¹⁴ l_p ≈ 10⁻²¹ m"; with 1 ly "≈ 0.2 fermi"; with 10⁵ ly "≈ 10⁻¹⁴ m".
  - **[checked]** The exact form, (3r₀l_P²/(8πf⁴))^{1/3}, gives 1.46×10⁻²¹ m, 3.1×10⁻¹⁶ m and 1.43×10⁻¹⁴ m. FR's "0.2 fm" is 0.31 fm.
- **Many fields:** "to get r₀ ≈ 1 m, we would need… 10⁶² fields". **[checked]** 1.5×10⁶².
- **Fewster & Roman (2005)**, null-projected QI on VKD-type wormholes ([gr-qc/0507013](https://arxiv.org/abs/gr-qc/0507013)): "r₀ ≲ 10⁵ l_p/√(1−b₀′); to approach 10²⁰ Planck lengths (≈ 1 fermi) one needs 1−b₀′ ≤ 10⁻³⁰". **[checked]** 1.6×10⁻¹⁵ m.

**Maldacena, Milekhin & Popov (2018), "Traversable wormholes in four dimensions"** [GR+QFT, SPEC as to realisation] — [arXiv:1807.04726](https://arxiv.org/abs/1807.04726)
- "a solution of an Einstein Maxwell theory plus charged massless fermions. The fermions give rise to a negative Casimir-like energy". It is "a long wormhole that does not lead to causality violations".
- Formulas: "ℓ = 16 r_e³/(G_N q) = 16π^{3/2}q²l_p/g³" and "E_min = −G_N q²/(256 r_e³)". It embeds in the Standard Model only "by making its overall size small compared to the electroweak scale".
- Traversal time is of order the light-crossing time r_e, while ℓ ≫ r_e outside.

**Maldacena & Milekhin (2021), "Humanly traversable wormholes", PRD 103, 066007** [SPEC; Randall–Sundrum assumptions] — [arXiv:2008.06618](https://arxiv.org/abs/2008.06618) §4
- Tidal limit "20 g > a ~ size/r_e² → r_e > 1.5×10⁷ m ~ .05 s", for a 0.5 m traveller.
- With R₅ = 50 µm, "r_e ~ .05 s, c₂ ~ 7×10⁷², ℓ ~ 3×10³ ly, γ = ℓ/r_e ~ 2×10¹², E_bin ~ −5×10⁹ kg".
- **[checked]** 1.514×10⁷ m; 0.0505 s; c₂ = 7.08×10⁷²; ℓ = 3,171 ly; γ = 1.98×10¹²; E_bin = −5.19×10⁹ kg (all within 6%).
- Proper traversal time πr_e/c = 0.16 s (abstract: "less than a second").
- Mouth mass r_e c²/G = 2.0×10³⁴ kg ≈ 1.0×10⁴ M_⊙. For r_e = R₅ the paper says "M ~ 10²³ kg"; I get 6.7×10²² kg.
- No exotic matter is added. The negative energy is the 4D Casimir energy of the fermions, which is classical in 5D (§3).
- MM's own fermion-only estimate requires "N_f > 10⁵²", which is too many. That is why they use the 5D model.

### Inferences
- The library's wormhole functions need no correction.
- Useful additions are short:
  - the thin-shell σ(a, M) and ϑ(a, M);
  - the QI throat bound, √(3/4π) l_P/f²;
  - the band bound, (3r₀l_P²/8πf⁴)^{1/3};
  - a Fewster–Eveson QI with a Gaussian or general sampling.
- MM 2021 is the only "human" wormhole in the mainstream literature. Its inputs are ~10⁴ M_⊙ magnetically charged black holes and a 50 µm extra dimension, so it is a consistency demonstration, not an engineering proposal.

### Gaps
- Morris & Thorne (1988, Am. J. Phys.) is not on arXiv and was not opened. Their "5×10⁴¹ dyn/cm²" quote stays as given in the DIRD review.
- I did not extract Fewster's later curved-spacetime QEIs (e.g. arXiv:1208.5399, which was not downloaded). The flat-space bounds above apply only on scales well below the curvature radius, which FR96 make explicit.

---

## Q4. Claimed experiments: White's Eagleworks interferometer (2011/2013) and White et al. 2021 (EPJ C); critiques; other lab claims

### Takeaway
No experiment has reported a measured spacetime distortion from an engineered warp or wormhole metric.
- **The 2011–2014 Eagleworks interferometer:**
  - It aimed at a 10⁻⁷ metric perturbation over 1 cm.
  - Standard GR predicts ~10⁻⁴⁵ from the capacitor energy.
  - White's own model needs a v_s ≈ 0.03c bubble, which costs ~3×10²⁰ kg of negative energy.
  - No result was published. Lee & Cleaver's critique holds, although its own R_S arithmetic is off by ~10²².
- **White et al. 2021** is a numerical, qualitative shape comparison with no quantitative measure and no data. Matching Casimir magnitudes to Alcubierre's equation gives a "bubble" speed of ~10⁻²² m/s.
- **White's 2013 energy chart** cannot be reproduced from the Alcubierre metric. It sits 6–9 orders of magnitude below the exact integral.

### Cited Findings
**White (2011), "Warp Field Mechanics 101", NASA JSC, 100 Year Starship Symposium** [FRINGE in its "boost" reading] — [NTRS 20110015936](https://ntrs.nasa.gov/citations/20110015936) (PDF read)
- Canonical form: ds² = (v_s²f² − 1)[dt − v_s f/(v_s²f²−1) dx]² − … The paper defines "φ = ½ ln|1 − v_s²f²|" and "γ = cosh(φ)".
- Claim: "the field's boost acts on the initial velocity as a scalar multiplier… ⟨v_eff⟩ = γ v_i". Example: "0.1c… boost field with a value of 100… apparent speed of 10c… Alpha Centauri in 0.43 years".
- Interferometer: "633 nm He-Ne laser… small (~1 cm) spherical region… measure space perturbations down to ~1 part in 10,000,000… a 1 cm diameter toroidal test article (something as simple as a very high-voltage capacitor ring) with a boost on the order of 1.0000001".
- **[checked, both ways]**
  - **White's own model.** 4.37 ly/10 = 0.437 yr. A boost of 100 needs v_s = 200 in his own formula, since γ = (u+1/u)/2 with u = √|1−v_s²|. A boost of 1.0000001 needs v_s = 0.0299c (γ ≈ 1 + v⁴/8). An Alcubierre bubble of 0.5 cm radius at 0.0299c with wall ≈ radius needs **−3.4×10²⁰ kg (−3.1×10³⁷ J)** of negative energy.
  - **Standard GR.** A generous capacitor field (Lee–Cleaver's 10⁶ V/m in a 0.1 m sphere, 0.019 J) produces 2GM/(c²r) = **3×10⁻⁴⁵**, about 38 orders below the 10⁻⁷ target. The phase sensitivity of 10⁻⁷ over 1 cm at 633 nm is 0.0099 rad per pass.
  - **The "γ·v_i" rule is not a GR result.** In Alcubierre's metric the interior moves at v_s whatever the craft's prior velocity. [FRINGE]
- White concedes the drive "violates all three energy conditions (strong, weak, and dominant)".

**White (2013), "Warp Field Mechanics 102: Energy Optimization"** (slides) — [NTRS 20130011213](https://ntrs.nasa.gov/citations/20130011213) (PDF read; page 16 rendered and read as an image)
- Chart "Exotic Mass Warp Requirements, 10m diameter, v_apparent = 10c", plotted against "Shell Thickness Fraction (2*R/S)" from 10⁻⁵ to 1. Reference lines: "Jupiter 1.9×10²⁷", "Earth 6.0×10²⁴", down to Voyager.
  - The dU/dt = 0 (plain Alcubierre/Chung–Freese) curve falls from ≈1×10²⁷ kg at 10⁻⁵ to ≈1×10²⁰ kg at 1.
  - Curves with "hyperspace" bulk velocity dU/dt → 0.999999999c fall to about Voyager mass at the thick end.
- The slides say "Changing topology greatly reduces the energy required" and "Oscillate the bubble intensity". They also describe the "White-Juday Warp Field Interferometer… perturbs optical index by 1 part in 10,000,000".
- **[mismatch, computed in White's own metric]**
  - For R = 5 m and v = 10, Alcubierre's exact integral has a floor over all wall thicknesses of **3.7×10²⁸ kg** (σR ≈ 1.2). The capacitance bound with a flat r₀ = 5 m is 5.6×10²⁸ kg.
  - The chart's thick-end value of ~10²⁰ kg is **~10⁸·⁶ below** that floor.
  - Reading the x-axis as Δ/R, PF Eq. 28 gives 5.6×10³³ kg at 10⁻⁵, against ~1.5×10²⁷ on the chart (~10⁶·⁶ below).
  - The axis "2*R/S" is not defined on the slide, and the Chung–Freese extra-dimension reduction has no stress-energy derivation in the slides.
  - Status: not reproducible. White's qualitative point, that thicker walls lower the peak density as σ², is correct GR.

**Eagleworks status** — [NTRS 20140005831](https://ntrs.nasa.gov/citations/20140005831) (2014 intern poster)
- The warp-field interferometry project was still building hardware: "5mm diameter open-air etalon mirrors", a planned "solenoid to induce the needed magnetic field", and an enclosure "to prevent air currents". No result is reported.
- Secondary sources say the results were "inconclusive" — [Wikipedia: Harold G. White](https://en.wikipedia.org/wiki/Harold_G._White) (search snippet only, not opened).

**Lee & Cleaver (2014), "The inability of the White–Juday warp field interferometer to spectrally resolve spacetime distortions"** — [arXiv:1407.7772](https://arxiv.org/abs/1407.7772)
- They assume "a generous 10⁶ V·m⁻¹ static electric field created by a 100 kV potential difference across a 0.1 m air gap": "u_E = 4.4 J·m⁻³", and 8πGu/c⁴ "is 9.2×10⁻⁴³ m⁻²".
- They also give "R_S is 1.92×10⁻²³ m, R_E is 8.76×10⁻¹² m… 1.81×10⁻⁶ arcseconds, which is 3 orders of magnitude below… 10⁻³ arcseconds".
- Air refraction from a 1 °C parcel boundary could "produce occasional spurious interference fringes".
- **[checked]** u_E = 4.43 J/m³ and 9.19×10⁻⁴³ m⁻². **[mismatch]** The field energy in a 0.1 m sphere gives R_S = 3.1×10⁻⁴⁶ m, or 2.3×10⁻⁴⁴ m for a 0.1 m × 10 m cylinder. Their 1.92×10⁻²³ m is ~10²¹–10²³ too large, so their conclusion holds by an even wider margin.

**White, Vera, Han, Bruccoleri & MacArthur (2021), "Worldline numerics applied to custom Casimir geometry generates unanticipated intersection with Alcubierre warp metric", Eur. Phys. J. C 81 (7), 677, doi:10.1140/epjc/s10052-021-09484-z** [SPEC; numerical only]
- Sources: abstract and metadata from the [SCOAP3 record](https://repo.scoap3.org/records/63808); full text through the [Springer article](https://link.springer.com/article/10.1140/epjc/s10052-021-09484-z) via a WebFetch summary. The quotes below come from that summary.
- Abstract:
  - "a micro/nano-scale structure has been discovered that predicts negative energy density distribution that closely matches requirements for the Alcubierre metric";
  - "qualitatively quite similar to a two-dimensional representation";
  - "a toy model consisting of a 1 µm diameter sphere centrally located in a 4 µm diameter cylinder… correlates well";
  - "chip-scale experiments might be explored… a real, albeit humble, warp bubble".
- Method: worldline numerics for a **massless scalar with Dirichlet boundaries**, with no material dispersion or temperature. Grids are 35×35, 50³ and 100×100 with 2,000-loop ensembles. The pillar "enhances negative energy density by a factor of 3–5".
- Sec. 1 gives a transient field "E = √(ħcπ²/(1080ε₀d⁴))", "~0.7 mV" for d = 4 µm.
- A conditional proposed test (Fig. 10b): conductor transit time through the sphere against free space. The paper reports **no correlation coefficient and no measurement**. Data availability: "This manuscript has no associated data."
- Internal inconsistency in the fetched summary: Fig. 3 is described with "sphere radius 4 µm", while the abstract says a 1 µm diameter sphere.
- **[checked, both ways]**
  - **The paper's own formula.** E = 357 V/m. E·d/2 = 0.71 mV matches "~0.7 mV". The 1080 = 3×360 assigns one third of the Casimir energy |u| = π²ħc/(720d⁴) = 1.7×10⁻⁶ J/m³ to ε₀E²/2.
  - **Alcubierre's equation.** A 1 µm bubble (R = 0.5 µm, wall ≈ 0.5 µm) needs a peak |ρ_E| = 5.2×10⁵⁴ v² J/m³. Setting this equal to the cavity's 1.7×10⁻⁶ J/m³ gives v = 5.7×10⁻³¹ c = **1.7×10⁻²² m/s**. Because T⁰⁰ ∝ v²f′², a shape match carries no information about speed. Note also that the Casimir |u| is in the vacuum and the Alcubierre requirement is a source, not a field strength.
- Critiques: I found no peer-reviewed rebuttal.
  - Press criticism (e.g. a [Gigazine summary of Ethan Siegel](https://gigazine.net/gsc_news/en/20211217-warp-bubble)) is secondary. The original was not opened.
  - Searches turned up no follow-up measurement by the Limitless Space Institute up to Oct 2026.

**Other "laboratory" claims**
- I found no peer-reviewed laboratory claim of a warp bubble or wormhole.
- Related but distinct work covered elsewhere:
  - the EmDrive and Q-thruster at Eagleworks, refuted as thermal artefacts (`gravity_programs_public_record.md` §6);
  - squeezed-vacuum sub-vacuum noise (Riek et al. 2015/2017, DIRD review);
  - optical analogue horizons (Philbin 2008, already in the library's `horizon_temperature`).
- Clough, Dietrich & Khan (2024) compute the gravitational waveform of a collapsing warp bubble as a possible observational signature; no detection is claimed (DIRD review).

### Inferences
- Status labels:
  - Eagleworks WFI: unpublished, no result, physically undetectable by ≥38 orders in GR [FRINGE as a test].
  - WFM 102 energy chart: not reproducible from Alcubierre [FRINGE].
  - White et al. 2021: a qualitative numerical analogy; GR content nil until a magnitude is computed [SPEC].
  - Lee–Cleaver: conclusion correct, arithmetic partly wrong.
- These belong in the library only as documented negative examples, e.g. a doc test that the capacitor perturbation is ~10⁻⁴⁵. They should not become physics functions.

### Gaps
- I could not open the White 2021 PDF directly (Springer and SCOAP3 bot walls). The full-text details come from a WebFetch summary.
- JBIS 66, 242 (2013), the journal version of WFM 102, was not obtained, so the chart's axis definition is unknown.
- No primary Eagleworks WFI data or result report was found on NTRS.

---

## Q5. A general tool: metric → stress-energy → energy conditions for any 3+1 metric — algorithm and test cases

### Takeaway
Two algorithms both work and are validated.
- **(a) Shift-only metrics** (unit lapse, flat slices, i.e. Alcubierre, Natário, Lentz, Fell–Heisenberg): SSV's closed forms need only the Jacobian of the shift and the time derivative of K.
- **(b) Any metric:** a nested 4th-order finite-difference Einstein tensor, the Eulerian tetrad, and Hawking–Ellis eigenvalues. At type I points this gives exact NEC/WEC/DEC/SEC without sampling observers. At other points, sample null directions.

This is what Warp Factory (MATLAB, 1 m grids) and warpax (JAX autodiff) do. A Rust version can use dual numbers or the same finite differences.

### Cited Findings
- **Warp Factory's pipeline** ([2404.03095](https://arxiv.org/abs/2404.03095) §3):
  - g_μν on a 4D grid;
  - "fourth-order central finite difference… 1-meter grid spacing";
  - Ricci, then the Einstein tensor, then T = (c⁴/8πG)G;
  - energy conditions by contracting with observer fields ("100 samples" of orientation and "10 samples" of velocity in Fuchs 2024). The Alcubierre figure uses n_observers = 1000.
- **warpax** ([2602.18023](https://arxiv.org/abs/2602.18023)): forward-mode autodiff from ADM variables, Hawking–Ellis classification, and an algebraic type-I check. It reports that single-frame (Eulerian) checks underestimate violations.
- **My implementation** (`me_scripts/metric_to_energy_conditions.py`, ~150 lines, numpy). Validated results:
  - **Schwarzschild** (isotropic, M = 1): max|T_ab| = 6.5×10⁻¹² at r = 3 and 1.7×10⁻¹¹ at r = 6, against curvature 0.11 and 0.028.
  - **Alcubierre** (v = 1, R = 1, σ = 8), Eulerian ρ against −(1/32π)v²f′²(y²+z²)/r²:
    - side wall: −1.59155×10⁻¹ (= −1/2π); NEC minimum −0.796; SEC_E −0.398;
    - 45° point: −7.958×10⁻²; NEC −0.762;
    - front axis: ρ ≈ 0 but NEC −0.677;
    - inside and outside points agree to 6+ digits.
  - **SSV ρ+p̄** matches at all six points: −0.37136 side, −0.39792 at 45°, −0.42441 on axis.
  - **Warp Factory case in SI** (comoving metric exactly as coded, v = 0.1, R = 300 m, σ = 0.015 m⁻¹): ρ(0, 300, 0) = −6.7751×10³⁵ J/m³, equal to the analytic value. The NEC minimum there is −1.48×10³⁷ J/m³.
  - **Natário** (v = 1, R = 1, σ = 4): ρ at (0.6, 0.85, 0) = −4.05272×10⁻², and at (1.0, 0.1, 0) = −1.18582×10⁻¹, both equal to Natário's formula. NEC minima are −1.23 and −0.55.
  - **Morris–Thorne** b = r₀²/r: ρ = p_r = −1/(8πr⁴) and ρ+p_r = −3.8376×10⁻² at r = 1.2. The NEC minimum is −3.8374×10⁻².
  - **Hawking–Ellis** (`hawking_ellis_check.out`):
    - Morris–Thorne is type I, with min(ρ+p_i) equal to the sampled NEC and SEC ρ+Σp = 0 as expected.
    - The Alcubierre side wall is type I: rest-frame ρ = −0.149, min(ρ+p_i) = −0.616. The sampled Eulerian-normalised NEC is −0.796; magnitudes depend on the null normalization, signs do not.
    - The Alcubierre front axis and a Natário wall point are not type I. They need sampling, consistent with warpax's ~4% non-type-I points.
- **Numerical cautions.**
  - Nested FD with step h has truncation (hσ)⁴ and rounding ε/h². I used h = 10⁻³ in G = c = 1 units for σ ~ 4–8, and h = 0.5 m for the SI case with a 67 m wall.
  - WEC minima over boosted observers scale as γ² wherever the NEC fails (−39.4 at 0.99c here). warpax quotes e^{2ζ}.

### Inferences
- For Rust the shift-only closed form is the cheap, exact path, using SSV's sign: ρ = (c⁴/16πG)(K² − trK²) and ρ+p̄ = (c⁴/24πG)(−2ℒ_nK + K² − 3trK²) with K_ij = v_(i,j). Total Eulerian energy, as a mass, is −(c²/32πG)∫|∇×v|²d³x in kg (v in units of c, ∇ in m⁻¹).
- The general FD evaluator is a natural home for type-I classification, using `nalgebra` if it is already a dependency, or a 4×4 eigen-solve otherwise.
- Test cases with published or exact values are listed in Q6.

### Gaps
- I did not test a non-unit-lapse metric (BM "modified time") or the Fuchs shell. Those need the TOV-plus-smoothing construction, and Warp Factory publishes no arrays. Validation would be against figures only.

---

## Q6. Placement: what to add to `gravitics::spacetime`, with test values

### Takeaway
The minimal, high-value additions are:
- two Alcubierre closed forms (local density, minimum energy) and a doc correction to `alcubierre_energy`;
- thin-shell wormhole surface stresses;
- QI throat and band bounds plus the Fewster–Eveson Gaussian QI;
- the SSV shift-drive density, optionally with a general finite-difference stress-energy evaluator.

Fringe items (White's chart, the WFI, the Lentz and Fell–Heisenberg positivity) should be doc examples labelled with their status, not functions.

### Cited Findings
**Catalogue (metric → energy → energy conditions → number to reproduce → experiment status)**

| Item | Closed form (SI) | Energy conditions | Number to reproduce | Experiment |
|---|---|---|---|---|
| Alcubierre ρ_E | −(c⁴/32πG)v²f′²(y²+z²)/r² ([Alcubierre](https://arxiv.org/abs/gr-qc/0009013)) | WEC, DEC, SEC, NEC violated | −1/(2π) geometric at (0,R,0), σR = 8; −6.775×10³⁵ J/m³ for WF case | none |
| Alcubierre total, linear wall | −(v²/12)(R²/Δ+Δ/12)c²/G (PF Eq. 28; in lib) | — | 6.94×10⁶² kg (R = 100, Δ = 100 l_P) | — |
| Alcubierre total, tanh wall | −(v²/36)R²σ c²/G (σR ≫ 1; derived here) | — | 0.376 M_⊙ (R = 100, σ = 2) | — |
| Minimum for flat interior r₀ | −(v²/12) r₀ c²/G ([BM](https://arxiv.org/abs/2102.06824) App.) | — | 1.122×10²⁸ kg (r₀ = 100 m, v = 1) | — |
| Van Den Broeck | T⁰⁰ = (1/8π)(B′²/B⁴ − 2B″/B³ − 4B′/(B³r)) | WEC violated (region II) | E_II,− = −1.38×10³⁰, E_II,+ = 4.87×10³⁰ kg; QI needs α ≈ 10³⁴ | — |
| Natário | ρ = −(v²c⁴/8πG)[3f′²cos²θ + (f′+rf″/2)²sin²θ] | ρ ≤ 0 everywhere; NEC violated | −0.0405272 geometric at (0.6, 0.85, 0) with v = 1, R = 1, σ = 4; total ≈ −(v²/180)R⁴σ³c²/G | — |
| Krasnikov tube | ER Eq. T | WEC violated | T_tt·ε² = −0.04222 (n = 10, δ = 0.01) | — |
| Shift-only drives | SSV ρ, ρ+p̄, ∫ρ d³x/c² = −(c²/32πG)∫ω² kg ([SSV](https://arxiv.org/abs/2105.03079)) | NEC violated (theorem) | ρ+p̄ = −0.37136 (Alcubierre side wall) | — |
| Morris–Thorne throat | τ₀ = c⁴/(8πG b₀²) (in lib) | NEC violated at throat | 4.815×10⁴⁰ Pa at 10 m | — |
| Thin-shell wormhole | σ = −(c⁴/2πGa)√(1−2GM/ac²); ϑ = −(c⁴/4πGa)(1−GM/ac²)/√(…) ([Visser](https://arxiv.org/abs/0809.0927)) | WEC, NEC violated | −1.926×10⁴³ J/m²; −2.693×10²⁷ kg total (a = 1 m) | — |
| Visser cube | −c⁴/(8G) (in lib) | exotic edges | −1.513×10⁴³ J/m | — |
| QI throat bound | r₀ ≤ √(3/4π) l_P/f² ([FR96](https://arxiv.org/abs/gr-qc/9510071)) | — | 7.90×10⁻³² m (f = 0.01) | — |
| QI band bound | a₀ ≤ (3r₀l_P²/(8πf⁴))^{1/3} | — | 1.46×10⁻²¹ m (r₀ = 1 m) | — |
| FE QI, Gaussian | −3ħc/(64π²(cτ)⁴) ([FE98](https://arxiv.org/abs/gr-qc/9805024)) | — | −1.86×10⁻²⁶ J/m³ at 1 ns (half of the lib's −3.72×10⁻²⁶) | squeezed-light QI-type tests (DIRD review) |
| MM human wormhole | r_e = c√(size/20g); E_bin = −(π²/16)R₅²c²/(G r_e ln(r_e/R₅)) ([MM21](https://arxiv.org/abs/2008.06618)) | no added exotic matter (5D Casimir) | 1.514×10⁷ m; −5.19×10⁹ kg | — |
| Lentz / FH "positive" | Eulerian-only ρ ≥ 0 | NEC/WEC violated ([SSV](https://arxiv.org/abs/2105.03079), [C&R 2025](https://arxiv.org/abs/2511.18251), FH's own text) | FH 9.25×10⁴³ J = 5.2×10⁻⁴ M_⊙c² | — |
| Fuchs shell | TOV shell + interior shift ([2405.02709](https://arxiv.org/abs/2405.02709)) | all satisfied (numerical, sampled observers) | M = 4.49×10²⁷ kg = 2.365 M_J; 2GM/c²R₂ = 0.333 | — |
| White WFI / WFM 102 | — | — | needs v_s = 0.030c → −3.4×10²⁰ kg in White's model; GR capacitor 3×10⁻⁴⁵ | no result published; Lee–Cleaver: undetectable |
| White 2021 Casimir | E = √(ħcπ²/1080ε₀d⁴) | — | 357 V/m, 0.71 mV at 4 µm; equivalent Alcubierre speed 1.7×10⁻²² m/s | proposed only, no data |

**Recommended Rust additions** (names follow the module's style, SI units, `metrology` constants). All test values are in `me_scripts/*.out`.
1. `alcubierre_energy_density(speed, radius, sigma, [x, y, z]) -> f64` (J/m³). Tests:
   - the WF case gives −6.7751×10³⁵ at (0, 300, 0);
   - in G = c = 1 with v = 1, R = 1, σ = 8, the point (0, 1, 0) gives exactly −1/(2π), since f′(R) = −σ/2.
2. `alcubierre_minimum_energy(r0, speed) -> f64` = −(v²/12) r₀ c²/G. Test: −1.122×10²⁸ kg at 100 m, v = 1. Doc note: BM's f = min(r₀/r, 1) attains it, and Alcubierre's σR = 8 profile is 2.68× above it.
3. Doc fix on `alcubierre_energy` (no code change):
   - Alcubierre's tanh wall gives −(v²/36)R²σc²/G, which is 2/3 of Eq. 28 under PF's Δ ≈ 2/σ matching. Test: 0.3762 M_⊙ at R = 100, σ = 2.
   - Lobo–Visser's "−v²R²σ" is an estimate 36× too large.
4. `thin_shell_wormhole(radius, mass) -> (f64, f64)`, giving (σ in J/m², surface pressure in N/m). Tests: σ(1 m, 0) = −1.926×10⁴³ J/m²; 4πa²σ/c² = −2.693×10²⁷ kg.
5. `wormhole_throat_bound(f)` = √(3/4π) l_P/f², and `wormhole_band_bound(r0, f)` = (3r₀l_P²/(8πf⁴))^{1/3}. Tests: 7.897×10⁻³² m; 1.461×10⁻²¹ m (r₀ = 1 m), 1.434×10⁻¹⁴ m (10⁵ ly).
6. `quantum_inequality_gaussian(tau)` = −3ħc/(64π²(cτ)⁴). Test: the ratio to `quantum_inequality` is exactly 1/2. Doc note on `quantum_inequality`: Fewster–Eveson tighten the Lorentzian constant to 27/2048, which is 9/64.
7. `shift_drive_density(jacobian: [[f64;3];3]) -> f64` = (c⁴/16πG)(K² − trK²), with K = sym(J). Optionally add `shift_drive_null_sum(jacobian, lie_n_trace_k)` for ρ+p̄. Tests: Alcubierre and Natário values above. ∫ρ against −(1/32π)∫ω² for Natário agrees to 10⁻⁹.
8. Optional, larger: `stress_energy(metric: impl Fn([f64;4]) -> [[f64;4];4], x: [f64;4]) -> [[f64;4];4]` using nested 4th-order finite differences, plus `energy_conditions(t, metric)` → Eulerian ρ, NEC minimum, and Hawking–Ellis type I (ρ, p_i). Tests: Schwarzschild |T| < 10⁻¹⁰; the Morris–Thorne and Alcubierre values above.
9. Doc examples only, not functions:
   - Van Den Broeck: E_II and the α = 10³⁴ finding;
   - Natário's R⁴σ³ cost;
   - MM 2021 parameters;
   - White's interferometer (3×10⁻⁴⁵ against 10⁻⁷) and the WFM 102 mismatch;
   - the Fell–Heisenberg ratio (5.2×10⁻⁴, not 10⁻⁴).

### Inferences
- Items 1–7 are each one to three lines, consistent with the module's existing style of closed forms with published test values.
- Item 8 is the only substantial new code. It is the general tool the brief asked about, and it subsumes the per-metric density functions for testing.
- No dependency is needed. Note that T^a_b is not symmetric: η·T is symmetric only in an orthonormal frame, and the Hawking–Ellis eigenproblem is the generalised η-symmetric one, so a plain symmetric Jacobi solver is not enough. A hand-written 4×4 routine is ~40 lines; the sampled NEC on a Fibonacci sphere needs none.
- Everything here is pure GR/QFT bookkeeping, and nothing has military application.

### Gaps
- The Fuchs shell needs its smoothed TOV profile, and Warp Factory publishes no arrays, so there are no pointwise test values for it. Only M and the compactness are checkable.
- `gravitics` depends only on `geopotential`, `gravitational` and `metrology`, with no `nalgebra` (checked in `crates/gravitics/Cargo.toml`). Item 8's 4×4 eigen-solve would need a small hand-written Jacobi or QR routine, or should be deferred.
