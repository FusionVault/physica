# Puthoff, "Advanced Space Propulsion Based on Vacuum (Spacetime Metric) Engineering" (JBIS 63, 82–89, 2010; arXiv:1204.2184; DIRD 15): equation-by-equation audit, PV versus GR, and `gravitics::spacetime` proposals

**Files beside this note.** Every computed number below comes from [puthoff_metric_checks.py](puthoff_metric_checks.py), and its output is in [puthoff_metric_checks.out.txt](puthoff_metric_checks.out.txt), cited as "§A–§J". It needs only numpy and scipy, and was run with `python3 -P` from the scratchpad.

**Downloads.** Each download is in its own `dl_*` folder and was read only as text or images:

| Folder | Contents |
|---|---|
| `dl_arxiv_1204.2184` | The paper (JBIS typeset PDF) |
| `dl_abs_1204.2184` | Its arXiv abstract page |
| `dl_dird15` | The DIRD 15 PDF |
| `dl_arxiv_gr-qc_9909037` | Puthoff 2002, *Found. Phys.* |
| `dl_arxiv_1012.5264` | Puthoff & Little, the JBIS 2002 text |
| `dl_arxiv_astro-ph_0107316v1`, `dl_abs_astro-ph_0107316` | The withdrawn 2001 version with Ibison |
| `dl_arxiv_physics_0403064` | Puthoff, Davis & Maccone, the Levi-Civita paper |
| `dl_arxiv_gr-qc_9702026` | Pfenning & Ford 1997, TeX source |

**Status labels used throughout:**
- **[EST]** Established physics or an exact identity.
- **[SPEC]** A consistent model, but untested or needing exotic inputs.
- **[FRINGE]** An interpretation without a derivation, usually the UAP-style "craft" narratives.

## Q1. What exactly is the metric, and how does the paper parameterise "engineering" it?

### Takeaway
The 2010 paper contains **no K at all**. It works with a generic static, diagonal metric and only g₀₀ and g₁₁. It "engineers" the metric by assumption: it posits g₀₀ > 1 and |g₁₁| < 1 around a craft. Its only physical parameterisation is an effective vacuum refractive index, n = √(−g₁₁/g₀₀) (eq. 8), presented as a manipulation of μ₀ and ε₀.

The polarizable-vacuum (PV) dictionary appears only by citation (refs [16] and [19]):
- In Puthoff's PV papers, g₀₀ = 1/K and g₁₁ = −K.
- So **n = K, not √K**.
- Remote observers see a light speed of c/K. The locally measured speed is always c.
- K = exp(2GM/rc²) is the vacuum solution.

### Cited Findings
**Bibliographic facts**
- The arXiv record gives JBIS vol. 63, pp. 82–89 (2010), "8 pages, 3 figures", v1 submitted 3 Feb 2012. — [arXiv:1204.2184](https://arxiv.org/abs/1204.2184)
- The typeset PDF ends "(Received 21 October 2010; 21 November 2010)". Ref. 13 was accessed 26 Nov 2010. — [paper, p. 89](https://arxiv.org/abs/1204.2184)
- DIRD 15 is dated 29 March 2010, with ICOD 1 December 2009 and number DIA-08-1003-015. It is "one in a series of advanced technology reports produced in FY 2009" under the AAWSA program. The author's name is redacted. It has 22 references. — [DIRD 15](https://documents2.theblackvault.com/documents/dia/AAWSAP-DIRDs/DIRD_15-DIRD_Advanced_Space_Propulsion_Based_on_Vacuum_Spacetime_Metric_Engineer.pdf)

**JBIS-only material.** The JBIS text adds the following, which a grep of the DIRD OCR does not find. The rest of the DIRD text is the same apart from wording.
- Refs [23] Leonhardt & Philbin and [24] Deardorff, Haisch, Maccabee & Puthoff, "Inflation Theory Implications for Extraterrestrial Visitation", JBIS 58 (2005).
- The "missing time / Rip van Winkle" sentence.
- The "20 minute … only a few minutes" example.
- The UV and soft-X-ray hazard sentences.
- "Accelerated aging of, say, plants in the area of a landed craft". — [paper §4.1, refs](https://arxiv.org/abs/1204.2184) versus [DIRD 15 pp. 6–12](https://documents2.theblackvault.com/documents/dia/AAWSAP-DIRDs/DIRD_15-DIRD_Advanced_Space_Propulsion_Based_on_Vacuum_Spacetime_Metric_Engineer.pdf)

**Every numbered equation** (signature +,−,−,−) — [paper pp. 83–85](https://arxiv.org/abs/1204.2184):
- **(1)** ds² = g_μν dx^μ dx^ν.
- **(2)** ds² = c²dt² − (dx²+dy²+dz²), with g₀₀ = 1 and g₁₁ = g₂₂ = g₃₃ = −1.
- **(3)** ds² = c²dt² − dr² − r²dθ² − r²sin²θ dφ².
- **(4)** "Schwarzschild-type solution" [ref. 10, Logunov & Mestvirishvili p. 76]: ds² = [(1−Gm/rc²)/(1+Gm/rc²)]c²dt² − [(1−Gm/rc²)/(1+Gm/rc²)]⁻¹dr² − (1+Gm/rc²)² r²(dθ²+sin²θdφ²).
- **(5)** Reissner–Nordström [ref. 11]: add +[Q²G/4πε₀c⁴]/[r²(1+Gm/rc²)²] inside both the g₀₀ bracket and the inverted g₁₁ bracket. The paper notes that the charge term has the opposite sign to the mass term, which is "electrogravitic repulsion" [ref. 12].
- **Unnumbered:** proper intervals √(−g_μν) dx^μ and √g₀₀ dt.
- **(6)** v_L^i = √(−g₁₁)dr/(√g₀₀ dt) = c. Local light speed is always c (ref. 13, Klauber).
- **(7)** v_L^e = dr/dt = √(g₀₀/−g₁₁) c, the remote, coordinate speed.
- **(8)** v_L^e = c/n, with n = √(−g₁₁/g₀₀).
- **Unnumbered:** c = 1/√(μ₀ε₀), with μ₀ = 4π×10⁻⁷ H/m and ε₀ = 8.854×10⁻¹² F/m. "The generation of an effective refractive index n = √(−g₁₁/g₀₀) ≠ 1 by technological means can … be interpreted as manipulation of the vacuum parameters μ₀ and ε₀", with the THεμ formalism (ref. 20) as the framework.

**Field strengths.** "The estimated electric or magnetic field strengths required to generate a given refractive index change are given by standard GR theory (the Levi-Civita Effect) and can be found in [19]". No number is given in the paper. — [paper §3.4](https://arxiv.org/abs/1204.2184)

**Model independence.** The approach is "model-independent, i.e., does not depend on knowledge of the specific mechanisms … only that a technology exists that can control and manipulate (i.e., engineer) the spacetime metric". The paper also says that "What is required for implementation is to determine appropriate sources for their generation". — [paper §2, §4.6, §5](https://arxiv.org/abs/1204.2184)

**The PV dictionary (Puthoff 2002)** — [Puthoff 2002, gr-qc/9909037](https://arxiv.org/abs/gr-qc/9909037):

| Eq. | Statement |
|---|---|
| (3) | D = Kε₀E, with K the "dielectric constant of the vacuum" |
| (6) | ε₀ → Kε₀ and μ₀ → Kμ₀ (needed for constant α and the weak equivalence principle) |
| (7) | c → c/K |
| (8) | E = E₀/√K |
| (9) | m = m₀K^{3/2} |
| (10) | ω = ω₀/√K |
| (11) | Δt = Δt₀√K |
| (13) | Δr = Δr₀/√K |
| (15) | The locally measured light speed "renormalizes … to the value c" |
| (17)/(19) | ds² = (1/K)c²dt² − K(dx²+dy²+dz²), so g₀₀ = 1/K and g₁₁ = −K, "an isotropic coordinate system" |
| (20)/(62) | K = exp(2GM/rc²) |
| (54)–(55) | K Lagrangian with λ = c⁴/32πG and f(K) = 1/K² ("as indicated by Dicke") |
| (59) | Field equation sourced by mass, EM energy density and the K-field energy density |

### Inferences
- **Eqs (4) and (5) are exactly Schwarzschild and Reissner–Nordström** in the harmonic radial coordinate r = R − Gm/c², where R is the areal radius. The residuals are below 2×10⁻¹⁵ (§A). The metric is not isotropic: −g₁₁ = (1+x)/(1−x) radially but (1+x)² tangentially, with x = Gm/rc². So the paper's n from eq. (8) is the *radial* index of one coordinate choice. Klauber and Boonserm (refs 13 and 17) are cited precisely because n is coordinate- and direction-dependent. **[EST]**
- **Translating the paper into PV is exact.** Put g₀₀ = 1/K and g₁₁ = −K into eqs (6)–(8) and Table 1 and you recover Puthoff 2002 eqs (7)–(13) term by term (§B). The answers to the brief's three options are:
  - **n = K**, not n = √K.
  - **Remote light speed c/K.** The *local* speed stays c, so "local light speed c/K" is a mislabel of the coordinate speed.
  - **"Engineering" in PV means K < 1.** **[EST, as a mapping]**
- **The paper never parameterises the source.** In GR terms it chooses g_μν and defers T_μν = (c⁴/8πG)G_μν to "the relevant literature [22]". Q3 shows what that T_μν must be.

### Gaps
- Logunov & Mestvirishvili (refs 10/11) and Klauber (ref. 13) were not opened. Eqs (4)–(5) were checked against Schwarzschild and Reissner–Nordström directly, not against the book.

---

## Q2. Table 1: every entry, its formula, values, and whether it is an identity (PV versus GR)

### Takeaway
Table 1 has seven rows, and **it contains no numbers**: every entry is a formula in g₀₀ and g₁₁ plus a qualitative label.
- **Rows 1–5** (time, frequency, energy, spatial measure, light speed) are correct kinematic identities for static observers compared via the timelike Killing vector. Rows 4–5 are coordinate-dependent.
- **Row 6** (mass, m → (−g₁₁/√g₀₀)m) is also an identity, but only of coordinate bookkeeping: it is the coefficient of ½v² in the slow-motion Lagrangian. It is not an invariant, and GR's conserved energy of a body at rest moves the *other* way, as √g₀₀.
- **Row 7** ("gravitational force f(g₀₀,g₁₁)") is qualitative.
- In PV variables the table is exactly Puthoff 2002's PV table.
- **The right-hand ("Spacetime-Engineered Metric") column is exactly the exterior of a negative-mass body.** This holds both in GR, as eq. (4) with m → −|m|, and in PV, as K = e^{−2G|M|/rc²}.

### Cited Findings
**Table 1** [paper p. 84](https://arxiv.org/abs/1204.2184), "as interpreted by a remote (unaltered spacetime) observer". The two columns are headed "Typical Stellar Mass (g₀₀ < 1, |g₁₁| > 1)" and "Spacetime-Engineered Metric (g₀₀ > 1, |g₁₁| < 1)":

| # | Variable | Formula | Stellar-mass column | Engineered column |
|---|---|---|---|---|
| 1 | Time interval | Δt → Δt/√g₀₀ | clocks "run slower" | clocks "run faster" |
| 2 | Frequency | ω → ω√g₀₀ | "red shift toward lower frequencies" | "blueshift toward higher frequencies" |
| 3 | Energy | E → E√g₀₀ | "energy states lowered" | "energy states raised" |
| 4 | Spatial measure | Δr → Δr/√(−g₁₁) | rulers "shrink" | rulers "expand" |
| 5 | Velocity of light | v_L = c → c√(g₀₀/−g₁₁) | "effective v_L < c" | "effective v_L > c" |
| 6 | Mass | m = E/c² → (−g₁₁/√g₀₀)m | "effective mass increases" | "effective mass decreases" |
| 7 | Gravitational "force" | f(g₀₀, g₁₁) | "gravitational" | "antigravitational" |

**Statements in the text around the table** — [paper §§3.1–4.3](https://arxiv.org/abs/1204.2184):
- §3.1 gives a numerical illustration of row 1. An interval "say, ten seconds" remote is judged locally as √g₀₀ dt, "say, 5 seconds".
- §4.3 states that "velocities in general in different coordinate systems scale as does the velocity of light, i.e., v → √(g₀₀/−g₁₁) v".
- §3.6 says antigravity arises if "the signs of the coefficients g₀₀ and g₁₁ could be reversed".

**The PV table** — [Puthoff 2002, final Table](https://arxiv.org/abs/gr-qc/9909037); [Puthoff & Little Tables 1–2](https://arxiv.org/abs/1012.5264):
- Puthoff 2002, for K ≥ 1 (a typical mass): v_L = c/K; m = m₀K^{3/2} ("effective mass increases"); ω = ω₀/√K; Δt = Δt₀√K; E = E₀/√K; L = L₀/√K; F_K ∝ ∇K.
- Puthoff & Little give the same rows twice: Table 1 for K ≥ 1 and Table 2 for "K ≤ 1 (engineered metric)", with "velocity of light > c", "effective mass decreases", "blueshift", "clocks run faster", "higher energy states" and "objects expand".

**Computed checks:**
- **PV substitution.** With g₀₀ = 1/K and g₁₁ = −K, Table 1 gives Δt×√K, ω×K^{−1/2}, E×K^{−1/2}, Δr×K^{−1/2}, v_L = c/K, n = K and m×K^{3/2}. This matches Puthoff 2002 for K = 4, 1, 0.5 and 0.05 (§B).
- **Row 6 is a Lagrangian coefficient.** For L = −mc²√(g₀₀ − (−g₁₁)v²/c²), the coefficient d²L/dv² equals −g₁₁/√g₀₀ to 10⁻⁶ at (0.7, −1.6) and (1.3, −0.8) (§B).
- **Series in x = GM/rc², isotropic coordinates (§C):**

  | Quantity | GR | PV |
  |---|---|---|
  | g₀₀ | 1 − 2x + 2x² − 1.5x³ | 1 − 2x + 2x² − 1.333x³ |
  | −g_ii | 1 + 2x + 1.5x² | 1 + 2x + 2x² |
  | n | 1 + 2x + 1.75x² | 1 + 2x + 2x² |
  | Mass factor | 1 + 3x + 4x² | 1 + 3x + 4.5x² |
  | Clock √g₀₀ | 1 − x + 0.5x² − 0.25x³ | 1 − x + 0.5x² − 0.1667x³ |

  So the PPN parameters are β = γ = 1 in both, as Puthoff 2002 eqs (42)–(46) state for the orders shown. The first differences are at x² in g_ii and at x³ in g₀₀.
- **At a given areal radius R**, the surface redshift factor √g₀₀ is:
  - identical for the Sun, where the difference is below float resolution at GM/Rc² = 2.12×10⁻⁶;
  - 0.20% lower in PV at GM/Rc² = 0.1723 (an illustrative 1.4 M☉, 12 km neutron star) (§E).

### Inferences
- **Rows 1–3 [EST]** are the gravitational redshift in its √g₀₀ forms, measured as noted in Q5. GR and PV agree to O(x²) and differ at O(x³). The paper's left-column labels are correct.
- **Rows 4–5 [EST, coordinate-dependent].** A remote observer has no instrument that reads "coordinate length" or "coordinate light speed" directly. The invariant observables behind them are:
  - Shapiro delay, whose sign follows n − 1;
  - lensing, a deflection toward high n;
  - angular size and redshift.
- **Row 6 [EST as bookkeeping; SPEC as physics].** It is the ½v² coefficient of the coordinate Lagrangian in whatever coordinates are used.
  - In eq. (4) coordinates it is anisotropic beyond first order. Its logarithm is 3x + x³ radially but 3x − x² + x³ tangentially (my expansion).
  - Local inertia is unchanged in both theories: by the equivalence principle in GR, and by PV 2002 eq. (15)-type renormalisation.
  - GR's invariant energy at infinity of a body at rest is mc²√g₀₀. It *decreases* near a mass (binding energy), opposite to the table's "effective mass increases".
  - Puthoff's own dynamics, Puthoff 2002 eq. (57), use m₀K^{3/2} as the coefficient of the coordinate velocity, so the row is internally consistent in PV.
- **Row 7 [EST, qualitative].** For slow bodies the static field is −c²∇ln√g₀₀, so only g₀₀ matters. g₁₁ enters for light, through the factor (1+γ). §3.6's "reverse the signs of g₀₀ and g₁₁" cannot be literal: a sign flip of g₀₀ makes t spacelike. The meaning has to be the sign of g₀₀ − 1 and of −g₁₁ − 1.
- **The engineered column is a negative mass. [EST math; SPEC/FRINGE physics]**
  - With m → −|m| in eq. (4): g₀₀ = (1+x)/(1−x) > 1 and −g₁₁ = (1−x)/(1+x) < 1. With K = e^{−2G|M|/rc²} in PV: K < 1.
  - So every right-column entry, including "antigravitational", is what a negative-mass Schwarzschild exterior already gives. In PV, g₀₀·(−g₁₁) = 1 identically. So "clocks faster" forces "rulers expand" and "n < 1".
- **Internal inconsistencies in the narrative sections [FRINGE as interpretation]:**
  - **§4.2 "large craft with spacious interior appearing … small".** This is the *shrink* effect, −g₁₁ > 1, i.e. the stellar-mass column. PV cannot combine it with g₀₀ > 1, because g₀₀·|g₁₁| = 1.
  - **§4.4 "light trapped within", "high-refractive-index optical fiber".** These need n > 1. The engineered column has n < 1, which acts as a diverging lens: light is deflected *away*, by 4G|M|/(bc²). That is 0.83″ at b = 10 m for the "1 ppm clock" craft of Q3 (§J).

### Gaps
- No observable in the paper is tied to a specific coordinate system, so rows 4–6 cannot be "tested" as written. Each must first be recast as an invariant: Shapiro delay, lensing or Killing energy.

---

## Q3. Energy densities and field strengths needed to change K (or g₀₀) appreciably, versus GR

### Takeaway
**The paper gives no energy density, field strength or propulsion number.** It says only that "the required energy densities predicted by present theory exceed by many orders of magnitude values achievable", and points to ref. [19].

Recomputed here:
1. **A theorem rules out the engineered column with ordinary sources.**
   - In GR, a static, asymptotically flat spacetime with g₀₀ > 1 anywhere must violate the strong energy condition. Ordinary matter and EM fields satisfy it.
   - In PV, K < 1 anywhere with a regular source needs negative mass-energy density at K's minimum.
2. **Weak field: GR and PV give the same required negative mass.** It is |M| = (g₀₀−1)Rc²/(2G), 6.7×10⁹ kg for a 10⁻¹⁸ clock speed-up at R = 10 m. Strong field: the two theories diverge by up to about 70×.
3. **The EM (Levi-Civita) route raises K above 1**, the wrong sign. Its size is 10⁻³⁸ for 1 T over 1 m.
4. **PV's charged-body route to K < 1** (Puthoff & Little App. B) works only for super-extremal charge, Q/M > 8.6×10⁻¹¹ C/kg. Even then the effect at 10 m is about 10⁻²⁴, and the source interior still needs negative energy.

### Cited Findings
**What the papers say:**
- "With regard to perturbation of the spacetime metric, the required energy densities predicted by present theory exceed by many orders of magnitude values achievable with existing engineering techniques." — [paper §1](https://arxiv.org/abs/1204.2184). The same sentence is in [Puthoff & Little §1](https://arxiv.org/abs/1012.5264).
- **Levi-Civita effect in PV**, the paper's ref. [19] — [Puthoff, Davis & Maccone, physics/0403064](https://arxiv.org/abs/physics/0403064):
  - (1) ∇²√K = −(√K/4λ)[½(B²/Kμ₀ + Kε₀E²) − (λ/K²)(∇K)²].
  - (3)/(4) K = α cos²βz with αβ² = 4πGB²/(μ₀c⁴).
  - (6) K ≈ 1 + (2πGB²/μ₀c⁴)[(L/2)² − z²].
  - (7) c′/c ≈ 1/K.
  - (8) τ ≈ (L/c)[1 + 2πGB²L²/(3μ₀c⁴)].
  - Electric case: (10) K = α sech²βz; (13) K ≈ 1 + (2πε₀GE²/c⁴)[(L/2)² − z²]; (15) τ ≈ (L/c)[1 + 2πε₀GE²L²/(3c⁴)].
  - (B-2) a = c²/(B√(4πG/μ₀)) ≈ 3.484×10¹⁸/B m.
  - "For even the strongest fields of interest (e.g., those of pulsars, B ~ 10⁹ Tesla) perturbation of the spacetime metric is sufficiently small".
  - Appendix B says the Levi-Civita metric "possess[es] none of the required characteristics of a traversable wormhole".
- **The charged PV solution** — [Puthoff 2002 eqs (63)–(71)](https://arxiv.org/abs/gr-qc/9909037); [Puthoff & Little App. B, eqs (B4)–(B6)](https://arxiv.org/abs/1012.5264):
  - √K = cosh(√(a²−b²)/r) + [a/√(a²−b²)] sinh(√(a²−b²)/r), with a = GM/c² and b² = Q²G/(4πε₀c⁴).
  - Weak field: g₀₀ = 1 − 2a/r + (2 + b²/a²)(a/r)², matching Reissner–Nordström eq. (67).
  - Puthoff & Little: "For b² > a², however, the hyperbolic solutions turn trigonometric, and K can take on values K < 1." This is the only *mechanism* for K < 1 offered in the PV papers.
- **Field energy densities have opposite signs.** "The energy densities of the EM fields and the K variable enter into Eq. (59) with opposite signs. Therefore, EM field effects can counteract gravitational field effects." — [Puthoff 2002 §IV](https://arxiv.org/abs/gr-qc/9909037)

**Computed values (§F, §G, §H, §J).**

*Levi-Civita constant.* c²/√(4πG/μ₀) = 3.4789×10¹⁸ T·m against the paper's 3.484×10¹⁸. It reproduces as 3.4837×10¹⁸ with c = 3×10⁸ m/s, so the paper rounded c.

*Levi-Civita PV, eq. (6) at the centre:*

| Field | Region | K(0) − 1 | Extra light transit time (eq. 8) |
|---|---|---|---|
| B = 1 T | L = 1 m | 1.03×10⁻³⁸ | 4.6×10⁻⁴⁷ s |
| B = 45 T | L = 1 m | 2.09×10⁻³⁵ | — |
| B = 10⁹ T (the paper's pulsar value) | L = 10 km | 1.03×10⁻¹² | — |
| E = 3×10⁶ V/m | L = 1 m | 1.0×10⁻⁴² | — |
| E = 1.32×10¹⁸ V/m (Schwinger) | L = 1 m | 2.0×10⁻¹⁹ | — |

All of these have **K > 1**: EM energy slows light, like mass.

*Charged route:*
- b > a ⇔ Q/M > √(4πε₀G) = 8.6175×10⁻¹¹ C/kg. This already exists as `gem::extremal_charge_to_mass()`.
- For b/a = 2:
  - K < 1 inside r = 1.654a.
  - √K → 0 at isotropic r = 0.827a. There g₀₀ → ∞ and the areal radius → 0, a naked singularity.
  - In GR Reissner–Nordström, g₀₀ > 1 inside the areal radius R = b²/2a = 2a, with a singularity at R = 0.
  - For the electron, R = kQ²/(2Mc²) = 1.409×10⁻¹⁵ m = r_e/2.
- A 10⁴ kg craft would need Q = 1.41×10⁶ C for K < 1 out to 10 m, giving E(10 m) = 1.27×10¹⁴ V/m. That is below the Schwinger field of 1.32×10¹⁸ V/m, but the metric change is only |g₀₀ − 1| ≈ 1.5×10⁻²⁴.

*Negative mass for a static region of areal radius R = 10 m with g₀₀ > 1 (§G).* GR uses the exterior Schwarzschild solution with M < 0. PV uses K = e^{−2G|M|/rc²} with R = r√K.

| Target | √g₀₀ | GR \|M\| | PV \|M\| | Mean \|ρc²\| (J/m³) | Surface repulsion (m/s²) | Entry speed needed |
|---|---|---|---|---|---|---|
| Optical-clock-detectable, g₀₀ − 1 = 10⁻¹⁸ | 1 + 5×10⁻¹⁹ | 6.73×10⁹ kg | 6.73×10⁹ kg | 1.4×10²³ | 4.5×10⁻³ | v ≥ c√(g₀₀−1) = 0.30 m/s |
| 1 ppm clock speed-up | 1 + 10⁻⁶ | 1.35×10²² kg | 1.35×10²² kg | 2.9×10³⁵ | 9.0×10⁹ | v ≥ 424 km/s |
| "Plants age" 10% faster | 1.1 | 1.41×10²⁷ kg | 1.41×10²⁷ kg | 3.0×10⁴⁰ | 8.6×10¹⁴ | v ≥ 0.417c |
| §4.1 "20 min vs a few min" (5 min) | 4 | 1.01×10²⁹ kg | 7.47×10²⁸ kg | 2.2×10⁴² | 1.7×10¹⁶ | v ≥ 0.968c |
| §4.1 "20 min vs a few min" (3 min) | 6.67 | 2.93×10²⁹ kg | 1.70×10²⁹ kg | 6.3×10⁴² | 2.9×10¹⁶ | v ≥ 0.989c |
| §4.1 300 K Wien peak (9.66 µm) → 550 nm | 17.6 | 2.07×10³⁰ kg (1.04 M☉) | 6.78×10²⁹ kg | 4.4×10⁴³ | 7.9×10¹⁶ | v ≥ 0.99838c |
| §4.1 "external sounds … redshift below the auditory range", 20 kHz → 20 Hz | 1000 | 6.7×10³³ kg | 9.3×10³¹ kg | 1.4×10⁴⁷ | 4.5×10¹⁸ | — |

Reference scales for the table:
- Casimir |u| at 1 µm is 4.33×10⁻⁴ J/m³. The 10⁻¹⁸ case needs 3.3×10²⁶ times that.
- Nuclear matter is about 2×10³⁴ J/m³ (an assumed round input).
- The weak-field central scale is u = δc⁴/(4πGL²): 9.6×10⁴⁰ J/m³ for δ = 1, L = 10 m. Here c⁴/G = 1.2103×10⁴⁴ J/m.

*§3.1's own example* (10 s → 5 s, so g₀₀ = 0.25, a slowed clock with positive mass):
- GR: areal R = (8/3)GM/c², 1.33× the Schwarzschild radius.
- PV: isotropic r = 2GM/(c² ln 4) = 1.443 GM/c², areal radius 2.885 GM/c².

### Inferences
- **The strong-energy-condition theorem [EST in GR].** For a static metric the lapse N = √g₀₀ obeys D²N = 4πN(ρ + p₁ + p₂ + p₃), with G = c = 1. This is the Tolman–Komar identity, cited from standard GR and not from an opened source.
  - With ρ + Σp ≥ 0 everywhere, N is subharmonic. So N ≤ 1, its value at infinity, everywhere.
  - Maxwell fields have ρ + Σp = 2ρ_EM ≥ 0. So **no arrangement of ordinary matter, charges or fields gives g₀₀ > 1 anywhere.**
  - Reissner–Nordström's "electrogravitic repulsion" region (g₀₀ > 1 inside R = kQ²/(2Mc²)) is therefore only reachable if the charged body's interior violates the strong energy condition. Equivalently, the field energy outside r exceeds the total mass-energy, so the enclosed mass is negative.
- **The PV analogue [SPEC model; the argument is mine].** Take PV 2002 eq. (59) in the static case. At an interior minimum of K, ∇K = 0, so ∇²√K = −(√K/4λ)(ρ_m c²/K + u_EM) must be ≥ 0. That requires ρ_m < 0 there. In a vacuum with only EM fields, √K cannot have a positive minimum (eq. 64b), so it falls monotonically to the singular K = 0. **So "K < 1 by charge" relocates the negative energy into the source; it does not remove it.** GR and PV agree on this.
- **GR versus PV for the craft column.**
  - Weak field: identical masses, the β = γ = 1 regime.
  - Strong field: PV needs less negative mass, by 1.4× at √g₀₀ = 4, 3× at 17.6 and 72× at 1000, because PV's g₀₀ is exponential in M.
  - Both are astrophysical masses of *negative* energy packed inside 10 m.
- **The repulsive barrier contradicts the "hardened craft" story [FRINGE claim; EST physics].** It follows from g₀₀ alone, so GR and PV agree.
  - Killing-energy conservation means anything arriving from outside must have γ∞ ≥ √g₀₀ to reach a static point of the region.
  - Water, sound and air at ordinary speeds cannot reach a craft with appreciable g₀₀ − 1. They are turned back by a surface repulsion of 10¹⁰–10¹⁸ m/s².
  - So the §4.1 claims that the craft "could … impact water at high velocities without apparent deleterious effects" and that "external sounds could … redshift below the auditory range" do not fit a metric with g₀₀ > 1: the water and sound never arrive.
  - The §4.1 numbers are also mutually inconsistent. The 20-min/few-min example implies √g₀₀ ≈ 4–7. Silencing all sound needs about 1000. Visible glow from 300 K needs about 17.6.
- **Confinement (zero exterior field) [SPEC].** Hiding the region from outside gravimeters needs a ± mass double layer with dipole density δc²/(8πG). That is 5.4×10⁷ kg/m for δ = 10⁻¹⁸ and 5.4×10¹⁹ kg/m for δ = 10⁻⁶ (§G). The inner layer must be negative mass.

### Gaps
- The full PV field equation was transcribed from OCR text with its square roots lost. My use of eq. (59)/(64b) assumes the reading ∇²√K on the left-hand side. This is confirmed by eq. (64b) reproducing the exponential solution and by my ODE residual of about 10⁻¹⁰ for eq. (65) (§F).
- I did not open Davis 2009 (ref. [22]), where the paper says "the options for meeting such requirements can be found". So I cannot say whether it contains an energy estimate specific to this paper.

---

## Q4. Warp drives (Alcubierre) and wormholes (Morris–Thorne) in the paper, against Pfenning & Ford and the quantum inequalities

### Takeaway
The paper's warp and wormhole content is qualitative. It has no metric, no energy and no quantum inequality. Its qualitative statements about Alcubierre agree with Pfenning & Ford:
- no local superluminal motion;
- no time dilation for the ship;
- negative energy needed in a thin shell.

But it says only that the requirements "are unlikely to be met in the foreseeable future", whereas Pfenning & Ford conclude "physically unattainable". Wormholes appear only as citations [4–6]. **Consistent, but it omits the quantitative no-go.**

### Cited Findings
**What the paper says** — [paper §4.6](https://arxiv.org/abs/1204.2184):
- Alcubierre's metric allows "arbitrarily large apparent faster-than-light speed … without violating the local speed-of-light constraint".
- "the proper (experienced) acceleration along the spaceship's path would be zero, and … the spaceship would suffer no time dilation".
- It needs "a thin-shell or surface layer of exotic matter, i.e., a quantum field having negative energy and/or negative pressure".
- "the technical requirements for such are unlikely to be met in the foreseeable future [22]".
- "arbitrary spacetime metrics … can in principle be postulated. What is required for implementation is to determine appropriate sources".

Wormholes appear only in §4.3: "This is the basis for discussion of warp drives and wormholes in the GR literature [2-6]".

**Pfenning & Ford** — [Pfenning & Ford, gr-qc/9702026](https://arxiv.org/abs/gr-qc/9702026), TeX source, equation numbers by count:
- Alcubierre's spacetime: "the proper time of the space going traveler's trip is identical to that of the elapsed time on earth … the spaceship never locally travels faster than the speed of light".
- Eq. (9), the quantum inequality: (τ₀/π)∫⟨T_μν u^μu^ν⟩dτ/(τ²+τ₀²) ≥ −3/(32π²τ₀⁴).
- Eq. (23): Δ ≤ 10² v_b L_Planck for α = 1/10.
- Eq. (28): E = −(v_b²/12)(R²/Δ + Δ/12).
- Text: "E ≤ −6.2×10⁷⁰ v_b L_Planck ∼ −6.2×10⁶⁵ v_b grams", which is about −3×10²⁰ M_galaxy v_b. A Δ = 1 m wall would be "on the order of a quarter of a solar mass … yet still not attainable".
- The abstract: "physically unattainable".
- Negative energy is similarly confined to a thin band at a wormhole throat (Ford & Roman 1996) and for the Krasnikov metric (Everett & Roman 1997).

**Recomputed P&F eq. (28)** for R = 100 m and Δ = 100 L_P (§I):
- E = −5.156×10³⁵ m = −3.19×10⁷⁰ L_P = −6.94×10⁶² kg (−6.94×10⁶⁵ g).
- The gram figure is 12% above the paper's, and the L_P figure is 1.9× below its "6.2×10⁷⁰ L_P".
- For Δ = 1 m: −0.564 M☉, not "a quarter".
- This confirms the repo's existing check and the `alcubierre_energy` doc comment.

**Exponential-metric wormhole.** "the absence of horizons — it is certainly clear that this geometry does not represent a black hole … instead one is dealing with a traversable wormhole". — [Boonserm et al. 2018, arXiv:1805.03781](https://arxiv.org/abs/1805.03781), abstract.
- My computation: the areal radius r·e^{GM/rc²} has its minimum e·GM/c² at r = GM/c². That is 4013.9 m for 1 M☉ (§E).

### Inferences
- **The paper's Alcubierre statements are correct [EST math]. Its feasibility wording [SPEC] understates the published bound by many orders of magnitude.** The paper cites Davis (ref. 22), not Pfenning & Ford, and states no quantum-inequality bound.
- **Wormholes.** The Morris–Thorne throat tension c⁴/(8πG b₀²) (4.8×10⁴⁰ Pa at b₀ = 10 m, from the crate's `wormhole_throat_tension` test and the repo review) and the Ford–Roman thin-band result are the relevant GR numbers. Both are already in the crate or the repo review. Nothing in the 2010 paper adds a formula.
- **The PV irony [SPEC].** Taken to strong field, PV turns every sufficiently compact *positive* mass into a horizonless traversable wormhole (throat e·GM/c²) with no negative-energy matter. The "exotic" stress in GR terms is supplied by the K-field. This is PV's most distinctive strong-field prediction, and it is testable (Q5). The 2010 paper does not mention it.

### Gaps
- Davis 2009 (ref. 22), Alcubierre 1994, Morris–Thorne 1988, Krasnikov and Everett–Roman were not opened this session. The statements about them are quoted from Pfenning & Ford's text or the existing repo review.

---

## Q5. What can be tested, by what measurement, and how large is the effect?

### Takeaway
Three groups of predictions, each with its effect size:
1. **PV as a gravity theory** is testable and already partly constrained.
   - It equals GR at 1PN, so all classical tests pass.
   - It differs at 2PN: 0.730 µas of extra light bending at the solar limb.
   - Strong field: a shadow 4.6% larger, no horizon, and a neutron-star redshift 0.2% lower at a fixed areal radius.
   - As published it has g₀ᵢ = 0. That implies no frame dragging, against Gravity Probe B's −37.2 ± 7.2 mas/yr (5.2σ from zero).
2. **Metric-engineering effects in the lab** are 10⁻³⁸ to 10⁻¹⁹ (EM routes), with the wrong sign.
3. **The craft signatures** are testable in principle, but only by observing astrophysical amounts of negative mass. The cheapest detector is a gravimeter: a 1 µGal instrument at 10 m resolves a monopole g₀₀ shift of 2.2×10⁻²⁴, about 10⁶ below optical-clock reach.

### Cited Findings
**Classical tests in PV** — [Puthoff 2002 §III](https://arxiv.org/abs/gr-qc/9909037):
- Redshift (eq. 22) "verified … to a precision of 5%" with the sodium D1 line (Brault).
- Pound–Rebka: "a measurement accuracy of Δω/ω ~10⁻¹⁵ for a height h = 22.5 meters".
- Bending (eq. 29) "For the sun this gives 1.75 seconds of arc".
- Perihelion (eq. 40) "43 seconds of arc per century".
- Note: "important topics such as gravitational radiation and frame-dragging effects were not addressed".
- Footnote 18: for "increasingly larger departures from flatness … the two approaches … begin to diverge".

**Gravity Probe B** measured a frame-dragging drift rate of "−37.2 ± 7.2 mas/yr" against a GR prediction of "−39.2 mas/yr". The geodetic rate was −6,601.8 ± 18.3 against −6,606.1 mas/yr. — [Everitt et al., arXiv:1105.3456](https://arxiv.org/abs/1105.3456)

**Levi-Civita "testable modifications for strong-field conditions"** — [Puthoff, Davis & Maccone §1](https://arxiv.org/abs/physics/0403064): the PV approach "reproduces results predicted by GR for standard (weak-field) astrophysical conditions while posing testable modifications for strong-field conditions".

**Computed effect sizes:**

| Test | GR | PV | Difference | Source |
|---|---|---|---|---|
| Light bending, 2nd-order coefficient | 15π/4 | 4π | 0.730 µas at the solar limb (first order 1.7512″) | numerical, 11.7809732 and 12.5663713 (§D; exact to 6×10⁻⁸) |
| Light bending, 3rd-order coefficient | 128/3 = 42.667 (numerical 42.66707) | 48 (numerical 48.0004) | — | numerical (§D). Corrected 2026-10-27: the first quadrature lost three digits near the turning point and gave 42.646 and 47.97; the cancellation-free form recovers 128/3 to 10⁻⁵ and confirms the report's 48 |
| Critical impact parameter | 3√3 GM/c² | 2e GM/c² | PV shadow deviation δ = +0.0463 | §E |
| Photon sphere (areal) | 3GM/c² | 3.297GM/c² | — | §E |
| Surface redshift factor at GM/Rc² = 0.1723 | — | — | −0.20% in PV | §E |
| Pound–Rebka gh/c² | 2.456×10⁻¹⁵ | identical | — | §J |
| Solar surface redshift | −2.1225×10⁻⁶ | identical | — | §J |

The PV deflection decomposition π(2 + 2γ − β + 3δ₂/4) uses δ₂ = 1 for GR and 4/3 for PV. That formula is my consistency check, not taken from a source.

### Inferences
- **The 2PN light-bending difference [SPEC, testable in principle].** 0.73 µas at the solar limb is below any solar-limb astrometry I have a source for. The exact strong-deflection integral is already in `exact_deflection`.
- **The shadow [SPEC, testable now].** PV's δ = +0.046 versus the EHT Sgr A* fiducial δ of −0.08 (VLTI prior) and −0.04 (Keck prior), which are the crate's test values from EHT Paper VI. PV lies 0.09–0.13 above them. Whether that is excluded depends on the EHT error bars, which I did not open (Gap). PV also predicts no horizon, so no event-horizon-based arguments apply. That is a separate, qualitative test.
- **Frame dragging [inference; PV as published: FRINGE-level tension].** The PV metric (eq. 17) has no g₀ᵢ, and the K Lagrangian is scalar. So a rotating Earth produces no Lense–Thirring term unless PV is extended. GP-B excludes zero at 5.2σ. Puthoff himself flags this as unaddressed, so the honest label is "PV is incomplete, not refuted", while any purely scalar completion is disfavoured.
- **Lab metric engineering [EST bound].**
  - Even the Schwinger field over 1 m gives K − 1 = 2×10⁻¹⁹, about 1/5 of today's best fractional clock resolution of 10⁻¹⁸ (an assumed round figure), and with K > 1.
  - Static fields, Casimir cavities or any SEC-satisfying source cannot give g₀₀ > 1 at all, by the theorem in Q3.
  - A Casimir cavity's negative field energy (4.3×10⁻⁴ J/m³ at 1 µm) is outweighed by its plates' positive rest mass. Its net g₀₀ shift is ordinary attraction.
- **Craft signatures [FRINGE, but testable].** Any real g₀₀ > 1 region with an unconfined monopole produces outward gravity c²δ/(2R): 4.5 mm/s² at R = 10 m for δ = 10⁻¹⁸. That exceeds any confinement-free stealth claim, and it is a clean observable for this repo's quantum gravimeters. A gravimeter at 1 µGal resolves δ ≈ 2×10⁻²⁴ at 10 m.

  The other observables:
  - Blueshifted thermal emission at T_obs = √g₀₀ T, with intensity rising as g₀₀².
  - Shapiro time *advance*.
  - Diverging-lens distortion of the background, 0.83″ at 10 m for |M| = 1.35×10²² kg.
  - A repulsive barrier, with entry needing γ ≥ √g₀₀.

  All follow from g₀₀ and g₁₁ alone, identically in GR and PV in the weak field.

### Gaps
- I did not open the EHT Sgr A* Paper VI error bars, LIGO/Virgo polarisation tests (PV's scalar sector would add a breathing mode), or current astrometric precision. So the strength of the shadow and gravitational-wave constraints on PV is not quantified here.
- Puthoff's 1996 *Physics Essays* paper (ref. 3) and Deardorff et al. 2005 (ref. 24, the source of the "craft" framing) were not opened.

---

## Q6. Relation to Puthoff 2002 (*Found. Phys.*) and Puthoff, Little & Ibison (JBIS 55, 137, 2002)

### Takeaway
The 2010 paper is the 2002 PV tables re-expressed in general metric coefficients:
- every K-formula is replaced by its (g₀₀, g₁₁) form;
- harmonic-coordinate Schwarzschild and Reissner–Nordström (eqs 4–5) are used instead of isotropic eqs (41)/(66);
- the "K ≤ 1 engineered" table becomes the right-hand column.

The introduction and conclusions reuse JBIS 2002 text nearly verbatim. Dropped from 2002:
- the propellantless-propulsion constraints (centre-of-mass and hidden momentum);
- the zero-point-energy, gravity and inertia material;
- the only K < 1 mechanism (the super-extremal charged solution).

Added: the UAP-style "craft signature" narratives.

### Cited Findings
**Shared equations and rows.** The PV-to-2010 map is: v_L = c/K → c√(g₀₀/−g₁₁); m = m₀K^{3/2} → (−g₁₁/√g₀₀)m; ω = ω₀/√K → ω√g₀₀; Δt = Δt₀√K → Δt/√g₀₀; E = E₀/√K → E√g₀₀; L = L₀/√K → Δr/√(−g₁₁); F_K ∝ ∇K → f(g₀₀, g₁₁). The "Mahajan, Qadir & Valanju" electrogravitic-repulsion reference is PV 2002's ref. 20 and 2010's ref. 12. — [Puthoff 2002 Table, eqs 7–13, 58](https://arxiv.org/abs/gr-qc/9909037); [paper Table 1](https://arxiv.org/abs/1204.2184)

**Shared text** — [Puthoff & Little](https://arxiv.org/abs/1012.5264):
- The T. D. Lee "vacuum engineering" quote.
- The "flyleaf" sentence.
- The "required energy densities exceed by many orders of magnitude" sentence.
- The "'Blue Sky'" phrase.
- The closing sentence "sample calculations … indicate the direction of potentially useful trends … constrained only by what is achievable practically in an engineering sense".

All appear in JBIS 2002 §§1 and 4, and again in 2010 §§1, 2 and 5. JBIS 2002 §4 explicitly "utilize[s] … the polarizable vacuum (PV) representation" with "K ≈ 1 + 2GM/rc² > 1" (Table 1) and "K ≤ 1 (engineered metric)" (Table 2).

**JBIS 2002 content absent in 2010** — [Puthoff & Little §2–3, App. A–B](https://arxiv.org/abs/1012.5264):
- §2.1: the centre of mass "cannot change its position if not acted upon by outside forces".
- App. A: "hidden momentum" cancels static E×H field momentum.
- ZPE gravity (Sakharov) and ZPE inertia (Haisch–Rueda–Puthoff).
- Casimir.
- App. B: eqs (B1)–(B6) and "K can take on values K < 1".

**Authorship history.** The JBIS 2002 paper was first posted with Ibison as astro-ph/0107316 (v1 17 Jul 2001, v2 28 Nov 2001). It was withdrawn as v3 on 26 Oct 2010 by Ibison, with the comment "Author Ibison does not subscribe to some of the speculations in this document". It was re-posted by Puthoff & Little as arXiv:1012.5264 on 23 Dec 2010. — [astro-ph/0107316 abstract page](https://arxiv.org/abs/astro-ph/0107316); [arXiv:1012.5264](https://arxiv.org/abs/1012.5264)

**PV's own claim of equivalence:**
- "to the order of expansion shown, the two metric tensors are identical" (eqs 42–46).
- For charge, eqs (67)–(68) versus (70)–(71).
- Footnote 18: PV agrees with Einstein–Yilmaz ("exponential solutions") and diverges from GR in strong field. — [Puthoff 2002](https://arxiv.org/abs/gr-qc/9909037)

### Inferences
- The 2010 "model-independent" framing removes the one place where PV made contact with a mechanism, the charged K < 1 solution. That solution is super-extremal and, as shown in Q3, needs negative energy in the source anyway.
- **Correction to the repo review** (`/home/dan/code/FusionVault/quantum-gravimeter/research/topics/fringe-physics-and-dirds/notes/dirds_gravity_and_metric_engineering.md`):
  1. Its Q2 entry "(6)–(8) Coordinate light speed v_L = c√(g₀₀/(−g₁₁))" conflates eqs (6) and (7). Eq. (6) is the *local* speed, = c. Eq. (7) is the coordinate speed.
  2. Its Q3 Puthoff-2002 entry gives the perihelion formula Δθ = 6π(GMm₀/Lc)² as "(39)". In the paper that is **eq. (40)**; eq. (39) is γ ≈ 1 − 3(GMm₀/ℓc)².
  3. Its "do not implement" verdict on Table 1 holds for the narratives. But two parts are worth code:
     - the mass row as the Lagrangian coefficient;
     - the right column as a negative-mass exterior, with its required mass, repulsion and entry threshold (Q7).
  4. Everything else in that review checks out: eq. (4) in harmonic coordinates, the Table 1 formulas, and "K = exp(2GM/rc²) is not in DIRD 15".

### Gaps
- The Kluwer 2002 PV chapter and Dicke 1957/1961 were not opened.
- The Levi-Civita paper's GRG version (vol. 37, 483) was not compared with the arXiv text.

---

## Q7. What should go in `gravitics::spacetime`, and what is already there

### Takeaway
**Already present:**
- `Metric::PolarizableVacuum` (K = exp(2GM/rc²)).
- `light_bending` (4π versus 15π/4).
- `exact_deflection`, with PV capture at 2e.
- `alcubierre_energy` (P&F eq. 28), `quantum_inequality` (P&F eq. 9), the Casimir functions, `wormhole_throat_tension`, `wormhole_edge_tension`, `horizon_temperature`.
- `shadow_diameter`, GR only.
- `gem::extremal_charge_to_mass` (the PV b > a threshold).

**Worth adding (small, each with a source-quoted test):**
1. Isotropic metric components on `Metric`.
2. Three Table 1 helpers.
3. The negative-mass requirement for a given g₀₀ (GR and PV).
4. The PV Levi-Civita delay.
5. The PV critical impact parameter, for shadows.
6. Optionally, the charged PV solution.

**Do not add:** the §4 craft narratives. Keep them as document checks only.

### Cited Findings
All test values below are quoted from opened sources: [Puthoff 2002](https://arxiv.org/abs/gr-qc/9909037), [Puthoff, Davis & Maccone](https://arxiv.org/abs/physics/0403064), [paper](https://arxiv.org/abs/1204.2184) and [Pfenning & Ford](https://arxiv.org/abs/gr-qc/9702026). Computed cross-checks are marked §.

| Proposed item | Formula | Test value and source | Status |
|---|---|---|---|
| `impl Metric { fn g00(self, mass, r) }` | GR ((1−x/2)/(1+x/2))², PV e^{−2x}, x = GM/rc² (isotropic) | Both = 1 − 2x + 2x² (PV 2002 eqs 42/45); differ at x³: −1.5 versus −4/3 (§C). Sun: √g₀₀ − 1 = −GM/Rc² = −2.12×10⁻⁶ (eq. 22, "5%") | EST (GR) / SPEC (PV) |
| `impl Metric { fn spatial(self, mass, r) }` (= −g_ii) | GR (1+x/2)⁴, PV e^{2x} | Both 1 + 2x (eqs 43/46); x² coefficient 3/2 versus 2 (§C) | as above |
| `fn refractive_index(g00, spatial) -> f64` | √(spatial/g00) (paper eq. 8); remote light speed c/n (eq. 7) | PV: n = K exactly (PV 2002 eqs 7, 19); n_PV/n_GR − 1 ≈ x²/4 (§C) | EST identity, coordinate-dependent |
| `fn redshift_factor(g00) -> f64` | √g₀₀ (Table 1 rows 1–3) | Paper §3.1: 10 s → 5 s ⇔ g₀₀ = 0.25. Pound–Rebka gh/c² = 2.46×10⁻¹⁵ at h = 22.5 m (PV 2002 §III.A; §J) | EST |
| `fn effective_mass_factor(g00, spatial) -> f64` | spatial/√g₀₀ (paper §3.5), the ½v² coefficient of the coordinate Lagrangian | PV: K^{3/2} (PV 2002 eq. 9); d²L/dv² check (§B) | EST bookkeeping; doc must say "not invariant; GR's Killing energy goes as √g₀₀" |
| `fn engineered_mass(excess, radius, metric) -> f64` | Negative mass giving g₀₀ = 1 + `excess` at areal `radius`. GR −excess·Rc²/(2G). PV y = ln(1+excess)/2, −y R e^y c²/G | No published value (the paper gives none). Identity tests: weak field GR = PV; √g₀₀ = 17.56 at 10 m → 2.07×10³⁰ versus 6.78×10²⁹ kg (§G) | Math EST; physics SPEC/FRINGE. Doc: needs SEC violation; exterior gravity c²·excess/(2R); entry needs γ ≥ √g₀₀ |
| `fn levi_civita_delay(energy_density, length) -> f64` | Extra light transit (L/c)(4πG u L²/(3c⁴)), u = B²/2μ₀ or ε₀E²/2 (PDM eqs 8, 15); equivalently K(0) − 1 = πGuL²/c⁴ (eqs 6, 13) | Eq. B-2: a = c²/(B√(4πG/μ₀)) ≈ 3.484×10¹⁸/B m (reproduced with c = 3×10⁸); identity K(0) − 1 = L²/(8a²) | EST at this order (PV = GR weak field) |
| `impl Metric { fn critical_impact(self) -> f64 }` (in GM/c²) | GR 3√3, PV 2e (both already inside `exact_deflection`); `shadow_diameter` could take a `Metric` | Crate's existing test: capture just inside 2e; PV δ = 2e/(3√3) − 1 = +0.0463 (§E) versus EHT fiducial −0.08/−0.04 (existing crate test values) | SPEC, testable |
| (optional) `fn charged_pv_dielectric(mass, charge, r) -> f64` | √K = cosh(k/r) + (a/k)sinh(k/r), k = √(a²−b²); trigonometric continuation for b > a (PV 2002 eq. 65; JBIS 2002 eq. B6) | Weak field g₀₀ = 1 − 2a/r + (2 + b²/a²)(a/r)² (PV 2002 eq. 70 = RN eq. 67); ODE (64b) residual about 10⁻¹⁰ (§F); K < 1 iff Q/M > `gem::extremal_charge_to_mass()` = 8.6175×10⁻¹¹ C/kg | SPEC (PV-specific strong field) |

### Inferences
- **Minimal set.** `Metric::g00` and `Metric::spatial`, plus the three free Table-1 helpers, cover the whole paper as executable identities. `engineered_mass` turns the "craft" column into a number with a GR and PV comparison, which is what the user's "compute both ways" rule asks for.
  - The existing private `photon` constants in `exact_deflection` (1/3 and 1/2) can become `Metric::critical_impact` without changing behaviour. That is a refactor, not new physics.
- **Doc labels.**
  - Table-1 helpers: [EST identity; coordinate-dependent rows 4–6].
  - `engineered_mass`: [EST math; requires negative energy (SEC theorem); FRINGE as a craft explanation].
  - PV items: [SPEC: matches GR at 1PN; differs at 2PN and in strong field; no frame dragging as published (GP-B 5.2σ)].
- **Leave out:** "hardened materials", "plant aging", "missing time", "cloaking", UV hazards and "large interior appearing small". These are narratives without a metric or source. Q2/Q3 show several are self-inconsistent: n < 1 cannot trap light, the repulsive barrier means no slow water impact, and PV forbids g₀₀ > 1 together with |g₁₁| > 1.
- Per the repo rule "rigs are examples", any "craft" calibration (the Q3 table) belongs in an example or test, not in library code.

### Gaps
- None of the proposed functions have a test value *published in the 2010 paper itself*, because it prints no numbers apart from μ₀, ε₀, "10 s / 5 s" and "20 min / few min". Their tests rest on Puthoff 2002, Puthoff–Davis–Maccone, Pfenning & Ford, and identities.
- Shapiro advance through a negative-mass region and the double-layer formula δc²/(8πG) are derived here but not taken from a source. Treat them as candidates pending a literature check.
