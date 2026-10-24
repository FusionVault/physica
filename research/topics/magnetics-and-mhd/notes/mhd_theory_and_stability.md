# MHD theory: equations, waves, shocks, equilibria, stability limits, analytic flows, dynamos and turbulence, as testable formulas

Conventions. SI units, CODATA 2018 constants, μ₀ = 1.25663706212×10⁻⁶ H m⁻¹. Temperatures in eV mean k_BT/e. "Computed" means I recomputed the value in `mhd_theory_checks.py` in this folder. It uses only the standard library and checks itself with asserts. I ran it as `python3 -I <path>/mhd_theory_checks.py` from the session scratchpad, outside this folder, and every assert passes. "Published" values are quoted as the source gives them. "Identity" means the check compares a formula against an exact limit of itself. Labels: **[mainstream]**, **[contested]**, **[refuted]**.

Scope. These notes do not re-derive what the `plasmoid` crate already has: Bennett, Pease–Braginskii, Lee scalings, Ostriker, Taylor/BFM/spheromak eigenvalues, FRC, Lundquist number, Sweet–Parker and plasmoid reconnection, radiation losses, Bosch–Hale and Lawson. They only point to it where a new formula plugs in. Example: the Kadomtsev test in Q4 runs on the Bennett profile from `pinch::bennett_*`, and the tearing growth rate in Q4 takes `relaxation::lundquist` as its input.

Source access. The NRL Plasma Formulary (nrl.navy.mil) returned HTTP 403, and the Glasgow and IPGP hosts of the Christensen papers refused connections. Where these notes quote "NRL-style" practical constants, I computed them from SI and say so. PDFs were read with `pdftotext` from the WebFetch cache.

## Q1. Ideal and resistive MHD parameters: Alfvén and magnetosonic speeds, β, pressure and tension, Rm, Ha, Λ, diffusion time, frozen-in flux

### Takeaway
Every Q1 quantity is a one-line closed form, and each has a published test value. Two are exact: PlasmaPy's docstring examples, v_A = 43 173.870 m/s for protons and 21 664.18 m/s for He²⁺ at B = 0.014 T and n = 5×10¹⁹ m⁻³. The rest are observational: coronal loop v_A of 759 and 537 km/s, active-region v_A of 2000–9400 km/s, and solar-wind β at 1 au. The script reproduces all of them. The practical constants 2.18×10¹¹ (v_A, cgs) and 2.28×10⁷ cm (d_i) follow from SI to 0.1%. Frozen-in flux holds when Rm = μ₀σUL ≫ 1. Laboratory sodium dynamos are hard because μ₀σ ≈ 10 s m⁻² for sodium, so Rm = 100 needs UL = 10 m² s⁻¹.

### Cited Findings

**Formula sheet [all mainstream]**

| # | Quantity | Exact form (SI) | Source | Published test value | Computed |
|---|---|---|---|---|---|
| 1 | Alfvén speed | v_A = B/√(μ₀ρ), ρ = n_i m_i + n_e m_e | [PlasmaPy `Alfven_speed` docs](https://plasmapy-rocco8773.readthedocs.io/en/latest/api/plasmapy.formulary.speeds.Alfven_speed.html) | B = 0.014 T, n = 5×10¹⁹ m⁻³, ρ = n(m_p+m_e) → **43 173.870 m/s**. He-4 2+ at the same n → **21 664.18 m/s** | 43 173.87 and 21 664.18 m/s (rel. err < 10⁻⁶) |
| 2 | Practical form | v_A = 2.1812×10¹⁶ B[T]/√(μ n[m⁻³]) m/s = 2.1812×10¹¹ B[G]/√(μ n[cm⁻³]) cm/s, where μ = m_i/m_p | Computed from SI. The NRL Formulary prints 2.18×10¹¹ in cgs, but I could not re-read it (403) | — | 2.1812×10¹¹ |
| 3 | Corona, typical loop | Row 1 with n = 10⁹ cm⁻³, B = 11 G, ρ = 1.67×10⁻¹⁵ g cm⁻³, β = 0.01 | [arXiv:2605.15057](https://arxiv.org/pdf/2605.15057) (search-result summary) | v_A ≈ **759 km/s** outside the loop, **537 km/s** inside (density ×2) | 759.3 and 536.9 km/s |
| 4 | Corona, active-region loops | Row 1 | [arXiv:2106.10884](https://arxiv.org/pdf/2106.10884) (search-result summary) | mean log n_e = 9.2, mean B ≈ 110 G, v_A ≈ **2000–9400 km/s**, mean β ≈ 0.002 | 6027 km/s at the sample means |
| 5 | Fast and slow magnetosonic speeds | v²_f,s = ½[(v_A²+c_s²) ± √((v_A²+c_s²)² − 4v_A²c_s²cos²θ)], c_s² = γp/ρ | standard ideal-MHD dispersion relation (Freidberg, Goedbloed & Poedts; not re-read this session) | — | Identities: θ = 0 gives (max, min)(v_A, c_s). θ = 90° gives √(v_A²+c_s²) and 0 |
| 6 | Plasma beta, magnetic pressure | β = 2μ₀p/B², p_B = B²/2μ₀. Tension force = (B·∇)B/μ₀ | standard. Form as used in [Crews & Turner 2026, arXiv:2608.00291](https://arxiv.org/pdf/2608.00291) (local β = p/p_B, p_B = B_θ²/2μ₀) | — | p_B(1 T) = 397.887 kPa. p_B(12.2 T, SPARC) = 59.2 MPa |
| 7 | Solar wind β at 1 au | Row 6 | [Wilson et al. 2018 ApJS 236, 41, arXiv:1802.08585](https://arxiv.org/pdf/1802.08585) | ~10 years of Wind data, mean (median): β_p,tot = **1.79 (1.05)**, β_e,tot = **2.31 (1.09)**. T_p = 12.7 (8.6) eV, T_e = 12.2 (11.9) eV. A 2023 erratum exists (ApJS 269, 62) and I did not check it | — |
| 8 | β in extreme solar-wind intervals | Row 6 with p = n k_B T_i | [Chen et al. 2014 GRL 41, 8081, Table 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC4459186/) | Low-β set: B = 13.3 nT, n_i = 3.9 cm⁻³, T_i = 1.64 eV, published β_i = **0.0122 ± 0.0059**. High-β set: B = 1.14 nT, n_i = 17.1 cm⁻³, T_i = 4.8 eV, published β_i = **23.2 ± 7.4** | 0.0146 and 25.4 from the table means. Both are inside 1σ; a mean of ratios is not the ratio of means |
| 9 | Magnetic Reynolds number | Rm = μ₀σUL = UL/η, with η = 1/(μ₀σ) | [Stefani, Gailitis & Gerbeth 2008, ZAMM, arXiv:0807.0299](https://arxiv.org/abs/0807.0299), Eq. (7) | Critical Rm for dynamo flows is "usually in the range of 10¹…10³". For sodium μ₀σ ≈ **10 s/m²**, so Rm = 100 needs UL = 10 m²/s, more than 1 m³ of sodium, and at least 100 kW | σ = 7.96×10⁶ S/m, η = 0.10 m²/s |
| 10 | Magnetic diffusion time | t_d = μ₀σl² (generic). For the slowest dipole mode of a sphere of radius r: τ = r²/(π²η) | t_d: [Stefani et al. 2008](https://arxiv.org/abs/0807.0299). Sphere: [Christensen & Tilgner 2004](https://www.maths.gla.ac.uk/~rs/res/B/DynamoScaling/2004ChristensenTilgner.pdf) (search-result summary; host unreachable) | — | Earth core, r = 3.48×10⁶ m, σ = 10⁶ S/m (assumed): **49 kyr**. Illustrative only |
| 11 | Hartmann number | Ha = BL√(σ/(ρν)). Ha² = Lorentz/viscous force ratio | [Bühler et al., KIT MHD group report](https://publikationen.bibliothek.kit.edu/1000097959/38980113) | DCLL blanket duct: L = 0.0808 m (half-width), simulated at **500 ≤ Ha ≤ 8000** | — |
| 12 | Magnetic Prandtl, Lundquist | Pm = νμ₀σ. S = Blσ√(μ₀/ρ) = Ha·√Pm | [Stefani et al. 2008](https://arxiv.org/abs/0807.0299) (they write "S =: Blσ µ0/ρ = Ha Pm"; the √ placement is garbled in the extracted text). Lundquist is already in `plasmoid::relaxation::lundquist` | Earth's core Pm ~ **10⁻⁶** ([Christensen 2009 lecture](https://cpaess.ucar.edu/sites/default/files/heliophysics/resources/presentations/2009_Christensen_planmagdynamo.pdf), slide 20) | Identity: S = Ha·√Pm holds when S = μ₀σv_A L |
| 13 | Elsasser number | Λ = σB²/(2ρΩ), the Lorentz/Coriolis ratio. Some authors drop the 2 | [Christensen 2009 lecture](https://cpaess.ucar.edu/sites/default/files/heliophysics/resources/presentations/2009_Christensen_planmagdynamo.pdf), slides 20 and 33 | Earth: **Λ = 0.3–5**. Models 0.03–100. Magnetostrophic rule Λ ~ 1, i.e. B²/2μ₀ ≈ ρηΩ | Λ = 2.49 for σ = 10⁶ S/m, B = 2 mT, ρ = 1.1×10⁴ kg/m³ (assumed inputs) |
| 14 | Frozen-in flux (Alfvén's theorem) | For Rm → ∞, ∂B/∂t = ∇×(u×B) and the flux through any material surface is constant. The field grows when the kinematic time l/v is shorter than t_d, i.e. Rm > Rm_c | [Stefani et al. 2008](https://arxiv.org/abs/0807.0299), Eqs. (6)–(7) | — | — |

**Practical test cases (what a calculator should reproduce)**
- Earth's outer core: Rm ≈ **10³**, E ≈ 10⁻¹⁴, Pm ≈ 10⁻⁶, Ro ≈ 5×10⁻⁶, Λ ≈ 0.3–5. The rms field at the top of the core in harmonic degrees 1–13 is **0.39 mT (3.9 G)**. The internal field is "~ 1 - 4 mT ?" (the question mark is in the original) — [Christensen 2009 lecture, slides 5, 20](https://cpaess.ucar.edu/sites/default/files/heliophysics/resources/presentations/2009_Christensen_planmagdynamo.pdf). Rm ≈ 1000 is also quoted in [arXiv:1307.6048](https://arxiv.org/pdf/1307.6048) (search-result summary).
- Torsional-wave constraint on the field in the core: travel time ~6 yr and B_s ≳ 2 mT, attributed to Gillet et al. 2010 and 2015 — [Jault, 1st Swarm Science Meeting proceedings](https://ftp.space.dtu.dk/data/magnetic-satellites/Swarm/SCARF/1st_Swarm_Intl_Sci_Mtg_Nantes/papers/16-Jault-Proc.pdf) (search-result summary; I did not read Gillet 2010 itself).
- Present-day dipole decay: the axial dipole is ~9% weaker than in 1840, falling ~16 nT/yr. Most of the recent decay comes from meridional flux advection by a core gyre, with diffusion a smaller share — [Finlay, Aubert & Gillet 2016, Nat. Commun. 7, 10422](https://www.nature.com/articles/ncomms10422).
- Tokamak core, ITER inductive reference (Q3 sources): n = 1.01×10²⁰ m⁻³, D–T (2.5 u), B = 5.3 T. **Computed** v_A = 7.30×10⁶ m/s and d_i = 3.6 cm.

### Inferences
- A Rust `mhd::params` module needs about 12 functions, each a one-liner. Test it against rows 1, 3, 4 and 8 exactly as the script does. The PlasmaPy docstring values are the cleanest regression targets, because PlasmaPy publishes the exact density convention ρ = n_i m_i + Z n_i m_e. Bake that convention in, or a 10⁻⁵-level mismatch appears for He²⁺.
- Rm and S both measure "advection or Alfvén crossing versus diffusion". S = Rm evaluated at u = v_A, so the existing `plasmoid::relaxation::lundquist` can be reused once the new module exposes `magnetic_diffusivity(σ)`.
- β_p ≈ β_e ≈ 1 in the median solar wind (row 7). Neither the low-β (cold plasma) nor the high-β limit of the magnetosonic formula applies there. Use the full row-5 form.

### Gaps
- I could not re-read the NRL Plasma Formulary 2019 or 2023 (HTTP 403 from nrl.navy.mil; a PSFC mirror failed on TLS). Its printed constants are therefore checked only against my own SI recomputation (2.1812×10¹¹ and 2.2771×10⁷).
- I retrieved no primary source for Earth's core conductivity (σ ~ 10⁶ S/m), outer-core density or flow speed (≈ 0.5 mm/s). The Earth Rm, Λ and decay-time values in the script use assumed inputs and are labelled illustrative. The published ranges they are compared against (Rm ~ 10³, Λ = 0.3–5) come from the Christensen lecture.
- I found no published numeric test for the fast or slow magnetosonic speed. Only identity checks are in the script.
- I retrieved no primary source for liquid sodium's viscosity, so I give no sodium Pm or Ha test case. The ~10⁻⁶ Pm figure is for Earth's core (row 12).

## Q2. MHD shocks and discontinuities: Rankine–Hugoniot jumps for a perpendicular shock

### Takeaway
For a perpendicular shock (B tangential), the density compression X = ρ₂/ρ₁ = B₂/B₁ = u₁/u₂ is the positive root of a quadratic in the upstream Alfvén Mach number M_A and β. I derived the quadratic from the conservation laws, and it matches the published form. In the script it satisfies the momentum and energy jumps to 10⁻¹², tends to gas dynamics as β → ∞, and reaches the strong-shock limit X → (γ+1)/(γ−1) = 4. Two worked bow-shock-like examples, both **computed**: M_A = 8, β = 1, γ = 5/3 gives X = 3.615 and a downstream β of 6.24. At the Wind-median β ≈ 2.1, X = 3.48.

### Cited Findings
- **Perpendicular-shock jump condition [mainstream]**, in the published form: 2(2−γ)X² + γX{2 + β[2 + (γ−1)M_s²]} − γ(γ+1)βM_s² = 0, where M_s is the sonic Mach number and β = 2μ₀p₁/B₁² — [Trinchieri et al. 2003 (Stephan's Quintet), arXiv:astro-ph/0302590](https://arxiv.org/pdf/astro-ph/0302590) (from the search-result summary). With M_A²/M_s² = γβ/2 (stated in the same search set), it is equivalent to my derivation, which uses M_A directly:
  **(2−γ)X² + [γ(1+β) + (γ−1)M_A²] X − (γ+1)M_A² = 0** (computed). Derivation: take mass flux ρu, flux freezing uB, momentum ρu² + p + B²/2μ₀ and energy ρu(u²/2 + γp/((γ−1)ρ)) + uB²/μ₀, then divide out the trivial root X = 1.
- The general oblique case is a cubic in X with coefficients in M_A, β and θ_Bn. Its constant term is ∝ cos²θ_Bn, so at θ_Bn = 90° it factors into X times the quadratic above — [search-result summary of CME-shock papers, e.g. arXiv:1503.04141](https://arxiv.org/pdf/1503.04141). The factorisation step is the summary's own check, not a statement in the sources.
- **Earth's bow shock, observed:** in the MMS bow-shock database (machine-learning compiled), Alfvén Mach numbers are mostly up to about 40. A few crossings reach ~150 when the upstream field is very weak, with large uncertainties. Measured magnetic compression ratios follow the Rankine–Hugoniot trend: more perpendicular shocks compress more, and the ratio tends to **4** at large M_A — [Lalti et al. 2022, arXiv:2203.04680](https://arxiv.org/pdf/2203.04680) (from the search-result summary; the full PDF exceeded the fetch limit).
- Upstream state for a worked example: median solar-wind β_p ≈ 1.05 and β_e ≈ 1.09 at 1 au — [Wilson et al. 2018](https://arxiv.org/pdf/1802.08585). I used β ≈ 2.1 as an illustrative total (a sum of medians, not the median of the sum).

**Computed values** (`perp_shock_X`, `conservation_residual`)

| M_A | β₁ | γ | X = ρ₂/ρ₁ = B₂/B₁ | Check |
|---|---|---|---|---|
| 8 | 1.0 | 5/3 | **3.6155** | Momentum and energy flux residuals < 10⁻¹² |
| 8 | 2.1 | 5/3 | **3.4834** | — |
| 10⁵ | 1.0 | 5/3 | → 4.0000 | Strong-shock limit (γ+1)/(γ−1) |
| M_s = 3 at β = 10⁹ | — | 5/3 | = (γ+1)M_s²/((γ−1)M_s²+2) = 3.0 | Gas-dynamic limit, rel. err < 10⁻⁶ |

### Inferences
- A Rust `shock::perpendicular(m_a, beta, gamma) -> f64` is a stable closed form. Use the quadratic with the "+√" root. Guard the input so that the fast Mach number M_f = M_A/√(1 + γβ/2) is above 1; otherwise no compressive shock exists. Downstream p₂ comes from the momentum jump, and the script returns it.
- Real bow-shock compression above 4 (reported in an unreviewed student MMS analysis turned up by the same search) would be a kinetic, non-steady or measurement effect. It is not a failure of the formula, which assumes a steady, planar, single-fluid shock with isotropic pressure.
- Oblique shocks need the cubic. The script does not implement it.

### Gaps
- I did not retrieve a primary bow-shock event with upstream M_A, β, θ_Bn and a measured density jump. A real event-level test value is therefore missing; the worked numbers are illustrative.
- I did not obtain Priest's or Goedbloed & Poedts' textbook derivation. The published form above comes from a search-result summary of one paper.
- Switch-on and switch-off shocks, intermediate (rotational) discontinuities, and the coplanarity theorem were not covered.

## Q3. Equilibria: Grad–Shafranov and Solov'ev; safety factor and Kruskal–Shafranov; Troyon, Greenwald, IPB98(y,2), with ITER and SPARC tests

### Takeaway
The ITER Physics Basis edge-q formula reproduces three published values to within 1%: ITER q95 = 3.0, SPARC Uckan q* = 3.05, and SPARC "≈ 3.6" with the aspect-ratio factor. The Greenwald limit reproduces ITER's ⟨n_e⟩ = 1.01×10²⁰ m⁻³ at f_G = 0.85, and gives SPARC 3.15×10²⁰ at f_G = 0.37. IPB98(y,2) cannot be tested forward from the retrieved sources, because neither gives the net heating power. Inverted, it needs P_net = 83.8 MW to give ITER's published 3.7 s, and P_loss = 29 MW to give SPARC's 0.77 s. Both are physically plausible. The Cerfon–Freidberg Solov'ev solution, solved in the script for an ITER-like shape, matches the target boundary to 1.2×10⁻⁵ of ψ_axis.

### Cited Findings

**Formula sheet**

| # | Quantity | Form | Source | Published test value | Computed |
|---|---|---|---|---|---|
| 1 | Grad–Shafranov [mainstream] | Δ*ψ ≡ R∂_R(R⁻¹∂_Rψ) + ∂_z²ψ = −μ₀R² dp/dψ − F dF/dψ. Normalised (x = R/R₀, y = z/R₀) with Solov'ev profiles: **Δ*ψ = (1−A)x² + A** | Cerfon & Freidberg, Phys. Plasmas 17, 032502 (2010), cited via [JPP reference list](https://resolve.cambridge.org/core/journals/journal-of-plasma-physics/article/solution-of-gradshafranov-equation-by-the-method-of-fundamental-solutions/8906F24CB0B95B8ED1C17B55B01FDD59) and [APS DPP 2009 abstract](https://archive.aps.org/dpp/2009/up8/99). The paper itself was not retrieved (paywall) | — | Finite-difference residual of the particular solution < 10⁻⁵ for A = 0, −0.155, 0.7 |
| 2 | Solov'ev / Cerfon–Freidberg solution [mainstream, exact] | ψ = x⁴/8 + A(x² ln x/2 − x⁴/8) + Σ₁⁷ cᵢψᵢ. ψ₁ = 1, ψ₂ = x², ψ₃ = y² − x² ln x, ψ₄ = x⁴ − 4x²y², ψ₅ = 2y⁴ − 9y²x² + 3x⁴ ln x − 12x²y² ln x, ψ₆ = x⁶ − 12x⁴y² + 8x²y⁴, ψ₇ = 8y⁶ − 140y⁴x² + 75y²x⁴ − 15x⁶ ln x + 180x⁴y² ln x − 120x²y⁴ ln x. Up-down symmetric conditions: ψ = 0 at (1±ε, 0) and (1−δε, κε); ψ_x = 0 at the top point; ψ_yy = −N₁,₂ψ_x at the equatorial points; ψ_xx = −N₃ψ_y at the top. N₁ = −(1+α)²/(εκ²), N₂ = (1−α)²/(εκ²), N₃ = −κ/(ε cos²α), with α = arcsin δ | Cerfon & Freidberg 2010, as above. The polynomials and N's are from memory and were **not re-read**, but each ψᵢ is verified numerically to satisfy Δ*ψᵢ = 0, and each N is verified against the boundary curve x = 1 + ε cos(τ + α sin τ), y = εκ sin τ | Follow-up work (Solov'ev profiles force an edge current jump; quadratic-profile extension) — [Cerfon, Guazzotto & Freidberg, EPS 2010 P4.115](https://lac913.epfl.ch/epsppd3/2010/pdf/P4.115.pdf) | ITER-like ε = 0.32, κ = 1.7, δ = 0.33, A = −0.155: c₁…c₇ = 6.665e−2, −1.955e−1, −5.111e−2, −4.597e−2, 5.532e−3, −5.531e−3, −1.480e−4. Axis at x = 1.0512 (Shafranov shift 0.051R₀). max\|ψ_boundary\|/\|ψ_axis\| = **1.16×10⁻⁵** |
| 3 | Kruskal–Shafranov [mainstream] | Cylinder of length L (period 2πR₀): q(a) = 2πaB_z/(LB_θ). External kink unstable for q(a) < 1 | [arXiv:physics/0411089](https://arxiv.org/pdf/physics/0411089); [Oz et al. 2011, PPPL](https://w3.pppl.gov/%7Ehji/paper/Oz11.pdf) (line-tied ropes kink at q_a = 1); [Fitzpatrick, UT Austin fusion notes](https://farside.ph.utexas.edu/teaching/plasma1/Fusionhtml/node13.html): q_a > 1 is the ideal-kink criterion, but tokamaks operate safely only for q_a ≳ 3 (all from search-result summaries) | Spheromak formation by gun-driven kink "agree[s] quantitatively with the Kruskal-Shafranov limit" (already in the plasmoid notes, arXiv:physics/0304104) | Practical constant: 2π/(μ₀·10⁶) = 5.0000, so q_cyl = 5a²B/(RI[MA]) |
| 4 | Edge q with shaping (ITER Physics Basis) [mainstream, empirical fit to equilibria] | q* (Uckan) = (5a²B/(RI))·[1 + κ²(1 + 2δ² − 1.2δ³)]/2. q95 = q*·(1.17 − 0.65ε)/(1 − ε²)² | Formula from memory of ITER Physics Basis 1999. Use as named in [SPARC overview, Creely et al. 2020, JPP 86, 865860502](https://www.cambridge.org/core/journals/journal-of-plasma-physics/article/overview-of-the-sparc-tokamak/DD3C44ECD26F5EACC554811764EF9FF0) | ITER reference **q95 = 3.0** ([Sips 2005, arXiv:physics/0410263, Table I](https://arxiv.org/abs/physics/0410263)). SPARC: "q* = 3.05" (Uckan definition, "corresponds to q95≈3.4 on SPARC"), and "the formula from the ITER Physics Basis gives q* ≈ 3.6" | ITER (a 2.0, R 6.2, B 5.3, I 15, κ 1.7, δ 0.33): **2.996**. SPARC (a 0.57, R 1.85, B 12.2, I 8.7, κ 1.75, δ 0.45): Uckan **3.059**, IPB **3.621** |
| 5 | Greenwald density limit [mainstream, empirical] | n_G[10²⁰ m⁻³] = I_p[MA]/(πa²[m²]) | [Sips 2005](https://arxiv.org/abs/physics/0410263) | ITER: ⟨n_e⟩ = **10.1×10¹⁹ m⁻³ at ⟨n_e⟩/n_GW = 0.85**. SPARC: "a Greenwald fraction of only 0.37" ([Creely 2020](https://www.cambridge.org/core/journals/journal-of-plasma-physics/article/overview-of-the-sparc-tokamak/DD3C44ECD26F5EACC554811764EF9FF0)), with ⟨n_e⟩ ≈ 3×10²⁰ m⁻³ in the abstract | n_G(ITER) = 1.194×10²⁰ → 1.015×10²⁰ at 0.85. n_G(SPARC) = 8.52×10²⁰ → 3.15×10²⁰ at 0.37 |
| 6 | Troyon β limit [mainstream; ideal-MHD numerical fit] | β_t[%] ≤ β_N·I[MA]/(a[m]B[T]), with β_N ≈ 2.8 without a wall | Troyon et al. 1984, PPCF 26, 209 (original not retrieved). β_N ≤ 2.8 no-wall; ~4–6 with a close wall plus feedback and rotation, "not sufficiently robust to be considered 'standard'" — [arXiv:2008.10485](https://arxiv.org/pdf/2008.10485) (search-result summary) | ITER reference **β_N = 1.8** ([Sips 2005](https://arxiv.org/abs/physics/0410263)). SPARC "roughly 1.0" ([Creely 2020](https://www.cambridge.org/core/journals/journal-of-plasma-physics/article/overview-of-the-sparc-tokamak/DD3C44ECD26F5EACC554811764EF9FF0)) | ITER β_t = 2.55% (Troyon ceiling 3.96%). SPARC β_t = 1.25% (ceiling 3.50%) |
| 7 | IPB98(y,2) H-mode confinement [mainstream, empirical; revision in progress] | τ_E = 0.0562 H₉₈ I_p^0.93 B_T^0.15 P^−0.69 ⟨n_e,19⟩^0.41 R₀^1.97 M^0.19 κ_a^0.78 ε^0.58. Units: s, MA, T, MW, 10¹⁹ m⁻³ line-average, m, amu. κ_a = V/(2π²R₀a²). P = net heating power (input minus radiation) | [Sips 2005, Eq. (1)](https://arxiv.org/abs/physics/0410263); same form in [Creely 2020, Eq. 4.1](https://www.cambridge.org/core/journals/journal-of-plasma-physics/article/overview-of-the-sparc-tokamak/DD3C44ECD26F5EACC554811764EF9FF0) (P_loss = power across the separatrix by transport, excluding radiation; κ₉₅ ≈ κ_a = 1.75) | ITER: **3.7 s** "for 40 MW additional heating", one technical s.d. ±14%, 95% interval ±28% ([Sips](https://arxiv.org/abs/physics/0410263)). SPARC: **τ_E ≈ 0.77 s**, Q ≈ 11, P_fus ≈ 140 MW at H₉₈ = 1; Q > 2 even at H₉₈ = 0.7 ([Creely](https://www.cambridge.org/core/journals/journal-of-plasma-physics/article/overview-of-the-sparc-tokamak/DD3C44ECD26F5EACC554811764EF9FF0)) | Inverse: ITER needs P_net = **83.8 MW** (κ_a = 1.7 and M = 2.5 assumed). SPARC needs P_loss = **29.0 MW** (n from f_G·n_G, M = 2.5 assumed) |

**Device parameters used** (published)
- ITER inductive reference — [Sips 2005, Table I](https://arxiv.org/abs/physics/0410263): R₀ 6.2 m, a 2.0 m, B_T 5.3 T, I_p 15 MA, q95 3.0, H₉₈ 1.0, β_N 1.8, ⟨n_e⟩ 10.1×10¹⁹ m⁻³ (0.85 n_GW), ⟨T_i⟩ 8.0 keV, ⟨T_e⟩ 8.8 keV, P_NB 33 MW + P_RF 7 MW, P_fus 400 MW, Q = 10, 400 s burn. Hybrid: 12 MA, q95 4.1, H 1.2, β_N 2.1, Q = 11. Reversed shear: 9 MA, q95 5.3, H 1.57, β_N 2.95, Q = 6. ITER's plasma volume is "more than 800 m³".
- SPARC — [Creely et al. 2020](https://www.cambridge.org/core/journals/journal-of-plasma-physics/article/overview-of-the-sparc-tokamak/DD3C44ECD26F5EACC554811764EF9FF0): B₀ 12.2 T, R₀ 1.85 m, a 0.57 m, ε 0.31, κ_sep 1.75 and δ_sep 0.54 (Fig. 1 caption, as extracted), κ₉₅ ≈ 1.75, δ₉₅ ≈ 0.45, I_p 8.7 MA, up to 25 MW ICRH, ~1 MW ohmic, 10 s flattop, ⟨n_e⟩ ≈ 3×10²⁰ m⁻³, ⟨T_e⟩ ≈ 7 keV, ~7 MW/m³ fusion power density.
- IPB98 revision: the ITPA "standard set" now has 5718 entries from 18 machines — [IAEA FEC contribution](https://conferences.iaea.org/indico/event/151/contributions/6399/) (search-result summary).

### Inferences
- The q95 formula is the strongest Q3 test. It reproduces three independent published numbers from two machines, so a Rust `q95(a, r, b, ip, kappa, delta)` can be shipped with these as unit tests. SPARC's 3.4 is presumably from an equilibrium calculation, not from this formula. Do not test against it.
- IPB98(y,2) needs P_net or P_loss, which neither retrieved source tabulates. A Rust test can still pin the formula's arithmetic by checking the inverse powers above (83.8 MW and 29.0 MW, both inside the physically allowed band of aux + alpha + ohmic power). It can also check that the ratio of τ_E between two inputs equals the product of the power-law ratios.
- The Cerfon–Freidberg solve is a 7×7 linear system. The boundary misfit of 10⁻⁵ is strong evidence that my recalled polynomials are right: a wrong ψᵢ would fail the Δ*ψᵢ = 0 check, and a wrong N would fail the curvature check. A Rust crate should still compare the cᵢ with the published paper before release.
- Kruskal–Shafranov ties the tokamak to the `plasmoid` crate. The spheromak and gun-driven kink results already in the plasmoid notes use the same q(a) = 1 threshold. In a tokamak, q95 ≈ 3 is the operating margin above q_a = 2, the external-kink disruption limit that tokamak notes usually quote. q_a = 2 is my inference and is not sourced here.

### Gaps
- The Cerfon & Freidberg 2010 paper was not retrieved, so the published cᵢ for the ITER-like case (and the A value itself) are unverified. The script's coefficients are self-consistent but have no external check.
- Neither source tabulates the ITER Q = 10 net heating power or radiated power, nor SPARC's P_loss (Creely's Table 2 was not in the extracted text), so I have no forward IPB98 test.
- I did not retrieve JET or DIII-D record parameters (β_N records, Greenwald-fraction records), nor the original Troyon 1984 and Greenwald 1988 papers.

## Q4. Instabilities: Z-pinch kink and sausage, magnetic Rayleigh–Taylor, interchange, tearing (FKR), resistive wall modes

### Takeaway
Three instability results have closed forms with published numbers, and the script reproduces all three. (1) The Kadomtsev energy-principle criteria on the Bennett profile give the published marginal m = 1 surface r* = √3 r_p with β* = 1/3. They also show that the Bennett profile is m = 0 unstable everywhere for γ = 5/3, because its current-concentration index α = 2 exceeds γ. (2) The magnetic Rayleigh–Taylor cut-off λ_c = 4πB²cos²θ/(μ₀gΔρ) reproduces Khomenko et al.'s 38 km at θ = 89° and 155 km at θ = 88°, for a 10 G prominence–corona interface. (3) The FKR constant-ψ tearing coefficient [Γ(1/4)/(2πΓ(3/4))]^{4/5} = 0.5474 matches the quoted ≈0.55. Resistive wall modes have no closed form in the sources retrieved.

### Cited Findings
- **Kadomtsev m = 0 (sausage/interchange) criterion [mainstream]:** −d ln p/d ln r ≤ 4γ/(2 + γβ), with local β = p/p_B and p_B = B_θ²/2μ₀. Equivalently, d/dr ln(p/⟨j_z⟩^γ) ≤ 0, where ⟨j_z⟩ is the area-averaged current density. For the polytropic diffuse family p/p₀ = (⟨j_z⟩/j₀)^α, the column is stable iff **α ≤ γ**. The marginal α = γ profile is the "Kadomtsev pinch". The Bennett profile p/p₀ = (1 + (r/r_p)²)⁻² has α = 2 and is m = 0 unstable for γ = 5/3 — [Crews & Turner 2026, arXiv:2608.00291](https://arxiv.org/pdf/2608.00291), Eqs. (15)–(19).
- **Kink (m ≥ 1) criterion [mainstream]:** −d ln p/d ln r ≤ m²/β. For m = 1 it becomes d(r p_B)/dr < 0, which "no profile with a free axis" satisfies everywhere. Compressibility gives no stabilisation, so the kink "is often the one that causes disruptions". Bennett: **β* = 1/3, unstable core out to r* = √3 r_p**. A sharp pinch has β* = 0. "Hard-core" profiles (p_B ∝ 1/r around a rigid conductor) are marginal with β* = 1. A conducting wall at r_w stabilises only k < k_c, with k_c r_p ~ (r_w/r_p)^−1.1 — [Crews & Turner 2026](https://arxiv.org/pdf/2608.00291), Eqs. (20)–(24).
- Sheared axial flow: "trans-Alfvénic sheared flow is shown to stabilize the MHD kink but also to excite shear-driven instabilities", and these shear-driven modes dominate — [Crews & Turner 2026, abstract](https://arxiv.org/pdf/2608.00291). [contested: an earlier sub-Alfvénic kink-suppression result was attributed, in a later Comment cited there, to an axis boundary-condition artefact]
- **Magnetic Rayleigh–Taylor [mainstream]:** in Gaussian form λ_c = B₀²cos²θ/((ρ₂−ρ₁)g). Published runs use B₀ = 10 G, ρ_cor = 3.7×10⁻¹² and ρ_chrom/prom = 2.9×10⁻¹⁰ kg m⁻³. At θ = 90° there is no cut-off; **θ = 89° → λ_c ≈ 38 km; θ = 88° → λ_c ≈ 155 km** — [Khomenko et al. 2014, arXiv:1310.7016](https://ar5iv.labs.arxiv.org/html/1310.7016). In SI the same result is λ_c = 4πB²cos²θ/(μ₀gΔρ), from the dispersion relation γ² = Agk − 2k²B²cos²θ/(μ₀(ρ_h+ρ_l)) with Atwood number A (computed). The peak growth rate sits at λ_max = 2λ_c — [Stone & Gardiner 2007, arXiv:0709.0452](https://ar5iv.arxiv.org/html/0709.0452) (search-result summary). With an inclined field the cut-off wavenumber disappears — [Vickers, Ballai & Erdélyi 2020, arXiv:1911.12734](https://arxiv.org/pdf/1911.12734) (search-result summary).
- **Tearing mode, constant-ψ FKR regime [mainstream]:** γ = [Γ(1/4)/(2πΓ(3/4))]^{4/5} Δ′^{4/5}/(τ_H^{2/5} τ_R^{3/5}). Normalised, γ̂ = 0.547 Δ′^{4/5} S^{−3/5}, unstable whenever Δ′ > 0. Layer width δ̂ ~ (γ̂/S)^{1/4}, so δ̂ ~ S^{−2/5} at Δ′ ~ 1. Valid for Δ′ ≪ S^{1/3} and γ̂ ≪ S^{−1/3} — [Fitzpatrick, Tearing Modes, Eqs. 9.59–9.66](https://farside.ph.utexas.edu/teaching/plasma/lectures/node127.html). The η^{3/5}Δ′^{4/5} scaling is also in [UW–Madison MHD lecture 30](https://magnetohydrodynamics.physics.wisc.edu/lecture30.html) (search-result summary). Another normalisation (oblique modes) uses [Γ(1/4)/(πΓ(3/4))]^{4/5} ≈ 0.95, a factor 2^{4/5} different ([arXiv:1201.0313](https://arxiv.org/pdf/1201.0313), search-result summary). Do not mix the two.
- **Computed:** coefficient = **0.5474**. At S = 10⁶ and Δ′a = 1, γτ = 1.375×10⁻⁴, layer width 4.0×10⁻³ a, and constant-ψ holds while Δ′a ≪ 100. At S = 10⁸: γτ = 8.68×10⁻⁶.
- **Interchange in a pinch:** this is the m = 0 Kadomtsev criterion above. I retrieved no source for the toroidal Suydam or Mercier forms (see Gaps).
- **Resistive wall modes:** none of the retrieved sources gives a closed-form growth rate. The only quantitative statement found is the wall-plus-feedback β_N range of 4–6 above the no-wall Troyon value ([arXiv:2008.10485](https://arxiv.org/pdf/2008.10485), search-result summary).

### Inferences
- The Bennett results are an analytic chain that a Rust test can run end to end. Force balance on p = (1+x²)⁻² gives B_θ²/2μ₀ = x²/(1+x²)², hence β = 1/x² (I derived this by hand and the script checks it). Then −d ln p/d ln r = 4x²/(1+x²), and the m = 1 criterion 4x²/(1+x²) ≤ x² becomes x ≥ √3. That exactly matches the published β* = 1/3, r* = √3 r_p, and it plugs into `plasmoid::pinch::bennett_*`.
- The m = 0 result for Bennett, 4x²/(1+x²) ≤ 4γx²/(2x² + γ), reduces to γ ≥ 2. A Bennett pinch is therefore sausage-stable only for γ ≥ 2, consistent with the published "α ≤ γ" rule since α = 2 for Bennett.
- The RT test needs g = 274 m s⁻², the standard solar surface gravity. Khomenko et al. do not state g, so it is my input. The 38.8 km and 155.3 km results agree with the published "≈ 38" and "≈ 155" km.
- Tearing growth needs Δ′ from an outer ideal solution, which depends on the current profile. A Rust crate can expose `fkr_growth(delta_prime_a, s, tau_a)` and take S from `plasmoid::relaxation::lundquist`. The plasmoid instability (already in the crate) is the S^{1/3}-limit continuation, because constant-ψ fails once Δ′ ≳ S^{1/3}.

### Gaps
- I found no published numeric kink or sausage growth rate (γa/v_A against ka) for a sharp or Bennett Z-pinch. Crews & Turner plot them (Figs. 1–4), but numbers were not extractable.
- Suydam (cylindrical interchange, D_S < 1/4) and Mercier criteria, and the resistive interchange, were not sourced this session.
- Resistive wall mode: no closed-form γτ_w formula or tested numbers retrieved.
- The non-constant-ψ (m = 1 resistive kink, γ ∝ S^{−1/3}) formula sits on Fitzpatrick's next page, which I did not fetch.

## Q5. Analytic MHD flows: Hartmann flow and its pressure drop; duct and pipe flow in liquid-metal blankets

### Takeaway
Hartmann flow between plates at y = ±a, with B normal to them, has the exact profile u/ū = Ha(cosh Ha − cosh(Ha·y/a))/(Ha cosh Ha − sinh Ha), with flat cores and Hartmann layers of thickness a/Ha. The pressure gradient is −dp/dx = σūB²·tanh Ha/(Ha − tanh Ha) for insulating walls (≈ σūB²/Ha) and σūB²·Ha/(Ha − tanh Ha) for perfectly conducting walls (≈ σūB²). I derived both and the script checks them, including the Poiseuille limit. For blanket ducts with conducting walls, the inertialess high-Ha result −dp/dx = σūB²·c/(1+c) reproduces KIT's published non-dimensional gradient: 0.0678 at c = 0.0727. The total pressure drop scales as ūB².

### Cited Findings
- Hartmann number Ha = BL√(σ/(ρν)) and Re = u₀L/ν, with L = half the duct size. KIT's 3D DCLL blanket simulations at **500 ≤ Ha ≤ 8000** show "the total pressure drop scales as p ~ u0B²", "the typical scaling for inertialess MHD flows in strong magnetic fields". Bends and U-turns do not add much to the total — [Bühler, Brinkmann, Klüber, Köhly, Mistrangelo, KIT MHD group report](https://publikationen.bibliothek.kit.edu/1000097959/38980113). [mainstream]
- **Circular-pipe test value** (MEKKA facility, NaK, 50 °C): wall conductance parameter c = (σ_w/σ)(R₀² − R²)/(R₀² + R²) = **0.0727** for the thick-walled pipe, and c_FCI = 0.00476 for the flow-channel-insert steel sheet. The fully developed gradients, normalised by σu₀B² (the report says "normalized by u0B²L"), are **0.0678** (pipe) and **0.00508** (infinitely long FCI), a predicted reduction of **13.3**. One uninsulated gap adds p₃D = 0.12, equal to the drop over l₃D = 23.6 of fully developed FCI flow — [same KIT report](https://publikationen.bibliothek.kit.edu/1000097959/38980113). [mainstream]
- The induced currents close through the walls, and the resulting Lorentz force is most of the pressure drop. Thin, resistive walls or insulating coatings cut it, and even a finite-resistance coating reduces it "by orders of magnitude". Hunt (1965) and Walker (1981) are the standard rectangular-duct references, and Chang & Lundgren (1961) treat thin walls of any cross-section — [KIT/UKAEA blanket MHD literature, search-result summaries](https://publikationen.bibliothek.kit.edu/270043415/3834101). A 2024 parametric study covers 0.01 ≤ c ≤ 0.5 and finds 3D effects grow as c falls — [UKAEA, Fusion Eng. Des. 2024](https://scientific-publications.ukaea.uk/wp-content/uploads/1-s2.0-S0920379624005386-main-1.pdf) (search-result summary).
- **Code-verification benchmarks:** the Smolentsev et al. 2015 V&V suite starts with 2D fully developed duct flow, i.e. Shercliff (insulating Hartmann walls) and Hunt (conducting Hartmann walls) flow. COMSOL verification reached Ha = 15 000 with flow-rate errors of ~0.2–0.4% (0.3596% at Ha = 15 000) — [Yan, Smolentsev, Abdou, COMSOL 2017](https://www.comsol.com/paper/validation-of-comsol-multiphysics-for-magnetohydrodynamics-mhd-flows-in-fusion-applications-50552) (search-result summary); also [FreeMHD, arXiv:2409.08950](https://arxiv.org/pdf/2409.08950).
- **Computed (Hartmann plates):**

| Ha | u_centre/ū | −dp/dx/(σūB²), insulating | −dp/dx/(σūB²), conducting |
|---|---|---|---|
| 0.01 | 1.5000 (Poiseuille) | 3.0×10⁴ (= 3/Ha²) | 3.0×10⁴ |
| 10 | 1.1110 | 0.1111 | 1.111 |
| 100 | 1.0101 | 0.0101 (≈ 1/(Ha−1)) | 1.010 |

### Inferences
- A Rust `flow::hartmann_profile(ha, eta)` and `flow::hartmann_gradient(ha, walls)` are exact, and four identities test them: mean = 1, the Poiseuille limit, the 1/(Ha−1) insulating asymptote, and the → 1 conducting asymptote. `flow::pipe_gradient(c) = c/(1+c)` is tested by KIT's 0.0678.
- The FCI value 0.00508 is not c/(1+c) = 0.00474. The likely cause is that the fluid inside the FCI has a smaller radius than the pipe-based scale L, so its mean velocity is higher. A ratio (L/R_FCI)² ≈ 1.072 would fit, but that is my inference; the report gives R_FCI only in a figure.
- These are inertialess, fully developed, uniform-B results. They fail at bends, manifolds, field gradients (3D current loops add the p₃D-type drops) and FCI gaps, and wherever buoyancy or MHD turbulence sets in.

### Gaps
- Shercliff and Hunt rectangular-duct series solutions, and their benchmark gradients at Ha = 500, 5000 and 15 000, were not retrieved as numbers (only the COMSOL error percentages).
- I retrieved no PbLi or NaK property data (σ, ρ, ν), so no dimensional blanket pressure-drop example in MPa is given.
- The Walker thin-wall rectangular formula c/(1 + a/(3b) + c) was not sourced and is left out.
- I did not source the Hartmann-layer laminar–turbulent transition criterion (commonly quoted as Re/Ha of a few hundred).

## Q6. Dynamo theory: kinematic thresholds, laboratory dynamos, αΩ dynamo number, geodynamo numbers, solar cycle, Christensen energy-flux scaling

### Takeaway
Kinematic dynamos need Rm above a geometry-dependent threshold. Analytic and laboratory values are: Roberts flow 3.5–3.9 (based on k), Ponomarenko 17.7, optimised Riga profiles 14.3 (convective) and 17.6 (absolute), VKS2 31–35 against numerical predictions of 48–133, and small-scale dynamos ~35 at Pm = 1. The Karlsruhe threshold is set by the α-number C = μ₀σα⊥R, measured at 8.4–9.3 against 8.12 predicted. For the plane-wave αΩ dynamo, D = αS/(η_T²k³) must exceed **2**, with cycle frequency ω = η_T k² at onset. The script recovers both values to 10⁻³ from the full dispersion relation. The geodynamo has Rm ~ 10³ and Λ ~ 0.3–5. Christensen's energy-flux law is B² = 1.2μ₀ f_ohm ρ^{1/3}(Fq)^{2/3}, fitted to dipolar simulations. It fits Earth, Jupiter and rapidly rotating low-mass stars, but I could not retrieve the planet-by-planet inputs needed for an Earth test.

### Cited Findings

**Laboratory and analytic kinematic thresholds [mainstream; experimentally confirmed]**

| Dynamo | Threshold (published) | Published facts | Source |
|---|---|---|---|
| Ponomarenko / Riga 1D theory | Rm_c = **17.7** (convective instability; Rm defined on the cylinder radius R). A back-flow makes it absolute | Gailitis & Freibergs 1976, 1980 | [Stefani et al. 2008](https://arxiv.org/abs/0807.0299); [Brandenburg & Subramanian 2005](https://arxiv.org/abs/astro-ph/0405052) |
| Riga (2000) | Bessel-function profiles: **12.0** (convective), **14.7** (absolute). Measured, optimised profiles: **14.3** and **17.6** | ~200 kW motor limit. First self-excitation in November 1999 at a sodium temperature of 210 °C; eigenfield at **1.326 Hz** at 2150 rpm. In July 2000 onset came above **≈1840 rpm**, followed by saturation. Hydraulic power ~200 kW with ~20 kW excess in saturation | [Stefani et al. 2008](https://arxiv.org/abs/0807.0299); [Gailitis et al. review, arXiv:1801.01749](https://ar5iv.labs.arxiv.org/html/1801.01749) |
| Karlsruhe (Dec 1999) | C = μ₀σα⊥R: measured C_crit = **8.4–9.3** against **8.12** from mean-field theory | 52 spin generators, radius 0.85 m, height 0.7 m. Self-excitation at V̇_C = 115 m³/h with V̇_H raised from 95 to 107 m³/h. Saturated at **≈7 mT** after ~120 s | [Stefani et al. 2008](https://arxiv.org/abs/0807.0299) |
| VKS2 (2006, soft-iron impellers) | Rm_c ≈ **31–35** against numerical predictions of 48–133. The field came out axisymmetric, "apparently at odds with Cowling's theorem" | 150 l sodium, 300 kW. Mean field ~40 G at ~30% above threshold. Reversals seen under asymmetric forcing | [Stefani et al. 2008](https://arxiv.org/abs/0807.0299); [Monchaux et al. 2007 PRL 98, 044502, arXiv:physics/0701075](https://arxiv.org/abs/physics/0701075) (threshold "about 30", 40 G, from the search-result summary) |
| Roberts flow (periodic) | R_crit = U₀/(ηk_f) ≈ **3.90** at k_z = k_x = k_y ≈ 0.71k_f. Based on k_z, the minimum is ≈ **3.49** at k_z ≈ 0.34k_f | — | [Brandenburg & Subramanian 2005, §4](https://arxiv.org/abs/astro-ph/0405052) |
| Small-scale (fluctuation) dynamo | Rm_crit ≈ **35** at Pm = 1. Kazantsev model gives 30–60; Kolmogorov WKBJ estimate ~60 | Helical large-scale dynamos need only ~1.2 | [Brandenburg & Subramanian 2005, §5](https://arxiv.org/abs/astro-ph/0405052) |
| Planetary rule of thumb | Rm = UR_c/λ > 50 | — | [Christensen 2009 lecture](https://cpaess.ucar.edu/sites/default/files/heliophysics/resources/presentations/2009_Christensen_planmagdynamo.pdf) |

**Mean-field dynamo numbers [mainstream]**
- α² dynamo: λ± = −η_T k² ± |αk|, unstable for 0 < k < α/η_T ≡ k_crit — [Brandenburg & Subramanian 2005, Eqs. 6.29–6.34](https://arxiv.org/abs/astro-ph/0405052).
- α²Ω dispersion: λ± = −η_T k² ± (α²k² − iαSk_z)^{1/2}. In the αΩ approximation, Re λ ≈ −η_T k² ± |½αSk_z|^{1/2} and ω_cyc ≈ |½αSk_z|^{1/2}. The dynamo-wave phase speed at marginality is c = ±η_T k_z, and the poloidal–toroidal phase shift is ±π/4 or ±3π/4 depending on the sign of the shear — [same, Eqs. 6.38–6.42, Table 6.2](https://arxiv.org/abs/astro-ph/0405052). **Computed:** marginality gives **D = αS/(η_T²k³) = 2**, with ω_cyc = η_T k² at onset (bisection on Eq. 6.38 with αk/S = 10⁻⁴ gives D_c = 2.0000 and ω = 1.0000 η_T k²).
- Spherical shells use two numbers, C_α = α₀R/η_T and C_Ω = ΔΩR²/η_T. In galactic discs D = C_ΩC_α, with C_Ω ≈ −10 and C_α ≈ 1, and D_crit ~ **6–10** depending on the α(z) profile — [same, §6.5.4–6.6](https://arxiv.org/abs/astro-ph/0405052).
- Solar cycle: ω_cyc = 2π/P_cyc, with P_cyc ≈ **11 years** for the Sun. The poloidal field reverses every 11 years — [same, §2](https://arxiv.org/abs/astro-ph/0405052).

**Geodynamo numbers [mainstream]** — [Christensen 2009 lecture](https://cpaess.ucar.edu/sites/default/files/heliophysics/resources/presentations/2009_Christensen_planmagdynamo.pdf), slides 19–20 (Earth versus model values): Ra* ≈ 5000× critical (?), E ≈ 10⁻¹⁴ (models ≥ 10⁻⁶), Pr 0.1–1, Pm ≈ 10⁻⁶ (models 0.06–20), Re ≈ 10⁸, **Rm ≈ 10³** (models 40–3000), Ro ≈ 5×10⁻⁶, **Λ = σB²/(2ρΩ) ≈ 0.3–5** (models 0.03–100). The dipole free-decay time is the longest decay time of a stagnant conductor, r²/(π²η) — [Christensen & Tilgner 2004](https://www.maths.gla.ac.uk/~rs/res/B/DynamoScaling/2004ChristensenTilgner.pdf) (search-result summary).

**Christensen, Holzwarth & Reiners 2009 energy-flux scaling [mainstream; prefactor and rotation-independence debated]**
- Hypothesis: "The magnetic energy density depends on thermodynamically available energy flux… The field strength is independent of rotation rate, conductivity, viscosity". Form: B²/2μ₀ ~ f_ohm ρ^{1/3}(L q_c/H_T)^{2/3}, where q_c is the convected heat flux, H_T = c_p/(αg) the temperature scale height, L the radial length scale, and f_ohm = ohmic/total dissipation. Fit to dipole-dominated models: **E_m* = 0.63 f_o (Fq*)^{2/3}**, or in dimensional form **B² = 1.2 μ₀ f_o ρ^{1/3}(Fq)^{2/3}**, where F is a thermodynamic efficiency factor. The fit holds for both Pm > 5 and Pm < 0.2. The planetary comparison assumes a ratio of 4–15 between the total internal field and the dipole field at the core–mantle boundary — [Christensen 2009 lecture, slides 34–36](https://cpaess.ucar.edu/sites/default/files/heliophysics/resources/presentations/2009_Christensen_planmagdynamo.pdf). The original paper is [Nature 457, 167 (2009)](https://www.nature.com/articles/nature07626) (not read).
- The law is "consistent with the maximum magnetic fields found in rapidly rotating low-mass stars, and some planetary fields like Earth and Jupiter". Stellar fields saturate when rotation is fast, "perhaps associated to the rate where Ro ≈ 0.1". Saturated pre-main-sequence stars and M dwarfs have fields of a few kG. Rapidly rotating evolved brown dwarfs are predicted at several hundred G to kG — [Reiners 2012, Living Rev. Sol. Phys., arXiv:1203.0241, §7.1](https://arxiv.org/abs/1203.0241).
- [contested in part] Schrinner 2012 argues the law is not truly independent of rotation, because the rotational dependence is absorbed into f_ohm when it is set to 1 — [arXiv:1212.6910](https://web3.arxiv.org/pdf/1212.6910) (search-result summary).
- **Computed, illustrative inverse:** with f_ohm = 1 and ρ = 1.1×10⁴ kg/m³ (assumed), B_rms = 1 mT needs Fq = 5.2 mW m⁻² and 4 mT needs 330 mW m⁻².

### Inferences
- A Rust `dynamo` module can carry `rm`, `elsasser`, `alpha2_kcrit`, `alpha_omega_marginal_d() = 2` and `cycle_frequency(eta, k)` with exact tests. It can also hold a table of published laboratory thresholds as constants, with doc-links, for "is this flow a dynamo?" checks.
- The prefactors 0.63 (dimensionless) and 1.2 (dimensional) are consistent if E_m* is B²/2μ₀ in units of ρ^{1/3}(Fq)^{2/3}, since 1.2/2 = 0.60 ≈ 0.63. That is my inference: the lecture slide is rounded and I did not read the normalisation in the paper.
- Fq is very sensitive to B (Fq ∝ B³), so an "Earth test" depends mostly on the assumed internal field strength (1–4 mT, itself uncertain by the lecture's own "?"). Treat the Christensen law as an order-of-magnitude predictor, which is how the paper presents it (log–log fits across ~10 decades).
- VKS shows that laboratory thresholds depend on boundary materials (soft iron cut Rm_c from 48–133 to 31–35). Kinematic predictions are not reliable to better than a factor of ~2 unless the boundaries are modelled.

### Gaps
- Christensen 2009's actual Earth and Jupiter inputs (q, F, ρ, f_ohm) and predicted values were not retrieved. The Glasgow and IPGP hosts refused connections, Nature is paywalled, and the 2025 giant-planet review (arXiv:2512.09725) exceeded the fetch size limit. No planet-level test value is available.
- I have no primary-source value for the electrical conductivity of the core (the 10⁶ S/m used is an assumption). The 49 kyr dipole decay time is therefore illustrative.
- I retrieved no solar-dynamo parameters (α, ΔΩ, η_T at the tachocline) that would give a published solar dynamo number.
- The Riga Rm based on a measured flow rate and the sodium σ(T) used in the experiments were not in the extracted text.

## Q7. MHD turbulence and Hall MHD: Goldreich–Sridhar critical balance, ion inertial length, when Hall terms matter

### Takeaway
Goldreich–Sridhar (1995) critical balance, k∥v_A ~ k⊥u⊥, gives E(k⊥) ∝ k⊥^{−5/3} and k∥ ∝ k⊥^{2/3}L^{−1/3}. Eddies become more elongated along B at small scales. The theory is phenomenological, holds for balanced (zero cross-helicity) turbulence, and is mainstream. Ideal MHD fails at the larger of d_i = c/ω_pi = v_A/Ω_i and ρ_i. Chen et al. (2014) showed that the solar-wind spectral break follows d_i at β_i ≪ 1 and ρ_i at β_i ≫ 1. Using their Table 1 means, f_b = v_sw/(2π max(d_i, ρ_i)) reproduces the published break frequencies: 0.635 Hz predicted against 0.58 ± 0.19 Hz, and 0.223 Hz against 0.223 ± 0.068 Hz.

### Cited Findings
- **Critical balance [mainstream, phenomenological]:** Goldreich & Sridhar 1995, ApJ 438, 763, set the cascade time (k⊥v_l)⁻¹ equal to the Alfvén period (k∥v_A)⁻¹. This gives l∥ ~ v_A ε^{−1/3} l⊥^{2/3}, i.e. k∥ ∝ k⊥^{2/3}, and E(k⊥) ∝ k⊥^{−5/3} — [Cho & Vishniac 2000, arXiv:astro-ph/0003403](https://arxiv.org/pdf/astro-ph/0003403) and related reviews (search-result summaries). Caveats: the theory "applies only to turbulence driven by nonlinear interactions among oppositely directed Alfvén modes carrying equal energies", and critical balance "remains a physically reliable hypothesis not strictly derived from basic principles". In mode decomposition, fast modes are isotropic (Iroshnikov–Kraichnan-like), while Alfvén and slow modes are anisotropic — same search set, e.g. [arXiv:1404.3691](https://arxiv.org/pdf/1404.3691). A recent overview is [Schekochihin 2022, "MHD turbulence: a biased review", arXiv:2010.00699](https://arxiv.org/pdf/2010.00699) (exceeded the fetch size; not read).
- **Ion-scale break [mainstream]:** ρ_i = v_th⊥i/Ω_i and d_i = c/ω_pi = v_A/Ω_i. The spectrum is close to −5/3 below ~0.3 Hz in the spacecraft frame. The break is defined where the local slope lies halfway between −5/3 and −8/3. Conclusion: "the break occurs at d_i [for β⊥i ≪ 1] and for β⊥i≫1 at ρ_i". Table 1 means: low β, B = 13.3 nT, n_i = 3.9 cm⁻³, v_sw = 460 km/s, T_i = 1.64 eV, **f_b = 0.58 ± 0.19 Hz** (18 intervals). High β, B = 1.14 nT, n_i = 17.1 cm⁻³, v_sw = 389 km/s, T_i = 4.8 eV, **f_b = 0.223 ± 0.068 Hz** (12 intervals) — [Chen et al. 2014, GRL 41, 8081](https://pmc.ncbi.nlm.nih.gov/articles/PMC4459186/). A 2D hybrid follow-up confirms that the break follows "the largest of the two" scales — [arXiv:1610.05158](https://ar5iv.labs.arxiv.org/html/1610.05158) (search-result summary).
- **Inertial lengths, exact test values:** PlasmaPy `inertial_length(5 m⁻³, "He+")` ≈ **2.03×10⁸ m** and `("e-")` ≈ **2.38×10⁶ m** — [PlasmaPy docs](https://plasmapy.readthedocs.io/en/stable/api/plasmapy.formulary.lengths.inertial_length.html) (values from the search-result summary of that page). **Computed:** 2.030×10⁸ m and 2.377×10⁶ m. The proton practical form is d_i = 2.277×10⁷ cm/√n[cm⁻³] (computed).
- **Computed scales:** solar wind, low-β set: d_i = 115 km, ρ_i = 14 km → f_b = 0.635 Hz. High-β set: d_i = 55 km, ρ_i = 278 km → f_b = 0.223 Hz (Taylor hypothesis f = v_sw/(2πl), with v_th = √(2k_BT/m)). ITER core: d_i = 3.6 cm. The GS anisotropy k⊥/k∥ = (k⊥L)^{1/3} is 100 at k⊥L = 10⁶.

### Inferences
- **When Hall MHD matters:** whenever the gradient scale of interest is ≲ d_i (at β ≲ 1) or ≲ ρ_i (at β ≳ 1). Examples: solar-wind turbulence below ~100 km, reconnection layers (Sweet–Parker width L·S^{−1/2} from `plasmoid::relaxation`, compared against d_i), and dense-plasma-focus pinch columns. In a tokamak core d_i is a few cm, much smaller than a, so global MHD is fine but tearing layers and island widths can reach d_i. A Rust `hall_regime(l, n, t, b)` returning the ratio l/max(d_i, ρ_i) is a direct calculator.
- The Chen match depends on the v_th = √(2k_BT/m) convention and on using one in 2π. The authors define the break by slope and work with spectra, so the 0.223 Hz exact agreement is partly fortunate. The low-β case (0.635 against 0.58 ± 0.19) is the more representative check.
- The GS spectral index (−5/3) and the Boldyrev dynamic-alignment index (−3/2) differ by an amount comparable to observational scatter in the perpendicular spectrum. The anisotropy k∥ ∝ k⊥^{2/3} is the more distinctive GS prediction.

### Gaps
- Boldyrev 2006 (dynamic alignment, E ∝ k⊥^{−3/2}) and the "tearing-mediated" range below the alignment scale were not read this session (the Schekochihin review was too large to fetch). I list them as [contested] alternatives, without a source-backed number.
- I found no published numerical test of the GS anisotropy prefactor (only the scaling).
- I retrieved no Hall-reconnection rate (~0.1 v_A) source, so the link to the plasmoid crate's reconnection rates is not quantified here.

## Q8. Where each formula is used in practice and where it breaks

### Takeaway
Most of these formulas share two assumptions: single-fluid, collisional MHD at scales above d_i and ρ_i, and steady or linear behaviour. The test values above sit inside those regimes. Each formula fails at a predictable edge: kinetic scales (Q7), finite Rm or S (Q1, Q4), non-planar or kinetic shocks (Q2), shaping and profile effects (Q3), 3D and inertial effects (Q5), and boundary-material effects (Q6, VKS).

### Cited Findings

| Formula | Used for | Breaks when | Evidence |
|---|---|---|---|
| v_A, β, magnetosonic speeds (Q1) | Solar corona seismology, wave heating, solar wind, tokamak τ_A | Scale ≲ d_i or ρ_i (dispersive kinetic Alfvén and whistler branches); collisionless plasma with anisotropic pressure | Spectral break at max(d_i, ρ_i) — [Chen et al. 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC4459186/). Solar wind is "typically considered a collisionless plasma"; wave–particle effects ≥ collisions — [Wilson et al. 2018](https://arxiv.org/pdf/1802.08585) |
| Frozen-in flux, Rm (Q1) | Dynamos, reconnection onset, liquid-metal experiments | Rm ≲ Rm_c (flux diffuses). Turbulent flows show an effective conductivity reduction, ~1% at Rm ≈ 1 in experiments; reversal statistics suggest ~3× in the Earth's core | [Stefani et al. 2008](https://arxiv.org/abs/0807.0299) |
| Perpendicular RH jump (Q2) | Bow shock, CME and interplanetary shocks, termination shocks | Oblique θ_Bn (needs the cubic); M_f < 1; kinetic, non-stationary high-M_A shocks; measured density ratios above 4 | [MMS bow-shock database, arXiv:2203.04680](https://arxiv.org/pdf/2203.04680) |
| Solov'ev GS (Q3) | Analytic equilibria for code verification, scoping studies | Solov'ev profiles force a current jump at the edge (unphysical); real p(ψ) and F(ψ) need numerical GS | [Cerfon, Guazzotto & Freidberg, EPS 2010](https://lac913.epfl.ch/epsppd3/2010/pdf/P4.115.pdf) |
| q95, Greenwald, Troyon, IPB98 (Q3) | Tokamak design (ITER, SPARC operating points) | Empirical fits: IPB98 carries a ±14% (1σ) / ±28% (95%) uncertainty for ITER and is being refitted (18 machines). Troyon rises to β_N ~4–6 with a wall plus feedback. Greenwald is a fit, not a hard limit (exceeding it is not sourced here) | [Sips 2005](https://arxiv.org/abs/physics/0410263); [IAEA ITPA contribution](https://conferences.iaea.org/indico/event/151/contributions/6399/); [arXiv:2008.10485](https://arxiv.org/pdf/2008.10485) |
| Kadomtsev, kink criteria (Q4) | Z-pinch design (sheared-flow pinches), pinch profile tailoring | Ideal-MHD energy principle only. Finite Larmor radius (gyroviscosity) reportedly stabilises marginal modes (search-result summary; the exact source was not pinned down). Sheared flow stabilises the kink but drives new shear modes | [Crews & Turner 2026](https://arxiv.org/pdf/2608.00291) |
| Magnetic RT (Q4) | Prominences, MagLIF and Z-pinch liners, supernova remnants | Inclined fields remove the cut-off. Partial ionisation and two-fluid effects | [Vickers et al. 2020](https://arxiv.org/pdf/1911.12734); [Khomenko et al. 2014](https://ar5iv.labs.arxiv.org/html/1310.7016) |
| FKR tearing (Q4) | Tokamak NTM/tearing onset, current-sheet stability | Δ′ ≳ S^{1/3} (non-constant-ψ, then the plasmoid regime); two-fluid layers when layer width < d_i | [Fitzpatrick](https://farside.ph.utexas.edu/teaching/plasma/lectures/node127.html) |
| Hartmann / c/(1+c) (Q5) | Blanket duct sizing, MHD pump and flowmeter design | 3D effects at gaps, bends and field gradients; FCI leakage; buoyancy or turbulence | [KIT report](https://publikationen.bibliothek.kit.edu/1000097959/38980113) |
| Rm_c, D_crit (Q6) | Dynamo experiment design, stellar and galactic dynamo models | Kinematic only (no saturation). Thresholds depend on boundary conditions (VKS 31–35 against 48–133 predicted) | [Stefani et al. 2008](https://arxiv.org/abs/0807.0299) |
| Christensen law (Q6) | Exoplanet and brown-dwarf radio-emission predictions | Needs fast rotation (dipolar regime, local Rossby ≲ 0.1). f_ohm and F are uncertain. Rotation dependence is debated | [Reiners 2012](https://arxiv.org/abs/1203.0241); [Schrinner 2012](https://web3.arxiv.org/pdf/1212.6910) |
| GS critical balance (Q7) | Solar wind, ISM, cosmic-ray scattering models | Imbalanced turbulence; scales below max(d_i, ρ_i); the −3/2 alternative | [Cho & Vishniac 2000](https://arxiv.org/pdf/astro-ph/0003403); [Chen et al. 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC4459186/) |

### Inferences
- Proposed crate layout (MIT, SI, clean-room from the cited forms):
  - `mhd-core`: Q1 parameters, magnetosonic speeds, Hall-regime ratio.
  - `mhd-shock`: perpendicular RH.
  - `mhd-equilibrium`: Solov'ev / Cerfon–Freidberg solver, q95, Greenwald, Troyon, IPB98.
  - `mhd-stability`: Kadomtsev, Kruskal–Shafranov, RT λ_c, FKR.
  - `mhd-flow`: Hartmann, pipe c/(1+c).
  - `mhd-dynamo`: Rm thresholds, αΩ.

  Each test in `mhd_theory_checks.py` maps one-to-one to a Rust unit test. The `plasmoid` crate supplies Bennett, Lundquist and Sweet–Parker, so these modules should depend on it rather than duplicate it.
- Status summary. [mainstream, exact or derived]: Q1 forms, RH, Solov'ev, Kadomtsev, magnetic RT, FKR, Hartmann, α²/αΩ dispersion. [mainstream, empirical fits with stated scatter]: IPB98, Troyon, Greenwald, q95 shaping formula, Christensen law. [contested]: VKS axisymmetric dynamo mechanism (open), Christensen's rotation independence, GS versus alignment spectra, sheared-flow kink stabilisation claims. I found nothing in this scope that is [refuted].

### Gaps
- Relativistic Alfvén-speed correction (v_A/√(1 + v_A²/c²)), Hall-MHD dispersion, and two-fluid tearing were not sourced this session.
- I did not quantify the experimental Greenwald-exceedance or Troyon-exceedance records (e.g. DIII-D).
