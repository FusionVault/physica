# Light, lensing and cosmology gaps: methods, primary sources and test values

Scope: the five gaps in `gravitics::spacetime`, `gravitational::{cosmology, lensing}` and `physica::analogies`. For each gap the notes give the method, the primary source, test values and where the code goes. Every published number below was read in a source that was opened, as cited. Every recomputed number comes from a script in `llc_gaps/scripts/`, next to these notes, with its saved `.out` file. The downloads are in their own `llc_gaps/dl_*` subfolders. Units are GM/c² = 1 for the shadow, SI elsewhere, and km/s/Mpc for H₀.

Scripts: `kerr_shadow.py`, `kerr_scan.py`, `jupiter_deflection.py`, `cosmology_checks.py`, `hogg_identities.py`, `finite_source.py`, `lens_on_limb.py` and `nfw_and_analogies.py`. Each has a `.out` file. Run them with python3, numpy and scipy. astropy is not installed, so its test values are reproduced independently.

## Gap 1: the Kerr shadow at any inclination, and the size measures

### Takeaway
Bardeen's parametric curve is exact and cheap. Each spherical photon orbit radius r gives (ξ, η), then a point (α, β). The curve gives the full outline at any inclination, and the existing pole-on function is its i → 0 limit. The size measure has to be named. The mean, median and areal radii all fall at most **7.08 %** below √27, exactly 2(1+√2)/(3√3) − 1 at a = 1 seen pole-on. The horizontal width falls up to **13.4 %** (to 9 M at a = 1, i = 90°). The EHT's "up to ~7.5 %" is a rounded envelope of the mean-radius measure.

### Cited Findings
- **Bardeen's curve, as restated by Cunha & Herdeiro (Sect. "The Kerr shadow").** A spherical photon orbit at radius r has η_CH ≡ L/E = (r³ + a²r + Ma² − 3Mr²)/(a(M − r)) and χ ≡ Q/E² = r²/(r² − a²)·(3r² + a² − η_CH²) ≥ 0. A static observer at infinity at polar angle θ_o sees the edge at x = −η_CH/sin θ_o, y = ±√(χ + a²cos²θ_o − η_CH²/tan²θ_o). Bardeen 1973 is cited for this. — [Cunha & Herdeiro 2018](https://arxiv.org/abs/1801.00860)
- **Notation pitfall.** Cunha & Herdeiro's η is Bardeen's ξ = L_z/E, and their χ is Bardeen's η = Q/E². Johannsen uses ξ = L_z/E and η = C/E², with x′ = −ξ/sin i and y′ = ±√(η + a²cos²i − ξ²cot²i) "similar to the ones for the Kerr metric (Bardeen 1973)". The y′ axis lies in the plane of the spin axis and the image is centred on the black hole. — [Johannsen 2013](https://arxiv.org/abs/1501.02814); [Cunha & Herdeiro 2018](https://arxiv.org/abs/1801.00860)
- **The curve as an explicit y(x).** The root r of r³ − 3Mr² + a(a+η)r + Ma(a−η) = 0 has a closed Viète form on three branches, with 𝒜 = M² − a(η+a)/3 and ℬ = M(M² − a²)|𝒜|^(−3/2). The paper notes that "y(2a) = 3√3 M for θ_o = π/2". — [Cunha & Herdeiro 2018](https://arxiv.org/abs/1801.00860)
- **The JP2010 size definitions.**
  - Displacement: D ≡ |x′_max + x′_min|/2.
  - Average radius: ⟨R⟩ ≡ ∫R dα/∫dα, with R ≡ √((x′ − D)² + y′²) and **tan α = y′/x′, the angle measured about the origin (the black hole)**.
  - Diameter: L ≡ 2⟨R⟩.
  - Asymmetry: A ≡ 2√(∫(R − ⟨R⟩)² dα/∫dα).
  - Fits: D = 2a sin i(1 − 0.41ε sin²i), plus 0.12a² sin i for a ≳ 0.7M; A/M = [0.84ε + 0.36(a/M)³] sin^{3/2} i.
  - "L = 10.4M" for 0 ≤ a ≤ 0.4, with a "systematic uncertainty of only ≃2%" at a = 0.9M.
  — [Johannsen & Psaltis 2010](https://arxiv.org/abs/1005.1931)
- **Johannsen 2013's semi-analytic fits.** He evaluates the same definitions semi-analytically, L = (2/π)∫⟨R̄⟩H(r)dr with H = |d arctan(y′/x′)/dr|, and fits 17,113 rings with L_K = L₁ + L₂cos(n₁i + φ) and D_K = D₁ sin(n₂i). The coefficients:
  - L₁ = 10.3907 − 0.520134a^2.19644 − 0.109423a^11.5509
  - D₁ = 2a + 0.188324a^3.29029 + …
  - n₁, φ, L₂ and n₂ as printed in the paper.

  The quoted accuracy is "<0.09% for a ≤ 0.997M" for L, and "<1.1% for a ≤ 0.97M" for D. For 0 ≤ |a| ≲ 0.77M, "10M ≲ L ≲ 10.4M". — [Johannsen 2013](https://arxiv.org/abs/1501.02814)
- **Chan, Psaltis & Özel 2013** use the same D, ⟨R⟩ and A, from GRay ray tracing with the image plane at r = 1000M. Their fits:
  - ⟨R⟩ ≃ R₀ + R₁cos(2.14i − 22.9°), with R₀ = (5.2 − 0.209a + 0.445a² − 0.567a³)M and R₁ = [0.24 − 3.3/((a − 0.9017)² + 0.059)]×10⁻³M.
  - A ≃ A₀ sin^n i, with A₀ = (0.332a³ + 0.176a^21.7 + 0.0756a^195)M and n = 1.55(1 − a)^−0.022 + 1.3(1 − a)^0.98.
  - Arguments in degrees.
  - "the ring retains a highly circular shape … significant asymmetries appear only at a ≳ 0.99".
  — [Chan et al. 2013](https://arxiv.org/abs/1303.5057)
- **Takahashi 2004** defines the shift x_c = (x_max + x_min)/2. For an optically thin medium, "the maximum size of a black hole shadow … concentrates on a narrow range … from ∼9.6 to ∼10.4", with r_g = Gm/c². His disc-occulted shadows are a different quantity. — [Takahashi 2004](https://arxiv.org/abs/astro-ph/0405099)
- **EHT, M87 Paper VI.** The shadow diameter "is predicted to be ≈9.6–10.4 GM/c²" (Takahashi 2004; Johannsen & Psaltis 2010), a "mean diameter that lies in the very small range (9.6–10.4) GM/c²". Figure 21 shades "the maximum ±4% discrepancy". At M87's inclination "the maximum asymmetry is 2% for all values of black hole spin (e.g., Chan et al. 2013)". — [EHT 2019 Paper VI](https://arxiv.org/abs/1906.11243)
- **EHT, Sgr A\* Paper VI.**
  - δ ≡ d_sh/d_sh,Sch − 1, with d_sh,Sch = 6√3 θ_g.
  - "Black holes with non-zero spin observed at different inclinations can have shadow sizes that are smaller by up to ∼7.5% from this limit (Takahashi 2004; Chan 2013)", so δ ∈ [−0.075, 0] is consistent with Kerr.
  - The radius "ranges from ∼4.8GMc⁻² to ∼5.2GMc⁻²".
  - For Kerr "we use the analytic approximation derived in Chan (2013)".
  - In the metric tests, d̃_metric "is the median shadow diameter, i.e., the locus of critical impact parameters", and size is defined "as its median radius".
  — [EHT 2022 Sgr A\* Paper VI](https://arxiv.org/abs/2311.09484)

### Inferences
- **Method (recommended).** In M = 1 units, use Bardeen's form, which is equivalent to Cunha & Herdeiro's. I checked numerically that their χ equals Bardeen's η on the photon orbits.

  ξ(r) = [r²(3 − r) − a²(r + 1)]/[a(r − 1)],  η(r) = r³[4a² − r(r − 3)²]/[a²(r − 1)²],
  α = −ξ/sin i,  β = ±√(η + a²cos²i − ξ²cot²i).

  - r runs over [r₁, r₂], the roots of β² = 0. Bracket them in [r_pro, r₀] and [r₀, r_retro], where r₀ is the existing ξ = 0 root of r³ − 3r² + a²r + a² = 0 (β² > 0 there for every i). The endpoints are the equatorial photon orbits r_pro/retro = 2{1 + cos[(2/3)arccos(∓a)]}.
  - At i = 90° exactly, r₁ = r_pro and r₂ = r_retro.
  - a = 0 needs its own branch (ξ has a in the denominator). It is the circle √27.
  - Cluster the samples at both ends with r = r₁ + (r₂ − r₁)(1 − cos t)/2, because β has square-root endpoints.
  - For the averages, the trapezoid rule over the parametrised curve converges to 1×10⁻¹¹ at 10⁵ points (`kerr_shadow.out`). Gauss–Legendre in t would need far fewer.
- **Recomputed test values** (`kerr_shadow.out`, `kerr_scan.out`), in units of GM/c²:

  | a | i | width α_max−α_min | height 2β_max | D | ⟨R⟩ (JP, angle about origin) | L = 2⟨R⟩ | Johannsen fit L | Chan ⟨R⟩ | median R | areal R | A |
  |---|---|---|---|---|---|---|---|---|---|---|---|
  | 0.5 | 90° | 10.234422 | 10.392305 | 1.020945 | 5.157025 | 10.314050 | 10.313523 | 5.150373 | 5.157138 | 5.157355 | 0.056866 |
  | 0.9 | 17° | 9.806789 | 9.890003 | 0.600227 | 4.924470 | 9.848941 | 9.848923 | 4.904852 | 4.924465 | 4.924498 | 0.029827 |
  | **0.9** | **90°** | **9.676741** | **10.392305** | **1.993949** | **5.027489** | **10.054979** | 10.056937 | 5.013799 | 5.028450 | 5.033226 | 0.278187 |
  | 0.998 | 90° | 9.107554 | 10.392305 | 2.442889 | 4.925690 | 9.851379 | 9.862710 | 4.918355 | 4.923983 | 4.938024 | 0.553937 |

  At a = 0.9 and i = 90°:
  - The exact mean diameter is 10.054979, against Johannsen's fit 10.056937, which is +0.019 % and inside his stated 0.09 %. This gives δ = −0.0325.
  - The displacement D = 1.993949, against Johannsen's fit 1.982213 (−0.59 %, inside his 1.1 %) and JP2010's 2a sin i + 0.12a² sin i = 1.8972 (−4.9 %).
  - Chan's ⟨R⟩ is 0.27 % low.
  - The asymmetry A = 0.2782, against JP2010's 0.36a³ = 0.262 and Chan's A₀ = 0.260.
  - Measuring the angle about the centre (D, 0) instead of the origin gives ⟨R⟩ = 5.031316 (+0.08 %). This is a convention trap.
- **Exact identities for doc tests**, all reproduced to round-off (`kerr_shadow.out`):
  1. **Pole-on.** As i → 0 the outline becomes a circle of radius √(η(r₀) + a²) = √(8r₀³/(r₀² − 1)), the existing `kerr_pole_on_shadow`. ⟨R⟩ at i = 10⁻⁴ matches it to 1×10⁻⁹ for a = 0.3, 0.9 and 0.999: 5.16974818, 4.91606338 and 4.82942567.
  2. **Edge-on, every spin.** r = 3 gives ξ = −2a and η = 27 exactly. The edge passes through (α, β) = (2a, ±3√3), the paper's "y(2a) = 3√3 M", and r = 3 is a stationary point of η. So the vertical diameter at i = 90° is **6√3 for every a**: the table's height column is 10.392305 at every spin.
  3. **Edge-on ends.** α_min = −ξ(r_pro) and α_max = −ξ(r_retro). At a = 0.9: −2.844421403 and 6.832319230. As a → 1 they tend to −2 and 7: width 9, D = 5/2.
  4. **Range over all a and i.** Mean, median and areal radius all run from **2(1+√2) = 4.828427** (a = 1, i → 0, where r₀ = 1+√2) to √27 = 5.196152. So δ_min = 2(1+√2)/(3√3) − 1 = **−0.070769**, while the half-width reaches −0.13374.
- **The 7.5 % claim.** It is defensible only as a rounded envelope of the mean-type diameter: Takahashi's "∼9.6" lower end against 6√3 = 10.39 gives −7.6 %. The exact Kerr floor for the mean, median and areal measures is −7.08 %. A shadow "diameter along an axis" can shrink 13.4 %. The doc comment on `shadow_diameter` should name the measure and give −0.0708 as exact. The EHT's −0.075 can stay as their published bound.
- **Where it goes.** In `gravitics::spacetime`, next to `kerr_pole_on_shadow`:
  - `kerr_shadow_edge(a, inclination, r) -> Option<[f64; 2]>`, returning α and β, and `None` where β² < 0.
  - `kerr_shadow_size(a, inclination) -> (mean_radius, displacement)` for the JP2010/Johannsen definitions.

  Doc tests: identities 1–4, plus the a = 0.9, i = 90° row checked against Johannsen's fit to 0.09 %. `shadow_deviation` already uses the EHT δ.

### Gaps
- Bardeen 1973 (Les Houches) was not opened; its formulas were taken from Cunha & Herdeiro and Johannsen 2013.
- The EHT's "median radius" algorithm (Medeiros et al. 2020; Psaltis et al. 2020) was not opened. My median, taken about (D, 0) over uniform angle, is one reasonable reading: it differs from ⟨R⟩ by 2×10⁻⁴ at a = 0.9, i = 90°.
- No published high-precision table of the exact mean diameter was found. Johannsen's fit, accurate to 0.09 %, is the most precise published comparison.

## Gap 2: multipole deflection at any orientation, and Jupiter's J₂ from primary sources

### Takeaway
Crosta & Mignard give the general-orientation J₂ deflection in closed form. The library's Born integral reproduces it to 2×10⁻¹². Its modulus is (1+γ)(2GMJ₂R²/c²b³)sin²ϑ with a spin-2 orientation pattern, so the deflection vanishes for light travelling along the spin axis. For the full zonal field on an equatorial ray there is an exact series: (4GM/c²b)[1 + J₂(R/b)² − J₄(R/b)⁴ + J₆(R/b)⁶ − …]. The odd J's give a north–south component of the same alternating form. On Jupiter's J₂: **14696.572 is Iess et al. 2018** (confirmed in the paper's own table). **14696.5735 (unnormalised, permanent tide included) is Durante 2020** (PDS SHADR). **14696.5063 is exactly Durante's value minus its 6.72×10⁻⁸ permanent tide.** **14696.514 ± 0.272 is attributed to Bolton et al. 2017 (Science), not Folkner 2017**, in the secondary paper that quotes it.

### Cited Findings
- **The deflection integral.** ΔΦ = (1+γ)∫∇⊥U dl. With χ the star–observer–planet angle, the general result is ΔΦ = ΔΦ₁n + ΔΦ₂m, with the χ-dependent factors given in the paper's Eqs. "defqu1-chi" and "defqu2-chi". For χ ≪ 1 (the infinite line):

  ΔΦ₁ = (2(1+γ)M/b)[1 + J₂(R²/b²)(1 − 2(n·z)² − (t·z)²)]
  ΔΦ₂ = (4(1+γ)MJ₂R²/b³)(m·z)(n·z)

  Here M = GM/c², t is the ray direction, n is the radial unit vector and m = t × n is the orthoradial one. — [Crosta & Mignard 2006](https://arxiv.org/abs/astro-ph/0512359)
- **Appendix B** states that "the modulo of the quadrupole deflection … will depend only on sin²ϑ, … the angle between the z axis of Jupiter and the direction from the star to the observer", with components ∝ cos 2α and sin 2α. — [Crosta & Mignard 2006](https://arxiv.org/abs/astro-ph/0512359)
- **Three typos in Crosta & Mignard.** Eq. "B-quad-mod" prints R³/b³, which is dimensionally wrong; it should be R². Eq. "B-mod1" prints −2(t·z)², but the main formula has −(t·z)². Eq. "defqu2-gam-eps" carries a stray extra γ factor. Implement from Eqs. "defqu1" and "defqu2". — [Crosta & Mignard 2006](https://arxiv.org/abs/astro-ph/0512359)
- **Crosta & Mignard Table 1**, grazing rays in Gaia geometry, post-Newtonian / quadrupole: Sun 1.″75 / ∼1 µas; Earth 574 / 0.6 µas; Jupiter 16280 / 239 µas; Saturn 5772 / 94 µas; Uranus 2081 / 25 µas; Neptune 2535 / 9 µas. Table 2 simulation parameters: R_J = 7.13×10⁷ m, GM_J/c² = 1.41 m, J₂ = 0.014736. "Using J₄ ∼ −0.0006 … this results for a grazing ray into a deflection of 10 µas". — [Crosta & Mignard 2006](https://arxiv.org/abs/astro-ph/0512359)
- **The Durante2020 SHADR file** (PDS: `jnogrv_0002/DATA/SHADR/GRV_JUGR_SHA_JUP_U_DUR2020_20200203_V01.TAB` and `.LBL`).
  - Header: reference radius 71492 km; GM = 126686534.1 ± 8.4 km³/s²; degree 12, order 2; normalisation-state column = 0, which the label defines as "0 coefficients are unnormalized".
  - C₂₀ = −1.46965735×10⁻² ± 1.7×10⁻⁹. Uncertainties are "3x formal (3-sigma)".
  - Zonal C_n0: C₃₀ = 4.5e-8, C₄₀ = 5.866085e-4, C₅₀ = 7.23e-8, C₆₀ = −3.42007e-5, C₇₀ = −1.2e-7, C₈₀ = 2.422e-6, C₉₀ = 1.13e-7, C₁₀,₀ = −1.81e-7, C₁₁,₀ = −1.6e-8, C₁₂,₀ = −6.2e-8.
  - The label: "The spherical harmonic coefficients are UNNORMALIZED … The value for J2 includes a permanent tidal contribution effect estimated to be 6.72e-8 from interior model prediction". The pole and Love numbers are also given.
  — [PDS TAB](https://atmos.nmsu.edu/PDS/data/jnogrv_0002/DATA/SHADR/GRV_JUGR_SHA_JUP_U_DUR2020_20200203_V01.TAB); [PDS LBL](https://atmos.nmsu.edu/PDS/data/jnogrv_0002/DATA/SHADR/GRV_JUGR_SHA_JUP_U_DUR2020_20200203_V01.LBL); [PDS model note](https://pdsatmospheres.nmsu.edu/data_and_services/atmospheres_data/JUNO/logs/The%20Durante2020%20Gravity%20Field%20Model.pdf)
- **The SHADR file-name convention** is GRV_JUGR_SHA_bbb_z_nnnn_YYYYMMDD_Vvv, where z = 'U' means UN-NORMALIZED and 'N' means NORMALIZED. The label calls the degree index "m" and the order "n", the reverse of the usual. — [Juno Gravity RDR SIS](https://atmos.nmsu.edu/data_and_services/atmospheres_data/JUNO/logs/JUNO_GRAV_RDR_SIS.PDF)
- **Iess et al. 2018's own table**, file nature25776-st2 on CaltechAUTHORS. Columns are Value / Uncertainty / Consider uncertainty, ×10⁻⁶:
  - J₂ = 14696.572 / 0.014 / 0.024
  - J₃ = −0.042 / 0.010
  - J₄ = −586.609 / 0.004
  - J₅ = −0.069
  - J₆ = 34.198 / 0.009
  - J₇ = 0.124
  - J₈ = −2.426
  - J₉ = −0.106
  - J₁₀ = 0.172
  - k₂₂ = 0.625
  - pole α = 268.0570°, δ = 64.4973°

  The caption, sigma level and reference radius are not in the image. — [Iess et al. 2018, CaltechAUTHORS record](https://authors.library.caltech.edu/records/x4fzt-7v779)
- **Iorio's quote** of J₂ = 14696.572×10⁻⁶ cites 2018Natur.555..220I, next to μ = 1.26713×10¹⁷ m³ s⁻² from IERS 2010. — [Iorio 2023](https://arxiv.org/abs/2310.02834)
- **Debras & Chabrier's quote** of J₂×10⁶ = 14696.514 ± 0.272 carries footnote (e), which is **Bolton et al. 2017, Science 356, 821**. Their GM_J = 126 686 533 ± 2 ×10⁹ m³ s⁻² carries footnote (b), Folkner et al. 2017. — [Debras & Chabrier 2017](https://arxiv.org/abs/1710.03969)

### Inferences
- **The Born integral matches the closed form** (`jupiter_deflection.out`, DUR2020 values, b = R unless stated):
  - Monopole at the limb, 4GM/(c²R) = 7.886632×10⁻⁸ rad = **16267.3 µas**. Table 1 gives 16280; Crosta & Mignard's own Table 2 parameters give 16316.
  - J₂ part at the limb with z ⟂ ray: **239.07 µas** (Table 1: 239).
  - Equatorial ray: |α| = 16506.4211 µas, the quadrupole adding +239.07426.
  - Ray over the pole (n ∥ z): 16028.2726 µas, the quadrupole subtracting 239.07426.
  - Ray along the spin axis (t ∥ z): exactly the monopole, 16267.3469 µas, with zero quadrupole.
  - Tilted ray: the closest point at 45° latitude, travelling at 30° to the equator, gives b = 0.9354R and quadrupole 219.06954 µas. That equals 4GMJ₂R²sin²ϑ/(c²b³) with sin²ϑ = 0.75.
  - In every case |Born − Crosta|/|quad| ≤ 2.1×10⁻¹².
- **Exact equatorial-ray series** for any zonal field, derived here and checked numerically. Use P_{2k}(0) and ∫(b² + l²)^{−(n+1)/2}dl = b^{−n}√π Γ(n/2)/Γ((n+1)/2). Then
  α_in-plane = (4GM/c²b)[1 + Σ_k (−1)^{k+1} J_{2k}(R/b)^{2k}],  α_z = (4GM/c²b) Σ_k (−1)^{k+1} J_{2k+1}(R/b)^{2k+1}.
  - At b = R: in-plane 16516.561365 µas from both the Born integral and the series. The terms are J₂ 239.0743, −J₄ 9.5426 (Crosta's "10 µas"), +J₆ 0.5564 and −J₈ 0.0394. α_z = +0.004495 µas.
  - At b = 2R: 8163.860339 µas, also from both methods.

  This is a high-precision doc test for `multipole_deflection` on the full DUR2020 zonal field.
- **Normalisation.** `geopotential::Coefficients` stores normalised C̄. Convert with C̄_n0 = C_n0/√(2n + 1); for Jupiter C̄₂₀ = −J₂/√5 = −6.5725074765×10⁻³. The existing test already divides by √5.
- **Resolving the three J₂ values.**
  - 14696.572 is Iess et al. 2018, primary and now confirmed, ± 0.014 (consider 0.024).
  - 14696.5735 ± 0.0017 (3σ) is Durante 2020. 14696.5735 − 0.0672 = **14696.5063** exactly, so the unconfirmed 14696.5063 is the Durante J₂ with the stated 6.72×10⁻⁸ permanent tide removed.
  - 14696.514 ± 0.272 is the early two-orbit value, credited to Bolton et al. 2017 in the secondary source.
  - Recommendation: replace the test's secondary citation with the DUR2020 PDS product: GM = 1.266865341×10¹⁷ m³/s², R = 7.1492×10⁷ m, unnormalised J_n, all 3σ. State that J₂ includes the permanent tide.
- **Second-order terms.** The second-order monopole term (15π/4)(GM/c²R)² is only 0.001 µas at Jupiter's limb, so the first-order Born deflection is complete at µas level.
- **Where it goes.** Tests only, in `gravitics::spacetime`:
  - Crosta's three orientation cases (polar, axial, tilted with the sin²ϑ modulus).
  - The DUR2020 zonal series at b = R and b = 2R.
  - Earth's 574 / 0.6 µas from the Earth coefficients already in `physica::analogies` (Gap 5).

  No new function is needed. A DUR2020 constant could live in the test, keeping the existing pattern.

### Gaps
- Klioner 1991 (Sov. Astron. 35, 523) was not opened.
- The Iess et al. 2018 main text and its Table 1 are behind the Nature paywall. The confirmation above comes from the paper's own table image on CaltechAUTHORS, which has no caption.
- Durante et al. 2020 (GRL) and Folkner et al. 2017 (GRL) could not be opened: Wiley and the Sapienza/Unibo repositories return a Cloudflare challenge, which was not bypassed. Folkner's own J₂ and its sigma level remain unverified.
- Saturn's GM, J₂ and R (Iess et al. 2019) were not opened, so Table 1's 94 µas was not recomputed.
- No newer Juno SHADR than DUR2020 is in PDS volume jnogrv_0002.

## Gap 3: curved, dark-energy and radiation cosmology distances

### Takeaway
Integrate D_C = D_H∫dz/E(z) once, with E² = Ω_m(1+z)³ + Ω_k(1+z)² + Ω_de·f_de(z) + Ω_γ(1+z)⁴[1 + ρ_ν/ρ_γ(z)]. Then apply S_k for curvature, and take D_A12 = S_k(χ₂ − χ₁)/(1 + z₂) for all Ω_k: **not Hogg's Eq. 19, which fails for Ω_k < 0**. Use the CPL form for w₀wₐ and the exact Fermi–Dirac integral for massive neutrinos. This reproduces every external-code value in astropy's tests to their printed digits (icosmo to 5×10⁻⁸, Mathematica w₀wₐ to 3×10⁻⁶, Mathematica neutrino densities to 3×10⁻⁶). The Planck18 preset recovers 100θ* = 1.04124 against Planck's 1.04119 ± 0.00029. Without radiation it gives 1.03639, which is 16σ off.

### Cited Findings
- **Hogg's distance formulas.**
  - D_H = c/H₀ (Eq. 4); Ω_M + Ω_Λ + Ω_k = 1 (Eq. 7).
  - E(z) = √(Ω_M(1+z)³ + Ω_k(1+z)² + Ω_Λ) (Eq. 14); D_C = D_H∫₀^z dz′/E(z′) (Eq. 15).
  - D_M = D_H Ω_k^{−1/2} sinh(√Ω_k D_C/D_H) for Ω_k > 0, D_C for Ω_k = 0, and D_H|Ω_k|^{−1/2} sin(√|Ω_k| D_C/D_H) for Ω_k < 0 (Eq. 16).
  - For Ω_Λ = 0 the closed form is D_M = D_H 2[2 − Ω_M(1 − z) − (2 − Ω_M)√(1 + Ω_M z)]/(Ω_M²(1 + z)) (Eq. 17); D_A = D_M/(1 + z) (Eq. 18).
  - D_A12 = [D_M2√(1 + Ω_k D_M1²/D_H²) − D_M1√(1 + Ω_k D_M2²/D_H²)]/(1 + z₂) (Eq. 19) holds "for Ω_k ≥ 0" and is "not correct for Ω_k < 0 (Phillip Helbig, 1998, private communication)".
  - Hogg has no radiation term.
  — [Hogg 1999](https://arxiv.org/abs/astro-ph/9905116)
- **External-code values in astropy's tests.**
  - test_flat_z1, FlatLambdaCDM(H0 = 70, Om0 = 0.27, Tcmb0 = 0), z = 1, from the Wright, Kempner and iCosmos calculators: D_C = 3364.5, 3364.8, 3364.7988 Mpc; D_A = 1682.3, 1682.4, 1682.3994; D_L = 6729.2, 6729.6, 6729.5976.
  - Kayser, Helbig & Schramm 1997 code: FlatLambdaCDM(100, 0.3) D_M(1, 2) = 1313.2232194828466 Mpc; FlatLambdaCDM(70.4, 0.272) D_A(1, 2) = 646.22968662822018 Mpc.
  - icosmo tables in data/cosmo_{flat,open,closed}.ecsv, all h = 0.7: (Ω_m, Ω_Λ) = (0.3, 0.7), (0.3, 0.1) and (2, 0.1), for z = 0 to 4.66. At z = 1, D_M = 3303.8288, 2979.346 and 2152.7954 Mpc.
  - Mathematica, w₀waCDM(H0 = 70, Om0 = 0.2, Ode0 = 0.8, w0 = −1.1, wa = 0.2): D_L(0.2, 0.4, 0.9, 1.2) = 1004.0, 2268.62, 6265.76, 9061.84 Mpc (rtol 10⁻⁴). de_density_scale(0.5, 1.5) = 0.9246310669529021, 0.9184087000251957. Also (0.3, 0.7, −0.9, 0): 971.667, 2141.67, 5685.96, 8107.41; and (0.3, 0.7, −0.9, −0.5): 974.087, 2157.08, 5783.92, 8274.08.
  - wCDM(70, 0.3, 0.6, w0 = −0.5) de_density_scale at z = 0.1, 0.2, 0.5, 1.5, 2.5: 1.15369, 1.31453, 1.83712, 3.95285, 6.5479.
  - Massive neutrinos, "compare … against the exact formula (eq 24/25 of Komatsu et al. 2011) computed using Mathematica". The ratio ρ_ν/ρ_γ/(0.2271·N_eff) at z = 0, 1, 2, 10, 1000, with Tcmb = 3 K and three equal masses:
    - 100 eV: 171969, 85984.5, 57323, 15633.5, 171.801
    - 0.25 eV: 429.924, 214.964, 143.312, 39.1005, 1.11086
    - 0.01 eV: 17.2347, 8.67345, 5.84348, 1.90671, 1.00021
    - mixed [0, 0.01, 0.25] eV with N_eff = 3.04: 149.386233, 74.87915, 50.0518, 14.002403, 1.03702333

    astropy's own fit "is only good to ∼0.5%".
  - Other non-flat z₁z₂ values in astropy's tests carry no stated source, for example LambdaCDM(70.4, 0.2, 0.5) D_A(1, 2) = 620.1175337852428 and LambdaCDM(100, 2, 1) D_A(1, 2) = 228.42914659246014.
  — [astropy test_lambdacdm.py](https://github.com/astropy/astropy/blob/main/astropy/cosmology/_src/tests/flrw/test_lambdacdm.py); [test_w0wacdm.py](https://github.com/astropy/astropy/blob/main/astropy/cosmology/_src/tests/flrw/test_w0wacdm.py); [test_w0cdm.py](https://github.com/astropy/astropy/blob/main/astropy/cosmology/_src/tests/flrw/test_w0cdm.py); [icosmo tables](https://github.com/astropy/astropy/tree/main/astropy/cosmology/_src/tests/flrw/data)
- **Planck 2018 Table 2**, the column "TT,TE,EE+lowE+lensing+BAO", 68 % limits:
  - Ω_bh² = 0.02242 ± 0.00014; Ω_ch² = 0.11933 ± 0.00091; 100θ_MC = 1.04101 ± 0.00029; τ = 0.0561 ± 0.0071; ln(10¹⁰A_s) = 3.047 ± 0.014; n_s = 0.9665 ± 0.0038.
  - H₀ = 67.66 ± 0.42; Ω_Λ = 0.6889 ± 0.0056; Ω_m = 0.3111 ± 0.0056; Ω_mh² = 0.14240 ± 0.00087; σ₈ = 0.8102 ± 0.0060.
  - z_re = 7.82 ± 0.71; age = 13.787 ± 0.020 Gyr; z* = 1089.80 ± 0.21; r* = 144.57 ± 0.22 Mpc; 100θ* = 1.04119 ± 0.00029; r_drag = 147.21 ± 0.23 Mpc.
  - The neighbouring column, without BAO: H₀ = 67.36 ± 0.54 and Ω_m = 0.3153 ± 0.0073. The abstract gives 67.4 ± 0.5.
  - "Ω_m includes the contribution from one neutrino with a mass of 0.06 eV".
  — [Planck 2018 VI](https://arxiv.org/abs/1807.06209)
- **astropy's Planck18 realization**: FlatLambdaCDM with H0 = 67.66, Om0 = 0.30966, Tcmb0 = 2.7255 K, Neff = 3.046, m_nu = [0, 0, 0.06] eV and Ob0 = 0.04897. Its reference: "Planck Collaboration 2018 … (Paper VI), Table 2 (TT, TE, EE + lowE + lensing + BAO)". — [astropy Planck18.ecsv](https://github.com/astropy/astropy/blob/main/astropy/cosmology/data/Planck18.ecsv)

### Inferences
- **Recomputed** (`cosmology_checks.out`, `hogg_identities.out`):
  - icosmo flat, open and closed: worst relative difference 5.4×10⁻⁸, 4.9×10⁻⁸ and 4.6×10⁻⁸ over 28 rows each, which is the tables' printed precision.
  - test_flat_z1: D_C = 3364.7987, D_A = 1682.3994, D_L = 6729.5974.
  - Kayser et al.: D_A(1, 2) = 646.2296866282 (agrees to 1×10⁻¹³); D_M(1, 2) = 1313.2232880, 5×10⁻⁸ from the quoted value.
  - w₀wₐ Mathematica: within 3×10⁻⁶ in every case. de_density_scale matches to all digits.
  - Neutrinos: the exact Fermi integral F(y) = ∫x²√(x² + y²)/(eˣ + 1)dx, normalised by 7π⁴/120, reproduces all Mathematica values to ≤3×10⁻⁶. Komatsu's fit, (1 + (Ay)^1.83)^{1/1.83} with A = 180ζ(3)/(7π⁴), errs by up to 2.7×10⁻³ (0.01 eV, z = 2).
  - Mattig's Eq. 17 against the integral: ≤7×10⁻¹⁶ for Ω_m = 0.3, 1 and 2.
  - Hogg Eq. 19 against S_k(χ₂ − χ₁)/(1 + z₂): they agree while √|Ω_k|χ < π/2 and diverge beyond it. For Ω_m = 3, Ω_Λ = 0, z₁ = 2, z₂ = 50, the S_k form gives 40.5883 Mpc and Hogg 30.0828 Mpc.
  - The crate's RXJ1131 doc test, 2334.06, is reproduced (2334.0558).
- **The Planck18 preset**, using astropy's parameters and exact neutrinos:
  - Ω_γ = 5.402015×10⁻⁵, Ω_ν = 1.439709×10⁻³, Ω_de = 0.688846, and Ω_m + Ω_ν = 0.31110 (Planck's 0.3111).
  - D_M(z* = 1089.80) = 13884.41 Mpc, so 100 r*/D_M = **1.04124** against Table 2's 1.04119 ± 0.00029 (0.17σ).
  - The same model with no radiation gives 1.03639.
  - D_C(1) = 3395.6346, D_A(1) = 1697.8173, D_L(1) = 6791.2692 Mpc (exact neutrinos). The Komatsu fit gives 3395.6343, 1697.8171 and 6791.2686.
- **Method for the crate.**
  - Replace the flat-only E(z) in `gravitational::cosmology` with a `Cosmology { h0, omega_m, omega_de, w0, wa, t_cmb, n_eff, m_nu: [f64; 3] }` struct, and a `PLANCK18` const holding astropy's parameters with Om0 excluding the massive neutrino.
  - Methods: `e(z)`, `comoving(z1, z2)` (χ in units of D_H), `transverse(z1, z2) = D_H S_k(χ)`, `angular_diameter(z1, z2)`, `luminosity(z)`, `time_delay(zd, zs)`.
  - Keep the existing free functions as the flat, radiation-free special case, so the RXJ1131 and H0LiCOW tests stand.
  - The crate's `integrate_log` (log-grid Simpson) suits it: integrate in ln(1 + z), which also handles z* ≈ 1090.
  - Doc tests: one icosmo open row and one closed row, the Mathematica w₀wₐ D_L, the Komatsu 0.25 eV ratio 429.924, Planck18's θ*, and the Mattig closed form.
- **Conventions and pitfalls.**
  - Ω_k = 1 − Ω_m − Ω_de − Ω_r. Positive Ω_k means open (sinh).
  - Planck's Ω_m includes the 0.06 eV neutrino; astropy's Om0 does not (0.30966 + 0.00144 = 0.3111).
  - astropy's Tcmb0 = 0 switches all radiation off.
  - N_eff is shared equally over three species, with T_ν = (4/11)^{1/3} T_cmb and 0.22710731766 = (7/8)(4/11)^{4/3} per species.
  - CPL dark energy: ρ_de ∝ (1+z)^{3(1+w₀+wₐ)} exp(−3wₐz/(1+z)).

### Gaps
- Komatsu et al. 2011 (Eqs. 24–26) was not opened; its content is taken from astropy's test comments, and the exact integral was verified against those Mathematica numbers instead.
- The sources of astropy's other non-flat D_M(z₁, z₂) values are not stated. I reproduce them to ≤6×10⁻⁷, but they are regression values, not external ones.

## Gap 4: lensing normalisation, finite sources and the H₀ registry

### Takeaway
- NFW: δ_c = (200/3)c³/[ln(1+c) − c/(1+c)], ρ_c = 3H²(z)/(8πG) and Σ_crit = c²D_s/(4πG D_d D_ds). For the concentration, use Dutton & Macciò's z-dependent Planck relation, with Duffy's WMAP5 relation as the alternative.
- Finite sources: Lee et al.'s exact one-dimensional integral for a uniform disc beats Gould's A(u)B₀(u/ρ) except for tiny ρ. Exact anchors: A(u = 0) = √(1 + 4/ρ²); A(u = ρ) in closed form, verified here; B₀(1) = 4/π.
- The H₀ registry should become method-tagged, with a common percentile convention. The new entries:
  - Planck 67.36 ± 0.54
  - SH0ES 73.04 ± 1.04, later 73.17 ± 0.86 and 73.49 ± 0.93
  - Local Distance Network 73.50 ± 0.81
  - CCHP 70.39 ± 1.22 ± 1.33 ± 0.70
  - DESI DR1 68.53 ± 0.80 and DR2 68.51 ± 0.58 (BAO + BBN)
  - GW170817 at 74.0 +16/−8 as median and percentiles, or 70.0 +12/−8 as MAP.

### Cited Findings
- **Wright & Brainerd's NFW normalisation.**
  - ρ(r) = δ_cρ_c/[(r/r_s)(1 + r/r_s)²], with ρ_c = 3H²(z)/(8πG) and r_s = r₂₀₀/c.
  - δ_c = (200/3)c³/[ln(1+c) − c/(1+c)] and M₂₀₀ = (800π/3)ρ_c r₂₀₀³.
  - Σ_c ≡ (c²/4πG)·D_s/(D_d D_ds), and κ = Σ/Σ_c.
  — [Wright & Brainerd 2000](https://arxiv.org/abs/astro-ph/9908213)
- **Duffy et al. 2008.**
  - c = A(M/M_pivot)^B(1 + z)^C, with M_pivot = 2×10¹² h⁻¹M☉. Table 1, NFW, full sample, z = 0–2: (A₂₀₀, B₂₀₀, C₂₀₀) = (5.71, −0.084, −0.47); z = 0: (5.74, −0.097, 0).
  - Virial: (7.85, −0.081, −0.71). Mean density: (10.14, −0.081, −1.01). Relaxed, 200: (6.71, −0.091, −0.44).
  - WMAP5 cosmology. The errors are 1σ from bootstrap.
  — [Duffy et al. 2008](https://arxiv.org/abs/0804.2486)
- **Dutton & Macciò 2014.**
  - log₁₀c = a + b log₁₀(M/[10¹² h⁻¹M☉]).
  - c₂₀₀: b = −0.101 + 0.026z, a = 0.520 + (0.905 − 0.520)exp(−0.617z^1.21).
  - c_vir: b = −0.097 + 0.024z, a = 0.537 + (1.025 − 0.537)exp(−0.718z^1.08).
  - At z = 0: a = 0.905 ± 0.001, b = −0.101 ± 0.001 (200); 1.025 and −0.097 (vir).
  - Planck cosmology: Ω_m = 0.3175, h = 0.671, σ₈ = 0.8344, n = 0.9624, Ω_b = 0.0490; Δ_vir(0) ≃ 104.2.
  — [Dutton & Macciò 2014](https://arxiv.org/abs/1402.7073)
- **Yoo et al. 2004.**
  - Gould's approximation: A_uni(u|ρ) ≃ A(u)B₀(u/ρ).
  - B₀(z) = (4/π)zE(k, z) with k = min(z⁻¹, 1), in Gradshteyn–Ryzhik notation.
  - Second order: A_uni = A(u)B₀(z)[1 + (ρ²/8)Q(z)], with Q = (1/3)[7 − 8z² − 4(1 − z²)F/E]. "We find numerically that −0.38 ≤ Q(z) ≤ 1, where the limits are saturated at z = 0.97 and z = 0".
  - Linear limb darkening: A ≃ A(u)[B₀ − ΓB₁], with B₁ = B₀ − (3/2)zρW₁. Limits: B₀ → 2z and B₁ → (2 − 3π/4)z as z → 0, and B₀ → 1, B₁ → 0 as z → ∞. Witt & Mao give "an exact evaluation … in closed, albeit cumbersome, form".
  — [Yoo et al. 2004](https://arxiv.org/abs/astro-ph/0309302)
- **Lee et al. 2009's exact uniform-source integral.** A_fs(u; ρ) = (1/πρ²)∫₀^π [u₂√(u₂² + 4) − u₁√(u₁² + 4)]dϑ, with u₂ = u cos ϑ + √(ρ² − u² sin²ϑ) and u₁ = 0 for u ≤ ρ. For u > ρ, u₁ = u cos ϑ − √(…) and the integral runs only to ϑ ≤ arcsin(ρ/u).
  - The Gould approximation, written as A_PS(u)(4u/πρ)E(ϑ_max, u/ρ), is "accurate only for a high-amplification event", with deviation ">2.5% for larger source".
  - Simpson's rule with n = 10 is within 0.5 % of Witt & Mao.
  - Witt & Mao's form "shows singularity when … u is similar to … ρ", so "Witt & Mao therefore derived a separate solution for the case of u = ρ".
  — [Lee et al. 2009](https://arxiv.org/abs/0901.1316)
- **H₀ primary values**, in km/s/Mpc:
  - **Planck 2018**, TT,TE,EE+lowE+lensing, flat ΛCDM: 67.36 ± 0.54. With BAO: 67.66 ± 0.42 (Table 2). — [Planck 2018 VI](https://arxiv.org/abs/1807.06209)
  - **SH0ES 2022**, Cepheid–SN Ia baseline including systematics: 73.04 ± 1.04. Cepheid + TRGB: 72.53 ± 0.99. ApJL, doi 10.3847/2041-8213/ac5c5b. — [Riess et al. 2022](https://arxiv.org/abs/2112.04510)
  - **SH0ES with the SMC as a fourth anchor**: 73.17 ± 0.86 (ApJ 973, 30, 2024; abstract read). — [Breuval et al. 2024](https://arxiv.org/abs/2404.08038)
  - **SH0ES JWST cross-check**: 72.6 ± 2.0 combining all JWST methods; the HST set gives 73.2 ± 0.9 (abstract read). — [Riess et al. 2024](https://arxiv.org/abs/2408.11770)
  - **SH0ES "Perfect Host"**: JWST + HST Cepheids 73.49 ± 0.93; with TRGB 73.18 ± 0.88 (abstract read). — [Riess et al. 2025](https://arxiv.org/abs/2509.01667)
  - **Local Distance Network consensus**: baseline 73.50 ± 0.81, which "differs by 7.1σ from flat ΛCDM with Planck+SPT+ACT and 5.0σ with BBN+BAO (DESI2)" (A&A 708, A166, 2026; abstract read). — [H0DN 2025](https://arxiv.org/abs/2510.23823)
  - **CCHP**: 70.39 ± 1.22 (stat) ± 1.33 (sys) ± 0.70 (σ_SN) from TRGB, 24 SN calibrators from HST + JWST. JWST-only TRGB: 68.81 ± 1.79 ± 1.32. JAGB: 67.80 ± 2.17 ± 1.64. — [Freedman et al. 2025](https://arxiv.org/abs/2408.06153)
  - **DESI DR1** (JCAP 02 (2025) 021): BAO + BBN 68.53 ± 0.80; + θ*: 68.52 ± 0.62. — [DESI 2024 VI](https://arxiv.org/abs/2404.03002)
  - **DESI DR2** (PRD 112, 083515): BAO + BBN 68.51 ± 0.58 in the equation and table, but **68.50 ± 0.58 in the conclusions**; + θ*: 68.45 ± 0.47; DESI + CMB: 68.17 ± 0.28. — [DESI DR2 II](https://arxiv.org/abs/2503.14738)
  - **GW170817** (Nature 551, 85), Table 1: 68.3 % symmetric 74.0 +16.0/−8.0; 68.3 % MAP **70.0 +12.0/−8.0**, the abstract's "maximum a posteriori and 68% credible interval". — [Abbott et al. 2017](https://arxiv.org/abs/1710.05835)
  - **TDCOSMO 2025** (A&A 704, A63): "constraints from **8** strongly lensed quasars". With Pantheon+, flat ΛCDM: 71.6 +3.9/−3.3. — [TDCOSMO 2025](https://arxiv.org/abs/2506.03023)
- **The existing crate text**: `H0_REGISTRY` names the TDCOSMO entry "the same six, new kinematics", and the docs say TDCOSMO 2025 "reanalyses the same six lenses". `KAPPA`'s doc says "θ_E = κMπ_rel". — [lensing.rs](/home/dan/code/FusionVault/quantum-gravimeter/crates/gravitational/src/lensing.rs)

### Inferences
- **Recomputed** (`finite_source.out`, `lens_on_limb.out`, `nfw_and_analogies.out`):
  - **B₀.** In SciPy, E(φ, k) = sin φ R_F(cos²φ, 1 − k²sin²φ, 1) − (k²/3)sin³φ R_D(…) with φ = arcsin min(1, 1/z) and modulus z. It needs no special case at z > 1; the closed form (4z/π)[zE(1/z²) − (z − 1/z)K(1/z²)] agrees to 2×10⁻¹⁵.
    - B₀(1) = 4/π = 1.2732395447351628.
    - B₀(0.1, 0.5, 0.9, 1.1, 2, 5) = 0.1994990586, 0.9342154577, 1.3426659198, 1.1667320471, 1.0346316184, 1.0050766018.
    - The maximum is 1.3430270547 at z = 0.90891.
    - Q runs from −0.3844 at z = 0.967 to 1.000000 at z → 0, confirming Yoo's quoted bounds.
  - **Exact uniform disc (Lee).**
    - A(0) = √(1 + 4/ρ²) to round-off at ρ = 0.01, 0.05, 0.5 and 1. For example, ρ = 0.5 gives 4.123105625617661.
    - **Lens on the limb:** A(ρ; ρ) = (1/π)[2/ρ + ((1 + ρ²)/ρ²)(π/2 + arcsin((ρ² − 1)/(ρ² + 1)))], matching to ≤1.4×10⁻¹² for ρ = 0.001 to 3. Values: ρ = 1 gives 1 + 2/π = 1.636619772367581; ρ = 0.5 gives 2.749075721239495; ρ = 0.05 gives 25.486000954805.
    - Gould's A·B₀ errs by −1.48×10⁻² at u = ρ/2, ρ = 0.5, and −6.5×10⁻⁶ at ρ = 0.01. Yoo's second order cuts the ρ = 0.05 error from 1.6×10⁻⁴ to 1.3×10⁻⁷.
    - Lee's Fig. 2 peaks: A(0.1; 0.5) = 4.086267 and A(0.01; 0.05) = 39.609821.
  - **Limb darkening.** B₁(0.05)/0.05 = −0.3545, tending to 2 − 3π/4 = −0.35619. B₁(0.5) = −0.09661963 and B₁(2) = 0.00735173, from a direct 2D quadrature of W₁.
  - **NFW.** δ_c(4, 5, 10) = 5271.147547, 8694.810191, 44778.661516.
    - Duffy c₂₀₀ at (M, z) = (10¹² h⁻¹M☉, 0) is 6.05233, against Dutton 8.03526 (c_vir 10.59254). At (10¹⁴, 0.5): Duffy 3.39750, Dutton 4.35381.
    - A 10¹⁴ h⁻¹M☉ halo at z_d = 0.3 in Planck18, with z_s = 1: ρ_c(z_d) = 1.180190×10⁻²⁶ kg/m³, r₂₀₀ = 1003.88 kpc, Dutton c₂₀₀ = 4.6457, δ_c = 7361.49, r_s = 216.09 kpc.
    - Σ_crit = 5.749602 kg/m² = 2753.170 M☉/pc², and κ_s = 2r_sδ_cρ_c/Σ_crit = 0.20151.
    - The mass identity 4πρ_cδ_c r_s³[ln(1+c) − c/(1+c)] = M₂₀₀ holds to 9×10⁻¹⁵.
- **Method and placement in `gravitational::lensing`.**
  - `sigma_crit(dd, ds, dds)`, `nfw_delta_c(c)`, `concentration_dutton(m200_h, z)` and Duffy as data.
  - `finite_source(u, rho)`, Lee's exact integral. Substitute u sin ϑ = ρ sin t for u > ρ to remove the square-root endpoint, and split at π/2 when u = ρ. This needs only √ and the existing quadrature, so `gravitational` does not need a `geopotential` dependency.
  - B₀ and B₁ belong in `physica` (Gap 5), where Carlson's R_F and R_D from `geopotential::special` are available.
  - Doc tests: A(0) = √(1 + 4/ρ²), the u = ρ closed form, the a ≫ ρ limit → Paczyński, and Lee's ≤0.5 % for Simpson n = 10.
- **H₀ registry.**
  - Widen `Hubble` with a `method` enum: `Cmb`, `Bao`, `DistanceLadder`, `TimeDelay(Assumption)` and `Siren`. Keep `value` as (median, +, −).
  - Use GW170817's **74.0 +16.0/−8.0**, the symmetric median interval, to match the struct's percentile convention, or note that 70.0 +12/−8 is a MAP with the highest-density interval.
  - For CCHP, store 70.39 with the quadrature total of 1.94 (computed), and keep the three components in `data`.
  - Fix the TDCOSMO 2025 text, which now says six lenses; the source says eight.
  - Suggested doc test: SH0ES and Planck differ by (73.04 − 67.36)/√(1.04² + 0.54²) = 4.85σ (computed). Riess 2022's text says "5σ".

### Gaps
- Witt & Mao 1994 (ApJ 430, 505) and Gould 1994 (ApJ 421, L71) are not on arXiv and were not opened. The u = ρ closed form above is my recollection of Witt & Mao, confirmed numerically against Lee's exact integral, not quoted from them.
- The Duffy 2008 arXiv v4 table may include the later erratum values; I did not check the MNRAS erratum.
- Of the four later SH0ES papers and H0DN, only the abstracts were read.

## Gap 5: the Phase 6 analogies 41i and 41j, and others the new modules enable

### Takeaway
41i should assert that every named temperature is one function of κ. Hawking's is exactly Unruh's at the horizon's surface gravity, c⁴/(4GM), and Gibbons–Hawking's is Unruh's at a = cH. The first law c²(dS/dM)⁻¹ = T_H also holds, with T_H·M = ħc³/(8πGk_B), equal to 6.170260×10⁻⁸ K M☉. 41j should assert θ_E² = κMπ_rel through the lens equation, and light_bending(M☉, AU) = κM☉[1 + (15π/16)GM☉/(c²AU)]. **The existing `KAPPA` doc's "θ_E = κMπ_rel" is wrong** (it should be θ_E² = κMπ_rel); `lens_mass` itself is right.

### Cited Findings
- **The Phase 6 plan.**
  - 41i: "One horizon temperature | gravitics::spacetime | ħκ/(2πk_B) for κ = a/c, c³/(4GM), H and a fibre's α | the fibre's 972.5 K beside the existing Unruh and Hawking tests".
  - 41j: "One lensing constant | gravitational::lensing ↔ gravitics::spacetime::light_bending | κ = 4GM☉/(c²·AU) | 8.144 mas per solar mass".
  - 41g: "One ring kernel … 4K(m)/√((a+ρ)² + z²) … K(0.5) = 1.854074677301372".
  — [Phase 6 plan](/home/dan/code/FusionVault/quantum-gravimeter/research/topics/cross-domain-frontiers/report.md)
- **The code.**
  - `horizon_temperature(kappa) = ħκ/(2πk_B)`. `unruh_temperature(a)` calls it with a/C, `hawking_temperature(m)` with c³/(4Gm), and `gibbons_hawking_temperature(h)` with h. The fibre doc test asserts 972.5 K at κ = 8×10¹⁴ s⁻¹ (Philbin et al. 2008, "reaches 10³ K").
  - `light_bending = 4x + c₂x²` with x = GM/(c²b) and c₂ = 15π/4 in GR.
  - `KAPPA = 4G/(c²AU)` and `lens_mass = θ_E/(κπ_E)`.
  — [spacetime.rs](/home/dan/code/FusionVault/quantum-gravimeter/crates/gravitics/src/spacetime.rs); [lensing.rs](/home/dan/code/FusionVault/quantum-gravimeter/crates/gravitational/src/lensing.rs)
- **The Sun's limb deflection** is 1.″75 (Table 1), and **Earth's** is 574 µas post-Newtonian with 0.6 µas from the quadrupole. — [Crosta & Mignard 2006](https://arxiv.org/abs/astro-ph/0512359)
- **Yoo's B₀ is a mean of 1/distance over the source disc**, B₀(z) ≡ zρW₀[z|ρ; x⁻¹]. — [Yoo et al. 2004](https://arxiv.org/abs/astro-ph/0309302)
- **The quadrupole deflection has a spin-2 pattern**, with components ∝ cos 2α and sin 2α. — [Crosta & Mignard 2006](https://arxiv.org/abs/astro-ph/0512359)

### Inferences
- **41i doc test** (`nfw_and_analogies.out`), in `physica` or in `gravitics::spacetime` itself, since only one crate is involved:
  - (i) `hawking_temperature(m) == horizon_temperature(c³/(4Gm))` and `== unruh_temperature(c⁴/(4Gm))`, to 0 ulp in my check.
  - (ii) `gibbons_hawking_temperature(h) == unruh_temperature(c·h)`.
  - (iii) T_H·M constant: 6.170260×10⁻⁸ K at 1 M☉ (computed with CODATA ħ and k_B, G, and GM☉ from IAU).
  - (iv) The first law, c²/(dS/dM) = T_H with S = k_B c³A/(4Għ) and A = 16πG²M²/c⁴, to 7×10⁻¹² by central difference.
  - (v) The fibre's 972.5 K already exists.
  - Optional cross-crate tie: `gibbons_hawking_temperature(PLANCK18.h0)` = **2.665597×10⁻³⁰ K** (computed), linking the new `cosmology::PLANCK18`.
  - Unruh needs 2.466083×10²⁰ m/s² for 1 K (computed).
- **41j doc test.**
  - (i) `KAPPA * SOLAR_MASS` = 8.143853 mas.
  - (ii) light_bending(M☉, AU, GR)/(KAPPA·M☉) − 1 = (15π/16)GM☉/(c²AU) = 2.907140×10⁻⁸. Assert the ratio of these two to 10⁻⁶; rounding limits it.
  - (iii) Solve the lens equation θ = α(θD_l)·D_ls/D_s with `light_bending` at Galactic distances (D_l = 4 kpc, D_s = 8 kpc, M = 0.3 M☉, π_rel = AU(1/D_l − 1/D_s)). It returns θ_E = √(KAPPA·M·π_rel) = 0.552625 mas, to about 10⁻⁹ from the second-order term.
  - (iv) KAPPA·M☉·AU/R☉ = 1.75119″, Table 1's 1.″75, with R☉ = 695700 km.
  - Fix the `KAPPA` doc to read θ_E² = κMπ_rel.
- **New analogies these modules enable**, with recomputed values:
  1. **"Earth's quadrupole in starlight"** (`geopotential` ↔ `gravitics::multipole_deflection`). The Earth coefficients already in `physica::analogies` (GM 3.986004415×10¹⁴, R 6378136.3 m, C̄₂₀ = −4.84169495×10⁻⁴, so J₂ = 1.082635903×10⁻³) give 573.70 µas and 0.6211 µas at the limb. Table 1 gives 574 and 0.6.
  2. **"The finite source is a disc of rings"** (`gravitational::lensing` ↔ `geopotential::special::ring`, extending 41g). B₀(z) = z⟨1/d⟩ over the unit disc, which is z times the in-plane Newtonian potential of a uniform disc of unit mass and radius. As a ring integral, (1/π)∫₀¹ r·4K(m)/(r + z) dr with m = 4rz/(r + z)². This matches the Carlson B₀ to 2×10⁻¹² at z = 0.25, 0.5, 2 and 5. Exact anchors: the potential at the rim is 4GM/(πa), so B₀(1) = 4/π; at the centre it is 2GM/a, so B₀ → 2z.
  3. **"One spin-2 pattern."** The J₂ deflection's modulus ∝ sin²ϑ, with orientation 2φ (Crosta Appendix B). Rotating the projected spin axis by ψ rotates the deflection by 2ψ, the same pattern as a gravitational wave's +/× antenna response. This needs an antenna-pattern function that `gravitational` does not yet have, so it is lower priority.
  4. **"The Einstein ring holds π R_E² Σ_crit"** (`lensing::KAPPA` ↔ the new `sigma_crit` with `cosmology` distances). M = π(θ_E D_l)²Σ_crit exactly; computed to 1×10⁻¹⁶. Microlensing's κ and galaxy lensing's Σ_crit are one constant.
  5. **"Radiation sets the CMB angle"** (`cosmology::PLANCK18` ↔ `metrology` CODATA σ_SB). Ω_γ from T_cmb = 2.7255 K gives 100θ* = 1.04124, against Planck's 1.04119 ± 0.00029.
  6. **"One spherical-harmonic engine, now Jupiter"**, extending the existing first analogy. The DUR2020 SHADR coefficients, normalised by 1/√(2n+1), go into `geopotential::Coefficients` and then `multipole_deflection`, which gives the 16516.561365 µas equatorial-limb series.

### Gaps
- No published source was opened for the Hawking temperature of one solar mass, 6.17×10⁻⁸ K, or for the Unruh acceleration of 1 K. These are recomputed from exact constants, so treat them as identities, not as published test values.
- The spin-2 analogy needs a new antenna-pattern function. I did not check whether `gravitational::memory` contains one.

## Conventions and pitfalls across the gaps

### Takeaway
Most failures in these gaps would be convention mismatches, not algorithm errors. Each convention below is pinned to its source.

### Cited Findings
- **Shadow angles.** JP2010, Johannsen and Chan measure the averaging angle about the black hole (tan α = y′/x′) but measure R from the displaced centre (D, 0). — [Johannsen & Psaltis 2010](https://arxiv.org/abs/1005.1931); [Chan et al. 2013](https://arxiv.org/abs/1303.5057)
- **Shadow size measures.** The EHT's δ uses the "median shadow diameter". M87 Paper VI quotes a "mean diameter" range. — [EHT Sgr A\* VI](https://arxiv.org/abs/2311.09484); [EHT M87 VI](https://arxiv.org/abs/1906.11243)
- **SHADR coefficients** are unnormalised (flag 0, 'U' in the file name), GM is in km³/s², and uncertainties are 3σ. The label swaps the m/n naming. — [DUR2020 LBL](https://atmos.nmsu.edu/PDS/data/jnogrv_0002/DATA/SHADR/GRV_JUGR_SHA_JUP_U_DUR2020_20200203_V01.LBL); [SIS](https://atmos.nmsu.edu/data_and_services/atmospheres_data/JUNO/logs/JUNO_GRAV_RDR_SIS.PDF)
- **Hogg's D_A12 is limited to Ω_k ≥ 0.** — [Hogg 1999](https://arxiv.org/abs/astro-ph/9905116)
- **Planck's Ω_m includes the 0.06 eV neutrino.** — [Planck 2018 VI](https://arxiv.org/abs/1807.06209)
- **GW170817's headline H₀** is a MAP with the minimal 68 % interval, not a median. — [Abbott et al. 2017](https://arxiv.org/abs/1710.05835)

### Inferences
- **Shadow measures to give.** Give ⟨R⟩ (JP angle about the origin) as the primary measure. Document the alternatives at a = 0.9, i = 90°: angle about the centre 5.031316, median 5.028450, areal 5.033226.
- **Curvature sign.** Ω_k > 0 is open (sinh) and Ω_k < 0 is closed (sin). Always take D_A12 from S_k(χ₂ − χ₁).
- **Units.** Keep the crate's SI (m, s⁻¹) and `H100`. Convert SHADR's km and km³ s⁻² at the boundary.
- **Normalisation.** C̄_n0 = C_n0/√(2n + 1) = −J_n/√(2n + 1).
- **The permanent tide.** State whether J₂ includes it: 14696.5735 includes the 0.0672×10⁻⁶ tide; 14696.5063 excludes it.
- **Lensing constant.** θ_E² = κMπ_rel, not θ_E = κMπ_rel.
- **Registry convention.** Store (median, +, −) consistently. MAP, HPD and stat/sys splits go in a note.

### Gaps
- None beyond those listed under each gap.
