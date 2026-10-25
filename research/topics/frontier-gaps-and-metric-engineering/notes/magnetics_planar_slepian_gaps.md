# Magnetics, planar transforms and Slepian functions: four gaps in the physica crates

Scope: the four Phase 1 gaps listed on FusionVault/physica issue #2 ("Not done in Phase 1"): items 3
(magnetic anomalies), 4 (planar transforms), 5 (Euler) and 7 (Slepian). Each section below answers the
three key questions for one gap: the most accurate testable method, the published values or synthetic
benchmarks that test it, and what mature codes do. Every published number is quoted from a source opened
in this session (equation or table named). Numbers marked "recomputed" or "ours" come from the scripts in
`gaps_scripts/` next to this file (python3 + numpy/scipy; outputs saved as `*.out`). Downloads are in
`gaps_downloads/<source>/`, read only; no downloaded code was run.

Crate conventions used throughout (read from the code, not changed): `magnetics::anomaly::body` uses
`B = (μ0/4π)(T − I tr T)·M` with `T` the Hessian of the unit-density potential `U = ∫dV/r`, so the
demagnetising field is `H = T·M/4π` inside and out and `N = −T/4π` inside an ellipsoid;
`geopotential::planar` has x north, y east, z down, power-of-two grids, `Θ_u = u_z + i(u_x k_x + u_y k_y)/|k|`;
`geopotential::slepian` uses 4π-normalised `P̄lm` with `(1/4π)∫(P̄lm cos mλ)² dΩ = 1`.

## Gap 1: induced plus remanent magnetisation and self-demagnetisation (method, test values, mature codes)

### Takeaway
The exact treatment for any ellipsoid (sphere, spheroid, triaxial, any orientation) is
`M = (I + χN)⁻¹(χH0 + M_r)` with `H0 = F/μ0` and `N` the body's constant internal demagnetising tensor
(Takahashi & Oliveira 2017, eqs 47–50; Clark, Saul & Emerson 1986 eq. 38 via that paper; Harmonica's
`get_magnetisation`); one Carlson integral, `N_c = (abc/3)·R_D(a², b², c²)`, gives every factor with no
shape branches and reproduces Chen, Brug & Goldfarb's (1991) Table I spheroid column and Takahashi &
Oliveira's Warrego and confocal numbers. Non-ellipsoids (prisms, finite cylinders) have non-uniform `M`;
the accurate route is a cell (collocation) integral equation built from the crate's own prism Hessians.

### Cited Findings
#### Formulas and primary sources
- Resultant magnetisation of an ellipsoid with isotropic susceptibility: `M = Λ(χH0 + M_R)`, `Λ = V(I + χÑ†)⁻¹Vᵀ` (eqs 49–50); anisotropic: `M = Λ(KH0 + M_R)`, `Λ = V(I + K̃Ñ†)⁻¹Vᵀ` (eqs 47–48), "consistent with that given by Clark et al. (1986, Eq. 38)". The remanent term is inside the same `Λ`, so remanence is demagnetised too. With `χ ≪ 1`, `M̆ = χH0 + M_R` (eq. 51), "the classical equation ... (Blakely, 1996, p. 89)". Along an axis, `P_i = ϑχH0/(1 + χñ_ii)` (eq. 63) — [Takahashi & Oliveira 2017, GMD 10, 3591](https://gmd.copernicus.org/articles/10/3591/2017/)
- Demagnetising factors as printed: prolate (`m = a/b > 1`) `ñ11 = [m/√(m²−1)·ln(m + √(m²−1)) − 1]/(m² − 1)`, `ñ22 = ñ33 = (1 − ñ11)/2` (eqs 30–31); oblate eq. 32 (arccos form; Emerson et al. 1985 used arctan); triaxial `a > b > c` through incomplete elliptic integrals F, E (eqs 25–29); `ñ11 + ñ22 + ñ33 = 1`; Stoner (1945) derivations — [Takahashi & Oliveira 2017](https://gmd.copernicus.org/articles/10/3591/2017/)
- Harmonica's implementation (commit a9cd442, 2026-10-16): triaxial `n_xx = abc/(√(a²−c²)(a²−b²))·(F(φ,k) − E(φ,k))`, `φ = arccos(c/a)`, `k = (a²−b²)/(a²−c²)` (scipy parameter), `n_zz = −abc/(√(a²−c²)(b²−c²))·E + b²/(b²−c²)`, `n_yy = 1 − n_xx − n_zz`; oblate `n_zz = [1 − q/√(1−q²)·arccos q]/(1 − q²)`, `q = c/b`; Clark et al. (1986) sign convention `ΔH = −N·M`; `get_magnetisation` solves `(I + N·χ)M = M_r + χH0` — [harmonica `_forward/ellipsoids/magnetic.py`](https://github.com/fatiando/harmonica/blob/a9cd44203308fd5db58acf74ac9ec3b7d1860136/src/harmonica/_forward/ellipsoids/magnetic.py)
- Any-ratio integral form: `N = ½(b/a)(c/a)∫0^∞ ds/[(s + c²/a²)R(s)]`, `R(s) = √((s+1)(s+b²/a²)(s+c²/a²))` along c (eqs Nellipsoid, Rs; Landau & Lifshitz 4.5, 4.25); Osborn (1945) "gives analytical solutions ... via the differences of incomplete elliptic integrals" — [Prozorov & Kogan 2018, arXiv:1712.06037 TeX](https://arxiv.org/abs/1712.06037)
- Finite cylinder, axial, χ = 0 (uniform M): `L_s = [2μ0n²/(3l²)]{(a²+l²)^½[l²F(k_s) + (a²−l²)E(k_s)] − a³}` (eq. 6), `k_s² = a²/(a²+l²)` (eq. 7), magnetometric `N_m = 1 − 2lL_s/(μ0πa²)` (eq. 15, n = 1 turn over length 2l); `N_m` and fluxmetric `N_f` defined by volume and midplane averages (eqs 4–5); cylinder `γ = l/a` — [Chen, Brug & Goldfarb 1991, IEEE Trans. Magn. 27, 3601 (NIST copy)](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=905398)
- Königsberger ratio `Q = M_rem/M_ind = M_rem/(χH)`; Q > 1 means remanence dominates; original reference Koenigsberger (1938), Terr. Magn. Atmos. Electr. 43, 299 — [Wikipedia, Koenigsberger ratio](https://en.wikipedia.org/wiki/Koenigsberger_ratio) (secondary; primary not opened)
- When to include self-demagnetisation: `‖δM‖/‖M‖ ≤ χ ñ_max` (eqs 58–59), so `χ_max = ε/ñ_max` (eq. 60), against the customary 0.1 SI — [Takahashi & Oliveira 2017](https://gmd.copernicus.org/articles/10/3591/2017/)
- Cell ("method of subsections") treatment of high susceptibility for non-ellipsoids is in Eskola & Tervo (1980), as listed in Takahashi & Oliveira's references — [GMD 10, 3591, reference list](https://gmd.copernicus.org/articles/10/3591/2017/) (paper itself not opened)
- Prisms: Aharoni (1998, J. Appl. Phys. 83, 3432) gives magnetometric factors of rectangular prisms; Prozorov & Kogan state the sum-to-one rule "is true only for ellipsoids" for their *effective* (diamagnetic, χ = −1) factors, and offer only an approximate `N ≈ 4ab/(4ab + 3c(a+b))` for cuboids — [Prozorov & Kogan 2018](https://arxiv.org/abs/1712.06037)

#### Published test values (all recomputed in `gaps_scripts/gap1_demag.py`, output `gap1_demag.out`)
- Chen et al. 1991 Table I, ellipsoid of revolution `N` (axis 2l, transverse 2a, γ = l/a): γ = 0.01 → 0.9845, 0.1 → 0.8608, 0.5 → 0.5272, 1 → 0.3333, 2 → 0.1736, 5 → 0.05582, 10 → 0.02029, 100 → 0.0004299. Recomputed (GMD eqs 30/32, Landau integral and Carlson form agree to 1e-15): 0.984490, 0.860804, 0.527200, 0.333333, 0.173564, 0.055821, 0.020286, 0.000430 — every printed digit matches — [Chen et al. 1991, Table I](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=905398)
- Chen et al. 1991 Table I, cylinder `N_m(0)`: 0.9650, 0.7967, 0.4745, 0.3116, 0.1819, 0.07991, 0.04119, 0.004232 for the same γ. Recomputed from eqs 6, 7, 15: 0.965040, 0.796676, 0.474490, 0.311577, 0.181864, 0.079907, 0.041193, 0.004232 — all match — [Chen et al. 1991](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=905398)
- Warrego orebody model (Table 2: a = 490.7 m, b = 69.7 m, c = 30.0 m, χ = 1.69 SI, B0 = (32610, 0, 39450) nT): "χ2 = 0.116 SI ... calculated by using Eq. (60) with a target relative error ε = 8%"; anomaly from ≈ −71 to ≈ 482 nT (≈ 553 nT peak-to-peak); with χ1 = 0.1 SI, `‖δM‖/‖M‖ ≈ 0.7 %` and ≈ 0.2 nT (≈ 0.6 %) anomaly difference. Recomputed: `ñ33 = 0.689521`, `χ_max = 0.08/0.6895 = 0.11602` — [Takahashi & Oliveira 2017, §4.3, Table 2](https://gmd.copernicus.org/articles/10/3591/2017/)
- Confocal ellipsoids (Table 1): Ellipsoid 1 = (900, 500, 100) m, χ = 1.2; Ellipsoid 2 = (≈1676.31, 1500, ≈1417.74) m with `u = 2×10⁶`, χ ≈ 0.014; volume ≈ 79 times larger, susceptibility ≈ 85 times smaller; field along a gives identical anomalies (eqs 63–64). Recomputed: a′ = 1676.31, c′ = 1417.74, χ′ = 0.01415, V2/V1 = 79.22, χ/χ′ = 84.8 — [Takahashi & Oliveira 2017, §4.2](https://gmd.copernicus.org/articles/10/3591/2017/)
- Sphere (factors "tend to 1/3", Stoner 1945; eq. 63 with ñ = 1/3): `M/(χH0) = 1/(1 + χ/3)` = 0.967742 at χ = 0.1, 0.75 at χ = 1, 0.230769 at χ = 10; the anomaly scales by the same factor, so ignoring demagnetisation overestimates a sphere's anomaly by χ/3 (3.3 % at 0.1 SI) — [Takahashi & Oliveira 2017](https://gmd.copernicus.org/articles/10/3591/2017/)
- Induced magnetisation for the crate's unit check: χ = 0.01 in F = 50 000 nT gives `χF/μ0 = 0.3979 A/m` (ours, from `H0 = F/μ0`).

#### Mature implementations
- Harmonica `ellipsoid_magnetic` (cites [Clark1986] and [Takahashi2018]) handles susceptibility tensors and remanence, internal and external demagnetising tensors (external via GMD eqs 34–37 with a stated sign difference), sphere/prolate/oblate/triaxial branches — [harmonica magnetic.py](https://github.com/fatiando/harmonica/blob/a9cd44203308fd5db58acf74ac9ec3b7d1860136/src/harmonica/_forward/ellipsoids/magnetic.py)
- Harmonica's tests are property tests, not published numbers: diagonal internal tensor, sum of factors = 1, demagnetisation lowers |M|, ellipsoid vs sphere and vs dipole at χ = 0.001, triaxial vs prolate/oblate limits — [harmonica test/ellipsoids/test_magnetic.py](https://github.com/fatiando/harmonica/blob/a9cd44203308fd5db58acf74ac9ec3b7d1860136/test/ellipsoids/test_magnetic.py)
- The GMD paper's code was Python in Fatiando a Terra (Zenodo doi:10.5281/zenodo.996479; github.com/pinga-lab/magnetic-ellipsoid) — [Takahashi & Oliveira 2017, §3](https://gmd.copernicus.org/articles/10/3591/2017/)

### Inferences
- **One formula, no branches.** `N_c = (abc/3)·R_D(a², b², c²)` follows from the Landau integral and Carlson's `R_D(x,y,z) = (3/2)∫dt/[(t+z)^{3/2}√((t+x)(t+y))]`; verified to 1e-15 against the integral for (3,2,1), (490.7,69.7,30), (1,1,10) and against all eight Chen spheroid entries (`gap1_demag.out`). `geopotential::special::carlson_rd` already exists, so `ellipsoid_factors(a, b, c) -> [f64; 3]` is three calls, with the sphere and both spheroid limits handled by `R_D` itself (Harmonica needs four branches and near-sphere tolerances).
- **Where:** `magnetics::anomaly`:
  - `magnetisation(chi, h0, remanent, n) -> [f64; 3]`: one 3×3 solve of `(I + χN)M = χH0 + M_r` (accept a 3×3 `chi` too, for anisotropy).
  - `ellipsoid_factors(a, b, c)`.
  - `koenigsberger(m_r, chi, h0)`.
  - Known `N` already implied by the crate: sphere `I/3`, infinite cylinder `(I − uuᵀ)/2`.
  - Feed `M` to the existing `sphere`, `cylinder` or `body`.
  About 40 lines plus tests. An ellipsoid's *external* field needs its gravity Hessian (confocal λ and Carlson forms). That is not in `geopotential::bodies` and is a separate, larger item.
- **Units:** the crate's ambient field is in tesla (`total_field(f, b)`), so `H0 = F/μ0`. Takahashi & Oliveira give H0 ≈ 18.7 A/m for their Table 1 field.
- **Non-ellipsoids, verified by `gap1_cells.py`:** with only prism Hessians (the crate has `prism`), the collocation system `(I − χ K)M = χH0`, `K_ij = T_j(c_i)/4π`, gives:
  - a single cube: exactly `1/(1 + χ/3)`;
  - staircase spheres, mean `M/(χH0)`: 0.75297, 0.75291 and 0.75193 at χ = 1 for 136, 552 and 1472 cubes (exact 0.75); 0.2403, 0.2396 and 0.2367 at χ = 10 (exact 0.2308).
  It converges slowly because of the staircase surface, and the 3n×3n dense solve costs O(n³). That is a reasonable `magnetics` item for χ ≳ 0.1 with fewer than about 2000 cells. For a single prism, the magnetometric factor (volume average of `−T/4π`) is a cheap approximation.
- **Risks:**
  1. Is remanence demagnetised? GMD eq. 49 and Harmonica both apply `(I + χN)⁻¹` to `M_r`. A lab-measured NRM of a high-χ sample is already the demagnetised value, so the API must say which one it takes.
  2. Q is ambiguous at high χ (intrinsic vs apparent).
  3. For prisms and cylinders, `−T/4π` at the centre is not the magnetometric factor. The cube centre gives 1/3, but Chen's cylinder γ = 1 has `N_m` = 0.3116.

### Gaps
- Osborn (1945, Phys. Rev. 67, 351), Stoner (1945), Clark, Saul & Emerson (1986, Explor. Geophys. 17, 189), Emerson et al. (1985), Blakely (1995), Aharoni (1998) and Eskola & Tervo (1980) were not opened: APS, Taylor & Francis and AIP paywalls, or a book. Chen et al.'s (1991) Table I stands in for Osborn's spheroid tables. The 2024 MDPI *Magnetism* update of Osborn's dataset (Kiss & Bakonyi) returned HTTP 403.
- No published closed-form prism value was opened, so a prism test can only use exact identities: a cube is 1/3, and the magnetometric factors sum to 1 because `tr T = −4π` inside any body. Aharoni's numbers remain unchecked.
- No published "anomaly of a high-susceptibility sphere" table was found beyond the exact `1/(1 + χ/3)` scaling.
- Takahashi & Oliveira's 0.7 % error at χ = 0.1 depends on their orientation angles (ε, ζ, η conventions); it was not recomputed.

## Gap 2: FFT padding, tapering and stabilised reduction to the pole (method, test values, mature codes)

### Takeaway
Harmonica pads by default: 25 % per side by edge (constant end-value) extension. On synthetic dipoles that is 2–12× better than no padding or zero padding whenever the field reaches the grid edge. A 50 % edge pad (which takes a power-of-two crate grid to the next power of two) with a cosine taper does best. Mirror padding helps only for narrow pads.

At low inclination, no Fourier-domain filter both stabilises the RTP and keeps amplitudes:
- the pseudo-inclination (Grant & Dodds / MAGMAP style, Seequent's formula) under-corrects N–S amplitudes by design;
- Wiener/Tikhonov damping needs a noise level;
- a damped equivalent-dipole layer (Silva 1986; Li & Oldenburg 2001 family) was the most accurate in our I = 10° benchmark (0.26 vs 0.38 Wiener and 0.84 plain, rms relative, with 1 % noise);
- reduction to the equator is stable and stays at the noise floor.

No published RTP error numbers for a 10° dipole were found; all numbers here are ours.

### Cited Findings
#### Padding: what mature codes do
- Harmonica `apply_filter`: "By default, use a padding width of 25% of each grid dimension. Fedi et al. (2012; doi:10.1111/j.1365-246X.2011.05259.x) suggest a padding of 100% but that seems exaggerated." Default mode is `"edge"`; the pad is removed after the inverse FFT; every transform (`upward_continuation`, `derivative_upward`, `reduction_to_pole`, …) takes `pad=True` — [harmonica filters/_filters.py](https://github.com/fatiando/harmonica/blob/a9cd44203308fd5db58acf74ac9ec3b7d1860136/src/harmonica/filters/_filters.py)
- Harmonica's `pad`/`unpad` are modified from xrft (MIT) and wrap `xarray.DataArray.pad`, offering numpy's modes constant, edge, linear_ramp (to `end_values`), maximum/mean/median/minimum, reflect, symmetric (with even/odd `reflect_type`) and wrap; the coordinates are extrapolated and the width stored as an attribute — [harmonica filters/_padding.py](https://github.com/fatiando/harmonica/blob/a9cd44203308fd5db58acf74ac9ec3b7d1860136/src/harmonica/filters/_padding.py)
- The Euler-inversion paper's notebook computes the upward derivative after `xrft.pad(..., mode="linear_ramp")` by half the grid size on each side, then unpads — [compgeolab/euler-inversion, code/synthetic-proof-of-concept.ipynb](https://github.com/compgeolab/euler-inversion)

#### Reduction to the pole: kernels and stabilisations
- Standard RTP kernel `1/(Θ_m Θ_f)`, `Θ = u_z + i(u_e k_e + u_n k_n)/|k|`, zero at k = 0; Harmonica applies it with no stabilisation — [harmonica `reduction_to_pole_kernel`](https://github.com/fatiando/harmonica/blob/a9cd44203308fd5db58acf74ac9ec3b7d1860136/src/harmonica/filters/_filters.py)
- Pseudo-inclination RTP (Seequent Oasis montaj REDP, Equation198): `L(θ) = [sin I − i cos I cos(D−θ)]² / ([sin²Ia + cos²Ia cos²(D−θ)]·[sin²I + cos²I cos²(D−θ)])`. Ia is the "Inclination for amplitude correction (never less than I)"; "if |Ia| < |I|, Ia = I"; it should be "large enough to counteract instability but small enough to minimize amplitude distortion"; θ is the wavenumber azimuth. The automatic rule (Equation199) is `Ia = min(90, max(0, int((0.008 I² − 1.71|I| + 80)/2.5)·2.5))`; the page cites no publication — [Seequent REDP help](https://help.seequent.com/Oasismontaj/2026.1/Content/gxhelp/fft2con/fft2con_redp.htm)
- Reduction to the equator (REDE, Equation193): the same denominator, with numerator `[sin I − i cos I cos(D−θ)]²·(−cos²(D−θ))`; "Default = ±20°" for Ia — [Seequent REDE help](https://help.seequent.com/Oasismontaj/2026.1/Content/gxhelp/fft2con/fft2con_rede.htm)
- Literature map (secondary, from an SBGf abstract), with the equation (1) RTP filter in polar wavenumber coordinates following Silva (1986):
  - Hansen & Pawlowski (1989): Wiener filtering;
  - Mendonça & Silva (1993): a truncated Taylor series;
  - Swain (2000): a pseudo-inclination in the denominator;
  - Li & Oldenburg (2001): "inverted the Fourier transform of the observed magnetic data with explicit regularization";
  - Silva (1986): equivalent sources;
  - Li et al. (2014): a regularised equivalent layer with positivity;
  - MacLeod et al. (2003): the total gradient (analytic signal) at low latitude.
  At I ≈ 0° the RTP "generates elongated artifacts in the direction of the declination" — [Melo & Barbosa 2019, SBGf](https://sbgf.org.br/mysbgf/eventos/expanded_abstracts/16th_CISBGf/Reduction%20to%20the%20pole%20in%20Fourier%20domain%20%e2%80%93%20good%20and%20bad%20filtering%20of%20real%20data%20in%20Brazil.pdf)
- Harmonica also lists a convolutional equivalent layer (Takahashi, Oliveira & Barbosa 2020, Geophysics 85(6), G129) in its references, the scalable regular-grid form of equivalent sources — [harmonica doc/references.rst](https://github.com/fatiando/harmonica/blob/a9cd44203308fd5db58acf74ac9ec3b7d1860136/doc/references.rst)

#### Synthetic benchmarks (ours, scripts `gap2_rtp_padding.py`, `gap2_rtp_low_inclination.py`, `gap2_rte_noise.py`)

The Python mirror first reproduces the crate's pinned values: upward continuation 6.924e-4 rms (crate 6.9e-4) and RTP at I = 45° 0.0114 (crate 0.011).

Padding, rms relative error on the original 128×128 grid at 10 m:

| case | none | zero 25 % | edge 25 % (Harmonica) | mirror 25 % | edge 50 % + taper |
|---|---|---|---|---|---|
| up 30 m, dipole 50 m deep (decays) | 6.9e-4 | 1.7e-3 | 3.9e-4 | 1.5e-3 | 2.5e-4 |
| up 30 m, dipole 300 m deep (edges carry signal) | 1.8e-2 | 1.6e-2 | 1.7e-3 | 2.5e-3 | 1.1e-3 |
| up 30 m, dipole 100 m deep at 1/8 from an edge | 2.8e-2 | 1.6e-2 | 2.3e-3 | 3.5e-2 | 1.8e-3 |
| RTP I = 45°, dipole 300 m deep | 0.33 | 0.18 | 0.082 | 0.099 | 0.042 |

- Mirror (`symmetric`) at 50–100 % is no better than none in the edge cases (1.8e-2, 3.5e-2, 0.27).

Low-inclination RTP, one induced dipole 100 m deep, 48×48 grid at 25 m, 25 % edge pad, rms relative error against the true pole anomaly:
- Noise is 1 % of the peak. The "floor" is rms(noise)/rms(truth): 0.070 at I = 10°.
- `*` marks a damping chosen as the best of a scan, an upper bound on any practical choice.

| I, D | plain | pseudo-I, Ia = I+10° | pseudo-I, Seequent auto Ia | Wiener* | equivalent layer* |
|---|---|---|---|---|---|
| 45°, 0 clean / noisy | 0.0103 / 0.139 | 0.146 / 0.189 | = plain (Ia = I) | 0.0103 / 0.138 | 0.0012 / 0.081 |
| 20°, 0 | 0.024 / 0.261 | 0.248 / 0.296 | 0.43 / 0.444 (Ia = 47.5°) | 0.023 / 0.235 | 0.020 / 0.151 |
| 10°, 0 | 0.092 / 0.838 | 0.29 / 0.403 | 0.564 / 0.571 (Ia = 62.5°) | 0.069 / 0.38 | 0.016 / 0.264 |
| 10°, 30° | 0.095 / 0.682 | 0.305 / 0.395 | 0.566 / 0.571 | 0.095 / 0.384 | 0.011 / 0.246 |
| 5°, 0 | 0.361 / 3.14 | 0.301 / 0.582 | 0.593 / 0.599 (Ia = 70°) | 0.141 / 0.424 | 0.0057 / 0.32 |
| 2°, 0 | 2.24 / 17.8 | 0.294 / 0.841 | 0.603 / 0.608 (Ia = 75°) | 0.183 / 0.433 | 0.0020 / 0.336 |

- Reduction to the equator against the true I = 0 anomaly (same grid and noise): at I = 10°, exact `Θ_eq²/Θ_f²` gives 0.0084 clean and 0.105 noisy (floor 0.115); Seequent REDE with Ia = 20° gives 0.049 and 0.107. At I = 5°: 0.0081 / 0.119 (floor 0.123).

### Inferences
- **Padding: method and where.** Add `Grid::padded(&self, values, mode) -> (Grid, Vec<f64>)` and a crop to `geopotential::planar`, or a `pad: bool` on `filter`, matching Harmonica:
  - mode edge (each edge row/column repeated);
  - width up to the next power of two, so 50 % per side on the crate's power-of-two grids (128 → 256 exactly);
  - an optional cosine (Tukey) roll-off over the pad.
  It costs 4× the FFT work and needs no new dependency. Zero padding should not be the default: it was *worse* than none (1.7e-3 vs 6.9e-4) for a field that decays inside the grid, because it puts a step at the boundary.
- **RTP: method and where.** Three additions to `planar`, in order of value per line:
  1. `reduce_to_equator(values, m, f)`: stable at any inclination, about 10 lines; the kernel is `Θ_eq,m Θ_eq,f / (Θ_m Θ_f)`, with the horizontal unit vectors at the same declinations.
  2. A damped pole kernel `conj(t)/(|t|² + ε)` with `t = Θ_mΘ_f`: one optional ε argument on `pole_kernel`. ε ≈ 0.02–0.03 was best for 1 % noise; choose it by GCV or discrepancy in practice.
  3. The Seequent pseudo-inclination, as an `Ia` argument. It is a reproducible industry convention, but in our test it biases amplitudes (0.146 rms at I = 45° with Ia = 55° on clean data), so document it as a display filter.
  The most accurate option, the damped equivalent-dipole layer, is a dense `n_cells³` solve. In our test it is O(2304³) via one eigendecomposition. A 128×128 grid needs the convolutional (BTTB + FFT + conjugate-gradient) form; that is a separate, larger item and could live with `gravimetry::local` or a new `planar::equivalent_layer`.
- **Benchmark caveat.** The equivalent-layer "source" is the same family as the truth (induced dipoles), which flatters it. Its clean-data advantage partly comes from avoiding the periodic edge, not only from the regularisation.
- **Test values for the crate** (our numbers, reproducible from the scripts; they lock behaviour rather than check against literature):
  - padded continuation 1.7e-3 vs unpadded 1.8e-2 for the 300 m dipole;
  - RTE at I = 10° below 1 % on clean data;
  - damped RTP at I = 10° with 1 % noise below 0.4, against 0.84 plain.
  Exact identities that need no benchmark: RTE∘(inverse RTE) = identity; the damped kernel → plain as ε → 0; the pseudo-I kernel equals plain when Ia = I.

### Gaps
- Blakely (1995, ch. 12), Grant & Dodds (1972, unpublished MAGMAP notes), Hansen & Pawlowski (1989), Mendonça & Silva (1993), Swain (2000), Li & Oldenburg (2001) and Li (2008, The Leading Edge) were not opened: book or SEG/CSIRO paywalls. Their formulas and published comparisons are not quoted.
- Fedi et al. (2012) was seen only through Harmonica's comment.
- Seequent's pages do not give the provenance of the auto-Ia rule.
- No published numeric RTP-error benchmark for a dipole at I = 10° was found. All such numbers above are ours.

## Gap 3: 3-D gridded and moving-window Euler, structural-index estimation, extended Euler (method, test values, mature codes)

### Takeaway
Go from the crate's profile solver to the 3-D Reid et al. (1990) least squares in two steps: add the y column, then a moving window ranked by the depth variance. The most accurate published single-window method is Uieda et al.'s (2025, GJI) Euler inversion:
- it is a total-least-squares Gauss–Newton fit that also predicts the data;
- it picks the structural index by minimum weighted misfit;
- it is now Harmonica's `EulerInversion`.

We reproduced its proof-of-concept numbers to the metre: Euler deconvolution (14626, 11865, −1553) m and Euler inversion (15045, 12029, −2663) m, against the paper's (15045, 12028, −2663) and a true depth of 3000 m. Extended Euler (Mushayandebvu et al. 2001) adds dip and susceptibility for 2-D contacts and sheets; its equations were not opened.

### Cited Findings
#### Formulas and algorithms
- Euler's equation (3-D, up positive): `(e_i − e0)∂f/∂e + (n_i − n0)∂f/∂n + (u_i − u0)∂f/∂u = η(b − f_i)`. Structural index table "after Reid & Thurston (2014)":
  - point or sphere: SI 3 (magnetic), 2 (gravity);
  - line, cylinder or thin-bed fault: 2 and 1;
  - thin sheet edge, thin sill or thin dyke: 1 and 0.
  Harmonica `EulerDeconvolution` solves the 4-parameter least squares and stores `covariance_ = χ²(AᵀA)⁻¹`; it "Assumes a single data window"; SI 0 is not supported — [harmonica `_euler_methods.py`](https://github.com/fatiando/harmonica/blob/a9cd44203308fd5db58acf74ac9ec3b7d1860136/src/harmonica/_euler_methods.py)
- Deconvolution as printed: `A p = c` with rows `[−∂x f_i, −∂y f_i, −∂z f_i, −η]` and `c_i = −x_i∂x f_i − y_i∂y f_i − z_i∂z f_i − ηf_i`; `p = (AᵀA)⁻¹Aᵀc` (eq. deconv-p); `C_p = σ̂0²(AᵀA)⁻¹`, `σ̂0² = ‖c − Ap‖²/(N − 4)` (eq. deconv-cov). The estimate "is valid only if the contents of the Jacobian matrix A are assumed to be error-free" — [Uieda et al. 2025, GJI 241(3) 1536, doi:10.1093/gji/ggaf114 (paper TeX in the repo)](https://github.com/compgeolab/euler-inversion)
- Euler inversion (Algorithm 1):
  - start from `p0` = the deconvolution solution and `d0 = 0.9 d_obs`;
  - each Gauss–Newton step updates both parameters and predicted data, with merit `‖W^½ r‖ + ν‖e‖`;
  - if the merit rises, undo the step and stop; also stop when the relative merit change falls below δ.
  Algorithm 2: run for each η and choose the η of minimum `WRMSE = √([d_o − d_η]ᵀW[d_o − d_η])`. Algorithm 3: overlapping windows; for each η, sort the solutions by the estimated variance of z_o and keep the first γ·M_η. Default weights (1, 0.1, 0.1, 0.025) for (f, ∂e, ∂n, ∂u) — [Uieda et al. 2025](https://github.com/compgeolab/euler-inversion)
- Harmonica `EulerInversion` defaults: `structural_index=(1, 2, 3)` (picks the best by data misfit), `max_iterations=20`, `tol=0.1`, `euler_misfit_balance=0.1`. Each step uses a diagonal `Q⁻¹`, so the cost per iteration is O(N) plus 4×4 inverses. Unlike the paper's code it does not undo a step that raises the merit — [harmonica `_euler_methods.py`](https://github.com/fatiando/harmonica/blob/a9cd44203308fd5db58acf74ac9ec3b7d1860136/src/harmonica/_euler_methods.py); [paper code/euler.py](https://github.com/compgeolab/euler-inversion)
- Moving and expanding windows in legacy Fatiando:
  - `EulerDeconvMW`: a grid of window centres of a given size; per window, the uncertainty is `√diag(inv(Hessian))[0:3]` and solutions are ranked by its norm, keeping the best `keep` fraction (default 0.2).
  - `EulerDeconvEW`: expanding windows about a centre, keeping the single best.
  The equation is formulated after Reid et al. (1990), Geophysics 55(1), 80–91 — [fatiando/gravmag/euler.py](https://github.com/fatiando/fatiando/blob/master/fatiando/gravmag/euler.py)
- "The usual strategy for multiple sources, introduced by Reid et al. (1990), is to use a moving window"; FitzGerald et al. (2004) review selection criteria; Silva et al. (2001) showed solutions are biased according to the window scheme; for estimating SI "we refer the reader to Barbosa et al. (1999) and Melo et al. (2013)" — [Uieda, Oliveira & Barbosa 2014, TLE 33(4) 448 (source manuscript.md)](https://github.com/pinga-lab/paper-tle-euler-tutorial)
- SI estimation: "Barbosa et al. (1999) proposed that the minimum correlation between base level estimates and magnetic anomaly gives the correct structural index, Melo et al. (2013) used this approach in 3D data"; Melo & Barbosa find that without interfering anomalies "the minimum variation of base level estimates indicates the correct structural index" — [Melo & Barbosa 2017, SBGf](https://sbgf.org.br/mysbgf/eventos/expanded_abstracts/15th_CISBGf/session/POTENTIAL%20FIELDS%20METHODS/Correct%20structural%20index%20defined%20by%20base%20level%20estimates%20in%20Euler%20deconvolution%20.pdf)
- Extended Euler: the conventional Euler equation solved jointly with "the rotational constraint equation"; "Solving both equations jointly (extended Euler) gives distance, depth, dip and susceptibility", assuming no remanence. A solution is kept if the two depth estimates differ by less than a user "Max % error". Cites Mushayandebvu, van Driel, Reid & Fairhead (2001), Geophysics 66(3), 814–823 — [Seequent GM-SYS Extended Euler help](https://help.seequent.com/GMSYS/2026.1/Content/ss/gmsys_2d/gx/pdepth_gxs/exeuler_gx.htm)
- Mushayandebvu et al. (2004) exploit the presence of 2-D sources (SI 1 magnetic / 0 gravity) in a window to estimate them correctly; finite-difference Euler (Gerovska et al. 2005) estimates non-integer SI — [Uieda et al. 2025, Introduction](https://github.com/compgeolab/euler-inversion)

#### Published test values
- Proof of concept:
  - setup: dipole 5e11 A m² at (15000, 12000, −3000) m, I = −30°, D = 15°; grid [0, 25e3]×[0, 20e3] at 300 m, height 800 m; 10 nT noise with `default_rng(42)`; base level 100 nT;
  - Euler inversion: (15045, 12028, −2663) m, base 93 nT, in 6 iterations;
  - Euler deconvolution: (14626, 11865, −1553) m, base 94 nT.

  Recomputed in `gaps_scripts/gap3_euler.py` (output `gap3_euler.out`) from the equations, with the notebook's data recipe:
  - deconvolution: (14626, 11865, −1553) m, base 93 nT;
  - inversion: (15045, 12029, −2663) m, base 91 nT, 6 iterations.

  The positions match to ≤ 1 m; the base levels differ by 1–2 nT, attributed to derivative implementation details (xarray/xrft) — [Uieda et al. 2025, `paper/variables/synthetic-proof-of-concept.tex`](https://github.com/compgeolab/euler-inversion)
- Other published test setups:
  - noise sweep: 0–40 nT in 0.2 nT steps, dipole at (15000, 11000, −5000) m, 500 m spacing; upward-coordinate errors at the largest noise (40 nT) of 4628 m (EI with unit weights), 4252 m (ED) and 2128 m (EI with weights 1, 0.1, 0.1, 0.025), as defined in the notebook `code/synthetic-noise-levels.ipynb`;
  - moving windows: 10 sources, I = −30°, D = −20°, 50 nT noise, base 1000 nT, regional gradient (0.02, −0.03) nT/m, 10 km windows with 5 km step, keep γ = 0.3 (ED), 0.35 (finite-difference ED), 0.25 (EI).

  These depend on the authors' random draws and Harmonica's derivative code, so they make good regression targets only with that recipe — [Uieda et al. 2025, `paper/variables/*.tex`](https://github.com/compgeolab/euler-inversion)
- Harmonica's tests:
  - Euler deconvolution with FFT derivatives of a noise-free dipole (1e14 A m², I = −40°, D = 15°, at (10e3, 15e3, −10e3) m, 500 m grid over [−100e3, 100e3]×[−80e3, 80e3]): location and base level within `rtol=0.01`;
  - with analytic point-mass derivatives: `atol=1e-3, rtol=1e-3`;
  - Euler inversion with 20 nT noise: horizontal within 1 %, depth 5 %, base 1 %, and SI 3 (dipole) or 2 (point mass) selected.

  — [harmonica `test/test_euler_methods.py`](https://github.com/fatiando/harmonica/blob/a9cd44203308fd5db58acf74ac9ec3b7d1860136/test/test_euler_methods.py)
- Ours, with noise-free analytic derivatives (central differences, h = 1 m), source at (10, 8, −2) km:
  - deconvolution at the true SI is exact to the derivative step: depth error 4.8e-4 m (dipole), 1.9e-4 m (point mass), 1.4e-4 m (line mass at azimuth 30°);
  - inversion WRMSE picks the right index each time: dipole 5.2e-3 / 1.9e-3 / 2.4e-9 for SI 1/2/3 → SI 3; point mass → SI 2; line mass → SI 1.

### Inferences
- **What to add, in `geopotential::planar` beside `euler`:**
  1. `euler3(points: &[[f64; 3]], f, fx, fy, fz, n) -> ([f64; 4], [[f64; 4]; 4])`. This is the profile solver with a y column, the same column scaling and Cholesky, and `σ̂0²(AᵀA)⁻¹` returned so windows can be ranked (≈ 25 lines).
  2. `euler_windows(grid, values, fx, fy, fz, n, size, step, keep)`, following Algorithm 3 / `EulerDeconvMW`: rank by `var(z0)`.
  3. `euler_inversion(...)` with SI selection over {1, 2, 3}, following Algorithm 1 with the paper's undo-on-merit-increase rule. That is about 80 lines with no new linear algebra beyond 4×4 inverses.
- **Derivatives:** use `filter` with `i k_x`, `i k_y` and `|k|` (sign for z down) on the Gap 2 padded grid. The paper itself used a half-grid linear-ramp pad for the vertical derivative and finite differences horizontally.
- **Test plan:**
  - the exact recovery above (analytic derivatives);
  - the Uieda et al. proof-of-concept as a regression: needs the same noise vector, which numpy's `default_rng(42)` cannot provide in Rust, so store the 68×84 noise array or use the noise-free variant;
  - the SI selection table.
- **Risks:**
  1. Sign conventions: the crate's z is down; Harmonica's and the paper's are up.
  2. 2-D sources make `AᵀA` singular along strike. Our line-mass case needed a minimum-norm solve; the crate's Cholesky would panic. Mushayandebvu et al. (2004) constrain this.
  3. Euler deconvolution's depth bias with noisy derivatives is large: −1553 vs −3000 m in the paper's test.
  4. Window selection rules strongly shape the solution clouds (Silva et al. 2001, via the tutorial).

### Gaps
- Reid et al. (1990), Barbosa, Silva & Medeiros (1999), Mushayandebvu et al. (2001, 2004), Reid & Thurston (2014) and Melo et al. (2013) were not opened (SEG paywall). Their exact equations are not quoted: notably extended Euler's rotational equation and Reid's acceptance criteria. Barbosa et al. are described only through Melo & Barbosa (2017).
- The published Euler-inversion results for the noise sweep and windows were not recomputed.

## Gap 4: Slepian functions of arbitrary regions, basin averages and leakage (method, test values, mature codes)

### Takeaway
Build the concentration matrix `D_lm,l′m′ = (1/4π)∫_R P̄lm P̄l′m′ T_m T_m′ dΩ` the way Simons' `kernelc.m` does: Gauss–Legendre in colatitude, exact longitude integrals over the region's longitude intervals at each node. Improve it by splitting at vertex and extremum colatitudes and clustering the nodes there. Our prototype:
- reproduces a rotated 40° cap's eigenvalues to 8e-15 at L = 18 (rotation invariance);
- gives `tr D = (L+1)²A/4π` to 1e-13;
- gives exact areas for polygons.

The basin average is the first row of the same matrix: Simons' `plm2avg`, and Swenson & Wahr (2002) eq. 17 with exact kernel coefficients eq. 18. Swenson & Wahr's Gaussian and Lagrange-optimal kernels (eqs 30–34, 44) and Harig & Simons' (2012) truncated Slepian basis (L = 60, 0.5° buffer, N = 20 for Greenland) are the leakage treatments.

### Cited Findings
#### Formulas
- `D_lm,l′m′ = ∫_R Y_lm Y_l′m′ dΩ` (eq. Dlmlmpdef, unit-normalised Y), `Dg = λg` (eq. eigen1), `1 > λ1 ≥ … ≥ λ_(L+1)² > 0`, and `N = Σλ = tr D = ∫_R D(r̂, r̂) dΩ = (L+1)²A/(4π)` (eq. shannon). Rounded N "roughly separate[s]" the λ > 0.5 functions. Cutoff wavelengths 6200, 3200, 2200, 1600 km for L = 6, 12, 18, 24. The commuting tridiagonal (Grünbaum) matrix is derived only for the axisymmetric polar cap (§"Concentration within an axisymmetric polar cap") — [Simons, Dahlen & Wieczorek 2006, SIAM Rev. 48, 504; arXiv:math/0408424 TeX](https://arxiv.org/abs/math/0408424)
- Harig & Simons restate D (eq. S4), the eigenproblem (S3), the Shannon number (S6) and the truncated expansion `d(r̂) ≈ Σ_{α≤N} d_α g_α(r̂)` for r̂ ∈ R (S7) — [Harig & Simons 2012, PNAS 109, 19934, SI (Princeton reprint)](https://geoweb.princeton.edu/people/simons/PDF/reprints/PNAS-2012.pdf)
- Swenson & Wahr:
  - eq. 2: `Δσ = (aρ_E/3)Σ(2l+1)/(1+k_l) P̄lm(cos θ)(ΔC cos mφ + ΔS sin mφ)`;
  - eq. 15: `ϑ` = 1 inside the basin, 0 outside;
  - eq. 16: `Δσ̄_region = (1/Ω_region)∫Δσ ϑ dΩ`;
  - eq. 17: `Δσ̄_region = (aρ_E/(3Ω_region))Σ(2l+1)/(1+k_l)(ϑ^c_lm ΔC_lm + ϑ^s_lm ΔS_lm)`;
  - eq. 18: `ϑ = (1/4π)Σ P̄lm(ϑ^c cos mφ + ϑ^s sin mφ)`;
  - eq. 8: `K_l = aρ_E(2l+1)/(3(1+k_l))`;
  - eq. 27: approximate kernel `Δσ̃ = Σ K_l/Ω_region (W^c ΔC + W^s ΔS)`;
  - Gaussian (eqs 30–34): `W(γ) = (b/2π)·exp[−b(1−cos γ)]/(1−e^{−2b})`, `b = ln 2/(1 − cos(r_½/a))`, `W_lm = 2πW_l ϑ_lm`, and the recursion `W0 = 1/2π`, `W1 = (1/2π)[(1+e^{−2b})/(1−e^{−2b}) − 1/b]`, `W_{l+1} = −(2l+1)/b·W_l + W_{l−1}`;
  - leakage, eq. 42: `var(ε_lkg) = ∫(W̄ − ϑ)²dΩ/∫ϑ²dΩ = (1/(4πΩ_region))Σ[(W^c − ϑ^c)² + (W^s − ϑ^s)²]`;
  - Lagrange-optimal kernel, eq. 44: `W_lm = [1 + λK_l²B_l²/(2l+1)]⁻¹ ϑ_lm`, with λ the unique positive root of eq. 46.

  — [Swenson & Wahr 2002, JGR 107(B9) 2193, doi:10.1029/2001JB000576](https://doi.org/10.1029/2001JB000576) (equations read from the publisher's equation images, saved in `gaps_downloads/jgr_swenson2002/eqs*.png`)

#### Mature implementations
- Simons' `kernelc.m` (slepian_alpha, a4bf0c1):
  - nodes: `nGL = max(ngl, 2·Lmax)` Gauss–Legendre points (default ngl = 200) in cos θ between the region's northern and southern colatitudes;
  - at each node, longitude intervals from the boundary crossings (`dphregion` → `phicurve`);
  - longitude integrals done exactly (`coscos`, `sinsin`, `sincos`), with Legendre products accumulated;
  - the result divided by 4π, so `Klmlmp(1)` is the fractional area;
  - "NOT FOR POLAR PATCHES! ... NOT WITHOUT MODIFICATIONS FOR REGIONS CONTAINING THE NORTH POLE OR THE SOUTH POLE";
  - an 'alternative' path masks a fine grid (`geoboxcap`, oversampling factor 2⁵) and does masked harmonic transforms;
  - `localization.m` diagonalises the kernel.

  `dphregion` treats coastlines "as occupying flat Cartesian geometry (which is, admittedly, a bit of a limitation)": edges straight in (lon, lat) — [slepian_alpha kernelc.m](https://github.com/csdms-contrib/slepian_alpha/blob/a4bf0c1e2e9354c828b51a76b9a4649994a07090/kernelc.m), [dphregion.m](https://github.com/csdms-contrib/slepian_alpha/blob/a4bf0c1e2e9354c828b51a76b9a4649994a07090/dphregion.m)
- Region outlines (`REGIONS/australia.m` etc.) come from `regselect` on a coastline file in `$IFILES`, not shipped in the repo — [slepian_alpha REGIONS](https://github.com/csdms-contrib/slepian_alpha/tree/a4bf0c1e2e9354c828b51a76b9a4649994a07090/REGIONS)
- Basin averaging in `plm2avg.m`: the "integration vector (also the first row/column from Klmlmp in KERNELC)" times the coefficients gives the integral, and `A = theINT/miniK(1)` the average. `integratebasis.m` integrates Slepian functions over the region in fractional sphere area — [slepian_delta plm2avg.m](https://github.com/csdms-contrib/slepian_delta/blob/0e8354306fbf41f2cddd836c04d9c9958b623a59/plm2avg.m), [integratebasis.m](https://github.com/csdms-contrib/slepian_delta/blob/0e8354306fbf41f2cddd836c04d9c9958b623a59/integratebasis.m)
- slepian_bravo holds linear inverse problems with scalar Slepian functions (Simons & Dahlen 2006, GJI 166, 1039) — [slepian_bravo README](https://github.com/csdms-contrib/slepian_bravo/blob/5f5f4e24aa2fe3443821f4ca2c1f06beee52f241/README.md)
- SHTOOLS `ComputeDMap(dh_mask, lmax, sampling, degrees)`:
  - the region is a Driscoll–Healy 0/1 mask, and each `Y_lm × mask` is transformed once to give a row of D;
  - the elements are approximate: "n should be chosen such that further increases in n do not alter the returned eigenvalues", with a bandwidth of about 5× lmax suggested;
  - it prints the pixel area over D(0,0) as a check.

  — [SHTOOLS ComputeDMap](https://shtools.github.io/SHTOOLS/pycomputedmap.html)
- Harig & Simons's procedure:
  - Slepian basis at L = 60 "matches the bandwidth of the GRACE data products"; coastlines buffered by 0.5° (chosen by synthetic recovery);
  - truncation at the Shannon number leaves 20 functions, the 12th "still concentrated to λ = 86.9%";
  - GIA (Paulson et al.) projected onto the same basis and removed; load Love numbers used;
  - mass is the sum over the 20 functions of coefficient × regional integral;
  - result −199.7 ± 6.3 Gt/y (2σ), acceleration −8.68 ± 4.1 Gt/y².

  — [Harig & Simons 2012](https://geoweb.princeton.edu/people/simons/PDF/reprints/PNAS-2012.pdf)

#### Published test values
- Table "tebel" (continents, area A/4π in %, rounded Shannon N for L = 6, 12, 18, 24):

  | region | A/4π | N at L = 6, 12, 18, 24 |
  |---|---|---|
  | Greenland | 0.44 % | 0, 1, 2, 3 |
  | Australia | 1.50 % | 1, 3, 5, 9 |
  | South America | 3.50 % | 2, 6, 13, 22 |
  | North America | 3.98 % | 2, 7, 14, 25 |
  | Africa | 5.78 % | 3, 10, 21, 36 |
  | Eurasia | 9.98 % | 5, 17, 36, 62 |

  At L = 18: Australia `λ5 = 0.607`, `λ12 = 0.049`; North America `λ11 = 0.642`, `λ12 = 0.596`; the 40° cap gives N = 42.

  Recomputed `(L+1)²A/4π` (`gap4_slepian.out`): e.g. Africa 2.83, 9.77, 20.87, 36.12 and Eurasia 4.89, 16.87, 36.03, 62.38. All 24 printed integers are the rounded values — [Simons, Dahlen & Wieczorek 2006, Table tebel and §"Continental concentration"](https://arxiv.org/abs/math/0408424)
- Greenland with a 0.5° buffer at L = 60: "The rounded Shannon number is n = 20"; λ12 = 0.869. Our inference: this needs A/4π between 0.524 % and 0.551 % (19.5/3721 to 20.5/3721), against SDW's 0.44 % unbuffered — [Harig & Simons 2012, Fig. S2 caption and SI text](https://geoweb.princeton.edu/people/simons/PDF/reprints/PNAS-2012.pdf)
- Swenson & Wahr results:
  - to reach 1 cm satellite error, a 500 km disc needs a Gaussian roll-off of 125 km and a 200 km disc 215 km;
  - for six North American basins of 2.2×10⁵–1.4×10⁶ km², at 200 km signal correlation and 1 cm satellite error, the leakage ratio is ≈ 0.4–0.1 with Gaussian kernels and ≈ 0.2–0.05 with the Lagrange kernel.

  These depend on preliminary JPL error degree amplitudes B_l (personal communication) and HYDRO1K basin masks, so they cannot be reproduced — [Swenson & Wahr 2002, §5–6](https://doi.org/10.1029/2001JB000576)

#### Our prototype checks (`gaps_scripts/gap4_slepian.py`, output `gap4_slepian.out`)
- Polar 40° cap at L = 18 through the general code reproduces the crate's pinned values: trace 42.228978 (Shannon 42.228978), λ42 = 0.5792, 42 eigenvalues above ½, λ1–3 within 3e-8 of 1.
- The same cap centred at colatitude 60°, longitude 30°: max |λ_rotated − λ_polar| = 3.9e-3 with 24 nodes per panel, 7.7e-15 with 48, 8.9e-15 with 96; trace error ≤ 1.6e-13.
- Polygons, with kernelc's lon-lat edge convention:
  - a lat-lon box (10–50°E, 20°S–15°N): `D00 = A/4π = 0.033379954912677` exactly, `tr D/((L+1)²A/4π) − 1 = −7.8e-16`;
  - a triangle ((0,0), (40,10), (15,50)): `4πD00 = 0.2601573187`, equal to an independent area sum; trace identity to 7.9e-13.
- Basin average from the first row: the mean of `P̄20` over the 40° polar cap is 1.512553001033004, against the closed form `√5·c(1−c²)/(2(1−c))` = 1.512553001033005 (c = cos 40°).
- Swenson & Wahr recursion vs eq. 33 by quadrature:
  - r_½ = 200 km: agrees to 4.4e-10·W0 up to l = 120;
  - r_½ = 500 km: the forward recursion drifts, by more than 1e-3·W0 at l = 118, while 2πW0 = 1 in both.

### Inferences
- **Method for the crate (`geopotential::slepian`):**
  - `region(lmax, intervals: impl Fn(f64) -> Vec<(f64, f64)>, breaks: &[f64], nodes) -> Vec<Slepian>`, where `Slepian` holds a full `(l, m, cos/sin)` coefficient vector rather than one order;
  - helpers `cap(colat, lon, radius)` with analytic intervals, and `polygon(&[(lon, lat)])` with lon-lat edges as in kernelc, breaks at vertex colatitudes;
  - reuse `Legendre`, `gauss_legendre`, and the exact longitude integrals: cos·cos, sin·sin and sin·cos antiderivatives with zero-frequency cases.
  A polar region works if its intervals return the full [0, 2π] near the pole; kernelc forbids this case. Great-circle edges would need only a different crossing formula.
- **Basin average (gravimetry's "not included yet"):**
  - `basin_kernel(lmax, region) -> Coefficients` holds `ϑ̂_lm = (1/4π)∫_R P̄lm{cos, sin} dΩ`, the D row for (0,0);
  - `average = Σ factor_n(C̄_nm ϑ̂c_nm + S̄_nm ϑ̂s_nm)/ϑ̂_00` with `Functional::WaterHeight { love }`'s factor `Rρ_E(2n+1)/(3ρ_w(1+k′_n))`. This is Swenson & Wahr eq. 17 in the crate's normalisation (their ϑ^c = 4πϑ̂, Ω_region = 4πϑ̂_00), checked by the `P̄20` cap test;
  - Gaussian smoothing `W_lm = 2πW_l ϑ̂_lm`: compute W_l by quadrature of eq. 33 (or stop the recursion where it drifts) rather than by forward recursion alone.
- **`Coefficients` conversion:** a unit vector c in the crate's 4π normalisation is already a set of `C̄nm`/`S̄nm` with `(1/4π)∫g² = 1`, so `Slepian::to_coefficients(gm, radius) -> Coefficients` is `Coefficients::zeros` plus `set(l, m, c, s)`. The projection of a field is the coefficient dot product (Parseval). Harig & Simons's regional mass is `Σ_{α≤N} d_α·4π(g_α·ϑ̂)`.
- **Cost and risks:**
  - Assembling D costs about K·(L+1)⁴/4 multiply-adds for K colatitude nodes: about 10¹⁰ for L = 60 and K ≈ 3000, seconds in Rust. The single kernel row costs K·(L+1)², which is trivial.
  - Diagonalisation is the bottleneck. The crate's cyclic Jacobi is O(n³) per sweep at n = 3721, too slow; it needs Householder tridiagonalisation + implicit QL, or subspace iteration for the top ~N eigenpairs.
  - There is no commuting tridiagonal matrix for arbitrary regions, so eigenvectors inside the clusters near 1 and 0 are ill-determined. Basin averages need no eigenvectors.
  - Coastline and buffer data are needed for real tests.
  - Edge convention (lon-lat vs great circle) and dateline wrap must be documented.
  - Polygons with hundreds of vertices make hundreds of panels. Without splitting, convergence is algebraic, as kernelc accepts with 200 nodes.

### Gaps
- No public outline data: Simons' coastline file in `$IFILES` is not in the repos. So the published continental eigenvalues (Australia λ5 = 0.607, North America λ11 = 0.642) and Greenland's N = 20 at L = 60 could not be recomputed; only Shannon-number arithmetic and exact identities were.
- Swenson & Wahr's numeric results depend on unpublished error spectra and HYDRO1K masks.
- The PNAS main text was read from the authors' reprint, not the PNAS page.
- SHTOOLS was read through its documentation page only, not its Fortran source.
