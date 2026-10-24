# Plasmoids and self-organising magnetised plasmas: Taylor relaxation, FRCs, spheromaks, reconnection, Helion and TAE

Conventions. SI units, CODATA 2018 constants, μ₀ = 1.25663706212×10⁻⁶ H m⁻¹, 1 keV = 1.160452×10⁷ K, so 100 million °C is 8.62 keV.
"computed" means I recomputed the value in `plasmoid_checks.py` in this folder. Run it with `python3 -I plasmoid_checks.py`. It uses only the standard library and checks its results with asserts against the published values.
"Published" values are quoted exactly as the source gives them. Labels: **[mainstream]**, **[contested]**, **[refuted]**.
For company items, [contested] means the claim comes from the company or its press releases and has had no independent peer-reviewed check. It does not mean anyone has shown the claim to be false.
Scope limits. Fusion reactivities ⟨σv⟩ for D–³He and p–¹¹B, dense plasma focus and Z-pinch plasmoids (including Lerner's "Focus Fusion") are left to the sibling notes in this folder.

## Q1. Self-organisation: Taylor relaxation, helicity, Woltjer's theorem and spheromak eigenvalues

### Takeaway
Magnetic helicity K = ∫A·B dV is exactly conserved in ideal MHD. In a weakly resistive, turbulent plasma it is still nearly conserved, while magnetic energy decays much faster. The state of minimum energy at fixed K is the linear force-free field ∇×B = λB with a single constant λ (Woltjer 1958). Taylor (1974) argued that real plasmas relax to that state. Everything a calculator needs then reduces to Bessel-function eigenvalues. The reversed-field pinch reverses at λa = 2.405 (Θ = 1.20). The helical state takes over at λa = 3.11 (Θ = 1.555). A cylindrical spheromak has (λa)² = 3.8317² + (πa/h)², and a spherical one has λa = 4.4934. A spheromak in a cylinder becomes tilt-unstable for h/a above about 1.67.

### Cited Findings

**Formula sheet**

| # | Quantity | Exact form (SI) | Source | Published test value | Computed |
|---|---|---|---|---|---|
| 1 | Magnetic helicity [Wb²] | K = ∫_V A·B dV. It is gauge-invariant when B·n̂ = 0 on ∂V. | Woltjer 1958 PNAS 44, 489, [doi:10.1073/pnas.44.6.489](https://doi.org/10.1073/pnas.44.6.489); Taylor 1986 RMP 58, 741, [doi:10.1103/RevModPhys.58.741](https://doi.org/10.1103/RevModPhys.58.741) (both cited as the origin; I did not re-read them in this session) | — | — |
| 2 | Woltjer minimum-energy state [mainstream] | Minimise W = ∫B²/(2μ₀) dV at fixed K. This gives ∇×B = λB with constant λ [m⁻¹]. Then W = λK/(2μ₀), because ∫B² = ∫B·∇×A = λ∫A·B for a closed flux conserver. | Woltjer 1958 (above) | — | For K = 10⁻³ Wb² and λ = 10 m⁻¹, W = 3.979 kJ (example only) |
| 3 | Taylor relaxation hypothesis [mainstream, with known limits] | A turbulent plasma with small resistivity conserves only the global K, approximately, and relaxes to the state in row 2. | Taylor 1974 PRL 33, 1139, [doi:10.1103/PhysRevLett.33.1139](https://doi.org/10.1103/PhysRevLett.33.1139); Taylor 1986 RMP (above) | — | — |
| 4 | Bessel-function model (BFM) of the cylindrical Taylor state | B_z = B₀J₀(μr), B_θ = B₀J₁(μr). Pinch parameter Θ = B_θ(a)/⟨B_z⟩ = μa/2. Reversal parameter F = B_z(a)/⟨B_z⟩ = μa·J₀(μa)/(2J₁(μa)). | standard algebra from row 2 | Field reversal (F = 0) occurs at the first zero of J₀ | μa = 2.4048, so **Θ = 1.2024** |
| 5 | Onset of the helical (m = 1) Taylor state | Above a critical μa the minimum-energy state is no longer axisymmetric | [Aspect-ratio scaling of single helical states in RFPs, arXiv:1509.07305](https://arxiv.org/pdf/1509.07305) (from the search-result summary) | The state is axisymmetric for μa below about **3.11**, and the m = 1 solution appears at μa = 3.11. Above μa ≈ 3.176 the first unstable mode is m = 1 with ka ≈ 0.86. | At μa = 3.11: **Θ = 1.555**, F = −1.545 |
| 6 | Cylindrical spheromak eigenvalue (flux conserver of radius a and length h, lowest axisymmetric mode) | (λa)² = j₁,₁² + (πa/h)², where j₁,₁ = 3.831706 is the first zero of J₁ | Chandrasekhar–Kendall/Taylor solution, standard (form as in Taylor 1986) | — | λa = 4.9550 at h/a = 1; 4.2686 at h/a = 1.67; 3.9722 at h/a = 3 |
| 7 | Tilt limit of a cylindrical spheromak | The spheromak tilts when the lowest m = 1 eigenvalue falls below the axisymmetric one | [APS DPP 2008 abstract, "Spheromak aspect-ratio effects on poloidal flux amplification"](https://archive.aps.org/dpp/2008/tp6/91) | "the tilt-mode stability limit (1.67) for an isolated spheromak". NIMROD runs were taken to L/R = 1.6. | The crossing therefore lies at λa ≈ 4.27 |
| 8 | Spherical spheromak (Woltjer/Rosenbluth–Bussac) | B ∝ spherical Bessel j₁(λr). B_r(a) = 0 requires j₁(λa) = 0, which is the same as tan x = x. | standard (I did not retrieve a published table) | — | **λa = 4.4934** |

**Measured spheromak performance (for calibration)** [mainstream]
- CTX, 1990, solid flux conserver: Ti0 = Te0 = 0.18 keV, n = 4.5×10¹⁹ m⁻³, τE* = 0.2 ms, n·T·τ = 1.6×10¹⁵ keV m⁻³ s — [Wurzel & Hsu, arXiv:2105.10954, Table 7](https://ar5iv.labs.arxiv.org/html/2105.10954)
- SSPX, 2007, shot 17524: Te0 = 0.5 keV (Ti inferred from Te), n = 2.25×10²⁰ m⁻³, τE* = 1 ms, n·T·τ = 1.1×10¹⁷ keV m⁻³ s — [same table](https://ar5iv.labs.arxiv.org/html/2105.10954)
- In coaxial-gun spheromak formation, a kink instability amplifies the poloidal flux. Its onset "agree[s] quantitatively with the Kruskal-Shafranov limit" — [arXiv:physics/0304104](https://arxiv.org/pdf/physics/0304104)
- In a high-aspect-ratio flux rope (L/R ≈ 11) the spheromak tilts and forms "a twisted double-helix Taylor state" — [arXiv:1310.2228](https://arxiv.org/pdf/1310.2228)

### Inferences
- A Rust calculator needs only J₀, J₁, the zeros 2.404826 and 3.831706, the root 4.493409 of tan x = x, and a root finder. The script uses a power series and bisection. All values are checked to 10⁻⁹.
- The BFM F–Θ curve is the mainstream benchmark. Real RFPs depart from it: experiments run at F ≈ −0.2 near Θ ≈ 1.6. The pure BFM gives F = −1.55 at Θ = 1.555, so measured profiles are only partially relaxed. This is why "modified BFM" fits exist. The departure is known qualitatively; I did not retrieve a sourced F–Θ dataset.
- The tilt limit in row 7 can be checked as a calculator test. Tabulate λa(h/a) for the m = 0 mode (row 6) and confirm the published 1.67 corresponds to λa = 4.27. The m = 1 eigenvalue in a finite cylinder needs a mixed-mode solution and was not computed here.

### Gaps
- I did not read Taylor 1986 RMP, Bellan's *Spheromaks*, Rosenbluth & Bussac 1979 or Finn, Manheimer & Ott 1981 in this session, so the published m = 1 eigenvalue curve and derivation of the 1.67 limit are unretrieved. The 1.67 number comes from a conference abstract.
- The resistive helicity decay rate dK/dt = −2∫ηJ·B dV and the gun helicity-injection rate dK/dt ≈ 2V_gunΦ_gun are standard results. I did not retrieve primary sources for them, so I list them here only as leads.

## Q2. Field-reversed configurations: equilibrium, the s-parameter, stability, confinement scalings, merging FRCs

### Takeaway
Combining radial pressure balance with axial force balance gives the profile-independent result ⟨β⟩ = 1 − x_s²/2, with x_s = r_s/r_c (Tuszewski 1988 Eq. 10, "perhaps the most useful relation in FRC physics"). Typical x_s ≈ 0.4–0.6 gives ⟨β⟩ ≈ 0.82–0.92. Stability and confinement are organised by s, the number of ion gyroradii between the field null and the separatrix. Experiments up to 1988 had s ≲ 2, and a reactor needs s ~ 20–40. Particle confinement roughly follows the R²/ρ_ie yardstick. Flux confinement is 3–5 times worse than classical. Merging two supersonic FRCs heats the plasma (Ti + Te > 0.5 keV in C-2) and amplifies flux more than tenfold. TAE's 2025 Norm result skips merging entirely.

### Cited Findings

**Equilibrium and stability formula sheet** (all from [Tuszewski, "Field reversed configurations", Nucl. Fusion 28, 2033 (1988), open PDF on Zenodo](https://zenodo.org/records/1235740), [doi:10.1088/0029-5515/28/11/008](https://doi.org/10.1088/0029-5515/28/11/008)) [mainstream]

| # | Quantity | Exact form | Tuszewski ref. | Published test value | Computed |
|---|---|---|---|---|---|
| 1 | Radial pressure balance | p_M = p + B²/(2μ₀). p_M is the peak pressure at the field null, equal to B_e²/(2μ₀). | Eq. (8) | — | — |
| 2 | Field-null radius | r_s = √2·R (R is the null, r_s the separatrix radius) | Eq. (9) | Fig. 23 confirms it experimentally | — |
| 3 | Average beta | ⟨β⟩ = (2/r_s²)∫₀^{r_s}(p/p_M) r dr = **1 − x_s²/2**, x_s = r_s/r_c | Eq. (10), first derived by Barnes et al. | "volume averaged beta values between 0.5 and 1" | x_s = 0.3/0.4/0.5/0.6/0.7 gives ⟨β⟩ = 0.955/0.920/**0.875**/0.820/0.755 |
| 4 | Poloidal-flux bounds | (x_s/√2)² < φ/(πR²B_e) < x_s/√2. The low-flux (LFSB) and high-flux (HFSB) sharp-boundary profiles are the limits. | Sec. 3.1.1 | — | x_s = 0.5 gives 0.125 < φ/(πR²B_e) < 0.354 |
| 5 | s-parameter | s = ∫_R^{r_s} r dr/(r_s ρ_i) ("approximate number of ion gyroradii between R and r_s"). For uniform T_i this equals s = φ/(2π r_s ρ_ie B_e), where ρ_ie is the ion gyroradius in the external field. | Eqs (11)–(12) | — | — |
| 6 | s bounds | With S* = R/ρ_ie: S*x_s²/(4√2) < s < S*x_s/4 | Sec. 3.1.1 | — | x_s = 0.5, S* = 10 gives 0.442 < s < 1.250. The script checks these bounds against rows 4 and 5. |
| 7 | Rigid-rotor profile | B = B_e tanh[K(r²/R² − 1)], with K chosen to satisfy Eq. (10). It gives flux close to the HFSB value. | Sec. 3.1.1 | — | — |
| 8 | 2-D equilibrium limit on separatrix beta | β_s < (1 − x_s²)² (Spencer et al.) | Sec. 3.1.2 | — | — |
| 9 | Adiabatic flux compression of elongated FRCs | T ∝ B_e^{4/5}, n ∝ B_e^{6/5}, from nT ∝ B², nV = const and nT ∝ V^{−5/3}. Full profile-dependent laws are in Table V. | Sec. 3.3 | Flux compression alone gives net heating efficiencies of "a few per cent". Wall compression keeps x_s constant and can be efficient. The Tor-liner device reached a total neutron yield "about 2 × 10⁸". | — |
| 10 | Stability status 1988 | Experiments have "s ≲ 2", while "s ~ 20–40 might be required for FRC reactors that use D–T fuel". Kinetic tilt theory predicts a "gradual transition to observable unstability for s > 3". LSX was designed for s ~ 8. | Sec. 1.2.2 | — | — |
| 11 | n = 2 rotational instability | Threshold rotation α ~ 1–1.2 (neutral-beam spectroscopy). It was eliminated by weak multipole fields. | Secs 1.2.1, 4.2 | — | — |
| 12 | Lifetime vs Alfvén time | "The configuration lifetimes exceed the Alfvén transit times by at most a factor of 100" | Sec. 1.2.2 | — | — |

**Confinement scalings** (Tuszewski 1988, Table VII, τ_N ∝ ℓ^a x_s^b r_c^c n^d T^e, using pressure balance to eliminate B_e) [mainstream as empirical fits; extrapolation contested]

| Scaling | a | b | c | d | e |
|---|---|---|---|---|---|
| R²/ρ_ie (empirical yardstick, Hamasaki & Krall) | 0 | 2 | 2 | 0.5 | 0 |
| x_s R²/ρ_ie (ℓ/r_c)^{3/4} (Ishimura) | 0.75 | 3 | 1.25 | 0.5 | 0 |
| TRX-1 | 0 | — | 2.2 | 0.6 | 0.3 |
| FRX-C/T | 0.2 | — | 1.4 | 0.5 | 0.4 |
| TRX | 0 | — | 3.6 | 0.9 | 0 |
| LHD (theory) | 0.6 | 3.3 | 2 | 0.5 | −0.7 |
| Bohm | 0 | 2 | 2 | 0.5 | −0.5 |
| Classical | 0 | 2 | 2 | 0 | 1.5 |

The OCR of the TRX-1, FRX-C/T and TRX rows lost the x_s exponent, so those cells show "—".
- Tuszewski: the gap between empirical and classical scalings is "approximately τ_N/τ_N(classical) ~ n^0.5/T^1.5". Flux confinement shows "moderate flux anomaly factors, mostly in the range 3–5" over classical Spitzer diffusivity. He uses η⊥/μ₀ (m² s⁻¹) = [T_e(eV)/100]^{−3/2} — [Tuszewski 1988 Secs 5.3.1–5.3.2](https://zenodo.org/records/1235740). Computed check: the NRL Spitzer η⊥/μ₀ at 100 eV with lnΛ = 12 is 0.98 m² s⁻¹, which matches the shorthand.
- The later review is L. C. Steinhauer, "Review of field-reversed configurations", Phys. Plasmas 18, 070501 (2011), [doi:10.1063/1.3613680](https://doi.org/10.1063/1.3613680) (cited by [arXiv:2606.15211](https://arxiv.org/pdf/2606.15211)). I could not read its content; see Gaps.

**Measured FRCs, from the independent ARPA-E compilation** ([Wurzel & Hsu, arXiv:2105.10954, Table 7](https://ar5iv.labs.arxiv.org/html/2105.10954); * marks a peak value inferred from the volume average) [mainstream]

| Device | Year | Ti0 (keV) | n (m⁻³) | τE* (s) | n·T·τ (keV m⁻³ s) |
|---|---|---|---|---|---|
| LSX (s ≈ 2) | 1993 | 0.547 (Te 0.253) | 1.30×10²¹* | 1×10⁻⁴ | 7.1×10¹⁶ |
| FRX-L | 2003 | 0.09 | 4×10²² | 3.3×10⁻⁶ | 1.2×10¹⁶ |
| TCS / TCSU | 2005/2008 | 0.025 / 0.1 | 6.5×10¹⁸ / 1.3×10¹⁹ | 4×10⁻⁵ / 7.5×10⁻⁵ | 6.5×10¹² / 9.7×10¹³ |
| C-2U | 2017 | 0.68* | 2.47×10¹⁹ | 2.4×10⁻⁴ | 4.0×10¹⁵ |
| C-2W (shot 118340) | 2021 | 3.5* | 1.30×10¹⁹ | 5×10⁻³ | 2.3×10¹⁷ (computed 2.275×10¹⁷) |

**Merging and colliding FRCs** [mainstream]
- C-2 (Tri Alpha): colliding and merging two θ-pinch-formed plasmoids gave "over tenfold" poloidal-flux amplification and "total temperature (Ti+Te) exceeding 0.5 keV", with flux confinement "approaching classical values" — [Binderbauer et al., PRL 105, 045003 (2010), PubMed abstract](https://pubmed.ncbi.nlm.nih.gov/20867853). Lifetimes exceeded 2 ms. Temperatures and fluxes "depend strongly on the speed of the initial individual plasmoids" — [APS DPP 2010 abstract](https://archive.aps.org/dpp/2010/di3/1).
- MSNW Inductive Plasma Accelerator (Slough et al., the forerunner of Helion): FRCs were ejected at "~250 km/s" into a central chamber, where they collided, merged and were magnetically compressed — [APS DPP 2007 abstract](https://meetings-archive.aps.org/dpp/2007/to4/14/).
- Computed: at 250 km/s the kinetic energy is 652 eV per deuteron and 326 eV per proton. This matches the ~0.5 keV total temperatures after merging if most of the directed energy thermalises.
- TAE Norm (2025): an FRC was formed by neutral-beam injection alone, with no θ-pinch sections and no merging. The change "from open field line to fully a formed field-reversed configuration is complete within ~ 10 ms". Separatrix radius 0.4 m, length 2 m, trapped flux about 6 mWb, plasma current 300–350 kA, total energy 9 kJ — [Roche et al., Nat. Commun. 16, 3487 (2025), PMC11993635](https://pmc.ncbi.nlm.nih.gov/articles/PMC11993635).

### Inferences
- Computed consistency check on Norm (order of magnitude only). A long current sheet with B(axis) ≈ −B_e gives B_e ≈ μ₀I/(2ℓ) = 0.11 T. Then φ/(πR²B_e) = 0.217, and the flux bounds (row 4) place x_s between 0.31 and 0.66, so ⟨β⟩ is between 0.78 and 0.95. This brackets the paper's "typical average β ~90%". The stored energy 1.5·⟨β⟩·p_M·V ≈ 6.5 kJ, using the cylinder volume as an upper bound, is the same order as the published 9 kJ. The published value includes anisotropic fast-ion energy. The point is that the published numbers are mutually consistent with Eq. (10). It does not measure B_e.
- An FRC calculator should take (B_e, r_c, r_s, ℓ_s, n, T_i, T_e) and return ⟨β⟩, R, the flux range, ρ_ie, the s range, the Alfvén time and the R²/ρ_ie yardstick. Watch the gyroradius convention: authors differ by √2 in the thermal speed (√(2T/m) vs √(T/m)). Fix it explicitly in code and state it in tests.
- Tuszewski's flux-compression law also tests merging-then-compression schemes such as Helion's (see Q4). Going from T₀ to T₁ needs B_e to rise by (T₁/T₀)^{5/4}.

### Gaps
- I could not read Steinhauer 2011. The search for the s/E (elongation-normalised) stability criterion failed with a tool error, so I have no sourced number for the s/E boundary.
- Tuszewski Fig. 35 gives τ_N against R²/ρ_ie, but I could not read a numeric prefactor from the OCR. Only exponents are given above.
- I did not retrieve TS-3/TS-4 (Tokyo) spheromak and FRC merging data.

## Q3. Magnetic reconnection: Sweet–Parker, the plasmoid instability and Lundquist numbers

### Takeaway
The Lundquist number is S = μ₀Lv_A/η. Sweet–Parker reconnection is S^{−1/2} slow. That is 1.8% of v_A in MRX (S ≈ 3×10³) but 2×10⁻⁷ in the solar corona (S ≈ 10¹³). It would take about 200 days to release a flare's flux, against about 6 minutes at a rate of 0.01. Above S_c ~ 10⁴, Sweet–Parker sheets break into plasmoid chains. The linear growth rate scales as S^{1/4}/τ_A and the number of plasmoids as S^{3/8} (Loureiro et al. 2007). The global rate then becomes roughly independent of S at about S_c^{−1/2} ~ 0.01 (Uzdensky et al. 2010). Bhattacharjee et al. 2009 measured a rate "weakly dependent" on S, about ten times Sweet–Parker.

### Cited Findings

**Formula sheet** [mainstream]

| # | Quantity | Exact form | Source | Published test value | Computed |
|---|---|---|---|---|---|
| 1 | Lundquist number | S = μ₀L_CS v_A/η, where v_A = B_R/√(μ₀ρ) and L_CS = εL is the sheet half-length (0 ≤ ε ≤ 1/2) | [Ji & Daughton 2011, arXiv:1109.0756](https://ar5iv.labs.arxiv.org/html/1109.0756) | — | — |
| 2 | Spitzer resistivity | η⊥ = 1.03×10⁻² Z lnΛ T_e^{−3/2} Ω cm, with T_e in eV; η∥ = η⊥/1.96 | [NRL Plasma Formulary](https://www.nrl.navy.mil/News-Media/Publications/NRL-Plasma-Formulary/) | — | η∥ = 5.26×10⁻⁵ Z lnΛ T_e^{−3/2} Ω m |
| 3 | Sweet–Parker | δ_SP = L/S^{1/2}. The rate v_in/v_A decays as S^{−1/2} and τ_SP = τ_A S^{1/2}. | [Bhattacharjee et al. 2009, arXiv:0906.5599](https://ar5iv.labs.arxiv.org/html/0906.5599) | Their runs at S = 3.14×10⁴, 2.51×10⁵ and 6.28×10⁵ give S^{−1/2} = 5.6×10⁻³, 2×10⁻³ and 1.26×10⁻³ | 5.64×10⁻³, 2.00×10⁻³, 1.26×10⁻³ |
| 4 | Plasmoid linear growth | γ_max ~ S^{1/4} v_A/L; plasmoid number ~ S^{3/8} | [Loureiro, Schekochihin & Cowley, Phys. Plasmas 14, 100703 (2007), arXiv:astro-ph/0703631](https://arxiv.org/abs/astro-ph/0703631) | — | For S = 10¹³: γτ_A ~ 2.2×10³ and N ~ 10⁵ |
| 5 | Same scalings from the tearing dispersion relation | γ_L,max ~ S^{1/4}, κ_max = k_max L ~ S^{3/8}. The number of plasmoids in the linear regime scales like κ_max. | [Bhattacharjee et al. 2009](https://ar5iv.labs.arxiv.org/html/0906.5599) | — | — |
| 6 | Critical Lundquist number | Plasmoid onset at S_c | [Ji & Daughton 2011](https://ar5iv.labs.arxiv.org/html/1109.0756); [Bhattacharjee et al. 2009](https://ar5iv.labs.arxiv.org/html/0906.5599) | Ji & Daughton: S_c ~ 10⁴. Bhattacharjee: "approximately 3×10⁴ in the present study". | — |
| 7 | Plasmoid-regime rate | The rate on the deepest Sweet–Parker level is about 1/√S_c | [Ji & Daughton 2011](https://ar5iv.labs.arxiv.org/html/1109.0756); [Uzdensky, Loureiro & Schekochihin, PRL 105, 235002 (2010), arXiv:1008.3330](https://arxiv.org/abs/1008.3330) | Ji & Daughton: "~ 0.01". Uzdensky: "The global reconnection rate is shown to be independent of the Lundquist number". Bhattacharjee: "weakly dependent", "larger than the Sweet-Parker rate by nearly an order of magnitude". | S_c = 10⁴ gives 0.010; S_c = 3×10⁴ gives 0.0058 |
| 8 | Plasmoid flux distribution | f(ψ) ∝ ψ⁻² | [Uzdensky et al. 2010](https://arxiv.org/abs/1008.3330) | "an inverse square law" | — |
| 9 | Collisionless transition | λ = L/d_i (or L/ρ_s). Secondary islands appear in collisionless sheets above λ_c ~ 50. | [Ji & Daughton 2011](https://ar5iv.labs.arxiv.org/html/1109.0756) | λ_c ~ 50 | — |

**Lundquist numbers by system** (from [Ji & Daughton 2011 Table 1](https://ar5iv.labs.arxiv.org/html/1109.0756); the ar5iv conversion can garble exponents, so check the PDF before treating these as final) [mainstream]

| System | L (m) | T_e (eV) | n_e (m⁻³) | B_T (T) | S (published) | λ (published) | Computed S |
|---|---|---|---|---|---|---|---|
| MRX (ε = 1/4, B_R = 0.3B_T) | 0.8 | 10 | 1×10¹⁹ | 0.1 | 3×10³ | 1.5×10² | **3.13×10³** (lnΛ = 10); Sweet–Parker rate 0.018 |
| MST | 1.0 | 1.3×10³ | 9×10¹⁸ | 0.5 | 3×10⁶ | 62 | — |
| TFTR | 0.9 | 1.3×10⁴ | 1×10²⁰ | 5.6 | 1×10⁸ | 2.3×10² | — |
| ITER | 4 | 2×10⁴ | 1×10²⁰ | 5.3 | 6×10⁸ | 5×10² | — |
| Magnetopause | 6×10⁷ | 300 | 1×10⁷ | 5×10⁻⁸ | 6×10¹³ | 9×10² | — |
| Magnetotail (B_R = B_T) | 6×10⁸ | 600 | 3×10⁵ | 2×10⁻⁸ | 4×10¹⁵ | 1.3×10³ | **2.8×10¹⁵** (lnΛ = 30) |
| Solar wind | 2×10¹⁰ | 10 | 7×10⁶ | 7×10⁻⁹ | 3×10¹² | 2×10⁵ | — |
| Solar corona (B_R = 0.1B_T) | 1×10⁷ | 200 | 1×10¹⁵ | 2×10⁻² | 1×10¹³ | 4×10⁷ | **2.3×10¹³** (lnΛ = 20) |
| Solar chromosphere | 1×10⁷ | 0.5 | 1×10¹⁷ | 2×10⁻² | 1×10⁸ | 3×10⁸ | — |
| Crab nebula flares | 1×10¹⁴ | 130 | 10⁶ | 10⁻⁷ | 5×10²⁰ | 2×10¹¹ | — |
| Interstellar medium | 5×10¹⁹ | 1 | 10⁵ | 5×10⁻¹⁰ | 2×10²⁰ | 1×10¹⁴ | — |

- Yamada, Kulsrud & Ji, "Magnetic reconnection", Rev. Mod. Phys. 82, 603 (2010), [doi:10.1103/RevModPhys.82.603](https://doi.org/10.1103/RevModPhys.82.603), is the standard review. I did not read it in this session.

### Inferences
- Computed timescales for the solar corona (L_CS = 5×10⁶ m, B_R = 2 mT): v_A = 1.38×10⁶ m s⁻¹ and τ_A = 3.6 s. τ_SP = τ_A S^{1/2} ≈ 1.75×10⁷ s, about 200 days. At a plasmoid-regime rate of 0.01, τ ≈ 360 s. Only the second is compatible with flares lasting minutes to hours, which is the standard reason S^{−1/2} Sweet–Parker cannot explain flares.
- My Spitzer-based S agrees with Ji & Daughton to within a factor of 2.5 in every row checked. The remaining spread comes from lnΛ and the parallel-vs-perpendicular η choice, so a Rust implementation should expose both as parameters.
- MRX at S ≈ 3×10³ lies below S_c ~ 10⁴, so it should show single-X-line Sweet–Parker-like behaviour (rate ~0.018). The corona, magnetotail and ISM lie far above S_c and far above λ_c ~ 50, so they sit in the multiple-X-line collisionless or plasmoid regime of the Ji–Daughton phase diagram.
- Relevance to FRC fusion: C-2 merging (Q2) and Helion merging are deliberately fast, driven reconnection events. The 2010 C-2 paper used a 2-D resistive MHD code ("LamyRidge") to reproduce the merging, per the [APS 2010 abstract](https://archive.aps.org/dpp/2010/di3/1).

### Gaps
- I did not retrieve Petschek's rate (~π/(8 ln S)) or the GEM-challenge collisionless rate (~0.1) from a primary source in this session. Ji & Daughton's wording on collisionless rates came through garbled and is not quoted.
- I did not retrieve MRX's measured confirmation of generalised Sweet–Parker scaling (Ji et al. 1998).

## Q4. Helion Energy: colliding-FRC pulsed approach, D–³He, inductive energy recovery, results and assessments to 2026

### Takeaway
Helion forms two FRCs, accelerates them into each other, merges them and compresses the result magnetically in under 1 ms. It then recovers energy inductively as the plasma pushes back on the coils. Company-reported milestones: Trenta reached 100 M°C (8.6 keV) in 2021. Polaris began operating at the end of 2024, switched to D–T in January 2026 and reported 150 M°C (12.9 keV). Helion reports a tenfold rise in neutron signal and 14.1 MeV neutron and 3.56 MeV alpha signatures. None of this is peer-reviewed, and no absolute yield or Q has been published. The ">95% round-trip" recovery figure belongs to the 2014 subscale Grande system, not to a fusion machine. Wurzel & Hsu's gain formula shows why recovery is the crux: with 95% efficiency in both directions, net electricity needs only Q_sci > 0.108, but at 80% it needs Q_sci > 0.56. In 2026 Helion removed the 2028 Orion power date from its website.

### Cited Findings

**Machine history (company-reported unless noted)** [contested: company data, not peer-reviewed]
- Venti (5th prototype, 2018): "7 T", ion temperature "2 keV", triple product "6.4 × 10¹⁸ keV·s/m³" — [Wikipedia, Helion Energy](https://en.wikipedia.org/wiki/Helion_Energy)
- Trenta (6th, 2021): "100 million degrees C" after "a 16-month test cycle with more than 10,000 pulses". Compression fields "exceeded 8 T", ion temperatures "surpassed 8 keV", electron temperatures "exceeded 1 keV", ion densities up to "3 × 10²² ions/m³", confinement "up to 0.5 ms" — [Wikipedia, Helion Energy](https://en.wikipedia.org/wiki/Helion_Energy)
- Grande (4th, 2014): "over 95% round-trip efficiency for over 1 million pulses", described as a subscale pulsed magnetic system — [Wikipedia, Helion Energy](https://en.wikipedia.org/wiki/Helion_Energy)
- Polaris (7th): "25% larger than Trenta". Pulse rate rises "from one pulse every 10 minutes to one pulse per second for short periods". Built in late 2024 with a Washington State Department of Health licence — [Wikipedia, Helion Energy](https://en.wikipedia.org/wiki/Helion_Energy). It "began operating at the end of 2024" and began using D–T "in January 2026" — [Helion newsroom, 13 Feb 2026](https://www.helionenergy.com/newsroom/helion-achieves-new-fusion-energy-milestones)
- Polaris temperature: "150MºC". A review by Dr. Ryan McBride quoted on the page cites "temperatures exceeding 13 keV or 150 million degrees Celsius". The page does not say whether this is the ion or electron temperature, and gives no field, density, yield or recovery number — [Helion newsroom](https://www.helionenergy.com/newsroom/helion-achieves-new-fusion-energy-milestones). Computed: 150 M°C = 12.93 keV, consistent with "exceeding 13 keV" after rounding.
- Polaris compression sequence: "Merged plasmas start at about 10 million to 20 million degrees" and are compressed to 150 M°C "in less than a millisecond". The long-term goal is 200 M°C, which Kirtley calls "that optimal sweet spot". Asked about scientific breakeven, Kirtley "demurred" — [TechCrunch, 13 Feb 2026](https://techcrunch.com/2026/02/13/fusion-startup-helion-hits-blistering-temps-as-it-races-toward-2028-deadline/)
- D–T measurement details: the first D–T pulse gave "a ten‑fold increase in neutron signal over typical D‑D pulses". The organic scintillators saturated. A diamond detector array saw 14.1 MeV peaks, Cu activation foils showed activity, and CR-39 recorded 3.56 MeV alphas. There were five more D–T pulses. "exact yield numbers are still being refined" — [Helion blog, How we conducted and measured D-T fusion](https://www.helionenergy.com/blog/how-we-conducted-and-measured-d-t-fusion)

**Energy recovery and gain**
- Helion's own description: expanding plasma induces current in the coils, which Helion compares to "regenerative braking". The condition it states is "If 95% of the input energy after the pulse is recovered (ηelec)". Net electricity is then "possible even when the scientific energy gain (Qsci) is below 1". FRC beta "approach[es] a ratio of 1"; each pulse lasts "less than a millisecond" — [Helion blog, How to make fusion electricity without ignition](https://www.helionenergy.com/blog/how-to-make-fusion-electricity-without-ignition) [contested]
- Helion also writes that "even at circuit efficiencies above 90%, there are still losses" and that it measures them with baseline and no-fusion pulses — [Helion blog, electrical diagnostics](https://www.helionenergy.com/blog/measuring-electricity-production-from-fusion-electrical-diagnostics) (from the search-result summary) [contested]
- Gain formula [mainstream]: Q_eng = η_elec·η_E·(Q_sci + 1) − 1 (Eq. 27). Here η_E is the efficiency from input electricity to applied heating, η_elec the efficiency from output power to electricity, and Q_sci = P_F/P_ext — [Wurzel & Hsu, arXiv:2105.10954](https://ar5iv.labs.arxiv.org/html/2105.10954)
- Computed break-even for Q_eng = 0, with η_elec = η_E = η: Q_sci > 1/η² − 1. That is **0.108** at η = 0.95, 0.235 at 0.90, 0.562 at 0.80 and 3.00 at 0.50.

**Fuel** [mainstream physics]
- Computed Q-values from AME2016 atomic masses: D + ³He → ⁴He + p, **18.353 MeV** (α 3.69 MeV, p 14.66 MeV). D + T, 17.589 MeV. D + D → T + p, 4.033 MeV. D + D → ³He + n, 3.269 MeV.
- Helion's fuel is D–³He, with "only 5%" of energy released as fast neutrons. ³He is to be bred from D–D side reactions — [Wikipedia, Helion Energy](https://en.wikipedia.org/wiki/Helion_Energy). Kirtley says ³He production has reached "very high efficiencies", without numbers — [TechCrunch](https://techcrunch.com/2026/02/13/fusion-startup-helion-hits-blistering-temps-as-it-races-toward-2028-deadline/) [contested]

**Commercial timeline**
- Microsoft PPA, 10 May 2023: "50-MWe to Microsoft starting in 2028". Nucor, Oct 2023: "500-MWe fusion power plant", targeting 2030 — [Wikipedia, Helion Energy](https://en.wikipedia.org/wiki/Helion_Energy)
- Orion construction "began in July 2025 in Malaga, Wash." — [Helion newsroom](https://www.helionenergy.com/newsroom/helion-achieves-new-fusion-energy-milestones). Wikipedia describes it as a planned 50 MW plant, "Expected to be twice the size of Polaris" — [Wikipedia](https://en.wikipedia.org/wiki/Helion_Energy)
- Sept 2026: Helion "removed a prediction from its website that Orion ... will begin producing electricity by 2028". The Polaris goal shifted from "net electricity" to "demonstrating electricity from fusion" — [Axios Pro, 8 Sep 2026](https://www.axios.com/pro/climate-deals/2026/09/08/helion-energy-fusion-net-electricity-date-year). The page returned HTTP 403, so this wording comes from the search-result summary only. [contested → timeline slipping]
- Funding: Series F of $425 M in Jan 2025 at a $5.4 B valuation, and Series G of $465 M led by Thrive Capital at $15.5 B — [Wikipedia](https://en.wikipedia.org/wiki/Helion_Energy)

**Independent assessments**
- The 2018 MITRE/JASON report for ARPA-E named the key challenge as "whether they can simultaneously achieve sufficiently high compression while maintaining plasma stability" — via [Wikipedia, Helion Energy](https://en.wikipedia.org/wiki/Helion_Energy)
- Daniel Jassby (retired PPPL) wrote in APS *Physics & Society* (Apr 2019) that such start-ups practise "voodoo fusion" and that "almost all have apparently never produced a single D–D fusion reaction" — via [Wikipedia](https://en.wikipedia.org/wiki/Helion_Energy). Helion's 2026 neutron reports, if confirmed, supersede the D–D part for Helion. [contested]
- MIT Technology Review, 2023: Paul Wilson (UW–Madison) said a 2028 plant would "surprise" him and that the key question is "whether they can create a large enough pulse to create enough energy". Jessica Lovering cited "two big unproven steps" — [MIT Technology Review](https://www.technologyreview.com/2023/05/10/1072812/this-startup-says-its-first-fusion-plant-is-five-years-away-experts-doubt-it/)
- TechCrunch's 2026 coverage contains no independent expert commentary — [TechCrunch](https://techcrunch.com/2026/02/13/fusion-startup-helion-hits-blistering-temps-as-it-races-toward-2028-deadline/). The table of ARPA-E Lawson-parameter data I retrieved has no Helion rows — [Wurzel & Hsu Table 7](https://ar5iv.labs.arxiv.org/html/2105.10954).

### Inferences
- Computed compression check, using Tuszewski's adiabatic law T ∝ B_e^{4/5}. Heating 10–20 M°C to 150 M°C by flux compression alone needs B_e to rise **12×–30×**, with density up 21×–58×. If the final field is around 10 T (Trenta reported >8 T), the merged plasmoid would start at roughly 0.3–0.8 T. This is an illustration, since real compression mixes flux and wall compression and is not purely adiabatic.
- Computed check on Trenta: the reported maxima n = 3×10²² m⁻³, Ti = 8 keV and Te = 1 keV give n(Ti + Te) = 4.3×10⁷ Pa. That needs B_e ≥ 10.4 T at β = 1, but 8 T gives only p_M = 2.5×10⁷ Pa. So the reported maxima were not simultaneous. Their product (1.2×10²⁰ keV m⁻³ s) is not a valid triple product and should not be compared with C-2W's 2.3×10¹⁷.
- D–³He at 13–17 keV gives a reactivity well below D–T at the same temperature (see the sibling reactivity notes). The 200 M°C "sweet spot" claim therefore rests on Helion's pulsed, high-β, direct-recovery economics, not on reactivity at its peak. This needs the sibling ⟨σv⟩ tables to quantify.
- The decisive unpublished number is the measured round-trip efficiency on Polaris pulses that contain plasma. Grande's 95% was a circuit without fusion plasma. The η_elec·η_E product enters the break-even Q_sci quadratically: 0.95² gives 0.108, while 0.90² gives 0.235.

### Gaps
- No peer-reviewed Polaris paper, absolute D–T neutron yield, Q_sci, field strength, density or measured recovery efficiency was found up to 24 Oct 2026.
- I could not retrieve the peer-reviewed MSNW/Helion merging paper (Slough et al., Nucl. Fusion 2011) or its numbers.
- I could not find whether the "13 keV" in the McBride review is an ion or electron temperature, or how it was measured.
- Microsoft PPA pricing and penalty terms are not public in the sources I found.

## Q5. TAE Technologies: C-2U, C-2W "Norman", Norm, Copernicus, Da Vinci, the beam-driven FRC and p–¹¹B plans

### Takeaway
TAE's beam-driven FRC programme has peer-reviewed, company-authored results. C-2 (2010) merged FRCs to Ti + Te > 0.5 keV. C-2U (2015) used 10+ MW of 15 keV neutral beams to sustain FRCs for 5+ ms. C-2W/Norman reached Te > 500 eV and Ttot > 3 keV for up to 30 ms, and the best independent-compilation triple product is 2.3×10¹⁷ keV m⁻³ s. That is about 10⁻⁴ of the D–T ignition level of about 3×10²¹. In 2025 TAE reported in Nature Communications that it formed FRCs with neutral beams alone ("Norm"). It then cancelled Copernicus (17 Nov 2025) and targets a p–¹¹B "Da Vinci" plant in the early 2030s. A $6 B all-stock merger with Trump Media & Technology Group, announced in Dec 2025, had not closed as of the last filings I found (S-4, Sept 2026).

### Cited Findings
- C-2 [mainstream]: peak Ti 400 eV, Te 150 eV, n = 1×10¹⁹ m⁻³ and 10⁹ neutrons per second for 3 ms — [Wikipedia, TAE Technologies](https://en.wikipedia.org/wiki/TAE_Technologies). Merging gave Ti + Te > 0.5 keV — [Binderbauer et al. PRL 2010](https://pubmed.ncbi.nlm.nih.gov/20867853)
- C-2U [mainstream]: NBI rose from about 4 MW to "10+ MW (15 keV hydrogen)". FRCs were "sustained for up to 5+ ms", limited "by hardware and electric supply constraints", and diamagnetism in the best shots lasted "over 11 ms" — [Gota et al., IAEA EX/P3-41 poster (TAE-hosted)](https://tae.com/wp-content/uploads/fliphtml5/9/files/basic-html/index.html); [APS DPP 2016 abstract](https://meetings-archive.aps.org/dpp/2016/cp10/62/). Wikipedia: "lasting 5 milliseconds with no sign of decay" — [Wikipedia](https://en.wikipedia.org/wiki/TAE_Technologies)
- C-2W/Norman [mainstream]: "Te > 500 eV; total electron and ion temperature, T_tot > 3 keV", sustained "for up to 30 ms", limited by NBI pulse length. Google machine-learning optimisation is credited — [Gota et al., Nucl. Fusion 2021 (PDF)](https://sun.ps.uci.edu/zlin/bib/gota21.pdf). The PDF was too large to fetch, so this wording comes from the search-result abstract. The independent compilation lists C-2W shots at Ti* = 0.6–3.5 keV, n = 1.3–2.08×10¹⁹ m⁻³ and τE* = 1.2–5 ms. The 2020–2021 rows cite a "conference talk, citable reference pending" — [Wurzel & Hsu Table 7](https://ar5iv.labs.arxiv.org/html/2105.10954)
- Norm, NBI-only formation [mainstream: peer-reviewed, company-authored]: "8 neutral beams with 15 keV energy, that can deliver up to 13 MW", of which "Typically about 8 MW" is absorbed. The seed plasma comes from an E×B discharge in 1–4 ms. Steady-state equilibrium lasts 21–30 ms, and the NBI-sustained phase "up to 40 ms", limited by power-supply storage. The dataset covers 4,671 shots in the reversal statistics and more than 50,000 Norman shots — [Roche et al., Nat. Commun. 16, 3487 (2025)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11993635); [APS DPP 2025 abstract](https://meetings-archive.aps.org/dpp/2025/gt03/5)
- Roadmap (TAE, 17 Nov 2025) [contested]: Norm "reduces reactor length and complexity by up to 50 percent". Norm is "being upgraded to achieve 100 million degrees Celsius". Copernicus is unnecessary. Commercial fusion from the hydrogen–boron design is targeted for the "early 2030s" — [TAE press release](https://tae.com/tae-shortens-device-roadmap-prepares-for-commercial-era/)
- Copernicus as planned (now cancelled): about $200 M, a target of about 100 M°C to validate D–T-equivalent conditions, and net energy "around 2025". Da Vinci targets "3 billion °C" and a prototype commercial reactor by 2030 — [Wikipedia, TAE Technologies](https://en.wikipedia.org/wiki/TAE_Technologies) [contested]. Computed: 3×10⁹ °C = 258.5 keV, and 100 M°C = 8.62 keV.
- p–¹¹B [mainstream physics]: computed Q = 8.682 MeV (p + ¹¹B → 3 ⁴He). TAE collaborators made the "first measurement of p-¹¹B fusion in magnetically confined plasma" at LHD in 2023 — [Wikipedia](https://en.wikipedia.org/wiki/TAE_Technologies)
- Theoretical gain for a 100 MW FRC [contested]: Q of 35 for D–T, 3 for D–³He, 2.7 for p–¹¹B and 4.3 for polarised p–¹¹B — [Wikipedia, TAE Technologies](https://en.wikipedia.org/wiki/TAE_Technologies)
- TMTG merger: announced 18 Dec 2025, all-stock, valued above $6 B, with each side owning about 50%. Devin Nunes and Michl Binderbauer would be co-CEOs. The parties planned to site and start a 50-MWe plant in 2026 — [ANS Nuclear Newswire](https://www.ans.org/news/2025-12-19/article-7632/trump-media-to-merge-with-fusion-startup-tae-technologies-in-6b-deal/). A 10 June 2026 update set closing for "the fourth quarter of 2026 or sooner" and dropped the Truth Social spin-off — [SEC 8-K exhibit](https://www.sec.gov/Archives/edgar/data/0001849635/000114036126024735/ef20075853_ex99-1.htm) (from the search summary). A Form S-4 was filed in Sept 2026 — [SEC S-4](https://www.sec.gov/Archives/edgar/data/0001849635/000143774926031534/djt20260922_s4.htm). Sources disagree on the cash: ANS says $200 M "at signing" plus $100 M on S-4 filing, while the SEC-filing summary says "$200 million in cash at closing".

### Inferences
- Computed: C-2W shot 118340 gives n·T·τ = 2.275×10¹⁷ keV m⁻³ s, matching the published 2.3×10¹⁷. That is 7.6×10⁻⁵ of the commonly quoted D–T ignition value of about 3×10²¹ ([Lawson criterion](https://en.wikipedia.org/wiki/Lawson_criterion)). p–¹¹B needs far more.
- The beam-driven FRC's main physics advance is sustainment, from 5 ms (C-2U) to 30–40 ms (C-2W/Norm), limited by power supplies rather than instabilities. The open question is energy confinement scaling with Te, which TAE plots but whose values I did not retrieve.
- Wikipedia's estimated Q of 2.7 for a p–¹¹B FRC, below 35 for D–T, shows that aneutronic operation trades gain for low neutron output. Bremsstrahlung at about 250 keV is the standard objection; see the sibling aneutronic notes.

### Gaps
- I did not retrieve the numeric τE–Te scaling from Gota et al. 2017, 2019 or 2021, nor Norm's field, density and Te values (the Roche paper does not state them).
- I could not confirm the merger's closing or vote results after the Sept 2026 S-4.
- I found no independent (non-TAE) assessment of the Da Vinci timeline.

## Q6. Other plasmoid claims: ball lightning as plasmoids, lab analogues, and contested plasmoid-energy claims

### Takeaway
A self-confined magnetised-plasma ball lightning that stores the reported energies is ruled out by the virial theorem (Finkelstein & Rubinstein 1964). Atmospheric pressure caps the magnetic plus thermal energy at 3p₀V, about 0.30 MJ m⁻³ or 1.3 kJ for a 20 cm ball. Estimates from reports are 1–10 MJ m⁻³, or 10–100 kJ. Chemical stores fit easily: oxidising about 120 g of silicon per m³ supplies 3.9 MJ m⁻³. The best laboratory analogue is the silicon-arc ball of Paiva et al. (PRL 2007): 1–4 cm, up to 8 s. The 2012 Qinghai spectrum (Cen et al., PRL 2014) showed Si, Fe and Ca lines from soil, supporting the silicon model. SAFIRE/Aureon claims of over-unity plasma output and transmutation have no peer-reviewed data.

### Cited Findings
- Natural ball lightning statistics [mainstream observation]: diameter "1–100 cm ..., most commonly 10–20 cm", lifetime "from one second to over a minute" (Charman 1972) — [Wikipedia, Ball lightning](https://en.wikipedia.org/wiki/Ball_lightning)
- Qinghai 2012 event (Cen, Yuan & Xue, PRL 2014) [mainstream, single event]: distance "900 m", video "1.64 seconds", diameter "5 m", travel "about 15 m" at "8.6 m/s". Neutral Si, Ca, Fe, N and O lines were seen, with light and O/N emission oscillating at about 100 Hz — [Wikipedia, Ball lightning](https://en.wikipedia.org/wiki/Ball_lightning)
- Paiva et al., PRL 98, 048501 (2007) [mainstream lab analogue]: arcs in pure Si produced luminous balls "with lifetimes on the order of seconds", demonstrating "the role of vaporization and oxidation of Si, as proposed by the Abrahamson-Dinniss theory" — [ADS abstract](https://ui.adsabs.harvard.edu/abs/2007PhRvL..98d8501P). Physics Today reported the balls as 1–4 cm across, with lifetimes up to 8 s and speeds of 5–30 cm/s, made at room temperature, atmospheric pressure and 70% relative humidity — [Physics Today](https://physicstoday.aip.org/news/ball-lightning-in-the-lab)
- Follow-on estimate by the same group: "mean energy density of 3.9 MJ m⁻³" — [SBPMat/ICAM 2009 abstract](https://www.sbpmat.org.br/icam2009dir/submission/palestrante/arquivos/PBB4.pdf). I saw only the search snippet, and the file attribution is uncertain. [contested]
- Other lab analogues [mainstream]: a plasma toroid lasting "about 200 milliseconds in a partial atmosphere". Atmospheric plasmoids over water last "up to only 350 microseconds" in the UV–NIR. Ohtsuki & Ofuruton (Nature 1991) made microwave "plasma fireballs" with "a 2.45 GHz, 5 kW (maximum power)" source. Jerby & Dikhtyar's microwave balls contain nanoparticles of "an average radius of 25 nm" — [Wikipedia, Ball lightning](https://en.wikipedia.org/wiki/Ball_lightning)
- Energy estimates [contested]: eyewitness-based totals are "in the order of 10 or 100 kJ", with a few reports up to 10 MJ. Stakhanov estimated about 1–10 J cm⁻³, and Bychkov et al. more than 1 kJ cm⁻³ for the most energetic cases. These come from search summaries of [arXiv:1102.2829](https://ar5iv.arxiv.org/html/1102.2829) and [arXiv:1007.4377](https://arxiv.org/pdf/1007.4377v52), and I could not tell which figure comes from which paper. Fedosin's electron–ion model gives "maximum diameter of 34 cm", charge "about 10 microcoulombs" and energy "reaching 11 kilojoules" — [Wikipedia](https://en.wikipedia.org/wiki/Ball_lightning)
- Virial theorem [mainstream; refutes the energetic self-confined plasmoid]: "confinement by self-field alone is inconsistent with conservation laws". Including air pressure "gives an upper bound to the stored energy", and "This upper bound is much less than the energies reported for some occurrences" — [Finkelstein & Rubinstein, Phys. Rev. 135, A390 (1964)](https://link.aps.org/doi/10.1103/PhysRev.135.A390). A 1981 Nature comment argues the virial objection stands for continuum plasmoid models — [Nature 290, 166](https://www.nature.com/articles/290166a0.pdf). A minority view says the bound does not apply to an open, field-fed ball — [arXiv:1007.4377](https://arxiv.org/pdf/1007.4377v52). [contested minority]
- SAFIRE/Aureon [contested; no data]: a secondhand forum claim says the "plasma reactor can produce about 14 times more energy than is needed to operate" it — [New Mars forum](https://newmars.com/forums/viewtopic.php?pid=218670). Aureon's newsletter says it "successfully transmuted thorium into many stable non-radioactive daughter elements", validated by an unnamed third-party lab — [Electrogenesis Substack](https://electrogenesis.substack.com/p/aureon-energy-the-safire-project). I found no peer-reviewed paper and no replication.

### Inferences
- Computed virial caps with p₀ = 101325 Pa: 3p₀ = 0.304 MJ m⁻³. A 4 cm ball can hold 10.2 J, a 20 cm ball 1.27 kJ and a 5 m ball 19.9 MJ. The largest field that 1 atm can confine is 0.505 T, while storing 3.9 MJ m⁻³ magnetically would take 3.13 T. A 20 cm ball at 1–10 MJ m⁻³ holds 4.2–42 kJ, 3–33 times the cap. So "BL is a self-confined magnetised plasmoid" is **[refuted]** for typical reported energies. "BL is a chemically fuelled aerosol, possibly with a short-lived plasma core" is **[mainstream]**, though still not proven.
- Computed chemical store: Si + O₂ → SiO₂ releases 911 kJ mol⁻¹ (standard [enthalpy of formation of SiO₂](https://en.wikipedia.org/wiki/Silicon_dioxide)). That is 7.55×10¹⁰ J per m³ of solid Si, so 3.9 MJ m⁻³ needs a Si volume fraction of only 5.2×10⁻⁵ (120 g m⁻³). The Paiva balls' 3.9 MJ m⁻³ exceeds the virial cap by about 13×, which is consistent with stored chemical energy, since the virial bound limits only field and plasma-thermal energy.
- Test for any over-unity plasmoid claim such as SAFIRE: ask for closed-boundary calorimetry, meaning total electrical input integrated over the run against heat out, with the error stated. Without that, a ratio like "14×" cannot be tested. I found no such data.

### Gaps
- I did not read the full Cen et al. 2014 PRL; its numbers here come via Wikipedia. Its exact temperature estimate was not retrieved.
- The Gatchina water-discharge plasmoid studies (Egorov & Stepanov; Versteegh et al.) were not found in the sources reached.
- I found no peer-reviewed SAFIRE/Aureon paper. I did not research other contested plasmoid claims (Shoulders' "EVOs", Chukanov) or Brilliant Light Power.
