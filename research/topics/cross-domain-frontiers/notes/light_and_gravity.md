# Light in gravitational fields, and gravity measured with light (2022–2026)

Every number marked "recomputed" was reproduced by `light_and_gravity_checks.py` in this folder
(standard library only; `python3 -I light_and_gravity_checks.py`, 0 failures). Constants are the
ones in `crates/metrology` plus IAU 2015 nominal R☉, R⊕, GM_J, R_J. Items older than 2022 are marked
**[foundational]** and are included only because a 2022–2026 result needs them. Status labels:
established; recent and replicated; recent and unreplicated; contested.

What the crates already have (read from the four READMEs and the source):
`gravitational::propagation::point_lens_geometric` (point-lens image magnifications and delay),
`lens_parameter` and `wave_optics_peak` (wave optics); `gravitics::spacetime::light_bending`
(point-mass deflection to second order, GR or Puthoff's vacuum, no PPN γ);
`gravitics::spacetime::horizon_temperature`, `unruh_temperature`, `hawking_temperature`;
`gravitics::gem` and `gravitics::superconductors::ring_frame_dragging` (frame dragging of matter,
not light); `gravitational::detectors::clock_strain` (clocks as wave detectors, not geodesy);
`gravitics::anomalies` (LLR used for the shielding bound only). Nothing in `geopotential` handles
the centrifugal potential, W₀, the normal potential U₀ or clock readings; nothing anywhere computes
cosmological distances, Shapiro delays, shadows, PPN γ/β/η or microlensing light curves.

## 1. Gravitational lensing beyond the point lens

### Takeaway
Four new lensing pieces have clean published tests: the singular isothermal sphere (Euclid's 2025
NGC 6505 ring reproduces its measured velocity dispersion to 0.6σ), NFW convergence, microlensing
light curves and masses (the isolated black hole OGLE-2011-BLG-0462 at 7.15 M☉ is reproduced to
0.1%), and time-delay distances (the six H0LiCOW lenses give back H₀ = 73.3 to 0.1%). All four need
one new shared piece, a flat-ΛCDM distance integral. TDCOSMO 2025 has moved the lensing H₀ down to
71.6 (+3.9/−3.3); the tension with Planck is no longer significant from lensing alone.

### Cited Findings
- **L0. Cosmological distances (prerequisite) [established; foundational].**
  D_M(z) = (c/H₀)∫₀ᶻ dz′/√(Ωm(1+z′)³ + 1 − Ωm) in flat ΛCDM, D_A = D_M/(1+z),
  D_ls = (D_M(z_s) − D_M(z_l))/(1+z_s). Test: through L3 (H0LiCOW) and L1 (NGC 6505) below. Crate:
  new `gravitational::cosmology` (the crate already defines the h = 1 Hubble rate in `lib.rs`);
  uses `metrology::{C, MPC}`. Links gravitational ↔ metrology. Source of the D_Δt definition, eq. 2:
  — [Wong et al. 2020, H0LiCOW XIII, arXiv:1907.04869](https://arxiv.org/abs/1907.04869)
- **L1. Singular isothermal sphere [established; foundational formula, 2025 test].**
  θ_E = 4π(σ/c)² D_ls/D_s; mass inside the ring M(<θ_E) = (c²/4G)(D_l D_s/D_ls)θ_E²;
  Σ_cr = c² D_s/(4πG D_l D_ls); two-image delay Δt = (D_Δt/2c)(θ_A² − θ_B²). Formula source:
  [Narayan & Bartelmann 1996, arXiv:astro-ph/9606001](https://arxiv.org/abs/astro-ph/9606001).
  Published test: Euclid's first lens, NGC 6505: z_l = 0.042, z_s = 0.406, θ_E = 2.5″ = 2.1 kpc,
  σ* = 303 ± 15 km/s, dark-matter fraction inside θ_E (11.1 +5.4/−3.5)% — [O'Riordan et al. 2025,
  A&A 694, A145, arXiv:2502.06505](https://arxiv.org/abs/2502.06505). Recomputed: SIS σ from θ_E
  = 313 km/s (0.64σ from the measured 303; this ratio is independent of H₀); θ_E = 2.07 kpc at
  H₀ = 70, Ωm = 0.3 (paper 2.1); M(<θ_E) = 1.48×10¹¹ M☉ (scales as 1/H₀; the paper's own value not
  found). Crate: new `gravitational::lensing`. Status: recent and unreplicated (single 2025 lens
  model; the ring itself is unambiguous).
- **L2. NFW lens convergence [established; foundational].** Σ(x) = 2r_s δ_c ρ_c f(x) with
  f(x) = [1 − (2/√(1−x²)) artanh√((1−x)/(1+x))]/(x²−1) for x < 1, the same with
  (2/√(x²−1)) arctan√((x−1)/(1+x)) for x > 1, and f(1) = 1/3; x = R/r_s. Test: continuity at x = 1
  and slope −0.4 there (recomputed). Implementation note found while checking: both closed forms
  cancel catastrophically near x = 1 (relative error ≈ 10⁻¹⁶/|x−1|², 3×10⁻⁵ at |x−1| = 10⁻⁶), so a
  port needs a series for |x−1| < 10⁻³ — [Wright & Brainerd 2000, ApJ 534, 34,
  arXiv:astro-ph/9908213](https://arxiv.org/abs/astro-ph/9908213). Crate: `gravitational::lensing`;
  uses the NFW halo, which could also feed `gravitics::modified` (rotation curves). Links
  gravitational ↔ gravitics.
- **L3. Time-delay distance and H₀ [contested (H₀ tension)].** Fermat potential
  τ(θ) = (D_Δt/c)[(θ−β)²/2 − ψ(θ)], D_Δt ≡ (1+z_d)D_d D_s/D_ds, so at fixed Ωm, H₀ ∝ 1/D_Δt —
  [Wong et al. 2020](https://arxiv.org/abs/1907.04869). Published inputs (Table 2, Mpc):
  B1608+656 (z 0.6304→1.394) 5156 +296/−236; RXJ1131−1231 (0.295→0.654) 2096 +98/−83; HE0435−1223
  (0.4546→1.693) 2707 +183/−168; SDSS1206+4332 (0.745→1.789) 5769 +589/−471; WFI2033−4723
  (0.6575→1.662) 4784 +399/−248; PG1115+080 (0.311→1.722) 1470 +137/−127; combined flat-ΛCDM
  H₀ = 73.3 +1.7/−1.8 km/s/Mpc — [Wong et al. 2020](https://arxiv.org/abs/1907.04869). Recomputed:
  per-lens H₀ = 69.7, 78.0, 71.3, 69.8, 71.4, 79.4; inverse-variance mean 73.28 ± 1.8 (D_Δt only,
  Ωm = 0.3, errors symmetrised). Crate: `gravitational::lensing` + `gravitational::cosmology`.
- **L3 newer values.** TDCOSMO 2025 (eight lensed quasars, new JWST/Keck/VLT stellar kinematics,
  blinded): H₀ = 71.6 +3.9/−3.3 km/s/Mpc with Pantheon+ Ωm; adding SLACS and SL2S gives 4.6%
  precision — [TDCOSMO 2025, A&A 704, A63, arXiv:2506.03023](https://arxiv.org/abs/2506.03023).
  Status: recent and unreplicated.
- **L3 lensed supernovae (same D_Δt algorithm, cluster lenses).** SN Refsdal: H₀ = 64.8 +4.4/−4.3
  (eight cluster models), 66.6 +4.1/−3.3 (two best models) — [Kelly et al. 2023, Science,
  via Carnegie](https://carnegiescience.edu/constraints-hubble-constant-supernova-refsdals-reappearance);
  delay SX–S1 = 376.0 +5.6/−5.5 days — [Kelly et al. 2023, ApJ, "The Magnificent Five Images"](https://doaj.org/article/a3c3b6986c3b4987a332829f14dac143).
  SN H0pe (z = 1.78 behind PLCK G165.7+67.0 at z = 0.35): H₀ = 71.8 +9.2/−8.1 from delays alone,
  75.7 +8.1/−5.5 using magnifications — [Pascale et al. 2025, ApJ, arXiv:2403.18902](https://arxiv.org/abs/2403.18902).
  Status: recent and unreplicated (one supernova each).
- **L4. Microlensing of real events [Paczyński formula established; the event recent and
  replicated].** u(t) = √(u₀² + ((t−t₀)/t_E)²); A(u) = (u²+2)/(u√(u²+4)) [this total is the sum of
  |μ±| already in `point_lens_geometric`; the time series is new]; astrometric shift
  δθ = θ_E u/(u²+2), maximal θ_E/(2√2) at u = √2; θ_E² = κMπ_rel, κ = 4G/(c²·AU) = 8.144 mas/M☉,
  M = θ_E/(κπ_E), π_rel = θ_Eπ_E, t_E = θ_E/μ_rel. Published (OGLE-2011-BLG-0462, HST astrometry over
  11 yr): M = 7.15 ± 0.83 M☉, D_L = 1.52 ± 0.15 kpc, D_S = 7.1 ± 1.7 kpc, θ_E = 5.496 ± 0.310 mas,
  π_E = 0.0944 ± 0.0124, t_E = 273.1 ± 6.3 d, u₀ = 0.00271, μ_rel = (−2.907, +6.751) mas/yr —
  [Sahu et al. 2025, arXiv:2503.07820](https://arxiv.org/html/2503.07820v2); first analysis
  M = 7.1 ± 1.3 M☉ at 1.58 ± 0.18 kpc — [Sahu et al. 2022, arXiv:2201.13296](https://arxiv.org/abs/2201.13296).
  Recomputed: κ = 8.1439; M = 7.149 M☉; D_L = 1.516 kpc; t_E = 273.1 d; peak A = 369; maximum
  astrometric shift 1.94 mas. Crate: `gravitational::lensing` (reuses `point_lens_geometric`).
- **L5. Lensing by a body with a spherical-harmonic field [established in the Born approximation;
  the cross-crate algorithm is a proposal].** α = ((1+γ)/c²)∫∇⊥U dl along the unperturbed ray, and
  the Shapiro delay Δt = ((1+γ)/c³)∫U dl; for a point mass α = 2(1+γ)GM/(c²b). For a ray in the
  equator of an oblate body with its axis perpendicular to the ray, U = GM/r + GMJ₂R²/(2r³) gives
  α = (4GM/c²b)(1 + J₂R²/b²). Recomputed: a numerical Born integral matches this to 10⁻⁶; Jupiter's
  J₂ (Juno, J₂ = 14 696.5×10⁻⁶ at R = 71 492 km; corrected 2026-10-25: Iess et al. give 14 696.572, and 14 696.5063 is Durante 2020 without its permanent tide — [Iess et al. 2018, Nature 555, 220](https://doi.org/10.1038/nature25776))
  adds 239 µas at the limb to the 16.27 mas monopole; Earth's J₂ adds 0.62 µas to 0.574 mas.
  Full-precision Gaia data were forecast to reach γ to ~10⁻⁶ using Jupiter — [Crosta & Mignard 2006,
  arXiv:astro-ph/0512359](https://arxiv.org/abs/astro-ph/0512359). Crate: `gravitics::spacetime`
  (or a new `gravitics::ppn`) taking any `geopotential::Coefficients`, integrating
  `geopotential` partials along the ray. **Links geopotential ↔ gravitics ↔ metrology.**

### Inferences
- One `cosmology` module unlocks L1, L3 and the lensed-supernova numbers; it is the single
  highest-leverage addition for this question.
- L5 is the cleanest "lensing that uses the gravity field" item: the same `.gfc` coefficients the
  crate already reads (Earth, or a Juno Jupiter model) would give µas-level deflections and Shapiro
  delays, and the J₂ term has an exact closed form to test against.
- The microlensing constant κ is exactly the solar deflection at b = 1 AU expressed per solar mass
  (4GM☉/(c²·AU) = 8.144 mas, twice the 4.07 mas at 90° elongation in P1). That identity could join
  the `physica` analogies table as a doc test.
- The lensing H₀ values now straddle the Planck–SH0ES gap (64.8 to 75.7), so the crate should carry
  them in a registry with status labels, as `gravitational::hf` does for limits.

### Gaps
- The NGC 6505 paper's own mass inside θ_E and its assumed cosmology were not found in the abstract;
  only the H₀-independent σ check is a strict test.
- An independent 2022 analysis of OGLE-2011-BLG-0462 (Lam et al.) reportedly found a lower mass;
  I did not fetch it, so the "replicated" label rests on the two Sahu-team analyses.
- Per-lens TDCOSMO 2025 D_Δt values were not extracted; the L3 test uses H0LiCOW 2020.

## 2. Light deflection and Shapiro delay as PPN tests

### Takeaway
Cassini's γ − 1 = (2.1 ± 2.3)×10⁻⁵ is still the best single measurement; the 2024 Living Review
combines planetary ephemerides to γ − 1 = (0.55 ± 1.35)×10⁻⁵ and β − 1 = (−0.45 ± 1.75)×10⁻⁵. No
Gaia light-deflection measurement of γ has been published. The 2021 Double Pulsar timing gives a
strong-field Shapiro test: range r = 6.162(21) µs against T☉m_B = 6.151 µs (0.5σ).

### Cited Findings
- **P1. Deflection with PPN γ and elongation [established; foundational].**
  δθ = ((1+γ)/2)(4GM/(c²d))((1+cos χ)/2), d = r sin χ, i.e. (1+γ)(GM/(c²r))cot(χ/2). Published
  checks: "at 90° from the Sun, the deflection is still 4 milliarcseconds"; VLBI 2009:
  γ − 1 = (−1.6 ± 1.5)×10⁻⁴; VLBA 2009 (3C273/3C279): (−2 ± 3)×10⁻⁴; Jupiter VLBI 1991: ~300 µas
  seen to ~50% — [Will 2014, Living Rev. Relativ. 17, 4, arXiv:1403.7377](https://arxiv.org/abs/1403.7377).
  Recomputed: solar limb 1.7512″; 90° elongation 4.072 mas; Jupiter limb 16.27 mas; Earth limb
  0.574 mas. Already covered in part: `gravitics::spacetime::light_bending` (point mass, b given,
  γ = 1). New: γ and the observer geometry. Crate: `gravitics::spacetime` or new `gravitics::ppn`.
- **P2. Shapiro delay and Cassini [established].** Round trip
  δt = 2(1+γ)(GM/c³)ln[(r⊕ + x⊕·n)(r_e − x_e·n)/d²], near the Sun
  δt ≈ ½(1+γ)[240 − 20 ln(d²/r)] µs (d in R☉, r in AU); Cassini 2002: b = 1.6 R☉, spacecraft at
  8.43 AU, X and Ka band, γ − 1 = (2.1 ± 2.3)×10⁻⁵ — [Will 2014, eqs. 62–63](https://arxiv.org/abs/1403.7377);
  primary: [Bertotti, Iess & Tortora 2003, Nature 425, 374](https://doi.org/10.1038/nature01997).
  Recomputed: exact log form 262.4 µs vs Will's approximation 263.8 µs; the Doppler observable
  y = 2(1+γ)(GM/c³)(db/dt)/b is ≈ 5×10⁻¹⁰ at conjunction (with db/dt taken as Earth's 29.8 km/s;
  order of magnitude only). Crate: new `gravitics::ppn` (uses `metrology::GM_SUN`). The same
  integral is L5's ∫U dl, so P2 and L5 share code.
- **P3. Planetary-ephemeris PPN bounds, 2024 [recent and replicated (several ephemerides)].**
  (β − 1) = (−0.45 ± 1.75)×10⁻⁵ and (γ − 1) = (0.55 ± 1.35)×10⁻⁵ from the common overlap of
  published ephemeris intervals; table entries include INPOP20a PF β−1 = (1.9 ± 6.28)×10⁻⁵,
  γ−1 = (2.64 ± 3.44)×10⁻⁵, J₂☉ = (2.165 ± 0.12)×10⁻⁷, and DE440 J₂☉ = (2.25 ± 0.09)×10⁻⁷ —
  [Fienga & Minazzoli 2024, Living Rev. Relativ. 27, 1, arXiv:2303.01821](https://arxiv.org/abs/2303.01821).
  Note: the combined value is an overlap of intervals, not one fit. Crate: registry in
  `gravitics::ppn`.
- **P4. Double Pulsar Shapiro delay [recent and unreplicated (one system)].**
  Δ_S = −2r ln Λ_u, Λ_u = 1 − e_T cos u − s[sin ω(cos u − e_T) + (1−e_T²)^½ cos ω sin u],
  r = Gm_B/c³ = T☉m_B, s = sin i, z_s ≡ −ln(1 − s); with lensing and retardation corrections
  −2r ln(Λ_u + δΛ_len + δΛ_ret). Published: z_s = 9.65(15), s = 0.999936, r = 6.162(21) µs,
  m_A = 1.338185 M☉, m_B = 1.248868 M☉, ω̇ = 16.899323(13) °/yr, Ṗ_b = −1.247920(78)×10⁻¹², GR's
  quadrupole formula confirmed to 1.3×10⁻⁴ (95%), and aberrational light bending detected —
  [Kramer et al. 2021, PRX 11, 041050, arXiv:2112.06795](https://arxiv.org/abs/2112.06795).
  Recomputed: T☉m_B = 6.1513 µs (0.51σ); s = 0.999936; peak delay ≈ 2r z_s = 119 µs. Crate:
  `gravitational::energy` (already has binary-pulsar decay; Shapiro r, s are the missing
  post-Keplerian pair). Links gravitational ↔ gravitics::ppn.
- **P5. Gaia [no result].** The only 2022 Gaia-related γ paper bounds the parallax-zero-point bias
  from γ's uncertainty at < 0.2 µas; it is not a γ measurement — [arXiv:2204.02789](https://arxiv.org/abs/2204.02789).

### Inferences
- P1, P2 and L5 are one integral (∫U dl and ∫∇⊥U dl) with a γ factor; a single PPN ray module in
  `gravitics` would serve the Sun (point mass), Jupiter and Earth (spherical harmonics from
  `geopotential`) and binary pulsars.
- Cassini remains the anchor; ephemeris combinations improve on it only by pooling, so the crate's
  registry should flag P3 as "combined".

### Gaps
- No 2022–2026 VLBI γ result was found; the latest VLBI numbers here are 2009–2011 via Will 2014.
- No BepiColombo superior-conjunction (SORT) γ result was found.
- The exact peak Doppler shift Bertotti et al. report was not extracted.

## 3. Optical clocks as gravity sensors (chronometric geodesy)

### Takeaway
A clock reading is a geopotential-number difference: Δν/ν = (C_B − C_A)/c² with C = W₀ − W, so
1 m²/s² = 1.113×10⁻¹⁷ and 1 cm ≈ 1.09×10⁻¹⁸. The 2022–2024 results span 1 mm (JILA, 0.5σ from g/c²),
1 cm (Wisconsin, 0.6σ), 450 m (Skytree) and 457 km of fibre (0.85σ from levelling). To compare with
the crate, `geopotential` must add the centrifugal potential (1.2×10⁻¹² of a clock rate at the
equator), the normal potential U₀ and the conventional W₀; all three have exact published tests.

### Cited Findings
- **C1. Clock reading ↔ geopotential number [established].** Δν/ν = (1+α)ΔU/c² (redshift
  parameter α = 0 in GR) — [Delva et al. 2018, eq. 1, arXiv:1812.03711](https://arxiv.org/abs/1812.03711).
  Geodetic form: Δν/ν = (C_B − C_A)/c², C(P) = W₀ − W(P), W = V + ½ω²(x² + y²), with the IHRS
  conventional W₀ = 62 636 853.4 m²/s² (IAG Resolution 1, 2015), replacing the IERS 2010 value
  62 636 856.0, a change of "up to −2.6 m²/s² (about 27 cm)" — [DGFI-TUM on IAG 2015](https://www.dgfi.tum.de/en/news-archive/na_w0/);
  [Sánchez et al. 2016, J. Geod. 90, 815](https://doi.org/10.1007/s00190-016-0913-x).
  L_G = W₀/c² = 6.969290134×10⁻¹⁰ (IAU 2000 B1.9, fixed from the 1998 W₀ = 62 636 856.0) —
  [IERS Conventions 2010, TN36](https://www.iers.org/TN36). Recomputed: L_G to 10⁻⁹; the W₀ change
  = 26.5 cm. Crate: new `geopotential::chronometric` (`clock_ratio(c_a, c_b)`,
  `geopotential_number(w, w0)`, `centrifugal_potential(r, lat)`); constants `OMEGA_EARTH`,
  `W0_IHRS`, `L_G` in `metrology`. **Links geopotential ↔ metrology ↔ gravitational::detectors.**
- **C2. GRS80 normal potential [established; foundational].** U₀ = (GM/E)arctan(e′) + ω²a²/3,
  E = √(a² − b²), e′ = E/b; GRS80: a = 6 378 137 m, GM = 3.986005×10¹⁴, ω = 7.292115×10⁻⁵ rad/s,
  b = 6 356 752.3141 m, U₀ = 62 636 860.850 m²/s² — [Moritz 2000, J. Geod. 74, 128](https://doi.org/10.1007/s001900050278).
  Recomputed to 0.01 m²/s²; U₀ − W₀(IHRS) = 7.45 m²/s² (≈ 76 cm); equatorial centrifugal potential
  108 159.5 m²/s². Crate: `geopotential::chronometric` (a clock compares against W, not the V the
  crate synthesises today).
- **C3. JILA millimetre redshift [recent and replicated].** Linear frequency gradient across one
  ~100 000-atom ⁸⁷Sr sample: −9.8(2.3)×10⁻²⁰/mm after systematics (raw −1.00(12)×10⁻¹⁹/mm, 14 runs
  over 10 days), predicted −1.09×10⁻¹⁹/mm; two-region synchronous comparison −1.28(27)×10⁻¹⁹/mm;
  instability 4.4×10⁻¹⁸/√τ, 7.6×10⁻²¹ over 92 h — [Bothwell et al. 2022, Nature 602, 420,
  arXiv:2109.12238](https://arxiv.org/abs/2109.12238) (numbers from the arXiv text; DOI
  10.1038/s41586-021-04349-7). Recomputed: g/c² = 1.090×10⁻¹⁹/mm at g = 9.796; measured is 0.48σ
  from it.
- **C4. Wisconsin centimetre redshift (independent replication) [recent and replicated].**
  Five ensembles over 1 cm, blinded: [−12.4 ± 0.7(stat) ± 2.5(sys)]×10⁻¹⁹/cm vs expected
  −10.9×10⁻¹⁹/cm — [Zheng et al. 2023, Nat. Commun. 14, arXiv:2207.07145](https://arxiv.org/abs/2207.07145).
  Recomputed: 0.58σ.
- **C5. Tokyo Skytree [foundational, 2020; restated in 2025].** Δν ≈ 21.18 Hz on ν ≈ 429.228 THz;
  ḡΔh/c² = (49 337.1 ± 1.4)×10⁻¹⁸; laser-ranged Δh = 452.596 ± 0.013 m; GNSS levelling
  452.650 ± 0.039 m; α = (1.4 ± 9.1)×10⁻⁵; RIKEN–UTokyo (15 km, 30 km fibre):
  Δν/ν = (−1652.9 ± 5.9)×10⁻¹⁸, α = (2.9 ± 3.6)×10⁻³ — [Shinkai, Takamoto & Katori 2025,
  arXiv:2502.06104](https://arxiv.org/html/2502.06104v1). Recomputed: 21.18/ν_Sr = 4.934×10⁻¹⁴
  (within 1.5×10⁻⁴ of the quoted term); implied ḡ = 9.797 m/s².
- **C6. 457 km fibre chronometric levelling [recent and unreplicated].** Two Sr lattice clocks:
  ΔW = 3918.1(2.6) m²/s² vs geodetic 3915.88(0.30) — [PR Applied 21, L061001 (2024), arXiv:2309.14953](https://arxiv.org/abs/2309.14953)
  (values as recorded in this repository's `research/topics/gravimetry-frontiers/report.md`).
  Recomputed: Δν/ν = 4.3595×10⁻¹⁴; agreement 0.85σ. This is the direct test for
  `geopotential::chronometric` against a `.gfc` model plus C2.
- **C7. Networks and space [recent].** Ten clocks in six countries compared for 45 days, 38 ratios
  — [Lindvall et al. 2025, Optica 12, 843, arXiv:2505.06763](https://arxiv.org/abs/2505.06763).
  ACES is on the ISS (launched 21 April 2025); no flight redshift result was found; a 2026 study
  is simulation only — [ESA ACES](https://www.esa.int/Science_Exploration/Human_and_Robotic_Exploration/ACES_Atomic_Clock_Ensemble_in_Space);
  [ACES overview, arXiv:2411.02912](https://arxiv.org/abs/2411.02912).

### Inferences
- The minimal, fully testable chain is: `.gfc` model → V (exists) → + centrifugal → W → C = W₀ − W →
  Δν/ν = ΔC/c² → compare with C6's 3918.1 m²/s² and C5's 49 337.1×10⁻¹⁸. Every link has a published
  number.
- `gravimetry::sim` could carry clocks on simulated orbits (see R1) to give proper time as a
  second observable next to ranging, joining gravimetry, geopotential and metrology.

### Gaps
- The 457 km paper's site coordinates were not extracted, so C6 cannot yet be recomputed from a
  geopotential model, only from its published ΔW.
- The published Nature version of Bothwell et al. may differ slightly from the arXiv numbers used.
- No published ACES in-orbit redshift value as of this search.

## 4. Gravitational redshift tests with eccentric Galileo satellites and later

### Takeaway
The 2018 Galileo tests remain the best: α = (0.19 ± 2.48)×10⁻⁵, 5.6× better than Gravity Probe A;
ground clocks (Skytree, α = (1.4 ± 9.1)×10⁻⁵) are close behind; ACES has no result yet. Redshift of
light from the Sun (633.1 m/s expected, 638 ± 6 measured) and from S2 at Sgr A* (≈ 200 km/s,
f = 0.90 ± 0.09 ± 0.15) both recompute exactly from GM, r and v.

### Cited Findings
- **R1. Eccentric-orbit clock redshift [established; two independent groups].**
  τ = ∫[1 − v²/(2c²) − (U_E + U_T)/c²]dt with U_E = GM/r[1 + J₂R²(1 − 3cos²θ)/(2r²)] and lunisolar
  tides; e = 0.162; 1008 days; "a linear and a periodic term, which amounts to 400 ns peak-to-peak";
  α = (+0.19 ± 2.48)×10⁻⁵ — [Delva et al. 2018, PRL 121, 231101, arXiv:1812.03711](https://arxiv.org/abs/1812.03711);
  companion result with uncertainty reduced more than 4× vs GP-A — [Herrmann et al. 2018, PRL 121,
  231102](https://journals.aps.org/prl/issues/121/23). Kepler form of the periodic gravitational
  part: Δτ = (√(GMa)/c²)e sin E; the total GNSS clock correction is twice that,
  F e√a sin E with F = −2√μ/c² = −4.442807633×10⁻¹⁰ s/m^½ (IS-GPS-200 [foundational]).
  Recomputed: with a = 27 979 km (assumed 37 orbits in 20 sidereal days, consistent with the
  paper's "around 13 hours"), peak-to-peak 2√(GMa)e/c² = 381 ns vs the paper's 400 ns; F to 10⁻⁸.
  Crate: `geopotential::chronometric::orbit_proper_time` integrating V from `Coefficients` along an
  orbit from `gravimetry::sim`. **Links geopotential ↔ gravimetry ↔ metrology.**
- **R2. Solar redshift [established effect; 2020 measurement].** v = GM☉/(cR☉) − GM☉/(c·AU) −
  GM⊕/(cR⊕); expected ≈ 633.1 m/s; measured 638 ± 6 m/s (weak Fe lines), 639 ± 14 m/s (15 strong
  lines) with HARPS and a laser frequency comb on Moon spectra — [González Hernández et al. 2020,
  A&A 643, A146, arXiv:2009.10558](https://arxiv.org/abs/2009.10558). Recomputed: 633.14 m/s; 0.81σ.
  Crate: `gravitics::ppn` (redshift tests).
- **R3. S2 at Sgr A* [recent and replicated].** z = GM/(c²r_p) + v_p²/(2c²) at pericentre;
  Schwarzschild precession Δφ = 6πGM/(c²a(1−e²)) per orbit. Published: f = 0.90 ± 0.09(stat) ±
  0.15(sys), v_p ≈ 7650 km/s, z ≈ 200 km/s / c, pericentre 120 AU (May 2018) —
  [GRAVITY 2018, arXiv:1807.09409](https://arxiv.org/abs/1807.09409); f_SP = 1.10 ± 0.19,
  "12′ per orbital period" — [GRAVITY 2020, arXiv:2004.07187](https://arxiv.org/abs/2004.07187);
  f_SP = 1.135 ± 0.110 (≈ 10σ), M = 4.2996 ± 0.0118 ×10⁶ M☉, R₀ = 8275.9 ± 8.6 pc,
  a = 0.12502″, e = 0.88444, extended mass < ~1200 M☉ — [GRAVITY 2024, A&A 692, A242,
  arXiv:2409.12261](https://arxiv.org/html/2409.12261). Recomputed from the 2024 elements: P = 16.05
  yr; Δφ = 12.21′; r_p = 119.6 AU; v_p = 7753 km/s; z·c = 206.7 km/s. Crate: `gravitics::ppn`;
  the same 6πGM/(c²a(1−e²)) is the periastron advance of P4. Links gravitics ↔ gravitational.
- **R4. Gravity Probe A comparison.** α ≈ 1.4×10⁻⁴ (1976), superseded 5.6× by Galileo —
  [Delva et al. 2018](https://arxiv.org/abs/1812.03711).

### Inferences
- Every redshift test here (Galileo, Skytree, Sun, S2) is one function of (U, v) with an α factor;
  a single `redshift(u_a, v_a, u_b, v_b, alpha)` with a registry of α bounds would cover them.
- R1 is the strongest cross-crate case: an orbit simulator (gravimetry), a field (geopotential)
  and constants (metrology) produce a measured nanosecond-level signal.

### Gaps
- Herrmann et al.'s α value was not fetched; (4.5 ± 3.1)×10⁻⁵ appears only in a slide quoted by
  a search summary.
- The Galileo semi-major axis was inferred from the period, not taken from the paper.
- Keck's independent S2 redshift (Do et al. 2019) was not fetched.

## 5. Frame dragging and light; black-hole shadows and photon rings

### Takeaway
The Schwarzschild shadow diameter 6√3 θ_g (θ_g = GM/(c²D)) reproduces the EHT: Sgr A* 53.3 µas
(VLTI prior) or 51.1 µas (Keck) against a 51.8 ± 2.3 µas ring and δ = −0.08 / −0.04; M87* 39.7 µas
against a 42 ± 3 µas ring (2017) and 43.9 ± 0.6 µas (2017–2021). Kerr shrinks the shadow to at
least 2(1+√2) ≈ 4.83 M in radius (pole-on, a = 1), matching EHT's quoted 4.8–5.2 M range. Light's
own frame-dragging signatures (gravitational Faraday rotation, gravitational spin Hall) are
unmeasured and the Faraday coefficient is disputed.

### Cited Findings
- **B1. Shadow size and the EHT deviation parameter [recent; M87* replicated over years, Sgr A*
  unreplicated].** d̂_m = α_c(1+δ)·6√3 θ_g, δ ≡ d_sh/d_sh,Sch − 1. Priors θ_g = 5.125 ± 0.009 ±
  0.020 µas (VLTI; M = 4.297×10⁶ M☉, R₀ = 8277 pc) and 4.92 ± 0.03 ± 0.01 µas (Keck;
  M = 3.951×10⁶ M☉, R₀ = 7935 pc); δ = −0.08 +0.09/−0.09 (VLTI) and −0.04 +0.09/−0.10 (Keck);
  Table 1 shadow ≈ 48–50 µas by method; Kerr shadow radius ∼4.8–5.2 GM/c²; image size within ~10%
  of Kerr — [EHT Sgr A* Paper VI, ApJL 930, L17, arXiv:2311.09484](https://arxiv.org/abs/2311.09484).
  Ring 51.8 ± 2.3 µas — [EHT Sgr A* Paper I, ApJL 930, L12, arXiv:2311.08680](https://arxiv.org/abs/2311.08680).
  M87*: ring 42 ± 3 µas, M = (6.5 ± 0.7)×10⁹ M☉ — [EHT M87* Paper I 2019, arXiv:1906.11238](https://arxiv.org/abs/1906.11238);
  persistent 43.9 ± 0.6 µas across 2017, 2018, 2021 — [EHT 2025, arXiv:2509.24593](https://arxiv.org/abs/2509.24593).
  Recomputed: θ_g from GRAVITY 2024 M and R₀ = 5.128 µas; 6√3θ_g = 53.3 (VLTI), 51.1 (Keck) µas;
  (1+δ)6√3θ_g = 49.0 and 49.1 µas (inside Table 1's range); M87* θ_g = 3.82 µas (EHT 3.8 ± 0.4) and
  ring/θ_g = 11.0–11.5. Crate: `gravitics::spacetime` (`shadow_diameter(m, d)`,
  `deviation(d_ring, alpha_c, theta_g)`).
- **B2. Kerr photon orbits and the pole-on shadow [established theory; foundational].** Equatorial
  circular photon orbits r_ph = 2M{1 + cos[(2/3)arccos(∓a)]} ([Bardeen, Press & Teukolsky 1972,
  ApJ 178, 347](https://doi.org/10.1086/151796)); spherical orbits ξ = [r²(3M−r) − a²(r+M)]/[a(r−M)],
  η = r³[4a²M − r(r−3M)²]/[a²(r−M)²]; a polar observer sees ξ = 0, radius √(η + a²) ([Teo 2003,
  GRG 35, 1909](https://doi.org/10.1023/A:1026286607562)). Recomputed: √27 at a → 0; exactly
  2(1+√2) = 4.828 M at a = 1; 5.121 M at a = 0.5; 4.884 M at a = 0.94; prograde/retrograde photon
  orbits 1 M and 4 M at a = 1. Crate: `gravitics::spacetime`; shares Bardeen's orbit algebra with
  `gravitational::inspiral::f_isco`. Links gravitics ↔ gravitational.
- **B3. Photon-ring subrings [theory established; detection contested].** Subrings become
  "exponentially narrower but weaker with increasing orbit number" with universal signatures on long
  baselines — [Johnson et al. 2020, Sci. Adv. 6, eaaz1310, arXiv:1907.04329](https://arxiv.org/abs/1907.04329).
  The demagnification per half orbit is e^γ with Lyapunov exponent γ; recomputed for Schwarzschild
  from λ = Ω = 1/(3√3 M): γ = λ·(π/Ω) = π exactly, so each subring is 23.14× thinner and fainter.
  Crate: `gravitics::spacetime`.
- **B4. Gravitational Faraday (Rytov–Skrotskii) rotation [contested theory; unmeasured].** The
  weak-field rotation of the polarisation plane by a spinning lens is a purely gravitomagnetic
  effect and its detectability has been studied — [Sereno 2005, arXiv:astro-ph/0410015](https://arxiv.org/abs/astro-ph/0410015);
  in magnetised Kerr spacetimes the weak-field rotation scales as aM²/R³ — [arXiv:2106.03520](https://arxiv.org/abs/2106.03520);
  standard reference [Ishihara, Takahashi & Tomimatsu 1988, PRD 38, 472](https://doi.org/10.1103/PhysRevD.38.472).
  A search summary reported that the numerical coefficient is still disputed between authors; I
  could not confirm a single coefficient. Crate: if added, `gravitics::gem` (it is the light-side
  partner of the existing gravitomagnetic field).
- **B5. Gravitational spin Hall effect of light [recent theory; unmeasured].** Polarisation-dependent
  ray equations from a covariant WKB analysis, with numerical Schwarzschild examples —
  [Oancea et al. 2020, PRD 102, 024075](https://collaborate.princeton.edu/en/publications/gravitational-spin-hall-effect-of-light/).
  No closed-form weak-field estimate was verified.

### Inferences
- B1 and B2 make a compact, fully tested addition: one function for 6√3θ_g, one for the Kerr
  pole-on radius, one for equatorial photon orbits; the EHT numbers and the exact 2(1+√2) are the
  tests.
- The crate's existing gravitomagnetic machinery (`gravitics::gem`) is for matter; B4 is the only
  light-side frame-dragging effect with a formula, and it should not be added until a coefficient is
  pinned from the primary source.

### Gaps
- The Sgr A* shadow value 48.7 ± 7.0 µas often quoted from Paper I was not found in the fetched
  abstracts; the Paper VI table values are used instead.
- The claimed n = 1 photon-ring detection in M87* and its rebuttal were not fetched.
- No primary-source coefficient for gravitational Faraday rotation was verified.

## 6. Analogue gravity with light

### Takeaway
The fibre-optical horizon gives a Hawking temperature k_BT = ħα/(2π) with α = −(1/δn)∂δn/∂τ, which
is exactly the crate's existing `horizon_temperature(kappa)` with κ = α; Philbin et al.'s ~10³ K
estimate recomputes to 973 K. Drori et al. 2019 observed stimulated Hawking emission at 231 nm for
a 1450 nm probe; a 2026 Nature paper from the same group reports its backreaction. The Doppler rule
ω′ = γ(1 − nu/c)ω explains why the UV partner has negative co-moving frequency; bulk silica
reproduces the sign and the ~10⁻² group–phase mismatch, but the measured wavelengths need the
fibre's own dispersion.

### Cited Findings
- **A1. Fibre horizon temperature [established theory; foundational].** Horizon where the pulse
  speed u equals the probe group velocity c/(n + ω dn/dω); co-moving frequency ω′ = (1 − nu/c)ω is
  conserved; linearising 1 − nu/c = α′τ gives k_BT′ = ħα′/(2π), α′ = −(u/c)∂δn/∂τ, and in the lab
  frame k_BT = ħα/(2π), α = −(1/δn)∂δn/∂τ; "the Hawking temperature (8) reaches 10³ K" for a shock
  steepness of 8×10¹⁴ Hz — [Philbin et al. 2008, Science 319, 1367, arXiv:0711.4796](https://arxiv.org/abs/0711.4796).
  Recomputed: 972.5 K. Already covered: `gravitics::spacetime::horizon_temperature(alpha)`. New:
  `fibre_surface_gravity(delta_n, d_delta_n_dtau)` and `comoving_frequency(omega, n, u)`. Links
  gravitics ↔ physica analogies (the same horizon formula as Unruh, Hawking, Gibbons–Hawking and the
  existing Schwinger–Unruh identity).
- **A2. Stimulated Hawking radiation in a fibre [unreplicated by another group; outside the window,
  2019].** ω′ = γ(1 − nu/c)ω with γ⁻² = 1 − u²/c²; Bogoliubov b̂± = αâ± + βâ†∓; 8 fs pump at
  800 nm in a 7 mm PCF (NKT NL-1.5-590), pump group velocity fitted to an 818.9 nm carrier; δn ≈ 10⁻³
  while group and phase velocity differ by ≈ 10⁻²; stimulated negative-frequency peak at 231 nm for
  a 1450 nm probe; 41 000 extra counts/s for 2×10⁷ probe photons, i.e. 2×10⁻³ spontaneous partners/s
  at ~10⁻² detection efficiency; horizon wavelength measured 1551 nm vs 1613 nm predicted; signal
  linear in probe power up to ≈1.5% of pump peak — [Drori et al. 2019, PRL 122, 010404,
  arXiv:1808.09244](https://arxiv.org/abs/1808.09244). Recomputed with fused-silica Sellmeier
  ([Malitson 1965, JOSA 55, 1205](https://doi.org/10.1364/JOSA.55.001205)): at 818.9 nm n = 1.4530,
  n_g = 1.4666, (n_g − n)/n = 0.0093 (paper ≈ 10⁻²); ω′/ω = +0.0146 at 1450 nm and −0.0361 at 231 nm
  (signs as the paper states). Crate: `gravitics::spacetime` (analogue horizons) — no physical
  gravity, so it belongs beside the horizon temperatures.
- **A3. 2026 follow-up [recent and unreplicated].** Procopio, Aguero-Santacruz, Bermudez and
  Leonhardt report a simple, direct process for fibre Hawking radiation (not the cascaded process
  previously assumed) and measured its backreaction; Nature, online 1 July 2026 —
  [arXiv:2607.01118](https://arxiv.org/abs/2607.01118). No numbers in the abstract.

### Inferences
- A1 adds almost no code: the crate already has the temperature; only α from a δn profile and the
  Doppler rule are new, and both have published numbers.
- A physica analogy doc test could assert that the fibre, Unruh, Hawking and Gibbons–Hawking
  temperatures are one `horizon_temperature` call with four different κ.

### Gaps
- The PCF's measured n(ω) is not public in the paper, so the 231 nm partner wavelength for a
  1450 nm probe cannot be recomputed; bulk silica only checks signs and magnitudes.
- Spontaneous optical Hawking radiation has not been reported; condensate analogues
  (Steinhauer group) were out of scope and not fetched.

## 7. Lunar laser ranging tests

### Takeaway
LLR's latest global solution gives Δ(m_g/m_i)_EM = (−2.1 ± 2.4)×10⁻¹⁴, Ġ/G = (−5.0 ± 9.6)×10⁻¹⁵ /yr,
β − 1 = (6.2 ± 7.2)×10⁻⁵ and γ − 1 = (1.7 ± 1.6)×10⁻⁴; the implied Nordtvedt η ≈ 5×10⁻⁵ corresponds
to a 0.6 mm range signal. A 2023 PRL bounds active/passive mass inequality at 3.9×10⁻¹⁴. No newer
LLR solution (2024–2026) was found.

### Cited Findings
- **M1. Nordtvedt effect [established null result].** m_p/m = 1 − η_N E_g/m with
  η_N = 4β − γ − 3 − (10/3)ξ − α₁ + (2/3)α₂ − (2/3)ζ₁ − (1/3)ζ₂; E_g/m = 4.6×10⁻¹⁰ (Earth),
  0.2×10⁻¹⁰ (Moon), 3.6×10⁻⁶ (Sun); range signal δr = 13.1 η_N cos(ω₀ − ω_s)t m; older bound
  |η_N| = (4.4 ± 4.5)×10⁻⁴ — [Will 2014, eqs. 66–67](https://arxiv.org/abs/1403.7377). Fully
  conservative form η = 4(β − 1) − (γ − 1) − α₁ − α₂ — [Fienga & Minazzoli 2024, eq. 126](https://arxiv.org/abs/2303.01821).
  Crate: `gravitics::ppn` (`nordtvedt_eta(beta, gamma)`, `llr_range_signal(eta)`); the
  self-energies E_g/m could come from a density model in `gravitics::anomalies`
  (which already holds column densities). Links gravitics::anomalies ↔ gravitics::ppn.
- **M2. Latest LLR solution [established, 2021 (just before the window)].**
  Ġ/G₀ = (−5.0 ± 9.6)×10⁻¹⁵ yr⁻¹; Δ(m_g/m_i)_EM = (−2.1 ± 2.4)×10⁻¹⁴; β − 1 = (6.2 ± 7.2)×10⁻⁵;
  γ − 1 = (1.7 ± 1.6)×10⁻⁴ — [Biskupek, Müller & Torre 2021, Universe 7, 34, arXiv:2012.12032](https://arxiv.org/abs/2012.12032);
  the 2024 review calls it more accurate than planetary ephemerides and quotes "a deduced η value of
  5×10⁻⁵" — [Fienga & Minazzoli 2024](https://arxiv.org/abs/2303.01821). Recomputed:
  η = 2.1×10⁻¹⁴/(4.6 − 0.2)×10⁻¹⁰ = 4.8×10⁻⁵ ± 5.5×10⁻⁵ (matches the review's 5×10⁻⁵); range signal
  13.1η = 0.63 mm; η from the same fit's β, γ: (0.8 ± 3.3)×10⁻⁴.
- **M3. Active vs passive mass [recent and unreplicated].** LLR 1970–2022 bounds the Al/Fe
  active–passive mass difference at 3.9×10⁻¹⁴, about 100× better than Bartlett & Van Buren 1986 —
  [Singh et al. 2023, PRL 131, 021401, arXiv:2212.09407](https://arxiv.org/abs/2212.09407).
- **M4. de Sitter precession of the lunar orbit [established].** Ω_G = (γ + ½)v × ∇U; "about 2
  arcseconds per century", measured to ~0.6% with LLR — [Will 2014, eq. 70](https://arxiv.org/abs/1403.7377).
  Recomputed: (3/2)(GM☉/(c²·AU))n⊕ = 1.919″/century = 19.2 mas/yr. Crate: `gravitics::ppn`;
  Gravity Probe B's 6.63″/yr geodetic term is the same formula in Earth orbit and would link to
  `gravitics::superconductors` (frame dragging), Will quoting 0.3% (geodetic) and 20% (frame
  dragging) for GP-B.
- **M5. Ġ/G from ephemerides [recent, secondary].** A 2025 review quotes |Ġ/G| in the
  (2–6)×10⁻¹⁴ yr⁻¹ range from ephemerides and LLR combinations and projects 4–10× gains from
  next-generation retroreflectors — [arXiv:2509.05910](https://arxiv.org/abs/2509.05910).

### Inferences
- M1–M4 and P1–P3 form one `gravitics::ppn` module: γ, β, η, α (redshift) and Ġ/G with formulas
  and a status-labelled registry, mirroring how `gravitics::anomalies` handles the CODATA G values.
- LLR already appears in `gravitics::anomalies` for the shielding bound; the same data's η and Ġ/G
  would sit naturally beside it.

### Gaps
- A 2023 ILRS workshop slide reportedly lists η = (−0.2 ± 1.1)×10⁻⁴; it is a presentation, not a
  paper, and the value was seen only in a search summary.
- No 2024–2026 LLR global solution was found.
