# Geomagnetism, space weather and astrophysical magnetism: formulas and published test values for verified Rust calculators

Check script: `research/topics/magnetics-and-mhd/notes/geomagnetism_checks.py` (stdlib only). It was run from the
session scratchpad, outside the notes folder:
`python3 -I .../notes/geomagnetism_checks.py /home/dan/code/geomagnetic/WMM2025COF/WMM2025_TestValues.txt`.
Result: **all checks passed** (38 published-value checks: 12 WMM2025 Table 6 points × 15 elements, all 100 rows of
NOAA's high-precision file, poles, eccentric dipole, SAA minimum, ATNF pulsars, RM constant, Schwinger field and
others). Without the optional path it checks 6 embedded high-precision rows instead of 100. The script embeds
WMM2025.COF exactly as distributed, plus the IGRF-14 DGRF-2020, IGRF-2025 and SV-2025-30 coefficients from the
IGRF-14 paper's Table 2.

Local primary copies already on this machine (the user's `geomagnetic` aurora-globe project, which has a TypeScript
WMM evaluator and an `npm run validate` against NOAA test values): `/home/dan/code/geomagnetic/WMM2025COF/`
(README, WMM2025.COF = WMM.COF, WMM2025_TestValues.txt/.pdf, dated 11 Dec 2024) and the full WMM2025 technical report
`/home/dan/code/geomagnetic/noaa_71569_DS1.pdf` (109 pp., created Aug 2025). The IGRF-14 paper PDF was fetched from
NORA during this research.

---

## 1. IGRF-14 and WMM2025: the potential, normalization, sign conventions, coordinates, secular variation (and how to run them on the `geopotential` engine)

### Takeaway
Both models use V = a Σₙ (a/r)ⁿ⁺¹ Σₘ [gₙᵐ cos mλ + hₙᵐ sin mλ] P̆ₙᵐ(sin φ′) with a = 6 371 200 m, Schmidt
semi-normalized P̆ₙᵐ with no Condon–Shortley phase, geocentric (r, φ′, λ), and linear secular variation. The 4π
fully normalized functions of the `geopotential` crate are P̄ₙₘ = √(2n+1) P̆ₙᵐ, so the gravity engine evaluates the
main field exactly if it is given GM := a², R := a, C̄ₙₘ = gₙᵐ/√(2n+1), S̄ₙₘ = hₙᵐ/√(2n+1), C̄₀₀ = 0, and the output
is mapped as X′ = −∇V_N, Y′ = −∇V_E, Z′ = +∇V_U. The check script shows the two routes agree to 6×10⁻¹¹ nT and
reproduce every WMM2025 test value.

### Cited Findings
- WMM main field is "the negative spatial gradient of a scalar potential" in geocentric (λ, φ′, r), expanded to
  N = 12 (WMM) or N = 133 (WMMHR), with a = 6,371,200 m, "Schmidt semi-normalized associated Legendre functions",
  and the definition the report says is used in geodesy and geomagnetism (Heiskanen & Moritz 1967 eq. 1–60;
  Langel 1987 eq. 8). The report adds that its P̆ relate to the Abramowitz & Stegun / Gradshteyn & Ryzhik functions
  by a stated factor (the equation is an image in the PDF) — [Chulliat et al. 2025, WMM2025 Technical Report §1.2](https://doi.org/10.25923/prbc-s316)
- IGRF-14 uses the same expansion with a = 6371.2 km "chosen to approximate the mean Earth radius", P̆ₙᵐ(cos θ)
  "Schmidt semi-normalized" (citing Winch et al. 2005), N = 10 up to and including 1995 and 13 afterwards —
  [Beggan et al. 2026, IGRF-14, EPS 78:127](https://doi.org/10.1186/s40623-025-02360-0)
- WMM procedure (§1.2): (1) height above MSL → height above WGS 84 via EGM96 (effect "of the order of 1 nT or less",
  optional); (2) geodetic (φ, h) → geocentric (φ′, r) with WGS 84 A = 6 378 137 m, 1/f = 298.257223563,
  R_c = A/√(1 − e² sin²φ), p = (R_c + h) cos φ, z = (R_c(1 − e²) + h) sin φ, r = √(p² + z²), φ′ = asin(z/r);
  (3) g(t) = g(t₀) + (t − t₀) ġ with t₀ = 2025.0 and t in decimal years; (4) X′, Y′, Z′; (5) rotate to the
  ellipsoidal frame; (6) H, F, I, D and their rates, with atan2 giving D ∈ (−π, π] and I ∈ [−π/2, π/2] —
  [WMM2025 Technical Report §1.2](https://doi.org/10.25923/prbc-s316)
- IGRF time dependence is piecewise linear between 5-year epochs, g(t) = g(T_t) + (t − T_t) ġ(T_t), with
  ġ(T_t) = [g(T_t + 5) − g(T_t)]/5 for T_t < 2025 and the predictive SV for 2025–2030 —
  [Beggan et al. 2026](https://doi.org/10.1186/s40623-025-02360-0)
- Sign conventions: declination east of true north positive; inclination positive when the north-seeking end
  dips into the ground; F = √(X² + Y² + Z²) — [Beggan et al. 2026](https://doi.org/10.1186/s40623-025-02360-0).
  Table 1 of the WMM report gives positive senses X north, Y east, Z down, I down, D east/clockwise and
  surface ranges X −17 000…43 000, Y −18 000…17 000, Z −67 000…62 000, H 0…43 000, F 23 000…67 000 nT —
  [WMM2025 Technical Report Table 1](https://doi.org/10.25923/prbc-s316)
- Grid variation: MIL-PRF-89500B defines two auxiliary polar angles, grid variation north (GVN) and south (GVS), from
  D, longitude and geodetic latitude (the defining equations are images in the PDF); Table 6 GV is "with respect to
  the Grid North of the Universal Polar Stereographic projection" — [WMM2025 Technical Report §1.1.2, Table 6](https://doi.org/10.25923/prbc-s316)
- WMM singularity: X′, Y′, Z′ are referred to a north–east–down frame undefined exactly at a geographic pole; the
  report fixes it by assigning the pole a longitude and gives the polar field at r = 6 371 200 m on 1 Jan 2025
  (WMM: X′, Y′, Z′ = 1709.5, 417.8, 56517.5 nT at the North Pole with λ = 0°) — [WMM2025 Technical Report §1.4](https://doi.org/10.25923/prbc-s316)
- The `geopotential` crate in this repo stores "fully normalized (4π, no Condon–Shortley phase) Stokes coefficients"
  and its `Partials::acceleration_neu(r, lat)` returns ∇V in (north, east, up) as
  `[dlat/r, dlon/(r cos lat), dr]` with geocentric latitude; `second_partials` is "finite at the poles" —
  `crates/geopotential/src/coefficients.rs`, `crates/geopotential/src/synthesis.rs`, `crates/geopotential/src/legendre.rs`
  (local source).

**Verified mapping (check script §1c, two independent recursions):**

| Gravity engine (`geopotential`) | Magnetic main field |
|---|---|
| `gm` | a² = 6 371 200² (units become nT·m², the engine does not care) |
| `radius` | a = 6 371 200 m |
| C̄₀₀ | 0 (no monopole) |
| C̄ₙₘ, S̄ₙₘ | gₙᵐ/√(2n+1), hₙᵐ/√(2n+1) (hₙ⁰ = 0) |
| `partials(r, φ′, λ).acceleration_neu(r, φ′)` = (N, E, U) | X′ = −N, Y′ = −E, Z′ = +U (since B = −∇V, Z is down) |
| same with ġ, ḣ as coefficients | Ẋ′, Ẏ′, Ż′ (synthesis is linear in the coefficients) |
| `second_partials` | magnetic gradient tensor ∂Bᵢ/∂xⱼ = −∂²V/∂xᵢ∂xⱼ |

Then rotate with ψ = φ′ − φ: X = X′ cos ψ − Z′ sin ψ, Y = Y′, Z = X′ sin ψ + Z′ cos ψ (same for the rates);
H = √(X²+Y²), F = √(H²+Z²), D = atan2(Y, X), I = atan2(Z, H); Ḣ = (XẊ+YẎ)/H, Ḟ = (XẊ+YẎ+ZŻ)/F,
Ḋ = (XẎ−YẊ)/H², İ = (HŻ−ZḢ)/F²; GV = D − λ (north polar), D + λ (south polar), wrapped to (−180°, 180°]. This
form reproduces every GV in Table 6 (all at ±80°); the 55° latitude cut-off used in the script is an assumption, not a
cited value.

Gaussian conversions: 1 nT = 1 γ = 10⁻⁵ G; 1 G = 10⁻⁴ T; a moment of 1 A m² = 10³ G cm³ (emu); the IGRF moment
formula uses μ₀ = 4π×10⁻⁷ H/m ([Beggan et al. 2026](https://doi.org/10.1186/s40623-025-02360-0)); the CODATA 2022 μ₀
differs by ~10⁻¹⁰ relative, which is irrelevant here.

### Inferences
- No (−1)ᵐ correction is needed between WMM/IGRF and the crate: neither includes the Condon–Shortley phase. A
  library whose associated Legendre functions include it (common in maths libraries following Abramowitz & Stegun,
  which the WMM report singles out as differently defined) would flip the sign of every odd-m term. The
  check script's Schmidt recursion and 4π recursion both omit it and agree to 6×10⁻¹¹ nT.
- The crate's `.gfc` reader is not usable for WMM.COF, igrf14coeffs.txt or SHC files; a small parser that fills
  `Coefficients::zeros(a², a, nmax)` with gₙᵐ/√(2n+1) is all that is needed. The `TideSystem` field is meaningless
  for magnetics and should be left `Unknown`.
- `acceleration_neu` is undefined exactly at a pole (east has no direction), which matches the WMM §1.4 caveat; use
  the report's pole convention if a polar value is needed.

### Gaps
- The exact Schmidt definition and the A&S relation in the WMM report are images that did not extract; the
  definition used here (P̆ₙᵐ = √(2(n−m)!/(n+m)!) Pₙ,ₘ for m > 0, Pₙ for m = 0) is verified numerically against the
  test values rather than read from the page.

---

## 2. Published test values, coefficient files, formats and terms of use

### Takeaway
WMM2025 ships 12 report test points (Table 6, with GV and rates) and a 100-row high-precision file; this script
reproduces all Table 6 values within the published rounding (worst 0.050 nT, 0.005°) and all 100 file rows to
≤ 0.0008 nT. IGRF-14 has no official test-point table; its paper's SAA minimum, pole positions and dipole statements
are the checkable values. NOAA states the WMM is public domain; the IGRF-14 paper is CC BY 4.0 and says the
coefficients and code are "freely available for use".

### Cited Findings
**WMM2025 Table 6 (exact, from the report; heights above the WGS 84 ellipsoid; "double precision arithmetic.
Single precision arithmetic can cause differences of up to 0.1 nT")** — [WMM2025 Technical Report Table 6](https://doi.org/10.25923/prbc-s316)

| Date | h km | Lat | Lon | X | Y | Z | H | F | I° | D° | GV° |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 2025.0 | 0 | 80 | 0 | 6521.6 | 145.9 | 54791.5 | 6523.2 | 55178.5 | 83.21 | 1.28 | 1.28 |
| 2025.0 | 0 | 0 | 120 | 39677.8 | −109.6 | −10580.2 | 39677.9 | 41064.3 | −14.93 | −0.16 | |
| 2025.0 | 0 | −80 | 240 | 6117.5 | 15751.9 | −52022.5 | 16898.1 | 54698.2 | −72.00 | 68.78 | −51.22 |
| 2025.0 | 100 | 80 | 0 | 6216.0 | 92.4 | 52598.8 | 6216.7 | 52964.9 | 83.26 | 0.85 | 0.85 |
| 2025.0 | 100 | 0 | 120 | 37688.6 | −96.2 | −10152.1 | 37688.7 | 39032.1 | −15.08 | −0.15 | |
| 2025.0 | 100 | −80 | 240 | 5907.6 | 14780.3 | −49540.7 | 15917.1 | 52035.0 | −72.19 | 68.21 | −51.79 |
| 2027.5 | 0 | 80 | 0 | 6500.8 | 294.5 | 54869.4 | 6507.5 | 55253.9 | 83.24 | 2.59 | 2.59 |
| 2027.5 | 0 | 0 | 120 | 39701.6 | −167.4 | −10381.8 | 39702.0 | 41036.9 | −14.65 | −0.24 | |
| 2027.5 | 0 | −80 | 240 | 6200.7 | 15730.3 | −51783.7 | 16908.3 | 54474.2 | −71.92 | 68.49 | −51.51 |
| 2027.5 | 100 | 80 | 0 | 6196.7 | 233.8 | 52670.5 | 6201.1 | 53034.3 | 83.29 | 2.16 | 2.16 |
| 2027.5 | 100 | 0 | 120 | 37711.5 | −148.7 | −9969.8 | 37711.8 | 39007.4 | −14.81 | −0.23 | |
| 2027.5 | 100 | −80 | 240 | 5984.0 | 14760.1 | −49317.7 | 15927.0 | 51825.7 | −72.10 | 67.93 | −52.07 |

Rates (nT/yr; °/yr), same rows in order:

| Date | h | Lat | Lon | Ẋ | Ẏ | Ż | Ḣ | Ḟ | İ | Ḋ |
|---|---|---|---|---|---|---|---|---|---|---|
| 2025.0 | 0 | 80 | 0 | −8.3 | 59.5 | 31.1 | −7.0 | 30.1 | 0.01 | 0.52 |
| 2025.0 | 0 | 0 | 120 | 9.5 | −23.1 | 79.4 | 9.6 | −11.2 | 0.11 | −0.03 |
| 2025.0 | 0 | −80 | 240 | 33.3 | −8.6 | 95.5 | 4.0 | −89.6 | 0.03 | −0.12 |
| 2025.0 | 100 | 80 | 0 | −7.7 | 56.5 | 28.7 | −6.9 | 27.6 | 0.01 | 0.52 |
| 2025.0 | 100 | 0 | 120 | 9.2 | −21.0 | 72.9 | 9.2 | −10.0 | 0.11 | −0.03 |
| 2025.0 | 100 | −80 | 240 | 30.6 | −8.0 | 89.2 | 3.9 | −83.8 | 0.03 | −0.11 |
| 2027.5 | 0 | 80 | 0 | −8.3 | 59.5 | 31.1 | −5.6 | 30.3 | 0.01 | 0.53 |
| 2027.5 | 0 | 0 | 120 | 9.5 | −23.1 | 79.4 | 9.6 | −10.7 | 0.11 | −0.03 |
| 2027.5 | 0 | −80 | 240 | 33.3 | −8.6 | 95.5 | 4.2 | −89.5 | 0.04 | −0.12 |
| 2027.5 | 100 | 80 | 0 | −7.7 | 56.5 | 28.7 | −5.6 | 27.8 | 0.01 | 0.52 |
| 2027.5 | 100 | 0 | 120 | 9.2 | −21.0 | 72.9 | 9.3 | −9.7 | 0.11 | −0.03 |
| 2027.5 | 100 | −80 | 240 | 30.6 | −8.0 | 89.2 | 4.0 | −83.7 | 0.03 | −0.11 |

- Report Table 7 gives the same 12 points for WMMHR2025 (e.g. 2025.0, 0 km, 80°N, 0°: X 6517.4, Y 144.8,
  Z 54701.3, F 55088.3 nT, D 1.27°) — [WMM2025 Technical Report Table 7](https://doi.org/10.25923/prbc-s316).
  Not checked (WMMHR coefficients to degree 133 not embedded).
- **High-precision file** `WMM2025_TestValues.txt` (in WMM2025COF.zip): 100 rows, 2025.0–2029.5, fields
  year, altitude km (HAE), geodetic lat, lon, D, I (2 decimals), H, X, Y, Z, F (6 decimals) and seven rates. First
  rows, verbatim: `2025.000000 28 89 -121 -99.77 88.47 1504.298146 -255.388723 -1482.460628 56194.288771 56214.419888`;
  `2025.000000 48 80 -96 -29.91 87.77 2164.285547 1875.982280 -1079.269389 55623.044051 55665.134163`;
  `2025.500000 6 -36 -137 20.28 -52.11 25353.230658 23781.930678 8786.698927 -32577.518648 41280.516301` —
  NOAA WMM2025COF package (local copy `/home/dan/code/geomagnetic/WMM2025COF/WMM2025_TestValues.txt`; distributed via
  [NCEI WMM page](https://www.ncei.noaa.gov/products/world-magnetic-model) as "WMM2025 Test Values")
- **Script result:** Table 6 worst |computed − published| per element: X 0.048, Y 0.044, Z 0.049, H 0.050,
  F 0.045 nT; I, D, GV 0.005°; rates ≤ 0.050 nT/yr and 0.005°/yr — all within rounding. 100-row file: max 0.00072 nT
  in X, H, Y, Z, F and ≤ 0.005° in D, I. The 0.0007 nT (10⁻⁸ relative) residual is present even at 2025.0, so it is
  not a time-handling issue; its source in NOAA's generator is unknown.
- **WMM.COF format**: header `    2025.0            WMM-2025     11/13/2024` (epoch, model name, release date),
  then 90 rows `n m gₙᵐ hₙᵐ ġₙᵐ ḣₙᵐ` (nT, nT/yr; hₙ⁰ = 0), terminated by two lines of 9s. First rows:
  `1 0 -29351.8 0.0 12.0 0.0`; `1 1 -1410.8 4545.4 9.7 -21.5`; `2 0 -2556.6 0.0 -11.6 0.0`;
  `2 1 2951.1 -3133.6 -5.2 -27.7`; `2 2 1649.3 -815.1 -8.0 -12.1` — local WMM2025.COF; identical to report Table 2
  ([WMM2025 Technical Report](https://doi.org/10.25923/prbc-s316)). Final WMM coefficients are the average of the
  NCEI and BGS candidates, rounded to one decimal (WMM) / four decimals (WMMHR) — same report §2.
- WMMHR2025 extends to degree 133; degrees 16–133 come from the MF7 crustal model and their SV is zero —
  [WMM2025 Technical Report §2](https://doi.org/10.25923/prbc-s316)
- **IGRF-14 file format and values**: `igrf14coeffs.txt` holds Schmidt semi-normalised coefficients in nT, degrees
  1–13, epochs 1900.0–2025.0 then a "2025-30" SV column in nT/yr — [igrf14coeffs.txt](https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/igrf14coeffs.txt).
  Low-degree rows read from that file (DGRF 2020.0 | IGRF 2025.0 | SV): g₁⁰ −29403.41 | −29350.0 | 12.6;
  g₁¹ −1451.37 | −1410.3 | 10.0; h₁¹ 4653.35 | 4545.5 | −21.5; g₂⁰ −2499.78 | −2556.2 | −11.2;
  g₂¹ 2981.96 | 2950.9 | −5.3; h₂¹ −2991.72 | −3133.6 | −27.3; g₂² 1676.85 | 1648.7 | −8.3;
  h₂² −734.62 | −814.2 | −11.1; g₁⁰ at 1900.0 −31543 and 2000.0 −29619.4. These match the paper's Table 2 row for
  row (two independent extractions) — [Beggan et al. 2026, Table 2](https://doi.org/10.1186/s40623-025-02360-0)
- IGRF-14 = DGRF-2020 + IGRF-2025 (both to degree 13) + SV-2025-2030 (to degree 8); valid 1900.0–2030.0; adopted
  November 2024 by IAGA Division V Working Group V-MOD; candidates from 19 teams / 27 institutes (14 DGRF, 15 IGRF,
  18 SV candidates) — [Beggan et al. 2026](https://doi.org/10.1186/s40623-025-02360-0); evaluation paper
  [Beggan et al. 2026, EPS 78:126](https://doi.org/10.1186/s40623-026-02382-2); special issue
  [EPS IGRF-14 collection](https://www.springeropen.com/collections/igrf14)
- IGRF-14 distribution: ASCII coefficients and calculators at [NCEI IGRF](https://www.ncei.noaa.gov/products/international-geomagnetic-reference-field),
  Zenodo coefficients `doi:10.5281/zenodo.14218973`, candidate models and evaluation code at
  `github.com/IAGA-VMOD/IGRF14eval`, Python software (Laundal et al. 2024) `doi:10.5281/zenodo.14231854`, BGS and
  Kyoto online calculators — [Beggan et al. 2026 §6](https://doi.org/10.1186/s40623-025-02360-0)
- WMM2025 and WMMHR2025 were released 17 Dec 2024; NCEI: "The WMM source code is in the public domain and not
  licensed or under copyright", "The information and software may be used freely by the public", with a 17 U.S.C.
  403 notice requirement for third-party copyrighted works incorporating it — [NCEI WMM page](https://www.ncei.noaa.gov/products/world-magnetic-model).
  The report is "Approved for public release, NGA-U-2025-01497"; dataset DOI 10.25921/aqfd-sd83 —
  [WMM2025 Technical Report](https://doi.org/10.25923/prbc-s316)
- IGRF-14 paper is CC BY 4.0; "The IGRF-14 coefficients, model code and maps are freely available for use" —
  [Beggan et al. 2026](https://doi.org/10.1186/s40623-025-02360-0)
- NOAA's December 2025 State of the Geomagnetic Field report: global RMS errors of WMM and WMMHR rose by less than
  1% over the year and stay well below the specification — [NOAA repository 72269](https://repository.library.noaa.gov/view/noaa/72269)
  (search summary). The January/February 2026 NOAA "updated model" news is the annual High Definition Geomagnetic
  Model (HDGM) 2026, not a WMM update — [NCEI node 8611](https://www.ncei.noaa.gov/node/8611); [repost](https://www.guiceoffshore.com/updated-noaa-model-of-earths-geomagnetic-field-provides-more-accurate-and-safer-energy-exploration/)
- WMM error budget (Table 16, global RMS, combined error at 2025.0, WMM2025): X 133, Y 81, Z 128, H 130, F 124 nT,
  I 0.18°, D 0.36°, GV 0.70°; at 2030.0: X 140, Y 98, Z 154, H 138, F 138 nT, D 0.40°. Military specification
  MIL-PRF-89500B: H 200 nT, F 280 nT, I and D 1.00°. Crustal-field omission (127/78/125 nT in X/Y/Z) dominates —
  [WMM2025 Technical Report Table 16](https://doi.org/10.25923/prbc-s316)
- Blackout Zone where H < 2000 nT; Caution Zone where 2000 ≤ H < 6000 nT — [WMM2025 Technical Report §1.8](https://doi.org/10.25923/prbc-s316)
- Vertical validity (Table 17): e.g. F valid 0–3,500 km at G0 activity, 0–1,900 km at G2–G4, "Not advised" at G5;
  X "Exceeds" spec at G2 and above — [WMM2025 Technical Report Table 17](https://doi.org/10.25923/prbc-s316)

### Inferences
- **Discrepancy:** the README inside the WMM2025COF zip ("Date December 17, 2024") still quotes sample test values
  for 2020.0/2020.5 (e.g. `2020.0 66 14 143 0.12 13.08 34916.9 ...`), which are WMM2020-era numbers outside the
  WMM2025 validity; the WMM2025 file gives `2025.000000 66 14 143 -0.19 12.82 35003.64 ...` at the same place. Use the
  report Table 6 and the TestValues file, not the README sample.
- **Discrepancy:** the IAGA V-MOD page `ngdc.noaa.gov/IAGA/vmod/igrf.html` still showed only IGRF-13 when fetched
  ([page](https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html)), while `coeffs/igrf14coeffs.txt` on the same host is live.
  The NCEI IGRF product page and Zenodo are the paper's recommended sources.
- WMM2025 and IGRF-14 at 2025.0 differ by at most 28 nT in F on a 10°×15° surface grid (computed), consistent
  with both being averages of overlapping candidate models but truncated at different degrees (12 vs 13).
- A Rust crate can redistribute WMM2025.COF verbatim under NOAA's public-domain statement (with attribution,
  per the 17 U.S.C. 403 note) and IGRF-14 coefficients under the paper's "freely available" statement; MIT for the
  code is compatible with both.

### Gaps
- No official IGRF-14 point-value test table was found; the IGRF online calculators (NCEI, BGS, Kyoto) could be used
  to generate one but are calculators, not published references.
- BGS's own terms for the WMM (BGS is co-producer) and the Zenodo IGRF-14 record's licence field were not checked.

---

## 3. Dipole, poles, field decline, South Atlantic Anomaly, Lowes–Mauersberger spectrum and the CMB

### Takeaway
From IGRF-14, the dipole moment M = 4π a³ √(g₁⁰² + g₁¹² + h₁¹²)/μ₀ is 7.690×10²² A m² at 2025.0 (7.708×10²² in
2020, 8.321×10²² in 1900, 7.669×10²² forecast for 2030), i.e. falling ~4.8% per century now and ~6.1% per century on
average since 1900. The north geomagnetic pole is at 80.79°N geocentric (80.85°N geodetic), 72.76°W, tilt 9.21°;
the IGRF-14 SAA surface minimum (22,071 nT at 26.22°S, 60.03°W) and all model pole positions are reproduced. The
paper's 92.31% dipole share of Lowes power is not reproduced by the standard definition (92.86%).

### Cited Findings
- IGRF-14 defines M(t) = (4π/μ₀) a³ √(g₁⁰² + g₁¹² + h₁¹²) and the axial moment M_A = 4π a³ |g₁⁰|/μ₀; the dipole
  "is forecast to continue to fall in strength over the next five years"; Fig. 5 plots M in ZA m² (no values in the
  text) — [Beggan et al. 2026](https://doi.org/10.1186/s40623-025-02360-0)
- Geomagnetic (centred-dipole) poles are where the dipole axis meets "a sphere of mean Earth radius (6371.2 km)";
  dip poles are where I = 90° on the WGS-84 ellipsoid, found by gradient descent on a 0.01° grid "within 1 km
  accuracy" — [Beggan et al. 2026](https://doi.org/10.1186/s40623-025-02360-0)
- IGRF-14 Table 3, 2025.0: north dip pole 85.73°N, 138.65°E; south dip pole 63.86°S, 135.07°E; north geomagnetic
  pole 80.85°N, 72.76°W; south 80.85°S, 107.24°E. 2030.0: north dip 84.73°N, 125.94°E; south dip 63.66°S, 134.27°E;
  geomagnetic 81.05°N, 72.96°W. 1900.0: north dip 70.46°N, 96.19°W; geomagnetic 78.68°N, 68.79°W. 2020.0: north dip
  86.49°N, 162.80°E; geomagnetic 80.65°N, 72.68°W — [Beggan et al. 2026, Table 3](https://doi.org/10.1186/s40623-025-02360-0)
- North dip pole speed peaked at 56 km/yr in 2002.5, slowed to 40 km/yr by 2012.5, forecast 31 km/yr in 2027.5
  (~100 m/day); the paper says the two dip poles "will briefly have the same longitude (around 136°W) in early 2026"
  — [Beggan et al. 2026](https://doi.org/10.1186/s40623-025-02360-0)
- WMM2025, 2025.0: geomagnetic north pole 72.76°W, 80.79°N geocentric (80.85°N geodetic), tilt 9.21°; dip poles
  139.30°E/85.76°N (text) or 139.29°E (Table 4) and 135.08°E/63.85°S; WMMHR2025 (N = 15) 139.38°E/85.77°N and
  135.07°E/63.86°S; north dip pole drift ~35 km/yr, south ~9 km/yr over 2025–2030; eccentric dipole (first eight
  Gauss coefficients, Langel 1987 p. 386) at (r, φ′, λ) = (605 km, 22.76°N, 135.34°E) — [WMM2025 Technical Report §1.6, Table 4](https://doi.org/10.25923/prbc-s316)
- SAA (IGRF-14): deepened ~150 nT in five years, moving west ~20 km/yr; minimum surface F at 2025.0 is 22,071 nT at
  26.22°S, 60.03°W (NE Argentina); a second low offshore South Africa keeps growing. Largest dF/dt −128 nT/yr in North
  America, +104 nT/yr in the Indian Ocean — [Beggan et al. 2026](https://doi.org/10.1186/s40623-025-02360-0)
- SAA (WMM): F below 60% of average; Swarm GVO trends 2019–2024 show ~406 nT decrease near the Caribbean; at 450 km
  altitude the 2025 minimum is above 18,600 nT and an 18,600 nT contour appears by 2030 (deepening, broadening, a
  possible two-lobed structure tied to two reversed-flux patches at the CMB) — [WMM2025 Technical Report §1.7](https://doi.org/10.25923/prbc-s316)
- Lowes–Mauersberger spectrum (Mauersberger 1956; Lowes 1966) is the report's spectral quality tool; the WMMHR power
  falls from degree 1 to 15 (core) and flattens from 16 to 133 (crust), i.e. "degrees one to 15 ... dominated by
  sources in the Earth's core, while degrees 16 to 133 are dominated by" crustal sources — [WMM2025 Technical Report §2.3, Fig. 8](https://doi.org/10.25923/prbc-s316)
- IGRF-14: "the dipole terms g₁⁰, g₁¹, h₁¹ account for 92.31% of the power in the main geomagnetic field, a slight
  fall from 92.52% in 2020" (Lowes spectrum at the surface, 2025) — [Beggan et al. 2026](https://doi.org/10.1186/s40623-025-02360-0)
- LCS-1 lithospheric model spans degrees 14–185 (CHAMP 2006–2010 + Swarm A/C 2014–2016 gradients) — [search summary of LCS-1 documentation](https://uknowledge.uky.edu/ees_facpub/11); CHAOS-8 merges its static
  internal field with LCS-1 above degree 25 — [DTU geomagnetic field models](https://geomag.space.dtu.dk/english/models)

**Script results (computed from the cited coefficients; "ok" = matches the published value):**

| Quantity | Computed | Published | |
|---|---|---|---|
| WMM2025 geomagnetic N pole lat, geocentric / geodetic | 80.7898° / 80.8504° | 80.79° / 80.85° | ok |
| WMM2025 geomagnetic N pole lon, tilt | −72.7567°, 9.2102° | 72.76°W, 9.21° | ok |
| IGRF-14 2025 / 2020 geomagnetic N pole (geodetic lat, lon) | 80.8500°, −72.7628° / 80.6491°, −72.6774° | 80.85, −72.76 / 80.65, −72.68 | ok |
| WMM2025 eccentric dipole | 604.6 km, 22.756°N, 135.336°E | 605 km, 22.76°N, 135.34°E | ok |
| WMM2025 north / south dip pole | 85.7621°N 139.295°E / 63.8512°S 135.078°E | 85.76 139.29(.30) / 63.85 135.08 | ok |
| IGRF-14 2025 north / south dip pole | 85.7273°N 138.616°E / 63.858°S 135.063°E | 85.73 138.65 / 63.86 135.07 | ok, 0.41 km / 0.39 km apart |
| IGRF-14 2025 SAA min F at h = 0 | 22,071.4 nT at 26.2231°S, 60.0285°W | 22,071 nT, 26.22°S, 60.03°W | ok |
| WMM2025 SAA min F at h = 0 | 22,079.2 nT at 26.09°S, 60.14°W | — | computed only |
| IGRF-14 M 1900 / 2020 / 2025 / 2030 | 8.321 / 7.708 / 7.690 / 7.669 ×10²² A m² | (Fig. 5 only) | computed |
| IGRF-14 M_A 2025; B₀ 2025 | 7.591×10²² A m²; 29,733.4 nT | — | computed |
| Dipole share of Lowes power 2025 / 2020 | 92.86% / 93.04% | 92.31% / 92.52% | **not reproduced** |

IGRF-14 2025 Lowes spectrum (computed, nT²): surface R₁…R₁₃ = 1.768×10⁹, 8.53×10⁷, 3.90×10⁷, 9.02×10⁶, 2.06×10⁶,
3.16×10⁵, 1.62×10⁵, 2.58×10⁴, 1.61×10⁴, 3.47×10³, 750, 222, 128; downward-continued to r = 3485 km with
Rₙ(r) = (n+1)(a/r)²ⁿ⁺⁴ Σₘ(g² + h²): 6.60×10¹⁰, 1.07×10¹⁰, 1.63×10¹⁰, 1.26×10¹⁰, 9.6×10⁹, 4.9×10⁹, 8.4×10⁹, 4.5×10⁹,
9.4×10⁹, 6.7×10⁹, 4.9×10⁹, 4.8×10⁹, 9.3×10⁹ (roughly white for n ≥ 2). RMS field n = 1–13: 43,636 nT at the surface,
409,869 nT at the CMB.

### Inferences
- **Discrepancy:** IGRF-14 Table 3 lists the geomagnetic pole latitude as 80.85°, which is the geodetic latitude of
  the 80.79° geocentric point; the text describes the pole on a sphere. The WMM report states both explicitly.
- **Discrepancy:** WMM2025 north dip pole longitude is 139.30°E in the report text and 139.29°E in Table 4; the
  computed 139.295°E explains both roundings.
- **Likely typo:** the IGRF-14 paper's "same longitude (around 136°W) in early 2026" should be east: computing the
  IGRF-14 dip poles gives north 135.42°E and south 134.91°E at 2026.0, crossing near 134.9°E around 2026.2 (both poles
  are in the eastern hemisphere throughout 2025–2030 per Table 3).
- **Unresolved:** the 92.31% / 92.52% dipole shares are ~0.55 points below the standard Lowes ratio
  R₁/Σₙ₌₁¹³ Rₙ (92.86% / 93.04%); weights n, 2n+1, 1 or the WGS-84 radius do not give the paper's numbers either.
  The paper's definition is in the companion evaluation paper (Beggan et al. 2026, EPS 78:126), not checked.
- The brief's range "≈7.7–8.0×10²² A m²" brackets 1960s–2010s values; at 2025.0 the IGRF-14 moment is
  7.69×10²² A m² and still falling.
- Core–crust crossover: the surface IGRF spectrum falls to ~10² nT² by n = 12–13 while crustal models (LCS-1, MF7)
  begin at n ≈ 14–16, consistent with the WMM report's degree-15/16 split; the classic "n ≈ 14" figure is supported
  by that, though no source with an exact crossover degree was fetched.
- CMB continuation is the same (R/r)ⁿ⁺¹ upward/downward continuation the gravity engine already does; for B it is
  (a/c)ⁿ⁺² per coefficient. The CMB radius (3485 km here) is an assumption, not a fetched value.

### Gaps
- No primary source for the numerical decline rate (%/century) was fetched; the values above are computed from
  IGRF-14 Table 2.
- No source with Langel & Estes' crossover degree, or the exact CMB radius used by IGRF/CHAOS teams, was fetched.

---

## 4. Magnetic missions, model families and the link to GRACE-style gravity

### Takeaway
Swarm (since Nov 2013) is the backbone of WMM2025, IGRF-14 and CHAOS-8; IGRF-14 and CHAOS-8 are the first to use
MSS-1 (launched 21 May 2023, 41° inclination), and CHAOS-8 also uses CSES; CSES-02 launched 14 June 2025. The
gravity link is two-way: the same SH machinery evaluates both, and Swarm's GPS data also give monthly gravity fields
(through COST-G) that bridged the 2017–2018 GRACE/GRACE-FO gap at degrees below ~12.

### Cited Findings
- WMM2025 uses Swarm Level-1b data (version 0605 and later), "freely provided by ESA"; magnetospheric fields were
  handled with CHAOS-7.18 parameterized by the RC index; MF7 supplies WMMHR degrees 16–133 —
  [WMM2025 Technical Report §2](https://doi.org/10.25923/prbc-s316)
- IGRF-14 used Swarm, ground observatories, Cryosat-2, CSES and, for the first time, MSS-1 data —
  [Beggan et al. 2026](https://doi.org/10.1186/s40623-025-02360-0)
- MSS-1: launched 21 May 2023; inclination 41°, covering all local times up to ±41° within two months; altitude
  ~450 km in one source and 500 km in a meeting listing; data released to the community on 1 Aug 2024 after
  calibration ended 31 Jul 2024; a lithospheric model from MSS-1 + Swarm-A + CHAMP used CHAOS-8 to remove other
  sources — [EPP 2025 lithospheric paper](https://eppcgs.org/article/doi/10.26464/epp2025058);
  [EPP special-issue introduction](https://www.eppcgs.org/en/article/pdf/preview/10.26464/epp2025001.pdf);
  [Jiang et al. 2024 initial MSS-1 field model](https://ftp.space.dtu.dk/pub/cfinl/publications/2024_Jiang_Finlay_Olsen_Toeffner-Clausen_Zhang.pdf);
  [RAS MSS-1 meeting programme](https://upgrade.ras.ac.uk/sites/default/files/2023-10/Final%20programme-MSS-1%20meeting%281%29.pdf)
- CSES-01 launched Feb 2018; CSES-02 (Zhangheng 1-02) launched 14 June 2025 07:56 UTC from Jiuquan; its
  High-Precision Magnetometer (vector + scalar) was switched on 18 June and works normally; sun-synchronous, ~500 km;
  CNSA with ASI/INFN — [NSSC news](https://english.nssc.ac.cn/news/202506/t20250630_1046452.html); [CSES Mission](https://en.wikipedia.org/wiki/CSES_Mission)
- CHAOS-8 (Kloss, Finlay, Olsen, Tøffner-Clausen, Gillet, Grayver; EPS vol. 78, 2026): time-dependent internal field to
  degree 20, static internal field merged with LCS-1 above degree 25, magnetospheric field with induced counterpart,
  co-estimated ionospheric climatology, new temporal regularization from geodynamo secular-acceleration statistics;
  data from Swarm, MSS-1, CSES (first time in CHAOS), CryoSat-2, CHAMP, SAC-C, Ørsted and observatories. CHAOS-8.3
  (Zenodo) uses Swarm L1b to 27 May 2025, MSS-1A to 20 May 2025, CryoSat-2 to 27 May 2025 —
  [CHAOS-8 Zenodo](https://zenodo.org/records/15806240); [DTU models](https://geomag.space.dtu.dk/english/models);
  [EPS 2026 entry](https://doaj.org/article/a7d1b5ea44c7496fa53048cfd3e2bbbd). The Zenodo record gives a 1997–2025
  span, the paper 1999–2025.
- CHAOS-7 (Finlay et al. 2020) is the predecessor used by WMM2025 — [CHAOS-7 paper](https://ar5iv.labs.arxiv.org/html/2010.11313)
- Kalmag (GFZ): an EPS paper compares the Kalmag lithospheric model's spectrum and degree correlation with WDMAM
  and LCS-1 — [Kalmag figure, EPS 10.1186/s40623-022-01692-5](https://earth-planets-space.springeropen.com/articles/10.1186/s40623-022-01692-5/figures/7)
- Swarm gravity: monthly GPS-derived gravity fields from Swarm (data since Dec 2013) combined by COST-G to bridge the
  GRACE/GRACE-FO gap (July 2017–May 2018); the signal is largely limited to degrees below 12; ~4 cm EWH agreement with
  GRACE over land at 750 km smoothing and ~1 cm over ocean at 3000 km — [Teixeira da Encarnação et al. 2020, ESSD 12:1385](https://essd.copernicus.org/articles/12/1385/2020/)

### Inferences
- One Rust engine can serve both: gravity (ICGEM `.gfc`, 4π-normalized) and magnetics (Schmidt, mapped by
  √(2n+1)). Magnetic models need far lower degree (WMM 12, IGRF 13, CHAOS-8 core 20, LCS-1 185, WMMHR 133) than the
  crate's 2700 limit, so performance is not a concern.
- Swarm is a natural shared test case for this repo: its GPS gravity sits in the same COST-G system the project's
  GRACE work uses, and its magnetic data define the field models above.
- MSS-1's 41° inclination fills local-time coverage that polar Swarm orbits sample slowly, which is why core-field
  and ionospheric modellers adopted it quickly; this is inferred from the cited coverage statements.

### Gaps
- Kalmag's degree range and release details, LCS-1's original GJI paper, and Swarm/CHAMP mission fact sheets were not
  fetched directly; LCS-1 details come from a search summary of secondary pages.
- No primary data on the Swarm magnetometer error budget or on CHAOS-8 resolution limits beyond the degrees above.

---

## 5. Space weather: magnetopause, Parker spiral, Dst and DPS, storms, GICs, plasmapause and L-shells

### Takeaway
Pressure balance with the factor-2 field compression gives r_mp = (2B₀²/(μ₀ρv²))^(1/6) R_E ≈ 10.1 R_E for
n = 5 cm⁻³, v = 400 km/s (8.0 R_E without compression); the Parker angle at 1 au is ≈ 47° for 400 km/s; DPS gives
3.84×10¹³ J of ring-current energy per nT of Dst (3.8×10¹⁵ J at −100 nT, 1.6×10¹⁶ J for May 2024, 3.3–6.8×10¹⁶ J for
Carrington estimates); the plane-wave geoelectric field is |E| = √(ω/(μ₀σ)) |B| (not √(ωμ₀/σ) B, which needs H).

### Cited Findings
- Chapman–Ferraro balance: (ρv²)_sw ≈ (4B(r)²/2μ₀)_m, the factor 4 because the field just inside a planar
  magnetopause is about twice the dipole value; with B(r) = μ₀m/(4πr³), r ≈ (μ₀m²/(8π²ρv²))^(1/6); typical subsolar
  distances 6–15 R_E — [Magnetopause](https://en.wikipedia.org/wiki/Magnetopause)
- Shue et al. (1998): r₀ = (10.22 + 1.29 tanh(0.184(B_z + 8.14))) D_p^(−1/6.6) (nPa, nT, R_E), flaring
  α = (0.58 − 0.007B_z)(1 + 0.024 ln D_p); accurate to ~1 R_E for average driving but overestimates standoff when
  compressed below 10 R_E; one set of slides misprints 0.1849 — [EGUsphere 2025 preprint](https://egusphere.copernicus.org/preprints/2025/egusphere-2025-4530/egusphere-2025-4530.pdf);
  [EGU2020 presentation](https://presentations.copernicus.org/EGU2020/EGU2020-21977_presentation.pdf); [Shue et al. 1998 (ADS)](https://ui.adsabs.harvard.edu/abs/1998JGR...10317691S)
- IMF at Earth: nominally ~45° to the Sun–Earth line; 1–37 nT, averaging about 6 nT — [Interplanetary magnetic field](https://en.wikipedia.org/wiki/Interplanetary_magnetic_field)
- DPS theorem: −μ ΔB_C = 2U_K (μ = Earth's dipole moment, U_K = total kinetic energy of trapped plasma), independent
  of spatial distribution and spectrum; −Dst is a proxy for ΔB_C when Dst < 0 after removing induced currents —
  [Principles of Heliophysics](https://arxiv.org/pdf/1910.14022). Tests: AMPTE/CCE data from 80 storms support DPS on
  the nightside; ring-current ions supply on average about half of the Dst depression (Turner et al. 2001) — same
  search summary.
- May 2024 (Gannon) storm: Dst minimum −412 nT at 02:00 UT on 11 May (−406 nT in other reports); G5, Kp 9; largest
  Dst and most negative SYM-H since March 1989 — [May 2024 solar storms](https://en.wikipedia.org/wiki/May_2024_solar_storms);
  [arXiv 2508.16780](https://arxiv.org/pdf/2508.16780); [arXiv 2507.07009](https://arxiv.org/pdf/2507.07009)
- Carrington 1859 Dst estimates: −850 nT (Siscoe et al. 2006, hourly Colaba averages), −1760 nT (Tsurutani et al.
  2003; Love & Mursula 2024 could not reproduce it), −1150 nT (Gonzalez et al. 2011), −850 to −1150 nT (Cliver &
  Dietrich 2013), median −964 nT with 68% interval [855, 1087] nT (Love et al. 2024) — [SWSC 2026 article](https://www.swsc-journal.org/articles/swsc/full_html/2026/01/swsc250084/swsc250084.html);
  [NSF PAR copy](https://par.nsf.gov/servlets/purl/10480345)
- March 1989: storm arrived 01:27 UT 13 March; James Bay network offline in under 90 s; outage 9 hours; Dst −589 nT;
  Kp 9, Ap 246 — [March 1989 geomagnetic storm](https://en.wikipedia.org/wiki/March_1989_geomagnetic_storm);
  six million people without power for nine hours; GICs of tens to hundreds of amperes recorded in networks —
  [Geomagnetically induced current](https://en.wikipedia.org/wiki/Geomagnetically_induced_current). Boteler (2019):
  two CMEs (X4.5 flare 10 March, M7.3 on 12 March); the second sudden commencement coincided with the substorm that hit
  Hydro-Québec; the blackout came before Dst peaked — [Boteler 2019, Space Weather 17:1427](https://doi.org/10.1029/2019SW002278) (abstract via search)
- NERC TPL-007 benchmark: peak 8 V/km at 60°N geomagnetic latitude over a resistive (layered, 3–20,000 Ω m) ground
  model, scaled by latitude and conductivity factors, using the March 1989 Ottawa record; the 1-in-100-year estimate
  spanned 3–8 V/km; a supplemental event uses 12 V/km — [NERC webinar](https://nercstg.nerc.com/pa/Stand/Project201303GeomagneticDisturbanceMitigation/GMD_Stage2_Industry_Webinar_04242014.pdf);
  [Régie de l'énergie filing](https://www.regie-energie.qc.ca/fr/participants/dossiers/R-4192-2022/doc/R-4192-2022-A-0006-Autres-Autre-2022_09_08.pdf);
  [ISO-NE](https://www.iso-ne.com/static-assets/documents/2022/05/gmd_2026_na_april_2022_pac_rev1.pdf)
- Plasmapause: L_pp = 5.6 − 0.46 Kp_max (max Kp over the preceding 24 h), valid roughly 2.25 < L < 8, 00–15 MLT
  (Carpenter & Anderson 1992, JGR 97:1097, eq. 7); Moldwin et al. 2002 use a 12 h window —
  [rbamlib plasmapause models](https://rbamlib.readthedocs.io/dev/models/lpp.html)
- Dipole field line r = L cos²λ (planetary radii, geomagnetic latitude λ); Van Allen belts roughly L = 1.5–2.5 and
  4–6; plasmapause typically near L ≈ 5 — [L-shell](https://en.wikipedia.org/wiki/L-shell)

**Script results (computed, B₀ = 29,733 nT from IGRF-14 2025):**

| Quantity | Value |
|---|---|
| Dynamic pressure, n = 5 cm⁻³ protons, v = 400 km/s | 1.338 nPa |
| Chapman–Ferraro r_mp with compression factor 2 (field) / 4 (pressure) | 10.08 R_E |
| … without compression | 8.00 R_E (ratio 4^(1/6) = 1.260) |
| Wikipedia form with M_E = 7.690×10²² A m² | 10.08 R_E (identical) |
| Shue 1998 r₀ at B_z = 0, 1.338 nPa / at B_z = −20 nT, 20 nPa | 10.90 / 5.69 R_E |
| Parker angle tan ψ = Ωr/v, Ω = 2π/25.38 d, 1 au: 400 / 800 km/s | 46.98° / 28.18° |
| DPS: ΔB = −μ₀E/(2πB₀R_E³) → energy per −1 nT | 3.845×10¹³ J |
| Ring-current energy at Dst −100 / −412 / −589 / −850 / −1760 nT | 3.8×10¹⁵ / 1.58×10¹⁶ / 2.26×10¹⁶ / 3.27×10¹⁶ / 6.77×10¹⁶ J |
| L_pp at Kp_max 2 / 4 / 6 | 4.68 / 3.76 / 2.84 |
| Dipole: equatorial B at L = 4; invariant latitude of L = 4 | 464.6 nT; 60.0° |
| Dipole: B = B₀ L⁻³ √(1+3 sin²λ)/cos⁶λ at L = 4, λ = 30°; r there | 1456.8 nT; 3.0 R_E |
| Plane wave, 500 nT, T = 300 s, σ = 10⁻³ S/m | 2.04 V/km (skin depth 276 km) |
| Plane wave, 1000 nT, T = 300 s, σ = 10⁻⁴ S/m | 12.9 V/km (872 km) |
| Plane wave, 1000 nT, T = 60 s, σ = 10⁻⁴ S/m | 28.9 V/km (390 km) |
| MT identity ρ_a = 0.2 T |E/B|² (E mV/km, B nT) | recovers 1/σ exactly (check passes) |

### Inferences
- **Formula correction:** with B in tesla, the uniform half-space (Cagniard) result is E = Z H = √(iωμ₀/σ) B/μ₀, so
  |E| = √(ω/(μ₀σ)) |B| with E leading B by 45°. The brief's E = √(ωμ₀/σ) B is right only if "B" means H = B/μ₀.
  The check script verifies the equivalent magnetotelluric identity ρ_a = 0.2 T |E/B|².
- The plane-wave examples show that kV-scale line voltages follow from a few V/km over hundreds of km of resistive
  shield, the regime of the 1989 Hydro-Québec event and NERC's 8 V/km benchmark; layered-earth and 3-D effects need
  a full MT model rather than this half-space formula.
- The Carpenter & Anderson formula gives L_pp < 2 for Kp 9 (May 2024, March 1989), outside its fitted range; extreme
  storms need a different model.
- The 25.38-day sidereal solar rotation used for the Parker angle is an assumed standard value (not fetched); the
  result ≈ 47° is consistent with the cited "approximately 45°".

### Gaps
- Boteler (2019) full text (dB/dt and geoelectric values at Ottawa and in Québec) was behind a 403; no fetched
  source gives the 1989 dB/dt in nT/min.
- No source fetched for the full Dst index definition (Kyoto WDC) or the pressure-correction Dst* = Dst − b√P + c.
- No primary source for radiation-belt particle physics (adiabatic invariants, drift periods) was fetched.

---

## 6. Astrophysical fields: pulsars and magnetars, the Schwinger field, Faraday rotation, galactic and cluster fields, Biermann battery

### Takeaway
B_s = √(3c³I/(8π²R⁶)) √(PṖ) = 3.2×10¹⁹ √(PṖ) G (3.2×10¹⁵ √(PṖ) T) for I = 10⁴⁵ g cm², R = 10 km; it reproduces the ATNF
v2.6.1 surface fields of the Crab (3.79×10¹² G), SGR 1806−20 (2.06×10¹⁵ G, 47 B_Schwinger) and 1E 1547−5408
(2.22×10¹⁴ G). B_Schwinger = m_e²c²/(eħ) = 4.414×10⁹ T. The RM constant is e³/(8π²ε₀m_e²c³) = 2.631×10⁻¹³ rad m⁻² T⁻¹
= 0.8119 rad m⁻² per (cm⁻³ μG pc).

### Cited Findings
- Minimum surface field of a canonical pulsar (uniform sphere, M ≈ 1.4 M_⊙, R ≈ 10 km, I ≈ 10⁴⁵ g cm²):
  B > (3c³I/(8π²R⁶))^(1/2) (PṖ)^(1/2) → B(G) > 3.2×10¹⁹ (PṖ/s²)^(1/2), "sometimes called the characteristic magnetic
  field"; Crab P = 0.033 s, Ṗ = 10⁻¹²·⁴ gives B > 4×10¹² G — [NRAO Essential Radio Astronomy ch. 6](https://www.cv.nrao.edu/~sransom/web/Ch6.html)
- ATNF Pulsar Catalogue v2.6.1 (queried 2026-10-24): J0534+2200 P0 = 0.0333924123 s, P1 = 4.20972×10⁻¹³,
  Ė = 4.5×10³⁸ erg/s, BSURF = 3.79×10¹² G; J1808−2024 (SGR 1806−20) P0 = 7.55592 s, P1 = 5.49×10⁻¹⁰, Ė = 5.0×10³⁴,
  BSURF = 2.06×10¹⁵ G; J1550−5418 P0 = 2.06983302 s, P1 = 2.318×10⁻¹¹, Ė = 1.0×10³⁵, BSURF = 2.22×10¹⁴ G —
  [ATNF Pulsar Catalogue](https://www.atnf.csiro.au/research/pulsar/psrcat/) (cite Manchester et al. 2005,
  [AJ 129:1993](https://arxiv.org/abs/astro-ph/0412641)). The help page defines BSURF only as "Surface magnetic flux
  density (Gauss)" (formula shown as an image).
- McGill Magnetar Catalog (Olausen & Kaspi 2014): SGR 1806−20 P = 7.547728 s, Ṗ ≈ 49.5×10⁻¹¹, B ≈ 20×10¹⁴ G,
  age 0.24 kyr; older timing (P = 7.47 s, Ṗ = 8.3×10⁻¹¹) gave 8×10¹⁴ G; Ṗ ≈ 4.95×10⁻¹⁰ from Woods et al. 2007 —
  [McGill catalog](https://arxiv.org/pdf/1309.4167); [arXiv 2101.07923](https://arxiv.org/pdf/2101.07923)
- B_c = m_e²c²/(q_eħ) ≃ 4.41×10⁹ T; E_c = m_e²c³/(q_eħ) ≃ 1.32×10¹⁸ V/m — [Schwinger limit](https://en.wikipedia.org/wiki/Schwinger_limit)
- RM = e³/(8π²ε₀m²c³) ∫ n_e B_∥ ds ≈ (2.62×10⁻¹³ T⁻¹) ∫ n_e B_∥ ds (SI); Gaussian RM = e³/(2πm²c⁴) ∫ n_e B_∥ ds —
  [Faraday effect](https://en.wikipedia.org/wiki/Faraday_effect)
- Spiral galaxies: isotropic turbulent fields 20–30 μG in arms and bars, 50–100 μG in central starbursts, ordered
  fields 10–15 μG in interarm regions; large-scale reversals like the Milky Way's not yet seen in external galaxies —
  [Beck 2016, A&ARv 24:4](https://arxiv.org/abs/1509.04522)
- Milky Way, solar neighbourhood: regular ≈ 2 μG, total ≈ 6 μG, isotropic random ≈ 3–4 μG (Haverkorn 2014); a
  2026 review gives 6 ± 2 μG total near the Sun (Beck 2001 synchrotron), ~10 ± 3 μG at 3 kpc from the centre, and a
  pulsar-RM coherent field of 1.4 ± 0.2 μG, noting equipartition holds only on large scales —
  [Haverkorn 2014](https://arxiv.org/pdf/1406.0283); [arXiv 2606.02149](https://arxiv.org/pdf/2606.02149) (search summary)
- Galaxy clusters: "typical field strengths of order 1 microGauss with high areal filling factors out to Mpc
  radii"; 10–40 μG in cooling-flow cores; measured by RM of embedded/background sources, synchrotron halos/relics and
  inverse-Compton X-rays — [Carilli & Taylor 2002, ARA&A 40:319](https://arxiv.org/abs/astro-ph/0110655)
- Biermann battery (Biermann 1950, Z. Naturforsch. 5:65): ∂B/∂t = ∇×(v×B) + (∇T_e × ∇n_e)/(e n_e), non-zero when
  the density and temperature gradients are not collinear, so a field grows from zero; a candidate seed-field
  mechanism in protogalaxies — [Biermann battery](https://en.wikipedia.org/wiki/Biermann_battery)

**Script results:** prefactor √(3c³I/(8π²R⁶)) = 3.1996×10¹⁹ G (ok vs 3.2×10¹⁹); Crab 3.794×10¹² G, Ė 4.46×10³⁸ erg/s;
SGR 1806−20 2.061×10¹⁵ G = 2.06×10¹¹ T = 46.7 B_c, Ė 5.02×10³⁴; J1550−5418 2.217×10¹⁴ G = 5.0 B_c, Ė 1.03×10³⁵; Crab
0.086 B_c (all ok vs ATNF); B_c = 4.41401×10⁹ T (ok); RM constant 2.6312×10⁻¹³ (Wikipedia 2.62×10⁻¹³, within 0.5%)
and 0.81190 rad m⁻² in practical units (ok vs 0.812); B = √(2μ₀u) for u = 1 eV cm⁻³ is 6.35 μG.

### Inferences
- **Minor discrepancy:** the RM constant computed from CODATA 2022 is 2.631×10⁻¹³ rad m⁻² T⁻¹ (giving 0.8119);
  Wikipedia's 2.62×10⁻¹³ is slightly low. Use the computed value.
- The "3.2×10¹⁹" is the equatorial field of an orthogonal vacuum dipole; polar fields are twice that, and
  force-free magnetosphere spin-down formulas change the prefactor by order unity. ATNF's BSURF matches the vacuum
  3.2×10¹⁹ √(PṖ) exactly for all three test pulsars.
- SGR 1806−20's field depends on the timing epoch (8×10¹⁴ to 2.4×10¹⁵ G across papers); ATNF v2.6.1's 2.06×10¹⁵ G
  is a reproducible test value because P0 and P1 are listed with it.
- Equipartition with 1 eV cm⁻³ gives ~6 μG, matching the local total field, but the 1 eV cm⁻³ energy density is an
  assumed round number, not a fetched value. SI Biermann: ∂B/∂t = (k_B/e)(∇T_e × ∇n_e)/n_e with T_e in kelvin;
  Gaussian multiplies by c.

### Gaps
- Lorimer & Kramer's handbook was not fetched; the NRAO text stands in for it.
- No fetched numeric example of a Biermann seed field.

---

## 7. Planetary magnetic fields: Jupiter, Saturn, Mercury, Ganymede

### Takeaway
Juno's JRM33 gives Jupiter a dipole of 4.177 G R_J³ (a = 71,492 km), tilted 10.25° toward λ_III = 196.38°,
1.53×10²⁷ A m² ≈ 19,850 M_⊕ ("~20,000 times"). Saturn's Cassini 11+ g₁⁰ = 21,141 nT (R_S = 60,268 km, tilt < 0.007°)
gives 4.63×10²⁵ A m² ≈ 602 M_⊕. Mercury: 190 nT R_M³ (2.76×10¹⁹ A m²), offset 479 km north. Ganymede: 719 nT
equatorial (1.31×10²⁰ A m²).

### Cited Findings
- JRM33 (Connerney et al. 2022, JGR Planets, doi 10.1029/2021JE007055): 32 of the first 33 orbits, degree-30
  expansion, well constrained through degree 13 (information to ~18), dynamo radius ~0.81 R_J; "The degree 1
  coefficients of the JRM33 model describe a dipole with moment M = 4.177 G, offset from the rotation axis by
  θd = 10.25° towards System III longitude of ϕd = 196.38°"; vs JRM09: ΔM = 0.007 G, Δθd = 0.06°, Δϕd = 0.23°, and
  Δg₁⁰ = 0.18%, Δg₁¹ = 0.26%, Δh₁¹ = 1.7%; Jupiter's dipole moment is "some ∼20,000 times that of Earth's";
  1 R_J = 71,492 km — [Connerney et al. 2022 (Caltech record and PDF)](https://authors.library.caltech.edu/records/dds0z-2vx98)
- JRM09 (Connerney et al. 2018) is degree/order 10 from the first nine orbits — [EGU2018 abstract](https://meetingorganizer.copernicus.org/EGU2018/EGU2018-7868.pdf?pdf);
  the non-dipolar field is concentrated in the northern hemisphere — [Moore et al. 2018, Nature](https://ideas.repec.org/a/nat/nature/v561y2018i7721d10.1038_s41586-018-0468-5.html)
- NASA fact sheets (older models): Jupiter GSFC-O6 dipole 4.30 G Rj³, tilt 9.4°, offset 0.119 Rj, Rj = 71,398 km,
  surface 4.0–13.0 G — [Jupiter fact sheet](https://nssdc.gsfc.nasa.gov/planetary/factsheet/jupiterfact.html); Saturn
  GSFC-Z3 0.215 G Rs³, tilt 0.0°, offset 0.038 Rs, Rs = 60,330 km, surface 0.18–0.84 G —
  [Saturn fact sheet](https://nssdc.gsfc.nasa.gov/planetary/factsheet/saturnfact.html); Mercury "MESSENGER Model 2010"
  0.002 G Rh³, offset 0.17 Rh, Rh = 2,440 km, surface 0.0025–0.007 G — [Mercury fact sheet](https://nssdc.gsfc.nasa.gov/planetary/factsheet/mercuryfact.html)
- Saturn: Cassini 11 (Dougherty et al. 2018) g₁⁰ = 21140, g₂⁰ = 1581, g₃⁰ = 2260 nT; Cassini 11+ (Cao et al. 2020)
  g₁⁰ = 21141, g₂⁰ = 1583, g₃⁰ = 2262 nT; reference 1-bar spheroid a = 60,268 km; dipole tilt < 0.007° (25.2″);
  magnetic equator ~0.0468 R_S north (2820 ± 26 km vs 2808.5 ± 12 km in Dougherty 2018) —
  [Cao et al. 2020, Icarus, arXiv 1911.06952](https://arxiv.org/abs/1911.06952) (Table 5); [Dougherty et al. 2018](https://resolver.caltech.edu/CaltechAUTHORS:20181005-093938794)
- Mercury (MESSENGER): Anderson et al. 2012: low-altitude equator crossings offset 479 ± 6 km north (486 ± 74 km at
  high altitude); offset dipole with southward moment 190 nT R_M³ fits with RMS < 14 nT; Anderson et al. 2011:
  195 ± 10 nT R_M³, 484 ± 11 km, tilt < 3° — [Anderson et al. 2012](https://scholars.unh.edu/cmerg/107);
  [EGU2012 abstract](https://meetingorganizer.copernicus.org/EGU2012/EGU2012-9245.pdf); [BepiColombo review](https://arxiv.org/pdf/2202.13243)
- Ganymede (Kivelson et al. 2002): permanent dipole with equatorial surface field 719 nT, tilted 176° from the spin
  axis, though the authors favour an inductive-response fit; the 1996 discovery gave ~750 nT —
  [Kivelson et al. 2002 (NTRS)](https://ntrs.nasa.gov/citations/20020044825); irfpy default moment −1.31×10²⁰ for
  719 nT — [irfpy Ganymede field](https://irfpy.irf.se/projects/pep/api/api_irfpy.jupsci.ganymede_magfield.html)

**Script results (M = 4πa³B₀/μ₀):** Earth 2025 7.690×10²² A m²; Jupiter 1.526×10²⁷ (19,849 M_⊕; ok vs "~20,000");
Saturn 4.628×10²⁵ (602 M_⊕); Mercury 2.760×10¹⁹ (3.6×10⁻⁴ M_⊕); Ganymede (R = 2631.2 km) 1.310×10²⁰ (ok vs 1.31×10²⁰).

### Inferences
- JRM09's dipole is 4.170 G (from JRM33's 4.177 G and the quoted ΔM = 0.007 G) and its tilt 10.31° (10.25° + 0.06°).
- The NASA fact-sheet values are pre-Juno/pre-Grand-Finale models with different reference radii (71,398 km,
  60,330 km); tests should use JRM33 and Cassini 11+ with their own radii (71,492 km and 60,268 km). The Saturn moment
  in G R_S³ from Cassini 11+ is 0.2114 G at 60,268 km, close to the 0.215 G sheet value.
- Mercury's 190 vs 195 nT R_M³ and 479 vs 484 km pairs come from the 2012 and 2011 Anderson papers respectively and
  should not be mixed.
- The JRM33 g₁⁰, g₁¹, h₁¹ values themselves are in the machine-readable supplement and PDS (JRM33.mod), which would let
  the same Schmidt engine evaluate Jupiter's field (degree 30); the Ganymede radius 2631.2 km is assumed.

### Gaps
- JRM33's individual Gauss coefficients and a Juno published test value (field at a given point) were not obtained.
- Uranus and Neptune (Voyager) and Mercury's BepiColombo 2026 orbit-insertion results were outside what was fetched.
