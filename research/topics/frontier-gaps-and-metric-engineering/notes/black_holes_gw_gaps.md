# Five gaps in the black-hole and gravitational-wave code: methods, sources, test values

Scope: `crates/gravitational/src/{inspiral,ringdown,remnant,energy,tides}.rs` and
`crates/gravitics/src/spacetime.rs`. Every test value below is quoted from a source I opened (named
with equation/table), or is marked **recomputed** (my Python, scripts in `bh_gw_scripts/`, run with
`python3 -P`). Downloads are in `bh_gw_sources/<id>/` (arXiv TeX sources, two PDFs, Berti's archived
data page and notebook, the `qnm` README and source read as text). No downloaded code was executed.

Scripts (all in `bh_gw_scripts/`):
- `gap1_kerr_isco_frequency.py`: Kerr ISCO frequency, limits, Isoyama's redshift identity, f_GW table.
- `gap2_leaver_kerr.py`: Leaver radial continued fraction + Cook–Zalutskiy spectral angular solver;
  checks against Berti's Schwarzschild file, the `qnm` README and BCW06 tables.
- `gap2_angular_cf_check.py`, `gap2_cf_only_robustness.py`: Leaver's angular continued fraction
  against the spectral method, and why it fails without care.
- `gap2_all_modes_from_fits.py`, `gap2_fit_errors.py` (+ `gap2_fit_errors.out`): all 63 BCW06 modes
  polished from the fits; true errors of the fits.
- `gap2_inverse_iteration.py`: the angular constant by Rayleigh-quotient inverse iteration (no eigen-solver).
- `gap3_remnant_precessing.py`: HBR16 generic-spin formula, BMR12 generic E_rad, Kerr κ and first law.
- `gap4_shapiro_nlo.py`: lensing and retardation terms for the Double Pulsar, with an independent
  vector check of KS99 Eq. 130.
- `gap5_love_tov.py`, `gap5_checks2.py`, `gap5_checks3.py`: TOV + Love-number solver, polytropes,
  incompressible stars, Read et al. piecewise polytropes (SLy, AP4), cancellation and series checks.

## Gap 1: Kerr ISCO gravitational-wave frequency

### Takeaway
Use M Ω = 1/((r_ISCO/M)^{3/2} + a) with the library's signed-spin `kerr_isco` radius, and
f_GW = Ω/π. This is printed in Isoyama et al. 2014, Eq. 3, which cites BPT72. It needs three lines in
`inspiral.rs`, with exact tests at a = 0, ±1. Do not add the conservative self-force shift: it is an
O(μ/M) correction meant for extreme mass ratios. At comparable masses it would be a 25–35% "correction"
far outside its range of validity.

### Cited Findings
- Isoyama et al. 2014 print the test-mass ISCO radius in BPT form, r/M = 3 + Z₂ ∓ [(3 − Z₁)(3 + Z₁ + 2Z₂)]^{1/2}, with Z₁, Z₂ as in the library. They give the orbital frequency in t as M Ω⁽⁰⁾_ISCO = [(r⁽⁰⁾_ISCO/M)^{3/2} + q]^{-1}, with q = S/M² > 0 prograde and < 0 retrograde (their Eq. 3; a commented-out line in the TeX shows the ± form ±M^{1/2}/(r^{3/2} ± q M^{3/2})) — [Isoyama et al. 2014, arXiv:1404.6133](https://arxiv.org/abs/1404.6133)
- Isoyama et al. also give the circular-orbit redshift z₀² = (1 − aΩ)[1 + aΩ − 3(MΩ)^{2/3}(1 − aΩ)^{1/3}]. The ISCO condition is z₀''(Ω_ISCO) = 0, and E₀ = z₀ − Ω z₀' — [Isoyama et al. 2014](https://arxiv.org/abs/1404.6133)
- They parametrise the conservative self-force shift as (M + μ)Ω_ISCO = MΩ⁽⁰⁾_ISCO(q)[1 + η C_Ω(q) + O(η²)], with η := μ/M. Table I gives C_Ω = 1.245568(2) at q = 0.1, 1.239927(4) at 0.3, 1.249234(6) at 0.5, 1.29319(1) at 0.7, 1.38157(3) at 0.9, and 1.286093(9) at −0.5, 1.314016(9) at −0.9. The text gives C_Ω(0) = 1.25101539 ± 4×10⁻⁸, with a minimum near q ≈ 0.3 — [Isoyama et al. 2014](https://arxiv.org/abs/1404.6133)
- The same paper says the sharp ISCO is replaced by a transition regime for a self-gravitating body. A useful ISCO survives only for the conservative dynamics — [Isoyama et al. 2014](https://arxiv.org/abs/1404.6133)
- Barack & Sago 2009 give Δr_ISCO = −3.269(±0.003)μ and ΔΩ_ISCO/Ω_ISCO = 0.4870(±0.0006) μ/M (their Eq. 300), and state that these are Lorenz-gauge values — [Barack & Sago 2009, arXiv:0902.0573](https://arxiv.org/abs/0902.0573)
- **Recomputed** (`gap1_kerr_isco_frequency.py`). The redshift identity holds at every spin tested (z₀'' and E₀' vanish to finite-difference precision ~1e-5), and E₀(Ω_ISCO) equals √(1 − 2/(3r)) to 1e-10. Exact limits: MΩ = 6^{-3/2} = 0.068041381743977 (a = 0), 1/2 (a = +1, r = 1), 1/26 (a = −1, r = 9). Values: a = 0.9 gives r = 2.3208830418, MΩ = 0.225441708835; a = 0.99 gives r = 1.4544979381, MΩ = 0.364409991055. For M = 10 M☉, f_GW = 439.704 Hz (a = 0), 1456.873 Hz (0.9), 2354.929 Hz (0.99), 3231.153 Hz (1), 248.550 Hz (−1). The constant-free ratio f(a)/f(0) is 3.31330 (0.9), 5.35571 (0.99), 6^{3/2}/2 = 7.348469228 (1) and 6^{3/2}/26 = 0.565266864 (−1) — [script](bh_gw_scripts/gap1_kerr_isco_frequency.py)

### Inferences
- BS09's 0.4870 is not comparable with C_Ω. With η = μ/M, (M+μ)Ω = MΩ₀(1 + μ/M)(1 + 0.4870 μ/M) implies C_Ω ≈ 1.487, but Isoyama et al. give 1.2510. The 0.236 difference is presumably the Lorenz-gauge time normalisation. I did not open the papers that correct it. Any C_Ω test should use Isoyama's numbers only.
- Implementation: `pub fn kerr_f_isco(m: f64, a: f64) -> f64 { let (r, _, _) = kerr_isco(a); C.powi(3) / (PI * G * m * (r.powf(1.5) + a)) }`. Tests: equals `f_isco(m)` at a = 0; ratio to `f_isco` is 6^{3/2}/2 at a = 1 and 6^{3/2}/26 at a = −1. The existing test at inspiral.rs:139 already pins the a = 0 case of the same formula. Cost: about 5 lines plus 5 test lines. Risk: the a = 1 limit relies on `cbrt(0)`, which is fine. The convention f_GW = 2 f_orb (the l = m = 2 harmonic) must be documented.
- `inspiral::horizon` uses `f_isco(m1 + m2)` as the upper cut-off. Leave it: switching it to a Kerr cut-off would need a remnant or effective spin, and nobody asked for that.

### Gaps
- BPT 1972 (ApJ 178, 347) is not on arXiv, and I did not open it. The formulas are quoted as reprinted by Isoyama et al. and HBR16.
- I found no primary source printing r_ISCO or Ω_ISCO at a = 0.9 or 0.99. A web search found none: [search hits included arXiv:1503.07060](https://arxiv.org/abs/1503.07060v2) on spinning particles, which is not the geodesic case. The a = 0.9/0.99 values above are recomputed only.

## Gap 2: Exact Kerr quasinormal modes

### Takeaway
The most accurate route that a few hundred lines of Rust can hold is Leaver's radial continued fraction, in the Cook–Zalutskiy confluent-Heun form with Nollert's tail. Pair it with the Cook–Zalutskiy spectral (spherical-harmonic) eigenproblem for the angular constant A, re-solved at every trial ω, and start the complex secant from the library's own BCW06 fit.

My Python version of exactly this reproduces several sources:
- Berti's Schwarzschild l = 2 overtones n = 0–5 to ≤ 8×10⁻¹⁶;
- the `qnm` package's a = 0.68 value to 4×10⁻¹⁶;
- all 693 BCW06 table entries to their 4 decimals, apart from two misprints that I traced to BCW06 itself.

It converges for all 63 modes up to j = 0.99 with no spin tracking. Leaver's angular continued fraction without inversion picks the wrong l in 36 of 315 cases, so use the spectral matrix, with inverse iteration in place of a full eigensolver.

### Cited Findings
- Cook & Zalutskiy 2014 solve the radial Teukolsky equation with an independent high-precision Leaver continued fraction and the angular equation with a spectral eigenvalue method. They used 24-digit (up to 32) arithmetic, absolute tolerance 1e-12, and conservatively claim at least ±1e-10 (typically ~1e-12). The matrix size is at least 4 and the continued fraction length at least 300, both chosen dynamically. The inversion number is n. Steps in a are ≤ 1e-3, down to 1e-9 near a = 1 — [Cook & Zalutskiy 2014, arXiv:1410.7698](https://arxiv.org/abs/1410.7698)
- Radial equations, M = 1, s = −2 (Cook & Zalutskiy Sec. II C). Each step:
  - Horizon frequencies: r± = 1 ± √(1−a²), σ± = (2ω r± − m a)/(r₊ − r₋).
  - Characteristic exponents, choice (ζ₊, ξ₋, η₊): ζ = iω, ξ = −s − iσ₊, η = −iσ₋.
  - Heun parameters: p = (r₊ − r₋)ζ/2; α = 1 + s + ξ + η − 2ζ + s iω/ζ; γ = 1 + s + 2η; δ = 1 + s + 2ξ.
  - Accessory parameter: σ = A + a²ω² − 8ω² + p(2α + γ − δ) + (1 + s − (γ+δ)/2)(s + (γ+δ)/2).
  - D-coefficients: D₀ = δ, D₁ = 4p − 2α + γ − δ − 2, D₂ = 2α − γ + 2, D₃ = α(4p − δ) − σ, D₄ = α(α − γ + 1).
  - Recurrence: αₙ = n² + (D₀+1)n + D₀, βₙ = −2n² + (D₁+2)n + D₃, γₙ = n² + (D₂−3)n + D₄ − D₂ + 2.
  - QNM condition: Cf(n;N) = 0, the n-th inversion of β₀ − α₀γ₁/(β₁ − α₁γ₂/(β₂ − …)), with tail r_N ≈ 1 + u₁/√N + u₂/N + u₃/N^{3/2}.
  - Tail terms: u₁ = ±√(−4p) (minimal branch), u₂ = −(8p − 4α + 2γ + 2δ + 3)/4, u₃ = [32p(2p − 4α + γ + 3δ + 4) + 4(γ+δ)(γ+δ−2) + 16σ + 3]/(32u₁).

  — [Cook & Zalutskiy 2014](https://arxiv.org/abs/1410.7698)
- Angular equation (Cook & Zalutskiy Eqs. 51–55). The pentadiagonal matrix is built from:
  - recurrence coefficients F_{sℓm} = √(((ℓ+1)²−m²)/((2ℓ+3)(2ℓ+1)))·√(((ℓ+1)²−s²)/(ℓ+1)²), G_{sℓm} = √((ℓ²−m²)/(4ℓ²−1))·√((ℓ²−s²)/ℓ²) (0 at ℓ = 0), and H_{sℓm} = −ms/(ℓ(ℓ+1));
  - products 𝒜 = F_ℓF_{ℓ+1}, 𝒟 = F_ℓ(H_{ℓ+1}+H_ℓ), ℬ = F_ℓG_{ℓ+1} + G_ℓF_{ℓ−1} + H_ℓ², ℰ = G_ℓ(H_{ℓ−1}+H_ℓ), 𝒞 = G_ℓG_{ℓ−1};
  - diagonal ℓ'(ℓ'+1) − s(s+1) − c²ℬ + 2csH, and off-diagonals −c²𝒜 (ℓ' = ℓ−2), −c²𝒟 + 2csF (ℓ−1), −c²ℰ + 2csG (ℓ+1), −c²𝒞 (ℓ+2), with c = aω.

  A is the eigenvalue continuous with ℓ(ℓ+1) − s(s+1) — [Cook & Zalutskiy 2014](https://arxiv.org/abs/1410.7698)
- Cook & Zalutskiy found that alternating radial and angular solves can converge slowly or fail. Their fix re-solves the angular problem at every ω used in the Newton finite differences — [Cook & Zalutskiy 2014](https://arxiv.org/abs/1410.7698)
- Stein's `qnm` package (JOSS 2019) implements the Leaver radial solver with the Cook–Zalutskiy spectral angular sector, root-polishing a single function E_r(ω) = E_r(ω, A(ω)) in 2-D. It interpolates initial guesses from a cached table — [Stein 2019, arXiv:1908.10377](https://arxiv.org/abs/1908.10377). Its source (read, not run) uses the characteristic exponents (ζ₊, ξ₋, η₊). A docstring there says the inner-horizon choice does not affect the minimal solution — [qnm radial.py](https://github.com/duetosymmetry/qnm/blob/master/qnm/radial.py)
- `qnm` README test value: (s, l, m, n) = (−2, 2, 2, 0) at a = 0.68 gives Mω = 0.5239751042900845 − 0.08151262363119974i. The same number is asserted in its test_qnm.py — [qnm README](https://github.com/duetosymmetry/qnm)
- Berti's Schwarzschild file s2l2.dat lists Re, Im (units 2M = 1), error and overtone n. Halved to Mω:

  | n | Mω |
  |---|---|
  | 0 | 0.3736716844180419 − 0.0889623156889357i (error 1.7e-13) |
  | 1 | 0.3467109968791634 − 0.2739148752912349i |
  | 2 | 0.3010534546123664 − 0.4782769832230718i |
  | 3 | 0.2515049621855907 − 0.7051482024334955i |
  | 4 | 0.2075145798130661 − 0.9468448908663515i |
  | 5 | 0.1692994030930433 − 1.195608054135847i |
  | 8 | ≈ −2i (Re 4.6e-16, Im 3.9999999999996 in 2M units: the algebraically special mode) |

  — [Berti ringdown page, archived 2024-03-18](https://web.archive.org/web/20240318065031/https://pages.jh.edu/~eberti2/ringdown/) and its [s2l2.dat](https://web.archive.org/web/20240318065031/https://pages.jh.edu/~eberti2/ringdown/s2l2.dat). This confirms the user's 0.37367168 − 0.08896232i.
- Berti's page warns that his Kerr frequencies are unreliable very close to extremality, pointing to arXiv:1212.3271, 1307.8086 and 1410.7698. His Kerr data are SharePoint links — [Berti ringdown page (archive)](https://web.archive.org/web/20240318065031/https://pages.jh.edu/~eberti2/ringdown/)
- Berti's Leaver notebook (Leaver units M = 1/2) gives Leaver's angular continued fraction:
  - k₁ = |m−s|/2, k₂ = |m+s|/2;
  - αₙ = −2(n+1)(n+2k₁+1);
  - βₙ = n(n−1) + 2n(k₁+k₂+1−2aω) − [2aω(2k₁+s+1) − (k₁+k₂)(k₁+k₂+1)] − [a²ω² + s(s+1) + A];
  - γₙ = 2aω(n+k₁+k₂+s).

  The notebook uses 130 terms and seven radial/angular iterations at a/M = 0.4 — [QNMLeaverKerr.nb (archive)](https://web.archive.org/web/20170210200510/http://www.phy.olemiss.edu:80/~berti/ringdown/QNMLeaverKerr.nb)
- BCW06 Tables II–IV give Mω to 4 decimals at j = 0, 0.1, …, 0.9, 0.98. Their fits are F = f₁ + f₂(1−j)^{f₃} and Q = q₁ + q₂(1−j)^{q₃}, with maximum errors over 0 ≤ j ≤ 0.99 in Tables VIII–X — [Berti, Cardoso & Will 2006, gr-qc/0512160](https://arxiv.org/abs/gr-qc/0512160) (PDF Tables II, IV, X read)
- BCW06 misprints found by the exact solver, checked in the PDF:
  - Table II (2, −2, 0) at j = 0.60 prints Im = .0890; the exact value is 0.0888917, i.e. .0889.
  - Table IV (4, 2, 0) at j = 0.00 prints Im = .0941; the Schwarzschild value is 0.0941640, i.e. .0942, as the other eight m columns print.
  - Table II (2, −2, 0) at j = 0.98 prints ".2927-.0881" with a stray minus sign.

  The library's `TABLES` string copies the first two faithfully — [BCW06 PDF](https://arxiv.org/abs/gr-qc/0512160); [script](bh_gw_scripts/gap2_all_modes_from_fits.py)
- **Recomputed, matching the sources above.**
  - Berti n = 0–5 agree to ≤ 8e-16; n = 6 to 4e-13 and n = 7 to 1e-10 with N = 1500.
  - qnm a = 0.68 agrees to 3.9e-16.
  - (2,2,0) at a = 0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 0.98, 0.99: 0.387017538366 − 0.088705699028i; 0.419526681764 − 0.087729271894i; 0.464123025976 − 0.085638834988i; 0.532600243551 − 0.080792873155i; 0.586016974909 − 0.075629552356i; 0.671614272132 − 0.064869235876i; 0.825429477482 − 0.038630197156i; 0.870892658736 − 0.029390424219i.
  - Angular constants at a = 0.9 and 0.99: A = 2.1098195941 + 0.2111235789i and 1.1196449289 + 0.1179539201i.
  - Others: (2,2,1) at a = 0.9: 0.6676575508 − 0.1952520671i. (2,1,0) at 0.9: 0.5162914042 − 0.0698043456i. (2,−2,0) at 0.9: 0.2972442136 − 0.0882806593i. (3,3,0) at 0.9: 1.0446370928 − 0.0654629044i. (4,4,0) at 0.9: 1.4104161198 − 0.0661559808i.

  All round to the BCW06 table entries — [script](bh_gw_scripts/gap2_leaver_kerr.py)
- **Recomputed.** The spectral and Leaver angular continued fractions agree to ≤ 7e-15 at a = 0.4–0.99 when started from the right A. Started blindly from ℓ(ℓ+1) − s(s+1) with no inversion, Leaver's angular continued fraction converged to the ℓ = 2 root for 36 of 315 cases, all (3,−2,n), (4,−2,n), (4,−1,n) and (4,2,n), so ω jumped to the wrong mode. The spectral "closest eigenvalue" choice never failed — [script](bh_gw_scripts/gap2_cf_only_robustness.py)
- **Recomputed: true errors of the BCW06 fits.** Exact Leaver values were taken at j = 0, 0.05, …, 0.95, 0.98, 0.99 for all 63 modes.
  - Every frequency fit's maximum error equals the stated F-% to the printed 2 decimals, which independently confirms the solver.
  - 13 quality-factor fits exceed the stated Q-% by more than 0.02 points. In nine l = 4 modes, (4,2,0–2), (4,1,0–2) and (4,0,0–2), the excess is exactly +2.00 or +3.00 points: stated 0.26 against true 2.26, 0.12 against 3.11, 0.31 against 3.31, and so on.
  - The others: (2,−1,1) is 3.22% at j = 0.6 against a stated 0.67; (3,3,0) 0.59 against 0.42; (4,4,1) 1.02 against 0.77; (4,4,2) 1.09 against 0.10.

  — [output](bh_gw_scripts/gap2_fit_errors.out)

### Inferences
- The exact +2.00/+3.00 excesses make it likely that BCW06 Table X's Q-error column lost a leading digit for l = 4, m = 0–2. The fits themselves are probably as intended, and the library's "eleven fits miss by up to 3.3%" is mostly a misprint in the stated error, not a bad fit. I checked that the PDF prints the same numbers as the library. The (2,−1,1) case (a peak at j = 0.6) may instead be a coefficient misprint, e.g. q₃ = 13.1683.
- Rust design, a ponytail-sized version (roughly 150–250 lines plus tests, no new dependency):
  1. a 40-line `Copy` complex type (add, sub, mul, div, sqrt);
  2. `fn leaver_cf(ω, a, l, m, n, A, N) -> Complex` as above, with Nollert's tail through u₃ and N ≈ 300–800;
  3. `fn sep_const(c, l, m, guess) -> Complex`, by inverse iteration with shift `guess` on the (l_max − l_min + 1)² ≈ 20×20 pentadiagonal matrix (banded complex Gaussian elimination, a few iterations, eigenvalue from the bilinear Rayleigh quotient vᵀMv/vᵀv, valid because the matrix is complex symmetric). **Recomputed:** 8 Rayleigh-quotient iterations started from v = e_l and λ = l(l+1) − 2 reproduce numpy's closest eigenvalue to ≤ 3.6e-14 for (2,2,0) a = 0.99, (3,−2,0) 0.9, (4,2,1) 0.98, (4,−1,2) 0.7 and (2,−2,0) 0.98, including the l = 3, 4 cases where Leaver's uninverted angular continued fraction failed ([script](bh_gw_scripts/gap2_inverse_iteration.py));
  4. complex secant on f(ω) = Cf(n;N; ω, A(aω)), carrying A from the previous evaluation, as Cook & Zalutskiy do;
  5. initial ω from the existing `qnm()` fit, ω₀ = F − iF/(2Q).

  This needs no spin tracking: all 63 modes × 11 spins converged from the fit guesses, and 315 more up to 0.99. It could live in `ringdown.rs` as `qnm_exact(l, m, n, j) -> Option<(f64, f64)>`, returning the same (F, Q) pair as `qnm`.
- Tests to add:
  - Berti Schwarzschild n = 0–3 (|Δ| < 1e-12);
  - the `qnm` README value at 0.68 (< 1e-12);
  - all BCW06 table entries within 6e-5, after correcting the two misprints in the test string;
  - Cf symmetry under ω → −ω*, m → −m, A → A*. Cook & Zalutskiy give A_{s,ℓ,−m}(c) = A_{s,ℓ,m}(−c) and A*(c) = A(c*). `qnm`'s test_qnm.py asserts the mirror mode −ω* at −m.
- Pitfalls:
  - (i) the sign of u₁ (take Re u₁ < 0, the minimal solution);
  - (ii) BCW06's m < 0 are the positive-Re ω counter-rotating modes, i.e. track the root with Re ω > 0 for negative m;
  - (iii) overtones n ≥ 5 at l = 2 approach the negative imaginary axis and the algebraically special −2i (Cook & Zalutskiy discuss this at length), so test only n ≤ 3 and expect N to grow;
  - (iv) j > 0.99: zero-damped modes need N growing like 1/√(1−a²) and tiny secant steps, and Berti flags his own data there;
  - (v) Leaver's angular continued fraction needs inversion number l − l_min to pick the right root; the spectral matrix avoids this;
  - (vi) units: Leaver and Berti use 2M = 1 (ω doubled, a halved), Cook and `qnm` use M = 1;
  - (vii) the sign convention e^{−iωt} gives Im ω < 0.
- Interpolating BCW06's tables is strictly worse: only 4 decimals at 11 spins. The fits already match the tables at that level.

### Gaps
- Leaver 1985 (Proc. R. Soc. A 402, 285) is not on arXiv, and I did not open it. Its formulas are taken from Cook & Zalutskiy's re-derivation and Berti's notebook.
- I could not download Cook's high-precision Kerr data (Zenodo record 2650358, KerrQNM_nn.h5): Zenodo returned 403 for unusual traffic, and the Wake Forest mirror is behind a Cloudflare challenge, which I did not attempt. Berti's Kerr tables are on SharePoint links I did not fetch. So the only primary 8+ digit Kerr test value I hold is the `qnm` README one. The other Kerr digits above are recomputed, cross-checked only to 4 decimals (BCW06) and by the internal identities. The `qnm` precomputed cache is a pickle, which I deliberately did not download or load.

## Gap 3: Precessing remnant spin, generic radiated energy, Kerr Hawking temperature

### Takeaway
HBR16's generic formula is a 20-line addition to `remnant.rs` that reuses the existing aligned pieces. BMR12's generic radiated energy needs no new code: it is the existing `radiated_energy` with a_i replaced by the projections |a_i| cos θ_i. The NR surrogate NRSur7dq4Remnant is about 10× more accurate. But it is a Gaussian-process model trained on 1528 simulations and needs its data file, so it is out of scope for a few hundred lines. The Kerr temperature is a 3-line function in `gravitics::spacetime` with exact rational tests.

### Cited Findings
- HBR16 generic spin. The remnant spin vector is a_fin = a_tot + ℓν, with a_tot = (a₁ + a₂q²)/(1+q)², q = M₂/M₁ ≤ 1 and ν = q/(1+q)². The magnitude of ℓ is |ℓ| = |L_ISCO(a_eff) − 2a_tot(E_ISCO(a_eff) − 1) + Σ_{i≤n_M} Σ_{j≤n_J} k_ij ν^{1+i} a_eff^j| (Eq. ell; the TeX writes the upper j limit as n_S). The projections that enter are a_tot = (|a₁|cos β + |a₂|cos γ q²)/(1+q)² and a_eff = a_tot + ξν(|a₁|cos β + |a₂|cos γ) (Eqs. atot, aeff). The magnitude is then |a_fin| = (1+q)^{-2}[|a₁|² + |a₂|²q⁴ + 2|a₁||a₂|q² cos α + 2(|a₁|cos β + |a₂|q² cos γ)|ℓ| q + |ℓ|² q²]^{1/2} (Eq. afinnorm). Here α is the angle between the spins and β, γ the angles of the spins to L̂, all at the initial separation (Eq. angles) — [Hofmann, Barausse & Rezzolla 2016, arXiv:1605.01938](https://arxiv.org/abs/1605.01938)
- HBR16 angle correction: Θ* = 2 arctan[(1 + ε_Θ) tan(Θ/2)] ≈ Θ + ε_Θ sin Θ, with ε_α ≈ 0 and ε_β = ε_γ ≈ 0.024 (fitted with the second form). It cuts the generic-spin residual bias from μ ≈ −0.005 (σ ≈ 0.007) to μ ≈ −0.001 (σ ≈ 0.006). The final-spin direction is cos θ_fin = Ĵ(r_in)·L̂(r_in) (Eq. dir), within 4° (6°) of NR in 64% (78%) of 157 cases and 18° at worst — [HBR16](https://arxiv.org/abs/1605.01938)
- HBR16 accuracy: rms σ ≈ 0.002 aligned and 0.006 generic over 619 NR simulations (abstract). On a blind test with 165 precessing GaTech runs, the n_M = 1, n_J = 2 formula has μ ≈ −0.004 (−0.0005 with adjusted angles) and σ ≈ 3.3×10⁻⁴ (3.5×10⁻⁴) as printed. Twin NR simulations differ by δa_fin ≈ 0.002 — [HBR16](https://arxiv.org/abs/1605.01938)
- BMR12 generic radiated energy: E_rad/M = [1 − Ẽ_ISCO(ã)]ν + 4ν²[4p₀ + 16p₁ã(ã+1) + Ẽ_ISCO(ã) − 1] (Eq. final_eq1). Here ã = L̂·(S₁+S₂)/M² = (|a₁|cos β + q²|a₂|cos γ)/(1+q)², with p₀ = 0.04827 ± 0.00039 and p₁ = 0.01707 ± 0.00032 fitted to equal-mass aligned runs only. Equal-mass misaligned NR data agree with residuals ≲ 1% and errors ≲ 10% of E_rad — [Barausse, Morozova & Rezzolla 2012, arXiv:1206.3803](https://arxiv.org/abs/1206.3803)
- NRSur7dq4Remnant was trained on 1528 simulations with q ≤ 4 and χ₁,₂ ≤ 0.8 in generic directions, using Gaussian-process regression. Its 95th-percentile out-of-sample errors are ≈ 5×10⁻⁴ M (mass), ≈ 2×10⁻³ (spin magnitude), ≈ 4×10⁻³ rad (spin direction) and ≈ 4×10⁻⁴ c (kick). It is at least an order of magnitude better than HBMR (HBR16 + BMR12 with n_M = n_J = 3), UIB, HL, HLZ and CLZM fits, even with LIGO's in-plane-spin post-processing — [Varma et al. 2019, arXiv:1905.09300](https://arxiv.org/abs/1905.09300)
- Jiménez-Forteza et al. 2017 (UIB) cover non-precessing binaries only, fitted to 427 NR runs. The final 3-D fit (lin+quad+mix) has RMSE 6.10×10⁻⁴ (9 coefficients, 298 points) for the final spin and 2.62×10⁻⁴ (8 coefficients) for E_rad (their summary tables) — [Jiménez-Forteza et al. 2017, arXiv:1611.00332](https://arxiv.org/abs/1611.00332)
- Kerr surface gravity: κ± = (r± − r∓)/(2(r±² + a²)) (Townsend ch. 4). An exercise gives κ = √(M⁴ − J²)/(2M(M² + √(M⁴ − J²))), which equals the user's c³√(1−χ²)/(2GM(1+√(1−χ²))) — [Townsend, Black Holes lecture notes, gr-qc/9707012](https://arxiv.org/abs/gr-qc/9707012). T = κ/(2π) in G = c = ħ = k = 1, and the first law reads δM = (κ/2π)δS + ΩδJ — [Wald, Living Rev. Rel. 4, gr-qc/9912119](https://arxiv.org/abs/gr-qc/9912119)
- **Recomputed** (`gap3_remnant_precessing.py`):
  - the generic formula at β, γ ∈ {0, π} equals |aligned Eq. 7| to 1.1e-16, and gives 0.68646 for q = 1 without spin;
  - with q = 1 and in-plane spins of 0.5 (or 0.9), |a_fin| = √((A/2)² + (ℓν)²) = 0.730566 (0.820809) when they are parallel and exactly ℓν = 0.686460 when they are antiparallel;
  - BMR12 at q = 1, |a| = 0.6, β = γ = 60° equals the aligned a = 0.3 value, 0.0600483; at aligned maximal spin it gives p₀ + 3p₁ = 0.09948;
  - κ(χ)/κ(0) = 2√(1−χ²)/(1+√(1−χ²)), which is exactly 8/9 at χ = 0.6, 3/4 at χ = 0.8 and 0 at χ = 1. Townsend's two forms agree, and dM/dA|_J = κ/(8π) = 0.0082883527 at χ = 0.7;
  - in SI, T_Schw(1 M☉) = 6.170074e-8 K (the library asserts 6.17e-8), so T(1 M☉, χ = 0.8) = 4.627555e-8 K.

  — [script](bh_gw_scripts/gap3_remnant_precessing.py)

### Inferences
- `remnant.rs`: add `final_spin_precessing(m1, m2, a1, a2, tilt1, tilt2, phi12) -> f64`. Compute cos α = cos β cos γ + sin β sin γ cos φ₁₂ from the in-plane azimuth difference φ₁₂. Order the masses so that q ≤ 1 (the existing `final_spin` is symmetric, but the generic formula is written with q = M₂/M₁ ≤ 1 and must swap the tilts with the masses). Reuse `kerr_isco` and the k_ij constants. About 25 lines.
  - Tests: reduces to `final_spin(...).abs()` for tilts 0/π; equal mass without spin gives 0.68646; the in-plane quadrature identities above; symmetric under swapping the bodies.
  - Optional: the ε = 0.024 map. It is a published fit refinement; leave it out unless asked.
- For BMR12, only the docs of `radiated_energy` need to say that a_i may be the projections |a_i|cos θ_i (BMR12's ã). No code.
- Risk: HBR16 itself says the angles are taken at the initial separation, and that "transitional precession" (spins nearly antialigned) breaks the direction formula. Expect σ ~ 0.006 in spin, much worse than NRSur7dq4Remnant's 2×10⁻³ at the 95th percentile.
- `gravitics::spacetime`: `pub fn kerr_hawking_temperature(m: f64, chi: f64) -> f64 { let s = (1.0 - chi * chi).sqrt(); horizon_temperature(C.powi(3) * s / (2.0 * G * m * (1.0 + s))) }`. Tests: equals `hawking_temperature` at χ = 0; ratio exactly 3/4 at χ = 0.8 and 8/9 at χ = 0.6; 0 at χ = 1. Note that `horizon_temperature` takes κ in s⁻¹ (κ_SI/c), consistent with the user's c³/(GM) form.

### Gaps
- Bardeen, Carter & Hawking 1973 and Hawking 1975, the primary sources for κ and T, are not on arXiv and I did not open them. The formula is cited from Townsend's notes and Wald's review.
- HBR16 gives no numeric generic-spin test cases, so the tests above are identities. NRSur7dq4Remnant evaluations (surfinBH) would give numbers, but that needs the package and its data, which I did not run.

## Gap 4: NLO Shapiro delay (lensing and retardation) for the Double Pulsar

### Takeaway
Kramer et al. replace Λ_u by Λ_u + δΛ_len + δΛ_ret, giving Δ_S = −2r ln(Λ_u + q_NLO(δΛ_len + δΛ_ret)). Here δΛ_len = 2rc/a_R = 2rs/(x + x_B), from Klioner & Zschocke's Eq. 73, and δΛ_ret is the velocity term inside the logarithm of Kopeikin & Schäfer's Eq. 130.

For the Double Pulsar, the lensing term lowers the delay by 0.757 µs at conjunction. Retardation is an antisymmetric ±1.1 µs within about ±11° of conjunction, and the related aberration-bending term is ±1.3 µs. Kramer fits q_NLO = 1.15(13).

### Cited Findings
- Kramer et al. 2021: Δ_S^{LO} = −2r ln Λ_u, with Λ_u = 1 − e_T cos u − s[sin ω(cos u − e_T) + (1 − e_T²)^{1/2} cos ω sin u]. Lensing enters as Λ_u → Λ_u + δΛ_u^len with δΛ_u^len = 2rc/a_R, an adapted version of Klioner & Zschocke Eq. (73) that they say has a few per cent error. Lensing reduces the computed propagation time, a consequence of Fermat's principle. With retardation the delay is Δ_S = −2r ln(Λ_u + δΛ^len + δΛ^ret), δΛ^ret being of order x/P_b and taken from Eq. (130) of Kopeikin & Schäfer. The velocity-dependent prefactor of the log is omitted (< 0.015 µs), and a_R = c(x + x_B)/s — [Kramer et al. 2021, arXiv:2112.06795](https://arxiv.org/abs/2112.06795)
- Kramer et al.: lensing alone is not yet observable because of its strong covariance with s; it shifts s by < 0.5σ. The effects are strongest near conjunction, where the eclipse (±1°) degrades timing. Lensing corrections to the retardation are below 0.02 µs post-fit — [Kramer et al. 2021](https://arxiv.org/abs/2112.06795)
- Kramer et al.: the aberration with signal bending is Δ_A = 𝒜[sin ψ + e_T sin ω] + 𝒟 cos ψ/Λ_u, with 𝒟 = (P/π) r/(x + x_B) and 𝒜 = P x/(P_b(1−e_T²)^{1/2} sin² i) ≃ 3.65 µs. q_NLO multiplies δΛ^ret, δΛ^len and 𝒟 together because the effects cannot be separated. The fit gives q_NLO = 1.15(13), "a ∼8σ measurement" of NLO propagation. Scaling only the deflection gives 1.26(24); scaling only the retardation gives 1.32(28) in the text equation but 1.32(24) in the summary table, a conflict inside the paper — [Kramer et al. 2021](https://arxiv.org/abs/2112.06795) (Table 3 and Sec. 6)
- Kramer et al. Table 3 inputs used: P_b = 0.1022515592973(10) d, x = 1.415028603(92) s, e_T = 0.087777023(61), ω₀ = 204.753686(47)°, z_s = 9.65(15), r = 6.162(21) µs, ν = 44.05406864196281 Hz, m_A = 1.338185, m_B = 1.248868 M☉, and R = m_A/m_B = x_B/x + O(c⁻⁴) — [Kramer et al. 2021](https://arxiv.org/abs/2112.06795)
- Klioner & Zschocke 2010, Eq. (73) ("elegant form"): cτ = R + (1+γ)m log[(x_obs + x_src + R + (1+γ)m)/(x_obs + x_src − R + (1+γ)m)] + O(m²/d) + O(m³), with m = GM/c². Their Eq. (72) is the explicit post-post-Newtonian form it repackages — [Klioner & Zschocke 2010, arXiv:1001.2133](https://arxiv.org/abs/1001.2133) (equation number checked in the arXiv PDF)
- Kopeikin & Schäfer 1999, Eq. (130): Δ_S(T) = −(2Gm_c/c³) ln{1 − e cos u − sin i[sin ω(cos u − e) + (1−e²)^{1/2} cos ω sin u] + (2π/sin i)(x/P_b)(m_p/m_c) e sin u − (2π sin i/(1−e²)^{1/2})(x/P_b)(m_p/m_c)[sin ω(cos u − e) + (1−e²)^{1/2} cos ω sin u] × [e cos ω + ((cos u − e) cos ω − (1−e²)^{1/2} sin ω sin u)/(1 − e cos u)]}. The velocity terms are of order 10⁻³–10⁻⁴. Their Eq. (119) is the expanded vector form −2m ln(r + k·r) + 2m(k·v_c) ln(r + k·r) + 2m[v_c·r − (k·r)(k·v_c)]/(r + k·r) — [Kopeikin & Schäfer 1999, gr-qc/9902030](https://arxiv.org/abs/gr-qc/9902030) (equation number checked in the arXiv PDF)
- **Recomputed** (`gap4_shapiro_nlo.py`):
  - 𝒜 = 3.6503 µs, reproducing Kramer's 3.65 µs;
  - x_B = 1.516229 lt-s, a_R = 878 826 km, δΛ_len = 4.2041×10⁻⁶ against Λ_u,min = 6.637×10⁻⁵ (1 − s = 6.44×10⁻⁵); peak LO delay 118.56 µs;
  - lensing changes the delay by −0.7569 µs at superior conjunction and −0.027 ns on the far side;
  - retardation amplitude 2πx_B/P_b = 1.0784×10⁻³; it changes the delay by −1.082/+1.189 µs (extremes within ±0.2 rad of eccentric anomaly around conjunction); with lensing, −1.493/+0.836 µs;
  - aberration NLO: 𝒟 = 15.189 ns, and 𝒟 cos ψ/Λ_u reaches ±1.30 µs.
  - Identities: at e = 0, δΛ_ret = −(2πx_B/P_b) s sin ψ cos ψ, with ψ = ω + u (residual 4e-19); at superior conjunction δΛ_ret ∝ (1/s − s), i.e. −1.1×10⁻⁸ here; Eq. 130 equals the Eq. 119 vector form to machine precision.
  - An independent construction, moving B to where it is when the signal passes (Δt = −k·r/c), agrees to first order in v. It leaves a constant ≈ v²/2 ≈ 5.5×10⁻⁷ offset, beyond KS99's order.

  — [script](bh_gw_scripts/gap4_shapiro_nlo.py)

### Inferences
- `energy.rs`: add `shapiro_delay_nlo(r, s, e, omega, u, x, x_b, pb, q_nlo) -> f64` returning −2r ln(Λ_u + q_nlo(δΛ_len + δΛ_ret)), with x m_p/m_c written as x_b (pulsar A as p, B as c). Keep `shapiro_delay` unchanged: it has one caller, the module test. About 20 lines.
  - Tests: q_nlo = 0 reproduces `shapiro_delay`; the e = 0 identity; δΛ_ret → 0 at conjunction as s → 1; the Double Pulsar lensing change −0.757 µs at conjunction to 1e-3 (a recomputed regression value).
  - Optionally add the 𝒟 cos ψ/Λ_u aberration term, as `aberration_bending`, since Kramer scale it with the same q_NLO.
- Pitfalls:
  - (i) KS99's text calls the focal parameter p = a_R(1−e²)^{1/2}; their Eq. 129 actually uses √(GM/p) with p = a_R(1−e²), the standard form. Implement Eq. 130 as printed, which is self-consistent.
  - (ii) KS99 uses k = −K = (0,0,−1), the propagation direction, and r = pulsar minus companion. The sign of the retardation is tied to that; my vector check confirms the printed signs with DD's Λ_u.
  - (iii) the lensing term is an approximation valid when the impact parameter is well above the Einstein radius (R_E ≈ 2500 km per Kramer), which holds outside eclipse;
  - (iv) q_NLO is not a PK parameter; Kramer stress that it only measures the relevance of the NLO terms.

### Gaps
- Kramer's text and table disagree on q_NLO[retardation] (1.32(28) or 1.32(24)). I did not check the published PRX version.
- I did not verify Kramer's "few %" accuracy claim for the KZ10 lensing term against an exact null-geodesic integration.

## Gap 5: Relativistic tidal Love number k₂

### Takeaway
Integrate the TOV equations together with the first-order Riccati equation for y = rH'/H (Postnikov et al.'s form of Hinderer's equation). Start from y(0) = 2, apply the surface jump y → y − 4πR³ε₋/M wherever the density is discontinuous, and evaluate k₂ with Hinderer's corrected Eq. 23. Below C ≈ 0.02, use a verified small-C series, because Eq. 23 cancels catastrophically there.

Embed Read et al.'s piecewise-polytrope parameters (4 numbers per EOS plus a 4-piece SLy crust). My Python version reproduces:
- DN09's closed form k₂(4/9) = 0.0017103 to 16 digits;
- DN09's γ = 2 μ-polytrope fit to ≤ 3×10⁻⁴;
- Read's SLy and AP4 R_1.4 and M_max to ≤ 0.05%;
- Hinderer et al. 2010's SLy k₂(1.4 M☉) = 0.0763 to 0.1%.

### Cited Findings
- Hinderer 2008 (arXiv version, which carries the erratum on its first pages) gives the perturbation equation H'' + H'[2/r + e^λ(2m/r² + 4πr(p − ρ))] + H[−6e^λ/r² + 4πe^λ(5ρ + 9p + (ρ+p)/(dp/dρ)) − ν'²] = 0. Near the centre H ≈ a₀r²[1 − (2π/7)(5ρ_c + 9p_c + (ρ_c+p_c)/(dp/dρ)_c)r²]. The corrected Eq. (23) is k₂ = (8C⁵/5)(1−2C)²[2 + 2C(y−1) − y] × {2C(6 − 3y + 3C(5y−8)) + 4C³[13 − 11y + C(3y−2) + 2C²(1+y)] + 3(1−2C)²[2 − y + 2C(y−1)] ln(1−2C)}⁻¹. The erratum also corrects Eq. (20), the exterior solution, and Table 1 — [Hinderer 2008 + erratum, arXiv:0711.2420](https://arxiv.org/abs/0711.2420)
- Erratum Table 1 (k₂ against M/R): n = 1 gives 0.2599 (1e-5), 0.122 (0.1), 0.0776 (0.15), 0.0459 (0.2), and a second "0.2" row reading 0.0253, evidently 0.25. n = 0.5 gives 0.4491, 0.251, 0.173, 0.095, 0.0569; n = 0.3 gives 0.5511, 0.294, 0.221, 0.119; n = 0.7 and 1.2 rows are also listed. Hinderer writes the polytrope as P = Kρ^{1+1/n}, with ρ the density in T_αβ, i.e. the energy density — [Hinderer 2008 erratum](https://arxiv.org/abs/0711.2420)
- Damour & Nagar 2009 say the published Hinderer paper had errors that invalidate its conclusion of a mild C-dependence, and that the erratum did not correct her Eq. (27) or Fig. 2 — [Damour & Nagar 2009, arXiv:0906.0096](https://arxiv.org/abs/0906.0096)
- DN09 define the e-polytrope (p = κe^γ) and the μ-polytrope (p = κμ^γ, e = μ + p/(γ−1)), and for γ = 2 fit k_ℓ = k_ℓ^N(1−2c)² Σ_{n=0}^4 a_n c^n. For ℓ = 2 their Table I gives a₀..a₄ = 0.9991, −2.9287, −1.1373, 14.0013, −50.9711; for ℓ = 3 and 4 too. The maximum compactness is 0.2145 (μ) and 0.265 (e) for γ = 2, and k₂^N(γ = 2) = 0.25991 — [DN09](https://arxiv.org/abs/0906.0096)
- DN09 incompressible star. The (e+p)de/dp term carries a δ-function at the surface, so y_out = y_in − 3, the Newtonian k_ℓ = 3/(4(ℓ−1)) (k₂ = 3/4), and at c_max = 4/9 the interior solution is analytic with y_out = 3ℓ − 4. This gives k₂(4/9) = 4096/(10935(308 − 81 ln 3)) ≃ 0.0017103. Their realistic EOS use linear interpolation of tables, to avoid spurious sound-speed oscillations — [DN09](https://arxiv.org/abs/0906.0096)
- Postnikov, Prakash & Lattimer 2010 rewrite the problem as the Riccati equation r y' + y² + y e^λ[1 + 4πr²(p − ρ)] + r²Q = 0, with Q = 4πe^λ(5ρ + 9p + (ρ+p)/c_s²) − 6e^λ/r² − ν'², y(0) = 2, and e^λ = (1−2m/r)⁻¹, ν' = 2e^λ(m + 4πpr³)/r². They warn of severe cancellations in Eq. 23 at small β and give a Taylor series for β < 0.1. A density jump gives y(r_d+) = y(r_d−) − 3Δρ/ρ̃, with ρ̃ = m(r_d)/(4πr_d³/3). The incompressible case gives y_R = −1 and k₂ = 3/4; n = 1 Newtonian gives y_R = (π²−9)/3 and k₂ = (15−π²)/(2π²). Polytropes with n < 1 need care because dρ/dh diverges at the surface. They report agreement with Hinderer, DN09 and Binnington & Poisson "to within numerical accuracy" — [Postnikov et al. 2010, arXiv:1004.5098](https://arxiv.org/abs/1004.5098)
- Hinderer et al. 2010 integrate H and β = dH/dr with the TOV equations from H = a₀r², β = 2a₀r. For a finite surface density, y = Rβ(R)/H(R) − 4πR³ε₋/m. Their Table I (1.4 M☉, tabulated EOS) gives SLy R = 11.74 km, m/R = 0.176, k₂ = 0.0763, λ = 1.70×10³⁶ g cm² s²; AP3 12.09, 0.171, 0.0858; AP1 9.36, 0.221, 0.0512; MPA1 12.47, 0.166, 0.0924; and so on — [Hinderer, Lackey, Lang & Read 2010, arXiv:0911.3535](https://arxiv.org/abs/0911.3535)
- Read et al. 2009 piecewise polytropes: p = K_iρ^{Γ_i} in pieces, ε = (1 + a_i)ρ + K_iρ^{Γ_i}/(Γ_i − 1), with a_i set by continuity. In c = 1 form, cgs needs p → p/c² and K → K/c², with ε and ρ in g/cm³. The core dividing densities are fixed at ρ₁ = 10^{14.7} and ρ₂ = 10^{15.0} g/cm³, with p₁ = p(ρ₁). The first core piece is extended down until it meets the SLy crust.
  - Crust (Table II, K_i, Γ_i, ρ_i): (6.80110e-9, 1.58425, 2.44034e7), (1.06186e-6, 1.28733, 3.78358e11), (5.32697e1, 0.62223, 2.62780e12), (3.99874e-8, 1.35692, up to the match density).
  - Core (Table III, log p₁ [dyne/cm²], Γ₁, Γ₂, Γ₃): SLy 34.384, 3.005, 2.988, 2.851; AP4 34.269, 2.830, 3.445, 3.348; and so on for 34 EOS.
  - Table III also lists observables from the tabulated EOS and the fit errors: SLy M_max 2.049 (+0.02%), R_1.4 11.736 km (−0.21%); AP4 M_max 2.213 (−1.08%), R_1.4 11.428 km (−0.90%).

  — [Read, Lackey, Owen & Friedman 2009, arXiv:0812.2163](https://arxiv.org/abs/0812.2163)
- **Recomputed** (`gap5_love_tov.py`, `gap5_checks2.py`, `gap5_checks3.py`):
  - k₂(C = 4/9, y = 2) = 0.0017103005961490 in floats and in 60-digit Decimal, equal to DN09's closed form.
  - Incompressible TOV + jump gives k₂ = 0.556257, 0.400305, 0.183019, 0.063206, 0.010347 and 0.001758 at C = 0.05, 0.1, 0.2, 0.3, 0.4, 0.444.
  - γ = 2 μ-polytrope: k₂ = 0.17919, 0.11710, 0.07075, 0.03729 at C = 0.05, 0.10, 0.15, 0.20, against the DN09 fit's 0.17921, 0.11706, 0.07077, 0.03727. The γ = 2 e-polytrope gives 0.12321, 0.07837, 0.04589 at C = 0.1, 0.15, 0.2, matching Hinderer's n = 1 rows (0.122, 0.0776, 0.0459) to ≤ 1%. So Hinderer's ρ is the energy density.
  - Hinderer's n = 0.5 rows do not match either convention: e-polytrope 0.2284/0.1531/0.0968 against 0.251/0.173/0.095, stable to 1e-4 under the surface stopping criterion. Her n = 0.3 rows are off by up to 6%.
  - Newtonian Lane–Emden + Riccati gives k₂ = 0.25991 (n = 1), 0.44915 (0.5), 0.36264 (0.7), 0.20619 (1.2), matching her M/R = 1e-5 rows. n = 0.3 gives 0.55315 against her 0.5511; my n < 1 surface handling is suspect there.
  - Read SLy piecewise polytrope: M_max = 2.0484 against Read's fit value 2.049 × 1.0002 = 2.0494; R_1.4 = 11.705 km against 11.736 × (1 − 0.0021) = 11.711 km; C = 0.1766, k₂ = 0.07620, Λ_1.4 = 295.6. Hinderer 2010's tabulated SLy gives k₂ = 0.0763, and Λ = (2/3)k₂/C⁵ = 300.4 with C from R = 11.74 km.
  - Read AP4: M_max 2.1879 against 2.1891; R_1.4 11.319 km against 11.325; k₂ = 0.07548, Λ_1.4 = 247.6.
  - The rounded crust constants leave a 1.2×10⁻⁴ relative pressure jump at the piece boundaries.

  — [script](bh_gw_scripts/gap5_love_tov.py), [checks](bh_gw_scripts/gap5_checks2.py), [checks](bh_gw_scripts/gap5_checks3.py)
- **Recomputed: Eq. 23 cancellation and the Postnikov series.** In double precision Eq. 23 has relative error ~1 at C ≤ 1e-4, 1e-3 at C = 1e-3, 1e-7 at C = 0.01 and 2e-11 at C = 0.05.
  - The true expansion is k₂ = ((1−2C)²/2) Σ b_k C^k, with b₀ = (2−y)/(y+3), b₁ = (y² + **6**y − 6)/(y+3)², b₂ = (y³ + 34y² − 8y + 12)/(7(y+3)³) and b₃ = (y⁴ + 62y³ + 84y² + 48y + 36)/(7(y+3)⁴). These were verified against 60-digit evaluations at y = −1, 0, 0.5, 1 and (π²−9)/3; the residual scales as C⁴.
  - Postnikov's printed series has −6y in the β term (it gives −0.707 against the true −0.386 at y = 0.29). It labels the b₃ term as β⁴, and its "β⁵" term is about half the true b₄ (0.0560 against 0.1124).

  — [checks](bh_gw_scripts/gap5_checks2.py)

### Inferences
- `tides.rs` design (about 150–200 lines):
  - (a) `pub fn love_number(compactness, y) -> f64`: Eq. 23 for C ≥ 0.02, and the 4-term series b₀..b₃ below that (error ≲ 3×10⁻⁸ at the switch on both sides);
  - (b) a fixed-step RK4 (or adaptive RK45) integration of (m, p, y) in r, or better in enthalpy h (Read et al. and Postnikov et al. use h, which handles the surface cleanly), from y = 2 at the centre, with the jump y −= 4πR³ε₋/M for a finite surface density and the same jump at internal discontinuities using the enclosed mean density;
  - (c) a small `PiecewisePolytrope` holding (log p₁, Γ₁, Γ₂, Γ₃) and the SLy crust, with SLy and AP4 constants.

  Then Λ = (2/3)k₂/C⁵ via the existing `lambda`.
- Tests:
  - the DN09 closed form k₂(4/9, 2) (exact, 1e-15);
  - incompressible TOV + jump against 3/4 in the C → 0 series limit, and smooth approach to 0.00171 at 4/9;
  - k₂ → (15−π²)/(2π²) as C → 0 for an n = 1 polytrope (via the series), complementing the existing Clairaut–Radau test;
  - the DN09 γ = 2 μ-polytrope fit at C = 0.05–0.2 within 1e-3 relative;
  - SLy and AP4 piecewise polytropes: R_1.4 within 0.2% of Read's fit-adjusted 11.711 and 11.325 km, M_max within 0.2% of 2.0494 and 2.1891;
  - SLy k₂(1.4) within 1% of Hinderer et al.'s 0.0763 (their value is for the tabulated EOS).
- Pitfalls:
  - (i) Read's K_i are in the c = 1 form: p/c² = Kρ^Γ, with ρ and p/c² in g/cm³, although the Table II caption says p is in dyne/cm². Using K with p in dyne/cm² gives nonsense crust pressures.
  - (ii) e- versus μ-polytropes differ by 5–20% in k₂ at C = 0.1–0.2. DN09's fits are μ-polytropes, Hinderer's table is e-polytropes.
  - (iii) Postnikov's jump equation has the density difference with the wrong sign in its middle expression; use y_out = y_in − 3Δρ/ρ̃ with Δρ = inside minus outside, which is DN09's −3 for the incompressible surface.
  - (iv) Do not use Postnikov's printed small-β series or Hinderer's n ≠ 1 table rows as test values.
  - (v) n < 1 polytropes have a divergent (ρ+p)/c_s² at the surface. It is integrable, but integration in r near the surface converges slowly, so integrate in enthalpy or ln h.
  - (vi) For tabulated EOS, DN09 interpolate linearly to avoid spurious sound-speed oscillations, but the dp/dρ entering Q must be consistent with the interpolation.
- EOS data: the piecewise-polytrope constants are a handful of published numbers (Read et al. Tables II–III), small enough to embed as constants with a citation. That avoids redistributing tabulated EOS files.

### Gaps
- I did not open any source tabulating Λ_1.4 for SLy or APR4. The 295.6 (SLy) and 247.6 (AP4) above are recomputed from Read's piecewise polytropes only.
- I did not check the licence terms of tabulated EOS repositories (CompOSE, LORENE, Lattimer's tables), or of the Read et al. tables beyond their being published numbers. Embedding about 20 published constants is very likely fine, but this is an inference, not a checked licence.
- Binnington & Poisson 2009, cited by Postnikov, was not opened.
- Hinderer's n = 0.5 and 0.3 finite-C rows disagree with my solver by up to 13% (n = 0.5) and 6% (n = 0.3). I could not settle whether that is her table or my n < 1 surface treatment: the n = 0.5 values are stable under my stopping criterion, but I did not cross-check them with an enthalpy-variable integration.

## Cross-cutting: what to implement, what to test, what to watch

### Takeaway
All five gaps fit the "few hundred lines of Rust" budget, and only the QNM solver needs real care.

| Gap | Change | Size |
|---|---|---|
| 1 | Kerr f_ISCO | 3 lines |
| 2 | Leaver + spectral QNM solver | ~200 lines incl. a complex type |
| 3 | HBR16 generic final spin; Kerr temperature | ~25 + 3 lines |
| 4 | NLO Shapiro | ~20 lines |
| 5 | TOV + Love solver with piecewise polytropes | ~150–200 lines |

Each has exact identities plus at least one published high-precision value, with two exceptions. Gap 2's only 8+ digit published Kerr value is the `qnm` README a = 0.68. Gap 5 has no published Λ_1.4 values that I opened.

### Cited Findings
- Exact identities available as tests, all verified in the scripts:
  - Gap 1: MΩ_ISCO = 6^{-3/2}, 1/2 and 1/26 at a = 0, +1, −1, and Isoyama's z₀''(Ω_ISCO) = 0 — [Isoyama et al. 2014](https://arxiv.org/abs/1404.6133);
  - Gap 2: the Schwarzschild overtone digits in Berti's s2l2.dat — [Berti data (archive)](https://web.archive.org/web/20240318065031/https://pages.jh.edu/~eberti2/ringdown/s2l2.dat);
  - Gap 3: κ ratios 3/4 at χ = 0.8 and 8/9 at χ = 0.6 — [Townsend gr-qc/9707012](https://arxiv.org/abs/gr-qc/9707012);
  - Gap 4: δΛ_ret at e = 0 equals −(2πx_B/P_b) s sin ψ cos ψ — [Kopeikin & Schäfer 1999](https://arxiv.org/abs/gr-qc/9902030);
  - Gap 5: k₂(4/9) = 4096/(10935(308 − 81 ln 3)) — [Damour & Nagar 2009](https://arxiv.org/abs/0906.0096).
- Errata and misprints that affect tests:
  - BCW06 Table II (2,−2,0) j = 0.6 and Table IV (4,2,0) j = 0 last digits, and the likely dropped leading digit in Table X's Q-% column — [BCW06](https://arxiv.org/abs/gr-qc/0512160);
  - Hinderer 2008 Eqs. 20 and 23 and Table 1 (corrected in the erratum; DN09 note that Eq. 27 and Fig. 2 remained wrong) — [arXiv:0711.2420](https://arxiv.org/abs/0711.2420), [DN09](https://arxiv.org/abs/0906.0096);
  - Postnikov's small-β series and the sign in the middle of the jump equation — [arXiv:1004.5098](https://arxiv.org/abs/1004.5098);
  - KS99's focal-parameter typo — [gr-qc/9902030](https://arxiv.org/abs/gr-qc/9902030);
  - Kramer's q_NLO[retardation] error bar, 28 against 24 — [arXiv:2112.06795](https://arxiv.org/abs/2112.06795).

### Inferences
- Order of work by value per line: Gap 1 and the Kerr temperature (trivial); Gap 4 (small, Double Pulsar regression values in hand); Gap 3 (small, identities only); Gap 5 (medium; the TOV integrator is reusable); Gap 2 last (largest, but it replaces fits whose quality-factor errors reach about 3%).
- The repository's issue workflow (branch per issue, PR closes it) suggests one issue per gap. The BCW06 misprints and the Hinderer table discrepancy are novel enough to comment on the relevant issue.

### Gaps
- Primary sources not opened, all off arXiv: BPT 1972, Leaver 1985, Bardeen–Carter–Hawking 1973, Hawking 1975.
- Data not obtained: Cook's Zenodo HDF5 (403), Berti's SharePoint Kerr tables, and the NRSur7dq4Remnant data and surfinBH evaluations.

## Implemented in the crates (2026-10-27, #4)

Items 43–49 are in `gravitational`, and 47 is in `gravitics`. Building them corrected these research values:
- **Love numbers (item 49):** this script's solver stopped at p = 10⁻¹⁴p_c and lost the outer 0.95 m of the star. Converged, SLy gives R₁.₄ = 11.706 km and k₂ = 0.07617, and AP4 gives 11.3197 km and k₂ = 0.07545. Λ does not change: 295.6 and 247.6, or 295.63 and 247.61 with metrology's solar mass. The DN09 agreement is 3×10⁻⁴ absolute but up to 5.5×10⁻⁴ relative. The series switch at C = 0.02 is 3×10⁻⁸ absolute.
- **Hinderer's n = 0.5 and 0.3 rows:** these are not caused by n < 1 surface handling. An independent enthalpy-variable solver agrees with this script to 10⁻⁴, and Hinderer's rows still differ by up to 13%, with errors that change sign, so the fault is in her table.
- **Shapiro retardation (item 48):** the factor-49 cross-check moved the companion to its passage position. Kopeikin & Schäfer define r with both bodies at emission time. An exact delay past a uniformly moving mass confirms their Eq. 130 to first order, with a misfit that falls as λ².
- **Kerr QNMs (item 44):** the depth must scale as 1/√(1 − j²) from j = 0, not only above 0.99. Five BCW06 entries are rounded up from just above the half-unit, so 6×10⁻⁵ is the right tolerance. The solver covers the fits' 63 modes (n ≤ 2).
