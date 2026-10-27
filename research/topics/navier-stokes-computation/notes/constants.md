# Every explicit constant in "Velocity-Field Detection of Computation in Forced Navier–Stokes Flows", recomputed

Scope: the explicit constants, bounds and inequalities in OpenAI's preprint (27 September 2026, Apache-2.0, not peer reviewed), issue #53. Each one is recomputed independently, and the note records whether it holds, how tight it is, and anything wrong. The theorems themselves (Theorem 1.1, Corollary 6.1) are **mainstream mathematics (preprint)**. Nothing here re-proves them; it checks the arithmetic they rest on.

Source: the manuscript PDF and its extracted text, read from the session scratchpad (`dl_ns/`, untrusted data, not committed). Equation numbers are the paper's.

Labels: **(paper)** marks a value printed in the manuscript; **(computed)** marks a value produced by `constants_checks.py`.

Checks: `notes/constants_checks.py`, stdlib only, run as `python3 -I constants_checks.py`. It makes 86 top-level assertions (several loop over hundreds of cases), runs in about 4 s, and exits 1 on any failure. All pass.

**Verdict in one line:** every inequality the paper asserts holds. None fails. Most are loose by factors of 10² to 10⁷, which is harmless for the proofs. Two are tight: eq. 18's 1/50, and §4's nonhalting clearance of exactly 2d. One constant the paper leaves implicit, the bump's derivative bound C, comes out at 41. The paper-scale parameters are far too extreme to simulate, so the planned numerical example must use relaxed parameters.

| Item | Paper | Verdict | Computed truth | Slack |
|---|---|---|---|---|
| Shell sum Σ_k e^{−\|Z+k\|²/8}, Lemma 3.3 | ≤ 1 + 8·1056 = 8449 < 10⁴ | holds | sup = 8π = 25.13 (at Z = 0) | 336× (vs 8449) |
| e^{−1/32} ≤ 32/33; Σ j(32/33)^j | 1056 | holds, exact | Σ j e^{−j/32} = 1023.92 | 3% |
| Heat bound, d = 1/32, eq. 14 | 8·10⁴/(4πd²) = 6.52·10⁶ < H = 10⁹ | holds | sup K = 1/(πe d²) = 119.91 | H/sup = 8.3·10⁶ |
| Heat bound, d = 1/16, eq. 14 | 1.63·10⁶ < H = 10⁸ | holds | sup K = 29.98 | H/sup = 3.3·10⁶ |
| Heat bound, τ ≥ 1 | ≤ 10⁴/(4π) = 796 | holds | K(1, 0) = 1 + 4e^{−4π²} ≈ 1 | 796× |
| M₀ = Σ 4b_n², eq. 13 | 4·10⁻¹²/(4R² − 1) < 10⁻¹¹ | holds, exact | ≤ 6.35·10⁻¹⁴ (R ≥ 4) | 157× |
| ∫g_n, eq. 11 | ≤ 4b_n² | holds | = 2.25 b_n² | 1.8× |
| 2δ_nD_n, eq. 18 | < 1/50 < 1/16 | holds | → 1/50 as D_n → ∞ | **tight** at 1/50 |
| Support radius √2 b_nRⁿ | = √2·10⁻⁶·2⁻ⁿ < 1/32 | holds, exact | ≤ 7.07·10⁻⁷ | 44 000× |
| 2HM₀ + 1/16, eq. 20 | < ½ | holds | ≤ 0.06263 (d = 1/32) | 8× |
| Nonhalting clearance, §4 | ≥ 2d = 1/16 | holds | exactly 1/16 | **tight** (≥ allows it) |
| §4 second pair | 7/24 > 1/4 | holds, exact | only ≥ 1/8 is needed | — |
| \|P″\|, eq. 37 | ≤ 360 | holds | 10/√3 = 5.774 | 62× |
| \|A∗″\| | ≤ 1440 | holds | 40/√3 = 23.09 | 62× |
| ‖∆φ_R‖∞ | ≤ 2880R⁻² ≤ 3000R⁻² | holds | 43.49R⁻² | 66× |
| T_n/R_n², eq. 39 | (16Λ/R₀)(K₀D(D/Λ)ⁿ + 2Λ⁻ⁿ) | exact identity | — | — |
| δ, eq. 39 | < 1/48 + 1/48 = 1/24 | holds (160 parameter sets) | ≈ (ν/(1+ν))/48 | detection needs only δ < ½ |
| R_{n+1} ≥ 2R_n | Λ ≥ 4 | holds | Λ = 4b² ≥ 16 | 8× |
| y′ ≤ bB_n − 1, §5.1 | left row | holds | attained | **tight** (as stated) |
| b(r) + b(1 − r) | ≥ e⁻² | holds | min = 2e⁻² at r = ½ | 2× |
| Bump derivative bound C, eq. 11 | "an effective integer C ≥ 1" | not given | ‖∇g‖ = 4, ‖D²g‖ = 40.46, so C = 41 | — |

---

## 1. Smooth step, bump and the bump constant C (§3.1, §3.4, §5.5)

### Takeaway
All the claimed properties hold. The paper never states C. The smallest admissible integer is 41, set by the Hessian.

### Cited Findings
- (paper) b(r) = e^{−1/r} for r > 0, else 0; θ(r) = b(r)/(b(r) + b(1 − r)); θ is 0 for r ≤ 0, 1 for r ≥ 1, flat at both ends.
- (paper) g(Z) = Π_j θ(2(Z_j + 1))θ(2(1 − Z_j)) "lies in [0, 1], equals one near zero, and is supported in [−1, 1]²".
- (paper) "Let an effective integer C ≥ 1 bound its first and second derivative norms."
- (paper) b(r) + b(1 − r) ≥ e⁻², and s^m e^{−s} ≤ (m + j)! s^{−j}.

### Inferences
- (computed) θ(r) + θ(1 − r) = 1, θ is nondecreasing, ∫₀¹θ = ½, and max θ′ = θ′(½) = 2.
- (computed) min over ℝ of b(r) + b(1 − r) is 2e⁻² = 0.2707, at r = ½. The paper's e⁻² holds with a factor of 2 to spare.
- (computed) g = 1 exactly on [−½, ½]², which is what "near zero" means here. g = 0 outside [−1, 1]², and ∫g = 1.5² = 2.25.
- (computed) max ‖∇g‖ = 4 and max ‖D²g‖_op = 40.46 (max |θ″| = 9.84, so |g₁″| ≤ 4·9.84 = 39.36), so **C = 41**. This matters for implementation, because D_n is proportional to C.
- (computed) The tail estimate s^m e^{−s} ≤ (m + j)! s^{−j} holds on the sampled grid. It is the Taylor bound e^s ≥ s^{m+j}/(m + j)!.

### Gaps
- C = 41 comes from a fine grid search (refined to step 5·10⁻⁵), not an interval-arithmetic certificate. Its margin to 41 is 0.54.

## 2. Theorem 3.1 scales: eqs. 10, 11, 13, 16, 18 and the support radius

### Takeaway
All hold, most by wide margins. Eq. 18's 1/50 is the only tight one. The mass bound 10⁻¹¹ is 157 times larger than needed.

### Cited Findings
- (paper) b_n = 10⁻⁶(2R)⁻ⁿ, D_n = 2Cb_n⁻²(1 + nL)R^{3n}, δ_n = 1/(100(1 + D_n)), t_n = 2n, R = 4^L with L ≥ 1 (eqs. 8, 10).
- (paper) ∫g_n ≤ 4b_n² (eq. 11); M₀ = Σ4b_n² = 4·10⁻¹²/(4R² − 1) < 10⁻¹¹ (eq. 13).
- (paper) ‖DΨ_s^{±1}‖ ≤ e^{Ln} ≤ Rⁿ, ‖D²Ψ_s^{±1}‖ ≤ nLe^{3Ln} ≤ nLR^{3n} (eq. 16).
- (paper) ‖w_new − w⁰‖∞ ≤ 2δ_nD_n < 1/50 < 1/16 (eq. 18); √2 b_nRⁿ = √2·10⁻⁶·2⁻ⁿ < 1/32 ≤ d.

### Inferences
- (computed) The geometric sum is exact. With x = 1/(4R²), 4·10⁻¹²·x/(1 − x) = 4·10⁻¹²/(4R² − 1), checked in exact rationals for L = 1…5. For R ≥ 4, M₀ ≤ 4·10⁻¹²/63 = 6.35·10⁻¹⁴. The claim < 10⁻¹¹ holds even for R = 1, so L ≥ 1 is not needed for it.
- (computed) ∫g_n = 2.25 b_n², so eq. 11's 4b_n² is loose by 1.78×.
- (computed) 2δ_nD_n = 2D_n/(100(1 + D_n)): 1/100 at D = 1, 1/55 at D = 10, and 1/50 to 15 digits at the real D_n. The strict < 1/50 always holds, but this bound has no slack. The 1/16 used later has 3× slack.
- (computed) The Hessian trace 2C(b_n⁻²R^{2n} + b_n⁻¹nLR^{3n}) ≤ D_n for L = 1, 2, 5 and n = 1…7. The wait 2 − 2δ_n > 1 always holds.
- (computed) In eq. 16 the variation-of-constants integral is exactly e^{2Ls} − e^{Ls}. The printed nLe^{3Ln} is valid but loose: e^{2Ln} ≤ R^{2n} would already do, and D_n could be cut to about 2Cb_n⁻²(1 + nL)R^{2n}. The proof does not need this.
- (computed) The support radius is at most √2·10⁻⁶/2 = 7.07·10⁻⁷, against 1/32, a margin of 44 000×.
- (computed) The real scales are extreme. At L = 1 and C = 41: b₁ = 1.25·10⁻⁷, D₁ = 6.72·10¹⁷, δ₁ = 1.49·10⁻²⁰. Each later block shrinks b_n by 8 and δ_n by about 4096.

### Gaps
- The bounds ‖DΨ‖ ≤ e^{Ln} were taken as the standard Gronwall argument and not re-derived for a specific V.

## 3. Lemma 3.3: the shell bound and the heat-kernel bounds (eq. 14)

### Takeaway
Both (d, H) pairs hold by six to seven orders of magnitude. The true supremum of the periodised unit-torus kernel over {dist ≥ d, τ ≤ 1} is 1/(πe d²), about 120 at d = 1/32 and 30 at d = 1/16. Over τ ≥ 1 it is 1 + 4e^{−4π²}, essentially 1.

### Cited Findings
- (paper) K(τ, Z) = (4πτ)⁻¹ Σ_k e^{−|Z+k|²/(4τ)}. For Z ∈ [−½, ½]², the shell |k|∞ = j "has 8j points and |Z + k| ≥ j/2".
- (paper) Σ_k e^{−|Z+k|²/8} ≤ 1 + 8Σ_j j e^{−j/32} < 10⁴, using e^{−1/32} ≤ 32/33 and Σ j q^j = q/(1 − q)², "an upper bound 1 + 8·1056 < 10⁴".
- (paper) For τ ≤ 1, "τ⁻¹e^{−d²/(8τ)} ≤ 8/d²", so K ≤ 8·10⁴/(4πd²), "less than 10⁹ for d = 1/32 and less than 10⁸ for d = 1/16". The maximum principle covers τ ≥ 1.

### Inferences
- (computed) e^{1/32} ≥ 1 + 1/32 = 33/32, so e^{−1/32} ≤ 32/33. In exact rationals, (32/33)/(1/33)² = 1056, and 1 + 8·1056 = 8449 < 10⁴. All exact.
- (computed) The shell has 8j points (checked to j = 29). In fact |Z + k| ≥ j − ½ ≥ j/2, so the paper uses the weaker half.
- (computed) The true sum Σ_j j e^{−j/32} is 1023.92, so 1056 is only 3% loose. The looseness is upstream: e^{−j²/32} was replaced by e^{−j/32}, and j − ½ by j/2. The actual supremum of the shell sum is (Σ_n e^{−n²/8})² = 8π(1 + O(e^{−8π²})) = 25.13, at Z = 0, against the bound 8449.
- (computed) sup_τ τ⁻¹e^{−d²/(8τ)} = 8/(e d²). The paper's 8/d² drops the 1/e.
- (computed) Paper bounds: 6.519·10⁶ at d = 1/32 and 1.630·10⁶ at d = 1/16.
- (computed) True suprema, found by golden-section search over τ ∈ [10⁻⁷, 1] on the circle |Z| = d (41 angles), with an off-circle grid confirming the circle is extremal:
  - d = 1/32: 119.9101 = 1/(πe d²), at τ = d²/4. H/sup = 8.34·10⁶, and the paper's intermediate bound is 5.4·10⁴ times the truth.
  - d = 1/16: 29.9775. H/sup = 3.34·10⁶.
- (computed) τ ≥ 1: the image sum equals the Fourier series Π_j Σ_m e^{−4π²m²τ}cos(2πmZ_j) to 10⁻¹², so sup K(1, ·) = K(1, 0) = 1 + 4e^{−4π²} + … ≈ 1 + 2.9·10⁻¹⁷. The paper's shell estimate gives 796. The kernel has unit mass (checked to 10⁻¹²).
- (computed) Consequently eq. 15, w(t_n) ≤ HM₀, also holds with H replaced by about 1, because the wait before every block exceeds 1.

### Gaps
- None for the arithmetic. The heat-kernel uniqueness argument ("uniqueness identifies it as the heat kernel") is standard and was not re-proved.

## 4. Eq. 20 and the §4 detector geometry

### Takeaway
Eq. 20 holds 8× below the threshold. The torus detector's nonhalting clearance is exactly 2d, which the hypothesis "≥ 2d" allows. That is the one place where the geometry has no slack.

### Cited Findings
- (paper) w ≤ 2HM₀ + 1/16 < ½ on E, "for both pairs in (9)".
- (paper) E = {1/32 < y < 1/8}, and the nonhalting trajectory satisfies 3/16 ≤ y ≤ 7/8; "The circle distance from this corridor to E is at least 1/16". Then d = 1/32 is used.
- (paper) The second pair: "min{2/3 − 3/8, 1 + 1/8 − 5/6} = 7/24 > 1/4".

### Inferences
- (computed) Using only M₀ < 10⁻¹¹: 2·10⁹·10⁻¹¹ + 1/16 = 0.0825 < ½, and 0.0645 for the second pair. Using M₀ ≤ 4·10⁻¹²/63: 0.06263 and 0.06251. The 1/16 diffusion allowance dominates.
- (computed) The corridor-to-E distance is min{3/16 − 1/8, 1 + 1/32 − 7/8} = min{1/16, 5/32} = 1/16 = 2d exactly. That is an infimum against the open set E. The theorem's hypothesis is "≥ 2d", so it holds, but a planar processor whose corridor dips below 3/16 by any amount would break it.
- (computed) 7/24 is exact on both sides. The relevant requirement is ≥ 2d = 1/8 for d = 1/16. The paper's "> 1/4" is true but not the condition used.

### Gaps
- The corridor bounds 3/16 ≤ y ≤ 7/8 and p = (¼, ¼) come from the companion preprint [12, Theorem 4.1], which was not opened.

## 5. Injection height, and the planar derivative kernel

### Takeaway
The time-derivative factor integrates to exactly 1, which gives the height-one injection. The derivative-kernel constant of Lemma 5.1 has a closed form.

### Cited Findings
- (paper) h = Σ δ_n⁻¹θ′((t − t_n)/δ_n)g_n (eq. 12), and w⁰ = θ((t − t_n)/δ_n)g_n during injection (eq. 17).
- (paper) "‖∂_ℓH_s‖₁ ≤ C_ν s^{−1/2}" with H_s = (4πνs)⁻¹e^{−|X|²/(4νs)}; ∫g dX = θ′(t) (eq. 32).

### Inferences
- (computed) ∫ over [t_n, t_n + δ_n] of δ_n⁻¹θ′((t − t_n)/δ_n) dt = θ(1) − θ(0) = 1. Simpson's rule gives 1 to 10⁻¹⁰ at δ = 10⁻² and 10⁻⁶. Since a = 0 during injection, w⁰_t = h and w⁰ reaches g_n, of height 1.
- (computed) ‖∂_ℓH_s‖₁ = (πνs)^{−1/2} exactly (∂₁H = −X₁H/(2νs), with E|X₁| = 2√(νs/π)). Checked at ν = 0.01, 1 and 7. So C_ν = (πν)^{−1/2} per direction.
- (computed) ∫θ′ = 1, so eq. 32's spatial factor integrates to 1, and M(t) = θ(t) (eq. 34).

### Gaps
- None.

## 6. Eq. 37: the quintic cutoff φ_R

### Takeaway
Every bound holds with about 62× slack. The true sup of |∆φ_R| is 43.49R⁻², against 2880R⁻² printed and C∗ = 3000. The two transition intervals of A∗ are disjoint, as claimed. The cutoff is exactly C² and not C³.

### Cited Findings
- (paper) P(s) = 10s³ − 15s⁴ + 6s⁵ on [0, 1], extended by 0 and 1, "is C², takes values in [0, 1], and satisfies |P″| ≤ 360". (The layout text extraction drops the exponent of s³ onto the next line; `pdftotext -raw` reads 10s³.)
- (paper) A∗(s) = P(2(1 + s))P(2(1 − s)) and φ_R(a) = A∗(a₁/R)A∗(a₂/R), which equals 1 on |a|∞ ≤ R/2 and 0 for |a|∞ ≥ R. "The transition intervals of the two factors in A∗ are disjoint; on each, the other factor equals one. Consequently |A∗″| ≤ 1440", and ‖∆φ_R‖∞ ≤ 2880R⁻² ≤ C∗R⁻².

### Inferences
- (computed) P(0) = P′(0) = P″(0) = 0, P(1) = 1, P′(1) = P″(1) = 0, and P′ = 30s²(1 − s)² ≥ 0, so P ∈ [0, 1] (exact rationals). P‴(0) = P‴(1) = 60, so the extension is C² and not C³.
- (computed) P″ = 60s(1 − s)(1 − 2s), with max |P″| = 10/√3 = 5.774 at s = ½ ± 1/(2√3). The paper's 360 bounds the coefficients term by term (60 + 180 + 120).
- (computed) The transitions are s ∈ [−1, −½] for the first factor and [½, 1] for the second. They are disjoint, and on each the other factor is identically 1 (the cross term P′P′ vanishes everywhere). So A∗″ = 4P″(·), max |A∗″| = 40/√3 = 23.09, consistent with 4·360 = 1440.
- (computed) sup |∆φ_R|·R² = 43.49, attained where both coordinates are in transition: 8|P″(u)P(u)| at u ≈ 0.809. The paper's 2·1440 is the sum of the two separate maxima. The true value is below the 2·23.09 = 46.19 that the same argument would give with exact P″.

### Gaps
- None.

## 7. Eqs. 25 and 39: the expanding-array scales and the loss budget

### Takeaway
The T_n/R_n² identity is exact, and δ < 1/24 holds for every tested (ν, b, N, m). The first half-budget is ≈ (ν/(1+ν))/48, so it approaches its 1/48 limit as ν → ∞. The second is negligible, below 10⁻¹⁰. The detection argument needs only δ < ½, so the 1/24 is a deliberate 12× margin.

### Cited Findings
- (paper) C∗ = 3000, Λ = 4D, R₀ = 1024C∗(1 + ν)Λ(K₀D + 2), R_n = R₀Λⁿ, S_n = 16R_n, T_n = S_{n+1}(K_{n+1} + 2), t₀ = 1, t_{n+1} = t_n + T_n (eq. 25); K_n = NB_n² = K₀Dⁿ, D = b² (eq. 23).
- (paper) δ = νC∗(R₀⁻² + Σ T_n/R_n²) < 1/24 (eq. 39), with T_n/R_n² = (16Λ/R₀)(K₀D(D/Λ)ⁿ + 2Λ⁻ⁿ). "Here D/Λ = 1/4 and Λ ≥ 4, so each geometric series has sum at most 4/3", giving νC∗Σ ≤ (ν/(1 + ν))·64/(3·1024) < 1/48. Also νC∗/R₀² ≤ ν/(1024²C∗(1 + ν)²) < 1/48.
- (paper) "Rn+1 ≥ 2Rn", "In particular Rn, Tn ≥ 1, and tn ≥ n + 1."

### Inferences
- (computed) The identity T_n/R_n² = 16ΛR_n(K₀D^{n+1} + 2)/R_n² = (16Λ/R₀)(K₀D(D/Λ)ⁿ + 2Λ⁻ⁿ) holds for independent Λ and D. It was checked in exact rationals at 200 random points for n = 0…11.
- (computed) In closed form, ΣT_n/R_n² = (16Λ/R₀)(4K₀D/3 + 2Λ/(Λ − 1)). This matches 60-term partial sums to 10⁻³⁰ relative.
- (computed) Over 160 parameter sets (b ∈ {2, 3, 5, 16}, N ∈ {2, 5, 100, 10⁴}, m ∈ {1, 3}, ν ∈ {10⁻⁶, 10⁻², 1, 10², 10⁶}), each half is < 1/48 and δ < 1/24, in exact rationals. The worst δ is 0.020833, at ν = 10⁶, within 0.1% of (ν/(1 + ν))/48.
- (computed) 64/(3·1024) = 1/48 exactly. So the first strict inequality comes from ν/(1 + ν) < 1, and from Σ Λ⁻ⁿ = Λ/(Λ − 1) ≤ 16/15 < 4/3.
- (computed) Λ = 4b² ≥ 16, so R_{n+1} = ΛR_n ≥ 16R_n. The paper's "≥ 2" and "Λ ≥ 4" are true but weak.
- (computed) The smallest instance (b = 2, N = 2, m = 1, ν = 1) gives K₀ = 8, D = 4, Λ = 16 and R₀ = 3 342 336 000.
- (computed) What δ must actually satisfy: event (3) needs 1 − δ > ½ and δ < ½. The "more than 3/4 of the unit mass" remark in §5.4 needs δ < ¼. Both have wide slack against 1/24.

### Gaps
- The parameter sweep is finite. The general proof is the two lines above, and the identity is exact for all parameters.

## 8. Eqs. 22–24: stack encoding, inverses and the address map

### Takeaway
All updates, inverses, range bounds and the address map are correct. Injectivity was confirmed exhaustively on 72 random tables that satisfy the compiler's incoming-rule property. The left-row bound y′ ≤ bB_n − 1 is attained.

### Cited Findings
- (paper) Update table (eq. 22): R: x′ = bx + d_β, y′ = ⌊y/b⌋. N: x′ = x, y′ = y − d_α + d_β. L: x′ = ⌊x/b⌋, y′ = d_γ + bd_β + b²⌊y/b⌋ with d_γ = x mod b.
- (paper) Inverses: ((x′ − d_β)/b, by′ + d_α), (x′, y′ − d_β + d_α) and (bx′ + d_γ, b⌊y′/b²⌋ + d_α). B_n = b^{m+n}; "y′ ≤ (b − 1) + b(b − 1) + b²(B_n/b − 1) = bB_n − 1"; address k = 1 + x + B_ny + B_n²(i − 1), 1 ≤ k ≤ K_n (eq. 24).

### Inferences
- (computed) For b ∈ {2, 3}, N ∈ {2, 3, 4}, levels n = 0, 1 and 6 random tables each, every configuration with a defined rule was enumerated. Each table has one move per target state, no rules out of the terminal state, and distinct (i′, β) pairs. In every case:
  - x′, y′ < B_{n+1};
  - the written digit is recovered from x′ mod b, y′ mod b or ⌊y′/b⌋ mod b by row;
  - (i′, β) identifies the rule, and the printed inverse returns (x, y);
  - k ∈ [1, K_n], k′ ∈ [1, K_{n+1}], and k → k′ is injective;
  - integer division recovers (x, y, i) from k.
- (computed) The left-row maximum equals bB_n − 1 for b = 2, 3, 4. The arithmetic is exact and the bound is attained.
- (computed) Negative control: two rules sharing (i′, β) map different configurations to the same image. The incoming-rule property from [12, Lemma 3.1] is necessary.

### Gaps
- The compiler that guarantees the incoming-rule property is in [12] and was not opened.

## 9. The gate (eq. 28) and the center paths (eqs. 26–27)

### Takeaway
W is divergence free, equals G on |ξ|∞ ≤ 2R, and vanishes beyond 3R. The paths start and end where claimed, stay S_n = 16R_n apart in one coordinate, and nonterminal routes stay at or below −S_n.

### Cited Findings
- (paper) χ(a) = Π_ℓ θ(a_ℓ + 3)θ(3 − a_ℓ) and W = ∇⊥[χ(ξ/R)(G₁ξ₂ − G₂ξ₁)] with ∇⊥ = (∂₂, −∂₁). It "is divergence free, equals G when |ξ|∞ ≤ 2R_n, and vanishes when |ξ|∞ ≥ 3R_n".
- (paper) Θ(r) = θ(2r − ½) "makes its transition between 1/4 and 3/4"; the paths are eq. 26 with Y_k = −(k + 2)S_n; c_{n,k,2} ≤ −S_n for nonterminal targets (eq. 27); |Y_k| ≤ S_n(K_n + 2) ≤ T_n.

### Inferences
- (computed) θ(a + 3)θ(3 − a) = 1 exactly for |a| ≤ 2 and 0 for |a| ≥ 3, so χ = 1 on |a|∞ ≤ 2 and χ = 0 on |a|∞ ≥ 3.
- (computed) The analytic W at R = 1.7 and G = (0.7, −2.2) on a 6561-point grid: W = G to 10⁻¹⁵ on the plateau, and W = 0 exactly beyond 3R. The central-difference divergence is at most 3.6·10⁻⁹, which is truncation error, since div ∇⊥ψ ≡ 0.
- (computed) Θ is 0 for r ≤ ¼ and 1 for r ≥ ¾.
- (computed) For random injective tables (b = 2 with 6 gates, b = 3 with 18 gates), sampled at 601 times across the stage:
  - c(t_n) = (S_nk, −S_n) and c(t_{n+1}) = (S_{n+1}k′, s_{i′}S_{n+1});
  - every pair of centers is at least S_n apart in some coordinate, which exceeds the 6R_n needed for disjoint gate supports;
  - eq. 27 holds.
- (computed) |Y_k| ≤ S_n(K_n + 2) ≤ T_n holds.

### Gaps
- Eq. 30's derivative bounds C_j were not computed. They are finite, but their values are not stated in the paper.

## 10. Eq. 21 viscosity rescaling and the eq. 5–6 reduction

### Takeaway
Both are exact identities, verified by exact polynomial algebra with ν as a symbol.

### Cited Findings
- (paper) U_ν = (νa(νt, Y), w(νt, Y)) and f_ν = (ν²F(νt, Y), νh(νt, Y)) (eq. 21). "Every horizontal equation acquires the factor ν², and every vertical term the factor ν."
- (paper) With F = ∂_ta + (a·∇)a − ν∆a and w solving eq. 6, U = (a, w) with p = 0 solves eq. 1 (eqs. 5–6).

### Inferences
- (computed) A polynomial stream function ψ(t, x, y) and scalar w(t, x, y) were chosen, both vanishing at t = 0, and F and h were formed at viscosity 1. The variables were (ν, t, x, y) with rational coefficients. The residuals of the viscosity-ν equations for (U_ν, f_ν), horizontal and vertical, are the zero polynomial, and div U_ν = 0. Since ν is a symbol, this holds for every ν.
- (computed) Negative control: f_ν = (νF(νt), …) leaves a nonzero residual.
- (computed) The eq. 5–6 reduction leaves a zero residual for (a, w), and div_Y a = 0.

### Gaps
- No sympy is installed, so the check uses a small exact polynomial class in the script instead of a CAS.

## 11. Eq. 8's R = 4^L, and why eq. 16 holds

### Takeaway
Eq. 16 needs e^{Ln} ≤ Rⁿ, that is e^L ≤ R. Eq. 8 defines **R = 4^L**: in the PDF the "L" sits raised and smaller (yMax 384.8 against the "4"'s 389.9). Since e < 4, e^L ≤ 4^L for every L, so eq. 16 holds. Plain-text extraction flattens the superscript to "R = 4L", under which e^L ≤ 4L would fail from L = 3 (e³ = 20.09 > 12). The build first reported that as a slip; the PDF's glyph positions show it is not one. *(Computed here; the PDF was checked with pdftotext's bounding boxes.)*

### Inferences
- Values that depend on R (b_n, D_n, M₀) use R = 4^L; at L = 1 the two readings coincide.

## 12. What this changes for the implementation

### Takeaway
Nothing in the paper needs correcting before implementation. The checks change five implementation choices.

### Inferences
- **C = 41.** The paper leaves C implicit. The `navier_stokes_computation` example should use C = 41, or compute and assert it.
- **The paper-scale torus construction cannot be simulated.** It has bump widths b₁ = 1.25·10⁻⁷ and δ₁ ≈ 10⁻²⁰, shrinking by 8 and 4096 per block. The §5 radii start at R₀ ≈ 3.3·10⁹. The planned small numerical run must use relaxed parameters and assert the same inequalities with the true constants:
  - sup K = 1/(πe d²) ≈ 120, not H;
  - the actual Lipschitz bound of the toy planar field, not 4^L;
  - ∫g = 2.25b².

  For example, a trivial V with b_n ≈ 0.02·2⁻ⁿ satisfies √2b_n < 1/32 and 2·120·M₀ + 1/16 < ½. That run checks the mechanism, not the paper's literal numbers.
- **Floating-point schedules.**
  - Even for the smallest table, T_n overflows u128 after about 14 levels and R_n after about 24, so use exact rationals or log-scale f64 for display.
  - In f64, assert δ < 1/24, not each half < 1/48. The first half is (ν/(1+ν))/48 times (4K₀D/3 + 2Λ/(Λ−1))/((4/3)(K₀D + 2)) ≤ 1. When ν and K₀D are both large, that is within f64 rounding of 1/48, and a strict comparison can fail spuriously.
- **φ_R is only C².** P‴ jumps at both ends. Test ∆φ_R analytically, not by spectral differentiation, which would ring. Assert the paper's ≤ 2880R⁻², and optionally the true 43.49R⁻².
- **Tests can be sharper than the paper.**
  - `torus_heat_kernel`: assert K ≤ 1/(πe d²)·(1 + 10⁻⁹) on dist ≥ d, τ ≤ 1, and K ≤ 1 + 10⁻¹⁵ for τ ≥ 1, as well as the eq. 14 bounds.
  - The R² derivative kernel: assert ‖∂_ℓH_s‖₁ = (πνs)^{−1/2} exactly.
  - `forcing` and the rescaling: the polynomial identities above give exact test cases.
- **Notation clashes to avoid in code.**
  - g is both the §3 bump and the §5 vertical source (eq. 31).
  - b is the smooth step b(r), the widths b_n and the alphabet size |Σ|.
  - D is both D_n (the Laplacian bound) and D = b².
  - δ_n (block duration) differs from δ (loss budget), and θ_j from θ.

### Gaps
- The companion preprints [10–12] were not opened. They supply the planar processor's corridor, the derivative bounds (8) and the injective compiler, and are out of scope for this note.
- The analytic steps (Lemma 3.2's Galerkin argument, Lemma 5.1's mild solution, the uniqueness energy estimate eq. 7) are standard. They were read, but not re-proved here.
