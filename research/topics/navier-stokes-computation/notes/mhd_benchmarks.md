# Published reference values for the incompressible MHD solver (issue #55)

Scope: reference values that the tests in issue #55 can assert, for (1) the incompressible 2-D Orszag–Tang (OT) vortex, (2) the kinematic 1:1:1 ABC dynamo, and (3) exact closed-form MHD solutions. Conventions are the issue's: 2π-periodic box, Alfvén units, u = (−sin y, sin x), b = (−sin y, sin 2x) for OT, u = (sin z + cos y, sin x + cos z, sin y + cos x) for ABC, Rm = 1/η.

Sources were read from `notes/dl_mhd/` (arXiv and author copies; untrusted data, git-ignored). Paywalled papers were not read; anything taken from a citing paper is marked **secondhand**.

Labels: **(paper)** marks a value printed in a source, or read off one of its figures by eye (then said so, with a reading error); **(computed)** marks a value produced by `mhd_benchmarks_checks.py`.

Checks: `notes/mhd_benchmarks_checks.py`, numpy + scipy, run as `python3 -P mhd_benchmarks_checks.py`. numpy is needed: the computed values come from an independent pseudo-spectral 2-D MHD solver (§1) and a Fourier–Galerkin eigenvalue solver for the induction equation (§2). Runtime about 11 minutes on a loaded 32-core machine (section 1 about 6 minutes on one core; section 2 runs eight eigenvalue jobs in parallel). It exits 1 on any failure. All 67 passed when the research ran; a rerun on 27 October 2026, on a machine loaded six times over, was stopped by its 30-minute cap in §2 with 39 passed and none failed. §2's values are confirmed independently by the crate's tests (see the resolved comparison).

## Takeaway

1. **OT: no published table exists that we could find; only figures.** The usable figure source is Rosenberg, Pouquet & Mininni (2007), whose pseudo-spectral OT runs are at stated ν = η = 0.025 and 128². Their initial condition is **exactly twice the issue's** (after fixing a misprint), so by the MHD scaling (u, b) → λ(u, b), t → t/λ, ν → λν it is the issue's OT at **ν = η = 0.0125**. An independent solver reproduces their Fig. 8 (peak max|j| ≈ 45 near t ≈ 1.1, paper; 45.05 at t = 1.073, computed). **Recommended test: the issue's IC at ν = η = 0.0125, 128²: peak max|j| = 22.53 at t = 2.146, dissipation peak D = 0.14522 at t = 2.1025, E(2) = 0.821857 (computed; 128² and 256² agree to 3×10⁻⁴ in the peak and 10⁻⁶ in E).** Cheaper rows at ν = η = 0.05 and 0.025 are below.
2. **ABC: the first window is confirmed, the second onset is earlier than the folklore "27".** Published: 8.9 ≤ Rm ≤ 17.5 and Rm ≥ 27 (Galloway & Frisch 1986, secondhand via Teyssier et al. 2006 and Archontis et al. 2003); Dorch (2000) prints 8.5–17.5. Computed: σ = 0 at **Rm = 8.934** and **Rm = 17.554**, and the second window opens at **Rm = 25.55** (Bouya & Dormy's Fig. 1 puts that zero crossing near 25 too). The mode crossing that Bouya & Dormy place at Rm₁ ∈ [24.05, 24.10] (paper) is reproduced (computed: between 24.05 and 24.10), with the pulsation jumping 0.5385 → 0.1308 (paper: "ω ≃ 0.53" → "ω ≃ 0.13").
3. **The first-window eigenvalue is a complex pair, not real and not a degenerate real pair** (paper: Bouya & Dormy, "the dominant eigenvalues are complex, leading to oscillations of the energy"). Computed at Rm = 12: **λ = 0.00807 ± 0.6177 i** for the field amplitude, so the magnetic energy grows at 2σ = 0.0161 with a superposed oscillation of period π/ω = 5.09. Bouya & Dormy's Fig. 1 gives σ_max ≈ 0.008 near Rm ≈ 11 (paper, read by eye).
4. **A fast dynamo test with a real growth rate:** the Roberts-type member A = 0, B = C = 1 of the same family. Computed: real, doubly degenerate λ, onset at Rm = 6.215, and λ = 0.08289 at Rm = 10. No published number for exactly this box and normalisation was found; Alexakis (2011) says qualitatively that these 2½-D flows are the first in the ABC family to dynamo in the 2π box, consistent with the computed onset.
5. **Exact solutions beyond Walén waves and Beltrami pairs:** the Alfvénic (Elsässer) states u = ±b with no mean field (Chandrasekhar 1956) and, in the kinematic problem, b = u for any Beltrami u (λ = −η|k|² exactly; computed).

## 1. Orszag–Tang vortex

### Cited findings (paper)

**Rosenberg, Pouquet & Mininni (2007)**, "Adaptive mesh refinement with spectral accuracy for magnetohydrodynamics in two space dimensions", [New J. Phys. 9, 304](https://doi.org/10.1088/1367-2630/9/8/304) (checked on Crossref); arXiv:physics/0703024 (`dl_mhd/physics_0703024.pdf`). Spectral-element code GASpAR against "a well–characterized pseudo–spectral code" (their refs. 13, 14).
- IC as printed (§4): Ψ = 2α[cos(2πx) + cos(2πy)], a_z = 2cos(2πx) + cos(2πy), α = 1, on [0, 1]²; "E_V and E_M, both equal to 2"; "initial correlation coefficient ρ̃_c is equal to 41%", ρ̃_c = 2u·b/(u² + b²).
- Dissipation D_T = −ν⟨ω²⟩ − η⟨j²⟩ (eq. 22).
- §5, first series: "ν = η = 0.025 (R_v = 80π)", pseudo-spectral comparison run at 128²; "the peak in total dissipation is reached near t = 1.1".
- Fig. 8 (128² series), read by eye from the black pseudo-spectral curve: J_max(0) ≈ 6; peak ≈ 45 at t ≈ 1.1 (inset: ±1, ±0.03); ≈ 29 over t = 1.5–1.6; minimum ≈ 10 near t ≈ 2.35; ≈ 23 at t = 2.5.
- Fig. 9 (256² series, ν not printed), by eye: ≈ 71 at t ≈ 0.95, peak ≈ 76.5 at t ≈ 1.1, ≈ 63 at t ≈ 1.65. Fig. 10: total energy 4 at t = 0 falling to ≈ 2.7 at t = 2.5; ⟨ω² + j²⟩/2 (as labelled) ≈ 14 at t = 0, peak ≈ 82 near t = 1.1.
- §4: enstrophy maximum "at t ≈ 0.85" for the fiducial R_v = 10³π run (512² equivalent), Fig. 7: J_max peak ≈ 180 (by eye). Too fine for our tests.
- **Only figures; no table of J_max or D_T values.**

**Biskamp & Welter (1989)**, Phys. Fluids B 1, 1964 (doi 10.1063/1.859060): paywalled, not read. Its IC, **secondhand** from Flint & Vahala (arXiv:1706.01424, eqs. 24–25): u = U₀[sin(y + 0.5) x̂ − sin(x + 1.4) ŷ], B = B₀[sin(y + 4.1) x̂ − 2 sin(2x + 2.3) ŷ], a phase-shifted OT, not the issue's IC. The citing paper compares only with B&W's Fig. 20 (dissipation rate), qualitatively.

**Politano, Pouquet & Sulem (1989)**, Phys. Fluids B 1, 2330; **Orszag & Tang (1979)**, J. Fluid Mech. 90, 129; **Dellar (2002)**, J. Comput. Phys. (citation from memory, not checked): paywalled or not found open; not read. No secondhand numbers for them were found.

### Inferences (computed)

**The printed a_z is a misprint, and the paper's IC is 2× the issue's.**
- The literal a_z = 2cos x + cos y (2π units) gives max|j|(0) = 3, not Fig. 8's 6; ruled out.
- a_z = 2cos x + cos 2y and a_z = cos 2x + 2cos y both give max|j|(0) = 6, E(0) = 4 (E_V = E_M = 2), and identical E, D and max|j| histories (to 4×10⁻¹⁴). The second is exactly 2 × (the issue's ψ = cos x + cos y, a = cos y + ½cos 2x).
- The paper's units are effectively a 2π box with these amplitudes: run that way, Fig. 8 is reproduced (table below). R_v = 80π at ν = 0.025 means U₀L₀ = 2π.
- ρ̃_c: every reading gives 2⟨u·b⟩/⟨u² + b²⟩ = 0.5, not the printed 41%. Not reconciled; it does not affect the dynamics.

**Fig. 8 reproduced** (amplitude-2 IC, ν = η = 0.025, 128², 2/3 rule, IF-RK4, dt = 2.5×10⁻³):

| Quantity | Fig. 8 / text (paper) | Computed |
|---|---|---|
| max\|j\| at t = 0 | ≈ 6 | 6 |
| peak max\|j\| | ≈ 45 at t ≈ 1.1 | 45.05 at t = 1.073 |
| max\|j\| at t = 1.5, 1.6 | ≈ 29 | 29.01, 29.49 |
| max\|j\| at t = 2.5 | ≈ 23 | 23.50 |
| peak of D | "near t = 1.1" | 1.1617 at t = 1.051 |

**Fig. 9's unprinted viscosity is ν = η = 0.01** (R_v = 200π): that run gives local maxima 71.52 at t = 0.954, 76.51 at 1.114 and 63.50 at 1.671, and a peak ⟨ω²⟩ + ⟨j²⟩ = 83.27, matching the figure. ν = 0.0125 gives 61.43, 67.97 and 54.07, which does not. Fig. 10's "/2" in the enstrophy label is inconsistent with its t = 0 value: ⟨ω²⟩ + ⟨j²⟩ = 4 + 10 = 14 without the /2.

**Reference values in the issue's units** (u = (−sin y, sin x), b = (−sin y, sin 2x); E = ½⟨u² + b²⟩, D = ν⟨ω²⟩ + η⟨j²⟩ = −dE/dt, ⟨·⟩ the box mean, max|j| over grid points). At t = 0: E = 1, ⟨ω²⟩ = 1, ⟨j²⟩ = 5/2, max|j| = 3.

| ν = η | Grid | Peak max\|j\| | Other local maxima of max\|j\| | Peak D | E(1) | E(2) | E(4) |
|---|---|---|---|---|---|---|---|
| 0.05 | 128² (64²: 7.5138) | 7.5144 at t = 2.0226 | 4.853 at 1.150 | 0.22001 at 1.1494 | 0.813019 | 0.599845 | 0.354156 |
| 0.025 | 128² (64²: 13.669) | 13.7375 at t = 2.0884 | 7.855 at 1.297; 9.804 at 1.757 | 0.17857 at 1.9549 | 0.895126 | 0.726973 | 0.499248 |
| **0.0125** | 128² | **22.527 at t = 2.1458** | 18.15 at 1.835; 14.75 at 3.229 | **0.14522 at 2.1025** | 0.943781 | 0.821857 | 0.618952 |
| 0.0125 | 256² | 22.534 at t = 2.1459 | 18.27 at 1.839; 14.82 at 3.228 | 0.14522 at 2.1025 | 0.943781 | 0.821857 | 0.618952 |

For ν = η = 0.0125 also E(3) = 0.697867 and E(5) = 0.559084. The ν = 0.0125 row is Rosenberg et al.'s Fig. 8 in the issue's units (times ×2, j ÷2): 2 × 22.527 = 45.05.

Notes for the test author:
- **Energies converge far faster than max|j|.** At ν = 0.0125, E agrees between 128² and 256² to 1.5×10⁻⁷ over t ≤ 5, but grid-point max|j| differs by up to ~0.5 % at some times (8.538 vs 8.551 at t = 1) because the maximum falls between grid points. The peak itself sits on a symmetric grid point and agrees to 3×10⁻⁴. This is the `gradient_sup` caveat of the issue, seen in practice.
- Assert E(t) and the peak of D to ~10⁻⁵; assert the peak of max|j| to ~10⁻³ at 128²; avoid asserting max|j| at arbitrary times.
- The global maximum of max|j| is a close contest between local maxima at ν = 0.05 (4.85 vs 7.51, safe) and ν = 0.0125 (18.15 vs 22.53, safe); at ν = 0.025 there are three. Assert the global peak, not "the first peak".

## 2. Kinematic ABC dynamo, A = B = C = 1

Conventions in all sources below: 2π-periodic box, field of the same periodicity as the flow, η = 1/Rm with the issue's u (amplitude 1 per Beltrami wave, ⟨u²⟩ = 3). Growth rates are of the **field** (b ∝ e^{λt}); the energy rate is 2σ.

### Cited findings (paper)

- **Teyssier, Fromang & Dormy (2006)**, J. Comput. Phys.; arXiv:astro-ph/0601715, §4.3: "at least two regions of instability ..., one for 8.9 ≤ Rm ≤ 17.5 and a second for Rm ≥ 27 (see Galloway and Frisch, 1986)". Same u at 1:1:1 (eq. 54); "Rm = (UL)/η" with U "some norm of the velocity", and their windows coincide with Galloway & Frisch's. Their finite-volume growth rates agree with Galloway & Frisch's spectral ones in Fig. 8 (figure only). At Rm = 159: growth rate 0.055 (256³ Cartesian) and 0.052 (AMR). Out of reach for our grids.
- **Archontis, Dorch & Nordlund (2003)**, A&A 397, 393; arXiv:astro-ph/0204208, §2: "two windows of dynamo action ...: A closed one in the interval Rem = 8.9–17.5, and an open one for magnetic Reynolds numbers above 27 (Galloway et al. 1984)"; §3: the growth rate "reaches its maximum value at Rem = 12 ... a maximum growth rate of the magnetic field in this window (Galanti et al. 1993)" (secondhand).
- **Dorch (2000)**, Physica Scripta 61, 717; arXiv:astro-ph/0007466: "Rem = 8.5 – 17.5 ([9] = Arnold & Korkina 1983) and one starting at Rem ≈ 27 ([10], [11] = Galloway & Frisch 1984, 1986)". **The lower edge differs (8.5 vs 8.9) between citing papers.** Dorch also gives the first-window seed b = (sin z − cos y, sin x − cos z, sin y − cos x).
- **Bouya & Dormy (2013)**, Phys. Fluids 25, 037103; arXiv:1206.5186 and the published version (`bouya_pof.pdf`); u as eq. 2, η = Rm⁻¹, spectral, N = 64 at small Rm:
  - mode crossing Rm₁ ∈ [24.05, 24.10], "in the stable window reported by Galloway and Frisch [1986] near Rm = 20";
  - "the dominant eigenvalues are complex, leading to oscillations of the energy"; at Rm₁ the pulsation "jumps discontinuously from ω ≃ 0.53, corresponding to the first window ..., to ω ≃ 0.13 corresponding to the second window";
  - eigenvalue coalescence Rm₂ ∈ [215.0, 215.4]; above it λ is real;
  - Fig. 1 (σ = Re λ vs Rm), read by eye (±0.001 in σ, ±3 % in Rm): σ ≈ −0.008 at Rm ≈ 7.7; +0.0035 at 9.3; **maximum ≈ 0.008 at Rm ≈ 11.2**; +0.006 at 13.5; +0.0016 at 16.1; −0.003 at 19.4; ≈ −0.005 at 24; zero crossing of the second-window branch at **Rm ≈ 25**; ≈ +0.02 at Rm ≈ 33.
  - **No table**: Fig. 1 and Fig. 3 (ω vs Rm) only.
- **Alexakis (2011)**, Phys. Rev. E 84, 026321; arXiv:1105.3692. Different normalisation: R_M = U/(ηk_u), U = ⟨u·u⟩^{1/2} (eqs. 3–5), so for the issue's 1:1:1 flow (U = √3, k_u = 1) **R_M = √3 Rm**, and γ/(Uk_u) = σ/√3. 2π box (k_uL = 2π), 32³ for R_M ≲ 20, figures only. Statements for 1:1:1: dynamo at R_M = 20 and 25; not a dynamo at R_M = 33.3 ("the no-dynamo window"); dynamo again at R_M = 100. The 2½-D flows (one coefficient zero, the other two equal) "result first in dynamo".
- **Galloway & Frisch (1984, 1986)**, GAFD 29, 13 and 36, 53; **Arnold & Korkina (1983)**; **Galanti, Sulem & Pouquet (1992)**, GAFD 66, 183; **Lau & Finn (1993)**; **Jones & Gilbert (2014)**, GAFD 108(1): paywalled, not read. The window edges above are secondhand from them.

### Inferences (computed)

Fourier–Galerkin operator for ∂ₜb = ∇×(u×b) + η∆b on |k_x|, |k_y|, |k_z| ≤ K, Leray-projected, k = 0 excluded, rightmost eigenvalues by ARPACK. K = 8 ↔ 17 modes a side (a 32³ grid with the 2/3 rule keeps |k_i| ≤ 10; 16³ keeps ≤ 5).

| Rm | λ (leading, field) | next | K used (convergence) |
|---|---|---|---|
| 5 | −0.08479 ± 0.61614 i | −0.2 (b = u) | 5, 8, 9 agree |
| 8 | −0.00823 ± 0.62761 i | −0.125 (b = u) | 9 |
| 9 | +0.000439 ± 0.62713 i | −0.1111 | 8 = 9 to 10⁻⁷ |
| 10 | +0.00529 ± 0.62518 i | −0.1 | 9 |
| **12** | **+0.00807 ± 0.61773 i** | −0.08333 (b = u); −0.11426 ± 0.11063 i | 8–12 agree; K = 5: 0.008069 ± 0.61854 i |
| 15 | +0.00441 ± 0.60015 i | −0.0667 (b = u); −0.0679 ± 0.1275 i | 9 |
| 17.5 | +0.000084 ± 0.58237 i | −0.0431 ± 0.1323 i | 10 = 11 to 2×10⁻⁷ (K = 8: 0.000065) |
| 20 | −0.00303 ± 0.56435 i | −0.02545 ± 0.13331 i | 9 (K = 8: −0.00309, −0.02529; K = 5: −0.0109) |
| 24.05 | −0.00548 ± 0.53853 i | −0.00566 ± 0.13082 i | 9 |
| 24.10 | −0.00536 ± 0.13080 i | −0.00547 ± 0.538 i (first-window branch) | 10 (K = 9: −0.00547) |
| 26 | +0.00154 ± 0.12871 i | — | 10 (K = 8: 0.00216; 9: 0.00136) |
| 27 | +0.00471 ± 0.12747 i | — | 11 (K = 8: 0.00549; 9: 0.00452; 10: 0.00475) |
| 30 | +0.0125 ± 0.1233 i | — | 9 only |

- **Window edges:** σ = 0 at Rm = 8.934 (K = 8) and 17.554 (K = 10) in the first window, and at Rm = 25.547 (K = 10; K = 11 moves it by about +0.02) for the second. Published "8.9" and "17.5" are confirmed to the precision they are printed; Dorch's "8.5" is not (σ(8.5) = −0.0033). The "27" for the second onset is not: σ(26) = +0.0015 at K = 10. Bouya & Dormy's Fig. 1 agrees with the earlier onset.
- **The maximum in the first window** is σ = 0.00807 near Rm = 12, matching both "maximum value at Rem = 12" (Archontis et al., citing Galanti et al.) and Bouya & Dormy's Fig. 1 (≈ 0.008).
- **Mode crossing:** the leading pulsation is 0.5385 at Rm = 24.05 and 0.1308 at Rm = 24.10, so Rm₁ ∈ [24.05, 24.10] is reproduced. The second-window branch converges slowly in K (it alternates with the parity of K); assert it only at K ≥ 10.
- **Alexakis's statements, mapped** (R_M = √3 Rm): R_M = 14.3 → Rm = 8.26, σ < 0; R_M = 20 → 11.55, σ > 0; R_M = 25 → 14.43, σ > 0; R_M = 33.3 → 19.23, σ < 0. All four consistent.
- **The b = u mode** is an exact eigenmode with λ = −η (Beltrami: u×u = 0). It is the second eigenvalue for Rm ≲ 15. A solver test: start from b = u and check e^{−ηt} decay to round-off.
- **Roberts-type flow** (A = 0, B = C = 1, u = (cos y, sin x, sin y + cos x)): real, doubly degenerate leading eigenvalue (the flow is z-independent). σ = −0.47474 (Rm = 2), −0.24964 (3), −0.05622 (5), −0.00820 (6), **onset Rm = 6.215**, +0.02539 (7), +0.08289 (10), +0.1374 (20), at K = 8. K = 6 is within 10⁻⁵ for Rm ≤ 10, and K = 10 matches to 3×10⁻⁵ at Rm = 20. Growth is non-oscillatory, so a 16³ run to t ≈ 100 at Rm = 10 measures σ directly.

### Comparison with the solver's runs (resolved)

`magnetics::navier_stokes`'s kinematic solver (`Torus::induction`), from a mean-free solenoidal seed, against the computed eigenvalues. σ is the field rate d ln‖b‖/dt; the early probes reported that same quantity, so there is no factor of 2.

| Flow, Rm | Grid | Solver σ | Solver π/ω | Eigen-solver |
|---|---|---|---|---|
| ABC, 12 | 16³ | +0.008069 | 5.0791 | K = 5: 0.008069 ± 0.61854 i (π/ω = 5.0791) |
| ABC, 12 | 32³ | +0.008070 | 5.0857 | converged: 0.00807 ± 0.61773 i (π/ω = 5.0857) |
| ABC, 5 | 16³ | −0.084789 | 5.0988 | −0.08479 ± 0.61614 i (π/ω = 5.0988) |
| ABC, b = u, 5 | 16³ | −0.200000 | — | −η exactly |
| Roberts, 10 | 32 × 32 × 4 | +0.082890 | — | +0.08289 |
| Roberts, 10 | 16 × 16 × 4 | +0.082718 | — | (0.2 % low: in-plane truncation) |
| Roberts, 5 | 16 × 16 × 4 | −0.05622 | — | −0.05622 |

- **The early "wander" was the estimator, not the solver.** ln‖b‖ oscillates by about ±0.5 in rate with period π/ω ≈ 5.1. Fixed windows that are not whole multiples of it leave an error of order (amplitude)/(window), about 0.005, which is the size of the drift seen. The slope of ln‖b‖ between the energy's peaks, a whole number of periods apart, recovers σ to 10⁻⁶ by t ≈ 80.
- **Both solvers agree:** the eigen-solver here (Fourier–Galerkin, numpy) and the time-stepping pseudo-spectral solver are independent implementations. The 16³ grid (modes to 5 per axis) reproduces the K = 5 truncation's eigenvalue to all printed digits.
- **The same holds for Orszag–Tang (§1).** On 64² and 128² the solver gives E(2) = 0.5998446, the dissipation peak 0.22001 at t ≈ 1.150 and max|j| 7.5143 at t ≈ 2.02 for ν = η = 0.05. For ν = η = 0.0125 it gives max|j| 22.527 at t ≈ 2.145, the dissipation peak 0.14522 at t ≈ 2.100 and E(1…5) to 10⁻⁷, all as in the tables above.
- These comparisons are the tests `orszag_tang_vortex`, `orszag_tang_against_rosenberg_pouquet_mininni_fig_8`, `roberts_flow_dynamo` and `abc_flow_dynamo` in `crates/magnetics/src/navier_stokes.rs`.

## 3. Exact solutions worth testing (brief)

- **Alfvénic states without a mean field** (Chandrasekhar 1956, "On the stability of the simplest solution of the equations of hydromagnetics", PNAS 42, 273; citation checked, text not read). With z± = u ± b, the equations are ∂ₜz± + (z∓·∇)z± = −∇P + ν∆z± when ν = η. So z⁻ = 0 (u = b, any solenoidal field) is an exact nonlinear solution: steady if ideal, heat-equation decay otherwise. The Walén wave about B₀ is the same statement shifted by the mean field. Test: a random solenoidal u = b stays u = b to round-off, and each mode decays as e^{−νk²t}.
- **Kinematic Beltrami eigenmode** (computed above): with u frozen and Beltrami, b = u decays as e^{−η|k|²t} exactly.
- Nothing else with a closed form for the periodic box was found that is not a special case of these or of the issue's existing Alfvén, free-decay and Beltrami-pair tests.

## Gaps

- **OT tables:** none found. Biskamp & Welter 1989, Politano, Pouquet & Sulem 1989, Orszag & Tang 1979 and Dellar 2002 are paywalled and were not read. Queries tried: "incompressible Orszag-Tang vortex pseudo-spectral benchmark maximum current density time viscosity resistivity table"; "incompressible MHD solver verification Orszag-Tang maximum current table t=1 spectral reference value ν=η=0.01" (extended); "Politano Pouquet Sulem 1989 ... Orszag-Tang current maximum time"; "Politano Pouquet Sulem 1989 Physics of Fluids B 1 2330 pdf hal"; "Biskamp Welter 1989 ... Orszag-Tang initial condition current maximum". The GHOST, Dedalus, Snoopy and TURBO code papers were not searched one by one. The only quantitative anchors are Rosenberg et al.'s figures, reproduced here.
- **Rosenberg et al.'s 41 % correlation** is not reproduced by any reading of their IC (0.5 computed).
- **ABC growth-rate tables:** none open. Galloway & Frisch 1986 (the source of 8.9 / 17.5 / 27) and Jones & Gilbert 2014 are paywalled; the Exeter repository copies refused automated access. The published first-window growth rates used here are read off Bouya & Dormy's Fig. 1. Queries tried: "ABC flow dynamo 1:1:1 8.9 17.5 magnetic Reynolds number window Arnold Korkina Galloway Frisch"; "ABC kinematic dynamo A=B=C=1 growth rate Rm = 12 first window maximum growth rate value"; "1:1:1 ABC flow kinematic dynamo eigenvalue table growth rate Rm 10 12 15 symmetry class Podvigina OR Jones OR Gilbert arXiv" (extended); "Jones Gilbert Dynamo action in the ABC flows using symmetries arXiv".
- **The second onset** is "27" in every citing paper, but both this computation (25.55) and Bouya & Dormy's Fig. 1 (≈ 25) put it lower. The original Galloway & Frisch 1986 value could not be checked.
- **Roberts flow:** no published growth rate for this exact flow, box and normalisation was found. Devlen, Brandenburg & Mitra 2013 (arXiv:1212.2626) treat Roberts IV, a different flow. G. O. Roberts 1972 (Phil. Trans. R. Soc. A; citation from memory, not checked) was not read. Queries tried: "Roberts flow kinematic dynamo growth rate table vertical wavenumber magnetic Reynolds number numerical values periodic"; "Roberts flow kinematic dynamo critical magnetic Reynolds number value k_z growth rate Plunian Rädler subharmonic arXiv".
