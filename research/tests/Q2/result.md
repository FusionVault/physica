# Q2: are the stated errors of COST-G's monthly fields right?

**Date:** 2026-10-24. **Status:** complete. This test came out of building the Kalman smoother
(M3).

![Are COST-G's stated errors right?](../../../figures/kalman.svg)

## Result

Each COST-G file states an error for every coefficient. Fitted to 24 years of fields, the data
vary more than those errors allow:

| Degrees | Scale needed on the stated errors (κ) | Signal to noise | Withheld months: rms of normalised errors |
|---|---|---|---|
| 2–10 | 1.25 | 9.6 | 1.05 |
| 11–30 | 2.0 | 2.8 | 1.03 |
| 31–60 | 2.5 | 0.4 | 1.08 |
| 1–60, all | | | **1.07** |

**What this means:**
- **The errors are about right only at the longest wavelengths.** Above degree 10, the month-to-month
  scatter is about 2 to 2.5 times what the files state. Anyone who weights months or propagates
  errors with these sigmas is too confident at the scales of river basins and smaller.
- **Scaled, they work.** With κ applied, the smoother predicts every tenth month, withheld from
  both the seasonal fit and the smoother, with normalised errors of rms 1.07. That is close to the
  ideal of 1, so its uncertainty maps can be trusted.
- **Above degree 30, a month is mostly noise** (signal-to-noise 0.4). That is why maps are
  smoothed over 300 km, and why the smoother pulls these coefficients strongly toward their
  seasonal cycle.

**Caveats:**
- **κ measures everything the model does not explain.** That includes any real signal that is not
  a seasonal cycle, a trend or a slowly wandering Ornstein–Uhlenbeck process. At high degrees the
  real signal is small, so there κ measures noise.
- **Coefficients are treated as independent.** Real GRACE errors are correlated between
  coefficients (the stripes), so κ describes each coefficient's variance, not those correlations.
- **Low degrees are special.** C̄20 and C̄30 come from laser ranging but keep the field file's
  errors, which overstate them (κ = 0.3 at degree 2). Degree 1 has no stated error, so it gets a
  robust estimate from its month-to-month changes.

## Method

The code is `gravviz kalman` ([`apps/gravviz/src/kalman.rs`](../../../apps/gravviz/src/kalman.rs)),
with the smoother in [`gravimetry::kalman`](../../../crates/gravimetry/src/kalman.rs).

```bash
cargo run --release -p gravviz -- kalman data/costg data/baseline data/TN-14_C30_C20_GSFC_SLR.txt data/TN-13_GEOC_CSR_RL0603.txt figures viewer/data
```

1. **Coefficients.** Take every coefficient of degree 1 to 60 against the 2004–2009 mean, after
   the TN-14 and TN-13 replacements, from all 259 months, 2002 to 2026.
2. **Seasonal fit.** Fit an intercept, a trend, and annual, semi-annual and 161-day terms by
   weighted least squares. What is left is the signal plus noise.
3. **Signal and noise model.** The signal is an Ornstein–Uhlenbeck process with size `sd` and
   correlation time `tau`. The noise is the stated error times κ.
4. **Fit per degree.** Choose `sd`, `tau` and κ for each degree, shared by all its orders, by
   maximum likelihood over a grid of candidates (12 × 8 × 11). The likelihood is exact and comes
   from the Kalman filter.
5. **Smooth.** The Rauch–Tung–Striebel smoother gives every calendar month, including the 35 with
   no field, a value and an uncertainty.
6. **Cross-validate.** Withhold every tenth observed month from both the seasonal fit and the
   smoother, predict it, and compute (observed − predicted) / √(prediction variance + (κ·σ)²).

**How the smoother is checked:**
- It matches the exact Gaussian-process solution (by brute-force linear algebra) to 1e-10, for
  the mean, the uncertainty and the likelihood.
- The fit recovers the settings of simulated data exactly.

**Products:** gap-free monthly maps of water height (300 km Gaussian, no destriping, because
the smoother suppresses the stripes itself), their uncertainty, and the mass flow between them.
All three are in the GRACE Mass Atlas under "Every month, gaps filled (Kalman)".
