# Q0: what makes a COST-G month noisy?

**Date:** 2026-10-24. **Status:** complete, as a control for W1 and E1.

![Noise in COST-G months](../../../figures/q0_noise.svg)

## Result

| Effect | Size | q | Meaning |
|---|---|---|---|
| Solar flux (F10.7) on degrees 61–90 | +30% per σ [+19, +42] | < 0.001 | Short-wavelength noise rises with solar activity, with mission phases held fixed. Our search found no earlier measurement of this. |
| Wide-pointing mode on the C₂₀ offset | +198% [+38, +544] | 0.03 | GRACE-FO's C₂₀ disagrees with laser ranging three times more since July 2023. The maps replace C₂₀ from laser ranging anyway. |
| GRACE transplant on the C₂₀ offset | +249% [+33, +813] | 0.04 | Matches the literature: C₂₀ and C₃₀ are sensitive to the accelerometer transplant. |
| GRACE transplant on degrees 61–90 | −35% | 0.04 | Lower, not higher (see below). |
| Wide-pointing on degrees 61–90 | −55% | < 0.001 | Lower, not higher (see below). |
| Anything on the open-ocean noise at 300 km | F10.7 +6% (p = 0.08); transplant +40% (q = 0.10) | ≥ 0.10 | **At the scale W1 and E1 use, no factor reaches q < 0.1.** |
| Geomagnetic activity (Ap), centre weighting | — | ≥ 0.17 | No effect anywhere. |

**What this means for the other tests:**

- **Solar maximum doesn't explain W1 or E1.** It overlaps the 2023–24 El Niño and the 2025
  Kamchatka window. It affects the short-wavelength coefficients, but not the 300 km noise those
  tests use.
- **The lower high-degree amplitudes** in the transplant and wide-pointing periods probably reflect
  how COST-G's contributors regularise the short wavelengths, not cleaner data. The metric measures
  what reaches the user.

**The noisiest months:**
- January–February 2015: GRACE battery management.
- September 2004: the near-repeat ground track.
- April 2023: the start of the wide-pointing transition.

`figures/q0_noise.csv` gives every month a noise weight, (median RMS / RMS)², for E1 and later tests.

## Method

The code is `gravviz q0` ([`apps/gravviz/src/q0.rs`](../../../apps/gravviz/src/q0.rs)). The metrics and
the model were fixed in the code before the first run.

```bash
cargo run --release -p gravviz -- q0 data/costg data/baseline data/TN-14_C30_C20_GSFC_SLR.txt data/TN-13_GEOC_CSR_RL0603.txt data/indices/Kp_ap_Ap_SN_F107_since_1932.txt figures
```

**Noise metrics** (259 months):

1. **Open-ocean residual RMS**, mm of water.
   - Field: degree ≤ 60, 300 km Gaussian, no destriping.
   - Cells: the 4,850 cells of the 2° grid that are within 60° of the equator and at least 3 cells
     from any land.
   - Each cell's trend, annual, semi-annual and 161-day terms are removed over the whole record.
2. **Degrees 61–90** against the 2004–2009 mean, mm of geoid.
3. **|C₂₀(COST-G) − C₂₀(laser ranging)|**, after removing each mission's median.

**Model:** ln(metric) regressed on these terms:
- mission-phase dummies: GRACE from 2011; GRACE transplant from September 2016; GRACE-FO;
  wide-pointing from July 2023;
- the standardised span means of daily F10.7 and ln Ap (GFZ Potsdam, CC BY 4.0);
- the effective number of contributing centres, 1/Σw², parsed from each file's `CMMNT` header.

Errors are Newey–West with 12 lags. Benjamini–Hochberg runs over the 15 tests of interest (5
factors × 3 metrics).

## Limits

- **F10.7 overlaps the mission phases.** The 2011–2015 solar maximum overlaps GRACE's late period,
  and the 2023–2025 maximum overlaps wide-pointing mode. The dummies absorb what they share, so the
  F10.7 effect is identified from variation within each phase.
- **The open-ocean metric includes real ocean-bottom-pressure variability**, so it bounds the noise
  from above.
