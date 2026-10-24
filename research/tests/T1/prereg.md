# T1 pre-registration: the Hunga Tonga Lamb wave in GRACE-FO laser ranging

**Registered:** 2026-10-24, as the git tag `prereg/T1-2026-10-24`.
**Data status:** no GRACE-FO Level-1B or AOD1B data for January 2022 had been downloaded or seen
by this project when this was registered. No published study was found that looks for the event
in GRACE-FO's inter-satellite ranging (KBR or LRI); see
[the atmospheric-waves research note](../research_notes/Gravitational%20waves%20and%20gravity%20beams/atmospheric_gravity_waves_from_satellites.md).

## Question

The eruption of Hunga Tonga–Hunga Ha'apai on 15 January 2022 at 04:14:45 UTC (20.546° S,
175.390° W) sent a Lamb wave around the world at 318.2 ± 6 m/s
([Wright et al. 2022](https://doi.org/10.1038/s41586-022-05012-5)). Far from the volcano,
barometers recorded a pressure jump of 1.3–2 hPa lasting 25–50 minutes
([Tien Shan station, arXiv 2303.04362](https://arxiv.org/abs/2303.04362);
[485 stations in China, Frontiers 2022](https://doi.org/10.3389/feart.2022.931545)).

A 1 hPa pulse is 10.2 kg/m² of air. Does GRACE-FO's laser ranging interferometer (LRI) feel the
gravity of that moving mass?

## Prediction

The forward model is `gravviz tonga` at this tag ([`apps/gravviz/src/tonga.rs`](../apps/gravviz/src/tonga.rs),
[`gravimetry::atmosphere`](../crates/gravimetry/src/atmosphere.rs)):

- The front is frozen while the pair crosses it: the pair moves at 7.08 km/s over the ground and
  the front at 0.318 km/s.
- Locally the front is a straight strip of surface mass with a Gaussian profile. Each satellite
  feels its horizontal pull, and the range acceleration is the difference of the along-track
  pulls: `ρ̈ = sin φ·[g_x(x_lead) − g_x(x_trail)]`, where φ is the angle between the track and the
  front.

![The Tonga Lamb wave as GRACE-FO would see it](../figures/tonga_lamb_wave.svg)

| Per hPa, one crossing | σ = 1.25 min | 2.5 min | 5 min | 10 min |
|---|---|---|---|---|
| Peak ρ̈ at 90°, nm/s² | 0.07 | 0.13 | 0.24 | 0.39 |
| SNR at 90° | 3.0 | 5.9 | 11.2 | 19.2 |
| SNR at 45° | 1.6 | 3.2 | 6.1 | 10.3 |
| SNR at 20° | 0.4 | 0.8 | 1.6 | 2.7 |

- **How the SNR is computed:**
  - It is a matched filter against the nominal pair noise, SuperSTAR accelerometers plus the
    GRACE-FO LRI ([Zingerle et al. 2024](https://doi.org/10.1029/2024EA003630)).
  - It scales linearly with the pulse amplitude.
  - 90% of SNR² lies between 1.2 and 7.6 mHz.
- **The ocean adds to the signal:**
  - The pulse travels faster than ocean long waves (√(gH) ≈ 200 m/s), so the sea rises under it
    instead of sinking.
  - This adds 0.63× (4 km deep) to 0.94× (5 km) of the air's own mass
    ([Proudman 1929](https://doi.org/10.1111/j.2044-8312.1929.tb05286.x)).
- **What this means:**
  - Against instrument noise alone, a perpendicular crossing of a 1.5 hPa, σ ≈ 5–10 min pulse
    over deep ocean gives SNR ≈ 25–55.
  - The real limit is the background: weather and ocean signals at degrees 7–43 are removed only
    by AOD1B, so the test has to survive model residuals, not instrument noise.

## Data

| Input | Version | Notes |
|---|---|---|
| GRACE-FO Level-1B | JPL RL04, `GRACEFO_L1B_ASCII_GRAV_JPL_RL04` (PO.DAAC), 13–17 January 2022 | LRI1B, GNV1B, ACT1B, SCA1B, KBR1B. The LRI files for 17 and 18 January are empty (3 kB in the catalogue), so 17 January is KBR only |
| AOD1B | GFZ RL06, `GRACE_AOD1B_GRAV_GFZ_RL06` (PO.DAAC), `AOD1B_2022-01_06` | RL07 is at GFZ only, which was unreachable on 2026-10-24. RL06 uses operational ECMWF analyses every 3 h |
| Monthly field | COST-G RL02.1, January 2022 | already in this project |
| Static field | GOCO06s, degrees 91–300 | already in this project |
| Ocean tides | FES2014b as in the COST-G standards, or EOT20 if FES2014b's coefficients are not openly available | chosen on availability before any post-eruption residual is computed, and recorded in the results |

The SHA-256 of every input is recorded in the results when the pipeline is frozen.

## Blinding

The residual pipeline is written after this registration. It is developed and checked only on
the residuals before the eruption: 13 January 00:00 UTC to 15 January 04:14 UTC.

- **Freezing:** when the pipeline, the noise model and the null are done, the code is tagged
  `prereg/T1-frozen`.
- **No peeking:** no residual after the eruption is plotted, filtered or summarised before that
  tag.
- **Disclosure:** any deviation from this document is listed in the results with its reason.

## Method

1. **Residual range acceleration.**
   - Start from the LRI range acceleration with its light-time and antenna-offset corrections
     (LRI1B).
   - Subtract the modelled range acceleration at the GNV1B positions and velocities:
     - gravity of the January 2022 COST-G field (to degree 90), with GOCO06s for degrees 91–300;
     - solid-Earth and pole tides (IERS 2010) and ocean tides;
     - AOD1B RL06 (the `GLO` field, linear in time between the 3-hourly epochs);
     - the Sun and Moon;
     - the non-gravitational accelerations (ACT1B rotated with SCA1B);
     - the velocity term `(|Δv|² − ρ̇²)/ρ`.
   - Per revolution, remove a least-squares fit of a bias, a trend, and 1/rev and 2/rev sines and
     cosines.
   - Band-pass 1–10 mHz with a zero-phase 4th-order Butterworth filter run forward and backward.
     The same filter is applied to every template.
2. **Crossings.**
   - The front's angular distance from the volcano is `θ_f(t) = 318.2·(t − t₀)/R`, with R =
     6371 km. Past the antipode the front returns, so it is folded into [0, π].
   - A crossing is a time at which a satellite's angular distance from the volcano equals θ_f.
   - Crossings with θ < 10° (near the volcano) or θ > 170° (where the wave refocuses at the
     antipode) are excluded.
   - The local crossing angle φ comes from the GNV1B ground track.
3. **Templates.** For each crossing:
   - The template is the forward model above at the local φ. Its amplitude follows the spreading
     shape `√(1/sin θ)`, normalised to 1 at θ = 90°, and it uses the atmosphere only.
   - Pulse widths are σ = 2.5, 5 and 10 min. The 1.25 min width is dropped because its SNR per
     crossing is too low.
   - `T_true` is the sum over all crossings in the window.
4. **AOD1B's copy of the pulse.**
   - ECMWF's analyses may hold the pulse. If so, AOD1B holds it at the 3-hourly epochs, and the
     linear interpolation turns it into two faded copies at the wrong places.
   - `T_aod` is the same pulse placed at the front positions of the 3-hourly epochs and weighted
     by the interpolation weights.
   - The residual is fitted as `r = A·T_true − B·T_aod + noise` by generalised least squares.
   - The noise spectrum is estimated (Welch, 6 h segments) from the residuals before the
     eruption only.
   - `B = 0` means AOD1B did not carry the pulse, and `B = A` means it carried it fully.
5. **Statistic.**
   - `Z = Â/σ_Â`, the largest of the three widths.
   - The window is the first 24 h after the eruption: 15 January 04:14:45 to 16 January 04:14:45
     UTC.
   - Â is in units of the template at 90°, i.e. hPa.
6. **Null.** The same procedure, widths and window are run with fake eruptions:
   - The source is the same volcano, with origin times every 10 minutes from 13 January 00:00 to
     14 January 04:14 UTC. That gives 170 fakes whose whole 24 h window lies before the real
     eruption.
   - To vary the geometry, 30 random fake sources over the oceans are added, with origin times
     drawn from the same span.
7. **Decision.**
   - **Detection:** `Z` exceeds the 99th percentile of the null `Z` values, and Â > 0 (a positive
     pressure pulse).
   - The p-value is the rank of `Z` among the null values.

## Reported either way

- **Â and B̂ with 68% and 95% intervals**, in hPa. They are compared with the 1.3–2 hPa recorded
  far away, times the ocean factor 1.63–1.94.
- **If no detection:** the 95% upper limit on A.
  - An upper limit below 1.3 hPa would contradict the forward model.
  - B̂ then says whether AOD1B removed the pulse.
- **Secondary results, not part of the decision:**
  - per-crossing Z for crossings at φ > 45°;
  - the second passage (24–44 h);
  - KBR on 17 January.

## Known limits

- **Flat-Earth template:**
  - At 490 km, horizontal distances are about 7% larger than on the ground.
  - The line of sight tilts about 1° from the track.
  - This scales A by a few percent and does not change Z.
- **GRACE-D accelerometer data:** these are transplanted from GRACE-C, which raises the
  non-gravitational noise above the nominal SuperSTAR curve.
- **Pulse shape:** the real pulse is not a Gaussian. It has a positive jump followed by a trough.
  A shape mismatch lowers Z; it does not create a false detection, because the null uses the same
  templates.
- **Background model:** AOD1B RL06, not RL07. RL06's operational ECMWF analyses and RL07's ERA5
  differ in how much of the pulse they resolve; `B` absorbs this.
