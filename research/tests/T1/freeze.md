# T1: the pipeline as frozen

**Frozen:** 2026-10-24, as the git tag `prereg/T1-frozen`, before any residual after the eruption
was computed. Everything below was developed and checked on 13 January 00:00 to 15 January
04:15 GPS only; the code drops every later record unless run with `--unblind`.

The registration is [`T1-2026-10-24.md`](T1-2026-10-24.md). The code is `gravviz t1-noise` and
`gravviz t1-test` ([`apps/gravviz/src/t1.rs`](../apps/gravviz/src/t1.rs)). To rerun:

```bash
cargo run --release -p gravviz -- t1-test data/gracefo_l1b data/gracefo_l1b/aod1b data/GOCO06s.gfc data/costg/GSM-2_2022001-2022031_GRFO_COSTG_BF01_0201.gfc data/tides/GOT4.10c/grids_oceantide figures --unblind 2022-01-13 2022-01-14 2022-01-15 2022-01-16
```

## Inputs

| Input | File | SHA-256 |
|---|---|---|
| GRACE-FO Level-1B RL04, 13–17 January 2022 | `gracefo_1B_2022-01-1?_RL04.ascii.{LRI,noLRI}.tgz` (PO.DAAC) | in `data/gracefo_l1b/SHA256SUMS`, copied below |
| AOD1B RL06, January 2022 | `AOD1B_2022-01_06.tgz` (PO.DAAC) | as above |
| GOT4.10c ocean tides | `got4.10c.tar.gz` (NASA GSFC) | `2431decebe847414925289e4904bce09709bf668ce3c888700cc884a1236e038` |
| COST-G January 2022, GOCO06s | already in the project | — |

## Deviations from the registration, and why

1. **Ocean-tide model: GOT4.10c** (Ray 2013; eight constituents, spherical harmonics to degree
   120 with load Love numbers, nodal corrections). The registration named FES2014b, or EOT20 if
   FES2014b was not openly available. FES2014b's coefficients need a login; EOT20 is distributed
   only as netCDF-4/HDF5, which this project cannot read without a heavy new dependency. GOT4.10c
   is openly downloadable in plain text. Chosen on availability. Its phases were checked on the
   pre-eruption residual (regressing it on each species' predicted signal and that signal shifted
   by 90°): diurnal 0.94 in phase and 0.02 in quadrature, semidiurnal 1.04 and −0.04. The check
   caught a 180° sign error in the diurnal arguments, fixed before this freeze.
2. **Orbits: dynamic orbits fitted to GNV1B.** The registration's model used GNV1B positions and
   velocities directly. Their velocities put noise into the velocity term (the 0.25 m/s²
   centripetal term must cancel to 10⁻¹⁰), so each satellite's orbit is integrated through every
   revolution (RK4, 2 s steps, the full gravity model, rotating-frame terms, nine empirical
   accelerations) and fitted to GNV1B positions every 10 s. Median fit 6.6 cm. The band-passed
   residual fell from 4.8 to 3.4×10⁻⁹ m/s² before tides.
3. **Moving-front templates.** The forward model assumed the front frozen while the pair crosses
   it. Each template now uses the pair's ground speed across the front relative to the front's own
   318.2 m/s, which reverses past the antipode (at 10° crossings the frozen front was off by about
   26%). AOD1B's 3-hourly copies do not move.
4. **Noise spectrum:** Welch with 10,800-sample (6 h) Hann segments zero-padded to 16,384, as
   registered; the band for the inner products is 0.7–14 mHz, where the band-pass passes more
   than about 1% of the power.
5. **Monthly field to degree 90, as registered.** Truncating it at 40 or 20 (GOCO06s above) gave
   no improvement on the pre-eruption residual (2.26 and 2.41 against 2.25×10⁻⁹ m/s²).

## The frozen null, before the eruption

Computed on 13 January 00:00 to 15 January 04:15 GPS with this pipeline (`figures/t1_null.txt`;
the `--unblind` run recomputes it with the tagged code):

- **Null:** 200 fake eruptions (170 from Tonga every 10 minutes from 13 January 00:00, 30 from
  random ocean points). Null Z: median 0.18, 1st–99th percentile −2.06 to 2.34. The **detection
  threshold is Z > 2.34**; the null is now close to a unit Gaussian (2.33).
- **Sensitivity:** the amplitude's σ is 0.77 hPa in a median window (5–95%: 0.43–1.65 hPa), with
  15 crossings. A 1.5 hPa pulse amplified 1.6–1.9× by the ocean would give Z ≈ 3.
- **Noise:** band-passed residual 2.2×10⁻⁹ m/s² with ocean tides and GNV1B orbits (it was 4.8 before
  tides and fitted orbits).

## Input checksums (SHA-256)

```
304e516f6aafef10cf6511d4f2db7e4d9d4a6e79900040fc90639eb82dd0374c  AOD1B_2022-01_06.tgz
60ad6694ad1cbbb4414153e59fd0077df466ef78c56227147ab68995a7b25160  gracefo_1B_2022-01-13_RL04.ascii.LRI.tgz
e84128587e1cb4d7aba40361f6109e720d0ab3ee15b7eac93c36c7b8e3476c6d  gracefo_1B_2022-01-13_RL04.ascii.noLRI.tgz
c6123e6ff3afb21ed8dd623c1b5d8954a08ddc6cb86f716ef2cc064508d0f17c  gracefo_1B_2022-01-14_RL04.ascii.LRI.tgz
04610c10070fb7cd56f3557aa6192ef39ecb5f042562664c2fe36cdb7ae67330  gracefo_1B_2022-01-14_RL04.ascii.noLRI.tgz
e9139b2c4cbad95597f73ec5358c8f133c3f951dae371eb62573372efedfe332  gracefo_1B_2022-01-15_RL04.ascii.LRI.tgz
b4ad3ea786683e635ba8a3f77813b2d83c92d0e40cf0262736c59068c675931a  gracefo_1B_2022-01-15_RL04.ascii.noLRI.tgz
f41a55d72c838021a93852de4434f435643e0e98dece4e742bd72a858bfddbc9  gracefo_1B_2022-01-16_RL04.ascii.LRI.tgz
b877c7e9eab6585f6c8b4e32fde0b95b684b12491e9b6359ae9d1b52c51c2222  gracefo_1B_2022-01-16_RL04.ascii.noLRI.tgz
84062514956453178e6421f146170f61a8f4992b513a2f08fc8ed7a6eb226ad3  gracefo_1B_2022-01-17_RL04.ascii.LRI.tgz
beef244f1613a47429bbadd75c5aa8ead92314a0f1092f5c3c62310d8a7ecce2  gracefo_1B_2022-01-17_RL04.ascii.noLRI.tgz
2431decebe847414925289e4904bce09709bf668ce3c888700cc884a1236e038  got4.10c.tar.gz
```
