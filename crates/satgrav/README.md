# satgrav

Instrument and mission models for satellite gravimetry, from electrostatic to cold-atom (quantum)
sensors.

- **Noise models:** `Accelerometer` and `Ranging` amplitude spectral densities (Zingerle et al.
  2024, eqs. 1–8), plus `pair_noise` for the line-of-sight noise of a GRACE-type pair.
- **`AtomInterferometer`:** light-pulse atom interferometers in microgravity. It covers phase,
  shot noise, and accelerometer and gradiometer noise. These are checked against the published
  CARIOQA, GRICE and Trimeche et al. (2019) figures.
- **`gradiometry_errors`:** closed-form formal coefficient errors of a gradiometer mission,
  checked against brute-force least squares.
- **`noise`:** an FFT, noise coloured to any spectrum, and whitening.
- **`sim`:** a closed-loop simulation of a GRACE-type pair. Two satellites fly circular orbits
  over a rotating Earth. The observable is their line-of-sight gravity difference, estimated by
  whitened least squares, and the output is the actual errors with a χ² consistency check.

```rust
use satgrav::{Accelerometer, Ranging, pair_noise};

// GRACE-FO class at 1 mHz: SuperSTAR accelerometers plus the GRACE-FO laser interferometer.
let n = pair_noise(Accelerometer::SuperStar, Ranging::GraceFo, 1e-3);
assert!((n - 3.5e-10).abs() < 0.05e-10);
```

**Not modelled yet:**
- time-variable gravity, and so temporal aliasing;
- orbit dynamics;
- attitude errors;
- multi-axis gradiometers and polar gaps.

`hybrid_asd` is an inverse-variance stand-in for an electrostatic plus cold-atom hybrid.

Licensed under either of MIT or Apache-2.0 at your option.
