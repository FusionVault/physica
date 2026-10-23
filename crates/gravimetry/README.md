# gravimetry

GRACE-style processing of monthly gravity fields into water-storage anomalies, on top of
[`geopotential`](https://crates.io/crates/geopotential).

- **`love_numbers`:** elastic load Love numbers (Wahr et al. 1998, Table 1), interpolated.
- **`gaussian_weights`, `smooth`:** Gaussian averaging (Jekeli 1981).
- **`destripe`:** Swenson & Wahr (2006) north–south stripe removal, ported from gravity-toolkit.
- **`mean`, `anomaly`:** baseline means and anomalies, rescaling between `GM`/radius conventions.

The usual chain follows GRACE Tellus (Landerer & Swenson 2012):
1. Replace C̄20/C̄30 and add degree 1, using `gravity-field`'s TN-14/TN-13 readers.
2. Subtract a baseline mean.
3. Truncate to degree 60.
4. Destripe.
5. Apply a 300 km Gaussian.
6. Synthesize `Functional::WaterHeight`.

```rust
use geopotential::{Coefficients, Functional};

let (gm, r) = (3.986004415e14, 6378136.3);
let mut month = Coefficients::zeros(gm, r, 60);
month.set(2, 0, 1e-10, 0.0);
let baseline = gravimetry::mean(&[Coefficients::zeros(gm, r, 60)]).unwrap();
let w = gravimetry::gaussian_weights(300e3, r, 60);
let dewh = gravimetry::smooth(&gravimetry::anomaly(&month, &baseline), &w);
let love = gravimetry::love_numbers(60).unwrap();
let pole = dewh.grid(&Functional::WaterHeight { love: &love }, r, &[1.5707963], &[0.0])[0];
assert!(pole > 0.01 && pole < 0.02); // metres of water
```

**Not included yet:**
- glacial isostatic adjustment;
- DDK filters;
- basin averages and leakage corrections.

Licensed under either of MIT or Apache-2.0 at your option. `destripe` is ported from
gravity-toolkit (MIT, Copyright (c) 2020 Tyler C Sutterley).
