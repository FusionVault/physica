# plasmoid

Fusion and self-organising plasma physics, computed and tested against published numbers.

```rust
use plasmoid::fusion::{Reaction, bosch_hale};

// D–T at 10 keV, Bosch & Hale (1992): 1.136×10⁻²² m³/s.
let sv = bosch_hale(Reaction::DT, 10.0).unwrap();
assert!((sv / 1.1362e-22 - 1.0).abs() < 1e-3);
```

Part of [quantum-gravimeter](https://github.com/FusionVault/quantum-gravimeter).
