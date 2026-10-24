# Formula audit, 2026-10-24

An independent audit re-derived every formula in `gwsens`, `satgrav::gw`, `gravitics` (then
`gravity-claims`), `gravimetry::atmosphere`, `gravviz tonga` and `geo5d`'s frame-dragging code from
the primary sources, and recomputed every test value in Python from the sources rather than from
the Rust. About 90 formulas and assertions were confirmed. These were not, and are now fixed:

## Errors

| Where | Problem | Fix |
|---|---|---|
| `satgrav::gw::pair_strain` | `h = x/x₀` (Ross et al. 2020's convention) is half the strain: a wave changes a length by `δL = ½hL`, the convention the rest of the module uses | `h = 2x/x₀`; 1 nm/√Hz over 220 km is 9.09×10⁻¹⁵ (satgrav 0.3.0) |
| `gwsens::hf::REGISTRY` | The Holometer's 2.4×10⁻¹⁹ Hz⁻¹ᐟ² is not in its paper, which reports better than 10⁻²¹/√Hz from 1 to 13 MHz | 10⁻²¹; MAGE now carries its measured 1.35×10⁻²⁰ rather than the predicted 6.6×10⁻²¹; the LSW and CAST bands follow Ejlli et al.; a test pins every entry |
| `gwsens::inspiral::horizon` | NaN when the binary merges below `f_lo` | returns 0, tested |
| `gravitics::superconductors` test | Tajmar's ρ*/ρ for niobium is 1.97×10⁻⁶ (λ_L = 39 nm), not 3.94×10⁻⁶; his field is twice that, `B_g = 2ωρ*/ρ = 3.9×10⁻⁶ω` | both asserted at 39 nm |

## Conventions and units

- `gwsens::lab::beam_tidal` now returns s⁻² (`−4GP/(c³ρ²)`); the paper's `R₀ₓ₀ₓ` is in m⁻².
- `gwsens::hf::cavity_noise` is now one-sided like the rest of the crate (twice the Living Review's
  two-sided Eq. 113).
- `gravitics::propulsion::mach_effect_thrust` includes the 4π the authors say SI needs; the test
  uses their 46 g stack (printed prefactor 1.87×10²¹, "~8 nN" without the 4π, about 0.1 µN with
  it). Fearn & Woodward reported "a couple of microNewtons" for this device, not 100 µN.
- `gravitics::hfgw`: the Li group's conversion rate is 8× the standard one; the shot-noise time is
  3.9×10¹⁷ s at their rate and 3.1×10¹⁸ s at the standard one. Both are computed from JASON's
  argument; neither number is printed in JASON.
- `satgrav::gw::clock_strain` returns Kolkowitz et al.'s σ_min, a fractional frequency, now
  documented as such.
- `geo5d`: the Gravity Probe B check now uses its orbit (a = 7027.4 km) and the published
  predictions, 6606.1 mas/yr and 39.2 mas/yr (the latter projected on the guide star) (geo5d 0.5.1).
- Doc notes added: the aLIGO fit's minimum (1.800×10⁻²⁴ at 228 Hz), the inclination average for
  sky-averaged curves, the field-of-view correction for magnet limits, NANOGrav's Ω computed from
  its A, and Pfenning & Ford's text against their Eq. 28.
- `gravviz tonga`: crossings now use the pair's speed relative to the moving front; at 10° the
  frozen-front template was off by about 26%.

## Tests tightened

The stochastic-background curve is pinned to √(2T); LISA's confusion noise, the wave flux and the
AION-10 strain are pinned to recomputed values; the ceiling test uses ΔN_eff = 0.2.
