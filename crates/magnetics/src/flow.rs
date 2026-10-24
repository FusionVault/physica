//! Conducting liquids in a magnetic field: the Hartmann flat profile and the pressure drop that
//! dominates liquid-metal blanket design.
//!
//! A field across a channel flattens the flow into a core with thin Hartmann layers of thickness
//! `a/Ha` at the walls, and the pressure needed to push the fluid grows with Ha (and with the
//! walls' conductance), which is why fusion blankets use insulating flow-channel inserts.
//!
//! ```
//! use magnetics::flow::pipe_gradient;
//!
//! // KIT's thin-walled pipe at wall conductance ratio c = 0.0727: 0.0678 σūB².
//! assert!((pipe_gradient(0.0727) - 0.0678).abs() < 1e-4);
//! ```

/// Hartmann velocity profile `u/ū` at position `eta` across the channel (−1 … 1) for Hartmann
/// number `ha`: `Ha(cosh Ha − cosh Ha·η)/(Ha cosh Ha − sinh Ha)`, evaluated without overflow.
pub fn hartmann_profile(ha: f64, eta: f64) -> f64 {
    // cosh(Ha·η)/cosh(Ha) = e^(Ha(|η|−1))(1 + e^(−2Ha|η|))/(1 + e^(−2Ha)).
    let e = eta.abs();
    let ratio = (ha * (e - 1.0)).exp() * (1.0 + (-2.0 * ha * e).exp()) / (1.0 + (-2.0 * ha).exp());
    ha * (1.0 - ratio) / (ha - ha.tanh())
}

/// Pressure gradient of Hartmann flow between plates, in units of `σūB²`: `tanh Ha/(Ha − tanh Ha)`
/// with insulating walls, `Ha/(Ha − tanh Ha)` with perfectly conducting ones.
pub fn hartmann_gradient(ha: f64, conducting_walls: bool) -> f64 {
    let t = ha.tanh();
    if conducting_walls { ha / (ha - t) } else { t / (ha - t) }
}

/// Inertialess pressure gradient of a thin-walled conducting pipe, in units of `σūB²`, for wall
/// conductance ratio `c = σ_w t/(σa)`: `c/(1 + c)`.
pub fn pipe_gradient(c: f64) -> f64 {
    c / (1.0 + c)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hartmann_limits() {
        assert!((hartmann_profile(10.0, 0.0) - 1.111_010_222).abs() < 1e-9);
        assert!((hartmann_profile(100.0, 0.0) - 1.010_101).abs() < 1e-6);
        assert!(hartmann_profile(10.0, 1.0).abs() < 1e-12);
        // Ha → 0 recovers Poiseuille: 3/2 at the centre, and a gradient of 3/Ha².
        assert!((hartmann_profile(1e-3, 0.0) - 1.5).abs() < 1e-5);
        assert!((hartmann_gradient(1e-2, false) * 1e-4 - 3.0).abs() < 1e-4);
        assert!((hartmann_gradient(10.0, false) - 0.111_111_11).abs() < 1e-7);
        assert!((hartmann_gradient(10.0, true) - 1.111_111_11).abs() < 1e-7);
        // Very large Ha, insulating: 1/(Ha − 1), and no overflow in the profile.
        assert!((hartmann_gradient(1e4, false) - 1.0 / 9999.0).abs() < 1e-12);
        assert!((hartmann_profile(1e4, 0.5) - 1e4 / 9999.0).abs() < 1e-9);
    }

    #[test]
    fn profile_has_unit_mean() {
        let n = 20_000;
        let mean: f64 = (0..n).map(|i| hartmann_profile(30.0, -1.0 + 2.0 * (i as f64 + 0.5) / n as f64)).sum::<f64>() / n as f64;
        assert!((mean - 1.0).abs() < 1e-6);
    }
}
