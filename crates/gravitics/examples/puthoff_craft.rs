//! Puthoff's "engineered metric" craft (Puthoff 2010, JBIS 63, 82; arXiv:1204.2184), its Table 1
//! column taken at its word. Each narrative in its §4 sets how fast clocks run in a 10 m region,
//! √g₀₀ > 1. Here, for each one, is the static mass that metric needs in general relativity and in
//! the polarizable vacuum (it is negative in both), the field that would bend a gravimeter, and the
//! electromagnetic route the paper cites. The paper prints no numbers of its own; these are the
//! frontier-gaps research's, reproduced by the library.
//!
//! ```text
//! cargo run -p gravitics --example puthoff_craft
//! ```

use gravitics::spacetime::{Metric, engineered_mass, levi_civita_delay, redshift_factor};
use metrology::{C, MU_0, SOLAR_MASS};

/// The craft's radius in the research's examples, m.
const RADIUS: f64 = 10.0;

/// `g₀₀ − 1` for clocks running `rate = √g₀₀` times fast.
fn excess_for_rate(rate: f64) -> f64 {
    (rate - 1.0) * (rate + 1.0)
}

fn main() {
    // (narrative, g₀₀ − 1 inside), from the paper's §4 as the research reads it. The excess is
    // carried directly: √(1 + 10⁻¹⁸) is 1.0 in f64.
    let scenarios = [
        ("optical-clock detectable, g00 - 1 = 1e-18", 1e-18),
        ("1 ppm clock speed-up", excess_for_rate(1.0 + 1e-6)),
        ("plants age 10% faster", excess_for_rate(1.1)),
        ("'20 min inside vs a few (5) outside'", excess_for_rate(4.0)),
        ("'20 min vs a few (3)'", excess_for_rate(20.0 / 3.0)),
        ("300 K Wien peak 9.66 um seen at 550 nm", excess_for_rate(9.66 / 0.55)),
        ("20 kHz sound redshifted to 20 Hz", excess_for_rate(1000.0)),
    ];
    println!("Static mass that a 10 m region with clocks running sqrt(g00) times fast requires:");
    for (what, excess) in scenarios {
        let gr = engineered_mass(excess, RADIUS, Metric::GeneralRelativity);
        let pv = engineered_mass(excess, RADIUS, Metric::PolarizableVacuum);
        println!("  {what:42} GR {gr:+.3e} kg ({:.2e} Msun), PV {pv:+.3e} kg", gr / SOLAR_MASS);
    }
    // The cleanest detector: an unconfined region pulls outward with c²δ/(2R).
    let pull = C * C * 1e-18 / (2.0 * RADIUS);
    println!("A g00 - 1 = 1e-18 region pushes out with {:.2} mm/s^2 at 10 m; a 1 uGal gravimeter sees g00 - 1 = {:.1e}", pull * 1e3, 2.0 * 1e-8 * RADIUS / (C * C));
    // The route the paper cites, the Levi-Civita effect of a field's own energy, moves the metric the
    // wrong way and by almost nothing: K − 1 = 3cτ/(4L) for 1 T over 1 m.
    let u = 1.0 / (2.0 * MU_0);
    let k = 3.0 * C * levi_civita_delay(u, 1.0) / 4.0;
    println!("1 T over 1 m changes the vacuum's index by {k:.2e}, against the 1e-18 the mildest narrative needs");
    println!("Every narrative needs negative mass; no field known to the paper's references supplies any.");
    assert!(redshift_factor(0.25) == 0.5);
}

#[test]
fn every_narrative_needs_negative_mass() {
    main();
    // The research's §G table (puthoff_metric_checks.py), GR then PV, kg.
    for (excess, gr, pv) in [
        (1e-18, 6.733e9, 6.733e9),
        (excess_for_rate(1.0 + 1e-6), 1.347e22, 1.347e22),
        (excess_for_rate(1.1), 1.414e27, 1.412e27),
        (excess_for_rate(4.0), 1.010e29, 7.467e28),
        (excess_for_rate(20.0 / 3.0), 2.925e29, 1.703e29),
        (excess_for_rate(9.66 / 0.55), 2.070e30, 6.777e29),
        (excess_for_rate(1000.0), 6.733e33, 9.302e31),
    ] {
        let m = |metric| -engineered_mass(excess, RADIUS, metric);
        assert!((m(Metric::GeneralRelativity) / gr - 1.0).abs() < 1e-3, "GR at {excess}");
        assert!((m(Metric::PolarizableVacuum) / pv - 1.0).abs() < 1e-3, "PV at {excess}");
    }
    // 4.49 mm/s² outward at 10 m; 1 µGal resolves 2.2×10⁻²⁴; the field route gives 1.03×10⁻³⁸.
    assert!((C * C * 1e-18 / (2.0 * RADIUS) - 4.494e-3).abs() < 1e-6);
    assert!((2.0 * 1e-8 * RADIUS / (C * C) / 2.225e-24 - 1.0).abs() < 1e-3);
    assert!((3.0 * C * levi_civita_delay(1.0 / (2.0 * MU_0), 1.0) / 4.0 / 1.033e-38 - 1.0).abs() < 1e-3);
}
