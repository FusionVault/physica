//! A tokamak design point checked against the empirical limits it has to respect: the edge safety
//! factor, the Greenwald density limit, the Troyon beta limit and the IPB98(y,2) confinement
//! scaling, for ITER and SPARC ([Creely et al. 2020](https://doi.org/10.1017/S0022377820001257)).
//!
//! ```text
//! cargo run -p plasmoid --example sparc
//! ```

use plasmoid::tokamak::{ITER, SPARC, Tokamak};

fn main() {
    for (name, t, beta_n) in [("ITER", ITER, 1.8), ("SPARC", SPARC, 1.0)] {
        println!("{name}: R = {} m, a = {} m, B = {} T, I = {} MA", t.major_radius, t.minor_radius, t.field, t.current);
        println!("  q* {:.2}, q95 {:.2} (kink limit is 1)", t.q_star(), t.q95());
        println!("  Greenwald limit {:.2e} m⁻³", t.greenwald());
        println!("  beta at β_N = {beta_n}: {:.2}% (Troyon ceiling at 2.8: {:.2}%)", t.troyon_beta(beta_n) * 100.0, t.troyon_beta(2.8) * 100.0);
    }
    // The confinement scaling, inverted: ITER's published 3.7 s needs this much net heating.
    let power = (10..2000).map(|k| f64::from(k) / 10.0).find(|&p| ITER.ipb98(10.1, p, 2.5, 1.7, 1.0) < 3.7);
    println!("ITER's 3.7 s requires about {:.0} MW of net heating", power.unwrap_or_default());
    assert!((ITER.q95() - 3.0).abs() < 0.01 && (SPARC.q_star() - 3.05).abs() < 0.01);
    assert!((power.unwrap_or_default() - 83.8).abs() < 0.5);
    let _ = Tokamak { ..ITER };
}

#[test]
fn design_points_match_the_published_limits() {
    main();
}
