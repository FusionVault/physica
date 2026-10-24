//! The Hulse–Taylor binary pulsar's orbit shrinks as it radiates gravitational waves. From its
//! published orbital elements, the Peters–Mathews formula predicts the period's decay, which 35
//! years of timing measure ([Weisberg & Huang 2016](https://arxiv.org/abs/1606.02744)).
//!
//! ```text
//! cargo run -p gravitational --example binary_pulsar
//! ```

use gravitational::energy::{peters_enhancement, peters_pbdot};
use metrology::SOLAR_MASS;

fn main() {
    // PSR B1913+16: masses (M☉), orbital period (days) and eccentricity.
    let (mp, mc, pb_days, e) = (1.438, 1.390, 0.322997448918, 0.6171340);
    let pbdot = peters_pbdot(mp * SOLAR_MASS, mc * SOLAR_MASS, pb_days * 86_400.0, e);
    println!("eccentricity boosts the radiated power {:.2}×", peters_enhancement(e));
    println!("predicted dPb/dt = {pbdot:.4e} s/s (general relativity, published: -2.40263e-12)");
    println!("measured, after the Galactic correction: -2.398e-12 ± 0.004e-12");
    println!("the orbit's period shortens by {:.1} µs per year", -pbdot * 365.25 * 86_400.0 * 1e6);
    // The masses are quoted to 10⁻³, which limits the agreement to about 2×10⁻⁴.
    assert!((pbdot / -2.40263e-12 - 1.0).abs() < 1e-3);
}

#[test]
fn reproduces_the_published_decay() {
    main();
}
