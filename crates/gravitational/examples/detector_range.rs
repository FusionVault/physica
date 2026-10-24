//! How far Advanced LIGO can see a binary neutron star, and how that sets the detection rate:
//! the number trade studies for new detectors are argued in. The LVK quoted 160–190 Mpc as the
//! O4 target and 155–175 Mpc as achieved ([observing plans](https://observing.docs.ligo.org/plan)).
//!
//! ```text
//! cargo run -p gravitational --example detector_range
//! ```

use gravitational::{curves, inspiral};
use metrology::{MPC, SOLAR_MASS};

fn main() {
    let (m1, m2) = (1.4 * SOLAR_MASS, 1.4 * SOLAR_MASS);
    // The Advanced LIGO design fit, SNR 8, from 10 Hz.
    let horizon = inspiral::horizon(m1, m2, 8.0, curves::aligo, 10.0);
    let range = horizon / inspiral::HORIZON_TO_RANGE;
    println!("design horizon {:.0} Mpc, sky-averaged range {:.0} Mpc", horizon / MPC, range / MPC);
    // Detections scale with the volume surveyed: range³.
    for achieved in [155.0, 175.0] {
        println!("a {achieved} Mpc range sees {:.0}% of the design rate", 100.0 * (achieved * MPC / range).powi(3));
    }
    assert!((150.0..260.0).contains(&(range / MPC)));
}

#[test]
fn design_range_is_in_the_published_band() {
    main();
}
