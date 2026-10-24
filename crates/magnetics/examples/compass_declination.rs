//! A compass correction from the World Magnetic Model, the navigation standard of NATO, the IHO
//! and every phone: magnetic declination and its drift at a few places, and why drillers need it
//! to 0.1° ([BGS](https://geomag.bgs.ac.uk/data_service/directionaldrilling/ifr.html)).
//!
//! ```text
//! cargo run -p magnetics --example compass_declination
//! ```

use magnetics::main_field::wmm2025;

fn main() {
    let wmm = wmm2025();
    for (place, lat, lon) in [("London", 51.5074, -0.1278), ("New York", 40.7128, -74.0060), ("Sydney", -33.8688, 151.2093), ("Singapore", 1.3521, 103.8198)] {
        let e = wmm.elements(2026.5, f64::to_radians(lat), f64::to_radians(lon), 0.0);
        println!("{place:<10} D {:+6.2}° ({:+.2}°/yr), I {:+6.2}°, F {:.0} nT", e.declination.to_degrees(), e.declination_rate.to_degrees(), e.inclination.to_degrees(), e.f);
    }
    // A wellbore steered 0.1° off over a 5 km horizontal reach ends up this far to the side:
    println!("0.1° of declination error over 5 km: {:.1} m at the bottom hole", 5000.0 * 0.1f64.to_radians().tan());

    // NOAA's published test point (WMM2025 report, Table 6): 2025.0, sea level, 80°N 0°E.
    let e = wmm.elements(2025.0, 80f64.to_radians(), 0.0, 0.0);
    assert!((e.declination.to_degrees() - 1.28).abs() < 0.005 && (e.f - 55178.5).abs() < 0.05);
}

#[test]
fn declination_at_the_official_test_point() {
    main();
}
