//! Designing an Earth-observation orbit with the gravity field: Earth's flattening (J2) turns an
//! orbit's plane, and at the right inclination it turns once a year, keeping the satellite over
//! every place at the same local time. Landsat 8 flies at 705 km and 98.2°
//! ([USGS](https://www.usgs.gov/landsat-missions/landsat-8)).
//!
//! ```text
//! cargo run -p geopotential --example sun_synchronous_orbit
//! ```

use core::f64::consts::TAU;

use geopotential::read_gfc;

fn main() {
    // GOCO06s to degree 4 (the full model is from ICGEM).
    let field = read_gfc(include_str!("../tests/data/goco06s_n4.gfc")).expect("bundled model").model.base;
    let j2 = -5f64.sqrt() * field.c(2, 0);
    println!("J2 = {j2:.6e} (GOCO06s)");
    // Nodal precession of a circular orbit: dΩ/dt = −(3/2) n J2 (R/a)² cos i.
    let year_rate = TAU / (365.2422 * 86_400.0);
    let inclination = |altitude: f64| {
        let a = 6_378_137.0 + altitude;
        let n = (field.gm / a.powi(3)).sqrt();
        (-year_rate / (1.5 * n * j2 * (field.radius / a).powi(2))).acos()
    };
    for (name, altitude) in [("Landsat 8", 705e3), ("GRACE-FO height", 490e3)] {
        println!("{name:<16} at {:.0} km: sun-synchronous at {:.2}°", altitude / 1e3, inclination(altitude).to_degrees());
    }
    assert!((inclination(705e3).to_degrees() - 98.2).abs() < 0.05);
}

#[test]
fn landsat_inclination() {
    main();
}
