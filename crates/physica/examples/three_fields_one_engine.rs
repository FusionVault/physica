//! Three fields of the Earth through one interface: gravity, the main magnetic field and the
//! gravitomagnetic field of its spin, evaluated along a Swarm-like orbit by the same loop. The
//! magnetic field is checked against NOAA's WMM2025 test value at 100 km (Table 6 of the
//! WMM2025 report).
//!
//! ```text
//! cargo run -p physica --example three_fields_one_engine
//! ```

use geopotential::{Coefficients, Field, geocentric};
use gravitics::gem::GravitomagneticDipole;
use magnetics::main_field::wmm2025;

fn main() {
    let mut gravity = Coefficients::zeros(3.986_004_415e14, 6_378_136.3, 2);
    gravity.set(0, 0, 1.0, 0.0);
    gravity.set(2, 0, -4.841_694_95e-4, 0.0); // GOCO06s C̄20
    let magnetism = wmm2025().at(2025.0);
    let spin = GravitomagneticDipole { j: 5.86e33 };
    let fields: [(&str, &str, &dyn Field); 3] = [("gravity", "m/s²", &gravity), ("main field", "nT", &magnetism), ("gravitomagnetic", "m/s²", &spin)];

    for lat in [-60.0f64, 0.0, 60.0] {
        let (r, lat_c) = geocentric(lat.to_radians(), 450e3);
        let values: Vec<String> = fields.iter().map(|(name, unit, f)| format!("{name} {:.4e} {unit}", f.magnitude(r, lat_c, 0.0))).collect();
        println!("{lat:+5.0}° at 450 km: {}", values.join(", "));
    }

    // WMM2025 Table 6: 2025.0, 100 km, 80°N 0°E, F = 52 964.9 nT. Magnitude is frame-free, so the
    // geocentric vector gives it directly.
    let (r, lat_c) = geocentric(80f64.to_radians(), 100e3);
    let f = magnetism.magnitude(r, lat_c, 0.0);
    println!("main field at 80°N, 100 km: {f:.1} nT (NOAA: 52964.9 nT)");
    assert!((f - 52_964.9).abs() < 0.05);
}

#[test]
fn wmm_table_6_through_the_trait() {
    main();
}
