//! Earth's dipole: where its axis points and how fast it is weakening. The WMM2025 coefficients
//! put the north geomagnetic pole at 80.79°N, 72.76°W, and the dipole moment falls about 5% per
//! century.
//!
//! ```text
//! cargo run -p magnetics --example magnetic_pole
//! ```

use magnetics::main_field::wmm2025;

fn main() {
    let wmm = wmm2025();
    for year in [2025.0, 2027.5, 2030.0] {
        let (lat, lon) = wmm.geomagnetic_pole(year);
        println!("{year}: north geomagnetic pole {:.2}°N {:.2}°W, dipole {:.4e} A m²", lat.to_degrees(), -lon.to_degrees(), wmm.dipole_moment(year));
    }
    let decline = 1.0 - wmm.dipole_moment(2030.0) / wmm.dipole_moment(2025.0);
    println!("the dipole weakens {:.1}% per century", decline * 20.0 * 100.0);
    let (lat, lon) = wmm.geomagnetic_pole(2025.0);
    assert!((lat.to_degrees() - 80.79).abs() < 0.005 && (lon.to_degrees() + 72.76).abs() < 0.005);
    assert!((4.0..6.0).contains(&(decline * 2000.0)));
}

#[test]
fn pole_and_decline() {
    main();
}
