//! A geomagnetic storm in numbers: the energy its ring current holds, the electric field it
//! induces in the ground, and how that compares with the benchmark grid operators plan for.
//!
//! ```text
//! cargo run -p magnetics --example storm
//! ```

use magnetics::main_field::wmm2025;
use magnetics::space_weather::{NERC_BENCHMARK, NERC_SUPPLEMENTAL, geoelectric_halfspace, ring_current_energy};

fn main() {
    let moment = wmm2025().dipole_moment(2025.0);
    println!("Earth's dipole moment in 2025: {moment:.3e} A m²");
    for (name, dst_nt) in [("quiet day", -20.0), ("May 2024 (Gannon)", -412.0), ("March 1989 (Québec)", -589.0), ("Carrington, low estimate", -850.0)] {
        println!("{name:<26} Dst {dst_nt:>7.0} nT → ring current {:.2e} J", ring_current_energy(dst_nt * 1e-9, moment));
    }

    // The induced ground field depends on how resistive the rock is.
    println!("\n{:>14} {:>12} {:>14}", "ground σ/(S/m)", "E/(V/km)", "skin depth/km");
    for sigma in [1e-2, 1e-3, 1e-4] {
        let (e, depth) = geoelectric_halfspace(1000e-9, 300.0, sigma);
        println!("{sigma:>14.0e} {:>12.2} {:>14.0}", e * 1e3, depth / 1e3);
    }
    let (e, _) = geoelectric_halfspace(1000e-9, 300.0, 1e-4);
    println!("\nNERC plans for {:.0} V/km (supplemental {:.0} V/km); resistive ground here gives {:.0} V/km",
        NERC_BENCHMARK * 1e3, NERC_SUPPLEMENTAL * 1e3, e * 1e3);
    assert!(e > NERC_BENCHMARK);
    assert!((ring_current_energy(-412e-9, moment) / 1.58e16 - 1.0).abs() < 0.02);
}

#[test]
fn storm_energies_and_ground_fields() {
    main();
}
