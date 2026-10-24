//! The Hunga Tonga eruption's pressure wave as a mass signal for satellite gravimetry: the Lamb
//! wave ran at 318 m/s ([Wright et al. 2022](https://doi.org/10.1038/s41586-022-05012-5)), which is
//! the speed of sound at 252 K. Over the deep Pacific it outran the ocean's long waves, so the sea
//! rose under it instead of sinking, adding to the air's mass rather than cancelling it.
//!
//! ```text
//! cargo run -p gravimetry --example tonga_pressure_wave
//! ```

use gravimetry::atmosphere::{forced_ocean_response, pressure_to_mass, sound_speed};

fn main() {
    let lamb = sound_speed(252.0);
    println!("sound at 252 K: {lamb:.1} m/s (the observed Lamb speed was 314–318 m/s)");
    // A 1 hPa pulse over a 4 km deep ocean.
    let air = pressure_to_mass(100.0);
    let sea = forced_ocean_response(lamb, 4000.0); // sea-level change per inverse-barometer centimetre
    let ocean = sea * air; // kg/m² of seawater (1 cm of water = 10 kg/m²)
    println!("1 hPa of air is {air:.1} kg/m²; the ocean answers with {ocean:+.1} kg/m² (response {sea:+.3})");
    println!("total mass anomaly per hPa: {:.1} kg/m², {:.0}% more than the air alone", air + ocean, 100.0 * ocean / air);
    assert!((lamb - 318.2).abs() < 0.1 && (sea - 0.6328).abs() < 1e-3);
}

#[test]
fn the_ocean_amplifies_the_wave() {
    main();
}
