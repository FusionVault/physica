//! What the ice sheets' mass loss, weighed from orbit, means for sea level: the JPL mascon trend
//! for Greenland ([NASA SVS](https://svs.gsfc.nasa.gov/webapps/sealevelgraphs/Updated/data/greenland_mass_200204_202404.txt))
//! and IMBIE's 30-year total ([Otosaka et al. 2023](https://repository.library.noaa.gov/view/noaa/69781)).
//!
//! ```text
//! cargo run -p gravimetry --example ice_and_sea_level
//! ```

use geopotential::RHO_WATER;
use gravimetry::sea_level_equivalent;

const GT: f64 = 1e12; // kg

fn main() {
    // Greenland, April 2002 to April 2024: −267.49 ± 21 Gt/yr.
    let rate = sea_level_equivalent(267.49 * GT);
    println!("Greenland adds {:.2} ± {:.2} mm/yr to global sea level", rate * 1e3, sea_level_equivalent(21.0 * GT) * 1e3);
    // Spread over the ice sheet's 1.71 million km², the same loss is a water-height change of:
    println!("…a {:.0} mm/yr fall in equivalent water height over the ice sheet", 267.49 * GT / (RHO_WATER * 1.71e12) * 1e3);
    assert!((rate * 1e3 - 0.739).abs() < 0.001);

    // IMBIE: both ice sheets raised sea level 21.0 ± 1.9 mm over 1992–2020.
    let mass = 21.0e-3 * RHO_WATER * gravimetry::OCEAN_AREA / GT;
    println!("IMBIE's 21.0 mm is {mass:.0} Gt of ice");
    assert!((mass - 7598.0).abs() < 1.0);
}

#[test]
fn sea_level_from_mass() {
    main();
}
