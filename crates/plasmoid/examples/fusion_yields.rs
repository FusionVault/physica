//! From a fusion shot's energy to the reactions and fuel behind it: NIF's 7 April 2025 record and
//! JET's final 69 MJ pulse ([LLNL](https://annual.llnl.gov/fy-2025/national-ignition-facility-2025);
//! [UK Government, 2024](https://www.gov.uk/government/news/jets-final-tritium-experiments-yield-new-fusion-energy-record)).
//!
//! ```text
//! cargo run -p plasmoid --example fusion_yields
//! ```

use metrology::{C, KEV, MEV};
use plasmoid::aneutronic::Q_DT;
use plasmoid::fusion::{M_D, M_T};

fn main() {
    let per_reaction = Q_DT * MEV; // 17.59 MeV, J
    // NIF N250406: 2.08 MJ of laser light in, 8.6 MJ of fusion out.
    let (laser, yield_nif) = (2.08e6, 8.6e6);
    let reactions = yield_nif / per_reaction;
    println!("NIF: target gain {:.2}, {reactions:.2e} D–T reactions (one neutron each)", yield_nif / laser);
    assert!((yield_nif / laser - 4.13).abs() < 0.01 && (reactions / 3.05e18 - 1.0).abs() < 0.01);

    // JET, 2023: 69.26 MJ in 5.2 s from 0.21 mg of fuel.
    let burned = 69.26e6 / per_reaction * (M_D + M_T) * KEV / (C * C); // kg of deuterium and tritium
    println!("JET: {:.1} MW for 5.2 s burned {:.3} mg of D–T (published 0.21 mg)", 69.26 / 5.2, burned * 1e6);
    assert!((burned * 1e6 - 0.205).abs() < 0.002);
}

#[test]
fn yields_match_the_published_shots() {
    main();
}
