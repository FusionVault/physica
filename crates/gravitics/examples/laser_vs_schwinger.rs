//! Planning a strong-field experiment: how far the most intense laser ever fired is from the
//! field at which the vacuum breaks down into electron–positron pairs, and what an electron beam
//! does to close the gap (LUXE).
//!
//! ```text
//! cargo run -p gravitics --example laser_vs_schwinger
//! ```

use gravitics::strong_field::{E_SCHWINGER, Intensity, chi_e, peak_field, schwinger_exponent};

fn main() {
    // CoReLS, 2021: 1.1×10²³ W/cm² at 0.8 µm, linearly polarised.
    let record = peak_field(1.1e27, Intensity::Linear);
    println!("record laser field {record:.2e} V/m = 1/{:.0} of the Schwinger field {E_SCHWINGER:.3e} V/m", E_SCHWINGER / record);
    println!("pair production there is suppressed by e^−{:.0}", schwinger_exponent(record));
    // A 16.5 GeV electron meeting a₀ = 7.9 light head-on sees the critical field in its own frame.
    let chi = chi_e(16.5, 7.9, 1.55, 17.2f64.to_radians());
    println!("LUXE: an electron at 16.5 GeV reaches χ = {chi:.2} (≥ 1 is the strong-field regime)");
    assert!((E_SCHWINGER / record / 1454.0 - 1.0).abs() < 1e-3 && chi > 1.0);
}

#[test]
fn the_gap_is_three_orders() {
    main();
}
