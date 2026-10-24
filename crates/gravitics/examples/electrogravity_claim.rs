//! A pulsed high-voltage "gravity" claim, computed three ways: the classical force that acts in
//! the same apparatus, what a real gravitational field would do, and the experiment that tells
//! them apart.
//!
//! ```text
//! cargo run -p gravitics --example electrogravity_claim
//! ```

use gravitics::electrogravity::{acceleration, dielectrophoretic_force, field_angular_momentum, mass_independence, rotation_torque};

fn main() {
    // Setup 1: a 2 cm square of 0.15 mm plastic between electrodes at 40 kV, where the field
    // reaches ~10⁷ V/m and falls off over centimetres.
    let volume = 0.02 * 0.02 * 0.15e-3;
    let f = dielectrophoretic_force(volume, 2.2, 2e16);
    println!("dielectrophoretic force on the sheet: {:.2} mN, the weight of {:.0} mg", f * 1e3, f / metrology::G0 * 1e6);
    println!("  it follows ∇E², so it pulls toward whichever electrode is live, whatever its sign");
    println!("  acceleration of the sheet: {:.1} m/s²", acceleration(f, volume, 950.0));

    // Setup 2: the ball between open coils. Their 10 µA makes only 7×10⁻¹¹ T, so the angular
    // momentum stored in the crossed fields is negligible; the polarisation torque is not.
    let l = field_angular_momentum(2.5e5, 7.4e-11, 0.05, 1.6e-3);
    let torque = rotation_torque(4.447e-3, 2.3, 2.5e5);
    println!("field angular momentum in the tank: {l:.1e} kg m²/s, too little to turn the ball");
    println!("polarisation torque on the ball, at most {torque:.1e} N m, enough to wind a cotton thread hundreds of turns");

    // The discriminator: same shape, different density.
    let ratio = mass_independence(950.0, 19_300.0);
    println!("polyethylene against tungsten of the same shape:");
    println!("  a classical force gives the lighter body {ratio:.1}× the acceleration");
    println!("  a gravitational field gives both the same, so the ratio would be 1.0");
    assert!((ratio - 20.3).abs() < 0.1);
    assert!(f > 1e-3); // milli-newtons, far above any claimed gravitational effect
}

#[test]
fn the_classical_force_dominates() {
    main();
}
