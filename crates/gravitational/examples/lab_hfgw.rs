//! Could a laboratory make a gravitational wave anyone could detect? The quadrupole formula says
//! what an energetic, asymmetric implosion radiates, and the high-frequency detector registry
//! says what the best instruments can hear. The gap is the answer.
//!
//! ```text
//! cargo run -p gravitational --example lab_hfgw
//! ```

use gravitational::hf::{REGISTRY, gertsenshtein_probability};
use gravitational::laboratory::{burst_energy, quadrupole_strain};
use gravitational::quanta::graviton_count;

fn main() {
    // An inertial-fusion implosion: ~20 kJ of shell kinetic energy stagnating in 100 ps.
    let (energy, tau) = (2e4, 1e-10);
    println!("{:>10} {:>14} {:>14}", "asymmetry", "h at 1 m", "energy/J");
    for eps in [0.001, 0.01, 0.1] {
        let i_ddot = 2.0 * eps * energy;
        println!("{eps:>10.3} {:>14.2e} {:>14.2e}", quadrupole_strain(i_ddot, 1.0), burst_energy(i_ddot, tau));
    }
    let h = quadrupole_strain(2.0 * 0.01 * energy, 1.0);

    // The best published limits near those frequencies.
    let best = REGISTRY.iter().filter(|l| l.band.0 > 1e9 && l.quantity.starts_with("h₀")).min_by(|a, b| a.value.total_cmp(&b.value));
    if let Some(l) = best {
        println!("best limit near GHz: {} at h₀ = {:.1e}", l.name, l.value);
        println!("the laboratory strain is {:.0e} times smaller", h / l.value);
    }

    // The other route: converting a laser pulse into gravitons in a strong static field.
    let p = gertsenshtein_probability(20.0, 2.0);
    let photons = 4.0 / (metrology::H * metrology::C / 0.8e-6); // a 4 J pulse at 0.8 µm
    println!("Gertsenshtein conversion in 20 T over 2 m: P = {p:.2e}");
    println!("  a 4 J pulse ({photons:.1e} photons) yields {:.1e} gravitons", photons * p);
    println!("  for comparison, a detectable wave holds {:.1e} gravitons", graviton_count(1e-21, 150.0));
    assert!(h / 3.9e-21 < 1e-20); // about twenty orders short
}

#[test]
fn the_gap_is_twenty_orders() {
    main();
}
