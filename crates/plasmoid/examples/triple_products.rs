//! How close each fusion machine has come to ignition: its triple product nTτ against the Lawson
//! requirement at its own temperature. The record values are Wurzel & Hsu's
//! ([Phys. Plasmas 32, 112106, 2025](https://arxiv.org/abs/2505.03834)), the ignition curve the
//! crate's Bosch–Hale D–T reactivity with bremsstrahlung.
//!
//! ```text
//! cargo run -p plasmoid --example triple_products
//! ```

use plasmoid::fusion::lawson;

fn main() {
    // Machine and shot, ion temperature (keV), density (m⁻³), energy confinement time (s), and the
    // published triple product (keV s m⁻³).
    let shots: [(&str, f64, f64, f64, f64); 6] = [
        ("TFTR 83546 (1995)", 43.0, 6.6e19, 0.28, 7.9e20),
        ("JT-60U E31872 (1998)", 16.8, 4.8e19, 0.69, 5.6e20),
        ("JET 99972 (2021)", 11.0, 8.5e19, 0.23, 2.2e20),
        ("W7-X (2017)", 3.5, 8.0e19, 0.22, 6.2e19),
        ("FuZE (2022)", 2.4, 4.9e23, 1e-6, 1.2e18),
        ("C-2W 122588 (2023)", 1.0, 1.69e19, 2e-3, 3.4e16),
    ];
    println!("{:<22} {:>12} {:>16}", "machine", "nTτ", "of ignition");
    for (name, t, n, tau, published) in shots {
        let ntt = n * t * tau;
        assert!((ntt / published - 1.0).abs() < 0.05, "{name}: {ntt:e} vs {published:e}");
        // Below 4.3 keV bremsstrahlung outruns D–T heating: no density and time ignite it.
        match lawson(t, f64::INFINITY) {
            Some(needed) => println!("{name:<22} {ntt:>12.2e} {:>15.1}%", 100.0 * ntt / (needed * t)),
            None => println!("{name:<22} {ntt:>12.2e} {:>16}", "too cold"),
        }
    }
    // The minimum of the ignition curve: about 3×10²¹ keV s m⁻³ near 14 keV.
    let best = (40..=600).map(|k| f64::from(k) * 0.05).filter_map(|t| lawson(t, f64::INFINITY).map(|n| n * t)).fold(f64::MAX, f64::min);
    println!("ignition needs at least {best:.2e} keV s m⁻³");
    assert!((best / 2.97e21 - 1.0).abs() < 0.01);
}

#[test]
fn triple_products_multiply_back() {
    main();
}
