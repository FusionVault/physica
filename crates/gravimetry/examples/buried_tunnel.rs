//! Finding a void from the surface: the gravity and gradient of a buried tunnel, and how long a
//! quantum gravimeter must average to see it. The geometry is the Birmingham survey's
//! ([Stray et al. 2022](https://doi.org/10.1038/s41586-021-04315-3)), which found a tunnel at
//! about 1.9 m depth.
//!
//! ```text
//! cargo run -p gravimetry --example buried_tunnel
//! ```

use gravimetry::local::{EOTVOS, MICROGAL, averaging_time, cylinder};

fn main() {
    let (radius, depth, contrast) = (1.0, 1.9, -2000.0); // an air-filled tunnel in soil
    println!("{:>8} {:>12} {:>12}", "offset/m", "Δg/µGal", "Γuu/E");
    for offset in [0.0, 1.0, 2.0, 3.0, 5.0] {
        let (gz, gradient) = cylinder(radius, depth, contrast, offset);
        println!("{offset:>8.1} {:>12.1} {:>12.1}", gz / MICROGAL, gradient / EOTVOS);
    }
    let (peak, _) = cylinder(radius, depth, contrast, 0.0);
    println!("peak anomaly {:.0} µGal", peak / MICROGAL);

    // An Exail-class gravimeter at 500 nm/s²/√Hz needs this long per station for a 10:1 signal.
    let target = peak.abs() / 10.0;
    println!("to resolve it 10:1 ({:.0} nm/s²): {:.0} s per station", target * 1e9, averaging_time(500e-9, target));
    assert!((peak / MICROGAL + 44.1).abs() < 0.1);
    // The gradient changes sign away from the tunnel, which is what localises it.
    assert!(cylinder(radius, depth, contrast, 0.0).1 < 0.0 && cylinder(radius, depth, contrast, 3.0).1 > 0.0);
}

#[test]
fn the_tunnel_anomaly() {
    main();
}
