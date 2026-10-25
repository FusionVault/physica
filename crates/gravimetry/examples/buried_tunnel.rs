//! Finding a void from the surface: the gravity and gradient of a buried tunnel, and how long a
//! quantum gravimeter or gradiometer must average to see it. The geometry and the gradiometer are
//! the Birmingham survey's ([Stray et al. 2022](https://doi.org/10.1038/s41586-021-04315-3)),
//! which found a tunnel at about 1.9 m depth.
//!
//! ```text
//! cargo run -p gravimetry --example buried_tunnel
//! ```

use gravimetry::atom::AtomInterferometer;
use gravimetry::local::{EOTVOS, MICROGAL, averaging_time, cylinder};

fn main() {
    let (radius, depth, contrast) = (1.0, 1.9, -2000.0); // an air-filled tunnel in soil

    // Stray et al.'s gradiometer: two single-diffraction interferometers 1 m apart vertically,
    // T = 85 ms, 10⁵ atoms at 0.7 Hz. It reads the difference in gravity between its sensors.
    let stray = AtomInterferometer { t: 0.085, contrast: 1.0, atoms: 1e5, cycle: 1.0 / 0.7, diffraction: 1 };
    let baseline = 1.0;
    let difference = |offset: f64| (cylinder(radius, depth, contrast, offset).0 - cylinder(radius, depth + baseline, contrast, offset).0) / baseline;

    println!("{:>8} {:>10} {:>10} {:>12} {:>12}", "offset/m", "Δg/µGal", "Γuu/E", "Δg/1 m in E", "phase/mrad");
    for offset in [0.0, 1.0, 2.0, 3.0, 5.0] {
        let (gz, gradient) = cylinder(radius, depth, contrast, offset);
        let d = difference(offset);
        println!("{offset:>8.1} {:>10.1} {:>10.1} {:>12.1} {:>12.2}", gz / MICROGAL, gradient / EOTVOS, d / EOTVOS, stray.phase(d * baseline) * 1e3);
    }
    let (peak, _) = cylinder(radius, depth, contrast, 0.0);
    println!("peak anomaly {:.0} µGal", peak / MICROGAL);

    // An Exail-class gravimeter at 500 nm/s²/√Hz needs this long per station for a 10:1 signal.
    let target = peak.abs() / 10.0;
    println!("gravimeter, to resolve it 10:1 ({:.0} nm/s²): {:.0} s per station", target * 1e9, averaging_time(500e-9, target));

    // The gradiometer at its projection limit, and at the 466 E/√Hz it reached outdoors.
    let target = difference(0.0).abs() / 10.0;
    for (what, asd) in [("projection limit", stray.gradient_asd(baseline)), ("outdoors", 466.0 * EOTVOS)] {
        println!("gradiometer, {what} {:.0} E/√Hz: {:.0} s per station", asd / EOTVOS, averaging_time(asd, target));
    }

    assert!((peak / MICROGAL + 44.1).abs() < 0.1);
    // The gradient changes sign away from the tunnel, which is what localises it.
    assert!(cylinder(radius, depth, contrast, 0.0).1 < 0.0 && cylinder(radius, depth, contrast, 3.0).1 > 0.0);
    // Over a line mass the difference across the baseline is the point gradient at the geometric
    // mean of the two depths, 2Gλ/(z₁z₂): about Stray et al.'s 150 E and 17.5 mrad.
    let mid = cylinder(radius, (depth * (depth + baseline)).sqrt(), contrast, 0.0).1;
    assert!((difference(0.0) / mid - 1.0).abs() < 1e-12);
    assert!((difference(0.0) / EOTVOS + 152.2).abs() < 0.1 && (stray.phase(difference(0.0)) * 1e3 + 17.7).abs() < 0.1);
}

#[test]
fn the_tunnel_anomaly() {
    main();
}
