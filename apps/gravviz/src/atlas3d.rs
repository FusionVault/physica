//! The data behind the 3D atlas: every number the page draws, computed by the crates.
//!
//! - geopotential: the GOCO06s geoid on a 2° grid;
//! - satgrav and the T1 pipeline: GRACE-FO's orbits for the 24 hours after the Tonga eruption;
//! - gravimetry: the Lamb wave's speed;
//! - geo5d: LAGEOS's frame-dragging node rate;
//! - gwsens and satgrav: detector noise curves and published limits;
//! - plasmoid: the spheromak eigenvalue, Bennett and Ostriker numbers;
//! - gravitics: warp-bubble energies, horizon temperatures and the Schwinger field.

use std::fmt::Write;
use std::path::Path;

use geopotential::Functional;
use gravity_field::read_gfc;

use crate::{Result, t1};

fn list(v: impl IntoIterator<Item = f64>, digits: usize) -> String {
    let items: Vec<String> = v.into_iter().map(|x| format!("{x:.digits$}")).collect();
    format!("[{}]", items.join(","))
}

fn sci(v: impl IntoIterator<Item = f64>) -> String {
    let items: Vec<String> = v.into_iter().map(|x| format!("{x:.4e}")).collect();
    format!("[{}]", items.join(","))
}

pub fn export(goco: &Path, l1b_dir: &Path, out: &Path) -> Result<()> {
    let mut j = String::from("{");

    // Geoid (GRS80 normal field removed), 2° cells, south to north, west to east.
    let mut field = read_gfc(&std::fs::read_to_string(goco)?)?.model.base;
    for (n, c) in [(2, -4.841_667_749_85e-4), (4, 7.903_037_335_11e-7), (6, -1.687_249_611_51e-9), (8, 3.460_524_683_94e-12)] {
        field.set(n, 0, field.c(n, 0) - c, 0.0);
    }
    for (n, m) in [(0, 0), (1, 0), (1, 1)] {
        field.set(n, m, 0.0, 0.0);
    }
    let field = field.truncated(180);
    let (lats, lons) = crate::grid_axes(2.0);
    let geoid = field.grid(&Functional::GeoidHeight, field.radius, &lats, &lons);
    write!(j, "\"geoid\":{{\"step\":2,\"rows\":{},\"cols\":{},\"values\":{}}},", lats.len(), lons.len(), list(geoid.iter().copied(), 1))?;

    // Coastlines (lon, lat) for the globe.
    let coast: Vec<String> = crate::coastlines().iter().map(|ring| {
        let pts: Vec<String> = ring.iter().map(|(lon, lat)| format!("[{lon:.2},{lat:.2}]")).collect();
        format!("[{}]", pts.join(","))
    }).collect();
    write!(j, "\"coast\":[{}],", coast.join(","))?;

    // GRACE-FO, every 30 s for 24 hours from the eruption (positions only).
    t1::UNBLIND.store(true, std::sync::atomic::Ordering::Relaxed);
    let orbit = t1::pair_orbit(l1b_dir, &["2022-01-15", "2022-01-16"], t1::ERUPTION, t1::ERUPTION + 86_400.0, 30.0)?;
    type Pick = fn(&(f64, t1::V3, t1::V3)) -> t1::V3;
    let coords = |pick: Pick| -> String {
        let pts: Vec<String> = orbit.iter().map(|o| { let p = pick(o); format!("[{:.1},{:.1},{:.1}]", p[0], p[1], p[2]) }).collect();
        format!("[{}]", pts.join(","))
    };
    write!(j, "\"orbit\":{{\"t0\":{},\"dt\":30,\"t\":{},\"c\":{},\"d\":{}}},", t1::ERUPTION,
        list(orbit.iter().map(|o| o.0 - t1::ERUPTION), 0), coords(|o| o.1), coords(|o| o.2))?;
    write!(j, "\"tonga\":{{\"lat\":-20.546,\"lon\":-175.390,\"speed\":318.2,\"sound_252k\":{:.1}}},", gravimetry::atmosphere::sound_speed(252.0))?;

    // LAGEOS frame dragging (geo5d).
    let node = geo5d::relativity::lense_thirring_node_rate(geo5d::constants::EARTH_ANGULAR_MOMENTUM, 12_270e3, 0.0045) * geo5d::relativity::MAS_PER_YEAR;
    write!(j, "\"lageos\":{{\"a_km\":12270,\"e\":0.0045,\"inc_deg\":109.84,\"node_mas_yr\":{node:.2}}},")?;

    // Gravitational-wave sensitivity (gwsens, satgrav): √S_n on a log grid, and published points.
    let freqs: Vec<f64> = (0..=240).map(|k| 10f64.powf(-9.0 + 18.0 * f64::from(k) / 240.0)).collect();
    let curve = |f: &dyn Fn(f64) -> f64| sci(freqs.iter().map(|&x| { let v = f(x); if v.is_finite() && v > 0.0 { v.sqrt() } else { 0.0 } }));
    write!(j, "\"gw\":{{\"freqs\":{},\"aligo\":{},\"et\":{},\"lisa\":{},", sci(freqs.iter().copied()),
        curve(&gwsens::curves::aligo), curve(&gwsens::curves::einstein_telescope), curve(&|f| gwsens::curves::lisa(f, gwsens::curves::Confusion::FourYears)))?;
    let aion = satgrav::gw::AtomGwDetector { baseline: 2000.0, t: 5.0, lmt: 40_000.0, phase_noise: 3e-6, k: satgrav::gw::K_SR87, diamonds: 1.0 };
    let points: Vec<String> = [
        ("GRACE-FO laser link as a detector", 1e-2, satgrav::gw::pair_strain(1e-9, 2.2e5), "satgrav"),
        ("AION-km atom interferometer", 0.1, aion.peak_strain(), "satgrav"),
    ].iter().map(|(n, f, h, src)| format!("{{\"name\":\"{n}\",\"f\":{f:e},\"h\":{h:.3e},\"src\":\"{src}\"}}"))
        .chain(gwsens::hf::REGISTRY.iter().filter(|l| l.quantity.starts_with("√S_h")).map(|l| format!("{{\"name\":\"{}\",\"f\":{:e},\"h\":{:.3e},\"src\":\"gwsens\"}}", l.name, (l.band.0 * l.band.1).sqrt(), l.value)))
        .collect();
    write!(j, "\"points\":[{}]}},", points.join(","))?;

    // Plasma (plasmoid) and spacetime (gravitics).
    let bennett = plasmoid::pinch::bennett_temperature(1e6, 1e20, 1.0);
    let ostriker = plasmoid::pinch::ostriker_line_mass(10.0, 2.3) * 3.0857e16 / 1.989e30;
    write!(j, "\"plasma\":{{\"sphere_eigenvalue\":{:.6},\"j11\":{:.6},\"bennett_ev_1ma\":{bennett:.0},\"ostriker_msun_pc\":{ostriker:.1},\"pease_braginskii_ma\":{:.2},\"magnetic_over_gravity\":{:.2e}}},",
        plasmoid::relaxation::SPHERE_EIGENVALUE, plasmoid::relaxation::J11,
        plasmoid::pinch::pease_braginskii(10.0, 1.0, plasmoid::pinch::Resistivity::Perpendicular) / 1e6,
        plasmoid::pinch::magnetic_over_gravitational(1e6, plasmoid::pinch::bennett_line_density(1e6, 2000.0) * 3.343_583_772e-27))?;
    let sun = 1.989e30;
    write!(j, "\"spacetime\":{{\"warp_100m_planck_kg\":{:.3e},\"warp_100m_1m_msun\":{:.3},\"schwinger_v_m\":{:.4e},\"unruh_k_per_ms2\":{:.4e},\"hawking_sun_k\":{:.3e},\"quantum_inequality_1ns\":{:.3e}}}",
        -gravitics::spacetime::alcubierre_energy(100.0, 100.0 * gravitics::PLANCK_LENGTH, 1.0),
        -gravitics::spacetime::alcubierre_energy(100.0, 1.0, 1.0) / sun,
        gravitics::strong_field::E_SCHWINGER, gravitics::spacetime::unruh_temperature(1.0),
        gravitics::spacetime::hawking_temperature(sun), gravitics::spacetime::quantum_inequality(1e-9))?;
    j.push('}');
    std::fs::create_dir_all(out.parent().unwrap_or(Path::new(".")))?;
    std::fs::write(out, &j)?;
    eprintln!("wrote {} ({} kB, {} orbit points)", out.display(), j.len() / 1024, orbit.len());
    // The page: the viewer's template with this data inlined.
    let page = include_str!("../../../viewer/orrery.html").replace("__ATLAS_DATA__", &j);
    let html = out.with_file_name("orrery.html");
    std::fs::write(&html, page)?;
    eprintln!("wrote {}", html.display());
    Ok(())
}
