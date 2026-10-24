//! Figures for `gravity-claims`: how far off each claim is, the energy a warp bubble needs, and
//! thrust per watt against the photon rocket. Every number comes from a crate function the crate's
//! tests check; the claimed values the crate does not hold are constants here, with sources.

use std::fmt::Write;

use gravity_claims::hfgw::{VACUUM_BREAKDOWN, shot_noise_time, surface_field};
use gravity_claims::propulsion::{ballistic_thrust, ion_wind_thrust, mach_effect_thrust, photon_thrust, shawyer_thrust};
use gravity_claims::spacetime::alcubierre_energy;
use gravity_claims::superconductors::{PODKLETNOV_CLAIM, TAJMAR_CLAIM, WEIGHT_BOUNDS, marx_bank, pendulum_impulse, ring_frame_dragging};
use gravity_claims::{AMU, C, E_CHARGE, PLANCK_LENGTH};

use crate::svg::{esc, open, pow10};

const LI_CLAIMED_TIME: f64 = 3e3; // s, Li et al. (PRD 2009)
const EAGLEWORKS: f64 = 1.2e-3 / 1e3; // N/W claimed (White et al. 2017)
const WOODWARD_REPORTED: f64 = 100e-6; // N at ~30 W (Monette et al. 2021, of Woodward's traces)
const DRESDEN_EMDRIVE: f64 = 30e-9 / 11.0; // N/W upper limit (Tajmar et al. 2022)
const SOLAR_MASS: f64 = 1.989e30;

pub fn figures() -> Vec<(&'static str, String)> {
    vec![("claims_gap.svg", gap()), ("warp_energy.svg", warp()), ("thrust_per_watt.svg", thrust())]
}

/// `3.6×10²¹`, or a plain number between 1 and 1000.
fn sci(v: f64) -> String {
    let e = v.abs().log10().floor() as i32;
    if (0..3).contains(&e) {
        return format!("{v:.0}");
    }
    format!("{:.1}×{}", v / 10f64.powi(e), pow10(e))
}

fn gap() -> String {
    let swing = pendulum_impulse(0.8, 0.142, 0.0185, 1e-4);
    let strictest = WEIGHT_BOUNDS.iter().map(|b| b.1).fold(f64::INFINITY, f64::min);
    let woodward = mach_effect_thrust(2.0 * std::f64::consts::PI * 32_350.0, 0.042, 3.2e-10, 4e-11, 0.019, 100.0, 7.9e3);
    let rows = [
        ("Tajmar: frame dragging near a spinning superconductor", "claim ÷ general relativity", TAJMAR_CLAIM / ring_frame_dragging(0.349, 0.075, 0.069)),
        ("Li–Baker: detecting a GHz gravitational wave", "shot-noise time ÷ claimed time", shot_noise_time(2.56e-18, 2.0) / LI_CLAIMED_TIME),
        ("Podkletnov: impulse-generator pendulum", "momentum ÷ the momentum of radiation with that energy", 0.0185 * swing.speed * C / swing.energy),
        ("Pais: charged-disc wave source", "surface field ÷ the breakdown limit his patent states", surface_field(5e4) / VACUUM_BREAKDOWN),
        ("Podkletnov: weight loss above a spinning disc", "claim ÷ strictest replication bound", PODKLETNOV_CLAIM / strictest),
        ("Woodward: Mach-effect thruster", "reported thrust ÷ his own formula", WOODWARD_REPORTED / woodward),
        ("EmDrive (Eagleworks)", "claimed thrust per watt ÷ a photon rocket", EAGLEWORKS / photon_thrust(1.0)),
        ("Podkletnov: impulse-generator charge", "charge needed ÷ what its Marx bank holds", 0.1 / marx_bank(25e-9, 20.0, 1e5).1),
    ];
    let (w, left, right, top, row) = (960.0, 360.0, 880.0, 120.0, 52.0);
    let h = top + row * rows.len() as f64 + 56.0;
    let x = |v: f64| left + (right - left) * v.log10() / 22.0;
    let mut s = open(w, h, "How far off are the claims?");
    let _ = write!(s, r#"<text x="24" y="34" class="title">How far off are the claims?</text>
<text x="24" y="56" class="subtitle">Each claim's own numbers, divided by what general relativity, its replications or its own formula allow. Log scale: each gridline is 10⁵.</text>
"#);
    for e in (0..=20).step_by(5) {
        let xx = x(10f64.powi(e));
        let _ = writeln!(s, r#"<line x1="{xx:.1}" y1="{:.1}" x2="{xx:.1}" y2="{:.1}" class="grid"/><text x="{xx:.1}" y="{:.1}" class="tick" text-anchor="middle">{}</text>"#,
            top - 12.0, h - 44.0, top - 20.0, if e == 0 { "1".into() } else { pow10(e) });
    }
    for (i, (name, against, gap)) in rows.iter().enumerate() {
        let y = top + row * i as f64;
        let end = x(*gap);
        let _ = writeln!(s, r#"<text x="24" y="{:.1}" style="fill: var(--ink); font-weight: 600">{}</text><text x="24" y="{:.1}" class="tick">{}</text>"#,
            y + 14.0, esc(name), y + 31.0, esc(against));
        let _ = writeln!(s, r#"<rect x="{left}" y="{:.1}" width="{:.1}" height="14" rx="4" fill="var(--s1)"><title>{}: {}</title></rect><text x="{:.1}" y="{:.1}" style="fill: var(--ink)">{}</text>"#,
            y + 9.0, (end - left).max(2.0), esc(name), sci(*gap), end + 8.0, y + 21.0, sci(*gap));
    }
    let _ = writeln!(s, r#"<line x1="{left}" y1="{:.1}" x2="{left}" y2="{:.1}" class="axis"/>"#, top - 12.0, h - 44.0);
    let _ = writeln!(s, r#"<text x="24" y="{:.1}" class="tick">From the gravity-claims crate, whose tests pin every number. Sources in the research report.</text>"#, h - 16.0);
    s.push_str("</svg>\n");
    s
}

fn warp() -> String {
    let (w, h, left, right, top, bottom) = (960.0, 560.0, 92.0, 900.0, 96.0, 500.0);
    let (x_lo, x_hi, y_lo, y_hi) = (-35.0, 3.0, 20.0, 70.0); // log10 of wall (m) and mass (kg)
    let x = |d: f64| left + (right - left) * (d.log10() - x_lo) / (x_hi - x_lo);
    let y = |m: f64| bottom - (bottom - top) * (m.log10() - y_lo) / (y_hi - y_lo);
    let mut s = open(w, h, "The energy a warp bubble needs");
    let _ = write!(s, r#"<text x="24" y="34" class="title">The energy a warp bubble needs</text>
<text x="24" y="56" class="subtitle">Negative energy, as a mass, of an Alcubierre bubble moving at light speed, against its wall thickness (Pfenning &amp; Ford 1997).</text>
"#);
    // The walls quantum field theory allows at v = c: Δ ≲ 100 Planck lengths.
    let qi = x(100.0 * PLANCK_LENGTH);
    let _ = writeln!(s, r#"<rect x="{left}" y="{top}" width="{:.1}" height="{:.1}" fill="var(--s3)" opacity="0.15"/><text x="{:.1}" y="{:.1}" class="tick">walls the quantum</text><text x="{:.1}" y="{:.1}" class="tick">inequality allows</text>"#,
        qi - left, bottom - top, left + 6.0, bottom - 26.0, left + 6.0, bottom - 12.0);
    for e in (20..=70).step_by(10) {
        let yy = y(10f64.powi(e));
        let _ = writeln!(s, r#"<line x1="{left}" y1="{yy:.1}" x2="{right}" y2="{yy:.1}" class="grid"/><text x="{:.1}" y="{:.1}" class="tick" text-anchor="end">{} kg</text>"#, left - 8.0, yy + 4.0, pow10(e));
    }
    for (e, label) in [(-35, "Planck length"), (-15, "proton"), (-10, "atom"), (-6, "µm"), (0, "1 m"), (3, "1 km")] {
        let xx = x(10f64.powi(e));
        let _ = writeln!(s, r#"<line x1="{xx:.1}" y1="{top}" x2="{xx:.1}" y2="{bottom}" class="grid"/><text x="{xx:.1}" y="{:.1}" class="tick" text-anchor="middle">{label}</text>"#, bottom + 18.0);
    }
    let _ = writeln!(s, r#"<line x1="{left}" y1="{bottom}" x2="{right}" y2="{bottom}" class="axis"/><text x="{:.1}" y="{:.1}" class="tick" text-anchor="middle">wall thickness</text>"#, (left + right) / 2.0, bottom + 38.0);
    // Reference masses.
    // Labels go where the curves are not: big masses at the right, small ones at the left.
    for (m, label) in [(5.97e24, "Earth"), (SOLAR_MASS, "Sun"), (2e42, "Milky Way, ≈10¹² Suns"), (1.5e53, "ordinary matter in the observable universe, ≈10⁵³ kg")] {
        let yy = y(m);
        let (tx, anchor) = if m > 1e40 { (right - 4.0, "end") } else { (qi + 8.0, "start") };
        let _ = writeln!(s, r#"<line x1="{left}" y1="{yy:.1}" x2="{right}" y2="{yy:.1}" stroke="var(--muted)" stroke-dasharray="4 4"/><text x="{tx:.1}" y="{:.1}" class="tick" text-anchor="{anchor}">{}</text>"#, yy - 6.0, esc(label));
    }
    let walls: Vec<f64> = (0..=200).map(|i| 10f64.powf(x_lo + (x_hi - x_lo) * i as f64 / 200.0)).collect();
    let series = [(1.0, "s2", "1 m bubble"), (100.0, "s1", "100 m bubble"), (1e4, "s4", "10 km bubble")];
    let mut lx = 24.0;
    for (r, class, label) in series {
        let d: String = walls.iter().map(|&wall| (wall, -alcubierre_energy(r, wall, 1.0)))
            .filter(|&(_, m)| m.log10() <= y_hi && m.log10() >= y_lo)
            .enumerate().map(|(j, (wall, m))| format!("{}{:.1},{:.1}", if j == 0 { 'M' } else { 'L' }, x(wall), y(m))).collect();
        let _ = writeln!(s, r#"<path d="{d}" class="line {class}"><title>{label}: E ∝ R²/Δ</title></path>"#);
        let _ = writeln!(s, r#"<line x1="{lx:.1}" y1="78" x2="{:.1}" y2="78" class="line {class}"/><text x="{:.1}" y="82">{label}</text>"#, lx + 16.0, lx + 22.0);
        lx += 130.0;
    }
    // The one point the quantum inequality allows for the 100 m bubble.
    let (wall, m) = (100.0 * PLANCK_LENGTH, -alcubierre_energy(100.0, 100.0 * PLANCK_LENGTH, 1.0));
    let _ = writeln!(s, r#"<circle cx="{:.1}" cy="{:.1}" r="5" fill="var(--s1)" stroke="var(--surface)" stroke-width="2"><title>100 m bubble with a 100 Planck-length wall: {} kg</title></circle><text x="{:.1}" y="{:.1}" style="fill: var(--ink)">{} kg: {} times all ordinary matter</text>"#,
        x(wall), y(m), sci(m), x(wall) + 10.0, y(m) - 8.0, sci(m), sci(m / 1.5e53));
    s.push_str("</svg>\n");
    s
}

fn thrust() -> String {
    let n2 = 28.0134 * AMU;
    // ARL's second lifter: 1.5 mA across 6 cm at 30 kV (45 W).
    let lifter_air = ion_wind_thrust(1.5e-3, 0.06, 2.5e-4) / 45.0;
    let lifter_vacuum = ballistic_thrust(1.5e-3, 3e4, n2, E_CHARGE) / 45.0;
    let rows = [
        ("Lifter in air (ion wind)", "standard", lifter_air),
        ("Same ions in vacuum (an ion thruster)", "standard", lifter_vacuum),
        ("Woodward thruster, reported (~30 W)", "claim", WOODWARD_REPORTED / 30.0),
        ("Shawyer's formula at Dresden's settings", "claim", shawyer_thrust(11.0, 741.0, 0.56) / 11.0),
        ("EmDrive, Eagleworks", "claim", EAGLEWORKS),
        ("Photon rocket, P/c", "standard", photon_thrust(1.0)),
        ("EmDrive, Dresden upper limit", "measured", DRESDEN_EMDRIVE),
    ];
    let (w, left, right, top, row) = (960.0, 340.0, 900.0, 120.0, 40.0);
    let h = top + row * rows.len() as f64 + 60.0;
    let (lo, hi) = (-9.0, -1.0); // log10 N/W
    let x = |v: f64| left + (right - left) * (v.log10() - lo) / (hi - lo);
    let mut s = open(w, h, "Thrust per watt");
    let _ = write!(s, r#"<text x="24" y="34" class="title">Thrust per watt</text>
<text x="24" y="56" class="subtitle">"Propellantless" thrusters against light and ion wind, in newtons per watt. Nothing that pushes on nothing beats light.</text>
"#);
    let mut lx = 24.0;
    for (class, label) in [("s1", "standard physics"), ("s2", "claimed"), ("s3", "measured limit")] {
        let _ = writeln!(s, r#"<circle cx="{:.1}" cy="80" r="5" fill="var(--{class})"/><text x="{:.1}" y="84">{label}</text>"#, lx + 5.0, lx + 16.0);
        lx += 150.0;
    }
    for e in (lo as i32)..=(hi as i32) {
        let xx = x(10f64.powi(e));
        let _ = writeln!(s, r#"<line x1="{xx:.1}" y1="{:.1}" x2="{xx:.1}" y2="{:.1}" class="grid"/><text x="{xx:.1}" y="{:.1}" class="tick" text-anchor="middle">{}</text>"#,
            top - 10.0, h - 50.0, h - 32.0, pow10(e));
    }
    let photon = x(photon_thrust(1.0));
    let _ = writeln!(s, r#"<line x1="{photon:.1}" y1="{:.1}" x2="{photon:.1}" y2="{:.1}" stroke="var(--ink2)" stroke-dasharray="4 4"/>"#, top - 10.0, h - 50.0);
    for (i, (name, kind, v)) in rows.iter().enumerate() {
        let y = top + row * i as f64 + 10.0;
        let class = match *kind { "standard" => "s1", "claim" => "s2", _ => "s3" };
        let _ = writeln!(s, r#"<text x="24" y="{:.1}" style="fill: var(--ink)">{}</text><line x1="{left}" y1="{y:.1}" x2="{right}" y2="{y:.1}" class="grid"/>"#, y + 4.0, esc(name));
        let (tx, anchor) = if x(*v) > right - 120.0 { (x(*v) - 12.0, "end") } else { (x(*v) + 12.0, "start") };
        let _ = writeln!(s, r#"<circle cx="{:.1}" cy="{y:.1}" r="6" fill="var(--{class})" stroke="var(--surface)" stroke-width="2"><title>{}: {} N/W</title></circle><text x="{tx:.1}" y="{:.1}" text-anchor="{anchor}">{} N/W</text>"#,
            x(*v), esc(name), sci(*v), y + 4.0, sci(*v));
    }
    s.push_str("</svg>\n");
    s
}

#[cfg(test)]
mod tests {
    #[test]
    fn scientific_labels() {
        assert_eq!(super::sci(3.65e21), "3.6×10²¹");
        assert_eq!(super::sci(360.0), "360");
        assert_eq!(super::sci(2.7e-6), "2.7×10⁻⁶");
    }
}
