//! Space weather from the solar wind to the power grid: where the magnetosphere stands, how much
//! energy a storm stores, how far the plasmasphere reaches, and the electric field a storm
//! induces in the ground.
//!
//! Planners use these numbers directly. NERC's TPL-007 standard sizes US grid defences for a
//! once-in-100-years peak geoelectric field of 8 V/km at 60° geomagnetic latitude
//! ([NERC](https://www.nerc.com/pa/Stand/Project201303GeomagneticDisturbanceMitigation/common_questions_June12.pdf));
//! the May 2024 storm reached Dst −412 nT ([arXiv:2407.07665](https://arxiv.org/pdf/2407.07665)).
//!
//! ```
//! use magnetics::space_weather::{dynamic_pressure, standoff};
//!
//! // A typical solar wind, 5 protons per cm³ at 400 km/s, against Earth's 29 733 nT dipole:
//! // the magnetopause stands 10.1 Earth radii upstream.
//! let r = standoff(29_733.4e-9, dynamic_pressure(5e6, 400e3), true);
//! assert!((r - 10.08).abs() < 0.01);
//! ```

use core::f64::consts::TAU;

use metrology::{M_P, MU_0};

/// The Sun's sidereal rotation rate at the equator (Carrington period 25.38 days), rad/s.
pub const SOLAR_ROTATION: f64 = TAU / (25.38 * 86_400.0);
/// NERC TPL-007's benchmark peak geoelectric field, V/m (8 V/km), before the latitude and ground
/// scaling factors.
pub const NERC_BENCHMARK: f64 = 8e-3;
/// NERC TPL-007's supplemental (localised) peak geoelectric field, V/m (12 V/km).
pub const NERC_SUPPLEMENTAL: f64 = 12e-3;

/// Solar-wind dynamic pressure `n m_p v²`, Pa, for proton density `n` (m⁻³) and speed `v` (m/s).
pub fn dynamic_pressure(n: f64, v: f64) -> f64 {
    n * M_P * v * v
}

/// Magnetopause standoff distance in planetary radii, from pressure balance between the
/// planet's equatorial surface dipole field `b0` (T) and the wind's dynamic pressure `p` (Pa):
/// `(2B₀²/(μ₀p))^(1/6)` with the factor-2 compression of the field at the boundary, or
/// `(B₀²/(2μ₀p))^(1/6)` without it (Chapman–Ferraro). The two differ by `4^(1/6)` = 1.26.
pub fn standoff(b0: f64, p: f64, compressed: bool) -> f64 {
    let k = if compressed { 2.0 } else { 0.5 };
    (k * b0 * b0 / (MU_0 * p)).powf(1.0 / 6.0)
}

/// The Shue et al. (1998) empirical subsolar magnetopause, Earth radii, for the interplanetary
/// field's north component `bz_nt` (nT) and dynamic pressure `p_npa` (nPa).
pub fn shue1998(bz_nt: f64, p_npa: f64) -> f64 {
    (10.22 + 1.29 * (0.184 * (bz_nt + 8.14)).tanh()) * p_npa.powf(-1.0 / 6.6)
}

/// Angle (radians) between the interplanetary magnetic field and the radial direction at
/// distance `r` (m) for solar-wind speed `v` (m/s): the Parker spiral, `atan(Ωr/v)`.
pub fn parker_angle(v: f64, r: f64) -> f64 {
    (SOLAR_ROTATION * r / v).atan()
}

/// Energy (J) of the ring current behind a storm depression `dst` (T) on a planet of dipole
/// moment `moment` (A m²): the Dessler–Parker–Sckopke relation, `M|ΔB|/2`.
pub fn ring_current_energy(dst: f64, moment: f64) -> f64 {
    moment * dst.abs() / 2.0
}

/// The plasmapause L-shell from the maximum Kp of the previous day (Carpenter & Anderson 1992):
/// `5.6 − 0.46 Kp`. `None` outside the fit's range, 2.25 < L < 8.
pub fn plasmapause(kp_max: f64) -> Option<f64> {
    let l = 5.6 - 0.46 * kp_max;
    (l > 2.25 && l < 8.0).then_some(l)
}

/// A dipole's field strength (same units as `b0`, the surface equatorial field) on the L-shell
/// `l` at magnetic latitude `mlat` (radians): `B₀L⁻³√(1 + 3 sin²λ)/cos⁶λ`, at radius `L cos²λ`.
pub fn dipole_field(b0: f64, l: f64, mlat: f64) -> f64 {
    b0 / l.powi(3) * (1.0 + 3.0 * mlat.sin().powi(2)).sqrt() / mlat.cos().powi(6)
}

/// Invariant latitude (radians) where the field line of L-shell `l` meets the surface: `acos L^(−½)`.
pub fn invariant_latitude(l: f64) -> f64 {
    l.powf(-0.5).acos()
}

/// Geoelectric field amplitude (V/m) and skin depth (m) induced over a uniform ground of
/// conductivity `sigma` (S/m) by a magnetic variation of amplitude `b` (T) and period `period`
/// (s): `|E| = √(ω/(μ₀σ))|B|`, leading B by 45°, and `δ = √(2/(ωμ₀σ))`. (The form `√(ωμ₀/σ)·B`
/// sometimes quoted holds only if "B" means `H = B/μ₀`.)
pub fn geoelectric_halfspace(b: f64, period: f64, sigma: f64) -> (f64, f64) {
    let omega = TAU / period;
    ((omega / (MU_0 * sigma)).sqrt() * b, (2.0 / (omega * MU_0 * sigma)).sqrt())
}

/// Apparent resistivity (Ω m) from the amplitudes of the electric field `e` (V/m) and magnetic
/// field `b` (T) at period `period` (s): `μ₀|E|²/(ω|B|²)`, the magnetotelluric inverse of
/// [`geoelectric_halfspace`]. In practical units this is `0.2T(E[mV/km]/B[nT])²`.
pub fn apparent_resistivity(e: f64, b: f64, period: f64) -> f64 {
    MU_0 * e * e / (TAU / period * b * b)
}

/// The radius (in L units, planetary radii) of a dipole field line of L-shell `l` at magnetic
/// latitude `mlat` (radians): `L cos²λ`.
pub fn field_line_radius(l: f64, mlat: f64) -> f64 {
    l * mlat.cos().powi(2)
}

#[cfg(test)]
mod tests {
    use super::*;
    use metrology::AU;

    fn close(a: f64, b: f64, rel: f64) -> bool {
        (a / b - 1.0).abs() < rel
    }

    #[test]
    fn magnetopause() {
        let p = dynamic_pressure(5e6, 400e3);
        assert!(close(p, 1.338_098e-9, 1e-6));
        let b0 = 29_733.4e-9; // IGRF-14 2025 dipole
        assert!(close(standoff(b0, p, true), 10.084_09, 1e-6) && close(standoff(b0, p, false), 8.003_75, 1e-6));
        assert!(close(standoff(b0, p, true) / standoff(b0, p, false), 4f64.powf(1.0 / 6.0), 1e-12));
        assert!(close(shue1998(0.0, 1.338), 10.895_67, 1e-5) && close(shue1998(-20.0, 20.0), 5.692_45, 1e-5));
    }

    #[test]
    fn heliosphere_and_storms() {
        assert!(close(parker_angle(400e3, AU).to_degrees(), 46.979_98, 1e-6));
        assert!(close(parker_angle(800e3, AU).to_degrees(), 28.182_80, 1e-6));
        // DPS with IGRF-14's 2025 moment: 3.845×10¹³ J per nT; May 2024 (−412 nT), March 1989 (−589 nT).
        let m = 7.690e22;
        assert!(close(ring_current_energy(-1e-9, m), 3.845e13, 1e-12));
        assert!(close(ring_current_energy(-412e-9, m), 1.584e16, 1e-3) && close(ring_current_energy(-589e-9, m), 2.265e16, 1e-3));
    }

    #[test]
    fn inner_magnetosphere() {
        assert!(close(plasmapause(2.0).unwrap(), 4.68, 1e-12) && close(plasmapause(6.0).unwrap(), 2.84, 1e-12));
        assert!(plasmapause(9.0).is_none());
        let b0 = 29_733.4e-9;
        assert!(close(dipole_field(b0, 4.0, 0.0), 464.6e-9, 1e-4));
        assert!(close(dipole_field(b0, 4.0, 30f64.to_radians()), 1456.8e-9, 1e-4));
        assert!(close(field_line_radius(4.0, 30f64.to_radians()), 3.0, 1e-12));
        assert!(close(invariant_latitude(4.0).to_degrees(), 60.0, 1e-12));
    }

    #[test]
    fn ground_electric_field() {
        let (e, d) = geoelectric_halfspace(500e-9, 300.0, 1e-3);
        assert!(close(e * 1e3, 2.04, 3e-3) && close(d / 1e3, 276.0, 2e-3));
        // A 1000 nT, 300 s variation over resistive (10⁻⁴ S/m) ground: 12.9 V/km, NERC's order.
        let (e, d) = geoelectric_halfspace(1000e-9, 300.0, 1e-4);
        assert!(close(e * 1e3, 12.9, 3e-3) && close(d / 1e3, 872.0, 2e-3) && e > NERC_BENCHMARK);
        // The inverse recovers the ground resistivity, and the practical form agrees.
        let (e, _) = geoelectric_halfspace(100e-9, 100.0, 1e-3);
        assert!(close(apparent_resistivity(e, 100e-9, 100.0), 1e3, 1e-12));
        assert!(close(apparent_resistivity(e, 100e-9, 100.0), 0.2 * 100.0 * (e * 1e6 / 100.0).powi(2), 1e-6));
    }
}
