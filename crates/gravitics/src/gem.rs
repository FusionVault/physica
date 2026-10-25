//! Gravitoelectromagnetism: linearised general relativity written like Maxwell's equations. A
//! spinning mass has a gravitomagnetic dipole field, exactly as a spinning charge has a magnetic
//! one, and it drags gyroscopes and orbits (Gravity Probe B, LAGEOS and LARES).
//!
//! The analogy's factors of 2 and 4 are conventions, with one physical exception: gravity is
//! spin-2, so a mass couples to the gravitomagnetic field twice as strongly as a charge to the
//! magnetic one. Three common conventions ([`Convention`]) put that factor in different places
//! and give the same force. The field outside a spinning body has a scalar potential, so
//! [`GravitomagneticDipole`] is a [`Field`] beside gravity and the geomagnetic main field.
//!
//! Blackett's 1947 hypothesis, that every rotating body is magnetised in proportion to its spin
//! (`μ = (β/2)√(4πε₀G)·J`), is here as a calculator: its constant β should be the same for every
//! body, and the measured values span a factor of 10⁶, with Mars below 2.4×10⁻⁵.
//!
//! ```
//! use geopotential::Field;
//! use gravitics::gem::GravitomagneticDipole;
//!
//! // Earth (J = 5.86×10³³ kg m²/s): 5.04×10⁻⁶ m/s² at the equator, twice that at the pole.
//! let earth = GravitomagneticDipole { j: 5.86e33 };
//! let r = 6.371e6;
//! assert!((earth.magnitude(r, 0.0, 0.0) / 5.045e-6 - 1.0).abs() < 1e-3);
//! assert!((earth.magnitude(r, 1.5707, 0.0) / earth.magnitude(r, 0.0, 0.0) - 2.0).abs() < 1e-6);
//! ```

use core::f64::consts::PI;

use geopotential::Field;
use metrology::{C, EPS_0, G};

/// `Gm_e²/(ħc)` = 1.75×10⁻⁴⁵: the gravitational fine-structure constant of the electron, and the
/// factor by which a condensate's gravitomagnetic response is suppressed below its magnetic one
/// (DeWitt 1966). It is why spinning a superconductor changes no weight that a balance can see.
pub const ELECTRON_GRAVITATIONAL_COUPLING: f64 = G * metrology::M_E * metrology::M_E / (metrology::HBAR * C);

/// The gravitational analogue of the vacuum permittivity, `1/(4πG)`, kg s² m⁻³.
pub const GM_PERMITTIVITY: f64 = 1.0 / (4.0 * PI * G);
/// The gravitational analogue of the vacuum permeability, `4πG/c²`, m/kg; with
/// [`GM_PERMITTIVITY`], `1/√(εμ) = c`.
pub const GM_PERMEABILITY: f64 = 4.0 * PI * G / (C * C);

/// Where a gravitoelectromagnetic convention puts the spin-2 factor.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Convention {
    /// Mashhoon: field `B = (G/c)[3(J·r̂)r̂ − J]/r³` (m/s²), force per mass `−2(v/c)×B`.
    Mashhoon,
    /// Ruggiero & Tartaglia: field `B_G = −2B` (Mashhoon's B), force per mass `(v/c)×B_G`.
    RuggieroTartaglia,
    /// The Maxwell-like SI form: field `B_g = GJ/(2c²r³)` (s⁻¹), force per mass `4v×B_g`.
    MaxwellSi,
}

/// The equatorial gravitomagnetic field magnitude of a body of angular momentum `j` (kg m²/s) at
/// distance `r` (m), in `convention`'s own units.
pub fn equatorial_field(j: f64, r: f64, convention: Convention) -> f64 {
    let mashhoon = G * j / (C * r.powi(3));
    match convention {
        Convention::Mashhoon => mashhoon,
        Convention::RuggieroTartaglia => 2.0 * mashhoon,
        Convention::MaxwellSi => G * j / (2.0 * C * C * r.powi(3)),
    }
}

/// The magnitude of the gravitomagnetic force per unit mass (m/s²) on a body moving at `v` (m/s)
/// across a field `b` given in `convention`'s units. Every convention gives the same number.
pub fn force_per_mass(v: f64, b: f64, convention: Convention) -> f64 {
    match convention {
        Convention::Mashhoon => 2.0 * v * b / C,
        Convention::RuggieroTartaglia => v * b / C,
        Convention::MaxwellSi => 4.0 * v * b,
    }
}

/// The gravitomagnetic dipole field of a body spinning about the polar axis with angular
/// momentum `j` (kg m²/s), in Mashhoon's convention: potential `ψ = (G/c)J sin φ/r²` (m²/s) and
/// field `B = −∇ψ` (m/s²). A gyroscope precesses at `|B|/c`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GravitomagneticDipole {
    /// Angular momentum, kg m²/s.
    pub j: f64,
}

impl Field for GravitomagneticDipole {
    fn potential(&self, r: f64, lat: f64, _lon: f64) -> f64 {
        G * self.j / C * lat.sin() / (r * r)
    }

    fn vector(&self, r: f64, lat: f64, _lon: f64) -> [f64; 3] {
        let k = G * self.j / (C * r.powi(3));
        [-k * lat.cos(), 0.0, 2.0 * k * lat.sin()]
    }
}

/// Radians per second to milliarcseconds per year, how precession rates are quoted.
pub const MAS_PER_YEAR: f64 = metrology::YEAR * 180.0 / PI * 3.6e6;

/// The Lense–Thirring drag of an orbit's node, `2GJ/(c²a³(1 − e²)^(3/2))` (rad/s), about a body
/// of spin angular momentum `j` (kg m²/s), for semi-major axis `a` (m) and eccentricity `e`: about
/// 31 mas/yr for LAGEOS and 118 for LARES (Ciufolini et al. 2013), measured to about 0.1% by
/// LARES-2 and LAGEOS. Prograde for every inclination.
pub fn lense_thirring_node_rate(j: f64, a: f64, e: f64) -> f64 {
    2.0 * G * j / (C * C * a.powi(3) * (1.0 - e * e).powf(1.5))
}

/// The Lense–Thirring drag of the perigee, `−6GJ cos i/(c²a³(1 − e²)^(3/2))` (rad/s), at
/// `inclination` (rad): three times the node rate, backwards, for a prograde orbit.
pub fn lense_thirring_perigee_rate(j: f64, a: f64, e: f64, inclination: f64) -> f64 {
    -3.0 * inclination.cos() * lense_thirring_node_rate(j, a, e)
}

/// The frame-dragging precession (rad/s) of a gyroscope on a circular polar orbit of radius `a`
/// (m), averaged over the orbit: `GJ/(2c²a³)`, half the equatorial [`GravitomagneticDipole`] field
/// over `c`, since Schiff's `(G/c²r³)[3(J·r̂)r̂ − J]` averages to `J/2` round a polar orbit. Gravity
/// Probe B's 39 mas/yr.
pub fn gyroscope_frame_dragging(j: f64, a: f64) -> f64 {
    equatorial_field(j, a, Convention::Mashhoon) / (2.0 * C)
}

/// `√(4πε₀G)`, C/kg: the charge-to-mass ratio at which electric repulsion balances gravity. It
/// is the extremal Reissner–Nordström ratio, the weak-gravity-conjecture scale and the constant in
/// Blackett's hypothesis.
pub fn extremal_charge_to_mass() -> f64 {
    (4.0 * PI * EPS_0 * G).sqrt()
}

/// A particle's charge-to-mass ratio over the extremal one: `q/(m√(4πε₀G))` (2×10²¹ for the
/// electron, which is why electromagnetism dwarfs gravity between particles).
pub fn wgc_ratio(q: f64, m: f64) -> f64 {
    q / (m * extremal_charge_to_mass())
}

/// Blackett's constant `β = 2μ/(√(4πε₀G)·J)` for a body of magnetic moment `mu` (A m²) and spin
/// angular momentum `j` (kg m²/s). The hypothesis predicts one β for every body. Status: refuted.
pub fn blackett_beta(mu: f64, j: f64) -> f64 {
    2.0 * mu / (extremal_charge_to_mass() * j)
}

/// The magnetic moment (A m²) Blackett's hypothesis predicts for spin `j` and constant `beta`:
/// `(β/2)√(4πε₀G)·J`.
pub fn blackett_moment(j: f64, beta: f64) -> f64 {
    beta / 2.0 * extremal_charge_to_mass() * j
}

#[cfg(test)]
mod tests {
    use super::*;
    use metrology::{E_CHARGE, M_E, M_P};

    fn close(a: f64, b: f64, rel: f64) -> bool {
        (a / b - 1.0).abs() < rel
    }

    #[test]
    fn earth_dipole() {
        let earth = GravitomagneticDipole { j: 5.86e33 };
        let r = 6.371e6;
        assert!(close(earth.magnitude(r, 0.0, 0.0), 5.044_984e-6, 1e-6));
        assert!(close(earth.magnitude(r, 0.0, 0.0) / C, 1.682_83e-14, 1e-5)); // gyroscope precession, rad/s
        // The field is minus the gradient of the potential (central differences in r and latitude).
        let (lat, h) = (0.4, 1e-5);
        let up = -(earth.potential(r * (1.0 + h), lat, 0.0) - earth.potential(r * (1.0 - h), lat, 0.0)) / (2.0 * r * h);
        let north = -(earth.potential(r, lat + h, 0.0) - earth.potential(r, lat - h, 0.0)) / (2.0 * h * r);
        let [n, e, u] = earth.vector(r, lat, 0.0);
        assert!(close(u, up, 1e-8) && close(n, north, 1e-8) && e == 0.0);
    }

    #[test]
    fn conventions_agree() {
        // A satellite at 7.9 km/s across Earth's equatorial field: 2.66×10⁻¹⁰ m/s² in all three.
        let (j, r, v) = (5.86e33, 6.371e6, 7.9e3);
        let forces = [Convention::Mashhoon, Convention::RuggieroTartaglia, Convention::MaxwellSi].map(|c| force_per_mass(v, equatorial_field(j, r, c), c));
        assert!(close(forces[0], 2.658_864e-10, 1e-6));
        assert!(close(forces[1], forces[0], 1e-12) && close(forces[2], forces[0], 1e-12));
        assert!(close(1.0 / (GM_PERMITTIVITY * GM_PERMEABILITY).sqrt(), C, 1e-12));
        assert!(close(GM_PERMITTIVITY, 1.192_297e9, 1e-6) && close(GM_PERMEABILITY, 9.331_99e-27, 1e-5));
    }

    #[test]
    fn frame_dragging_of_orbits_and_gyroscopes() {
        use metrology::{EARTH_ANGULAR_MOMENTUM as J, GM_EARTH};
        // Ciufolini et al. (2013, arXiv:1309.1699): LAGEOS at 5900 km altitude, "approximately
        // 31 mas/yr" (30.67 for a = 12 270 km, e = 0.0045); LARES at a = 7820 km, e = 0.0007,
        // "approximately 118". The perigee drags back three times as fast at i = 0.
        let lageos = lense_thirring_node_rate(J, 12_270e3, 0.0045) * MAS_PER_YEAR;
        assert!(close(lageos, 30.67, 2e-4) && close(lageos, 31.0, 0.02));
        assert!(close(lense_thirring_node_rate(J, 7820e3, 0.0007) * MAS_PER_YEAR, 118.0, 0.01));
        let i = 109.8f64.to_radians();
        assert!(close(lense_thirring_perigee_rate(J, 12_270e3, 0.0045, i) / lense_thirring_node_rate(J, 12_270e3, 0.0045), -3.0 * i.cos(), 1e-12));
        // Gravity Probe B (Everitt et al. 2011): a 642 km polar orbit, and the predicted
        // frame-dragging rate towards IM Pegasi (declination 16.84°) of 39.2 mas/yr, with its
        // geodetic 6606.1; the bare formulas over the mean radius give both to 0.5%.
        let a = 6371e3 + 642e3;
        let dragging = gyroscope_frame_dragging(J, a) * MAS_PER_YEAR;
        assert!(close(dragging * 16.84f64.to_radians().cos(), 39.2, 5e-3), "{dragging}");
        assert!(close(gravitational::energy::geodetic_precession(GM_EARTH, a) * MAS_PER_YEAR, 6606.1, 6e-3));
        // Half the dipole's equatorial field over c, and 1.68×10⁻¹⁴ rad/s at the equator's surface.
        let earth = GravitomagneticDipole { j: J };
        assert!(close(gyroscope_frame_dragging(J, a), 0.5 * earth.magnitude(a, 0.0, 0.0) / C, 1e-12));
    }

    #[test]
    fn blackett_and_extremal_ratio() {
        assert!(close(extremal_charge_to_mass(), 8.617_52e-11, 1e-5));
        assert!(close(ELECTRON_GRAVITATIONAL_COUPLING, 1.751_81e-45, 1e-5));
        // The electromagnetic fine-structure constant is 42 orders larger.
        assert!(close(metrology::ALPHA / ELECTRON_GRAVITATIONAL_COUPLING, 4.165_57e42, 1e-4));
        assert!(close(wgc_ratio(E_CHARGE, M_E), 2.040_98e21, 1e-5) && close(wgc_ratio(E_CHARGE, M_P), 1.111_55e18, 1e-5));
        // Earth: μ = 7.91×10²² A m², J = 5.87×10³³ gives β = 0.31; IGRF-14's 2025 moment, 0.30.
        assert!(close(blackett_beta(7.91e22, 5.87e33), 0.312_742, 1e-5));
        assert!(close(blackett_beta(7.690e22, 5.87e33), 0.304_044, 1e-5));
        assert!(close(blackett_moment(5.87e33, blackett_beta(7.91e22, 5.87e33)), 7.91e22, 1e-12));
    }
}
