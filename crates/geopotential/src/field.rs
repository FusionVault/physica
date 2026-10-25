//! The one interface every potential field in the crate family shares: gravity from
//! [`Coefficients`], the main magnetic field (the `magnetics` crate), and the gravitomagnetic
//! field of a spinning body (`gravitics`).

use crate::Coefficients;

/// Geocentric radius (m) and latitude (radians) of a point at geodetic latitude `lat` (radians)
/// and height `height` (m) above the WGS84 ellipsoid: where to evaluate a [`Field`] for a point
/// given by GNSS or a map. [`NormalGravity::geocentric`](crate::normal::NormalGravity::geocentric)
/// does the same for any ellipsoid.
pub fn geocentric(lat: f64, height: f64) -> (f64, f64) {
    crate::normal::WGS84.geocentric(lat, height)
}

/// A static field with a scalar potential `V`, evaluated at geocentric radius `r` (m), latitude
/// `lat` and longitude `lon` (radians).
pub trait Field {
    /// The scalar potential `V` at the point, in the field's own units.
    fn potential(&self, r: f64, lat: f64, lon: f64) -> f64;

    /// The field vector in local (north, east, up) components, in the field's own units: `∇V`
    /// for gravity (m/s²), `−∇V` for a magnetic field (nT for the geomagnetic models). Undefined
    /// exactly at a pole, where east has no direction.
    fn vector(&self, r: f64, lat: f64, lon: f64) -> [f64; 3];

    /// The field's magnitude at the point.
    fn magnitude(&self, r: f64, lat: f64, lon: f64) -> f64 {
        let [n, e, u] = self.vector(r, lat, lon);
        (n * n + e * e + u * u).sqrt()
    }
}

impl Field for Coefficients {
    fn potential(&self, r: f64, lat: f64, lon: f64) -> f64 {
        self.partials(r, lat, lon).v
    }

    fn vector(&self, r: f64, lat: f64, lon: f64) -> [f64; 3] {
        self.partials(r, lat, lon).acceleration_neu(r, lat)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use metrology::{WGS84_A, WGS84_F};

    #[test]
    fn geodetic_to_geocentric() {
        let (r, lat) = geocentric(0.0, 0.0);
        assert!((r - WGS84_A).abs() < 1e-9 && lat == 0.0);
        let (r, lat) = geocentric(core::f64::consts::FRAC_PI_2, 0.0);
        assert!((r - WGS84_A * (1.0 - WGS84_F)).abs() < 1e-6 && (lat - core::f64::consts::FRAC_PI_2).abs() < 1e-12);
        // On the ellipsoid tan φ′ = (1 − e²) tan φ: 11.55′ less than geodetic at 45°.
        let (_, lat) = geocentric(45f64.to_radians(), 0.0);
        assert!((lat.tan() - (1.0 - WGS84_F * (2.0 - WGS84_F))).abs() < 1e-15);
        assert!(((45.0 - lat.to_degrees()) * 60.0 - 11.546).abs() < 1e-3);
    }

    #[test]
    fn a_point_mass_through_the_trait() {
        let mut earth = Coefficients::zeros(3.986_004_418e14, 6_378_137.0, 0);
        earth.set(0, 0, 1.0, 0.0);
        let r = 7_000e3;
        let field: &dyn Field = &earth;
        assert!((field.potential(r, 0.3, 1.0) - earth.gm / r).abs() < 1e-6);
        let [north, east, up] = field.vector(r, 0.3, 1.0);
        assert!(north.abs() < 1e-12 && east.abs() < 1e-12 && (up + earth.gm / (r * r)).abs() < 1e-12);
        assert!((field.magnitude(r, 0.3, 1.0) - earth.gm / (r * r)).abs() < 1e-12);
    }
}
