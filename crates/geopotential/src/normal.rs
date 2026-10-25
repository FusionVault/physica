//! The normal gravity field: the gravity of a level ellipsoid of revolution, the reference that
//! gravity anomalies, geoid heights and physical heights are measured against.
//!
//! Four constants define a level ellipsoid: its semi-major axis `a`, flattening `f`, geocentric
//! gravitational constant `GM` and rotation rate `ω`. Everything else follows in closed form
//! (Heiskanen & Moritz 1967, §2.7–2.10; [Moritz 2000](https://doi.org/10.1007/s001900050278)):
//! normal gravity at the equator and the pole, Somigliana's constant, the potential `U₀` on the
//! ellipsoid and the zonal harmonics `J₂ₙ` of the normal field. So [`GRS80`] and [`WGS84`] are
//! just their four defining constants, and the tests reproduce every derived constant their
//! defining documents print.
//!
//! ```
//! use geopotential::normal::{GRS80, WGS84};
//!
//! // Somigliana's formula: 9.7803 m/s² at the equator and 9.8322 at the pole.
//! assert!((GRS80.gamma(0.0) - 9.780_326_771_5).abs() < 1e-10);
//! assert!((GRS80.gamma(90f64.to_radians()) - 9.832_186_368_5).abs() < 1e-9);
//! // WGS84 differs from GRS80 by 0.14 mGal, through GM; a kilometre up, normal gravity is
//! // 308.5 mGal less.
//! assert!(((GRS80.gamma(0.0) - WGS84.gamma(0.0)) * 1e5 - 0.1436).abs() < 1e-3);
//! let lat = 45f64.to_radians();
//! assert!(((GRS80.gamma(lat) - GRS80.gamma_at(lat, 1e3)) * 1e5 - 308.5).abs() < 0.1);
//! ```

use core::f64::consts::FRAC_PI_4;

use metrology::{EARTH_ROTATION, GM_EARTH, GRS80_F, GRS80_GM, WGS84_A, WGS84_F};

use crate::Coefficients;

/// A level ellipsoid of revolution: an equipotential surface of its own gravity, which its four
/// defining constants fix completely.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NormalGravity {
    /// Semi-major axis `a`, m.
    pub a: f64,
    /// Flattening `f`.
    pub f: f64,
    /// Geocentric gravitational constant `GM`, m³/s².
    pub gm: f64,
    /// Angular velocity `ω`, rad/s.
    pub omega: f64,
}

/// The Geodetic Reference System 1980 (Moritz 2000), under most national datums.
pub const GRS80: NormalGravity = NormalGravity { a: WGS84_A, f: GRS80_F, gm: GRS80_GM, omega: EARTH_ROTATION };

/// WGS84 (NIMA TR8350.2), the GPS ellipsoid. It differs from GRS80 in the ninth digit of `f` and
/// in `GM`, which is the current IERS value.
pub const WGS84: NormalGravity = NormalGravity { a: WGS84_A, f: WGS84_F, gm: GM_EARTH, omega: EARTH_ROTATION };

impl NormalGravity {
    /// The semi-minor axis `b = a(1 − f)`, m.
    pub fn b(&self) -> f64 {
        self.a * (1.0 - self.f)
    }

    /// The first eccentricity squared, `e² = f(2 − f)`.
    pub fn e2(&self) -> f64 {
        self.f * (2.0 - self.f)
    }

    /// The second eccentricity `e′ = E/b`, with `E = √(a² − b²)` the linear eccentricity.
    fn e_prime(&self) -> f64 {
        (self.a * self.a - self.b() * self.b()).sqrt() / self.b()
    }

    /// `m = ω²a²b/GM`, the centrifugal-to-gravitational ratio at the equator.
    pub fn m(&self) -> f64 {
        self.omega * self.omega * self.a * self.a * self.b() / self.gm
    }

    /// Heiskanen & Moritz's `q₀` and `q₀′` of the second eccentricity (their Eqs. 2-58 and 2-67).
    fn q(&self) -> (f64, f64) {
        let e = self.e_prime();
        let at = e.atan();
        (0.5 * ((1.0 + 3.0 / (e * e)) * at - 3.0 / e), 3.0 * (1.0 + 1.0 / (e * e)) * (1.0 - at / e) - 1.0)
    }

    /// Normal gravity at the equator, `γ_e = (GM/ab)(1 − m − m e′q₀′/(6q₀))`, m/s².
    pub fn gamma_equator(&self) -> f64 {
        let (q0, q0p) = self.q();
        let m = self.m();
        self.gm / (self.a * self.b()) * (1.0 - m - m * self.e_prime() * q0p / (6.0 * q0))
    }

    /// Normal gravity at the poles, `γ_p = (GM/a²)(1 + m e′q₀′/(3q₀))`, m/s².
    pub fn gamma_pole(&self) -> f64 {
        let (q0, q0p) = self.q();
        self.gm / (self.a * self.a) * (1.0 + self.m() * self.e_prime() * q0p / (3.0 * q0))
    }

    /// Somigliana's constant `k = bγ_p/(aγ_e) − 1`.
    pub fn k(&self) -> f64 {
        self.b() * self.gamma_pole() / (self.a * self.gamma_equator()) - 1.0
    }

    /// The normal potential on the ellipsoid, `U₀ = (GM/E) arctan e′ + ω²a²/3`, m²/s². The geoid's
    /// `W₀` is compared with it in [`bruns`].
    pub fn u0(&self) -> f64 {
        let e = self.e_prime();
        self.gm / (e * self.b()) * e.atan() + self.omega * self.omega * self.a * self.a / 3.0
    }

    /// The dynamical form factor `J₂ = (e²/3)(1 − 2me′/(15q₀))`. GRS80 defined `J₂` and derived
    /// `f`; this is the inverse, and returns its 108 263×10⁻⁸ to 13 digits.
    pub fn j2(&self) -> f64 {
        let (q0, _) = self.q();
        self.e2() / 3.0 * (1.0 - 2.0 * self.m() * self.e_prime() / (15.0 * q0))
    }

    /// The zonal harmonic `J₂ₙ` of the normal potential, `(−1)ⁿ⁺¹ 3e²ⁿ/((2n+1)(2n+3)) ·
    /// (1 − n + 5nJ₂/e²)` (Heiskanen & Moritz, Eq. 2-92); `J₀ = −1` by the same formula is not
    /// meant, so `n ≥ 1`.
    pub fn zonal(&self, n: usize) -> f64 {
        let (e2, nf) = (self.e2(), n as f64);
        let sign = if n % 2 == 1 { 1.0 } else { -1.0 };
        sign * 3.0 * e2.powi(n as i32) / ((2.0 * nf + 1.0) * (2.0 * nf + 3.0)) * (1.0 - nf + 5.0 * nf * self.j2() / e2)
    }

    /// The normal potential as fully normalized coefficients to degree `nmax`: `C̄₂ₙ,₀ =
    /// −J₂ₙ/√(4n+1)`, nothing else. Synthesized with the centrifugal term
    /// ([`Partials::gravity_neu`](crate::Partials::gravity_neu)) it is Somigliana's gravity on the
    /// ellipsoid to 10⁻¹⁰ m/s² by degree 20, and the reference field a disturbing potential is
    /// taken from.
    pub fn coefficients(&self, nmax: usize) -> Coefficients {
        let mut c = Coefficients::zeros(self.gm, self.a, nmax);
        c.set(0, 0, 1.0, 0.0);
        for n in 1..=nmax / 2 {
            c.set(2 * n, 0, -self.zonal(n) / (4.0 * n as f64 + 1.0).sqrt(), 0.0);
        }
        c
    }

    /// Normal gravity on the ellipsoid at geodetic latitude `lat` (rad): Somigliana's closed form
    /// `γ = γ_e(1 + k sin²φ)/√(1 − e² sin²φ)`, m/s².
    pub fn gamma(&self, lat: f64) -> f64 {
        let s2 = lat.sin().powi(2);
        self.gamma_equator() * (1.0 + self.k() * s2) / (1.0 - self.e2() * s2).sqrt()
    }

    /// Normal gravity at geodetic latitude `lat` (rad) and ellipsoidal height `h` (m), to second
    /// order in `h/a` (Heiskanen & Moritz, Eq. 2-123): `γ(1 − 2(1 + f + m − 2f sin²φ)h/a +
    /// 3h²/a²)`. The free-air gradient, 0.3086 mGal/m, is its first term.
    pub fn gamma_at(&self, lat: f64, h: f64) -> f64 {
        let x = h / self.a;
        self.gamma(lat) * (1.0 - 2.0 * self.height_coefficient(lat) * x + 3.0 * x * x)
    }

    /// The mean of normal gravity between the ellipsoid and height `h` (m) along the normal,
    /// `γ(1 − (1 + f + m − 2f sin²φ)h/a + h²/a²)`: what a normal height divides by.
    pub fn mean_gamma(&self, lat: f64, h: f64) -> f64 {
        let x = h / self.a;
        self.gamma(lat) * (1.0 - self.height_coefficient(lat) * x + x * x)
    }

    fn height_coefficient(&self, lat: f64) -> f64 {
        1.0 + self.f + self.m() - 2.0 * self.f * lat.sin().powi(2)
    }

    /// Geocentric radius (m) and latitude (rad) of the point at geodetic latitude `lat` (rad) and
    /// height `h` (m) above this ellipsoid.
    pub fn geocentric(&self, lat: f64, h: f64) -> (f64, f64) {
        let e2 = self.e2();
        let (s, c) = lat.sin_cos();
        let n = self.a / (1.0 - e2 * s * s).sqrt();
        let (p, z) = ((n + h) * c, (n * (1.0 - e2) + h) * s);
        (p.hypot(z), z.atan2(p))
    }

    /// The dynamic height `C/γ₄₅` (m) of geopotential number `c = W₀ − W` (m²/s²): every point of
    /// one level surface gets the same number, which is why it is used along canals and lakes,
    /// though it is not a length.
    pub fn dynamic_height(&self, c: f64) -> f64 {
        c / self.gamma(FRAC_PI_4)
    }

    /// The normal height `C/γ̄` (m) of geopotential number `c` at latitude `lat`, with `γ̄` the
    /// [`mean_gamma`](Self::mean_gamma) up to that height, solved by iteration: the height above
    /// the quasigeoid, which needs no assumption about the rock below. The orthometric height
    /// divides by the mean *actual* gravity along the plumb line instead.
    pub fn normal_height(&self, c: f64, lat: f64) -> f64 {
        let mut h = c / self.gamma(lat);
        for _ in 0..8 {
            let next = c / self.mean_gamma(lat, h);
            if (next - h).abs() < 1e-12 {
                return next;
            }
            h = next;
        }
        h
    }
}

/// Bruns' formula with its zero-degree term: the height (m) of the level surface `W = w0` above
/// the ellipsoid of normal potential `u0`, from the disturbing potential `t = W − U` (m²/s²) and
/// normal gravity `gamma` there: `(T − (W₀ − U₀))/γ`. Leaving the second term out puts the geoid
/// about half a metre wrong; with [`metrology::W0_IHRS`] and GRS80's `U₀` it is +0.76 m.
pub fn bruns(t: f64, gamma: f64, w0: f64, u0: f64) -> f64 {
    (t - (w0 - u0)) / gamma
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Field;
    use metrology::W0_IHRS;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    #[test]
    fn grs80_derived_constants() {
        // Moritz (2000), "Geodetic Reference System 1980", derived constants to their digits.
        assert!(close(GRS80.e2(), 0.006_694_380_022_90, 5e-15) && close(GRS80.m(), 0.003_449_786_003_08, 5e-15));
        assert!(close(GRS80.gamma_equator(), 9.780_326_771_5, 5e-11) && close(GRS80.gamma_pole(), 9.832_186_368_5, 5e-11));
        assert!(close(GRS80.k(), 0.001_931_851_353, 5e-13) && close(GRS80.u0(), 62_636_860.850, 5e-4));
        // J₂ was GRS80's defining constant; it returns from the derived f to 13 digits.
        assert!(close(GRS80.j2(), 108_263e-8, 5e-16));
        assert!(close(GRS80.zonal(2), -0.000_002_370_912_22, 5e-15) && close(GRS80.zonal(3), 0.000_000_006_083_47, 5e-15));
        assert!(close(GRS80.zonal(4), -0.000_000_000_014_27, 5e-15));
        assert!(close(GRS80.gamma(FRAC_PI_4), 9.806_199_203, 5e-10));
        assert!(close(GRS80.gamma_pole(), GRS80.gamma(core::f64::consts::FRAC_PI_2), 1e-14));
    }

    #[test]
    fn wgs84_derived_constants() {
        // NIMA TR8350.2 (2000), Table 3.3: from its four defining constants, J₂ is derived.
        assert!(close(WGS84.e2(), 6.694_379_990_14e-3, 5e-15) && close(WGS84.m(), 0.003_449_786_506_84, 5e-15));
        assert!(close(WGS84.gamma_equator(), 9.780_325_335_9, 5e-11) && close(WGS84.gamma_pole(), 9.832_184_937_8, 1e-10));
        assert!(close(WGS84.k(), 0.001_931_852_652_41, 5e-15) && close(WGS84.u0(), 62_636_851.714_6, 5e-5));
        assert!(close(WGS84.j2(), 1.082_629_821_31e-3, 5e-15));
        // Its even zonals as the WGS84 document prints them (C̄₂₀ = −0.484166774985×10⁻³, C̄₄₀ =
        // 0.790303733511×10⁻⁶, C̄₆₀ = −0.168724961151×10⁻⁸, C̄₈₀ = 0.346052468394×10⁻¹¹); GRS80's
        // C̄₂₀, from its own J₂, is −0.484166854896×10⁻³. (geo5d's "GRS80 normal field" carried
        // WGS84's, with C̄₈₀ ten times too large.)
        let wgs84 = WGS84.coefficients(8);
        assert!(close(wgs84.c(2, 0), -0.484_166_774_985e-3, 5e-16) && close(wgs84.c(4, 0), 0.790_303_733_511e-6, 5e-18));
        assert!(close(wgs84.c(6, 0), -0.168_724_961_151e-8, 5e-20) && close(wgs84.c(8, 0), 0.346_052_468_394e-11, 5e-22));
        assert!(close(GRS80.coefficients(2).c(2, 0), -0.484_166_854_896e-3, 5e-16));
    }

    #[test]
    fn the_synthesized_normal_field_is_somiglianas() {
        // The even zonals to degree 20, with the centrifugal term, give Somigliana's gravity on
        // the ellipsoid to 10⁻¹⁰ m/s², and the potential U₀ to 10⁻⁵ m²/s², at every latitude.
        for normal in [GRS80, WGS84] {
            let field = normal.coefficients(20);
            for deg in [0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 89.9] {
                let lat = f64::to_radians(deg);
                let (r, psi) = normal.geocentric(lat, 0.0);
                let p = field.partials(r, psi, 0.3);
                let [n, e, u] = p.gravity_neu(r, psi, normal.omega);
                let g = (n * n + e * e + u * u).sqrt();
                assert!(close(g, normal.gamma(lat), 1e-10), "{deg}°: {g} against {}", normal.gamma(lat));
                assert!(close(p.geopotential(r, psi, normal.omega), normal.u0(), 1e-5), "{deg}°");
            }
            // Degree 8, as the normal fields are usually carried, is good to 2×10⁻⁸ m/s².
            let (r, psi) = normal.geocentric(0.7, 0.0);
            let [n, e, u] = normal.coefficients(8).partials(r, psi, 0.0).gravity_neu(r, psi, normal.omega);
            assert!(close((n * n + e * e + u * u).sqrt(), normal.gamma(0.7), 2e-8));
        }
    }

    #[test]
    fn gravity_above_the_ellipsoid() {
        // The free-air gradient: 0.3086 mGal/m, 0.30855 at 45° exactly, and 3h²/a² takes 0.07 mGal
        // back by 1 km. The mean to 1 km is Simpson's rule on γ(h), exact for a quadratic.
        let lat = FRAC_PI_4;
        let gradient = (GRS80.gamma(lat) - GRS80.gamma_at(lat, 1.0)) * 1e5;
        assert!(close(gradient, 0.308_55, 1e-5) && close(gradient, 0.3086, 1e-4), "{gradient}");
        assert!(close(GRS80.gamma(lat) - GRS80.gamma_at(lat, 1e3), 308.48e-5, 1e-7));
        let simpson = (GRS80.gamma(lat) + 4.0 * GRS80.gamma_at(lat, 500.0) + GRS80.gamma_at(lat, 1e3)) / 6.0;
        assert!(close(GRS80.mean_gamma(lat, 1e3), simpson, 1e-12));
        // Geodetic to geocentric: WGS84 through the crate's `geocentric`, and GRS80 within 1 mm.
        let (r, psi) = GRS80.geocentric(lat, 100.0);
        let (r2, psi2) = crate::geocentric(lat, 100.0);
        assert!(close(r, r2, 1e-3) && close(psi, psi2, 1e-10) && (45.0 - psi.to_degrees()) * 60.0 > 11.5);
    }

    #[test]
    fn heights_from_geopotential_numbers() {
        // C = 1000 γ₄₅ is 1000 m of dynamic height anywhere; the normal height is 0.157 m more
        // at 45°, because the mean normal gravity up to 1 km is 0.154 mGal·… less. Zero at the
        // quasigeoid.
        let c = 1000.0 * GRS80.gamma(FRAC_PI_4);
        assert!(close(GRS80.dynamic_height(c), 1000.0, 1e-9));
        let h = GRS80.normal_height(c, FRAC_PI_4);
        assert!(close(h * GRS80.mean_gamma(FRAC_PI_4, h), c, 1e-9) && close(h, 1000.157, 1e-3), "{h}");
        assert_eq!(GRS80.normal_height(0.0, 0.3), 0.0);
        // Bruns: the normal field's own geoid is the ellipsoid when W₀ = U₀; a W₀ 9.8 m²/s² lower
        // lies 1 m higher, so the IHRS W₀ puts the geoid 0.76 m above GRS80's U₀ surface.
        let gamma = GRS80.gamma(0.5);
        assert_eq!(bruns(0.0, gamma, GRS80.u0(), GRS80.u0()), 0.0);
        assert!(close(bruns(0.0, gamma, GRS80.u0() - 9.8, GRS80.u0()), 1.0, 0.002));
        assert!(close(bruns(0.0, GRS80.gamma(0.0), W0_IHRS, GRS80.u0()), 0.762, 1e-3));
    }

    #[test]
    fn deflections_of_the_vertical() {
        // A disturbing potential of one degree-3 zonal term: ξ is minus its meridian slope over
        // γr (central differences), 0.17″, and η is zero. A tesseral term gives η ≠ 0, minus its
        // slope along the parallel.
        let (lat, lon) = (0.5, 0.7);
        let (r, psi) = GRS80.geocentric(lat, 0.0);
        let gamma = GRS80.gamma(lat);
        let mut t = Coefficients::zeros(GRS80.gm, GRS80.a, 3);
        t.set(3, 0, 1e-6, 0.0);
        let (xi, eta) = t.partials(r, psi, lon).deflections(r, psi, gamma);
        let h = 1e-5;
        let slope = (t.potential(r, psi + h, lon) - t.potential(r, psi - h, lon)) / (2.0 * h);
        assert!(close(xi, -slope / (gamma * r), 1e-9 * xi.abs()) && xi.abs() > 1e-7 && eta == 0.0, "{xi} {eta}");
        t.set(3, 0, 0.0, 0.0);
        t.set(2, 1, 1e-6, 1e-6);
        let (_, eta) = t.partials(r, psi, lon).deflections(r, psi, gamma);
        let slope = (t.potential(r, psi, lon + h) - t.potential(r, psi, lon - h)) / (2.0 * h);
        assert!(close(eta, -slope / (gamma * r * psi.cos()), 1e-9 * eta.abs()) && eta.abs() > 1e-7, "{eta}");
    }
}
