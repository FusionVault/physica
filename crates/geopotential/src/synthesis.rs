use crate::coefficients::Coefficients;
use crate::functional::Functional;
use crate::legendre::{Legendre, idx};

/// The gravitational potential and its first partial derivatives at one point.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Partials {
    /// `V`, m²/s².
    pub v: f64,
    /// `∂V/∂r`, m/s².
    pub dr: f64,
    /// `∂V/∂φ` (geocentric latitude), m²/s² per radian.
    pub dlat: f64,
    /// `∂V/∂λ`, m²/s² per radian.
    pub dlon: f64,
}

impl Partials {
    /// The gravitational acceleration `∇V` in local (north, east, up) components, m/s². Undefined
    /// exactly at a pole, where east has no direction.
    pub fn acceleration_neu(&self, r: f64, lat: f64) -> [f64; 3] {
        [self.dlat / r, self.dlon / (r * lat.cos()), self.dr]
    }
}

impl Coefficients {
    /// `V` and its partial derivatives at geocentric radius `r` (m), latitude `lat` and longitude
    /// `lon` (radians). Finite at the poles.
    ///
    /// # Panics
    /// If the degree exceeds [`MAX_DEGREE`](crate::MAX_DEGREE).
    // shortcut: builds a fresh Legendre table per call; reuse one when evaluating many points
    // (orbits) if this shows up in profiles.
    pub fn partials(&self, r: f64, lat: f64, lon: f64) -> Partials {
        let nmax = self.max_degree();
        let mut leg = Legendre::new(nmax);
        leg.compute_with_derivative(lat);
        let (cos_m, sin_m): (Vec<f64>, Vec<f64>) =
            (0..=nmax).map(|m| ((m as f64 * lon).cos(), (m as f64 * lon).sin())).unzip();
        let q = self.radius / r;
        let (mut v, mut vr, mut vlat, mut vlon) = (0.0, 0.0, 0.0, 0.0);
        let mut qn = 1.0;
        for n in 0..=nmax {
            let (mut sv, mut slat, mut slon) = (0.0, 0.0, 0.0);
            for m in 0..=n {
                let k = idx(n, m);
                let (c, s) = (self.c[k], self.s[k]);
                let even = c * cos_m[m] + s * sin_m[m];
                sv += leg.p[k] * even;
                slat += leg.dp[k] * even;
                slon += leg.p[k] * m as f64 * (s * cos_m[m] - c * sin_m[m]);
            }
            v += qn * sv;
            vr -= (n + 1) as f64 * qn * sv;
            vlat += qn * slat;
            vlon += qn * slon;
            qn *= q;
        }
        let gmr = self.gm / r;
        Partials { v: gmr * v, dr: gmr / r * vr, dlat: gmr * vlat, dlon: gmr * vlon }
    }

    /// A functional on the grid `lats × lons` (radians) at radius `r`, row-major by latitude:
    /// `out[i * lons.len() + j]` is at `(lats[i], lons[j])`.
    ///
    /// # Panics
    /// If the degree exceeds [`MAX_DEGREE`](crate::MAX_DEGREE).
    // shortcut: direct O(nmax · lons) longitude sum per row; switch to a real FFT (realfft) when
    // grids reach high degree and fine resolution.
    pub fn grid(&self, functional: &Functional, r: f64, lats: &[f64], lons: &[f64]) -> Vec<f64> {
        let nmax = self.max_degree();
        let factor: Vec<f64> =
            (0..=nmax).map(|n| functional.factor(n, self.gm, self.radius, r)).collect();
        let mut leg = Legendre::new(nmax);
        let (mut a, mut b) = (vec![0.0; nmax + 1], vec![0.0; nmax + 1]);
        let mut out = Vec::with_capacity(lats.len() * lons.len());
        for &lat in lats {
            leg.compute(lat);
            a.fill(0.0);
            b.fill(0.0);
            for (n, &f) in factor.iter().enumerate() {
                for m in 0..=n {
                    let k = idx(n, m);
                    let w = f * leg.p[k];
                    a[m] += w * self.c[k];
                    b[m] += w * self.s[k];
                }
            }
            out.extend(lons.iter().map(|&lon| {
                (0..=nmax)
                    .map(|m| {
                        let (s, c) = (m as f64 * lon).sin_cos();
                        a[m] * c + b[m] * s
                    })
                    .sum::<f64>()
            }));
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const GM: f64 = 3.986004415e14;
    const R: f64 = 6378136.3;

    /// A Kaula-like field (|C̄nm| ~ 1e-5/n²) with deterministic pseudo-random signs and sizes.
    fn field(nmax: usize) -> Coefficients {
        let mut seed: u64 = 0x2545_F491_4F6C_DD1D;
        let mut rand = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            (seed >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0
        };
        let mut f = Coefficients::zeros(GM, R, nmax);
        f.set(0, 0, 1.0, 0.0);
        for n in 2..=nmax {
            let kaula = 1e-5 / (n * n) as f64;
            for m in 0..=n {
                let s = if m == 0 { 0.0 } else { kaula * rand() };
                f.set(n, m, kaula * rand(), s);
            }
        }
        f
    }

    #[test]
    fn monopole_is_gm_over_r() {
        let mut f = Coefficients::zeros(GM, R, 0);
        f.set(0, 0, 1.0, 0.0);
        let r = 7.0e6;
        let p = f.partials(r, 0.4, 1.2);
        assert!((p.v - GM / r).abs() < 1e-8);
        assert!((p.dr + GM / (r * r)).abs() < 1e-14);
        assert_eq!((p.dlat, p.dlon), (0.0, 0.0));
    }

    #[test]
    fn partials_match_finite_differences() {
        let mut f = field(30);
        f.set(0, 0, 0.0, 0.0); // without GM/r the differences keep more significant digits
        let (r, lat, lon) = (6.9e6, 0.7, -2.1);
        let p = f.partials(r, lat, lon);
        let v = |r, lat, lon| f.partials(r, lat, lon).v;
        let (hr, ha) = (1.0, 1e-6);
        let dr = (v(r + hr, lat, lon) - v(r - hr, lat, lon)) / (2.0 * hr);
        let dlat = (v(r, lat + ha, lon) - v(r, lat - ha, lon)) / (2.0 * ha);
        let dlon = (v(r, lat, lon + ha) - v(r, lat, lon - ha)) / (2.0 * ha);
        assert!((p.dr - dr).abs() < 1e-8 * p.dr.abs(), "{} vs {dr}", p.dr);
        assert!((p.dlat - dlat).abs() < 1e-5 * p.dlat.abs(), "{} vs {dlat}", p.dlat);
        assert!((p.dlon - dlon).abs() < 1e-5 * p.dlon.abs(), "{} vs {dlon}", p.dlon);
    }

    /// Outside the masses V is harmonic, so the trace of its Hessian (here by central differences
    /// of the analytic gradient in Cartesian coordinates) vanishes.
    #[test]
    fn laplace_trace_vanishes() {
        let f = field(40);
        let grad = |x: [f64; 3]| {
            let r = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
            let (lat, lon) = ((x[2] / r).asin(), x[1].atan2(x[0]));
            let [n, e, u] = f.partials(r, lat, lon).acceleration_neu(r, lat);
            let ((sp, cp), (sl, cl)) = (lat.sin_cos(), lon.sin_cos());
            [
                -sp * cl * n - sl * e + cp * cl * u,
                -sp * sl * n + cl * e + cp * sl * u,
                cp * n + sp * u,
            ]
        };
        let x0 = [3.1e6, -4.4e6, 4.0e6];
        let h = 1.0;
        let mut trace = 0.0;
        let mut scale = 0.0f64;
        for i in 0..3 {
            let (mut hi, mut lo) = (x0, x0);
            hi[i] += h;
            lo[i] -= h;
            let d = (grad(hi)[i] - grad(lo)[i]) / (2.0 * h);
            trace += d;
            scale = scale.max(d.abs());
        }
        assert!(trace.abs() < 1e-7 * scale, "trace {trace:e} vs diagonal {scale:e}");
    }

    /// The grid functionals equal the matching combinations of the point partials.
    #[test]
    fn grid_functionals_match_partials() {
        let mut f = field(25);
        f.set(0, 0, 0.0, 0.0); // a disturbing potential has no monopole
        let r = 6.6e6;
        let (lats, lons) = ([-1.2, 0.1, 1.5], [0.0, 2.0, 4.5, 6.2]);
        let grid = |func| f.grid(&func, r, &lats, &lons);
        let (pot, geoid, dist, anom) = (
            grid(Functional::Potential),
            grid(Functional::GeoidHeight),
            grid(Functional::GravityDisturbance),
            grid(Functional::GravityAnomaly),
        );
        // Tolerances are relative to the largest value, since single points can cancel to ~0.
        let max_abs = |v: &[f64]| v.iter().fold(0.0f64, |a, x| a.max(x.abs()));
        let close = |a: f64, b: f64, grid: &[f64]| (a - b).abs() <= 1e-12 * max_abs(grid);
        for (i, &lat) in lats.iter().enumerate() {
            for (j, &lon) in lons.iter().enumerate() {
                let k = i * lons.len() + j;
                let p = f.partials(r, lat, lon);
                assert!(close(pot[k], p.v, &pot));
                assert!(close(geoid[k], p.v / (GM / (r * r)), &geoid));
                assert!(close(dist[k], -p.dr, &dist));
                assert!(close(anom[k], -p.dr - 2.0 * p.v / r, &anom));
            }
        }
        // Radial gradient against a finite difference of the disturbance in r.
        let h = 1.0;
        let up = f.grid(&Functional::GravityDisturbance, r + h, &lats, &lons);
        let down = f.grid(&Functional::GravityDisturbance, r - h, &lats, &lons);
        let trr = grid(Functional::RadialGradient);
        for k in 0..trr.len() {
            let fd = -(up[k] - down[k]) / (2.0 * h);
            assert!((trr[k] - fd).abs() < 1e-6 * max_abs(&trr), "{} vs {fd}", trr[k]);
        }
    }

    #[test]
    fn rescaling_leaves_the_field_unchanged() {
        let f = field(20);
        let g = f.rescaled(3.986004418e14, 6378137.0);
        let (a, b) = (f.partials(7.0e6, 0.3, 1.1), g.partials(7.0e6, 0.3, 1.1));
        assert!((a.v - b.v).abs() < 1e-14 * a.v.abs() && (a.dlat - b.dlat).abs() < 1e-12 * a.dlat.abs());
    }

    #[test]
    fn water_height_degree_two() {
        // A GRACE-sized ΔC̄20 of 1e-10 is about 19 mm of water at the pole (k'_2 = −0.303,
        // P̄20(pole) = √5).
        let love = [0.0, 0.027, -0.303];
        let mut d = Coefficients::zeros(GM, R, 2);
        d.set(2, 0, 1e-10, 0.0);
        let ewh = d.grid(&Functional::WaterHeight { love: &love }, R, &[core::f64::consts::FRAC_PI_2], &[0.0])[0];
        let expect = R * 5517.0 / 3000.0 * 5.0 / 0.697 * 1e-10 * 5f64.sqrt();
        assert!((ewh - expect).abs() < 1e-15 && (ewh - 0.0188).abs() < 1e-4, "{ewh}");
    }
}
