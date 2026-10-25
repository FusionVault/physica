use crate::coefficients::Coefficients;
use crate::functional::Functional;
use crate::legendre::{Legendre, gauss_legendre, idx};

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

    /// The geopotential `W = V + ½ω²r²cos²φ` (m²/s²) of a body rotating at `omega` (rad/s): the
    /// potential a clock on it senses and a level surface is a surface of.
    pub fn geopotential(&self, r: f64, lat: f64, omega: f64) -> f64 {
        self.v + 0.5 * (omega * r * lat.cos()).powi(2)
    }

    /// Gravity `∇W`, the acceleration plus the centrifugal term `ω²r cos φ` directed away from the
    /// axis, in local (north, east, up) components, m/s².
    pub fn gravity_neu(&self, r: f64, lat: f64, omega: f64) -> [f64; 3] {
        let [n, e, u] = self.acceleration_neu(r, lat);
        let (s, c) = lat.sin_cos();
        let w2r = omega * omega * r;
        [n - w2r * c * s, e, u + w2r * c * c]
    }

    /// The deflections of the vertical `(ξ, η)` (north and east, rad) from the partials of a
    /// disturbing potential `T` and normal gravity `gamma` there, in the spherical approximation
    /// (Heiskanen & Moritz, Eq. 2-204): `ξ = −(∂T/∂φ)/(γr)`, `η = −(∂T/∂λ)/(γr cos φ)`. Multiply by
    /// 206 265 for arcseconds.
    pub fn deflections(&self, r: f64, lat: f64, gamma: f64) -> (f64, f64) {
        (-self.dlat / (gamma * r), -self.dlon / (gamma * r * lat.cos()))
    }
}

/// The potential with its first and second partial derivatives in `(r, φ, λ)` at one point: what
/// a gravity-gradient tensor is built from.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SecondPartials {
    /// `V`, `∂V/∂r`, `∂V/∂φ`, `∂V/∂λ`.
    pub first: Partials,
    /// `∂²V/∂r²`.
    pub drr: f64,
    /// `∂²V/∂φ²`.
    pub dlatlat: f64,
    /// `∂²V/∂λ²`.
    pub dlonlon: f64,
    /// `∂²V/∂r∂φ`.
    pub drlat: f64,
    /// `∂²V/∂r∂λ`.
    pub drlon: f64,
    /// `∂²V/∂φ∂λ`.
    pub dlatlon: f64,
}

impl SecondPartials {
    /// The gravity-gradient (tidal) tensor `Γᵢⱼ = ∂²V/∂xᵢ∂xⱼ` in the local frame at geocentric
    /// radius `r` (m) and latitude `lat` (radians), s⁻², as the symmetric matrix
    /// `[[Γnn, Γne, Γnu], [Γne, Γee, Γeu], [Γnu, Γeu, Γuu]]` with x = north, y = east, z = up.
    /// Multiply by 1e9 for Eötvös. Its trace vanishes in free space (Laplace's equation), and a
    /// gradiometer measures it directly: a rigid vehicle's own acceleration is common to both
    /// proof masses and cancels, leaving only this.
    ///
    /// For a point mass it is `diag(−GM/r³, −GM/r³, +2GM/r³)`. Undefined exactly at a pole.
    pub fn tensor_neu(&self, r: f64, lat: f64) -> [[f64; 3]; 3] {
        let (sec, tan) = (1.0 / lat.cos(), lat.tan());
        let (vr, vlat, vlon) = (self.first.dr, self.first.dlat, self.first.dlon);
        let nn = vr / r + self.dlatlat / (r * r);
        let ee = vr / r - tan * vlat / (r * r) + self.dlonlon * sec * sec / (r * r);
        let uu = self.drr;
        let nu = (self.drlat - vlat / r) / r;
        let eu = (self.drlon - vlon / r) * sec / r;
        let ne = (self.dlatlon + tan * vlon) * sec / (r * r);
        [[nn, ne, nu], [ne, ee, eu], [nu, eu, uu]]
    }
}

impl Coefficients {
    /// `V` with its first and second partial derivatives at geocentric radius `r` (m), latitude
    /// `lat` and longitude `lon` (radians). Finite at the poles.
    ///
    /// # Panics
    /// If the degree exceeds [`MAX_DEGREE`](crate::MAX_DEGREE).
    pub fn second_partials(&self, r: f64, lat: f64, lon: f64) -> SecondPartials {
        let nmax = self.max_degree();
        let mut leg = Legendre::new(nmax);
        leg.compute_with_second_derivative(lat);
        let (cos_m, sin_m): (Vec<f64>, Vec<f64>) =
            (0..=nmax).map(|m| ((m as f64 * lon).cos(), (m as f64 * lon).sin())).unzip();
        let q = self.radius / r;
        let mut acc = [0.0f64; 10]; // v, r, φ, λ, rr, φφ, λλ, rφ, rλ, φλ
        let mut qn = 1.0;
        for n in 0..=nmax {
            let mut sum = [0.0f64; 6]; // P·e, P'·e, P''·e, P·de, P·d²e, P'·de
            for m in 0..=n {
                let k = idx(n, m);
                let (c, s) = (self.c[k], self.s[k]);
                let mf = m as f64;
                let even = c * cos_m[m] + s * sin_m[m];
                let deven = mf * (s * cos_m[m] - c * sin_m[m]);
                sum[0] += leg.p[k] * even;
                sum[1] += leg.dp[k] * even;
                sum[2] += leg.ddp[k] * even;
                sum[3] += leg.p[k] * deven;
                sum[4] -= leg.p[k] * mf * mf * even;
                sum[5] += leg.dp[k] * deven;
            }
            let (n1, n12) = ((n + 1) as f64, ((n + 1) * (n + 2)) as f64);
            let terms = [sum[0], -n1 * sum[0], sum[1], sum[3], n12 * sum[0], sum[2], sum[4], -n1 * sum[1], -n1 * sum[3], sum[5]];
            for (a, t) in acc.iter_mut().zip(terms) {
                *a += qn * t;
            }
            qn *= q;
        }
        let gmr = self.gm / r;
        SecondPartials {
            first: Partials { v: gmr * acc[0], dr: gmr / r * acc[1], dlat: gmr * acc[2], dlon: gmr * acc[3] },
            drr: gmr / (r * r) * acc[4],
            dlatlat: gmr * acc[5],
            dlonlon: gmr * acc[6],
            drlat: gmr / r * acc[7],
            drlon: gmr / r * acc[8],
            dlatlon: gmr * acc[9],
        }
    }

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

/// The grid that [`Coefficients::analysis`] reads to degree `nmax`: `nmax + 1` Gauss–Legendre
/// latitudes (radians, south to north) and `2·nmax + 2` longitudes `2πj/(2·nmax + 2)`.
pub fn gauss_grid(nmax: usize) -> (Vec<f64>, Vec<f64>) {
    let lats = gauss_legendre(nmax + 1).0.iter().map(|x| x.asin()).collect();
    let nlon = 2 * nmax + 2;
    (lats, (0..nlon).map(|j| core::f64::consts::TAU * j as f64 / nlon as f64).collect())
}

impl Coefficients {
    /// The coefficients to degree `nmax` of a functional sampled on [`gauss_grid`]`(nmax)` at
    /// radius `r`, row-major by latitude as [`grid`](Self::grid) writes it: its exact inverse
    /// for a field of that degree. The latitude integral is Gauss–Legendre quadrature, exact for
    /// the polynomial `P̄nm P̄n′m` of degree `n + n′ ≤ 2·nmax`; the longitude sum is exact for
    /// orders up to `nmax`. A degree whose factor is zero, such as degree 1 of a gravity anomaly,
    /// is left at zero.
    ///
    /// # Panics
    /// If `values` is not one value per grid point.
    pub fn analysis(functional: &Functional, r: f64, values: &[f64], nmax: usize, gm: f64, radius: f64) -> Coefficients {
        let (x, w) = gauss_legendre(nmax + 1);
        let nlon = 2 * nmax + 2;
        assert_eq!(values.len(), x.len() * nlon, "one value per point of gauss_grid(nmax)");
        let trig: Vec<(f64, f64)> =
            (0..=nmax).flat_map(|m| (0..nlon).map(move |j| (m as f64 * core::f64::consts::TAU * j as f64 / nlon as f64).sin_cos())).collect();
        let mut out = Coefficients::zeros(gm, radius, nmax);
        let mut leg = Legendre::new(nmax);
        for (i, (&xi, &wi)) in x.iter().zip(&w).enumerate() {
            leg.compute(xi.asin());
            let row = &values[i * nlon..(i + 1) * nlon];
            for m in 0..=nmax {
                let (mut a, mut b) = (0.0, 0.0);
                for (f, &(sin, cos)) in row.iter().zip(&trig[m * nlon..(m + 1) * nlon]) {
                    a += f * cos;
                    b += f * sin;
                }
                // (1/4π)∫∫ f P̄nm (cos, sin) dΩ, with dλ = 2π/nlon.
                let k = wi / (2.0 * nlon as f64);
                for n in m..=nmax {
                    let p = k * leg.p[idx(n, m)];
                    out.c[idx(n, m)] += p * a;
                    out.s[idx(n, m)] += p * b;
                }
            }
        }
        for n in 0..=nmax {
            let factor = functional.factor(n, gm, radius, r);
            for m in 0..=n {
                let (c, s) = (out.c[idx(n, m)], out.s[idx(n, m)]);
                let (c, s) = if factor == 0.0 { (0.0, 0.0) } else { (c / factor, s / factor) };
                out.set(n, m, c, s);
            }
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
    fn second_partials_match_differences_and_laplace() {
        let mut f = field(25);
        f.set(0, 0, 0.0, 0.0);
        for (r, lat, lon) in [(6.9e6, 0.6, -1.1), (7.2e6, -1.3, 2.9), (6.7e6, 1.569, 0.4)] {
            let sp = f.second_partials(r, lat, lon);
            let (a, b) = (sp.first, f.partials(r, lat, lon));
            for (x, y) in [(a.v, b.v), (a.dr, b.dr), (a.dlat, b.dlat), (a.dlon, b.dlon)] {
                assert!((x - y).abs() <= 1e-13 * y.abs().max(1e-30), "{x} vs {y}"); // same sums, other order
            }
            let (hr, ha) = (1.0, 1e-6);
            let p = |r, lat, lon| f.partials(r, lat, lon);
            let close = |a: f64, b: f64| (a - b).abs() <= 1e-5 * (a.abs() + b.abs() + 1e-30);
            assert!(close(sp.drr, (p(r + hr, lat, lon).dr - p(r - hr, lat, lon).dr) / (2.0 * hr)));
            assert!(close(sp.dlatlat, (p(r, lat + ha, lon).dlat - p(r, lat - ha, lon).dlat) / (2.0 * ha)));
            assert!(close(sp.dlonlon, (p(r, lat, lon + ha).dlon - p(r, lat, lon - ha).dlon) / (2.0 * ha)));
            assert!(close(sp.drlat, (p(r, lat + ha, lon).dr - p(r, lat - ha, lon).dr) / (2.0 * ha)));
            assert!(close(sp.drlon, (p(r, lat, lon + ha).dr - p(r, lat, lon - ha).dr) / (2.0 * ha)));
            assert!(close(sp.dlatlon, (p(r, lat, lon + ha).dlat - p(r, lat, lon - ha).dlat) / (2.0 * ha)));
            // Laplace in spherical coordinates: V_rr + 2V_r/r + (V_φφ − tanφ·V_φ + V_λλ/cos²φ)/r² = 0.
            let lap = sp.drr + 2.0 * sp.first.dr / r
                + (sp.dlatlat - lat.tan() * sp.first.dlat + sp.dlonlon / lat.cos().powi(2)) / (r * r);
            assert!(lap.abs() < 1e-9 * sp.drr.abs(), "Laplacian {lap:e} vs {:e}", sp.drr);
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

    /// The gradient tensor against finite differences of the potential in a true local Cartesian
    /// frame, and against the point-mass form.
    #[test]
    fn gradient_tensor_in_the_local_frame() {
        let (gm, radius) = (3.986_004_415e14, 6_378_136.3);
        let mut field = Coefficients::zeros(gm, radius, 4);
        field.set(0, 0, 1.0, 0.0);
        field.set(2, 0, -4.841_695e-4, 0.0);
        field.set(2, 2, 2.439_37e-6, -1.400_31e-6);
        field.set(3, 1, 2.030_45e-6, 2.482_42e-7);
        field.set(4, 4, -1.880_28e-7, 3.088_01e-7);

        let (r, lat, lon) = (7_000e3, 0.4, 1.1);
        let got = field.second_partials(r, lat, lon).tensor_neu(r, lat);

        // Displace in the local frame, convert back to spherical, and evaluate V there.
        let (sp, cp) = lat.sin_cos();
        let (sl, cl) = lon.sin_cos();
        let p0 = [r * cp * cl, r * cp * sl, r * sp];
        let axes = [[-sp * cl, -sp * sl, cp], [-sl, cl, 0.0], [cp * cl, cp * sl, sp]];
        let v_at = |d: [f64; 3]| {
            let p: Vec<f64> = (0..3).map(|k| p0[k] + (0..3).map(|a| d[a] * axes[a][k]).sum::<f64>()).collect();
            let rr = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
            field.partials(rr, (p[2] / rr).asin(), p[1].atan2(p[0])).v
        };
        // h balances round-off (∝1/h², from differencing a 10⁷-sized potential) against
        // truncation (∝h²); the residual noise floor is about 1e-14 s⁻².
        let h = 400.0;
        for i in 0..3 {
            for j in 0..3 {
                let step = |a: f64, b: f64| { let mut d = [0.0; 3]; d[i] += a; d[j] += b; v_at(d) };
                let fd = (step(h, h) - step(h, -h) - step(-h, h) + step(-h, -h)) / (4.0 * h * h);
                assert!((got[i][j] - fd).abs() < 1e-13, "({i},{j}): {} vs {fd}", got[i][j]);
            }
        }
        // Laplace: the trace vanishes.
        assert!((0..3).map(|k| got[k][k]).sum::<f64>().abs() < 1e-16);

        // A point mass: diag(-GM/r^3, -GM/r^3, +2GM/r^3), with no off-diagonal terms.
        let mut mono = Coefficients::zeros(gm, radius, 0);
        mono.set(0, 0, 1.0, 0.0);
        let t = mono.second_partials(r, lat, lon).tensor_neu(r, lat);
        let g3 = gm / r.powi(3);
        assert!((t[0][0] + g3).abs() < 1e-18 && (t[1][1] + g3).abs() < 1e-18 && (t[2][2] - 2.0 * g3).abs() < 1e-18);
        assert!(t[0][1].abs() < 1e-20 && t[0][2].abs() < 1e-20 && t[1][2].abs() < 1e-20);
        // 1 Eotvos is 1e-9 s^-2: the radial gradient 2GM/r^3 at 7000 km is about 2330 E.
        assert!((2.0 * g3 * 1e9 - 2325.0).abs() < 5.0);
    }

    #[test]
    fn analysis_inverts_the_grid() {
        // Unit-size random coefficients to degree 16, the notes' case: synthesis on the Gauss grid
        // and back recovers them to 9×10⁻¹⁵ (the notes' 1.6×10⁻¹⁴). A gravity disturbance at 400 km
        // round-trips too.
        let nmax = 16;
        let mut seed: u64 = 7;
        let mut rand = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            (seed >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0
        };
        let mut f = Coefficients::zeros(GM, R, nmax);
        for n in 0..=nmax {
            for m in 0..=n {
                let s = if m == 0 { 0.0 } else { rand() };
                f.set(n, m, rand(), s);
            }
        }
        let (lats, lons) = gauss_grid(nmax);
        for (functional, r, tol) in [(Functional::Potential, R, 2e-14), (Functional::GravityDisturbance, R + 400e3, 2e-14)] {
            let values = f.grid(&functional, r, &lats, &lons);
            let back = Coefficients::analysis(&functional, r, &values, nmax, GM, R);
            let err = (0..=nmax).flat_map(|n| (0..=n).map(move |m| (n, m))).map(|(n, m)| (back.c(n, m) - f.c(n, m)).abs().max((back.s(n, m) - f.s(n, m)).abs())).fold(0.0, f64::max);
            assert!(err < tol, "{err:e}");
        }
    }
}
