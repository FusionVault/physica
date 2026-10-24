//! From instrument noise to gravity-field errors: formal errors of spherical-harmonic
//! coefficients for a gravity-gradiometer mission, and the spatial resolution of a degree.

/// The radial gravity gradient `V_zz` that a unit fully normalized coefficient of degree `n`
/// produces at radius `r`: `λ_n = GM/r³·(R/r)ⁿ·(n+1)(n+2)`, s⁻².
pub fn radial_gradient_transfer(n: usize, gm: f64, radius: f64, r: f64) -> f64 {
    gm / r.powi(3) * (radius / r).powi(n as i32) * ((n + 1) * (n + 2)) as f64
}

/// Formal error `σ_n` of every fully normalized coefficient of degree `n = 0..=nmax` (the same
/// for all orders) from a mission measuring the radial gradient `V_zz` at altitude `altitude` (m)
/// for `duration` seconds, with noise density `asd(f)` in s⁻² Hz⁻¹ᐟ²:
///
/// `σ_n = s(f_n) / (λ_n·√(2·duration))`, with `f_n = n / T_orbit`.
///
/// Least squares over `duration/Δt` samples of white noise with variance `s²/(2Δt)` and uniform
/// global coverage gives exactly this, independent of `Δt` (the test checks it against a
/// brute-force solution). Coloured noise is read at the degree's along-track frequency `n·f_orb`,
/// the usual heuristic. Not modelled: polar gaps, the other tensor components, attitude and
/// orbit errors, regularization.
pub fn gradiometry_errors(
    gm: f64,
    radius: f64,
    altitude: f64,
    duration: f64,
    nmax: usize,
    asd: impl Fn(f64) -> f64,
) -> Vec<f64> {
    let r = radius + altitude;
    let orbit = crate::missions::orbit_period(gm, radius, altitude);
    (0..=nmax)
        .map(|n| {
            let f = (n.max(1) as f64) / orbit;
            asd(f) / (radial_gradient_transfer(n, gm, radius, r) * (2.0 * duration).sqrt())
        })
        .collect()
}

/// Half-wavelength resolution of degree `n` on a sphere of radius `radius`, `πR/n` (m).
pub fn resolution(n: usize, radius: f64) -> f64 {
    core::f64::consts::PI * radius / n as f64
}

#[cfg(test)]
mod tests {
    use super::*;
    use geopotential::Legendre;

    #[test]
    fn resolution_of_zingerle_gradiometry_degrees() {
        // Zingerle et al. 2024: maximum degrees 223 and 263 correspond to ~90 and ~76 km.
        assert!((resolution(223, 6_378_137.0) / 1e3 - 89.9).abs() < 0.1);
        assert!((resolution(263, 6_378_137.0) / 1e3 - 76.2).abs() < 0.1);
    }

    /// Brute-force least squares on a Fibonacci-lattice sphere (near-uniform coverage) gives the
    /// same formal errors as the closed form, σ_n = σ_obs/(λ_n·√M).
    #[test]
    fn closed_form_matches_least_squares() {
        let (gm, radius, r, nmax, m_obs) = (3.986004415e14, 6378136.3, 6.6e6, 6usize, 20_000usize);
        let cols: Vec<(usize, usize, bool)> = (0..=nmax)
            .flat_map(|n| (0..=n).flat_map(move |m| [(n, m, false), (n, m, true)].into_iter().filter(move |c| !(c.1 == 0 && c.2))))
            .collect();
        let k = cols.len();
        let mut normal = vec![0.0; k * k];
        let mut leg = Legendre::new(nmax);
        let golden = core::f64::consts::PI * (3.0 - 5f64.sqrt());
        let mut row = vec![0.0; k];
        for i in 0..m_obs {
            let lat = (1.0 - 2.0 * (i as f64 + 0.5) / m_obs as f64).asin();
            let lon = golden * i as f64;
            leg.compute(lat);
            for (j, &(n, m, sine)) in cols.iter().enumerate() {
                let trig = if sine { (m as f64 * lon).sin() } else { (m as f64 * lon).cos() };
                row[j] = radial_gradient_transfer(n, gm, radius, r) * leg.p(n, m) * trig;
            }
            for a in 0..k {
                for b in 0..=a {
                    normal[a * k + b] += row[a] * row[b];
                }
            }
        }
        // Invert the symmetric normal matrix (Gauss–Jordan; it is small and well conditioned).
        let mut inv = vec![0.0; k * k];
        for a in 0..k {
            inv[a * k + a] = 1.0;
            for b in 0..a {
                normal[b * k + a] = normal[a * k + b];
            }
        }
        for p in 0..k {
            let d = normal[p * k + p];
            for c in 0..k {
                normal[p * k + c] /= d;
                inv[p * k + c] /= d;
            }
            for q in (0..k).filter(|&q| q != p) {
                let f = normal[q * k + p];
                for c in 0..k {
                    normal[q * k + c] -= f * normal[p * k + c];
                    inv[q * k + c] -= f * inv[p * k + c];
                }
            }
        }
        for (j, &(n, ..)) in cols.iter().enumerate() {
            let ls = inv[j * k + j].sqrt();
            let closed = 1.0 / (radial_gradient_transfer(n, gm, radius, r) * (m_obs as f64).sqrt());
            assert!((ls / closed - 1.0).abs() < 0.01, "{:?}: {ls:e} vs {closed:e}", cols[j]);
        }
    }

    #[test]
    fn errors_scale_with_noise_and_duration() {
        let gm = 3.986004415e14;
        let a = gradiometry_errors(gm, 6378136.3, 250e3, 1e7, 100, |_| 1e-11);
        let b = gradiometry_errors(gm, 6378136.3, 250e3, 4e7, 100, |_| 2e-11);
        assert!(a.iter().zip(&b).all(|(x, y)| (y / x - 1.0).abs() < 1e-12));
    }
}
