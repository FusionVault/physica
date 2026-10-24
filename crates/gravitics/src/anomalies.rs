//! Gravity's open anomalies, each set against the measurement that tests it: the scatter of the
//! measured Newtonian constant, the Earth-flyby speed changes, gravitational shielding and the
//! Pioneer acceleration.
//!
//! The pattern is that single-apparatus anomalies shrink when a second laboratory rebuilds the
//! apparatus. The sixteen G measurements behind CODATA 2022 span 551 ppm and disagree with their
//! own error bars by a Birge ratio of 3.6, so CODATA expands the uncertainty 3.9-fold
//! ([CODATA 2022](https://arxiv.org/abs/2409.03787)). Anderson's flyby formula fitted six flybys
//! and then predicted +6 mm/s for Juno, which showed 0 ± 0.8 mm/s
//! ([Wikipedia: Flyby anomaly](https://en.wikipedia.org/wiki/Flyby_anomaly)). Majorana's
//! shielding coefficient would hide 4.6% of the Sun's pull behind the Earth; lunar laser ranging
//! bounds it ten orders lower. Pioneer's acceleration decayed with the spacecraft's plutonium,
//! which identifies it as thermal recoil ([Turyshev et al. 2012](https://arxiv.org/abs/1204.2507)).
//!
//! ```
//! use gravitics::anomalies::{CODATA_2022_G, CODATA_2022_G_CORRELATIONS, birge_ratio, gls_mean};
//!
//! // Rebuild CODATA 2022's G from its sixteen inputs (units of 10⁻¹¹ m³ kg⁻¹ s⁻²).
//! let n = CODATA_2022_G.len();
//! let mut cov = vec![vec![0.0; n]; n];
//! for (i, row) in CODATA_2022_G.iter().enumerate() {
//!     cov[i][i] = row.2 * row.2;
//! }
//! let at = |id: &str| CODATA_2022_G.iter().position(|r| r.0 == id).unwrap();
//! for &(a, b, r) in &CODATA_2022_G_CORRELATIONS {
//!     let (i, j) = (at(a), at(b));
//!     cov[i][j] = r * CODATA_2022_G[i].2 * CODATA_2022_G[j].2;
//!     cov[j][i] = cov[i][j];
//! }
//! let values: Vec<f64> = CODATA_2022_G.iter().map(|r| r.1).collect();
//! let (mean, _, chi2) = gls_mean(&values, &cov).unwrap();
//! assert!((mean - 6.674_30).abs() < 1e-6);
//! assert!((birge_ratio(chi2, n - 1) - 3.61).abs() < 0.01);
//! ```

/// The sixteen measurements of G in the CODATA 2022 adjustment: identifier, value and standard
/// uncertainty, in 10⁻¹¹ m³ kg⁻¹ s⁻² ([CODATA 2022](https://arxiv.org/abs/2409.03787)).
pub const CODATA_2022_G: [(&str, f64, f64); 16] = [
    ("NIST-82", 6.672_48, 0.000_43),
    ("TR&D-96", 6.672_9, 0.000_5),
    ("LANL-97", 6.673_98, 0.000_70),
    ("UWash-00", 6.674_255, 0.000_092),
    ("BIPM-01", 6.675_59, 0.000_27),
    ("UWup-02", 6.674_22, 0.000_98),
    ("MSL-03", 6.673_87, 0.000_27),
    ("HUST-05", 6.672_22, 0.000_87),
    ("UZur-06", 6.674_25, 0.000_12),
    ("HUST-09", 6.673_49, 0.000_18),
    ("BIPM-14", 6.675_54, 0.000_16),
    ("LENS-14", 6.671_91, 0.000_99),
    ("UCI-14", 6.674_35, 0.000_13),
    ("HUSTT-18", 6.674_184, 0.000_078),
    ("HUSTA-18", 6.674_484, 0.000_077),
    ("JILA-18", 6.672_60, 0.000_25),
];

/// The correlation coefficients between pairs of [`CODATA_2022_G`] inputs.
pub const CODATA_2022_G_CORRELATIONS: [(&str, &str, f64); 3] =
    [("NIST-82", "LANL-97", 0.351), ("HUST-05", "HUST-09", 0.134), ("HUST-09", "HUSTT-18", 0.068)];

/// Majorana's 1920 gravitational absorption coefficient, 6.7×10⁻¹² cm²/g, in m²/kg.
pub const MAJORANA_ABSORPTION: f64 = 6.7e-13;
/// Lunar laser ranging's bound on gravitational absorption (Eckhardt 1990), m²/kg.
pub const LLR_ABSORPTION_BOUND: f64 = 1e-22;

/// Generalised-least-squares mean of correlated `values` with covariance matrix `cov`:
/// `(mean, uncertainty, χ²)`, with `Ḡ = 1ᵀC⁻¹y/1ᵀC⁻¹1`, `u = (1ᵀC⁻¹1)^(−½)` and `χ² = rᵀC⁻¹r`.
/// `None` if `cov` is not positive definite or the sizes differ.
pub fn gls_mean(values: &[f64], cov: &[Vec<f64>]) -> Option<(f64, f64, f64)> {
    let n = values.len();
    if n == 0 || cov.len() != n || cov.iter().any(|row| row.len() != n) {
        return None;
    }
    // Cholesky factor L of C, then C⁻¹b by forward and back substitution.
    let mut l = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..=i {
            let s = cov[i][j] - (0..j).map(|k| l[i][k] * l[j][k]).sum::<f64>();
            if i == j {
                if s <= 0.0 {
                    return None;
                }
                l[i][i] = s.sqrt();
            } else {
                l[i][j] = s / l[j][j];
            }
        }
    }
    let solve = |b: &[f64]| {
        let mut y = vec![0.0; n];
        for i in 0..n {
            y[i] = (b[i] - (0..i).map(|k| l[i][k] * y[k]).sum::<f64>()) / l[i][i];
        }
        for i in (0..n).rev() {
            y[i] = (y[i] - (i + 1..n).map(|k| l[k][i] * y[k]).sum::<f64>()) / l[i][i];
        }
        y
    };
    let w = solve(&vec![1.0; n]);
    let total: f64 = w.iter().sum();
    let mean = w.iter().zip(values).map(|(w, y)| w * y).sum::<f64>() / total;
    let r: Vec<f64> = values.iter().map(|y| y - mean).collect();
    let chi2 = r.iter().zip(solve(&r)).map(|(a, b)| a * b).sum();
    Some((mean, total.sqrt().recip(), chi2))
}

/// The Birge ratio `√(χ²/ν)`: how far a set of measurements disagrees with its stated errors (1
/// means they agree). CODATA expands the uncertainty by a factor k (3.9 for G in 2022), which
/// divides χ² by k².
pub fn birge_ratio(chi2: f64, dof: usize) -> f64 {
    (chi2 / dof as f64).sqrt()
}

/// Anderson et al.'s 2008 flyby formula, `ΔV = K V∞ (cos δᵢ − cos δₒ)` with
/// `K = 2ω⊕R⊕/c = 3.099×10⁻⁶`, for asymptotic speed `v_inf` (m/s) and incoming and outgoing
/// declinations (radians). It reproduces the six fitted flybys and failed out of sample: Juno was
/// predicted at +6 mm/s and observed at 0 ± 0.8 mm/s. A labelled reproduction, not physics.
pub fn flyby_anderson(v_inf: f64, dec_in: f64, dec_out: f64) -> f64 {
    const K: f64 = 2.0 * 7.292_115e-5 * 6_371e3 / metrology::C;
    K * v_inf * (dec_in.cos() - dec_out.cos())
}

/// Column density `2Rρ` (kg/m²) through the centre of a body of `radius` (m) and mean `density`.
pub fn column_density(radius: f64, density: f64) -> f64 {
    2.0 * radius * density
}

/// The fraction of gravity a column `column` (kg/m²) would absorb with absorption coefficient
/// `h` (m²/kg): `1 − e^(−hΣ)`. With [`MAJORANA_ABSORPTION`] the Earth would hide 4.6% of the Sun's
/// pull at eclipse; [`LLR_ABSORPTION_BOUND`] caps it at 7×10⁻¹².
pub fn shielding_loss(h: f64, column: f64) -> f64 {
    -(-h * column).exp_m1()
}

/// The Pioneer anomaly as thermal recoil: an acceleration `a0` at the start that halves every
/// `tau`, after time `t` (any unit, the same for both): `a₀ 2^(−t/τ)`.
pub fn pioneer_thermal(a0: f64, tau: f64, t: f64) -> f64 {
    a0 * (-t / tau).exp2()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, rel: f64) -> bool {
        (a / b - 1.0).abs() < rel
    }

    fn codata() -> (Vec<f64>, Vec<Vec<f64>>) {
        let n = CODATA_2022_G.len();
        let mut cov = vec![vec![0.0; n]; n];
        for (i, row) in CODATA_2022_G.iter().enumerate() {
            cov[i][i] = row.2 * row.2;
        }
        let at = |id: &str| CODATA_2022_G.iter().position(|r| r.0 == id).unwrap();
        for &(a, b, r) in &CODATA_2022_G_CORRELATIONS {
            let (i, j) = (at(a), at(b));
            cov[i][j] = r * CODATA_2022_G[i].2 * CODATA_2022_G[j].2;
            cov[j][i] = cov[i][j];
        }
        (CODATA_2022_G.iter().map(|r| r.1).collect(), cov)
    }

    #[test]
    fn codata_g() {
        let (y, cov) = codata();
        let (mean, u, chi2) = gls_mean(&y, &cov).unwrap();
        assert!((mean - 6.674_300_4).abs() < 1e-7 && close(u, 3.789_13e-5, 1e-5) && close(chi2, 195.679, 1e-5));
        assert!(close(birge_ratio(chi2, 15), 3.611_82, 1e-5));
        // Expanded by 3.9: 22 ppm and χ² of 12.9.
        assert!(close(3.9 * u / mean * 1e6, 22.14, 1e-3) && close(chi2 / 3.9f64.powi(2), 12.865, 1e-3));
        for (id, z) in [("BIPM-14", 7.7475), ("JILA-18", -6.8016), ("BIPM-01", 4.7763)] {
            let row = CODATA_2022_G.iter().find(|r| r.0 == id).unwrap();
            assert!(close((row.1 - mean) / row.2, z, 1e-4), "{id}");
        }
        let (lo, hi) = y.iter().fold((f64::MAX, f64::MIN), |(a, b), &v| (a.min(v), b.max(v)));
        assert!(close((hi - lo) / mean * 1e6, 551.37, 1e-4));
        assert!(gls_mean(&[1.0], &[vec![0.0]]).is_none() && gls_mean(&[1.0, 2.0], &[vec![1.0]]).is_none());
    }

    #[test]
    fn flybys() {
        // Anderson et al. 2008's predictions, mm/s: GLL-I, GLL-II, NEAR, Cassini, Rosetta, MESSENGER.
        let deg = f64::to_radians;
        for (v, din, dout, mm) in [(8949.0, -12.52, -34.15, 4.12296), (8877.0, -34.26, -4.87, -4.67440), (6851.0, -20.76, -71.96, 13.2794),
            (16010.0, -12.92, -4.99, -1.06817), (3863.0, -2.81, -34.29, 2.06651), (4056.0, 31.44, -31.92, 0.0553082)]
        {
            assert!(close(flyby_anderson(v, deg(din), deg(dout)) * 1e3, mm, 1e-5), "{v}");
        }
    }

    #[test]
    fn shielding_and_pioneer() {
        let (earth, moon) = (column_density(6_371e3, 5514.0), column_density(1_737.4e3, 3344.0));
        assert!(close(earth, 7.02594e10, 1e-5) && close(moon, 1.16197e10, 1e-5));
        assert!(close(shielding_loss(MAJORANA_ABSORPTION, earth), 4.59830e-2, 1e-5));
        assert!(close(shielding_loss(MAJORANA_ABSORPTION, moon), 7.755e-3, 1e-3));
        assert!(close(shielding_loss(LLR_ABSORPTION_BOUND, earth), 7.02594e-12, 1e-5));
        // The LLR bound at a solar eclipse: Σ☾ × the Sun's 5.93×10⁻³ m/s² is 6.9×10⁻¹⁵ m/s².
        assert!(close(shielding_loss(LLR_ABSORPTION_BOUND, moon) * 5.9301e-3, 6.9e-15, 0.01));
        assert!(close(pioneer_thermal(10.1e-10, 28.8, 20.0), 6.24127e-10, 1e-5));
    }
}
