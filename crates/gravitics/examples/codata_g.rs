//! Why Newton's constant is the least well known: CODATA 2022 rebuilt from its sixteen
//! measurements. Their weighted mean is 6.674 30, but they disagree with their own error bars by
//! a Birge ratio of 3.6, so CODATA expands the uncertainty 3.9-fold
//! ([CODATA 2022](https://arxiv.org/abs/2409.03787)).
//!
//! ```text
//! cargo run -p gravitics --example codata_g
//! ```

use gravitics::anomalies::{CODATA_2022_G, CODATA_2022_G_CORRELATIONS, birge_ratio, gls_mean};

fn main() {
    let n = CODATA_2022_G.len();
    let mut cov = vec![vec![0.0; n]; n];
    for (i, row) in CODATA_2022_G.iter().enumerate() {
        cov[i][i] = row.2 * row.2;
    }
    let at = |id: &str| CODATA_2022_G.iter().position(|r| r.0 == id).expect("listed");
    for &(a, b, r) in &CODATA_2022_G_CORRELATIONS {
        let (i, j) = (at(a), at(b));
        cov[i][j] = r * CODATA_2022_G[i].2 * CODATA_2022_G[j].2;
        cov[j][i] = cov[i][j];
    }
    let values: Vec<f64> = CODATA_2022_G.iter().map(|r| r.1).collect();
    let (mean, u, chi2) = gls_mean(&values, &cov).expect("positive definite");
    let birge = birge_ratio(chi2, n - 1);
    println!("G = {mean:.6} ± {:.6} ×10⁻¹¹ (stated errors), χ² = {chi2:.1} for {} degrees of freedom", u, n - 1);
    println!("Birge ratio {birge:.2}: expanded by 3.9, u = {:.1} ppm", 3.9 * u / mean * 1e6);
    for row in &CODATA_2022_G {
        println!("  {:<9} {:+6.2} σ from the mean", row.0, (row.1 - mean) / row.2);
    }
    assert!((mean - 6.674_30).abs() < 1e-6 && (birge - 3.61).abs() < 0.01);
}

#[test]
fn rebuilds_codata() {
    main();
}
