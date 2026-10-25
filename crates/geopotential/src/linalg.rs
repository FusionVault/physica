//! Dense linear algebra that the potential-field crates share: pivoted Gaussian elimination for a
//! general system, and Cholesky for the symmetric positive-definite normal equations of least
//! squares.
//!
//! ```
//! use geopotential::linalg::{cholesky, cholesky_solve, gauss};
//!
//! // x + 2y = 5, 3x + 4y = 6 → (−4, 4.5).
//! let x = gauss(vec![vec![1.0, 2.0, 5.0], vec![3.0, 4.0, 6.0]]);
//! assert!((x[0] + 4.0).abs() < 1e-14 && (x[1] - 4.5).abs() < 1e-14);
//!
//! // The textbook positive-definite matrix [[4, 12, −16], [12, 37, −43], [−16, −43, 98]]
//! // has L = [[2], [6, 1], [−8, 5, 3]].
//! let mut a = vec![4.0, 0.0, 0.0, 12.0, 37.0, 0.0, -16.0, -43.0, 98.0];
//! cholesky(&mut a, 3);
//! assert_eq!([a[0], a[3], a[4], a[6], a[7], a[8]], [2.0, 6.0, 1.0, -8.0, 5.0, 3.0]);
//! let x = cholesky_solve(&a, 3, &[1.0, 2.0, 3.0]);
//! assert!((4.0 * x[0] + 12.0 * x[1] - 16.0 * x[2] - 1.0).abs() < 1e-10);
//! ```

/// Solve a dense linear system by Gaussian elimination with partial pivoting. Each row carries its
/// right-hand side as its last element: `rows[i] = [a_i0, …, a_i(n−1), b_i]`. For matrices with no
/// symmetry, such as an electrode collocation matrix.
pub fn gauss(mut p: Vec<Vec<f64>>) -> Vec<f64> {
    let n = p.len();
    for col in 0..n {
        let pivot = (col..n).max_by(|&a, &b| p[a][col].abs().total_cmp(&p[b][col].abs())).unwrap();
        p.swap(col, pivot);
        let pivot_row = p[col].clone();
        for row in &mut p[col + 1..] {
            let f = row[col] / pivot_row[col];
            for (x, y) in row[col..].iter_mut().zip(&pivot_row[col..]) {
                *x -= f * y;
            }
        }
    }
    let mut q = vec![0.0; n];
    for row in (0..n).rev() {
        let s: f64 = (row + 1..n).map(|k| p[row][k] * q[k]).sum();
        q[row] = (p[row][n] - s) / p[row][row];
    }
    q
}

/// In-place Cholesky factor `L` (with `A = LLᵀ`) of a symmetric positive-definite `n×n` matrix
/// stored row-major in its lower triangle, as least-squares normal equations are.
///
/// # Panics
/// If the matrix is not positive definite.
pub fn cholesky(a: &mut [f64], n: usize) {
    for j in 0..n {
        let mut d = a[j * n + j];
        for p in 0..j {
            d -= a[j * n + p] * a[j * n + p];
        }
        assert!(d > 0.0, "normal matrix not positive definite at {j}");
        let d = d.sqrt();
        a[j * n + j] = d;
        for i in j + 1..n {
            let mut s = a[i * n + j];
            for p in 0..j {
                s -= a[i * n + p] * a[j * n + p];
            }
            a[i * n + j] = s / d;
        }
    }
}

/// Solve `LLᵀx = b` with the factor [`cholesky`] left in `l`.
pub fn cholesky_solve(l: &[f64], n: usize, b: &[f64]) -> Vec<f64> {
    let mut z = b.to_vec();
    for i in 0..n {
        for p in 0..i {
            z[i] -= l[i * n + p] * z[p];
        }
        z[i] /= l[i * n + i];
    }
    for i in (0..n).rev() {
        for p in i + 1..n {
            z[i] -= l[p * n + i] * z[p];
        }
        z[i] /= l[i * n + i];
    }
    z
}

/// Diagonal of N⁻¹ = L⁻ᵀL⁻¹: the squared norms of the columns of L⁻¹.
pub fn cholesky_inverse_diagonal(l: &[f64], n: usize) -> Vec<f64> {
    let mut diag = vec![0.0; n];
    let mut col = vec![0.0; n];
    for j in 0..n {
        // Solve L·c = e_j.
        col.iter_mut().for_each(|v| *v = 0.0);
        col[j] = 1.0 / l[j * n + j];
        for i in j + 1..n {
            let mut s = 0.0;
            for p in j..i {
                s -= l[i * n + p] * col[p];
            }
            col[i] = s / l[i * n + i];
        }
        // Column j of L⁻¹; (N⁻¹)_jj is its squared norm.
        diag[j] = col[j..].iter().map(|c| c * c).sum();
    }
    diag
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gauss_needs_its_pivot() {
        // A zero on the diagonal: without row exchange this fails.
        let x = gauss(vec![vec![0.0, 1.0, 2.0], vec![1.0, 0.0, 3.0]]);
        assert_eq!(x, vec![3.0, 2.0]);
        // A 6×6 non-symmetric system with a known solution: residual at round-off.
        let n = 6;
        let a: Vec<Vec<f64>> = (0..n).map(|i| (0..n).map(|j| 1.0 / (1.0 + i as f64 + 2.0 * j as f64) + if i == j { 1.0 } else { 0.0 }).collect()).collect();
        let truth: Vec<f64> = (0..n).map(|i| i as f64 - 2.5).collect();
        let rows = a.iter().map(|r| {
            let mut row = r.clone();
            row.push(r.iter().zip(&truth).map(|(x, y)| x * y).sum());
            row
        });
        let x = gauss(rows.collect());
        assert!(x.iter().zip(&truth).all(|(a, b)| (a - b).abs() < 1e-12), "{x:?}");
    }

    #[test]
    fn cholesky_inverse_diagonal_matches_the_inverse() {
        // det = (2·1·3)² = 36, so by cofactors the inverse's diagonal is 1777/36, 136/36 and 4/36;
        // and A·A⁻¹ = I column by column.
        let full = [[4.0, 12.0, -16.0], [12.0, 37.0, -43.0], [-16.0, -43.0, 98.0]];
        let mut a = vec![4.0, 0.0, 0.0, 12.0, 37.0, 0.0, -16.0, -43.0, 98.0];
        cholesky(&mut a, 3);
        let diag = cholesky_inverse_diagonal(&a, 3);
        for (d, exact) in diag.iter().zip([1777.0 / 36.0, 136.0 / 36.0, 4.0 / 36.0]) {
            assert!((d / exact - 1.0).abs() < 1e-12, "{d} vs {exact}");
        }
        for (j, d) in diag.iter().enumerate() {
            let mut e = [0.0; 3];
            e[j] = 1.0;
            let col = cholesky_solve(&a, 3, &e);
            assert!((col[j] - d).abs() < 1e-10 * d.abs(), "{j}: {} vs {d}", col[j]);
            for (i, row) in full.iter().enumerate() {
                let v: f64 = row.iter().zip(&col).map(|(x, y)| x * y).sum();
                assert!((v - e[i]).abs() < 1e-10);
            }
        }
    }
}
