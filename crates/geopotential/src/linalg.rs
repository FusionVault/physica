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

/// Eigenvalues and eigenvectors of a symmetric `n×n` matrix `a` (row-major). Eigenvalue `k` is
/// `values[k]`, with its unit eigenvector in column `k` of `vectors` (row-major), in no particular
/// order. Up to `n = 32` by cyclic Jacobi rotations, which keep small eigenvalues accurate to the
/// matrix's round-off; above that by [`householder_ql`], which costs one `O(n³)` instead of one
/// per sweep.
pub fn symmetric_eigen(a: &[f64], n: usize) -> (Vec<f64>, Vec<f64>) {
    if n > 32 { householder_ql(a, n) } else { jacobi(a, n) }
}

/// [`symmetric_eigen`] by Householder reduction to tridiagonal form and the implicit QL
/// iteration with Wilkinson shifts (Golub & Van Loan, *Matrix Computations*, Algorithms 8.3.1
/// and 8.3.3; EISPACK's `tred2` and `tql2`), `(4/3)n³` operations for the reduction and `O(n²)`
/// per eigenvalue after it. Eigenvalues are accurate to round-off times the matrix's norm.
pub fn householder_ql(a: &[f64], n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut a = a.to_vec();
    let mut q: Vec<f64> = (0..n * n).map(|i| if i / n == i % n { 1.0 } else { 0.0 }).collect();
    // Householder reflections H = I − 2vvᵀ, each clearing one column below the subdiagonal:
    // A ← HAH as A − 2vqᵀ − 2qvᵀ with p = Av, q = p − (vᵀp)v, and Q ← QH.
    for k in 0..n.saturating_sub(2) {
        let norm = (k + 1..n).map(|i| a[i * n + k].powi(2)).sum::<f64>().sqrt();
        if norm == 0.0 {
            continue;
        }
        let mut v = vec![0.0; n];
        let x0 = a[(k + 1) * n + k];
        let alpha = if x0 >= 0.0 { -norm } else { norm };
        for i in k + 1..n {
            v[i] = a[i * n + k];
        }
        v[k + 1] -= alpha;
        let vn = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        v.iter_mut().for_each(|x| *x /= vn);
        let p: Vec<f64> = (0..n).map(|i| (k + 1..n).map(|j| a[i * n + j] * v[j]).sum()).collect();
        let vp: f64 = (k + 1..n).map(|i| v[i] * p[i]).sum();
        let qv: Vec<f64> = (0..n).map(|i| p[i] - vp * v[i]).collect();
        for i in 0..n {
            for j in 0..n {
                a[i * n + j] -= 2.0 * (v[i] * qv[j] + qv[i] * v[j]);
            }
        }
        for i in 0..n {
            let dot: f64 = (k + 1..n).map(|j| q[i * n + j] * v[j]).sum();
            for j in k + 1..n {
                q[i * n + j] -= 2.0 * dot * v[j];
            }
        }
    }
    let mut d: Vec<f64> = (0..n).map(|i| a[i * n + i]).collect();
    let mut e: Vec<f64> = (0..n).map(|i| if i + 1 < n { a[(i + 1) * n + i] } else { 0.0 }).collect();
    // Implicit QL on the tridiagonal (d, e), with the rotations applied to Q's columns.
    for l in 0..n {
        for _ in 0..100 {
            let mut m = l;
            while m + 1 < n && e[m].abs() > f64::EPSILON * (d[m].abs() + d[m + 1].abs()) {
                m += 1;
            }
            if m == l {
                break;
            }
            let mut g = (d[l + 1] - d[l]) / (2.0 * e[l]);
            let r = g.hypot(1.0);
            g = d[m] - d[l] + e[l] / (g + if g >= 0.0 { r } else { -r });
            let (mut s, mut c, mut p) = (1.0, 1.0, 0.0);
            let mut i = m;
            let mut underflow = false;
            while i > l {
                i -= 1;
                let (f, b) = (s * e[i], c * e[i]);
                let r = f.hypot(g);
                e[i + 1] = r;
                if r == 0.0 {
                    d[i + 1] -= p;
                    e[m] = 0.0;
                    underflow = true;
                    break;
                }
                (s, c) = (f / r, g / r);
                g = d[i + 1] - p;
                let r = (d[i] - g) * s + 2.0 * c * b;
                p = s * r;
                d[i + 1] = g + p;
                g = c * r - b;
                for row in q.chunks_mut(n) {
                    let (x, y) = (row[i], row[i + 1]);
                    (row[i], row[i + 1]) = (c * x - s * y, s * x + c * y);
                }
            }
            if underflow {
                continue;
            }
            d[l] -= p;
            e[l] = g;
            e[m] = 0.0;
        }
    }
    (d, q)
}

/// [`symmetric_eigen`] by cyclic Jacobi rotations.
fn jacobi(a: &[f64], n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut a = a.to_vec();
    let mut v: Vec<f64> = (0..n * n).map(|i| if i / n == i % n { 1.0 } else { 0.0 }).collect();
    let scale = a.iter().map(|x| x * x).sum::<f64>();
    for _ in 0..100 {
        let off: f64 = (0..n).flat_map(|p| (p + 1..n).map(move |q| (p, q))).map(|(p, q)| a[p * n + q].powi(2)).sum();
        if off <= 1e-34 * scale {
            break;
        }
        for p in 0..n {
            for q in p + 1..n {
                if a[p * n + q] == 0.0 {
                    continue;
                }
                // The rotation that zeroes a_pq: t = tan φ, the smaller root of t² + 2θt − 1 = 0.
                let theta = (a[q * n + q] - a[p * n + p]) / (2.0 * a[p * n + q]);
                let t = theta.signum() / (theta.abs() + theta.hypot(1.0));
                let c = 1.0 / t.hypot(1.0);
                let s = t * c;
                for k in 0..n {
                    let (kp, kq) = (a[k * n + p], a[k * n + q]);
                    (a[k * n + p], a[k * n + q]) = (c * kp - s * kq, s * kp + c * kq);
                }
                for k in 0..n {
                    let (pk, qk) = (a[p * n + k], a[q * n + k]);
                    (a[p * n + k], a[q * n + k]) = (c * pk - s * qk, s * pk + c * qk);
                    let (kp, kq) = (v[k * n + p], v[k * n + q]);
                    (v[k * n + p], v[k * n + q]) = (c * kp - s * kq, s * kp + c * kq);
                }
            }
        }
    }
    ((0..n).map(|k| a[k * n + k]).collect(), v)
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

    #[test]
    fn householder_ql_agrees_with_jacobi() {
        // A 60×60 symmetric matrix with eigenvalues from 10⁻⁶ to 1 (so clustered small ones):
        // both solvers give A v = λ v and orthonormal columns at round-off, the same spectrum to
        // 10⁻¹⁴ of the norm, and the small eigenvalues to 10⁻¹³ of themselves.
        let n = 60;
        let mut state = 3u64;
        let mut next = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (state >> 11) as f64 / (1u64 << 53) as f64 - 0.5
        };
        let b: Vec<f64> = (0..n * n).map(|_| next()).collect();
        // A = B D Bᵀ with D the chosen eigenvalues and B orthogonalised by Gram–Schmidt.
        let mut qm = b.clone();
        for j in 0..n {
            for k in 0..j {
                let dot: f64 = (0..n).map(|i| qm[i * n + j] * qm[i * n + k]).sum();
                (0..n).for_each(|i| qm[i * n + j] -= dot * qm[i * n + k]);
            }
            let norm = (0..n).map(|i| qm[i * n + j].powi(2)).sum::<f64>().sqrt();
            (0..n).for_each(|i| qm[i * n + j] /= norm);
        }
        let lambda: Vec<f64> = (0..n).map(|k| 10f64.powf(-6.0 * k as f64 / (n - 1) as f64)).collect();
        let a: Vec<f64> = (0..n * n).map(|ij| (0..n).map(|k| qm[(ij / n) * n + k] * lambda[k] * qm[(ij % n) * n + k]).sum()).collect();
        for (values, vectors) in [householder_ql(&a, n), jacobi(&a, n)] {
            for k in 0..n {
                for i in 0..n {
                    let av: f64 = (0..n).map(|j| a[i * n + j] * vectors[j * n + k]).sum();
                    assert!((av - values[k] * vectors[i * n + k]).abs() < 1e-14, "A v = λ v");
                }
                let dot: f64 = (0..n).map(|i| vectors[i * n + k] * vectors[i * n + (k + 1) % n]).sum();
                assert!(dot.abs() < 1e-13 && ((0..n).map(|i| vectors[i * n + k].powi(2)).sum::<f64>() - 1.0).abs() < 1e-13);
            }
            let mut got = values.clone();
            got.sort_by(|a, b| b.total_cmp(a));
            for (g, l) in got.iter().zip(&lambda) {
                assert!((g - l).abs() < 1e-13 * l.max(1e-2), "{g} vs {l}");
            }
        }
    }

    #[test]
    fn symmetric_eigen_of_the_discrete_laplacian() {
        // tridiag(−1, 2, −1) of size n has eigenvalues 2 − 2cos(kπ/(n+1)) in closed form.
        let n: usize = 12;
        let a: Vec<f64> = (0..n * n).map(|i| match (i / n).abs_diff(i % n) { 0 => 2.0, 1 => -1.0, _ => 0.0 }).collect();
        let (mut values, vectors) = symmetric_eigen(&a, n);
        for k in 0..n {
            // A v = λ v, and the columns are orthonormal.
            for i in 0..n {
                let av: f64 = (0..n).map(|j| a[i * n + j] * vectors[j * n + k]).sum();
                assert!((av - values[k] * vectors[i * n + k]).abs() < 1e-14);
            }
            for l in 0..n {
                let dot: f64 = (0..n).map(|i| vectors[i * n + k] * vectors[i * n + l]).sum();
                assert!((dot - if k == l { 1.0 } else { 0.0 }).abs() < 1e-14);
            }
        }
        values.sort_by(f64::total_cmp);
        for (k, v) in values.iter().enumerate() {
            let want = 2.0 - 2.0 * ((k + 1) as f64 * core::f64::consts::PI / (n + 1) as f64).cos();
            assert!((v - want).abs() < 1e-14, "λ{k} = {v} vs {want}");
        }
    }
}
