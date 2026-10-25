"""Gap 2: eigen-solver-free spectral angular step for a Rust port: shifted inverse iteration on the
Cook-Zalutskiy matrix with the bilinear Rayleigh quotient v^T M v / v^T v (matrix is complex symmetric).
Compared with numpy's eigenvalue closest to the shift, at converged QNMs of several modes and spins."""
import numpy as np
from gap2_leaver_kerr import F, G, H, sep_const, solve

def matrix(c, l, m, s=-2, extra=25):
    lmin = max(abs(m), abs(s)); ls = np.arange(lmin, l + extra + 1); n = len(ls)
    M = np.zeros((n, n), complex)
    for i, L in enumerate(ls):
        for j, Lp in enumerate(ls):
            if Lp == L - 2: M[i, j] = -c*c*F(s, Lp, m)*F(s, Lp+1, m)
            elif Lp == L - 1: M[i, j] = -c*c*F(s, Lp, m)*(H(s, Lp+1, m)+H(s, Lp, m)) + 2*c*s*F(s, Lp, m)
            elif Lp == L: M[i, j] = Lp*(Lp+1) - s*(s+1) - c*c*(F(s, Lp, m)*G(s, Lp+1, m) + G(s, Lp, m)*F(s, Lp-1, m) + H(s, Lp, m)**2) + 2*c*s*H(s, Lp, m)
            elif Lp == L + 1: M[i, j] = -c*c*G(s, Lp, m)*(H(s, Lp-1, m)+H(s, Lp, m)) + 2*c*s*G(s, Lp, m)
            elif Lp == L + 2: M[i, j] = -c*c*G(s, Lp, m)*G(s, Lp-1, m)
    return M, ls

def inv_iter(M, shift, ls, l, iters=8):
    v = (ls == l).astype(complex)       # start from the spherical harmonic l
    lam = shift
    for _ in range(iters):
        w = np.linalg.solve(M - lam * np.eye(len(v)), v)   # Rust: banded Gaussian elimination
        v = w / np.sqrt(w @ w)
        lam = (v @ M @ v) / (v @ v)                       # bilinear Rayleigh quotient (Rayleigh-quotient iteration)
    return lam

worst = 0
for (l, m, n, a, w0) in [(2, 2, 0, 0.99, 0.87 - 0.03j), (3, -2, 0, 0.9, 0.519 - 0.09j), (4, 2, 1, 0.98, 1.2 - 0.25j), (4, -1, 2, 0.7, 0.75 - 0.47j), (2, -2, 0, 0.98, 0.29 - 0.088j)]:
    w, A = solve(a, l, m, n, w0, l*(l+1) - 2 + 0j, N=800)
    M, ls = matrix(a * w, l, m)
    A_ii = inv_iter(M, l*(l+1) - 2.0 + 0j, ls, l)
    worst = max(worst, abs(A_ii - A))
    print(f"({l},{m},{n}) a={a}: omega={w:.8f}  A(numpy closest)={A:.10f}  A(inverse iteration from l(l+1)-2)={A_ii:.10f}")
print("max diff", worst)
