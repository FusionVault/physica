"""Gap 2: Leaver's angular continued fraction (coefficients as in Berti's QNMLeaverKerr.nb, after
Leaver 1985 Eqs. 18-27) versus the Cook-Zalutskiy spectral eigenvalue, at the converged QNMs.
Shows that a Rust port can avoid a complex eigen-solver by root-finding A with a scalar secant."""
import numpy as np
from gap2_leaver_kerr import sep_const, track, A_closest

def ang_cf(A, c, s, m, N=300):
    k1, k2 = abs(m - s) / 2, abs(m + s) / 2
    al = lambda n: -2 * (n + 1) * (n + 2 * k1 + 1)
    be = lambda n: n * (n - 1) + 2 * n * (k1 + k2 + 1 - 2 * c) - (2 * c * (2 * k1 + s + 1) - (k1 + k2) * (k1 + k2 + 1)) - (c * c + s * (s + 1) + A)
    ga = lambda n: 2 * c * (n + k1 + k2 + s)
    r = 0.0
    for n in range(N, 0, -1):
        r = ga(n) / (be(n) - al(n) * r)
    return be(0) / al(0) - r   # Berti: Leaver33Ang = beta0/alpha0 - Leaver31Ang

def A_from_cf(c, s, m, A0):
    x0, x1 = A0, A0 * (1 + 1e-6) + 1e-6
    f0, f1 = ang_cf(x0, c, s, m), ang_cf(x1, c, s, m)
    for _ in range(60):
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0); x0, f0, x1 = x1, f1, x2; f1 = ang_cf(x1, c, s, m)
        if abs(x1 - x0) < 1e-14: break
    return x1

w0 = 0.37367168441804 - 0.08896231568894j
res = track(2, 2, 0, [0.4, 0.7, 0.9, 0.99], w0)
for a, (w, A) in res.items():
    Acf = A_from_cf(a * w, -2, 2, A)
    print(f"a={a}: A_spectral = {A:.12f}, A_LeaverCF = {Acf:.12f}, |diff| = {abs(A-Acf):.1e}")
print("Berti notebook a/M=0.4 (2M=1: ar=0.2): his iterate 0.879684276903489/2 =", 0.879684276903489/2, " ours:", res[0.4][0])
