"""Gap 2: eigen-solver-free variant for a Rust port: radial Leaver CF + Leaver angular CF, A found by
a scalar secant started at l(l+1) - s(s+1) every time; omega started at the BCW06 fit. Compared with
the spectral-angular solver for all 63 modes at high spins (no tracking in either)."""
import re, sys
import numpy as np
from gap2_leaver_kerr import solve, cf
from gap2_angular_cf_check import A_from_cf

src = open(sys.argv[1]).read()
fits = {}
for row in re.finditer(r"\(\((\d), (-?\d), (\d)\), \[([-\d.e]+), ([-\d.e]+), ([-\d.e]+)\], \[([-\d.e]+), ([-\d.e]+), ([-\d.e]+)\]", src):
    fits[(int(row[1]), int(row[2]), int(row[3]))] = [float(row[i]) for i in range(4, 10)]

def solve_cf_only(a, l, m, n, w, N=800):
    A0 = l * (l + 1) - 2.0 + 0j
    f = lambda w: cf(w, a, l, m, n, A_from_cf(a * w, -2, m, A0), N)
    w0, w1 = w, w * (1 + 1e-6) + 1e-7
    f0, f1 = f(w0), f(w1)
    for _ in range(100):
        w2 = w1 - f1 * (w1 - w0) / (f1 - f0); w0, f0, w1 = w1, f1, w2; f1 = f(w1)
        if abs(w1 - w0) < 1e-14: break
    return w1

worst, bad = 0, []
for (l, m, n), (f1, f2, f3, q1, q2, q3) in sorted(fits.items()):
    for j in [0.3, 0.7, 0.9, 0.98, 0.99]:
        F = f1 + f2 * (1 - j) ** f3; Q = q1 + q2 * (1 - j) ** q3
        g = F - 1j * F / (2 * Q)
        ws, _ = solve(j, l, m, n, g, l * (l + 1) - 2 + 0j, N=800)
        wc = solve_cf_only(j, l, m, n, g)
        d = abs(ws - wc); worst = max(worst, d)
        relfit = abs(ws - g) / abs(g)
        if d > 1e-9 or relfit > 0.05:
            bad.append(((l, m, n), j, ws, wc, relfit))
print(f"max |spectral - CF-only| over 315 cases: {worst:.1e}")
for b in bad: print("  check:", b)
