"""Gap 2: true errors of the BCW06 fits (library ringdown.rs FITS) against exact Leaver QNMs on
j = 0, 0.05, ..., 0.95, 0.98, 0.99; compared with the maxima BCW06 state (Tables VIII-X)."""
import re, sys
import numpy as np
from gap2_leaver_kerr import solve
src = open(sys.argv[1]).read()
rows = re.findall(r"\(\((\d), (-?\d), (\d)\), \[([-\d.e]+), ([-\d.e]+), ([-\d.e]+)\], \[([-\d.e]+), ([-\d.e]+), ([-\d.e]+)\], ([\d.]+), ([\d.]+)\)", src)
js = list(np.arange(0, 0.951, 0.05)) + [0.98, 0.99]
over = []
for r in rows:
    l, m, n = int(r[0]), int(r[1]), int(r[2]); f = list(map(float, r[3:6])); q = list(map(float, r[6:9])); fe, qe = float(r[9]), float(r[10])
    eF = eQ = 0; jF = jQ = 0
    for j in js:
        F = f[0] + f[1] * (1 - j) ** f[2]; Q = q[0] + q[1] * (1 - j) ** q[2]
        w, _ = solve(j, l, m, n, F - 1j * F / (2 * Q), l * (l + 1) - 2 + 0j, N=800)
        Fx, Qx = w.real, w.real / (-2 * w.imag)
        if abs(F / Fx - 1) * 100 > eF: eF, jF = abs(F / Fx - 1) * 100, j
        if abs(Q / Qx - 1) * 100 > eQ: eQ, jQ = abs(Q / Qx - 1) * 100, j
    flag = " <-- Q over stated" if eQ > qe else ""
    flag += " <-- F over stated" if eF > fe else ""
    print(f"({l},{m:+d},{n}) F max {eF:5.2f}% at j={jF:.2f} (stated {fe})  Q max {eQ:5.2f}% at j={jQ:.2f} (stated {qe}){flag}")
