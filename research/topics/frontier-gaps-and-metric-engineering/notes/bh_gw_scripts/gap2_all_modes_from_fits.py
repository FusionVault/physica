"""Gap 2: polish every BCW06 mode (l=2..4, all m, n=0..2) from the BCW06 fit as initial guess, at the
eleven spins of BCW06 Tables II-IV, and compare with those tables (4 decimals). No spin tracking.
FITS and TABLES are parsed from the library's own ringdown.rs (transcribed from BCW06)."""
import re, sys
import numpy as np
from gap2_leaver_kerr import solve

src = open(sys.argv[1]).read()
fits = {}
for l, m, n, f1, f2, f3 in re.findall(r"\(\((\d), (-?\d), (\d)\), \[([-\d.]+), ([-\d.]+), ([-\d.]+)\], \[([-\d.]+), ([-\d.]+), ([-\d.]+)\]".replace(r", \[([-\d.]+), ([-\d.]+), ([-\d.]+)\]$", ""), src)[:0]:
    pass
for row in re.finditer(r"\(\((\d), (-?\d), (\d)\), \[([-\d.e]+), ([-\d.e]+), ([-\d.e]+)\], \[([-\d.e]+), ([-\d.e]+), ([-\d.e]+)\]", src):
    l, m, n = int(row[1]), int(row[2]), int(row[3])
    fits[(l, m, n)] = [float(row[i]) for i in range(4, 10)]
tab = re.search(r'const TABLES: &str = "\\\n(.*?)";', src, re.S).group(1)
table = {}
for line in tab.strip().splitlines():
    parts = line.split()
    l, n, j = int(parts[0]), int(parts[1]), float(parts[2])
    for k, pair in enumerate(parts[3:]):
        re_, im_ = map(float, pair.split(","))
        table[(l, l - k, n, j)] = re_ - 1j * im_
print(len(fits), "fits,", len(table), "table entries")

worst, fails = 0, []
for (l, m, n), (f1, f2, f3, q1, q2, q3) in sorted(fits.items()):
    for j in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.98]:
        F = f1 + f2 * (1 - j) ** f3; Q = q1 + q2 * (1 - j) ** q3
        guess = F - 1j * F / (2 * Q)
        w, A = solve(j, l, m, n, guess, l * (l + 1) - 2 + 0j, N=800)
        ref = table[(l, m, n, j)]
        d = max(abs(w.real - ref.real), abs(w.imag - ref.imag))
        worst = max(worst, d)
        if d > 6e-5:
            fails.append(((l, m, n), j, w, ref))
print(f"worst |Leaver - table| over all {len(fits)*11} cases: {worst:.1e}")
for f in fails:
    print("  mismatch", f)
