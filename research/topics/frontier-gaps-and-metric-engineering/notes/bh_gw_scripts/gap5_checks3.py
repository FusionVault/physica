"""Gap 5: (1) verify the corrected small-C series coefficients at several y; (2) n = 0.5 surface
sensitivity of the relativistic solver (stop pressure), to judge Hinderer's erratum Table 1 rows."""
from decimal import Decimal as D
import numpy as np
import gap5_love_tov as g
from gap5_checks2 import k2_exact

for y in [-1.0, 0.0, 0.5, 1.0, 2.0]:
    Y = y + 3
    b = [(2 - y) / Y, (y*y + 6*y - 6) / Y**2, (y**3 + 34*y**2 - 8*y + 12) / (7 * Y**3), (y**4 + 62*y**3 + 84*y**2 + 48*y + 36) / (7 * Y**4)]
    worst = 0
    for C in [1e-4, 3e-4, 1e-3]:
        ex = float(k2_exact(C, y))
        ser = (1 - 2*C)**2 / 2 * sum(bk * C**k for k, bk in enumerate(b))
        worst = max(worst, abs(ser / ex - 1) / C**4)
    print(f"y={y}: max |series/exact - 1| / C^4 = {worst:.3f}  (O(1) => b0..b3 correct)")

# n=0.5 e-polytrope, vary stopping pressure
G = 3.0
e_of = lambda p: p ** (1 / G)
dedp = lambda p: e_of(p) / (G * p) if p > 0 else 0.0
for stop in [1e-10, 1e-14, 1e-18]:
    orig = g.star.__code__
    def star(eos_e, eos_cs2inv, pc, r0=1e-6, stop=stop):
        from scipy.integrate import solve_ivp
        PI = np.pi
        def rhs(r, u):
            m, p, y = u; p = max(p, 0.0); e = eos_e(p); el = 1 / (1 - 2 * m / r)
            dp = -(e + p) * (m + 4 * PI * r**3 * p) / (r * (r - 2 * m)); nup = 2 * el * (m + 4 * PI * p * r**3) / r**2
            Q = 4 * PI * el * (5 * e + 9 * p + (e + p) * eos_cs2inv(p)) - 6 * el / r**2 - nup**2
            return [4 * PI * r * r * e, dp, -(y * y + y * el * (1 + 4 * PI * r * r * (p - e)) + r * r * Q) / r]
        ev = lambda r, u: u[1] - pc * stop; ev.terminal = True; ev.direction = -1
        s = solve_ivp(rhs, [r0, 1e3], [4 * PI / 3 * eos_e(pc) * r0**3, pc, 2.0], events=ev, rtol=1e-12, atol=1e-40, method="LSODA")
        return s.y[0, -1], s.t[-1], s.y[2, -1]
    from scipy.optimize import brentq
    for Ct in [0.01, 0.1, 0.2]:
        f = lambda lp: (lambda M, R, y: M / R)(*star(e_of, dedp, np.exp(lp))) - Ct
        lp = brentq(f, np.log(1e-8), np.log(5.0), xtol=1e-12)
        M, R, y = star(e_of, dedp, np.exp(lp))
        print(f"  n=0.5 e-polytrope stop={stop:g}: C={M/R:.4f} k2={g.k2_of(M/R, y):.5f}")
