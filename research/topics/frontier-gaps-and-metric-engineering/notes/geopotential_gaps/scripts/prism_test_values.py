"""Candidate regression values for elongated prisms (recomputed here, not published): reference by
8- and 10-point split Gauss-Legendre (two orders must agree), with the crate's current error."""
import numpy as np
from prism_checks import gl_split, crate_prism
cases = [("needle 1x1x100", [-0.5, -0.5, -50], [0.5, 0.5, 50], [300.0, 400.0, 1200.0]),
         ("needle 1x1x100", [-0.5, -0.5, -50], [0.5, 0.5, 50], [30.0, 40.0, 120.0]),
         ("plate 100x100x1", [-50, -50, -0.5], [50, 50, 0.5], [120.0, -90.0, 300.0])]
for name, lo, hi, p in cases:
    lo, hi, p = (np.array(a, float) for a in (lo, hi, p))
    a, b = gl_split(p, lo, hi, 4.0, 8), gl_split(p, lo, hi, 6.0, 10)
    old = crate_prism(p, lo, hi)
    print(f"{name} at {p.tolist()}: V = {b[0]:.16e} (8 vs 10 pt rel {abs(a[0] / b[0] - 1):.0e}); crate today rel {abs(old[0] / b[0] - 1):.0e}")
    print(f"   grad V = [{b[1][0]:.16e}, {b[1][1]:.16e}, {b[1][2]:.16e}] (spread {np.max(np.abs(a[1] - b[1])) / np.linalg.norm(b[1]):.0e}; crate today {np.linalg.norm(old[1] - b[1]) / np.linalg.norm(b[1]):.0e})")
    print(f"   T_zz = {b[2][2, 2]:.16e}, T_xy = {b[2][0, 1]:.16e} (spread {np.max(np.abs(a[2] - b[2])) / np.max(np.abs(b[2])):.0e}; crate today {np.max(np.abs(old[2] - b[2])) / np.max(np.abs(b[2])):.0e})")
    # far-field identity: through the quadrupole, V = M/R + M sum C_i X_i^2 / R^5, C_x = (2a^2 - b^2 - c^2)/6 (half-dims)
    h = 0.5 * (hi - lo); M = np.prod(hi - lo); R = np.linalg.norm(p)
    C = np.array([(2 * h[0]**2 - h[1]**2 - h[2]**2), (2 * h[1]**2 - h[0]**2 - h[2]**2), (2 * h[2]**2 - h[0]**2 - h[1]**2)]) / 6
    v2 = M / R + M * np.dot(C, p * p) / R**5
    print(f"   monopole+quadrupole {v2:.16e}, rel residual {abs(v2 / b[0] - 1):.1e} (O((h/R)^4) = {(h.max() / R)**4:.1e})")
