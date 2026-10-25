"""Scan spin and inclination: the range of each Kerr shadow size measure relative to sqrt(27)."""
import numpy as np
from kerr_shadow import measures, pole_radius, xi_eta

s27 = 27 ** 0.5
res = {k: (9, None, -9, None) for k in ["meanR_origin", "medianR", "areal", "halfwidth", "halfheight"]}
for a in [0.0, 0.3, 0.6, 0.9, 0.99, 0.999, 1 - 1e-6]:
    for ideg in [0.5, 10, 30, 50, 70, 90]:
        m = measures(a, np.radians(ideg), 40001)
        m["halfwidth"], m["halfheight"] = m["width"] / 2, m["height"] / 2
        for k in res:
            lo, alo, hi, ahi = res[k]
            d = m[k] / s27 - 1
            if d < lo: lo, alo = d, (a, ideg)
            if d > hi: hi, ahi = d, (a, ideg)
            res[k] = (lo, alo, hi, ahi)
for k, (lo, alo, hi, ahi) in res.items():
    print(f"{k:13s}: delta from {lo:+.5f} at (a,i)={alo} to {hi:+.5f} at {ahi}")
r = pole_radius(1.0)
print("pole-on a=1: r =", r, " radius sqrt(8r^3/(r^2-1)) =", np.sqrt(8 * r**3 / (r * r - 1)), " 2(1+sqrt2) =", 2 * (1 + 2 ** 0.5),
      " delta =", 2 * (1 + 2 ** 0.5) / s27 - 1)
