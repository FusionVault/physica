"""Cost and accuracy of the proposed prism algorithm at tolerance 1e-13 vs 1e-14, near a plate and a needle."""
import numpy as np
from prism_recursive import prism, reference
from prism_checks import errs
for name, lo, hi, pt in [("plate 100x100x1, 0.3 above the face", [-50, -50, -0.5], [50, 50, 0.5], [3.1, -7.2, 0.8]),
                         ("needle 1x1x1e4, 2 from the side", [-0.5, -0.5, -5e3], [0.5, 0.5, 5e3], [2.5, 0.3, 17.0]),
                         ("needle 1x1x100, 0.3 from the side", [-0.5, -0.5, -50], [0.5, 0.5, 50], [0.8, 0.1, 3.0])]:
    lo, hi, pt = np.array(lo, float), np.array(hi, float), np.array(pt, float)
    ref = reference(pt, lo, hi)
    for tol in [1e-14, 1e-13, 1e-12]:
        out = prism(pt, lo, hi, tol=tol)
        print(f"{name}: tol {tol:.0e}: max rel err {max(errs(out, ref)):.1e}, {out[3]} kernel evaluations")
