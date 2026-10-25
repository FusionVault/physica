"""Regression values on the Eros mesh shipped with polyhedral-gravity-model (unitless, 14744 faces):
analytic edge/face sums against the tetrahedral Gauss-Jacobi cubature (n = 10). Recomputed, not published."""
import math, sys
import numpy as np
from polyhedron_checks import polyhedron, tet_cubature, mesh_check
read = lambda path: np.array([l.split()[1:4] for l in open(path) if l.strip() and not l.startswith("#")][1:], float)
ev, ef = read(sys.argv[1]), read(sys.argv[2]).astype(np.int64)
print(mesh_check(ev, ef))
for p in [(2.0, 0.0, 0.0), (0.0, 0.0, 5.0), (20.0, 30.0, 40.0)]:
    p = np.array(p)
    a, c = polyhedron(p, ev, ef), tet_cubature(p, ev, ef, 10, True)
    print(f"p = {p.tolist()}: V analytic {a[0]:.15e}  GJ10 {c[0]:.15e}  rel {abs(a[0] / c[0] - 1):.0e}")
    print(f"   grad V analytic {a[1]}  GJ10 {c[1]}")
print("trace at origin:", np.trace(polyhedron(np.zeros(3), ev, ef)[2]), " -4pi =", -4 * math.pi)
