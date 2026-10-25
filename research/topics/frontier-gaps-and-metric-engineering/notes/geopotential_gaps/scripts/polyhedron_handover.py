"""Where to hand the polyhedron over to the tetrahedral Gauss-Jacobi cubature, and how many points:
max rel error in V (math.fsum of per-face terms, so the test's own summation is exact enough) for the
cube mesh (reference: prism split GL) and Eros (reference: GJ n=10), against distance in body radii."""
import math, sys
import numpy as np
from polyhedron_checks import polyhedron, duffy_rule, verts, faces
from prism_checks import gl_split

def tet_v(p, V, F, n):
    A = V[np.unique(F)].mean(axis=0)
    B, C, D = (V[F[:, k]] for k in range(3))
    vol6 = np.einsum("ij,ij->i", B - A, np.cross(C - A, D - A))
    U, Vv, W, Wt = duffy_rule(n, True)
    X = A + U[None, :, None] * ((B - A)[:, None, :] + Vv[None, :, None] * ((C - B)[:, None, :] + W[None, :, None] * (D - C)[:, None, :]))
    r = np.linalg.norm(X - p, axis=2)
    return math.fsum((vol6[:, None] * Wt[None, :] / r).ravel())

rng = np.random.default_rng(5)
dirs = rng.normal(size=(6, 3)); dirs /= np.linalg.norm(dirs, axis=1)[:, None]
lo, hi = np.array([-0.5] * 3), np.array([0.5] * 3)
rad = math.sqrt(3) / 2
print("cube mesh (12 faces), distance in circumradii; analytic = crate edge/face sums")
for k in [2, 5, 10, 20, 50]:
    row = {}
    for d in dirs:
        p = k * rad * d
        ref = gl_split(p, lo, hi, 4.0, 8)[0]
        row["analytic"] = max(row.get("analytic", 0), abs(polyhedron(p, verts, faces)[0] / ref - 1))
        for n in [3, 4, 5, 6, 8]:
            row[f"GJ{n}"] = max(row.get(f"GJ{n}", 0), abs(tet_v(p, verts, faces, n) / ref - 1))
    print(f"  {k:3d}: " + "  ".join(f"{a} {b:.0e}" for a, b in row.items()))
if len(sys.argv) > 2:
    read = lambda path: np.array([l.split()[1:4] for l in open(path) if l.strip() and not l.startswith("#")][1:], float)
    ev, ef = read(sys.argv[1]), read(sys.argv[2]).astype(np.int64)
    c = ev[np.unique(ef)].mean(axis=0)
    size = np.max(np.linalg.norm(ev[np.unique(ef)] - c, axis=1))
    print("Eros (14744 faces), distance in max radii from the vertex centroid")
    for k in [2, 5, 10, 20, 50, 200]:
        row = {}
        for d in dirs[:3]:
            p = c + k * size * d
            ref = tet_v(p, ev, ef, 10)
            row["analytic"] = max(row.get("analytic", 0), abs(polyhedron(p, ev, ef)[0] / ref - 1))
            for n in [3, 4, 5, 6, 8]:
                row[f"GJ{n}"] = max(row.get(f"GJ{n}", 0), abs(tet_v(p, ev, ef, n) / ref - 1))
        print(f"  {k:3d}: " + "  ".join(f"{a} {b:.0e}" for a, b in row.items()))
