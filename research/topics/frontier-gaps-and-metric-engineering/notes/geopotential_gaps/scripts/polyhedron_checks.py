"""Gap 2: the polyhedron's far field and mesh checks.
Ports geopotential::bodies::polyhedron (Werner & Scheeres 1997 edge/face sums), measures its far-field
loss on a cube mesh and on the Eros mesh shipped with ESA's polyhedral-gravity-model (data only),
tries a log1p fix for L_e, a Duffy-collapsed Gauss-Legendre cubature over tetrahedra (face + centroid)
for the far field, and an O(F log F) directed-edge closure/orientation check plus signed volume.
Run: python3 polyhedron_checks.py <path to Eros.node> <path to Eros.face>"""
import math
import sys
import numpy as np
from prism_checks import gl_split, corner_sums

EPS = np.finfo(float).eps


def polyhedron(p, verts, faces, l_form="crate"):
    """Vectorised port of the crate's polyhedron(); l_form 'crate' = ln((s+e)/gap), 'log1p' = log1p(2e/gap)."""
    P = verts[faces]  # F x 3 x 3
    r = P - p
    l = np.linalg.norm(r, axis=2)
    n = np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0])
    n /= np.linalg.norm(n, axis=1)[:, None]
    triple = np.einsum("ij,ij->i", r[:, 0], np.cross(r[:, 1], r[:, 2]))
    triple = np.where(triple == 0, -0.0, triple)
    den = (l[:, 0] * l[:, 1] * l[:, 2] + l[:, 0] * np.einsum("ij,ij->i", r[:, 1], r[:, 2])
           + l[:, 1] * np.einsum("ij,ij->i", r[:, 2], r[:, 0]) + l[:, 2] * np.einsum("ij,ij->i", r[:, 0], r[:, 1]))
    omega = 2 * np.arctan2(triple, den)
    rf = np.einsum("ij,ij->i", n, r[:, 0])
    v = -0.5 * np.sum(rf * rf * omega)
    g = np.sum(n * (rf * omega)[:, None], axis=0)
    t = -np.einsum("i,ij,ik->jk", omega, n, n)
    for k in range(3):
        i, j = k, (k + 1) % 3
        edge = P[:, j] - P[:, i]
        e = np.linalg.norm(edge, axis=1)
        ne = np.cross(edge, n)
        ne /= np.linalg.norm(ne, axis=1)[:, None]
        s = l[:, i] + l[:, j]
        ab = np.einsum("ij,ij->i", r[:, i], r[:, j])
        c = np.cross(r[:, i], r[:, j])
        gap = np.where(ab < 0, 2 * np.einsum("ij,ij->i", c, c) / ((l[:, i] * l[:, j] - ab) * (s + e)), s - e)
        le = np.log((s + e) / gap) if l_form == "crate" else np.log1p(2 * e / gap)
        ner = np.einsum("ij,ij->i", ne, r[:, i])
        v += 0.5 * np.sum(np.einsum("ij,ij->i", n, r[:, i]) * ner * le)
        g -= np.sum(n * (ner * le)[:, None], axis=0)
        t += np.einsum("i,ij,ik->jk", le, n, ne)
    t = 0.5 * (t + t.T)
    return v, g, t


def duffy_rule(n, jacobi=False):
    """Nodes (u, v, w) in [0,1]^3 and weights for int_T f = 6 Vol sum w f(x(u,v,w)),
    x = A + u[(B-A) + v((C-B) + w(D-C))], Jacobian u^2 v."""
    if jacobi:
        from scipy.special import roots_jacobi
        xu, wu = roots_jacobi(n, 0, 2)  # weight (1+x)^2 on [-1,1]
        xv, wv = roots_jacobi(n, 0, 1)
        u, wu = (xu + 1) / 2, wu / 8
        vv, wv = (xv + 1) / 2, wv / 4
    else:
        x, w = np.polynomial.legendre.leggauss(n)
        u, vv = (x + 1) / 2, (x + 1) / 2
        wu, wv = w / 2 * u ** 2, w / 2 * vv
    x, w = np.polynomial.legendre.leggauss(n)
    ww, xw = w / 2, (x + 1) / 2
    U, V, W = np.meshgrid(u, vv, xw, indexing="ij")
    Wt = np.einsum("i,j,k->ijk", wu, wv, ww)
    return U.ravel(), V.ravel(), W.ravel(), Wt.ravel()


def tet_cubature(p, verts, faces, n=4, jacobi=False):
    """Far field: signed tetrahedra (centroid A, face B C D) by a collapsed Gauss rule."""
    A = verts.mean(axis=0)
    B, C, D = (verts[faces[:, k]] for k in range(3))
    vol6 = np.einsum("ij,ij->i", B - A, np.cross(C - A, D - A))  # 6 x signed volume
    U, V, W, Wt = duffy_rule(n, jacobi)
    # points: F x Q x 3
    X = A + U[None, :, None] * ((B - A)[:, None, :] + V[None, :, None] * ((C - B)[:, None, :] + W[None, :, None] * (D - C)[:, None, :]))
    d = X - p
    r2 = np.sum(d * d, axis=2)
    r = np.sqrt(r2)
    wq = vol6[:, None] * Wt[None, :]
    v = np.sum(wq / r)
    g = np.einsum("fq,fqa->a", wq / (r2 * r), d)
    t = np.einsum("fq,fqa,fqb->ab", wq * 3 / (r2 * r2 * r), d, d) - np.eye(3) * np.sum(wq / (r2 * r))
    return v, g, t, len(faces) * len(Wt)


def mesh_check(verts, faces):
    """Closed and consistently oriented iff every directed edge (i, j) occurs exactly once and its reverse
    (j, i) exactly once; outward iff the signed volume sum (v0 . v1 x v2)/6 is positive."""
    edges = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
    if np.any(edges[:, 0] == edges[:, 1]):
        return "degenerate face (repeated vertex)"
    keys = edges[:, 0].astype(np.int64) * len(verts) + edges[:, 1]
    rev = edges[:, 1].astype(np.int64) * len(verts) + edges[:, 0]
    uniq, cnt = np.unique(keys, return_counts=True)
    if np.any(cnt > 1):
        return f"{int(np.sum(cnt > 1))} directed edges used twice: inconsistent orientation or non-manifold"
    missing = ~np.isin(rev, keys)
    if np.any(missing):
        return f"{int(np.sum(missing))} edges without a partner: the mesh is open"
    vol = np.sum(np.einsum("ij,ij->i", verts[faces[:, 0]], np.cross(verts[faces[:, 1]], verts[faces[:, 2]]))) / 6
    return f"closed and consistent, volume {vol:+.15g} ({'outward' if vol > 0 else 'INWARD: flip all faces'})"


def rel(a, b):
    return abs(a / b - 1)


# Cube mesh as in the crate's tests
verts = np.array([[(i & 1) - 0.5, ((i >> 1) & 1) - 0.5, ((i >> 2) & 1) - 0.5] for i in range(8)], float)
quads = [[0, 4, 6, 2], [1, 3, 7, 5], [0, 1, 5, 4], [2, 6, 7, 3], [0, 2, 3, 1], [4, 5, 7, 6]]
faces = np.array([t for a, b, c, d in quads for t in ([a, b, c], [a, c, d])])
lo, hi = np.array([-0.5] * 3), np.array([0.5] * 3)
def main():
    print("mesh check, cube:", mesh_check(verts, faces))
    print("mesh check, inverted cube:", mesh_check(verts, faces[:, [1, 0, 2]]))
    bad = faces.copy(); bad[3] = bad[3][[1, 0, 2]]
    print("mesh check, one face flipped:", mesh_check(verts, bad))
    print("mesh check, one face missing:", mesh_check(verts, faces[1:]))
    pin = np.array([0.1, 0.2, -0.05])
    print("inverted cube negates V:", polyhedron(pin, verts, faces[:, [1, 0, 2]])[0], -polyhedron(pin, verts, faces)[0],
          " trace inside", np.trace(polyhedron(pin, verts, faces[:, [1, 0, 2]])[2]), "(+4pi when inverted)")
    print("one flipped face, V error at an inside point:", polyhedron(pin, verts, bad)[0] / polyhedron(pin, verts, faces)[0] - 1)

    print("\npolyhedral-gravity-model cube [-1,1]^3 reference rows (G = 6.67430e-11, rho = 1):")
    G = 6.67430e-11
    v2 = verts * 2
    for pt, vref, gref in [([-5, -5, -5], 6.16562723600017E-11, [4.11088282102974e-12] * 3),
                           ([-5, -4.5, -5], 6.37068727767111E-11, [4.53491850099316e-12, 4.08129150394744e-12, 4.53491850099315e-12])]:
        v, g, _ = polyhedron(np.array(pt, float), v2, faces)
        vp, gp, _ = corner_sums(np.array(pt, float), -np.ones(3), np.ones(3))
        print(f"  {pt}: polyhedron G*V {G * v:.15e}  prism {G * vp:.15e}  file {vref:.15e}  rel {rel(G * v, vref):.1e};"
              f" g {G * g}  file {gref}")

    print("\nFar field, cube mesh: max rel error over 6 directions (reference: prism split GL)")
    rng = np.random.default_rng(3)
    dirs = rng.normal(size=(6, 3)); dirs /= np.linalg.norm(dirs, axis=1)[:, None]
    for k in [3, 10, 30, 100, 300, 1000]:
        out = {}
        for d in dirs:
            pt = k * d
            ref = gl_split(pt, lo, hi, 4.0, 8)[:3]
            for name, res in [("crate ln", polyhedron(pt, verts, faces)), ("log1p", polyhedron(pt, verts, faces, "log1p")),
                              ("tet GL n=4", tet_cubature(pt, verts, faces, 4)), ("tet GJ n=4", tet_cubature(pt, verts, faces, 4, True))]:
                e = max(rel(res[0], ref[0]), np.linalg.norm(res[1] - ref[1]) / np.linalg.norm(ref[1]),
                        np.max(np.abs(res[2] - ref[2])) / np.max(np.abs(ref[2])))
                out[name] = max(out.get(name, 0), e)
        print(f"  {k:5d} sizes: " + "  ".join(f"{a} {b:.1e}" for a, b in out.items()))

    if len(sys.argv) > 2:
        def read(path, cols):
            rows = [l.split() for l in open(path) if l.strip() and not l.startswith("#")]
            return np.array([r[1:1 + cols] for r in rows[1:]], dtype=float)
        ev = read(sys.argv[1], 3)
        ef = read(sys.argv[2], 3).astype(np.int64)
        print(f"\nEros mesh: {len(ev)} vertices listed, {len(np.unique(ef))} used, {len(ef)} faces")
        print("  mesh check:", mesh_check(ev, ef))
        c = ev[np.unique(ef)].mean(axis=0)
        size = np.max(np.linalg.norm(ev[np.unique(ef)] - c, axis=1))
        print(f"  max radius from vertex centroid {size:.6f}")
        tr_in = np.trace(polyhedron(c, ev, ef)[2])
        print(f"  trace at the vertex centroid {tr_in:.15f} (-4pi = {-4 * math.pi:.15f})")
        for k in [2, 10, 30, 100, 1000]:
            out = {}
            for d in dirs[:3]:
                pt = c + k * size * d
                ref = tet_cubature(pt, ev, ef, 8, True)
                for name, res in [("crate ln", polyhedron(pt, ev, ef)), ("log1p", polyhedron(pt, ev, ef, "log1p")),
                                  ("tet GJ n=4", tet_cubature(pt, ev, ef, 4, True)), ("tet GJ n=6", tet_cubature(pt, ev, ef, 6, True))]:
                    e = max(rel(res[0], ref[0]), np.linalg.norm(res[1] - ref[1]) / np.linalg.norm(ref[1]),
                            np.max(np.abs(res[2] - ref[2])) / np.max(np.abs(ref[2])))
                    out[name] = max(out.get(name, 0), e)
            print(f"  {k:5d} radii: " + "  ".join(f"{a} {b:.1e}" for a, b in out.items()))
        pt = c + 3 * size * dirs[0]
        print("  Eros at 3 radii, V by analytic / tet GJ n=8 / n=10:", polyhedron(pt, ev, ef)[0], tet_cubature(pt, ev, ef, 8, True)[0],
              tet_cubature(pt, ev, ef, 10, True)[0])


if __name__ == "__main__":
    main()
