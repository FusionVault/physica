"""Gap 1, proposed algorithm: per-axis Gauss-Legendre orders from the Bernstein ellipse of 1/r when
that is affordable, else bisect an elongated prism along its longest side (only for points outside
the closed box, so no split plane passes through p), else the corner sums.
Tested on a cube, a 1x1x100 and a 1x1x10^4 needle and a 100x100x1 plate, near and far.
Run: python3 prism_recursive.py"""
import math
import numpy as np
from prism_checks import corner_sums, gl, dist_to_box, gl_split, errs, EPS

BUDGET = 512  # GL points allowed before bisecting instead


def orders(p, lo, hi, tol=EPS):
    h = 0.5 * (hi - lo)
    d = dist_to_box(p, lo, hi)
    if d == 0.0:
        return None
    n = []
    for a in range(3):
        rho = (d + math.sqrt(d * d + h[a] ** 2)) / h[a]
        n.append(max(1, int(math.ceil(math.log(1 / tol) / (2 * math.log(rho))))))
    return n


def prism(p, lo, hi, count=None, tol=1e-14):
    if count is None:
        count = [0]
    L = hi - lo
    n = orders(p, lo, hi)
    if n is not None and n[0] * n[1] * n[2] <= BUDGET:
        v, g, t, k = gl(p, lo, hi, tuple(n))
        count[0] += k
        return v, g, t, count[0]
    v, g, t = corner_sums(p, lo, hi)
    count[0] += 8
    a = int(np.argmax(L))
    if n is not None and L[a] > 2 * L.min():
        # loss of the corner sums: eps times the largest corner term over |V|
        big = max(np.linalg.norm(np.where([(c >> k) & 1 for k in range(3)], hi, lo) - p) for c in range(8))
        if EPS * big * big * (math.pi / 2 + abs(math.log(big))) > tol * abs(v):
            mid = 0.5 * (lo[a] + hi[a])
            hi1, lo2 = hi.copy(), lo.copy()
            hi1[a], lo2[a] = mid, mid
            v1, g1, t1, _ = prism(p, lo, hi1, count, tol)
            v2, g2, t2, _ = prism(p, lo2, hi, count, tol)
            return v1 + v2, g1 + g2, t1 + t2, count[0]
    return v, g, t, count[0]


def reference(p, lo, hi):
    """Long-double corner sums over pieces no longer than twice the shortest side (absolute loss
    ~1e-19 r^2 per piece) when p is within one longest side of the box, else 8-point split GL."""
    L = hi - lo
    if dist_to_box(p, lo, hi) > L.max():
        return gl_split(p, lo, hi, 4.0, 8)[:3]
    n = [int(math.ceil(L[a] / (2 * L.min()))) for a in range(3)]
    v, g, t = np.longdouble(0), np.zeros(3, np.longdouble), np.zeros((3, 3), np.longdouble)
    edges = [np.linspace(lo[a], hi[a], n[a] + 1) for a in range(3)]
    for i in range(n[0]):
        for j in range(n[1]):
            for k in range(n[2]):
                a = np.array([edges[0][i], edges[1][j], edges[2][k]])
                b = np.array([edges[0][i + 1], edges[1][j + 1], edges[2][k + 1]])
                if np.linalg.norm(p - 0.5 * (a + b)) > 10 * np.linalg.norm(b - a):
                    dv, dg, dt, _ = gl(p, a, b, (8, 8, 8))  # far piece: error ~(1/40)^16
                else:
                    dv, dg, dt = corner_sums(p, a, b, np.longdouble)
                v, g, t = v + dv, g + dg, t + dt
    return float(v), g.astype(float), t.astype(float)


def main():
    rng = np.random.default_rng(11)
    dirs = rng.normal(size=(8, 3))
    dirs /= np.linalg.norm(dirs, axis=1)[:, None]
    bodies = {"cube": ([-0.5] * 3, [0.5] * 3), "needle 1x1x100": ([-0.5, -0.5, -50], [0.5, 0.5, 50]),
              "plate 100x100x1": ([-50, -50, -0.5], [50, 50, 0.5]), "needle 1x1x1e4": ([-0.5, -0.5, -5e3], [0.5, 0.5, 5e3])}
    for name, (lo, hi) in bodies.items():
        lo, hi = np.array(lo, float), np.array(hi, float)
        rb = 0.5 * np.linalg.norm(hi - lo)
        print(name)
        gaps = [0.01, 0.3, 2.0] if "1e4" not in name else [0.3, 2.0]
        # near: points at a gap g outside the surface, on rays from the centre
        for gap in gaps:
            e, cnt, pts = np.zeros(3), 0, 0
            for d in dirs[: (3 if "1e4" in name else 8)]:
                # walk out along d until the distance to the box equals gap
                s0, s1 = 0.0, 2 * rb + gap
                for _ in range(80):
                    s = 0.5 * (s0 + s1)
                    if dist_to_box(0.5 * (lo + hi) + s * d, lo, hi) < gap:
                        s0 = s
                    else:
                        s1 = s
                pt = 0.5 * (lo + hi) + s1 * d
                out = prism(pt, lo, hi)
                e = np.maximum(e, errs(out, reference(pt, lo, hi)))
                cnt = max(cnt, out[3])
            old = max(errs((*corner_sums(pt, lo, hi), 8), reference(pt, lo, hi)))
            print(f"  gap {gap:5.2f} from surface: new V {e[0]:.0e} g {e[1]:.0e} T {e[2]:.0e} (<= {cnt} evals); whole-box corner sums at last point {old:.0e}")
        for k in [1.5, 3, 10, 100]:
            e, cnt = np.zeros(3), 0
            for d in dirs[: (3 if "1e4" in name else 8)]:
                pt = 0.5 * (lo + hi) + k * rb * d
                if dist_to_box(pt, lo, hi) == 0:
                    continue
                out = prism(pt, lo, hi)
                e = np.maximum(e, errs(out, reference(pt, lo, hi)))
                cnt = max(cnt, out[3])
            print(f"  {k:5.1f} R_B from centre: new V {e[0]:.0e} g {e[1]:.0e} T {e[2]:.0e} (<= {cnt} evals)")


if __name__ == "__main__":
    main()
