"""Gap 1: long thin prisms. Ports geopotential::bodies::prism (Nagy et al. 2000 corner sums with the
crate's safe ln/atan) and its 4x4x4 Gauss-Legendre far field, measures the error of each against
an independent reference, and compares hand-over strategies by accuracy and cost.
Also reproduces Fukushima (2019, GJI 220, 610) Table 4 rigorous values.
Run: python3 prism_checks.py   (numpy only; reference = extended precision or high-order split GL)."""
import math
import numpy as np

EPS = np.finfo(float).eps


def corner_sums(p, lo, hi, dt=np.float64):
    """V, grad V, Hessian at unit G rho by the crate's corner sums, in dtype dt."""
    p, lo, hi = (np.asarray(a, dtype=dt) for a in (p, lo, hi))
    v = dt(0)
    g = np.zeros(3, dtype=dt)
    t = np.zeros((3, 3), dtype=dt)
    for corner in range(8):
        up = [(corner >> a) & 1 == 1 for a in range(3)]
        s = dt(1) if bin(corner).count("1") % 2 == 1 else dt(-1)
        x = np.array([(hi[a] if up[a] else lo[a]) - p[a] for a in range(3)], dtype=dt)
        r = np.sqrt(np.sum(x * x))
        ln = np.zeros(3, dtype=dt)
        at = np.zeros(3, dtype=dt)
        for a in range(3):
            b, c = (a + 1) % 3, (a + 2) % 3
            if p[a] >= hi[a]:
                ln[a] = -np.log(r - x[a])
            elif x[a] >= 0:
                ln[a] = np.log(x[a] + r)
            else:
                ln[a] = np.log((x[b] * x[b] + x[c] * x[c]) / (r - x[a]))
            num = x[b] * x[c]
            at[a] = dt(0) if num == 0 else np.arctan(num / (x[a] * r))
        for a in range(3):
            b, c = (a + 1) % 3, (a + 2) % 3
            v += s * (x[b] * x[c] * ln[a] - dt(0.5) * x[a] * x[a] * at[a])
            g[a] -= s * (x[b] * ln[c] + x[c] * ln[b] - x[a] * at[a])
            t[a, a] -= s * at[a]
            t[b, c] += s * ln[a]
            t[c, b] = t[b, c]
    return v, g, t


_GL = {}


def gl_nodes(n):
    if n not in _GL:
        _GL[n] = np.polynomial.legendre.leggauss(n)
    return _GL[n]


def gl(p, lo, hi, n=(4, 4, 4)):
    """Tensor Gauss-Legendre cubature of 1/r and its derivatives; returns (V, g, T, evaluations)."""
    p, lo, hi = (np.asarray(a, dtype=float) for a in (p, lo, hi))
    half, mid = 0.5 * (hi - lo), 0.5 * (hi + lo)
    (x0, w0), (x1, w1), (x2, w2) = (gl_nodes(k) for k in n)
    X, Y, Z = np.meshgrid(mid[0] + half[0] * x0, mid[1] + half[1] * x1, mid[2] + half[2] * x2, indexing="ij")
    W = np.einsum("i,j,k->ijk", w0, w1, w2) * np.prod(half)
    d = np.stack([X - p[0], Y - p[1], Z - p[2]])
    r2 = np.sum(d * d, axis=0)
    r = np.sqrt(r2)
    v = np.sum(W / r)
    g = np.array([np.sum(W * d[a] / (r2 * r)) for a in range(3)])
    t = np.array([[np.sum(W * (3 * d[a] * d[b] - (r2 if a == b else 0)) / (r2 * r2 * r)) for b in range(3)] for a in range(3)])
    return v, g, t, W.size


def dist_to_box(p, lo, hi):
    q = np.clip(p, lo, hi)
    return float(np.linalg.norm(np.asarray(p) - q))


def gl_split(p, lo, hi, ratio, n=4):
    """Uieda et al. (2016)-style adaptive discretisation: halve along every axis whose length exceeds
    (distance from p to the piece's centre)/ratio, then an n-point rule per axis on each piece."""
    stack = [(np.asarray(lo, float), np.asarray(hi, float))]
    v, g, t, cnt = 0.0, np.zeros(3), np.zeros((3, 3)), 0
    while stack:
        a, b = stack.pop()
        dist = np.linalg.norm(np.asarray(p) - 0.5 * (a + b))
        split = [(b[k] - a[k]) * ratio > dist for k in range(3)]
        if any(split):
            pieces = [(a, b)]
            for k in range(3):
                if split[k]:
                    new = []
                    for (c, d) in pieces:
                        m = 0.5 * (c[k] + d[k])
                        d1, c2 = d.copy(), c.copy()
                        d1[k], c2[k] = m, m
                        new += [(c, d1), (c2, d)]
                    pieces = new
            stack += pieces
        else:
            dv, dg, dt_, k = gl(p, a, b, (n, n, n))
            v, g, t, cnt = v + dv, g + dg, t + dt_, cnt + k
    return v, g, t, cnt


def crate_prism(p, lo, hi, far=10.0):
    """The crate today: corner sums, or 4x4x4 GL beyond `far` body diagonals from the centre."""
    lo, hi = np.asarray(lo, float), np.asarray(hi, float)
    diag = np.linalg.norm(hi - lo)
    dist = np.linalg.norm(np.asarray(p) - 0.5 * (lo + hi))
    if dist > far * diag:
        return gl(p, lo, hi)
    v, g, t = corner_sums(p, lo, hi)
    return v, g, t, 8


def reference(p, lo, hi):
    """Extended-precision corner sums near the body; 8-point split GL (ratio 4) further out."""
    lo, hi = np.asarray(lo, float), np.asarray(hi, float)
    rb = 0.5 * np.linalg.norm(hi - lo)
    dist = np.linalg.norm(np.asarray(p) - 0.5 * (lo + hi))
    if dist < 1.2 * rb:
        v, g, t = corner_sums(p, lo, hi, np.longdouble)
        return float(v), g.astype(float), t.astype(float)
    v, g, t, _ = gl_split(p, lo, hi, ratio=4.0, n=8)
    return v, g, t


def errs(res, ref):
    v, g, t = res[:3]
    vr, gr, tr = ref
    return abs(v / vr - 1), np.linalg.norm(g - gr) / np.linalg.norm(gr), np.max(np.abs(t - tr)) / np.max(np.abs(tr))


def main():
    rng = np.random.default_rng(7)
    dirs = rng.normal(size=(12, 3))
    dirs /= np.linalg.norm(dirs, axis=1)[:, None]

    print("=== Fukushima (2019) Table 4: half-dims (2.5, 3.1, 5.0), point (30, 40, 50), G = rho = 1 ===")
    lo, hi = np.array([-2.5, -3.1, -5.0]), np.array([2.5, 3.1, 5.0])
    p = np.array([30.0, 40.0, 50.0])
    vref, gref, tref = gl_split(p, lo, hi, ratio=8.0, n=10)[:3]
    vl, gl_, tl = corner_sums(p, lo, hi, np.longdouble)
    vd, gd, td = corner_sums(p, lo, hi)
    print(f"V   split-GL {vref:.16e}  long-double corners {float(vl):.16e}  double corners {vd:.16e}  table +4.385413035403429E+00")
    for a, tab in zip(range(3), [-2.638376163531892e-02, -3.515460560153808e-02, -4.380838400289370e-02]):
        print(f"g_{'xyz'[a]} split-GL {gref[a]:.16e}  long-double {float(gl_[a]):.16e}  table {tab:.15e}")
    print(f"Gxx split-GL {tref[0, 0]:.16e}  long-double {float(tl[0, 0]):.16e}  table -4.019417444564296E-04")
    print(f"M/R order 0: {310 / np.linalg.norm(p):.16e} (table order 0 +4.384062043356595E+00)")
    qc = [(2 * 2.5**2 - 3.1**2 - 5**2) / 6, (2 * 3.1**2 - 2.5**2 - 5**2) / 6, (2 * 5**2 - 2.5**2 - 3.1**2) / 6]
    print("quadrupole coefficients (2a²-b²-c²)/6 etc.:", qc, " (table C200 -3.685, C020 -2.005, C002 +5.69)")
    R = np.linalg.norm(p)
    v2 = 310 * (1 / R + (qc[0] * p[0]**2 + qc[1] * p[1]**2 + qc[2] * p[2]**2) / R**5)
    print(f"through order 2: {v2:.16e}")

    print("\n=== Corner-sum loss against distance (max over 12 directions), relative errors V, |g|, T ===")
    for name, (lo, hi) in {"cube 1x1x1": ([-0.5] * 3, [0.5] * 3), "1x1x10": ([-0.5, -0.5, -5], [0.5, 0.5, 5]),
                           "1x1x100 needle": ([-0.5, -0.5, -50], [0.5, 0.5, 50])}.items():
        lo, hi = np.array(lo, float), np.array(hi, float)
        diag = np.linalg.norm(hi - lo)
        abc = np.prod(hi - lo)
        print(name, f"diagonal {diag:.3f}, abc {abc}")
        for k in [1, 3, 10, 30, 100]:
            e = np.zeros(3)
            for d in dirs:
                pt = 0.5 * (lo + hi) + k * diag * d
                e = np.maximum(e, errs(corner_sums(pt, lo, hi), reference(pt, lo, hi)))
            r = k * diag
            print(f"  {k:4d} diagonals: corner sums V {e[0]:.1e} g {e[1]:.1e} T {e[2]:.1e};  eps r^3/abc = {EPS * r**3 / abc:.1e}")

    print("\n=== Corner sums next to a long needle (point 2 from the axis, at mid-length) ===")
    for L in [10.0, 100.0, 1e4]:
        lo, hi = np.array([-0.5, -0.5, -L / 2]), np.array([0.5, 0.5, L / 2])
        pt = np.array([2.0, 0.3, 0.17])
        ref = corner_sums(pt, lo, hi, np.longdouble)
        ref = (float(ref[0]), ref[1].astype(float), ref[2].astype(float))
        e = errs(corner_sums(pt, lo, hi), ref)
        print(f"  L = {L:g}: V {e[0]:.1e} g {e[1]:.1e} T {e[2]:.1e}")

    print("\n=== Strategies for the 1x1x100 needle: max rel error over 12 directions, and kernel evaluations ===")
    lo, hi = np.array([-0.5, -0.5, -50.0]), np.array([0.5, 0.5, 50.0])
    diag = np.linalg.norm(hi - lo)


    def near_cubes(p, lo, hi):
        """Split along each axis into pieces no longer than twice the shortest side, each by crate_prism."""
        L = hi - lo
        n = [int(math.ceil(L[a] / (2 * L.min()))) for a in range(3)]
        v, g, t, cnt = 0.0, np.zeros(3), np.zeros((3, 3)), 0
        for i in range(n[0]):
            for j in range(n[1]):
                for k in range(n[2]):
                    a = lo + L * np.array([i, j, k]) / n
                    b = lo + L * np.array([i + 1, j + 1, k + 1]) / n
                    dv, dg, dt_, c = crate_prism(p, a, b)
                    v, g, t, cnt = v + dv, g + dg, t + dt_, cnt + c
        return v, g, t, cnt


    def bernstein_gl(p, lo, hi, tol=EPS):
        """Per-axis GL order from the Bernstein ellipse of 1/r: rho_a = (d + sqrt(d^2 + h_a^2))/h_a with
        d the distance from p to the box (conservative), n_a = ceil(ln(1/tol)/(2 ln rho_a))."""
        h = 0.5 * (hi - lo)
        d = dist_to_box(p, lo, hi)
        n = []
        for a in range(3):
            rho = (d + math.sqrt(d * d + h[a] ** 2)) / h[a]
            n.append(max(2, int(math.ceil(math.log(1 / tol) / (2 * math.log(rho))))))
        return gl(p, lo, hi, tuple(n))


    def hybrid(p, lo, hi):
        """Corner sums where their loss estimate 4 eps r^3/abc is below the Bernstein GL cost cap,
        else Bernstein-ordered GL."""
        r = np.linalg.norm(p - 0.5 * (lo + hi))
        if 4 * EPS * r**3 / np.prod(hi - lo) < 1e-14 or dist_to_box(p, lo, hi) < 0.25 * np.max(hi - lo):
            v, g, t = corner_sums(p, lo, hi)
            return v, g, t, 8
        return bernstein_gl(p, lo, hi)


    strategies = {"crate (10 diag, 4^3 GL)": crate_prism, "near-cubes x crate": near_cubes,
                  "split GL ratio 8, n=4": lambda p, lo, hi: gl_split(p, lo, hi, 8.0, 4),
                  "Bernstein per-axis GL": bernstein_gl, "hybrid corner/Bernstein": hybrid}
    for k in [0.6, 1, 2, 3, 5, 9.99, 10.01, 20, 50]:
        row = f"{k:6.2f} diag:"
        refs = []
        for d in dirs:
            pt = 0.5 * (lo + hi) + k * diag * d
            if dist_to_box(pt, lo, hi) < 1e-9:
                continue
            refs.append((pt, reference(pt, lo, hi)))
        for name, f in strategies.items():
            e, cnt = np.zeros(3), 0
            for pt, ref in refs:
                out = f(pt, lo, hi)
                e = np.maximum(e, errs(out, ref))
                cnt = max(cnt, out[3])
            row += f" | {name}: {max(e):.0e} ({cnt})"
        print(row)


if __name__ == "__main__":
    main()
