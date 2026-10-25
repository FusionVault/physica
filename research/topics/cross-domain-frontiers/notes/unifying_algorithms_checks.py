"""Benchmarks for the shared algorithms in unifying_algorithms_across_crates.md.

Standard library only. Run as:  python3 -I unifying_algorithms_checks.py
Each check prints PASS/FAIL with the computed value, the reference and its source tag.
"""
import cmath
import math
import random

G = 6.67430e-11
C = 299_792_458.0
MU0 = 1.25663706212e-6
EPS0 = 8.8541878128e-12
QE = 1.602176634e-19
ME = 9.1093837015e-31
MP = 1.67262192369e-27
GM_SUN = 1.32712440018e20
R_SUN = 6.957e8
M_SUN = GM_SUN / G
AU = 1.495978707e11
PC = 3.0856775814913673e16
RAD2AS = 180 / math.pi * 3600

FAILS = []


def check(name, got, want, rel=None, abs_=None):
    if rel is not None:
        ok = abs(got - want) <= rel * abs(want)
    else:
        ok = abs(got - want) <= abs_
    print(f"{'PASS' if ok else 'FAIL'}  {name}: got {got!r}, want {want!r}")
    if not ok:
        FAILS.append(name)


def info(name, value):
    print(f"INFO  {name}: {value}")


# ---------------------------------------------------------------- shared numerics
def gauss_legendre(n):
    """Nodes and weights on [-1, 1] by Newton on P_n."""
    xs, ws = [], []
    for i in range(1, n + 1):
        x = math.cos(math.pi * (i - 0.25) / (n + 0.5))
        for _ in range(100):
            p0, p1 = 1.0, x
            for k in range(2, n + 1):
                p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
            dp = n * (x * p1 - p0) / (x * x - 1)
            dx = p1 / dp
            x -= dx
            if abs(dx) < 1e-16:
                break
        xs.append(x)
        ws.append(2 / ((1 - x * x) * dp * dp))
    return xs, ws


GL64 = gauss_legendre(64)


def integrate(f, a, b, pieces=8, rule=GL64):
    xs, ws = rule
    h = (b - a) / pieces
    s = 0.0
    for p in range(pieces):
        lo = a + p * h
        s += sum(w * f(lo + (x + 1) * h / 2) for x, w in zip(xs, ws)) * h / 2
    return s


def cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def norm(a):
    return math.sqrt(dot(a, a))


def sub(a, b):
    return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]


def rk4(f, y, t, h):
    k1 = f(t, y)
    k2 = f(t + h / 2, [a + h / 2 * b for a, b in zip(y, k1)])
    k3 = f(t + h / 2, [a + h / 2 * b for a, b in zip(y, k2)])
    k4 = f(t + h, [a + h * b for a, b in zip(y, k3)])
    return [a + h / 6 * (b + 2 * c + 2 * d + e) for a, b, c, d, e in zip(y, k1, k2, k3, k4)]


print("=" * 72, "\nA. Right rectangular prism (Nagy, Papp & Benedek 2000) and Poisson's relation")


def prism_terms(p, box):
    """Unit-density (G*rho = 1) potential, attraction and Hessian of a prism at point p.
    box = (x1, x2, y1, y2, z1, z2). Returns V, g (grad V), T (Hessian of V)."""
    V = 0.0
    g = [0.0, 0.0, 0.0]
    T = [[0.0] * 3 for _ in range(3)]
    xs, ys, zs = box[0:2], box[2:4], box[4:6]
    for i, X0 in enumerate(xs):
        for j, Y0 in enumerate(ys):
            for k, Z0 in enumerate(zs):
                s = (-1) ** (i + j + k + 1)  # + at the upper corner (i=j=k=1)
                X, Y, Z = X0 - p[0], Y0 - p[1], Z0 - p[2]
                R = math.sqrt(X * X + Y * Y + Z * Z)
                lx, ly, lz = math.log(X + R), math.log(Y + R), math.log(Z + R)
                ax, ay, az = math.atan(Y * Z / (X * R)), math.atan(Z * X / (Y * R)), math.atan(X * Y / (Z * R))
                V += s * (X * Y * lz + Y * Z * lx + Z * X * ly - X * X / 2 * ax - Y * Y / 2 * ay - Z * Z / 2 * az)
                # d/dp = -d/dX of the primitive
                g[0] -= s * (Y * lz + Z * ly - X * ax)
                g[1] -= s * (Z * lx + X * lz - Y * ay)
                g[2] -= s * (X * ly + Y * lx - Z * az)
                T[0][0] -= s * ax
                T[1][1] -= s * ay
                T[2][2] -= s * az
                T[0][1] += s * lz
                T[0][2] += s * ly
                T[1][2] += s * lx
    T[1][0], T[2][0], T[2][1] = T[0][1], T[0][2], T[1][2]
    return V, g, T


# A1: far field of a 1 m cube (rho = 1000) tends to a point mass; a cube has no quadrupole.
rho = 1000.0
cube = (-0.5, 0.5, -0.5, 0.5, -0.5, 0.5)
for dist in (10.0, 100.0, 1000.0):
    p = [dist * 0.6, dist * 0.0, dist * 0.8]
    V, gv, T = prism_terms(p, cube)
    gz = G * rho * gv[2]
    want = -G * rho * 1.0 * p[2] / dist**3  # attraction toward the cube at the origin
    info(f"prism far field at {dist:g} m: relative error of g_z vs point mass", f"{abs(gz / want - 1):.2e}")
V, gv, T = prism_terms([6.0, 0.0, 8.0], cube)
check("prism g_z vs point mass at 10 m (cube: hexadecapole only, ~(a/r)^4)", G * rho * gv[2], -G * rho * 8.0 / 1000.0, rel=1e-4)

# A2: Bouguer slab limit, 2*pi*G*rho*h; the textbook 0.1119 mGal/m for 2670 kg/m3.
h, L = 100.0, 2e6
V, gv, T = prism_terms([0.0, 0.0, 1e-3], (-L, L, -L, L, -h, 0.0))
slab = -G * 2670.0 * gv[2]
check("prism -> Bouguer plate, g per metre of 2670 kg/m3 (mGal/m)", slab / h * 1e5, 0.1119, rel=1e-3)
check("prism -> Bouguer plate vs 2*pi*G*rho*h", slab, 2 * math.pi * G * 2670.0 * h, rel=1e-4)

# A3: Laplace outside (trace 0) and Poisson inside (trace -4 pi G rho).
V, gv, T = prism_terms([0.3, 0.7, 2.1], (-1, 1.2, -0.4, 0.9, -2.0, -0.3))
check("prism tensor trace outside (Laplace)", T[0][0] + T[1][1] + T[2][2], 0.0, abs_=1e-12)
V, gv, T = prism_terms([0.1, 0.2, -0.3], (-1, 1.2, -0.4, 0.9, -2.0, 0.5))
check("prism tensor trace inside / (-4 pi) (Poisson)", (T[0][0] + T[1][1] + T[2][2]) / (-4 * math.pi), 1.0, rel=1e-13)
# The potential at the centre of a unit cube: a closed number, checked by brute-force quadrature.
V0, _, Tc = prism_terms([1e-9, 2e-9, 3e-9], cube)
corner = integrate(lambda x: integrate(lambda y: integrate(lambda z: 1 / math.sqrt(x * x + y * y + z * z), 0, 0.5, 4), 0, 0.5, 4), 0, 0.5, 4)
check("unit-cube potential at its centre (Nagy formula vs 3-D quadrature)", V0, 8 * corner, rel=1e-6)
info("unit-cube centre potential  int 1/r dV", f"{V0:.12f}")

# A4: Poisson's relation, B = (mu0/4pi) * Hessian(int 1/r dV) . M  (B outside; mu0*H inside).
M = [0.0, 0.0, 1e5]  # A/m
for name, pt, vol in (("far field of a magnetised cube vs its dipole", [3.0, 4.0, 12.0], 1.0),):
    _, _, T = prism_terms(pt, cube)
    B = [MU0 / (4 * math.pi) * sum(T[i][j] * M[j] for j in range(3)) for i in range(3)]
    r = norm(pt)
    m = [x * vol for x in M]
    rh = [x / r for x in pt]
    Bd = [MU0 / (4 * math.pi) * (3 * dot(m, rh) * rh[i] - m[i]) / r**3 for i in range(3)]
    check(name + " (B_z)", B[2], Bd[2], rel=1e-4)
_, _, T = prism_terms([1e-9, 2e-9, 3e-9], cube)
check("cube centre: demagnetising factor N_zz = -Hessian_zz/(4 pi) = 1/3", -T[2][2] / (4 * math.pi), 1 / 3, rel=1e-8)

print("=" * 72, "\nA5. Polyhedron (Werner & Scheeres 1997) on the triangulated cube, against the prism")


def polyhedron(p, verts, faces):
    """Unit-density potential and attraction, and the solid-angle sum, of a closed polyhedron
    with outward counter-clockwise triangles."""
    U = 0.0
    g = [0.0, 0.0, 0.0]
    wsum = 0.0
    edges = {}
    for f in faces:
        a, b, c = (verts[i] for i in f)
        n = cross(sub(b, a), sub(c, a))
        nn = norm(n)
        n = [x / nn for x in n]
        r1, r2, r3 = sub(a, p), sub(b, p), sub(c, p)
        l1, l2, l3 = norm(r1), norm(r2), norm(r3)
        w = 2 * math.atan2(dot(r1, cross(r2, r3)), l1 * l2 * l3 + l1 * dot(r2, r3) + l2 * dot(r3, r1) + l3 * dot(r1, r2))
        wsum += w
        rf = dot(n, r1)
        U -= 0.5 * rf * rf * w
        for k in range(3):
            g[k] += n[k] * rf * w
        for i0, i1 in ((f[0], f[1]), (f[1], f[2]), (f[2], f[0])):
            e = sub(verts[i1], verts[i0])
            ne = cross(e, n)
            le = norm(ne)
            ne = [x / le for x in ne]
            edges.setdefault((min(i0, i1), max(i0, i1)), []).append((n, ne))
    for (i0, i1), pair in edges.items():
        ra, rb = sub(verts[i0], p), sub(verts[i1], p)
        la, lb, le = norm(ra), norm(rb), norm(sub(verts[i1], verts[i0]))
        Le = math.log((la + lb + le) / (la + lb - le))
        # E_e = n_A ne_A^T + n_B ne_B^T ; contribution r.E.r L and -E r L
        Er = [sum(nA[i] * neA[j] * ra[j] for (nA, neA) in pair for j in range(3)) for i in range(3)]
        U += 0.5 * dot(ra, Er) * Le
        for k in range(3):
            g[k] -= Er[k] * Le
    return U, g, wsum


cv = [[x, y, z] for x in (-0.5, 0.5) for y in (-0.5, 0.5) for z in (-0.5, 0.5)]
idx = {tuple(v): i for i, v in enumerate(cv)}


def q(x, y, z):
    return idx[(x, y, z)]


h_ = 0.5
quads = [  # outward CCW
    [q(h_, -h_, -h_), q(h_, h_, -h_), q(h_, h_, h_), q(h_, -h_, h_)],
    [q(-h_, -h_, -h_), q(-h_, -h_, h_), q(-h_, h_, h_), q(-h_, h_, -h_)],
    [q(-h_, h_, -h_), q(-h_, h_, h_), q(h_, h_, h_), q(h_, h_, -h_)],
    [q(-h_, -h_, -h_), q(h_, -h_, -h_), q(h_, -h_, h_), q(-h_, -h_, h_)],
    [q(-h_, -h_, h_), q(h_, -h_, h_), q(h_, h_, h_), q(-h_, h_, h_)],
    [q(-h_, -h_, -h_), q(-h_, h_, -h_), q(h_, h_, -h_), q(h_, -h_, -h_)],
]
tris = [t for a, b, c, d in quads for t in ((a, b, c), (a, c, d))]
pt = [0.9, -0.4, 1.3]
Up, gp, wout = polyhedron(pt, cv, tris)
Vq, gq, _ = prism_terms(pt, cube)
check("polyhedron potential = prism potential (outside)", Up, Vq, rel=1e-12)
check("polyhedron g_z = prism g_z", gp[2], gq[2], rel=1e-12)
check("solid-angle sum outside = 0", wout, 0.0, abs_=1e-12)
_, _, win = polyhedron([0.1, 0.2, -0.05], cv, tris)
check("solid-angle sum inside = 4 pi (inside/outside test, Laplacian = -G rho * sum)", win, 4 * math.pi, rel=1e-13)

print("=" * 72, "\nB. Boris pusher (Boris 1970; Qin et al. 2013) against RK4")


def boris_step(x, v, qm, E, B, dt, exact_phase=False):
    vm = [v[i] + qm * E[i] * dt / 2 for i in range(3)]
    t = [qm * B[i] * dt / 2 for i in range(3)]
    if exact_phase:  # Zenitani & Umeda 2018: t -> tan(Omega dt/2) b
        tm = norm(t)
        if tm > 0:
            t = [ti * math.tan(tm) / tm for ti in t]
    t2 = dot(t, t)
    s = [2 * ti / (1 + t2) for ti in t]
    vp = [a + b for a, b in zip(vm, cross(vm, t))]
    vplus = [a + b for a, b in zip(vm, cross(vp, s))]
    vn = [vplus[i] + qm * E[i] * dt / 2 for i in range(3)]
    return [x[i] + vn[i] * dt for i in range(3)], vn


# B1: uniform B, Omega = 1, Omega*dt = 0.1, 10^5 steps (~1592 gyro-periods).
dt, nsteps = 0.1, 100_000
x, v = [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]
for _ in range(nsteps):
    x, v = boris_step(x, v, 1.0, [0, 0, 0], [0, 0, 1.0], dt)
check("Boris |v| after 1e5 steps in uniform B (exact energy)", norm(v), 1.0, abs_=1e-13)


def lorentz_rhs(t, y):
    vx, vy, vz = y[3:]
    return [vx, vy, vz, vy * 1.0, -vx * 1.0, 0.0]


y = [1.0, 0.0, 0.0, 0.0, -1.0, 0.0]
for _ in range(nsteps):
    y = rk4(lorentz_rhs, y, 0, dt)
erk = abs(math.hypot(y[3], y[4]) ** 2 - 1)
info("RK4 energy error after 1e5 steps, Omega dt = 0.1", f"{erk:.3e} (predicted 1e5*(h^6/72) = {1e5 * dt**6 / 72:.3e})")
check("RK4 secular energy loss = N h^6/72 (amplification |R(ih)|^2 = 1 - h^6/72 + ...)", erk, 1e5 * dt**6 / 72, rel=2e-2)
# B2: gyrophase: one Boris step rotates by 2 atan(Omega dt/2), not Omega dt.
x, v = [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]
x, v = boris_step(x, v, 1.0, [0, 0, 0], [0, 0, 1.0], dt)
check("Boris rotation per step = 2 atan(Omega dt / 2)", math.atan2(-v[0], -v[1]), 2 * math.atan(dt / 2), rel=1e-14)
x, v = boris_step([1.0, 0.0, 0.0], [0.0, -1.0, 0.0], 1.0, [0, 0, 0], [0, 0, 1.0], dt, exact_phase=True)
check("Zenitani-Umeda corrected rotation per step = Omega dt", math.atan2(-v[0], -v[1]), dt, rel=1e-14)
info("Boris phase lag per gyro-period at Omega dt = 0.1 (rad)", f"{2 * math.pi / dt * (dt - 2 * math.atan(dt / 2)):.4e}  (~ 2 pi (Omega dt)^2 / 12 = {2 * math.pi * dt**2 / 12:.4e})")
# B3: E x B drift reproduced exactly (non-relativistic Boris).
x, v = [0.0, 0.0, 0.0], [0.0, -0.1, 0.0]
for _ in range(10_000):
    x, v = boris_step(x, v, 1.0, [0.1, 0.0, 0.0], [0, 0, 1.0], dt)
check("Boris: a particle at the E x B velocity stays there (10^4 steps)", abs(v[0]) + abs(v[1] + 0.1), 0.0, abs_=1e-14)
x, v = [0.0, 0.0, 0.0], [0.3, 0.0, 0.0]
ysum, steps = 0.0, 0
Ef = [0.1, 0.0, 0.0]
xs0 = x
for k in range(62_832):  # ~1000 gyro-periods at dt = 0.1
    x, v = boris_step(x, v, 1.0, Ef, [0, 0, 1.0], dt)
info("Boris E x B drift: mean v_y over ~1000 periods", f"{x[1] / (62_832 * dt):.6f} (exact -E/B = -0.1)")
check("Boris E x B drift with gyration (mean v_y; residual = gyro-orbit end offset / time)", x[1] / (62_832 * dt), -0.1, rel=2e-4)

# B4: a 100 keV proton in Earth's dipole at L = 3, equatorial pitch angle 30 degrees: bounce
# period from Boris against the guiding-centre integral T(alpha) (Schulz & Lanzerotti 1974).
B0, RE = 3.12e-5, 6.371e6


def dipole(r):
    rr = norm(r)
    rh = [c / rr for c in r]
    f = B0 * (RE / rr) ** 3
    return [f * (3 * rh[2] * rh[i] - (1.0 if i == 2 else 0.0)) for i in range(3)]


def bounce_T(alpha):
    s2 = math.sin(alpha) ** 2
    lo, hi = 0.0, math.pi / 2
    for _ in range(200):
        mid = (lo + hi) / 2
        if math.cos(mid) ** 6 - s2 * math.sqrt(1 + 3 * math.sin(mid) ** 2) > 0:
            lo = mid
        else:
            hi = mid
    lm = lo

    def f(ph):
        lam = lm * math.sin(ph)
        sl = math.sin(lam)
        val = 1 - s2 * math.sqrt(1 + 3 * sl * sl) / math.cos(lam) ** 6
        return math.cos(lam) * math.sqrt(1 + 3 * sl * sl) / math.sqrt(max(val, 1e-300)) * lm * math.cos(ph)

    return integrate(f, 0.0, math.pi / 2, 16)


check("T(alpha -> 0) = 1 + ln(2+sqrt3)/(2 sqrt3)", bounce_T(1e-7), 1 + math.log(2 + math.sqrt(3)) / (2 * math.sqrt(3)), rel=1e-4)
check("T(alpha -> 90 deg) = pi/(3 sqrt2), the small-oscillation limit", bounce_T(math.pi / 2 - 1e-3), math.pi / (3 * math.sqrt(2)), rel=1e-4)
Lsh, alpha, ek = 3.0, math.radians(30.0), 100e3 * QE
vmag = math.sqrt(2 * ek / MP)
qm = QE / MP
x = [Lsh * RE, 0.0, 0.0]
v = [0.0, vmag * math.sin(alpha), -vmag * math.cos(alpha)]
Beq = B0 / Lsh**3
dt = 2 * math.pi / (qm * Beq / math.sin(alpha) ** 2) / 60
t, crossings = 0.0, []
prevz = x[2]
while len(crossings) < 4:
    x, v = boris_step(x, v, qm, [0, 0, 0], dipole(x), dt)
    t += dt
    if prevz < 0 <= x[2]:
        crossings.append(t - dt * x[2] / (x[2] - prevz))
    prevz = x[2]
tb_boris = (crossings[-1] - crossings[0]) / (len(crossings) - 1)
tb_gc = 4 * Lsh * RE / vmag * bounce_T(alpha)
info("dipole bounce period, Boris vs guiding centre (s)", f"{tb_boris:.4f} vs {tb_gc:.4f}")
check("Boris bounce period in a dipole vs guiding-centre T(alpha)", tb_boris, tb_gc, rel=3e-3)
check("Boris keeps |v| in the dipole (energy, magnetic field only)", norm(v), vmag, rel=1e-12)

print("=" * 72, "\nC. Kepler (Hairer, Lubich & Wanner, e = 0.6) and Mercury's perihelion")
e = 0.6


def kepler_acc(qv):
    r3 = (qv[0] ** 2 + qv[1] ** 2) ** 1.5
    return [-qv[0] / r3, -qv[1] / r3]


def H(qv, pv):
    return 0.5 * (pv[0] ** 2 + pv[1] ** 2) - 1 / math.hypot(*qv)


def run(method, h, periods):
    qv, pv = [1 - e, 0.0], [0.0, math.sqrt((1 + e) / (1 - e))]
    n = round(periods * 2 * math.pi / h)
    maxdh = 0.0
    if method == "verlet":
        a = kepler_acc(qv)
        for _ in range(n):
            ph = [pv[i] + h / 2 * a[i] for i in range(2)]
            qv = [qv[i] + h * ph[i] for i in range(2)]
            a = kepler_acc(qv)
            pv = [ph[i] + h / 2 * a[i] for i in range(2)]
            maxdh = max(maxdh, abs(H(qv, pv) + 0.5))
    elif method == "yoshida4":
        w1 = 1 / (2 - 2 ** (1 / 3))
        w0 = -(2 ** (1 / 3)) * w1
        for _ in range(n):
            for w in (w1, w0, w1):
                a = kepler_acc(qv)
                ph = [pv[i] + w * h / 2 * a[i] for i in range(2)]
                qv = [qv[i] + w * h * ph[i] for i in range(2)]
                a = kepler_acc(qv)
                pv = [ph[i] + w * h / 2 * a[i] for i in range(2)]
            maxdh = max(maxdh, abs(H(qv, pv) + 0.5))
    else:
        y = qv + pv
        f = lambda t, y: [y[2], y[3]] + kepler_acc(y[:2])
        for _ in range(n):
            y = rk4(f, y, 0, h)
            maxdh = max(maxdh, abs(H(y[:2], y[2:]) + 0.5))
        qv, pv = y[:2], y[2:]
    Lz = qv[0] * pv[1] - qv[1] * pv[0]
    return maxdh, abs(H(qv, pv) + 0.5), abs(Lz - math.sqrt(1 - e * e)), math.hypot(qv[0] - (1 - e), qv[1])


for m_, hh in (("verlet", 0.01), ("yoshida4", 0.01), ("rk4", 0.01)):
    for per in (10, 100):
        mx, fin, dl, dq = run(m_, hh, per)
        info(f"Kepler e=0.6 {m_:8s} h={hh} {per:4d} periods", f"max|dH| {mx:.2e}  final|dH| {fin:.2e}  |dL| {dl:.1e}  |q-q0| {dq:.2e}")
mx10, fin10, _, _ = run("rk4", 0.01, 10)
mx100, fin100, _, _ = run("rk4", 0.01, 100)
check("RK4 energy error grows ~linearly (x10 periods -> ~x10)", fin100 / fin10, 10.0, rel=0.2)
mxv10, _, _, _ = run("verlet", 0.01, 10)
mxv100, _, dLv, _ = run("verlet", 0.01, 100)
check("Stormer-Verlet energy error stays bounded (max over 100 periods / over 10)", mxv100 / mxv10, 1.0, rel=0.05)
check("Stormer-Verlet conserves angular momentum to rounding", dLv, 0.0, abs_=1e-11)

# Mercury: integrate the Schwarzschild orbit equation u'' + u = mu/h^2 + 3 mu u^2 / c^2 for one
# radial period; the extra angle is 6 pi mu / (c^2 a (1 - e^2)), 42.98"/century.
a_m, e_m, P_m = 5.7909050e10, 0.205630, 87.9691
p_m = a_m * (1 - e_m**2)
u0 = 1 / (a_m * (1 - e_m))
hh2 = GM_SUN * p_m  # Newtonian h^2 (the GR shift at this order is independent of the h choice)
f = lambda phi, y: [y[1], GM_SUN / hh2 - y[0] + 3 * GM_SUN * y[0] ** 2 / C**2]
y, phi, dphi = [u0, 0.0], 0.0, 2 * math.pi / 200_000
prev = y
while True:
    yn = rk4(f, y, phi, dphi)
    if phi > math.pi and y[1] > 0 >= yn[1]:
        frac = y[1] / (y[1] - yn[1])
        lo_, hi_ = 0.0, dphi
        for _ in range(60):  # bisection on u' = 0 inside the last step
            mid = (lo_ + hi_) / 2
            if rk4(f, y, phi, mid)[1] > 0:
                lo_ = mid
            else:
                hi_ = mid
        phi_peri = phi + lo_
        break
    y, phi = yn, phi + dphi
shift = phi_peri - 2 * math.pi
century = shift * 36525 / P_m * RAD2AS
check("Mercury perihelion advance from the orbit equation (arcsec/century)", century, 42.98, rel=2e-3)
check("... equals 6 pi GM/(c^2 a (1-e^2)) per orbit", shift, 6 * math.pi * GM_SUN / (C**2 * p_m), rel=1e-4)

print("=" * 72, "\nD. Rays: one Hamiltonian integrator for curved spacetime and plasma")


def trace_ray(n2_and_grad, x, k, s_end, h_of):
    """Integrate dx/ds = k, dk/ds = grad(n^2)/2 (H = (k^2 - n^2)/2 = 0) with RK4."""
    y = x + k
    s = 0.0

    def rhs(_, y):
        _, gr = n2_and_grad(y[:3])
        return [y[3], y[4], y[5], gr[0] / 2, gr[1] / 2, gr[2] / 2]

    while s < s_end:
        hstep = min(h_of(y), s_end - s)
        y = rk4(rhs, y, s, hstep)
        s += hstep
    return y


rs = 2 * GM_SUN / C**2 / R_SUN  # in solar radii


def gordon_isotropic(p):  # Schwarzschild in isotropic coordinates: n = (1 + rs/4r)^3 / (1 - rs/4r)
    r = norm(p)
    a = rs / (4 * r)
    n = (1 + a) ** 3 / (1 - a)
    dn_dr = (-(a / r)) * (3 * (1 + a) ** 2 * (1 - a) + (1 + a) ** 3) / (1 - a) ** 2
    return n * n, [2 * n * dn_dr * c / r for c in p]


Lr = 2e4
y = trace_ray(gordon_isotropic, [-Lr, 1.0, 0.0], [gordon_isotropic([-Lr, 1.0, 0.0])[0] ** 0.5, 0.0, 0.0], 2 * Lr,
              lambda y: max(0.02, 0.02 * abs(y[0])))
defl = -math.atan2(y[4], y[3])
# finite path: add the two missing tails 2 * (2GM/c^2 b) * (1 - L/sqrt(L^2 + b^2))
tails = 2 * rs * (1 - Lr / math.hypot(Lr, 1.0))
check("ray through Gordon index of the Sun: deflection at the limb (arcsec)", (defl + tails) * RAD2AS, 1.7512, rel=2e-4)


def linear_plasma(p):  # omega_p^2/omega^2 = z/H for z > 0, H = 1
    z = p[2]
    return (1 - max(z, 0.0)), [0.0, 0.0, -1.0 if z > 0 else 0.0]


th0 = math.radians(40)
y = trace_ray(linear_plasma, [0.0, 0.0, 1e-12], [math.sin(th0), 0.0, math.cos(th0)], 4 * math.cos(th0), lambda y: 1e-3)
check("ray in a linear plasma layer: ground range 2H sin(2 theta0) (projectile analogue)", y[0], 2 * math.sin(2 * th0), rel=1e-6)

k_tec = QE**2 / (8 * math.pi**2 * EPS0 * ME)
check("ionospheric constant e^2/(8 pi^2 eps0 m_e) (m^3/s^2)", k_tec, 40.3, rel=2e-3)
info("ionospheric constant, unrounded", f"{k_tec:.4f}")
check("L1 range delay per TECU (m)", k_tec * 1e16 / 1575.42e6**2, 0.1624, rel=1e-3)
# exact plasma index integrated through a Chapman layer vs 40.3 TEC/f^2
NmF2, hm, Hs = 1e12, 300e3, 60e3


def chap(hh):
    zz = (hh - hm) / Hs
    return NmF2 * math.exp(0.5 * (1 - zz - math.exp(-zz)))


f1 = 1575.42e6
fp2 = lambda n_: n_ * QE**2 / (EPS0 * ME) / (2 * math.pi) ** 2
phase_adv = integrate(lambda hh: 1 - math.sqrt(1 - fp2(chap(hh)) / f1**2), 0.0, 2000e3, 40)
group_del = integrate(lambda hh: 1 / math.sqrt(1 - fp2(chap(hh)) / f1**2) - 1, 0.0, 2000e3, 40)
tec = integrate(chap, 0.0, 2000e3, 40)
check("Chapman layer: phase advance = 40.3 TEC/f^2 to first order", phase_adv, k_tec * tec / f1**2, rel=1e-4)
check("Chapman layer: group delay = +40.3 TEC/f^2 to first order", group_del, k_tec * tec / f1**2, rel=1e-4)
info("Chapman layer TEC (TECU) and L1 group delay (m)", f"{tec / 1e16:.2f} TECU, {group_del:.3f} m")
# parabolic layer: group height h' = h_b + (y_m/2)(f/f_c) ln((1+f/f_c)/(1-f/f_c))
hb, ym = 200e3, 100e3
for ratio in (0.5, 0.9):
    zr = ym * (1 - math.sqrt(1 - ratio**2))  # reflection height above base
    # group index 1/sqrt(1 - X), X = ratio^-2 * (1 - ((z - ym)/ym)^2) ... substitute z = zr sin^2-like
    X = lambda z: (1 - ((z - ym) / ym) ** 2) / ratio**2
    num = integrate(lambda t: 1 / math.sqrt(max(1 - X(zr * (1 - t * t)), 1e-300)) * 2 * zr * t, 0.0, 1.0, 16)
    want = ym / 2 * ratio * math.log((1 + ratio) / (1 - ratio))
    check(f"parabolic layer virtual height above base at f/fc = {ratio}", num, want, rel=1e-6)

# Strong deflection: exact Schwarzschild (Darwin) against Bozza 2002, units G = c = M = 1.
def darwin(b):
    lo, hi = 0.0, 1 / 3  # closest-approach u0 solves 1/b^2 = u0^2 - 2 u0^3
    for _ in range(200):
        mid = (lo + hi) / 2
        if mid * mid - 2 * mid**3 < 1 / b**2:
            lo = mid
        else:
            hi = mid
    u0 = lo
    # 1/b^2 - u^2 + 2u^3 = u0 t^2 Q(u) with u = u0(1 - t^2): divide out the root analytically
    Q = lambda u: (u + u0) - 2 * (u * u + u * u0 + u0 * u0)
    g = lambda t: 2 * math.sqrt(u0) / math.sqrt(Q(u0 * (1 - t * t)))
    # graded pieces toward t = 0 where the near-double root makes the integrand peak
    edges = [0.0] + [10.0 ** (-k) for k in range(12, 0, -1)] + [1.0]
    tot = sum(integrate(g, a_, b_, 2) for a_, b_ in zip(edges, edges[1:]))
    return 2 * tot - math.pi


bc = 3 * math.sqrt(3)
bbar = math.log(216 * (7 - 4 * math.sqrt(3))) - math.pi
check("Bozza b-bar = ln[216(7-4 sqrt3)] - pi", bbar, -0.4002, abs_=1e-4)
for eps in (1e-2, 1e-3, 1e-4):
    ex = darwin(bc * (1 + eps))
    info(f"Schwarzschild deflection at b/bc - 1 = {eps:g}", f"exact {ex:.6f} rad, Bozza {-math.log(eps) + bbar:.6f}, diff {ex + math.log(eps) - bbar:+.2e}")
check("exact deflection -> Bozza strong-field limit at b/bc-1 = 1e-4", darwin(bc * 1.0001), -math.log(1e-4) + bbar, abs_=2e-3)
bw = 2e4
check("exact deflection weak field: 4/b + 15 pi/(4 b^2) (the crate's light_bending)", darwin(bw), 4 / bw + 15 * math.pi / (4 * bw * bw), rel=1e-7)
theta_m87 = 2 * math.sqrt(27) * G * 6.5e9 * M_SUN / (C**2 * 16.8e6 * PC) * RAD2AS * 1e6
info("M87* shadow diameter for M = 6.5e9 Msun, D = 16.8 Mpc (uas)", f"{theta_m87:.1f} (EHT ring 42 +/- 3)")
theta_sgr = 2 * math.sqrt(27) * G * 4.0e6 * M_SUN / (C**2 * 8.0e3 * PC) * RAD2AS * 1e6
info("Sgr A* shadow diameter for M = 4.0e6 Msun, D = 8.0 kpc (uas)", f"{theta_sgr:.1f} (EHT 2022 shadow 48.7 +/- 7.0, ring 51.8 +/- 2.3)")
rph = lambda a, sgn: 2 * (1 + math.cos(2 / 3 * math.acos(-sgn * a)))
check("Kerr photon orbit, a = 0: 3M", rph(0.0, 1), 3.0, rel=1e-14)
check("Kerr photon orbit, a = M prograde: M", rph(1.0, 1), 1.0, rel=1e-14)
check("Kerr photon orbit, a = M retrograde: 4M", rph(1.0, -1), 4.0, rel=1e-14)
shapiro = 2 * GM_SUN / C**3 * math.log(4 * AU * 1.524 * AU / R_SUN**2) * 2
info("Shapiro round-trip delay Earth-Mars at superior conjunction, grazing (us)", f"{shapiro * 1e6:.1f}")

print("=" * 72, "\nE. Special functions")


# E1: Faddeeva w(z) = exp(-z^2) erfc(-iz): Taylor series for |z| < 3, Laplace continued fraction beyond.
# (A check only. In the crate use Weideman 1994 or Zaghloul & Ali 2011, uniform to ~1e-14.)
def faddeeva(z):
    if abs(z) < 3:  # the series loses ~|z|^2/ln(10) digits: 4 at |z| = 3
        s, term, n = 0j, 1 + 0j, 0
        # w(z) = sum (iz)^n / Gamma(n/2 + 1)
        while True:
            t = (1j * z) ** n / math.gamma(n / 2 + 1)
            s += t
            if n > 10 and abs(t) < 1e-17 * abs(s):
                return s
            n += 1
    if z.imag < 0:  # reflection w(z) = 2 exp(-z^2) - w(-z)
        return 2 * cmath.exp(-z * z) - faddeeva(-z)
    cf = 0j
    for k in range(400, 0, -1):
        cf = (k / 2) / (z - cf)
    return 1j / math.sqrt(math.pi) / (z - cf)


check("w(i) = e erfc(1)", faddeeva(1j).real, math.e * math.erfc(1.0), rel=1e-13)
check("Re w(1) = exp(-1)", faddeeva(1 + 0j).real, math.exp(-1), rel=1e-13)
check("Im w(1) = (2/sqrt pi) Dawson(1), F(1) = 0.5380795069127684", faddeeva(1 + 0j).imag * math.sqrt(math.pi) / 2, 0.5380795069127684, rel=1e-12)
zz = 1.79 + 2.39j  # |z| just under 3, where both methods are valid
cfv = 0j
for k in range(400, 0, -1):
    cfv = (k / 2) / (zz - cfv)
check("w by series vs continued fraction at 1.79+2.39i", abs(faddeeva(zz) - 1j / math.sqrt(math.pi) / (zz - cfv)), 0.0, abs_=1e-11)
Zf = lambda zeta: 1j * math.sqrt(math.pi) * faddeeva(zeta)


def langmuir_root(kld, w0):
    """Root of 1 + (1 + zeta Z(zeta))/(k lambda_D)^2 = 0, zeta = omega/(sqrt2 k v_th), omega/omega_p."""
    D = lambda w: 1 + (1 + (w / (math.sqrt(2) * kld)) * Zf(w / (math.sqrt(2) * kld))) / kld**2
    w = w0
    for _ in range(60):
        dw = 1e-7
        step = D(w) / ((D(w + dw) - D(w - dw)) / (2 * dw))
        w -= step
        if abs(step) < 1e-14:
            break
    return w


wl = langmuir_root(0.5, 1.4 - 0.15j)
info("Langmuir wave root at k lambda_D = 0.5 (omega/omega_p)", f"{wl.real:.4f} {wl.imag:+.4f}i")
check("Landau damping root, real part 1.4156", wl.real, 1.4156, abs_=1e-4)
check("Landau damping root, damping -0.1534 (commonly quoted 0.1533)", wl.imag, -0.15336, abs_=2e-4)
# Kinetic Jeans instability: same Z, opposite sign: k^2/kJ^2 = 1 - sqrt(pi) y e^{y^2} erfc(y), gamma = sqrt2 k sigma y
target = 0.25
lo, hi = 0.0, 10.0
for _ in range(200):
    mid = (lo + hi) / 2
    if 1 - math.sqrt(math.pi) * mid * math.exp(mid * mid) * math.erfc(mid) > target:
        lo = mid
    else:
        hi = mid
yJ = lo
Dj = 1 - (1 / target) * (1 + (1j * yJ) * Zf(1j * yJ))
check("kinetic Jeans root at k = kJ/2 satisfies the Landau dispersion with G for -e^2", abs(Dj), 0.0, abs_=1e-10)
info("kinetic Jeans growth at k = kJ/2: gamma/(sqrt2 k sigma)", f"{yJ:.6f}")


# E2: Carlson symmetric integrals (duplication), checked against Carlson (1995) Section 3.
def RF(x, y, z):
    for _ in range(200):
        lam = math.sqrt(x * y) + math.sqrt(y * z) + math.sqrt(z * x)
        x, y, z = (x + lam) / 4, (y + lam) / 4, (z + lam) / 4
        A = (x + y + z) / 3
        if max(abs(A - x), abs(A - y), abs(A - z)) < 1e-4 * abs(A):
            break
    X, Y = 1 - x / A, 1 - y / A
    Z = -X - Y
    E2, E3 = X * Y - Z * Z, X * Y * Z
    return (1 - E2 / 10 + E3 / 14 + E2 * E2 / 24 - 3 * E2 * E3 / 44) / math.sqrt(A)


def RD(x, y, z):
    s, fac = 0.0, 1.0
    for _ in range(200):
        lam = math.sqrt(x * y) + math.sqrt(y * z) + math.sqrt(z * x)
        s += fac / (math.sqrt(z) * (z + lam))
        fac /= 4
        x, y, z = (x + lam) / 4, (y + lam) / 4, (z + lam) / 4
        A = (x + y + 3 * z) / 5
        if max(abs(A - x), abs(A - y), abs(A - z)) < 1e-4 * abs(A):
            break
    X, Y = 1 - x / A, 1 - y / A
    Z = -(X + Y) / 3
    E2 = X * Y - 6 * Z * Z
    E3 = (3 * X * Y - 8 * Z * Z) * Z
    E4 = 3 * (X * Y - Z * Z) * Z * Z
    E5 = X * Y * Z**3
    return 3 * s + fac * (1 - 3 * E2 / 14 + E3 / 6 + 9 * E2 * E2 / 88 - 3 * E4 / 22 - 9 * E2 * E3 / 52 + 3 * E5 / 26) / (A * math.sqrt(A))


check("RF(1,2,0) (Carlson 1995)", RF(1, 2, 0), 1.3110287771461, rel=1e-12)
check("RF(2,3,4)", RF(2, 3, 4), 0.58408284167715, rel=1e-12)
check("RD(0,2,1)", RD(0, 2, 1), 1.7972103521034, rel=1e-12)
check("RD(2,3,4)", RD(2, 3, 4), 0.16510527294261, rel=1e-12)
K = lambda m: RF(0, 1 - m, 1)
Ee = lambda m: RF(0, 1 - m, 1) - m / 3 * RD(0, 1 - m, 1)
check("K(0.5) = RF(0, 1-m, 1) = 1.854074677301372 (A&S 17.1; the electrodes AGM test)", K(0.5), 1.854074677301372, rel=1e-14)
check("E(0.5) = RF - (m/3) RD = 1.350643881047675", Ee(0.5), 1.350643881047675, rel=1e-14)


def cel(kc, p, a, b):  # Bulirsch 1969, in the form of Derby & Olbert 2010
    k, em = abs(kc), 1.0
    if p > 0:
        pp = math.sqrt(p)
        bb = b / pp
        aa = a
    else:
        f = kc * kc
        q_ = 1 - f
        g = 1 - p
        f = f - p
        q_ = q_ * (b - a * p)
        pp = math.sqrt(f / g)
        aa = (a - b) / g
        bb = -q_ / (g * g * pp) + aa * pp
    f = aa
    aa = aa + bb / pp
    g = k / pp
    bb = 2 * (bb + f * g)
    pp = g + pp
    g = em
    em = k + em
    kk = k
    while abs(g - k) > g * 1e-12:
        k = 2 * math.sqrt(kk)
        kk = k * em
        f = aa
        aa = aa + bb / pp
        g = kk / pp
        bb = 2 * (bb + f * g)
        pp = g + pp
        g = em
        em = k + em
    return math.pi / 2 * (bb + aa * em) / (em * (em + pp))


kc = math.sqrt(0.5)
check("cel(kc,1,1,1) = K", cel(kc, 1, 1, 1), K(0.5), rel=1e-13)
check("cel(kc,1,1,kc^2) = E", cel(kc, 1, 1, kc * kc), Ee(0.5), rel=1e-13)


# Current loop off axis (K, E) against Biot-Savart summation, and the on-axis formula.
def loop_B(I, a, rho, z):
    s = (a + rho) ** 2 + z * z
    d = (a - rho) ** 2 + z * z
    m = 4 * a * rho / s
    k_, e_ = K(m), Ee(m)
    c0 = MU0 * I / (2 * math.pi * math.sqrt(s))
    bz = c0 * (k_ + (a * a - rho * rho - z * z) / d * e_)
    br = c0 * z / rho * (-k_ + (a * a + rho * rho + z * z) / d * e_) if rho > 0 else 0.0
    return br, bz


def loop_biot(I, a, p, n=20000):
    B = [0.0, 0.0, 0.0]
    for i in range(n):
        t = 2 * math.pi * (i + 0.5) / n
        src = [a * math.cos(t), a * math.sin(t), 0.0]
        dl = [-a * math.sin(t) * 2 * math.pi / n, a * math.cos(t) * 2 * math.pi / n, 0.0]
        r = sub(p, src)
        rr = norm(r)
        c_ = cross(dl, r)
        for k in range(3):
            B[k] += MU0 * I / (4 * math.pi) * c_[k] / rr**3
    return B


br, bz = loop_B(1.0, 0.04, 0.023, 0.037)
bb_ = loop_biot(1.0, 0.04, [0.023, 0.0, 0.037])
check("current loop B_z (K,E) vs Biot-Savart", bz, bb_[2], rel=1e-10)
check("current loop B_rho (K,E) vs Biot-Savart", br, bb_[0], rel=1e-10)
check("current loop on axis = magnetics::loop_bz", loop_B(1.0, 0.04, 0.0, 0.05)[1], MU0 * 0.04**2 / (2 * (0.04**2 + 0.05**2) ** 1.5), rel=1e-14)
# Near the axis the textbook B_rho cancels catastrophically (Ortner et al. 2023).
rho_s = 1e-7 * 0.04
br_text = loop_B(1.0, 0.04, rho_s, 0.02)[0]
br_series = 3 * MU0 * 1.0 * 0.04**2 * 0.02 * rho_s / (4 * (0.04**2 + 0.02**2) ** 2.5)  # leading term: -(rho/2) dBz/dz
info("loop B_rho at rho = 1e-7 a: textbook K,E formula vs near-axis series", f"{br_text:.6e} vs {br_series:.6e} (relative error {abs(br_text / br_series - 1):.1e})")


def solenoid_DO(n_turns_per_m, I, a, half, rho, z):  # Derby & Olbert 2010, via cel
    B0 = MU0 * n_turns_per_m * I / math.pi
    zp, zm = z + half, z - half
    out = []
    for zs in (zp, zm):
        den = math.sqrt(zs * zs + (rho + a) ** 2)
        al, be = a / den, zs / den
        ga = (a - rho) / (a + rho)
        kk = math.sqrt((zs * zs + (a - rho) ** 2) / (zs * zs + (a + rho) ** 2))
        out.append((al * cel(kk, 1, 1, -1), be * cel(kk, ga * ga, 1, ga)))
    return B0 * (out[0][0] - out[1][0]), B0 * a / (a + rho) * (out[0][1] - out[1][1])


nt, Ls_, aa_ = 1000.0, 0.2, 0.03
bz_ax = MU0 * nt * 1.0 / 2 * ((0.05 + 0.1) / math.hypot(0.05 + 0.1, aa_) - (0.05 - 0.1) / math.hypot(0.05 - 0.1, aa_))
check("Derby-Olbert solenoid on axis (rho -> 0) = magnetics::solenoid_bz", solenoid_DO(nt, 1.0, aa_, 0.1, 1e-12, 0.05)[1], bz_ax, rel=1e-9)
# off axis against a stack of loops (midpoint rule over the length)
nl = 4000
sb_r = sum(loop_B(nt * Ls_ / nl, aa_, 0.02, 0.05 - (-0.1 + (i + 0.5) * Ls_ / nl))[0] for i in range(nl))
sb_z = sum(loop_B(nt * Ls_ / nl, aa_, 0.02, 0.05 - (-0.1 + (i + 0.5) * Ls_ / nl))[1] for i in range(nl))
do_r, do_z = solenoid_DO(nt, 1.0, aa_, 0.1, 0.02, 0.05)
check("Derby-Olbert solenoid off axis B_z vs summed loops", do_z, sb_z, rel=1e-6)
check("Derby-Olbert solenoid off axis B_rho vs summed loops", do_r, sb_r, rel=1e-5)


# E3: Bessel. plasmoid::relaxation::bessel_j uses Simpson on Bessel's integral (2000 panels).
def bessel_simpson(n, x, steps=2000):
    h = math.pi / steps
    f = lambda t: math.cos(n * t - x * math.sin(t))
    inner = sum((4 if k % 2 else 2) * f(k * h) for k in range(1, steps))
    return (f(0) + inner + f(math.pi)) * h / 3 / math.pi


check("J0(1) (DLMF/A&S 9.1) via the crate's method", bessel_simpson(0, 1.0), 0.7651976865579666, rel=1e-14)
check("J1(1)", bessel_simpson(1, 1.0), 0.4400505857449335, rel=1e-14)
check("J0 at its first zero 2.404825557695773", bessel_simpson(0, 2.404825557695773), 0.0, abs_=1e-15)
info("crate-method J0(50) vs A&S 0.05581232766925182", f"{bessel_simpson(0, 50.0):.16f}")

# E4: Gaussian smoothing weights (Jekeli 1981 / Wahr et al. 1998) are W_n = i_n(b)/i_0(b),
# modified spherical Bessel functions: forward recursion is unstable; backward is not.
R_E = 6_378_136.3


def weights_forward(radius, nmax, floor=True):
    b = math.log(2) / (1 - math.cos(radius / R_E))
    e2 = math.exp(-2 * b)
    w = [1.0, (1 + e2) / (1 - e2) - 1 / b]
    for n in range(1, nmax):
        nxt = -(2 * n + 1) / b * w[n] + w[n - 1]
        if floor and (w[n] <= 1e-10 or not (1e-10 <= nxt <= w[n])):
            nxt = 1e-10
        w.append(nxt)
    return w[: nmax + 1]


def weights_backward(radius, nmax):
    b = math.log(2) / (1 - math.cos(radius / R_E))
    rho_ = 0.0
    ratios = [0.0] * (nmax + 1)
    for n in range(nmax + 400, 0, -1):  # i_n/i_{n-1} = 1/((2n+1)/b + i_{n+1}/i_n)
        rho_ = 1 / ((2 * n + 1) / b + rho_)
        if n <= nmax:
            ratios[n] = rho_
    w = [1.0]
    for n in range(1, nmax + 1):
        w.append(w[-1] * ratios[n])
    return w


# direct quadrature of the Legendre transform of exp(-b(1 - cos g)) for a few degrees
def legendre_p(n, x):
    p0, p1 = 1.0, x
    if n == 0:
        return 1.0
    for k in range(2, n + 1):
        p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
    return p1


def w_quad(radius, n):
    b = math.log(2) / (1 - math.cos(radius / R_E))
    # the kernel is concentrated near x = 1: integrate on [1 - 40/b, 1]
    f0 = lambda x: math.exp(-b * (1 - x))
    lo = max(-1.0, 1 - 40 / b)
    return integrate(lambda x: f0(x) * legendre_p(n, x), lo, 1.0, 40) / integrate(f0, lo, 1.0, 40)


wb = weights_backward(300e3, 60)
check("W_30 at 300 km: backward recursion = Legendre transform of the Gaussian", wb[30], w_quad(300e3, 30), rel=1e-10)
for radius, nmax in ((300e3, 60), (300e3, 96), (500e3, 120), (1000e3, 96), (100e3, 180)):
    fw = weights_forward(radius, nmax, floor=False)
    fwf = weights_forward(radius, nmax)
    bw_ = weights_backward(radius, nmax)
    err = max(abs(a - c) for a, c in zip(fw, bw_))
    errf = max(abs(a - c) for a, c in zip(fwf, bw_))
    info(f"Gaussian weights {radius / 1e3:.0f} km to n={nmax}: forward error (no floor / with crate floor)", f"{err:.1e} / {errf:.1e}")

print("=" * 72, "\nF. Spherical harmonics: Gauss-Legendre analysis is exact; Slepian polar cap")


def pbar_all(L, x):
    """Fully normalised (4 pi) associated Legendre functions, no Condon-Shortley phase."""
    s = math.sqrt(max(0.0, 1 - x * x))
    P = [[0.0] * (L + 1) for _ in range(L + 1)]
    P[0][0] = 1.0
    for m in range(1, L + 1):
        P[m][m] = (math.sqrt(3.0) if m == 1 else math.sqrt((2 * m + 1) / (2 * m))) * s * P[m - 1][m - 1]
    for m in range(0, L):
        P[m + 1][m] = math.sqrt(2 * m + 3) * x * P[m][m]
    for m in range(0, L + 1):
        for n in range(m + 2, L + 1):
            a_ = math.sqrt((2 * n - 1) * (2 * n + 1) / ((n - m) * (n + m)))
            b_ = math.sqrt((2 * n + 1) * (n + m - 1) * (n - m - 1) / ((n - m) * (n + m) * (2 * n - 3)))
            P[n][m] = a_ * x * P[n - 1][m] - b_ * P[n - 2][m]
    return P


Lsh_ = 16
random.seed(7)
coef = {(n, m): (random.gauss(0, 1), random.gauss(0, 1) if m else 0.0) for n in range(Lsh_ + 1) for m in range(n + 1)}
xg, wg = gauss_legendre(Lsh_ + 1)
nlon = 2 * Lsh_ + 2
grid = []
for x in xg:
    P = pbar_all(Lsh_, x)
    row = []
    for j in range(nlon):
        lam = 2 * math.pi * j / nlon
        row.append(sum(P[n][m] * (cc * math.cos(m * lam) + ss * math.sin(m * lam)) for (n, m), (cc, ss) in coef.items()))
    grid.append(row)
maxerr = 0.0
Ps = [pbar_all(Lsh_, x) for x in xg]
for (n, m), (cc, ss) in coef.items():
    ca = sa = 0.0
    for i, x in enumerate(xg):
        fc = sum(grid[i][j] * math.cos(m * 2 * math.pi * j / nlon) for j in range(nlon)) * 2 * math.pi / nlon
        fs = sum(grid[i][j] * math.sin(m * 2 * math.pi * j / nlon) for j in range(nlon)) * 2 * math.pi / nlon
        ca += wg[i] * Ps[i][n][m] * fc / (4 * math.pi)
        sa += wg[i] * Ps[i][n][m] * fs / (4 * math.pi)
    maxerr = max(maxerr, abs(ca - cc), abs(sa - ss))
check(f"Gauss-Legendre analysis recovers a degree-{Lsh_} field exactly ((L+1) x (2L+2) grid)", maxerr, 0.0, abs_=1e-12)


def jacobi_eig(A):
    n = len(A)
    A = [row[:] for row in A]
    for _ in range(100):
        off = sum(A[i][j] ** 2 for i in range(n) for j in range(n) if i != j)
        if off < 1e-28:
            break
        for p in range(n):
            for q_ in range(p + 1, n):
                if abs(A[p][q_]) < 1e-300:
                    continue
                th = (A[q_][q_] - A[p][p]) / (2 * A[p][q_])
                t = (1 if th >= 0 else -1) / (abs(th) + math.sqrt(th * th + 1))
                c_ = 1 / math.sqrt(t * t + 1)
                s_ = t * c_
                for k in range(n):
                    akp, akq = A[k][p], A[k][q_]
                    A[k][p], A[k][q_] = c_ * akp - s_ * akq, s_ * akp + c_ * akq
                for k in range(n):
                    apk, aqk = A[p][k], A[q_][k]
                    A[p][k], A[q_][k] = c_ * apk - s_ * aqk, s_ * apk + c_ * aqk
    return sorted((A[i][i] for i in range(n)), reverse=True)


Lc, cap = 18, math.radians(40)
xq, wq_ = gauss_legendre(80)
x0 = math.cos(cap)
nodes = [(x0 + (1 - x0) * (x + 1) / 2, w * (1 - x0) / 2) for x, w in zip(xq, wq_)]
Pn = [pbar_all(Lc, x) for x, _ in nodes]
eigs, trace = [], 0.0
for m in range(Lc + 1):
    ls = list(range(m, Lc + 1))
    fac = 0.5 if m == 0 else 0.25  # (1/4pi) * (2pi or pi)
    D = [[fac * sum(w * Pn[i][l1][m] * Pn[i][l2][m] for i, (_, w) in enumerate(nodes)) for l2 in ls] for l1 in ls]
    ev = jacobi_eig(D)
    mult = 1 if m == 0 else 2
    eigs += ev * mult
    trace += mult * sum(D[i][i] for i in range(len(ls)))
shannon = (Lc + 1) ** 2 * (1 - math.cos(cap)) / 2
check("Slepian polar cap: sum of eigenvalues = Shannon number (L+1)^2 (1-cos T)/2", sum(eigs), shannon, rel=1e-12)
eigs.sort(reverse=True)
info(f"Slepian cap L={Lc}, Theta=40 deg: N = {shannon:.3f}; lambda_1..3, lambda_N, count > 0.5",
     f"{eigs[0]:.6f} {eigs[1]:.6f} {eigs[2]:.6f}, {eigs[round(shannon) - 1]:.3f}, {sum(1 for e_ in eigs if e_ > 0.5)}")

print("=" * 72, "\nG. Multipole truncation bound and Barnes-Hut")
random.seed(1)
src = []
while len(src) < 60:
    p = [random.uniform(-1, 1) for _ in range(3)]
    if norm(p) <= 1:
        src.append((p, random.uniform(0.5, 1.5)))
Qtot = sum(qq for _, qq in src)
target = [1.2, -2.0, 1.9]
rt = norm(target)
exact = sum(qq / norm(sub(target, p)) for p, qq in src)
ok = True
for pmax in range(0, 13):
    approx = 0.0
    for p, qq in src:
        rp = norm(p)
        cg = dot(p, target) / (rp * rt) if rp > 0 else 1.0
        approx += qq * sum(rp**n / rt ** (n + 1) * legendre_p(n, cg) for n in range(pmax + 1))
    bound = Qtot / (rt - 1.0) * (1.0 / rt) ** (pmax + 1)
    ok &= abs(approx - exact) <= bound
check("multipole error <= Q/(r-a) (a/r)^(p+1) for p = 0..12 (Greengard & Rokhlin 1987)", 1.0 if ok else 0.0, 1.0, abs_=0)

random.seed(3)
Nb = 1000
bodies = [[random.random() for _ in range(3)] for _ in range(Nb)]
mass = 1.0 / Nb


class Node:
    __slots__ = ("c", "h", "m", "com", "kids", "idx")

    def __init__(self, c, h):
        self.c, self.h, self.m, self.com, self.kids, self.idx = c, h, 0.0, [0.0, 0.0, 0.0], None, None


def insert(node, i):
    p = bodies[i]
    if node.kids is None and node.idx is None and node.m == 0.0:
        node.idx = i
    else:
        if node.kids is None:
            node.kids = [None] * 8
            j, node.idx = node.idx, None
            put(node, j)
        put(node, i)
    node.com = [(node.com[k] * node.m + p[k] * mass) / (node.m + mass) for k in range(3)]
    node.m += mass


def put(node, i):
    p = bodies[i]
    o = (p[0] > node.c[0]) + 2 * (p[1] > node.c[1]) + 4 * (p[2] > node.c[2])
    if node.kids[o] is None:
        hh = node.h / 2
        node.kids[o] = Node([node.c[k] + (hh if (o >> k) & 1 else -hh) for k in range(3)], hh)
    insert(node.kids[o], i)


root = Node([0.5, 0.5, 0.5], 0.5)
for i in range(Nb):
    insert(root, i)


def force(i, node, theta, count):
    if node is None or node.m == 0.0 or node.idx == i:
        return [0.0, 0.0, 0.0]
    d = sub(node.com, bodies[i])
    r = norm(d)
    if node.kids is None or (2 * node.h) / r < theta:
        count[0] += 1
        return [node.m * x / r**3 for x in d]
    tot = [0.0, 0.0, 0.0]
    for kid in node.kids:
        f_ = force(i, kid, theta, count)
        tot = [a + b for a, b in zip(tot, f_)]
    return tot


direct = []
for i in range(Nb):
    tot = [0.0, 0.0, 0.0]
    for j in range(Nb):
        if i != j:
            d = sub(bodies[j], bodies[i])
            r = norm(d)
            tot = [tot[k] + mass * d[k] / r**3 for k in range(3)]
    direct.append(tot)
for theta in (0.3, 0.5, 0.7, 1.0):
    errs, count = [], [0]
    for i in range(Nb):
        f_ = force(i, root, theta, count)
        errs.append(norm(sub(f_, direct[i])) / norm(direct[i]))
    errs.sort()
    info(f"Barnes-Hut monopole, N={Nb} uniform cube, theta={theta}", f"median rel. force error {errs[Nb // 2]:.1e}, 99th pct {errs[int(0.99 * Nb)]:.1e}, interactions/body {count[0] / Nb:.0f} (direct {Nb - 1})")

print("=" * 72, "\nH. Dimensionless numbers: one table, exact identities")
sigma, nu, rho_f, Bf, Lf, Uf, Om = 1e7, 7e-7, 930.0, 0.1, 1.0, 1.0, 10.0  # liquid sodium
eta = 1 / (MU0 * sigma)
vA = Bf / math.sqrt(MU0 * rho_f)
Rm = Uf * Lf / eta
S = vA * Lf / eta
Pm = nu / eta
Ha = Bf * Lf * math.sqrt(sigma / (rho_f * nu))
Re = Uf * Lf / nu
Elsasser = sigma * Bf**2 / (rho_f * Om)
Em = eta / (Om * Lf**2)
check("Hartmann = Lundquist / sqrt(magnetic Prandtl)", Ha, S / math.sqrt(Pm), rel=1e-13)
check("Rm = Re * Pm", Rm, Re * Pm, rel=1e-13)
check("Elsasser = Lundquist^2 * magnetic Ekman", Elsasser, S**2 * Em, rel=1e-13)
check("Lundquist = Rm at U = v_A", S, vA * Lf * MU0 * sigma, rel=1e-13)
ratio = 4 * math.pi * EPS0 * G * MP**2 / QE**2
check("(k_Jeans/k_Debye)^2 = G m_p^2 / (e^2/4 pi eps0) for hydrogen", 4 * math.pi * G * MP**2 / (QE**2 / EPS0), ratio, rel=1e-13)
info("G m_p^2 / (e^2 / 4 pi eps0)", f"{ratio:.4e}")

print("=" * 72)
print("FAILED:" if FAILS else "ALL CHECKS PASSED", *FAILS, sep="\n  ")
