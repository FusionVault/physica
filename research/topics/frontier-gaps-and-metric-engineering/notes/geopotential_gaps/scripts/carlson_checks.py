"""Gap 3: Carlson's R_C and R_J Cauchy principal values, and the incomplete Legendre integrals
F(phi|m), E(phi|m), Pi(n; phi|m) built on them, checked against Carlson (1995), DLMF ch. 19,
A&S tables 17.5/17.6/17.9, scipy (Carlson-based elliprf/rd/rj/rc and Cephes ellipkinc/ellipeinc)
and direct quadrature. Run: python3 -I carlson_checks.py
The R_* ports follow geopotential::special line by line (Carlson 1995, r = eps)."""
import math
import numpy as np
from scipy import special, integrate, optimize

EPS = np.finfo(float).eps


def rf(x, y, z):
    a0 = (x + y + z) / 3.0
    q = (3.0 * EPS) ** (-1 / 6) * max(abs(a0 - x), abs(a0 - y), abs(a0 - z))
    v, a, f = [x, y, z], a0, 1.0
    while f * q >= abs(a):
        sx, sy, sz = (math.sqrt(w) for w in v)
        lam = sx * sy + sx * sz + sy * sz
        v, a, f = [(w + lam) / 4 for w in v], (a + lam) / 4, f / 4
    dx, dy = (a0 - x) * f / a, (a0 - y) * f / a
    dz = -dx - dy
    e2, e3 = dx * dy - dz * dz, dx * dy * dz
    return (1 - e2 / 10 + e3 / 14 + e2 * e2 / 24 - 3 * e2 * e3 / 44) / math.sqrt(a)


def series(e2, e3, e4, e5):
    return 1 - 3 * e2 / 14 + e3 / 6 + 9 * e2 * e2 / 88 - 3 * e4 / 22 - 9 * e2 * e3 / 52 + 3 * e5 / 26


def rd(x, y, z):
    a0 = (x + y + 3 * z) / 5
    q = (EPS / 4) ** (-1 / 6) * max(abs(a0 - x), abs(a0 - y), abs(a0 - z))
    v, a, f, s = [x, y, z], a0, 1.0, 0.0
    while f * q >= abs(a):
        sx, sy, sz = (math.sqrt(w) for w in v)
        lam = sx * sy + sx * sz + sy * sz
        s += f / (sz * (v[2] + lam))
        v, a, f = [(w + lam) / 4 for w in v], (a + lam) / 4, f / 4
    dx, dy = (a0 - x) * f / a, (a0 - y) * f / a
    dz = -(dx + dy) / 3
    e2, e3 = dx * dy - 6 * dz * dz, (3 * dx * dy - 8 * dz * dz) * dz
    e4, e5 = 3 * (dx * dy - dz * dz) * dz * dz, dx * dy * dz ** 3
    return f / (a * math.sqrt(a)) * series(e2, e3, e4, e5) + 3 * s


def rc_pos(x, y):
    """R_C for x >= 0, y > 0 (the crate's carlson_rc)."""
    a0 = (x + 2 * y) / 3
    q = (3 * EPS) ** (-1 / 8) * abs(a0 - x)
    v, a, f = [x, y], a0, 1.0
    while f * q >= abs(a):
        lam = 2 * math.sqrt(v[0]) * math.sqrt(v[1]) + v[1]
        v, a, f = [(w + lam) / 4 for w in v], (a + lam) / 4, f / 4
    s = (y - a0) * f / a
    return (1 + s * s * (3 / 10 + s * (1 / 7 + s * (3 / 8 + s * (9 / 22 + s * (159 / 208 + s * 9 / 8)))))) / math.sqrt(a)


def rc(x, y):
    """R_C with the Cauchy principal value for y < 0: Carlson (1995) eq. (2.13a),
    R_C(x, -y) = sqrt(x/(x+y)) R_C(x+y, y), y > 0 (DLMF 19.2.20 is the same identity)."""
    if y > 0:
        return rc_pos(x, y)
    return math.sqrt(x / (x - y)) * rc_pos(x - y, -y)


def rj_pos(x, y, z, p):
    a0 = (x + y + z + 2 * p) / 5
    delta = (p - x) * (p - y) * (p - z)
    q = (EPS / 4) ** (-1 / 6) * max(abs(a0 - w) for w in (x, y, z, p))
    v, a, f, s = [x, y, z, p], a0, 1.0, 0.0
    while f * q >= abs(a):
        sx, sy, sz, sp = (math.sqrt(w) for w in v)
        lam = sx * sy + sx * sz + sy * sz
        d = (sp + sx) * (sp + sy) * (sp + sz)
        s += f / d * rc_pos(1.0, 1.0 + f ** 3 * delta / (d * d))
        v, a, f = [(w + lam) / 4 for w in v], (a + lam) / 4, f / 4
    dx, dy, dz = ((a0 - w) * f / a for w in (x, y, z))
    dp = -(dx + dy + dz) / 2
    e2 = dx * dy + dx * dz + dy * dz - 3 * dp * dp
    e3 = dx * dy * dz + 2 * e2 * dp + 4 * dp ** 3
    e4 = (2 * dx * dy * dz + e2 * dp + 3 * dp ** 3) * dp
    e5 = dx * dy * dz * dp * dp
    return f / (a * math.sqrt(a)) * series(e2, e3, e4, e5) + 6 * s


def rj_carlson(x, y, z, p):
    """R_J; for p < 0 the principal value by Carlson (1995) eq. (2.23a) with y the middle value:
    (y+q) R_J(x,y,z,-q) = (p-y) R_J(x,y,z,p) - 3 R_F(x,y,z) + 3 sqrt(xyz/(xz+pq)) R_C(xz+pq, pq),
    p - y = (z-y)(y-x)/(y+q)."""
    if p > 0:
        return rj_pos(x, y, z, p)
    x, y, z = sorted((x, y, z))  # middle value in y makes (z-y)(y-x) >= 0
    q = -p
    pp = y + (z - y) * (y - x) / (y + q)
    return ((pp - y) * rj_pos(x, y, z, pp) - 3 * rf(x, y, z)
            + 3 * math.sqrt(x * y * z / (x * z + pp * q)) * rc_pos(x * z + pp * q, pp * q)) / (y + q)


def rj_dlmf(x, y, z, p):
    """R_J principal value by DLMF 19.20.14-15 with z the largest value:
    (q+z) R_J(x,y,z,-q) = (p-z) R_J(x,y,z,p) - 3 R_F + 3 sqrt(xyz/(xy+pq)) R_C(xy+pq, pq),
    p = (z(x+y+q) - xy)/(z+q)."""
    if p > 0:
        return rj_pos(x, y, z, p)
    x, y, z = sorted((x, y, z))
    q = -p
    pp = (z * (x + y + q) - x * y) / (z + q)
    return ((pp - z) * rj_pos(x, y, z, pp) - 3 * rf(x, y, z)
            + 3 * math.sqrt(x * y * z / (x * y + pp * q)) * rc_pos(x * y + pp * q, pp * q)) / (q + z)


def rj_quad(x, y, z, p):
    """(3/2) PV int_0^inf dt / ((t+p) sqrt((t+x)(t+y)(t+z))) by QUADPACK's Cauchy weight."""
    g = lambda t: 1.0 / math.sqrt((t + x) * (t + y) * (t + z))
    if p > 0:
        return 1.5 * integrate.quad(lambda t: g(t) / (t + p), 0, np.inf, epsabs=0, epsrel=1e-13, limit=500)[0]
    q = -p
    a = integrate.quad(g, 0, 2 * q, weight='cauchy', wvar=q, epsabs=0, epsrel=1e-13, limit=500)[0]
    # tail t = 2q + u/(1-u)
    b = integrate.quad(lambda t: g(t) / (t - q), 2 * q, np.inf, epsabs=0, epsrel=1e-13, limit=500)[0]
    return 1.5 * (a + b)


def rel(a, b):
    return abs(a / b - 1)


print("=== Carlson (1995) section 3 check values, principal values ===")
cases = [("R_C(1/4,-2) = ln2/3", rc(0.25, -2.0), 0.23104906018665, math.log(2) / 3),
         ("R_J(2,3,4,-0.5)", rj_carlson(2, 3, 4, -0.5), 0.24723819703052, None),
         ("R_J(2,3,4,-5)", rj_carlson(2, 3, 4, -5.0), -0.12711230042964, None)]
for name, got, table, exact in cases:
    print(f"{name:22s} ours {got:.15f}  table {table:.14f}  rel {rel(got, table):.1e}"
          + (f"  exact {exact:.16f} rel {rel(got, exact):.1e}" if exact else ""))
for (x, y, z, p) in [(2, 3, 4, -0.5), (2, 3, 4, -5.0), (0.5, 1, 7, -3.0), (0.1, 1, 2, -0.7)]:
    a, b, c, d = rj_carlson(x, y, z, p), rj_dlmf(x, y, z, p), special.elliprj(x, y, z, p), rj_quad(x, y, z, p)
    print(f"R_J({x},{y},{z},{p}): Carlson2.23a {a:.16f} DLMF19.20.14 {b:.16f} scipy {c:.16f} quadPV {d:.16f}"
          f"  max rel spread {max(rel(a, b), rel(a, c), rel(a, d)):.1e}")
for (x, y) in [(0.25, -2.0), (1.0, -1.0), (3.0, -0.1), (0.0, -2.0)]:
    print(f"R_C({x},{y}) ours {rc(x, y):.16f} scipy {special.elliprc(x, y):.16f}")
# Where the principal value changes sign (Carlson: near p = -1.2552) and the digits lost near it.
root = optimize.brentq(lambda p: rj_carlson(2, 3, 4, p), -2, -1, xtol=1e-15)
print("  quadPV at the root and at Carlson -1.2552:", rj_quad(2, 3, 4, root), rj_quad(2, 3, 4, -1.2552), " ours at -1.2552:", rj_carlson(2, 3, 4, -1.2552))
print(f"R_J(2,3,4,p) = 0 at p = {root:.10f} (Carlson: approximately -1.2552)")
for p in [root * (1 + 1e-6), -1.25, -1.3]:
    a, c = rj_carlson(2, 3, 4, p), special.elliprj(2, 3, 4, p)
    print(f"  p = {p:.9f}: ours {a:.6e} scipy {c:.6e} rel diff {rel(a, c):.1e}")
# DLMF 19.20.8 closed forms for the principal value
y, q = 2.0, 3.0
print("DLMF 19.20.8 R_J(0,y,y,-q) = -3pi/(2 sqrt(y)(y+q)):", rj_carlson(0, y, y, -q), -3 * math.pi / (2 * math.sqrt(y) * (y + q)))
print("DLMF 19.20.8 R_J(x,y,y,p) = 3(R_C(x,y)-R_C(x,p))/(p-y), x=1,y=2,p=-3:",
      rj_carlson(1, 2, 2, -3.0), 3 * (rc(1, 2) - rc(1, -3)) / (-3 - 2))

print("\n=== Incomplete Legendre integrals through Carlson ===")


def F(phi, m):
    """F(phi|m) = int_0^phi (1 - m sin^2)^(-1/2), DLMF 19.25.5 in the form sin(phi) R_F(cos^2, 1-m sin^2, 1),
    with DLMF 19.2.10 for |phi| > pi/2. m > 1 is real only for m sin^2 phi <= 1, which also rules out
    |phi| > pi/2 (returns nan)."""
    j = round(phi / math.pi)
    r = phi - j * math.pi  # |r| <= pi/2
    s, c = math.sin(r), math.cos(r)
    d2 = 1 - m * s * s
    if d2 < 0 or (j != 0 and m >= 1):
        return math.nan
    val = s * rf(c * c, d2, 1.0)
    return val + 2 * j * rf(0.0, 1 - m, 1.0) if j else val


def E(phi, m):
    """E(phi|m): DLMF 19.25.9 for m <= 0, 19.25.10 for 0 < m <= 1 (no cancellation near m = 1,
    DLMF 19.36(i), Reinsch & Raab 2000), 19.25.11 for 1 < m <= csc^2 phi; DLMF 19.2.10 beyond pi/2."""
    j = round(phi / math.pi)
    r = phi - j * math.pi
    sgn = -1.0 if r < 0 else 1.0
    s, c = abs(math.sin(r)), math.cos(r)
    # 1 - m s^2 without cancellation near m = 1, phi = pi/2 (Boost's (c-1) + (1-k^2))
    x, y = c * c, (c * c + (1 - m) * s * s if 0 < m <= 1 else 1 - m * s * s)
    if y < 0 or (j != 0 and m > 1):
        return math.nan
    if s == 0:
        val = 0.0
    elif m <= 0:
        val = s * rf(x, y, 1) - m / 3 * s ** 3 * rd(x, y, 1)
    elif m == 1:
        val = s  # DLMF 19.6.9
    elif m < 1:
        val = (1 - m) * s * rf(x, y, 1) + m * (1 - m) / 3 * s ** 3 * rd(x, 1, y) + m * s * c / math.sqrt(y)
    else:
        val = -(1 - m) / 3 * s ** 3 * rd(y, 1, x) + s * math.sqrt(y) / c
    val *= sgn
    if j:
        val += 2 * j * (rf(0.0, 1 - m, 1.0) - m / 3 * rd(0.0, 1 - m, 1.0))
    return val


def E_naive(phi, m):
    s, c = math.sin(phi), math.cos(phi)
    return s * rf(c * c, 1 - m * s * s, 1) - m / 3 * s ** 3 * rd(c * c, 1 - m * s * s, 1)


def Pi(n, phi, m):
    """Pi(n; phi|m) = int_0^phi dθ / ((1 - n sin^2θ) sqrt(1 - m sin^2θ)) (DLMF 19.2.7 / A&S 17.2.14 sign),
    = s R_F(x,y,1) + (n/3) s^3 R_J(x, y, 1, 1 - n s^2)  (DLMF 19.25.14 by homogeneity);
    1 - n s^2 < 0 is the Cauchy principal value (R_J principal value); DLMF 19.2.10-style
    periodicity with the complete Pi beyond pi/2."""
    j = round(phi / math.pi)
    r = phi - j * math.pi
    s, c = math.sin(r), math.cos(r)
    x, y, p = c * c, 1 - m * s * s, 1 - n * s * s
    val = s * rf(x, y, 1) + n / 3 * s ** 3 * rj_carlson(x, y, 1.0, p) if s != 0 else 0.0
    if j:
        val += 2 * j * (rf(0.0, 1 - m, 1.0) + n / 3 * rj_carlson(0.0, 1 - m, 1.0, 1 - n))
    return val


def Pi_dlmf_19_25_16(n, phi, m):
    """Pi for alpha^2 = n > c = csc^2 phi by DLMF 19.25.16 (no negative R_J argument needed)."""
    c = 1 / math.sin(phi) ** 2
    w2 = m / n
    return (-w2 / 3 * rj_pos(c - 1, c - m, c, c - w2)
            + math.sqrt((c - 1) * (c - m) / ((n - 1) * (1 - w2))) * rc(c * (n - 1) * (1 - w2), (n - c) * (c - w2)))


def quad_pv(f_regular, theta0, phi):
    """PV of int_0^phi f(θ)/(θ - θ0) dθ with QUADPACK's Cauchy weight."""
    return integrate.quad(f_regular, 0, phi, weight='cauchy', wvar=theta0, epsabs=0, epsrel=1e-13, limit=500)[0]


def Pi_quad(n, phi, m):
    g = lambda t: 1 / ((1 - n * math.sin(t) ** 2) * math.sqrt(1 - m * math.sin(t) ** 2))
    if n * math.sin(min(abs(phi), math.pi / 2)) ** 2 <= 1 or abs(phi) > math.pi / 2:
        if n <= 1:
            return integrate.quad(g, 0, phi, epsabs=0, epsrel=1e-13, limit=500)[0]
    t0 = math.asin(1 / math.sqrt(n))  # pole of 1/(1 - n sin^2) in (0, pi/2)
    # g(t) = h(t)/(t - t0) with h regular at t0
    h = lambda t: (t - t0) * g(t) if t != t0 else -1 / (n * math.sin(2 * t0) * math.sqrt(1 - m / n))
    return quad_pv(h, t0, phi)


d = math.radians
print("A&S Table 17.5 F(phi\\alpha), m = sin^2 alpha:")
for phi, al, tab in [(30, 45, 0.53562273), (30, 60, 0.54222911), (30, 90, 0.54930614), (20, 40, 0.35198869)]:
    m = math.sin(d(al)) ** 2
    v = F(d(phi), m)
    print(f"  F({phi}°\\{al}°) ours {v:.10f}  A&S {tab:.8f}  diff {v - tab:.1e}  scipy {special.ellipkinc(d(phi), m):.15f}")
print("A&S Table 17.6 E(phi\\alpha):")
for phi, al, tab in [(30, 45, 0.51204932), (30, 60, 0.50609207), (90, 45, 1.35064388), (70, 60, 1.02663689),
                     (75, 85, 0.96992212), (80, 88, 0.98568915), (85, 88, 0.99748392)]:
    m = math.sin(d(al)) ** 2
    v, vn = E(d(phi), m), E_naive(d(phi), m)
    print(f"  E({phi}°\\{al}°) ours {v:.10f}  A&S {tab:.8f}  diff {v - tab:.1e}  scipy {special.ellipeinc(d(phi), m):.15f}"
          f"  19.25.10-vs-19.25.9 {v - vn:.1e}")
print("A&S Table 17.9 Pi(n; phi\\alpha):")
for n, phi, al, tab in [(0.5, 60, 45, 1.38218), (0.3, 75, 60, 1.93879), (0.6, 45, 30, 0.91527), (0.5, 90, 45, 2.70129),
                        (0.1, 90, 0, 1.65576), (0.0, 75, 90, 2.02759)]:
    m = math.sin(d(al)) ** 2
    v = Pi(n, d(phi), m)
    print(f"  Pi({n};{phi}°\\{al}°) ours {v:.10f}  A&S {tab:.5f}  diff {v - tab:.1e}  quad {Pi_quad(n, d(phi), m):.12f}")

print("Exact identities (DLMF 19.6, 19.2.10):")
phi, m = 1.2, 0.7
print("  F(phi,1) = gd^-1(phi) = atanh(sin phi):", F(0.7, 1.0), math.atanh(math.sin(0.7)))
print("  E(phi,1) = sin phi:", E(0.7, 1.0), math.sin(0.7))
print("  F(phi,0) = phi; E(phi,0) = phi:", F(1.3, 0.0), E(1.3, 0.0))
print("  Pi(phi,1,0) = tan phi (DLMF 19.6.11):", Pi(1.0, 1.1, 0.0), math.tan(1.1))
print("  Pi(phi,n,0) = R_C(c-1, c-n) (DLMF 19.6.12), n=3, phi=1.1 (PV):", Pi(3.0, 1.1, 0.0),
      rc(1 / math.sin(1.1) ** 2 - 1, 1 / math.sin(1.1) ** 2 - 3))
Kc, Ec = special.ellipk(m), special.ellipe(m)
for jj, ph in [(1, 0.4), (2, -0.3), (-1, 1.0)]:
    print(f"  19.2.10 F({jj}pi+{ph}) ours {F(jj * math.pi + ph, m):.15f} vs 2jK+F {2 * jj * Kc + F(ph, m):.15f} scipy {special.ellipkinc(jj * math.pi + ph, m):.15f}")
    print(f"  19.2.10 E({jj}pi+{ph}) ours {E(jj * math.pi + ph, m):.15f} vs 2jE+E {2 * jj * Ec + E(ph, m):.15f} scipy {special.ellipeinc(jj * math.pi + ph, m):.15f}")
print("  F, E odd in phi:", F(-0.8, 0.6) + F(0.8, 0.6), E(-0.8, 0.6) + E(0.8, 0.6))
# m > 1: DLMF 19.7.4 reciprocal modulus, F(phi|m) = F(beta|1/m)/sqrt(m), sin beta = sqrt(m) sin phi
m, phi = 2.5, 0.5
beta = math.asin(math.sqrt(m) * math.sin(phi))
print(f"  m>1 F: ours {F(phi, m):.15f}  19.7.4 {F(beta, 1 / m) / math.sqrt(m):.15f}  scipy {special.ellipkinc(phi, m):.15f}")
print(f"  m>1 E: ours {E(phi, m):.15f}  19.7.4 {(E(beta, 1 / m) - (1 - 1 / m) * F(beta, 1 / m)) * math.sqrt(m):.15f}  scipy {special.ellipeinc(phi, m):.15f}")
print("  m>1 beyond asin(1/sqrt m) (complex):", F(1.0, m), E(1.0, m))
# m < 0: DLMF 19.7.5 imaginary modulus is not needed; compare with quadrature
for m in [-3.0, -50.0]:
    fq = integrate.quad(lambda t: 1 / math.sqrt(1 - m * math.sin(t) ** 2), 0, 1.4, epsabs=0, epsrel=1e-13)[0]
    eq = integrate.quad(lambda t: math.sqrt(1 - m * math.sin(t) ** 2), 0, 1.4, epsabs=0, epsrel=1e-13)[0]
    print(f"  m={m}: F rel {rel(F(1.4, m), fq):.1e}  E rel {rel(E(1.4, m), eq):.1e}")
print("Pi with n > 1 (principal value) vs QUADPACK Cauchy weight and DLMF 19.25.16 / 19.6.5:")
for n, phi, m in [(2.0, 1.2, 0.5), (4.0, 0.9, 0.3), (1.5, 1.4, 0.9), (10.0, 0.5, 0.2)]:
    a, b = Pi(n, phi, m), Pi_quad(n, phi, m)
    cdl = Pi_dlmf_19_25_16(n, phi, m)
    print(f"  Pi({n};{phi}|{m}) ours {a:.15f} quadPV {b:.15f} DLMF19.25.16 {cdl:.15f} rel {rel(a, b):.1e} {rel(a, cdl):.1e}")
for n, m in [(2.0, 0.5), (5.0, 0.9)]:
    print(f"  complete Pi({n}|{m}) ours {Pi(n, math.pi / 2, m):.15f}  DLMF19.6.5 K - Pi(m/n) {special.ellipk(m) - Pi(m / n, math.pi / 2, m):.15f}")
print("  periodicity with n>1: Pi(2; pi+0.3|0.5)", Pi(2.0, math.pi + 0.3, 0.5), " 2Pi(2|0.5)+Pi(2;0.3|0.5)",
      2 * Pi(2.0, math.pi / 2, 0.5) + Pi(2.0, 0.3, 0.5))
