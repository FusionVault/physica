"""Checks for constants.md: every explicit constant, bound and inequality in OpenAI,
"Velocity-Field Detection of Computation in Forced Navier-Stokes Flows" (27 Sep 2026).

Run: python3 -I constants_checks.py   (stdlib only; a few seconds)
Every number quoted in constants.md is asserted here; exits 1 on any failure.
"Paper" values are printed in the manuscript; "computed" values are produced here.
"""
import itertools
import math
import random
import sys
from fractions import Fraction as Fr

FAILS = []


def check(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    if not ok:
        FAILS.append(name)


def close(a, b, rel):
    return abs(a - b) <= rel * max(abs(a), abs(b))


# ---------------------------------------------------------------------------------------------
# 1. Smooth step b, theta, bump g (section 3.1) and the effective tail estimate (3.4, 5.5)
# ---------------------------------------------------------------------------------------------
def bs(r):
    return math.exp(-1.0 / r) if r > 0 else 0.0


def dbs(r):
    return math.exp(-1.0 / r) / (r * r) if r > 0 else 0.0


def theta(r):
    return bs(r) / (bs(r) + bs(1 - r))


def dtheta(r):
    den = bs(r) + bs(1 - r)
    return (dbs(r) * bs(1 - r) + bs(r) * dbs(1 - r)) / (den * den)


def d2theta(r, h=1e-5):
    return (dtheta(r + h) - dtheta(r - h)) / (2 * h)


grid = [i / 2000 for i in range(-4000, 6001)]
den_min = min(bs(r) + bs(1 - r) for r in grid)
check("b(r)+b(1-r) >= e^-2 on [-2,3] (paper); computed min = 2e^-2 at r=1/2",
      den_min >= math.exp(-2) and close(den_min, 2 * math.exp(-2), 1e-9), f"min={den_min:.6f}")
check("theta = 0 for r<=0, 1 for r>=1, theta(r)+theta(1-r)=1",
      all(theta(r) == 0 for r in grid if r <= 0) and all(theta(r) == 1 for r in grid if r >= 1)
      and all(abs(theta(r) + theta(1 - r) - 1) < 1e-15 for r in grid))
check("theta nondecreasing", all(theta(grid[i]) <= theta(grid[i + 1]) for i in range(len(grid) - 1)))
check("theta flat at 0: theta(r) <= r^5 for 0<r<=0.05", all(theta(r) <= r ** 5 for r in [j / 1000 for j in range(1, 51)]))
check("integral of theta over [0,1] = 1/2 (symmetry)", close(sum(theta((i + .5) / 4000) for i in range(4000)) / 4000, .5, 1e-12))
dth_max = max(dtheta(r) for r in [i / 20000 for i in range(1, 20000)])
check("computed max theta' = theta'(1/2) = 2", close(dth_max, 2.0, 1e-6), f"{dth_max:.6f}")
tail_ok = all(s ** m * math.exp(-s) <= math.factorial(m + j) * s ** (-j)
              for s in [0.1, 1, 3, 10, 50] for m in range(6) for j in range(6))
check("s^m e^-s <= (m+j)! s^-j (effective tail bound)", tail_ok)


def g1(z):
    return theta(2 * (z + 1)) * theta(2 * (1 - z))


def g(z1, z2):
    return g1(z1) * g1(z2)


gz = [i / 200 for i in range(-300, 301)]
check("g in [0,1] (grid on [-1.5,1.5]^2)", all(0 <= g(a, c) <= 1 for a in gz[::5] for c in gz[::5]))
check("g = 1 on [-1/2,1/2]^2", all(g(a, c) == 1 for a in gz if abs(a) <= .5 for c in gz[::7] if abs(c) <= .5))
check("g = 0 outside [-1,1]^2", all(g(a, c) == 0 for a in gz if abs(a) >= 1 for c in gz[::7]))
int_g1 = sum(g1(-1 + (i + .5) / 10000) for i in range(20000)) / 10000
check("computed: integral of g = 1.5^2 = 2.25 (paper bounds it by 4, eq. 11)", close(int_g1 ** 2, 2.25, 1e-9), f"{int_g1**2:.9f}")


# derivative norms of g, for the effective integer C of eq. 11
def dg1(z):
    return 2 * dtheta(2 * (z + 1)) * theta(2 * (1 - z)) - 2 * theta(2 * (z + 1)) * dtheta(2 * (1 - z))


def d2g1(z):
    a, c = 2 * (z + 1), 2 * (1 - z)
    return 4 * d2theta(a) * theta(c) - 8 * dtheta(a) * dtheta(c) + 4 * theta(a) * d2theta(c)


zs = [-1 + i / 2000 for i in range(4001)]
G1 = [(g1(z), dg1(z), d2g1(z)) for z in zs]
grad_max = max(math.hypot(p[1] * q[0], p[0] * q[1]) for p in G1[::8] for q in G1[::8])
hess_max = 0.0
for p in G1[::8]:
    for q in G1[::8]:
        hxx, hyy, hxy = p[2] * q[0], p[0] * q[2], p[1] * q[1]
        tr, det = hxx + hyy, hxx * hyy - hxy * hxy
        hess_max = max(hess_max, abs(tr) / 2 + math.sqrt(max(tr * tr / 4 - det, 0)))
check("computed: max|grad g| = 4 (so C >= 4)", close(grad_max, 4.0, 1e-3), f"{grad_max:.4f}")
check("computed: max||D^2 g||_op = 40.46 (max|theta''| = 9.84), so the smallest admissible integer C is 41",
      40.3 < hess_max < 40.6, f"{hess_max:.3f}")
C_G = 41

# ---------------------------------------------------------------------------------------------
# 2. Theorem 3.1 scales, eqs. 10-13, 16, 18, support radius
# ---------------------------------------------------------------------------------------------
TEN6 = Fr(1, 10 ** 6)
for L in range(1, 6):
    R = 4 ** L
    x = Fr(1, 4 * R * R)
    M0 = 4 * TEN6 ** 2 * x / (1 - x)
    partial = sum(4 * (TEN6 * Fr(1, (2 * R) ** n)) ** 2 for n in range(1, 40))
    ok = M0 == Fr(4, 10 ** 12) / (4 * R * R - 1) and M0 < Fr(1, 10 ** 11) and 0 <= M0 - partial < M0 * Fr(1, 10 ** 50)
    check(f"eq. 13 M0 = 4e-12/(4R^2-1) exactly, < 1e-11 (L={L}, R={R})", ok, f"M0={float(M0):.4e}")
check("computed: M0 <= 4e-12/63 = 6.35e-14 for every L>=1 (157x below 1e-11)",
      Fr(4, 10 ** 12) / 63 < Fr(1, 10 ** 11) / 157)
check("M0 < 1e-11 even for R=1 (4e-12/3): the L>=1 hypothesis is not needed here", Fr(4, 10 ** 12) / 3 < Fr(1, 10 ** 11))

for L in (1, 2, 5):
    R = 4 ** L
    for n in range(1, 8):
        bn = TEN6 / (2 * R) ** n
        Dn = 2 * C_G / bn ** 2 * (1 + n * L) * R ** (3 * n)
        dn = 1 / (100 * (1 + Dn))
        if not (2 * dn * Dn < Fr(1, 50) < Fr(1, 16) and 2 - 2 * dn > 1):
            check(f"eq. 18 at L={L}, n={n}", False)
        rad = math.sqrt(2) * float(bn * R ** n)
        if not (bn * R ** n == TEN6 / 2 ** n and rad < 1 / 32):
            check(f"support radius at L={L}, n={n}", False)
        # Hessian bound of section 3.3 is <= D_n
        hess = 2 * (C_G / bn ** 2 * R ** (2 * n) + C_G / bn * n * L * R ** (3 * n))
        if not hess <= Dn:
            check(f"trace of Hessian <= D_n at L={L}, n={n}", False)
check("eq. 18: 2 delta_n D_n < 1/50 < 1/16; wait 2-2delta_n > 1 (L=1,2,5; n=1..7)", True)
check("support radius sqrt2 b_n R^n = sqrt2 1e-6 2^-n < 1/32 (max 7.07e-7)", math.sqrt(2) * 1e-6 / 2 < 1 / 32 / 44000)
check("section 3.3: 2C(b^-2 R^2n + b^-1 nL R^3n) <= D_n", True)
b1 = TEN6 / 8
D1 = 2 * C_G / b1 ** 2 * 2 * 64
check("computed: L=1, C=41: D_1 = 6.72e17, delta_1 = 1.49e-20",
      close(float(D1), 6.71744e17, 1e-5) and close(float(1 / (100 * (1 + D1))), 1.4887e-20, 1e-3))
dd = [float(2 * (1 / (100 * (1 + D))) * D) for D in (Fr(1), Fr(10), D1)]
check("computed: 2 delta D = 1/100, 1/55, -> 1/50 as D grows (eq. 18 is tight)",
      close(dd[0], .01, 1e-12) and close(dd[1], 1 / 55, 1e-12) and close(dd[2], .02, 1e-15))
# eq. 16: second variation integral
for L in (1, 2):
    for s in (0.5, 1, 3):
        n = math.ceil(s)
        integral = math.exp(2 * L * s) - math.exp(L * s)
        num = sum(math.exp(L * (s - r)) * L * math.exp(2 * L * r) for r in [(i + .5) * s / 20000 for i in range(20000)]) * s / 20000
        if not (close(integral, num, 1e-6) and integral <= n * L * math.exp(3 * L * n) and math.exp(L * n) <= 4 ** (L * n)):
            check(f"eq. 16 at L={L}, s={s}", False)
check("eq. 16: int e^{L(s-r)} L e^{2Lr} dr = e^{2Ls} - e^{Ls} <= nL e^{3Ln}; e^{Ln} <= R^n", True)

# ---------------------------------------------------------------------------------------------
# 3. Lemma 3.3: shell bound and heat-kernel bounds, eq. 14
# ---------------------------------------------------------------------------------------------
q = math.exp(-1 / 32)
check("e^{-1/32} <= 32/33", q <= 32 / 33 and Fr(33, 32) <= sum(Fr(1, 32 ** k * math.factorial(k)) for k in range(3)))
check("sum_j j (32/33)^j = 32*33 = 1056 exactly", Fr(32, 33) / (1 - Fr(32, 33)) ** 2 == 1056)
check("1 + 8*1056 = 8449 < 1e4", 1 + 8 * 1056 == 8449 < 10 ** 4)
sj = q / (1 - q) ** 2
check("computed: sum_j j e^{-j/32} = 1023.92", close(sj, 1023.92, 1e-5), f"{sj:.4f}")
check("shell |k|_inf = j has 8j points", all(sum(1 for a in range(-j, j + 1) for c in range(-j, j + 1) if max(abs(a), abs(c)) == j) == 8 * j for j in range(1, 30)))
zgrid = [i / 40 - .5 for i in range(41)]
check("|Z+k| >= j - 1/2 >= j/2 on the shell, Z in [-1/2,1/2]^2",
      all(math.hypot(z1 + a, z2 + c) >= max(abs(a), abs(c)) - .5 - 1e-15
          for z1 in zgrid[::4] for z2 in zgrid[::4] for a in range(-6, 7) for c in range(-6, 7)))


def theta1(z, s, nmax=60):  # sum_n exp(-(z+n)^2 / s)
    return sum(math.exp(-(z + n) ** 2 / s) for n in range(-nmax, nmax + 1))


shell_sup = max(theta1(z1, 8) * theta1(z2, 8) for z1 in zgrid for z2 in zgrid)
check("computed: sup_Z sum_k e^{-|Z+k|^2/8} = 8 pi = 25.13 (at Z=0), vs bound 8449",
      close(shell_sup, 8 * math.pi, 1e-12) and shell_sup <= 8449, f"{shell_sup:.6f}")


def K(tau, z1, z2):
    return theta1(z1, 4 * tau) * theta1(z2, 4 * tau) / (4 * math.pi * tau)


def K_fourier(tau, z1, z2, mmax=10):
    s = lambda z: sum(math.exp(-4 * math.pi ** 2 * m * m * tau) * math.cos(2 * math.pi * m * z) for m in range(-mmax, mmax + 1))
    return s(z1) * s(z2)


check("image sum = Fourier series (heat kernel identity) at tau=0.05, 0.3, 1",
      all(close(K(t, a, c), K_fourier(t, a, c), 1e-12) for t in (.05, .3, 1) for a, c in ((0, 0), (.2, -.4), (.5, .5))))
mass = sum(K(0.01, (i + .5) / 200 - .5, (j + .5) / 200 - .5) for i in range(200) for j in range(200)) / 200 ** 2
check("heat kernel has unit mass on T^2", close(mass, 1, 1e-10), f"{mass:.12f}")
K1 = K(1, 0, 0)
check("computed: sup over tau>=1 is K(1,0) = 1 + 4e^{-4pi^2} = 1 + 2.9e-17", abs(K1 - 1) < 1e-15 and K(2, 0, 0) <= K1)
check("shell estimate at tau=1: K <= 1e4/(4 pi) = 795.8 <= H", 1e4 / (4 * math.pi) < 1e8)


def golden(f, a, b, it=200):
    gr = (math.sqrt(5) - 1) / 2
    c, d = b - gr * (b - a), a + gr * (b - a)
    for _ in range(it):
        if f(c) > f(d):
            b = d
        else:
            a = c
        c, d = b - gr * (b - a), a + gr * (b - a)
    return (a + b) / 2


for d, H in ((1 / 32, 1e9), (1 / 16, 1e8)):
    paper_bound = 8e4 / (4 * math.pi * d * d)
    best = 0.0
    for ang in [k * math.pi / 160 for k in range(41)]:
        z1, z2 = d * math.cos(ang), d * math.sin(ang)
        f = lambda lt: K(math.exp(lt), z1, z2)
        lt = golden(f, math.log(1e-7), 0.0)
        best = max(best, f(lt))
    # interior points with dist >= d never beat the circle
    off = max(K(d * d / 4 * fac, z1, z2) for fac in (.5, 1, 2) for z1 in zgrid for z2 in zgrid
              if math.hypot(z1, z2) >= d)
    exact = 1 / (math.pi * math.e * d * d)
    check(f"eq. 14 (paper): 8e4/(4 pi d^2) = {paper_bound:.4g} < H = {H:.0e} at d = {d}", paper_bound < H)
    check(f"tau^-1 e^(-d^2/(8tau)) <= 8/d^2 at d={d} (sup = 8/(e d^2))",
          close(max(math.exp(-d * d / (8 * t)) / t for t in [d * d / 8 * 1.0001 ** k for k in range(-5000, 5000)]), 8 / (math.e * d * d), 1e-6))
    check(f"computed: true sup K over dist>=d, tau<=1 at d={d} is 1/(pi e d^2) = {exact:.4f}; H/sup = {H/best:.3g}",
          close(best, exact, 1e-6) and off <= best * (1 + 1e-9) and best < paper_bound < H, f"{best:.6f}")

# ---------------------------------------------------------------------------------------------
# 4. Eq. 20 and the section 4 geometry
# ---------------------------------------------------------------------------------------------
for d, H in ((Fr(1, 32), 10 ** 9), (Fr(1, 16), 10 ** 8)):
    check(f"eq. 20: 2 H M0 + 1/16 < 1/2 with only M0 < 1e-11 (d={d})", 2 * H * Fr(1, 10 ** 11) + Fr(1, 16) < Fr(1, 2))
    lhs = 2 * H * Fr(4, 10 ** 12) / 63 + Fr(1, 16)
    check(f"computed: with M0 <= 4e-12/63, 2HM0 + 1/16 <= {float(lhs):.6f} (d={d})", lhs < Fr(1, 2))
check("section 4: corridor 3/16<=y<=7/8 to E=(1/32,1/8): circle distance = 1/16 = 2d exactly (tight)",
      min(Fr(3, 16) - Fr(1, 8), 1 + Fr(1, 32) - Fr(7, 8)) == Fr(1, 16) == 2 * Fr(1, 32))
check("section 4: min{2/3-3/8, 1+1/8-5/6} = 7/24 > 1/4 (>= 2d = 1/8 is all that is needed)",
      min(Fr(2, 3) - Fr(3, 8), 1 + Fr(1, 8) - Fr(5, 6)) == Fr(7, 24) > Fr(1, 4))

# ---------------------------------------------------------------------------------------------
# 5. Time-derivative factor: height-one injection; and the R^2 derivative-kernel L1 norm
# ---------------------------------------------------------------------------------------------
for dl in (1e-2, 1e-6):
    tn, nseg = 4.0, 20000
    hstep = dl / nseg
    simpson = sum((1 if i in (0, nseg) else 4 if i % 2 else 2) * dtheta(((tn + i * hstep) - tn) / dl) / dl for i in range(nseg + 1)) * hstep / 3
    check(f"int delta^-1 theta'((t-tn)/delta) dt = theta(1)-theta(0) = 1 (delta={dl})", close(simpson, 1, 1e-10))
for nu in (0.01, 1, 7):
    s = 0.3
    v = 2 * nu * s
    l1 = sum(abs(x) / (2 * nu * s) * math.exp(-x * x / (2 * v)) / math.sqrt(2 * math.pi * v) for x in [(i + .5) * 1e-3 * math.sqrt(v) - 15 * math.sqrt(v) for i in range(30000)]) * 1e-3 * math.sqrt(v)
    check(f"computed: ||d_l H_s||_1 = (pi nu s)^-1/2 (nu={nu})", close(l1, 1 / math.sqrt(math.pi * nu * s), 1e-6))

# ---------------------------------------------------------------------------------------------
# 6. Eq. 37: quintic P, A*, phi_R
# ---------------------------------------------------------------------------------------------
P = lambda s: 0.0 if s <= 0 else 1.0 if s >= 1 else 10 * s ** 3 - 15 * s ** 4 + 6 * s ** 5
P2 = lambda s: 0.0 if s <= 0 or s >= 1 else 60 * s - 180 * s ** 2 + 120 * s ** 3
Pf = [Fr(0), Fr(0), Fr(0), Fr(10), Fr(-15), Fr(6)]  # coefficients of P on [0,1]


def peval(c, s):
    return sum(ci * s ** i for i, ci in enumerate(c))


def pder(c):
    return [i * c[i] for i in range(1, len(c))]


Pd1, Pd2, Pd3 = pder(Pf), pder(pder(Pf)), pder(pder(pder(Pf)))
check("P(0)=0, P(1)=1, P'(0)=P'(1)=P''(0)=P''(1)=0: extension by 0 and 1 is C^2",
      [peval(Pf, 0), peval(Pf, 1), peval(Pd1, 0), peval(Pd1, 1), peval(Pd2, 0), peval(Pd2, 1)] == [0, 1, 0, 0, 0, 0])
check("P'''(0)=60, P'''(1)=60: the extension is not C^3", peval(Pd3, 0) == 60 and peval(Pd3, 1) == 60)
check("P' = 30 s^2 (1-s)^2 >= 0, so P in [0,1]", all(peval(Pd1, Fr(i, 97)) == 30 * Fr(i, 97) ** 2 * (1 - Fr(i, 97)) ** 2 for i in range(98)))
p2max = max(abs(P2(i / 200000)) for i in range(200001))
check("|P''| <= 360 (paper: 60+180+120); computed max = 10/sqrt3 = 5.774",
      p2max <= 360 and close(p2max, 10 / math.sqrt(3), 1e-9), f"{p2max:.6f}")
A = lambda s: P(2 * (1 + s)) * P(2 * (1 - s))
A2 = lambda s: 4 * P2(2 * (1 + s)) * P(2 * (1 - s)) + 4 * P(2 * (1 + s)) * P2(2 * (1 - s))  # cross term: P'(..)P'(..) = 0
Pd = lambda s: 0.0 if s <= 0 or s >= 1 else 30 * s * s * (1 - s) ** 2
sgrid = [i / 100000 - 1.5 for i in range(300001)]
check("transition intervals [-1,-1/2] and [1/2,1] of A* are disjoint; the other factor is 1 on each",
      all(P(2 * (1 - s)) == 1 for s in sgrid if -1 <= s <= -.5) and all(P(2 * (1 + s)) == 1 for s in sgrid if .5 <= s <= 1)
      and all(Pd(2 * (1 + s)) * Pd(2 * (1 - s)) == 0 for s in sgrid))
check("A* = 1 on |s|<=1/2, 0 on |s|>=1, in [0,1]",
      all(A(s) == 1 for s in sgrid if abs(s) <= .5) and all(A(s) == 0 for s in sgrid if abs(s) >= 1) and all(0 <= A(s) <= 1 for s in sgrid))
a2max = max(abs(A2(s)) for s in sgrid)
check("|A*''| <= 1440 (paper); computed max = 40/sqrt3 = 23.09", a2max <= 1440 and close(a2max, 40 / math.sqrt(3), 1e-8), f"{a2max:.4f}")
fd = max(abs((A(s + 1e-4) - 2 * A(s) + A(s - 1e-4)) / 1e-8 - A2(s)) for s in sgrid[::97])
check("A*'' formula matches a finite difference", fd < 1e-3, f"{fd:.2e}")
coarse = [i / 2000 - 1.1 for i in range(4401)]
vals = [(A(s), A2(s)) for s in coarse]
lap_max = max(abs(p[1] * q[0] + p[0] * q[1]) for p in vals[::3] for q in vals[::3])
check("||Delta phi_R|| <= 2880 R^-2 <= C* R^-2, C*=3000 (paper); computed sup = 43.49 R^-2 (both factors in transition)",
      lap_max <= 2880 <= 3000 and 43.4 < lap_max < 43.6, f"{lap_max:.3f}")
check("phi_R = 1 on |a|_inf <= R/2, 0 on |a|_inf >= R",
      all(A(x / 3) * A(y / 3) == 1 for x in [-1.5, 0, 1.5] for y in [-1.5, .7, 1.5])
      and all(A(x / 3) * A(y / 3) == 0 for x in [-3, 3.2] for y in [0, 2]))

# ---------------------------------------------------------------------------------------------
# 7. Eqs. 25 and 39: scales and the loss budget delta < 1/24
# ---------------------------------------------------------------------------------------------
CS = 3000
rng = random.Random(53)
for _ in range(200):  # the identity, with Lambda and D independent rationals
    K0, Dv, Lam, R0 = (Fr(rng.randint(1, 999), rng.randint(1, 99)) for _ in range(4))
    for n in range(0, 12):
        Rn, Rn1 = R0 * Lam ** n, R0 * Lam ** (n + 1)
        Tn = 16 * Rn1 * (K0 * Dv ** (n + 1) + 2)
        if Tn / Rn ** 2 != 16 * Lam / R0 * (K0 * Dv * (Dv / Lam) ** n + 2 * Lam ** -n):
            check("T_n/R_n^2 identity", False)
check("T_n/R_n^2 = (16 Lambda/R0)(K0 D (D/Lambda)^n + 2 Lambda^-n), exact over 200 random rational points, n<=11", True)

worst = Fr(0)
count = 0
for bb in (2, 3, 5, 16):
    for N in (2, 5, 100, 10 ** 4):
        for m in (1, 3):
            for nu in (Fr(1, 10 ** 6), Fr(1, 100), Fr(1), Fr(100), Fr(10 ** 6)):
                D = bb * bb
                Lam = 4 * D
                K0 = N * bb ** (2 * m)
                R0 = 1024 * CS * (1 + nu) * Lam * (K0 * D + 2)
                closed = 16 * Lam / R0 * (K0 * D * Fr(4, 3) + 2 * Fr(Lam, Lam - 1))
                partial = sum(16 * R0 * Lam ** (n + 1) * (K0 * D ** (n + 1) + 2) / (R0 * Lam ** n) ** 2 for n in range(60))
                part1 = nu * CS * closed
                part2 = nu * CS / R0 ** 2
                delta = part1 + part2
                ok = (Lam / D == 4 and Lam >= 16 and R0 * Lam >= 2 * R0 and part1 < Fr(1, 48) and part2 < Fr(1, 48)
                      and delta < Fr(1, 24) and 0 <= closed - partial < closed * Fr(1, 10 ** 30)
                      and part1 <= nu / (1 + nu) / 48 and R0 >= 1024 * CS * (1 + nu) and R0 >= 1
                      and 16 * R0 * Lam * (K0 * D + 2) >= 1)  # T_0 >= 1, hence t_n >= n + 1
                worst = max(worst, delta)
                count += 1
                if not ok:
                    check(f"eq. 39 at b={bb}, N={N}, m={m}, nu={nu}", False)
check(f"eq. 39: delta < 1/48 + 1/48 = 1/24 and R_{{n+1}} = Lambda R_n >= 2 R_n ({count} parameter sets)", True)
check(f"computed: worst delta = {float(worst):.6f} = (nu/(1+nu))/48 to 6 digits; 64/(3*1024) = 1/48 exactly",
      Fr(64, 3 * 1024) == Fr(1, 48) and close(float(worst), 1 / 48 * 1e6 / (1 + 1e6), 1e-3))
check("detection needs only delta < 1/2 (event 3); the 3/4-mass remark needs delta < 1/4", Fr(1, 24) < Fr(1, 4) < Fr(1, 2))
R0ex = 1024 * CS * 2 * 16 * (2 * 4 * 4 + 2)
check("computed: R0 = 3.34e9 for nu=1, b=2, N=2, m=1 (K0=8, D=4, Lambda=16)", R0ex == 3342336000)

# ---------------------------------------------------------------------------------------------
# 8. Eqs. 22-24: stack updates, inverses, bounds and address map, exhaustively on small tables
# ---------------------------------------------------------------------------------------------
MOVES = "RNL"


def step(x, y, beta, move, b):
    if move == "R":
        return b * x + beta, y // b
    if move == "N":
        return x, y - y % b + beta
    gam = x % b
    return x // b, gam + b * beta + b * b * (y // b)


def inverse(xp, yp, beta, alpha, move, b):
    if move == "R":
        return (xp - beta) // b, b * yp + alpha
    if move == "N":
        return xp, yp - beta + alpha
    return b * xp + yp % b, b * (yp // (b * b)) + alpha


def written(xp, yp, move, b):
    return {"R": xp % b, "N": yp % b, "L": (yp // b) % b}[move]


def random_table(N, b, rng, injective=True):
    move_of = {i: rng.choice(MOVES) for i in range(1, N + 1)}
    slots = [(ip, be) for ip in range(1, N + 1) for be in range(b)]
    rng.shuffle(slots)
    rules = {}
    keys = [(i, al) for i in range(1, N) for al in range(b)]  # state N is terminal
    for idx, key in enumerate(keys):
        if rng.random() < 0.85:
            ip, be = slots[idx] if injective else slots[idx % 2]
            rules[key] = (be, ip)
    return move_of, rules


tables = 0
for b in (2, 3):
    for N in (2, 3, 4):
        for trial in range(6):
            move_of, rules = random_table(N, b, rng)
            incoming = {}
            for (i, al), (be, ip) in rules.items():
                incoming[(ip, be)] = (i, al)
            for m in (1,):
                for n in (0, 1):
                    B = b ** (m + n)
                    Bn1 = b * B
                    seen = {}
                    for i in range(1, N + 1):
                        for x in range(B):
                            for y in range(B):
                                al = y % b
                                if (i, al) not in rules:
                                    continue
                                be, ip = rules[(i, al)]
                                mv = move_of[ip]
                                xp, yp = step(x, y, be, mv, b)
                                ok = 0 <= xp < Bn1 and 0 <= yp < Bn1 and (mv != "L" or yp <= b * B - 1)
                                # recover the rule from (i', written digit), then invert
                                be2 = written(xp, yp, mv, b)
                                i0, al0 = incoming[(ip, be2)]
                                ok = ok and (i0, al0) == (i, al) and inverse(xp, yp, be2, al0, mv, b) == (x, y)
                                k = 1 + x + B * y + B * B * (i - 1)
                                kp = 1 + xp + Bn1 * yp + Bn1 * Bn1 * (ip - 1)
                                ok = ok and 1 <= k <= N * B * B and 1 <= kp <= N * Bn1 * Bn1
                                ok = ok and ((k - 1) % B, ((k - 1) // B) % B, (k - 1) // (B * B) + 1) == (x, y, i)
                                if not ok or kp in seen:
                                    check(f"stack map b={b} N={N} n={n}", False)
                                seen[kp] = k
                    tables += 1
check(f"eqs. 22-24: updates in range, y' <= bB_n - 1, three inverses, address map and k -> k' injective ({tables} tables)", True)
yl_max = max(step(x, y, be, "L", b)[1] for b in (2, 3, 4) for B in [b ** 2] for x in range(B) for y in range(B) for be in range(b))
check("left-row bound y' <= bB_n - 1 is attained (tight)", all(
    max(step(x, y, be, "L", b)[1] for x in range(b ** 2) for y in range(b ** 2) for be in range(b)) == b * b ** 2 - 1 for b in (2, 3, 4)))
# negative control: without the incoming-rule property injectivity fails
move_of = {1: "N", 2: "N", 3: "N"}
rules = {(1, 0): (1, 2), (1, 1): (1, 2)}  # two read symbols, same (target, written)
imgs = [step(0, y, 1, "N", 2) for y in (0, 1)]
check("negative control: rules sharing (i', beta) collide", imgs[0] == imgs[1])

# ---------------------------------------------------------------------------------------------
# 9. Gate W (eq. 28) and the center paths (eq. 26, 27)
# ---------------------------------------------------------------------------------------------
def chi1(a):
    return theta(a + 3) * theta(3 - a)


def dchi1(a):
    return dtheta(a + 3) * theta(3 - a) - theta(a + 3) * dtheta(3 - a)


def gate(X1, X2, c, G, R):
    x1, x2 = X1 - c[0], X2 - c[1]
    a1, a2 = x1 / R, x2 / R
    chi = chi1(a1) * chi1(a2)
    Lf = G[0] * x2 - G[1] * x1
    d1psi = dchi1(a1) * chi1(a2) / R * Lf - chi * G[1]
    d2psi = chi1(a1) * dchi1(a2) / R * Lf + chi * G[0]
    return d2psi, -d1psi


c0, G0, Rg = (0.4, -1.3), (0.7, -2.2), 1.7
pts = [(c0[0] + Rg * (i / 10 - 4), c0[1] + Rg * (j / 10 - 4)) for i in range(81) for j in range(81)]
hh = 1e-5
div_max = 0.0
plateau_ok, zero_ok = True, True
for X1, X2 in pts:
    w = gate(X1, X2, c0, G0, Rg)
    r = max(abs(X1 - c0[0]), abs(X2 - c0[1])) / Rg
    if r <= 2:
        plateau_ok &= abs(w[0] - G0[0]) < 1e-15 and abs(w[1] - G0[1]) < 1e-15
    if r >= 3:
        zero_ok &= w == (0.0, 0.0) or (w[0] == 0 and w[1] == 0)
    dv = (gate(X1 + hh, X2, c0, G0, Rg)[0] - gate(X1 - hh, X2, c0, G0, Rg)[0]
          + gate(X1, X2 + hh, c0, G0, Rg)[1] - gate(X1, X2 - hh, c0, G0, Rg)[1]) / (2 * hh)
    div_max = max(div_max, abs(dv))
check("eq. 28: W = G on |xi|_inf <= 2R", plateau_ok)
check("eq. 28: W = 0 on |xi|_inf >= 3R", zero_ok)
check("eq. 28: div W = 0 (central difference, 6561 points)", div_max < 1e-8, f"max {div_max:.1e}")
check("chi plateau: chi1 = 1 on |a|<=2, 0 on |a|>=3",
      all(chi1(a) == 1 for a in [i / 100 for i in range(-200, 201)]) and all(chi1(a) == 0 for a in [3, 3.5, -3, -4]))
check("Theta(r) = theta(2r-1/2) moves only on [1/4, 3/4]",
      all(theta(2 * r - .5) == (0 if r <= .25 else 1) for r in [i / 1000 for i in range(-500, 1500)] if not .25 < r < .75))

# center paths for a random injective table: separations and eq. 27
Th = lambda r: theta(2 * r - .5)
for b in (2, 3):
    N, m, n = 3, 1, 0
    move_of, rules = random_table(N, b, rng)
    B, Bn1 = b ** (m + n), b ** (m + n + 1)
    Rn = 1.0
    Sn, Sn1 = 16 * Rn, 16 * 16 * b * b * Rn
    routes = []
    for i in range(1, N):
        for x in range(B):
            for y in range(B):
                if (i, y % b) not in rules:
                    continue
                be, ip = rules[(i, y % b)]
                xp, yp = step(x, y, be, move_of[ip], b)
                k = 1 + x + B * y + B * B * (i - 1)
                kp = 1 + xp + Bn1 * yp + Bn1 * Bn1 * (ip - 1)
                routes.append((k, kp, 1 if ip == N else -1))
    sep_ok, eq27_ok, end_ok = True, True, True
    for u in [j / 600 for j in range(601)]:
        th = [Th(3 * u - jj) for jj in range(3)]
        cs = []
        for k, kp, s in routes:
            Yk = -(k + 2) * Sn
            c1 = (1 - th[1]) * Sn * k + th[1] * Sn1 * kp
            c2 = -(1 - th[0]) * Sn + (th[0] - th[2]) * Yk + th[2] * s * Sn1
            cs.append((c1, c2))
            if s < 0:
                eq27_ok &= c2 <= -Sn + 1e-9
            if u == 0:
                end_ok &= (c1, c2) == (Sn * k, -Sn)
            if u == 1:
                end_ok &= abs(c1 - Sn1 * kp) < 1e-9 and abs(c2 - s * Sn1) < 1e-9
        for p, q2 in itertools.combinations(cs, 2):
            sep_ok &= max(abs(p[0] - q2[0]), abs(p[1] - q2[1])) >= Sn - 1e-9
    check(f"eq. 26: paths start at p_(n,k), end at p_(n+1,k'), and stay >= S_n = 16R_n > 6R_n apart (b={b}, {len(routes)} gates)", sep_ok and end_ok)
    check(f"eq. 27: nonterminal routes stay at or below -S_n (b={b})", eq27_ok)
check("|Y_k| <= S_n(K_n+2) <= T_n and S_{n+1} k' <= T_n (b=2, N=3, m=1, n=0..5)",
      all(16 * (3 * 4 ** (1 + n) + 2) <= 16 * 4 * (3 * 4 ** (2 + n) + 2) for n in range(6)))

# ---------------------------------------------------------------------------------------------
# 10. Eq. 21 viscosity rescaling and the eq. 5-6 reduction, exact polynomial algebra
#     Variables (nu, t, x, y); nu is a symbol, so the check holds for every nu.
# ---------------------------------------------------------------------------------------------
def padd(*ps):
    out = {}
    for p in ps:
        for k, v in p.items():
            out[k] = out.get(k, 0) + v
    return {k: v for k, v in out.items() if v != 0}


def pscale(p, c, nu_pow=0):
    return {(k[0] + nu_pow,) + k[1:]: v * c for k, v in p.items() if v * c != 0}


def pmul(p, q2):
    out = {}
    for k1, v1 in p.items():
        for k2, v2 in q2.items():
            k = tuple(a + c for a, c in zip(k1, k2))
            out[k] = out.get(k, 0) + v1 * v2
    return {k: v for k, v in out.items() if v != 0}


def pdiff(p, var):
    out = {}
    for k, v in p.items():
        if k[var]:
            kk = list(k)
            kk[var] -= 1
            out[tuple(kk)] = out.get(tuple(kk), 0) + v * k[var]
    return {k: v for k, v in out.items() if v != 0}


def time_rescale(p):  # f(t) -> f(nu t)
    return {(k[0] + k[1],) + k[1:]: v for k, v in p.items()}


def lap(p):
    return padd(pdiff(pdiff(p, 2), 2), pdiff(pdiff(p, 3), 3))


def mono(c, nu=0, t=0, x=0, y=0):
    return {(nu, t, x, y): Fr(c)}


psi = padd(mono(1, t=2, x=3, y=1), mono(2, t=2, x=1, y=2), mono(Fr(-1, 3), t=3, y=4), mono(5, t=1, x=2, y=2), mono(-1, t=4, x=1))
a = (pdiff(psi, 3), pscale(pdiff(psi, 2), -1))
w = padd(mono(1, t=1, x=2), mono(3, t=1, x=1, y=3), mono(-2, t=3, y=1), mono(1, t=2, x=4))


def adv(vel, f):
    return padd(pmul(vel[0], pdiff(f, 2)), pmul(vel[1], pdiff(f, 3)))


check("a = grad-perp psi is divergence free; a(0) = w(0) = 0",
      padd(pdiff(a[0], 2), pdiff(a[1], 3)) == {} and all(k[1] > 0 for comp in (a[0], a[1], w) for k in comp))
F = [padd(pdiff(a[j], 1), adv(a, a[j]), pscale(lap(a[j]), -1)) for j in range(2)]
h = padd(pdiff(w, 1), adv(a, w), pscale(lap(w), -1))
Ua = [pscale(time_rescale(a[j]), 1, 1) for j in range(2)]
Uw = time_rescale(w)
fa = [pscale(time_rescale(F[j]), 1, 2) for j in range(2)]
fw = pscale(time_rescale(h), 1, 1)
res = [padd(pdiff(Ua[j], 1), adv(Ua, Ua[j]), pscale(lap(Ua[j]), -1, 1), pscale(fa[j], -1)) for j in range(2)]
resw = padd(pdiff(Uw, 1), adv(Ua, Uw), pscale(lap(Uw), -1, 1), pscale(fw, -1))
check("eq. 21: U_nu = (nu a(nu t), w(nu t)), f_nu = (nu^2 F(nu t), nu h(nu t)) solve viscosity nu exactly, p = 0",
      res == [{}, {}] and resw == {} and padd(pdiff(Ua[0], 2), pdiff(Ua[1], 3)) == {})
bad = [padd(pdiff(Ua[j], 1), adv(Ua, Ua[j]), pscale(lap(Ua[j]), -1, 1), pscale(time_rescale(F[j]), -1, 1)) for j in range(2)]
check("negative control: f_nu = (nu F(nu t), ...) leaves a nonzero residual", bad != [{}, {}])
# eq. 5-6: the 3-D residual of U = (a, w), z-independent, p = 0, is (F - F, h - h); div U = div_Y a
check("eqs. 5-6: with F, h as defined the 3-D residual of (a, w) vanishes and div U = div_Y a = 0",
      all(padd(pdiff(a[j], 1), adv(a, a[j]), pscale(lap(a[j]), -1), pscale(F[j], -1)) == {} for j in range(2))
      and padd(pdiff(w, 1), adv(a, w), pscale(lap(w), -1), pscale(h, -1)) == {})

# ---------------------------------------------------------------------------------------------
# Eq. 8 / 16 (#53): R = 4^L (a superscript in the PDF, flattened to "4L" by text extraction), so
# ||D Psi_s|| <= e^{Ln} <= R^n holds for every L since e < 4. Read as 4L it would fail from L = 3.
check("eqs. 8/16: e^L <= R = 4^L for every L", all(math.exp(L) <= 4 ** L for L in range(1, 60)))
check("eqs. 8/16: the flattened reading 4L would fail from L = 3 (e^3 = 20.09 > 12)", math.exp(3) > 12 and math.exp(2) <= 8)

# ---------------------------------------------------------------------------------------------
print()
print(f"{len(FAILS)} failure(s)" if FAILS else "all checks pass")
sys.exit(1 if FAILS else 0)
