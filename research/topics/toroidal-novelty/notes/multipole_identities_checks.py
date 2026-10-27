#!/usr/bin/env python3
"""Numerical checks for multipole_identities.md (toroidal-family moments of closed currents).

Run with numpy, isolated from the working directory:  python3 -P multipole_identities_checks.py
Every identity is asserted; the script prints one line per check and exits 1 if any fails.

Test currents (all closed, ∇·J = 0):
  * a 24-turn poloidal winding on a torus (polygon), shifted and tilted off the origin;
  * a trefoil knot (polygon);
  * 120 random closed polygons (generic closed currents, for the rank count);
  * J = (∇×)^n K for a random sum of vector Gaussians K, by FFT on a 96³ grid (n = 1..6).
Polygon moments use Gauss-Legendre nodes on each straight segment, so moments of polynomials
are exact to rounding; the identities then hold to ~1e-13.
SI-like, no 1/c: T = (1/10)∫[(r·J)r − 2r²J] dV (Dubovik-Tugushev, Ederer-Spaldin).
"""
import math
import sys

import numpy as np

FAILED = []


def check(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    if not ok:
        FAILED.append(name)


def rel(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float(np.max(np.abs(a - b)) / max(np.max(np.abs(b)), 1e-300))


# ---------------------------------------------------------------- current representations
GL_X, GL_W = np.polynomial.legendre.leggauss(8)  # exact for polynomials of degree <= 15


def polygon(pts, current=1.0):
    """Closed polygon (vertices pts, closing edge implied) -> nodes X and weights W = I dl dw."""
    p = np.asarray(pts, float)
    q = np.roll(p, -1, axis=0)
    t = 0.5 * (GL_X + 1.0)
    X = (p[:, None, :] * (1 - t)[None, :, None] + q[:, None, :] * t[None, :, None]).reshape(-1, 3)
    W = (current * (q - p)[:, None, :] * (0.5 * GL_W)[None, :, None]).reshape(-1, 3)
    return X, W


def join(*cs):
    return np.vstack([c[0] for c in cs]), np.vstack([c[1] for c in cs])


def torus_winding(turns=24, m=40, big_r=1.0, a=0.3, current=1.0, centre=(0, 0, 0), axis=(0, 0, 1)):
    """`turns` poloidal circles (m-gons) on a torus about `axis` through `centre`."""
    axis = np.asarray(axis, float) / np.linalg.norm(axis)
    e1 = np.cross(axis, [1.0, 0, 0] if abs(axis[0]) < 0.9 else [0, 1.0, 0])
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(axis, e1)
    loops = []
    for j in range(turns):
        ph = 2 * math.pi * (j + 0.5) / turns
        rho_hat = math.cos(ph) * e1 + math.sin(ph) * e2
        th = 2 * math.pi * np.arange(m) / m
        pts = (np.asarray(centre) + np.outer(big_r + a * np.cos(th), rho_hat)
               + np.outer(a * np.sin(th), axis))
        loops.append(polygon(pts, current))
    return join(*loops)


def trefoil(n=400, scale=0.5, shift=(0.3, -0.2, 0.1)):
    t = 2 * math.pi * np.arange(n) / n
    pts = scale * np.stack([np.sin(t) + 2 * np.sin(2 * t), np.cos(t) - 2 * np.cos(2 * t), -np.sin(3 * t)], 1)
    return polygon(pts + np.asarray(shift))


def mom(c, f):
    """∫ f(r) J dV for a scalar weight f given as a function of X (n,3) -> (n,)."""
    X, W = c
    return (f(X)[:, None] * W).sum(0)


def r2(X):
    return (X * X).sum(1)


def rdotJ(c):
    return (c[0] * c[1]).sum(1)


def toroidal_dubovik(c):
    X, W = c
    return 0.1 * (rdotJ(c)[:, None] * X - 2 * r2(X)[:, None] * W).sum(0)


# ---------------------------------------------------------------- FFT curl fields
N_GRID, L_BOX = 96, 6.0
xs = (np.arange(N_GRID) - N_GRID // 2) * (2 * L_BOX / N_GRID)
GX, GY, GZ = np.meshgrid(xs, xs, xs, indexing="ij")
GR = np.stack([GX, GY, GZ], -1)
DV = (2 * L_BOX / N_GRID) ** 3
kk = 2 * math.pi * np.fft.fftfreq(N_GRID, d=2 * L_BOX / N_GRID)
KX, KY, KZ = np.meshgrid(kk, kk, kk, indexing="ij")


def curl(F):
    Fh = [np.fft.fftn(F[..., i]) for i in range(3)]
    k = (KX, KY, KZ)
    out = []
    for i in range(3):
        j, l = (i + 1) % 3, (i + 2) % 3
        out.append(np.real(np.fft.ifftn(1j * (k[j] * Fh[l] - k[l] * Fh[j]))))
    return np.stack(out, -1)


rng = np.random.default_rng(20261027)
K_FIELD = np.zeros(GR.shape)
for _ in range(4):
    cen, amp = rng.uniform(-1, 1, 3), rng.normal(size=3)
    K_FIELD += amp * np.exp(-((GR - cen) ** 2).sum(-1) / (2 * 0.6 ** 2))[..., None]
GRID_R2 = (GR * GR).sum(-1)


def gmom(F, f):
    return (f[..., None] * F).sum((0, 1, 2)) * DV


def gcross(F):
    return np.cross(GR, F)


# ================================================================ Q1/Q2: one-term reductions
solenoid = torus_winding(centre=(0.4, -0.3, 0.25), axis=(0.3, 0.2, 1.0))
knot = trefoil()
CURRENTS = {"tilted torus winding": solenoid, "trefoil knot": knot}
FIELDS = {}
M1 = K_FIELD  # J = curl M1, M1 localized
J1 = curl(M1)
FIELDS["random curl field"] = (J1, M1)

for name, c in CURRENTS.items():
    X, W = c
    T = toroidal_dubovik(c)
    half = 0.5 * (rdotJ(c)[:, None] * X).sum(0)
    quarter = -0.25 * mom(c, r2)
    sixth = (np.cross(X, np.cross(X, W))).sum(0) / 6
    planes = 0.25 * (np.cross(X, np.cross(X, W))).sum(0)
    check(f"(b) ½∫(r·J)r = −¼∫r²J = (1/6)∫r×(r×J) = Dubovik T, {name}",
          max(rel(half, T), rel(quarter, T), rel(sixth, T)) < 1e-12, f"max rel {max(rel(half, T), rel(quarter, T), rel(sixth, T)):.1e}")
    check(f"(a) Planes ¼∫r×(r×J) = (3/2)T, {name}", rel(planes, 1.5 * T) < 1e-12, f"ratio {planes @ T / (T @ T):.12f}")
    for p in range(4):
        lhs = ((r2(X) ** p * rdotJ(c))[:, None] * X).sum(0)
        rhs = -mom(c, lambda X: r2(X) ** (p + 1)) / (2 * p + 2)
        sc = np.abs((r2(X) ** p * np.abs(rdotJ(c)))[:, None] * np.abs(X)).sum()
        check(f"radii ∫r^{2*p}(r·J)r = −∫r^{2*p+2}J/{2*p+2}, {name}", np.max(np.abs(lhs - rhs)) < 1e-12 * sc)
        s = (r2(X) ** p * rdotJ(c)).sum()
        check(f"scalar ∫r^{2*p}(r·J) = 0, {name}", abs(s) < 1e-12 * np.abs(r2(X) ** p * rdotJ(c)).sum())

J, M = FIELDS["random curl field"]
T = 0.1 * gmom(GR * (GR * J).sum(-1, keepdims=True) - 2 * GRID_R2[..., None] * J, np.ones(GRID_R2.shape))
check("(b) ½∫r×M = Dubovik T for J = ∇×M (random curl field)",
      rel(0.5 * gmom(gcross(M), np.ones(GRID_R2.shape)), T) < 1e-9)
check("(b) −¼∫r²J = Dubovik T (random curl field)", rel(-0.25 * gmom(J, GRID_R2), T) < 1e-9)
check("(a) ¼∫r×(r×J) = (3/2)T (random curl field)",
      rel(0.25 * gmom(gcross(gcross(J)), np.ones(GRID_R2.shape)), 1.5 * T) < 1e-9)

# General order: ∫r^{2p+2} J·∇h = −(2p+2)∫r^{2p}(r·J) h for harmonic h = Re[(a·r)^L], a·a = 0.
for name, c in CURRENTS.items():
    X, W = c
    for L in range(1, 6):
        u = rng.normal(size=3); u /= np.linalg.norm(u)
        v = np.cross(u, rng.normal(size=3)); v /= np.linalg.norm(v)
        a = u + 1j * v
        ar = X @ a
        h = np.real(ar ** L)
        gradh = np.real(L * ar[:, None] ** (L - 1) * a[None, :])
        for p in range(3):
            lhs = (r2(X) ** (p + 1) * (W * gradh).sum(1)).sum()
            rhs = -(2 * p + 2) * (r2(X) ** p * rdotJ(c) * h).sum()
            sc = np.abs(r2(X) ** p * rdotJ(c) * h).sum()
            if not abs(lhs - rhs) < 1e-11 * sc:
                check(f"spherical one-term reduction L={L} p={p}, {name}", False, f"{lhs} vs {rhs}")
                break
        else:
            continue
        break
    else:
        check(f"spherical one-term reduction (harmonic weights, L = 1..5, p = 0..2), {name}", True)

# Cartesian toroidal quadrupole and octupole (Dumitriu-Vrejoiu 2005; Gurvitz et al. 2019 Table 2).
def stf3(S):
    tr = np.einsum("ijj->i", S)
    d = np.eye(3)
    return S - 0.2 * (np.einsum("ij,k->ijk", d, tr) + np.einsum("ik,j->ijk", d, tr) + np.einsum("jk,i->ijk", d, tr))


def tq_gurvitz(X, W):
    rj, rr = (X * W).sum(1), r2(X)
    t = 4 * np.einsum("n,nj,nk->jk", rj, X, X) + 2 * (rj * rr).sum() * np.eye(3) \
        - 5 * (np.einsum("n,nj,nk->jk", rr, X, W) + np.einsum("n,nj,nk->jk", rr, W, X))
    return t / 42


def to_gurvitz(X, W):
    rj, rr = (X * W).sum(1), r2(X)
    t = 35 * np.einsum("n,nj,nk,nl->jkl", rj, X, X, X)
    t -= 20 * (np.einsum("n,nl,nk,nj->jkl", rr, W, X, X) + np.einsum("n,nk,nl,nj->jkl", rr, W, X, X)
               + np.einsum("n,nj,nl,nk->jkl", rr, W, X, X))
    vec = (rj * rr)[:, None] * X + 4 * (rr ** 2)[:, None] * W
    vp = vec.sum(0)
    d = np.eye(3)
    t += np.einsum("lj,kp,p->jkl", d, d, vp) + np.einsum("lk,jp,p->jkl", d, d, vp) + np.einsum("lp,jk,p->jkl", d, d, vp)
    return t / 300


open_rand = (rng.normal(size=(500, 3)), rng.normal(size=(500, 3)))  # not divergence-free
check("Gurvitz T^(Qe) (1/42) is traceless for any J", abs(np.trace(tq_gurvitz(*open_rand))) < 1e-10 * np.abs(tq_gurvitz(*open_rand)).max())
check("Gurvitz T^(Oe) (1/300) is traceless for any J",
      np.abs(np.einsum("jjl->l", to_gurvitz(*open_rand))).max() < 1e-10 * np.abs(to_gurvitz(*open_rand)).max())
for name, c in CURRENTS.items():
    X, W = c
    rj = rdotJ(c)
    q2 = np.einsum("n,nj,nk->jk", rj, X, X)
    check(f"static T^(Qe) = ⅓∫(r·J)x x = −⅙∫r²(xJ+Jx), {name}",
          max(rel(tq_gurvitz(X, W), q2 / 3),
              rel(tq_gurvitz(X, W), -(np.einsum("n,nj,nk->jk", r2(X), X, W) + np.einsum("n,nj,nk->jk", r2(X), W, X)) / 6)) < 1e-12)
    q3 = stf3(np.einsum("n,nj,nk,nl->jkl", rj, X, X, X))
    check(f"static T^(Oe) = ¼ STF∫(r·J)x x x, {name}", rel(to_gurvitz(X, W), q3 / 4) < 1e-12)

# Krynytskyi & Rovenchak (arXiv:1804.10494, eqs. 56-57), any order n: T = ∫ρ_t x^{⊗n},
# ρ_t = −div{2(J·r)r − (n+3)r²J}/(2(n+1)(2n+3)). Contracted with a null a^{⊗n} (its STF part):
# static value (1/(n+1))∫(r·J)(a·r)^n, the Cartesian form of Dubovik-Tugushev (B.4) / Szmytkowski (C.13).
for name, c in CURRENTS.items():
    X, W = c
    worst = 0.0
    for n in range(1, 7):
        u = rng.normal(size=3); u /= np.linalg.norm(u)
        v = np.cross(u, rng.normal(size=3)); v /= np.linalg.norm(v)
        a = u + 1j * v
        ar, aJ, rj = X @ a, W @ a, rdotJ(c)
        kr = (2 * rj * n * ar ** n - (n + 3) * r2(X) * n * ar ** (n - 1) * aJ).sum() / (2 * (n + 1) * (2 * n + 3))
        one = (rj * ar ** n).sum() / (n + 1)
        worst = max(worst, abs(kr - one) / np.abs(rj * ar ** n).sum())
    check(f"Krynytskyi-Rovenchak T^(n), n = 1..6, static STF part = (1/(n+1))∫(r·J)x^n, {name}", worst < 1e-12, f"{worst:.1e}")

# Rank count: order-N moments ∫x^{⊗N}⊗J of closed currents span N(N+2) dims (open: 3·C(N+2,2)).
loops = [polygon(rng.uniform(-1, 1, (7, 3))) for _ in range(120)]
segs = [(rng.uniform(-1, 1, (1, 3)), rng.uniform(-1, 1, (1, 3))) for _ in range(120)]


def moment_vec(c, n):
    X, W = c
    t = W
    for _ in range(n):
        t = np.einsum("n...,ni->n...i", t, X)
    return t.sum(0).ravel()


for n in range(1, 5):
    def rank(cs):
        s = np.linalg.svd(np.array([moment_vec(c, n) for c in cs]), compute_uv=False)
        return int((s > 1e-9 * s[0]).sum())
    rc, ro = rank(loops), rank(segs)
    check(f"counting: order {n} moments of closed currents have rank {n*(n+2)} (open {3*(n+1)*(n+2)//2})",
          rc == n * (n + 2) and ro == 3 * (n + 1) * (n + 2) // 2, f"closed {rc}, open {ro}")

# ================================================================ hierarchy J = (∇×)^n K
def h_formula(n, J):
    k = (n - 1) // 2
    if n % 2:  # n = 2k+1: axial, ∫r^{2k}(r×J)
        c = (-1) ** k * 3 / (2 ** (k + 1) * math.factorial(k) * math.prod(range(1, 2 * k + 4, 2)))
        return c * gmom(gcross(J), GRID_R2 ** k)
    k = (n - 2) // 2  # n = 2k+2: polar, ∫r^{2k+2}J
    return (-1) ** (k + 1) * 1.5 / math.factorial(2 * k + 3) * gmom(J, GRID_R2 ** (k + 1))


def naive(n, J):  # 2^-n ∫(r×)^n J, the iterated-density (Planes / Prosandeev-Bellaiche) convention
    F = J
    for _ in range(n):
        F = gcross(F)
    return gmom(F, np.ones(GRID_R2.shape)) / 2 ** n


def ratio_formula(n):
    k = (n - 1) // 2 if n % 2 else (n - 2) // 2
    if n % 2:
        return math.factorial(k) * math.prod(range(1, 2 * k + 4, 2)) / (3 * 2 ** k)
    return (2 * k + 3) * math.factorial(2 * k + 3) / (3 * 2 ** (2 * k + 1) * (2 * k + 2))


hK = gmom(K_FIELD, np.ones(GRID_R2.shape))
J = K_FIELD
coeffs = {1: "1/2", 2: "-1/4", 3: "-1/20", 4: "1/80", 5: "1/560", 6: "-1/3360"}
for n in range(1, 7):
    J = curl(J)
    hf = h_formula(n, J)
    check(f"hierarchy h_{n} = ∫K = ({coeffs[n]})∫{'r^%d(r×J)' % (n - 1) if n % 2 else 'r^%dJ' % n}", rel(hf, hK) < 1e-7, f"rel {rel(hf, hK):.1e}")
    nv = naive(n, J)
    check(f"iterated ½r× convention / h_{n} = {ratio_formula(n):g}", rel(nv, ratio_formula(n) * hK) < 1e-7,
          f"measured {nv @ hK / (hK @ hK):.9f}")

# A rig: 16 small tori on a ring, axes tangential. h_3 from the one-term formula = ½Σ c_i × T_i.
ring = []
centres, tvecs = [], []
for i in range(16):
    ph = 2 * math.pi * i / 16
    cen = 2.0 * np.array([math.cos(ph), math.sin(ph), 0.0])
    ax = np.array([-math.sin(ph), math.cos(ph), 0.0])
    c = torus_winding(turns=24, m=40, big_r=0.3, a=0.1, centre=cen, axis=ax)
    ring.append(c)
    centres.append(cen)
    Xc, Wc = c
    tvecs.append(toroidal_dubovik((Xc - cen, Wc)))
rig = join(*ring)
h3_rig = 0.5 * np.cross(np.array(centres), np.array(tvecs)).sum(0)
h3_formula = -(r2(rig[0])[:, None] * np.cross(rig[0], rig[1])).sum(0) / 20
check("rig of 16 tori: −(1/20)∫r²(r×J) = ½Σ cᵢ×Tᵢ (the crate's hypertoroidal_moment)", rel(h3_formula, h3_rig) < 1e-10,
      f"rel {rel(h3_formula, h3_rig):.1e}")

# ================================================================ Q3: Planes eq. 28 vs eq. 29
# Coulomb-gauge A = ∫J/|R−r| (μ0/4π dropped) of an untilted torus winding at the origin
# (m = 0, no magnetic quadrupole): the leading term is ∇(t·∇)(1/R) with t = Dubovik T, not 3/2 T.
sol0 = torus_winding(turns=24, m=40)
T0 = toroidal_dubovik(sol0)
for Rd in (50.0, 100.0):
    Rv = Rd * np.array([0.48, 0.36, 0.8])
    A = (sol0[1] / np.linalg.norm(Rv - sol0[0], axis=1)[:, None]).sum(0)
    lead = 3 * (T0 @ Rv) * Rv / Rd ** 5 - T0 / Rd ** 3
    fit = A @ lead / (lead @ lead)
    check(f"Planes eq. 28 needs t = T: far A fits ∇(t·∇)(1/R) with t/T = {fit:.6f} at R = {Rd:g}", abs(fit - 1) < 2e-3 and abs(fit - 1.5) > 0.4)

# ================================================================ published numbers
# Luk'yanchuk et al. (arXiv:2406.19728) model Hopfion in a sphere of radius R = 1:
# T = ½∫r×P = 0.96 αt R⁴ and G = ¼∫r×(r×P) = 1.19 αp R⁵ ("Prosandeev and Bellaiche, 2008").
nr = 600
r = (np.arange(nr) + 0.5) / nr
th = (np.arange(nr) + 0.5) / nr * math.pi
rr, tt = np.meshgrid(r, th, indexing="ij")
dV = 2 * math.pi * rr ** 2 * np.sin(tt) * (1 / nr) * (math.pi / nr)
eta = np.sin(math.pi * rr / 2) * np.sin(tt)
T_hopf = 0.5 * (rr * np.sin(tt) * eta * dV).sum()
# P_pol = ∇×(ζ φ̂), ζ = sin(πr)/r · sin θ: P_r = 2 sin(πr) cos θ / r², P_θ = −π cos(πr) sin θ / r
Pr = 2 * np.sin(math.pi * rr) * np.cos(tt) / rr ** 2
Pt = -math.pi * np.cos(math.pi * rr) * np.sin(tt) / rr
Pz = Pr * np.cos(tt) - Pt * np.sin(tt)
G_hopf = 0.25 * ((rr * Pr * rr * np.cos(tt) - rr ** 2 * Pz) * dV).sum()
h2_hopf = -0.25 * (rr ** 2 * Pz * dV).sum()
check(f"Luk'yanchuk Hopfion T = 0.96 αt R⁴ reproduced ({T_hopf:.4f})", round(T_hopf, 2) == 0.96)
check(f"Luk'yanchuk Hopfion G = 1.19 αp R⁵ reproduced ({G_hopf:.4f})", round(G_hopf, 2) == 1.19)
check(f"Hopfion G / (−¼∫r²P) = 3/2 (solenoidal P): {G_hopf / h2_hopf:.6f}", abs(G_hopf / h2_hopf - 1.5) < 1e-4)

# Nemkov et al. 2018 eq. 13: two nested poloidal tori with zero toroidal dipole have no static
# toroidal octupole iff R² − a²/4 is equal for both; one-term form ∫(r·J) r³P₃(cos θ).
def octupole(c):
    X = c[0]
    z = X[:, 2]
    return (rdotJ(c) * (5 * z ** 3 - 3 * z * r2(X)) / 2).sum()


def nested(R_out, a_out):
    """Inner torus (R = 1, a = 0.3, I = 1) plus an outer one whose current zeroes the total T."""
    inner = torus_winding(turns=48, m=120, big_r=1.0, a=0.3)
    outer = torus_winding(turns=48, m=120, big_r=R_out, a=a_out)
    i2 = -toroidal_dubovik(inner)[2] / toroidal_dubovik(outer)[2]
    both = join(inner, (outer[0], i2 * outer[1]))
    return octupole(both) / abs(octupole(inner)), toroidal_dubovik(both)[2], to_gurvitz(*both)[2, 2, 2]


a_out = 0.6
R_match = math.sqrt(1.0 - 0.3 ** 2 / 4 + a_out ** 2 / 4)
o3m, t1m, gm = nested(R_match, a_out)
o3x, _, gx = nested(R_match + 0.1, a_out)
check(f"Nemkov eq. 13: matched radii leave no toroidal octupole (rel {o3m:.1e}, Gurvitz T^(Oe)_zzz ratio {gm / gx:.1e}); mismatched rel {o3x:.2e}",
      abs(o3m) < 1e-3 and abs(gm / gx) < 1e-2 and abs(t1m) < 1e-12 and abs(o3x) > 1e-2)

# Gurvitz et al. 2019 Table 2 prints T^(2e) with 1/10; with their eq. (23), p + (ik/c)T + (ik³/c)T^(2e),
# Alaee et al.'s exact dipole −(1/iω){∫J j0 + (k²/2)∫[3(r·J)r − r²J] j2/(kr)²} needs 1/280.
def j0(x):
    return np.sin(x) / x


def j2_over_x2(x):
    return ((3 / x ** 3 - 1 / x) * np.sin(x) - 3 * np.cos(x) / x ** 2) / x ** 2


X, W = knot
rn = np.sqrt(r2(X))
F = lambda k: (j0(k * rn)[:, None] * W).sum(0) + 0.5 * k ** 2 * (
    (j2_over_x2(k * rn))[:, None] * (3 * rdotJ(knot)[:, None] * X - r2(X)[:, None] * W)).sum(0)
Tk = toroidal_dubovik(knot)
X4 = ((3 * r2(X) ** 2)[:, None] * W - 2 * (r2(X) * rdotJ(knot))[:, None] * X).sum(0) / 280
g = {k: (F(k) - k ** 2 * Tk) / k ** 4 for k in (0.2, 0.1)}  # F(0) = ∫J = 0
c4 = (4 * g[0.1] - g[0.2]) / 3  # Richardson: remove the O(k²) remainder
check(f"Alaee p: k² coefficient is Dubovik T (1/10) and k⁴ coefficient is (1/280)∫[3r⁴J − 2r²(r·J)r] (rel {rel(c4, X4):.1e})",
      rel(c4, X4) < 1e-4 and rel(c4, 28 * X4) > 0.9)
check("static (1/280)∫[3r⁴J − 2r²(r·J)r] = (1/80)∫r⁴J = h_4", rel(X4, mom(knot, lambda X: r2(X) ** 2) / 80) < 1e-12)

# Generating functions (closed J): the exact dipole brackets are one-term Bessel moments,
# P(k) = (3/2)∫J j0(kr) = (3/2)k∫(r·J) r̂ j1(kr) (ρ = 0 form of the exact electric coefficient) and
# M(k) = (3/2)∫(r×J) j1(kr)/(kr); their Taylor coefficients are the ladder h_2, h_4, … and h_1, h_3, ….
def j1(x):
    return np.sin(x) / x ** 2 - np.cos(x) / x


for name, c in CURRENTS.items():
    X, W = c
    rn = np.sqrt(r2(X))
    Fc = lambda k: (j0(k * rn)[:, None] * W).sum(0) + 0.5 * k ** 2 * (
        (j2_over_x2(k * rn))[:, None] * (3 * rdotJ(c)[:, None] * X - r2(X)[:, None] * W)).sum(0)
    worst = 0.0
    for k in (0.7, 1.3, 2.9):
        one = 1.5 * (j0(k * rn)[:, None] * W).sum(0)
        jac = 1.5 * k * ((j1(k * rn) * rdotJ(c) / rn)[:, None] * X).sum(0)
        worst = max(worst, rel(Fc(k), one), rel(jac, one))
    check(f"exact electric-dipole bracket = (3/2)∫J j0(kr) = (3/2)k∫(r·J)r̂ j1(kr) at kR ~ 1-5, {name}", worst < 1e-10,
          f"max rel {worst:.1e}")

print(f"\n{'ALL PASS' if not FAILED else str(len(FAILED)) + ' FAILED'}")
sys.exit(1 if FAILED else 0)
