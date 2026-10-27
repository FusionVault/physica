"""Checks for observability.md: what can and cannot be measured about toroidal dipoles,
how sharp the Mie a1 zero is, and exact properties of supertoroidal pulses.

Run: python3 -P observability_checks.py   (needs numpy; about 1 minute)
Every number quoted in observability.md is asserted here; exits 1 on any failure.
Units: SI where stated; the pulse integrals use c = mu0 = eps0 = 1 and lengths in units of q1.
"""
import cmath
import math
import sys

import numpy as np

C = 299792458.0
FAILS = []


def check(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))
    if not ok:
        FAILS.append(name)


def close(a, b, rel):
    return abs(a - b) <= rel * max(abs(a), abs(b))


# ---------------------------------------------------------------------------------------------
# 1. Single frequency: a non-radiating family J = (curl curl - k^2) F carries any (p, T).
#    F = z_hat * phi(r) with phi a smooth bump of radius a. For a radial phi,
#    J = (z u'(r)/r) r_vec + z_hat w(r),  u = phi'/r,  w = -phi'' - phi'/r - k^2 phi.
# ---------------------------------------------------------------------------------------------
GL_X, GL_W = np.polynomial.legendre.leggauss(400)


def radial_quad(f, a):
    r = 0.5 * a * (GL_X + 1.0)
    return 0.5 * a * np.sum(GL_W * f(r))


def bump(kind, a):
    """phi, phi', phi'' for phi = (1-s)^4 (kind 0) or s (1-s)^4 (kind 1), s = r^2/a^2."""
    def d(r):
        s = r * r / a**2
        ds, d2s = 2 * r / a**2, 2 / a**2
        if kind == 0:
            g, g1, g2 = (1 - s) ** 4, -4 * (1 - s) ** 3, 12 * (1 - s) ** 2
        else:
            g = s * (1 - s) ** 4
            g1 = (1 - s) ** 4 - 4 * s * (1 - s) ** 3
            g2 = -8 * (1 - s) ** 3 + 12 * s * (1 - s) ** 2
        return g, g1 * ds, g2 * ds * ds + g1 * d2s
    return d


def sph_j0(x):
    return np.sinc(x / np.pi)


def sph_j2_over_x2(x):
    # j2(x)/x^2, with its series near 0 (limit 1/15)
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    small = np.abs(x) < 1e-2
    xs = x[small]
    out[small] = 1 / 15 - xs**2 / 210 + xs**4 / 7560
    xl = x[~small]
    out[~small] = ((3 / xl**2 - 1) * np.sin(xl) / xl - 3 * np.cos(xl) / xl**2) / xl**2
    return out


def nr_moments(kind, a, k):
    """Return (int J_z, T_z, P_exact_z [the bracket of Alaee's exact dipole], Phi0, Phi2)."""
    d = bump(kind, a)

    def parts(r):
        phi, p1, p2 = d(r)
        u = np.where(r > 0, p1 / np.where(r > 0, r, 1), 0)
        # u' = phi''/r - phi'/r^2
        up = np.where(r > 0, p2 / np.where(r > 0, r, 1) - p1 / np.where(r > 0, r * r, 1), 0)
        w = -p2 - u - k * k * phi
        return phi, u, up, w

    four_pi = 4 * math.pi
    jz = four_pi * radial_quad(lambda r: r * r * (r * parts(r)[2] / 3 + parts(r)[3]), a)
    tz = four_pi * radial_quad(lambda r: r * r * (-(r**3) * parts(r)[2] / 3 - 5 * r * r * parts(r)[3] / 3), a) / 10
    pex = four_pi * radial_quad(
        lambda r: r * r * ((r * parts(r)[2] / 3 + parts(r)[3]) * sph_j0(k * r)
                           + 0.5 * k * k * (2 * r**3 * parts(r)[2] / 3) * sph_j2_over_x2(k * r) * 1.0),
        a)
    phi0 = four_pi * radial_quad(lambda r: r * r * parts(r)[0], a)
    phi2 = four_pi * radial_quad(lambda r: r**4 * parts(r)[0], a)
    return jz, tz, pex, phi0, phi2


print("== 1. single-frequency non-radiating family")
a = 1.0
for kind in (0, 1):
    for ka in (0.3, 1.0, 2.5):
        k = ka / a
        jz, tz, pex, phi0, phi2 = nr_moments(kind, a, k)
        # exact electric dipole bracket: (k^2/2)*[...]*j2/(kr)^2 has (kr)^2 in the denominator -> k^2 * r^2 * j2(kr)/(kr)^2
        check(f"exact a_E1 = 0 (kind {kind}, ka {ka})", abs(pex) < 1e-12 * abs(jz), f"{abs(pex)/abs(jz):.1e}")
        check(f"int J = -k^2 Phi0 (kind {kind}, ka {ka})", close(jz, -k * k * phi0, 1e-12))
        check(f"T = Phi0 + k^2 Phi2/6 (kind {kind}, ka {ka})", close(tz, phi0 + k * k * phi2 / 6, 1e-12))
        # p = (i/w) int J and (ik/c)T, in units omega = c k: p_eff/(i/w) = int J + k^2 T
        resid = abs(jz + k * k * tz) / abs(jz)
        check(f"p + (ik/c)T = O(k^4): residual k^2<r^2>/6 (kind {kind}, ka {ka})",
              close(resid, k * k * phi2 / phi0 / 6, 1e-10), f"{resid:.4g}")

# Two bumps reach any (int J, T) pair: e.g. p = 0 with T = 1 (an exactly non-radiating "pure toroid").
k = 1.0
m0 = nr_moments(0, a, k)
m1 = nr_moments(1, a, k)
M = np.array([[m0[0], m1[0]], [m0[1], m1[1]]])
c0, c1 = np.linalg.solve(M, [0.0, 1.0])
check("NR source with int J = 0, T = 1 exists (ka = 1)", abs(c0 * m0[2] + c1 * m1[2]) < 1e-12 and close(c0 * m0[1] + c1 * m1[1], 1.0, 1e-12))

# Reception = radiation (reciprocity): int J . eps e^{i k n.r} = 0 for transverse eps; but a uniform
# (non-vacuum) field E = z_hat couples with int J_z = -k^2 Phi0.
def plane_wave_coupling(kind, a, k, nvec, eps, nr=80, nmu=80, nph=64):
    d = bump(kind, a)
    xr, wr = np.polynomial.legendre.leggauss(nr)
    xm, wm = np.polynomial.legendre.leggauss(nmu)
    r = 0.5 * a * (xr + 1)
    wr = 0.5 * a * wr
    ph = 2 * np.pi * (np.arange(nph) + 0.5) / nph
    R, MU, PH = np.meshgrid(r, xm, ph, indexing="ij")
    W = (wr[:, None, None] * wm[None, :, None]) * (2 * np.pi / nph) * R**2
    st = np.sqrt(1 - MU**2)
    X, Y, Z = R * st * np.cos(PH), R * st * np.sin(PH), R * MU
    phi, p1, p2 = d(R)
    up = p2 / R - p1 / R**2
    w = -p2 - p1 / R - k * k * phi
    fac = Z * up / R
    J = [fac * X, fac * Y, fac * Z + w]
    ph_ = np.exp(1j * k * (nvec[0] * X + nvec[1] * Y + nvec[2] * Z))
    return sum(eps[i] * np.sum(W * J[i] * ph_) for i in range(3)), np.sum(W * J[2])


th = 0.7
nvec = (math.sin(th), 0.0, math.cos(th))
for eps in [(math.cos(th), 0.0, -math.sin(th)), (0.0, 1.0, 0.0)]:
    cpl, jz_num = plane_wave_coupling(0, a, 1.3, nvec, eps)
    check(f"NR source does not couple to a plane wave, eps={tuple(round(e, 3) for e in eps)}", abs(cpl) < 1e-10 * abs(jz_num), f"{abs(cpl)/abs(jz_num):.1e}")
check("but couples to a uniform field: int J_z = -k^2 Phi0", close(jz_num, -1.3**2 * nr_moments(0, a, 1.3)[3], 1e-6))

# ---------------------------------------------------------------------------------------------
# 2. Many frequencies: a fixed pattern j(r) with known drive is identified; a pattern that moves at
#    O(k^2) is not. Pattern: 24-turn toroidal solenoid (R = 1, a = 0.3, 1 A) plus a 0.2 m z-segment
#    on the axis carrying 1 A (an open wire, so int j != 0).
# ---------------------------------------------------------------------------------------------
print("== 2. band data and the fixed-pattern hypothesis")


def solenoid(turns=24, i=1.0, big_r=1.0, aa=0.3, m=200):
    segs = []
    for j in range(turns):
        ph = 2 * math.pi * (j + 0.5) / turns
        def pt(kk):
            t = 2 * math.pi * kk / m
            rho = big_r + aa * math.cos(t)
            return np.array([rho * math.cos(ph), rho * math.sin(ph), aa * math.sin(t)])
        for kk in range(m):
            p, q = pt(kk), pt(kk + 1)
            segs.append((0.5 * (p + q), i * (q - p)))
    return segs


def zseg(i=1.0, ell=0.2, m=400):
    zs = -ell / 2 + ell * (np.arange(m) + 0.5) / m
    return [(np.array([0.0, 0.0, z]), np.array([0.0, 0.0, i * ell / m])) for z in zs]


def phi_exact(segs, k):
    """z-component of int J j0(kr) + (k^2/2) int [3(r.J)r - r^2 J] j2(kr)/(kr)^2 (= -i w a_E1-type bracket)."""
    out = 0.0
    for r, idl in segs:
        rr = float(np.linalg.norm(r))
        x = np.array([k * rr])
        out += idl[2] * sph_j0(x)[0] + 0.5 * k * k * (3 * np.dot(r, idl) * r[2] - rr * rr * idl[2]) * sph_j2_over_x2(x)[0]
    return out


def tor(segs):
    t = np.zeros(3)
    for r, idl in segs:
        t += 0.1 * (np.dot(r, idl) * r - 2 * np.dot(r, r) * idl)
    return t


sol, seg = solenoid(), zseg()
both = sol + seg
T_sol, T_seg, T_both = tor(sol)[2], tor(seg)[2], tor(both)[2]
check("solenoid |T| = N I V/(4 pi) = 3.393 A m^3 (to the 200-gon's 1.6e-4)", close(abs(T_sol), 24 * math.pi * 0.3**2 / 2, 2e-4), f"{T_sol:.5f}")
check("segment T = -I l^3/120 (midpoint rule)", close(T_seg, -0.2**3 / 120, 1e-5))
ks = [0.05, 0.10, 0.15]
vals = [phi_exact(both, kk) for kk in ks]
A = np.array([[1, kk**2, kk**4] for kk in ks])
c = np.linalg.solve(A, vals)
check("fixed pattern: band fit gives int j = 0.2 A m", close(c[0], 0.2, 1e-7), f"{c[0]:.10f}")
check("fixed pattern: band fit gives T of the whole pattern", close(c[1], T_both, 1e-4), f"{c[1]:.6f} vs {T_both:.6f}")


def phi_alt(k):
    # the segment alone, its current rising as I(k) = 1 + k^2 T_sol/l: no toroid at all
    return phi_exact(zseg(i=1.0 + k * k * T_sol / 0.2), k)


mis = lambda kk: abs(phi_alt(kk) - phi_exact(both, kk)) / abs(phi_exact(both, kk) - phi_exact(both, 0.0))
d1, d2 = mis(0.025), mis(0.05)
# mismatch relative to the k^2 (toroidal) part of a_E1 itself
check("segment with O(k^2) current change mimics the toroid: mismatch 3e-4 of the toroidal term at kR = 0.05", d2 < 3.5e-4, f"{d2:.2e}")
check("and the mismatch grows as k^2 relative to that term (ratio 4 when k doubles)", close(d2 / d1, 4.0, 0.05), f"{d2/d1:.3f}")
check("while its toroidal moment is only the segment's", close(tor(zseg(i=1.0 + 0.05**2 * T_sol / 0.2))[2], T_seg * (1 + 0.05**2 * T_sol / 0.2), 1e-9))

# An overlapping probe does see T: a closed current in a uniform external current density J0 z_hat
# (A_ext = -mu0 J0 rho^2/4 z_hat) has interaction energy sum(I dl . A_ext) = mu0 J0 T_z exactly.
U_int = sum(-0.25 * (r[0] ** 2 + r[1] ** 2) * idl[2] for r, idl in sol)  # per unit mu0 J0
check("static toroid in a uniform current density: U = mu0 T . J_ext", close(U_int, T_sol, 1e-9), f"{U_int:.6f} vs {T_sol:.6f}")

# ---------------------------------------------------------------------------------------------
# 3. The Mie a1 zero: sensitivity and tolerance.
# ---------------------------------------------------------------------------------------------
print("== 3. Mie a1 zero")


def psi1(z): return cmath.sin(z) / z - cmath.cos(z)
def dpsi1(z): return cmath.cos(z) / z - cmath.sin(z) / z**2 + cmath.sin(z)
def chi1(z): return -cmath.cos(z) / z - cmath.sin(z)
def dchi1(z): return cmath.sin(z) / z + cmath.cos(z) / z**2 - cmath.cos(z)
def psi2(z): return (3 / z**2 - 1) * cmath.sin(z) - 3 * cmath.cos(z) / z
def chi2(z): return (3 / z**2 - 1) * cmath.cos(z) + 3 * cmath.sin(z) / z
def dpsi2(z): return psi1(z) - 2 * psi2(z) / z
def dchi2(z): return chi1(z) - 2 * chi2(z) / z


def mie_a(m, x, order=1):
    ps, dps, ch, dch = (psi1, dpsi1, chi1, dchi1) if order == 1 else (psi2, dpsi2, chi2, dchi2)
    xi, dxi = ps(x) + 1j * ch(x), dps(x) + 1j * dch(x)
    return (m * ps(m * x) * dps(x) - ps(x) * dps(m * x)) / (m * ps(m * x) * dxi - xi * dps(m * x))


def mie_b1(m, x):
    xi, dxi = psi1(x) + 1j * chi1(x), dpsi1(x) + 1j * dchi1(x)
    return (psi1(m * x) * dpsi1(x) - m * psi1(x) * dpsi1(m * x)) / (psi1(m * x) * dxi - m * xi * dpsi1(m * x))


def zero_d(n, lam):
    """Mirror of magnetics::toroidal::mie_a1_zero_diameter (same scan and bisection)."""
    ps = lambda z: math.sin(z) / z - math.cos(z)
    dp = lambda z: math.cos(z) / z - math.sin(z) / z**2 + math.sin(z)
    g = lambda x: n * ps(n * x) * dp(x) - ps(x) * dp(n * x)
    step = 0.01 / max(n, 1.0)
    lo = step
    while g(lo) * g(lo + step) > 0:
        lo += step
        if lo > 100:
            return float("nan")
    hi = lo + step
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if g(lo) * g(mid) <= 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi) * lam / math.pi


def bisect(f, lo, hi):
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if f(lo) * f(mid) <= 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


# Values printed by `cargo run` of a scratch binary calling magnetics::toroidal::mie_a1_zero_diameter
RUST_MIE = {3.9999: 204.033898480970, 4.0: 204.028352425795, 4.0001: 204.022806691133,
            4.077: 199.850804267955, 4.112: 198.011206711165, 4.045: 201.564702724423}
for n, v in RUST_MIE.items():
    check(f"Python mirror = Rust mie_a1_zero_diameter(n = {n}, 550)", abs(zero_d(n, 550.0) - v) < 1e-9)

lam = 550.0
d0 = zero_d(4.0, lam)
x0 = math.pi * d0 / lam
check("d0 = 204.028 nm (Miroshnichenko et al. ~204 nm)", abs(d0 - 204.028) < 1e-3, f"{d0:.4f}")
check("a1 vanishes there", abs(mie_a(4.0, x0)) < 1e-13)
dd_dn = (RUST_MIE[4.0001] - RUST_MIE[3.9999]) / 2e-4
check("dd/dn = -55.46 nm per unit index", abs(dd_dn + 55.46) < 0.01, f"{dd_dn:.3f}")
check("non-dispersive dd/dlambda = x0/pi = 0.3710", abs(d0 / lam - 0.37096) < 1e-5)
# Silicon (Green 2008, 300 K): n(540) = 4.112, n(550) = 4.077 + 0.027968i, n(560) = 4.045.
dd_dl_si = (RUST_MIE[4.045] * 560 / 550 - RUST_MIE[4.112] * 540 / 550) / 20
check("with Si dispersion dd/dlambda = 0.541", abs(dd_dl_si - 0.541) < 1e-3, f"{dd_dl_si:.4f}")

F = lambda d: abs(mie_a(4.0, math.pi * d / lam)) ** 2 - 1e-3
dplus, dminus = bisect(F, d0, d0 + 5) - d0, bisect(F, d0 - 5, d0) - d0
check("|a1|^2 < 1e-3 needs d within -1.020/+1.082 nm", abs(dplus - 1.082) < 1e-3 and abs(dminus + 1.020) < 1e-3, f"{dminus:.4f}, {dplus:.4f}")
Fn = lambda n: abs(mie_a(n, x0)) ** 2 - 1e-3
nplus, nminus = bisect(Fn, 4.0, 4.1) - 4.0, bisect(Fn, 3.9, 4.0) - 4.0
check("|a1|^2 < 1e-3 needs n within -0.0184/+0.0195", abs(nplus - 0.01951) < 1e-4 and abs(nminus + 0.01840) < 1e-4, f"{nminus:.5f}, {nplus:.5f}")
Fl = lambda L: abs(mie_a(4.0, math.pi * d0 / L)) ** 2 - 1e-3
lplus, lminus = bisect(Fl, lam, lam + 10) - lam, bisect(Fl, lam - 10, lam) - lam
check("|a1|^2 < 1e-3 needs lambda within -2.90/+2.76 nm (n fixed)", abs(lplus - 2.763) < 0.01 and abs(lminus + 2.901) < 0.01, f"{lminus:.3f}, {lplus:.3f}")
check("at the a1 zero |b1|^2 = 0.0762: the plane-wave far field is not zero", abs(abs(mie_b1(4.0, x0)) ** 2 - 0.0762) < 1e-4)
a2sq = abs(mie_a(4.0, x0, 2)) ** 2
check("at the a1 zero |a2|^2 = 0.0018 (residual under radially polarised light)", abs(a2sq - 0.00184) < 1e-5, f"{a2sq:.5f}")

# Absorption floor: min over d of |a1|^2 -> |da1/dn|^2 kappa^2 as kappa -> 0.
da_dn = abs((mie_a(4.0 + 1e-6, x0) - mie_a(4.0 - 1e-6, x0)) / 2e-6)
check("|da1/dn| = 1.6705 at the zero", abs(da_dn - 1.6705) < 1e-3, f"{da_dn:.4f}")
for kap in (1e-3, 3e-3):
    floor = min(abs(mie_a(4.0 + 1j * kap, math.pi * d / lam)) ** 2 for d in np.linspace(d0 - 1, d0 + 1, 4001))
    check(f"loss floor = |da1/dn|^2 kappa^2 (kappa = {kap})", close(floor, da_dn**2 * kap**2, 0.02), f"{floor:.3e} vs {da_dn**2*kap**2:.3e}")
si = 4.077 + 0.027968j
ds = np.linspace(195, 205, 10001)
vals = [abs(mie_a(si, math.pi * d / lam)) ** 2 for d in ds]
i = int(np.argmin(vals))
check("crystalline Si at 550 nm: floor 1.91e-3 > 1e-3, at d = 199.9 nm", abs(vals[i] - 1.906e-3) < 2e-6 and abs(ds[i] - 199.9) < 0.05, f"{vals[i]:.4e} at {ds[i]:.2f}")
kmax = math.sqrt(1e-3) / da_dn
check("a 1e-3 floor needs kappa < 0.019 (at n = 4)", abs(kmax - 0.019) < 5e-4, f"{kmax:.4f}")

# ---------------------------------------------------------------------------------------------
# 4. Supertoroidal pulses f = f0/[(q1 + i tau)(s + q2)^alpha], TE fields from the Hertz potential.
# ---------------------------------------------------------------------------------------------
print("== 4. supertoroidal pulses")
rng = np.random.default_rng(1)
rr, zz, tt = rng.uniform(0, 50, 20000), rng.uniform(-50, 50, 20000), rng.uniform(-50, 50, 20000)
for q1, q2 in ((1.0, 20.0), (0.3, 0.01)):
    s = rr**2 / (q1 + 1j * (zz - tt)) - 1j * (zz + tt)
    check(f"Re(s + q2) >= q2 (q1 = {q1}, q2 = {q2})", np.min((s + q2).real) >= q2 - 1e-12)


def f_si(r, z, t, q1, q2, al):
    tau, sig = z - C * t, z + C * t
    A = q1 + 1j * tau
    return 1 / (A * (r * r / A - 1j * sig + q2) ** al)


RUST_STLP = [((0.3, 0.2, 1e-9, 0.5), (2.205916646180394e-1, 2.473834112263143e-2)),
             ((2.0, -1.5, -3e-9, 20.0), (-3.791221069632554e-28, -3.435365403528666e-29)),
             ((5.0, 3.0, 0.0, 2.5), (1.016660961346763e-4, -1.614383353647747e-5))]
for (r, z, t, al), (re, im) in RUST_STLP:
    v = f_si(r, z, t, 1.0, 20.0, al)
    check(f"principal power = Rust supertoroidal_generating_fn at alpha = {al}", abs(v - complex(re, im)) < 1e-12 * abs(v))


def fields(r, z, t, q1, q2, al):
    """E_phi = d_t d_r f, H_r = d_r d_z f, H_z = -lap_perp f (c = mu0 = 1, f0 = 1)."""
    A, Q = q1 + 1j * (z - t), q2 - 1j * (z + t)
    D = r * r + A * Q
    E = -2 * al * 1j * r * ((al + 1) * A ** (al - 1) * (q1 + q2 - 2j * t) / D ** (al + 2) - (al - 1) * A ** (al - 2) / D ** (al + 1))
    Hr = 2 * al * 1j * r * ((al + 1) * A ** (al - 1) * (q2 - q1 - 2j * z) / D ** (al + 2) - (al - 1) * A ** (al - 2) / D ** (al + 1))
    Hz = 4 * al * A ** (al - 1) * (A * Q - al * r * r) / D ** (al + 2)
    return E, Hr, Hz


def hz_shen_eq6(r, z, t, q1, q2, al):
    A, Q = q1 + 1j * (z - t), q2 - 1j * (z + t)
    D = r * r + A * Q
    return -4 * al * (A ** (al - 1) * (r * r - al * A * Q) / D ** (al + 2) + (al - 1) * A ** (al - 2) * Q / D ** (al + 1))


# our fields come from f by finite differences
h = 1e-4
def f1(r, z, t, al): return 1 / ((1 + 1j * (z - t)) * (r * r / (1 + 1j * (z - t)) - 1j * (z + t) + 20.0) ** al)
for al in (1.0, 5.0):
    r, z, t = 0.7, 0.4, 0.3
    E, Hr, Hz = fields(r, z, t, 1.0, 20.0, al)
    fr = lambda r, z, t: (f1(r + h, z, t, al) - f1(r - h, z, t, al)) / (2 * h)
    E_fd = (fr(r, z, t + h) - fr(r, z, t - h)) / (2 * h)
    Hr_fd = (fr(r, z + h, t) - fr(r, z - h, t)) / (2 * h)
    lap = (f1(r + h, z, t, al) - 2 * f1(r, z, t, al) + f1(r - h, z, t, al)) / h**2 + fr(r, z, t) / r
    check(f"E_phi, H_r, H_z are the Hertz-potential derivatives (alpha = {al})",
          abs(E - E_fd) < 1e-5 * abs(E) and abs(Hr - Hr_fd) < 1e-5 * abs(Hr) and abs(Hz + lap) < 1e-4 * abs(Hz))
    divs = []
    for hz in (lambda *x: fields(*x)[2], hz_shen_eq6):
        dr = ((r + h) * fields(r + h, z, t, 1.0, 20.0, al)[1] - (r - h) * fields(r - h, z, t, 1.0, 20.0, al)[1]) / (2 * h * r)
        dz = (hz(r, z + h, t, 1.0, 20.0, al) - hz(r, z - h, t, 1.0, 20.0, al)) / (2 * h)
        divs.append(abs(dr + dz) / abs(dr))
    if al == 1.0:
        check("alpha = 1: Shen et al. eq. (6) agrees and div H = 0", divs[0] < 1e-6 and divs[1] < 1e-6)
    else:
        check("alpha = 5: our H_z gives div H = 0, Shen et al. eq. (6) as printed does not", divs[0] < 1e-6 and divs[1] > 1.0, f"{divs[0]:.1e} vs {divs[1]:.1f}")


def energy_closed(q1, q2, al):
    """Energy of the complex field, int (|E|^2 + |H|^2)/2 (c = mu0 = 1); the physical Re-field has half."""
    return math.pi**1.5 * math.gamma(al + 0.5) / math.gamma(al) * (al * q1 + q2) / (q1**3 * q2 ** (2 * al + 1))


def energy_momentum(q1, q2, al, t=0.0, nu=800, nR=4000):
    s = np.linspace(math.log(1e-4), math.log(1e6), nR)
    R, dsx = np.exp(s), s[1] - s[0]
    x, w = np.polynomial.legendre.leggauss(nu)
    U = P = 0.0
    for xi, wi in zip(x, w):
        E, Hr, Hz = fields(R * math.sqrt(1 - xi * xi), R * xi, t, q1, q2, al)
        U += wi * np.sum(0.5 * (abs(E) ** 2 + abs(Hr) ** 2 + abs(Hz) ** 2) * R**3) * dsx * 2 * math.pi
        P += wi * np.sum(-(E * np.conj(Hr)).real * R**3) * dsx * 2 * math.pi
    return U, P


for al, q1, q2, t, tol in ((1.0, 1.0, 3.0, 0.0, 1e-10), (2.5, 1.0, 3.0, 0.0, 1e-10), (2.5, 1.0, 3.0, 4.0, 1e-10),
                           (3.0, 1.0, 3.0, 0.0, 1e-10), (5.0, 1.0, 20.0, 0.0, 1e-10)):
    U, P = energy_momentum(q1, q2, al, t)
    beta = (q2 - al * q1) / (q2 + al * q1)
    check(f"U closed form (alpha {al}, q2/q1 {q2/q1:g}, t {t})", close(U, energy_closed(q1, q2, al), tol), f"{U/energy_closed(q1, q2, al)-1:.1e}")
    check(f"c P_z / U = (q2 - alpha q1)/(q2 + alpha q1) = {beta:.4f}", abs(P / U - beta) < 1e-10, f"{P/U:.12f}")
U, _ = energy_momentum(1.0, 3.0, 0.5, nu=3000)
check("alpha = 0.5 (below Shen's alpha >= 1): energy finite and equal to the closed form", close(U, energy_closed(1.0, 3.0, 0.5), 1e-6), f"{U/energy_closed(1.0,3.0,0.5)-1:.1e}")

# Boost covariance: f(r, lam tau', sigma'/lam; q1, q2) = lam^(alpha-1) f(r, tau', sigma'; q1/lam, lam q2)
def f_ts(r, tau, sig, q1, q2, al):
    A = q1 + 1j * tau
    return 1 / (A * (r * r / A - 1j * sig + q2) ** al)
for al, lmb in ((2.5, 1.7), (20.0, 0.6)):
    r, tp, sp = 1.3, 0.4, -2.2
    lhs = f_ts(r, lmb * tp, sp / lmb, 1.0, 20.0, al)
    rhs = lmb ** (al - 1) * f_ts(r, tp, sp, 1.0 / lmb, 20.0 * lmb, al)
    check(f"boost maps (q1, q2) -> (q1/lambda, lambda q2) (alpha {al})", abs(lhs - rhs) < 1e-12 * abs(lhs))

# On-axis H_z = 4 alpha/(A^2 Q^(alpha+1)): at fixed tau it falls as (2t)^-(alpha+1).
for al in (1.0, 5.0):
    tau0 = 0.3
    v1 = abs(fields(0.0, tau0 + 1e4, 1e4, 1.0, 20.0, al)[2])
    v2 = abs(fields(0.0, tau0 + 2e4, 2e4, 1.0, 20.0, al)[2])
    check(f"on-axis H_z decay exponent alpha + 1 = {al+1:g}", abs(math.log(v1 / v2, 2) - (al + 1)) < 1e-3)

# Whole-plane polarity p = (cos b(inf) - cos b(0))/2 of Re H in transverse planes at t = 0.
def polarity(al, z, q1=1.0, q2=20.0):
    r = np.logspace(-6, 8, 4001)
    A, Q = q1 + 1j * z, q2 - 1j * z
    D = r * r + A * Q
    lA, lD = np.log(A), np.log(D)
    off = -(al + 1) * lD.real
    t1 = np.exp((al - 1) * lA - (al + 2) * lD - off)
    t2 = np.exp((al - 2) * lA - (al + 1) * lD - off)
    hr = (2 * al * 1j * r * ((al + 1) * t1 * (q2 - q1 - 2j * z) - (al - 1) * t2)).real
    hz = (4 * al * t1 * (A * Q - al * r * r)).real
    cb = hz / np.hypot(hr, hz)
    return 0.5 * (cb[-1] - cb[0])
expect = {(1.0, 0.5): -1.0, (1.0, 2.0): 0.0, (2.0, 3.0): 0.0, (5.0, 1.0): -0.5, (5.0, 3.0): 0.5,
          (20.0, 14.0): 0.5, (20.0, 10.0): -0.5, (20.0, 6.0): 0.5, (20.0, 3.5): -0.5}
for (al, z), pv in expect.items():
    got = polarity(al, z)
    check(f"whole-plane polarity at alpha {al:g}, z = {z:g} q1 is {pv:+g}", abs(got - pv) < 1e-5, f"{got:+.6f}")

# On-axis zeros of Re H_z at t = 0 for z in [0, 80 q1], q2 = 20 q1 (Shen et al. Fig. 6a): our H_z drops
# from 4 to 2 zeros between alpha = 5.5 and 5.6, the feature they report; eq. (6) as printed does not.
zax = np.linspace(1e-6, 80, 800001)
Aax, Qax = 1 + 1j * zax, 20 - 1j * zax
def nzeros(v):
    s = np.sign(v)
    return int(np.sum(s[1:] * s[:-1] < 0))
ours = {al: nzeros((4 * al / (Aax**2 * Qax ** (al + 1))).real) for al in (5.5, 5.6)}
shen = {al: nzeros((-4 * al * (-al / (Aax**2 * Qax ** (al + 1)) + (al - 1) / (Aax**3 * Qax**al))).real) for al in (5.5, 5.6)}
check("on-axis zeros: ours 4 -> 2 at alpha 5.5 -> 5.6 (Shen Fig. 6 drop); printed eq. (6) gives 3 -> 3",
      ours == {5.5: 4, 5.6: 2} and shen == {5.5: 3, 5.6: 3}, f"{ours} {shen}")

print()
if FAILS:
    print(f"{len(FAILS)} FAILED: " + "; ".join(FAILS))
    sys.exit(1)
print("all checks passed")
