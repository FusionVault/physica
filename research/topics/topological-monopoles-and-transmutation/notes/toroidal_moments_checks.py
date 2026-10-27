#!/usr/bin/env python3
"""Checks for toroidal_moments.md. Standard library only; run with `python3 -I`.

Every computed number is asserted against a published value (source in the
comment) with a stated tolerance. Exit code is non-zero on any failure.
Conventions: SI units, time dependence exp(-i*omega*t), toroidal dipole
T = (1/10) * Integral[(r.J) r - 2 r^2 J] dV (SI; Gaussian texts divide by c).
"""
import cmath
import math
import sys

# CODATA 2018 exact / recommended constants
H = 6.62607015e-34
E = 1.602176634e-19
HBAR = H / (2 * math.pi)
C = 299792458.0
MU0 = 1.25663706212e-6
EPS0 = 1.0 / (MU0 * C * C)

failures = []


def check(name, got, want, rel=None, abs_=None):
    err = abs(got - want)
    tol = abs_ if abs_ is not None else rel * abs(want)
    ok = err <= tol
    print(f"{'PASS' if ok else 'FAIL'}  {name}: got {got:.6g}, want {want:.6g} (tol {tol:.3g})")
    if not ok:
        failures.append(name)


# --- vector helpers ---------------------------------------------------------
def add(a, b): return (a[0] + b[0], a[1] + b[1], a[2] + b[2])
def sub(a, b): return (a[0] - b[0], a[1] - b[1], a[2] - b[2])
def sc(a, s): return (a[0] * s, a[1] * s, a[2] * s)
def dot(a, b): return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
def cross(a, b): return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def toroidal_solenoid(n_turns, current, R, a, m_seg):
    """Segments (midpoint, I*dl) of n_turns circular turns (minor radius a) on a torus (major radius R)."""
    segs = []
    for j in range(n_turns):
        ph = 2 * math.pi * (j + 0.5) / n_turns
        u = (math.cos(ph), math.sin(ph), 0.0)
        pts = [add(sc(u, R + a * math.cos(2 * math.pi * k / m_seg)), (0.0, 0.0, a * math.sin(2 * math.pi * k / m_seg)))
               for k in range(m_seg + 1)]
        for k in range(m_seg):
            mid = sc(add(pts[k], pts[k + 1]), 0.5)
            segs.append((mid, sc(sub(pts[k + 1], pts[k]), current)))
    return segs


def toroidal_dipole(segs):
    """T = (1/10) sum I[(r.dl) r - 2 r^2 dl]  (SI, A m^3)."""
    t = (0.0, 0.0, 0.0)
    for r, idl in segs:
        t = add(t, sub(sc(r, dot(r, idl)), sc(idl, 2 * dot(r, r))))
    return sc(t, 0.1)


def magnetic_dipole(segs):
    m = (0.0, 0.0, 0.0)
    for r, idl in segs:
        m = add(m, cross(r, idl))
    return sc(m, 0.5)


# 1. Toroidal dipole of a toroidal solenoid: |T| = N I V / (4 pi) = N I pi R a^2 / 2
#    Nemkov, Basharin & Fedotov, PRA 98, 023858 (2018), App. B: magnitude c^-1 pi^2 j0 R rho = I V/(4 pi c)
#    (Gaussian; SI drops the 1/c).
N_T, I_T, R_T = 24, 1.0, 1.0
for a_over_R in (0.05, 0.5):
    segs = toroidal_solenoid(N_T, I_T, R_T, a_over_R * R_T, 200)
    T = toroidal_dipole(segs)
    want = N_T * I_T * math.pi * R_T * (a_over_R * R_T) ** 2 / 2
    check(f"toroid |T| (a/R={a_over_R}) vs N I V/(4 pi)", abs(T[2]), want, rel=1e-3)
    check(f"toroid net magnetic dipole (a/R={a_over_R}) is zero", math.sqrt(dot(magnetic_dipole(segs), magnetic_dipole(segs))), 0.0, abs_=1e-9)

# 1b. Same T from (1/2) sum r_c x m over turns (the identity T = (1/2) Integral r x M dV for J = curl M)
a_T = 0.3
segs = toroidal_solenoid(N_T, I_T, R_T, a_T, 200)
T = toroidal_dipole(segs)
half_rxm = (0.0, 0.0, 0.0)
for j in range(N_T):
    ph = 2 * math.pi * (j + 0.5) / N_T
    centre = (R_T * math.cos(ph), R_T * math.sin(ph), 0.0)
    m_turn = magnetic_dipole([(sub(r, centre), idl) for r, idl in segs[j * 200:(j + 1) * 200]])
    half_rxm = add(half_rxm, sc(cross(centre, m_turn), 0.5))
check("T(1/10 formula) == (1/2) sum r x m", T[2], half_rxm[2], rel=1e-3)


# 2. Radiation: far-field integral F(n) = sum I dl exp(-i k n.r); transverse part -> -i w (p + (ik/c) T).
#    Pure toroid (p = 0): F_perp = k^2 T_perp, power P_T = mu0 w^6 |T|^2 / (12 pi c^5)  (w^6 vs w^4 for p, m).
#    Published: the toroidal far field repeats the electric-dipole pattern scaled by a factor w^2
#    (Miroshnichenko et al., Nat. Commun. 6, 8069 (2015)); anapole d + i k tau = 0 (Nemkov et al. 2018, Gaussian).
def far_field(segs, k, theta):
    n = (math.sin(theta), 0.0, math.cos(theta))
    fx = fy = fz = 0j
    for r, idl in segs:
        ph = cmath.exp(-1j * k * dot(n, r))
        fx += idl[0] * ph; fy += idl[1] * ph; fz += idl[2] * ph
    return n, (fx, fy, fz)


def perp_sq(n, f):  # |n x F|^2 for complex F
    cx = (n[1] * f[2] - n[2] * f[1], n[2] * f[0] - n[0] * f[2], n[0] * f[1] - n[1] * f[0])
    return sum(abs(z) ** 2 for z in cx)


def radiated_power(segs, k, n_theta=90, extra_F=None):
    w = k * C
    tot = 0.0
    for i in range(n_theta):
        th = math.pi * (i + 0.5) / n_theta
        n, f = far_field(segs, k, th)
        if extra_F is not None:
            f = tuple(f[q] + extra_F[q] for q in range(3))
        tot += perp_sq(n, f) * math.sin(th) * (math.pi / n_theta) * 2 * math.pi  # azimuthal symmetry
    return MU0 * C * k * k / (32 * math.pi ** 2) * tot, w


segs = toroidal_solenoid(N_T, I_T, R_T, a_T, 60)
T = toroidal_dipole(segs)
k1 = 0.05 / R_T
n, f = far_field(segs, k1, math.pi / 2)
check("F_perp(theta=90deg) = k^2 T_z", f[2].real, k1 * k1 * T[2], rel=1e-2)
P1, w1 = radiated_power(segs, k1)
check("toroid power vs mu0 w^6 T^2/(12 pi c^5)", P1, MU0 * w1 ** 6 * T[2] ** 2 / (12 * math.pi * C ** 5), rel=1e-2)
P2, _ = radiated_power(segs, 2 * k1)
check("toroid power scales as w^6 (P(2w)/P(w) = 64)", P2 / P1, 64.0, rel=2e-2)
# Anapole: add a point electric dipole p = -(i k / c) T; its F = -i w p.  Residual power must vanish.
p_z = -1j * k1 / C * T[2]
F_p = (0j, 0j, -1j * k1 * C * p_z)
P_ana, _ = radiated_power(segs, k1, extra_F=F_p)
check("anapole p = -(ik/c)T: residual/uncancelled power < 1e-3", P_ana / P1, 0.0, abs_=1e-3)
# Wrong sign doubles the amplitude -> 4x power (shows the sign matters with exp(-i w t)).
P_bad, _ = radiated_power(segs, k1, extra_F=tuple(-x for x in F_p))
check("p = +(ik/c)T gives 4x power", P_bad / P1, 4.0, rel=1e-2)


# 3. Mie electric-dipole zero (a1 = 0) of a sphere, n = 4, lambda = 550 nm: "diameter ~204 nm"
#    (Miroshnichenko et al. 2015, Fig. 2: P_sph = 0 at ~204 nm).
def psi1(z): return math.sin(z) / z - math.cos(z)
def dpsi1(z): return math.cos(z) / z - math.sin(z) / z ** 2 + math.sin(z)


def a1_numerator(d_nm, m=4.0, lam=550.0):
    x = math.pi * d_nm / lam
    return m * psi1(m * x) * dpsi1(x) - psi1(x) * dpsi1(m * x)


lo, hi = 180.0, 230.0
assert a1_numerator(lo) * a1_numerator(hi) < 0, "no a1 zero bracketed"
for _ in range(60):
    mid = 0.5 * (lo + hi)
    if a1_numerator(lo) * a1_numerator(mid) <= 0: hi = mid
    else: lo = mid
check("Mie a1 = 0 diameter, n=4, 550 nm [nm]", 0.5 * (lo + hi), 204.0, abs_=3.0)


# 4. Aharonov-Bohm / Tonomura: flux quantum h/2e and AB phase e Phi / hbar = pi for one h/2e.
#    Osakabe et al. PRA 34, 815 (1986): relative phase "0 or pi due to ... flux quantization in units of h/2e".
PHI0 = H / (2 * E)
check("flux quantum h/2e [Wb] (CODATA 2.067833848e-15)", PHI0, 2.067833848e-15, rel=1e-9)
check("AB phase for one h/2e = pi", E * PHI0 / HBAR, math.pi, rel=1e-12)
check("AB phase for two h/2e = 2 pi (indistinguishable from 0)", (E * 2 * PHI0 / HBAR) % (2 * math.pi), 0.0, abs_=1e-9)
# Dirac monopole: e g = 2 pi hbar n with g the total flux -> g_D = h/e = 2 Phi0 (webers); as a
# pole strength in A m (B = mu0 q_m / (4 pi r^2)) it is h/(mu0 e).
check("Dirac flux h/e = 2 x (h/2e)", H / E, 2 * PHI0, rel=1e-12)

# 5. Puthoff US 5,845,220 quotes a Josephson frequency of 483.2 MHz/uV; 2e/h (CODATA) = 483.5978 MHz/uV.
kj = 2 * E / H / 1e6 / 1e6  # MHz per microvolt
check("Josephson 2e/h [MHz/uV] vs patent's 483.2 (0.08% low)", kj, 483.2, rel=1e-3)

# 6. Lockheed Martin US 9,502,202: mean free path ~9.3 mm, collision time ~100 ns, speed "about 100 km/s".
check("US 9,502,202 speed = 9.3 mm / 100 ns [m/s]", 9.3e-3 / 100e-9, 1.0e5, rel=0.1)

# 7. 133Cs anapole: kappa = 0.127 +- 0.019, "differing from zero by ~7 sigma" (Haxton, Science 275, 1753 (1997)).
check("Cs anapole kappa / sigma", 0.127 / 0.019, 7.0, abs_=0.5)

# 8. MEG scale (Hamalainen et al., RMP 65, 413 (1993)).
#    "a dipole moment of 10 nAm would correspond to 40 mm^2 of active cortex" at 250 nA/mm^2 x 1 mm.
check("cortex area for 10 nAm at 250 nA/mm^2, 1 mm [mm^2]", (10.0 * 1000.0) / (250.0 * 1.0), 40.0, rel=1e-12)  # nA*mm / (nA/mm)
#    Radial field of a 10 nAm tangential dipole 30 mm below the surface of a 100 mm sphere. The radial
#    component is unaffected by volume currents in a spherical conductor, so Biot-Savart of Q is exact.
Q = (10e-9, 0.0, 0.0); r0 = (0.0, 0.0, 0.070)
bmax = 0.0
for i in range(1, 400):
    th = 0.6 * i / 400
    for sgn in (1, -1):  # field extrema lie in the plane perpendicular to Q
        rr = (0.0, sgn * 0.100 * math.sin(th), 0.100 * math.cos(th)); d = sub(rr, r0)
        br = MU0 / (4 * math.pi) * dot(cross(Q, d), sc(rr, 1 / 0.100)) / dot(d, d) ** 1.5
        bmax = max(bmax, abs(br))
print(f"info  MEG peak radial field of 10 nAm at 30 mm depth: {bmax*1e15:.0f} fT "
      "(published: signals 'typically 50-500 fT')")
check("MEG peak field inside the published 50-500 fT band [fT] (band check)", min(max(bmax * 1e15, 50), 500), bmax * 1e15, rel=1e-12)

# 8b. Illustration (no published value): cortical currents arranged as a toroid. Loops of total
#     current-dipole strength 10 nAm closing on radius 5 mm (m = Q rho / 2), ring radius 20 mm (T = R m / 2).
m_c = 10e-9 * 5e-3 / 2
T_c = 0.020 * m_c / 2
w_a = 2 * math.pi * 10.0
print(f"info  cortical toroid: m = {m_c:.2g} A m^2, T = {T_c:.2g} A m^3, static field outside = 0, "
      f"radiated at 10 Hz = {MU0 * w_a**6 * T_c**2 / (12 * math.pi * C**5):.2g} W")

# 9. Thin-wire heating. IACS: 100% conductivity = 58 MS/m -> 1.7241e-8 ohm m at 20 C.
RHO_CU = 1 / 58e6
check("IACS resistivity [ohm m]", RHO_CU, 1.7241e-8, rel=1e-4)
d_w = 8e-6
r_per_m = RHO_CU / (math.pi * (d_w / 2) ** 2)
print(f"info  8 um copper: {r_per_m:.1f} ohm/m")
for v in (230.0, 10e3, 100e3):
    i_w = 1000.0 / v
    print(f"info  1 kW at {v/1e3:g} kV through 8 um Cu: I = {i_w:.3g} A, J = {i_w/(math.pi*(d_w/2)**2):.3g} A/m^2, "
          f"I^2R = {i_w**2*r_per_m:.3g} W per metre")

if failures:
    print(f"\n{len(failures)} FAILED: {failures}")
    sys.exit(1)
print("\nall checks passed")
