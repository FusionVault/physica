#!/usr/bin/env python3
"""Checks for monopole_toroidal.md. Standard library only; run with `python3 -I`.

Every derived number in the note is recomputed here and asserted, against a
published value where one exists (source in the comment) or against an
independent route (closed form vs. brute-force sum or quadrature). Exits 1 on
any failure.

Conventions: SI, time dependence exp(-i omega t), magnetic charge in A m
(div B = mu0 rho_m, -curl E = dB/dt + mu0 J_m), the crate's toroidal dipole
T = (1/10) Int[(r.J) r - 2 r^2 J] dV and anapole p = -(ik/c) T.
Constants are the ones in crates/metrology (CODATA 2018), so numbers match
the magnetics crate's tests.
"""
import cmath
import math
import random
import sys

H = 6.62607015e-34
E = 1.602176634e-19
HBAR = H / (2 * math.pi)
C = 299792458.0
MU0 = 1.25663706212e-6
EPS0 = 1.0 / (MU0 * C * C)
MU_B = 9.2740100783e-24
ALPHA = 7.2973525693e-3
M_E = 9.1093837015e-31

G_D = H / (MU0 * E)          # Dirac charge, A m (magnetics::monopole::dirac_charge(1))
DIRAC_FLUX = H / E           # mu0 g_D = h/e
FLUX_QUANTUM = H / (2 * E)   # h/2e

failures = []
count = 0


def check(name, got, want, rel=None, abs_=None):
    global count
    count += 1
    err = abs(got - want)
    tol = abs_ if abs_ is not None else rel * abs(want)
    ok = err <= tol
    print(f"{'PASS' if ok else 'FAIL'}  {name}: got {got:.6g}, want {want:.6g} (tol {tol:.3g})")
    if not ok:
        failures.append(name)


def add(a, b): return tuple(x + y for x, y in zip(a, b))
def sub(a, b): return tuple(x - y for x, y in zip(a, b))
def sc(a, s): return tuple(x * s for x in a)
def dot(a, b): return sum(x * y for x, y in zip(a, b))
def cross(a, b): return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])
def norm(a): return math.sqrt(dot(a, a))
def cabs2(a): return sum(abs(x) ** 2 for x in a)


# =============================================================================
# 1. Dyonic multipoles: Ampere and Gilbert pictures of one oscillating magnetised
#    torus radiate identically, which fixes the electric dipole of a magnetic
#    current loop at  p = -d_m / c^2,  d_m = (1/2) Int r x J_m dV.
# =============================================================================
print("\n# 1. Dyonic doublet, Ampere/Gilbert radiation, duality invariance")

N_TURNS, I0, R0, A0 = 24, 1.0, 1.0, 0.3
K = 0.05
OMEGA = K * C


def toroidal_solenoid(n_turns, current, big_r, a, m_seg):
    segs = []
    for j in range(n_turns):
        ph = 2 * math.pi * (j + 0.5) / n_turns
        u = (math.cos(ph), math.sin(ph), 0.0)
        pts = [add(sc(u, big_r + a * math.cos(2 * math.pi * k / m_seg)), (0.0, 0.0, a * math.sin(2 * math.pi * k / m_seg)))
               for k in range(m_seg + 1)]
        for k in range(m_seg):
            segs.append((sc(add(pts[k], pts[k + 1]), 0.5), sc(sub(pts[k + 1], pts[k]), current)))
    return segs


def toroidal_dipole(segs):
    t = (0.0, 0.0, 0.0)
    for r, idl in segs:
        t = add(t, sub(sc(r, dot(r, idl)), sc(idl, 2 * dot(r, r))))
    return sc(t, 0.1)


# This winding sense with +I gives T along -z (magnetisation along -phi_hat), so the
# Ampere picture of a core magnetised along +phi_hat carries -I.
segs = toroidal_solenoid(N_TURNS, -I0, R0, A0, 40)
T_A = toroidal_dipole(segs)[2]
T_IDEAL = N_TURNS * I0 * math.pi * R0 * A0 ** 2 / 2       # magnetics::toroidal::toroidal_solenoid_moment
check("Ampere toroid T (40-gon turns) vs N I V/(4 pi)", T_A, T_IDEAL, rel=5e-3)

# Gilbert picture of the same magnet: core magnetisation M = N I/(2 pi rho) phi_hat,
# so rho_m = -div M = 0 and the only source is J_m = dM/dt = -i omega M (azimuthal ring).
GL4 = [(-0.8611363115940526, 0.3478548451374538), (-0.3399810435848563, 0.6521451548625461),
       (0.3399810435848563, 0.6521451548625461), (0.8611363115940526, 0.3478548451374538)]
core = []          # (r, M dV)
NTH, NPH = 16, 48
for xs, ws in GL4:
    s = 0.5 * A0 * (xs + 1)
    for it in range(NTH):
        th = 2 * math.pi * (it + 0.5) / NTH
        for ip in range(NPH):
            ph = 2 * math.pi * (ip + 0.5) / NPH
            rho = R0 + s * math.cos(th)
            r = (rho * math.cos(ph), rho * math.sin(ph), s * math.sin(th))
            # M dV = [N I/(2 pi rho)] * s rho ds dth dph: the rho cancels.
            mdv = N_TURNS * I0 / (2 * math.pi) * s * (0.5 * A0 * ws) * (2 * math.pi / NTH) * (2 * math.pi / NPH)
            core.append((r, sc((-math.sin(ph), math.cos(ph), 0.0), mdv)))
t_core = 0.5 * sum(cross(r, m)[2] for r, m in core)
check("Gilbert core t = (1/2) Int r x M dV = N I V/(4 pi)", t_core, T_IDEAL, rel=1e-9)
mag_total = [sum(m[q] for _, m in core) for q in range(3)]
check("Gilbert core net moment Int M dV = 0", math.sqrt(sum(x * x for x in mag_total)), 0.0, abs_=1e-9)

# d_m = (1/2) Int r x J_m = -i omega t, so -d_m/c^2 = (i omega/c^2) t = (ik/c) t: the same
# electric dipole the crate's p_eff = p + (ik/c)T gives the Ampere toroid.
d_m = -1j * OMEGA * t_core
check("-d_m/c^2 equals (ik/c) t (imaginary part)", (-d_m / C ** 2).imag, (K / C * t_core), rel=1e-12)

# Far field: E ~ (i mu0 omega/4 pi r) e^{ikr} G,  G = (n x F_e) x n - (1/c) n x F_m,
# F = Sum J dV e^{-ik n.r}.  P = mu0 c k^2/(32 pi^2) Oint |G|^2 dOmega.  The J_m term is the
# xi = -pi/2 duality image of the electric-current far field (E' = -cB0, J_e0 = J_m/c).
NT, NP = 40, 16
dirs = []
for it in range(NT):
    th = math.pi * (it + 0.5) / NT
    for ip in range(NP):
        ph = 2 * math.pi * (ip + 0.5) / NP
        dirs.append(((math.sin(th) * math.cos(ph), math.sin(th) * math.sin(ph), math.cos(th)),
                     math.sin(th) * (math.pi / NT) * (2 * math.pi / NP)))


def ff(elems, n):
    f = [0j, 0j, 0j]
    for r, jdv in elems:
        ph = cmath.exp(-1j * K * dot(n, r))
        for q in range(3):
            f[q] += jdv[q] * ph
    return f


FA = [ff(segs, n) for n, _ in dirs]                                         # electric current, real I dl
FG = [[-1j * OMEGA * x for x in ff(core, n)] for n, _ in dirs]              # magnetic current J_m dV


def power(ce=0.0, cm=0.0, p_z=0j, cm_from_e=None):
    """Power of ce*(Ampere toroid) + cm*(Gilbert core, J_m) + electric dipole p_z at the origin."""
    tot = 0.0
    for (n, w), fa, fg in zip(dirs, FA, FG):
        fe = [ce * x for x in fa]
        fe[2] += -1j * OMEGA * p_z
        fm = [cm * x for x in fg]
        t1 = cross(cross(n, fe), n)
        t2 = cross(n, fm)
        g = [t1[q] - t2[q] / C for q in range(3)]
        tot += cabs2(g) * w
    return MU0 * C * K * K / (32 * math.pi ** 2) * tot


P_dip = MU0 * OMEGA ** 6 * T_IDEAL ** 2 / (12 * math.pi * C ** 5)            # toroidal_radiated_power
P_A, P_G = power(ce=1.0), power(cm=1.0)
check("Ampere toroid power vs mu0 w^6 T^2/(12 pi c^5)", P_A, P_dip, rel=2e-2)
check("Gilbert (magnetic-current) core power vs same", P_G, P_dip, rel=2e-2)
check("Ampere toroid + Gilbert core of -M: radiation cancels (frac. of P_A)", power(ce=1.0, cm=-1.0) / P_A, 0.0, abs_=1e-3)
check("Ampere + Gilbert of +M doubles the field: 4x power", power(ce=1.0, cm=1.0) / P_A, 4.0, rel=2e-2)
p_ana = -1j * K / C * T_IDEAL
check("Gilbert core + crate anapole p = -(ik/c)T: power fraction", power(cm=1.0, p_z=p_ana) / P_G, 0.0, abs_=1e-3)
check("Gilbert core + opposite p: 4x power", power(cm=1.0, p_z=-p_ana) / P_G, 4.0, rel=2e-2)

# Duality rotation (c J_e, J_m) -> rotated by xi leaves the power unchanged direction by direction:
# G' = cos(xi) G + sin(xi) n x G, and Re[G* . (n x G)] = 0.
xi = 0.7
tot0 = tot1 = 0.0
for (n, w), fa, fg in zip(dirs, FA, FG):
    fe = [x + 0.3 * y / C for x, y in zip(fa, fg)]
    fe[2] += -1j * OMEGA * (0.4 * p_ana)
    fm = [0.5 * y - 0.2 * C * x for x, y in zip(fa, fg)]
    fe2 = [math.cos(xi) * a + math.sin(xi) * b / C for a, b in zip(fe, fm)]
    fm2 = [-math.sin(xi) * C * a + math.cos(xi) * b for a, b in zip(fe, fm)]
    for (e_, m_), acc in (((fe, fm), 0), ((fe2, fm2), 1)):
        t1 = cross(cross(n, e_), n)
        t2 = cross(n, m_)
        v = cabs2([t1[q] - t2[q] / C for q in range(3)]) * w
        if acc == 0:
            tot0 += v
        else:
            tot1 += v
check("duality rotation xi = 0.7 leaves radiated power invariant", tot1 / tot0, 1.0, rel=1e-12)

# =============================================================================
# 2. Dirac strings, toroidal solenoids, composite and fractional charges
# =============================================================================
print("\n# 2. Dirac string = closed thin toroid of flux h/e; AB phases")
check("g_D = h/(mu0 e) (crate dirac_charge(1), 3.291e-9 A m)", G_D, 3.291e-9, rel=2e-4)
check("mu0 g_D = h/e = 2 Phi_0", MU0 * G_D / FLUX_QUANTUM, 2.0, rel=1e-15)


def ab(q, flux):
    return q * flux / HBAR


check("electron round h/e: 2 pi (invisible)", ab(E, DIRAC_FLUX), 2 * math.pi, rel=1e-15)
check("Cooper pair (2e) round h/e: 4 pi", ab(2 * E, DIRAC_FLUX), 4 * math.pi, rel=1e-15)
check("electron round h/2e: pi (Osakabe 1986 '0 or pi')", ab(E, FLUX_QUANTUM), math.pi, rel=1e-15)
check("charge e/3 round h/e: 2 pi/3 mod 2 pi (visible)", math.fmod(ab(E / 3, DIRAC_FLUX), 2 * math.pi), 2 * math.pi / 3, rel=1e-12)
check("charge 2e/3 round h/e: 4 pi/3 mod 2 pi (visible)", math.fmod(ab(2 * E / 3, DIRAC_FLUX), 2 * math.pi), 4 * math.pi / 3, rel=1e-12)
check("e/3 round 3 h/e (pole 3 g_D): 2 pi, invisible again", ab(E / 3, 3 * DIRAC_FLUX), 2 * math.pi, rel=1e-12)
# A neutron (u d d) whose three quarks all encircle the line: total charge 0, phase 0.
check("neutron udd all encircling: sum of quark AB phases", sum(ab(q * E, DIRAC_FLUX) for q in (2 / 3, -1 / 3, -1 / 3)), 0.0, abs_=1e-12)

# Dual AB (Dowling, Williams & Franson 1999): a pole g circling electric flux Psi = Int D.dA
# (coulombs) gains mu0 g Psi/hbar in magnitude; for g = g_D that is 2 pi Psi/e.
check("dual AB: g_D round Psi = e gives 2 pi", MU0 * G_D * E / HBAR, 2 * math.pi, rel=1e-15)
# Polarisation quantum e R/Omega (R normal to a cell face, Omega = |R| A_face): the flux
# through one face changes by e, the dual AB phase by exactly 2 pi.
a_cell = 4.0e-10
dP = E * a_cell / a_cell ** 3
check("polarisation quantum x face area = e", dP * a_cell ** 2 / E, 1.0, rel=1e-15)
check("... shifts g_D's dual AB phase by 2 pi", MU0 * G_D * dP * a_cell ** 2 / HBAR, 2 * math.pi, rel=1e-12)

# Thin toroidal solenoid: T = Phi * (pi R^2)/mu0 as a/R -> 0 (exact flux of a circular-section
# ideal toroid: Phi = mu0 N I (R - sqrt(R^2 - a^2))).
for a_over_r, tol in ((0.01, 3e-5), (0.1, 3e-3)):
    a = a_over_r * R0
    phi = MU0 * N_TURNS * I0 * (R0 - math.sqrt(R0 ** 2 - a ** 2))
    t_crate = N_TURNS * I0 * math.pi * R0 * a ** 2 / 2
    check(f"T_crate / (Phi pi R^2/mu0) at a/R = {a_over_r}", t_crate / (phi * math.pi * R0 ** 2 / MU0), 1.0, rel=tol)
# So a closed Dirac string of n units round a loop of vector area A has "T" = n g_D A.
area_1nm = math.pi * 1e-9 ** 2
check("closed Dirac string loop, r = 1 nm: g_D pi r^2 (A m^3)", G_D * area_1nm, 1.0339e-26, rel=1e-4)

# Charge-transport identity: a magnetic charge Q carried round any closed polygon leaves a
# magnetisation line Q dl whose toroidal moment (1/2) Int r x M dV equals Q times the loop's
# vector area, whatever the origin.
random.seed(7)
poly = [(random.uniform(-1, 1), random.uniform(-1, 1), random.uniform(-1, 1)) for _ in range(9)]
Q_test = 2.7
for shift in ((0.0, 0.0, 0.0), (5.0, -3.0, 2.0)):
    pts = [add(p, shift) for p in poly]
    tm = (0.0, 0.0, 0.0)
    va = (0.0, 0.0, 0.0)
    for i in range(len(pts)):
        p, q = pts[i], pts[(i + 1) % len(pts)]
        tm = add(tm, sc(cross(sc(add(p, q), 0.5), sc(sub(q, p), Q_test)), 0.5))
        va = add(va, sc(cross(p, q), 0.5))
    for c_ in range(3):
        check(f"Delta t = Q x vector area, component {c_}, origin shift {shift[0]}", tm[c_], Q_test * va[c_], rel=1e-12, abs_=None)

# =============================================================================
# 3. Spin ice (Dy2Ti2O7): string flux, hexagon-loop toroidal moments
# =============================================================================
print("\n# 3. Spin ice")
MU = 10 * MU_B
R_NN = 3.54e-10                       # the crate's `a`: nearest-neighbour Dy-Dy distance
A_D = math.sqrt(1.5) * R_NN           # diamond bond, tetrahedron centre to centre
Q = 2 * MU / A_D                      # magnetics::monopole::spin_ice_charge(10 mu_B, 3.54 A)
check("spin_ice_charge(10 mu_B, 3.54 A) = 4.278e-13 A m (crate test)", Q, 4.278e-13, rel=1e-3)
A_CUBIC = 2 * math.sqrt(2) * R_NN
check("cubic constant 2 sqrt2 r_nn vs Dhar et al. 'a = 10.0(1) A'", A_CUBIC * 1e10, 10.0, abs_=0.1)
check("g_D/Q = 7692 (crate test)", G_D / Q, 7692.0, rel=1e-3)
phi_string = 2 * math.pi * Q / G_D    # e mu0 Q / hbar
check("AB step across a spin-ice Dirac string 2 pi Q/g_D vs Dhar et al. '~0.8 mrad'", phi_string * 1e3, 0.8, rel=5e-2)
check("... = 40 pi r_e/a_d (Castelnovo eq. 4 form)", phi_string, 40 * math.pi * (ALPHA * HBAR / (M_E * C)) / A_D, rel=1e-6)
check("Dhar et al. fitted q_m = (4.4 +- 0.3)e-13 A m is within 1 sigma", abs(Q - 4.4e-13) / 0.3e-13, 0.0, abs_=1.0)

# Diamond lattice of tetrahedron centres; find a hexagonal ring through the origin.
fcc = [(0, 0, 0), (0, 0.5, 0.5), (0.5, 0, 0.5), (0.5, 0.5, 0)]
sites = []
for i in range(-2, 3):
    for j in range(-2, 3):
        for k in range(-2, 3):
            for b in fcc:
                for basis in ((0, 0, 0), (0.25, 0.25, 0.25)):
                    sites.append(tuple(A_CUBIC * (x + y + z) for x, y, z in zip((i, j, k), b, basis)))
bond = A_CUBIC * math.sqrt(3) / 4
check("diamond bond sqrt3 a_c/4 = sqrt(3/2) r_nn", bond, A_D, rel=1e-12)
origin = min(sites, key=lambda s: dot(s, s))
nbrs = {}
near = [s for s in sites if dot(sub(s, origin), sub(s, origin)) < (3 * bond) ** 2]
for s in near:
    nbrs[s] = [t for t in near if abs(norm(sub(s, t)) - bond) < 1e-3 * bond]


def find_ring():
    def dfs(path):
        if len(path) == 6:
            return path if origin in nbrs[path[-1]] else None
        for t in nbrs[path[-1]]:
            if t not in path:
                r = dfs(path + [t])
                if r:
                    return r
        return None
    return dfs([origin])


ring = find_ring()
va = (0.0, 0.0, 0.0)
for i in range(6):
    va = add(va, sc(cross(ring[i], ring[(i + 1) % 6]), 0.5))
area = norm(va)
check("vector area of the diamond chair hexagon = (4/sqrt3) a_d^2", area, 4 / math.sqrt(3) * A_D ** 2, rel=1e-9)
axis = sc(va, 1 / area)
check("... along a <111> axis (|component| = 1/sqrt3)", abs(axis[0]), 1 / math.sqrt(3), rel=1e-9)
# Spins on the bond midpoints, head to tail along the ring (a flippable loop).
spins = []
for i in range(6):
    p, q = ring[i], ring[(i + 1) % 6]
    spins.append((sc(add(p, q), 0.5), sc(sub(q, p), MU / A_D)))
cen = sc(tuple(sum(s[0][q] for s in spins) for q in range(3)), 1 / 6)
radii = [norm(sub(s[0], cen)) for s in spins]
check("spin sites form a hexagon of circumradius r_nn", max(radii), R_NN, rel=1e-9)
check("... and are coplanar (max offset from the plane / r_nn)", max(abs(dot(sub(s[0], cen), axis)) for s in spins) / R_NN, 0.0, abs_=1e-9)
t_loop = dot(sc(tuple(sum(cross(sub(r, cen), m)[q] for r, m in spins) for q in range(3)), 0.5), axis)
check("toroidal moment of a flippable hexagon = (Q/2) |A| = 2 sqrt2 mu r_nn", abs(t_loop), 2 * math.sqrt(2) * MU * R_NN, rel=1e-9)
check("... in mu_B A (100.1)", abs(t_loop) / (MU_B * 1e-10), 100.13, rel=1e-3)
check("loop flip changes t by Q |A| (200.3 mu_B A)", 2 * abs(t_loop), Q * area, rel=1e-9)
# The same change from the monopole current: J_m dt = Q dl along the six bonds.
t_jm = 0.5 * sum(dot(cross(sub(sc(add(ring[i], ring[(i + 1) % 6]), 0.5), cen), sc(sub(ring[(i + 1) % 6], ring[i]), Q)), axis) for i in range(6))
check("(1/2) Int r x J_m dt over the monopole's loop = flip change", abs(t_jm), 2 * abs(t_loop), rel=1e-9)

# Electron phase along the hexagon's axis, from the point dipoles' A = mu0 m x R/(4 pi R^3):
# Int A.axis dl over the infinite axis, by quadrature, against the closed form
# (mu0/2 pi) Sum (rho_hat x m).axis/rho per dipole.
def a_axis_integral(spins_):
    tot = 0.0
    nq = 4000
    for iq in range(nq):
        u = -math.pi / 2 + math.pi * (iq + 0.5) / nq      # l = r_nn tan u
        l = R_NN * math.tan(u)
        dl = R_NN / math.cos(u) ** 2 * (math.pi / nq)
        x = add(cen, sc(axis, l))
        a_vec = (0.0, 0.0, 0.0)
        for r, m in spins_:
            d = sub(x, r)
            a_vec = add(a_vec, sc(cross(m, d), MU0 / (4 * math.pi) / norm(d) ** 3))
        tot += dot(a_vec, axis) * dl
    return tot


ia = a_axis_integral(spins)
closed = sum(MU0 / (2 * math.pi) * dot(cross(sc(sub(r, cen), 1 / R_NN), m), axis) / R_NN for r, m in spins)
check("axis integral of A: quadrature vs closed form", ia, closed, rel=1e-6)
check("closed form = (mu0/2pi) 4 sqrt2 mu/r_nn", abs(closed), MU0 / (2 * math.pi) * 4 * math.sqrt(2) * MU / R_NN, rel=1e-9)
phase_axis = E * abs(ia) / HBAR
phase_cont = math.pi * Q / G_D                        # continuum line of moment Q/2 per length
check("continuum loop phase pi Q/g_D (mrad)", phase_cont * 1e3, 0.4084, rel=1e-3)
check("discrete on-axis phase / continuum = 4 sqrt2/(2 pi) / sqrt(2/3)", phase_axis / phase_cont, 4 * math.sqrt(2) / (2 * math.pi) / math.sqrt(2 / 3), rel=1e-6)
check("discrete on-axis phase (mrad)", phase_axis * 1e3, 0.4503, rel=1e-3)
# Ring of linear magnetic material (M = chi H) round a wire carrying I: t = chi I V/(4 pi), exact.
chi, i_w = 0.37, 2.0
t_ring = 0.5 * sum(cross(r, sc(m, chi * i_w / (N_TURNS * I0)))[2] for r, m in core)   # core M scaled to chi I/(2 pi rho)
check("ring core round a wire: t = chi I V/(4 pi)", t_ring, chi * i_w * 2 * math.pi ** 2 * R0 * A0 ** 2 / (4 * math.pi), rel=1e-9)

# =============================================================================
# 4. Toroidal 'Aharonov-Casher': a moving point toroid in a static E field
# =============================================================================
print("\n# 4. Moving toroidal moment in a static E field")
# Rest frame: B' = -v x E/c^2, U' = -T.curl B'. Identity:
# curl(v x E) = v div E - (v.grad)E, so L = -U' = -mu0 rho T.v + (1/c^2) d(T.E)/dt.
QB, SIG = 1e-12, 0.2      # Gaussian blob of charge QB, width SIG (m), at the origin


def e_field(x):
    r = norm(x)
    if r < 1e-12:
        return (0.0, 0.0, 0.0)
    qin = QB * (math.erf(r / (math.sqrt(2) * SIG)) - math.sqrt(2 / math.pi) * r / SIG * math.exp(-r * r / (2 * SIG * SIG)))
    return sc(x, qin / (4 * math.pi * EPS0 * r ** 3))


def rho(x):
    return QB / (2 * math.pi * SIG * SIG) ** 1.5 * math.exp(-dot(x, x) / (2 * SIG * SIG))


T_V = (0.3, -0.5, 0.8)
V_V = (1.0, 2.0, -0.5)
X0 = (0.1, 0.15, -0.05)
hh = 1e-4


def grad(f, x):
    return tuple((f(add(x, sc(ei, hh))) - f(sub(x, sc(ei, hh)))) / (2 * hh) for ei in ((1, 0, 0), (0, 1, 0), (0, 0, 1)))


def curl_vxe(x):
    comp = lambda q: (lambda y: cross(V_V, e_field(y))[q])
    g = [grad(comp(q), x) for q in range(3)]
    return (g[2][1] - g[1][2], g[0][2] - g[2][0], g[1][0] - g[0][1])


lhs = -dot(T_V, sc(curl_vxe(X0), -1 / C ** 2))           # U' = -T.curl(B')
div_e = sum(grad(lambda y, q=q: e_field(y)[q], X0)[q] for q in range(3))
check("Gauss law in the blob: div E = rho/eps0", div_e, rho(X0) / EPS0, rel=1e-5)
rhs = (dot(T_V, V_V) * rho(X0) / EPS0 - dot(V_V, grad(lambda y: dot(T_V, e_field(y)), X0))) / C ** 2
check("U' = [(T.v) rho/eps0 - (v.grad)(T.E)]/c^2", lhs, rhs, rel=1e-5)


def phase_loop(radius, terms):
    """Oint L dt/hbar round a circle in the xy plane, speed |v| = 1 m/s along the path."""
    tot = mag = 0.0
    n = 2000
    for i in range(n):
        a = 2 * math.pi * (i + 0.5) / n
        x = (radius * math.cos(a), radius * math.sin(a), 0.05)
        tang = (-math.sin(a), math.cos(a), 0.0)
        dl = 2 * math.pi * radius / n
        if terms == "contact":
            d = -MU0 * rho(x) * dot(T_V, tang) * dl
        else:   # (1/c^2) grad(T.E).dl, a total derivative
            d = dot(grad(lambda y: dot(T_V, e_field(y)), x), tang) * dl / C ** 2
        tot += d
        mag += abs(d)
    return tot / HBAR, mag / HBAR


g_tot, g_mag = phase_loop(0.15, "grad")
check("total-derivative term round a loop through the blob / its |integrand| sum", g_tot / g_mag, 0.0, abs_=1e-6)
# Two-arm interferometer: arm 1 straight through the blob centre along x, arm 2 far away
# (rho = 0 there). Phase difference -(mu0/hbar) T_x Int rho dx = -(mu0/hbar) T_x QB/(2 pi SIG^2).
nx = 4000
arm = sum(-MU0 * rho((-2.0 + 4.0 * (i + 0.5) / nx, 0.0, 0.0)) * T_V[0] * (4.0 / nx) for i in range(nx)) / HBAR
check("contact phase on an arm through the blob = -(mu0/hbar) T_x Q/(2 pi sigma^2)", arm, -MU0 * T_V[0] * QB / (2 * math.pi * SIG ** 2) / HBAR, rel=1e-9)
c_out, c_mag = phase_loop(3.0, "contact")
check("contact phase on a loop at r = 3 m, where rho ~ 0, / the arm's", c_out / arm, 0.0, abs_=1e-30)
# Reciprocity with AB: for T = Phi A/mu0, mu0 rho T.L = Phi * (charge swept through the hole).
phi_t, a_t, l_t, rho_t = FLUX_QUANTUM, 1e-12, 1e-3, 5.0
T_thin = phi_t * a_t / MU0
check("mu0 rho T L = Phi (rho A L): toroidal AC = AB of the swept charge", MU0 * rho_t * T_thin * l_t, phi_t * rho_t * a_t * l_t, rel=1e-15)
check("one electron through an h/2e toroid's hole: pi", MU0 * T_thin * (E / (a_t * l_t)) * l_t / HBAR, math.pi, rel=1e-12)

print(f"\n{count - len(failures)}/{count} checks passed")
if failures:
    print("FAILED:", ", ".join(failures))
    sys.exit(1)
