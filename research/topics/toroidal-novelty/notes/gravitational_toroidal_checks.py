"""Checks for gravitational_toroidal.md: toroidal mass currents, anapoles and the
gravitomagnetic Aharonov-Bohm phase in linearised gravity.

Needs numpy (user site). Run as:  python3 -P gravitational_toroidal_checks.py  (-P keeps this folder off sys.path)
Every check prints computed vs expected and asserts a tolerance; exits 1 on any failure.

Sections
  A  gravitomagnetic AB phase round a toroidal mass current (Q3)
  B  the gravitational toroidal dipole cannot radiate; the EM one can (Q2)
  C  the first gravitational toroidal radiation is l = 2 and needs the stresses (Q2)
  D  a gravitational anapole: an exactly non-radiating conserved source (Q1)
  E  the toroidal term's share of the 1PN binary flux (Q2)
"""
import sys
from math import pi, sqrt

import numpy as np

# Constants as in crates/metrology (SI 2019, CODATA 2018)
C = 299_792_458.0
G = 6.674_30e-11
H = 6.626_070_15e-34
HBAR = H / (2 * pi)
E_CHARGE = 1.602_176_634e-19
M_E = 9.109_383_701_5e-31
AMU = 1.660_539_066_60e-27
M_RB87 = 86.909_180_529 * AMU      # as in crates/magnetics/src/monopole.rs
MU_0 = 1.256_637_062_12e-6
GM_PERMEABILITY = 4 * pi * G / C**2  # gravitics::gem::GM_PERMEABILITY

n_fail = 0


def check(name, got, want, rel=None, abs_=None):
    global n_fail
    if rel is not None:
        ok = abs(got / want - 1) <= rel if want != 0 else abs(got) <= rel
        tol = f"rel {rel:g}"
    else:
        ok = abs(got - want) <= abs_
        tol = f"abs {abs_:g}"
    print(f"[{'ok' if ok else 'FAIL'}] {name}: got {got:.6g}, want {want:.6g} ({tol})")
    if not ok:
        n_fail += 1


def check_true(name, cond, detail=""):
    global n_fail
    print(f"[{'ok' if cond else 'FAIL'}] {name} {detail}")
    if not cond:
        n_fail += 1


# ---------------------------------------------------------------------------------------------
# A. Gravitomagnetic AB phase round a toroidal mass current
#
# Slow particle of mass m: L ⊃ m c h_0i v^i, so Δφ = (m c/ħ)∮h_0i dx^i. Harmonic gauge, static:
# h_0i = −(4G/c³) a_i with a(x) = ∫ρv/|x−x'| d³x'. So |Δφ| = (4Gm/(ħc²)) Γ, Γ = ∮a·dl.
# ∇×∇×a = 4πρv, so b = ∇×a obeys Ampère with 4π: inside an ideal toroidal winding of total
# poloidal mass flux Ṁ (kg/s through the hole), b_φ = 2Ṁ/s, and Γ = 4πṀ(R − √(R²−a²)).
#   Δφ = 16πGmṀ(R − √(R² − a²))/(ħc²)
# ---------------------------------------------------------------------------------------------
print("== A. gravitomagnetic Aharonov-Bohm phase of a toroidal mass current")


def gm_ab_phase(m, mdot, major, minor):
    return 16 * pi * G * m * mdot * (major - sqrt(major**2 - minor**2)) / (HBAR * C**2)


def em_ab_phase_toroid(q, n_i, major, minor):
    # magnetics::toroidal::aharonov_bohm_phase(q, Φ) with Φ = μ0 N I (R − √(R²−a²))
    return q * MU_0 * n_i * (major - sqrt(major**2 - minor**2)) / HBAR


def toroid_segments(turns, mdot_per_turn, big_r, a, m):
    """Same winding as crates/magnetics/src/toroidal.rs tests: (midpoint, I dl) per segment."""
    mids, idl = [], []
    th = 2 * pi * np.arange(m + 1) / m
    for j in range(turns):
        ph = 2 * pi * (j + 0.5) / turns
        rho = big_r + a * np.cos(th)
        p = np.stack([rho * np.cos(ph), rho * np.sin(ph), a * np.sin(th)], axis=1)
        mids.append(0.5 * (p[1:] + p[:-1]))
        idl.append(mdot_per_turn * (p[1:] - p[:-1]))
    return np.concatenate(mids), np.concatenate(idl)


def loop(centre, radius, normal_axis, m):
    """Circle as (midpoints, dl); normal_axis 'y' puts it in the xz-plane."""
    t = 2 * pi * np.arange(m + 1) / m
    if normal_axis == "y":
        p = np.stack([centre[0] + radius * np.cos(t), np.full_like(t, centre[1]), centre[2] + radius * np.sin(t)], axis=1)
    else:  # 'z': horizontal circle
        p = np.stack([centre[0] + radius * np.cos(t), centre[1] + radius * np.sin(t), np.full_like(t, centre[2])], axis=1)
    return 0.5 * (p[1:] + p[:-1]), p[1:] - p[:-1]


def neumann(src, dst):
    """Γ = Σ_src Σ_dst (I dl · dl_C)/|x − x'|, the mass-current analogue of Neumann's formula."""
    (xs, ids), (xc, dlc) = src, dst
    tot = 0.0
    for k in range(0, len(xc), 50):
        d = np.linalg.norm(xc[k:k + 50, None, :] - xs[None, :, :], axis=2)
        tot += np.sum((dlc[k:k + 50] @ ids.T) / d)
    return tot


turns, R_, a_ = 96, 1.0, 0.3
src = toroid_segments(turns, 1.0, R_, a_, 120)     # 1 kg/s per turn: Ṁ = 96 kg/s
mdot = turns * 1.0
want_gamma = 4 * pi * mdot * (R_ - sqrt(R_**2 - a_**2))
g1 = abs(neumann(src, loop((R_, 0.0, 0.0), 0.45, "y", 400)))
g2 = abs(neumann(src, loop((R_, 0.0, 0.1), 0.6, "y", 400)))   # a different linking loop
g0 = abs(neumann(src, loop((0.0, 0.0, 2.0), 3.0, "z", 400)))  # does not link the tube
check("A1 Γ for a loop linking the tube = 4πṀ(R − √(R²−a²))", g1, want_gamma, rel=5e-3)
check("A1 a second linking loop gives the same Γ (topological)", g2, want_gamma, rel=5e-3)
check_true("A1 a loop not linking the tube sees Γ ≈ 0", g0 < 1e-3 * want_gamma, f"(Γ/Γ_link = {g0 / want_gamma:.2e})")

# Thin-torus limit in terms of the toroidal moment T_g = ṀπRa²/2 (toroidal_solenoid_moment form)
a_thin = 0.01
t_g = mdot * pi * R_ * a_thin**2 / 2
check("A2 thin torus: Δφ → 16GmT_g/(ħc²R²)", gm_ab_phase(M_RB87, mdot, R_, a_thin),
      16 * G * M_RB87 * t_g / (HBAR * C**2 * R_**2), rel=1e-4)

# Infinite rotating cylinder limit vs Mo & Modesto (arXiv:2608.15587) eq. 4.28: 8πGmμωρ0²/(ħc²).
# Shell of linear density μ, angular velocity ω, radius ρ0: mass current per length K = μω/(2π);
# a torus of large R with Ṁ = 2πR K and minor radius ρ0 tends to the cylinder.
mu_lin, omega, rho0, big = 2.0, 50.0, 0.1, 1e4
check("A3 cylinder limit equals Mo & Modesto eq. 4.28", gm_ab_phase(M_RB87, 2 * pi * big * mu_lin * omega / (2 * pi), big, rho0),
      8 * pi * G * M_RB87 * mu_lin * omega * rho0**2 / (HBAR * C**2), rel=1e-6)

# The naive analogue (μ0 → 4πG/c², q → m, NI → Ṁ) is 4× too small: gravity's spin-2 factor in a
# potential. Lynden-Bell & Katz (arXiv:1205.1616, eqs. 1.5, 1.7): curl H = −16πJ against EM's 4πj.
naive = M_RB87 * GM_PERMEABILITY * mdot * (R_ - sqrt(R_**2 - a_**2)) / HBAR
check("A4 gravity / naive analogue = 4 (LB&K 16π vs 4π)", gm_ab_phase(M_RB87, mdot, R_, a_) / naive, 16 * pi / (4 * pi), rel=1e-12)

# A5. Sizes. (i) The electron mass current of the toroid in magnetics::toroidal's doctest
# (24 turns of 1 A, R = 1 m, a = 0.3 m): Ṁ = 24 m_e/e.
mdot_e = 24 * 1.0 * M_E / E_CHARGE
phi_e_rb = gm_ab_phase(M_RB87, mdot_e, 1.0, 0.3)
phi_e_em = em_ab_phase_toroid(E_CHARGE, 24.0, 1.0, 0.3)
phi_e_gm = gm_ab_phase(M_E, mdot_e, 1.0, 0.3)
print(f"    electron mass current of the 24 A-turn toroid: {mdot_e:.4e} kg/s")
print(f"    GM AB phase, ⁸⁷Rb: {phi_e_rb:.4e} rad;  electron EM AB phase: {phi_e_em:.4e} rad;  electron GM phase: {phi_e_gm:.4e} rad")
check("A5 electron: GM/EM phase ratio = 16πG m_e²/(μ0 e² c²)",
      phi_e_gm / phi_e_em, 16 * pi * G * M_E**2 / (MU_0 * E_CHARGE**2 * C**2), rel=1e-12)
# = 4α_G/α with α_G = Gm_e²/(ħc) = 1.75181e-45 (gravitics::gem test) and α = 7.2973525693e-3
check("A5 that ratio is 4α_G/α", phi_e_gm / phi_e_em, 4 * 1.751_81e-45 / 7.297_352_5693e-3, rel=1e-4)

# (ii) N spinning rings round a torus (total mass M, rim speed v, ring radius a, torus radius R),
# thin limit: Δφ = 4 G m M v a/(ħ c² R). Inputs below are ASSUMED design values, not sourced.
def rings_phase(m, m_tot, v, a, r):
    return 4 * G * m * m_tot * v * a / (HBAR * C**2 * r)


# check rings formula against the continuum one: N rings of mass M/N, radius a, speed v carry
# mass current M v/(2π a N) each, Ṁ = M v/(2π a)
m_tot, v_rim = 1.25, 100.0
check("A6 ring formula = continuum formula (thin)", rings_phase(M_RB87, m_tot, v_rim, 0.01, 1.0),
      gm_ab_phase(M_RB87, m_tot * v_rim / (2 * pi * 0.01), 1.0, 0.01), rel=1e-4)
sigma_shot = 1 / (0.8 * sqrt(1e5))    # CARIOQA per-shot phase noise, gravimetry::atom test (3.95e-3)
check("A7 CARIOQA per-shot phase noise (crate value 3.95e-3)", sigma_shot, 3.95e-3, rel=0.01)
for label, mt, v in [("1.25 kg (Overstreet's source mass) at 100 m/s", 1.25, 100.0),
                     ("1 t at 1 km/s", 1e3, 1e3)]:
    p = rings_phase(M_RB87, mt, v, 0.5, 1.0)   # a/R = 1/2
    shots = (sigma_shot / p) ** 2
    print(f"    {label}, a/R = 0.5: Δφ = {p:.3e} rad; shots for SNR 1: {shots:.2e} ({shots * 3.5 / 3.156e7:.2e} yr at 3.5 s)")
p_big = rings_phase(M_RB87, 1e3, 1e3, 0.5, 1.0)
check("A8 1 t at 1 km/s, a/R = 1/2: Δφ", p_big, 2.0325e-12, rel=1e-3)
# Earth-rotation Sagnac phase over 1 m², 2mΩA/ħ, for scale (Ω⊕ = 7.292115e-5 rad/s, IERS)
sagnac = 2 * M_RB87 * 7.292_115e-5 * 1.0 / HBAR
print(f"    Earth Sagnac phase over 1 m²: {sagnac:.3e} rad; ratio GM/Sagnac = {p_big / sagnac:.2e}")
check("A9 GM phase / Earth Sagnac (1 m²)", p_big / sagnac, 1.0e-17, rel=0.1)
# Overstreet et al. 2022 scalar AB rate, gravimetry::atom::aharonov_bohm_rate: 1.63 rad/s
check("A10 scalar AB rate, Rb at 7 cm from 1.25 kg", M_RB87 * G * 1.25 / (HBAR * 0.07), 1.631, rel=1e-3)

# ---------------------------------------------------------------------------------------------
# Grid tools for B–D (units G = c = 1 inside these sections; e^(−iωt), ∂_t → −iω, k = ω)
# ---------------------------------------------------------------------------------------------
n, Lb = 80, 4.0
h = 2 * Lb / n
ax = -Lb + h * np.arange(n)
X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
XS = [X, Y, Z]
R2 = X**2 + Y**2 + Z**2
dV = h**3
q1 = 2 * pi * np.fft.fftfreq(n, d=h)
QX, QY, QZ = np.meshgrid(q1, q1, q1, indexing="ij")
QS = [QX, QY, QZ]


def d(f, i):
    return np.fft.ifftn(1j * QS[i] * np.fft.fftn(f))


def lap(f):
    return np.fft.ifftn(-(QX**2 + QY**2 + QZ**2) * np.fft.fftn(f))


def integ(f):
    return np.sum(f) * dV


def xhat2(i, j):
    return XS[i] * XS[j] - (R2 / 3 if i == j else 0.0)


def xhat3(i, j, k):
    t = XS[i] * XS[j] * XS[k]
    t = t - R2 / 5 * ((XS[k] if i == j else 0) + (XS[j] if i == k else 0) + (XS[i] if j == k else 0))
    return t


def blob(c, w):
    return np.exp(-((X - c[0]) ** 2 + (Y - c[1]) ** 2 + (Z - c[2]) ** 2) / w**2)


rng = np.random.default_rng(20261027)


def random_sym_T(nblobs=4, w=0.45):
    T = [[0j * X for _ in range(3)] for _ in range(3)]
    for _ in range(nblobs):
        c = rng.uniform(-0.6, 0.6, 3)
        b = blob(c, w)
        m = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
        m = m + m.T
        for i in range(3):
            for j in range(3):
                T[i][j] = T[i][j] + m[i, j] * b
    return T


# ---------------------------------------------------------------------------------------------
# B. l = 1. Blanchet LRR eq. 138 (1PN mass moments, σ = (T00 + Tii)/c², σi = T0i/c) at l = 1:
#   I_i = ∫x_i σ + (1/10c²)∂t²∫r² x_i σ − (6/5c²)∂t∫x̂_ij σ_j.
# For a conserved source the O(1/c²) part must vanish (no dipole radiation). With T_ij given,
# Σ_i = ∂_j T_ij/(iω), ρ = ∂_k Σ_k/(iω). The O(1/c²) part is
#   ∫x_i T_kk − (ω²/10)∫r² x_i ρ + (6iω/5)∫x̂_ij Σ_j.
# EM: p_i = ∫x_i ρ + (1/10c²)∂t²∫r² x_i ρ − (3/10c²)∂t∫x̂_ij J_j, and the c⁻² part is −Ṫ_i/c² with
# Dubovik's T_i = (1/10)∫[(x·J)x_i − 2r²J_i]: it does not vanish.
# ---------------------------------------------------------------------------------------------
print("== B. the toroidal dipole: gravity cannot radiate it, EM can")
omega = 0.7
T = random_sym_T()
Sig = [sum(d(T[i][j], j) for j in range(3)) / (1j * omega) for i in range(3)]
rho = sum(d(Sig[k], k) for k in range(3)) / (1j * omega)
Tkk = T[0][0] + T[1][1] + T[2][2]
worst, worst_no_stress = 0.0, 0.0
for i in range(3):
    t1 = integ(XS[i] * Tkk)
    t2 = -(omega**2 / 10) * integ(R2 * XS[i] * rho)
    t3 = (6j * omega / 5) * integ(sum(xhat2(i, j) * Sig[j] for j in range(3)))
    scale = abs(t1) + abs(t2) + abs(t3)
    worst = max(worst, abs(t1 + t2 + t3) / scale)
    worst_no_stress = max(worst_no_stress, abs(t2 + t3) / scale)
    # the leading term ∫x_i ρ vanishes too (momentum conservation)
    lead = integ(XS[i] * rho)
    check_true(f"B1 ∫x_{i}ρ = 0 (relative)", abs(lead) / integ(np.abs(XS[i] * rho)) < 1e-10, f"({abs(lead) / integ(np.abs(XS[i] * rho)):.1e})")
check_true("B2 gravity: O(1/c²) dipole terms cancel", worst < 1e-10, f"(worst residual / scale = {worst:.1e})")
check_true("B3 ... but only with the stress T_kk in σ", worst_no_stress > 0.05, f"(without it: {worst_no_stress:.2f})")

# EM: random compact current, ρ_e = ∂·J/(iω)
J = [sum((rng.normal() + 1j * rng.normal()) * blob(rng.uniform(-0.6, 0.6, 3), 0.45) for _ in range(4)) for _ in range(3)]
rho_e = sum(d(J[k], k) for k in range(3)) / (1j * omega)
xJ = X * J[0] + Y * J[1] + Z * J[2]
for i in range(3):
    corr = -(omega**2 / 10) * integ(R2 * XS[i] * rho_e) + (3j * omega / 10) * integ(sum(xhat2(i, j) * J[j] for j in range(3)))
    t_dub = 0.1 * integ(xJ * XS[i] - 2 * R2 * J[i])
    check(f"B4 EM: O(1/c²) dipole correction = −Ṫ_{i} = iωT_{i} (re)", corr.real, (1j * omega * t_dub).real, rel=1e-8)
check_true("B5 EM toroidal dipole of a generic current is nonzero", abs(t_dub) > 1e-3)

# ---------------------------------------------------------------------------------------------
# C. l = 2 for a closed (divergence-free) axisymmetric poloidal mass current, ρ_osc = 0.
# Construction: U = ∇×(h(s,z)(−y, x, 0)) is poloidal and divergence-free; T_ij = ½(x_iU_j + x_jU_i).
# Then ∂_jT_ij = ½(x·∇U_i + 4U_i) and ∂_i∂_jT_ij = 0, so Σ = ∂_jT_ij/(iω) is a closed poloidal mass
# current and ρ_osc = 0. All compact. (T_ij = ∂_(iV_j) with ∇·V = 0 also gives ρ_osc = 0, but its
# l = 2 moments vanish identically, so it starts at higher order.)
# Prediction (Blanchet eq. 138, l = 2, ρ static): I_ij = (1/c²)[∫x̂_ij T_kk + (4/3)iω A_ij],
# A_ij = ∫(x·Σ)x̂_ij — the toroidal term — equivalently −(1/c²)[(1/3)∫x̂_ij T_kk + (8/3)S_ij],
# S_ij = ∫x_k T_k<i x_j>. Current moments vanish for a poloidal axisymmetric current.
# Exact linearised power: P = (Gω²/4πc⁵)∮Λ_ij,kl T̃_ij T̃*_kl dΩ, T̃(k) = ∫T e^(−ik·x).
# ---------------------------------------------------------------------------------------------
print("== C. l = 2 toroidal radiation needs the stresses")


def poloidal_T(cells, w=0.6):
    """cells: list of (z0, sign), each a Gaussian vortex-ring cell on the axis (a closed poloidal
    mass current, like Hill's vortex). Returns T_ij and ∂_jT_ij (= iωΣ_i)."""
    hfun = sum(sg * np.exp(-(X**2 + Y**2 + (Z - z0) ** 2) / w**2) for z0, sg in cells)
    Phi = [-Y * hfun, X * hfun, 0 * X]
    U = [d(Phi[2], 1) - d(Phi[1], 2), d(Phi[0], 2) - d(Phi[2], 0), d(Phi[1], 0) - d(Phi[0], 1)]
    Tm = [[0.5 * (XS[i] * U[j] + XS[j] * U[i]) for j in range(3)] for i in range(3)]
    divT = [sum(d(Tm[i][j], j) for j in range(3)) for i in range(3)]
    return Tm, divT


def ttproj_power(Tm, k, n_mu=24, n_phi=1):
    """∮ Λ T̃ T̃* dΩ (Gauss–Legendre in cos θ; n_phi = 1 assumes axisymmetry)."""
    mus, wts = np.polynomial.legendre.leggauss(n_mu)
    tot = 0.0
    for mu, wt in zip(mus, wts):
        st = sqrt(1 - mu * mu)
        for ip in range(n_phi):
            ph = 2 * pi * (ip + 0.5) / n_phi
            nv = np.array([st * np.cos(ph), st * np.sin(ph), mu])
            ph_fac = np.exp(-1j * k * (nv[0] * X + nv[1] * Y + nv[2] * Z))
            Tt = np.array([[integ(Tm[i][j] * ph_fac) for j in range(3)] for i in range(3)])
            P = np.eye(3) - np.outer(nv, nv)
            A = P @ Tt @ P
            A = A - 0.5 * P * np.trace(A)
            tot += wt * (2 * pi / n_phi) * np.sum(np.abs(A) ** 2)
    return tot


def power_exact(Tm, k, **kw):
    return k**2 / (4 * pi) * ttproj_power(Tm, k, **kw)


def quad_power(I, k):
    return (1 / 5) * 0.5 * k**6 * np.sum(np.abs(I) ** 2)


# divergence check of the construction
Tm, LV = poloidal_T([(0.6, 1.0), (-0.6, -1.0)])   # two opposed cells: toroidal dipoles ±T at z = ±0.6
divS = sum(d(LV[k], k) for k in range(3))
check_true("C0 Σ is divergence-free (ρ_osc = 0)", np.max(np.abs(divS)) < 1e-8 * np.max(np.abs(LV[2])), f"({np.max(np.abs(divS)) / np.max(np.abs(LV[2])):.1e})")
Uphi = -np.sin(np.arctan2(Y, X)) * LV[0] + np.cos(np.arctan2(Y, X)) * LV[1]
check_true("C0 Σ is poloidal (no azimuthal part)", np.max(np.abs(Uphi)) < 1e-8 * np.max(np.abs(LV[2])))


def l2_moments(Tm, LV, k):
    Sg = [LV[i] / (1j * k) for i in range(3)]   # Σ = ∂_jT_ij/(iω)
    tkk = Tm[0][0] + Tm[1][1] + Tm[2][2]
    xS = X * Sg[0] + Y * Sg[1] + Z * Sg[2]
    I_t = np.zeros((3, 3), complex)
    I_s = np.zeros((3, 3), complex)
    I_alt = np.zeros((3, 3), complex)
    for i in range(3):
        for j in range(3):
            a_ij = integ(xS * xhat2(i, j))
            m_kk = integ(xhat2(i, j) * tkk)
            s_raw = lambda p, q: integ(sum(XS[kk] * Tm[kk][p] for kk in range(3)) * XS[q])
            s_ij = 0.5 * (s_raw(i, j) + s_raw(j, i)) - (sum(s_raw(l, l) for l in range(3)) / 3 if i == j else 0)
            I_t[i, j] = (4 / 3) * 1j * k * a_ij
            I_s[i, j] = m_kk
            I_alt[i, j] = -(1 / 3) * m_kk - (8 / 3) * s_ij
    return I_t, I_s, I_alt


for kk in [0.05, 0.1]:
    I_t, I_s, I_alt = l2_moments(Tm, LV, kk)
    I_full = I_t + I_s
    check_true(f"C1 k = {kk}: the two forms of I_ij agree (Ȧ = ∫x̂T_kk + 2S)", np.sum(np.abs(I_full - I_alt)) < 1e-8 * np.sum(np.abs(I_full)), f"({np.sum(np.abs(I_full - I_alt)) / np.sum(np.abs(I_full)):.1e})")
    pe = power_exact(Tm, kk)
    check(f"C2 k = {kk}: exact power / l = 2 prediction (toroidal + stress)", pe / quad_power(I_full, kk), 1.0, rel=0.05 if kk > 0.07 else 0.015)
    print(f"    toroidal-only prediction / exact = {quad_power(I_t, kk) / pe:.3f};  stress-only / exact = {quad_power(I_s, kk) / pe:.3f}")
    check_true(f"C3 k = {kk}: toroidal term alone is wrong by > 2×", abs(quad_power(I_t, kk) / pe - 1) > 1.0)
# fixed T (fixed V): P ∝ ω⁶; fixed mass current Σ ∝ V/ω: P ∝ ω⁸
p1, p2 = power_exact(Tm, 0.03), power_exact(Tm, 0.06)
check("C4 pair of opposed toroids at fixed stress: P(2ω)/P(ω) → 2⁶", p2 / p1, 64.0, rel=0.01)
print(f"    so at fixed mass current (V ∝ ω) the ratio is {p2 / p1 * 4:.1f} ≈ 2⁸")
# single toroid (z-odd): I_ij = 0 by symmetry, leading is l = 3 at 1/c²: ω⁸ at fixed stress
Tm1, LV1 = poloidal_T([(0.0, 1.0)])
I_t1, I_s1, _ = l2_moments(Tm1, LV1, 0.05)
check_true("C5 single toroid: l = 2 moment vanishes", np.sum(np.abs(I_t1 + I_s1)) < 1e-8 * np.sum(np.abs(I_t + I_s)))
# This centred cell, with the stress T = x_(iU_j), radiates nothing at any frequency: U is built from
# a radial Gaussian, so T̃_ij(q) is made of q_i, q_j and δ_ij terms only (a longitudinal stress)
# and the TT projection kills it. Add a divergence-free stress S_ij = ε_ikz ε_jmz ∂_k∂_m g: Σ and ρ
# are unchanged, the same mass current now radiates, at l = 2 with I_ij = ∫x̂_ij S_kk (A_ij = 0).
p_pair = power_exact(Tm, 0.05)
check_true("C6 the centred cell (same kind of closed mass current) radiates nothing", power_exact(Tm1, 0.05) < 1e-15 * p_pair, f"(P/P_pair = {power_exact(Tm1, 0.05) / p_pair:.1e})")
Sg1 = [LV1[i] / 1j for i in range(3)]   # Σ at ω = 1
t_dip = 0.5 * integ((X * Sg1[0] + Y * Sg1[1] + Z * Sg1[2]) * Z)   # toroidal dipole ½∫(x·Σ)z
t_scale = 0.5 * integ(np.abs((X * Sg1[0] + Y * Sg1[1] + Z * Sg1[2]) * Z))
print(f"    centred cell: toroidal dipole T_z = ½∫(x·Σ)z = {t_dip.real:.4g}{t_dip.imag:+.4g}i (scale {t_scale:.4g})")
check_true("C6 ... although its toroidal dipole is nonzero (EM: a lone T radiates)", abs(t_dip) > 0.1 * t_scale)
g = blob((0.0, 0.0, 0.3), 0.6)
Sd = [[(-d(d(g, 1), 1) if (i, j) == (0, 0) else -d(d(g, 0), 0) if (i, j) == (1, 1) else d(d(g, 0), 1) if (i, j) in [(0, 1), (1, 0)] else 0 * g) for j in range(3)] for i in range(3)]
divSd = max(np.max(np.abs(sum(d(Sd[i][j], j) for j in range(3)))) for i in range(3))
check_true("C7 S_ij is divergence-free", divSd < 1e-10 * np.max(np.abs(Sd[0][0])))
Tm2 = [[Tm1[i][j] + Sd[i][j] for j in range(3)] for i in range(3)]
skk = Sd[0][0] + Sd[1][1]
I_S = np.array([[integ(xhat2(i, j) * skk) for j in range(3)] for i in range(3)])
for kk in [0.05, 0.1]:
    check(f"C8 k = {kk}: same mass current plus S radiates (G/5c⁵)⟨I⃛²⟩, I = ∫x̂S_kk", power_exact(Tm2, kk) / quad_power(I_S, kk), 1.0, rel=0.02 if kk > 0.07 else 0.005)

# normalisation check of power_exact on a generic source with oscillating density: P → Newtonian
Tg = random_sym_T()
kq = 0.02
I_N = np.array([[-(2 / kq**2) * integ(Tg[i][j] - ((Tg[0][0] + Tg[1][1] + Tg[2][2]) / 3 if i == j else 0)) for j in range(3)] for i in range(3)])
check("C9 generic source: exact power → (G/5c⁵)⟨I⃛²⟩ with I = −(2/ω²)∫T_<ij>", power_exact(Tg, kq, n_mu=16, n_phi=16) / quad_power(I_N, kq), 1.0, rel=2e-3)

# ---------------------------------------------------------------------------------------------
# D. A gravitational anapole (Devaney–Wolf construction). T_ij = (∇² + k²)F_ij with compact F;
# Σ, ρ by conservation. Then T̃_ij = (k² − q²)F̃_ij vanishes on the shell q = k: no radiation, and
# the outgoing solution of (∇² + k²)u = T is u = F, zero outside the source (h̄ = −16πG F/c⁴).
# Its Newtonian quadrupole −(2/ω²)∫T_<ij> = −2∫F_<ij> (c = 1) is not zero: the 1/c² toroidal and
# stress terms cancel it.
# ---------------------------------------------------------------------------------------------
print("== D. gravitational anapole: an exactly non-radiating conserved source")
kD = 0.5
F = random_sym_T(nblobs=3, w=0.5)
TD = [[lap(F[i][j]) + kD**2 * F[i][j] for j in range(3)] for i in range(3)]
p_ana = power_exact(TD, kD, n_mu=12, n_phi=12)
p_ref = power_exact(F, kD, n_mu=12, n_phi=12)   # F itself as a source, same size and shape
check_true("D1 anapole source radiates nothing", p_ana < 1e-12 * p_ref, f"(P/P_ref = {p_ana / p_ref:.1e})")
I_ND = np.array([[-(2 / kD**2) * integ(TD[i][j] - ((TD[0][0] + TD[1][1] + TD[2][2]) / 3 if i == j else 0)) for j in range(3)] for i in range(3)])
I_F = np.array([[-2 * integ(F[i][j] - ((F[0][0] + F[1][1] + F[2][2]) / 3 if i == j else 0)) for j in range(3)] for i in range(3)])
check("D2 its Newtonian quadrupole is −2∫F_<ij> (nonzero)", np.sum(np.abs(I_ND)), np.sum(np.abs(I_F)), rel=1e-8)
check_true("D2 ... and would radiate on its own", quad_power(I_ND, kD) > 1e-3 * p_ref)
# exterior field: u(x0) = −∫ e^(ik|x0−x|)/(4π|x0−x|) T(x) d³x at x0 = (6, 0, 0)
x0 = np.array([6.0, 0.0, 0.0])
dist = np.sqrt((X - x0[0]) ** 2 + (Y - x0[1]) ** 2 + (Z - x0[2]) ** 2)
Gk = np.exp(1j * kD * dist) / (4 * pi * dist)
u_ext = abs(integ(Gk * TD[0][1]))
u_scale = integ(np.abs(Gk * TD[0][1]))
check_true("D3 exterior field of the anapole vanishes (|∫G T| / ∫|G T|)", u_ext < 1e-8 * u_scale, f"({u_ext / u_scale:.1e})")

# ---------------------------------------------------------------------------------------------
# E. 1PN flux of a circular binary from Blanchet–Damour's point-particle quadrupole
#   I_ij = Σ m_A[1 + (3/2)v_A²/c² − Gm_B/(r c²)] x̂_A + (1/14c²) d²/dt² Σ m_A r_A² x̂_A
#          − (20/21c²) d/dt Σ m_A v_A^k x̂_A^ijk      (the last is the toroidal-type term)
# plus the Newtonian mass octupole and current quadrupole, with F = (G/c⁵)[(1/5)I⁽³⁾² +
# (1/c²)((1/189)I⁽⁴⁾² + (16/45)J⁽³⁾²)], harmonic-coordinate Kepler ω² = (Gm/r³)[1 + (−3+ν)γ].
# Target: F = (32/5)ν²x⁵[1 + (−1247/336 − 35ν/12)x] (Blanchet LRR, opened, dl_grav/1310.1528).
# ---------------------------------------------------------------------------------------------
print("== E. the toroidal term's share of the 1PN binary flux")


def stf2(a):
    return 0.5 * (a + np.swapaxes(a, -1, -2)) - np.eye(3) * (np.trace(a, axis1=-2, axis2=-1)[..., None, None] / 3)


def stf3(x):
    r2 = np.sum(x * x, axis=-1)
    t = np.einsum("ti,tj,tk->tijk", x, x, x)
    e = np.eye(3)
    t -= (r2[:, None, None, None] / 5) * (np.einsum("ij,tk->tijk", e, x) + np.einsum("ik,tj->tijk", e, x) + np.einsum("jk,ti->tijk", e, x))
    return t


def tdiff(f, w, order):
    """Spectral time derivative of samples over one orbital period."""
    nt = f.shape[0]
    freq = 2 * pi * np.fft.fftfreq(nt, d=(2 * pi / w) / nt)
    sh = (nt,) + (1,) * (f.ndim - 1)
    return np.real(np.fft.ifft((1j * freq.reshape(sh)) ** order * np.fft.fft(f, axis=0), axis=0))


def flux_ratio(nu, x, toroidal=True):
    m, gm, w = 1.0, 1.0, 1.0   # G = m = ω = 1; c set by x = (Gmω/c³)^(2/3)
    cl = (gm * w) ** (1 / 3) / sqrt(x)
    # 1PN Kepler: solve ω² = (Gm/r³)(1 + (−3+ν)γ), γ = Gm/(r c²), for r
    r = (gm / w**2) ** (1 / 3)
    for _ in range(50):
        r = (gm / w**2 * (1 + (-3 + nu) * gm / (r * cl**2))) ** (1 / 3)
    dm = sqrt(1 - 4 * nu)
    m1, m2 = m * (1 + dm) / 2, m * (1 - dm) / 2
    nt = 64
    t = (2 * pi / w) * np.arange(nt) / nt
    xr = np.stack([r * np.cos(w * t), r * np.sin(w * t), 0 * t], axis=1)
    vr = np.stack([-r * w * np.sin(w * t), r * w * np.cos(w * t), 0 * t], axis=1)
    bodies = [(m1, m2 / m * xr, m2 / m * vr, m2), (m2, -m1 / m * xr, -m1 / m * vr, m1)]
    I = np.zeros((nt, 3, 3))
    r2x = np.zeros((nt, 3, 3))
    vx3 = np.zeros((nt, 3, 3))
    I3 = np.zeros((nt, 3, 3, 3))
    Jc = np.zeros((nt, 3, 3))
    eps = np.zeros((3, 3, 3))
    eps[0, 1, 2] = eps[1, 2, 0] = eps[2, 0, 1] = 1
    eps[0, 2, 1] = eps[2, 1, 0] = eps[1, 0, 2] = -1
    for ma, xa, va, mb in bodies:
        xx = stf2(np.einsum("ti,tj->tij", xa, xa))
        v2 = np.sum(va * va, axis=1)
        ra2 = np.sum(xa * xa, axis=1)
        I += ma * (1 + 1.5 * v2 / cl**2 - gm * mb / m / (r * cl**2))[:, None, None] * xx
        r2x += ma * ra2[:, None, None] * xx
        vx3 += ma * np.einsum("tk,tijk->tij", va, stf3(xa))
        I3 += ma * stf3(xa)
        Jc += ma * stf2(np.einsum("abi,ta,tb,tj->tij", eps, xa, va, xa))
    I = I + tdiff(r2x, w, 2) / (14 * cl**2)
    if toroidal:
        I = I - (20 / 21) * tdiff(vx3, w, 1) / cl**2
    f = (1 / 5) * np.mean(np.sum(tdiff(I, w, 3) ** 2, axis=(1, 2)))
    f += (1 / 189) * np.mean(np.sum(tdiff(I3, w, 4) ** 2, axis=(1, 2, 3))) / cl**2
    f += (16 / 45) * np.mean(np.sum(tdiff(Jc, w, 3) ** 2, axis=(1, 2))) / cl**2
    f /= cl**5
    return f / (32 / 5 * nu**2 * x**5 * cl**5)   # over (32/5)(c⁵/G)ν²x⁵


def c1(nu, toroidal=True):
    x = 1e-4
    f1, f2 = flux_ratio(nu, x, toroidal), flux_ratio(nu, 2 * x, toroidal)
    return (4 * f1 - f2 - 3) / (2 * x)


for nu in [0.25, 0.1, 0.01]:
    full = c1(nu)
    check(f"E1 ν = {nu}: 1PN flux coefficient = −1247/336 − 35ν/12", full, -1247 / 336 - 35 * nu / 12, rel=1e-5)
    share = full - c1(nu, toroidal=False)
    check(f"E2 ν = {nu}: share of the −(20/21) toroidal-type term = −(32/21)(1 − 3ν)", share, -(32 / 21) * (1 - 3 * nu), rel=1e-5)
print(f"    equal masses: toroidal share {-(32 / 21) * 0.25:.4f} of {-1247 / 336 - 35 / 48:.4f}; test-mass limit {-(32 / 21):.4f} of {-1247 / 336:.4f}")

print()
if n_fail:
    print(f"{n_fail} check(s) FAILED")
    sys.exit(1)
print("all checks passed")
