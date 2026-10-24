"""Recomputed numbers for mhd_theory_and_stability.md. Run from outside this folder:
    python3 -I /path/to/mhd_theory_checks.py
Standard library only. Every assert compares against a published value (tolerance noted) or,
where marked "identity", against an exact limit of the formula itself."""
import cmath
import math

# CODATA 2018
MU0 = 1.25663706212e-6      # H/m
E = 1.602176634e-19         # C
MP = 1.67262192369e-27      # kg
ME = 9.1093837015e-31       # kg
M_ALPHA = 6.6446573357e-27  # kg (He-4 nucleus)
U = 1.66053906660e-27       # kg
EPS0 = 8.8541878128e-12     # F/m
C = 299792458.0             # m/s
KB = 1.380649e-23           # J/K


def close(x, ref, rel):
    return abs(x - ref) <= rel * abs(ref)


# ---------------------------------------------------------------- Q1 basic MHD numbers
def alfven(b, rho):
    return b / math.sqrt(MU0 * rho)


def magnetosonic(va, cs, theta):
    """(fast, slow) phase speeds at angle theta between k and B."""
    s = va * va + cs * cs
    d = math.sqrt(s * s - 4 * va * va * cs * cs * math.cos(theta) ** 2)
    return math.sqrt((s + d) / 2), math.sqrt(max(s - d, 0.0) / 2)


def beta(n, t_ev, b):
    return 2 * MU0 * n * t_ev * E / b ** 2


def inertial_length(n, m, z=1):
    return C / math.sqrt(n * (z * E) ** 2 / (EPS0 * m))


def gyroradius(t_ev, b, m=MP, z=1):
    return math.sqrt(2 * t_ev * E / m) / (z * E * b / m)


print("== Q1 basic MHD ==")
# PlasmaPy Alfven_speed docstring: B = 0.014 T, n = 5e19 m^-3, rho = n (m_p + m_e) -> 43173.870 m/s;
# He-4 2+ (rho = n m_alpha + 2 n m_e) -> 21664.18 m/s
va_p = alfven(0.014, 5e19 * (MP + ME))
va_he = alfven(0.014, 5e19 * (M_ALPHA + 2 * ME))
print(f"v_A proton {va_p:.2f} m/s, He2+ {va_he:.2f} m/s")
assert close(va_p, 43173.870, 1e-6) and close(va_he, 21664.18, 1e-5)
# Coronal loop, arXiv:2605.15057: B = 11 G, rho = 1.67e-15 g/cm^3 -> 759 km/s outside; 2x denser inside -> 537 km/s
va_out, va_in = alfven(11e-4, 1.67e-12), alfven(11e-4, 2 * 1.67e-12)
print(f"corona v_A {va_out/1e3:.1f} km/s outside, {va_in/1e3:.1f} km/s inside")
assert close(va_out, 759e3, 2e-3) and close(va_in, 537e3, 2e-3)
# Active-region loops, arXiv:2106.10884: log n_e = 9.2, mean B ~110 G; v_A range 2000-9400 km/s
va_ar = alfven(110e-4, 10 ** 9.2 * 1e6 * MP)
print(f"active-region loop v_A {va_ar/1e3:.0f} km/s (published range 2000-9400)")
assert 2000e3 < va_ar < 9400e3
# practical-unit constants recomputed from SI (the NRL Formulary prints 2.18e11 and 2.28e7 in cgs; not re-read here)
k_va = 1 / math.sqrt(MU0 * MP)                     # v_A = k_va * B[T] / sqrt(mu n[m^-3]) m/s
k_di = C / math.sqrt(1e6 * E ** 2 / (EPS0 * MP))   # d_i = k_di / sqrt(n[cm^-3]) m, protons
print(f"v_A = {k_va:.4e} B[T]/sqrt(mu n[m^-3]) m/s = {k_va*1e-4/1e3*100:.4e} B[G]/sqrt(mu n[cm^-3]) cm/s; "
      f"d_i = {k_di*100:.4e}/sqrt(n[cm^-3]) cm")
assert close(k_va * 1e-4 / 1e3 * 100, 2.18e11, 1e-3) and close(k_di * 100, 2.28e7, 2e-3)
# magnetic pressure, identity
print(f"B^2/2mu0 at 1 T = {1/(2*MU0)/1e3:.3f} kPa; at SPARC 12.2 T = {12.2**2/(2*MU0)/1e6:.1f} MPa")
# fast/slow identities: theta = 0 -> (max, min)(vA, cs); theta = 90 deg -> sqrt(vA^2+cs^2), 0
f0, s0 = magnetosonic(3.0, 1.0, 0.0)
f90, s90 = magnetosonic(3.0, 1.0, math.pi / 2)
assert close(f0, 3.0, 1e-12) and close(s0, 1.0, 1e-12) and close(f90, math.sqrt(10), 1e-12) and s90 == 0.0
# Practical tokamak q constant: 2 pi/(mu0 * 1e6) = 5 when I is in MA (q_cyl = 5 a^2 B/(R I))
assert close(2 * math.pi / (MU0 * 1e6), 5.0, 1e-9)
# Sodium, Stefani et al. 2008: mu0*sigma ~ 10 s/m^2, so Rm = 100 needs U L = 10 m^2/s
print(f"sodium: sigma = {10/MU0:.3g} S/m, eta = {1/10:.2f} m^2/s; Rm=100 needs UL = {100/10:.0f} m^2/s")

# ---------------------------------------------------------------- Q2 perpendicular MHD shock
def perp_shock_X(ma, beta1, g=5 / 3):
    """Density compression X = rho2/rho1 for a perpendicular shock. beta1 = 2 mu0 p1 / B1^2, M_A = u1/vA1."""
    a, b, c = 2 - g, g * (beta1 + 1) + (g - 1) * ma ** 2, -(g + 1) * ma ** 2
    return (-b + math.sqrt(b * b - 4 * a * c)) / (2 * a)


def conservation_residual(ma, beta1, X, g=5 / 3):
    """Normalised rho1 = B1 = mu0 = 1 (vA1 = 1). Returns relative residuals of momentum and energy flux."""
    u1, p1, b1 = ma, beta1 / 2, 1.0
    u2, b2 = u1 / X, b1 * X
    p2 = p1 + u1 ** 2 * (1 - 1 / X) + (b1 ** 2 - b2 ** 2) / 2
    mom = (u1 ** 2 + p1 + b1 ** 2 / 2, X * u2 ** 2 + p2 + b2 ** 2 / 2)
    en = (u1 * (u1 ** 2 / 2 + g / (g - 1) * p1) + u1 * b1 ** 2,
          X * u2 * (u2 ** 2 / 2 + g / (g - 1) * p2 / X) + u2 * b2 ** 2)
    return abs(mom[0] - mom[1]) / mom[0], abs(en[0] - en[1]) / en[0], p2


print("\n== Q2 perpendicular shock ==")
X = perp_shock_X(8.0, 1.0)
rm, re, p2 = conservation_residual(8.0, 1.0, X)
print(f"M_A = 8, beta = 1, gamma = 5/3: X = rho2/rho1 = B2/B1 = {X:.4f}; downstream beta = {2*p2/X**2:.3f}")
assert rm < 1e-12 and re < 1e-12          # identity: the quadratic satisfies the jump conditions exactly
assert close(X, 3.6155, 1e-4)
X2 = perp_shock_X(8.0, 2.1)
print(f"M_A = 8, beta = 2.1 (Wind median beta_e + beta_p): X = {X2:.4f}")
# beta -> infinity at fixed sonic Mach M recovers gas dynamics X = (g+1)M^2/((g-1)M^2 + 2)
g, Ms, bb = 5 / 3, 3.0, 1e9
ma_big = Ms * math.sqrt(g * bb / 2)        # M_A^2 / M_s^2 = g beta / 2
assert close(perp_shock_X(ma_big, bb), (g + 1) * Ms ** 2 / ((g - 1) * Ms ** 2 + 2), 1e-6)
assert close(perp_shock_X(1e5, 1.0), 4.0, 1e-4)   # strong-shock limit (g+1)/(g-1) = 4 (published asymptote)
# fast-mode Mach number must exceed 1: M_A = 1.5, beta = 1 -> M_f = 1.5/sqrt(1 + g beta/2) = 1.15 > 1, X > 1
assert perp_shock_X(1.5, 1.0) > 1.0

# ---------------------------------------------------------------- Q3 tokamak equilibria and limits
def q95_ipb(a, r, b, ip, kappa, delta):
    """ITER Physics Basis edge q (I in MA). Returns (q_star Uckan form without the eps factor, q95)."""
    eps = a / r
    qstar = 5 * a * a * b / (r * ip) * (1 + kappa ** 2 * (1 + 2 * delta ** 2 - 1.2 * delta ** 3)) / 2
    return qstar, qstar * (1.17 - 0.65 * eps) / (1 - eps ** 2) ** 2


def greenwald(ip_ma, a):
    return ip_ma / (math.pi * a * a) * 1e20


def ipb98(ip, b, n19, p, r, kappa_a, eps, m, h=1.0):
    return 0.0562 * h * ip ** 0.93 * b ** 0.15 * n19 ** 0.41 * p ** -0.69 * r ** 1.97 * kappa_a ** 0.78 * eps ** 0.58 * m ** 0.19


print("\n== Q3 equilibria ==")
_, q_iter = q95_ipb(2.0, 6.2, 5.3, 15.0, 1.7, 0.33)
qs_sparc, q_sparc = q95_ipb(0.57, 1.85, 12.2, 8.7, 1.75, 0.45)
print(f"ITER q95 = {q_iter:.3f} (published 3.0); SPARC q* = {qs_sparc:.3f} (published 3.05), IPB q95 = {q_sparc:.3f} (published ~3.6)")
assert close(q_iter, 3.0, 0.01) and close(qs_sparc, 3.05, 0.01) and close(q_sparc, 3.6, 0.01)
ng_iter, ng_sparc = greenwald(15, 2.0), greenwald(8.7, 0.57)
print(f"n_G ITER = {ng_iter:.4g} m^-3 -> 0.85 n_G = {0.85*ng_iter:.4g} (published <n_e> = 1.01e20)")
assert close(0.85 * ng_iter, 1.01e20, 0.01)
print(f"n_G SPARC = {ng_sparc:.4g} m^-3 -> f_G = 0.37 gives {0.37*ng_sparc:.3g} (published <n_e> ~ 3e20)")
assert close(0.37 * ng_sparc, 3e20, 0.06)
for name, bn, ip, a, b in (("ITER", 1.8, 15, 2.0, 5.3), ("SPARC", 1.0, 8.7, 0.57, 12.2)):
    bt = bn * ip / (a * b)
    print(f"{name}: beta_N = {bn} -> beta_t = {bt:.2f} %; Troyon (g = 2.8) ceiling {2.8*ip/(a*b):.2f} %")
# IPB98(y,2): ITER published tau_E = 3.7 s (kappa_a = 1.7, M = 2.5, n19 = 10.1); the net power is not in the source.
tau_iter_at = lambda p: ipb98(15, 5.3, 10.1, p, 6.2, 1.7, 2.0 / 6.2, 2.5)
p_iter = math.exp(math.log(tau_iter_at(1.0) / 3.7) / 0.69)
print(f"ITER: tau_E = 3.7 s needs P_net = {p_iter:.1f} MW (40 MW aux + ~80 MW alpha - radiation)")
assert 40 < p_iter < 120
tau_sparc_at = lambda p: ipb98(8.7, 12.2, 0.37 * ng_sparc / 1e19, p, 1.85, 1.75, 0.57 / 1.85, 2.5)
p_sparc = math.exp(math.log(tau_sparc_at(1.0) / 0.77) / 0.69)
print(f"SPARC: tau_E = 0.77 s needs P_loss = {p_sparc:.1f} MW (P_alpha ~ 28 MW + up to 25 MW ICRH + ~1 MW ohmic)")
assert 15 < p_sparc < 55


# Grad-Shafranov, Solov'ev profiles in Cerfon-Freidberg normalisation: Delta* psi = (1 - A) x^2 + A
def dstar(f, x, y, h=1e-3):
    fxx = (f(x + h, y) - 2 * f(x, y) + f(x - h, y)) / h ** 2
    fyy = (f(x, y + h) - 2 * f(x, y) + f(x, y - h)) / h ** 2
    fx = (f(x + h, y) - f(x - h, y)) / (2 * h)
    return fxx - fx / x + fyy


L = math.log
HOM = [
    lambda x, y: 1.0,
    lambda x, y: x * x,
    lambda x, y: y * y - x * x * L(x),
    lambda x, y: x ** 4 - 4 * x * x * y * y,
    lambda x, y: 2 * y ** 4 - 9 * y * y * x * x + 3 * x ** 4 * L(x) - 12 * x * x * y * y * L(x),
    lambda x, y: x ** 6 - 12 * x ** 4 * y * y + 8 * x * x * y ** 4,
    lambda x, y: (8 * y ** 6 - 140 * y ** 4 * x * x + 75 * y * y * x ** 4 - 15 * x ** 6 * L(x)
                  + 180 * x ** 4 * y * y * L(x) - 120 * x * x * y ** 4 * L(x)),
]


def psi_p(A):
    return lambda x, y: x ** 4 / 8 + A * (x * x * L(x) / 2 - x ** 4 / 8)


for i, f in enumerate(HOM):
    for (x, y) in ((0.8, 0.1), (1.2, -0.3), (1.0, 0.5)):
        assert abs(dstar(f, x, y)) < 1e-3 * max(1.0, abs(f(x, y))), f"psi_{i+1} not homogeneous"
for A in (0.0, -0.155, 0.7):
    for (x, y) in ((0.8, 0.1), (1.2, -0.3)):
        assert abs(dstar(psi_p(A), x, y) - ((1 - A) * x * x + A)) < 1e-5


def d1(f, x, y, ax, h=1e-4):
    return ((f(x + h, y) - f(x - h, y)) if ax == "x" else (f(x, y + h) - f(x, y - h))) / (2 * h)


def d2(f, x, y, ax, h=1e-3):
    return ((f(x + h, y) - 2 * f(x, y) + f(x - h, y)) if ax == "x"
            else (f(x, y + h) - 2 * f(x, y) + f(x, y - h))) / h ** 2


def solve(m, rhs):
    n = len(rhs)
    a = [row[:] + [rhs[i]] for i, row in enumerate(m)]
    for k in range(n):
        p = max(range(k, n), key=lambda r: abs(a[r][k]))
        a[k], a[p] = a[p], a[k]
        for r in range(k + 1, n):
            fct = a[r][k] / a[k][k]
            for c in range(k, n + 1):
                a[r][c] -= fct * a[k][c]
    xs = [0.0] * n
    for k in range(n - 1, -1, -1):
        xs[k] = (a[k][n] - sum(a[k][c] * xs[c] for c in range(k + 1, n))) / a[k][k]
    return xs


def solovev(eps, kappa, delta, A):
    """Up-down symmetric Cerfon-Freidberg equilibrium. Boundary x = 1 + eps cos(t + alpha sin t), y = eps kappa sin t."""
    al = math.asin(delta)
    n1, n2, n3 = -(1 + al) ** 2 / (eps * kappa ** 2), (1 - al) ** 2 / (eps * kappa ** 2), -kappa / (eps * math.cos(al) ** 2)
    out, inn, top = (1 + eps, 0.0), (1 - eps, 0.0), (1 - delta * eps, kappa * eps)
    conds = [
        lambda f: f(*out), lambda f: f(*inn), lambda f: f(*top), lambda f: d1(f, *top, "x"),
        lambda f: d2(f, *out, "y") + n1 * d1(f, *out, "x"),
        lambda f: d2(f, *inn, "y") + n2 * d1(f, *inn, "x"),
        lambda f: d2(f, *top, "x") + n3 * d1(f, *top, "y"),
    ]
    pp = psi_p(A)
    coef = solve([[cond(h) for h in HOM] for cond in conds], [-cond(pp) for cond in conds])
    psi = lambda x, y: pp(x, y) + sum(c * h(x, y) for c, h in zip(coef, HOM))
    return psi, coef, (n1, n2, n3), al


eps, kap, dlt, A = 0.32, 1.7, 0.33, -0.155
psi, coef, (n1, n2, n3), al = solovev(eps, kap, dlt, A)
# curvature constants N1..N3 against the parametric boundary (identity)
bx = lambda t: 1 + eps * math.cos(t + al * math.sin(t))
by = lambda t: eps * kap * math.sin(t)
h = 1e-4
for t0, nref, kind in ((0.0, n1, "xyy"), (math.pi, n2, "xyy"), (math.pi / 2, n3, "yxx")):
    xp, yp = (bx(t0 + h) - bx(t0 - h)) / (2 * h), (by(t0 + h) - by(t0 - h)) / (2 * h)
    xpp, ypp = (bx(t0 + h) - 2 * bx(t0) + bx(t0 - h)) / h ** 2, (by(t0 + h) - 2 * by(t0) + by(t0 - h)) / h ** 2
    curv = (xpp * yp - xp * ypp) / yp ** 3 if kind == "xyy" else (ypp * xp - yp * xpp) / xp ** 3
    assert close(curv, nref, 1e-5)
axis = min((psi(1 - eps + 2 * eps * k / 400, 0.0), 1 - eps + 2 * eps * k / 400) for k in range(1, 400))
edge = max(abs(psi(bx(2 * math.pi * k / 360), by(2 * math.pi * k / 360))) for k in range(360))
print(f"Solov'ev ITER-like (eps .32, kappa 1.7, delta .33, A -.155): psi_axis = {axis[0]:.4e} at x = {axis[1]:.4f}; "
      f"max |psi| on target boundary / |psi_axis| = {edge/abs(axis[0]):.2e}")
print("  coefficients c1..c7 =", " ".join(f"{c:.5e}" for c in coef))
assert edge / abs(axis[0]) < 0.02 and 1 < axis[1] < 1 + eps   # closed surface, Shafranov-shifted axis

# ---------------------------------------------------------------- Q4 instabilities
print("\n== Q4 instabilities ==")
# Kadomtsev criteria on the Bennett profile p = (1 + x^2)^-2 (x = r/r_p). Force balance gives beta = 1/x^2.
for x in (0.5, 1.0, math.sqrt(3), 3.0):
    lhs = 4 * x * x / (1 + x * x)                  # -dln p / dln r
    beta_loc = 1 / x ** 2
    m1_ok = lhs <= 1 / beta_loc + 1e-12            # m >= 1: -dlnp/dlnr <= m^2/beta
    m0_ok = lhs <= 4 * g / (2 + g * beta_loc)      # m = 0: -dlnp/dlnr <= 4 gamma/(2 + gamma beta)
    print(f"x = {x:.3f}: beta = {beta_loc:.3f}, m=1 stable {m1_ok}, m=0 stable {m0_ok}")
    assert m0_ok is False
# m = 1 marginal surface: 4/(1+x^2) = 1 -> x* = sqrt(3), beta* = 1/3 (published, arXiv:2608.00291)
xs = math.sqrt(4 - 1)
assert close(xs, math.sqrt(3), 1e-12) and close(1 / xs ** 2, 1 / 3, 1e-12)


def rt_lambda_c(b, theta, rho_h, rho_l, g_acc):
    """Critical (shortest unstable) wavelength, SI. theta = angle between B and the wavevector k, both in the
    interface plane (theta = 90 deg: no tension, no cut-off), as in Khomenko et al."""
    return 4 * math.pi * (b * math.cos(theta)) ** 2 / (MU0 * g_acc * (rho_h - rho_l))


# Khomenko et al. 2014 (arXiv:1310.7016): B = 10 G, rho = 2.9e-10 / 3.7e-12 kg m^-3, theta = 89 deg -> ~38 km, 88 deg -> ~155 km.
# g_sun = 274 m s^-2 (not stated in that paper; standard solar surface gravity).
lc89 = rt_lambda_c(1e-3, math.radians(89), 2.9e-10, 3.7e-12, 274.0)
lc88 = rt_lambda_c(1e-3, math.radians(88), 2.9e-10, 3.7e-12, 274.0)
print(f"magnetic RT lambda_c = {lc89/1e3:.1f} km (89 deg), {lc88/1e3:.1f} km (88 deg); published 38 and 155 km")
assert close(lc89, 38e3, 0.03) and close(lc88, 155e3, 0.01)
# fastest-growing wavelength = 2 lambda_c (identity of gamma^2 = A g k - 2 k^2 B^2 cos^2/(mu0 (rho_h+rho_l)))
Agk = lambda k, bc=1e-3 * math.cos(math.radians(89)): ((2.9e-10 - 3.7e-12) / (2.9e-10 + 3.7e-12) * 274 * k
                                                       - 2 * k * k * bc * bc / (MU0 * (2.9e-10 + 3.7e-12)))
kc = 2 * math.pi / lc89
kgrid = [kc * j / 2000 for j in range(1, 2000)]
kmax = max(kgrid, key=Agk)
assert close(2 * math.pi / kmax, 2 * lc89, 2e-3)

fkr = (math.gamma(0.25) / (2 * math.pi * math.gamma(0.75))) ** 0.8
print(f"FKR constant-psi coefficient [Gamma(1/4)/(2 pi Gamma(3/4))]^(4/5) = {fkr:.4f} (quoted ~0.55)")
assert close(fkr, 0.55, 0.01)
for S in (1e6, 1e8):
    print(f"  S = {S:.0e}, Delta'a = 1: gamma tau = {fkr * S ** -0.6:.3e}, layer width ~ S^-2/5 = {S ** -0.4:.2e}, "
          f"constant-psi valid while Delta'a << S^1/3 = {S ** (1/3):.0f}")

# ---------------------------------------------------------------- Q5 Hartmann flow and duct pressure drop
def hartmann_profile(ha, eta):
    """u(eta)/u_mean between plates at eta = y/a = +-1, field normal to plates."""
    return ha * (math.cosh(ha) - math.cosh(ha * eta)) / (ha * math.cosh(ha) - math.sinh(ha))


def hartmann_dpdx(ha, insulating=True):
    """-dp/dx / (sigma u_mean B^2). Insulating walls: net current zero; perfectly conducting walls: E = 0."""
    t = math.tanh(ha)
    return t / (ha - t) if insulating else ha / (ha - t)


print("\n== Q5 Hartmann and duct flow ==")
for ha in (1e-2, 10.0, 100.0):
    mean = sum(hartmann_profile(ha, -1 + (k + 0.5) / 5000) for k in range(10000)) / 10000
    assert close(mean, 1.0, 1e-6)
    print(f"Ha = {ha:g}: u_centre/u_mean = {hartmann_profile(ha, 0):.4f}, "
          f"-dp/dx/(sigma u B^2) insulating {hartmann_dpdx(ha):.4g}, conducting {hartmann_dpdx(ha, False):.4g}")
# limits (identity): Ha -> 0 gives Poiseuille u_c/u_mean = 1.5 and -dp/dx = 3 mu u/a^2 = 3 sigma B^2 u / Ha^2
assert close(hartmann_profile(1e-3, 0), 1.5, 1e-5) and close(hartmann_dpdx(1e-3) * 1e-6, 3.0, 1e-4)
assert close(hartmann_dpdx(1e4), 1 / (1e4 - 1), 1e-9)      # insulating: wall friction only, ~ 1/Ha
# Circular pipe, thin/finite conducting wall, Ha >> 1: -dp/dx = sigma u B^2 c/(1+c).
# KIT (Buhler et al.): c = 0.0727 -> published nondimensional gradient 0.0678
print(f"pipe c = 0.0727: c/(1+c) = {0.0727/1.0727:.4f} (published 0.0678)")
assert close(0.0727 / 1.0727, 0.0678, 1e-3)

# ---------------------------------------------------------------- Q6 dynamos
print("\n== Q6 dynamos ==")


def alpha_omega_re(D, alpha):
    """Re of the growing branch of Eq. 6.38 (Brandenburg & Subramanian 2005), eta_T = k = 1, S = D/alpha."""
    return (-1 + cmath.sqrt(alpha ** 2 - 1j * alpha * (D / alpha))).real


lo, hi = 0.1, 10.0
for _ in range(200):
    mid = 0.5 * (lo + hi)
    lo, hi = (mid, hi) if alpha_omega_re(mid, 1e-4) < 0 else (lo, mid)
Dc = 0.5 * (lo + hi)
w_c = abs(cmath.sqrt(-1j * Dc).imag)
print(f"alpha-Omega plane wave: D_crit = {Dc:.5f} (analytic 2), cycle frequency at onset = {w_c:.4f} eta k^2 (analytic 1)")
assert close(Dc, 2.0, 1e-3) and close(w_c, 1.0, 1e-3)
# Christensen et al. 2009 in dimensional form B^2 = 1.2 mu0 f_ohm rho^(1/3) (F q)^(2/3) (Christensen 2009 lecture).
# Inverse, illustrative only: flux F q needed for a 1-4 mT rms field in a rho = 1.1e4 kg m^-3 core with f_ohm = 1.
for b_core in (1e-3, 4e-3):
    fq = (b_core ** 2 / (1.2 * MU0 * 1.1e4 ** (1 / 3))) ** 1.5
    print(f"  B_rms = {b_core*1e3:.0f} mT needs F q = {fq*1e3:.2f} mW m^-2")
# Earth illustrative numbers (inputs NOT from a retrieved source; see notes): sigma 1e6 S/m, U 5e-4 m/s, D 2.26e6 m
sig, Uc, Dsh, rc, rho_c, Om = 1e6, 5e-4, 2.26e6, 3.48e6, 1.1e4, 7.292e-5
print(f"  Earth: Rm = {MU0*sig*Uc*Dsh:.0f} (lecture: ~1e3), dipole free decay r^2 mu0 sigma/pi^2 = "
      f"{MU0*sig*rc**2/math.pi**2/3.156e7/1e3:.0f} kyr, Elsasser sigma B^2/(2 rho Omega) at 2 mT = "
      f"{sig*(2e-3)**2/(2*rho_c*Om):.2f} (lecture: 0.3-5)")
assert 300 < MU0 * sig * Uc * Dsh < 3000 and 0.3 < sig * (2e-3) ** 2 / (2 * rho_c * Om) < 5

# ---------------------------------------------------------------- Q7 Hall scales and turbulence
print("\n== Q7 Hall MHD and turbulence ==")
# PlasmaPy inertial_length docstring: n = 5 m^-3, He+ -> 2.03e8 m; e- -> 2.38e6 m
d_he = inertial_length(5.0, 4.00260325413 * U - ME)
d_e = inertial_length(5.0, ME)
print(f"d(He+, 5 m^-3) = {d_he:.4g} m, d_e = {d_e:.4g} m")
assert close(d_he, 2.03e8, 3e-3) and close(d_e, 2.38e6, 3e-3)
# identity d_i = v_A / Omega_i
n_t, b_t = 5e6, 5e-9
assert close(inertial_length(n_t, MP), alfven(b_t, n_t * MP) / (E * b_t / MP), 1e-12)
# Chen et al. 2014 (GRL 41, 8081), Table 1 means; break at the larger of d_i and rho_i.
# Taylor hypothesis f_b ~ v_sw / (2 pi l_break). Published f_b = 0.58 +- 0.19 Hz (low beta), 0.223 +- 0.068 Hz (high beta)
for label, b_nt, n_cc, v_kms, ti, f_pub, sd, beta_pub in (("low beta", 13.3, 3.9, 460, 1.64, 0.58, 0.19, 0.0122),
                                                       ("high beta", 1.14, 17.1, 389, 4.8, 0.223, 0.068, 23.2)):
    di, ri = inertial_length(n_cc * 1e6, MP), gyroradius(ti, b_nt * 1e-9)
    fb = v_kms * 1e3 / (2 * math.pi * max(di, ri))
    bi = beta(n_cc * 1e6, ti, b_nt * 1e-9)
    print(f"{label}: d_i = {di/1e3:.0f} km, rho_i = {ri/1e3:.0f} km, predicted f_b = {fb:.3f} Hz "
          f"(published {f_pub} +- {sd}); beta_i from means = {bi:.4g} (published mean {beta_pub})")
    assert abs(fb - f_pub) < sd
# tokamak-core Hall scale, illustrative: ITER n = 1.01e20, D-T (2.5 u)
print(f"ITER core: d_i = {inertial_length(1.01e20, 2.5*U)*100:.1f} cm, v_A = {alfven(5.3, 1.01e20*2.5*U)/1e6:.2f} Mm/s")
# Goldreich-Sridhar: k_par = k_perp^(2/3) L^(-1/3) -> anisotropy k_perp/k_par = (k_perp L)^(1/3)
print(f"GS95 anisotropy at k_perp L = 1e6: k_perp/k_par = {(1e6) ** (1/3):.0f}")
print("\nall checks passed")
