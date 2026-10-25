"""Gap 5: relativistic Love number k2 = TOV + Riccati equation for y = rH'/H (Postnikov et al. 2010
Eq. 7-9, Hinderer 2008 Eqs. 15/23 as corrected in the erratum), surface jump y -> y - 4 pi R^3 e_-/M
(Damour & Nagar 2009 Eq. 'minus_3'; Hinderer et al. 2010 after Eq. 5).

Units G = c = 1, lengths in km for EOS work. Checks:
  * Newtonian n = 1: k2 -> (15 - pi^2)/(2 pi^2) as C -> 0 (Hinderer 2008; Postnikov Eq. n=1)
  * incompressible: k2 -> 3/4 as C -> 0; k2(4/9) = 4096/(10935 (308 - 81 ln 3)) = 0.0017103 (DN09)
  * gamma = 2 mu-polytrope vs DN09 Table I fit; Hinderer 2008 erratum Table 1 (n = 1, 0.5, ...)
  * SLy and AP4 piecewise polytropes (Read et al. 2009 Tables II-III) vs Read's R_1.4, M_max and
    Hinderer et al. 2010 Table I (k2, R at 1.4 Msun for tabulated SLy)
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

PI = np.pi


def k2_of(C, y):
    """Hinderer 2008 erratum Eq. (23)."""
    num = 8 * C**5 / 5 * (1 - 2 * C) ** 2 * (2 + 2 * C * (y - 1) - y)
    den = (2 * C * (6 - 3 * y + 3 * C * (5 * y - 8)) + 4 * C**3 * (13 - 11 * y + C * (3 * y - 2) + 2 * C**2 * (1 + y))
           + 3 * (1 - 2 * C) ** 2 * (2 - y + 2 * C * (y - 1)) * np.log(1 - 2 * C))
    return num / den


def k2_series(C, y):
    """Postnikov Eq. (10) small-C series, as printed (its 4th term is labelled beta^4; see notes)."""
    Y = y + 3
    return (1 - 2 * C) ** 2 / 2 * ((2 - y) / Y + (y * y - 6 * y - 6) / Y**2 * C + (y**3 + 34 * y**2 - 8 * y + 12) / (7 * Y**3) * C**2)


def star(eos_e, eos_cs2inv, pc, r0=1e-6, jump=0.0):
    """Integrate m, p, y outward in r (G=c=1). eos_e(p) -> energy density; eos_cs2inv(p) -> de/dp.
    Returns M, R, y_R (with surface jump if the surface density is finite)."""
    def rhs(r, u):
        m, p, y = u
        p = max(p, 0.0)
        e = eos_e(p)
        el = 1 / (1 - 2 * m / r)
        dp = -(e + p) * (m + 4 * PI * r**3 * p) / (r * (r - 2 * m))
        nup = 2 * el * (m + 4 * PI * p * r**3) / r**2
        Q = 4 * PI * el * (5 * e + 9 * p + (e + p) * eos_cs2inv(p)) - 6 * el / r**2 - nup**2
        dy = -(y * y + y * el * (1 + 4 * PI * r * r * (p - e)) + r * r * Q) / r
        return [4 * PI * r * r * e, dp, dy]

    surf = lambda r, u: u[1] - pc * 1e-14
    surf.terminal, surf.direction = True, -1
    e0 = eos_e(pc)
    u0 = [4 * PI / 3 * e0 * r0**3, pc, 2.0]
    sol = solve_ivp(rhs, [r0, 1e3], u0, events=surf, rtol=1e-11, atol=1e-30, method="DOP853", max_step=0.05)
    m, p, y = sol.y[:, -1]
    R = sol.t[-1]
    if jump:  # finite surface energy density e_- = jump
        y -= 4 * PI * R**3 * jump / m
    return m, R, y


# ---------------- 1. incompressible star (Schwarzschild interior), r in units where e0 = 1 ------------
def incompressible(C):
    """Exact TOV pressure for uniform density, y integrated in r, then y_out = y_in - 3."""
    e0 = 1.0
    R = np.sqrt(3 * C / (4 * PI * e0))
    M = C * R
    def p_of(r):
        a = np.sqrt(1 - 2 * M * r * r / R**3); b = np.sqrt(1 - 2 * M / R)
        return e0 * (a - b) / (3 * b - a)
    def rhs(r, u):
        y = u[0]; m = 4 * PI / 3 * e0 * r**3; p = p_of(r); el = 1 / (1 - 2 * m / r)
        nup = 2 * el * (m + 4 * PI * p * r**3) / r**2
        Q = 4 * PI * el * (5 * e0 + 9 * p) - 6 * el / r**2 - nup**2   # (e+p)de/dp = 0 inside
        return [-(y * y + y * el * (1 + 4 * PI * r * r * (p - e0)) + r * r * Q) / r]
    s = solve_ivp(rhs, [1e-8 * R, R], [2.0], rtol=1e-12, atol=1e-14, method="DOP853")
    return k2_of(C, s.y[0, -1] - 3.0)


print("incompressible: k2(C->0) =", incompressible(1e-4), "(3/4); k2(0.44) =", incompressible(0.44),
      "; DN09 k2(4/9) =", 4096 / (10935 * (308 - 81 * np.log(3))))
for C in [0.05, 0.1, 0.2, 0.3, 0.4, 0.444]:
    print(f"   incompressible C={C}: k2={incompressible(C):.6f}")

# ---------------- 2. polytropes ----------------
def polytrope(kind, Gam, C_target):
    """kind 'mu': p = K mu^G, e = mu + p/(G-1) (DN09 mu-polytrope); kind 'e': p = K e^G. K = 1."""
    K = 1.0
    if kind == "mu":
        mu = lambda p: (p / K) ** (1 / Gam)
        e_of = lambda p: mu(p) + p / (Gam - 1)
        dedp = lambda p: (mu(p) / (Gam * p) + 1 / (Gam - 1)) if p > 0 else 0.0  # de/dp = de/dmu / (dp/dmu)
    else:
        e_of = lambda p: (p / K) ** (1 / Gam)
        dedp = lambda p: (e_of(p) / (Gam * p)) if p > 0 else 0.0
    def comp(pc):
        M, R, y = star(e_of, dedp, pc)
        return M / R, k2_of(M / R, y)
    pc = brentq(lambda lp: comp(np.exp(lp))[0] - C_target, np.log(1e-8), np.log(0.5), xtol=1e-12)
    return comp(np.exp(pc))


k2N = (15 - PI**2) / (2 * PI**2)
print("\nNewtonian n=1 limit k2N =", k2N)
a = [0.9991, -2.9287, -1.1373, 14.0013, -50.9711]   # DN09 Table I, l = 2
for C in [1e-3, 0.05, 0.1, 0.15, 0.2]:
    Cm, k_mu = polytrope("mu", 2.0, C)
    Ce, k_e = polytrope("e", 2.0, C)
    fit = k2N * (1 - 2 * C) ** 2 * sum(an * C**i for i, an in enumerate(a))
    print(f"  Gamma=2 C={C}: mu-polytrope k2={k_mu:.5f} (DN09 fit {fit:.5f}); e-polytrope k2={k_e:.5f}")
print("  Hinderer 2008 erratum Table 1, n=1: C=0.1 0.122, C=0.15 0.0776, C=0.2 0.0459")
for n, rows in [(0.5, [(0.1, 0.251), (0.15, 0.173), (0.2, 0.095)]), (0.3, [(0.1, 0.294), (0.2, 0.119)])]:
    G = 1 + 1 / n
    for C, ref in rows:
        print(f"  n={n} C={C}: mu {polytrope('mu', G, C)[1]:.4f}  e {polytrope('e', G, C)[1]:.4f}  Hinderer {ref}")

# ---------------- 3. Read et al. 2009 piecewise polytropes ----------------
c_cgs, G_cgs, Msun = 2.99792458e10, 6.67430e-8, 1.98847e33
geo = G_cgs / c_cgs**2 * 1e10           # g/cm^3 -> km^-2
crust = [(6.80110e-09, 1.58425, 2.44034e+07), (1.06186e-06, 1.28733, 3.78358e+11), (5.32697e+01, 0.62223, 2.62780e+12), (3.99874e-08, 1.35692, None)]


def pwp(logp1, G1, G2, G3):
    rho1, rho2 = 10**14.7, 10**15.0
    p1 = 10**logp1 / c_cgs**2              # p/c^2 in g/cm^3 (Read: K, p in c = 1 form)
    K1 = p1 / rho1**G1; K2 = p1 / rho1**G2; K3 = K2 * rho2**(G2 - G3)
    Kc, Gc, _ = crust[3]
    rho0 = (Kc / K1) ** (1 / (G1 - Gc))
    pieces = crust[:3] + [(Kc, Gc, rho0), (K1, G1, rho1), (K2, G2, rho2), (K3, G3, np.inf)]
    # energy density continuity: e = (1+a_i) rho + K rho^G/(G-1)
    lo, avals, e_lo = 0.0, [], 0.0
    for i, (K, G, hi) in enumerate(pieces):
        a = 0.0 if i == 0 else e_lo / lo - 1 - K / (G - 1) * lo ** (G - 1)
        avals.append(a)
        e_lo = (1 + a) * hi + K / (G - 1) * hi**G if np.isfinite(hi) else None
        lo = hi
    pb = [K * hi**G for (K, G, hi) in pieces[:-1]]   # pressure at the upper bound of each piece

    def piece(p):
        return np.searchsorted(pb, p)

    def e_of(p):   # p and e in km^-2
        if p <= 0:
            return 0.0
        pg = p / geo
        i = piece(pg); K, G, _ = pieces[i]
        rho = (pg / K) ** (1 / G)
        return ((1 + avals[i]) * rho + K / (G - 1) * rho**G) * geo

    def dedp(p):
        if p <= 0:
            return 0.0
        pg = p / geo
        i = piece(pg); K, G, _ = pieces[i]
        rho = (pg / K) ** (1 / G)
        return ((1 + avals[i]) + K * G / (G - 1) * rho ** (G - 1)) / (K * G * rho ** (G - 1))

    cont = [(pieces[i][0] * pieces[i][2] ** pieces[i][1]) / (pieces[i + 1][0] * pieces[i][2] ** pieces[i + 1][1]) - 1 for i in range(len(pieces) - 1)]
    return e_of, dedp, cont


def ns(eos, mass):
    e_of, dedp, _ = eos
    def m_of(lpc):
        M, R, y = star(e_of, dedp, np.exp(lpc))
        return M, R, y
    def mass_msun(lpc):
        return m_of(lpc)[0] / (G_cgs * Msun / c_cgs**2 / 1e5)
    lps = np.linspace(np.log(1e-5), np.log(3e-3), 60)
    ms = [mass_msun(l) for l in lps]
    imax = int(np.argmax(ms))
    from scipy.optimize import minimize_scalar
    r = minimize_scalar(lambda l: -mass_msun(l), bracket=(lps[imax - 1], lps[imax], lps[imax + 1]), tol=1e-10)
    Mmax = -r.fun
    l14 = brentq(lambda l: mass_msun(l) - mass, lps[0], lps[imax], xtol=1e-13)
    M, R, y = m_of(l14)
    C = M / R
    k2 = k2_of(C, y)
    return Mmax, R, C, k2, 2 / 3 * k2 / C**5


print("\nRead et al. 2009 piecewise polytropes, 1.4 Msun:")
for name, prm, read_R, read_Rpct, read_M, read_Mpct, hind in [
        ("SLy", (34.384, 3.005, 2.988, 2.851), 11.736, -0.21, 2.049, 0.02, (11.74, 0.176, 0.0763)),
        ("AP4", (34.269, 2.830, 3.445, 3.348), 11.428, -0.90, 2.213, -1.08, None)]:
    eos = pwp(*prm)
    print(f"  {name}: pressure continuity at boundaries (rel): {max(abs(c) for c in eos[2]):.2e}")
    Mmax, R, C, k2, Lam = ns(eos, 1.4)
    print(f"  {name}: M_max = {Mmax:.4f} (Read fit: {read_M*(1+read_Mpct/100):.4f}); R_1.4 = {R:.3f} km (Read fit: {read_R*(1+read_Rpct/100):.3f}, table {read_R})"
          f"; C = {C:.4f}; k2 = {k2:.5f}; Lambda_1.4 = {Lam:.1f}" + (f"   [Hinderer+2010 tabulated: R {hind[0]}, C {hind[1]}, k2 {hind[2]}, Lambda {2/3*hind[2]/hind[1]**5:.0f}]" if hind else ""))
