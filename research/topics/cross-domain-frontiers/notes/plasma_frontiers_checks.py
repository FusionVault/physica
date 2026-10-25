"""Recompute the published numbers behind plasma_frontiers_and_kinetic_algorithms.md.

Standard library only. Run: python3 -I plasma_frontiers_checks.py
Each check prints the recomputed value beside the published one and asserts a tolerance.
Constants are CODATA 2018, the values the `metrology` crate uses.
"""
import cmath
from decimal import Decimal, getcontext
from math import pi, sqrt, log, exp, erfc, cos, tan

E = 1.602176634e-19
ME = 9.1093837015e-31
MP = 1.67262192369e-27
AMU = 1.66053906660e-27
EPS0 = 8.8541878128e-12
KB = 1.380649e-23
HBAR = 1.054571817e-34
KEV = 1e3 * E
M_D, M_T = 2.013553212745 * AMU, 3.01550071621 * AMU  # nuclear masses, as in plasmoid::fusion


def check(name, got, want, rtol=None, atol=None):
    tol = atol if atol is not None else rtol * abs(want)
    ok = abs(got - want) <= tol
    print(f"{'ok ' if ok else 'BAD'} {name}: {got:.10g} (published {want:.10g})")
    assert ok, name


# ---------------------------------------------------------------- 1. Faddeeva w(z)
def w_weideman(z, n=32):
    """Weideman (1994) SIAM J. Numer. Anal. 31, 1497, Eq. (38)-(39): rational series in
    Z = (L + iz)/(L - iz), L = sqrt(N/sqrt 2), coefficients a_n from a 2M-point cosine
    transform of f(t) = exp(-t^2)(L^2 + t^2) at t = L tan(theta/2), M = 2N. Im z >= 0."""
    m = 2 * n
    L = sqrt(n / sqrt(2.0))
    ks = range(-m + 1, m)
    f = {k: exp(-(L * tan(k * pi / (2 * m))) ** 2) * (L * L + (L * tan(k * pi / (2 * m))) ** 2) for k in ks}
    a = [sum(f[k] * cos(j * k * pi / m) for k in ks) / (2 * m) for j in range(1, n + 1)]
    zz = (L + 1j * z) / (L - 1j * z)
    p = 0j
    for c in reversed(a):  # Horner: sum a_j Z^(j-1)
        p = p * zz + c
    return 2 * p / (L - 1j * z) ** 2 + 1 / (sqrt(pi) * (L - 1j * z))


def w_cf(z, terms=200):
    """Laplace continued fraction w = (i/sqrt pi)/(z - (1/2)/(z - 1/(z - (3/2)/(z - ...)))),
    accurate for |z| >> 1 with Im z > 0."""
    t = 0j
    for k in range(terms, 0, -1):
        t = (k / 2) / (z - t)
    return 1j / (sqrt(pi) * (z - t))


getcontext().prec = 90


def _dpi():
    def atan_inv(x):  # atan(1/x)
        x = Decimal(x)
        term, s, k, x2 = 1 / x, 1 / x, 1, x * x
        while True:
            term = -term / x2
            k += 2
            d = term / k
            if abs(d) < Decimal(10) ** -(getcontext().prec + 2):
                return s
            s += d
    return 16 * atan_inv(5) - 4 * atan_inv(239)


DPI = _dpi()


def w_series(x, y):
    """Reference: w(z) = sum (iz)^n / Gamma(n/2 + 1) in 90-digit decimal arithmetic.
    x, y are strings; good for |z| <~ 7."""
    x, y = Decimal(x), Decimal(y)
    pr, pi_ = -y, x  # iz = -y + ix
    tr, ti = Decimal(1), Decimal(0)  # (iz)^n
    g = [Decimal(1), DPI.sqrt() / 2]  # Gamma(1), Gamma(3/2)
    sr, si = Decimal(1), Decimal(0)
    n = 0
    while True:
        n += 1
        tr, ti = tr * pr - ti * pi_, tr * pi_ + ti * pr
        if n >= 2:
            g.append(Decimal(n) / 2 * g[n - 2])
        dr, di = tr / g[n], ti / g[n]
        sr += dr
        si += di
        if n > 20 and abs(dr) + abs(di) < Decimal(10) ** -60:
            return complex(float(sr), float(si))


# Zaghloul & Ali (2011) ACM TOMS 38(2) 15, Table 2, Mathematica erfi-based reference values.
ZA = [
    ("6.3e-2", "1e-20", 9.960388660702479e-01, 7.090008726353683e-02),
    ("6.3e-2", "1e-2", 9.849424862549036e-01, 6.965909657459020e-02),
    ("6.3e-2", "10", 5.613881832823887e-02, 3.502232333332985e-04),
    ("6.3", "1e-20", 5.792460778844102e-18, 9.072765968412736e-02),
    ("6.3", "1e-6", 1.478934493028413e-08, 9.072765968412492e-02),
    ("6.3", "1e-2", 1.478930389133942e-04, 9.072741516349275e-02),
    ("6.3", "10", 4.040671157393860e-02, 2.527577277549421e-02),
    ("1", "1e-20", 3.678794411714423e-01, 6.071577058413937e-01),
    ("5.5", "1e-14", 7.307386729528773e-14, 1.043674364367812e-01),
]
print("== Faddeeva w(z) against Zaghloul & Ali 2011 Table 2")
for xs, ys, v, l in ZA:
    ref = w_series(xs, ys)
    check(f"series Re w({xs}+{ys}i)", ref.real, v, rtol=1e-13)
    check(f"series Im w({xs}+{ys}i)", ref.imag, l, rtol=1e-13)
    wd = w_weideman(complex(float(xs), float(ys)), 32)
    # Weideman N=32 is accurate in absolute terms (~1e-14 of |w|) but not where Re w << |w|.
    check(f"Weideman32 |w - ref| at {xs}+{ys}i", abs(wd - complex(v, l)), 0.0, atol=2e-13)
# Far field: continued fraction against the table.
for x, y, v, l in [(630, 10, 1.421137820009847e-05, 8.953145713915760e-04),
                   (630, 1e-20, 1.421495882582394e-26, 8.955401496757104e-04),
                   (6.3e-2, 200, 2.820912377324508e-03, 8.885651855627418e-07)]:
    r = w_cf(complex(x, y))
    check(f"CF Re w({x}+{y}i)", r.real, v, rtol=1e-12)
    check(f"CF Im w({x}+{y}i)", r.imag, l, rtol=1e-12)
# Weideman fails in relative terms near the real axis where Re w = exp(-x^2) is tiny:
wd = w_weideman(complex(6.3, 1e-20), 32)
print(f"   Weideman32 Re w(6.3+1e-20 i) = {wd.real:.3e} against 5.792e-18: relative error "
      f"{abs(wd.real / 5.792460778844102e-18 - 1):.1e}")
# Closed forms: w(iy) = erfcx(y); Re w(x) = exp(-x^2); Im w(x) = 2 F(x)/sqrt(pi), F(1) = 0.5380795069127684
check("w(2i) = exp(4) erfc(2)", w_weideman(2j, 32).real, exp(4) * erfc(2), rtol=1e-13)
check("Dawson F(1) = (sqrt pi/2) Im w(1)", sqrt(pi) / 2 * w_series("1", "0").imag, 0.5380795069127684, rtol=1e-14)


def w(z, n=40):
    """Faddeeva w on the whole plane: w(z) = 2 exp(-z^2) - w(-z) below the real axis."""
    return w_weideman(z, n) if z.imag >= 0 else 2 * cmath.exp(-z * z) - w_weideman(-z, n)


def Zf(zeta):
    """Plasma dispersion function Z(zeta) = i sqrt(pi) w(zeta) (Fried & Conte 1961), analytic
    continuation of the Landau contour included."""
    return 1j * sqrt(pi) * w(zeta)


print("== Plasma dispersion function against PlasmaPy docstring values")
check("Z(0) imag = sqrt(pi)", Zf(0j).imag, 1.77245385, rtol=1e-8)
for z, zr in [(1 + 1j, -0.36905845 + 0.54014504j), (0.3 + 0j, -0.56526333 + 1.61990085j),
              (0.7 + 2.3j, -0.09995023 + 0.37685142j)]:
    check(f"Z({z}) real", Zf(z).real, zr.real, atol=1e-8)  # docstring digits are truncated
    check(f"Z({z}) imag", Zf(z).imag, zr.imag, atol=1e-8)
# Z'(zeta) = -2[1 + zeta Z(zeta)] against a centred difference
h = 1e-5
z0 = 0.8 - 0.3j
check("Z' identity", abs((Zf(z0 + h) - Zf(z0 - h)) / (2 * h) - (-2 * (1 + z0 * Zf(z0)))), 0.0, atol=1e-8)

# ---------------------------------------------------------------- 2. Landau damping
def langmuir_root(K, guess):
    """Root of 1 + (1/K^2)[1 + xi Z(xi)] = 0, xi = omega/(sqrt2 K), omega in units of omega_pe,
    K = k lambda_De, v_t = sqrt(T/m)."""
    om = guess
    for _ in range(60):
        xi = om / (sqrt(2) * K)
        Z = Zf(xi)
        D = 1 + (1 + xi * Z) / K**2
        dZ = -2 * (1 + xi * Z)
        dD = (Z + xi * dZ) / K**2 / (sqrt(2) * K)
        step = D / dD
        om -= step
        if abs(step) < 1e-14:
            return om
    raise RuntimeError("no convergence")


print("== Langmuir waves: exact kinetic root against the Vlasov benchmark k lambda_D = 0.5")
om = langmuir_root(0.5, 1.4 - 0.15j)
check("omega_r(k lambda_D = 0.5)", om.real, 1.4156, atol=1e-4)
check("gamma(k lambda_D = 0.5)", -om.imag, 0.1533, atol=1e-4)
for K in (0.1, 0.2, 0.3, 0.4, 0.5):
    om = langmuir_root(K, complex(sqrt(1 + 3 * K * K), -0.01))
    bg = sqrt(1 + 3 * K * K)
    weak = sqrt(pi / 8) / K**3 * exp(-1 / (2 * K * K) - 1.5)
    print(f"   K={K}: omega = {om.real:.6f} - {-om.imag:.6e} i | Bohm-Gross {bg:.6f} | "
          f"weak-damping formula {weak:.6e}")


# Weak damping is below Weideman's ~3e-14 absolute error: at K = 0.1 the root's imaginary part is
# noise. Use the exact Im Z(xi) = sqrt(pi) exp(-xi^2) on the real axis instead (perturbative rate).
K = 0.1
om_r = langmuir_root(K, complex(sqrt(1 + 3 * K * K), -0.01)).real
xi = om_r / (sqrt(2) * K)
im_D = sqrt(pi) * xi * exp(-xi * xi) / K**2
d = 1e-6
dre = ((1 + ((om_r + d) / (sqrt(2) * K)) * Zf(complex((om_r + d) / (sqrt(2) * K), 0)).real)
       - (1 + ((om_r - d) / (sqrt(2) * K)) * Zf(complex((om_r - d) / (sqrt(2) * K), 0)).real)) / (2 * d) / K**2
check("K = 0.1 perturbative gamma against the weak-damping formula", im_D / dre,
      sqrt(pi / 8) / K**3 * exp(-1 / (2 * K * K) - 1.5), rtol=0.05)


def ion_acoustic_root(K, tau, mu, guess):
    """Electrons + singly charged ions, n_i = n_e: 1 + (1/K^2)[1 + xi_e Z(xi_e)] +
    (tau/K^2)[1 + xi_i Z(xi_i)] = 0, Omega = omega/omega_pi, xi_e = Omega/(sqrt2 K sqrt mu),
    xi_i = Omega sqrt(tau)/(sqrt2 K), tau = Te/Ti, mu = mi/me."""
    om = guess
    for _ in range(100):
        def D(o):
            xe, xi = o / (sqrt(2) * K * sqrt(mu)), o * sqrt(tau) / (sqrt(2) * K)
            return 1 + (1 + xe * Zf(xe)) / K**2 + tau * (1 + xi * Zf(xi)) / K**2
        d = 1e-7 * max(abs(om), 1e-3)
        step = D(om) / ((D(om + d) - D(om - d)) / (2 * d))
        om -= step
        if abs(step) < 1e-13 * abs(om):
            return om
    raise RuntimeError("no convergence")


print("== Ion-acoustic waves, hydrogen, K = k lambda_De = 0.1")
mu = MP / ME
for tau in (3.0, 5.0, 10.0, 20.0):
    fluid = 0.1 / sqrt(1 + 0.01) * sqrt(1 + 3 / tau)
    om = ion_acoustic_root(0.1, tau, mu, complex(fluid, -0.01 * fluid))
    approx = sqrt(pi / 8) * (sqrt(1 / mu) + tau**1.5 * exp(-tau / 2 - 1.5))
    print(f"   Te/Ti={tau:>4}: omega/omega_pi = {om.real:.6f}, -gamma/omega = {-om.imag / om.real:.4f} "
          f"(textbook small-K estimate {approx:.4f})")

# ---------------------------------------------------------------- 3. Debye, Coulomb log, Spitzer
def debye(n, t_ev):
    return sqrt(EPS0 * t_ev * E / (n * E * E))


print("== Coulomb logarithm and resistivity against PlasmaPy (n_e = 1e19 m^-3, T = 1e6 K, e-p)")
T = KB * 1e6
lam_d = sqrt(EPS0 * T / (1e19 * E * E))
mu_ep = ME * MP / (ME + MP)


def lnlam_classical(V):
    b_perp = E * E / (4 * pi * EPS0 * mu_ep * V * V)
    b_db = HBAR / (2 * mu_ep * V)
    return log(lam_d / max(b_perp, b_db)), b_perp


V_th = sqrt(2 * T / mu_ep)
ln_th, b_perp = lnlam_classical(V_th)
check("PlasmaPy classical lnL (thermal V)", ln_th, 14.545527556, rtol=1e-8)
check("PlasmaPy classical lnL (V = 1e6 m/s)", lnlam_classical(1e6)[0], 11.363478378, rtol=1e-8)
nu = 1e19 * pi * (2 * b_perp) ** 2 * V_th * ln_th  # PlasmaPy's collision_frequency convention
check("PlasmaPy Spitzer_resistivity = mu nu/(n e^2)", mu_ep * nu / (1e19 * E * E), 2.4915e-6, rtol=1e-3)  # 3e-4 residual untraced
t_ev = T / E
nrl_ei = 24 - log(sqrt(1e13) / t_ev)  # NRL, T_e > 10 Z^2 eV, n in cm^-3
print(f"   NRL electron-ion lnL at the same point = {nrl_ei:.4f} (PlasmaPy 14.5455)")
eta_perp = 1.03e-4 * nrl_ei * t_ev**-1.5
print(f"   NRL eta_perp = {eta_perp:.4e} Ohm m, eta_par = eta_perp/1.96 = {eta_perp / 1.96:.4e}; "
      f"PlasmaPy's 2.49e-6 is {2.4915e-6 / eta_perp:.2f} x eta_perp")
# NRL's 1.03e-4 follows from Braginskii tau_e = 3.44e5 T^1.5/(n lnL) (n cm^-3): eta_perp = me/(n e^2 tau_e)
tau_e_coef = 6 * sqrt(2) * pi**1.5 * EPS0**2 * sqrt(ME) * E**1.5 / E**4  # SI, T in eV, n in m^-3
check("Braginskii tau_e coefficient (cm^-3 units)", tau_e_coef / 1e6, 3.44e5, rtol=2e-3)
check("NRL eta_perp coefficient", ME / (E * E * tau_e_coef), 1.03e-4, rtol=5e-3)
check("Bohm D_B = T/(16 B), NRL 6.25e6 cm^2/s at 1 eV, 1 G", 1.0 / (16 * 1e-4) * 1e4, 6.25e6, rtol=1e-12)
# Plasma parameter and coupling at the PlasmaPy point
N_D = 4 * pi / 3 * 1e19 * lam_d**3
print(f"   lambda_D = {lam_d:.5e} m, N_D = {N_D:.4e}, ln(12 pi n lambda_D^3) = {log(9 * N_D):.4f}")

# Epperlein-Haines alpha_par, beta_par as tabulated by Sadler, Walsh & Li 2021 Table I
ZS = [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 20, 30, 60]
ALPHA = [0.5061, 0.4295, 0.3950, 0.3750, 0.3618, 0.3524, 0.3454, 0.3399, 0.3319, 0.3263, 0.3221, 0.3144, 0.3081, 0.3015]
BETA = [0.7029, 0.9054, 1.0180, 1.0923, 1.1456, 1.1861, 1.2180, 1.2439, 1.2834, 1.3121, 1.3341, 1.3770, 1.4139, 1.4547]
lorentz = 3 * pi / 32
print("== Spitzer: alpha_par(Z) against the Lorentz limit 3 pi/32 and Spitzer-Harm gamma_E")
for z, sh in [(1, 0.582), (2, 0.683), (4, 0.785)]:  # Spitzer & Harm 1953 gamma_E (recalled)
    check(f"(3pi/32)/alpha_par(Z={z}) = gamma_E", lorentz / ALPHA[ZS.index(z)], sh, rtol=6e-3)
check("alpha_par(60) -> 3pi/32", ALPHA[-1], lorentz, rtol=0.03)
check("beta_par(60) -> 3/2", BETA[-1], 1.5, rtol=0.04)
check("1/alpha_par(1) against NRL's eta_perp/eta_par = 1.96", 1 / ALPHA[0], 1.96, rtol=0.01)


def gamma_wedge(chi, b):  # Sadler, Walsh & Li 2021 Eq. (14), Z = 1 row
    return (chi + BETA[0] * chi**2) / (b[0] + b[1] * chi + b[2] * chi**2 + chi**3)


print(f"   cross-gradient Nernst gamma_wedge(chi=1, Z=1) = {gamma_wedge(1, (0.5589, 1.0599, 2.1643)):.5f}")

# ---------------------------------------------------------------- 4. ICF: NIF N210808
print("== NIF N210808 (Abu-Shawareb et al. 2022 PRL 129, 075001, Table I)")
GBAR = 1e14  # Pa
for name, p, V, Y, El, Ehs_pub in [("N210207", 351, 3.3, 0.170, 1.93, 17.4), ("N210307", 353, 2.7, 0.145, 1.909, 14.6),
                                   ("N210808", 569, 6.4, 1.37, 1.917, 55.0)]:
    Ehs = 1.5 * p * GBAR * V * 1e5 * 1e-18 / 1e3  # kJ
    check(f"{name} E_hs = (3/2) p V", Ehs, Ehs_pub, rtol=0.03)
    ratio = Y * 1e3 / 5 / (2 * Ehs_pub)
    print(f"   {name}: E_alpha/(2 E_hs) = {ratio:.3f} (ignition > 1.4), G_target = {Y / El:.3f}")
check("N210808 G_target", 1.37 / 1.917, 0.72, atol=0.006)  # 0.715 from the rounded table
check("N210207 G_target", 0.170 / 1.93, 0.089, atol=0.0012)  # yields in Table I are rounded
m_bar = (M_D + M_T) / 2
check("c_DT = 3 k/m_bar (GJ/(g keV))", 3 * KEV / (m_bar * 1e3) / 1e9, 0.115, rtol=0.01)
check("Q_alpha coefficient 3.5 MeV/(4 m_bar^2)", 3.5e3 * KEV / (4 * (m_bar * 1e3) ** 2) / 1e9, 8.2e24, rtol=0.025)
check("Q_B coefficient 5.34e-37/m_bar^2", 5.34e-37 / m_bar**2 * 1e3 * 1e-12, 3.1e7, rtol=0.02)
# Hot-spot state from p, T_DD and V: rho R and Coulomb coupling
p, t_dd, V = 569 * GBAR, 8.94, 6.4e5 * 1e-18
n_i = p / (2 * t_dd * KEV)
rho = n_i * m_bar
R = (3 * V / (4 * pi)) ** (1 / 3)
check("N210808 rho R_hs from p/(2T_DD), V", rho * R / 10, 0.44, rtol=0.05)  # kg/m^2 -> g/cm^2
a_ws = (3 / (4 * pi * n_i)) ** (1 / 3)
gam = E * E / (4 * pi * EPS0 * a_ws * t_dd * KEV)
print(f"   hot spot: n_i = {n_i:.3e} m^-3, rho = {rho / 1e3:.1f} g/cc, R = {R * 1e6:.1f} um, Gamma = {gam:.2e}, "
      f"lambda_D = {debye(n_i, t_dd * 1e3) * 1e9:.2f} nm")


def bosch_hale_dt(t):  # Bosch & Hale 1992, as in plasmoid::fusion (m^3/s)
    bg, mr = 34.3827, 1_124_656.0
    c = [1.17302e-9, 1.51361e-2, 7.51886e-2, 4.60643e-3, 1.35000e-2, -1.06750e-4, 1.36600e-5]
    th = t / (1 - t * (c[1] + t * (c[3] + t * c[5])) / (1 + t * (c[2] + t * (c[4] + t * c[6]))))
    xi = (bg * bg / (4 * th)) ** (1 / 3)
    return 1e-6 * c[0] * th * sqrt(xi / (mr * t**3)) * exp(-3 * xi)


Ysim = (n_i / 2) ** 2 * bosch_hale_dt(t_dd) * V * 89e-12 * 17.59e3 * KEV
print(f"   0-D yield n_D n_T <sv>(T_DD) V tau_BW x 17.59 MeV = {Ysim / 1e6:.2f} MJ (measured 1.37 MJ)")
# ITFX (Patel 2020, as Eq. 5-6 of the PRL): Y_amp = exp(0.9 ITFX^0.47); chi_noalpha ~ ITFX_noalpha^0.34
for itfx in (0.3, 1.0, 2.0):
    noa = itfx * exp(-0.9 * itfx**0.47) + 0.007 * itfx**1.65
    print(f"   ITFX_alpha = {itfx}: ITFX_no-alpha = {noa:.4f}, Y_amp = {exp(0.9 * itfx**0.47):.3f}, "
          f"chi_no-alpha = {noa**0.34:.3f}")

print("== NIF ignition shots: target gain = yield / laser energy (LLNL)")
for date, el, y, g in [("2022-12-05", 2.05, 3.15, 1.5), ("2023-07-30", 2.05, 3.88, None), ("2023-10-08", 1.9, 2.4, None),
                       ("2023-10-30", 2.2, 3.4, None), ("2024-02-12", 2.2, 5.2, None), ("2024-11-18", 2.2, 4.1, None),
                       ("2025-02-23", 2.05, 5.0, 2.44), ("2025-04-07", 2.08, 8.6, 4.13), ("2025-10-01", 2.065, 3.5, 1.74)]:
    print(f"   {date}: G = {y / el:.3f}" + (f" (published {g})" if g else ""))
check("N221204 gain (PRL 132, 065102: 3.1 MJ / 2.05 MJ ~ 1.5)", 3.1 / 2.05, 1.5, atol=0.02)
check("2025-02-23 gain", 5.0 / 2.05, 2.44, atol=0.005)
check("2025-04-07 gain", 8.6 / 2.08, 4.13, atol=0.01)
print(f"   2025-10-01: 3.5/2.065 = {3.5 / 2.065:.3f} against the published 1.74: the page's yield or gain is rounded")
print(f"   2026-06-20: 7.9 MJ at gain ~3.8 implies {7.9 / 3.8:.2f} MJ of laser energy")

# ---------------------------------------------------------------- 5. Long-pulse energy turnover
print("== Energy turnover = heating power x duration")
for name, gj, s in [("W7-X Feb 2023", 1.3, 480), ("W7-X 2025", 1.8, 360), ("WEST Feb 2025 #61299", 2.6, 1337)]:
    print(f"   {name}: {gj} GJ over {s} s = {gj * 1e3 / s:.2f} MW mean")

# ---------------------------------------------------------------- 6. Ionization potential depression
print("== IPD: Stewart-Pyatt (high-density form) and Ecker-Kroll, Perez-Callejo et al. 2024 PRE 109, 045204")
K_E2 = E / (4 * pi * EPS0)  # e^2/(4 pi eps0) in eV m


def ipd_sp(z, n_e):
    return 1.5 * (z + 1) * K_E2 * (4 * pi * n_e / (3 * (z + 1))) ** (1 / 3)


def ipd_ek(z, n_e, n_i):
    return (z + 1) * K_E2 * (4 * pi * (n_e + n_i) / 3) ** (1 / 3)


n_mg = 1.738e3 / (24.305 * AMU)  # solid magnesium, ions per m^3
for zbar in (8.0, 9.0, 10.0):
    sp, ek = ipd_sp(9, zbar * n_mg), ipd_ek(9, zbar * n_mg, n_mg)
    print(f"   Li-like Mg (Z=9), Zbar={zbar}: SP {sp:.1f} eV, EK {ek:.1f} eV, EK/SP {ek / sp:.3f}")
check("EK/SP model ratio at Zbar = 9 against the fitted 1.30/0.88",
      ipd_ek(9, 9 * n_mg, n_mg) / ipd_sp(9, 9 * n_mg), 1.30 / 0.88, rtol=0.05)

# ---------------------------------------------------------------- 7. Thomas-Fermi mean ionization (More fit)
def zbar_more(Z, A, rho_gcc, t_ev):
    """More (1985) fit to Thomas-Fermi average-atom ionization; rho in g/cm^3, T in eV.
    Coefficients as in LANL Multi-BGK src/zBar.c (More 1985 not fetched)."""
    a1, a2, a3, a4 = 0.003323, 0.9718, 9.26148e-5, 3.10165
    b0, b1, b2, c1, c2 = -1.7630, 1.43175, 0.31546, -0.366667, 0.983333
    alpha, beta = 14.3139, 0.6624
    R = rho_gcc / (Z * A)
    T0 = t_ev / Z ** (4 / 3)
    Tf = T0 / (1 + T0)
    Af = a1 * T0**a2 + a3 * T0**a4
    B = -exp(b0 + b1 * Tf + b2 * Tf**7)
    C = c1 * Tf + c2
    Q1 = Af * R**B
    Q = (R**C + Q1**C) ** (1 / C)
    x = alpha * Q**beta
    return Z * x / (1 + x + sqrt(1 + 2 * x))


print("== More's Thomas-Fermi fit, solid aluminium (2.70 g/cc)")
for t in (0.0, 1.0, 10.0, 100.0, 1000.0, 1e4):
    print(f"   T = {t:>7} eV: Zbar = {zbar_more(13, 26.98, 2.70, t):.3f}")
assert zbar_more(13, 26.98, 2.70, 1e5) > 12.9

# ---------------------------------------------------------------- 8. Jeans twin of the Langmuir root
print("== Kinetic Jeans modes: 1 - (1/kappa^2)[1 + xi Z(xi)] = 0, kappa = k/k_J, xi = omega/(sqrt2 k sigma)")


def jeans_growth(kappa):
    """Unstable branch, omega = i gamma: kappa^2 = 1 - sqrt(pi) y erfcx(y), gamma = sqrt2 kappa y
    in units of omega_J = sqrt(4 pi G rho) = k_J sigma. Bisection on y > 0."""
    f = lambda y: 1 - sqrt(pi) * y * w(complex(0, y)).real - kappa**2
    lo, hi = 1e-12, 1e6
    for _ in range(200):
        mid = sqrt(lo * hi)
        lo, hi = (mid, hi) if f(mid) > 0 else (lo, mid)
    return sqrt(2) * kappa * lo


for kappa in (0.01, 0.5, 0.9, 0.99):
    print(f"   k/k_J = {kappa}: growth rate = {jeans_growth(kappa):.6f} omega_J")
check("cold limit gamma -> omega_J as k -> 0", jeans_growth(1e-3), 1.0, rtol=2e-3)
check("marginal at k = k_J (Jeans 1902)", jeans_growth(0.9999), 0.0, atol=2e-3)
# Gravitational plasma frequency identity used by both crates: lambda_D = v_te/omega_pe
n0, t0 = 1e19, 100.0
omega_pe = sqrt(n0 * E * E / (EPS0 * ME))
check("lambda_D = sqrt(T/m)/omega_pe", debye(n0, t0), sqrt(t0 * E / ME) / omega_pe, rtol=1e-14)
print("all checks passed")
