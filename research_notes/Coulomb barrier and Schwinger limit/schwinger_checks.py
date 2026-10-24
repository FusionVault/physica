"""Strong-field QED check numbers. Pure stdlib. Run: python3 -I sfqed_checks.py
Constants: CODATA 2022 (e, h, c exact in SI).
Every printed line is a "computed" value for the notes; published values are compared as ratios.
"""
from math import pi, sqrt, exp, cos, radians, asin, log

# --- constants (CODATA 2022) ---
c = 299_792_458.0
e = 1.602_176_634e-19
h = 6.626_070_15e-34
hbar = h / (2 * pi)
m_e = 9.109_383_7139e-31
alpha = 7.297_352_5643e-3
eps0 = 8.854_187_8188e-12
mu0 = 1 / (eps0 * c**2)
mc2_eV = m_e * c**2 / e            # 510998.95 eV
lamC = hbar / (m_e * c)            # reduced Compton wavelength

def line(k, v, pub=None, unit=""):
    s = f"{k:<58s} {v:.6g} {unit}"
    if pub is not None:
        s += f"   published {pub:.4g}  ratio {v/pub:.4f}"
    print(s)

print("== 1. Critical fields ==")
E_S = m_e**2 * c**3 / (e * hbar)
B_S = m_e**2 * c**2 / (e * hbar)
line("E_S = m^2 c^3/(e hbar)  [V/m]", E_S, 1.3e18)
line("E_S [V/cm]", E_S / 100, 1.3e16)
line("B_S = m^2 c^2/(e hbar)  [T]", B_S, 4.4e9)
line("B_S [G]", B_S * 1e4, 4.4e13)
I_half = 0.5 * eps0 * c * E_S**2 / 1e4      # W/cm^2, I = eps0 c E^2 / 2
I_full = eps0 * c * E_S**2 / 1e4            # W/cm^2, I = eps0 c E^2 (= Gaussian c E^2/4pi)
line("I_S = eps0 c E_S^2 (Di Piazza I=F^2/4pi) [W/cm^2]", I_full, 4.6e29)
line("I_S = eps0 c E_S^2 / 2 (cycle-avg linear) [W/cm^2]", I_half, 2.3e29)
line("reduced Compton wavelength [m]", lamC)
line("work over lamC at E_S / mc^2 (should be 1)", e * E_S * lamC / (m_e * c**2))

print("\n== 2. Schwinger rate, constant E (spinor QED) ==")
pref = c / (4 * pi**3 * lamC**4)  # = (eE_S)^2/(4 pi^3 hbar^2 c)  [m^-3 s^-1]
line("prefactor (eE_S)^2/(4 pi^3 hbar^2 c) [m^-3 s^-1]", pref)
line("same, via SI e,hbar,c", (e * E_S)**2 / (4 * pi**3 * hbar**2 * c))
line("same, as alpha eps0 E_S^2/(pi^2 hbar)", alpha * eps0 * E_S**2 / (pi**2 * hbar))
line("Melissinos (8') alpha E_c^2/(pi hbar) [cm^-3 s^-1] (HL field)", alpha * eps0 * E_S**2 / (pi * hbar) * 1e-6, 3.45e50)

def w_sum(x, nmax=2000):     # vacuum-decay rate (2 Im L / hbar), Schwinger 1951
    return pref * x * x * sum(exp(-n * pi / x) / n**2 for n in range(1, nmax))
def w_pair(x):               # mean pair-number rate (Nikishov 1970; Cohen & McGady 2008): n=1 term
    return pref * x * x * exp(-pi / x)
def w_scalar(x, nmax=2000):  # scalar QED: (eE)^2/(8 pi^3) sum (-1)^(n+1)/n^2 exp(-n pi/x)
    return pref / 2 * x * x * sum((-1)**(n + 1) * exp(-n * pi / x) / n**2 for n in range(1, nmax))

print(f"{'E/E_S':>6} {'exp(-pi/x)':>11} {'Gamma_pair m^-3s^-1':>20} {'w_sum m^-3s^-1':>15} {'w_sum/G_pair':>12} {'G_pair cm^-3s^-1':>17} {'pairs in (1um)^3 x 1um/c':>24}")
for x in (1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01):
    gp, ws = w_pair(x), w_sum(x)
    n_lam = gp * (1e-6)**4 / c
    print(f"{x:6.2f} {exp(-pi/x):11.4e} {gp:20.4e} {ws:15.4e} {ws/gp:12.6f} {gp*1e-6:17.4e} {n_lam:24.4e}")
line("sum at x=1: zeta-like ratio w_sum/Gamma_pair", w_sum(1) / w_pair(1))
line("scalar/spinor (sum) at x=0.1", w_scalar(0.1) / w_sum(0.1))
# field for one pair in (1 um)^3 during 1 um/c (illustrative focal volume) - bisection
lo, hi = 0.005, 0.5
for _ in range(200):
    mid = (lo + hi) / 2
    if w_pair(mid) * (1e-6)**4 / c > 1: hi = mid
    else: lo = mid
line("E/E_S for 1 pair in (1um)^3 over 1um/c", mid)
line("  corresponding I = eps0 c E^2/2 [W/cm^2]", I_half * mid**2)

print("\n== 3. Laser parameters ==")
def a0_lin(I_Wcm2, lam_um):   # peak a0, linear pol, I = eps0 c E0^2 / 2
    E0 = sqrt(2 * I_Wcm2 * 1e4 / (eps0 * c))
    return e * E0 * lam_um * 1e-6 / (2 * pi * m_e * c**2)
def a0_dp(I_Wcm2, lam_um):    # Di Piazza convention I = E0^2/4pi (Gaussian) = eps0 c E0^2 (SI)
    return a0_lin(I_Wcm2, lam_um) / sqrt(2)
line("a0 coefficient, lin, per sqrt(1e18 W/cm2) per um", a0_lin(1e18, 1.0), 0.855)
line("a0 coefficient, Di Piazza conv., per sqrt(1e20) per um", a0_dp(1e20, 1.0), 6.0)
hw800 = h * c / 800e-9 / e
line("photon energy at 800 nm [eV]", hw800, 1.55)
line("Di Piazza xi0 = 7.5 sqrt(I/1e20)/w[eV] at 1.55 eV", a0_dp(1e20, 0.8), 7.5 / 1.55)

def chi_e(E_GeV, a0, hw_eV, theta_deg=0.0):  # theta = angle from head-on; plane-wave peak
    gamma = E_GeV * 1e9 / mc2_eV
    beta = sqrt(1 - 1 / gamma**2)
    return (1 + beta * cos(radians(theta_deg))) * gamma * a0 * hw_eV / mc2_eV
line("Di Piazza chi0 coefficient (1 GeV, 1e20 W/cm2, 800nm, head-on)", chi_e(1, a0_dp(1e20, 0.8), hw800), 5.9e-2)

print("-- CoReLS 2021 (lambda assumed 0.8 um)")
line("a0 at 1.1e23 W/cm2, lin convention", a0_lin(1.1e23, 0.8))
line("a0 at 1.1e23 W/cm2, Di Piazza convention", a0_dp(1.1e23, 0.8))
line("E0 peak at 1.1e23 (lin) / E_S", sqrt(2 * 1.1e27 / (eps0 * c)) / E_S)
line("I_S(eps0 c E^2/2) / 1.1e23", I_half / 1.1e23)

print("-- LUXE (16.5 GeV, 800 nm, 17.2 deg crossing)")
for I, xi_pub, chi_pub in ((1.33e20, 7.9, 1.5), (1.2e21, 23.6, 4.45)):
    a = a0_lin(I, 0.8)
    line(f"  xi at {I:.3g} W/cm2", a, xi_pub)
    line(f"  chi_e(16.5 GeV, xi_pub={xi_pub}, 17.2 deg)", chi_e(16.5, xi_pub, hw800, 17.2), chi_pub)

print("-- E144 (46.6 GeV, 527 nm, 17 deg)")
hw527 = h * c / 527e-9 / e
gam144 = 46.6e9 / mc2_eV
eta144 = (1 + cos(radians(17))) * gam144 * hw527 / mc2_eV
line("  photon energy 527 nm [eV]", hw527)
line("  eta = k.p/m^2 (ESPPU 'eta')", eta144, 0.83)
line("  chi = xi*eta with xi=0.36 (ESPPU point)", 0.36 * eta144, 0.30)
line("  LUXE-quoted E144 chi (xi=0.4)", 0.4 * eta144, 0.25)
# threshold: gamma(29 GeV) + n laser photons -> e+e-, s = 2 w1 n w2 (1+cos th)
nthr = 2 * mc2_eV**2 / (29e9 * hw527 * (1 + cos(radians(17))))
line("  min laser photons for 29 GeV photon (17 deg)", nthr)
line("  E144 positron signal significance 106/14 [sigma]", 106 / 14)

print("-- E-320 (10 GeV, 800 nm, ~30 deg)")
line("  chi for a0=sqrt(10)", chi_e(10, sqrt(10), hw800, 30), 0.5)
line("  chi for a0=sqrt(20) (if a0^2 is cycle-avg)", chi_e(10, sqrt(20), hw800, 30), 0.5)
line("  ESPPU E320 point chi = 2*0.15", 2 * 0.15)
line("  ESPPU E320II point chi = 16*0.15", 16 * 0.15)
line("  ESPPU NA63 point chi = 20*0.1", 20 * 0.1)

print("-- ESPPU: 5 GeV head-on with 1e23 W/cm2, 800 nm -> chi ~ 10")
line("  lin a0", chi_e(5, a0_lin(1e23, 0.8), hw800), 10)
line("  Di Piazza a0", chi_e(5, a0_dp(1e23, 0.8), hw800), 10)

print("-- ESPPU Table 1 xi reproduction: FWHM 2 um spot holding 60% of energy, Gaussian in time")
def I_peak(E_J, tau_fs, fwhm_um=2.0, frac=0.6):
    P = frac * E_J * 2 * sqrt(log(2) / pi) / (tau_fs * 1e-15)
    w0 = fwhm_um * 1e-4 / sqrt(2 * log(2))   # cm, 1/e^2 intensity radius
    return 2 * P / (pi * w0**2)
for name, EJ, tau, xi_pub, lam in (("FACET-II", 0.3, 30, 6, 0.8), ("GEMINI", 15, 45, 32, 0.8),
                                    ("ZEUS", 75, 25, 100, 0.8), ("APOLLON", 70, 22, 100, 0.8),
                                    ("CoReLS", 80, 20, 110, 0.8), ("ELI-NP", 200, 20, 180, 0.8),
                                    ("ELI-BL (1.057um?)", 1500, 150, 230, 1.057), ("NSF OPAL (per beam)", 500, 20, 280, 0.8),
                                    ("VULCAN2020 (1.053um?)", 400, 20, 330, 1.053), ("SEL as tabulated 25 kJ", 25000, 25, 1800, 0.8),
                                    ("SEL if 100 PW = 2.5 kJ/25 fs", 2500, 25, 1800, 0.8)):
    I = I_peak(EJ, tau)
    print(f"  {name:<30s} I={I:9.3e} W/cm2  xi_DP={a0_dp(I, lam):7.1f}  xi_lin={a0_lin(I, lam):7.1f}  pub {xi_pub}  ratio_DP {a0_dp(I, lam)/xi_pub:.3f}")
print("  chi per unit xi for 10 GeV head-on at 800 nm:", f"{chi_e(10, 1, hw800):.4f}")

print("-- Ritus-Narozhny")
chi_RN = alpha**-1.5
line("  chi where alpha chi^(2/3) = 1", chi_RN, 1600)
gam_RN = chi_RN * mc2_eV / (2 * 100 * hw800)
line("  electron energy for chi=1600 head-on, a0=100, 800nm [TeV]", gam_RN * mc2_eV / 1e12)
line("  alpha chi^(2/3) at chi=1", alpha)
line("  alpha chi^(2/3) at chi=4.45 (LUXE ph-1)", alpha * 4.45**(2 / 3))

print("-- Keldysh / dynamically assisted (Schutzhold, Gies, Dunne 2008)")
line("  Keldysh gamma = (hw/mc^2)(E_S/E), 800nm, E=0.1E_S", hw800 / mc2_eV / 0.1)
def A_inst(g):   # in units of m^2/(qE); valid g >= pi/2, else pi
    if g <= pi / 2: return pi
    return 2 * asin(pi / (2 * g)) + pi / (2 * g * g) * sqrt(4 * g * g - pi * pi)
line("  A_inst(pi/2)/pi", A_inst(pi / 2) / pi, 1.0)
line("  1 - A_inst(pi)/pi (published 'about 40%')", 1 - A_inst(pi) / pi, 0.40)
line("  A_inst(100)*100/(2pi) -> 1 (1/gamma tail)", A_inst(100) * 100 / (2 * pi), 1.0)
line("  hbar*omega of assisting pulse for gamma=pi at E=0.1E_S [keV]", pi * 0.1 * mc2_eV / 1e3)
line("  exp(-pi*10) (Schwinger exponent at 0.1 E_S)", exp(-10 * pi))
line("  exp(-A_inst(pi)*10) (assisted, gamma=pi)", exp(-10 * A_inst(pi)))

print("\n== 4. Vacuum birefringence (Euler-Heisenberg, weak field) ==")
A_e = 2 * alpha**2 * hbar**3 / (45 * mu0 * m_e**4 * c**5)
line("A_e = 2 alpha^2 hbar^3 / (45 mu0 m^4 c^5) [T^-2]", A_e, 1.32e-24)
line("Delta n_u = 3 A_e [T^-2]", 3 * A_e, 4e-24)
line("Delta n at 1 T", 3 * A_e * 1.0)
line("Delta n at 2.5 T", 3 * A_e * 2.5**2, 2.5e-23)
line("n_par - 1 = 3.5 A_e B^2 at 2.5 T", 3.5 * A_e * 6.25)
line("n_perp - 1 = 2 A_e B^2 at 2.5 T", 2 * A_e * 6.25)
line("PVLAS: sigma(Delta n)/Delta n_QED at 2.5 T", 17e-23 / (3 * A_e * 6.25), 7)
line("PVLAS: (meas - QED)/sigma", (12e-23 - 3 * A_e * 6.25) / 17e-23)
line("PVLAS per T^2: 12e-23/6.25", 12e-23 / 6.25)
line("PVLAS per T^2: 17e-23/6.25", 17e-23 / 6.25)
line("EH Lagrangian coeff 2 alpha^2 eps0^2 hbar^3/(45 m^4 c^5) [m^3/J]", 2 * alpha**2 * eps0**2 * hbar**3 / (45 * m_e**4 * c**5))
line("B/B_S at 1e14 G (magnetar scale)", 1e10 / B_S)
line("STAR cos4phi significance 16.8/2.5 [sigma]", 16.8 / 2.5)

print("\n== 5. Supercritical Coulomb fields ==")
line("point-nucleus Dirac 1s collapse Z = 1/alpha", 1 / alpha)
for Z in (82, 92):
    line(f"  point-nucleus 1s binding Z={Z} [keV]", mc2_eV * (1 - sqrt(1 - (Z * alpha)**2)) / 1e3)
line("U+U total Z", 184)
v = 2 * 35e-15 / 0.25e-20
line("  v implied by tau=2R/v, R=35 fm, tau=0.25e-20 s [units of c]", v / c, 0.1)
ke = e**2 / (4 * pi * eps0) / e / 1e6 * 1e15   # MeV fm
line("  e^2/(4 pi eps0) [MeV fm]", ke, 1.44)
Ecm = 92 * 92 * ke / 16.5
line("  U+U Coulomb energy at R=16.5 fm [MeV]", Ecm)
line("  lab energy per nucleon for head-on U+U to 16.5 fm [MeV/u]", 2 * Ecm / 238)

print("\n== 6. Graphene Schwinger analogue ==")
def gamma_2d(E, vF, species=4):  # per area per time; Allor-Cohen-McGady: 4 species -> (qE)^1.5 / (pi^2 hbar^1.5 v^0.5)
    return species / 4 * (e * E)**1.5 / (pi**2 * hbar**1.5 * sqrt(vF))
for E in (1e4, 1e5, 1e6):
    for vF in (1e6, 3e4):
        print(f"  E={E:.0e} V/m vF={vF:.0e} m/s: rate {gamma_2d(E, vF):.4e} m^-2 s^-1 = {gamma_2d(E, vF)*1e-4:.4e} cm^-2 s^-1")
line("  n = j_m/(e vF), j_m = 0.05 mA/um, vF=1e6 [cm^-2]", 0.05e-3 / 1e-6 / (e * 1e6) * 1e-4, 3e10)
line("  j_c = n e vF, n=0.4e16 m^-2 [mA/um]", 0.4e16 * e * 1e6 * 1e3 * 1e-6)

print("\n== 7. NA63 effective crystal fields from chi = gamma E/E_S ==")
for name, EGeV, chi in (("Si <111>, 178.2 GeV e+", 178.2, 1.4), ("Ge <110>, 200 GeV e-", 200, 2.4)):
    g = EGeV * 1e9 / mc2_eV
    line(f"  {name}: E_lab = chi E_S / gamma [V/m]", chi * E_S / g)

print("\n== 8. Radiation reaction helpers (ESPPU 2025 eqs. 4-5) ==")
r_e = e**2 / (4 * pi * eps0 * m_e * c**2)
line("tau0 = 2 r_e/(3c) [s]", 2 * r_e / (3 * c), 6.4e-24)
g1 = lambda x: (1 + 12 * x + 31 * x**2 + 3.7 * x**3) ** (-4 / 9)
g2 = lambda x: (1 + 4.8 * (1 + x) * log(1 + 1.7 * x) + 2.44 * x**2) ** (-2 / 3)
for x in (0.1, 0.3, 1.0, 3.0):
    print(f"  Gaunt g(chi={x}): approx1 {g1(x):.4f}  approx2 {g2(x):.4f}")
for x in (0.1, 0.3, 1.0):
    print(f"  nonlinear BW quasi-static suppression exp(-8/(3 chi)) at chi={x}: {exp(-8/(3*x)):.4e}")
