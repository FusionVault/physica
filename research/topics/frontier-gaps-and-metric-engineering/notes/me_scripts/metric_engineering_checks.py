"""Recompute the published numbers in the metric-engineering literature (1994-2026).

Run: python3 -I metric_engineering_checks.py   (numpy + scipy only; reads no downloaded files)
Every line prints: what, computed value, published value (with source), and a verdict.
Geometric units (G = c = 1) inside the formulas; SI conversions are explicit.
"""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar, brentq

c, G, hbar, eps0, k_B = 299792458.0, 6.67430e-11, 1.054571817e-34, 8.8541878128e-12, 1.380649e-23
lP = np.sqrt(hbar * G / c**3)            # 1.616e-35 m
mP = np.sqrt(hbar * c / G)               # 2.176e-8 kg
Msun, MJ, ly, g0 = 1.98847e30, 1.89813e27, 9.4607304725808e15, 9.80665
C2G = c**2 / G                           # kg per metre of geometric length (1.3466e27)
C4G = c**4 / G                           # J per metre (1.2103e44)


def show(name, got, want=None, src="", tol=None):
    if want is None:
        print(f"{name:58s} {got: .4g}")
        return
    ratio = got / want
    ok = "" if tol is None else ("  OK" if abs(ratio - 1) <= tol else "  MISMATCH")
    print(f"{name:58s} {got: .4g}  vs published {want: .4g} [{src}]  ratio {ratio:.3g}{ok}")


def f_alc(r, R, s):
    return (np.tanh(s * (r + R)) - np.tanh(s * (r - R))) / (2 * np.tanh(s * R))


def df_alc(r, R, s):
    return s * (1 / np.cosh(s * (r + R))**2 - 1 / np.cosh(s * (r - R))**2) / (2 * np.tanh(s * R))


def I_alc(R, s):
    """int_0^inf r^2 f'^2 dr for Alcubierre's tanh wall (geometric length^1)."""
    pts = [max(R - 40 / s, 0), R, R + 40 / s]
    return sum(quad(lambda r: r * r * df_alc(r, R, s)**2, a, b, limit=400, epsabs=0, epsrel=1e-11)[0]
               for a, b in zip([0] + pts[:-1], pts)) + quad(lambda r: r * r * df_alc(r, R, s)**2, pts[-1], np.inf)[0]


print("=== 1. Alcubierre 1994 / Pfenning & Ford 1997 / Lobo & Visser 2004 ===")
v = 1.0
E_pf = lambda R, D, v: -(v**2 / 12) * (R**2 / D + D / 12)          # PF Eq. (28), geometric
show("PF Eq.28, R=100 m, Delta=100 v lP (kg, per v)", -E_pf(100, 100 * lP, 1) * C2G, 6.2e62, "PF text '6.2e65 g'", 0.15)
show("  same, in Planck lengths (per v)", -E_pf(100, 100 * lP, 1) / lP, 6.2e70, "PF text '6.2e70 L_Planck'", 0.15)
show("  same / Milky Way 2e45 g = 2e42 kg", -E_pf(100, 100 * lP, 1) * C2G / 2e42, 3e20, "PF text '3e20 M_galaxy'", 0.2)
show("PF Eq.28, R=100 m, Delta=1 m (solar masses)", -E_pf(100, 1, 1) * C2G / Msun, 0.25, "PF text 'a quarter'", 0.2)
# exact tanh wall versus PF's piecewise-linear wall, using PF's own matching Delta(sigma)
R = 100.0
for s in (2.0, 20.0, 200.0):
    D_pf = (1 + np.tanh(s * R)**2)**2 / (2 * s * np.tanh(s * R))
    exact = (v**2 / 12) * I_alc(R, s)
    show(f"  tanh wall sigma={s:g}: exact/PF Eq.28 (same Delta)", exact / -E_pf(R, D_pf, 1))
    show(f"  tanh wall sigma={s:g}: exact/(v^2 R^2 sigma/36)", exact / (R * R * s / 36))
    show(f"  tanh wall sigma={s:g}: Lobo-Visser 'v^2 R^2 sigma' / exact", (R * R * s) / exact)
# minimum over wall thickness for fixed R: E_min = k (v^2/12) R c^2/G
res = minimize_scalar(lambda ls: I_alc(1.0, np.exp(ls)), bounds=(np.log(0.05), np.log(50)), method="bounded")
k_min, s_min = res.fun, np.exp(res.x)
show("min over sigma of int r^2 f'^2 dr / R  (tanh, R=1)", k_min)
show("  at sigma*R =", s_min)
show("Bobrick-Martire optimum f=min(r0/r,1): int r^2 f'^2 / r0", 1.0)
show("  tanh minimum / BM optimum (BM claim: 'about a factor of 3')", k_min / 1.0)
show("  tanh at sigma*R=8 (Alcubierre's plotted example) / BM optimum", I_alc(1.0, 8.0) / 1.0)

print("\n=== 2. White 2011 (WFM101) / 2013 (WFM102) / White-Juday interferometer ===")
show("Alpha Centauri 4.37 ly at 10c (years)", 4.37 / 10, 0.43, "WFM101 '0.43 years'", 0.03)
# boost gamma = cosh(0.5 ln|1 - v^2 f^2|) = 1.0000001 inside the bubble (f = 1)
vs = brentq(lambda x: np.cosh(0.5 * np.log(abs(1 - x * x))) - (1 + 1e-7), 1e-4, 0.5)
show("WFM101 boost 1.0000001 needs v_s (units of c)", vs)
Ewfi = (vs**2 / 12) * I_alc(0.005, 1 / 0.005 * 1.0) * C2G     # 1 cm sphere, wall ~ radius
show("  Alcubierre negative mass for that 1 cm bubble (kg)", -Ewfi)
show("  as energy (J)", -Ewfi * c**2)
phase = 2 * np.pi * 1e-7 * 0.01 / 633e-9
show("  phase shift of 1e-7 over 1 cm at 633 nm (rad, single pass)", phase)
# WFM102 chart: 10 m diameter, v_apparent = 10c. Alcubierre energy floor for R = 5 m, v = 10:
Rw = 5.0
show("WFM102 R=5 m v=10: tanh minimum over all walls (kg)", (100 / 12) * k_min * Rw * C2G)
show("  chart dU/dt=0 curve at shell fraction 1 (read off slide 16)", 1e20)
show("  PF Eq.28 at Delta = 1e-5 R (kg)", -E_pf(Rw, 1e-5 * Rw, 10) * C2G)
show("  chart dU/dt=0 curve at fraction 1e-5 (read off slide 16)", 1.5e27)
# Lee & Cleaver 2014 arithmetic
uE = eps0 * (1e5)**2 / (2 * 0.1**2)
show("Lee-Cleaver u_E for 100 kV over 0.1 m (J/m^3)", uE, 4.4, "arXiv:1407.7772 Eq.1", 0.02)
show("  8 pi G u/c^4 (m^-2)", 8 * np.pi * G * uE / c**4, 9.2e-43, "Eq.2", 0.02)
M_sph = uE * 4 / 3 * np.pi * 0.1**3 / c**2
show("  R_S of that field energy in a 0.1 m-radius sphere (m)", 2 * G * M_sph / c**2, 1.92e-23, "text", None)
M_cyl = uE * np.pi * 0.1**2 * 10 / c**2
show("  R_S if the field filled a 0.1 m x 10 m cylinder (m)", 2 * G * M_cyl / c**2, 1.92e-23, "text", None)
show("  dimensionless metric perturbation 2GM/(c^2 r), sphere", 2 * G * M_sph / c**2 / 0.1)

print("\n=== 3. White et al. 2021 EPJC 81:677 (Casimir 'warp bubble') ===")
d = 4e-6
Ef = np.sqrt(hbar * c * np.pi**2 / (1080 * eps0 * d**4))
show("E = sqrt(hbar c pi^2/(1080 eps0 d^4)) at d = 4 um (V/m)", Ef)
show("  times d/2 (V)  [paper '~0.7 mV']", Ef * d / 2, 0.7e-3, "EPJC Sec.1 via Springer HTML", 0.05)
show("  times d   (V)", Ef * d)
u4 = np.pi**2 * hbar * c / (720 * d**4)
show("Casimir |u| between ideal plates 4 um apart (J/m^3)", u4)
show("  eps0 E^2/2 / |u|  (shows the 1080 = 3 x 360 split)", eps0 * Ef**2 / 2 / u4)
# Alcubierre requirement for a 1 um bubble at v = 1e-? : what speed could |u4| buy?
Rb, sb = 0.5e-6, 1 / 0.25e-6
peak = C4G / (32 * np.pi) * max(df_alc(r, Rb, sb)**2 for r in np.linspace(0, 2 * Rb, 20001))
show("Alcubierre peak |rho| for R=0.5 um, wall 0.5 um, v=c (J/m^3)", peak)
show("  speed whose peak |rho| equals the 4 um Casimir |u| (units of c)", np.sqrt(u4 / peak))
show("  that speed in m/s", np.sqrt(u4 / peak) * c)

print("\n=== 4. Van Den Broeck 1999 ===")
alpha, Dt, Rt, Rv, n = 1e17, 1e-15, 1e-15, 3e-15, 80
show("E_IV = (1/12)((R+D/2)^2/D + D/12), D=100 lP (kg)", ((Rv + 50 * lP)**2 / (100 * lP) + 100 * lP / 12) / 12 * C2G,
     6.3e29, "gr-qc/9905084 Sec.3", 0.05)
show("inner pocket diameter 2 alpha Rt (m)", 2 * alpha * Rt, 200, "'200 m'", 1e-9)


def vdb(w):
    """B and T^00 (units 1/Dt^2) at w = (Rt+Dt-r)/Dt, from VdB Eq. for T^00."""
    g = n * w**(n - 1) - (n - 1) * w**n
    g1 = n * (n - 1) * (w**(n - 2) - w**(n - 1))
    g2 = n * (n - 1) * ((n - 2) * w**(n - 3) - (n - 1) * w**(n - 2))
    B = alpha * g + 1
    Br, Brr, r = -alpha * g1, alpha * g2, 2 - w        # d/dr = -(1/Dt) d/dw; r in units of Dt (Rt = Dt)
    T = (Br**2 / B**4 - 2 * Brr / B**3 - 4 * Br / (B**3 * r)) / (8 * np.pi)
    R1212 = Br**2 / B**4 - Brr / B**3 - Br / (B**3 * r)
    return B, T, R1212


ws = np.linspace(1e-6, 1 - 1e-9, 2_000_001)
B, T, Rr = vdb(ws)
neg = T < 0
w0 = brentq(lambda w: vdb(w)[1], 0.9, 0.9999)
i = np.argmin(T)
show("T^00 changes sign at w =", w0, 0.981, "VdB", 0.005)
show("T^00 most negative at w =", ws[i], 0.349, "VdB", 0.02)
show("  its value x Dt^2", T[i], -4.9e2, "VdB", 0.05)
show("  B there", B[i])
dens = 4 * np.pi * (2 - ws)**2 * B**3 * T                   # E = Dt * int dens dw (geometric)
show("E_II,- (kg)", Dt * np.trapezoid(np.where(neg, dens, 0), ws) * C2G, -1.4e30, "VdB", 0.05)
show("E_II,+ (kg)", Dt * np.trapezoid(np.where(neg, 0, dens), ws) * C2G, 4.9e30, "VdB", 0.05)
j = np.argmax(np.abs(Rr))
show("max |R_1212| at w =", ws[j], 0.348, "VdB", 0.02)
show("r_c,min / Dt = 1/sqrt|R_1212|", 1 / np.sqrt(abs(Rr[j])), 1 / 72.5, "VdB 'Dt/72.5'", 0.05)
show("r_c,min (m) with Dt = 1e-15 m", Dt / np.sqrt(abs(Rr[j])), 1.4e-34, "VdB", None)
show("  Dt implied by the paper's r_c,min = 1.4e-34 m", 1.4e-34 * 72.5)
tau0 = 0.1 * 1.4e-34
show("QI RHS 3 hbar/(32 pi^2 c tau0^4), tau0 = 0.1 r_c (kg/m^3)", 3 * hbar / (32 * np.pi**2 * c * tau0**4), 9.2e94, "VdB", 0.1)
show("QI LHS 4.9e2/Dt^2 x c^2/G with Dt = 1e-15 m (kg/m^3)", 4.9e2 / Dt**2 * C2G, 6.6e93, "VdB", None)
show("QI LHS with Dt = 1e-32 m (kg/m^3)", 4.9e2 / 1e-32**2 * C2G, 6.6e93, "VdB", 0.05)

print("\n=== 5. Natario 2002 zero-expansion drive (f: 0 inside -> 1/2 outside) ===")


def nat_energy(R, s, v=1.0):
    """E = -(v^2/2) int f'^2 r^2 - (v^2/3) int (f' + r f''/2)^2 r^2, with f = (1 - f_alc)/2 (geometric)."""
    h = 1e-4 / s
    fp = lambda r: -0.5 * df_alc(r, R, s)
    fpp = lambda r: -0.5 * (df_alc(r + h, R, s) - df_alc(r - h, R, s)) / (2 * h)
    lo, hi = max(R - 40 / s, 0), R + 40 / s
    a = quad(lambda r: fp(r)**2 * r * r, lo, hi, limit=400)[0]
    b = quad(lambda r: (fp(r) + r * fpp(r) / 2)**2 * r * r, lo, hi, limit=400)[0]
    vort = quad(lambda r: (4 * fp(r) + r * fpp(r))**2 * r * r, lo, hi, limit=400)[0]   # SSV: -(1/32pi) int w^2
    return -(v**2 / 2) * a - (v**2 / 3) * b, -(v**2 / 12) * vort


for R_, s_ in ((100, 2.0), (100, 20.0), (10, 2.0)):
    En, Ev = nat_energy(R_, s_)
    show(f"Natario R={R_} sigma={s_}: E (M_sun, v=1)", En * C2G / Msun)
    show("  SSV vorticity identity -(1/32pi) int w^2 / direct", Ev / En)
    show("  Natario / Alcubierre (same R, sigma)", En / (-(1 / 12) * I_alc(R_, s_)))
    show("  E / (v^2 R^4 sigma^3)", En / (R_**4 * s_**3))

print("\n=== 6. Krasnikov tube (Everett & Roman 1997) ===")


def T_tt(n, delta, eps=1.0):
    k = lambda rho: 1 - 0.5 * (1 - delta / 2) * (np.tanh(2 * (2 * (n * eps - rho) / eps - 1)) + 1)
    rho, h = (n - 1) * eps, 1e-4 * eps
    k0, k1, k2 = k(rho), (k(rho + h) - k(rho - h)) / (2 * h), (k(rho + h) - 2 * k(rho) + k(rho - h)) / h**2
    return (-4 * (1 + k0) / rho * k1 + 3 * k1**2 - 4 * (1 + k0) * k2) / (32 * np.pi * (1 + k0)**2)


fit = lambda n, d: (0.271 * (3.504 - 4.034 * n - 0.008 * d + 0.017 * n * d - 0.872 * d * d + n * d * d)) / (8 * np.pi * (n - 1) * (1.018 + 0.491 * d)**2)
for n_, d_ in ((10, 0.01), (100, 0.01), (10, 1.0), (3, 0.5)):
    show(f"T_tt end cap, n={n_}, delta={d_} (x eps^2): Eq.T vs ER fit", T_tt(n_, d_), fit(n_, d_), "ER Eq.Tttji", 0.02)
show("T_tt x 8 pi eps^2 at n=1000, delta->0", T_tt(1000, 1e-6) * 8 * np.pi, -1.0, "ER '~ -1/(8 pi eps^2)'", 0.1)
show("QI: eps <~ lP/sigma^2 at sigma=0.01 (m)", lP / 0.01**2, 1e-31, "ER '1e4 lP ~ 1e-31 m'", 0.7)
E_kt = lambda rho, D, eps, a: a * rho * D / (4 * eps)               # T*V = (1/8pi eps^2)(2 pi a eps rho D)
show("tube energy alpha rho D/eps, 1 m x 1 m, eps=100 lP (Planck masses / alpha)", 1 / (100 * lP) * 1 / lP, 1e68, "ER", None)
show("  with the 1/4 the formula T V actually gives (/alpha)", E_kt(1 / lP, 1 / lP, 100, 1), 1e68, "ER", None)
show("  ER's rounded version in grams (/alpha)", 1 / (100 * lP) / lP * mP * 1e3, 1e63, "ER '1e63 g'", None)
show("  galaxies for alpha=0.01, D = 4e16 m (M_gal = 1e12 Msun)",
     0.01 * (1 / lP) * (4e16 / lP) / 100 * mP / (1e12 * Msun), 1e32, "ER '1e32 M_galaxy'", None)
show("  eta for QI-safe tube, tau0 = eps = 1 cm: lP^2 eps^2/tau0^4", lP**2 * 0.01**2 / 0.01**4, 1e-66, "ER", None)

print("\n=== 7. 'Positive-energy' solitons ===")
show("Lentz R^2/w c^2/G for R=100, w=1 (M_sun)", 100**2 / 1 * C2G / Msun)
show("  C implied by '(few)x1e-1 M_sun' (taking 0.3)", 0.3 / (100**2 * C2G / Msun))
show("Fell-Heisenberg E_sun = M_sun c^2 (J)", Msun * c**2, 1.78e47, "FH Sec.Numerics", 0.01)
show("  E_total/E_sun (paper: 'four orders of magnitude smaller')", 9.25e43 / (Msun * c**2))
show("  E_total as mass (kg)", 9.25e43 / c**2)
show("  rho_max 3.2e26 kg/m^3 / nuclear density 2.3e17", 3.2e26 / 2.3e17)

print("\n=== 8. Wormholes: Morris-Thorne, Visser, Ford-Roman 1996 ===")
show("MT throat tension c^4/(8 pi G b0^2), b0 = 10 m (Pa)", C4G / (8 * np.pi * 100), 4.82e40, "library test", 5e-3)
show("Visser cube edge tension c^4/(8G) (J/m)", C4G / 8, 1.513e43, "library test", 1e-3)
show("  12 edges of a 1 m cube as mass (kg)", 12 * C4G / 8 / c**2)
sig = lambda a, M=0.0: -C4G / (2 * np.pi * a) * np.sqrt(1 - 2 * G * M / (a * c * c))        # J/m^2
show("thin-shell sigma, a = 1 m, M = 0 (J/m^2)", sig(1.0))
show("  total shell mass 4 pi a^2 sigma/c^2 (kg)", 4 * np.pi * sig(1.0) / c**2)
show("  same for a = 1 m around M = 0.1 a c^2/G (kg)", 4 * np.pi * sig(1.0, 0.1 * C2G) / c**2)
f_ = 0.01
show("FR b=r0^2/r: r0 <= sqrt(3/(4pi)) lP/f^2 at f=0.01 (m)", np.sqrt(3 / (4 * np.pi)) * lP / f_**2, 1e-31, "FR '1e4 lP = 1e-31 m'", 0.3)
show("  FR's approximation lP/(2 f^2) (m)", lP / (2 * f_**2))
a0 = lambda r0: (3 * r0 / (8 * np.pi * f_**4) * lP**2)**(1 / 3)
show("FR absurdly benign band a0, r0 = 1 m (m)", a0(1.0), 1e-21, "FR '1e14 lP ~ 1e-21 m'", 0.6)
show("  r0 = 1 ly (m)", a0(ly), 0.2e-15, "FR '2e19 lP ~ 0.2 fm'", None)
show("  r0 = 1e5 ly (m)", a0(1e5 * ly), 1e-14, "FR '1e21 lP ~ 1e-14 m'", 0.6)
show("  FR's (r0/(8 f^4 lP))^(1/3) lP, r0 = 1 m (m)", (1 / (8 * f_**4 * lP))**(1 / 3) * lP)
show("fields N for r0 = 1 m: (2 f^2 r0/lP)^2", (2 * f_**2 / lP)**2, 1e62, "FR '1e62 fields'", 0.6)

print("\n=== 9. Maldacena & Milekhin 2021 'humanly traversable' ===")
re = c * np.sqrt(0.5 / (20 * g0))
show("r_e from 20 g tidal limit, 0.5 m body (m)", re, 1.5e7, "MM Eq.remin", 0.02)
show("  as light time (s)", re / c, 0.05, "MM", 0.02)
R5 = 5e-5
L = np.log(re / R5)
show("c2 = 4 pi R5 r_e/(lP^2 sqrt(log))", 4 * np.pi * R5 * re / (lP**2 * np.sqrt(L)), 7e72, "MM Eq.Params", 0.05)
ell = 4 / np.pi * re**2 / R5 * np.sqrt(L)
show("ell = (4/pi)(r_e^2/R5) sqrt(log) (ly)", ell / ly, 3e3, "MM", 0.05)
show("gamma = ell/r_e", ell / re, 2e12, "MM", 0.05)
show("E_bin = -(pi^2/16) R5^2 c^2/(G r_e log) (kg)", -(np.pi**2 / 16) * R5**2 * C2G / (re * L), -5e9, "MM", 0.05)
show("proper traversal time pi r_e/c (s)", np.pi * re / c)
show("mass of extremal BH with r_e = R5: r_e c^2/G (kg)", R5 * C2G, 1e23, "MM '~1e23 kg'", None)
show("mass of the human-sized mouth: r_e c^2/G (M_sun)", re * C2G / Msun)

print("\n=== 10. Quantum inequalities (Ford-Roman; Fewster-Eveson) ===")
lor = quad(lambda t: (2 * t * t - 1)**2 / (t * t + 1)**5, -np.inf, np.inf)[0] / np.pi    # f = sqrt(Lorentzian)
show("FE bound (1/16pi^2) int f''^2, Lorentzian weight, tau0=1", lor / (16 * np.pi**2), 3 / (32 * np.pi**2), "FR C", 1e-9)
gau = quad(lambda t: ((t * t - 1) * np.exp(-t * t / 2))**2, -np.inf, np.inf)[0] / np.sqrt(np.pi)
show("FE bound, Gaussian weight exp(-t^2)/sqrt(pi), tau0=1", gau / (16 * np.pi**2), 3 / (64 * np.pi**2), "this check", 1e-9)
show("Lorentzian QI at tau0 = 1 ns (J/m^3)", -3 * hbar * c / (32 * np.pi**2 * (c * 1e-9)**4), -3.72e-26, "library test", 5e-3)

print("\n=== 11. Warp Factory Alcubierre test case (Helmerich et al. 2024 Fig.) ===")
rr = np.linspace(250, 350, 100001)
pk = C4G / (32 * np.pi) * 0.1**2 * np.max(df_alc(rr, 300, 0.015)**2)
show("peak |T00| v=0.1, R=300, sigma=0.015 (J/m^3)", pk, 2.7e36, "WF Fig. colour bar ~ -2.7e36", None)
show("  ring at r (m)", rr[np.argmax(df_alc(rr, 300, 0.015)**2)])
show("  total Eulerian energy (kg)", -(0.01 / 12) * I_alc(300, 0.015) * C2G)
show("Fuchs shell M/M_J", 4.49e27 / MJ, 2.365, "Fuchs 2024", 2e-3)
show("  2GM/(c^2 R2), R2 = 20 m", 2 * G * 4.49e27 / (c * c * 20))

print("\n=== 12. Extra closed forms used in the notes ===")
show("PF wall bound (3/4) sqrt(3/pi) v/alpha^2 at alpha=0.1 (Planck lengths, per v)", 0.75 * np.sqrt(3 / np.pi) / 0.01, 1e2, "PF Eq. wall_thickness '1e2 v L_P'", None)
show("Capacitance bound (v^2/12) r0 c^2/G, r0 = 100 m, v = 1 (kg)", 100 / 12 * C2G)
show("  in M_sun", 100 / 12 * C2G / Msun)
show("  in Jupiter masses", 100 / 12 * C2G / MJ)
show("tanh thin-wall (v^2/36) R^2 sigma c^2/G, R=100, sigma=2 (M_sun)", 1e4 * 2 / 36 * C2G / Msun)
show("Natario thin-wall (v^2/180) R^4 sigma^3 c^2/G, R=100, sigma=2 (M_sun)", 1e8 * 8 / 180 * C2G / Msun)
show("Fewster-Roman 2005: r0 <~ 1e5 lP/sqrt(1-b0'), 1-b0' = 1e-30 (m)", 1e5 * lP / 1e-15, 1e-15, "FR05 '1e20 lP ~ 1 fermi'", 0.7)
vs10 = brentq(lambda x: np.cosh(0.5 * np.log(abs(1 - x * x))) - 100, 1.5, 1e4)
show("WFM101 '10c' = boost 100 x 0.1c: bubble v_s in White's own formula", vs10)
show("  energy ratio to a v_s = 1 bubble of the same shape (v_s^2)", vs10**2)
show("Alcubierre rho at (0, R, 0), v=1, sigma R = 8 (geometric, exact -1/(2 pi))", -df_alc(1.0, 1.0, 8.0)**2 / (32 * np.pi), -1 / (2 * np.pi), "f'(R) = -sigma/2", 1e-6)
