"""Checks for Puthoff, "Advanced Space Propulsion Based on Vacuum (Spacetime Metric) Engineering"
(JBIS 63, 82-89, 2010; arXiv:1204.2184; DIRD 15) and the PV papers it rests on
(Puthoff 2002, gr-qc/9909037; Puthoff & Little, arXiv:1012.5264 = JBIS 55, 137, 2002;
Puthoff, Davis & Maccone, physics/0403064 = GRG 37, 483, 2005).

Run: python3 -P puthoff_metric_checks.py  (numpy + scipy from user site; -P keeps this folder off sys.path). Output is saved beside it as
puthoff_metric_checks.out.txt. Units SI; CODATA values as in crates/metrology.
"""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

C, G, HBAR = 299_792_458.0, 6.674_30e-11, 1.054_571_817e-34
EPS0, MU0, QE, ME, KB = 8.854_187_812_8e-12, 1.256_637_062_12e-6, 1.602_176_634e-19, 9.109_383_701_5e-31, 1.380_649e-23
MSUN = 1.327_124_4e20 / G
LP = 1.616_255e-35
KE = 1 / (4 * np.pi * EPS0)
rg_per_kg = G / C**2  # m per kg

def head(s):
    print("\n== " + s)

def taylor(f, n=6, r=0.05, N=128):
    """Taylor coefficients of an analytic f at 0 from a Cauchy integral on |x| = r (FFT)."""
    z = r * np.exp(2j * np.pi * np.arange(N) / N)
    a = np.fft.fft(f(z)) / N
    return np.real(a[:n] / r ** np.arange(n))

# ---------------------------------------------------------------- A. Eq. (4) and (5) of the paper
head("A. Eq.(4)/(5): Schwarzschild and Reissner-Nordstrom in harmonic radial coordinate r = R - m")
m, q2 = 1.0, 0.3  # geometric units, m = Gm/c^2, q2 = Q^2 G/(4 pi eps0 c^4)
for r in [1.5, 3.0, 10.0, 100.0]:
    x = m / r
    R = r + m
    g00_4 = (1 - x) / (1 + x)
    g00_s = 1 - 2 * m / R
    grr_4 = 1 / g00_4               # paper: [(1-x)/(1+x)]^-1 dr^2, and dR = dr
    grr_s = 1 / (1 - 2 * m / R)
    gth_4 = (1 + x) ** 2 * r**2     # = R^2
    g00_5 = (1 - x) / (1 + x) + q2 / (r**2 * (1 + x) ** 2)
    g00_rn = 1 - 2 * m / R + q2 / R**2
    print(f"r={r:6.1f}: g00 eq4-Schw {g00_4-g00_s:+.1e}  grr {grr_4-grr_s:+.1e}  R^2 {gth_4-R**2:+.1e}  eq5-RN {g00_5-g00_rn:+.1e}")

# ---------------------------------------------------------------- B. Table 1 in PV variables
head("B. Table 1 with the PV metric g00 = 1/K, g11 = -K (Puthoff 2002 eq. 19)")
for K in [4.0, 1.0 + 1e-9, 0.5, 0.05]:
    g00, g11 = 1 / K, -K
    t = 1 / np.sqrt(g00); w = np.sqrt(g00); E = np.sqrt(g00); L = 1 / np.sqrt(-g11)
    vL = np.sqrt(g00 / -g11); n = np.sqrt(-g11 / g00); mfac = -g11 / np.sqrt(g00)
    ok = np.allclose([t, w, E, L, vL, n, mfac], [np.sqrt(K), K**-0.5, K**-0.5, K**-0.5, 1 / K, K, K**1.5])
    print(f"K={K:<10g} dt x{t:.4g} omega x{w:.4g} E x{E:.4g} dr x{L:.4g} vL/c {vL:.4g} n {n:.4g} m x{mfac:.4g}  == PV-2002 eqs 7-13: {ok}")

# Table 1 mass entry = coefficient of v^2/2 in L = -m c^2 sqrt(g00 - (-g11) v^2/c^2) (slow motion)
for g00, g11 in [(0.7, -1.6), (1.3, -0.8)]:
    Lag = lambda v: -np.sqrt(g00 - (-g11) * v**2)  # m = c = 1
    h = 1e-4
    coef = (Lag(h) - 2 * Lag(0) + Lag(-h)) / h**2
    print(f"g00={g00}, g11={g11}: d2L/dv2 = {coef:.6f}  vs  -g11/sqrt(g00) = {-g11/np.sqrt(g00):.6f}")
# The same coordinate bookkeeping gives local/coordinate acceleration sqrt(-g11)/g00 (= K^1.5 in PV).

# ---------------------------------------------------------------- C. Series: GR (isotropic) vs PV
head("C. Taylor coefficients in x = GM/(r c^2), isotropic coordinates")
gr = {
    "g00": lambda x: ((1 - x / 2) / (1 + x / 2)) ** 2,
    "-gii": lambda x: (1 + x / 2) ** 4,
    "n = sqrt(-gii/g00)": lambda x: (1 + x / 2) ** 3 / (1 - x / 2),
    "mass -gii/sqrt(g00)": lambda x: (1 + x / 2) ** 5 / (1 - x / 2),
    "clock sqrt(g00)": lambda x: (1 - x / 2) / (1 + x / 2),
}
pv = {
    "g00": lambda x: np.exp(-2 * x),
    "-gii": lambda x: np.exp(2 * x),
    "n = sqrt(-gii/g00)": lambda x: np.exp(2 * x),
    "mass -gii/sqrt(g00)": lambda x: np.exp(3 * x),
    "clock sqrt(g00)": lambda x: np.exp(-x),
}
for k in gr:
    a, b = taylor(gr[k], 5), taylor(pv[k], 5)
    print(f"{k:22s} GR {np.round(a, 6)}  PV {np.round(b, 6)}")
print("=> PPN beta = gamma = 1 in both; first difference: -gii at x^2 (GR 3/2, PV 2), g00 at x^3 (GR -3/2, PV -4/3)")

# ---------------------------------------------------------------- D. Light bending from n(r)
head("D. Deflection by an isotropic index n(u), u = M/r; alpha = 2 int_0^u0 b du/sqrt(n^2 - b^2 u^2) - pi")
def deflect(n, b):
    u0 = brentq(lambda u: n(u) - b * u, 0.5 / b, 1.5 / b, xtol=1e-300, rtol=1e-15)  # b >= 100 here
    f = lambda s: 2 * u0 * s * b / np.sqrt(max(n(u0 * (1 - s * s)) ** 2 - (b * u0 * (1 - s * s)) ** 2, 1e-300))
    val, _ = quad(f, 0, 1, epsabs=1e-15, epsrel=1e-14, limit=200)
    return 2 * val - np.pi
nGR = lambda u: (1 + u / 2) ** 3 / (1 - u / 2)
nPV = lambda u: np.exp(2 * u)
eps = np.array([2e-3, 3e-3, 4e-3, 5e-3, 6e-3, 8e-3, 1e-2])
for name, n in [("GR", nGR), ("PV", nPV)]:
    resid = np.array([(deflect(n, 1 / e) - 4 * e) / e**2 for e in eps])
    c4, c3, c2 = np.polyfit(eps, resid, 2)
    c2x = 15 * np.pi / 4 if name == "GR" else 4 * np.pi
    c4b, c3b = np.polyfit(eps, (resid - c2x) / eps, 1)
    print(f"{name}: c2 = {c2:.6f} (15pi/4 = {15*np.pi/4:.6f}, 4pi = {4*np.pi:.6f}); c3 = {c3b:.3f} given exact c2 (128/3 = {128/3:.3f})")
print("Epstein-Shapiro-type check: c2 = pi(2 + 2gamma - beta + 3 delta/4) with delta = (2/3)*(x^2 coeff of -gii) (my check, not from a source):"
      f" GR {np.pi*(2+2-1+0.75*1):.6f}, PV {np.pi*(2+2-1+0.75*4/3):.6f}")
arcsec = 180 / np.pi * 3600
x_limb = G * MSUN / (C**2 * 6.957e8)
print(f"Solar limb: first order {4*x_limb*arcsec:.4f} arcsec; PV-GR second-order difference {(np.pi/4)*x_limb**2*arcsec*1e6:.3f} microarcsec")

# ---------------------------------------------------------------- E. Strong field: photon sphere, throat, shadow
head("E. Strong field (geometric units M = 1)")
bc_pv, bc_gr = 2 * np.e, 3 * np.sqrt(3)
print(f"critical impact parameter: PV 2e = {bc_pv:.5f}, GR 3sqrt3 = {bc_gr:.5f}; PV shadow deviation delta = {bc_pv/bc_gr-1:+.4f}")
r_throat = 1.0  # min of r*exp(1/r)
print(f"PV areal radius r e^(1/r) minimum at r = {r_throat}: e*M = {np.e:.5f} M -> {np.e*G*MSUN/C**2:.1f} m for 1 Msun (no horizon; Boonserm et al. 2018)")
print(f"PV photon sphere: isotropic r = 2M, areal {2*np.exp(0.5):.4f} M; GR areal 3M")
# Surface redshift at a given AREAL radius R (observable with radius + redshift): GR sqrt(1-2y), PV exp(-m/r) with R = r exp(m/r)
for y in [2.12e-6, 0.1723]:  # Sun's surface; illustrative neutron star 1.4 Msun at R = 12 km
    r = brentq(lambda r: r * np.exp(1 / r) - 1 / y, 1.0 + 1e-12, 1e9)
    gr_, pv_ = np.sqrt(1 - 2 * y), np.exp(-1 / r)
    print(f"GM/(R c^2) = {y}: sqrt(g00) GR {gr_:.10f}, PV {pv_:.10f}, PV/GR - 1 = {pv_/gr_-1:+.2e} (series -y^3/6 = {-y**3/6:+.2e})")

# ---------------------------------------------------------------- F. Charged PV solution vs RN
head("F. PV charged solution (Puthoff 2002 eq. 65; JBIS 2002 eq. B6) and K < 1")
def sqrtK(r, a, b):
    if a > b:
        k = np.sqrt(a * a - b * b); return np.cosh(k / r) + a / k * np.sinh(k / r)
    k = np.sqrt(b * b - a * a); return np.cos(k / r) + a / k * np.sin(k / r)
for a, b in [(1.0, 0.5), (1.0, 2.0)]:
    r, h = 7.0, 1e-3
    y = lambda r: sqrtK(r, a, b)
    d1 = (y(r + h) - y(r - h)) / (2 * h); d2 = (y(r + h) - 2 * y(r) + y(r - h)) / h**2
    lhs = d2 + 2 / r * d1; rhs = (d1**2 - b * b / r**4) / y(r)
    xs = 1e-3; r_big = a / xs
    g00 = 1 / sqrtK(r_big, a, b) ** 2
    pred = 1 - 2 * xs + (2 + b * b / (a * a)) * xs**2
    print(f"a={a}, b={b}: ODE (64b) residual {lhs-rhs:+.1e}; weak-field g00 vs eq.(70) {g00-pred:+.1e} (O(x^3) = {xs**3:.0e})")
a, b = 1.0, 2.0
k = np.sqrt(b * b - a * a)
r0 = k / (np.pi - np.arctan(k / a))  # first zero of cos(k/r) + (a/k) sin(k/r) coming in from infinity
print(f"b/a = 2: sqrtK -> 0 (g00 = 1/K -> infinity; areal radius r sqrt(K) -> 0) at isotropic r = {r0:.4f} a")
r_K1 = brentq(lambda r: sqrtK(r, a, b) - 1, r0 * 1.0001, 100)
print(f"   K < 1 (g00 > 1) inside r = {r_K1:.4f} a (areal = isotropic there); GR RN: g00 > 1 inside areal R = b^2/(2a) = {b*b/(2*a):.3f} a, singular at R = 0")
qm = np.sqrt(4 * np.pi * EPS0 * G)
print(f"b > a  <=>  Q/M > sqrt(4 pi eps0 G) = {qm:.4e} C/kg (gem::extremal_charge_to_mass)")
re = KE * QE**2 / (ME * C**2)
print(f"RN g00 > 1 for areal R < k Q^2/(2 M c^2): electron {re/2:.4e} m = r_e/2 (r_e = {re:.4e} m)")
Mcraft, Rcraft = 1e4, 10.0
Qcraft = np.sqrt(2 * Mcraft * C**2 * Rcraft / KE)
print(f"craft M = 1e4 kg: K < 1 out to 10 m needs k Q^2/(2 r) = M c^2 -> Q = {Qcraft:.3e} C, E(10 m) = {KE*Qcraft/Rcraft**2:.3e} V/m;"
      f" Schwinger field {ME**2*C**3/(QE*HBAR):.3e} V/m")
print(f"   but the metric change there is only |g00 - 1| ~ 2GM/(rc^2) = {2*G*Mcraft/(Rcraft*C**2):.2e}")

# ---------------------------------------------------------------- G. What it takes to make g00 > 1 (Table 1 right column)
head("G. Negative mass needed for a static region of areal radius R = 10 m with redshift factor sqrt(g00) > 1")
R = 10.0
def m_gr(d):  # d = g00 - 1; exterior Schwarzschild with M < 0: g00 = 1 + 2G|M|/(R c^2)
    return d * R * C**2 / (2 * G)
def m_pv(d):  # K = exp(-2G|M|/(r c^2)), areal R = r sqrt(K)
    y = np.log1p(d) / 2; r = R * np.exp(y); return y * r * C**2 / G
wien = 2.897_771_955e-3
cases = [
    ("optical-clock detectable, g00-1 = 1e-18", 1e-18),
    ("1 ppm clock speed-up", 2e-6),
    ("plants age 10% faster (sqrt g00 = 1.1)", 1.1**2 - 1),
    ("'20 min vs a few (5) min' (sqrt g00 = 4)", 15.0),
    ("'20 min vs a few (3) min' (sqrt g00 = 6.67)", (20 / 3) ** 2 - 1),
    (f"300 K Wien peak {wien/300*1e6:.2f} um -> 550 nm (sqrt g00 = {wien/300/550e-9:.1f})", (wien / 300 / 550e-9) ** 2 - 1),
    ("20 kHz sound -> 20 Hz (sqrt g00 = 1000)", 1e6 - 1),
]
for label, d in cases:
    g00 = 1 + d
    mg, mp = m_gr(d), m_pv(d)
    rho = mg * C**2 / (4 / 3 * np.pi * R**3)
    gam = np.sqrt(g00)  # gamma_inf - 1 ~ d/2 for small d
    acc = C**2 * d / (2 * R) / np.sqrt(g00)  # outward proper acceleration of a static observer at R (GR)
    print(f"{label:55s} GR |M| = {mg:.3e} kg ({mg/MSUN:.2e} Msun), PV |M| = {mp:.3e} kg; "
          f"mean |rho c^2| = {rho:.2e} J/m3; entry needs gamma_inf >= {gam:.6g}; surface repulsion {acc:.2e} m/s2")
print("(Casimir energy density at 1 um: %.2e J/m3; nuclear matter ~2.3e17 kg/m3 -> %.1e J/m3)" % (np.pi**2 * HBAR * C / (720 * 1e-24), 2.3e17 * C**2))
print("Section 3.1 example (10 s -> 5 s, g00 = 0.25): GR areal R = 8/3 GM/c^2 =", f"{8/3:.4f}", "; PV isotropic r = 2/ln4 =", f"{2/np.log(4):.4f}", ", areal", f"{2/np.log(4)*2:.4f} GM/c^2")
for d in [1e-18, 1e-6]:
    print(f"monopole: g00 shift {d:.0e} at R = 10 m implies exterior gravity c^2 d/(2R) = {C**2*d/(2*R):.2e} m/s2;"
          f" confined (zero net mass) needs a +/- mass double layer of dipole density d c^2/(8 pi G) = {d*C**2/(8*np.pi*G):.2e} kg/m")
print(f"a 1e-8 m/s2 (1 microGal) gravimeter at 10 m resolves a monopole g00 shift of 2gR/c^2 = {2*1e-8*R/C**2:.1e}")
print(f"Weak-field scale: u = delta c^4/(4 pi G L^2) for central |g00-1| = delta, uniform sphere L: delta=1, L=10 m -> {C**4/(4*np.pi*G*100):.2e} J/m3; c^4/G = {C**4/G:.4e} J/m")

# ---------------------------------------------------------------- H. Levi-Civita effect (paper ref. [19])
head("H. Levi-Civita effect in PV (Puthoff, Davis & Maccone eqs 6, 8, 13, 15, B-2)")
print(f"B-2: c^2/sqrt(4 pi G/mu0) = {C**2/np.sqrt(4*np.pi*G/MU0):.4e} T m (paper: 3.484e18)")
for B, L in [(1.0, 1.0), (45.0, 1.0), (1e9, 1e4)]:
    dK = 2 * np.pi * G * B**2 / (MU0 * C**4) * (L / 2) ** 2
    extra = L / C * 2 * np.pi * G * B**2 * L**2 / (3 * MU0 * C**4)
    print(f"B = {B:g} T over L = {L:g} m: K(0) - 1 = {dK:.3e}; extra light transit time {extra:.3e} s")
for E, L in [(3e6, 1.0), (ME**2 * C**3 / (QE * HBAR), 1.0)]:
    dK = 2 * np.pi * EPS0 * G * E**2 / C**4 * (L / 2) ** 2
    print(f"E = {E:.3g} V/m over L = {L:g} m: K(0) - 1 = {dK:.3e}")
print("Sign: K > 1 (light slower) -- EM energy gravitates like mass; it cannot make K < 1.")

# ---------------------------------------------------------------- I. Pfenning & Ford numbers (paper ref. [22] context)
head("I. Alcubierre total energy, Pfenning & Ford eq. (28), R = 100 m, wall Delta = 100 v L_P (eq. 23)")
Rb, v = 100.0, 1.0
Dl = 100 * v * LP
E_len = -(v**2 / 12) * (Rb**2 / Dl + Dl / 12)  # metres (G = c = 1)
print(f"E = {E_len:.3e} m = {E_len/LP:.3e} L_P = {E_len*C**2/G:.3e} kg = {E_len*C**2/G*1e3:.3e} g (P&F: -6.2e70 L_P ~ -6.2e65 g)")
print(f"Delta = 1 m: {-(1/12)*(Rb**2/1+1/12)*C**2/G/MSUN:.3f} Msun (P&F: 'a quarter of a solar mass')")

# ---------------------------------------------------------------- J. Small checks quoted in the notes
head("J. Small checks")
print(f"Levi-Civita B-2 with c = 3e8: {9e16/np.sqrt(4*np.pi*G/MU0):.4e} T m (paper 3.484e18)")
for gam in [1.1, 4.0, 20 / 3, 17.5623, 1000.0]:
    print(f"gamma_inf >= {gam}: v >= {np.sqrt(1-1/gam**2):.6f} c")
print(f"needed mean |rho c^2| for a 1e-18 shift at 10 m over Casimir |u| at 1 um: {1.44e23/(np.pi**2*HBAR*C/(720*1e-24)):.2e}")
print(f"GP-B frame dragging: 37.2/7.2 = {37.2/7.2:.2f} sigma from zero")
print(f"Sun surface redshift -GM/(R c^2) = {-G*MSUN/(6.957e8*C**2):.4e}; Pound-Rebka g h/c^2 (g = 9.81, h = 22.5 m) = {9.81*22.5/C**2:.4e}")
print(f"negative-mass lensing for the 1 ppm craft (|M| = 1.347e22 kg) at b = 10 m: {4*G*1.347e22/(C**2*10)*180/np.pi*3600:.2f} arcsec (away from the craft)")
