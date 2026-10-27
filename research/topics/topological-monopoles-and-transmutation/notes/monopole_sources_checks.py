"""Recompute the primary-source formulas for magnetics::monopole (issue #24). Stdlib only.
Run: python3 -I monopole_sources_checks.py   (exits 1 on any failure)"""
import math, sys

FAIL = []
def check(name, got, want, rtol):
    ok = abs(got - want) <= rtol * abs(want)
    print(f"{'ok  ' if ok else 'FAIL'} {name}: {got:.6g} (want {want:.6g}, rtol {rtol:g})")
    if not ok: FAIL.append(name)

def bisect(f, lo, hi, n=200):
    flo = f(lo)
    for _ in range(n):
        mid = 0.5 * (lo + hi); fm = f(mid)
        if (fm > 0) == (flo > 0): lo, flo = mid, fm
        else: hi = mid
    return 0.5 * (lo + hi)

print("== 1. Monopole-nucleus binding ==")
# Olaussen & Sollie, NPB 255 (1985) 465, eqs 4.16a,b and 4.17 (hbar = c = 1):
# X = g_M |q|, q = Z e g = nZ/2, g_M = A mu/(S Z) (mu in nuclear magnetons).
def os_lams(Z, A, S, mu, n=1):
    q = abs(n) * Z / 2; X = A * mu / (S * Z) * q
    l0 = S * X - q
    r = math.sqrt((q + 1)**2 + (2 * S - 1 - q) * X + X * X / 4)
    lp = (S - 0.5) * X - 2 * q - 1 + r
    lm = (S - 0.5) * X - 2 * q - 1 - r
    return q, X, l0, lp, lm
# Olaussen-Sollie Table (|eg| = 1/2): S, mu, |q|, X, lambda0, lambda+, lambda-
os_table = {"23Na": (11, 23, 1.5, 2.2175, 5.5, 17.0, 20.0, 12.4, -2.4),
            "27Al": (13, 27, 2.5, 3.6414, 6.5, 19.7, 42.8, 35.6, 15.2),
            "55Mn": (25, 55, 2.5, 3.444, 12.5, 37.9, 82.3, 64.6, 35.0),
            "85Rb": (37, 85, 2.5, 1.3524, 18.5, 23.0, 39.0, 21.4, -5.4),
            "135Ba": (56, 135, 1.5, 0.8365, 28.0, 37.6, 28.4, -4.7, -34.1),
            "138La": (57, 138, 5.0, 3.707, 28.5, 51.2, 227.5, 195.4, 149.4)}
for k, (Z, A, S, mu, q_, X_, l0_, lp_, lm_) in os_table.items():
    q, X, l0, lp, lm = os_lams(Z, A, S, mu)
    check(f"O-S {k} X", X, X_, 3e-3)
    check(f"O-S {k} lambda0", l0, l0_, 4e-3)
    check(f"O-S {k} lambda+", lp, lp_, max(0.1 / abs(lp_), 3e-3))
    check(f"O-S {k} lambda-", lm, lm_, max(0.1 / abs(lm_), 3e-3))

# Milton gamma-hat = (A/Z) mu/mu_N = S g_M  => X = gamma_hat * l / S, l = |q|
def lam_pm_ghat(gh, S, l):
    X = gh * l / S; r = math.sqrt((l + 1)**2 + (2 * S - 1 - l) * X + X * X / 4)
    return (S - 0.5) * X - 2 * l - 1 + r, (S - 0.5) * X - 2 * l - 1 - r
def thresholds(S, l):
    gm = bisect(lambda g: lam_pm_ghat(g, S, l)[0] - 0.25, 0.0, 50.0)
    gp = bisect(lambda g: lam_pm_ghat(g, S, l)[1] - 0.25, 0.0, 50.0)
    return gm, gp
def gc32(l): s = math.sqrt(33 + 32 * l); return (3 / (4 * l) * (6 + 4 * l - s), 3 / (4 * l) * (6 + 4 * l + s))
def gc52(l): s = math.sqrt(1161 + 1296 * l + 64 * l * l); return ((36 + 28 * l - s) / (12 * l), (36 + 28 * l + s) / (12 * l))
def gc1(l): return 3 / (4 * l) * (3 + 16 * l + 16 * l * l) / (9 + 4 * l)
for l in (2, 4, 6.5, 20):
    a, b = thresholds(1.5, l); c, d = gc32(l)
    check(f"Milton 9.12 vs O-S 4.16b root, S=3/2 l={l}: gc-", c, a, 1e-9)
    check(f"Milton 9.12 vs O-S 4.16b root, S=3/2 l={l}: gc+", d, b, 1e-9)
    a, b = thresholds(2.5, l); c, d = gc52(l)
    check(f"Milton 9.13 vs O-S 4.16b root, S=5/2 l={l}: gc-", c, a, 1e-9)
    check(f"Milton 9.13 vs O-S 4.16b root, S=5/2 l={l}: gc+", d, b, 1e-9)
    a = bisect(lambda g: lam_pm_ghat(g, 1, l)[0] - 0.25, 0.0, 50.0)
    check(f"Milton 9.10 vs O-S 4.16b root, S=1 l={l}", gc1(l), a, 1e-9)
check("9.12 at l=2 (Milton 1.557)", gc32(2)[0], 1.557, 3e-4); check("9.12 at l=2 (Milton 8.943)", gc32(2)[1], 8.943, 3e-4)
check("9.13 at l=6.5 (Milton 1.374)", gc52(6.5)[0], 1.374, 3e-4); check("9.13 at l=6.5 (Milton 4.216)", gc52(6.5)[1], 4.216, 3e-4)
print(f"     limits l->inf: 9.12 -> {gc32(1e12)[0]:.4f}, {gc32(1e12)[1]:.4f}; 9.13 -> {gc52(1e12)[0]:.4f}, {gc52(1e12)[1]:.4f}")
# lowest state lambda0 > 1/4  <=>  gamma_hat > 1 + 1/(4l)  (Milton 9.7), any S
for S in (0.5, 1.5, 2.5):
    l = 6.5; g = 1 + 1 / (4 * l)
    check(f"lambda0 = S X - l at gamma_hat = 1 + 1/4l is 1/4 (S={S})", S * (g * l / S) - l, 0.25, 1e-12)

# Bracci & Fiorentini NPB 232 (1984) 236, spin 1/2, n = 1: mu = e K S/m (m = A m_p) => K = A mu~;
# lambda_I = |q| - K/2 (2.19), bound if lambda_I < -1/4; lambda_II = (j^2-q^2) - sqrt(j^2 - |q| K + K^2/4) (2.27), j = J + 1/2
def bf_lI(Z, A, mu): return Z / 2 - A * mu / 2
def bf_lII(Z, A, mu, j): q = Z / 2; K = A * mu; return (j * j - q * q) - math.sqrt(j * j - q * K + K * K / 4)
for nm, Z, A, mu, want in (("1H", 1, 1, 2.79, 0.90), ("13C", 6, 13, 0.70, 1.57), ("19F", 9, 19, 2.63, 20.47), ("31P", 15, 31, 1.13, 10.04), ("3H (table prints 2.97: typo, lambda_II row confirms K)", 1, 3, 2.98, 3.97)):
    check(f"B-F Table I -lambda_I {nm}", -bf_lI(Z, A, mu), want, 0.02)
for nm, Z, A, mu, want in (("3He", 2, 3, -2.13, 1.54), ("113Cd", 48, 113, -0.62, 10.48), ("3H", 1, 3, 2.98, 2.21)):
    q = Z / 2; check(f"B-F Table I -lambda_II {nm} (j = q+1)", -bf_lII(Z, A, mu, q + 1), want, 0.02)
# Milton 9.8 |kappa_hat| > (J^2+J-l^2)/l reproduces the lambda_II = -1/4 threshold (n = 1, K = 2 l gamma_hat)
for l, J in ((1, 1.5), (24, 24.5), (5, 7.5)):
    j = J + 0.5; f = lambda gh: (j*j - l*l) - math.sqrt(j*j - l*(2*l*gh) + (2*l*gh)**2/4) + 0.25
    kc = (J*J + J - l*l) / l
    check(f"Milton 9.8 kappa_c vs B-F 2.27 root (+), l={l} J={J}", bisect(f, 1.0, 1e3) - 1, kc, 1e-9)
    check(f"Milton 9.8 kappa_c vs B-F 2.27 root (-), l={l} J={J}", 1 - bisect(f, -1e3, 1.0), kc, 1e-9)
# Milton 9.1 (hbar = 1, Gaussian): E = -(mu q^2/2)[N + 1/2 + ((j+1/2)^2 - m'^2)^(1/2)]^-2
Ry = 13.605693122990
def E91(Z, N, j, m): return -Z * Z * Ry / (N + 0.5 + math.sqrt((j + 0.5)**2 - m * m))**2
check("9.1, m'=0, j=0: hydrogen", E91(1, 0, 0, 0), -Ry, 1e-12)
check("9.1, m'=1/2, j=1/2, N=0 [eV]", E91(1, 0, 0.5, 0.5), -7.29134, 1e-5)

print("\n== 2. Stopping power ==")
NA = 6.02214076e23; alpha = 7.2973525643e-3; c_cm = 2.99792458e10; hbar = 1.054571817e-34; me = 9.1093837139e-31
e2_eVcm = 1.43996448e-7; mec2_eV = 0.51099895069e6
def vF_cm_s(rho, Z, A): Ne = rho / A * NA * Z * 1e6; return hbar / me * (3 * math.pi**2 * Ne)**(1/3) * 100
check("v_F Si 4.3 g/cm3 (P&S 3.9e8)", vF_cm_s(4.3, 14, 28.0855), 3.9e8, 0.01)
def slow(beta, rho, Z, A, n=1, spin=1.37):  # Gaussian eq 12, MeV/cm
    Ne = rho / A * NA * Z; vF = vF_cm_s(rho, Z, A)
    g2e2 = (n / (2 * alpha))**2 * e2_eVcm**2
    return spin * 2 * math.pi * Ne * g2e2 * beta / (mec2_eV * vF / c_cm) * (math.log(2 * vF / (alpha * c_cm)) - 0.5) / 1e6
for rho in (2.33, 4.3):
    print(f"     eq12 Si rho={rho}: x1 {slow(1e-3, rho, 14, 28.0855, spin=1):.2f}, x1.37 {slow(1e-3, rho, 14, 28.0855):.2f} MeV/cm at beta=1e-3; Milton 8.11 (A-K eq 60): 45 MeV/cm")
check("eq12 x1.37 Si 2.33 beta=1e-3 [MeV/cm] (repo 61)", slow(1e-3, 2.33, 14, 28.0855), 61, 0.02)
# fast regime, Ahlen 1978 as printed by Milton 8.10 / P&S 11; PDG Bethe K = 4 pi N_A r_e^2 m_e c^2
Kb = 0.307075  # MeV cm^2/mol
def fast(beta, rho, ZA, I_eV, n=1, delta=0.0):
    KYG = 0.406 if n == 1 else 0.346; B = {1: 0.248, 2: 0.672, 3: 1.022, 6: 1.685}[n]
    g = 1 / math.sqrt(1 - beta * beta)
    L = math.log(2 * mec2_eV * beta**2 * g**2 / I_eV) + KYG / 2 - delta / 2 - 0.5 - B
    return Kb * ZA * rho * (n / (2 * alpha))**2 * L  # MeV/cm
print(f"     (g_D/e)^2 = {(1/(2*alpha))**2:.2f}")
for b in (0.1, 0.2, 0.5, 0.8):
    print(f"     eq11 Si 2.33 I=173 eV delta=0 beta={b}: {fast(b, 2.33, 14/28.0855, 173.0)/1e3:.4f} GeV/cm")

print("\n== 3. PDG 2025 rev. 94 ==")
h = 6.62607015e-34; e = 1.602176634e-19; mu0 = 1.25663706127e-6; c = 299792458.0
gD = h / (mu0 * e)
check("94.2 Q_M^D = 2pi/e in HL natural units => Q_M/e = 1/(2 alpha) == g_D/(e c)", 2 * math.pi / (4 * math.pi * alpha), gD / (e * c), 1e-9)
GeV = 1e9 * e
E_gain = gD * 3e-10 * 3.0857e19  # 3 uG over 1 kpc
print(f"     g_D x 3 uG x 1 kpc = {E_gain/GeV:.4g} GeV; 1e11 GeV over {1e11*GeV/(gD*3e-10)/3.0857e19:.4f} kpc")
def v_mag(M_GeV): return 1.0 if M_GeV <= 1e11 else 1e-3 * math.sqrt(1e17 / M_GeV)
for M in (1e13, 1e17):
    print(f"     94.8: M={M:.0e} GeV v/c={v_mag(M):.3g}, (1/2)Mv^2 = {0.5*M*v_mag(M)**2:.3g} GeV")
def parker(M): return 1e-15 if M <= 1e17 else 1e-15 * M / 1e17
def seed(M): return (M / 1e17 + 3e-6) * 1e-16
for M, wp, ws in ((1e15, 1e-15, 1.0003e-18), (1e17, 1e-15, 1.000003e-16), (1e19, 1e-13, 1.00000003e-14)):
    check(f"94.9 Parker M={M:.0e}", parker(M), wp, 1e-12); check(f"94.10 seed M={M:.0e}", seed(M), ws, 1e-12)
print(f"     seed bound floor (M -> 0): {seed(0):.3g}; seed = Parker/10 at M = 1e17: {seed(1e17)/parker(1e17):.5f}")
# Super-K 2012 eq 9 from eqs 7, 8 with eps = 1
R_sun = 7.0e10; t_sun = 4.6e9 * 3.15576e7; NM = 7.1e17
def F_sk(beta): return NM / (math.pi * R_sun**2 * (1 + (2e-3 / beta)**2) * 4 * math.pi * t_sun)
print(f"     Super-K eq 8 small-beta form (beta_esc/beta)^2: F(1e-3) = {NM/(math.pi*R_sun**2*(2e-3/1e-3)**2*4*math.pi*t_sun):.3g}")
print(f"     Super-K eqs 7-8 (eps=1): F(1e-3) = {F_sk(1e-3):.3g}; F(1e-4) = {F_sk(1e-4):.3g}; F(1e-4)/F(1e-3) = {F_sk(1e-4)/F_sk(1e-3):.4f} (beta^2 law: 0.01)")
def superk(beta, sigma0_mb=1.0, fpi=0.5): return 6.3e-24 * (beta / 1e-3)**2 / (sigma0_mb * fpi / 0.5)
check("Super-K eq 9 (1e-4, 1 mb, 0.5)", superk(1e-4), 6.3e-26, 1e-12)
# IceCube eqs 4-5
def sigma_cat(s0, beta, beta0=0.0, gam=0.0): return s0 / beta * ((beta / beta0)**gam if beta < beta0 else 1.0)
def lam_cat(s0, beta, n, beta0=0.0, gam=0.0): return 1 / (sigma_cat(s0, beta, beta0, gam) * n)
n_H = 2 * 0.917 / 18.015 * NA
print(f"     ice H density {n_H:.4g}/cm3; lambda for sigma_cat 1.7e-22: {1/(1.7e-22*n_H)*10:.3f} mm; 1.7e-23: {1/(1.7e-23*n_H):.3f} cm")
check("IceCube: sigma0 = 1.7e-25 cm2 at beta=1e-3 -> sigma_cat 1.7e-22", sigma_cat(1.7e-25, 1e-3), 1.7e-22, 1e-12)

print("\n== 4. Lambda_MSbar from alpha_s (4-loop, PDG eq 9.5 truncated at b3) ==")
z3 = 1.2020569031595942
def bcoef(nf):
    b0 = (33 - 2 * nf) / (12 * math.pi); b1 = (153 - 19 * nf) / (24 * math.pi**2)
    b2 = (2857 - 5033 / 9 * nf + 325 / 27 * nf**2) / (128 * math.pi**3)
    b3 = (149753 / 6 + 3564 * z3 - (1078361 / 162 + 6508 / 27 * z3) * nf + (50065 / 162 + 6472 / 81 * z3) * nf**2 + 1093 / 729 * nf**3) / (256 * math.pi**4)
    return b0, b1, b2, b3
def as4(mu, Lam, nf):
    b0, b1, b2, b3 = bcoef(nf); t = math.log(mu * mu / (Lam * Lam)); l = math.log(t)
    return 1 / (b0 * t) * (1 - b1 * l / (b0**2 * t) + (b1**2 * (l * l - l - 1) + b0 * b2) / (b0**4 * t * t)
                          + (b1**3 * (-2 * l**3 + 5 * l * l + 4 * l - 1) - 6 * b0 * b2 * b1 * l + b0**2 * b3) / (2 * b0**6 * t**3))
mZ = 91.1876
L5_2018 = bisect(lambda L: as4(mZ, L, 5) - 0.1181, 0.05, 0.5)
check("Lambda5 for alpha_s=0.1181 (PDG 2018 eq 9.24b: 210 MeV)", L5_2018 * 1e3, 210, 0.01)
L5_2025 = bisect(lambda L: as4(91.1880, L, 5) - 0.1180, 0.05, 0.5)
print(f"     Lambda5 for PDG 2025 alpha_s = 0.1180: {L5_2025*1e3:.1f} MeV (derived, not printed by PDG 2025)")
print("\nFAILED: " + ", ".join(FAIL) if FAIL else "\nall checks passed")

print("\n== SI forms for Rust (alpha-free: g^2 e^2 -> n^2 (hbar c)^2 / 4) ==")
hbarc = hbar * c; me_c2 = me * c * c
def slow_si(beta, n_e, v_f, n=1, spin=1.37):  # J/m; n_e in m^-3, v_f in m/s
    return spin * math.pi / 2 * n_e * n * n * hbarc**2 * beta / (me * c * v_f) * (math.log(2 * v_f / (alpha * c)) - 0.5)
def fast_si(beta, n_e, I_J, n=1, delta=0.0):
    KYG = 0.406 if n == 1 else 0.346; B = {1: 0.248, 2: 0.672, 3: 1.022, 6: 1.685}[n]
    g2 = 1 / (1 - beta * beta)
    return math.pi * n_e * n * n * hbarc**2 / me_c2 * (math.log(2 * me_c2 * beta**2 * g2 / I_J) + KYG / 2 - delta / 2 - 0.5 - B)
ne_si = 2.33 / 28.0855 * NA * 14 * 1e6
vf = vF_cm_s(2.33, 14, 28.0855) / 100
MeVcm = 1e6 * e / 1e-2
check("SI slow == Gaussian slow (Si 2.33, beta 1e-3) [MeV/cm]", slow_si(1e-3, ne_si, vf) / MeVcm, slow(1e-3, 2.33, 14, 28.0855), 1e-6)
check("SI fast == Gaussian fast (Si 2.33, beta 0.5) [MeV/cm]", fast_si(0.5, ne_si, 173 * e) / MeVcm, fast(0.5, 2.33, 14 / 28.0855, 173.0), 1e-5)
print(f"     Si 2.33: n_e = {ne_si:.6g} m^-3, v_F = {vf:.6g} m/s; slow(1e-3) = {slow_si(1e-3, ne_si, vf):.6g} J/m; fast(0.5) = {fast_si(0.5, ne_si, 173*e):.6g} J/m")
sys.exit(1 if FAIL else 0)
