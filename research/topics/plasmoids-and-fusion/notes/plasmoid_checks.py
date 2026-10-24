"""Recomputed numbers for plasmoids_frc_reconnection.md. Run: python3 -I plasmoid_checks.py
Standard library only. Every assert compares against a published value (tolerance noted)."""
import math

# CODATA 2018
MU0 = 1.25663706212e-6      # H/m
E = 1.602176634e-19         # C
KB_EV = 8.617333262e-5      # eV/K
MP = 1.67262192369e-27      # kg
MD = 3.3435837724e-27       # kg (deuteron)
U_MEV = 931.49410242        # MeV per u
P_ATM = 101325.0            # Pa


def close(x, ref, rel):
    return abs(x - ref) <= rel * abs(ref)


# ---------- Bessel functions by power series (fine for x < ~10) ----------
def jn(n, x):
    s, k = 0.0, 0
    while k < 60:
        s += (-1) ** k * (x / 2) ** (2 * k + n) / (math.factorial(k) * math.factorial(k + n))
        k += 1
    return s


def bisect(f, a, b, tol=1e-12):
    fa = f(a)
    for _ in range(200):
        m = 0.5 * (a + b)
        fm = f(m)
        if fa * fm <= 0:
            b = m
        else:
            a, fa = m, fm
        if b - a < tol:
            break
    return 0.5 * (a + b)


j01 = bisect(lambda x: jn(0, x), 2.0, 3.0)          # first zero of J0
j11 = bisect(lambda x: jn(1, x), 3.0, 4.5)          # first zero of J1
sph = bisect(lambda x: math.tan(x) - x, 4.3, 4.6)   # first root of tan x = x  (j1(x)=0, spherical Bessel)
assert close(j01, 2.404825557695773, 1e-9) and close(j11, 3.831705970207512, 1e-9)
assert close(sph, 4.493409457909064, 1e-9)

print("== Q1 Taylor states ==")
print(f"J0 zero 2.4048 -> BFM field reversal at mu*a = {j01:.4f}, Theta = mu*a/2 = {j01/2:.4f}")


def bfm(mua):
    """Bessel-function model: Theta = B_theta(a)/<B_z>, F = B_z(a)/<B_z>."""
    return mua / 2, mua * jn(0, mua) / (2 * jn(1, mua))


th, F = bfm(3.11)    # mu*a = 3.11: onset of the m=1 helical Taylor state (arXiv:1509.07305)
print(f"BFM at mu*a = 3.11: Theta = {th:.3f}, F = {F:.3f}")
assert close(th, 1.555, 1e-3)
th0, F0 = bfm(j01)
assert abs(F0) < 1e-9


def spheromak_cyl(h_over_a):
    """Lowest axisymmetric Taylor state in a cylinder radius a, length h: (lam a)^2 = j11^2 + (pi a/h)^2."""
    return math.sqrt(j11 ** 2 + (math.pi / h_over_a) ** 2)


for ha in (1.0, 1.67, 3.0):
    print(f"cylindrical spheromak h/a = {ha}: lambda*a = {spheromak_cyl(ha):.4f}")
print(f"spherical spheromak (flux conserver radius a): lambda*a = {sph:.4f}")
# Taylor state energy in a closed flux conserver: W = lambda*K/(2 mu0)
K, lam = 1e-3, 10.0    # Wb^2, 1/m (example)
print(f"example W = lam K/(2 mu0) for K=1e-3 Wb^2, lam=10/m: {lam*K/(2*MU0):.3e} J")

print("\n== Q2 FRC equilibrium ==")


def beta_avg(xs):
    return 1 - xs ** 2 / 2          # Tuszewski 1988 Eq. (10), Barnes et al.


for xs in (0.3, 0.4, 0.5, 0.6, 0.7):
    print(f"xs = {xs}: <beta> = {beta_avg(xs):.3f}; flux bounds phi/(pi R^2 Be) in "
          f"[{xs**2/2:.3f}, {xs/math.sqrt(2):.3f}]")
assert close(beta_avg(0.5), 0.875, 1e-12)


def s_bounds(xs, S_star):
    """S* = R/rho_ie; Tuszewski 1988: S* xs^2/(4 sqrt2) < s < S* xs/4."""
    return S_star * xs ** 2 / (4 * math.sqrt(2)), S_star * xs / 4


lo, hi = s_bounds(0.5, 10.0)
print(f"s bounds for xs=0.5, R/rho_ie=10: {lo:.3f} < s < {hi:.3f}")
# consistency: the s bounds follow from s = phi/(2 pi rs rho_ie Be) and the flux bounds with rs = sqrt2 R
assert close(lo, (0.5 ** 2 / 2) * 10 / (2 * math.sqrt(2)), 1e-12)
assert close(hi, (0.5 / math.sqrt(2)) * 10 / (2 * math.sqrt(2)), 1e-12)

# Norm (Roche et al. 2025): rs = 0.4 m, length 2 m, flux ~6 mWb, current ~300-350 kA.
# shortcut: Be from a long uniform current sheet with B(axis) ~ -Be, so I/l ~ 2 Be/mu0; order-of-magnitude only.
rs, ell, phi, I = 0.4, 2.0, 6e-3, 350e3
Be = MU0 * I / (2 * ell)
R = rs / math.sqrt(2)
ratio = phi / (math.pi * R ** 2 * Be)
xs_min, xs_max = ratio * math.sqrt(2), math.sqrt(2 * ratio)
print(f"Norm check: Be ~ {Be:.3f} T, phi/(pi R^2 Be) = {ratio:.3f} -> xs in [{xs_min:.2f}, {xs_max:.2f}], "
      f"<beta> in [{beta_avg(xs_max):.2f}, {beta_avg(xs_min):.2f}]")
pM = Be ** 2 / (2 * MU0)
Vcyl = math.pi * rs ** 2 * ell
print(f"Norm check: pM = {pM:.0f} Pa; 1.5*<beta>*pM*V(cylinder, upper bound on volume) ~ "
      f"{1.5*0.9*pM*Vcyl/1e3:.1f} kJ vs published total 9 kJ")

# Adiabatic flux compression of an elongated FRC: T ~ Be^(4/5), n ~ Be^(6/5) (Tuszewski 1988 Sec 3.3)
for T0 in (10e6, 20e6):
    bratio = (150e6 / T0) ** (5 / 4)
    print(f"Helion: {T0/1e6:.0f} -> 150 MdegC by adiabatic flux compression alone needs Be x{bratio:.1f}, "
          f"n x{bratio**1.2:.0f}")

# Trenta separately reported maxima: n=3e22 m^-3, Ti=8 keV, Te=1 keV, B=8 T (Wikipedia)
p_needed = 3e22 * (8e3 + 1e3) * E
print(f"Trenta: n(Ti+Te) = {p_needed:.3e} Pa needs Be >= {math.sqrt(2*MU0*p_needed):.1f} T at beta=1 "
      f"(reported 8 T gives pM = {64/(2*MU0):.3e} Pa)")

print("\n== Q3 reconnection ==")


def spitzer_par(Te_eV, lnL, Z=1.0):
    """NRL formulary: eta_perp = 1.03e-2 Z lnL T^-3/2 ohm cm; eta_par = eta_perp/1.96 -> ohm m."""
    return 1.03e-4 * Z * lnL * Te_eV ** -1.5 / 1.96


def lundquist(L, n, Te, B, lnL, A=1.0):
    vA = B / math.sqrt(MU0 * n * A * MP)
    return MU0 * L * vA / spitzer_par(Te, lnL), vA


# Ji & Daughton 2011 Table 1 rows (L_CS = eps*L, B_R = fraction*B_T as their notes say)
rows = {
    "MRX":        (0.25 * 0.8, 1e19, 10, 0.3 * 0.1, 10, 3e3),
    "corona":     (0.5 * 1e7, 1e15, 200, 0.1 * 2e-2, 20, 1e13),
    "magnetotail": (0.5 * 6e8, 3e5, 600, 2e-8, 30, 4e15),
}
for k, (L, n, Te, B, lnL, Spub) in rows.items():
    S, vA = lundquist(L, n, Te, B, lnL)
    tA = L / vA
    print(f"{k}: vA = {vA:.3e} m/s, tau_A = {tA:.3e} s, S = {S:.2e} (published {Spub:.0e}); "
          f"SP rate S^-1/2 = {S**-0.5:.2e}; tau_SP = {tA*S**0.5:.2e} s; tau at rate 0.01 = {tA/0.01:.2e} s; "
          f"Loureiro N~S^3/8 = {S**0.375:.1e}; gamma*tau_A~S^1/4 = {S**0.25:.1e}")
    assert Spub / 3 < S < Spub * 3, k     # lnL and eta convention differ from the paper; factor-3 band
# Tuszewski 1988 Sec 5.3.2 uses eta_perp/mu0 (m^2/s) = [Te(eV)/100]^-3/2 for flux-decay anomaly factors
eperp = 1.03e-4 * 12 * 100 ** -1.5 / MU0
print(f"Spitzer eta_perp/mu0 at Te=100 eV, lnL=12: {eperp:.2f} m^2/s (Tuszewski shorthand: 1.0)")
assert close(eperp, 1.0, 0.05)
Sc = 1e4
print(f"plasmoid regime rate ~ Sc^-1/2 = {Sc**-0.5:.3f} for Sc = 1e4; Bhattacharjee Sc = 3e4 -> {3e4**-0.5:.4f}")
for S in (3.14e4, 2.51e5, 6.28e5):   # Bhattacharjee et al. 2009 Table 1: S^-1/2 = 5.6e-3, 2e-3, 1.26e-3
    print(f"  S = {S:.2e}: S^-1/2 = {S**-0.5:.2e}")
assert close((3.14e4) ** -0.5, 5.6e-3, 0.01) and close((2.51e5) ** -0.5, 2.0e-3, 0.01)
assert close((6.28e5) ** -0.5, 1.26e-3, 0.01)

print("\n== Q4/Q5 fusion fuels and machine numbers ==")
m = {"n": 1.00866491595, "H1": 1.00782503223, "H2": 2.01410177812, "H3": 3.01604928132,
     "He3": 3.01602932265, "He4": 4.00260325413, "B11": 11.00930517}   # atomic masses, u (AME2016)
Q = {
    "D+T -> He4+n": m["H2"] + m["H3"] - m["He4"] - m["n"],
    "D+He3 -> He4+p": m["H2"] + m["He3"] - m["He4"] - m["H1"],
    "D+D -> T+p": 2 * m["H2"] - m["H3"] - m["H1"],
    "D+D -> He3+n": 2 * m["H2"] - m["He3"] - m["n"],
    "p+B11 -> 3He4": m["H1"] + m["B11"] - 3 * m["He4"],
}
for k, du in Q.items():
    print(f"{k}: Q = {du*U_MEV:.4f} MeV")
assert close(Q["D+T -> He4+n"] * U_MEV, 17.589, 1e-3)
assert close(Q["D+He3 -> He4+p"] * U_MEV, 18.353, 1e-3)
assert close(Q["p+B11 -> 3He4"] * U_MEV, 8.68, 2e-3)
assert close(Q["D+D -> T+p"] * U_MEV, 4.03, 2e-3) and close(Q["D+D -> He3+n"] * U_MEV, 3.27, 2e-3)
# two-body split, nonrelativistic: light product gets Q*M_heavy/(M_heavy+M_light)
qdhe = Q["D+He3 -> He4+p"] * U_MEV
print(f"D-He3 split: alpha {qdhe*m['H1']/(m['H1']+m['He4']):.2f} MeV, proton {qdhe*m['He4']/(m['H1']+m['He4']):.2f} MeV")

for label, TC in (("Trenta/Norm 100 MdegC", 100e6), ("Polaris 150 MdegC", 150e6),
                  ("Helion goal 200 MdegC", 200e6), ("Da Vinci 3 GdegC", 3e9)):
    print(f"{label}: {(TC+273.15)*KB_EV/1e3:.2f} keV")
assert close(150e6 * KB_EV / 1e3, 12.93, 1e-3)

for v in (250e3,):
    print(f"merging FRCs at {v/1e3:.0f} km/s each: KE per deuteron = {0.5*MD*v*v/E:.0f} eV, "
          f"per proton = {0.5*MP*v*v/E:.0f} eV")


def q_eng(q_sci, eta_elec, eta_E):
    """Wurzel & Hsu 2022 Eq. (27)."""
    return eta_elec * eta_E * (q_sci + 1) - 1


for eta in (0.95, 0.90, 0.80, 0.50):
    print(f"eta_elec = eta_E = {eta}: Q_eng > 0 needs Q_sci > {1/(eta*eta)-1:.3f}")
assert abs(q_eng(1 / 0.95 ** 2 - 1, 0.95, 0.95)) < 1e-12

# Triple products (Wurzel & Hsu Table 7): C-2W 2021 shot 118340
tp = 3.5 * 1.3e19 * 0.005
print(f"C-2W 118340: n T tau = {tp:.3e} keV m^-3 s (published 2.3e17); fraction of 3e21 = {tp/3e21:.1e}")
assert close(tp, 2.3e17, 0.02)
print(f"Trenta separate maxima product 8 keV*3e22*5e-4 s = {8*3e22*5e-4:.1e} (not a simultaneous triple product)")

print("\n== Q6 ball lightning ==")
for d in (0.04, 0.2, 5.0):
    V = math.pi * d ** 3 / 6
    print(f"d = {d} m: V = {V:.3e} m^3; virial cap 3 p0 V = {3*P_ATM*V:.3e} J; "
          f"at 1-10 MJ/m^3: {1e6*V:.2e}-{1e7*V:.2e} J; at 3.9 MJ/m^3: {3.9e6*V:.2e} J")
print(f"virial energy-density cap 3 p0 = {3*P_ATM/1e6:.3f} MJ/m^3; max B held by 1 atm = "
      f"{math.sqrt(2*MU0*P_ATM):.3f} T; B storing 3.9 MJ/m^3 = {math.sqrt(2*MU0*3.9e6):.2f} T")
# Si + O2 -> SiO2, dH_f = -911 kJ/mol; Si 2329 kg/m^3, 28.0855 g/mol
u_si = 911e3 * 2329 / 28.0855e-3
print(f"Si oxidation energy per m^3 of solid Si = {u_si:.2e} J/m^3; Si volume fraction for 3.9 MJ/m^3 = "
      f"{3.9e6/u_si:.1e} ({3.9e6/u_si*2329*1e3:.0f} g Si per m^3)")
print("\nall asserts passed")
