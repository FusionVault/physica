"""Recompute key numbers for z_pinch_and_plasma_focus.md (stdlib only).

Run: python3 -I z_pinch_checks.py   -- every assert is a self-check against a published value.
"""
import math

MU0 = 4e-7 * math.pi          # H/m (pre-2019 exact value; new SI differs by 1e-10)
E = 1.602176634e-19           # C, also J/eV
KB = 1.380649e-23             # J/K
G = 6.67430e-11               # m^3 kg^-1 s^-2
MP = 1.67262192e-27           # kg
MD = 3.3435837768e-27         # kg (deuteron)
ME = 9.1093837e-31            # kg
MH = 1.6735575e-27            # kg (H atom)
HBAR = 1.054571817e-34
MSUN = 1.98847e30
PC = 3.0856775814913673e16
LY = 9.4607304725808e15


def close(a, b, rel):
    return abs(a - b) <= rel * abs(b)


# ---------- 1. Bennett relation  mu0 I^2/(8 pi) = N k (Te+Ti) ----------
def bennett_line_density(I, Tsum_eV):
    """N (ions per m) for current I (A) and Te+Ti (eV), z=1."""
    return MU0 * I**2 / (8 * math.pi * Tsum_eV * E)


def bennett_T(I, N, z=1):
    """Bennett temperature T_B (eV) with T_e=T_i, line density N (m^-1), charge z."""
    return MU0 * I**2 / (8 * math.pi * N * (1 + z)) / E


# Klir 2005 practical form: T_B[keV] = 3.12 I[MA]^2 / ((1+z) N[1e18 cm^-1])
tb = bennett_T(1e6, 1e20) / 1e3          # 1 MA, 1e18 cm^-1 = 1e20 m^-1
assert close(tb * 2, 3.12, 0.005), tb
print(f"Bennett coefficient (keV, MA, 1e18/cm): {tb*2:.3f}  [published 3.12]")
print(f"N for 1 MA, Te+Ti=2 keV: {bennett_line_density(1e6, 2000):.3e} m^-1")

# FuZE 2019 (Zhang et al. PRL 122, 135001): I~200 kA, n~1e17 cm^-3, a~0.3 cm, T 1-2 keV published
for label, I, n, a in [("FuZE 2019", 200e3, 1e23, 3e-3)]:
    N = n * math.pi * a**2               # Bennett profile: N = pi a^2 n0
    T = bennett_T(I, N)
    print(f"{label}: N={N:.2e}/m  Bennett T_e=T_i={T/1e3:.2f} keV  [published 1-2 keV]")
    B = MU0 * I / (2 * math.pi * a)
    rho = n * MP                         # 80% H / 20% D fill: proton mass is the right order
    vA = B / math.sqrt(MU0 * rho)
    tau = a / vA
    print(f"  B(a)={B:.1f} T  vA={vA/1e3:.0f} km/s  a/vA={tau*1e9:.1f} ns  "
          f"16 us quiescent = {16e-6/tau:.0f} growth times")
    print(f"  shear threshold 0.1 k vA at ka=1: dVz/dr={0.1*vA/a:.2e} 1/s (dV across a = {0.1*vA/1e3:.0f} km/s)")
    fuze_tau = tau

# ---------- 2. Pease-Braginskii current ----------
# Klir 2005 eq 2.20: I_PB = (pi/mu0) sqrt(48 lnL /(sigma0 A)) (1+z)/z ; parabolic n profile,
# P_brems = A (kTe)^1/2 * (z^3 n_i^2) [i.e. z^2 n_e n_i], sigma = sigma0 (kTe)^3/2 /(z lnL).
def i_pb(lnL=10.0, z=1, perpendicular=True):
    c_sig = 9.71e3 if perpendicular else 1.96 * 9.71e3   # NRL: eta_perp = 1.03e-4 z lnL / T_eV^1.5 ohm m
    sigma0 = c_sig / E**1.5                              # S m^-1 J^-3/2
    A = 1.69e-38 / math.sqrt(E)                          # NRL bremsstrahlung, W m^3 J^-1/2
    return (math.pi / MU0) * math.sqrt(48 * lnL / (sigma0 * A)) * (1 + z) / z


ipb = i_pb()
assert close(ipb, 1.4e6, 0.05), ipb
print(f"I_PB (lnL=10, z=1, sigma_perp): {ipb/1e6:.2f} MA [published ~1.4 MA]; "
      f"sigma_par: {i_pb(perpendicular=False)/1e6:.2f} MA; lnL=20: {i_pb(20)/1e6:.2f} MA")

# ---------- 3. Sharp-boundary Z-pinch m=0 / m=1 growth (Kruskal-Schwarzschild geometry, derived) ----------
def bessel_i(m, x, n=2000):
    h = math.pi / n
    s = sum((1 if k in (0, n) else 4 if k % 2 else 2) * math.exp(x * math.cos(k * h)) * math.cos(m * k * h)
            for k in range(n + 1))
    return s * h / 3 / math.pi


def bessel_k(m, x, tmax=12.0, n=4000):
    h = tmax / n
    s = sum((1 if k in (0, n) else 4 if k % 2 else 2) * math.exp(-x * math.cosh(k * h)) * math.cosh(m * k * h)
            for k in range(n + 1))
    return s * h / 3


assert close(bessel_i(0, 1), 1.266066, 1e-5) and close(bessel_i(1, 1), 0.565159, 1e-5)
assert close(bessel_k(0, 1), 0.421024, 1e-5) and close(bessel_k(1, 1), 0.601907, 1e-5)


def growth(m, ka):
    """(gamma a / vA)^2 = ka I_m'/I_m * [1 + m^2 K_m/(ka K_m')]; vA from edge field."""
    ip = (bessel_i(abs(m - 1), ka) + bessel_i(m + 1, ka)) / 2
    kp = -(bessel_k(abs(m - 1), ka) + bessel_k(m + 1, ka)) / 2
    return ka * ip / bessel_i(m, ka) * (1 + m * m * bessel_k(m, ka) / (ka * kp))


for m in (0, 1, 2):
    row = "  ".join(f"ka={x}:{math.sqrt(g) if (g := growth(m, x)) > 0 else 0:.3f}" for x in (0.5, 1, 2, 5))
    print(f"gamma a/vA, m={m}: {row}")
assert close(math.sqrt(growth(0, 1.0)), math.sqrt(0.565159 / 1.266066), 1e-4)

# ---------- 4. Zap pressures ----------
for label, n, T in [("FuZE-3 2025 (n=5.18e24, Te=1 keV)", 5.18e24, 1000),
                    ("Goyon 2024 (4.9e17 cm^-3, 2.25 keV)", 4.9e23, 2250)]:
    print(f"{label}: p_e = {n*T*E/1e6:.0f} MPa")
assert close(5.18e24 * 1000 * E, 830e6, 0.01)

# ---------- 5. Lee model: dynamic resistance and Table I of Lee 2009 APL 95 151503 ----------
c = 2.0
dLdt = 2e-7 * math.log(c) * 1e5          # dz/dt = 10 cm/us
DR0 = dLdt / 2
assert close(dLdt, 1.4e-2, 0.02) and close(DR0, 7e-3, 0.02)
V0, L0 = 30e3, 30e-9
table = [(0.45, 1e-6, 152), (4.5, 10e-6, 447), (45, 100e-6, 1156), (135, 300e-6, 1676),
         (450, 1000e-6, 2321), (1080, 2400e-6, 2781), (4500, 10000e-6, 3407), (45000, 100000e-6, 4209)]
for E0, C0, ipk_pub in table:
    Z0 = math.sqrt(L0 / C0)
    ipk = V0 / (Z0 + 7e-3 + 0.1 * Z0) / 1e3
    assert close(0.5 * C0 * V0**2 / 1e3, E0, 0.02)
    if E0 < 45000:  # paper rounds Ztotal; its 45 MJ row (4209 kA) is inconsistent with its own Ztotal=7.6 mOhm
        assert close(ipk, ipk_pub, 0.015), (E0, ipk)
    print(f"  E0={E0} kJ Z0={Z0*1e3:.2f} mOhm Ipeak={ipk:.0f} kA [published {ipk_pub}]")
print(f"Lee Table I reproduced; asymptote V0/DR0 = {V0/7e-3/1e6:.2f} MA")


def yn_peak(I_MA):
    return 1.8e10 * I_MA**3.8


def yn_pinch(I_MA):
    return 3.2e11 * I_MA**4.5


print(f"Lee fit Yn(Ipeak=1 MA)={yn_peak(1):.2e}; Yn(Ipeak=1.1 MA)={yn_peak(1.1):.2e} "
      f"vs FF-1 2016 measured 2.5e11 -> ratio {2.5e11/yn_peak(1.1):.1f}")


def snowplow_speed(I, a, b, rho0, fm=1.0, fc=1.0):
    """Steady axial run-down speed from rho0 fm pi(b^2-a^2) v^2 = mu0 (fc I)^2 ln(b/a)/(4 pi)."""
    cc = b / a
    return fc / math.sqrt(fm) * math.sqrt(MU0 * math.log(cc) / (4 * math.pi**2 * (cc**2 - 1) * rho0)) * I / a


rho_ff1 = 1.6e24 * MD                    # FF-1 fill 1.6e18 atoms/cm^3 (Lerner 2023)
for fm, fc in ((1, 1), (0.1, 0.7)):
    v = snowplow_speed(1e6, 0.028, 0.05, rho_ff1, fm, fc)
    print(f"FF-1 snowplow speed fm={fm}, fc={fc}: {v/1e4:.1f} cm/us  [Lee: ~10 cm/us typical]")

# ---------- 6. LPPFusion / Lerner numbers ----------
ecyc = HBAR * E * 1e5 / ME / E           # hbar omega_ce in eV per GG (1 GG = 1e5 T)
assert close(ecyc, 11.6, 0.005)
for ion, mass_ratio in (("p", MP / ME), ("alpha", 6.6446573e-27 / ME), ("11B", 11.009305 * 1.66053907e-27 / ME)):
    print(f"QMF threshold B for Ei=300 keV, {ion}: {300e3/(mass_ratio*ecyc):.2f} GG")
Bc = 4 * 3 * 6 * (MP / ME) * 1e5         # Bc = 4 z mu (M/m) B, pB11 z=3 mu=6, B=100 kG
print(f"Lerner Eq.2 pB11, B=100 kG: Bc={Bc/1e9:.1f} GG [published ~10 GG]")
for fuel, mu, z in (("D", 2, 1), ("pB11", 6, 3)):
    print(f"Lerner Eq.9 ion-heating current {fuel}: I > {5.0/(mu**(1/6)*z**(5/6)):.2f} MA")
ntt = 3.5e19 * 40e-9 * 240                # cm^-3 * s * keV
print(f"n tau T = {ntt:.2e} keV s cm^-3 [published 3.4e14]")
e_dd = (3.269 + 4.033) * 1e6 * E          # both DD branches per neutron
print(f"2.5e11 DD neutrons -> {2.5e11*e_dd:.2f} J (both branches) / {2.5e11*3.269e6*E:.2f} J (n branch) "
      f"[published 0.2 J]; eta = 0.2 J/60 kJ = {0.2/60e3:.2e}")
print(f"2013 ion beam 120 kA x 3 MV x 5 ns = {120e3*3e6*5e-9/1e3:.1f} kJ [published 2 kJ]")

# ---------- 7. MagLIF ----------
n_dt = 2e3 / (17.59e6 * E)
print(f"2 kJ DT-equivalent = {n_dt:.2e} DT reactions = {n_dt/1.1e13:.0f} x the 1.1e13 DD neutrons")
print(f"2026 record 2.33e13 scaled at the same ratio: ~{2.33e13/1.1e13*2:.1f} kJ DT-eq (rough; ignores Ti)")

# ---------- 8. Pinch-gravity analogy ----------
cs2 = KB * 10 / (2.3 * MH)
mline = 2 * cs2 / G
print(f"Ostriker critical line mass T=10 K mu=2.3: {mline*PC/MSUN:.1f} Msun/pc [published 16.6]")
assert close(mline * PC / MSUN, 16.6, 0.01)


def ostriker_mass(rhoc, cs2, rmax_factor=1e4, n=200000):
    """Integrate rho_c/(1+r^2/(8H^2))^2 2 pi r dr, H^2 = cs^2/(4 pi G rho_c)."""
    H2 = cs2 / (4 * math.pi * G * rhoc)
    rmax = rmax_factor * math.sqrt(H2)
    h = rmax / n
    return sum(rhoc / (1 + ((k + 0.5) * h)**2 / (8 * H2))**2 * 2 * math.pi * (k + 0.5) * h * h for k in range(n))


for rhoc in (1e-17, 1e-15):
    assert close(ostriker_mass(rhoc, cs2), mline, 1e-3)
print("Ostriker profile integrates to 2 cs^2/G for any rho_c (checked 2 values)")

# Mapping: G m^2 <-> mu0 e^2 u^2/(4 pi), cs^2 <-> k(Te+Ti)/m  =>  Bennett.
u, Tsum = 3e5, 2000.0
N_map = 2 * (Tsum * E / MD) / (MU0 * E**2 * u**2 / (4 * math.pi * MD**2)) / MD
I_map = E * u * N_map
assert close(MU0 * I_map**2 / (8 * math.pi), N_map * Tsum * E, 1e-12)
print(f"Analogy check: u={u/1e3:.0f} km/s, Te+Ti=2 keV -> N={N_map:.3e}/m, I={I_map/1e3:.0f} kA satisfies Bennett")
u_eq = MP * math.sqrt(4 * math.pi * G / MU0) / E
print(f"Drift at which magnetic attraction equals gravity for protons: {u_eq:.2e} m/s")
N1 = bennett_line_density(1e6, 2000)
ratio = MU0 * 1e6**2 / (4 * math.pi) / (G * (N1 * MD)**2)
print(f"1 MA D pinch, Te+Ti=2 keV: line mass {N1*MD:.2e} kg/m; magnetic/gravitational binding = {ratio:.1e}; "
      f"line mass gravity would need = {2*Tsum*E/MD/G:.2e} kg/m")

# ---------- 9. MARAUDER / plasmoid kinematics ----------
for m_mg, v_kms in ((1, 1000), (2, 3000), (1, 10000)):
    print(f"MARAUDER KE {m_mg} mg at {v_kms} km/s = {0.5*m_mg*1e-6*(v_kms*1e3)**2/1e6:.1f} MJ")
a_g = 1e11 * 9.80665
print(f"1e11 g to 1000 km/s: distance {1e12/(2*a_g):.2f} m, time {1e6/a_g*1e6:.2f} us")
print(f"Bostick plasmoids 1e7-2e7 cm/s = {1e5/1e3:.0f}-{2e5/1e3:.0f} km/s")

# ---------- 10. Lockheed T4B / TX energy and beta checks ----------
for label, n, Ti, Te, B, V, W_pub in (("T4B", 5e19, 75, 250, 0.1, 0.2, 1170), ("TX", 5e20, 9.6e3, 12.6e3, 2.3, 16.3, 51e6)):
    p = n * (Ti + Te) * E
    print(f"{label}: p={p:.3e} Pa beta={p/(B**2/(2*MU0)):.2f} thermal W={1.5*p*V:.3g} J [published {W_pub:.3g} J]")

# ---------- 11. Peratt galaxy filaments ----------
I, d = 1e18, 3e5 * LY
print(f"Peratt filaments 1e18 A at 300,000 ly: B={MU0*I/(2*math.pi*d)*1e10:.1f} microgauss, "
      f"F/L={MU0*I**2/(2*math.pi*d):.2e} N/m")
print("all checks passed")
