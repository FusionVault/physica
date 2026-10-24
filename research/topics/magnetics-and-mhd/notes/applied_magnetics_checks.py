"""Recompute the test values in applied_magnetics_engineering.md (standard library only).

Run from outside this folder:  python3 -I /path/to/applied_magnetics_checks.py
Every `check(...)` compares a recomputed number with a published one; a failure raises.
Lines tagged ILLUSTRATION use assumed inputs (no published test case) and are not asserts.
"""
import math

MU0 = 4e-7 * math.pi          # H/m; 2019 SI value 1.25663706212e-6 differs by 5e-10 relative
G0 = 9.80665                  # m/s^2 standard gravity (used in Isp definitions)
MP = 1.67262192e-27           # kg
HBAR = 1.054571817e-34        # J s
PSI = 6894.757                # Pa
LB = 0.45359237               # kg
FT = 0.3048                   # m
KN = 1852 / 3600              # m/s per knot
TNT = 4.184e6                 # J/kg (convention)

FAILS = []


def check(label, got, want, rel):
    ok = abs(got - want) <= rel * abs(want)
    print(f"[{'ok' if ok else 'FAIL'}] {label}: computed {got:.6g}, published {want:.6g} (tol {rel:.1%})")
    if not ok:
        FAILS.append(label)


def illus(label, value, unit=""):
    print(f"[ILLUSTRATION] {label}: {value:.6g} {unit}")


def pmag(B):
    """Magnetic pressure / energy density B^2/(2 mu0), Pa = J/m^3."""
    return B * B / (2 * MU0)


# ---------------------------------------------------------------- 1. record-field arithmetic
print("\n# 1. Record magnets: component sums and Earth-field ratios")
check("Hefei hybrid 2022: 34.22 T resistive + 11 T SC", 34.22 + 11.0, 45.22, 1e-9)
check("NHMFL hybrid (MagLab page): 33.5 + 11.5 T", 33.5 + 11.5, 45.0, 1e-9)
check("NHMFL hybrid (1996 design): 31 + 14 T", 31 + 14, 45.0, 1e-9)
check("LBC3 2019: 14.4 T HTS in 31.1 T resistive", 14.4 + 31.1, 45.5, 1e-9)
check("LBC9 2025: 17.6 T in '31 T' background (rounded; 31.1 T gives 48.7)", 17.6 + 31.1, 48.7, 1e-9)
check("LBC9 stress >150,000 psi in GPa", 150_000 * PSI / 1e9, 1.034, 0.001)
check("LANL 100.75 T ~ 2,000,000 x Earth (50 uT)", 100.75 / 50e-6, 2.0e6, 0.01)
check("Hefei 42.02 T ~ 800,000 x Earth (52.5 uT implied)", 42.02 / 52.5e-6, 8.0e5, 0.01)
check("SECUF 35.6 T ~ 700,000 x Earth (51 uT implied)", 35.6 / 50.9e-6, 7.0e5, 0.01)

print("\n# 2. Magnetic pressure B^2/2mu0 at record fields (MPa)")
for label, B in [("TFMC 20.1 T", 20.1), ("32 T all-SC", 32), ("Hefei 45.22 T", 45.22),
                 ("LBC9 48.7 T", 48.7), ("LANL 100.75 T", 100.75), ("ISSP STC 300 T", 300),
                 ("ISSP EMFC 1200 T", 1200), ("VNIIEF 2800 T", 2800)]:
    print(f"    {label:18s} p = {pmag(B)/1e6:12.1f} MPa")
check("p(45 T) = 806 MPa", pmag(45) / 1e6, 805.7, 0.001)
check("p(100 T) = 3.98 GPa", pmag(100) / 1e9, 3.979, 0.001)
# thin winding: B = mu0 J t  =>  pressure on winding J t B/2 equals B^2/2mu0
J, t = 1.26e9, 5e-3
B = MU0 * J * t
check("thin-winding identity J t B/2 == B^2/2mu0", J * t * B / 2, pmag(B), 1e-12)

# ---------------------------------------------------------------- 3. CFS/MIT TFMC 2021
print("\n# 3. CFS/MIT TFMC (Sept 2021)")
I_tf, B_tf, E_tf, M_tf = 40.5e3, 20.1, 110e6, 9270.0
check("Lorentz line load I*B = 815 kN/m", I_tf * B_tf / 1e3, 815, 0.002)
L_tf = 2 * E_tf / I_tf**2
illus("inferred inductance L = 2E/I^2 (110 MJ at 40.5 kA; not a published value)", L_tf, "H")
illus("TFMC stored energy per coil mass", E_tf / M_tf / 1e3, "kJ/kg")
illus("virial bound sigma/rho at 1 GPa, steel 7900 kg/m3", 1e9 / 7900 / 1e3, "kJ/kg")

# ---------------------------------------------------------------- 4. coil formulas
print("\n# 4. Solenoid, Helmholtz, Fabry factor")


def loop_Bz(R, z, I=1.0, N=1):
    return MU0 * N * I * R**2 / (2 * (R**2 + z**2) ** 1.5)


def solenoid_Bz(z, R, L, n, I):
    """On-axis field of a thin finite solenoid, n turns/m, centred at 0, length L."""
    a, b = z + L / 2, z - L / 2
    return MU0 * n * I / 2 * (a / math.hypot(a, R) - b / math.hypot(b, R))


def helmholtz_B0(R, N, I):
    return (4 / 5) ** 1.5 * MU0 * N * I / R


check("Helmholtz coefficient (4/5)^1.5", (4 / 5) ** 1.5, 0.715542, 1e-5)
R = 0.1
pair = lambda z: loop_Bz(R, z - R / 2) + loop_Bz(R, z + R / 2)
check("Helmholtz centre field = sum of two loops", pair(0), helmholtz_B0(R, 1, 1), 1e-12)
h = 1e-4
d2 = (pair(h) - 2 * pair(0) + pair(-h)) / h**2
check("Helmholtz d2B/dz2 at centre ~ 0 (relative to B/R^2)", abs(d2) * R**2 / pair(0) + 1, 1.0, 1e-5)
check("finite solenoid -> mu0 n I when L >> R", solenoid_Bz(0, 0.01, 100.0, 1000, 1.0), MU0 * 1000, 1e-6)
check("finite solenoid end field -> half of centre (long)", solenoid_Bz(50.0, 0.01, 100.0, 1000, 1.0),
      MU0 * 1000 / 2, 1e-4)


def fabry_G_SI(alpha, beta):
    """B0 = mu0 * G * sqrt(P*lam/(rho*a1)) for uniform J; alpha = a2/a1, beta = half-length/a1."""
    F = beta * math.log((alpha + math.hypot(alpha, beta)) / (1 + math.hypot(1, beta)))
    return F / math.sqrt(2 * math.pi * beta * (alpha**2 - 1))


best = max((fabry_G_SI(1 + i * 0.02, 0.05 + j * 0.02), 1 + i * 0.02, 0.05 + j * 0.02)
           for i in range(1, 400) for j in range(0, 400))
G_cgs = best[0] * 4 * math.pi / 10   # Montgomery practical units: B[G] = G sqrt(P[W] lam/(rho[ohm cm] a1[cm]))
print(f"    uniform-J Fabry factor max: G_SI={best[0]:.4f}, G_Montgomery={G_cgs:.4f} at alpha={best[1]:.2f}, beta={best[2]:.2f}")
check("Fabry max (uniform J) ~0.179 in Montgomery units (literature value, from memory)", G_cgs, 0.179, 0.005)
# what a Bitter coil needs: B0 for P = 26.9 MW, copper rho 2e-8, lam 0.8, a1 = 16 mm, G_SI=0.142
illus("uniform-J 26.9 MW, a1 16 mm, rho 2e-8, lam 0.8 (ideal upper bound)",
      MU0 * best[0] * math.sqrt(26.9e6 * 0.8 / (2e-8 * 0.016)), "T")

# ---------------------------------------------------------------- 5. REBCO scaling (Senatore 2016)
print("\n# 5. REBCO Jc(T,B) = Jc00 exp(-T/T*) B^-alpha  (Senatore et al. SuST 29 014002)")


def rebco_jc_ratio(T1, B1, T2, B2, Tstar=25.0, alpha=0.65):
    return math.exp(-(T2 - T1) / Tstar) * (B2 / B1) ** (-alpha)


illus("Jc(4.2 K)/Jc(20 K) at fixed B, T*=25 K", rebco_jc_ratio(20, 10, 4.2, 10), "x")
illus("Jc(20 T)/Jc(12 T) at fixed T, alpha=0.65", rebco_jc_ratio(20, 12, 20, 20), "x")
check("temperature lift factor independent of start T (exp(-dT/T*))",
      rebco_jc_ratio(10, 5, 20, 5) / rebco_jc_ratio(30, 5, 40, 5), 1.0, 1e-12)

# ---------------------------------------------------------------- 6. SMES virial bound
print("\n# 6. SMES: virial theorem M_struct >= rho E / sigma")
check("Tixador 100 MPa, 8000 kg/m3 -> ~10 kJ/kg (sigma/rho=12.5)", 100e6 / 8000 / 1e3, 12.5, 1e-9)
check("CERN 2025 abstract: 400 MPa / 8000 -> 50 kJ/kg (upper end)", 400e6 / 8000 / 1e3, 50, 1e-9)
check("BPA 1983 SMES 30 MJ = 8.333 kWh", 30e6 / 3.6e6, 8.333, 1e-3)
illus("BPA SMES full-power duration 30 MJ / 10 MW", 30e6 / 10e6, "s")
illus("Bullard-Haugan 30 Wh/kg in kJ/kg", 30 * 3.6, "kJ/kg")
illus("structure stress needed for 30 Wh/kg at rho 8000 (virial)", 30 * 3600 * 8000 / 1e6, "MPa")

# ---------------------------------------------------------------- 7. diamagnetic levitation
print("\n# 7. Berry & Geim 1997 (Eur. J. Phys. 18, 307)")
chi_w, rho_w, g_bg = 8.8e-6, 1000.0, 9.81
BdB = MU0 * rho_w * g_bg / chi_w
check("B dB/dz for water = 1400.9 T^2/m", BdB, 1400.9, 1e-4)
B0_pred = math.sqrt(BdB / (0.63 * 8.15))
check("required centre field with B(zi)=0.63 B0, grad -8.15 B0 T/m -> 16.5 T", B0_pred, 16.5, 0.003)
m_dip = chi_w / MU0 * (4 / 3) * math.pi * 0.01**3 * 16.5
check("induced loop current for 10 mm frog ~1.5 A (using B0)", m_dip / (math.pi * 0.01**2), 1.5, 0.05)
check("induced field |chi| B0 ~ 1.5 gauss", chi_w * 16.5 / 1e-4, 1.5, 0.05)

# ---------------------------------------------------------------- 8. Inductrack / maglev
print("\n# 8. Inductrack (Post, LLNL patent US6629503) and maglev")
check("L/D=30 at 50 km/h, drag 6500 N -> lift ~ 20,000 kg car", 6500 * 30 / G0, 20000, 0.01)
check("drag power 6500 N at 50 km/h ~ 90 kW", 6500 * 50 / 3.6 / 1e3, 90, 0.01)
vt1 = 500 / 300      # L/D = v / v_t  (L/D = kvL/R)
vt2 = 50 / 30
check("transition speed from L/D 300 @500 km/h equals that from L/D 30 @50 km/h", vt1, vt2, 1e-9)
illus("implied Inductrack transition speed v_t = R/(kL)", vt1, "km/h")
illus("SCMaglev 603 km/h record (21 Apr 2015)", 603 / 3.6, "m/s")

# ---------------------------------------------------------------- 9. MHD generator
print("\n# 9. Faraday MHD generator  P = sigma u^2 B^2 k (1-k)")


def mhd_power_density(sigma, u, B, k):
    return sigma * u**2 * B**2 * k * (1 - k)


kmax = max((mhd_power_density(1, 1, 1, i / 1000), i / 1000) for i in range(1001))[1]
check("power density maximal at load factor k = 0.5", kmax, 0.5, 1e-9)
illus("sigma 10 S/m, u 1000 m/s, B 5 T, k 0.5 (assumed seeded-combustion values)",
      mhd_power_density(10, 1000, 5, 0.5) / 1e6, "MW/m^3")
illus("Avco Mk V: 32 MWe / 550 MWth", 32 / 550 * 100, "% enthalpy->electric")
illus("Avco Mk V per Wikipedia: 35 MW gross - 8 MW magnet", 35 - 8, "MW net")
illus("U-25: 20.4 MWe / 250 MWth combustor", 20.4 / 250 * 100, "%")

# ---------------------------------------------------------------- 10. MHD seawater thrust, flowmeter
print("\n# 10. Seawater MHD thrust (Yamato-1) and EM flowmeter")
SIG_SEA = 42.914e-3 * 100    # C(35,15,0) = 42.914 mS/cm -> S/m
check("seawater C(35,15,0) in S/m", SIG_SEA, 4.2914, 1e-6)
illus("Yamato-1: effective electrode length F/(I B) = 8000 N /(4000 A x 4 T)", 8000 / (4000 * 4), "m")


def mhd_thruster_eta_e(sigma, u, B, f):
    """electrical efficiency of a Faraday thruster: f = J B thrust density (N/m^3), u duct velocity."""
    return sigma * u * B**2 / (sigma * u * B**2 + f)


for Bt in (4, 10, 20):
    illus(f"eta_e at sigma 4.29, u 5 m/s, B {Bt} T, f 2e4 N/m^3", mhd_thruster_eta_e(SIG_SEA, 5, Bt, 2e4))
check("Yamato 8 kn in m/s", 8 * KN, 4.116, 1e-3)
check("Yamato 15 km/h in knots", 15 / 3.6 / KN, 8.1, 0.003)
illus("flowmeter U = B D v: 5 mT, 0.1 m, 1 m/s", 5e-3 * 0.1 * 1.0 * 1e3, "mV")

# ---------------------------------------------------------------- 11. railgun, skin depth
print("\n# 11. Railgun F = 1/2 L' I^2, skin depth")
check("NavWeaps saboted 44 lb round at nominal 2500 m/s, 'about 63 MJ'", 0.5 * 44 * LB * 2500**2 / 1e6, 63, 0.015)
check("8,268 fps = 2520 m/s", 8268 * FT, 2520, 0.001)
ke08 = 0.5 * 7 * LB * (8268 * FT) ** 2 / 1e6
check("Jan 2008 shot: 7 lb at 8,268 fps vs published 10.64 MJ (expect ~5% gap)", ke08, 10.64, 0.06)
illus("mass implied by 33 MJ at 2.5 km/s", 2 * 33e6 / 2500**2, "kg")
Lp = 0.5e-6
illus("current for 33 MJ over 10 m of rail, L'=0.5 uH/m: I = sqrt(2F/L')",
      math.sqrt(2 * (33e6 / 10) / Lp) / 1e6, "MA")


def skin_depth(rho, f, mur=1.0):
    return math.sqrt(2 * rho / (2 * math.pi * f * MU0 * mur))


check("Cu 60 Hz skin depth 8.42 mm (sigma 5.96e7 S/m)", skin_depth(1 / 5.96e7, 60) * 1e3, 8.42, 0.002)
check("Cu 50 Hz skin depth 9.22 mm", skin_depth(1 / 5.96e7, 50) * 1e3, 9.22, 0.002)
check("Cu 1 MHz skin depth 65.2 um", skin_depth(1 / 5.96e7, 1e6) * 1e6, 65.2, 0.002)
check("Cu 60 Hz 'about 8.5 mm' with sigma 58.5e6", skin_depth(1 / 58.5e6, 60) * 1e3, 8.5, 0.002)

# ---------------------------------------------------------------- 12. flux compression
print("\n# 12. Flux compression B = lambda B0 (r0/r)^2")


def flux_comp(B0, r0, r, lam=1.0):
    return lam * B0 * (r0 / r) ** 2


illus("ideal: 3 T seed, r0/r = 25", flux_comp(3, 25, 1), "T")
illus("flux retained to turn 3 T, r0/r=25 into 2800 T", 2800 / flux_comp(3, 25, 1))
check("VNIIEF 170 kg TNT ~ 680 MJ (4.0 MJ/kg; 4.184 gives 711)", 170 * TNT / 1e6, 680, 0.05)
illus("MagLIF BR 0.4 MG cm in T m", 0.4e6 * 1e-4 * 1e-2, "T m")
illus("MagLIF: 10 T seed in 2.3 mm liner -> r for BR = 0.4 T m (flux conserved)",
      10 * 2.3e-3**2 / 0.4 * 1e6, "um")

# ---------------------------------------------------------------- 13. electric propulsion
print("\n# 13. Thrusters: eta = T g0 Isp / (2 P)")


def eta(T, Isp, P):
    return T * G0 * Isp / (2 * P)


check("SPT-100B: 83 mN, 1600 s, 1350 W -> 48%", eta(0.083, 1600, 1350), 0.48, 0.01)
check("SPT-100: 80 mN, 1600 s, 1350 W -> 47% (2-figure rounding)", eta(0.080, 1600, 1350), 0.47, 0.015)
illus("AEPS: 580 mN, 2600 s, 12.5 kW", eta(0.58, 2600, 12.5e3))
illus("X3: thrust-to-power at 5.4 N / 99 kW", 5.4 / 99, "N/kW")
isp_lo, isp_hi = (2 * e * 99e3 / (G0 * 5.4) for e in (0.54, 0.67))
illus("X3 Isp implied by eta 0.54-0.67 at 5.4 N/99 kW", isp_lo, f"- {isp_hi:.0f} s (published range 1800-2650 s)")
check("VX-200 2011: 5.7 N, 5000 s, 200 kW -> 72%", eta(5.7, 5000, 200e3), 0.72, 0.03)
illus("VX-200 JPP 2014 72% at 4900 s, 200 kW -> implied thrust", 2 * 0.72 * 200e3 / (G0 * 4900), "N")
check("Zubrin: 200 MW at 20 kg/kW = 4000 t", 200e3 * 20 / 1e3, 4000, 1e-9)


def maecker(I, ra_rc):
    return MU0 * I**2 / (4 * math.pi) * (math.log(ra_rc) + 0.75)


illus("Maecker self-field MPD thrust, 20 kA, ra/rc = 5", maecker(20e3, 5), "N")

# ---------------------------------------------------------------- 14. sails
print("\n# 14. Magnetic sail and E-sail")
n_sw, u_sw = 5e6, 450e3
rho_sw = n_sw * MP
illus("solar-wind dynamic pressure rho u^2 at 1 AU (n 5/cc, 450 km/s)", rho_sw * u_sw**2 * 1e9, "nPa")
for cso in (0.5, 1, 2):
    illus(f"magnetopause field B_mp = u sqrt(rho mu0 C_so), C_so={cso}", u_sw * math.sqrt(rho_sw * MU0 * cso) * 1e9, "nT")
F_ms = 1.0 * 0.5 * rho_sw * u_sw**2 * math.pi * 100e3**2
illus("magsail drag C_d=1, R_mp=100 km", F_ms, "N")
illus("same after Freeland's /3.1 correction", F_ms / 3.1, "N")
illus("E-sail tether for 1 N at 0.5 mN/km", 1 / 0.5e-3, "km")
illus("E-sail tether mass for 1 N at 10 g/km", 2000 * 10e-3, "kg")
illus("E-sail thrust per length at 5 AU (1/r scaling)", 500 / 5, "nN/m")

# ---------------------------------------------------------------- 15. MRI, sensors, magnetocaloric
print("\n# 15. MRI, magnetometers, magnetocaloric")
GAMMA_P = 42.577478461      # MHz/T, CODATA 2022 gamma_p/2pi
GAMMA_P_SH = 42.57638543    # MHz/T, shielded proton in H2O sphere 25 C
for Bm in (1.5, 3.0, 7.0, 11.7):
    print(f"    Larmor at {Bm:4.1f} T: {GAMMA_P*Bm:9.4f} MHz (water-shielded {GAMMA_P_SH*Bm:9.4f})")
check("3 T Larmor 127.73 MHz", GAMMA_P * 3, 127.732, 1e-5)
check("IEC partial-body normal mode 10-8R at R=0.5", 10 - 8 * 0.5, 6.0, 1e-9)
eps = (0.16e-15) ** 2 * 0.45e-6 / (2 * MU0) / HBAR
check("SERF 160 aT/rtHz in 0.45 cm^3 -> energy resolution 44 hbar", eps, 44, 0.02)
illus("Gd magnetocaloric dT at 7 T (3.6 K at 1 T, ~B^0.7)", 3.6 * 7**0.7, "K")

print()
if FAILS:
    raise SystemExit(f"{len(FAILS)} check(s) failed: {FAILS}")
print("all published-value checks passed")
