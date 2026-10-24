"""Recompute key numbers for alternative_gravity_theories.md (run: python3 -I altgrav_checks.py)."""
import math

c = 299_792_458.0
G = 6.67430e-11
hbar = 1.054571817e-34
h = 2 * math.pi * hbar
eV = 1.602176634e-19
hc_eVm = h * c / eV                     # 1.23984e-6 eV m
GM_sun = 1.32712440018e20
R_sun = 6.957e8
AU = 1.495978707e11
Mpc = 3.0856775814913673e22
ARCSEC = math.pi / 180 / 3600
m_e, m_n = 9.1093837015e-31, 1.67492749804e-27

def out(name, val, fmt="{:.4g}"):
    print(f"{name:62s} {fmt.format(val)}")

print("== PPN / solar system ==")
# Brans-Dicke: gamma = (1+w)/(2+w)  =>  gamma-1 = -1/(2+w)
bd_gamma = lambda w: (1 + w) / (2 + w)
lo = (2.1 - 2 * 2.3) * 1e-5            # Cassini central - 2 sigma
w_min = -1 / lo - 2
out("Brans-Dicke omega bound from Cassini (2 sigma lower edge)", w_min)
assert abs(bd_gamma(40000) - 1 + 1 / 40002) < 1e-15
out("Brans-Dicke gamma-1 at omega=40000", bd_gamma(40000) - 1)
out("Brans-Dicke Nordtvedt eta = 4b-g-3 = 1/(2+w) at w=40000", 4 * 1 - bd_gamma(40000) - 3)

defl = 4 * GM_sun / (c**2 * R_sun) / ARCSEC
out("Light deflection at solar limb, GR (arcsec) [(1+g)/2 x this]", defl)

# Mercury perihelion (GR): 6 pi GM / (c^2 a (1-e^2)) per orbit
a, e, P_days = 5.7909050e10, 0.205630, 87.9691
dphi = 6 * math.pi * GM_sun / (c**2 * a * (1 - e**2))
per_century = dphi * (36525 / P_days) / ARCSEC
out("Mercury GR perihelion advance (arcsec/century)", per_century)
ppn_fac = lambda g, b: (2 + 2 * g - b) / 3
out("  x (2+2g-b)/3 for Brans-Dicke w=40000", per_century * ppn_fac(bd_gamma(40000), 1))
out("  x (2+2g-b)/3 for vDVZ-like gamma=1/2, beta=1 (no Vainshtein)", per_century * ppn_fac(0.5, 1))

# Shapiro round-trip excess delay, Cassini 2002 conjunction geometry (r_E=1 AU, r_C=8.43 AU, b=1.6 R_sun)
r1, r2, b = AU, 8.43 * AU, 1.6 * R_sun
shapiro = 2 * (1 + 1) * GM_sun / c**3 * math.log(4 * r1 * r2 / b**2)
out("Cassini-geometry round-trip Shapiro delay, GR (microseconds)", shapiro * 1e6)

print("\n== Graviton mass <-> Compton wavelength lambda_g = h/(m_g c) ==")
for label, m in [("GWTC-3 v1 (2021) 1.27e-23 eV", 1.27e-23), ("GWTC-3 v3 / PRD 2025 2.42e-23 eV", 2.42e-23),
                 ("GWTC-4.0 II (2026) 1.92e-23 eV", 1.92e-23), ("INPOP19a 90% 3.16e-23 eV", 3.16e-23),
                 ("INPOP17b 90% 6.76e-23 eV", 6.76e-23)]:
    out(f"lambda_g (km) for {label}", hc_eVm / m / 1e3)
assert abs(hc_eVm / 3.16e-23 / 1e3 - 3.93e13) / 3.93e13 < 0.01   # matches Bernus et al. 2020 abstract
for lam_km in (1.2e14, 2.2e14):
    out(f"m_g (eV) for Will 2018 Mars bound lambda_g={lam_km:.1e} km", hc_eVm / (lam_km * 1e3))
lam = hc_eVm / 1.92e-23
out("1 - v_g/c at f=100 Hz for m_g=1.92e-23 eV  [0.5 (c/(lam f))^2]", 0.5 * (c / (lam * 100)) ** 2)
out("Yukawa factor exp(-r/lambda) at r=1 AU for m_g=1.92e-23 eV", math.exp(-AU / lam))

print("\n== GW170817 speed bound ==")
D = 26 * Mpc                      # conservative lower 90% distance used in ApJL 848 L13
out("upper: +1.74 s * c / D", 1.74 * c / D)
out("lower: (1.74 - 10) s * c / D  (GW emitted up to 10 s before gamma)", (1.74 - 10) * c / D)

print("\n== MOND / emergent gravity ==")
H0 = 67.4e3 / Mpc
a0 = 1.2e-10
out("c H0 (m/s^2), H0=67.4", c * H0)
out("c H0 / (2 pi)", c * H0 / (2 * math.pi))
out("c H0 / 6  (Verlinde point-mass MOND-like scale)", c * H0 / 6)
out("MOND radius sqrt(GM_sun/a0) (AU)", math.sqrt(GM_sun / a0) / AU)
M = 1e11 * 1.98847e30
out("BTFR V_flat = (G M a0)^(1/4) for M_b=1e11 M_sun (km/s)", (G * M * a0) ** 0.25 / 1e3)
rar = lambda gbar, gd=1.20e-10: gbar / (1 - math.exp(-math.sqrt(gbar / gd)))
for gbar in (1e-12, 1.2e-10, 1e-8):
    out(f"RAR g_obs for g_bar={gbar:.1e} (m/s^2)", rar(gbar))
assert abs(rar(1e-14) / math.sqrt(1e-14 * 1.2e-10) - 1) < 0.01   # deep-MOND limit sqrt(g_bar a0)

print("\n== Mach / Sciama ==")
rho_c = 3 * H0**2 / (8 * math.pi * G)
out("rho_crit (kg/m^3), H0=67.4", rho_c)
out("G rho_c / H0^2 = 3/(8 pi)", G * rho_c / H0**2)
out("4 pi G rho_c T^2 (T=1/H0) = 3/2", 4 * math.pi * G * rho_c / H0**2)
R_H = c / H0
M_H = 4 / 3 * math.pi * rho_c * R_H**3
out("G M_H / (R_H c^2) for Hubble sphere at rho_c = 1/2", G * M_H / (R_H * c**2))

# Rodal: GR Machian term suppressed by (U_N/c^2)(omega L/c)^2; example lab device 0.1 kg, L=1 cm, 50 kHz
Mdev, L, f = 0.1, 0.01, 50e3
wl = (2 * math.pi * f * L / c) ** 2
out("(omega L/c)^2 for L=1 cm, f=50 kHz", wl)
out("(U_N/c^2)(omega L/c)^2 with U_N = G M/L, M=0.1 kg", G * Mdev / L / c**2 * wl)

# Assis/Weber: inertial mass shift of a charge q inside a shell at potential U: dm = qU/(3c^2)
out("Weber/Assis electron dm/m inside shell at U=10 kV  [eU/(3 m_e c^2)]", eV * 1e4 / (3 * m_e * c**2))

print("\n== Einstein-Cartan: density where spin-torsion ~ mass energy, rho = m^2 c^4/(8 pi G hbar^2) ==")
for lab, m in (("electron", m_e), ("neutron", m_n)):
    out(f"rho_EC {lab} (g/cm^3)", m**2 * c**4 / (8 * math.pi * G * hbar**2) * 1e-3)

print("\n== Weber-like gravity: numerically integrated perihelion advance ==")
def weber_precession(mu, gam, GM=1.0, a=1.0, e=0.2, cc=100.0, steps=200000):
    """F = -(GM m/r^2)(1 - mu rdot^2 + gam r rddot) r_hat, mu,gam in units of 1/cc^2.
    Returns perihelion advance per orbit (rad) from RK4 in (r, rdot, phi)."""
    mu, gam = mu / cc**2, gam / cc**2
    rp = a * (1 - e)
    hh = math.sqrt(GM * a * (1 - e**2))          # Newtonian angular momentum at perihelion
    def deriv(s):
        r, v, ph = s
        phidot = hh / r**2
        rdd = (r * phidot**2 - GM / r**2 * (1 - mu * v**2)) / (1 + gam * GM / r)
        return (v, rdd, phidot)
    s = (rp, 0.0, 0.0)
    T = 2 * math.pi * math.sqrt(a**3 / GM)
    dt = 1.3 * T / steps
    prev_v = 0.0
    passed_apo = False
    for _ in range(steps):
        k1 = deriv(s)
        k2 = deriv(tuple(x + dt / 2 * k for x, k in zip(s, k1)))
        k3 = deriv(tuple(x + dt / 2 * k for x, k in zip(s, k2)))
        k4 = deriv(tuple(x + dt * k for x, k in zip(s, k3)))
        new = tuple(x + dt / 6 * (p + 2 * q + 2 * r_ + w) for x, p, q, r_, w in zip(s, k1, k2, k3, k4))
        if new[1] < 0: passed_apo = True
        if passed_apo and prev_v < 0 <= new[1]:      # rdot crosses 0 upward -> perihelion
            frac = -prev_v / (new[1] - prev_v)        # linear interpolation in phi
            ph = s[2] + frac * (new[2] - s[2])
            return ph - 2 * math.pi
        prev_v = new[1]
        s = new
    raise RuntimeError("no perihelion found")

gr = 6 * math.pi * 1.0 / (100.0**2 * 1.0 * (1 - 0.2**2))
for lab, mu, gam in (("Tisserand mu=1,gam=2", 1, 2), ("Weber-original mu=1/2,gam=1", 0.5, 1),
                     ("Assis mu=3,gam=6", 3, 6), ("Raguza mu=9,gam=6", 9, 6), ("mu=0,gam=6", 0, 6)):
    out(f"Weber-like {lab}: advance / GR advance", weber_precession(mu, gam) / gr)
assert abs(weber_precession(3, 6) / gr - 1) < 0.03

print("\n== Vacuum energy, Casimir, quantum inequality ==")
cas_E = lambda d: math.pi**2 * hbar * c / (720 * d**3)     # J/m^2, ideal plates
cas_P = lambda d: math.pi**2 * hbar * c / (240 * d**4)     # Pa
for d in (1e-6, 1e-7, 1e-8):
    out(f"Casimir |E|/A at d={d:.0e} m (J/m^2)", cas_E(d))
    out(f"Casimir pressure at d={d:.0e} m (Pa)", cas_P(d))
out("Forward one-shot: 1 m^2 pair collapsed 1 um -> 10 nm (J)", cas_E(1e-8) - cas_E(1e-6))
# Moddel 2021: 1.4 pW max power, 0.02 um^2, 33 nm cavity
Pm, Am = 1.4e-12, 0.02e-12
out("Moddel 2021 claimed power density (W/m^2)", Pm / Am)
E4h = Pm * 4 * 3600
out("Moddel: energy at 1.4 pW for 4 h (J)", E4h)
Ecav = cas_E(33e-9) * Am
out("Ideal Casimir energy of a 33 nm, 0.02 um^2 cavity (J)", Ecav)
out("  ratio 4-h output / one-shot cavity Casimir energy", E4h / Ecav)
qi = lambda tau: 3 * hbar / (32 * math.pi**2 * c**3 * tau**4)   # Ford-Roman 4D massless scalar, J/m^3
for tau in (1e-9, 1e-6, 1.0):
    out(f"Ford-Roman QI: max sampled negative energy density, tau={tau:.0e} s (J/m^3)", qi(tau))
rho_L = 0.6847 * 3 * (67.36e3 / Mpc) ** 2 / (8 * math.pi * G) * c**2
out("Dark-energy density rho_Lambda c^2 (J/m^3), Planck 2018", rho_L)
out("  in 1 km^3 (J)", rho_L * 1e9)
out("  in Earth's volume (J)", rho_L * 4 / 3 * math.pi * 6.371e6**3)
g0 = 9.80665
out("Lift 1 kg by 1 m (J)", g0)
out("Fall 1 kg from infinity to Earth surface GM/R (MJ)", 3.986004418e14 / 6.371e6 / 1e6)
