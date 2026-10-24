"""Recompute the key numbers for gravity_unruh_hawking_links.md (SI, CODATA 2018 exact/recommended)."""
import math

hbar = 1.054571817e-34; c = 299792458.0; kB = 1.380649e-23; e = 1.602176634e-19
G = 6.67430e-11; me = 9.1093837015e-31; mp = 1.67262192369e-27; mu = 1.66053906660e-27
eps0 = 8.8541878128e-12; alpha = 7.2973525693e-3; h = 2 * math.pi * hbar
Msun = 1.9884e30; Rsun = 6.957e8; Lsun = 3.828e26; yr = 3.15576e7; g0 = 9.80665
keV = 1e3 * e; MeV = 1e6 * e
mRb87 = 86.909180527 * mu

out = {}
def show(k, v, fmt="{:.4g}"):
    out[k] = v; print(f"{k:58s} {fmt.format(v)}")

# --- Unruh ---
TU = lambda a: hbar * a / (2 * math.pi * c * kB)
show("Unruh coefficient hbar/(2 pi c kB) [K s^2/m]", TU(1.0))
show("a for T_U = 1 K [m/s^2]", 1 / TU(1.0))
show("a for T_U = 1 K [g0]", 1 / TU(1.0) / g0)
assert abs(1 / TU(1.0) - 2.47e20) / 2.47e20 < 0.01   # Wikipedia: 2.47e20 m/s^2 ~ 1 K

# --- Schwinger field, Unruh at Schwinger field ---
ES = me**2 * c**3 / (e * hbar)
show("Schwinger field E_S = m^2c^3/(e hbar) [V/m]", ES)
aS = e * ES / me
show("a = eE_S/m = m c^3/hbar [m/s^2]", aS)
show("T_U at E_S [K]", TU(aS))
show("k T_U at E_S [keV] (= m c^2/2pi)", kB * TU(aS) / keV)
# identity: pi E_S/E == m c^2 / (2 kB T_U(eE/m)) for any E
for E in (1e16, 1e17, ES):
    lhs = math.pi * ES / E; rhs = me * c**2 / (2 * kB * TU(e * E / me))
    assert abs(lhs - rhs) / lhs < 1e-12
show("Schwinger exponent at E=E_S: exp(-pi)", math.exp(-math.pi))
# Schwinger rate, spin-1/2, 3+1D: Gamma = (eE)^2/(4 pi^3 hbar^2 c) sum n^-2 exp(-n pi E_S/E)
def schwinger_rate(E, nmax=200):
    pre = (e * E)**2 / (4 * math.pi**3 * hbar**2 * c)
    return pre * sum(n**-2 * math.exp(-n * math.pi * ES / E) for n in range(1, nmax + 1))
show("Schwinger rate at E_S [pairs m^-3 s^-1]", schwinger_rate(ES))
show("Schwinger rate at 0.1 E_S [pairs m^-3 s^-1]", schwinger_rate(0.1 * ES))
# record laser 1.1e23 W/cm^2 -> peak E
I = 1.1e23 * 1e4
Elas = math.sqrt(2 * I / (c * eps0))
show("Peak E of 1.1e23 W/cm^2 laser [V/m]", Elas)
show("E_S / E_laser", ES / Elas)
show("Schwinger exponent pi E_S/E at record laser", math.pi * ES / Elas)
# field energy at E_S and its gravity over a 1 m sphere
u = 0.5 * eps0 * ES**2
show("EM energy density at E_S [J/m^3]", u)
show("equivalent mass density u/c^2 [kg/m^3]", u / c**2)
for r in (1.0, 1e3):
    M = u / c**2 * 4 / 3 * math.pi * r**3
    show(f"GM/(r c^2) for E_S-filled sphere r={r:g} m", G * M / (r * c**2))

# --- storage rings / linac (Crispino et al. 2008 Sec. IV) ---
show("Sokolov-Ternov P0 = 8/(5 sqrt3)", 8 / (5 * math.sqrt(3)))
def tau_ST(E_GeV, rho):
    gam = E_GeV * 1e9 * e / (me * c**2)
    return 8 * me**2 * c * rho**3 / (5 * math.sqrt(3) * alpha * hbar**2 * gam**5)
show("tau_ST SPEAR 3.7 GeV rho=12.7 m [min] (bend-only)", tau_ST(3.7, 12.7) / 60)
show("tau_ST LEP 45.6 GeV rho=3096 m [h] (bend-only)", tau_ST(45.6, 3096) / 3600)
gL = 45.6e9 * e / (me * c**2); aL = gL**2 * c**2 / 3096
show("LEP 45.6 GeV proper accel gamma^2 c^2/rho [m/s^2]", aL)
show("LEP T_U [K]", TU(aL))
show("LEP circular T_eff (DeltaE>>a): x pi/sqrt3 [K]", TU(aL) * math.pi / math.sqrt(3))
a10 = e * 10e6 / me
show("Linac 10 MV/m: a = eE/m [m/s^2]", a10)
show("Linac 10 MV/m: a [g0]", a10 / g0)
show("Linac 10 MV/m: T_U [K]", TU(a10))
show("Linac  1 MV/m: T_U [K]", TU(a10 / 10))
show("Crispino's a=2e16 g0 -> T_U [K]", TU(2e16 * g0))
# Rogers 1988 Penning trap: v=0.6c, B=15 T
v = 0.6 * c; gam = 1 / math.sqrt(1 - 0.36); w0 = e * 15 / (gam * me); aR = gam**2 * v * w0
show("Rogers: omega0 [1/s]", w0); show("Rogers: a [g0]", aR / g0); show("Rogers: T_U [K]", TU(aR))
# Chen-Tajima 1999 estimate
w0l = 2e15; a0 = 100
pre = 9 / math.pi**2 * hbar * w0l / (me * c**2) * a0
show("Chen-Tajima dI_U/dI_L with ln", pre * math.log(a0 / math.pi))
show("Chen-Tajima dI_U/dI_L with log10", pre * math.log10(a0 / math.pi))
aCT = 2 * a0 * w0l * c
show("Chen-Tajima peak proper accel 2 a0 omega0 c [m/s^2]", aCT)
show("Chen-Tajima T_U [K]", TU(aCT))

# --- Hawking ---
TH = lambda M: hbar * c**3 / (8 * math.pi * G * M * kB)
show("T_H(1 Msun) [K]", TH(Msun))
show("T_H(1 Msun) [nK]", TH(Msun) * 1e9)
tev = lambda M: 5120 * math.pi * G**2 * M**3 / (hbar * c**4)
show("t_ev photon-only 5120piG^2M^3/(hbar c^4), 1 Msun [yr]", tev(Msun) / yr)
assert abs(tev(Msun) / yr / 2.140e67 - 1) < 0.03  # Wikipedia 2.140e67 yr; we get 2.095e67 (2% off, unexplained)
show("P = hbar c^6/(15360 pi G^2 M^2), 1 Msun [W]", hbar * c**6 / (15360 * math.pi * G**2 * Msun**2))
# Page 1976: P = 2e-4 hbar c^6 G^-2 M^-2 (massless species, M >> 1e17 g)
kP = 2.0e-4
tauP = lambda M: G**2 * M**3 / (3 * kP * hbar * c**4)
show("Page lifetime coefficient [s per g^3]", tauP(1e-3))
show("Page lifetime 1 Msun [yr]", tauP(Msun) / yr)
show("ratio photon-only formula / Page", tev(Msun) / tauP(Msun))
show("T_H(5.00e11 kg) [K]", TH(5.00e11))
show("k T_H(5.00e11 kg) [MeV]", kB * TH(5.00e11) / MeV)
show("Page massless-only lifetime 5e11 kg [Gyr]", tauP(5e11) / yr / 1e9)
TCMB = 2.7255
show("M with T_H = T_CMB [kg]", hbar * c**3 / (8 * math.pi * G * kB * TCMB))
show("r_s of 1 Msun [km]", 2 * G * Msun / c**2 / 1e3)

# --- charged black holes ---
k4 = 4 * math.pi * eps0
Mcrit = c**4 / (math.sqrt(k4 * G**3) * ES)
show("M where extremal RN horizon field = E_S [kg]", Mcrit)
show("  ... [Msun]", Mcrit / Msun)
Ehor = lambda M: c**4 / (math.sqrt(k4) * G**1.5 * M)
show("extremal horizon field, 10 Msun [V/m]", Ehor(10 * Msun))
show("extremal horizon field, 1e17 g [V/m]", Ehor(1e14))

# --- de Sitter ---
Mpc = 3.0856775814913673e22
H0 = 67.4e3 / Mpc; OmL = 1 - 0.315
TGH = lambda H: hbar * H / (2 * math.pi * kB)
show("H0 [1/s]", H0)
show("T_GH(H0) [K]", TGH(H0))
show("T_GH(H0 sqrt(OmL)) asymptotic dS [K]", TGH(H0 * math.sqrt(OmL)))
# Cai-Kim dS2: T_dS = T_U + sqrt(T_U^2 + T_GH^2); flat limit 2 T_U; E=0 limit T_GH
def N_dS(m, qE, H):  # natural units hbar=c=kB=1, Eq. (30) of arXiv:1407.4569
    gam = math.sqrt((qE / H)**2 + m**2 - H**2 / 4)
    Smu = 2 * math.pi * gam / H; Slam = 2 * math.pi * qE / H**2
    return (math.exp(-(Smu - Slam)) + math.exp(-2 * Smu)) / (1 - math.exp(-2 * Smu))
m = 1.0
for qE in (0.5, 2.0):
    flat = math.exp(-math.pi * m**2 / qE)
    small_H = N_dS(m, qE, 1e-3)
    show(f"N_dS(m=1,qE={qE},H=1e-3)/Schwinger e^(-pi m^2/qE)", small_H / flat)
show("N_dS(m=1,E=0,H=0.2) / exp(-2 pi m/H)", N_dS(1.0, 0.0, 0.2) / math.exp(-2 * math.pi / 0.2))

# --- analogue BH (Steinhauer 2019): kB T_H = 0.35 nK = 0.12 m c_out^2 ---
T = 0.35e-9
cout = math.sqrt(kB * T / (0.12 * mRb87))
show("Steinhauer: implied c_out [mm/s]", cout * 1e3)
kappa = 2 * math.pi * kB * T / hbar
show("Steinhauer: surface gravity kappa=2pi kB T/hbar [1/s]", kappa)
show("T_H(1 Msun)/T_analogue", TH(Msun) / T)

# --- graphene Klein-Schwinger (Schmitt et al 2023) ---
Dl = 0.2 * e; vF = 1e6
show("graphene E_S = Delta^2/(e hbar vF), Delta=0.2 eV [V/m]", Dl**2 / (e * hbar * vF))
show("G0 = 1.20 x 4 e^2/h [mS]", 1.20 * 4 * e**2 / h * 1e3)
show("c/vF", c / vF)

# --- HRP zero-point-field inertia ---
wP = math.sqrt(c**5 / (hbar * G))
show("Planck angular frequency [1/s]", wP)
Gz = lambda mass, wc: 2 * math.pi * mass * c**2 / (hbar * wc**2)  # m = Gz hbar wc^2/(2 pi c^2)
show("HRP Gamma_z for electron, wc = Planck [s]", Gz(me, wP))
show("HRP Gamma_z for proton, wc = Planck [s]", Gz(mp, wP))
show("Planck time [s]", math.sqrt(hbar * G / c**5))
uZPF = hbar * wP**4 / (8 * math.pi**2 * c**3)
show("ZPF energy density, Planck cutoff [J/m^3]", uZPF)
rhoc = 3 * H0**2 / (8 * math.pi * G)
uL = OmL * rhoc * c**2
show("observed dark-energy density [J/m^3]", uL)
show("log10 ratio ZPF(Planck)/observed", math.log10(uZPF / uL))
# cutoff that would match observed vacuum energy
wmatch = (uL * 8 * math.pi**2 * c**3 / hbar) ** 0.25
show("cutoff matching rho_Lambda: hbar w [meV]", hbar * wmatch / e * 1e3)

# --- Sun: hydrostatic estimates ---
Pc_ssm = 2.477e16; Tc_ssm = 1.571e7; rhoc_ssm = 1.622e5
show("GM^2/R^4 [Pa]", G * Msun**2 / Rsun**4)
show("uniform-density P_c = 3GM^2/(8 pi R^4) [Pa]", 3 * G * Msun**2 / (8 * math.pi * Rsun**4))
# Lane-Emden n=3 by RK4
def lane_emden(n, h=1e-4):
    xi = 1e-6; th = 1 - xi**2 / 6; dth = -xi / 3
    f = lambda x, y, z: (z, -max(y, 0)**n - 2 * z / x)
    while th > 0:
        k1 = f(xi, th, dth); k2 = f(xi + h/2, th + h/2*k1[0], dth + h/2*k1[1])
        k3 = f(xi + h/2, th + h/2*k2[0], dth + h/2*k2[1]); k4 = f(xi + h, th + h*k3[0], dth + h*k3[1])
        th_n = th + h/6*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0]); dth_n = dth + h/6*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1])
        if th_n <= 0:
            fr = th / (th - th_n); return xi + fr * h, dth + fr * (dth_n - dth)
        xi += h; th, dth = th_n, dth_n
xi1, dth1 = lane_emden(3)
show("Lane-Emden n=3 xi1", xi1); show("  -xi1^2 theta'(xi1)", -xi1**2 * dth1)
assert abs(xi1 - 6.89685) < 1e-3
rhobar = Msun / (4 / 3 * math.pi * Rsun**3)
rc_ratio = -xi1 / (3 * dth1)
show("  rho_c/rho_bar", rc_ratio); show("  rho_c n=3 [kg/m^3]", rc_ratio * rhobar)
Pc3 = G * Msun**2 / Rsun**4 / (4 * math.pi * (3 + 1) * dth1**2)
show("  P_c n=3 [Pa]", Pc3); show("  P_c SSM / P_c n=3", Pc_ssm / Pc3)
show("T ~ G M mu m_u/(kB R), mu=0.61 [K]", G * Msun * 0.61 * mu / (kB * Rsun))
show("ideal-gas T_c from SSM P_c, rho_c, mu_c=0.85 [K]", Pc_ssm * 0.85 * mu / (rhoc_ssm * kB))
show("Kelvin-Helmholtz GM^2/(2RL) [Myr]", G * Msun**2 / (2 * Rsun * Lsun) / yr / 1e6)

# --- Coulomb barrier and Gamow peak for p+p at T_c ---
kT = kB * Tc_ssm
show("kT_c [keV]", kT / keV)
EC = e**2 / (k4 * 2e-15)
show("pp Coulomb barrier at r=2 fm [keV]", EC / keV)
show("barrier / kT", EC / kT)
mr = mp / 2
EG = 2 * mr * c**2 * (math.pi * alpha)**2
show("Gamow energy E_G pp [keV]", EG / keV)
E0 = (EG * kT**2 / 4) ** (1 / 3)
show("Gamow peak E0 [keV]", E0 / keV)
show("Gamow width 4 sqrt(E0 kT/3) [keV]", 4 * math.sqrt(E0 * kT / 3) / keV)
show("exp(-3 E0/kT)", math.exp(-3 * E0 / kT))
S0 = 4.01e-22 * keV * 1e-28  # keV b -> J m^2
dS = 11.2 / MeV               # S'/S
def sigv(T):
    kT = kB * T; N = 20000; Emax = 60 * keV; s = 0.0
    for i in range(1, N + 1):
        E = Emax * i / N
        s += S0 * (1 + dS * E) * math.exp(-E / kT - math.sqrt(EG / E)) * Emax / N
    return math.sqrt(8 / (math.pi * mr)) * kT**-1.5 * s
sv = sigv(Tc_ssm)
show("<sigma v>_pp at T_c [m^3/s] (no screening)", sv)
Xc = 0.34
npc = Xc * rhoc_ssm / (mp)
show("n_p centre (X=0.34) [m^-3]", npc)
rpp = npc**2 * sv / 2
show("pp reaction rate density [m^-3 s^-1]", rpp)
show("proton lifetime vs pp 1/(n_p<sv>) [Gyr]", 1 / (npc * sv) / yr / 1e9)
show("power density, 13.1 MeV per pp (ppI) [W/m^3]", rpp * 13.1 * MeV)
show("Lawson nT tau for Sun centre (n_e*kT*tauKH) [keV s/m^3]",
     (rhoc_ssm / mu * (Xc + 0.5 * (1 - Xc))) * kT / keV * G * Msun**2 / (2 * Rsun * Lsun))

# --- extra ties ---
show("T_U at a = g0 (Earth surface) [K]", TU(g0))
show("T_U at a = 274 m/s^2 (Sun surface) [K]", TU(274.0))
Pdrain = schwinger_rate(ES) * 2 * me * c**2
show("pair rest-mass power drain at E_S [W/m^3]", Pdrain)
show("time to drain field energy at E_S (rest mass only) [s]", u / Pdrain)
show("BH mass with T_H = m_e c^2/(2 pi kB) [kg]", hbar * c**3 / (8 * math.pi * G * kB * TU(aS)))
