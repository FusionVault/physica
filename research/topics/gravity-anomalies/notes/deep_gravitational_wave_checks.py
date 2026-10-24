"""Recompute test values for deep GW physics notes. Pure stdlib. Run: python3 -I gwdeep.py"""
import math, cmath

G = 6.67430e-11; c = 299792458.0; hbar = 1.054571817e-34; hP = 6.62607015e-34
kB = 1.380649e-23; eV = 1.602176634e-19; mu0 = 1.25663706212e-6
GMsun = 1.32712440018e20; Msun = GMsun / G; Tsun = GMsun / c**3
Mpc = 3.0856775814913673e22; kpc = Mpc / 1e3; pc = Mpc / 1e6
yr = 365.25 * 86400; erg = 1e-7
me = 9.1093837015e-31; mp = 1.67262192369e-27; alpha = 7.2973525693e-3
out = {}
def p(k, v): out[k] = v; print(f"{k:60s} {v:.6g}")

# 1 Planck luminosity and event energetics
LP = c**5 / G; p("L_P = c^5/G [W]", LP); p("c^5/(4G) [W]", LP / 4)
p("Msun c^2 [J]", Msun * c**2)
p("GW150914 E_rad 3.0 Msun c^2 [J]", 3.0 * Msun * c**2)
L150914 = 3.5e56 * erg
p("GW150914 L_peak 3.5e56 erg/s [W]", L150914)
p("GW150914 200 Msun c^2/s [W]", 200 * Msun * c**2)
p("GW150914 L_peak / L_P", L150914 / LP)
p("GW190521 L_peak 3.7e56 erg/s / L_P", 3.7e56 * erg / LP)
p("GW170817 E_rad > 0.025 Msun c^2 [J]", 0.025 * Msun * c**2)
D150914 = 410 * Mpc
p("GW150914 peak flux L/(4 pi D^2) [W/m^2]", L150914 / (4 * math.pi * D150914**2))
p("GW150914 fluence E/(4 pi D^2) [J/m^2]", 3.0 * Msun * c**2 / (4 * math.pi * D150914**2))
hdot = 2 * math.pi * 250 * 1e-21
p("PRL estimate c^3 hdot^2/(16 pi G), h=1e-21,f=250 Hz [W/m^2]", c**3 * hdot**2 / (16 * math.pi * G))
# Flanagan-Hughes Eq 5.41: T^tz = c^3 h0^2 w^2 /(32 pi G)
w = 2 * math.pi * 1e3
p("Flanagan-Hughes flux h0=1e-22, f=1kHz [mW/m^2]", 1e3 * c**3 * (1e-22)**2 * w**2 / (32 * math.pi * G))

# 2 graviton counts
nu = 2 * math.pi * 150; h0 = 1e-21
p("Tobar N = h0^2 c^5/(32 pi G hbar nu^2), h0=1e-21, 150 Hz", h0**2 * c**5 / (32 * math.pi * G * hbar * nu**2))
def n_dB(h, f):  # Carney Eq 1 in SI: rho/(hbar w) * (c/f)^3, rho = c^2 h^2 w^2/(32 pi G)
    w = 2 * math.pi * f; rho = c**2 * h**2 * w**2 / (32 * math.pi * G)
    return rho / (hbar * w) * (c / f)**3
p("Carney n*lambda_dB^3, h=1e-22, f=1 kHz", n_dB(1e-22, 1e3))
p("Carney n*lambda_dB^3, h=7e-28, f=1e18 Hz", n_dB(7e-28, 1e18))
# frequency at which a strain h has one quantum per de Broglie cell
p("f where n*lambda^3=1 for h=1e-30 [Hz]", 1e3 * math.sqrt(n_dB(1e-30, 1e3)))

# 3 Tobar et al. 2024 rates
def gamma_stim(vs, M, h, l=1): return vs**2 * M * h**2 / (4 * l**2 * math.pi**3 * hbar)
p("Tobar Gamma_stim Al 1800 kg, vs=5.4 km/s, h=5e-22 [Hz]", gamma_stim(5400, 1800, 5e-22))
p("Tobar Gamma_stim Al 1800 kg, vs=5.1 km/s, h=5e-22 [Hz]", gamma_stim(5100, 1800, 5e-22))
p("Tobar Gamma_spon Nb rho=8570 vs=5e3 L=1 R=0.5 [Hz]", 8 * math.pi * G * 8570 * 5e3**4 * 0.5**2 / (1.0 * c**5))
def tobar_mass(vs, h0, Mc, f):
    w = 2 * math.pi * f
    return (24 * math.pi**2 / 5) * hbar / (h0**2 * vs**2) * (G * Mc * Msun / (2 * c**3))**(5 / 3) * w**(8 / 3)
p("Tobar M_opt Be vs=13 km/s, h0=2e-22, Mc=1.19, 100 Hz [kg]", tobar_mass(13e3, 2e-22, 1.19, 100))
p("Tobar M_opt Al vs=5.4 km/s, h0=2e-22, Mc=1.19, 150 Hz [kg]", tobar_mass(5.4e3, 2e-22, 1.19, 150))

# 4 Rothman-Boughn / Dyson arithmetic (Gaussian units where noted)
lP2 = G * hbar / c**3
p("Planck length^2 [cm^2]", lP2 * 1e4)
p("(G m_e/c^2)^2 gravito-Compton [cm^2]", (G * me / c**2)**2 * 1e4)
alpha_g = G * mp**2 / (hbar * c)
p("Jupiter-mass bound (alpha/alpha_g)^1.5 m_p [kg]", (alpha / alpha_g)**1.5 * mp)
p("max ionization sigma = 3072 pi/(5*32) lP^2 [cm^2] (ak=1)", 3072 * math.pi / (5 * 32) * lP2 * 1e4)
Gcgs = G * 1e3; ccgs = c * 100
Lmix = 2 * ccgs**2 / (math.sqrt(Gcgs) * 1.0)  # cm, B = 1 G
p("Gertsenshtein mixing length 2c^2/(sqrt(G) B), B=1 G [Mpc]", Lmix / 100 / Mpc)
Dgc = 8 * kpc * 100
p("P photon->graviton, B=5 uG, D=8 kpc", math.sin(math.sqrt(Gcgs) * 5e-6 * Dgc / (2 * ccgs**2))**2)
# Pshirkov Eq 1 check from SI 4 pi G B^2 L^2/(mu0 c^4)
p("4piGB^2L^2/(mu0 c^4), B=1 G, L=1 cm (Pshirkov 8.3e-50)", 4 * math.pi * G * (1e-4)**2 * (0.01)**2 / (mu0 * c**4))
p("Pshirkov P_g->gamma B=1uG, L=1 kpc, vacuum coherence", 8.3e-50 * (1e-6)**2 * (kpc * 100)**2)
p("Pshirkov plasma L_coh, f=1e10 Hz, n_e=1e-6 cm^-3 [Mpc]", 3e14 / 1e-6 / 100 / Mpc)
p("plasma frequency 8.98 kHz sqrt(n_e), n_e=0.03 cm^-3 [Hz]", 8980 * math.sqrt(0.03))

# 5 speed of gravity
D17 = 26 * Mpc
p("dt c / D, dt=1.74 s, D=26 Mpc (upper bound)", 1.74 * c / D17)
p("(10-1.74) c/D, intrinsic delay 10 s (lower bound)", (10 - 1.74) * c / D17)
p("H0 bright siren v_H/d = 3017/43.8 [km/s/Mpc]", 3017 / 43.8)

# 6 graviton mass
for m in (1.92e-23, 2.23e-23, 2.42e-23, 1.27e-23, 1.01e-24):
    p(f"Compton wavelength h/(m c) for m_g={m:g} eV [m]", hP * c / (m * eV))
m = 1.92e-23; f = 100.0; D = 1e3 * Mpc
p("massive-graviton delay (D/2c)(mc^2/hf)^2, D=1 Gpc, 100 Hz [s]", D / (2 * c) * (m * eV / (hP * f))**2)

# 7 memory
def h_mem_0pn(m1, m2, D, x, theta=math.pi / 2):
    M = (m1 + m2) * GMsun / c**2; eta = m1 * m2 / (m1 + m2)**2
    s2 = math.sin(theta)**2; c2 = math.cos(theta)**2
    return 2 * eta * M * x / D * s2 * (17 + c2) / 96
p("eta M/R GW150914 (36+29, 410 Mpc)", 36 * 29 / 65**2 * 65 * GMsun / c**2 / D150914)
p("Favata 0PN memory at x=1/6, edge-on, GW150914", h_mem_0pn(36, 29, D150914, 1 / 6))
p("Favata 0PN memory at x=1/6, edge-on, 1e9+1e9 Msun at 1 Gpc", h_mem_0pn(1e9, 1e9, 1e3 * Mpc, 1 / 6))
p("PTA memory ramp B=1/2, h0=1e-14, 5 yr [us]", 0.5 * 1e-14 * 5 * yr * 1e6)

# 8 lensing (Takahashi & Nakamura 2003)
def w_lens(ML, f): return 8 * math.pi * G * ML * Msun * f / c**3
p("w = 8 pi G M f/c^3, M=100 Msun, f=100 Hz", w_lens(100, 100))
p("M_L with w=1 at f=100 Hz [Msun]", 1 / w_lens(1, 100))
wv = 1.0
p("|F|^2 at y=0, w=1: pi w/(1-exp(-pi w))", math.pi * wv / (1 - math.exp(-math.pi * wv)))
y = 1.0
p("mu+ at y=1", 0.5 + (y**2 + 2) / (2 * y * math.sqrt(y**2 + 4)))
p("mu- at y=1", 0.5 - (y**2 + 2) / (2 * y * math.sqrt(y**2 + 4)))
s = math.sqrt(y**2 + 4)
dt = 4 * 1e8 * Tsun * (y * s / 2 + math.log((s + y) / (s - y)))
p("time delay y=1, M_Lz=1e8 Msun [s]", dt)
p("4 G M/c^3, M=1e8 Msun [s]", 4 * 1e8 * Tsun)

# 9 absorption / dispersion by matter
for name, eta in (("water 1e-3 Pa s", 1e-3), ("air 1.8e-5 Pa s", 1.8e-5)):
    g = 16 * math.pi * G * eta / c**2
    p(f"Hawking/Baym damping 16 pi G eta/c^2, {name} [1/s]", g)
    p(f"  absorption length c/gamma, {name} [Gpc]", c / g / (1e3 * Mpc))
p("Baym 1/(omega tau_U), omega=1e3, tau_U=1e18 s", 1 / (1e3 * 1e18))
p("Weinberg neutrino damping amplitude factor sqrt(1-0.356)", math.sqrt(1 - 0.356))

# 10 resonant bar (Aguiar Eq 4.1) and GW170817 energy deposit
def sigma_bar(M, v, L, R, pois, n=1, theta=math.pi / 2):
    return 8 / math.pi * M / n**2 * (G / c) * (v / c)**2 * math.sin(theta)**4 * (1 + pois * (1 - 2 * pois) * (n * math.pi * R / L)**2 / 2)
Sig = sigma_bar(2300, 5400, 3.0, 0.6, 0.345)
p("bar integrated cross-section LSU 2300 kg [m^2 Hz]", Sig)
p("bar physical area / cross-section (2R*L)", 2 * 0.6 * 3.0 / Sig)
def dEdf(Mc_sun, f):  # Newtonian inspiral spectrum, J/Hz
    return (math.pi**(2 / 3) / 3) * G**(2 / 3) * (Mc_sun * Msun)**(5 / 3) * f**(-1 / 3)
D17b = 40 * Mpc; f0 = 900.0
Ff = dEdf(1.188, f0) / (4 * math.pi * D17b**2)
p("GW170817 spectral fluence at 900 Hz, isotropic [J/m^2/Hz]", Ff)
Eabs = Sig * Ff
p("energy absorbed by LSU bar (optimal orientation) [J]", Eabs)
p("  in quanta hbar*omega at 900 Hz", Eabs / (hbar * 2 * math.pi * f0))
p("  ratio to k_B*(0.1 K)", Eabs / (kB * 0.1))
p("  momentum E/c [kg m/s]", Eabs / c)
p("  bar velocity kick E/(cM) [m/s]", Eabs / (c * 2300))
p("tidal accel 0.5 w^2 h L, h=1e-21, 250 Hz, L=1 m [m/s^2]", 0.5 * (2 * math.pi * 250)**2 * 1e-21)
p("LIGO arm change h L/2, h=1e-21, L=4 km [m]", 1e-21 * 4000 / 2)
p("GW150914-like at 8 kpc: strain", 1e-21 * 410e3 / 8)
p("GW150914-like at 8 kpc: peak flux [W/m^2]", L150914 / (4 * math.pi * (8 * kpc)**2))
p("  ratio to solar constant 1361 W/m^2", L150914 / (4 * math.pi * (8 * kpc)**2) / 1361)
def h_from_lum(wp, f, r):  # luminosity wp*c^5/G -> sinusoidal h0 via T^tz = c^3 h0^2 w^2/(32 pi G)
    return 2 * math.sqrt(2) * math.sqrt(wp) * c / (2 * math.pi * f * r)
p("h0 bound from L = 1e-3 c^5/G, 250 Hz, 410 Mpc", h_from_lum(1e-3, 250, D150914))
p("h0 bound from L = c^5/G at one wavelength (r = lambda)", h_from_lum(1.0, 1.0, c / 1.0))
p("kick v/c GW200129 1542 km/s", 1542e3 / c)

# 11 inflation
As = math.exp(3.044) * 1e-10
Mpl = math.sqrt(hbar * c / (8 * math.pi * G)) * c**2 / eV / 1e9  # reduced Planck mass in GeV
p("A_s from ln(1e10 As)=3.044", As); p("reduced M_Pl [GeV]", Mpl)
for r in (0.036, 0.034, 0.032, 0.01, 0.001):
    V14 = (1.5 * math.pi**2 * As * r)**0.25 * Mpl
    H = math.pi * Mpl * math.sqrt(As * r / 2)
    p(f"V^1/4 r={r} [GeV]", V14); p(f"H_inf r={r} [GeV]", H)
Orad = 2.47e-5 * (1 + 3.046 * (7 / 8) * (4 / 11)**(4 / 3))
p("Omega_rad h^2 (photons 2.47e-5 + 3.046 nu)", Orad)
p("Omega_GW h^2 = Orad r As/24, r=0.036 (no g* or nu damping)", Orad * 0.036 * As / 24)

# 12 binary pulsars (Peters-Mathews)
def pbdot(mp_, mc_, Pb_days, e):
    fe = (1 + 73 / 24 * e**2 + 37 / 96 * e**4) / (1 - e**2)**3.5
    Pb = Pb_days * 86400
    return -192 * math.pi / 5 * Tsun**(5 / 3) * (Pb / (2 * math.pi))**(-5 / 3) * fe * mp_ * mc_ * (mp_ + mc_)**(-1 / 3)
p("Pbdot GR B1913+16 (1.438,1.390,0.322997448918 d,0.6171340)", pbdot(1.438, 1.390, 0.322997448918, 0.6171340))
p("Pbdot GR J0737-3039A/B (1.338185,1.248868,0.1022515592973 d,0.087777023)", pbdot(1.338185, 1.248868, 0.1022515592973, 0.087777023))
p("Peters f(e) Hulse-Taylor", (1 + 73 / 24 * 0.617134**2 + 37 / 96 * 0.617134**4) / (1 - 0.617134**2)**3.5)

# 13 multipole coefficients (Blanchet LRR Table 1)
def dfact(n): return math.prod(range(n, 0, -2))
for l in (2, 3, 4):
    a = (l + 1) * (l + 2) / ((l - 1) * l * math.factorial(l) * dfact(2 * l + 1))
    b = 4 * l * (l + 2) / ((l - 1) * math.factorial(l + 1) * dfact(2 * l + 1))
    p(f"energy-flux a_{l} (mass)", a); p(f"energy-flux b_{l} (current)", b)
    p(f"  1/a_{l}", 1 / a); p(f"  1/b_{l}", 1 / b)

# 14 Dyson 1963 eq (11) check
Mcgs = 2e33; R = 2e9; Gc = 6.7e-8; cc = 3e10
V = math.sqrt(Gc * Mcgs / (4 * R))
p("Dyson V = sqrt(GM/4R) [cm/s]", V)
p("Dyson W = 128 V^10/(5 G c^5) [erg/s]", 128 * V**10 / (5 * Gc * cc**5))
a_sep = 2 * R
p("Peters equal-mass (64/5)G^4M^5/(c^5 a^5) [erg/s]", 64 / 5 * Gc**4 * Mcgs**5 / (cc**5 * a_sep**5))
p("(2048/5) V^10/(G c^5) [erg/s]", 2048 / 5 * V**10 / (Gc * cc**5))
p("Dyson period 2 pi R/V [s]", 2 * math.pi * R / V)

# self-checks against published numbers (loose tolerances)
def near(k, ref, tol): assert abs(out[k] / ref - 1) < tol, (k, out[k], ref)
near("Flanagan-Hughes flux h0=1e-22, f=1kHz [mW/m^2]", 1.5, 0.1)
near("Tobar N = h0^2 c^5/(32 pi G hbar nu^2), h0=1e-21, 150 Hz", 4e36, 0.25)
near("Tobar Gamma_stim Al 1800 kg, vs=5.4 km/s, h=5e-22 [Hz]", 1.0, 0.2)
near("Carney n*lambda_dB^3, h=1e-22, f=1 kHz", 2e35, 0.25)
near("4piGB^2L^2/(mu0 c^4), B=1 G, L=1 cm (Pshirkov 8.3e-50)", 8.3e-50, 0.02)
near("Pbdot GR B1913+16 (1.438,1.390,0.322997448918 d,0.6171340)", -2.40263e-12, 1e-3)  # masses quoted to 1e-3 only
near("Dyson W = 128 V^10/(5 G c^5) [erg/s]", 2e37, 0.15)
near("GW150914 200 Msun c^2/s [W]", 3.5e49, 0.05)
near("bar integrated cross-section LSU 2300 kg [m^2 Hz]", 4e-25, 0.1)
print("self-checks passed")

# 15 angular- and linear-momentum flux coefficients (Blanchet LRR Table 1), l=2
l = 2
ap = (l + 1) * (l + 2) / ((l - 1) * math.factorial(l) * dfact(2 * l + 1))
app = 2 * (l + 2) * (l + 3) / (l * math.factorial(l + 1) * dfact(2 * l + 3))
cpp = 8 * (l + 2) / ((l - 1) * math.factorial(l + 1) * dfact(2 * l + 1))
p("angular-momentum a'_2 (expect 2/5)", ap); p("linear-momentum a''_2 (expect 2/63)", app); p("linear-momentum c''_2 (expect 16/45)", cpp)
assert abs(ap - 0.4) < 1e-12 and abs(app - 2 / 63) < 1e-12 and abs(cpp - 16 / 45) < 1e-12
# circular binary 0PN flux in terms of x: (32/5)(c^5/G) eta^2 x^5 ; equal mass at ISCO x=1/6
p("0PN flux equal mass at x=1/6, in units of c^5/G", 32 / 5 * 0.25**2 * (1 / 6)**5)
print("momentum-flux checks passed")

# 16 extra numbers quoted in the notes
p("l_P^2 / sigma_Thomson (graviton vs photon cross-section)", lP2 / 6.6524587321e-29)
p("thrust per watt of a collimated GW beam 1/c [N/W]", 1 / c)
p("tidal accel per metre, GW150914-like at 8 kpc [m/s^2]", out["tidal accel 0.5 w^2 h L, h=1e-21, 250 Hz, L=1 m [m/s^2]"] * 410e3 / 8)
g_mantle = 16 * math.pi * G * 1e21 / c**2
p("misapplied mantle viscosity 1e21 Pa s: gamma [1/s]", g_mantle); p("  c/gamma [m] (wrong regime)", c / g_mantle)
