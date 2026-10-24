"""Checks for the review of research/external/2026-10-24-cross-crate-frontier-regimes.md.
Standard library only. Run from outside the repository: python3 -I <path>."""
# Batch 1, item 1: can a site-dependent (geomagnetic or geopotential) bias explain the G scatter?
from math import pi, sqrt, radians, exp
G = [("NIST-82",6.67248,0.00043),("TR&D-96",6.6729,0.0005),("LANL-97",6.67398,0.00070),
     ("UWash-00",6.674255,0.000092),("BIPM-01",6.67559,0.00027),("UWup-02",6.67422,0.00098),
     ("MSL-03",6.67387,0.00027),("HUST-05",6.67222,0.00087),("UZur-06",6.67425,0.00012),
     ("HUST-09",6.67349,0.00018),("BIPM-14",6.67554,0.00016),("LENS-14",6.67191,0.00099),
     ("UCI-14",6.67435,0.00013),("HUSTT-18",6.674184,0.000078),("HUSTA-18",6.674484,0.000077),
     ("JILA-18",6.67260,0.00025)]
mean = 6.674300396
d = {n:(v,u) for n,v,u in G}
def ppm(x): return (x/mean)*1e6
# Same site, different apparatus: HUST (Wuhan) has four inputs.
hust = [d[k] for k in ("HUST-05","HUST-09","HUSTT-18","HUSTA-18")]
vals=[v for v,_ in hust]
print("HUST spread (ppm):", ppm(max(vals)-min(vals)))
t,a = d["HUSTT-18"], d["HUSTA-18"]
print("HUSTT-18 vs HUSTA-18, same lab, same year: diff ppm", ppm(a[0]-t[0]), "sigma", (a[0]-t[0])/(t[1]**2+a[1]**2)**0.5)
b1,b2 = d["BIPM-01"], d["BIPM-14"]
print("BIPM-01 vs BIPM-14, same site:", ppm(b1[0]-b2[0]), "ppm,", (b1[0]-b2[0])/(b1[1]**2+b2[1]**2)**0.5, "sigma")
# Magnetic force on a susceptible test mass in Earth's own field: F = chi V/mu0 * B dB/dx.
mu0 = 1.25663706212e-6
B = 5e-5; dBdx = 3*B/6.371e6                 # dipole gradient scale 3B/r
chi_au, V = -3.4e-5, 0.2/19300               # 0.2 kg gold sphere
F_mag = abs(chi_au)*V/mu0*B*dBdx
F_grav = 6.6743e-11*10*0.2/0.1**2            # 10 kg source at 10 cm
print("Earth-field gradient force on 0.2 kg gold (N):", F_mag, " Cavendish signal (N):", F_grav, " ratio:", F_mag/F_grav)
# What lab field gradient would bias G by the CODATA spread (2.2e-5)?
need = 2.2e-5*F_grav/(abs(chi_au)*V/mu0)    # needed B*dB/dx, T^2/m
print("B dB/dx needed for a 2.2e-5 bias (T^2/m):", need, "e.g. B = 1e-4 T with dB/dx (T/m):", need/1e-4)

# Gravitomagnetic force on a test mass at 1 mm/s against the Cavendish acceleration.
Bg = 6.6743e-11*5.86e33/299792458.0/6.371e6**3
print("gravitomagnetic: 2 v B/c at 1 mm/s (m/s^2):", 2*1e-3/299792458.0*Bg, " fraction of 6.7e-8:", 2*1e-3/299792458.0*Bg/(F_grav/0.2))

# ---- Item 2: lab high-frequency gravitational waves ----
G=6.67430e-11; C=299792458.0; HBAR=1.054571817e-34; MU0=1.25663706212e-6; H=6.62607015e-34
k = G/C**4
print("G/c^4 (s^2/(kg m)):", k, " 8 pi G/c^4:", 8*pi*k)
# NIF implosion: shell kinetic energy ~ 20 kJ, hot-spot stagnation in ~100 ps, P2 asymmetry eps.
E_kin, tau = 2e4, 1e-10
for eps in (0.01, 0.1):
    i_ddot = 2*eps*E_kin                          # second time derivative of the quadrupole, kg m^2/s^2 (= J)
    h1m = 2*k*i_ddot/1.0
    e_gw = G/(5*C**5)*(i_ddot/tau)**2*tau         # (G/5c^5) <I'''^2> over tau
    f = 1/tau
    print(f"NIF eps={eps}: h at 1 m {h1m:.3e}; energy radiated {e_gw:.3e} J; gravitons at {f:.0e} Hz: {e_gw/(H*f):.3e}")
# The fusion yield itself (8.6 MJ), emitted with a small anisotropy eps over the burn:
for eps in (0.01,):
    print(f"8.6 MJ neutron burst, eps={eps}: memory-like h at 1 m {4*k*eps*8.6e6:.3e}")
# Gertsenshtein conversion in a SPARC-like 20 T field over 2 m, and a CoReLS-like pulse through it.
P = 4*pi*G*20**2*2**2/(MU0*C**4)
print("P(gamma->g), 20 T over 2 m:", P)
photons = 4.0/(H*C/0.8e-6)                         # a 4 J pulse at 0.8 um
print("photons in a 4 J, 0.8 um pulse:", photons, " gravitons produced over 2 m of 20 T:", photons*P)
# A 150 GA, 1 ns Z-pinch (the batch's upper current): magnetic energy per metre mu0 I^2/(4 pi) ln(b/a) ~ mu0 I^2/(4 pi) for ln = 1
I = 150e9
print("magnetic energy per metre at 150 GA (J/m, ln=1):", MU0*I**2/(4*pi))

# ---- Item 3: buried bodies, sensor averaging, DOV in navigation ----
G=6.67430e-11; g0=9.80665; R=6.371e6
EOTVOS=1e-9
# A buried horizontal cylinder (a tunnel), radius a, depth z of its axis, density contrast drho:
def cylinder(a, z, drho, x):
    lam = drho*pi*a*a
    gz = 2*G*lam*z/(x*x+z*z)
    gzz = 2*G*lam*(z*z-x*x)/(x*x+z*z)**2
    gxx = -gzz                                  # 2-D Laplace: gxx + gzz = 0
    return gz, gzz, gxx
gz, gzz, gxx = cylinder(1.0, 1.9, -2000.0, 0.0)
print("tunnel a=1 m, axis 1.9 m, drho=-2000: dg_z (uGal)", gz/1e-8, " G_zz (E)", gzz/EOTVOS)
# A buried sphere: dg_z and G_zz above the centre.
def sphere(a, z, drho, x):
    m = drho*4/3*pi*a**3
    r2 = x*x+z*z
    return G*m*z/r2**1.5, G*m*(2*z*z-x*x)/r2**2.5
s = sphere(5.0, 20.0, -2700.0, 0.0)
print("sphere a=5 m at 20 m, drho=-2700: dg_z (uGal)", s[0]/1e-8, " G_zz (E)", s[1]/EOTVOS)
# Averaging a white-noise gravimeter: 500 nm/s^2/sqrt(Hz) (Exail AQG datasheet)
asd = 500e-9
for t in (90, 600, 2400):
    print(f"AQG averaged {t} s: {asd/sqrt(t)*1e9:.1f} nm/s^2")
print("time to 2 uGal (20 nm/s^2):", (asd/20e-9)**2, "s")
# Equivalence principle: one accelerometer reads a - g; a gradiometer's two proof masses a baseline L apart
# read the same vehicle acceleration and differ by the tidal tensor times L.
L = 0.5
print("gradiometer: 1 E over 0.5 m is a differential acceleration (m/s^2):", 1*EOTVOS*L)
# Deflection of the vertical theta: an inertial navigator mistakes g*theta for horizontal acceleration;
# Schuler tuning bounds the position error at peak 2 R theta (steady R theta).
for arcsec in (1, 5, 10):
    th = radians(arcsec/3600)
    print(f"DOV {arcsec}\": horizontal gravity error {g0*th*1e6:.1f} um/s^2 (= {g0*th/1e-8:.0f} uGal); Schuler peak position error {2*R*th:.0f} m")
print("Schuler period (min):", 2*pi*sqrt(R/g0)/60)

# ---- Item 4: PEGS timing, AGW growth, ionospheric dynamo ----
G=6.67430e-11; C=299792458.0; g0=9.80665
# PEGS: a point mass change dM at depth, the prompt (speed-of-light) gravity change at range r.
# Monopole term dg ~ G dM / r^2; the actual coseismic signal is quadrupolar and smaller.
# Tohoku Mw 9.07: M0 = 4.59e22 N m. A representative mass redistribution: dM ~ M0/(c_rupture^2 * ...);
# use the published PEGS amplitude scale instead of a toy dM.
M0 = 4.59e22
# Vallee et al 2017 and Juhel et al 2018 report PEGS of ~ 0.1-1 nm/s^2 at ~1000 km, tens of seconds
# before P. Here just show the monopole scale for a plausible dM and the arrival-time gap.
for dM, r in ((1e13, 1e6), (1e14, 1e6)):
    print(f"monopole dg for dM={dM:.0e} kg at r={r:.0e} m: {G*dM/r**2/1e-9:.3f} nm/s^2")
vp, c_lamb = 7e3, 318.0
r = 1e6
print(f"at {r/1e3:.0f} km: light travel {r/C*1e3:.3f} ms; P-wave {r/vp:.0f} s; Lamb wave {r/c_lamb:.0f} s")
print(f"PEGS lead over P-wave at 1000 km: {r/vp - r/C:.0f} s")
# Acoustic-gravity wave amplitude growth: energy flux ~ rho v^2 c conserved; rho falls exponentially
# with scale height H ~ 7 km, so v grows as exp(z/2H).
H = 7e3
for z in (100e3, 300e3):
    print(f"AGW velocity amplification by z={z/1e3:.0f} km (H=7 km): factor {exp(z/(2*H)):.3e}")
# Ionospheric dynamo: neutral wind u across B drives J = sigma (u x B); the dynamo EMF field is u*B.
# Hunga Tonga Lamb wave: surface particle velocity v0 ~ p/(rho0 c) for a 780 Pa pulse.
p0, rho0 = 780.0, 1.2
v0 = p0/(rho0*c_lamb)
print(f"Lamb-wave surface particle velocity v0 = p/(rho0 c): {v0:.2f} m/s")
# At 150 km the neutral velocity is amplified; take u ~ v0 * exp(z/2H) capped by dissipation, B ~ 3e-5 T.
B = 3e-5
u150 = v0*exp(100e3/(2*H))   # use 100 km growth as a conservative cap
E_dyn = u150*B
print(f"dynamo field u*B at ~120-150 km (u~{u150:.0f} m/s, B=30000 nT): {E_dyn*1e3:.2f} mV/m")
# Alfven speed in the F-region: B/sqrt(mu0 rho), rho ~ n_i m_i, n_i ~ 1e12 m^-3 of O+ (16 amu)
mu0=1.25663706212e-6; amu=1.66053907e-27
rho_iono = 1e12*16*amu
vA = B/sqrt(mu0*rho_iono)
print(f"F-region Alfven speed (n=1e12 O+, B=30000 nT): {vA/1e3:.0f} km/s; transit to conjugate point ~ few 1000 km:", f"{2e6/vA:.0f} s")

# ---- Item 5: condensate suppression, mid-band, atom phases ----
G=6.67430e-11; C=299792458.0; HBAR=1.054571817e-34; H=6.62607015e-34; ME=9.1093837015e-31; E=1.602176634e-19
# DeWitt gravito-magnetic coupling suppression for a superconductor: G m_e^2/(hbar c).
supp = G*ME**2/(HBAR*C)
print("G m_e^2/(hbar c):", supp)
# Why Podkletnov-scale weight change is unobservable: the GEM permeability is 4 pi G/c^2.
print("gravito-magnetic permeability 4 pi G/c^2 (m/kg):", 4*pi*G/C**2)
# Ultralight dark matter: Compton frequency f = m c^2 / h for a scalar of mass m (eV).
for m_eV in (1e-15, 1e-14, 1e-13):
    f = m_eV*E/H
    print(f"ULDM m={m_eV:.0e} eV: Compton frequency {f:.3e} Hz (mid-band {0.01} - 10 Hz?) ")
# The mid-band blind spot between LISA (<~0.1 Hz) and LIGO (>~10 Hz) maps to masses:
for f in (0.01, 1.0, 10.0):
    print(f"f={f} Hz <-> scalar mass {H*f/E:.3e} eV")
# Atom-interferometer gravimeter phase: dphi = k_eff g T^2 (Kasevich-Chu). Rosi 2014: keff=2*2pi/780nm, T=0.16s
keff = 2*2*pi/780.241e-9; T=0.160
print("Rosi k_eff T^2 (rad per m/s^2):", keff*T**2, " phase for g0:", keff*T**2*9.80665)
# Gravitational Aharonov-Bohm: phase from the potential around a mass M at distance r for Rb at 7 cm (Overstreet 2022)
mRb = 86.909*1.66053907e-27
print("Overstreet-scale phase rate m*GM/(r hbar), M=1.25 kg r=0.07 m (rad/s):", mRb*G*1.25/0.07/HBAR)
