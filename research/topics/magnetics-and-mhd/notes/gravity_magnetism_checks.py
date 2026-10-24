"""Recompute the test values in gravity_magnetism_links_and_claims.md. SI units; Gaussian/HL conversions noted.
Run: python3 -I gravity_magnetism_checks.py   (standard library only; asserts are the self-test)."""
from math import pi, sqrt, cos, radians
G = 6.67430e-11; c = 2.99792458e8; hbar = 1.054571817e-34; h = 2*pi*hbar
e = 1.602176634e-19; me = 9.1093837015e-31; mp = 1.67262192e-27
mu0 = 1.25663706212e-6; eps0 = 8.8541878128e-12; alpha = 1/137.035999
muB = 9.2740100783e-24; g0 = 9.80665; Msun = 1.98892e30; yr = 3.15576e7
MAS = 180/pi*3600e3            # mas per rad
eV = e; GeV = 1e9*eV; hbarc_eVm = hbar*c/eV   # 1.9733e-7 eV m

def p(name, val, unit=""):
    print(f"{name:66s} {val:.4g} {unit}")

def close(a, b, rel):
    assert abs(a/b - 1) < rel, (a, b)

print("== 1. Gravitoelectromagnetism: conventions and Earth ==")
M_E, R_E, w_E = 5.9722e24, 6.3781e6, 7.292115e-5
J_E = 8.034e37*w_E                  # C = 0.3307 M R^2 (assumed IERS-type value)
p("Earth angular momentum J", J_E, "kg m^2/s (Biemond table: 5.87e33)")
B_M_eq = G*J_E/(c*R_E**3)            # Mashhoon B = curl A, A = (G/c) J x r / r^3  [m/s^2]
p("Mashhoon B at equator, G J/(c R^3)", B_M_eq, "m/s^2")
p("Mashhoon B at pole, 2 G J/(c R^3)", 2*B_M_eq, "m/s^2")
p("Gyro precession Omega_P = B/c at equator", B_M_eq/c, "rad/s")
B_RT_eq = -2*B_M_eq                  # Ruggiero-Tartaglia eq.(24) dipole = -2 x curl of their eq.(10) A
B_W_eq = G*J_E/(2*c**2*R_E**3)       # Wikipedia-SI convention, units 1/s, force m(E + 4 v x B)
p("Wikipedia-SI B_g at equator, G J/(2 c^2 R^3)", B_W_eq, "1/s")
v = 7.9e3                            # orbital speed
f_M = 2*(v/c)*B_M_eq; f_RT = (v/c)*abs(B_RT_eq); f_W = 4*v*B_W_eq
p("v x B force/mass, Mashhoon 2(v/c)B", f_M, "m/s^2"); p("  R&T (v/c)B_G", f_RT, "m/s^2"); p("  Wikipedia 4 v B_g", f_W, "m/s^2")
close(f_M, f_RT, 1e-12); close(f_M, f_W, 1e-12)
p("Wikipedia uniform-sphere Earth formula 2 pi R g/(5 c^2 T) / exact", (2*pi*R_E*G*M_E/R_E**2/(5*c**2*2*pi/w_E))/B_W_eq, "(=0.4/0.3307)")
mu_g = 4*pi*G/c**2; eps_g = 1/(4*pi*G)
p("gravitomagnetic permeability mu_g = 4 pi G/c^2 (force factor 4)", mu_g, "m/kg")
p("  same with B'=4B_g (force factor 1): 16 pi G/c^2", 4*mu_g, "m/kg")
p("gravitoelectric permittivity 1/(4 pi G)", eps_g, "kg s^2/m^3"); close(1/sqrt(mu_g*eps_g), c, 1e-12)
p("hbar*Omega_P, spin-1/2 at pole (Mashhoon ~1e-29 eV)", hbar*2*B_M_eq/c/eV, "eV")
p("hbar*Omega_Earth (Mashhoon ~1e-19 eV)", hbar*w_E/eV, "eV")
def lt_node(a, ecc): return 2*G*J_E/(c**2*a**3*(1-ecc**2)**1.5)*yr*MAS
sats = {"LAGEOS": (12270e3, 0.0045), "LAGEOS 2": (12163e3, 0.014), "LARES": (7821e3, 0.0008), "LARES 2": (12266e3, 0.0001)}
for k, (a, ecc) in sats.items(): p(f"Lense-Thirring node rate {k} (elements assumed)", lt_node(a, ecc), "mas/yr")
s = lt_node(*sats["LAGEOS"]) + lt_node(*sats["LARES 2"])
p("LAGEOS + LARES 2 node sum (Nature 2026 reports 61.3)", s, "mas/yr"); close(s, 61.3, 0.01)
a_gpb = R_E + 642e3
fd = G*J_E/(2*c**2*a_gpb**3)*yr*MAS*cos(radians(16.84))   # polar-orbit average, projected for IM Peg dec
p("GP-B frame dragging, polar orbit, IM Peg projection (GR 39.2)", fd, "mas/yr"); close(fd, 39.2, 0.02)
p("GP-B measured 37.2 +- 7.2: (meas-GR)/sigma", (37.2-39.2)/7.2)

print("\n== 2. Blackett / Schuster-Wilson ==")
k_g = sqrt(4*pi*eps0*G)      # 'gravitational charge-to-mass', = RN extremal Q/M
p("sqrt(4 pi eps0 G) (C/kg)", k_g)
def beta_SI(mu, J): return 2*mu/(k_g*J)        # mu = (beta/2) sqrt(4 pi eps0 G) J  == Gaussian mu = beta sqrt(G) J/(2c)
Gcgs, ccgs = 6.67430e-8, 2.99792458e10
def beta_cgs(mu, J): return mu/(sqrt(Gcgs)*J/(2*ccgs))
close(beta_SI(7.91e22, 5.87e33), beta_cgs(7.91e25, 5.87e40), 1e-3)   # emu -> A m^2 is 1e-3, g cm^2/s -> SI is 1e-7
# Biemond (viXra 2006.0221) Table 1 values, Gaussian (G cm^3, g cm^2/s)
bodies = {"Earth": (7.91e25, 5.87e40), "Sun": (1.35e33, 1.57e48), "78 Vir": (2.1e36, 4.2e50), "Jupiter": (1.35e30, 4.29e45),
          "Saturn": (4.54e28, 7.71e44), "Venus": (2.99e22, 1.81e38), "Moon": (2.26e21, 2.36e36),
          "96 isolated pulsars (mean)": (2.6e30, 1.0e46), "3 binary MSPs (mean)": (4.5e27, 1.6e47), "10 isolated WDs (mean)": (1.6e34, 1.2e47)}
bs = {}
for k, (mu, J) in bodies.items():
    bs[k] = beta_cgs(mu, J); p(f"beta {k}", bs[k])
I_mars = 0.4*0.925*6.4e26*(3.39e8)**2; J_mars = I_mars*2*pi/(1.03*86400)
p("Mars: beta upper limit from mu < 2e20 G cm^3 (Campanelli: 2e-5)", beta_cgs(2e20, J_mars))
p("spread max/min of beta over the table", max(bs.values())/min(bs.values()))
p("electron: beta (Blackett, g=2 spin) = 2 e/(m_e sqrt(4 pi eps0 G))", 2*e/(me*k_g))
p("electron: Barrow-Gibbons beta = N^1/2 = e/(m_e sqrt(4 pi eps0 G))", e/(me*k_g))
sig, rho, R = 1510e6, 19250.0, 1.0      # tungsten sphere, R = 1 m, at break-up spin
Bmax_G = 1e-13*sqrt(sig/1e6)*sqrt(rho/1e3)*R   # Campanelli eq.5 scaling, gauss
p("Campanelli eq.5: max Blackett field, W sphere R=1 m, beta=1e-5", Bmax_G*1e-4, "T")

print("\n== 3. Magnetic monopoles ==")
qD = h/(mu0*e)
p("Dirac charge g_D = h/(mu0 e) (SI, A m)", qD, "A m"); p("  = e/(2 alpha) in Gaussian, in units of e", 1/(2*alpha))
p("flux of one Dirac monopole mu0 g_D = h/e = 2 Phi0", mu0*qD, "Wb")
F_cab = 1/(151*86400*20*2*pi)
p("Cabrera 1982: 1 event / (151 d x 20 cm^2 x 2 pi sr) (PRL 6.1e-10)", F_cab, "cm^-2 s^-1 sr^-1"); close(F_cab, 6.1e-10, 0.01)
for name, lim in (("Parker 1e-15", 1e-15), ("MACRO 1.4e-16", 1.4e-16), ("IceCube 2022 2e-19", 2e-19)):
    p(f"Cabrera flux / {name}", F_cab/lim)
def parker(M_GeV): return 1e-15*max(1.0, M_GeV/1e17)
def parker_seed(M_GeV): return (M_GeV/1e17 + 3e-6)*1e-16
for M in (1e10, 1e15, 1e17, 1e19): p(f"Parker (PDG 94.9) / seed bound (94.10) at M={M:.0e} GeV", parker(M), f"/ {parker_seed(M):.3g}")
p("GUT monopole mass ~ M_X/alpha_GUT, M_X=1e16 GeV, alpha=1/40", 1e16*40, "GeV (PDG: 1e17-1e18)")
qm_ice = 2*10*muB/(sqrt(1.5)*3.54e-10)
p("spin-ice monopole q_m = 2 mu/a_d (Dy2Ti2O7, mu=10 muB)", qm_ice, "A m"); p("  q_D / q_m (CMS 2008: ~8000)", qD/qm_ice)

print("\n== 4. Axions: ADMX haloscope, CAST helioscope, and Gertsenshtein in the same magnets ==")
T_nat = sqrt(1/mu0/eV*hbarc_eVm**3)   # 1 T in HL natural units, eV^2
p("1 tesla in Heaviside-Lorentz natural units", T_nat, "eV^2")
def g_axion(f_Hz, ggam):             # g_agg = (ggam alpha/pi) m_a / Lambda^2, Lambda = 78 MeV
    m = h*f_Hz/eV; return ggam*alpha/pi*m/(78e6)**2, m      # eV^-1
def haloscope_W(V_L, B_T, C, ggam, rho_GeVcc, f_Hz, Q):
    g, m = g_axion(f_Hz, ggam)
    rho = rho_GeVcc*1e9*(hbarc_eVm*100)**3                 # eV^4
    Vn = V_L*1e-3/hbarc_eVm**3                             # eV^-3
    P = g**2*rho*(B_T*T_nat)**2*Vn*C*Q/m                   # eV^2
    return P*eV*(eV/hbar)
P_admx = haloscope_W(136, 6.8, 0.4, 0.97, 0.45, 650e6, 50000)
p("ADMX 2018 reference: P = g^2 (rho/m) B^2 V C Q (paper 1.9e-22 W)", P_admx, "W"); close(P_admx, 1.9e-22, 0.05)
p("  KSVZ g_agg at 650 MHz", g_axion(650e6, 0.97)[0]*1e9, "GeV^-1"); p("  DFSZ power at same settings", haloscope_W(136, 6.8, 0.4, 0.36, 0.45, 650e6, 50000), "W")
p("  photons/s for KSVZ at 650 MHz", P_admx/(h*650e6), "1/s")
B, L = 9.0, 9.26
for yrlab, g in (("CAST 2017", 0.66e-10), ("CAST 2024", 0.58e-10)):
    p(f"{yrlab} P(a->gamma) = (g B L/2)^2 at limit g={g:.2g}/GeV", (g/1e9*B*T_nat*L/hbarc_eVm/2)**2)
P_gw = 4*pi*G*B**2*L**2/(mu0*c**4)
p("Gertsenshtein P(g->gamma) = 4 pi G B^2 L^2/(mu0 c^4), CAST magnet", P_gw)
MPl = sqrt(hbar*c/G)*c**2/GeV
p("equivalent axion coupling of the graviton sqrt(16 pi G) = sqrt(16 pi)/M_Pl", sqrt(16*pi)/MPl, "GeV^-1")
close((sqrt(16*pi)/MPl/1e9*B*T_nat*L/hbarc_eVm/2)**2, P_gw, 1e-3)
p("CAST-2024 coupling / graviton-equivalent coupling", 0.58e-10/(sqrt(16*pi)/MPl))

print("\n== 5. Magnetism in strong gravity ==")
def bz_SI(Phi_Wb, Om, k=1/(6*pi)): return k*Phi_Wb**2*Om**2/(mu0*c)   # = (k/4pi) Phi_G^2 Om^2/c (Gaussian) = k Phi_HL^2 Om^2 (HL)
Phi_G = 1e33; Om = 1e-5
close(bz_SI(Phi_G*1e-8, Om)*1e7, (1/(6*pi))/(4*pi)*Phi_G**2*Om**2/ccgs, 1e-9)  # SI vs Gaussian
M87 = 6.5e9*Msun; rg = G*M87/c**2
for a in (0.5, 0.9, 0.94):
    rH = rg*(1+sqrt(1-a*a)); OmH = a*c/(2*rH)
    for BG in (1, 30):
        Phi = 2*pi*rH**2*BG*1e-4
        p(f"M87* a={a}: B_H={BG} G -> BZ power (k=1/6pi)", bz_SI(Phi, OmH)*1e7, "erg/s")
    Breq = sqrt(1e42*1e-7/bz_SI(2*pi*rH**2*1.0, OmH))
    p(f"  B at horizon needed for 1e42 erg/s (EHT V minimum)", Breq*1e4, "G")
for Md in (3e-4, 2e-3): p(f"M87* Mdot c^2 at {Md:g} Msun/yr", Md*Msun/yr*c**2*1e7, "erg/s")
Bp, Rns, Mns = 1.6e11, 1e4, 1.4*Msun
Eext = pi*Bp**2*Rns**3/(3*mu0)
p("SGR 1806-20 external dipole energy pi B^2 R^3/(3 mu0) (2e47 erg)", Eext*1e7, "erg"); close(Eext*1e7, 2e47, 0.1)
p("2004 giant flare 2e46 erg / that", 2e46/(Eext*1e7))
p("B_p / Schwinger B_c = m_e^2 c^2/(e hbar)", Bp/(me**2*c**2/(e*hbar)))
p("virial max field: (B^2/2mu0)(4pi/3)R^3 = (3/5) G M^2/R", sqrt(0.6*G*Mns**2/Rns*2*mu0/(4*pi/3*Rns**3)), "T")
Esat = 4e50*1e-7; Vrem = 4/3*pi*(15e3)**3
p("merger: rms B for 4e50 erg in r=15 km (Kiuchi 2015)", sqrt(2*mu0*Esat/Vrem)*1e4, "G")
p("spacetime curvature radius of a 1 T field c^2 sqrt(mu0/(4 pi G))/B", c**2*sqrt(mu0/(4*pi*G)), "m")
p("  for magnetar B_p", c**2*sqrt(mu0/(4*pi*G))/Bp, "m")

print("\n== 6. Diamagnetic levitation ==")
def BdB(rho, chi): return mu0*rho*g0/abs(chi)
p("water, chi=-8.8e-6, g=9.81 (Berry-Geim eq.27: 1400.9)", mu0*1000*9.81/8.8e-6, "T^2/m"); close(mu0*1000*9.81/8.8e-6, 1400.9, 1e-3)
p("water, chi=-9.05e-6", BdB(1000, -9.05e-6), "T^2/m")
p("  chi_cgs = chi_SI/(4 pi)", -9.05e-6/(4*pi))
p("Bitter magnet: B0 from B=0.63 B0, dB/dz=-8.15 B0 /m (paper 16.5 T)", sqrt(1400.9/(0.63*8.15)), "T")
p("induced B inside frog |chi| B0", 8.8e-6*16.5*1e4, "G (paper ~1.5 G)")
p("pyrolytic graphite chi_perp=-4.5e-4, rho=2200", BdB(2200, -4.5e-4), "T^2/m")
p("  NdFeB surface 0.5 T falling to 0 over 4 mm gives ~", 0.5*0.5/4e-3, "T^2/m")
p("field energy density at 16.5 T", 16.5**2/(2*mu0), "J/m^3")

print("\n== 7. Fringe claims vs conservation ==")
Vmag = 0.040*0.0254*0.0381
p("MEG magnet volume (US 6,362,718)", Vmag, "m^3")
for lab, u in (("Br=1.2 T: Br^2/2mu0", 1.2**2/(2*mu0)), ("(BH)max 300 kJ/m^3", 3e5)):
    p(f"  stored energy, {lab}", u*Vmag, "J")
Pin, Pout = 14.0, 48.0
p("  claimed COP per output coil 48 W/14 W", Pout/Pin); p("  excess power, 2 coils", 2*Pout-Pin, "W")
p("  seconds to drain a 22 J magnet at that excess", 1.2**2/(2*mu0)*Vmag/(2*Pout-Pin), "s")
m_rem = 225.0; Vrem_RG = m_rem/7500.0
p("Roshchin-Godin magnet volume (225 kg, rho 7500 assumed)", Vrem_RG, "m^3")
p("  stored field energy Br^2/2mu0 V, Br=0.85 T", 0.85**2/(2*mu0)*Vrem_RG, "J")
p("  same at the paper's '150' read as 150 kJ/m^3", 150e3*Vrem_RG, "J")
KE = 0.5*115*(550*2*pi/60*0.5)**2
p("  roller-ring kinetic energy, 115 kg at r=0.5 m, 550 rpm", KE, "J"); p("  runtime at 7 kW", KE/7e3, "s")
F = 0.35*350*g0
p("  claimed lift 35% of 350 kg", F, "N"); p("  photon-rocket power F c", F*c, "W")
A = pi*0.5**2
p("  ideal actuator-disk power for F over 1 m disc, air", F**1.5/sqrt(2*1.2*A), "W")
p("  induced air speed sqrt(F/(2 rho A))", sqrt(F/(2*1.2*A)), "m/s")
m_dip = 0.85/mu0*Vrem_RG
p("  dipole field of the magnets at 15 m (claimed 0.05 T)", mu0*2*m_dip/(4*pi*15**3), "T")
Vwalls = sum(2*pi*r*0.065*2.5 for r in [0.55*i for i in range(1, 28)])
p("  'magnetic walls' energy: 0.05 T in shells to 15 m, h=2.5 m", 0.05**2/(2*mu0)*Vwalls, "J")
rb = (3*27.2/(4*pi*7870))**(1/3); need = 27.2*g0/(pi*rb**2)
p("Hutchison 60-lb iron ball: pressure needed under it", need, "Pa")
p("  electrostatic pressure at air breakdown eps0 E^2/2, 3 MV/m", eps0*(3e6)**2/2, "Pa")
p("  magnetic pressure B^2/2mu0 equal to it needs B", sqrt(2*mu0*need), "T")

print("\n== 8. Unification numbers ==")
R_KK = sqrt(16*pi*G*4*pi*eps0)*hbar/(c*e)
p("Kaluza-Klein radius sqrt(16 pi G) hbar/(e c) [Gaussian form]", R_KK, "m")
p("  2 pi R (Klein 1926: 0.8e-30 cm)", 2*pi*R_KK*100, "cm"); close(R_KK, sqrt(16*pi/alpha)*sqrt(hbar*G/c**3), 1e-9)
p("  bound |Rdot/R| if alpha ~ 1/R^2, from adot/a = 1.8(2.5)e-19/yr (2 sigma)", (1.8+2*2.5)*1e-19/2, "1/yr")
p("WGC: electron e/(m sqrt(4 pi eps0 G))", e/(me*k_g)); p("WGC: proton", e/(mp*k_g))
p("Aharonov-Bohm: h/e", h/e, "Wb"); p("  h/2e (Tonomura 1986 half-quantum)", h/(2*e), "Wb")
print("\nall asserts passed")
