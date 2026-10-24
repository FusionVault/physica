"""Reproduce published numbers for superconductor-gravity claims and compare with GR. SI units."""
from math import pi, sqrt, cos, atan, log10, tanh
G = 6.67430e-11; c = 2.99792458e8; hbar = 1.054571817e-34; h = 2*pi*hbar
me = 9.1093837015e-31; e = 1.602176634e-19; mu0 = 1.25663706212e-6; eps0 = 8.8541878128e-12
g0 = 9.80665

def p(name, val, unit=""):
    print(f"{name:62s} {val:.4g} {unit}")

print("== 1. Podkletnov 1997 disc (cond-mat/9701074) ==")
Ro, Ri, t, rho_ybco = 0.1375, 0.040, 0.010, 6380.0  # rho: YBCO X-ray density ~6.38 g/cm3 (assumed)
V = pi*(Ro**2 - Ri**2)*t; M = rho_ybco*V
w = 5000*2*pi/60
I = 0.5*M*(Ro**2 + Ri**2)
p("disc volume", V, "m^3"); p("disc mass (assumed density)", M, "kg")
p("omega at 5000 rpm", w, "rad/s"); p("rim speed", w*Ro, "m/s"); p("rim centripetal accel / g", w**2*Ro/g0)
p("rotational KE at 5000 rpm", 0.5*I*w**2, "J"); p("KE/c^2 (mass-energy)", 0.5*I*w**2/c**2, "kg")
sigma = rho_ybco*t
for z in (0.025, 0.5, 3.0):  # on-axis Newtonian field of a uniform annulus (adds to weight)
    gz = 2*pi*G*sigma*((z/sqrt(z**2+Ri**2)) - (z/sqrt(z**2+Ro**2)))
    p(f"GR/Newton: annulus pull at z={z} m on axis (dg/g, sign=+weight)", gz/g0)
# sample sits above the annulus, not the hole: use point above mid-annulus via thin-disc max 2 pi G sigma
p("upper bound 2*pi*G*sigma/g (infinite sheet)", 2*pi*G*sigma/g0)
p("claimed max weight change (decrease)", -0.021)
p("|claim| / GR upper bound", 0.021/(2*pi*G*sigma/g0))
# Modanese-type shadow geometry, uniform Earth: fraction of g from cone of half-angle th = 1 - cos^3(th)
for hh in (1, 2, 3, 4, 6, 8, 10, 3.0/Ro):
    th = atan(1/hh)
    p(f"shadow fraction, uniform Earth, h={hh:.3g} disc radii", 1 - cos(th)**3)

print("\n== 2. Impulse gravity generator (physics/0108005) ==")
l, m_s = 0.800, 0.0185
tab1 = [(500,56.5,2.0,3.6),(750,91.3,5.2,9.5),(1000,110.4,7.7,13.9),(1250,123.0,9.5,17.3),(1500,131.6,10.9,19.8),(1750,137.6,11.9,21.7),(2000,142.0,12.7,23.1)]
tab2 = [(500,40.0,1.0,1.8),(750,70.9,3.2,5.7),(1000,85.3,4.6,8.3),(1250,94.6,5.6,10.2),(1500,100.8,6.4,11.6),(1750,104.7,6.9,12.5),(2000,107.1,7.2,13.1)]
for name, tab in (("Emitter 2 (Table 1)", tab1), ("Emitter 1 (Table 2)", tab2)):
    print(name)
    for V_kV, dl, dh_pub, dE_pub in tab:
        d = dl/1000
        dh = l - sqrt(l**2 - d**2)
        dE = m_s*g0*dh
        v = sqrt(2*g0*dh)  # = F*t, eq. (1)
        print(f"  {V_kV:5d} kV  dh={dh*1e3:6.2f} mm (pub {dh_pub})  dE={dE*1e4:6.2f}e-4 J (pub {dE_pub})  v=Ft={v:.3f} m/s  F(t=1e-4s)/g={v/1e-4/g0:7.0f}")
v = sqrt(2*g0*(l - sqrt(l**2 - 0.142**2)))
p("momentum of 18.5 g pendulum at 2 MV", m_s*v, "kg m/s"); p("E/c for same energy", m_s*g0*(l-sqrt(l**2-0.142**2))/c, "kg m/s")
p("p / (E/c)", m_s*v/(m_s*g0*(l-sqrt(l**2-0.142**2))/c))
Cst, n = 25e-9, 20
p("Marx bank stored energy at 100 kV/stage", 0.5*n*Cst*(100e3)**2, "J")
p("Marx erected capacitance", Cst/n, "F"); p("Marx charge at 2 MV", Cst/n*2e6, "C")
p("paper's Q = I*t (1e4 A x 1e-5 s)", 1e4*1e-5, "C"); p("paper's U=QV/2 with Q=0.1 C, V=2 MV", 0.1*2e6/2, "J")
dm = 2500/c**2
p("GR: Newtonian pull of 2.5 kJ mass-energy at 6 m", G*dm/36, "m/s^2")
p("GR: delta-v over 1e-4 s", G*dm/36*1e-4, "m/s"); p("claimed v / GR delta-v", v/(G*dm/36*1e-4))
p("Newtonian pull of a 1000 kg apparatus at 6 m", G*1000/36, "m/s^2"); p("claimed F (t=1e-4 s) / that", (v/1e-4)/(G*1000/36))
p("64c claim: delay over 1211 m", 1211/(64*c), "s"); p("light time over 1211 m", 1211/c, "s")

print("\n== 3. Modanese numbers (hep-th/9505094; physics/0108005) ==")
n_pairs = 1e26  # 1e20 cm^-3
p("pair rest-energy density m n c^2 (n=1e20 cm^-3)", 2*me*c**2*n_pairs, "J/m^3")
rho_L = 5.9e-27  # kg/m^3, observed dark-energy density (Planck 2018 approx)
p("observed vacuum energy density rho_L c^2", rho_L*c**2, "J/m^3")
p("ratio pair energy density / vacuum energy", 2*me*c**2*n_pairs/(rho_L*c**2))
# natural units check of 1996 numbers
mpair_cm = 2*me*c/hbar/100  # cm^-1
p("Cooper-pair mass in natural units", mpair_cm, "cm^-1 (paper: ~1e10)")
phi0 = sqrt(1e20/mpair_cm); p("phi0 = sqrt(n/m)", phi0, "cm^-1 (paper: ~1e5)")
p("m^2 phi0^2 = m n", mpair_cm*1e20, "cm^-4 (paper: ~1e30)")
Lam_obs = 1.1e-52  # m^-2
p("observed Lambda in cm^-2", Lam_obs*1e-4, "cm^-2 (paper requires ~1e-36 at 1e-8..1e-6 cm)")
p("orders by which required local Lambda exceeds observed", log10(1e-36/(Lam_obs*1e-4)))
# 2001 GL estimates: L2 ~ hbar^2/(2m) n/xi^2 (amplitude-gradient term), L1 = beta n^2 / 2
mstar = 2*me; ns = 1e27; xi = 1e-9; kappa = 1e2; estar = 2*e
L2 = hbar**2/(2*mstar)*ns/xi**2
beta = 2*mu0*hbar**2*estar**2*kappa**2/mstar**2  # from kappa^2 = m^2 beta/(2 mu0 hbar^2 e^2)
p("L2 ~ hbar^2 n/(2 m xi^2)", L2, "J/m^3 (paper: ~1e7)")
p("L1 = beta n^2/2 with paper's kappa relation (e*=2e)", 0.5*beta*ns**2, "J/m^3 (paper: ~1e4)")
p("orders by which L2 exceeds observed vacuum energy", log10(L2/(rho_L*c**2)))

print("\n== 4. Tajmar / de Matos ==")
coef = me**2/(mu0*e**2)
p("m_e^2/(mu0 e^2) = rho* lambda_L^2", coef, "kg/m")
p("... / g0 (Tajmar eq.6 coefficient 2.624e-18)", coef/g0)
rhoNb, rhoPb = 8570.0, 11340.0
for lam in (27.6e-9, 39e-9, 41.6e-9):
    p(f"Nb rho*/rho with lambda_L={lam*1e9:.1f} nm", coef/lam**2/rhoNb)
p("Tate dm/m = (1.000084-0.999992)/0.999992", (1.000084-0.999992)/0.999992)
p("Tate significance (92 ppm / 21 ppm)", 92/21)
p("B_g/omega needed = 2 dm/m (paper 1.84e-4)", 2*9.2e-5)
r_s = 0.036
p("pred. coupling g/omegadot, r=3.6cm, rho*/rho=1.74e-6", r_s*1.74e-6/g0, "s^2 (paper -0.7e-8)")
p("pred. coupling from Tate dm/m, r=3.6 cm", r_s*9.2e-5/g0, "s^2 (paper -37e-8)")
p("measured: 100 ug / 1500 rad/s^2", 100e-6/1500, "s^2 (paper -6.6e-8)")
k_meas = 9.64e-8*g0/r_s
p("measured k = Omega_eff/omega from -9.64e-8 s^2 at r=3.6 cm", k_meas, "(paper ~2.6e-5)")
Ro, Ri, H = 0.075, 0.069, 0.015
Mring = rhoNb*pi*(Ro**2-Ri**2)*H
p("Nb ring mass 150x138x15 mm", Mring, "kg")
k_GR = 4*G*Mring/(c**2*(Ro+Ri))  # Tajmar 0707.3806 eq.(4), = 2GM/(c^2 R)
p("GR k = Omega_LT/omega at ring centre = 2GM/(c^2 Rbar)", k_GR)
p("ratio k_meas/k_GR", k_meas/k_GR); p("log10", log10(k_meas/k_GR))
p("ratio 2006 B_g convention (2k)/k_GR", 2*k_meas/k_GR)
a_GR = r_s*k_GR*1500
p("GR tangential accel at r=3.6cm, omegadot=1500", a_GR, "m/s^2"); p("claimed 100 ug / GR", 100e-6*g0/a_GR)
p("Tajmar 2007 gyro: gamma for Bg/omega=3.2e-8 with this ring", 3.2e-8/k_GR, "(paper 1.2e18 for Nb+Al)")
p("2011 limit 4e-11 / k_GR", 4e-11/k_GR); p("2011 limit vs de Matos prediction 3.95e-6", 3.95e-6/4e-11)
p("2011 limit vs 2007 gyro claim 3.2e-8", 3.2e-8/4e-11)
# equivalent GR source density for the claim
MR = k_meas*c**2/(2*G); p("GR M/R needed for k_meas", MR, "kg/m")
p("equivalent density in the ring volume", MR*0.072/(pi*(Ro**2-Ri**2)*H), "kg/m^3")
p("equivalent density for gamma=1e18", 1e18*rhoNb, "kg/m^3")
J_E, R_E = 5.86e33, 6.371e6
p("Earth Omega_LT at equator GJ/(c^2 R^3)", G*J_E/(c**2*R_E**3), "rad/s")

print("\n== 5. London moment, Tate, DeWitt, gravitational penetration depth ==")
p("London moment B/omega = 2 m_e/e", 2*me/e, "T s/rad")
p("B at 650 rad/s", 2*me/e*650, "T")
p("GR frame-drag correction to B at 650 rad/s for Nb ring", 2*me/e*k_GR*650, "T")
p("fractional GR correction to measured Cooper-pair mass", k_GR)
p("Hoang 2020 Nb gap significance 18/2.1", 18/2.1)
M, v_over_d = 1.0, 300*2*pi/2
I_dw = 8*pi*eps0*G*me*M*v_over_d/e
p("DeWitt current I ~ 8 pi eps0 G m M v/(e d), M=1kg, v/d=300*pi", I_dw, "A (K&W: ~1e-29)")
lamNb = 39e-9; ns_Nb = 2*me/(mu0*(2*e)**2*lamNb**2)  # pair density
lam_g = sqrt(c**2/(4*pi*G*(2*me)*ns_Nb))
p("Nb Cooper-pair density from lambda_L=39 nm", ns_Nb, "m^-3")
p("Cooper-pair mass density rho*", 2*me*ns_Nb, "kg/m^3")
p("gravitational penetration depth sqrt(c^2/(4 pi G rho*))", lam_g, "m")
p("lambda_g/lambda_e", lam_g/lamNb, "(Ummarino: ~1e21)")
p("lambda_g for neutron-star density 1e17 kg/m^3", sqrt(c**2/(4*pi*G*1e17)), "m (Lano/Modanese: ~12 km)")

print("\n== 6. Li & Torr PRD 43, 457 ==")
mu_g = 4*pi*G/c**2
p("mu_g = 4 pi G/c^2", mu_g, "m/kg")
p("induced B_g per tesla: mu_g (m/q) B0 / mu0", mu_g*(me/e)/mu0, "s^-1 per T")
p("residual B from Earth B_g~1e-14 s^-1: (m/q) B_g", me/e*1e-14, "T")

print("\n== 7. Ummarino & Gallerati ==")
def gstar(lam, kap, Hc, D): return lam*kap*me*g0/(sqrt(2)*e*Hc*D)
p("YBCO 77 K g*", gstar(3.3e-7, 94.4, 0.2, 3.2e-4), "(paper 2e-11)")
p("Pb 6.26 K g*", gstar(7.8e-8, 0.48, 0.018, 1.0), "(paper 8.2e-17)")
p("YBCO tau = lambda^2/D", (3.3e-7)**2/3.2e-4, "s (paper 3.4e-10)")
p("YBCO eta = 4 pi mu0 sigma D", 4*pi*mu0*(1/4e-7)*3.2e-4, "(paper 1.27e-2)")
p("Pb eta = 4 pi mu0 sigma D", 4*pi*mu0*(1/2.5e-9)*1.0, "(paper 6.6e3)")
p("YBCO recovery time eta*tau (eq.93)", 1.27e-2*3.4e-10, "s")
lam_b = 500e-9/sqrt((109-105)/109)
p("Bi2212 105 K tau", lam_b**2/1e-3, "s (paper 6.8e-9)"); p("Bi2212 eta = mu0 sigma D", mu0*(1/3.6e-6)*1e-3, "(paper 3.5e-4)")
for ell in (1e-9, 300e-9):
    om = 2*me*g0*ell/hbar
    p(f"Josephson omega=2 m g l/hbar, l={ell*1e9:.0f} nm", om, "s^-1"); p("   period", 2*pi/om, "s"); p("   V_g = m g l/e", me*g0*ell/e, "V")
eta_b = mu0*(1/3.6e-6)*1e-3
p("Bi2223 Ez = 4 pi Bc D/(eta L), L=0.15 m", 4*pi*0.32*1e-3/(eta_b*0.15), "V/m (paper ~77)")
p("Schiff-Barnhill field m g/e (mainstream equilibrium)", me*g0/e, "V/m")
