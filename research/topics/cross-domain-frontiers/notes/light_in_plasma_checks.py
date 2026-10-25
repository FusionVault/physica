"""Recompute the published numbers behind light_in_plasma_and_magnetic_fields.md.

Standard library only. Run: python3 -I light_in_plasma_checks.py
Each check prints the recomputed value beside the published one and asserts a tolerance.
"""
from math import pi, sqrt, log, asinh, degrees

# CODATA 2018/2022 (SI)
E = 1.602176634e-19      # C
ME = 9.1093837139e-31    # kg
MP = 1.67262192595e-27   # kg
C = 299792458.0          # m/s
EPS0 = 8.8541878188e-12  # F/m
HBAR = 1.054571817e-34   # J s
H = 2 * pi * HBAR
RE = E**2 / (4 * pi * EPS0 * ME * C**2)  # classical electron radius
MEC2_MEV = ME * C**2 / E / 1e6
MPC2_MEV = MP * C**2 / E / 1e6
PC = 3.0856775814913673e16  # m
EV = E


def check(name, got, want, rtol):
    ok = abs(got - want) <= rtol * abs(want)
    print(f"{'ok ' if ok else 'BAD'} {name}: {got:.5g} (published {want:.5g}, tol {rtol:.0%})")
    assert ok, name


def n_c(lam):  # critical density (m^-3) for wavelength lam (m)
    w = 2 * pi * C / lam
    return EPS0 * ME * w * w / E**2


def omega_p(n):  # electron plasma frequency (rad/s), n in m^-3
    return sqrt(n * E**2 / (EPS0 * ME))


# 1. Relativistic self-focusing: P_c = 8 pi eps0 me^2 c^5 / e^2 * (w/wp)^2 (Gaussian 2c(mc^2/e)^2)
PC_UNIT = 8 * pi * EPS0 * ME**2 * C**5 / E**2
check("P_c coefficient (GW), Sprangle 1987 / Esarey 2009", PC_UNIT / 1e9, 17.4, 0.005)
check("n_c at 1 um (cm^-3)", n_c(1e-6) / 1e6, 1.115e21, 0.002)

# 2. Lu et al. 2007 PRSTAB 10, 061301, Eqs. 1-7, 10
def lu(P, n_cm3, lam_um, tau_fwhm, pc_gw=17.0):
    n = n_cm3 * 1e6
    lam = lam_um * 1e-6
    r = n_c(lam) / n                          # (w0/wp)^2
    Pc = pc_gw * 1e9 * r
    a0 = 2 * (P / Pc) ** (1 / 3)              # Eq. 2
    kp = omega_p(n) / C
    R = 2 * sqrt(a0) / kp                     # Eq. 1
    Ld = 2 / 3 * r * R                        # Eq. 4
    Letch = r * C * tau_fwhm                  # Eq. 3
    dE = 2 / 3 * MEC2_MEV * r * a0 / 1e3      # Eq. 6 first line (GeV)
    dE_pr = 1.7 * (P / 100e12) ** (1 / 3) * (1e18 / n_cm3) ** (2 / 3) * (0.8 / lam_um) ** (4 / 3)  # Eq. 6
    N = 8 / 15 / (2 * pi / lam * RE) * sqrt(P / (PC_UNIT / 2))  # Eq. 10, m^2c^5/e^2 -> 4pi eps0 me^2 c^5/e^2
    return dict(a0=a0, R=R, Ld=Ld, Letch=Letch, dE=dE, dE_pr=dE_pr, N=N, Pc=Pc)

s = lu(200e12, 1.5e18, 0.8, 30e-15)
check("Lu sample 200 TW, 1.5e18: dE Eq.6 (GeV)", s["dE_pr"], 1.6, 0.03)
check("Lu sample: dephasing L_d (cm)", s["Ld"] * 100, 1.31, 0.08)
check("Lu sample: pump depletion L_pd (cm)", s["Letch"] * 100, 0.96, 0.10)
check("Lu Eq.10: N at 100 TW, 0.8 um", lu(100e12, 1e18, 0.8, 30e-15)["N"], 2.5e9, 0.05)
# Eq. 6 two forms agree (coefficient 1.7)
q = lu(100e12, 1e18, 0.8, 30e-15)
check("Lu Eq.6: (2/3)mc^2 (w0/wp)^2 a0 vs 1.7 GeV practical form", q["dE"], q["dE_pr"], 0.06)

# Picksley et al. 2024 PRL 133, 255001: 500 TW, 1e17 cm^-3, 815 nm, 40 fs, w0 = 53 um, a0 ~ 2.2, 9.2 GeV
p = lu(500e12, 1e17, 0.815, 40e-15)
check("Picksley: Lu matched spot w0 (um) vs 53 um used", p["R"] * 1e6, 53.0, 0.05)
check("Picksley: Lu a0 vs ~2.2 reported", p["a0"], 2.2, 0.10)
I_w0 = 2 * 500e12 / (pi * (53e-6) ** 2)          # W/m^2, Gaussian peak
a0_I = 0.855 * 0.815 * sqrt(I_w0 / 1e22)          # a0 = 0.855 lam[um] sqrt(I/1e18 W cm^-2)
check("Picksley: a0 from 500 TW into w0 = 53 um", a0_I, 2.2, 0.10)
print(f"    Picksley: Lu Eq.6 energy {p['dE_pr']:.1f} GeV (measured 9.2; their simulated optimum >13)")
print(f"    Picksley: L_d {p['Ld']*100:.0f} cm, L_pd {p['Letch']*100:.0f} cm (channel 30 cm)")
# Gonsalves et al. 2019 PRL 122, 084801: 850 TW, ~3.0e17 cm^-3, 7.8 GeV
g = lu(850e12, 3.0e17, 0.815, 35e-15)
check("Gonsalves 2019: Lu Eq.6 (GeV) vs 7.8 measured", g["dE_pr"], 7.8, 0.06)

# Lu's field falls linearly to zero over L_d, so a stage of length L < L_d gains dE * (2x - x^2), x = L/L_d
def truncated(d, L):
    x = min(L / d["Ld"], 1.0)
    return d["dE_pr"] * (2 * x - x * x)

print(f"    Picksley truncated: to L_pd = {truncated(p, p['Letch']):.1f} GeV, to 30 cm = {truncated(p, 0.30):.1f} GeV (measured 9.2)")
print(f"    Gonsalves: L_d {g['Ld']*100:.0f} cm vs 20 cm capillary; truncated gain {truncated(g, 0.20):.1f} GeV (measured 7.8)")

# 3. Plasma-channel matched depth dn_c = 1/(pi r_e w0^2) (Lu 2007; Esarey 2009)
print(f"    channel depth for w_m = 37 um: {1/(pi*RE*(37e-6)**2)/1e6:.3g} cm^-3")
check("dn_c/n_p = 4/(kp w0)^2 identity", 1 / (pi * RE * (53e-6) ** 2) / 1e23,
      4 / (omega_p(1e23) / C * 53e-6) ** 2, 1e-9)

# 4. Ponderomotive energies
w1 = 2 * pi * C / 1e-6
E0 = sqrt(2 * 1e18 / (EPS0 * C))                 # 1e14 W/cm^2 linear, peak field
Up = E**2 * E0**2 / (4 * ME * w1**2) / EV
check("Up = e^2E^2/(4 m w^2) at 1e14 W/cm^2, 1 um (eV)", Up, 9.33, 0.005)
a0_137 = sqrt(1.37e18 / 1.37e18)
check("Wilks T_h at I lam^2 = 1.37e18 (a0 = 1, MeV)", MEC2_MEV * (sqrt(1 + a0_137**2 / 2) - 1), 0.115, 0.01)
check("a0 coefficient 0.855 (I = eps0 c E0^2/2)",
      E * sqrt(2 * 1e22 / (EPS0 * C)) / (ME * C * w1), 0.855, 0.002)

# 5. Hole boring (Robinson 2009; Macchi 2013 RMP Eqs. 28-30)
def hb(Pi):
    b = sqrt(Pi) / (1 + sqrt(Pi))
    E_reflect = 2 * MPC2_MEV * Pi / (1 + 2 * sqrt(Pi))
    return b, E_reflect

for Pi in (1e-4, 0.1, 1.0):
    b, Er = hb(Pi)
    beta_i = 2 * b / (1 + b * b)                 # ion bounced at 2 v_hb (relativistic addition)
    check(f"HB Eq.30 = reflected-ion kinetic energy, Pi={Pi}", Er, MPC2_MEV * (1 / sqrt(1 - beta_i**2) - 1), 1e-9)
    P_rad = 2 * (1 - b) / (1 + b)               # Eq. 28 LHS / (I/c), with n_i m_i c^3 = I/Pi
    g_hb = 1 / sqrt(1 - b * b)
    check(f"HB Eq.28 momentum balance, Pi={Pi}", P_rad, (1 / Pi) * g_hb * 2 * g_hb * b * b, 1e-9)
# Palmer et al. 2011 PRL 106 014801 (via Macchi 2013): 10 um, 6e15 W/cm^2 CP, a0 ~ 0.5, n_e ~ few n_c
I_cm2 = 6e15
a0_pal = sqrt(I_cm2 * 1e4 / (ME * C**3 * n_c(10e-6)))   # I = a0^2 me c^3 n_c (circular)
check("Palmer 2011 a0 from 6e15 W/cm^2 at 10 um", a0_pal, 0.5, 0.08)
for k in (1, 2, 4):
    Pi = a0_pal**2 * (ME / MP) / k
    print(f"    Palmer HB at n_e = {k} n_c: E_max {hb(Pi)[1]:.3f} MeV (observed up to 1.2 MeV)")

# 6. Light sail (Macchi 2013 Eqs. 33-36): F = 2 int I dt / (sigma c^2)
F = 2 * 1e12 / (1000 * 10e-9 * C**2)               # 1e8 J/cm^2, 1 g/cm^3, 10 nm
check("Light sail practical F = 2.2 F_1e8/(rho_1 l_10)", F, 2.2, 0.02)
beta = ((1 + F) ** 2 - 1) / ((1 + F) ** 2 + 1)
check("Light sail Eq.35 = mp c^2 (gamma - 1)", MPC2_MEV * F**2 / (2 * (F + 1)), MPC2_MEV * (1 / sqrt(1 - beta**2) - 1), 1e-9)
check("Light sail Eq.36 eta = 2b/(1+b) = 1 - 1/(F+1)^2", 2 * beta / (1 + beta), 1 - 1 / (F + 1) ** 2, 1e-12)

# 7. TNSA, Mora 2003 / Macchi Eq. 24: u_f = 2 c_s ln(tau + sqrt(tau^2+1)) => E = 2 Z T_h ln^2(...)
Th, Z, A, tau = 1.0, 1, 1, 10.0
cs2 = Z * Th / (A * MPC2_MEV)                       # c_s^2 in units of c^2 (non-relativistic ion)
uf = 2 * sqrt(cs2) * asinh(tau)
check("Mora: (1/2) m u_f^2 = 2 Z T_h asinh^2(tau)", 0.5 * A * MPC2_MEV * uf**2, 2 * Z * Th * asinh(tau) ** 2, 1e-12)

# 8. Faraday rotation / rotation measure
K_RM = E**3 / (8 * pi**2 * EPS0 * ME**2 * C**3)     # rad m^-2 per (m^-3 T m)
check("RM SI coefficient e^3/(8 pi^2 eps0 me^2 c^3)", K_RM, 2.63e-13, 0.003)
k_astro = K_RM * 1e6 * 1e-10 * PC                   # cm^-3, uG, pc
check("RM astro coefficient 0.812 rad m^-2/(cm^-3 uG pc)", k_astro, 0.812, 0.002)
check("<B_par> = 1.232 RM/DM uG", 1 / k_astro, 1.232, 0.002)
K_CM = E**4 / (16 * pi**3 * EPS0 * ME**3 * C**4)    # Cotton-Mouton plasma phase coefficient
check("Plasma Cotton-Mouton coefficient 2.45e-11 rad T^-2 m^-2", K_CM, 2.45e-11, 0.005)
# ATNF catalogue v2.6.1 (DM pc cm^-3, RM rad m^-2)
for name, dm, rm in (("B0329+54", 26.7641, -64.33), ("B0531+21", 56.77, -45.44),
                     ("B0833-45", 67.771, 31.38), ("B1937+21", 71.0, 11.0)):
    print(f"    {name}: <B_par> = {rm/(k_astro*dm):+.3f} uG")
# FRB 20121102A (Michilli et al. 2018 Nature 553, 182): RM_src = 1.46e5, 1.33e5; z = 0.193
z = 0.193
print(f"    FRB 20121102A observed-frame RM = {1.46e5/(1+z)**2:.4g} and {1.33e5/(1+z)**2:.4g} rad m^-2")
check("FRB 20121102A 2016 -> 2023 (30,755) drop of ~70% (arXiv 2507.15790)", 1 - 30755 / (1.46e5 / (1 + z) ** 2), 0.70, 0.02)
# Faraday rotation angle of a 1 m (21 cm) wave
print(f"    rotation at 21 cm through RM=53.8 (CHIME median |RM|): {degrees(53.8*0.21**2):.0f} deg")

# 9. Dispersion delay and the plasma photon mass
K_DM = E**2 / (8 * pi**2 * EPS0 * ME * C)           # s Hz^2 per m^-2
check("k_DM = 4.1488 ms GHz^2 / (pc cm^-3)", K_DM * PC * 1e6 / 1e18 * 1e3, 4.148808, 1e-4)
ne, d, f = 0.03e6, 1000 * PC, 1.4e9
fp = omega_p(ne) / (2 * pi)
m_eff_ev = H * fp / EV                               # m_eff c^2 = h f_p = hbar w_p
t_plasma = K_DM * ne * d / f**2
t_mg = d / C * 0.5 * (m_eff_ev * EV / (H * f)) ** 2  # gravitational::propagation::massive_graviton_delay
check("plasma DM delay == massive-graviton delay with m c^2 = hbar w_p", t_mg, t_plasma, 1e-12)
check("f_p at 0.03 cm^-3 (kHz) vs gravitational::hf doc 1.56 kHz", fp / 1e3, 1.56, 0.01)
print(f"    photon effective mass at 0.03 cm^-3: {m_eff_ev:.3g} eV = {m_eff_ev*EV/C**2:.3g} kg (FRB Proca bound 4.8e-51 kg)")

# 10. Magneto-optics in matter
print(f"    TGG (-40 rad/T/m at 1064 nm): B L for 45 deg = {pi/4/40:.4f} T m")
V_he, kcm_he = 3.87e-5, 2.19e-16                     # Cadene et al. 2013 PRA 88 043815 (per atm)
A_e = 2 * (E**2 / (4 * pi * EPS0 * HBAR * C)) ** 2 * HBAR**3 / (45 * 4e-7 * pi * ME**4 * C**5)
dn_vac = 3 * A_e                                     # gravitics::strong_field::birefringence at 1 T
check("QED vacuum Delta n_u = 3 A_e (T^-2)", dn_vac, 3.97e-24, 0.01)
print(f"    He pressure whose Cotton-Mouton equals QED vacuum: {dn_vac/kcm_he:.2g} atm = {dn_vac/kcm_he*101325:.2g} Pa")

# 11. Magnetised plasma cut-offs (Stix), ITER from plasmoid::tokamak::ITER (B = 5.3 T)
B, ne = 5.3, 1e20
fce = E * B / (2 * pi * ME)
fpe = omega_p(ne) / (2 * pi)
fR = 0.5 * (fce + sqrt(fce**2 + 4 * fpe**2))
fUH = sqrt(fce**2 + fpe**2)
check("f_ce per tesla (GHz)", E / (2 * pi * ME) / 1e9, 27.99, 0.001)
print(f"    ITER 5.3 T, 1e20: f_ce {fce/1e9:.1f} GHz, f_pe {fpe/1e9:.1f} GHz, R-cutoff {fR/1e9:.1f} GHz, UH {fUH/1e9:.1f} GHz")
lam = 118.8e-6
print(f"    ITER-like Cotton-Mouton phase, 118.8 um, 4 m chord: {degrees(K_CM*lam**3*ne*B**2*4):.0f} deg")

# 12. Flux compression of the NIF magnetised implosion (Moody et al. 2022 PRL 129 195002)
print(f"    26 T -> 4.9 kT needs radial convergence {sqrt(4900/26):.1f} under B r^2 = const")

# 13. Relativistic transparency threshold n_e < gamma n_c, gamma = sqrt(1 + a0^2/2)
for a0 in (1, 10, 30):
    print(f"    a0 = {a0}: relativistic critical density {sqrt(1+a0*a0/2):.2f} n_c")
print("all checks passed")
