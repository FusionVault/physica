"""Checks for fringe Coulomb-barrier / Schwinger claims. Run: python3 -I barrier_checks.py
Constants CODATA 2018. Inputs marked RECALLED were not re-fetched this session."""
import math

e = 1.602176634e-19; hbar = 1.054571817e-34; c = 2.99792458e8; eps0 = 8.8541878128e-12
mu0 = 1.25663706212e-6; me = 9.1093837015e-31; alpha = 7.2973525693e-3; kB = 1.380649e-23
MeV = 1e6 * e; NA = 6.02214076e23
mec2_keV = 510.99895; md_c2_keV = 1875612.9  # deuteron rest energy
mp_c2 = 938.27209; mn_c2 = 939.56542  # MeV

def hdr(s): print(f"\n== {s} ==")

# 1. D-D neutrons per watt
hdr("1 D-D yield per watt")
Q_n, Q_p, Q_g = 3.269, 4.033, 23.847  # MeV, standard Q values
Qbar = 0.5 * (Q_n + Q_p)
R_per_W = 1 / (Qbar * MeV)
n_per_W = 0.5 * R_per_W
print(f"mean Q = {Qbar:.3f} MeV; reactions/s/W = {R_per_W:.3e}; neutrons/s/W = {n_per_W:.3e}")
He_per_W = 1 / (Q_g * MeV)
print(f"if all heat from d+d->4He: 4He/s/W = {He_per_W:.3e} (and the same number of 23.8 MeV gammas if the standard channel)")
# Fleischmann-Pons Table 1 + p.306: 0.4x10 cm rod at 64 mA/cm2: 1.751 W excess, neutrons ~4e4 /s
P_fp, Nn_fp = 1.751, 4e4
print(f"F-P 0.4x10 cm @64 mA/cm2: expected n/s = {P_fp*n_per_W:.3e}; claimed 4e4; deficit = {P_fp*n_per_W/Nn_fp:.2e}")
print(f"F-P 0.4x10 cm @512 mA/cm2 (26.8 W): expected n/s = {26.8*n_per_W:.3e}")
# dose at 1 m; RECALLED fluence-to-ambient-dose coefficient ~4e-10 Sv cm2 for ~2.5 MeV neutrons (ICRP 74 order)
phi = P_fp * n_per_W / (4 * math.pi * 100.0**2)
print(f"flux at 1 m = {phi:.3e} n/cm2/s; dose rate ~ {phi*4e-10*3600:.1f} Sv/h (RECALLED coefficient)")
# energy per Pd atom for 4 MJ/cm3 (F-P p.305); Pd 12.02 g/cm3, 106.42 g/mol
nPd = 12.02 / 106.42 * NA
print(f"Pd atoms/cm3 = {nPd:.3e}; 4 MJ/cm3 = {4e6/nPd/e:.0f} eV per Pd atom")
V_rod = math.pi * 0.2**2 * 10
rate_pair_heat = P_fp * R_per_W / (nPd * V_rod)
rate_pair_n = 2 * Nn_fp / (nPd * V_rod)
print(f"rod volume {V_rod:.3f} cm3; per-pair rate needed for heat = {rate_pair_heat:.2e} /s; for the neutrons = {rate_pair_n:.2e} /s (F-P quote 1e-19)")
KN = 3e-64
print(f"vs Koonin-Nauenberg D2 3e-64/s: heat needs x{rate_pair_heat/KN:.1e} ({math.log10(rate_pair_heat/KN):.0f} orders)")

# 2. Gamow factor and screening
hdr("2 Gamow / screening")
mu_keV = md_c2_keV / 2
EG = 2 * mu_keV * (math.pi * alpha) ** 2  # keV, Z1=Z2=1
print(f"d+d Gamow energy E_G = {EG:.1f} keV; 2*pi*eta = {math.sqrt(EG):.2f}/sqrt(E[keV])")

def enh(E, U):  # sigma(E+U)/sigma(E) with S const and the 1/E geometric factor kept at E
    return math.exp(math.sqrt(EG / E) - math.sqrt(EG / (E + U)))

for E in (1.0, 2.5, 5.0, 10.0):
    print(f"E_cm={E:5.1f} keV: f(U=0.025)={enh(E,0.025):.3f} f(0.3)={enh(E,0.3):.3f} f(0.6)={enh(E,0.6):.3g}")
S0 = 53.7 + 55.6  # keV b, RECALLED Bosch-Hale S(0) for the n and p branches
kT_keV = kB * 300 / e / 1e3
def log10_sigma_thermal(U):  # barn, pair energy kT(300 K), model sigma = S/E exp(-sqrt(EG/(E+U)))
    return math.log10(S0 / kT_keV) - math.sqrt(EG / (kT_keV + U)) / math.log(10)
for U in (0.0, 0.027, 0.1, 0.3, 0.6, 1.0):
    print(f"300 K, U_e={U*1e3:6.0f} eV: sigma ~ 10^{log10_sigma_thermal(U):.1f} b")
# rate per gram for D/Pd=1 at 300 K
v_th = math.sqrt(2 * kT_keV * 1e3 * e / (md_c2_keV * 1e3 * e / c**2 / 2))  # relative speed for pair energy kT
nD_m3 = nPd * 1e6; N_g = nPd / 12.02
def rate_per_g(sig_b): return 0.5 * nD_m3 * sig_b * 1e-28 * v_th * N_g  # (1/2) n^2 sigma v per m3, divided by density
print(f"v_rel(300 K) = {v_th:.3e} m/s")
r13 = rate_per_g(1e-13)
print(f"MARRS check: sigma=1e-13 b -> {r13:.2e} reactions/s/g = {r13/R_per_W:.2f} W/g (slide: >=1 W/g)")
def U_needed(target_rate_g):
    sig = target_rate_g / rate_per_g(1.0)
    x = math.log(S0 / kT_keV / sig)  # sqrt(EG/(E+U))
    return EG / x**2 - kT_keV
fp_heat_g = P_fp * R_per_W / (V_rod * 12.02)
fp_n_g = 2 * Nn_fp / (V_rod * 12.02)
for label, r in (("F-P heat 1.751 W", fp_heat_g), ("F-P neutrons 4e4/s", fp_n_g), ("1 W/g", R_per_W),
                 ("MARRS SOA 0.1/s/g", 0.1), ("MARRS Ph2 1e6/s/g", 1e6)):
    print(f"{label:20s}: {r:.2e} reactions/s/g -> U_e needed ~ {U_needed(r)*1e3:.0f} eV (thermal model)")
# Karahadian >1e18 enhancement at sub-keV: what U_e?
for E in (0.05, 0.1, 0.2):
    lo, hi = 0.0, 1e3
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if enh(E, mid) < 1e18: lo = mid
        else: hi = mid
    print(f"enhancement 1e18 at E_cm={E} keV needs U_e ~ {lo*1e3:.0f} eV")
# Debye check (slides: R_D = 69 sqrt(T/(n_eff rho_a)) m)
R_D = math.sqrt(eps0 * kB * 300 / (e**2 * 1 * 6e28))
print(f"Debye R_D(300 K, 6e28) = {R_D:.2e} m; U_D = {e/(4*math.pi*eps0*R_D):.0f} eV")

# 3. DARPA MARRS metrics in W/g
hdr("3 MARRS metrics")
for lab, r in (("SOA", 0.1), ("Phase 1", 1e3), ("Phase 2", 1e6), ("heat-detectable", 1e9)):
    print(f"{lab:16s} {r:.0e} reactions/s/g = {r/R_per_W:.1e} W/g")
print(f"F-P heat claim per gram = {fp_heat_g:.2e} reactions/s/g = x{fp_heat_g/0.1:.1e} the MARRS SOA")

# 4. Pais plasma compression fusion device
hdr("4 Pais PCFD")
k = 1 / (4 * math.pi * eps0)
for r in (0.15, 1.0):
    print(f"1 C on sphere r={r} m: E_surf={k/r**2:.2e} V/m, V={k/r:.2e} V, U={k/(2*r):.2e} J")
R, w = 0.15, 1e4
B_shell = (2 / 3) * mu0 * 1.0 * w / (4 * math.pi * R)
print(f"spinning 1 C shell R={R} m, w={w:.0e} rad/s: B_inside = {B_shell:.2e} T; claimed 1e6 T is x{1e6/B_shell:.1e}")
Bc = 1e6
pB = Bc**2 / (2 * mu0)
Vsph = 4 / 3 * math.pi * 0.15**3
print(f"B=1e6 T: pressure {pB:.2e} Pa ({pB/101325:.1e} atm); energy in 0.3 m sphere {pB*Vsph:.2e} J = {pB*Vsph/4.184e15:.2f} Mt TNT")
print(f"gain ~ B^3: B=1e6 T -> 1e18 implies G=(B/1 T)^3; at the 30 T the patent cites: {30**3:.1e}")
T_keV_DD = 232e6 / 1.160452e7
n = 1e20
print(f"232 MC = {T_keV_DD:.1f} keV; Lawson 3e21 keV s/m3 at n=1e20 -> tau >= {3e21/(n*T_keV_DD):.2f} s")
p_pl = 2 * n * T_keV_DD * 1e3 * e
print(f"plasma pressure {p_pl:.2e} Pa -> beta=1 field {math.sqrt(2*mu0*p_pl):.2f} T")
def sigv_DD(T_keV, N=20000):  # Maxwellian <sigma v>, constant S0 (RECALLED), m^3/s
    mu_kg = md_c2_keV * 1e3 * e / c**2 / 2
    tot = 0.0; Emax = 60 * T_keV; dE = Emax / N
    for i in range(1, N + 1):
        E = i * dE
        sig = S0 / E * math.exp(-math.sqrt(EG / E)) * 1e-28  # m^2
        tot += sig * E * math.exp(-E / T_keV) * dE
    return math.sqrt(8 / (math.pi * mu_kg)) / (T_keV * 1e3 * e) ** 1.5 * tot * (1e3 * e) ** 2  # keV^2 -> J^2
sv = sigv_DD(T_keV_DD)
Pden = 0.5 * n**2 * sv * Qbar * MeV
print(f"<sigma v>_DD(20 keV, const S) = {sv:.2e} m3/s; P_fus(n=1e20) = {Pden:.2e} W/m3; in 0.3 m sphere {Pden*Vsph:.0f} W")
n_GW = math.sqrt(2 * 1e9 / Vsph / (sv * Qbar * MeV))
p_GW = 2*n_GW*T_keV_DD*1e3*e
print(f"1 GW in 0.3 m sphere at 20 keV needs n = {n_GW:.1e} m^-3, pressure {p_GW:.1e} Pa, beta=1 field {math.sqrt(2*mu0*p_GW):.0f} T")

# 5. Schwinger / vacuum polarization
hdr("5 Schwinger and Euler-Heisenberg")
ES = me**2 * c**3 / (e * hbar); BS = ES / c
IS_half = 0.5 * c * eps0 * ES**2 / 1e4; IS_full = c * eps0 * ES**2 / 1e4
print(f"E_S = {ES:.4e} V/m, B_S = {BS:.4e} T, I_S = {IS_half:.2e} (cε0E²/2) or {IS_full:.2e} (cε0E²) W/cm2")
I_rec = 1.1e23
E_rec = math.sqrt(2 * I_rec * 1e4 / (c * eps0))
print(f"record laser 1.1e23 W/cm2 -> E0 = {E_rec:.2e} V/m = {E_rec/ES:.2e} E_S; intensity short by {IS_half/I_rec:.1e}")
def schwinger_rate(E):  # pairs m^-3 s^-1, leading term, constant field
    pref = (e * E) ** 2 / (4 * math.pi**3 * hbar**2 * c)
    return math.log10(pref) - math.pi * ES / E / math.log(10)
for lab, E in (("vacuum breakdown 1e9", 1e9), ("record laser", E_rec), ("Pais sigma=5e4 C/m2", 5e4/eps0),
               ("0.1 E_S", 0.1 * ES), ("E_S", ES)):
    print(f"{lab:22s} E={E:.2e}: log10(rate m^-3 s^-1) = {schwinger_rate(E):.4g}; EH fractional correction = {alpha/(45*math.pi)*(E/ES)**2:.2e}")
print(f"Pais 1e6 T: B/B_S = {1e6/BS:.2e}; EH fraction = {alpha/(45*math.pi)*(1e6/BS)**2:.2e}")
wC = me * c**2 / hbar
rho_zp = hbar * wC**4 / (8 * math.pi**2 * c**3)
print(f"ZPE density to Compton cutoff = {rho_zp:.3e} J/m3; eps0 E_S^2/2 = {0.5*eps0*ES**2:.3e}; ratio = {0.5*eps0*ES**2/rho_zp:.1f} (pi/alpha = {math.pi/alpha:.1f})")
nu_pais = (eps0 * ES**2 * 8 * math.pi**2 * c**3 / hbar) ** 0.25
print(f"Pais eq.4 nu = {nu_pais:.2e} (= {nu_pais/wC:.2f} omega_C)")
A_e = 2 * alpha**2 * hbar**3 / (45 * mu0 * me**4 * c**5)
print(f"PVLAS A_e = {A_e:.3e} T^-2; Delta n = 3 A_e B^2 at 2.5 T = {3*A_e*2.5**2:.2e} (PVLAS: 2.5e-23 predicted, (12±17)e-23 measured)")

# 6. Taleyarkhan 2002
hdr("6 bubble fusion")
for nr in (4e4, 8e4, 7e5):
    print(f"{nr:.0e} n/s -> D-D power {2*nr*Qbar*MeV:.2e} W")
print(f"tritium-implied / measured neutron ratio = {7e5/8e4:.1f} to {7e5/4e4:.1f} (D-D predicts ~1)")

# 7. Lugano
hdr("7 Lugano E-Cat")
E_net = 1618194 * 3600
print(f"net 1,618,194 Wh = {E_net/1e6:.0f} MJ; per 1 g fuel {E_net/1e-3/1e6:.1e} MJ/kg (petrol ~46 MJ/kg -> x{E_net/1e-3/46e6:.1e})")
n7 = 0.011 / 7.016 * NA
print(f"p+7Li->2 4He on 0.011 g 7Li: {n7*17.35*MeV/1e6:.0f} MJ")
print(f"MFMP: File-4 excess 102.9 W (COP 1.13), File-8 289.8 W (COP 1.32); Lugano File-4 net 1658 W -> x{1658/102.9:.0f} overstatement")

# 8. Widom-Larsen
hdr("8 Widom-Larsen")
beta0 = (mn_c2 - mp_c2) / (mec2_keV / 1e3)
print(f"threshold e+p->n+nu: (mn-mp-me)c2 = {mn_c2-mp_c2-mec2_keV/1e3:.3f} MeV; beta0 = (mn-mp)/me = {beta0:.3f}")
print(f"W-L eq.5 at beta=20: {1.86e-3*(20-2.59)**2:.3f} Hz; Ciuchi et al. 1.8e-3 Hz -> ratio {1.86e-3*(20-2.59)**2/1.8e-3:.0f}")

# 9. Hydrino (Mills: W_q = 13.6 eV / q^2)
hdr("9 hydrino")
for q in (0.5, 0.25):
    print(f"q={q:.4f}: binding {13.6/q**2:.1f} eV; release from H(1) {13.6/q**2-13.6:.1f} eV")

# self-checks against published numbers
assert abs(R_per_W - 1.71e12) / 1.71e12 < 0.01
assert abs(e / (4 * math.pi * eps0 * R_D) - 300) < 20  # slides: ~300 eV
assert abs(3 * A_e * 2.5**2 - 2.5e-23) / 2.5e-23 < 0.05  # PVLAS QED prediction
assert 0.5 < r13 / R_per_W < 5  # MARRS ">=1 W/g for 1e-13 b"
assert 1e-24 < sv < 1e-23  # DD reactivity at 20 keV is of order a few 1e-24 m3/s
print("\nself-checks passed")
