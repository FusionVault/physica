"""Checks for propulsion_and_energy.md (Greenyer TeslaTech 2026 talk: ARV, oxygen fuel,
relic neutrinos, Casimir, HFGW).

Standard library only. Run:  python3 -I propulsion_and_energy_checks.py
Every computed number is asserted against a published value (cited) or against the value
quoted in the notes (marked "notes"), with a stated relative tolerance. Exit code 1 on failure.
"""
import math
import sys

FAIL = []


def check(name, got, want, rel):
    ok = abs(got / want - 1.0) <= rel
    print(f"{'ok  ' if ok else 'FAIL'} {name}: {got:.6g} (want {want:.6g}, tol {rel:g})")
    if not ok:
        FAIL.append(name)


# CODATA 2018 / SI exact
C = 299_792_458.0
HBAR = 1.054_571_817e-34
H = 6.626_070_15e-34
E = 1.602_176_634e-19
KB = 1.380_649e-23
G = 6.674_30e-11
EPS0 = 8.854_187_812_8e-12
MU0 = 1.256_637_062_12e-6
ME = 9.109_383_701_5e-31
U_KG = 1.660_539_066_60e-27
U_KEV = 931_494.102_42
MEV = 1e6 * E
G0 = 9.806_65
DAY = 86_400.0
MPC = 3.085_677_581_491_367e22

# ---------------------------------------------------------------------------
# 1. Oxygen as fuel: AME2020 atomic mass excesses (keV), Wang et al. 2021, mass_1.mas20
# ---------------------------------------------------------------------------
D = {
    "n": 8071.31806, "1H": 7288.971064, "2H": 13135.722895, "3H": 14949.81090, "3He": 14931.21888, "4He": 2424.91587, "12C": 0.0,
    "16O": -4737.00217, "17O": -808.76421, "22Ne": -8024.716, "28Si": -21492.79711,
    "31P": -24440.54442, "31S": -19042.531, "32S": -26015.53714, "56Fe": -60607.163,
}
A = {"n": 1, "1H": 1, "2H": 2, "3H": 3, "3He": 3, "4He": 4, "12C": 12, "16O": 16, "17O": 17, "22Ne": 22, "28Si": 28,
     "31P": 31, "31S": 31, "32S": 32, "56Fe": 56}


def q_kev(lhs, rhs):
    """Q (keV) from atomic mass excesses; electrons balance when Z balances."""
    return sum(D[x] for x in lhs) - sum(D[x] for x in rhs)


def atomic_mass_u(x):
    return A[x] + D[x] / U_KEV


q_s32 = q_kev(["16O", "16O"], ["32S"])
q_si_a = q_kev(["16O", "16O"], ["28Si", "4He"])
q_p31 = q_kev(["16O", "16O"], ["31P", "1H"])
q_s31 = q_kev(["16O", "16O"], ["31S", "n"])
check("16O+16O -> 32S Q (MeV) [notes; AME2020]", q_s32 / 1e3, 16.5415, 1e-4)
check("16O+16O -> 28Si+a Q (MeV) [notes]", q_si_a / 1e3, 9.5939, 1e-4)
check("16O+16O -> 31P+p Q (MeV) [notes]", q_p31 / 1e3, 7.6776, 1e-4)
check("16O+16O -> 31S+n Q (MeV) [notes]", q_s31 / 1e3, 1.4972, 1e-3)
# The coordinator's AME2020 screen of the user's audit values
q_ohsawa = q_kev(["17O"] * 4, ["56Fe", "12C"])
q_huang = q_kev(["17O"] * 2, ["12C", "22Ne"])
check("4 17O -> 56Fe+12C Q (MeV) [coordinator screen 57.37; audit says 49.47]", q_ohsawa / 1e3, 57.37, 1e-3)
check("2 17O -> 12C+22Ne Q (MeV) [coordinator screen 6.41; audit says 17.13]", q_huang / 1e3, 6.41, 1e-3)
check("12C+16O -> 28Si Q (MeV) [audit 16.76]", q_kev(["12C", "16O"], ["28Si"]) / 1e3, 16.76, 1e-3)

m16 = atomic_mass_u("16O")
# energy fractions eps = Q / (rest energy of the reactants)
eps_oo = q_s32 / (2 * m16 * U_KEV)
check("16O+16O -> 32S mass fraction released [notes]", eps_oo, 5.551e-4, 1e-3)
e_oo = eps_oo * C * C
check("16O+16O -> 32S specific energy (J/kg) [notes]", e_oo, 4.989e13, 1e-3)
e_ohsawa = q_ohsawa * 1e3 * E / (4 * atomic_mass_u("17O") * U_KG)
check("4 17O -> Fe+C specific energy (J/kg of 17O) [notes]", e_ohsawa, 8.14e13, 2e-3)
# chemistry: H2 + 1/2 O2 -> H2O(g), dHf = -241.826 kJ/mol (CODATA, NIST WebBook)
M_H2O, M_O, M_H2 = 18.01528e-3, 15.9994e-3, 2.01588e-3
e_chem_mix = 241.826e3 / M_H2O
e_chem_o = 241.826e3 / M_O
check("H2/O2 to H2O(g) per kg propellant (J/kg) [notes]", e_chem_mix, 1.3424e7, 1e-3)
check("H2/O2 per kg of oxygen (J/kg O) [notes]", e_chem_o, 1.5115e7, 1e-3)
eps_chem = e_chem_mix / C**2
check("chemistry mass fraction [notes]", eps_chem, 1.4936e-10, 1e-3)
# complete conversion; the talk's '938 MeV per proton, 16 per oxygen'
check("c^2 (J/kg)", C * C, 8.98755e16, 1e-5)
talk_baryon = 16 * 938.272_088_16 * MEV / (16 * U_KG)          # what the audit code does
true_full = C * C
check("talk/audit '16 x m_p c^2 per 16 u' over E=mc^2 of 16O [notes: exceeds by 0.76%]",
      16 * 938.272_088_16 / (m16 * U_KEV / 1e3), 1.00766, 1e-4)
check("talk/audit baryon-decay specific energy (J/kg) [notes]", talk_baryon, 9.0530e16, 1e-4)
# D-D: the two branches from AME2020, taken as equally likely
q_dd_n = q_kev(["2H", "2H"], ["3He", "n"]) / 1e3
q_dd_p = q_kev(["2H", "2H"], ["3H", "1H"]) / 1e3
check("D+D -> 3He+n Q (MeV)", q_dd_n, 3.2689, 1e-4)
check("D+D -> T+p Q (MeV)", q_dd_p, 4.0327, 1e-4)
dd_mean = 0.5 * (q_dd_n + q_dd_p)
check("m_p c^2 / mean D-D Q (per reaction) [notes: 257, not 166]", 938.272 / dd_mean, 257.0, 2e-3)
check("full conversion / D-D per unit mass [notes: ~1028]", 931.494 / (dd_mean / 4.0282), 1028.0, 1e-2)

# ---------------------------------------------------------------------------
# 2. Rocket accounting: craft of 5000 kg (audit) carrying 950 kg O (talk)
# ---------------------------------------------------------------------------
M0, MO = 5000.0, 950.0
R = M0 / (M0 - MO)


def u_exhaust(eps):
    """Rocket-frame exhaust speed if fraction eps of propellant rest mass becomes exhaust KE."""
    return C * math.sqrt(eps * (2 - eps))


def dv_rel(u, ratio):
    """Relativistic rocket: dv = c tanh((u/c) ln R) (Ackeret 1946)."""
    return C * math.tanh(u / C * math.log(ratio))


check("photon rocket dv/c = (R^2-1)/(R^2+1)", dv_rel(C, R) / C, (R * R - 1) / (R * R + 1), 1e-12)
check("photon rocket dv/c with 950 kg of 5000 kg [notes]", dv_rel(C, R) / C, 0.2077, 1e-3)
check("O+O fusion exhaust speed (m/s) [notes]", u_exhaust(eps_oo), 9.988e6, 1e-3)
check("O+O fusion dv (m/s) [notes]", dv_rel(u_exhaust(eps_oo), R), 2.104e6, 2e-3)
r_chem = (M0 - MO + MO * (1 + M_H2 / M_O)) / (M0 - MO)  # add 119.7 kg H2
check("H2/O2 ideal exhaust speed (m/s) [notes]", u_exhaust(eps_chem), 5181.0, 1e-3)
check("H2/O2 dv with 950 kg O + 120 kg H2 (m/s) [notes]", dv_rel(u_exhaust(eps_chem), r_chem), 1214.0, 2e-3)
# constant proper acceleration with a perfect photon rocket: tau = (c/a) ln R
tau_1g = C / G0 * math.log(R) / DAY
check("photon rocket at 1 g: proper time (days) [notes]", tau_1g, 74.5, 2e-3)
a_28 = C * math.log(R) / (28 * DAY)
check("acceleration if the 950 kg lasts 28 days as a photon rocket (g) [notes]", a_28 / G0, 2.66, 5e-3)
p_28 = MO * C * C / (28 * DAY)
check("photon-rocket power to burn 950 kg in 28 days (W) [notes]", p_28, 3.529e13, 1e-3)
p_1kg = C * C / (31.5 * DAY)
check("power, 1 kg fully converted over 31.5 days (W) [notes]", p_1kg, 3.302e10, 1e-3)
check("photon thrust from it (N) [notes]", p_1kg / C, 110.2, 2e-3)
check("photon thrust from 1 kW (N) [repo propulsion::photon_thrust test]", 1e3 / C, 3.336e-6, 1e-3)

# ---------------------------------------------------------------------------
# 3. The user's audit simulator, with its own numbers
# ---------------------------------------------------------------------------
mdot = 1.0 / (33.0 * DAY)
v_eff = 90e3 * 3.5
p_audit = mdot * talk_baryon
f_audit = p_audit / v_eff
check("audit burn rate (kg/s) [audit 1 kg/33 d]", mdot, 3.507e-7, 1e-3)
check("audit power, baryon pathway (W) [notes]", p_audit, 3.176e10, 1e-3)
check("audit v_eff = f R0 (m/s)", v_eff, 3.15e5, 1e-12)
check("audit thrust F = P/v_eff (N) [notes]", f_audit, 1.008e5, 1e-3)
check("audit acceleration on 5000 kg (g) [notes]", f_audit / M0 / G0, 2.056, 1e-3)
check("audit thrust over photon bound P/c = c/v_eff [notes]", f_audit / (p_audit / C), 951.7, 1e-3)
check("max thrust from mdot if all of it left at c (N) [notes]", mdot * C, 105.1, 2e-3)
p_audit_fe = mdot * e_ohsawa
check("audit power, corrected 4 17O pathway (W) [notes]", p_audit_fe, 2.856e7, 2e-3)
check("audit thrust, corrected 4 17O pathway (N) [notes]", p_audit_fe / v_eff, 90.66, 2e-3)
check("tank life at the audit burn rate (years) [notes]", MO / mdot / (365.25 * DAY), 85.83, 1e-3)
check("Newtonian integrator passes c after (days) [notes]", C / (f_audit / M0) / DAY, 172.2, 2e-3)
# air-breathing ideal actuator disk: P = F^1.5 / sqrt(2 rho A)
rho_air, A_disk = 1.225, math.pi * 3.5**2
f_hover = M0 * G0
p_hover = f_hover**1.5 / math.sqrt(2 * rho_air * A_disk)
check("ideal power to hover 5000 kg on a 3.5 m-radius disk in air (W) [notes]", p_hover, 1.119e6, 2e-3)
check("ideal air thrust from the audit's 31.8 GW (N) [notes]", (2 * rho_air * A_disk * p_audit**2) ** (1 / 3), 4.57e7, 3e-3)
# audit's displacement current: 1 uF, 2.5 MV, 250 ns
i_disp = 1e-6 * 2.5e6 / 250e-9
check("audit displacement current C dV/dt (A) [notes]", i_disp, 1.0e7, 1e-12)
check("its stored energy 1/2 C V^2 (J) [notes]", 0.5 * 1e-6 * 2.5e6**2, 3.125e6, 1e-12)
check("its power if dumped in 250 ns (W) [notes]", 0.5 * 1e-6 * 2.5e6**2 / 250e-9, 1.25e13, 1e-12)
check("plate area per metre of gap for 1 uF (m^2) [notes]", 1e-6 / EPS0, 1.129e5, 1e-3)

# ---------------------------------------------------------------------------
# 4. The 48 plates at 2.5 MV, 30-150 kHz, and 7.3 TA
# ---------------------------------------------------------------------------
A_ann = math.pi * (3.5**2 - 0.5**2)
gap = 1.0
E_field = 2.5e6 / gap
cap = EPS0 * A_ann / gap
check("annulus 0.5-3.5 m area (m^2)", A_ann, 37.70, 1e-3)
check("electrostatic pressure eps0 E^2/2 at 2.5 MV/m (Pa) [notes]", 0.5 * EPS0 * E_field**2, 27.67, 1e-3)
check("force on 37.7 m^2 (N) [notes]", 0.5 * EPS0 * E_field**2 * A_ann, 1043.0, 2e-3)
check("capacitance (F) [notes]", cap, 3.338e-10, 1e-3)
check("stored energy (J) [notes]", 0.5 * cap * 2.5e6**2, 1043.0, 2e-3)
check("displacement current amplitude at 90 kHz (A) [notes]", cap * 2.5e6 * 2 * math.pi * 90e3, 471.9, 2e-3)
I_T = 7.3e12
B35 = MU0 * I_T / (2 * math.pi * 3.5)
check("B of 7.3 TA at 3.5 m (T) [notes]", B35, 4.171e5, 1e-3)
check("magnetic pressure there (Pa) [notes]", B35**2 / (2 * MU0), 6.924e16, 2e-3)
check("over the 1200 T record field (Nakamura et al. 2018)", B35 / 1200.0, 347.6, 2e-3)
L_ring = MU0 * 3.5 * (math.log(8 * 3.5 / 0.4) - 2.0)
U_ring = 0.5 * L_ring * I_T**2
check("ring inductance R=3.5 m, a=0.4 m (H) [notes]", L_ring, 9.889e-6, 1e-3)
check("magnetic energy of 7.3 TA ring (J) [notes]", U_ring, 2.635e20, 2e-3)
check("that over 950 kg x c^2 [notes]", U_ring / (MO * C * C), 3.086, 3e-3)
check("electron flow of 7.3 TA (kg/s of electrons) [notes]", I_T / E * ME, 41.5, 2e-3)
check("Alfven current 4 pi eps0 m c^3/e (A) (Alfven 1939)", 4 * math.pi * EPS0 * ME * C**3 / E, 1.7045e4, 1e-3)
check("beta-gamma for 7.3 TA beam below Alfven limit [notes]", I_T / (4 * math.pi * EPS0 * ME * C**3 / E), 4.283e8, 2e-3)

# ---------------------------------------------------------------------------
# 5. Relic neutrinos
# ---------------------------------------------------------------------------
T_CMB = 2.7255                       # Fixsen 2009, as used by Planck 2018
T_NU = (4 / 11) ** (1 / 3) * T_CMB
zeta3 = 1.202056903159594
n_gamma = 2 * zeta3 / math.pi**2 * (KB * T_CMB / (HBAR * C)) ** 3 / 1e6   # cm^-3
n_nu_flav = 3 / 11 * n_gamma
check("T_nu (K) [PDG 2024 ch. 26: ~1.9 K]", T_NU, 1.9, 0.03)
check("T_nu (K) [notes, exact (4/11)^(1/3) T_CMB]", T_NU, 1.9454, 1e-4)
check("n_gamma (cm^-3) [notes]", n_gamma, 410.7, 1e-3)
check("n_nu, all 3 flavours nu+nubar, instant decoupling (cm^-3) [notes 336]", 3 * n_nu_flav, 336.0, 1e-3)
check("vs PDG 2024 n_nu = 339.5 cm^-3 (with QED/decoupling corrections)", 3 * n_nu_flav, 339.5, 0.012)
omega_nu_desi = 0.064 / 93.14
check("Omega_nu h^2 at DESI DR2 limit 0.064 eV [notes]", omega_nu_desi, 6.871e-4, 1e-3)
check("Omega_nu/Omega_b at DESI limit (Omega_b h^2 = 0.0224) [notes]", omega_nu_desi / 0.0224, 0.03068, 2e-3)
check("Omega_nu/Omega_b at Planck limit 0.12 eV [notes]", 0.12 / 93.14 / 0.0224, 0.05752, 2e-3)
check("sum m_nu for 'four to five times' (4.5x) baryons (eV) [notes]", 4.5 * 0.0224 * 93.14, 9.389, 1e-3)
check("that over DESI limit [notes]", 4.5 * 0.0224 * 93.14 / 0.064, 146.7, 2e-3)
check("Omega_c h^2 / Omega_b h^2 (Planck 2018 0.120/0.0224) [notes: 5.4]", 0.120 / 0.0224, 5.357, 1e-3)
# rho_crit / h^2 from 3 H^2 / (8 pi G) and the neutrino energy density
h100 = 1e5 / MPC
rho_c = 3 * h100**2 / (8 * math.pi * G)
check("rho_crit/h^2 (kg/m^3)", rho_c, 1.8785e-26, 1e-3)
u_nu = omega_nu_desi * rho_c * C**2
check("relic-neutrino energy density at 0.064 eV (J/m^3) [notes]", u_nu, 1.160e-12, 2e-3)
v_cmb = 370e3   # input: order of the CMB-dipole speed of the Sun
check("ram pressure rho_nu v^2 at 370 km/s, total absorption (Pa) [notes]", u_nu * (v_cmb / C) ** 2, 1.765e-18, 3e-3)
# Duda, Gelmini & Nussinov 2001: largest accelerations ~ f x 1e-26 cm/s^2, f <~ 17
a_dgn = 17 * 1e-26 * 1e-2
check("g over the largest relic-neutrino acceleration [notes]", G0 / a_dgn, 5.769e27, 3e-3)
# Pauli (Tremaine-Gunn type) ceiling on a bound population: n <= g (4 pi/3) (m v / h)^3
m_ev, v_esc = 0.05, 550e3
p_ev = m_ev * v_esc / C
hc_ev_cm = H * C / E * 100
n_max = 2 * 4 * math.pi / 3 * (p_ev / hc_ev_cm) ** 3
check("Pauli ceiling, 0.05 eV neutrinos bound at 550 km/s (cm^-3) [notes]", n_max, 3.39, 5e-3)

# ---------------------------------------------------------------------------
# 6. Casimir and vacuum energy
# ---------------------------------------------------------------------------
def casimir_p(d):
    return math.pi**2 * HBAR * C / (240 * d**4)


def casimir_e_area(d):
    return math.pi**2 * HBAR * C / (720 * d**3)


check("Casimir pressure 1 um (Pa) [notes]", casimir_p(1e-6), 1.300e-3, 2e-3)
check("Casimir pressure 10 nm (Pa) [notes]", casimir_p(1e-8), 1.300e5, 2e-3)
check("Casimir energy per area at 10 nm (J/m^2) [notes]", casimir_e_area(1e-8), 4.334e-4, 2e-3)
check("Casimir energy per area at 1 um (J/m^2) [notes]", casimir_e_area(1e-6), 4.334e-10, 2e-3)
# area of 10 nm-gap plates needed to equal 1 kg of H2/O2 chemistry, one-shot
check("area to match 1 kg of H2/O2, one-shot at 10 nm (km^2) [notes]", e_chem_mix / casimir_e_area(1e-8) / 1e6, 3.097e4, 3e-3)

# ---------------------------------------------------------------------------
# 7. Gravitational waves from the craft (repo gravitational::laboratory formulas)
# ---------------------------------------------------------------------------
def rod_power(m, l, w):
    return 2 / 45 * G * m * m * l**4 * w**6 / C**5


p_rod = rod_power(5000.0, 7.0, 2 * math.pi * 10)
check("GW power, 5000 kg 7 m rod at 10 rev/s (W) [notes]", p_rod, 4.524e-33, 3e-3)
# field energy of the plates (1043 J) as mass, oscillating across 7 m at 90 kHz: P ~ (G/c^5)(M L^2 w^3)^2
m_field = 1043.0 / C**2
w90 = 2 * math.pi * 90e3
p_field = G / C**5 * (m_field * 3.5**2 * w90**3) ** 2
check("GW power, plate field energy pulsing at 90 kHz (W, order) [notes]", p_field, 1.82e-44, 0.05)
check("G/c^5 (W^-1)", G / C**5, 2.7562e-53, 1e-3)

print()
if FAIL:
    print(f"{len(FAIL)} FAILED:", *FAIL, sep="\n  ")
    sys.exit(1)
print("all checks passed")
