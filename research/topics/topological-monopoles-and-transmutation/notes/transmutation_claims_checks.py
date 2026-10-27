"""Checks for transmutation_claims.md: nuclear accounting of the transmutation claims in
Greenyer's talk (TeslaTech, 15 Aug 2026) and in the user's audit of it.

Standard library only. Run: python3 -I transmutation_claims_checks.py  (exit 1 on failure).

Masses: AME2020, file mass_1.mas20.txt (Wang, Huang, Kondev, Audi, Naimi, Chinese Physics C 45,
030003 (2021), doi:10.1088/1674-1137/abddaf), downloaded from
https://www-nds.iaea.org/amdc/ame2020/mass_1.mas20.txt. Each entry below is (Z, A, mass excess in
keV, atomic mass in micro-u), copied from that file; the two columns are checked against each
other so a transcription slip in either one fails.
Q-values are from atomic masses: Q = sum(Delta_in) - sum(Delta_out), with Z and A conserved;
electrons then balance automatically for EC/beta-minus; each emitted positron costs 2 m_e c^2.
"""

import math
import sys

U_KEV = 931_494.102_42        # 1 u in keV (CODATA 2018)
ME_KEV = 510.998_950          # electron, keV
MP_MEV = 938.272_088_16       # proton, MeV (CODATA 2018)
MEV_J = 1.602_176_634e-13
U_KG = 1.660_539_066_60e-27
NA = 6.022_140_76e23
ALPHA = 7.297_352_5693e-3
KB_EV = 8.617_333_262e-5      # eV/K
E2_MEV_FM = 1.439_964_548     # e^2/(4 pi eps0), MeV fm
C2 = 299_792_458.0 ** 2       # J/kg
DAY = 86_400.0

AME = {
    "n": (0, 1, 8071.31806, 1008664.91590),
    "H1": (1, 1, 7288.971064, 1007825.031898),
    "H2": (1, 2, 13135.722895, 2014101.777844),
    "H3": (1, 3, 14949.81090, 3016049.28132),
    "He3": (2, 3, 14931.21888, 3016029.32197),
    "He4": (2, 4, 2424.91587, 4002603.25413),
    "Be8": (4, 8, 4941.672, 8005305.102),
    "C12": (6, 12, 0.0, 12000000.0),
    "O16": (8, 16, -4737.00217, 15994914.61926),
    "O17": (8, 17, -808.76421, 16999131.75595),
    "F17": (9, 17, 1951.701, 17002095.237),
    "Ne20": (10, 20, -7041.93217, 19992440.17525),
    "Ne22": (10, 22, -8024.716, 21991385.113),
    "Na23": (11, 23, -9529.85352, 22989769.28195),
    "Mg24": (12, 24, -13933.578, 23985041.689),
    "Al27": (13, 27, -17196.864, 26981538.408),
    "Si27": (14, 27, -12384.505, 26986704.687),
    "Si28": (14, 28, -21492.79711, 27976926.53442),
    "P30": (15, 30, -20200.856, 29978313.490),
    "P31": (15, 31, -24440.54442, 30973761.99768),
    "S31": (16, 31, -19042.531, 30979557.002),
    "Si30": (14, 30, -24432.962, 29973770.137),
    "S32": (16, 32, -26015.53714, 31972071.17354),
    "K39": (19, 39, -33807.19535, 38963706.48482),
    "Ca40": (20, 40, -34846.402, 39962590.850),
    "Ca44": (20, 44, -41468.730, 43955481.489),
    "Ca46": (20, 46, -43139.610, 45953687.726),
    "Ca48": (20, 48, -44224.868, 47952522.654),
    "Fe56": (26, 56, -60607.163, 55934935.537),
    "Co56": (27, 56, -56040.518, 55939838.032),
    "Ni56": (28, 56, -53907.649, 55942127.761),
    "Cu63": (29, 63, -65579.868, 62929597.119),
    "Cu65": (29, 65, -67263.677, 64927789.476),
    "W180": (74, 180, -49636.243, 179946713.304),
    "W182": (74, 182, -48246.144, 181948205.636),
    "W184": (74, 184, -45705.453, 183950933.180),
    "Xe132": (54, 132, -89278.97451, 131904155.08346),
    "Xe134": (54, 134, -88125.83443, 133905393.030),
    "Xe136": (54, 136, -86429.170, 135907214.474),
    "Xe140": (54, 140, -72986.461, 139921645.814),
    "Xe144": (54, 144, -56872.301, 143938945.076),
    "Pb208": (82, 208, -21748.519, 207976652.005),
    "Au197": (79, 197, -31139.751, 196966570.103),
}

failures = []


def check(name, got, want, tol, rel=False):
    err = abs(got - want) / abs(want) if rel else abs(got - want)
    ok = err <= tol
    print(f"{'ok ' if ok else 'FAIL'} {name}: {got:.6g} (want {want:.6g}, tol {tol:g}{' rel' if rel else ''})")
    if not ok:
        failures.append(name)


def mass_u(k):
    return AME[k][3] * 1e-6


def q_mev(lhs, rhs, positrons=0, weak=0):
    """Q (MeV) of a reaction given as lists of AME keys, from atomic masses.

    A must balance. Z must balance after `positrons` beta-plus decays and `weak` further charge
    changes by electron capture (weak > 0) or beta-minus (weak < 0); atomic masses carry the
    electrons for EC and beta-minus, and each positron costs 2 m_e c^2 (kinetic Q only)."""
    a = lambda side: sum(AME[k][1] for k in side)
    z = lambda side: sum(AME[k][0] for k in side)
    assert a(lhs) == a(rhs), (lhs, rhs, "A")
    assert z(lhs) - positrons - weak == z(rhs), (lhs, rhs, "Z")
    return (sum(AME[k][2] for k in lhs) - sum(AME[k][2] for k in rhs) - 2 * ME_KEV * positrons) / 1e3


# 0. The embedded table: mass excess and atomic mass columns agree (AME internal consistency).
worst = max(abs((AME[k][3] * 1e-6 - AME[k][1]) * U_KEV - AME[k][2]) for k in AME)
check("AME2020 table: max |(M - A)u - Delta| (keV)", worst, 0.0, 0.01)

# 1. Published Q-values reproduced (tests of q_mev itself).
check("D+D -> 3He+n (fusion.rs 3.268 MeV)", q_mev(["H2", "H2"], ["He3", "n"]), 3.268, 0.001)
check("D+D -> T+p (fusion.rs 4.033 MeV)", q_mev(["H2", "H2"], ["H3", "H1"]), 4.033, 0.001)
oo = ["O16", "O16"]
# Published: Woosley, Heger & Weaver RMP 74 (2002) and Clayton (1983), as listed on Wikipedia's
# "Oxygen-burning process"; tolerance 3 keV for the rounding of older mass tables, except the
# 31S + n channel, where AME2020 gives 1.497 MeV against the listed 1.459 (38 keV; see the note).
for prods, want, tol in [(["Si28", "He4"], 9.593, 0.003), (["P31", "H1"], 7.676, 0.003),
                         (["S31", "n"], 1.459, 0.04), (["Si30", "H1", "H1"], 0.381, 0.003),
                         (["P30", "H2"], -2.409, 0.003), (["S32"], 16.539, 0.003),
                         (["Mg24", "He4", "He4"], -0.393, 0.003)]:
    check(f"16O+16O -> {'+'.join(prods)}", q_mev(oo, prods), want, tol)
check("16O+16O -> 31S+n, AME2020 value (this note)", q_mev(oo, ["S31", "n"]), 1.4972, 0.0005)

# 2. Reactions on the slide and in the talk.
q_co = q_mev(["C12", "O16"], ["Si28"])
check("12C+16O -> 28Si (audit 16.76)", q_co, 16.76, 0.01)
check("12C+16O -> 24Mg+a (this note)", q_mev(["C12", "O16"], ["Mg24", "He4"]), 6.771, 0.001)
check("12C+16O -> 27Al+p (this note)", q_mev(["C12", "O16"], ["Al27", "H1"]), 5.171, 0.001)
check("12C+16O -> 27Si+n (this note)", q_mev(["C12", "O16"], ["Si27", "n"]), -0.424, 0.001)
q_sisi = q_mev(["Si28", "Si28"], ["Ni56"])
check("28Si+28Si -> 56Ni (audit 10.92, mislabelled as ->56Fe)", q_sisi, 10.92, 0.005)
q_ni_fe = q_mev(["Ni56"], ["Fe56"], weak=2)  # two EC, atomic masses
check("56Ni -> 56Fe by 2 EC (this note)", q_ni_fe, 6.6995, 0.0005)
check("56Ni -> 56Co EC (this note)", q_mev(["Ni56"], ["Co56"], weak=1), 2.1329, 0.0005)
check("56Co -> 56Fe EC (this note)", q_mev(["Co56"], ["Fe56"], weak=1), 4.5666, 0.0005)
q_sisi_fe_beta = q_mev(["Si28", "Si28"], ["Fe56"], positrons=2)
check("2 28Si -> 56Fe + 2e+ + 2nu, kinetic (coordinator 15.58)", q_sisi_fe_beta, 15.58, 0.01)
q_sisi_fe = q_sisi + q_ni_fe
check("2 28Si -> 56Fe total incl. annihilation (this note)", q_sisi_fe, 17.622, 0.001)
q_ohsawa_ni = 2 * q_co + q_sisi
check("2(12C+16O) -> 56Ni (audit 44.44 'per 56Fe')", q_ohsawa_ni, 44.43, 0.01)
check("2(12C+16O) -> 56Fe (this note)", q_ohsawa_ni + q_ni_fe, 51.133, 0.001)
q_4o17 = q_mev(["O17"] * 4, ["C12", "Fe56"])
check("4 17O -> 12C + 56Fe (coordinator 57.37; audit says 49.47)", q_4o17, 57.37, 0.005)
q_huang = q_mev(["O17", "O17"], ["C12", "Ne22"])
check("17O+17O -> 12C+22Ne (coordinator 6.41; audit says 17.13)", q_huang, 6.41, 0.005)
q_pe = q_mev(["H1", "O16"], ["O17"], weak=1)         # p + e- + 16O -> 17O + nu (Huang eq. 7)
check("p + e- + 16O -> 17O + nu (this note)", q_pe, 3.3607, 0.0005)
check("16O + p -> 17F + gamma (this note)", q_mev(["O16", "H1"], ["F17"]), 0.6003, 0.0005)
check("Huang chain per 22Ne: 2(p+e+16O) + 17O+17O (this note)", 2 * q_pe + q_huang, 13.13, 0.01)
check("Kervran 23Na+16O -> 39K (this note)", q_mev(["Na23", "O16"], ["K39"]), 19.540, 0.001)
check("Kervran 39K+1H -> 40Ca (this note)", q_mev(["K39", "H1"], ["Ca40"]), 8.328, 0.001)
check("Kervran 24Mg+16O -> 40Ca (this note)", q_mev(["Mg24", "O16"], ["Ca40"]), 16.175, 0.001)
# Tungsten -> calcium + xenon: only these pairs of stable isotopes balance A and Z.
for w, ca, xe, want in [("W184", "Ca48", "Xe136", 84.949), ("W182", "Ca48", "Xe134", 84.105),
                        ("W182", "Ca46", "Xe136", 81.323), ("W180", "Ca48", "Xe132", 83.868),
                        ("W180", "Ca40", "Xe140", 58.197), ("W184", "Ca40", "Xe144", 46.013)]:
    check(f"{w} -> {ca} + {xe} (this note)", q_mev([w], [ca, xe]), want, 0.002)
check("Parkhomov 28Si+12C -> 40Ca (this note)", q_mev(["Si28", "C12"], ["Ca40"]), 13.354, 0.001)
check("8Be -> 2 4He: 8Be unbound by (MeV, this note)", q_mev(["Be8"], ["He4", "He4"]), 0.0918, 0.0002)
# Oxygen burning: neutrons per watt from the 5% 31S+n branch, mean Q weighted by the listed branching.
q_oo_mean = 0.34 * 9.593 + 0.56 * 7.676 + 0.05 * 1.459 + 0.05 * -2.409
check("16O+16O neutrons/s per watt (this note)", 0.05 / (q_oo_mean * MEV_J), 4.16e10, 5e-3, rel=True)
# Copper -> lead / gold (Adamenko): nucleon budget, charge balanced with free p and alphas.
check("4 63Cu -> 208Pb + 5a + 24p, strong only (this note)",
      q_mev(["Cu63"] * 4, ["Pb208"] + ["He4"] * 5 + ["H1"] * 24), -427.63, 0.01)
# Least endothermic: the spare nucleons leave as alphas, with 12 EC turning protons into neutrons.
q_cu_pb = q_mev(["Cu63"] * 4, ["Pb208"] + ["He4"] * 11, weak=12)
check("4 63Cu -> 208Pb + 11a via 12 EC (this note)", q_cu_pb, -267.2, 0.1)
check("3 63Cu + 65Cu -> 197Au + 10a + 17p, strong only (this note)",
      q_mev(["Cu63"] * 3 + ["Cu65"], ["Au197"] + ["He4"] * 10 + ["H1"] * 17), -381.03, 0.01)
check("3 63Cu + 65Cu -> 197Au + 14a + n via 9 EC (this note)",
      q_mev(["Cu63"] * 3 + ["Cu65"], ["Au197"] + ["He4"] * 14 + ["n"], weak=9), -274.9, 0.1)
# Lead -> carbon (Matsumoto): 208Pb -> 17 12C + 4He needs 22 neutrons turned into protons
# (beta-minus, atomic masses carry the electrons): Q = Delta(208Pb) - 17 Delta(12C) - Delta(4He).
q_pb_c = q_mev(["Pb208"], ["C12"] * 17 + ["He4"], weak=-22)
check("208Pb -> 17 12C + 4He via 22 beta- (this note)", q_pb_c, -24.173, 0.001)
# Wendt & Irion 1922: 184W -> 46 4He needs 18 beta- (Z 74 -> 92).
q_w_he = q_mev(["W184"], ["He4"] * 46, weak=-18)
check("184W -> 46 4He via 18 beta- (this note)", q_w_he, -157.25, 0.01)

# 3. Wendt & Irion, JACS 44, 1887 (1922): 0.713 mg W fully to helium "almost exactly 4 cc";
#    explosion energy from their 0.1 uF condenser at the 31.5 kV used in most shots.
mol_w = 0.713e-3 / 183.84
check("Wendt-Irion: He volume at STP if all W -> He (paper ~4 cc)", mol_w * 46 * 22_414.0, 4.0, 0.03, rel=True)
e_shot = 0.5 * 0.1e-6 * 31.5e3 ** 2
e_need = mol_w * NA * -q_w_he * MEV_J
check("Wendt-Irion: 1/2 C V^2 per shot (J, this note)", e_shot, 49.6, 0.1)
check("Wendt-Irion: energy to split 0.713 mg W (MJ, this note)", e_need / 1e6, 58.8, 0.1)

# 4. Coulomb barriers and tunnelling (same formulas as plasmoid::fusion).
def mr_mev(a, b):
    ma, mb = mass_u(a) * U_KEV / 1e3, mass_u(b) * U_KEV / 1e3
    return ma * mb / (ma + mb)


def gamow_energy_mev(a, b):
    return 2 * mr_mev(a, b) * (math.pi * ALPHA * AME[a][0] * AME[b][0]) ** 2


# Reproduce the plasmoid::fusion doctest: D-T Gamow peak at 10 keV is 30.9 keV.
eg_dt = 2 * (1875.6129424 * 2808.9211325 / (1875.6129424 + 2808.9211325)) * (math.pi * ALPHA) ** 2
check("D-T Gamow peak at kT=10 keV (fusion.rs 30.9 keV)", (eg_dt * 1e3 * 10.0 ** 2 / 4) ** (1 / 3), 30.9, 0.1)
# Barrier: 12C+16O measured V_B ~7.8 MeV (Goasduff et al., arXiv:1107.1124) -> point-Coulomb radius.
check("12C+16O point-Coulomb radius for V_B=7.8 MeV (fm, this note)", E2_MEV_FM * 48 / 7.8, 8.862, 0.001)
for a, b, eg_want in [("C12", "O16", 15_468.0), ("O16", "O16", 32_074.0), ("O17", "O17", 34_088.0),
                      ("Si28", "Si28", 526_166.0), ("Na23", "O16", 71_520.0)]:
    eg = gamow_energy_mev(a, b)
    check(f"Gamow energy {a}+{b} (MeV, this note)", eg, eg_want, 2e-3, rel=True)
# log10 of the Maxwell-averaged tunnelling factor exp(-3 E0/kT), E0 = (E_G kT^2/4)^(1/3).
def log10_thermal(a, b, t_k):
    kt = KB_EV * t_k * 1e-6
    return -3 * (gamow_energy_mev(a, b) / (4 * kt)) ** (1 / 3) / math.log(10)


check("12C+16O at 5000 K: log10 exp(-3E0/kT) (this note)", log10_thermal("C12", "O16", 5000), -2707.5, 1.0)
check("12C+16O at 1e4 K (this note)", log10_thermal("C12", "O16", 1e4), -2149.0, 1.0)
check("16O+16O at 1e4 K (this note)", log10_thermal("O16", "O16", 1e4), -2740.4, 1.0)
check("28Si+28Si at 5000 K (this note)", log10_thermal("Si28", "Si28", 5000), -8773.0, 1.0)
check("23Na+16O at 310 K (this note)", log10_thermal("Na23", "O16", 310), -11397.0, 1.0)
# Even with 1 keV of screening (far beyond measured U_e of a few hundred eV) at E = 1 keV:
check("12C+16O tunnelling at E=1 keV: log10 exp(-sqrt(E_G/E)) (this note)",
      -math.sqrt(gamow_energy_mev("C12", "O16") / 1e-3) / math.log(10), -1708.0, 2.0)

# 5. Srinivasan: "4 tonnes of excess metal per day for 11 weeks", 11 MW furnace.
fe_atoms_day = 4e6 / mass_u("Fe56") * NA   # alpha-conjugate routes make pure 56Fe (55.935 u)
p_si = fe_atoms_day * q_sisi * MEV_J / DAY            # if Fe came from 28Si+28Si -> 56Ni
p_co = fe_atoms_day * (q_ohsawa_ni + q_ni_fe) * MEV_J / DAY   # all the way from C+O
check("4 t Fe/day: atoms/day (this note)", fe_atoms_day, 4.3066e28, 1e-4, rel=True)
check("4 t Fe/day via Si+Si: power (W, this note)", p_si, 8.722e11, 1e-3, rel=True)
check("4 t Fe/day via 2(C+O): power (W, this note)", p_co, 4.083e12, 1e-3, rel=True)
check("ratio to an 11 MW furnace, Si+Si (this note)", p_si / 11e6, 7.929e4, 1e-3, rel=True)
check("4 t over 11 weeks instead: Si+Si power (W, this note)", p_si / 77, 1.1327e10, 1e-3, rel=True)
# Steady-state 56Ni/56Co activity if the iron passed through 56Ni (production = decay rate).
check("steady-state 56Ni activity (Ci, this note)", fe_atoms_day / DAY / 3.7e10, 1.3472e13, 1e-3, rel=True)
# Etemadi et al. Heliyon 9, e13956 (2023): FeSi75 at 8.46-9.18 MWh/t -> an 11 MW furnace's output.
check("11 MW furnace FeSi75 output at 8.46 MWh/t (t/day)", 11 * 24 / 8.46, 31.2, 0.05)
check("11 MW furnace FeSi75 output at 9.18 MWh/t (t/day)", 11 * 24 / 9.18, 28.76, 0.05)
# 17O budget for 4 17O -> 12C + 56Fe at 4 t/day (natural 17O 0.038%, NIST).
o_total_t = 4 * 4e6 / mass_u("Fe56") / 0.00038 * 15.999 / 1e6
check("oxygen needed per day for 4 t Fe via 4 17O (tonnes, this note)", o_total_t, 12_043.0, 1e-3, rel=True)

# 6. Adamenko: 300 J (talk) per shot at 100% conversion to the endothermic Cu -> Pb.
pb_atoms = 300 / (-q_cu_pb * MEV_J)
check("Pb atoms per 300 J shot at most (this note)", pb_atoms, 7.006e12, 2e-3, rel=True)
check("Pb mass per 300 J shot at most (ng, this note)", pb_atoms * 207.98 / NA * 1e9, 2.42, 0.01)
d_um = (6 * pb_atoms * 207.98 / NA / 11.34 / math.pi) ** (1 / 3) * 1e4
check("diameter of that much Pb as a sphere (um, rho 11.34, this note)", d_um, 7.4, 0.05)

# 7. Sundaresan & Bockris / Ohsawa: energy and activity per microgram of 56Fe made from 2(C+O).
fe_ug = 1e-6 / 55.935 * NA
check("2(C+O) -> 56Fe: heat per ug Fe (kJ, this note)", fe_ug * (q_ohsawa_ni + q_ni_fe) * MEV_J / 1e3, 88.2, 0.1)
check("56Ni activity per ug if made as 56Ni (Ci, t1/2 6 d, this note)",
      fe_ug * math.log(2) / (6 * DAY) / 3.7e10, 0.389, 0.002)

# 8. Monopole "938 MeV per proton ... 16 times for oxygen ... 166 times D-D".
e16 = 16 * MP_MEV
e16_true = mass_u("O16") * U_KEV / 1e3
check("16 m_p c^2 (audit 15012.35 MeV)", e16, 15_012.35, 0.01)
check("M(16O atom) c^2, the ceiling (MeV, this note)", e16_true, 14_899.17, 0.01)
check("audit specific energy / c^2 (this note: exceeds E=mc^2)",
      e16 * MEV_J / (16 * U_KG) / C2, 1.00728, 1e-5)
dd = (3.268 + 4.033) / 2
check("m_p / D-D mean Q per reaction (this note)", MP_MEV / dd, 257.0, 0.5)
check("m_p / D-D per nucleon (this note)", MP_MEV / (dd / 4), 1028.0, 1.0)
check("(m_p/16) / (D-D per nucleon): one proton per 16O (this note)", MP_MEV / 16 / (dd / 4), 64.25, 0.05)
check("m_p / (D+D->4He 23.85 MeV / 4) (this note)",
      MP_MEV / (q_mev(["H2", "H2"], ["He4"]) / 4), 157.4, 0.2)

# 9. The audit's specific energies per kg of oxygen, corrected.
def per_kg_o(q, n_o, key):
    return q * MEV_J / (n_o * mass_u(key) * U_KG)


check("4 17O -> 12C+56Fe: J per kg 17O (this note)", per_kg_o(q_4o17, 4, "O17"), 8.133e13, 2e-3, rel=True)
check("2(C+O) -> 56Fe: J per kg 16O (this note)", per_kg_o(q_ohsawa_ni + q_ni_fe, 2, "O16"), 1.540e14, 2e-3, rel=True)
check("17O+17O -> 12C+22Ne: J per kg 17O (this note)", per_kg_o(q_huang, 2, "O17"), 1.818e13, 2e-3, rel=True)
# What the audit's own specific_energy_joules_per_kg_o2 returns (Q / (n_O * A * u)).
for name, q, n, a, want in [("Ohsawa", 49.47, 4, 17, 7.021e13), ("Srinivasan", 44.44, 2, 16, 1.340e14),
                            ("Huang", 17.13, 2, 17, 4.862e13), ("monopole", 16 * MP_MEV, 1, 16, 9.053e16)]:
    check(f"audit formula, {name} (J/kg, this note)", q * MEV_J / (n * a * U_KG), want, 2e-3, rel=True)
check("whole-atom annihilation of 16O: J/kg = c^2", per_kg_o(e16_true, 1, "O16"), C2, 1e-9, rel=True)
check("1 kg per 33 days at c^2: power (GW, this note)", C2 / (33 * DAY) / 1e9, 31.52, 0.02)

# 10. O2 molar susceptibility: CRC +3449e-6 cm^3/mol (CGS) vs Moldover et al. PRA 78 (2008),
#     (42.92 +/- 0.06)e-9 m^3/mol at 293.15 K. SI = 4 pi 1e-6 x CGS.
check("O2: CRC 3449e-6 cgs in SI vs Moldover 42.92e-9", 3449e-6 * 4 * math.pi * 1e-6, 42.92e-9, 0.015, rel=True)
# Wolfram ElementData (periodictable.com), SI m^3/mol: Tb 2.161385e-6, Dy 8.85625e-7, O2 4.27184e-8.
check("Tb / O2 per mole (this note)", 2.161385e-6 / 4.27184e-8, 50.6, 0.1)
check("Dy / O2 per mole (this note)", 8.85625e-7 / 4.27184e-8, 20.73, 0.02)
# Curie law: chi ~ 1/T, so the audit's "at STP" (273.15 K) value is 7% above the 293 K one.
check("O2 at 273.15 K from 3449e-6 at 293 K, Curie (cgs, this note)", 3449e-6 * 293 / 273.15, 3699.6e-6, 1e-3, rel=True)

# 11. Isotopic signatures (NIST isotopic compositions): natural 22Ne/20Ne and the Fe pattern.
check("natural 22Ne/20Ne (NIST 0.0925/0.9048)", 0.0925 / 0.9048, 0.1022, 0.0002)

# 12. Shoulders: an unneutralised EV of 1e11 electrons, 1 um across (US 5,018,180).
q_ev = 1e11 * 1.602_176_634e-19
check("EV surface potential (MV, this note)", 8.987_551_79e9 * q_ev / 0.5e-6 / 1e6, 288.0, 0.5)

if failures:
    print(f"\n{len(failures)} FAILED: {failures}")
    sys.exit(1)
print("\nall checks passed")
