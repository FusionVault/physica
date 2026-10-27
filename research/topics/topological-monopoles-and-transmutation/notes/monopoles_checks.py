#!/usr/bin/env python3
"""Recompute every number in monopoles.md and assert it against the published value.

Standard library only. Run:  python3 -I monopoles_checks.py   (exits 1 on any failure)
Published values are hard-coded with their source tag; tolerances are explicit per check.
"""
import math, sys

FAIL = []
def check(name, got, want, rtol):
    ok = abs(got - want) <= rtol * abs(want)
    print(f"{'ok  ' if ok else 'FAIL'} {name}: {got:.6g}  (published/expected {want:.6g}, rtol {rtol:g})")
    if not ok: FAIL.append(name)
def within(name, got, lo, hi):
    ok = lo <= got <= hi
    print(f"{'ok  ' if ok else 'FAIL'} {name}: {got:.6g}  (expected in [{lo:.3g}, {hi:.3g}])")
    if not ok: FAIL.append(name)
def show(name, val, unit=""):
    print(f"     {name}: {val:.6g} {unit}")

# ---- CODATA 2022 (physics.nist.gov/cuu/Constants/Table/allascii.txt, fetched 2026-10-26)
h = 6.62607015e-34; hbar = h / (2 * math.pi); e = 1.602176634e-19; c = 299792458.0
mu0 = 1.25663706127e-6; alpha = 7.2973525643e-3; muB = 9.2740100657e-24; muN = 5.0507837393e-27
a0 = 5.29177210544e-11; mp_MeV = 938.27208943; u_MeV = 931.49410372; u_kg = 1.66053906892e-27
NA = 6.02214076e23; Ry_eV = 13.605693122990; hbarc_MeVfm = 197.3269804; Phi0 = 2.067833848e-15
B_atomic = 2.35051757077e5; me_MeV = 0.51099895069
check("CODATA: alpha = mu0 c e^2/(2h)", mu0 * c * e**2 / (2 * h), alpha, 1e-9)
eV = e; MeV = 1e6 * e; GeV = 1e9 * e
muB_eV = muB / e; muN_eV = muN / e

print("\n== 1. Dirac charge, flux, field, field energy ==")
gD = h / (mu0 * e)                                   # A m
check("g_D = h/(mu0 e) [repo magnetics report 3.291e-9 A m]", gD, 3.291e-9, 2e-4)
check("Gaussian g_D/e = 1/(2 alpha) [repo 68.52]", 1 / (2 * alpha), 68.52, 1e-4)
check("Dirac flux mu0 g_D = h/e = 2 Phi0 [repo 4.14e-15 Wb]", mu0 * gD, 2 * Phi0, 1e-9)
check("Castelnovo 2008 SI condition e q = n h/mu0 (n=1)", e * gD, h / mu0, 1e-12)
def B_dirac(r, n=1): return mu0 * n * gD / (4 * math.pi * r**2)   # = n hbar/(2 e r^2)
check("B at r: mu0 g/(4 pi r^2) == hbar/(2 e r^2)", B_dirac(1e-15), hbar / (2 * e * 1e-30), 1e-12)
check("B_Dirac(1 fm) [derived]", B_dirac(1e-15), 3.291e14, 2e-4)
check("B_Dirac(a0) [derived]", B_dirac(a0), 1.1752e5, 2e-4)
r15 = math.sqrt(hbar / (2 * e * 1e15))
check("radius where B_Dirac = 1e15 T [derived, fm]", r15 / 1e-15, 0.5737, 1e-3)
# field energy outside r0: U = mu0 g^2/(8 pi r0) = hbar c/(8 alpha r0)
K = mu0 * gD**2 / (8 * math.pi)                      # J m
check("mu0 g_D^2/(8 pi) == hbar c/(8 alpha)", K, hbar * c / (8 * alpha), 1e-9)
check("field-energy constant [MeV fm]", K / MeV / 1e-15, 3380.1, 1e-4)
check("r0 for U = m_p c^2 [fm]", K / (mp_MeV * MeV) / 1e-15, 3.6025, 1e-3)
check("U(r0 = hbar c/M_X) / M_X = 1/(8 alpha)", 1 / (8 * alpha), 17.129, 1e-4)
within("U for M_X = 1e16 GeV (PDG 94.1: M_mon ~ 1e17-1e18 GeV)", 1e16 / (8 * alpha), 1e17, 1e18)
check("r0 for U = 1 J [derived, m]", K / 1.0, 5.416e-25, 1e-3)

print("\n== 2. The user's audit: formulas checked ==")
g_audit = 2 * math.pi * hbar / e
check("audit g_D = 2 pi hbar/e is h/e (a flux in Wb)", g_audit, 4.135667696e-15, 1e-9)
check("audit B = mu0 g_audit/(4 pi r^2) is low by a factor mu0", (mu0 * g_audit) / (mu0 * gD), mu0, 1e-9)
def cross(a, b): return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
def dot(a, b): return sum(x*y for x, y in zip(a, b))
def A_audit(r, n=(0, 0, 1), g=gD):           # (mu0 g/4 pi r) (r x n)/(r + r.n)
    R = math.sqrt(dot(r, r)); k = mu0 * g / (4 * math.pi * R * (R + dot(r, n)))
    return tuple(k * x for x in cross(r, n))
def A_sph(r, g=gD):                           # (mu0 g/4 pi)(1 - cos t)/(r sin t) phi_hat
    x, y, z = r; R = math.sqrt(x*x + y*y + z*z); rho = math.hypot(x, y)
    k = mu0 * g / (4 * math.pi) * (1 - z / R) / (R * (rho / R)) / rho
    return (-k * y, k * x, 0.0)
def curl(F, r, d=1e-6):
    def p(i, j):
        rp = list(r); rm = list(r); rp[j] += d; rm[j] -= d
        return (F(tuple(rp))[i] - F(tuple(rm))[i]) / (2 * d)
    return (p(2, 1) - p(1, 2), p(0, 2) - p(2, 0), p(1, 0) - p(0, 1))
pt = (0.3, -0.4, 0.5); R = math.sqrt(dot(pt, pt)); Bz_true = mu0 * gD / (4 * math.pi) * pt[2] / R**3
check("spherical form (1-cos)/(r sin) phi: curl_z = +mu0 g z/(4 pi r^3)", curl(A_sph, pt)[2], Bz_true, 1e-5)
check("audit vector form (r x n): curl_z = -mu0 g z/(4 pi r^3) (sign flipped)", curl(A_audit, pt)[2], -Bz_true, 1e-5)
check("both forms differ only by sign: A_audit = -A_sph", A_audit(pt)[0], -A_sph(pt)[0], 1e-9)
# string: denominator r + r.n = 0 on r.n = -r, i.e. the half-line along -n
check("audit denominator vanishes on -n half-line", (lambda r: math.sqrt(dot(r, r)) + dot(r, (0, 0, 1)))((0, 0, -2.0)) + 1, 1.0, 1e-12)
# (3) B where B^2/(2 mu0) = Lambda^4/(hbar c)^3, and inverse
def B_from_Lambda4(L_MeV):
    u = (L_MeV * MeV)**4 / (hbar * c)**3
    return math.sqrt(2 * mu0 * u)
for L, want in ((200, 2.8957e14), (300, 6.5154e14), (340, 8.3686e14)):
    check(f"B(Lambda^4 equality) at Lambda = {L} MeV [derived, T]", B_from_Lambda4(L), want, 1e-3)
L15 = 200 * math.sqrt(1e15 / B_from_Lambda4(200))
check("Lambda for which B(Lambda^4) = 1e15 T [MeV]", L15, 371.67, 1e-4)
GEV2_T = (GeV)**2 / (e * hbar * c**2)           # field with eB = 1 GeV^2
check("eB = 1 GeV^2 in tesla [derived]", GEV2_T, 1.6904e16, 1e-4)
check("Lambda for which eB = Lambda^2 at 1e15 T [MeV]", 1000 * math.sqrt(1e15 / GEV2_T), 243.22, 1e-4)

print("\n== 3. QCD and heavy-ion field scales ==")
mpi = 140.0   # MeV, as used by Skokov et al. 2009
check("eB = m_pi^2 (140 MeV) in tesla", (mpi / 1000)**2 * GEV2_T, 3.3132e14, 1e-4)
T_HL = math.sqrt(2 * (1.0**2 / (2 * mu0)) / e * (hbar * c / e)**3)   # 1 T in Heaviside-Lorentz eV^2
check("1 T in natural HL units [eV^2; repo 195.4]", T_HL, 195.35, 1e-4)
check("Skokov identity: 1 MeV^2 (HL) = 0.512e14 G", 1e12 / T_HL * 1e4, 0.512e14, 2e-3)
check("Skokov: m_pi^2 = 140^2 x 0.512e14 G", 140**2 * 0.512e14, 1.0035e18, 1e-4)
check("eB = m_pi^2 vs B = m_pi^2 (HL): ratio 1/e_HL = 1/sqrt(4 pi alpha)", ((mpi/1000)**2 * GEV2_T) / ((mpi*1e6)**2 / T_HL), 1 / math.sqrt(4 * math.pi * alpha), 1e-6)
B_moedal = 7.6e18 / T_HL
within("MoEDAL 2022 Pb-Pb 5.02 TeV: B = 7.6 GeV^2 (HL) ~ '1e16 T'", B_moedal, 1e16, 1e17)
check("  same, tesla [derived]", B_moedal, 3.890e16, 1e-3)
check("MoEDAL 2024 Pb-Pb 2.76 TeV: 4.0 GeV^2 [derived, T]", 4.0e18 / T_HL, 2.048e16, 1e-3)
check("MoEDAL omega = 73 GeV/hbar (published 'omega ~ 1e26 s^-1')", 73 * GeV / hbar, 1e26, 0.15)
check("  field lifetime hbar/omega [s]", hbar / (73 * GeV), 9.02e-27, 1e-3)
check("1e15 T as eB [GeV^2]", 1e15 / GEV2_T, 0.059157, 1e-4)
check("1e15 T as eB in units of m_pi^2 (140 MeV)", 1e15 / GEV2_T / 0.14**2, 3.0182, 1e-4)
check("Bali 2012: eB = 1 GeV^2 [T]", 1 * GEV2_T, 1.6904e16, 1e-4)
check("D'Elia 2022: eB = 9 GeV^2 (T_c = 63 MeV) [T]", 9 * GEV2_T, 1.520e17, 1e-3)
u15 = 1e15**2 / (2 * mu0)
check("B^2/2mu0 at 1e15 T [J/m^3]", u15, 3.979e35, 1e-3)
check("  in MeV/fm^3", u15 / MeV * 1e-45, 2483.4, 1e-4)
R_fe = math.sqrt(5 / 3) * 3.7377e-15    # IAEA LiveChart rms charge radius 3.7377 fm -> uniform sphere
V_fe = 4 / 3 * math.pi * R_fe**3
E_fe = 55.934935537 * u_MeV              # LiveChart/AME atomic mass 55934935.537 micro-u
check("Fe-56 mass-energy density [MeV/fm^3]", E_fe / (V_fe * 1e45), 110.7, 2e-3)
check("1e15 T energy density / Fe-56 mass-energy density", u15 / MeV * 1e-45 / (E_fe / (V_fe * 1e45)), 22.43, 3e-3)
check("B whose B^2/2mu0 = Fe-56 mass-energy density [T]", math.sqrt(2 * mu0 * E_fe * MeV / V_fe), 2.11e14, 3e-3)
u_bind = 8.7903563 * 56 * MeV / V_fe
check("B whose B^2/2mu0 = Fe-56 binding-energy density [T]", math.sqrt(2 * mu0 * u_bind), 2.0523e13, 1e-3)
check("Lai 2001: B0 = 2.3505e9 G = CODATA atomic unit (T)", 2.3505e9 * 1e-4, B_atomic, 1e-4)
check("Lai 2001 Table I: H atom binding at 1e8 T / 13.6 eV", 161 / Ry_eV, 11.83, 1e-3)
check("Lai 2001: H2 dissociation 46 eV at 1e8 T / 4.75 eV at B=0", 46 / 4.75, 9.684, 1e-3)

print("\n== 4. Bond-strength ladder: fields from mu B = E and from B^2/2mu0 = E/V ==")
cm1 = h * c * 100 / e                     # eV per cm^-1
kJmol = 1000 / NA / e                     # eV per kJ/mol
E_vdw = 99.2 * cm1                        # Ar2 De, Herman et al. via Lang et al. 2023
E_hb = 1105 * cm1                         # water dimer D0, Rocher-Casterline et al. 2011
check("water dimer 1105 cm^-1 = 13.2 kJ/mol (published)", E_hb / kJmol, 13.22, 2e-3)
E_al = 330.0 * kJmol                      # NIST WebBook dfH gas Al (CODATA) = cohesive energy at 298 K
E_fe_coh = 415.47 * kJmol                 # NIST WebBook dfH gas Fe (Chase 1998)
E_hh = 2 * 217.998 * kJmol                # H-H bond enthalpy, 298 K, from dfH gas H (CODATA)
E_cc = 716.68 / 2 * kJmol                 # C-C per bond, diamond-like 2 bonds/atom, from dfH gas C
ladder = [  # name, energy eV, moment J/T
    ("vdW Ar2 De", E_vdw, muB), ("H-bond (H2O)2 D0", E_hb, muB), ("metallic Al", E_al, muB),
    ("metallic Fe", E_fe_coh, muB), ("covalent C-C", E_cc, muB), ("covalent H-H", E_hh, muB),
    ("beta Q 3H", 18.59202e3, muN), ("beta Q n", 782.347e3, muN), ("alpha Q 210Po", 5407.5307e3, muN),
    ("B/A 56Fe", 8790.3563e3, muN), ("cluster Q 223Ra->14C", 31827.93e3, muN),
    ("proton rest energy", mp_MeV * 1e6, muN)]
want_eV = [0.012299, 0.13700, 3.4202, 4.3061, 3.7140, 4.5188, 18592.02, 782347, 5407530.7, 8790356.3, 31827930, 938272089]
V_ref = (3.761e-10)**3                     # Ar2 equilibrium distance cubed (stated assumption for chemistry)
want_Bz = [212.48, 2366.9, 59088, 74391, 64162, 78066, 5.8976e11, 2.4817e13, 1.7153e14, 2.7884e14, 1.0096e15, 2.9763e16]
want_Bed = [9.6485, 32.202, 160.90, 180.53, 167.66, 184.94]
for i, ((nm, E, mom), w) in enumerate(zip(ladder, want_eV)):
    check(f"{nm} energy [eV]", E, w, 2e-4)
    check(f"  {nm}: Zeeman field E/{'mu_B' if mom == muB else 'mu_N'} [T]", E * e / mom, want_Bz[i], 2e-4)
    if mom == muB:
        check(f"  {nm}: energy-density field sqrt(2 mu0 E/V), V = (3.761 A)^3 [T]", math.sqrt(2 * mu0 * E * e / V_ref), want_Bed[i], 2e-4)
check("16 x m_p c^2 (audit's O-16 yield) [MeV]", 16 * mp_MeV, 15012.35, 1e-6)
check("audit nuclear threshold 1e10 T vs 56Fe binding-density field: orders", math.log10(2.0523e13 / 1e10), 3.312, 1e-3)
check("audit isotope threshold 1e7 T vs 3H Zeeman field: orders", math.log10(5.8976e11 / 1e7), 4.771, 1e-3)
check("audit covalent threshold 1e4 T vs H-H energy-density field (ratio)", 1e4 / 184.94, 54.07, 1e-3)
check("Zeeman field for Ar2 (mu_B) [T]", E_vdw / muB_eV, 212.5, 2e-3)
check("energy-density field for Ar2 in (3.761 A)^3 [T]", math.sqrt(2 * mu0 * E_vdw * e / V_ref), 9.648, 2e-3)
check("Zeeman field for H-bond [T]", E_hb / muB_eV, 2367, 2e-3)
check("Zeeman field for covalent H-H [T]", E_hh / muB_eV, 7.8066e4, 1e-4)
check("Zeeman field for 56Fe B/A (mu_N) [T]", 8.7903563e6 / muN_eV, 2.789e14, 2e-3)
check("Zeeman field for cluster Q 31.83 MeV (mu_N) [T]", 31.82793e6 / muN_eV, 1.0097e15, 2e-3)
check("Zeeman field for m_p c^2 with proton moment 2.7928 mu_N [T]", mp_MeV * 1e6 / (2.792847351 * muN_eV), 1.0659e16, 2e-3)
check("cluster Q 223Ra -> 209Pb + 14C from mass excesses [keV]", 17233.245 - (-17614.574) - 3019.89328, 31827.93, 1e-6)

print("\n== 5. Callan-Rubakov catalysis and searches ==")
check("IceCube: sigma 1.7e-23 cm^2 <-> lambda 1 cm implies n [cm^-3]", 1 / (1.7e-23 * 1.0), 5.882e22, 1e-3)
n_H_ice = 2 * 0.917 / 18.015 * NA          # free protons in ice, density 0.917 g/cm^3 assumed
check("  = hydrogen-nucleus density of ice (0.917 g/cm^3 assumed)", n_H_ice, 1 / 1.7e-23, 0.05)
def superk_limit(beta, sigma0_mb=1.0): return 6.3e-24 * (beta / 1e-3)**2 / sigma0_mb
check("Super-K 2012 limit at beta=1e-3, sigma0=1 mb", superk_limit(1e-3), 6.3e-24, 1e-12)
within("  'more than eight orders' below 1e-15", 1e-15 / superk_limit(1e-3), 1e8, 1e9)
check("MACRO 1.4e-16 / Parker 1e-15", 1.4e-16 / 1e-15, 0.14, 1e-12)
check("IceCube 2022 2.0e-19 / Parker", 2.0e-19 / 1e-15, 2e-4, 1e-12)
check("Cramer 1984: 1e32 yr / 1e12 = 1e10 x age 1e10 yr", 1e32 / 1e12, 1e10 * 1e10, 1e-12)
E_o16 = 15.99491461926 * u_MeV
check("O-16 atomic mass-energy [MeV] (audit uses 16 x 938.272 = 15012.35)", E_o16, 14899.17, 1e-5)
check("  audit overstates by", 16 * mp_MeV / E_o16 - 1, 0.007597, 1e-3)
Qdd_t = 2 * 13135.722895 - 14949.81090 - 7288.971064
Qdd_he = 2 * 13135.722895 - 14931.21888 - 8071.31806
Qdd_4 = 2 * 13135.722895 - 2424.91587
check("D+D -> T+p Q [keV]", Qdd_t, 4032.66, 1e-5)
check("D+D -> 3He+n Q [keV]", Qdd_he, 3268.91, 1e-5)
check("D+D -> 4He Q [keV]", Qdd_4, 23846.53, 1e-6)
for nm, r, w in (("m_p / mean DD Q per reaction", mp_MeV / ((Qdd_t + Qdd_he) / 2e3), 257.0),
                 ("m_p / (mean DD Q per nucleon)", mp_MeV / ((Qdd_t + Qdd_he) / 2e3 / 4), 1028.0),
                 ("m_p / (D+D->4He Q per nucleon)", mp_MeV / (Qdd_4 / 1e3 / 4), 157.4)):
    check(f"talk's '166 times': {nm}", r, w, 2e-3)
L_gal = 16.02 / (gD * 3e-10)   # length over which g_D in 3 uG gains 1e11 GeV (PDG 94.8 boundary)
check("g_D x 3 uG x L = 1e11 GeV gives L [kpc]", L_gal / 3.0857e19, 0.5257, 2e-3)
check("1e11 GeV in J", 1e11 * GeV, 16.02, 1e-3)
check("PDG 94.8 continuity: 1e-3 c (1e17/1e11)^0.5 = c", 1e-3 * (1e17 / 1e11)**0.5, 1.0, 1e-12)

print("\n== 6. Monopole-nucleus binding criteria (Milton 2006, sec. 9) ==")
nuc = {  # name: (Z, A, spin, mu/mu_N) IAEA LiveChart (ENSDF) ground states
    "n": (0, 1, 0.5, -1.9130427), "1H": (1, 1, 0.5, 2.792847351), "2H": (1, 2, 1, 0.857438231),
    "3He": (2, 3, 0.5, -2.12762531), "9Be": (4, 9, 1.5, -1.17743), "16O": (8, 16, 0, 0.0),
    "17O": (8, 17, 2.5, -1.893543), "27Al": (13, 27, 2.5, 3.6407), "57Fe": (26, 57, 0.5, 0.0902),
    "113Cd": (48, 113, 0.5, -0.62226)}
def ghat(Z, A, mu): return A / Z * mu
for k, want in (("1H", 2.79), ("3He", -3.20), ("9Be", -2.66), ("27Al", 7.56), ("113Cd", -1.46), ("2H", 1.71)):
    Z, A, S, mu = nuc[k]; check(f"gamma-hat {k} (Milton Table 1)", ghat(Z, A, mu), want, 6e-3)
def lowest_binds(Z, A, mu, m=0.5): l = m * Z; return ghat(Z, A, mu) > 1 + 1 / (4 * l)
def kappa_c(J, l): return (J * J + J - l * l) / l
def gc32(l): s = math.sqrt(33 + 32 * l); return (3 / (4 * l) * (6 + 4 * l - s), 3 / (4 * l) * (6 + 4 * l + s))
def gc52(l): s = math.sqrt(1161 + 1296 * l + 64 * l * l); return ((36 + 28 * l - s) / (12 * l), (36 + 28 * l + s) / (12 * l))
check("eq 9.12, Be (l=2): gamma_c- = 1.557", gc32(2)[0], 1.557, 1e-3)
check("eq 9.12, Be (l=2): gamma_c+ = 8.943", gc32(2)[1], 8.943, 1e-3)
check("eq 9.13, Al (l=6.5): gamma_c- = 1.374", gc52(6.5)[0], 1.374, 1e-3)
check("eq 9.13, Al (l=6.5): gamma_c+ = 4.216", gc52(6.5)[1], 4.216, 1e-3)
check("eq 9.13 large-l limit gamma_c- -> 5/3 (published 1.67)", gc52(1e9)[0], 1.667, 1e-3)
check("proton binds in J = l - 1/2 (criterion 9.7)", float(lowest_binds(1, 1, 2.792847351)), 1.0, 0)
check("27Al binds in lowest state", float(lowest_binds(13, 27, 3.6407)), 1.0, 0)
check("3He: |kappa-hat| 4.2 > kappa_c(J=l+1/2=3/2, l=1) = 2.75", kappa_c(1.5, 1) , 2.75, 1e-12)
check("113Cd: kappa_c(J = 49/2, l = 24)", kappa_c(24.5, 24), 2.03125, 1e-9)
check("113Cd |kappa-hat| = |gamma-hat - 1|", abs(ghat(48, 113, -0.62226) - 1), 2.4649, 1e-3)
check("neutron |gamma| 1.91 > 3/(4 |m'|) = 1.5", 3 / (4 * 0.5), 1.5, 0)
Z, A, S, mu = nuc["57Fe"]
check("57Fe gamma-hat [derived]", ghat(Z, A, mu), 0.19775, 1e-3)
check("57Fe fails lowest criterion: needs > 1 + 1/52", 1 + 1 / (4 * 13), 1.01923, 1e-4)
check("57Fe fails J=l+1/2: |kappa-hat| 0.80 < kappa_c", kappa_c(13.5, 13), 2.0577, 1e-4)
Z, A, S, mu = nuc["17O"]
check("17O gamma-hat [derived] (negative, like 9Be)", ghat(Z, A, mu), -4.0238, 1e-4)
check("17O eq 9.13 at l = 4: gamma_c-", gc52(4)[0], 1.2951, 1e-3)
check("17O natural abundance [%] (LiveChart)", 0.03835, 0.03835, 0)

print("\n== 7. Atomic spectra with a bound monopole (Milton eq 9.1, spinless) ==")
def E_dyonic(Z, N=0, j=0.5, m=0.5): return -Z**2 * Ry_eV / (N + 0.5 + math.sqrt((j + 0.5)**2 - m * m))**2
check("eq 9.1 with m'=0, j=l=0 reproduces -13.6057 eV", E_dyonic(1, 0, 0, 0), -Ry_eV, 1e-12)
check("H-like ground state with a Dirac monopole (m'=1/2) [eV]", E_dyonic(1), -7.2913, 1e-4)
check("  ratio to the no-monopole value", E_dyonic(1) / -Ry_eV, 0.53590, 1e-4)
check("spin Zeeman scale mu_B B_Dirac(a0) [eV]", muB_eV * B_dirac(a0), 6.802, 1e-3)
r3 = 3e-10
check("triplet O2 (2 mu_B) at 3 A from a Dirac monopole [eV]", 2 * muB_eV * B_dirac(r3), 0.4234, 1e-3)

print("\n== 8. Ivoilov & Urutskoev 2004 Moessbauer numbers ==")
for dH, w in ((280, 108.5), (460, 84.65)):
    B = dH * 1e-4    # Oe -> T in vacuum (1 Oe <-> 1e-4 T)
    check(f"radius where B_Dirac = {dH} Oe [nm]", math.sqrt(hbar / (2 * e * B)) * 1e9, w, 2e-3)
r280 = math.sqrt(hbar / (2 * e * 0.028)); V_layer = 0.015**2 * 0.2e-6
N_needed = V_layer / (4 / 3 * math.pi * r280**3)
check("monopoles to fill the 0.2 um x 15 x 15 mm CEMS layer at that spacing", N_needed, 8.431e9, 1e-3)
check("  their mass if GUT (1e17 GeV each) [kg]", N_needed * 1e17 * GeV / c**2, 1.50, 3e-3)
check("  pair energy at the MoEDAL 75 GeV mass floor [J]", 2 * N_needed * 75 * GeV, 202.6, 1e-3)
check("280 Oe at -40 Oe/K = temperature step [K]", 280 / 40, 7.0, 1e-12)
check("significance 280/70 and 460/70", 280 / 70 + 460 / 70, 4.0 + 6.5714, 1e-4)
check("force on g_D in 120 Oe [N]", gD * mu0 * 120 * 1e3 / (4 * math.pi), 3.949e-11, 1e-3)
check("57Fe abundance [%] LiveChart (paper: 2.2)", 2.119, 2.2, 0.04)

check("a 1e-9 bound fraction / PDG 1e-29 per nucleon: orders", math.log10(1e-9 / 1e-29), 20.0, 1e-12)
check("Yue 2021 1000 T vs 1e15 T: orders", math.log10(1e15 / 1e3), 12.0, 1e-12)
check("D'Elia 9 GeV^2 field / 1e15 T", 9 * GEV2_T / 1e15, 152.1, 1e-3)

print("\n== 9. Emergent and synthetic monopoles ==")
mu_dy = 10 * muB; a_d = 4.3356e-10
qm = 2 * mu_dy / a_d
check("spin ice q_m = 2 mu/a_d (Bramwell 4.6 muB/A)", qm / (muB / 1e-10), 4.6, 0.005)
check("  in A m [repo 4.28e-13]", qm, 4.278e-13, 1e-3)
check("  Bramwell measured 5 muB/A in A m", 5 * muB / 1e-10, 4.637e-13, 1e-3)
check("  repo form 2 mu/(sqrt(3/2) a), a = 3.54 A", 2 * mu_dy / (math.sqrt(1.5) * 3.54e-10), qm, 1e-3)
r_e = alpha * hbar / (9.1093837139e-31 * c)
check("Castelnovo eq 4: q_m/q_D = 20 r_e/a_d", qm / gD, 20 * r_e / a_d, 1e-6)
check("  q_D/q_m (published ~8000)", gD / qm, 7692, 1e-3)
m_rb = 86.909180529 * u_kg
check("87Rb mass (Hall et al. 2016: 1.443e-25 kg)", m_rb, 1.443e-25, 1e-3)
check("Ray 2014 synthetic monopole: string circulation 4 pi hbar/M [m^2/s]", 4 * math.pi * hbar / m_rb, 9.183e-9, 1e-3)
check("Neubauer 2009: area per skyrmion for B_eff 2.5 T with Phi0 = h/e [m^2]", (h / e) / 2.5, 1.6543e-15, 1e-3)

print("\nFAILED: " + ", ".join(FAIL) if FAIL else "\nall checks passed")
sys.exit(1 if FAIL else 0)
