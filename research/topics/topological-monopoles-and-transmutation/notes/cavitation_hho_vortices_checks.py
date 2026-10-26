"""Checks for cavitation_hho_vortices.md: bubble collapse, D-D fusion at sonoluminescence
temperatures, oxyhydrogen thermochemistry, Ohmasa and Kladov patent numbers, vortex
pressures, dust-devil charge, Helmholtz and plate acoustics.

Standard library only. Run: python3 -I cavitation_hho_vortices_checks.py
Every computed number is compared with a published value (source in the comment) or, where
no published value exists, with the derivation it rests on; inferences are asserted as the
inequality the notes state. Exit status is non-zero on any failure.
"""
import math
import sys

FAIL = []


def check(name, got, want, rtol):
    ok = abs(got / want - 1.0) <= rtol
    print(f"{'ok  ' if ok else 'FAIL'} {name}: {got:.6g} vs {want:.6g} (rtol {rtol:g})")
    if not ok:
        FAIL.append(name)


def holds(name, cond, detail):
    print(f"{'ok  ' if cond else 'FAIL'} {name}: {detail}")
    if not cond:
        FAIL.append(name)


# CODATA 2018 exact / recommended constants
K_B = 1.380649e-23
E = 1.602176634e-19
N_A = 6.02214076e23
R = K_B * N_A
F = E * N_A
EPS0 = 8.8541878128e-12
AMU = 1.66053906660e-27
ALPHA = 7.2973525693e-3
H_C_CM = 1.98644586e-23  # h c in J cm (1 cm^-1 = 1.986e-23 J)

# ---------------------------------------------------------------- 1. Rayleigh collapse
# t_c = tau R0 sqrt(rho/dp), tau = sqrt(3/2) * B(5/6, 1/2) / 3.
# Klotz, Phys. Fluids 25, 082109 (2013), arXiv:1206.5753: "0.91468".
beta = math.gamma(5 / 6) * math.gamma(0.5) / math.gamma(4 / 3)
tau = math.sqrt(1.5) * beta / 3
check("Rayleigh collapse constant", tau, 0.91468, 1e-5)


def rayleigh_plesset_collapse(r0=1.0, rho=1.0, dp=1.0, steps=200000):
    """Empty cavity, inviscid, no surface tension: R R'' + 1.5 R'^2 = -dp/rho.
    Integrate in the energy form R'^2 = (2dp/3rho)((r0/R)^3 - 1) over R."""
    # t = int dR / |R'| from 0 to r0; substitute R = r0 (1 - s^2) to remove the r0 singularity
    total = 0.0
    n = steps
    for i in range(n):
        s = (i + 0.5) / n
        r = r0 * (1 - s * s)
        drds = 2 * r0 * s
        v = math.sqrt(2 * dp / (3 * rho) * ((r0 / r) ** 3 - 1))
        total += drds / v / n
    return total


check("Rayleigh-Plesset integration = tau", rayleigh_plesset_collapse(), 0.91468, 1e-3)

# Plesset & Chapman, JFM 47, 283 (1971): jet ~130 m/s (bubble touching wall) and ~170 m/s
# (half a radius off) for dp/rho = 1e6 cm^2/s^2 ~ 1 atm / water density.
check("Plesset-Chapman scaling: 1 atm/rho_water vs 1e6 cm^2/s^2", 101325 / 1000.0 / 100.0, 1.0, 0.02)
k_touch, k_half = 130 / math.sqrt(100), 170 / math.sqrt(100)
check("jet coefficient (touching)", k_touch, 13.0, 1e-9)
check("jet coefficient (gamma=1.5)", k_half, 17.0, 1e-9)

# Water hammer p = rho c v. Philipp & Lauterborn, JFM 361, 75 (1998): 83 m/s -> "about 0.1 GPa".
p_wh = 998.0 * 1483.0 * 83.0
check("water-hammer pressure at 83 m/s (GPa)", p_wh / 1e9, 0.1, 0.30)

# Huang et al., Sci. Rep. 14, 214 (2024): ~21 bar water, pressure spikes >720 bar.
# Inference: a vapour void collapsing under 21 bar reaches wall speeds whose water hammer exceeds
# 720 bar once R0/R ~ 2 (Rayleigh energy form), so condensation-induced water hammer suffices.
dp_huang, rho_w, c_w = 21e5, 1000.0, 1480.0
u_wall = math.sqrt(2 * dp_huang / (3 * rho_w) * (2 ** 3 - 1))
holds("Huang: 21-bar collapse water hammer exceeds 720 bar at R0/R=2",
      rho_w * c_w * u_wall > 720e5, f"U={u_wall:.0f} m/s, p_wh={rho_w*c_w*u_wall/1e5:.0f} bar")
dv_720 = 720e5 / (rho_w * c_w)
holds("Huang: 720 bar needs only ~49 m/s of arrested flow", 40 < dv_720 < 60, f"dv={dv_720:.1f} m/s")

# ---------------------------------------------------------------- 2. Sonoluminescence temperatures
# Flannigan & Suslick, Nature 434, 52 (2005): Ar emission temperatures 4,000-15,000 K;
# Ar emitting states >13 eV "cannot be thermally populated".
kt_ev = K_B * 15000 / E
check("kT at 15,000 K (eV)", kt_ev, 1.2926, 1e-3)
boltz = math.exp(-13.0 / kt_ev)
holds("13 eV Boltzmann factor at 15,000 K is small", boltz < 1e-4, f"exp(-13/kT) = {boltz:.2e}")
# Ideal adiabatic compression of Ar (gamma=5/3): T = T0 (R0/R)^(3(gamma-1)); ratio for 15,000 K
g = 5 / 3
ratio = (15000 / 300) ** (1 / (3 * (g - 1)))
check("Ar compression ratio for 15,000 K from 300 K", ratio, math.sqrt(50), 1e-12)
check("adiabatic T round trip", 300 * ratio ** (3 * (g - 1)), 15000, 1e-12)

# ---------------------------------------------------------------- 3. D-D fusion (port of plasmoid::fusion)
M_D = 1_875_612.9424  # keV, as crates/plasmoid/src/fusion.rs
MR_DD = M_D / 2


B_G_DD = 31.3970


def bh_s_ddn(e):
    """Bosch & Hale (1992) S-factor Pade fit for D(d,n)3He, keV mb; e in keV (centre of mass)."""
    a = [5.3701e4, 3.3027e2, -1.2706e-1, 2.9327e-5, -2.5151e-9]
    return a[0] + e * (a[1] + e * (a[2] + e * (a[3] + e * a[4])))


def bh_sigma_ddn(e):
    """Cross section, mb: S(E) / (E exp(B_G/sqrt E))."""
    return bh_s_ddn(e) / (e * math.exp(B_G_DD / math.sqrt(e)))


def bh_reactivity_ddn(t):
    """Bosch & Hale (1992) reactivity form, m^3/s, t in keV (0.2-100)."""
    b_g, mr = 31.3970, 937_814.0
    c = [5.43360e-12, 5.85778e-3, 7.68222e-3, 0.0, -2.96400e-6, 0.0, 0.0]
    theta = t / (1 - t * (c[1] + t * (c[3] + t * c[5])) / (1 + t * (c[2] + t * (c[4] + t * c[6]))))
    xi = (b_g ** 2 / (4 * theta)) ** (1 / 3)
    return 1e-6 * c[0] * theta * math.sqrt(xi / (mr * t ** 3)) * math.exp(-3 * xi)


def maxwell_reactivity(s_kev_mb, kt, mr_kev, e0, width, b_g=B_G_DD):
    """<sigma v> = sqrt(8/(pi m)) (kT)^-3/2 int S(E) exp(-E/kT - B_G/sqrt E) dE, Simpson's rule
    about the Gamow peak (sigma E = S e^(-B_G/sqrt E) keeps the exponent in one place). m^3/s."""
    m = mr_kev * 1e3 * E / (299792458.0 ** 2)  # kg
    lo, hi = max(e0 - 8 * width, e0 / 10), e0 + 8 * width
    n = 4000
    h = (hi - lo) / n
    acc = 0.0
    for i in range(n + 1):
        en = lo + i * h
        w = 1 if i in (0, n) else (4 if i % 2 else 2)
        acc += w * s_kev_mb(en) * math.exp(-en / kt - b_g / math.sqrt(en))
    integral = acc * h / 3 * (1e-31 * 1e3 * E) * (1e3 * E)  # (keV mb -> J m^2) x (dE keV -> J)
    return math.sqrt(8 / (math.pi * m)) * (kt * 1e3 * E) ** -1.5 * integral


def gamow_peak(e_g, kt):
    e0 = (e_g * kt * kt / 4) ** (1 / 3)
    return e0, 4 * math.sqrt(e0 * kt / 3)


# Repository test value (crates/plasmoid/src/fusion.rs): sigma_DDn(50 keV) = 16.487 mb
check("Bosch-Hale sigma D(d,n) at 50 keV (mb)", bh_sigma_ddn(50.0), 16.487, 1e-3)
e_g = 2 * MR_DD * (math.pi * ALPHA) ** 2
check("D-D Gamow energy (keV) = B_G^2", e_g, 31.3970 ** 2, 2e-4)
e0, wid = gamow_peak(e_g, 1.0)
sv_num = maxwell_reactivity(bh_s_ddn, 1.0, MR_DD, e0, wid)
check("Maxwellian <sv> D(d,n) at 1 keV: integral vs Bosch-Hale fit", sv_num, bh_reactivity_ddn(1.0), 0.03)

kt_sl = K_B * 15000 / E / 1e3  # keV
e0_sl, wid_sl = gamow_peak(e_g, kt_sl)
sv_sl = maxwell_reactivity(bh_s_ddn, kt_sl, MR_DD, e0_sl, wid_sl)
tau_g = 3 * e0_sl / kt_sl
print(f"     D-D at 15,000 K: Gamow peak {e0_sl*1e3:.1f} eV, exponent {tau_g:.1f}, <sv> = {sv_sl:.2e} m^3/s")
check("Gamow exponent 3E0/kT = 3(E_G/4kT)^(1/3) at 15,000 K", tau_g, 3 * (e_g / (4 * kt_sl)) ** (1 / 3), 1e-12)
holds("Gamow peak at 15,000 K is ~74 eV, exponent ~173", 73 < e0_sl * 1e3 < 76 and 170 < tau_g < 175, f"{e0_sl*1e3:.1f} eV, {tau_g:.1f}")
holds("D(d,n) reactivity at 15,000 K is below 1e-90 m^3/s", sv_sl < 1e-90, f"{sv_sl:.2e}")
# Generous hot spot: n = 1e29 m^-3 (above liquid D density), r = 1 um, 1 ns -> fusions per collapse
n_d, vol, dt = 1e29, 4 / 3 * math.pi * (1e-6) ** 3, 1e-9
per_collapse = 0.5 * n_d ** 2 * sv_sl * vol * dt
holds("fusions per collapse at 15,000 K < 1e-60", per_collapse < 1e-60, f"{per_collapse:.1e}")
# Taleyarkhan et al., Science 295, 1868 (2002) abstract: simulated 1e6-1e7 K. At 1e7 K (0.862 keV):
sv_1e7 = bh_reactivity_ddn(K_B * 1e7 / E / 1e3)
print(f"     D(d,n) <sv> at 1e7 K = {sv_1e7:.2e} m^3/s; same hot spot -> {0.5*n_d**2*sv_1e7*vol*dt:.2e} per collapse")
holds("ratio of reactivities 1e7 K vs 15,000 K exceeds 1e64", sv_1e7 / sv_sl > 1e64, f"{sv_1e7/sv_sl:.1e}")
holds("the claimed 1e7 K hot spot would give >1 fusion per collapse (so the claim is testable)",
      0.5 * n_d ** 2 * sv_1e7 * vol * dt > 1, f"{0.5*n_d**2*sv_1e7*vol*dt:.0f} per collapse")

# ---------------------------------------------------------------- 4. Electrolysis and oxyhydrogen
# NIST WebBook (CODATA, Cox, Wagman & Medvedev): dfH(H2O,l) = -285.830 kJ/mol, dfH(H2O,g) = -241.826,
# S(H2O,l) = 69.95, S(H2) = 130.680, S(O2) = 205.152 J/(mol K).
dh_l, dh_g = -285.830e3, -241.826e3
ds = 69.95 - 130.680 - 0.5 * 205.152
dg = dh_l - 298.15 * ds
check("reversible cell voltage (V), textbook 1.229", -dg / (2 * F), 1.229, 5e-4)
check("thermoneutral voltage (V), textbook 1.481", -dh_l / (2 * F), 1.481, 5e-4)
check("H2 higher heating value (MJ/kg), from NIST dfH", -dh_l / 2.01588e-3 / 1e6, 141.79, 1e-3)
check("H2 lower heating value (MJ/kg), from NIST dfH", -dh_g / 2.01588e-3 / 1e6, 119.96, 1e-3)

# Ideal-gas molar volume at 0 C, 1 atm; CODATA 2018: 22.413969 L/mol
vm0 = R * 273.15 / 101325 * 1e3
check("molar volume at STP (L/mol)", vm0, 22.413969, 1e-6)
# "Implosion": 2 H2 + O2 (3 mol gas) -> 2 H2O (liquid, 36.03 g ~ 36 cm^3)
shrink = 3 * vm0 / (2 * 18.015 / 1.0 / 1e3)
holds("volume shrinks ~1,860x when the burnt gas condenses", 1800 < shrink < 1900, f"{shrink:.0f}x")

# Ohmasa, US 7,459,071 B2, first embodiment: 2.0 V, 5 A/dm^2, 50 anodes, 1,000 L/h of mixed gas;
# cell 270 x 1660 x 390 mm (H). Assume litres at 0 C, 1 atm (the patent does not say).
mol_h2_per_h = 1000 * (2 / 3) / vm0
i_needed = mol_h2_per_h * 2 * F / 3600
area_dm2 = i_needed / 5.0
plate_max_dm2 = 2.7 * 3.9  # one face of a plate filling the cell cross-section
holds("Ohmasa: Faraday current for 1,000 L/h fits on 50 anodes at 5 A/dm^2",
      area_dm2 < 50 * plate_max_dm2, f"I = {i_needed:.0f} A, {area_dm2:.0f} dm^2 of {50*plate_max_dm2:.0f} available")
check("Ohmasa: electrical input at 2.0 V (kW)", 2.0 * i_needed / 1e3, 3.19, 3e-3)
eta_hhv = (-dh_l / (2 * F)) / 2.0
check("Ohmasa: HHV efficiency at 2.0 V (100% Faradaic)", eta_hhv, 1.481 / 2.0, 5e-4)
holds("Ohmasa: gas energy < electrical input (COP < 1)", eta_hhv < 1, f"{eta_hhv:.3f}")

# Fuel-cell bound: electrical work per mol H2 <= -dG = 2 F E_rev.
check("max fuel-cell work per mol H2 (kJ)", -dg / 1e3, 237.14, 1e-3)

# Ohmasa, US 2016/0145521 A1: "44.8 MJ/m3 (10,700 Kcal/m3)" with CO2 20% / OHMASA-GAS 80%,
# "43.9 MJ/m3 (10,500 Kcal/m3)" with 40/60.
check("kcal->MJ 10,700", 10700 * 4.1868e-3, 44.8, 1e-3)
check("kcal->MJ 10,500", 10500 * 4.1868e-3, 43.9, 2e-3)
hhv_h2_m3 = -dh_l / vm0 * 1e3 / 1e6  # MJ per m3 H2 at STP
check("H2 HHV per m^3 at 0 C, 1 atm (MJ)", hhv_h2_m3, 12.75, 1e-3)
for co2, claim in ((0.2, 44.8), (0.4, 43.9)):
    for x_h2, label in ((2 / 3, "2:1 gas"), (0.70, "70% H2, CA 2,522,593 upper bound")):
        bound = (1 - co2) * x_h2 * hhv_h2_m3
        holds(f"CO2 {co2:.0%}: claimed {claim} MJ/m3 exceeds chemical bound ({label})",
              claim / bound > 4, f"bound {bound:.2f} MJ/m3, claim/bound = {claim/bound:.1f}")

# ---------------------------------------------------------------- 5. Adiabatic flame temperature
# GRI-Mech 3.0 NASA-7 polynomials (thermo30.dat), [high a1..a7, low a1..a7], Tmid 1000 K.
NASA = {
    "O":   ([2.56942078E+00, -8.59741137E-05, 4.19484589E-08, -1.00177799E-11, 1.22833691E-15, 2.92175791E+04, 4.78433864E+00],
            [3.16826710E+00, -3.27931884E-03, 6.64306396E-06, -6.12806624E-09, 2.11265971E-12, 2.91222592E+04, 2.05193346E+00]),
    "O2":  ([3.28253784E+00, 1.48308754E-03, -7.57966669E-07, 2.09470555E-10, -2.16717794E-14, -1.08845772E+03, 5.45323129E+00],
            [3.78245636E+00, -2.99673416E-03, 9.84730201E-06, -9.68129509E-09, 3.24372837E-12, -1.06394356E+03, 3.65767573E+00]),
    "H":   ([2.50000001E+00, -2.30842973E-11, 1.61561948E-14, -4.73515235E-18, 4.98197357E-22, 2.54736599E+04, -4.46682914E-01],
            [2.50000000E+00, 7.05332819E-13, -1.99591964E-15, 2.30081632E-18, -9.27732332E-22, 2.54736599E+04, -4.46682853E-01]),
    "H2":  ([3.33727920E+00, -4.94024731E-05, 4.99456778E-07, -1.79566394E-10, 2.00255376E-14, -9.50158922E+02, -3.20502331E+00],
            [2.34433112E+00, 7.98052075E-03, -1.94781510E-05, 2.01572094E-08, -7.37611761E-12, -9.17935173E+02, 6.83010238E-01]),
    "OH":  ([3.09288767E+00, 5.48429716E-04, 1.26505228E-07, -8.79461556E-11, 1.17412376E-14, 3.85865700E+03, 4.47669610E+00],
            [3.99201543E+00, -2.40131752E-03, 4.61793841E-06, -3.88113333E-09, 1.36411470E-12, 3.61508056E+03, -1.03925458E-01]),
    "H2O": ([3.03399249E+00, 2.17691804E-03, -1.64072518E-07, -9.70419870E-11, 1.68200992E-14, -3.00042971E+04, 4.96677010E+00],
            [4.19864056E+00, -2.03643410E-03, 6.52040211E-06, -5.48797062E-09, 1.77197817E-12, -3.02937267E+04, -8.49032208E-01]),
}
ATOMS = {"O": (0, 1), "O2": (0, 2), "H": (1, 0), "H2": (2, 0), "OH": (1, 1), "H2O": (2, 1)}  # (H, O)
R_GRI = 8.314510  # Marzouk (2023) value


def coeffs(sp, t):
    hi, lo = NASA[sp]
    return hi if t >= 1000 else lo


def h_rt(sp, t):
    a = coeffs(sp, t)
    return a[0] + a[1] * t / 2 + a[2] * t * t / 3 + a[3] * t ** 3 / 4 + a[4] * t ** 4 / 5 + a[5] / t


def s_r(sp, t):
    a = coeffs(sp, t)
    return a[0] * math.log(t) + a[1] * t + a[2] * t * t / 2 + a[3] * t ** 3 / 3 + a[4] * t ** 4 / 4 + a[6]


def h(sp, t):
    return h_rt(sp, t) * R_GRI * t


check("GRI-Mech dfH(H2O,g) at 298.15 K vs NIST -241.826 kJ/mol",
      (h("H2O", 298.15) - h("H2", 298.15) - 0.5 * h("O2", 298.15)) / 1e3, -241.826, 1e-3)


def equilibrium(t, p_atm=1.0):
    """Mole numbers per (2 H2 + O2) at t, p: species H2O H2 O2 OH H O, from formation
    constants relative to H2 and O2 (1 bar standard state)."""
    p = p_atm * 1.01325  # bar
    lnk = {}
    for sp, (nh, no) in ATOMS.items():
        g = h_rt(sp, t) - s_r(sp, t)
        g_ref = nh / 2 * (h_rt("H2", t) - s_r("H2", t)) + no / 2 * (h_rt("O2", t) - s_r("O2", t))
        lnk[sp] = -(g - g_ref)

    def parts(x, y):  # x = ln p_H2, y = ln p_O2
        return {sp: math.exp(lnk[sp] + nh / 2 * x + no / 2 * y) for sp, (nh, no) in ATOMS.items()}

    x, y = math.log(0.1 * p), math.log(0.05 * p)
    for _ in range(200):
        def resid(x, y):
            pp = parts(x, y)
            tot = sum(pp.values())
            hb = sum(ATOMS[s][0] * v for s, v in pp.items())
            ob = sum(ATOMS[s][1] * v for s, v in pp.items())
            return math.log(tot / p), math.log(hb / (2 * ob))
        f1, f2 = resid(x, y)
        if abs(f1) + abs(f2) < 1e-13:
            break
        d = 1e-7
        a11 = (resid(x + d, y)[0] - f1) / d
        a12 = (resid(x, y + d)[0] - f1) / d
        a21 = (resid(x + d, y)[1] - f2) / d
        a22 = (resid(x, y + d)[1] - f2) / d
        det = a11 * a22 - a12 * a21
        dx, dy = (-f1 * a22 + f2 * a12) / det, (-a11 * f2 + a21 * f1) / det
        x, y = x + max(-2, min(2, dx)), y + max(-2, min(2, dy))
    pp = parts(x, y)
    ob = sum(ATOMS[s][1] * v for s, v in pp.items())
    return {s: v * 2 / ob for s, v in pp.items()}  # scale to 2 O atoms


def bisect(fn, a, b):
    fa = fn(a)
    for _ in range(100):
        m = 0.5 * (a + b)
        fm = fn(m)
        if (fm > 0) == (fa > 0):
            a, fa = m, fm
        else:
            b = m
    return 0.5 * (a + b)


h_react = 2 * h("H2", 298.15) + h("O2", 298.15)
aft_eq = bisect(lambda t: sum(n * h(s, t) for s, n in equilibrium(t).items()) - h_react, 2000, 4000)
check("oxyhydrogen equilibrium AFT vs Marzouk (2023) GRI-Mech/Cantera 3,076.92 K", aft_eq, 3076.92, 2e-3)
check("oxyhydrogen equilibrium AFT vs Marzouk (2023) CEARUN 3,074.51 K", aft_eq, 3074.51, 3e-3)
aft_cc = bisect(lambda t: 2 * h("H2O", t) - h_react, 3000, 6000)
check("oxyhydrogen complete-combustion AFT vs Marzouk GRI-Mech Excel 4,890.21 K", aft_cc, 4890.21, 3e-3)
eq = equilibrium(aft_eq)
tot = sum(eq.values())
print("     equilibrium mole fractions at AFT:", ", ".join(f"{s} {v/tot:.3f}" for s, v in eq.items()))
check("equilibrium products keep H/O = 2", sum(ATOMS[k][0] * v for k, v in eq.items()) / sum(ATOMS[k][1] * v for k, v in eq.items()), 2.0, 1e-9)
holds("about 42% of the products are dissociated at the AFT", 0.40 < 1 - eq["H2O"] / tot < 0.44, f"{1-eq['H2O']/tot:.3f}")

# ---------------------------------------------------------------- 6. Kladov RU 2054604 C1
# Rule: P1 = (0.3-0.7)(P2+P3), P2+P3-P1 = (1-10) sigma. Example 1: tap water 20 C, sigma ~0.35 MPa,
# P3 ~ 0, P1 = 0.15 MPa, P2 = 0.5 MPa; N1 = 6 kW, N2 = 6.5 kW, "coefficient (k) of 1.08".
sigma_t, p3 = 0.35, 0.0
p2 = sigma_t / (1 - 0.3)
p1 = 0.3 * (p2 + p3)
check("Kladov P2 (MPa)", p2, 0.5, 1e-9)
check("Kladov P1 (MPa)", p1, 0.15, 1e-9)
check("Kladov k = N2/N1", 6.5 / 6.0, 1.08, 5e-3)
# Kladov: collapse energy "several orders of magnitude higher than the binding energy of nucleons".
# Mechanical energy density available ~ (P1+P2) per unit liquid volume, per water molecule:
e_per_molecule_ev = (p1 + p2) * 1e6 * (18.015e-3 / 998 / N_A) / E
holds("Kladov: pressure work per water molecule is ~1e-4 eV", 1e-5 < e_per_molecule_ev < 1e-3, f"{e_per_molecule_ev:.2e} eV")
holds("Kladov: even 15,000 K hot spots are ~1e-7 of 8 MeV/nucleon", kt_ev / 8e6 < 1e-6, f"{kt_ev/8e6:.1e}")

# ---------------------------------------------------------------- 7. Ohmasa vibrator: cavitation number
# US 7,459,071: blade tip amplitude 0.1-15.0 mm at 200-1,000 per minute.
v_tip = 2 * math.pi * (1000 / 60) * 15e-3
p_v_assumed = 20e3  # Pa, above water's vapour pressure at the patent's 55 C (assumption)
sigma_c = (101325 - p_v_assumed) / (0.5 * 1000 * v_tip ** 2)
holds("Ohmasa blade: peak tip speed < 2 m/s", v_tip < 2, f"{v_tip:.2f} m/s")
holds("Ohmasa blade: cavitation number >> 1 (no hydrodynamic inception)", sigma_c > 10, f"sigma = {sigma_c:.0f}")

# ---------------------------------------------------------------- 8. Tornadoes, H-bonds, Rankine vortex
# Rocher-Casterline et al., J. Chem. Phys. 134, 211101 (2011): D0 = 1105 cm^-1 (13.2 kJ/mol).
d0 = 1105 * H_C_CM * N_A / 1e3
check("water dimer D0 (kJ/mol)", d0, 13.2, 2e-3)
# Wurman et al., BAMS 88, 31 (2007): winds up to 132 m/s.
v_t = 132.0
ke = 0.5 * 18.015e-3 * v_t ** 2 / 1e3  # kJ/mol of bulk-flow kinetic energy
holds("bulk KE of water moving at 132 m/s is ~1% of one H-bond", ke / d0 < 0.02, f"{ke:.3f} kJ/mol = {ke/d0:.3%} of D0")
holds("thermal RT at 298 K exceeds that KE >10x", R * 298.15 / 1e3 / ke > 10, f"RT/KE = {R*298.15/1e3/ke:.0f}")
pv_work = 1e4 * 18.07e-6 / 1e3  # 100 hPa on 1 mol liquid water, kJ/mol
holds("100 hPa pressure work on liquid water << H-bond", pv_work / d0 < 1e-4, f"{pv_work/d0:.1e} of D0")


def rankine_deficit(vmax, rho, rc=1.0, n=200000, rmax=2000.0):
    """Integrate cyclostrophic dp/dr = rho v^2/r for a Rankine vortex from 0 to rmax."""
    acc = 0.0
    r_prev = 0.0
    for i in range(1, n + 1):
        r = rmax * (i / n) ** 3  # dense near the core
        rm = 0.5 * (r + r_prev)
        v = vmax * rm / rc if rm < rc else vmax * rc / rm
        acc += rho * v * v / rm * (r - r_prev)
        r_prev = r
    return acc


rho_air = 1.15
check("Rankine core deficit integral = rho vmax^2", rankine_deficit(v_t, rho_air), rho_air * v_t ** 2, 1e-3)
check("Rankine deficit at 132 m/s (hPa), about twice the 100 hPa Manchester record", rho_air * v_t ** 2 / 100, 2 * 100, 0.01)
check("dynamic pressure at 132 m/s (hPa)", 0.5 * rho_air * v_t ** 2 / 100, 100.0, 0.01)
check("oxyhydrogen gas volume 3 mol at STP (L)", 3 * vm0, 67.2, 1e-3)
check("Rankine speed for a 100 hPa deficit (m/s)", math.sqrt(1e4 / rho_air), 93.0, 0.01)

# ---------------------------------------------------------------- 9. Dust-devil charge (Farrell et al. 2004)
# E = 4.3 kV/m (saturated) at r ~ 50 m, diameter ~10 m; n = 4 pi eps0 r^2 E / (e V/2) ~ 1e6 el/cc
e_dd, r_dd, d_dd = 4.3e3, 50.0, 10.0
for h_dd in (100.0, 200.0):
    vol_dd = math.pi * (d_dd / 2) ** 2 * h_dd
    n_dd = 4 * math.pi * EPS0 * r_dd ** 2 * e_dd / (E * vol_dd / 2) / 1e6
    holds(f"dust-devil charge density ~1e6 /cc (h={h_dd:.0f} m)", 2e5 < n_dd < 5e6, f"{n_dd:.2e} /cc")
holds("4.3 kV/m is 10-20% of Mars breakdown 20 kV/m (Farrell: '10-20%')", 0.1 <= 4.3 / 20 <= 0.22, f"{4.3/20:.3f}")
holds("0.8 MV total vertical potential implies a 100-200 m column", 100 <= 0.8e6 / e_dd <= 200, f"h = {0.8e6/e_dd:.0f} m")

# ---------------------------------------------------------------- 10. Helmholtz resonator
# Monteiro et al., arXiv:1805.04014: f = (c/2pi) sqrt(A/(V L')), L' = L + 1.5a; neck d = 29.0 mm,
# L = 70 mm, V = 435 cm^3; reference sound speed 344 m/s at 22.3 C; wavelengths 1-1.5 m.
c_air = math.sqrt(1.4 * R * (273.15 + 22.3) / 28.9647e-3)
check("speed of sound at 22.3 C (m/s)", c_air, 344.0, 3e-3)
a_n = 29.0e-3 / 2
area = math.pi * a_n ** 2
l_eff = 70e-3 + 1.5 * a_n


def helmholtz(c, area, vol, l_eff):
    return c / (2 * math.pi) * math.sqrt(area / (vol * l_eff))


f_full = helmholtz(344.0, area, 435e-6, l_eff)
check("Helmholtz: full-bottle wavelength (m) vs paper's upper 1.5 m", 344.0 / f_full, 1.5, 0.05)
# volume giving a 1.5 m and 1.0 m wavelength must lie inside the 435 cm^3 bottle
v_15 = area / (l_eff * (2 * math.pi * (344 / 1.5) / 344) ** 2)
v_10 = area / (l_eff * (2 * math.pi * (344 / 1.0) / 344) ** 2)
holds("Helmholtz: 1-1.5 m wavelengths need air volumes inside the bottle", 0 < v_10 < v_15 < 435e-6,
      f"f(435 cm^3) = {f_full:.0f} Hz; V(1.5 m) = {v_15*1e6:.0f} cm^3, V(1.0 m) = {v_10*1e6:.0f} cm^3")

# ---------------------------------------------------------------- 11. Clamped circular plate (Chladni)
# Leissa, Vibration of Plates, NASA SP-160 (1969), Table 2.1: lambda^2 = w a^2 sqrt(rho h / D)
# (n nodal diameters, s nodal circles): (0,0) 10.2158 (scan reads 10.2168), (1,0) 21.26,
# (2,0) 34.88, (0,1) 39.771.


def bessel_j(n, x):
    return sum((-1) ** k * (x / 2) ** (2 * k + n) / (math.factorial(k) * math.factorial(k + n)) for k in range(40))


def bessel_i(n, x):
    return sum((x / 2) ** (2 * k + n) / (math.factorial(k) * math.factorial(k + n)) for k in range(40))


def clamped_roots(n, count):
    fn = lambda lam: bessel_j(n, lam) * bessel_i(n + 1, lam) + bessel_i(n, lam) * bessel_j(n + 1, lam)
    roots, lam, step = [], 0.5, 0.01
    while len(roots) < count:
        if (fn(lam) > 0) != (fn(lam + step) > 0):
            roots.append(bisect(fn, lam, lam + step))
        lam += step
    return roots


r0 = clamped_roots(0, 2)
check("clamped plate (0,0) lambda^2", r0[0] ** 2, 10.2158, 2e-4)
check("clamped plate (0,1) lambda^2", r0[1] ** 2, 39.771, 2e-4)
check("clamped plate (1,0) lambda^2", clamped_roots(1, 1)[0] ** 2, 21.26, 1e-3)
check("clamped plate (2,0) lambda^2", clamped_roots(2, 1)[0] ** 2, 34.88, 1e-3)

print()
if FAIL:
    print(f"{len(FAIL)} check(s) failed: {FAIL}")
    sys.exit(1)
print("all checks passed")
