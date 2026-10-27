"""Checks for strange_radiation.md. Run: python3 -I strange_radiation_checks.py
Standard library only. Constants CODATA 2018 (exact SI values where defined).
Every published value below was read this session from the source named beside it."""
import math, sys

h = 6.62607015e-34; hbar = h / (2 * math.pi); e = 1.602176634e-19; c = 2.99792458e8
mu0 = 1.25663706212e-6; me = 9.1093837015e-31; alpha = 7.2973525693e-3; kB = 1.380649e-23
NA = 6.02214076e23; u = 1.66053906660e-27; muN_eV = 3.15245125844e-8  # nuclear magneton, eV/T
a0 = 5.29177210903e-11; MeV = 1e6 * e; R_gas = 8.314462618
fails = 0

def check(name, got, want, rtol):
    global fails
    ok = abs(got - want) <= rtol * abs(want)
    fails += not ok
    print(f"{'ok  ' if ok else 'FAIL'} {name}: got {got:.6g}, published {want:.6g} (tol {rtol:.1%})")

def hdr(s): print(f"\n== {s} ==")

# 1. Dirac charge (PDG 2024 rev. 95, Eq. 95.2: Q_E Q_M = 2*pi in natural units)
hdr("1 Dirac magnetic charge")
def dirac_charge_si(n=1):  # A*m; flux mu0*g = n h/e
    return n * h / (mu0 * e)
gD = dirac_charge_si()
check("g_D [A m] (de Lacheze-Murel 2022 footnote: 3.29e-9 SI)", gD, 3.29e-9, 0.005)
check("g_D/(e c) = 1/(2 alpha) (Mavromatos-Mitsou Eq. 11: 68.5)", gD / (e * c), 68.5, 0.001)
check("(g_D/e)^2 ionisation ratio (Mavromatos-Mitsou: ~4,700)", (gD / (e * c)) ** 2, 4700, 0.01)

def energy_gain_J(B_T, L_m, n=1):  # work done on a pole by a uniform field
    return dirac_charge_si(n) * B_T * L_m
check("energy gain 1 G over 0.1 m [keV] (de Lacheze-Murel 2022: 205 keV)",
      energy_gain_J(1e-4, 0.1) / e / 1e3, 205, 0.005)
print(f"   in Earth's field 50 uT: {energy_gain_J(50e-6, 1.0)/MeV:.2f} MeV per metre")

def monopole_field_T(r_m, n=1):  # Coulomb-like B = mu0 g / (4 pi r^2)
    return mu0 * dirac_charge_si(n) / (4 * math.pi * r_m ** 2)
for lab, r in (("1 fm", 1e-15), ("10 fm", 1e-14), ("Fe nearest neighbour 2.48 A", 2.48e-10), ("Bohr radius", a0)):
    print(f"   B at {lab}: {monopole_field_T(r):.3e} T")
r33 = math.sqrt(mu0 * gD / (4 * math.pi * 33.0)); print(f"   radius where B = 33 T (alpha-Fe hyperfine field): {r33*1e9:.2f} nm")
Bc = me ** 2 * c ** 2 / (e * hbar); rc = math.sqrt(mu0 * gD / (4 * math.pi * Bc))
print(f"   radius where B = m_e^2c^2/(e hbar) = {Bc:.3e} T: {rc*1e15:.0f} fm")
def pole_moment_energy_eV(mu_nuclear, r_m):  # |mu . B| of a nuclear moment at distance r from a Dirac pole
    return mu_nuclear * muN_eV * monopole_field_T(r_m)
print(f"   57Fe ground moment 0.0902 muN at 10 fm from a pole: |U| = {pole_moment_energy_eV(0.0902,1e-14)/1e3:.1f} keV")

# 2. 57Fe Mossbauer: moments from IAEA LiveChart (ENSDF), 14.4129 keV; Preston-Hanna-Heberle 1962 abstract
hdr("2 57Fe hyperfine sextet")
mu_g, I_g, mu_x, I_x, Eg = 0.0902, 0.5, -0.15531, 1.5, 14.4129e3
check("mu1/mu0 from ENSDF moments vs Preston 1962 (-1.715 +- 0.004)", mu_x / mu_g, -1.715, 0.006)
def line_velocities_mm_s(B_T):  # six allowed dm=0,+-1 lines, pure magnetic splitting, no isomer shift
    gg, gx = mu_g / I_g, mu_x / I_x
    out = []
    for mx, mg in ((1.5, .5), (.5, .5), (-.5, .5), (.5, -.5), (-.5, -.5), (-1.5, -.5)):
        dE = (-gx * mx + gg * mg) * muN_eV * B_T
        out.append(dE / Eg * c * 1e3)
    return sorted(out)
v = line_velocities_mm_s(33.0)
print("   lines at 33.0 T [mm/s]:", " ".join(f"{x:+.3f}" for x in v), f"; outer split {v[-1]-v[0]:.3f}")
# Urutskoev et al. 2002 p.723: standard 330.42 kG, N +0.24 kG, S -0.29 kG, error 0.012 kG, "relative change ~8e-4"
H0, dN, dS = 330.42, 0.24, -0.29
check("mean |dH|/H (Urutskoev 2002: ~8e-4)", (abs(dN) + abs(dS)) / 2 / H0, 8e-4, 0.10)
print(f"   outer-line shift for +0.24 kG: {(v[-1]-v[0])/2*dN/H0*1e3:.2f} um/s; for -0.29 kG: {(v[-1]-v[0])/2*dS/H0*1e3:.2f} um/s")
dHdT = -40.0  # Oe/K, quoted by Ivoilov & Urutskoev 2004 from Preston 1962
for lab, dH in (("2002 N", 240), ("2002 S", -290), ("2004 No.6(N) 280+-70", 280), ("2004 No.5(S) 460+-70", 460)):
    print(f"   {lab} Oe = temperature difference of {abs(dH/dHdT):.1f} K at -40 Oe/K")
def effective_field(B_hf, b_par, b_perp):  # |B_hf + b| for b split along/across the magnetisation
    return math.hypot(B_hf + b_par, b_perp)
# self-check: a field across the magnetisation shifts |B| upward for either sign
assert effective_field(33.0, 0, 1.26) > 33.0 and effective_field(33.0, 0, -1.26) > 33.0
b_need = math.sqrt(2 * 33.0 * 0.024)
sigma = 2 * b_need / mu0  # A m per m^2 of a uniform pole sheet giving b_need on one side
print(f"   transverse field giving +0.024 T at second order: {b_need:.2f} T; pole sheet {sigma/gD/1e4:.2e} Dirac poles/cm^2")

# 3. second-order Doppler coefficient, Debye model (theta_D = 400 K, Preston 1962) at 308 K (CEMS detector at 35 C)
hdr("3 second-order Doppler")
def debye_cv_ratio(T, theta, n=2000):  # C_v / (3 k_B) per atom
    x = theta / T; s = 0.0
    for i in range(1, n + 1):
        t = (i - 0.5) * x / n
        s += t ** 4 * math.exp(t) / (math.expm1(t)) ** 2
    return 3 / x ** 3 * s * x / n
def sod_coeff_um_s_K(T, theta, A=56.935):  # |d(delta)/dT| = C_v/(2 M c) in velocity units
    return 3 * kB * debye_cv_ratio(T, theta) / (2 * A * u * c) * 1e6
check("|d delta/dT| at 308 K (Ivoilov-Urutskoev 2004: 0.65 um/s/K)", sod_coeff_um_s_K(308, 400), 0.65, 0.05)

# 4. Monopole stopping power (Patrizii & Spurio 2015, Eqs. 11-12)
hdr("4 monopole energy loss")
def fermi_velocity_cm_s(rho_g_cm3, Z, A):
    Ne = rho_g_cm3 / A * NA * Z * 1e6  # m^-3
    return hbar / me * (3 * math.pi ** 2 * Ne) ** (1 / 3) * 100
check("v_F silicon at 4.3 g/cm3 (Patrizii-Spurio: 3.9e8 cm/s)", fermi_velocity_cm_s(4.3, 14, 28.086), 3.9e8, 0.01)
def dEdx_slow_GeV_cm(beta, rho, Z, A, n=1, spin_factor=1.37):  # Ahlen-Kinoshita, 1e-3 < beta < 1e-2, Gaussian units
    Ne = rho / A * NA * Z  # cm^-3
    e2 = 1.43996448e-7  # e^2 in eV cm
    g2e2 = (n / (2 * alpha)) ** 2 * e2 ** 2  # g^2 e^2 with g = n e/(2 alpha)
    vF = fermi_velocity_cm_s(rho, Z, A)
    mevF = 0.51099895e6 * vF / (c * 100)  # m_e c v_F in eV
    L = math.log(2 * vF / (alpha * c * 100)) - 0.5  # 2 m_e v_F a0 / hbar = 2 v_F/(alpha c)
    return spin_factor * 2 * math.pi * Ne * g2e2 * beta / mevF * L / 1e9
# shortcut: free-electron-gas model counts all Z electrons; order of magnitude only for AgBr
for lab, rho, Z, A in (("Si 2.33", 2.33, 14, 28.086), ("AgBr 6.47", 6.47, 82, 187.77)):
    s = dEdx_slow_GeV_cm(1e-3, rho, Z, A)
    print(f"   beta=1e-3 in {lab} g/cm3: {s*1e3:.0f} MeV/cm -> {s*0.1*1e3:.0f}-{s*0.3*1e3:.0f} MeV over 1-3 mm (claimed 700-1000 MeV)")
lo = dEdx_slow_GeV_cm(1e-3, 6.47, 82, 187.77) * 0.3; hi = dEdx_slow_GeV_cm(1e-3, 2.33, 14, 28.086) * 0.1
print(f"   claimed 0.7-1 GeV is {0.7/lo:.0f}-{1.0/hi:.0f} x the slow-monopole loss over the track")
# Urutskoev 2002: beta = 1e-3 from track + pulse; scintillator time-of-flight 20-40 m/s
print(f"   beta=1e-3 is {1e-3*c:.0f} m/s, i.e. {1e-3*c/40:.0f}-{1e-3*c/20:.0f} x the 20-40 m/s time-of-flight value")

# 5. Energy budget of the 2002 claims
hdr("5 energy budget")
W_bank = 50e3  # J, Urutskoev 2002 p.702
for N in (1e19, 1e20):
    print(f"   {N:.0e} transformations x 0.1 MeV bound = {N*0.1*MeV/1e3:.0f} kJ vs {W_bank/1e3:.0f} kJ bank; 48Ti mass {N*48e-3/NA*1e6:.1f} mg")
def ti48_after_burnout(f, a48=0.7372):  # all lost Ti assumed to be 48Ti (Urutskoev 2002 Fig. 10 assumption)
    return (a48 - f) / (1 - f)
print(f"   48Ti abundance after 4% burnout: {ti48_after_burnout(0.04):.4f}; burnout giving 65%: {(0.7372-0.65)/(1-0.65):.3f}")
check("natural Ti abundances sum (IAEA LiveChart)", 8.25 + 7.44 + 73.72 + 5.41 + 5.18, 100.0, 0.001)
# Priem et al. 2008 Table 3: electrodes 18-8 steel (Fe 72, Cr 20, Ni 8); shot 2 Fe 5.16 Cr 1.1 Ni 0.52; shot 3 Fe 6.6 Cr 1.3 Ni 0.67
def steel_fingerprint(fe, cr, ni, steel=(72.0, 20.0, 8.0)):
    k = fe / steel[0]
    return k, cr / (k * steel[1]), ni / (k * steel[2])
for lab, fe, cr, ni in (("shot 2", 5.16, 1.1, 0.52), ("shot 3", 6.6, 1.3, 0.67)):
    k, rc_, rn = steel_fingerprint(fe, cr, ni)
    print(f"   Priem {lab}: steel fraction {k:.3f}; Cr obs/expected {rc_:.2f}; Ni obs/expected {rn:.2f}")
    assert 0.6 < rc_ < 1.2 and 0.8 < rn < 1.2

# Klimov 2024 JCMNS 38 p.290: measured COP 1.62-2.82 (exps 3-7), 0.73-0.80 (exps 1-2); authors add x1.5-2 for unmeasured losses
print(f"   Klimov 2024: measured COP 1.62-2.82 becomes {1.62*1.5:.2f}-{2.82*2:.2f} only after the unmeasured x1.5-2 loss factor")
# Klimov 2016 JCMNS 19 p.160: Li7/Li6 13.192 -> 18.996
print(f"   Klimov 2016 Li7/Li6 change: x{18.996/13.192:.3f}")

# 6. Decay rates
hdr("6 decay-rate changes")
check("7Be in C60 vs Be metal (Ohtsuki 2004: 0.83%)", (53.12 - 52.68) / 53.12 * 100, 0.83, 0.01)
print(f"   bare/neutral 187Re half-life ratio (Bosch 1996: 42 Gyr vs 32.9 yr): {42e9/32.9:.2e}")
def step_after_suppression(T_half_d, dt_d):  # fractional activity excess after decays stop for dt (atoms not used up)
    return math.expm1(math.log(2) * dt_d / T_half_d)
for lab, T in (("60Co (1925.28 d, IAEA)", 1925.28), ("137Cs (30.08 y, IAEA)", 30.08 * 365.25)):
    print(f"   2-day total suppression of {lab}: permanent step {step_after_suppression(T, 2):.2e} vs Pomme 2016 limit 6e-6 to 8e-5")

print(f"   60Co step / Pomme limits: {step_after_suppression(1925.28,2)/8e-5:.0f}-{step_after_suppression(1925.28,2)/6e-6:.0f} x")

# 7. Air fluorescence brightness
hdr("7 air fluorescence")
Y337 = 5.61  # photons/MeV, 337 nm band, 1013 hPa, 293 K (AIRFLY, Ave et al. 2013)
rho_air = 101325 * 28.965e-3 / (R_gas * 293.15)
def column_rayleigh(dose_Gy_h, path_m, yield_per_MeV=Y337):  # isotropic emission, column rate / 1e10
    rate = dose_Gy_h / 3600 * rho_air / MeV * yield_per_MeV  # photons s^-1 m^-3
    return rate * path_m / 1e10
for D in (1, 10, 100):
    print(f"   {D:>3} Gy/h over 100 m of air: 337-nm band alone {column_rayleigh(D, 100)/1e3:.2f} kR (IBC I = 1 kR, IBC III = 100 kR)")
check("AIRFLY yield used", Y337, 5.61, 0.001)

print(f"\n{'ALL CHECKS PASSED' if not fails else f'{fails} CHECK(S) FAILED'}")
sys.exit(1 if fails else 0)
