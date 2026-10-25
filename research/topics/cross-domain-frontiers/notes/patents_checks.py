"""Recompute the numbers stated in patents (2020-2026) against standard physics.

Run: python3 -I patents_checks.py
Each block prints the patent's own number, the standard calculation, and the ratio.
Inputs marked ASSUMED are not stated in the patent; everything else is quoted from it.
"""
from math import pi, sqrt, asin, degrees, log, exp

G = 6.67430e-11
C = 299_792_458.0
MU0 = 1.25663706212e-6
EPS0 = 8.8541878128e-12
KB = 1.380649e-23
H = 6.62607015e-34
HBAR = H / (2 * pi)
E = 1.602176634e-19
MUB = 9.2740100783e-24
R_GAS = 8.314462618
INCH = 0.0254
TORR = 133.322368


def head(t):
    print("\n== " + t)


# ---------------------------------------------------------------- fringe propulsion
head("US 11,511,891 B2 (Aurigema & Buhler): asymmetrical electrostatic pressure")
gap = 0.25 * INCH
v = 25e3
e_field = v / gap
p_es = 0.5 * EPS0 * e_field**2            # standard electrostatic (Maxwell) pressure
area8, area15 = 8 * INCH**2, 15 * INCH**2
print(f"patent: eps0*E^2 at 1 MV/m 'order of 1 Pa'; standard eps0E^2/2 = {0.5*EPS0*1e12:.2f} Pa, eps0E^2 = {EPS0*1e12:.2f} Pa")
print(f"field across 0.25 in at 25 kV = {e_field/1e6:.2f} MV/m, pressure = {p_es:.1f} Pa")
print(f"one-sided pressure force on 8 in^2 = {p_es*area8*1e3:.0f} mN; measured 0.237 mN = {0.237e-3/(p_es*area8):.1e} of it")
print(f"area ratio 15/8 = {area15/area8:.3f}; measured 421/237 = {421/237:.3f}")
print("standard: closed-surface integral of the Maxwell stress of a static, neutral, isolated system is 0 (field ~1/r^3);")
print("          net thrust of a static field is 0, so the bound is a systematic budget, not a thrust/W")
m_dev = (0.030, 0.040)                      # The Debrief: 30-40 g articles, '1 g' claimed, '10 mN' best stack
print(f"secondary (Debrief): 1 g on 30-40 g needs {m_dev[0]*9.80665*1e3:.0f}-{m_dev[1]*9.80665*1e3:.0f} mN vs stated best 10 mN")

head("US 11,961,666 B2 (Purvis): pulsed e-field propulsion")
f_claim, p_claim = 7.72, 144.0
print(f"patent: {f_claim} N from {p_claim} W = {p_claim/f_claim:.2f} W/N (patent 18.65)")
photon = p_claim / C                        # gravitics::propulsion::photon_thrust
print(f"photon-rocket ceiling P/c = {photon:.2e} N; claim / ceiling = {f_claim/photon:.2e}")
l_h, r_ohm = 1e-6, 1.0
print(f"own circuit: envelope decay 2L/R = {2*l_h/r_ohm*1e6:.1f} us (patent 2 us); 4L/(R^2 C)=2 -> C = {4*l_h/(2*r_ohm**2)*1e6:.1f} uF")
print(f"own impulse field mu0/4pi * Idot = {1e-7*1e7:.1f} V/m per unit geometry factor (patent: Idot ~1e7 A/s)")
for name, f, p in [("Biefeld-Brown", 50e-3, 40.0), ("Cannae", 9e-3, 10.5), ("Eagleworks", 0.091e-3, 17.0)]:
    print(f"  prior art as cited, {name}: {f/p*1e3:.3f} mN/W = {f/p/(1/C):.1e} x photon rocket")

head("US 11,149,719 B2 (Von Bargen): EdDrive internal-mass drive")
d, dv = 36 * INCH, 2.045 * 0.44704
print(f"patent: 36 in in 1 s = {d:.4f} m; 2.045 mph = {dv:.4f} m/s; 1 kg x that / 1 s = {1.0*dv:.4f} N (patent 914.4 mN)")
print(f"914.4 mN / 12 s cycle = {914.4/12:.1f} mN (patent 76.2)")
m_host, m_proj, stroke = 1.0, 0.5, 0.30     # projectile 50% of host (patent); stroke ASSUMED 0.30 m
print(f"standard: centre of mass fixed in free space; hull excursion = m_p/(M+m_p) x stroke = {m_proj/(m_host+m_proj)*stroke:.3f} m, returns each cycle")
print("on a table any net drift is friction asymmetry (Dean-drive mechanism); thrust in vacuum = 0")

# ---------------------------------------------------------------- gravitational waves
head("US 12,494,696 B2 (Unspace): semi-coherent gravitational-wave beam")
rho_w = 19_300.0
r, t = 0.2, 0.1                              # patent upper sizes (m)
v_tip = 90.88                                # patent upper tip speed
omega = v_tip / r
f_gw = omega / pi                            # twice the rotation frequency
m_bar = rho_w * (2 * r) * t * t              # ASSUMED most asymmetric shape: a tungsten bar 2R x T x T
p_gw = 2 / 45 * G * m_bar**2 * (2 * r)**4 * omega**6 / C**5   # gravitational::laboratory::spinning_rod_power
print(f"tip 90.88 m/s at R = 0.2 m: omega = {omega:.0f} rad/s ({omega*60/(2*pi):.0f} RPM), f_GW = {f_gw:.0f} Hz")
print(f"bar mass {m_bar:.1f} kg; quadrupole power = {p_gw:.2e} W; energy in 1 s = {p_gw:.2e} J vs claimed > 6e-11 J (ratio {6e-11/p_gw:.1e})")
omega_max = 32_000 * 2 * pi / 60
p_max = 2 / 45 * G * m_bar**2 * (2 * r)**4 * omega_max**6 / C**5
print(f"at 32,000 RPM: tip {omega_max*r:.0f} m/s, hoop stress rho v^2 = {rho_w*(omega_max*r)**2/1e9:.1f} GPa (tungsten UTS ~1-1.5 GPa); P = {p_max:.2e} W")
print("an axisymmetric disc has no changing quadrupole: P = 0 exactly")
for f in (500.0, 1700.0):
    print(f"lambda = c/f at {f:.0f} Hz = {C/f:.2e} m (patent 2e5-6e5 m)")
lam, w0 = C / 1700.0, 1.5e-3
print(f"Gaussian beam, waist 1.5 mm at lambda {lam:.2e} m: Rayleigh range = {pi*w0**2/lam:.1e} m, divergence lambda/(pi w0) = {lam/(pi*w0):.1e} rad -> no beam")
h_1km = 2 * G / (C**4 * 1e3) * (m_bar * (2 * r)**2 / 12) * (2 * omega)**2 * 2
print(f"strain at 1 km (r = 1000 m, far zone needs r >> lambda, so this overstates it): h ~ {h_1km:.1e} vs LIGO ~1e-23")

# ---------------------------------------------------------------- light: laser-plasma
head("US 11,438,999 B2 (Leemans et al., UC): frozen-gas plasma channels")
lam_l, w = 800e-9, 25e-6
print(f"patent: Rayleigh range 'a few mm' for 800 nm, 25 um spot; pi w^2/lambda = {pi*w**2/lam_l*1e3:.2f} mm")
n0 = 1e17                                     # cm^-3, patent wall density
print(f"cold wave-breaking field at 1e17 cm^-3: E0 = 96 sqrt(n0) V/m = {96*sqrt(n0)/1e9:.1f} GV/m (Esarey 2009)")


def deposition_temp(p_bar, t_tp, p_tp, dh_sub):
    """Clausius-Clapeyron on the sublimation curve from the triple point."""
    return 1 / (1 / t_tp - R_GAS / dh_sub * log(p_bar / p_tp))


p10psi = 10 * 0.0689476
t_co2 = deposition_temp(p10psi, 216.58, 5.185, 26.1e3)   # CO2 triple point; dH_sub ~26 kJ/mol (approx)
t_n2o = deposition_temp(p10psi, 182.34, 0.878, 23.0e3)   # N2O triple point; dH_sub ~23 kJ/mol (approx)
print(f"CO2 at 10 psi deposits below {t_co2-273.15:.0f} C; patent says '-50 C or lower'; -50 C is above the triple point (-56.6 C): no solid at low pressure")
print(f"N2O at 10 psi deposits below {t_n2o-273.15:.0f} C; patent says '-155 C or lower' (consistent, conservative)")
for p in (5, 30):
    n_h2 = p * TORR / (KB * 293.15) * 1e-6
    print(f"second gas {p} torr H2 at 293 K: {n_h2:.1e} molecules/cm^3 -> n_e {2*n_h2:.1e} cm^-3 if fully ionized")

# ---------------------------------------------------------------- magnetics
head("US 11,543,474 B1 (Sandia & QuSpin): optically pumped gradient magnetometer")
g_j, g_i = 2.00233113, -0.000995141           # 87Rb, nuclear g in Bohr-magneton units
g_f2 = g_j / 4 + 3 * g_i / 4
g_f1 = -g_j / 4 + 5 * g_i / 4
coef = (2 * g_f2 - 1 * g_f1) * MUB / H * 1e-9
print(f"|1,1>-|2,2> low-field shift = {coef:.3f} Hz/nT (patent 3 x 7 = 21 Hz/nT); F=2 ratio {g_f2*MUB/H*1e-9:.3f} Hz/nT (patent 7)")
nu_hf = 6_834_682_610.904
print(f"cell spacing = microwave wavelength c/nu_hf = {C/nu_hf*100:.2f} cm (patent 4.4 cm)")
beta, delta, gamma = 547.9, 0.52, -0.0013    # patent Table 1 (Vanier 1982), N2
for t_c in (0.0, 100.0):
    shift = beta + delta * t_c + gamma * t_c**2
    print(f"N2 shift at {t_c:.0f} C = {shift:.1f} Hz/torr; 30 vs 15 torr offset = {15*shift/1e3:.2f} kHz (patent measured 7 kHz at ~100 C)")
p_eq = 15 * 273.15 / 293.15                  # ASSUMED fill at 20 C, coefficient referred to 0 C density
print(f"if fill pressures are at 20 C: {p_eq:.2f} torr-equivalent -> {p_eq*(beta+delta*100+gamma*1e4)/1e3:.2f} kHz")

head("US 11,852,701 B2: fibre-coupled NV-diamond magnetometer")
p_fl, lam_fl = (0.3e-3, 1e-3), 700e-9        # fluorescence power from patent; ~700 nm ASSUMED band centre
rates = [p / (H * C / lam_fl) for p in p_fl]
print(f"photon rate 0.3-1 mW: {rates[0]:.1e}-{rates[1]:.1e} /s; R x 12 ns = {rates[0]*12e-9:.1e}-{rates[1]*12e-9:.1e} NV (patent ~1e8; a lower bound, unit collection)")
gamma_nv = 28.024e9                          # Hz/T
lw, contrast = 1e6, 0.01                     # ASSUMED ensemble ODMR linewidth and contrast
for rate in rates:
    eta = 4 / (3 * sqrt(3)) * lw / (gamma_nv * contrast * sqrt(rate))
    print(f"CW-ODMR shot-noise limit (Barry et al. 2020) at R={rate:.1e}: {eta*1e12:.0f} pT/sqrt(Hz) (patent: sub-nT/sqrt(Hz))")
vol = 2e-9                                    # 2 mm^3, the patent's largest sample
print(f"energy resolution at 1 nT/sqrt(Hz), 2 mm^3: {(1e-9)**2*vol/(2*MU0*HBAR):.1e} hbar (bound is ~hbar)")

# ---------------------------------------------------------------- plasma
head("US 11,758,640 B2 (Zap Energy): sheared-flow Z-pinch")
molecules = TORR * 1e-3 / (KB * 273.15)
print(f"1 Torr-L at 273 K = {molecules:.3e} molecules (patent 2.5e19; ratio {2.5e19/molecules:.2f})")
print(f"500 vs 10 Torr-L in equal time = {500/10:.0f}x rate (patent 1/50 time)")
i_p, t_tot_ev, n_e = 200e3, 2e3, 1e23       # secondary (OSTI PRL): ~200 kA, Te+Ti ~1-2 keV each -> 2 keV ASSUMED sum, 1e17 cm^-3
n_line = MU0 * i_p**2 / (8 * pi * t_tot_ev * E)   # plasmoid::pinch::bennett_line_density
print(f"Bennett line density at 200 kA, Te+Ti = 2 keV: {n_line:.2e} /m -> radius at 1e23 m^-3 = {sqrt(n_line/(pi*n_e))*1e3:.1f} mm")

head("US 11,195,627 B2 (TAE): tunable neutral beams")
print(f"15 -> 40 keV at constant current: injected power x{40/15:.2f} (patent: up to 2x heating power, after losses)")

# ---------------------------------------------------------------- gravity sensing
head("US 11,269,111 B2 (Birmingham / Delta G): gravity gradiometer")
print(f"patent tilt asin(1/15) = {degrees(asin(1/15)):.2f} deg (patent 3.8)")
k_eff = 4 * pi / 780.241e-9
t_fall = sqrt(2 * 0.15 / 9.81)
t_pulse = t_fall / 2                          # dropped atoms: 2T fits the 15 cm fall
sens = k_eff * t_pulse**2                     # gravimetry::atom::phase_sensitivity
print(f"15 cm drop: fall {t_fall*1e3:.0f} ms, T = {t_pulse*1e3:.0f} ms, k_eff T^2 = {sens:.3e} rad/(m/s^2)")
n_at, c_con, base, t_cyc = 1e6, 0.5, 1.0, 0.3  # ASSUMED atoms, contrast; baseline 1 m (patent 0.8-1.2 m); cycle 300 ms (patent)
dphi = 1 / (c_con * sqrt(n_at))
dgrad = sqrt(2) * dphi / sens / base / 1e-9   # two clouds, uncorrelated shot noise
asd = dgrad * sqrt(t_cyc)
print(f"shot-noise gradient per shot {dgrad:.1f} E, ASD {asd:.1f} E/sqrt(Hz), time to 1 E = {asd**2:.0f} s")
print(f"Rb D2 natural linewidth 1/(2 pi 26.24 ns) = {1/(2*pi*26.24e-9)/1e6:.2f} MHz (patent 6 MHz)")
print(f"Rb hyperfine 6.835 GHz (patent 6.8 GHz)")
