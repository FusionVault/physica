"""Recompute the published numbers behind light_and_gravity.md.

Standard library only. Run as:  python3 -I light_and_gravity_checks.py
Every check prints computed vs published and asserts the stated tolerance.
"""
from math import pi, sqrt, log, exp, atan, cos, acos, sin, tan, atanh, atan2

# Constants as in crates/metrology (SI 2019, CODATA 2018, IAU 2015 B3, IERS 2010)
C = 299_792_458.0
G = 6.674_30e-11
HBAR = 1.054_571_817e-34
K_B = 1.380_649e-23
GM_SUN = 1.327_124_4e20          # IAU 2015 nominal
GM_EARTH = 3.986_004_418e14
AU = 149_597_870_700.0
PARSEC = 3.085_677_581_491_367e16
MPC = 1e6 * PARSEC
YEAR = 365.25 * 86_400.0
DAY = 86_400.0
R_SUN = 6.957e8                   # IAU 2015 B3 nominal solar radius
R_EARTH = 6.3781e6                # IAU 2015 B3 nominal equatorial radius
GM_JUP = 1.266_865_3e17           # IAU 2015 B3 nominal
R_JUP = 7.1492e7
MSUN = GM_SUN / G
MAS = pi / 180 / 3600 / 1000      # rad per milliarcsecond
UAS = MAS / 1000
ARCSEC = MAS * 1000
T_SUN = GM_SUN / C**3             # 4.925490947 µs

n_fail = 0


def check(name, got, want, rel=None, abs_=None):
    """Print and assert |got-want| within rel (relative) or abs_ (absolute)."""
    global n_fail
    if rel is not None:
        ok = abs(got / want - 1) <= rel
        tol = f"rel {rel:g}"
    else:
        ok = abs(got - want) <= abs_
        tol = f"abs {abs_:g}"
    flag = "ok " if ok else "FAIL"
    if not ok:
        n_fail += 1
    print(f"[{flag}] {name}: computed {got:.6g}, published {want:.6g} ({tol})")


def sigma(a, b, ua, ub=0.0):
    return abs(a - b) / sqrt(ua * ua + ub * ub)


def simpson(f, a, b, n=2000):
    h = (b - a) / n
    s = f(a) + f(b) + sum((4 if i % 2 else 2) * f(a + i * h) for i in range(1, n))
    return s * h / 3


# ---------- cosmology (flat LCDM, radiation neglected) ----------
def comoving(z, h0, om):
    """Transverse comoving distance D_M (m) in flat LCDM, H0 in km/s/Mpc."""
    dh = C / (h0 * 1e3 / MPC)
    return dh * simpson(lambda x: 1 / sqrt(om * (1 + x) ** 3 + 1 - om), 0.0, z)


def ang_dists(zl, zs, h0, om):
    dl, ds = comoving(zl, h0, om), comoving(zs, h0, om)
    return dl / (1 + zl), ds / (1 + zs), (ds - dl) / (1 + zs)


print("\n== 1. Optical clocks and chronometric geodesy ==")
# Bothwell et al. 2022 (arXiv:2109.12238): g/c^2 per mm
g_boulder = 9.796   # m/s^2, assumption (Boulder local g); paper quotes -1.09e-19/mm
check("JILA g/c^2 per mm", g_boulder / C**2 * 1e-3, 1.09e-19, rel=0.005)
print(f"      measured -9.8(2.3)e-20/mm is {sigma(9.8e-20, 1.09e-19, 2.3e-20):.2f} sigma from g/c^2")
# Zheng et al. 2023 Nat Commun: -12.4 ± 0.7 ± 2.5 e-19/cm vs -10.9e-19/cm
check("Wisconsin g/c^2 per cm", 9.80 / C**2 * 1e-2, 10.9e-19, rel=0.005)
print(f"      measured -12.4(0.7)(2.5)e-19/cm is {sigma(12.4, 10.9, 0.7, 2.5):.2f} sigma from g/c^2")
# Tokyo Skytree (Takamoto et al. 2020; numbers from Shinkai-Takamoto-Katori arXiv:2502.06104)
nu_sr = 429_228_004_229_873.0
frac = 21.18 / nu_sr
check("Skytree 21.18 Hz / nu_Sr vs gΔh/c^2 = 49337.1e-18", frac, 49_337.1e-18, rel=2e-4)
gbar = 49_337.1e-18 * C**2 / 452.596
check("Skytree implied mean g (m/s^2)", gbar, 9.797, rel=1e-3)
# 457 km fibre (PR Applied 21, L061001, 2024): ΔW = 3918.1(2.6) vs 3915.88(0.30) m^2/s^2
check("457 km link Δν/ν = ΔW/c^2", 3918.1 / C**2, 4.3595e-14, rel=1e-4)
print(f"      clock vs geodesy: {sigma(3918.1, 3915.88, 2.6, 0.30):.2f} sigma; 1 m^2/s^2 = {1/C**2:.3e} fractional")
check("1 cm height on Earth as Δν/ν", 9.81 * 0.01 / C**2, 1.09e-18, rel=0.005)
# L_G = W0/c^2 (IAU 2000 B1.9 defines L_G = 6.969290134e-10 from W0 = 62636856.0)
check("L_G = W0(IERS 2010)/c^2", 62_636_856.0 / C**2, 6.969_290_134e-10, rel=1e-9)
dW0 = 62_636_856.0 - 62_636_853.4   # IERS 2010 vs IAG 2015 conventional W0
check("W0 update in cm of height (≈27 cm per Sánchez et al.)", dW0 / 9.80 * 100, 26.5, rel=0.03)
# GRS80 normal potential U0 = (GM/E) arctan(e') + ω²a²/3 (Moritz 2000)
a, gm, om_e, b = 6_378_137.0, 3.986_005e14, 7.292_115e-5, 6_356_752.3141
E = sqrt(a * a - b * b)
U0 = gm / E * atan(E / b) + om_e**2 * a * a / 3
check("GRS80 U0 (m^2/s^2)", U0, 62_636_860.850, abs_=0.01)
print(f"      U0 - W0(IAG 2015) = {U0 - 62_636_853.4:.2f} m^2/s^2 (~{(U0-62_636_853.4)/9.8*100:.0f} cm)")
phi_c_eq = 0.5 * om_e**2 * a * a
print(f"      centrifugal potential at the equator {phi_c_eq:.1f} m^2/s^2 = {phi_c_eq/C**2:.3e} fractional")
# Galileo eccentric satellites (Delva et al. 2018): e = 0.162, period ~13 h, '400 ns peak-to-peak'
P = 20 * 86_164.0905 / 37        # 37 orbits per 20 sidereal days, assumption (ground-track repeat)
a_gal = (GM_EARTH * P * P / (4 * pi * pi)) ** (1 / 3)
pp = 2 * sqrt(GM_EARTH * a_gal) * 0.162 / C**2
check("Galileo redshift periodic term peak-to-peak (s)", pp, 400e-9, rel=0.06)
print(f"      a = {a_gal/1e3:.0f} km, apogee-perigee = {2*a_gal*0.162/1e3:.0f} km (paper: 'some 8500 km')")
F = -2 * sqrt(3.986_005e14) / C**2
check("GNSS eccentricity constant F (s/m^0.5)", F, -4.442_807_633e-10, rel=1e-8)
# Solar gravitational redshift seen from Earth (González Hernández et al. 2020): ~633.1 m/s
v_sun = GM_SUN / (C * R_SUN) - GM_SUN / (C * AU) - GM_EARTH / (C * R_EARTH)
check("solar GRS at Earth (m/s)", v_sun, 633.1, rel=2e-4)
print(f"      measured 638 ± 6 m/s is {sigma(638, v_sun, 6):.2f} sigma")

print("\n== 2. Light deflection, Shapiro delay, PPN ==")
check("solar-limb deflection (arcsec)", 4 * GM_SUN / (C**2 * R_SUN) / ARCSEC, 1.7517, rel=1e-3)


def defl_sun(chi, gamma=1.0, r=AU):
    """Deflection of a star at elongation chi from the Sun seen from distance r (Will 2014 §4.1.1)."""
    return (1 + gamma) / 2 * 4 * GM_SUN / (C**2 * r * sin(chi)) * (1 + cos(chi)) / 2


check("deflection at 90° from the Sun (mas, Will: 'still 4')", defl_sun(pi / 2) / MAS, 4.07, rel=0.005)
print(f"      Jupiter limb {4*GM_JUP/(C**2*R_JUP)/MAS:.2f} mas; Earth limb {4*GM_EARTH/(C**2*R_EARTH)/MAS:.3f} mas")
# Cassini 2002: b = 1.6 R_sun, spacecraft at 8.43 AU. Round-trip Shapiro delay
b, re = 1.6 * R_SUN, 8.43 * AU
dt_exact = 2 * 2 * T_SUN * log(4 * AU * re / b**2)
dt_will = 240 - 20 * log(1.6**2 / 8.43)
check("Cassini round-trip Shapiro delay exact vs Will eq. 63 (µs)", dt_exact * 1e6, dt_will, rel=0.01)
y_peak = 4 * T_SUN / b * 29.8e3   # db/dt ≈ Earth's orbital speed (assumption)
print(f"      Doppler y = 2(1+γ)(GM/c^3)(db/dt)/b ≈ {y_peak:.1e} near conjunction (order of magnitude)")
# Double pulsar (Kramer et al. 2021, PRX 11, 041050)
check("double pulsar Shapiro range r = T_sun m_B (µs)", T_SUN * 1.248868 * 1e6, 6.162, rel=0.004)
print(f"      {sigma(6.162, T_SUN*1.248868e6, 0.021):.2f} sigma; s = 1 - exp(-9.65) = {1-exp(-9.65):.6f} (paper 0.999936)")
print(f"      peak Shapiro delay ≈ 2 r z_s = {2*6.162*9.65:.0f} µs")
# Nordtvedt: Δ(mg/mi)_EM = -η[(Eg/m)_E - (Eg/m)_M], Eg/m = 4.6e-10 (Earth), 0.2e-10 (Moon) (Will 2014)
d_selfE = 4.6e-10 - 0.2e-10
eta = 2.1e-14 / d_selfE
check("η from Biskupek et al. 2021 Δ(mg/mi) (Fienga&Minazzoli: '5e-5')", eta, 5e-5, rel=0.1)
print(f"      ±{2.4e-14/d_selfE:.1e}; range signal 13.1 η m = {13.1*eta*1e3:.2f} mm")
eta_ppn = 4 * 6.2e-5 - 1.7e-4
print(f"      η = 4(β-1) - (γ-1) from the same LLR fit: {eta_ppn:.1e} ± {sqrt((4*7.2e-5)**2+(1.6e-4)**2):.1e}")
n_e = 2 * pi / YEAR
omega_ds = 1.5 * GM_SUN / (C**2 * AU) * n_e
check("de Sitter precession of the lunar orbit (arcsec/century, Will: 'about 2')", omega_ds * 100 * YEAR / ARCSEC, 1.92, rel=0.01)

# Deflection by any field: α = (1+γ)/c² ∫ ∇⊥U dl along the straight ray (Born approximation).
# Ray in the equatorial plane of an oblate body, spin axis ⟂ ray: U = GM/r + GM J2 R²/(2r³) in that plane.
def born_deflection(gm, j2, r_eq, b, gamma=1.0, n=4000):
    """∂U/∂b integrated along x = b tanθ; returns the deflection (rad)."""
    def dudb(theta):
        x = b * tan(theta)
        r = sqrt(x * x + b * b)
        g = gm * b / r**3 + 1.5 * gm * j2 * r_eq**2 * b / r**5   # -∂U/∂b
        return g * b / cos(theta) ** 2                              # dx = b sec²θ dθ
    return (1 + gamma) / C**2 * simpson(dudb, -pi / 2 + 1e-9, pi / 2 - 1e-9, n)


J2_JUP = 14_696.5063e-6   # Durante et al. 2020's 14696.5735 without its 6.72e-8 permanent tide; Iess et al. 2018 give 14696.572 (corrected 2026-10-25)
d_num = born_deflection(GM_JUP, J2_JUP, R_JUP, R_JUP)
d_ana = 4 * GM_JUP / (C**2 * R_JUP) * (1 + J2_JUP)
check("Jupiter limb deflection with J2, numerical Born integral vs 4GM/(c²b)(1+J2 R²/b²)", d_num, d_ana, rel=1e-6)
print(f"      J2 part at Jupiter's limb = {4*GM_JUP/(C**2*R_JUP)*J2_JUP/UAS:.0f} µas; "
      f"Earth's J2 part at its limb = {4*GM_EARTH/(C**2*R_EARTH)*1.0826e-3/UAS:.2f} µas")

print("\n== 3. Lensing ==")
kappa = 4 * G / (C**2 * AU) * MSUN / MAS
check("microlensing κ = 4G/(c^2 AU) (mas/M_sun)", kappa, 8.144, rel=1e-3)
# OGLE-2011-BLG-0462 (Sahu et al. 2025, arXiv:2503.07820)
thE, piE, DS = 5.496, 0.0944, 7.1
M = thE / (kappa * piE)
check("OGLE-2011-BLG-0462 mass (M_sun)", M, 7.15, rel=0.01)
pi_rel = thE * piE                 # mas
DL = 1 / (pi_rel + 1 / DS)          # kpc
check("OGLE-2011-BLG-0462 lens distance (kpc)", DL, 1.52, rel=0.01)
mu = sqrt(2.907**2 + 6.751**2)
check("t_E = θ_E/μ_rel (days)", thE / mu * YEAR / DAY, 273.1, rel=0.01)


def pacz(u):
    """Total point-lens magnification (Paczyński 1986)."""
    return (u * u + 2) / (u * sqrt(u * u + 4))


print(f"      peak A(u0 = 0.00271) = {pacz(0.00271):.0f}; max astrometric shift θ_E/(2√2) = {thE/(2*sqrt(2)):.2f} mas")
check("A(u=1) = 3/√5", pacz(1.0), 3 / sqrt(5), rel=1e-12)
# Euclid NGC 6505 ring (O'Riordan et al. 2025): θ_E = 2.5", z_l = 0.042, z_s = 0.406, σ* = 303 ± 15 km/s
zl, zs, th = 0.042, 0.406, 2.5 * ARCSEC
Dl, Ds, Dls = ang_dists(zl, zs, 70.0, 0.3)
sig_sis = C * sqrt(th * Ds / (4 * pi * Dls)) / 1e3
print(f"      NGC 6505: SIS σ from θ_E = {sig_sis:.0f} km/s vs measured 303 ± 15 ({sigma(sig_sis, 303, 15):.2f} sigma; H0-independent)")
check("NGC 6505 θ_E in kpc (paper 2.1; H0=70, Ωm=0.3)", Dl * th / (1e3 * PARSEC), 2.1, rel=0.03)
Mein = C**2 / (4 * G) * Dl * Ds / Dls * th**2 / MSUN
print(f"      projected mass inside θ_E = {Mein:.3e} M_sun (H0=70; scales as 1/H0)")
# H0LiCOW (Wong et al. 2020, Table 2): D_Δt and redshifts; H0 ∝ 1/D_Δt at fixed Ωm
lenses = {  # name: (z_d, z_s, D_dt, +err, -err) in Mpc
    "B1608+656": (0.6304, 1.394, 5156, 296, 236),
    "RXJ1131-1231": (0.295, 0.654, 2096, 98, 83),
    "HE0435-1223": (0.4546, 1.693, 2707, 183, 168),
    "SDSS1206+4332": (0.745, 1.789, 5769, 589, 471),
    "WFI2033-4723": (0.6575, 1.662, 4784, 399, 248),
    "PG1115+080": (0.311, 1.722, 1470, 137, 127),
}
wsum = hsum = 0.0
for name, (zd, zsrc, ddt, ep, em) in lenses.items():
    dd, ds_, dds = ang_dists(zd, zsrc, 70.0, 0.3)
    ddt70 = (1 + zd) * dd * ds_ / dds / MPC
    h = 70.0 * ddt70 / ddt
    s = h * 0.5 * (ep + em) / ddt
    wsum += 1 / s**2
    hsum += h / s**2
    print(f"      {name:14s} D_Δt(H0=70) = {ddt70:6.0f} Mpc -> H0 = {h:5.1f} ± {s:4.1f}")
check("H0LiCOW inverse-variance H0 (paper 73.3 +1.7 -1.8)", hsum / wsum, 73.3, rel=0.02)
print(f"      ± {1/sqrt(wsum):.1f} (symmetrised errors, D_Δt only, Ωm fixed at 0.3)")


def nfw_f(x):
    """Wright & Brainerd 2000 eq. 11 shape: Σ(x) = 2 r_s δ_c ρ_c f(x)."""
    if x < 1:
        return (1 - 2 / sqrt(1 - x * x) * atanh(sqrt((1 - x) / (1 + x)))) / (x * x - 1)
    if x > 1:
        return (1 - 2 / sqrt(x * x - 1) * atan(sqrt((x - 1) / (1 + x)))) / (x * x - 1)
    return 1 / 3


# shortcut: closed forms cancel catastrophically near x = 1 (rel. error ~1e-16/|x-1|²); a Rust port
# should switch to a series for |x - 1| < 1e-3. Here we only check continuity from both sides.
# x<1 and x>1 branches meet at 1/3 with one slope: their mean is 1/3 to O(ε²), their slope ≈ -0.4
check("NFW Σ shape: mean of both branches at 1 ± 1e-3", 0.5 * (nfw_f(1 - 1e-3) + nfw_f(1 + 1e-3)), 1 / 3, rel=1e-5)
check("NFW Σ shape: slope across x = 1", (nfw_f(1 + 1e-3) - nfw_f(1 - 1e-3)) / 2e-3, -0.4, rel=1e-3)

print("\n== 4. Black-hole shadows and photon rings ==")
th_vlti, th_keck = 5.125, 4.92   # µas, EHT Sgr A* Paper VI priors
check("Sgr A* θg from GRAVITY 2024 M, R0 (µas)",
      G * 4.2996e6 * MSUN / (C**2 * 8275.9 * PARSEC) / UAS, 5.125, rel=0.002)
for lab, tg, d in (("VLTI", th_vlti, -0.08), ("Keck", th_keck, -0.04)):
    print(f"      Sgr A* {lab}: 6√3 θg = {6*sqrt(3)*tg:.1f} µas; (1+δ)·6√3 θg = {(1+d)*6*sqrt(3)*tg:.1f} µas "
          f"(EHT VI Table 1 shadow ≈ 48–50 µas)")
tg87 = G * 6.5e9 * MSUN / (C**2 * 16.8 * MPC) / UAS
check("M87* θg for 6.5e9 M_sun at 16.8 Mpc (µas, EHT 2019: 3.8 ± 0.4)", tg87, 3.8, rel=0.02)
print(f"      M87* Schwarzschild shadow 6√3 θg = {6*sqrt(3)*tg87:.1f} µas; ring 42 ± 3 (2017), 43.9 ± 0.6 (2017–2021)")
print(f"      ring/θg: {42/tg87:.1f} (2017), {43.9/tg87:.1f} (2017–2021)")


def kerr_polar_shadow(a):
    """Shadow radius (units of M) seen along the spin axis: ξ = 0 spherical photon orbit."""
    lo, hi = 1.0 + sqrt(2) - 1e-12, 3.0
    f = lambda r: r**3 - 3 * r * r + a * a * r + a * a
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(lo) * f(mid) <= 0:
            hi = mid
        else:
            lo = mid
    r = 0.5 * (lo + hi)
    eta = r**3 * (4 * a * a - r * (r - 3) ** 2) / (a * a * (r - 1) ** 2)
    return sqrt(eta + a * a)


check("Kerr a→0 polar shadow radius = √27", kerr_polar_shadow(1e-4), sqrt(27), rel=1e-6)
check("Kerr a=1 polar shadow radius = 2(1+√2) (EHT VI: Kerr range ~4.8–5.2)", kerr_polar_shadow(1.0),
      2 * (1 + sqrt(2)), rel=1e-6)
print(f"      a = 0.5: {kerr_polar_shadow(0.5):.3f} M; a = 0.94: {kerr_polar_shadow(0.94):.3f} M")
r_ph = lambda a, s: 2 * (1 + cos(2 / 3 * acos(-s * a)))
check("equatorial prograde photon orbit a=1 (M)", r_ph(1.0, 1), 1.0, abs_=1e-12)
check("equatorial retrograde photon orbit a=1 (M)", r_ph(1.0, -1), 4.0, rel=1e-12)
lam = 1 / (3 * sqrt(3))          # Lyapunov exponent per unit coordinate time (c = G = M = 1)
t_half = pi * 3 * sqrt(3)        # half an orbit at Ω = 1/(3√3)
check("Schwarzschild photon-ring Lyapunov exponent γ = λ·T/2", lam * t_half, pi, rel=1e-12)
print(f"      each subring demagnified by e^γ = {exp(pi):.2f}")

print("\n== 5. S2 at Sgr A* (GRAVITY 2018, 2020, 2024) ==")
M_bh = 4.2996e6 * MSUN
a_s2 = 0.12502 * ARCSEC * 8275.9 * PARSEC
e_s2 = 0.88444
gm = G * M_bh
check("S2 period from Kepler (yr, paper: 16)", 2 * pi * sqrt(a_s2**3 / gm) / YEAR, 16.05, rel=0.01)
check("S2 Schwarzschild precession per orbit (arcmin, paper: 12')",
      6 * pi * gm / (C**2 * a_s2 * (1 - e_s2**2)) / ARCSEC / 60, 12.0, rel=0.03)
rp = a_s2 * (1 - e_s2)
check("S2 pericentre (AU, paper: 120)", rp / AU, 120, rel=0.01)
vp = sqrt(gm * (2 / rp - 1 / a_s2))
check("S2 pericentre speed (km/s, paper: ~7650)", vp / 1e3, 7650, rel=0.02)
z = gm / (C**2 * rp) + vp**2 / (2 * C**2)
check("S2 redshift at pericentre, gravitational + transverse Doppler (km/s, paper: ~200)", z * C / 1e3, 200, rel=0.05)

print("\n== 6. Optical analogue horizons ==")
alpha = 8e14
check("Philbin 2008 fibre Hawking T = ħα/(2π k_B) at α = 8e14 s^-1 (paper: 10^3 K)",
      HBAR * alpha / (2 * pi * K_B), 1e3, rel=0.05)


def n_silica(lam_um):
    """Fused silica, Malitson 1965 Sellmeier, valid 0.21–3.71 µm."""
    l2 = lam_um * lam_um
    return sqrt(1 + 0.6961663 * l2 / (l2 - 0.0684043**2) + 0.4079426 * l2 / (l2 - 0.1162414**2)
                + 0.8974794 * l2 / (l2 - 9.896161**2))


def ng_silica(lam_um, h=1e-5):
    return n_silica(lam_um) - lam_um * (n_silica(lam_um + h) - n_silica(lam_um - h)) / (2 * h)


ng_pump = ng_silica(0.8189)
print(f"      bulk silica at the pump (818.9 nm): n = {n_silica(0.8189):.4f}, n_g = {ng_pump:.4f}, "
      f"(n_g - n)/n = {(ng_pump-n_silica(0.8189))/n_silica(0.8189):.4f} (Drori note 41: ≈1e-2)")
for lam_nm in (1450, 231):
    s = 1 - n_silica(lam_nm / 1e3) / ng_pump
    print(f"      ω'/ω = 1 - n u/c at {lam_nm} nm = {s:+.4f} ({'positive' if s > 0 else 'negative'} co-moving frequency)")

print(f"\n{n_fail} failed")
assert n_fail == 0
