"""Recompute the published numbers behind gravity_and_gravimetry_frontiers.md.

Standard library only. Run as:  python3 -I gravity_frontiers_checks.py
Each check prints the computed value beside the published one and asserts the agreement the
notes claim. Sources are named on each block; DOIs and arXiv IDs are in the notes.
"""
from datetime import datetime
from math import acos, asin, atan, cos, exp, log, pi, radians, sin, sqrt

G = 6.674_30e-11
C = 299_792_458.0
GM_SUN = 1.327_124_400_18e20
T_SUN = GM_SUN / C**3  # GM_sun/c^3, s
L_SUN = GM_SUN / C**2  # GM_sun/c^2, m
GM_MOON = 4.902_800e12
GM_EARTH = 3.986_004_418e14
U = 1.660_539_066_60e-27  # atomic mass unit, kg
ALPHA = 7.297_352_569_3e-3
E2_4PIEPS0 = 2.307_077_552e-28  # e^2/(4 pi eps0), J m
MICROGAL = 1e-8
EOTVOS = 1e-9


def close(x, y, rel):
    return abs(x - y) <= rel * abs(y)


def show(label, got, want=None):
    print(f"{label:58s} {got:.6g}" + ("" if want is None else f"   (published {want})"))


# ---------------------------------------------------------------------------------------------
# 1. Longman (1959) lunisolar gravity tide, test vector from LongmanTide (Leeman, MIT licence).
# ---------------------------------------------------------------------------------------------
def longman(lat, lon_east, alt, when):
    """Return (moon, sun, total) vertical tidal acceleration in mGal, times 1 + h2 - 1.5 k2."""
    days = (when - datetime(1899, 12, 31, 12)).total_seconds() / 86400.0
    T = days / 36525.0  # Julian centuries from 1899-12-31 12:00
    t0 = when.hour + when.minute / 60 + when.second / 3600  # UT hour
    mu, M, S = 6.673e-8, 7.3537e25, 1.993e33  # cgs: G, moon and sun masses
    e, m = 0.05490, 0.074804  # lunar eccentricity; ratio of mean motions sun/moon
    c, c1 = 3.84402e10, 1.495e13  # mean distances moon and sun, cm
    h2, k2 = 0.612, 0.303  # Longman's Love numbers
    a = 6.378270e8  # equatorial radius, cm
    i = 0.08979719  # inclination of lunar orbit to ecliptic, rad
    omega = radians(23.452)  # obliquity
    L = -lon_east  # Longman counts west longitude positive
    lam = radians(lat)
    H = alt * 100.0
    s = 4.72000889397 + 8399.70927456 * T + 3.45575191895e-05 * T**2 + 3.49065850399e-08 * T**3
    p = 5.83515162814 + 71.0180412089 * T + 0.000180108282532 * T**2 + 1.74532925199e-07 * T**3
    h = 4.88162798259 + 628.331950894 * T + 5.23598775598e-06 * T**2
    N = 4.52360161181 - 33.757146295 * T + 3.6264063347e-05 * T**2 + 3.39369576777e-08 * T**3
    I = acos(cos(omega) * cos(i) - sin(omega) * sin(i) * cos(N))
    nu = asin(sin(i) * sin(N) / sin(I))
    t = radians(15.0 * (t0 - 12.0) - L)
    chi = t + h - nu
    cos_al = cos(N) * cos(nu) + sin(N) * sin(nu) * cos(omega)
    sin_al = sin(omega) * sin(N) / sin(I)
    xi = N - 2 * atan(sin_al / (1 + cos_al))
    sigma = s - xi
    l = (sigma + 2 * e * sin(s - p) + 1.25 * e * e * sin(2 * (s - p))
         + 3.75 * m * e * sin(s - 2 * h + p) + 11 / 8 * m * m * sin(2 * (s - h)))
    p1 = 4.90822941839 + 0.0300025492114 * T + 7.85398163397e-06 * T**2 + 5.3329504922e-08 * T**3
    e1 = 0.01675104 - 0.00004180 * T - 0.000000126 * T**2
    chi1 = t + h
    l1 = h + 2 * e1 * sin(h - p1)
    cth = sin(lam) * sin(I) * sin(l) + cos(lam) * (cos(I / 2)**2 * cos(l - chi) + sin(I / 2)**2 * cos(l + chi))
    cph = sin(lam) * sin(omega) * sin(l1) + cos(lam) * (cos(omega / 2)**2 * cos(l1 - chi1) + sin(omega / 2)**2 * cos(l1 + chi1))
    r = a / sqrt(1 + 0.006738 * sin(lam)**2) + H
    ap, ap1 = 1 / (c * (1 - e * e)), 1 / (c1 * (1 - e1 * e1))
    d = 1 / (1 / c + ap * e * cos(s - p) + ap * e * e * cos(2 * (s - p))
             + 15 / 8 * ap * m * e * cos(s - 2 * h + p) + ap * m * m * cos(2 * (s - h)))
    D = 1 / (1 / c1 + ap1 * e1 * cos(h - p1))
    gm = mu * M * r / d**3 * (3 * cth**2 - 1) + 1.5 * mu * M * r * r / d**4 * (5 * cth**3 - 3 * cth)
    gs = mu * S * r / D**3 * (3 * cph**2 - 1)
    f = (1 + h2 - 1.5 * k2) * 1e3  # gal -> mGal, times the gravimetric factor
    return gm * f, gs * f, (gm + gs) * f


gm, gs, gt = longman(40.7914, 282.1414, 370.0, datetime(2015, 4, 23))
show("Longman moon, State College 2015-04-23 00:00 UT, mGal", gm, 0.0324029651226)
show("Longman sun, mGal", gs, -0.0288682178454)
show("Longman total, mGal", gt, 0.00353474727722)
assert abs(gm - 0.0324029651226) < 1e-9 and abs(gs + 0.0288682178454) < 1e-9 and abs(gt - 0.00353474727722) < 1e-9

# ---------------------------------------------------------------------------------------------
# 2. Love numbers and tidal factors (IERS Conventions 2010, Tables 6.3 and ch. 7 step 1).
# ---------------------------------------------------------------------------------------------
h2, l2, k2 = 0.6078, 0.0847, 0.30190  # IERS ch.7 h(0), l(0); Table 6.3 anelastic Re k20
h3, k3 = 0.292, 0.093


def delta(n, h, k):  # gravimetric factor of degree n
    return 1 + 2 * h / n - (n + 1) * k / n


show("delta_2 (IERS h2 0.6078, k20 0.30190)", delta(2, h2, k2), "≈1.16")
show("delta_2 (Longman h2 0.612, k2 0.303)", delta(2, 0.612, 0.303))
show("delta_3 (h3 0.292, k3 0.093)", delta(3, h3, k3))
show("gamma_2 = 1 + k - h (tilt / levelling factor)", 1 + k2 - h2, "≈0.69 (Mäkinen 2021 Eq. 52)")
assert close(delta(2, h2, k2), 1.15495, 1e-5) and close(1 + k2 - h2, 0.6941, 1e-4)

# Rigid-Earth peak vertical tide (degree 2, body overhead), mean distances.
R = 6_371_000.0
d_moon, d_sun = 3.844e8, 1.495_978_707e11
g_moon = 2 * GM_MOON * R / d_moon**3
g_sun = 2 * GM_SUN * R / d_sun**3
show("lunar peak tide, rigid, µGal", g_moon / MICROGAL)
show("solar peak tide, rigid, µGal", g_sun / MICROGAL)
show("solar/lunar ratio", g_sun / g_moon, "≈0.46")
assert close(g_sun / g_moon, 0.46, 0.01)
show("lunar tidal gradient 2GM/d^3 at Earth, E", 2 * GM_MOON / d_moon**3 / EOTVOS)

# ---------------------------------------------------------------------------------------------
# 3. Permanent tide (Mäkinen 2021, J Geod 95:106, Eqs. 18, 20, A-4, A-8).
# ---------------------------------------------------------------------------------------------
A = -2.9166  # m^2 s^-2, W_T = A (r/a)^2 (sin^2 phi - 1/3)
show("A' = 2A/3", 2 * A / 3, -1.9444)
show("A'' = A'/sqrt5", 2 * A / 3 / sqrt(5), -0.86956)
assert abs(2 * A / 3 + 1.9444) < 1e-4 and abs(2 * A / 3 / sqrt(5) + 0.86956) < 1e-5
a = 6_378_137.0
gT = lambda phi: -2 * A / a * (sin(phi)**2 - 1 / 3) / MICROGAL  # spherical approximation
show("g_T equator, µGal (spherical)", gT(0.0), -30.49)
show("g_T pole, µGal (spherical)", gT(pi / 2), "60.77 = -30.49+90.95+0.31")
assert abs(gT(0.0) + 30.49) < 0.01 and abs(gT(pi / 2) - 60.77) < 0.3
gamma_e = 9.780_326_771_5  # GRS80 normal gravity at the equator
N_eq = A * (0 - 1 / 3) / gamma_e
show("N_mean - N_zero at equator, m", N_eq, "0.099 (Ekman 1989 via A-8)")
assert abs(N_eq - 0.099) < 0.001
corr_eq = (0.605 - 0.30190) * 99.0
show("Ekman tide-free -> zero-tide height corr at equator, mm", corr_eq, "(h2-k2)·99 mm")

# ---------------------------------------------------------------------------------------------
# 4. W0 and clock-measured potential differences (IAG 2015 Res. 1; IAU 2000 B1.9).
# ---------------------------------------------------------------------------------------------
W0_IHRS = 62_636_853.4
L_G = 6.969_290_134e-10
W0_TT = L_G * C**2
show("W0 implied by L_G (TT definition), m^2/s^2", W0_TT, "62 636 856.0")
show("W0(TT) - W0(IHRS), m^2/s^2", W0_TT - W0_IHRS)
show("  as a height, m", (W0_TT - W0_IHRS) / 9.80)
assert abs(W0_TT - 62_636_856.0) < 0.1
# GRS80 (Moritz 2000) Somigliana normal gravity; gamma_e, k, e^2 must reproduce gamma_p.
GE, KS, E2 = 9.780_326_771_5, 0.001_931_851_353, 0.006_694_380_022_90
somigliana = lambda phi: GE * (1 + KS * sin(phi)**2) / sqrt(1 - E2 * sin(phi)**2)
show("GRS80 normal gravity at the pole, m/s^2", somigliana(pi / 2), 9.8321863685)
assert abs(somigliana(pi / 2) - 9.832_186_368_5) < 1e-9
U0_GRS80 = 62_636_860.850
show("U0(GRS80) - W0(IHRS) as a height, m", (U0_GRS80 - W0_IHRS) / 9.80)
dW = 3918.1  # 457 km Sr link, PR Applied 21, L061001 (2024)
show("3918.1 m^2/s^2 as fractional frequency", dW / C**2)
show("  as a height at g = 9.807, m", dW / 9.807)

# ---------------------------------------------------------------------------------------------
# 5. Equivalence principle: Eötvös parameter and MICROSCOPE coupling bounds.
# ---------------------------------------------------------------------------------------------
eta_quad = sqrt(2.3**2 + 1.5**2)
show("MICROSCOPE eta uncertainty (quadrature), 1e-15", eta_quad, "≈2.7")
g_orbit = GM_EARTH / (6_371e3 + 710e3)**2
show("  differential acceleration it resolves at 710 km, m/s^2", eta_quad * 1e-15 * g_orbit)
show("Asenbaum 2020 eta uncertainty (quadrature), 1e-12", sqrt(1.8**2 + 3.4**2))
ratio = E2_4PIEPS0 / (G * U**2)
show("alpha/(G u^2) (Fayet 2025 Eq. 7)", ratio, "1.2536e36")
assert close(ratio, 1.2536e36, 2e-4)
d_bl = 0.00079 - 0.05704
show("Delta((B-L)/A_r) Ti-Pt = Delta(B/A_r) - Delta(L/A_r)", d_bl, -0.05625)
coef_bl = 3.623e34
show("implied (B-L)/A_r of Earth", coef_bl / (ratio * abs(d_bl)))
for spin, bound, pub in (("spin-1", 4.5e-15, 1.1e-25), ("spin-0", 6.5e-15, 1.3e-25)):
    eps = sqrt(bound / coef_bl)
    g_bl = eps * sqrt(4 * pi * ALPHA)
    show(f"|g_B-L| bound, {spin}", g_bl, pub)
    assert close(g_bl, pub, 0.05)

# ---------------------------------------------------------------------------------------------
# 6. Atom gradiometer field survey (Stray et al. 2022, Nature 602, 590).
# ---------------------------------------------------------------------------------------------
k_eff = 4 * pi / 780.241e-9
T, Lb = 0.085, 1.0
phase = k_eff * 150 * EOTVOS * Lb * T**2
show("phase of 150 E over 1 m at T = 85 ms, mrad", phase * 1e3, 17.5)
assert close(phase * 1e3, 17.5, 0.01)
qpn = sqrt(2) / sqrt(1e5) / (k_eff * Lb * T**2) / sqrt(0.7)
show("QPN gradient ASD, 1e5 atoms, C = 1, 0.7 Hz, E/rtHz", qpn / EOTVOS, "≈44")
assert close(qpn / EOTVOS, 44, 0.1)
show("averaging time to 20 E at 466 E/rtHz, s", (466 / 20)**2, "'within 10 min'")
lam = -1800 * 4.0  # kg/m line density of a 2 m x 2 m void, soil contrast -1.8 g/cm^3
for h0 in (0.0, 0.2):
    z1, z2 = 1.89 + h0, 2.89 + h0
    dg = 2 * G * lam * (1 / z1 - 1 / z2)
    show(f"line-mass Δg over 1 m baseline, lower sensor {h0} m up, E", dg / EOTVOS, "≈150 (modelled peak)")

# ---------------------------------------------------------------------------------------------
# 7. Kerr ringdown (Berti, Cardoso & Will 2006, Table VIII; values Table II).
# ---------------------------------------------------------------------------------------------
QNM = {  # (l, m, n): (f1, f2, f3, q1, q2, q3)
    (2, 2, 0): (1.5251, -1.1568, 0.1292, 0.7000, 1.4187, -0.4990),
    (2, 2, 1): (1.3673, -1.0260, 0.1628, 0.1000, 0.5436, -0.4731),
    (2, 1, 0): (0.6000, -0.2339, 0.4175, -0.3000, 2.3561, -0.2277),
    (3, 3, 0): (1.8956, -1.3043, 0.1818, 0.9000, 2.3430, -0.4810),
    (4, 4, 0): (2.3000, -1.5056, 0.2244, 1.1929, 3.1191, -0.4825),
}


def qnm(lmn, j):
    f1, f2, f3, q1, q2, q3 = QNM[lmn]
    F = f1 + f2 * (1 - j)**f3
    Q = q1 + q2 * (1 - j)**q3
    return F, F / (2 * Q), Q  # M omega_R, M omega_I, quality factor


MAX_ERR = {(2, 2, 0): (0.0185, 0.0088), (2, 2, 1): (0.0156, 0.0169)}
table2 = {((2, 2, 0), 0.0): (0.3737, 0.0890), ((2, 2, 0), 0.7): (0.5326, 0.0808),
          ((2, 2, 0), 0.9): (0.6716, 0.0649), ((2, 2, 1), 0.0): (0.3467, 0.2739),
          ((2, 2, 1), 0.7): (0.5212, 0.2442)}
for (lmn, j), (wr, wi) in table2.items():
    F, WI, _ = qnm(lmn, j)
    show(f"BCW fit {lmn} j={j}: M omega_R", F, wr)
    show(f"BCW fit {lmn} j={j}: M omega_I", WI, wi)
    fmax, qmax = MAX_ERR[lmn]  # the paper's maximum fit errors over j in [0, 0.99]
    assert close(F, wr, fmax + 0.001) and close(F / (2 * WI), wr / (2 * wi), qmax + 0.002)


def ringdown(m_det_solar, j, lmn=(2, 2, 0)):
    F, _, Q = qnm(lmn, j)
    f = F / (2 * pi * m_det_solar * T_SUN)
    return f, Q / (pi * f)


f, tau = ringdown(68.0, 0.67)
show("GW150914 f_220 from M_f = 68 (detector), a_f = 0.67, Hz", f, "251 ± 8")
show("GW150914 tau_220, ms", tau * 1e3, "4.0 ± 0.3")
assert abs(f - 251) < 8 and abs(tau * 1e3 - 4.0) < 0.3
f, tau = ringdown(62.7 * 1.09, 0.68)
show("GW250114 f_220 from 62.7 M_sun (1+z = 1.09), 0.68, Hz", f, "not found")
show("GW250114 tau_220, ms", tau * 1e3)

# ---------------------------------------------------------------------------------------------
# 8. Remnant spin (Hofmann, Barausse & Rezzolla 2016) and mass (Barausse, Morozova & Rezzolla 2012).
# ---------------------------------------------------------------------------------------------
def r_isco(a):
    z1 = 1 + (1 - a * a)**(1 / 3) * ((1 + a)**(1 / 3) + (1 - a)**(1 / 3))
    z2 = sqrt(3 * a * a + z1 * z1)
    sgn = (a > 0) - (a < 0)
    return 3 + z2 - sgn * sqrt((3 - z1) * (3 + z1 + 2 * z2))


def e_isco(a):
    return sqrt(1 - 2 / (3 * r_isco(a)))


def l_isco(a):
    return 2 / (3 * sqrt(3)) * (1 + 2 * sqrt(3 * r_isco(a) - 2))


K = {(0, 1): 3.39221, (0, 2): 4.48865, (0, 3): -5.77101, (0, 4): -13.0459,
     (1, 0): 35.1278, (1, 1): -72.9336, (1, 2): -86.0036, (1, 3): 93.7371, (1, 4): 200.975,
     (2, 0): -146.822, (2, 1): 387.184, (2, 2): 447.009, (2, 3): -467.383, (2, 4): -884.339,
     (3, 0): 223.911, (3, 1): -648.502, (3, 2): -697.177, (3, 3): 753.738, (3, 4): 1166.89}
XI = 0.474046
# k00 is fixed by Eq. 11: sqrt(3)/2 + sum_i k_i0 / 4^(2+i) = 0.68646 (Scheel et al. 2009).
K[(0, 0)] = 16 * (0.68646 - sqrt(3) / 2 - sum(K[(i, 0)] / 4**(2 + i) for i in (1, 2, 3)))
show("HBR16 k00 (nM=3, nJ=4)", K[(0, 0)], "≈ -5.9 (footnote 7)")


def final_spin(m1, m2, a1, a2):
    """Aligned spins; a1 belongs to the heavier body m1."""
    q = m2 / m1
    nu = q / (1 + q)**2
    atot = (a1 + a2 * q * q) / (1 + q)**2
    aeff = atot + XI * nu * (a1 + a2)
    s = sum(k * nu**(2 + i) * aeff**j for (i, j), k in K.items())
    return atot + nu * (l_isco(aeff) - 2 * atot * (e_isco(aeff) - 1)) + s


def radiated(m1, m2, a1, a2, p0=0.04827, p1=0.01707):
    q = m2 / m1
    nu = q / (1 + q)**2
    at = (a1 + q * q * a2) / (1 + q)**2
    E = e_isco(at)
    return (1 - E) * nu + 4 * nu * nu * (4 * p0 + 16 * p1 * at * (at + 1) + E - 1)


show("HBR16 a_f, q = 1, no spin", final_spin(1, 1, 0, 0), 0.68646)
assert abs(final_spin(1, 1, 0, 0) - 0.68646) < 1e-9
show("BMR12 M_f/M, q = 1, no spin", 1 - radiated(1, 1, 0, 0), "0.95173 by construction (p0)")
af = final_spin(33.6, 32.2, 0, 0)
mf = (33.6 + 32.2) * (1 - radiated(33.6, 32.2, 0, 0))
show("GW250114 a_f from 33.6 + 32.2 M_sun, no spin", af, "0.68 ± 0.01")
show("GW250114 M_f, M_sun", mf, "62.7 (+1.0 -1.1)")
assert abs(af - 0.68) < 0.01 and abs(mf - 62.7) < 1.1
show("test-particle limit a_f(q=1e-4, a1=0.9) ≈ a1", final_spin(1, 1e-4, 0.9, 0))


# ---------------------------------------------------------------------------------------------
# 9. Hawking area law with Kerr areas.
# ---------------------------------------------------------------------------------------------
def kerr_area(m_solar, chi):
    rg = m_solar * L_SUN
    return 8 * pi * rg * rg * (1 + sqrt(1 - chi * chi))


Ai = kerr_area(33.6, 0) + kerr_area(32.2, 0)
Af = kerr_area(62.7, 0.68)
show("GW250114 initial area (spins 0), km^2 source frame", Ai / 1e6, "≈240 000 (Ars Technica)")
show("GW250114 final area, km^2 source frame", Af / 1e6, "≈400 000 (Ars Technica)")
show("A_f / A_i", Af / Ai)
assert close(Ai / 1e6, 240_000, 0.02)
q1 = kerr_area(0.9516, 0.6865) / (2 * kerr_area(0.5, 0))
show("A_f/A_i, equal-mass non-spinning (M_f 0.9516, a_f 0.6865)", q1)

# ---------------------------------------------------------------------------------------------
# 10. Tidal deformability (Flanagan & Hinderer 2008; GW170817 PRX 9, 011001).
# ---------------------------------------------------------------------------------------------
def lambda_tilde(m1, m2, l1, l2):
    M = m1 + m2
    return 16 / 13 * ((m1 + 12 * m2) * m1**4 * l1 + (m2 + 12 * m1) * m2**4 * l2) / M**5


assert abs(lambda_tilde(1.4, 1.4, 190, 190) - 190) < 1e-9
show("Λ̃(1.46, 1.27, 200, 450)", lambda_tilde(1.46, 1.27, 200, 450))
show("GW170817 Λ̃ HPD upper edge 300 + 420", 300 + 420, "720")
k2_n1 = (15 - pi**2) / (2 * pi**2)
show("Newtonian k2 of an n = 1 polytrope (15 - pi^2)/(2 pi^2)", k2_n1, "≈0.26")
comp = 1.4 * L_SUN / 11.9e3
show("compactness of 1.4 M_sun, 11.9 km", comp)
show("k2 that gives Λ = 190 at that compactness", 1.5 * 190 * comp**5)
print("\nall checks passed")
