"""Gap 1: demagnetising factors and induced-plus-remanent magnetisation.

Checks
- spheroid factors (Takahashi & Oliveira 2017 GMD eqs 30-33) against the ellipsoid column N of
  Chen, Brug & Goldfarb 1991 (IEEE Trans. Magn. 27, 3601) Table I;
- cylinder magnetometric N_m(chi=0) from Chen et al. eqs (6), (7), (15) against Table I;
- the triaxial factors (GMD eqs 25-27, as coded in Harmonica) for the Warrego model of GMD Table 2,
  chi_max = eps / n_max (GMD eq 60) with eps = 8 % -> 0.116 SI;
- the confocal ellipsoids of GMD Table 1 (eq 64): chi' ~ 0.014, volume ratio ~ 79, chi ratio ~ 85;
- the high-susceptibility sphere M = chi/(1 + chi/3) H0.
"""
import numpy as np
from scipy.special import ellipk, ellipe, ellipkinc, ellipeinc
from scipy.integrate import quad

def prolate(m):  # a/b = m > 1, factor along the long axis (GMD eq 30)
    s = np.sqrt(m * m - 1)
    n11 = (m / s * np.log(m + s) - 1) / (m * m - 1)
    return n11, (1 - n11) / 2

def oblate_axial(q):  # q = c/a < 1 along the short (symmetry) axis (Harmonica's form of GMD eq 32)
    return (1 - q / np.sqrt(1 - q * q) * np.arccos(q)) / (1 - q * q)

def landau(a, b, c):  # Prozorov & Kogan 2018 eq (Nellipsoid): factor along c, any ratios
    k = max(a, b, c); a, b, c = a / k, b / k, c / k  # scale-free; quad needs order-1 numbers
    f = lambda s: 0.5 * a * b * c / ((s + c * c) * np.sqrt((s + a * a) * (s + b * b) * (s + c * c)))
    return quad(f, 0, np.inf, limit=400, epsabs=1e-14, epsrel=1e-13)[0]

def triaxial(a, b, c):  # a > b > c, Harmonica _demag_tensor_triaxial_internal (GMD eqs 25-27)
    phi, k = np.arccos(c / a), (a * a - b * b) / (a * a - c * c)
    F, E = ellipkinc(phi, k), ellipeinc(phi, k)
    r = np.sqrt(a * a - c * c)
    nxx = a * b * c / (r * (a * a - b * b)) * (F - E)
    nzz = -a * b * c / (r * (b * b - c * c)) * E + b * b / (b * b - c * c)
    return nxx, 1 - nxx - nzz, nzz

def spheroid_axial(gamma):  # axis 2l, transverse 2a, gamma = l/a: Chen Table I column N
    if gamma > 1: return prolate(gamma)[0]
    if gamma < 1: return oblate_axial(gamma)
    return 1 / 3

def cylinder_nm(gamma):  # Chen 1991 eqs (6), (7), (15) with n = 1, a = 1, l = gamma; mu0 cancels
    a, l = 1.0, gamma
    m = a * a / (a * a + l * l)  # k_s^2
    Ls = 2 / (3 * l * l) * (np.sqrt(a * a + l * l) * (l * l * ellipk(m) + (a * a - l * l) * ellipe(m)) - a ** 3)
    return 1 - 2 * l * Ls / (np.pi * a * a)

# Chen, Brug & Goldfarb 1991, Table I: gamma, N_m(0), N (ellipsoid)
table = [(0.01, 0.9650, 0.9845), (0.1, 0.7967, 0.8608), (0.5, 0.4745, 0.5272), (1.0, 0.3116, 0.3333),
         (2.0, 0.1819, 0.1736), (5.0, 0.07991, 0.05582), (10.0, 0.04119, 0.02029), (100.0, 0.004232, 0.0004299)]
print("gamma   N_spheroid(ours)  Chen N   | N_m cyl(ours)  Chen N_m | Landau integral")
for g, nm, ne in table:
    ns, nc = spheroid_axial(g), cylinder_nm(g)
    nl = landau(1.0, 1.0, g)
    print(f"{g:7.2f}  {ns:.6f}  {ne:<9} | {nc:.6f}  {nm:<9} | {nl:.6f}")
    assert abs(ns - ne) < 0.6 * 10 ** -(len(str(ne).split('.')[1]) if ne < 1 else 4) + 1e-12, (g, ns, ne)
    assert abs(nl - ns) < 1e-10
    assert abs(nc - nm) <= 0.5 * 10 ** -(len(str(nm).split('.')[1])) + 1e-12, (g, nc, nm)

# Triaxial: sum rule and limits
for abc in [(3.0, 2.0, 1.0), (1000.0, 700.0, 200.0), (490.7, 69.7, 30.0)]:
    n = triaxial(*abc)
    nl = [landau(abc[1], abc[2], abc[0]), landau(abc[0], abc[2], abc[1]), landau(*abc)]
    print("triaxial", abc, ["%.6f" % x for x in n], "sum", sum(n), "landau", ["%.6f" % x for x in nl])
    assert max(abs(x - y) for x, y in zip(n, nl)) < 1e-9
n_max = triaxial(490.7, 69.7, 30.0)[2]
print("Warrego n33 =", n_max, " chi_max(eps=8%) =", 0.08 / n_max, "(GMD: 0.116 SI)")

# Confocal pair, GMD Table 1 / eq 64, field along a
a, b, c, chi, u = 900.0, 500.0, 100.0, 1.2, 2e6
a2, b2, c2 = np.sqrt(a * a + u), np.sqrt(b * b + u), np.sqrt(c * c + u)
v1, v2 = 4 / 3 * np.pi * a * b * c, 4 / 3 * np.pi * a2 * b2 * c2
H0 = 1.0
P = v1 * chi * H0 / (1 + chi * triaxial(a, b, c)[0])
chi2 = P / (v2 * H0 - triaxial(a2, b2, c2)[0] * P)
print(f"confocal: a'={a2:.2f} b'={b2:.2f} c'={c2:.2f}  chi'={chi2:.5f}  V2/V1={v2/v1:.2f}  chi/chi'={chi/chi2:.1f}")

# High-susceptibility sphere: M = chi H0 / (1 + chi/3), so the anomaly scales by 1/(1 + chi/3)
for chi in [0.01, 0.1, 1.0, 1.69, 10.0]:
    print(f"sphere chi={chi}: M/(chi H0) = {1/(1+chi/3):.6f}  (error of ignoring demag {chi/3*100:.3f} %)")

# One formula for every ellipsoid, any axis order, sphere and spheroids included (no branches):
# N_c = (abc/3) R_D(a^2, b^2, c^2), from the Landau integral with Carlson's R_D (geopotential::special
# already has carlson_rd). Checked against the Landau integral and Chen's Table I.
from scipy.special import elliprd
def carlson_factor(a, b, c):
    k = max(a, b, c); a, b, c = a / k, b / k, c / k
    return a * b * c / 3 * elliprd(a * a, b * b, c * c)
worst = 0.0
for abc in [(3.0, 2.0, 1.0), (490.7, 69.7, 30.0), (1.0, 1.0, 1.0), (1.0, 1.0, 10.0), (1.0, 1.0, 0.1), (900.0, 500.0, 100.0)]:
    a, b, c = abc
    n = [carlson_factor(b, c, a), carlson_factor(c, a, b), carlson_factor(a, b, c)]
    worst = max(worst, abs(sum(n) - 1), abs(n[2] - landau(*abc)))
for g, nm, ne in table:
    worst_chen = abs(carlson_factor(1.0, 1.0, g) - ne)
    assert worst_chen <= 0.5 * 10 ** -(len(str(ne).split('.')[1])) + 1e-12
print(f"Carlson R_D form: max |sum - 1| or |N - Landau| = {worst:.1e}; all Chen Table I spheroid N reproduced to the printed digits")
print(f"induced magnetisation for chi = 0.01 in F = 50 000 nT: chi F/mu0 = {0.01 * 50e-6 / (4e-7 * np.pi):.4f} A/m")
