"""Gap 3: HBR16 generic-spin final spin (arXiv:1605.01938 Eqs. ell, atot, aeff, afinnorm, angle_map,
eq:dir) and BMR12 generic radiated energy (arXiv:1206.3803 Eq. final_eq1 with a~ of the generic
section); Kerr surface gravity (Townsend gr-qc/9707012 ch. 4 and exercise) and the first-law check.

Coefficients: HBR16 Table I top block (n_M = 1, n_J = 2), k00 from the equal-mass constraint, as the
library's remnant.rs does; BMR12 p0, p1."""
import numpy as np

K01, K02, K10, K11, K12, XI = -1.2019, -1.20764, 3.79245, 1.18385, 4.90494, 0.41616
EQ = 0.68646
K00 = 16 * (EQ - np.sqrt(3) / 2) - K10 / 4
P0, P1 = 0.04827, 0.01707


def isco(a):
    z1 = 1 + np.cbrt(1 - a * a) * (np.cbrt(1 + a) + np.cbrt(1 - a))
    z2 = np.sqrt(3 * a * a + z1 * z1)
    r = 3 + z2 - np.sign(a) * np.sqrt((3 - z1) * (3 + z1 + 2 * z2))
    return r, np.sqrt(1 - 2 / (3 * r)), 2 / (3 * np.sqrt(3)) * (1 + 2 * np.sqrt(3 * r - 2))


def final_spin_aligned(m1, m2, a1, a2):   # the library's Eq. 7 (remnant.rs)
    m = m1 + m2; nu = m1 * m2 / m**2; at = (a1 * m1**2 + a2 * m2**2) / m**2
    ae = at + XI * nu * (a1 + a2)
    _, e, l = isco(ae)
    fit = nu**2 * (K00 + K01 * ae + K02 * ae**2) + nu**3 * (K10 + K11 * ae + K12 * ae**2)
    return at + nu * (l - 2 * at * (e - 1)) + fit


def final_spin_generic(m1, m2, A1, A2, beta, gamma, alpha, eps=0.0):
    """|a_fin| for spin magnitudes A1, A2 at angles beta, gamma to L and alpha between them.
    q = M2/M1 <= 1 is required by HBR16's convention; eps is the angle-map epsilon (0.024)."""
    if m2 > m1:
        m1, m2, A1, A2, beta, gamma = m2, m1, A2, A1, gamma, beta
    amap = lambda t, e: 2 * np.arctan((1 + e) * np.tan(t / 2)) if e else t
    beta, gamma = amap(beta, eps), amap(gamma, eps)
    q = m2 / m1; nu = q / (1 + q) ** 2
    at = (A1 * np.cos(beta) + A2 * np.cos(gamma) * q * q) / (1 + q) ** 2
    ae = at + XI * nu * (A1 * np.cos(beta) + A2 * np.cos(gamma))
    _, e, l = isco(ae)
    ell = abs(l - 2 * at * (e - 1) + nu * (K00 + K01 * ae + K02 * ae**2) + nu**2 * (K10 + K11 * ae + K12 * ae**2))
    s = A1**2 + A2**2 * q**4 + 2 * A1 * A2 * q * q * np.cos(alpha) + 2 * (A1 * np.cos(beta) + A2 * q * q * np.cos(gamma)) * ell * q + ell**2 * q * q
    return np.sqrt(s) / (1 + q) ** 2, ell


def erad(m1, m2, A1, A2, beta, gamma):
    m = m1 + m2; nu = m1 * m2 / m**2
    a = (A1 * np.cos(beta) * m1**2 + A2 * np.cos(gamma) * m2**2) / m**2
    _, e, _ = isco(a)
    return (1 - e) * nu + 4 * nu * nu * (4 * P0 + 16 * P1 * a * (a + 1) + e - 1)


print("k00 =", K00)
# 1. aligned limit
worst = 0
for (m1, m2, a1, a2) in [(1, 1, 0, 0), (36, 29, 0.3, -0.2), (1, 0.25, 0.9, 0.5), (1, 0.5, -0.8, 0.6), (1, 1, 0.99, 0.99)]:
    ali = final_spin_aligned(m1, m2, a1, a2)
    gen, _ = final_spin_generic(m1, m2, abs(a1), abs(a2), 0 if a1 >= 0 else np.pi, 0 if a2 >= 0 else np.pi,
                                0 if a1 * a2 >= 0 else np.pi)
    worst = max(worst, abs(abs(ali) - gen))
print(f"1. generic formula at beta, gamma in {{0, pi}} vs aligned Eq. 7: max |diff| = {worst:.1e}")
print("2. equal mass, no spin:", final_spin_generic(1, 1, 0, 0, 0, 0, 0)[0], "(0.68646)")
# 3. in-plane spins: the spins add in quadrature with ell*nu along L (beta = gamma = pi/2)
for A in [0.5, 0.9]:
    g, ell = final_spin_generic(1, 1, A, A, np.pi / 2, np.pi / 2, 0.0)
    print(f"3. q=1, both spins {A} in-plane and parallel: |a_fin| = {g:.6f} = sqrt((A/2)^2 + (ell nu)^2) = {np.hypot(A/2, ell/4):.6f}")
    g2, _ = final_spin_generic(1, 1, A, A, np.pi / 2, np.pi / 2, np.pi)
    print(f"   ... and antiparallel in-plane: |a_fin| = {g2:.6f} = ell nu = {ell/4:.6f}")
# 4. epsilon map
print("4. angle map with eps=0.024 at beta=60deg:", np.degrees(2*np.arctan(1.024*np.tan(np.radians(30)))), "deg; small-eps form", 60 + np.degrees(0.024*np.sin(np.radians(60))))
# 5. BMR12 generic: projections enter only through a~
print("5. BMR12 E_rad q=1 a=0.6 at beta=gamma=60deg equals aligned a=0.3:", erad(1, 1, 0.6, 0.6, np.pi/3, np.pi/3), erad(1, 1, 0.3, 0.3, 0, 0))
print("   BMR12 q=1 aligned maximal: ", erad(1, 1, 1, 1, 0, 0), " (p0 + 3 p1 =", P0 + 3 * P1, ")")

# 6. Kerr surface gravity, G = c = M = 1: kappa = sqrt(1-chi^2)/(2(1+sqrt(1-chi^2)))
kap = lambda chi: np.sqrt(1 - chi * chi) / (2 * (1 + np.sqrt(1 - chi * chi)))
for chi in [0.0, 0.6, 0.8, 1.0]:
    print(f"6. chi={chi}: kappa/kappa_Schw = {kap(chi)/0.25:.12f}")
# Townsend's two forms agree: (r+ - r-)/(2(r+^2+a^2)) and sqrt(M^4-J^2)/(2M(M^2+sqrt(M^4-J^2)))
chi = 0.7; rp, rm = 1 + np.sqrt(1 - chi**2), 1 - np.sqrt(1 - chi**2)
print("   Townsend forms:", (rp - rm) / (2 * (rp**2 + chi**2)), np.sqrt(1 - chi**2) / (2 * (1 + np.sqrt(1 - chi**2))), kap(chi))
# first law: dM = kappa/(8 pi) dA + Omega_H dJ with A = 8 pi (M^2 + sqrt(M^4 - J^2)); check dM/dA at fixed J
def M_of_AJ(A, J):
    from scipy.optimize import brentq
    return brentq(lambda M: 8 * np.pi * (M * M + np.sqrt(max(M**4 - J * J, 0))) - A, np.sqrt(J) * (1 + 1e-12), 10)
J = 0.7; A0 = 8 * np.pi * (1 + np.sqrt(1 - J * J)); h = 1e-6
dMdA = (M_of_AJ(A0 + h, J) - M_of_AJ(A0 - h, J)) / (2 * h)
print(f"   first law at chi=0.7: dM/dA|_J = {dMdA:.10f}, kappa/(8 pi) = {kap(0.7)/(8*np.pi):.10f}")
# SI test value: T(1 Msun, chi = 0.8) = 0.75 * T_Schw
hbar, kB, G, c, Msun = 1.054571817e-34, 1.380649e-23, 6.67430e-11, 299792458.0, 1.98847e30
TS = hbar * c**3 / (8 * np.pi * G * Msun * kB)
print(f"   T_Schw(1 Msun) = {TS:.6e} K; T(chi=0.8) = {0.75*TS:.6e} K")
