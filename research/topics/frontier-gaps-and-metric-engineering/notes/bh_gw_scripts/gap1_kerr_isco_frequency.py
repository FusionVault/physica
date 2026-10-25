"""Gap 1: Kerr ISCO orbital frequency M Omega = 1/(r^{3/2} + a) (Isoyama et al. 2014, arXiv:1404.6133,
Eq. 3, after Bardeen, Press & Teukolsky 1972), signed spin a, and checks:
  * exact limits a = 0 (6^{-3/2}), a = +1 (1/2), a = -1 (1/26)
  * Isoyama's redshift z0(Omega)^2 = (1 - a Omega)[1 + a Omega - 3 Omega^{2/3} (1 - a Omega)^{1/3}]
    has z0''(Omega_isco) = 0 (their ISCO condition), and E0 = z0 - Omega z0' has E0'(Omega_isco) = 0
  * Isoyama Table I C_Omega gives the O(eta) conservative self-force shift (M+mu) Omega = M Omega0 (1 + eta C)
"""
import numpy as np

G, c, Msun = 6.67430e-11, 299792458.0, 1.98847e30


def r_isco(a):
    z1 = 1 + np.cbrt(1 - a * a) * (np.cbrt(1 + a) + np.cbrt(1 - a))
    z2 = np.sqrt(3 * a * a + z1 * z1)
    return 3 + z2 - np.sign(a) * np.sqrt((3 - z1) * (3 + z1 + 2 * z2))


def m_omega(a):
    return 1 / (r_isco(a) ** 1.5 + a)


for a, exact in [(0.0, 6 ** -1.5), (1.0, 0.5), (-1.0, 1 / 26)]:
    print(f"a={a:+}: M Omega = {m_omega(a):.15f}, exact {exact:.15f}")


def z0(W, a):
    return np.sqrt((1 - a * W) * (1 + a * W - 3 * W ** (2 / 3) * (1 - a * W) ** (1 / 3)))


for a in [-0.9, -0.5, 0.0, 0.5, 0.9, 0.99]:
    W = m_omega(a); h = 1e-4 * W
    d2 = (z0(W + h, a) - 2 * z0(W, a) + z0(W - h, a)) / h**2
    E = lambda w: z0(w, a) - w * (z0(w + 1e-7, a) - z0(w - 1e-7, a)) / 2e-7
    dE = (E(W + h) - E(W - h)) / (2 * h)
    r = r_isco(a)
    E_bpt = np.sqrt(1 - 2 / (3 * r))
    print(f"a={a:+.2f}: r={r:.10f} M Omega={W:.12f} z0''={d2:+.1e} E0'={dE:+.1e} E0(Omega)={E(W):.10f} vs sqrt(1-2/(3r))={E_bpt:.10f}")

print("\nGW frequency f = Omega/pi for M = 10 Msun (Hz):")
for a in [-1.0, -0.9, 0.0, 0.5, 0.7, 0.9, 0.99, 0.998, 1.0]:
    f = m_omega(a) * c**3 / (np.pi * G * 10 * Msun)
    print(f"  a={a:+.3f}: r_isco={r_isco(a):.6f}  M Omega={m_omega(a):.8f}  f_GW={f:9.3f} Hz  ratio to Schwarzschild {m_omega(a)/6**-1.5:.5f}")

# Isoyama Table I (q = spin) and C_Omega(0) from text
C = {0.0: 1.25101539, 0.1: 1.245568, 0.3: 1.239927, 0.5: 1.249234, 0.7: 1.29319, 0.9: 1.38157, -0.5: 1.286093, -0.9: 1.314016}
print("\nconservative GSF shift (M+mu)Omega/(M Omega0) - 1 = eta C_Omega, eta = mu/M:")
for a, Cw in C.items():
    print(f"  a={a:+.1f}: C={Cw}; eta=1e-3 -> {1e-3*Cw:.2e}; eta=0.1 -> {0.1*Cw:.3f}")
