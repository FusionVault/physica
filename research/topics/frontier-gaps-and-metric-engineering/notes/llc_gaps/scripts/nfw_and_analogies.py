"""NFW normalisation (Wright & Brainerd 2000), concentration-mass relations (Duffy et al. 2008 Table 1;
Dutton & Maccio 2014 Eqs. 7-10), Sigma_crit, and the cross-crate analogy identities.
Constants: CODATA 2018 / IAU 2015 nominal, as in metrology."""
import math
import sys
import numpy as np
from scipy.integrate import quad
from scipy.special import ellipk

sys.path.insert(0, ".")
from cosmology_checks import Cosmo  # noqa: E402
from finite_source import B0  # noqa: E402

G, C = 6.67430e-11, 299792458.0
HBAR, KB = 1.054571817e-34, 1.380649e-23
GM_SUN = 1.3271244e20
M_SUN = GM_SUN / G
AU = 149597870700.0
PC = AU * 648000 / math.pi
MPC = 1e6 * PC
R_SUN = 6.957e8  # IAU 2015 nominal solar radius
MAS = math.pi / 180 / 3600 * 1e-3
UAS = MAS * 1e-3


def delta_c(c):
    return 200 / 3 * c ** 3 / (math.log(1 + c) - c / (1 + c))


def duffy_c200(M_h, z, A=5.71, B=-0.084, Cz=-0.47):
    """Duffy et al. 2008 Table 1, NFW, full sample, z = 0-2, M200 in h^-1 Msun, pivot 2e12."""
    return A * (M_h / 2e12) ** B * (1 + z) ** Cz


def dutton_c200(M_h, z):
    """Dutton & Maccio 2014 Eqs. 7-8 with log10 c = a + b log10(M/1e12 h^-1 Msun)."""
    b = -0.101 + 0.026 * z
    a = 0.520 + (0.905 - 0.520) * math.exp(-0.617 * z ** 1.21)
    return 10 ** (a + b * math.log10(M_h / 1e12))


def dutton_cvir(M_h, z):
    b = -0.097 + 0.024 * z
    a = 0.537 + (1.025 - 0.537) * math.exp(-0.718 * z ** 1.08)
    return 10 ** (a + b * math.log10(M_h / 1e12))


if __name__ == "__main__":
    print("== delta_c(c) = (200/3) c^3/(ln(1+c) - c/(1+c))  (Wright & Brainerd Eq. 2)")
    for c in [4.0, 5.0, 10.0]:
        print(f"  c={c}: delta_c = {delta_c(c):.6f}")
    print("== concentration-mass relations (M in h^-1 Msun)")
    for M in [1e12, 2e12, 1e14, 1e15]:
        for z in [0.0, 0.5, 1.0]:
            print(f"  M={M:.0e} z={z}: Duffy c200 = {duffy_c200(M, z):.5f}  Dutton c200 = {dutton_c200(M, z):.5f}  Dutton cvir = {dutton_cvir(M, z):.5f}")
    print("== an NFW halo in Planck18: M200 = 1e14 Msun/h at z_d = 0.3, source at z_s = 1.0")
    p = Cosmo(67.66, 0.30966, Tcmb0=2.7255, Neff=3.046, m_nu=(0, 0, 0.06))
    h = 0.6766
    zd, zs = 0.3, 1.0
    Hz = p.H0 * p.E(zd) * 1e3 / MPC
    rho_c = 3 * Hz ** 2 / (8 * math.pi * G)
    M200 = 1e14 / h * M_SUN
    r200 = (3 * M200 / (800 * math.pi * rho_c)) ** (1 / 3)
    c = dutton_c200(1e14, zd)
    rs = r200 / c
    dc = delta_c(c)
    # mass identity: 4 pi rho_c delta_c rs^3 [ln(1+c) - c/(1+c)] = M200
    m_check = 4 * math.pi * rho_c * dc * rs ** 3 * (math.log(1 + c) - c / (1 + c))
    dd, ds, dds = p.DA(zd) * MPC, p.DA(zs) * MPC, p.DA(zd, zs) * MPC
    sig_c = C ** 2 / (4 * math.pi * G) * ds / (dd * dds)
    print(f"  rho_c(z_d) = {rho_c:.6e} kg/m^3, r200 = {r200/MPC*1e3:.4f} kpc, c200 = {c:.4f}, delta_c = {dc:.3f}, r_s = {rs/MPC*1e3:.3f} kpc")
    print(f"  mass identity M(<r200)/M200 - 1 = {m_check/M200 - 1:.1e}")
    print(f"  Sigma_crit = {sig_c:.6f} kg/m^2 = {sig_c / M_SUN * PC**2:.3f} Msun/pc^2; kappa_s = 2 r_s delta_c rho_c/Sigma_crit = {2*rs*dc*rho_c/sig_c:.5f}")

    print("\n== 41j one lensing constant")
    kappa = 4 * G / (C * C * AU)
    print(f"  kappa M_sun = {kappa*M_SUN/MAS:.6f} mas")
    x = GM_SUN / (C * C * AU)
    lb = 4 * x + 15 * math.pi / 4 * x * x
    print(f"  light_bending(M_sun, 1 AU)/(kappa M_sun) - 1 = {lb/(kappa*M_SUN)-1:.6e}; (15 pi/16) GM/(c^2 AU) = {15*math.pi/16*x:.6e}")
    print(f"  solar limb: kappa M_sun AU/R_sun = {kappa*M_SUN*AU/R_SUN/MAS/1e3:.5f} arcsec (Crosta & Mignard Table 1: 1.75\")")
    # Einstein ring: theta = alpha(theta D_l) D_ls/D_s with Euclidean Galactic distances
    Dl, Ds, M = 4e3 * PC, 8e3 * PC, 0.3 * M_SUN
    pi_rel = AU * (1 / Dl - 1 / Ds)
    # the point-lens Einstein radius from the first-order deflection 4GM/(c^2 b) at b = theta D_l
    th_E = math.sqrt(4 * G * M / C ** 2 * (Ds - Dl) / (Dl * Ds))
    print(f"  theta_E = {th_E/MAS:.6f} mas; sqrt(kappa M pi_rel) = {math.sqrt(kappa*M*pi_rel)/MAS:.6f} mas; ratio-1 = {th_E/math.sqrt(kappa*M*pi_rel)-1:.1e}")
    print(f"  mass in the ring: pi (theta_E D_l)^2 Sigma_crit / M - 1 = {math.pi*(th_E*Dl)**2*C**2/(4*math.pi*G)*Ds/(Dl*(Ds-Dl))/M - 1:.1e}")

    print("\n== Earth's quadrupole deflection (Crosta & Mignard Table 1: Earth 574 uas, Q 0.6 uas)")
    GM_E, R_E, C20bar = 3.986004415e14, 6378136.3, -4.84169495e-4  # analogies.rs values
    J2E = -C20bar * math.sqrt(5)
    mono = 4 * GM_E / (C * C * R_E)
    print(f"  J2 = {J2E:.9e}; monopole at the limb {mono/UAS:.2f} uas; quadrupole {mono*J2E/UAS:.4f} uas")

    print("\n== finite-source B0 as the in-plane potential of a uniform disc (unit mass, unit radius, G = 1)")
    def disc_potential(z):
        # <1/d> over the unit disc: (1/pi) int_0^1 r dr int_0^2pi dth / sqrt(r^2 + z^2 - 2 r z cos th)
        # inner angular integral is the ring kernel 4 K(m)/(r+z) with m = 4rz/(r+z)^2
        f = lambda r: r * 4 * ellipk(4 * r * z / (r + z) ** 2) / (r + z) / math.pi
        pts = [z] if 0 < z < 1 else None
        return quad(f, 0, 1, points=pts, epsabs=0, epsrel=1e-12, limit=400)[0]
    for z in [0.25, 0.5, 2.0, 5.0]:
        print(f"  z={z}: B0(z) = {B0(z):.12f}; z * <1/d>_disc (ring-kernel integral) = {z*disc_potential(z):.12f}")
    print(f"  rim: B0(1) = 4/pi = {4/math.pi:.15f} (the rim potential of a uniform disc 4GM/(pi a)); centre: B0(z)/z -> 2 (2GM/a)")

    print("\n== 41i one horizon temperature")
    T = lambda kap: HBAR * kap / (2 * math.pi * KB)
    Th = lambda m: T(C ** 3 / (4 * G * m))
    print(f"  Hawking T(1 Msun) = {Th(M_SUN):.6e} K; T*M = hbar c^3/(8 pi G k_B) = {HBAR*C**3/(8*math.pi*G*KB)/M_SUN:.6e} K Msun")
    print(f"  Unruh: 1 K needs a = {2*math.pi*KB*C/HBAR:.6e} m/s^2")
    H0 = 67.66e3 / MPC
    print(f"  Gibbons-Hawking at the Planck18 H0 = 67.66: {T(H0):.6e} K")
    # first law: T_H = c^2 (dS/dM)^-1 with S = k_B c^3 A/(4 G hbar), A = 16 pi G^2 M^2/c^4
    S = lambda m: KB * C ** 3 * 16 * math.pi * G * G * m * m / C ** 4 / (4 * G * HBAR)
    m = 10 * M_SUN
    dSdM = (S(m * (1 + 1e-6)) - S(m * (1 - 1e-6))) / (2e-6 * m)
    print(f"  first law: c^2/(dS/dM) / T_H - 1 = {C*C/dSdM/Th(m) - 1:.1e}")
    # Schwarzschild surface gravity as Unruh acceleration at the horizon: a = c^4/(4GM)
    print(f"  T_H(M) / T_U(c^4/(4GM)) - 1 = {Th(m) / T(C**4/(4*G*m)/C) - 1:.1e}")
