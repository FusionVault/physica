"""Finite-source point-lens magnification: Lee et al. 2009 (arXiv:0901.1316) Eq. A_fs (exact, uniform
disc), Yoo et al. 2004 (astro-ph/0309302) B0(z) = (4z/pi) E(arcsin min(1, 1/z) | modulus z), Q(z),
and limb-darkened B1 by direct 2D integration of Yoo's W_1.
"""
import math
import numpy as np
from scipy.special import elliprf, elliprd, ellipk, ellipe
from scipy.integrate import quad, dblquad


def ell_e_inc(phi, k):
    """E(phi, k) by Carlson: valid for k sin(phi) <= 1, any k."""
    s, c = math.sin(phi), math.cos(phi)
    y = max(1 - (k * s) ** 2, 0.0)
    return s * elliprf(c * c, y, 1) - (k * k) * s ** 3 / 3 * elliprd(c * c, y, 1)


def ell_f_inc(phi, k):
    s, c = math.sin(phi), math.cos(phi)
    return s * elliprf(c * c, max(1 - (k * s) ** 2, 0.0), 1)


def B0(z):
    if z == 0:
        return 0.0
    phi = math.asin(min(1.0, 1.0 / z))
    return 4 * z / math.pi * ell_e_inc(phi, z)


def Q(z):
    phi = math.asin(min(1.0, 1.0 / z))
    return (7 - 8 * z * z - 4 * (1 - z * z) * ell_f_inc(phi, z) / ell_e_inc(phi, z)) / 3


def A_ps(u):
    return (u * u + 2) / (u * math.sqrt(u * u + 4))


def A_uni(u, rho):
    """Lee et al. 2009 Eq. A_fs: (1/(pi rho^2)) int_0^pi [u2 sqrt(u2^2+4) - u1 sqrt(u1^2+4)] dtheta."""
    g = lambda x: x * math.sqrt(x * x + 4)
    if u <= rho:
        f = lambda t: g(u * math.cos(t) + math.sqrt(rho * rho - (u * math.sin(t)) ** 2))
        val = quad(f, 0, math.pi, epsabs=0, epsrel=1e-13, limit=200)[0]
    else:
        tmax = math.asin(rho / u)

        def f(t):
            d = math.sqrt(max(rho * rho - (u * math.sin(t)) ** 2, 0.0))
            return g(u * math.cos(t) + d) - g(u * math.cos(t) - d)
        val = quad(f, 0, tmax, epsabs=0, epsrel=1e-13, limit=200)[0]
    return val / (math.pi * rho * rho)


def B1(z):
    """Yoo Eq. b1def: B0 - (3/2) z rho W_1[z|rho; 1/x], with W_1 = (1/pi) int dtheta int r dr sqrt(1-r^2)/sqrt(r^2+z^2-2rz cos)."""
    def w1():
        # integrate over the unit disc in polar coords about the source centre; 1/x singular at the lens point
        f = lambda r, t: r * math.sqrt(1 - r * r) / math.sqrt(max(r * r + z * z - 2 * r * z * math.cos(t), 1e-300))
        return 2 * dblquad(f, 0, math.pi, 0, 1, epsabs=1e-12, epsrel=1e-10)[0] / math.pi
    # z rho W_1[z|rho; x^-1] with x = rho * dist: the rho cancels
    return B0(z) - 1.5 * z * w1()


if __name__ == "__main__":
    print("B0 limits: B0(1e-6)/(2e-6) =", B0(1e-6) / 2e-6, " B0(1) =", B0(1.0), "(4/pi =", 4 / math.pi, ")  B0(1e4) =", B0(1e4))
    # B0 for z > 1 by the reciprocal-modulus closed form: z E(1/z) - (z - 1/z) K(1/z), parameter 1/z^2
    for z in [1.5, 3.0]:
        m = 1 / z ** 2
        closed = 4 * z / math.pi * (z * ellipe(m) - (z - 1 / z) * ellipk(m))
        print(f"  z={z}: Carlson {B0(z):.15f}  complete-integral form {closed:.15f}")
    print("B0 table:", {z: round(B0(z), 10) for z in [0.1, 0.5, 0.9, 1.0, 1.1, 2.0, 5.0]})
    zs = np.linspace(0.01, 3, 30001)
    zs = zs[np.abs(zs - 1) > 1e-9]
    q = np.array([Q(z) for z in zs])
    i = q.argmin()
    print(f"Q(z): min {q.min():.4f} at z={zs[i]:.3f} (Yoo: -0.38 at z=0.97), Q(0.0001)={Q(1e-4):.6f} (Yoo: max 1 at z=0)")
    zb = np.linspace(0.5, 1.5, 100001)
    bb = np.array([B0(z) for z in zb])
    print(f"max B0 = {bb.max():.10f} at z = {zb[bb.argmax()]:.5f}")
    print("\nExact uniform disc (Lee 2009) vs identities and Gould's A(u) B0(u/rho):")
    for rho in [0.01, 0.05, 0.5, 1.0]:
        a0 = A_uni(0.0, rho)
        print(f"  rho={rho}: A(u=0) = {a0:.15f}  sqrt(1+4/rho^2) = {math.sqrt(1 + 4 / rho ** 2):.15f}")
        for zz in [0.5, 1.0, 2.0]:
            u = zz * rho
            ex = A_uni(u, rho)
            gould = A_ps(u) * B0(zz)
            print(f"     u={zz}rho: exact {ex:.12f}; Gould A*B0 {gould:.12f} ({gould/ex-1:+.2e}); 2nd order "
                  f"{gould*(1+rho*rho/8*Q(zz)) if zz != 1.0 else float('nan'):.12f}")
    # Lee et al. Fig 2 cases: u0=0.1 rho=0.5; u0=0.01 rho=0.05: peak values
    print(f"  Lee Fig.2 peaks: A(u=0.1, rho=0.5) = {A_uni(0.1, 0.5):.6f};  A(u=0.01, rho=0.05) = {A_uni(0.01, 0.05):.6f}")
    print("\nLimb darkening (Yoo): B1(z->0)/z -> 2 - 3 pi/4 =", 2 - 3 * math.pi / 4)
    for z in [0.05, 0.5, 2.0]:
        print(f"  B1({z}) = {B1(z):.8f}  (B1/z = {B1(z)/z:.6f})")
