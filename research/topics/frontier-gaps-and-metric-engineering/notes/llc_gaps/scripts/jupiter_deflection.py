"""Born-integral light deflection by a zonal gravity field, checked against
Crosta & Mignard 2006 (astro-ph/0512359) Eqs. defqu1-defqu2 and Appendix B, and against the exact
equatorial-ray zonal series derived here. Jupiter values from the Durante2020 PDS SHADR file
(unnormalised C_n0 = -J_n), GM and R from its header.

Reads ../dl_pds_shadr/GRV_JUGR_SHA_JUP_U_DUR2020_20200203_V01.TAB as plain text.
"""
import math
import sys
import numpy as np
from numpy.polynomial import legendre as L

C = 299792458.0
MAS = math.pi / 180 / 3600 * 1e-3
UAS = MAS * 1e-3


def read_shadr(path):
    lines = open(path).read().splitlines()
    head = [x.strip() for x in lines[0].split(",")]
    R, GM, sGM, deg, order, norm = float(head[0]) * 1e3, float(head[1]) * 1e9, float(head[2]) * 1e9, int(head[3]), int(head[4]), int(head[5])
    J = {}
    for ln in lines[1:]:
        p = [x.strip() for x in ln.split(",")]
        if len(p) >= 6 and p[1] == "0":
            J[int(p[0])] = (-float(p[2]), float(p[4]))  # J_n = -C_n0 (unnormalised)
    return dict(R=R, GM=GM, sGM=sGM, deg=deg, order=order, norm=norm, J=J)


def grad_U(p, GM, R, J):
    """Gradient of U = (GM/r)[1 - sum J_n (R/r)^n P_n(z/r)], body frame, z = spin axis."""
    x, y, z = p
    r = math.sqrt(x * x + y * y + z * z)
    s = z / r
    g = -GM / r ** 3 * np.array(p)
    for n, jn in J.items():
        Pn = L.legval(s, [0] * n + [1])
        dPn = L.legval(s, L.legder([0] * n + [1]))
        # V_n = -GM J_n R^n P_n(s) / r^(n+1); grad s = (e_z - s r_hat)/r
        f = -GM * jn * R ** n
        rhat = np.array(p) / r
        dV_dr = f * Pn * (-(n + 1)) / r ** (n + 2)
        grad_s = (np.array([0, 0, 1.0]) - s * rhat) / r
        g = g + dV_dr * rhat + f * dPn / r ** (n + 1) * grad_s
    return g


def born(point, t, GM, R, J, gamma=1.0, nodes=400):
    """(1+gamma)/c^2 int grad_perp U dl along the infinite line point + l t."""
    t = np.array(t, float) / np.linalg.norm(t)
    p0 = np.array(point, float)
    p0 = p0 - np.dot(p0, t) * t
    b = np.linalg.norm(p0)
    x, w = L.leggauss(nodes)
    th = math.pi / 2 * x
    s = np.zeros(3)
    for th_i, w_i in zip(th, w):
        p = p0 + b * math.tan(th_i) * t
        g = grad_U(p, GM, R, J)
        s += math.pi / 2 * w_i * b / math.cos(th_i) ** 2 * (g - np.dot(g, t) * t)
    return (1 + gamma) * s / C ** 2, p0, b


def crosta(p0, t, b, GM, R, J2, gamma=1.0):
    """Crosta & Mignard Eqs. defqu1-defqu2 (chi -> 0); n points from the ray toward the planet."""
    z = np.array([0, 0, 1.0])
    n = -p0 / b
    m = np.cross(t, n)
    M = GM / C ** 2
    d1 = 2 * (1 + gamma) * M / b * (1 + J2 * R * R / b / b * (1 - 2 * np.dot(n, z) ** 2 - np.dot(t, z) ** 2))
    d2 = 4 * (1 + gamma) * M * J2 * R * R / b ** 3 * np.dot(m, z) * np.dot(n, z)
    return d1 * n + d2 * m


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "../dl_pds_shadr/GRV_JUGR_SHA_JUP_U_DUR2020_20200203_V01.TAB"
    f = read_shadr(path)
    R, GM, J = f["R"], f["GM"], f["J"]
    print(f"DUR2020 header: R={R} m GM={GM:.4e} m^3/s^2 (+-{f['sGM']:.1e}, 3 sigma) degree={f['deg']} order={f['order']} normalisation flag={f['norm']} (0 = unnormalised)")
    for n in sorted(J):
        print(f"  J{n} = {J[n][0]*1e6:+.7f}e-6 +- {J[n][1]*1e6:.4f}e-6 (3 sigma)")
    J2 = J[2][0]
    print(f"  J2 without the 6.72e-8 permanent tide: {(J2 - 6.72e-8)*1e6:.4f}e-6; normalised C20bar = -J2/sqrt5 = {-J2/math.sqrt(5):.10e}")
    mono = 4 * GM / (C * C * R)
    print(f"\nmonopole at the limb 4GM/(c^2 R) = {mono:.6e} rad = {mono/UAS:.1f} uas  (Crosta & Mignard Table 1: 16280)")
    print(f"J2 at the limb, z perpendicular to the ray: {mono*J2/UAS:.2f} uas (Table 1: 239)")
    # Crosta's own parameter set (Table 2): M = 1.41 m, R = 7.13e7 m, J2 = 0.014736
    m_cm, R_cm, J2_cm = 1.41, 7.13e7, 0.014736
    print(f"Crosta Table 2 parameters: monopole {4*m_cm/R_cm/UAS:.0f} uas, quadrupole {4*m_cm/R_cm*J2_cm/UAS:.1f} uas")

    zonal2 = {2: J2}
    print("\n== J2-only field: Born integral vs Crosta & Mignard, b = R, various orientations")
    cases = {
        "equatorial ray (t.z=0, n.z=0)": ([R, 0, 0], [0, 1, 0]),
        "ray over the pole (n || z)": ([0, 0, R], [0, 1, 0]),
        "ray along the spin axis (t || z)": ([R, 0, 0], [0, 0, 1]),
        "tilted 30 deg, closest point at 45 deg lat": ([R * math.cos(math.radians(45)), 0, R * math.sin(math.radians(45))],
                                                     [0, math.cos(math.radians(30)), math.sin(math.radians(30))]),
    }
    for name, (pt, t) in cases.items():
        a, p0, b = born(pt, t, GM, R, zonal2)
        tt = np.array(t, float) / np.linalg.norm(t)
        cr = crosta(p0, tt, b, GM, R, J2)
        mono_v = 4 * GM / (C * C * b) * (-p0 / b)
        q_b, q_c = a - mono_v, cr - mono_v
        sin2 = 1 - tt[2] ** 2
        print(f"  {name}: |alpha|={np.linalg.norm(a)/UAS:.4f} uas; quad Born {np.linalg.norm(q_b)/UAS:.5f} uas, Crosta {np.linalg.norm(q_c)/UAS:.5f};"
              f" |dQ| = 4GM J2 R^2 sin^2/(c^2 b^3) = {mono*J2*sin2*(R/b)**3/UAS:.5f} (b/R={b/R:.4f}, sin^2={sin2:.3f}); |Born-Crosta|/|quad| = {np.linalg.norm(q_b-q_c)/max(np.linalg.norm(q_c),1e-30):.1e}")
    print("\n== full DUR2020 zonal field, equatorial ray at b = R and b = 2R: Born vs exact series")
    full = {n: v[0] for n, v in J.items()}
    for bfac in [1.0, 2.0]:
        b = bfac * R
        a, p0, _ = born([b, 0, 0], [0, 1, 0], GM, R, full, nodes=800)
        even = 1 + sum((-1) ** (k + 1) * full.get(2 * k, 0) * (R / b) ** (2 * k) for k in range(1, 7))
        odd = sum((-1) ** (k + 1) * full.get(2 * k + 1, 0) * (R / b) ** (2 * k + 1) for k in range(1, 6))
        m0 = 4 * GM / (C * C * b)
        print(f"  b={bfac}R: Born in-plane {-a[0]/UAS:.6f} uas, series {m0*even/UAS:.6f}; Born z {a[2]/UAS:+.6f} uas, series {m0*odd/UAS:+.6f}")
        print(f"     terms (uas): J2 {m0*full[2]*(R/b)**2/UAS:.4f}, -J4 {-m0*full[4]*(R/b)**4/UAS:.4f}, +J6 {m0*full[6]*(R/b)**6/UAS:.4f}, -J8 {-m0*full[8]*(R/b)**8/UAS:.4f}")
    # Born integral of a J2 field: GR second-order term for scale
    print(f"\n2nd-order monopole (15 pi/4)(GM/c^2R)^2 at the limb = {15*math.pi/4*(GM/C**2/R)**2/UAS:.3f} uas")
