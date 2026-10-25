"""Curved w0wa-CDM distances with radiation and massive neutrinos (Hogg 1999 + CPL dark energy +
Komatsu et al. 2011 neutrinos), checked against the external-code values in astropy's tests.

Reads the icosmo .ecsv tables downloaded into ../dl_astropy (plain-text parsing only).
"""
import math
import sys
import numpy as np
from scipy.integrate import quad

C_KMS = 299792.458
SIGMA_SB, G, C_SI = 5.670374419e-8, 6.67430e-11, 299792458.0
MPC_KM = 3.0856775814913673e19
KB_EV = 8.617333262e-5
ZETA3 = 1.2020569031595942
NU_PREFAC = 7 / 8 * (4 / 11) ** (4 / 3)  # 0.2271073...


def fermi_ratio(y):
    """Exact F(y)/F(0), F(y) = int x^2 sqrt(x^2+y^2)/(e^x+1) dx (Komatsu et al. 2011 Eqs. 24-25)."""
    if y == 0:
        return 1.0
    f = lambda x: x * x * math.sqrt(x * x + y * y) / (math.exp(x) + 1.0) if x < 700 else 0.0
    # split at the peak region; integrand decays as e^-x
    val = quad(f, 0, 50, epsabs=0, epsrel=1e-13, limit=400)[0] + quad(f, 50, 700, epsabs=0, epsrel=1e-13, limit=400)[0]
    return val / (7 * math.pi ** 4 / 120)


def komatsu_fit(y):
    """Komatsu et al. 2011 Eq. 26 fit used by astropy: (1 + (A y)^p)^(1/p), p = 1.83."""
    A = 180 * ZETA3 / (7 * math.pi ** 4)
    return (1 + (A * y) ** 1.83) ** (1 / 1.83)


class Cosmo:
    def __init__(self, H0, Om0, Ode0=None, w0=-1.0, wa=0.0, Tcmb0=0.0, Neff=3.046, m_nu=(0.0, 0.0, 0.0),
                 exact_nu=True):
        self.H0, self.Om0, self.w0, self.wa = H0, Om0, w0, wa
        h0_si = H0 / MPC_KM
        rho_crit = 3 * h0_si ** 2 / (8 * math.pi * G)
        self.Ogamma0 = 4 * SIGMA_SB / C_SI ** 3 * Tcmb0 ** 4 / rho_crit if Tcmb0 > 0 else 0.0
        self.Neff, self.m_nu, self.exact = Neff, list(m_nu), exact_nu
        self.Tnu = (4 / 11) ** (1 / 3) * Tcmb0
        self.y0 = [m / (KB_EV * self.Tnu) if Tcmb0 > 0 else 0.0 for m in self.m_nu]
        self.Onu0 = self.Ogamma0 * self.nu_rel(0.0)
        if Ode0 is None:  # flat
            Ode0 = 1 - Om0 - self.Ogamma0 - self.Onu0
        self.Ode0 = Ode0
        self.Ok0 = 1 - Om0 - Ode0 - self.Ogamma0 - self.Onu0

    def nu_rel(self, z):
        """rho_nu / rho_gamma, astropy convention: Neff shared equally among the species."""
        if self.Ogamma0 == 0:
            return 0.0
        f = fermi_ratio if self.exact else komatsu_fit
        n = len(self.m_nu)
        return NU_PREFAC * self.Neff / n * sum(f(y / (1 + z)) for y in self.y0)

    def de_scale(self, z):
        return (1 + z) ** (3 * (1 + self.w0 + self.wa)) * math.exp(-3 * self.wa * z / (1 + z))

    def E(self, z):
        zp = 1 + z
        rad = self.Ogamma0 * (1 + self.nu_rel(z)) * zp ** 4 if self.Ogamma0 else 0.0
        return math.sqrt(self.Om0 * zp ** 3 + self.Ok0 * zp ** 2 + self.Ode0 * self.de_scale(z) + rad)

    def DH(self):
        return C_KMS / self.H0

    def chi(self, z1, z2):
        """line-of-sight comoving distance / D_H between z1 and z2 (dimensionless)."""
        if z2 <= 1e4:
            return quad(lambda z: 1 / self.E(z), z1, z2, epsabs=0, epsrel=1e-12, limit=500)[0]
        # integrate in ln(1+z) at high z
        return quad(lambda u: math.exp(u) / self.E(math.exp(u) - 1), math.log1p(z1), math.log1p(z2), epsabs=0,
                    epsrel=1e-12, limit=500)[0]

    def Sk(self, x):
        ok = self.Ok0
        if abs(ok) < 1e-15:
            return x
        s = math.sqrt(abs(ok))
        return math.sinh(s * x) / s if ok > 0 else math.sin(s * x) / s

    def DM(self, z1, z2=None):
        if z2 is None:
            z1, z2 = 0.0, z1
        return self.DH() * self.Sk(self.chi(z1, z2))

    def DA(self, z1, z2=None):
        if z2 is None:
            z1, z2 = 0.0, z1
        return self.DM(z1, z2) / (1 + z2)

    def DL(self, z):
        return self.DM(z) * (1 + z)

    def DA_hogg(self, z1, z2):
        """Hogg Eq. 19, valid only for Ok >= 0."""
        dm1, dm2, dh = self.DM(z1), self.DM(z2), self.DH()
        ok = self.Ok0
        return (dm2 * math.sqrt(1 + ok * dm1 ** 2 / dh ** 2) - dm1 * math.sqrt(1 + ok * dm2 ** 2 / dh ** 2)) / (1 + z2)


def read_ecsv(path):
    meta, rows = {}, []
    for line in open(path):
        if line.startswith("# - {") and ":" in line and "name" not in line:
            k, v = line[5:].strip().rstrip("}").split(":", 1)
            meta[k.strip()] = v.strip()
        elif not line.startswith("#") and line.strip() and not line.startswith("redshift"):
            rows.append([float(x) for x in line.split()])
    return meta, np.array(rows)


def rel(a, b):
    return a / b - 1


if __name__ == "__main__":
    dl = sys.argv[1] if len(sys.argv) > 1 else "../dl_astropy"
    print("== icosmo tables (astropy test_flat_open_closed_icosmo)")
    for name in ["cosmo_flat", "cosmo_open", "cosmo_closed"]:
        meta, t = read_ecsv(f"{dl}/{name}.ecsv")
        c = Cosmo(100 * float(meta["h"]), float(meta["Om"]), float(meta["Ol"]))
        worst = max(max(abs(rel(c.DM(z), dm)), abs(rel(c.DA(z), da)), abs(rel(c.DL(z), dl_))) for z, dm, da, dl_ in t[1:])
        print(f"  {name}: Om={meta['Om']} Ol={meta['Ol']} Ok={c.Ok0:+.2f}: worst rel diff over {len(t)-1} rows = {worst:.2e};"
              f" z=1: DM={c.DM(1.0):.4f} (table {t[t[:,0]==1.0][0,1]})")
    print("== test_flat_z1: FlatLambdaCDM(70, 0.27, Tcmb0=0), z=1 vs Wright/Kempner/iCosmos")
    c = Cosmo(70, 0.27)
    print(f"  DC={c.DM(1):.4f} (3364.5, 3364.8, 3364.7988)  DA={c.DA(1):.4f} (1682.3994)  DL={c.DL(1):.4f} (6729.5976)")
    print("== Kayser, Helbig & Schramm 1997 code values quoted in astropy tests")
    c = Cosmo(100, 0.3)
    print(f"  Flat(100,0.3) DM(1,2) = {c.DM(1,2):.10f} (1313.2232194828466)")
    c = Cosmo(70.4, 0.272)
    print(f"  Flat(70.4,0.272) DA(1,2) = {c.DA(1,2):.10f} (646.22968662822018)")
    print("== astropy non-flat z1z2 values (astropy's own, source not stated)")
    c = Cosmo(100, 0.3, 0.5)
    z1, z2 = (0, 0, 2, 0.5, 1), (2, 1, 1, 2.5, 1.1)
    exp = (3535.931375645655, 2226.430046551708, -1208.6817970036532, 2595.567367601969, 151.36592003406884)
    print("  LambdaCDM(100,0.3,0.5) DM(z1,z2):", [f"{rel(c.DM(a, b), e):+.1e}" for a, b, e in zip(z1, z2, exp)])
    c = Cosmo(100, 1.0, 0.2)
    exp = (-281.31602666724865, 0.0, 248.58093707820436, 843.9331377460543, 1618.6104987686672, 2287.5626543279927)
    print("  LambdaCDM(100,1.0,0.2) DM(0.1,z2):", [f"{(c.DM(0.1, b) - e):+.1e}" for b, e in zip((0, 0.1, 0.2, 0.5, 1.1, 2), exp)])
    c = Cosmo(70.4, 0.2, 0.5)
    print(f"  LambdaCDM(70.4,0.2,0.5) DA(1,2) = {c.DA(1,2):.10f} (620.1175337852428); Hogg Eq.19 gives {c.DA_hogg(1,2):.10f}")
    c = Cosmo(100, 2, 1)
    print(f"  LambdaCDM(100,2,1) Ok={c.Ok0} DA(1,2) = {c.DA(1,2):.10f} (228.42914659246014); Hogg Eq.19 (invalid, Ok<0) gives {c.DA_hogg(1,2):.6f}")
    print("== w0wa vs Mathematica (astropy test_varyde_lumdist_mathematica), rtol 1e-4")
    zs = [0.2, 0.4, 0.9, 1.2]
    for (om, ode, w0, wa, exp) in [(0.2, 0.8, -1.1, 0.2, [1004.0, 2268.62, 6265.76, 9061.84]),
                                   (0.3, 0.7, -0.9, 0.0, [971.667, 2141.67, 5685.96, 8107.41]),
                                   (0.3, 0.7, -0.9, -0.5, [974.087, 2157.08, 5783.92, 8274.08])]:
        c = Cosmo(70, om, ode, w0, wa)
        print(f"  Om={om} w0={w0} wa={wa}:", [f"{c.DL(z):.3f}/{e} ({rel(c.DL(z), e):+.1e})" for z, e in zip(zs, exp)])
    c = Cosmo(70, 0.2, 0.8, -1.1, 0.2)
    print("  de_density_scale(0.5, 1.5):", c.de_scale(0.5), c.de_scale(1.5), "(0.9246310669529021, 0.9184087000251957)")
    c = Cosmo(70, 0.3, 0.6, -0.5)
    print("  wCDM w=-0.5 de scale z=[.1,.2,.5,1.5,2.5]:", [round(c.de_scale(z), 5) for z in [0.1, 0.2, 0.5, 1.5, 2.5]],
          "(1.15369, 1.31453, 1.83712, 3.95285, 6.5479)")
    print("== massive neutrinos: exact Komatsu Eqs. 24-25 vs astropy's Mathematica values (rho_nu/rho_gamma/(prefac Neff))")
    zt = [0.0, 1.0, 2.0, 10.0, 1000.0]
    for m, H0, om, exp in [(100.0, 75, 0.25, [171969, 85984.5, 57323, 15633.5, 171.801]),
                           (0.25, 75, 0.25, [429.924, 214.964, 143.312, 39.1005, 1.11086]),
                           (0.01, 80, 0.30, [17.2347, 8.67345, 5.84348, 1.90671, 1.00021])]:
        for exact in (True, False):
            c = Cosmo(H0, om, Tcmb0=3.0, Neff=3, m_nu=(m, m, m), exact_nu=exact)
            got = [c.nu_rel(z) / (NU_PREFAC * 3) for z in zt]
            print(f"  m={m} eV {'exact' if exact else 'fit  '}:", [f"{g:.6g}({rel(g, e):+.1e})" for g, e in zip(got, exp)])
    c = Cosmo(80, 0.30, Tcmb0=3.0, Neff=3.04, m_nu=(0.0, 0.01, 0.25))
    got = [c.nu_rel(z) / (NU_PREFAC * 3.04) for z in zt]
    print("  mixed [0,0.01,0.25] Neff=3.04 exact:", [f"{g:.6f}" for g in got], "(149.386233, 74.87915, 50.0518, 14.002403, 1.03702333)")
    c = Cosmo(80, 0.30, Tcmb0=3.0, Neff=3, m_nu=(0.01,) * 3, exact_nu=False)
    print("  efunc fit m=0.01:", c.E(1.0), c.E(2.0), "(1.76225893, 2.97022048 astropy, rtol 1e-4)")
    c = Cosmo(80, 0.30, Tcmb0=3.0, Neff=3, m_nu=(0.01,) * 3, exact_nu=True)
    print("  efunc exact m=0.01:", c.E(1.0), c.E(2.0))
    print("== Planck18 preset (astropy: H0=67.66, Om0=0.30966, Tcmb0=2.7255, Neff=3.046, m_nu=[0,0,0.06], Ob0=0.04897)")
    for exact in (True, False):
        p = Cosmo(67.66, 0.30966, Tcmb0=2.7255, Neff=3.046, m_nu=(0, 0, 0.06), exact_nu=exact)
        zstar, rstar = 1089.80, 144.57
        dm = p.DM(zstar)
        print(f"  exact_nu={exact}: Ogamma0={p.Ogamma0:.6e} Onu0={p.Onu0:.6e} Ode0={p.Ode0:.6f} Om0+Onu0={p.Om0+p.Onu0:.5f} (Planck Om 0.3111)")
        print(f"    DM(z*=1089.80)={dm:.3f} Mpc; 100 theta* = 100 r*/DM = {100*rstar/dm:.5f} (Table 2: 1.04119 +- 0.00029)")
        print(f"    DC(1)={p.DM(1):.4f}  DA(1)={p.DA(1):.4f}  DL(1)={p.DL(1):.4f} Mpc")
    p0 = Cosmo(67.66, 0.30966 + 0.0014)  # no radiation at all
    print(f"  without radiation: DM(z*) = {p0.DM(1089.80):.1f} Mpc -> 100 theta* = {100*144.57/p0.DM(1089.80):.5f}")
    # flat LCDM no radiation vs crate doc-test: RXJ1131 D_dt at H0=70, Om=0.3 = 2334.06 Mpc
    c = Cosmo(70, 0.3)
    dd, ds, dds = c.DA(0.295), c.DA(0.654), c.DA(0.295, 0.654)
    print("== crate cross-check: RXJ1131 D_dt =", (1 + 0.295) * dd * ds / dds, "(crate doc test 2334.06)")
