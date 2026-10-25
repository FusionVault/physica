"""Kerr shadow outline at any inclination (Bardeen 1973, as given by Cunha & Herdeiro 2018
Eqs. eq_eta, eq_chi, eq_shadow) and the size measures used in the literature.

Units: M = G = c = 1. Bardeen's notation: xi = L_z/E, eta = Q/E^2 (Cunha & Herdeiro call these
eta and chi). Image coordinates alpha = -xi/sin i, beta = +-sqrt(eta + a^2 cos^2 i - xi^2 cot^2 i).
"""
import numpy as np


def xi_eta(r, a):
    xi = (r * r * (3 - r) - a * a * (r + 1)) / (a * (r - 1))
    eta = r ** 3 * (4 * a * a - r * (r - 3) ** 2) / (a * a * (r - 1) ** 2)
    return xi, eta


def beta2(r, a, inc):
    xi, eta = xi_eta(r, a)
    return eta + a * a * np.cos(inc) ** 2 - xi * xi / np.tan(inc) ** 2


def bisect(f, lo, hi, n=200):
    flo = f(lo)
    for _ in range(n):
        mid = 0.5 * (lo + hi)
        if (f(mid) > 0) == (flo > 0):
            lo, flo = mid, f(mid)
        else:
            hi = mid
    return 0.5 * (lo + hi)


def photon_radii(a):
    # equatorial prograde / retrograde circular photon orbits (Bardeen, Press & Teukolsky 1972)
    rp = 2 * (1 + np.cos(2 / 3 * np.arccos(-abs(a))))
    rm = 2 * (1 + np.cos(2 / 3 * np.arccos(abs(a))))
    return rp, rm


def pole_radius(a):
    # largest root of r^3 - 3r^2 + a^2 r + a^2 = 0 (xi = 0), Newton from r = 3 (as in the crate)
    r = 3.0
    for _ in range(100):
        r -= (r ** 3 - 3 * r * r + a * a * r + a * a) / (3 * r * r - 6 * r + a * a)
    return r


def outline(a, inc, n=400001):
    """Upper half of the shadow edge, (alpha, beta), clustered toward both ends."""
    if a == 0:
        t = np.linspace(0, np.pi, n)
        return 27 ** 0.5 * np.cos(t), 27 ** 0.5 * np.sin(t)
    rp, rm = photon_radii(a)
    r0 = pole_radius(a)
    f = lambda r: beta2(r, a, inc)
    if abs(inc - np.pi / 2) < 1e-15:
        r1, r2 = rp, rm
    else:
        r1, r2 = bisect(f, rp + 1e-15, r0), bisect(f, r0, rm)
    t = np.linspace(0, np.pi, n)
    r = r1 + (r2 - r1) * (1 - np.cos(t)) / 2
    xi, eta = xi_eta(r, a)
    al = -xi / np.sin(inc)
    b2 = np.maximum(beta2(r, a, inc), 0.0)
    return al, np.sqrt(b2)


def measures(a, inc, n=400001):
    al, be = outline(a, inc, n)
    xmin, xmax = al.min(), al.max()
    D = (xmax + xmin) / 2  # Takahashi 2004 x_c; JP2010 use |.|
    width = xmax - xmin
    height = 2 * be.max()
    # JP2010 / Johannsen 2013 / Chan 2013: R from (D, 0), angle tan(alpha) = y'/x' about the ORIGIN
    th0 = np.unwrap(np.arctan2(be, al))
    R = np.hypot(al - D, be)
    trap = lambda y, x: float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(x)))
    span0 = abs(th0[-1] - th0[0])  # should be pi
    meanR_origin = abs(trap(R, th0)) / span0
    A_origin = 2 * np.sqrt(abs(trap((R - meanR_origin) ** 2, th0)) / span0)
    # the same average with the angle about the centre (D, 0)
    thc = np.unwrap(np.arctan2(be, al - D))
    meanR_centre = abs(trap(R, thc)) / abs(thc[-1] - thc[0])
    # median radius about the centre, over uniformly spaced angle
    grid = np.linspace(min(thc[0], thc[-1]), max(thc[0], thc[-1]), 200001)
    order = np.argsort(thc)
    Rg = np.interp(grid, thc[order], R[order])
    medR = float(np.median(Rg))
    # area-equivalent radius (Kocherlakota et al. 2021 use the areal radius)
    area = 2 * abs(trap(be, al))
    return dict(xmin=xmin, xmax=xmax, D=D, width=width, height=height, meanR_origin=meanR_origin,
                A=A_origin, meanR_centre=meanR_centre, medianR=medR, areal=np.sqrt(area / np.pi), span0=span0)


def johannsen2013_fit(i, a):
    """Johannsen 2013 (arXiv:1501.02814) Eqs. L_Kerrfit, D_Kerrfit; i in radians, a in M."""
    L1 = 10.3907 - 0.520134 * a ** 2.19644 - 0.109423 * a ** 11.5509
    L2 = (0.141644 * a ** 2.06874 + 0.077576 * a ** 2.80628 - 0.305764 * a ** 4.44562 + 0.603316 * a ** 6.58314
          - 0.801615 * a ** 9.97291 + 0.892182 * a ** 15.1866 - 0.508699 * a ** 18.8103)
    n1 = -2.00026 - 0.0424525 * a ** 6.43112 + 0.0503142 * a ** 36.5515 + 0.118102 * a ** 153.450
    ph = 3.14246 + 0.0970934 * a ** 5.40872 + 0.0964881 * a ** 13.3813 - 0.0990121 * a ** 196.868
    D1 = 2 * a + 0.188324 * a ** 3.29029 + 0.153158 * a ** 10.4479 + 0.0982701 * a ** 48.1825 + 0.0500847 * a ** 328.089
    n2 = (1 - 0.0149883 * a ** 0.656102 + 0.0188647 * a ** 0.729558 + 0.0370128 * a ** 2.49086 + 0.147758 * a ** 30.2625
          - 0.163190 * a ** 31.1614 - 0.0148721 * a ** 259.142)
    return L1 + L2 * np.cos(n1 * i + ph), D1 * np.sin(n2 * i)


def chan2013_fit(i_deg, a):
    """Chan, Psaltis & Ozel 2013 (arXiv:1303.5057) Eq. Rapprox; arguments in degrees."""
    R0 = 5.2 - 0.209 * a + 0.445 * a * a - 0.567 * a ** 3
    R1 = (0.24 - 3.3 / ((a - 0.9017) ** 2 + 0.059)) * 1e-3
    return R0 + R1 * np.cos(np.radians(2.14 * i_deg - 22.9))


if __name__ == "__main__":
    s27 = 27 ** 0.5
    # 1. pole-on limit reproduces kerr_pole_on_shadow: sqrt(8 r^3/(r^2-1)) at xi = 0
    for a in [0.3, 0.9, 0.999]:
        r0 = pole_radius(a)
        _, eta0 = xi_eta(r0, a)
        assert abs(np.sqrt(eta0 + a * a) / np.sqrt(8 * r0 ** 3 / (r0 * r0 - 1)) - 1) < 1e-13
        m = measures(a, 1e-4, 20001)
        print(f"a={a} i->0: meanR {m['meanR_origin']:.10f} vs sqrt(eta+a^2) {np.sqrt(eta0+a*a):.10f}")
    # 2. Cunha & Herdeiro: y(x = 2a) = 3 sqrt 3 at theta_o = pi/2 (r = 3 exactly)
    for a in [0.5, 0.9, 0.99]:
        xi, eta = xi_eta(3.0, a)
        print(f"a={a}: r=3 -> xi={xi:.15f} (=-2a), eta={eta:.15f} (=27)")
    # 3. edge-on extremes from the equatorial photon orbits
    for a in [0.5, 0.9, 0.998, 1 - 1e-9]:
        rp, rm = photon_radii(a)
        xp, _ = xi_eta(rp, a)
        xm, _ = xi_eta(rm, a)
        print(f"a={a}: r_pro={rp:.12f} r_retro={rm:.12f} alpha_min={-xp:.12f} alpha_max={-xm:.12f} "
              f"width={xp - xm:.12f} D={-(xp + xm) / 2:.12f}")
    # 4. table of measures
    print("\n a     i    width    height   D        <R>_JP   2<R>_JP   L_fit(J13)  D_fit(J13) <R>Chan  medianR  arealR  A")
    for a in [0.0, 0.5, 0.9, 0.998]:
        for ideg in [17, 30, 60, 90]:
            inc = np.radians(ideg)
            m = measures(a, inc)
            L, Df = johannsen2013_fit(inc, a)
            print(f"{a:5.3f} {ideg:3d} {m['width']:.6f} {m['height']:.6f} {m['D']:.6f} {m['meanR_origin']:.6f} "
                  f"{2*m['meanR_origin']:.6f} {L:.6f} {Df:.6f} {chan2013_fit(ideg, a):.6f} {m['medianR']:.6f} "
                  f"{m['areal']:.6f} {m['A']:.6f}")
    # convergence check
    m1, m2 = measures(0.9, np.pi / 2, 100001), measures(0.9, np.pi / 2, 400001)
    print("\nconvergence a=0.9 i=90:", m1['meanR_origin'] - m2['meanR_origin'], m1['A'] - m2['A'])
    m = measures(0.9, np.pi / 2, 1600001)
    print("a=0.9 i=90 (1.6M pts):", {k: round(v, 10) for k, v in m.items()})
    # 5. how small can the Kerr shadow be? delta = L/(6 sqrt 3) - 1 for each measure
    print("\nfractional deviation from 6 sqrt 3 at a -> 1, i = 90 deg:")
    m = measures(1 - 1e-6, np.pi / 2, 1600001)
    for k in ['meanR_origin', 'meanR_centre', 'medianR', 'areal']:
        print(f"  {k}: R = {m[k]:.6f}, delta = {m[k] / s27 - 1:+.5f}")
    print(f"  width/2: {m['width']/2:.6f} delta={m['width']/2/s27-1:+.5f}; height/2: {m['height']/2:.6f} delta={m['height']/2/s27-1:+.5f}")
