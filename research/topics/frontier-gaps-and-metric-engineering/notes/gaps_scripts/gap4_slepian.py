"""Gap 4: Slepian functions of an arbitrary region, and basin averages.

Concentration matrix D_(lm)(l'm') = (1/4pi) Int_R Pbar_lm Pbar_l'm' T_m T_m' dOmega (4pi-normalised
real harmonics, T = cos m phi or sin m phi, as geopotential::slepian), the definition of Simons,
Dahlen & Wieczorek (2006) eq. (Dlmlmpdef) in the crate's normalisation. Integrated as Simons'
kernelc.m does it: Gauss-Legendre in colatitude, exact longitude integrals over the region's
longitude intervals at each node; here improved by splitting the colatitude range at every vertex
or extremum (panels) and clustering nodes at panel ends (theta = a + (b-a)(1-cos t)/2).
Checks: polar cap vs the crate's pinned values; a rotated cap must have the same eigenvalues
(rotation invariance); trace = (L+1)^2 A/4pi (SDW eq. shannon); box and polygon areas; the basin
average of Pbar_20 over a cap; Swenson & Wahr (2002) Gaussian W_l, eq. 34 against eq. 33.
"""
import numpy as np

def legendre_table(L, theta):
    """Pbar_lm(cos theta), 4pi-normalised, shape (L+1, L+1, len(theta))."""
    x, s = np.cos(theta), np.sin(theta)
    P = np.zeros((L + 1, L + 1, len(theta)))
    P[0, 0] = 1.0
    for m in range(1, L + 1):
        P[m, m] = (np.sqrt(3.0) if m == 1 else np.sqrt((2 * m + 1) / (2 * m))) * s * P[m - 1, m - 1]
    for m in range(0, L):
        P[m + 1, m] = np.sqrt(2 * m + 3) * x * P[m, m]
        for l in range(m + 2, L + 1):
            a = np.sqrt((2 * l - 1) * (2 * l + 1) / ((l - m) * (l + m)))
            b = np.sqrt((2 * l + 1) * (l + m - 1) * (l - m - 1) / ((l - m) * (l + m) * (2 * l - 3)))
            P[l, m] = a * x * P[l - 1, m] - b * P[l - 2, m]
    return P

def basis(L):
    return [(l, m, t) for l in range(L + 1) for m in range(l + 1) for t in ((0,) if m == 0 else (0, 1))]

def lon_integrals(L, intervals):
    """I[t, m, t', m'] = sum over intervals of Int T_tm T_t'm' dphi, t = 0 cos, 1 sin."""
    m = np.arange(L + 1)
    I = np.zeros((2, L + 1, 2, L + 1))
    def S(f, p):  # antiderivative of cos(f p)
        return np.where(f == 0, p, np.sin(f * p) / np.where(f == 0, 1, f))
    def C(f, p):  # antiderivative of sin(f p)
        return np.where(f == 0, 0.0, -np.cos(f * p) / np.where(f == 0, 1, f))
    a, b = m[:, None], m[None, :]
    for lo, hi in intervals:
        d = lambda F, f: F(f, hi) - F(f, lo)
        I[0, :, 0, :] += 0.5 * (d(S, a - b) + d(S, a + b))
        I[1, :, 1, :] += 0.5 * (d(S, a - b) - d(S, a + b))
        I[1, :, 0, :] += 0.5 * (d(C, a + b) + d(C, a - b))  # sin(a) cos(b)
    I[0, :, 1, :] = I[1, :, 0, :].T
    return I

def concentration(L, intervals_at, breaks, n_gl=64):
    B = basis(L); idx = {k: i for i, k in enumerate(B)}
    D = np.zeros((len(B), len(B)))
    t, w = np.polynomial.legendre.leggauss(n_gl)
    t = (t + 1) * np.pi / 2; w = w * np.pi / 2
    for a, b in zip(breaks[:-1], breaks[1:]):
        th = a + (b - a) * (1 - np.cos(t)) / 2
        wt = w * (b - a) * np.sin(t) / 2 * np.sin(th)
        P = legendre_table(L, th)
        for k in range(len(th)):
            I = lon_integrals(L, intervals_at(th[k]))
            if not I.any():
                continue
            v = np.array([P[l, m, k] for (l, m, tt) in B])
            Ik = np.array([[I[t1, m1, t2, m2] for (_, m2, t2) in B] for (_, m1, t1) in B])
            D += wt[k] * np.outer(v, v) * Ik
    return D / (4 * np.pi), B

def cap_intervals(theta_c, phi_c, radius):
    def f(th):
        q = (np.cos(radius) - np.cos(th) * np.cos(theta_c)) / max(np.sin(th) * np.sin(theta_c), 1e-300)
        if q <= -1: return [(0.0, 2 * np.pi)]
        if q >= 1: return []
        dl = np.arccos(q); return [(phi_c - dl, phi_c + dl)]
    return f

def polygon_intervals(lonlat_deg):
    """kernelc/dphregion convention: edges straight in (lon, lat); crossings paired west to east."""
    V = np.radians(np.asarray(lonlat_deg, float)); lon, colat = V[:, 0], np.pi / 2 - V[:, 1]
    def f(th):
        xs = []
        for i in range(len(V)):
            j = (i + 1) % len(V)
            c0, c1 = colat[i], colat[j]
            if (c0 - th) * (c1 - th) < 0:
                xs.append(lon[i] + (lon[j] - lon[i]) * (th - c0) / (c1 - c0))
        xs.sort()
        return list(zip(xs[0::2], xs[1::2]))
    return f, sorted(set(colat))

L, TH = 18, np.radians(40.0)
# 1. polar cap: the crate's pinned values (Shannon 42.229, lambda_42 = 0.579, 42 above 1/2)
D, B = concentration(L, lambda th: [(0.0, 2 * np.pi)] if th < TH else [], [0.0, TH], 40)
lam_polar = np.sort(np.linalg.eigvalsh(D))[::-1]
print(f"polar cap: trace {np.trace(D):.6f} (Shannon {(L+1)**2*(1-np.cos(TH))/2:.6f}); lambda_42 {lam_polar[41]:.4f}; >1/2: {(lam_polar > 0.5).sum()}; lambda_1..3 {lam_polar[:3]}")
# 2. the same cap centred at colatitude 60 deg, longitude 30 deg
tc, pc = np.radians(60.0), np.radians(30.0)
for n_gl in (24, 48, 96):
    Dr, _ = concentration(L, cap_intervals(tc, pc, TH), [tc - TH, tc + TH], n_gl)
    lam_rot = np.sort(np.linalg.eigvalsh(Dr))[::-1]
    print(f"rotated cap, {n_gl} nodes per panel: max |lambda - polar| {np.abs(lam_rot - lam_polar).max():.2e}; trace error {np.trace(Dr) - (L+1)**2*(1-np.cos(TH))/2:.2e}")
# 3. a lat-lon box and a triangle as polygons: D_00 = A/4pi, trace = (L+1)^2 A/4pi
box = [(10, -20), (50, -20), (50, 15), (10, 15)]
f, br = polygon_intervals(box)
Db, _ = concentration(12, f, br, 32)
A = (np.sin(np.radians(15)) - np.sin(np.radians(-20))) * np.radians(40)
print(f"box: D00 {Db[0,0]:.15f} vs A/4pi {A/4/np.pi:.15f}; trace/(13^2 A/4pi) - 1 = {np.trace(Db)/(169*A/4/np.pi)-1:.1e}")
tri = [(0, 0), (40, 10), (15, 50)]
f, br = polygon_intervals(tri)
Dt, _ = concentration(12, f, br, 32)
# reference area of the plate-carree triangle: Int cos(lat) dlat dlon, by fine midpoint sum
lat = np.radians(np.linspace(0, 50, 200001)); dl = lat[1] - lat[0]
width = np.array([sum(b - a for a, b in f(np.pi / 2 - la)) for la in lat[:-1] + dl / 2])
At = np.sum(width * np.cos(lat[:-1] + dl / 2)) * dl
print(f"triangle: D00*4pi {Dt[0,0]*4*np.pi:.10f} vs area {At:.10f}; trace/(13^2 A/4pi) - 1 = {np.trace(Dt)/(169*At/4/np.pi)-1:.1e}")
# 4. basin average of Pbar_20 over the polar cap: (4pi/A) theta_lm . f_lm vs closed form
c = np.cos(TH)
row = D[0]  # theta_lm = D[(0,0),(l,m,t)] since Pbar_00 = 1
i20 = B.index((2, 0, 0))
avg = row[i20] / row[0]
exact = 2 * np.pi * np.sqrt(5) * (-(c ** 3 - c) / 2) / (2 * np.pi * (1 - c))
print(f"basin average of Pbar_20 over the 40 deg cap: {avg:.15f} vs closed form {exact:.15f}")
# 5. Swenson & Wahr (2002) eqs 30-34: W_l by recursion vs the integral, r_1/2 = 200 km
a = 6371.0
for r_half in (200.0, 500.0):
    b = np.log(2) / (1 - np.cos(r_half / a))
    W = [1 / (2 * np.pi), (1 / (2 * np.pi)) * ((1 + np.exp(-2 * b)) / (1 - np.exp(-2 * b)) - 1 / b)]
    for l in range(1, 120):
        W.append(-(2 * l + 1) / b * W[l] + W[l - 1])
    g, wg = np.polynomial.legendre.leggauss(4000)  # x = cos gamma on [-1, 1]
    kern = b / (2 * np.pi) * np.exp(-b * (1 - g)) / (1 - np.exp(-2 * b))
    from numpy.polynomial import legendre as Lg
    Wd = [np.sum(wg * kern * Lg.legval(g, [0] * l + [1])) for l in range(121)]  # eq 33 with Pbar_l0 = sqrt(2l+1) P_l
    worst = max(abs(W[l] - Wd[l]) / abs(Wd[0]) for l in range(121))
    first_bad = next((l for l in range(121) if abs(W[l] - Wd[l]) > 1e-3 * abs(Wd[0])), None)
    print(f"Gaussian r1/2={r_half:.0f} km: 2pi W_0 = {2*np.pi*W[0]:.3f}; recursion vs integral max |dW|/W_0 to l=120: {worst:.1e}; first l off by >1e-3 W_0: {first_bad}")
# 6. Table of Simons, Dahlen & Wieczorek (2006): rounded (L+1)^2 A/4pi against the printed N
table = {"Greenland": (0.44, [0, 1, 2, 3]), "Australia": (1.50, [1, 3, 5, 9]), "South America": (3.50, [2, 6, 13, 22]),
         "North America": (3.98, [2, 7, 14, 25]), "Africa": (5.78, [3, 10, 21, 36]), "Eurasia": (9.98, [5, 17, 36, 62])}
for name, (pct, Ns) in table.items():
    ours = [(Lx + 1) ** 2 * pct / 100 for Lx in (6, 12, 18, 24)]
    print(f"  {name:14s} {pct:5.2f}%: (L+1)^2 A/4pi = " + ", ".join(f"{o:.2f}" for o in ours) + f"  printed {Ns}")
print(f"Harig & Simons (2012): N = 20 at L = 60 needs A/4pi in [{19.5/3721*100:.3f}, {20.5/3721*100:.3f}] %")
