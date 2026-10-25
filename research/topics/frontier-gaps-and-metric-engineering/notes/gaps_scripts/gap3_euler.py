"""Gap 3: 3-D Euler deconvolution (Reid et al. 1990) and Euler inversion (Uieda et al. 2025).

A. Re-implements, from the equations in Uieda et al. (2025, GJI 241, 1536; paper source
   content.tex, Algorithms 1-3; code/euler.py read, not run), the paper's proof-of-concept test and
   compares with its published numbers (paper/variables/synthetic-proof-of-concept.tex):
   true (15000, 12000, -3000) m, base 100 nT; Euler inversion (15045, 12028, -2663) m, 93 nT;
   Euler deconvolution (14626, 11865, -1553) m, 94 nT.
   The data recipe follows the paper's notebook: verde grid [0, 25e3, 0, 20e3] at 300 m and 800 m
   height, dipole 5e11 A m^2 at I = -30, D = 15, 10 nT noise from default_rng(42), base 100 nT;
   x/y derivatives by xarray.differentiate (numpy.gradient, edge_order 1), the up derivative by FFT
   after a linear ramp to zero over half the grid, then Harmonica's own 25 % edge pad.
B. Noise-free analytic derivatives: Euler deconvolution is exact; the structural index is found
   by the minimum weighted misfit of Euler inversion (Algorithm 2).
"""
import numpy as np

MU0_4PI = 1e-7

def grid_coords(w, e, s, n, spacing):  # verde.grid_coordinates, adjust="spacing"
    ne = int(round((e - w) / spacing)) + 1; nn = int(round((n - s) / spacing)) + 1
    return np.meshgrid(np.linspace(w, e, ne), np.linspace(s, n, nn))  # (northing, easting) arrays

def angles_to_vec(intensity, inc, dec):  # harmonica.magnetic_angles_to_vec -> (e, n, u)
    i, d = np.radians(inc), np.radians(dec)
    return np.array([intensity * np.cos(i) * np.sin(d), intensity * np.cos(i) * np.cos(d), -intensity * np.sin(i)])

def dipole_b(e, n, u, src, m):  # nT, (e, n, u) components
    r = np.array([e - src[0], n - src[1], u - src[2]])
    r2 = (r ** 2).sum(0); rm = sum(r[i] * m[i] for i in range(3))
    return [1e9 * MU0_4PI * (3 * rm * r[i] / r2 ** 2.5 - m[i] / r2 ** 1.5) for i in range(3)]

def fft_up_deriv(f, de, dn):
    """Harmonica derivative_upward (kernel -|k|, default 25 % edge pad) on a grid already ramped."""
    pe, pn = int(0.25 * f.shape[1]), int(0.25 * f.shape[0])
    g = np.pad(f, ((pn, pn), (pe, pe)), mode="edge")
    kn = 2 * np.pi * np.fft.fftfreq(g.shape[0], dn)[:, None]; ke = 2 * np.pi * np.fft.fftfreq(g.shape[1], de)[None, :]
    out = np.real(np.fft.ifft2(np.fft.fft2(g) * -np.hypot(kn, ke)))
    return out[pn:pn + f.shape[0], pe:pe + f.shape[1]]

def euler_deconvolution(coords, data, si):
    e, n, u = coords; f, fe, fn, fu = data
    A = np.column_stack([fe, fn, fu, np.full(f.size, float(si))])
    c = e * fe + n * fn + u * fu + si * f
    H = A.T @ A; p = np.linalg.lstsq(A, c, rcond=1e-12)[0]  # min-norm for 2-D sources
    s2 = np.sum((c - A @ p) ** 2) / (f.size - 4)
    return p, s2 * np.linalg.pinv(H)

def euler_inversion(coords, data, si, weights=(1, 0.1, 0.1, 0.025), max_it=20, tol=0.1, nu=0.1):
    e, n, u = coords; N = data[0].size
    dobs = np.concatenate(data); dpred = 0.9 * dobs
    p = euler_deconvolution(coords, data, si)[0]
    w = np.repeat(np.asarray(weights, float), N); winv = [np.full(N, 1 / x) for x in weights]
    def euler_eq(d, p):
        f, fe, fn, fu = np.split(d, 4)
        return (e - p[0]) * fe + (n - p[1]) * fn + (u - p[2]) * fu + si * (f - p[3])
    def merit(d, p):
        return np.linalg.norm((dobs - d) * w) + nu * np.linalg.norm(euler_eq(d, p)), np.linalg.norm((dobs - d) * w)
    m0, _ = merit(dpred, p); history = [p.copy()]; cof = None
    for _ in range(max_it):
        f, fe, fn, fu = np.split(dpred, 4)
        A = -np.column_stack([fe, fn, fu, np.full(N, float(si))])
        B = [np.full(N, float(si)), e - p[0], n - p[1], u - p[2]]  # dEuler/d(f, fe, fn, fu)
        r = dobs - dpred
        Qinv = 1 / sum(b * b * wi for b, wi in zip(B, winv))
        Br = sum(b * ri for b, ri in zip(B, np.split(r, 4)))
        eu = euler_eq(dpred, p)
        cof_new = np.linalg.pinv(A.T @ (Qinv[:, None] * A), rcond=1e-12)
        dp = -cof_new @ (A.T @ (Qinv * (eu + Br)))
        lam = Qinv * (Br + eu + A @ dp)
        dd = r - np.concatenate([wi * b * lam for b, wi in zip(B, winv)])
        m1, _ = merit(dpred + dd, p + dp)
        if m1 > m0:
            break
        p, dpred, cof = p + dp, dpred + dd, cof_new
        history.append(p.copy())
        if abs(m0 - m1) / m0 < tol:
            m0 = m1; break
        m0 = m1
    wrmse = np.sqrt(np.sum(w * (dobs - dpred) ** 2) / dobs.size)
    return p, history, wrmse

# ---------------- A. the paper's proof of concept ----------------
E, Nn = grid_coords(0, 25e3, 0, 20e3, 300.0); U = np.full_like(E, 800.0)
inc, dec = -30.0, 15.0; src = (15e3, 12e3, -3e3)
b = dipole_b(E, Nn, U, src, angles_to_vec(5e11, inc, dec))
fhat = angles_to_vec(1, inc, dec)
tfa = sum(b[i] * fhat[i] for i in range(3))
tfa = tfa + np.random.default_rng(42).normal(0, 10, size=tfa.shape) + 100.0
de, dn = E[0, 1] - E[0, 0], Nn[1, 0] - Nn[0, 0]
f_e = np.gradient(tfa, de, axis=1); f_n = np.gradient(tfa, dn, axis=0)
pn, pe = tfa.shape[0] // 2, tfa.shape[1] // 2
ramp = np.pad(tfa, ((pn, pn), (pe, pe)), mode="linear_ramp", end_values=0)
f_u = fft_up_deriv(ramp, de, dn)[pn:pn + tfa.shape[0], pe:pe + tfa.shape[1]]
coords = (E.ravel(), Nn.ravel(), U.ravel()); data = (tfa.ravel(), f_e.ravel(), f_n.ravel(), f_u.ravel())
p_ed, _ = euler_deconvolution(coords, data, 3)
p_ei, hist, _ = euler_inversion(coords, data, 3)
print(f"grid {tfa.shape} spacing ({de:.1f}, {dn:.1f}) m")
print(f"Euler deconvolution: ({p_ed[0]:.0f}, {p_ed[1]:.0f}, {p_ed[2]:.0f}) m, base {p_ed[3]:.0f} nT   [paper: 14626, 11865, -1553, 94]")
print(f"Euler inversion:     ({p_ei[0]:.0f}, {p_ei[1]:.0f}, {p_ei[2]:.0f}) m, base {p_ei[3]:.0f} nT   [paper: 15045, 12028, -2663, 93; {len(hist)-1} iterations, paper 6]")

# ---------------- B. analytic derivatives, structural index by misfit ----------------
def numeric_grad(fun, e, n, u, h=1.0):
    return [(fun(e + h, n, u) - fun(e - h, n, u)) / (2 * h), (fun(e, n + h, u) - fun(e, n - h, u)) / (2 * h), (fun(e, n, u + h) - fun(e, n, u - h)) / (2 * h)]
srcB = (10e3, 8e3, -2e3)
sources = {
    "dipole (SI 3)": lambda e, n, u: sum(c * fc for c, fc in zip(dipole_b(e, n, u, srcB, angles_to_vec(1e11, 60, 0)), angles_to_vec(1, 60, 0))),
    "point mass gz (SI 2)": lambda e, n, u: 1e13 * (u - srcB[2]) / ((e - srcB[0]) ** 2 + (n - srcB[1]) ** 2 + (u - srcB[2]) ** 2) ** 1.5,
    # an infinite line through srcB at azimuth 30 deg: Euler is rank-deficient along strike
    "line mass, azimuth 30 deg, gz (SI 1)": lambda e, n, u: 1e9 * (u - srcB[2]) / (((e - srcB[0]) * np.cos(np.radians(30)) - (n - srcB[1]) * np.sin(np.radians(30))) ** 2 + (u - srcB[2]) ** 2),
}
EB, NB = grid_coords(0, 20e3, 0, 16e3, 250.0); UB = np.full_like(EB, 500.0)
cb = (EB.ravel(), NB.ravel(), UB.ravel())
for name, fun in sources.items():
    f0 = fun(EB, NB, UB) + 50.0
    g = numeric_grad(fun, EB, NB, UB)
    d = (f0.ravel(), g[0].ravel(), g[1].ravel(), g[2].ravel())
    true_si = int(name.split("SI ")[1][0])
    p, _ = euler_deconvolution(cb, d, true_si)
    fits = {si: euler_inversion(cb, d, si) for si in (1, 2, 3)}
    best = min(fits, key=lambda s: fits[s][2])
    perp = (p[0] - srcB[0]) * np.cos(np.radians(30)) - (p[1] - srcB[1]) * np.sin(np.radians(30)) if "line" in name else np.hypot(p[0] - srcB[0], p[1] - srcB[1])
    print(f"{name}: ED at true SI: horizontal error {abs(perp):.2e} m, depth error {p[2]-srcB[2]:.2e} m, base {p[3]-50:.1e}; "
          f"EI WRMSE by SI " + ", ".join(f"{s}: {fits[s][2]:.3g}" for s in (1, 2, 3)) + f" -> SI {best}")
