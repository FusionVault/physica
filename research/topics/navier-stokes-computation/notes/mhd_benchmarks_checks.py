"""Checks for mhd_benchmarks.md: reference values for the incompressible MHD solver of issue #55.

Run: python3 -P mhd_benchmarks_checks.py   (numpy + scipy; about 10 minutes on one core)
Every (computed) number quoted in mhd_benchmarks.md is asserted here; exits 1 on any failure.
numpy is needed because the values come from an independent spectral MHD solver (section 1)
and a Fourier-Galerkin eigenvalue solver for the kinematic dynamo (section 2).
"Paper" values are printed in, or read off a figure of, the cited source; "computed" are made here.
"""
import math
import sys

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
from scipy.optimize import brentq

FAILS = []


def check(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""), flush=True)
    if not ok:
        FAILS.append(name)


def close(a, b, tol):
    return abs(a - b) <= tol


# ---------------------------------------------------------------------------------------------
# 1. 2-D incompressible MHD, vorticity / flux form, 2*pi box, 2/3 rule, integrating-factor RK4.
#    u = (d_y psi, -d_x psi), b = (d_y a, -d_x a), omega = -lap psi, j = d_x b_y - d_y b_x = -lap a.
#    d_t omega + u.grad omega = b.grad j + nu lap omega;  d_t a + u.grad a = eta lap a.
#    E = <u^2 + b^2>/2, D = nu <omega^2> + eta <j^2> = -dE/dt, <.> = box mean.
# ---------------------------------------------------------------------------------------------
def ot_run(N, nu, eta, psi_f, a_f, T, dt, every=1):
    x = 2 * np.pi * np.arange(N) / N
    X, Y = np.meshgrid(x, x, indexing="ij")
    KX, KY = np.meshgrid(np.fft.fftfreq(N, 1 / N), np.fft.rfftfreq(N, 1 / N), indexing="ij")
    K2 = KX**2 + KY**2
    deal = (np.abs(KX) < N / 3) & (np.abs(KY) < N / 3)
    F = np.fft.rfft2
    B = lambda g: np.fft.irfft2(g, s=(N, N))
    inv = np.where(K2 > 0, 1 / np.where(K2 > 0, K2, 1), 0)
    wh = K2 * F(psi_f(X, Y)) * deal
    ah = F(a_f(X, Y)) * deal
    wt = np.ones_like(K2)
    wt[:, 1:-1] = 2  # rfft half-plane weights (N even)

    def mean2(gh):
        return float((wt * np.abs(gh) ** 2).sum() / N**4)

    def rhs(wh, ah):
        ph = wh * inv
        ux, uy = B(1j * KY * ph), B(-1j * KX * ph)
        bx, by = B(1j * KY * ah), B(-1j * KX * ah)
        jh = K2 * ah
        nw = -(ux * B(1j * KX * wh) + uy * B(1j * KY * wh)) + bx * B(1j * KX * jh) + by * B(1j * KY * jh)
        na = -(ux * B(1j * KX * ah) + uy * B(1j * KY * ah))
        return F(nw) * deal, F(na) * deal

    Ew, Ea = np.exp(-nu * K2 * dt / 2), np.exp(-eta * K2 * dt / 2)
    out = []

    def diag(t, wh, ah):
        jh = K2 * ah
        Ev, Em = 0.5 * mean2(np.sqrt(K2) * wh * inv), 0.5 * mean2(np.sqrt(K2) * ah)
        Z, J2 = mean2(wh), mean2(jh)
        out.append((t, Ev, Em, nu * Z + eta * J2, float(np.abs(B(jh)).max()), Z, J2))

    diag(0.0, wh, ah)
    for s in range(int(round(T / dt))):
        k1w, k1a = rhs(wh, ah)
        k2w, k2a = rhs(Ew * (wh + dt / 2 * k1w), Ea * (ah + dt / 2 * k1a))
        k3w, k3a = rhs(Ew * wh + dt / 2 * k2w, Ea * ah + dt / 2 * k2a)
        k4w, k4a = rhs(Ew * Ew * wh + dt * Ew * k3w, Ea * Ea * ah + dt * Ea * k3a)
        wh = Ew * Ew * wh + dt / 6 * (Ew * Ew * k1w + 2 * Ew * (k2w + k3w) + k4w)
        ah = Ea * Ea * ah + dt / 6 * (Ea * Ea * k1a + 2 * Ea * (k2a + k3a) + k4a)
        if (s + 1) % every == 0:
            diag((s + 1) * dt, wh, ah)
    return np.array(out)  # columns: t, Ev, Em, D, max|j| (grid), <omega^2>, <j^2>


def local_maxima(o, col):
    """(t, value) of each interior local maximum of column col, refined by a parabola."""
    t, y, r = o[:, 0], o[:, col], []
    for i in range(1, len(y) - 1):
        if y[i] > y[i - 1] and y[i] >= y[i + 1]:
            d = (y[i - 1] - y[i + 1]) / (2 * (y[i - 1] - 2 * y[i] + y[i + 1]))
            r.append((t[i] + d * (t[1] - t[0]), y[i] - 0.25 * (y[i - 1] - y[i + 1]) * d))
    return r


def at(o, t):
    return o[np.argmin(abs(o[:, 0] - t))]


cos = np.cos
OT_PSI = lambda X, Y: cos(X) + cos(Y)          # u = (-sin y, sin x)
OT_A = lambda X, Y: cos(Y) + 0.5 * cos(2 * X)  # b = (-sin y, sin 2x)

print("== 1a. Orszag-Tang initial condition (the issue's normalisation)")
o = ot_run(64, 0.0, 0.0, OT_PSI, OT_A, 0.0, 1.0)
t0 = o[0]
check("E(0) = 1, E_v = E_m = 1/2", close(t0[1], 0.5, 1e-12) and close(t0[2], 0.5, 1e-12))
check("<omega^2>(0) = 1, <j^2>(0) = 5/2", close(t0[5], 1, 1e-12) and close(t0[6], 2.5, 1e-12))
check("max|j|(0) = 3 (j = cos y + 2 cos 2x)", close(t0[4], 3, 1e-12))
# cross-helicity: u.b = sin^2 y + sin x sin 2x -> <u.b> = 1/2; 2<u.b>/<u^2+b^2> = 1/2
check("2<u.b>/<u^2+b^2> = 0.5 (Rosenberg et al. print 41% for their IC; not reproduced)", True, "0.5")

print("== 1b. Rosenberg, Pouquet & Mininni 2007 IC: which reading of the printed a_z")
P2 = lambda X, Y: 2 * (cos(X) + cos(Y))
A_lit = lambda X, Y: 2 * cos(X) + cos(Y)          # as printed: max|j|(0) = 3
A_v1 = lambda X, Y: 2 * cos(X) + cos(2 * Y)       # max|j|(0) = 6
A_v2 = lambda X, Y: cos(2 * X) + 2 * cos(Y)       # = 2 x the issue's IC; max|j|(0) = 6
oL = ot_run(64, 0.025, 0.025, P2, A_lit, 0.0, 1.0)
check("printed a_z = 2cos x + cos y gives max|j|(0) = 3, not Fig. 8's ~6 -> misprint", close(oL[0, 4], 3, 1e-12))
o1 = ot_run(64, 0.025, 0.025, P2, A_v1, 1.0, 2.5e-3, every=4)
o2 = ot_run(64, 0.025, 0.025, P2, A_v2, 1.0, 2.5e-3, every=4)
check("readings 2cos x + cos 2y and cos 2x + 2cos y give identical histories",
      np.abs(o1[:, 1:] - o2[:, 1:]).max() < 1e-9, f"max diff {np.abs(o1[:, 1:] - o2[:, 1:]).max():.1e}")
check("both have E(0) = 4 (paper: E_V = E_M = 2) and max|j|(0) = 6",
      close(o1[0, 1] + o1[0, 2], 4, 1e-12) and close(o1[0, 4], 6, 1e-12))

print("== 1c. Rosenberg et al. Fig. 8 (nu = eta = 0.025, 128^2 pseudo-spectral), reproduced")
oR = ot_run(128, 0.025, 0.025, P2, A_v2, 2.5, 2.5e-3, every=2)
pk = max(local_maxima(oR, 4), key=lambda p: p[1])
check("peak max|j| = 45.05 at t = 1.073 (computed)", close(pk[1], 45.05, 0.01) and close(pk[0], 1.073, 0.002),
      f"{pk[1]:.3f} at {pk[0]:.4f}")
# figure-read values (paper, Fig. 8 black curve), read by eye at the stated tolerances
check("Fig. 8: peak ~45 near t ~1.1 (+-2, +-0.07)", close(pk[1], 45, 2) and close(pk[0], 1.1, 0.07))
check("Fig. 8: ~29 at t = 1.5-1.6 (+-1.5)", close(at(oR, 1.5)[4], 29, 1.5) and close(at(oR, 1.6)[4], 29, 1.5),
      f"{at(oR, 1.5)[4]:.2f}, {at(oR, 1.6)[4]:.2f}")
check("Fig. 8: ~23 at t = 2.5 (+-1.5)", close(at(oR, 2.5)[4], 23, 1.5), f"{at(oR, 2.5)[4]:.2f}")
dpk = max(local_maxima(oR, 3), key=lambda p: p[1])
check("text: 'peak in total dissipation is reached near t = 1.1' -> computed 1.051", close(dpk[0], 1.051, 0.002),
      f"D peak {dpk[1]:.4f} at {dpk[0]:.4f}")

print("== 1d. the issue's IC at nu = eta = 0.0125 (Fig. 8 in the issue's units: t x2, j /2)")
oA = ot_run(128, 0.0125, 0.0125, OT_PSI, OT_A, 5.0, 2.5e-3, every=2)
oB = ot_run(256, 0.0125, 0.0125, OT_PSI, OT_A, 5.0, 1.25e-3, every=4)
mA, mB = local_maxima(oA, 4), local_maxima(oB, 4)
pA, pB = max(mA, key=lambda p: p[1]), max(mB, key=lambda p: p[1])
check("peak max|j| = 22.527 at t = 2.1458 (128^2)", close(pA[1], 22.527, 2e-3) and close(pA[0], 2.1458, 5e-4),
      f"{pA[1]:.4f} at {pA[0]:.4f}")
check("peak max|j| = 22.534 at t = 2.1459 (256^2)", close(pB[1], 22.534, 2e-3) and close(pB[0], 2.1459, 5e-4),
      f"{pB[1]:.4f} at {pB[0]:.4f}")
check("scaling identity: 2 x 128^2 peak = amplitude-2 run's peak (45.05) to 1e-3",
      close(2 * pA[1], pk[1], 1e-3 * pk[1]) and close(pA[0] / 2, pk[0], 1e-3), f"{2 * pA[1]:.4f} vs {pk[1]:.4f}")
check("other local maxima: 18.15 at 1.835 and 14.75 at 3.229 (128^2)",
      any(close(v, 18.15, 0.01) and close(t, 1.835, 0.003) for t, v in mA)
      and any(close(v, 14.75, 0.01) and close(t, 3.229, 0.003) for t, v in mA))
dA = max(local_maxima(oA, 3), key=lambda p: p[1])
check("dissipation peak D = 0.14522 at t = 2.1025", close(dA[1], 0.14522, 1e-5) and close(dA[0], 2.1025, 5e-4),
      f"{dA[1]:.5f} at {dA[0]:.4f}")
for t, E in ((1, 0.943781), (2, 0.821857), (3, 0.697867), (4, 0.618952), (5, 0.559084)):
    check(f"E({t}) = {E} (128^2 and 256^2 agree to 1e-6)",
          close(at(oA, t)[1] + at(oA, t)[2], E, 1e-6) and close(at(oB, t)[1] + at(oB, t)[2], E, 1e-6))
eb = np.abs(np.interp(oB[:, 0], oA[:, 0], oA[:, 1] + oA[:, 2]) - oB[:, 1] - oB[:, 2]).max()
check("energy-identity proxy: 128^2 and 256^2 energy histories differ < 1e-6", eb < 1e-6, f"{eb:.1e}")

print("== 1e. the issue's IC at nu = eta = 0.05 and 0.025, 128^2 (64^2 for the resolution check)")
REF = {  # nu: (peak max|j| (t, value), D peak (t, value), {t: E})
    0.05: ((2.0226, 7.5144), (1.1494, 0.22001), {1: 0.813019, 2: 0.599845, 4: 0.354156}),
    0.025: ((2.0884, 13.7375), (1.9549, 0.17857), {1: 0.895126, 2: 0.726973, 4: 0.499248}),
}
for nu, (jp, dp, Es) in REF.items():
    o = ot_run(128, nu, nu, OT_PSI, OT_A, 4.0, 2.5e-3, every=4)
    o64 = ot_run(64, nu, nu, OT_PSI, OT_A, 4.0, 2.5e-3, every=4)
    p = max(local_maxima(o, 4), key=lambda q: q[1])
    p64 = max(local_maxima(o64, 4), key=lambda q: q[1])
    d = max(local_maxima(o, 3), key=lambda q: q[1])
    check(f"nu={nu}: peak max|j| = {jp[1]} at t = {jp[0]}", close(p[1], jp[1], 2e-4) and close(p[0], jp[0], 2e-4),
          f"{p[1]:.5f} at {p[0]:.5f}; 64^2 gives {p64[1]:.5f}")
    check(f"nu={nu}: D peak = {dp[1]} at t = {dp[0]}", close(d[1], dp[1], 1e-5) and close(d[0], dp[0], 2e-4))
    for t, E in Es.items():
        check(f"nu={nu}: E({t}) = {E}", close(at(o, t)[1] + at(o, t)[2], E, 1e-6))

print("== 1f. Rosenberg et al. Fig. 9 (256^2 series): which viscosity")
for nu, want in ((0.01, ((0.954, 71.52), (1.114, 76.51), (1.671, 63.50))), (0.0125, ((0.945, 61.4), (1.103, 67.97), (1.66, 54.07)))):
    o = ot_run(256, nu, nu, P2, A_v2, 1.8, 1.25e-3, every=4)
    m = local_maxima(o, 4)
    ok = all(any(close(t, wt_, 0.01) and close(v, wv, 0.1) for t, v in m) for wt_, wv in want)
    check(f"amplitude-2 IC, nu = eta = {nu}: max|j| local maxima {want}", ok,
          str([(round(float(t), 3), round(float(v), 2)) for t, v in m]))
    if nu == 0.01:
        check("Fig. 9 (read by eye): ~71 at ~0.95, ~76.5 at ~1.1, ~63 at ~1.65 -> matches nu = 0.01 (+-1.5)",
              all(any(close(t, a, 0.06) and close(v, b, 1.5) for t, v in m) for a, b in ((0.95, 71), (1.1, 76.5), (1.65, 63))))
        check("Fig. 10: peak <omega^2>+<j^2> ~82 near t ~1.1 (+-2)", close((o[:, 5] + o[:, 6]).max(), 82, 2),
              f"{(o[:, 5] + o[:, 6]).max():.2f}")


# ---------------------------------------------------------------------------------------------
# 2. Kinematic dynamo, d_t b = curl(u x b) + eta lap b, 2*pi box, Rm = 1/eta, mean field excluded.
#    Fourier-Galerkin on |k_x|,|k_y|,|k_z| <= K (k != 0), Leray-projected; gradient modes pushed
#    to -10 so they never lead.  Eigenvalue lambda = sigma + i omega of the field (energy: 2 sigma).
#    u = (A sin z + C cos y, B sin x + A cos z, C sin y + B cos x): A = B = C = 1 is the issue's ABC;
#    A = 0, B = C = 1 is a Roberts-type flow (z-independent).
# ---------------------------------------------------------------------------------------------
def cross_mat(a):
    z = np.zeros(a.shape[0], a.dtype)
    return np.stack([np.stack([z, -a[:, 2], a[:, 1]], -1), np.stack([a[:, 2], z, -a[:, 0]], -1),
                     np.stack([-a[:, 1], a[:, 0], z], -1)], 1)


def uhat(A, B, C):
    return {(0, 0, 1): (-0.5j * A, 0.5 * A, 0), (0, 0, -1): (0.5j * A, 0.5 * A, 0),
            (0, 1, 0): (0.5 * C, 0, -0.5j * C), (0, -1, 0): (0.5 * C, 0, 0.5j * C),
            (1, 0, 0): (0, -0.5j * B, 0.5 * B), (-1, 0, 0): (0, 0.5j * B, 0.5 * B)}


def induction_matrix(K, eta, abc=(1, 1, 1)):
    r = np.arange(-K, K + 1)
    ks = np.stack([g.ravel() for g in np.meshgrid(r, r, r, indexing="ij")], 1)
    ks = ks[np.any(ks != 0, 1)]
    n = len(ks)
    idx = -np.ones((2 * K + 1,) * 3, int)
    idx[ks[:, 0] + K, ks[:, 1] + K, ks[:, 2] + K] = np.arange(n)
    k2 = (ks**2).sum(1).astype(float)
    P = np.eye(3)[None] - ks[:, :, None] * ks[:, None, :] / k2[:, None, None]
    rows, cols, vals = [], [], []

    def add(ri, ci, blk):
        for a in range(3):
            for b in range(3):
                rows.append(3 * ri + a)
                cols.append(3 * ci + b)
                vals.append(blk[:, a, b])

    add(np.arange(n), np.arange(n), -eta * k2[:, None, None] * P - 10 * (np.eye(3)[None] - P))
    for q, uq in uhat(*abc).items():
        if not any(uq):
            continue
        src = ks - np.array(q)
        ok = np.all(np.abs(src) <= K, 1) & np.any(src != 0, 1)
        ri, ci = np.nonzero(ok)[0], idx[src[ok, 0] + K, src[ok, 1] + K, src[ok, 2] + K]
        M = cross_mat(1j * ks[ok].astype(complex)) @ cross_mat(np.tile(np.array(uq, complex), (len(ri), 1)))
        add(ri, ci, M @ P[ci])
    A = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(3 * n, 3 * n))
    return A, ks


def leading(K, Rm, abc=(1, 1, 1), nev=8):
    A, _ = induction_matrix(K, 1 / Rm, abc)
    w = sla.eigs(A, k=nev, which="LR", ncv=40, tol=1e-12, maxiter=20000, return_eigenvectors=False)
    return w[np.argsort(-w.real)]


print("== 2a. operator sanity")
x = np.array([0.3, 1.7, 4.1])
u_pt = np.array([math.sin(x[2]) + math.cos(x[1]), math.sin(x[0]) + math.cos(x[2]), math.sin(x[1]) + math.cos(x[0])])
s = sum(np.array(v) * np.exp(1j * np.dot(q, x)) for q, v in uhat(1, 1, 1).items())
check("Fourier coefficients reproduce u = (sin z + cos y, sin x + cos z, sin y + cos x)", np.abs(s - u_pt).max() < 1e-14)
A, ks = induction_matrix(6, 1 / 12)
uvec = np.zeros(A.shape[0], complex)
for q, v in uhat(1, 1, 1).items():
    i = np.nonzero(np.all(ks == np.array(q), 1))[0][0]
    uvec[3 * i:3 * i + 3] = v
check("b = u (Beltrami) is an exact eigenmode with lambda = -eta = -1/12", np.abs(A @ uvec + uvec / 12).max() < 1e-14)


def _job(args):
    K, Rm, abc = args
    return args, leading(K, Rm, abc)


def run_jobs(jobs):
    from multiprocessing import Pool
    with Pool(min(8, len(jobs))) as p:
        return dict(p.map(_job, jobs))


R111, ROB = (1, 1, 1), (0, 1, 1)
# brentq (xtol 2e-3) gave roots 8.934 (K = 8), 17.554 (K = 10), 25.547 (K = 10); brackets of +-0.01 to 0.05:
EDGES = [(8, 8.92), (8, 8.95), (10, 17.54), (10, 17.57), (10, 25.5), (10, 25.6)]
if __name__ == "__main__":
    jobs = [(5, 5, R111), (8, 5, R111), (5, 12, R111), (8, 12, R111), (8, 20, R111), (5, 20, R111), (8, 8.5, R111),
            (8, 8.26, R111), (8, 11.55, R111), (8, 14.43, R111), (8, 19.23, R111),
            (9, 24.05, R111), (9, 24.10, R111), (10, 26, R111), (10, 27, R111)]
    jobs += [(K, Rm, R111) for K, Rm in EDGES]
    jobs += [(8, Rm, ROB) for Rm in (2, 3, 5, 6, 7, 10, 20)] + [(8, 6.205, ROB), (8, 6.225, ROB)]
    L = run_jobs(jobs)
    lam = lambda K, Rm, abc=R111: L[(K, Rm, abc)]

    print("== 2b. 1:1:1 ABC: leading eigenvalues (field amplitude; energy rate is 2 sigma)")
    w = lam(8, 12)
    check("Rm = 12: lambda = 0.00807 +- 0.61773 i (complex pair, not real)",
          close(w[0].real, 0.00807, 1e-5) and close(abs(w[0].imag), 0.61773, 1e-5) and close(w[1], w[0].conjugate(), 1e-9),
          f"{w[0]:.6f}")
    check("Rm = 12: next eigenvalue -1/12 (the b = u mode)", close(w[2].real, -1 / 12, 1e-9))
    w5 = lam(5, 12)
    check("Rm = 12 at K = 5 (16^3 with 2/3 rule): 0.008069 +- 0.61854 i",
          close(w5[0].real, 0.008069, 2e-6) and close(abs(w5[0].imag), 0.61854, 2e-5), f"{w5[0]:.6f}")
    for K in (5, 8):
        w = lam(K, 5)
        check(f"Rm = 5, K = {K}: lambda = -0.08479 +- 0.61614 i, next -0.2",
              close(w[0].real, -0.08479, 1e-5) and close(abs(w[0].imag), 0.61614, 1e-5) and close(w[2].real, -0.2, 1e-9))
    w = lam(8, 20)
    # K = 9 (scan, in the notes) gives -0.003031 +- 0.56435 i and -0.025445 +- 0.133307 i
    check("Rm = 20, K = 8: -0.00309 +- 0.5643 i, then -0.02529 +- 0.1334 i",
          close(w[0].real, -0.00309, 5e-6) and close(abs(w[0].imag), 0.5643, 1e-4)
          and close(w[2].real, -0.02529, 1e-5) and close(abs(w[2].imag), 0.1334, 1e-4), f"{w[0]:.6f}, {w[2]:.6f}")
    check("Rm = 20, K = 5 (16^3): under-resolved, sigma = -0.0109", close(lam(5, 20)[0].real, -0.0109, 1e-4))
    check("Rm = 8.5: sigma = -0.0033 < 0 (Dorch's printed lower edge 8.5 is not an edge)",
          close(lam(8, 8.5)[0].real, -0.0033, 5e-5))

    print("== 2c. window edges and the mode crossing")
    (Ka, a1), (Kb, a2), (Kc, b1), (Kd, b2), (Ke, c1), (Kf, c2) = EDGES
    check(f"first window opens between Rm = {a1} and {a2}", lam(Ka, a1)[0].real < 0 < lam(Kb, a2)[0].real)
    check(f"first window closes between Rm = {b1} and {b2}", lam(Kc, b1)[0].real > 0 > lam(Kd, b2)[0].real)
    check(f"second window opens between Rm = {c1} and {c2}", lam(Ke, c1)[0].real < 0 < lam(Kf, c2)[0].real)
    check("published 8.9 < Rm < 17.5 lies inside the computed first window", a2 <= 8.95 and b1 >= 17.5)
    check("published second onset 'Rm >= 27' is late: sigma(26) = +0.0015 at K = 10",
          close(lam(10, 26)[0].real, 0.00154, 3e-5), f"{lam(10, 26)[0].real:.6f}")
    check("Rm = 27, K = 10: sigma = 0.00475, omega = 0.1275",
          close(lam(10, 27)[0].real, 0.00475, 3e-5) and close(abs(lam(10, 27)[0].imag), 0.1275, 2e-4))
    w1, w2 = lam(9, 24.05)[0], lam(9, 24.10)[0]
    check("Bouya & Dormy Rm1 in [24.05, 24.10]: leading omega jumps 0.5385 -> 0.1308",
          close(abs(w1.imag), 0.5385, 2e-4) and close(abs(w2.imag), 0.1308, 2e-4), f"{w1:.6f} -> {w2:.6f}")

    print("== 2d. Alexakis (2011) statements, R_M = sqrt(3) Rm")
    check("U = <u.u>^(1/2) = sqrt(3), k_u = 1 for the issue's flow", True, "each Beltrami pair has <.^2> = 1")
    for RM, sign in ((14.3, -1), (20, 1), (25, 1), (33.3, -1)):
        Rm = round(RM / math.sqrt(3), 2)
        check(f"R_M = {RM} -> Rm = {Rm}: dynamo {'yes' if sign > 0 else 'no'}", sign * lam(8, Rm)[0].real > 0)

    print("== 2e. Roberts-type flow A = 0, B = C = 1: real, doubly degenerate")
    for Rm, s in ((2, -0.47474), (3, -0.24964), (5, -0.05622), (6, -0.00820), (7, 0.02539), (10, 0.08289), (20, 0.13738)):
        w = lam(8, Rm, ROB)
        check(f"Rm = {Rm}: lambda = {s} (real, double)",
              close(w[0].real, s, 1e-5) and abs(w[0].imag) < 1e-9 and close(w[0], w[1], 1e-9), f"{w[0].real:.6f}")
    check("onset between Rm = 6.205 and 6.225 (root 6.215)", lam(8, 6.205, ROB)[0].real < 0 < lam(8, 6.225, ROB)[0].real)

    print()
    print(f"{len(FAILS)} failure(s)" if FAILS else "all checks pass")
    sys.exit(1 if FAILS else 0)
