"""A general metric -> stress-energy -> energy-condition evaluator (the Warp Factory recipe, pointwise),
validated on Schwarzschild, Alcubierre, Natario, Morris-Thorne, and the SSV closed forms.

Run: python3 metric_to_energy_conditions.py   (numpy only; reads no downloaded files)
Units G = c = 1 unless stated. Signature (-,+,+,+). T = G_{mu nu}/(8 pi).
"""
import numpy as np

H = 1e-3  # finite-difference step (4th-order central differences, nested once)


def d4(fun, X, a, h=None):
    h = H if h is None else h
    e = np.zeros(4); e[a] = h
    return (8 * (fun(X + e) - fun(X - e)) - (fun(X + 2 * e) - fun(X - 2 * e))) / (12 * h)


def christoffel(g, X):
    gi = np.linalg.inv(g(X))
    dg = np.array([d4(g, X, a) for a in range(4)])          # dg[a, m, n] = d_a g_mn
    s = dg.transpose(1, 0, 2) + dg.transpose(1, 2, 0) - dg  # [s, m, n] = d_m g_sn + d_n g_sm - d_s g_mn
    return 0.5 * np.einsum("ls,smn->lmn", gi, s)


def einstein(g, X):
    G4 = christoffel(g, X)
    dG = np.array([d4(lambda Y: christoffel(g, Y), X, a) for a in range(4)])   # dG[a, l, m, n]
    Ric = (np.einsum("llmn->mn", dG) - np.einsum("nllm->mn", dG)
           + np.einsum("lls,smn->mn", G4, G4) - np.einsum("lns,sml->mn", G4, G4))
    gx = g(X)
    return Ric - 0.5 * np.einsum("mn,mn", np.linalg.inv(gx), Ric) * gx, gx


def eulerian_frame(gx):
    gi = np.linalg.inv(gx)
    alpha = np.sqrt(-1 / gi[0, 0])
    e = [-alpha * gi[:, 0]]                                  # n^mu = g^{mu nu} n_nu with n_nu = (-alpha,0,0,0)
    for i in (1, 2, 3):                                      # coordinate d_i are tangent to t = const: Gram-Schmidt
        u = np.zeros(4); u[i] = 1.0
        for b in e[1:]:
            u -= (u @ gx @ b) * b
        e.append(u / np.sqrt(u @ gx @ u))
    return np.array(e)


SPHERE = (lambda k: np.stack([np.sqrt(1 - (1 - 2 * (k + .5) / 6000)**2) * np.cos(np.pi * (1 + 5**.5) * k),
                              np.sqrt(1 - (1 - 2 * (k + .5) / 6000)**2) * np.sin(np.pi * (1 + 5**.5) * k),
                              1 - 2 * (k + .5) / 6000], 1))(np.arange(6000))


def conditions(g, X):
    """Eulerian rho, mean pressure, flux, and minima of NEC (over null k with k^0=1), WEC (boosts up to 0.99c), SEC_E, DEC_E."""
    Gmn, gx = einstein(g, X)
    E = eulerian_frame(gx)
    T = E @ (Gmn / (8 * np.pi)) @ E.T                        # orthonormal components T_{a b}
    rho, flux, S = T[0, 0], T[0, 1:], T[1:, 1:]
    nec = np.min(rho + 2 * SPHERE @ flux + np.einsum("ki,ij,kj->k", SPHERE, S, SPHERE))
    wec = min(np.min((rho + 2 * b * SPHERE @ flux + b * b * np.einsum("ki,ij,kj->k", SPHERE, S, SPHERE)) / (1 - b * b))
              for b in np.linspace(0, 0.99, 34))
    return dict(rho=rho, pbar=np.trace(S) / 3, flux=flux, nec=nec, wec=wec,
                sec=0.5 * (rho + np.trace(S)), dec=rho - np.linalg.norm(flux), T=T, Tcoord=Gmn / (8 * np.pi))


# ---------- metrics ----------
def f_alc(r, R, s):
    return (np.tanh(s * (r + R)) - np.tanh(s * (r - R))) / (2 * np.tanh(s * R))


def df_alc(r, R, s):
    return s * (1 / np.cosh(s * (r + R))**2 - 1 / np.cosh(s * (r - R))**2) / (2 * np.tanh(s * R))


def shift_metric(vfield):
    """ds^2 = -dt^2 + sum (dx^i - v^i dt)^2 (SSV/Natario convention: shift beta^i = -v^i, lapse 1, flat slices)."""
    def g(X):
        vv = vfield(X)
        m = np.eye(4); m[0, 0] = -1 + vv @ vv; m[0, 1:] = m[1:, 0] = -vv
        return m
    return g


def alcubierre(v, R, s):
    return shift_metric(lambda X: np.array([v * f_alc(np.sqrt((X[1] - v * X[0])**2 + X[2]**2 + X[3]**2), R, s), 0, 0]))


def alcubierre_comoving(v, R, s):
    """Warp Factory's metricGet_AlcubierreComoving: g_0x = beta_x = v(1 - f), static."""
    def g(X):
        b = v * (1 - f_alc(np.sqrt(X[1]**2 + X[2]**2 + X[3]**2), R, s))
        m = np.eye(4); m[0, 0] = -1 + b * b; m[0, 1] = m[1, 0] = b
        return m
    return g


def natario_X(v, R, s):
    """Natario's zero-expansion field (polar axis x): X = -2 v f cos(th) e_r + v (2f + r f') sin(th) e_th, f = (1 - f_alc)/2."""
    def X(Y):
        x, y, z = Y[1:]
        r = np.sqrt(x * x + y * y + z * z); rho = np.sqrt(y * y + z * z)
        f, fp = 0.5 * (1 - f_alc(r, R, s)), -0.5 * df_alc(r, R, s)
        ct, st = x / r, rho / r
        er = np.array([x, y, z]) / r
        eth = np.array([-st, ct * y / rho, ct * z / rho])
        return -2 * v * f * ct * er + v * (2 * f + r * fp) * st * eth
    return X


def natario_rho(v, R, s, x, y, z):
    r = np.sqrt(x * x + y * y + z * z); h = 1e-5
    fp = -0.5 * df_alc(r, R, s); fpp = -0.5 * (df_alc(r + h, R, s) - df_alc(r - h, R, s)) / (2 * h)
    ct2 = (x / r)**2
    return -v * v / (8 * np.pi) * (3 * fp * fp * ct2 + (fp + r * fpp / 2)**2 * (1 - ct2))


def schwarzschild_isotropic(M):
    def g(X):
        r = np.linalg.norm(X[1:]); a = M / (2 * r)
        return np.diag([-((1 - a) / (1 + a))**2] + [(1 + a)**4] * 3)
    return g


def morris_thorne(r0):
    """Phi = 0, b = r0^2/r, upper sheet r > r0, Cartesian: g_ij = delta_ij + (1/(1-b/r) - 1) x_i x_j / r^2."""
    def g(X):
        xyz = X[1:]; r = np.linalg.norm(xyz)
        m = np.eye(4); m[0, 0] = -1
        m[1:, 1:] += (1 / (1 - r0 * r0 / r**2) - 1) * np.outer(xyz, xyz) / r**2
        return m
    return g


def ssv_rho_plus_pbar(vfield, X, sign=+1):
    """SSV 2022 Eq. (nec-1): rho + pbar = (1/24 pi)(-2 L_n K + K^2 - 3 tr K^2), K_ij = sign * v_(i,j), L_n = d_t + v^i d_i."""
    def Kmat(Y):
        J = np.array([d4(lambda Z: vfield(Z), Y, a) for a in (1, 2, 3)])   # J[j, i] = d_j v^i
        return sign * 0.5 * (J + J.T)
    K = Kmat(X); trK = lambda Y: np.trace(Kmat(Y))
    LnK = d4(trK, X, 0) + vfield(X) @ np.array([d4(trK, X, a) for a in (1, 2, 3)])
    return (-2 * LnK + np.trace(K)**2 - 3 * np.trace(K @ K)) / (24 * np.pi), (np.trace(K)**2 - np.trace(K @ K)) / (16 * np.pi)


if __name__ == "__main__":
    print("--- Schwarzschild (isotropic, M = 1): T should vanish ---")
    for r in (3.0, 6.0):
        c = conditions(schwarzschild_isotropic(1.0), np.array([0, r / np.sqrt(3)] * 1 + [r / np.sqrt(3)] * 2))
        print(f"r = {r}: max |T_ab| = {np.max(np.abs(c['T'])):.2e} (curvature scale 1/r^2 = {1/r**2:.2e})")

    v, R, s = 1.0, 1.0, 8.0
    print(f"\n--- Alcubierre (lab frame, t = 0), v = {v}, R = {R}, sigma = {s} ---")
    vf = lambda X: np.array([v * f_alc(np.sqrt((X[1] - v * X[0])**2 + X[2]**2 + X[3]**2), R, s), 0, 0])
    for name, P in (("wall, side (theta=90)", (0, 0, 1.0, 0)), ("wall, 45 deg", (0, .7071, .7071, 0)),
                    ("wall, front on axis", (0, 1.0, 1e-3, 0)), ("wall, rear on axis", (0, -1.0, 1e-3, 0)),
                    ("inside", (0, 0.2, 0.1, 0)), ("outside", (0, 0, 1.6, 0))):
        X = np.array(P, float); r = np.linalg.norm(X[1:])
        c = conditions(alcubierre(v, R, s), X)
        ana = -v * v / (32 * np.pi) * df_alc(r, R, s)**2 * (X[2]**2 + X[3]**2) / r**2
        rp_ssv, rho_ssv = ssv_rho_plus_pbar(vf, X, +1)
        rp_ssv_m, _ = ssv_rho_plus_pbar(vf, X, -1)
        print(f"{name:22s} rho = {c['rho']: .5e} (analytic {ana: .5e}, SSV {rho_ssv: .5e}); "
              f"rho+pbar = {c['rho'] + c['pbar']: .4e} (SSV K=+v_(i,j): {rp_ssv: .4e}, K=-v_(i,j): {rp_ssv_m: .4e}); "
              f"NEC min {c['nec']: .3e}; WEC min(<=0.99c) {c['wec']: .3e}; SEC_E {c['sec']: .3e}; DEC_E {c['dec']: .3e}")

    print("\n--- Warp Factory test case in SI: comoving Alcubierre v = 0.1, R = 300 m, sigma = 0.015 /m ---")
    C4G = 299792458.0**4 / 6.67430e-11
    H = 0.5  # metres: wall width ~ 1/sigma = 67 m
    c = conditions(alcubierre_comoving(0.1, 300.0, 0.015), np.array([0, 0, 300.0, 0]))
    print(f"rho at (0, 300, 0) = {c['rho'] * C4G: .4e} J/m^3 ; analytic {-C4G/(32*np.pi)*0.01*df_alc(300.0,300.0,0.015)**2: .4e}; "
          f"NEC min {c['nec'] * C4G: .3e} J/m^3")

    H = 1e-3
    v, R, s = 1.0, 1.0, 4.0
    print(f"\n--- Natario zero expansion (static), v = {v}, R = {R}, sigma = {s} ---")
    gN = shift_metric(natario_X(v, R, s))
    for P in ((0, 0.6, 0.85, 0.0), (0, 1.0, 0.1, 0.0), (0, 0.1, 1.0, 0.0), (0, -0.7, 0.6, 0.3)):
        X = np.array(P, float)
        c = conditions(gN, X)
        print(f"x,y,z = {P[1:]}: rho = {c['rho']: .5e} (Natario formula {natario_rho(v, R, s, *P[1:]): .5e}); "
              f"expansion-free check; NEC min {c['nec']: .3e}; WEC min {c['wec']: .3e}")

    print("\n--- Morris-Thorne b = r0^2/r, r0 = 1 ---")
    for r in (1.2, 2.0):
        X = np.array([0, r / np.sqrt(2), r / np.sqrt(2), 0.0])
        c = conditions(morris_thorne(1.0), X)
        er = np.r_[0, X[1:] / r] * np.sqrt(1 - 1 / r**2)  # unit radial vector: g_rr = 1/(1 - b/r)
        pr = er @ c["Tcoord"] @ er
        print(f"r = {r}: rho = {c['rho']: .5e} (analytic {-1/(8*np.pi*r**4): .5e}); p_r = {pr: .5e} "
              f"(analytic {-1/(8*np.pi*r**4): .5e}); rho+p_r = {c['rho'] + pr: .4e} (analytic {-2/(8*np.pi*r**4): .4e}); NEC min {c['nec']: .4e}")
