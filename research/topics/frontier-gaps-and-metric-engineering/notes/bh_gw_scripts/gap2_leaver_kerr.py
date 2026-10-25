"""Gap 2: Kerr QNMs by Leaver's continued fraction (radial) + Cook-Zalutskiy spectral method (angular).

Written from the equations printed in Cook & Zalutskiy 2014 (arXiv:1410.7698): radial Teukolsky
equation -> confluent Heun form with (zeta_+, xi_-, eta_+), D_0..D_4 (their Eq. 31), the 3-term
recurrence alpha_n, beta_n, gamma_n, the n-th inversion Cf(n;N) and the Nollert tail r_N; angular
equation -> pentadiagonal matrix of their Eqs. (52)-(55). Units M = 1, s = -2.

Checks: Berti's Schwarzschild table s2l2.dat (2M = 1 units), qnm README (a = 0.68), BCW06 Tables II-IV.
Reads one downloaded data file (plain numbers) as data only.
"""
import sys
import numpy as np

S = -2


# ---------- angular: Cook & Zalutskiy Eqs. (51)-(55) ----------
def F(s, l, m):
    return np.sqrt(((l + 1) ** 2 - m * m) / ((2 * l + 3) * (2 * l + 1))) * np.sqrt(((l + 1) ** 2 - s * s) / (l + 1) ** 2)


def G(s, l, m):
    return 0.0 if l == 0 else np.sqrt((l * l - m * m) / (4 * l * l - 1)) * np.sqrt((l * l - s * s) / (l * l))


def H(s, l, m):
    return 0.0 if (l == 0 or s == 0) else -m * s / (l * (l + 1))


def sep_const(c, l, m, s=S, lmax_extra=25):
    """Eigenvalue A_{slm}(c) of the spectral matrix, the one continuous with l(l+1)-s(s+1) at c=0."""
    lmin = max(abs(m), abs(s))
    ls = np.arange(lmin, l + lmax_extra + 1)
    n = len(ls)
    M = np.zeros((n, n), complex)
    for i, L in enumerate(ls):
        for j, Lp in enumerate(ls):
            if Lp == L - 2:
                M[i, j] = -c * c * F(s, Lp, m) * F(s, Lp + 1, m)
            elif Lp == L - 1:
                M[i, j] = -c * c * F(s, Lp, m) * (H(s, Lp + 1, m) + H(s, Lp, m)) + 2 * c * s * F(s, Lp, m)
            elif Lp == L:
                B = F(s, Lp, m) * G(s, Lp + 1, m) + G(s, Lp, m) * F(s, Lp - 1, m) + H(s, Lp, m) ** 2
                M[i, j] = Lp * (Lp + 1) - s * (s + 1) - c * c * B + 2 * c * s * H(s, Lp, m)
            elif Lp == L + 1:
                M[i, j] = -c * c * G(s, Lp, m) * (H(s, Lp - 1, m) + H(s, Lp, m)) + 2 * c * s * G(s, Lp, m)
            elif Lp == L + 2:
                M[i, j] = -c * c * G(s, Lp, m) * G(s, Lp - 1, m)
    ev = np.linalg.eigvals(M)
    return ev


def A_closest(c, l, m, A_guess):
    ev = sep_const(c, l, m)
    return ev[np.argmin(abs(ev - A_guess))]


# ---------- radial: Cook & Zalutskiy Sec. II C 1 ----------
def cf(omega, a, l, m, n_inv, A, N=600, s=S):
    root = np.sqrt(1 - a * a)
    rp, rm = 1 + root, 1 - root
    sp = (2 * omega * rp - m * a) / (rp - rm)
    sm = (2 * omega * rm - m * a) / (rp - rm)
    zeta = 1j * omega                      # zeta_+
    xi = -s - 1j * sp                      # xi_-
    eta = -1j * sm                         # eta_+
    p = (rp - rm) * zeta / 2
    al = 1 + s + xi + eta - 2 * zeta + s * 1j * omega / zeta
    ga = 1 + s + 2 * eta
    de = 1 + s + 2 * xi
    sig = A + a * a * omega * omega - 8 * omega * omega + p * (2 * al + ga - de) + (1 + s - (ga + de) / 2) * (s + (ga + de) / 2)
    D0, D1, D2, D3, D4 = de, 4 * p - 2 * al + ga - de - 2, 2 * al - ga + 2, al * (4 * p - de) - sig, al * (al - ga + 1)
    k = np.arange(N + 1)
    alpha = k * k + (D0 + 1) * k + D0
    beta = -2 * k * k + (D1 + 2) * k + D3
    gamma = k * k + (D2 - 3) * k + D4 - D2 + 2
    # Nollert tail r_N = 1 + u1/sqrt(N) + u2/N + u3/N^1.5, minimal branch Re(u1) < 0
    u1 = np.sqrt(-4 * p + 0j)
    if u1.real > 0:
        u1 = -u1
    u2 = -(8 * p - 4 * al + 2 * ga + 2 * de + 3) / 4
    u3 = (32 * p * (2 * p - 4 * al + ga + 3 * de + 4) + 4 * (ga + de) * (ga + de - 2) + 16 * sig + 3) / (32 * u1)
    rN = 1 + u1 / np.sqrt(N) + u2 / N + u3 / N ** 1.5
    e = beta[N] + alpha[N] * rN
    for j in range(N - 1, n_inv, -1):
        e = beta[j] - alpha[j] * gamma[j + 1] / e
    Q = alpha[n_inv] * gamma[n_inv + 1] / e
    P = 0
    if n_inv > 0:
        d = beta[0]
        for j in range(1, n_inv):
            d = beta[j] - alpha[j - 1] * gamma[j] / d
        P = alpha[n_inv - 1] * gamma[n_inv] / d
    return beta[n_inv] - P - Q


def solve(a, l, m, n, omega_guess, A_guess, N=600, tol=1e-14):
    """Complex secant on f(omega) = Cf(n; N; omega, A(a omega)), A re-solved at every omega (Cook's fix)."""
    state = {"A": A_guess}

    def f(w):
        state["A"] = A_closest(a * w, l, m, state["A"])
        return cf(w, a, l, m, n, state["A"], N)

    w0, w1 = omega_guess, omega_guess * (1 + 1e-6) + 1e-7
    f0, f1 = f(w0), f(w1)
    for _ in range(100):
        w2 = w1 - f1 * (w1 - w0) / (f1 - f0)
        w0, f0, w1 = w1, f1, w2
        f1 = f(w1)
        if abs(w1 - w0) < tol:
            break
    return w1, state["A"]


def track(l, m, n, a_targets, omega0, N=600, da=0.01):
    """Follow a mode from a = 0 (Schwarzschild omega0) to each target spin."""
    out = {}
    A = l * (l + 1) - S * (S + 1) + 0j
    w, a = omega0, 0.0
    w, A = solve(0.0, l, m, n, w, A, N)
    hist = [(0.0, w)]
    for at in sorted(a_targets):
        while a < at - 1e-12:
            step = min(da, at - a)
            a_new = a + step
            if len(hist) >= 3:  # quadratic extrapolation (as in Cook & Zalutskiy) for the next guess
                (x0, y0), (x1, y1), (x2, y2) = hist[-3:]
                g = (y0 * (a_new - x1) * (a_new - x2) / ((x0 - x1) * (x0 - x2)) + y1 * (a_new - x0) * (a_new - x2) / ((x1 - x0) * (x1 - x2))
                     + y2 * (a_new - x0) * (a_new - x1) / ((x2 - x0) * (x2 - x1)))
            else:
                g = w
            w, A = solve(a_new, l, m, n, g, A, N)
            a = a_new
            hist.append((a, w))
        out[at] = (w, A)
    return out


if __name__ == "__main__":
    s2l2 = sys.argv[1]  # Berti's s2l2.dat (Schwarzschild, 2M = 1): Re, Im, error, n
    berti = np.loadtxt(s2l2)
    print("Schwarzschild l=2 vs Berti s2l2.dat (M omega = (2M omega)/2):")
    for n in range(8):
        ref = (berti[n, 0] - 1j * berti[n, 1]) / 2
        w, A = solve(0.0, 2, 2, n, ref * (1 + 2e-3), 4 + 0j, N=1500)
        print(f"  n={n}: {w.real:.14f} {w.imag:+.14f}i   ref {ref.real:.14f} {ref.imag:+.14f}i   |diff| {abs(w-ref):.1e}")
    for l, nn in [(3, 0), (4, 0)]:
        pass

    # Kerr, (2,2,0) and others vs qnm README and BCW06 tables
    w220 = (berti[0, 0] - 1j * berti[0, 1]) / 2
    spins = [0.1, 0.3, 0.5, 0.68, 0.7, 0.8, 0.9, 0.98, 0.99]
    res = track(2, 2, 0, spins, w220)
    print("(l,m,n) = (2,2,0):")
    for a in spins:
        w, A = res[a]
        print(f"  a={a}: M omega = {w.real:.12f} {w.imag:+.12f}i   A = {A.real:.10f} {A.imag:+.10f}i  Q = {w.real/(-2*w.imag):.8f}")
    qnm_readme = 0.5239751042900845 - 0.08151262363119974j
    print(f"  qnm README a=0.68: |diff| = {abs(res[0.68][0]-qnm_readme):.2e}")

    for (l, m, n) in [(2, 2, 1), (2, 1, 0), (2, -2, 0), (3, 3, 0), (4, 4, 0)]:
        if l == 2:
            w0 = (berti[n, 0] - 1j * berti[n, 1]) / 2
        else:
            w0 = {3: 0.5994 - 0.0927j, 4: 0.8092 - 0.0942j}[l]  # BCW06 Tables III/IV j=0, refined by the solver
        r = track(l, m, n, [0.5, 0.7, 0.9, 0.98], w0)
        print(f"({l},{m},{n}):", "  ".join(f"a={a}: {r[a][0].real:.10f}{r[a][0].imag:+.10f}i" for a in [0.5, 0.7, 0.9, 0.98]))
