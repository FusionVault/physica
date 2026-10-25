"""Gap 5 follow-ups: (a) k2(4/9, y=2) closed form vs DN09; (b) cancellation in Eq. 23 at small C and
Postnikov's series (checked against the exact formula in 60-digit Decimal arithmetic); (c) Newtonian
k2 for n = 0.5 and n = 1 polytropes from Lane-Emden + Newtonian Riccati (Postnikov Eq. 11)."""
from decimal import Decimal as D, getcontext
import numpy as np
from scipy.integrate import solve_ivp
from gap5_love_tov import k2_of
getcontext().prec = 60

def ln1m(x):  # ln(1 - x) by series, Decimal
    s, t, n = D(0), D(1), 1
    while True:
        t *= x; term = t / n; s -= term; n += 1
        if abs(term) < D(10) ** -58: return s

def k2_exact(C, y):
    C, y = D(C), D(y)
    num = D(8) * C**5 / 5 * (1 - 2 * C) ** 2 * (2 + 2 * C * (y - 1) - y)
    den = (2 * C * (6 - 3 * y + 3 * C * (5 * y - 8)) + 4 * C**3 * (13 - 11 * y + C * (3 * y - 2) + 2 * C**2 * (1 + y))
           + 3 * (1 - 2 * C) ** 2 * (2 - y + 2 * C * (y - 1)) * ln1m(2 * C))
    return num / den

print("(a) k2(4/9, y=2): float", k2_of(4/9, 2.0), " Decimal", k2_exact(D(4)/9, 2), " DN09:", 4096/(10935*(308-81*np.log(3))))
print("(b) cancellation: float Eq.23 vs 60-digit exact, y = (pi^2-9)/3 (n=1 Newtonian):")
y = (np.pi**2 - 9) / 3
for C in [1e-5, 1e-4, 1e-3, 1e-2, 0.05]:
    ex = float(k2_exact(C, y))
    Y = y + 3
    t = [(2 - y) / Y, (y*y - 6*y - 6) / Y**2, (y**3 + 34*y**2 - 8*y + 12) / (7 * Y**3),
         (y**4 + 62*y**3 + 84*y**2 + 48*y + 36) / (7 * Y**4),
         5 / 294 * (5*y**5 + 490*y**4 + 1272*y**3 + 1884*y**2 + 1476*y + 648) / Y**5]
    ser_as_printed = (1 - 2*C)**2 / 2 * (t[0] + t[1]*C + t[2]*C**2 + t[3]*C**4 + t[4]*C**5)
    ser_beta3 = (1 - 2*C)**2 / 2 * (t[0] + t[1]*C + t[2]*C**2 + t[3]*C**3 + t[4]*C**4)
    print(f"   C={C:g}: exact {ex:.12f}  float rel err {k2_of(C, y)/ex-1:+.1e}  series(as printed, C^4,C^5) {ser_as_printed/ex-1:+.1e}  series(C^3,C^4) {ser_beta3/ex-1:+.1e}")

# fit the true series coefficients from the exact formula: k2/((1-2C)^2/2) = sum b_k C^k
Cs = [D(k) / 10**6 for k in range(1, 8)]
import itertools
vals = [k2_exact(c, y) * 2 / (1 - 2 * c) ** 2 for c in Cs]
# solve Vandermonde in Decimal (7x7)
n = len(Cs); A = [[c**j for j in range(n)] + [v] for c, v in zip(Cs, vals)]
for i in range(n):
    p = max(range(i, n), key=lambda r: abs(A[r][i])); A[i], A[p] = A[p], A[i]
    for r in range(n):
        if r != i:
            f = A[r][i] / A[i][i]; A[r] = [a - f * b for a, b in zip(A[r], A[i])]
b = [A[i][n] / A[i][i] for i in range(n)]
Y = y + 3
print("   true series coefficients b0..b4:", [f"{float(x):.10f}" for x in b[:5]])
print("   Postnikov terms t0..t4         :", [f"{x:.10f}" for x in [(2 - y) / Y, (y*y - 6*y - 6) / Y**2, (y**3 + 34*y**2 - 8*y + 12) / (7 * Y**3), (y**4 + 62*y**3 + 84*y**2 + 48*y + 36) / (7 * Y**4), 5/294*(5*y**5+490*y**4+1272*y**3+1884*y**2+1476*y+648)/Y**5]])

# (c) Newtonian polytrope k2: Lane-Emden theta'' = -theta^n - 2 theta'/xi; y Riccati: xi y' + y^2 + y - 6 + xi^2 * (4 pi rho / (dp/drho)) ...
def newt_k2(n):
    # In Lane-Emden units 4 pi G rho/(dp/drho) * r^2 = xi^2 * theta^n / ((n+1)/n * theta) * ... = xi^2 * n * theta^(n-1) / 1  (see notes)
    def rhs(x, u):
        th, dth, yy = u
        th = max(th, 1e-300)
        return [dth, -th**n - 2 * dth / x, (-(yy * yy + yy - 6) - x * x * n * th ** (n - 1)) / x]
    ev = lambda x, u: u[0]; ev.terminal = True; ev.direction = -1
    x0 = 1e-6
    s = solve_ivp(rhs, [x0, 20], [1 - x0**2 / 6, -x0 / 3, 2.0], events=ev, rtol=1e-12, atol=1e-14, method="DOP853")
    yR = s.y[2, -1]
    return (2 - yR) / (2 * (3 + yR)), s.t[-1]
for n in [1.0, 0.5, 0.3, 0.7, 1.2]:
    k, xi1 = newt_k2(n)
    print(f"(c) Newtonian n={n}: k2 = {k:.5f}, xi_1 = {xi1:.5f}")
print("    Hinderer 2008 erratum Table 1 at M/R = 1e-5: n=0.3 0.5511, 0.5 0.4491, 0.7 0.3626, 1.0 0.2599, 1.2 0.2062")
