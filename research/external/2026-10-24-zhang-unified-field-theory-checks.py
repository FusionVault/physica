"""Checks for the review of research/external/2026-10-24-zhang-unified-field-theory.md.
Standard library only. Run from outside the repository: python3 -I <path>. Takes a few seconds."""
from math import pi, sqrt, log, sin, e as EULER
EPS0 = 8.8541878128e-12; MU0 = 1.25663706212e-6; C = 299792458.0; G = 6.67430e-11; H = 6.62607015e-34
QE = 1.602176634e-19; KB = 1.380649e-23; G0 = 9.80665; GM_EARTH = 3.986004418e14
K = 1 / (4 * pi * EPS0)

def elliptic(m):
    a, b, c = 1.0, sqrt(1 - m), sqrt(m); s, w = 0.5 * m, 0.5
    for _ in range(64):
        if c <= 2.2e-16 * a: break
        a, b, c = 0.5 * (a + b), sqrt(a * b), 0.5 * (a - b); w *= 2; s += w * c * c
    k = pi / (2 * a); return k, k * (1 - s)

def ring_phi(q, R, z0, rho, z):
    s = (R + rho) ** 2 + (z - z0) ** 2; return K * q * 2 / pi * elliptic(4 * R * rho / s)[0] / sqrt(s)

def ring_E(q, R, z0, rho, z):
    dz = z - z0; s = (R + rho) ** 2 + dz * dz; d = (R - rho) ** 2 + dz * dz; k, e = elliptic(4 * R * rho / s)
    ez = K * q * 2 / pi * dz * e / (sqrt(s) * d)
    er = 0.0 if rho == 0 else K * q / (pi * rho * sqrt(s)) * (k - (R * R - rho * rho + dz * dz) / d * e)
    return er, ez

def solve(rings, wire, pots):
    n = len(rings)
    P = [[(log(8 * Ri / wire) / (4 * pi * pi * EPS0 * Ri) if i == j else ring_phi(1.0, Rj, zj, Ri, zi))
          for j, (Rj, zj) in enumerate(rings)] + [pots[i]] for i, (Ri, zi) in enumerate(rings)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(P[r][c])); P[c], P[p] = P[p], P[c]
        for r in range(c + 1, n):
            f = P[r][c] / P[c][c]; Pr, Pc = P[r], P[c]
            for k in range(c, n + 1): Pr[k] -= f * Pc[k]
    q = [0.0] * n
    for r in range(n - 1, -1, -1):
        q[r] = (P[r][n] - sum(P[r][k] * q[k] for k in range(r + 1, n))) / P[r][r]
    return q

def field(rings, q, rho, z):
    er = ez = 0.0
    for (R, z0), qi in zip(rings, q):
        a, b = ring_E(qi, R, z0, rho, z); er += a; ez += b
    return er, ez

print("== Ring field against direct summation (1 nC, R = 4 cm, at rho 2.3 cm, 2.7 cm above)")
n = 20000; ex = ez = 0.0
for i in range(n):
    t = 2 * pi * (i + 0.5) / n; dx, dy, dz = 0.023 - 0.04 * __import__("math").cos(t), -0.04 * sin(t), 0.027
    d3 = (dx * dx + dy * dy + dz * dz) ** 1.5; ex += K * 1e-9 / n * dx / d3; ez += K * 1e-9 / n * dz / d3
print("  closed form", ring_E(1e-9, 0.04, 0.01, 0.023, 0.037), " summed", (ex, ez))

print("== 1. US 12,417,870: 12-turn coils (R 3.5 cm, 6.5-17.5 cm from the centre), top at 50 kV")
coil = [(0.035, 0.065 + 0.01 * i) for i in range(12)]; rings = coil + [(R, -z) for R, z in coil]
q = solve(rings, 0.75e-3, [50e3] * 12 + [0.0] * 12)
E0 = field(rings, q, 0.0, 0.0)[1]
e2 = lambda z: sum(x * x for x in field(rings, q, 0.0, z))
grad = (e2(1e-4) - e2(-1e-4)) / 2e-4; U = 0.5 * sum(q[:12]) * 50e3
print(f"  charge live {sum(q[:12])*1e9:.2f} nC, grounded {sum(q[12:])*1e9:.2f} nC, energy {U*1e3:.3f} mJ")
print(f"  E at the ball {E0:.4e} V/m, dE2/dz {grad:.4e} V^2/m^3")
m, eps = 0.35e-3, 2.3; a = (3 * m / 950 / (4 * pi)) ** (1 / 3); Kcm = (eps - 1) / (eps + 2); W = m * G0
f_dep = 2 * pi * EPS0 * a ** 3 * Kcm * grad; qs = 12 * pi * EPS0 * a * a * abs(E0) * eps / (eps + 2)
tau = 2 * pi * EPS0 * a ** 3 * E0 ** 2 * (1 - Kcm); I = 0.4 * m * a * a
print(f"  ball radius {a*1e3:.3f} mm; DEP {f_dep*1e6:.3f} uN = {f_dep/W:.2e} W; Pauthenier {qs*1e12:.1f} pC -> {qs*abs(E0)/W:.4f} W")
for T in (10, 30, 60):
    kap = 4 * pi * pi * I / T ** 2; print(f"  torque ceiling {tau:.3e} N m; twist period {T} s -> {tau/kap/2/pi:.0f} turns")
B = MU0 * 100 * 10e-6 / 2 * ((0.12 + 0.06) / sqrt((0.12 + 0.06) ** 2 + 0.035 ** 2) - (0.12 - 0.06) / sqrt((0.12 - 0.06) ** 2 + 0.035 ** 2))
print(f"  B at the ball from 10 uA: {B:.3e} T = {B/50e-6:.1e} of Earth's field")
L = EPS0 * 0.05 * abs(E0) * B * pi * 0.05 ** 2 * 0.2
print(f"  field angular momentum over the tank {L:.3e} kg m^2/s -> spin {L/I:.2e} rad/s")
print(f"  GR: field energy as mass {U/C**2:.3e} kg, pull at 6 cm {G*U/C**2/0.06**2:.3e} m/s^2")
vd = 10e-6 / (8.47e28 * QE * pi * 0.75e-3 ** 2)
print(f"  own terms: drift velocity {vd:.3e} m/s; at 30 kHz: {2*pi*3e4*vd:.2e} m/s^2, {vd/(2*pi*3e4):.1e} m")

print("== 2. The gas in a 10 cm tank (air, A = 15/(cm Torr), B = 365 V/(cm Torr), gamma = 0.01)")
TC = 133.322368 * 0.01; A_, B_ = 15 / TC, 365 / TC; lnln = log(log(1 + 1 / 0.01))
pdmin = EULER * log(1 + 1 / 0.01) / A_
print(f"  Paschen minimum {B_*pdmin:.2f} V at {pdmin:.4f} Pa m ({pdmin/TC:.3f} Torr cm); no breakdown below {log(101)/A_/0.1:.2f} Pa over 10 cm")
for p in (1e5, 1e3, 100, 10, 4.5, 1, 0.1, 0.01):
    den = log(A_ * p * 0.1) - lnln; vb = B_ * p * 0.1 / den if den > 0 else float("inf")
    lam = KB * 293.15 / (sqrt(2) * pi * 3.711e-10 ** 2 * p)
    print(f"  {p:8g} Pa: mean free path {lam:.3e} m, Townsend breakdown {vb:.4g} V")
print(f"  Peek onset, 0.75 mm wire: {3e6*(1+0.0301/sqrt(0.75e-3)):.4e} V/m; surface field at 50 kV ~ {50e3/(0.75e-3*log(0.2/0.75e-3)):.2e} V/m")

print("== 3. The nylon ring on a balance: wire 1.5 mm, 50 mm gap, 55 kV, ring r 20-50 mm, 50 mm tall")
def weighing(top, bottom, rin=0.020, rout=0.050, R=0.75e-3, wire=0.75e-3, n=200, step=1.5e-3):
    up = [(R, 0.025 + step / 2 + step * i) for i in range(n)]; rings = up + [(r, -z) for r, z in up]
    q = solve(rings, wire, [top] * n + [bottom] * n)
    s, N = 0.0, 60; dr = (rout - rin) / N
    for i in range(N + 1):
        rho = rin + dr * i; w = 0.5 if i in (0, N) else 1.0
        e_t = sum(x * x for x in field(rings, q, rho, 0.025)); e_b = sum(x * x for x in field(rings, q, rho, -0.025))
        s += w * (e_t - e_b) * 2 * pi * rho * dr
    return 0.5 * EPS0 * s, q[n // 2] / step  # upward force per unit (eps_r - 1), line charge mid-wire
for name, t, b in (("top live", 55e3, 0.0), ("bottom live", 0.0, 55e3), ("+-27.5 kV", 27.5e3, -27.5e3)):
    f, lam = weighing(t, b)
    print(f"  {name:12s} weight change {-f*2.6/G0*1e3:+.3f} g (dilute) / {-f*0.722/G0*1e3:+.3f} g (normal field); line charge {lam*1e9:.0f} nC/m")
print(f"  isolated finite wire for comparison: 2 pi eps0 V/ln(L/a) = {2*pi*EPS0*55e3/log(0.3/0.75e-3)*1e9:.0f} nC/m")
f, _ = weighing(55e3, 0.0, rin=0.011)
print(f"  ring reaching to 11 mm (5 mm from the cable): {-f*2.6/G0*1e3:+.3f} / {-f*0.722/G0*1e3:+.3f} g")
cage1, _ = weighing(55e3, 0.0, n=100, step=3e-3)
cage7, _ = weighing(55e3, 0.0, rin=0.025, R=0.013, n=100, step=3e-3)
print(f"  seven live cables as a cage of radius 13 mm, ring from 25 mm: {cage7/cage1:.2f} x the single cable")
fwd, rev = (-0.83 - 1.20) / 2, (0.76 + 1.08) / 2; fwd7, rev7 = (-6.17 - 8.36) / 2, (5.47 + 7.63) / 2
print(f"  reported: odd part {(fwd-rev)/2:+.4f} g, even part {(fwd+rev)/2:+.4f} g; seven cables odd {(fwd7-rev7)/2:+.4f} g (x{(fwd7-rev7)/(fwd-rev):.2f})")
ring_mass = 1140 * pi * (0.05 ** 2 - 0.02 ** 2) * 0.05
print(f"  ring mass {ring_mass:.3f} kg; claimed field {1e-3*G0/ring_mass:.4f} m/s^2 per gram ({1e-3/ring_mass*1e6:.0f} ppm of g)")
print(f"  own terms: 30 uA ramped in 1 ms -> carrier acceleration {30e-6/1e-3/(8.47e28*QE*pi*0.75e-3**2):.3e} m/s^2; steady DC -> 0")
q_needed = 0.9675e-3 * G0 / 2e5; area = 2 * pi * 0.05 * 0.05 + 2 * pi * 0.02 * 0.05 + 2 * pi * (0.05 ** 2 - 0.02 ** 2)
print(f"  a net charge instead: {q_needed*1e9:.0f} nC at 2e5 V/m, {q_needed/area:.2e} C/m^2 = {q_needed/area/(EPS0*3e6):.2f} of air breakdown")

print("== 4. The air gap of a 0.5 T core: objects on cotton threads")
bgrad = 0.5 * 0.5 / 0.018
for name, vol, mass, chi in (("leaf (water)", 0.038 * 0.014 * 0.30e-3, 125e-6, -9.04e-6), ("pork skin (water)", 0.020 * 0.010 * 3e-3, 0.5e-3, -9.04e-6),
                             ("paper (cellulose ~ water)", 0.037 * 0.0125 * 0.11e-3, 50e-6, -9.04e-6)):
    F = abs(chi) * vol * bgrad / MU0; print(f"  {name:26s} susceptibility force {F:.2e} N = {F/(mass*G0):.3f} of its weight")
delta = sqrt(2 * 2.65e-8 / (2 * pi * 50 * MU0)); ratio = MU0 * (1 / 2.65e-8) * 2 * pi * 50 * 0.002 * 0.012 / 4
print(f"  aluminium at 50 Hz: skin depth {delta*1e3:.1f} mm; induced/applied field ~ {ratio:.3f}; repulsion ~ {ratio**2*0.04*0.012*0.002*bgrad/MU0:.2e} N vs weight {2e-3*G0:.2e} N")

print("== 5. The Michelson fringe shift: Maxwell stress")
print(f"  eps0 E^2/2 at 1e6 V/m: {0.5*EPS0*1e12:.3f} Pa; on 1 cm^2 {0.5*EPS0*1e12*1e-4:.2e} N; GR path change ~ {2*G*U/C**4/0.05*0.1:.1e} m")

print("== 6. Zhang's E_theta/e_r = A x R / c^2 is the kink field ratio")
qq, acc, r, th = 1e-9, 1e12, 1.0, pi / 2
print(f"  E_rad = {qq*acc*sin(th)/(4*pi*EPS0*C*C*r):.6e} V/m; e_r * a r sin(th)/c^2 = {K*qq/r**2*acc*r*sin(th)/C**2:.6e} V/m")

print("== 7. The five-constant equation is Kepler's third law")
mass = GM_EARTH * (1 + 1 / 81.30057) / G; nu = mass * C * C / H; r, T = 384399e3, 27.321661 * 86400
print(f"  Moon: 4 pi^2 r^3 c^2/(T^2 G h nu) = {4*pi*pi*r**3*C*C/(T*T*G*H*nu):.5f}  (nu = {nu:.3e} Hz)")
print("  GW170817 bounds the speed of gravity: -3e-15 < (v - c)/c < 7e-16 (gravitational::propagation::speed_bounds)")
