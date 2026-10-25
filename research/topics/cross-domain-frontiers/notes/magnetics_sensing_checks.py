#!/usr/bin/env python3
"""Checks behind magnetics_and_quantum_sensing.md. Standard library only; run as `python3 -I`.

Every `check` compares a recomputed value with a published one (source in the notes); every
`ident` is an exact identity the proposed Rust function must satisfy; `show` prints a value that
has no published counterpart. Reads two trusted files from this repository: the bundled
WMM2025.COF and the IGRF-14 block of the earlier geomagnetism_checks.py.
"""
import cmath
import os
import re
from math import acos, asin, atan2, cos, degrees, exp, log, pi, radians, sin, sqrt

MU0 = 1.25663706212e-6        # CODATA 2018, as in the metrology crate
HBAR = 1.054571817e-34; H = 6.62607015e-34; QE = 1.602176634e-19; MUB = 9.2740100783e-24
G = 6.67430e-11
CM = MU0 / (4 * pi)           # magnetic constant C_m of Blakely (1995) in SI, 1e-7 T m/A
A_REF = 6371200.0
WGS_A = 6378137.0; WGS_F = 1 / 298.257223563; WGS_E2 = WGS_F * (2 - WGS_F)
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../.."))

FAILS = []
def check(name, got, pub, tol, unit=""):
    ok = abs(got - pub) <= tol
    if not ok:
        FAILS.append(name)
    print(f"{'ok ' if ok else 'BAD'} {name:<66s} {got:13.6g} {unit:<7s} pub {pub:<11.6g} diff {got - pub:+.3g}")

def ident(name, err, tol=1e-9):
    ok = abs(err) <= tol
    if not ok:
        FAILS.append(name)
    print(f"{'ok ' if ok else 'BAD'} {name:<66s} |err| {abs(err):.2e} (tol {tol:g})")

def show(name, v, unit=""):
    print(f"    {name:<66s} {v:13.6g} {unit}")

# ============================================================== 1. quantum magnetometry
print("== 1. Energy resolution, NV shot noise, spin-projection noise ==")
def energy_resolution(db, volume, t=1.0):
    """E_R/hbar = dB^2 V T / (2 mu0 hbar); db in T/sqrt(Hz) with t = 1 s, or single-shot T with t."""
    return db * db * volume * t / (2 * MU0 * HBAR)
check("Dang et al. 2010: 160 aT/sqrt(Hz), 0.45 cm3 (existing crate test)", energy_resolution(160e-18, 0.45e-6), 44, 0.6, "hbar")
check("Palacios Alvarez et al. 2022 BEC: 72 fT, 1091 um3, 3.5 s", energy_resolution(72e-15, 1091e-18, 3.5), 0.075, 0.002, "hbar")
show("Kominis et al. 2003 SERF: 0.54 fT/sqrt(Hz) in 0.3 cm3", energy_resolution(0.54e-15, 0.3e-6), "hbar")
show("Sekiguchi et al. 2024 NV CW-ODMR: 9.4 pT/sqrt(Hz) in 0.03 mm3", energy_resolution(9.4e-12, 0.03e-9), "hbar")

gamma_nv = 2.003 * MUB / H                         # Hz/T, g_e = 2.003 (Barry et al. 2020)
show("NV electron gyromagnetic ratio g mu_B / h (g = 2.003)", gamma_nv / 1e9, "GHz/T")
check("Sekiguchi 2024: slope 332 pA/Hz x gamma -> A/T", 332e-12 * gamma_nv, 9.3, 0.05, "A/T")
# shot noise of the photocurrent through a lock-in: sqrt(2 x 2 q I) (their Eq. 2), divided by the slope
check("Sekiguchi 2024 photon-shot-noise limit sqrt(4 q I)/S, 6.4 mA, 9.3 A/T", sqrt(4 * QE * 6.4e-3) / 9.3 * 1e12, 6.9, 0.05, "pT/rtHz")
eta_cw = 4 / (3 * sqrt(3)) / gamma_nv * 0.48e6 / (0.03 * sqrt(6.4e-3 / QE))
show("Barry et al. 2020 Eq. 16 with dnu 0.48 MHz, C 3 %, R = I/q (NOT 6.9)", eta_cw * 1e12, "pT/rtHz")

def spin_projection(gamma, n_atoms, t2, t=1.0):
    """dB = 1/(gamma sqrt(N T2 t)), gamma in rad s^-1 T^-1 (Budker & Romalis 2007 Eq. 1)."""
    return 1 / (gamma * sqrt(n_atoms * t2 * t))
ge = 2 * pi * 2.0023 * MUB / H                     # free-electron gamma, rad/(s T)
for q in (1, 6):                                   # q = 6: SERF slowing-down factor, I = 3/2, low polarisation (assumed)
    t2 = 1 / ((ge / q) ** 2 * 1e14 * (1e-17) ** 2)  # T2 that puts N = 1e14 (1 cm3 at 1e14 cm^-3) at 1e-17 T/rtHz
    show(f"T2 implied by Budker-Romalis 1e-17 T/rtHz per cm3, q = {q}", t2 * 1e3, "ms")
ident("spin_projection inverse: dB(N, T2) recovers 1e-17 T/rtHz", spin_projection(ge, 1e14, 1 / (ge ** 2 * 1e14 * 1e-34)) - 1e-17, 1e-27)
check("Lucivero et al. 2022: 14 fT/rtHz over 0.2 cm baseline", 14 / 0.2, 70, 0.01, "fT/cm/rtHz")
# gravimetry::local::averaging_time(asd, target) = (asd/target)^2 applied to magnetometers
show("Sekiguchi 2024: white-noise time to 0.3 pT from 9.4 pT/rtHz (paper: Allan, 'a few thousand s')", (9.4 / 0.3) ** 2, "s")
ident("Silani et al. 2023: 70 fT s^1/2 averaged 1 h is below the published 2 fT", 0.0 if 70 / sqrt(3600) < 2 else 1.0)
show("   the white-noise value after 1 h", 70 / sqrt(3600), "fT")

# atom interferometer <-> Ramsey magnetometer: both 1/(C sqrt N) per shot over a scale factor
K_EFF = 4 * pi / 780.241e-9
check("gravimetry::atom doc: k_eff T^2 at T = 0.16 s (Rosi et al.)", K_EFF * 0.16 ** 2 / 1e5, 4.12, 0.005, "1e5 rad/(m/s2)")
gamma_rb = 2 * pi * (2.0023 / 4) * MUB / H         # 87Rb F = 2: g_F = g_s/(2I+1) ~ 1/2 (nuclear term neglected)
show("87Rb F=2 gamma/2pi (g_F = 1/2)", gamma_rb / 2 / pi / 1e9, "GHz/T")
n_at, c_at, tau = 1e6, 1.0, 0.32
show("same N=1e6, C=1: atom gravimeter da per shot, T=0.16 s, double diffr.", 1 / (c_at * sqrt(n_at)) / (2 * K_EFF * 0.16 ** 2), "m/s2")
show("same N=1e6, C=1: 87Rb Ramsey dB per shot, tau = 2T = 0.32 s", 1 / (c_at * sqrt(n_at)) / (gamma_rb * tau), "T")

# ============================================================== 2. aeromagnetic compensation and MagNav
print("\n== 2. Tolles-Lawson (Gnadt, Wollaber & Nielsen 2022) and MagNav arithmetic ==")
def solve(a, b):
    """Gaussian elimination with partial pivoting; a is n x n (list of lists), b length n."""
    n = len(b); m = [row[:] + [b[i]] for i, row in enumerate(a)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(m[r][c])); m[c], m[p] = m[p], m[c]
        for r in range(c + 1, n):
            f = m[r][c] / m[c][c]
            for k in range(c, n + 1):
                m[r][k] -= f * m[c][k]
    x = [0.0] * n
    for r in range(n - 1, -1, -1):
        x[r] = (m[r][n] - sum(m[r][k] * x[k] for k in range(r + 1, n))) / m[r][r]
    return x

def lstsq(rows, y, ridge=0.0):
    n = len(rows[0])
    ata = [[sum(r[i] * r[j] for r in rows) + (ridge if i == j else 0.0) for j in range(n)] for i in range(n)]
    aty = [sum(r[i] * yy for r, yy in zip(rows, y)) for i in range(n)]
    return solve(ata, aty)

def rot_body_from_ned(roll, pitch, yaw):
    cr, sr, cp, sp, cy, sy = cos(roll), sin(roll), cos(pitch), sin(pitch), cos(yaw), sin(yaw)
    return [[cp * cy, cp * sy, -sp],
            [sr * sp * cy - cr * sy, sr * sp * sy + cr * cy, sr * cp],
            [cr * sp * cy + sr * sy, cr * sp * sy - sr * cy, cr * cp]]

def matvec(m, v):
    return [sum(m[i][k] * v[k] for k in range(3)) for i in range(3)]

def tl_row(b, bdot_hat):
    """18 Tolles-Lawson terms: 3 permanent, 6 induced (|B| c_i c_j, i <= j), 9 eddy (|B| c_i dc_j/dt)."""
    f = sqrt(sum(x * x for x in b)); c = [x / f for x in b]
    ind = [f * c[i] * c[j] for i in range(3) for j in range(i, 3)]
    eddy = [f * c[i] * bdot_hat[j] for i in range(3) for j in range(3)]
    return c + ind + eddy

# calibration box: roll, pitch, yaw manoeuvres (periods 5-7 s) on four headings, 10 Hz
dt, Be_ned = 0.1, [17000.0, -3000.0, 50000.0]      # nT, an Ottawa-like field (illustrative)
samples = []
for leg in range(4):
    for k in range(600):
        t = k * dt
        samples.append((radians(10) * sin(2 * pi * t / 5), radians(5) * sin(2 * pi * t / 7),
                        radians(90 * leg) + radians(5) * sin(2 * pi * t / 6)))
Be_body = [matvec(rot_body_from_ned(*s), Be_ned) for s in samples]
def cosines(vs):
    return [[x / sqrt(sum(y * y for y in v)) for x in v] for v in vs]
def dcos(cs, i):
    j0, j1 = max(i - 1, 0), min(i + 1, len(cs) - 1)
    return [(cs[j1][k] - cs[j0][k]) / ((j1 - j0) * dt) for k in range(3)]

# (a) the linear algebra: y = |B_e| + A beta exactly, beta recovered by least squares
beta_true = [120, -80, 40, 0.004, 0.002, -0.003, 0.006, 0.001, -0.002] + [0.5 * (-1) ** i for i in range(9)]
cs = cosines(Be_body)
A = [tl_row(b, dcos(cs, i)) for i, b in enumerate(Be_body)]
y = [sum(a * bb for a, bb in zip(row, beta_true)) for row in A]
beta = lstsq(A, y)
ident("TL: 18 coefficients recovered from noiseless map-based calibration", max(abs(a - b) for a, b in zip(beta, beta_true)), 1e-6)
# (b) the physics: a true aircraft field a + b B_e (b symmetric), scalar |B_e + B_a|, TL fit, residual
a_p = [150.0, -60.0, 90.0]; b_ind = [[0.004, 0.001, -0.0015], [0.001, 0.006, 0.0005], [-0.0015, 0.0005, -0.002]]
Ba = [[a_p[i] + sum(b_ind[i][k] * v[k] for k in range(3)) for i in range(3)] for v in Be_body]
Bt = [[u + w for u, w in zip(v, ba)] for v, ba in zip(Be_body, Ba)]
F_e = sqrt(sum(x * x for x in Be_ned)); scalar = [sqrt(sum(x * x for x in v)) for v in Bt]
cst = cosines(Bt); At = [tl_row(b, dcos(cst, i)) for i, b in enumerate(Bt)]
bt = lstsq(At, [s - F_e for s in scalar])
res = [s - sum(r * c for r, c in zip(row, bt)) - F_e for s, row in zip(scalar, At)]
show("TL: rms aircraft field |B_a| in this simulation", sqrt(sum(sum(x * x for x in v) for v in Ba) / len(Ba)), "nT")
show("TL: rms raw scalar error |B_t| - |B_e|", sqrt(sum((s - F_e) ** 2 for s in scalar) / len(scalar)), "nT")
show("TL: rms residual after 18-term fit", sqrt(sum(r * r for r in res) / len(res)), "nT")
lin = [s - sum(x * c for x, c in zip(ba, cc)) - F_e for s, ba, cc in zip(scalar, Ba, cst)]
sec = [s - sum(x * c for x, c in zip(ba, cc)) + (sum(x * x for x in ba) - sum(x * c for x, c in zip(ba, cc)) ** 2) / (2 * s) - F_e
       for s, ba, cc in zip(scalar, Ba, cst)]
show("TL Eq. 12 first-order error with the true B_a (max)", max(abs(v) for v in lin), "nT")
show("TL Eq. 40 second-order error with the true B_a (max)", max(abs(v) for v in sec), "nT")
check("Q-CTRL 2025: 22 m is 0.006 % of the flight -> distance", 22 / 6e-5 / 1e3, 366.7, 0.1, "km")
check("Canciani 2016/2020: unaided INS 230 m over MagNav 13.1 m DRMS", 230 / 13.1, 17.6, 0.05, "x")
show("MagNav CRLB, 1 axis: sigma_x = sigma_B/|grad T|, 1 nT over 10 nT/km", 1.0 / (10 / 1e3), "m")

# ============================================================== 3. joint gravity-magnetic models
print("\n== 3. Poisson's relation, magnetic sphere/cylinder, Euler, RTP, continuation ==")
def grav_sphere(a, depth, drho, x):          # gravimetry::local::sphere
    dm = 4 / 3 * pi * a ** 3 * drho; s2 = x * x + depth * depth
    return G * dm * depth / s2 ** 1.5, G * dm * (2 * depth ** 2 - x * x) / s2 ** 2.5
def grav_cylinder(a, depth, drho, x):        # gravimetry::local::cylinder
    lam = pi * a * a * drho; s2 = x * x + depth * depth
    return 2 * G * lam * depth / s2, 2 * G * lam * (depth ** 2 - x * x) / s2 ** 2
def unit(inc, dec):                          # x north, y east, z down
    return (cos(inc) * cos(dec), cos(inc) * sin(dec), sin(inc))
def mag_sphere(a, depth, mag, inc_m, dec_m, inc_f, dec_f, x, y=0.0, z=0.0):
    """Total-field anomaly (T) of a uniformly magnetised sphere = dipole m = (4/3)pi a^3 M at (0,0,depth)."""
    m = 4 / 3 * pi * a ** 3 * mag; mh, fh = unit(inc_m, dec_m), unit(inc_f, dec_f)
    r = (x, y, z - depth); rr = sqrt(sum(v * v for v in r)); mr = sum(p * q for p, q in zip(mh, r))
    b = [CM * m * (3 * mr * r[i] / rr ** 5 - mh[i] / rr ** 3) for i in range(3)]
    return sum(p * q for p, q in zip(b, fh))
def mag_cylinder_pole(a, depth, mag, x):
    """2-D horizontal cylinder, vertical magnetisation and field: 2 C_m M pi a^2 (z^2 - x^2)/(x^2+z^2)^2."""
    s2 = x * x + depth * depth
    return 2 * CM * mag * pi * a * a * (depth ** 2 - x * x) / s2 ** 2
check("sphere at the pole: 2 C_m m / z^3, a 5 m, z 20 m, M 1 A/m", mag_sphere(5, 20, 1, pi / 2, 0, pi / 2, 0, 0) * 1e9, 13.09, 0.005, "nT")
for x in (0.0, 7.0, 25.0):
    gz, guu = grav_sphere(5, 20, 2700, x)
    ident(f"Poisson: sphere dT(pole) = C_m M/(G rho) Guu at x = {x} m", mag_sphere(5, 20, 1, pi / 2, 0, pi / 2, 0, x) - CM * 1 / (G * 2700) * guu, 1e-20)
    gz, guu = grav_cylinder(3, 15, 2700, x)
    ident(f"Poisson: cylinder dT(pole) = C_m M/(G rho) Guu at x = {x} m", mag_cylinder_pole(3, 15, 1, x) - CM / (G * 2700) * guu, 1e-20)
gz, guu = grav_sphere(5, 20, -2700, 0)
check("gravimetry::local doc value: 5 m void at 20 m, Guu", guu / 1e-9, -23.6, 0.05, "E")
show("same void if it also removes M = 1 A/m: Poisson dT at the pole", CM * (-1) / (G * -2700) * guu * 1e9 * -1, "nT")
# general direction: dT = (C_m M/(G rho)) d_f d_m U, U = G dM / r, by central differences
def U(px, py, pz, depth=20.0, a=5.0, rho=2700.0):
    return G * 4 / 3 * pi * a ** 3 * rho / sqrt(px * px + py * py + (pz - depth) ** 2)
def dd_u(p, mh, fh, h=1e-3):
    f = lambda s, t: U(*[p[i] + s * mh[i] + t * fh[i] for i in range(3)])
    return (f(h, h) - f(h, -h) - f(-h, h) + f(-h, -h)) / (4 * h * h)
im, dm_, i_f, d_f = radians(60), radians(10), radians(-35), radians(-20)
worst = 0.0
for p in ((0.0, 0.0, 0.0), (12.0, -5.0, 0.0), (-30.0, 8.0, -3.0)):
    poisson = CM * 1 / (G * 2700) * dd_u(p, unit(im, dm_), unit(i_f, d_f))
    direct = mag_sphere(5, 20, 1, im, dm_, i_f, d_f, *p)
    worst = max(worst, abs(poisson - direct) / abs(direct))
ident("Poisson, any inclination/declination: relative error (finite differences)", worst, 1e-5)
check("sphere at the magnetic equator, directly above: -C_m m/z^3", mag_sphere(5, 20, 1, 0, 0, 0, 0, 0) * 1e9, -6.545, 0.001, "nT")
# MDR by regression of dT on Guu (Chandler et al. 1981 moving window), single window here
xs = [i - 50.0 for i in range(101)]
gu = [grav_sphere(5, 20, 2700, x)[1] for x in xs]; dt_ = [mag_sphere(5, 20, 2.5, pi / 2, 0, pi / 2, 0, x) for x in xs]
mg, md = sum(gu) / len(gu), sum(dt_) / len(dt_)
slope = sum((u - mg) * (d - md) for u, d in zip(gu, dt_)) / sum((u - mg) ** 2 for u in gu)
ident("Poisson regression recovers M/rho = 2.5/2700 A m^2/kg", slope * G / CM - 2.5 / 2700, 1e-12)

# Euler homogeneity: (x-x0) df/dx + (z-z0) df/dz = -N f
def euler_n(fun, x, z, h=1e-4):
    fx = (fun(x + h, z) - fun(x - h, z)) / (2 * h); fz = (fun(x, z + h) - fun(x, z - h)) / (2 * h)
    return -(x * fx + (z - 20.0) * fz) / fun(x, z)
cases = (("gravity sphere g_z", lambda x, z: grav_sphere(5, 20 - z, 2700, x)[0], 2),
         ("gravity cylinder g_z", lambda x, z: grav_cylinder(5, 20 - z, 2700, x)[0], 1),
         ("magnetic sphere dT (pole)", lambda x, z: mag_sphere(5, 20, 1, pi / 2, 0, pi / 2, 0, x, 0, z), 3),
         ("magnetic sphere dT (I = 60 deg)", lambda x, z: mag_sphere(5, 20, 1, radians(60), 0, radians(60), 0, x, 0, z), 3),
         ("magnetic cylinder dT (pole)", lambda x, z: mag_cylinder_pole(5, 20 - z, 1, x), 2))
for name, fun, n in cases:
    ident(f"Euler structural index {n}: {name}", euler_n(fun, 9.0, 0.0) - n, 1e-5)
# Euler solve: unknowns x0, z0, background; N known
def euler_solve(fun, n, xs_, h=1e-4):
    rows, rhs = [], []
    for x in xs_:
        f = fun(x, 0.0); fx = (fun(x + h, 0.0) - fun(x - h, 0.0)) / (2 * h); fz = (fun(x, h) - fun(x, -h)) / (2 * h)
        rows.append([fx, fz, n]); rhs.append(x * fx + 0.0 * fz + n * f)
    return lstsq(rows, rhs)
sol = euler_solve(lambda x, z: mag_sphere(5, 20, 1, radians(60), 0, radians(60), 0, x - 3.0, 0, z), 3, [i * 2.0 - 40 for i in range(41)])
ident("Euler deconvolution: dipole depth recovered (20 m)", sol[1] - 20.0, 1e-4)
ident("Euler deconvolution: dipole offset recovered (3 m)", sol[0] - 3.0, 1e-4)

# ---- FFT: reduction to the pole and upward continuation on a 128 x 128 grid
def fft(a, inv=False):
    n = len(a)
    if n == 1:
        return a[:]
    ev, od = fft(a[0::2], inv), fft(a[1::2], inv); s = 1 if inv else -1
    out = [0j] * n
    for k in range(n // 2):
        t = cmath.exp(s * 2j * pi * k / n) * od[k]; out[k] = ev[k] + t; out[k + n // 2] = ev[k] - t
    return out
def fft2(g, inv=False):
    rows = [fft(r, inv) for r in g]; cols = [fft(list(c), inv) for c in zip(*rows)]
    out = [list(r) for r in zip(*cols)]
    if inv:
        nn = len(g) * len(g[0]); out = [[v / nn for v in r] for r in out]
    return out
N, DX, ZD = 128, 10.0, 50.0
ks = [2 * pi * (i if i < N // 2 else i - N) / (N * DX) for i in range(N)]
def grid(fun):
    return [[complex(fun((i - N // 2) * DX, (j - N // 2) * DX)) for j in range(N)] for i in range(N)]  # i: north x, j: east y
inc, dec = radians(45), radians(0)
obs = grid(lambda x, y: mag_sphere(10, ZD, 1, inc, dec, inc, dec, x, y))
pole = grid(lambda x, y: mag_sphere(10, ZD, 1, pi / 2, 0, pi / 2, 0, x, y))
def theta(u, kx, ky):
    k = sqrt(kx * kx + ky * ky)
    return u[2] + 1j * (u[0] * kx + u[1] * ky) / k
def rms_rel(a, b):
    num = sum(abs(a[i][j] - b[i][j]) ** 2 for i in range(N) for j in range(N)); den = sum(abs(b[i][j]) ** 2 for i in range(N) for j in range(N))
    return sqrt(num / den)
spec = fft2(obs); u_m = u_f = unit(inc, dec)
for sign in (1, -1):                                       # pick the sign convention that the FFT kernel implies
    rtp_spec = [[0j if (i == 0 and j == 0) else spec[i][j] / (theta(u_m, sign * ks[i], sign * ks[j]) * theta(u_f, sign * ks[i], sign * ks[j]))
                 for j in range(N)] for i in range(N)]
    err = rms_rel(fft2(rtp_spec, True), pole)
    if sign == 1:
        ident("RTP at I = 45 deg (Blakely sign, numpy-style FFT): rms error vs pole", err, 0.02)
    else:
        show("RTP with the opposite kernel sign (wrong): rms error", err)
show("RTP filter amplification along declination at I = 10 deg: 1/sin^2 I", 1 / sin(radians(10)) ** 2)
ident("RTP filter is 1 at the pole", abs(1 / (theta(unit(pi / 2, 0), 0.3, -0.2) ** 2) - 1), 1e-15)
hup = 30.0
pspec = fft2(pole)
up = fft2([[pspec[i][j] * exp(-sqrt(ks[i] ** 2 + ks[j] ** 2) * hup) for j in range(N)] for i in range(N)], True)
pole_hi = grid(lambda x, y: mag_sphere(10, ZD + hup, 1, pi / 2, 0, pi / 2, 0, x, y))
show("upward continuation by exp(-|k| h), h = 30 m: rms error vs analytic", rms_rel(up, pole_hi))
# pseudo-gravity (Baranov 1957): F[g_z] = (G rho/(C_m M)) F[dT_pole]/|k|; compare peak minus corner
ps = fft2([[0j if (i == 0 and j == 0) else pspec[i][j] / sqrt(ks[i] ** 2 + ks[j] ** 2) * G * 2700 / (CM * 1) for j in range(N)] for i in range(N)], True)
gz_true = grid(lambda x, y: grav_sphere(10, ZD, 2700, sqrt(x * x + y * y))[0])
show("pseudo-gravity: (peak - corner) relative error vs gravimetry::local::sphere",
     ((ps[N // 2][N // 2] - ps[0][0]).real - (gz_true[N // 2][N // 2] - gz_true[0][0]).real) / (gz_true[N // 2][N // 2] - gz_true[0][0]).real)

# cross-gradient (Gallardo & Meju 2003): t = dm1/dx dm2/dz - dm1/dz dm2/dx
def xgrad(m1, m2, x, z, h=1e-5):
    d = lambda f, a, b: ((f(x + h, z) - f(x - h, z)) / (2 * h), (f(x, z + h) - f(x, z - h)) / (2 * h))
    a1, b1 = d(m1, x, z); a2, b2 = d(m2, x, z)
    return a1 * b2 - b1 * a2
m1 = lambda x, z: exp(-((x - 1) ** 2 + (z - 2) ** 2)); m2 = lambda x, z: 3 * m1(x, z) + m1(x, z) ** 2
ident("cross-gradient is zero for m2 = f(m1)", max(abs(xgrad(m1, m2, x, z)) for x in (0.2, 1.5) for z in (1.0, 2.7)), 1e-8)
show("cross-gradient for a shifted body (non-zero)", xgrad(m1, lambda x, z: m1(x - 0.5, z), 0.2, 1.0))

# ============================================================== 4. geomagnetism from WMM2025 and IGRF-14
print("\n== 4. SV timescale, source-radius fit, SAA drift, core waves ==")
def read_cof(path):
    g, h, gd, hd = {}, {}, {}, {}
    with open(path) as f:
        for ln in f.read().split("\n")[1:]:
            p = ln.split()
            if len(p) != 6:
                continue
            k = (int(p[0]), int(p[1])); g[k], h[k], gd[k], hd[k] = map(float, p[2:])
    return g, h, gd, hd
WMM = read_cof(os.path.join(ROOT, "crates/magnetics/data/WMM2025.COF"))
NMAX = 12
def power(g, h, n):
    return (n + 1) * sum(g[(n, m)] ** 2 + h[(n, m)] ** 2 for m in range(n + 1))
R = {n: power(WMM[0], WMM[1], n) for n in range(1, NMAX + 1)}; S = {n: power(WMM[2], WMM[3], n) for n in range(1, NMAX + 1)}
tau = {n: sqrt(R[n] / S[n]) for n in R}
for n in (1, 2, 3, 6, 12):
    show(f"WMM2025 tau_n = sqrt(R_n/S_n), n = {n}", tau[n], "yr")
tau_sv = exp(sum(log(n * tau[n]) for n in range(2, NMAX + 1)) / (NMAX - 1))
check("Lhuillier et al. 2011: tau_SV in tau_n = tau_SV/n, n = 2..12 (WMM2025)", tau_sv, 415, 55, "yr")
def fit_radius(weight):
    xs_ = list(range(2, NMAX + 1)); ys = [log(R[n] / weight(n)) for n in xs_]
    mx, my = sum(xs_) / len(xs_), sum(ys) / len(ys)
    sl = sum((x - mx) * (yy - my) for x, yy in zip(xs_, ys)) / sum((x - mx) ** 2 for x in xs_)
    return A_REF * exp(sl / 2) / 1e3                  # R_n ~ K w(n) (c/a)^(2n+4): slope = 2 ln(c/a)
show("source radius from R_n = K (c/a)^(2n+4), n = 2..12", fit_radius(lambda n: 1.0), "km")
show("source radius from R_n = K (n+1/2) (c/a)^(2n+4), n = 2..12", fit_radius(lambda n: n + 0.5), "km")
show("source radius from R_n = K (n+1/2)/(n(n+1)) (c/a)^(2n+4), n = 2..12", fit_radius(lambda n: (n + 0.5) / (n * (n + 1))), "km")
print("    (Voorhies et al. 2002 published 3512 +- 64 km from other data; seismic CMB 3480 km)")

def schmidt(nmax, t, u):
    P = [[0.0] * (nmax + 1) for _ in range(nmax + 1)]; D = [[0.0] * (nmax + 1) for _ in range(nmax + 1)]
    P[0][0] = 1.0
    for n in range(1, nmax + 1):
        P[n][n] = u * P[n - 1][n - 1] * (1.0 if n == 1 else sqrt((2 * n - 1) / (2 * n)))
        for m in range(n):
            p2 = P[n - 2][m] if n >= 2 else 0.0
            P[n][m] = ((2 * n - 1) * t * P[n - 1][m] - sqrt((n - 1) ** 2 - m * m) * p2) / sqrt(n * n - m * m)
    for n in range(1, nmax + 1):
        for m in range(n + 1):
            D[n][m] = (-n * t * P[n][m] + (sqrt(n * n - m * m) * P[n - 1][m] if m < n else 0.0)) / u
    return P, D
def total_field(g, h, lat_d, lon_d, h_m=0.0):
    phi = radians(lat_d); rc = WGS_A / sqrt(1 - WGS_E2 * sin(phi) ** 2)
    p = (rc + h_m) * cos(phi); z = (rc * (1 - WGS_E2) + h_m) * sin(phi); r = sqrt(p * p + z * z); phic = asin(z / r)
    nmax = max(n for n, _ in g); P, D = schmidt(nmax, sin(phic), cos(phic)); lam = radians(lon_d)
    X = Y = Z = 0.0
    for n in range(1, nmax + 1):
        q = (A_REF / r) ** (n + 2)
        for m in range(n + 1):
            gc = g[(n, m)] * cos(m * lam) + h[(n, m)] * sin(m * lam)
            X -= q * gc * D[n][m]; Y += q * m * (g[(n, m)] * sin(m * lam) - h[(n, m)] * cos(m * lam)) * P[n][m]; Z -= (n + 1) * q * gc * P[n][m]
    return sqrt(X * X + (Y / cos(phic)) ** 2 + Z * Z)
def at(model, t):
    g, h, gd, hd = model
    return {k: g[k] + (t - 2025.0) * gd[k] for k in g}, {k: h[k] + (t - 2025.0) * hd[k] for k in h}
check("WMM2025 Table 6: F at (2025.0, 0 km, 80N, 0E)", total_field(*at(WMM, 2025.0), 80, 0), 55178.5, 0.1, "nT")
def local_min(g, h, la, lo):
    f = total_field(g, h, la, lo); step = 1.0
    while step > 1e-5:
        moved = False
        for dla, dlo in ((step, 0), (-step, 0), (0, step), (0, -step)):
            v = total_field(g, h, la + dla, lo + dlo)
            if v < f:
                f, la, lo, moved = v, la + dla, lo + dlo, True
        if not moved:
            step /= 2
    return f, la, lo
def km_between(la1, lo1, la2, lo2):
    a1, a2, dl = radians(la1), radians(la2), radians(lo2 - lo1)
    return 6371.0 * acos(min(1.0, sin(a1) * sin(a2) + cos(a1) * cos(a2) * cos(dl)))
def igrf14():
    txt = open(os.path.join(ROOT, "research/topics/magnetics-and-mhd/notes/geomagnetism_checks.py")).read()
    block = re.search(r'IGRF14 = """(.*?)"""', txt, re.S).group(1)
    g20, h20, g25, h25 = {}, {}, {}, {}
    for ln in block.split("\n"):
        p = ln.split()
        if len(p) < 5:
            continue
        k = (int(p[1]), int(p[2]))
        (g20 if p[0] == "g" else h20)[k] = float(p[3]); (g25 if p[0] == "g" else h25)[k] = float(p[4])
    for d in (h20, h25):
        for k in g25:
            d.setdefault(k, 0.0)
    return (g20, h20), (g25, h25)
(G20, H20), (G25, H25) = igrf14()
f25, la25, lo25 = local_min(G25, H25, -26.0, -60.0)
check("IGRF-14 2025 SAA minimum F (existing report value)", f25, 22071, 0.6, "nT")
f20, la20, lo20 = local_min(G20, H20, -26.0, -58.0)
show("IGRF-14 2020 SAA minimum F", f20, "nT"); show("   at latitude", la20, "deg"); show("   at longitude", lo20, "deg")
show("IGRF-14 2020 -> 2025 SAA-minimum drift", km_between(la20, lo20, la25, lo25) / 5, "km/yr")
show("IGRF-14 2020 -> 2025 SAA-minimum westward longitude rate", (lo25 - lo20) / 5, "deg/yr")
show("IGRF-14 2020 -> 2025 SAA-minimum intensity rate", (f25 - f20) / 5, "nT/yr")
gw0, hw0 = at(WMM, 2025.0); gw5, hw5 = at(WMM, 2030.0)
fa, laa, loa = local_min(gw0, hw0, -26.0, -60.0); fb, lab, lob = local_min(gw5, hw5, -26.0, -60.0)
show("WMM2025 SV forecast 2025 -> 2030 SAA-minimum drift", km_between(laa, loa, lab, lob) / 5, "km/yr")
show("WMM2025 SV forecast 2025 -> 2030 SAA-minimum intensity rate", (fb - fa) / 5, "nT/yr")
# is there a second (eastern) local minimum south-west of Africa at degree 13?
cand = []
for la in range(-60, -14, 2):
    for lo in range(-30, 41, 2):
        f0 = total_field(G25, H25, la, lo)
        if all(f0 < total_field(G25, H25, la + a, lo + b) for a in (-2, 0, 2) for b in (-2, 0, 2) if (a, b) != (0, 0)):
            cand.append((f0, la, lo))
show("IGRF-14 2025: grid local minima of F between 30W and 40E (count)", len(cand))
for f0, la, lo in cand:
    fe25, lae25, loe25 = local_min(G25, H25, la, lo); fe20, lae20, loe20 = local_min(G20, H20, la, lo)
    show(f"   eastern local minimum 2025 near ({la}, {lo}): F", fe25, "nT")
    show("      at latitude, longitude", lae25, f"deg, {loe25:.3f} deg")
    show("      2020 -> 2025 intensity rate", (fe25 - fe20) / 5, "nT/yr")
    show("      2020 -> 2025 drift", km_between(lae20, loe20, lae25, loe25) / 5, "km/yr")
    show("      eastern/western intensity-rate ratio", (fe25 - fe20) / (f25 - f20))
ESA = (24000 - 22000) / 50
show("ESA Swarm: SAA minimum fell 24000 -> 22000 nT over 1970-2020", -ESA, "nT/yr")
show("Gao et al. 2025: western minimum -290 nT 2015-2024.5", -290 / 9.5, "nT/yr")
show("Gao et al. 2025: eastern minimum -730 nT 2015-2024.5", -730 / 9.5, "nT/yr")

rho_core = 1.1e4                                   # kg/m3, assumed mean outer-core density
for bs in (2e-3, 4e-3):
    va = bs / sqrt(MU0 * rho_core)
    show(f"torsional Alfven speed, B_s = {bs * 1e3:.0f} mT (mhd::alfven_speed)", va * 3.15576e7 / 1e3, "km/yr")
    show(f"   travel time across 1220-3480 km at that speed", (3480 - 1220) / (va * 3.15576e7 / 1e3), "yr")
lam_mc = 1500 * 7
show("MC wave: wavelength = phase speed x period, 1500 km/yr x 7 yr", lam_mc, "km")
show("   azimuthal wavenumber at the CMB, 2 pi 3480 km / lambda", 2 * pi * 3480 / lam_mc)

# ============================================================== 5. high-field magnets (Hahn et al. 2019) and records
print("\n== 5. No-insulation REBCO insert LBC3 (Hahn et al. 2019) and field records ==")
L, RC, KM = 50.4e-3, 47.1e-3, 59.7e-3
check("NI charging time constant L/R_c", L / RC, 1.07, 0.005, "s")
check("turn-to-turn leakage current V/R_c, 9.45 mV", 9.45e-3 / RC, 0.2, 0.005, "A")
check("its centre field, leakage x magnet constant", 9.45e-3 / RC * KM * 1e3, 12, 0.5, "mT")
check("calculated centre field 31.1 T + 59.7 mT/A x 245.3 A", 31.1 + KM * 245.3, 45.7, 0.05, "T")
check("screening-current deficit, calculated minus measured 45.5 T", 31.1 + KM * 245.3 - 45.5, 0.2, 0.05, "T")
check("conductor current density 245.3 A / (4.02 mm x 0.043 mm)", 245.3 / (4.02 * 0.043), 1420, 5, "A/mm2")
check("winding current density 2717 x 245.3 A / (10 mm x 53.1 mm)", 2717 * 245.3 / (10 * 53.1), 1260, 10, "A/mm2")
check("turns: 12 pancakes x 226.4", 12 * 226.4, 2717, 0.5)
s_in, s_out = 45.5 * 1420e6 * 7e-3 / 1e6, 31.1 * 1420e6 * 17e-3 / 1e6
show("hoop stress BJR at the inner radius, 45.5 T, 7 mm", s_in, "MPa"); show("BJR at the outer radius, 31.1 T background, 17 mm", s_out, "MPa")
ident("published peak 691 MPa lies inside the BJR envelope", 0.0 if s_in <= 691 <= s_out else 1.0)
a1, a2, hb = 7e-3, 17e-3, 53.1e-3 / 2
k_thick = MU0 * 2717 / (2 * (a2 - a1)) * log((a2 + sqrt(a2 ** 2 + hb ** 2)) / (a1 + sqrt(a1 ** 2 + hb ** 2)))
show("thick-solenoid centre constant (uniform J) vs published 59.7 mT/A", k_thick * 1e3, "mT/A")
am, cdep = (a1 + a2) / 2, a2 - a1
wheeler = 31.6 * 2717 ** 2 * am ** 2 / (6 * am + 9 * 2 * hb + 10 * cdep) * 1e-6
show("Wheeler multilayer inductance vs published 50.4 mH", wheeler * 1e3, "mH")
show("insert self-energy 0.5 L I^2 at 245.3 A", 0.5 * L * 245.3 ** 2, "J")
show("all-superconducting 35.6 T (2026) over NHMFL 32 T (2017)", 35.6 / 32 - 1, "fraction")
show("resistive 42.02 T over 41.4 T (2017)", 42.02 - 41.4, "T")
show("resistive record: B^2 per MW, 42.02 T at 32.3 MW", 42.02 ** 2 / 32.3, "T^2/MW")

print(f"\n{len(FAILS)} failed: {FAILS}" if FAILS else "\nall checks and identities pass")
