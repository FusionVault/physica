"""Polarizable-vacuum (exponential) metric against GR: every observable computed both ways.

Units G = c = M = 1 unless stated. PV: ds^2 = -e^{-2/r}dt^2 + e^{2/r}(dr^2 + r^2 dΩ^2) (isotropic r).
GR: Schwarzschild, areal R (or isotropic where said). Run: python3 -I pv_metric_tests.py
"""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, curve_fit

E = np.e
PI = np.pi
out = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    out.append(s)


# ---------------------------------------------------------------- 1. series of the metrics
say("== 1. PN expansions (U = M/r isotropic, U_A = M/R areal)")
# PV: -g00 = e^{-2U}; g_ij = e^{2U}.  GR isotropic: -g00 = ((1-U/2)/(1+U/2))^2, g_ij = (1+U/2)^4.
# Numerical Taylor coefficients by polynomial fit on tiny U (checks the hand algebra).
U = np.linspace(-0.02, 0.02, 81)
def coeffs(y, deg=6):
    return np.polynomial.polynomial.polyfit(U, y, deg)
for name, g00, gij in [("PV", np.exp(-2 * U), np.exp(2 * U)),
                       ("GR", ((1 - U / 2) / (1 + U / 2)) ** 2, (1 + U / 2) ** 4)]:
    c0, c1 = coeffs(g00), coeffs(gij)
    say(f"{name}: -g00 = 1 {c0[1]:+.6f}U {c0[2]:+.6f}U^2 {c0[3]:+.6f}U^3;  g_ij = 1 {c1[1]:+.6f}U {c1[2]:+.6f}U^2 {c1[3]:+.6f}U^3")
# PPN-2: -g00 = 1 - 2U + 2βU^2 + ...;  g_ij = 1 + 2γU + (3/2)δ U^2  ->  δ_PV = 4/3, δ_GR = 1
say("2PN parameter δ (g_ij = 1+2γU+(3/2)δU^2): PV 4/3, GR 1;  β = γ = 1 both")

# areal-coordinate tt series of PV (Psaltis/EHT κ convention: g_tt = -1 + 2/R - 2κ1/R^2 + 2κ2/R^3 - 2κ3/R^4)
def pv_u_of_UA(UA):  # u = M/r solves u e^{-u} = U_A, outer branch (Lambert W0)
    return brentq(lambda u: u * np.exp(-u) - UA, 0.0, 1.0, xtol=1e-17, rtol=8.9e-16)
UAs = np.linspace(1e-4, 0.03, 60)
mgtt = np.array([np.exp(-2 * pv_u_of_UA(x)) for x in UAs])
c = np.polynomial.polynomial.polyfit(UAs, mgtt - (1 - 2 * UAs), 7)
say(f"PV areal: -g_tt - (1-2U_A) ≈ {c[2]:+.2e}U_A^2 {c[3]:+.6f}U_A^3 {c[4]:+.6f}U_A^4 {c[5]:+.5f}U_A^5")
say(f"  -> κ1 = 0, κ2 = {-c[3]/2:.6f} (exact 1/6), κ3 = {c[4]/2:.6f}")
zeta = 1 / 6
say(f"Psaltis 2020 truncated-2PN shadow 3√3(1+ζ/9) with ζ=κ2=1/6: {np.sqrt(27)*(1+zeta/9):.4f} "
    f"({zeta/9*100:.2f}% larger)  vs exact 2e = {2*E:.4f} ({(2*E/np.sqrt(27)-1)*100:.3f}% larger)")

# ---------------------------------------------------------------- 2. exact deflection
# α = 2∫_0^{u0} du/√F − π, u = M/r.  PV (isotropic): F = e^{4u}/b^2 − u^2.  GR (areal): F = 1/b^2 − u^2(1−2u).
def turning(b, metric):
    g = (lambda u: u * u * np.exp(-4 * u)) if metric == "PV" else (lambda u: u * u * (1 - 2 * u))
    photon = 0.5 if metric == "PV" else 1 / 3
    return brentq(lambda u: b * b * g(u) - 1.0, 0.0, photon, xtol=1e-300, rtol=8.9e-16, maxiter=500)

GL_X, GL_W = np.polynomial.legendre.leggauss(16)
def deflection(b, metric):
    u0 = turning(b, metric)
    def q(u):
        if metric == "GR":
            return (u + u0) - 2 * (u * u + u * u0 + u0 * u0)
        d = u - u0
        slope = np.where(d == 0, 4.0, np.expm1(4 * d) / np.where(d == 0, 1, d))
        return (u + u0) - np.exp(4 * u0) * slope / (b * b)
    s, upper = 0.0, 1.0
    for k in range(1, 50):  # dyadic panels toward t = 0 (the near-photon-sphere peak)
        lower = 0.0 if k == 49 else 0.5 ** k
        mid, half = 0.5 * (upper + lower), 0.5 * (upper - lower)
        t = mid + half * GL_X
        s += half * np.sum(GL_W * 2 * np.sqrt(u0) / np.sqrt(q(u0 * (1 - t * t))))
        upper = lower
    return 2 * s - PI

say("\n== 2. Light bending")
# check against the series in the repo docstring: GR = 4x + 15π/4 x^2 + 128/3 x^3 (+ 3465π/64 x^4)
xs = np.linspace(2e-4, 8e-3, 40)
for metric in ("GR", "PV"):
    al = np.array([deflection(1 / x, metric) for x in xs])
    cfit = np.polynomial.polynomial.polyfit(xs, al, 7)
    say(f"{metric}: α = {cfit[1]:.8f}x + {cfit[2]:.6f}x^2 + {cfit[3]:.4f}x^3 + {cfit[4]:.2f}x^4  "
        f"(c2/π = {cfit[2]/PI:.6f})")
say(f"GR exact series values: 15π/4 = {15*PI/4:.6f}, 128/3 = {128/3:.4f}, 3465π/64 = {3465*PI/64:.2f}; PV 4π = {4*PI:.6f}")
# medium picture: n = 1 + 2u + n2 u^2 ->  α2 = π(n2 + 2): PV n2 = 2 -> 4π ; GR isotropic n2 = 7/4 -> 15π/4
m_sun, R_sun = 1476.625, 6.957e8   # GM/c^2 (m), solar radius (m)
x = m_sun / R_sun
uas = PI / 180 / 3600 * 1e-6
say(f"Sun limb: 1PN 4x = {4*x/uas/1e6:.6f} arcsec; GR 2PN = {15*PI/4*x*x/uas:.3f} µas; PV 2PN = {4*PI*x*x/uas:.3f} µas; "
    f"PV−GR = {(PI/4)*x*x/uas:.3f} µas (= 1/15 = {1/15*100:.2f}% of the 2PN term)")
say(f"LATOR (proposed): 2PN term to ~2e-3 -> PV/GR separation {(1/15)/2e-3:.0f}σ if that accuracy is reached")

say("\n== 3. Strong deflection (Bozza 2002 formulas)")
def bozza(A, Ap, App, B, C, Cp, Cpp, Ainv, xm):
    ym, Cm = A(xm), C(xm)
    um = np.sqrt(Cm / ym)
    beta_m = Cm * (1 - ym) ** 2 * (Cpp(xm) * ym - Cm * App(xm)) / (2 * ym ** 2 * Cp(xm) ** 2)
    R = lambda z: 2 * np.sqrt(B(Ainv((1 - ym) * z + ym)) * ((1 - ym) * z + ym)) / (C(Ainv((1 - ym) * z + ym)) * Ap(Ainv((1 - ym) * z + ym))) * (1 - ym) * np.sqrt(Cm)
    f = lambda z: 1 / np.sqrt(ym - ((1 - ym) * z + ym) * Cm / C(Ainv((1 - ym) * z + ym)))
    R0 = R(0.0)
    g = lambda z: R(z) * f(z) - R0 / (np.sqrt(beta_m) * z)
    bR = quad(g, 1e-12, 1, limit=400, epsabs=1e-13, epsrel=1e-13)[0]
    abar = R0 / (2 * np.sqrt(beta_m))
    bbar = -PI + bR + abar * np.log(2 * beta_m / ym)
    return abar, bbar, um, bR, beta_m

# Schwarzschild, areal x, M = 1
gr = bozza(lambda x: 1 - 2 / x, lambda x: 2 / x ** 2, lambda x: -4 / x ** 3, lambda x: 1 / (1 - 2 / x),
           lambda x: x * x, lambda x: 2 * x, lambda x: 2.0, lambda y: 2 / (1 - y), 3.0)
pv = bozza(lambda x: np.exp(-2 / x), lambda x: 2 / x ** 2 * np.exp(-2 / x), lambda x: (4 / x ** 4 - 4 / x ** 3) * np.exp(-2 / x),
           lambda x: np.exp(2 / x), lambda x: x * x * np.exp(2 / x), lambda x: (2 * x - 2) * np.exp(2 / x),
           lambda x: (2 - 4 / x + 4 / x ** 2) * np.exp(2 / x), lambda y: -2 / np.log(y), 2.0)
say(f"GR: ā = {gr[0]:.10f}, b̄ = {gr[1]:.10f} (Bozza: 1, −0.4002; exact ln[216(7−4√3)]−π = {np.log(216*(7-4*np.sqrt(3)))-PI:.10f}), u_m = {gr[2]:.6f}, b_R = {gr[3]:.6f} (Bozza 0.9496)")
say(f"PV: ā = {pv[0]:.10f}, b̄ = {pv[1]:.10f}, u_m = {pv[2]:.6f} (2e = {2*E:.6f}), b_R = {pv[3]:.8f}, β_m = {pv[4]:.6f} (e(1−1/e)^2 = {E*(1-1/E)**2:.6f})")
say(f"   PV b̄ closed part: −π + b_R + ln[2(e−1)^2] with ln[2(e−1)^2] = {np.log(2*(E-1)**2):.8f}")
# direct check: α(b) + ā ln(b/b_c − 1) -> b̄ as b -> b_c
for metric, bc, ab in (("GR", np.sqrt(27), gr), ("PV", 2 * E, pv)):
    eps = np.array([1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6, 1e-6])
    res = np.array([deflection(bc * (1 + e), metric) + np.log(e) for e in eps])
    # residual ~ b̄ + c ε ln ε + d ε: least squares
    Afit = np.vstack([np.ones_like(eps), eps * np.log(eps), eps]).T
    sol = np.linalg.lstsq(Afit, res, rcond=None)[0]
    say(f"{metric}: α+ln(ε) at ε=1e-6: {res[-1]:.8f}; extrapolated b̄ = {sol[0]:.8f} (Bozza formula {ab[1]:.8f})")
# images / photon-ring subrings: impact parameter where α = nπ ; ratio of successive offsets -> e^{π/ā}
def image(angle, metric, bc):
    lo, hi = bc * (1 + 1e-14), 2 * bc
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if deflection(mid, metric) > angle:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)
for metric, bc in (("GR", np.sqrt(27)), ("PV", 2 * E)):
    b4, b5 = image(4 * PI, metric, bc), image(5 * PI, metric, bc)
    say(f"{metric}: (b_4π−b_c)/(b_5π−b_c) = {(b4-bc)/(b5-bc):.6f}  (e^π = {np.exp(PI):.6f});  "
        f"b(α=2π)/b_c − 1 = {image(2*PI, metric, bc)/bc-1:.6e}; b(α=3π)/b_c − 1 = {image(3*PI, metric, bc)/bc-1:.6e}")
# Bozza observables for Sgr A* (VLTI mass/distance used by EHT): θ∞, s, r (mag)
GMsun_c2 = 1476.625
pc = 3.0856775814913673e16
for lab, Msol, Dpc in (("VLTI", 4.297e6, 8277.0), ("Keck", 3.951e6, 7935.0)):
    thg = Msol * GMsun_c2 / (Dpc * pc) / uas
    for metric, (ab, bb, um) in (("GR", gr[:3]), ("PV", pv[:3])):
        thinf = um * thg
        s = thinf * np.exp((bb - 2 * PI) / ab)
        rmag = 2.5 * np.log10(np.exp(2 * PI / ab))
        say(f"Sgr A* {lab} θ_g={thg:.4f} µas {metric}: shadow diameter {2*thinf:.3f} µas, θ∞={thinf:.3f} µas, s={s*1e3:.3f} nas, r={rmag:.4f} mag")

# ---------------------------------------------------------------- 4. EHT tension
say("\n== 4. Shadow size vs EHT")
d_pv = 2 * E / np.sqrt(27) - 1
say(f"δ_PV = 2e/√27 − 1 = {d_pv:.6f}; Kerr range −0.075..0 (EHT); MT wormhole Φ=−M/R (EHT's example): r_sh = e = {E:.4f}, δ = {E/np.sqrt(27)-1:.4f}")
rows = [("eht-imaging VLTI GRMHD (fiducial)", -0.08, 0.09, 0.09), ("eht-imaging Keck GRMHD (fiducial)", -0.04, 0.09, 0.10),
        ("eht-imaging VLTI an.Kerr", -0.05, 0.09, 0.11), ("eht-imaging Keck an.Kerr", -0.02, 0.11, 0.11),
        ("eht-imaging VLTI an.non-Kerr", -0.07, 0.10, 0.09), ("eht-imaging Keck an.non-Kerr", -0.02, 0.10, 0.09),
        ("SMILI VLTI GRMHD", -0.10, 0.12, 0.10), ("SMILI Keck GRMHD", -0.06, 0.13, 0.10),
        ("DIFMAP VLTI GRMHD", -0.12, 0.10, 0.08), ("DIFMAP Keck GRMHD", -0.08, 0.09, 0.09),
        ("mG-ring VLTI GRMHD", -0.17, 0.11, 0.10), ("mG-ring Keck GRMHD", -0.13, 0.11, 0.11),
        ("M87* 2017 stellar-dynamics prior", -0.01, 0.17, 0.17)]
for lab, mu, up, lo in rows:
    z_pv = (d_pv - mu) / up
    z_gr = (0 - mu) / up
    say(f"{lab:38s} δ={mu:+.2f}+{up:.2f}/−{lo:.2f}: PV at {z_pv:.2f}σ, Schwarzschild at {z_gr:.2f}σ, "
        f"Gaussian likelihood ratio GR:PV = {np.exp(0.5*(z_pv**2-z_gr**2)):.2f}")
say(f"Separation PV−Schwarzschild in units of the fiducial σ=0.09: {d_pv/0.09:.2f}σ; σ needed for 3σ: {d_pv/3:.4f}; for 5σ: {d_pv/5:.4f}")

# ---------------------------------------------------------------- 5. circular orbits, ISCO, efficiency
say("\n== 5. Circular orbits")
def circ(Af, Apf, Cf, Cpf, r):
    Om = np.sqrt(Apf(r) / Cpf(r))
    ut = 1 / np.sqrt(Af(r) - Cf(r) * Om ** 2)
    return Om, Af(r) * ut, Cf(r) * Om * ut   # Ω, E, L per unit mass
PVf = (lambda r: np.exp(-2 / r), lambda r: 2 / r ** 2 * np.exp(-2 / r), lambda r: r * r * np.exp(2 / r), lambda r: (2 * r - 2) * np.exp(2 / r))
GRf = (lambda r: 1 - 2 / r, lambda r: 2 / r ** 2, lambda r: r * r, lambda r: 2 * r)
r_pv, r_gr = 3 + np.sqrt(5), 6.0
for lab, f, r, areal in (("PV", PVf, r_pv, r_pv * np.exp(1 / r_pv)), ("GR", GRf, r_gr, r_gr)):
    Om, En, L = circ(*f, r)
    say(f"{lab} ISCO r={r:.6f} (areal {areal:.6f}): MΩ = {Om:.7f}, E = {En:.7f}, efficiency 1−E = {(1-En)*100:.4f}%, L = {L:.6f}")
Om_pv, Om_gr = circ(*PVf, r_pv)[0], circ(*GRf, r_gr)[0]
say(f"Ω_ISCO PV/GR − 1 = {(Om_pv/Om_gr-1)*100:.3f}%  (Makukov & Mychelkin 2018 quote 0.0633326 and 6.92%)")
say(f"f_GW(ISCO)=Ω/π for 10 M_sun: GR {Om_gr/PI/(10*4.925490947e-6):.1f} Hz, PV {Om_pv/PI/(10*4.925490947e-6):.1f} Hz")
# max redshift of light from a static emitter at the photon sphere / ISCO (static clocks)
for lab, A, r in (("PV ISCO static", PVf[0], r_pv), ("GR ISCO static", GRf[0], r_gr), ("PV photon sphere static", PVf[0], 2.0), ("GR photon sphere static", GRf[0], 3.0)):
    say(f"{lab}: 1+z = {1/np.sqrt(A(r)):.5f}")
# near-circular periastron advance K = Ω_φ/Ω_r vs x = (MΩ)^{2/3}  (gauge invariant)
SECOND = {  # (A'', C'', g_rr) analytic
    "PV": (lambda r: (4 / r ** 4 - 4 / r ** 3) * np.exp(-2 / r), lambda r: (2 - 4 / r + 4 / r ** 2) * np.exp(2 / r), lambda r: np.exp(2 / r)),
    "GR": (lambda r: -4 / r ** 3, lambda r: 2.0 + 0 * r, lambda r: 1 / (1 - 2 / r)),
}
def K_of_r(lab, r):
    f = PVf if lab == "PV" else GRf
    Af, Apf, Cf, Cpf = f
    App, Cpp, grr = SECOND[lab]
    Om, En, L = circ(*f, r)
    # g_rr ṙ^2 = W = E^2/A − L^2/C − 1; at a circular orbit W = W' = 0 and ω_r^2 = −W''/(2 g_rr) per proper time
    W2 = En ** 2 * (2 * Apf(r) ** 2 / Af(r) ** 3 - App(r) / Af(r) ** 2) - L ** 2 * (2 * Cpf(r) ** 2 / Cf(r) ** 3 - Cpp(r) / Cf(r) ** 2)
    Om_r = np.sqrt(-W2 / (2 * grr(r))) * Af(r) / En
    return Om / Om_r, Om ** (2 / 3)
for lab in ("GR", "PV"):
    rs = np.geomspace(80, 4e3, 40)
    Ks, xs_ = np.array([K_of_r(lab, r) for r in rs]).T
    cK = np.polynomial.polynomial.polyfit(xs_, Ks - 1 - 3 * xs_, [2, 3, 4, 5])
    say(f"{lab}: K = Ω_φ/Ω_r = 1 + 3x + {cK[2]:.4f}x^2 + {cK[3]:.2f}x^3 ...   (GR exact (1−6x)^(−1/2): 27/2, 135/2)")
for lab, rr in (("GR", [7.0, 10.0, 20.0]), ("PV", [6.0, 9.0, 19.0])):
    say(f"{lab}: " + ", ".join(f"x={K_of_r(lab, r)[1]:.5f}: K={K_of_r(lab, r)[0]:.5f}" for r in rr)
        + (f"; GR check (1−6x)^−1/2 at x=1/10: {(1-0.6)**-0.5:.5f}" if lab == "GR" else ""))
# null-orbit Lyapunov exponent λ (per unit t) and γ = λ π b_c (subring demagnification exponent per half orbit)
for lab, rm, bc in (("GR", 3.0, np.sqrt(27)), ("PV", 2.0, 2 * E)):
    f = PVf if lab == "PV" else GRf
    Af, Apf, Cf, Cpf = f
    App, Cpp, grr = SECOND[lab]
    En, L = 1.0, bc
    W2 = En ** 2 * (2 * Apf(rm) ** 2 / Af(rm) ** 3 - App(rm) / Af(rm) ** 2) - L ** 2 * (2 * Cpf(rm) ** 2 / Cf(rm) ** 3 - Cpp(rm) / Cf(rm) ** 2)
    lam = np.sqrt(W2 / (2 * grr(rm))) * Af(rm) / En
    say(f"{lab}: photon orbit λ = {lam:.6f}/M, Ω_c = {1/bc:.6f}/M, γ = λπ/Ω_c = {lam*PI*bc:.6f} (π = {PI:.6f}); "
        f"eikonal l=2,n=0: Mω ≈ {2.5/bc:.4f} − {0.5*lam:.4f}i")

# ---------------------------------------------------------------- 6. S2 precession at fixed E, L
say("\n== 6. S2 periapsis advance per orbit, same E and L in both metrics")
def precession(f, grr, En, L):
    Af, Apf, Cf, Cpf = f
    W = lambda r: En ** 2 / Af(r) - L ** 2 / Cf(r) - 1
    # turning points: bracket around Newtonian guesses
    rs = np.geomspace(50, 1e6, 20000)
    w = W(rs)
    idx = np.where(np.diff(np.sign(w)))[0]
    r1 = brentq(W, rs[idx[0]], rs[idx[0] + 1], xtol=1e-12, rtol=8.9e-16)
    r2 = brentq(W, rs[idx[-1]], rs[idx[-1] + 1], xtol=1e-12, rtol=8.9e-16)
    a, d = 0.5 * (r1 + r2), 0.5 * (r2 - r1)
    integrand = lambda ch: (L / Cf(a - d * np.cos(ch))) * np.sqrt(grr(a - d * np.cos(ch))) * d * np.sin(ch) / np.sqrt(np.maximum(W(a - d * np.cos(ch)), 1e-300))
    val = quad(integrand, 0, PI, limit=500, epsabs=1e-14, epsrel=1e-14)[0]
    return 2 * val - 2 * PI
# S2 (GRAVITY 2020 numbers quoted in EHT Sgr A* VI): a = 125 mas at 8246.7 pc, e = 0.8846 (eccentricity from GRAVITY 2020 not opened: use 0.88466)
M_S2 = 4.261e6 * GMsun_c2
a_m = 0.125 / 206264.806 * 8246.7 * pc
e_S2 = 0.88466
p = a_m * (1 - e_S2 ** 2) / M_S2   # semi-latus rectum in M
# Newtonian E, L from p, e; then identical (E, L) in both metrics
L = np.sqrt(p) * 1.0
En = np.sqrt(1 - (1 - e_S2 ** 2) / p)
dgr = precession(GRf, lambda r: 1 / (1 - 2 / r), En, L)
dpv = precession(PVf, lambda r: np.exp(2 / r), En, L)
say(f"S2: p = {p:.1f} M; Δφ_GR = {dgr:.6e} rad = {dgr*206264.8/60:.3f} arcmin; Δφ_PV = {dpv:.6e}; PV−GR = {dpv-dgr:.3e} rad "
    f"= {(dpv-dgr)/dgr:.2e} of the total; 1PN 6π/p = {6*PI/p:.6e}")
say(f"   measured f_SP = 1.10 ± 0.19 -> σ(Δφ) ≈ {0.19*6*PI/p:.2e} rad; PV−GR is {(dpv-dgr)/(0.19*6*PI/p):.1e} σ")

# ---------------------------------------------------------------- 7. Shapiro delay at 2PN
say("\n== 7. Shapiro delay at 2PN")
# Optical radius w = n r (= b at the turning point): T = ∫ G(w) w dw/√(w²−b²), G = n dr/dw = 1/(1 + r (ln n)').
# G−1: PV (2/r)/(1−2/r); GR isotropic: −s/(1+s), s = −x[3/(1+x)+1/(1−x)], x = 1/2r.
def r_of_w(w, metric):
    if metric == "PV":
        return brentq(lambda r: r * np.exp(2 / r) - w, 2.0, w, xtol=1e-300, rtol=8.9e-16)
    n = lambda r: (1 + 0.5 / r) ** 3 / (1 - 0.5 / r)
    return brentq(lambda r: n(r) * r - w, (2 + np.sqrt(3)) / 2, w, xtol=1e-300, rtol=8.9e-16)
def Gm1(w, metric):
    r = r_of_w(w, metric)
    if metric == "PV":
        return (2 / r) / (1 - 2 / r)
    x = 0.5 / r
    s = -x * (3 / (1 + x) + 1 / (1 - x))
    return -s / (1 + s)
def dT_fixed_b(b):  # full pass, endpoints at infinity, same b
    f = lambda eta: (Gm1(b * np.cosh(eta), "PV") - Gm1(b * np.cosh(eta), "GR")) * b * np.cosh(eta)
    return 2 * quad(f, 0, 40, limit=400, epsabs=0, epsrel=1e-10)[0]
for bR in (1.0, 2.0, 10.0):
    b = bR * R_sun / m_sun
    dTb = dT_fixed_b(b)
    dal = deflection(b, "PV") - deflection(b, "GR")
    dT_ends = dTb - b * dal  # fixed emitter and receiver (Legendre/envelope: δT|Θ = δT|b − b δΘ|b)
    to_ps = m_sun / 299792458 * 1e12
    say(f"Sun b = {bR:4.1f} R_sun: ΔT|b = {dTb*to_ps:.3f} ps (π/2 M²/b = {PI/2/b*to_ps:.3f}); Δα = {dal:.4e} (π/4 M²/b² = {PI/4/b**2:.4e}); "
        f"fixed endpoints ΔT = {dT_ends*to_ps:.3f} ps (π/4 M²/b = {PI/4/b*to_ps:.3f} ps)")
b16 = 1.6 * R_sun / m_sun
say(f"Cassini geometry b = 1.6 R_sun: fixed-endpoint ΔT = {PI/4/b16*m_sun/299792458*1e12:.3f} ps one way; Doppler-equivalent δγ = (π/8)M/b = {PI/8/b16:.2e}; "
    f"deflection-equivalent δγ at the limb (π/8)M/R_sun = {PI/8*m_sun/R_sun:.2e}")
say("   GR Shapiro one-way at the limb ≈ 2M(1+γ)ln(4 r_E r_p/b²)/c: for r_E = 1 AU, r_p = 10 AU, "
    f"{2*m_sun*2*np.log(4*1.496e11*1.496e12/R_sun**2)/299792458*1e6:.1f} µs")

# ---------------------------------------------------------------- 8. surface redshift at given areal radius
say("\n== 8. Static redshift at the same areal radius R (GM/c^2R = U_A)")
for UA in (1 / 1400.0, 0.1, 0.172, 0.25, 0.3):
    u = pv_u_of_UA(UA)
    z_pv = np.exp(u) - 1
    z_gr = 1 / np.sqrt(1 - 2 * UA) - 1
    say(f"U_A = {UA:.5f}: z_GR = {z_gr:.6e}, z_PV = {z_pv:.6e}, PV−GR = {z_pv-z_gr:.3e} (U_A^3/6 = {UA**3/6:.3e})")

# ---------------------------------------------------------------- 9. Binary-pulsar radiation (Ibison 2003 eq 65)
say("\n== 9. Binary pulsar: PV/GR radiated power, Ibison (astro-ph/0302273) eq 65")
ratio = lambda e: 0.5 * (51 * e ** 4 + 396 * e ** 2 + 128) / (37 * e ** 4 + 292 * e ** 2 + 96)
for lab, e, meas, sig in (("circular", 0.0, None, None), ("PSR B1913+16", 0.6171340, 0.9983, 0.0016), ("PSR J0737-3039A/B", 0.0878, 1.0, 1.3e-4 / 1.96)):
    r_ = ratio(e)
    extra = f"; measured {meas} ± {sig:.1e} -> PV off by {(meas-r_)/sig:.0f}σ" if meas else ""
    say(f"{lab}: e = {e}: P_PV/P_GR = {r_:.5f}{extra}")
# check of eq 64 normalisation: f(1, e) is Peters–Mathews
pm = lambda e: (1 + 73 / 24 * e ** 2 + 37 / 96 * e ** 4) / (1 - e ** 2) ** 3.5
fl = lambda lam, e: ((227 + 69 * lam) * e ** 4 + 4 * (443 + 141 * lam) * e ** 2 + 192 * (3 + lam)) / (768 * (1 - e ** 2) ** 3.5)
say(f"eq 64 check: f(1,0.6171) = {fl(1,0.6171):.6f}, Peters–Mathews = {pm(0.6171):.6f}; f(−1/3,e)/f(1,e) at e=0.6171: {fl(-1/3,0.6171)/fl(1,0.6171):.5f}")

# ---------------------------------------------------------------- 10. QNMs of test fields (Boonserm et al. 2018 potentials)
say("\n== 10. Test-field ringing (fixed metric, 'inverse Cowling'); time-domain + fit")
def tortoise_grid(metric, rstar):
    # returns r(r*) on a uniform r* grid
    if metric == "GR":
        # r* = r + 2 ln(r/2 − 1); solve for L = ln(r − 2) (monotonic), return r − 2 to keep 1 − 2/r exact
        rho = np.empty_like(rstar)
        for i, s in enumerate(rstar):
            Lr = brentq(lambda L: 2 + np.exp(L) + 2 * (L - np.log(2)) - s, -2000, 20)
            rho[i] = np.exp(Lr)
        return rho
    # PV: r* = ∫_2^r e^{2/s} ds  (r* = 0 at the photon sphere)
    rg = np.concatenate([np.geomspace(0.12, 2.0, 40000)[:-1], np.geomspace(2.0, 4000, 40000)])
    w = np.exp(2 / rg)
    cum = np.concatenate([[0], np.cumsum(0.5 * (w[1:] + w[:-1]) * np.diff(rg))])
    cum -= np.interp(2.0, rg, cum)
    return np.interp(rstar, cum, rg)
def potential(metric, s, l, r):
    if metric == "GR":  # here r holds ρ = r − 2
        rho, r = r, r + 2
        return (rho / r) * (l * (l + 1) / r ** 2 + (1 - s * s) * 2 / r ** 3)
    return np.exp(-4 / r) * (l * (l + 1) + (s * s - 1) * ((1 - 1 / r) ** 2 - 1)) / r ** 2
def ring(metric, s, l, L=700.0, h=0.02, T=260.0, obs=30.0):
    xs_ = np.arange(-L, L + h / 2, h)
    r = tortoise_grid(metric, xs_)
    V = potential(metric, s, l, r)
    V[~np.isfinite(V)] = 0.0
    psi = np.exp(-((xs_ - 10.0) ** 2) / 8.0)
    prev = psi.copy()          # zero initial velocity
    dt = 0.5 * h
    iobs = np.argmin(abs(xs_ - obs))
    ts, sig = [], []
    for n in range(int(T / dt)):
        lap = np.zeros_like(psi)
        lap[1:-1] = (psi[2:] - 2 * psi[1:-1] + psi[:-2]) / h ** 2
        new = 2 * psi - prev + dt * dt * (lap - V * psi)
        prev, psi = psi, new
        ts.append((n + 1) * dt)
        sig.append(psi[iobs])
    return np.array(ts), np.array(sig), V, r, xs_
def fit_mode(t, y, t0, t1):
    m = (t > t0) & (t < t1)
    model = lambda t, A, wi, wr, ph: A * np.exp(-wi * (t - t0)) * np.cos(wr * (t - t0) + ph)
    best = None
    for wr0 in np.linspace(0.08, 0.8, 19):
        try:
            p, _ = curve_fit(model, t[m], y[m], p0=[np.max(abs(y[m])), 0.09, wr0, 0.0], maxfev=20000)
            res = np.sum((model(t[m], *p) - y[m]) ** 2)
            if p[2] > 0 and (best is None or res < best[1]):
                best = (p, res)
        except Exception:
            pass
    return best[0][2], best[0][1]
for metric, s, l in (("GR", 2, 2), ("PV", 2, 2), ("GR", 0, 2), ("PV", 0, 2), ("GR", 1, 1), ("PV", 1, 1)):
    t, y, V, r, xs_ = ring(metric, s, l)
    pk = np.argmax(np.where(np.isfinite(V), V, 0))
    # count local maxima of V (double barrier -> trapped modes / echoes)
    inner = V[1:-1]
    nmax = int(np.sum((inner > V[:-2]) & (inner > V[2:]) & (inner > 1e-6)))
    wr, wi = fit_mode(t, y, 70.0, 140.0)
    wr2, wi2 = fit_mode(t, y, 55.0, 110.0)
    rp = r[pk] + 2 if metric == "GR" else r[pk]
    say(f"{metric} s={s} l={l}: Mω ≈ {wr:.4f} − {wi:.4f}i (window 55–110: {wr2:.4f} − {wi2:.4f}i); V peak {V[pk]:.5f} at r = {rp:.4f}; local maxima of V: {nmax}")

open(__file__.replace(".py", "_output.txt"), "w").write("\n".join(out) + "\n")
