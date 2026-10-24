"""Reproduce numbers quoted in the gravity DIRDs and the follow-up literature.
Pure Python (no numpy). CODATA-ish constants; DIRD values used where the DIRD states them."""
import math

G, c, hbar, h = 6.67430e-11, 2.99792458e8, 1.054571817e-34, 6.62607015e-34
lP = math.sqrt(hbar * G / c**3)          # Planck length
mP = math.sqrt(hbar * c / G)             # Planck mass
Msun, MJ, ME, RE, Rsun, AU = 1.98847e30, 1.89813e27, 5.9722e24, 6.3781e6, 6.957e8, 1.495978707e11
LY = 9.4607e15
pr = lambda k, v: print(f"{k:58s} {v}")

print("== DIRD 15 Puthoff: harmonic-coordinate Schwarzschild (eq 4) == areal Schwarzschild")
m_geo = G * Msun / c**2
for r in (10 * m_geo, 3 * m_geo):
    R = r + m_geo                                     # areal radius
    g00_dird = (1 - m_geo / r) / (1 + m_geo / r)
    g00_areal = 1 - 2 * m_geo / R
    assert abs(g00_dird - g00_areal) < 1e-12
pr("eq(4) g00 == 1-2GM/(Rc^2) with R=r+GM/c^2", "OK")

print("== DIRD 13 Obousy: eq 2.3, Table 2, eq 2.4 ==")
G13, c13 = 6.673e-11, 3e8
k = c13**4 / G13
pr("c^4/G (DIRD: 1.21e44)", f"{k:.3e}")
R, sigma = 50.0, 1e3
for v in (1e-5, 1e-4, 0.01, 0.5, 1, 2, 10, 100):
    pr(f"  E_warp(v={v}) = -(c^4/G) v^2 R^2 sigma [J]", f"{-k*v*v*R*R*sigma:.3e}")
Mship, Delta = 1e6, 1e-3
vmax = math.sqrt(G13 / c13**2 * Mship * Delta / R**2)
pr("G/c^2 (DIRD: 7.41e-28)", f"{G13/c13**2:.3e}")
pr("eq 2.4 v_max (DIRD: 1.72e-14)", f"{vmax:.3e}  -> {vmax*c13:.2e} m/s")
Rv1 = math.sqrt(G13 / c13**2 * Mship * Delta)
pr("R that gives v_max = 1 from eq 2.4 (DIRD says 3.37e24 m)", f"{Rv1:.3e} m")
pr("  the DIRD figure equals R^2/(G M/c^2) =", f"{R**2/(G13/c13**2*Mship):.3e}  (dimensionally wrong)")
H0 = 70e3 / 3.0857e22
pr("H0 [1/s]", f"{H0:.3e}")
pr("Hubble rate for v=c across 1 m: c/1m / H0 (DIRD: 1e26)", f"{c/1.0/H0:.2e}")
rhoL = 0.69 * 3 * H0**2 * c**2 / (8 * math.pi * G)
pr("dark-energy density rho_L c^2 [J/m^3]", f"{rhoL:.2e}")
pr("x 1e52 (DIRD eq 7.5: 1e42 J/m^3)", f"{rhoL*1e52:.2e}")
pr("Table 3 at v=1 if 100 m^3 x 1e42 J/m^3 (DIRD: 1e42 J)", f"{100*1e42:.1e}")

print("== DIRD 18 Davis wormholes ==")
G18, c18, hb18 = 6.673e-11, 3e8, 1.055e-34
pref = (math.pi**3 / 30) ** 0.25 * math.sqrt(math.sqrt(hb18 * G18 / c18**3))
pr("eq 8 prefactor (DIRD: 4.05e-18 m^1/2)", f"{pref:.3e}")
for r in (1e3, 1.50e11):
    d = pref * math.sqrt(r)
    rho = -math.pi**2 * hb18 * c18 / (720 * d**4)
    pr(f"  r_throat={r:.2e} m: d, rho_Casimir", f"{d:.3e} m, {rho:.3e} J/m^3")
pr("eq 9 c^2/G (DIRD: 1.35e27 kg per m)", f"{c18**2/G18:.3e}")
for lam in (1000, 100, 10, 1, 0.1, 0.01):
    M = lam * c18**2 / G18
    pr(f"  lambda={lam} m: M_wh/MJ, M_wh/ME", f"{M/1.90e27:.3g} MJ, {M/5.98e24:.3g} ME")
# Morris-Thorne throat tension tau = c^4/(8 pi G b0^2)
for b0 in (10.0, 1.0):
    tau = c**4 / (8 * math.pi * G * b0**2)
    pr(f"  MT throat tension b0={b0} m [Pa] (MT: 5e41 dyn/cm^2 at 10 m)", f"{tau:.3e} Pa = {tau*10:.2e} dyn/cm^2")
# thin-shell spherical throat a: sigma = -c^4/(2 pi G a)
a = 1.0
pr("thin-shell sphere a=1 m: sigma [J/m^2], total E [J], E/c^2 [kg]",
   f"{-c**4/(2*math.pi*G*a):.3e}, {-2*c**4*a/G:.3e}, {-2*c**2*a/G:.3e}")
# gravitational squeezing Table 1 / eq 7
rsE = 2 * G * ME / c**2
lam = 8 * math.pi * rsE
pr("Earth r_s, lambda=8 pi r_s (DIRD: 8.87e-3 m, 0.23 m)", f"{rsE:.3e}, {lam:.3f}")
pr("rho=-2 pi^2 hbar c/lambda^4 (DIRD -2.23e-22)", f"{-2*math.pi**2*hbar*c/lam**4:.3e}")
pr("  same with h instead of hbar", f"{-2*math.pi**2*h*c/lam**4:.3e}")
lc = math.sqrt(RE**3 * c**2 / (G * ME))
pr("Ford lambda_c at Earth surface (DIRD 2.42e11 m)", f"{lc:.3e}")
pr("  rho=-2 pi^2 hbar c/lambda_c^4 (DIRD -1.82e-70)", f"{-2*math.pi**2*hbar*c/lc**4:.3e}")

print("== DIRD 19 antigravity ==")
pr("GM_E/R_E [MJ/kg] (DIRD 62.5)", f"{6.673e-11*5.972e24/6.378e6/1e6:.2f}")
pr("eta0=16 pi G/c^2 (DIRD 3.73e-26 m/kg)", f"{16*math.pi*6.673e-11/9e16:.3e}")
pr("gamma0=1/(4 pi G) (DIRD 1.19e9)", f"{1/(4*math.pi*6.673e-11):.3e}")
pr("G rho/c^2 for rho=5500 (DIRD 4.08e-24 m^-2)", f"{6.673e-11*5500/9e16:.3e}")

print("== DIRD 21 HFGW (Baker) ==")
pr("32G/(5c^5) (Baker jerk-formula constant 1.76e-52)", f"{32*G/(5*c**5):.3e}")
P = 1.76e-52 * (2 * 10 * 4e8 / 2e-10) ** 2
pr("P for r=10 m, df=4e8 N, dt=2e-10 s (DIRD 2.8e-13 W)", f"{P:.2e}")
for name, conv in (("F=c^3 w^2 A^2/(8 pi G)", 8), ("F=c^3 w^2 A^2/(16 pi G)", 16), ("F=c^3 w^2 A^2/(32 pi G)", 32)):
    kA = math.sqrt(conv * math.pi * G / c**3) / (2 * math.pi)
    pr(f"  A = k sqrt(F)/nu, k from {name}", f"{kA:.3e}")
A1 = 1.28e-18 * math.sqrt(7e-11) / 5e9
A1_typo = 1.28e-18 * math.sqrt(7e-11 / 5e9)
pr("A for F=7e-11, nu=5e9 with sqrt(F)/nu (DIRD 2e-33)", f"{A1:.2e}")
pr("  with eq (2) as printed sqrt(F/nu)", f"{A1_typo:.2e}")
A2 = 1.28e-18 * math.sqrt(3.0) / 1.2e14
pr("A for F=3 W/m^2, nu=1.2e14 (DIRD 1.8e-32)", f"{A2:.2e}")
B, S, N = 1e5, 3.0, 1e-8
pr("Shannon B log2(1+S/N) (DIRD 1.9e6 bps)", f"{B*math.log2(1+S/N):.3e}")
pr("  B ln(1+S/N)", f"{B*math.log(1+S/N):.3e}")
pr("received flux 1.48e14/(7e6)^2 (DIRD 3 W/m^2)", f"{1.48e14/7e6**2:.2f}")

print("== Pfenning-Ford 1997 ==")
Rb = 100.0
for vb in (1.0,):
    Dl = 100 * vb * lP
    E_len = (vb**2 / 12) * (Rb**2 / Dl + Dl / 12)          # geometrised length
    pr("E = (v^2/12)(R^2/Delta) with Delta=100 v lP, R=100 m [kg]", f"{E_len*c**2/G:.2e}  (paper: 6.2e62 v kg)")
    pr("   same in Planck masses", f"{E_len*c**2/G/mP:.2e}  (paper: 6.2e70)")
    E1 = (vb**2 / 12) * (Rb**2 / 1.0)
    pr("Delta = 1 m: E [Msun] (paper: ~1/4 Msun)", f"{E1*c**2/G/Msun:.3f}")
pr("Lentz E~C v^2 R^2/w, R=100, w=1, C=1 [Msun]", f"{1e4*c**2/G/Msun:.2f} (paper: few x 0.1 Msun, C~O(1))")

print("== Ford-Roman QI (massless scalar, Lorentzian sampling) ==")
for t0 in (1e-15, 1e-9, 1.0):
    pr(f"  rho_min = -3 hbar c/(32 pi^2 (c t0)^4), t0={t0} s [J/m^3]", f"{-3*hbar*c/(32*math.pi**2*(c*t0)**4):.3e}")

print("== Fuchs et al 2024 shell ==")
M = 4.49e27
pr("M/MJ (paper 2.365)", f"{M/MJ:.3f}")
pr("2GM/(c^2 R2), R2=20 m", f"{2*G*M/(c**2*20):.3f}   r_s = {2*G*M/c**2:.2f} m")

print("== Casimir (ideal plates, T=0) ==")
for d in (1e-6, 1e-7, 1e-8):
    pr(f"  d={d:g} m: u=-pi^2 hbar c/(720 d^4), P=-pi^2 hbar c/(240 d^4)",
       f"{-math.pi**2*hbar*c/(720*d**4):.3e} J/m^3, {-math.pi**2*hbar*c/(240*d**4):.3e} Pa")
Rs, d = 100e-6, 100e-9
pr("sphere-plate PFA F=-pi^3 hbar c R/(360 d^3), R=100um,d=100nm", f"{-math.pi**3*hbar*c*Rs/(360*d**3):.3e} N")
# weight of Casimir energy (Fulling et al.): F = |E| g / c^2, 1 cm^2 plates at 10 nm
E = -math.pi**2 * hbar * c / (720 * (1e-8)**3) * 1e-4
pr("Casimir energy of 1 cm^2 plates at d=10 nm [J]; weight |E| g/c^2 [N]", f"{E:.3e}, {abs(E)*9.81/c**2:.3e}")

print("== PV / exponential metric vs Schwarzschild ==")
pr("light bending 4GM/(c^2 R_sun) [arcsec]", f"{4*G*Msun/(c**2*Rsun)*206264.806:.4f}")
a_m, e_m, P_m = 5.7909e10, 0.20563, 87.969
dphi = 6 * math.pi * G * Msun / (c**2 * a_m * (1 - e_m**2))
pr("Mercury perihelion 6 pi GM/(c^2 a(1-e^2)) [arcsec/century]", f"{dphi*206264.806*36525/P_m:.2f}")
pr("Pound-Rebka g h/c^2, h=22.5 m", f"{9.81*22.5/c**2:.3e}")

def deflection(n, x):
    """exact bending for isotropic metric with index n(u) (u = M/r, M=1), closest approach u0=x.
    phi = 2 int_0^{u0} b du / sqrt(n(u)^2 - b^2 u^2), b = n(u0)/u0; substitute u = u0 (1 - t^2)."""
    u0 = x
    b = n(u0) / u0
    def f(t):
        u = u0 * (1 - t * t)
        return 2 * b * (2 * u0 * t) / math.sqrt(n(u) ** 2 - (b * u) ** 2)
    # composite 5-point Gauss-Legendre (open rule, so the removable 0/0 at t=0 is never evaluated)
    xs = (0.0, -0.5384693101056831, 0.5384693101056831, -0.9061798459386640, 0.9061798459386640)
    ws = (0.5688888888888889, 0.4786286704993665, 0.4786286704993665, 0.2369268850561891, 0.2369268850561891)
    N = 400
    hstep, s = 1.0 / N, 0.0
    for i in range(N):
        mid = (i + 0.5) * hstep
        s += sum(w * f(mid + 0.5 * hstep * xi) for xi, w in zip(xs, ws))
    return s * 0.5 * hstep - math.pi

n_exp = lambda u: math.exp(2 * u)                                 # sqrt(B/A) = e^{2M/r}
n_sch = lambda u: (1 + u / 2) ** 3 / (1 - u / 2)                  # isotropic Schwarzschild
for x in (1e-2, 3e-3):
    de, ds = deflection(n_exp, x), deflection(n_sch, x)
    # in terms of M/b: b = n(u0)/u0  -> m = M/b
    me, ms = x / n_exp(x), x / n_sch(x)
    pr(f"  r0=M/{1/x:.0f}: (alpha-4m)/m^2 exp-metric (expect 4pi={4*math.pi:.4f})", f"{(de-4*me)/me**2:.4f}")
    pr(f"  r0=M/{1/x:.0f}: (alpha-4m)/m^2 Schwarzschild (expect 15pi/4={15*math.pi/4:.4f})", f"{(ds-4*ms)/ms**2:.4f}")
def c2(n, xs=(3e-3, 1e-3)):
    """(alpha-4m)/m^2 = c2 + c3 m: extrapolate linearly to m -> 0"""
    pts = [(x / n(x), (deflection(n, x) - 4 * x / n(x)) / (x / n(x)) ** 2) for x in xs]
    (m1, y1), (m2, y2) = pts
    return y2 - (y1 - y2) / (m1 - m2) * m2
ce, cs = c2(n_exp), c2(n_sch)
pr("extrapolated c2: exp-metric (4pi) / Schwarzschild (15pi/4)", f"{ce:.4f} / {cs:.4f}")
assert abs(ce - 4 * math.pi) < 2e-2 and abs(cs - 15 * math.pi / 4) < 2e-2
mb = G * Msun / (c**2 * Rsun)
pr("Sun-grazing 2nd-order difference (pi/4)(M/b)^2 [micro-arcsec]", f"{math.pi/4*mb**2*206264.806e6:.3f}")
pr("exponential metric throat: r=GM/c^2, areal radius e*GM/c^2 for 1 Msun [m]", f"{math.e*G*Msun/c**2:.1f}")

print("== DIRD 24 ZPE density with Planck cutoff ==")
wP = math.sqrt(c**5 / (hbar * G))
pr("omega_P = sqrt(c^5/(hbar G)) [1/s]", f"{wP:.3e}")
pr("rho = hbar w_P^4/(8 pi^2 c^3)  (integral of hbar w^3/(2pi^2c^3))", f"{hbar*wP**4/(8*math.pi**2*c**3):.2e} J/m^3")
pr("rho if nu_P := omega_P is used as a frequency (2pi^2 hbar nu^4/c^3)", f"{2*math.pi**2*hbar*wP**4/c**3:.2e} J/m^3 (DIRD ~1e113)")

print("== Visser 1989 cubic wormhole ==")
pr("edge tension mu = T = -c^4/(8G) [J/m] (paper -1.52e43)", f"{-c**4/(8*G):.3e}")
pr("  in Planck units (paper -1/8 m_P/l_P): mu/(m_P c^2/l_P)", f"{-c**4/(8*G)/(mP*c**2/lP):.3f}")
pr("  total edge energy, 1 m cube (12 edges) [J], [kg]", f"{-12*c**4/(8*G):.3e}, {-12*c**2/(8*G):.3e}")
