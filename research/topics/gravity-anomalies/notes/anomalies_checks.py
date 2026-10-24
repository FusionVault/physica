# Recompute key numbers for open_gravity_anomalies.md. Inputs are published values (sources in the notes).
import math

def solve(A, b):  # Gauss-Jordan, small dense systems
    n = len(A); M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(M[r][c])); M[c], M[p] = M[p], M[c]
        for r in range(n):
            if r != c:
                f = M[r][c] / M[c][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return [M[i][n] / M[i][i] for i in range(n)]

print("== 1. G: CODATA 2018/2022 Table XXX (16 inputs, 1e-11 m^3 kg^-1 s^-2) ==")
G = [("NIST-82",6.67248,0.00043),("TR&D-96",6.6729,0.0005),("LANL-97",6.67398,0.00070),
     ("UWash-00",6.674255,0.000092),("BIPM-01",6.67559,0.00027),("UWup-02",6.67422,0.00098),
     ("MSL-03",6.67387,0.00027),("HUST-05",6.67222,0.00087),("UZur-06",6.67425,0.00012),
     ("HUST-09",6.67349,0.00018),("BIPM-14",6.67554,0.00016),("LENS-14",6.67191,0.00099),
     ("UCI-14",6.67435,0.00013),("HUSTT-18",6.674184,0.000078),("HUSTA-18",6.674484,0.000077),
     ("JILA-18",6.67260,0.00025)]
idx = {n:i for i,(n,_,_) in enumerate(G)}
N = len(G); C = [[0.0]*N for _ in range(N)]
for i,(_,_,u) in enumerate(G): C[i][i] = u*u
for a,b,r in [("NIST-82","LANL-97",0.351),("HUST-05","HUST-09",0.134),("HUST-09","HUSTT-18",0.068)]:
    i,j = idx[a],idx[b]; C[i][j] = C[j][i] = r*G[i][2]*G[j][2]
y = [v for _,v,_ in G]
w = solve(C, [1.0]*N)                       # C^-1 * 1
Gm = sum(wi*yi for wi,yi in zip(w,y))/sum(w)
uG = 1/math.sqrt(sum(w))
r = [yi-Gm for yi in y]; chi2 = sum(ri*zi for ri,zi in zip(r, solve(C, r)))
print(f"GLS mean = {Gm:.6f}, u(unexpanded) = {uG:.6f} ({uG/Gm*1e6:.1f} ppm)")
print(f"chi2 = {chi2:.1f} for nu = {N-1}; Birge ratio sqrt(chi2/nu) = {math.sqrt(chi2/(N-1)):.2f}")
print(f"expanded x3.9: u = {3.9*uG:.5f} ({3.9*uG/Gm*1e6:.1f} ppm); chi2/3.9^2 = {chi2/3.9**2:.1f}")
for n in ("BIPM-14","JILA-18","BIPM-01"):
    i = idx[n]; print(f"  normalized residual {n}: {(y[i]-Gm)/G[i][2]:+.2f}")
print(f"spread max-min = {(max(y)-min(y))/Gm*1e6:.0f} ppm ({max(y)-min(y):.5f})")
t,a = G[idx['HUSTT-18']], G[idx['HUSTA-18']]
print(f"HUST TOS vs AAF: diff {(a[1]-t[1])/Gm*1e6:.0f} ppm = {(a[1]-t[1])/math.hypot(t[2],a[2]):.1f} sigma")
print(f"Quinn 2013 6.67545(18) -> 2014 6.67554(16): +{(6.67554-6.67545)/6.67545*1e6:.1f} ppm")
for lab,v,u in [("NIST-2026 (Metrologia abstract)",6.67366,0.00020),("NIST-2026 (APS Physics news)",6.67387,57e-6*6.67387)]:
    print(f"{lab}: u_r = {u/v*1e6:.0f} ppm; vs BIPM-14 {(v-6.67554)/6.67554*1e6:+.0f} ppm "
          f"= {(v-6.67554)/math.hypot(u,0.00016):+.1f} sigma; vs CODATA {(v-6.67430)/6.67430*1e6:+.0f} ppm")
print(f"Rosi 2014 vs CODATA2022: {(6.67191-6.67430)/6.67430*1e6:+.0f} ppm")

print("\n== 2. Flyby anomaly, Anderson et al. 2008 ==")
wE, RE, c = 7.292115e-5, 6371e3, 299792458.0
K = 2*wE*RE/c; print(f"K = 2 wE RE / c = {K:.4e}")
fly = [("GLL-I",8.949,-12.52,-34.15,3.92,0.3,4.12),("GLL-II",8.877,-34.26,-4.87,-4.6,1,-4.67),
       ("NEAR",6.851,-20.76,-71.96,13.46,0.01,13.28),("Cassini",16.01,-12.92,-4.99,-2,1,-1.07),
       ("Rosetta-I",3.863,-2.81,-34.29,1.80,0.03,2.07),("MESSENGER",4.056,31.44,-31.92,0.02,0.01,0.06)]
for n,v,di,do,obs,s,pred in fly:
    p = K*v*1e6*(math.cos(math.radians(di))-math.cos(math.radians(do)))  # km/s*1e6 -> mm/s
    print(f"  {n:10s} V={v:6.3f} km/s  dV_pred={p:+7.3f} mm/s (paper {pred:+.2f}), obs {obs:+.2f}+-{s}  dV/V obs={obs/(v*1e6):+.3e}")

print("\n== 3. Gravitational absorption (Majorana): F = F_N exp(-h * int rho dl) ==")
col_E = 2*6371e3*5514.0; col_M = 2*1737.4e3*3344.0
print(f"Earth diametral column {col_E:.3e} kg/m^2, Moon {col_M:.3e} kg/m^2")
for lab,h in [("Majorana 1920 h=6.7e-12 cm^2/g",6.7e-13),("Zurich lab bound 4.3e-15",4.3e-15),
              ("Eckhardt 1990 LLR <1e-22",1e-22),("Williams+ 2009 (3+-5)e-22, 1sig upper 8e-22",8e-22)]:
    print(f"  {lab}: Earth-through fraction 1-exp(-h*col) = {1-math.exp(-h*col_E):.3e}; Moon = {1-math.exp(-h*col_M):.3e}")
gS = 1.32712440018e20/1.495978707e11**2
print(f"Sun's gravity at 1 au = {gS:.4e} m/s^2")
h_wang = (7.0e-8/gS)/col_M
print(f"Wang 2000 7.0e-8 m/s^2 read as lunar shielding of Sun -> h = {h_wang:.2e} m^2/kg = {h_wang/1e-22:.1e} x LLR bound")
print(f"LLR-bound eclipse signal: {1e-22*col_M*gS:.1e} m/s^2 = {1e-22*col_M*gS*1e8:.1e} uGal")

print("\n== 4. Pioneer ==")
aP = 8.74e-10; H0 = 70e3/3.0856775814913673e22
print(f"cH0 (H0=70) = {c*H0:.3e} m/s^2; aP/cH0 = {aP/(c*H0):.2f}")
print(f"thermal fit a0=7.4e-10 / aP = {7.4/8.74:.2f}; Doppler a0 10.1e-10 tau 28.8 yr -> at t0+20yr {10.1*2**(-20/28.8):.2f}e-10")
m = 241.0  # kg, approx Pioneer mass used only for scale (assumption)
print(f"power needed for aP via collimated photons on {m} kg: P = m a c = {m*aP*c:.0f} W")

print("\n== 5. MOND / RAR (McGaugh+2016): g_obs = g_bar / (1 - exp(-sqrt(g_bar/g+))) ==")
gd = 1.20e-10
print(f"cH0/(2pi) = {c*H0/(2*math.pi):.3e}")
for gb in (1e-12,1e-11,1e-10,1e-9,1e-8):
    go = gb/(1-math.exp(-math.sqrt(gb/gd)))
    print(f"  g_bar={gb:.0e}: g_obs={go:.4e}  ratio={go/gb:.4f}  deep-MOND sqrt(g_bar a0)={math.sqrt(gb*gd):.4e}")
M = 1.0e11*1.98847e30; G0 = 6.67430e-11
print(f"BTFR v^4 = G M a0 for M=1e11 Msun: v = {(G0*M*gd)**0.25/1e3:.1f} km/s; MOND radius sqrt(GM/a0) = {math.sqrt(G0*M/gd)/3.0857e19:.1f} kpc")
Msun2 = 2*1.98847e30
r0 = math.sqrt(G0*Msun2/gd)/1.495978707e11
print(f"wide binary (2 Msun total) separation where g_N = a0: {r0:.0f} au")

print("\n== 6. Wide binaries ==")
print(f"Chae 2023: 10^(sqrt2*0.109) = {10**(math.sqrt(2)*0.109):.3f}; 10^(sqrt2*0.034) = {10**(math.sqrt(2)*0.034):.3f}")
for G_,lab in [(0.085,"Chae 2025 MOND regime"),(0.063,"Chae 2025 MOND+transition"),(0.102,"Chae+ 2026 <10^-9.1"),(0.132,"Chae+ 2026 <10^-9.5"),(0.07,"MOND pred (Chae)")]:
    print(f"  Gamma={G_}: gamma_g = 10^(2 Gamma) = {10**(2*G_):.3f}  ({lab})")
print(f"Banik+2024: MOND expects ~20% velocity boost -> gravity boost {1.2**2:.2f}")

print("\n== 7. Antimatter / EP ==")
u = math.hypot(0.13,0.16); print(f"ALPHA-g a_g = 0.75 +- {u:.3f} g; (1-0.75)/u = {(1-0.75)/u:.2f} sigma from +1g; (0.75+1)/u = {(0.75+1)/u:.1f} sigma from -1g")
u = math.hypot(2.3,1.5); print(f"MICROSCOPE eta = -1.5 +- {u:.2f} e-15 (quadrature); |eta|<{abs(-1.5)+2*u:.1f}e-15 approx 2-sigma envelope")
print(f"Eot-Wash 2008 Be-Ti: Delta a = eta * g_N toward source; 1.8e-13*9.8/ (check vs 3.1e-15 m/s^2) -> implied g_eff {3.1e-15/1.8e-13:.4f} m/s^2 (horizontal component of Earth's pull at Seattle)")

print("\n== 8. Tate Cooper-pair mass ==")
meas, th = 1.000084, 0.999992
d = (meas-th)/th
print(f"discrepancy (m*-m_th)/m_th = {d*1e6:.0f} ppm = {(meas-th)/21e-6:.1f} sigma; Tajmar-de Matos Bg = 2*omega*{d:.3e} = {2*d:.3e} omega")
print(f"London moment B = 2 m omega / e : per rad/s B = {2*9.1093837e-31/1.602176634e-19:.4e} T")
