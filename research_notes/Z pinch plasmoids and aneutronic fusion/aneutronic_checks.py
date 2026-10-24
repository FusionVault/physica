"""Recompute the key numbers for aneutronic_fusion_direct_conversion.md.

Standard library only. Run: python3 -I aneutronic_checks.py
SI units unless a name says otherwise. CODATA 2018 constants; NIST atomic masses (AME2020).
Every assert ties a computed number to a published one (tolerance stated in the line).
"""
import math

# --- constants (CODATA 2018) ---
c = 299792458.0; e = 1.602176634e-19; me = 9.1093837015e-31; eps0 = 8.8541878128e-12
alpha = 7.2973525693e-3; re = 2.8179403262e-15; mu0 = 1.25663706212e-6
u_MeV = 931.49410242; u_kg = 1.66053906660e-27
keV = 1e3 * e; MeV = 1e6 * e; mec2_keV = 510.99895

def show(k, v, fmt="{:.4g}"):
    print(f"{k:66s} {fmt.format(v)}")

def close(a, b, rel):
    return abs(a - b) <= rel * abs(b)

# --- 1. Q-values from atomic masses (NIST/AME2020, u); electron counts balance in every line ---
M = dict(n=1.00866491595, H1=1.00782503223, H2=2.01410177812, H3=3.0160492779,
         He3=3.0160293201, He4=4.00260325413, Li6=6.0151228874, Li7=7.0160034366,
         B11=11.00930536, N14=14.00307400443,
         C11=11.0114336)  # shortcut: C11 from AME2020 mass excess 10650.3 keV (recalled, not fetched)
def Q(lhs, rhs):
    return (sum(M[x] for x in lhs) - sum(M[x] for x in rhs)) * u_MeV
reactions = {
    "p + 11B -> 3 alpha":        (["H1", "B11"], ["He4"] * 3, 8.68),
    "D + 3He -> alpha + p":      (["H2", "He3"], ["He4", "H1"], 18.35),
    "3He + 3He -> alpha + 2p":   (["He3", "He3"], ["He4", "H1", "H1"], 12.86),
    "p + 6Li -> alpha + 3He":    (["H1", "Li6"], ["He4", "He3"], 4.02),
    "D + D -> T + p":            (["H2", "H2"], ["H3", "H1"], 4.03),
    "D + D -> 3He + n":          (["H2", "H2"], ["He3", "n"], 3.27),
    "D + T -> alpha + n":        (["H2", "H3"], ["He4", "n"], 17.59),
    "p + 7Li -> 2 alpha":        (["H1", "Li7"], ["He4", "He4"], 17.35),
    "D + 6Li -> 2 alpha":        (["H2", "Li6"], ["He4", "He4"], 22.37),
    "alpha + 11B -> 14N + n":    (["He4", "B11"], ["N14", "n"], 0.158),
    "p + 11B -> 11C + n":        (["H1", "B11"], ["C11", "n"], -2.765),
}
Qv = {}
for name, (l, r, pub) in reactions.items():
    Qv[name] = Q(l, r); show(f"Q [{name}] [MeV]  (published {pub})", Qv[name], "{:.4f}")
    assert close(Qv[name], pub, 0.01) or abs(Qv[name] - pub) < 0.01
# two-body product split (zero initial KE, non-relativistic): E_a = Q m_b/(m_a+m_b)
def split(q, ma, mb):
    return q * mb / (ma + mb), q * ma / (ma + mb)
Ep, Ea = split(Qv["D + 3He -> alpha + p"], M["H1"], M["He4"])
show("D-3He: proton energy [MeV] (published 14.7)", Ep); show("D-3He: alpha energy [MeV] (published 3.6)", Ea)
assert close(Ep, 14.7, 0.01) and close(Ea, 3.6, 0.03)
Ea6, Eh6 = split(Qv["p + 6Li -> alpha + 3He"], M["He4"], M["He3"])
show("p-6Li: alpha / 3He energies [MeV] (published 1.7 / 2.3)", Ea6); show("   3He", Eh6)
Et, Ep2 = split(Qv["D + D -> T + p"], M["H3"], M["H1"])
Eh, En = split(Qv["D + D -> 3He + n"], M["He3"], M["n"])
show("D-D(p): T / p energies [MeV] (1.01 / 3.02)", Et); show("   p", Ep2)
show("D-D(n): 3He / n energies [MeV] (0.82 / 2.45)", Eh); show("   n", En)
assert close(En, 2.45, 0.01) and close(Ep2, 3.02, 0.01)
show("p-11B: mean alpha energy Q/3 [MeV]", Qv["p + 11B -> 3 alpha"] / 3)

# --- 2. Gamow energy E_G = 2 m_r c^2 (pi alpha Z1 Z2)^2 ; B_G = sqrt(E_G) ---
def mr_keV(a, b):
    return M[a] * M[b] / (M[a] + M[b]) * u_MeV * 1e3
def EG(a, b, z1, z2):
    return 2 * mr_keV(a, b) * (math.pi * alpha * z1 * z2) ** 2
for (a, b, z1, z2, lab) in [("H1", "B11", 1, 5, "p-11B (NS: E_G=22589 keV, mrc2=859526)"),
                            ("H2", "He3", 1, 2, "D-3He (BH: B_G=68.7508, mrc2=1124572)"),
                            ("H2", "H3", 1, 1, "D-T  (BH: B_G=34.3827)"),
                            ("H2", "H2", 1, 1, "D-D  (BH: B_G=31.3970)"),
                            ("He3", "He3", 2, 2, "3He-3He"), ("H1", "Li6", 1, 3, "p-6Li")]:
    show(f"E_G [keV] {lab}", EG(a, b, z1, z2), "{:.1f}"); show("   sqrt(E_G) [keV^0.5]; m_r c^2 [keV]", math.sqrt(EG(a, b, z1, z2)), "{:.4f}")
    show("   m_r c^2 [keV] (atomic masses; BH/NS use nuclear)", mr_keV(a, b), "{:.0f}")
assert close(EG("H1", "B11", 1, 5), 22589, 0.002)
assert close(math.sqrt(EG("H2", "He3", 1, 2)), 68.7508, 0.002)

# --- 3. Bosch-Hale / Nevins-Swain / Tentori-Belloni reactivity form ---
def peres(T, BG2, mrc2, C):
    """<sigma v> = C1 theta sqrt(xi/(mrc2 T^3)) exp(-3 xi); T in keV; units of C1 set the output."""
    C1, C2, C3, C4, C5, C6, C7 = C
    th = T / (1 - T * (C2 + T * (C4 + T * C6)) / (1 + T * (C3 + T * (C5 + T * C7))))
    xi = (BG2 / (4 * th)) ** (1 / 3)
    return C1 * th * math.sqrt(xi / (mrc2 * T ** 3)) * math.exp(-3 * xi)

# p-11B: Nevins & Swain 2000 (as tabulated in Tentori & Belloni 2023 Table 2; output m^3/s)
NS = (4.4467e-14, -5.9357e-2, 2.0165e-1, 1.0404e-3, 2.7621e-3, -9.1653e-6, 9.8305e-7)
TB = (9.7827e-14, -5.1610e-2, 1.3240e-1, 3.8446e-4, 1.2499e-3, -5.6715e-6, 1.1615e-6)
def sv_NS(T, resonance=True):
    v = peres(T, 22589.0, 859526.0, NS)
    if resonance:  # NS 148 keV narrow-resonance term, cm^3/s -> m^3/s
        v += 5.41e-15 * T ** -1.5 * math.exp(-148.0 / T) * 1e-6
    return v
def sv_TB(T):
    return peres(T, 22589.0, 859526.0, TB)
print("\n p-11B reactivity [m^3/s]   T[keV]   NS(+res)   NS(no res)   TB    TB/NS-1")
for T in (70, 100, 150, 200, 250, 300, 400, 500):
    a, b, t = sv_NS(T), sv_NS(T, False), sv_TB(T)
    print(f"   {T:6.0f} {a:11.4e} {b:11.4e} {t:11.4e} {t / a - 1:8.3f}")
# published (Ahmad et al. 2026, Frontiers, Fig. 1 range 75-500 keV): NS 2.455e-23 / 3.741e-22, TB 2.777e-23 / 5.619e-22
# they match the NS fit WITHOUT the 148 keV resonance term, at 75 and 500 keV, to 4 digits:
assert close(sv_NS(75, False), 2.455e-23, 0.001) and close(sv_NS(500, False), 3.741e-22, 0.001)
assert close(sv_TB(75), 2.777e-23, 0.001) and close(sv_TB(500), 5.619e-22, 0.001)
# Tentori & Belloni text: +12% @100 keV, >20% @200, 30% @300, 50% @500
assert close(sv_TB(100) / sv_NS(100) - 1, 0.12, 0.15) and close(sv_TB(500) / sv_NS(500) - 1, 0.50, 0.05)
assert close(sv_TB(300) / sv_NS(300) - 1, 0.30, 0.10) and sv_TB(200) / sv_NS(200) - 1 > 0.20
Tpk = max(range(50, 501), key=sv_TB)
show("TB fit: T of max <sigma v> inside 50-500 keV [keV] (edge => still rising)", Tpk)

# Bosch & Hale 1992 reactivity coefficients (cm^3/s).
# shortcut: C1..C7 recalled from BH Table VII (paper not fetched); checked below by integrating the
# BH cross section (coefficients read from WarpX source) - agreement validates both.
BH = {
    "D3He": (68.7508 ** 2, 1124572.0, (5.51036e-10, 6.41918e-3, -2.02896e-3, -1.91080e-5, 1.35776e-4, 0, 0)),
    "DT":   (34.3827 ** 2, 1124656.0, (1.17302e-9, 1.51361e-2, 7.51886e-2, 4.60643e-3, 1.35000e-2, -1.06750e-4, 1.36600e-5)),
    "DDn":  (31.3970 ** 2, 937814.0, (5.43360e-12, 5.85778e-3, 7.68222e-3, 0, -2.96400e-6, 0, 0)),
    "DDp":  (31.3970 ** 2, 937814.0, (5.65718e-12, 3.41267e-3, 1.99167e-3, 0, 1.05060e-5, 0, 0)),
}
def sv_BH(r, T):
    BG2, mrc2, C = BH[r]
    return peres(T, BG2, mrc2, C)  # cm^3/s
# BH cross sections, S(E) Pade, E = CM energy in keV, sigma in mb (coefficients as in WarpX source)
XS = {
    "D3He": (68.7508, (5.7501e6, 2.5226e3, 4.5566e1, 0, 0), (-3.1995e-3, -8.5530e-6, 5.9014e-8, 0), 900.0),
    "DT":   (34.3827, (6.927e4, 7.454e8, 2.050e6, 5.2002e4, 0), (6.38e1, -9.95e-1, 6.981e-5, 1.728e-4), 550.0),
    "DDn":  (31.3970, (5.3701e4, 3.3027e2, -1.2706e-1, 2.9327e-5, -2.5151e-9), (0, 0, 0, 0), 4900.0),
    "DDp":  (31.3970, (5.5576e4, 2.1054e2, -3.2638e-2, 1.4987e-6, 1.8181e-10), (0, 0, 0, 0), 4900.0),
}
def sigma_mb(r, E):
    BG, A, B, _ = XS[r]
    S = A[0] + E * (A[1] + E * (A[2] + E * (A[3] + E * A[4])))
    S /= 1 + E * (B[0] + E * (B[1] + E * (B[2] + E * B[3])))
    return S / (E * math.exp(BG / math.sqrt(E)))
def sv_num(r, T, mrc2):
    """Maxwellian <sigma v> [cm^3/s] by Simpson integration of sigma(E) E exp(-E/T) to the fit's E_max."""
    Emax = min(XS[r][3], 40 * T); N = 20000; h = Emax / N
    s = 0.0
    for i in range(1, N + 1):
        E = i * h; w = 1 if i == N else (4 if i % 2 else 2)
        s += w * sigma_mb(r, E) * E * math.exp(-E / T)
    integral = s * h / 3 * 1e-27 * keV ** 2          # mb keV^2 -> cm^2 J^2
    mu = mrc2 * keV / c ** 2
    return math.sqrt(8 / (math.pi * mu)) * (T * keV) ** -1.5 * integral * 100  # m/s -> cm/s
print("\n Bosch-Hale <sigma v> [cm^3/s]: fit vs numerical integral of BH sigma(E)")
for r, Ts in (("D3He", (10, 20, 50, 100)), ("DT", (5, 10, 20, 50)), ("DDn", (10, 50, 100)), ("DDp", (10, 50, 100))):
    for T in Ts:
        f, n = sv_BH(r, T), sv_num(r, T, BH[r][1])
        print(f"   {r:5s} T={T:5.0f} keV  fit {f:10.4e}  integral {n:10.4e}  ratio {f / n:.4f}")
        assert close(f, n, 0.05), (r, T)  # BH quote ~2% max error for each fit; 5% bounds their sum
show("D-T <sigma v> at 10 keV [cm^3/s] (BH 1992 table 1.136e-16; recalled)", sv_BH("DT", 10))

# --- 4. D-3He neutron fraction from D-D side reactions ---
def d3he_neutron_fraction(T, xD, burn_T=False):
    """xD = n_D/(n_D+n_He3). Fraction of fusion power carried by neutrons. burn_T: bred T burns at once."""
    nD, nH = xD, 1 - xD
    p_dhe = nD * nH * sv_BH("D3He", T) * Qv["D + 3He -> alpha + p"]
    p_ddn = 0.5 * nD * nD * sv_BH("DDn", T) * Qv["D + D -> 3He + n"]
    p_ddp = 0.5 * nD * nD * sv_BH("DDp", T) * Qv["D + D -> T + p"]
    pn = 0.5 * nD * nD * sv_BH("DDn", T) * En
    tot = p_dhe + p_ddn + p_ddp
    if burn_T:
        tot += 0.5 * nD * nD * sv_BH("DDp", T) * Qv["D + T -> alpha + n"]
        pn += 0.5 * nD * nD * sv_BH("DDp", T) * 14.07
    return pn / tot
print("\n D-3He neutron power fraction  T[keV]  xD  (no T burn / bred T burns)")
for T in (50, 70, 100):
    for xD in (0.5, 0.25):
        print(f"   T={T:4.0f}  xD={xD:.2f}  {d3he_neutron_fraction(T, xD):.4f}  {d3he_neutron_fraction(T, xD, True):.4f}")

# --- 5. Bremsstrahlung and synchrotron ---
W0 = alpha * re ** 2 * me * c ** 3                  # e^6/(hbar m c^2) in SI form, W m^3
show("W0 = alpha r_e^2 m_e c^3 [W m^3]", W0)
gB = 2 * math.sqrt(3) / math.pi
coef = 16 / 3 * math.sqrt(2 * math.pi / 3) * gB
show("non-rel. coefficient (16/3)sqrt(2pi/3)(2sqrt3/pi)  (Putvinski: 8.511)", coef)
assert close(coef, 8.511, 1e-3)
Cbr = coef * W0 / math.sqrt(mec2_keV)              # W m^3 keV^-1/2
show("P_br/(Z_eff n_e^2 sqrt(T_keV)) [W m^3]  (NRL 1.69e-32 W cm^3 eV^-1/2 = 5.344e-37)", Cbr)
assert close(Cbr, 1.69e-32 * 1e-6 * math.sqrt(1e3), 0.005)
def Pbr_putvinski(ne, Te, Zeff):
    """Putvinski et al. 2019 Eq. 7 (fit to QED results, <1%): W/m^3, Te keV."""
    t = Te / mec2_keV
    return coef * W0 * ne ** 2 * math.sqrt(t) * (Zeff * (1 + 1.78 * t ** 1.34) + 2.12 * t * (1 + 1.1 * t + t ** 2 - 1.25 * t ** 2.5))
def Pbr_rider(ne, Te, Zeff, pref=5.172e-37):
    """Rider form as quoted by Wang et al. 2026 Eq. 10 (5.172e-43 MW m^3 => 5.172e-37 W m^3, Te keV)."""
    t = Te / mec2_keV
    return pref * ne ** 2 * math.sqrt(Te) * (Zeff * (1 + 0.7936 * t + 1.874 * t ** 2) + 3 / math.sqrt(2) * t)
show("Rider prefactor / NRL prefactor", 5.172e-37 / Cbr)
for Te in (50, 150, 300):
    show(f"Putvinski/Rider brems ratio, Z_eff=3, Te={Te} keV (same prefactor)", Pbr_putvinski(1, Te, 3) / Pbr_rider(1, Te, 3, Cbr))
    show(f"relativistic enhancement Putvinski(T)/nonrel, Z_eff=3, Te={Te}", Pbr_putvinski(1, Te, 3) / (Cbr * 3 * math.sqrt(Te)))

Csync = e ** 4 * keV / (3 * math.pi * eps0 * me ** 3 * c ** 3)   # Larmor, <v_perp^2>=2kT/m
show("P_sync/(B^2 n_e T_keV) [W m^3 T^-2 keV^-1] (NRL 6.21e-28 cgs-eV = 6.21e-17)", Csync)
assert close(Csync, 6.21e-17, 0.005)

# --- 6. p-11B: fusion vs bremsstrahlung ---
EpB = Qv["p + 11B -> 3 alpha"] * MeV
def pB_ratio(Ti, Te, fB, sv=sv_TB, brem=Pbr_putvinski):
    """fB = n_B/(n_p+n_B). Returns P_fus/P_br (density-independent)."""
    npp, nB = 1 - fB, fB
    ne = npp + 5 * nB; Zeff = (npp + 25 * nB) / ne
    return npp * nB * sv(Ti) * EpB / brem(ne, Te, Zeff)
best = max((pB_ratio(300, 300, f / 1000), f / 1000) for f in range(20, 500))
show("Te=Ti=300 keV, TB: max P_fus/P_br over fB", best[0]); show("   at fB = n_B/n_i", best[1])
print("\n p-11B P_fus/P_br, fB=0.15 (Putvinski optimum)   Ti   Te/Ti   NS   TB")
for Ti in (150, 200, 300, 400, 500):
    for r in (1.0, 0.5, 0.33, 0.25):
        print(f"   Ti={Ti:4.0f} Te/Ti={r:4.2f}  NS {pB_ratio(Ti, r * Ti, 0.15, sv_NS):6.3f}  TB {pB_ratio(Ti, r * Ti, 0.15):6.3f}")
for sv, lab in ((sv_NS, "NS"), (sv_TB, "TB")):
    lo = next((T for T in range(70, 501) if pB_ratio(T, T, 0.15, sv) >= 1), None)
    show(f"{lab}: lowest T (Te=Ti, fB=0.15) with P_fus >= P_br [keV] (None = never)", lo if lo else float("nan"))
show("Rider-form check: P_br/P_fus at Ti=300, Te=Ti, fB=0.1, NS", 1 / pB_ratio(300, 300, 0.1, sv_NS, lambda n, T, Z: Pbr_rider(n, T, Z)))
show("Liu 2024 check: P_br/P_fus at Ti=300, Te=23.33, fB=0.1, NS (published 0.5 at Rider's sigma)", 1 / pB_ratio(300, 23.33, 0.1, sv_NS, lambda n, T, Z: Pbr_rider(n, T, Z)))
show("   same with TB (published 0.27 at 2016 sigma)", 1 / pB_ratio(300, 23.33, 0.1, sv_TB, lambda n, T, Z: Pbr_rider(n, T, Z)))
# synchrotron vs brems at fixed beta: B^2 = 2 mu0 n (Ti_i + Te)/beta
def sync_over_brem(Ti, Te, fB, beta, ne=1e20):
    npp, nB = (1 - fB), fB; ne_rel = npp + 5 * nB; Zeff = (npp + 25 * nB) / ne_rel
    ni = ne / ne_rel
    B2 = 2 * mu0 * (ni * Ti + ne * Te) * keV / beta
    return Csync * B2 * ne * Te / Pbr_putvinski(ne, Te, Zeff)
for beta in (1.0, 0.5, 0.1):
    show(f"p-11B Ti=300, Te=150, fB=0.15: P_sync(no reflection)/P_br at beta={beta}", sync_over_brem(300, 150, 0.15, beta))
show("   field needed at beta=1, n_e=1e20 [T]", math.sqrt(2 * mu0 * (1e20 / 1.6 * 300 + 1e20 * 150) * keV))

# --- 7. Direct conversion: N-stage electrostatic collector on a spread-energy ion flux ---
def eta_stages(N, shift, Emax=40.0, n=4000):
    """Ion flux f(E) ~ (E-shift) exp(-(E-shift)) for E>shift (units of T). Collectors at V1<..<VN; an ion is
    caught by the highest V <= E and returns qV. Coordinate search for the best voltages. Returns max efficiency."""
    h = (Emax - shift) / n
    Es = [shift + (i + 0.5) * h for i in range(n)]
    f = [(E - shift) * math.exp(-(E - shift)) for E in Es]
    tot = sum(fi * E for fi, E in zip(f, Es))
    def eff(V):
        s = 0.0
        for fi, E in zip(f, Es):
            caught = max((v for v in V if v <= E), default=0.0)
            s += fi * caught
        return s / tot
    V = [shift + (k + 1) * 1.0 for k in range(N)]
    step = 0.5
    while step > 0.01:
        improved = False
        for k in range(N):
            for d in (-step, step):
                W = sorted(V[:k] + [max(0.0, V[k] + d)] + V[k + 1:])
                if eff(W) > eff(V) + 1e-9:
                    V, improved = W, True
        if not improved:
            step /= 2
    return eff(V)
print("\n ideal N-stage collector efficiency for flux (E-s)exp(-(E-s)), s = potential drop / T")
for s in (0.0, 3.0, 5.0):
    print("   s=%.0f  " % s + "  ".join(f"N={N}:{eta_stages(N, s, n=800):.3f}" for N in (1, 2, 4, 8)))
eta_mono = 1.0  # monoenergetic beam with one matched collector: no energy-spread loss (ideal)

# --- 8. Plant arithmetic ---
def Qeng_mcf(Qsci, eta_c, eta_d):
    """Steady magnetic: heat P_h absorbed; out = eta_c (P_fus + P_h); in = P_h/eta_d."""
    return eta_c * eta_d * (Qsci + 1)
def Qeng_icf(G, eta_c, eta_d):
    """Driver energy not recovered: out = eta_c G E_d; in = E_d/eta_d."""
    return eta_c * eta_d * G
show("Wikipedia/Entler: Q_E=1/(eta_heat f_recirc eta_elec (1-f_ch)), 0.7,0.2,0.4,0.2 (pub ~22)", 1 / (0.7 * 0.2 * 0.4 * 0.8))
assert close(1 / (0.7 * 0.2 * 0.4 * 0.8), 22, 0.02)
for eta_c, eta_d in ((0.40, 0.70), (0.40, 0.50), (0.60, 0.70), (0.80, 0.90)):
    show(f"MCF: Q_sci for Q_eng=1 at eta_c={eta_c}, eta_d={eta_d}", 1 / (eta_c * eta_d) - 1)
show("NIF 2022: Q_sci = 3.15/2.05", 3.15 / 2.05)
show("NIF 2022: yield/wall-plug (3.15 MJ / ~400 MJ)", 3.15 / 400)
show("NIF 2022: Q_eng if yield -> electricity at 0.4", 0.4 * 3.15 / 400)
show("NIF 2025: Q_sci = 8.6/2.08", 8.6 / 2.08)
show("NIF 2025: Q_eng at 0.4, assuming ~400 MJ wall plug (assumption)", 0.4 * 8.6 / 400)
show("Mehlhorn HB11 plant: fusion power 280 MJ x 5 Hz [MW]", 280 * 5)
show("   gross electric at 40% [MW]", 0.4 * 1400); show("   recirculating = 2 MJ x 5 Hz / 0.07 [MW]", 2 * 5 / 0.07)
show("   Q_eng (published 560/140 = 4)", Qeng_icf(140, 0.4, 0.07))
show("   Q_eng with 0.8 direct conversion", Qeng_icf(140, 0.8, 0.07))
show("   target gain for Q_eng=1 at 0.4, 7%", 1 / (0.4 * 0.07))
mB = 280e6 / EpB * M["B11"] * u_kg
show("   11B burned per 280 MJ pulse [mg] (published 3.525)", mB * 1e6)
# Helion-style pulsed recovery: E_out = eta_rec E_in + eta_dc f_ch Q_sci E_in
for eta_rec in (0.95, 0.90, 0.80, 0.60):
    show(f"pulsed inductive: Q_sci for Q_eng=1, eta_rec={eta_rec}, eta_dc=0.9, f_ch=0.95", (1 - eta_rec) / (0.9 * 0.95))
# p-11B with brems: charged share to DEC, X-rays to thermal
for fx in (0.5, 0.8):
    show(f"p-11B plant: eta_c with X-ray fraction {fx} of output (0.8 DEC, 0.4 thermal)", (1 - fx) * 0.8 + fx * 0.4)

# --- 9. Laser p-11B shots: fusion energy out / laser energy in (isotropic emission assumed) ---
def gain_from_alpha_per_sr(n_sr, E_laser_J):
    return 4 * math.pi * n_sr / 3 * EpB / E_laser_J
show("Margarone 2020 LFEX: 5e9 /sr, 1.4 kJ -> E_fus/E_laser (isotropic)", gain_from_alpha_per_sr(5e9, 1.4e3))
show("arXiv:2407.06729 PALS B18H22: 1.7e9 /sr, 550 J", gain_from_alpha_per_sr(1.7e9, 550))
show("Giuffrida 2020 PALS: 1e11 alphas total, ~600 J (energy assumed)", 1e11 / 3 * EpB / 600)
show("beam fusion: thick-target 6.56e-5 reactions per 1 MeV proton -> E_fus/E_p", 6.56e-5 * 8.68)
print("\nall asserts passed")
