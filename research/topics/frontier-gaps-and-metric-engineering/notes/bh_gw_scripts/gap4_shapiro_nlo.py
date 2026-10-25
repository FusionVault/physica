"""Gap 4: sizes of the NLO Shapiro terms for the Double Pulsar.

Inputs: Kramer et al. 2021 (arXiv:2112.06795) Table 3 (P_b, x, e_T, omega_0, r, z_s, nu, m_A, m_B).
Formulas: Kramer Eqs. ShaLO, Lambda_u, ShapRet, aberrationBend, aberrD, aberrationA;
lensing dLambda = 2rc/a_R with a_R = c(x+x_B)/s (Kramer, after Klioner & Zschocke 2010 Eq. 73);
retardation dLambda from Kopeikin & Schafer 1999 Eq. 130.
"""
import numpy as np

day = 86400.0
Pb = 0.1022515592973 * day
x = 1.415028603          # light-s
e = 0.087777023
om = np.radians(204.753686)
r = 6.162e-6             # s, measured Shapiro range
zs = 9.65
s = 1 - np.exp(-zs)
nu = 44.05406864196281
P = 1 / nu
mA, mB = 1.338185, 1.248868
xB = x * mA / mB          # R = m_A/m_B = x_B/x (Kramer, Sec. 6 end)

# Kramer Eq. (aberrationA): A = P x /(P_b sqrt(1-e^2) sin^2 i) ~ 3.65 us
A = P * x / (Pb * np.sqrt(1 - e**2) * s**2)
print(f"s = {s:.6f}, x_B = {xB:.6f} lt-s, A = {A*1e6:.4f} us (Kramer: 3.65 us)")

# lensing: dLambda_len = 2 r c / a_R = 2 r s /(x + x_B)
dlen = 2 * r * s / (x + xB)
print(f"a_R = {(x+xB)/s*2.99792458e8/1e3:.1f} km, dLambda_len = {dlen:.4e}")

# aberration NLO coefficient D = (P/pi) r/(x + x_B)
D = P / np.pi * r / (x + xB)
print(f"D = {D*1e9:.4f} ns")


def lam(u, om, s=s):
    return 1 - e * np.cos(u) - s * (np.sin(om) * (np.cos(u) - e) + np.sqrt(1 - e**2) * np.cos(om) * np.sin(u))


def dret(u, om, s=s):
    """Kopeikin & Schafer 1999 Eq. (130), the velocity terms in the log argument (x m_p/m_c = x_B)."""
    k = 2 * np.pi * xB / Pb
    br = np.sin(om) * (np.cos(u) - e) + np.sqrt(1 - e**2) * np.cos(om) * np.sin(u)
    last = e * np.cos(om) + ((np.cos(u) - e) * np.cos(om) - np.sqrt(1 - e**2) * np.sin(om) * np.sin(u)) / (1 - e * np.cos(u))
    return k / s * e * np.sin(u) - k * s / np.sqrt(1 - e**2) * br * last


# superior conjunction: omega + f = pi/2; find u there
f_sc = np.pi / 2 - om
u_sc = 2 * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(f_sc / 2))
u = u_sc + np.linspace(-0.2, 0.2, 400001)   # +-11 deg around conjunction
L = lam(u, om)
lo = -2 * r * np.log(L)
dS_len = -2 * r * np.log(L + dlen) - lo
dS_ret = -2 * r * np.log(L + dret(u, om)) - lo
dS_both = -2 * r * np.log(L + dlen + dret(u, om)) - lo
f = 2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(u / 2))
psi = om + f
dA = D * np.cos(psi) / L
print(f"Lambda_u min = {L.min():.4e} (1-s = {1-s:.4e}); peak LO delay = {lo.max()*1e6:.3f} us")
print(f"k_ret = 2 pi x_B/P_b = {2*np.pi*xB/Pb:.4e}; dLambda_ret range near conj: {dret(u,om).min():.3e} .. {dret(u,om).max():.3e}")
print(f"lensing:     delay change at conj {dS_len[np.argmin(L)]*1e6:+.4f} us, extreme {dS_len.min()*1e6:+.4f} us")
print(f"retardation: delay change extremes {dS_ret.min()*1e6:+.4f} / {dS_ret.max()*1e6:+.4f} us")
print(f"both:        extremes {dS_both.min()*1e6:+.4f} / {dS_both.max()*1e6:+.4f} us")
print(f"aberration NLO D cos(psi)/Lambda extremes {dA.min()*1e6:+.4f} / {dA.max()*1e6:+.4f} us")

# away from conjunction (orbit average size)
uu = np.linspace(0, 2 * np.pi, 20001)
LL = lam(uu, om)
print(f"lensing far from conj (u_sc+pi): {(-2*r*np.log(LL+dlen)+2*r*np.log(LL))[np.argmax(LL)]*1e9:+.3f} ns")

# circular-orbit identity check: e = 0 => dret = -k s sin(psi) cos(psi)
e0 = e
e = 0.0
uu = np.linspace(0, 2 * np.pi, 7)
chk = dret(uu, om) + 2 * np.pi * xB / Pb * s * np.sin(om + uu) * np.cos(om + uu)
e = e0
print("e=0 identity residual", np.abs(chk).max())

# Independent check of KS99 Eq. (130): build vectors (KS99 Eqs. 120-125, Omega = 0, k = (0,0,-1)),
# move B to where it is when the signal passes it (dt = -k.r/c), and recompute r' + k.r'.
inc = np.arcsin(s)
lv = np.array([1.0, 0, 0]); mv = np.array([0, np.cos(inc), np.sin(inc)]); k = np.array([0, 0, -1.0])
Pv = lv * np.cos(om) + mv * np.sin(om); Qv = -lv * np.sin(om) + mv * np.cos(om)
aR = (x + xB) / s; nb = 2 * np.pi / Pb; mp_over_M = mA / (mA + mB)
worst = 0
for uu in u_sc + np.linspace(-0.1, 0.1, 41):
    ff = 2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(uu / 2))
    rr = aR * (1 - e * np.cos(uu))
    rv = rr * (Pv * np.cos(ff) + Qv * np.sin(ff))
    rdot = nb * aR / np.sqrt(1 - e**2) * (-Pv * np.sin(ff) + Qv * (np.cos(ff) + e))
    vc = -mp_over_M * rdot
    lam_vec = (rr + k @ rv) / aR
    assert abs(lam_vec - lam(uu, om)) < 1e-12, (lam_vec, lam(uu, om))
    rp = rv + vc * (k @ rv)            # r - v_c dt with dt = -k.r
    direct = (np.linalg.norm(rp) + k @ rp) / aR - lam_vec
    worst = max(worst, abs(direct - dret(uu, om)) / abs(dret(uu, om)))
print(f"KS99 Eq.130 vs retarded-position construction: max rel diff {worst:.2e}")
