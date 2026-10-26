"""Checks for ball_lightning_and_clusters.md. Standard library only; run with `python3 -I`.

Each check recomputes a number from its source's own inputs and asserts it against the
published value with a stated tolerance. "derived" checks have no published value and are
asserted against an independent algebraic form of the same quantity.
"""
import math
import sys

E = 1.602176634e-19      # C
ME = 9.1093837015e-31    # kg
C = 299_792_458.0        # m/s
EPS0 = 8.8541878128e-12  # F/m
K = 1 / (4 * math.pi * EPS0)
NA = 6.02214076e23
R_GAS = 8.314462618
SIGMA = 5.670374419e-8
P0 = 101_325.0
MC2_KEV = 510.99895

fails = []


def check(name, got, want, rel):
    ok = abs(got / want - 1) <= rel
    print(f"{'ok  ' if ok else 'FAIL'} {name}: {got:.6g} vs {want:.6g} (tol {rel:.2g})")
    if not ok:
        fails.append(name)


def sphere_volume(d):
    return math.pi * d**3 / 6


# --- Cen, Yuan & Xue 2014, PRL 112, 035001 -----------------------------------------------
pixel = 4e-6 * 902.0 / 3.3e-3                      # pixel pitch x distance / focal length
check("Cen: pixel footprint ~1.1 m", pixel, 1.1, 0.01)
check("Cen: 2D speed 8.6 m/s (9 px over 20-1160 ms)", 9 * pixel / (1.160 - 0.020), 8.6, 0.01)
check("Cen: sound ranging low end 0.898 km", 2.64 * 340.51, 898.0, 0.002)
check("Cen: sound ranging high end 0.9 km", 2.64 * 341.15, 900.0, 0.002)
check("Cen: 99.4 Hz flicker = 2 x 50 Hz mains", 99.4, 2 * 50.0, 0.01)

# --- Virial cap 3 p0 V (Finkelstein & Rubinstein 1964; crates/plasmoid relaxation.rs) ----
def virial_cap(d, p0=P0):
    return 3 * p0 * sphere_volume(d)

check("virial cap, 20 cm ball (repo test value)", virial_cap(0.2), 1.273e3, 1e-3)
check("virial cap, 5 m ball (Cen apparent size; repo report)", virial_cap(5.0), 19.9e6, 0.01)

# --- Abrahamson & Dinniss 2000, Nature 403, 519 ---------------------------------------------
DH_SIO2 = 910.7e3      # J/mol, standard enthalpy of formation of SiO2 (magnitude)
M_SI = 28.0855e-3      # kg/mol
e_chem = 0.100 * DH_SIO2 / M_SI                   # 100 g/m^3 Si loading
check("A&D: chemical energy density 3.3 MJ/m^3", e_chem, 3.3e6, 0.03)
check("A&D: 300 mm ball at 79.0 kW/m^3 -> ~1,100 W", 79.0e3 * sphere_volume(0.3), 1100.0, 0.02)
arr = lambda t1, t2: math.exp(189e3 / R_GAS * (1 / t1 - 1 / t2))
check("A&D: Arrhenius k(1400)/k(1300) vs lifetimes 10 s / 3 s", arr(1300, 1400), 10 / 3, 0.06)
# loose: the 1,200 K lifetime also depends on heat loss, so only order of magnitude.
check("A&D: Arrhenius k(1300)/k(1200) vs lifetimes 30 s / 10 s (loose)", arr(1200, 1300), 30 / 10, 0.5)
check("derived: A&D store / virial cap", e_chem / (3 * P0), e_chem / virial_cap(1.0) * sphere_volume(1.0), 1e-12)

# --- Mesyats 1995, Physics-Uspekhi 38, 567: ectons -------------------------------------------
check("Mesyats: j = (h/t_d)^1/2 = 1e9 A/cm^2 at h=1e9, t_d=1 ns", math.sqrt(1e9 / 1e-9), 1e9, 1e-9)
n_ecton = 1.6 * 10e-9 / E                          # i_th = 1.6 A for 10 ns
check("Mesyats: 1.6 A x 10 ns = 1e11 electrons (portion 1e11-1e12)", n_ecton, 1e11, 0.01)
q_ev = 1e11 * E                                    # Shoulders' typical EV charge
m_eroded = 0.6e-4 * q_ev                           # g, gamma_m = 0.6e-4 g/C
check("Mesyats: energy to explode the eroded mass at 1e4 J/g ~ 1e-8 J", 1e4 * m_eroded, 1e-8, 0.05)

# --- Shoulders' EV (US 5,018,180): 1e11 electrons in ~1 um, ions <= 1e-5 per electron -------
r = 0.5e-6
u_sphere = 0.6 * K * q_ev**2 / r                   # uniform sphere
check("derived: EV Coulomb self-energy 3kQ^2/5R = (3/5) Q V_surface", u_sphere, 0.6 * q_ev * (K * q_ev / r), 1e-12)
check("derived: EV self-energy ~2.77 J", u_sphere, 2.77, 0.01)
check("derived: EV surface potential ~288 MV", K * q_ev / r, 2.88e8, 0.01)
check("derived: EV surface field ~5.8e14 V/m", K * q_ev / r**2, 5.76e14, 0.01)
n_e = 1e11 / sphere_volume(1e-6)
n_cu = 8960 / 63.546e-3 * NA                       # one conduction electron per Cu atom
check("derived: EV electron density / Cu conduction density ~2.25", n_e / n_cu, 2.25, 0.01)
wp = math.sqrt(n_e * E**2 / (EPS0 * ME))
check("derived: explosion time sqrt(3)/w_p = sqrt(m R^3 / kQe)", math.sqrt(3) / wp, math.sqrt(ME * r**3 / (K * q_ev * E)), 1e-9)
check("derived: explosion time ~7e-17 s", math.sqrt(3) / wp, 7.04e-17, 0.01)
gamma = 1 / math.sqrt(1 - 0.1**2)
check("derived: KE per electron at 0.1c ~2.57 keV", (gamma - 1) * MC2_KEV, 2.57, 0.01)
check("derived: EV kinetic energy at 0.1c ~41 uJ", 1e11 * (gamma - 1) * MC2_KEV * 1e3 * E, 4.12e-5, 0.01)
check("derived: Budker 1/gamma^2 at 0.1c = 1 - beta^2", 1 / gamma**2, 0.99, 1e-12)
check("derived: Budker shortfall vs 1e-5 ions per electron", (1 / gamma**2) / 1e-5, 9.9e4, 1e-9)
b_brill = math.sqrt(2 * ME * n_e / EPS0)
check("derived: Brillouin field to hold n_e, eps0 B^2/2m = n", EPS0 * b_brill**2 / (2 * ME), n_e, 1e-12)
check("derived: Brillouin field ~2.0e5 T", b_brill, 1.98e5, 0.01)
check("derived: EV self-energy / ecton energy ~2.9e8", u_sphere / (1e4 * m_eroded), 2.9e8, 0.02)

# --- Iron: Wikipedia infobox values; "FeO2" ------------------------------------------------
M_FE, M_O = 55.845, 15.999
L_FE = 13.81e3 / (M_FE * 1e-3)                     # J/kg
check("derived: Fe latent heat 13.81 kJ/mol = 247 kJ/kg", L_FE, 2.473e5, 0.002)
fe_frac = {f: M_FE * n / (M_FE * n + M_O * o) for f, (n, o) in
           {"FeO": (1, 1), "Fe3O4": (3, 4), "Fe2O3": (2, 3), "FeO2": (1, 2)}.items()}
check("derived: Fe mass fraction in FeO 77.7%", fe_frac["FeO"], 0.777, 0.002)
check("derived: Fe mass fraction in Fe3O4 72.4%", fe_frac["Fe3O4"], 0.724, 0.002)
check("derived: Fe mass fraction in Fe2O3 69.9%", fe_frac["Fe2O3"], 0.699, 0.002)
check("derived: Fe mass fraction in FeO2 63.6%", fe_frac["FeO2"], 0.636, 0.002)
x_max = (0.07 / M_O) / (0.93 / M_FE)               # Adamenko 93 wt% Fe, rest all O at most
check("derived: 93 wt% Fe allows at most FeO_0.26", x_max, 0.263, 0.01)

RHO_FE = 6980.0
def solidify_radiative(d, t=1811.0, emissivity=1.0):
    """Lower bound on latent-heat removal time (s) by radiation alone, lumped droplet."""
    return RHO_FE * L_FE * d / (6 * emissivity * SIGMA * t**4)

def solidify_conductive(d, k_gas=0.067, dt=1500.0):
    """Latent-heat removal time (s) by still-gas conduction, Nu = 2 (k_gas assumed)."""
    return RHO_FE * L_FE * d**2 / (12 * k_gas * dt)

check("derived: 50 um Fe drop, radiation-only solidification ~24 ms", solidify_radiative(50e-6), 0.0236, 0.01)
check("derived: 50 um Fe drop, conduction Nu=2 ~3.6 ms", solidify_conductive(50e-6), 3.58e-3, 0.01)
check("derived: 10 um drop conduction scales as d^2",
      solidify_conductive(10e-6) / solidify_conductive(100e-6), 0.01, 1e-12)

# --- Wendt & Irion: energy per W atom at Mesyats's 1e4 J/g electrical-explosion density -------
ev_per_atom = 1e4 * 183.84 / 96485.33212           # J/g x g/mol / (J/mol per eV)
check("derived: 1e4 J/g in tungsten ~19 eV per atom", ev_per_atom, 19.05, 0.01)

if fails:
    print(f"\n{len(fails)} check(s) failed: {fails}")
    sys.exit(1)
print("\nall checks passed")
