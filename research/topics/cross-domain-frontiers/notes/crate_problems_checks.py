# Checks of four problems in the published crates that the report adds to the notes' findings.
# Standard library only: python3 -I crate_problems_checks.py. Constants as in crates/metrology.
import math

C, G, GM_SUN = 299_792_458.0, 6.674_30e-11, 1.327_124_4e20
M_E, M_P, MU_0 = 9.109_383_701_5e-31, 1.672_621_923_69e-27, 1.256_637_062_12e-6
R_SUN = 6.957e8  # IAU 2015 nominal
ARCSEC = math.pi / 180 / 3600


def check(label, got, ok):
    print(f"{'ok ' if ok else 'BAD'} {label}: {got:.4g}")
    assert ok, label


# 1. gravitics::electrodes::ring_field (Zypman form, crate's AGM elliptic) near the axis, against
#    the near-axis series E_rho = -(rho/2) dE_z/dz (div E = 0). Ring and offset from the crate test.
def elliptic(m):
    a, b, c = 1.0, math.sqrt(1 - m), math.sqrt(m)
    s, w = 0.5 * m, 0.5
    for _ in range(64):
        if c <= 2.220446049250313e-16 * a:
            break
        a, b, c = 0.5 * (a + b), math.sqrt(a * b), 0.5 * (a - b)
        w *= 2.0
        s += w * c * c
    k = math.pi / (2 * a)
    return k, k * (1 - s)


def ring_e_rho(r, rho, dz):  # charge q = 1 with 1/(4 pi eps0) = 1, as in the crate up to K
    s, d = (r + rho) ** 2 + dz * dz, (r - rho) ** 2 + dz * dz
    k, e = elliptic(4 * r * rho / s)
    return 1 / (math.pi * rho * math.sqrt(s)) * (k - (r * r - rho * rho + dz * dz) / d * e)


def series_e_rho(r, rho, dz):
    return -(rho / 2) * (r * r - 2 * dz * dz) / (r * r + dz * dz) ** 2.5


r, dz = 0.04, 0.01
err = {f: abs(ring_e_rho(r, f * r, dz) / series_e_rho(r, f * r, dz) - 1) for f in (1e-3, 1e-5, 1e-7)}
check("ring E_rho at rho = 1e-3 a, relative error (series truncation only)", err[1e-3], err[1e-3] < 1e-5)
check("ring E_rho at rho = 1e-5 a, relative error (cancellation)", err[1e-5], 1e-7 < err[1e-5] < 1e-4)
check("ring E_rho at rho = 1e-7 a, relative error (cancellation)", err[1e-7], err[1e-7] > 1e-2)

# 2. gravitational::hf::plasma_index doc: "Gordon optical metric ... n = 1/sqrt(g00)". In isotropic
#    Schwarzschild coordinates the optical index is sqrt(B/A) = (1 + rs/4r)^3/(1 - rs/4r), so
#    n - 1 -> 2GM/(rc^2); 1/sqrt(g00) = (1 + rs/4r)/(1 - rs/4r) gives n - 1 -> GM/(rc^2). A ray with
#    n - 1 = k/r is deflected by 2k/b, so the doc's index bends light by half the GR value.
full = 2 * (2 * GM_SUN / C**2) / R_SUN / ARCSEC
half = 2 * (GM_SUN / C**2) / R_SUN / ARCSEC
check("limb deflection from n = sqrt(B/A), arcsec", full, abs(full - 1.7512) < 1e-3)
check("limb deflection from the doc's n = 1/sqrt(g00), arcsec", half, abs(half - 0.8756) < 1e-3)

# 3. Spitzer coefficient: plasmoid::pinch uses sigma0 = 9.71e3 (eta = 1/9.71e3), relaxation::lundquist
#    uses 1.03e-4; the Braginskii derivation in plasma_frontiers_checks.py gives 1.0314e-4.
pinch, relax, derived = 1 / 9.71e3, 1.03e-4, 1.0314e-4
check("pinch vs lundquist Spitzer coefficient, relative gap", relax / pinch - 1, abs(relax / pinch - 1) < 2e-4)
check("derived vs pinch Spitzer coefficient, relative gap", derived / pinch - 1, 1e-3 < derived / pinch - 1 < 2e-3)

# 4. Alfven speed: relaxation::lundquist uses rho = n A m_p; magnetics::mhd::alfven_speed with
#    mass_density adds the electrons (n (m_p + m_e) for hydrogen).
gap = math.sqrt(1 + M_E / M_P) - 1
check("lundquist v_A over mhd::alfven_speed - 1, hydrogen", gap, abs(gap - 2.72e-4) < 1e-6)

print("all checks passed")
