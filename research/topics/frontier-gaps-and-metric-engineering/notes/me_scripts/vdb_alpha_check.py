"""Van Den Broeck 1999: which alpha reproduces the paper's T00 peak (w=0.349, -4.9e2/Dt^2) and r_c = Dt/72.5?"""
import numpy as np
for alpha in (1e17, 1e34):
    for n in (60, 80, 100):
        w = np.linspace(1e-6, 1 - 1e-9, 2_000_001)
        g = n * w**(n - 1) - (n - 1) * w**n
        g1 = n * (n - 1) * (w**(n - 2) - w**(n - 1))
        g2 = n * (n - 1) * ((n - 2) * w**(n - 3) - (n - 1) * w**(n - 2))
        B, Br, Brr, r = alpha * g + 1, -alpha * g1, alpha * g2, 2 - w
        T = (Br**2 / B**4 - 2 * Brr / B**3 - 4 * Br / (B**3 * r)) / (8 * np.pi)
        R1212 = Br**2 / B**4 - Brr / B**3 - Br / (B**3 * r)
        i, j = np.argmin(T), np.argmax(abs(R1212))
        dens = 4 * np.pi * r**2 * B**3 * T
        Em, Ep = np.trapezoid(np.where(T < 0, dens, 0), w), np.trapezoid(np.where(T < 0, 0, dens), w)
        print(f"alpha={alpha:.0e} n={n}: T00 min {T[i]:.4g}/Dt^2 at w={w[i]:.4f}; max|R1212| at w={w[j]:.4f}, r_c=Dt/{np.sqrt(abs(R1212[j])):.2f}; "
              f"E-/(alpha Dt)={Em/alpha:.4g}, E+/(alpha Dt)={Ep/alpha:.4g}")

# QI test (VdB Eq. QInumbers) with the parameters the paper lists (alpha=1e17, Dt=1e-15 m) versus alpha=1e34, Dt=1e-32 m
hbar, c, G = 1.054571817e-34, 299792458.0, 6.67430e-11
for Dt, Tmin, rc_div in ((1e-15, 175.9, 43.39), (1e-32, 492.0, 72.46)):
    tau0 = 0.1 * Dt / rc_div
    lhs = Tmin / Dt**2 * c**2 / G
    rhs = 3 * hbar / (32 * np.pi**2 * c * tau0**4)
    print(f"Dt={Dt:.0e} m: |T00| = {lhs:.3g} kg/m^3, QI allows {rhs:.3g} kg/m^3 -> ratio {lhs/rhs:.3g} ({'violated' if lhs > rhs else 'satisfied'})")
