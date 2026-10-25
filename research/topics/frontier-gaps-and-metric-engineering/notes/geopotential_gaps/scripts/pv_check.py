"""Independent principal-value check of Pi(1.5; 1.4 | 0.9) by singularity subtraction (no QAWC)."""
import math
from scipy import integrate
n, phi, m = 1.5, 1.4, 0.9
g = lambda t: 1 / ((1 - n * math.sin(t) ** 2) * math.sqrt(1 - m * math.sin(t) ** 2))
t0 = math.asin(1 / math.sqrt(n))
c = -1 / (n * math.sin(2 * t0) * math.sqrt(1 - m / n))  # residue-like coefficient: g ~ c/(t - t0)
reg = lambda t: g(t) - c / (t - t0) if abs(t - t0) > 1e-6 else 0.0
val = integrate.quad(reg, 0, phi, points=[t0], epsabs=0, epsrel=1e-12, limit=1000)[0] + c * math.log((phi - t0) / t0)
print("Pi(1.5;1.4|0.9) by subtraction:", val)
