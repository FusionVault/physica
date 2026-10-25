"""The uniform-disc magnification with the lens on the limb (u = rho): Lee et al. 2009's exact integral
against the closed form (1/pi)[2/rho + (1+rho^2)/rho^2 (pi/2 + arcsin((rho^2-1)/(rho^2+1)))]
(recalled as Witt & Mao 1994's u = rho case; that paper was not opened)."""
import math
from finite_source import A_uni

for rho in [0.001, 0.01, 0.05, 0.5, 1.0, 3.0]:
    closed = (2 / rho + (1 + rho * rho) / rho ** 2 * (math.pi / 2 + math.asin((rho * rho - 1) / (rho * rho + 1)))) / math.pi
    ex = A_uni(rho, rho)
    print(f"rho={rho}: Lee integral {ex:.15f}  closed form {closed:.15f}  rel {ex/closed-1:+.1e}")
