"""Hogg 1999 Eq. 17 (Mattig, Omega_Lambda = 0) as an exact test of the curved D_M, and a case where
Hogg's D_A12 formula (valid only for Omega_k >= 0) fails while S_k(chi2 - chi1)/(1+z2) holds."""
import math
from cosmology_checks import Cosmo

for om in [0.3, 1.0, 2.0]:
    c = Cosmo(70, om, 0.0)
    for z in [0.5, 2.0, 5.0]:
        mattig = c.DH() * 2 * (2 - om * (1 - z) - (2 - om) * math.sqrt(1 + om * z)) / (om * om * (1 + z))
        print(f"Om={om} OL=0 Ok={c.Ok0:+.1f} z={z}: integral D_M={c.DM(z):.10f}  Mattig={mattig:.10f}  rel={c.DM(z)/mattig-1:+.1e}")
# closed universe where sqrt|Ok| chi exceeds pi/2: Hogg's D_A12 picks the wrong branch
c = Cosmo(70, 3.0, 0.0)
for z1, z2 in [(0.5, 1.0), (1.0, 10.0), (2.0, 50.0)]:
    x = math.sqrt(abs(c.Ok0)) * c.chi(0, z2)
    print(f"Om=3 Ok={c.Ok0:+.0f}: z1={z1} z2={z2}: sqrt|Ok| chi2 = {x:.3f} rad; S_k form D_A12 = {c.DA(z1, z2):.4f} Mpc; Hogg Eq.19 = {c.DA_hogg(z1, z2):.4f} Mpc")
