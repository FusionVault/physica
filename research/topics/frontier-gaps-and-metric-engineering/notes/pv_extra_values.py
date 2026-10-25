"""Extra test values for pv_metric_tests.md (run: python3 -P pv_extra_values.py)."""
import numpy as np
src = open(__file__.replace("pv_extra_values.py", "pv_metric_tests.py")).read()
exec(src[:src.index('say("\\n== 2. Light bending")')])   # imports, helpers, deflection()
for b in (6.0, 7.0, 10.0, 20.0, 100.0):
    print(f"b = {b:6.1f} M: α_GR = {deflection(b,'GR'):.12f}  α_PV = {deflection(b,'PV'):.12f}")
for eps in (1e-3, 1e-6):
    print(f"b = b_c(1+{eps:g}): α_GR = {deflection(np.sqrt(27)*(1+eps),'GR'):.10f}, α_PV = {deflection(2*np.e*(1+eps),'PV'):.10f}")
# Mercury: near-circular extra precession 2π(K_PV − K_GR) ≈ 2π x^2, x = GM/(c^2 a)
x = 1476.625 / 5.7909e10
per_orbit = 2 * np.pi * x * x
print(f"Mercury x = {x:.4e}: extra per orbit {per_orbit:.3e} rad; per century {per_orbit*100/0.240846*206264.806e6:.3f} µas")
# Vagnozzi et al. 2023 combined Sgr A* δ, 1σ interval −0.125..0.005 -> Gaussian mean −0.06, σ 0.065
d_pv = 2 * np.e / np.sqrt(27) - 1
z_pv, z_gr = (d_pv + 0.06) / 0.065, 0.06 / 0.065
print(f"combined Sgr A*: PV {z_pv:.2f}σ, Schwarzschild {z_gr:.2f}σ, LR GR:PV = {np.exp((z_pv**2 - z_gr**2)/2):.2f}; PV r_sh = {2*np.e:.4f} vs 1σ upper 5.22, 2σ upper 5.56")
