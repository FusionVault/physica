"""Window-stability check of the time-domain test-field QNMs in pv_metric_tests.py (run: python3 -P pv_qnm_check.py)."""
import numpy as np
from scipy.optimize import brentq, curve_fit
src = open(__file__.replace("pv_qnm_check.py", "pv_metric_tests.py")).read()
exec(src[src.index("def tortoise_grid"):src.index("for metric, s, l in")])
for metric, s, l in (("GR", 2, 2), ("PV", 2, 2), ("GR", 1, 1), ("PV", 1, 1)):
    t, y, V, r, xs_ = ring(metric, s, l, L=900.0, h=0.02, T=330.0, obs=30.0)
    fits = [fit_mode(t, y, a, a + 70) for a in (60, 70, 80, 90, 100)]
    print(metric, s, l, " ".join(f"[{a}]{wr:.4f}-{wi:.4f}i" for a, (wr, wi) in zip((60, 70, 80, 90, 100), fits)))
    # late-time envelope: |psi| at t = 150, 200, 250, 300
    print("   |psi| at t=150,200,250,300:", " ".join(f"{abs(y[np.argmin(abs(t-T0))]):.2e}" for T0 in (150, 200, 250, 300)))
