"""DLMF 19.36(i): E(phi|m) by 19.25.9 loses digits for m near 1 and phi > pi/4 (Reinsch & Raab 2000);
19.25.10 does not. Reference: E = int_0^phi sqrt(1 - m sin^2) by adaptive quadrature in extended precision
is not available, so use scipy.special.ellipeinc (Cephes) and the m = 1 limit sin(phi) as checks."""
import contextlib, io, math
with contextlib.redirect_stdout(io.StringIO()):
    from carlson_checks import E, E_naive
from scipy import special, integrate
for m in [0.99, 1 - 1e-6, 1 - 1e-10, 1 - 1e-14]:
    for phi in [1.0, 1.5, 1.5707]:
        q = integrate.quad(lambda t: math.sqrt(1 - m * math.sin(t) ** 2), 0, phi, epsabs=0, epsrel=1e-13, limit=500)[0]
        print(f"m = 1-{1 - m:.0e}, phi = {phi}: 19.25.10 rel {E(phi, m) / q - 1:+.1e}  19.25.9 rel {E_naive(phi, m) / q - 1:+.1e}"
              f"  Cephes rel {special.ellipeinc(phi, m) / q - 1:+.1e}")
