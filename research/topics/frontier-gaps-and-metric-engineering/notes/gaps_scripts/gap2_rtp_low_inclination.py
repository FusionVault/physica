"""Gap 2b: RTP of one induced dipole at low inclination, every method on the same grid.

48 x 48 grid at 25 m (1200 m square, near the crate test grid), dipole 100 m below the centre,
induced (m = f), declination 0 unless stated. Truth: the same dipole at the pole. FFT methods use
25 % edge padding (Harmonica's default). Errors are rms(estimate - truth)/rms(truth) over the grid.
Noise: Gaussian, 1 % of the anomaly's peak, seed 7. 'floor' is rms(noise)/rms(truth), what a
perfect noise-passing filter would leave. Damping parameters are oracle-chosen (best of a scan):
an upper bound on what a practical choice (GCV, L-curve, discrepancy) achieves.
"""
import numpy as np
from gap2_rtp_padding import unit, dipole, rms_rel, theta, rtp_pseudo, rtp_wiener, magmap_ia

n, h, depth = 48, 25.0, 100.0
xs, ys = np.meshgrid(np.arange(n) * h, np.arange(n) * h, indexing="ij")
c0 = (n // 2) * h
pole = np.array([0.0, 0.0, 1.0])
truth = dipole(depth, pole, pole, c0, c0, xs, ys)
P = n // 4

def fft_apply(v, kern):
    vp = np.pad(v, P, mode="edge")
    m = vp.shape[0]
    kx = 2 * np.pi * np.fft.fftfreq(m, h)[:, None] * np.ones((1, m))
    ky = 2 * np.pi * np.fft.fftfreq(m, h)[None, :] * np.ones((m, 1))
    return np.real(np.fft.ifft2(np.fft.fft2(vp) * kern(kx, ky)))[P:P + n, P:P + n]

def plain(f):
    def k(kx, ky):
        t = theta(f, kx, ky) ** 2; out = 1 / t; out[0, 0] = 0; return out
    return k

px, py = xs.ravel(), ys.ravel()
DXM, DYM = px[:, None] - px[None, :], py[:, None] - py[None, :]
def eqs_factory(f, layer):
    """Damped equivalent-dipole layer: one eigendecomposition of G^T G, then any damping."""
    G = dipole(layer, f, f, 0.0, 0.0, DXM, DYM); Gp = dipole(layer, pole, pole, 0.0, 0.0, DXM, DYM)
    e, V = np.linalg.eigh(G.T @ G)
    s2 = e.max()
    def solve(data, damp):
        q = V @ ((V.T @ (G.T @ data.ravel())) / (e + damp ** 2 * s2))
        return (Gp @ q).reshape(n, n)
    return solve

rng = np.random.default_rng(7)
noise_unit = rng.standard_normal((n, n))
for inc, dec in [(45, 0), (20, 0), (10, 0), (10, 30), (5, 0), (2, 0)]:
    f = unit(inc, dec)
    clean = dipole(depth, f, f, c0, c0, xs, ys)
    noise = 0.01 * np.abs(clean).max() * noise_unit
    print(f"I={inc:>2} D={dec:>2}  floor {np.sqrt(np.sum(noise**2)/np.sum(truth**2)):.3f}")
    eqs = eqs_factory(f, 2.5 * h)
    for label, d in (("clean", clean), ("noisy", clean + noise)):
        r = {"plain": rms_rel(fft_apply(d, plain(f)), truth)}
        for ia in sorted({magmap_ia(inc), min(90, inc + 10), min(90, inc + 20)}):
            r[f"pseudoI Ia={ia:g}"] = rms_rel(fft_apply(d, rtp_pseudo(inc, dec, ia)), truth)
        r["wiener*"] = min(rms_rel(fft_apply(d, rtp_wiener(f, f, e)), truth) for e in 10.0 ** np.arange(-8, 0.5, 0.25))
        r["eqs*"] = min(rms_rel(eqs(d, dmp), truth) for dmp in 10.0 ** np.arange(-6, -0.5, 0.25))
        print(f"   {label}: " + ", ".join(f"{k} {v:.3g}" for k, v in r.items()))
