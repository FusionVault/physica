"""Gap 2c: reduction to the equator with 1 % noise, I = 10 and 5 deg (grid and noise as gap2b)."""
import numpy as np
from gap2_rtp_padding import unit, dipole, rms_rel, theta, rte_pseudo
n, h, depth = 48, 25.0, 100.0
xs, ys = np.meshgrid(np.arange(n) * h, np.arange(n) * h, indexing="ij"); c0 = (n // 2) * h
P = n // 4
def fft_apply(v, kern):
    vp = np.pad(v, P, mode="edge"); m = vp.shape[0]
    kx = 2 * np.pi * np.fft.fftfreq(m, h)[:, None] * np.ones((1, m)); ky = 2 * np.pi * np.fft.fftfreq(m, h)[None, :] * np.ones((m, 1))
    return np.real(np.fft.ifft2(np.fft.fft2(vp) * kern(kx, ky)))[P:P + n, P:P + n]
noise_unit = np.random.default_rng(7).standard_normal((n, n))
eq = unit(0, 0); truth = dipole(depth, eq, eq, c0, c0, xs, ys)
for inc in (10, 5):
    f = unit(inc, 0); clean = dipole(depth, f, f, c0, c0, xs, ys)
    noise = 0.01 * np.abs(clean).max() * noise_unit
    exact = lambda kx, ky: np.where(np.hypot(kx, ky) == 0, 0, theta(eq, kx, ky) ** 2 / theta(f, kx, ky) ** 2)
    for lab, d in (("clean", clean), ("noisy", clean + noise)):
        print(f"I={inc} {lab}: floor {np.sqrt(np.sum(noise**2)/np.sum(truth**2)):.3f}; exact RTE {rms_rel(fft_apply(d, exact), truth):.3g}; Seequent REDE Ia=20 {rms_rel(fft_apply(d, rte_pseudo(inc, 0, 20.0)), truth):.3g}")
