"""Gap 2: FFT padding and stabilised reduction to the pole, on synthetic dipoles.

Mirrors geopotential::planar (128 x 128 grid at 10 m, x north = axis 0, y east, z down) and first
reproduces its pinned values (continuation 6.9e-4 rms, RTP 1.1 % at I = 45 deg), then measures
(a) padding modes and widths, (b) RTP at low inclination: plain, pseudo-inclination (Seequent /
Geosoft MAGMAP formula), Wiener/Tikhonov damping, a damped equivalent-dipole layer, and the
reduction to the equator. All numbers printed here are this script's own synthetic results.
"""
import numpy as np

NX = NY = 128
DX = 10.0
X, Y = np.meshgrid(np.arange(NX) * DX, np.arange(NY) * DX, indexing="ij")
KX = 2 * np.pi * np.fft.fftfreq(NX, DX)[:, None] * np.ones((1, NY))
KY = 2 * np.pi * np.fft.fftfreq(NY, DX)[None, :] * np.ones((NX, 1))

def unit(inc, dec):
    i, d = np.radians(inc), np.radians(dec)
    return np.array([np.cos(i) * np.cos(d), np.cos(i) * np.sin(d), np.sin(i)])

def dipole(depth, m, f, x0=640.0, y0=640.0, X=X, Y=Y, z=0.0):
    d = [X - x0, Y - y0, np.full_like(X, z - depth)]
    r2 = d[0] ** 2 + d[1] ** 2 + d[2] ** 2
    dm = sum(d[i] * m[i] for i in range(3)); df = sum(d[i] * f[i] for i in range(3))
    return (3 * dm * df - r2 * np.dot(m, f)) / r2 ** 2.5

def rms_rel(a, b):
    return np.sqrt(np.sum((a - b) ** 2) / np.sum(b ** 2))

def theta(u, kx, ky):
    k = np.hypot(kx, ky); k = np.where(k == 0, 1, k)
    return u[2] + 1j * (u[0] * kx + u[1] * ky) / k

def apply(values, kernel_fn, pad=0, mode="none", taper=False):
    """Pad (numpy.pad modes), optionally cosine-taper the pad, filter, crop."""
    v = values
    if pad:
        v = np.pad(values, pad, mode="constant" if mode == "zero" else mode)
        if taper:  # Tukey-like cosine roll-off over the padded margin
            n0, n1 = v.shape
            w = lambda n: np.concatenate([0.5 - 0.5 * np.cos(np.pi * np.arange(pad) / pad), np.ones(n - 2 * pad), 0.5 + 0.5 * np.cos(np.pi * np.arange(1, pad + 1) / pad)])
            v = v * np.outer(w(n0), w(n1))
    n0, n1 = v.shape
    kx = 2 * np.pi * np.fft.fftfreq(n0, DX)[:, None] * np.ones((1, n1))
    ky = 2 * np.pi * np.fft.fftfreq(n1, DX)[None, :] * np.ones((n0, 1))
    out = np.real(np.fft.ifft2(np.fft.fft2(v) * kernel_fn(kx, ky)))
    return out[pad:pad + values.shape[0], pad:pad + values.shape[1]] if pad else out

def up(h):
    return lambda kx, ky: np.exp(-np.hypot(kx, ky) * h)

def rtp(m, f):
    def k(kx, ky):
        t = theta(m, kx, ky) * theta(f, kx, ky)
        out = 1 / t; out[0, 0] = 0
        return out
    return k

def magmap_ia(inc):  # Seequent Oasis montaj REDP help, Ia auto rule (Equation199)
    ia = min(90.0, max(0.0, int((0.008 * inc ** 2 - 1.71 * abs(inc) + 80) / 2.5) * 2.5))
    return inc if abs(ia) < abs(inc) else ia

def rtp_pseudo(inc, dec, ia):  # Seequent REDP L(theta), induced magnetisation (Equation198)
    def k(kx, ky):
        th = np.arctan2(ky, kx)  # wavenumber azimuth from north (x) towards east (y)
        c = np.cos(np.radians(dec) - th)
        I, A = np.radians(inc), np.radians(ia)
        num = (np.sin(I) - 1j * np.cos(I) * c) ** 2
        den = (np.sin(A) ** 2 + np.cos(A) ** 2 * c ** 2) * (np.sin(I) ** 2 + np.cos(I) ** 2 * c ** 2)
        out = num / den; out[0, 0] = 0
        return out
    return k

def rtp_wiener(m, f, eps):
    def k(kx, ky):
        t = theta(m, kx, ky) * theta(f, kx, ky)
        out = np.conj(t) / (np.abs(t) ** 2 + eps); out[0, 0] = 0
        return out
    return k

def rte_pseudo(inc, dec, ia):  # Seequent REDE L(theta) (Equation193)
    def k(kx, ky):
        th = np.arctan2(ky, kx)
        c = np.cos(np.radians(dec) - th)
        I, A = np.radians(inc), np.radians(ia)
        num = (np.sin(I) - 1j * np.cos(I) * c) ** 2 * (-c ** 2)
        den = (np.sin(A) ** 2 + np.cos(A) ** 2 * c ** 2) * (np.sin(I) ** 2 + np.cos(I) ** 2 * c ** 2)
        out = num / den; out[0, 0] = 0
        return out
    return k


if __name__ == "__main__":
    # ---- 1. reproduce the crate's pinned values -------------------------------------------------
    pole = np.array([0.0, 0.0, 1.0])
    e_up = rms_rel(apply(dipole(50, pole, pole), up(30)), dipole(80, pole, pole))
    t45 = unit(45, 0)
    e45 = rms_rel(apply(dipole(50, t45, t45), rtp(t45, t45)), dipole(50, pole, pole))
    print(f"crate check: upward 30 m {e_up:.3e} (crate 6.9e-4); RTP I=45 {e45:.4f} (crate 0.011)")

    # ---- 2. padding ------------------------------------------------------------------------------
    print("\npadding: rms relative error on the original 128x128 grid")
    cases = {
        "up 30 m, dipole 50 m deep (decays in grid)": (dipole(50, pole, pole), up(30), dipole(80, pole, pole)),
        "up 30 m, dipole 300 m deep (edges carry signal)": (dipole(300, pole, pole), up(30), dipole(330, pole, pole)),
        "up 30 m, dipole 100 m deep at 1/8 from an edge": (dipole(100, pole, pole, x0=160.0), up(30), dipole(130, pole, pole, x0=160.0)),
        "RTP I=45, dipole 300 m deep": (dipole(300, t45, t45), rtp(t45, t45), dipole(300, pole, pole)),
    }
    for name, (field, kern, truth) in cases.items():
        row = [f"none {rms_rel(apply(field, kern), truth):.2e}"]
        for frac in (0.25, 0.5, 1.0):
            p = int(frac * NX)
            for mode in ("zero", "edge", "symmetric", "reflect"):
                row.append(f"{mode}{int(frac*100)} {rms_rel(apply(field, kern, p, mode), truth):.2e}")
            row.append(f"edge{int(frac*100)}+taper {rms_rel(apply(field, kern, p, 'edge', True), truth):.2e}")
            row.append(f"symm{int(frac*100)}+taper {rms_rel(apply(field, kern, p, 'symmetric', True), truth):.2e}")
        print(" ", name); print("    " + "; ".join(row))

    # ---- 3. reduction to the pole at low inclination ---------------------------------------------
    rng = np.random.default_rng(1)
    print("\nRTP of a dipole 50 m deep, induced, rms relative error against the true pole anomaly")
    for inc, dec in [(45, 0), (20, 0), (10, 0), (10, 30), (5, 0)]:
        f = unit(inc, dec); data = dipole(50, f, f); truth = dipole(50, pole, pole)
        noise = 0.01 * np.abs(data).max() * rng.standard_normal(data.shape)
        res = []
        for label, d in (("clean", data), ("1% noise", data + noise)):
            plain = rms_rel(apply(d, rtp(f, f)), truth)
            ia = magmap_ia(inc)
            pseudo = rms_rel(apply(d, rtp_pseudo(inc, dec, ia)), truth)
            best = min((rms_rel(apply(d, rtp_wiener(f, f, e)), truth), e) for e in 10.0 ** np.arange(-8, 0.5, 0.25))
            res.append(f"{label}: plain {plain:.3g}, pseudo-I(Ia={ia:g}) {pseudo:.3g}, Wiener best {best[0]:.3g} (eps={best[1]:.1e})")
        print(f"  I={inc:>2} D={dec:>2}: " + " | ".join(res))

    # The equivalent-dipole layer comparison is in gap2_rtp_low_inclination.py (same grid for all methods).

    # Reduction to the equator: compare with the anomaly the same dipole would have at I = 0.
    print("\nRTE (Seequent REDE, Ia = 20 deg default) against the true equator anomaly")
    for inc in (10, 20, 45):
        f = unit(inc, 0); data = dipole(50, f, f)
        eq = unit(0, 0); truth = dipole(50, eq, eq)
        rte = apply(data, rte_pseudo(inc, 0, max(20.0, inc)))
        exact = apply(data, lambda kx, ky: (lambda t: np.where(np.hypot(kx, ky) == 0, 0, t))(theta(eq, kx, ky) ** 2 / (theta(f, kx, ky) ** 2)))
        print(f"  I={inc}: RTE(Ia) {rms_rel(rte, truth):.3g}; exact ratio Theta_eq^2/Theta_f^2 {rms_rel(exact, truth):.3g}")
