//! Time series with a given noise spectrum, and whitening of series by that spectrum: the two
//! frequency-domain steps of a closed-loop simulation.

use core::f64::consts::TAU;

/// In-place radix-2 FFT of `(re, im)`; `inverse` uses `e^{+i…}` and divides by the length.
///
/// # Panics
/// If the length is not a power of two.
pub fn fft(re: &mut [f64], im: &mut [f64], inverse: bool) {
    let n = re.len();
    assert!(n.is_power_of_two() && im.len() == n, "FFT length must be a power of two");
    let mut j = 0;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j |= bit;
        if i < j {
            re.swap(i, j);
            im.swap(i, j);
        }
    }
    let sign = if inverse { 1.0 } else { -1.0 };
    let mut len = 2;
    while len <= n {
        let (wr, wi) = ((TAU / len as f64).cos(), sign * (TAU / len as f64).sin());
        for start in (0..n).step_by(len) {
            let (mut cr, mut ci) = (1.0, 0.0);
            for k in 0..len / 2 {
                let (a, b) = (start + k, start + k + len / 2);
                let (tr, ti) = (re[b] * cr - im[b] * ci, re[b] * ci + im[b] * cr);
                (re[b], im[b]) = (re[a] - tr, im[a] - ti);
                (re[a], im[a]) = (re[a] + tr, im[a] + ti);
                (cr, ci) = (cr * wr - ci * wi, cr * wi + ci * wr);
            }
        }
        len <<= 1;
    }
    if inverse {
        let k = 1.0 / n as f64;
        re.iter_mut().chain(im.iter_mut()).for_each(|x| *x *= k);
    }
}

/// Multiply the spectrum of `x` (sampled every `dt` s) by `gain(f)`, zero-padded to a power of
/// two; the DC term is dropped. Shared by noise colouring and whitening.
fn filter(x: &mut [f64], dt: f64, gain: impl Fn(f64) -> f64) {
    let n = x.len().next_power_of_two();
    let (mut re, mut im) = (vec![0.0; n], vec![0.0; n]);
    re[..x.len()].copy_from_slice(x);
    fft(&mut re, &mut im, false);
    for k in 0..n {
        let f = k.min(n - k) as f64 / (n as f64 * dt);
        let g = if k == 0 { 0.0 } else { gain(f) };
        re[k] *= g;
        im[k] *= g;
    }
    fft(&mut re, &mut im, true);
    x.copy_from_slice(&re[..x.len()]);
}

/// `n` samples, every `dt` seconds, of Gaussian noise whose one-sided amplitude spectral density
/// is `asd(f)` (units per √Hz). Deterministic for a given `seed`.
pub fn colored_noise(asd: impl Fn(f64) -> f64, dt: f64, n: usize, seed: u64) -> Vec<f64> {
    let mut x = gaussian(n, seed);
    // Unit-variance white noise has one-sided PSD 2·dt; rescale it to asd².
    let k = 1.0 / (2.0 * dt).sqrt();
    filter(&mut x, dt, |f| asd(f) * k);
    x
}

/// Divide the spectrum of `x` by `asd(f)`: noise with that ASD becomes white with one-sided PSD 1,
/// i.e. per-sample variance `1/(2·dt)`. Apply it to observations and to every design-matrix column.
pub fn whiten(x: &mut [f64], dt: f64, asd: impl Fn(f64) -> f64) {
    filter(x, dt, |f| 1.0 / asd(f));
}

/// Standard normal samples (xorshift64* and Box–Muller).
fn gaussian(n: usize, seed: u64) -> Vec<f64> {
    // xorshift is linear over GF(2), so nearby seeds give dependent streams (seed 3 = seed 1 ⊕
    // seed 2): scramble the seed with SplitMix64 first.
    let mut z = seed.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    let mut s = (z ^ (z >> 31)).max(1);
    let mut uniform = move || {
        s ^= s >> 12;
        s ^= s << 25;
        s ^= s >> 27;
        ((s.wrapping_mul(0x2545_F491_4F6C_DD1D) >> 11) as f64 + 0.5) / (1u64 << 53) as f64
    };
    (0..n.div_ceil(2))
        .flat_map(|_| {
            let (u1, u2) = (uniform(), uniform());
            let r = (-2.0 * u1.ln()).sqrt();
            [r * (TAU * u2).cos(), r * (TAU * u2).sin()]
        })
        .take(n)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fft_matches_a_direct_dft_and_inverts() {
        let n = 64;
        let x: Vec<f64> = (0..n).map(|i| ((i * 7 % 13) as f64 - 6.0) * 0.3).collect();
        let (mut re, mut im) = (x.clone(), vec![0.0; n]);
        fft(&mut re, &mut im, false);
        for k in 0..n {
            let (mut dr, mut di) = (0.0, 0.0);
            for (t, v) in x.iter().enumerate() {
                let a = -TAU * (k * t) as f64 / n as f64;
                dr += v * a.cos();
                di += v * a.sin();
            }
            assert!((re[k] - dr).abs() < 1e-10 && (im[k] - di).abs() < 1e-10);
        }
        fft(&mut re, &mut im, true);
        assert!(re.iter().zip(&x).all(|(a, b)| (a - b).abs() < 1e-12) && im.iter().all(|v| v.abs() < 1e-12));
    }

    #[test]
    fn colored_noise_has_the_requested_spectrum() {
        let dt = 5.0;
        let n = 1 << 16;
        // White: per-sample std = asd·√(1/(2dt)).
        let w = colored_noise(|_| 1e-10, dt, n, 7);
        let std = (w.iter().map(|v| v * v).sum::<f64>() / n as f64).sqrt();
        assert!((std / (1e-10 / (2.0 * dt).sqrt()) - 1.0).abs() < 0.02, "{std:e}");
        // Coloured (∝ 1/f): band powers follow ∫asd²df = ln(f2/f1) in each band.
        let x = colored_noise(|f| 1e-10 / f.sqrt(), dt, n, 11);
        let (mut re, mut im) = (x.clone(), vec![0.0; n]);
        fft(&mut re, &mut im, false);
        let band = |f1: f64, f2: f64| -> f64 {
            (1..n / 2).filter(|&k| { let f = k as f64 / (n as f64 * dt); f >= f1 && f < f2 })
                .map(|k| 2.0 * (re[k] * re[k] + im[k] * im[k]) * dt / n as f64).sum::<f64>() / (n as f64 * dt)
        };
        for (f1, f2) in [(1e-3f64, 4e-3f64), (1e-2, 4e-2)] {
            let expect = 1e-20 * (f2 / f1).ln();
            assert!((band(f1, f2) / expect - 1.0).abs() < 0.1, "{f1}: {:e} vs {expect:e}", band(f1, f2));
        }
    }

    #[test]
    fn whitening_undoes_colouring() {
        let dt = 10.0;
        let asd = |f: f64| 1e-10 * (1.0 + 5e-3 / f).sqrt();
        let mut x = colored_noise(asd, dt, 1 << 15, 3);
        whiten(&mut x, dt, asd);
        let var = x.iter().map(|v| v * v).sum::<f64>() / x.len() as f64;
        assert!((var / (1.0 / (2.0 * dt)) - 1.0).abs() < 0.03, "{var}");
    }
}
