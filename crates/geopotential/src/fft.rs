//! The fast Fourier transform, for gridded fields and time series alike.

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
}
