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

/// In-place FFT of a row-major array of shape `dims` (the last index fastest), along each axis in
/// turn with [`fft`]; `inverse` divides by the total length.
///
/// # Panics
/// If a dimension is not a power of two or the lengths differ from the shape's.
pub fn fftn(re: &mut [f64], im: &mut [f64], dims: &[usize], inverse: bool) {
    let n: usize = dims.iter().product();
    assert!(re.len() == n && im.len() == n, "one value per grid point");
    let mut stride = n;
    for &len in dims {
        stride /= len;
        if stride == 1 {
            for (r, i) in re.chunks_mut(len).zip(im.chunks_mut(len)) {
                fft(r, i, inverse);
            }
            continue;
        }
        let (mut lr, mut li) = (vec![0.0; len], vec![0.0; len]);
        for start in (0..n).step_by(len * stride).flat_map(|o| o..o + stride) {
            for k in 0..len {
                (lr[k], li[k]) = (re[start + k * stride], im[start + k * stride]);
            }
            fft(&mut lr, &mut li, inverse);
            for k in 0..len {
                (re[start + k * stride], im[start + k * stride]) = (lr[k], li[k]);
            }
        }
    }
}

/// The FFT's wavenumber (rad per unit length) at index `i` of `n` samples `d` apart, negative
/// from `n/2` up.
pub fn wavenumber(i: usize, n: usize, d: f64) -> f64 {
    let i = if i < n / 2 { i as f64 } else { i as f64 - n as f64 };
    TAU * i / (n as f64 * d)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fftn_matches_a_direct_dft_and_inverts() {
        // A 4 × 2 × 8 array against the triple sum, then back.
        let dims = [4, 2, 8];
        let x: Vec<f64> = (0..64).map(|i| ((i * 11 % 17) as f64 - 8.0) * 0.2).collect();
        let (mut re, mut im) = (x.clone(), vec![0.0; 64]);
        fftn(&mut re, &mut im, &dims, false);
        let split = |i: usize| [i / 16, i / 8 % 2, i % 8];
        for k in 0..64 {
            let (mut dr, mut di) = (0.0, 0.0);
            for (t, v) in x.iter().enumerate() {
                let a = -TAU * (0..3).map(|d| (split(k)[d] * split(t)[d]) as f64 / dims[d] as f64).sum::<f64>();
                dr += v * a.cos();
                di += v * a.sin();
            }
            assert!((re[k] - dr).abs() < 1e-12 && (im[k] - di).abs() < 1e-12);
        }
        fftn(&mut re, &mut im, &dims, true);
        assert!(re.iter().zip(&x).all(|(a, b)| (a - b).abs() < 1e-14) && im.iter().all(|v| v.abs() < 1e-14));
        let pi = core::f64::consts::PI;
        assert_eq!([wavenumber(3, 8, 0.5), wavenumber(4, 8, 0.5), wavenumber(7, 8, 0.5)], [1.5 * pi, -2.0 * pi, -0.5 * pi]);
    }

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
