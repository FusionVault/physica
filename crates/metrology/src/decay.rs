//! Radioactive decay: the one-line laws other crates share. A source in steady state, made as
//! fast as it decays, has an activity equal to its production rate, so that case needs no
//! function.

use core::f64::consts::LN_2;

/// One curie, Bq (exact): the activity of about one gram of ²²⁶Ra.
pub const CURIE: f64 = 3.7e10;

/// Activity `N ln2 / t½`, Bq, of `n` nuclei with half-life `half_life` (s).
pub fn activity(n: f64, half_life: f64) -> f64 {
    n * LN_2 / half_life
}

/// The permanent relative step in a source's activity after its decay stops for `duration`
/// (same units as `half_life`): `exp(ln2·Δt/t½) − 1`. Afterwards it is that much stronger than
/// a source that kept decaying, so claimed suppressions can be checked against the long-term
/// stability of routine activity measurements.
pub fn activity_step_after_suppression(half_life: f64, duration: f64) -> f64 {
    (LN_2 * duration / half_life).exp_m1()
}

/// Relative half-life change `(T₁ − T₂)/T₁` between two environments, as measured for
/// electron capture, whose rate depends on the atom's own electrons.
pub fn ec_half_life_shift(t1: f64, t2: f64) -> f64 {
    (t1 - t2) / t1
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{AMU, N_A, YEAR};

    #[test]
    fn radium_is_about_one_curie_per_gram() {
        // 1 g of ²²⁶Ra (226.025 u, t½ 1600 yr) is 0.988 Ci, 3.66×10¹⁰ Bq.
        let n = 1e-3 / (226.025_410 * AMU);
        let a = activity(n, 1600.0 * YEAR);
        assert!((a / 3.657e10 - 1.0).abs() < 1e-3 && (a / CURIE - 0.988).abs() < 1e-3);
        assert!((n / (N_A / 226.025_410) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn a_two_day_stop_of_cobalt_60_leaves_a_permanent_step() {
        // ⁶⁰Co, t½ = 1925.28 d (IAEA): 7.2×10⁻⁴, 9–120 times the 6×10⁻⁶–8×10⁻⁵ limits on
        // annual modulation from 14 laboratories (Pommé et al. 2016).
        let step = activity_step_after_suppression(1925.28, 2.0);
        assert!((step / 7.2e-4 - 1.0).abs() < 1e-2 && step / 8e-5 > 9.0 && step / 6e-6 < 121.0);
        assert!((activity_step_after_suppression(5.0, 5.0) - 1.0).abs() < 1e-15); // one half-life: twice as strong
    }

    #[test]
    fn beryllium_7_inside_c60() {
        // Ohtsuki et al., PRL 93, 112501 (2004): 52.68 d in C₆₀ against 53.12 d in Be metal.
        assert!((ec_half_life_shift(53.12, 52.68) - 0.0083).abs() < 1e-4);
    }
}
