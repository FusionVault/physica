//! Toroidal moments, anapoles, supertoroidal pulses and the Aharonov–Bohm phase: the third family
//! of multipoles, next to the electric and magnetic ones, and what it does and does not do.
//!
//! A current wound poloidally on a torus has no electric or magnetic dipole, but it has a toroidal
//! dipole `T = (1/10)∫[(r·J)r − 2r²J] dV` (SI, A m³; Gaussian texts put `1/c` in front). Static,
//! it makes no field and no scalar potential outside itself, so it neither feels nor exerts a force
//! on external static fields, and it cannot be one of [`geopotential::Field`]'s Gauss
//! coefficients. Oscillating as `e^(−iωt)`, it radiates exactly like an electric dipole
//! `p_eff = p + (ik/c)T`, so its power rises as `ω⁶` where a dipole's rises as `ω⁴`, and
//! `p = −(ik/c)T` cancels the far field: the anapole. Outside an ideal toroidal solenoid `B = 0`
//! but `∮A·dl = Φ`, and an electron going through the hole gains the phase `eΦ/ħ`, which Tonomura's
//! group measured as 0 or π because the trapped flux comes in units of `h/2e`. The phase is not an
//! energy: the potentials alone do no work. Status of all of this: **mainstream**.
//!
//! ```
//! use magnetics::toroidal::{aharonov_bohm_phase, toroidal_radiated_power, toroidal_solenoid_moment};
//! use metrology::{E_CHARGE, FLUX_QUANTUM};
//!
//! // One superconducting flux quantum h/2e through the hole shifts an electron by π.
//! assert!((aharonov_bohm_phase(E_CHARGE, FLUX_QUANTUM) - core::f64::consts::PI).abs() < 1e-12);
//! // 24 turns of 1 A, minor radius 0.3 m, on a 1 m torus: T = NIV/4π = 3.393 A m³.
//! assert!((toroidal_solenoid_moment(24.0, 1.0, 1.0, 0.3) - 3.3929).abs() < 1e-4);
//! // A cortical "toroid" of 10 nAm loops (T = 2.5×10⁻¹³ A m³) radiates 5×10⁻⁶⁵ W at 10 Hz.
//! let p = toroidal_radiated_power(2.5e-13, core::f64::consts::TAU * 10.0);
//! assert!((p / 5.293e-65 - 1.0).abs() < 1e-3);
//! ```
//!
//! Sources: Dubovik & Cheshkov, [Fiz. Elem. Chast. At. Yadra 5, 791 (1974)](https://www1.jinr.ru/Archive/Pepan/1974-v5/v-5-3/6.htm)
//! (eq. 37, `E_lm = Q̇_lm + k²T_lm`); Dubovik & Tugushev, [Phys. Rep. 187, 145 (1990)](https://doi.org/10.1016/0370-1573(90)90042-Z)
//! (not opened); Miroshnichenko et al., [Nat. Commun. 6, 8069 (2015)](https://doi.org/10.1038/ncomms9069)
//! (eqs. 5, 8, 10–12); Fernandez-Corbaton, Nanz & Rockstuhl, [arXiv:1507.00755](https://arxiv.org/abs/1507.00755)
//! (eqs. 5–6, the SI `1/10` form for `e^(−iωt)`); Nemkov, Basharin & Fedotov,
//! [Phys. Rev. A 98, 023858 (2018)](https://arxiv.org/abs/1804.07342) (eq. 18, App. B); Yang &
//! Bozhevolnyi, [arXiv:1810.01221](https://arxiv.org/abs/1810.01221) (eq. 27); Shen et al.,
//! [Nat. Commun. 12, 5891 (2021)](https://doi.org/10.1038/s41467-021-26037-w); Planes, Castán &
//! Saxena, [arXiv:1504.04958](https://arxiv.org/abs/1504.04958); Osakabe et al.,
//! [Phys. Rev. A 34, 815 (1986)](https://doi.org/10.1103/PhysRevA.34.815).

use core::f64::consts::PI;

use metrology::{C, HBAR, MU_0};

pub(crate) fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub(crate) fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

/// Toroidal dipole (A m³) of a closed circuit given as segments `(r, I·dl)`, position (m) and
/// current element (A m): `T = (1/10)Σ[(r·I dl)r − 2r² I dl]`, the SI form of Dubovik's moment.
/// For a closed current it equals `½Σ(r·I dl)r` and `½∫r × M dV` for any `M` with `∇×M = J`.
/// Status: **mainstream**.
pub fn toroidal_dipole(segs: &[([f64; 3], [f64; 3])]) -> [f64; 3] {
    let mut t = [0.0; 3];
    for &(r, idl) in segs {
        let (rj, rr) = (dot(r, idl), dot(r, r));
        for k in 0..3 {
            t[k] += 0.1 * (rj * r[k] - 2.0 * rr * idl[k]);
        }
    }
    t
}

/// Toroidal dipole (A m³) of `n` turns carrying `i` (A) wound on a torus of major radius `major`
/// and minor radius `minor` (m): `NIπRa²/2 = NIV/(4π)`, along the axis. Nemkov et al. (2018,
/// App. B) give `IV/(4πc)` in Gaussian units. Status: **mainstream**.
pub fn toroidal_solenoid_moment(n: f64, i: f64, major: f64, minor: f64) -> f64 {
    n * i * PI * major * minor * minor / 2.0
}

/// The electric dipole (C m) that radiates like `p` and `t` together, one Cartesian component of
/// each as `[re, im]`, at wavenumber `k` (1/m), for time dependence `e^(−iωt)`:
/// `p_eff = p + (ik/c)T`. The anapole is `p = −(ik/c)T`; for `e^(+iωt)` the sign of `i` flips.
/// The SI form of Nemkov et al.'s Gaussian `d + ikτ = 0` and Miroshnichenko et al.'s eq. 12.
/// Outside the source only this sum is observable: no external measurement at one frequency
/// separates `p` from `T` (Fernandez-Corbaton, Nanz & Rockstuhl 2017). Status: **mainstream**.
pub fn effective_dipole(p: [f64; 2], t: [f64; 2], k: f64) -> [f64; 2] {
    [p[0] - k / C * t[1], p[1] + k / C * t[0]]
}

/// Power (W) radiated by a toroidal dipole of amplitude `t` (A m³) oscillating at `omega`
/// (rad/s): `μ₀ω⁶T²/(12πc⁵)`, the electric-dipole power `μ₀ω⁴|p|²/(12πc)` of `p = (ik/c)T`.
/// The `ω⁶` law assumes a current pattern that keeps its shape across frequency; a pattern that
/// deforms at `O(k²)` can mimic `T` exactly. Status: **mainstream**.
pub fn toroidal_radiated_power(t: f64, omega: f64) -> f64 {
    MU_0 * omega.powi(6) * t * t / (12.0 * PI * C.powi(5))
}

/// The smallest diameter (same unit as `wavelength`) of a sphere of relative refractive index `n`
/// at which the Mie electric-dipole coefficient `a₁` vanishes, the anapole of Miroshnichenko et
/// al. (2015): the first root of `nψ₁(nx)ψ₁′(x) − ψ₁(x)ψ₁′(nx)`, `x = πd/λ`, with the
/// Riccati–Bessel `ψ₁(z) = sin z/z − cos z`. NaN if there is none below `x = 100`.
/// Status: **mainstream**.
pub fn mie_a1_zero_diameter(n: f64, wavelength: f64) -> f64 {
    let psi = |z: f64| z.sin() / z - z.cos();
    let dpsi = |z: f64| z.cos() / z - z.sin() / (z * z) + z.sin();
    let g = |x: f64| n * psi(n * x) * dpsi(x) - psi(x) * dpsi(n * x);
    // Steps much finer than the period π/n of ψ₁(nx), so no pair of roots is skipped.
    let step = 0.01 / n.max(1.0);
    let mut lo = step;
    while g(lo) * g(lo + step) > 0.0 {
        lo += step;
        if lo > 100.0 {
            return f64::NAN;
        }
    }
    let mut hi = lo + step;
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        if g(lo) * g(mid) <= 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    0.5 * (lo + hi) * wavelength / PI
}

/// The next moment up from point moments `moments` at `centres` (m): `½Σ rᵢ × vᵢ`. Magnetic
/// dipoles in give the toroidal dipole; toroidal dipoles in give the hypertoroidal moment of
/// Prosandeev & Bellaiche (2009), and so on, the parity alternating with each order (polar,
/// axial). It depends on the origin unless the lower moments sum to zero.
///
/// Planes et al. (2015, eq. 29) write the toroidal moment as `¼∫r × (r × J) dV`, which is 3/2 of
/// this for any closed current; `½∫r × M dV` needs a true magnetisation (`∇×M = J`), not
/// `½r × J`. Status: **mainstream** (definition); the golden-ratio radii and "n ≥ 3" rules of
/// fractal toroids are **unverified** (in no source found).
pub fn hypertoroidal_moment(centres: &[[f64; 3]], moments: &[[f64; 3]]) -> [f64; 3] {
    let mut h = [0.0; 3];
    for (&r, &v) in centres.iter().zip(moments) {
        let c = cross(r, v);
        for k in 0..3 {
            h[k] += 0.5 * c[k];
        }
    }
    h
}

/// Aharonov–Bohm phase (rad) of a charge `charge` (C) taken round a flux `flux` (Wb): `qΦ/ħ`.
/// One `h/2e` gives π for an electron, so Tonomura's superconductor-clad toroid showed 0 or π.
/// Status: **mainstream**.
pub fn aharonov_bohm_phase(charge: f64, flux: f64) -> f64 {
    charge * flux / HBAR
}

/// EMF (V) round any loop threading a flux changing at `dflux_dt` (Wb/s): `∮E·dl = −dΦ/dt`.
/// So `E` cannot vanish everywhere on a loop round a toroid whose flux changes, which is what
/// "E = B = 0 signalling" with potentials alone (US 5,845,220) would need. Status: Faraday's law,
/// **mainstream**; the patent's claim **contested**.
pub fn loop_emf(dflux_dt: f64) -> f64 {
    -dflux_dt
}

/// The supertoroidal generating function `[re, im]`, with `f₀ = 1`, at cylindrical radius `r` and
/// height `z` (m) and time `t` (s): `f = 1/[(q₁ + iτ)(s + q₂)^α]`, `s = r²/(q₁ + iτ) − iσ`,
/// `τ = z − ct`, `σ = z + ct` (Shen et al. 2021). `α = 1` is the Hellwarth–Nouchi pulse
/// `1/[r² + (q₁ + iτ)(q₂ − iσ)]`. The pulse's energy is finite for every `α > 0` (Shen et al.'s
/// `α ≥ 1` is sufficient, not necessary; derived in the toroidal-novelty research, #23). Any `g(s)/(q₁ + iτ)` solves the
/// wave equation, and `Re(s + q₂) > 0` for `q₁, q₂ > 0`, so the principal power used here is
/// analytic everywhere. Status: **mainstream**.
pub fn supertoroidal_generating_fn(r: f64, z: f64, t: f64, q1: f64, q2: f64, alpha: f64) -> [f64; 2] {
    let (tau, sigma) = (z - C * t, z + C * t);
    let a2 = q1 * q1 + tau * tau;
    let w = [r * r * q1 / a2 + q2, -r * r * tau / a2 - sigma];
    let mag = 1.0 / (a2.sqrt() * w[0].hypot(w[1]).powf(alpha));
    let phase = -(tau.atan2(q1) + alpha * w[1].atan2(w[0]));
    [mag * phase.cos(), mag * phase.sin()]
}

#[cfg(test)]
mod tests {
    use super::*;
    use metrology::{DIRAC_FLUX, E_CHARGE, FLUX_QUANTUM};

    fn close(a: f64, b: f64, rel: f64) -> bool {
        (a / b - 1.0).abs() < rel
    }

    /// `turns` circular turns of `m` segments each, minor radius `a`, on a torus of major radius
    /// `big_r` about the z axis, each turn carrying `i`: the check script's toroidal solenoid.
    fn solenoid(turns: usize, i: f64, big_r: f64, a: f64, m: usize) -> Vec<([f64; 3], [f64; 3])> {
        let mut segs = Vec::new();
        for j in 0..turns {
            let ph = 2.0 * PI * (j as f64 + 0.5) / turns as f64;
            let pt = |k: usize| {
                let th = 2.0 * PI * k as f64 / m as f64;
                let rho = big_r + a * th.cos();
                [rho * ph.cos(), rho * ph.sin(), a * th.sin()]
            };
            for k in 0..m {
                let (p, q) = (pt(k), pt(k + 1));
                segs.push(([0.5 * (p[0] + q[0]), 0.5 * (p[1] + q[1]), 0.5 * (p[2] + q[2])], [i * (q[0] - p[0]), i * (q[1] - p[1]), i * (q[2] - p[2])]));
            }
        }
        segs
    }

    #[test]
    fn toroidal_solenoid_moment_matches_the_segment_sum() {
        // Script: |T| = N I V/(4π) to 1e-3 (it finds 1.6×10⁻⁴, the polygon's chord error).
        for a in [0.05, 0.5] {
            let t = toroidal_dipole(&solenoid(24, 1.0, 1.0, a, 200));
            assert!(close(t[2].abs(), toroidal_solenoid_moment(24.0, 1.0, 1.0, a), 1e-3), "a = {a}");
            assert!(t[0].abs() < 1e-12 && t[1].abs() < 1e-12);
        }
        // For a closed current the 1/10 form equals ½Σ(r·I dl)r, and Planes et al.'s ¼Σr×(r×I dl)
        // is 3/2 of it.
        let segs = solenoid(24, 1.0, 1.0, 0.3, 200);
        let t = toroidal_dipole(&segs)[2];
        let half: f64 = segs.iter().map(|&(r, idl)| 0.5 * dot(r, idl) * r[2]).sum();
        let planes: f64 = segs.iter().map(|&(r, idl)| 0.25 * cross(r, cross(r, idl))[2]).sum();
        assert!(close(half, t, 1e-12) && close(planes, 1.5 * t, 1e-12));
    }

    #[test]
    fn ring_identity_and_hypertoroids() {
        // ½Σ rᵢ × mᵢ over the turns' centres and magnetic moments is the toroidal dipole.
        let (turns, m) = (24, 200);
        let segs = solenoid(turns, 1.0, 1.0, 0.3, m);
        let mut centres = Vec::new();
        let mut moments = Vec::new();
        for (j, turn) in segs.chunks(m).enumerate() {
            let ph = 2.0 * PI * (j as f64 + 0.5) / turns as f64;
            let c = [ph.cos(), ph.sin(), 0.0];
            let mut mj = [0.0; 3];
            for &(r, idl) in turn {
                let x = cross([r[0] - c[0], r[1] - c[1], r[2] - c[2]], idl);
                (0..3).for_each(|k| mj[k] += 0.5 * x[k]);
            }
            centres.push(c);
            moments.push(mj);
        }
        assert!(close(hypertoroidal_moment(&centres, &moments)[2], toroidal_dipole(&segs)[2], 1e-12));
        // A ring of N azimuthal toroidal dipoles T at radius ρ: ½NρT along the axis.
        let ring: Vec<_> = (0..8).map(|j| f64::from(j) * PI / 4.0).collect();
        let c: Vec<_> = ring.iter().map(|p| [2.0 * p.cos(), 2.0 * p.sin(), 0.0]).collect();
        let t: Vec<_> = ring.iter().map(|p| [-3.0 * p.sin(), 3.0 * p.cos(), 0.0]).collect();
        let h = hypertoroidal_moment(&c, &t);
        assert!(close(h[2], 0.5 * 8.0 * 2.0 * 3.0, 1e-12) && h[0].abs() < 1e-12 && h[1].abs() < 1e-12);
    }

    /// Far-field power (W) of the segments plus a point dipole `p_z` (`[re, im]`) at the origin,
    /// at wavenumber `k`: `μ₀ck²/(32π²)∮|n × F|² dΩ`, `F = Σ I dl e^(−ik n·r) − iω p`.
    fn far_field_power(segs: &[([f64; 3], [f64; 3])], k: f64, p_z: [f64; 2]) -> f64 {
        let (nt, np) = (60, 24);
        let w = k * C;
        let mut total = 0.0;
        for it in 0..nt {
            let th = PI * (it as f64 + 0.5) / nt as f64;
            for ip in 0..np {
                let ph = 2.0 * PI * (ip as f64 + 0.5) / np as f64;
                let n = [th.sin() * ph.cos(), th.sin() * ph.sin(), th.cos()];
                let (mut re, mut im) = ([0.0; 3], [0.0; 3]);
                for &(r, idl) in segs {
                    let a = -k * dot(n, r);
                    (0..3).for_each(|q| {
                        re[q] += idl[q] * a.cos();
                        im[q] += idl[q] * a.sin();
                    });
                }
                re[2] += w * p_z[1];
                im[2] -= w * p_z[0];
                let sq = dot(cross(n, re), cross(n, re)) + dot(cross(n, im), cross(n, im));
                total += sq * th.sin() * (PI / nt as f64) * (2.0 * PI / np as f64);
            }
        }
        MU_0 * C * k * k / (32.0 * PI * PI) * total
    }

    #[test]
    fn radiation_and_the_anapole() {
        // Script: 24 turns, a/R = 0.3, kR = 0.05; power to 1e-2 (it finds 5.5×10⁻⁴).
        let segs = solenoid(24, 1.0, 1.0, 0.3, 60);
        let t = toroidal_dipole(&segs)[2];
        let k = 0.05;
        let p1 = far_field_power(&segs, k, [0.0; 2]);
        assert!(close(p1, toroidal_radiated_power(t, k * C), 1e-2));
        // ω⁶: exactly 64 from the formula, 63.9 from the numerical far field (script 63.8955).
        assert!(close(toroidal_radiated_power(t, 2.0 * k * C) / toroidal_radiated_power(t, k * C), 64.0, 1e-12));
        assert!(close(far_field_power(&segs, 2.0 * k, [0.0; 2]) / p1, 64.0, 1e-2));
        // The anapole p = −(ik/c)T leaves < 1e-3 of the power (script 7.9×10⁻⁸); p_eff is zero.
        let ana = [0.0, -k / C * t];
        assert_eq!(effective_dipole(ana, [t, 0.0], k), [0.0, 0.0]);
        assert!(far_field_power(&segs, k, ana) / p1 < 1e-3);
        // The opposite sign doubles p_eff, so 4× the power (script 4.001): the sign of
        // (ik/c)T is right for e^(−iωt).
        let wrong = [0.0, k / C * t];
        let pe = effective_dipole(wrong, [t, 0.0], k);
        let dipole_power = MU_0 * (k * C).powi(4) * (pe[0] * pe[0] + pe[1] * pe[1]) / (12.0 * PI * C);
        assert!(close(dipole_power / p1, 4.0, 1e-2));
        assert!(close(far_field_power(&segs, k, wrong) / p1, 4.0, 1e-2));
    }

    #[test]
    fn silicon_sphere_anapole() {
        // n = 4 at 550 nm: "the diameter ∼204 nm" (Miroshnichenko et al. 2015; script 204.028).
        let d = mie_a1_zero_diameter(4.0, 550.0);
        assert!((d - 204.0).abs() < 3.0 && (d - 204.028).abs() < 1e-3, "{d}");
    }

    #[test]
    fn aharonov_bohm_and_faraday() {
        // Osakabe et al. 1986: "0 or π" for flux in units of h/2e; h/e (Dirac) is 2π.
        assert!(close(aharonov_bohm_phase(E_CHARGE, FLUX_QUANTUM), PI, 1e-12));
        assert!(close(aharonov_bohm_phase(E_CHARGE, DIRAC_FLUX), 2.0 * PI, 1e-12));
        assert_eq!(loop_emf(2.5), -2.5);
    }

    #[test]
    fn supertoroidal_pulses_solve_the_wave_equation() {
        let (q1, q2) = (1.0, 20.0);
        // α = 1 is Hellwarth–Nouchi's 1/[r² + (q₁ + iτ)(q₂ − iσ)].
        for (r, z, t) in [(0.3, 0.2, 1e-9), (2.0, -1.5, -3e-9), (0.0, 0.0, 0.0)] {
            let f = supertoroidal_generating_fn(r, z, t, q1, q2, 1.0);
            let (tau, sigma) = (z - C * t, z + C * t);
            let d = [r * r + q1 * q2 + tau * sigma, q2 * tau - q1 * sigma];
            let dd = d[0] * d[0] + d[1] * d[1];
            assert!((f[0] - d[0] / dd).hypot(f[1] + d[1] / dd) < 1e-12 / dd.sqrt());
        }
        // (∇² − ∂²/c²∂t²)f = 0 by central differences, against the size of its terms.
        for alpha in [1.0, 2.5, 5.0] {
            for (r, z, t) in [(0.3, 0.2, 1e-9), (2.0, -1.5, -3e-9), (1.0, 4.0, 2e-8)] {
                let f = |r: f64, z: f64, t: f64| supertoroidal_generating_fn(r, z, t, q1, q2, alpha);
                let h = 1e-3;
                let ht = h / C;
                for q in 0..2 {
                    let f0 = f(r, z, t)[q];
                    let frr = (f(r + h, z, t)[q] - 2.0 * f0 + f(r - h, z, t)[q]) / (h * h);
                    let fr = (f(r + h, z, t)[q] - f(r - h, z, t)[q]) / (2.0 * h);
                    let fzz = (f(r, z + h, t)[q] - 2.0 * f0 + f(r, z - h, t)[q]) / (h * h);
                    let ftt = (f(r, z, t + ht)[q] - 2.0 * f0 + f(r, z, t - ht)[q]) / (ht * ht * C * C);
                    let box_f = frr + fr / r + fzz - ftt;
                    let scale = frr.abs() + (fr / r).abs() + fzz.abs() + ftt.abs();
                    assert!(box_f.abs() < 1e-5 * scale, "α = {alpha}, ({r}, {z}, {t}): {box_f} vs {scale}");
                }
            }
        }
    }
}
