//! The construction in OpenAI's *Velocity-Field Detection of Computation in Forced Navier–Stokes
//! Flows* (27 September 2026), at a size a computer can hold: the integer-stack machine encoding
//! (eqs. 22–24), both scale schedules (eqs. 10 and 25) with their inequalities, and a run of the
//! torus mechanism (§3) through [`magnetics::navier_stokes`]'s solver.
//!
//! ```text
//! cargo run -p magnetics --example navier_stokes_computation
//! ```
//!
//! The paper's own scales cannot be simulated: its first bump is 1.25×10⁻⁷ wide and its first
//! stirring burst lasts ~10⁻²⁰. So the run uses a relaxed bump and a uniform translation for the
//! planar processor, and checks that the solver keeps the bounds the proof relies on. It does not
//! and cannot check the halting theorem.

use std::collections::HashSet;

use magnetics::navier_stokes::{Torus, bump, bump_derivatives, smooth_step, smooth_step_derivatives, torus_heat_kernel};

/// A head move.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Move {
    Left,
    Stay,
    Right,
}

/// One rule of a one-head machine: in `state` reading digit `read`, write `write`, move, then
/// enter `next`.
#[derive(Clone, Copy, Debug)]
struct Rule {
    state: usize,
    read: u64,
    write: u64,
    step: Move,
    next: usize,
}

/// A four-state, two-symbol table with the compiler's properties: state 4 is terminal; each
/// target state fixes its incoming move; target and written symbol identify the rule. From a
/// blank tape it halts after three steps.
const TABLE: [Rule; 6] = [
    Rule { state: 1, read: 0, write: 1, step: Move::Right, next: 2 },
    Rule { state: 1, read: 1, write: 0, step: Move::Right, next: 2 },
    Rule { state: 2, read: 0, write: 0, step: Move::Left, next: 3 },
    Rule { state: 2, read: 1, write: 1, step: Move::Left, next: 3 },
    Rule { state: 3, read: 0, write: 1, step: Move::Stay, next: 1 },
    Rule { state: 3, read: 1, write: 0, step: Move::Stay, next: 4 },
];
const STATES: usize = 4;
const SYMBOLS: u64 = 2;
/// Grid points a side for the torus run, and its relaxed bump's half-width (the paper's first
/// is 1.25×10⁻⁷). 32² keeps the example's test near two seconds in debug; it represents the bump
/// to a few per cent (64² takes four times the steps on four times the points).
const N: usize = 32;
const WIDTH: f64 = 0.25;

/// Eq. 22: the stacks `(x, y)` after `rule`. `x` holds the tape left of the head, nearest cell
/// lowest; `y` the head cell and the tape to its right, digits base `SYMBOLS`.
fn update(rule: &Rule, x: u64, y: u64) -> (u64, u64) {
    let b = SYMBOLS;
    match rule.step {
        Move::Right => (b * x + rule.write, y / b),
        Move::Stay => (x, y - rule.read + rule.write),
        Move::Left => (x / b, x % b + b * rule.write + b * b * (y / b)),
    }
}

/// The inverse of eq. 22 from the new state and stacks: the move comes from the state, the
/// written digit from the stacks, then the rule from both; `None` if no rule leads there.
fn undo(next: usize, x: u64, y: u64) -> Option<(usize, u64, u64)> {
    let b = SYMBOLS;
    let step = TABLE.iter().find(|r| r.next == next)?.step;
    let written = match step {
        Move::Right => x % b,
        Move::Stay => y % b,
        Move::Left => y / b % b,
    };
    let rule = TABLE.iter().find(|r| r.next == next && r.write == written)?;
    let (a, beta) = (rule.read, rule.write);
    Some(match step {
        Move::Right => (rule.state, (x - beta) / b, b * y + a),
        Move::Stay => (rule.state, x, y - beta + a),
        Move::Left => (rule.state, b * x + y % b, b * (y / (b * b)) + a),
    })
}

/// Eq. 24: the address `k = 1 + x + B y + B²(i − 1)` of a configuration at level width `B`.
fn address(i: usize, x: u64, y: u64, big_b: u64) -> u64 {
    1 + x + big_b * y + big_b * big_b * (i as u64 - 1)
}

/// The configuration `(i, x, y)` at address `k`, by integer division.
fn configuration(k: u64, big_b: u64) -> (usize, u64, u64) {
    let k = k - 1;
    ((k / (big_b * big_b)) as usize + 1, k % big_b, k / big_b % big_b)
}

fn encoding() {
    // The compiler's properties of the table.
    for next in 1..=STATES {
        let incoming: Vec<&Rule> = TABLE.iter().filter(|r| r.next == next).collect();
        assert!(incoming.windows(2).all(|w| w[0].step == w[1].step), "state {next}: one move");
        assert!(incoming.iter().map(|r| r.write).collect::<HashSet<_>>().len() == incoming.len(), "state {next}: injective");
    }
    assert!(TABLE.iter().all(|r| r.state != STATES));
    // Every configuration with stacks below B_n = b^(m+n) (m = 1, the blank tape's stacks are 0)
    // goes to one below B_{n+1}, with y′ ≤ bB_n − 1, distinct sources to distinct targets, and
    // the inverse recovers the source.
    let m = 1;
    let mut checked = 0;
    for n in 1..=4u32 {
        let (big_b, next_b) = (SYMBOLS.pow(m + n), SYMBOLS.pow(m + n + 1));
        let mut targets = HashSet::new();
        for i in 1..STATES {
            for x in 0..big_b {
                for y in 0..big_b {
                    let rule = TABLE.iter().find(|r| r.state == i && r.read == y % SYMBOLS).expect("every non-terminal pair has a rule");
                    let (x2, y2) = update(rule, x, y);
                    assert!(x2 < next_b && y2 < next_b && y2 < SYMBOLS * big_b);
                    let k = address(rule.next, x2, y2, next_b);
                    assert!((1..=STATES as u64 * next_b * next_b).contains(&k) && configuration(k, next_b) == (rule.next, x2, y2));
                    assert!(targets.insert(k), "two sources reach address {k}");
                    assert_eq!(undo(rule.next, x2, y2), Some((i, x, y)));
                    checked += 1;
                }
            }
        }
    }
    // Run it from the blank tape: 1 → 2 → 3 → 4, the terminal state.
    let (mut i, mut x, mut y, mut steps) = (1, 0, 0, 0);
    while i != STATES {
        let rule = TABLE.iter().find(|r| r.state == i && r.read == y % SYMBOLS).expect("defined");
        (x, y) = update(rule, x, y);
        (i, steps) = (rule.next, steps + 1);
    }
    println!("Integer stacks (eqs. 22–24): {checked} configurations on levels 1–4 map injectively to");
    println!("  the next level's addresses, each undone by the inverse rows; the table halts in {steps} steps.");
    assert_eq!((steps, x, y), (3, 0, 0));
}

/// The bump's derivative constant `C` (§3.1), which the paper leaves unstated: the ceiling of
/// the larger of `max‖∇g‖ = 4` and `max‖D²g‖ ≈ 40.46` (found by the constants check of #53).
fn bump_constant() -> f64 {
    let (mut grad, mut hess) = (0f64, 0f64);
    for i in 0..=800 {
        for j in 0..=800 {
            let z = [-1.0 + i as f64 / 400.0, -1.0 + j as f64 / 400.0];
            let (_, g, h) = bump_derivatives(z);
            grad = grad.max(g[0].hypot(g[1]));
            let (mean, half) = ((h[0][0] + h[1][1]) / 2.0, ((h[0][0] - h[1][1]) / 2.0).hypot(h[0][1]));
            hess = hess.max(mean.abs() + half);
        }
    }
    assert!((grad - 4.0).abs() < 1e-12 && (hess - 40.46).abs() < 0.01, "{grad} {hess}");
    grad.max(hess).ceil()
}

fn torus_schedule(c: f64) {
    // Eq. 10 with R = 4L, in logarithms (D_n reaches 10^(10³) within a few hundred blocks).
    println!("\nTorus schedule (eq. 10), C = {c}:");
    for l in [1.0f64, 2.0, 5.0] {
        let r = 4.0 * l;
        let ln_b = |n: f64| (1e-6f64).ln() - n * (2.0 * r).ln();
        let ln_d = |n: f64| (2.0 * c).ln() - 2.0 * ln_b(n) + (1.0 + n * l).ln() + 3.0 * n * r.ln();
        // M₀ = Σ4b_n² = 4×10⁻¹²/(4R² − 1) < 10⁻¹¹; the bump's true mass is 2.25b_n².
        let m0 = 4e-12 / (4.0 * r * r - 1.0);
        let partial: f64 = (1..200).map(|n| 4.0 * (2.0 * ln_b(n as f64)).exp()).sum();
        assert!(m0 < 1e-11 && (partial / m0 - 1.0).abs() < 1e-12);
        for h in [1e9, 1e8] {
            assert!(2.0 * h * m0 + 1.0 / 16.0 < 0.5, "eq. 20");
        }
        for n in [1.0, 2.0, 10.0, 100.0] {
            // 2δ_nD_n = (1/50)·D/(1 + D), short of 1/50 by 1/(50(1 + D)): exact, but a tie in f64.
            let margin_log10 = -(50f64.ln() + ln_d(n)) / 10f64.ln();
            // √2 b_n Rⁿ = √2×10⁻⁶ 2⁻ⁿ < 1/32 (eq. 19's separation).
            let spread = (2f64.sqrt().ln() + ln_b(n) + n * r.ln()).exp();
            assert!(spread < 1.0 / 32.0 && (spread / (2f64.sqrt() * 1e-6 * 0.5f64.powf(n)) - 1.0).abs() < 1e-9);
            if n == 1.0 {
                println!("  L = {l}: b₁ = {:.3e}, D₁ = 10^{:.1}, δ₁ = 10^{:.1}, 1/50 − 2δ₁D₁ = 10^{margin_log10:.1}, M₀ = {m0:.2e}",
                    ln_b(1.0).exp(), ln_d(1.0) / 10f64.ln(), -(100f64.ln() + ln_d(1.0)) / 10f64.ln());
            }
        }
        // Eq. 16 needs e^(Ln) ≤ Rⁿ, that is e^L ≤ R = 4L, which holds for L = 1, 2 only; R = ⌈e^L⌉
        // (or any R ≥ max(4L, e^L)) repairs it and changes nothing else in §3.
        let holds = l.exp() <= r;
        println!("    eq. 16's e^L ≤ R = 4L: {:.2} ≤ {r}: {holds}", l.exp());
        assert_eq!(holds, l < 3.0);
    }
}

fn array_schedule() {
    // Eq. 25 for the table above: N = 4 states, b = 2 symbols, m = 1.
    let (n_states, base, m) = (STATES as f64, SYMBOLS as f64, 1);
    let (c_star, growth) = (3000.0, base * base); // the paper's C∗ and D = b²
    let lambda = 4.0 * growth;
    let k0 = n_states * base.powi(2 * m);
    assert!(2880.0 <= c_star && growth / lambda == 0.25 && lambda >= 4.0);
    println!("\nExpanding-array schedule (eq. 25), N = {n_states}, b = {base}, m = {m}:");
    for nu in [1e-3, 1.0, 1e3] {
        let r0 = 1024.0 * c_star * (1.0 + nu) * lambda * (k0 * growth + 2.0);
        // T_n/R_n² = (16Λ/R₀)(K₀D(D/Λ)ⁿ + 2Λ⁻ⁿ), small even when T_n and R_n overflow.
        let ratio = |n: i32| 16.0 * lambda / r0 * (k0 * growth * (growth / lambda).powi(n) + 2.0 * lambda.powi(-n));
        // Check that form against the definitions while they fit: R_n = R₀Λⁿ, S_n = 16R_n,
        // K_n = K₀Dⁿ, T_n = S_{n+1}(K_{n+1} + 2), t_{n+1} = t_n + T_n; also R_{n+1} ≥ 2R_n, t_n ≥ n + 1.
        let mut t = 1.0;
        for n in 0..20 {
            let rn = r0 * lambda.powi(n);
            let tn = 16.0 * r0 * lambda.powi(n + 1) * (k0 * growth.powi(n + 1) + 2.0);
            assert!((tn / (rn * rn) / ratio(n) - 1.0).abs() < 1e-12 && rn >= 1.0 && tn >= 1.0 && t >= n as f64 + 1.0);
            assert!(r0 * lambda.powi(n + 1) >= 2.0 * rn);
            t += tn;
        }
        let sum: f64 = (0..400).map(ratio).sum();
        let delta = nu * c_star * (1.0 / (r0 * r0) + sum);
        // The paper's halves: νC∗ΣT_n/R_n² ≤ ν/(48(1 + ν)) and νC∗/R₀² ≤ ν/(1024²C∗(1 + ν)²).
        assert!(nu * c_star * sum <= nu / (48.0 * (1.0 + nu)) * (1.0 + 1e-12));
        assert!(nu * c_star / (r0 * r0) <= nu / (1024f64.powi(2) * c_star * (1.0 + nu).powi(2)) * (1.0 + 1e-12));
        assert!(delta < 1.0 / 24.0, "eq. 39");
        println!("  ν = {nu:e}: R₀ = {r0:.3e}, δ = {delta:.3e} < 1/24 (bound ν/48(1+ν) = {:.3e})", nu / (48.0 * (1.0 + nu)));
    }
}

/// The torus mechanism (§3) at unit viscosity on a 64² unit torus, with a uniform translation
/// for the planar processor: `shift` per period, smoothly in phase, zero near integer phases.
/// Returns the largest `w` seen in the detector band `E` and its blocks.
fn torus_run(shift: [f64; 2], blocks: usize, band: (f64, f64), c: f64) -> Vec<f64> {
    let torus = Torus { n: [N, N], length: 1.0 };
    let p = [0.25, 0.25];
    // A translation has DΨ = I and D²Ψ = 0, so eq. 10's Hessian bound is D = 2C/width² and
    // δ = 1/(100(1 + D)).
    let width = WIDTH;
    let hessian_bound = 2.0 * c / (width * width);
    let delta = 1.0 / (100.0 * (1.0 + hessian_bound));
    let wrap = |x: f64| x - x.round();
    let profile = |centre: [f64; 2]| torus.sample(|[x, y]| bump([wrap(x - centre[0]) / width, wrap(y - centre[1]) / width]));
    let g = profile(p);
    let in_band: Vec<bool> = (0..torus.size()).map(|i| (band.0..band.1).contains(&torus.position(i)[1])).collect();
    let detector = |w: &[f64]| w.iter().zip(&in_band).filter(|(_, b)| **b).map(|(v, _)| *v).fold(f64::MIN, f64::max);
    let zero = || [vec![0.0; torus.size()], vec![0.0; torus.size()]];
    // Θ(r) = θ(2r − ½) moves the phase between ¼ and ¾ of each period.
    let phase_step = |r: f64| smooth_step(2.0 * r - 0.5);
    let phase_rate = |r: f64| 2.0 * smooth_step_derivatives(2.0 * r - 0.5)[1];
    // Eq. 18: ‖w_new − w₀‖∞ ≤ 2δD < 1/50 through injection and stirring, for the new part
    // w_new = w − w_old. The grid holds only the 2/3-rule truncation P_J w₀ of the reference, so
    // the bound is checked against P_J w₀; the truncation's own error ‖P_J w₀ − w₀‖∞ is reported.
    let max_diff = |a: &[f64], b: &[f64]| a.iter().zip(b).fold(0f64, |m, (x, y)| m.max((x - y).abs()));
    let check = |w: &[f64], old: &[f64], reference: &[f64], worst: &mut [f64; 2]| {
        let new: Vec<f64> = w.iter().zip(old).map(|(a, b)| a - b).collect();
        let truncated = torus.truncate(reference);
        let (galerkin, grid) = (max_diff(&new, &truncated), max_diff(&truncated, reference));
        assert!(galerkin <= 2.0 * delta * hessian_bound, "eq. 18: {galerkin}");
        *worst = [worst[0].max(galerkin), worst[1].max(grid)];
    };
    let (mut w, mut mass, mut seen, mut worst) = (vec![0.0; torus.size()], 0.0, Vec::new(), [0.0; 2]);
    for n in 1..=blocks {
        let t_n = 2.0 * n as f64;
        // Before the block, the old mass is spread out: w ≤ K(1, 0)·mass (eq. 15).
        let bound = torus_heat_kernel(1.0, [0.0, 0.0]) * mass;
        assert!(w.iter().all(|&v| v >= -1e-12 && v <= bound + 1e-12), "eq. 15");
        let old = w.clone();
        // Injection: h = δ⁻¹θ′((t − t_n)/δ)g over [t_n, t_n + δ].
        w = torus.advection_diffusion(&w, 1.0, (t_n, t_n + delta), 32, |t| (zero(), {
            let rate = smooth_step_derivatives((t - t_n) / delta)[1] / delta;
            g.iter().map(|v| v * rate).collect()
        }));
        check(&w, &old, &g, &mut worst);
        // Stirring: n periods in time δ, a = (n/δ)ΣΘ′(s − j)·shift with s = n(t − t_n − δ)/δ.
        let speed = |t: f64| {
            let s = n as f64 * (t - t_n - delta) / delta;
            (0..n).map(|j| phase_rate(s - j as f64)).sum::<f64>() * n as f64 / delta
        };
        let samples = 8 * n;
        let k_max = core::f64::consts::TAU * (N / 3) as f64;
        for k in 0..samples {
            let (t0, t1) = (t_n + delta * (1.0 + k as f64 / samples as f64), t_n + delta * (1.0 + (k + 1) as f64 / samples as f64));
            // Steps that keep |a|·k_max·dt ≤ 1 over this eighth of a period (none needed where
            // the phase collars stop the flow).
            let fastest = (0..=16).map(|j| speed(t0 + (t1 - t0) * j as f64 / 16.0)).fold(0.0, f64::max);
            let steps = (fastest * shift[0].hypot(shift[1]) * k_max * (t1 - t0)).ceil() as usize;
            w = torus.advection_diffusion(&w, 1.0, (t0, t1), steps.max(1), |t| ([vec![shift[0] * speed(t); N * N], vec![shift[1] * speed(t); N * N]], vec![0.0; N * N]));
            // The inviscid reference w₀ = g∘Ψ⁻¹ is the bump moved by shift·ΣΘ(s − j).
            let s = n as f64 * (t1 - t_n - delta) / delta;
            let moved: f64 = (0..n).map(|j| phase_step(s - j as f64)).sum();
            let reference = profile([p[0] + moved * shift[0], p[1] + moved * shift[1]]);
            check(&w, &old, &reference, &mut worst);
            seen.push((n, detector(&w)));
        }
        // The heat-only wait to t_{n+1}, sampled; the integrating factor makes each step exact.
        mass += torus.mean(&g);
        for k in 0..4 {
            let (t0, t1) = (t_n + 2.0 * delta + k as f64 * (2.0 - 2.0 * delta) / 4.0, t_n + 2.0 * delta + (k + 1) as f64 * (2.0 - 2.0 * delta) / 4.0);
            w = torus.advection_diffusion(&w, 1.0, (t0, t1), 1, |_| (zero(), vec![0.0; N * N]));
            seen.push((n, detector(&w)));
        }
        // The mass is Σ∫g, to the time quadrature of the injection (32 steps across θ′).
        assert!((torus.mean(&w) / mass - 1.0).abs() < 1e-6, "mass {} {mass}", torus.mean(&w));
    }
    println!("    δ = {delta:.2e}; ‖w_new − P_J w₀‖∞ ≤ {:.2e} (eq. 18: 2δD = {:.3}); ‖P_J w₀ − w₀‖∞ = {:.2e}; mass {mass:.5} = {blocks} × 2.25·width²",
        worst[0], 2.0 * delta * hessian_bound, worst[1]);
    (1..=blocks).map(|n| seen.iter().filter(|s| s.0 == n).map(|s| s.1).fold(f64::MIN, f64::max)).collect()
}

fn torus_mechanism(c: f64) {
    println!("\nTorus mechanism (§3), relaxed: {N}² unit torus, ν = 1, bump half-width {WIDTH} at p = (¼, ¼),");
    println!("  the processor a uniform translation, detector band E = 47/64 < y < 49/64:");
    let band = (47.0 / 64.0, 49.0 / 64.0);
    // "Halting": ¼ up per period puts Ψ₂(p) = (¼, ¾) in E. Block 1 runs one period; block 2
    // runs two, and the fresh bump crosses ½ (the proof's > 15/16).
    let halting = torus_run([0.0, 0.25], 2, band, c);
    println!("  orbit into E, max u₃ in E per block: {:.4}, {:.4}", halting[0], halting[1]);
    assert!(halting[1] > 15.0 / 16.0);
    // "Not halting": ¼ sideways per period keeps the bump's support (y ≤ ½) d = 15/64 from E.
    // Eq. 20's bound with the relaxed mass M = 2.25·width² per block: old mass ≤ K(1, 0)M, the
    // fresh bump ≤ H_d M with H_d = sup K(τ ≤ 1, dist ≥ d) ≥ K(1, 0), and the burst error ≤ 1/50.
    // The relaxed bumps do not shrink from block to block as the paper's bₙ do, so the old mass
    // grows by M each block and the bound stays below ½ for two blocks only.
    let blocks = 2;
    let running = torus_run([0.25, 0.0], blocks, band, c);
    let d = (band.0 - 0.25 - WIDTH).min(1.0 + 0.25 - WIDTH - band.1);
    let h_d = (1..400)
        .map(|i| torus_heat_kernel(i as f64 / 400.0, [d, 0.0]).max(torus_heat_kernel(i as f64 / 400.0, [d / 2f64.sqrt(), d / 2f64.sqrt()])))
        .fold(0.0, f64::max);
    let mass = 2.25 * WIDTH * WIDTH;
    let bound = torus_heat_kernel(1.0, [0.0, 0.0]) * (blocks - 1) as f64 * mass + h_d * mass + 1.0 / 50.0;
    println!("  sideways orbit, max u₃ in E per block: {:.2e}, {:.2e} (eq. 20 analogue, H_d = {h_d:.3}: {bound:.3})", running[0], running[1]);
    assert!(running.iter().all(|&v| v < bound) && bound < 0.5);
}

fn main() {
    println!("What this shows, and what it does not.");
    println!("  OpenAI's preprint builds smooth forces for 3-D Navier–Stokes from rest whose solution has");
    println!("  u₃ > ½ in a fixed band exactly when a given Turing machine halts, so that event is");
    println!("  undecidable. No run can show undecidability. This checks the finite parts: the machine's");
    println!("  integer-stack encoding, the scale schedules and their inequalities with the paper's");
    println!("  constants, and the scalar mechanism on a grid with relaxed scales.\n");
    encoding();
    let c = bump_constant();
    torus_schedule(c);
    array_schedule();
    torus_mechanism(c);
    println!("\nThe paper's own scales (b₁ ~ 10⁻⁷, δ₁ ~ 10⁻²⁰) are far beyond any grid, and the planar");
    println!("processor of its companion preprint is not built here; the halting equivalence is a theorem.");
}

#[test]
fn navier_stokes_computation() {
    main();
}
