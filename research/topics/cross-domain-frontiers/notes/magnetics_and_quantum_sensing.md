# Magnetic sensing, geomagnetism and joint gravity–magnetic methods, 2022–2026: a testable catalogue for `magnetics` and its neighbours

All recomputed numbers come from `magnetics_sensing_checks.py`, which sits next to these notes. It uses only the standard library and runs with `python3 -I`. It reads two trusted repository files: `crates/magnetics/data/WMM2025.COF` and the IGRF-14 block of `research/topics/magnetics-and-mhd/notes/geomagnetism_checks.py`. In its output, `ok` means a published number was reproduced, `ident` lines are exact identities that a Rust test can assert, and `show` lines are values with no published counterpart. **All 44 checks and identities pass.** The script runs in under a second.

Notation used for every item: **F** formula; **S** primary source; **T** published test value, with inputs and units; **C** proposed crate and module; **X** whether it links crates; **L** status label. "Already covered" means the function exists in a crate's "What is inside" table. "Already proposed" means it appears in `research/topics/magnetics-and-mhd/report.md` but has not shipped.

The following are left out on purpose:
- magnet devices, magnetars and other astrophysical magnetism, and frontier propulsion, all on hold;
- military uses, including the F-16 MagNav flights and magnetic anomaly detection (MAD) of submarines;
- the US defence recommendation for WMMHR.

The Tolles–Lawson model began as a MAD tool. Only its civil use, compensating survey and navigation aircraft, is catalogued here.

---

## 1. Quantum magnetometry and its fundamental limits, compared with atom-interferometer gravimeters

### Takeaway
The useful frontier items are four: the spin-projection formula, the NV-ensemble shot-noise and CW-ODMR formulas, the energy-resolution bound and its two 2022–2025 violations, and the white-noise averaging law. Each has a reproducible published number. The headline numbers:
- `magnetics::magnets::energy_resolution` already reproduces the BEC magnetometer's 0.075 ħ (0.0747 ħ here).
- The 2024 NV record's photon-shot-noise limit of 6.9 pT/√Hz is reproduced to 0.2% (6.886).
- The projection-noise algebra is the same as `gravimetry::atom`'s 1/(C√N), which gives a clean `physica::analogies` item.

### Cited Findings

**Q1. `magnetics::sensing::spin_projection_limit(gamma, n_atoms, t2, t)`** (new module `sensing`)
- **F:** δB ≈ (1/γ)·1/√(N·T₂·t). In the source's form, δB ≃ ħ/(g μ_B)·1/√(N τ T), ignoring factors of order unity. Here γ = gμ_B/ħ (rad s⁻¹ T⁻¹), N is the atom number, τ = T₂ is the coherence time and T is the measurement time. Foundational (2007). — [Budker & Romalis 2007, arXiv:physics/0611246, Nature Phys. doi:10.1038/nphys566](https://arxiv.org/abs/physics/0611246)
- **T:** The source puts the potassium SERF fundamental limit at "about 10⁻¹⁷ V^(−1/2) T/√Hz" (V in cm³), with alkali density "about 10¹⁴ cm⁻³" in SERF. Demonstrated SERF sensitivity "exceeds 10⁻¹⁵ T/√Hz", and the projected limit is "below 10⁻¹⁷ T/√Hz". — [Budker & Romalis 2007](https://arxiv.org/html/physics/0611246v1)
  - The formula reaches 10⁻¹⁷ T/√Hz for N = 10¹⁴ in 1 cm³ if T₂ = 3.2 ms, using the free-electron γ.
  - With the SERF slowing-down factor q = 6 (I = 3/2, low polarisation; assumed, not sourced) the required T₂ is 116 ms.
  - Testable identity: `spin_projection_limit(γ, 1e14, 1/(γ²·1e14·1e-34)) = 1e-17` exactly.
- **T (record):** Kominis et al. 2003 measured 0.54 fT/√Hz in 0.3 cm³ and projected a fundamental limit below 0.01 fT/√Hz. — [Kominis et al. 2003, Nature 422:596, doi:10.1038/nature01484](https://collaborate.princeton.edu/en/publications/a-subfemtotesla-multichannel-atomic-magnetometer/)
  - The energy resolution of that result is 330 ħ.
- **C:** `magnetics::sensing`. **X:** Yes, through the `physica::analogies` item Q6. **L:** established.

**Q2. `magnetics::sensing::nv_ramsey_sensitivity(n, tau, t2star, p, contrast, n_avg, t_init, t_read)`**
- **F (Barry et al. Eq. 9):** spin-projection limit η_sp ≈ ħ/(Δm_s g_e μ_B)·1/√(Nτ).
- **F (Eq. 14):** full Ramsey sensitivity η = ħ/(Δm_s g_e μ_B)·(1/√(Nτ))·e^{(τ/T₂*)^p}·√(1 + 1/(C² n_avg))·√((t_I + τ + t_R)/τ).
- **F (Eq. 11):** readout parameter σ_R = √(1 + 2(a+b)/(a−b)²). The best pulsed ensemble readouts are about 100× from the projection limit. — [Barry et al. 2020, Rev. Mod. Phys., arXiv:1903.08176](https://ar5iv.labs.arxiv.org/html/1903.08176)
- **T:** No worked numeric example was retrieved, so only identities can be tested:
  - At τ ≪ T₂*, C²n_avg → ∞ and t_I = t_R = 0, Eq. 14 reduces to Eq. 9.
  - γ_NV/2π = g μ_B/h = 28.03 GHz/T at g = 2.003. — [Barry et al. 2020](https://ar5iv.labs.arxiv.org/html/1903.08176)
- **C:** `magnetics::sensing`. **X:** with `gravimetry::atom` (the same 1/√N structure). **L:** established.

**Q3. `magnetics::sensing::photocurrent_shot_field(i_dc, slope, lockin)` and `cw_odmr_sensitivity(linewidth, contrast, rate)`**
- **F (photocurrent form, used by Sekiguchi et al.):** η = √(2·κ·q_e·I_dc)/S. Here S = γ_e·dĨ/dδ is the field slope in A/T, and κ = 2 for lock-in demodulation (their Eq. 2: ⟨n⟩² + 2×2q_eI_fl).
- **T:** I_dc = 6.4 mA and S = 9.3 A/T give **6.886 pT/√Hz**, against the published "photon-shot-noise-limited sensitivity of 6.9 pT/√Hz". The slope conversion 332 pA/Hz × 28.03 GHz/T = 9.307 A/T reproduces the published 9.3 A/T. — [Sekiguchi et al. 2024, Phys. Rev. Applied 21, 064010, arXiv:2309.04093](https://arxiv.org/html/2309.04093)
- **Record:** the measured 9.4 ± 0.1 pT/√Hz (5–100 Hz) is the best CW-ODMR NV-ensemble value reported. The paper also reports:
  - T₂* ≈ 2 µs, from about 500 ppm ¹³C and 5.6 ppm N;
  - an excitation volume of 0.03 mm³ and a 3% contrast;
  - 0.48 MHz fitted linewidth, 0.39 W of light and a 0.3 pT floor from the Allan deviation.

  — [Sekiguchi et al. 2024](https://arxiv.org/abs/2309.04093)
- **F (Barry et al. Eq. 16):** η_CW = (4/(3√3))·(h/(g_e μ_B))·Δν/(C√R). — [Barry et al. 2020](https://ar5iv.labs.arxiv.org/html/1903.08176)
  - With Sekiguchi's Δν = 0.48 MHz, C = 3% and R = I/q_e it gives **2.20 pT/√Hz, not 6.9**. The lock-in factor √2 and the different linewidth definitions account for part of the gap.
  - **Do not test Eq. 16 against this paper.** Ship it with identity tests only.
- **Energy resolution:** 9.4 pT/√Hz in 0.03 mm³ gives E_R = 1.0×10⁷ ħ (computed).
- **C:** `magnetics::sensing`. **X:** no. **L:** formula established; record recent and unreplicated (one group).

**Q4. Energy resolution E_R = δB²·V·T/(2μ₀). Already covered: `magnetics::magnets::energy_resolution`.** Add these test vectors:
- **Below ħ with a BEC:** 72(8) fT single-shot sensitivity over V = 1091(30) µm³ for 3.5 s gives E_R = 0.075(16) ħ. — [Palacios Alvarez et al. 2022, PNAS (Aalto repository copy)](https://acris.aalto.fi/ws/portalfiles/portal/79762019/Single_domain_Bose_condensate_magnetometer_achieves_energy_resolution_per_bandwidth_below_.pdf)
  - Recomputed: **0.0747 ħ**.
- **The bound itself:** several dc technologies approach E_R = ħ but did not surpass it before 2020. Model calculations for NV, SQUID and alkali-vapour sensors all predict limits near ħ. — [Mitchell & Palacios Alvarez 2020, Rev. Mod. Phys. 92, 021001, arXiv:1905.00618](https://www.arxiv.org/pdf/1905.00618)
- **Below ħ with a levitated ferromagnet:** a hard ferromagnet levitated above a superconductor at cryogenic temperature reached about 0.064 ħ. Below 10⁻³ ħ is projected. — [Ahrens, Ji, Budker, Timberlake, Ulbricht & Vinante 2025, PRL 134, 110801](https://eprints.soton.ac.uk/506233); theory in [Vinante et al. 2021, arXiv:2104.14425](https://arxiv.org/pdf/2104.14425)
  - The field noise and volume were not retrieved, so this is not yet a test vector.
- **Existing test, for reference:** Dang et al. 2010, 160 aT/√Hz in 0.45 cm³ → 43.5 ħ (published 44 ħ). — [magnetics README](https://crates.io/crates/magnetics)
- **L:** sub-ħ is recent and replicated, on two independent platforms. The "ħ limit" is contested as a universal bound.

**Q5. Earth-field scalar OPMs for MagNav and surveying (records, no new function)**
- **Pulsed ⁸⁷Rb multipass gradiometer:**
  - 14 fT/√Hz differential sensitivity at Earth-scale fields, which is 70 fT/cm/√Hz over a 0.2 cm baseline (recomputed 70.0);
  - common-mode rejection above 10⁴, compared against a Cramér–Rao bound that includes photon shot noise, spin noise and diffusion.

  — [Lucivero et al. 2022, Phys. Rev. Applied 18, L021001, arXiv:2112.09004](https://arxiv.org/abs/2112.09004v1); [Princeton record](https://collaborate.princeton.edu/en/publications/femtotesla-nearly-quantum-noise-limited-pulsed-gradiometer-at-ear/)
- **NV RF record:**
  - about 70 fT·s^1/2 at 0.35 MHz, using ferrite flux concentrators with ~300× gain (70 kHz–3.6 MHz);
  - "below 2 fT after 1 hour".

  The white-noise law gives 70/√3600 = **1.17 fT**, consistent with that statement. — [Silani et al. 2023, arXiv:2302.12401](https://ar5iv.labs.arxiv.org/html/2302.12401)
- **Averaging law. Already covered: `gravimetry::local::averaging_time(asd, target) = (asd/target)²`.** It applies unchanged to magnetometers:
  - 9.4 pT/√Hz needs 982 s of white-noise averaging to reach 0.3 pT. Sekiguchi reports 0.3 pT from an Allan deviation over "a few thousand seconds". The longer time is consistent with excess low-frequency noise; this is an inference.
- **L:** recent and unreplicated (each a single group).

**Q6. `physica::analogies`: projection noise in the atom gravimeter and the Ramsey magnetometer**
- **F:**
  - `gravimetry::atom`: δa = 1/(C√N)/(2 k_eff T²), double diffraction. This is already covered.
  - Magnetometer: δB = 1/(C√N)/(γ τ).
  - The two differ only in the scale factor: k_eff T² in rad per m s⁻², against γτ in rad per T.
- **T:**
  - The crate's own doc value k_eff T² = 4.12×10⁵ rad/(m s⁻²) at T = 0.16 s (Rosi et al.) recomputes to 4.123×10⁵.
  - ⁸⁷Rb F = 2 has γ/2π = 7.006 GHz/T (g_F = g_s/4; nuclear term neglected).
  - With N = 10⁶, C = 1: δa = 1.21×10⁻⁹ m s⁻² per shot at T = 0.16 s, and δB = 7.1×10⁻¹⁴ T per shot at τ = 0.32 s (computed, illustrative).
  — [gravimetry README](https://crates.io/crates/gravimetry) (atom module, Lévèque et al. 2021, arXiv:2011.03382)
- **C:** `physica::analogies`. **X:** yes, `gravimetry` and `magnetics`. **L:** established.

### Inferences
- The SERF form of the energy resolution explains why E_R ≈ ħ appears across technologies. With δB from Q1 and T₂ limited by spin-destruction collisions (1/T₂ = n σ_sd v̄), E_R = σ_sd v̄/(2μ₀ħγ²) does not depend on density. This is an inference; σ_sd and v̄ were not sourced. It would make a good doc paragraph, but not a test.
- The two sub-ħ platforms (BEC 2022, levitated ferromagnet 2025) bypass the single-spin argument. A `sensing` module should therefore present ħ as a benchmark, not a bound.
- For `physica`, the cross-crate fact worth stating is that one formula serves both sensors, `averaging_time(asd, target)`, and that both noise floors share the 1/(C√N) per-shot form.

### Gaps
- No 2022–2026 SERF record below Dang et al.'s 160 aT/√Hz (2010) was found. One search returned only a 2017 theory projection of 0.018 aT/√Hz for a K-Rb-He hybrid ([arXiv:1702.07755](https://arxiv.org/pdf/1702.07755)).
- No paper retrieved here gives a full (N, T₂, δB) triple that tests Q1 beyond order of magnitude.
- Field noise and volume for Ahrens et al. 2025: the Southampton page returned 403.
- SQUID records and the SQUID energy-resolution formula ε = S_Φ/(2L) were not researched with sources. Barry et al. state only that SQUIDs and SERF sensors reach the fT/√Hz level.
- The SERF slowing-down factor q = 6 is from memory.

---

## 2. Magnetic anomaly maps and navigation: MagNav and Tolles–Lawson compensation

### Takeaway
The full 18-coefficient Tolles–Lawson model, its calibration and its second-order correction are published in closed form (Gnadt, Wollaber & Nielsen 2022). On a synthetic calibration the script recovers all 18 coefficients to 2×10⁻⁷. The real flight data (SGL, Ottawa, 2020) and the published error figures are:
- tail-stinger baseline RMSE of 45–59 nT;
- 13.1 m DRMS position error;
- 22 m final error in Q-CTRL's 2025 quantum MagNav, up to 46× better than a strategic-grade INS.

These figures give arithmetic tests now, and dataset tests if the user downloads the challenge data.

### Cited Findings

**N1. `magnetics::aeromag::tolles_lawson_row(b_vec, b_vec_prev, dt)` → `[f64; 18]`, `fit(rows, y, ridge)`, `compensate(scalar, row, beta)`** (new module `aeromag`)
- **F, the linearised scalar:** |B_e| ≈ |B_t| − B_a·B̂_t (Eq. 12), where B̂_t = B_t/|B_t| are the direction cosines from a vector magnetometer.
- **F, the aircraft field:** B_a = a + bB_t + cḂ_t (Eq. 23), with three parts:
  - permanent a: 3 terms (Eq. 24);
  - induced b: 9 terms, but b is symmetric, so 6 are independent (Eq. 25);
  - eddy c: 9 terms (Eq. 26).

  That gives 18 coefficients β (Eq. 28).
- **F, the regression row:** δ = [B̂ (3); |B|·B̂_iB̂_j for i ≤ j (6); |B|·B̂_i·dB̂_j/dt (9)] (Eqs. 29–30).
- **F, calibration:**
  - The Earth field is removed by band-passing both sides, bpf(B_scalar) = bpf(A)β (Eq. 34), with a 0.1–0.9 Hz Butterworth passband.
  - β is solved by least squares, or by ridge regression (A_fᵀA_f + λI)⁻¹A_fᵀy (Eq. 36) with λ = 0.00025–0.025.
  - Compensation is then |B_e| = B_scalar − Aβ (Eq. 37).
  - When a good map is available, map-based calibration needs no band-pass.
- **F, second order:** |B_e| ≈ |B_t| − B_a·B̂_t + |B_a × B̂_t|²/(2|B_t|) (Eq. 40). Ignoring the second-order term "typically results in less than 0.05% error".
- **F, the aircraft vector field:** given in §4 of the paper, with β₅/2, β₆/2 and β₈/2 off the diagonal.

  — [Gnadt, Wollaber & Nielsen 2022, arXiv:2212.09899](https://arxiv.org/abs/2212.09899)
- **T, synthetic (script):**
  - **Map-based calibration recovers all 18 β to 1.9×10⁻⁷**. This is an identity test.
  - In a physical simulation the aircraft field (permanent + symmetric induced) was 130 nT rms. It caused a raw scalar error of 33.4 nT rms. After the 18-term fit the residual is 0.002 nT.
  - With the true aircraft field, the first-order error of Eq. 12 peaks at 0.32 nT, and Eq. 40 cuts it to 4.9×10⁻⁴ nT.
- **T, published flight numbers. They need the MagNav dataset; downloading it is the user's call.**
  - The Tolles–Lawson baseline gives 58.85 nT RMSE on Flight 1003 and 45.13 nT on Flight 1007.
  - The best LTC-CfC network gives 18.20 and 19.14 nT.

  — [NeurIPS 2023 slides (from search snippet, not fetched)](https://neurips.cc/media/neurips-2023/Slides/76255.pdf). A physics-informed LTC model reports up to a 64% cut in compensation RMSE. — [SandboxAQ publication page](https://pub.sandboxaq.com/publications/physics-informed-calibration-of-aeromagnetic-compensation-in-magnetic-navigation-systems-using-liquid-time-constant-networks)
- **The dataset:** collected by Sander Geophysics near Ottawa from a Cessna Grand Caravan, with:
  - five split-beam caesium OPMs, Mag 1 on the tail stinger at x = −12.01 m, y = 0, z = 1.37 m from the front seat rail;
  - four fluxgates;
  - flights 1002–1007, at 400 m and 800 m survey altitudes.

  The Mag 1 compensation "exactly matches" SGL's professional compensation; Mags 3 and 4 "would not be suitable". The core field is igrf = mag_1_dc − mag_1_igrf. — [Gnadt et al. 2022/2023, MagNav challenge, arXiv:2007.12158](https://arxiv.org/abs/2007.12158)
- **C:** `magnetics::aeromag`.
- **X:** Yes:
  - it needs `main_field` (IGRF/WMM core removal) and a band-pass filter;
  - `gravimetry::noise` already has the FFT;
  - the navigation filter can reuse `gravimetry::kalman`.
- **L:** established (1950 model; 2022 derivation). The machine-learning variants are recent and replicated on one shared dataset.

**N2. `magnetics::aeromag::position_crlb(sigma_b, grad)`, the map-matching bound**
- **F:** For one scalar measurement with noise σ_B, the Fisher information on position is J = ∇T∇Tᵀ/σ_B². Along the gradient, σ_x ≥ σ_B/|∇T|; a series of fixes adds its J's. SandboxAQ uses the CRLB as the map-quality metric and pairs it with an EKF variant that has unimodal and multimodal measurement noise. — [SandboxAQ, DASC 2024, publication page](https://pub.sandboxaq.com/publications/lower-bounds-on-magnetic-navigation-performance-as-a-function-of-magnetic-anomaly-map-quality)
- **T:** No published numeric CRLB was retrieved. As an illustration, 1 nT over a 10 nT/km gradient gives 100 m (computed).
- **C:** `magnetics::aeromag`. **X:** `geopotential` (gradient from the engine or planar grid), `gravimetry::kalman`. **L:** established estimation theory; 2024 application.

**N3. Published MagNav position errors (arithmetic tests, plus a status record)**
- **Q-CTRL "Ironstone Opal":**
  - best final position error 22 m, which is 0.006% of the flight distance;
  - up to 46× better than a velocity-aided strategic-grade INS, and at least 11× across conditions;
  - flights up to 19,000 ft; the parameters are learned online, without calibration manoeuvres;
  - the ground vehicle's bounded error was 7× lower than the INS.

  — [Muradoglu et al. 2025, arXiv:2504.08167](https://arxiv.org/abs/2504.08167)
  - The implied flight length is 22/6×10⁻⁵ = **366.7 km** (computed).
  - The magnetometers are proprietary and no sensitivity is given.
  - The aircraft was a Cessna 208B near Griffith, Australia, per secondary press, not verified. — [Q-CTRL press coverage](https://innovationaus.com/?p=66343)
- **Canciani (AFIT) civil-data result:** 13.1 m DRMS against 230 m for the unaided INS, a ratio of **17.6×** (computed). The filter was a marginalised particle filter with OPC magnetometers. — [Canciani 2016 dissertation](https://scholar.afit.edu/etd/251); [DSIAC 2020 briefing](https://dsiac.dtic.mil/wp-content/uploads/2020/03/DSIACWebinar_MagNav_Canciani-2.pdf)
  - Military follow-on flights are deliberately left out.
- **MagNav.jl** implements the marginalised particle filter for "highly non-linear, non-Gaussian measurement" map matching. — [MagNav.jl docs](https://docs.juliahub.com/MagNav/3z1LY/1.1.5/nav)
- **A 2026 preprint** reports fully in-flight calibration, with no calibration flight, at accuracy comparable to offline-trained models on the challenge data. — [arXiv:2603.08265](https://arxiv.org/abs/2603.08265) (search snippet only)
- **L:** Q-CTRL is recent and unreplicated (a preprint from a single company). Canciani's result is foundational (2016).

**N4. `geopotential::planar::upward_continue(grid, dx, dy, h)`: multiply the 2-D spectrum by e^{−|k|h}**
- **F:** F[T](k; z+h) = F[T](k; z)·e^{−|k|h}, for both gravity and magnetic grids. Blakely 1995, ch. 12 ([book DOI, not fetched](https://doi.org/10.1017/CBO9780511549816)).
- **Use in MagNav:** navigation systems ingest maps as a "stack of grids" upward-continued from the survey grid, typically by FFT. Characterising the uncertainty of the map gradient is described as critical. — [EGU26-14369](https://meetingorganizer.copernicus.org/EGU26/EGU26-14369.html)
- **T:** On a 128×128 grid at 10 m spacing, a pole anomaly of a 10 m sphere at 50 m depth, continued up 30 m, matches the analytic anomaly at 80 m depth to **6.9×10⁻⁴ rms** (script).
- **C:** `geopotential::planar` (new), so that `gravimetry` and `magnetics` share it. Alternatively, reuse `gravimetry::noise`'s FFT. **X:** yes, all three. **L:** established.

### Inferences
- `aeromag` needs no new numerics beyond an 18×18 least-squares solve and a band-pass filter. The FFT exists in `gravimetry::noise`. The crate dependency direction needs a decision: either `magnetics` depends on `gravimetry`, or the FFT moves to `geopotential`.
- For a Rust test that runs without the dataset, use the synthetic recovery identity and Eq. 40's error reduction. Keep the 58.85 and 45.13 nT RMSE numbers as `#[ignore]` tests gated on the user's download.

### Gaps
- No fetched source for the NeurIPS 2023 RMSE table; the numbers come from a search snippet.
- The eddy-current column definition is ambiguous in the PDF text: |B|·B̂_i·dB̂_j/dt, or |Ḃ|·… The script uses the MagNav.jl-style |B|·B̂_i·dB̂_j/dt. Confirm against `create_TL_A` in MagNav.jl before shipping.
- No published CRLB numbers were found.
- Q-CTRL's absolute ground error was not retrieved. A "180 m over 18 km" figure appears only in an unverified press summary.

---

## 3. Joint gravity–magnetic methods

### Takeaway
Poisson's relation turns the existing `gravimetry::local::{sphere, cylinder}` vertical gradient Γ_uu into the reduced-to-pole magnetic anomaly of the same body, exactly: ΔT_pole = (C_m M/(Gρ))·Γ_uu. The script confirms it to 10⁻²⁴ T at three offsets, and to 2.6×10⁻⁸ relative error at arbitrary inclination and declination. The same set of cross-crate items also includes:
- magnetic sphere and cylinder models, the magnetic twins of `gravimetry::local`;
- an MDR regression;
- Euler deconvolution with structural indices 1–3, all verified;
- reduction to the pole and pseudo-gravity, verified numerically to 1.1% and 0.13%;
- the cross-gradient constraint.

### Cited Findings

**J1. `magnetics::anomaly::sphere(a, depth, magnetization, inc_m, dec_m, inc_f, dec_f, x, y)` → ΔT (T)**
- **F:** m = (4/3)πa³M. B = C_m[3(m·r)r/r⁵ − m/r³], with C_m = μ₀/4π = 10⁻⁷ T m/A, and ΔT = B·F̂.
- Directions: m̂ = (cos I cos D, cos I sin D, sin I) with x north, y east, z down. Blakely 1995 ([doi, not fetched](https://doi.org/10.1017/CBO9780511549816)); [Nabighian et al. 2005, Geophysics 70(6):33ND, doi:10.1190/1.2133784](https://pubs.usgs.gov/publication/70027609)
- **T (computed; no published numeric example retrieved):**
  - At the pole: ΔT = 2C_m m/z³. With a = 5 m, z = 20 m and M = 1 A/m this is **13.09 nT**.
  - At the magnetic equator directly above: ΔT = −C_m m/z³ = **−6.545 nT**.
- **C:** `magnetics::anomaly` (new). **X:** yes, it is the twin of `gravimetry::local::sphere`. **L:** established.

**J2. `magnetics::anomaly::cylinder(a, depth, magnetization, …, x)`, a 2-D horizontal cylinder (line dipole)**
- **F:** m′ = Mπa² per unit length; B = 2C_m[2(m′·r̂)r̂ − m′]/r². At the pole, ΔT = 2C_m Mπa²(z² − x²)/(x² + z²)².
- **T:** Exact identity against `gravimetry::local::cylinder` through J3 (script: error 0 to 3×10⁻²⁴ T).
- **C:** `magnetics::anomaly`. **X:** with `gravimetry::local::cylinder`. **L:** established.

**J3. `physica::analogies::poisson(grav_gradient, m_over_rho)`, Poisson's relation**
- **F:** For a body with uniform density contrast ρ and magnetisation M, the magnetic potential is V = −(C_m M/(Gρ))·∂U/∂m̂. Hence B = (C_m M/(Gρ))·∇(m̂·∇U) and ΔT = (C_m M/(Gρ))·∂_f∂_m U.
- At the pole, with field and magnetisation vertical, **ΔT = (C_m M/(Gρ))·Γ_uu**. Γ_uu is exactly the second output of `gravimetry::local::{sphere, cylinder}`.
- The relation needs a uniform density and magnetisation contrast, but holds for any shape and location. Joint use is mostly to estimate the magnetisation-to-density ratio (MDR); Garland 1951, Chandler et al. 1981 and Chandler & Malek 1991 are the classic references. — [Mendonça 2003, SBGf (review of the method)](https://sbgf.org.br/mysbgf/eventos/expanded_abstracts/8th_CISBGf/arq_251.pdf); [analytic-signal MDR variant](https://jesphys.ut.ac.ir/article_28442.html?lang=en)
- **T (script):**
  - Sphere and cylinder at x = 0, 7 and 25 m: |ΔT − (C_m M/Gρ)Γ_uu| ≤ 3.3×10⁻²⁴ T.
  - General direction (I_m 60°, D_m 10°, I_f −35°, D_f −20°), using finite differences of U = GΔM/r: relative error 2.6×10⁻⁸.
  - The crate's own doc value for a 5 m void at 20 m depth, Γ_uu = −23.6 E, recomputes to −23.589 E. If the void also removes M = 1 A/m, the same Γ_uu predicts a 13.09 nT magnetic low.
- **C:** `physica::analogies`, or a dev-dependency test in `magnetics::anomaly`. **X:** yes, `gravimetry` and `magnetics`. **L:** established.

**J4. `magnetics::anomaly::mdr_regression(dt, gzz)`, moving-window Poisson analysis**
- **F:** Regress ΔT_pole on Γ_uu (or ΔT on the vertical derivative of gravity) within a window. The slope is C_m M/(Gρ), so MDR = slope·G/C_m (A m² kg⁻¹). The total gradient of the magnetic field is preferred because it is "almost independent on the direction of the total-magnetization". — [Tyrrhenian seamounts application (Univ. Naples record)](https://www.iris.unina.it/handle/11588/971233); [Doo et al. 2009 (NCU)](https://ir.lib.ncu.edu.tw/handle/987654321/27577)
- **T:** Synthetic: recovers M/ρ = 2.5/2700 exactly (0 error, script). No published field MDR value with inputs was retrieved.
- **C:** `magnetics::anomaly`. **X:** with `gravimetry`. **L:** established.

**J5. `geopotential::planar::euler(fx, fz, f, x, n)` → (x₀, z₀, background), with a `structural_index` table**
- **F:** Euler's homogeneity equation (x − x₀)∂f/∂x + (z − z₀)∂f/∂z = −N(f − B), solved by least squares in a window. The structural index N must be an integer to keep its physical meaning (Reid & Thurston 2014). — [SBGf abstract on structural index](https://sbgf.org.br/mysbgf/eventos/expanded_abstracts/15th_CISBGf/session/POTENTIAL%20FIELDS%20METHODS/Correct%20structural%20index%20defined%20by%20base%20level%20estimates%20in%20Euler%20deconvolution%20.pdf)
- **T (exact, script):** the indices, verified to ≤ 5×10⁻⁹ by finite differences:

  | Body | Field | N |
  |---|---|---|
  | sphere | gravity g_z | 2 |
  | horizontal cylinder | gravity g_z | 1 |
  | sphere | magnetic ΔT, at the pole and at I = 60° | 3 |
  | horizontal cylinder | magnetic ΔT | 2 |

  The inverse problem recovers a dipole's 20 m depth to 1.2×10⁻⁹ m and its 3 m offset to 7×10⁻¹¹ m.
- **C:** `geopotential::planar` (shared). **X:** `gravimetry::local`, `magnetics::anomaly`. **L:** established.

**J6. `geopotential::planar::reduce_to_pole(grid, inc_m, dec_m, inc_f, dec_f)` and `pseudo_gravity(grid, m_over_rho, …)`**
- **F:** F[ΔT_r] = F[ΔT]/(Θ_mΘ_f), with Θ = m̂_z + i(m̂_x k_x + m̂_y k_y)/|k|. Pseudo-gravity (Baranov 1957) is F[g_z] = (Gρ/(C_m M))·F[ΔT_pole]/|k|. Blakely 1995, ch. 12 (equation numbers not checked) ([doi](https://doi.org/10.1017/CBO9780511549816)).
- **T (script, 128×128 grid, I = 45°, sphere a = 10 m at 50 m):**
  - The RTP output matches the analytic pole anomaly to **1.1% rms**. The opposite sign convention in Θ gives 168%, so the test pins the sign for an FFT with e^{−ik·x}, the numpy-style convention.
  - At the pole the filter is 1 (to 10⁻¹⁶).
  - Low-latitude instability: along the declination direction the filter amplitude reaches 1/sin²I, which is **33.2 at I = 10°**.
  - Pseudo-gravity matches `gravimetry::local::sphere` g_z (peak minus corner) to **−0.13%**.
- **C:** `geopotential::planar`. **X:** `gravimetry` and `magnetics`. **L:** established. Low-latitude RTP remains an active problem (see Gaps).

**J7. `physica::joint::cross_gradient(m1, m2)`**
- **F:** t = ∇m₁ × ∇m₂, which is zero wherever the two models share structure. It is used as the coupling term in joint gravity–magnetic inversion.
- **S:** introduced by Gallardo & Meju 2003, GRL 30(13):1658, [doi:10.1029/2003GL017370](https://doi.org/10.1029/2003GL017370); JGR 2004, [doi:10.1029/2003JB002716](https://doi.org/10.1029/2003JB002716).
- **T:** exact identity t = 0 for m₂ = 3m₁ + m₁² (1.5×10⁻¹¹, script). A shifted body gives t ≠ 0 (−0.026).
- **C:** `physica`. **X:** yes. **L:** established.

### Inferences
- Poisson's relation is the strongest cross-crate test in this catalogue. It needs no published magnetic number: every magnetic body model can be tested against the already-tested gravity models in `gravimetry::local`. The relation is also an exact identity, unlike the potential-theory analogies already in `physica::analogies`.
- `gravimetry::local` returns Γ_uu already. A magnetic sphere at the pole is therefore one line: `CM * m / (G * rho) * sphere(a, z, rho, x).1`. The general-direction model still needs the dipole formula J1.
- Upward continuation, Euler, RTP and pseudo-gravity all act on planar grids of either field. Putting them in one shared `planar` module avoids writing them twice.

### Gaps
- No published numeric test value was retrieved for a magnetic sphere, cylinder or prism; Blakely's book was not fetched. The Bhattacharyya (1964) prism model, the 3-D twin of a gravity prism, was not sourced. Test it by its far-field limit (the dipole J1) and by J3.
- No 2022–2026 primary paper on MDR or Poisson analysis was retrieved. The DOIs for Chandler et al. 1981, Reid et al. 1990, Reid & Thurston 2014 and Baranov 1957 were not verified. Their values are verified here by the identities instead.
- No recent paper on low-latitude RTP methods was looked up.

---

## 4. Geomagnetism 2022–2026: jerks, secular acceleration, the South Atlantic Anomaly, core flow, the lithospheric field

### Takeaway
The best new tests come from the field models the crate already parses. IGRF-14's 2020→2025 South Atlantic Anomaly minimum:
- moves **21.0 km/yr** west (published: ~20 km/yr);
- deepens by **−30.5 nT/yr** (published: −290 nT over 2015–2024.5, i.e. −30.5 nT/yr).

The eastern lobe south-west of Africa deepens **2.9×** faster (−87.4 nT/yr) and moves only 8.4 km/yr; published: "2–3 times greater … no apparent motion". WMM2025 gives a secular-variation timescale τ_SV = 380 yr against the published 415 yr, and a Voorhies-type spectral fit gives a source radius of 3541–3589 km against the published 3512 ± 64 km. Two results have only order-of-magnitude checks, through the existing `mhd::alfven_speed`: the 7-yr magneto-Coriolis waves (2022) and the 2019–2020 jerk.

### Cited Findings

**G1. `magnetics::main_field::saa_minima(model, t)` → `Vec<(F, lat, lon)>` and `saa_drift(model_a, t_a, model_b, t_b)`. Already proposed: `minimum_intensity`; the second lobe and the drift are new.**
- **F:** A local-minimum search of F at h = 0, as in the report's pattern search, run separately from seeds near (−26°, −60°) and (−40°, 0°). Drift is the great-circle distance divided by Δt; the intensity rate is ΔF/Δt.
- **T (script, IGRF-14):**
  - **Western minimum:**
    - 2020: 22,223.8 nT at (−26.267°, −58.978°).
    - 2025: 22,071.4 nT at (−26.22°, −60.03°), the report's existing value reproduced to 0.4 nT.
    - Rates: **20.98 km/yr, −0.210°/yr in longitude, −30.49 nT/yr**.
  - **WMM2025 SV forecast, 2025→2030:** 21.05 km/yr and −29.31 nT/yr.
  - **Eastern lobe:**
    - 2025: 23,461.8 nT at (−40.70°, −0.72°).
    - Rates: **−87.4 nT/yr and 8.4 km/yr**; ratio to the western rate 2.87.
- **Published:**
  - The SAA moved west at about 20 km/yr over 1970–2020, its minimum fell from about 24,000 to 22,000 nT, and a second minimum has emerged south-west of Africa. — [ESA Swarm](https://www.esa.int/Applications/Observing_the_Earth/FutureEO/Swarm/Swarm_probes_weakening_of_Earth_s_magnetic_field)
  - From the MSCM model (Swarm + CSES + MSS-1, 2014.0–2024.5): the western minimum fell 290 nT since 2015 and moved west at 20 km/yr; the eastern minimum fell 730 nT ("2–3 times greater") with "no apparent motion". — [Gao, Wang, Livermore, Rogers & Liu, EGU25-12183](https://meetingorganizer.copernicus.org/EGU25/EGU25-12183.html)
- **C:** `magnetics::main_field`. **X:** `geopotential` engine. **L:** recent and replicated (ESA/Swarm, MSCM and IGRF-14 agree).

**G2. `magnetics::main_field::sv_timescale(model, n)` and `fit_tau_sv(model, n_range)`**
- **F:** τ_n = √(R_n/S_n), where R_n = (n+1)Σ(g² + h²) and S_n = (n+1)Σ(ġ² + ḣ²) are the Lowes spectra of the field and of its secular variation. The one-parameter law is τ_n = τ_SV/n, fitted here as τ_SV = exp(mean ln(nτ_n)) over n = 2–12.
- **S:** the law holds for the non-dipole field within statistical scatter, with τ_SV = 415 yr from recent satellite models; the error bars are garbled in the extracted abstract, likely +45/−55. — [Lhuillier, Fournier, Hulot & Aubert 2011, GRL 38, L09306, doi:10.1029/2011GL047356](https://www.geophysik.uni-muenchen.de/en/research/publications/1774)
- The scaling with degree is still debated. — [Tsang & Jones 2024, GJI](https://www.mas.ncl.ac.uk/~nyt27/paper/ggae234.pdf)
- **T (script, WMM2025):**
  - **τ_SV = 380 yr** (n = 2–12), against 415 yr published. This is within −55 if that lower bar is right.
  - τ₁ = 1124 yr, τ₂ = 158 yr, τ₆ = 70.5 yr, τ₁₂ = 17.2 yr.
- **C:** `magnetics::main_field`. It reuses the proposed `lowes_spectrum` and `geopotential::Coefficients::degree_variances` (already covered) applied to the SV coefficients. **X:** `geopotential`. **L:** established; scaling contested.

**G3. `magnetics::main_field::source_radius(model, n_range)`, a spectral depth to the core**
- **F:** Fit ln R_n = ln K + ln w(n) + (2n+4) ln(c/a), so that c = a·exp(slope/2). The core-source spectral shape involves K(n + ½) times (c/a)^(2n+4) (Voorhies 2004; exact weight not confirmed).
- **S:** Voorhies, Sabaka & Purucker 2002 found a core radius of 3512 ± 64 km (seismic 3480 km), and a crustal shell at 6367 ± 14 km. — [NASA NTRS](https://ntrs.nasa.gov/citations/20020046450)
- **T (script, WMM2025, n = 2–12):** the fitted radius depends strongly on the weight w(n):

  | Weight w(n) | Radius (km) |
  |---|---|
  | 1 | 3322 |
  | n + ½ | 3079 |
  | (n + ½)/(n(n+1)) | **3589** (3541 for n = 3–12) |

  Only the last is within ~1σ of 3512 ± 64.
- **C:** `magnetics::main_field`. **X:** `geopotential`, and the `physica::analogies` continuation item. **L:** established. Read Voorhies's paper for the exact weight before writing the test.

**G4. `magnetics::main_field::jerk_fit(t, sv)` → (t_jerk, ΔSA), a two-segment linear fit of the secular variation**
- **F:** A jerk is a V-shaped change in a secular-variation series, i.e. a step in secular acceleration. Find the break t_j that minimises the residual of a two-segment fit; ΔSA is the difference of the two slopes (nT/yr²).
- **S:** A new jerk around 2019–2020 is seen at observatories and in Swarm virtual observatories, and CHAOS shows it is global. A pulse in global-average secular acceleration at the core–mantle boundary appeared around mid-2018. — [Pavón-Carrasco et al. 2021, EPS 73:175, doi:10.1186/s40623-021-01504-2](https://www.doi.org/10.1186/s40623-021-01504-2)
- From CHAOS-7.2, minima in the acceleration of the SAA's area line up with jerks, and a radial-field SA pulse has its maximum in 2016.2 and minimum in 2017.5. — [Frontiers in Earth Science 2020, doi:10.3389/feart.2020.607049](https://www.frontiersin.org/journals/earth-science/articles/10.3389/feart.2020.607049/full)
- **T:** no observatory ΔSA amplitude was retrieved (the Springer page needs a cookie login). Ship with a synthetic V-shaped series (exact recovery) only.
- **C:** `magnetics::main_field` or a small `secular` module. **X:** can reuse `gravimetry::kalman` for gappy series. **L:** recent; one primary paper retrieved, with several data types.

**G5. Core waves. Existing `magnetics::mhd::alfven_speed` with new test vectors.**
- **Magneto-Coriolis waves:**
  - Quasi-geostrophic magneto-Coriolis waves with a period of ~7 yr travel westward along the core–mantle boundary at up to ~1500 km/yr, with flow amplitude ~3 km/yr. They are seen in 1999–2021 satellite data (Ørsted, CHAMP, Swarm). — [Gillet, Gerick, Jault, Schwaiger, Aubert & Istas 2022, PNAS 119, e2115258119 (ESA summary)](https://www.esa.int/Applications/Observing_the_Earth/FutureEO/Swarm/Swarm_unveils_magnetic_waves_deep_down); [IPGP](https://www.ipgp.fr/~aubert/a-new-category-of-waves-identified-in-earth-s-core.html)
  - No axisymmetric background field tested supports modes consistent with recent observations; some non-axisymmetric fields do. — [Gerick & Livermore 2024, arXiv:2403.03011](https://arxiv.org/abs/2403.03011)
- **Core flow:** the planetary gyre flows at roughly 15–50 km/yr, with ~7-yr waves at low latitudes. — [Finlay, Gillet, Aubert, Livermore & Jault 2023, Nat. Rev. Earth Environ. 4:377, doi:10.1038/s43017-023-00425-w](https://orbit.dtu.dk/en/publications/gyres-jets-and-waves-in-the-earths-core/)
- **Torsional waves:** inverted for 6–8 yr periods (Gillet et al. 2010, Nature 465:74, doi:10.1038/nature09010). — [IDEAS record](https://ideas.repec.org/a/nat/nature/v465y2010i7294d10.1038_nature09010.html)
  - A lecture summary gives "~6 yrs traveltime: Bs ≳ 2 mT". — [GdR slides](https://mhd.ens.fr/GdR/GDR_Vercors.pdf)
  - Maffei & Jackson give propagation from the inner core to the CMB in about 4 yr. — [2013 abstract](https://mhd.ens.fr/GDRE/GDRE2013/ABS/abstract_Maffei.pdf)
- **T (script):** with an assumed mean ρ = 1.1×10⁴ kg/m³:
  - B_s = 2 mT → V_A = 537 km/yr → **4.2 yr** across 1220–3480 km. This is consistent with "about 4 yr" but not with the slide's 6 yr.
  - 4 mT → 2.1 yr.
  - The magneto-Coriolis wavelength is 1500 km/yr × 7 yr = 10,500 km, i.e. azimuthal wavenumber ≈ 2.1 at the CMB (inference).
- **C:** test vectors on `magnetics::mhd::alfven_speed`. **X:** `mhd` and `main_field`. **L:**
  - torsional-wave field strength: contested (2 against 4 mT, 4 against 6 yr);
  - magneto-Coriolis waves: recent, with their interpretation debated.

**G6. The lithospheric field: WMMHR2025 and LCS-1 through the existing engine (no new numerics)**
- **WMMHR2025:**
  - core field and SV for n = 1–15, and crust for n = 16–133 from MF7 (last two years of CHAMP);
  - ~300 km resolution at the equator, against 3300 km for WMM;
  - file `WMMHR.COF`, released 17 Dec 2024, valid to 31 Dec 2029. Test values ship inside the coefficient zip, behind NCEI's survey page.

  — [NCEI WMMHR page](https://www.ncei.noaa.gov/products/world-magnetic-model-high-resolution); [CIRES, DOI 10.25921/qb1c-vn52](https://geomag.colorado.edu/wmm2025-and-wmmhr2025)
- **LCS-1:**
  - built from CHAMP 2006–2010 and Swarm 2014–2016 gradient data;
  - degree 14–185, on a 0.25° grid;
  - GJI 211:1461 (2017), so foundational. — [search summary; LCS-1 use in UKY study](https://uknowledge.uky.edu/ees_facpub/11)
- **C:** `magnetics::main_field::read_cof`, already covered for degree 12. It needs degree up to 133, and `geopotential::Legendre`, which is stable to degree 2700 (already covered). **X:** `geopotential`.
- **L:** established. This is already listed under "What stays out" in the magnetics report until test values are obtained. The test values exist in the NCEI zip, and downloading them is the user's call.

### Inferences
- G1 is the best new geomagnetic test. Two independent published figures (20 km/yr and −290 nT over ~9.5 yr) are reproduced to about 5% and 0.1% by IGRF-14 coefficients that are already embedded in the repository. The eastern-lobe ratio of 2.87 falls inside the published "2–3 times". The time base "since 2015", taken as 9.5 yr to 2024.5, is my reading of the abstract.
- G2 and G3 turn the proposed `lowes_spectrum` into two physically meaningful numbers, the SV timescale and the source depth, with published comparison values. They would make good README examples for `main_field`.
- Swarm's own core-field products (CHAOS-7/8) were not needed. IGRF-14 and WMM2025 carry enough secular variation for these tests.

### Gaps
- The exact τ_SV error bars (Lhuillier 2011) and Voorhies's exact spectral weight need reading the papers.
- Observatory jerk amplitudes (ΔSA, nT/yr²) for 2014, 2017 and 2020 were not retrieved. The source radius and density for the torsional-wave test are assumed values.
- CHAOS-7 or CHAOS-8 numerical point values were not retrieved, and neither were Finlay et al. 2020's specific numbers. No 2022–2026 lithospheric degree spectrum in nT² was retrieved.
- The 4 mT core-field figure attributed to Gillet et al. 2010 was not confirmed.

---

## 5. Superconducting and high-field magnets 2022–2026

### Takeaway
`magnetics::magnets` lacks four simple engineering formulas: hoop stress σ = BJR, the no-insulation (NI) coil time constant L/R_c, the bypass current V/R_c, and a thick-solenoid centre constant. Hahn et al.'s 45.5 T paper (2019, foundational) supplies published values for all four, reproduced here:
- 1.070 s against 1.07 s;
- 0.2006 A and 11.98 mT against 0.2 A and 12 mT;
- 45.74 T against 45.7 T.

The 2024–2026 records come only from press releases:
- resistive: 42.02 T at 32.3 MW (CHMFL, Sept 2024);
- all-superconducting user magnet: 35.6 T in a 35 mm bore (SECUF, Jan 2026); a 35.1 T result held for 30 min (ASIPP, Sept 2025).

### Cited Findings

**M1. `magnetics::magnets::hoop_stress(b, j, r)` = B·J·R.** Distinct from the existing `lorentz_line_load` (I·B, per unit length).
- **S:** LBC3 was "dry-wound" with self-supporting turns, "with a hoop stress of σθ = BJR"; peak magnetic stress was estimated at 691 MPa and total strain at 0.38%. — [Hahn et al. 2019, Nature 570:496, doi:10.1038/s41586-019-1293-1](https://ideas.repec.org/a/nat/nature/v570y2019i7762d10.1038_s41586-019-1293-1.html) (text extracted from the [NHMFL copy](https://public.magnet.fsu.edu/AToth/NSF%20publication%20access/2019/09%20Newly%20Entered%20Pubs%20September%202019/2019%20cited/2019_45.5-tesla%20direct-current%20magnetic%20field.pdf))
- **T:** J = 1420 A/mm² (published; recomputed 245.3 A/(4.02 × 0.043 mm²) = 1419.1).
  - BJR at the inner radius (45.5 T, 7 mm): 452 MPa.
  - BJR at the outer radius (31.1 T background, 17 mm): 751 MPa.
  - The published 691 MPa peak lies inside this envelope (script). Exact reproduction needs B(R) across the winding.
- **C:** `magnetics::magnets`. **X:** with `plasmoid` (magnet stress limits; the existing `virial_specific_energy`). **L:** established.

**M2. `magnetics::magnets::ni_time_constant(l, r_c)` = L/R_c and `ni_bypass_current(v, r_c)` = V/R_c**
- **T (all published, Hahn et al. 2019 Table 1 and text):**
  - **L = 50.4 mH, R_c = 47.1 mΩ → 1.0701 s** (published 1.07 s).
  - Coil voltage **9.45 mV → leakage 0.2006 A** (published 0.2 A).
  - × magnet constant 59.7 mT/A → **11.98 mT** (published "centre field of 12 mT").
- **C:** `magnetics::magnets`. **X:** no. **L:** established. The NI technique is recent and replicated: CFS's no-insulation, no-twist (NINT) coils also use it. — [SPARC TFMC special section](https://snf.ieeecsc.org/post/announcement/tas-special-section-sparc-tfmc-program)

**M3. Screening-current field deficit, as a published number (no model)**
- **T:** The design constant gives a calculated field of 31.1 T + 59.7 mT/A × 245.3 A = **45.74 T** (published "45.7 T"). The measured field was 45.5 T. The 0.2 T deficit is attributed to "diamagnetic screening currents". — [Hahn et al. 2019](https://ideas.repec.org/a/nat/nature/v570y2019i7762d10.1038_s41586-019-1293-1.html)
- **C:** a doc example on `magnets::solenoid_bz`. **L:** established.

**M4. `magnetics::magnets::thick_solenoid_center(n, i, a1, a2, half_length)`. Already proposed as `fabry_factor` in the report.**
- **F:** B₀ = μ₀ N I/(2(a₂ − a₁))·ln[(a₂ + √(a₂² + b²))/(a₁ + √(a₁² + b²))], assuming uniform winding current density.
- **T:** LBC3 has N = 2717 (12 pancakes × 226.4 turns, recomputed 2716.8), a₁ = 7 mm, a₂ = 17 mm and b = 26.55 mm. This gives **58.45 mT/A against the published design constant of 59.7 mT/A (−2.1%)**.
- **Inductance check:** Wheeler's multilayer formula gives 51.69 mH against the published 50.4 mH (+2.6%).
- **Recommended tolerance:** 3%. The report's fabry_factor still had no test value; this supplies one.
- **C:** `magnetics::magnets`. **L:** established.

**M5. The record table (doc constants, arithmetic tests only)**
- **45.5 T, test coil, 2019:** LBC3 added 14.4 T inside a 31.1 T, 18 MW resistive magnet with a 50 mm warm bore. Earlier LBCs reached 40 T and 42.5 T. The prior 45 T hybrid (since 2000) uses a 31 MW, 33.6 T resistive insert inside 11.4 T LTS coils. Quench propagation in the NI coil was 1.75 m/s measured and 1.46 m/s simulated. — [Hahn et al. 2019](https://ideas.repec.org/a/nat/nature/v570y2019i7762d10.1038_s41586-019-1293-1.html)
- **Resistive record:** 42.02 T at 32.3 MW, CHMFL Hefei, 22 Sept 2024, beating NHMFL's 41.4 T (2017). — [CAS newsroom](https://english.cas.cn/newsroom/research_news/phys/202409/t20240923_690300.shtml); [HFIPS](https://english.hf.cas.cn/nr/headline/202409/t20240922_690291.html)
  - Computed: +0.62 T, and B²/P = 54.7 T²/MW.
- **All-superconducting user magnet:** 35.6 T in a 35 mm usable bore, SECUF Beijing, 24 Jan 2026. It breaks NHMFL's 32 T record (+11.25%, computed). — [IOP CAS](https://english.iop.cas.cn/news/202604/t20260408_1155395.html)
- **All-superconducting, ASIPP:** 35.1 T held for 30 min, 27 Sept 2025. Some outlets call it a world record, others a national record. — [TechNode](https://technode.com/2025/09/29/chinas-all-superconducting-magnet-reaches-field-700000-times-stronger-than-earths/)
- **CFS/MIT TFMC (2021; papers 2024):** 20.1 T and 40.5 kA, already covered by `lorentz_line_load`'s 815 kN/m test. — [IEEE CSC](https://snf.ieeecsc.org/post/announcement/tas-special-section-sparc-tfmc-program)
- **L:**
  - 45.5 T: established (peer-reviewed);
  - 42.02 T: recent, facility press release;
  - 35.6 T and 35.1 T: recent and unreplicated, press releases with no peer-reviewed numbers retrieved.

### Inferences
- M1, M2 and M4 together let one LBC3 example in `magnets` reproduce a whole published magnet: the field, current density, stress envelope, time constant and leakage. This mirrors how the CFS coil is already used.
- Records belong in docs with their category (test coil, user magnet, resistive, hybrid). The report's existing "record conflicts" row (35.1 against 35.5 T) is now superseded by 35.6 T (SECUF) and 35.1 T (ASIPP), which are different magnets.

### Gaps
- The CHMFL 45.22 T hybrid record (2022) was not verified in this session.
- No source with inputs was found for a simple screening-current model, such as Brandt–Indenbom thin-strip magnetisation (M = J_c w/4).
- No inputs were found for the adiabatic normal-zone propagation formula v = (J/C)√(ρk/(T_s − T₀)), for an LTS or HTS coil.
- The peak stress cannot be reproduced exactly without the field map B(R).
- The SECUF and ASIPP magnets lack published LTS/HTS splits, currents and stresses.
- CFS TFMC stored energy and quench data are in the open-access IEEE TAS papers, which were not read.
