# Light in plasmas and magnetic fields: laser–plasma physics and magneto-optics (2022–2026)

Every number marked "recomputed" was rerun from its stated inputs, using standard-library Python
(`python3 -I light_in_plasma_checks.py`, in this directory, CODATA 2018/2022 constants). All its
asserts pass. Crate dependencies matter for where an item can go: `plasmoid` and `gravitational`
depend only on `metrology`, `magnetics` on `geopotential` and `metrology`, `gravitics` also on
`gravitational`, and only `physica` sees every crate (read from each `Cargo.toml`). So a
"links two crates" item either takes the other crate's output as a plain argument, or lives as a
cross-check test or example in `physica`.

Already in the crates (from the "What is inside" tables and the source):
`gravitics::strong_field::{a0, chi_e, peak_field, critical_intensity, birefringence,
euler_heisenberg_fraction}` (with the `Intensity::{Linear, Peak}` conventions);
`gravitational::hf::{plasma_frequency, critical_density, plasma_index}`;
`gravitational::propagation::{speed_deficit, massive_graviton_delay}`;
`plasmoid::aneutronic::{synchrotron, bremsstrahlung}`; `plasmoid::pinch::cyclotron_energy`;
`plasmoid::tokamak::ITER`; `magnetics::mhd::{inertial_length, gyroradius}`;
`magnetics::magnets::{flux_compression, larmor_frequency, skin_depth}`. Nothing in any crate
covers ponderomotive energy, wakefield scalings, ion acceleration, Faraday rotation, rotation
measure, dispersion measure, Verdet constants or Cotton–Mouton effects (grep for `faraday`,
`verdet`, `ponderomot`, `wakefield`, `rotation measure` found nothing).

## Laser–plasma basics with recent benchmarks (ponderomotive force, self-focusing, LWFA scalings, channel guiding, inverse bremsstrahlung)

### Takeaway
Lu et al.'s 2007 blowout-regime scalings are short, closed-form and testable. Recomputed from
the papers' own inputs, Eq. 6 gives 7.6 GeV against Gonsalves 2019's measured 7.8 GeV. For
Picksley 2024 it gives a 13.2 GeV ceiling against 9.2 GeV measured, and the gap closes to about
7.5–10 GeV once the stage is cut at the pump-depletion length and at the 30 cm channel. Lu's
matched-spot rule (k_pR = 2√a₀) gives the 52 µm focal spot BELLA actually used. These belong in
a new `plasmoid::laser` module that takes `a₀` from `gravitics::strong_field::a0` as an argument.

### Cited Findings
**L1. Relativistic self-focusing critical power** (established; foundational, Sprangle et al. 1987).
- Formula: `P_c = 8π ε₀ m_e² c⁵ / e² · (ω₀/ω_p)² = 17.4 (ω₀/ω_p)² GW = 17.4 (n_c/n_e) GW`, with `n_c = ε₀ m_e ω₀²/e²`. Lu et al. use `P_c = 17 ω₀²/ω_p² [GW]` — [Lu et al. 2007, PRSTAB 10, 061301, arXiv:physics/0612227](https://arxiv.org/abs/physics/0612227).
- Recomputed: the coefficient `8π ε₀ m_e² c⁵/e²` = 17.42 GW, and `n_c(1 µm)` = 1.115×10²¹ cm⁻³ — [light_in_plasma_checks.py](light_in_plasma_checks.py).
- Crate: `plasmoid::laser::critical_power(n_e, wavelength)`. Links: `gravitational::hf::critical_density` gives the same `n_c` from a frequency (already covered: `critical_density`), but `plasmoid` cannot import `gravitational`, so either duplicate the one-liner or add a `physica` test that the two agree.

**L2. Lu et al. 2007 blowout-regime design equations** (established 2007; recent benchmarks replicated at BELLA in 2019 and 2024). Equation numbers are from the arXiv text — [Lu et al. 2007](https://arxiv.org/abs/physics/0612227):
- Eq. 1, matched spot: `k_p R ≃ k_p w₀ = 2√a₀`. Eq. 2: `a₀ ≃ 2 (P/P_c)^{1/3}`.
- Eq. 3, pump depletion (etching): `L_etch ≃ (c/v_etch) cτ_FWHM ≃ (ω₀²/ω_p²) cτ_FWHM`, with `v_etch ≃ c ω_p²/ω₀²`.
- Eq. 4, dephasing: `L_d ≃ c/(c − v_φ) R ≃ (2/3)(ω₀²/ω_p²) R`, with `v_φ ≃ c[1 − 3ω_p²/(2ω₀²)]`.
- Eq. 6, energy gain: `ΔE ≃ (2/3) m_e c² (ω₀/ω_p)² a₀`, i.e. `ΔE[GeV] ≃ 1.7 (P[TW]/100)^{1/3} (10¹⁸/n_p[cm⁻³])^{2/3} (0.8/λ₀[µm])^{4/3}`. Eq. 7: `ΔE[GeV] ≃ 3.8 (P/P_c)^{−2/3} (P[TW]/100)`.
- Eq. 10, charge: `N ≃ (8/15)/(k₀ r_e) · √(P/(m²c⁵/e²)) ≃ 2.5×10⁹ (λ₀[µm]/0.8) √(P[TW]/100)`.
- The condition for reaching dephasing is `L_etch > L_d ⇒ cτ_FWHM > 2R/3`. The paper's own test case is 200 TW, 30 fs, 0.8 µm, n = 1.5×10¹⁸ cm⁻³: 1.5 GeV in PIC against "≃1.6 GeV" from Eq. 6, with `L_d ≃ 1.31 cm > 0.96 cm ≃ L_pd` — [Lu et al. 2007](https://arxiv.org/abs/physics/0612227).
- Recomputed from those inputs: Eq. 6 gives 1.63 GeV (matches). `L_d` = 1.40 cm against the printed 1.31 (7% high) and `L_pd` = 1.04 cm against 0.96 (9% high), probably from rounding of `R` or `τ` in the paper. Eq. 10 at 100 TW gives N = 2.58×10⁹ against the printed 2.5×10⁹. The exact form `(2/3)mc²(ω₀/ω_p)²a₀` with a₀ from Eq. 2 gives 1.78 GeV where the practical form gives 1.7 — [light_in_plasma_checks.py](light_in_plasma_checks.py).
- Crate: `plasmoid::laser::{matched_a0, matched_spot, dephasing_length, depletion_length, lu_energy_gain, lu_charge}`. Links: `gravitics::strong_field::a0` (the a₀ from an intensity), `gravitational::hf::plasma_frequency`.

**L3. Benchmark: BELLA 2024, 9.2 GeV in a 30 cm hydrogen HOFI channel** (recent; same group's 2019 record replicated and extended).
- Published inputs: 500 TW, 21.3 ± 3 J, ~40 fs FWHM, λ = 815 nm, w₀ = 53 ± 1 µm, peak a₀ ≈ 2.2, n₀ ≈ 1×10¹⁷ cm⁻³, 30 cm channel, low-power matched spot w_m ≈ 37 µm. Outputs: quasi-monoenergetic peaks up to 9.2 GeV with spreads as low as 3%, charge extending beyond 10 GeV, average field ≈ 30 GV/m. Simulations of an optimised pulse reach 13.0–13.1 GeV — [Picksley et al. 2024, PRL 133, 255001, arXiv:2408.00740](https://arxiv.org/html/2408.00740v1); DOI [10.1103/PhysRevLett.133.255001](https://link.aps.org/accepted/10.1103/PhysRevLett.133.255001). Physics Today notes that 9.2 GeV is not the wakefield record — [Physics Today](https://physicstoday.aip.org/news/plasma-channel-guides-electrons-to-10-gev).
- Recomputed with Lu: matched spot `2√a₀/k_p` = 52.2 µm against the 53 µm used. a₀ from Eq. 2 is 2.41, and from 500 TW focused to w₀ = 53 µm (`a₀ = 0.855 λ[µm] √(I/10¹⁸ W cm⁻²)`) it is 2.35, against ≈2.2 reported. Eq. 6 gives 13.2 GeV, which is the paper's own simulated optimum. With `L_d` = 58 cm, `L_pd` = 20 cm and the channel at 30 cm, the dephasing condition `cτ > 2R/3` fails (12 µm < 35 µm) — [light_in_plasma_checks.py](light_in_plasma_checks.py).
- Testable statement: `lu_energy_gain(500 TW, 1e17, 0.815 µm)` = 13.2 GeV (±5%), and `matched_spot` = 52 µm against 53 µm.

**L4. Benchmark: BELLA 2019, 7.8 GeV in a 20 cm laser-heated capillary** (recent before 2022; replicated by L3).
- Published: peak power 0.85 PW guided over 15 diffraction lengths (20 cm); 7.8 GeV with 5 pC, up to 62 pC in 6 GeV peaks, 0.2 mrad divergence. The channel was deepened by laser inverse-bremsstrahlung heating of a capillary discharge — [Gonsalves et al. 2019, PRL 122, 084801](https://physics.aps.org/featured-article-pdf/10.1103/PhysRevLett.122.084801). The on-axis density "≈3.0×10¹⁷ cm⁻³" comes from a search snippet of the same PRL, not a full read. The 3.4×10¹⁷ figure often quoted was **not** confirmed — [search summary of PRL 122, 084801](https://physics.aps.org/featured-article-pdf/10.1103/PhysRevLett.122.084801).
- Recomputed: Eq. 6 at 850 TW, 3.0×10¹⁷ cm⁻³ and 815 nm gives 7.55 GeV against 7.8 measured (3% low). There `L_d` ≈ 15 cm is shorter than the 20 cm capillary, so the stage reaches dephasing and the full Eq. 6 gain applies — [light_in_plasma_checks.py](light_in_plasma_checks.py).

**L5. Benchmark: Texas Petawatt 2024, 10 GeV from a 10 cm nanoparticle-assisted cell** (recent and unreplicated).
- 130 ± 10 J on target; 340 pC, 10 ± 1.86 GeV bunches with 3.4 GeV rms convolved spread and 0.9 mrad rms divergence; nanoparticles ablated into helium aid injection — [Aniculaesei et al. 2024, Matter Radiat. Extremes 9, 014001, DOI 10.1063/5.0161687](https://www.researching.cn/articles/OJ5f6bc391bfb566b7). Picksley et al. cite it as "~10 GeV using 118 J" — [arXiv:2408.00740](https://arxiv.org/html/2408.00740v1). The two papers give different laser energies, 130 ± 10 J and 118 J.
- Not used as a Lu test here because the plasma density and pulse duration were not obtained (see Gaps).

**L6. Plasma-channel guiding: the matched channel depth** (established; Lu 2007 citing Esarey/Sprangle).
- Formula: `Δn_c = 1/(π r_e w₀²)`, equivalently `Δn_c/n_p ≃ 4/(k_p w₀)²`. Self-guiding needs `a₀²/8 ≳ 4/(k_p w₀)²`, i.e. `P ≳ P_c`. With a channel, `Δn_c/n_p ≃ 1/a₀` at matching — [Lu et al. 2007](https://arxiv.org/abs/physics/0612227).
- Recomputed: for w_m = 37 µm (Picksley's matched spot), `Δn_c` = 8.25×10¹⁶ cm⁻³; the two forms agree to 10⁻⁹ — [light_in_plasma_checks.py](light_in_plasma_checks.py). No published channel depth was found to test against (Gaps).
- Crate: `plasmoid::laser::matched_channel_depth(w0)`.

**L7. Ponderomotive (quiver) energy and hot-electron temperature** (established).
- Non-relativistic quiver energy `U_p = e²E₀²/(4m_eω²)`, which is 9.33×10⁻¹⁴ I[W cm⁻²] λ²[µm] eV. Recomputed: 9.337 eV at 10¹⁴ W cm⁻², 1 µm — [light_in_plasma_checks.py](light_in_plasma_checks.py).
- Relativistic ponderomotive force `f_p = −m_e c² ∇(1 + ⟨a²⟩)^{1/2}`, Eq. 4, and ponderomotive ("Wilks") energy `E_p = m_e c²(√(1 + a₀²/2) − 1)`, Eq. 6, for linear polarisation with `a₀ = 0.85 (Iλ²/10¹⁸ W cm⁻² µm²)^{1/2}` — [Macchi, Borghesi & Passoni 2013, Rev. Mod. Phys. 85, 751, arXiv:1302.1775](https://arxiv.org/abs/1302.1775). Recomputed: 0.115 MeV at a₀ = 1, and the coefficient 0.855 from `I = ε₀cE₀²/2` — [light_in_plasma_checks.py](light_in_plasma_checks.py).
- Crate: `plasmoid::laser::{quiver_energy, ponderomotive_energy}`, taking `a₀`. Link: the `Intensity` convention in `gravitics::strong_field` (Linear means `I = ε₀cE₀²/2`, matching the 0.855 coefficient).

**L8. Relativistic transparency threshold** (established; at the heart of the 2024 150 MeV record).
- `n_e < γ n_c` with `γ = √(1 + a₀²/2)` (linear polarisation), or for a thin foil `a₀ > ζ = π (n_e/n_c)(ℓ/λ)` (Macchi Eq. 5; `a₀ = ζ` is the "optimal" light-sail thickness) — [Macchi et al. 2013](https://arxiv.org/abs/1302.1775). Ziegler et al. rely on "relativistically induced transparency" at the main pulse — [Ziegler et al. 2024](https://www.nature.com/articles/s41567-024-02505-0).
- Recomputed: γ n_c = 1.22, 7.14 and 21.2 n_c at a₀ = 1, 10 and 30 — [light_in_plasma_checks.py](light_in_plasma_checks.py).
- Crate: `plasmoid::laser::relativistic_critical_density(n_c, a0)`, which is `gravitational::hf::critical_density` times γ (a cross-check in `physica`).

### Inferences
- Lu's Eq. 6 assumes the electrons run the full dephasing length. Its field falls linearly to zero over `L_d`, so a stage of length L < L_d gains `ΔE · (2x − x²)` with x = L/L_d. For Picksley that gives 7.5 GeV if the pulse is spent at `L_pd` = 20 cm and 10.1 GeV if it lasts the whole 30 cm. The measured 9.2 GeV falls between, which supports the authors' point that the pulse shape limited the gain. This truncated-stage form is my derivation from Lu's linear-field assumption and is not in the paper — [light_in_plasma_checks.py](light_in_plasma_checks.py).
- Gonsalves (stage longer than `L_d`) agrees with Eq. 6 to 3%; Picksley (stage shorter than `L_d`) needs the truncation. Together they are a clean pair of tests: one where the formula should hold and one where it must be truncated.
- `plasmoid::laser` should take `a₀` and `n_e/n_c` as inputs rather than recomputing them, so `gravitics::strong_field::a0` and `gravitational::hf::critical_density` stay the single sources. `physica` can then carry a test that wires them together.

### Gaps
- Inverse-bremsstrahlung absorption: the source for the coefficient was not obtained. The NRL Plasma Formulary 2023 page exists but returned HTTP 403 — [NRL](https://www.nrl.navy.mil/ppd/content/nrl-plasma-formulary). The textbook form `κ_ib = (ν_ei(n_c)/c)(n_e/n_c)²/√(1 − n_e/n_c)` (Kruer 1988) is left unverified in this session, and no published absorption number was found to test it against. Gonsalves 2019 uses the process (laser heating of the capillary), but its numbers were not obtained.
- The Aniculaesei 2024 density and pulse duration were not obtained, so that 10 GeV point was not run through Lu.
- No published channel depth Δn_c was found for the Picksley HOFI channel. The HOFI channel is not purely parabolic either.
- Lu's printed `L_d` and `L_pd` (1.31 and 0.96 cm) are 7–9% below recomputation from the stated inputs. The cause was not found; set test tolerances at 10%.

## Laser-driven ion acceleration: TNSA and RPA scalings, hole boring, and the 2022–2025 records

### Takeaway
The hole-boring and light-sail formulas (Macchi 2013 Eqs. 28–36) are exact closed forms, and
their internal identities were all verified. Mora's TNSA front velocity is also closed form, but
its cut-off grows without bound with time and needs an empirical stopping time. The 2024–2025
records above 100 MeV (150 MeV at DRACO, >110 MeV at CoReLS) come from hybrid
transparency/RPA/Coulomb regimes. No single scaling reproduces them, so they belong in a records
registry, not in a formula test.

### Cited Findings
**I1. Hole-boring (piston) velocity and ion energy** (established theory; experiment "less clear" in solids).
- Momentum balance in the piston frame (Eq. 28): `(2I/c)(1 − β)/(1 + β) = n_i γ_hb (2m_i γ_hb v_hb) v_hb`. Solution (Eq. 29): `v_hb/c = √Π/(1 + √Π)`, with `Π = I/(m_i n_i c³) = (Z/A)(n_c/n_e)(m_e/m_p) a₀²`. Maximum energy per nucleon of the reflected ions (Eq. 30): `E_max = 2 m_p c² Π/(1 + 2√Π)`, which reduces to `2m_p c² Π` for Π ≪ 1. "Experimental evidence of HB acceleration in solid targets is less clear at present" — [Macchi et al. 2013, RMP 85, 751](https://arxiv.org/abs/1302.1775).
- Recomputed identities (to 10⁻⁹): Eq. 30 equals the kinetic energy of an ion bounced at relativistic velocity `2β/(1 + β²)`, and Eq. 29 satisfies Eq. 28, for Π = 10⁻⁴, 0.1 and 1 — [light_in_plasma_checks.py](light_in_plasma_checks.py).
- Published test (Palmer et al. 2011, PRL 106, 014801, via Macchi): 10 µm CO₂ laser, ~6×10¹⁵ W cm⁻² circularly polarised, "a₀ ≃ 0.5", hydrogen jet at a few n_c, protons up to 1.2 MeV with narrow spread. Energies "were actually higher than expected taking the vacuum laser intensity", which Macchi attributes to self-focusing — [Macchi et al. 2013](https://arxiv.org/abs/1302.1775).
- Recomputed: a₀ = 0.47 from `I = a₀² m_e c³ n_c` (circular; the same as `Intensity::Peak`). Eq. 30 then gives 0.22, 0.11 and 0.055 MeV at n_e = 1, 2 and 4 n_c, a factor of 5–20 below the 1.2 MeV observed. This reproduces the known shortfall — [light_in_plasma_checks.py](light_in_plasma_checks.py).
- Recent: a time-resolved Doppler-spectroscopy measurement of hole boring at PHELIX confirmed an analytical (non-relativistic) hole-boring model with arbitrary density and intensity profiles — [Nature Communications, DOI 10.1038/s41467-021-27363-9 (PMC8636483)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8636483/). Its numbers were not extracted.
- Crate: `plasmoid::laser::{hole_boring_pi, hole_boring_velocity, hole_boring_ion_energy}`. Link: `gravitics::strong_field::a0`. It is related to `gravitics::propulsion::photon_thrust`: the β → 0 limit of the moving-mirror pressure `(2I/c)(1 − β)/(1 + β)` is the 2I/c of a perfect reflector (not checked against that function's exact form). Status: established (theory); contested or partial (quantitative agreement with experiment).

**I2. Light-sail RPA** (established theory; recent experiments mixed with other mechanisms).
- Eqs. 31–36: `F(t) = (2/σc²) ∫ I dt′` with `σ = m_i n_i ℓ`; `β = [(1 + F)² − 1]/[(1 + F)² + 1]`; `E = m_p c² F²/[2(F + 1)]`; `η = 2β/(1 + β) = 1 − 1/(F + 1)²`; `γ(t) = sinh(u) + 1/(4 sinh u)` with `u = (1/3) asinh(3Ωt + 2)`, `Ω = (Z m_e a₀²/(A m_p ζ)) ω`. Practical units: `F_∞ = 2.2 F_{1e8} ρ₁⁻¹ ℓ₁₀⁻¹`, where `F_{1e8}` is fluence in 10⁸ J cm⁻², ρ₁ = density in g cm⁻³, ℓ₁₀ = thickness in units of 10 nm — [Macchi et al. 2013](https://arxiv.org/abs/1302.1775).
- Recomputed: `F` = 2.225 for 10⁸ J cm⁻² on a 10 nm, 1 g cm⁻³ foil (matches 2.2). Eq. 35 equals `m_p c²(γ − 1)` with Eq. 34's β (720 MeV per nucleon at F = 2.225). Eq. 36's two forms agree — [light_in_plasma_checks.py](light_in_plasma_checks.py).
- Crate: `plasmoid::laser::{sail_fluence, sail_velocity, sail_energy, sail_efficiency}`. Link: `gravitics::propulsion::photon_thrust`; it is the same relativistic mirror physics as a laser-driven light sail.

**I3. TNSA, Mora's isothermal expansion** (established; foundational, Mora 2003 PRL 90, 185002).
- Ion-front velocity `u_f(t) ≃ 2c_s ln(τ + √(τ² + 1))` and field `E(t) ≃ √(2/e) E₀/√(τ² + 1)`, with `τ = ω_pi t/√(2e)` (e is Euler's number), Eqs. 23–24. "The major drawback ... the cut-off energy diverges logarithmically with time" — [Macchi et al. 2013](https://arxiv.org/abs/1302.1775). With `c_s² = Z T_h/m_i`, `E_max = 2 Z T_h [asinh τ]²`, and `T_h` from L7.
- Recomputed: `(1/2) m_i u_f² = 2 Z T_h asinh²(τ)` exactly — [light_in_plasma_checks.py](light_in_plasma_checks.py).
- Crate: `plasmoid::laser::tnsa_max_energy(t_h, z, tau)`. A published TNSA cut-off with all inputs, needed to test it, was not obtained (Gaps).

**I4. Proton-energy records, 2024–2025** (recent and unreplicated individually; >100 MeV now reached at two facilities by different routes).
- DRACO-PW (HZDR): "a spectrally separated high-energy component of up to 150 MeV" from solid-density plastic foils; a controlled pre-pulse makes the foil relativistically transparent at the main pulse ("cascaded acceleration regimes"; 3D PIC) — [Ziegler et al. 2024, Nature Physics 20, 1211, DOI 10.1038/s41567-024-02505-0](https://www.nature.com/articles/s41567-024-02505-0). The facility previously reached ~80 MeV, and the authors expect ">200 MeV" to be possible (a prediction) — [chemeurope/HZDR release](https://chemeurope.com/en/news/1183445/milestone-in-plasma-acceleration-new-energy-record-for-next-generation-proton-accelerators).
- CoReLS (multi-PW, Korea): cut-off ">110 MeV" from a conduction-restricted nanometre polymer foil of finite lateral size, attributed to radiation pressure plus Coulomb repulsion, with fields above 10¹³ V m⁻¹ — [Shou et al. 2025, Nat. Commun. 16, DOI 10.1038/s41467-025-56667-3](https://www.nature.com/articles/s41467-025-56667-3).
- Laser energy, intensity and thickness for these shots were not obtained, so neither record was run through I1–I3.
- Crate: a `plasmoid::laser::PROTON_RECORDS` registry (facility, year, MeV, mechanism, DOI, status), in the style of `gravitational::hf::REGISTRY`.

### Inferences
- The two >100 MeV results use transparency (L8) and Coulomb-assisted RPA. Testing them would need the target's areal density and the laser fluence, which then fix `ζ` and `F` in I2. The light-sail energy at the measured fluence is the right comparison, not TNSA.
- For a small, safe, testable module, I1 and I2 are closed forms with exact identities. I3 is the classic baseline, and its stopping time must be an explicit argument (not hidden), because it diverges.

### Gaps
- No published TNSA shot with T_h, n_e, Z and acceleration time was obtained, so Fuchs et al. 2006's empirical `t_acc ≈ 1.3 τ_L` was not verified in this session.
- Laser parameters (energy on target, duration, spot) for Ziegler 2024 and Shou 2025 were not extracted. The full texts were not read.
- No 2022–2026 relativistic hole-boring velocity measurement with published inputs was found. The PHELIX Doppler result is non-relativistic, and its numbers were not extracted.

## Faraday rotation and rotation measure: pulsars, fast radio bursts and the DM–RM field estimate

### Takeaway
`RM = (e³/8π²ε₀m_e²c³) ∫ n_e B∥ dl` = 2.631×10⁻¹³ in SI = 0.8119 rad m⁻² per (cm⁻³ µG pc), and
`⟨B∥⟩ = 1.232 RM/DM µG`. All three were recomputed from CODATA. The same coefficient and
polarisation algebra serve pulsars, FRBs and tokamak polarimeters. The published tests are a
pulsar catalogue (ATNF), FRB 20121102A's source-frame RM of 1.46×10⁵ rad m⁻² falling to
30,755 rad m⁻² by 2023 (a 70% drop), FRB 20190520B's two sign reversals, and the CHIME
non-repeater median |RM| of 53.8 rad m⁻².

### Cited Findings
**R1. Faraday rotation angle and the RM coefficient** (established).
- `Δψ = RM λ²`, with `RM = K ∫ n_e B∥ dl` and `K = e³/(8π² ε₀ m_e² c³)`. Recomputed: K = 2.631×10⁻¹³ rad T⁻¹ (SI: n_e in m⁻³, B in T, l in m), which converts to 0.8119 rad m⁻² for n_e in cm⁻³, B in µG and l in pc. The 2.63×10⁻¹³ coefficient is the standard tokamak-polarimetry form (with λ² in m²) — [light_in_plasma_checks.py](light_in_plasma_checks.py).
- In tokamaks the Faraday angle scales as λ² ∫ n_e B_pol dl and the Cotton–Mouton phase as λ³ ∫ n_e B_tor² dl. At 118.8 µm both are large and must be separated; a 57.2 µm design gives a Faraday rotation of about 25–30° — [NIFS annual report 2011](https://www.nifs.ac.jp/report/annrep11/pdf/144.pdf). Z-cut quartz windows add negligible rotation at 118.8 µm and 0.3 T — [Plasma Fusion Res. 13, 1405112 (2018)](https://www.jspf.or.jp/PFR/PFR_articles/pfr2018/pfr2018_13-1405112.html).
- Crate: `magnetics::optics::{rotation_measure, faraday_angle}`. Links: `plasmoid::tokamak` (polarimetry of a machine, in a `physica` example), `magnetics::main_field` (ionospheric Faraday rotation of satellite signals). Status: established.

**R2. Plasma Cotton–Mouton phase** (established).
- `φ_CM = (e⁴/(16π³ ε₀ m_e³ c⁴)) λ³ ∫ n_e B⊥² dl`; recomputed coefficient 2.457×10⁻¹¹ rad T⁻² m⁻² (SI) — [light_in_plasma_checks.py](light_in_plasma_checks.py). The λ³ B_tor² scaling is confirmed by the ITER polarimetry literature — [NIFS 2011](https://www.nifs.ac.jp/report/annrep11/pdf/144.pdf).
- Recomputed demonstration (not a published test): using `plasmoid::tokamak::ITER`'s 5.3 T with n_e = 10²⁰ m⁻³, a 4 m chord and λ = 118.8 µm gives φ_CM ≈ 27°. This shows why ITER moved toward shorter wavelengths — [light_in_plasma_checks.py](light_in_plasma_checks.py).
- Crate: `magnetics::optics::cotton_mouton_phase`. Link: `plasmoid::tokamak::ITER`.

**R3. DM–RM mean parallel field** (established).
- `⟨B∥⟩ = RM/(0.812 DM) = 1.232 RM/DM µG` (RM in rad m⁻², DM in pc cm⁻³); recomputed coefficient 1.2317 — [light_in_plasma_checks.py](light_in_plasma_checks.py).
- Published inputs, from the ATNF catalogue v2.6.1 — [ATNF psrcat](https://www.atnf.csiro.au/research/pulsar/psrcat/):

  | Pulsar | DM (pc cm⁻³) | RM (rad m⁻²) | ⟨B∥⟩ recomputed |
  |---|---|---|---|
  | J0332+5434 (B0329+54) | 26.7641 | −64.33 ± 0.06 | −2.960 µG |
  | J0534+2200 (Crab) | 56.77 | −45.44 ± 0.08 | −0.986 µG |
  | J0835−4510 (Vela) | 67.771 | +31.38 ± 0.01 | +0.570 µG |
  | J1939+2134 (B1937+21) | 71.0 | +11.0 ± 0.3 | +0.191 µG |

- FAST GPPS: Faraday rotation for 134 faint halo pulsars plus 311 more RMs; halo field ≈ 2 µG with scale height ≥ 2.7 ± 0.3 kpc; disk-field reversals in the Perseus arm and beyond — [Xu, Han, Wang & Yan 2022, Sci. China PMA 65, 129704, arXiv:2211.11302](https://arxiv.org/abs/2211.11302).
- Crate: `magnetics::optics::mean_parallel_field(rm, dm)`. Link: `gravitational::propagation`'s dispersion (see P1), which supplies DM from a measured delay.

**R4. FRB 20121102A: the extreme, decaying RM** (recent and replicated across Arecibo/GBT and FAST).
- Source-frame RM +1.46×10⁵ and +1.33×10⁵ rad m⁻², 7 months apart; z = 0.193 — [Michilli et al. 2018, Nature 553, 182, arXiv:1801.03965](https://arxiv.org/abs/1801.03965). Recomputed observed-frame values (`RM_obs = RM_src/(1 + z)²`): 1.026×10⁵ and 9.35×10⁴ rad m⁻² — [light_in_plasma_checks.py](light_in_plasma_checks.py).
- FAST, 3 March 2023: RM = +30,755 (−15/+10) rad m⁻² by RM synthesis, continuing the 2016–2019 decline — [ATel #15980](https://www.astronomerstelegram.org/?read=15980). This is "a 70% decrease compared to that from December 2016", with a decay slope of −28.5 (+2.0/−1.8) rad m⁻² day⁻¹ and σ_RM down 13% — [arXiv:2507.15790](https://arxiv.org/abs/2507.15790). Recomputed: 1 − 30,755/1.026×10⁵ = 0.700 — [light_in_plasma_checks.py](light_in_plasma_checks.py).
- The ATel's "1.02×10⁵ rad m⁻² per year" decline rate is inconsistent with its own values and is probably a typo — [ATel #15980](https://www.astronomerstelegram.org/?read=15980).

**R5. FRB 20190520B: two RM sign reversals in 17 months** (recent; two telescopes in one study).
- "The amount of Faraday rotation is highly variable", it "changes sign twice" over seventeen months, and this is interpreted as a turbulent magnetised screen 10⁻⁵–100 pc from the source — [Anna-Thomas et al. 2023, Science 380, 599, arXiv:2202.11112](https://arxiv.org/abs/2202.11112). The press figure of roughly ±10⁴ rad m⁻² was not confirmed from the paper — [Newswise](https://www.newswise.com/articles/view/792285/?sc=c128).

**R6. CHIME/FRB non-repeater RM distribution** (recent; catalogue).
- 128 non-repeaters (400–800 MHz): 89 linearly polarised above 6σ, 29 unpolarised, 10 contaminated; |RM| after Milky Way subtraction 0.5–1160 rad m⁻², median 53.8 rad m⁻²; median linear polarisation 63%; "most ... consistent with Milky Way-like host galaxies" — [Pandhi et al. 2024, ApJ 968, 50, DOI 10.3847/1538-4357/ad40aa, arXiv:2401.17378](https://arxiv.org/abs/2401.17378).
- Crate: a `magnetics::optics::FRB_RM` registry (source, RM, DM, z, epoch, reference, status), next to R3.

### Inferences
- One module, `magnetics::optics`, covers tokamak polarimetry (λ² and λ³ coefficients), pulsar fields and FRB environments with the same two constants. A `physica` example can run ITER's polarimeter numbers from `plasmoid::tokamak::ITER` and a pulsar's ⟨B∥⟩ from the same function.
- A lower bound on the host field follows from `B∥ ≥ 1.232 RM_src/DM_host`. It was not computed here because FRB 20121102A's DM_host was not obtained.

### Gaps
- The full Anna-Thomas et al. RM values and DM were not read (the abstract is qualitative).
- No Han et al. 2018 ApJS or 2023–2025 FAST paper with a fitted Galactic field model and published pulsar inputs was read in full; the halo value comes from the 2022 abstract.
- No published ITER Faraday angle at 118.8 µm with its line-integral inputs was found. The 25–30° figure at 57.2 µm lacks stated inputs, so it is an order-of-magnitude check only.

## Magneto-optics in matter: Verdet constants, Cotton–Mouton and Voigt effects, MOKE, inverse Faraday effect

### Takeaway
The Faraday rotation `θ = V B L` and the Cotton–Mouton birefringence `Δn = k_CM B²` have clean,
published test values. Helium at 1064 nm has V = (3.87 ± 0.12)×10⁻⁵ rad T⁻¹ m⁻¹ atm⁻¹ and
k_CM = (2.19 ± 0.12)×10⁻¹⁶ T⁻² atm⁻¹, agreeing with ab-initio theory. TGG at 1064 nm is about
−38 to −40 rad T⁻¹ m⁻¹ (secondary sources). Matched against `gravitics::strong_field::birefringence`,
helium's Cotton–Mouton equals the QED vacuum's at 1.8×10⁻⁸ atm (1.8 mPa), a sharp cross-crate
test of the residual-gas requirement in PVLAS-type experiments. The inverse Faraday effect in
plasmas has kilotesla predictions but no 2022–2026 kilotesla measurement was found.

### Cited Findings
**M1. Verdet rotation** (established).
- `θ_F = V B∥ L`; for a gas, `V = π k_F/λ` with `Δn_F = k_F B∥`, and in a cavity of finesse F, `Θ_F = (2F/π) V B∥ L_B` — [Cadène et al. 2013, PRA 88, 043815, arXiv:1307.6394](https://ar5iv.labs.arxiv.org/html/1307.6394).
- Helium, 1064 nm: V = (3.87 ± 0.12)×10⁻⁵ atm⁻¹ rad T⁻¹ m⁻¹ at 294 ± 1 K, or (4.17 ± 0.13)×10⁻⁵ normalised to 273.15 K. Theory: 4.06×10⁻⁵ (Ekström 2005) and (4.09 ± 0.02)×10⁻⁵ (Savukov 2012) — [Cadène et al. 2013](https://ar5iv.labs.arxiv.org/html/1307.6394).
- TGG, 1064 nm: −38 rad T⁻¹ m⁻¹ (−131 at 632 nm) — [Wikipedia, TGG](https://en.wikipedia.org/wiki/Terbium_gallium_garnet); 38.5 rad T⁻¹ m⁻¹ for TGG and 36 for KTF — [Northrop Grumman SYNOPTICS datasheet](https://cdn.northropgrumman.com/-/media/AprimoConnector/SPACE/Products/TGG-and-KTF-Faraday-Crystals-Specs-Northrop-Grumman-SYNOPTICS.pdf). The primary measurement compilation (Vojna, Slezák, Lucianetti & Mocek 2019, Appl. Sci. 9, 3160, DOI 10.3390/app9153160) returned HTTP 403 and was not read.
- Recomputed: a 45° isolator in TGG at −40 rad T⁻¹ m⁻¹ needs B L = 0.0196 T m, e.g. 1 T over 19.6 mm — [light_in_plasma_checks.py](light_in_plasma_checks.py).
- Crate: `magnetics::optics::faraday_rotation(verdet, b, length)` with a small table of V(λ, T). Status: established.

**M2. Cotton–Mouton (Voigt) birefringence in gases, and the QED vacuum** (established; links `gravitics`).
- `Δn_CM = k_CM B⊥²`; helium at 1064 nm: k_CM = (2.19 ± 0.12)×10⁻¹⁶ T⁻² atm⁻¹ at 293 K, or (2.35 ± 0.13)×10⁻¹⁶ at 273.15 K, against the theoretical 2.4036×10⁻¹⁶ (Coriani et al. 1999), within 1σ — [Cadène et al. 2013](https://ar5iv.labs.arxiv.org/html/1307.6394).
- Already covered: `gravitics::strong_field::birefringence` (`Δn = 3A_e B²`, the vacuum's Cotton–Mouton). Recomputed 3A_e = 3.97×10⁻²⁴ T⁻². Helium's Cotton–Mouton equals this at p = 1.8×10⁻⁸ atm = 1.8 mPa — [light_in_plasma_checks.py](light_in_plasma_checks.py).
- Crate: `magnetics::optics::cotton_mouton(k_cm, b, pressure)`. Link: `gravitics::strong_field::birefringence` (a `physica` test computing the residual-gas pressure at which gas mimics QED). Status: established.

**M3. Inverse Faraday effect (IFE) in plasmas** (established phenomenon; kilotesla regime recent and unreplicated, simulation only).
- IFE: "a quasi-static axial magnetic field that is self-generated when a circularly polarized beam propagates in a plasma"; 3D PIC with unconventional polarisation states — [Liu, Jia & Zheng 2023, Matter Radiat. Extremes 8, 014405](https://www.researching.cn/articles/OJ447086f19a3348ae). Spin- and orbital-angular-momentum beams are predicted to make kilotesla, multi-picosecond axial fields over hundreds of microns — [Longman & Fedosejevs, arXiv:2111.13243](https://arxiv.org/abs/2111.13243). Twisted multi-beam schemes are also predicted to give multi-kT fields — [arXiv:2312.15298](https://arxiv.org/html/2312.15298v1).
- Measurements: fields "in excess of 7 MG" (700 T) in underdense helium at up to 10¹⁹ W cm⁻² — [Najmudin et al., HAL hal-00520546](https://hal-polytechnique.archives-ouvertes.fr/hal-00520546). Experiments with circular polarisation remain scarce — [arXiv:2405.14375](https://arxiv.org/html/2405.14375v1).
- Crate: none proposed until a closed-form B_IFE with a published test value is obtained (Gaps).

### Inferences
- M2 joins `magnetics` (gas Cotton–Mouton) to `gravitics` (QED vacuum birefringence) with one comparable quantity, Δn per T². It is the most direct "light in a magnetic field" cross-crate test available.
- The Voigt effect in atomic vapours and the Cotton–Mouton effect in liquids and gases are the same quadratic-in-B transverse birefringence. One function covers both, with the material constant as an input.

### Gaps
- MOKE: no source with a testable published Kerr rotation and its inputs was retrieved. Nothing is proposed.
- TGG temperature dependence: paramagnetic V ∝ 1/(T − T₀) is the expected form, but no fitted TGG constants at 1064 nm were obtained. The primary Vojna et al. 2019 paper could not be read (HTTP 403).
- No closed-form IFE field with a published benchmark was found. The 7 MG measurement is pre-2022, and its formula was not extracted.

## Photon propagation in a plasma: effective photon mass, group delay, whistlers, synchrotron and cyclotron emission

### Takeaway
The cold-plasma group delay `Δt = k_DM DM/f²`, with `k_DM = e²/(8π²ε₀m_ec)` = 4.1488 ms GHz²
per (pc cm⁻³), is exactly `gravitational::propagation::massive_graviton_delay` evaluated at the
plasma's effective photon mass `m c² = ħω_p`. For a uniform medium the two agree to 10⁻¹². This
identity is a free cross-module test. It also explains why FRB photon-mass bounds (4.7–4.9×10⁻⁵¹ kg,
2025) must model DM: interstellar plasma alone gives an effective mass of 1.15×10⁻⁴⁷ kg.

### Cited Findings
**P1. Dispersion delay and the photon's plasma mass** (established; links two modules).
- `Δt = k_DM DM (f₁⁻² − f₂⁻²)`, with `k_DM = e²/(8π² ε₀ m_e c)`. Recomputed: 4.148808 ms GHz² cm³ pc⁻¹. The effective photon mass is `m_eff = ħω_p/c²` — [light_in_plasma_checks.py](light_in_plasma_checks.py).
- Identity, recomputed: for n_e = 0.03 cm⁻³ over 1 kpc at 1.4 GHz, the plasma delay of 63.50 ms equals `(d/c)·½(m_eff c²/hf)²`, which is `gravitational::propagation::massive_graviton_delay` with m = h f_p, to 10⁻¹² — [light_in_plasma_checks.py](light_in_plasma_checks.py). `f_p(0.03 cm⁻³)` = 1.555 kHz, matching the doc comment of `gravitational::hf::plasma_frequency` (1.56 kHz) — [crates/gravitational/src/hf.rs](../../../../crates/gravitational/src/hf.rs).
- ISM effective mass, recomputed: ħω_p = 6.43×10⁻¹² eV = 1.15×10⁻⁴⁷ kg at 0.03 cm⁻³ — [light_in_plasma_checks.py](light_in_plasma_checks.py).
- Already covered: `gravitational::hf::{plasma_frequency, plasma_index}`. New: `gravitational::propagation::plasma_delay(dm, f)` and `photon_plasma_mass(n_e)`. Status: established.

**P2. Photon rest-mass bounds from FRBs** (recent and unreplicated; method-dependent).
- 1σ upper limits of 4.83, 4.71 and 4.86×10⁻⁵¹ kg (ΛCDM, wCDM, w₀wₐCDM) from localised FRBs with an improved extragalactic DM distribution; "accepted for publication in PRD" (v2, 23 Jan 2026) — [Zhang, Liu, Yu & Wu, arXiv:2511.14163](https://arxiv.org/abs/2511.14163). Pulsar-timing plus dedispersed-FRB bound: 9.52×10⁻⁴⁶ kg (5.34×10⁻¹⁰ eV/c²) — [arXiv:2403.06422](https://arxiv.org/html/2403.06422v3).
- Crate: add these to a `gravitational::propagation` bounds registry next to the graviton-mass bounds. The photon-mass limit plays the same role as `massive_graviton_delay`'s graviton bounds.

**P3. Magnetised-plasma cut-offs and cyclotron frequency** (established; links `magnetics` field to `gravitational::hf`).
- `f_ce = eB/(2π m_e)` = 27.99 GHz T⁻¹; R cut-off `f_R = ½[f_ce + √(f_ce² + 4f_pe²)]`; upper hybrid `f_UH = √(f_ce² + f_pe²)`. Recomputed for ITER (5.3 T from `plasmoid::tokamak::ITER`, n_e = 10²⁰ m⁻³): f_ce 148.4 GHz, f_pe 89.8 GHz, f_R 190.6 GHz, f_UH 173.4 GHz — [light_in_plasma_checks.py](light_in_plasma_checks.py). These formulas were not checked against a fetched source this session (textbook Stix forms). Already covered: `plasmoid::pinch::cyclotron_energy` (ħω_ce in eV), `plasmoid::aneutronic::synchrotron` (NRL power density `6.21×10⁻¹⁷ B² n_e T_e`).
- Crate: `magnetics::optics::{cyclotron_frequency, r_cutoff, upper_hybrid}`, taking `f_pe` from `gravitational::hf::plasma_frequency` as an argument.

**P4. Magnetised laser-driven implosions: light, plasma and field together** (recent and replicated at NIF in 2022 and 2024).
- A 26 T axial field on a room-temperature, indirectly driven D₂ capsule raised the ion temperature by ~40% and the neutron yield by ~3.2×, with an estimated compressed core field of ~4.9 kT — [Moody et al. 2022, PRL 129, 195002](https://link.aps.org/accepted/10.1103/PhysRevLett.129.195002). Follow-ups: 28 T gave +28% T_i and 2.19× yield; 12 T gave +33% and 2.5× — [arXiv:2503.00186](https://arxiv.org/html/2503.00186v1).
- Recomputed: 26 T → 4.9 kT implies a radial convergence of 13.7 under `B r² = const` — [light_in_plasma_checks.py](light_in_plasma_checks.py). Already covered: `magnetics::magnets::flux_compression(b0, r0, r, lambda)`.
- Crate: a `physica` example joining `magnetics::magnets::flux_compression` with `plasmoid` yields. The convergence ratio is my inference, not a published number.

### Inferences
- Because plasma dispersion and a Proca mass both give an f⁻² delay, the photon-mass bound depends on separating the two through their different redshift dependence. The 4.8×10⁻⁵¹ kg bound is below the effective plasma mass of the mean intergalactic medium. This is consistent only because the bound applies to the residual after the DM is modelled. A test in `gravitational` can show the degeneracy explicitly by computing both delays at one frequency.
- `gravitational::hf::plasma_index` already gives the phase index. The group index `n_g = 1/n` (cold, unmagnetised) gives the same delay to first order and would make P1 self-checking inside one crate.

### Gaps
- Whistler dispersion (Eckersley `t = D f^{−1/2}`, `D = (1/2c)∫ f_pe/√f_ce ds`): no source or published dispersion value was fetched this session, so it is not proposed with a test.
- No 2022–2026 synchrotron or cyclotron benchmark beyond what `plasmoid::aneutronic` already tests was gathered.

## What links these to existing crates (plasmoid, magnetics, gravitational::hf, gravitics::strong_field)

### Takeaway
Three new modules would cover the catalogue with the least duplication: `plasmoid::laser`
(L1–L8, I1–I4), `magnetics::optics` (R1–R6, M1–M2, P3) and two functions in
`gravitational::propagation` (P1–P2). The cross-crate tests belong in `physica`, because only
`physica` depends on every crate.

### Cited Findings
- Crate dependency graph: `plasmoid` and `gravitational` depend on `metrology` only; `magnetics` on `geopotential` and `metrology`; `gravitics` on `geopotential`, `gravitational` and `metrology`; `physica` on all seven — [crates/*/Cargo.toml](../../../../crates/).
- Existing functions to reuse, not duplicate: `gravitics::strong_field::a0` (with `Intensity::Linear`/`Peak`), `birefringence`; `gravitational::hf::{plasma_frequency, critical_density, plasma_index}`; `gravitational::propagation::massive_graviton_delay`; `magnetics::magnets::flux_compression`; `plasmoid::tokamak::ITER`; `plasmoid::aneutronic::synchrotron`; `plasmoid::pinch::cyclotron_energy` — [crate READMEs and sources](../../../../crates/).
- Cross-crate test list, each with its number from the script — [light_in_plasma_checks.py](light_in_plasma_checks.py):
  1. `plasmoid::laser::lu_energy_gain` with a₀ from `gravitics::strong_field::a0`: Gonsalves 7.55 GeV against 7.8; Picksley matched spot 52 µm against 53 µm.
  2. `gravitational::propagation::massive_graviton_delay(h f_p, f, d)` against `plasma_delay(DM, f)`: equal to 10⁻¹².
  3. `magnetics::optics::cotton_mouton(He)` against `gravitics::strong_field::birefringence`: equal at 1.8 mPa.
  4. `magnetics::optics::cotton_mouton_phase` with `plasmoid::tokamak::ITER`: ≈ 27° at 118.8 µm over a 4 m chord.
  5. `plasmoid::laser::relativistic_critical_density` against `gravitational::hf::critical_density` × √(1 + a₀²/2).
  6. `magnetics::optics::mean_parallel_field` on ATNF pulsars: −2.960 µG for B0329+54.

### Inferences
- The highest-value additions, in order: (1) `magnetics::optics` RM/DM/Verdet/Cotton–Mouton, because they are exact constants with catalogue tests and join three crates; (2) `plasmoid::laser` Lu scalings plus hole boring and light sail, with closed forms, BELLA benchmarks and exact identities; (3) `gravitational::propagation::plasma_delay`, a one-liner with an exact identity test against existing code.
- Proton and FRB records fit the registry pattern already used by `gravitational::hf::REGISTRY`, with a status field.

### Gaps
- Whether `metrology` should host `n_c`, `ω_p` and `r_e` so that `plasmoid` and `magnetics` stop duplicating them is a design choice for the user; it was not assessed.
