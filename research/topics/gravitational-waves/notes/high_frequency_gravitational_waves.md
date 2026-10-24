# High-frequency gravitational waves (HFGW, ~10 kHz to THz and above): sources, detection concepts, limits, and formulas (2023–2026 literature)

Scope note. The backbone source is the 2025 update of the Living Review: N. Aggarwal, O. D. Aguiar, D. Blas, ..., V. Domcke, S. Ellis, G. Franciolini, J. Kopp, S. M. Lee, ..., G. White, "Challenges and Opportunities of Gravitational Wave Searches above 10 kHz", arXiv:2501.11723 (v2 dated 3 Feb 2026), published in Living Reviews in Relativity 28, article 10 (3 Nov 2025), open access, CC BY. Below it is cited as "Review 2025" with section and equation numbers, taken from the arXiv v2 text extracted in this session. It succeeds the 2021 version ("Challenges and opportunities of gravitational-wave searches at MHz to GHz frequencies", Living Rev. Relativ. 24, 4 (2021), arXiv:2011.12414). I did not open the 2021 version in this session, so I cite it only where the 2025 text refers to it. Equations in Review 2025 use natural units (c = ħ = G = 1, Heaviside-Lorentz electromagnetism) unless noted. SI conversions marked "derived" are my own, checked numerically in this session.

## Q1. Which sources produce HFGWs, what are their predicted h(f), h_c(f) or Ω_GW(f), and what is the BBN/ΔN_eff bound?

### Takeaway
Every cosmological (stochastic) HFGW source is capped by the dark-radiation bound Ω_GW,0 h² < 1.1×10⁻⁶ (ΔN_eff < 0.18). At high f that is a tiny strain: h_c ≲ 1.3×10⁻²¹ (1 Hz/f), i.e. about 10⁻³⁰ at 1 GHz. The only sources that can locally exceed this are late-Universe coherent or transient ones: mergers and close encounters of light primordial black holes (PBHs), boson clouds from black-hole superradiance, and GW memory. These are not bounded by ΔN_eff, but their strain falls off as 1/D and the expected nearby event rates are small. The 2025 review lists every source in its Table 1 (late Universe) and Table 2 (stochastic).

### Cited Findings

**BBN/CMB bound (Review 2025 §4.2, Eqs. 82–87)**
- Extra radiation is counted as Δρ_rad = (π²/30)(7/4)(4/11)^{4/3} ΔN_eff T⁴. The GW bound is ρ_GW/ρ_γ|_{T≲MeV} ≤ (7/8)(4/11)^{4/3} ΔN_eff^max. Redshifted to today: (ρ_GW h²/ρ_c)_0 ≤ Ω_rad,0 h² (7/8)(4/11)^{4/3} ΔN_eff^max = 5.6×10⁻⁶ ΔN_eff^max, with Ω_rad,0 h² = (ρ_γ/ρ_c)_0 h² = 2.47×10⁻⁵ (photons only). The SM value is N_eff = 3.043. — [Review 2025 §4.2](https://arxiv.org/abs/2501.11723); N_eff from [Cielo et al. 2023, arXiv:2306.05460](https://arxiv.org/abs/2306.05460)
- BBN alone gives ΔN_eff < 0.2 (95% CL). The combined CMB+BBN analysis gives ΔN_eff < 0.18 (95% CL), hence **Ω_GW,0 h² < 1.1×10⁻⁶** (Eq. 86). This applies to backgrounds already present before BBN (f ≥ 1.5×10⁻¹² Hz). With CMB-only limits and homogeneous initial conditions it tightens to **Ω_GW,0 h² < 2.9×10⁻⁷** (Eq. 87), valid down to f ≳ 10⁻¹⁵ Hz. If future CMB experiments reach the theory uncertainty on N_eff (~10⁻³), the bound would become Ω h² < 5.6×10⁻⁹. — [Review 2025 §4.2](https://arxiv.org/abs/2501.11723), combined analysis [Yeh, Shelton, Olive, Fields 2022, arXiv:2207.13133](https://arxiv.org/abs/2207.13133)
- The bound limits the *integrated* GW energy density inside the horizon. It acts as a bound on Ω_GW(f) itself except for spectra with very narrow peaks (Δf ≪ f). — [Review 2025 §4.2](https://arxiv.org/abs/2501.11723)
- "No proposal above the LIGO–VIRGO–KAGRA band currently reaches below the cosmological bound of Ω_GW h² ≲ 10⁻⁶". So cosmological sources "seem currently out of reach", and current sensitivities mostly correspond to Ω_GW ≫ 1, which should be read as a probe of *local* GW overdensities. — [Review 2025 §3](https://arxiv.org/abs/2501.11723)

**Late-Universe coherent and transient sources (Review 2025 §4.1, Table 1)**
- *PBH inspirals.* f_ISCO = 4400 Hz (M_⊙/M) with M = m₁+m₂ (Eq. 46). The band 10⁴–10¹⁵ Hz maps to M ≈ 10⁻¹²–10⁻¹ M_⊙. 0-PN Fourier-domain strain (Eq. 47): h̃(f) = (5/24)^{1/2} π^{-2/3} D⁻¹ (GM_c)^{5/6} f^{-7/6} e^{iψ} Q(θ,φ,ϕ), where M_c = (m₁m₂)^{3/5}/(m₁+m₂)^{1/5}. For equal masses this is ≈ 2×10⁻³⁷ s (kpc/D)(m_PBH/10⁻¹² M_⊙)^{5/6}(f/GHz)^{-7/6}, equivalently ≈ 10⁻³¹ s (kpc/D)(m_PBH/10⁻⁵ M_⊙)^{5/6}(f/GHz)^{-7/6} (Table 1). — [Review 2025 §4.1.2](https://arxiv.org/abs/2501.11723)
- PBH chirp timing. ḟ = (96/5) π^{8/3} (GM_c/c³)^{5/3} f^{11/3} ≃ 4.6×10¹¹ Hz² (m_PBH/10⁻⁹ M_⊙)^{5/3}(f/GHz)^{11/3} (Eq. 49). N_cycles = f²/ḟ ≃ 2.2×10⁶ (f/GHz)^{-5/3}(m_PBH/10⁻⁹ M_⊙)^{-5/3} (Eq. 48). Time to coalescence τ(f) ≈ 83 s (m_PBH/10⁻¹² M_⊙)^{-5/3}(f/GHz)^{-8/3} (Eq. 50). — [Review 2025 §4.1.2](https://arxiv.org/abs/2501.11723)
- PBH merger rate for early (primordial) binaries (Eq. 51): dR/(d ln m₁ d ln m₂) = 0.0038 kpc⁻³ yr⁻¹ × f_PBH^{53/37} (t/t₀)^{-34/37} (M/10⁻¹² M_⊙)^{-32/37} [m₁m₂/(m₁+m₂)²]^{-34/37} S(M, f_PBH, ψ) ρ(m₁)ρ(m₂). It is valid for f_PBH ≲ 0.1, and S is a suppression factor. Local enhancement uses an NFW Milky Way halo with ρ_DM(r_⊙) = 7.9×10⁻³ M_⊙/pc³, r_⊙ ≃ 8.0 kpc and r₀ = 15.6 kpc, giving δ(r) between 1 and 2×10⁵ (Eqs. 53–54). Distance reach for an optimally oriented source at fixed SNR (Eq. 56, as read from the plain-text extraction): d = (5/24)^{1/2} π^{-2/3} (GM_c)^{5/6} Δt^{-(k-1)/4} {(2/SNR) ∫_{f_min}^{f_max} df [f^{-7/3}/S_h^noise(f)]^k}^{1/(2k)}, with k = 1 for linear and k = 2 for quadratic detectors and f_max = min(f_ISCO, detector maximum). Check the exponent placement against the PDF. — [Review 2025 §4.1.2](https://arxiv.org/abs/2501.11723), building on [Franciolini, Maharana, Muia 2022, arXiv:2205.02153](https://arxiv.org/abs/2205.02153)
- *SGWB from PBH binaries.* It is a late-time background (z ~ O(10)) and so is **not** subject to the ΔN_eff bound, but for realistic f_PBH ≈ 1 its amplitude still lies below the BBN level. The inspiral tail goes as Ω ∝ f^{2/3}, h_c ∝ f^{-2/3}. Table 1 gives h_c,sto ≈ 5×10⁻³¹ (f_ISCO/GHz)^{-1.07}. The background is expected to be strongly anisotropic (Poisson). — [Review 2025 §4.1.2, Table 1](https://arxiv.org/abs/2501.11723)
- *PBH hyperbolic encounters.* Peak frequency f_p = (1/2π)√(GM(e+1)/r_p³) ≃ 1.6 GHz (10⁻⁵ M_⊙/M)(R_S/r_p)^{3/2}√((e+1)/2) (Eq. 59). Peak strain h_p = 3.6×10⁻²⁵ [4q/(1+q)²][G(e)/G(1)](M/10⁻⁵ M_⊙)^{5/3}(f_p/1.6 GHz)^{2/3}(1 Mpc/D), with G(e) = e + 2/(e+1)^{1/3} (Eq. 60). Burst duration ~1/f_p (Eq. 61). — [Review 2025 §4.1.2](https://arxiv.org/abs/2501.11723)
- *GW memory from PBH mergers.* h(f) ~ 1.2×10⁻²⁴ s (f_ISCO/f)(M/10⁻⁵ M_⊙)(D/kpc)⁻¹ (Eq. 69). It extends to f ≪ f_ISCO. — [Review 2025 §4.1.3](https://arxiv.org/abs/2501.11723)
- *Black-hole superradiance (boson clouds around light PBHs).* Nearly monochromatic emission at f_GW ~ ω_R/π ~ 5 MHz (μ/10⁻⁸ eV) ~ 5 GHz (10⁻⁵ M_⊙/m_PBH)(Gm_PBHμ/0.1) (Eq. 73). Peak strain for scalars h_S ≈ 5×10⁻³⁰ (m_PBH/10⁻⁵ M_⊙)(α/0.1)⁷((χ_i−χ_f)/0.5)(kpc/D) (Eq. 74). For vectors and tensors h_{V,T} ≈ 10⁻²⁶ (same scalings, with (α/0.1)⁵) (Eq. 75). Here α ≡ Gm_PBHμ. Signal duration τ_S ≈ 0.13 yr (m/10⁻⁵ M_⊙)(0.1/α)¹⁵(0.5/Δχ) and τ_{V,T} ≈ 0.17 s (m/10⁻⁵ M_⊙)(0.1/α)¹¹(0.5/Δχ) (Eqs. 76–77). Frequency drift ḟ/f² ≃ 3×10⁻²⁰ (α/0.1)¹⁷ (Eq. 78). Bosenova bursts (QCD axion) occur at f_bn ≈ 30 MHz (16/c_bn)(m/10⁻⁵ M_⊙)⁻¹(...)², with h(f) ≈ 10⁻²⁷ s × (scalings) (Eqs. 79–80). — [Review 2025 §4.1.5](https://arxiv.org/abs/2501.11723), citing [Brito, Cardoso, Pani, arXiv:1501.06570](https://arxiv.org/abs/1501.06570)
- *Known astrophysics.* NS-merger post-merger emission is ≲ (1–5) kHz, with h(f) ≲ 10⁻²⁶ s. A QCD first-order phase transition (FOPT) inside an NS merger would peak at f ≈ 0.6 MHz (v_w = 0.1, τ = 1 ms) with h_c,sto ≃ 1.5×10⁻²⁴ v_f² (100 Mpc/D) (Eqs. 36–39). Photon-sphere Gertsenshtein conversion around SMBH accretion disks gives f ≃ 3.3×10¹⁹ Hz and Ω_GW ≃ 2×10⁻¹² ξ (ξ ≲ 1). The thermal plasma of the Sun gives a spectrum over 10¹²–10¹⁹ Hz peaking near 10¹⁸ Hz, with h_c ~ 10⁻⁴². — [Review 2025 §4.1.1, Table 1](https://arxiv.org/abs/2501.11723); MHz NS-merger FOPT from [Casalderrey-Solana et al., arXiv:2210.03171](https://arxiv.org/abs/2210.03171); SMBH disks from [Saito, Soda, Yoshino 2021, arXiv:2106.05552](https://arxiv.org/abs/2106.05552)

**Early-Universe stochastic sources (Review 2025 §4.2, Table 2)**
- *Inflation.* The scale-invariant vacuum signal is Ω_GW,0 ≲ 10⁻¹⁶. The extra-species (axion–gauge field) mechanism gives Ω_GW,0 ≃ 3.6×10⁻⁹ Ω_rad,0 (H⁴/M_Pl⁴) e^{4πξ}/ξ⁶, with ξ = φ̇/(2f_aH) ≲ 4.7 in the perturbative regime, reaching about 10⁻¹⁰ (Eq. 89). The broken-spatial-reparametrization EFT gives P_T = (2H²/(π²M_Pl²c_T³))(k/k*)^{n_T} with n_T = −2ε + (2/3)m_h²/H² (Eq. 91). Second-order induced GWs give Ω_GW,0 ~ Ω_rad,0 P_ζ², i.e. ~10⁻⁹ for P_ζ ≲ 10^{-2.5}. Mode–frequency map: ln(f/10⁻¹⁸ Hz) ≃ N_CMB − N (Eq. 88). — [Review 2025 §4.2.1](https://arxiv.org/abs/2501.11723)
- *Preheating.* Lattice results: for V ∝ φ², f ≃ 10⁸–10⁹ Hz and Ω_GW,0 ≃ 10⁻¹²–10⁻¹¹ (q* ∈ 10⁴–10⁶); for V ∝ φ⁴, f ≃ 10⁷–10⁸ Hz and Ω ≃ 10⁻¹³–10⁻¹¹. Gauge-field (axial) preheating can reach Ω ~ 10⁻⁷–10⁻⁶. Oscillons: f ≃ 10⁶ Hz × X × (m/10¹² GeV)^{5/6} with X ≈ 10–10³ (Eq. 96, as extracted), and Ω_GW,0 ≃ 10⁻¹³–10⁻¹⁰. Table 2 summary: (10⁶–10⁹) Hz, Ω ≲ 10⁻¹⁰, h_c ≲ 10⁻²⁹ (MHz/f). — [Review 2025 §4.2.2](https://arxiv.org/abs/2501.11723), lattice parameter study [Figueroa & Torrenti 2017, arXiv:1707.04533](https://arxiv.org/abs/1707.04533)
- *Thermal background, the "cosmic gravitational microwave background" (CGMB).* It is an SM-guaranteed signal peaking near 100 GHz. Ω_GW,0(f) ≃ (1440√10/2π²)[g*s(T_end)^{1/3}/g*s(T_max)^{5/6}] Ω_rad,0 (f³/T₀³)(T_max/M_Pl) η̂(T_max, ...) (Eq. 97). The peak is at **f_peak ≈ 79.8 GHz (106.75/g*s(T_max))^{1/3}** (Eq. 98). The amplitude is linear in T_max and approaches the dark-radiation bound only for T_max ~ few×10¹⁹ GeV. CMB tensor limits require T_max < 6.6×10¹⁵ GeV (106.75/g*s)^{1/4} for slow-roll inflation with instantaneous reheating. Table 2: f_peak ~ (10–100) GHz, Ω ≲ 10⁻⁶. Double-graviton corrections were added in 2024. — [Review 2025 §4.2.3](https://arxiv.org/abs/2501.11723); [Ringwald, Schütte-Engel, Tamarit 2021, arXiv:2011.04731](https://arxiv.org/abs/2011.04731); [Ghiglieri, Laine, Schütte-Engel, Speranza 2024, arXiv:2401.08766](https://arxiv.org/abs/2401.08766)
- *First-order phase transitions (sound waves).* f_peak ≃ 260 MHz × (H*R*)⁻¹ (T*/10¹⁵ GeV)(g*/100)^{1/6} (Eq. 99), with 1 ≲ (H*R*)⁻¹ ≲ 10⁴ and R* ≈ (8π)^{1/3}v_w/β. Ω_GW,0(f) ≃ 2H*R*[1 − (1+2H*τ_v)^{-1/2}](../../../../research_notes/Gravitational%20waves%20and%20gravity%20beams/100/g%2A)^{1/3} K² Ω̃_GW Ω_rad,0 S(f/f_peak), with S(x) = x³[7/(4+3x²)]^{7/2} (Eqs. 100–101). Ω̃_GW ≈ 0.058 (also quoted as O(10⁻²)), K = κα/(1+α), and κ ≈ α/(0.73 + 0.083√α + α) for near-luminal walls. Result: Ω_GW,0 ≲ 10⁻⁷, reachable only with strong supercooling. The frequency is capped at ~10² GHz. Caveat: the plain-text extraction renders τ_v as "4R*/(3K)", but the standard form is τ_v ≈ R*/Ū_f with Ū_f ≈ √(3K/4); check the PDF. — [Review 2025 §4.2.4](https://arxiv.org/abs/2501.11723), [Hindmarsh et al. 2017, arXiv:1704.05871](https://arxiv.org/abs/1704.05871)
- *Cosmic strings (Nambu–Goto loops).* Plateau **Ω_GW,0^plateau ≈ 8.04 Ω_rad,0 √(Gμ/Γ)** with Γ ≃ 50 (Eq. 102). The high-frequency cutoff is set by the network formation temperature; T_max ≲ 10¹⁶ GeV implies f_Δ ≲ 10⁹ Hz (the extracted text writes "GeV", evidently a typo for Hz). Bounds: CMB Gμ ≲ 10⁻⁷, PTA Gμ ≲ 10⁻¹¹ (the Fig. 4 description says < 10⁻¹⁰; the review is internally inconsistent here). Metastable strings evade the PTA bound and are limited only by BBN, Gμ ≲ 10⁻⁴ (up to ~10⁻³ above the LIGO band per the Fig. 4 text), with the plateau spanning 10³–10⁹ Hz. Global strings: Ω ≃ 4×10⁴ (Gμ)² Ω_rad,0 (Eq. 103). Gauge textures: f ~ 10¹¹ Hz (v/M_Pl), Ω ~ 2×10⁻⁴ (v/M_Pl)⁴ (Eq. 104). — [Review 2025 §4.2.5](https://arxiv.org/abs/2501.11723); plateau formula from [Auclair et al. 2020, arXiv:1909.00819](https://arxiv.org/abs/1909.00819); metastable strings [Buchmuller, Domcke, Schmitz 2020, arXiv:2009.10649](https://arxiv.org/abs/2009.10649)
- *Evaporating PBHs (Hawking gravitons).* The spectrum lies at 10¹³–10²² Hz. Maximum Ω_GW,0 ≈ 10^{-7.5} for PBHs evaporating just before BBN (m_i ≲ 10⁹ g), and ≈ 10^{-6.5} for 10⁹–10¹⁷ g with the peak at 10¹⁸–10²² Hz. GUT-scale PBHs peak near 10¹⁵ Hz with Ω ~ 10⁻⁸. — [Review 2025 §4.2.6](https://arxiv.org/abs/2501.11723); [Anantua, Easther, Giblin 2009, arXiv:0812.0825](https://arxiv.org/abs/0812.0825); [Dong, Kinney, Stojkovic 2016, arXiv:1511.05642](https://arxiv.org/abs/1511.05642)
- *Lab generation is hopeless.* A 50 cm piezo crystal at 300 MHz gives h ≈ 10⁻³⁸ at 10 m. A 10⁴ J/m³ cavity gives h ≈ 10⁻⁴³ (10 m/r) at 10¹⁰ Hz. LEP-2 gives h ≈ 10⁻⁴³ at 100 km. Counter-propagating lasers give h ≈ 5.2×10⁻³⁸ (λ/10 µm)²(τ/1 ps)(I/10²³ W cm⁻²)(10 cm/r) (Eq. 105). — [Review 2025 §4.3](https://arxiv.org/abs/2501.11723)

### Inferences
- Converted to strain with Maggiore/Moore normalization (see Q4), the BBN bound Ω h² < 1.1×10⁻⁶ gives h_c < 1.33×10⁻²¹ (1 Hz/f): 1.3×10⁻²⁶ at 100 kHz, 1.3×10⁻²⁷ at 1 MHz, 1.3×10⁻³⁰ at 1 GHz, 1.3×10⁻³³ at 1 THz. The CMB-only bound (2.9×10⁻⁷) gives 6.8×10⁻²² (1 Hz/f). (Derived numerically.) Using the Review 2025 h_c normalization, divide these by √2.
- The superradiance annihilation frequency is f_GW ≈ 2μc²/h = 4.84×10¹⁴ Hz (μ/1 eV), i.e. 4.8 MHz at μ = 10⁻⁸ eV. This matches the review's Eq. 73 normalization (derived; it neglects O(α²) binding corrections).
- For a Rust signal-model crate, a natural split is: (a) deterministic waveforms (PBH 0-PN chirp with ḟ, N_cycles and τ; hyperbolic-burst peak scaling; monochromatic superradiance with ḟ drift), all with 1/D amplitude scaling; (b) stochastic Ω_GW(f) templates (PT sound-wave S(x), string plateau with cutoff, CGMB via Eq. 97 with η̂ tabulated from Ringwald et al.), each clipped against the ΔN_eff bound.

### Gaps
- I did not extract η̂(T, k) for the CGMB. It must be taken from Ringwald et al. 2021, Ghiglieri & Laine 2015 or Ghiglieri et al. 2024 (tabulated rates), so I have no closed-form fit.
- 2023–2025 primary papers on UHF cosmic strings (for example Servant & Simakachorn) were not opened. Cosmic-string numbers above come from the review only.
- The review flags the Nambu–Goto versus field-theory string debate (~20 years, unresolved) as changing the predicted amplitude by large factors. I found no quantitative resolution.
- Eq. 96 (oscillon frequency) exponent placement was ambiguous in the text extraction. Verify against the PDF.

## Q2. Main detection concepts, best published sensitivities, and the formulas behind them

### Takeaway
No concept, built or proposed, reaches realistic sources: they fall "short by at least several orders of magnitude" (Review 2025 §6).
- **Published GW limits** exist from: the Holometer (1–13 MHz); the 100 MHz synchronous-recycling interferometers; QUEST (13–80 MHz); graviton–magnon reinterpretations (8–14 GHz); recast axion haloscope data (ABRACADABRA-10cm at 10 kHz–5 MHz; CAPP-12T at 5.3 GHz); light-shining-through-a-wall (LSW) and helioscope photon data (optical and X-ray); and astrophysical conversion (neutron stars, Galactic field, CMB/radio).
- **Under construction or R&D, no GW limits yet:** the Levitated Sensor Detector (LSD) and MAGO/SRF heterodyne cavities.
- **Mechanical vs electromagnetic:** mechanical coupling (bars, BAW, levitated sensors, MAGO, magnetic Weber bars) favors kHz–GHz; magnetic (Gertsenshtein) conversion favors ≳ GHz.

### Cited Findings

**General detector formalism (Review 2025 §2.3)**
- Noise-equivalent strain PSD: S_h^noise(f) = S_n̄(f)|T_n(f)|²/|T_h(f)|² (Eq. 16), with T_h the strain→readout transfer function. Optimal SNR is SNR_lin = [2Δt ∫₀^∞ df S_h/S_h^noise]^{1/2} for observables linear in h (Eq. 22). For quadratic observables (e.g. photon counting), SNR_quad = [2Δt ∫ df (S_h/S_h^noise)²]^{1/2} (Eq. 23). — [Review 2025 §2.3](https://arxiv.org/abs/2501.11723)
- Cross-correlating N linear detectors gives SNR ≃ [N(N−1) t_int ∫ df Γ²(f)(S_h/S_h^noise)²]^{1/2} (Eq. 24), with overlap reduction function Γ. In Ω form: SNR ≃ (3H₀²/4π²)[N(N−1) t_int ∫ df Γ²Ω_GW²/(f⁶ S_h^noise²)]^{1/2} (Eq. 26). — [Review 2025 §2.3.1](https://arxiv.org/abs/2501.11723)
- Persistent monochromatic signal: h ≳ SNR × (S_h^noise(f_g)/t_int)^{1/2} for linear detectors (Eq. 31). For quadratic detectors it is h ≳ SNR^{1/2}(S_h^noise/t_int)^{1/2} (Eq. 32), or SNR^{1/2}(S_h^noise)^{1/2}(Δf_det/t_int)^{1/4} when the detector response is broader than 1/t_int (Eq. 33). Transients: S_h^lin ≳ SNR² S_h^noise/(ΔfΔt) and S_h^quad ≳ SNR² S_h^noise (ΔfΔt)^{-1/2} (Eqs. 27–28); burst strain h ~ √(S_h/τ) (Eq. 29). — [Review 2025 §2.3.2–2.3.3](https://arxiv.org/abs/2501.11723)
- Resonant scanning gains nothing at large Q for a stochastic background: the Q gain is cancelled by less time per bin. Running at a fixed frequency gives Ω sensitivity ∝ Q^{-1/2}. Sky-averaged sensitivity for wide-field detectors is about 10× worse than the optimal-orientation curves usually plotted. For narrow-field LSW experiments the loss is much larger. — [Review 2025 §3](https://arxiv.org/abs/2501.11723)

**Interferometers and resonant masses (Review 2025 §5.1)**
- 0.75 m synchronous-recycling interferometers at 100 MHz (0.5 W Nd:YAG, kHz bandwidth): S_n^{1/2} ~ 10⁻¹⁶ Hz^{-1/2}, h_c,sto ≲ 7×10⁻¹⁴. A design study shows 10⁻²⁰ Hz^{-1/2} is possible at 100 MHz with 2 kHz bandwidth, but the approach is not competitive above 1 GHz. — [Review 2025 §5.1.2](https://arxiv.org/abs/2501.11723); [Akutsu et al. 2008, arXiv:0803.4094](https://arxiv.org/abs/0803.4094); [Nishizawa et al. 2008, arXiv:0710.1944](https://arxiv.org/abs/0710.1944). Ejlli et al. 2019 instead quote h_c ≤ 1.4×10⁻¹² at 100 MHz for the same experiment — [Ejlli et al. 2019](https://arxiv.org/abs/1908.00232). The two numbers conflict.
- Holometer (Fermilab): two co-located 40 m (39 m arm) power-recycled Michelsons. Cross-correlated sensitivity S_n^{1/2} ≃ 10⁻²¹ Hz^{-1/2} over 1–13 MHz from a 10³ h dataset. Ejlli et al. quote h_c < 2.5×10⁻¹⁸ at 1 MHz falling to 2.4×10⁻¹⁹ at 13 MHz (3σ). A 704 h analysis found no harmonic sources (cosmic-string loops, eccentric binaries) at 1–25 MHz. — [Review 2025 §5.1.2](https://arxiv.org/abs/2501.11723); [Chou et al. (Holometer) 2017, arXiv:1611.05560](https://arxiv.org/abs/1611.05560); [Martinez & Kamai 2020, arXiv:2010.06118](https://arxiv.org/abs/2010.06118); [Ejlli et al. 2019](https://arxiv.org/abs/1908.00232)
- QUEST (Cardiff): twin co-located power-recycled Michelsons designed for shot-noise-limited 2×10⁻¹⁹ m/√Hz over 1–200 MHz. A 10⁴ s coincident run gives strain ≈ 3×10⁻²⁰ Hz^{-1/2} at 40 MHz and limits on correlated length fluctuations over 13–80 MHz. These are described as "the first broadband constraints for a stochastic gravitational wave background at these frequencies". Plans: 200 MHz bandwidth and 100× better sensitivity. Published PRL 2025, DOI 10.1103/61j9-cjkk. — [Patra et al., arXiv:2410.09175](https://arxiv.org/abs/2410.09175); [Review 2025 §5.1.2](https://arxiv.org/abs/2501.11723)
- Spheres and bars: MiniGRAIL S_n^{1/2} ≃ 1.5×10⁻²⁰ at 2942.9 Hz (5 K). Schenberg 1.1×10⁻¹⁹ at ~3.2 kHz. AURIGA ≃ 10⁻²⁰ near 900 Hz. Resonant masses are unlikely to go far beyond kHz. — [Review 2025 §5.1.3](https://arxiv.org/abs/2501.11723)
- The quantum Cramér–Rao bound links stored optical energy to the bandwidth × sensitivity product, so LIGO-level broadband MHz interferometers are judged not viable. — [Review 2025 §5.1.2](https://arxiv.org/abs/2501.11723)

**Levitated sensors (Review 2025 §5.2.1)**
- Concept: a dielectric particle or disk is trapped at an antinode of the cavity standing wave. A GW changes the cavity length, displacing the particle from the antinode, and the trap frequency is tuned to the GW. The sensitivity is set by thermal noise of the trapped particle and ultimately by photon-recoil heating, so it *improves* with frequency. — [Arvanitaki & Geraci 2013, arXiv:1207.5320](https://arxiv.org/abs/1207.5320); [Review 2025 §5.2.1](https://arxiv.org/abs/2501.11723)
- LSD: a 1 m prototype Michelson with Fabry–Pérot arms is under construction at Northwestern, targeting S_n^{1/2} ≲ 10⁻¹⁹ Hz^{-1/2} at 10 kHz and ~10⁻²¹ at 100 kHz. A 10 m room-temperature instrument is projected at ~10⁻²⁰ (10 kHz) and ~10⁻²² (100 kHz); a cryogenic 100 m instrument at more than 10× better. The target band is 10–100 kHz (original concept 50–300 kHz), and the main physics target is superradiance. — [Aggarwal et al., arXiv:2010.13157](https://arxiv.org/abs/2010.13157); [Review 2025 §5.2.1, Table 6](https://arxiv.org/abs/2501.11723)
- Hardware milestones: optical trapping of NaYF₄ hexagonal prisms (PRL 129, 053604, 2022) and a 10× increase in levitated mass (SPIE 2024). As of an April 2026 APS abstract the LSD is "still in development", with an expected thermal-noise-limited h ~ 10⁻²¹–10⁻²³. **No GW limits from levitated sensors have been published** as of this search. — [Winstone et al., arXiv:2204.10843](https://arxiv.org/abs/2204.10843); [SPIE 2024 proceeding](https://proceedings.spiedigitallibrary.org/conference-proceedings-of-spie/PC12902/PC129020K/Detecting-high-frequency-gravitational-waves-with-optically-levitated-high-aspect/10.1117/12.3013033.full); [APS 2026 abstract](https://meetings-archive.aps.org/smt/2026/apr-v91/7/)
- A general-relativistic treatment (PRD, accepted) finds the strain response depends strongly on trap position and is maximized near the input mirror. — [PRD accepted, DOI 10.1103/mhz8-h774](https://journals.aps.org/prd/accepted/10.1103/mhz8-h774)
- Related concepts: levitated membranes with Q > 10¹⁰–10¹² (DESY R&D). A superfluid-⁴He detector claims S_n^{1/2} < 8×10⁻¹⁹ at ~100 kHz. A superconducting levitated sphere projects S_n^{1/2} ~ 10⁻¹⁹ Hz^{-1/2} over 10 kHz–1 MHz. — [Review 2025 §5.2.1, §5.2.4](https://arxiv.org/abs/2501.11723); [Vadakkumbatt et al. 2021, arXiv:2107.00120](https://arxiv.org/abs/2107.00120); [Carney et al. 2024, arXiv:2408.01583](https://arxiv.org/abs/2408.01583)

**Bulk acoustic wave (BAW) resonators (Review 2025 §5.2.2)**
- Cryogenic quartz BAW resonators reach Q up to 8×10⁹ on longitudinal overtones spanning 5–700 MHz (1–700 MHz quoted overall), with piezoelectric readout through SQUIDs. The estimated sensitivity is S_n^{1/2} ≃ 2×10⁻²² Hz^{-1/2}, depending on mode geometry. The 2014 sensitivity paper carries a 2023 erratum (PRD 108, 129901). — [Goryachev & Tobar 2014, arXiv:1410.2334](https://arxiv.org/abs/1410.2334); [Review 2025 §5.2.2](https://arxiv.org/abs/2501.11723)
- MAGE (UWA): first stage uses two near-identical quartz BAWs, with sensitivity down to 6.6×10⁻²¹ strain/√Hz in several narrow MHz bands. The second detector allows local events to be vetoed, and a stated goal is to identify the source of the earlier rare events. A Milan program (BAUSCIA) and MAGE aim for networks of O(10) BAWs with O(100) modes. — [Campbell, Goryachev, Tobar 2023, Sci. Rep. 13, 10638, arXiv:2307.00715](https://arxiv.org/abs/2307.00715); [Review 2025 §5.2.2](https://arxiv.org/abs/2501.11723)

**Deformation of microwave cavities: MAGO / heterodyne SRF (Review 2025 §5.2.3)**
- Principle: two nearly degenerate modes ω₁, ω₂ with power pumped into mode 1. A GW at ω_g transfers power into mode 2 when ||ω₂−ω₁| − ω_g| ≲ Δω₂. Below a GHz the coupling is dominated by wall deformation (the mechanical channel) rather than direct EM coupling. — [Review 2025 §5.2.3](https://arxiv.org/abs/2501.11723); original MAGO [Ballantini et al., arXiv:gr-qc/0502054](https://arxiv.org/abs/gr-qc/0502054)
- MAGO 2.0 theory: noise-equivalent S_n^{1/2} ~ 10⁻²¹ Hz^{-1/2} for 100 kHz ≲ 2πf ≲ GHz with resonant (narrow) operation. With an overcoupled signal mode (broadband) the projection is S_n^{1/2} ≲ 10⁻¹⁸ Hz^{-1/2} over two decades around 100 kHz. — [Berlin et al. 2023, PRD 108, 084058, arXiv:2303.01518](https://arxiv.org/abs/2303.01518); [Review 2025 §5.2.3](https://arxiv.org/abs/2501.11723)
- Experimental revival (DESY/UHH/Fermilab-SQMS):
  - First characterisation of the 20-year-old spherical 2-cell niobium MAGO cavity (CQG 42, 2025).
  - A Jan 2026 paper describes a seismically isolated 1.8 K cryostat for a planned 10 kHz–100 MHz search.
  - A July 2026 paper reports cryogenic tests at ~2 K at Fermilab and DESY, with ~11 kHz mode splitting after room-temperature tuning, but **mechanical Q factors "well below commonly assumed theoretical values"**.
  - The first goal is a weak bound in the 10–100 kHz range.
  - GravNet (INFN Frascati, IFAE, Bonn, Mainz) pursues networks.
  - Sources: [Fischer et al., arXiv:2411.18346](https://arxiv.org/abs/2411.18346); [Wenskat et al., arXiv:2601.18719](https://arxiv.org/abs/2601.18719); [Dokuyucu et al., arXiv:2607.02476](https://arxiv.org/abs/2607.02476); [Review 2025 §5.2.3](https://arxiv.org/abs/2501.11723)

**Magnetic Weber bars (Review 2025 §5.2.4)**
- A persistent-mode superconducting magnet acts as a resonant mass. Coil deformation changes the field by O(hB₀), read by a pickup loop and SQUID. The ADMX-EFR MRI magnet gives an estimated broadband S_n^{1/2} ~ 10⁻²⁰ Hz^{-1/2} from a few kHz to ~10 MHz, with a peak of ~10⁻²² at ~kHz on a mechanical resonance. DMRadio-GUT is also considered. — [Domcke, Ellis, Rodd 2024, arXiv:2408.01483](https://arxiv.org/abs/2408.01483); [Review 2025 §5.2.4](https://arxiv.org/abs/2501.11723)

**Electromagnetic oscillators: microwave cavities and LC haloscopes (Review 2025 §5.3)**
- In curved spacetime the GW acts as an effective current j_eff^μ = ∂_ν(−½hF₀^{μν} + h^ν_α F₀^{μα} − h^μ_α F₀^{να}) (Eq. 108). This is analogous to the axion current (∂_νa)F̃^{μν}, so axion experiments double as HFGW detectors. Surface currents on the boundary of the magnetized volume must be included. — [Review 2025 §5.3](https://arxiv.org/abs/2501.11723); [Domcke, Garcia-Cely, Lee, Rodd 2024, JHEP 03 (2024) 128, arXiv:2306.03125](https://arxiv.org/abs/2306.03125)
- Microwave cavity:
  - Induced field Ẽ_h ≃ ηQ(2πfL)B₀h̃ (Eq. 109), with η ≃ 0.1.
  - Signal power P_sig ≃ ½QL⁵(2πf)³B₀²η²S_hΔf (Eq. 110).
  - Thermal noise P_noise ≃ k_BT_sysΔf.
  - Result: **S_h^noise ≃ 2k_BT_sys/[(2πf)³η²QB₀²L⁵]** (Eq. 113, natural units).
  - The cavity literature quotes h₀ ~ 10⁻²²–10⁻²¹ at O(GHz) — [Review 2025 §5.3.1](https://arxiv.org/abs/2501.11723); [Berlin et al. 2022, PRD 105, 116011, arXiv:2112.11465](https://arxiv.org/abs/2112.11465)
- Benchmark parameters (Review 2025 Table 3; f [GHz], Q, B₀ [T], L [m], T_sys [K]):
  - ADMX (0.65–1.02, 8×10⁴, 7.5, 0.51, 0.6)
  - HAYSTAC (5.6–5.8, 3×10⁴, 9, 0.13, 0.13)
  - CAPP (1.6–1.65, 4×10⁴, 7.3, 0.15, 1.2)
  - ORGAN (15–16 and 26–27, 10⁴, 11.5, 0.023, 5.3)
  - SQMS (1–2, 10⁶, 5, 0.46, 1)
  - cubic cavities 1/2/3 (0.1/1/10 GHz)
  - RADES-BabyIAXO (0.25–0.33 and 2.5–3.4, 10⁵, 2, 0.5–5, 4.6)
  - ABRACADABRA (10⁻⁴–0.002, 1, 1, 0.096, 0.5)
  - SHAFT (3×10⁻⁶–0.003, 1, 1.51, 0.046, 4.2)
  - ADMX SLIC (0.043, 3×10³, 7, 0.2, 20)
  - BASE (4×10⁻⁴, 4×10⁴, 1.85, 0.025, 5.7)
  - WISPLC (3×10⁻⁵–0.005, 10⁴, 14, 0.29, 4)
  - DMRadio-m³ (0.005–0.2, 10⁵, 4, 1.3, 0.02)
  - DMRadio-GUT (10⁻⁴–0.03, 2×10⁷, 16, 2.2, 0.01)
  - — [Review 2025 Table 3](https://arxiv.org/abs/2501.11723)
- Low-mass LC haloscopes:
  - Induced flux Φ̃_h ≃ ηB₀(2πfL)²L²h̃ (Eq. 114).
  - Resonant readout: **S_h^noise ≃ 2T_sys/[(2πf)⁵η²B₀²L⁷Q]** (Eq. 117).
  - Broadband SQUID-limited: S_h^noise ≃ (10⁻⁶Φ₀)²Hz⁻¹/[2π⁴η²f⁴B₀²L⁷L_sq] (Eq. 118), with α = 1/√2 and L_sq ≃ 1 nH.
  - Recast axion sensitivities for coherent GWs: h₀ ~ 10⁻⁹ (2 MHz, ABRA), ~10⁻¹⁶ (40 MHz, ADMX SLIC), ~10⁻¹⁵ (6 MHz, WISPLC), ~10⁻²¹ (100 MHz, DMRadio-m³).
  - Sources: [Domcke, Garcia-Cely, Rodd 2022, PRL 129, 041101, arXiv:2202.00695](https://arxiv.org/abs/2202.00695); [Review 2025 §5.3.2](https://arxiv.org/abs/2501.11723)
- Proposals with mixed status: plasma cavities (ALPHA) project h₀ ~ 3×10⁻²⁴–10⁻²² at 10–50 GHz. Claims that the AC–DC cross-term gives extra sensitivity are **rejected** by the review, which says the linear-signal sensitivity is independent of the background field. — [Review 2025 §5.3.1](https://arxiv.org/abs/2501.11723)

**Photon (re)generation: LSW, helioscopes, dielectric haloscopes (Review 2025 §5.4)**
- Induced field Ẽ_h, B̃_h ≃ 2πfLB₀h̃ for 2πfL ≫ 1 (Eq. 119). Then ⟨S⟩ ≃ ∫df (2πfLB)²S_h (Eq. 120), S_P,sig ≃ A(2πfLB)²S_h (Eq. 121), and the photon-count PSD is S_N,sig ≃ 2πεfAL²B²Δt S_h (Eq. 122).
- Single-photon (dark-count Γ_D) limit: **S_h^noise ≃ (Γ_D/Δf)/[2εA(2πf)L²B²]** (Eq. 125). Thermal or heterodyne readout: **S_h^noise ≃ k_BT_sys/(8π²f²B²L²A)** (Eq. 126). Detection also requires N_γ,sig ≥ 1, which dominates for short transients such as PBH mergers. Low-frequency waveguide cutoff near f ~ L/(4πA) (Eq. 123, as extracted). — [Review 2025 §5.4.1](https://arxiv.org/abs/2501.11723)
- Table 4 parameters (f, B₀, L, A, Γ_D, ε):
  - OSQAR II: (0.3–1)×10⁶ GHz, 9 T, 14.3 m, 5×10⁻⁴ m², 1.14×10⁻³ Hz, 0.9
  - ALPS II: (0.3–1)×10⁶ GHz, 5.3 T, 106 m, 2×10⁻³ m², 10⁻⁶ Hz, 0.75
  - JURA: 13 T, 960 m, 8×10⁻³ m², 10⁻⁶ Hz, 1
  - CAST: (0.5–1.6)×10⁹ GHz, 9 T, 9.26 m, 2.9×10⁻³ m², 1.5×10⁻⁴ Hz, 0.7
  - BabyIAXO: (0.25–2)×10⁹ GHz, 2.5 T, 10 m, 0.77 m², 10⁻³ Hz, 1
  - IAXO: 2.5 T, 20 m, 3.08 m², 10⁻⁴ Hz, 1
  - — [Review 2025 Table 4](https://arxiv.org/abs/2501.11723)
- A IAXO-type magnet in the GHz range could reach h_c ≲ 10⁻²² with heterodyne receivers or ≲ 10⁻²⁵ with single-photon detectors near 4×10¹⁰ Hz. — [Ringwald, Schütte-Engel, Tamarit 2021, arXiv:2011.04731](https://arxiv.org/abs/2011.04731), via [Review 2025 §5.4.1](https://arxiv.org/abs/2501.11723)
- Dielectric haloscopes (MADMAX), fully resonant mode: S_h^{noise 1/2} ~ 10⁻²² Hz^{-1/2} (1 m²/L²)(10 T/B₀)(43/N_d)(10 GHz/f) with N_d ≲ 50 disks (Eq. 127). A hybrid configuration is favored. DALI is an alternative (6–60 GHz). — [Domcke, Ellis, Kopp 2024, arXiv:2409.06462](https://arxiv.org/abs/2409.06462); [Review 2025 §5.4.2](https://arxiv.org/abs/2501.11723)

**Other EM, spin, atomic and quantum concepts (Review 2025 §5.5, §5.7)**
- Graviton–magnon resonance (a GW couples to spins of a ferromagnetic insulator; readout via a cavity-coupled magnon mode): reinterpreted ferromagnetic-haloscope data give S_n^{1/2} ≲ 7.6×10⁻²² Hz^{-1/2} at 14 GHz and ≃ 1.2×10⁻²⁰ at 8.2 GHz. Expressed as stochastic h_c these are ≤ 9.1×10⁻¹⁷ (14 GHz) and ≤ 1.1×10⁻¹⁵ (8.2 GHz), as quoted by Ejlli et al. Larger samples (size ~ λ_GW) could reach h₀ ~ 10⁻²⁰. — [Ito, Ikeda, Miuchi, Soda 2020, arXiv:1903.04843](https://arxiv.org/abs/1903.04843); [Ito & Soda 2023, arXiv:2212.04094](https://arxiv.org/abs/2212.04094); [Review 2025 §5.7.2](https://arxiv.org/abs/2501.11723); [Ejlli et al. 2019](https://arxiv.org/abs/1908.00232)
- One-electron quantum cyclotron (Penning trap, the "electron proposal"): the transition rate grows with the spatial size of the electron wavefunction, a GW-specific effect. Existing dark-photon setups give h₀ ~ 3.8×10⁻²⁰ near 100 GHz. With large-wavefunction state preparation h₀ ≲ 6.9×10⁻²³, scannable over 20–200 GHz. — [Ito & Kitano 2024, JCAP 04 (2024) 068, arXiv:2309.02992](https://arxiv.org/abs/2309.02992); [Review 2025 §5.7.4](https://arxiv.org/abs/2501.11723)
- Rydberg atoms (EIT, heterodyne with a "local laser"): h₀ ≲ few×10⁻²⁰ near 30 GHz, ≲ 10⁻²³ if quantum-noise limited, and potentially 10³× better with entanglement. — [Kanno, Soda, Taniguchi 2023, arXiv:2311.03890](https://arxiv.org/abs/2311.03890); [Review 2025 §5.7.5](https://arxiv.org/abs/2501.11723)
- Optical frequency modulation with atomic-clock techniques and an "optical rectifier": h₀ ≲ 10⁻¹⁷–10⁻²¹ over 10 kHz–10 GHz for 1 s integration, under optimistic assumptions. — [Bringmann, Domcke, Fuchs, Kopp 2023, PRD 108, L061303, arXiv:2304.10579](https://arxiv.org/abs/2304.10579)
- Resonant polarization rotation: the Cruise–Ingley waveguide ring gives S_n^{1/2} ≲ 10⁻¹⁴ Hz^{-1/2} at 100 MHz with two correlated detectors (h_c ≤ 1.4×10⁻¹⁰ per Ejlli et al.). An ALPS II optical-cavity polarimetric version is proposed for 0.1 MHz–0.1 THz. — [Review 2025 §5.5.3](https://arxiv.org/abs/2501.11723); [García-Cely, Marsili, Ringwald, Spector 2025, arXiv:2501.08382](https://arxiv.org/abs/2501.08382)
- High-energy pulsed lasers (inverse Gertsenshtein at 2× the laser frequency): h₀ ≲ 10⁻²⁰ over 10¹³–10¹⁹ Hz, down to 10⁻²⁶ with next-generation lasers. Atomic-electric-field conversion targets 10²⁰–10²⁴ Hz using neutrino detectors such as JUNO. — [Review 2025 §5.5.1–5.5.2](https://arxiv.org/abs/2501.11723)
- Neutron proposals: I found **no** neutron-interferometer HFGW detector. A neutron displacement-noise-free interferometer (PRD 2022) targets ~0.1 Hz, not high frequencies. An ultracold-neutron gravitational bouncing-state "reservoir" (Unnikrishnan, IJMPD 2024) is a single theoretical proposal. — [Neutron DFI, arXiv:2112.11982](https://arxiv.org/abs/2112.11982); [EBSCO record of the UCN proposal](https://www.ebsco.com/articles/science/f2701231-1511-5983-b14a-e3c60bd52fea/a-novel-wideband-neutron-earth-resonant-quantum-detector-for-gravitational-waves)

**Astrophysical and cosmological conversion (Review 2025 §5.6)**
- Figure of merit: P_{h→γ} ~ (BL)²/M_Pl² (Eq. 128). Cosmological fields give BL ≃ 10¹³ G km (B/nG)(L/Gpc), comparable to neutron stars at 10¹³ G km (B/10¹² G)(L/10 km). Plasma effective photon mass suppresses the conversion. — [Review 2025 §5.6](https://arxiv.org/abs/2501.11723)
- Recasts of radio and CMB data:
  - ARCADE 2: h_c,sto < 10⁻²⁴ (strongest cosmic fields) to < 10⁻¹⁴ (weakest) at 3–30 GHz.
  - EDGES at ~78 MHz: h_c < 10⁻¹² to < 10⁻²¹ (the extraction pairs 10⁻¹² with the strongest fields; the pairing should be checked against the review).
  - The 10-order spread reflects unknown primordial magnetic fields.
  - Sources: [Domcke & Garcia-Cely 2021, PRL 126, 021104, arXiv:2006.01161](https://arxiv.org/abs/2006.01161); [Review 2025 §5.6](https://arxiv.org/abs/2501.11723)
- Neutron stars:
  - Non-resonant conversion (Ito, Kohri, Nakayama): h_c ≲ 10⁻¹⁴–10⁻¹⁸ at 0.1–1 GHz and ≲ 10⁻¹⁶–10⁻²⁶ over 10¹³–10²⁷ Hz. The review calls this modeling "rudimentary".
  - Resonant conversion: McDonald & Ellis.
  - Population constraints: Dandoy et al.
  - Sources: [Ito, Kohri, Nakayama 2024, PRD 109, 063026, arXiv:2305.13984](https://arxiv.org/abs/2305.13984); [McDonald & Ellis 2024, PRD 110, 103003, arXiv:2406.18634](https://arxiv.org/abs/2406.18634); [Dandoy et al. 2024, JCAP 12 (2024) 023, arXiv:2402.14092](https://arxiv.org/abs/2402.14092)
- Galactic magnetic field: diffuse-photon constraints over ~10⁹–10²⁶ Hz (f ≳ 1 MHz), 1–2 orders stronger than earlier astrophysical bounds. The Galactic limits are the only astrophysical ones that translate to Ω_GW h² < 10¹⁵. LHAASO gives similar sensitivity at 100 TeV–PeV photon energies. Planetary magnetospheres (an Earth-orbit satellite and Juno at Jupiter) also give limits. — [Lella, Calore, Carenza, Mirizzi 2024, PRD 110, 083042, arXiv:2406.17853](https://arxiv.org/abs/2406.17853); [Ito, Kohri, Nakayama 2024, PTEP, arXiv:2309.14765](https://arxiv.org/abs/2309.14765); [Ramazanov et al. 2023, arXiv:2304.11222](https://arxiv.org/abs/2304.11222); [Liu, Ren, Zhang 2024, PRL 132, 131402, arXiv:2305.01832](https://arxiv.org/abs/2305.01832); [Review 2025 §3, §5.6](https://arxiv.org/abs/2501.11723)

**Contested or speculative concepts, flagged in the literature**
- Li–Baker "heterodyne enhancement of magnetic conversion" (seeding the conversion volume with a microwave beam): concerns were raised by Woods (2012) and by the 2008 JASON report (Eardley). The claims rely on non-existent technology and have no experimental results. — [Review 2025 §5.5.4](https://arxiv.org/abs/2501.11723); [JASON report](https://irp.fas.org/agency/dod/jason/gravwaves.pdf)
- Superfluid and superconducting "impedance-matched" antennas (Chiao, Anandan): h₀ ~ 10⁻³⁰ at 10 GHz is speculated. Kiefer & Weber (2005) cast doubt on the underlying assumptions. — [Review 2025 §5.7.1](https://arxiv.org/abs/2501.11723); [Kiefer & Weber, arXiv:gr-qc/0408010](https://arxiv.org/abs/gr-qc/0408010)

### Inferences
- SI form of the cavity sensitivity (Eq. 113), with μ₀ and c reinstated (derived): S_h^noise ≃ 2k_BT_sys μ₀c²/[ω³η²QB₀²L⁵]. For the ADMX benchmark (f = 0.8 GHz, η = 0.1, Q = 8×10⁴, B = 7.5 T, L = 0.51 m, T_sys = 0.6 K) this gives √S_h ≈ 1×10⁻²² Hz^{-1/2}. With t_int ~ 100 s that is h₀ ~ 10⁻²³, the same order as the review's "h₀ ~ 10⁻²²–10⁻²¹ at GHz" once O(1) factors and duty cycles are included. The other natural-unit formulas (Eqs. 117, 118, 125, 126) need the same care: each B² comes with 1/μ₀, and powers of c and ħ follow from dimensions.
- Linear detectors (cavities, interferometers, BAW, levitated sensors) beat photon counters for short transients. Photon counters win only where single-photon dark counts are tiny (optical/X-ray LSW) and the signal lasts long enough to yield ≥ 1 photon.
- For implementation, each detector reduces to a few parameters plus one closed-form S_h^noise(f) (Tables 3–5 with Eqs. 113–127), which is easy to encode as data-driven structs. The official curve data are on Zenodo (see Q5).

### Gaps
- I did not obtain the full text of QUEST, the Holometer or the MAGE runs for exact frequency-resolved limit curves. Only headline numbers are given.
- Exact numerical limits from the 2026 precision-spectroscopy paper (Budker et al.) were not in the abstract.
- The review does not resolve whether its plotted S_n values from older papers were converted from one-sided to its two-sided convention. A factor-√2 ambiguity in comparisons is possible.

## Q3. The (inverse) Gertsenshtein conversion probability in SI units, and how strain maps to photon flux

### Takeaway
In vacuum, for a magnetic field B transverse to the GW direction over a length L ≫ λ, the GW-to-photon energy conversion probability is P = 4πGB²L² in natural Heaviside-Lorentz units. Equivalently P = (κBL/2)² with κ = √(16πG). In SI this is P = 4πGB²L²/(μ₀c⁴) ≈ 8.3×10⁻³⁸ (B/1 T)²(L/1 m)². Multiplying the GW energy flux by P makes G cancel. The resulting photon power and rate are simple closed forms. Ejlli et al.'s stochastic limit formula, h_c^min ≃ 1.6×10⁻¹⁶ (1 T/B)(1 m/L)√[(Ṅ/1 Hz)(1 m²/A)(1 Hz/Δf)/ε], follows exactly from them.

### Cited Findings
- Exact vacuum solution: for a +-polarized GW along x crossing a uniform B₀ẑ region of length L, the outgoing EM energy-momentum tensor is ⟨T^{μν}⟩ = ½h₊²B₀²ω²L² (in the lightlike forward component). Dividing the Poynting flux by the GW flux h₊²ω²/(8πG) gives **P_{x>L/2} = 4πGB₀²L²** for ωL ≫ 1, with essentially no backward emission and the same result for × polarization (Eq. 140). Equivalent methods are the effective current j_eff (effective magnetization M_i = −h_ij B₀^j) and the S-matrix/Feynman route of Raffelt & Stodolsky. — [Review 2025 Appendix A](https://arxiv.org/abs/2501.11723)
- Schematic form P_{h→γ} ~ (BL)²/M_Pl² (Eq. 128). — [Review 2025 §5.6](https://arxiv.org/abs/2501.11723)
- Induced field E_h, B_h ≃ 2πfLB₀h (Eq. 119, an order-of-magnitude "≃"). Power per area ⟨S⟩ ≃ ∫df(2πfLB)²S_h (Eq. 120). Photon number N_γ = εP_sig Δt/(2πf) (Eq. 122, natural units). — [Review 2025 §5.4.1](https://arxiv.org/abs/2501.11723)
- Ejlli et al. (Lorentz-Heaviside natural units): photon energy flux from a stochastic background Φ_γ ≃ ∫ B²z²h_c²(ω) ω/4 dω (their Eq. 2, vacuum, Δz ≪ oscillation length). Requiring it to exceed the CCD background gives h_c² ≥ 4N(ω,t)/[B²L²A ε] with N = N_exp/Δω (Eq. 3). In practical units, h_c^min ≃ 1.6×10⁻¹⁶ (1 T/B)(1 m/L)√[(N_exp/1 Hz)(1 m²/A)(1 Hz/Δf)(1/ε)] (Eq. 4). Results: excluded bands 2.7–14×10¹⁴ Hz and 5–12×10¹⁸ Hz, down to h_c^min ≈ 6×10⁻²⁶ and ≈ 5×10⁻²⁸ (95% CL). — [Ejlli, Ejlli, Cruise, Pisano, Grote 2019, EPJC 79, 1032, arXiv:1908.00232](https://arxiv.org/abs/1908.00232)
- **Correction (contested):** interpreted as limits on an *isotropic* background, the LSW and helioscope numbers must be degraded by the instruments' tiny field of view, giving h_c,sto ≃ 4×10⁻²³ (optical, OSQAR II) and ≃ 8×10⁻²⁶ (X-ray, CAST). — [Liu, Ren, Zhang 2024, arXiv:2305.01832](https://arxiv.org/abs/2305.01832), as reported in [Review 2025 §5.4.1](https://arxiv.org/abs/2501.11723)
- The GW energy density is ρ_g = ⟨ḣ₊² + ḣ×²⟩/(32πG) (Review Eq. 7, natural units; multiply by c² for SI energy density). — [Review 2025 §2.2.1](https://arxiv.org/abs/2501.11723)

### Inferences
These are derived in this session and checked numerically. They hold in vacuum with B ⊥ k, ωL ≫ 1 and no waveguide or plasma effects.
- **Conversion probability (SI):** P = 4πG B²L²/(μ₀c⁴) = 8.26×10⁻³⁸ (B/1 T)²(L/1 m)². Examples: CAST (9 T, 9.26 m) gives P ≈ 5.7×10⁻³⁴; ALPS II (5.3 T, 106 m) gives P ≈ 2.6×10⁻³².
- **Monochromatic linearly polarized wave, h(t) = h₀ cos ωt:**
  - GW intensity F = c³ω²h₀²/(32πG).
  - Converted EM power P_γ = A·F·P = A B²L²ω²h₀²/(8μ₀c), with G cancelling. The equivalent induced field amplitude is E = cB_ind = (ωL/2)B h₀; the review's Eq. 119 drops this factor ½.
  - Photon rate Ṅ_γ = εP_γ/(ħω) = ε A B²L²ω h₀²/(8μ₀cħ) ≈ 2.0×10³¹ s⁻¹ × ε (h₀)²(f/1 Hz)(B/1 T)²(L/1 m)²(A/1 m²).
  - Example: ALPS II at f = 3×10¹⁴ Hz with h₀ = 10⁻²⁵ gives ~0.04 photons/s before ε and field-of-view losses.
- **Stochastic background (Maggiore/Moore h_c):** dρ_GW/dω = ωh_c²/(16πG) × c² in SI, so dΦ_γ/dω = B²L²ωh_c²/(4μ₀c). This exactly reproduces Ejlli's Eq. 2, so Ejlli's h_c uses the Maggiore/Moore normalization. The photon rate is Ṅ_γ = εAB²L²h_c²Δω/(4μ₀cħ) for flat h_c, i.e. 3.95×10³¹ s⁻¹ × ε h_c²(Δf/1 Hz)(B/1 T)²(L/1 m)²(A/1 m²). Inverting for Ṅ = 1 s⁻¹ gives exactly Ejlli's 1.6×10⁻¹⁶ prefactor.
- **With an effective photon mass (plasma or gas) or finite coherence,** the standard two-level mixing result (by analogy with axions, Raffelt & Stodolsky 1988) replaces L² with L²[sin(qL/2)/(qL/2)]², where q = |k_γ − k_g| ≈ ω_pl²/(2ω c). This is not re-verified from a primary source in this session; treat it as the standard form to check against Ejlli et al. Appendix A.
- In a Rust crate these are three short pure functions (probability, coherent photon rate, stochastic photon rate) plus the Ejlli limit inversion. They can be unit-tested against 8.26×10⁻³⁸ and 1.59×10⁻¹⁶.

### Gaps
- I did not verify the exact mixing-matrix expressions (Ejlli Eq. A.7, Δ_x,y) for media with refractive index or QED vacuum birefringence. These matter for long or gas-filled magnets and astrophysical plasma.
- The field-of-view correction factor (Liu et al. 2024) was not extracted as a formula, only its resulting numbers.

## Q4. Standard relations between h_c, S_h and Ω_GW(f), with H₀ conventions

### Takeaway
Use the Maggiore / Moore–Cole–Berry convention with one-sided PSD:
- H₀²Ω_GW(f) = (2π²/3) f³ S_h(f) = (2π²/3) f² h_c²(f), with h_c² = f S_h.
- Numerically, h_c = 1.263×10⁻¹⁸ (1 Hz/f) √(h²Ω_GW), with H₀ = 100h km s⁻¹ Mpc⁻¹ = 3.2408×10⁻¹⁸ h s⁻¹.

Warning: the 2025 Living Review switched to **two-sided** PSDs with a 4π angular normalization and eₐ·eₐ = 1 polarization tensors. Its Ω = (4π²/3H₀²) f³ S_h and h_c,sto = √(f S_h) = (H₀/2πf)√(3Ω) give an h_c that is **√2 smaller** for the same Ω. Quote Ω h² rather than h_c to avoid both this ambiguity and the H₀ dependence.

### Cited Findings
- Moore, Cole, Berry 2015, Eq. 29: H₀²Ω_GW(f) = (8πG/3c²) f² S_E(f) = (2π²/3) f³ S_h(f) = (2π²/3) f²[h_c(f)]² = (8π²/3) f⁴|h̃(f)|². Here S_h is the **one-sided** PSD and ρ_c = 3H₀²/(8πG). Their characteristic strain is defined as [h_c(f)]² = 4f²|h̃(f)|² (Eq. 17). For inspirals h_c = √(2f²/ḟ) h₀ (Eq. 35). For a monochromatic source the "characteristic strain" is the amplitude times the square root of the number of periods observed. — [Moore, Cole, Berry 2015, CQG 32, 015014, arXiv:1408.0740](https://arxiv.org/abs/1408.0740)
- Review 2025 conventions:
  - Fourier transform q(f) = ∫dt e^{−2πift}q(t) (Eq. 1).
  - Two-sided PSD ⟨q(f)q*(f′)⟩ = S_q(f)δ(f−f′) (Eq. 2); one-sided S_q^{(1)} = 2S_q (footnote 4).
  - Polarization tensors e⁺_ij = (ûᵢûⱼ − v̂ᵢv̂ⱼ)/√2, so eᵃ_ij eᵃ_ij = 2 summed over both polarizations (Eq. 5).
  - Isotropic background: ⟨h_a h_a′*⟩ = (1/4π)S_h^a δ(f−f′)δ²(Ω̂)δ_aa′ (Eq. 8), with S_h = ½ΣS_h^a.
  - Ω_GW(f) = (4π²/3H₀²)|f|³S_h(|f|) (Eq. 9); h_c,sto ≡ √(f S_h) (Eq. 12); h_c,sto = (H₀/2πf)(3Ω_GW)^{1/2} (Eq. 13).
  - Monochromatic: S_h = ½h₀²[δ(f−f_g) + δ(f+f_g)] (Eq. 15).
  - Stated equivalence: S_h = 4πS_h^{Allen–Romano} = ½S_h^{Maggiore} = S_h^{Thrane–Romano, Moore} (Eq. 11).
  - Source: [Review 2025 §2.1–2.2](https://arxiv.org/abs/2501.11723)
- The review says its convention switch from one-sided (first version) to two-sided PSDs keeps the conversion formulas "formally identical". It also explicitly warns that most h_c definitions assume matched filtering of known-phase signals in broadband interferometers and "should not be directly applied" to resonant or quadratic detectors (§2.4, footnote 10). — [Review 2025 §2.1, §2.4](https://arxiv.org/abs/2501.11723)
- H₀: the review uses H₀ = h × 100 km s⁻¹ Mpc⁻¹, assumes h = 0.7 when a number is needed, and acknowledges the early/late-Universe H₀ tension. — [Review 2025 footnote 21](https://arxiv.org/abs/2501.11723)

### Inferences
- **Internal inconsistency to flag (my reading):** Review Eq. 9 (prefactor 4π²/3) and Moore Eq. 29 (prefactor 2π²/3) agree only if S_h^{Review} = ½S_h^{Moore}, as expected for two-sided versus one-sided. Review Eq. 11, however, states S_h^{Review} = S_h^{Moore}. Either Eq. 11 counts the polarization-normalization factor differently or it is a typo. In practice h_c^{Review} = h_c^{Moore}/√2 for a given Ω, and S_h^{Review} = S_h^{Moore}/2. Encode both conventions explicitly, e.g. an enum, and test round-trips.
- Numerical constants (derived): H₀/h = 3.2408×10⁻¹⁸ s⁻¹.
  - Moore/Maggiore: h_c = 1.2634×10⁻¹⁸ (1 Hz/f)√(h²Ω).
  - Review 2025: h_c,sto = 8.934×10⁻¹⁹ (1 Hz/f)√(h²Ω).
  - Inverse (Moore): h²Ω = 6.27×10³⁵ (f/Hz)² h_c² = 6.27×10³⁵ (f/Hz)³ (S_h^{one-sided}/Hz⁻¹).
- Exact chain for a Rust implementation in the Moore convention: S_h^{(1)}(f) = h_c²/f; Ω_GW = 2π²f³S_h^{(1)}/(3H₀²). Monochromatic: S_h^{(1)} = ½h₀² × δ(f−f_g) for h = h₀cos(2πf_gt), and one-polarization intensity F = c³(2πf)²h₀²/(32πG). Both conventions give the same h_c to Ω mapping in "h²Ω" units only if the prefactor matches the chosen h_c definition.
- Converting a detector's noise to an "Ω sensitivity" is meaningful only with an integration time and bandwidth via the power-law-integrated (PLS) recipe (Review Eqs. 23, 26, 35). The naive Ω_n(f) = 2π²f³S_n/(3H₀²) is the per-frequency noise-equivalent energy density, not a detection threshold.

### Gaps
- I did not open Maggiore's Phys. Rep. (gr-qc/9909001) in this session. The "1.263×10⁻¹⁸" prefactor is derived from Moore Eq. 29, which matches it.
- The Review's Eq. 11 inconsistency could not be resolved without contacting the authors or checking the typeset journal version.

## Q5. Numbers for plotting the HFGW landscape: published curves and anchor points

### Takeaway
The authoritative machine-readable source is the HFGWPlotter dataset by Muia, Ringwald and Tamarit (Zenodo, June 2025), which backs the Review 2025 figures: noise-equivalent strain curves (Figs. 1–2) and stochastic Ω sensitivities and signals (Fig. 4). The review deliberately stopped drawing a single "h_c vs f" plot with sources and detectors overlaid, because it hides integration-time and bandwidth variables. Use noise-equivalent √S_n(f) for detectors and h(f), h₀ or Ω h² for sources.

### Cited Findings
- HFGWPlotter datasets: "High frequency gravitational wave plotter: Noise-equivalent strain", DOI 10.5281/zenodo.15720342, and "...: Stochastic signals", DOI 10.5281/zenodo.15720443. — [Review 2025 refs. 19–20](https://arxiv.org/abs/2501.11723); [Zenodo noise curves](https://doi.org/10.5281/zenodo.15720342); [Zenodo stochastic](https://doi.org/10.5281/zenodo.15720443)
- Fig. 4 of the review uses PLS curves for broadband detectors (1 yr, SNR_thr = 10, Ω = Ω₀(f/f*)^α templates; Eq. 35). Resonant detectors are shown at a fixed frequency for 1 yr. For photon counters (OSQAR, ALPS, CAST, IAXO) the curve is set by "≥ 1 photon in 1 yr". The horizontal BBN line is Ω h² ≈ 10⁻⁶. Fig. 3 shows PBH distance reach for broadband detectors (SNR 10, equal mass, ISCO cutoff) against the 1 event/yr distance for f_PBH = 1 and 0.1%. — [Review 2025 §3](https://arxiv.org/abs/2501.11723)
- Footnote 10: the first review's combined h_c plot "contains many hidden variables ... which may lead to misleading conclusions". — [Review 2025 §3](https://arxiv.org/abs/2501.11723)

**Detector anchor points.** Each entry gives: quantity; frequency; status; sources.
- **Published limits or measured sensitivities**
  - Holometer: √S_n ≃ 10⁻²¹ Hz^{-1/2} (cross-correlated); h_c < 2.5×10⁻¹⁸ (1 MHz) to 2.4×10⁻¹⁹ (13 MHz); 1–13 MHz; published. — [arXiv:1611.05560](https://arxiv.org/abs/1611.05560), [Review 2025](https://arxiv.org/abs/2501.11723), [Ejlli 2019](https://arxiv.org/abs/1908.00232)
  - QUEST: √S_n ≈ 3×10⁻²⁰ Hz^{-1/2} at 40 MHz; limits over 13–80 MHz; published (PRL 2025). — [arXiv:2410.09175](https://arxiv.org/abs/2410.09175)
  - 100 MHz interferometers: √S_n ~ 10⁻¹⁶, h_c ≲ 7×10⁻¹⁴ (Review) or 1.4×10⁻¹² (Ejlli); 100 MHz; published. — [arXiv:0803.4094](https://arxiv.org/abs/0803.4094)
  - Cruise–Ingley polarization ring: √S_n ≲ 10⁻¹⁴, h_c ≤ 1.4×10⁻¹⁰; 100 MHz; prototype. — [Review 2025](https://arxiv.org/abs/2501.11723), [Ejlli 2019](https://arxiv.org/abs/1908.00232)
  - Magnon reinterpretation: √S_n ≲ 7.6×10⁻²² at 14 GHz and 1.2×10⁻²⁰ at 8.2 GHz (h_c ≤ 9.1×10⁻¹⁷ and 1.1×10⁻¹⁵); recast. — [arXiv:1903.04843](https://arxiv.org/abs/1903.04843)
  - CAPP-12T reanalysis: h₀ ≈ 3.9×10⁻²¹ (90% CL, best sky regions); 2 MHz span at 5.311 GHz; published reanalysis (late 2025). — [arXiv:2511.17817](https://arxiv.org/abs/2511.17817)
  - ABRACADABRA-10cm: √S_n ≃ 8.4×10⁻⁹ Hz^{-1/2} (f/1 MHz)⁻²; transient (PBH) sensitivity h ~ 4×10⁻⁵ near 0.2 MHz; 10 kHz–5 MHz; first HFGW-dedicated haloscope search, Dec 2023 data. — [arXiv:2505.02821](https://arxiv.org/abs/2505.02821)
  - OSQAR / ALPS (optical): h_c^min ≈ 6×10⁻²⁶ (Ejlli), or ≃ 4×10⁻²³ after field-of-view correction; 2.7–14×10¹⁴ Hz; recast. — [arXiv:1908.00232](https://arxiv.org/abs/1908.00232), [arXiv:2305.01832](https://arxiv.org/abs/2305.01832)
  - CAST (X-ray): h_c^min ≈ 5×10⁻²⁸ (Ejlli), or ≃ 8×10⁻²⁶ after field-of-view correction; 5–12×10¹⁸ Hz; recast. — same sources
  - MiniGRAIL / Schenberg: √S_n ≃ 1.5×10⁻²⁰ at 2.94 kHz / 1.1×10⁻¹⁹ at 3.2 kHz; published. — [Review 2025](https://arxiv.org/abs/2501.11723)
  - MAGE (BAW): √S_n down to 6.6×10⁻²¹ in narrow MHz bands; operating. — [arXiv:2307.00715](https://arxiv.org/abs/2307.00715)
- **Projections**
  - LSD 1 m prototype: √S_n ~ 10⁻¹⁹ (10 kHz) and ~10⁻²¹ (100 kHz); 10 m: 10⁻²⁰ and 10⁻²²; 10–100 kHz; under construction. — [arXiv:2010.13157](https://arxiv.org/abs/2010.13157)
  - MAGO 2.0: √S_n ~ 10⁻²¹ (resonant), ≲ 10⁻¹⁸ (broadband); ~10 kHz–GHz; R&D. — [arXiv:2303.01518](https://arxiv.org/abs/2303.01518)
  - BAW estimate: √S_n ≃ 2×10⁻²²; 1–700 MHz; R&D. — [arXiv:1410.2334](https://arxiv.org/abs/1410.2334)
  - Magnetic Weber bar (MRI magnet): √S_n ~ 10⁻²⁰ broadband, ~10⁻²² at kHz; kHz–10 MHz; proposal. — [arXiv:2408.01483](https://arxiv.org/abs/2408.01483)
  - Levitated superconducting sphere: √S_n ~ 10⁻¹⁹; 10 kHz–1 MHz; proposal. — [arXiv:2408.01583](https://arxiv.org/abs/2408.01583)
  - LC haloscopes (recast or projected, coherent signals): h₀ ~ 10⁻⁹ (2 MHz) to 10⁻²¹ (100 MHz, DMRadio-m³). — [arXiv:2202.00695](https://arxiv.org/abs/2202.00695)
  - IAXO in the GHz band: h_c ≲ 10⁻²² (heterodyne), ≲ 10⁻²⁵ (single-photon); ~40 GHz; proposal. — [arXiv:2011.04731](https://arxiv.org/abs/2011.04731)
  - MADMAX resonant: √S_n ~ 10⁻²² × (scalings); ~10 GHz; R&D. — [arXiv:2409.06462](https://arxiv.org/abs/2409.06462)
  - Quantum cyclotron: h₀ ~ 3.8×10⁻²⁰, down to 6.9×10⁻²³; 20–200 GHz; proposal. — [arXiv:2309.02992](https://arxiv.org/abs/2309.02992)
  - Rydberg atoms: h₀ ≲ few×10⁻²⁰; ~30 GHz; proposal. — [arXiv:2311.03890](https://arxiv.org/abs/2311.03890)
  - NEMO: √S_n ≃ 10⁻²⁴ Hz^{-1/2}; 1–2.5 kHz; R&D. — [Review 2025](https://arxiv.org/abs/2501.11723)
- **Astrophysical (isotropic SGWB, model-dependent)**
  - ARCADE 2: h_c < 10⁻²⁴ to 10⁻¹⁴; 3–30 GHz. EDGES: h_c < 10⁻¹² to 10⁻²¹; 78 MHz. — [arXiv:2006.01161](https://arxiv.org/abs/2006.01161)
  - Neutron stars (non-resonant): h_c ≲ 10⁻¹⁴–10⁻¹⁸ (0.1–1 GHz); ≲ 10⁻¹⁶–10⁻²⁶ (10¹³–10²⁷ Hz). — [arXiv:2305.13984](https://arxiv.org/abs/2305.13984)

**Source anchor points (upper envelopes, Review 2025 Table 2 and Fig. 4)**
- Ω ≲ 10⁻¹⁶ (inflation, vacuum); ≲ 10⁻¹⁰ (extra species, 10⁵–10⁸ Hz; broken spatial reparametrization); ≲ 10⁻⁸ (second-order induced; phase transitions ≲ 10⁹ Hz; GUT-PBH evaporation).
- Preheating and oscillons: ≲ 10⁻¹⁰ at 10⁶–10⁹ Hz.
- CGMB: ≲ 10⁻⁶ at the 10–100 GHz peak, with h_c ≲ 10⁻³¹ (MHz/f).
- Table 2 strain envelopes: inflation h_c ≲ 10⁻³² (MHz/f); extra species and preheating ≲ 10⁻²⁹ (MHz/f); phase transitions ≲ 10⁻²⁸ (MHz/f). These use the review's h_c convention.
- — [Review 2025 Table 2](https://arxiv.org/abs/2501.11723)

### Inferences
Derived with the Moore convention h²Ω = 6.27×10³⁵ f²h_c².
- Holometer h_c < 2.4×10⁻¹⁹ at 13 MHz corresponds to Ω h² ~ 6×10¹², and 2.5×10⁻¹⁸ at 1 MHz to ~4×10¹².
- QUEST √S_n = 3×10⁻²⁰ at 40 MHz corresponds to Ω h² ~ 4×10¹⁹ (one-sided) or ~7×10¹⁹ (review two-sided).
- The 100 MHz interferometer gives Ω h² ~ 3×10²⁵ using the review's h_c = 7×10⁻¹⁴, or ~10²⁸ using Ejlli's 1.4×10⁻¹².
- These are 18–31 orders of magnitude above the BBN line, which is the main message for a landscape plot.
- A good landscape figure has two panels: (i) √S_n(f) for detectors with resonant curves dashed, against coherent-source h₀ and transient h(f) at reference distances (1 kpc for superradiance and PBH mergers, the d_yr distance from Fig. 7); (ii) Ω h²(f) with PLS curves and the BBN band. Drawing from the HFGWPlotter Zenodo files avoids digitizing figures.

### Gaps
- I did not download or inspect the Zenodo files, so their format and column conventions (one- or two-sided) are unverified.
- Exact frequency-resolved limit curves for QUEST, MAGE, CAPP-12T and ABRA-10cm were not extracted; only headline points are given.

## Q6. The experimental frontier in 2024–2026: new results, the BAW "rare events", and new proposals

### Takeaway
2024–2026 brought the first *dedicated* HFGW data analyses with laboratory hardware:
- QUEST's 13–80 MHz broadband limits (PRL 2025);
- ABRACADABRA-10cm's 10 kHz–5 MHz search (2025);
- the CAPP-12T reanalysis at 5.3 GHz (2025);
- precision-spectroscopy constraints at 100 kHz–100 MHz (July 2026 preprint);
- the first cryogenic characterisation of the MAGO SRF cavity (2025–2026).

There are no levitated-sensor GW limits yet; the LSD is still under construction in 2026. The 2021 BAW "rare events" are not regarded as GWs, following Lasky & Thrane (2021) and the 2025 review. A 2026 quantum-limit analysis argues that only kHz–MHz stochastic backgrounds are detectable even in principle.

### Cited Findings
- **BAW rare events.**
  - 153 days over two runs. Event 1 was on 12 May 2019 in the 5.506 MHz mode (3B) only; the 8.392 MHz mode (5C) saw nothing. Event 2 was on 27 Nov 2019, seen in both the 5.506 MHz (3B) and 4.993 MHz (3C) modes. The UWA press release described the events as possibly HFGWs (from PBHs or dark-matter clouds) while the origin stayed open. — [Goryachev et al. 2021, PRL 127, 071102, arXiv:2102.05859](https://arxiv.org/abs/2102.05859); [UWA press release](https://www.uwa.edu.au/news/Article/2021/August/World-first-detector-reports-rare-events)
  - **Disputed by Lasky & Thrane (PRD 104, 103017, 2021).** A stochastic interpretation implies Ω_GW ≈ 10⁸. Any high-frequency source loud enough would produce GW memory seen by LIGO/Virgo as a burst with SNR ≳ 10⁶, which was not observed. — [Lasky & Thrane, arXiv:2110.13319](https://arxiv.org/abs/2110.13319)
  - Review 2025: "given their strength they are not considered to be viable GW candidates". MAGE was built partly to identify the origin. — [Review 2025 §5.2.2](https://arxiv.org/abs/2501.11723); [MAGE, arXiv:2307.00715](https://arxiv.org/abs/2307.00715)
- **QUEST (PRL 2025):** described above in Q2. — [arXiv:2410.09175](https://arxiv.org/abs/2410.09175)
- **ABRACADABRA-10cm (arXiv:2505.02821, v2 Sept 2026):** a semicircular figure-8 pickup coupling to h₊ and h×. Six days of data (4–9 Dec 2023, T = 518,459 s). Covers 10 kHz–5 MHz, "the first limit at frequencies above those probed by LIGO O4a". √S_n ≃ 8.4×10⁻⁹ Hz^{-1/2}(f/1 MHz)⁻², matching theory. A PBH transient search over 24,975 segments excludes sources only within 46 km, i.e. h ~ 10⁻⁴. — [Pappas et al., arXiv:2505.02821](https://arxiv.org/abs/2505.02821)
- **CAPP-12T reanalysis (Nov 2025):** 2 MHz span at 5.311 GHz with no rescan candidates. h₀ ≈ 3.9×10⁻²¹ (90% CL) in the best sky regions. In a superradiance benchmark it excludes ~1.22×10⁻⁶ M_⊙ black holes within ~O(10⁻²) AU of Earth. — [Kim et al., arXiv:2511.17817](https://arxiv.org/abs/2511.17817)
- **Precision spectroscopy (July 2026):** GW-induced frequency modulation of light in laser cavities, analyzed with optical spectroscopy, constrains "much of" 100 kHz–100 MHz "for the first time". The authors project 8 orders of magnitude improvement and extension to ≥ 1 GHz. — [Budker, Domcke, Kopp, Tretiak, arXiv:2607.12617](https://arxiv.org/abs/2607.12617)
- **MAGO / SRF:** CQG 2025 first characterisation, then a Jan 2026 search-design paper (10 kHz–100 MHz, 1.8 K), then July 2026 cryogenic tests (11 kHz splitting, mechanical Q well below the assumed values, an unfavorable surprise for projections). — [arXiv:2411.18346](https://arxiv.org/abs/2411.18346); [arXiv:2601.18719](https://arxiv.org/abs/2601.18719); [arXiv:2607.02476](https://arxiv.org/abs/2607.02476)
- **New concepts 2023–2025:**
  - magnetic Weber bars, i.e. axion magnets as resonant masses — [arXiv:2408.01483](https://arxiv.org/abs/2408.01483)
  - superconducting levitated detector — [arXiv:2408.01583](https://arxiv.org/abs/2408.01583)
  - dielectric haloscopes (MADMAX/DALI) — [arXiv:2409.06462](https://arxiv.org/abs/2409.06462)
  - optical-cavity polarimetry at ALPS II — [arXiv:2501.08382](https://arxiv.org/abs/2501.08382)
  - quantum cyclotron — [arXiv:2309.02992](https://arxiv.org/abs/2309.02992)
  - Rydberg atoms — [arXiv:2311.03890](https://arxiv.org/abs/2311.03890)
  - optical frequency modulation — [arXiv:2304.10579](https://arxiv.org/abs/2304.10579)
  - cubic SRF cavity for polarization and direction — [Navarro et al. 2024, arXiv:2312.02270](https://arxiv.org/abs/2312.02270)
  - astrophysical conversion: NS (resonant/non-resonant/population), Galactic, planetary — Q2 above.
- **Fundamental quantum limit (PRD 113, 044020, Feb 2026):** a universal quantum bound across detection schemes. Stochastic backgrounds at kHz–MHz are "in principle observable", while higher-frequency backgrounds lie below the quantum limit. — [Guo, Miao, Wang, Yang, Zhou, arXiv:2501.18146](https://arxiv.org/abs/2501.18146)
- **Community:** the UHF-GW Initiative (http://www.ctc.cam.ac.uk/activities/UHF-GW.php) grew out of workshops at ICTP Trieste 2019, online 2021 and CERN 2023, whose summary is the Living Review. The review states "None of the detection concepts listed in this report currently reach the sensitivity needed to probe realistic sources". It nonetheless sees "clear development paths" toward √S_n ≃ 10⁻³⁸ Hz^{-1/2} via optical-frequency magnetic conversion; this is a speculative extrapolation. — [Review 2025 §1, §6](https://arxiv.org/abs/2501.11723); [Living Rev. Relativ. record (DESY)](https://bib-pubdb1.desy.de/record/646132)

### Inferences
- The realistic near-term discovery targets are nearby coherent sources: superradiance around light PBHs at ≲ AU to kpc distances, and very close PBH mergers. They are probed by linear, resonant, long-integration searches (cavity and haloscope reanalyses, LSD, BAW arrays, MAGO). Cosmological backgrounds above ~MHz are below both the BBN line and, per Guo et al. 2026, the quantum limit.
- The MAGO mechanical-Q result and the ABRA-10cm "matches theory but h ~ 10⁻⁴ transient reach" result show that 2024–2026 hardware is still validating noise models, not approaching source-level sensitivity.

### Gaps
- I found no published GW limits from levitated sensors, BAUSCIA, GravNet or DMRadio as of Oct 2026. Absence of evidence in my searches, not proof.
- The Goryachev et al. full text (event energies and alternative explanations considered, such as cosmic rays or internal stress release) was not retrieved. The abstract lists none of these.
- Numerical limits of the July 2026 spectroscopy paper and the full CAPP-12T sky map were not extracted.
- I did not check for MAGE or ORGAN GW results published after 2023, or for HAYSTAC/ADMX GW reanalyses beyond CAPP-12T.
