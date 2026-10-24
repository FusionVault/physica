# Atom-interferometer and quantum-sensor gravitational-wave detectors, and their links to gravimetry (2023–2026)

Scope: notes for extending the Rust crate `satgrav` (which already models cold-atom accelerometers and
gradiometers with phase `2·k_eff·a·T²` and shot noise `1/(C√N)`) with GW strain sensitivity and the
noise sources shared with gravimetry. "Demonstrated" means measured in hardware; everything else is a
design target or projection. Numbers marked *(computed here)* are my own evaluations of the cited
formulas, not values from the sources.

## 1. Status and design of atom-interferometer GW detectors (AION, MAGIS-100, ZAIGA, MIGA, ELGAR, AEDGE/AEDGE+, space, TVLBAI roadmap)

### Takeaway
As of late 2026 no atom-interferometer (AI) GW detector is operating. The two 100 m-class instruments are
MAGIS-100 (Fermilab; its laser lab was completed in January 2026, commissioning expected by 2028) and the
proposed AICE (a 100 m AI in CERN's PX46 shaft, letter of intent in 2025). AION-10 (Oxford, 10 m) has a
2025 technical design report and is aimed at dark matter, not GWs. The TVLBAI roadmap targets km-scale
detectors from 2035. Demonstrated large momentum transfer (LMT) is about 400ħk, while GW-grade designs
assume n = 5 000–40 000. That gap of one to two orders of magnitude in LMT, together with the phase-noise
goals, separates today's hardware from the targets.

### Cited Findings

**Design parameter tables (all are targets or projections)**

| Scenario | L [m] | T [s] | δφ_noise [rad/√Hz] | LMT n (or n_p^max) | Other | Source |
|---|---|---|---|---|---|---|
| AION-10 initial | 10 | 1.4 | 1e-3 | 100 | Sr, single-photon clock | [Badurina et al. 2020, JCAP 05 011, arXiv:1911.11755, Table 1](https://arxiv.org/abs/1911.11755) |
| AION-10 goal | 10 | 1.4 | 1e-4 | 1 000 | | same |
| AION-100 initial | 100 | 1.4 | 1e-4 | 1 000 | | same |
| AION-100 goal | 100 | 1.4 | 1e-5 | 40 000 | | same |
| AION-km | 2 000 | 5 | 0.3e-5 | 40 000 | | same |
| AICE (100 m, CERN PX46) | 100 | 1.4 | 1e-5 | 5 494 | Δz = 50 m between 2 sources, N = 1e8 atoms, Δt = 1.5 s, T_int = 1e8 s | [Baynham et al. 2025, AICE LoI, arXiv:2509.11867, Table 1](https://arxiv.org/abs/2509.11867) |
| Earth-km (benchmark) | 2 000 | T_int = 5 | 0.3e-5 | n_p^max = 40 000 | | [El-Neaj et al. (AEDGE) 2020, EPJ Quantum Technol. 7:6, arXiv:1908.00802, Table 1](https://arxiv.org/abs/1908.00802) |
| AEDGE (space, 2 satellites in MEO) | 4.4e7 | T_int = 300 | 1e-5 | n_p^max = 1 000 | n_p = 2Q(2n−1)+1, Δx = ħk_eff T/m < 90 cm (clouds inside the spacecraft) | same |
| Graham 2016 space example | 4.4e7 | T_max = 300 | δφ_a = 1e-5 | n_max = 1e3 | 50 cm telescope, ~10 pK atoms | [Graham, Hogan, Kasevich, Rajendran 2016, PRD 94, 104022, arXiv:1606.01860](https://arxiv.org/abs/1606.01860) |

- In AION Table 1, "T" is defined as giving an AI duration of 2T, and T = 1.4 s corresponds to a 10 m AI in launch mode. AION projections assume 5 years of observation split across 10 log-spaced resonance frequencies, and use power-law-integrated (PI) curves at SNR = 10 — [Badurina et al. 2020](https://arxiv.org/abs/1911.11755)
- AEDGE+ is a variant with atom clouds about 100 m in size outside the spacecraft. That allows larger LMT and lower-frequency reach, but its technical feasibility was not assessed — [Badurina et al. 2021, Phil. Trans. R. Soc. A 380:20210060, arXiv:2108.02468](https://arxiv.org/abs/2108.02468); [AEDGE paper](https://arxiv.org/abs/1908.00802)
- AEDGE was submitted to ESA's Voyage 2050 call. An ESA senior committee received it favourably but called for a preparatory technology programme. I found no 2024–2026 ESA decision — [AEDGE summary page](https://www.emergentmind.com/topics/aedge) (secondary); [AEDGE white paper](https://www.cosmos.esa.int/documents/1866264/3219248/BuchmuellerO_White-Paper-AEDGE-submitted-02082019.pdf/678dff63-4773-3590-65a8-0494f5b09f8a?t=1565184625388)
- Badurina et al. 2021 state that AION-10 "is not expected to detect GWs". Their focus is AION-100, AION-km and AEDGE. AION-100 sites under consideration are Boulby, Daresbury (UK) and CERN PX46. The 2021 study uses the full sky-averaged antenna response instead of the optimal response used earlier — [arXiv:2108.02468](https://arxiv.org/abs/2108.02468)

**Project status, 2025–2026**
- **MAGIS-100 (Fermilab, NuMI shaft, 100 m)**: a 100 m magnetically shielded vacuum tower, three Sr atom sources, and a laser room at the top. Laser-lab construction began in 2023 and was complete by January 2026. Equipment is being moved in, and commissioning is expected by 2028. Science goals are ultralight dark matter, quantum mechanics at macroscopic scales, and GW-technology demonstration; one report says it is not expected to be sensitive enough to detect known GW sources — [Fermilab news tag MAGIS-100](https://news.fnal.gov/tag/magis-100/); [Jan 2026 repost of Fermilab release](https://sciencesprings.wordpress.com/2026/01/10/from-the-does-fermi-national-accelerator-laboratory-fermilab-completes-laser-lab-construction-for-worlds-largest-vertical-atom-interferometer/); [Interesting Engineering](https://interestingengineering.com/space/laser-lab-for-worlds-largest-vertical-interferometer). Design paper: [Abe et al. 2021, Quantum Sci. Technol. 6, 044003, arXiv:2104.02835](https://arxiv.org/abs/2104.02835) (its PDF exceeded the fetch limit, so its parameter table was not extracted).
- **AION-10 (Oxford)**: a technical design report was posted in August 2025. It describes a 10 m vertical tower with two ultracold-Sr sources in UHV, active vibration isolation and magnetic shielding. Analysis shows critical optics within spec 97% of the time; this is an analysis, not measured performance — [Bongs et al. 2025, AION-10 TDR, arXiv:2508.03491](https://arxiv.org/abs/2508.03491)
- **AION prototype (demonstrated)**: a tabletop single-photon long-baseline-style AI on the ⁸⁷Sr clock transition operates at the standard quantum limit, even with added laser phase noise. This shows the laser-noise rejection of the single-photon gradiometer scheme — [AION-REPORT/2025-02, arXiv:2504.09158](https://arxiv.org/abs/2504.09158) (seen via search snippet only)
- **AICE (CERN PX46)**: PX46 is about 143 m from the SX4 surface building to LHC level, with a 10.1 m internal diameter. The LoI is CERN-LHCC-2025-008. Site preparation would run during LS3 and the AI would operate during Run 4. Experiment cost is estimated at 30–50 MCHF; the technical services prefer installation work in H2 2028. Peak sensitivity band is 50 mHz–10 Hz. Magnetic noise must be ≤ 100 nT/√Hz, with < 100 pT/√Hz desirable. Seismic noise measured at PX46 lies between the NHNM and NLNM. The design has up to 5–10 atom sources along the tube — [arXiv:2509.11867](https://arxiv.org/abs/2509.11867); implementation study [Arduini et al., CERN-PBC-Report-2025-004, arXiv:2508.09694](https://arxiv.org/abs/2508.09694) (cited, not fetched)
- **MIGA (France, LSBB Rustrel, about 300 m underground)**: a horizontal two-arm detector with a 150 m maximum gradiometer length, under construction in 2025 — [TVLBAI ESPP input 2025, arXiv:2503.21366](https://arxiv.org/abs/2503.21366); [Junca et al. 2019, PRD 99 104026, arXiv:1902.05337](https://arxiv.org/abs/1902.05337); design [Canuel et al. 2018, Sci. Rep. 8, 14064](https://arxiv.org/abs/1703.02490)
- **ZAIGA (China, near Wuhan)**: stage 1 (a 240 m vertical shaft and a horizontal gallery longer than 1 km) is funded and under construction, with completion scheduled for 2027. Phase 2 (2027–2035) would add a triangle of 1 km galleries for GW and dark-matter searches — [TVLBAI 2023 workshop summary, arXiv:2310.08183](https://arxiv.org/abs/2310.08183) (via search summary); design [Zhan et al. 2020, IJMPD 29, 1940005, arXiv:1903.09288](https://arxiv.org/abs/1903.09288) (not fetched)
- **ELGAR (Europe)**: a proposal (not funded) for a multi-km, two-arm horizontal detector. It is reported as 32 km arms and a peak strain of 3.3e-22/√Hz at 1.7 Hz, from a secondary search summary not checked against the paper — [TVLBAI ESPP input](https://arxiv.org/abs/2503.21366); [Canuel et al. 2020, CQG 37, 225017, arXiv:1911.03701](https://arxiv.org/abs/1911.03701)
- **10 m vertical prototypes**: Stanford, and VLBAI in Hannover — [arXiv:2503.21366](https://arxiv.org/abs/2503.21366)

**TVLBAI roadmap**
- 1st workshop (CERN, 2023): the goal is a roadmap for one or more km-scale detectors operating by the mid-2030s — [Abend et al. 2024, AVS Quantum Sci. 6, 024701, DOI 10.1116/5.0185291, arXiv:2310.08183](https://arxiv.org/abs/2310.08183). 2nd workshop: [Abdalla et al. 2025, EPJ Quantum Technol. 12, 42](https://arxiv.org/abs/2509.11867) (citation from the AICE reference list). 3rd workshop: Leibniz University Hannover ([indico 1522217](https://indico.cern.ch/event/1522217/)).
- The proto-collaboration MoU has been signed by 55 institutions in 20 countries, plus 3 observers. Each workshop had more than 200 participants. AICE is endorsed as a "priority demonstrator" — [arXiv:2509.11867, App. A](https://arxiv.org/abs/2509.11867)
- Timeline: 2025–2029 for detailed design, prototyping and small-facility tests (MAGIS-100 and MIGA under construction, a possible 100 m AI in PX46 during LS3); 2029–2034 for construction and commissioning of one or more full-scale interferometers; 2035 onwards for full operation — [TVLBAI ESPP input, arXiv:2503.21366](https://arxiv.org/abs/2503.21366)
- Cost estimates: an initial 100 m ULDM-capable detector at about 30 MCHF, plus about 20 MCHF to reach ultimate sensitivity (25% contingency included); km-scale about 10× more. Infrastructure: 2–10 atom sources at 100 m and up to 100 at 1 km; ≤ 100 kW (1 MW) for 100 m (1 km) atom sources — [arXiv:2503.21366](https://arxiv.org/abs/2503.21366)

**Demonstrated LMT and atom-optics performance (vs targets)**
- Floquet atom optics on the Sr ¹S₀–³P₁ (689 nm) line: momentum separation above 400ħk, pulse efficiency above 99.4% — [Wilkason et al. 2022, PRL 129, 183202, arXiv:2205.06965](https://arxiv.org/abs/2205.06965)
- First LMT clock-transition (single-photon) AI: up to 141ħk — [Rudolph et al. 2020, PRL, arXiv:1910.05459](https://arxiv.org/abs/1910.05459)
- The AION collaboration "has demonstrated LMTs using several hundred laser pulses", in a paper still in preparation as of September 2025 — [AICE LoI footnote 2](https://arxiv.org/abs/2509.11867)
- In 2016, n_max = 1e3 pulses was "an order of magnitude above the current state-of-the-art for a discrete pulse LMT sequence" — [Graham et al. 2016](https://arxiv.org/abs/1606.01860)

### Inferences
- Peak broadband strain ASD, h ≈ δφ/(2·n·k·L) with k = 2π/698 nm *(computed here, optimal orientation, no sky averaging)*: AION-10 initial 5.6e-14; AION-100 initial 5.6e-17; AICE (L = Δz = 50 m) 2.0e-18; AION-100 goal 1.4e-19; AION-km 2.1e-21 (all at f = 1/(2T)); AEDGE with n ≈ 250 about 5e-23/√Hz near 3 mHz. These follow the published order of magnitude and are useful as unit tests for a `satgrav` GW module.
- AICE's δφ = 1e-5 rad/√Hz with N = 1e8 atoms per 1.5 s cycle is about 12× below the single-cloud SQL (√1.5/√1e8 ≈ 1.2e-4 rad/√Hz) *(computed here)*. It therefore implies about 20 dB of spin squeezing, concurrent interleaved clouds, or higher flux. The LoI's text mentions squeezing and more atoms as routes. A crate should store δφ as an input, not derive it from N alone, or it should expose a squeezing factor.
- The AI programme has the same technology stack as satgrav's CARIOQA/GRICE-type sensors (cold atoms, LMT, long T, magnetic shielding), but uses Sr single-photon clock transitions rather than Rb Raman. So `K_EFF_RB87` needs a sibling constant for Sr-87 698 nm with k = ω_A/c (single photon, not 4π/λ).

### Gaps
- MAGIS-100 parameter table (n_LMT, atom number, δφ target, strain curve): the PDF exceeded the fetch limit. A 2025–2026 MAGIS-100 strain curve was not found.
- No 2025–2026 primary update was found for ZAIGA, MIGA or ELGAR. ELGAR numbers are unverified.
- No ESA decision or roadmap for AEDGE or other space AI GW missions dated 2024–2026 was found. I found no reliable source on newer Chinese space AI proposals (e.g. AIGSO/AMIGO are mentioned only by name in [Shimoda et al.](https://arxiv.org/abs/1812.01835)).
- TVLBAI 2nd-workshop summary contents were not read (citation only).

## 2. Strain-sensitivity formula for a two-interferometer gradiometer GW detector

### Takeaway
For a single-photon (clock-transition) gradiometer of baseline L, the GW phase response is Graham et al.
2016 Eq. 9. In broadband mode at low frequency it reduces to Δφ ≈ 2·n·k·h·L·sin²(ωT/2), and on resonance
with Q diamonds the peak is 2·Q·n·k·h·L. The strain ASD is h_n(f) = δφ_noise / |Δφ/h|. This is the same
Mach–Zehnder acceleration transfer function satgrav already uses, driven by a GW tidal acceleration
Δa = ½ω²hL. So the GW mode can be added as a thin layer over the existing gradiometer model.

### Cited Findings

**Response function (Graham et al. 2016, PRD 94, 104022, arXiv:1606.01860)** — [source](https://arxiv.org/abs/1606.01860)
- GW: h(t) = h cos(ω(t − z/c) + φ₀). Light-arrival-time difference across the gradiometer (Eq. 7): Δτ(t) = L/c + (h/2ω)[sin(ω(t + L/c) + φ₀) − sin(ωt + φ₀)].
- Response amplitude (Eq. 9), with ΔΦ_grad(t₀) = Δφ·cos(ωt₀ + φ₀):

  `Δφ = k_eff·h·L · [sin(ωQT)/cos(ωT/2)] · sinc(ωnL/(2c)) · sin(ωT/2 − ω(n−1)L/(2c))`

  where k_eff ≡ n·ω_A/c for an n-pulse LMT sequence on an optical transition of angular frequency ω_A, T is the pulse spacing (each diamond lasts 2T), and Q is the number of diamonds.
- Resonance at ω_r ≡ π/T, bandwidth ~ω_r/Q. On resonance (Eq. 10): `Δφ_res = 2·Q·k_eff·h·L·sinc(ω_r n L/(2c))·cos(ω_r (n−1) L/(2c))`, which is ≈ 2Qk_eff hL for ω_r ≪ c/(nL).
- Constraints: total pulses n_p = 2Q(2n−1)+1 ≤ n_max; total time 2TQ ≤ T_max; LMT duration 2(n−1)L/c < T. Swept-resonance mode improves strain sensitivity to long-lived sources by about √Q.
- Stochastic sensitivity for two cross-correlated detectors (Eq. 11): `Ω_GW(f) = [π c² f³ / (ρ_c G |γ(x₁,x₂,f)|)] · √(2/(τ_int Δf)) · 1.645 · h_n²(f)`, with γ ≤ 8π/5 and Δf ≈ f_r/Q. So Ω_GW ∝ Q^(−3/2) at fixed f_r.
- Technical-noise requirements set equal to atom shot noise δφ_a = 1e-5 rad/√Hz (space design, L = 4.4e7 m):
  - Photon shot noise of the laser link: telescope diameter `d = 50 cm·(L/4.4e7 m)^(1/2)·(1 W/P_t)^(1/5)·((f_R/n_p)/(2 Hz/200))^(1/5)·(1e-5/δφ_a)^(2/5)` (Eq. 12, exponents as best reconstructed from the extracted text), with d_min ∼ L^(2/5) n^(1/5).
  - Wavefront aberration: δλ = (λ/300)·(200/(nQ))·(Λ/1 cm)·(δφ_a/1e-5)·(1 nm/√Hz / δx) (Eq. 13).
  - Laser frequency noise: δω = 2π·10 Hz/√Hz·(200/(nQ))·(1 ns/t_d)·(δφ_a/1e-5) (Eq. 14).
  - Platform acceleration noise couples through the relative velocity Δv between the two interferometers.
- The same Eq. 2.6 form (written 2QnkhL with k = ω_A/c) is used by AION, which optimises n and Q per frequency under the n_p and T_max constraints — [Badurina et al. 2020](https://arxiv.org/abs/1911.11755)

**Gradiometer strain PSD including GGN and detection noise (MIGA; Junca et al. 2019, PRD 99, 104026, arXiv:1902.05337)** — [source](https://arxiv.org/abs/1902.05337)
- Differential phase PSD (Eq. 4): S_ψ(d,ω) = (2nk)²·[d²·S_h(ω)/4 + S_Δx(d,ω)]·|ωG(ω)|² + 2S_ε(ω), with sensitivity function |ωG(ω)| = 4 sin²(ωT/2) (Eq. 5) and S_Δx = S_Δa/ω⁴ (Eq. 6), where S_Δa is the PSD of the differential gravity acceleration projected on the baseline (this is the GGN).
- Strain sensitivity (SNR = 1, Eq. 7): `S_h(ω) = 4·S_Δa(d,ω)/(ω⁴ d²) + S_ε(ω)/(2·(2nk)²·d²·sin⁴(ωT/2))`.
- Differential-acceleration PSD from a spatial correlation C (Eq. 9): S_Δa(d,ω) = 2C(0,ω) − C(d,ω) − C(−d,ω).
- The GW-equivalent of a differential acceleration (from Eq. 7) is Δa = ω²hd/2, so h_n = 2√S_Δa/(ω²d).

**How sensitivity curves are drawn**
- AION/AEDGE plot (i) ASD-like strain curves h_n(f) from δφ_noise and the response; (ii) monochromatic Ω_GW h² sensitivity; (iii) power-law-integrated (PI) curves for stochastic backgrounds (SNR = 10, 5 yr); (iv) characteristic strain. The resonant-mode envelope is made of stair steps set by the n_max and T_max constraints — [Badurina et al. 2020](https://arxiv.org/abs/1911.11755); [Graham et al. 2016](https://arxiv.org/abs/1606.01860); [Badurina et al. 2021](https://arxiv.org/abs/2108.02468)
- For a single-axis (single-baseline) detector, averaging over sky and polarisation reduces power response by 4/15. Ross et al. multiply by 15/4 when converting a single-axis strain ASD to an Ω_GW limit. The single-axis antenna pattern is h_xx = h₊(cos²θcos²φ − sin²φ) + 2h×·cosθ·sinφ·cosφ — [Ross et al. 2020, arXiv:2002.02044](https://arxiv.org/abs/2002.02044)

### Inferences
- Link to satgrav's existing model: for a light-pulse AI with effective wavenumber k_eff, the acceleration-to-phase transfer is k_eff·(4/ω²)·sin²(ωT/2). That gives k_eff·a·T² at DC; the MZ double-diffraction case in satgrav uses 2k_eff. Feeding in the GW tidal differential acceleration Δa = ½ω²hL gives Δφ = 2·k_eff·h·L·sin²(ωT/2), exactly Graham's low-frequency, Q = 1 limit of Eq. 9 *(derived here; consistent with MIGA Eqs. 4–7)*. So a `satgrav` GW module needs only: (a) `h_to_delta_a(f, L) = ½(2πf)²L`; (b) the existing MZ transfer function; (c) optional finite-light-time factors sinc(ωnL/2c)·sin(ωT/2 − ω(n−1)L/2c) and the resonant factor sin(ωQT)/cos(ωT/2); (d) a sky-average factor (4/15 in power) for single-baseline detectors.
- The same formula converts any gradiometer's gradient noise √S_Γ (s⁻²/√Hz) to strain: h_n = 2√S_Γ/ω². For a GOCE-class 10 mE/√Hz gradiometer that gives h ≈ 5e-9 (10 mHz), 5e-11 (0.1 Hz), 5e-13 (1 Hz) per √Hz *(computed here)*. That is 8–12 orders of magnitude short of LISA or AION-km. A CARIOQA/GRICE-type space gradiometer is not a GW detector. The GW-relevant quantity is baseline length, which is why AEDGE uses L = 4.4e7 m between satellites rather than a gradiometer arm.
- For a Raman two-photon gradiometer (Rb), laser phase noise does not cancel over long baselines. Single-photon Sr transitions are what make km–40 000 km baselines viable ([Graham et al. 2016](https://arxiv.org/abs/1606.01860), introduction). This is why satgrav's Rb Raman parameters cannot simply be scaled to L = km.

### Gaps
- Graham et al. 2013 (PRL 110, 171102, arXiv:1206.0818), the original single-photon proposal, was not fetched; its formula is reproduced through the 2016 paper's Q = 1 limit.
- The exact sky-averaged response function used by Badurina et al. 2021 is not given in the fetched text (it refers to earlier work).
- The power-law-integrated curve method (Thrane & Romano 2013) is cited by AION but was not fetched; its formula is not reproduced here.

## 3. Gravity-gradient (Newtonian) noise: mechanisms, formulas, and links to geophysical gravimetry

### Takeaway
Gravity-gradient noise (GGN, also called Newtonian noise or NN) is the same physics as gravimetric forward
modelling of density changes. Seismic Rayleigh waves, P and S body waves, and atmospheric infrasound move
density around and produce a gravity acceleration δa that no shielding can block. For an AI gradiometer at
low frequency, the GGN strain is independent of baseline. Mitigation is by Wiener filtering from seismometer
and strainmeter arrays (residual ≈ 1/(N·SNR²)), multi-gradiometer arrays, and going underground (factor
e^(−kh)). Terrestrial AIs are expected to be GGN-limited below about 0.1 Hz.

### Cited Findings

**Rayleigh-wave GGN at depth h (MIGA formalism, following Harms)** — [Junca et al. 2019, arXiv:1902.05337](https://arxiv.org/abs/1902.05337)
- Gravity potential at depth h (Eq. 11): δU(r₀,t) = 2πGρ₀·(ξ_z(ω)/…)·e^(i(k_R·r₀−ωt))·[−2e^(−h q_z^S) + (1 + q_z^P/q_z^S)e^(−h k_R)]·(…), with q_z^{P,S} = √(k_R² − k_{P,S}²).
- Horizontal acceleration on a test mass at r₀ = X·e_x − h·e_z (Eq. 15): δa_x(X,t) = κ_R(k_R)·ξ_z(ω)·e^(−iωt)·e^(i k_R X cosφ)·i·cosφ, with `κ_R(k_R) = 2πGρ₀γ[e^(−h k_R) + b(e^(−h q_z^P) − e^(−h k_R))]`. Here γ and b are Rayleigh-wave factors depending on c_R, α (P speed) and β (S speed).
- Isotropic Rayleigh field correlation (Eq. 17): C_R(d,ω) = κ_R²·S_ξz(ω)·⟨cos²φ·e^(i k_R d cosφ)⟩_φ = (κ_R²S_ξz/2)·[J₀(k_R d) − J₂(k_R d)]. The ½ is reconstructed so that Eq. 9 reproduces Eq. 19. GGN (Eq. 18): `S_Δa(d,ω) = κ_R²(k_R)·S_ξz(ω)·χ_R(k_R,d)`, with χ_R = 1 − J₀(k_R d) + J₂(k_R d) (Eq. 19) ≈ (3/8)k_R²d² for f ≪ c_R/(2πd) (Eq. 20).
- Consequence: with Eq. 7, the seismic GGN strain limit is independent of gradiometer length at low frequency. MIGA example: d = 150 m, c_R ≈ 2.4 km/s (valid for f ≪ 2.5 Hz); LSBB parameters α ≈ 4.66 km/s, β ≈ 2.61 km/s, ρ₀ = 2500 kg/m³, depth 300 m. The seismic GGN is below the initial MIGA shot noise (2-photon, 1 mrad/s detection noise), and only becomes visible in the mHz range for the upgraded iMIGA (2×100 photon LMT, 0.1 mrad/s).

**Atmospheric infrasound GGN** — [Junca et al. 2019](https://arxiv.org/abs/1902.05337)
- Adiabatic sound wave: δρ/ρ₀ = δp/(γ_a p₀). Acceleration along the gradiometer (Eq. 25): δa_x = κ_I(k_I)·δp(ω)·e^(−iωt)·e^(i k_I X sinθ cosφ)·e^(−h k_I sinθ)·sinθ·i·cosφ. Here κ_I is printed as "4πGρ γp0 kI" in the extracted text; dimensional analysis gives `κ_I = 4πGρ₀/(γ_a p₀ k_I)`, with ρ₀ the air density, γ_a the adiabatic index and p₀ the mean pressure.
- GGN (Eqs. 27–28): S_Δa = κ_I²·S_δp(ω)·χ_I, with χ_I = ⟨2(1 − cos(k_I d sinθ cosφ))·e^(−2hk_I sinθ)·sin²θ·cos²φ⟩_(θ,φ), and χ_I ∝ k_I²d² for k_I d ≪ 1 and 2hk_I ≪ 1. Unlike the seismic case, around MIGA's peak (about 2 Hz) the full expression is needed and infrasound GGN depends on d. MIGA extrapolates the pressure PSD at low frequency as S(δp;ω) = S₀ω^(−2.2).

**Body-wave GGN in a cavity and Wiener-filter cancellation (Einstein Telescope context)** — [Badaracco & Harms 2019, arXiv:1903.07936](https://arxiv.org/abs/1903.07936)
- General formula (Eq. 1): δa(r₀,t) = −G∫dV ρ(r)·(ξ(r,t)·∇₀)·(r − r₀)/|r − r₀|³.
- Small cavity in an infinite homogeneous medium (Eq. 7): `δa(ω) = (4/3)πGρ₀·(2ξ^P(ω) − ξ^S(ω))`. With P fraction p (S(ξ^P) = p·S(ξ), S(ξ^S) = (1−p)·S(ξ)), the test-mass PSD is (Eq. 9) `C_NN = (4πGρ₀/3)²·S(ξ;ω)·(3p+1)`.
- Wiener residual (Eq. 2): `R(ω) = 1 − C_SN†·C_SS⁻¹·C_SN / C_NN`. For N ideal sensors of signal-to-noise ratio SNR (Eq. 3): R_min ≈ 1/(N·SNR²).
- Results for p = 1/3 and SNR = 15: 6 sensors give √R ≈ 0.43; √R < 0.1 needs N > 14 (about 15 sensors for 10× amplitude reduction); three-axis sensors are not significantly better than single-axis.
- 2025 update: "fusion" arrays of seismometers plus DAS-type fibre strainmeters. Six borehole seismometers plus 14 in-tunnel strainmeters give √R = 0.125, matching 20 borehole seismometers. Unconstrained 20-sensor optimum: √R = 0.082. Strain-only arrays reach only √R = 0.63 for pure S-wave fields. Gravimeters and tiltmeters are not used — [Ophardt, Badaracco, Isleif 2025, arXiv:2512.13554](https://arxiv.org/abs/2512.13554). Related: [joint body-wave optimisation for ET, arXiv:2310.05709](https://arxiv.org/abs/2310.05709); [cost-efficient robust ET arrays, arXiv:2606.21110](https://arxiv.org/abs/2606.21110); [full-tensor detector NN cancellation, arXiv:1504.04724](https://arxiv.org/abs/1504.04724)

**GGN in AI projects**
- AION models GGN with the USGS New Low-Noise Model (NLNM) and simple depth scaling, plus mitigation from N interferometers: S_h^GGN ∝ 1/(√N·L) (the √ placement is ambiguous in the extracted text). The seismic campaigns at Fermilab NuMI (MAGIS), LSBB (MIGA) and CERN (surface and LHC tunnel) are "in broad agreement" with the NLNM. Above about 0.1 Hz GGN is assumed fully mitigable; below about 0.1 Hz GGN "will start to dominate". In space, the GGN analogue comes from asteroids (below about 1e-6 Hz) and unresolved galactic and extragalactic binaries — [Badurina et al. 2021](https://arxiv.org/abs/2108.02468)
- AICE plots GGN for the NHNM, the NLNM, and PX46-measured seismic data (sandstone/molasse; glacial till gives a similar result). It plans multiple sources (5–10) to control GGN, "particularly in the presence of different rock strata", and possible active mitigation with "purpose-built instrumentation for seismic and weather monitoring" — [arXiv:2509.11867](https://arxiv.org/abs/2509.11867)
- Multi-gradiometer arrays to average GGN are the basis of ELGAR's 2D geometry — [Canuel et al. 2020](https://arxiv.org/abs/1911.03701) (not fetched)

### Inferences
- Combining MIGA Eqs. 7 and 20 gives a closed-form low-frequency Rayleigh GGN strain for an AI gradiometer: S_h^GGN ≈ (3/2)·κ_R²·k_R²·S_ξz/ω⁴ = (3/2)·κ_R²·S_ξz/(ω²c_R²), independent of d *(derived here)*. At the surface (h = 0, b-term ~0) κ_R ≈ 2πGρ₀γ. This is a compact function to add to satgrav, taking the NLNM/NHNM as input S_ξz.
- For a laser-interferometer arm with four test masses that are uncorrelated at separations L ≫ λ_R, and an isotropic field (⟨cos²φ⟩ = ½), the same κ_R gives S_h^NN ≈ 2κ_R²·S_ξz/(ω⁴L²) *(derived here from the MIGA single-mass expression; check against Harms' Living Review before use)*.
- Link to gravimetry: NN subtraction is gravimetric forward modelling of measured density changes (seismic displacement or pressure), followed by subtraction. The infrasound κ_I term is the dynamic, wavenumber-dependent generalisation of the barometric admittance correction applied to superconducting gravimeter records. Wiener filtering a gravimeter against a barometer array is the same estimator as Eq. 2 above. A satgrav "environmental correction" module could share code between (a) ground-gravimeter barometric and seismic corrections and (b) a GW-detector NN budget.
- For satellites the terrestrial GGN problem is replaced by time-variable gravity (tides, atmosphere/ocean de-aliasing, hydrology), which GRACE-FO models and removes below about 30 mHz ([Ross et al. 2020](https://arxiv.org/abs/2002.02044), their Fig. 1 note). That is satgrav's existing signal, and it is "noise" for a GW search.

### Gaps
- Harms, "Terrestrial gravity fluctuations", Living Rev. Relativ. 22:6 (2019), [doi:10.1007/s41114-019-0022-2](https://doi.org/10.1007/s41114-019-0022-2): the arXiv PDF (1507.05850) exceeded the 10 MB fetch limit and the Springer page redirects to a login. Exact equation numbers, the γ(ν) and b expressions, and Harms' surface-detector strain formula were therefore not verified. Use the MIGA formulas above, or retrieve the review manually.
- Saulson 1984 (PRD 30, 732) was not retrieved; its estimate is not reproduced.
- The exact infrasound low-frequency prefactor for χ_I was illegible in the extracted text.
- No quantitative LIGO/ET NN curve values were extracted.

## 4. Prompt elasto-gravity signals (PEGS) from earthquakes

### Takeaway
A rupture changes the mass distribution immediately, and the gravity change propagates at light speed, so
it arrives before the P wave. On the ground the measured PEGS is the gravity perturbation partly cancelled
by the ground acceleration it induces. Observed amplitudes are about 1 nm/s² (0.1 µGal) for Tohoku Mw 9.1
at 1000–2000 km. PEGS were first reported on the Kamioka superconducting gravimeter plus F-net (2016),
clearly on broadband networks in 2017, and since 2022 they have been extracted by deep learning for Mw ≳ 7.8.
The first operational PEGS tsunami-warning implementation was in Peru in 2025. No satellite detection of
PEGS has been published. GOCE has detected the later atmospheric infrasound from Tohoku.

### Cited Findings
- **Theory**: an analytic solution for gravity perturbations from a point double-couple in an infinite, homogeneous, non-self-gravitating elastic medium, valid between rupture onset and the P arrival. It was validated against spectral-element simulations of strike-slip and dip-slip half-space ruptures. Detection within 10 s of onset is "severely challenged" by microseism noise but may be possible with future gravity strainmeters — [Harms, Ampuero, Barsuglia, Chassande-Mottin, Montagner, Somala, Whiting 2015, GJI 201(3), 1416–1425, doi:10.1093/gji/ggv090](https://authors.library.caltech.edu/records/hs6t9-kcm59)
- **First report (Tohoku 2011, Mw 9.1)**: pre-P signals in the superconducting gravimeter (SG) at Kamioka (about 510 km from the epicentre) and five nearby F-net broadband seismometers. The theory-expected amplitude was "a fraction of µGal at ∼65 s" after onset — [Montagner et al. 2016, Nat. Commun. 7, 13349, doi:10.1038/ncomms13349](https://pubmed.ncbi.nlm.nih.gov/27874858/) (via search summaries)
- **Clear detection**: gravity perturbations are "particularly well observed" with broadband seismometers at 1000–2000 km from Tohoku. The figure scale bar is 1 nm/s² for stations at about 1390–3044 km — [Vallée et al. 2017, Science 358, 1164–1168, doi:10.1126/science.aao0746](https://www.science.org/doi/10.1126/science.aao0746) (via search summary)
- **Dense-array reanalysis (Japan)**: at Kamioka, the noise level A_N fell from 70 to 0.4 nm/s² after filtering, and Montagner et al. 2016 "could not visually detect a clear signal in the time domain" — [Kimura et al. 2019, Earth Planets Space 71:27, doi:10.1186/s40623-019-1006-x](https://link.springer.com/article/10.1186/s40623-019-1006-x) (via search summary)
- **Review**: [PEGS and their potential use in modern seismology, EPSL 2020](https://www.sciencedirect.com/science/article/abs/pii/S0012821X20300935)
- **Deep learning (PEGSNet)**: instantaneous tracking of earthquake growth from PEGS — [Licciardi et al. 2022, Nature, doi:10.1038/s41586-022-04672-7](https://www.nature.com/articles/s41586-022-04672-7). The follow-on work gives magnitude and focal mechanism for Mw ≥ 7.8 events about 2 minutes after origin, using the Alaska and Western Canada broadband networks — [Commun. Earth Environ. 2024, doi:10.1038/s43247-024-01725-9](https://www.nature.com/articles/s43247-024-01725-9) (via search summary)
- **2025–2026**: PEGSGraph, a graph neural network that uses the whole station network once 30 stations detect a P wave, updating magnitude every 10 s — [Hourcade et al. 2025, JGR: Machine Learning and Computation, doi:10.1029/2024JH000360](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2024JH000360). First implementation of a PEGS-based tsunami warning system, in Peru — [GJI 244(1), ggaf419](https://academic.oup.com/gji/article/244/1/ggaf419/8300379); [Temblor summary](https://temblor.net/earthquake-insights/upgrading-tsunami-warning-systems-for-faster-and-more-accurate-alerts-16975/). Sensitivity of PEGS to 3-D Earth structure — [GJI 245(1), ggaf519](https://academic.oup.com/gji/article/245/1/ggaf519/8382458)
- **Detectability with future gravity gradiometers and strainmeters**: conventional inertial sensors (gravimeters, seismometers) can "almost" never detect prompt gravity from M < ~8 within a few tens of seconds; future gravity strainmeters could detect M > 7 — [Juhel et al. 2018, JGR Solid Earth, doi:10.1029/2018JB016698](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018JB016698); [Shimoda et al. 2021, GJI 224, 533](https://academic.oup.com/gji/article/224/1/533/5920617) (via search summary). A TOBA-class gradiometer at 1e-15/√Hz at 0.1 Hz could detect an Mw 7 at about 100 km, about 10 s earlier than the current warning system, with two detectors 75 km apart for location — [Takano et al. 2024, arXiv:2412.01323](https://arxiv.org/abs/2412.01323)
- **Satellites**: GOCE (about 270 km altitude) detected the Tohoku *infrasound* (not PEGS): up to 11% air-density perturbation and 1.35e-7 m/s² vertical platform acceleration, first seen about 27 minutes after the earthquake. Detection above about 300 km is expected to be hard because infrasound fades with altitude — [Garcia et al. 2013, GRL, "GOCE: the first seismometer in orbit around the Earth"](https://hal-insu.archives-ouvertes.fr/insu-02563854v1) (via search summary). GRACE-FO's laser ranging interferometer (LRI) "captures gravitational signals as small as 0.1 nm/s² at 490 km altitude" (about 10× better than KBR) — [Ghobadi-Far et al. 2020, GRL, doi:10.1029/2020GL089445](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020GL089445)
- No PEGS study of the 29 July 2025 Kamchatka Mw 8.8 earthquake was found. The Kamchatka arXiv papers found concern waveform cross-correlation at IMS arrays — [arXiv:2606.17060](https://arxiv.org/pdf/2606.17060), [arXiv:2607.16226](https://arxiv.org/pdf/2607.16226)

### Inferences
- Order of magnitude: δψ ~ G·I₂[M]/r³ and δg ~ G·I₂[M]/r⁴, with I₂ = ∫u·M(t−u)du (dimensional form consistent with Harms 2015's "compact" solution; exact angular factors not verified). For an instant step M₀ this gives δg ~ G·M₀·t²/(2r⁴): for Mw 9.1 (M₀ ≈ 5.6e22 N·m) at 1000 km just before the P arrival (t ≈ 125 s), δg ~ 3e-8 m/s² *(computed here)*. This is an upper bound. Tohoku's moment release spreads over about 150 s, so M(t) is much smaller at early times, and on the ground the induced acceleration cancels much of the perturbation. Together these explain the observed ~1 nm/s².
- Satellites are in free fall and do not suffer the ground-acceleration cancellation; they see the gravity perturbation itself. However (a) at 490 km altitude the short-wavelength part is attenuated; (b) GRACE-FO measures the *difference* along a 220 km baseline; (c) the signal lasts minutes, and the pair travels about 7.6 km/s, so a useful geometry happens only by chance; (d) the LRI floor of about 0.1 nm/s² is quoted for its gravity-difference analysis bandwidth, not for a minute-long transient. A PEGS detection with GRACE-FO is plausible only for Mw ≳ 9 directly beneath a fortunate overpass, and would be marginal. The coseismic *static* step, already used in the repo's E1 Kamchatka work, is the practical satellite observable. Treat this as a hypothesis for a pre-registered test, not a known result.
- A satgrav PEGS module could reuse the AI gradiometer noise model. For a gravity-gradient (TOBA/AI) sensor, the PEGS observable is the gradient δΓ ~ δg/r, and for an SG it is δg minus the induced ground acceleration.

### Gaps
- The explicit Harms 2015 early-time formula (with radiation-pattern factors and the I₂ integral) could not be retrieved: the IPGP host refused connections and the Caltech S3 link timed out. Its exact form should be checked before implementation.
- Vallée et al. 2017 full text (exact amplitude-vs-distance values, the cancellation formula δg − ü) was not fetched.
- No published satellite PEGS detection or feasibility study was found.

## 5. Torsion-bar antennas and superconducting/atom gravity gradiometers as low-frequency GW detectors

### Takeaway
TOBA responds to the cross-polarised tidal field through the bar's dynamic quadrupole:
θ(f) = (q/2I)·f²/(f₀²(1+iφ) − f²)·h(f). Above the torsion resonance it is a free-mass gravity-gradient
sensor. The demonstrated sensitivity is 1e-8 (2010) and 1e-10/√Hz (2015). Phase-III targets 1e-15/√Hz at
0.1 Hz, mainly for earthquake early warning and Newtonian-noise studies, and the final 10 m TOBA targets
1e-19. SOGRO (superconducting tensor gradiometer) targets 1e-19–1e-21/√Hz at 0.1–10 Hz but is a design
only. Any gradiometer converts to strain via h = 2√S_Γ/ω².

### Cited Findings
- **TOBA equations** — [Takano et al. 2024, "TOrsion-Bar Antenna: A Ground-Based Detector for Low-Frequency Gravity Gradient Measurement", arXiv:2412.01323](https://arxiv.org/abs/2412.01323):
  - GW torque (Eq. 9): N_GW = −∂U_GW/∂θ = ¼·ḧ_ij·q^ij; in polarisations (Eq. 15): N_GW = (ḧ₊/2)q₊ + (ḧ×/2)q×.
  - Dynamic quadrupole (Eq. 10): q^ij = ∫dV ρ(ξ^i w^j + ξ^j w^i − ⅔δ^ij ξ^k w_k).
  - Equation of motion (Eq. 16): I·θ̈ + κθ = N_GW; resonance (Eq. 18) f₀ = (1/2π)√(κ/I).
  - Transfer function (Eqs. 19–21): θ(f) = Σ_A H_A(f)h_A(f), with `H_A(f) = (q_A/2I)·f²/(f₀²(1+iφ_rot) − f²)`. For a bar along x (Eq. 25): H₊ = 0 and H× ≈ ½·f²/(f₀²(1+iφ_rot) − f²), with I = M(l_x² + l_y²)/12 and q× = M(l_x² − l_y²)/12. Terrestrial gravity gradients cannot be distinguished from GW effects.
  - **Demonstrated**: Phase-I (2010) about 1e-8/√Hz at 0.1 Hz; Phase-II (2015) about 1e-10/√Hz at 3–8 Hz; seismic cross-coupling 5e-6 rad/m at 0.1 Hz (2017); cryogenic suspension to 6.1 K in 10 days (2020); cryogenic monolithic interferometer 3.6e-14 m/√Hz at 0.1 Hz at 12 K (2024).
  - **Targets**: Phase-III TOBA about 1e-15/√Hz at 0.1 Hz with 30 cm test masses; final TOBA 1e-19/√Hz at 0.1 Hz with a 10 m detector.
- The earlier overview gives the same target structure. The final TOBA is limited by quantum noise (shot and radiation-pressure) with 10 m cryogenic bars at 4 K and 10 W laser power. 1e-15/√Hz at 0.1 Hz is enough for geophysics (Mw 7 early warning) and for testing NN models and their cancellation, since NN reaches about 1e-15/√Hz near 0.1 Hz. Phase-II rotation was read with Fabry–Perot cavities — [Shimoda et al. 2019, arXiv:1812.01835](https://arxiv.org/abs/1812.01835)
- Original proposal: [Ando, Ishidoshiro, Yamamoto, Yagi, Kokuyama, Tsubono, Takamori 2010, PRL 105, 161101](https://toyama-kaken.devel.elsevierpure.com/en/publications/torsion-bar-antenna-for-low-frequency-gravitational-wave-observat) (no arXiv version found)
- **SOGRO (Superconducting Omni-directional Gravitational Radiation Observatory)**: six magnetically levitated superconducting test masses on orthogonal axes with SQUID differential readout. Target strain 1e-19–1e-21/√Hz at 0.1–10 Hz; seismic common-mode rejection of 1 part in 1e9; NN reduction via the full-tensor response. This is a design concept with no built instrument found — [Paik et al., APS April Meeting 2018 abstract D04.3](https://meetings-archive.aps.org/apr/2018/d04/3) (conference abstract, not peer reviewed)
- **Atom gradiometer**: MIGA Eq. 7 above gives S_h for an AI gradiometer, including GGN and detection noise — [Junca et al. 2019](https://arxiv.org/abs/1902.05337)

### Inferences
- Above f₀, TOBA's response tends to H× → −½, so θ ≈ −h×/2. Equivalently, the bar's rotation measures the gravity-gradient component Γ_xy(t) = −½ḧ_xy, consistent with the gradiometer conversion h = 2√S_Γ/ω² used in Section 2 *(derived here)*. TOBA, SOGRO, MIGA and satgrav's satellite gradiometers can all share one `GradientNoise → StrainNoise` conversion, with sensor-specific responses below their resonances.
- Order of magnitude: Phase-III's 1e-15/√Hz at 0.1 Hz corresponds to √S_Γ ≈ ½·1e-15·(2π·0.1)² ≈ 2e-16 s⁻²/√Hz ≈ 2e-7 E/√Hz *(computed here)*. That is about 5 × 10⁴ times better than GOCE's 10 mE/√Hz, which shows why ground torsion bars, not space gradiometers, are the PEGS/NN instruments.

### Gaps
- Paik et al. 2016 (CQG 33, 075003, the SOGRO design paper) was not fetched; the SOGRO numbers above come only from conference abstracts.
- Superconducting-gravimeter-based GW limits (e.g. Coughlin & Harms bounds, cited by [Ross et al.](https://arxiv.org/abs/2002.02044) as refs [10, 12]) were not retrieved.

## 6. Space clocks and satellite ranging as GW detectors (optical-lattice clocks, LISA, GRACE-FO)

### Takeaway
LISA was adopted by ESA on 25 January 2024 (2.5 Gm arms, launch planned for 2035). Its standard analytic
sensitivity is Robson, Cornish and Liu's S_n(f), with explicit P_OMS, P_acc and a galactic-confusion fit.
Optical-lattice clocks on two drag-free satellites sharing one laser link see s = δν/ν = h·sin(πfd/c).
With 7e6 Sr atoms that gives about 1e-20/√Hz in fractional frequency (a target 4 orders beyond the
clock stability demonstrated in 2016). GRACE-FO's LRI, read as strain over 220 km, constrains
Ω_GW < 3.3e7 at 100 mHz, with a best characteristic strain of 2e-15 at 0.5 Hz. It is a real but
astrophysically useless bound.

### Cited Findings

**LISA status** — [ESA, "Capturing the ripples of spacetime: LISA gets go-ahead"](https://www.esa.int/Science_Exploration/Space_Science/LISA/Capturing_the_ripples_of_spacetime_LISA_gets_go-ahead); [ESA LISA factsheet](https://www.esa.int/Science_Exploration/Space_Science/LISA_factsheet)
- ESA's Science Programme Committee adopted LISA on 25 January 2024. It had been selected as L3 on 20 June 2017. Construction was to begin in 2025 once a prime contractor was chosen. Launch is planned for 2035 on Ariane 6. The design has three spacecraft with 2.5 million km arms, and the band is 0.1 mHz–1 Hz.

**LISA analytic sensitivity curve** — [Robson, Cornish & Liu 2019, CQG 36, 105011, arXiv:1803.01944](https://arxiv.org/abs/1803.01944)
- Full form (Eq. 13 plus confusion noise, Eq. 1):
  `S_n(f) = (10/(3L²))·[P_OMS(f) + 2(1 + cos²(f/f*))·P_acc(f)/(2πf)⁴]·(1 + (6/10)(f/f*)²) + S_c(f)`,
  with L = 2.5e9 m, f* = c/(2πL) = 19.09 mHz.
- `P_OMS = (1.5e-11 m)²·(1 + (2 mHz/f)⁴) Hz⁻¹` (Eq. 10); `P_acc = (3e-15 m s⁻²)²·(1 + (0.4 mHz/f)²)·(1 + (f/8 mHz)⁴) Hz⁻¹` (Eq. 11). These come from the LISA "Science Requirements / Sensitivity Curves" note LISA-LCST-SGS-TN-001.
- Sky-averaged transfer function fit (Eq. 9): R(f) ≈ (3/10)/(1 + 0.6(f/f*)²). The 3/10 sums LISA's two low-frequency channels; many papers quote 3/20 for a single channel.
- Galactic confusion (Eq. 14): `S_c(f) = A·f^(−7/3)·exp(−f^α + β f sin(κf))·[1 + tanh(γ(f_k − f))] Hz⁻¹`, with A = 9e-45 and:

  | T_obs | α | β | κ | γ | f_k [Hz] |
  |---|---|---|---|---|---|
  | 6 mo | 0.133 | 243 | 482 | 917 | 0.00258 |
  | 1 yr | 0.171 | 292 | 1020 | 1680 | 0.00215 |
  | 2 yr | 0.165 | 299 | 611 | 1340 | 0.00173 |
  | 4 yr | 0.138 | −221 | 521 | 1680 | 0.00113 |

- Check values for √S_n without S_c *(computed here)*: 1.28e-19 (1 mHz), 1.74e-20 (3 mHz), 1.20e-20 (10 mHz), 4.6e-20 (0.1 Hz) per √Hz.

**Optical-lattice-clock GW detector** — [Kolkowitz, Pikovski, Langellier, Lukin, Walsworth, Ye 2016, PRD 94, 124043, arXiv:1606.01859](https://arxiv.org/abs/1606.01859)
- Two drag-free satellites, each with a Sr lattice clock, with the laser in B phase-locked to light from A (one shared clock laser). Signal (Eq. 1): `s ≡ δν/ν = h·sin(π f_GW d/c)`, so s = h at d = λ_GW/2.
- The clock transfer function scales as T_ν ∝ f²·T_φ, i.e. ∝ f² below c/2d and flat above (where an optical interferometer's T_φ ∝ 1/f). The chosen d = 5e10 m gives a peak-sensitivity floor at c/2d ≈ 3 mHz; a LISA-like 5e9 m works above about 30 mHz.
- Shot-noise-limited fractional frequency (Eq. 2): σ_min(τ) = δν_min/ν (per √τ). The extracted text loses the root; the form `σ_min = √Δ_A/(ν√(2πτN))` reproduces their number *(checked here)*. With ⁸⁷Sr at ν = 430 THz, Δ_A = 1 mHz (T_max = 160 s) and N = 7e6, σ_min = 1.1e-20/√Hz, a 4-orders-of-magnitude improvement over the clock stability demonstrated in 2016. Dynamical-decoupling (DD) sequences give narrow-band sensitivity to monochromatic GWs.
- A link-related estimate of h_min ≈ 3.7e-23 for a continuous GW at 3 mHz–10 Hz (with ΔL = 30 mHz, η = 0.5, ≤ 1 day averaging) appears in the paper. Its context (laser-link/photon-number requirement) was only partially extracted.

**GRACE-FO as a GW detector** — [Ross, Hagedorn, Shaw, Lockwood, Iritani, Lee, Venkateswara, Gundlach 2020, arXiv:2002.02044](https://arxiv.org/abs/2002.02044)
- Strain from the published LRI displacement spectrum: h̃_xx(f) = x̃(f)/x₀ with x₀ = 220 km. Noise below 30 mHz is terrestrial gravity, not instrument noise. Single-axis antenna pattern as in Section 2. The orbit modulates the antenna pattern at f_orb = 1/94.5 min (φ(t) = 2πf_orb t + φ₀).
- Stochastic limit (Eq. 3): Ω_GW(f) = (prefactor/H₀²)·f³·h̃_meas(f)²·(15/4), with H₀ = 70.3 km/s/Mpc and 95th-percentile bins of 20 points. The prefactor reads "32π³/3" in the extracted text; the standard one-sided-PSD conversion is 2π²f³S_h/(3H₀²), so check the original before coding. Result: **Ω_GW < 3.3e7 at 100 mHz**, "many orders of magnitude" above expected sources.
- Best characteristic-strain sensitivity is 2e-15 at 0.5 Hz. With a dedicated analysis (30 mHz–5 Hz) GRACE-FO might see local neutron-star binaries, sub-galactic stellar-mass BH binaries, or Milky Way IMBH mergers. A LISA Pathfinder analogue (proof masses 376 mm apart, acceleration noise converted to strain) is given in the appendix.
- The LRI is LISA-heritage technology, but its science case is Earth gravity — [NASA GRACE-FO LRI page](https://gracefo.jpl.nasa.gov/resources/14/laser-ranging-interferometer/)

### Inferences
- For satgrav, LISA-RCL is a self-contained closed form (two noise terms, a transfer factor and a 4-row lookup table). It is the natural reference curve to put next to satgrav's GRACE-FO/NGGM and AI-gradiometer curves on one strain axis. A unit test against the check values above would catch transcription errors.
- The GRACE-FO → strain mapping h = x/x₀ (displacement ASD over separation) also applies to NGGM/MAGIC-type pairs. Their LRI noise and 100–200 km separations put them around 1e-15–1e-14/√Hz *(order-of-magnitude inference)*. They remain useless as GW detectors, but a "strain view" is a cheap, honest cross-check for the crate.
- Clocks and AIs share the single-photon Sr clock transition. Kolkowitz's σ_min and AION/AEDGE's δφ_noise are both atom-projection-noise limits (∝ 1/√N). A shared `AtomEnsemble {N, contrast, coherence}` type could feed both the existing satgrav AI model and a clock-GW model.

### Gaps
- The LISA Definition Study Report (Red Book, ESA-SCI-DIR-RP-002, 2024; arXiv:2402.07571) and any updated 2024–2026 requirement curve were not retrieved. RCL is based on the older LISA-LCST-SGS-TN-001 levels.
- No 2023–2026 optical-clock GW proposal newer than Kolkowitz 2016 was found or verified (e.g. Loeb & Maoz 2015, or space-clock missions).
- The publication venue of Ross et al. 2020 was not confirmed (arXiv v3, May 2020). No newer GRACE-FO GW analysis was found.
- No published analysis was found of NGGM/MAGIC (Next Generation Gravity Mission, ESA–NASA) as a GW constraint.
