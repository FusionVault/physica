# Li–Baker HFGW detector and generator mathematics, the JASON critique, and the Pais Navy patents

Conventions in these notes:
- **Units.** SI unless marked "Gaussian" (JASON writes in Gaussian/cgs). Symbols: ħω is photon or graviton energy, Z0 = μ0c ≈ 376.7 Ω, A or h is the GW strain amplitude, B̂ is a static field, ψ0 is the Gaussian-beam (GB) peak field.
- **Status labels.** **[mainstream]** agrees with standard GR and electrodynamics. **[contested]** is disputed in the literature with no settled refutation. **[refuted]** has a published refutation, named in each case.
- **"Our check".** A number recomputed in Python (`python3 -I`, CODATA constants) from the authors' own formula. The ratio is our value divided by the published value. Every "our check" line is an inference by the note author, not a published claim.
- **Downloads.** All primary PDFs were read in full text (pdftotext) from arXiv, the authors' own website, Google Patents and an archive.org copy of the JASON report. The live FAS copy sits behind a bot challenge, which was not bypassed.

## Q1. Which papers did Li, Baker et al. publish? Key equations, parameters and predicted numbers

### Takeaway
Every Li-group detector paper from 2003 to 2023 computes the same quantity. They call it the first-order "perturbative photon flux" (PPF): the cross term between the inverse-Gertsenshtein field the GW induces in a static B field and a synchro-resonant 5–10 GHz Gaussian beam. The standard parameter set (10 W beam, W0 = 5 cm, 3 T, 6 m, 5 GHz, h = 10⁻³⁰) yields about 8×10² signal photons/s with a claimed integration time of 10³–10⁵ s. Their arithmetic reproduces from their own formulas to within a few percent. The physics interpretation is what JASON refutes (Q2).

### Cited Findings
**Papers (peer-reviewed unless noted)**
- **Li, Tang, Shi, Phys. Rev. D 67, 104008 (2003), arXiv:gr-qc/0306092, DOI 10.1103/PhysRevD.67.104008.** "Electromagnetic response of a Gaussian beam to high-frequency relic gravitational waves in quintessential inflationary models." A GB in a static B field under the "synchroresonance" condition ω_e = ω_g. For ν_g = 10¹⁰ Hz and h = 10⁻³⁰ the abstract predicts a PPF of about 10³ s⁻¹ through a 10⁻² m² region. — [arXiv:gr-qc/0306092](https://arxiv.org/abs/gr-qc/0306092)
- **Li, Baker, Chen (2006), arXiv:gr-qc/0604109 (preprint; the Living Review cites it as an arXiv e-print).**
  - Parameters: P = 10 W, ψ0 ≈ 1.26×10³ V/m, W0 = 0.05 m, B̂_y = 3 T, Δs ≈ 10⁻² m² (0 ≤ y ≤ W0, 0 ≤ z ≤ 0.3 m), l = 0.3 m, ν = 10¹⁰ Hz.
  - Strains: A⊗ = 2.00×10⁻³¹ (QIM) and 6.32×10⁻³¹ (PBBS).
  - Predicted PPF: about 4.04×10² s⁻¹ and 1.27×10³ s⁻¹ respectively. The maximum is at x = ±3.16 cm. At x = 0 the PPF is 1.15×10³ s⁻¹ while N_x^(0) = 0.
  - Source: [arXiv:gr-qc/0604109](https://arxiv.org/abs/gr-qc/0604109)
- **Li, Baker, Fang, Stephenson, Chen, Eur. Phys. J. C 56, 407–423 (2008), arXiv:0806.1989.** This is the "Li–Baker Chinese HFGW detector" paper and JASON's reference [11]. Source for every item below: [arXiv:0806.1989](https://arxiv.org/abs/0806.1989)
  - **Gaussian beam, eq. (27):**
    ψ = ψ0/√(1+(z/f)²) · exp(−r²/W²) · exp{i[(k_e z − ω_e t) − tan⁻¹(z/f) + k_e r²/(2R) + δ]}
    with r² = x²+y², f = πW0²/λ_e, W = W0[1+(z/f)²]^½ and R = z + f²/z. The static field is B̂_y^(0) for −l1 ≤ z ≤ l2 (eq. 30).
  - **Background photon flux (BPF) density, eq. (31):**
    n_x^(0) = (1/ħω_e)(1/μ0)⟨Ẽ_y^(0) B̃_z^(0)⟩ = (1/(2μ0ħω_e)) Re{ψ_y* (i/ω_e)(∂ψ_x/∂y − ∂ψ_y/∂x)} = f_x^(0) exp(−2r²/W²)
    in photons s⁻¹ m⁻². Because Ẽ_y^(0) = 0 at x = 0, n_x^(0)|_{x=0} = 0.
  - **PPF density, eq. (59):**
    n_x^(1) = (c/ħω_e)⟨T^01⟩_{ω_e=ω_g} = (1/(2μ0ħω_e)) Re{Ẽ_y^(1)* (i/ω_e)(∂ψ_x/∂y − ∂ψ_y/∂x)}_{ω_e=ω_g}
    The leading term is ∝ A⊗ B̂_y ψ0 k_g y (z+l1)/(4μ0 [1+(z/f)²]^½ (z+f²/z)). The full expression includes the quasi-probability integrals F1(y) = ∫exp(−y²/W²)cos(k_g y²/2R)dy and F2(y) = ∫exp(−y²/W²)sin(k_g y²/2R)dy (eq. 60).
  - **Signal, eq. (61):** N_x^(1) = ∬_Δs n_x^(1)|_{x=0} dy dz.
  - **Stated parameters:**
    - (1) P = 10 W, giving ψ0 ≈ 1.26×10³ V/m for W0 = 0.05 m.
    - (2) B̂_y = 3 T.
    - (3) Δs ≈ 10⁻² m².
    - (4) l = l1 + l2 = 0.3 m and 6 m.
    - (5) ν_e = ν_g = 5 GHz (λ_g = 0.06 m).
  - **Table I (Δν = 1 Hz, l = 6 m):** N_x^(1)(x=0) is 8.21×10¹⁰, 8.21×10⁸, 8.21×10⁶, 8.21×10⁴, 8.21×10², 8.21 and 0 s⁻¹ for h_rms = 10⁻²², 10⁻²⁴, …, 10⁻³² and 10⁻³⁴. N_x^(0)(x=0) = 0. At x = 3.5 cm, N_x^(0) = 1.24×10²² s⁻¹ and N_x^(1) = 3.54×10² s⁻¹ (for h = 10⁻³⁰).
  - **Graviton flux:** a plane HFGW with h = 10⁻³⁰ at 5 GHz carries N_g = 3.77×10¹⁶ s⁻¹ through the beam waist (W0 = 5 cm).
- **Li, Yang, Fang, Baker, Stephenson, Wen, Phys. Rev. D 80, 064013 (2009), arXiv:0909.4118, DOI 10.1103/PhysRevD.80.064013.** "Signal photon flux and background noise in a coupling electromagnetic detecting system for HFGWs." Source for every item below: [arXiv:0909.4118](https://arxiv.org/abs/0909.4118)
  - **Pure inverse-Gertsenshtein ("G-effect") part.** Parameters (eq. 6): B̂_y = 10 T, L = 10 m, ν = 5 GHz, k ≈ 100 m⁻¹, hν = 3.3×10⁻²⁴ J, Δs = 0.01 m².
    - Eq. (3): Ẽ^(1) ≈ A B̂_y k_g c z exp[i(k_g z − ω_g t)]
    - Eq. (5): u_em = (1/μ0)(A B̂_y k_g L)² c
    - Eq. (7): U_em^(2) = u_em Δs ≈ 2.3×10⁻⁴⁰ W at h = 10⁻³⁰
    - Eq. (8): N_γ^(2) ≈ 7.0×10⁻¹⁷ s⁻¹
    - Eq. (9): U_gw = (c³/8πG) ω² A² Δs ≈ 1.6×10⁻⁷ W
    - Eq. (10): N_g ≈ 4.8×10¹⁶ s⁻¹
    - Eq. (11): conversion ratio ≈ 1.4×10⁻³³
    - Eq. (12): Δt ≈ 1.4×10¹⁶ s
    - Eq. (13): h = 10⁻²⁶ gives 7×10⁻⁹ s⁻¹ (1.4×10⁸ s); h = 10⁻²⁴ gives 7×10⁻⁵ s⁻¹ (1.4×10⁴ s).
    - The authors' own conclusion: a pure inverse-G-effect lab detector "would not be useful" in the microwave band.
  - **Plane EM wave plus static field.**
    - Eq. (14): 10 W over 0.01 m² gives |E_x^(0)| ≈ 8.7×10² V/m.
    - Eq. (15): N^(0) ≈ 3.0×10²⁴ s⁻¹.
    - Eq. (17): U_z^(1) = Re[(1/μ0c) E^(1)* E_x^(0)] Δs ≈ 6.9×10⁻²⁰ W.
    - Eq. (18): N_z^(1) ≈ 2.1×10⁴ s⁻¹.
    - Eq. (20): ratio N_z^(1)/N_z^(0) ≈ 7.0×10⁻²¹.
    - Eq. (23): visibility condition N_z^(1)(Δt)^½ > (N_z^(0))^½, so at h = 10⁻²⁶ one needs Δt > 6.8×10⁷ s.
  - **Interference term.** Eqs. (32)–(35) write the cross term as n^(1) = 2(Ṅ0 Ṅ_GW)^½ cos δ, and the total count as N_d = Ṅ0Δt + 2(Ṅ0Ṅ_GW)^½ cos δ·Δt + Ṅ_GWΔt. This is the same expression as JASON eq. (3-32).
  - **Gaussian-beam scheme.**
    - Eq. (50): N_x^(1)Δt > (N_x^(0)Δt)^½
    - Eq. (51): Δt > N_x^(0)/(N_x^(1))² = Δt_min
    - Eqs. (52)–(54): N_x^(1) = N_x,max^(1) exp(−x²/W²); N_x^(0) ∝ x exp(−2x²/W²); N_x(other)^(0) ∝ exp(−2x²/W²)
    - Eq. (56): Δt_min = (x N_x,max^(0) + |N_x(other)^(0)|_max)/|N_x^(1)|²_max, with N_x,max^(0) ≈ 1.2×10²² s⁻¹
    - Eq. (57), at x = 0: Δt_min = |N_x(other)^(0)|_max / |N_x^(1)|²_max
    - Eq. (58), for B̂ = 3 T, L = 6 m, P = 10 W, h = 10⁻³⁰: N_x^(1) ≈ 8.2×10² s⁻¹. Then Δt_min ≈ 3.0×10³ s if N_other < 2.1×10⁹ s⁻¹, or ≈ 3.0×10⁵ s (3.5 days) if N_other < 2.1×10¹¹ s⁻¹.
  - **Off-axis receiving surfaces** (eqs. 62–65), for h = 10⁻²⁶:
    - x = 2 cm: N_x^(0) ≈ 1.7×10²⁰ s⁻¹, N_x^(1) ≈ 7.0×10⁶ s⁻¹, Δt_min ≈ 3.5×10⁶ s.
    - x = 10 cm: N_x^(0) ≈ 4.0×10¹⁷ s⁻¹, N_x^(1) ≈ 1.5×10⁵ s⁻¹, Δt_min ≈ 1.2×10⁷ s.
  - **Abstract:** the required accumulation time is about 10³–10⁵ s for h_rms ≈ 10⁻²⁶–10⁻³⁰ at 5 GHz in a 1 Hz bandwidth.
- **Li Xin, Wang Sai, Wen Hao, Chinese Physics C 40(8), 085101 (2016), arXiv:1509.02443, DOI 10.1088/1674-1137/40/8/085101.**
  - Detector parameters: B̂_y = 10 T, l1 = l2 = 0.3 m, ψ0 ≈ 1.26×10³ V/m, W0 = 0.05 m, Δs = 10⁻² m², evaluated at x = 3.5 cm.
  - As the tensor index n_t runs from −6.25×10⁻³ to 0.87: h_t(5 GHz) runs from 10⁻³⁶ to 10⁻²⁵, N_x1 from 1.40×10⁻⁴ to 2.85×10⁷ s⁻¹, and N_x1/BPF from 10⁻²⁸ to 10⁻¹⁶.
  - Sources: [arXiv:1509.02443](https://arxiv.org/abs/1509.02443); [CPC](https://cpc.ihep.ac.cn/article/doi/10.1088/1674-1137/40/8/085101)
- **Li F.-Y., Wen, Fang, Li D., Zhang T.-J., Eur. Phys. J. C 80, 879 (2020), arXiv:1712.00766, DOI 10.1140/epjc/s10052-020-08429-2 ("3DSR" paper).**
  - Gives analytic EM perturbations from all six possible GW polarizations in the "three-dimensional synchro-resonance system" (3DSR).
  - Claims all six are separable and detectable. The 3DSR band (10⁸–10¹² Hz) overlaps FAST (7×10⁷–3×10⁹ Hz).
  - Eqs. (28)–(29) carry a factor ½: Ẽ_y^(1) = (i/2) A⊗ B̂_y^(0) k_g c z exp[i(k_g z − ω_g t)].
  - Sources: [arXiv:1712.00766](https://arxiv.org/abs/1712.00766); [EPJC](https://link.springer.com/10.1140/epjc/s10052-020-08429-2)
- **Wu M., Li J., Jiang Q., Chinese Physics C 47(10), 105104 (2023), DOI 10.1088/1674-1137/ace9c3.** "Simulation of HFGW detection using modulated Gaussian beam" (a Li–Baker-type detector).
  - Parameters: ψ0 = 1260 (≈10 W), W0 = 0.06 m, f_e = f_g = 10⁹ Hz, h = 10⁻³¹.
  - Results: PPF peak about 2700 photons; BPF peak 1.57×10²⁷ photons; SNR 259 at x = 0.08 m and 320 at x = 0.064 m.
  - It argues the BPF "can be eliminated" at special transverse positions. It cites no JASON report.
  - Source: [CPC 2023](https://cpc.ihep.ac.cn/article/doi/10.1088/1674-1137/ace9c3)
- **Other Li-group items.**
  - Li & Yang, Chin. Phys. Lett. 21, 2113–2116 (2004), arXiv:gr-qc/0410060.
  - Woods, "Diffraction from embedded reflectors in Li-Baker HFGW detector," Physics Procedia 38, 66–76 (2012).
  - Both as listed in the Living Review bibliography (refs 538 and 542): [arXiv:2501.11723](https://arxiv.org/abs/2501.11723)
- **No JASON response in the later papers.** A full-text search of the EPJC 2008, PRD 2009, CPC 2016 and EPJC 2020 texts finds zero occurrences of "JASON". The EPJC 2020 "Eardley" hit is a different paper (Eardley, Lee, Lightman et al., PRL 30, 884, 1973). — [arXiv:0909.4118](https://arxiv.org/abs/0909.4118); [arXiv:1712.00766](https://arxiv.org/abs/1712.00766)

**Baker's AIP Conference Proceedings (STAIF) and related conference papers.** All are as listed in Baker's own reference lists: [Baker SPESIF-2009 refs](https://www.drrobertbaker.com/docs/HFGW%27s%20Power%20Requirements.pdf); [laser paper header](https://www.drrobertbaker.com/docs/AIP;%20HFGW%20Laser%20Generator.pdf); [Physics Procedia 2012 header](https://inspirehep.net/files/6bf245c473964d406575a245362a0ede)
- Baker, Woods, Li, "Piezoelectric-crystal-resonator HFGW generation and synchro-resonance detection," AIP Conf. Proc. 813, 1280–1300 (STAIF-2006).
- Baker, Li, Li R., "Ultra-high-intensity lasers for gravitational wave generation and detection," AIP Conf. Proc. 813, 1249–1258 (STAIF-2006).
- Fontana & Baker, "Generation of gravitational waves with nuclear reactions," AIP Conf. Proc. 813, 1352–1358 (2006).
- Baker, Stephenson, Li, "Analyses of the frequency and intensity of laboratory generated HFGWs," AIP Conf. Proc. 969, 1036–1044 (STAIF-2008).
- Baker, Stephenson, Li, "Proposed ultra-high sensitivity HFGW detector," AIP Conf. Proc. 969, 1045–1054 (STAIF-2008). This is JASON's reference [10].
- Baker, "Input power requirements for HFGW generators," and Baker & Black, "Radiation pattern for a multiple-element HFGW generator," both SPESIF-2009.
- Baker & Baker, "Gravitational wave generator apparatus," Physics Procedia 38, 288–297 (2012).

### Inferences
**Our checks of the Li-group arithmetic**
- **ψ0 for 10 W and W0 = 5 cm.** The standard GB relation P = πW0²ψ0²/(4Z0) gives ψ0 = 1.385×10³ V/m against the published 1.26×10³ (ratio 1.10). The paper does not state its convention.
- **Graviton flux through the waist.** (c³ω²A²/8πG)·πW0²/ħω = 3.76×10¹⁶ s⁻¹ (published 3.77×10¹⁶, ratio 0.997).
- **PRD 2009 eqs. (7)–(23) all reproduce within 0–4%** using the paper's own k = 100 m⁻¹ and hν = 3.3×10⁻²⁴ J:
  - U_em^(2) = 2.39×10⁻⁴⁰ W
  - N^(2) = 7.2×10⁻¹⁷ s⁻¹
  - U_gw = 1.585×10⁻⁷ W
  - N_g = 4.79×10¹⁶ s⁻¹
  - conversion 1.5×10⁻³³
  - E0 = 868 V/m
  - U1 = 6.91×10⁻²⁰ W
  - N1 = 2.09×10⁴ s⁻¹
  - Δt(h = 10⁻²⁶) = 6.9×10⁷ s
  - Δt_min(x = 0) = 3.1×10³ s and 3.1×10⁵ s
- **The published "signal" at x = 3.5 cm is the homodyne interference term.** Take the pure inverse-G photon rate for the EPJC parameters (3 T, 6 m, 0.01 m², 5 GHz, h = 10⁻³⁰, Li's field formula): Ṅ_GW = 2.56×10⁻¹⁸ s⁻¹. Then 2√(N_BPF(3.5 cm)·Ṅ_GW) = 2√(1.24×10²² × 2.56×10⁻¹⁸) = 3.57×10² s⁻¹, against Table I's 3.54×10² s⁻¹ (ratio 1.007).
  - This is exactly JASON's homodyne cross term 2(Ṅ0Ṅ_GW)^½ (eq. 3-32), with Ṅ0 equal to the local background flux.
  - Li et al. themselves call N_x^(1) "equivalent to" 2(Ṅ0Ṅ_GW)^½ in the plane-wave case (PRD 2009 §4-3).
- **The Li formulas use a GW energy flux of c³ω²A²/(8πG).** That is 4× the standard single-polarization flux c³ω²A²/(32πG) (2× the circular-polarization value). See Q5.
- **Rust test vectors** (published value, then our recomputation):
  - ψ0: 1.26e3 (ours 1.385e3)
  - N_g (waist): 3.77e16 (ours 3.758e16)
  - U_em^(2): 2.3e-40 W (ours 2.386e-40)
  - N_γ^(2): 7.0e-17 s⁻¹ (ours 7.229e-17)
  - conversion: 1.4e-33 (ours 1.505e-33)
  - E0: 8.7e2 V/m (ours 8.680e2)
  - U_z^(1): 6.9e-20 W (ours 6.907e-20)
  - N_z^(1): 2.1e4 (ours 2.093e4)
  - Δt (h = 1e-26): 6.8e7 s (ours 6.916e7)
  - N_x^(1)(3.5 cm), via 2√(N0·N_GW): 3.54e2 (ours 3.565e2)
  - Δt_min: 3.0e3 s (ours 3.123e3)

### Gaps
- Table I's N_x^(1)(x=0) = 8.21×10² s⁻¹ was not re-integrated from eqs. (59)–(61). That needs a numerical 2-D integral of the full eq. (59), which is printed with garbled layout in the arXiv text. The scaling ∝ h (linear) is exact by construction.
- The DOI of EPJC 56, 407 (2008) was not verified (the arXiv id is given instead). The 2025 EPJC "Laguerre–Gaussian beam HFGW source" paper (DOI 10.1140/epjc/s10052-025-14158-1) surfaced in search but was not read.
- I could not retrieve Woods et al., J. Mod. Phys. (2011), cited by Baker & Baker as giving sensitivity 10⁻³⁰–10⁻³² with SNR > 1500 on p. 511.

## Q2. What did JASON (JSR-08-506, 2008) conclude, what is wrong and by how much, and what later analyses exist?

### Takeaway
JASON (lead author D. Eardley, for ODNI) found that the Li–Baker sensitivity is too optimistic by more than 20 orders of magnitude. The error is not gravitation physics: the Li–Baker analysis treats the homodyne cross term 2E0·E_GW as a signal without counting the shot noise of the reference beam that must be present for that term to exist. Li et al. answered only implicitly (PRD 2009): they place the receiver at x = 0, where the transverse background Poynting flux vanishes. The 2025 Living Review records the concerns and notes no experimental results. Status: Li–Baker sensitivity claims **[refuted — JASON 2008; concerns noted by Aggarwal et al. 2025]**.

### Cited Findings
- **Report identity.** JSR-08-506, October 2008, "D. Eardley et al.", The MITRE Corporation / JASON, sponsored by the ODNI and DIA. Approved for public release. The study was requested by the National MASINT Committee staff of ODNI.
  - JASON received briefings on June 17–18, 2008 from Robert Baker, Gary Stephenson, Paul Murad, Patricia Walters, Ronald Pandolfi, Kevin Pollpeter and Mark Pesses.
  - It also collected about 1000 pages of material. A classified appendix exists.
  - Source: [JASON JSR-08-506 (archive.org copy of FAS PDF)](https://web.archive.org/web/2010/http://www.fas.org/irp/agency/dod/jason/gravwaves.pdf); original [FAS URL](https://irp.fas.org/agency/dod/jason/gravwaves.pdf)
- **Conclusions** (Executive Summary and §6). Source: [JASON](https://web.archive.org/web/2010/http://www.fas.org/irp/agency/dod/jason/gravwaves.pdf)
  - (1) The previous Li–Baker analysis is incorrect by many orders of magnitude.
  - (2) Detecting natural relic HFGW and artificial HFGW sources is infeasible in the foreseeable future.
  - (3) No foreign HFGW threat is credible: not communications, not radar or tomography, not propulsion, not any other use.
  - (4) In weak fields GR and its experimental basis are complete and reliable.
  - §6 says the homodyne-sensitivity estimates are "too optimistic by more than 20 orders of magnitude", and that this is an error of electrical engineering, not of gravitation physics.
- **JASON's generator physics (§3.3, Gaussian units).** Source: [JASON](https://web.archive.org/web/2010/http://www.fas.org/irp/agency/dod/jason/gravwaves.pdf)
  - (3-13) P_GW(out) = F P_EMW(in)
  - (3-14) F = 4πG B0² L²/c⁴ ≈ 10⁻³⁵ for B0 = 10⁵ G and L² = 10³ cm²
  - (3-15) P_GW = 4πG B0² U L/c³
  - (3-16) U = P_EMW L/c; P = 10 kW and L = 30 cm give U = 10⁻⁵ J
  - (3-17)–(3-19) Copper cavity with Ḣ ≈ 10⁶ W and Q ≈ 2×10³: U_max ≈ 0.2 J, giving max P_GW ≈ 2×10⁻²⁷ W ≈ 2×10⁻² graviton/s
  - (3-20)–(3-22) Superconducting cavity with field-emission limit E⊥ < 50 MV/m, u < 10⁷ erg cm⁻³, V ≈ 3×10³ cm³: P_GW ≈ 10⁻²³ W ≈ 10² graviton/s
  - (3-23)–(3-24) Mechanical limit u ≈ 10¹⁰ dyn/cm²: 10⁻²⁰ W ≈ 10⁵ graviton/s
  - (3-25) Container-explosion limit: 10⁻¹⁸ W ≈ 10⁷ graviton/s
  - (3-26)–(3-28) Target flux P_GW·b/(4πd²) with beaming factor b ≈ 10², below 30 gravitons cm⁻² s⁻¹ even for V = 10⁷ cm³
- **JASON's detector physics (§3.4).** Source: [JASON](https://web.archive.org/web/2010/http://www.fas.org/irp/agency/dod/jason/gravwaves.pdf)
  - (3-29) P_EM(out) = F P_GW(in), with F ≈ 10⁻³⁵. The best generator (10² graviton/s), b ≈ 10², d ≈ 10 m and a 10⁴ cm² detector give about 10⁻³⁴ photons/s.
  - (3-30) t̂_γ = 10³⁴ s.
  - Homodyne analysis:
    - (3-31) Ṅ_d ∝ |E_T|² = |E0|² + 2E0·E_GW + |E_GW|²
    - (3-32) N_d = Ṅ0t + 2(Ṅ0Ṅ_GW)^½ t cos δ + Ṅ_GW t
    - (3-33) (Ṅ0Ṅ_GW)^½ t > (Ṅ0t)^½
    - (3-34) t > 1/Ṅ_GW = t̂_γ, which is independent of Ṅ0. Heterodyning gives the same result within a factor 2 (Marcuse).
  - Cavity:
    - (3-36) Ṅ_γ ≈ Q/t̂_γ
    - (3-38) t1 ≈ (L/c · t̂_γ)^½ ≈ 10¹² s
    - (3-40) t̂_γ/Q ≈ 10³¹ s for copper (Q ≈ 2×10³)
    - (3-41) ≈ 10²³ s for superconducting walls (Q ≈ 10¹¹)
  - (3-11) Radar back-scatter R ≈ (4πGΔρ/ω²)² ≈ 10⁻⁵² for Δρ ≈ 1 g/cm³ and ω ≈ 10¹⁰ s⁻¹.
- **§4.1, the specific Li–Baker error.** Source: [JASON](https://web.archive.org/web/2010/http://www.fas.org/irp/agency/dod/jason/gravwaves.pdf)
  - The proposal claims about 10²⁵ times the sensitivity of (3-30)/(3-34).
  - That factor comes from a reference beam |E0| = 1.17×10³ V/m added to E_GW ≈ ĥB0ωL/c = 10⁻²² V/m (eq. 4-42, for ĥ = 3×10⁻³², B0 = 10⁵ G and ωL/c = 10). The cross term then exceeds |E_GW|² by 10²⁵.
  - The proposal neglects the |E0|² photon noise by claiming E0² never reaches the detector. JASON calls this "physically untenable": the interference term cannot be detected unless E0 is present at the same time.
  - A fractal membrane cannot separate exactly parallel E0 and E_GW beams.
- **§4.2.** Source: [JASON](https://web.archive.org/web/2010/http://www.fas.org/irp/agency/dod/jason/gravwaves.pdf)
  - Model the relic background as a 2 K graviton bath: about 10 nW lands on 10 cm × 10 cm.
  - The graviton-to-photon conversion coefficient is about 10⁻³⁸ (9 T), giving about 10⁻⁴⁶ W, or one photon per 10²⁴ s.
  - Thermal emission would require cooling to about 20 mK.
- **§4.3, fractal membrane.**
  - Stephenson's briefing claimed a Q contribution of 3.4×10²¹ and an antenna gain of 6.3×10⁴. JASON: essentially none of this is valid.
  - The membrane cannot be superconducting in 9 T, cannot reduce phase space, and is a linear lens bound by the diffraction limit.
  - Source: [JASON](https://web.archive.org/web/2010/http://www.fas.org/irp/agency/dod/jason/gravwaves.pdf)
- **§4.4, Stephenson's stored energy.**
  - Stephenson's assumption: 10¹³ J stored, from 10⁵ W × 10³ s × Q 10⁵.
  - JASON: 10¹³ J is 2.4 kt and would need 10⁵ MBar of containment.
  - A 10⁵ W source with Q = 10⁵ at 5 GHz holds only about 2 J (residence time about 2×10⁻⁵ s).
  - The implied field (about 1.5 GV/m) is about 50× the best accelerator RF fields (30 MV/m).
  - JASON also calls the gas-pressure calculation in reference [10] "grossly the wrong calculation," and says superconductors have no enhanced GW interaction.
  - Source: [JASON](https://web.archive.org/web/2010/http://www.fas.org/irp/agency/dod/jason/gravwaves.pdf)
- **§5, applications.**
  - 1 g propulsion by a GW focus 100 m ahead needs about 10³⁸ W, more than 10²⁵ times world power consumption.
  - Earth imaging at 10⁻⁴⁰ rad phase needs about 10⁸⁰ gravitons.
  - Baker's "$50 B timing benefit" is called nonsense.
  - Source: [JASON](https://web.archive.org/web/2010/http://www.fas.org/irp/agency/dod/jason/gravwaves.pdf)
- **Baker's rebuttal (GravWave "Q&A").**
  - It says JASON found GravWave wrong "by an order of magnitude". That is contradicted by JASON's own text ("many orders of magnitude", ">20 orders").
  - It argues the detector rests on the "Li-effect" (synchro-resonance plus Gertsenshtein), not the Gertsenshtein effect, and that JASON did not consult the presenters.
  - It cites Bekenstein (1973) GW recoil of 100–300 km/s as evidence for GW propulsion.
  - Source: [Baker, Q&A: JASON Report](https://www.drrobertbaker.com/docs/Q%20&%20A%20JASON.pdf)
- **Li et al.'s implicit answer (PRD 2009).** They accept JASON's plane-wave homodyne result (their eq. 23 matches JASON 3-33). They argue that in the GB geometry the PPF and BPF have different directions and decay rates, so at x = 0 the noise is only N_x(other) (scattering, diffraction, thermal), giving Δt_min = N_other/N1² (eqs. 50–58). — [arXiv:0909.4118](https://arxiv.org/abs/0909.4118)
- **Living Review, Aggarwal et al. (2025), arXiv:2501.11723, §5.5.4 "Heterodyne Enhancement of Magnetic Conversion".**
  - Refs. [538–541] (Li & Yang 2004; Li, Baker, Chen 2006; Baker, Stephenson, Li AIP 969 (2008); Li et al. PRD 80 064013) proposed seeding the conversion volume with a same-frequency EM wave.
  - Concerns were raised in [542] (Woods 2012) and [543] (Eardley 2008, i.e. JASON).
  - The claimed sensitivities rely on technology that does not exist, and no experimental results suggest feasibility.
  - The review gives the schematic conversion probability P_{h→γ} ∼ (BL)²/M_Pl² (eq. 128).
  - Source: [arXiv:2501.11723](https://arxiv.org/abs/2501.11723)
- **Eardley's role.** He is the named lead author ("D. Eardley et al.") on the report documentation page, and the Living Review cites the report as "D. Eardley, High frequency gravitational waves, 2008." — [JASON](https://web.archive.org/web/2010/http://www.fas.org/irp/agency/dod/jason/gravwaves.pdf); [arXiv:2501.11723](https://arxiv.org/abs/2501.11723)
- **Stephenson.** Gary V. Stephenson (Seculine Consulting) is a co-author of EPJC 2008, PRD 2009 and AIP 969. He gave JASON the briefing "The Standard Quantum Limit for the Li-Baker HFGW Detector" (JASON ref. [15]). Baker & Baker (2012) cite Stephenson for an improved standard quantum limit via superconducting resonance chambers. — [JASON](https://web.archive.org/web/2010/http://www.fas.org/irp/agency/dod/jason/gravwaves.pdf); [Baker & Baker 2012](https://inspirehep.net/files/6bf245c473964d406575a245362a0ede)

### Inferences
- **The exact point of disagreement.** Both sides agree on the classical fields and on the cross-term algebra. They disagree on the noise floor.
  - JASON's floor is the shot noise √(Ṅ0 t) of every reference-beam photon that carries the interference term. Then SNR = 2√(Ṅ_GW t), so t ≳ 1/Ṅ_GW whatever Ṅ0 is. This is the standard quantum limit for homodyne detection, a textbook result: a linear phase-insensitive or homodyne measurement cannot beat the photon number of the signal itself.
  - Li et al. take the floor as the transverse Poynting-flux component n_x^(0) crossing a chosen surface. That is zero at x = 0 by symmetry, although the field energy, including B̃_z^(0), which appears in their own cross term, is present there.
  - Our check above shows the published x = 3.5 cm "signal" equals 2√(N_BPF Ṅ_GW). That supports JASON's reading that the PPF is the homodyne term.
- **Size of the discrepancy for the EPJC/PRD parameters** (3 T, 6 m, 0.01 m², 5 GHz, h = 10⁻³⁰; our computation):
  - Ṅ_GW = 2.6×10⁻¹⁸ s⁻¹ (Li's field formula) or 3.2×10⁻¹⁹ s⁻¹ (standard ½hBωL field). JASON's floor is t ≥ 3.9×10¹⁷ to 3.1×10¹⁸ s.
  - Li et al. claim Δt_min = 3×10³ s: a ratio of about 1.3×10¹⁴ to 10¹⁵ in time or photon number.
  - At fixed t = 3×10³ s the JASON-limited strain is h ≥ 1.1×10⁻²³ (Li's field formula), about 7 orders of magnitude worse than 10⁻³⁰.
  - JASON's own "10²⁵" is the power ratio 2E0/E_GW for its parameter set. Our check gives 2.3×10²⁵.
- **JASON arithmetic slips** (our checks). None changes the conclusion, and all err in the proponents' favour:
  - E_GW for its stated parameters is 9×10⁻²² V/m, not 10⁻²² V/m.
  - u = E²/8π at 50 MV/m is 1.1×10⁵ erg cm⁻³, about 100× below the 10⁷ used.
  - Eq. (3-19) gives 1.9×10⁻²⁷ W, which is 1.8×10⁻³ graviton/s at ħω (ω = 10¹⁰ s⁻¹), not 2×10⁻² graviton/s.
  - The §4.4 stored energy is P·Q/ω = 0.32 J (JASON: about 2 J).
  - The §4.4 "E = 5×10⁵ esu/cm or 150 MV/cm or 1.5 GV/m" chain mixes units: 5×10⁵ statV/cm is 1.5×10¹⁰ V/m.
  - JASON's Gaussian F = 4πGB²L²/c⁴ is 4π higher than the first-order result F = GB²L²/c⁴ (see Q5).
- **Rust test vectors (JASON):**
  - F(10⁵ G, 10³ cm²) = 1.038e-35
  - U(10 kW, 30 cm) = 1.0e-5 J
  - U_max(1 MW, Q = 2e3, ω = 1e10) = 0.2 J
  - P_GW(eq. 3-19) = 1.87e-27 W
  - P_GW(eq. 3-24) = 2.8e-20 W
  - 2 K bath over 0.01 m² = 9.1e-9 W
  - wait time at 10⁻⁴⁶ W = 1.1e24 s
  - 10¹³ J = 2.39 kt

### Gaps
- The content of Woods (2012, Physics Procedia 38, 66–76), cited by the Living Review as a "concern", could not be retrieved. Note that Woods is also a Baker co-author on the generator papers.
- No standalone Stephenson SQL paper was found; only the JASON summary of his briefing.
- No peer-reviewed rebuttal of JASON by Li or Baker was found. PRD 2009 does not cite JASON.
- I did not fetch the FAS "Secrecy News" coverage (Dec 2008), which is behind a bot challenge.
- I did not find a published "independent analysis" beyond JASON, the Living Review paragraph and Woods (2012).

## Q3. Baker's HFGW generator power estimates (FBAR arrays, lasers, IR rings) against the standard quadrupole formula

### Takeaway
Baker's "jerk" formula P = 1.76×10⁻⁵²(2rΔf/Δt)² W uses the correct constant 32G/5c⁵. For a single FBAR pair it agrees with the standard quadrupole result to within a factor of about 2 (both about 10⁻³³ W). The large published numbers come from three multipliers:
- coherent N² summation over arrays far larger than the GW wavelength;
- a "superradiant" needle beam about 10¹¹ times narrower in solid angle than diffraction allows;
- a GW flux-to-strain relation 2× off.

Status: the single-element formula is **[mainstream]**. The N² array power and needle-beam fluxes are **[refuted/contested]**: JASON rejects lab HFGW applications generally, and the diffraction-limit contradiction is our inference.

### Cited Findings
- **Core formulas.**
  - Laser paper (STAIF-2006, AIP 813, eq. 1): "P = 2G(2δm)²r⁴ω⁶/45c⁵ = 1.76×10⁻⁵² (2rΔf/Δt)² W", presented as the quadrupole approximation P = (G/45c⁵)(d³I/dt³)².
  - Eq. 2: A = (8πG F_GW/c³ω²)^½ = 1.28×10⁻¹⁸ F_GW^½/ν_GW.
  - Source: [Baker, Li, Li 2006](https://www.drrobertbaker.com/docs/AIP;%20HFGW%20Laser%20Generator.pdf)
- **Same formulas elsewhere.** Physics Procedia 2012 eqs. (1) and (4) repeat them. Baker & Black (SPESIF-2009) eq. (9) gives P(r, Δf, ν) = 1.76×10⁻⁵²(2r ν_GW Δf)² W, attributed to Baker, Woods, Li (2006) eq. (2). — [Baker & Baker 2012](https://inspirehep.net/files/6bf245c473964d406575a245362a0ede); [Baker & Black 2009](https://www.drrobertbaker.com/docs/Analyses%20of%20HFGW%20Generators%20and%20Radiation%20Pattern.pdf)
- **Radiation pattern (Baker & Black).**
  - Single pair: half-angle 47.5°, cap 0.1358 m² on the unit-area sphere, flux at 1 m F(1.0) = 0.01146 P_i W m⁻².
  - N pairs: A_cap = A½(N=1)/N and F(1.0) = N²·0.01146 P_i, i.e. an "N² law".
  - Source: [Baker & Black 2009](https://www.drrobertbaker.com/docs/Analyses%20of%20HFGW%20Generators%20and%20Radiation%20Pattern.pdf)
- **FBAR proof-of-concept** (Baker, SPESIF-2009). Source: [Baker 2009](https://www.drrobertbaker.com/docs/HFGW%27s%20Power%20Requirements.pdf)
  - Hardware: 1.8×10⁸ cell-phone FBARs (100 µm with 10 µm margins) driven by 10,000 magnetrons, 19.8 km long.
  - Output: total HFGW power 0.066 W; flux 3.53 W m⁻² at one wavelength (6.1 cm) beyond the last element; A = 4.9×10⁻²⁸.
  - Staggered and sliced variants: 198 m, then about 198 cm, then about 2 cm (A = 4.9×10⁻²⁶ at 1 kW per band), or 20 cm with A = 4×10⁻²⁷.
  - Power: 20,000 magnetrons need about 25 MW. A nanotech version needs less than about 20 W.
  - Baker states the amplitude is ∝ N (not N²) once elements are sliced.
- **Woods IR ring** (as summarised by Baker 2009). Source: [Baker 2009](https://www.drrobertbaker.com/docs/HFGW%27s%20Power%20Requirements.pdf)
  - Pentane at 6.3×10²⁸ m⁻³ in a 100 m-radius monomode ring gives 2N = 4.85×10¹⁷ jerkable masses, with P_i = 4.62×10⁻¹⁶ W, and a ring flux of 1.29 W m⁻².
  - With 1 m rings ×100 per plate and 10⁶ stacked plates (n² law): 1.29×10⁸ W m⁻² and A = 1.21×10⁻²⁸.
- **Double-helix FBAR array** (Baker & Baker, Physics Procedia 38, 288–297, 2012). Source: [Baker & Baker 2012](https://inspirehep.net/files/6bf245c473964d406575a245362a0ede)
  - Geometry: 20 m long; r1 = 0.35 m, r2 = 0.05 m; volume 7.5 m³; FBARs 110×110×20 µm; N = 1.55×10¹³ pairs.
  - Drive: Δf = 2 N per FBAR, Δt = 4×10⁻¹⁰ s (2.5 GHz drive, GW at 5 GHz, λ_GW = 6 cm), 2r_rms = 0.5 m.
  - Coherent sum ("standard approach"): P = 1.76×10⁻⁵²(0.5 × 3.1×10¹³/4×10⁻¹⁰)² = 2.64×10⁻⁷ W.
  - Needle cap A_cap = 1.71/N m² at 1 m, so S(1 m) = (P/4)/(1.71/N) = 6×10⁵ W m⁻² and A = 1.98×10⁻²⁵.
  - At 1.3×10⁷ m: S = 3.55×10⁻⁹ W m⁻², A = 1.52×10⁻³².
  - "Conservative" (incoherent) case: NP = 1.7×10⁻²⁰ W, S(1 m) = 3.8×10⁻⁸ W m⁻², A = 5×10⁻³²; at Earth's diameter S = 1.33×10⁻²⁰ W m⁻², A = 3.8×10⁻³⁹.
  - Energizing: 1000 × 1 kW magnetrons, 3.2×10⁴ W m⁻², described as "32 gigawatts per square meter" of superradiant microwaves.
- **Laser generator** (STAIF-2006). Source: [Baker, Li, Li 2006](https://www.drrobertbaker.com/docs/AIP;%20HFGW%20Laser%20Generator.pdf)
  - SIOM pulses of 23 TW, 33.9 fs, 800 nm: 3.14×10¹⁸ photons per pulse, giving an impulsive force Δf = (1+R)(h/λ)(photons/s) = 1.5×10⁵ N.
  - P ≈ 1.76×10⁻⁵²(2 targets × 2r × Δf/Δt)² = 5.5×10⁻¹² W.
  - Peak flux 1.23×10⁻⁷ W m⁻², A = 4.5×10⁻²³ at 3 mm.
  - Detector: 2.5 mm-radius pulsed GB in 15 T, ψ0 = 1.8×10¹² V/m for 10¹⁴ W, claimed sensitive to about 10⁻³². Space-borne at lunar distance: about 2×10³ W and 10¹³–10¹⁴ W m⁻² at the focal spot.
- **Calibration claim (PSR 1913+16).** Separation 4.1×10⁹ m, Δf = 1.16×10³³ N and 7.2×10⁻⁵ Hz give 10.1×10²⁴ W, which Baker calls "very close" to the relativistic value. The flux at Earth (9.5×10¹⁵ m) is 2.26×10⁻⁸ W m⁻². — [Baker 2009](https://www.drrobertbaker.com/docs/HFGW%27s%20Power%20Requirements.pdf)
- **JASON's generic generator ceilings** apply to EM-driven sources (Gertsenshtein route): 2×10⁻²⁷ W (Cu cavity) to 10⁻¹⁸ W (explosion limit). Baker's rebuttal says no HFGW generator proposal uses the Gertsenshtein effect. — [JASON](https://web.archive.org/web/2010/http://www.fas.org/irp/agency/dod/jason/gravwaves.pdf); [Baker Q&A](https://www.drrobertbaker.com/docs/Q%20&%20A%20JASON.pdf)

### Inferences
**Our checks of Baker's formulas and numbers**
- **The constant.** 1.76×10⁻⁵² = 32G/(5c⁵) = 1.764×10⁻⁵² SI. For a circular binary, |df/dt| = m r ω³, so 1.76×10⁻⁵²(2rΔf/Δt)² equals the standard Peters result (32G/5c⁵)μ²a⁴ω⁶ **[mainstream]**.
- **The laser eq. (1) is internally inconsistent by a factor of 144.**
  - The left side 2G(2δm)²r⁴ω⁶/45c⁵ is the rotating-rod formula with rod length r. For a dumbbell it is (8/45)/(128/5) = 1/144 of the right side.
  - The two sides agree only if Δf/Δt = δm r ω³/12.
- **Single FBAR pair, standard quadrupole for linear oscillation.**
  - Two masses at ±r, each driven by force F0 sin ωt: I⃛_xx ≈ 4rḟ and Σ Q⃛ij² = (2/3) I⃛_xx², so P = (G/5c⁵)(2/3)(4r)²⟨ḟ²⟩ = (16G/15c⁵) r² F0² ω².
  - For r = 0.25 m, F0 = 2 N and ω = 2π×2.5 GHz: 1.8×10⁻³³ W. Baker's form gives 1.1×10⁻³³ W (ratio 0.61), so close.
- **Array totals.**
  - Fully coherent N² (valid only if all sources sit within ≲ λ/2 or are perfectly phased into one beam): standard 4.4×10⁻⁷ W against Baker's 2.64×10⁻⁷ W.
  - Incoherent N: 2.8×10⁻²⁰ W against Baker's "conservative" 1.7×10⁻²⁰ W.
  - The array (0.7 m × 20 m) is far larger than λ_GW = 6 cm, so the N² figure is an upper bound that a phased end-fire array cannot reach in total power.
- **The "superradiant" needle beam breaks the diffraction limit.**
  - A_cap = 1.71/N m² at 1 m implies a half-angle of 1.9×10⁻⁷ rad. Any wave obeying the free wave equation (JASON eq. 3-9) from an aperture D ≈ 0.7 m at λ = 6 cm spreads at least about 1.22λ/D ≈ 0.10 rad.
  - So the solid angle, and hence the flux, is overstated by about 3×10¹¹.
  - Even granting the coherent 4.4×10⁻⁷ W, the diffraction-limited flux at 1 m is about 1.3×10⁻⁵ W m⁻², or A ≈ 1.8×10⁻³⁰ (standard flux relation). Incoherent and isotropic it is 2.2×10⁻²¹ W m⁻², or A ≈ 2.4×10⁻³⁸.
  - Detecting h ≈ 10⁻³⁰ at 5 GHz is itself beyond the JASON limit (Q2).
- **Flux-to-strain relation.** Baker's A = (8πG F/c³ω²)^½ uses F = c³ω²A²/(8πG), 4× the standard single-polarization flux. For a given flux the standard A is 2× Baker's: the FBAR case at 3.53 W m⁻² gives A = 9.6×10⁻²⁸, not 4.8×10⁻²⁸.
  - The published values reproduce from Baker's relation: our 4.80×10⁻²⁸ (published 4.9×10⁻²⁸), 1.50×10⁻³² (1.52×10⁻³²), coefficient 1.256×10⁻¹⁸ (1.28×10⁻¹⁸).
- **Laser force.** Δf = 1.50×10⁵ N reproduces (ratio 0.997).
- **PSR 1913+16.** Peters–Mathews with m1 = 1.4398, m2 = 1.3886 M☉, a = 1.95×10⁹ m and e = 0.6171 gives 7.7×10²⁴ W. Baker's jerk form with r = 2.05×10⁹ m gives 2.1×10²⁵ W; we could not reproduce his 1.01×10²⁵ W exactly.
- **Rust test vectors:**
  - P_pair(Baker) = 1.100e-33 W
  - P_pair(std) = 1.813e-33 W
  - P_total(Baker coherent) = 2.643e-7 W
  - P_total(std coherent) = 4.357e-7 W
  - P_total(std incoherent) = 2.811e-20 W
  - A(3.53 W m⁻², 4.915 GHz) = 4.80e-28 (Baker) or 9.60e-28 (standard)

### Gaps
- I could not reproduce the 0.066 W total for the 1.8×10⁸-FBAR array from the documents read. It depends on Baker, Stephenson, Li (AIP 969, 2008) eq. (8), whose gravwave.com URL now returns no PDF.
- Woods & Baker (2005, STAIF) and Woods & Baker (2009, IR rings) were cited but not read.
- No source was found for a Baker "HFGW from rocket exhaust" estimate. I found no evidence that such a paper exists.
- No independent peer-reviewed analysis of the FBAR or double-helix generator was found. Our critique above is our own inference.

## Q4. The Pais patents: equations, numbers, Navy actions and tests, and critiques

### Takeaway
Pais's patents rest on one formula, S_max = f_G(σ²/ε0)×(velocity-like term). It reproduces his headline numbers (10²⁴–10³³ W m⁻², 10¹⁰ W of HFGW, 10³⁵ gravitons/s) only with surface charges 10²–10⁷ times above the vacuum breakdown limit that his own patent states. NAWCAD spent $508k (Oct 2016–Sep 2019) on the "High Energy Electromagnetic Field Generator". The spin test found no anomalous effect, and the Navy says the "Pais Effect" "could not be proven". No test of the HFGW generator itself was reported. Status: **[refuted]** in the sense of no supporting evidence and a null Navy test, with critiques from physicists (Willis, Webb) and the USPTO examiner's initial rejection.

### Cited Findings
- **US 10,322,827 B2 "High frequency gravitational wave generator"** (inventor Pais, assignee US Navy; Google Patents lists it as expired, fee related). Source for all items below: [US10322827B2](https://patents.google.com/patent/US10322827B2/en)
  - **Device:** a gas-filled shell with an electrically charged outer surface vibrated by microwave emitters, plus sound generators and acoustic resonant cavities. Mössbauer Co-57/Fe-57 doping is said to "further amplify" the EM flux. Microwave pulsing is at 300 MHz–300 GHz.
  - **Eq. 1:** S_max = f_G(σ²/ε0)[(R_v ν²) t_op], where f_G = 1 for a disc, σ is surface charge density, R_v is vibration amplitude, ν is "angular frequency of vibration in Hertz" and t_op is the operating time at maximum acceleration.
  - **Eq. 2:** S_max = f_G(σ²/ε0)[(R_v ν³) t_op²].
  - **Eq. 3:** E_T = m R_v² Ω² [exp(2Ωt)], with negative damping b/2m ≫ Ω0.
  - **Eq. 4:** P_GW(out) = [(4πG/c⁴)B0²L²] P_EM(in), taken from JASON (3-13)/(3-14).
  - **Claims:**
    - With all parameters other than ν "of unit order", S ≈ 10³³ W m⁻² at about 10⁷ Hz.
    - JASON's eq. 4 is "missing" a factor 1/μ0 ≈ 10⁶.
    - About 10¹⁰ W of HFGW for an input of 10³³ W m⁻² with B0 ≈ 10¹⁰ T. Via JASON eq. (3-25) this is said to equal 10³⁵ gravitons/s.
    - Focused HFGWs could "induce a curvature singularity".
    - JASON's conclusions "must be revisited, are incorrect".
- **US 10,144,532 B2 "Craft using an inertial mass reduction device"** (application 15/141,270, filed 2016-04-28). Source for all items below: [US10144532B2](https://patents.google.com/patent/US10144532B2/en)
  - **Eq. 1:** S_max = f_G(σ²/ε0)[R_r ω + R_v ν + v_R].
    - Example: σ = 50,000 C/m², a 2 m disc at 30,000 RPM gives about 10²⁴ W m⁻².
    - Adding 10⁹–10¹⁸ Hz vibration gives 10²⁴–10²⁸ W m⁻².
  - **Eq. 2:** with 10¹¹ Hz, 1 C/m² and t_op about the inverse of the vibration amplitude, the patent claims 10³³ W m⁻².
  - **Eq. 3:** 2γ(dm/dt)_pp c² = S_max A_S (pair-production mass rate).
  - **Eq. 4:** (S_max/c) = ħν_v⁴/(8π²c³). Setting the left side to ε0E² at the Schwinger field (about 10¹⁸ V/m) gives ν_v of order 10²² Hz.
  - **Eq. 5:** Hayasaka–Takeuchi gyroscope weight loss ΔW_R(ω) = −2×10⁻¹⁰ M r_eq ω kg m s⁻².
  - **Eq. 6:** ΔW_R(ν) = −f_ep M A_ν ν.
- **US 10,135,366 B2 (high-energy EM field generator for asteroid deflection).** Source for all items below: [US10135366B2](https://patents.google.com/patent/US10135366B2/en)
  - **Relativistic Eq. 1:** S_max = f_G γ²(σ²/ε0)[R_r ω + R_v ν + v_R], with f_G = 2/3 for a sphere, 1/6 for a hemisphere and 1 for a disc.
  - **Eq. 3, momentum balance:** (M_A + M_KI)Δv = k M_KI u_KI. The kinetic-impactor case gives 0.001 m/s. Replacing the right side with (S_max/c²)Vol_S, using S = 10²⁴ W m⁻² and Vol = 10³ m³, gives 1 m/s for a 10¹⁰ kg asteroid.
  - **Stated breakdown limits:** 3×10⁶ V/m in dry air, 10⁹ V/m in vacuum.
  - **Power:** a 1–10 MW pebble-bed reactor.
- **US 2019/0295733 A1 "Plasma compression fusion device".** Source: [US20190295733A1](https://patents.google.com/patent/US20190295733A1/en)
  - Counter-spinning, electrically charged conical "dynamic fusors" whose "concentrated magnetic energy flux" compresses and heats the plasma.
  - Lawson criterion (eq. 1, pressure × confinement time). Device size 0.3–2 m.
  - A secondary article reports claimed output of gigawatts to terawatts: [Interesting Engineering](https://interestingengineering.com/science/nikola-tesla-engineer-ufo-patents-us-navy)
- **US 2019/0058105 A1 "Piezoelectricity-induced room temperature superconductor" (the superconducting EM field generator).** Source: [US20190058105A1](https://patents.google.com/patent/US20190058105A1/en)
  - Eq. 1: London depth λ_L = [m_s/(μ0 n_s q_s²)]^½.
  - Eq. 2: the Meissner condition is written as [(μ0I/2πR) + (μ0 σ A_v ω_v² Δt)] ≥ (B_E)_MAX, claimed independent of T_c.
  - Eq. 3: α = e/(φ0 ε0 c).
  - Drive: vibration above 100 MHz.
- **USPTO prosecution of US 10,144,532.**
  - Examiner Philip Bonzell's rejection (dated in the article as March 30, 2018) said there is no such thing as a "repulsive EM energy field". He noted the scheme would need about 10⁹ T and 10¹⁸ V/m.
  - Dr. James Sheehy (CTO, Naval Aviation Enterprise) wrote that "China is already investing significantly in this area".
  - A notice of allowance followed on Oct. 31, 2018.
  - Source: [TWZ, Tingley & Rogoway, 2019](https://www.twz.com/28729/docs-show-navy-got-ufo-patent-granted-by-warning-of-similar-chinese-tech-advances)
- **NAWCAD official statement** (Timothy Boulay, NAWCAD Communications Director, to TWZ). Source: [TWZ, Tingley, updated Feb 1 2021](https://www.twz.com/39012/the-navy-finally-speaks-up-about-its-bizarre-ufo-patent-experiments)
  - Testing of the High Energy Electromagnetic Field Generator ran from Oct 2016 to Sep 2019 and cost $508,000 over three years, of which about $462k was salaries.
  - When testing ended in Sep 2019 the "Pais Effect" could not be proven.
  - No further research has been conducted, and the project has not transitioned anywhere.
  - Pais left NAWCAD in June 2019 for Strategic Systems Programs, then moved to the USAF.
- **FOIA test details.** Source: [TWZ, Tingley, updated Jan 26 2021](https://www.twz.com/38937/navy-ufo-patent-documents-talk-of-spacetime-modification-weapon-detail-experimental-testing)
  - Funding and effort: at least $466,810 over FY2017–2019, at least 10 technicians and engineers, about 1,600 hours (2018–2019).
  - Test articles: a coin-cell capacitor about 0.276 in diameter, spun to 100,000 RPM by a 2-inch air motor, with no means of vibrating the spindle. Also a 12-inch disc with cruciform piezoelectric elements.
  - The September 2018 spin test "did not measure any anomalous electromagnetic effects". The Aug 2019 NAVAIR STAIRS report attributes this to charge sitting at the capacitor's centre, and calls the effect "neither observed nor disproved".
  - An FOUO slide describes a "Spacetime Modification Weapon".
- **Internal review.** Source: [TWZ, Tingley, updated Dec 17 2020](https://www.twz.com/37134/emails-show-navys-ufo-patents-went-through-significant-internal-review-resulted-in-a-demo)
  - NAWCAD Technology Transfer Office reviews (2015), Invention Evaluation Board (2016–2017), and a USPTO interview (June 6, 2019) on the superconductor patent.
  - The editors of *Joule* required "compelling experimental validation" before considering Pais's superconductor paper.
- **Expert critiques** (qualitative only).
  - Carl Willis (UNM): "a classic case of pathological science".
  - Stephen Webb (Portsmouth): it is "puzzling" that the patents were awarded, and they do not make sense in terms of physics.
  - Charles Collett (Muhlenberg): possibly not outlandish in theory, but faces "significant engineering challenges".
  - Source: [Interesting Engineering, Ratner, Sep 24 2021](https://interestingengineering.com/science/nikola-tesla-engineer-ufo-patents-us-navy)

### Inferences
**Our checks of the Pais numbers**
- **US 10,144,532 example.** σ²/ε0 · Rω with σ = 5×10⁴ C/m², R = 2 m, ω = 3142 rad/s gives 1.77×10²⁴ W m⁻² (published "order 10²⁴"; reproduces).
  - The implied surface field E = σ/ε0 = 5.6×10¹⁵ V/m is 5.6×10⁶ times the 10⁹ V/m vacuum limit stated in his own US 10,135,366. It is still about 240× below the Schwinger field (1.32×10¹⁸ V/m), so even this example does not reach the "vacuum polarization" regime the patents invoke.
  - With σ capped at 10⁹ V/m (8.9×10⁻³ C/m²), the same formula gives 5.6×10¹⁰ W m⁻², 14 orders lower.
- **US 10,322,827 eq. 2.** With σ = R_v = t_op = f_G = 1 and ν = 10⁷ it gives 1.1×10³² W m⁻² (2.8×10³⁴ if ν is angular). The claimed "10³³" is reproducible to order of magnitude, but σ = 1 C/m² means E = 1.1×10¹¹ V/m, about 100× above breakdown.
- **Eq. 4 as Pais applies it.** (4πG/c⁴)(10¹⁰ T)²(1 m)²(10³³ W) = 1.04×10¹⁰ W, which reproduces "10¹⁰ W". 10¹⁰ W divided by JASON's implied 10⁻²⁵ J per graviton gives 10³⁵ s⁻¹, which also reproduces.
- **On the "missing 1/μ0".** JASON's formula is in Gaussian units (B0 = 10⁵ G, cm), so nothing is missing there. Pais is right that an SI version needs 1/μ0. The first-order SI result 4πGB²L²/(μ0c⁴) applied to his inputs gives 8.3×10¹⁵ W.
  - The inputs are unphysical for a laboratory: B = 10¹⁰ T has energy density B²/2μ0 = 4.0×10²⁵ J/m³ (his "10²⁵ J/m³"), equivalent to about 4×10⁸ kg/m³, with a magnetic pressure of about 4×10²⁵ Pa.
- **US 10,135,366 asteroid.** (10²⁴/c²)·10³/10¹⁰ = 1.11 m/s (reproduces "1 m/s"). Here S/c² is the field momentum density, a correct identity **[mainstream]**. But a bound, steadily spinning charged shell is a static field configuration that radiates nothing, so this momentum is not delivered to the asteroid (our inference).
- **Craft eq. 4.** ε0E_S² = ħν⁴/(8π²c³) gives ν = 4.2×10²¹ Hz (claimed order 10²², ratio 0.42).
  - The right side is the standard zero-point energy density up to an angular cutoff ν.
  - Equating it to ε0E² has no basis in Schwinger pair production, whose rate scales as exp(−πE_S/E).
- **Physics of S_max.** f_G(σ²/ε0)·v is about twice the electrostatic energy density σ²/2ε0 carried along at speed v.
  - With v = a·t_op it grows without bound as operating time increases. It has no c in it and no link to the mechanical input power, which violates energy conservation.
  - Radiated power from accelerated charge follows Larmor/dipole formulas (∝ q²a²/6πε0c³). A uniformly spinning axisymmetric charge distribution is a static field configuration and does not radiate at all.
  - Eq. 3's exponential growth is simply the work done by the negative-damping drive.
- **Superconductor eq. 2.** It is dimensionally consistent: μ0σv is the magnetic field of a moving charge sheet with v = A_vω_v²Δt. But the Meissner effect is a thermodynamic phase property, not an inequality between fields, so "independent of T_c" does not follow.

### Gaps
- No primary NAWCAD test report was read (only TWZ summaries of the FOIA documents). No measured numbers (fields, charges, voltages) were published in the articles.
- I found no test of the HFGW generator patent (US 10,322,827) itself, and no NAVAIR statement specific to it.
- I found no quantitative peer-reviewed physics critique of the Pais patents. Expert critiques in the press are qualitative.
- The exact grant date of US 10,322,827 was not captured; Google Patents was not queried for it.
- The 2015 Pais paper in Int. J. Space Science and Engineering 3(4), 312–317 was not read.

## Q5. What is genuinely calculable and correct, and where exactly do the proposals depart from mainstream physics?

### Takeaway
The inverse-Gertsenshtein pieces are standard physics: an induced field ∝ hBωL, a second-order photon rate ∝ h²B²L², a conversion probability of about 10⁻³³ for 10 T × 10 m, and the PRD 2009 admission that a pure G-effect lab detector is useless. They agree with mainstream physics to within a factor of about 2. So does the classical algebra of the first-order cross term. What departs from mainstream:
- (a) the noise model that turns the cross term into a background-free photon count (refuted by JASON);
- (b) the fractal-membrane "magic" (JASON);
- (c) Baker's coherent N² summation plus a sub-diffraction beam (our inference);
- (d) Pais's S_max formula and its inputs (no physical basis; null Navy test).

### Cited Findings
- **Standard schematic conversion.** P_{h→γ} ∼ (BL)²/M_Pl² (Living Review eq. 128). The figure of merit is BL: about 10¹³ G km × (B/nG)(L/Gpc) for cosmological fields. — [arXiv:2501.11723](https://arxiv.org/abs/2501.11723)
- **Li et al. on the pure inverse-G effect.** Conversion 1.4×10⁻³³ and 7×10⁻¹⁷ photons/s for 10 T, 10 m, 0.01 m², h = 10⁻³⁰. The authors conclude it is "not available" for lab detection. This agrees with JASON (3-29)/(3-30). — [arXiv:0909.4118](https://arxiv.org/abs/0909.4118); [JASON](https://web.archive.org/web/2010/http://www.fas.org/irp/agency/dod/jason/gravwaves.pdf)
- **The factor ½ in the induced field.**
  - The 3DSR paper writes the induced field with the factor ½: Ẽ^(1) = (i/2)A B̂ k_g c z. — [arXiv:1712.00766](https://arxiv.org/abs/1712.00766)
  - PRD 2009 eq. (3) writes Ẽ^(1) ≈ A B̂ k_g c z without the ½. — [arXiv:0909.4118](https://arxiv.org/abs/0909.4118)
- **The interference-term algebra is common to both sides.** Li et al.'s eqs. (32)–(35) are identical to JASON (3-31)/(3-32); the disagreement is only over which photons count as noise. — [arXiv:0909.4118](https://arxiv.org/abs/0909.4118); [JASON](https://web.archive.org/web/2010/http://www.fas.org/irp/agency/dod/jason/gravwaves.pdf)
- **No experiment has demonstrated the scheme.** The Living Review says no experimental results suggest feasibility. — [arXiv:2501.11723](https://arxiv.org/abs/2501.11723)

### Inferences
- **Our first-order derivation.** For TT h_+ = h cos(kz−ωt) in transverse B0, the effective current is j^x = −B0∂_z h_+. The resonant solution of □A_x is A_x = −(B0h/2) z cos(kz−ωt), so |E| = ½ h B0 ω z.
  - This gives, in Gaussian units, P_conv = I_EM/I_GW = G B0² L²/c⁴.
  - In SI: P_conv = 4πG B² L²/(μ0 c⁴), which is (BL)²/(2M̄_Pl²) in Heaviside–Lorentz natural units with the reduced Planck mass M̄_Pl² = 1/(8πG). This is consistent with the Living Review's schematic form.
  - For 10 T and 10 m: 8.3×10⁻³⁴ (Li: 1.4×10⁻³³, about 1.7× higher). For JASON's case: 8.3×10⁻³⁷ (JASON: 10⁻³⁵, 4π higher).
  - Li's PRD 2009 field is 2× ours and their GW flux is 4× the one-polarization standard, so the conversion ratio partly cancels.
  - Treat the polarization convention (A⊕ = A⊗ = A "circular") as the source of the factor-2 ambiguities. Expose it as a parameter in code.
- **Standard GW energy flux.** F = (c³/16πG)⟨ḣ_+² + ḣ_×²⟩, which is c³ω²h²/(32πG) for one linear polarization of amplitude h **[mainstream]**. Li and Baker both use c³ω²A²/(8πG).
- **Suggested Rust calculator split** (each function takes `convention: Li | Standard`):
  1. `gb_peak_field(P, W0)` gives √(4Z0P/πW0²). Expect the Li value 1.26e3 to differ by about 10%.
  2. `gw_flux(h, f)` is the standard flux, with Li/Baker = 4× standard.
  3. `gertsenshtein_field(h, B, f, L)` is ½hBωL (standard) or hBωL (Li PRD 2009).
  4. `second_order_rate(...)` is ⟨S⟩Δs/ħω. Test: 7.0e-17 s⁻¹ (Li) for 10 T, 10 m, 0.01 m², h = 1e-30.
  5. `cross_term_rate(N0, N_gw)` is 2√(N0·N_gw). Test: N0 = 1.24e22 gives 3.54e2.
  6. `li_dt_min(N_other, N1)` is N_other/N1². Test: 3.0e3 s.
  7. `jason_t_min(N_gw)` is 1/N_gw. Test: 3.9e17 s.
  8. `jason_cavity(H, Q, ω, B, L)` covers eqs. 3-17 to 3-19.
  9. `quadrupole_pair_linear(r, F0, ω)` is (16G/15c⁵)r²F0²ω², with `baker_jerk(r, Δf, Δt)` = (32G/5c⁵)(2rΔf/Δt)².
  10. `baker_strain(F, ν)` is 1.28e-18√F/ν, against the standard 2× value.
  11. `diffraction_cap(λ, D)` checks Baker's 1.71/N.
  12. `pais_smax(fG, σ, v)` = fGσ²/ε0·v, with a breakdown check (σ/ε0 ≤ 1e9 V/m).
  13. `pais_hfgw(B, L, P)` as written, plus the SI-correct version.
  - Every "published" test vector is listed in Q1–Q4 above.

### Gaps
- No independent experimental HFGW limits in the GHz band were checked against the Li–Baker targets. The Living Review has these (for example magnetic-conversion and polarimetry limits), but extracting them was outside this pass.
- The factor-2 polarization conventions in the Li-group papers were not resolved against a single authoritative source. Raffelt & Stodolsky (1988) graviton–photon mixing is the standard reference, but it was not fetched this session, so the SI first-order formula above rests on the derivation shown here.
