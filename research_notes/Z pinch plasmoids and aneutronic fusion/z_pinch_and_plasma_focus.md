# Z-pinch, dense plasma focus, plasmoids, plasma cosmology and the pinch–gravity analogy: calculable formulas with test values

Conventions: SI units unless stated; T in eV where marked (1 eV = 1.602176634e-19 J); μ0 = 4π×10⁻⁷ H/m. Every item is labelled **[mainstream]**, **[contested]** or **[refuted]**. "computed" means I recomputed it with `z_pinch_checks.py` (in this folder, run `python3 -I z_pinch_checks.py`, stdlib only; all asserts pass). Published values are quoted as given in the source.

## 1. Z-pinch equilibrium, radiative collapse, m=0/m=1 instabilities, sheared-flow stabilisation (Zap Energy), Sandia Z and MagLIF

### Takeaway
The Bennett relation and the Pease–Braginskii current (≈1.4 MA) are textbook results that you can reproduce exactly from Spitzer resistivity and bremsstrahlung. Sheared-flow stabilisation is peer-reviewed. FuZE has gone from 1.25×10⁵ neutrons per pulse (2019) to more than 10¹⁰ (2025, company-reported), with total pressure around 1.6 GPa. MagLIF on Z set a new D–D record in 2026 of 2.33×10¹³ neutrons, about twice the 2020 record of 1.1×10¹³.

### Cited Findings

**Bennett relation [mainstream]**
- Formula: μ0 I₀²/(8π) = k T_B N_i (1+z), where N_i = 2π∫n_i r dr (ions per metre) and T_B is the density-weighted temperature. It is an integral equilibrium "valid for any profiles of p(r) and I(r)". The derivation comes from μ0 I₀² = 16π² ∫₀ᴿ r p dr. Practical form: T_B ≈ 3.12 I₀²/((1+z)N_i) [keV; MA, 10¹⁸ cm⁻¹] — [Klír 2005 thesis, CTU Prague, §2.1.1, eqs 2.12–2.16](https://arxiv.org/pdf/physics/0703207). The original paper is Bennett, Phys. Rev. 45, 890 (1934), cited there.
- Test value (computed): 1 MA with N = 10²⁰ m⁻¹ and z = 1 gives T_B = 1.560 keV, so the coefficient is 3.121 against the published 3.12. 1 MA with T_e + T_i = 2 keV needs N = 1.560×10²⁰ m⁻¹.

**Pease–Braginskii current and radiative collapse [mainstream]**
- Formula (Klír eq. 2.20): I_PB = (π/μ0)·√(48 lnΛ/(σ0 A))·(1+z)/z. It assumes an optically thin plasma, uniform T and current density, bremsstrahlung only, and a parabolic profile n_i = n_max(1 − r²/R²). The Spitzer conductivity is σ = σ0(kT_e)^{3/2}/(z lnΛ), and P_BS = A_z(kT_e)^{1/2}∫n_i² dA. The thesis gives "I_PB ≈ 1.4 MA" for deuterium with lnΛ ≈ 10. Above I_PB, radiation exceeds Joule heating and the pinch collapses. The thesis notes that I_PB can be "substantially changed by anomalous resistivity, leakage currents, axial magnetic field, … opacity …" — [Klír 2005 §2.1.2](https://arxiv.org/pdf/physics/0703207)
- Pease and Braginskii predicted the effect independently in 1957. Experiments have rarely seen a clean collapse — [Akel et al., Nukleonika 2016 (per search summary)](https://yadda.icm.edu.pl/baztech/element/bwmeta1.element.baztech-1e6d33d9-2af6-49d8-a970-efe0c3bce656/c/akel_Experiments_and_simulations.pdf). In fibre pinches, rapid expansion was seen even though the polyethylene I_PB (600 kA) was exceeded — [Klír 2005](https://arxiv.org/pdf/physics/0703207)
- Test value (computed), using NRL constants:
  - Inputs: η⊥ = 1.03×10⁻⁴ z lnΛ T_eV^{-3/2} Ω·m, which gives σ0 = 9.71×10³/e^{3/2}; bremsstrahlung 1.69×10⁻³⁸ z² n_e n_i T_eV^{1/2} W/m³.
  - z = 1, lnΛ = 10, σ⊥: I_PB = 1.37 MA.
  - Using σ∥ instead (σ∥ = 1.96 σ⊥): 0.98 MA.
  - lnΛ = 20: 1.94 MA.
  - I_PB scales as √lnΛ, and it depends on which conductivity you use. A Z-pinch has J ⊥ B, so σ⊥ is the consistent choice. Klír's σ_Spitzer ≈ 7.7×10⁻³ T_K^{3/2}/(z lnΛ) S/m is the ⊥ value.

**Sausage (m=0) and kink (m=1) instabilities [mainstream]**
- Perturbations go as ξ ∝ exp(i(mθ + kz) + Γt). For a Bennett equilibrium in ideal MHD, m=0 has growth Re(Γ) ~ v_A/R with v_A ≈ B_max/√(μ0ρ). The ideal growth rate grows without bound with k. Finite-Larmor-radius effects make the maximum finite — [Klír 2005 §2.1.3, eqs 2.21–2.23](https://arxiv.org/pdf/physics/0703207)
- Sharp-boundary (Kruskal–Schwarzschild 1954 geometry) dispersion relation, **derived by me, not taken from a fetched source**. It covers an incompressible column with a surface current, no internal B_z and a vacuum B_θ = B_a a/r: (Γa/v_A)² = [ka I_m′(ka)/I_m(ka)]·[1 + m²K_m(ka)/(ka K_m′(ka))], with v_A = B_a/√(μ0ρ). m=0 is always unstable. m=1 is unstable at all k. m ≥ 2 is stable for ka < ~m². Check it against Freidberg, *Ideal MHD*, before shipping.
  - Test values (computed), Γa/v_A at ka = 0.5, 1, 2, 5:
    - m=0: 0.348, 0.668, 1.181, 2.114
    - m=1: 0.481, 0.714, 1.075, 1.942
    - m=2: 0, 0, 0, 1.242
  - Closed-form check: m=0 at ka = 1 gives √(I₁(1)/I₀(1)) = 0.668.

**Sheared-flow stabilisation [mainstream for the measurements; scaling to gain is unproven]**
- Threshold: theory predicts that a sheared axial flow dV_z/dr > 0.1 k V_A stabilises the pure Z-pinch. This is Shumlak & Hartman, PRL 1995, cited in the 2001 paper. ZaP saw "a stable period which is over 700 times the expected instability growth time in a static Z-pinch". The measured shear was above threshold while the plasma was stable and near zero afterwards. The authors say causality cannot be proven — [Shumlak, Golingo, Nelson, Den Hartog, PRL 87 (2001) via UW-Madison copy](https://plasma2.physics.wisc.edu/pub_files/journal/Shumlak170.pdf). ZaP plasmas were 100 cm long with a 1 cm radius — [Shumlak papers summary](https://wx1.ans.org/pubs/journals/fst/a_13407). A 2005 poster reported global stability for about 2000 growth times — [APS DPP 2005](https://meetings-archive.aps.org/dpp/2005/bp1/20/)
- FuZE 2019 (Zhang et al., PRL 122, 135001):
  - Neutron emission lasted about 5 µs at about 200 kA, inside a 16 µs quiescent period.
  - Fill: 20% D₂ / 80% H₂.
  - Yield: (1.25 ± 0.45)×10⁵ neutrons per pulse, scaling as the square of the D fraction.
  - T: 1–2 keV; n ≈ 10¹⁷ cm⁻³; pinch radius 0.3 cm.
  - Source: [ADS](https://ui.adsabs.harvard.edu/abs/2019PhRvL.122m5001Z)
- FuZE 2020: up to 400 kA, T_i 1–2 keV, n 1–2×10²³ m⁻³, Y_n ~10⁷ per pulse, with yield strongly dependent on pinch current — [APS DPP 2020 abstract](https://meetings-archive.aps.org/dpp/2020/cp20/14/)
- Levitt et al., PRL 132, 155101 (2024): Thomson scattering measured T_e = 1–3 keV on axis, 20 cm downstream of the nose cone. T_e peaks at the same time as the neutron pulse, the pinch current and the inductive voltage — [Zap update](https://www.zapenergy.com/updates/2024/04/elevated-electron-temp); [phys.org](https://phys.org/news/2024-04-electron-temperatures-small-scale-stabilized.html)
- Ryan et al., arXiv:2408.05171 (2024), on FuZE: −25 kV charge, mean current 370 kA, 4×10⁷ D–D neutrons per discharge. The deuteron beam energy is below 7.4 ± 5.6 (stat) ± 3.7 (syst) keV, so most yield is isotropic and thermonuclear. It later appeared in Nuclear Fusion (DOI 10.1088/1741-4326/ada8bf, per Zap) — [arXiv](https://arxiv.org/abs/2408.05171); [Zap 2025 highlights](https://legacy.zapenergy.com/blog/2025-highlights)
- FuZE-Q 2024: single-shot neutron production above 5×10⁹, "triple last year's" — [Zap 2024 review (company)](https://legacy.zapenergy.com/blog/a-look-back-at-2024)
- 2025, company-reported and not yet peer-reviewed:
  - More than 11,700 plasmas across five devices; neutron yield > 10¹⁰; density > 7×10²⁴ m⁻³; pressure 1.6 GPa; shear velocity > 200 km/s.
  - FuZE-3 (three electrodes, commissioned 2025) passed 1 GPa.
  - Century ran 1,080 shots in 3 h in a DOE milestone, at 39 kW average input. Century uses H or He fuel, so it produces no neutrons.
  - Source: [Zap 2025 highlights](https://legacy.zapenergy.com/blog/2025-highlights)
- FuZE-3 (APS DPP, Nov 2025): electron pressure up to 830 MPa, about 1.6 GPa total estimated, T_e > 1 keV, n_e 3–5×10²⁴ m⁻³, lasting about 1 µs. No journal paper is cited yet — [ScienceDaily](https://www.sciencedaily.com/releases/2025/11/251120002836.htm)
- Earlier ZaP-HD: about 1 keV at 10²⁴ m⁻³, quiescent for 50 µs — [PPPL colloquium abstract](https://www.pppl.gov/events/2023/pppl-colloquium-sheared-flow-stabilized-z-pinch-approach-fusion-energy)

**Sandia Z and MagLIF [mainstream]**
- Z can drive 27 MA (Slutz & Vesey PRL 108, 025003, 2012) with a current pulse of about 100 ns — [Wikipedia MagLIF](https://en.wikipedia.org/wiki/Magnetized_liner_inertial_fusion)
- Gomez et al., PRL 113, 155003 (2014): 10 T, 2.5 kJ Z-Beamlet laser, 19 MA, up to 2×10¹² D–D neutrons — [Wikipedia MagLIF](https://en.wikipedia.org/wiki/Magnetized_liner_inertial_fusion). Physics World described about 10¹² neutrons and about 35 million degrees — [Physics World](https://physicsworld.com/a/new-approach-to-fusion-delivers-copious-neutrons/). The two sources disagree on the exact yield.
- Gomez et al., PRL 125, 155002 (2020): 1.1×10¹³ D–D neutrons (about 2 kJ DT-equivalent), burn-averaged T_i = 3.1 keV. The changes were field 10.4 → 15.9 T, preheat 0.46 → 1.2 kJ and current 16 → 20 MA — [Wikipedia MagLIF](https://en.wikipedia.org/wiki/Magnetized_liner_inertial_fusion); [phys.org](https://phys.org/news/2020-10-neutron-sandia-labs-machine-fusion.html); [IAEA overview](https://conferences.iaea.org/event/214/contributions/17017/contribution.pdf)
- **2026 record:** Weis et al., Phys. Plasmas 33, 062711 (published 18 June 2026):
  - D–D yields: 2.11×10¹³ and 2.33×10¹³.
  - T_i: 3.0 and 3.3 keV; stagnation pressures 1.6 and 1.3 Gbar; Lawson parameter χ = 0.2 and 0.1.
  - Drive: AR9 Be liner with a dielectric coating, about 2.2 kJ preheat, 18.5 MA, 15 T.
  - Simulations suggest an ice layer could raise yield up to 2.5×.
  - Source: [AIP](https://pubs.aip.org/aip/pop/article/33/6/062711/3395915/Improving-the-stability-and-performance-of-MagLIF)
- A projected 30–50 kJ DT yield on Z is a projection, not a result — [newswise/Sandia](https://www.newswise.com/articles/record-neutron-numbers-at-sandia-labs-z-machine-fusion-experiments)

### Inferences
- FuZE 2019 Bennett consistency (computed): I = 200 kA, n0 = 10²³ m⁻³ and a = 3 mm give N = πa²n0 = 2.83×10¹⁸ m⁻¹ for a Bennett profile, so T_e = T_i = 2.21 keV. The published value is 1–2 keV, so the data are close to Bennett balance. Edge field B = 13.3 T, v_A ≈ 920 km/s and a/v_A = 3.3 ns. The 16 µs quiescent period is therefore about 4900 static growth times. That is the same order as the ZaP figures of more than 700 and about 2000.
- Shear threshold at ka = 1 for FuZE 2019 (computed): dV_z/dr > 3.1×10⁷ s⁻¹, about 92 km/s across the radius. Zap's 2025 "> 200 km/s" shear is above this scale.
- FuZE-3 pressure check (computed): n_e = 5.18×10²⁴ m⁻³ at 1 keV gives exactly 830 MPa. The "1.6 GPa total" assumes T_i ≈ T_e.
- MagLIF DT-equivalent (computed): 2 kJ / 17.59 MeV = 7.10×10¹⁴ DT reactions, about 65 per D–D neutron at about 3 keV. Scaling the 2026 2.33×10¹³ D–D the same way gives about 4.2 kJ DT-equivalent. This is my rough number, not a published one, and it ignores the T_i difference.
- Neutron yield rising steeply with pinch current, with 1.25×10⁵ at 200 kA, about 10⁷ at 400 kA and 4×10⁷ at 370 kA, is consistent with a high power of I. The FuZE data do not pin the exponent.

### Gaps
- I did not open Haines 2011 (PPCF 53, 093001) or Ryutov, Derzon & Matzen 2000 (RMP), so the Haines prefactor for I_PB is not checked. I used Klír's derivation instead.
- I did not find a published FuZE-Q peak current, bank energy or pinch radius. FuZE-3's current is also unknown. Zap's 2025 numbers are not peer-reviewed yet.
- Goyon et al. (Phys. Plasmas 2024; T_e up to 2.25 keV, n up to 4.9×10¹⁷ cm⁻³) appeared only in a search snippet, and I did not open the primary.
- I did not open or verify the Kadomtsev m=0 criterion (−d ln p/d ln r < 4γ/(2+γβ)), so it is not included as a sourced formula.
- I found no published Z stored energy in what I fetched. The Shumlak–Hartman 1995 PRL full text was not retrieved.

## 2. Dense plasma focus: Lee model scaling and saturation, LPPFusion / Focus Fusion plasmoids and beams, p-¹¹B status, Bostick 1956

### Takeaway
The Lee model gives reproducible scaling: Y_n ≈ 1.8×10¹⁰ I_peak^3.8 and 3.2×10¹¹ I_pinch^4.5 (MA). Yield saturates because the run-down "dynamic resistance" (about 7 mΩ) caps I_peak near V₀/DR₀. LPPFusion's FF-1/FF-2B results are peer-reviewed: confined ion energies above 200 keV, 2.5×10¹¹ neutrons and Z_eff ≈ 1.004. Their plasmoid model and the QMF effect remain [contested]. As of 2026 there is still no measured p-¹¹B fusion yield.

### Cited Findings

**Lee model and saturation [mainstream]**
- Dynamic resistance: tube inductance is L = (μ/2π) ln(c) z with c = b/a, so dL/dt = 2×10⁻⁷ (ln c) dz/dt. With dz/dt ≈ 10 cm/µs (about the same for all D₂ foci) and c = 2, dL/dt = 1.4×10⁻² H/s. The power into a changing L has a resistive part ½(dL/dt)I², which gives DR₀ = ½ dL/dt ≈ 7 mΩ. Only about 0.7 of the current drives the sheath. Lee models this with f_m and f_c, and including f_c lowers DR₀ by about 40% — [Lee, APL 95, 151503 (2009)](https://www.plasmafocus.net/IPFS/2009/09PP8APL.pdf)
- Saturation: as C₀ → ∞, Z₀ = √(L₀/C₀) → 0 and I_peak → V₀/DR₀. Y_n ~ E₀² at low energy falls to about E₀^0.8 towards 25 MJ, "deteriorating already at several hundred kilojoules". Relations: Y_n ~ I_pinch^4.5 and Y_n ~ I_peak^3.8 — [Lee 2009](https://www.plasmafocus.net/IPFS/2009/09PP8APL.pdf)
- Table I test vectors, with V₀ = 30 kV, L₀ = 30 nH, DR₀ = 7 mΩ, stray resistance 0.1Z₀ and I_peak = V₀/Z_total. Published (E₀ kJ, C₀ µF, I_peak kA): (0.45, 1, 152), (4.5, 10, 447), (45, 100, 1156), (135, 300, 1676), (450, 1000, 2321), (1080, 2400, 2781), (4500, 10⁴, 3407), (45000, 10⁵, 4209) — [Lee 2009](https://www.plasmafocus.net/IPFS/2009/09PP8APL.pdf)
- Fits with prefactors: Y_n = 3.2×10¹¹ I_pinch^4.5 for I_pinch 0.2–2.4 MA, and Y_n = 1.8×10¹⁰ I_peak^3.8 for I_peak 0.3–5.7 MA. Y_n ~ E₀^2.0 at tens of kJ and E₀^0.84 at MJ (up to 25 MJ). The method fits the computed current waveform to the measured one — [Saw & Lee summary (KUSET)](https://journals.ku.edu.np/kuset/article/view/413)
- Below 500 J (six machines): Y_n ~ E₀^2.1, I_peak^3.9, I_pinch^5.0 — [Singh, Lee, Saw (KUSET 2014)](https://nepjol.info/index.php/KUSET/article/view/63652)
- Experimental data fit Y_n ∝ I_p^3.5 against a simple-scaling I^4. Exponents above 4 are possible if T rises with I — [Hassan et al., ICFPPT 2011 (IAEA TECDOC)](https://www-pub.iaea.org/MTCD/Publications/PDF/TE-1713-CD/html/s10_4.htm)
- A 2026 arXiv preprint argues that suppressing filamentation could beat saturation. It is unreviewed, so treat it as [contested] — [arXiv:2606.27108](https://arxiv.org/pdf/2606.27108)

**LPPFusion / Focus Fusion (Lerner et al.)**
- Proposed mechanism [contested]: the current sheath breaks into filaments, which merge and kink into a plasmoid "only tens of microns across". As the magnetic field decays it induces an electric field that "generates a beam of electrons in one direction and a beam of ions in the other". Plasmoid electron densities of 10¹⁹–10²² cm⁻³ are cited from the literature — [Lerner et al., J. Fusion Energy 42:7 (2023)](https://www.lppfusion.com/storage/s10894-023-00345-z.pdf)
- Lerner model equations [contested], as printed. The PDF extraction lost the exponent signs; I fixed them using the paper's own worked examples.
  - Eq 1, core radius: r_c ≈ 1.32×10⁻³ (μz)^{2/3} r (cm).
  - Eq 2, core field: B_c = 4zμ(M/m)B (G).
  - Eq 3, core density: n_c ≈ 3.7×10¹⁰ μ²zI²/r² (I in A, r in cm).
  - Eq 4, yield: Y ∝ μ^2.75 I⁴ f(T_i).
  - Symbols: μ is the mean ion mass number, z the ion charge, M/m the proton/electron mass ratio, B the peak field at the cathode, and r the cathode radius.
  - The paper's example: B = 100 kG gives B_c = 10 GG for p-¹¹B — [Lerner 2023](https://www.lppfusion.com/storage/s10894-023-00345-z.pdf)
- Quantum magnetic field (QMF) effect [contested]:
  - Landau energy: E_b = (n + ½)ħω_ce ≈ 11.6 eV × B(GG).
  - Ion→electron energy transfer is suppressed when E_i < (M/m)E_b. For E_i = 300 keV this needs B > 14 GG (p), 3.5 GG (α) or 1.3 GG (¹¹B).
  - The transfer rate can drop "by a factor as large as 25", with x-ray emission reduced up to 5×.
  - 0-D simulations claim ignition with up to 80% burn-up, for r_min 14 µm and n_i up to 4×10²³ cm⁻³.
  - Source: [Lerner 2023](https://www.lppfusion.com/storage/s10894-023-00345-z.pdf)
- Viscous ion heating [contested]: T_i = 6.2×10⁴ z_eff^1.6 n_i^0.2 (lnΛ L_max B)^0.4 (Eq 7). An electron-beam ion-heating condition is I > 5.0/(μ^{1/6} z^{5/6}) MA (Eq 9). The impurity criterion requires f z_i² < 1.4 (z_eff < 1.2) at the end of run-down — [Lerner 2023](https://www.lppfusion.com/storage/s10894-023-00345-z.pdf)
- FF-1 hardware [mainstream]: 12 capacitors, 113 µF, 45 kV max, 115 kJ max (usually run with 8). Charge 24–40 kV, peak current "slightly over 1 MA", rise time 1.8 µs. Anode radius 2.8 cm, cathode radius 5 cm. Renamed FF-2B in 2019 when the electrodes were changed to Be — [Lerner 2023](https://www.lppfusion.com/storage/s10894-023-00345-z.pdf)
- 2012, Phys. Plasmas 19, 032704 [mainstream as a measurement; its interpretation is contested]: ions above 150 keV confined for about 7–30 ns, and up to 1.5×10¹¹ neutrons per pulse, attributed to confined ions rather than a beam — [NextBigFuture abstract reproduction](https://www.nextbigfuture.com/2012/03/new-cathode-design-for-dense-plasma.html); [Lerner 2023](https://www.lppfusion.com/storage/s10894-023-00345-z.pdf)
- 2017, Phys. Plasmas 24, 102708 [mainstream measurement; interpretation contested]:
  - Ion energy: best shot mean 240 ± 20 keV; ten-shot mean 124 keV.
  - Isotropic neutrons, with no cold dense component detected. X-ray-derived total density was 0.8–2.5× the hot-ion density.
  - Yield: record 2.5 ± 0.25×10¹¹ neutrons, 0.2 ± 0.02 J, on 23 May 2016, from 60 kJ input, so η = 3.3×10⁻⁶.
  - Plasmoid volume 1.2×10⁻³ cm³. n ≥ 3×10¹⁹ cm⁻³ (upper limit 4×10¹⁹), which is 20× the 1.6×10¹⁸ cm⁻³ fill. Confinement about 40 ns.
  - nτT = 3.4 ± 0.8×10²⁰ keV·s/m³.
  - Source: [Lerner 2023](https://www.lppfusion.com/storage/s10894-023-00345-z.pdf)
- Beams [mainstream measurement]: the 28 Feb 2013 ion beam had 120 kA peak, 3 MeV mean ion energy (time of flight), 5 ns and 2 kJ, with "4 kJ" including the opposite electron beam. Inductance rose about 10 nH in the best pinches, about 4–5 kJ into the pinch — [Lerner 2023](https://www.lppfusion.com/storage/s10894-023-00345-z.pdf)
- Plateau [mainstream]: DPF yield plateaus above about 1 MA. The best wall-plug efficiency is still Speed-2 (2001), 1 J from 80 kJ (1.25×10⁻⁵). LPP has not beaten its own 2016 yield — [Lerner 2023](https://www.lppfusion.com/storage/s10894-023-00345-z.pdf)
- Purity [company claim, published in JFE]:
  - Be deposition is 6 ± 1 µg per shot into the pinch, with 3.3 mg of gas in the sheath. That is ≤ 0.2% by mass, giving z_eff = 1.004 (1.012 worst case).
  - JET's best was 1.2.
  - LPP projects 25 J if the plasmoid shrinks from 250 µm to 50 µm.
  - Source: [Lerner 2023](https://www.lppfusion.com/storage/s10894-023-00345-z.pdf)
- p-¹¹B status [mainstream fact; the projections are contested]:
  - 2025: first "tiny pinch" with almost pure decaborane, showing a plasmoid "but not yet measurable fusion reactions".
  - Boron separated from hydrogen and deposited Be–B on the anode, which cracked it in August 2025.
  - The proposed fix is a slower run-down — [LPPFusion year-end summary (1 Jan 2026)](https://www.lppfusion.com/year-end-lppfusion-summary/); [LPP posts](https://www.lppfusion.com/results-confirm-theory-show-way-back-to-boron/)
- Other p-¹¹B: HB11/LFEX produced 1.4×10¹¹ alphas (66 mJ) from 6.5 MJ input, so η = 1.0×10⁻⁸ (as cited by Lerner) — [Lerner 2023](https://www.lppfusion.com/storage/s10894-023-00345-z.pdf)

**Bostick 1956 plasmoids [mainstream, historical]**
- Bostick, "Experimental Study of Ionized Matter Projected across a Magnetic Field", Phys. Rev. 104, 292–299 (1956). A plasma gun projected plasmoids across a field at up to 10⁷ cm/s, and up to 2×10⁷ cm/s in an Armour Research Foundation report. A perpendicular field halves the speed and stretches the plasmoid along B — [per search summary of Bostick reports](https://digital.library.unt.edu/ark:/67531/metadc714766/); [DTIC](https://apps.dtic.mil/sti/pdfs/AD1038383.pdf). Lerner credits the Bostick–Nardi model as the origin of the plasmoid picture — [Lerner 2023](https://www.lppfusion.com/storage/s10894-023-00345-z.pdf)

### Inferences
- Lee Table I (computed): all rows up to 4.5 MJ reproduce to within 1.2%. The paper rounds Z_total, so allow about 1.5%. **The 45 MJ row is internally inconsistent.** Z₀ = 0.55 mΩ gives Z_total = 7.6 mΩ, so I_peak = 30 kV/7.6 mΩ = 3946 kA, not 4209 kA. Don't use that row as a test vector. Asymptote: V₀/DR₀ = 4.29 MA.
- Snowplow run-down speed (derived/computed): from ρ₀f_mπ(b²−a²)v² = μ0(f_cI)² ln(b/a)/(4π), v = (f_c/√f_m)·√(μ0 ln c/(4π²(c²−1)ρ₀))·I/a. For FF-1 (a = 2.8 cm, b = 5 cm, 1 MA, D fill 1.6×10¹⁸ cm⁻³) this gives 4.5 cm/µs with f_m = f_c = 1, and 9.9 cm/µs with f_m = 0.1, f_c = 0.7. That matches Lee's "about 10 cm/µs".
- FF-1 against the Lee fit (computed): Y(I_peak = 1.1 MA) = 2.6×10¹⁰, while FF-1 measured 2.5×10¹¹, about 10× above the fit. This fits Lerner's claim that FF-1 beat other DPFs at its current.
- LPP number checks (computed):
  - QMF thresholds: 11.58 eV/GG; 14.11 GG (p), 3.55 GG (α), 1.29 GG (¹¹B).
  - B_c for p-¹¹B (z = 3, μ = 6) at 100 kG: 13.2 GG, against "10 GG".
  - Eq 9: 4.45 MA for D and 1.48 MA for p-¹¹B.
  - nτT: 3.5×10¹⁹ cm⁻³ × 40 ns × 240 keV = 3.36×10¹⁴ keV·s/cm³ ✓.
  - Ion beam: 120 kA × 3 MV × 5 ns = 1.8 kJ, against 2 kJ.
  - **0.2 J from 2.5×10¹¹ neutrons:** this is 0.13 J counting the n branch only, or 0.29 J counting both D–D branches. The published value falls between the two, so the convention is unclear.
- Two things are separate here: the I⁴ scaling and plateau are mainstream, but the claim that impurities and filament disruption cause the plateau, and that it can therefore be beaten, is contested. Lee attributes the plateau to circuit physics instead.

### Gaps
- I did not open the 2012 and 2017 Phys. Plasmas papers directly. The 2012 ion energy distribution and B estimates come only through Lerner 2023 and a blog abstract.
- I found no independent replication of LPP's ion-energy or QMF claims, and no peer-reviewed critique.
- The 2016 vs 2017 date for the 200 keV record differs between LPP pages.
- A March/May 2026 FF-2B reassembly is reported in LPP news posts (search summary only). I found no 2026 p-¹¹B yield.
- The "Lee & Saw J. Fusion Energy 27, 292 (2008)" source of the I^4.5 fit is cited in Lee 2009, but I did not open it.

## 3. Peratt's Z-pinch plasma cosmology (Birkeland currents, galaxy simulations, "Peratt instability", petroglyphs) and mainstream assessment

### Takeaway
Peratt ran 3-D electromagnetic particle-in-cell (PIC) simulations of two or more galaxy-scale current filaments (about 10¹⁸ A, about 300,000 ly apart), published in IEEE TPS in 1986. He claimed they reproduce double radio galaxies, spiral morphology and flat rotation curves without dark matter. Mainstream astrophysics rejects plasma cosmology as a cosmology (CMB and light-element abundances), and has not observed the required galactic Birkeland currents. The petroglyph Z-pinch aurora work is peer-reviewed in IEEE TPS but fringe.

### Cited Findings
- [contested] Peratt (LANL) and Alfvén outlined the "plasma universe". Peratt simulated colliding plasma clouds that he said mimic galaxy shapes (Peratt, Green & Nielson, PRL 1980). He proposed galaxies form when filaments "300,000 light years apart" carrying Birkeland currents of 10¹⁸ A join in a Z-pinch. He claimed jets resembling quasars without supermassive black holes, and flat rotation curves without dark matter — [Wikipedia: Plasma cosmology](https://en.wikipedia.org/wiki/Plasma_cosmology)
- [contested] Peratt, IEEE TPS 14, 639 (1986), Part I: the model is "the interaction between field-aligned current-conducting, galactic-dimensioned plasma sheets or filaments (Birkeland currents)". It uses "fully three-dimensional electromagnetic particle-in-cell simulations". The synchrotron emission was said to match double radio galaxies and quasars in power, morphology, spectra, polarisation and jets. Part II (IEEE TPS 14, 763) ran the interactions over 10⁹–5×10⁹ yr and produced elliptical, peculiar, barred and spiral forms, compared with H I distributions and rotation curves — [OpenAlex record / abstract](https://openalex.org/W2000605847); [AIP-hosted letter](https://aip.brightspotcdn.com/PTO.v38.i9.13_1.online.pdf)
- [refuted, as a cosmology] Mainstream view:
  - Large-scale Birkeland currents have not been observed.
  - The charge-neutrality length is far smaller than cosmological scales (Siegel & Fry 2006).
  - No plasma-cosmology model reproduces the CMB anisotropy spectrum (WMAP and Planck fit the Big Bang).
  - Alfvén–Klein nucleosynthesis would overproduce X- and γ-rays.
  - Peebles (1993): "there is no way that the results can be consistent with the isotropy of the cosmic microwave background radiation."
  - Source: [Wikipedia: Plasma cosmology](https://en.wikipedia.org/wiki/Plasma_cosmology). Peebles' criticism targets Alfvén–Klein cosmology, not the galaxy PIC runs specifically.
- [contested] Petroglyphs: Peratt, "Characteristics for the Occurrence of a High-Current Z-Pinch Aurora as Recorded in Antiquity", IEEE TPS 31(6), Dec 2003 (doi 10.1109/TPS.2003.820956). It classifies petroglyphs against laboratory instability morphologies and proposes an intense prehistoric aurora. Part II (IEEE TPS 35(4), 2007) compares 84 instability configurations, and a survey of "several million" field-of-view sites claims an axial-south alignment. There is no independent replication — [APS 2005 abstract](https://meetings-archive.aps.org/apr/2005/l1/30); [APS DPP 2006 abstract](https://archive.aps.org/dpp/2006/gp1/132)
- Peratt cited Bostick's 1950s laboratory plasmoid discharges as galaxy-like precursors — [Wikipedia: Plasma cosmology](https://en.wikipedia.org/wiki/Plasma_cosmology)

### Inferences
- Filament numbers (computed): two parallel 10¹⁸ A filaments at d = 300,000 ly = 2.84×10²¹ m give B = μ0I/(2πd) = 7.0×10⁻¹¹ T = 0.7 µG at the partner filament, and F/L = μ0I²/(2πd) = 7.05×10⁷ N/m. Those are µG-scale fields, the same order as observed galactic fields, so the field strength alone does not rule the model out. The objections are the missing current-closure observations and cosmological data.
- The simulations are real numerical work published in a peer-reviewed engineering journal. The mainstream response has mostly been neglect rather than formal rebuttal. "Refuted" applies to plasma cosmology as an alternative to ΛCDM, not to each PIC result.

### Gaps
- I found no peer-reviewed astrophysics rebuttal that addresses Peratt's 1986 galaxy PIC runs specifically. Critiques found were either general (cosmology) or informal forum posts, and I did not cite the forum posts.
- I found no simulation particle count, grid or scaling factors (the PIC-to-galaxy mapping), because I did not open the IEEE papers.
- I found no source using the term "Peratt instability". It seems to be used informally in Electric Universe circles for the diocotron/Z-pinch column sequences in the petroglyph papers. That is unverified.

## 4. The pinch–gravity analogy: Bennett pinch ≡ self-gravitating isothermal cylinder (Ostriker 1964)

### Takeaway
Under a uniform drift u the pinch and the gravitating cylinder are the same equilibrium: the Bennett pinch is Ostriker's isothermal cylinder with G m² ↔ μ0e²u²/(4π) and c_s² ↔ k(T_e+T_i)/m_i. Both have the profile n ∝ (1 + r²/r₀²)⁻² and a fixed line mass or line density set by temperature alone: 2c_s²/G on the gravity side, Bennett's N on the pinch side. Magnetic attraction in a 1 MA pinch is about 5×10²⁷ times stronger than the plasma's own gravity.

### Cited Findings
- [mainstream] Bennett: μ0I²/(8π) = N k(T_e+T_i), valid for any profile — [Klír 2005](https://arxiv.org/pdf/physics/0703207). The Bennett-profile vector potential is A_z = −A₀ ln(1 + (r/r_p)²), which gives n ∝ (1 + r²/r_p²)⁻² — [arXiv:2411.06674 (search summary)](https://arxiv.org/pdf/2411.06674)
- [mainstream] Stodółkiewicz (1963) and Ostriker (1964) found the isothermal filament independently, with line mass M_line = 2kT/(Gμm_H) = 2c_s²/G. Above this, cylinders collapse radially and fragment. With turbulence, M_line,vir = 2σ_tot²/G. For μ = 2.3 and T = 10 K, M_line = 16.6 M_⊙/pc, with 16.3 (T/10 K) quoted elsewhere — [arXiv:2307.11162](https://arxiv.org/pdf/2307.11162); [arXiv:2604.04501](https://arxiv.org/pdf/2604.04501); ≈16 M_⊙/pc × (T/10 K) — [André et al., arXiv:1309.7762](https://arxiv.org/pdf/1309.7762). In Aquila, filaments with λ ≳ 5×(2c_s²/G) host protoclusters — [arXiv:1309.7762](https://arxiv.org/pdf/1309.7762)
- Ostriker profile [mainstream; standard form]: ρ(r) = ρ_c/(1 + r²/(8H²))², with H² = c_s²/(4πGρ_c). The integrated line mass is 2c_s²/G, independent of ρ_c. I checked this numerically; the form comes from Ostriker 1964, ApJ 140, 1056, which I did not fetch.

### Inferences
- **Mapping (derived here, mainstream physics):**
  - Gravity: g(r) = 2Gμ(r)/r, so the force density is ρg = 2Gρμ(r)/r.
  - Pinch: B = μ0I(r)/(2πr) with J = enu, so JB = n(μ0e²u²/2π)N(r)/r.
  - Matching the two gives G m_i² ↔ μ0e²u²/(4π) and c_s² ↔ k(T_e+T_i)/m_i.
  - Result: line mass 2c_s²/G ↔ N = 8πk(T_e+T_i)/(μ0e²u²). With I = euN this becomes μ0I²/(8π) = Nk(T_e+T_i), which is Bennett exactly. Computed check: u = 300 km/s and T_e+T_i = 2 keV give N = 2.774×10¹⁸ m⁻¹ and I = 133 kA, which satisfies Bennett to 10⁻¹².
- **Virial forms:**
  - Gravity: 2K + W = 0 with K per length = ∫p dA = μc_s² and W = −Gμ² for any cylindrical distribution, giving μ = 2c_s²/G.
  - Magnetic: 2∫p dA = μ0I²/(4π), the Bennett relation, since 2Nk(T_e+T_i) = μ0I²/(4π).
  - The correspondence is Gμ² ↔ μ0I²/(4π). The full mapping needs uniform u. The integral relations hold for any profile.
- **Worked gravity example (computed):** T = 10 K, μ = 2.3 gives c_s = 0.189 km/s and M_line,crit = 1.075×10¹⁵ kg/m = 16.7 M_⊙/pc, against the published 16.6.
- **Worked pinch example (computed):** a D plasma at 1 MA with T_e+T_i = 2 keV has N = 1.56×10²⁰ m⁻¹, a line mass of 5.22×10⁻⁷ kg/m.
  - Magnetic/gravitational binding ratio: μ0I²/(4π) ÷ G(Nm_D)² = 5.5×10²⁷.
  - To confine the same plasma by its own gravity would take 2c_s²/G = 2.87×10²¹ kg/m.
- Equivalent drift (computed): two protons attract magnetically as strongly as gravitationally when their relative drift is u = m_p√(4πG/μ0)/e = 2.70×10⁻¹⁰ m/s. This is a neat Rust test constant.
- The analogy is formal. Both are cylinders of an isothermal gas held by a 1/r attractive force, and the sign works because parallel currents attract. It does not imply any coupling between pinches and gravity [refuted if claimed as a physical gravity effect]. I found no claim of that kind in the sources I read.

### Gaps
- I found no published paper that states the Bennett ↔ Ostriker correspondence explicitly. A search for it found nothing, so the mapping here is my own derivation. Ostriker 1964 itself was not fetched.

## 5. DoD and national-lab pinch and plasmoid programmes: AFRL Shiva Star and MARAUDER, Navy and DTRA, Lockheed CFR, directed-energy plasmoids

### Takeaway
On the public record:
- MARAUDER: Phillips Lab/AFRL formed and accelerated 0.1–2 mg compact toroids on the 9.4 MJ Shiva Star bank from late 1990. The goal was above 1000 km/s. 3,000 and 10,000 km/s were projections, and no measured final speed is published after 1995.
- Lockheed CFR: its public claims (100 MW in 2.1×3.0 m) were never backed by fusion data, and the effort was halted before 2021.
- Navy: the Pais "plasma compression fusion device" is a patent with no experimental data.

### Cited Findings
- [mainstream] Shiva Star/MARAUDER: CT experiments ran "on the 9.4 MJ Shiva Star fast capacitor bank at Phillips Laboratory (Kirtland AFB, New Mexico) since late 1990" with 0.1–1.0 mg rings — [Peterkin, Degnan et al. (DTIC ADA275568)](https://apps.dtic.mil/sti/pdfs/ADA275568.pdf). Degnan et al., Phys. Fluids B 5(8), 2938 (1993), DOI 10.1063/1.860681, used toroids of 0.5–2 mg. CTs form in 7–15 µs and "can be accelerated intact". The goal was to compress toroids to high density and field and accelerate them to high speed — [search summary of Degnan 1993](https://ilse.kobv.de/KobvIndexRecord/nat_lic_papers_NLZ219268614). An AF study's aim was "velocities in excess of 1,000 km/s" for a nonnuclear test system — [DTIC ADA639263](https://apps.dtic.mil/sti/tr/pdf/ADA639263.pdf)
- [mainstream for the hardware; the projections were not demonstrated] MARAUDER summary:
  - 1–2 mg in a 1 m diameter toroid, up to 10¹¹ g acceleration (1993 gun).
  - Projected 3,000 km/s (1995) and 10,000 km/s (2000).
  - Megajoule-class kinetic energy from a 5–10 MJ coaxial gun; one shot was compared to about 5 lb TNT (9.5 MJ).
  - No published fate after 1995. The related LLNL RACE (Ring Accelerator Experiment) dates to 1991.
  - Source: [Wikipedia: MARAUDER](https://en.wikipedia.org/wiki/MARAUDER)
- [contested; never demonstrated] Lockheed CFR:
  - Began 2010 at Skunk Works. In Oct 2014 it promised a test in under a year and a prototype within 5 years.
  - Claimed a 100 MW D–T unit of 2.1×3.0 m.
  - McGuire's 2015 concepts: 200 t (15 T, 1 m shielding) to 2,000 t (5 T, 2 m shielding). A 200 MW_th, 18×7 m, about 2,000 t concept was also described.
  - T4 (2015): T_e 20 eV, n_e 10¹⁶ m⁻³, under 1% ionisation, 3 kW.
  - T4B (design): 1 m × 2 m, 1 MW/25 keV neutral beams for 3 ms, n = 5×10¹⁹ m⁻³, β = 1, 0.1 T, 0.2 m³, 1,170 J, T_i 75 eV, T_e 250 eV.
  - TX (design): 7 m × 18 m, 320 MW gross, 5×10²⁰ m⁻³, 2.3 T, 16.3 m³, 51 MJ, T_i 9.6 keV, T_e 12.6 keV.
  - T5 was due online at the end of 2019. Aviation Week (2023) reported "Skunk Works Halted Nuclear Fusion Effort Before 2021".
  - Steven Cowley questioned the small size, noting that doubling size gives about 8× better confinement.
  - Source: [Wikipedia: Lockheed Martin CFR](https://en.wikipedia.org/wiki/Lockheed_Martin_Compact_Fusion_Reactor); [The War Zone 2019](https://www.twz.com/29074/skunk-works-exotic-fusion-reactor-program-moves-forward-with-larger-more-powerful-design)
- [contested; no data] Navy (NAWCAD) Pais "plasma compression fusion device": application filed 22 March 2018, published 26 September 2019. It claims a 0.3–2 m device producing gigawatt to terawatt output "with input power only". Experts say the claims need "considerable proof" — [TechXplore](https://techxplore.com/news/2019-10-patent-plasma-compression-fusion-device.html); [The War Zone](https://www.twz.com/30256/scientist-behind-the-navys-ufo-patents-has-now-filed-one-for-a-compact-fusion-reactor.md). Reports of an issued patent (US10620121B2) rest on an unvetted secondary source — [Notion analysis](https://supt.notion.site/Salvatore-Pais-Patent-Analysis-Advanced-Propulsion-Energy-Systems-22f8f11f547a80e19e6ddeec9abdd9da)

### Inferences
- MARAUDER kinematics (computed): KE = ½mv² gives 0.5 MJ for 1 mg at 1,000 km/s, 9.0 MJ for 2 mg at 3,000 km/s (which matches the "5 lb TNT ≈ 9.5 MJ" figure) and 50 MJ for 1 mg at 10,000 km/s. At 10¹¹ g, reaching 1,000 km/s takes 0.51 m and 1.02 µs, which is consistent with a metre-scale coaxial gun and a µs-scale bank. The 9 MJ case would need almost all of Shiva Star's 9.4 MJ, so the 3,000 km/s and 2 mg projection implies near-100% coupling. That is implausible, and it fits the absence of any published demonstration.
- Lockheed checks (computed): T4B gives p = 2.60 kPa against B²/2μ0 = 3.98 kPa, so β = 0.65, not 1. Its thermal energy is 1.5pV = 781 J against the 1,170 J quoted; the difference may be fast-ion energy. TX gives β = 0.84 and a thermal energy of 43.5 MJ against 51 MJ quoted. The design tables are internally approximate.
- Bostick-type plasmoids (100–200 km/s) and MARAUDER CTs (aimed at ≥ 1,000 km/s) are the same family of objects: magnetised plasma rings fired from coaxial guns. Only the CT programme had a declared directed-energy or weapons-effects purpose ("nonnuclear test system").

### Gaps
- I found nothing public on DTRA pinch or plasmoid work beyond the general knowledge that Z-pinch x-ray simulators are used for radiation-effects testing. No numbers or sources were retrieved, so nothing is claimed here.
- I found no Navy plasmoid weapons data. I found no declassified directed-energy plasmoid results with measured speeds or energies beyond the MARAUDER figures above.
- I found no MARAUDER measured final CT velocity. The Degnan 1993 full text was not opened.
- I found no RACE (LLNL) numbers.
- No Lockheed T5 results were ever published, and no Lockheed statement was found confirming the shutdown. That rests on the Aviation Week 2023 report cited by Wikipedia.
