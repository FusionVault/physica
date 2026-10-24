# Tunnelling crosses the barrier; the vacuum holds firm

Nothing in physics breaks the Coulomb barrier. Every fusion system that works, from the Sun's core to the National Ignition Facility (NIF), tunnels through it at the Gamow peak, 30–45 times below the barrier top, and the mainstream formulas for that tunnelling (Gamow factors, Bosch–Hale reactivities and the Wurzel–Hsu Lawson criterion) reproduce their published tables to within 0.3% *(check)*. NIF's record shot turned **2.08 MJ of laser light into 8.6 MJ, a target gain of 4.13**, and stars fuse only because gravity holds a 1.35 keV plasma together for billions of years. The levers that act on the barrier itself are small: electron screening multiplies rates by 1.2–3 at the 2.5–10 keV energies where it is measured, a muon catalyses **150 fusions where 270–340 would repay its cost**, and pyroelectric crystals, fusors and NASA's "lattice confinement" return 10⁻⁶ to 10⁻¹⁴ of their input as fusion energy. The Schwinger field, **1.32×10¹⁸ V/m**, has never been reached: the brightest laser falls **1,454× short in field and 2.1×10⁶ in intensity**, the strong-field regime has been entered only in boosted frames (χ ≈ 0.3 at SLAC E144, up to 2.4 in CERN NA63's crystals), and laboratory vacuum birefringence remains 7× below detection. Every fringe claim checked here fails by a margin that can be computed: Fleischmann and Pons saw **3.7×10⁷ times fewer neutrons** than D–D fusion of their heat requires, Pais's fusion device yields about 200 W on its own parameters rather than gigawatts, Widom–Larsen overstate their own rate 300×, and Rossi's Lugano COP falls from 3.2–3.6 to 1.13–1.32 once the emissivity is corrected. The one live contested frontier is narrow: the 300–700 eV screening measured in deuterated metals, which now underpins DARPA's 2026 MARRS programme even though its phase-2 goal of 10⁶ reactions s⁻¹ g⁻¹ is only 6×10⁻⁷ W/g. Gravity's tie to Schwinger physics is a correspondence, not a lever: the Schwinger exponent is exactly a Boltzmann factor at twice the Unruh temperature, every astrophysical black hole is colder than 62 nK, and confinement is gravity's only contribution to fusion. About forty calculators with published test vectors belong in the crate's new `fusion` and `strong_field` modules and its existing `spacetime` module.

*How to read this report.* **[mainstream]** means the result agrees with standard nuclear physics, QED or GR and has been reproduced. **[contested]** means it is published but disputed, with no settled refutation. **[refuted]** means it has a named refutation or failed replication, or is internally inconsistent by many orders of magnitude. Numbers marked *(check)* were recomputed for this report with the five scripts in the notes folder, all run with `python3 -I`: [coulomb_checks.py](notes/coulomb_checks.py), [solar_pp_checks.py](notes/solar_pp_checks.py), [schwinger_checks.py](notes/schwinger_checks.py), [fringe_barrier_checks.py](notes/fringe_barrier_checks.py) and [gravity_unruh_hawking_checks.py](notes/gravity_unruh_hawking_checks.py). They are our inferences, not published claims. "(snippet)" marks a fact seen only in a search-engine summary, which is a weaker source.

## Every working fusion system tunnels at three to four times kT

### The barrier is 30–45 times the temperature of a burning plasma

Mainstream fusion physics starts from one factorisation. The cross section is σ(E) = S(E)/E · exp(−2πη), where η = Z₁Z₂α√(μc²/2E) is the Sommerfeld parameter and the S-factor varies slowly for non-resonant reactions ([Solar Fusion II](https://arxiv.org/abs/1004.2318); [Solar Fusion III](https://arxiv.org/abs/2405.06470)) **[mainstream]**. Bosch and Hale write the exponent as B_G/√E, with Gamow constant B_G = παZ₁Z₂√(2m_rc²), so the Gamow energy is E_G = B_G² ([UWFDM-1268](https://fti.neep.wisc.edu/fti.neep.wisc.edu/pdf/fdm1268.pdf)). From CODATA 2018 masses the Gamow constants match Bosch–Hale's published values to 10⁻⁵ for D–T (34.3826 against 34.3827 keV^½) and D–D. For D–³He they differ by 0.017% (68.7626 against 68.7508), probably because of older masses ([MFE Formulary notebook](https://notebook.community/MFEFormulary/MFEFormulary/ipython/BH_ReactionRates)) *(check)*. Code that must reproduce the Bosch–Hale tables bit for bit should therefore carry their tabulated constants rather than derive them.

The height of the barrier depends on where the nuclei are taken to touch, so no single published number exists. A textbook problem treats D and T as touching 2 fm spheres and gets **360 keV** ([esaral](https://www.esaral.com/q/consider-the-d-t-reaction-deuterium-tritium-fusion-56156)). A university lecture gives about 0.36 MeV, and notes that the D–T cross section peaks near 110 keV in the lab frame, well below the barrier ([Lehigh ME362](https://www.lehigh.edu/~eus204/teaching/ME362/lectures/lecture03.pdf), via search summary). Using R = r₀(A₁^⅓ + A₂^⅓), the D–T barrier is 444 keV at r₀ = 1.2 fm and 381 keV at 1.4 fm, and p–¹¹B's is 1.6–1.9 MeV *(check)*. Wikipedia's "the energy necessary to overcome the Coulomb barrier is 0.1 MeV" ([Nuclear fusion](https://en.wikipedia.org/wiki/Nuclear_fusion)) is roughly where σ_DT becomes large, not a barrier height.

Plasmas fuse anyway because of tunnelling. At 10 keV the classical turning point for two singly charged nuclei is 144 fm, about 40 nuclear radii out, and the D–T tunnelling factor exp(−2πη) is 1.9×10⁻⁵. Each step up in charge costs many orders of magnitude: at the same energy D–³He gives 3.6×10⁻¹⁰ and p–¹¹B gives 2.3×10⁻²¹ *(check)*. Multiplying the rising tunnelling factor by the falling Maxwell tail gives the Gamow peak, E₀ = (E_G(kT)²/4)^⅓, with 1/e width Δ = 4√(E₀kT/3) ([Solar Fusion II, Eqs. 8–9](https://arxiv.org/abs/1004.2318)). Rauscher's approximations, E₀ ≈ 0.1220(Z₁²Z₂²μT₉²)^⅓ MeV and Δ ≈ 0.2368(Z₁²Z₂²μT₉⁵)^⅙ MeV, agree with the exact form to within 0.3% ([Rauscher 2010](https://elog.gsi.de/esr/E127/190925_135951/PhysRevC.81.045807.pdf), via search summary) *(check)*.

| Reaction | E_G (keV) | Barrier at r₀ = 1.2 / 1.4 fm (keV) | exp(−2πη) at 10 keV | kT | Gamow peak E₀ | 1/e width | exp(−3E₀/kT) |
|---|---|---|---|---|---|---|---|
| D–T | 1182.2 | 444 / 381 | 1.90×10⁻⁵ | 10 keV | 30.9 keV | 40.6 keV | 9.4×10⁻⁵ |
| D–D | 985.8 | 476 / 408 | 4.88×10⁻⁵ | 10 keV | 29.1 keV | 39.4 keV | 1.6×10⁻⁴ |
| D–³He | 4728.3 | 888 / 761 | 3.60×10⁻¹⁰ | 10 keV | 49.1 keV | 51.2 keV | 4.0×10⁻⁷ |
| p–¹¹B | 22589.1 | 1861 / 1595 | 2.28×10⁻²¹ | 300 keV | 798 keV | 1130 keV | 3.4×10⁻⁴ |
| p+p, solar centre | 493.1 | 720 (2 fm spheres) | — | 1.34 keV | 6.04 keV | 6.56 keV | 1.3×10⁻⁶ |

All values are *(check)*. The table is the whole of the mainstream answer in one place. The reactions that matter come from the Maxwell tail at three to four times kT, suppressed by 10⁻⁴ for D–T at 10 keV and by 10⁻⁶ for p+p in the Sun, whose "∼6 keV" Gamow peak in Solar Fusion II reproduces as 6.04 keV.

### Bosch–Hale and Wurzel–Hsu reproduce to their last published digit

The Bosch–Hale (1992) reactivity is a seven-coefficient fit:

```
θ    = T / [1 − T(C₂ + T(C₄ + TC₆)) / (1 + T(C₃ + T(C₅ + TC₇)))]
ξ    = (B_G² / 4θ)^⅓
⟨σv⟩ = C₁ θ √(ξ / (m_rc² T³)) e^(−3ξ)
```

([UWFDM-1268 Eqs. 2–5](https://fti.neep.wisc.edu/fti.neep.wisc.edu/pdf/fdm1268.pdf); [MFE Formulary notebook](https://notebook.community/MFEFormulary/MFEFormulary/ipython/BH_ReactionRates)) **[mainstream]**. With the notebook's exact coefficients it reproduces all 28 entries of UWFDM-1268's Table 6 to within 0.3% *(check)*. D–T at 10 keV gives 1.1362×10⁻¹⁶ cm³/s against the tabulated 1.14×10⁻¹⁶. The notebook's printed D(d,n) endpoints come back exactly, 1.6572186×10⁻³⁷ m³/s at 0.1 keV and 2.6816742×10⁻²³ m³/s at 100 keV, and D–T peaks at 8.9465×10⁻¹⁶ cm³/s at 66.7 keV, Wikipedia's "about 70 keV". The fit is valid from 0.2 to 100 keV, and from 0.5 to 190 keV for ³He(d,p) ([UWFDM-1268 Table 7](https://fti.neep.wisc.edu/fti.neep.wisc.edu/pdf/fdm1268.pdf)); the notebook's wider "T range" column is actually the cross-section fit range. The older Duane-based NRL Formulary values differ from Bosch–Hale by −20% for D–T and −67% for D–³He at 1 keV *(check)*, so Bosch–Hale should be the reference.

The companion S-factor fit gives cross sections ([Wikipedia: Nuclear fusion, Bosch–Hale table](https://en.wikipedia.org/wiki/Nuclear_fusion)). Its D–T cross section peaks at 5.07 b at 64.8 keV in the centre-of-mass frame, which is 108 keV deuteron lab energy, the lecture's "near 110 keV" *(check)*. D(d,n) at 50 keV centre-of-mass gives 16.49 mb, within 1% of the pyroelectric paper's "16.3 mb at 100 keV" lab energy ([Geuther, Danon & Saglime, PRL 96, 054803](https://pages.github.rpi.edu/danony/html/Papers/PhysRevLett_96_054803%20Nuclear%20Reactions%20Induced%20by%20a%20Pyroelectric%20Accelerator.pdf)) *(check)*. The two fits also check each other: a numerical Maxwellian average of the cross-section fit matches the reactivity fit to 0.2–1.2% for D–T and D–D at 2, 10 and 50 keV, and to 0.6–3.1% for D–³He *(check)*. That gives a Rust crate an independent integration test.

The Lawson criterion reproduces end to end. Wurzel and Hsu balance alpha heating against bremsstrahlung:

```
nTτ_E = 3T² / [ (f_c + 1/Q_fuel) ⟨σv⟩ ε_F / 4 − C_B √T ]
C_B = 5.34×10⁻³⁷ W m³ keV^−½
```

They state that "when T < 4.3 keV … ignition is not possible" ([Wurzel & Hsu, Phys. Plasmas 29, 062103](https://arxiv.org/abs/2105.10954)) **[mainstream]**. With f_c = 3.5/17.6 and Bosch–Hale reactivity, their published values come back as follows *(check)*:

| Quantity | Published | Recomputed |
|---|---|---|
| Ignition (Q_fuel = ∞), minimum nTτ_E | 2.9×10²¹ keV s m⁻³ at 14 keV | 2.968×10²¹ at 14.25 keV |
| Ignition, minimum nτ_E | 1.6×10²⁰ s m⁻³ at 26 keV | 1.567×10²⁰ at 26.0 keV |
| Q_fuel = 1, minimum nTτ_E | 4.6×10²⁰ at 14 keV | 4.677×10²⁰ at 13.65 keV |
| Q_fuel = 1, minimum nτ_E | 2.5×10¹⁹ at 26 keV | 2.527×10¹⁹ at 25.7 keV |
| Ideal ignition temperature | 4.3 keV | 4.305 keV |

Every row agrees to the published two significant figures, which makes this an ideal integration test. Other fuels are far harder: catalysed D–D needs 1.1×10²³ at 52 keV, D–³He 5.2×10²² at 68 keV, and p–¹¹B has no entry at all.

p–¹¹B has a Gamow energy 19 times larger than D–T's. Tentori and Belloni (2023) refit its reactivity in the same Bosch–Hale form and raise it above Nevins and Swain (2000) by "12% at 100 keV and 50% at 500 keV" ([Nucl. Fusion 63, 086001](https://iopscience.iop.org/article/10.1088/1741-4326/acda4b)) **[mainstream]**. Both statements reproduce, as +12.6% and +50.1% *(check)*, but the 100 keV figure only works with a Nevins–Swain resonance term recalled from memory rather than from a source. At 300 and 500 keV that term is negligible, so +29.6% and +50.2% are the robust test values. The new fit reaches 5.62×10⁻²² m³/s at 500 keV, 63% of D–T's peak reactivity, but only at 7.5 times D–T's peak temperature *(check)*.

### NIF's gain of 4.13 counts laser light, not the wall plug

The 1997 National Academies definition of ignition is fusion energy exceeding laser energy ([Wurzel & Hsu](https://arxiv.org/abs/2105.10954)). NIF crossed that line on 5 December 2022 and had done so eleven times by June 2026 ([LLNL, Achieving Fusion Ignition](https://lmf.llnl.gov/science/achieving-fusion-ignition)) **[mainstream]**:

| Shot | Laser energy on target | Fusion yield | Target gain |
|---|---|---|---|
| 5 Dec 2022 | 2.05 MJ | 3.15 MJ | 1.54 *(check)* |
| 30 Jul 2023 | 2.05 MJ | 3.88 MJ | 1.89 *(check)* |
| 23 Feb 2025 | 2.05 MJ | 5.0 MJ | 2.44 |
| 7 Apr 2025 (record) | 2.08 MJ | **8.6 ± 0.45 MJ** | **4.13** |
| 1 Oct 2025 | 2.065 MJ | 3.5 ± 0.17 MJ | 1.74 |
| 20 Jun 2026 (11th) | ≈2.08 MJ (implied) | 7.9 ± 0.4 MJ | ≈3.8 |

LLNL's FY2025 annual report gives the record shot's peak power as 456 TW, and lists the October 2025 yield as 3.6 MJ, which conflicts with the 3.5 MJ above ([LLNL FY2025](https://annual.llnl.gov/fy-2025/national-ignition-facility-2025)). LLNL's own experts estimate that a test facility would need a target gain near 15 to reach zero net electric demand ([LLNL annual report](https://annual.llnl.gov/sites/annual/files/2026-02/llnlannual25-3-national-ignition-facility.pdf), via search summary). Ignition is real, and it got there through heat and compression rather than through the barrier; its gain is counted against the laser light reaching the target, not against the electricity drawn from the grid.

## Screening, muons and beams change rates by factors, not by the 10⁴⁰ cold fusion needs

### Bound electrons shave 20 eV off the barrier in gas and 300 eV in metals

Electrons around the target nucleus lower the barrier by a screening energy U_e, which raises the laboratory cross section by f ≈ exp(πηU_e/E) for U_e ≪ E ([Solar Fusion II, Eq. 15](https://arxiv.org/abs/1004.2318); [Solar Fusion III, Eq. 25](https://arxiv.org/abs/2405.06470)) **[mainstream]**. For d+d on a gas target the data "were compatible with the adiabatic limit (U_e = 20 eV)" ([Solar Fusion III](https://arxiv.org/abs/2405.06470)). Conference slides quote the Greife et al. (1995) gas value both as 25 ± 5 eV ([Czerski slides](https://comex5.ifj.edu.pl/slides/czerski.pdf), via search summary) and as 15 ± 5 eV, in "fair agreement with adiabatic limit (Ue = 20–30 eV)" ([Edinburgh slides](https://indico.ph.ed.ac.uk/event/373/contributions/3866/attachments/1856/2896/electron-screening.pdf)), so the gas d+d value is best treated as 15–25 eV. The effect follows directly *(check)*:

| Case | 1 keV | 2.5 keV | 5 keV | 10 keV |
|---|---|---|---|---|
| d+d, U_e = 25 eV | 1.47 | 1.10 | 1.036 | 1.012 |
| d+d, U_e = 300 eV (metal-like) | 47 (first order: 111) | 2.99 | 1.50 | 1.16 |
| d+d, U_e = 600 eV | 717 | 7.6 | 2.17 | 1.33 |
| d+³He, U_e = 219 eV | 653 (first order: 1862) | — | 1.92 | 1.26 |

The plain entries apply the shift exactly inside the Gamow factor, P(E+U_e)/P(E); the bracketed ones use the first-order exponential. The two differ by more than 2× once U_e/E exceeds a few percent, so a calculator must expose both.

Gas-target values already exceed theory for some reactions. ³He(d,p)⁴He gives U_e = 219 ± 15 eV against an adiabatic 119 eV ([Solar Fusion II](https://arxiv.org/abs/1004.2318)) **[contested]**, and ³He(³He,2p)⁴He fits 305 ± 90 eV ([Solar Fusion III](https://arxiv.org/abs/2405.06470)) **[contested]**. Measured values "tend to exceed theory by a factor ~2", and the candidate explanations are low-energy stopping powers, Trojan-Horse constraints on the bare cross section and clustering ([Solar Fusion II](https://arxiv.org/abs/1004.2318)) **[contested]**.

Metals give ten times more. Raiola et al. (2004) found Cu at 470 eV and Pt at 670 eV, "factor ~20–25 higher" than gas, while Nd and Hf stayed below 30 eV ([Edinburgh slides](https://indico.ph.ed.ac.uk/event/373/contributions/3866/attachments/1856/2896/electron-screening.pdf)). Huke et al. (2008) measured Ta at 322 ± 15, Zr at 297 ± 8, Pd at 313 ± 2 and Al at 190 ± 15 eV. They traced the earlier "unusual changes" mainly to deuteron density dynamics and surface contamination, and "clearly excluded" the Debye–Hückel model ([Huke et al., PRC 78, 015803](https://arxiv.org/abs/0805.4538)). Dielectric-function theory gives only 112–136 eV against measured values near 300 eV ([Czerski et al., EPL 68, 363](https://www.newenergytimes.com/v2/library/2004/2004Czerski-The2Hdp3HReaction.pdf), snippet), and a 2024 review concludes that "a suitable theory has not yet been found" ([Vesić 2024](https://preprints.org/manuscript/202409.0669/v1), snippet). The DOE's 2004 panel viewed these beam-on-foil studies "favorably" but called them "peripheral" to cold fusion ([DOE 2004 report](https://lenr-canr.org/acrobat/DOEreportofth.pdf)). The data are **[mainstream]**, reproduced by several groups; the interpretation is **[contested]**.

None of this matters in a burning plasma, where the Gamow peak sits at 30–40 keV and even 300 eV changes σ by a few percent *(check)*. In the Sun, Salpeter's weak-screening factor exp(0.188Z₁Z₂ζρ^½T₆^−³ᐟ²) gives **f₀ ≈ 1.05** for p+p, using an assumed core composition (ρ = 150 g/cm³, X = 0.34) *(check)*. Field-theory analyses expect only about 4% corrections to it, though "controversies have not entirely died out" ([Solar Fusion II](https://arxiv.org/abs/1004.2318); [Solar Fusion III](https://arxiv.org/abs/2405.06470)) **[mainstream]**. Screening is a big lever only for low-energy beams on cold targets.

### A muon catalyses 150 fusions but costs the energy of 270 to 340

A muon is 207 times heavier than an electron, so it binds d and t into a molecule small enough to fuse at about 10¹² s⁻¹ ([Koonin & Nauenberg preprint](https://newenergytimes.com/v2/library/1989/1989KooninNauenberg.shtml)). The best measured yield is **150 ± 20 d–t fusions per muon**, from LAMPF in 1986 ([MuFusE, arXiv:2606.05333](https://arxiv.org/abs/2606.05333); [Jones, Nature 321, 127](https://www.nature.com/articles/321127a0)) **[mainstream]**. The ceiling is X_μ = 1/(ω_s + λ₀/λ_c), where ω_s is the chance that the muon sticks to the alpha particle, λ₀ = 1/τ_μ and λ_c is the cycle rate. Initial sticking is 0.912%, and stripping lowers the effective value to 0.35–0.58% ([Pahlavani & Motevalli, Acta Phys. Pol. B 40, 319](https://www.actaphys.uj.edu.pl/R/40/2/319/pdf)); measured high-density sticking is "typically 10–50% lower" than standard theory ([MuFusE](https://arxiv.org/abs/2606.05333)). Kelly, Hart and Rose take ω_s = 0.47% and t_c/τ = 0.18%, which gives X = 154. Their simulation puts the cost at **4.70 GeV per muon**, against 5–6 GeV in the literature ([Kelly et al., J. Phys. Energy 3, 035003](https://iopscience.iop.org/article/10.1088/2515-7655/abfb4b)); their optimised Q = 1.87 is **[contested]**, a single-group design study.

The arithmetic is unforgiving *(check)*. At 17.6 MeV per fusion, energy breakeven needs 267, 284 and 341 fusions at 4.70, 5 and 6 GeV per muon, against Pahlavani's "about 285" at 5 GeV, and converting heat back to electricity at 40% raises those to 668–852. One hundred and fifty fusions release 2.64 GeV, which is 25 times the muon's rest energy but only about half its production cost, and the sticking floor caps the yield at 1/ω_s ≈ 213–333. Muon-catalysed fusion is therefore **[mainstream]** science without a power case, unless reactivation improves or the cost per muon halves.

### Beam devices fuse for real and return 10⁻⁶ to 10⁻¹⁴ of their input

Every small device that genuinely fuses is a beam hitting a target:

| Device | Driver | Fusion rate | P_fus / P_in *(check)* | Label |
|---|---|---|---|---|
| Pyroelectric crystal, UCLA 2005 | > 100 keV, > 4 nA (0.4 mW) | ~1000 s⁻¹ | 1.5×10⁻⁶ (beam only) | mainstream |
| Pyroelectric, RPI replication | 200 kV, 0.55 nA inferred | ≈10,500 n per cooling phase | — | mainstream |
| Fusor, TU Eindhoven | 43 kV × 30 mA (1.29 kW) | < 10⁶ n s⁻¹ | 9×10⁻¹⁰ | mainstream |
| Hirsch fusor, 1967 | up to 150 kV | ~10⁸ s⁻¹ | ~6×10⁻⁵ W against kW | mainstream |
| NASA Glenn "lattice confinement" | 2.9 MeV × 15 mA (43.5 kW) | 1.5 ± 0.3×10³ n s⁻¹ | 4×10⁻¹⁴ | contested |

The pyroelectric abstract reports "> 100 keV and > 4 nA" and "a neutron flux over 400 times the background level" ([RePEc](https://ideas.repec.org/a/nat/nature/v434y2005i7037d10.1038_nature03575.html), numbers from snippet), and Wikipedia gives "around 1000 fusion reactions per second" from a tip field of "about 25 gigavolts per meter" ([Pyroelectric fusion](https://en.wikipedia.org/wiki/Pyroelectric_fusion)). The other rows come from the RPI replication ([Geuther et al.](https://pages.github.rpi.edu/danony/html/Papers/PhysRevLett_96_054803%20Nuclear%20Reactions%20Induced%20by%20a%20Pyroelectric%20Accelerator.pdf)), the Eindhoven fusor ([FuseNet](https://fusenet.eu/node/361)), Hirsch ([Chemie.de lexicon](https://www.chemie.de/lexikon/Farnsworth-Hirsch-Fusor.html)) and NASA ([Steinetz et al., NASA/TP-20205001616](https://ntrs.nasa.gov/api/citations/20205001616/downloads/TP-20205001616.pdf)). In each case the fast ion carries 10–200 keV and its thick-target yield is set by ∫σ(E)/(dE/dx)dE, so the Bosch–Hale cross sections suffice to estimate it; the RPI paper matches its measured yield to first order this way. Fusors fail on their loss channels. Grid losses "can be at least five orders of magnitude higher than the energy released", and Rider showed that non-Maxwellian plasmas lose energy to thermalisation and bremsstrahlung faster than fusion can replace it ([Wikipedia: Fusor](https://en.wikipedia.org/wiki/Fusor)).

NASA Glenn's experiment is a chain of steps. Bremsstrahlung photons above 2.226 MeV split deuterons, and the photoneutrons heat lattice deuterons to an average of 64 keV, which then fuse. Samples of ErD₃ and TiD₂, each holding 5×10²⁴ deuterium atoms, gave 1.5–1.6×10³ n s⁻¹ of 2.45 MeV neutrons, with unfueled Er and Ti controls run under the same conditions ([NASA/TP-20205001616](https://ntrs.nasa.gov/api/citations/20205001616/downloads/TP-20205001616.pdf); theory in [Pines et al., PRC 101, 044609](https://journals.aps.org/prc/abstract/10.1103/PhysRevC.101.044609)). Photoneutron-driven d–d fusion is expected conventional physics. What stays **[contested]** is the predicted 1.2×10³ n s⁻¹, which relies on enhanced screening, and the reading of extra 4–5 MeV neutrons as Oppenheimer–Phillips stripping. Wikipedia notes that the reactions were not self-sustaining and the electron source "energetically expensive", and gives a beam endpoint of about 2.2 MeV that conflicts with the paper's 2.9 MeV ([Lattice confinement fusion](https://en.wikipedia.org/wiki/Lattice_confinement_fusion)). No independent replication was found.

### A mole of D₂ waits 10³² years for one fusion

Koonin and Nauenberg calculated the baseline that any room-temperature claim must beat. Ground-state molecules fuse at p+d ≈ 10⁻⁵⁵ s⁻¹, d+d ≈ **3×10⁻⁶⁴ s⁻¹** and d+t ≈ 10⁻⁶⁸·⁹ s⁻¹. Reaching the 1989 claims would require "enhancements of the electron mass by factors of 5–10", for which they knew no plausible mechanism ([Koonin & Nauenberg preprint](https://newenergytimes.com/v2/library/1989/1989KooninNauenberg.shtml)) **[mainstream]**. At that d+d rate a mole of D₂ produces one fusion per 1.7×10³² years, and the muonic molecule is 10⁷⁵·⁵ times faster *(check)*. Metzler, Hunt and Galvanetto put the enhancement needed for observable fusion at ">40 orders of magnitude"; the baseline is **[mainstream]**, while their proposal to cascade known mechanisms worth up to 30 orders each is **[contested]** ([arXiv:2208.07245](https://arxiv.org/abs/2208.07245)).

## The Schwinger field sits 1,454 times above the strongest laser field

### The vacuum switches on within a factor of ten of E_S

The QED critical field is E_S = m_e²c³/(eħ) = **1.32329×10¹⁸ V/m**, and its magnetic equivalent B_S = E_S/c = **4.41401×10⁹ T** *(check)*. Over one Compton wavelength, such a field does work 2m_ec² on an electron ([Di Piazza et al., RMP 84, 1177](https://arxiv.org/abs/1111.3886)) **[mainstream]**. Di Piazza et al. quote 1.3×10¹⁶ V/cm and I_cr ≈ 4.6×10²⁹ W/cm². Their arXiv text prints the magnetic field as "4.4×10¹⁴ G", a typo or extraction artefact for 4.414×10¹³ G *(check)*. The critical intensity is 4.65×10²⁹ W/cm² in the ε₀cE_S² convention and 2.32×10²⁹ W/cm² cycle-averaged (ε₀cE_S²/2), and that unstated factor of 2 explains many apparent disagreements in the literature.

Schwinger's result is a series, and only its first term counts pairs. Cohen and McGady separate the vacuum-decay rate w, the full series, from the mean pair-production rate Γ, "just the first term of that series" ([arXiv:0807.1117](https://arxiv.org/abs/0807.1117)). In SI units, for a constant field and x = E/E_S:

```
Γ(E) = (αε₀E²/π²ħ) exp(−π/x)                 (pairs m⁻³ s⁻¹)
w(E) = (αε₀E²/π²ħ) Σₙ n⁻² exp(−nπ/x)         (vacuum decay)
prefactor at x = 1: c/(4π³ƛ_C⁴) = 1.08704×10⁵⁶ m⁻³ s⁻¹
```

The only published anchor for the prefactor is Melissinos's "αE_c²/πħ = 3.45×10⁵⁰ cm⁻³ s⁻¹". It recomputes as 3.415×10⁵⁰ (ratio 0.990), and Schwinger's constant-field prefactor is exactly that number divided by π ([Melissinos, hep-ph/9805507](https://arxiv.org/abs/hep-ph/9805507)) *(check)*. The full sum exceeds the first term by 1.1% at x = 1 and by less than 10⁻⁶ for x ≤ 0.2; scalar QED has half the prefactor.

| E/E_S | exp(−πE_S/E) | Γ (m⁻³ s⁻¹) | Pairs in (1 µm)³ during 1 µm/c |
|---|---|---|---|
| 1 | 4.32×10⁻² | 4.70×10⁵⁴ | 1.6×10²² |
| 0.2 | 1.51×10⁻⁷ | 6.55×10⁴⁷ | 2.2×10¹⁵ |
| **0.1** | 2.27×10⁻¹⁴ | **2.47×10⁴⁰** | 8.2×10⁷ |
| **0.05** | 5.16×10⁻²⁸ | **1.40×10²⁶** | 4.7×10⁻⁷ |
| **0.01** | 3.65×10⁻¹³⁷ | **3.97×10⁻⁸⁵** | 1.3×10⁻¹¹⁷ |
| 6.9×10⁻⁴ (record laser, as if static) | — | 10⁻¹⁹³³ | — |
| 4.3×10⁻³ (Pais's 5.6×10¹⁵ V/m) | — | 10⁻²⁶⁸ | — |
| 7.6×10⁻¹⁰ (vacuum breakdown, 10⁹ V/m) | — | 10^(−1.8×10⁹) | — |

All rows are *(check)*. One pair in a cubic micron within one light-crossing time needs 0.064 E_S, equivalent to 9.6×10²⁶ W/cm² for a static field. The table assumes a uniform static field: a single plane wave has both field invariants zero and makes no pairs at all, as Melissinos notes. The general-field Nikishov form, written with the invariants, appears in the notes only from memory and needs checking against Fedotov et al. (2023) before it is coded.

### The record laser reaches a₀ ≈ 230; χ above 1 needs GeV electrons

The intensity record is still **(1.1 ± 0.2)×10²³ W/cm²**, set by CoReLS in Korea in 2021 over 80 consecutive shots at 0.1 Hz with a 1.1 µm FWHM spot ([Yoon et al., Optica 8, 630, IBS repository](https://pr.ibs.re.kr/handle/8788114/10375), snippet) **[mainstream]**. At an assumed 0.8 µm that is a₀ = 227 in the linear convention or 160 in the "E²/4π" convention, E₀/E_S = 6.9×10⁻⁴, and 2.1×10⁶ below I_S *(check)*. Peak power has risen since then, but peak intensity has not. ELI-Beamlines' L4 ATON delivered 5.1 PW (786 J in 154 fs) in September 2025 and plans "to increase power toward 10 PW in 2026" ([ELI news](https://www.eli-laser.eu/news/eli-s-l4-aton-laser-achieves-5-petawatt/)). ZEUS ran its first 2 PW experiment in May 2025 ([Michigan Engineering](https://news.engin.umich.edu/2025/05/the-us-has-a-new-most-powerful-laser), snippet), and Shanghai's 100 PW-class Station of Extreme Light is still under construction ([SHINE circular](https://sari.cas.cn/kjdt/xshd/202507/t20250703_7879010.html), snippet).

Two dimensionless numbers matter: the laser strength a₀ = eE₀/(m_ecω), and χ = E_rest/E_S, the field in the electron's rest frame in units of E_S. For a collision at angle θ, χ_e = (1 + β cos θ)γa₀ħω/(m_ec²), which reproduces LUXE's design values to within 1.6% ([Yap for LUXE, arXiv:2302.04699](https://arxiv.org/abs/2302.04699)) *(check)*. The European Strategy (ESPPU 2025) input tabulates indicative a₀ ("ξ") values for the main facilities ([Sarri et al., arXiv:2504.02608](https://arxiv.org/abs/2504.02608)):

| Facility | Energy, pulse length | ξ (ESPPU, indicative) | Notes |
|---|---|---|---|
| CoReLS | 80 J, 20 fs | 110 | holds the intensity record |
| ZEUS / Apollon | 75 / 70 J, 25 / 22 fs | 100 | — |
| ELI-NP | 200 J, 20 fs | 180 | — |
| ELI-Beamlines | 1500 J, 150 fs | 230 | reproduces only at λ ≈ 1.05 µm *(check)* |
| NSF OPAL | 2 × 500 J, 20 fs | 280 | design |
| SEL | 25 kJ, 25 fs | 1800 | ≈540 for a 100 PW (2.5 kJ) design *(check)* |
| LUXE phase 0 / 1 | 40 / 350 TW with a 16.5 GeV beam | < 7.9 / < 23.6 | χ_e < 1.5 / < 4.45; laser and detectors still need funding |
| E-320, FACET-II | about 10 TW with a 10 GeV beam | about 1–5 | χ ≳ 0.5, preliminary |

The table's ξ values reproduce to 0.90–0.98 under its stated focus assumptions *(check)*. SEL's 1,800 reproduces only from the tabulated 25 kJ, about 1 EW, not from its public 100 PW design, so its a₀ is better quoted as roughly 500–800 depending on convention. Because χ grows with γ, a 10 GeV electron meeting an 800 nm laser head-on gains 0.119 in χ per unit of a₀, so a ZEUS-class ξ of 100 gives χ ≈ 12 before any radiation losses on the way into the focus *(check)*. ESPPU's own example of 5 GeV against 10²³ W/cm², which "would result in … χ ≈ 10", recomputes as 9.1 or 12.8 depending on the a₀ convention.

The far end of the scale is unknown territory. The Ritus–Narozhny conjecture says perturbative QED breaks down when αχ^⅔ ≈ 1, that is at χ ≈ 1604, which at a₀ = 100 would need a 1.35 TeV electron *(check)*. Bubble-chain corrections do scale as αχ^⅔, but other examples show the expansion converging, and ESPPU concludes that "it is simply not known how QED behaves there" ([ESPPU input](https://arxiv.org/abs/2504.02608); [Mironov, Meuren & Fedotov, arXiv:2003.06909](https://arxiv.org/abs/2003.06909)) **[contested]**. At LUXE phase 1 the parameter is only 0.020.

### E144 and NA63 entered the regime in boosted frames; E-320 and LUXE extend it

| Experiment | Probe | χ reached | Result | Label |
|---|---|---|---|---|
| SLAC E144 (1997) | 46.6 GeV e⁻ + terawatt 527 nm laser | ≈0.3 (ξ = 0.36) | 106 ± 14 positrons above background; at least 4 laser photons per pair | mainstream; the "Schwinger" framing is contested |
| CERN NA63 (2018) | 178.2 GeV e⁺ in Si ⟨111⟩ | ≤ 1.4 | quantum radiation reaction needed to fit the spectra | mainstream |
| CERN NA63 (2023) | 200 GeV e⁻ in 400 µm Ge ⟨110⟩ | ≤ 2.4 | trident yield agrees with LCFA over ~3 orders | mainstream |
| SLAC E-320 (2024–25) | 10 GeV e⁻ + about 10 TW laser at ~30° | ≳ 0.5 | nonlinear Compton scattering; energy losses to ≥ 5 GeV | contested (preliminary) |
| STAR, RHIC (2021) | ultra-peripheral Au+Au at 200 GeV | — | 6085 exclusive e⁺e⁻ pairs; cos 4Δφ modulation of 16.8 ± 2.5% | measurement mainstream; interpretation contested |

E144 saw "106±14 positrons above background" ([OSTI 491608](https://www.osti.gov/biblio/491608), snippet). Melissinos explains this as a two-step process: a backscattered photon of about 29 GeV absorbs at least four more laser photons, so the yield "varies as η^2n with n = 5.1 ± 0.2". The threshold recomputes as 3.91 photons, which rounds up to four ([Melissinos](https://arxiv.org/abs/hep-ph/9805507)) *(check)*. His description of it as Schwinger pair production is **[contested]**: with ξ < 1 the process was perturbative multiphoton absorption, not tunnelling, and his Schwinger-type estimate of four pairs per shot compares with 0.1 observed. At NA63, "only a full quantum theory of radiation reaction is capable of explaining the experimental results" ([Wistisen et al., Nat. Commun. 9, 795](https://arxiv.org/abs/1704.01080)), and the trident study "probes values for the strong field parameter χ up to 2.4" ([Nielsen et al., PRL 130, 071601](https://arxiv.org/abs/2211.02390)). In the lab frame those crystal fields are only 5–8×10¹² V/m, 4–6×10⁻⁶ of E_S; Lorentz factors of 3.5–3.9×10⁵ supply the rest *(check)*.

E-320's abstract reports a peak a₀² ≳ 10 and χ ≳ 0.5 at about 30° and marks the analysis preliminary ([ExHILP 2025](https://indico.eli-laser.eu/event/198/contributions/935/)). The numbers are consistent only if a₀² is a cycle-averaged quantity, since a₀ = √10 gives χ = 0.35 and a₀ = √20 gives 0.50 *(check)*. STAR's modulation is 6.7σ from zero *(check)*, and the collaboration says it "has implications for vacuum birefringence" ([STAR, arXiv:1910.12400](https://arxiv.org/abs/1910.12400); [Brandenburg et al. review](https://arxiv.org/abs/2208.14943)). But the photons are quasi-real Weizsäcker–Williams photons, and the birefringence is inferred from the polarisation dependence of the Breit–Wheeler process rather than measured as a refractive index. Only crystals have published, peer-reviewed χ > 1. LUXE phase 1, at up to 4.45, would be the first laser experiment clearly above 1; its beamline is EU-funded, but "significant investment is still required" for the laser and detectors ([ESPPU input](https://arxiv.org/abs/2504.02608)).

### Vacuum birefringence is still 7 times below laboratory reach

The weak-field Euler–Heisenberg Lagrangian makes a magnetised vacuum birefringent, with Δn = n∥ − n⊥ = 3A_eB² and A_e = 2α²ħ³/(45μ₀m_e⁴c⁵) = **1.3247×10⁻²⁴ T⁻²**. That gives 2.48×10⁻²³ at 2.5 T ([Di Piazza et al.](https://arxiv.org/abs/1111.3886)) *(check)* **[mainstream theory]**. PVLAS spent 25 years on it, and its final result is Δn = **(12 ± 17)×10⁻²³** at 2.5 T against Δn(QED) = 2.5×10⁻²³ ([Ejlli et al., Physics Reports 871](https://arxiv.org/abs/2005.12913)). The uncertainty is 6.8 times the predicted signal, and the central value sits 0.56σ from QED *(check)*. There is no laboratory detection.

The astrophysical evidence depends on models **[contested]**. For the magnetar 4U 0142+61, IXPE measured polarisation rising from 14 ± 1% at 2–4 keV to 41 ± 7% at 5.5–8 keV, with a ~90° swing of the polarisation angle in between, which the authors say lends "further support" to vacuum birefringence ([Taverna et al., Science 378, 646](https://arxiv.org/abs/2205.08898)). Lai reads the swing as mode conversion at the QED vacuum resonance in the magnetar's atmosphere ([arXiv:2209.13640](https://arxiv.org/abs/2209.13640), snippet). For 1RXS J1708, the small emitting regions "do not allow a robust conclusion" ([Zane et al.](https://arxiv.org/abs/2301.12919), snippet), and the earlier optical claim for RX J1856.5−3754 rests on a polarisation of 16.43 ± 5.26% ([Mignani et al.](https://arxiv.org/abs/1610.08323)). Magnetar fields of 10¹⁴ G are 2.3 B_S, so the weak-field formula does not even apply there *(check)*, and a calculator must refuse B approaching B_S.

### Supercritical nuclei: the theory is sound, the signals were refuted

For a point nucleus the Dirac 1s level collapses at Z = 1/α = 137. For extended nuclei it reaches the negative-energy continuum at about Z = 173, and the vacuum should then become charged ([Rafelski et al., arXiv:1604.08690](https://arxiv.org/abs/1604.08690)) **[mainstream theory]**. Uranium–uranium collisions reach Z₁ + Z₂ = 184, but only briefly. A collision just below the Coulomb barrier lasts about 0.25×10⁻²⁰ s, while positron emission takes "typically 100 times longer", so only about 1% of the expected yield appears ([Rafelski et al.](https://arxiv.org/abs/1604.08690)).

GSI's EPOS and ORANGE groups reported narrow positron lines at about 250–400 keV from 1983 onwards ([Schweppe et al., PRL 51, 2261, GSI repository](https://repository.gsi.de/record/50234)). Reading them as a new 1.8 MeV particle was ruled out by e⁺e⁻ (Bhabha) scattering ([arXiv:nucl-ex/9706009](https://arxiv.org/abs/nucl-ex/9706009), snippet). Argonne's APEX found "no evidence … for narrow structures" ([OSTI 196563](https://www.osti.gov/biblio/196563), snippet), although a comment argued it was less sensitive to extremely narrow lines, and a later double-ORANGE run saw no lines within its statistics ([arXiv:nucl-ex/0004009](https://arxiv.org/abs/nucl-ex/0004009), snippet). The consensus attributes the early signals to "highly system dependent nuclear excitations converting into pairs" ([Rafelski et al.](https://arxiv.org/abs/1604.08690)), so the vacuum-decay or new-particle reading is **[refuted]**. Maltsev et al. propose a cleaner test, pair probability against impact parameter, in which the supercritical signature is a probability that rises as the collision energy falls. They call it "feasible with the CRYRING facility at GSI/FAIR" ([PRL 123, 113401](https://arxiv.org/abs/1903.08546)), but no such run was found scheduled by October 2026.

### Assisted tunnelling and graphene lower the exponent only on paper or in analogues

Schützhold, Gies and Dunne showed that adding a weak, fast field to a strong, slow one shortens the tunnelling path. Once the combined Keldysh parameter γ = mω/(qE) exceeds π/2, the instanton action becomes:

```
A ≈ (m²/qE) [ 2 arcsin(π/2γ) + (π/2γ²)√(4γ² − π²) ]
```

At γ = π the action "is reduced by about 40%" ([PRL 101, 130404](https://arxiv.org/abs/0807.0754)) **[mainstream theory, not observed in QED]**. The formula reproduces: the reduction at γ = π is 39.1%, and the tail approaches 2π/γ to within 4×10⁻⁵. At 0.1 E_S it would raise exp(−10π) = 2.3×10⁻¹⁴ to 4.9×10⁻⁹ *(check)*, but reaching γ = π there needs a 160 keV assisting photon overlapping a field more than a hundred times stronger than any laser has made *(check)*.

Condensed matter tests the prefactor rather than the exponent. Allor, Cohen and McGady give the graphene pair rate per unit area as π⁻²ħ^−³ᐟ²v_F^−½(qE)^³ᐟ² ([arXiv:0708.1471](https://arxiv.org/abs/0708.1471)). Berdyugin et al. observed Schwinger-like electron–hole generation above the critical current j_c ≈ nev_F, "at a rate ∝ E^(3/2)", and their quoted n ≈ 3×10¹⁰ cm⁻² recomputes as 3.12×10¹⁰ ([Science 375, 430](https://arxiv.org/abs/2106.12609)) *(check)* **[mainstream as an analogue]**. Schmitt et al. measured a universal one-dimensional Schwinger conductance in ballistic graphene transistors, 0.18 ± 0.02 mS against a predicted 0.186 mS. Their effective critical field is ≲6×10⁷ V/m for a 0.2 eV gap, but that gap rests on an ansatz with "no full theory" ([Nat. Phys. 19, 830](https://arxiv.org/abs/2207.13400)) **[mainstream analogue; contested in detail]**. These devices mimic a QED effect "predicted to occur at enormous fields of 10¹⁸ V m⁻¹" without touching the QED vacuum.

## Fringe claims miss their mainstream bounds by 300 to 10⁵³

### Pais's fusion application was abandoned, and falls 10⁶–10⁹ short on its own numbers

US application 15/928,703, "Plasma Compression Fusion Device", was filed by the Navy on 22 March 2018 and published on 26 September 2019. It was never granted: the USPTO mailed non-final actions on 16 August and 11 December 2020, and the application was abandoned on 14 July 2021 "for failure to respond to an office action" ([US20190295733A1](https://patents.google.com/patent/US20190295733A1/en)). Press reports that the USPTO "has issued a patent" confused pre-grant publication with grant ([Atomic Insights](https://atomicinsights.com/?p=17049)). The companion IEEE Transactions on Plasma Science paper is a design study, with no experimental results found ([IEEE Xplore](https://ieeexplore.ieee.org/abstract/document/8871349/authors), snippet), and the words "Schwinger", "vacuum polarization" and "Pais effect" do not appear in the fusion application. Its stated parameters are a density of about 10²⁰ m⁻³, D–D at about 232 million °C, "at least one Coulomb" per dynamic fusor, a field B_MAX of order 10⁶ T, an energy gain of order 10¹⁸ "derived from Equation 2", and gigawatt-to-terawatt output from kW–MW input ([US20190295733A1](https://patents.google.com/patent/US20190295733A1/en)).

Every one of those fails on its own terms *(check)*. One coulomb on a 0.15 m-radius sphere gives a surface field of 4.0×10¹¹ V/m, which is 400 times the 10⁹ V/m vacuum limit Pais himself states ([US 10,135,366](https://patents.google.com/patent/US10135366B2/en)), and stores 3.0×10¹⁰ J; even a 1 m sphere sits 9× above that limit. Spinning that charge at 10⁴ rad/s makes 4.4 mT inside the shell, 2.2×10⁸ times below the claimed field; the patent's eq. 11 reaches 10⁶ T only through an unbounded (ωt)² factor after about 1.2 s of operation. A 10⁶ T field would push at 4×10¹⁷ Pa and hold 5.6×10¹⁵ J (1.3 Mt TNT) in a 0.3 m sphere, all from kW–MW of input. "Gain ~ B³" yields 10¹⁸ only if the gain is taken to be (B/1 T)³, and at the 30 T the patent itself cites the same rule gives 2.7×10⁴. Finally, 232 million °C is 20.0 keV, where the Bosch–Hale D–D reactivity is 5.0×10⁻²⁴ m³/s. At 10²⁰ m⁻³ that yields 1.46×10⁴ W/m³: **207 W in a 0.3 m sphere, 5×10⁶ to 5×10⁹ times below "GW–TW"**. The notes' constant-S estimate of 161 W was 29% low. One gigawatt in that sphere would need n = 2.2×10²³ m⁻³ held for seconds, and would emit about 8.5×10²⁰ neutrons per second.

The "vacuum polarization" side of Pais's work fails by more. His craft patent's surface field, σ/ε₀ = **5.6×10¹⁵ V/m** at σ = 5×10⁴ C/m² ([US 10,144,532](https://patents.google.com/patent/US10144532B2/en)), is 5.6×10⁶ times his own breakdown limit and still 240 times below E_S. At that field the Euler–Heisenberg correction to the field energy, (α/45π)(E/E_S)², is 9×10⁻¹⁰, and the Schwinger rate is 10⁻²⁶⁸ m⁻³ s⁻¹ *(check)*. The craft patent also sets ε₀E_S² equal to ħν⁴/(8π²c³), the zero-point energy density cut off at ν, and obtains ν ≈ 4.2×10²¹ Hz. That is a dimensional identity: ε₀E_S²/2 is exactly π/α = 430.5 times the zero-point density cut off at the Compton frequency, so the equation only returns ν = 5.4 ω_C. It is numerology, not a condition for pair creation *(check)*. The empirical record closes the matter. NAWCAD spent **$508,000** testing the related field generator from 2016 to 2019, and the "Pais Effect" "could not be proven" ([TWZ, NAWCAD statement](https://www.twz.com/39012/the-navy-finally-speaks-up-about-its-bizarre-ufo-patent-experiments)). No test of the fusion device is on record **[refuted]**.

### Fleischmann and Pons lack 37 million times the neutrons they reported

Fleischmann and Pons's 1989 Table 1 reports 1.751 W of excess heat from a 0.4×10 cm Pd rod at 64 mA/cm², and 26.8 W at 512 mA/cm² rescaled from a 1.25 cm electrode. They give a neutron flux "of the order 4×10⁴ s⁻¹" and tritium "to the extent of 1–2×10⁴ atoms s⁻¹", and say themselves that the heat "would require rates … in the range 10¹¹–10¹⁴" ([F–P 1989, NET copy](https://www.newenergytimes.com/v2/library/1989/1989Fleischmann-PrelimNote-SK-MARKUP.pdf)). Wikipedia's "about 4,000 neutrons per second" conflicts with the primary text ([Wikipedia: Cold fusion](https://en.wikipedia.org/wiki/Cold_fusion)). D–D fusion releases a mean 3.651 MeV, so 1 W corresponds to 1.71×10¹² reactions and 8.55×10¹¹ neutrons per second. The 1.751 W cell should therefore have emitted **1.50×10¹² n/s, 3.7×10⁷ times the claimed rate**. That is a flux of 1.2×10⁷ n cm⁻² s⁻¹ at 1 m, roughly 17 Sv/h with a recalled dose coefficient, or a lethal dose within half an hour. The heat needs a rate per D pair of 3.5×10⁻¹¹ s⁻¹, 1.2×10⁵³ times the Koonin–Nauenberg rate *(check)*.

The DOE 2004 review adds that the ⁴He + γ channel occurs "approximately only once for every 10⁷ D + D fusion reactions", and that the claim that its 24 MeV goes into lattice heat without gammas raised "serious concerns" ([DOE 2004 report](https://lenr-canr.org/acrobat/DOEreportofth.pdf)). Petrasso et al. found the Utah gamma peak narrower than the detector's resolution ([AIP press release](https://newenergytimes.com/v2/sr/ColdFusion1989/APS-Baltimore/19890515-AIP-Press-Release.shtml), snippet). D–D fusion as the source of the heat is **[refuted]**; ⁴He as gamma-free ash is **[contested]** at best. Both federal reviews reached the same place. In 1989 DOE accepted the ERAB panel's "skepticism regarding the scientific validity of cold fusion" and saw "no need to establish special cold fusion programs" ([DOE/ERAB transmittal](https://documents.theblackvault.com/documents/doe/osti-coldfusion/ERAB.pdf)). In 2004, eighteen reviewers split "approximately evenly" on excess heat. Most found the effects "not repeatable", with a magnitude that "has not increased in over a decade", and two-thirds did not find the nuclear evidence conclusive. ⁴He appeared in 5 of 16 heat-producing cells, "very close to … background", and the conclusions were "similar to those found in the 1989 review" ([DOE 2004 report](https://lenr-canr.org/acrobat/DOEreportofth.pdf)).

### Google, ARPA-E and DARPA turned cold fusion into beam-target physics

Google funded about 30 researchers with $10 million from 2015 to 2018 ([Physics World](https://physicsworld.com/a/cold-case-of-cold-fusion-reopened-by-google-funded-scientists)), producing 28 peer-reviewed papers ([IEEE Spectrum](https://spectrum.ieee.org/google-nuclear-energy)). Palladium exceeded the claimed 0.875 loading threshold "only twice" in several dozen tests. None of 420 nickel-powder samples "provided evidence of excess heat", and weeks of pulsed-plasma runs saw no enhanced tritium. The team "found no evidence of anomalous effects claimed by proponents of cold fusion" (Berlinguette et al., Nature 570, 45, 2019, via [Physics World](https://physicsworld.com/a/cold-case-of-cold-fusion-reopened-by-google-funded-scientists)) **[mainstream null result]**.

ARPA-E's 2023 exploratory topic made eight awards totalling **$9,956,549**, from $295,924 to $2,000,000, to Stanford, MIT, Amphionic, the Energetics Technology Center, LBNL, Texas Tech and two Michigan teams. Its director framed the question as "does this area show promise … Or can we conclusively show that it does not?" ([New Energy Times listing](https://news.newenergytimes.net/2023/03/01/u-s-department-of-energys-arpa-e-funds-low-energy-nuclear-reaction-research/); [ARPA-E FY2023 report](https://arpa-e.energy.gov/sites/default/files/2025-09/ARPA-E%20FY%202023%20Annual%20Report.pdf), snippet). No published confirmation of excess heat was found from any team. The visible outputs are measurements. LBNL and UC Davis report sub-2.5 keV D–D fusion in Pd and Ti hydrides "enhanced by over 10¹⁸ compared to theoretical bare-nucleus fusion rates", at rates they call "low" ([Karahadian et al., arXiv:2512.06212](https://arxiv.org/abs/2512.06212v1)) **[contested, not peer reviewed]**. Texas Tech showed that untreated stainless steel 316 holds 15–240 pmol of ⁴He, so helium as nuclear ash needs pretreated materials ([arXiv:2405.02307](https://arxiv.org/abs/2405.02307)) **[mainstream]**. At UBC, TRIUMF and LBNL, electrochemical loading in the "Thunderbird Reactor" raised beam-driven D–D rates in Pd by **15(2)%** ([Chen et al., Nature 644, 640](https://ideas.repec.org/a/nat/nature/v644y2025i8077d10.1038_s41586-025-09042-7.html)), in a device that returns "less than about 10⁻¹⁰" of the energy it consumes ([Physics World](https://physicsworld.com/?p=123517)) **[mainstream]**. The sources conflict on its beam energy: DARPA's slide says 15 keV, while Physics World's "megaelectronvolt" probably refers to the products.

DARPA's MARRS programme, with a proposers' day on 20 January 2026, builds on those results. It aims to raise solid-state rates "from ~0.1 (state-of-the-art) to 10⁶ fusion reactions/second/gram": 10³ at the end of an 18-month phase 1 and 10⁶ at the end of phase 2, with a design target of "Ue≥600 eV" ([DARPA MARRS slides](https://www.darpa.mil/sites/default/files/attachment/2026-01/marrs-program-proposers-day-presentation.pdf)). The slides' own detection floors are about 10⁹ s⁻¹ g⁻¹ by calorimetry and 0.1 s⁻¹ g⁻¹ by neutrons. In watts the state of the art is 5.8×10⁻¹⁴ W/g, and the phase 2 goal of 5.8×10⁻⁷ W/g is still 10³ below calorimetric visibility. Fleischmann and Pons's 1.39 W/cm³ (0.116 W/g) is 2×10¹² times the state of the art and 2×10⁵ times even the phase 2 goal *(check)*. MARRS is mainstream research resting on a **[contested]** premise, that screening of 600 eV or more can be engineered, and even complete success would leave the 1989 heat claims unsupported by 5–12 orders of magnitude.

Europe's CleanHME (grant 951974, €5.6 million, 17 partners, 2020–2025) claimed that "the total heat energy produced exceeds by many orders of magnitude the chemical energy" ([Polish Science](https://www.polishscience.pl/en/?p=141620), snippet), but CORDIS lists no results ([CORDIS 951974](https://cordis.europa.eu/project/id/951974)) **[contested]**. A €4 million companion, HERMES, is reported from the same funding round ([Science|Business](https://sciencebusiness.net/news/projects-touting-potential-new-hydrogen-energy-source-among-winners-latest-eu-research-awards), snippet). The Navy SPAWAR group attributed triple tracks in CR-39 to neutrons of at least 9.6 MeV ([Mosier-Boss et al. 2009](https://newenergytimes.com/v2/library/2008/2008BossTripleTracks.pdf), snippet), which Wikipedia notes cannot be validated without quantitative neutron analysis ([Wikipedia: Cold fusion](https://en.wikipedia.org/wiki/Cold_fusion)) **[contested]**.

### Widom–Larsen, Rossi, hydrinos and bubble fusion fail their own arithmetic

Widom and Larsen propose that "heavy" surface electrons convert protons into slow neutrons, because "No Coulomb barriers exist for the weak interaction neutron production" ([cond-mat/0505026](https://arxiv.org/abs/cond-mat/0505026)). Ciuchi et al. redid the rate. For β = m*/m_e = 20, the claimed 0.56 Hz becomes 1.8×10⁻³ Hz, "about 300 times smaller"; the claimed value would need β = 61, and Widom and Larsen's own estimator gives only β = 1.01–1.8, "below threshold" ([Ciuchi et al., EPJC 72, 2193](https://arxiv.org/abs/1209.6501)). Their eq. 5 reproduces as 0.564 Hz, a ratio of 313 *(check)* **[refuted]**.

Rossi's 2014 Lugano report claimed a COP of "about 3.2 and 3.6" over 32 days, 1.6 MWh net, and ⁶²Ni enrichment from 3.6% to 98.7%. Temperatures came from an infrared camera with literature emissivity, and Rossi himself handled the fuel ([Lugano report](https://www.mech.kth.se/~hanno/Lugano_report.pdf)). A replica calibrated by a pro-LENR group found COP 1.13 and 1.32, with true surface temperatures 300–360 °C lower than reported, and concluded "at least an order of magnitude lower excess heat" ([Valat et al., JCMNS 21, 81](https://jcmns.org/article/72418-investigations-of-the-lugano-hotcat-reactor/attachment/150447.pdf)). The claimed 5825 MJ from 1 g of fuel is 1.3×10⁵ times petrol, and burning every ⁷Li nucleus in the charge would supply only 45% of it *(check)*. Industrial Heat said it could not reproduce the results, and the resulting lawsuit settled in 2017 ([Wikipedia: Energy Catalyzer](https://en.wikipedia.org/wiki/Energy_Catalyzer)). The COP of 3.2–3.6 is **[refuted]**, the residual excess **[contested]**, and the 2019 E-Cat SK **[refuted/unvalidated]**.

Mills's hydrino states below the hydrogen ground state "are not square integrable", his wave equation is not Lorentz invariant, and there is "no theoretical support of the hydrino hypothesis" ([Rathke, New J. Phys. 7, 127](https://ar5iv.labs.arxiv.org/html/quant-ph/0505150)) **[refuted]**. The q = ½ state would bind at 54.4 eV, above Rathke's quoted stability bound of about 20 eV *(check)*. A company-associated 300 kW SunCell report has no independent peer-reviewed validation ([peak-oil.org summary](https://news.peak-oil.org/2020/01/validation-report-of-300000w-of-power-produced-by-the-suncell-brilliant-light-power/), snippet) **[contested]**.

Taleyarkhan et al. (Science 2002) claimed 4–8×10⁴ n/s from cavitated, chilled deuterated acetone. Their own tritium implied about 7×10⁵ n/s, 9–17 times the measured neutrons, where D–D predicts about 1:1 ([Taleyarkhan et al., NET copy](https://www.newenergytimes.com/v2/library/2002/2002-TAL-EvidenceForNuclearEmissions.shtml)). ORNL's Shapira and Saltmarsh found any correlated emission "at least 4 orders of magnitude too small" to explain the tritium ([arXiv:1307.3217](https://arxiv.org/pdf/1307.3217), snippet), and Naranjo matched the 2006 spectrum to ²⁵²Cf. Purdue substantiated misconduct on 2 of 12 allegations without ruling on the results' validity, and ONR debarred Taleyarkhan for 28 months ([Wikipedia: Bubble fusion](https://en.wikipedia.org/wiki/Bubble_fusion)). Even if true, the claimed power was 5×10⁻⁸ to 8×10⁻⁷ W *(check)* **[refuted]**. PLOS ONE retracted Holmlid's paper on laser-induced "mesons" in ultra-dense hydrogen, citing conflict with baryon-number conservation and decay constants that matched electronic pickup ([PMC retraction notice](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6386439/), snippet) **[refuted]**.

### Metal screening cannot be carried to room temperature

The one bridge proponents build from real data to cold fusion is the 300–700 eV screening measured in metals, and it collapses under its own arithmetic. Take a naive thermal model of Pd–D at 300 K, σ = (S/E)·exp(−√(E_G/(E + U_e))) with D/Pd = 1. Its recalled S-factors (53.7 and 55.6 keV b) match the Bosch–Hale A₁ coefficients in the mainstream table. The cross sections come out as follows *(check)*:

| U_e | σ at 300 K |
|---|---|
| 0 | 10⁻²⁶⁷⁵ b |
| 27 eV (gas) | 10⁻⁷⁶ b |
| 100 eV | 10⁻³⁶·⁵ b |
| 300 eV | 10⁻¹⁸·³ b |
| 600 eV | 10⁻¹¹ b |

The model reproduces MARRS's theory slide, "≥1 Watt/gram for 10⁻¹³ barn", giving 2.5 W/g. To match MARRS's measured state of the art it needs U_e = 168 eV; its phase 2 goal needs 270 eV and Fleischmann–Pons's heat 423 eV. Applied unchanged at thermal energies, the measured 300–700 eV would make palladium deuteride emit 2×10⁷ to 9×10¹⁵ reactions per second per gram, and neutron counters that see 0.1 would see that. So either keV-energy screening does not carry over to thermal deuterons as a constant shift, or the gas model is wrong for localised deuterons 2.9 Å apart. Either way, using the screening anomaly to "explain" cold-fusion heat is **[refuted]** by the neutron-rate bound. The Debye model behind "~300 eV" reproduces numerically (295 eV), but only with a screening radius of 0.05 Å that no degenerate metal supports, so it is **[refuted]** as an explanation, consistent with Huke's exclusion ([Huke et al.](https://arxiv.org/abs/0805.4538)).

The same arithmetic deflates the "10¹⁸" headline. An enhancement of 10¹⁸ over bare nuclei needs U_e of only 51 eV at a centre-of-mass energy of 50 eV, 195 eV at 100 eV and 990 eV at 200 eV *(check)*, so the figure mainly measures how small bare sub-keV cross sections are. Szczecin's ZrD₂ data fall 7 orders of magnitude and then plateau below about 2.5 keV, with U_e = 340 eV against "about 100 eV" for defect-free material, which the group calls the "first direct observation" of thermal fusion in ion tracks ([Czerski et al., arXiv:2409.02112](https://arxiv.org/abs/2409.02112v2)) **[contested]**.

### The ledger

| Claim | Published | Reproduces? | Mainstream bound | Gap | Label |
|---|---|---|---|---|---|
| Pais fusor power | GW–TW at 10²⁰ m⁻³, 20 keV | no derivation given | Bosch–Hale D–D: 207 W | 5×10⁶–5×10⁹ | refuted |
| Pais 1 C fusor | ≥ 1 C | — | his own 10⁹ V/m breakdown | 400× (4×10¹¹ V/m) | refuted |
| Pais 10⁶ T | eq. 11, μ₀σRω³t² | only via unbounded t² | spinning shell: 4.4 mT | 2.2×10⁸ | refuted |
| Pais vacuum polarisation | ε₀E_S² = ħν⁴/8π²c³ | yes; an identity (π/α) | Schwinger rate 10⁻²⁶⁸ at 5.6×10¹⁵ V/m | 240× in field | refuted |
| F–P heat as D–D | 1.751 W, 4×10⁴ n/s | — | 1.5×10¹² n/s expected | 3.7×10⁷ | refuted |
| F–P rate per pair | 10⁻¹⁹ s⁻¹ (neutrons) | 9.4×10⁻¹⁹ | Koonin–Nauenberg 3×10⁻⁶⁴ s⁻¹ | 10⁴⁵ (heat: 10⁵³) | refuted |
| ⁴He as gamma-free ash | — | — | ⁴He + γ in 1 of 10⁷ reactions | — | contested at best |
| Google 2015–2018 | no anomalies | — | — | — | mainstream (null) |
| MARRS goals | 0.1 → 10⁶ s⁻¹ g⁻¹ | 5.8×10⁻¹⁴ → 5.8×10⁻⁷ W/g | calorimetry floor 10⁹ s⁻¹ g⁻¹ | phase 2 is 10³ below it | contested premise |
| Karahadian 2025 | > 10¹⁸ over bare | — | needs U_e of only 51–990 eV | — | contested |
| Screening explains heat | 300–700 eV | — | would give 2×10⁷–9×10¹⁵ s⁻¹ g⁻¹, against 0.1 seen | 10⁸–10¹⁷ | refuted |
| Debye screening | ~300 eV | yes, 295 eV | R_D = 0.05 Å is unphysical | — | refuted |
| NASA LCF | 1.5×10³ n/s | prediction 1.2×10³ (with screening) | P_fus/P_beam = 4×10⁻¹⁴ | — | contested |
| Widom–Larsen | 0.56 Hz at β = 20 | yes, 0.564 | Ciuchi: 1.8×10⁻³ Hz | 313× | refuted |
| Lugano COP | 3.2–3.6 | yes, 3.13 per file | replica: COP 1.13–1.32 | 16× in net output | refuted (residual contested) |
| Hydrino | W₁/q² | 54.4 eV at q = ½ | not square integrable | — | refuted |
| Bubble fusion | 4–8×10⁴ n/s plus tritium | — | tritium implies 7×10⁵ n/s; ORNL ≥ 10⁴ | 9–17× internally | refuted |
| SPAWAR CR-39 | ≥ 9.6 MeV neutrons | — | no quantitative neutron analysis | — | contested |
| CleanHME heat | "many orders" above chemistry | — | no CORDIS results | — | contested |
| Holmlid H(0) mesons | — | — | baryon number; pickup | — | refuted (retracted) |
| E144 as "Schwinger pairs" | — | — | ξ < 1, perturbative multiphoton | — | contested framing |
| STAR vacuum birefringence | 16.8 ± 2.5% | 6.7σ | inferred from polarisation, not Δn | — | contested interpretation |
| GSI positron lines | 250–400 keV peaks | — | APEX, ORANGE null; Bhabha scattering | — | refuted |

## Gravity supplies confinement and a temperature, never a barrier breaker

### Schwinger's exponent is a Boltzmann factor at twice the Unruh temperature

A uniformly accelerated observer sees a thermal bath at T_U = ħa/(2πck_B), about 4.06×10⁻²¹ K per m/s², so 1 K needs 2.47×10²⁰ m/s² ([Wikipedia: Unruh effect](https://en.wikipedia.org/wiki/Unruh_effect); [Crispino, Higuchi & Matsas, RMP 80, 787](https://arxiv.org/abs/0710.5373)) **[mainstream]**. The standard review states that the effect "does not need experimental confirmation any more than free Quantum Field Theory does" ([Crispino et al.](https://arxiv.org/abs/0710.5373)). Set a = eE/m_e and the link to Schwinger is exact: πE_S/E ≡ m_ec²/(2k_BT_U), an identity that holds to 10⁻¹² *(check)*. The flat-space pair rate is thus a Boltzmann factor for one rest mass at **2T_U**, not at T_U. At E_S itself the electron's acceleration is mc³/ħ = 2.3×10²⁹ m/s², and its Unruh temperature is m_ec²/(2πk_B) = 9.4×10⁸ K, or 81 keV *(check)*.

Cai and Kim find that the Schwinger effect in de Sitter and anti-de Sitter space reduces to the flat formula "with T=2T_U", through an effective temperature T_dS = √(T_GH² + T_U²) + T_U ([arXiv:1407.4569](https://ar5iv.labs.arxiv.org/html/1407.4569)). Its first term matches the √(a² + H²)/2π temperature that an accelerated detector sees in de Sitter space ([Kim, arXiv:1512.03142](https://arxiv.org/abs/1512.03142)) **[mainstream]**. The cleanest published bridge between the effects is Nikishov and Ritus's detailed balance for photon emission in a constant field at β⁻¹ = a/2π ([Crispino et al.](https://arxiv.org/abs/0710.5373)) **[mainstream as a result, contested in interpretation]**, and Parentani and Massar connect Schwinger, Unruh and accelerated black holes ([INSPIRE 416584](https://inspirehep.net/literature/416584)) **[mainstream]**. The factor of 2 itself remains open. Volovik calls it a "factor of 2 problem" and argues that coherence between processes doubles the Unruh, de Sitter and Hawking temperatures alike ([arXiv:2206.02799](https://arxiv.org/abs/2206.02799); [arXiv:2409.03316](https://arxiv.org/abs/2409.03316)) **[contested]**, while Labun and Rafelski conjecture that pair production is captured by a g = 1 vacuum-fluctuation temperature ([arXiv:1203.6148](https://arxiv.org/abs/1203.6148)) **[contested]**. "Gravity as Schwinger" is therefore a correspondence between a Euclidean periodicity and an exponent, not a mechanism by which gravitational fields make pairs at useful rates. A calculator should expose `unruh_temperature` and a separate Schwinger effective temperature of 2T_U, and test the identity, so that the factor of 2 stays explicit.

### Every gravitational horizon is unmeasurably cold

One formula, T = ħκ/(2πk_B), with κ the surface gravity expressed as a rate, covers every case. Only charged particles in strong electromagnetic fields reach temperatures anyone could measure:

| Situation | Acceleration or surface gravity | Temperature | Source |
|---|---|---|---|
| Cosmological horizon, H₀ = 67.4 km/s/Mpc | H₀ = 2.18×10⁻¹⁸ s⁻¹ | 2.66×10⁻³⁰ K | check; [Planck 2018](https://arxiv.org/abs/1807.06209) |
| Earth's surface gravity (Unruh) | 9.81 m/s² | 4.0×10⁻²⁰ K | check |
| Sun's surface gravity (Unruh) | 274 m/s² | 1.1×10⁻¹⁸ K | check |
| Black hole of 1 M☉ (Hawking) | — | **61.7 nK** | check |
| BEC analogue black hole | κ = 288 s⁻¹ (check) | 0.35 nK | [Steinhauer 2019](https://arxiv.org/abs/1809.00913) |
| Electron in a 10 MV/m linac | 1.76×10¹⁸ m/s² | 7.1 mK | check |
| Penning-trap proposal (Rogers 1988) | 6×10¹⁹ g | 2.4 K (published "≈ 2 K") | [Crispino et al.](https://arxiv.org/abs/0710.5373) |
| 45.6 GeV electrons in LEP | 2.3×10²³ m/s² | ≈ 940 K | check |
| Laser-driven electron, a₀ = 100 (Chen–Tajima) | 1.2×10²⁶ m/s² | 4.9×10⁵ K | check |
| Electron at E_S | 2.3×10²⁹ m/s² | 9.4×10⁸ K (81 keV) | check |
| Primordial black hole evaporating now, 5.00×10¹¹ kg | — | 2.45×10¹¹ K (21 MeV) | check; [Wikipedia: Hawking radiation](https://en.wikipedia.org/wiki/Hawking_radiation) |

Wikipedia's "10⁻⁷ K" for a 3 M☉ hole is about 5 times too high, since the formula gives 20.6 nK ([Wikipedia: Hawking radiation](https://en.wikipedia.org/wiki/Hawking_radiation)), and a hole of 4.5×10²² kg, about 0.6 lunar masses, sits in equilibrium with the CMB, so heavier holes currently grow ([Fixsen 2009](https://arxiv.org/abs/0911.1955)) *(check)*. The photon-only evaporation time, 2.095×10⁶⁷ yr × (M/M☉)³, sits 2% below Wikipedia's 2.140×10⁶⁷ for reasons not identified. Page's total power of "2×10⁻⁴ℏc⁶G⁻²M⁻²" gives a lifetime 9.65 times shorter, because neutrinos dominate the emission ([Page, PRD 13, 198](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.13.198)) *(check)*. At the 21 MeV temperature of a hole evaporating today, the massless-only Page lifetime is 34.5 Gyr, so emission of muons and pions must be what brings it down to the age of the universe *(check)*. No Hawking radiation has been detected from astrophysical or primordial black holes ([Wikipedia: Hawking radiation](https://en.wikipedia.org/wiki/Hawking_radiation)) **[mainstream status]**. The cosmological horizon radiates at ħH/(2πk_B) ([Gibbons & Hawking, PRD 15, 2738](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.15.2738)), colder than anything else in the table by ten orders of magnitude.

Charged black holes are where gravity and the Schwinger field actually meet. The horizon field of an extremal Reissner–Nordström hole equals E_S at 1.06×10³⁶ kg, or 5.3×10⁵ M☉, so lighter extremal holes would have supercritical horizons and discharge by pair creation. That is consistent with Gibbons's statement that holes "more massive than 10¹⁷ g" lose charge by "a Schwinger type formula" ([Gibbons, CMP 44, 245](https://link.springer.com/doi/10.1007/BF01609829), abstract via search) *(check; our derivation)*. A horizon hot enough to emit electron pairs thermally needs T_H ≈ 9.4×10⁸ K, a black hole of about 1.3×10¹⁴ kg, and no laboratory gravitational field comes within 10²⁸ of the required acceleration *(check)*.

### Laboratory "Unruh detections" reinterpret QED, and analogues test only the kinematics

There is still no undisputed laboratory detection of the Unruh effect: "It is currently not clear whether the Unruh effect has actually been observed, since the claimed observations are disputed" ([Wikipedia: Unruh effect](https://en.wikipedia.org/wiki/Unruh_effect)). Sokolov–Ternov spin polarisation in storage rings reaches 8/(5√3) ≈ **92.38%** and can be described with an effective temperature, but circular motion is not a true thermal bath, and the "thermal" and exact polarisation formulas "show substantial differences" ([Crispino et al.](https://arxiv.org/abs/0710.5373)) **[mainstream that the analogy is imperfect]**. Akhmedov and Singleton's argument that the two effects "physically coincide" is **[contested]** ([arXiv:0705.2525](https://arxiv.org/abs/0705.2525)). The review's Bell–Leinaas linac example quotes 2×10¹⁶ g and 0.7 mK for a 10 MV/m field, but those two numbers are mutually consistent only for 1 MV/m; 10 MV/m gives 7.1 mK *(check)*. Chen and Tajima's laser proposal estimates Unruh radiation at "≈ 3×10⁻⁴" of the Larmor background, a figure that reproduces with log₁₀ (3.5×10⁻⁴) but not with ln (8.1×10⁻⁴), so the original must be checked before it is coded ([Chen & Tajima, PRL 83, 256](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.83.256)) *(check)* **[proposal]**. Lynch et al. claim "the first observation of acceleration-induced thermality in a non-analogue system" in positron-channelling data ([PRD 104, 025015](https://arxiv.org/abs/1903.00043)) **[contested]**. A driven BEC of about 10⁵ atoms reproduces Unruh's temperature and entropy as a quantum simulation ([Hu et al., Nat. Phys. 2019](https://arxiv.org/abs/1807.07504)) **[mainstream as a simulation]**, and the 2025 superradiant-mirror ([phys.org on Deswal et al.](https://phys.org/news/2025-11-faint-quantum-space.html)) and annular-Josephson-junction ([Hiroshima University](https://www.hiroshima-u.ac.jp/node/93889)) schemes remain **[proposals]**. Because T_U = 1 K needs 2.5×10¹⁹ g, only charged particles in strong fields get close, and their ordinary radiation is already fully described by inertial-frame QED. Any "detection" reinterprets existing physics rather than finding new physics.

Analogue horizons confirm the kinematics: horizons plus fluctuations give thermal pairs, whatever the medium. Steinhauer's group repeated a Rb-87 condensate experiment 7400 times and found the spectrum "thermal at the Hawking temperature implied by the analogue surface gravity", 0.35 nK, with no analogue firewall ([de Nova et al., Nature 569, 688](https://arxiv.org/abs/1809.00913)) **[mainstream]**; a real solar-mass hole is 176 times hotter *(check)*. The 2016 entanglement claim ([Steinhauer, Nat. Phys. 12, 959](https://arxiv.org/abs/1510.00621)) was questioned by Leonhardt ([arXiv:1609.03803](https://arxiv.org/abs/1609.03803)) and defended by Steinhauer ([arXiv:1609.09017](https://arxiv.org/abs/1609.09017)) **[contested]**, while stationarity was shown in 2021 by observing the system at six different times ([Kolobov et al., Nat. Phys. 17, 362](https://arxiv.org/abs/1910.09363)) **[mainstream]**. Water-wave and optical-fibre analogues showed stimulated, classical mode conversion ([Weinfurtner et al., PRL 106, 021302](https://arxiv.org/abs/1008.1911); [Drori et al., PRL 122, 010404](https://arxiv.org/abs/1808.09244)) **[mainstream]**. None of this says anything about engineering gravity: the analogue "metric" belongs to the fluid, and Einstein's equations play no part.

### Vacuum-inertia and strong-field gravity claims fail on QCD, cosmology and GR

Haisch, Rueda and Puthoff (1994) proposed that "inertia is an electromagnetic resistance" to the zero-point field (ZPF) as seen in accelerated frames ([PRA 49, 678](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.49.678)). Their mass formula, m = Γ_zħω_c²/(2πc²), has two free parameters by the authors' own account ([arXiv:physics/9807023](https://arxiv.org/abs/physics/9807023)). With a Planck cutoff the electron needs Γ_z = 1.4×10⁻⁶⁵ s, about 10⁻²² Planck times, so the parameter simply encodes the mass *(check)*. Woodward pointed out that hadron inertia sits in electrically neutral gluon fields ([Found. Phys. 31, 819](https://dx.doi.org/10.1023/A%3A1017500513005)), so the electromagnetic version is **[refuted]**; the generalised reply that adds gluon and W/Z fields is **[contested]** and makes no quantitative prediction ([Dobyns, Rueda & Haisch](https://arxiv.org/abs/gr-qc/0002069)). Two measurements leave such models little room. MICROSCOPE's η = [−1.5 ± 2.3 ± 1.5]×10⁻¹⁵ ([PRL 129, 121102](https://arxiv.org/abs/2209.15487)) forces any material-dependent mass mechanism to be tuned into universality, and Rainville et al. confirmed E = mc² in nuclear reactions to 4×10⁻⁷ ([NIST](https://www.nist.gov/node/399321)) **[mainstream]**. The 1998 paper's suggestion that "the energy released in fusion would be coming from the ZPF" is **[contested]** and makes no distinct prediction. Puthoff's 1989 ZPF account of gravity did not survive Cole, Rueda and Danley, who found it does not yield Newtonian gravity "at least not with this particular starting point" ([PRA 63, 054101](https://www.bu.edu/simulation/publications/lab/pdf/PRA2001.pdf); [Puthoff 1989](https://link.aps.org/doi/10.1103/PhysRevA.39.2333)) **[refuted]**.

Treating the vacuum as an energy reservoir fails on cosmology. Cut off at the Planck frequency, the ZPF holds 5.9×10¹¹¹ J/m³, against an observed dark-energy density of 5.3×10⁻¹⁰ J/m³: a gap of **10¹²¹**, with the matching cutoff at 6.7 meV *(check)*. That is the cosmological-constant problem **[mainstream]** ([Weinberg, RMP 61, 1](https://doi.org/10.1103/RevModPhys.61.1)), and it rules out any scheme that treats Planck-cutoff vacuum energy as a gravitating, usable reservoir. Strong-field "gravity control" fails on GR and QED together. A field at E_S stores 7.75×10²⁴ J/m³, equivalent to 8.6×10⁷ kg/m³, yet a 1 m sphere of it produces a potential of only GM/(rc²) = 2.7×10⁻¹⁹, and pair creation would drain the field in about 10⁻¹⁷ s *(check)*. Pais-style gravity control is **[refuted]** on both counts.

### The Sun burns at 10⁻⁶ Gamow suppression because gravity can wait

Gravity does not lower the barrier; it confines and compresses. NASA gives the solar centre as 1.571×10⁷ K, 1.622×10⁵ kg/m³ and 2.477×10¹⁶ Pa ([NASA Sun fact sheet](https://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html)), and models put the central fusion power at "about 276.5 watts/m³", "similar to an active compost heap" ([Wikipedia: Solar core](https://en.wikipedia.org/wiki/Solar_core)) **[mainstream]**. The p+p barrier at 2 fm is 720 keV, 532 times kT, and the Gamow peak at 6.09 keV is suppressed to 1.4×10⁻⁶ *(check)*. With S₁₁(0) = (4.01 ± 0.04)×10⁻²² keV b and its derivative ([Solar Fusion II](https://arxiv.org/abs/1004.2318)), the unscreened rate is ⟨σv⟩ = 1.6×10⁻⁴⁹ m³/s, which at X_c = 0.34 gives 183 W/m³, 0.66 of the model value *(check)*. Screening, the choice of T_c (the rate goes roughly as T⁴), the ppII and ppIII branches and CNO plausibly close the gap, but it was not reconciled. The implied proton lifetime of 6.0–6.5 Gyr also conflicts with Wikipedia's "about one billion years", perhaps through a per-proton against per-pair definition; that is unresolved.

Simple gravity gets the conditions right to within a factor of 2 *(check)*. An n = 3 polytrope gives ρ_c/ρ̄ = 54.18 and P_c = 1.245×10¹⁶ Pa, half the standard solar model's value; the virial temperature GMμm_u/(k_BR) is 1.40×10⁷ K; and the Kelvin–Helmholtz time of 15.7 Myr measures how long gravity holds the heat. The Sun can burn at all because of a trade. The p+p S-factor is 1.4×10²³ times smaller than D(d,p)'s 55.6 keV b (the Bosch–Hale A₁), because the first step is a weak interaction. The Sun pays for that with a "Lawson product" n_ekT_cτ_KH of 4.4×10⁴⁶ keV s m⁻³, 1.5×10²⁵ times the D–T ignition requirement *(check)*. Stellar gravity's Unruh temperature of 10⁻¹⁸ K contributes nothing, and there is no "quantum-gravity assist" to fusion anywhere short of a horizon.

## Conclusion

The three barriers in this report are one exponent in three guises: exp(−√(E_G/E)) for nuclei, exp(−πE_S/E) for the vacuum, and exp(−mc²/2k_BT) at a horizon. Every route that works raises the energy entering the exponent instead of lowering the barrier. Plasmas heat the fuel until the Gamow peak sits three or four kT out, NA63 and E144 used Lorentz factors of 10⁵ to turn laboratory fields of about 10⁻⁶ E_S into rest-frame fields near E_S, and stars wait for billions of years. The genuine barrier-lowering effects are small and bounded: a few hundred eV of screening, and a 40% cut in the exponent from assisted tunnelling that so far exists only in theory. The fringe claims instead promise to change a prefactor by 10⁶ to 10⁵³, or never compute the exponent at all, as Pais's eq. 4 does when it equates two energy densities. That gives the crate its design rule: every claim check reports the controlling exponent beside the claimed number, and every test asserts the ratio to the bound.

The contested frontier is narrow, quantitative and worth tracking. The screening puzzle is the one place where mainstream theory falls short of the data, by about 2× in gas and about 3× in metals, and the 2025 Nature result shows that electrochemistry can move beam-driven rates by 15%. A verdict here could change with LUXE phase 1 crossing χ = 1 with a laser, a PVLAS-class measurement reaching Δn = 2.5×10⁻²³, or impact-parameter-resolved U–U collisions at CRYRING. None of these would revive the heat claims, since even complete success for MARRS would leave Fleischmann–Pons 2×10⁵ short. A screening calculator that exposes both forms, beside Bosch–Hale with enforced validity ranges, lets this repository's code separate the real anomaly, a factor of 1.2–3 at keV energies, from its extrapolation, a factor of 10⁴⁰ at room temperature.

## What goes into the crates

The calculators belong in the repository's crate for exotic and frontier physics, currently [`crates/gravity-claims`](../../../crates/gravity-claims) and about to be renamed. Two modules are new. `fusion` covers the barrier, reactivities, Lawson, screening, muon catalysis, stellar confinement and the fusion claim checks. `strong_field` covers Schwinger rates, laser parameters, vacuum birefringence, supercritical charge and analogues. The horizon temperatures go into the existing `spacetime` module beside its Casimir and quantum-inequality functions, and the Pais checks reuse the existing `hfgw::surface_field` and `hfgw::VACUUM_BREAKDOWN`. As before, each function is written from the cited formula rather than translated from anyone's code, so the MIT licence holds, and every claim test asserts the published number and its ratio to the mainstream function beside it. A sensible build order is the `fusion` core (`gamow_energy`, `gamow_peak`, `bosch_hale`, `lawson`, `screening`), then the `strong_field` core (`schwinger_rate`, `a0`, `chi_e`, `birefringence`), then the horizons, then the claim checks.

Five conventions decide whether the tests pass:

1. **Constants.** The crate uses CODATA 2018 (m_e = 9.1093837015×10⁻³¹ kg), while `schwinger_checks.py` uses CODATA 2022. The resulting difference of about 3×10⁻⁸ in E_S sits well inside a 10⁻⁶ formula tolerance. New constants are needed for α, k_B, the nuclear masses in keV (p 938272.08816, d 1875612.9424, t 2808921.1325, ³He 2808391.6074) and τ_μ = 2.1969811 µs.
2. **Bosch–Hale constants.** The Bosch–Hale tests must use the tabulated B_G and m_rc², not CODATA-derived values, which differ by 0.017% for D–³He.
3. **Intensity and a₀.** These take an explicit enum: linear (ε₀cE²/2) or Gaussian "E²/4π" (ε₀cE²). Most 1.3–1.4× mismatches in the literature are this convention.
4. **Validity ranges.** Functions return `None` outside them: Bosch–Hale 0.2–100 keV (³He(d,p) 0.5–190 keV), Tentori–Belloni 70–500 keV, Δn only for B ≪ B_S, and the Schwinger rate only for a constant field.
5. **Tolerances.** Test formulas to 10⁻⁶ against the scripts, and published comparisons to their rounding: 0.3% for UWFDM Table 6, two significant figures for Wurzel–Hsu, 1.6% for LUXE χ_e, and 3% for Di Piazza's χ₀ coefficient, whose 2.7% gap is unexplained.

### `fusion`: mainstream calculators

| Function | Formula | Test vector (from the notes) | Label |
|---|---|---|---|
| `coulomb_barrier(z1, z2, r_fm)` → keV | 1.439965·Z₁Z₂/R (MeV fm) | (1, 1, 4.0 fm) → 360.0 keV (textbook 360 keV); D–T at R = 3.243 fm → 444 keV; p–¹¹B at R = 4.514 fm → 1595 keV | mainstream; the caller chooses the R convention |
| `sommerfeld(z1, z2, mr_kev, e_kev)` and `tunnelling(…)` = e^(−2πη) | η = Z₁Z₂α√(m_rc²/2E) | D–T, 10 keV, m_rc² = 1124647.3 keV → η = 1.7304, P = 1.8969×10⁻⁵; D–D, 1 keV → 4.9970, 2.3147×10⁻¹⁴; p–¹¹B, 10 keV → 7.5643, 2.2848×10⁻²¹ | mainstream |
| `gamow_energy(z1, z2, mr_kev)` | 2m_rc²(παZ₁Z₂)² | D–T 1182.2, D–D 985.8, D–³He 4728.3, p–¹¹B 22589.1, p–p 493.1 keV; √(D–T) = 34.3826 against Bosch–Hale's 34.3827 | mainstream |
| `gamow_peak(e_g, kt)` → (E₀, Δ) | E₀ = (E_G kT²/4)^⅓; Δ = 4√(E₀kT/3) | D–T at kT = 10 keV → 30.918, 40.607 keV; p+p at 1.3357 keV (1.55×10⁷ K) → 6.036, 6.557 keV; p–¹¹B at 300 keV → 798.0, 1130.0 keV; Rauscher's 0.1220/0.2368 forms within 0.3% | mainstream |
| `bosch_hale(reaction, t_kev)` → cm³/s | θ, ξ, C₁…C₇ with Bosch–Hale's B_G and m_rc² | all 28 values of UWFDM-1268 Table 6 within 0.3%; D–T at 10 keV → 1.1362×10⁻¹⁶; D(d,n) at 0.1 and 100 keV → 1.6572186×10⁻³⁷ and 2.6816742×10⁻²³ m³/s (exact MFE endpoints); D–T maximum 8.9465×10⁻¹⁶ at 66.7 keV | mainstream; 0.2–100 keV |
| `bosch_hale_sigma(reaction, e_kev)` → mb | S(E) = (A₁+E(A₂+…))/(1+E(B₁+…)); σ = S/(E·e^(B_G/√E)) | D–T at 10 keV → S = 1.4246×10⁷ keV mb, σ = 27.02 mb; D–T peak 5065.7 mb at 64.8 keV; D(d,n) and D(d,p) at 50 keV → 16.487 and 15.573 mb | mainstream; D–T low-energy range (≤ 550 keV) only |
| `pb11_reactivity(fit, t_kev)` → m³/s | Bosch–Hale form; E_G = 22589.1 keV, m_rc² = 859604.5 keV, P₁…P₇ for Nevins–Swain 2000 and Tentori–Belloni 2023 | Tentori–Belloni at 100, 300, 500 keV → 7.0650×10⁻²³, 4.3870×10⁻²², 5.6190×10⁻²²; TB/NS − 1 = +29.6% at 300 and +50.2% at 500 keV (published 30% and 50%) | mainstream; 70–500 keV |
| `lawson(t_kev, q_fuel)` → nτ_E | 3T/[(f_c + 1/Q)⟨σv⟩ε_F/4 − C_B√T], f_c = 3.5/17.6, C_B = 5.34×10⁻³⁷ W m³ keV^−½ | Q = ∞: minimum nTτ 2.968×10²¹ keV s m⁻³ at 14.25 keV, minimum nτ 1.567×10²⁰ s m⁻³ at 26.0 keV; Q = 1: 4.677×10²⁰ at 13.65 keV, 2.527×10¹⁹ at 25.7 keV; ignition at 4.305 keV (published 2.9×10²¹, 1.6×10²⁰, 4.6×10²⁰, 2.5×10¹⁹, 4.3 keV) | mainstream |
| `screening(z1, z2, mr_kev, e_kev, ue_kev, form)` | first order e^(πηU_e/E), or exact shift P(E+U_e)/P(E) | d+d, 25 eV at 10 keV → 1.0125 (both forms); d+d, 300 eV at 1 keV → 111.0 / 47.47; d+³He, 219 eV at 1 keV → 1862 / 653.5, at 5 keV → 1.961 / 1.920 | mainstream formula; metal U_e inputs contested |
| `salpeter(z1z2, zeta, rho_gcc, t6)` | e^(0.188·Z₁Z₂·ζ·√ρ·T₆^−³ᐟ²) | (1, 1.291, 150, 15.7) → 1.0489; T₆ = 15.5 → 1.0499 | mainstream (composition assumed) |
| `nonresonant_reactivity(z1, z2, mr_kev, s0, s1, kt, f0)` | Solar Fusion II Eq. 8 with S_eff = S(0)[1 + 5kT/36E₀] + S′(0)E₀[1 + 35kT/36E₀] | p+p at 1.55×10⁷ K, S(0) = 4.01×10⁻²² keV b, S′ = 4.49×10⁻²⁴ b, f₀ = 1.049 → S_eff = 4.4626×10⁻²² keV b, ⟨σv⟩ = 1.5929×10⁻⁴³ cm³/s (1.5973×10⁻⁴³ from the 1.301×10⁻¹⁴ form) | mainstream |
| `fusion_power_density(n1, n2, sigma_v, q, identical)` | n₁n₂⟨σv⟩Q/(1 + δ₁₂) | n = 10²⁰ m⁻³, D–D total 5.002×10⁻²⁴ m³/s at 20 keV, Q = 3.651 MeV → 1.463×10⁴ W/m³, or 207 W in a sphere of radius 0.15 m | mainstream |
| `muon_yield(sticking, cycle_rate)` and `muon_breakeven(e_mu_mev)` | 1/(ω_s + λ₀/λ_c), λ₀ = 1/τ_μ = 4.5517×10⁵ s⁻¹; E_μ/17.6 MeV | (0.0047, 2.5287×10⁸ s⁻¹) → 153.85 (Kelly ≈ 154; measured 150 ± 20); 1/ω_s = 110 at 0.912%; breakeven 267, 284, 341 at 4.70, 5, 6 GeV (Pahlavani "about 285" at 5 GeV) | mainstream |
| `DD_MEAN_Q` and `dd_per_watt()` | Q̄ = (3.268 + 4.033)/2 MeV = 3.651 MeV | 1.710×10¹² reactions s⁻¹ W⁻¹ and 8.548×10¹¹ n s⁻¹ W⁻¹ (MARRS: "1 Watt~10¹² reactions/second") | mainstream |
| `lane_emden(n)` → (ξ₁, −ξ₁²θ′, ρ_c/ρ̄) | RK4 integration of θ″ + (2/ξ)θ′ + θⁿ = 0 | n = 3 → 6.897, 2.018, 54.18; with GM²/[4π(n+1)θ′(ξ₁)²R⁴] the Sun's P_c = 1.245×10¹⁶ Pa (NASA 2.477×10¹⁶) | mainstream |
| `kelvin_helmholtz(m, r, l)` and `virial_temperature(m, r, mu)` | GM²/(2RL); GMμm_u/(k_BR) | Sun → 15.7 Myr; μ = 0.61 → 1.40×10⁷ K (NASA T_c 1.571×10⁷ K) | mainstream |

The `fusion` module also needs integration tests that combine functions rather than add API:

- The Maxwellian average of `bosch_hale_sigma` matches `bosch_hale` within 1.2% for D–T and D–D at 2, 10 and 50 keV, and within 3.1% for D–³He.
- The unscreened numerical p+p rate at 1.571×10⁷ K is 1.6×10⁻⁴⁹ m³/s, giving 183 W/m³ at ρ_c = 1.622×10⁵ kg/m³ and X = 0.34 (the model gives 276.5). The solar "Lawson product" is 4.39×10⁴⁶ keV s m⁻³.
- Device energy ratios: pyroelectric 1.46×10⁻⁶, Eindhoven fusor 9.07×10⁻¹⁰, NASA LCF 4.03×10⁻¹⁴.
- NIF target gains: 1.537, 1.893, 2.439 and 4.135.
- Koonin–Nauenberg's 3×10⁻⁶⁴ s⁻¹ gives one fusion per 1.75×10³² yr per mole.

### `strong_field`: mainstream calculators

| Function | Formula | Test vector (from the notes) | Label |
|---|---|---|---|
| `E_SCHWINGER`, `B_SCHWINGER`, `critical_intensity(conv)` | m_e²c³/(eħ); E_S/c; ε₀cE_S²·{1, ½} | 1.32329×10¹⁸ V/m (Di Piazza 1.3×10¹⁶ V/cm); 4.41401×10⁹ T; 4.64811×10²⁹ and 2.32406×10²⁹ W/cm² (Di Piazza 4.6×10²⁹) | mainstream |
| `schwinger_rate(e)` and `vacuum_decay_rate(e)` | (αε₀E²/π²ħ)·e^(−πE_S/E); same prefactor × Σ n⁻²e^(−nπE_S/E) | prefactor 1.08704×10⁵⁶ m⁻³ s⁻¹; decay/rate at E_S = 1.011016; E/E_S = 0.1, 0.05, 0.01 → 2.4688×10⁴⁰, 1.4017×10²⁶, 3.9683×10⁻⁸⁵; π × prefactor = 3.415×10⁵⁰ cm⁻³ s⁻¹ (Melissinos 3.45×10⁵⁰, ratio 0.990); scalar/spinor = 0.5 | mainstream; constant field only |
| `schwinger_exponent(e)` | πE_S/E | 4566 at 9.104×10¹⁴ V/m (the record laser); identity test with `spacetime::unruh_temperature` below | mainstream |
| `a0(intensity, lambda, conv)` | eE₀/(m_ecω) | 1.1×10²³ W/cm² at 0.8 µm → 226.84 (linear) and 160.40 (Gaussian); coefficients 0.85493 and 6.045 (published 0.855 and 6.0); LUXE 1.33×10²⁰ → 7.888 (7.9), 1.2×10²¹ → 23.69 (23.6) | mainstream |
| `chi_e(energy, a0, photon_ev, angle)` | (1 + β cos θ)·γ·a₀·ħω/(m_ec²) | LUXE (16.5 GeV, a₀ = 7.9, 1.55 eV, 17.2°) → 1.5127 (published 1.5); a₀ = 23.6 → 4.519 (4.45); E144 η (46.6 GeV, 2.3526 eV, 17°) → 0.8214 (ESPPU 0.83); 10 GeV head-on at 800 nm → 0.1187 per unit a₀ | mainstream |
| `ritus_narozhny(chi)` | αχ^⅔ | equals 1 at χ = 1604.18 (published ≈ 1600); 0.0197 at χ = 4.45 | contested regime |
| `assisted_action(gamma)` → Ã, with action = (E_S/E)·Ã | Ã = 2 arcsin(π/2γ) + (π/2γ²)√(4γ² − π²) for γ ≥ π/2, else π | Ã(π/2)/π = 1; 1 − Ã(π)/π = 0.391 (published "about 40%"); 100·Ã(100)/2π = 0.99996; at E = 0.1 E_S, e^(−10·Ã(π)) = 4.9088×10⁻⁹ against e^(−10π) = 2.2711×10⁻¹⁴ | mainstream theory, unobserved |
| `birefringence(b)` | Δn = 3A_eB², A_e = 2α²ħ³/(45μ₀m_e⁴c⁵) | A_e = 1.32467×10⁻²⁴ T⁻²; 2.5 T → 2.48375×10⁻²³ (PVLAS's QED value 2.5×10⁻²³); n∥ − 1 = 2.898×10⁻²³, n⊥ − 1 = 1.656×10⁻²³; returns `None` for B > 0.1 B_S | mainstream |
| `eh_fraction(e)` | (α/45π)(E/E_S)² | 0.1 E_S → 5.16×10⁻⁷; 9.10×10¹⁴ V/m → 2.44×10⁻¹¹ | mainstream |
| `breit_wheeler_suppression(chi)` | e^(−8/(3χ)) | χ = 0.1 → 2.6231×10⁻¹²; 0.3 → 1.3791×10⁻⁴; 1 → 6.9483×10⁻² | mainstream (quasi-static asymptote) |
| `dirac_1s_binding(z)` and `Z_CRITICAL_EXTENDED = 173` | m_ec²[1 − √(1 − (Zα)²)] | Z = 82 → 101.581 keV; Z = 92 → 132.28 keV; `None` above Z = 137.036 | mainstream (point nucleus; 173 is published) |
| `graphene_critical_field(gap_ev, v_f)` and `graphene_pair_rate(e, v_f)` | Δ²/(eħv_F); π⁻²ħ^−³ᐟ²v_F^−½(eE)^³ᐟ² (four species) | 0.2 eV, 10⁶ m/s → 6.077×10⁷ V/m (Schmitt ≲ 6×10⁷); 10⁴ and 10⁶ V/m at 10⁶ m/s → 6.0×10²⁴ and 6.0×10²⁷ m⁻² s⁻¹ | mainstream analogue |

These need tests only, with no public function:

- CoReLS gives E₀/E_S = 6.880×10⁻⁴ and I_S/I = 2.113×10⁶.
- One pair in (1 µm)³ within 1 µm/c needs 0.06443 E_S.
- NA63's lab-frame fields χE_S/γ are 5.312×10¹² V/m (Si) and 8.114×10¹² V/m (Ge).
- E144's photon threshold is 3.913.
- τ₀ = 2r_e/3c = 6.266×10⁻²⁴ s (ESPPU 6.4×10⁻²⁴).
- The two Gaunt approximations give g(1) = 0.1795 and 0.1811.
- Graphene: n = 3.121×10¹⁰ cm⁻² (published 3×10¹⁰), and G₀ = 1.20 × 4e²/h = 0.186 mS.
- U+U at 16.5 fm has a Coulomb energy of 738.7 MeV and needs 6.207 MeV/u head-on.
- Significances: PVLAS σ/QED = 6.84; STAR 6.72σ.

### `spacetime`: horizon additions

| Function | Formula | Test vector (from the notes) | Label |
|---|---|---|---|
| `horizon_temperature(kappa)`, with `unruh_temperature(a)` (κ = a/c), `hawking_temperature(m)` (κ = c³/4GM) and `gibbons_hawking_temperature(h)` (κ = H) | ħκ/(2πk_B) | 4.055×10⁻²¹ K per m s⁻² (Wikipedia 4.06×10⁻²¹); 1 K at 2.466×10²⁰ m/s² (2.47×10²⁰); g₀ → 3.977×10⁻²⁰ K; 1 M☉ → 6.17×10⁻⁸ K; 5×10¹¹ kg → 2.454×10¹¹ K; T_CMB at 4.502×10²² kg; H₀ = 67.4 → 2.655×10⁻³⁰ K; κ = 287.9 s⁻¹ → 0.35 nK (Steinhauer); Rogers' trap → 2.405 K (≈ 2 K) | mainstream |
| identity test (no new function) | `strong_field::schwinger_exponent(E)` = m_ec²/(2k_B·`unruh_temperature`(eE/m_e)) | holds to 10⁻¹² at three fields; T_U(E_S) = 9.438×10⁸ K = 81.33 keV | mainstream correspondence |
| `hawking_lifetime(m, model)` | photons only: 5120πG²M³/(ħc⁴); Page: G²M³/(3·2×10⁻⁴ħc⁴) | 1 M☉ → 2.095×10⁶⁷ yr (Wikipedia 2.140×10⁶⁷; 2% gap unexplained) and 2.171×10⁶⁶ yr; Page coefficient 8.716×10⁻²⁷ s g⁻³ | mainstream; the Page coefficient comes from the rounded abstract |
| `extremal_rn_field(m)` | c⁴/(√(4πε₀)·G^³ᐟ²·M) | equals E_S at 1.061×10³⁶ kg (5.34×10⁵ M☉); 10 M☉ → 7.063×10²² V/m | mainstream (our derivation) |
| `sokolov_ternov_time(gamma, rho)` and `SOKOLOV_TERNOV_P0` | 8m_e²cρ³/(5√3·αħ²γ⁵); P₀ = 8/(5√3) | P₀ = 0.92376; 3.7 GeV, ρ = 12.7 m → 4.857 min; 45.6 GeV, ρ = 3096 m → 4.125 h (bend-only, with assumed radii) | mainstream |
| `zpf_energy_density(omega_c)` | ħω_c⁴/(8π²c³) | Planck cutoff (1.855×10⁴³ s⁻¹) → 5.868×10¹¹¹ J/m³; Compton cutoff → 1.801×10²² J/m³, which is ε₀E_S²/2 divided by π/α = 430.5; matches ρ_Λ = 5.253×10⁻¹⁰ J/m³ at ħω = 6.678 meV | mainstream identity |
| `field_compactness(e, r)` | G·(ε₀E²/2c²)·(4πr³/3)/(rc²) | E_S at 1 m → 2.683×10⁻¹⁹; at 1 km → 2.683×10⁻¹³ | mainstream |

### Labelled claim checks

Each claim check sits in its topical module, and its doc comment opens with the label and names the bound it violates.

| Function or test | Formula | Test vector (from the notes) | Label and bound |
|---|---|---|---|
| `charged_sphere(q, r)` → (field, potential, energy) (`fusion`) | Q/(4πε₀r²), Q/(4πε₀r), Q²/(8πε₀r) | 1 C, 0.15 m → 3.99×10¹¹ V/m, 5.99×10¹⁰ V, 3.00×10¹⁰ J; 1 m → 8.99×10⁹ V/m | refuted; ÷ `hfgw::VACUUM_BREAKDOWN` = 399 |
| `spinning_shell_field(q, r, omega)` (`fusion`) | (2/3)μ₀Qω/(4πR) | 1 C, 0.15 m, 10⁴ rad/s → 4.44×10⁻³ T | refuted; the claimed 10⁶ T is 2.2×10⁸ times this |
| Pais power test | `fusion_power_density` with `bosch_hale` (D–D total at 20 keV) | 207 W; 1 GW ÷ 207 W = 4.8×10⁶ (1 TW: 4.8×10⁹); 1 GW needs n = 2.2×10²³ m⁻³; at 10⁶ T, B²/2μ₀ = 3.98×10¹⁷ Pa and 5.62×10¹⁵ J in a 0.3 m sphere | refuted |
| Pais vacuum test | `schwinger_rate`, `eh_fraction`, `zpf_energy_density`, `hfgw::surface_field` | `surface_field(5e4)` = 5.647×10¹⁵ V/m; log₁₀ Γ = −268.4; EH fraction 9.40×10⁻¹⁰; craft eq. 4 → ν = 4.21×10²¹ s⁻¹ = 5.42 ω_C | refuted |
| `dd_neutron_deficit(watts, n_claimed)` (`fusion`) | P · 8.548×10¹¹ against the claim | 1.751 W, 4×10⁴ → 1.497×10¹² n/s, deficit 3.74×10⁷; 26.8 W → 2.291×10¹³ n/s; rate per pair 3.50×10⁻¹¹ s⁻¹ = 1.2×10⁵³ × Koonin–Nauenberg | refuted (F–P heat as D–D) |
| `widom_larsen_rate(beta, beta0)` (`fusion`) | 1.86×10⁻³(β − β₀)² Hz | β = 20, β₀ = 2.59 → 0.564 Hz (published ≈ 0.56); ÷ Ciuchi's 1.8×10⁻³ Hz = 313 (published "about 300"); β₀ from (m_n − m_p)/m_e = 2.531 | refuted |
| `hrp_gamma_z(m, omega_c)` (`spacetime`) | 2πmc²/(ħω_c²) | electron at the Planck cutoff → 1.418×10⁻⁶⁵ s; proton → 2.603×10⁻⁶² s, against t_P = 5.391×10⁻⁴⁴ s | refuted (EM-only inertia model) |

These checks are a multiplication or two each, so they stay as documented assertions:

| Check | Values to assert | Label |
|---|---|---|
| MARRS in watts | 0.1, 10³, 10⁶ and 10⁹ s⁻¹ g⁻¹ → 5.8×10⁻¹⁴, 5.8×10⁻¹⁰, 5.8×10⁻⁷ and 5.8×10⁻⁴ W/g; F–P's 1.98×10¹¹ s⁻¹ g⁻¹ = 2.0×10¹² × the state of the art | contested premise |
| Bubble fusion | 4×10⁴, 8×10⁴ and 7×10⁵ n/s → 4.68×10⁻⁸, 9.36×10⁻⁸ and 8.19×10⁻⁷ W; tritium/neutron ratio 8.8–17.5 | refuted |
| Lugano | 5825 MJ per gram (1.3×10⁵ × petrol); full burn-up of 0.011 g of ⁷Li → 2625 MJ; replica File 4: 1658 W → 102.9 W (16×) | refuted (COP); contested (residual) |
| Hydrino | 13.6 eV/q² → 54.4 eV (q = ½), 217.6 eV (q = ¼) | refuted |
| Debye screening | R_D = 4.88×10⁻¹² m and U_D = 295 eV at 300 K and 6×10²⁸ m⁻³ | refuted as explanation |
| "10¹⁸" enhancement | needs U_e of 51, 195 and 992 eV at 50, 100 and 200 eV centre-of-mass | contested |
| F–P energy per atom | 4 MJ cm⁻³ = 367 eV per Pd atom | refuted as D–D |

### What stays out

Leave these out as physics: any model of cold-fusion heat; the Widom–Larsen mechanism, keeping only the rate reproduction above; hydrino energy levels as states; Pais's B_MAX and "gain ~ B³" relations and his ε₀E² = ħν⁴ condition, except as the claim tests above; Haisch–Rueda–Puthoff inertia as a mass model; every zero-point-energy device; and Holmlid's H(0). The naive thermal screening model of Pd–D is the notes' own construction, not a published model, so it belongs in a test comment at most, even though it reproduces MARRS's 1 W/g at 10⁻¹³ b (2.51 W/g).

Five mainstream items should also wait:

- the general-field Nikishov rate, until its factors are checked against Fedotov et al. (2023);
- strong-field Euler–Heisenberg above B_S, which magnetar work needs and the notes do not cover;
- the Tentori–Belloni 10–70 keV fit and the Bosch–Hale D–T high-energy S-fit, neither of which was transcribed;
- Melissinos's interpolating g(η) and the Chen–Tajima Unruh/Larmor ratio, until the ln or log₁₀ question is settled against the originals;
- any theory of metal screening (Debye, ion tracks, threshold resonance), since the review literature says none fits yet.

The calculator to build first is `screening` with both forms, beside `bosch_hale` with enforced validity ranges. That pair alone separates the measured anomaly from the room-temperature extrapolation that every cold-fusion argument depends on.
