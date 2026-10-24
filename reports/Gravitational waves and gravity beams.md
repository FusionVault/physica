# Chart gravitational waves honestly, from nanohertz to terahertz

The most useful thing these crates can gain from 2023–2026 gravitational-wave research is a permissively licensed **strain-and-sensitivity crate** (working name `gwsens`, unclaimed on crates.io on 24 October 2026). It would put LIGO, LISA, pulsar timing arrays, atom interferometers, GRACE-FO and the megahertz-to-terahertz detector concepts on one axis beside a hard cosmological ceiling. Alongside it goes a thin gravitational-wave layer in `satgrav` that reuses the atom-interferometer transfer function the crate already has. Gravity-beam generators are not worth building. High-frequency gravitational waves (HFGWs) are real physics with real new experiments: QUEST's first broadband 13–80 MHz limits, an ABRACADABRA-10cm search, a CAPP-12T reanalysis and a revived MAGO cavity. But the 2025 Living Review states that no concept, built or proposed, reaches any realistic source. Every early-Universe background is capped at **Ω_GW h² < 1.1×10⁻⁶**, about 10⁻³⁰ strain at 1 GHz, while the best laboratory limits sit **18–31 orders of magnitude** above that line. "Gravity beams" split three ways:

- **Established GR, but unmeasurable:** the gravitational near field of a light or particle beam, about 10⁻¹⁷ m/s² for a petawatt beam at 1 mm.
- **Serious but speculative:** collimated wave emission through electromagnetic conversion, with an optimistic strain of 10⁻³².
- **Fringe or refuted:** Podkletnov's gravity shielding, Li–Baker HFGW propulsion, Ning Li's "AC gravity" and the Tajmar gravitomagnetic anomalies.

Gravitomagnetism itself is real: LARES-2 measured Earth's frame dragging to about one part in a thousand in 2026, which makes it a natural, validated addition to `geo5d`. In Rust, waveforms were covered in October 2026 by a GPL-licensed crate, so the open gap is sensitivity curves, signal-to-noise (SNR) and background tooling, and HFGW curves; no tool for the last exists in any language. The project's real research edge is where gravimetry and gravitational-wave instruments share physics. Newtonian noise is gravimetric forward modelling, any gradiometer converts to strain by one line, and the "gravity waves" already present in GOCE and GRACE-FO data are atmospheric buoyancy waves. The 2022 Tonga eruption is an untested target in GRACE-FO ranging residuals.

## A cosmological ceiling caps every high-frequency background near 10⁻³⁰ strain

The reference for everything above 10 kHz is the 2025 update of the community Living Review, "Challenges and Opportunities of Gravitational Wave Searches above 10 kHz" (Aggarwal et al., Living Rev. Relativ. 28:10, v2 February 2026) ([Review 2025](https://arxiv.org/abs/2501.11723)). Its most important number is not a detector sensitivity but a bound. Gravitational waves present before Big Bang nucleosynthesis (BBN) count as extra radiation, so the combined CMB and BBN limit ΔN_eff < 0.18 caps their energy density today at **Ω_GW,0 h² < 1.1×10⁻⁶**. CMB-only limits tighten this to 2.9×10⁻⁷, and a future CMB measurement of N_eff to its theoretical precision would push it to about 5.6×10⁻⁹ ([Review 2025 §4.2](https://arxiv.org/abs/2501.11723); [Yeh et al. 2022](https://arxiv.org/abs/2207.13133)).

In strain units the bound is brutal. With the one-sided convention of Moore, Cole and Berry, **h_c = 1.263×10⁻¹⁸ (1 Hz/f) √(h²Ω)** ([Moore, Cole & Berry 2015](https://arxiv.org/abs/1408.0740)). The ceiling is then 1.3×10⁻²⁶ at 100 kHz, 1.3×10⁻²⁷ at 1 MHz, **1.3×10⁻³⁰ at 1 GHz** and 1.3×10⁻³³ at 1 THz (computed here). The review says plainly that "no proposal above the LIGO–VIRGO–KAGRA band currently reaches below the cosmological bound". It reads current sensitivities, which correspond to Ω_GW ≫ 1, only as probes of local overdensities ([Review 2025 §3](https://arxiv.org/abs/2501.11723)). A 2026 analysis goes further. It derives a universal quantum limit across detection schemes and finds that stochastic backgrounds are "in principle observable" only at kHz–MHz, while higher-frequency backgrounds lie below the quantum limit itself ([Guo et al., PRD 113, 044020](https://arxiv.org/abs/2501.18146)).

Within that ceiling the sources sort cleanly into what is physically guaranteed and what is hypothetical. The table gives the formulas a Rust crate can encode. Each is a short closed form that the code should clip against the ceiling.

| Source | Frequency and amplitude (Review 2025 numbering) | Is the source real? |
|---|---|---|
| Primordial black hole (PBH) inspiral | f_ISCO = 4400 Hz (M_⊙/M); 0-PN h̃ ∝ D⁻¹ M_c^(5/6) f^(−7/6); ḟ = (96/5)π^(8/3)(GM_c/c³)^(5/3) f^(11/3) (Eqs. 46–50) | GR is exact. Sub-solar PBHs are hypothetical. |
| Superradiant boson cloud around a light PBH | f ≈ 5 MHz (μ/10⁻⁸ eV); scalar h ≈ 5×10⁻³⁰ (m/10⁻⁵ M_⊙)(α/0.1)⁷(kpc/D) (Eqs. 73–74) | Needs an ultralight boson, which is hypothetical. This is the main target of near-term searches. |
| Cosmic gravitational microwave background (CGMB) | peak **79.8 GHz (106.75/g*s)^(1/3)**; Ω ∝ T_max/M_Pl (Eqs. 97–98) | Guaranteed by the Standard Model, but near the ceiling only if T_max ~ 10¹⁹ GeV. |
| First-order phase transition | f_peak ≃ 260 MHz (H*R*)⁻¹(T*/10¹⁵ GeV); shape S(x) = x³[7/(4+3x²)]^(7/2) (Eqs. 99–101) | Speculative. Ω ≲ 10⁻⁷ needs strong supercooling. |
| Cosmic strings | plateau Ω ≈ 8.04 Ω_rad √(Gμ/Γ), Γ ≈ 50 (Eq. 102) | Speculative. The amplitude is disputed by large factors (the Nambu–Goto versus field-theory debate). |
| Preheating, oscillons | 10⁶–10⁹ Hz, Ω ≲ 10⁻¹⁰ (Table 2) | Model-dependent. |
| Inflation (vacuum) | Ω ≲ 10⁻¹⁶ (Table 2) | Standard, but hopelessly small. |
| Known astrophysics | Neutron-star post-merger ≲ 5 kHz; solar thermal plasma h_c ~ 10⁻⁴² near 10¹⁸ Hz (Table 1) | Real, and negligible. |

All rows are from [Review 2025 §4](https://arxiv.org/abs/2501.11723). The superradiance frequency is f ≈ 4.84×10¹⁴ Hz per eV of boson mass (computed here, neglecting binding corrections), which matches the review's normalisation.

The physics of the waves is settled GR. What is speculative is whether any source exists in the band loudly enough. The only sources that can locally exceed the ceiling are late-Universe and coherent: PBH mergers, boson clouds and gravitational-wave memory. They escape the ΔN_eff bound because they form after BBN, but their strain falls as 1/D and their nearby event rates are small. A crate that ships these templates should label each by its source assumption, not just its formula.

A convention trap matters for any code that converts between Ω_GW, S_h and h_c. The 2025 review switched to **two-sided** power spectral densities with a 4π angular normalisation, so its h_c is **√2 smaller** than the Moore/Maggiore h_c for the same Ω. The two prefactors are 8.934×10⁻¹⁹ and 1.263×10⁻¹⁸ (computed here) ([Review 2025 Eqs. 9–13](https://arxiv.org/abs/2501.11723)). The review's own Eq. 11 states S_h^Review = S_h^Moore, which contradicts its Eq. 9 unless the polarisation normalisation is counted differently. The notes could not resolve this without the typeset journal version. The safe design is an explicit `Convention` enum, Ω h² as the quantity of record, and round-trip tests. The review also warns that its first edition's combined h_c plot "contains many hidden variables … which may lead to misleading conclusions", and it no longer overlays sources and detectors on a single h_c axis ([Review 2025 §3, footnote 10](https://arxiv.org/abs/2501.11723)).

## The first 2024–2026 hardware limits sit 18 to 31 orders above the ceiling

The period 2024–2026 produced the first dedicated HFGW analyses from laboratory hardware. Every one of them is a genuine measurement, and every one is far from any source:

- **QUEST** (Cardiff): twin co-located power-recycled Michelsons. It reached about **3×10⁻²⁰ Hz⁻¹ᐟ² at 40 MHz** and published "the first broadband constraints for a stochastic gravitational wave background" over 13–80 MHz ([Patra et al., PRL 2025](https://arxiv.org/abs/2410.09175)).
- **ABRACADABRA-10cm**, an axion haloscope repurposed with a figure-8 pickup: six days of December 2023 data over 10 kHz–5 MHz. Its noise matched theory. Its transient search excludes primordial-black-hole sources only **within 46 km**, a strain of order 10⁻⁴ ([Pappas et al.](https://arxiv.org/abs/2505.02821)).
- **CAPP-12T reanalysis:** h₀ ≈ **3.9×10⁻²¹** (90% CL) in a 2 MHz span at 5.311 GHz ([Kim et al. 2025](https://arxiv.org/abs/2511.17817)).
- **Laser-cavity spectroscopy** (July 2026 preprint): it constrains "much of" 100 kHz–100 MHz for the first time ([Budker et al.](https://arxiv.org/abs/2607.12617)).
- **MAGO**, a superconducting RF cavity from twenty years ago: revived by DESY and Fermilab. Its July 2026 cryogenic tests found mechanical Q factors "well below commonly assumed theoretical values", a direct hit on its sensitivity projections ([Dokuyucu et al.](https://arxiv.org/abs/2607.02476)).
- **Levitated Sensor Detector:** still "in development" in April 2026, with no published limits ([APS 2026 abstract](https://meetings-archive.aps.org/smt/2026/apr-v91/7/)).

Converted to energy density with the one-sided relation h²Ω = 6.27×10³⁵ (f/Hz)² h_c² (computed here), these limits sit absurdly high. The Holometer's 2.4×10⁻¹⁹ at 13 MHz is Ω h² ~ 6×10¹², and QUEST's 40 MHz floor is Ω h² ~ 4×10¹⁹. Against the 10⁻⁶ ceiling that is **18–31 orders of magnitude**. Any honest landscape plot must show this gap ([Holometer, Chou et al. 2017](https://arxiv.org/abs/1611.05560); [Ejlli et al. 2019](https://arxiv.org/abs/1908.00232)).

The most-publicised HFGW "detection" claim has not survived scrutiny. In 2021 a cryogenic quartz bulk-acoustic-wave (BAW) resonator at the University of Western Australia recorded two rare events at 5.5 MHz over 153 days. The press release framed them as possible gravitational waves ([Goryachev et al., PRL 127, 071102](https://arxiv.org/abs/2102.05859); [UWA](https://www.uwa.edu.au/news/Article/2021/August/World-first-detector-reports-rare-events)). Lasky and Thrane showed that a gravitational-wave reading would need Ω_GW ≈ 10⁸. Any source that loud would also have produced a memory burst in LIGO/Virgo with SNR ≳ 10⁶, and none was seen ([Lasky & Thrane 2021](https://arxiv.org/abs/2110.13319)). The 2025 review concludes the events "are not considered to be viable GW candidates". The two-detector MAGE follow-up, at 6.6×10⁻²¹ strain/√Hz in narrow bands, was built partly to find their true origin ([Campbell et al. 2023](https://arxiv.org/abs/2307.00715)).

### Magnetic conversion has a clean SI formula in which G cancels

The physics most worth encoding is the (inverse) Gertsenshtein effect. A gravitational wave crossing a transverse magnetic field B over a length L ≫ λ converts to photons with probability 4πGB²L² in natural Heaviside–Lorentz units ([Review 2025 Appendix A](https://arxiv.org/abs/2501.11723)). In SI this is **P = 4πG B²L²/(μ₀c⁴) = 8.26×10⁻³⁸ (B/1 T)²(L/1 m)²**, which gives 5.7×10⁻³⁴ for CAST and 2.6×10⁻³² for ALPS II (computed here).

Multiplying P by the wave's energy flux cancels G. For a monochromatic wave the photon rate becomes Ṅ = εAB²L²ωh₀²/(8μ₀cħ) = 1.98×10³¹ s⁻¹ × ε h₀² (f/Hz)(B/T)²(L/m)²(A/m²). For a flat stochastic background it is Ṅ = 3.95×10³¹ s⁻¹ × ε h_c² (Δf/Hz)(B/T)²(L/m)²(A/m²) (both derived and checked here). Inverting the stochastic rate for one photon per second gives exactly Ejlli et al.'s published prefactor, **h_c^min ≃ 1.6×10⁻¹⁶ (1 T/B)(1 m/L)√[(Ṅ/Hz)(m²/A)(Hz/Δf)/ε]**. That fixes their h_c as the one-sided Moore convention ([Ejlli et al. 2019](https://arxiv.org/abs/1908.00232)). The match is a ready-made unit test.

A contested correction belongs in the code as an explicit parameter rather than a silent default. Liu, Ren and Zhang argue that the light-shining-through-wall and helioscope limits must be degraded by the instruments' tiny field of view when read as limits on an isotropic background. On their reading the optical limit falls from 6×10⁻²⁶ to about 4×10⁻²³, and the X-ray limit from 5×10⁻²⁸ to 8×10⁻²⁶ ([Liu, Ren & Zhang 2024](https://arxiv.org/abs/2305.01832)).

The cavity and haloscope formulas from the review are equally compact once SI factors are restored (each B² brings a 1/μ₀). For a microwave cavity, **S_h ≃ 2k_BT_sys μ₀c²/[ω³η²QB₀²L⁵]**. With the ADMX benchmark (0.8 GHz, η = 0.1, Q = 8×10⁴, 7.5 T, 0.51 m, 0.6 K) this gives **√S_h ≈ 1×10⁻²² Hz⁻¹ᐟ²** (computed here). That matches the review's "h₀ ~ 10⁻²²–10⁻²¹ at GHz" once integration time is included ([Review 2025 Eq. 113, Table 3](https://arxiv.org/abs/2501.11723)). The resonant LC haloscope (Eq. 117), the broadband SQUID readout (Eq. 118) and photon counting (Eqs. 125–126) follow the same pattern. Their parameters for about fifteen named experiments are tabulated in the review's Tables 3 and 4. The machine-readable versions of the review's curves are the HFGWPlotter datasets on Zenodo ([noise-equivalent strain](https://doi.org/10.5281/zenodo.15720342); [stochastic signals](https://doi.org/10.5281/zenodo.15720443)). Their licence and column conventions were not verified.

Two detector-side claims fall outside the mainstream. The Li–Baker "heterodyne enhancement" detector rests on technology that does not exist and has no results, and both Woods and the 2008 JASON panel raised concerns. The review also rejects claims that an AC–DC cross-term adds sensitivity, stating that linear-signal sensitivity is independent of the background field ([Review 2025 §5.3.1, §5.5.4](https://arxiv.org/abs/2501.11723)). Chiao's superconducting "impedance-matched" antennas, speculated at h₀ ~ 10⁻³⁰, rest on assumptions Kiefer and Weber called into doubt ([Kiefer & Weber](https://arxiv.org/abs/gr-qc/0408010)).

## Gravity beams are real near fields, ten orders below any sensor

### Light and particle beams carry a tube of gravity, not a projected ray

The serious core of the "gravity beam" idea is ninety years old and fully established. Tolman, Ehrenfest and Podolsky showed in 1931 that parallel light rays do not attract while counter-propagating ones do ([INSPIRE](https://inspirehep.net/literature/46619)). Bonnor found exact Einstein solutions for steady light beams in 1969; they are plane-fronted waves ([Bonnor beam](https://en.wikipedia.org/wiki/Bonnor_beam)). Faraoni and Dumse explained the cancellation: for parallel beams the gravitoelectric pull is "exactly cancelled by a gravitomagnetic force" ([Faraoni & Dumse 1999](https://arxiv.org/abs/gr-qc/9811052)). Rätzel, Wilkens and Menzel worked out the linearised field of a laser pulse in 2016. A test mass at rest at distance ρ feels **a ≈ 4GP/(c³ρ)**, and only while the pulse passes. All effects travel on spherical shells expanding at c from the emission and absorption events ([Rätzel et al. 2016](https://arxiv.org/abs/1511.01023)).

The factor 4G, rather than the 2G of a Newtonian line mass, is the "light gravitates twice as strongly" result from the pressure term. For a petawatt at 1 mm the acceleration is **9.9×10⁻¹⁸ m/s²**. A 100 J femtosecond pulse at 1 mm delivers a velocity kick of about 10⁻³⁰ m/s (both computed here).

Spengler, Rätzel and Braun did the most careful feasibility study. Their curvature is R₀ₓ₀ₓ ≃ −4GP/(c⁵ρ²), and they stress that "there is no gravitational wave generated as this type of source is not quadrupolar in nature. Rather one can detect the gravitational near-field." The best candidate source is not a laser but the LHC beam, with 3.8×10¹² W average power and a 16 µm waist. An optimised oscillator reaches S/N ≈ 0.6 in a week. With the high-luminosity upgrade and tungsten test masses it reaches S/N ≈ 16, "neglecting all sources of technical noise". The authors call even that "a very difficult experiment" ([Spengler et al. 2022](https://arxiv.org/abs/2104.09209)).

Spin and orbital angular momentum add frame-dragging terms of about **10⁻²⁹ m/s²** for a petawatt beam, "too small to be detected with current technology" ([Schneiter et al. 2018](https://arxiv.org/abs/1804.08706)). The main 2023 twisted-light paper computes h ~ 1.9×10⁻³⁶ for a petawatt at 1.36 µm. Its authors judge the chance of detecting light's gravity "extremely low" ([Atonga et al.](https://arxiv.org/abs/2309.04191)). The genuine laboratory frontier is not beams at all. Brack et al. drove a high-Q beam through gravity alone from a transmitter vibrating at 42 Hz, measured picometre displacements, and matched Newtonian modelling to about 2% ([Nature Physics 2022, via APS Physics](https://physics.aps.org/articles/v15/121)). In short, a beam of light or matter does carry gravity, but as a logarithmic potential confined to the tube it occupies. Co-moving objects feel nothing, and nothing is projected ahead of the pulse.

### Laboratory wave generation fails by the factor G/c⁵

Making true gravitational *radiation* in a lab is worse than detecting a beam's near field. The quadrupole formula for a rigid body spun about a principal axis gives **P = (32/5)(G/c⁵)(I₁−I₂)²Ω⁶**. The classic Misner–Thorne–Wheeler steel beam (490 t, 20 m, spun to breakup at 28 rad/s) radiates **2.3×10⁻²⁹ W** at a wavelength of 34,000 km (computed here; parameters from [Chiao et al.](https://arxiv.org/abs/1712.08680)). Every table-top experiment below about 1 GHz therefore sits deep in the near field, where Newtonian coupling dominates radiation by powers of λ/r. A genuine "gravitational Hertz experiment", with generator and receiver in each other's far field, needs f ≳ 1–10 GHz. That leaves only electromagnetic schemes.

The best of these is Kolosnitsyn and Rudenko's Gertsenshtein generator: a 10 T, 100 m superconducting Fabry–Pérot cavity pumped at 30 GHz. It reaches h ~ 10⁻³²–10⁻³⁰ and about 14 detected photons per second. The authors call this a "potential opportunity" that assumes no losses and ignores instrumental noise ([Kolosnitsyn & Rudenko 2015](https://arxiv.org/abs/1504.06548)). The pump is about 5×10²⁶ photons/s at the same frequency as a signal of tens of photons per second. A real experiment would therefore need roughly 260–270 dB of isolation between generator and receiver, which the paper does not address (inference from their numbers).

Collimation obeys ordinary diffraction, θ ~ λ_GW/D. A beam therefore needs a coherent source much larger than its wavelength, which is again why every "beam" proposal clusters at GHz and above. Gravitational-wave lensing is established theory, but LIGO–Virgo–KAGRA found no confident lensed event in O3 ([LVK 2023](https://arxiv.org/abs/2304.08393)).

The "gravitational-wave laser" has a serious form and a fringe form. The serious form is astrophysical. Black-hole boson clouds supply enormous occupation numbers, so stimulated emission is plausible there ([Liu 2024](https://arxiv.org/abs/2401.16096)). Schützhold's 2025 PRL proposes detecting stimulated graviton emission or absorption by light, but he notes a semiclassical interaction gives the same phase, so it would not prove gravitons exist ([APS Physics Viewpoint](https://physics.aps.org/articles/v18/175)). Laboratory "gasers" built from high-Tc superconductors (Fontana) or from Chiao's dynamical-Casimir cavities rely on an unproven "gravitational permeativity". Chiao's 2003 microwave transducer test found no signal, set only an upper limit of 1.6×10⁻⁵ on conversion efficiency, and was never published in a journal ([Gallerati, Modanese & Ummarino 2022](https://arxiv.org/abs/2206.07574)).

### Gravitomagnetism is the measured version of "moving mass makes a field"

The velocity-dependent, magnetic-like part of gravity that cancels the attraction between parallel beams is real and measured. Gravity Probe B measured frame dragging of **−37.2 ± 7.2 mas/yr** against GR's −39.2, and a geodetic drift of −6601.8 ± 18.3 against −6606.1 ([Everitt et al. 2011](https://arxiv.org/abs/1105.3456)). In 2026 LARES-2 combined with LAGEOS and GRACE gravity fields measured Earth's Lense–Thirring node drag to about **one part in a thousand**, an order of magnitude better than before ([Ciufolini et al., Nature 2026](https://www.nature.com/articles/s41586-026-10715-0); [phys.org](https://phys.org/news/2026-07-validate-einstein.html)). Iorio disputes the error budget, arguing that J₂ biases enter through the other parameters, and the LARES team replies that his claims are "not reproducible in any independent analyses" ([Iorio, EPJC 2025](https://link.springer.com/article/10.1140/epjc/s10052-025-13964-x)).

The node precession is Ω̇_LT = 2GJ/(c²a³(1−e²)^(3/2)), which gives 30.7 mas/yr for LAGEOS (computed here). It uses the gravity-field coefficients the project already reads. That makes it the honest "gravity beam" module for `geo5d`: real gravitomagnetism, measured with the same GRACE fields the crates process.

### Five fringe claims share one failure pattern

The claims that inspired much popular "gravity beam" interest all report effects at the 10⁻⁴–10⁻² level, with no GR mechanism, that vanished as controls improved or were never reproduced:

- **Podkletnov's rotating YBCO disk** claimed weight losses of 0.05–0.5%. Hathaway et al. found nothing to 0.001%. NASA Marshall saw at most 3–5×10⁻⁵ cm/s², "significantly lower" than the claimed 2.1%, and it coincided with magnetic leakage ([Gallerati et al. 2022](https://arxiv.org/abs/2206.07574); [NASA NTRS](https://ntrs.nasa.gov/api/citations/19990064220)). His "impulse gravity generator" and its "64c" timing come only from the inventor and his co-author ([Podkletnov & Modanese](https://arxiv.org/abs/physics/0108005v2)).
- **Tajmar's gravitomagnetic London moment** claimed B_g ≈ 1.84×10⁻⁴ ω. GR frame dragging inside a 1 kg, 10 cm rotating shell is about 10⁻²⁶ ω, so the claim is 22 orders of magnitude above GR (computed here). A Canterbury ring-laser test was null, and Tajmar's own group later attributed the signals to acoustic noise ([TU Dresden](https://tu-dresden.de/ing/maschinenwesen/ilr/rfs/forschung/forschungsfelder/raumfahrtantriebe-und-neue-konzepte/breakthrough-propulsion-physics)).
- **Ning Li's superconductor theory** gave fields "too large by many orders of magnitude" ([Harris 1999](https://dx.doi.org/10.1023%2FA%3A1021621425670)).
- **The Li–Baker HFGW programme** was reviewed by the US government's JASON panel in 2008, which concluded that "No foreign threat in HFGW is credible" for communication, imaging or propulsion ([FAS summary](https://fas.org/publication/jason_study)).

Even the most optimistic peer-reviewed generator, at about 10⁻¹² W of gravitational waves, would give a thrust of P/c ≈ 3×10⁻²¹ N (inference).

| Class | Topics |
|---|---|
| **Established**: standard GR, measured or exact | GW propagation and polarisation; BBN/CMB ceiling; beam near fields (Tolman, Bonnor, Rätzel); quadrupole radiation; Gertsenshtein conversion; GW lensing theory; gravitomagnetism (GP-B, LARES-2); dynamic near-field gravity (Brack 2022); Newtonian noise; prompt elasto-gravity signals |
| **Speculative but serious**: mainstream journals, untested | PBH and superradiance sources; phase transitions and cosmic strings at MHz–GHz; every HFGW detector projection (LSD, MAGO, BAW arrays, magnetic Weber bars, Rydberg, quantum cyclotron); the LHC near-field experiment; Gertsenshtein Hertz experiments; stimulated graviton detection; astrophysical boson-cloud lasing; the "development path" to √S ~ 10⁻³⁸ |
| **Contested** | field-of-view correction to photon-regeneration limits; the LARES-2 error budget (Iorio) |
| **Fringe or refuted** | BAW events as GWs; Podkletnov shielding and impulse generator; Tajmar anomalies; Ning Li "AC gravity"; Li–Baker propulsion and imaging; lab gasers via "gravitational permeativity"; AC–DC cross-term sensitivity gains |

## Gravimetry and gravitational-wave instruments share their noise physics

### An atom gradiometer becomes a wave detector by feeding it a tidal acceleration

`satgrav` already models the Mach–Zehnder atom interferometer as Δφ = 2k_eff a T² with projection noise 1/(C√N). Extending it to gravitational waves needs almost no new physics. A wave of strain h produces a differential tidal acceleration **Δa = ½ω²hL** across a baseline L. Fed through the interferometer's acceleration transfer function k_eff(4/ω²)sin²(ωT/2), it gives **Δφ = 2·n·k·h·L·sin²(ωT/2)**. That is exactly the low-frequency, single-diamond limit of Graham et al.'s full response:

Δφ = k_eff hL [sin(ωQT)/cos(ωT/2)] sinc(ωnL/2c) sin(ωT/2 − ω(n−1)L/2c)

On resonance it peaks at 2Qk_eff hL ([Graham, Hogan, Kasevich & Rajendran 2016](https://arxiv.org/abs/1606.01860)). The MIGA paper gives the matching strain spectrum with gravity-gradient noise: S_h = 4S_Δa/(ω⁴d²) + S_ε/(2(2nk)²d² sin⁴(ωT/2)) ([Junca et al. 2019](https://arxiv.org/abs/1902.05337)). The one physics change is the wavenumber. Long-baseline detectors use the single-photon ⁸⁷Sr clock transition at 698 nm with k = ω_A/c, because Rb Raman laser phase noise does not cancel over long baselines. A sibling of `K_EFF_RB87` is required, and Rb parameters cannot simply be scaled to a kilometre.

The detectors themselves are still a decade away. **No atom-interferometer gravitational-wave detector operates in 2026**:

- **MAGIS-100** (Fermilab, 100 m): laser-lab construction finished in January 2026, with commissioning expected by 2028 ([Fermilab](https://news.fnal.gov/tag/magis-100/)).
- **AION-10** (Oxford, 10 m): it has a 2025 technical design report and "is not expected to detect GWs" ([AION-10 TDR](https://arxiv.org/abs/2508.03491); [Badurina et al. 2021](https://arxiv.org/abs/2108.02468)).
- **AICE**, a 100 m interferometer in CERN's PX46 shaft: a 2025 letter of intent, costed at 30–50 MCHF ([Baynham et al. 2025](https://arxiv.org/abs/2509.11867)).
- **TVLBAI roadmap:** km-scale detectors operating from 2035 ([TVLBAI 2025](https://arxiv.org/abs/2503.21366)).

The decisive hardware gap is large momentum transfer (LMT). The demonstrated figure is about **400ħk** ([Wilkason et al. 2022](https://arxiv.org/abs/2205.06965)), while the AICE and AION designs assume **5,000–40,000** ([AION, Badurina et al. 2020](https://arxiv.org/abs/1911.11755)).

Evaluating h ≈ δφ/(2nkL) on the published design tables gives clean unit tests (computed here, optimal orientation):

| Design | Peak strain (Hz⁻¹ᐟ²) |
|---|---|
| AION-10 initial | 5.6×10⁻¹⁴ |
| AION-100 initial | 5.6×10⁻¹⁷ |
| AICE | 2.0×10⁻¹⁸ |
| AION-100 goal | 1.4×10⁻¹⁹ |
| AION-km | 2.1×10⁻²¹ |

AICE's target of 10⁻⁵ rad/√Hz with 10⁸ atoms per 1.5 s cycle is about 12× below the single-cloud standard quantum limit. It implies roughly 20 dB of squeezing or interleaved clouds, so the crate should take δφ as an input or expose a squeezing factor rather than derive it from N.

The same algebra turns any gradiometer into a strain detector: **h = 2√S_Γ/ω²**. A GOCE-class 10 mE/√Hz gradiometer gives 5×10⁻⁹ at 10 mHz, 5×10⁻¹¹ at 0.1 Hz and 5×10⁻¹³ at 1 Hz (computed here). That is 8–12 orders of magnitude short of LISA or AION-km. A CARIOQA- or GRICE-type space gradiometer is not a gravitational-wave detector; what matters is baseline, which is why the AEDGE space design separates its atom clouds by 4.4×10⁷ m ([El-Neaj et al. 2020](https://arxiv.org/abs/1908.00802)).

GRACE-FO makes the point with real data. Reading its laser-ranging displacement over the 220 km separation as strain, h = x/x₀, gives Ω_GW < 3.3×10⁷ at 100 mHz and a best characteristic strain of 2×10⁻¹⁵ at 0.5 Hz. This is a real bound, "many orders of magnitude" above expected sources ([Ross et al. 2020](https://arxiv.org/abs/2002.02044)). For the strain view, LISA is the natural reference curve. ESA adopted it on 25 January 2024, with launch planned for 2035 ([ESA](https://www.esa.int/Science_Exploration/Space_Science/LISA/Capturing_the_ripples_of_spacetime_LISA_gets_go-ahead)). The Robson–Cornish–Liu analytic curve is two noise terms, a transfer factor and a four-row confusion table ([Robson et al. 2019](https://arxiv.org/abs/1803.01944)). Optical-lattice clocks on two satellites reach a shot-noise fractional frequency of 1.1×10⁻²⁰/√Hz with 7×10⁶ Sr atoms, four orders beyond the clock stability demonstrated in 2016 ([Kolkowitz et al. 2016](https://arxiv.org/abs/1606.01859)). Clocks and interferometers share the same atom-projection limit. A shared atom-ensemble type in `satgrav` could therefore feed both.

### Newtonian noise is gravimetric forward modelling, and earthquakes outrun their own P waves

Gravity-gradient or "Newtonian" noise in gravitational-wave detectors is the same physics as forward-modelling density changes in gravimetry. Seismic waves and infrasound move mass, and no shield can block the resulting gravity. For an isotropic Rayleigh-wave field at depth h, the MIGA formalism gives the differential acceleration PSD S_Δa = κ_R² S_ξz χ_R, with κ_R = 2πGρ₀γ[e^(−hk_R) + b(e^(−hq_z^P) − e^(−hk_R))] and χ_R = 1 − J₀(k_Rd) + J₂(k_Rd) ≈ (3/8)k_R²d² at low frequency ([Junca et al. 2019](https://arxiv.org/abs/1902.05337)). Combined with the strain formula, the low-frequency seismic strain noise is **S_h ≈ (3/2)κ_R²S_ξz/(ω²c_R²), independent of baseline** (derived here). This is a compact function that takes Peterson's low-noise or high-noise seismic model as input.

The mitigation algebra is equally transferable. For a small cavity in a homogeneous medium the noise is C_NN = (4πGρ₀/3)²S_ξ(3p+1), and a Wiener filter from N sensors at signal-to-noise ratio SNR leaves a residual of about 1/(N·SNR²). About 15 sensors give a tenfold amplitude reduction ([Badaracco & Harms 2019](https://arxiv.org/abs/1903.07936)), and 2025 "fusion" arrays of seismometers and fibre strainmeters match 20 borehole seismometers ([Ophardt et al. 2025](https://arxiv.org/abs/2512.13554)). The same estimator applied to a barometer array is the dynamic generalisation of the barometric admittance correction applied to superconducting gravimeter records. One code path can therefore serve both a ground-gravimeter correction and a detector noise budget.

Prompt elasto-gravity signals (PEGS) are the geophysical payoff. A rupture changes the mass distribution instantly, and the gravity change arrives at light speed, before the P wave. Signals of about **1 nm/s²** from Tohoku were seen clearly on broadband seismometers 1000–2000 km away ([Vallée et al. 2017](https://www.science.org/doi/10.1126/science.aao0746)). Since then deep learning has tracked earthquake growth from PEGS in real time ([Licciardi et al. 2022](https://www.nature.com/articles/s41586-022-04672-7)), and Peru implemented the first PEGS-based tsunami-warning system in 2025 ([GJI 244, ggaf419](https://academic.oup.com/gji/article/244/1/ggaf419/8300379)). A TOBA torsion bar at its Phase-III target of 10⁻¹⁵/√Hz at 0.1 Hz could detect an Mw 7 at 100 km about 10 s earlier than current warnings ([Takano et al. 2024](https://arxiv.org/abs/2412.01323)). That target equals √S_Γ ≈ 2×10⁻¹⁶ s⁻²/√Hz, about 5×10⁴ times better than GOCE (computed here). Torsion bars, not space gradiometers, are therefore the PEGS instruments.

No satellite PEGS detection has been published. A GRACE-FO detection is plausible only for Mw ≳ 9 under a fortunate overpass, and would be marginal. That is a hypothesis, not a finding, and the exact Harms 2015 early-time formula could not be retrieved to test it.

### The gravity waves already in our satellite data are atmospheric

The word "gravity waves" also names buoyancy waves in the atmosphere, and these are the only gravity-related waves the project's satellites have actually recorded. GOCE, at about 270 km, recorded Tohoku's infrasound directly in 2011: air-density changes of about **11%** and vertical accelerations of about **10⁻⁷ m/s²**, matching a forward model ([Physics Today](https://physicstoday.aip.org/news/a-seismometer-in-orbit); [ESA](https://www.esa.int/Applications/Observing_the_Earth/GOCE/GOCE_the_first_seismometer_in_orbit)). The 2022 Hunga Tonga eruption sent a Lamb wave around the globe at **318 m/s** ([Wright et al., Nature 2022](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9492537/)). GRACE-FO and Swarm-C accelerometer densities show three thermospheric waves at 200–450 m/s near 500 km ([Li et al. 2023](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2022GL102265)).

The cause of those waves is a live 2025 dispute, which makes it a good test case. MESORAC-HIAMCM modelling attributes them to secondary gravity waves and finds that the Lamb wave is much weaker at 510 km ([Li et al. 2025](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2024AV001470)). WACCM-X gives the Lamb wave about 25% of the disturbance ([Liu et al. 2023](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2023GL103682)). On the operational side, NASA JPL's GUARDIAN system flagged the 29 July 2025 Kamchatka tsunami in ionospheric GNSS data about **30 minutes before landfall** in Hawai'i ([Natural Hazards 2026](https://ideas.repec.org/a/spr/nathaz/v122y2026i9d10.1007_s11069-026-08151-4.html); [JPL](https://www.jpl.nasa.gov/news/see-nasas-guardian-catch-a-tsunami/)).

The governing physics is the textbook acoustic-gravity dispersion relation, k_z² = (ω² − ω_a²)/c² − k_h²(1 − N²/ω²). It has an isothermal acoustic cut-off ω_a = γg/(2c) and a Brunt–Väisälä frequency N = √(γ−1)·g/c ([Ann. Geophys. 2024](https://angeo.copernicus.org/articles/42/491/2024/); [arXiv:1705.01388](https://arxiv.org/pdf/1705.01388)). At 250 K these give c = 316.9 m/s, ω_a = 0.0217 rad/s (4.8 min) and N = 0.0196 rad/s (5.3 min) (computed here). That explains the 3–5 min co-seismic acoustic periods.

The untested question is whether the Tonga Lamb front shows up in GRACE-FO laser-ranging residuals on 15–16 January 2022. The AOD1B RL07 de-aliasing product is built from 3-hourly ERA5 ([GFZ AOD1B RL07](https://isdc-data.gfz.de/grace/DOCUMENTS/Level-1/GRACE_AOD1B_Product_Description_Document_for_RL07.pdf)), so a sharp pressure front lasting minutes would be under-sampled and would leak into those residuals. A detection would be new, and a null would be a publishable bound. No source has tested it, and it needs Level-1B data the crates do not yet read.

## Build one strain crate, one satgrav layer and one frame-dragging module

### Waveforms are covered in Rust, so build sensitivity tooling

The Rust gravitational-wave ecosystem changed in the week this research was done:

- **waveformrs** 0.2.0, published 22 October 2026 by LVK developer Colm Talbot, ships TaylorF2, IMRPhenomD and XPHM, NRTidal and NRSur7dq4, "checked against LALSimulation". It is **GPL-3.0** because much of it is translated from LALSuite ([crates.io](https://crates.io/crates/waveformrs)).
- **bilby_rust** (antenna response, sidereal time), **distmarg** (distance-marginalised likelihood) and **igwn-ligolw** (LIGO_LW XML) cover other pieces ([bilby_rust](https://github.com/ColmTalbot/bilby_rust); [distmarg](https://crates.io/crates/distmarg); [igwn-ligolw](https://crates.io/crates/igwn-ligolw)).
- **libgrav** has binary-parameter conversions but is not on crates.io ([libgrav](https://github.com/transientlunatic/libgrav)).

The gap is the sensitivity and forecast layer, the least-maintained part of the Python stack too: gwent's last release was 2023, hasasia's 2021 and pygwb's 2024 ([gwent](https://pypi.org/project/gwent/); [hasasia](https://pypi.org/project/hasasia/); [pygwb](https://pypi.org/project/pygwb/)). For HFGWs the only tooling is a web plotter at Mainz built on the review's data ([Shplot](https://incandenza-01.zdv.uni-mainz.de/Shplot/)). No Rust or WASM sensitivity-curve or SNR crate turned up. Writing every formula from the papers, rather than translating LAL code, keeps the result MIT like the rest of the workspace. waveformrs and LALSuite then serve as external test oracles, never as dependencies.

The design should follow `satgrav`'s existing pattern of passing noise curves as `impl Fn(f64) -> f64`. `gwsens` and `satgrav` then need no dependency on each other: a `satgrav` interferometer ASD is just a closure that `gwsens` integrates.

### Priority order

| # | Build | Home | Contents | Gate before the next step |
|---|---|---|---|---|
| 1 | Strain conventions, curves and SNR | new crate `gwsens` | `Convention` enum (one-sided Moore/Thrane–Romano versus two-sided LRR 2025); Ω ↔ S_h ↔ h_c with H₀ as a parameter; Sathyaprakash–Schutz aLIGO and ET fits; Robson LISA with confusion table; log-log interpolator for tabulated ASDs; inner product, optimal SNR, Newtonian horizon and range (÷2.264) | Round-trip conversions to 10⁻¹²; LISA, Hellings–Downs, BNS and NANOGrav fixtures below; Chen et al. horizons with the O4high file to about 1% |
| 2 | Stochastic background and pulsar timing | `gwsens::sgwb` | power laws (β = 2α + 2), network S_eff, power-law-integrated curves (β = −8…8), Hellings–Downs, γ = 3 − 2α | PI curve matches the Thrane–Romano reference code (LIGO-P1300115) |
| 3 | High-frequency module | `gwsens::hf` | BBN, CMB and future-CMB ceiling lines; Gertsenshtein probability, coherent and stochastic photon rates, Ejlli inversion, optional field-of-view factor flagged as contested; cavity, LC, SQUID and photon-counting noise in SI; detector registry with a `Status` tag (published limit, recast, prototype, projection, disputed) and a citation per entry; source templates clipped to the ceiling | Gertsenshtein, Ejlli, ADMX, f_ISCO, superradiance and S(1) fixtures below |
| 4 | Lab-gravity reality check | `gwsens::lab` | beam near-field acceleration, curvature and velocity kick; quadrupole luminosity for rods and rotors; GW flux F = c³ω²h²/(32πG) | Beam, LHC and MTW-rod fixtures below |
| 5 | GW layer for atom sensors | `satgrav::gw` | Sr-87 clock wavenumber; Δa = ½ω²hL; Graham Eq. 9 response with Q and finite light time; MIGA strain with S_Δa; h = 2√S_Γ/ω²; h = x/x₀ for satellite pairs; Rayleigh Newtonian-noise closed form; 4/15 sky average | AION/AICE strain table, GOCE, TOBA and clock fixtures below |
| 6 | Frame dragging | `geo5d` relativity | Lense–Thirring node rate; Schiff gyroscope precession; orbit averages | LAGEOS and GP-B fixtures below |
| 7 | Landscape figure | `apps/gravviz` | two panels: √S_n(f) for detectors against coherent sources, and Ω h²(f) with PI curves and the ceiling band | Reproduces the review's Fig. 1 and Fig. 4 shapes from Zenodo, after a licence check |
| 8 | Pre-registered test | gravimetry hypothesis log | acoustic-gravity dispersion helpers; Tonga Lamb-front timing against GRACE-FO laser-ranging residuals | 250 K dispersion fixtures; the 318 m/s Lamb speed; the protocol registered before the data are opened |

Items 1–4 need no network data and are pure formula work with closed-form tests. Item 5 is under a few hundred lines because it reuses `satgrav`'s transfer function. Item 8 is the only one that could produce a new scientific result, and the only one that needs new data readers.

### Formulas and test values

The values below were recomputed in double precision for this report from the formulas in the cited sources. Test the formulas to about 10⁻⁶ relative accuracy and the published comparisons to their stated tolerance.

| Quantity | Formula | Test value | Source |
|---|---|---|---|
| H₀ per unit h | 100 km s⁻¹ Mpc⁻¹ | 3.2408×10⁻¹⁸ s⁻¹ | standard |
| Strain from energy density, one-sided | h_c = √(3H₀²Ω/2π²)/f | 1.2634×10⁻¹⁸ (1 Hz/f)√(h²Ω); LRR 2025 convention ÷√2 = 8.934×10⁻¹⁹ | [Moore et al.](https://arxiv.org/abs/1408.0740); [Review 2025](https://arxiv.org/abs/2501.11723) |
| BBN/CMB ceiling | Ω h² < 1.1×10⁻⁶ | h_c < 1.33×10⁻³⁰ at 1 GHz; check 5.6×10⁻⁶ × ΔN_eff (gives 1.0×10⁻⁶ at 0.18) to 10% | [Review 2025 Eq. 86](https://arxiv.org/abs/2501.11723) |
| Gertsenshtein probability | 4πGB²L²/(μ₀c⁴) | 8.263×10⁻³⁸ (1 T, 1 m); CAST (9 T, 9.26 m) 5.74×10⁻³⁴; ALPS II (5.3 T, 106 m) 2.61×10⁻³² | [Review 2025 App. A](https://arxiv.org/abs/2501.11723) |
| Stochastic photon rate | εAB²L²h_c²·2πΔf/(4μ₀cħ) | 3.954×10³¹ s⁻¹ per unit; one photon/s gives h_c = 1.590×10⁻¹⁶ | [Ejlli et al. Eq. 4](https://arxiv.org/abs/1908.00232) |
| Coherent photon rate | εAB²L²ωh₀²/(8μ₀cħ) | 1.977×10³¹ s⁻¹ × εh₀²(f/Hz) | derived from the above |
| Cavity noise (SI) | 2k_BT_sys μ₀c²/(ω³η²QB²L⁵) | ADMX: √S_h = 9.7×10⁻²³ Hz⁻¹ᐟ² | [Review 2025 Eq. 113, Table 3](https://arxiv.org/abs/2501.11723) |
| PBH ISCO | c³/(6^(3/2)πGM) | 4397 Hz for 1 M_⊙ (review rounds to 4400) | [Review 2025 Eq. 46](https://arxiv.org/abs/2501.11723) |
| Superradiance line | f = 2μc²/h | 4.836×10¹⁴ Hz per eV | [Review 2025 Eq. 73](https://arxiv.org/abs/2501.11723) |
| Phase-transition shape | S(x) = x³[7/(4+3x²)]^(7/2) | S(1) = 1 exactly | [Review 2025 Eq. 101](https://arxiv.org/abs/2501.11723) |
| Beam near field | a = 4GP/(c³ρ) | 1 PW, 1 mm: 9.91×10⁻¹⁸ m/s²; LHC 3.8×10¹² W, 100 µm: 3.77×10⁻¹⁹ | [Rätzel et al.](https://arxiv.org/abs/1511.01023); [Spengler et al.](https://arxiv.org/abs/2104.09209) |
| Spinning rod | P = (2/45)GM²L⁴Ω⁶/c⁵ | 4.9×10⁵ kg, 20 m, 28 rad/s: 2.27×10⁻²⁹ W; λ = πc/Ω = 3.36×10⁷ m | [Chiao et al.](https://arxiv.org/abs/1712.08680) |
| Lense–Thirring node rate | 2GJ/(c²a³(1−e²)^(3/2)) | LAGEOS (J = 5.86×10³³, a = 12,270 km, e = 0.0045): 30.67 mas/yr | [LARES-2 context](https://www.nature.com/articles/s41586-026-10715-0) |
| GP-B orbit averages | geodetic (3/2)(GM)^(3/2)/(c²a^(5/2)); frame dragging GJ/(2c²a³) | a = 7013 km: 6638 and 41.1 mas/yr (full GR with oblateness: 6606 and 39.2) | [Everitt et al.](https://arxiv.org/abs/1105.3456) |
| Hellings–Downs | ζ = (3/2)x ln x − x/4 + 1/2, x = (1−cos ψ)/2 | 0.44250 (10°), 0.21163 (30°), −0.08236 (60°), −0.14486 (90°), −0.01114 (120°), 0.16971 (150°), 0.25 (180°) | [Thrane & Romano](https://arxiv.org/abs/1310.5300) |
| NANOGrav amplitude | Ω = (2π²/3H₀²)f²h_c² | A = 2.4×10⁻¹⁵ at 1/yr gives Ω h² = 3.62×10⁻⁹ | [NANOGrav 15-yr](https://arxiv.org/abs/2306.16213) |
| LISA instrument curve | Robson Eq. 13 | √S_n = 4.60×10⁻¹⁷, 1.28×10⁻¹⁹, 1.74×10⁻²⁰, 1.20×10⁻²⁰, 4.61×10⁻²⁰ at 0.1, 1, 3, 10, 100 mHz | [Robson et al.](https://arxiv.org/abs/1803.01944) |
| BNS 1.4 + 1.4 M_⊙ | M_c, τ₀ = (5/256)(πf)^(−8/3)(GM_c/c³)^(−5/3) | M_c = 1.21877 M_⊙; τ₀(20 Hz) = 157.9 s; ḟ(100 Hz) = 17.36 Hz/s; f_ISCO = 1570.4 Hz | [Sathyaprakash & Schutz](https://arxiv.org/abs/0903.0338) |
| Horizon | Newtonian TaylorF2 at SNR 8 | 449.4 Mpc for O4high with the tabulated file; about 496.5 Mpc with the analytic fit (smoke test only) | [Chen et al.](https://arxiv.org/abs/1709.08079) |
| Atom-interferometer strain | δφ/(2nkL), k = 2π/698 nm | AION-km 2.08×10⁻²¹; AICE 2.02×10⁻¹⁸; AION-100 initial 5.55×10⁻¹⁷ | [Badurina et al.](https://arxiv.org/abs/1911.11755); [AICE](https://arxiv.org/abs/2509.11867) |
| Gradiometer to strain | h = 2√S_Γ/ω² | 10 mE/√Hz: 5.07×10⁻⁹ (10 mHz), 5.07×10⁻¹¹ (0.1 Hz) | [Junca et al.](https://arxiv.org/abs/1902.05337) |
| Clock detector | √Δ_A/(ν√(2πτN)) | ν = 429.2 THz, Δ_A = 1 mHz, N = 7×10⁶: 1.11×10⁻²⁰/√Hz | [Kolkowitz et al.](https://arxiv.org/abs/1606.01859) |
| Acoustic-gravity cut-offs | c = √(γRT); ω_a = γg/2c; N = √(γ−1)g/c | 250 K: 316.9 m/s, 0.02167 rad/s, 0.01958 rad/s | [arXiv:1705.01388](https://arxiv.org/pdf/1705.01388) |

Several formulas carry known transcription risks:

- **Graham Eq. 12:** the telescope-size exponents were reconstructed from extracted text.
- **Infrasound coupling:** κ_I was restored by dimensional analysis.
- **Kolkowitz noise formula:** the square root was restored so that it reproduces the paper's number.
- **GRACE-FO conversion:** Ross et al.'s prefactor reads "32π³/3" in extraction against the standard 2π²/3.
- **Review Eqs. 56 and 96:** both have ambiguous exponent placement.

Each needs a check against the typeset PDF before it is coded. The CGMB spectrum needs tabulated rates from Ringwald et al. or Ghiglieri et al. 2024, because no closed-form fit was extracted.

### What to leave out

Leave out waveform models beyond the Newtonian amplitude. waveformrs already covers them, and its GPL licence would contaminate an MIT crate. Write a clean-room TaylorF2 only if a permissive WASM path is actually needed.

Leave out GWF frame-file reading, which is a large parser for little project value, and defer even the GWOSC HDF5 reader. O4a and O4b are public under CC BY 4.0 ([GWOSC O4b](https://gwosc.org/O4/O4b)), but nothing in this plan needs raw strain. Real-time detection pipelines are also out: GstLAL, PyCBC Live and MBTA own that space.

Leave out any model of a gravity-beam generator, shield, thruster or communicator. Podkletnov, Tajmar, Ning Li and Li–Baker belong in the documentation only, as a short table that sets each claimed effect against its GR value. The GR calculation is the useful part.

Do not plot the review's "development path" to √S ~ 10⁻³⁸ as a projection, do not list the 2021 BAW events as candidates, and do not draw a single h_c chart mixing sources and detectors. Defer a PEGS module until the Harms 2015 formula is in hand, and do not build a GNSS-TEC pipeline, which would duplicate VARION and GUARDIAN. Finally, no documentation should present CARIOQA-class gradiometers as gravitational-wave instruments.

## Conclusion

The honest answer and the useful software turn out to be the same artifact: one strain axis with a hard ceiling drawn on it. The user's enthusiasm for high-frequency waves and gravity beams is best served not by modelling devices that cannot work, but by a calculator that makes the gap visible. It would show, for example, a petawatt beam at 10⁻¹⁷ m/s², a lab rotor at 10⁻²⁹ W and the best MHz limit at Ω h² ~ 10¹⁹, all against the same 10⁻⁶ line. That reframing also turns a liability into a contribution. No tool in any language currently puts HFGW detectors, atom interferometers, LISA and satellite gravimetry on one axis with their status and provenance attached, and the few Python tools that come close have gone years without a release.

The deeper insight is that gravimetry and gravitational-wave detection are the same measurement read at different baselines. What disqualifies a gradiometer as a wave detector, its short arm, is what makes it a good mapper of mass. And what wave detectors must subtract as Newtonian noise, moving density, is exactly what gravimetry exists to measure. The crates' genuine research opportunities lie in that overlap: the shared transfer function in `satgrav`, frame dragging measured with the very GRACE fields the project reads, and the atmospheric waves sitting unexamined in GRACE-FO residuals. Not in the terahertz band. The events worth watching are the first levitated-sensor limits, MAGO's fight with its mechanical Q, MAGIS-100 commissioning in 2028 and the AICE decision. Each will move one row of the registry, and none is likely to move the ceiling.
