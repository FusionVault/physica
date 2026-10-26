# Toroidal moments, anapoles, super/hypertoroids, non-radiating sources, vector potentials and the Aharonov–Bohm effect, and the papers and patents Greenyer's TeslaTech 2026 talk builds on

Scope: the established physics, with formulas and published test values, plus an identification and verdict for each source the talk cites (transcript `/home/dan/Desktop/transcript_1793004282390.txt`, lines 236–330, 393–470, 549–650, 842–870, 893–902, 1026–1031, 1171–1222, 1751–1781, 1999–2008).
Conventions here: SI units, time dependence e^(−iωt). Every number below comes from a source I opened (or Crossref metadata); where a number came only from a search snippet it is marked "(snippet, not opened)". Downloaded texts are in `notes/dl_toroidal/`. Checks: `notes/toroidal_moments_checks.py` (standard library only, 21 assertions, runs in about 0.3 s, all pass).

Status labels: **mainstream** (published and reproduced), **contested** (published, disputed or never independently reproduced), **refuted** (contradicted by measurement or by an exact theorem), **unverified** (no source found, or source not opened).

## 1. The toroidal dipole moment, its radiation, the anapole, Zel'dovich and the ¹³³Cs anapole

### Takeaway
The toroidal dipole T = (1/10)∫[(r·J)r − 2r²J] dV (SI; Gaussian texts put 1/c in front) is a well-established third family of multipoles, next to the electric and magnetic ones. It radiates with the same pattern as an electric dipole, carried by an effective dipole p + (ik/c)T, so its own radiated power goes as ω⁶ where an electric or magnetic dipole goes as ω⁴. An "anapole" is the case p = −(ik/c)T, where the far field cancels. Zel'dovich proposed the anapole in a paper received in 1957 and published in 1958, not in 1954. The 1997 measurement the speaker probably means is the ¹³³Cs nuclear anapole moment (Wood et al.). It measured a static, parity-violating nuclear moment, not a "fractal toroidal moment".

### Cited Findings
- The standard review is Dubovik & Tugushev, "Toroid moments in electrodynamics and solid-state physics", Phys. Rep. 187, 145–202 (1990). Crossref-checked; full text not opened. — [Crossref DOI 10.1016/0370-1573(90)90042-Z](https://doi.org/10.1016/0370-1573(90)90042-Z)
- The other standard review is Papasimakis, Fedotov, Savinov, Raybould & Zheludev, "Electromagnetic toroidal excitations in matter and free space", Nat. Mater. 15, 263–271 (2016). Crossref-checked; the accepted-version PDF at Southampton returned an HTML block page, so the text was not read. — [doi:10.1038/nmat4563](https://doi.org/10.1038/nmat4563)
- The toroidal multipoles come from Dubovik & Cheshkov: Sov. Phys. JETP 24 (1967) and Sov. J. Part. Nucl. 5, 318–364 (1974). These are reference-list entries in the review opened next. — [Yang & Bozhevolnyi, arXiv:1810.01221](https://arxiv.org/abs/1810.01221)
- In the spherical-multipole expansion, the toroidal coefficient T_lm depends on the radial current (r·J ≠ 0), and electric and toroidal multipoles share the same vector harmonic Ψ_lm. The generalised anapole condition is Q_lm = −T_lm (that review's eq. 27). — [Yang & Bozhevolnyi 2018, arXiv:1810.01221](https://arxiv.org/abs/1810.01221)
- Nemkov, Basharin & Fedotov (PRA 98, 023858, 2018) state that a toroidal multipole can cancel a charge multipole when Q⁽⁰⁾_lm + ikT⁽⁰⁾_lm = 0, and that "at l = 1 it gives the familiar anapole condition d + ikτ = 0" (Gaussian units). They also state that mean-square-radius corrections "have exactly the same radiation pattern as their parent multipoles". — [arXiv:1804.07342](https://arxiv.org/abs/1804.07342); [Crossref](https://doi.org/10.1103/PhysRevA.98.023858)
- For poloidal current on a torus, the same paper gives a toroidal dipole "of c⁻¹π²j₀Rρ = IV/4πc, with V being the volume of the torus" (its App. B). — [arXiv:1804.07342](https://arxiv.org/abs/1804.07342)
- Miroshnichenko et al. (Nat. Commun. 6, 8069, 2015) on the radiation pattern: an oscillating toroidal dipole radiates "with a pattern fully repeating that of an electric dipole moment but being scaled by a factor ω²". — [PMC4560796](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4560796/); [doi:10.1038/ncomms9069](https://doi.org/10.1038/ncomms9069)
- The same paper on static toroids: a toroid with a constant poloidal current "generates no field outside, but may possess a non-zero potential". — [Miroshnichenko et al. 2015](https://doi.org/10.1038/ncomms9069)
- That paper reports the first optical anapole, in a silicon nanodisk 50 nm high with a diameter of 200–400 nm. — [Miroshnichenko et al. 2015](https://doi.org/10.1038/ncomms9069)
- For a sphere with refractive index n = 4 at a wavelength of 550 nm, "the diameter ∼204 nm" gives P_sph = 0 with P_car ≠ 0 (an anapole). Near "196 nm" the Cartesian dipole vanishes instead. — [Miroshnichenko et al. 2015](https://doi.org/10.1038/ncomms9069)
- The first toroidal dipolar response in a metamaterial: Kaelberer, Fedotov, Papasimakis, Tsai & Zheludev, Science 330, 1510–1512 (2010). Crossref title and metadata only. — [doi:10.1126/science.1197172](https://doi.org/10.1126/science.1197172)
- Zel'dovich's paper is cited as "ZelDovich Y B 1958 Electromagnetic interaction with parity violation Sov. Phys." That review says "In 1957, he pointed out…". — [Yang & Bozhevolnyi 2018](https://arxiv.org/abs/1810.01221)
- The JETP archive record gives Sov. Phys. JETP 6(6), 1184 (June 1958) and ZhETF 33(6), 1531, received 26 Sept 1957 (snippet, not opened). — [JETP record](https://jetp.ras.ru/cgi-bin/e/index/e/6/6/p1184?a=list)
- Haxton places the anapole after the 1957 discovery of parity violation: "Shortly thereafter… Vaks and Zeldovich independently noted that particles could then have parity-violating couplings to the electromagnetic field." — [Haxton, Science 275, 1753 (1997), arXiv:nucl-th/9704059](https://arxiv.org/abs/nucl-th/9704059)
- Wood, Bennett, Cho, Masterson, Roberts, Tanner & Wieman, "Measurement of Parity Nonconservation and an Anapole Moment in Cesium", Science 275, 1759–1763 (21 Mar 1997). Crossref-checked. — [doi:10.1126/science.275.5307.1759](https://doi.org/10.1126/science.275.5307.1759)
- Haxton on that result: "a seven-fold improvement to 0.35% in the ¹³³Cs results". The anapole parameter is "κ … is 0.127 ± 0.019, a result differing from zero by ∼7σ". The anapole moment grows as A^(2/3). — [Haxton 1997](https://arxiv.org/abs/nucl-th/9704059)
- Wood et al.'s own values, Im(E1pnc)/β = −1.5935(56) mV/cm and a nuclear-spin-dependent part of 0.077(11) mV/cm, come only from a search snippet. The abstract was not opened. — [INSPIRE record](https://inspirehep.net/literature/455239)

### Inferences
- For a divergence-free current I derived the toroidal part of ∫J_i x_j x_k dV as −2T_iδ_jk + T_jδ_ik + T_kδ_ij. So the second-order term of the far-field integral ∫J e^(−ik n·r) dV has transverse part k²T⊥, and the effective radiating dipole is **p_eff = p + (ik/c)T** (SI, e^(−iωt)). The radiated power follows:
  - P = μ₀ω⁴|p_eff|²/(12πc)
  - P_T = μ₀ω⁶|T|²/(12πc⁵) for a pure toroid.

  The check script verifies this on a discretised 24-turn toroidal solenoid at kR = 0.05:
  - far field: 0.03 % error
  - total power: 0.05 %
  - P(2ω)/P(ω) = 63.9 against 64
  - adding p = −(ik/c)T leaves 8×10⁻⁸ of the power; the opposite sign gives 4.001× the power.

  The ω² factor in the pattern quoted by Miroshnichenko et al. is the (k/c)² ratio of intensities.
- The script also computes the zero of the Mie coefficient a₁ for n = 4 at 550 nm: 204.0 nm, against the published ∼204 nm.
- The speaker's "Zel'dovich 1954" is **refuted** by the publication record: parity violation was found in 1957, and the paper was received in 1957 and published in 1958.
- "Verified in the West in 1997" fits Wood et al. 1997, which measured a static nuclear anapole moment (~7σ by Haxton's count). It is not a measurement of a "fractal" or "hypertoroidal" moment. The first dynamic toroidal dipole was seen in a metamaterial in 2010 (Kaelberer et al.).

### Gaps
- Dubovik & Tugushev 1990 and Papasimakis 2016 were not opened, so the 1/10 normalisation is checked here against Nemkov 2018's torus result, not against those reviews. The script finds that the 1/10 formula gives exactly N I V/(4π), the SI form of Nemkov's IV/(4πc).
- Radescu & Vaman were not located.
- Wood et al.'s full abstract was not opened.
- The claim that Dubovik found the toroidal moment in 1965 and put it in a 1967 Dubna thesis is unverified. The 1967 JETP paper with Cheshkov exists as a reference-list entry.

## 2. Fields of a toroidal solenoid, the Aharonov–Bohm phase and Tonomura 1986; do potentials without fields carry energy or force?

### Takeaway
An ideal toroidal solenoid confines B inside its winding and has B = 0 outside, while ∮A·dl around a loop through the hole equals the enclosed flux. An electron passing through the hole picks up the Aharonov–Bohm phase eΦ/ħ relative to one passing outside. Tonomura's group measured this with a 6 µm toroidal ferromagnet clad in superconducting niobium. Because the flux is quantised in units of h/2e, the phase came out as 0 or π. The vector potential outside does no work and exerts no classical force: the Lorentz force involves only E and B, and the AB effect is a phase shift without energy exchange.

### Cited Findings
- Aharonov & Bohm, "Significance of Electromagnetic Potentials in the Quantum Theory", Phys. Rev. 115, 485–491 (1959). Crossref-checked. — [doi:10.1103/PhysRev.115.485](https://doi.org/10.1103/PhysRev.115.485)
- Tonomura, Osakabe, Matsuda, Kawasaki, Endo, Yano & Yamada, PRL 56, 792–795 (1986). From the abstract:
  - "A toroidal ferromagnet was covered with a superconductor layer to confine the field, and further with a copper layer for complete shielding from the electron wave."
  - "The experiment gave direct evidence for flux quantization also."

  — [INSPIRE abstract / doi:10.1103/PhysRevLett.56.792](https://doi.org/10.1103/PhysRevLett.56.792)
- Osakabe et al., PRA 34, 815–822 (1986), abstract: "Below the transition temperature the relative phase shift was measured to be 0 or π due to the magnetic-flux quantization in units of h/2e." — [doi:10.1103/PhysRevA.34.815](https://doi.org/10.1103/PhysRevA.34.815)
- Hitachi's account of the experiment: "a doughnut-shaped (toroidal) ferromagnet six micrometers in diameter… covered… with a niobium superconductor". It was held at 5 K, and the fringes were displaced by "just half a fringe". — [Hitachi R&D page](https://www.hitachi.com/rd/research/materials/quantum/aharonov-bohm/index.html)
- Afanasiev & Dubovik, "Electromagnetic properties of a toroidal solenoid", J. Phys. A 25, 4869–4886 (1992). The abstract says the paper gives the multipole expansion and toroidal form factors. It also says "uniform rotation of the toroidal solenoid with a constant current… leads to the appearance of a magnetic field outside the solenoid", and that for non-uniform rotation "an electric field arises". — [IOP abstract](https://iopscience.iop.org/article/10.1088/0305-4470/25/18/020)
- Other AB-on-toroid theory: Lyuboshitz & Lyuboshitz, JETP 91, 673–681 (2000), on the momentum-transfer cross-section for AB scattering on a toroidal solenoid; Afanasiev, Physica Scripta 48, 385–392 (1993), on the toroidal Aharonov–Casher effect. Crossref metadata only. — [doi:10.1134/1.1326959](https://doi.org/10.1134/1.1326959); [doi:10.1088/0031-8949/48/4/001](https://doi.org/10.1088/0031-8949/48/4/001)
- The field conjugate to the magnetic toroidization is ∇×B. The interaction energy is ∫(∇×B)·T_M d³r, and the electric-toroidization energy vanishes because ∇×E = 0 in statics. — [Planes, Castán & Saxena, arXiv:1504.04958](https://arxiv.org/abs/1504.04958)
- The Lockheed patent (section 5) itself describes the AB effect as shifting phase "without exchanging energy". — [US 9,502,202](https://patents.google.com/patent/US9502202)

### Inferences
- From the script:
  - h/2e = 2.067833848×10⁻¹⁵ Wb.
  - The AB phase e(h/2e)/ħ is exactly π, so even multiples of h/2e give 0 and odd multiples give π. This is exactly Osakabe's "0 or π".
  - The Dirac flux h/e is 2 × h/2e, giving an AB phase of 2π, which no interference experiment can see.
- **Energy and force:**
  - A static toroid has no external E or B, so a charge outside feels no Lorentz force and gains no energy.
  - A static toroidal moment couples only to ∇×B, which in vacuum outside sources is μ₀ε₀∂E/∂t, zero in statics. So a static toroid neither feels nor exerts forces from or on external static fields.
  - What remains is (i) the AB phase, which is quantum, topological, gauge-invariant and not an energy, and (ii) radiation, if the toroid is driven in time (section 1).

  The talk's "no fields but strong vector potential" is **mainstream** as a description of an ideal static toroid. Any reading that it carries usable energy or force is **refuted** by the Lorentz force law and by the ∇×B coupling.
- Afanasiev & Dubovik 1992 show that fields appear outside only when the toroid moves (rotates). That is consistent with the above: the outside is field-free only for the static, ideal case.

### Gaps
- Not opened: the Tonomura PRL full text (copper-layer and niobium thicknesses, flux values); the full text of Afanasiev & Dubovik 1992; Feynman's static charge-plus-magnet angular-momentum discussion (the Feynman Lectures site returned 403).

## 3. Supertoroidal pulses, hypertoroidal moments, "fractal" toroids and skyrmionic textures: what is real

### Takeaway
Three real, separate things share the vocabulary:
- **Supertoroidal light pulses**: free-space solutions of Maxwell's equations (Shen, Hou, Papasimakis & Zheludev 2021). They have skyrmionic field textures and self-similar, "fractal-like" arrangements of singularities as the order parameter α increases.
- **Hypertoroidal moments**: a defined hierarchy in ferroelectrics and ferrotoroidics (Prosandeev & Bellaiche 2009; Thorner et al. PRB 2014), built as curls of curls of the toroidization.
- **Higher-order "supertoroidal" currents in matter**: iterated (fractal) windings of solenoidal currents (Afanasiev, Dubovik, Stepanovsky and others).

None of these is a monopole. None of them confines vacuum energy or links to matter disassembly, and none is described as forming in charge-separated media through hydrodynamic shear. Those parts are the speaker's extrapolation.

### Cited Findings
- Shen, Hou, Papasimakis & Zheludev, "Supertoroidal light pulses as electromagnetic skyrmions propagating in free space", Nat. Commun. 12, 5891 (2021). Crossref-checked.
  - They extend the Hellwarth–Nouchi pulse (PRE 54, 889, 1996). The generating function is f(r,t) = f₀/[(q₁ + iτ)(s + q₂)^α], with s = r²/(q₁ + iτ) − iσ, τ = z − ct, σ = z + ct, and finite energy requires α ≥ 1.
  - The pulses show "skyrmionic structure… multiple singularities in the Poynting vector maps and fractal-like distributions of energy backflow".
  - Their singular shells form "a self-similar pattern that seems infinitely repeated".

  — [PMC8501108](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8501108/); [doi:10.1038/s41467-021-26037-w](https://doi.org/10.1038/s41467-021-26037-w)
- The same paper: "higher order, supertoroidal modes in matter have been introduced in the form of the fractal iterations of solenoidal currents" (its refs. 42–47). Those references are:
  - Afanasiev & Stepanovsky, J. Phys. A 28, 4565 (1995)
  - Afanasiev & Dubovik, Phys. Part. Nucl. 29 (1998)
  - Afanasiev, J. Phys. D 34, 539 (2001)
  - Boardman et al. (SPIE 2005)
  - Gurvitz et al., Laser Photon. Rev. 13, 1800266 (2019)
  - Nemkov et al., PRA 98, 023858 (2018)

  — [Shen et al. 2021](https://doi.org/10.1038/s41467-021-26037-w)
- A follow-up (Shen, Papasimakis & Zheludev, arXiv:2204.05676) introduces non-diffracting supertoroidal pulses. Per its abstract, "the topological properties of the original pulses do not survive upon propagation" (snippet, not opened). — [arXiv:2204.05676](https://arxiv.org/abs/2204.05676)
- "Hypertoroidal moment in complex dipolar structures": Prosandeev & Bellaiche, J. Mater. Sci. 44, 5235–5248 (2009). Crossref only. — [doi:10.1007/s10853-009-3460-5](https://doi.org/10.1007/s10853-009-3460-5)
- "Axial hypertoroidal moment in a ferroelectric nanotorus: A way to switch local polarization": Thorner, Kiat, Bogicevic & Kornev, Phys. Rev. B 89, 220103 (2014). Crossref only. — [doi:10.1103/PhysRevB.89.220103](https://doi.org/10.1103/PhysRevB.89.220103)
- How hypertoroidal moments are defined:
  - Spiral configurations of toroidal moments are "characterized by divergence-free vectors expressed as the curl of higher order toroidal moments. These moments are denoted as hypertoroidal moments", with T_M → T_M + ∇×T_M⁽²⁾ and so on, "generalized to any order".
  - A transverse hypertoroidal moment in off-centre nanorings "is a polar vector".

  — [Planes, Castán & Saxena, arXiv:1504.04958](https://arxiv.org/abs/1504.04958)
- The "synthetic skyrmion in a BEC at Aalto" is Lee, Gheorghe, Tiurev, Ollikainen, Möttönen & Hall, "Synthetic electromagnetic knot in a three-dimensional skyrmion", Sci. Adv. 4, eaao3820 (2018). The work is from Amherst College together with Aalto University (QCD Labs) and Jyväskylä.
  - The texture was made "within an optically trapped and spin-polarized ferromagnetic ⁸⁷Rb BEC".
  - Its synthetic field lines form a Hopf fibration ("a Shankar skyrmion with a unit Hopf number").
  - The fields are synthetic gauge fields acting on neutral atoms, not real E or B.

  — [PMC5834309](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5834309/); [doi:10.1126/sciadv.aao3820](https://doi.org/10.1126/sciadv.aao3820)

### Inferences
- "Supertoroid (Russia) / hypertoroid (West)":
  - "Supertoroidal" is used by the Southampton group (UK) for both pulses and higher-order modes in matter.
  - "Hypertoroidal" is used in the ferroelectrics literature (US and France).
  - The Russian line (Afanasiev, Dubovik, Stepanovsky at JINR Dubna) is cited for the fractal iterations of solenoidal currents.

  So the iterated-toroid concept is **mainstream** as mathematics and electrodynamics, but the East/West naming split is the speaker's simplification.
- In none of these sources does a higher order make a structure "non-radiating" or stable. Nemkov et al. say the higher-order (MSR) terms radiate with the same pattern as their parents. A static ideal toroid of order 1 already has no outside field.
- Fractal nesting in the speaker's sense (golden-ratio scaling, at most 48 sublevels, swastika shapes, hollow spheres that collapse matter) appears in none of the opened sources: **unverified**.
- The Aalto/Amherst image shows a synthetic gauge-field texture in a cold-atom condensate. Comparing it with SEM images of discharge debris is the speaker's analogy, not a result of the paper.

### Gaps
- The PRB paper on hypertoroidal moments "published January 20", in 2023 or 2024, which the speaker says he predicted in September 2022, was not found. The only PRB hit with "hypertoroidal" in the title is Thorner et al. 2014. The candidate needs the slide or title.
- The full text of Prosandeev & Bellaiche was not opened. Its exact definition (my reading is a hypertoroidal moment h = ½Σ r × g of local toroidal moments g) is unconfirmed.

## 4. Afanasiev & Dubovik, and the "Physics of Plasmas paper with 1.36 kT replicating 1976 work in Japan"

### Takeaway
The Afanasiev & Dubovik work is real: J. Phys. A 25, 4869 (1992) on the toroidal solenoid, and Phys. Part. Nucl. 29 (1998) on remarkable charge–current configurations. The plasma paper is, per the user's audit and Crossref, Yue et al., Phys. Plasmas 28, 042303 (2021). It is a particle-in-cell (PIC) **simulation** in which a moving electron vortex forms an expanding magnetic ring "with an internal magnetic field up to 1000 Tesla". The abstract mentions no 1.36 kT, no 1976 Japanese experiment and no field converging to a point singularity.

### Cited Findings
- Yue, Chen, Geng, Yuan, Weng, Bulanov, … Mima, Sheng & Zhang, "Dynamics of moving electron vortices and magnetic ring in laser plasma interaction", Phys. Plasmas 28, 042303 (1 Apr 2021). Crossref-checked. — [doi:10.1063/5.0034098](https://doi.org/10.1063/5.0034098)
- From that abstract:
  - "In three dimensions, the structure appears as an expanding magnetic ring with an internal magnetic field up to 1000 Tesla."
  - "The moving velocity can be as large as a 0.2 c level".
  - More than 5 % of the laser energy is carried "as far as 50 μm away".
  - The method is "multi-dimensional Particle-in-Cell simulations".

  — [Strathprints 77046](https://strathprints.strath.ac.uk/77046)
- Afanasiev & Dubovik 1992, J. Phys. A 25, 4869 (see section 2). — [IOP](https://iopscience.iop.org/article/10.1088/0305-4470/25/18/020)
- For "Some remarkable charge–current configurations", Phys. Part. Nucl. 29 (1998), the sources disagree:
  - Crossref lists a single author (Afanasiev) and page 366.
  - Shen et al. 2021 cite it as Afanasiev & Dubovik, pp. 891–945.

  — [Crossref doi:10.1134/1.953084](https://doi.org/10.1134/1.953084); [Shen et al. ref. 43](https://doi.org/10.1038/s41467-021-26037-w)
- Experimental kT laser-driven fields exist in other work. Ehret et al. (Phys. Rev. E) report kilotesla-scale quasi-static magnetised plasma structures (snippet, not opened). — [MEPhI repository](https://openrepository.mephi.ru/handle/123456789/27655)

### Inferences
- If a 1.36 kT number appears in Yue et al., it is a colour-bar peak or a 2-D value from a simulation, not a measurement. The abstract's 3-D figure is ~1000 T. **Unverified**: the full text was not opened.
- The talk's question, "what is the field at the point singularity?", has no basis in the abstract. A PIC simulation resolves the ring core on its grid; no divergence to a point is reported.
- "Replicating 1976 work in Japan" is **unverified**. Laser-produced self-generated magnetic fields have a 1970s history, but no source linking Yue 2021 to a 1976 Japanese experiment was found.

### Gaps
- The Yue et al. full text (figures, 1.36 kT, 2-D versus 3-D peaks, laser parameters) is paywalled and was not opened.
- The "Moscow Nuclear Physics Institute" paper on the "superconductivity of the Dirac monopole" (the RCCNT&BL 2008 conference) and "Vishnevsky 2008" were not found:
  - RCCNT&BL-15 was held in Sochi in October 2008 ([New Energy Times](https://newenergytimes.com/v2/conferences/2008/RCCNTBL15/RCCNTBL-15-Announce.shtml)), but no proceedings were found.
  - The closest published ideas are Birbrair's superconducting ball-lightning preprint ([arXiv:physics/0102012](https://arxiv.org/abs/physics/0102012)) and Stephan's 2024 monopole–ball-lightning preprint ([arXiv:2411.00240](https://arxiv.org/abs/2411.00240)), both search-level only.

## 5. The patents: Lockheed Martin (coherent matter-wave beams), the "1988 Max Planck" patent, and the 1993 potentials-only communication patent

### Takeaway
- **Lockheed Martin:** US 9,502,202 B2 (Arman & Chase; filed 28 Dec 2011; granted 22 Nov 2016). It claims phase synchronisation of electron streams through the Aharonov–Bohm effect, modelled with the Kuramoto equations. It does not mention toroids or monopoles.
- **"1988 Max Planck" patent:** no such patent was found. The Max Planck patent in this field is US 6,476,383 (Esslinger, Hänsch & Bloch; priority 1999), a magnetic trap with continuous extraction of a matter wave, i.e. the atom laser. Lockheed cites it. A 1988 matter-wave patent in the same citation list is MIT's US 4,886,964.
- **1993 patent:** US 5,845,220 by Harold E. Puthoff (filed 23 Aug 1993, priority 31 May 1991, granted 1 Dec 1998). It claims communication by scalar and vector potentials "without including any electromagnetic field", using solenoids or toroidal coils plus driven plates, received by a Josephson junction. It does not use the words anapole or toroidal moment.

### Cited Findings
- **US 9,502,202 B2 (Lockheed Martin):**
  - Inventors Moe J. Arman and Charles Chase; application US 13/338,925, filed 28 Dec 2011; published as US 2013/0169157 on 4 Jul 2013; granted 22 Nov 2016.
  - Each beam-generating unit has a cathode–anode cavity producing a stream of charged particles (electrons as the example). A magnetic field orthogonal to the channels phase-synchronises the streams "through a vector potential", invoking the AB effect "without exchanging energy" and the Kuramoto model.
  - Numbers in the embodiment: ~100 gauss (1 T in the synchronisation calculation); a housing of about 10 µm × 1 µm × 0.1 µm; a 0.1 mm gap; electron mean free path ~9.3 mm; collision time ~100 ns; electron speed ~100 km/s; synchronisation in "a few hundred nanoseconds"; "at least one billion" particles.
  - No toroid or monopole is mentioned. Cited prior art includes US 6,476,383 B1 (Max Planck, priority 1999) and US 4,886,964 (MIT, 1988).

  — [Google Patents US9502202](https://patents.google.com/patent/US9502202)
- **US 6,476,383 (Max Planck Society):** "Device and method for generating and manipulating coherent matter waves", by Esslinger, Hänsch & Bloch. Filed 31 Aug 1999; granted 5 Nov 2002. It describes a trap of two quadrupole coils plus an Ioffe coil, a shield, and continuous extraction of matter waves (snippet, not opened). — [Google Patents US6476383](https://patents.google.com/patent/US6476383)
- **US 5,845,220 (Puthoff):**
  - Assignee EarthTech International. Filed 23 Aug 1993 as a continuation of an application of 31 May 1991; granted 1 Dec 1998.
  - Abstract: information is sent "by transmitting a signal comprising scalar and vector potentials without including ay [sic] electromagnetic field".
  - Geometries: a solenoid plus plates driven by a voltage that tracks dI/dt (FIGS. 2–3), and two coaxial toroidal coils plus voltage-driven screens (FIG. 7).
  - Receiver: a DC-biased Josephson junction.
  - It distinguishes itself from Gelinas because "the present invention employs a suppressed electromagnetic field (i.e. E=B=0)".
  - It quotes a Josephson frequency of 483.2 MHz/µV.

  — [Google Patents US5845220](https://patents.google.com/patent/US5845220)
- Later continuation filings by Puthoff and Eusebi (Quantcomm LLC), such as US 10,361,792 (filed 2016, granted 2019), describe "a communication system using vector and scalar potential" for sea water and plasma (snippet, not opened). — [Justia US 2018/0062765](https://patents.justia.com/patent/20180062765)

### Inferences
- **Lockheed:**
  - The patent's own numbers are consistent with each other: 9.3 mm / 100 ns = 93 km/s, against the stated ~100 km/s (script).
  - Physics check (my inference): an AB phase is the same for every electron that takes the same path, so it cannot pull independent electrons' random phases together.
  - By Liouville's theorem, a static magnetic field cannot raise a beam's phase-space density, and so cannot raise its degeneracy or coherence.
  - Kuramoto synchronisation needs a coupling between the oscillators. A shared vector potential is not such a coupling.

  Status: **unverified**. It is a granted patent, and I found no published demonstration. The talk's link between this patent and monopoles is not in the patent.
- **"1988 Max Planck":** probably US 6,476,383 (the 1999 atom laser) with the date misremembered, or a mix-up with MIT's 1988 US 4,886,964. Which patent's front page the slide shows is **unverified**.
- **Puthoff:** the talk's "1993" matches the filing date. "Nine months after he came back from Moscow" and "working with NASA" are **unverified**; no biography source was opened.
- The physics test for the Puthoff patent (my inference from Faraday's law):
  - If E = B = 0 throughout a simply connected region, the potentials there are a pure gauge and carry no gauge-invariant signal.
  - On any loop that threads a toroid with changing flux, ∮E·dl = −dΦ/dt ≠ 0, so E cannot vanish everywhere on that loop.

  So a time-varying signal "without any field" can be detected only where E ≠ 0 somewhere along the detector's circuit. Status: **contested**. The patent is granted, and no independent peer-reviewed demonstration was found.
- The patent's 483.2 MHz/µV is 0.08 % below 2e/h = 483.598 MHz/µV (script). It is a small numerical slip in the patent text, as reported by the page fetch.

### Gaps
- A 1988 Max Planck patent was not found in patent searches.
- The claims of US 6,476,383 were not opened.
- Whether any Puthoff–NASA work exists was not checked.

## 6. "Electromagnetic phantom", closed Poynting vectors and "vacuum currents" (Fryberger, Zvirblis)

### Takeaway
"Freiberger" is David Fryberger of SLAC. He wrote a vorton model of ball lightning (SLAC-PUB-6473, 1994), a small-particle model of the "cavity lights" seen in superconducting niobium cavities (SLAC-PUB-8405; NIM A 459, 29, 2001), and a ball-lightning model for the cavity lights (SLAC-PUB-13583, 2009) that uses symmetrised Maxwell equations with "vacuum currents". This is single-author speculative work with no independent confirmation (**unverified**). "Zvirblis 1995" (electromagnetic phantoms, closed Poynting vectors) was not found in any indexed source. Closed circulating Poynting flux around static charge-plus-magnet systems is standard electrodynamics, and it implies no energy capture.

### Cited Findings
- SLAC-PUB-8405, "A Small particle model as a possible explanation of recently reported cavity lights", NIM A 459, 29–43 (2001). The abstract describes luminous objects "inside superconducting niobium cavities under an electrical excitation equivalent to an average acceleration of ∼4 MV/m", some "executing elliptical orbits". It concludes that "it is highly unlikely that a small object orbiting in such a potential well can be made of any known material". — [INSPIRE](https://inspirehep.net/search?p=a%20Fryberger%20and%20t%20cavity%20lights)
- SLAC-PUB-13583 (2009), "A Ball Lightning Model as a Possible Explanation of Recently Reported Cavity Lights", and SLAC-PUB-6473 (1994), "A Model for ball lightning". INSPIRE metadata only. — [INSPIRE](https://inspirehep.net/search?p=a%20Fryberger%20and%20t%20ball%20lightning)
- The OSTI record of the 2009 report says the ball-lightning model uses symmetrised Maxwell equations with "vacuum currents" (snippet; OSTI refused the connection). — [OSTI 962083](https://www.osti.gov/biblio/962083)
- Fryberger's later arXiv work includes "On Magneticons and some related matters" (2018) and a proposal to search for magnetic charge at Belle II (with Sullivan, 2017). arXiv metadata only. — [arXiv:1807.09606](https://arxiv.org/abs/1807.09606)
- Supertoroidal pulses have "multiple singularities in the Poynting vector maps" and energy backflow. This is mainstream, with no energy capture involved. — [Shen et al. 2021](https://doi.org/10.1038/s41467-021-26037-w)
- In the anapole, "the energy is concentrated inside the particle" while far-field scattering vanishes. That is a driven, lossy resonance, not a source of energy. — [Miroshnichenko et al. 2015](https://doi.org/10.1038/ncomms9069)

### Inferences
- A closed (circulating) Poynting vector with zero divergence, as in a static charge next to a magnet, carries field momentum but no net energy flow (∮S·dA = 0). Collapsing one closed loop to "energise the other" is not a mechanism that appears in any opened source: **unverified**.
- That vacuum currents make baryons "lose their reference to the Dirac sea" is Fryberger's speculative vorton framework. Status: **unverified** (no independent tests found).

### Gaps
- Zvirblis (Zhvirblis) 1995 was not found in Crossref, arXiv or web search; the only hit was a Substack blog. It may be the popular-science magazine Khimiya i Zhizn, which is unconfirmed.
- The year and laboratory of the cavity-light observations were not in the abstracts opened.

## 7. fMRI, the brain's magnetic fields, and "the brain produces toroidal moments"

### Takeaway
fMRI's BOLD signal comes from blood oxygenation: paramagnetic deoxyhaemoglobin acts as an endogenous contrast agent. It does not image the brain's currents or magnetic fields. The brain's own fields, measured by MEG, are 10 fT to 1 pT, and typical source strengths are about 10 nAm. Neural current loops can have a toroidal component in principle, but a static toroidal moment has no outside field and the radiated power at brain frequencies is negligible (~10⁻⁶⁵ W). The claim that fMRI "proves the brain produces toroidal moments" is **refuted** as a reading of what fMRI measures. The added claims about imprinting information into relic-neutrino condensates are **unverified**, with no source.

### Cited Findings
- Ogawa, Lee, Kay & Tank, PNAS 87, 9868–9872 (1990): "Paramagnetic deoxyhemoglobin in venous blood is a naturally occurring contrast agent for magnetic resonance imaging". They defined "blood oxygenation level-dependent (BOLD) contrast". — [PMC55275](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC55275/); [doi:10.1073/pnas.87.24.9868](https://doi.org/10.1073/pnas.87.24.9868)
- Pauling & Coryell, "The Magnetic Properties and Structure of Hemoglobin, Oxyhemoglobin and Carbonmonoxyhemoglobin", PNAS 22, 210–216 (1936). This is the source of oxy-diamagnetism and deoxy-paramagnetism; Crossref metadata only, since PMC holds only a scan. — [doi:10.1073/pnas.22.4.210](https://doi.org/10.1073/pnas.22.4.210)
- Hämäläinen, Hari, Ilmoniemi, Knuutila & Lounasmaa, RMP 65, 413–497 (1993):
  - "the weak 10 fT–1 pT magnetic fields produced by electric currents flowing in neurons"
  - "Neuromagnetic signals are typically 50–500 fT"
  - "A typical strength of a dipole… is 10 nA m"
  - a single postsynaptic potential ≈ 20 fA m
  - "a dipole moment of 10 nAm would correspond to 40 mm² of active cortex" at 100–250 nA/mm²

  — [Aalto PDF](https://acris.aalto.fi/ws/portalfiles/portal/14383291/RevModPhys.65.413.pdf); [doi:10.1103/RevModPhys.65.413](https://doi.org/10.1103/RevModPhys.65.413)

### Inferences
- From the script:
  - The radial field of a 10 nAm tangential dipole 30 mm deep in a 100 mm sphere peaks at **355 fT**, inside the published 50–500 fT band. Volume currents do not change the radial component in a spherical conductor.
  - 10 nAm / (250 nA/mm² × 1 mm) = 40 mm², as published.
- An illustrative cortical "toroid" (not from any paper): 10 nAm loops closing on a 5 mm radius give m = 2.5×10⁻¹¹ A m². Arranged on a 20 mm ring they give T = 2.5×10⁻¹³ A m³, zero static field outside, and 5×10⁻⁶⁵ W radiated at 10 Hz.
- BOLD phase and susceptibility maps show where blood oxygenation changes. "Phase singularities" in such maps are features of the MR phase image, not of neural currents.

### Gaps
- Which fMRI figure the slide shows was not identified.
- The "magnetic loop known since the '60s" was not identified. It may refer to the cardiac or neural magnetic fields first measured in the 1960s.

## 8. "Kilowatts over an 8-micron wire" (Moray King 2012; Avramenko single-wire)

### Takeaway
Single-wire power transmission (the Avramenko plug; Strebkov) is a real but marginal line of work. It is an open-circuit, high-voltage, resonant arrangement in which the return current is capacitive (displacement current). Ohmic heating depends only on the current. At kilovolts and above, kilowatts in an 8 µm copper wire are ordinary physics. At mains voltage they would vaporise the wire. The claim is therefore not extraordinary in itself, and it is extraordinary only if it is made at low voltage with no heating. The voltage and the source were not found, so the status is **unverified**.

### Cited Findings
- Avramenko & Avramenko, US 6,104,107, "Method and apparatus for single line electrical transmission" (filed 11 Jan 1995; granted 15 Aug 2000). Snippet; patent page not opened. — [Google Patents US6104107](https://patents.google.com/patent/US6104107)
- A 2023 Energies review reportedly says that in July 1990 at the Moscow Power Engineering Institute, 1.3 kW was transmitted through a single tungsten wire 2.75 m long (snippet; the MDPI page returned 403, so the diameter is unknown). — [MDPI Energies 16, 5089](https://www.mdpi.com/1996-1073/16/13/5089)
- IACS: 100 % conductivity is 58 MS/m, i.e. 1.7241×10⁻⁸ Ω·m at 20 °C (snippet). — [Wikipedia: IACS](https://en.wikipedia.org/wiki/International_Annealed_Copper_Standard)

### Inferences
- From the script, for 8 µm copper (343 Ω/m), carrying 1 kW:

| Voltage | Current | Current density | Dissipation |
|---|---|---|---|
| 230 V | 4.35 A | 8.7×10¹⁰ A/m² | 6.5 kW per metre: the wire vaporises |
| 10 kV | 0.1 A | 2×10⁹ A/m² | 3.4 W/m: hot but plausible for a short run |
| 100 kV | 0.01 A | — | 0.034 W/m: trivial |

- The "mainstream" reading is a high-voltage, low-current resonant (Tesla-style) line, in which the wire carries the charging current of the far-end capacitance.

### Gaps
- Moray King's 2012 presentation (wire material, voltage, frequency, measured power) was not found.
- The diameter in Avramenko's experiment is not confirmed: 8 µm is the speaker's figure, and the review page was not opened.

## 9. Other references in the talk

### Takeaway
"The UA model", "the 1995 classified-then-declassified work with six yin-yang pairs" and "Holt 1979" were not researched in depth or not found. They are listed for the report writer as open items.

### Cited Findings
- The SLAC tribute says Fryberger died in May 2026 (snippet). — [SLAC AHRO](https://ahro.slac.stanford.edu/news/2026-06-03-memory-david-fryberger)

### Inferences
- None beyond the labels below.

### Gaps
- "UA paper/model" (transcript line ~2001): unidentified.
- "Vishnevsky 2008 flux model": unidentified (see section 4).
- The 1995 declassified work: unidentified.
- A.C. Holt 1979, "Field resonance propulsion concept" (NASA): not opened, and outside this scope.

## 10. The user's audit: its formulas checked

Source: `research/external/2026-10-26-greenyer-audit-and-crate.md`, an unverified user document. Each formula is checked against the opened sources and the check script.

### Takeaway
Correct as written:
- the toroidal-dipole formula and its implementation
- the Aharonov–Bohm phase eΦ/ħ

Need correction:
- **The anapole condition P + ikT = 0.** It holds only in Gaussian units, where T carries 1/c. With the audit's SI T it must read **p + (ik/c)T = 0** (for e^(−iωt); the sign of i flips for e^(+iωt)).
- **The Dirac-monopole field.** It mixes units.

Not found in the papers they cite (invented):
- the "Yue2021" and "Vishnevsky2008" field models
- the 1,360 T constant
- the golden-ratio and n ≥ 3 rules for the fractal toroidal moment

### Cited Findings
- **Toroidal dipole T⁽¹⁾ = (1/10)∫[r(r·J) − 2r²J] d³r**, implemented as 0.1·Σ[r(r·Idl) − 2r²Idl]. Correct. The script's version of this formula gives |T| = N I V/(4π) for a toroidal solenoid, the SI form of Nemkov's IV/(4πc). It also matches ½Σ r × m. — [Nemkov et al. 2018](https://arxiv.org/abs/1804.07342)
- **Anapole "P + ikT⁽¹⁾ = 0":** Nemkov et al. write d + ikτ = 0 in units where τ carries c⁻¹. — [Nemkov et al. 2018](https://arxiv.org/abs/1804.07342)
- **"Closed internal soliton loop" ∮S·dA = 0:** zero net radiated flux is what an anapole means. The near field is non-zero, with energy "concentrated inside the particle". No source calls it a soliton. — [Miroshnichenko et al. 2015](https://doi.org/10.1038/ncomms9069)
- **"Nested recursively to order n ≥ 3 the structure forms a non-radiating boundary":** in no source. An ideal static toroid of order 1 already has no external field ([Miroshnichenko 2015](https://doi.org/10.1038/ncomms9069)), and higher orders radiate with their parents' patterns ([Nemkov 2018](https://arxiv.org/abs/1804.07342)).
- **Recursion T⁽ᵏ⁾ = ½Σ r_k,i × T⁽ᵏ⁻¹⁾_i with R_k = R₀φ⁻ᵏ and N_k ∈ [2, 48]:** the ½Σ r × (lower moment) form matches the hypertoroidal construction in the ferroelectric literature. The review describes hypertoroidal moments as curls of higher-order toroidizations, polar and axial in turn ([Planes et al.](https://arxiv.org/abs/1504.04958); [Prosandeev & Bellaiche 2009](https://doi.org/10.1007/s10853-009-3460-5)). The golden-ratio radii and the 2–48 sub-tori come from the talk, not from any source.
- **AB phase Δφ = (e/ħ)∮A·dl = (e/ħ)Φ:** correct. The script gives e(h/2e)/ħ = π, matching Osakabe's "0 or π". — [Osakabe et al. 1986](https://doi.org/10.1103/PhysRevA.34.815)
- **Yue et al. 2021:** the paper exists and is Crossref-checked. The abstract gives "up to 1000 Tesla", from a simulation. It does not give 1,360 T or a "central dead-stop phase singularity". — [Strathprints](https://strathprints.strath.ac.uk/77046)
- **"Zel'dovich (1957)":** received in 1957, published 1958. Acceptable as "1957/1958". — [Yang & Bozhevolnyi](https://arxiv.org/abs/1810.01221)
- **"Afanasev & Dubovik (1992/1998)":** the 1992 paper is confirmed. The 1998 authorship is inconsistent between Crossref and Shen et al. — [IOP 1992](https://iopscience.iop.org/article/10.1088/0305-4470/25/18/020); [Crossref 1998](https://doi.org/10.1134/1.953084)

### Inferences
- **Fractal recursion:**
  - Its parity alternates (r × axial = polar, r × polar = axial), so the "order" changes the moment's symmetry class. That deserves a unit test.
  - The code's magnitude ½·N·r_k·T⁽ᵏ⁻¹⁾ is valid only if every r_k,i ⊥ T⁽ᵏ⁻¹⁾_i and all the cross products are parallel, i.e. an ideal symmetric ring.
  - The result depends on the choice of origin unless the lower-order net moments vanish.
- **Dirac monopole:**
  - The audit takes g_D = 2πħ/e = h/e (webers, a flux) but then writes B = (μ₀g/4π) r/r³. That needs g as a pole strength in A·m, equal to h/(μ₀e).
  - With g in webers, B = g r/(4πr³).
  - The Dirac string vector potential ∝ (1 − cos θ)/(r sin θ) φ̂ is the standard textbook form, not checked against an opened source.
  - The script checks h/e = 2 × (h/2e).
- **"Yue2021" field model:**
  - B₀ exp(−((ρ−R₀)² + z²)/a²) φ̂ + B₀ (R₀/√(ρ² + z² + ε²))^γ ẑ is **invented**: no such formula is in the abstract, and the full text was not opened.
  - The test's B_core > 10¹⁰ T is set by the arbitrary ε and γ, not by physics.
  - A central axial field ∝ r^(−γ) also violates ∇·B = 0 unless it is paired with radial components.
- **"Vishnevsky2008" field model** (monopole plus Beltrami helix ∇×B = αB): the source was not found, so the model is **unverified**. A Beltrami field is a standard force-free field. A monopole core contradicts ∇·B = 0 in classical electromagnetism unless magnetic charge is assumed.
- **Kelvin magnetophoresis** F = (χV/μ₀)(B·∇)B = (χV/2μ₀)∇B²: the standard form for small χ in a current-free region. It is outside this scope and not checked against a source.

### Gaps
- Yue et al. full text (for 1.36 kT).
- Vishnevsky 2008.
- Zvirblis 1995.

## 11. Calculators for the crates

How toroidal multipoles fit `geopotential::Field`:
- `geopotential::Field` (`crates/geopotential/src/field.rs`) is a static field with a scalar potential, evaluated outside its sources. It unifies gravity (`Coefficients`), the geomagnetic main field (`magnetics`) and the gravitomagnetic dipole (`gravitics`).
- A static toroidal multipole produces **no field and no scalar potential outside its sources**, so its `Field` contribution is identically zero. Gauss coefficients cannot represent it, and it should not implement `Field`.
- Toroidal moments belong as **source-moment calculators** (from currents), as **dynamic radiation** through p_eff = p + (ik/c)T, and in the **coupling to ∇×B**.
- Suggested home: a new `magnetics::toroidal` module, next to `magnetics::magnets` (loop and solenoid fields).
- The electromagnetic AB phase fits `magnetics`. `gravimetry::atom` already has the gravitational Aharonov–Bohm rate.
- Constants belong in `metrology`.

| Function signature | Formula | Primary source | Published test value | Target crate::module | Status |
|---|---|---|---|---|---|
| `toroidal_dipole(segs: &[([f64;3],[f64;3])]) -> [f64;3]` (position, I·dl) | T = (1/10)Σ I[(r·dl)r − 2r²dl] (SI) | Dubovik & Tugushev 1990; Nemkov et al. 2018 | torus: \|T\| = N I V/(4π) (script: 1e-3) | magnetics::toroidal | mainstream |
| `toroidal_solenoid_moment(n: f64, i: f64, major: f64, minor: f64) -> f64` | T = N I π R a²/2 = N I V/(4π) | Nemkov et al. 2018 App. B (IV/4πc, Gaussian) | identical to `toroidal_dipole` on a discretised torus | magnetics::toroidal | mainstream |
| `effective_dipole(p: Complex, t: Complex, k: f64) -> Complex` | p_eff = p + (ik/c)T, e^(−iωt) | Nemkov 2018 (d + ikτ = 0); Yang & Bozhevolnyi eq. 27 | p = −(ik/c)T leaves < 1e-3 of the power (script 8e-8) | magnetics::toroidal | mainstream |
| `toroidal_radiated_power(t: f64, omega: f64) -> f64` | P = μ₀ω⁶T²/(12πc⁵) | derived; pattern "scaled by ω²" (Miroshnichenko 2015) | numerical far-field power within 1e-2; P(2ω)/P(ω) = 64 | magnetics::toroidal | mainstream |
| `mie_a1_zero_diameter(n: f64, wavelength: f64) -> f64` | root of m ψ₁(mx)ψ₁′(x) − ψ₁(x)ψ₁′(mx), x = πd/λ | Miroshnichenko 2015 | n = 4, 550 nm → ∼204 nm (script 204.0, ±3) | magnetics::toroidal | mainstream |
| `hypertoroidal_moment(centres: &[[f64;3]], moments: &[[f64;3]]) -> [f64;3]` | ½Σ rᵢ × Tᵢ (parity alternates per order) | Prosandeev & Bellaiche 2009; Planes et al. 2015 | none opened; unit-test the ring identity against `toroidal_dipole` | magnetics::toroidal | mainstream (definition); golden-ratio radii unverified |
| `aharonov_bohm_phase(charge: f64, flux: f64) -> f64` | Δφ = qΦ/ħ | Aharonov & Bohm 1959; Osakabe et al. 1986 | Φ = h/2e → π ("0 or π") | magnetics (beside gravimetry::atom) | mainstream |
| `FLUX_QUANTUM`, `JOSEPHSON_CONSTANT`, `DIRAC_FLUX` | h/2e; 2e/h; h/e = 2Φ₀ | CODATA; Tonomura 1986 | 2.067833848e-15 Wb; 483.598 MHz/µV (patent: 483.2) | metrology | mainstream |
| `loop_emf(dflux_dt: f64) -> f64` | ∮E·dl = −dΦ/dt (judges "E = B = 0 signalling") | Faraday's law; used against US 5,845,220 | none (identity) | magnetics | mainstream; patent claim contested |
| `supertoroidal_generating_fn(r, z, t, q1, q2, alpha) -> Complex` | f = f₀/[(q₁ + iτ)(s + q₂)^α], s = r²/(q₁ + iτ) − iσ | Shen et al. 2021; Hellwarth & Nouchi 1996 | α = 1 → Hellwarth–Nouchi; test (∇² − ∂²/c²∂t²)f = 0 numerically | magnetics::toroidal (or web/ as a 3-D model) | mainstream |
| `meg_radial_field(q: [f64;3], r0: [f64;3], r: [f64;3]) -> f64` | B_r = μ₀/(4π)(Q × (r − r₀))·r̂/\|r − r₀\|³ (sphere: exact radial) | Hämäläinen et al. 1993 | 10 nAm, 30 mm deep, 100 mm sphere → 355 fT, inside "50–500 fT" | magnetics | mainstream |
| `wire_joule_per_length(p: f64, v: f64, d: f64, rho: f64) -> f64` | I = P/V; P′ = I²ρ/(πd²/4) | IACS 1/58 MS/m | ρ = 1.7241e-8 Ω·m; 8 µm Cu = 343 Ω/m | magnetics::magnets | mainstream; King claim unverified |
| `cs_anapole_kappa() -> (f64, f64)` (reference value) | κ = 0.127 ± 0.019 | Wood et al. 1997, via Haxton 1997 | ~7σ (script 6.7σ) | physica (claims/reference table) | mainstream (measured) |
| Yue2021 / Vishnevsky2008 / golden-ratio FTM models (audit) | as written in the audit | none found | none | physica (labelled claim models, examples only) | unverified (Yue formula not in the paper; Vishnevsky not found) |
