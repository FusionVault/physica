# Toroidal moments meet magnetic monopoles and duality: what is known, what is new, what could be proved

Scope: whether anything genuinely new and provable lies where toroidal moments, magnetic charge and electromagnetic duality meet. Four questions: (1) the multipole family with electric and magnetic sources, and a duality-covariant anapole; (2) Dirac strings as thin toroidal flux lines, and the AB phases of composite particles; (3) toroidal moments carried by spin-ice monopole currents; (4) a toroidal analogue of the Aharonov–Casher effect. Builds on `crates/magnetics/src/{toroidal,monopole}.rs` and `research/topics/topological-monopoles-and-transmutation/notes/{monopoles,toroidal_moments}.md`, which are not redone here.

Conventions: SI, time dependence e^(−iωt). Magnetic charge is in A·m, so ∇·B = µ₀ρ_m and −∇×E = ∂B/∂t + µ₀J_m. The toroidal dipole is the crate's T = (1/10)∫[(r·J)r − 2r²J] dV, with anapole p = −(ik/c)T. Constants are those in `crates/metrology` (CODATA 2018), so every number matches the crate tests. Checks: `notes/monopole_toroidal_checks.py` (standard library, `python3 -I`, 60 assertions, about 2 s, all pass, exits 1 on failure). Downloads are in `notes/dl_mono/` (git-ignored, untrusted). DOIs were resolved on api.crossref.org on 2026-10-27.

Labels for each claim: **known** (published; citation given), **new to our search** (derived here; no source found by the queries in §6, which is not proof that none exists), **uncertain** (source not opened, or our comparison is ambiguous).

## 1. Duality: the multipole family with electric and magnetic sources, and a duality-covariant anapole

### Takeaway
With magnetic sources there are four multipole families under parity P and time reversal T: electric (polar, T-even), magnetic (axial, T-odd), magnetic toroidal (polar, T-odd) and electric toroidal (axial, T-even). The electric toroidal family is exactly the toroidal moment of a magnetic current. In a material it is a vortex of electric polarisation, G = ½∫r × P. All of this is **known**. The radiating dipole amplitudes form a duality doublet (c𝒫, ℳ). The anapole condition "doublet = 0" is therefore duality-invariant, and the magnetic anapole m + (ik/c)T⁽ᵐ⁾ = 0 is its rotated half. Also **known** in substance. Nothing here is paper material.

### Cited Findings
- Source substitution chain, credited to Miller and quoted by Dubovik & Saha: ρ_e → j_e = c curl M → −div M = ρ_m → j_m = c curl P → −div P = ρ_e (Gaussian units). The axial toroid moment couples to the vortex electric field that a changing B induces. — Dubovik & Saha, [arXiv:cond-mat/9806057](https://arxiv.org/abs/cond-mat/9806057) (1998), §§1–2 (downloaded).
- Four fundamental multipoles under the space-time inversion group (their Table I): E (Q_lm, sources ρ_e and j_m), M (j_e), MT (j_e) and ET (j_m). The ET multipole is G_lm = [1/c(l+1)]∫[r·j_m]O_lm = ∫G·∇O_lm, with P = ∇×G. Putting j_m in place of j_e in the *magnetic* multipole gives back the *electric* Q_lm. — Hayami & Kusunose, J. Phys. Soc. Jpn. 87, 033709 (2018), [arXiv:1712.02927](https://arxiv.org/abs/1712.02927), eqs. 3–10 (downloaded).
- Magnetic anapole: the lowest-order magnetic anapole comes from interference between the magnetic toroidal dipole and the magnetic dipole. — Li & Crozier, "Origin of the anapole condition as revealed by a simple expansion beyond the toroidal multipole", [Phys. Rev. B 97, 245423 (2018)](https://doi.org/10.1103/PhysRevB.97.245423) (Crossref and search summary only).
- Duality relations among the AB, Aharonov–Casher (AC) and He–McKellar–Wilkens (HMW) phases, plus a fourth, the dual AB (DAB) phase of a monopole circling a line of electric dipoles. — Dowling, Williams & Franson, [PRL 83, 2486 (1999)](https://doi.org/10.1103/PhysRevLett.83.2486), [arXiv:quant-ph/9906016](https://arxiv.org/abs/quant-ph/9906016) (downloaded).
- Duality-invariant phase of a dyon circling a solenoid that holds both electric and magnetic flux. — Heras, "The quantum phase of a dyon", [Eur. Phys. J. Plus 138 (2023)](https://doi.org/10.1140/epjp/s13360-023-03914-5), [arXiv:1910.01117](https://arxiv.org/abs/1910.01117) (downloaded; abstract read).
- The quasistatic "magnetic and toroidal analogues of charge" expansion. — Krynytskyi & Rovenchak, [arXiv:1804.10494](https://arxiv.org/abs/1804.10494) (downloaded; no magnetic sources, there are "no magnetic charges hence B_p ≡ 0").
- Toroidal multipoles as higher-order terms of the exact electric multipoles. — Fernandez-Corbaton, Nanz & Rockstuhl, [arXiv:1507.00755](https://arxiv.org/abs/1507.00755) (already a crate source).

### Inferences
- **Doublet (derived; known in substance).** For a source with (ρ_e, J_e, ρ_m, J_m):
  - c𝒫 = c[p + (ik/c)T⁽ᵉ⁾ − d⁽ᵐ⁾/c²], with d⁽ᵐ⁾ = ½∫r × J_m dV (the electric dipole of a magnetic-current loop);
  - ℳ = m_A + m_G + (ik/c)T⁽ᵐ⁾, with m_A = ½∫r × J_e, m_G = ∫rρ_m dV, and T⁽ᵐ⁾ the 1/10 formula applied to J_m.

  The duality rotation (cJ_e, J_m) → (cJ_e cos ξ + J_m sin ξ, −cJ_e sin ξ + J_m cos ξ) rotates (c𝒫, ℳ) by ξ. Radiated power ∝ |c𝒫|² + |ℳ|² is invariant. So the joint condition 𝒫 = ℳ = 0 is duality-covariant, and ℳ = 0 with T⁽ᵐ⁾ is the magnetic anapole. **Known in substance**: Dubovik & Saha, Hayami & Kusunose, Li & Crozier. The explicit SI sign −d⁽ᵐ⁾/c² is **new to our search** only as a worked equation.
- **Sign check by two pictures of one magnet (new to our search as a test; physics known).** Describe an oscillating magnetisation M in two ways:
  - Ampère: J_e = ∇×M, which has toroidal moment t = ½∫r × M;
  - Gilbert: ρ_m = −∇·M and J_m = ∂M/∂t, so d⁽ᵐ⁾ = ṫ = −iωt.

  The doublet gives −d⁽ᵐ⁾/c² = (ik/c)t in the Gilbert picture, exactly the crate's p_eff term in the Ampère picture. The script checks this numerically with a 24-turn toroid at kR = 0.05, against its Gilbert twin: a continuous azimuthal ring of magnetic current with M = NI/(2πρ).
  - The two powers agree with µ₀ω⁶T²/(12πc⁵) to 0.9% and 0.05%.
  - Ampère(M) + Gilbert(−M) leaves 1.7×10⁻⁵ of the power, and Ampère(M) + Gilbert(+M) gives 4.017×.
  - The crate's anapole p = −(ik/c)T cancels the Gilbert ring to 7.9×10⁻⁸.
  - A duality rotation by ξ = 0.7 leaves the power of a mixed source unchanged to 10⁻¹².

  This fixes the relative sign of the magnetic-current term in the far field, G = (n × F_e) × n − (n × F_m)/c.
- **Static coupling (known).** A static electric toroidal moment couples to ∇×E = −∂B/∂t, not to any static field. That is Dubovik & Saha's "vortex electric field", and Planes et al.'s point (already in toroidal_moments.md §2) that the electric-toroidisation energy vanishes in statics. Dubovik & Saha's application of it, "aromagnetism" in organic liquids, is **uncertain/contested**: we have only their citation, not the experiment.

### Gaps
- Dubovik & Tugushev (Phys. Rep. 1990), Dubovik, Martsenyuk & Saha (PRE 61, 7087, 2000) and Li & Crozier 2018 were not opened. If a fully dyonic (ρ_m, J_m) exact multipole expansion is published, it is likely in one of them.

## 2. Monopoles inside a toroid: Dirac strings, h/e, composite and fractional charges

### Takeaway
A Dirac string *is* an infinitely thin solenoid carrying flux h/e = 2Φ₀ = µ₀g_D. A closed string is an infinitely thin toroidal solenoid of that flux. The AB phase of a charge q around it is 2π(q/e). That is invisible for every integer charge and visible (2π/3, 4π/3) for charge e/3 or 2e/3. This is why the minimum pole is 3g_D if fractionally charged particles are free, and why colour-magnetic charge is needed otherwise. All **known**. Our closed-form identity: as a "magnetisation", a closed n-unit Dirac string spanning vector area 𝒜 has toroidal moment n g_D 𝒜, which the crate's toroid formula reproduces in the thin limit (§2 Inferences). It is a gauge artefact for integer-charged matter: **new to our search**, but elementary.

### Cited Findings
- PDG 2025 review 94: a vector potential for the monopole field "must be singular along some line" (the Dirac string). Invisibility of the string for every pair of charges gives Q_E Q_M = 2π. Without colour-magnetic charge the minimum becomes 6π/e, because quark charges enter. — PDG via `topological-monopoles-and-transmutation/notes/monopoles.md` §1 (opened there).
- The string can be modelled as an infinitely thin solenoid returning the monopole's flux. Ferrell & Hopfield made the long-solenoid version of the argument (search summaries of [arXiv:2012.11941](https://arxiv.org/abs/2012.11941) and [arXiv:hep-ph/9911525](https://arxiv.org/abs/hep-ph/9911525); not opened; **uncertain** as to wording).
- The DAB phase is ϕ_DAB = −gΦ_E/ħc (Gaussian), for "a line of electric dipoles (electric flux tube)", and the HMW phase is −4πdλ_M/ħc. — Dowling et al. 1999 (downloaded, p. 2).
- In a crystalline insulator a Dirac monopole sees D = P as an effective dual magnetic field B̃ = 2πD (e = ħ = 1). The monopole's Berry phase is set by the polarisation, and the polarisation is defined only mod 1 per cell. — Song, He, Vishwanath & Wang, [Phys. Rev. Research 3, 023011 (2021)](https://doi.org/10.1103/PhysRevResearch.3.023011), [arXiv:1909.08637](https://arxiv.org/abs/1909.08637), Fig. 2 and §IV (downloaded).
- Electron AB phase 0 or π with flux in units of h/2e (Osakabe et al. 1986). — toroidal_moments.md §2.

### Inferences
All asserted in the script, §2:
- Phases round h/e: electron 2π, Cooper pair 4π, e/3 → 2π/3, 2e/3 → 4π/3 (mod 2π). A 3g_D string gives 2π again for e/3. A neutron whose three quarks all encircle the line gains Σq_iΦ/ħ = 0. A string threading *between* the quarks of a hadron would show the fractional phases, which is the physical content of PDG's 6π/e. **Known.**
- Electron round h/2e: π, the Tonomura/Osakabe value. **Known.**
- **DAB phase in SI (known).** A pole g circling electric flux Ψ = ∫D·dA (C) gains µ₀gΨ/ħ in magnitude, which is 2πΨ/e for g = g_D. Inside a ferroelectric polarisation vortex, D = P and E = 0 macroscopically, so the monopole's phase is 2π(∫P·dA)/e mod 2π. A shift of P by the polarisation quantum eR/Ω changes the flux through a cell face by exactly e, and the phase by exactly 2π.
  - So the Dirac condition and the polarisation quantum are one statement. This is Song et al.'s result in our units. I first derived it independently, then found it there: **known**.
  - Caution, from my own working: arguing "macroscopic E = 0, so no phase" is wrong. The phase is defined only mod 2π, and the integer count of dipoles cut by the spanning surface drops out of it while the fractional part does not.
- **Closed string = thin toroid (derived).** For an ideal circular-section toroid, Φ = µ₀NI(R − √(R² − a²)) exactly. The crate's T = NIπRa²/2 equals ΦπR²/µ₀ × [1 + O(a²/R²)]: ratios 0.999975 at a/R = 0.01 and 0.9975 at 0.1. So a thin flux loop has T = Φ𝒜/µ₀, and a closed n-unit Dirac string has T = n g_D 𝒜, which is 1.034×10⁻²⁶ A·m³ for a loop of radius 1 nm.
  - For integer-charged matter this "toroidal moment" is unobservable. The singular gauge transformation that removes the string is single-valued, so it also removes the contact coupling −T·(∇×B).
  - **New to our search** as a stated identity; trivially provable.

### Gaps
- Wu–Yang (1975) was not opened. Its fibre-bundle construction removes the string altogether, which makes the closed-string toroidal moment manifestly unphysical. That is consistent with the above.

## 3. Spin ice: toroidal moments of monopole currents

### Takeaway
Spin-ice strings are not quantised and are real, so the toroidal moment a monopole leaves behind is physical there. One exact, origin-independent identity covers both this case and the Dirac string: a magnetic charge Q carried round any closed path of vector area 𝒜 changes the toroidal moment ½∫r × M by Δt = Q𝒜.
- With the crate's `spin_ice_charge` (Q = 4.278×10⁻¹³ A·m for Dy₂Ti₂O₇), a flippable hexagon carries t = (Q/2)|𝒜| = 2√2 µ r_nn = 100.1 µ_B·Å.
- The smallest monopole loop (one hexagon flip) changes t by Q|𝒜| = 200.3 µ_B·Å.
- An electron crossing a single spin-ice Dirac string gains 2πQ/g_D = 0.817 mrad. Dhar et al.'s simulated holograms show "∼0.8 mrad", so that number is **known**.
- The identity Δt = Q𝒜, the hexagon values, and the reading of Dhar et al.'s "emergent electric field" phase regions as toroidal dipoles left by closed monopole loops are **new to our search**.

No published "toroidal susceptibility of a monopole gas" was found.

### Cited Findings
- Spin-ice charge q_m = 2µ/a_d ≈ q_D/8000. The strings are observable and not quantised. — Castelnovo, Moessner & Sondhi, Nature 451, 42 (2008), via monopoles.md §5 (opened there).
- Dirac strings seen by diffuse neutron scattering in Dy₂Ti₂O₇. — Morris et al., [Science 326, 411 (2009)](https://doi.org/10.1126/science.1178868) (Crossref only).
- Monopole charge and current measured ("magnetricity"); the monopole current is ∂M/∂t. — Bramwell et al., [Nature 461, 956 (2009)](https://doi.org/10.1038/nature08500) (Crossref; search summary).
- Holography simulations of a one-unit-cell Dy₂Ti₂O₇ film (cubic a = 10.0(1) Å) — Dhar, Jaubert, Cassidy, Shintake & Shannon, [arXiv:2112.01362](https://arxiv.org/abs/2112.01362), Figs. 4–5 (downloaded):
  - A polarised film with one monopole shows a phase jump across the string of about 0.8 mrad, and q_m = (4.4 ± 0.3)×10⁻¹³ A·m fitted at r = 2 nm.
  - In thermalised samples, pairs of monopoles that traverse a hexagon and annihilate leave bounded regions of phase difference ±0.4 mrad. The authors read a magnetic charge going round a loop as a "magnetic current loop", a source of electric field.
- Pollard & Zhu, J. Electron Microsc. (2013): AB phase imaging of reversal in artificial spin ice. — [abstract](https://digitalcommons.memphis.edu/facpubs/6974).
- Electric dipoles on spin-ice monopoles (3-in-1-out tetrahedra). — Khomskii, [Nat. Commun. 3, 904 (2012)](https://doi.org/10.1038/ncomms1904). Magnetoelectric response near the monopole liquid–gas transition. — Grams et al., [Nat. Commun. 5, 4853 (2014)](https://doi.org/10.1038/ncomms5853). Both Crossref only.
- Artificial "magneto-toroidal" and direct-kagome spin ices with ferrotoroidic phases — [arXiv:2404.19377](https://arxiv.org/abs/2404.19377), [arXiv:2505.12011](https://arxiv.org/abs/2505.12011) (search summaries only; **uncertain**). These are toroidal moments of nanomagnet vortices, not of monopole currents.
- "Monopolisation" in Spaldin et al., [PRB 88, 094429 (2013)](https://arxiv.org/abs/1306.5396), is the magnetoelectric monopole a = ∫r·µ dV. Its abstract does not mention spin ice.

### Inferences
All asserted in the script, §3:
- **Identity (new to our search; elementary).** Moving Q by dl flips magnetisation Q dl along the path (in spin ice, 2µ per bond of length a_d, and Q a_d = 2µ). So Δt = ½∮r × Q dl = Q𝒜, independent of origin because ∮dl = 0. It is the time integral of ½∫r × J_m, the same d⁽ᵐ⁾ = ṫ as in §1. The script checks it on a random 9-gon and again with the origin shifted.
- **Geometry.**
  - The crate's `a` = 3.54 Å is the nearest-neighbour distance r_nn. The cubic constant is 2√2 r_nn = 10.013 Å, matching Dhar's 10.0(1) Å.
  - The `spin_ice_charge` docstring calls `a` the "cubic constant", which is wrong; the formula and tests are right.
  - The diamond bond is √3a_c/4 = √(3/2) r_nn = a_d. The smallest monopole loop is a chair hexagon of six bonds, whose vector area lies along ⟨111⟩ with |𝒜| = (4/√3)a_d² = 4.341×10⁻¹⁹ m².
  - Its six spins sit on a planar regular hexagon of circumradius r_nn (checked from explicit lattice coordinates).
- **Numbers.**
  - A flippable (head-to-tail) hexagon: t = 2√2 µ r_nn = (Q/2)|𝒜| = 9.286×10⁻³² A·m³ = 100.1 µ_B·Å.
  - A flip changes it by Q|𝒜| = 200.3 µ_B·Å.
  - String AB step: 2πQ/g_D = 40πr_e/a_d = 0.8168 mrad. Dhar's 0.8 mrad agrees to 2%, and their 4.4(3)×10⁻¹³ A·m is within 0.4σ of the crate's Q.
- **Electron phase through a hexagon (derived; comparison uncertain).**
  - Continuum: a closed line of moment Q/2 per length links flux µ₀Q/2, which gives πQ/g_D = 0.408 mrad. Flipping the loop reverses it, a step of 0.817 mrad.
  - Point dipoles: the six point dipoles give ∫A·dl = (µ₀/2π)·4√2µ/r_nn along the hexagon's axis (quadrature and closed form agree to 10⁻⁶). That is 0.450 mrad, 1.103 times the continuum value, because the flux of discrete dipoles is not confined at this scale.
  - Dhar's "±0.4 mrad" regions are consistent with ±πQ/g_D about a mean. But we cannot tell from the text whether their ± refers to the two loop states, so this match is **uncertain**.
- **Reinterpretation (new to our search).** In the emergent gauge theory, Dhar et al.'s closed-loop regions are an emergent electric field. In ordinary electromagnetism they are the AB signature of a real static toroidal dipole t = Q𝒜 per loop: a flux ring µ₀Q linking the beam. While the monopole circulates, the same loop radiates as the electric dipole −d⁽ᵐ⁾/c² = (ik/c)t. Both descriptions are of one object.
- **"Toroidal susceptibility" (derived; elementary).** For a ring of any linear magnetic material (M = χH) threaded by a wire carrying I, t = ½∫r × M = χIV/(4π) exactly, the analogue of the crate's NIV/(4π). In spin ice the monopole current J_m = ∂M/∂t is what builds this azimuthal M. So the static toroidal response to ∇×H = J_free is just the transverse susceptibility. No distinct "monopole-gas toroidal susceptibility" was found in the literature. I see no reason one would differ from χ_transverse at q → 0 (**uncertain**; not derived from a spin-ice model).
- **Measurability.**
  - Per-loop: 0.4–0.8 mrad over an area of about 10 Å², near Dhar's quoted need of better than about 1 mrad phase resolution at sub-nm and sub-kelvin. Not demonstrated.
  - Macroscopic: t = χIV/(4π) is an ordinary magnetisation measurement.

### Gaps
- Morris 2009 and Bramwell 2009 full texts not opened. Dhar et al. not found in a journal (arXiv v2, Dec 2021). No spin-ice Monte Carlo was run here, so a q-dependent toroidal response of the Coulomb phase was not computed.

## 4. Aharonov–Casher and the toroid: is there a topological phase for a neutral particle with a toroidal moment?

### Takeaway
No. A point toroidal moment T moving at v through a static field gets the Lagrangian L = −µ₀ρ(r) T·v + (1/c²)d(T·E)/dt. The second term is a total derivative and drops out of every closed-loop phase. The first is a contact term with the local charge density. So the phase is −(µ₀/ħ)∫ρ T·dl: zero on any path through charge-free space, and not topological. For a thin flux loop (T = Φ𝒜/µ₀) it is exactly the AB phase of the charge swept through the loop's hole, Φq_through/ħ. That is the AB effect viewed from the solenoid's frame. The derivation is **new to our search**; the content is almost certainly what Afanasiev (1993) calls the "3-dimensional analogue of the Aharonov–Casher effect", which we could not open (**uncertain**). No unexploited configuration was found.

### Cited Findings
- AC: a magnetic dipole circling a line charge, ϕ_AC = 4πmλ_E/ħc (Gaussian), observed with neutrons and in atoms. — Aharonov & Casher, [PRL 53, 319 (1984)](https://doi.org/10.1103/PhysRevLett.53.319); Dowling et al. 1999.
- HMW: an electric dipole circling a line of monopoles. — He & McKellar, [PRA 47, 3424 (1993)](https://doi.org/10.1103/PhysRevA.47.3424); Wilkens, [PRL 72, 5 (1994)](https://doi.org/10.1103/PhysRevLett.72.5). It was measured by lithium atom interferometry, at most 27 mrad, with a 31% discrepancy from theory. — Lepoutre et al., [PRL 109, 120404 (2012)](https://doi.org/10.1103/PhysRevLett.109.120404) (search summary of the abstract).
- Afanasiev, "Toroidal solenoids in an electromagnetic field and toroidal Aharonov–Casher effect", [Phys. Scr. 48, 385 (1993)](https://doi.org/10.1088/0031-8949/48/4/001). The abstract says it derives quantum and classical equations of motion for toroidal moments in external fields, recovers a 3-D AC analogue, and predicts quantum scattering of toroidal moments on a Coulomb centre. Full text paywalled, not opened.
- A static toroid has no field outside, so it does not couple directly to static external fields (Southampton review, search summary). The toroidal moment couples to ∇×B (Planes et al.; toroidal_moments.md §2).
- Azimov & Ryndin argue that added electric-dipole or higher moments destroy the AC duality, and that AC is not purely topological (search summary of [arXiv:hep-ph/9710433](https://arxiv.org/abs/hep-ph/9710433); **uncertain**).

### Inferences
- **Derivation (new to our search).**
  1. Expanding U = −∫J·A to second order with the crate's T gives U = −T·(∇×B), using the toroidal part −2T_iδ_jk + T_jδ_ik + T_kδ_ij of ∫J_i x_j x_k from toroidal_moments.md §1.
  2. In the toroid's rest frame, B′ = −v × E/c², and ∇×(v × E) = v(∇·E) − (v·∇)E, so U′ = [(T·v)ρ/ε₀ − (v·∇)(T·E)]/c².
  3. The script checks this identity by finite differences inside a Gaussian charge blob (rel. 10⁻⁵). It checks that the total-derivative term integrates to zero round a loop through the blob (10⁻¹⁴ of its absolute sum). It checks that an interferometer arm through the blob centre gives −(µ₀/ħ)T_xQ/(2πσ²) (10⁻⁹), and a loop in charge-free space gives zero.
- **Reciprocity (new to our search as stated; physics known).** With T = Φ𝒜/µ₀, µ₀ρT·L = Φ(ρ𝒜L) = Φq_through. One electron through the hole of an h/2e thin toroid gives π, the Tonomura value seen from the other side.
- **Why it is not topological.** T couples only to the current density at its own location: µ₀J in statics, and µ₀J′ = −µ₀ρv in its rest frame. AC is topological because m × E is curl-free outside the line charge. The toroidal analogue has no field-free region in which the coupling survives.
- **Toroidal moment in a monopole field.** A static pole has ∇×B = 0 off the pole, so there is no coupling. A moving pole has J_m and a displacement current; T couples to the latter only through ∇×B = (1/c²)∂E/∂t, which is again local. **New to our search**, and a one-line consequence of U = −T·(∇×B).

### Gaps
- Afanasiev 1993 and his book were not opened, so we cannot say whether his "toroidal AC effect" is our contact term, the reciprocity statement, or something else.

## 5. Assessment, feasibility and novelty verdict

| # | Candidate | Status | Proof sketch | Numerical check (script) | Crate feasibility |
|---|---|---|---|---|---|
| 1 | Four multipole families; ET = toroidal moment of J_m = ½∫r × P | **known** (Miller via Dubovik & Saha 1998; Hayami & Kusunose 2018) | apply j_e → j_m to the MT moment | — | doc note only |
| 2 | Doublet (c𝒫, ℳ), 𝒫 ∋ −d⁽ᵐ⁾/c², ℳ ∋ (ik/c)T⁽ᵐ⁾; duality-covariant anapole | **known in substance** (duality of the textbook expansion; magnetic anapole: Li & Crozier 2018); explicit SI sign new to our search | rotate (cJ_e, J_m); Ampère and Gilbert pictures of one magnet must radiate alike | §1: 11 checks (cancellation to 1.7×10⁻⁵, anapole to 7.9×10⁻⁸, duality invariance 10⁻¹²) | possible: `toroidal::effective_dipole` could take a magnetic-current term; low value |
| 3 | Dirac string = thin toroid of h/e; phases 2π(q/e); 3g_D for free e/3 | **known** (Dirac; PDG 2025) | single-valuedness | §2: 9 checks | exists (`aharonov_bohm_phase`, `DIRAC_FLUX`) |
| 4 | Closed n-unit Dirac string has T = n g_D 𝒜 (gauge artefact) | **new to our search**, elementary | thin-toroid limit T = Φ𝒜/µ₀ | §2: ratios 0.999975, 0.9975; 1.0339×10⁻²⁶ A·m³ | not worth a function |
| 5 | DAB phase 2πΨ/e; polarisation quantum ↔ Dirac quantum | **known** (Dowling et al. 1999; Song et al. 2021) | D = P in a vortex; ΔP·A_face = e | §2: 3 checks | easy: `monopole::dual_aharonov_bohm_phase(g, psi) = µ₀gΨ/ħ`, tests 2π at (g_D, e) |
| 6 | Δt = Q𝒜 for any magnetic charge transported round a loop | **new to our search** (likely folklore) | Δ M = Q dl along the path; ½∮r × dl = 𝒜 | §2: 6 checks; §3: monopole-current form | `toroidal::toroidal_dipole` already does it if fed (r, Q dl) |
| 7 | Spin-ice string AB step 2πQ/g_D = 0.817 mrad | **known** (Dhar et al. 2021, "∼0.8 mrad") | flux µ₀Q | §3 | easy: `monopole::spin_ice_string_phase(moment, a)`, test 0.8168 mrad and Dhar's 0.8 |
| 8 | Hexagon toroidal moment 100.1 µ_B·Å; flip Δt = 200.3 µ_B·Å; reading Dhar's loop regions as toroidal dipoles | **new to our search** | explicit pyrochlore coordinates | §3: 9 checks | example, not library |
| 9 | On-axis discrete phase 0.450 mrad vs continuum 0.408 mrad | derived; match to Dhar's ±0.4 **uncertain** | point-dipole A integral | §3 | no |
| 10 | Toroidal susceptibility of a monopole gas | none published found; static t = χIV/(4π) elementary | M = χH round a wire | §3 | no |
| 11 | Toroidal "AC" phase −(µ₀/ħ)∫ρT·dl: local, = AB of swept charge | **new to our search** as derived; probably Afanasiev 1993 (**uncertain**) | Lorentz-boost U = −T·∇×B | §4: 7 checks | possible `toroidal::moving_toroid_phase(t, rho, path)`; low value |

**Novelty verdict (frank).** Nothing found here is a publishable new result in physics. The duality structure (1, 2), the string-as-solenoid and composite-charge phases (3), the dual AB phase and its tie to the polarisation quantum (5), and the spin-ice string phase (7) are all in the literature. The items new to our search (4, 6, 8, 11) are elementary identities, each provable in a few lines. Their likely status is folklore or a special case of known papers (Afanasiev 1993 for 11), not discoveries. The most a paper could be is a short pedagogical note, for example in Eur. J. Phys. or Am. J. Phys.: "Monopoles leave toroidal moments behind", unifying 4, 6, 8 and 11, with the Ampère/Gilbert radiation check of 2. Even that needs Afanasiev 1993, Dubovik & Tugushev 1990 and the published version of Dhar et al. read first, to make sure it is not already there.

**Strongest honest claim.** A magnetic charge Q carried round a closed path of vector area 𝒜 changes the toroidal moment by exactly Q𝒜, whatever the origin.
- For a Dirac monopole this is a closed h/e string, n g_D 𝒜, invisible to integer charges.
- In spin ice it is a real toroidal dipole: 200.3 µ_B·Å per hexagon flip in Dy₂Ti₂O₇, with an electron AB step of 2πQ/g_D = 0.82 mrad.
- Dhar et al.'s simulated "closed-loop" holography regions are these toroidal dipoles.
- While the loop is being traversed it radiates as the electric dipole −d⁽ᵐ⁾/c² = (ik/c)t, the same term the crate uses for Ampère toroids.

Proved here: by derivation, and by 60 numerical assertions. Not proved: that nobody has written it down.

**Crate notes (no crate code edited).**
- `magnetics::monopole::spin_ice_charge`'s docstring says `a` is the "cubic constant". The value used (3.54 Å) and the formula a_d = √(3/2)·a are for the nearest-neighbour distance; the cubic constant is 10.01 Å. The tests are correct; only the wording is wrong.
- Worth adding, if anything: `dual_aharonov_bohm_phase` and `spin_ice_string_phase`, both one-liners with published test values (2π; Dhar's 0.8 mrad).

## 6. Searches run (2026-10-27)

WebSearch queries, verbatim in substance:
1. toroidal moment magnetic current duality "magnetic toroidal" dual multipole expansion magnetic charges
2. Dowling Williams Franson "Maxwell duality, Lorentz invariance, and topological phase" Aharonov-Casher He-McKellar
3. Afanasiev "toroidal" "Aharonov-Casher" Physica Scripta 1993
4. spin ice monopole toroidal moment hexagon loop magnetoelectric
5. magnetic monopole Aharonov-Bohm phase polarization quantum Berry phase "polarization" Dirac quantization ferroelectric dual
6. anapole toroidal dipole neutral particle topological phase Aharonov-Casher analogue "toroidal moment" interferometry
7. "Dirac string" toroidal solenoid infinitely thin flux h/e invisible closed loop toroidal moment
8. "spin ice" "toroidal susceptibility" OR "toroidal moment" monopole current Dirac string pyrochlore
9. Lehmann Donnelly Derlet Heyderman Fiebig "magneto-toroidal" artificial spin ice toroidal moment susceptibility
10. He McKellar 1993 … Wilkens 1994 … Lepoutre HMW measurement 2012
11. Spaldin Fiebig Mostovoy 2008 toroidal moment magnetoelectric … E×H
12. "toroidal" "Aharonov-Casher" Afanasiev toroidal moment moving charged wire phase "Topological Effects in Quantum Mechanics"
13. anapole moment interaction Hamiltonian moving toroidal dipole electric field "v ×" Lorentz transformation …
14. "dual Aharonov-Bohm" monopole "polarization" ferroelectric vortex OR "electric flux tube" realization quasiparticle spin ice
15. "polarization quantum" "magnetic monopole" phase "modulo" Dirac quantization Vanderbilt Witten effect …
16. electron Aharonov-Bohm phase spin ice Dirac string flux "not quantized" electron holography …
17. "magnetic anapole" OR "dual anapole" magnetic dipole cancelled by magnetic toroidal dipole duality condition
18. moving magnetic monopole closed path leaves toroidal moment "vector area" string magnetization spin ice loop flip toroidal dipole

Queries 4, 6, 7, 8, 18 found nothing that combines toroidal moments with monopole currents or Dirac strings; that is the basis of every "new to our search" label above. Downloaded to `dl_mono/`: arXiv quant-ph/9906016, 1712.02927, 1804.10494, 1910.01117, cond-mat/9806057, 1909.08637, 2112.01362. Opened only as abstracts or Crossref records: Afanasiev 1993, Spaldin et al. 2008 and 2013, Khomskii 2012, Grams 2014, Morris 2009, Bramwell 2009, Li & Crozier 2018, Lepoutre 2012, Pollard & Zhu 2013.
