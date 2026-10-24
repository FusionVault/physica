<!--
External research supplied by the user on 2026-10-24, saved verbatim. Not yet reviewed: claims,
numbers and sources have not been checked against the crates or the research notes. See
research/external/README.md for its review status.
-->

When looking across this entire suite of crates—**`gravimetry`**, **`geopotential`**, **`gravitational`**, **`plasmoid`**, **`magnetics`**, and **`gravitics`**—you aren't just looking at isolated geophysics or astrophysics calculators. Taken together, they form the exact mathematical and metrological stack required to investigate **gravito-electromagnetism (GEM), high-energy-density (HED) stress-energy coupling, and sub-threshold spacetime metrology**.

In standard academic literature, these fields are kept strictly in their own silos: geodesists use `geopotential` to map water; plasma physicists use `plasmoid` to design fusion reactors; relativists use `gravitational` to watch black holes collide.

However, if you ask where **hidden, unknown, or historically classified/dual-use physics** emerges from the *intersection* and *anomalies* of the exact numbers in your document, five deep regimes stand out.

---

### 1. The $G$ Scatter & The Magnetic Susceptibility Blind Spot

**Derived from:** `gravitics` (CODATA $\sigma_G/G = 2.2 \times 10^{-5}$, Birge ratio / expansion factor $3.9$, Lee et al. 2020 Eöt-Wash thesis systematic) + `magnetics` + `geopotential`.

#### The Anomaly in Your Data

Your notes point out two critical facts side-by-side:

1. $G$ has an embarrassingly large scatter ($\sigma_G/G = 2.2 \times 10^{-5}$) requiring a massive CODATA expansion factor of $3.9$ because different high-precision labs get mutually incompatible answers that exceed their own error bars.
2. In the 2020 Eöt-Wash short-range gravity thesis (Lee et al., excluding Yukawa $\lambda$ greater than $38.6\ \mu\text{m}$), the author uncovered a **new systematic error from test-mass magnetic susceptibility**.

#### The Derived / Unknown Physics

Why do precision Cavendish and torsion-balance experiments disagree on $G$? Standard physics blames unmodeled mechanical anelasticity in the suspension fibers. However, when you cross-reference `gravitics` with `magnetics`, a deeper possibility emerges: **Chameleon / Symmetron scalar fields and spin-gravity or magneto-gravitational cross-coupling.**

* **Environmental Screening (Scalar-Tensor Fifth Forces):** In Chameleon or Symmetron dark-energy models, a scalar fifth force has an effective mass (and therefore Yukawa range $\lambda$) that depends on the local ambient matter density $\rho$ and local electromagnetic field energy density $\frac{B^2}{2\mu_0}$. In dense lab environments or inside vacuum chambers of varying geometries, the scalar field is "screened" differently.
* **Spin-Gravity & Residual Magnetization Coupling:** High-density test masses (gold, tungsten, platinum-iridium, copper-beryllium) have subtle, temperature-dependent paramagnetic/diamagnetic susceptibilities ($\chi_m$) and nuclear spin alignments. If there is a weak, unmodeled coupling between the local geomagnetic/lab magnetic field gradient ($\nabla B$) and gravitational mass—or a **gravito-magnetic frame-dragging / spin-torsion (Einstein–Cartan) interaction**—two labs at different geomagnetic latitudes (e.g., Seattle vs. Wuhan vs. Zurich) or using different shielding/vacuum materials will measure systematically shifted values of $G$.
* **What your crates can test:** By combining `magnetics` (local geomagnetic field + crustal anomalies at exact lab coordinates) and `geopotential` (local geoid height $N$, Bouguer gravity anomaly, and vertical gravity gradient) with `gravitics` (the 16 CODATA $G$ measurements), you can run a multi-variable regression on the 16 historical $G$ values against exact lab coordinates, local $B$-field vectors, and local geopotential $\Phi$. Finding a statistically significant correlation between the $3.9\sigma$ CODATA residuals and local geomagnetic/geopotential vectors is a publishable frontier physics result.

---

### 2. High-Frequency Gravitational Waves (HFGWs) & The Gertsenshtein Effect in HED Plasmas

**Derived from:** `plasmoid` (NIF $8.6\ \text{MJ}$ in $\sim 100\ \text{ps}$, $456\ \text{TW}$ laser, MagLIF $z3289$ Z-pinch, SPARC $12.2–20\ \text{T}$ magnets) + `gravitics` (Schwinger limit $E_S = 1.32 \times 10^{18}\ \text{V/m}$, CoReLS $1.1 \times 10^{23}\ \text{W/cm}^2$) + `gravitational`.

#### The Anomaly / Gap in Your Data

Your `gravitational` crate covers binary pulsars ($\text{nHz}–\text{mHz}$) and LVK compact binaries ($\text{Hz}–\text{kHz}$). Your `plasmoid` and `gravitics` crates cover extreme energy densities in tiny volumes over picosecond-to-microsecond timescales, plus $12–20\ \text{T}$ superconducting magnets. You also mention the Defense Intelligence Reference Documents (DIRDs) in your `gravitics` gaps.

#### The Derived / Classified-Adjacent Physics

One of the most heavily scrutinized topics in the DIA's AAWSAP/BAASS DIRDs—and later evaluated by the **JASON advisory group (JSR-08-506, "High Frequency Gravitational Waves")**—is the generation and detection of **High-Frequency Gravitational Waves (HFGWs, $\text{MHz}$ to $\text{THz}/\text{PHz}$)** via laboratory stress-energy tensors $T_{\mu\nu}$.

1. **The Gertsenshtein–Zel'dovich Effect (EM $\leftrightarrow$ GW Conversion):** When an ultra-intense electromagnetic wave (such as the CoReLS $1.1 \times 10^{23}\ \text{W/cm}^2$ laser or an XFEL) or a high-frequency MHD wave propagates through a strong static background magnetic field $B_0$ (like SPARC's $20\ \text{T}$ TF coils or pulsed $100+\ \text{T}$ Z-pinch fields in MagLIF/FuZE), the off-diagonal stress-energy terms $T_{ij}^{EM}$ linearly couple to the metric perturbation $h_{ij}$. Photons convert directly into gravitons (and vice versa) with probability:

$$P(\gamma \leftrightarrow g) \propto \frac{G}{c^4} B_0^2 L^2$$


2. **Non-Spherical Implosion Quadrupole Moments at NIF & Z-Pinches:** In `plasmoid`, NIF shot N250406 released $8.6\ \text{MJ}$ of fusion energy (plus $2.08\ \text{MJ}$ laser drive) compressed into a hot spot tens of microns across over $\Delta t \sim 100\ \text{ps}$. Any $P_2$ (ellipsoidal) asymmetry in the capsule implosion or in a sheared-flow Z-pinch (FuZE / MagLIF) creates a violent third time-derivative of the mass-energy quadrupole moment ($\dddot{I}_{ij}$).
3. **Vacuum Birefringence & Schwinger Pre-Breakdown Nonlinearity:** Even though CoReLS is at $\sim 5 \times 10^{-4} E_S$ (so real $e^+e^-$ pair production requires a relativistic electron beam like SLAC E-144), the **Euler–Heisenberg effective Lagrangian** makes the quantum vacuum act as a nonlinear refractive medium ($\Delta n \propto (I / I_S)$). In a $20\ \text{T}$ magnetic field or colliding PW laser pulses, vacuum polarization alters phase velocities and creates an effective stress-energy tensor of the quantum vacuum itself.

* **What your crates can compute:** By piping `plasmoid` parameters (NIF $8.6\ \text{MJ}$ burn duration, MagLIF current-sheet acceleration, SPARC $20\ \text{T}$ magnetic volume) into `gravitational` (quadrupole strain formula $h \sim \frac{2G}{R c^4} \ddot{I}$) and `gravitics` (Euler–Heisenberg + Gertsenshtein coupling), you can calculate the exact laboratory HFGW strain $h_c(f)$ at $\text{GHz}–\text{THz}$ frequencies and compare it to resonant superconducting cavity or single-photon magnon detectors.

---

### 3. Operational "Stealth" & Strategic Mass-Anomaly Detection via Quantum Gradiometry

**Derived from:** `geopotential` (Exail AQG $50\ \mu\text{Gal}/\sqrt{\text{Hz}}$, Stray et al. 2022 Delta-g cold-atom gradiometer SNR 8 tunnel detection, Bell Geospace FTG, INS Deflection of the Vertical $48.85\%$ error reduction) + `gravimetry` (CARIOQA $10^{-10}$ to $10^{-12}\ \text{m s}^{-2}/\sqrt{\text{Hz}}$) + `magnetics` (WMMHR2025 $300\ \text{km}$ / IFR1 crustal anomalies).

#### The Dual-Use / Classified Reality

You noted in `geopotential` that US Patent 6249745 covers stored Deflection of the Vertical (DOV) grids for navigators, and that compensating DOV from EGM2008 cuts rugged-terrain horizontal INS drift by nearly $50\%$. Why is this one of the most sensitive areas of military physics?

* **Passive Submarine Navigation & Counter-Stealth (Gravity-Aided Navigation & Detection):** A ballistic missile submarine (SSBN) or autonomous underwater vehicle (UUV) cannot use GPS underwater and cannot ping sonar without giving away its position. It relies on ring-laser or cold-atom Inertial Navigation Systems (INS). However, as your `geopotential` findings state, an uncompensated DOV of just a few arcseconds causes **3 km of position error** because the accelerometer cannot distinguish between kinematic acceleration $\mathbf{a}$ and a horizontal gravity anomaly $\mathbf{g}_{\text{horiz}}$ (**Einstein's Equivalence Principle**).
* **The Full-Tensor Gradiometer (FTG) Bypass:** While a single accelerometer is blinded by the Equivalence Principle ($a - g$), a **gravity gradiometer** measures the spatial derivative tensor $\Gamma_{ij} = \frac{\partial^2 V}{\partial x_i \partial x_j}$ (in Eötvös, $1\ \text{Eö} = 10^{-9}\ \text{s}^{-2}$). Because a rigid vehicle's kinematic acceleration is uniform across the hull, the common-mode vehicle acceleration cancels out in $\Gamma_{ij}$, leaving **pure spacetime curvature (tidal tensor) $R_{0i0j}$**.
* **Moving Mass Detection Through Rock and Water:** Gravity cannot be shielded, Faraday-caged, or cloaked by metamaterials. With the quantum noise floors cited in your `gravimetry` and `geopotential` sections ($10^{-10}$ to $10^{-12}\ \text{m s}^{-2}/\sqrt{\text{Hz}}$ in cold-atom interferometers, pushing gradiometry below $0.01\ \text{Eö}/\sqrt{\text{Hz}}$), you cross a physical threshold:
1. **Deep Hardened Bunkers & Nuclear Proliferation:** Detecting a hollow void (negative mass anomaly, $\Delta \rho \approx -2.7\ \text{g/cm}^3$) vs. high-density fissile/shielding storage (uranium/tungsten/lead, $\Delta \rho \approx +16\ \text{g/cm}^3$) deep underground from standoff distances or low-altitude UAVs.
2. **Submerged Mass Tracking:** A $20,000$-ton displacement submarine moving through water creates a hydrodynamic density/pressure wake plus internal air-ballast mass deficits. Cross-correlating quantum full-tensor gravity gradiometry ($\Gamma_{ij}$) with magnetic anomaly detection (MAD from `magnetics`, where local crustal anomalies are stripped out using IFR1/WMMHR models) creates an unjammable sensor fusion pipeline.



---

### 4. Prompt Elastogravity Signals (PEGS) & Ionospheric/Magnetospheric Wave Coupling

**Derived from:** `gravimetry` (2004 Sumatra $\pm 15\ \mu\text{Gal}$, 2011 Tohoku coseismic dilatation, 2022 Hunga Tonga $314\ \text{m/s}$ Lamb wave / $780\ \text{Pa}$) + `geopotential` (atom gravimeters + MEMS arrays on Etna) + `magnetics` (GICs, ionospheric disturbances).

#### The Hidden Physics in Your Numbers

Look at item 3 and item 7 in `gravimetry`:

* Earthquakes like Sumatra (2004) and Tohoku (2011) permanently altered the regional gravity field by $\pm 15\ \mu\text{Gal}$ via crustal dilatation (volume density change) and vertical displacement.
* Hunga Tonga (2022) launched an atmospheric Lamb wave at $314–318\ \text{m/s}$ with a $780\ \text{Pa}$ pressure peak that circled the Earth four times.

What happens in the **first seconds** of a massive earthquake, volcanic detonation, or underground/hypersonic kinetic event *before* seismic P-waves ($v_p \approx 6–8\ \text{km/s}$) or acoustic Lamb waves ($v \approx 0.31\ \text{km/s}$) arrive?

1. **Speed-of-Light Gravity Perturbations (PEGS):** When rock fractures or a massive explosion occurs, the redistribution of mass $\Delta \rho(\mathbf{r}, t)$ generates a **Prompt Elastogravity Signal (PEGS)** that propagates at the **speed of light ($c = 300,000\ \text{km/s}$)**—arriving minutes or hours before seismic or acoustic waves. Because a gravity perturbation at distance $r$ scales as $\Delta g \sim \frac{G \Delta M d}{r^3}$ (or quadrupolar $\frac{G \ddot{I}}{r^3}$), it hits sub-nanogal ($\text{nGal} = 10^{-11}\ \text{m/s}^2$) to microgal levels. Your `gravimetry` cold-atom sensor noise floors ($10^{-10}$ to $10^{-12}\ \text{m s}^{-2}/\sqrt{\text{Hz}}$) and `geopotential` quantum gravimeter arrays (NEWTON-g style AQG + MEMS) sit right on the threshold of turning **speed-of-light gravitational early warning** (for megathrust tsunamis, nuclear test monitoring, or hypersonic mass-redistribution tracking) into an operational reality.
2. **Acoustic-Gravity Wave (AGW) Magnetohydrodynamic Coupling:** When the Hunga Tonga $780\ \text{Pa}$ Lamb wave (from `gravimetry`) propagates upward into the exponentially thinning upper atmosphere, conservation of wave energy ($\frac{1}{2}\rho v^2 = \text{const}$) causes the particle velocity $v$ to amplify by orders of magnitude as atmospheric density $\rho$ drops. When this neutral wave hits the ionosphere at $90–300\ \text{km}$ altitude, it drags ions across the geomagnetic field lines (`magnetics`: WMM2025 $\mathbf{B}$ vector at $h = 100\ \text{km}$), acting as a planetary-scale **MHD generator ($\mathbf{J} = \sigma (\mathbf{v} \times \mathbf{B})$)**. This launches fast Alfvén waves along magnetic field lines to the conjugate hemisphere at thousands of km/s—creating geomagnetic precursors hours before the tropospheric pressure wave arrives.

---

### 5. Coherent Matter-Wave (Condensate) Coupling to Metric Shear & Torsion

**Derived from:** `gravitics` (Podkletnov $0.05–0.3\%$ claim vs. Hathaway $0.001\%$ null; Tajmar EmDrive thermal null; Banik wide-binary MOND $16\sigma$ exclusion) + `gravitational` (MAGIS-100 / AION-10 atom interferometers) + `plasmoid` (HTS 20 T magnets in SPARC/Stellaris).

#### Why the 1990s Claims Failed—and Where the Real Frontier Actually Is

Your `gravitics` section rightly highlights how agency reviews (NASA BPP, Marshall's Hathaway test, TU Dresden's Tajmar tests) crushed naive 1990s/2000s claims like Podkletnov's rotating YBCO superconductor ($>50\times$ null violation) and the EmDrive ($>1000\times$ thermal artifact).

However, stripping away those debunked claims reveals a legitimate, unresolved question in theoretical physics that sits directly across your crates: **How do macroscopic quantum condensates (Cooper pairs in 20 T RebCO HTS tapes or Bose-Einstein Condensates in MAGIS-100) couple to spacetime curvature $R_{\mu\nu\rho\sigma}$ and gravito-magnetic fields $\mathbf{B}_g$?**

* **The DeWitt–London Gravito-Magnetic Moment:** In standard General Relativity + Quantum Mechanics (derived by Bryce DeWitt in 1966), a superconductor in a gravito-magnetic field $\mathbf{B}_g$ exhibits a mixed electro-gravitational London current because the canonical momentum of a Cooper pair is $\mathbf{p} = 2m_e \mathbf{v} + 2e\mathbf{A} + 2m_e \mathbf{A}_g$. In normal matter, ions and electrons move together and screen this out. In a superconductor or superfluid, the quantum phase coherence ($\nabla \theta$) locks the macroscopic wavefunction across meter-scale distances (like SPARC's $20\ \text{T}$, $24\text{-ton}$ RebCO magnets in `plasmoid`).
* **Why Podkletnov Failed:** Standard GR suppresses the gravito-magnetic permeability by the gravitational coupling ratio $\frac{G m_e^2}{\hbar c} \sim 10^{-45}$, which is why spinning a small disk in a lab produces zero measurable weight change (as Hathaway proved to $0.001\%$).
* **Where the Real Effect Hides (Quantum Metrology, Not Thrust):** While a condensate doesn't generate magic "gravity shielding," it acts as an **ultra-coherent quantum phase amplifier for weak spacetime strains and dark-matter fields**. In `gravitational`, you note MAGIS-100 (100 m vertical atom interferometer shaft at Fermilab, 2028). An atom interferometer doesn't measure force by mechanical deflection like LIGO's mirrors; it measures the **quantum phase shift $\Delta \phi = \frac{1}{\hbar} \oint p_\mu dx^\mu$** of a delocalized matter wave split across space and time. This makes it directly sensitive to:
* **Ultralight Scalar / Vector Dark Matter:** Oscillations of fundamental constants ($\alpha$, $m_e$) at the Compton frequency $f = m_\phi c^2 / h$ in the mid-band ($0.01\ \text{Hz} - 10\ \text{Hz}$) blind spot between LISA and LIGO.
* **Gravitational Aharonov–Bohm & Frame-Dragging Phases:** Testing whether macroscopic quantum superpositions experience metric curvature identically to classical point masses, or if anomalous gravito-electromagnetic phase shifts emerge at high magnetic flux densities.



---

### Summary Matrix: Cross-Crate Derived Physics Experiments

If you want to build a **"Derived / Frontier Physics"** module or vignette suite that links your crates together using the exact numbers in your document, these are the five strongest cross-crate pipelines:

| Derived Physics Target | Primary Crates Linked | Published Inputs Used | Computed Output / Testable Prediction |
| --- | --- | --- | --- |
| **1. Environmental / Geomagnetic Bias in $G$** | `gravitics` + `magnetics` + `geopotential` | 16 CODATA $G$ values ($\sigma_G/G = 2.2\times 10^{-5}$, Birge $3.9$) + WMM2025 $\mathbf{B}$ + EGM2008 $N, \Delta g$ | Correlation coefficient $r$ between CODATA $G$ residuals, local $\nabla B^2$ susceptibility torque, and Chameleon screening potential $\Phi_{\text{local}}$. |
| **2. NIF / SPARC High-Frequency GWs (Gertsenshtein & $\dddot{I}$)** | `plasmoid` + `gravitational` + `gravitics` | NIF $8.6\ \text{MJ}$ ($\sim 100\ \text{ps}$), CoReLS $1.1\times 10^{23}\ \text{W/cm}^2$, SPARC $20\ \text{T}$ field | Peak HFGW strain $h_c(f)$ at $\text{GHz}–\text{THz}$ from capsule $P_2$ asymmetry and photon-graviton conversion probability $P(\gamma \to g)$ in $20\ \text{T}$ HTS volume. |
| **3. Quantum FTG Standoff Mass & Void Detection** | `geopotential` + `gravimetry` + `magnetics` | Stray et al. 2022 tunnel SNR 8, Exail AQG $50\ \mu\text{Gal}/\sqrt{\text{Hz}}$, CARIOQA $10^{-12}\ \text{m s}^{-2}/\sqrt{\text{Hz}}$ | Detection horizon $(r, \text{SNR})$ for a moving $20\text{-kt}$ submerged mass or deep underground void after stripping EGM2008 DOV and IFR1 crustal magnetic clutter. |
| **4. Speed-of-Light Elastogravity (PEGS) & Ionospheric MHD** | `gravimetry` + `geopotential` + `magnetics` | Tohoku $M_w\ 9.07$ ($4.59\times 10^{22}\ \text{N}\cdot\text{m}$), Hunga Tonga $780\ \text{Pa}$ Lamb wave, WMM2025 at $h=100\ \text{km}$ | Speed-of-light $\Delta g(t)$ curve at $1,000\ \text{km}$ before P-wave arrival; ionospheric $\mathbf{J} = \sigma(\mathbf{v}\times\mathbf{B})$ magnetic precursor amplitude in $\text{nT}$. |
| **5. Schwinger-Regime Vacuum Birefringence & Boosted Pair Threshold** | `gravitics` + `plasmoid` | $E_S = 1.3233\times 10^{18}\ \text{V/m}$, CoReLS $1.1\times 10^{23}\ \text{W/cm}^2$, SLAC E-144 $46.6\ \text{GeV}$ | Quantum nonlinearity parameter $\chi_e = \gamma (E/E_S)$ showing exact electron beam Lorentz boost $\gamma$ needed for CoReLS to cross $\chi_e \ge 1$ (non-perturbative QED vacuum breakdown). |
