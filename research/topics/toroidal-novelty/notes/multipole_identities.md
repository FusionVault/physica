# The algebra of toroidal and hypertoroidal moments of closed currents: what is known, what can be proved, and what is new

Scope: the reduction of toroidal-family moments (toroidal dipole, toroidal quadrupole and octupole, mean-square radii, hypertoroidal/supertoroidal ladders) to single integrals of J for closed currents, the conversion factors between published conventions, the 3/2 in Planes, Castán & Saxena, and numerical checks. Asked for issue #23's novelty investigation.

Conventions: SI without the Gaussian 1/c; "closed current" means J localised (compact support or fast decay) with ∇·J = 0. Dubovik's toroidal dipole is T = (1/10)∫[(r·J)r − 2r²J] dV. STF means the symmetric trace-free part.

Labels: **known** (with citation to a source opened here), **known (not opened)** (cited by an opened source), **new to our search** (not found in the sources and queries listed in §5; that is not proof of novelty), **uncertain**.

Downloads (untrusted, not executed) are in `notes/dl_identities/`. Checks: `notes/multipole_identities_checks.py`, numpy, run with `python3 -P`. It makes 61 assertions, runs in about 2 s, all pass, and exits 1 on any failure.

**Verdict in one line:** the reductions, including the general-order theorem, are textbook or published. The ladder normalisations are an elementary corollary that is new to our search. Two misprints are worth a note to the authors: Planes et al.'s 3/2, and Gurvitz et al.'s 1/10 for T⁽²ᵉ⁾. None of this is a research paper. At most it is a short pedagogical note or comment.

| Candidate | Verdict | Evidence |
|---|---|---|
| (b) T = ½∫(r·J)r = −¼∫r²J for closed J | **known**, textbook | Haxton, Liu & Ramsey-Musolf 2002 eqs. 4–6; Dumitriu & Vrejoiu 2005; Ederer & Spaldin 2007 eq. 5 |
| One-term form of the toroidal 2^L-pole at every L | **known** | Dubovik & Tugushev 1990 eq. B.4 (not opened, cited by Szmytkowski & Łukasik); Agre 2011 eq. 21; Szmytkowski & Łukasik 2016 App. C (proof plus a one-parameter family) |
| Same for all mean-square radii (r^{2p} weights) | one-line corollary of the same identity; explicit general statement **new to our search** | §2, Theorem B; checked to p = 3 |
| Counting: N(N+2) independent order-N moments | orders 1–2 **known** (Haxton et al.); general case is a Clebsch–Gordan exercise | checked numerically to N = 4 |
| Generating function | **known**: Jackson-type exact multipole coefficient with ρ = 0; Alaee et al. 2018 exact dipoles | the closed form (3/2)∫J j₀(kr) is a corollary (§2) |
| Ladder h_n = ∫K for J = (∇×)ⁿK, closed forms for all n | ladder definition **known** (Boardman & Marinov 2005); MSR ↔ high-order toroidal "equivalent up to a prefactor" **known** (Borovkov & Canós Valero 2021); explicit all-n static closed forms and ratios to the iterated-r× convention **new to our search** | §2, Theorem D; checked n = 1–6 |
| Planes et al. eq. 29/30: ¼∫r×(r×J) = (3/2)T | **erratum-level**, no published erratum found; no published numbers found to depend on it | §3 |
| Gurvitz et al. 2019 Table 2: T⁽²ᵉ⁾ with 1/10 | **uncertain → likely misprint**: their own eq. 23, with Alaee's exact dipole, needs 1/280 | §3, check |

## 1. Are (b) and its general-order versions known?

### Takeaway
Yes. (b) is stated, with the warning that it needs current conservation, in the nuclear-anapole literature. The general-order statement, that for a closed current every toroidal 2^L-pole is a single integral, is in Dubovik & Tugushev's 1990 review, according to Szmytkowski & Łukasik, who cite its eq. B.4. It is also in Agre's 2011 Physics–Uspekhi methodological note. Szmytkowski & Łukasik prove it, with a one-parameter family of equivalent forms. The mean-square-radius versions are not written out anywhere we found, but they follow from the same one-line identity.

### Cited Findings
- **Haxton, Liu & Ramsey-Musolf** (nucl-th/0109014, PRC 2002) define the anapole as a = (M²/6)∫x×(x×j) (eq. 4) and write: "Current conservation allows Eq. (4) to be rewritten as" a = −(M²/4)∫x²j (eq. 6). They add that this form "is obtained only after exploiting the constraints of current conservation". They count six constraints on ∫xᵢjⱼ, leaving three components (the magnetic dipole), and ten constraints on the 18 trilinear moments, leaving eight (the magnetic quadrupole plus the anapole). In Lorentz–Heaviside units a/M² is exactly Dubovik's T. — `dl_identities/nucl-th_0109014.txt` §II.A
- **Dumitriu & Vrejoiu** (physics/0503164) write the toroidal dipole as T = (1/6)∫ξ×(ξ×j) − (1/30)Ṗ_qqi, with Ṗ_qqi = ∫ξ²ξᵢρ̇. Using ∫ξ²ξᵢρ̇ = ∫[2(ξ·j)ξᵢ + ξ²jᵢ], they reach Dubovik's (1/10) form. They also give the toroidal quadrupole Tᵢₖ = (1/42)∫[4(ξ·j)ξᵢξₖ − 5ξ²(ξᵢjₖ + ξₖjᵢ) + 2ξ²(ξ·j)δᵢₖ]. — `physics_0503164.txt` p. 10
- **Ederer & Spaldin** (0706.1974, PRB 76, 214404) give t = (1/10c)∫[r(r·j) − 2r²j] (eq. 1) and, for j = c∇×µ, t = ½∫r×µ (eq. 5). For localised spins this gives t = ½Σ r_α × m_α (eq. 7), and under a shift of origin t → t + ½R₀×m. — `0706.1974.txt`
- **Fernandez-Corbaton, Nanz & Rockstuhl** (1507.00755), appendix: the same divergence-theorem step on J r² r_k is what produces the 1/10 (search snippet; the earlier toroidal-moments note opened this paper for eqs. 5–6).
- **Agre** (Phys. Usp. 54, 167, 2011): the static toroidal 2^l-pole is a single integral, T_lm ∝ ∫r^{l+1} j·Y_{l,l+1,m}, eq. 21. Dubovik's two-term expression, "on the face of it, does not coincide", but the two are "equivalent under current stationarity conditions". Static toroidal moments make no contribution to B. — `ufn112d.txt` pp. 170–171
- **Szmytkowski & Łukasik** (1602.00613), App. C. For a sourceless stationary current they prove that T_LM, defined with Y^{L+1}_{LM} (C.1), equals the Y^{L−1}_{LM} form (C.8) and a one-parameter family (C.9). The family includes T_LM = (1/(L+1))√(4π/(2L+1)) ∫r^L Y_LM (r·j) d³r (C.13), "given before in Refs. [40, Eq. (B.4)] and [25, Eq. (2.2)]" (Dubovik & Tugushev 1990; Mielewczyk & Szmytkowski 2006). It also includes forms in r×j (C.16, C.18). They show T₀₀ = ∫r·j = 0. — `1602.00613.txt` App. B–C
- **Krynytskyi & Rovenchak** (1804.10494) give toroidal moments of every order as T_{i₁…iₙ} = ∫ρ_t(n) x_{i₁}…x_{iₙ}, with a "toroidal charge" ρ_t(n) = −div{2(j·r)r − (n+3)r²j}/(2c(n+1)(2n+3)) (eqs. 56–57). This reproduces the 1/10 and 1/42 tensors (eqs. 59–60). They note that another reference uses 28 where they use 42 (the "Answer" remark). — `1804.10494.txt` §6
- **Góngora & Ley-Koo** (Rev. Mex. Fís. E 52, 188, 2006): the exact electric multipole source is ∫j_l(kr′)Y*_lm[ik r′·J − c(2 + r′·∇′)ρ] (their eq. 42 and its long-wavelength form, eq. 47). With ρ = 0 only the one-term (r·J) j_l weight is left. They note a missing factor of k in another reference's eq. 5b. — `gongora_leykoo2006.txt` pp. 193–194
- **Dubovik & Cheshkov 1974** (JINR archive scan, first 10 pages only): J⊥ = rot N + rot rot R, with N = xψ and R = xχ (eqs. 15–17), and the magnetic dipole's mean-square radius ρ̄² = ½∫ξ²[ξ×J]dξ. — `dubovik_cheshkov1974.pdf` pp. 797, 800

### Inferences
- (b) and its general-L version are known, and the general-L proof is in print (Szmytkowski & Łukasik App. C). Our script confirms the spherical statement in Cartesian form, ∫r^{2p+2}J·∇h = −(2p+2)∫r^{2p}(r·J)h for harmonic h, for L = 1–5 and p = 0–2 on a tilted toroidal winding and a trefoil knot, to rounding error.
- The mean-square-radius versions follow from the same identity with f = r^{2p+2}h. For L = 1 this reads ∫r^{2p}(r·J)r = −∫r^{2p+2}J/(2p+2), which the script checks for p = 0–3. Since Dubovik and Radescu & Vaman define the radii with two terms, each static radius is also a single integral. We found no source that writes this down for general p, but it is a corollary, not a result.

### Gaps
- Dubovik & Tugushev 1990 (paywalled), Radescu & Vaman 2002 (paywalled) and Papasimakis et al. 2016 (403 at Southampton) were not opened. Their content is known here only through Szmytkowski & Łukasik, Borovkov & Canós Valero and Krynytskyi & Rovenchak.
- Gray, Karl & Novikov, Am. J. Phys. 78, 936 (2010), "Magnetic multipolar contact fields: the anapole and related moments", is the most likely textbook-level source for general-order static anapole moments. It was not found in open text.
- Google Scholar and NASA ADS were not queried directly: ADS needs an API token and Scholar blocks automated access. Every search here went through a general web search engine (§5).

## 2. Is there a clean general theorem? Statement, proof and conversion factors

### Takeaway
Yes, but its core is textbook: a closed current's Fourier transform is transverse. The statement worth recording ties every published convention together:
- every toroidal-family moment of a closed current is one of two single integrals per channel;
- each published definition is a fixed rational multiple of it;
- the curl-ladder normalisation hₙ (current loop → toroidal solenoid → supertoroid) coincides exactly with the coefficients of the exact long-wavelength dipole expansion, at every order.

The explicit all-n closed forms are new to our search. The ingredients are not.

### Cited Findings
- **Alaee, Rockstuhl & Fernandez-Corbaton** (1701.00755, Opt. Commun. 407, 17, 2018), Table II: p_α = −(1/iω){∫J_α j₀(kr) + (k²/2)∫[3(r·J)r_α − r²J_α] j₂(kr)/(kr)²} and m_α = (3/2)∫(r×J)_α j₁(kr)/(kr). — `1701.00755.txt`
- **Borovkov & Canós Valero** (J. Phys. Conf. Ser. 2015, 012021, 2021):
  - From Alaee's p they give the n-th order electric toroidal dipole T⁽ⁿ⁾ ∝ a(n)∫(kr)^{2n}[(n+1)J − n(r·J)r/r²], with a(n) = (−1)ⁿ6(n+1)/(2n+3)! (eqs. 2–3). They use the normalisation p = p + (ik/c)T⁽¹⁾ + (ik³/c)T⁽²⁾ + … (eq. 5).
  - From Radescu & Vaman they give the n-th mean-square radius of the toroidal dipole as R⁽ⁿ⁾ = [1/(2c(n+1)(2n+5))]∫r^{2n}[(n+1)(r·J)r − (n+2)r²J] (eq. 10).
  - The relation is T⁽ⁿ⁺¹⁾ = i a(n) k^{2n+1} R⁽ⁿ⁾ (eq. 11): "equivalent up to a prefactor". They say the magnetic case extends directly. — `borovkov2021.txt`
- **Boardman & Marinov** (physics/0510155), the supertoroid ladder:
  - T₀ = m = ½∫r×dI (eq. 13);
  - T₁ = ½∫r×dm = NπdR²I/2 (eqs. 14, 16);
  - T₂ = ½∫r×dT₁ = π²D²NI(R₂² − R₁²)/16 (eqs. 17, 19);
  - "and so on". — `physics_0510155.txt` §3
- **Nemkov, Basharin & Fedotov 2018**: a third-order supertoroidal coil "represents a source of the 1st order MSR of a magnetic dipole" (Fig. 4c), and a fourth-order coil represents the first MSR of a toroidal dipole (Fig. 5c). Two nested poloidal tori with zero toroidal dipole lose their "usual multipoles of the same order" when R² − ρ²/4 = R₀² − ρ₀²/4 (eq. 13). — `dl_toroidal/1804.07342.txt` §III
- **Gurvitz et al.** (1810.04945, Laser Photonics Rev. 13, 1800266, 2019), Table 2 (read from the PDF page image):
  - T⁽ᵉ⁾ = (1/10)∫[(J·r)r − 2r²J];
  - T⁽²ᵉ⁾ = (1/10)∫[3r⁴J − 2r²(r·J)r];
  - T⁽Qᵉ⁾ (1/42, as Vrejoiu);
  - T⁽Oᵉ⁾ = (1/300)∫[35(r·J)rrr − 20r²(Jrr + …) + (δδ + …)((r·J)r²r + 4Jr⁴)];
  - T⁽ᵐ⁾ = −(iω/20)∫r²(r×J).

  Their eq. 23 has the combination p + (ik/c)T⁽ᵉ⁾ + (ik³/c)T⁽²ᵉ⁾. — `dl_toroidal/1810.04945.pdf` pp. 13, 15

### Inferences (statements and proofs; status of each)

**Lemma 0 (master identity; known, textbook).** For a closed J and any smooth f of polynomial growth, ∫J·∇f dV = ∫∇·(fJ) − ∫f∇·J = 0. With f = e^{−ik·r} this is k·J̃(k) = 0: the form factor is transverse. Every identity below is a choice of f.

**Theorem A ((b) and its radii; known for p = 0, trivial for p ≥ 1).** f = r^{2p+2}xⱼ gives ∫r^{2p}(r·J)r = −∫r^{2p+2}J/(2p+2). At p = 0 this gives T = ½∫(r·J)r = −¼∫r²J = (1/6)∫r×(r×J) = ½∫r×M for any localised M with ∇×M = J. The last step uses ∫r²∇×M = −2∫r×M. Any α∫(r·J)r + β∫r²J with β = α/2 − ¼ equals T; Dubovik (1/10, −1/5) and Haxton (1/6, −1/6) both satisfy this. Planes' ¼∫r×(r×J) has (α, β) = (¼, −¼) and does not. The scalar case f = r^{2p+2} gives ∫r^{2p}(r·J) = 0.

**Theorem B (one-term reduction at every order; known in spherical form, Cartesian statement new to our search).** f = r^{2p+2}h, with h a harmonic polynomial of degree L, gives
∫r^{2p+2} J·∇h = −(2p+2)∫r^{2p}(r·J) h.
In STF Cartesian form, with τ_{L,p} = STF∫r^{2p}(r·J)x^{⊗L} and σ_{L,p} = STF∫r^{2p+2}J⊗x^{⊗(L−1)}, this is σ_{L,p} = −(2p+2)/L · τ_{L,p}. For static closed currents, all the published irreducible toroidal 2^L-poles then reduce to
T^{(L)} = τ_{L,0}/(L+1).
This holds for Dubovik–Tugushev (via Szmytkowski–Łukasik C.13), Vrejoiu/Gurvitz (1/10, 1/42, 1/300 at L = 1, 2, 3) and Krynytskyi–Rovenchak (every n). For Krynytskyi–Rovenchak, contract with a null a^{⊗n}:
[1/(2(n+1)(2n+3))]·[2n + 2(n+3)]·∫(r·J)(a·r)ⁿ = ∫(r·J)(a·r)ⁿ/(n+1).
This gives, for example:
- the static toroidal quadrupole Tⱼₖ = ⅓∫(r·J)xⱼxₖ = −⅙∫r²(xⱼJₖ + xₖJⱼ);
- the static toroidal octupole = ¼ STF∫(r·J)xxx.

**Theorem C (counting; known for N ≤ 2, standard angular-momentum algebra).** Let N be the number of position factors. The order-N moments ∫J⊗x^{⊗N} span Sym^N ⊗ V = ⊕_l (D_{l+1} ⊕ D_l ⊕ D_{l−1}), with l = N, N−2, … Lemma 0 with f = x_{i₀}…x_{i_N} kills the totally symmetric part, Sym^{N+1} = D_{N+1} ⊕ D_{N−1} ⊕ … What survives:
- the axial channels L = N, N−2, … ≥ 1 once each (the magnetic 2^L-poles and their radii);
- the polar channels L = N−1, N−3, … ≥ 1 once each (the toroidal 2^L-poles and their radii, which are the τ_{L,p} with L + 2p + 1 = N).

In total that is N(N+2) components.
- Corollary: each channel has multiplicity 1, so any two definitions of a toroidal moment or radius in that channel differ by a scalar factor for closed currents. This is why every convention conversion is a single number.
- The script checks the rank on 120 random closed polygons: 3, 8, 15, 24 for N = 1–4. Open segments give 9, 18, 30, 45.

**Theorem D (the ladder; new to our search as explicit closed forms, ingredients known).** Let J = (∇×)ⁿK with K localised, and set hₙ = ∫K. K is unique up to ∇χ with χ localised, which does not change ∫K. Then:
- h_{2k+1} = (−1)^k · 3/(2^{k+1} k! (2k+3)!!) · ∫r^{2k}(r×J);
- h_{2k+2} = (−1)^{k+1} · (3/2)/(2k+3)! · ∫r^{2k+2}J.

That is:
- h₁ = ½∫r×J = m;
- h₂ = −¼∫r²J = T;
- h₃ = −(1/20)∫r²(r×J);
- h₄ = (1/80)∫r⁴J;
- h₅ = (1/560)∫r⁴(r×J);
- h₆ = −(1/3360)∫r⁶J.

Proof: two integrations by parts. The first, (B), is ∫r^{2p+2}∇×Y = −(2p+2)∫r^{2p}(r×Y). The second, (A), is ∫r×(∇×Y) = 2∫Y, and for p ≥ 1 and solenoidal Y, ∫r^{2p}(r×∇×Y) = (2p+3)∫r^{2p}Y, which uses Theorem A on Y. Apply them alternately down the chain X_{j−1} = ∇×X_j, from X_n = K to X_0 = J.

Consequences:
- **Ladder = exact-expansion coefficients.** For closed J, Alaee's dipole brackets are generating functions of the ladder:
  - P(k) := ∫J j₀ + (k²/2)∫[3(r·J)r − r²J] j₂/(kr)² = (3/2)∫J j₀(kr) = (3/2)k∫(r·J) r̂ j₁(kr) = Σ_{m≥1} k^{2m} h_{2m};
  - M(k) := (3/2)∫(r×J) j₁(kr)/(kr) = Σ_{m≥0} k^{2m} h_{2m+1}.

  The P(k) closed form is Lemma 0 with f = j₀(kr)xⱼ; the series follow by Taylor-expanding j₀ and j₁/x and applying Theorem A. So the supertoroid ladder ½∫r×d(lower) of Boardman & Marinov and Afanasiev, Borovkov's T⁽ⁿ⁾ (in the eq. 5 normalisation) and Gurvitz's T⁽ᵉ⁾ and T⁽²ᵉ⁾ (with 1/280) are the same numbers for static closed currents: Boardman–Marinov T_{n−1} = hₙ. Radescu–Vaman's static radii are R⁽ⁿ⁾ = −∫r^{2n+2}J/(4(n+1)), and h_{2n+2} = a(n)R⁽ⁿ⁾, which is the static content of Borovkov eq. 11.
- **The iterated-density convention.** Nₙ = 2⁻ⁿ∫(r×)ⁿJ is what Planes et al. use for the toroidal moment and Prosandeev & Bellaiche (as quoted by Luk'yanchuk et al.) for the hypertoroidal moment. It is not the ladder. Nₙ/hₙ:
  - for n = 2k+1: k!(2k+3)!!/(3·2^k);
  - for n = 2k+2: (2k+3)(2k+3)!/(3·2^{2k+1}(2k+2)).

  The ratios are 1, 3/2, 5/2, 25/4, 35/2 and 245/4 for n = 1–6. Planes' 3/2 is the n = 2 entry.
- **The crate.** `magnetics::toroidal::hypertoroidal_moment` (½Σcᵢ×vᵢ) is the ladder hₙ. For a rig of 16 poloidal tori on a ring, −(1/20)∫r²(r×J) equals ½Σcᵢ×Tᵢ to 1.6×10⁻¹⁴. That is exact for these sub-tori, which have no local m, magnetic quadrupole or magnetic radius.

**Conversion table (static, closed J; S_m := ∫r^{2m}J and A_m := ∫r^{2m}(r×J)).**

| Source | Quantity as published | Static value | Relation to the ladder |
|---|---|---|---|
| Dubovik–Tugushev; Ederer–Spaldin; Gurvitz T⁽ᵉ⁾; Krynytskyi–Rovenchak | (1/10)∫[(r·J)r − 2r²J] | −¼S₁ | h₂ |
| Haxton et al. 2002 (LH units) | a/M² = (1/6)∫x×(x×j) | −¼S₁ | h₂ |
| Dumitriu–Vrejoiu 2005 | (1/6)∫ξ×(ξ×j) − (1/30)Ṗ_qqi | −¼S₁ | h₂ |
| Planes et al. 2015, eqs. 29, 33 | ¼∫r×(r×J) | −⅜S₁ | (3/2)h₂ |
| Boardman–Marinov 2005 | T₁ = ½∫r×dm; T₂ = ½∫r×dT₁ | — | h₂; h₃ |
| Gurvitz T⁽ᵐ⁾ (magnitude) | (ω/20)∫r²(r×J) | — | \|h₃\| up to the iω/c bookkeeping |
| Dubovik–Cheshkov 1974 | ρ̄² = ½∫ξ²(ξ×J) | ½A₁ | −10h₃ |
| Gurvitz T⁽²ᵉ⁾ with 1/280; Borovkov T⁽²⁾ | (1/280)∫[3r⁴J − 2r²(r·J)r] | (1/80)S₂ | h₄ |
| Radescu–Vaman via Borovkov eq. 10 | R⁽ⁿ⁾ | −S_{n+1}/(4(n+1)) | h_{2n+2}/a(n) |
| Prosandeev–Bellaiche (via Luk'yanchuk) | T = ½∫r×P; G = ¼∫r×(r×P) | — | h₁[P]; (3/2)h₂[P] for ∇·P = 0 |
| Wikipedia (spins) | T = Σrᵢ×Sᵢ | — | 2× Ederer–Spaldin's ½Σr×m, with S in place of m |

### Gaps
- Theorem D's closed forms and ratio table were not found in any opened source. Borovkov & Canós Valero, Radescu & Vaman (not opened) or Gray–Karl–Novikov (not opened) may contain equivalent formulas in another basis. Treat it as **new to our search**, not as new.
- Jackson's exact multipole coefficient (Classical Electrodynamics, 3rd ed., §9.10) is cited from memory and was not re-opened. Góngora & Ley-Koo's eq. 42 is the opened source for the ρ = 0 one-term Bessel weight.

## 3. Planes et al.'s 3/2, and other normalisation slips

### Takeaway
The 3/2 is a genuine internal inconsistency in Planes, Castán & Saxena. Their eq. 28 is the toroidal term of the vector-potential expansion of their eq. 23. It holds with t = Dubovik's T, but their eq. 29 gives t = ¼∫r×(r×J) = (3/2)T. Eq. 30 arises from treating m(r) = ½r×J as a magnetisation, which it is not.

We found no erratum. The review cites no number computed with eq. 29 or 33, and its spin form (eq. 34, ½Σr×s) is the standard one, so no published value we found is affected.

The same ¼∫r×(r×·) appears as a definition in the ferroelectric hypertoroidal moment of Prosandeev & Bellaiche, as quoted by Luk'yanchuk et al. There it is a convention, internally consistent, and their published numbers reproduce. It differs from the ladder normalisation by exactly 3/2 when P is solenoidal.

A second slip: Gurvitz et al.'s Table 2 prints T⁽²ᵉ⁾ with 1/10. Their own eq. 23 combined with Alaee's exact dipole needs 1/280.

### Cited Findings
- **Planes, Castán & Saxena** (arXiv:1504.04958 v1, 20 Apr 2015; journal ref. "Multiferroic Materials 1 (2014) 9", DOI 10.2478/muma-2014-0002, which resolves to a Portico archive record). The equations as extracted:
  - eq. 23: A(R) = (µ₀/4π)Σ((−1)ⁿ/n!)∫J (r·∇)ⁿ(1/R);
  - eq. 28: A⁽²⁾_tor = ∇(t·∇)(1/R) + tδ(R);
  - eq. 29: t = ¼∫r×[r×J];
  - eq. 30: "Defining m(r) = ½[r×J(r)] as the distribution of magnetic moments", t = ½∫r×m;
  - eq. 33: point charges, ¼Σq r×(r×u);
  - eq. 34: ½Σr_α×s_α times a spin prefactor.

  Ref. [24] for the expansion is Raab & de Lange's book. Only one arXiv version exists. — `dl_toroidal/ferrotoroidic_review.txt` §3; arXiv abstract page
- **Ederer & Spaldin 2007** use ½∫r×µ with ∇×µ = j/c, the correct counterpart of eq. 30. — `0706.1974.txt` eqs. 1, 5, 7
- **Luk'yanchuk et al.**, "Topological ferroelectric chirality" (arXiv:2406.19728): "toroidal, T, and hypertoroidal, G, moments (Prosandeev and Bellaiche, 2008)" are defined by densities ½r×(P − P̄) and ¼r×r×(P − P̄) (eqs. 4–5). For the model Hopfion (eqs. 1–2, a divergenceless P in a sphere) they report T = 0.96 αₜR⁴e_z and G = 1.19 α_pR⁵e_z. — `2406.19728.txt`
- The Prosandeev–Bellaiche 2008 APS abstract (search snippet) describes "a double cross product of the local dipoles with their positions". The primary PRB 77, 060101(R) was not opened.

### Inferences
- The script computes the Coulomb-gauge A = ∫J/|R−r| of a 24-turn winding (m = 0, no magnetic quadrupole) at R = 50 and 100. It fits ∇(t·∇)(1/R) with t/T = 0.999829 and 0.999957, which approaches 1, not 1.5. So eq. 29 is inconsistent with eq. 28 by 3/2.
- The full second-order term is (µ₀/4π)[∇(T·∇)(1/R) + 4πTδ(R)], from ∫xᵢxₗJⱼ = δᵢⱼTₗ + δₗⱼTᵢ − 2δᵢₗTⱼ for closed J. Planes' "tδ(R)" therefore appears to lack the 4π. This is **uncertain**: it rests on text extraction, and the prefactor may sit outside the extracted line.
- Published numbers: none found that depend on eq. 29 or 33.
  - The review itself computes nothing with them.
  - Spin-based (eq. 34; Ederer & Spaldin) and dipole-based (g = ½Σr×p) toroidal moments use the ½r× form, which is unaffected.
  - The ¼ double cross product survives in ferroelectrics as the hypertoroidal *definition*. The script reproduces Luk'yanchuk's T = 0.9648 (published 0.96) and G = 1.1894 (published 1.19), and finds G/(−¼∫r²P) = 1.500005. In the ladder normalisation the Hopfion's moment would be 1.1894/1.5 ≈ 0.793 α_pR⁵. This is a choice of convention, not an error, but anyone comparing it with Dubovik-type moments needs the 3/2.
  - A citation search (who cites Planes eq. 29 and computes with it) was not possible with these tools.
- **Gurvitz T⁽²ᵉ⁾.** Expanding Alaee's p for a closed knot current, the script finds the k² coefficient equals Dubovik's T, and the k⁴ coefficient equals (1/280)∫[3r⁴J − 2r²(r·J)r] to 1.1×10⁻⁶ after Richardson extrapolation. Borovkov's static bracket is (n + 3/2)S_n; with a(2) = 18/7! and 7/2 it gives the same value, (1/80)S₂. The 1/10 in Gurvitz's Table 2 is therefore **likely a misprint** in the arXiv version. That is **uncertain** until the published LPR version is checked; the integrand itself is right.
- **Erratum or comment?** The 3/2 is worth a short note to the authors or a comment. The journal appears to be archived, so a correction may not be possible there. It does not warrant a paper.

### Gaps
- The primary Prosandeev & Bellaiche PRB 77, 060101(R) (2008) and J. Mater. Sci. 44, 5235 (2009) were not opened. Their normalisation (¼ or ½ per level, and any 1/N per cell) is known here only through Luk'yanchuk et al.
- Papasimakis et al. 2016 was not opened (403), so its hypertoroidal/higher-order normalisation is unknown.
- The published (LPR) version of Gurvitz et al. was not checked.

## 4. Testable consequences: the numerical checks

### Takeaway
Every identity in §2 holds to rounding on non-trivial closed currents: a tilted, shifted toroidal winding; a trefoil knot; 120 random closed polygons; and a random divergence-free field built as (∇×)ⁿ of random vector Gaussians, for n = 1–6. The checks also reproduce two published results from the one-term forms: Luk'yanchuk's Hopfion values and Nemkov's eq. 13. The identities need closure: the rank test shows that open currents carry the extra components which closure removes.

### Cited Findings
- Script `notes/multipole_identities_checks.py` (numpy; `python3 -P`; 61 assertions; exit 0). Polygon moments use 8-point Gauss–Legendre per segment, exact for polynomial weights up to degree 15. Curl fields use FFT derivatives on a 96³ grid in a ±6 box.

### Inferences (results as printed by the script)
- **(b) and Planes:** ½∫(r·J)r = −¼∫r²J = (1/6)∫r×(r×J) = T to ≤5.3×10⁻¹⁵ on the torus and knot. Planes/T = 1.500000000000 on both, and on the curl field to <10⁻⁹. On the curl field ½∫r×M = T as well.
- **Radii and scalars:** ∫r^{2p}(r·J)r = −∫r^{2p+2}J/(2p+2) and ∫r^{2p}(r·J) = 0 for p = 0–3.
- **One-term reductions:**
  - harmonic weights for L = 1–5, p = 0–2;
  - Gurvitz/Vrejoiu T⁽Qᵉ⁾ and T⁽Oᵉ⁾ are traceless for any J, and their static values are ⅓∫(r·J)xx and ¼STF∫(r·J)xxx;
  - Krynytskyi–Rovenchak's T⁽ⁿ⁾ matches (1/(n+1))∫(r·J)xⁿ for n = 1–6.
- **Counting:** ranks 3, 8, 15 and 24 (closed) against 9, 18, 30 and 45 (open).
- **Ladder:**
  - hₙ = ∫K for n = 1–6 to between 1.2×10⁻¹³ and 5.5×10⁻¹⁰;
  - the iterated-r× ratios measure 1.000000000, 1.500000000, 2.500000000, 6.250000000, 17.499999999 and 61.249999987;
  - the 16-torus rig gives h₃ = ½Σcᵢ×Tᵢ to 1.6×10⁻¹⁴.
- **Generating function:** Alaee's bracket = (3/2)∫J j₀(kr) = (3/2)k∫(r·J)r̂ j₁(kr) at k = 0.7, 1.3 and 2.9, to ≤9.4×10⁻¹⁵.
- **Nemkov eq. 13 from one line:**
  - For a poloidal winding (R, a), r·dl = −aR sin θ dθ. This gives T ∝ a²R and a static toroidal octupole ∝ a²R(R² − a²/4).
  - With T = 0, the octupole of two nested tori therefore cancels iff R₁² − a₁²/4 = R₂² − a₂²/4, which is Nemkov's eq. 13.
  - The script finds the octupole at 3.2×10⁻⁵ of a single torus's (the polygon error) when the radii match, against 0.222 when the outer major radius is 0.1 larger.
  - What remains is the first radius, τ_{1,1,z} = ∫r²(r·J)z = −πN₁I₁a₁²R₁·(5/4)(a₁² − a₂²) by the same algebra (not separately asserted).

### Gaps
- Nemkov's eq. 12 (the magnetic case) was not tested, because its surface-current prescription is ambiguous in the text.
- The trefoil and random polygons are generic closed currents, but not physically "knotted-field" currents, such as helicity-carrying volume fields. The curl-field test covers smooth volume currents.

## 5. Search log

Queries (WebSearch, general engine, 2026-10-27), with what each returned of use:
1. Gray Karl Novikov "Magnetic multipolar contact fields: the anapole and related moments" r^2 J toroidal: AJP 78, 936 citation only.
2. anapole moment static current "a = -π ∫ r^2 j" toroidal dipole equivalence closed current: Haxton et al. nucl-th/0109014.
3. Vrejoiu "Electromagnetic multipoles in Cartesian coordinates" toroidal moments static currents: physics/0310064, 0503164, 0503228, 0511088.
4. Radescu Vaman "toroid moments" arXiv exact calculation: Vrejoiu & Zus 0906.0503 and 0910.1313; Radescu–Vaman not on arXiv.
5. Prosandeev Bellaiche hypertoroidal moment definition (two queries): the APS 2008 abstract; Luk'yanchuk 2406.19728.
6. Spaldin Fiebig Mostovoy toroidal moment definition: Ederer & Spaldin 0706.1974.
7. Gongora Ley-Koo "Complete electromagnetic multipole expansion including toroidal moments": the RMF-E PDF.
8. Castellanos Panizo Rivas "Magnetostatic multipoles in Cartesian coordinates": record only.
9. static toroidal multipole moments closed currents reduce r^2 J (r·J) higher order: Szmytkowski & Łukasik 1602.00613; Agre 2011 (UFN).
10. Papasimakis … eprints: Southampton PDF (403).
11. Alaee Rockstuhl Fernandez-Corbaton exact multipole arXiv: 1701.00755.
12. Castán Planes Saxena toroidocaloric toroidal moment definition: 1205.1392 (no current-based definition).
13. "toroidal moment" "r × (r × j)" 1/6 OR 1/4 erratum Planes: nothing found.
14. hypertoroidal moments curl hierarchy mean-square radius equivalence: Boardman & Marinov physics/0510155; Borovkov & Canós Valero 2021.
15. exact multipole expansion closed current toroidal form factor Taylor j_l(kr) r·J: Krynytskyi & Rovenchak 1804.10494; Nemkov 2201.08333; Fernandez-Corbaton 1507.00752.

Also fetched: arXiv abstract page of 1504.04958 (one version, journal ref.); Wikipedia "Toroidal moment" (formula and spin convention); Dubovik & Cheshkov 1974 (JINR archive scan, 10 pages); IOP open-access PDF of Borovkov & Canós Valero.

Not available to these tools: Google Scholar and NASA ADS (no API token), Dubovik & Tugushev 1990, Radescu & Vaman 2002, Papasimakis 2016, Gray–Karl–Novikov 2010, Prosandeev & Bellaiche 2008/2009 primaries.

## 6. Novelty verdict

- **Paper-worthy result:** none. The general one-term theorem is published (Dubovik & Tugushev 1990 via Szmytkowski & Łukasik 2016 App. C; Agre 2011). The MSR ↔ high-order-toroidal equivalence is published (Borovkov & Canós Valero 2021), and so is the supertoroid ladder (Boardman & Marinov 2005; Nemkov et al. 2018).
- **Erratum or comment:**
  - Planes et al.'s ¼∫r×(r×J) = (3/2)T, inconsistent with their own eq. 28, with no published numbers found to depend on it;
  - Gurvitz et al.'s T⁽²ᵉ⁾ prefactor, 1/10 where 1/280 is needed (uncertain until the journal version is checked).
- **Textbook or pedagogical:** a unified note stating Theorems A–D with the conversion table and the ratio sequence Nₙ/hₙ = 1, 3/2, 5/2, 25/4, 35/2, 245/4, … would be a correct short Eur. J. Phys. or Am. J. Phys.-style note. Its only new-to-our-search content is the explicit closed forms for hₙ, the exact equality "supertoroid ladder = exact-dipole Taylor coefficients", and the closed form (3/2)∫J j₀(kr). All three are short corollaries.
- **Strongest honest claim:** "For closed currents, every published toroidal, hypertoroidal and mean-square-radius moment is a rational multiple of one of the single integrals ∫r^{2k}(r×J) or STF∫r^{2p}(r·J)x^{⊗L}. We give the multiples for eleven conventions, show that the supertoroid ladder ½∫r×d(lower) is exactly the Taylor series of Alaee's exact dipoles, and identify two normalisation misprints (Planes 2015, a factor 3/2; Gurvitz 2019, 1/10 for 1/280), all verified numerically."
