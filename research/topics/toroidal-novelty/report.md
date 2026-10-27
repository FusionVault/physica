# Toroidal moments hold one publishable result, about pulses, and two errata

This report asks whether toroidal moments, anapoles and the Aharonov–Bohm effect hold anything genuinely new that can be proved and written up. It was investigated after `magnetics::toroidal` was built (#23), along four independent lines. Each line searched the literature, proved every candidate, and checked it numerically with a gating script. All four scripts pass.

**The answer is mostly no.** Three of the four lines found only known physics. The one candidate that could stand as a short paper is a closed-form energy and momentum for the supertoroidal pulse family at every α. That claim holds only if three sources we could not open do not already contain it. The investigation also produced two corrections to published papers and four crate fixes.

## Verdicts

| Line | Notes | Checks | Verdict |
|---|---|---|---|
| Closed-current multipole identities | [multipole_identities.md](notes/multipole_identities.md) | 61, `-P` | **Known.** Szmytkowski & Łukasik 2016 prove that every toroidal-family moment of a closed current reduces to a single moment of J, at every order. They attribute it to Dubovik & Tugushev 1990. New to our search, but elementary: closed forms for every rung of the hypertoroidal ladder, the fact that the ladder is the Taylor series of the exact dipoles, and conversion factors between 11 conventions. |
| Gravitational analogues | [gravitational_toroidal.md](notes/gravitational_toroidal.md) | `-P` | **Mostly known.** The exterior of a toroidal mass current is known (Lynden-Bell & Katz), and so is the toroidal term in the mass moments (Blanchet eq. 138). New to our search, but quick consequences: a non-radiating, field-free "gravitational anapole"; the fact that a single closed vortex with longitudinal stress radiates nothing; that the stresses count as much as the mass current at quadrupole order; and a closed-form gravitomagnetic AB phase round a torus, exactly 4× the electromagnetic analogue. That phase is 2×10⁻¹² rad for a tonne of rings at 1 km/s, about 10¹⁸ interferometer shots. Damour & Iyer 1991 must be read first. |
| Monopoles and duality | [monopole_toroidal.md](notes/monopole_toroidal.md) | 60, `-I` | **Known.** The four multipole families with magnetic currents, the duality-covariant anapole, the Dirac string as an h/e line, and the dual AB phase are all published. New to our search, but small: a magnetic charge carried round a loop of vector area 𝒜 leaves a toroidal moment Q𝒜. In spin ice that is 200.3 µ_B·Å per hexagon flip. |
| Observability | [observability.md](notes/observability.md) | 91, `-P` | **The debate is settled** by Fernandez-Corbaton et al. 2017: at one frequency, no outside measurement separates p from T. New to our search, as a reconciliation: T is identifiable exactly when the current pattern is rigid and the drive known, and a pattern deforming at O(k²) mimics T exactly. **The strongest candidate:** the supertoroidal pulse's energy U = (π^{3/2}/2)μ₀f₀²[Γ(α+½)/Γ(α)](αq₁ + q₂)/(q₁³q₂^{2α+1}) and momentum cP_z/U = (q₂ − αq₁)/(q₂ + αq₁). Both are finite for every α > 0, at rest when q₂ = αq₁, covariant under Lorentz boosts, and match numerical integration to 10⁻¹²; at α = 1 they reproduce Lekner. |

## What could be written

1. **A short paper or Comment on supertoroidal pulses** (Shen et al., Nat. Commun. 12, 5891, 2021). It would contain:
   - the closed-form energy and momentum for every α;
   - the finding that "α ≥ 1 for finite energy" is unnecessary;
   - the finding that their printed H_z (eq. 6 and SI eq. C3) is wrong for α ≠ 1, since it fails ∇·H = 0, while their Fig. 6 matches the correct field;
   - the finding that the whole-plane polarity of their magnetic "skyrmions" is ±½ for α ∉ {1, 2}.

   **Before writing:** read Ziolkowski 1989, Hellwarth & Nouchi 1996 and Lekner's 2018 book; any of them may already give the general-α energy.
2. **Errata notes** to the authors of two papers:
   - Planes, Castán & Saxena 2015: their eq. 29/30 is 3/2 of the toroidal moment their own eq. 28 implies, because m = ½r×J is not a magnetisation.
   - Gurvitz et al. 2019: Table 2's 1/10 for the second-order toroidal term should be 1/280. This is uncertain until the journal version is read.
3. **A teaching note** collecting the conversion factors between toroidal-moment conventions and the closed forms of the hypertoroidal ladder (Am. J. Phys. or Eur. J. Phys.). It is worthwhile only if Dubovik & Tugushev 1990 does not already contain it.

Nothing here supports a research paper on new physics of toroidal moments or monopoles. The talk that started this line of research claimed links between the two (see [the topological-monopole report](../topological-monopoles-and-transmutation/report.md)); this investigation finds the real physics nearby, and it is known.

## What went into the crates

- `magnetics::monopole::spin_ice_charge`: `a` is the nearest-neighbour spacing (3.54 Å), not the cubic constant (10.01 Å). The formula and tests were already right.
- `magnetics::toroidal`: the supertoroidal energy is finite for every α > 0. `effective_dipole` now says that only p + (ik/c)T is observable from outside. `toroidal_radiated_power`'s ω⁶ law now says it assumes a rigid current pattern.
- **Candidates for later, each with a test value from the scripts:**
  - `supertoroidal_energy` and `supertoroidal_momentum_ratio`;
  - `supertoroidal_te_fields`, with the correct H_z;
  - `toroidal_interaction_energy`, μ₀T·J;
  - `gravitomagnetic_ab_phase`, with a `physica` analogy asserting the factor of 4;
  - `dual_aharonov_bohm_phase` and `spin_ice_string_phase`.

## Sources still to read

- Dubovik & Tugushev 1990
- Papasimakis et al. 2016
- Damour & Iyer 1991
- Ziolkowski 1989
- Hellwarth & Nouchi 1996
- Lekner 2018
- Nanz 2016
- Radescu & Vaman
- Afanasiev 1993
- the published version of Gurvitz et al. 2019

Most of these are paywalled. Each line's notes record its search queries, so the novelty claims can be repeated.
