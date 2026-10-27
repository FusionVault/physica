# 2026-10-27: OpenAI's "Velocity-Field Detection of Computation in Forced Navier–Stokes Flows" (link)

**Source:** <https://github.com/openai/math/blob/main/preprints/Velocity-Field-Detection-of-Computation-in-Forced-Navier-Stokes-Flows-September-27-2026/manuscript.pdf>

- OpenAI, 27 September 2026; 18 pages.
- Licensed Apache-2.0 (the `openai/math` repository).
- A preprint, not peer reviewed.
- The paper itself is not stored here; this is the link and a summary in our words.

**The user's request (verbatim):**

> OpenAI just recently released a bunch of preprints including navier stokes. Now, i want to implement navier stokes and all relating equations in the paper into our repo stack please (goes in whatever related crates). There are ALOT of other mathematical papers here too, which i want to go over eventually but i think navier stokes is the most important since it applies to magnetohydrodynamics. […] As mentioned please file an issue for it please specific to just this paper

**Summary.**

- **Theorem 1.1:** for 3-D incompressible Navier–Stokes started from rest (eq. 1), at any computable viscosity, there are smooth forces whose velocity field detects whether a Turing machine halts.
- **Two constructions:**
  - On the flat torus T³: a pointwise test, u₃ > ½ somewhere in 1/32 < y < 1/8. It uses short stirring bursts and long heat-only waits.
  - On R² × T: an integral test, ∫_{X₂>0} u₃ > ½. It uses an expanding array of translation gates, with every mixed derivative of the force bounded.
- **The common device** (eqs. 5–6): for a planar divergence-free field a(t, Y), U = (a, w) with pressure zero solves Navier–Stokes. This needs w to solve an advection–diffusion equation, and the horizontal force to be F = ∂ₜa + (a·∇)a − ν∆a.
- **The diffusion estimates** carry the proofs:
  - torus heat-kernel bounds (Lemma 3.3);
  - a mass bound below 10⁻¹¹ (eq. 13);
  - burst errors below 1/50 (eq. 18);
  - a moving-cutoff loss budget below 1/24 (eq. 39).
- **Corollary 6.1:** neither velocity event is decidable.
- **Related:** the companion preprints [10–12] (a fixed-particle test, slow clocks, eventually stationary forcing) supply the planar processor. They are separate papers.

**Review status:** tracked in [#53](https://github.com/FusionVault/physica/issues/53).
