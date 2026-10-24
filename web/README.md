# web

3D models of the physics in the crates, each one a single HTML file that runs in the browser.
Open the file directly; there is nothing to build and no server to run. three.js loads from a CDN,
so the first load needs a connection.

| Model | What it shows |
|---|---|
| [`gravity-patent.html`](gravity-patent.html) | US 12,417,870 rebuilt from its embodiment: a 0.35 g polyethylene ball on a cotton thread in a 10 cm vacuum tank, between two open coils, the upper one live. The coils' electric field is solved as 24 charged rings and drawn as field lines and a field slice. A panel gives every force on the ball beside its weight: the pull toward the live coil, the charge any glow would give it, and the polarisation torque against the thread. It also shows whether the gas in the tank breaks down at the pressure you set, the corona at the bare copper ends, the magnetic field of the charging current, the claim in its own terms, and general relativity's prediction. Switch the drive to ±V/2 or the ball to PTFE to run the tests that would settle it. |

The equations are the crates', ported line for line: `gravitics::electrodes` (the ring solver),
`gravitics::electrogravity` and `plasmoid::discharge`. At 50 kV, every number matches
`crates/physica/examples/gravity_patent.rs`, which `cargo test` asserts. The research behind the
model, and what it assumes beyond the patent, is in
[the review of Zhang's sources](../research/external/2026-10-24-zhang-unified-field-theory-review.md).
