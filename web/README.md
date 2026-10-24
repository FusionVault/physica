# web

3D models of the physics in the crates, each one a single HTML file that runs in the browser.
Open the file directly; there is nothing to build and no server to run. three.js loads from a CDN,
so the first load needs a connection.

| Model | What it shows |
|---|---|
| [`gravity-patent.html`](gravity-patent.html) | US 12,417,870 rebuilt from its embodiment and simulated: a 0.35 g polyethylene ball on a cotton thread in a 10 cm vacuum tank, between two open coils, the upper one live. The coils and leads are solved as 93 charged rings and wire segments, drawn as field lines and a field slice. The ball swings and turns under every force and torque the crates compute: its charge and charge patches, polarisation, the gas (drag, breakdown, a glow's clamping and draining), magnetism, the field's own gravity, and, when chosen, Zhang's vortex field with each of his three constants. Set the pressure, the drive, the ball, its charge and roundness, and the thread; Reset returns to the patent's settings. |

The physics is the crates', ported line for line: `gravitics::electrodes`,
`gravitics::electrogravity`, `plasmoid::discharge` and `magnetics`, put together exactly as
`crates/physica/examples/gravity_patent.rs` does. On load the page reruns the example's pinned
reference values (the field, the charges, and a 1 s trajectory) and shows whether they match;
`cargo test` pins the same values on the Rust side. The research behind the
model, and what it assumes beyond the patent, is in
[the review of Zhang's sources](../research/external/2026-10-24-zhang-unified-field-theory-review.md).
