# Animated visualisation of gravity and water-storage fields

Scope note: research done 24 October 2026. Library versions come from live npm registry, cdnjs API and crates.io queries made that day. Gzip sizes were measured by downloading the jsDelivr file and running `gzip -9`. The shared web-search budget ran out part-way through, so some planned lookups were not made. Each section's Gaps list says which. Facts about our own project are cited to files in this repo: `viewer/index.html`, `viewer/data/tws_index.json`, `crates/*`.

Project facts these notes are tailored to:
- The viewer is one 505-line `viewer/index.html`. It has CSS light/dark tokens, a month slider, play via `setInterval`, scrub via `requestAnimationFrame`, an orthographic globe and Equal Earth ([viewer/index.html](../../viewer/index.html)).
- Data: COST-G RL02.1 equivalent water height in mm, relative to the 2004–2009 mean. Degree 60, destriped, 300 km Gaussian, GIA not removed. A 2° grid of 180 × 90 cells, 259 months from 2002-04 to 2026-07. 35 months are missing across 19 gaps, the largest being 11 months (2017-06 → 2018-06, the GRACE→GRACE-FO gap). The grids ship as a 5.7 MB zlib+base64 JSON ([viewer/data/tws_index.json](../../viewer/data/tws_index.json); month gaps counted from that file).
- Crates:
  - `geopotential`: SH synthesis to degree 2700, grids of geoid/gravity/water height, degree variances and amplitudes.
  - `gravimetry`: Love numbers, Gaussian smoothing, destriping, anomalies.
  - `gravity-field`: ICGEM/TN-13/TN-14 readers.
  - `satgrav`: noise ASDs, `orbit_period`, `degree_at_frequency`, `sim::Pair::positions` (Earth-fixed positions of both GRACE-type satellites at time t), and a closed-loop `simulate` least-squares solution.

  Sources: `crates/*/Cargo.toml` and [crates/satgrav/src/sim.rs](../../crates/satgrav/src/sim.rs).

## Notable gravity/water animations (NASA SVS, JPL, GFZ GravIS, ESA GOCE, Earth Observatory): what makes them effective?

### Takeaway
The reference animations are simple. Each shows one variable with a two-hue brown/orange-to-green/blue "dry–wet" diverging key, a fixed projection, a visible date and a fixed frame rate (30 fps movies). They are effective through restraint and narrative framing, and through physical metaphor in the geoid "potato" (the factor-7000 exaggeration is always stated). None of the pages fetched documents how they interpolate between months or how they show uncertainty. That is the gap our viewer can fill honestly.

### Cited Findings
- NASA SVS #4415, "Terrestrial Water Storage Anomaly 2002–2015", uses GRACE TWSA from January 2002 to September 2015. It has a Robinson-projection frame set plus composite and date-overlay frame sets, and 4,993 frames at 1920×1080, shown as a 1080p 30 fps movie. Brown means less water than normal and green more, with a colourbar in cm. The page does not describe smoothing or interpolation — [NASA SVS 4415](https://svs.gsfc.nasa.gov/4415)
- The newer version, SVS #4476 (March 2015–March 2016), switches to orange for drier and blue for wetter — [NASA SVS 4415 page, "newer version" link](https://svs.gsfc.nasa.gov/4415); [SVS 4476](https://svs.gsfc.nasa.gov/4476)
- SVS has a separate California TWS animation for 2002–2017, with a newer version available — [SVS 30862](https://svs.gsfc.nasa.gov/30862)
- SVS #5574, "GRACE FO Soil Moisture Within Continental United States: Monitoring Drought":
  - Three separate percentile-based dryness-index layers: groundwater, root-zone soil moisture and surface soil moisture.
  - Updated weekly, covering roughly the year to the current week.
  - 3840×2160 movies at 30 fps, 1,000 frames per layer, made for NASA's Earth Information Center.
  - First published 2 March 2026, updated 23 October 2026. Scientists: M. Rodell and H. Kato Beaudoing.

  Source: [NASA SVS 5574](https://svs.gsfc.nasa.gov/5574/)
- JPL's "GRACE-FO: from range observations to global mass change" animation shows the twin satellites, with their separation changing as they pass over mass anomalies. It explains that the polar pair gathers enough data for a global map about every 30 days — [JPL GRACE resource 40](https://grace.jpl.nasa.gov/resources/40/grace-fo-from-range-observations-to-global-mass-change)
- JPL's "Himalaya plot" (PIA22507) draws the ground track with the inter-satellite range variation following the terrain below. It is an early, very effective "measurement → mass" explainer — [JPL PIA22507](https://www.jpl.nasa.gov/images/pia22507-recently-launched-twin-satellites-create-the-himalaya-plot/)
- GRACE-FO orbit: about 490 km initially, 475 km now, inclination 89.0° — [DLR GRACE-FO](https://www.dlr.de/en/rb/research-operation/missions/earth-observation-science/grace-follow-on)
- The JPL GRACE(-FO) Data Analysis Tool is for "quick-look comparisons". All data are interpolated to a 1×1° grid for computational reasons. The page does not describe its animation or colour features — [JPL DAT](https://grace.jpl.nasa.gov/data/data-analysis-tool/)
- GFZ GravIS (with AWI and TU Dresden) visualises Level-3 products:
  - Gridded mass anomalies and basin-average time series of TWS (non-glaciated land), ocean bottom pressure, and Antarctic/Greenland ice-mass change.
  - A prototype groundwater product from the EU G3P project.
  - Products update with each new monthly field.

  Sources: [GFZ GravIS page](https://www.gfz.de/en/section/global-geomonitoring-and-gravity-field/data-products-services/gravis-gravity-information-service); [gravis.gfz.de](http://gravis.gfz.de/home); [GravIS paper, ESSD 17:611, 2025](https://essd.copernicus.org/articles/17/611/2025/essd-17-611-2025.html)
- The GFZ "Potsdam Gravity Potato" visualises the EIGEN-6C combined model (LAGEOS, GRACE, GOCE, ground gravity and altimetry) — [Universe Today](https://www.universetoday.com/articles/the-potsdam-gravity-potato-shows-earths-gravity-variations)
- ESA's GOCE geoid video shows the geoid amplified by a factor of 7,000 (quoting GOCE mission manager R. Floberghagen). The real geoid ranges from about 100 m below the ellipsoid to 80 m above it (Rummel) — [Scientific American blog](https://www.scientificamerican.com/blog/observations/the-geoid-why-a-map-of-earths-gravity-yields-a-potato-shaped-planet/). Both of these are secondary reports; no primary GFZ/ESA page on the render was found.
- earth.nullschool.net, the benchmark for animated geophysical fields:
  - One field at a time, with a coloured overlay plus animated particles.
  - Keyboard time-stepping ("p" toggles animation), pressure levels, an HD toggle.
  - Colour scales credited to ColorBrewer, Kindlmann linear-luminance, MYCARTA and cubehelix. Inspired by the HINT.FM wind map.

  Source: [earth.nullschool.net/about](https://earth.nullschool.net/about.html)

### Inferences
- The things these pieces do well, and that our viewer should copy:
  - One quantity per view.
  - A two-hue diverging scale with a neutral centre, labelled in cm.
  - A persistent date stamp.
  - A fixed frame rate with one month per frame, or a few frames per month.
  - A narrative frame such as drought, groundwater loss or melting ice.
  - Explicit statement of any exaggeration (the ×7000).
- Arithmetic check: 100 m × 7000 = 700 km, about 11% of Earth's 6371 km radius. So the potato's bumps are about one tenth of the radius, which is why it reads as a potato.
- What the reference set lacks, and so where our viewer can differentiate:
  - Visible uncertainty.
  - An explicit marking of missing months, especially the 11-month 2017–18 gap.
  - Linked time series.
  - The spectral view (degrees and coefficients).
  - Interactive replay of the measurement principle (pair ranging, ground-track coverage).

### Gaps
- NASA Earth Observatory GRACE pieces were not reached; the search budget ran out.
- No SVS global GRACE-FO-era (2018–2026) TWS animation was confirmed. Searches found only GRACE-era global ones and the 2026 CONUS drought one.
- GravIS map and animation features (colour maps, time slider, projections) are not described on its landing page. The per-product pages were not fetched.
- No primary GFZ or ESA documentation of the potato rendering (lighting, mesh resolution) was found.

## Flow/particle techniques, and how to animate "mass transport" derived from gravity

### Takeaway
Particle advection comes in two forms. nullschool uses CPU and Canvas2D on a 1° grid. mapbox/webgl-wind uses a GPU fragment-shader state texture, with up to 1 M particles at 60 fps. Both are small enough to hand-write in a static page. Gravity gives no velocity, only mass change. A defensible "transport" field is the curl-free flux F = −∇χ with ∇²χ = ∂σ/∂t, solved exactly in spherical harmonics, and labelled as apparent, minimum-norm transport. It is physically meaningful mainly over the ocean (bottom pressure). Over land the sources are mostly vertical (precipitation minus evapotranspiration and runoff).

### Cited Findings
- **webgl-wind facts** (all from the [README](https://github.com/mapbox/webgl-wind) and source files):
  - Up to 1 million particles at 60 fps. ISC licence. No tagged releases. Wind data is pre-generated as PNG plus JSON.
  - Particle state is a square RGBA8 texture, one pixel per particle. The position is decoded as `(R/255 + B, G/255 + A)`, giving two bytes, about 16-bit fixed point, per axis.
  - The update fragment shader runs at `precision highp`. It does manual bilinear interpolation of the wind texture from 4 texels, with a comment that hardware filtering is "lower-res". Particles move with `velocity.x / cos(lat)` to correct EPSG:4326 distortion and wrap with `fract()`.
  - Particles respawn randomly with probability `dropRate + speed·dropRateBump`, defaulting to 0.003 + 0.01·speed, "to avoid degeneration".
  - Trails come from drawing onto a screen texture that is ping-ponged and faded by `fadeOpacity = 0.996` per frame. A `floor(255·c·opacity)/255` hack guarantees the trail fades to zero.
  - The default particle ramp is ColorBrewer Spectral, red–yellow–blue and rainbow-like.

  Sources: [update.frag.glsl](https://github.com/mapbox/webgl-wind/blob/master/src/shaders/update.frag.glsl); [screen.frag.glsl](https://github.com/mapbox/webgl-wind/blob/master/src/shaders/screen.frag.glsl); [index.js](https://github.com/mapbox/webgl-wind/blob/master/src/index.js)
- The explanatory Mapbox blog post ("How I built a wind map with WebGL") returned HTTP 403/404 from both blog URLs tried, so the notes rely on the source code instead — [mapbox blog URL](https://blog.mapbox.com/how-i-built-a-wind-map-with-webgl-b63022b5537f)
- **cambecc/earth** (nullschool's open-source predecessor, MIT):
  - An SVG map, with an HTML5 Canvas for particles and a second Canvas for the colour overlay.
  - GFS data at 1°, bilinearly interpolated in the browser, which the README calls costly.
  - Particle paths use finite-difference estimates of projection distortion at each pixel, so they are correct in any D3 projection.
  - A detached canvas with the globe's bounding sphere drawn into it serves as a mask.
  - GRIB2 is converted to JSON with grib2json. The README calls it a personal learning project, and it has 95 open issues.

  Source: [github.com/cambecc/earth](https://github.com/cambecc/earth)
- Line integral convolution (LIC) traces a streamline forward and backward from each pixel and averages a noise texture along it. It is computationally heavy but well suited to a GPU shader with a real-time WebGL display. Smooth noise such as simplex avoids aliasing — [US patent 11672465 (EGF visualisation)](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/11672465)
- A 2012 WebGL vector-field paper (ocean currents) used hierarchical integration with framebuffer passes. It cut cost "from linear to logarithmic" but warned that plain LIC may struggle on low-capacity devices — [EAFIT/Vicomtech 2012, SciTePress](https://www.scitepress.org/Papers/2012/38437/38437.pdf)
- Atmospheric water-balance work (the DWB project, EGU abstracts):
  - Vertically integrated moisture flux convergence equals P − E, which equals runoff plus storage change.
  - So GRACE storage change plus observed runoff can be used to solve for P − E, and the moisture-flux divergence gives an independent route to storage change.

  These are conference abstracts, not peer-reviewed papers, and the attribution of each statement to a specific abstract was not checked. Sources: [EGU2009-10219](https://meetingorganizer.copernicus.org/EGU2009/EGU2009-10219.pdf); [EGU2008-A-06976](https://meetings.copernicus.org/www.cosis.net/abstracts/EGU2008/06976/EGU2008-A-06976.pdf); [EGU2017-13632](https://meetingorganizer.copernicus.org/EGU2017/EGU2017-13632.pdf?pdf)
- Mobile GPUs often lack linear filtering of 32-bit float textures (details and percentages under the globe-rendering question). Manual bilinear sampling, as webgl-wind does, avoids depending on it — [web3dsurvey OES_texture_float_linear](https://web3dsurvey.com/webgl2/extensions/OES_texture_float_linear)

### Inferences
- **Mass-transport construction.** This is derived here, not from a source.
  - Let σ be surface mass density (kg m⁻², i.e. EWH × ρ_w) and Δσ/Δt the month-to-month tendency.
  - Assuming horizontal transport only, ∂σ/∂t + ∇·F = 0. The curl-free part F = −∇χ satisfies ∇²χ = ∂σ/∂t.
  - On the sphere this is diagonal in SH: χ_lm = −R²/(l(l+1)) · (∂σ/∂t)_lm for l ≥ 1. The l = 0 term must be zero, which holds for global mass conservation.
  - The grid of F (east and north components) comes from `geopotential`'s gradient and partials machinery applied to χ_lm.
- What the construction can and cannot show:
  - The rotational (divergence-free) part of any real flow is invisible to gravity.
  - Over land, most of ∂σ/∂t is fed vertically by P − E. A land "flux" is therefore the smoothest horizontal redistribution consistent with the tendency, not actual river or moisture flow.
  - Honest labels: "apparent (curl-free) transport implied by mass change". Over land, a better choice is to animate sources and sinks (convergence) as pulsing glyphs rather than streaks.
  - Ocean bottom pressure (a GravIS product) is the case where horizontal mass redistribution dominates and the flux is most meaningful.
- Implementation in our stack:
  - Particle state goes in RG32F or RGBA16F float textures. EXT_color_buffer_float is about 99.96% supported in WebGL2, so the RGBA8 byte-packing trick is unnecessary.
  - Flux on a 2° grid is a 180×90 RG16F texture.
  - On the globe, advect in lon/lat with the cos φ correction (as webgl-wind does), then project in the vertex shader. Wrap longitude, and respawn near the poles to avoid singular cos φ.
  - Trails use the fade-texture ping-pong. Trail opacity is a design knob for implied speed, and should be stated.
- LIC is the better static or print companion (an SVG figure via Rust). Particles are better for the live page. A cheap compromise is "animated LIC": phase-shifting the noise along streamlines.
- Speed scaling: flux magnitudes are tiny, kg m⁻¹ s⁻¹ scale. Normalise particle speed to a stated reference, e.g. "1 s of animation = 1 month", and show a speed legend. Otherwise the motion implies velocities that do not exist.

### Gaps
- No published example of a particle or flux animation derived from GRACE was found. The search budget ran out before the targeted search, so this may be novel or may simply not have been found.
- No current maintained WebGL LIC library was found; the sources are older papers.
- The Mapbox blog post explaining webgl-wind was unreachable (403/404).

## Globe rendering options in a static page: raw WebGL, three.js, deck.gl GlobeView, CesiumJS, MapLibre globe

### Takeaway
For a 2° (180×90) data field on a phone, hand-written WebGL2 is the best fit, optionally via twgl.js (~21 KB gz) or regl (~28 KB gz). A full-screen fragment shader inverse-projects each pixel (orthographic in closed form; Equal Earth by Newton iteration as d3 does). It samples a data texture and colours through a 1-D colour-map texture. three.js (~194 KB gz across two files) is worth it only for a 3D displaced-geoid mesh. deck.gl 9.4 (575 KB gz, GlobeView still experimental), MapLibre 6 (ESM-only, ~440 KB gz with worker, not on cdnjs) and CesiumJS (1.77 MB gz plus asset folders plus an ion token for imagery) are oversized. They also fit our CDN-only, tile-free constraint badly.

### Cited Findings
- **Library versions and sizes.** All rows are as of 24 Oct 2026. Versions come from the npm registry `/latest` and the cdnjs API. Gzip figures are `gzip -9` of the jsDelivr file.

  | Library | npm latest (date) | cdnjs | Main browser file (raw / gzip) | Licence |
  |---|---|---|---|---|
  | three | 0.186.1 (2026-09-24) | 0.186.1 | `build/three.module.min.js` 393 KB / 90 KB, plus `three.core.min.js` 416 KB / 104 KB | MIT |
  | deck.gl | 9.4.0 (2026-09-05) | not on cdnjs | `dist.min.js` 2.07 MB / 575 KB | MIT |
  | maplibre-gl | 6.15.0 (2026-10-18) | 6.15.0, CSS files only | `dist/maplibre-gl.mjs` 1.10 MB / 292 KB, plus worker `.mjs` 519 KB / 147 KB | BSD-3 |
  | cesium | 1.146.0 (2026-10-01) | 1.146.0 | `Build/Cesium/Cesium.js` 6.08 MB / 1.77 MB, plus Workers, Assets, Widgets | Apache-2.0 |
  | globe.gl | 2.46.2 (2026-08-22) | not on cdnjs | `globe.gl.min.js` 1.89 MB / 524 KB | MIT |
  | regl | 2.1.1 | 2.1.1 | `regl.min.js` 87 KB / 28 KB | MIT |
  | twgl.js | 7.0.0 | 7.0.1 | `twgl.min.js` 81 KB / 21.5 KB | MIT |
  | d3 | 7.9.0 | 7.9.0 | `d3.min.js` 280 KB / 92 KB | ISC |
  | @observablehq/plot | 0.6.17 | — | — | ISC |
  | weatherlayers-gl | 2026.9.1 | — | — | MPL-2.0 or commercial terms |

  Sources: [npm three](https://registry.npmjs.org/three/latest); [npm deck.gl](https://registry.npmjs.org/deck.gl/latest); [npm maplibre-gl](https://registry.npmjs.org/maplibre-gl/latest); [npm cesium](https://registry.npmjs.org/cesium/latest); [npm globe.gl](https://registry.npmjs.org/globe.gl/latest); [cdnjs API three.js](https://api.cdnjs.com/libraries/three.js?fields=version); [cdnjs API maplibre-gl files](https://api.cdnjs.com/libraries/maplibre-gl/6.15.0?fields=files); [jsDelivr three.module.min.js](https://cdn.jsdelivr.net/npm/three@0.186.1/build/three.module.min.js); [jsDelivr deck.gl dist.min.js](https://cdn.jsdelivr.net/npm/deck.gl@9.4.0/dist.min.js); [jsDelivr Cesium.js](https://cdn.jsdelivr.net/npm/cesium@1.146.0/Build/Cesium/Cesium.js)
- The jsDelivr `three.module.min.js` begins with an `import{Matrix3,Vector2,…}` statement, so it pulls in the separate core file. The real cost is roughly the sum of both files — [jsDelivr three.module.min.js](https://cdn.jsdelivr.net/npm/three@0.186.1/build/three.module.min.js)
- deck.gl GlobeView is marked experimental:
  - "No high-precision rendering at high zoom levels (> 12)".
  - Only `lnglat` coordinates.
  - HeatmapLayer and ContourLayer are unsupported. MaskExtension is unsupported.
  - TileLayer/MVTLayer support is experimental.
  - Mixing GlobeView and MapView has known issues.
  - The globe is see-through unless culled.
  - `resolution` defaults to 10° tessellation; `maxPitch` is 60.

  Source: [deck.gl GlobeView docs](https://deck.gl/docs/api-reference/core/globe-view)
- deck.gl 9.4 (5 Sep 2026) changes:
  - TerrainLayer, Tile3DLayer and TerrainExtension work on GlobeView.
  - New `project_common_position_to_flat` shader functions let custom layers sample textures in flat space while drawn on the globe.
  - GlobeController gains inertia, pitch and bearing.
  - WebGPU support across all layers is experimental and not advised for production.

  deck.gl 9.3 (13 Apr 2026) gave TimelineWidget controlled `time` and `playing` props with `onTimeChange` and `onPlayingChange` callbacks. Source: [deck.gl What's New](https://deck.gl/docs/whats-new)
- MapLibre GL JS:
  - Globe rendering arrived in 5.0 (January 2025). The globe is an "adaptive composite" view that changes to Mercator at deep zoom. Globe view is not in MapLibre Native. Sources: [MapLibre globe roadmap](https://maplibre.org/roadmap/maplibre-gl-js/globe-view/); [MapLibre newsletter Dec 2024](https://maplibre.org/news/2025-01-05-maplibre-newsletter-december-2024/)
  - v6.0.0 (22 Jul 2026) is ESM-only (`maplibre-gl.mjs`). The UMD and CSP bundles were removed, and the "worker [loads] as a real URL, so `worker-src blob:` is no longer required". Pages must use `<script type="module">`. Source: [MapLibre v6.0.0 release notes](https://github.com/maplibre/maplibre-gl-js/releases/tag/v6.0.0)
- CesiumJS quick start:
  - Loads `Cesium.js` and `widgets.css` from `cesium.com/downloads/cesiumjs/releases/1.146/…`. The NPM route needs the Workers, ThirdParty, Assets and Widgets directories served statically.
  - Global imagery, terrain and buildings come through a Cesium ion access token.
  - The page does not cover CSP or bundle size.

  Source: [Cesium quickstart](https://cesium.com/learn/cesiumjs-learn/cesiumjs-quickstart/)
- d3-geo's Equal Earth inverse uses Newton–Raphson with up to 12 iterations (`iterations = 12`, break when |δ| < ε²), with the A1..A4 polynomial constants. It ports directly to GLSL — [d3-geo equalEarth.js](https://github.com/d3/d3-geo/blob/main/src/projection/equalEarth.js). Our Rust `apps/gravviz/src/map.rs` already has `forward` and `inverse` with a round-trip test ([map.rs](../../apps/gravviz/src/map.rs)).
- Browser API support:
  - WebGL 2.0: 97.13% of global users (caniuse). Safari ≥ 15, iOS ≥ 15, Chrome ≥ 56, Firefox ≥ 51 — [caniuse webgl2](https://caniuse.com/webgl2)
  - WebGPU: 87.82% of global users. Chrome and Edge ≥ 113, Chrome Android and iOS Safari 26+ supported. Safari desktop 26 "partial". Firefox desktop listed as "disabled by default" — [caniuse webgpu](https://caniuse.com/webgpu)
- Float-texture precision:
  - OES_texture_float_linear (linear filtering of 32-bit float textures) is needed even in WebGL2. Surveyed support: about 74–75% on Android, 46–48% on iOS, about 99% on Windows. Sources: [web3dsurvey](https://web3dsurvey.com/webgl2/extensions/OES_texture_float_linear); [MDN OES_texture_float_linear](https://developer.mozilla.org/en-US/docs/Web/API/OES_texture_float_linear)
  - A 2017 Khronos list post reported that iOS claimed the extension but silently fell back to nearest — [Khronos public_webgl 2017-03](https://www.khronos.org/webgl/public-mailing-list/public_webgl/1703/msg00037.php)
  - EXT_color_buffer_float (render to float targets) is in 99.96% of WebGL2 browsers: Android 99.84%, iOS 100% — [web3dsurvey EXT_color_buffer_float](https://web3dsurvey.com/webgl2/extensions/EXT_color_buffer_float)
- Shader precision:
  - GLSL ES `mediump` is at least 16-bit with range ±2¹⁴. Some phones (e.g. iPhone 6+) really use 16 bits.
  - Some cheap phones lack `highp` in fragment shaders. Use `#ifdef GL_FRAGMENT_PRECISION_HIGH`, and test highp support via LINK_STATUS.
  - `getShaderPrecisionFormat` was unreliable on iPhone Safari as of 2020.

  Source: [WebGL Fundamentals: precision issues](https://webglfundamentals.org/webgl/lessons/webgl-precision-issues.html)

### Inferences
- **Recommended architecture (raw WebGL2, about 300 lines of JS and GLSL, no dependency):**
  1. Upload all months once as a `TEXTURE_2D_ARRAY` of R16F, 180×90×259 layers. That is about 4.2 M texels, about 8.4 MB of GPU memory, with values in mm. Half-float has an 11-bit significand, giving about 0.5 mm resolution near 1000 mm, which is well below GRACE noise. Alternatively use one 2D atlas.
  2. A full-screen triangle whose fragment shader computes (lon, lat) by inverse projection:
     - Orthographic: closed form.
     - Equal Earth: Newton iterations in `highp`. Discard pixels outside the projection outline.
  3. Bilinear-sample months k and k+1 manually. Mix by tween fraction t, unless either month is missing.
  4. Divide by the colour limit and look up a 256×1 RGBA8 texture holding Crameri vik (built from the published hex list).
  5. Coastlines go as a line strip in a second pass, or on a Canvas2D overlay as now.

  This keeps the existing light/dark CSS tokens. Only the land/ocean mask and the colour map's background colour need theming.
- Precision: at global scale float32 is ample, since the data resolution is 2° and a pixel is about 0.1° on a phone. The real risks are:
  - `mediump` in fragment shaders on low-end phones, which breaks Newton iterations and colour-limit maths. Force `highp`, or fall back to Canvas2D.
  - R32F textures with LINEAR filtering on iOS. Use R16F, which is filterable in WebGL2 (that is a GLES 3.0 property, not verified in a fetched source), or do manual bilinear filtering.
- A displaced 3D geoid mesh is three.js territory. The alternative is about 80 lines of raw WebGL with an icosphere and a vertex-shader displacement from a geoid texture. Raw WebGL saves about 190 KB gz.
- Avoid for this project:
  - Cesium: its ion token and cesium.com asset hosts conflict with the artifact's script-host allowlist, and it is the heaviest option.
  - MapLibre 6: ESM-only and not on cdnjs. Its module worker loaded "as a real URL" from a CDN may hit the same-origin restriction on `new Worker(url)`, which was not tested. It also needs a style and source plumbing we do not use.
  - deck.gl: a 575 KB bundle for one bitmap.
- deck.gl's TimelineWidget (v9.3) is a ready-made controlled play/scrub UI, but our existing slider already covers it.
- The WebGPU route (e.g. Rust `wgpu` 30.0.1 compiled to WASM) reaches only about 88% of users, and Firefox desktop is listed as disabled by default. Use WebGL2 as the baseline, with any WebGPU path optional.
- Rust-in-browser: compile `geopotential` to WASM (wasm-bindgen 0.2.129, web-sys 0.3.106) for on-device SH synthesis, e.g. degree build-up or a user-chosen functional. Rendering stays in JS/WebGL. A 2° grid to degree 60 is about 16 k points × 1.9 k coefficients, roughly 30 M multiply-adds per frame, which is feasible in WASM at interactive rates. Precomputing per-degree grids in Rust is simpler still. This is an estimate, not benchmarked.

### Gaps
- Whether MapLibre 6's module worker loads cross-origin from jsDelivr/unpkg was not tested. Its JS is not on cdnjs; cdnjs lists only CSS for 6.15.0.
- Firefox WebGPU: caniuse shows "disabled by default" for desktop through v161. The Mozilla release notes on Windows enablement could not be checked because the search budget was exhausted. Treat this as unresolved.
- deck.gl BitmapLayer behaviour on GlobeView is not mentioned in its docs.
- Cesium's CSP and worker requirements when loaded from cdnjs or jsDelivr, rather than cesium.com, are not documented on the quickstart page.
- R16F filterability in WebGL2 is stated from GLES 3.0 knowledge, not verified in a fetched source.

## Time animation design: tweening, seasonal phase/amplitude, anomaly pulses, small multiples, scrollytelling, trails, linked brushing, sonification

### Takeaway
User studies say animation is good for presenting and for engagement, but worse than static small multiples for analysis: slower, with more errors. So the viewer should pair any animation with a static, linked alternative (small multiples, phase map, time series). Tweening should interpolate data values, never colours, and must stop at real data gaps. Seasonal cycles are best shown statically with a cyclic map (romaO, four-phase) for phase and lightness for amplitude, with an optional "seasonal wave" animation synthesised from the fit.

### Cited Findings
- Robertson et al. 2008 (InfoVis, TVCG 10.1109/TVCG.2008.125) compared animation, static traces and small multiples:
  - Animation was "the fastest technique for presentation", and participants found it enjoyable and exciting, but it led to many errors.
  - For analysis, both static views were significantly faster than animation, and small multiples were the most accurate.

  Sources: [Robertson et al. 2008 (PDF)](https://faculty.cc.gatech.edu/~stasko/papers/infovis08-anim.pdf); [MSR page](https://www.microsoft.com/en-us/research/?p=145396)
- A 2019 study on mobile phones found that small multiples beat animation in 7 of 9 tasks for completion time. It notes that the 2008 study used large displays, which matters for our phone-width constraint — [mobile small-multiples study 2019, arXiv:1907.03919](https://arxiv.org/pdf/1907.03919)
- GRACE seasonal-cycle mapping in the literature (GLWS2.0) shows annual amplitude (mm) and phase as "month of maximum (1 = January)". It notes that phase mismatches of one to two months between GRACE and models are common — [GLWS2.0, arXiv:2211.17138](https://arxiv.org/pdf/2211.17138)
- Gaussian smoothing biases seasonal signals: at 800 km smoothing, annual amplitudes drop by about 25–40% and phases shift by up to 10° (basin scale) — [CU sea-level group bibliography entry](https://sealevel.colorado.edu/index.php/bibcite/reference/302). Our grids use 300 km ([tws_index.json](../../viewer/data/tws_index.json)).
- O'Neill, Rodell and Loomis (EGU 2026) argue that simplified seasonal representations recover mean annual TWS amplitude with much less attenuation and leakage — [EGU26-16019](https://meetingorganizer.copernicus.org/EGU26/EGU26-16019.html)
- Cyclic colour maps "are intended for the presentation of orientation values or angular phase data" (Kovesi). Problems with the HSV hue wheel:
  - Uneven contrast, with flat spots from cyan to yellow and red to magenta.
  - Light secondaries that create false anomalies.
  - A three-way red/green/blue split that does not match four-quadrant thinking such as compass points or seasons.

  Kovesi recommends four-colour cyclic maps, e.g. blue → darkened yellow → dark red → pink → blue, with smoothed lightness reversals — [Kovesi 2015, arXiv:1509.03700](https://arxiv.org/abs/1509.03700)
- Crameri's cyclic maps (e.g. romaO) are perceptually uniform, CVD-friendly and readable in B&W, contain no black, and are "with circular data only" — [Crameri cyclic colour maps](https://www.fabiocrameri.ch/cycliccolourmaps/)
- Scrollytelling: scrollama 3.2.0 is on cdnjs (npm last modified 2022) — [cdnjs API scrollama](https://api.cdnjs.com/libraries/scrollama?fields=version); [npm scrollama](https://registry.npmjs.org/scrollama/latest)
- Sonification:
  - Tone.js npm latest is 15.1.22 (modified 2026-10-22, about 79 KB gz `build/Tone.js`) — [npm tone](https://registry.npmjs.org/tone/latest); [jsDelivr Tone.js](https://cdn.jsdelivr.net/npm/tone@15.1.22/build/Tone.js)
  - cdnjs reports Tone version 15.5.56, which is inconsistent with npm "latest" and probably a prerelease tag — [cdnjs API tone](https://api.cdnjs.com/libraries/tone?fields=version)
  - A search for GRACE water-storage sonification found no prior example.

### Inferences
- **Tweening rules:**
  - Interpolate data (mm) per cell in the shader, then colour-map. Interpolating RGB produces off-map hues through the neutral centre.
  - Use linear interpolation in time between consecutive available months only. Across a gap (any missing month, especially the 11-month 2017-06 → 2018-06 gap) do not interpolate. Show the last frame faded or hatched with a "no data: GRACE → GRACE-FO gap" banner, and mark missing months on the slider track.
  - Use a time-proportional playhead: frames per real month, so a 3-month gap takes 3× as long, rather than frames per available field. That keeps seasonal rhythm honest.
  - Do not ease (ease-in/out) between monthly values; easing invents sub-monthly dynamics. At most, cross-fade over 15–25% of the inter-month interval to reduce flicker.
- **Seasonal phase/amplitude:**
  - Fit per cell (Rust, `gravimetry`): x(t) = a + b·t + A₁cos(2πt − φ₁) + A₂cos(4πt − φ₂) + residual.
  - Static key view: hue from romaO (or a four-phase cyclic map, so the four seasons are the quadrants) shows φ₁. Lightness or opacity shows A₁, which makes the map bivariate.
  - Animated view: render A₁cos(2πτ − φ₁) for a synthetic year τ looping at about 4–8 s. Seasonal water "waves" then visibly sweep, e.g. across the Amazon and the monsoon belt.
  - Caveat to display: phases and amplitudes are biased by the 300 km smoothing.
- **Anomaly pulse:** cells where the de-seasonalised residual exceeds k·σ (per-cell σ) get a brief ring or outline pulse of at most about 0.5 s. Use opacity, not scale (see the accessibility section). Never more than 3 flashes per second.
- **Trails and accumulation:**
  - A "cumulative change since 2002-04" mode that integrates the trend.
  - A "record" layer: cells reaching a new running minimum or maximum light up and stay lit, faded, as a trail.
  - A "time of emergence" map: the year the trend exceeded 2σ.
- **Small multiples:** a 12-column (month) × 25-row (year) grid of 90×45-pixel thumbnails from the same R16F texture array, drawn with one draw call per thumbnail via a scissor and viewport, or as a single instanced pass.
  - At phone width this is a scrolling list of year rows.
  - Tapping a cell sets the main map's month. This is the analysis-friendly complement to playback.
- **Linked brushing:**
  - Brushing a time range on the point or basin series sets the map to the mean (or trend) over the brush.
  - Brushing a region on the map (lasso or basin pick) updates an area-weighted series, using cos φ weights.
  - Hovering the series scrubs the map. Our viewer already has point time series, so the main additions are range-brush and region-average.
- **Scrollytelling:** three to five steps (global trend → Greenland/Antarctica → California → Middle East groundwater → the 2017–18 gap and Swarm). It can be built with scrollama, or with a plain IntersectionObserver and no dependency.
- **Sonification:** a basin series mapped to pitch, with loudness proportional to |anomaly| and a click at each January. A pure WebAudio OscillatorNode avoids Tone.js. Treat it as an accessibility and engagement extra, not an analysis tool.

### Gaps
- Harrower & Fabrikant (2008) and Tversky et al. (2002) on map-animation cognition and change blindness were not fetched; the search budget was exhausted.
- No source was found for the specific frame-rate or dwell-time recommendations for monthly geoscience animations. The values above are design judgement.
- No precedent for GRACE sonification was found.
- No Observable notebooks were found, because of the search limit.

## Spectral and abstract animations: SH build-up, rotating exaggerated geoid, coefficient triangles, error-spectrum races, ground tracks, aliasing

### Takeaway
Our crates already compute every quantity these explainers need: per-degree synthesis, degree variances, error triangles, orbit positions, degree-to-frequency mapping. The animations are therefore mostly a front-end task. The strongest storytelling ones are:
- an SH build-up (degree 2 → 60 or → 180) synchronised with a degree-spectrum cursor;
- an exaggeration slider on a rotating geoid, with the factor always labelled;
- a ground-track coverage accumulation driven by `satgrav::sim::Pair::positions`.

The literature's central message, that temporal aliasing rather than sensor noise limits GRACE-type missions, is the best "race chart" story.

### Cited Findings
- The `geopotential` crate offers normalised Legendre functions to degree 2700, potential and gradient synthesis, grids of geoid, gravity and water height, `degree_variances` and `degree_amplitudes`, and a `truncated` method — [crates/geopotential/Cargo.toml](../../crates/geopotential/Cargo.toml)
- `satgrav` provides:
  - `orbit_period` (circular orbit).
  - `degree_at_frequency` (n ≈ f·T, one cycle per revolution being degree 1).
  - `resolution(n) = πR/n`.
  - `sim::Pair::positions`: Earth-fixed positions of the leading and trailing satellites t seconds after A's northbound equator crossing at 0° longitude.
  - `simulate`: closed-loop least squares, with memory of about 1.2 GB at degree 40 for one month of 30 s samples.

  Sources: [crates/satgrav/src/lib.rs](../../crates/satgrav/src/lib.rs); [crates/satgrav/src/sim.rs](../../crates/satgrav/src/sim.rs); [crates/satgrav/src/errors.rs](../../crates/satgrav/src/errors.rs)
- Existing figures and data: `figures/noise_spectra.csv/.svg`, `figures/error_triangle_2026-07.svg`, `figures/gradiometry_errors.csv/.svg` and `figures/sst_closed_loop.csv/.svg` — [figures/](../../figures)
- Temporal aliasing dominates: in Zingerle et al. 2024 closed-loop simulations, sensor noise only shows once the time-variable signal is scaled below 10% (GRACE-FO-class sensors) or below 0.1% (a 10⁻¹² m s⁻² Hz⁻¹ᐟ² atom accelerometer) — [Zingerle et al. 2024](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024EA003630), as summarised in [our report](../../reports/Satellite%20gravimetry%20data%20and%20Rust%20tooling.md)
- The potato's exaggeration is ×7000, and the geoid spans about −100 m to +80 m — [Scientific American](https://www.scientificamerican.com/blog/observations/the-geoid-why-a-map-of-earths-gravity-yields-a-potato-shaped-planet/)
- The GRACE-FO pair collects a global map about every 30 days. The orbit is 490 → 475 km at 89° inclination. Sources: [JPL resource 40](https://grace.jpl.nasa.gov/resources/40/grace-fo-from-range-observations-to-global-mass-change); [DLR](https://www.dlr.de/en/rb/research-operation/missions/earth-observation-science/grace-follow-on)
- Swarm is the only gravity source spanning the 2017–2018 gap: 145 monthly fields, agreeing with GRACE to about 4 cm EWH over land at 750 km smoothing — [our report, citing Teixeira da Encarnação et al. 2020](https://essd.copernicus.org/articles/12/1385/2020/)

### Inferences
- **SH build-up:**
  - Precompute per-degree contribution grids g_l(lon, lat) for l = 2..60. That is 59 × 16,200 values, about 1 M, or about 2 MB as f16.
  - Animate the partial sum S_L = Σ_{l≤L} g_l with L ticking up about 4 per second.
  - Side panel: the degree-amplitude curve with a moving cursor at L and the half-wavelength πR/L in km.
  - For the static geoid, go to L ≈ 180–360 by computing in WASM or streaming precomputed grids.
  - The ringing at low L (Gibbs) is itself instructive, and it motivates the Gaussian filter. A second pass can animate the Gaussian weight curve while the map smooths.
- **Rotating exaggerated geoid:**
  - An icosphere of about 40 k vertices, displaced by N(lon, lat) × k in the vertex shader. Colour by N with vik or roma. Lambert shading reveals the bumps.
  - An exaggeration slider runs from k = 1 (a visually perfect sphere) to 7000 and beyond, always showing "×k" on screen.
  - Under reduced-motion, rotation is user-driven only.
  - A monthly-anomaly version (EWH displaced, k about 10⁵–10⁶) is possible but strongly needs the exaggeration label.
- **Coefficient-triangle animation:**
  - The degree (0–60) × order (−60..60, S on the left and C on the right) triangle of monthly ΔC/ΔS, or of their signal-to-error ratio, on a log-scaled diverging map.
  - Played through months, the north–south striping appears as high-order noise at fixed m patterns. Destriping can be shown as a before/after toggle.
  - It costs only 3,721 cells per month, which is trivial in Canvas2D.
- **Error-spectrum race chart:**
  - Degree-RMS curves (EWH mm on a log axis) for: signal; GRACE KBR; GRACE-FO LRI; CAI scenarios; aliasing with and without de-aliasing; multi-pair.
  - Animate as successive reveals or as a "race" in which curves morph when a parameter changes (altitude, separation, accelerometer floor).
  - The crossover degree (signal = error) moves, giving "achievable resolution πR/n_cross". This visualises the report's main message.
- **Ground-track coverage:**
  - Evaluate `Pair::positions` at 30–60 s steps for 30 days, about 43–86 k points.
  - Animate the ground track with a short fading trail, the inter-satellite line, and a coverage density raster (hits per 2° cell) that fills in.
  - Show the gaps of short sub-cycles (e.g. 3-day or 7-day) versus the full month, and the polar convergence from the 89° inclination.
  - A time-warp control is needed: one orbit takes about 94 min.
- **Aliasing explainer:**
  - Panel 1: a sine of period P (e.g. a tidal or sub-monthly signal) sampled at the orbit's revisit interval Δ. The samples trace a slow alias of period 1/|1/P − k/Δ|.
  - Panel 2: map frames from a synthetic field showing how a fast signal leaks into monthly solutions as stripes. The `simulate` output with and without the time-variable truth would supply this.
  - Pair it with the existing `sst_closed_loop` figure.
- **Mission-transition animation:** GRACE → gap (Swarm, fading to lower resolution and higher noise) → GRACE-FO. The resolution change is shown honestly, not hidden.

### Gaps
- The published alias periods of specific tides (e.g. S2, K1) for GRACE and GRACE-FO orbits were not retrieved; the search budget was exhausted. Compute them from the orbit's repeat or sub-cycle, or cite a primary paper, before publishing numbers.
- No existing web example of an animated SH build-up or coefficient-triangle animation for gravity was found. This is not evidence there are none.

## Accessibility and honesty: perceptually uniform diverging maps, avoiding rainbow, reduced motion, uncertainty, misleading interpolation

### Takeaway
Use Crameri's perceptually uniform, CVD-safe maps. For EWH anomaly use vik, or roma for a less saturated look; for phase use romaO. Smooth the diverging map's centre (Kovesi). Honour `prefers-reduced-motion` by starting paused and replacing rotation and scaling with fades. Give every auto-playing animation a pause control (WCAG 2.2.2, Level A, applies after 5 s). Show uncertainty explicitly: noise-floor hatching, or HOP-style animated realisations, which outperform error bars for comparisons.

### Cited Findings
- Crameri, Shephard & Heron 2020 (Nature Communications 11:5444):
  - Colour-vision deficiency affects 0.5% of women and 8% of men worldwide.
  - Rainbow is "the most-widely used colour map".
  - Visual error from non-uniform maps like jet "can diverge from an objective representation by more than seven percent". Turbo fails uniformity despite being perceptually ordered.
  - For data diverging about a central value, "a divergent or a multi-sequential colour map should be chosen".
  - Perceptual uniformity is checked with flat CIEDE2000 ΔE along the map.

  Source: [Crameri et al. 2020 (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7595127)
- Scientific colour maps:
  - Version 8.0.1, 5 Oct 2023, Zenodo DOI 10.5281/zenodo.1243862.
  - Perceptually uniform and ordered, readable by colour-vision-deficient people, and readable in B&W print.
  - Extra export formats include Adobe .ase and others.

  Source: [fabiocrameri.ch/colourmaps](https://www.fabiocrameri.ch/colourmaps/). The fetched page did not list the diverging map names: vik, roma, broc, cork, berlin, etc.
- Kovesi on diverging maps:
  - The common blue–white–red map has a lightness-gradient reversal at its centre that "induces the perception of a false feature". Smoothing the lightness (e.g. a Gaussian filter) removes it, at a small cost.
  - Values straddling the centre become nearly isoluminant.
  - His sine-on-ramp test image (features at 10% of range) reveals flat spots.

  Source: [Kovesi 2015](https://arxiv.org/abs/1509.03700)
- NASA SVS moved from brown/green (2002–2015 TWSA) to orange/blue (2015–16 TWSA) for dry/wet — [SVS 4415](https://svs.gsfc.nasa.gov/4415)
- `prefers-reduced-motion`:
  - Baseline "widely available" since January 2020.
  - `reduce` means the user wants motion "removed, reduced, or replaced".
  - "Scaling or panning large objects can be vestibular motion triggers". MDN's example replaces a scaling pulse with an opacity "dissolve".
  - The OS switches are on macOS, iOS, Windows 10/11, Android 9+, GNOME and KDE.

  Source: [MDN prefers-reduced-motion](https://developer.mozilla.org/en-US/docs/Web/CSS/@media/prefers-reduced-motion)
- WCAG 2.2 SC 2.2.2 Pause, Stop, Hide (Level A): moving content that starts automatically, lasts more than 5 s and is shown alongside other content needs a mechanism to pause, stop or hide it. A single control for all motion is recommended. Related criteria are 2.3.1 "Three Flashes or Below Threshold" and 2.3.3 "Animation from Interactions" — [W3C Understanding 2.2.2](https://www.w3.org/WAI/WCAG22/Understanding/pause-stop-hide.html)
- Hypothetical Outcome Plots (HOPs) "animate a finite set of individual draws" from a distribution. For two or three quantities, people judged far more accurately with HOPs than with error bars or violins; single-quantity accuracy was similar — [Hullman lab HOPs](https://www.hullmanlab.northwestern.edu/paper%20discrete-outcome/2015/01/01/hops.html)
- With HOPs, untrained viewers judged trends better than with static aggregates when evidence was weak — [Kale et al. 2018 InfoVis](https://mucollective.northwestern.edu/files/2018-HOPsTrends-InfoVis.pdf)
- NetHOPs: 51 experts' estimates were within 11% of ground truth on average — [Zhang, Adar & Hullman 2021](https://ar5iv.labs.arxiv.org/html/2108.09870)
- webgl-wind's default ramp is ColorBrewer Spectral (rainbow-like), so copying its defaults would import a non-recommended map — [webgl-wind index.js](https://github.com/mapbox/webgl-wind/blob/master/src/index.js)

### Inferences
- **Colour choices:**
  - EWH anomaly: vik (light centre) in light mode. In dark mode, a dark-centre diverging map (Crameri berlin, vanimo or lisbon) avoids a glaring white "zero" sea. If one map must serve both themes, roma is a fair compromise. Map names are from domain knowledge; not confirmed on a fetched page.
  - Trend: the same diverging map with a symmetric limit.
  - Amplitude: sequential batlow, or a lightness ramp.
  - Phase: romaO, or a Kovesi four-phase map so the seasons are quadrants. Bivariate phase × amplitude uses the cyclic hue at full amplitude, fading to the background at zero.
  - Uncertainty: grey hatch.
- **Bake the colour limit into the legend and keep it fixed during playback.** Auto-rescaling per frame makes quiet months look as dramatic as extreme ones. This is a common animation lie. Show a "values beyond ±L are clipped" tick, which the viewer already does with ≤/≥.
- **Reduced motion:**
  - If `reduce` is set: no autoplay; no globe auto-rotation or inertia; pulses become static outlines; particles become a static LIC or streamline image.
  - Month changes become instant swaps, or crossfades of at most 200 ms. Always honour the live `change` event, as the viewer already does for colour scheme.
- **Uncertainty in animation:**
  1. A noise-floor mask: hatch or desaturate cells where |anomaly| is below the monthly noise σ. The per-month σ can come from calibrated errors, or from the residual RMS over the ocean far from coasts.
  2. A HOP mode: cycle at about 2–4 frames per second through realisations (field plus a noise sample from `satgrav` simulated mission errors or formal covariance). Persistent features stay, and noise twinkles.
  3. A ribbon on the time series at ±1σ.

  Also state the leakage, GIA and smoothing caveats in the footer; the viewer already notes "GIA not removed".
- **Misleading interpolation:** interpolation in time and space both need disclosure.
  - Time: no tweening across missing months; mark gaps on the slider.
  - Space: label it "2° grid, 300 km smoothing", and do not render sub-grid detail with smooth shading that implies higher resolution. Bilinear interpolation is acceptable but should be described. Nearest (blocky) is the most honest at full zoom.
- **Flash safety:** keep any blink or pulse below 3 per second and confined to small areas, per WCAG 2.3.1.

### Gaps
- The full Crameri map list and the berlin/vanimo dark-centre suitability were not verified on a fetched page.
- The exact numeric flash thresholds in WCAG 2.3.1 (area, luminance) were not fetched.
- No study of uncertainty animation specifically for geophysical raster time series was found, owing to the search limit.

## Export: shareable GIF/MP4/APNG from canvas in the browser, or from Rust

### Takeaway
For frame-exact, smooth files, render offline. Step the month and tween deterministically, then encode each frame with WebCodecs, through Mediabunny's `CanvasSource` or the 2026 CCapture.js v2 rewrite, which does WebCodecs MP4 with a virtual clock. Use `gifenc` (about 4 KB gz) for GIF and UPNG.js for APNG. `canvas.captureStream()` with MediaRecorder is universally available but real-time only, and its format varies by browser. mp4-muxer is deprecated. In Rust, render frames with `tiny-skia` or `resvg` and write them with the `gif` or `png` (APNG) crates, or pipe them to an external encoder.

### Cited Findings
- `HTMLCanvasElement.captureStream(frameRate)` returns a MediaStream. A frameRate of 0 means frames are captured only on `requestFrame()`. It is Baseline "widely available" since January 2020. It throws SecurityError on a non-origin-clean canvas — [MDN captureStream](https://developer.mozilla.org/en-US/docs/Web/API/HTMLCanvasElement/captureStream)
- WebCodecs VideoEncoder:
  - Chromium ≥ 94, Firefox desktop ≥ 130, Safari ≥ 16.4 per the caniuse `VideoEncoder.close` table.
  - Firefox for Android is not supported.
  - Encoding codec support is narrower than decoding; use `VideoEncoder.isConfigSupported()`.

  Source: [caniuse VideoEncoder](https://caniuse.com/mdn-api_videoencoder_close). A blog claims full audio+video parity only in Safari 26, which conflicts with caniuse — [Digital Samba blog](https://www.digitalsamba.com/blog/webcodecs-api-explained)
- mp4-muxer 5.2.2: "This library is deprecated … in favor of Mediabunny" — [mp4-muxer README](https://cdn.jsdelivr.net/npm/mp4-muxer@5.2.2/README.md). The same applies to webm-muxer 5.1.4 (same author, same date) — [npm webm-muxer](https://registry.npmjs.org/webm-muxer/latest)
- Mediabunny:
  - Version 1.62.1 (npm, 2026-10-22, MPL-2.0). The full ESM bundle `mediabunny.min.mjs` is 704 KB raw / 181 KB gz, and the library is tree-shakeable.
  - Not on cdnjs; available on jsDelivr.
  - `CanvasSource(canvas, {codec, quality})` with `await add(timestamp, duration)` takes explicit timestamps and gives backpressure, which suits non-real-time rendering.

  Sources: [Mediabunny media sources guide](https://mediabunny.dev/guide/media-sources); [npm mediabunny](https://registry.npmjs.org/mediabunny/latest); [jsDelivr bundle](https://cdn.jsdelivr.net/npm/mediabunny@1.62.1/dist/bundles/mediabunny.min.mjs)
- CCapture.js 2.0.0 (npm, modified 2026-07-27, by Jaume Sanchez/spite):
  - "a ground-up rewrite … (ES modules, WebCodecs, GPU motion blur)".
  - A virtual clock advances one fixed step per captured frame, so slow-to-render frames still give a smooth file at an exact frame rate. `new CCapture({format:"mp4", framerate:60})`.
  - The UMD build `build/ccapture.umd.min.js` is 117 KB raw / 34 KB gz. Not on cdnjs.

  Sources: [CCapture v2 README](https://cdn.jsdelivr.net/npm/ccapture.js@2.0.0/README.md); [github.com/spite/ccapture.js](https://github.com/spite/ccapture.js)
- Other encoders:
  - gifenc 1.0.3: 9 KB raw / 3.9 KB gz ESM; jsDelivr only.
  - gif.js 0.2.0: on cdnjs.
  - UPNG.js 2.1.0, for APNG: on cdnjs.
  - h264-mp4-encoder 1.0.12: WASM, last modified 2022.

  Sources: [npm gifenc](https://registry.npmjs.org/gifenc/latest); [cdnjs gif.js](https://api.cdnjs.com/libraries/gif.js?fields=version); [cdnjs upng-js](https://api.cdnjs.com/libraries/upng-js?fields=version); [npm h264-mp4-encoder](https://registry.npmjs.org/h264-mp4-encoder/latest)
- Rust crates (crates.io, 24 Oct 2026):

  | Crate | Version | Updated | Note |
  |---|---|---|---|
  | `gif` | 0.14.2 | 2026-04-09 | |
  | `png` | 0.18.1 | 2026-02-14 | already a dependency of gravviz |
  | `image` | 0.25.10 | | |
  | `resvg` | 0.48.1 | 2026-08-02 | SVG → raster |
  | `tiny-skia` | 0.12.0 | | |
  | `openh264` | 0.9.9 | 2026-10-15 | H.264 bindings |
  | `rav1e` | 0.8.1 | | AV1 encoder |
  | `minimp4` | 0.1.2 | 2024 | |
  | `wgpu` | 30.0.1 | 2026-08-22 | |
  | `plotters` | 0.3.7 | 2024-09 | |

  Sources: [crates.io gif](https://crates.io/crates/gif); [crates.io png](https://crates.io/crates/png); [crates.io resvg](https://crates.io/crates/resvg); [crates.io openh264](https://crates.io/crates/openh264); [crates.io rav1e](https://crates.io/crates/rav1e); [crates.io wgpu](https://crates.io/crates/wgpu); [apps/gravviz/src/svg.rs `png_data_uri`](../../apps/gravviz/src/svg.rs)

### Inferences
- **In-page "Save clip" button (recommended path):**
  - Pause the live loop and render N frames off-screen at a fixed size (e.g. 1080×1080 for social, 1280×720 for slides).
  - Feed `CanvasSource.add(i/fps, 1/fps)`. Import Mediabunny from jsDelivr as ESM, selecting only the MP4 output and the canvas source to keep it small.
  - Offer a Blob download. Feature-detect `VideoEncoder`. If it is missing (Firefox Android), fall back to gifenc GIF at 480 px and 10–15 fps, or to MediaRecorder WebM in real time.
  - Burn in the date stamp, colour bar, units and source credit ("COST-G RL02.1, CC BY 4.0"), so the clip stays honest when shared out of context.
- GIF's 256-colour palette is adequate for a single continuous diverging map. Pre-build the palette from the colour map's 256 entries plus the UI colours, rather than per-frame quantisation; this avoids flicker. APNG via UPNG.js gives lossless quality for short loops, with larger files.
- **Rust path for reports:** the gravviz app already makes SVG and PNG data URIs. Extend it to render each month as RGBA (tiny-skia, or resvg rasterising the existing SVG) and write:
  - an animated GIF (`gif` crate, global palette from vik);
  - an APNG (`png` crate; its APNG encoding support is believed present but was not verified);
  - PNG frames for `ffmpeg -framerate 12 -i %04d.png -pix_fmt yuv420p out.mp4`, run outside the build.

  Avoid bundling openh264 or rav1e unless needed: heavy compiles, and the user's machine rules require capped builds.
- File-size budget: 259 months × 1080² at 12 fps is about 22 s of H.264. A few MB is plausible, but this was not measured.

### Gaps
- MediaRecorder's MIME types per browser (WebM in Chromium/Firefox, MP4 in Safari) are from general knowledge, as the MDN captureStream page did not cover them. Check with `MediaRecorder.isTypeSupported()`.
- APNG encoding support in the Rust `png` 0.18 crate was not verified from documentation.
- The Mediabunny tree-shaken size for MP4 + CanvasSource only was not measured. 181 KB gz is the full bundle.

## Concrete animation ideas tailored to our gravity data (what each needs)

### Takeaway
There are fourteen ideas. Each is buildable as a static page with raw WebGL2 or Canvas2D plus our Rust precomputation, and each states its data needs. The first five extend the existing viewer with little new data. Ideas 6–11 use spectral and mission outputs the crates already compute. Ideas 12–14 are explainers and outreach.

### Cited Findings
- The current data are 259 monthly 2° EWH grids in mm. Processing: degree 60, destriped, 300 km Gaussian, relative to 2004–2009. 35 months are missing, including 2017-06 → 2018-06 — [viewer/data/tws_index.json](../../viewer/data/tws_index.json)
- The crates provide degree spectra, per-degree truncation, grids of any functional, Gaussian smoothing and destriping, error propagation, orbit positions and closed-loop simulation — [crates/*/Cargo.toml](../../crates); [crates/satgrav/src/sim.rs](../../crates/satgrav/src/sim.rs)
- The effectiveness evidence and design rules behind these ideas are cited in the sections above: Robertson 2008, HOPs, Crameri 2020, Kovesi 2015, MDN, WCAG, webgl-wind.

### Inferences
Each idea lists its data needs, then its build and notes.

1. **Honest month tween.**
   - Data: the existing grids, converted to binary Int16 or Float16 (`.bin`, about 8.4 MB, or about 3–4 MB gzipped, instead of 5.7 MB base64 JSON), plus the missing-month list.
   - Build: R16F texture array; shader mix between adjacent months; gap frames hatched; time-proportional playhead.
2. **Seasonal "breathing" map.**
   - Data: per-cell annual and semi-annual amplitude and phase plus trend (Rust least squares over 259 months), as 4–6 float grids, about 400 KB.
   - Build: static romaO phase × amplitude key map; animated synthetic year A cos(2πτ − φ). The seasonal wave visibly travels.
3. **Trend reveal and "cumulative change since 2002".**
   - Data: per-cell linear trend and its σ; cumulative anomaly per month (already in the grids).
   - Build: the map accumulates; cells whose trend exceeds 2σ get an outline; a "time of emergence" year map.
4. **Anomaly pulse / extremes tracker.**
   - Data: de-seasonalised residuals and a per-cell σ.
   - Build: during playback, cells beyond ±2σ get a 0.5 s opacity ring and a running ticker of the largest anomalies. Under reduced motion it is a static outline. At most 3 pulses per second.
5. **Small multiples calendar plus linked brushing.**
   - Data: existing grids, and basin or region polygons for region averages.
   - Build: a year × month thumbnail grid; tap to set the month; brush on the series to show the range mean; lasso on the map for a weighted series.
6. **Apparent mass-transport particles (ocean first).**
   - Data: monthly SH coefficients of Δσ/Δt; Poisson solve χ_lm = −R²(Δσ/Δt)_lm / (l(l+1)); a grid of ∇χ (east and north) at 2°; ideally ocean bottom pressure (GravIS OBP, or our own with an ocean mask).
   - Build: webgl-wind-style particle state textures (float), a fade-trail ping-pong and globe projection. Label it "apparent curl-free transport".
7. **Spherical-harmonic build-up.**
   - Data: per-degree contribution grids for l = 2..60 (about 2 MB f16), or WASM synthesis; degree amplitudes for the side chart.
   - Build: the partial sum animates while a cursor moves on the degree spectrum, with a πR/L readout. A second act animates the Gaussian filter's weights.
8. **Rotating exaggerated geoid ("potato") with an exaggeration slider.**
   - Data: a static geoid grid (e.g. EIGEN-6C4 or EGM2008 truncated to about 180, at 1°; 360×180 floats), computed by `geopotential`.
   - Build: three.js or raw WebGL icosphere displacement, with "×k" always shown. Auto-rotation is off under reduced motion.
9. **Coefficient-triangle playback.**
   - Data: monthly ΔC/ΔS to degree and order 60 (3,721 per month), raw and destriped.
   - Build: the triangle animates with month, linked to the map. Stripes appear as order-dependent structure, with a destriping toggle.
10. **Error-spectrum race.**
    - Data: degree-RMS signal and error curves per scenario (`figures/noise_spectra.csv`, `sst_closed_loop.csv`, `gradiometry_errors.csv`).
    - Build: curves morph as altitude, separation, accelerometer floor or aliasing toggle change; the crossover degree is marked, with resolution in km.
11. **Uncertainty HOPs for a month.**
    - Data: K about 20 noise realisations per month, from `satgrav` simulated mission errors or formal covariance, synthesised to the grid. Or one error-degree-variance curve to draw random coefficients from.
    - Build: frames cycle at 2–4 per second, so robust features persist and noise twinkles. A static alternative shows the probability that |signal| > 0.
12. **GRACE pair ground-track coverage.**
    - Data: `Pair::positions` sampled every 30–60 s for 30 days, about 1 MB as f32 lon/lat pairs, plus a hit-count grid.
    - Build: a dot with a fading trail; the coverage raster fills in; scrubbing compares 1-day, 7-day and 30-day sampling; the 89° polar gap is visible.
13. **Aliasing explainer.**
    - Data: an analytic sine and sampling interval, plus `simulate` runs with and without the time-variable truth.
    - Build: panel 1 shows dots on a fast sine revealing a slow alias; panel 2 shows how stripes appear in the monthly map.
14. **Scrollytelling story with sonified basins.**
    - Data: the existing grids, 4–5 basin masks and series.
    - Build: IntersectionObserver or scrollama steps with fly-to on the globe; an optional WebAudio pitch track of the basin series (click on each January).
    - Export uses the same frame path as the "Save clip" button.

Shared plumbing for all 14:
- One WebGL2 field renderer: texture array, inverse projection and colour-map LUT.
- One time controller: real-month clock, gap awareness, reduced-motion aware, with the WCAG pause control.
- One linked-selection store (month, range, region).
- One export path (Mediabunny or CCapture, with a GIF fallback).

### Gaps
- No performance benchmarks were run, per the "web research only" constraint. Frame-rate claims for phones are extrapolated from webgl-wind's "1 M particles at 60 fps" claim and our small grid size.
- Ocean bottom pressure inputs (GravIS OBP grids, or an ocean mask with GAD restoration) were not checked for availability in our pipeline.
