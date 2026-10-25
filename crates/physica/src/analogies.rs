//! The same equations across gravity, electromagnetism and plasma, each stated as a test that
//! calls the crates on both sides.
//!
//! **One spherical-harmonic engine.** Earth's gravity and its main magnetic field are both
//! potential fields expanded in spherical harmonics, so one engine evaluates both through
//! [`geopotential::Field`]: gravity is `∇V`, the magnetic field `−∇V` with `GM = a²` and
//! `C̄ₙₘ = gₙᵐ/√(2n+1)`.
//!
//! ```
//! use geopotential::{Coefficients, Field};
//! use magnetics::main_field::wmm2025;
//!
//! let mut gravity = Coefficients::zeros(3.986_004_415e14, 6_378_136.3, 2);
//! gravity.set(0, 0, 1.0, 0.0);
//! gravity.set(2, 0, -4.841_694_95e-4, 0.0);
//! let magnetism = wmm2025().at(2025.0);
//! let fields: [&dyn Field; 2] = [&gravity, &magnetism];
//! let r = 6_371_200.0 + 450e3; // a Swarm-like orbit
//! let [g, b] = fields.map(|f| f.magnitude(r, 0.5, 1.0));
//! assert!((8.0..8.8).contains(&g)); // m/s²
//! assert!((20_000.0..60_000.0).contains(&b)); // nT
//! ```
//!
//! **Spinning dipoles.** A spinning charge and a spinning mass both make a dipole field, twice as
//! strong at the pole as at the equator: the geomagnetic axial dipole and the gravitomagnetic
//! field of Earth's spin.
//!
//! ```
//! use geopotential::{Coefficients, Field};
//! use gravitics::gem::GravitomagneticDipole;
//! use magnetics::main_field::{MainField, REFERENCE_RADIUS};
//!
//! let a = REFERENCE_RADIUS;
//! let mut axial = Coefficients::zeros(a * a, a, 1);
//! axial.set(1, 0, -29_351.8 / 3f64.sqrt(), 0.0); // WMM2025's g₁⁰, nT
//! let magnetic = MainField(axial);
//! let gravitomagnetic = GravitomagneticDipole { j: 5.86e33 };
//! let pole = core::f64::consts::FRAC_PI_2 - 1e-9;
//! for field in [&magnetic as &dyn Field, &gravitomagnetic] {
//!     assert!((field.magnitude(a, pole, 0.0) / field.magnitude(a, 0.0, 0.0) - 2.0).abs() < 1e-9);
//! }
//! ```
//!
//! **The spin-2 factor inside a spinning shell.** Inside a spinning charged shell the magnetic
//! field is `(2/3)(μ₀/4π)Qω/R`; inside a spinning massive shell, frames are dragged at
//! `4GMω/(3c²R)`. Replace `μ₀/4π` by `G/c²` and `Q` by `M`, and gravity comes out exactly twice as
//! strong: the spin-2 graviton couples to mass-currents with twice the strength.
//!
//! ```
//! use metrology::{C, G, MU_0};
//! use plasmoid::fusion::spinning_shell_field;
//!
//! let (m, r, omega) = (1e4, 1.0, 300.0);
//! let analogue = spinning_shell_field(m, r, omega) * (G / (C * C)) / (MU_0 / (4.0 * core::f64::consts::PI));
//! let frame_dragging = 4.0 * G * m * omega / (3.0 * C * C * r);
//! assert!((frame_dragging / analogue - 2.0).abs() < 1e-12);
//! assert!((frame_dragging / 2.970e-21 - 1.0).abs() < 1e-3); // rad/s
//! ```
//!
//! **One charge-to-mass constant.** `√(4πε₀G)` is where electric repulsion matches gravity
//! (the extremal black hole and the weak-gravity scale). The drift speed at which two protons'
//! magnetic attraction matches their gravity, times the proton's ratio to that constant, is the
//! speed of light exactly.
//!
//! ```
//! use gravitics::gem::wgc_ratio;
//! use metrology::{C, E_CHARGE, M_P};
//! use plasmoid::pinch::equal_attraction_drift;
//!
//! assert!((equal_attraction_drift() * wgc_ratio(E_CHARGE, M_P) / C - 1.0).abs() < 1e-9);
//! ```
//!
//! **Pinch and filament.** The Bennett relation that holds a Z-pinch together, `μ₀I²/4π`, is the
//! Ostriker relation that holds a self-gravitating gas filament together, `Gλ²`, with magnetism in
//! place of gravity. In a 1 MA pinch, magnetism wins by 10²⁷.
//!
//! ```
//! use plasmoid::pinch::{bennett_line_density, magnetic_over_gravitational};
//!
//! let line_mass = bennett_line_density(1e6, 2000.0) * 3.3436e-27; // deuterium at 2 keV
//! assert!(magnetic_over_gravitational(1e6, line_mass) > 1e27);
//! ```
//!
//! **Schwinger and Unruh.** The exponential suppression of pair creation in a field `E`,
//! `e^(−πE_S/E)`, is exactly a Boltzmann factor `e^(−m_ec²/(2k_BT_U))` at the Unruh temperature
//! of an electron accelerated by that field.
//!
//! ```
//! use gravitics::spacetime::unruh_temperature;
//! use gravitics::strong_field::{E_SCHWINGER, schwinger_exponent};
//! use metrology::{C, E_CHARGE, K_B, M_E};
//!
//! let field = 0.01 * E_SCHWINGER;
//! let t_u = unruh_temperature(E_CHARGE * field / M_E);
//! assert!((schwinger_exponent(field) / (M_E * C * C / (2.0 * K_B * t_u)) - 1.0).abs() < 1e-12);
//! ```
//!
//! **Gravity's vacuum constants.** The gravitational analogues of ε₀ and μ₀, `1/(4πG)` and
//! `4πG/c²`, give `1/√(εμ) = c`: gravitational waves and light travel at one speed.
//!
//! ```
//! use gravitics::gem::{GM_PERMEABILITY, GM_PERMITTIVITY};
//!
//! assert!((1.0 / (GM_PERMITTIVITY * GM_PERMEABILITY).sqrt() / metrology::C - 1.0).abs() < 1e-12);
//! ```
//!
//! **One plasma, three crates.** The plasma frequency and the Alfvén speed appear in
//! `gravitational` (waves meeting a plasma), `magnetics` (MHD) and `plasmoid` (reconnection), which
//! share no dependency but `metrology`, so each keeps its one-line formula. These tests hold the
//! copies together. The electron inertial length is `c` over the plasma's angular frequency,
//! exactly. The reconnection calculator's Alfvén speed counts only the ions' mass, `A·m_p`, so it
//! sits `√(1 + m_e/m_p) − 1` = 2.7×10⁻⁴ above MHD's, which counts the electrons too.
//!
//! ```
//! use gravitational::hf::plasma_frequency;
//! use magnetics::mhd::{alfven_speed, inertial_length, mass_density};
//! use metrology::{C, M_E, M_P};
//! use plasmoid::relaxation::lundquist;
//!
//! let n = 1e19;
//! let skin = inertial_length(n, M_E, 1.0);
//! assert!((skin * 2.0 * core::f64::consts::PI * plasma_frequency(n) / C - 1.0).abs() < 1e-14);
//!
//! let (_, ions_only) = lundquist(0.1, n, 10.0, 0.1, 10.0, 1.0);
//! let with_electrons = alfven_speed(0.1, mass_density(n, M_P, 1.0));
//! let gap = ions_only / with_electrons - 1.0;
//! assert!((gap - ((1.0 + M_E / M_P).sqrt() - 1.0)).abs() < 1e-14 && (gap - 2.72e-4).abs() < 1e-6);
//! ```
//!
//! **Poisson's relation.** A uniformly magnetised body's anomaly is its own gravity gradient. With
//! magnetisation `M` along `M̂` and density contrast `Δρ`, `B = (μ₀/4π)(M/(GΔρ))(Γ − I tr Γ)·M̂`,
//! so at the magnetic pole a total-field survey reads `(μ₀/4π)(M/(GΔρ))Γ_uu`. The sphere in
//! `gravimetry` and the sphere in `magnetics` are written independently and agree on the profile.
//!
//! ```
//! use core::f64::consts::PI;
//! use gravimetry::local::sphere as gravity;
//! use magnetics::anomaly::{sphere as magnetic, total_field};
//! use metrology::{G, MU_0};
//!
//! let (a, depth, rho, m) = (5.0, 20.0, 2700.0, 1.0);
//! for x in [0.0, 7.0, 25.0] {
//!     let (_, gamma_uu) = gravity(a, depth, rho, x);
//!     let b = magnetic([x, 0.0, 0.0], [0.0, 0.0, depth], a, [0.0, 0.0, m]); // z down, at the pole
//!     let poisson = MU_0 / (4.0 * PI) * m / (G * rho) * gamma_uu;
//!     assert!((b[2] / poisson - 1.0).abs() < 1e-14);
//!     // A magnetometer in a 60 000 nT field reads it to within B²/2F.
//!     assert!((total_field([0.0, 0.0, 6e-5], b) - poisson).abs() < 1e-12);
//! }
//! ```
