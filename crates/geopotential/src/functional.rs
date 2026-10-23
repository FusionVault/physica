/// Mean density of the Earth used by Wahr et al. (1998), kg/m³. gravity-toolkit instead derives
/// it from `GM` (≈5.51e3), which changes water heights by about 0.1%.
pub const RHO_EARTH: f64 = 5517.0;

/// Density of fresh water, kg/m³.
pub const RHO_WATER: f64 = 1000.0;

/// A quantity that is a per-degree multiple of the harmonic sum
/// `Σ_m P̄nm(sin φ)(C̄nm cos mλ + S̄nm sin mλ)` (spherical approximation).
///
/// For the geoid, gravity disturbance, anomaly and gradient, pass disturbing-potential
/// coefficients (model minus normal field) or a difference of two fields such as one GRACE month
/// minus a mean.
#[derive(Clone, Copy, Debug)]
pub enum Functional<'a> {
    /// Gravitational potential `V`, m²/s².
    Potential,
    /// Geoid height / height anomaly by Bruns' formula with `γ = GM/r²`, m.
    GeoidHeight,
    /// Gravity disturbance `−∂T/∂r`, m/s² (×1e5 for mGal).
    GravityDisturbance,
    /// Gravity anomaly `−∂T/∂r − 2T/r`, m/s².
    GravityAnomaly,
    /// Radial gravity gradient `∂²T/∂r²`, s⁻² (×1e9 for Eötvös).
    RadialGradient,
    /// Equivalent water height on the reference sphere (Wahr et al. 1998, eq. 14 divided by
    /// `ρ_w`), m. `love` holds the load Love numbers `k'_n`, indexed by degree; `r` is ignored.
    WaterHeight {
        /// Load Love numbers `k'_n` for every degree used.
        love: &'a [f64],
    },
}

impl Functional<'_> {
    /// The factor `κ_n` multiplying the degree-`n` harmonic sum, at radius `r`, for a field with
    /// the given `gm` and reference `radius`.
    ///
    /// # Panics
    /// For [`Functional::WaterHeight`], if `love` has no entry for degree `n`.
    pub fn factor(&self, n: usize, gm: f64, radius: f64, r: f64) -> f64 {
        let q = (radius / r).powi(n as i32);
        let nf = n as f64;
        match self {
            Functional::Potential => gm / r * q,
            Functional::GeoidHeight => r * q,
            Functional::GravityDisturbance => gm / (r * r) * (nf + 1.0) * q,
            Functional::GravityAnomaly => gm / (r * r) * (nf - 1.0) * q,
            Functional::RadialGradient => gm / (r * r * r) * (nf + 1.0) * (nf + 2.0) * q,
            Functional::WaterHeight { love } => {
                radius * RHO_EARTH / (3.0 * RHO_WATER) * (2.0 * nf + 1.0) / (1.0 + love[n])
            }
        }
    }
}
