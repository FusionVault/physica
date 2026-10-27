//! Nuclear accounting: the masses of nuclides, the energy a reaction releases or costs, which
//! isotopes occur in nature, and the bookkeeping that turns a claimed rate of transmutation into
//! the power and radiation it would need.
//!
//! Q-values come from atomic mass excesses (AME2020), so electrons are counted for electron
//! capture and β⁻ and each emitted positron costs `2m_ec²`. A reaction that does not balance its
//! mass number, or its charge after the stated weak conversions, is refused rather than given a
//! number. Isotope ratios, not the appearance of an element, are what separate transmutation from
//! contamination: an α-conjugate route makes pure ⁵⁶Fe, W → Ca + Xe makes ⁴⁸Ca-rich calcium.
//!
//! ```
//! use plasmoid::nuclear::{nuclide, q_value};
//!
//! // ¹²C + ¹⁶O → ²⁸Si releases 16.756 MeV; 4 ¹⁷O → ¹²C + ⁵⁶Fe would release 57.372 MeV.
//! let (c, o16, o17) = ((1, 6, 12), (1, 8, 16), (4, 8, 17));
//! assert!((q_value(&[c, o16], &[(1, 14, 28)], 0, 0).unwrap() - 16.756).abs() < 1e-3);
//! assert!((q_value(&[o17], &[c, (1, 26, 56)], 0, 0).unwrap() - 57.372).abs() < 1e-3);
//! // Copper cannot become lead without paying for it: 4 ⁶³Cu → ²⁰⁸Pb + 11α costs 267.2 MeV.
//! let q = q_value(&[(4, 29, 63)], &[(1, 82, 208), (11, 2, 4)], 0, 12).unwrap();
//! assert!((q + 267.2).abs() < 0.1);
//! assert!(nuclide(8, 16).is_some() && q_value(&[c], &[o16], 0, 0).is_err());
//! ```
//!
//! Sources: Wang et al., [AME2020](https://doi.org/10.1088/1674-1137/abddaf), Chin. Phys. C 45,
//! 030003 (2021); NIST, [Atomic Weights and Isotopic Compositions](https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl).

use metrology::{AMU, C, MEV};

use crate::fusion::{AMU_KEV, ME_KEV};

/// A nuclide's atomic mass excess `Δ = (M − A)u·c²`, keV (AME2020).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Nuclide {
    pub symbol: &'static str,
    pub z: u32,
    pub a: u32,
    pub mass_excess_kev: f64,
}

impl Nuclide {
    /// Atomic mass, u.
    pub fn atomic_mass_u(&self) -> f64 {
        f64::from(self.a) + self.mass_excess_kev / AMU_KEV
    }

    /// Rest energy of the whole atom, MeV: the most a reaction could ever release from it, as
    /// in monopole-catalysed decay (14 899.17 MeV for ¹⁶O, not 16 × 938.272).
    pub fn rest_energy_mev(&self) -> f64 {
        self.atomic_mass_u() * AMU_KEV / 1e3
    }
}

/// The AME2020 nuclides this crate's tests and the topological-monopole research use, from the
/// evaluation's `mass_1.mas20` (all measured, none estimated). Add rows from it as needed.
pub const AME2020: &[Nuclide] = &[
    Nuclide { symbol: "n", z: 0, a: 1, mass_excess_kev: 8071.31806 },
    Nuclide { symbol: "H", z: 1, a: 1, mass_excess_kev: 7288.971064 },
    Nuclide { symbol: "H", z: 1, a: 2, mass_excess_kev: 13135.722895 },
    Nuclide { symbol: "H", z: 1, a: 3, mass_excess_kev: 14949.8109 },
    Nuclide { symbol: "He", z: 2, a: 3, mass_excess_kev: 14931.21888 },
    Nuclide { symbol: "He", z: 2, a: 4, mass_excess_kev: 2424.91587 },
    Nuclide { symbol: "Li", z: 3, a: 6, mass_excess_kev: 14086.88044 },
    Nuclide { symbol: "Li", z: 3, a: 7, mass_excess_kev: 14907.10463 },
    Nuclide { symbol: "Be", z: 4, a: 8, mass_excess_kev: 4941.672 },
    Nuclide { symbol: "B", z: 5, a: 11, mass_excess_kev: 8667.708 },
    Nuclide { symbol: "C", z: 6, a: 11, mass_excess_kev: 10649.397 },
    Nuclide { symbol: "C", z: 6, a: 12, mass_excess_kev: 0.0 },
    Nuclide { symbol: "N", z: 7, a: 14, mass_excess_kev: 2863.41683 },
    Nuclide { symbol: "O", z: 8, a: 16, mass_excess_kev: -4737.00217 },
    Nuclide { symbol: "O", z: 8, a: 17, mass_excess_kev: -808.76421 },
    Nuclide { symbol: "F", z: 9, a: 17, mass_excess_kev: 1951.701 },
    Nuclide { symbol: "Ne", z: 10, a: 20, mass_excess_kev: -7041.93217 },
    Nuclide { symbol: "Ne", z: 10, a: 22, mass_excess_kev: -8024.716 },
    Nuclide { symbol: "Na", z: 11, a: 23, mass_excess_kev: -9529.85352 },
    Nuclide { symbol: "Mg", z: 12, a: 24, mass_excess_kev: -13933.578 },
    Nuclide { symbol: "Al", z: 13, a: 27, mass_excess_kev: -17196.864 },
    Nuclide { symbol: "Si", z: 14, a: 27, mass_excess_kev: -12384.505 },
    Nuclide { symbol: "Si", z: 14, a: 28, mass_excess_kev: -21492.79711 },
    Nuclide { symbol: "Si", z: 14, a: 30, mass_excess_kev: -24432.962 },
    Nuclide { symbol: "P", z: 15, a: 30, mass_excess_kev: -20200.856 },
    Nuclide { symbol: "P", z: 15, a: 31, mass_excess_kev: -24440.54442 },
    Nuclide { symbol: "S", z: 16, a: 31, mass_excess_kev: -19042.531 },
    Nuclide { symbol: "S", z: 16, a: 32, mass_excess_kev: -26015.53714 },
    Nuclide { symbol: "K", z: 19, a: 39, mass_excess_kev: -33807.19535 },
    Nuclide { symbol: "Ca", z: 20, a: 40, mass_excess_kev: -34846.402 },
    Nuclide { symbol: "Ca", z: 20, a: 44, mass_excess_kev: -41468.73 },
    Nuclide { symbol: "Ca", z: 20, a: 46, mass_excess_kev: -43139.61 },
    Nuclide { symbol: "Ca", z: 20, a: 48, mass_excess_kev: -44224.868 },
    Nuclide { symbol: "Fe", z: 26, a: 56, mass_excess_kev: -60607.163 },
    Nuclide { symbol: "Co", z: 27, a: 56, mass_excess_kev: -56040.518 },
    Nuclide { symbol: "Ni", z: 28, a: 56, mass_excess_kev: -53907.649 },
    Nuclide { symbol: "Cu", z: 29, a: 63, mass_excess_kev: -65579.868 },
    Nuclide { symbol: "Cu", z: 29, a: 65, mass_excess_kev: -67263.677 },
    Nuclide { symbol: "Xe", z: 54, a: 132, mass_excess_kev: -89278.97451 },
    Nuclide { symbol: "Xe", z: 54, a: 134, mass_excess_kev: -88125.83443 },
    Nuclide { symbol: "Xe", z: 54, a: 136, mass_excess_kev: -86429.17 },
    Nuclide { symbol: "Xe", z: 54, a: 140, mass_excess_kev: -72986.461 },
    Nuclide { symbol: "Xe", z: 54, a: 144, mass_excess_kev: -56872.301 },
    Nuclide { symbol: "W", z: 74, a: 180, mass_excess_kev: -49636.243 },
    Nuclide { symbol: "W", z: 74, a: 182, mass_excess_kev: -48246.144 },
    Nuclide { symbol: "W", z: 74, a: 184, mass_excess_kev: -45705.453 },
    Nuclide { symbol: "Au", z: 79, a: 197, mass_excess_kev: -31139.751 },
    Nuclide { symbol: "Pb", z: 82, a: 208, mass_excess_kev: -21748.519 },
];

/// The nuclide of charge `z` and mass number `a` in [`AME2020`].
pub fn nuclide(z: u32, a: u32) -> Option<&'static Nuclide> {
    AME2020.iter().find(|n| n.z == z && n.a == a)
}

/// Why a reaction was refused.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Imbalance {
    /// Mass number in minus out.
    MassNumber(i64),
    /// Charge in minus out, less the stated positrons and weak conversions.
    Charge(i64),
    /// A nuclide (z, a) not in [`AME2020`].
    Unknown(u32, u32),
}

/// Q-value (MeV, positive when energy is released) of a reaction whose sides list
/// `(count, z, a)`. `positrons` is the number of β⁺ decays on the way (each costs `2m_ec²` with
/// atomic masses), `weak` the further charge changes by electron capture (positive) or β⁻
/// (negative). The reaction must conserve A, and Z up to those conversions.
pub fn q_value(lhs: &[(u32, u32, u32)], rhs: &[(u32, u32, u32)], positrons: u32, weak: i32) -> Result<f64, Imbalance> {
    let side = |s: &[(u32, u32, u32)]| -> Result<(i64, i64, f64), Imbalance> {
        s.iter().try_fold((0, 0, 0.0), |(a, z, d), &(k, zz, aa)| {
            let n = nuclide(zz, aa).ok_or(Imbalance::Unknown(zz, aa))?;
            Ok((a + i64::from(k * aa), z + i64::from(k * zz), d + f64::from(k) * n.mass_excess_kev))
        })
    };
    let ((a_in, z_in, d_in), (a_out, z_out, d_out)) = (side(lhs)?, side(rhs)?);
    if a_in != a_out {
        return Err(Imbalance::MassNumber(a_in - a_out));
    }
    let dz = z_in - z_out - i64::from(positrons) - i64::from(weak);
    if dz != 0 {
        return Err(Imbalance::Charge(dz));
    }
    Ok((d_in - d_out - 2.0 * ME_KEV * f64::from(positrons)) / 1e3)
}

/// Energy per kilogram of fuel (J/kg) for a reaction releasing `q_mev` from fuel of total mass
/// `fuel_mass_u` (u): `Q/(m u)`. `None` above `c²`, which no reaction can exceed: the
/// "16 × 938 MeV from each ¹⁶O" of monopole catalysis claims is 1.0073 c².
pub fn specific_energy(q_mev: f64, fuel_mass_u: f64) -> Option<f64> {
    let e = q_mev * MEV / (fuel_mass_u * AMU);
    (e <= C * C * (1.0 + 1e-12)).then_some(e)
}

/// Power (W) a source must release to make `kg_per_s` of a product of atomic mass `mass_u` (u)
/// by a reaction of `q_mev` per product atom.
pub fn production_power(kg_per_s: f64, mass_u: f64, q_mev: f64) -> f64 {
    kg_per_s / (mass_u * AMU) * q_mev * MEV
}

/// The most product atoms `energy_j` (J) can make by an endothermic route costing `|q_mev|` each,
/// at perfect efficiency. Its inverse, `N·|Q|`, is the energy a claimed number of conversions
/// needs.
pub fn max_atoms(energy_j: f64, q_mev: f64) -> f64 {
    energy_j / (q_mev.abs() * MEV)
}

/// Natural isotopic compositions `(z, a, fraction)` of the 84 elements with any (NIST); an
/// isotope absent here does not occur in nature.
pub const ABUNDANCES: &[(u32, u32, f64)] = &[
    (1, 1, 0.999885), (1, 2, 0.000115), (2, 3, 1.34e-06), (2, 4, 0.99999866), (3, 6, 0.0759), (3, 7, 0.9241),
    (4, 9, 1.0), (5, 10, 0.199), (5, 11, 0.801), (6, 12, 0.9893), (6, 13, 0.0107), (7, 14, 0.99636),
    (7, 15, 0.00364), (8, 16, 0.99757), (8, 17, 0.00038), (8, 18, 0.00205), (9, 19, 1.0), (10, 20, 0.9048),
    (10, 21, 0.0027), (10, 22, 0.0925), (11, 23, 1.0), (12, 24, 0.7899), (12, 25, 0.1), (12, 26, 0.1101),
    (13, 27, 1.0), (14, 28, 0.92223), (14, 29, 0.04685), (14, 30, 0.03092), (15, 31, 1.0), (16, 32, 0.9499),
    (16, 33, 0.0075), (16, 34, 0.0425), (16, 36, 0.0001), (17, 35, 0.7576), (17, 37, 0.2424), (18, 36, 0.003336),
    (18, 38, 0.000629), (18, 40, 0.996035), (19, 39, 0.932581), (19, 40, 0.000117), (19, 41, 0.067302), (20, 40, 0.96941),
    (20, 42, 0.00647), (20, 43, 0.00135), (20, 44, 0.02086), (20, 46, 4e-05), (20, 48, 0.00187), (21, 45, 1.0),
    (22, 46, 0.0825), (22, 47, 0.0744), (22, 48, 0.7372), (22, 49, 0.0541), (22, 50, 0.0518), (23, 50, 0.0025),
    (23, 51, 0.9975), (24, 50, 0.04345), (24, 52, 0.83789), (24, 53, 0.09501), (24, 54, 0.02365), (25, 55, 1.0),
    (26, 54, 0.05845), (26, 56, 0.91754), (26, 57, 0.02119), (26, 58, 0.00282), (27, 59, 1.0), (28, 58, 0.68077),
    (28, 60, 0.26223), (28, 61, 0.011399), (28, 62, 0.036346), (28, 64, 0.009255), (29, 63, 0.6915), (29, 65, 0.3085),
    (30, 64, 0.4917), (30, 66, 0.2773), (30, 67, 0.0404), (30, 68, 0.1845), (30, 70, 0.0061), (31, 69, 0.60108),
    (31, 71, 0.39892), (32, 70, 0.2057), (32, 72, 0.2745), (32, 73, 0.0775), (32, 74, 0.365), (32, 76, 0.0773),
    (33, 75, 1.0), (34, 74, 0.0089), (34, 76, 0.0937), (34, 77, 0.0763), (34, 78, 0.2377), (34, 80, 0.4961),
    (34, 82, 0.0873), (35, 79, 0.5069), (35, 81, 0.4931), (36, 78, 0.00355), (36, 80, 0.02286), (36, 82, 0.11593),
    (36, 83, 0.115), (36, 84, 0.56987), (36, 86, 0.17279), (37, 85, 0.7217), (37, 87, 0.2783), (38, 84, 0.0056),
    (38, 86, 0.0986), (38, 87, 0.07), (38, 88, 0.8258), (39, 89, 1.0), (40, 90, 0.5145), (40, 91, 0.1122),
    (40, 92, 0.1715), (40, 94, 0.1738), (40, 96, 0.028), (41, 93, 1.0), (42, 92, 0.1453), (42, 94, 0.0915),
    (42, 95, 0.1584), (42, 96, 0.1667), (42, 97, 0.096), (42, 98, 0.2439), (42, 100, 0.0982), (44, 96, 0.0554),
    (44, 98, 0.0187), (44, 99, 0.1276), (44, 100, 0.126), (44, 101, 0.1706), (44, 102, 0.3155), (44, 104, 0.1862),
    (45, 103, 1.0), (46, 102, 0.0102), (46, 104, 0.1114), (46, 105, 0.2233), (46, 106, 0.2733), (46, 108, 0.2646),
    (46, 110, 0.1172), (47, 107, 0.51839), (47, 109, 0.48161), (48, 106, 0.0125), (48, 108, 0.0089), (48, 110, 0.1249),
    (48, 111, 0.128), (48, 112, 0.2413), (48, 113, 0.1222), (48, 114, 0.2873), (48, 116, 0.0749), (49, 113, 0.0429),
    (49, 115, 0.9571), (50, 112, 0.0097), (50, 114, 0.0066), (50, 115, 0.0034), (50, 116, 0.1454), (50, 117, 0.0768),
    (50, 118, 0.2422), (50, 119, 0.0859), (50, 120, 0.3258), (50, 122, 0.0463), (50, 124, 0.0579), (51, 121, 0.5721),
    (51, 123, 0.4279), (52, 120, 0.0009), (52, 122, 0.0255), (52, 123, 0.0089), (52, 124, 0.0474), (52, 125, 0.0707),
    (52, 126, 0.1884), (52, 128, 0.3174), (52, 130, 0.3408), (53, 127, 1.0), (54, 124, 0.000952), (54, 126, 0.00089),
    (54, 128, 0.019102), (54, 129, 0.264006), (54, 130, 0.04071), (54, 131, 0.212324), (54, 132, 0.269086), (54, 134, 0.104357),
    (54, 136, 0.088573), (55, 133, 1.0), (56, 130, 0.00106), (56, 132, 0.00101), (56, 134, 0.02417), (56, 135, 0.06592),
    (56, 136, 0.07854), (56, 137, 0.11232), (56, 138, 0.71698), (57, 138, 0.0008881), (57, 139, 0.9991119), (58, 136, 0.00185),
    (58, 138, 0.00251), (58, 140, 0.8845), (58, 142, 0.11114), (59, 141, 1.0), (60, 142, 0.27152), (60, 143, 0.12174),
    (60, 144, 0.23798), (60, 145, 0.08293), (60, 146, 0.17189), (60, 148, 0.05756), (60, 150, 0.05638), (62, 144, 0.0307),
    (62, 147, 0.1499), (62, 148, 0.1124), (62, 149, 0.1382), (62, 150, 0.0738), (62, 152, 0.2675), (62, 154, 0.2275),
    (63, 151, 0.4781), (63, 153, 0.5219), (64, 152, 0.002), (64, 154, 0.0218), (64, 155, 0.148), (64, 156, 0.2047),
    (64, 157, 0.1565), (64, 158, 0.2484), (64, 160, 0.2186), (65, 159, 1.0), (66, 156, 0.00056), (66, 158, 0.00095),
    (66, 160, 0.02329), (66, 161, 0.18889), (66, 162, 0.25475), (66, 163, 0.24896), (66, 164, 0.2826), (67, 165, 1.0),
    (68, 162, 0.00139), (68, 164, 0.01601), (68, 166, 0.33503), (68, 167, 0.22869), (68, 168, 0.26978), (68, 170, 0.1491),
    (69, 169, 1.0), (70, 168, 0.00123), (70, 170, 0.02982), (70, 171, 0.1409), (70, 172, 0.2168), (70, 173, 0.16103),
    (70, 174, 0.32026), (70, 176, 0.12996), (71, 175, 0.97401), (71, 176, 0.02599), (72, 174, 0.0016), (72, 176, 0.0526),
    (72, 177, 0.186), (72, 178, 0.2728), (72, 179, 0.1362), (72, 180, 0.3508), (73, 180, 0.0001201), (73, 181, 0.9998799),
    (74, 180, 0.0012), (74, 182, 0.265), (74, 183, 0.1431), (74, 184, 0.3064), (74, 186, 0.2843), (75, 185, 0.374),
    (75, 187, 0.626), (76, 184, 0.0002), (76, 186, 0.0159), (76, 187, 0.0196), (76, 188, 0.1324), (76, 189, 0.1615),
    (76, 190, 0.2626), (76, 192, 0.4078), (77, 191, 0.373), (77, 193, 0.627), (78, 190, 0.00012), (78, 192, 0.00782),
    (78, 194, 0.3286), (78, 195, 0.3378), (78, 196, 0.2521), (78, 198, 0.07356), (79, 197, 1.0), (80, 196, 0.0015),
    (80, 198, 0.0997), (80, 199, 0.1687), (80, 200, 0.231), (80, 201, 0.1318), (80, 202, 0.2986), (80, 204, 0.0687),
    (81, 203, 0.2952), (81, 205, 0.7048), (82, 204, 0.014), (82, 206, 0.241), (82, 207, 0.221), (82, 208, 0.524),
    (83, 209, 1.0), (90, 232, 1.0), (91, 231, 1.0), (92, 234, 5.4e-05), (92, 235, 0.007204), (92, 238, 0.992742),
];

/// The natural atom fraction of the isotope (z, a), or `None` if it does not occur in nature.
pub fn abundance(z: u32, a: u32) -> Option<f64> {
    ABUNDANCES.iter().find(|&&(zz, aa, _)| zz == z && aa == a).map(|&(_, _, f)| f)
}

/// The ways a parent (z, a) could split into a nucleus of charge `z1` and its partner with both
/// fragments occurring in nature, as `(a1, a − a1)`. ¹⁸⁴W → Ca + Xe allows only ⁴⁸Ca + ¹³⁶Xe, so
/// calcium made that way would be ⁴⁸Ca-rich, against natural calcium's 96.9% ⁴⁰Ca.
pub fn stable_splits(z: u32, a: u32, z1: u32) -> Vec<(u32, u32)> {
    ABUNDANCES
        .iter()
        .filter(|&&(zz, a1, _)| zz == z1 && a1 < a && abundance(z - z1, a - a1).is_some())
        .map(|&(_, a1, _)| (a1, a - a1))
        .collect()
}

/// How far a contaminant of known composition explains an observed mix of elements: with
/// `observed` and `source` in the same element order (any units), it scales the source to the
/// first element and returns `(k, ratios)`, the scale and each element's observed/expected. Ratios
/// near 1 point to the contaminant; the method of Priem et al. (2008) for the iron, chromium and
/// nickel left by exploded wires.
pub fn contamination_ratios(observed: &[f64], source: &[f64]) -> (f64, Vec<f64>) {
    let k = observed[0] / source[0];
    (k, observed.iter().zip(source).map(|(o, s)| o / (k * s)).collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, rel: f64) -> bool {
        (a / b - 1.0).abs() < rel
    }

    fn q(lhs: &[(u32, u32, u32)], rhs: &[(u32, u32, u32)]) -> f64 {
        q_value(lhs, rhs, 0, 0).unwrap()
    }

    const N: (u32, u32, u32) = (1, 0, 1);
    const H: (u32, u32, u32) = (1, 1, 1);
    const D: (u32, u32, u32) = (1, 1, 2);
    const T: (u32, u32, u32) = (1, 1, 3);
    const HE3: (u32, u32, u32) = (1, 2, 3);
    const HE4: (u32, u32, u32) = (1, 2, 4);
    const C12: (u32, u32, u32) = (1, 6, 12);
    const O16: (u32, u32, u32) = (1, 8, 16);

    #[test]
    fn the_table() {
        assert_eq!(AME2020.len(), 48);
        for (i, n) in AME2020.iter().enumerate() {
            assert!(AME2020[i + 1..].iter().all(|m| (m.z, m.a) != (n.z, n.a)), "{} {}", n.symbol, n.a);
        }
        // Mass excess and AME2020's printed atomic masses agree to 10⁻⁸ u (0.01 keV), through
        // CODATA 2018's u = 931 494.102 42 keV, the value AME2020 uses.
        for (z, a, m) in [(1, 1, 1.007_825_031_898), (8, 16, 15.994_914_619_26), (26, 56, 55.934_935_537), (82, 208, 207.976_652_005), (0, 1, 1.008_664_915_90)] {
            assert!((nuclide(z, a).unwrap().atomic_mass_u() - m).abs() < 1e-8, "{z} {a}");
        }
        assert_eq!(nuclide(6, 12).unwrap().mass_excess_kev, 0.0);
        assert!(nuclide(92, 238).is_none());
    }

    #[test]
    fn published_q_values() {
        // Bosch & Hale's D–D branches (fusion.rs 3.268, 4.033 MeV); oxygen burning (Woosley,
        // Heger & Weaver 2002): ²⁸Si+α 9.593, ³¹P+p 7.676, ³²S 16.539; carbon–oxygen 16.76.
        assert!((q(&[(2, 1, 2)], &[HE3, N]) - 3.268).abs() < 1e-3 && (q(&[(2, 1, 2)], &[T, H]) - 4.033).abs() < 1e-3);
        let oo = [(2, 8, 16)];
        assert!((q(&oo, &[(1, 14, 28), HE4]) - 9.593).abs() < 3e-3 && (q(&oo, &[(1, 15, 31), H]) - 7.676).abs() < 3e-3);
        assert!((q(&oo, &[(1, 16, 32)]) - 16.539).abs() < 3e-3 && (q(&[C12, O16], &[(1, 14, 28)]) - 16.76).abs() < 1e-2);
        // The aneutronic constants are AME2020's to their last digit.
        use crate::aneutronic::{Q_DDN, Q_DDP, Q_DHE3, Q_DT, Q_HE3HE3, Q_PB11, Q_PLI6};
        for (got, want) in [
            (q(&[H, (1, 5, 11)], &[(3, 2, 4)]), Q_PB11),
            (q(&[D, HE3], &[HE4, H]), Q_DHE3),
            (q(&[(2, 2, 3)], &[HE4, (2, 1, 1)]), Q_HE3HE3),
            (q(&[H, (1, 3, 6)], &[HE4, HE3]), Q_PLI6),
            (q(&[D, T], &[HE4, N]), Q_DT),
            (q(&[(2, 1, 2)], &[HE3, N]), Q_DDN),
            (q(&[(2, 1, 2)], &[T, H]), Q_DDP),
        ] {
            assert!((got - want).abs() < 5e-5, "{got} {want}");
        }
        // And fusion.rs's nuclear masses are the atomic masses less electrons and their binding
        // (13.6 eV, 79.0 eV), to 0.3 eV.
        use crate::fusion::{M_D, M_HE3, M_T};
        let nuclear = |z: u32, a: u32, binding_kev: f64| nuclide(z, a).unwrap().atomic_mass_u() * AMU_KEV - f64::from(z) * ME_KEV + binding_kev;
        assert!((nuclear(1, 2, 13.6e-3) - M_D).abs() < 3e-4 && (nuclear(1, 3, 13.6e-3) - M_T).abs() < 4e-4 && (nuclear(2, 3, 79.005e-3) - M_HE3).abs() < 3e-4);
        // The plasmoids report's other rows (NIST masses, so to 0.2 keV): p + ⁷Li → 2α 17.3462,
        // D + ⁶Li → 2α 22.3728, α + ¹¹B → ¹⁴N + n 0.1581.
        assert!((q(&[H, (1, 3, 7)], &[(2, 2, 4)]) - 17.3462).abs() < 3e-4 && (q(&[D, (1, 3, 6)], &[(2, 2, 4)]) - 22.3728).abs() < 3e-4);
        assert!((q(&[HE4, (1, 5, 11)], &[(1, 7, 14), N]) - 0.1581).abs() < 3e-4);
        // p + ¹¹B → ¹¹C + n, the neutron channel: −2.7640 MeV.
        assert!((q(&[H, (1, 5, 11)], &[(1, 6, 11), N]) + 2.7640).abs() < 1e-4);
    }

    #[test]
    fn the_talk_reactions() {
        // 2 ²⁸Si → ⁵⁶Ni 10.922; then two electron captures to ⁵⁶Fe, 6.700; the same with two
        // positrons instead, 15.578 kinetic.
        assert!((q(&[(2, 14, 28)], &[(1, 28, 56)]) - 10.922).abs() < 1e-3);
        assert!((q_value(&[(1, 28, 56)], &[(1, 26, 56)], 0, 2).unwrap() - 6.700).abs() < 1e-3);
        assert!((q_value(&[(2, 14, 28)], &[(1, 26, 56)], 2, 0).unwrap() - 15.578).abs() < 1e-3);
        // 4 ¹⁷O → ¹²C + ⁵⁶Fe 57.372 (the audit's 49.47 is wrong); 2 ¹⁷O → ¹²C + ²²Ne 6.407 (not
        // 17.13); ²³Na + ¹⁶O → ³⁹K 19.540; ¹⁸⁴W → ⁴⁸Ca + ¹³⁶Xe 84.949.
        assert!((q(&[(4, 8, 17)], &[C12, (1, 26, 56)]) - 57.372).abs() < 1e-3 && (q(&[(2, 8, 17)], &[C12, (1, 10, 22)]) - 6.407).abs() < 1e-3);
        assert!((q(&[(1, 11, 23), O16], &[(1, 19, 39)]) - 19.540).abs() < 1e-3 && (q(&[(1, 74, 184)], &[(1, 20, 48), (1, 54, 136)]) - 84.949).abs() < 1e-3);
        // Endothermic: 4 ⁶³Cu → ²⁰⁸Pb + 11α (12 captures) −267.245; ¹⁸⁴W → 46α (18 β⁻) −157.252.
        assert!((q_value(&[(4, 29, 63)], &[(1, 82, 208), (11, 2, 4)], 0, 12).unwrap() + 267.245).abs() < 1e-3);
        assert!((q_value(&[(1, 74, 184)], &[(46, 2, 4)], 0, -18).unwrap() + 157.252).abs() < 1e-3);
        // Unbalanced reactions are refused.
        assert_eq!(q_value(&[C12], &[O16], 0, 0), Err(Imbalance::MassNumber(-4)));
        assert_eq!(q_value(&[(1, 28, 56)], &[(1, 26, 56)], 0, 0), Err(Imbalance::Charge(2)));
        assert_eq!(q_value(&[(1, 92, 238)], &[(1, 92, 238)], 0, 0), Err(Imbalance::Unknown(92, 238)));
    }

    #[test]
    fn energy_bookkeeping() {
        // The whole ¹⁶O atom is 14 899.17 MeV, and that is c² per kilogram; 16 m_p is 1.0073 c², so
        // it is refused.
        let o16 = nuclide(8, 16).unwrap();
        assert!((o16.rest_energy_mev() - 14_899.17).abs() < 0.01);
        assert!(close(specific_energy(o16.rest_energy_mev(), o16.atomic_mass_u()).unwrap(), C * C, 1e-12));
        assert!(specific_energy(16.0 * 938.272_088, o16.atomic_mass_u()).is_none());
        // 4 ¹⁷O → Fe + C: 8.14×10¹³ J/kg.
        let fuel = 4.0 * nuclide(8, 17).unwrap().atomic_mass_u();
        assert!(close(specific_energy(57.372, fuel).unwrap(), 8.141e13, 1e-3));
        // Srinivasan's "4 tonnes of iron a day" made as ⁵⁶Fe: 8.722×10¹¹ W via Si+Si and
        // 4.083×10¹² W from C+O, against an 11 MW furnace.
        let (rate, fe) = (4000.0 / 86_400.0, nuclide(26, 56).unwrap().atomic_mass_u());
        assert!(close(production_power(rate, fe, 10.922_055), 8.722e11, 1e-3));
        assert!(close(production_power(rate, fe, 51.133_159), 4.083e12, 1e-3));
        // Adamenko: 300 J makes at most 7.0065×10¹² lead atoms (2.42 ng) at 267.245 MeV each.
        let n = max_atoms(300.0, -267.245);
        assert!(close(n, 7.0065e12, 1e-4) && close(n * 207.98 * AMU * 1e12, 2.42, 1e-3));
    }

    #[test]
    fn isotopes_in_nature() {
        // Every element's fractions sum to 1.
        for z in 1..=92 {
            let s: f64 = ABUNDANCES.iter().filter(|r| r.0 == z).map(|r| r.2).sum();
            assert!(s == 0.0 || (s - 1.0).abs() < 1e-9, "{z}: {s}");
        }
        assert_eq!(ABUNDANCES.len(), 288);
        assert_eq!(abundance(22, 48), Some(0.7372));
        assert_eq!(abundance(20, 40), Some(0.969_41));
        assert_eq!(abundance(8, 17), Some(0.000_38));
        assert!(abundance(28, 56).is_none());
        // ¹⁸⁴W can split into Ca + Xe only as ⁴⁸Ca + ¹³⁶Xe; ¹⁸³W and ¹⁸⁶W not at all.
        assert_eq!(stable_splits(74, 184, 20), vec![(48, 136)]);
        assert!(stable_splits(74, 183, 20).is_empty() && stable_splits(74, 186, 20).is_empty());
        // Priem et al. 2008: 18-8 steel (Fe 72, Cr 20, Ni 8) explains shot 2's Fe 5.16, Cr 1.1,
        // Ni 0.52 and shot 3's 6.6, 1.3, 0.67 at 0.71–0.91 of the expected Cr and Ni.
        let steel = [72.0, 20.0, 8.0];
        let (k2, r2) = contamination_ratios(&[5.16, 1.1, 0.52], &steel);
        let (k3, r3) = contamination_ratios(&[6.6, 1.3, 0.67], &steel);
        assert!(close(k2, 0.071_67, 1e-3) && close(r2[1], 0.767, 1e-3) && close(r2[2], 0.907, 1e-3) && r2[0] == 1.0);
        assert!(close(k3, 0.091_67, 1e-3) && close(r3[1], 0.709, 1e-3) && close(r3[2], 0.914, 1e-3));
    }
}
