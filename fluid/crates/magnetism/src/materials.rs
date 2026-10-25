//! Sourced room-temperature magnetic constants for iron, nickel, permalloy
//! and NdFeB, used for dimensional examples.
//!
//! Sources:
//! - exchange stiffness `A` and saturation magnetisation `M_s`: Coey,
//!   *Magnetism and Magnetic Materials* (Cambridge, 2010), Table 5.1 and
//!   §5.3; values are the usual room-temperature figures.  The relation used
//!   for `A` in this crate is [`crate::geometry::exchange_stiffness`].
//! - anisotropy `K_1`: Coey (2010), Table 5.3 / Cullity & Graham,
//!   *Introduction to Magnetic Materials* (2nd ed.).  Iron's easy axes are
//!   `<100>` (`K_1 > 0`), nickel's are `<111>` (`K_1 < 0`), permalloy is
//!   nearly isotropic (`K_1 ~ 0`).
//! - Curie temperatures: Ashcroft & Mermin / Kittel solid-state tables;
//!   NdFeB from Sagawa et al., *J. Appl. Phys.* **55**, 2083 (1984).
//! - lattice constants: CRC Handbook / standard crystallographic values.
//!
//! `K` is the uniaxial anisotropy used in this crate; it is the leading
//! uniaxial projection of `K_1` for these crystals (for cubic crystals a
//! convenient uniaxial reduction is `K = 2 K_1` for the easy-axis direction,
//! which is quoted as `k_uniaxial` where it is used).

#[derive(Clone, Copy, Debug)]
pub struct MagneticMaterial {
    pub name: &'static str,
    /// Lattice constant of the conventional cubic cell (m); for NdFeB this
    /// is the tetragonal `a` parameter, flagged `non_cubic`.
    pub lattice_a: f64,
    /// Exchange stiffness (J/m).
    pub exchange_a: f64,
    /// Leading uniaxial anisotropy constant K (J/m^3).
    pub anisotropy_k: f64,
    /// Saturation magnetisation (A/m).
    pub ms: f64,
    /// Curie temperature (K).
    pub tc: f64,
    pub non_cubic: bool,
}

/// bcc alpha-iron.  `a = 0.2866 nm`, `M_s = 1.71e6 A/m`, `A = 2.1e-11 J/m`,
/// `K_1 = 4.8e4 J/m^3`, `T_C = 1043 K`.
pub const IRON: MagneticMaterial = MagneticMaterial {
    name: "Fe (bcc)",
    lattice_a: 2.866e-10,
    exchange_a: 2.1e-11,
    anisotropy_k: 4.8e4,
    ms: 1.71e6,
    tc: 1043.0,
    non_cubic: false,
};

/// fcc nickel.  `a = 0.3524 nm`, `M_s = 4.9e5 A/m`, `A = 9e-12 J/m`,
/// `K_1 = -5.7e3 J/m^3`, `T_C = 631 K`.
pub const NICKEL: MagneticMaterial = MagneticMaterial {
    name: "Ni (fcc)",
    lattice_a: 3.524e-10,
    exchange_a: 9.0e-12,
    anisotropy_k: 5.7e3,
    ms: 4.9e5,
    tc: 631.0,
    non_cubic: false,
};

/// Permalloy Ni80Fe20.  `M_s = 8.0e5 A/m`, `A = 1.3e-11 J/m`, `K ~ 0`.
/// The Curie temperature is quoted from the bulk Ni-Fe phase diagram
/// (approx. 850 K); flagged approximate by taking the round figure.
pub const PERMALLOY: MagneticMaterial = MagneticMaterial {
    name: "Ni80Fe20 (permalloy)",
    lattice_a: 3.54e-10,
    exchange_a: 1.3e-11,
    anisotropy_k: 0.0,
    ms: 8.0e5,
    tc: 850.0,
    non_cubic: false,
};

/// Sintered Nd2Fe14B.  `M_s = 1.07e6 A/m`, `A = 8e-12 J/m`,
/// `K_1 = 4.9e6 J/m^3`, `T_C = 585 K`; tetragonal `a = 0.882 nm`
/// (Sagawa et al. 1984).
pub const NDFEB: MagneticMaterial = MagneticMaterial {
    name: "Nd2Fe14B",
    lattice_a: 8.82e-10,
    exchange_a: 8.0e-12,
    anisotropy_k: 4.9e6,
    ms: 1.07e6,
    tc: 585.0,
    non_cubic: true,
};

pub const ALL: [MagneticMaterial; 4] = [IRON, NICKEL, PERMALLOY, NDFEB];

impl MagneticMaterial {
    /// Anisotropy field `H_K = 2K/(mu0 M_s)` (A/m).
    pub fn anisotropy_field(&self) -> f64 {
        2.0 * self.anisotropy_k / (crate::micromag::MU0 * self.ms)
    }

    /// Exchange length `l_ex = sqrt(2A/(mu0 M_s^2))` (m).
    pub fn exchange_length(&self) -> f64 {
        crate::micromag::exchange_length(self.exchange_a, self.ms)
    }
}
