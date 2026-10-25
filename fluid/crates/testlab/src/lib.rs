//! Virtual materials-testing lab.
//!
//! For one material at every emergent level — atomistic, single-crystal,
//! polycrystal, continuum — this crate runs the standard material tests
//! (static tensile, ultrasonic, atomistic-to-continuum fold, and a 1D
//! plate-impact dynamic test) and records the measured value, a published
//! or exact reference, the error and the level. Where the fold between two
//! levels is not valid, the summary says so rather than hiding it.
//!
//! The organising idea is [the emergence ladder](../../../analysis/emergence_ladder.md):
//! a cheap model, a cheap validity indicator, and an unfold to a finer model
//! where the indicator fails. Here the rungs are
//!
//! | level | model | validity indicator |
//! |---|---|---|
//! | atomistic | LJ pair potential / measured single-crystal `C_ij` | — |
//! | single crystal | Cauchy–Born or measured anisotropic elasticity | `(ka)²` dispersion error |
//! | polycrystal | Voigt/Reuss/Hill isotropic average | grain size / wavelength `d/λ` |
//! | continuum shock | GPR elastic–plastic solid + NASG EOS | Hugoniot `Us(up)` and `σ_HEL` |
//!
//! Everything is SI (m, kg, s, Pa). See `analysis/testlab.md` for the
//! measured numbers and the honest misses.

pub mod atomistic;
pub mod dynamic;
pub mod elastic;
pub mod interfaces;
pub mod magnetic;
pub mod mechanical;
pub mod optical;
pub mod report;
pub mod thermal;
pub mod viscoelastic;

use materials::anisotropic::{cubic, M6};

/// The emergent level a measurement was made at.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Level {
    Atomistic,
    SingleCrystal,
    Polycrystal,
    Continuum,
    Macroscopic,
}

impl Level {
    pub fn label(self) -> &'static str {
        match self {
            Level::Atomistic => "atomistic",
            Level::SingleCrystal => "single-crystal",
            Level::Polycrystal => "polycrystal",
            Level::Continuum => "continuum",
            Level::Macroscopic => "macroscopic",
        }
    }
}

/// One standard-test measurement.
#[derive(Clone, Debug)]
pub struct Measurement {
    pub material: String,
    pub test: String,
    pub level: Level,
    pub quantity: String,
    pub value: f64,
    pub reference: f64,
    pub unit: &'static str,
    /// True when the fold to this level is certified for this test, false
    /// when the validity indicator says it must be unfolded.
    pub fold_valid: bool,
    pub note: String,
}

impl Measurement {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        material: &str,
        test: &str,
        level: Level,
        quantity: &str,
        value: f64,
        reference: f64,
        unit: &'static str,
        fold_valid: bool,
        note: &str,
    ) -> Self {
        Measurement {
            material: material.to_string(),
            test: test.to_string(),
            level,
            quantity: quantity.to_string(),
            value,
            reference,
            unit,
            fold_valid,
            note: note.to_string(),
        }
    }

    /// Signed error `value − reference`.
    pub fn error(&self) -> f64 {
        self.value - self.reference
    }

    /// Relative error, zero reference guarded.
    pub fn rel_error(&self) -> f64 {
        let s = self.reference.abs().max(1e-300);
        (self.value - self.reference) / s
    }
}

/// The atomistic model behind a material. For a real metal the honest
/// statement is that a central pair potential is inadequate; the bottom
/// rung is then the measured single-crystal stiffness.
#[derive(Clone, Copy, Debug)]
pub enum BottomRung {
    /// A fully self-consistent pair-potential crystal (LJ argon).
    PairPotential {
        potential: &'static str,
        sigma_m: f64,
        epsilon_j: f64,
        mass_kg: f64,
        /// Equilibrium fcc lattice constant at zero pressure, m.
        a0_m: f64,
        note: &'static str,
    },
    /// The measured single-crystal constants are the bottom rung; a pair
    /// potential would not reproduce them.
    MeasuredSingleCrystal { note: &'static str },
}

/// A material description gathering what each emergent level needs.
#[derive(Clone, Debug)]
pub struct Material {
    pub name: &'static str,
    /// Mass density, kg/m³.
    pub density: f64,
    /// Single-crystal Voigt stiffness in the crystal frame, Pa.
    pub c: M6,
    /// Cubic constants `(C11, C12, C44)`, Pa (convenience).
    pub c11: f64,
    pub c12: f64,
    pub c44: f64,
    pub bottom: BottomRung,
    /// Handbook polycrystal Young's modulus and Poisson ratio.
    pub young_ref: f64,
    pub poisson_ref: f64,
    pub elastic_source: &'static str,
    /// Polycrystal bulk yield stress, Pa.
    pub yield_stress: f64,
    pub yield_source: &'static str,
    /// Shock Hugoniot `Us = c0 + s·up`, c0 in m/s.
    pub hugoniot_c0: f64,
    pub hugoniot_s: f64,
    pub hugoniot_source: &'static str,
}

impl Material {
    /// Single-crystal Young's modulus along a crystal direction.
    pub fn young(&self, n: [f64; 3]) -> f64 {
        materials::anisotropic::youngs_modulus(&self.c, n)
    }

    /// Zener anisotropy ratio `2C44/(C11−C12)`.
    pub fn zener(&self) -> f64 {
        materials::anisotropic::zener_ratio(self.c11, self.c12, self.c44)
    }
}

/// Copper: measured single-crystal constants (Simmons & Wang 1971).
pub fn copper() -> Material {
    let (c11, c12, c44) = (168.4e9, 121.4e9, 75.4e9);
    Material {
        name: "Cu",
        density: 8960.0,
        c: cubic(c11, c12, c44),
        c11,
        c12,
        c44,
        bottom: BottomRung::MeasuredSingleCrystal {
            note: "no central pair potential reproduces copper's d-band \
                   bonding or positive stacking-fault energy; the measured \
                   single-crystal constants are the bottom rung",
        },
        young_ref: 128.0e9,
        poisson_ref: 0.34,
        elastic_source: "single crystal: Simmons & Wang 1971; polycrystal \
                         E, ν: Ashby & Jones, Engineering Materials / CRC",
        yield_stress: 70.0e6,
        yield_source: "annealed polycrystal Cu, Callister, Materials Science \
                       and Engineering, Table 6.1 (order of magnitude; temper \
                       dependent)",
        hugoniot_c0: 3940.0,
        hugoniot_s: 1.49,
        hugoniot_source: "Marsh, LASL Shock Hugoniot Data (1980), Cu fit",
    }
}

/// Aluminium: measured single-crystal constants (Simmons & Wang 1971).
pub fn aluminium() -> Material {
    let (c11, c12, c44) = (107.3e9, 60.9e9, 28.3e9);
    Material {
        name: "Al",
        density: 2700.0,
        c: cubic(c11, c12, c44),
        c11,
        c12,
        c44,
        bottom: BottomRung::MeasuredSingleCrystal {
            note: "as for Cu: a pair potential misses the metallic bonding; \
                   the measured single-crystal constants are the bottom rung",
        },
        young_ref: 70.0e9,
        poisson_ref: 0.35,
        elastic_source: "single crystal: Simmons & Wang 1971; polycrystal \
                         E, ν: Ashby & Jones / CRC",
        yield_stress: 35.0e6,
        yield_source: "annealed 1100 Al, Callister, Materials Science and \
                       Engineering, Table 6.1 (order of magnitude)",
        hugoniot_c0: 5350.0,
        hugoniot_s: 1.34,
        hugoniot_source: "Marsh, LASL Shock Hugoniot Data (1980), Al fit",
    }
}

/// Lennard-Jones argon: the fully self-consistent case, where the atomistic,
/// single-crystal and continuum rungs are all generated by the same pair
/// potential.
pub fn argon() -> Material {
    let (c, density, a0) = atomistic::argon_cauchy_born();
    // Cubic constants in SI from the reduced Cauchy–Born tensor.
    let c11 = c[0][0];
    let c12 = c[0][1];
    let c44 = c[3][3];
    Material {
        name: "Ar",
        density,
        c,
        c11,
        c12,
        c44,
        bottom: BottomRung::PairPotential {
            potential: "12-6 Lennard-Jones",
            sigma_m: md::argon::SIGMA_M,
            epsilon_j: md::argon::epsilon_j(),
            mass_kg: md::argon::mass_kg(),
            a0_m: a0,
            note: "the only fully self-consistent material here: the same pair \
                   potential gives the Cauchy–Born stiffness, the phonons, \
                   the MD state point and (as a closure) the continuum; \
                   central forces give the Cauchy relation C12 = C44",
        },
        young_ref: f64::NAN,
        poisson_ref: f64::NAN,
        elastic_source: "Cauchy–Born from the LJ pair potential (no handbook \
                         comparison; quantum solid argon is not a classical LJ \
                         crystal)",
        yield_stress: f64::NAN,
        yield_source: "not tested dynamically",
        hugoniot_c0: f64::NAN,
        hugoniot_s: f64::NAN,
        hugoniot_source: "not tested dynamically",
    }
}

/// The three materials of the report.
pub fn materials() -> Vec<Material> {
    vec![copper(), aluminium(), argon()]
}

/// Convenience: look up a material by name.
pub fn material(name: &str) -> Option<Material> {
    materials().into_iter().find(|m| m.name == name)
}
