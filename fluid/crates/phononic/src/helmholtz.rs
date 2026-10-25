//! Acoustic metamaterials from fluid-filled resonators: Helmholtz resonators,
//! side-branch-loaded waveguides and their effective-medium (folded) limit.
//!
//! The rung above the plain layered/duct description is the **lumped**
//! Helmholtz resonator. A cavity of volume `V` connected through a neck of
//! cross-section `S` and effective length `L_eff` is a mass–spring oscillator
//! with acoustic mass `M_a = rho L_eff / S` and acoustic compliance
//! `C_a = V / (rho c^2)`, resonating at
//!
//! ```text
//!   omega0 = c sqrt(S / (V L_eff)).
//! ```
//!
//! `L_eff = L + delta_out + delta_in` accounts for the fluid that moves outside
//! the geometric neck. For a circular neck of radius `a` the low-frequency end
//! corrections are `delta = (8/(3 pi)) a = 0.8488 a` for a flanged opening and
//! `delta = 0.6133 a` for an unflanged (free-space) opening (Rayleigh, *Theory
//! of Sound* II; Levine & Schwinger, J. Acoust. Soc. Am. **20**, 637 (1948);
//! Kinsler & Frey, *Fundamentals of Acoustics*, 4th ed., ch. 9). The flanged
//! value is the low-`ka` reactance of a baffled piston,
//! `X_rad = (8/(3 pi)) ka rho c`; the unflanged value is the Rayleigh open-pipe
//! correction.
//!
//! This module keeps both the cheap lumped model and the exact one, on the same
//! geometry:
//!
//! * **Lumped (fold):** [`Helmholtz::lumped_impedance`], a series mass and
//!   compliance. With thermoviscous losses it becomes
//!   [`Helmholtz::lumped_impedance_lossy`], using the Stinson low-reduced-
//!   frequency (LRF) narrow-tube model.
//! * **Exact (unfold):** [`Helmholtz::input_impedance_exact`] concatenates the
//!   neck and cavity as transmission-line ducts and terminates the cavity
//!   rigidly. The validity indicator is the neck and cavity length over the
//!   wavelength, `L/lambda`; the lumped model fails as `L` or `L_c` approaches
//!   `lambda/8`–`lambda/4`, where the cavity stops being a pure compliance.
//!
//! # Thermoviscous losses
//!
//! The narrow-tube model of Zwikker & Kosten (1949) and Stinson (J. Acoust.
//! Soc. Am. **89**, 550 (1991)) gives the complex density and compressibility
//!
//! ```text
//!   rho_c = rho0 / [1 - 2 J1(alpha a)/(alpha a J0(alpha a))],
//!   C_c   = (1/B0)[1 + (gamma-1) 2 J1(alpha_t a)/(alpha_t a J0(alpha_t a))],
//!   alpha^2 = -i omega rho0 / mu,   alpha_t = sqrt(Pr) alpha,
//! ```
//!
//! with `B0 = rho0 c0^2`, `Pr` the Prandtl number and `mu` the shear
//! viscosity. The characteristic impedance is `sqrt(rho_c/C_c)` and the
//! wavenumber is `omega sqrt(rho_c C_c)`. In the low-frequency (Poiseuille)
//! limit the neck series impedance tends to `8 mu L/(pi a^4) + i omega rho0
//! L/(pi a^2)`, which the tests check.
//!
//! # Side branches and the negative-bulk-modulus band
//!
//! A side branch of admittance `Y_b = 1/Z_b` on a waveguide of impedance `Z0`
//! transmits `t = 2 Z_b / (2 Z_b + Z0)` (pressure amplitude). A periodic array
//! of such branches, spacing `a`, has the exact Bloch relation
//! `cos(q a) = tr(M_cell)/2` with `M_cell` the duct segment times the shunt.
//! At long wavelength it folds to an effective medium with bulk modulus
//!
//! ```text
//!   1/B_eff = (1/B0) [1 - F omega0^2 / (omega^2 - omega0^2 + i Gamma omega)],
//!   F = V / (A a),
//! ```
//!
//! which is the Fang et al. form (Nature Materials **5**, 452 (2006)) with the
//! filling factor `F` fixed by the geometry (cavity volume over cell volume).
//! `B_eff < 0` for `omega` just above `omega0`: the negative-bulk-modulus band.
//! The folded model is valid while `k a = omega a/c` is small; the tests measure
//! the error against the exact Bloch wavevector and record where it grows.
//!
//! Filling the same cavity with water instead of air raises the fluid bulk
//! modulus by four orders of magnitude. Wall compliance then matters: a
//! spherical cavity in a solid of shear modulus `G` has an extra compliance
//! `3/(4G)` per unit volume (Lamé solution for a pressurised spherical cavity),
//! so `1/B_cav = 1/B0 + 3/(4G)`. For water in steel the wall term is ~2%,
//! for water in a polymer it is comparable to the fluid term.
//!
//! # The fold/unfold map
//!
//! | coarse model (fold) | indicator | unfold |
//! |---|---|---|
//! | lumped Helmholtz mass+compliance | `L/lambda`, `L_c/lambda` | exact duct transfer matrix |
//! | lossless lumped model | `s = a sqrt(omega rho/mu)`, `ka` | Stinson LRF, then unsteady |
//! | effective bulk modulus `B_eff` | `k a` | exact Bloch of the periodic array |
//! | rigid cavity compliance | `3 B0/(4 G)` | elastic-wall cavity compliance |

use crate::complex::C;
use crate::layered::Mat2;
use crate::materials::{Fluid, Solid};
use std::f64::consts::PI;

/// Dry air at 20 C, 1 atm. Density and sound speed from the CRC Handbook
/// (rho = 1.204 kg m^-3, c = 343.2 m s^-1); the rounded 343 m/s is used here.
pub const AIR: Fluid = Fluid {
    name: "air",
    rho: 1.204,
    c: 343.0,
    source: "CRC Handbook, dry air 20 C, 1 atm (rho=1.204, c=343.2)",
};

/// Transport properties needed by the thermoviscous narrow-tube model.
#[derive(Clone, Copy, Debug)]
pub struct Thermoviscous {
    /// Dynamic (shear) viscosity [Pa s].
    pub mu: f64,
    /// Prandtl number `mu c_p / kappa`.
    pub pr: f64,
    /// Ratio of specific heats, `c_p/c_v`; enters only the thermal term
    /// `(gamma - 1)`.
    pub gamma: f64,
}

/// Air transport at 20 C: mu = 1.825e-5 Pa s, Pr = 0.71, gamma = 1.4
/// (CRC Handbook / Kinsler & Frey, *Fundamentals of Acoustics*).
pub const AIR_TRANSPORT: Thermoviscous = Thermoviscous {
    mu: 1.825e-5,
    pr: 0.71,
    gamma: 1.4,
};

/// Water transport at 20 C: mu = 1.002e-3 Pa s, Pr = 7.0, gamma = 1.0
/// (CRC Handbook; gamma is 1 for a nearly incompressible liquid, so the
/// thermal term drops out).
pub const WATER_TRANSPORT: Thermoviscous = Thermoviscous {
    mu: 1.002e-3,
    pr: 7.0,
    gamma: 1.0,
};

// ---------------------------------------------------------------------------
// Complex Bessel functions J0, J1 (power series, with a large-argument limit)
// ---------------------------------------------------------------------------

/// Complex `J0(z)` by its power series. The series loses accuracy to
/// cancellation for `|z|` beyond a few tens; callers use [`bessel_ratio`],
/// which switches to the asymptotic ratio there.
fn bessel_j0(z: C) -> C {
    let mut sum = C::real(1.0);
    let mut term = C::real(1.0);
    let z2 = z * z;
    let mut k = 1.0f64;
    loop {
        term = term * (-z2) / (4.0 * k * k);
        sum = sum + term;
        if term.abs() <= 1e-17 * sum.abs() && k > 2.0 {
            break;
        }
        k += 1.0;
        if k > 500.0 {
            break;
        }
    }
    sum
}

/// Complex `J1(z)` by its power series.
fn bessel_j1(z: C) -> C {
    let mut sum = C::real(1.0);
    let mut term = C::real(1.0);
    let z2 = z * z;
    let mut k = 1.0f64;
    loop {
        term = term * (-z2) / (4.0 * k * (k + 1.0));
        sum = sum + term;
        if term.abs() <= 1e-17 * sum.abs() && k > 2.0 {
            break;
        }
        k += 1.0;
        if k > 500.0 {
            break;
        }
    }
    (z / 2.0) * sum
}

/// `2 J1(z) / (z J0(z))`, the boundary-layer factor that appears in the
/// Stinson density and compressibility. For large `|z|` along the passive
/// direction (`Im z < 0`) the Bessel ratio tends to `i` and the factor to
/// `-2 i / z`; the asymptotic value is used beyond `|z| = 60` (error < 1.5%,
/// far below the LRF model error there).
fn bessel_ratio(z: C) -> C {
    if z.abs() > 60.0 {
        C::real(-2.0) * C::i() / z
    } else {
        C::real(2.0) * bessel_j1(z) / (z * bessel_j0(z))
    }
}

impl Thermoviscous {
    /// Shear wave parameter `s = a sqrt(omega rho / mu)` (dimensionless).
    /// The LRF narrow-tube model requires `s` not large *and* the pressure to
    /// be uniform across the tube, `k0 a << 1`.
    pub fn shear_wave_number(&self, fluid: &Fluid, radius: f64, omega: f64) -> f64 {
        radius * (omega * fluid.rho / self.mu).sqrt()
    }

    /// Complex density `rho_c` of a circular tube of radius `a`.
    pub fn complex_density(&self, fluid: &Fluid, radius: f64, omega: f64) -> C {
        let alpha = (C::new(0.0, -1.0) * omega * fluid.rho / self.mu).sqrt();
        let z = alpha * radius;
        C::real(fluid.rho) / (C::real(1.0) - bessel_ratio(z))
    }

    /// Complex compressibility `C_c` of a circular tube of radius `a`.
    pub fn complex_compressibility(&self, fluid: &Fluid, radius: f64, omega: f64) -> C {
        let alpha = (C::new(0.0, -1.0) * omega * fluid.rho / self.mu).sqrt();
        let zt = alpha * (self.pr.sqrt() * radius);
        let b0 = fluid.bulk();
        let thermal = (C::real(1.0) + bessel_ratio(zt) * (self.gamma - 1.0)) / b0;
        thermal
    }

    /// Characteristic specific impedance `Z_c = sqrt(rho_c / C_c)` [Pa s/m].
    pub fn characteristic_impedance(&self, fluid: &Fluid, radius: f64, omega: f64) -> C {
        let rho_c = self.complex_density(fluid, radius, omega);
        let c_c = self.complex_compressibility(fluid, radius, omega);
        (rho_c / c_c).sqrt()
    }

    /// Complex wavenumber `k_c = omega sqrt(rho_c C_c)` [1/m].
    pub fn wavenumber(&self, fluid: &Fluid, radius: f64, omega: f64) -> C {
        let rho_c = self.complex_density(fluid, radius, omega);
        let c_c = self.complex_compressibility(fluid, radius, omega);
        (rho_c * c_c).sqrt() * omega
    }
}

// ---------------------------------------------------------------------------
// Neck end corrections
// ---------------------------------------------------------------------------

/// The open end of a neck: baffled (flanged) or free (unflanged).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EndCondition {
    /// Opening in an infinite baffle. `delta/a = 8/(3 pi) = 0.8488`.
    Flanged,
    /// Free-space opening. `delta/a = 0.6133` (Rayleigh; Levine & Schwinger).
    Unflanged,
}

impl EndCondition {
    /// Low-frequency end-correction coefficient `delta / a`.
    pub fn coefficient(self) -> f64 {
        match self {
            EndCondition::Flanged => 8.0 / (3.0 * PI),
            EndCondition::Unflanged => 0.6133,
        }
    }
}

// ---------------------------------------------------------------------------
// The Helmholtz resonator
// ---------------------------------------------------------------------------

/// A Helmholtz resonator: a cavity of volume `V` behind a circular neck of
/// radius `a` and length `L`.
#[derive(Clone, Copy, Debug)]
pub struct Helmholtz {
    /// Neck radius `a` [m].
    pub neck_radius: f64,
    /// Geometric neck length `L` [m].
    pub neck_length: f64,
    /// Cavity volume `V` [m^3].
    pub cavity_volume: f64,
    /// End condition at the mouth (opening to the exterior).
    pub mouth: EndCondition,
    /// End condition at the neck–cavity junction.
    pub junction: EndCondition,
}

impl Helmholtz {
    pub fn new(
        neck_radius: f64,
        neck_length: f64,
        cavity_volume: f64,
        mouth: EndCondition,
        junction: EndCondition,
    ) -> Self {
        assert!(neck_radius > 0.0 && neck_length > 0.0 && cavity_volume > 0.0);
        Self {
            neck_radius,
            neck_length,
            cavity_volume,
            mouth,
            junction,
        }
    }

    /// Neck cross-sectional area `S = pi a^2`.
    pub fn neck_area(&self) -> f64 {
        PI * self.neck_radius * self.neck_radius
    }

    /// Total end correction `delta_out + delta_in`.
    pub fn end_correction(&self) -> f64 {
        (self.mouth.coefficient() + self.junction.coefficient()) * self.neck_radius
    }

    /// Effective neck length `L_eff = L + delta`.
    pub fn effective_length(&self) -> f64 {
        self.neck_length + self.end_correction()
    }

    /// Lumped resonance `omega0 = c sqrt(S/(V L_eff))`.
    pub fn omega0(&self, fluid: &Fluid) -> f64 {
        fluid.c * (self.neck_area() / (self.cavity_volume * self.effective_length())).sqrt()
    }

    /// Lumped resonance frequency `f0 = omega0/(2 pi)`.
    pub fn frequency0(&self, fluid: &Fluid) -> f64 {
        self.omega0(fluid) / (2.0 * PI)
    }

    /// Acoustic mass of the neck `M_a = rho L_eff / S`.
    pub fn neck_acoustic_mass(&self, fluid: &Fluid) -> f64 {
        fluid.rho * self.effective_length() / self.neck_area()
    }

    /// Acoustic compliance of the cavity `C_a = V / (rho c^2) = V / B`.
    pub fn cavity_acoustic_compliance(&self, fluid: &Fluid) -> f64 {
        self.cavity_volume / fluid.bulk()
    }

    /// Lossless lumped input impedance `i omega M_a + 1/(i omega C_a)`.
    pub fn lumped_impedance(&self, fluid: &Fluid, omega: f64) -> C {
        let m = self.neck_acoustic_mass(fluid);
        let c = self.cavity_acoustic_compliance(fluid);
        C::i() * omega * m + C::real(1.0) / (C::i() * omega * c)
    }

    /// Thermoviscous lumped input impedance. The neck uses the Stinson
    /// complex density at the neck radius (`Z_neck = i omega rho_c L_eff / S`),
    /// the cavity uses the Stinson complex compressibility
    /// (`Z_cav = 1/(i omega C_c V)`).
    pub fn lumped_impedance_lossy(
        &self,
        fluid: &Fluid,
        transport: &Thermoviscous,
        omega: f64,
    ) -> C {
        let rho_c = transport.complex_density(fluid, self.neck_radius, omega);
        let m = rho_c * self.effective_length() / self.neck_area();
        let c_c = transport.complex_compressibility(fluid, self.neck_radius, omega);
        let c = c_c * self.cavity_volume;
        C::i() * omega * m + C::real(1.0) / (C::i() * omega * c)
    }

    /// Low-frequency radiation impedance at the mouth, per unit volume
    /// velocity. Specific impedance `rho c[(ka)^2/2 + i coeff ka]`, divided by
    /// `S`; the imaginary part is exactly `omega M_rad` with `M_rad = rho
    /// delta/S`. Valid for `ka << 1`.
    pub fn radiation_impedance_lf(&self, fluid: &Fluid, omega: f64) -> C {
        let ka = omega * self.neck_radius / fluid.c;
        let coeff = self.mouth.coefficient();
        let z_specific = C::new(0.5 * ka * ka, coeff * ka) * fluid.impedance();
        z_specific / self.neck_area()
    }

    /// Length of the cavity if it is modelled as a duct of area `area`.
    pub fn cavity_length(&self, area: f64) -> f64 {
        self.cavity_volume / area
    }

    /// Exact input impedance at the neck mouth for a rigidly closed cavity.
    ///
    /// The neck (length `L_eff`, area `S`) and the cavity (length `V/A_c`,
    /// area `A_c`) are concatenated as lossless transmission-line ducts and
    /// the back face is rigid (`U = 0`), so `Z_in = -m22/m21`.
    pub fn input_impedance_exact(&self, fluid: &Fluid, omega: f64, cavity_area: f64) -> C {
        let m = self.duct_matrix_exact(fluid, omega, cavity_area);
        C::real(-1.0) * m.m22 / m.m21
    }

    /// Exact Stinson (LRF) input impedance at the neck mouth, rigid back face.
    /// Each segment uses its own complex `Z_c` and `k_c` at its own radius.
    pub fn input_impedance_exact_lossy(
        &self,
        fluid: &Fluid,
        transport: &Thermoviscous,
        omega: f64,
        cavity_radius: f64,
    ) -> C {
        let z_n = transport.characteristic_impedance(fluid, self.neck_radius, omega)
            / self.neck_area();
        let k_n = transport.wavenumber(fluid, self.neck_radius, omega);
        let m_n = Mat2::layer(z_n, k_n, self.effective_length());
        let a_c = PI * cavity_radius * cavity_radius;
        let z_c = transport.characteristic_impedance(fluid, cavity_radius, omega) / a_c;
        let k_c = transport.wavenumber(fluid, cavity_radius, omega);
        let l_c = self.cavity_volume / a_c;
        let m_c = Mat2::layer(z_c, k_c, l_c);
        let m = m_c * m_n;
        C::real(-1.0) * m.m22 / m.m21
    }

    /// Total characteristic matrix of the exact (lossless) neck + cavity ducts.
    fn duct_matrix_exact(&self, fluid: &Fluid, omega: f64, cavity_area: f64) -> Mat2 {
        let z_n = C::real(fluid.impedance() / self.neck_area());
        let k = C::real(omega / fluid.c);
        let m_n = Mat2::layer(z_n, k, self.effective_length());
        let z_c = C::real(fluid.impedance() / cavity_area);
        let l_c = self.cavity_length(cavity_area);
        let m_c = Mat2::layer(z_c, k, l_c);
        let m = m_c * m_n;
        m
    }

    /// Lowest series resonance (`Im Z` cross from negative to positive) of an
    /// exact input impedance, searched in `[0.05, 20] * omega0`.
    pub fn resonance_exact<F>(&self, fluid: &Fluid, impedance: F) -> f64
    where
        F: Fn(f64) -> C,
    {
        let w0 = self.omega0(fluid);
        let (lo, hi) = (0.05 * w0, 20.0 * w0);
        let n = 200_000;
        let mut prev_w = lo;
        let mut prev = impedance(lo).im;
        for i in 1..=n {
            let w = lo + (hi - lo) * i as f64 / n as f64;
            let cur = impedance(w).im;
            if prev < 0.0 && cur >= 0.0 {
                // Linear interpolation of the reactance zero.
                let frac = -prev / (cur - prev);
                return prev_w + frac * (w - prev_w);
            }
            prev_w = w;
            prev = cur;
        }
        f64::NAN
    }
}

// ---------------------------------------------------------------------------
// Elastic-walled cavity compliance
// ---------------------------------------------------------------------------

/// Hydrostatic compliance of an elastic-walled cavity, per unit cavity volume,
/// added to the fluid compliance `1/B0`.
///
/// For a spherical cavity of radius `R` in an infinite isotropic solid of
/// shear modulus `G`, Lamé gives a radial displacement `u(R) = p R/(4G)` under
/// internal pressure `p`, so `Delta V/V = 3p/(4G)` and the wall contribution
/// to `1/B_cav` is `3/(4G)`. For a long cylindrical cavity in plane strain the
/// outer field gives `Delta V/V = p/G`, contribution `1/G`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CavityShape {
    Spherical,
    Cylindrical,
}

impl CavityShape {
    /// Wall compliance per unit volume, `1/B_wall`.
    pub fn wall_compliance(self, solid: &Solid) -> f64 {
        let g = solid.mu();
        match self {
            CavityShape::Spherical => 3.0 / (4.0 * g),
            CavityShape::Cylindrical => 1.0 / g,
        }
    }
}

/// Effective cavity bulk modulus with the elastic-wall compliance folded in,
/// `1/B_cav = 1/B0 + 1/B_wall`.
pub fn cavity_bulk_modulus_elastic(fluid: &Fluid, solid: &Solid, shape: CavityShape) -> f64 {
    1.0 / (1.0 / fluid.bulk() + shape.wall_compliance(solid))
}

/// Ratio of the wall compliance to the fluid compliance. Wall compliance
/// matters when this is not small compared with one.
pub fn wall_to_fluid_compliance(fluid: &Fluid, solid: &Solid, shape: CavityShape) -> f64 {
    shape.wall_compliance(solid) * fluid.bulk()
}

// ---------------------------------------------------------------------------
// Side branches and periodic arrays
// ---------------------------------------------------------------------------

/// Pressure transmission amplitude of a junction with a side branch of
/// impedance `Z_b` on a waveguide of impedance `Z0`: `t = 2 Z_b/(2 Z_b + Z0)`.
pub fn side_branch_transmission(z_b: C, z0: f64) -> C {
    C::real(2.0) * z_b / (C::real(2.0) * z_b + z0)
}

/// Transmission loss of one side branch, `TL = -20 log10 |t|` [dB].
pub fn side_branch_loss_db(z_b: C, z0: f64) -> f64 {
    -20.0 * side_branch_transmission(z_b, z0).abs().log10()
}

/// A waveguide of area `A` loaded by a periodic array of Helmholtz side
/// branches with spacing `a`.
#[derive(Clone, Copy, Debug)]
pub struct SideBranchArray {
    /// Waveguide cross-sectional area [m^2].
    pub area: f64,
    /// Spacing between resonators [m].
    pub spacing: f64,
    pub resonator: Helmholtz,
    /// Thermoviscous transport used by the lossy branch impedance.
    pub transport: Thermoviscous,
    /// Phenomenological loss rate `Gamma` in the effective-medium formula
    /// [rad/s]. Zero gives the lossless folded model.
    pub gamma_loss: f64,
}

impl SideBranchArray {
    pub fn new(area: f64, spacing: f64, resonator: Helmholtz, transport: Thermoviscous) -> Self {
        assert!(area > 0.0 && spacing > 0.0);
        Self {
            area,
            spacing,
            resonator,
            transport,
            gamma_loss: 0.0,
        }
    }

    /// Waveguide characteristic impedance `Z0 = rho c / A`.
    pub fn waveguide_impedance(&self, fluid: &Fluid) -> f64 {
        fluid.impedance() / self.area
    }

    /// Lumped side-branch input impedance (lossless).
    pub fn branch_impedance(&self, fluid: &Fluid, omega: f64) -> C {
        self.resonator.lumped_impedance(fluid, omega)
    }

    /// Lumped side-branch input impedance with thermoviscous losses.
    pub fn branch_impedance_lossy(&self, fluid: &Fluid, omega: f64) -> C {
        self.resonator
            .lumped_impedance_lossy(fluid, &self.transport, omega)
    }

    /// Transmission coefficient of one side branch.
    pub fn transmission_one(&self, fluid: &Fluid, omega: f64) -> C {
        side_branch_transmission(
            self.branch_impedance(fluid, omega),
            self.waveguide_impedance(fluid),
        )
    }

    /// Exact unit-cell transfer matrix: a duct segment of length `a` followed
    /// by the shunt admittance `Y_b = 1/Z_b`.
    pub fn cell_matrix(&self, fluid: &Fluid, omega: f64) -> Mat2 {
        let z0 = self.waveguide_impedance(fluid);
        let k = omega / fluid.c;
        let duct = Mat2::layer(C::real(z0), C::real(k), self.spacing);
        let y_b = C::real(1.0) / self.branch_impedance(fluid, omega);
        // Flow entering the dead-end branch is removed from the through-flow
        // (`U2 = U1 - Y_b p`), so the shunt matrix has `m21 = -Y_b`.
        let shunt = Mat2 {
            m11: C::real(1.0),
            m12: C::real(0.0),
            m21: C::real(-1.0) * y_b,
            m22: C::real(1.0),
        };
        duct * shunt
    }

    /// `cos(q a) = tr(M_cell)/2` of the exact periodic array.
    pub fn bloch_cos(&self, fluid: &Fluid, omega: f64) -> C {
        let m = self.cell_matrix(fluid, omega);
        (m.m11 + m.m22) * 0.5
    }

    /// Chain matrix of a finite array of `cells` resonators, with a duct
    /// segment of length `a` on both sides of each shunt.
    pub fn stack_matrix(&self, fluid: &Fluid, omega: f64, cells: usize) -> Mat2 {
        let mut m = Mat2::identity();
        for _ in 0..cells {
            m = self.cell_matrix(fluid, omega) * m;
        }
        m
    }

    /// Pressure transmission amplitude of a finite array between two
    /// waveguides of impedance `Z0`, `t = 2/(A + B/Z0 + C Z0 + D)`.
    pub fn stack_transmission(&self, fluid: &Fluid, omega: f64, cells: usize) -> C {
        let m = self.stack_matrix(fluid, omega, cells);
        let z0 = C::real(self.waveguide_impedance(fluid));
        C::real(2.0) / (m.m11 + m.m12 / z0 + m.m21 * z0 + m.m22)
    }

    /// Transmission loss of a finite array, `-20 log10 |t|` [dB].
    pub fn stack_loss_db(&self, fluid: &Fluid, omega: f64, cells: usize) -> f64 {
        -20.0 * self.stack_transmission(fluid, omega, cells).abs().log10()
    }

    /// Exact Bloch wavevector of the periodic array, reduced to the first
    /// Brillouin zone.
    pub fn bloch_wavevector(&self, fluid: &Fluid, omega: f64) -> C {
        let h = self.bloch_cos(fluid, omega);
        let a = self.spacing;
        if h.im.abs() < 1e-12 && h.re.abs() <= 1.0 {
            C::real(h.re.clamp(-1.0, 1.0).acos() / a)
        } else {
            let alpha = h.abs().max(1.0).acosh();
            let re_q = if h.re < 0.0 { PI / a } else { 0.0 };
            C::new(re_q, alpha / a)
        }
    }

    /// Exact array in a band gap (`cos(q a)` outside `[-1, 1]`).
    pub fn in_band_gap(&self, fluid: &Fluid, omega: f64) -> bool {
        let h = self.bloch_cos(fluid, omega);
        h.im.abs() < 1e-9 && h.re.abs() > 1.0
    }

    /// Filling factor `F = V/(A a)`, the fraction of the cell volume occupied
    /// by the cavity fluid.
    pub fn filling_factor(&self) -> f64 {
        self.resonator.cavity_volume / (self.area * self.spacing)
    }

    /// Fang effective bulk modulus of the array,
    /// `1/B_eff = (1/B0)[1 - F omega0^2/(omega^2 - omega0^2 + i Gamma omega)]`.
    pub fn effective_bulk_modulus(&self, fluid: &Fluid, omega: f64) -> C {
        let w0 = self.resonator.omega0(fluid);
        let denom = C::new(omega * omega - w0 * w0, self.gamma_loss * omega);
        let bracket = C::real(1.0)
            - C::real(self.filling_factor() * w0 * w0) / denom;
        C::real(fluid.bulk()) / bracket
    }

    /// Folded wavenumber `k_eff = omega sqrt(rho0/B_eff)` from the effective
    /// bulk modulus (density unchanged).
    pub fn effective_medium_wavenumber(&self, fluid: &Fluid, omega: f64) -> C {
        let b_eff = self.effective_bulk_modulus(fluid, omega);
        (C::real(fluid.rho) / b_eff).sqrt() * omega
    }

    /// Sub-wavelength validity indicator `k a = omega a / c`.
    pub fn ka(&self, fluid: &Fluid, omega: f64) -> f64 {
        omega * self.spacing / fluid.c
    }
}

/// Fang-form effective bulk modulus from bare scalars,
/// `B_eff = B0 [1 - F omega0^2/(omega^2 - omega0^2 + i Gamma omega)]^{-1}`.
pub fn fang_bulk_modulus(b0: f64, omega: f64, omega0: f64, filling: f64, gamma_loss: f64) -> C {
    let denom = C::new(omega * omega - omega0 * omega0, gamma_loss * omega);
    C::real(b0) / (C::real(1.0) - C::real(filling * omega0 * omega0) / denom)
}

// ---------------------------------------------------------------------------
// Double-negative (negative B and negative rho)
// ---------------------------------------------------------------------------

/// A membrane-loaded waveguide gives an effective density
/// `rho_eff = rho0 / [1 - F_m omega^2/(omega^2 - omega_m^2 + i gamma_m omega)]`,
/// negative above the membrane resonance. Combined with the Helmholtz
/// negative-bulk-modulus band this is the double-negative medium of Lee et al.
/// (Phys. Rev. Lett. **104**, 054301 (2010)).
#[derive(Clone, Copy, Debug)]
pub struct MembraneDensity {
    pub omega_m: f64,
    pub filling: f64,
    pub gamma: f64,
}

impl MembraneDensity {
    pub fn new(omega_m: f64, filling: f64, gamma: f64) -> Self {
        assert!(omega_m > 0.0 && filling > 0.0);
        Self {
            omega_m,
            filling,
            gamma,
        }
    }
    /// Effective density `rho_eff(omega)`.
    pub fn effective_density(&self, rho0: f64, omega: f64) -> C {
        let denom = C::new(omega * omega - self.omega_m * self.omega_m, self.gamma * omega);
        C::real(rho0)
            / (C::real(1.0) - C::real(self.filling * omega * omega) / denom)
    }
}

/// Effective index `n = sqrt(rho_eff / B_eff)` of a combined resonator/membrane
/// medium. `Re n < 0` with both `Re rho_eff < 0` and `Re B_eff < 0` is the
/// double-negative band. The branch is chosen with `Re n >= 0` unless both
/// effective parameters are negative, where the negative-real root is taken.
pub fn double_negative_index(rho_eff: C, b_eff: C) -> C {
    let n2 = rho_eff / b_eff;
    let n = n2.sqrt();
    if rho_eff.re < 0.0 && b_eff.re < 0.0 {
        -n
    } else {
        n
    }
}
