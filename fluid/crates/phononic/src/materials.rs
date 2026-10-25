//! Sourced elastic properties of common fluids and solids.
//!
//! Every number carries its source. The values are representative handbook
//! room-temperature figures; they are elastic constants, not a materials
//! database, and the tests state the sensitivity of the benchmarks to them.

/// A homogeneous isotropic elastic solid, described by density and the
/// longitudinal (P) and transverse (S) phase speeds.
#[derive(Clone, Copy, Debug)]
pub struct Solid {
    pub name: &'static str,
    /// Mass density [kg m^-3].
    pub rho: f64,
    /// Compressional (P) wave speed [m s^-1].
    pub c_p: f64,
    /// Shear (S) wave speed [m s^-1].
    pub c_s: f64,
    pub source: &'static str,
}

impl Solid {
    /// First Lamé parameter `lambda = rho (c_p^2 - 2 c_s^2)`.
    pub fn lambda(&self) -> f64 {
        self.rho * (self.c_p * self.c_p - 2.0 * self.c_s * self.c_s)
    }
    /// Second Lamé parameter (shear modulus) `mu = rho c_s^2`.
    pub fn mu(&self) -> f64 {
        self.rho * self.c_s * self.c_s
    }
    /// Bulk modulus `K = lambda + 2 mu / 3`.
    pub fn bulk(&self) -> f64 {
        self.lambda() + 2.0 * self.mu() / 3.0
    }
    /// Young's modulus `E = mu (3 lambda + 2 mu) / (lambda + mu)`.
    pub fn young(&self) -> f64 {
        let (l, m) = (self.lambda(), self.mu());
        m * (3.0 * l + 2.0 * m) / (l + m)
    }
    /// Longitudinal acoustic impedance `Z = rho c_p`.
    pub fn impedance_p(&self) -> f64 {
        self.rho * self.c_p
    }
}

/// A homogeneous inviscid fluid, described by density and sound speed.
#[derive(Clone, Copy, Debug)]
pub struct Fluid {
    pub name: &'static str,
    /// Mass density [kg m^-3].
    pub rho: f64,
    /// Sound speed [m s^-1].
    pub c: f64,
    pub source: &'static str,
}

impl Fluid {
    /// Bulk modulus `B = rho c^2`.
    pub fn bulk(&self) -> f64 {
        self.rho * self.c * self.c
    }
    /// Acoustic impedance `Z = rho c`.
    pub fn impedance(&self) -> f64 {
        self.rho * self.c
    }
}

/// Water at 20 C. Density from the CRC Handbook; sound speed from the
/// IAPWS-97 / CRC tabulated value at 20 C, 1 atm.
pub const WATER: Fluid = Fluid {
    name: "water",
    rho: 998.2,
    c: 1482.0,
    source: "CRC Handbook of Chemistry and Physics, 20 C, 1 atm (rho=998.2, c=1482)",
};

/// Mild/structural steel. Density, P and S speeds are the standard NDT
/// handbook values.
pub const STEEL: Solid = Solid {
    name: "steel",
    rho: 7850.0,
    c_p: 5960.0,
    c_s: 3230.0,
    source: "Kinsler et al., Fundamentals of Acoustics / NDT handbook steel (5940-5960, 3200-3230)",
};

/// Aluminium (wrought, 2024/6061 range). Handbook NDT values.
pub const ALUMINIUM: Solid = Solid {
    name: "aluminium",
    rho: 2700.0,
    c_p: 6420.0,
    c_s: 3040.0,
    source: "NDT handbook aluminium (rho=2700, c_p=6420, c_s=3040)",
};

/// PMMA (poly(methyl methacrylate)), a common stiff polymer. Sound speeds
/// from ultrasonic handbook tables; the low shear speed is the reason
/// polymer plates show strong coincidence transmission.
pub const PMMA: Solid = Solid {
    name: "PMMA",
    rho: 1180.0,
    c_p: 2750.0,
    c_s: 1330.0,
    source: "ultrasonic handbook PMMA (rho=1180, c_p=2750, c_s=1330)",
};

/// A representative room-temperature engineering epoxy. Sound speeds are
/// estimates typical of the epoxy matrix used in phononic-crystal papers;
/// they are low, giving a large contrast against steel or tungsten.
pub const EPOXY: Solid = Solid {
    name: "epoxy",
    rho: 1180.0,
    c_p: 2530.0,
    c_s: 1160.0,
    source: "estimate; representative epoxy (rho=1180, c_L=2530, c_T=1160) as in the phononic-crystal literature",
};
