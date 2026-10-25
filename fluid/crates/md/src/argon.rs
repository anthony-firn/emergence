//! Reduced-to-SI closures for argon and comparison with measured liquid
//! argon.
//!
//! The LJ parameters are the standard argon values used by the materials
//! crate: `sigma = 3.405e-10 m`, `epsilon/k_B = 119.8 K` (Rowley &
//! Nicholson 1972, via standard LJ tables), `m = 39.948 u`. With these,
//! the reduced units map as
//!
//! | quantity | unit |
//! |---|---|
//! | length | `sigma = 3.405e-10 m` |
//! | time | `sigma sqrt(m/epsilon) = 2.156e-12 s` |
//! | viscosity | `sqrt(m epsilon)/sigma^2 = 9.0346e-5 Pa s` |
//! | thermal conductivity | `k_B sqrt(epsilon/m)/sigma^2 = 1.8804e-2 W/(m K)` |
//! | self-diffusion | `sigma sqrt(epsilon/m) = 5.3767e-8 m^2/s` |
//! | pressure | `epsilon/sigma^3 = 4.1898e7 Pa` |
//!
//! Measured liquid-argon values near 90 K are taken from the transport
//! correlations of B. A. Younglove and H. J. M. Hanley, "The viscosity and
//! thermal conductivity coefficients of gaseous and liquid argon",
//! *J. Phys. Chem. Ref. Data* **15**, 1323 (1986), as evaluated by the
//! CoolProp library together with the argon reference equation of state of
//! Tegeler, Span & Wagner (1999). On the saturated liquid line:
//! `T = 90.0 K`, `rho = 1378.6 kg/m^3`, `mu = 2.4004e-4 Pa s`,
//! `k = 0.12451 W/(m K)`.

/// Boltzmann constant, exact SI value.
pub const KB: f64 = 1.380_649e-23;
/// Atomic mass unit, CODATA 2018.
pub const ATOMIC_MASS: f64 = 1.660_539_066_60e-27;

/// Argon Lennard-Jones length, metres.
pub const SIGMA_M: f64 = 3.405e-10;
/// Argon Lennard-Jones well depth in kelvin.
pub const EPSILON_OVER_KB_K: f64 = 119.8;
/// Argon molar mass in atomic mass units.
pub const MASS_U: f64 = 39.948;

/// Argon molecular mass, kg.
pub fn mass_kg() -> f64 {
    MASS_U * ATOMIC_MASS
}
/// Argon Lennard-Jones well depth, joules.
pub fn epsilon_j() -> f64 {
    EPSILON_OVER_KB_K * KB
}
/// Reduced time unit `sigma sqrt(m/epsilon)`, seconds.
pub fn time_unit() -> f64 {
    SIGMA_M * (mass_kg() / epsilon_j()).sqrt()
}
/// Reduced velocity unit `sqrt(epsilon/m)`, m/s.
pub fn velocity_unit() -> f64 {
    SIGMA_M / time_unit()
}
/// Reduced pressure unit `epsilon/sigma^3`, Pa.
pub fn pressure_unit() -> f64 {
    epsilon_j() / SIGMA_M.powi(3)
}
/// Reduced viscosity unit `sqrt(m epsilon)/sigma^2`, Pa s.
pub fn viscosity_unit() -> f64 {
    (mass_kg() * epsilon_j()).sqrt() / SIGMA_M.powi(2)
}
/// Reduced thermal-conductivity unit `k_B sqrt(epsilon/m)/sigma^2`,
/// W/(m K).
pub fn thermal_conductivity_unit() -> f64 {
    KB * (epsilon_j() / mass_kg()).sqrt() / SIGMA_M.powi(2)
}
/// Reduced self-diffusion unit `sigma sqrt(epsilon/m)`, m^2/s.
pub fn diffusion_unit() -> f64 {
    SIGMA_M * (epsilon_j() / mass_kg()).sqrt()
}

/// A set of reduced transport coefficients and their SI images.
#[derive(Clone, Copy, Debug)]
pub struct ArgonClosure {
    pub temperature_k: f64,
    pub density_kg_m3: f64,
    /// Reduced viscosity `eta*`.
    pub eta_star: f64,
    /// Reduced thermal conductivity `kappa*`.
    pub kappa_star: f64,
    /// Reduced self-diffusion `D*`.
    pub d_star: f64,
    /// SI viscosity, Pa s.
    pub viscosity_pa_s: f64,
    /// SI thermal conductivity, W/(m K).
    pub conductivity_w_m_k: f64,
    /// SI self-diffusion, m^2/s.
    pub diffusion_m2_s: f64,
}

/// Fold reduced LJ transport coefficients at `(t_star, rho_star)` into SI
/// argon closures.
pub fn argon_closure(
    t_star: f64,
    rho_star: f64,
    eta_star: f64,
    kappa_star: f64,
    d_star: f64,
) -> ArgonClosure {
    ArgonClosure {
        temperature_k: t_star * EPSILON_OVER_KB_K,
        density_kg_m3: rho_star * mass_kg() / SIGMA_M.powi(3),
        eta_star,
        kappa_star,
        d_star,
        viscosity_pa_s: eta_star * viscosity_unit(),
        conductivity_w_m_k: kappa_star * thermal_conductivity_unit(),
        diffusion_m2_s: d_star * diffusion_unit(),
    }
}

/// Measured saturated-liquid-argon transport at a fixed temperature, from
/// Younglove & Hanley (1986) via CoolProp / Tegeler et al. (1999).
#[derive(Clone, Copy, Debug)]
pub struct LiquidArgonReference {
    pub temperature_k: f64,
    pub density_kg_m3: f64,
    pub viscosity_pa_s: f64,
    pub conductivity_w_m_k: f64,
}

/// Saturated liquid argon at 90.0 K (Younglove & Hanley 1986, via CoolProp).
pub fn liquid_argon_90k() -> LiquidArgonReference {
    LiquidArgonReference {
        temperature_k: 90.0,
        density_kg_m3: 1378.6,
        viscosity_pa_s: 2.4004e-4,
        conductivity_w_m_k: 0.124_51,
    }
}

/// Saturated liquid argon at the temperature corresponding to the LJ
/// triple point `T* = 0.722` (86.50 K), used to compare with the LJ
/// reference state of Meier et al. (2004).
pub fn liquid_argon_at_lj_triple() -> LiquidArgonReference {
    LiquidArgonReference {
        temperature_k: 0.722 * EPSILON_OVER_KB_K,
        density_kg_m3: 1400.3,
        viscosity_pa_s: 2.6676e-4,
        conductivity_w_m_k: 0.129_64,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn units_have_the_right_magnitude() {
        assert!((time_unit() - 2.156_349e-12).abs() / 2.156_349e-12 < 1e-5);
        assert!((viscosity_unit() - 9.034_578e-5).abs() / 9.034_578e-5 < 1e-5);
        assert!((thermal_conductivity_unit() - 1.880_386e-2).abs() / 1.880_386e-2 < 1e-5);
        assert!((diffusion_unit() - 5.376_691e-8).abs() / 5.376_691e-8 < 1e-5);
        assert!((pressure_unit() - 4.189_756e7).abs() / 4.189_756e7 < 1e-5);
    }

    #[test]
    fn closure_round_trips() {
        let c = argon_closure(0.722, 0.8442, 3.258, 6.8, 0.032);
        assert!((c.viscosity_pa_s - 2.943_5e-4).abs() / 2.943_5e-4 < 1e-3);
        assert!((c.temperature_k - 86.4956).abs() < 1e-3);
    }
}
