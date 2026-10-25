//! The emergent boundaries of the light → heat → sound chain, as cheap
//! indicator functions.
//!
//! Each coarse model in this crate is valid in a regime, and each regime is
//! decided by a dimensionless number that costs almost nothing to evaluate:
//!
//! | coarse model | indicator | threshold | unfold |
//! |---|---|---|---|
//! | Beer–Lambert heating (no field solve) | `n d / L_c` (coherence) | 1 | [`photonics::multilayer`] transfer matrix |
//! | instantaneous thermal deposition | thermal confinement `tau_p / (delta^2/alpha)` | 1 | solve heat conduction during the pulse |
//! | instantaneous stress deposition | stress confinement `tau_p / (delta/c)` | 1 | solve the coupled acoustic problem during the pulse |
//! | conduction-only heat transport | Rayleigh number `Ra` | 1708 | unfold to the flow solver |
//!
//! The thresholds are exact linear-stability or asymptotic results, not tuned
//! constants; the tests exercise each at clear-cut values on both sides.
//!
//! Note: the code uses `delta_t` for a temperature difference and `delta` for
//! the optical absorption depth; they are different quantities.

/// Absorption depth `delta = 1 / mu_a` [m] for an absorption coefficient
/// `mu_a` [m^-1].
pub fn absorption_depth(mu_a: f64) -> f64 {
    1.0 / mu_a
}

/// Thermal-confinement ratio `tau_p / (delta^2 / alpha)`.
///
/// Below one, the pulse is shorter than the thermal relaxation time of the
/// absorption depth, so the deposited heat can be treated as instantaneous.
pub fn thermal_confinement_ratio(tau_pulse: f64, alpha: f64, delta: f64) -> f64 {
    tau_pulse * alpha / (delta * delta)
}

/// Whether the pulse is thermally confined (`ratio < 1`).
pub fn is_thermally_confined(tau_pulse: f64, alpha: f64, delta: f64) -> bool {
    thermal_confinement_ratio(tau_pulse, alpha, delta) < 1.0
}

/// Stress-confinement ratio `tau_p / (delta / c)`.
///
/// Below one, the pulse is shorter than the acoustic transit time across the
/// absorption depth, so the initial pressure is the local deposited energy
/// density (the usual `p0 = Gamma mu_a F` photoacoustic initial condition).
pub fn stress_confinement_ratio(tau_pulse: f64, c: f64, delta: f64) -> f64 {
    tau_pulse * c / delta
}

/// Whether the pulse is stress confined (`ratio < 1`).
pub fn is_stress_confined(tau_pulse: f64, c: f64, delta: f64) -> bool {
    stress_confinement_ratio(tau_pulse, c, delta) < 1.0
}

/// Which optical model the coherence indicator selects.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OpticalRegime {
    /// Fringes are resolved: use the coherent transfer matrix.
    Coherent,
    /// Phase is lost: use the incoherent intensity sum / Beer–Lambert.
    Incoherent,
}

/// Coherence indicator `n d / L_c` (optical thickness over coherence length).
pub fn coherence_indicator(n: f64, d: f64, coherence_length: f64) -> f64 {
    n * d / coherence_length
}

/// Coherent below `n d / L_c = 1`, incoherent above it.
pub fn optical_regime(n: f64, d: f64, coherence_length: f64) -> OpticalRegime {
    if coherence_indicator(n, d, coherence_length) < 1.0 {
        OpticalRegime::Coherent
    } else {
        OpticalRegime::Incoherent
    }
}

/// Critical Rayleigh number for the onset of steady Rayleigh–Bénard convection
/// in a horizontal layer between rigid (no-slip), isothermal boundaries.
/// Drazin & Reid, *Hydrodynamic Stability* (1981), and Chandrasekhar,
/// *Hydrodynamic and Hydromagnetic Stability* (1961): `Ra_c = 1707.762`.
pub const RA_CRIT_RAYLEIGH_BENARD: f64 = 1707.762;

/// Classical Rayleigh number `Ra = g beta dT L^3 / (nu alpha)` for a layer of
/// depth `L` held at a temperature difference `dT`.
pub fn rayleigh_benard(g: f64, beta: f64, delta_t: f64, l: f64, nu: f64, alpha: f64) -> f64 {
    g * beta * delta_t * l * l * l / (nu * alpha)
}

/// Rayleigh number from a heat flux `flux` [W m^-2] entering the base of a
/// layer, using the conduction temperature scale `dT = flux * L / k`.
pub fn rayleigh_from_flux(
    g: f64,
    beta: f64,
    flux: f64,
    l: f64,
    k: f64,
    nu: f64,
    alpha: f64,
) -> f64 {
    rayleigh_benard(g, beta, flux * l / k, l, nu, alpha)
}

/// Critical internal Rayleigh number for uniform volumetric heating of a
/// horizontal layer between rigid, isothermal boundaries. Sparrow, Goldstein
/// & Jonsson, "Thermal instability in a horizontal fluid layer: effect of
/// boundary conditions and non-linear temperature profile", J. Fluid Mech. 18,
/// 513 (1964), give `Ra_I,c = 2772` for the definition
/// `Ra_I = g beta q L^5 / (k nu alpha)`.
pub const RA_CRIT_INTERNAL: f64 = 2772.0;

/// Internal-heating Rayleigh number `Ra_I = g beta q L^5 / (k nu alpha)` for a
/// volumetric heat source `q` [W m^-3] distributed over a layer of depth `L`.
///
/// The fifth power of `L` follows from the conduction temperature scale
/// `q L^2 / k` raised to `L^3`.
pub fn rayleigh_internal(g: f64, beta: f64, q: f64, l: f64, k: f64, nu: f64, alpha: f64) -> f64 {
    g * beta * q * l.powi(5) / (k * nu * alpha)
}

/// Whether the layer is expected to convect, given the classical
/// Rayleigh–Bénard number.
pub fn convects_rayleigh_benard(ra: f64) -> bool {
    ra > RA_CRIT_RAYLEIGH_BENARD
}

/// Whether a volumetrically heated layer is expected to convect given the
/// internal Rayleigh number.
pub fn convects_internal(ra_i: f64) -> bool {
    ra_i > RA_CRIT_INTERNAL
}

/// Choose the convection indicator: surface-like heating when the absorption
/// depth is a small fraction of the layer (`mu_a L >~ 3`), volumetric
/// otherwise. Returns the Rayleigh number that should be compared with the
/// corresponding critical value.
pub fn absorbed_heating_rayleigh(
    g: f64,
    beta: f64,
    absorbed_flux: f64,
    q_vol: f64,
    l: f64,
    k: f64,
    nu: f64,
    alpha: f64,
    mu_a: f64,
) -> (f64, f64) {
    if mu_a * l >= 3.0 {
        (
            rayleigh_from_flux(g, beta, absorbed_flux, l, k, nu, alpha),
            RA_CRIT_RAYLEIGH_BENARD,
        )
    } else {
        (
            rayleigh_internal(g, beta, q_vol, l, k, nu, alpha),
            RA_CRIT_INTERNAL,
        )
    }
}
