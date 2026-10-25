//! Effective media: folding a sub-wavelength structure into a homogeneous one.
//!
//! This module is the *coarse* rung above the exact transfer matrix. When the
//! period of a composite is much smaller than the wavelength, the wave cannot
//! resolve the individual inclusions and the material acts as a homogeneous
//! medium. Two closures are provided.
//!
//! * **Maxwell Garnett** — isolated inclusions of permittivity `eps_i` in a
//!   host `eps_m`, volume fraction `f`:
//!   `(eps - eps_m)/(eps + 2 eps_m) = f (eps_i - eps_m)/(eps_i + 2 eps_m)`.
//!   Valid at low fill fraction.
//! * **Bruggeman** — a symmetric, self-consistent mixture that treats both
//!   constituents on an equal footing (percolation-capable):
//!   `f (eps_i - eps)/(eps_i + 2 eps) + (1-f)(eps_m - eps)/(eps_m + 2 eps) = 0`.
//!
//! The exact long-wavelength limit of a *layered* medium is known in closed
//! form and is anisotropic (form birefringence):
//!
//! ```text
//!   eps_parallel      = f eps1 + (1-f) eps2          (E in the plane)
//!   1/eps_perp        = f/eps1 + (1-f)/eps2          (E normal to the planes)
//! ```
//!
//! `anisotropic_slab_reflectance` evaluates the p-polarised response of the
//! homogenised uniaxial slab, so a finite stack can be compared with the
//! folded model. The validity indicator is the period-to-wavelength ratio
//! `Lambda/lambda`; the missing correction is `O((Lambda/lambda)^2)`. The
//! `homogenisation_error` envelope grows with that ratio (and with index
//! contrast) and marks where the fine structure must be unfolded back to
//! `multilayer`.
use crate::complex::C;
use crate::multilayer::{characteristic_matrix, Polarization};

/// Maxwell Garnett effective permittivity.
pub fn maxwell_garnett(host: C, inclusion: C, f: f64) -> C {
    assert!((0.0..=1.0).contains(&f));
    let alpha = (inclusion - host) / (inclusion + host * 2.0);
    host * ((C::real(1.0) + alpha * (2.0 * f)) / (C::real(1.0) - alpha * f))
}

/// Bruggeman effective permittivity (symmetric mixture).
pub fn bruggeman(eps1: C, eps2: C, f: f64) -> C {
    assert!((0.0..=1.0).contains(&f));
    // 2 eps^2 - b eps - eps1 eps2 = 0, b = (3f-1) eps1 + (2-3f) eps2.
    let b = eps1 * (3.0 * f - 1.0) + eps2 * (2.0 - 3.0 * f);
    let disc = b * b + eps1 * eps2 * 8.0;
    let root = disc.sqrt();
    let r_plus = (b + root) * 0.25;
    let r_minus = (b - root) * 0.25;
    // The passive physical root has non-negative imaginary part.
    if r_plus.im >= r_minus.im {
        r_plus
    } else {
        r_minus
    }
}

/// Exact long-wavelength in-plane permittivity of a two-component laminate.
pub fn layered_parallel(eps1: C, eps2: C, f: f64) -> C {
    assert!((0.0..=1.0).contains(&f));
    eps1 * f + eps2 * (1.0 - f)
}

/// Exact long-wavelength normal permittivity of a two-component laminate.
pub fn layered_perpendicular(eps1: C, eps2: C, f: f64) -> C {
    assert!((0.0..=1.0).contains(&f));
    (C::real(f) / eps1 + C::real(1.0 - f) / eps2).recip()
}

/// Homogenised uniaxial permittivity of a laminate.
#[derive(Clone, Copy, Debug)]
pub struct Uniaxial {
    pub eps_parallel: C,
    pub eps_perpendicular: C,
}

impl Uniaxial {
    pub fn from_layers(eps1: C, eps2: C, f: f64) -> Self {
        Self {
            eps_parallel: layered_parallel(eps1, eps2, f),
            eps_perpendicular: layered_perpendicular(eps1, eps2, f),
        }
    }
    pub fn is_birefringent(&self) -> bool {
        (self.eps_parallel - self.eps_perpendicular).abs() > 1e-12
    }
    /// Intrinsic birefringence `n_e - n_o` (real parts), useful as a diagnostic.
    pub fn birefringence(&self) -> f64 {
        self.eps_perpendicular.sqrt().re - self.eps_parallel.sqrt().re
    }
}

fn isotropic_admittance(n: f64, s: f64, pol: Polarization) -> C {
    let q = (C::real(n * n - s * s)).sqrt();
    match pol {
        Polarization::TE => q,
        Polarization::TM => C::real(n * n) / q,
    }
}

/// Reflectance of a homogeneous uniaxial slab of thickness `d` between
/// isotropic media `n0` and `ns`, at angle `theta0`.
pub fn anisotropic_slab_reflectance(
    medium: Uniaxial,
    d: f64,
    n0: f64,
    ns: f64,
    lambda: f64,
    theta0: f64,
    pol: Polarization,
) -> f64 {
    let s = n0 * theta0.sin();
    let (q, eta) = match pol {
        Polarization::TE => {
            let q = (medium.eps_parallel - C::real(s * s)).sqrt();
            (q, q)
        }
        Polarization::TM => {
            let q = (medium.eps_parallel
                - medium.eps_parallel / medium.eps_perpendicular * C::real(s * s))
            .sqrt();
            (q, medium.eps_parallel / q)
        }
    };
    let eta0 = isotropic_admittance(n0, s, pol);
    let etas = isotropic_admittance(ns, s, pol);
    let m = characteristic_matrix(q, eta, d, lambda);
    let denom = m.m11 - m.m12 * etas;
    let y = (m.m22 * etas - m.m21) / denom;
    let r = (eta0 - y) / (eta0 + y);
    r.abs2()
}

/// Leading-order homogenisation reflectance-error indicator for a laminate.
///
/// The measured normal- and oblique-incidence reflectance error of a laminate
/// against its folded uniaxial slab follows `err ~ 2 (Lambda/lambda)^2
/// (1 + contrast)`, where `contrast = |eps1 - eps2| / (|eps1| + |eps2|)`.
/// This is an empirical envelope extracted from the convergence benchmark
/// (`tests/optics.rs`), not a fitted physical constant: it is the leading
/// quadratic term that tells the caller when the fine structure must be
/// unfolded. Unfold when it exceeds the target reflectance accuracy.
pub fn homogenisation_error(period: f64, lambda: f64, eps1: C, eps2: C) -> f64 {
    let contrast = (eps1 - eps2).abs() / (eps1.abs() + eps2.abs()).max(1e-30);
    2.0 * (period / lambda).powi(2) * (1.0 + contrast)
}

/// Simple ratio indicator: unfold the fine structure when this exceeds the
/// tolerance implied by `homogenisation_error`.
pub fn period_over_wavelength(period: f64, lambda: f64) -> f64 {
    period / lambda
}

// ---------------------------------------------------------------------------
// Two-dimensional mixtures (the long-wavelength limit of a 2D crystal)
// ---------------------------------------------------------------------------

/// Exact long-wavelength in-plane permittivity of a 2D array of circular
/// inclusions for the electric field parallel to the inclusion axes (TM): the
/// field is continuous along the axes, so the permittivity is the area average
/// `eps_eff = f eps_rod + (1-f) eps_bg`. This is exact at every fill fraction.
pub fn area_average_2d(eps_rod: C, eps_bg: C, f: f64) -> C {
    assert!((0.0..=1.0).contains(&f));
    eps_rod * f + eps_bg * (1.0 - f)
}

/// Two-dimensional Maxwell Garnett permittivity for circular cylinders with
/// the field perpendicular to the axes (TE). The depolarisation factor of a
/// cylinder is `1/2`, so with `beta = (eps_i - eps_m)/(eps_i + eps_m)`,
/// `eps_eff = eps_m (1 + f beta)/(1 - f beta)`. Valid at low fill fraction; it
/// approaches `eps_i` at `f = 1` and is the standard closure used as the TE
/// bound.
pub fn maxwell_garnett_2d(host: C, inclusion: C, f: f64) -> C {
    assert!((0.0..=1.0).contains(&f));
    let beta = (inclusion - host) / (inclusion + host);
    host * ((C::real(1.0) + beta * f) / (C::real(1.0) - beta * f))
}

/// Effective index `sqrt(eps_eff)` of a 2D mixture, with the principal branch.
pub fn effective_index_2d(eps_eff: C) -> C {
    eps_eff.sqrt()
}

/// Leading-order homogenisation error indicator for a 2D crystal. Same form as
/// the 1D laminate indicator, with the lattice constant in place of the
/// period: `2 (a/lambda)^2 (1 + contrast)`,
/// `contrast = |eps1 - eps2| / (|eps1| + |eps2|)`. Empirical envelope from the
/// PWE convergence benchmark, used as a switch, not a certificate.
pub fn homogenisation_error_2d(a: f64, lambda: f64, eps1: C, eps2: C) -> f64 {
    let contrast = (eps1 - eps2).abs() / (eps1.abs() + eps2.abs()).max(1e-30);
    2.0 * (a / lambda).powi(2) * (1.0 + contrast)
}
