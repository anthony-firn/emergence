use crate::{A_RAD, C};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Regime {
    LocalEquilibrium,
    Diffusion,
    AngularPackets,
    WaveOpticsRequired,
}
/// This selector is an eligibility check, not a certificate. Beams and
/// boundaries can be non-isotropic even with a large cell optical depth.
pub struct Indicators {
    pub transport_depth: f64,
    pub mean_free_path_over_gradient_length: f64,
    pub reduced_flux: f64,
    pub absorption_relaxations: f64,
    pub lte_relative_error: f64,
    pub directional_structure: bool,
    pub wavelength_over_feature: f64,
}
impl Indicators {
    pub fn regime(&self) -> Regime {
        assert!([
            self.transport_depth,
            self.mean_free_path_over_gradient_length,
            self.reduced_flux,
            self.absorption_relaxations,
            self.lte_relative_error,
            self.wavelength_over_feature
        ]
        .iter()
        .all(|v| v.is_finite() && *v >= 0.0));
        if self.wavelength_over_feature > 0.01 {
            return Regime::WaveOpticsRequired;
        }
        if self.directional_structure
            || self.transport_depth < 10.0
            || self.reduced_flux > 0.1
            || self.mean_free_path_over_gradient_length > 0.1
        {
            return Regime::AngularPackets;
        }
        if self.absorption_relaxations >= 100.0 && self.lte_relative_error < 0.01 {
            Regime::LocalEquilibrium
        } else {
            Regime::Diffusion
        }
    }
}
#[derive(Clone, Copy, Debug)]
pub struct Exchange {
    pub temperature: f64,
    pub radiation: f64,
    pub matter_delta: f64,
}
/// Conservative backward-Euler LTE grey emission/absorption at fixed volume.
/// Cv is volumetric [J m^-3 K^-1], E [J m^-3], opacity [m^-1].
/// Frozen Cv and Planck-mean opacity: suitable for a split step away from
/// latent-heat plateaus; supply the actual enthalpy/EOS solve at a transition.
pub fn exchange(t: f64, cv: f64, e: f64, opacity: f64, dt: f64) -> Exchange {
    assert!(t >= 0.0 && t.is_finite() && cv > 0.0 && cv.is_finite() && e >= 0.0 && e.is_finite());
    assert!(opacity >= 0.0 && opacity.is_finite() && dt >= 0.0 && dt.is_finite());
    let s = C * opacity * dt;
    if s == 0.0 {
        return Exchange {
            temperature: t,
            radiation: e,
            matter_delta: 0.0,
        };
    }
    let w = if s.is_infinite() { 1.0 } else { s / (1.0 + s) };
    let trad = (e / A_RAD).sqrt().sqrt();
    let (mut lo, mut hi) = (t.min(trad), t.max(trad));
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        let residual = cv * (mid - t) + w * (A_RAD * mid.powi(4) - e);
        if residual > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let tn = 0.5 * (lo + hi);
    let en = ((1.0 - w) * e + w * A_RAD * tn.powi(4)).min(e + cv * t);
    let delta = e - en;
    Exchange {
        temperature: t + delta / cv,
        radiation: en,
        matter_delta: delta,
    }
}
/// Closed, uniform 1D grid, implicit radiation diffusion, O(n) solve.
/// k_transport = k_absorption + (1-g) k_scattering; units m^-1.
/// Requires optically thick cells. Vacuum/thin cells must use transport.
pub fn diffuse_closed(
    e: &mut [f64],
    k_transport: &[f64],
    dx: f64,
    dt: f64,
) -> Result<(), &'static str> {
    if e.is_empty()
        || e.len() != k_transport.len()
        || !dx.is_finite()
        || dx <= 0.0
        || !dt.is_finite()
        || dt < 0.0
    {
        return Err("invalid diffusion grid/time step");
    }
    if e.iter().any(|v| !v.is_finite() || *v < 0.0)
        || k_transport.iter().any(|k| !k.is_finite() || *k * dx < 10.0)
    {
        return Err(
            "diffusion requires positive radiation and transport optical depth >= 10 in every cell",
        );
    }
    let n = e.len();
    let mut diagonal = vec![1.0; n];
    let mut upper = vec![0.0; n];
    let mut lower = vec![0.0; n];
    for i in 0..n - 1 {
        // Harmonic diffusivity, equivalent to arithmetic opacity.
        let r = dt * C / (1.5 * (k_transport[i] + k_transport[i + 1]) * dx * dx);
        if !r.is_finite() {
            return Err("diffusion coefficient overflow");
        }
        diagonal[i] += r;
        diagonal[i + 1] += r;
        upper[i] = -r;
        lower[i + 1] = -r;
    }
    let mut rhs = e.to_vec();
    for i in 1..n {
        let r = lower[i] / diagonal[i - 1];
        diagonal[i] -= r * upper[i - 1];
        rhs[i] -= r * rhs[i - 1];
    }
    rhs[n - 1] /= diagonal[n - 1];
    for i in (0..n - 1).rev() {
        rhs[i] = (rhs[i] - upper[i] * rhs[i + 1]) / diagonal[i];
    }
    if rhs.iter().any(|v| !v.is_finite() || *v < 0.0) {
        return Err("diffusion solve exceeds floating point conditioning; reduce step");
    }
    e.copy_from_slice(&rhs);
    Ok(())
}
