//! The transfer-matrix method for planar stacks.
//!
//! A stack of isotropic layers is the exact local model of a thin film, a
//! Bragg mirror, an antireflection coating or a metallised mirror, at any
//! incidence angle and for both polarisations. Each layer is characterised by
//! its characteristic (Abelès) matrix, which relates the tangential electric
//! and magnetic fields at its two faces:
//!
//! ```text
//!   [ E ]         [ cos(delta)   i sin(delta) / eta ] [ E ]
//!   [   ]      =  [                                  ] [   ]     (forward into +z)
//!   [ H ]         [ i eta sin(delta)   cos(delta)    ] [ H ]
//! ```
//!
//! with `delta = 2 pi d q / lambda` and the tilted optical admittance
//! `eta = q` for TE (s) and `eta = N^2 / q` for TM (p). Here
//! `q = sqrt(N^2 - s^2)` is the normal wavevector, `s = N0 sin(theta0)` is the
//! conserved tangential wavevector, and `N = n + i kappa` is the complex
//! index. The matrix product in layer order gives the input admittance `Y`,
//! then `r = (eta0 - Y)/(eta0 + Y)`, `t = (1+r)/(M11 - M12 eta_sub)`.
//!
//! The exact flat-interface results (Fresnel, Brewster, total internal
//! reflection, the quarter-wave coating and mirror, energy conservation) are
//! the benchmarks. The physical sub-grid structure this rung resolves is the
//! interference of multiply reflected partial waves: it is a **spatial**
//! emergent boundary. When the layers are much thicker than the source
//! coherence length the coherent picture must be unfolded *down* to the
//! incoherent intensity sum (`crate::incoherent`); when the period is much
//! smaller than the wavelength the stack must be folded *up* into an
//! effective medium (`crate::effective_medium`). The cheap indicators for the
//! two hand-offs are the coherence-length-to-optical-thickness ratio and the
//! period-to-wavelength ratio.
use crate::complex::C;
use std::f64::consts::PI;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Polarization {
    /// TE / s: electric field parallel to the interfaces.
    TE,
    /// TM / p: magnetic field parallel to the interfaces.
    TM,
}

#[derive(Clone, Copy, Debug)]
pub struct Layer {
    /// Complex refractive index n + i kappa; kappa >= 0 for a passive medium.
    pub index: C,
    /// Complex relative permeability `mu = mu' + i mu''`. Use `C::real(1.0)`
    /// for a non-magnetic layer. `index^2 = eps * mu` for the physical branch.
    pub mu: C,
    /// Thickness in metres.
    pub thickness: f64,
}

impl Layer {
    pub fn new(n: f64, kappa: f64, thickness: f64) -> Self {
        assert!(n > 0.0 && n.is_finite() && thickness > 0.0 && thickness.is_finite());
        assert!(kappa.is_finite());
        Self {
            index: C::new(n, kappa),
            mu: C::real(1.0),
            thickness,
        }
    }
    /// Magnetic layer with complex index `n + i kappa` and complex relative
    /// permeability `mu_re + i mu_im`.
    pub fn magnetic(n: f64, kappa: f64, mu_re: f64, mu_im: f64, thickness: f64) -> Self {
        assert!(n.is_finite() && kappa.is_finite() && thickness > 0.0 && thickness.is_finite());
        Self {
            index: C::new(n, kappa),
            mu: C::new(mu_re, mu_im),
            thickness,
        }
    }
    /// Layer specified by relative permittivity and permeability. The
    /// refractive index uses the product of the principal square roots,
    /// `n = sqrt(eps) sqrt(mu)`, which selects the negative-index branch when
    /// both parts are negative (`eps = mu = -1` gives `n = -1`, not `+1`).
    pub fn from_eps_mu(eps: C, mu: C, thickness: f64) -> Self {
        assert!(thickness > 0.0 && thickness.is_finite());
        Self {
            index: eps.sqrt() * mu.sqrt(),
            mu,
            thickness,
        }
    }
    pub fn lossless(n: f64, thickness: f64) -> Self {
        Self::new(n, 0.0, thickness)
    }
    /// Quarter-wave optical thickness at `lambda` (design wavelength).
    pub fn quarter_wave(n: f64, lambda: f64) -> Self {
        Self::lossless(n, lambda / (4.0 * n))
    }
    /// Relative permittivity `eps = index^2 / mu` for the physical branch.
    pub fn epsilon(&self) -> C {
        self.index * self.index / self.mu
    }
    /// Index seen by the wave in this layer at the conserved tangential `s`.
    pub fn wavevector(&self, s: f64) -> C {
        normal_wavevector(self.index, s)
    }
}

/// Normal wavevector `q = N cos(theta) = N sqrt(1 - (s/N)^2)`, the physical
/// branch for a medium of signed refractive index `N`. Writing it as `N`
/// times the principal cosine (rather than `sqrt(N^2 - s^2)`) is what keeps
/// the sign of a negative-index medium: for `N < 0` the factor flips the
/// branch, so an evanescent wave grows in `+z` as it must for Pendry's
/// perfect lens. For a positive-index medium it reduces to
/// `sqrt(N^2 - s^2)` with `Im(q) >= 0` (absorption and total reflection both
/// decay along `+z`). Here `N = sqrt(eps) sqrt(mu)` is the physical branch and
/// `N^2 = eps mu`.
pub fn normal_wavevector(n: C, s: f64) -> C {
    n * (C::real(1.0) - C::real(s * s) / (n * n)).sqrt()
}

/// Tilted optical admittance. For a medium with permittivity `eps` and
/// permeability `mu`, `eta = q / mu` for TE and `eps / q = N^2/(mu q)` for TM;
/// setting `mu = 1` recovers the non-magnetic convention `eta = q` / `N^2/q`.
fn tilted_admittance(n: C, mu: C, q: C, pol: Polarization) -> C {
    match pol {
        Polarization::TE => q / mu,
        Polarization::TM => n * n / (mu * q),
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Mat2 {
    pub m11: C,
    pub m12: C,
    pub m21: C,
    pub m22: C,
}

impl Mat2 {
    pub fn identity() -> Self {
        Self {
            m11: C::real(1.0),
            m12: C::real(0.0),
            m21: C::real(0.0),
            m22: C::real(1.0),
        }
    }
    pub fn layer(delta: C, eta: C) -> Self {
        let i = C::i();
        let (cosd, sind) = complex_cos_sin(delta);
        Self {
            m11: cosd,
            m12: (i * sind) / eta,
            m21: i * eta * sind,
            m22: cosd,
        }
    }
    pub fn apply(self, v: (C, C)) -> (C, C) {
        (
            self.m11 * v.0 + self.m12 * v.1,
            self.m21 * v.0 + self.m22 * v.1,
        )
    }
}

impl std::ops::Mul for Mat2 {
    type Output = Self;
    fn mul(self, b: Self) -> Self {
        Self {
            m11: self.m11 * b.m11 + self.m12 * b.m21,
            m12: self.m11 * b.m12 + self.m12 * b.m22,
            m21: self.m21 * b.m11 + self.m22 * b.m21,
            m22: self.m21 * b.m12 + self.m22 * b.m22,
        }
    }
}

/// `(cos delta, sin delta)` for complex `delta`, via Euler's formula.
fn complex_cos_sin(delta: C) -> (C, C) {
    let e = (C::i() * delta).exp();
    let em = (C::i() * (-delta)).exp();
    ((e + em) * 0.5, (e - em) / (C::i() * 2.0))
}

/// Characteristic matrix from a normal wavevector `q` and tilted admittance
/// `eta`; both polarisations share this once `(q, eta)` are chosen.
pub fn characteristic_matrix(q: C, eta: C, thickness: f64, lambda: f64) -> Mat2 {
    let delta = q * (2.0 * PI * thickness / lambda);
    let (cosd, sind) = complex_cos_sin(delta);
    Mat2 {
        m11: cosd,
        m12: (C::i() * sind) / eta,
        m21: C::i() * eta * sind,
        m22: cosd,
    }
}

/// Characteristic matrix of one layer at conserved tangential wavevector `s`.
pub fn layer_matrix(
    index: C,
    mu: C,
    thickness: f64,
    lambda: f64,
    s: f64,
    pol: Polarization,
) -> Mat2 {
    let q = normal_wavevector(index, s);
    let eta = tilted_admittance(index, mu, q, pol);
    characteristic_matrix(q, eta, thickness, lambda)
}

/// A complete plane stack: incident medium, finite layers, exit medium.
#[derive(Clone, Debug)]
pub struct Stack {
    pub incident: C,
    pub layers: Vec<Layer>,
    pub exit: C,
}

/// Internal field state of one layer, kept for the heat-source profile.
#[derive(Clone, Copy, Debug)]
pub struct LayerField {
    pub z_start: f64,
    pub thickness: f64,
    pub index: C,
    pub mu: C,
    pub q: C,
    pub eta: C,
    /// Forward and backward electric-field amplitudes at `z_start`.
    pub a: C,
    pub b: C,
}

#[derive(Clone, Debug)]
pub struct Solution {
    /// Amplitude reflection coefficient (H-based admittance convention).
    pub r: C,
    /// Amplitude transmission coefficient.
    pub t: C,
    pub reflectance: f64,
    pub transmittance: f64,
    /// Fraction of incident power absorbed in the finite layers.
    pub absorptance: f64,
    /// Absorbed fraction in each layer (sums with `transmittance` to one).
    pub layer_absorptance: Vec<f64>,
    pub fields: Vec<LayerField>,
    pub incident_admittance: C,
    pub lambda: f64,
    /// Normalised incident `cos(theta)` factor `q0 = n0 cos(theta0)`.
    pub q0: f64,
}

impl Stack {
    pub fn new(incident: f64, layers: Vec<Layer>, exit: f64) -> Self {
        assert!(incident > 0.0 && exit > 0.0);
        Self {
            incident: C::real(incident),
            layers,
            exit: C::real(exit),
        }
    }
    pub fn with_complex(incident: C, layers: Vec<Layer>, exit: C) -> Self {
        Self {
            incident,
            layers,
            exit,
        }
    }
    pub fn single(incident: f64, layer: Layer, exit: f64) -> Self {
        Self::new(incident, vec![layer], exit)
    }

    /// Solve the stack at wavelength `lambda` (m) and incidence angle `theta`.
    pub fn solve(&self, lambda: f64, theta: f64, pol: Polarization) -> Solution {
        assert!(lambda > 0.0 && lambda.is_finite());
        assert!(theta.is_finite());
        self.solve_s(lambda, self.incident.re * theta.sin(), pol)
    }

    /// Solve the stack at a fixed conserved tangential wavevector
    /// `s = k_x / k0`. Unlike [`Stack::solve`] this accepts `s > 1`, which is
    /// the evanescent regime used for the perfect-lens and total-internal-
    /// reflection checks.
    pub fn solve_s(&self, lambda: f64, s: f64, pol: Polarization) -> Solution {
        assert!(lambda > 0.0 && lambda.is_finite());
        assert!(s.is_finite());
        let q0 = normal_wavevector(self.incident, s);
        let eta0 = tilted_admittance(self.incident, C::real(1.0), q0, pol);
        let q_sub = normal_wavevector(self.exit, s);
        let eta_sub = tilted_admittance(self.exit, C::real(1.0), q_sub, pol);

        let mut m = Mat2::identity();
        let n = self.layers.len();
        let mut fields: Vec<LayerField> = Vec::with_capacity(n);
        let mut matrices: Vec<Mat2> = Vec::with_capacity(n);
        let mut z = 0.0;
        for l in &self.layers {
            let q = l.wavevector(s);
            let eta = tilted_admittance(l.index, l.mu, q, pol);
            let mj = layer_matrix(l.index, l.mu, l.thickness, lambda, s, pol);
            matrices.push(mj);
            m = m * mj;
            fields.push(LayerField {
                z_start: z,
                thickness: l.thickness,
                index: l.index,
                mu: l.mu,
                q,
                eta,
                a: C::real(0.0),
                b: C::real(0.0),
            });
            z += l.thickness;
        }

        let denom = m.m11 - m.m12 * eta_sub;
        // At an input-admittance pole `denom = 0` the slab supports a bound
        // mode; for a driven wave the reflected amplitude is total and the
        // transmission ratio is indeterminate, so take `r = -1`, `t = 0`
        // rather than dividing by zero.
        let scale = m.m11.abs() + (m.m12 * eta_sub).abs();
        let (r, t) = if denom.abs() <= 1e-13 * scale.max(1e-300) {
            (C::real(-1.0), C::real(0.0))
        } else {
            let y = (m.m22 * eta_sub - m.m21) / denom;
            let r = (eta0 - y) / (eta0 + y);
            let t = (r + C::real(1.0)) / denom;
            (r, t)
        };
        let reflectance = r.abs2();
        let transmittance = if eta0.re.abs() > 0.0 {
            eta_sub.re / eta0.re * t.abs2()
        } else {
            0.0
        };

        // Back-propagate the (E, H) state to get each layer's field amplitudes.
        let mut state = (t, -eta_sub * t);
        let mut states = vec![state; n + 1];
        for j in (0..n).rev() {
            state = matrices[j].apply(state);
            states[j] = state;
        }
        let mut layer_absorptance = Vec::with_capacity(n);
        for (j, f) in fields.iter_mut().enumerate() {
            let (e_in, h_in) = states[j];
            let (e_out, h_out) = states[j + 1];
            let p_in = -((e_in * h_in.conj()).re);
            let p_out = -((e_out * h_out.conj()).re);
            layer_absorptance.push((p_in - p_out) / eta0.re);
            f.a = (e_in - h_in / f.eta) * 0.5;
            f.b = (e_in + h_in / f.eta) * 0.5;
        }
        let absorptance = 1.0 - reflectance - transmittance;
        Solution {
            r,
            t,
            reflectance,
            transmittance,
            absorptance,
            layer_absorptance,
            fields,
            incident_admittance: eta0,
            lambda,
            q0: q0.re,
        }
    }

    pub fn reflectance(&self, lambda: f64, theta: f64, pol: Polarization) -> f64 {
        self.solve(lambda, theta, pol).reflectance
    }
}

impl Solution {
    /// Electric field amplitude (relative to unit incident amplitude) at the
    /// depth `z` measured from the front interface, in the same polarisation.
    pub fn electric_field(&self, z: f64) -> C {
        for f in &self.fields {
            if z <= f.z_start + f.thickness + 1e-15 {
                let local = (z - f.z_start).max(0.0);
                let beta = f.q * (2.0 * PI * local / self.lambda);
                return f.a * (C::i() * beta).exp() + f.b * (C::i() * (-beta)).exp();
            }
        }
        C::real(0.0)
    }

    /// Time-averaged absorbed power density [W m^-3] at depth `z` for an
    /// incident irradiance `incident_irradiance` [W m^-2] (the component
    /// normal to the surface, i.e. the beam irradiance on the stack).
    ///
    /// `p = (2 pi / lambda) (Im(eps) / q0) I0 |E|^2` with `|E|` normalised to
    /// a unit incident amplitude. Integrating `p` over a layer recovers
    /// `I0 * layer_absorptance`, which is the coupling sink handed to a thermal
    /// solver. Only the electric loss `Im(eps)` is included; magnetic loss
    /// `Im(mu)` would add a term proportional to `|H|^2`.
    pub fn power_density(&self, z: f64, incident_irradiance: f64) -> f64 {
        let e2 = self.electric_field(z).abs2();
        let im_eps = self
            .fields
            .iter()
            .find(|f| z <= f.z_start + f.thickness + 1e-15)
            .map(|f| (f.index * f.index / f.mu).im)
            .unwrap_or(0.0);
        (2.0 * PI / self.lambda) * (im_eps / self.q0) * incident_irradiance * e2
    }

    /// Cross-section-averaged heat source in layer `j`: total absorbed power
    /// per unit area divided by the layer thickness [W m^-3].
    pub fn layer_heat_source(&self, j: usize, incident_irradiance: f64) -> f64 {
        incident_irradiance * self.layer_absorptance[j] / self.fields[j].thickness
    }

    /// Sampled volumetric heat source `(z, p)` across the whole stack, ready to
    /// be interpolated or cell-averaged onto a thermal grid. `samples_per_layer`
    /// points are placed at cell centres of each layer.
    pub fn power_density_profile(
        &self,
        incident_irradiance: f64,
        samples_per_layer: usize,
    ) -> Vec<(f64, f64)> {
        assert!(samples_per_layer > 0);
        let mut out = Vec::new();
        for f in &self.fields {
            for i in 0..samples_per_layer {
                let z = f.z_start + f.thickness * (i as f64 + 0.5) / samples_per_layer as f64;
                out.push((z, self.power_density(z, incident_irradiance)));
            }
        }
        out
    }
}

/// Exact Fresnel power reflectance of a single flat interface between two
/// non-magnetic media (mu = 1).
pub fn interface_reflectance(n0: C, n1: C, theta0: f64, pol: Polarization) -> f64 {
    let s = n0.re * theta0.sin();
    let q0 = normal_wavevector(n0, s);
    let q1 = normal_wavevector(n1, s);
    let eta0 = tilted_admittance(n0, C::real(1.0), q0, pol);
    let eta1 = tilted_admittance(n1, C::real(1.0), q1, pol);
    let r = (eta0 - eta1) / (eta0 + eta1);
    r.abs2()
}
