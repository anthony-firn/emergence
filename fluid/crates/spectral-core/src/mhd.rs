//! Quasi-static (low magnetic Reynolds number) magnetohydrodynamics for a
//! conducting liquid such as mercury or liquid sodium.
//!
//! Model ladder. Full Maxwell + Navier–Stokes reduces in two emergent steps:
//! 1. magnetoquasistatics, when the flow is slow against light
//!    (L ≪ c/f, charge relaxation ε/σ ≪ flow time — 1e-19 s for mercury):
//!    ∇·J = s, J = σ(E + u×B), E = −∇φ (+ a uniform applied field);
//! 2. low Rm = μ₀σUL ≪ 1: the field induced by J is negligible against the
//!    applied B, so B is a given field and the induction equation drops out.
//!
//! With uniform σ, ∇·J = s gives J in closed form on the torus:
//!   J = σ P[u×B] + ∇Δ⁻¹s + J̄,
//! where P is the Leray projector (the electrostatic potential exactly
//! cancels the gradient part of the motional EMF u×B), s the current source
//! density of the electrodes, and J̄ the mean current density (zero for an
//! open circuit). The Lorentz force per unit mass is J×B/ρ.
//!
//! For a uniform B the induced part is diagonal in Fourier space:
//!   P[(σ/ρ) P(û×B) × B] = −(σ/ρ)(k̂·B)² û
//! (Joule damping acts only on modes varying along B). That exact linear
//! rate can be folded into the integrating factor, removing the Hartmann
//! stiffness instead of resolving it with small steps.
//!
//! A permanent magnet with periodic magnetization M has B = μ₀(M + H),
//! H = −∇ψ, ∇·B = 0, hence B = μ₀ P[M] + B̄: its field is the Leray
//! projection of the magnetization.

use crate::fft::Fft3;
use crate::grid::Grid;
use crate::solver::Spec;
use num_complex::Complex64;
use rayon::prelude::*;

const I: Complex64 = Complex64 { re: 0.0, im: 1.0 };

#[derive(Clone, Debug)]
pub struct Mhd {
    /// σ/ρ: the induced current's force coefficient.
    pub sigma_over_rho: f64,
    /// 1/ρ: force coefficient of the applied (electrode / mean) currents.
    pub inv_rho: f64,
    /// Uniform applied field B̄.
    pub b0: [f64; 3],
    /// Additional non-uniform field (physical, N³ per component), e.g. from
    /// [`magnet_field`]. `None`: B = B̄.
    pub b_var: Option<[Vec<f64>; 3]>,
    /// Current source density s = ∇·J (physical, N³, zero mean): positive at
    /// the anode, negative at the cathode.
    pub source: Option<Vec<f64>>,
    /// Mean current density J̄ (zero: open circuit).
    pub j_mean: [f64; 3],
    /// Treat the Joule damping of a uniform field exactly in the integrating
    /// factor (IF-RK4 only; requires `b_var == None`).
    pub joule_in_if: bool,
    /// Fluid indicator (physical, N³; 1 in the liquid, 0 in solids). With
    /// marker walls the velocity inside a body is a fictitious smooth
    /// extension; solids are conductors at rest, so their motional EMF must
    /// be masked out. `None`: the whole box is liquid.
    pub fluid_mask: Option<Vec<f64>>,
}

impl Mhd {
    pub fn uniform(sigma_over_rho: f64, inv_rho: f64, b0: [f64; 3]) -> Self {
        Self { sigma_over_rho, inv_rho, b0, b_var: None, source: None, j_mean: [0.0; 3], joule_in_if: false, fluid_mask: None }
    }

    /// Exact Joule damping rate (σ/ρ)(k̂·B̄)² of mode k (uniform field).
    #[inline]
    pub fn joule_rate(&self, k: [f64; 3]) -> f64 {
        let k2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
        if k2 == 0.0 {
            return 0.0;
        }
        let kb = k[0] * self.b0[0] + k[1] * self.b0[1] + k[2] * self.b0[2];
        self.sigma_over_rho * kb * kb / k2
    }

    fn induced_in_if(&self) -> bool {
        self.joule_in_if && self.b_var.is_none() && self.fluid_mask.is_none()
    }

    /// Current density J (spectral) for the velocity û; the induced part is
    /// left out when it is integrated exactly.
    pub fn current(&self, g: &Grid, fft: &Fft3, uh: &[Spec; 3]) -> [Spec; 3] {
        let ns = g.spec_len();
        let mut j: [Spec; 3] = std::array::from_fn(|_| vec![Complex64::default(); ns]);
        if !self.induced_in_if() && self.sigma_over_rho != 0.0 {
            // w = u × B
            let w: [Spec; 3] = match (&self.b_var, &self.fluid_mask) {
                (None, None) => {
                    let b = self.b0;
                    let (a0, a1, a2) = (&uh[0], &uh[1], &uh[2]);
                    let mut w: [Spec; 3] = std::array::from_fn(|_| vec![Complex64::default(); ns]);
                    let [w0, w1, w2] = &mut w;
                    w0.par_iter_mut().zip(w1.par_iter_mut()).zip(w2.par_iter_mut()).enumerate().for_each(|(i, ((x, y), z))| {
                        *x = a1[i] * b[2] - a2[i] * b[1];
                        *y = a2[i] * b[0] - a0[i] * b[2];
                        *z = a0[i] * b[1] - a1[i] * b[0];
                    });
                    w
                }
                _ => {
                    let mut u: Vec<Vec<f64>> = (0..3).map(|c| to_phys(g, fft, &uh[c])).collect();
                    if let Some(m) = &self.fluid_mask {
                        for uc in u.iter_mut() {
                            uc.par_iter_mut().zip(m.par_iter()).for_each(|(v, w)| *v *= w);
                        }
                    }
                    let b = self.field(g);
                    let wp: [Vec<f64>; 3] = std::array::from_fn(|c| {
                        let (a, bb) = ((c + 1) % 3, (c + 2) % 3);
                        (0..g.real_len()).into_par_iter().map(|q| u[a][q] * b[bb][q] - u[bb][q] * b[a][q]).collect()
                    });
                    std::array::from_fn(|c| to_spec(g, fft, &wp[c]))
                }
            };
            let s = self.sigma_over_rho / self.inv_rho; // σ
            for c in 0..3 {
                j[c].par_iter_mut().zip(w[c].par_iter()).for_each(|(o, v)| *o = *v * s);
            }
            leray(g, &mut j);
        }
        if let Some(src) = &self.source {
            let sh = to_spec(g, fft, src);
            for c in 0..3 {
                j[c].par_iter_mut().enumerate().for_each(|(i, o)| {
                    let k = g.k(i);
                    let k2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
                    if k2 > 0.0 {
                        *o -= I * k[c] * sh[i] / k2; // ∇Δ⁻¹ s
                    }
                });
            }
        }
        for c in 0..3 {
            j[c][0] = Complex64::new(self.j_mean[c], 0.0);
        }
        j
    }

    /// Total field B (physical, N³ per component).
    pub fn field(&self, g: &Grid) -> [Vec<f64>; 3] {
        std::array::from_fn(|c| match &self.b_var {
            None => vec![self.b0[c]; g.real_len()],
            Some(bv) => bv[c].iter().map(|v| v + self.b0[c]).collect(),
        })
    }

    /// Add the Lorentz force J×B/ρ (spectral, unprojected) to `out`.
    pub fn add_force(&self, g: &Grid, fft: &Fft3, uh: &[Spec; 3], out: &mut [Spec; 3]) {
        let j = self.current(g, fft, uh);
        let ir = self.inv_rho;
        match &self.b_var {
            None => {
                let b = self.b0;
                for c in 0..3 {
                    let (a, bb) = ((c + 1) % 3, (c + 2) % 3);
                    let (ja, jb) = (&j[a], &j[bb]);
                    out[c].par_iter_mut().enumerate().for_each(|(i, o)| *o += (ja[i] * b[bb] - jb[i] * b[a]) * ir);
                }
            }
            Some(_) => {
                let b = self.field(g);
                let jp: Vec<Vec<f64>> = (0..3).map(|c| to_phys(g, fft, &j[c])).collect();
                for c in 0..3 {
                    let (a, bb) = ((c + 1) % 3, (c + 2) % 3);
                    let f: Vec<f64> = (0..g.real_len()).into_par_iter().map(|q| (jp[a][q] * b[bb][q] - jp[bb][q] * b[a][q]) * ir).collect();
                    let fh = to_spec(g, fft, &f);
                    out[c].par_iter_mut().zip(fh.par_iter()).for_each(|(o, v)| *o += *v);
                }
            }
        }
    }

    /// Joule dissipation rate ∫ |J|²/σ per unit mass (box mean of
    /// |J|²·(1/ρσ)), and the electric power input of the sources.
    pub fn joule_dissipation(&self, g: &Grid, fft: &Fft3, uh: &[Spec; 3]) -> f64 {
        let mut s = self.clone();
        s.joule_in_if = false;
        let j = s.current(g, fft, uh);
        let e: f64 = (0..3).map(|c| (0..g.spec_len()).map(|i| g.weight(i) * j[c][i].norm_sqr()).sum::<f64>()).sum();
        e * self.inv_rho * self.inv_rho / self.sigma_over_rho
    }
}

/// Field of a permanent magnet with periodic magnetization μ₀M (physical,
/// N³ per component, in field units): B = P[μ₀M], whose mean is μ₀M̄
/// (the mean of H vanishes on the torus). Any external uniform field goes
/// in B̄.
pub fn magnet_field(g: &Grid, mu0_m: &[Vec<f64>; 3]) -> [Vec<f64>; 3] {
    let fft = Fft3::new(g);
    let mut s: [Spec; 3] = std::array::from_fn(|c| to_spec(g, &fft, &mu0_m[c]));
    leray(g, &mut s);
    std::array::from_fn(|c| to_phys(g, &fft, &s[c]))
}

/// Smooth radial current between coaxial cylindrical electrodes along z:
/// current I (per box height 2π) enters on a ring of radius `r_in` and
/// leaves on a ring of radius `r_out`, each spread over a compact C^∞ bump
/// of half-width `w` placed on the solid side of the wall. Rings are
/// normalized discretely, so ∫ s = 0 exactly and each carries ±I.
pub fn coaxial_source(g: &Grid, center: [f64; 2], r_in: f64, r_out: f64, w: f64, current: f64) -> Vec<f64> {
    let n = g.n;
    let bump = |d: f64| if d.abs() < 1.0 { (-1.0 / (1.0 - d * d)).exp() } else { 0.0 };
    let ring = |rc: f64| -> Vec<f64> {
        (0..g.real_len())
            .into_par_iter()
            .map(|p| {
                let (i, j) = (p / (n * n), (p / n) % n);
                let r = ((g.x(i) - center[0]).powi(2) + (g.x(j) - center[1]).powi(2)).sqrt();
                bump((r - rc) / w)
            })
            .collect()
    };
    let (a, b) = (ring(r_in - w), ring(r_out + w));
    let dv = g.dx().powi(3);
    let (sa, sb) = (a.iter().sum::<f64>() * dv, b.iter().sum::<f64>() * dv);
    a.iter().zip(b.iter()).map(|(x, y)| current * (x / sa - y / sb)).collect()
}

fn leray(g: &Grid, v: &mut [Spec; 3]) {
    let [a, b, c] = v;
    a.par_iter_mut().zip(b.par_iter_mut()).zip(c.par_iter_mut()).enumerate().for_each(|(idx, ((x, y), z))| {
        let k = g.k(idx);
        let k2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
        if !g.kept(idx) {
            *x = Complex64::default();
            *y = Complex64::default();
            *z = Complex64::default();
            return;
        }
        if k2 == 0.0 {
            return;
        }
        let kd = (*x * k[0] + *y * k[1] + *z * k[2]) / k2;
        *x -= kd * k[0];
        *y -= kd * k[1];
        *z -= kd * k[2];
    });
}

fn to_phys(g: &Grid, fft: &Fft3, s: &Spec) -> Vec<f64> {
    let mut w = s.clone();
    let mut out = vec![0.0; g.real_len()];
    fft.inverse_inplace(&mut w, &mut out);
    out
}

fn to_spec(g: &Grid, fft: &Fft3, f: &[f64]) -> Spec {
    let mut s = vec![Complex64::default(); g.spec_len()];
    fft.forward(f, &mut s);
    s
}
