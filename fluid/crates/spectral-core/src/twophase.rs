//! Two-phase incompressible flow (e.g. water and nitrogen gas) with surface
//! tension and gravity: Cahn–Hilliard phase field coupled to variable-density
//! Navier–Stokes, periodic box, pseudo-spectral.
//!
//! φ ∈ [−1, 1] marks fluid 1 (φ = +1) and fluid 2 (φ = −1); ρ, μ are linear
//! in φ (clipped). The interface is diffuse, thickness ~η — itself an
//! emergent model of a molecularly thin surface, valid while η is resolved
//! and small against curvature radii (the sharp-interface limit).
//!
//!   ∂t φ + u·∇φ = M Δμ_c,        μ_c = λ(−Δφ + (φ³ − φ)/η²),
//!   ∂t u + u·∇u = −∇p/ρ + ∇·(μ(∇u + ∇uᵀ))/ρ + μ_c ∇φ/ρ + (1 − ρ̄/ρ) g,
//!   ∇·u = 0,                     σ = (2√2/3) λ/η.
//!
//! Gravity enters as (1 − ρ̄/ρ)g (the mean weight is carried by the mean
//! pressure gradient in a periodic box). By default, large density ratios are
//! handled with constant-coefficient operators (Dong & Shen, J. Comput. Phys.
//! 231, 2012): (1/ρ)∇p ≈ (1/ρ₀)∇p + (1/ρ − 1/ρ₀)∇p*, ρ₀ = min ρ, and the
//! viscous term is ν₀Δu implicit plus the remainder explicit. Every solve is
//! then a division in Fourier space; the pressure step is exactly the Leray
//! projection: u^{n+1} = P Ĝ / (γ₀/Δt + ν₀|k|²), p^{n+1} = −ρ₀ i k·Ĝ/|k|².
//! [`Projection::Variable`] instead solves the variable-coefficient pressure
//! equation ∇·(ρ⁻¹∇p) = ∇·u* by PCG (constant-coefficient FFT Poisson
//! preconditioner, pseudo-spectral dealiased operator), the exact
//! variable-density projection; see `analysis/twophase_and_mhd.md`.
//! Cahn–Hilliard uses linear stabilization (Shen & Yang). Time stepping:
//! BDF2 with second-order extrapolation (BDF1 on the first step).

use num_complex::Complex64;
use rayon::prelude::*;

use crate::fft::Fft3;
use crate::grid::Grid;
use crate::solver::{project, zeros_spec, Spec};

const I: Complex64 = Complex64::new(0.0, 1.0);

/// Pressure projection used by [`TwoPhaseSolver::step`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Projection {
    /// Dong & Shen (2012) constant-coefficient splitting with ρ₀ = min ρ:
    /// the pressure gradient is removed by an FFT division with weight 1/ρ₀.
    /// Full (1/ρ − 1/ρ₀)∇p* correction at the previous pressure. This is the
    /// default and is second order, but a capillary force whose two standard
    /// forms differ by a gradient leaves a spurious solenoidal residue
    /// (~−27% capillary frequency at density ratio 10).
    Constant,
    /// Variable-coefficient projection: solve
    /// `∇·(ρ⁻¹∇p) = ∇·u*` with PCG, the operator applied pseudo-spectrally
    /// (dealiased) and the constant-coefficient FFT Poisson solve used as
    /// preconditioner, then `u = (u* − ρ⁻¹∇p)/(γ₀/Δt + ν₀|k|²)`. Reduces to
    /// the constant-coefficient projection when ρ is uniform. A gradient
    /// added to the force is absorbed exactly, so the two capillary force
    /// forms agree.
    Variable,
}

#[derive(Clone, Copy, Debug)]
pub struct TwoPhaseParams {
    /// densities of fluid 1 (φ = +1) and fluid 2 (φ = −1)
    pub rho: [f64; 2],
    /// dynamic viscosities
    pub mu: [f64; 2],
    /// surface tension
    pub sigma: f64,
    /// gravitational acceleration vector
    pub gravity: [f64; 3],
    /// interface thickness parameter η
    pub eta: f64,
    /// Cahn–Hilliard mobility M
    pub mobility: f64,
    /// stabilization constant S (Shen–Yang), ≥ ~1
    pub stab: f64,
}

pub struct TwoPhaseSolver {
    pub grid: Grid,
    pub prm: TwoPhaseParams,
    pub t: f64,
    pub u: [Spec; 3],
    pub phi: Spec,
    u_prev: [Spec; 3],
    phi_prev: Spec,
    p: Spec,
    p_prev: Spec,
    steps: usize,
    fft: Fft3,
    k2: Vec<f64>,
    lambda: f64,
    rho0: f64,
    nu0: f64,
    /// Pressure projection (default: constant-coefficient Dong–Shen).
    pub projection: Projection,
    /// PCG tolerance (relative 2-norm of the residual), variable projection.
    pub vc_tol: f64,
    /// Maximum PCG iterations, variable projection.
    pub vc_max_iters: usize,
    /// CG iterations and final relative residual of the last variable-
    /// density projection (0 when the constant projection was used).
    pub last_vc_iters: usize,
    pub last_vc_resid: f64,
}

/// `(1/ρ) ∂_d p` in spectral space for a real density field `ρ` given by its
/// reciprocal `inv_rho` on the physical grid. The gradient is formed in
/// Fourier, multiplied by `1/ρ` pointwise in physical space, and the result
/// is transformed back and dealiased. Preserves Hermitian symmetry.
fn rho_grad(fft: &Fft3, g: &Grid, inv_rho: &[f64], p: &Spec, d: usize) -> Spec {
    let ns = g.spec_len();
    let nr = g.real_len();
    let mut gp: Spec = p.par_iter().enumerate().map(|(i, v)| I * g.k(i)[d] * *v).collect();
    let mut phys = vec![0.0; nr];
    fft.inverse_inplace(&mut gp, &mut phys);
    phys.par_iter_mut().zip(inv_rho.par_iter()).for_each(|(a, b)| *a *= *b);
    let mut fs = vec![Complex64::default(); ns];
    fft.forward(&phys, &mut fs);
    fs.par_iter_mut().enumerate().for_each(|(i, v)| {
        if !g.kept(i) {
            *v = Complex64::default();
        }
    });
    fs
}

/// `A p = −∇·(ρ⁻¹∇p)`, symmetric positive semi-definite on the zero-mean
/// subspace, applied pseudo-spectrally with the same dealiasing as the rest
/// of the solver.
fn vc_apply(fft: &Fft3, g: &Grid, inv_rho: &[f64], p: &Spec, out: &mut Spec) {
    out.par_iter_mut().for_each(|v| *v = Complex64::default());
    for d in 0..3 {
        let fs = rho_grad(fft, g, inv_rho, p, d);
        out.par_iter_mut().zip(fs.par_iter()).enumerate().for_each(|(i, (o, f))| {
            if g.kept(i) {
                *o -= I * g.k(i)[d] * f;
            }
        });
    }
}

impl TwoPhaseSolver {
    /// Solve `A p = −∇·u*` for the variable-density projection, `A` as in
    /// [`vc_apply`], by preconditioned conjugate gradients. The preconditioner
    /// is the constant-coefficient FFT Poisson solve with a reference density
    /// (the volume harmonic mean of ρ, which centres the preconditioned
    /// spectrum at one). `guess` seeds the iteration (warm start). Returns
    /// `(p, iterations, relative residual)`.
    fn vc_solve(&self, gh: &[Spec; 3], inv_rho: &[f64], guess: &Spec) -> (Spec, usize, f64) {
        let g = &self.grid;
        let ns = g.spec_len();
        let k2 = &self.k2;
        // b = −∇·u*, masked to the retained modes and with zero mean.
        let mut b = vec![Complex64::default(); ns];
        for i in 0..ns {
            if !g.kept(i) {
                continue;
            }
            let k = g.k(i);
            b[i] = -(I * (gh[0][i] * k[0] + gh[1][i] * k[1] + gh[2][i] * k[2]));
        }
        b[0] = Complex64::default();
        let rho_ref = 1.0 / (inv_rho.iter().sum::<f64>() / inv_rho.len() as f64);
        // Preconditioner: ρ_ref/|k|² (k = 0 and non-retained modes → 0).
        let pre = |r: &Spec| -> Spec {
            r.par_iter()
                .enumerate()
                .map(|(i, v)| {
                    if g.kept(i) && k2[i] > 0.0 {
                        *v * (rho_ref / k2[i])
                    } else {
                        Complex64::default()
                    }
                })
                .collect()
        };
        // Weighted real inner product (half-spectrum multiplicity).
        let wdot = |a: &Spec, b: &Spec| -> f64 {
            a.par_iter()
                .zip(b.par_iter())
                .enumerate()
                .map(|(i, (x, y))| g.weight(i) * (x.conj() * y).re)
                .sum::<f64>()
        };
        let mut x: Spec = guess
            .par_iter()
            .enumerate()
            .map(|(i, v)| if g.kept(i) { *v } else { Complex64::default() })
            .collect();
        x[0] = Complex64::default();
        let mut r = b.clone();
        let mut ap = vec![Complex64::default(); ns];
        vc_apply(&self.fft, g, inv_rho, &x, &mut ap);
        r.par_iter_mut().zip(ap.par_iter()).for_each(|(a, c)| *a -= *c);
        let bnorm = wdot(&b, &b).sqrt().max(1e-300);
        let mut resid = wdot(&r, &r).sqrt() / bnorm;
        let mut z = pre(&r);
        let mut p = z.clone();
        let mut rz = wdot(&r, &z);
        let mut iters = 0usize;
        while resid > self.vc_tol && iters < self.vc_max_iters {
            vc_apply(&self.fft, g, inv_rho, &p, &mut ap);
            let pap = wdot(&p, &ap);
            if pap <= 0.0 {
                break;
            }
            let alpha = rz / pap;
            x.par_iter_mut().zip(p.par_iter()).for_each(|(a, c)| *a += alpha * c);
            r.par_iter_mut().zip(ap.par_iter()).for_each(|(a, c)| *a -= alpha * c);
            iters += 1;
            resid = wdot(&r, &r).sqrt() / bnorm;
            if resid <= self.vc_tol {
                break;
            }
            z = pre(&r);
            let rz_new = wdot(&r, &z);
            let beta = rz_new / rz;
            p.par_iter_mut().zip(z.par_iter()).for_each(|(a, c)| *a = *c + beta * *a);
            rz = rz_new;
        }
        (x, iters, resid)
    }
    pub fn new(grid: Grid, prm: TwoPhaseParams, u_phys: [Vec<f64>; 3], phi_phys: Vec<f64>) -> Self {
        let fft = Fft3::new(&grid);
        let ns = grid.spec_len();
        let k2: Vec<f64> = (0..ns)
            .map(|i| {
                let k = grid.k(i);
                k[0] * k[0] + k[1] * k[1] + k[2] * k[2]
            })
            .collect();
        let mut u = zeros_spec(ns);
        for c in 0..3 {
            fft.forward(&u_phys[c], &mut u[c]);
        }
        project(&grid, &mut u);
        let mut phi = vec![Complex64::default(); ns];
        fft.forward(&phi_phys, &mut phi);
        for i in 0..ns {
            if !grid.kept(i) {
                phi[i] = Complex64::default();
            }
        }
        let lambda = 3.0 * prm.sigma * prm.eta / (2.0 * 2f64.sqrt());
        let rho0 = prm.rho[0].min(prm.rho[1]);
        let nu0 = (prm.mu[0] / prm.rho[0]).max(prm.mu[1] / prm.rho[1]);
        let projection = match std::env::var("VCPROJ").as_deref() {
            Ok("1") | Ok("variable") => Projection::Variable,
            _ => Projection::Constant,
        };
        Self {
            u_prev: u.clone(),
            phi_prev: phi.clone(),
            p: vec![Complex64::default(); ns],
            p_prev: vec![Complex64::default(); ns],
            u,
            phi,
            t: 0.0,
            steps: 0,
            fft,
            k2,
            lambda,
            rho0,
            nu0,
            projection,
            vc_tol: 1e-8,
            vc_max_iters: 400,
            last_vc_iters: 0,
            last_vc_resid: 0.0,
            prm,
            grid,
        }
    }

    /// Select the pressure projection (see [`Projection`]).
    pub fn set_projection(&mut self, projection: Projection) {
        self.projection = projection;
    }

    fn to_phys(&self, s: &Spec) -> Vec<f64> {
        let mut w = s.clone();
        let mut out = vec![0.0; self.grid.real_len()];
        self.fft.inverse_inplace(&mut w, &mut out);
        out
    }

    fn to_spec(&self, f: &[f64]) -> Spec {
        let mut s = vec![Complex64::default(); self.grid.spec_len()];
        self.fft.forward(f, &mut s);
        s
    }

    /// ∂_d of a spectral field, in physical space.
    fn deriv(&self, s: &Spec, d: usize) -> Vec<f64> {
        let g = &self.grid;
        let ds: Spec = s.par_iter().enumerate().map(|(i, v)| I * g.k(i)[d] * *v).collect();
        self.to_phys(&ds)
    }

    fn dealias(&self, s: &mut Spec) {
        let g = &self.grid;
        s.par_iter_mut().enumerate().for_each(|(i, v)| {
            if !g.kept(i) {
                *v = Complex64::default();
            }
        });
    }

    /// Mean density ρ̄ (exact: ρ is linear in φ and φ̄ is the k = 0 mode).
    pub fn mean_density(&self) -> f64 {
        let pb = self.phi[0].re.clamp(-1.0, 1.0);
        0.5 * (self.prm.rho[0] + self.prm.rho[1]) + 0.5 * (self.prm.rho[0] - self.prm.rho[1]) * pb
    }

    /// Advance one BDF2/EXT2 step of size dt.
    pub fn step(&mut self, dt: f64) {
        let g = self.grid.clone();
        let (ns, nr) = (g.spec_len(), g.real_len());
        let first = self.steps == 0;
        let (g0, a1, a2, b1, b2) = if first { (1.0, 1.0, 0.0, 1.0, 0.0) } else { (1.5, 2.0, -0.5, 2.0, -1.0) };
        let prm = self.prm;
        let (eta, m_ch, s_ch, lam) = (prm.eta, prm.mobility, prm.stab, self.lambda);
        let k2 = self.k2.clone();
        let vc = self.projection == Projection::Variable;

        // extrapolations
        let ext = |a: &Spec, b: &Spec| -> Spec { a.par_iter().zip(b.par_iter()).map(|(x, y)| *x * b1 + *y * b2).collect() };
        let phi_s = ext(&self.phi, &self.phi_prev);
        let u_s: [Spec; 3] = std::array::from_fn(|c| ext(&self.u[c], &self.u_prev[c]));
        let p_s = ext(&self.p, &self.p_prev);
        let usp: [Vec<f64>; 3] = std::array::from_fn(|c| self.to_phys(&u_s[c]));

        // ---- Cahn–Hilliard
        let phis_p = self.to_phys(&phi_s);
        let gphi: [Vec<f64>; 3] = std::array::from_fn(|d| self.deriv(&phi_s, d));
        let adv: Vec<f64> = (0..nr).into_par_iter().map(|q| usp[0][q] * gphi[0][q] + usp[1][q] * gphi[1][q] + usp[2][q] * gphi[2][q]).collect();
        let nonlin: Vec<f64> = phis_p.par_iter().map(|&f| (f * f * f - f) / (eta * eta) - s_ch / (eta * eta) * f).collect();
        let (adv_h, nl_h) = (self.to_spec(&adv), self.to_spec(&nonlin));
        let mut phi_new: Spec = (0..ns)
            .into_par_iter()
            .map(|i| {
                let hat = self.phi[i] * a1 + self.phi_prev[i] * a2;
                let rhs = hat / dt - adv_h[i] - nl_h[i] * (m_ch * lam * k2[i]);
                rhs / (g0 / dt + m_ch * lam * (k2[i] * k2[i] + s_ch / (eta * eta) * k2[i]))
            })
            .collect();
        self.dealias(&mut phi_new);

        // ---- properties and surface tension at t^{n+1}
        let phin_p = self.to_phys(&phi_new);
        let lap_phi: Spec = phi_new.par_iter().zip(k2.par_iter()).map(|(v, &kk)| *v * kk).collect(); // −Δφ
        let lap_p = self.to_phys(&lap_phi);
        let mu_c: Vec<f64> = (0..nr).into_par_iter().map(|q| lam * (lap_p[q] + (phin_p[q].powi(3) - phin_p[q]) / (eta * eta))).collect();
        let gphin: [Vec<f64>; 3] = std::array::from_fn(|d| self.deriv(&phi_new, d));
        // capillary force per volume: μ_c∇φ (default) or −φ∇μ_c (CAPFORM=1)
        let cap_alt = std::env::var("CAPFORM").map(|v| v == "1").unwrap_or(false);
        let gmu: Option<[Vec<f64>; 3]> = cap_alt.then(|| {
            let ms = self.to_spec(&mu_c);
            std::array::from_fn(|d| self.deriv(&ms, d))
        });
        let (r1, r2) = (prm.rho[0], prm.rho[1]);
        let (m1, m2) = (prm.mu[0], prm.mu[1]);
        let rho: Vec<f64> = phin_p.par_iter().map(|&f| {
            let c = f.clamp(-1.0, 1.0);
            0.5 * (r1 + r2) + 0.5 * (r1 - r2) * c
        }).collect();
        let mu: Vec<f64> = phin_p.par_iter().map(|&f| {
            let c = f.clamp(-1.0, 1.0);
            0.5 * (m1 + m2) + 0.5 * (m1 - m2) * c
        }).collect();
        let rho_bar = {
            let pb = phi_new[0].re.clamp(-1.0, 1.0);
            0.5 * (r1 + r2) + 0.5 * (r1 - r2) * pb
        };

        // ---- momentum: explicit terms G (physical part)
        let grad_u: [[Vec<f64>; 3]; 3] = std::array::from_fn(|i| std::array::from_fn(|j| self.deriv(&u_s[i], j))); // ∂_j u_i
        // The constant-coefficient splitting carries (1/ρ − 1/ρ₀)∇p* as an
        // explicit correction; the variable-coefficient projection solves for
        // the pressure with the correct weight and needs no such term.
        let grad_p: Option<[Vec<f64>; 3]> = (!vc).then(|| std::array::from_fn(|d| self.deriv(&p_s, d)));
        // ∇·(μ D), D_ij = ∂_j u_i + ∂_i u_j
        let mut divmd: [Vec<f64>; 3] = std::array::from_fn(|_| vec![0.0; nr]);
        for i in 0..3 {
            for j in 0..3 {
                let md: Vec<f64> = (0..nr).into_par_iter().map(|q| mu[q] * (grad_u[i][j][q] + grad_u[j][i][q])).collect();
                let mds = self.to_spec(&md);
                let dj = self.deriv(&mds, j);
                divmd[i].par_iter_mut().zip(dj.par_iter()).for_each(|(a, b)| *a += b);
            }
        }
        let inv_rho0 = 1.0 / self.rho0;
        let gv = prm.gravity;
        let gphys: [Vec<f64>; 3] = std::array::from_fn(|i| {
            (0..nr)
                .into_par_iter()
                .map(|q| {
                    let conv = usp[0][q] * grad_u[i][0][q] + usp[1][q] * grad_u[i][1][q] + usp[2][q] * grad_u[i][2][q];
                    let ir = 1.0 / rho[q];
                    let pcorr = match &grad_p {
                        None => 0.0,
                        Some(gp) => (ir - inv_rho0) * gp[i][q],
                    };
                    -conv - pcorr + ir * divmd[i][q] + ir * match &gmu { None => mu_c[q] * gphin[i][q], Some(gm) => -phin_p[q] * gm[i][q] } + (1.0 - rho_bar * ir) * gv[i]
                })
                .collect()
        });
        let nu0 = self.nu0;
        let mut gh: [Spec; 3] = std::array::from_fn(|c| {
            let ph = self.to_spec(&gphys[c]);
            (0..ns)
                .into_par_iter()
                .map(|i| {
                    let hat = self.u[c][i] * a1 + self.u_prev[c][i] * a2;
                    hat / dt + ph[i] + u_s[c][i] * (nu0 * k2[i]) // −ν₀Δu* moved to the explicit side
                })
                .collect()
        });
        for c in 0..3 {
            self.dealias(&mut gh[c]);
        }
        // pressure and projection
        let p_new: Spec;
        if vc {
            let inv_rho: Vec<f64> = rho.par_iter().map(|r| 1.0 / r).collect();
            // Extrapolate the previous two pressures as the initial guess
            // (BDF2-consistent); it costs nothing and cuts CG iterations.
            let guess: Spec = self
                .p
                .par_iter()
                .zip(self.p_prev.par_iter())
                .map(|(a, b)| *a * 2.0 - *b)
                .collect();
            let (p_sol, iters, resid) = self.vc_solve(&gh, &inv_rho, &guess);
            self.last_vc_iters = iters;
            self.last_vc_resid = resid;
            // u* ← u* − ρ⁻¹∇p, in physical space so the variable weight is exact.
            for d in 0..3 {
                let corr = rho_grad(&self.fft, &g, &inv_rho, &p_sol, d);
                gh[d].par_iter_mut().zip(corr.par_iter()).for_each(|(a, c)| *a -= *c);
            }
            p_new = p_sol;
        } else {
            self.last_vc_iters = 0;
            self.last_vc_resid = 0.0;
            p_new = (0..ns)
                .into_par_iter()
                .map(|i| {
                    if k2[i] == 0.0 {
                        return Complex64::default();
                    }
                    let k = g.k(i);
                    let kg = gh[0][i] * k[0] + gh[1][i] * k[1] + gh[2][i] * k[2];
                    -I * kg * (self.rho0 / k2[i])
                })
                .collect();
            // In a periodic box ∫ρu (not ∫u) is conserved and the buoyancy term
            // would drive a spurious volume-mean drift; the mean velocity is held
            // at zero (the projection zeroes k = 0), as is usual for periodic
            // two-phase simulations.
            project(&g, &mut gh);
        }
        for c in 0..3 {
            gh[c].par_iter_mut().zip(k2.par_iter()).for_each(|(v, &kk)| *v /= g0 / dt + nu0 * kk);
        }

        self.u_prev = std::mem::replace(&mut self.u, gh);
        self.phi_prev = std::mem::replace(&mut self.phi, phi_new);
        self.p_prev = std::mem::replace(&mut self.p, p_new);
        self.steps += 1;
        self.t += dt;
    }

    /// φ in physical space.
    pub fn phi_phys(&self) -> Vec<f64> {
        self.to_phys(&self.phi)
    }

    /// Pressure (spectral potential `p`) in physical space. In the variable-
    /// density path, and at a steady state, this is the physical pressure up to
    /// an additive constant; both paths store the same convention so a pressure
    /// *jump* is comparable across them.
    pub fn pressure_phys(&self) -> Vec<f64> {
        self.to_phys(&self.p)
    }

    /// Interface height along z (the φ = 0 crossing nearest `z_guess`) in
    /// the column (i, j), by linear interpolation on a 4× refined spectral
    /// sampling of the column.
    pub fn interface_height(&self, i: usize, j: usize, z_guess: f64) -> f64 {
        let g = &self.grid;
        let n = g.n;
        let x = [g.x(i), g.x(j)];
        let fine = 4 * n;
        let col: Vec<f64> = (0..fine)
            .map(|q| {
                let z = q as f64 * 2.0 * std::f64::consts::PI / fine as f64;
                // direct spectral evaluation along the column
                let mut v = 0.0;
                for idx in 0..g.spec_len() {
                    if self.phi[idx].norm_sqr() == 0.0 {
                        continue;
                    }
                    let k = g.k(idx);
                    v += g.weight(idx) * (self.phi[idx] * Complex64::from_polar(1.0, k[0] * x[0] + k[1] * x[1] + k[2] * z)).re;
                }
                v
            })
            .collect();
        let dz = 2.0 * std::f64::consts::PI / fine as f64;
        let mut best = (f64::INFINITY, z_guess);
        for q in 0..fine {
            let (a, b) = (col[q], col[(q + 1) % fine]);
            if a == 0.0 || a * b < 0.0 {
                let z = q as f64 * dz + dz * a / (a - b);
                let d = (z - z_guess).abs();
                if d < best.0 {
                    best = (d, z);
                }
            }
        }
        best.1
    }

    /// Kinetic energy ∫ ½ ρ |u|² (box mean).
    pub fn kinetic_energy(&self) -> f64 {
        let g = &self.grid;
        let phi = self.phi_phys();
        let u: Vec<Vec<f64>> = (0..3).map(|c| self.to_phys(&self.u[c])).collect();
        let (r1, r2) = (self.prm.rho[0], self.prm.rho[1]);
        (0..g.real_len())
            .map(|q| {
                let c = phi[q].clamp(-1.0, 1.0);
                let rho = 0.5 * (r1 + r2) + 0.5 * (r1 - r2) * c;
                0.5 * rho * (u[0][q].powi(2) + u[1][q].powi(2) + u[2][q].powi(2))
            })
            .sum::<f64>()
            / g.real_len() as f64
    }

    pub fn phase_mean(&self) -> f64 {
        self.phi[0].re
    }
}

#[cfg(test)]
mod vc_absorption_tests {
    use super::*;
    use crate::init;

    fn build(n: usize, etas: f64, ratio: f64) -> TwoPhaseSolver {
        let g = Grid::new(n);
        let eta = etas * g.dx();
        let prm = TwoPhaseParams {
            rho: [1.0, 1.0 / ratio],
            mu: [2e-4, 2e-4 * 15.0 / ratio],
            sigma: 1.0,
            gravity: [0.0; 3],
            eta,
            mobility: 1e-8,
            stab: 2.0,
        };
        let a0 = 0.03;
        let k = 2.0;
        let w = 2f64.sqrt() * eta;
        let (zlo, zhi) = (0.5 * std::f64::consts::PI, 1.5 * std::f64::consts::PI);
        let phi = init::from_fn(&g, |x, _, z| {
            let zeta = a0 * (k * x).cos();
            [((z - zlo - zeta) / w).tanh() - ((z - zhi - zeta) / w).tanh() - 1.0, 0.0, 0.0]
        })[0]
            .clone();
        TwoPhaseSolver::new(g.clone(), prm, init::from_fn(&g, |_, _, _| [0.0; 3]), phi)
    }

    /// Constant and variable projections must give measurably different
    /// velocities once ρ varies. This guards against the new path silently
    /// being a no-op.
    #[test]
    fn variable_and_constant_projections_differ() {
        let mut ca = build(24, 1.2, 10.0);
        let mut vv = build(24, 1.2, 10.0);
        ca.set_projection(Projection::Constant);
        vv.set_projection(Projection::Variable);
        let dt = 0.01;
        ca.step(dt);
        vv.step(dt);
        let g = &ca.grid;
        let (mut num, mut den) = (0.0, 0.0);
        for c in 0..3 {
            for i in 0..g.spec_len() {
                num += g.weight(i) * (ca.u[c][i] - vv.u[c][i]).norm_sqr();
                den += g.weight(i) * ca.u[c][i].norm_sqr();
            }
        }
        eprintln!("relative u difference (constant vs variable): {:.3e}", (num / den).sqrt());
        assert!((num / den).sqrt() > 1e-12, "projections are identical");
    }

    /// The two capillary force forms differ by `(1/ρ)∇(φ μ_c)`, a discrete
    /// gradient only up to the truncation of the product. Measure how much of
    /// that difference the variable projection can absorb, i.e. the relative
    /// norm of the residual `δgh − ρ⁻¹∇p`. A well-resolved interface should
    /// give ~0.
    #[test]
    fn variable_projection_absorbs_force_gradient() {
        for &(n, etas) in &[(24usize, 1.2f64), (24, 2.0), (24, 3.0)] {
            let s = build(n, etas, 10.0);
            let g = &s.grid;
            let nr = g.real_len();
            let ns = g.spec_len();
            let phi = s.phi.clone();
            let k2 = &s.k2;
            let lap_phi: Spec = phi.par_iter().zip(k2.par_iter()).map(|(v, &kk)| *v * kk).collect();
            let lap_p = s.to_phys(&lap_phi);
            let phin_p = s.to_phys(&phi);
            let eta2 = s.prm.eta * s.prm.eta;
            let mu_c: Vec<f64> = (0..nr)
                .map(|q| s.lambda * (lap_p[q] + (phin_p[q].powi(3) - phin_p[q]) / eta2))
                .collect();
            let gphin: [Vec<f64>; 3] = std::array::from_fn(|d| s.deriv(&phi, d));
            let ms = s.to_spec(&mu_c);
            let gm: [Vec<f64>; 3] = std::array::from_fn(|d| s.deriv(&ms, d));
            let (r1, r2) = (s.prm.rho[0], s.prm.rho[1]);
            let rho: Vec<f64> = phin_p
                .iter()
                .map(|&f| {
                    let c = f.clamp(-1.0, 1.0);
                    0.5 * (r1 + r2) + 0.5 * (r1 - r2) * c
                })
                .collect();
            let inv_rho: Vec<f64> = rho.iter().map(|r| 1.0 / r).collect();
            // δgh = (1/ρ)(μ_c∇φ + φ∇μ_c)  (= the difference of the two forms)
            let diff: [Vec<f64>; 3] = std::array::from_fn(|i| {
                (0..nr)
                    .map(|q| (mu_c[q] * gphin[i][q] + phin_p[q] * gm[i][q]) / rho[q])
                    .collect()
            });
            let mut dh: [Spec; 3] = std::array::from_fn(|i| s.to_spec(&diff[i]));
            for c in 0..3 {
                s.dealias(&mut dh[c]);
            }
            let zero = vec![Complex64::default(); ns];
            let (p_sol, iters, resid) = s.vc_solve(&dh, &inv_rho, &zero);
            let (mut rn, mut dn) = (0.0, 0.0);
            for d in 0..3 {
                let corr = rho_grad(&s.fft, g, &inv_rho, &p_sol, d);
                for i in 0..ns {
                    let rr = dh[d][i] - corr[i];
                    rn += g.weight(i) * rr.norm_sqr();
                    dn += g.weight(i) * dh[d][i].norm_sqr();
                }
            }
            eprintln!(
                "n={n} eta/dx={etas}: absorbed residual {:.3}% (cg {iters}, r {resid:.1e})",
                100.0 * (rn / dn).sqrt()
            );
        }
    }

    /// Solenoidal (non-pressure) part left in the default force μ_c∇φ by the
    /// constant-coefficient projection (remove k·k/k²) vs the variable one.
    #[test]
    fn solenoidal_residue_constant_vs_variable() {
        for &ratio in &[1.0f64, 10.0, 100.0] {
            let s = build(24, 1.2, ratio);
            let g = &s.grid;
            let nr = g.real_len();
            let ns = g.spec_len();
            let phi = s.phi.clone();
            let k2 = &s.k2;
            let lap_phi: Spec = phi.par_iter().zip(k2.par_iter()).map(|(v, &kk)| *v * kk).collect();
            let lap_p = s.to_phys(&lap_phi);
            let phin_p = s.to_phys(&phi);
            let eta2 = s.prm.eta * s.prm.eta;
            let mu_c: Vec<f64> = (0..nr)
                .map(|q| s.lambda * (lap_p[q] + (phin_p[q].powi(3) - phin_p[q]) / eta2))
                .collect();
            let gphin: [Vec<f64>; 3] = std::array::from_fn(|d| s.deriv(&phi, d));
            let (r1, r2) = (s.prm.rho[0], s.prm.rho[1]);
            let rho: Vec<f64> = phin_p
                .iter()
                .map(|&f| {
                    let c = f.clamp(-1.0, 1.0);
                    0.5 * (r1 + r2) + 0.5 * (r1 - r2) * c
                })
                .collect();
            let inv_rho: Vec<f64> = rho.iter().map(|r| 1.0 / r).collect();
            let fphys: [Vec<f64>; 3] = std::array::from_fn(|i| {
                (0..nr).map(|q| mu_c[q] * gphin[i][q] / rho[q]).collect()
            });
            let mut fh: [Spec; 3] = std::array::from_fn(|i| s.to_spec(&fphys[i]));
            for c in 0..3 {
                s.dealias(&mut fh[c]);
            }
            // constant projection: remove k(k·f)/k²
            let mut solc: [Spec; 3] = fh.clone();
            for i in 0..ns {
                if k2[i] == 0.0 {
                    continue;
                }
                let k = g.k(i);
                let kf = (fh[0][i] * k[0] + fh[1][i] * k[1] + fh[2][i] * k[2]) / k2[i];
                for d in 0..3 {
                    solc[d][i] -= kf * k[d];
                }
            }
            // variable projection
            let zero = vec![Complex64::default(); ns];
            let (p_sol, _, _) = s.vc_solve(&fh, &inv_rho, &zero);
            let mut solv: [Spec; 3] = fh.clone();
            for d in 0..3 {
                let corr = rho_grad(&s.fft, g, &inv_rho, &p_sol, d);
                for i in 0..ns {
                    solv[d][i] -= corr[i];
                }
            }
            let (mut nc, mut nv, mut nf) = (0.0, 0.0, 0.0);
            for d in 0..3 {
                for i in 0..ns {
                    nf += g.weight(i) * fh[d][i].norm_sqr();
                    nc += g.weight(i) * solc[d][i].norm_sqr();
                    nv += g.weight(i) * solv[d][i].norm_sqr();
                }
            }
            eprintln!(
                "ratio={ratio}: solenoidal part of μ_c∇φ: constant {:.3}% variable {:.3}%",
                100.0 * (nc / nf).sqrt(),
                100.0 * (nv / nf).sqrt()
            );
        }
    }

    /// Project the two force forms separately and compare the *projected*
    /// (solenoidal) fields. If the absorption test is right this should be at
    /// the aliasing level, not tens of percent.
    #[test]
    fn projected_force_forms_agree() {
        let s = build(24, 1.2, 10.0);
        let g = &s.grid;
        let nr = g.real_len();
        let ns = g.spec_len();
        let phi = s.phi.clone();
        let k2 = &s.k2;
        let lap_phi: Spec = phi.par_iter().zip(k2.par_iter()).map(|(v, &kk)| *v * kk).collect();
        let lap_p = s.to_phys(&lap_phi);
        let phin_p = s.to_phys(&phi);
        let eta2 = s.prm.eta * s.prm.eta;
        let mu_c: Vec<f64> = (0..nr)
            .map(|q| s.lambda * (lap_p[q] + (phin_p[q].powi(3) - phin_p[q]) / eta2))
            .collect();
        let gphin: [Vec<f64>; 3] = std::array::from_fn(|d| s.deriv(&phi, d));
        let ms = s.to_spec(&mu_c);
        let gm: [Vec<f64>; 3] = std::array::from_fn(|d| s.deriv(&ms, d));
        let (r1, r2) = (s.prm.rho[0], s.prm.rho[1]);
        let rho: Vec<f64> = phin_p
            .iter()
            .map(|&f| {
                let c = f.clamp(-1.0, 1.0);
                0.5 * (r1 + r2) + 0.5 * (r1 - r2) * c
            })
            .collect();
        let inv_rho: Vec<f64> = rho.iter().map(|r| 1.0 / r).collect();
        let fd: [Vec<f64>; 3] = std::array::from_fn(|i| {
            (0..nr).map(|q| mu_c[q] * gphin[i][q] / rho[q]).collect()
        });
        let fa: [Vec<f64>; 3] = std::array::from_fn(|i| {
            (0..nr).map(|q| -phin_p[q] * gm[i][q] / rho[q]).collect()
        });
        let mut fdh: [Spec; 3] = std::array::from_fn(|i| s.to_spec(&fd[i]));
        let mut fah: [Spec; 3] = std::array::from_fn(|i| s.to_spec(&fa[i]));
        for c in 0..3 {
            s.dealias(&mut fdh[c]);
            s.dealias(&mut fah[c]);
        }
        let zero = vec![Complex64::default(); ns];
        let (pd, _, _) = s.vc_solve(&fdh, &inv_rho, &zero);
        let (pa, _, _) = s.vc_solve(&fah, &inv_rho, &zero);
        let (mut nd, mut na) = (0.0, 0.0);
        for d in 0..3 {
            let cd = rho_grad(&s.fft, g, &inv_rho, &pd, d);
            let ca = rho_grad(&s.fft, g, &inv_rho, &pa, d);
            for i in 0..ns {
                let vd = fdh[d][i] - cd[i];
                let va = fah[d][i] - ca[i];
                nd += g.weight(i) * vd.norm_sqr();
                na += g.weight(i) * (vd - va).norm_sqr();
            }
        }
        eprintln!(
            "projected form difference relative to projected default: {:.3}",
            (na / nd).sqrt()
        );
        assert!((na / nd).sqrt() < 0.15, "projected force forms differ by {:.1}%", 100.0 * (na / nd).sqrt());
    }
}
