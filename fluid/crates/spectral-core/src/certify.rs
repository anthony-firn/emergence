//! A posteriori error certificate for unforced, unpenalized runs.
//!
//! Let v(t) be the numerical solution, made a continuous-in-time function by
//! cubic Hermite interpolation between steps (nodes vⁿ, tendencies Fⁿ), and
//! let R = ∂t v − νΔv + P(v·∇v) be its residual in the exact Navier–Stokes
//! equations (Leray-projected; the pressure absorbs the gradient part). For
//! any Leray–Hopf weak solution u with the same data, w = u − v satisfies
//! (weak–strong stability; all constants explicit):
//!
//!   d/dt ‖w‖² ≤ 2 μ(t) ‖w‖² + ‖R‖²_{H⁻¹} / (2ν),
//!   μ(t) = max_x ( −λ_min(sym ∇v(x, t)) ),
//!
//! from −∫ w·(∇v)w ≤ μ‖w‖² and |∫ R·w| ≤ ‖R‖_{H⁻¹}‖∇w‖ ≤ ν‖∇w‖² + ‖R‖²_{H⁻¹}/(4ν).
//! Norms are box means, ⟨·⟩ = (2π)⁻³ ∫ ·, consistently on both sides.
//!
//! The residual has two orthogonal parts:
//! * truncation: the modes of P(v × ω) that the 2/3 rule discards. v has
//!   modes |k_i| ≤ k_max, so v × ω has |k_i| ≤ 2k_max and is computed
//!   *exactly* on a padded grid of M ≥ 4k_max + 1 points per direction.
//! * time: R_time = v̇ − F(v) for the interpolant, F the discrete tendency.
//!   It vanishes at the nodes (v̇ = Fⁿ there), so it is sampled at the
//!   midpoint and integrated by Simpson's rule.
//!
//! Caveats, stated plainly: μ is the maximum over the padded grid's points
//! (a sampled estimate of the supremum, on a grid 1.5× finer than the
//! solver's), held at its maximum over each step's three samples; time
//! integrals use Simpson quadrature; residuals are evaluated in floating
//! point. Within those, ‖u − v‖ ≤ √y(t) holds for
//! every Leray–Hopf solution, so it also certifies closeness for the true
//! flow. The factor exp(2∫μ) makes the bound tight only for moderate
//! times / Reynolds numbers; for long turbulent runs it is honest but large.

use num_complex::Complex64;
use rayon::prelude::*;

use crate::fft::Fft3;
use crate::grid::{wavenumber, Grid};
use crate::solver::{project, zeros_spec, Solver, Spec};

const I: Complex64 = Complex64::new(0.0, 1.0);

/// Diagnostics of one certified step.
#[derive(Clone, Copy, Debug, Default)]
pub struct CertStep {
    pub t: f64,
    /// ‖P R_trunc‖_{H⁻¹} at the step's midpoint
    pub resid_trunc: f64,
    /// ‖R_time‖_{H⁻¹} at the step's midpoint
    pub resid_time: f64,
    /// sampled max of −λ_min(sym ∇v) over the step (endpoints and midpoint)
    pub mu: f64,
    /// certified bound on ‖u − v‖ (box-mean L²) at t
    pub bound: f64,
}

/// Exact evaluation of v × ω on a padded grid.
struct Padded {
    gp: Grid,
    fft: Fft3,
    /// padded spectral index → N-grid spectral index, for retained modes
    from_n: Vec<Option<usize>>,
    /// N-grid retained modes as (padded index, N index)
    pairs: Vec<(usize, usize)>,
    spec: [Spec; 3],
    phys_u: [Vec<f64>; 3],
    phys_w: [Vec<f64>; 3],
    tmp: Spec,
}

impl Padded {
    fn new(gn: &Grid) -> Self {
        let kmax = gn.kmax();
        // even, FFT-friendly, ≥ 4 k_max + 2 (3N/2 satisfies this for even N)
        let mut m = (3 * gn.n) / 2;
        if m < 4 * kmax + 2 {
            m = 4 * kmax + 2;
        }
        m += m % 2;
        let gp = Grid::new(m);
        let mut from_n = vec![None; gp.spec_len()];
        let mut pairs = Vec::new();
        for idx in 0..gp.spec_len() {
            let (i, j, l) = gp.spec_ijl(idx);
            let (kx, ky, kz) = (wavenumber(i, m), wavenumber(j, m), l as i64);
            let km = kmax as i64;
            if kx.abs() <= km && ky.abs() <= km && kz <= km {
                let wrap = |k: i64| k.rem_euclid(gn.n as i64) as usize;
                let ni = (wrap(kx) * gn.n + wrap(ky)) * gn.nh + l;
                from_n[idx] = Some(ni);
                pairs.push((idx, ni));
            }
        }
        let (ns, nr) = (gp.spec_len(), gp.real_len());
        Self {
            fft: Fft3::new(&gp),
            from_n,
            pairs,
            spec: zeros_spec(ns),
            phys_u: std::array::from_fn(|_| vec![0.0; nr]),
            phys_w: std::array::from_fn(|_| vec![0.0; nr]),
            tmp: vec![Complex64::default(); ns],
            gp,
        }
    }

    /// Returns (‖P(v×ω)|discarded‖²_{H⁻¹}, max_x −λ_min(sym ∇v)).
    fn evaluate(&mut self, uh: &[Spec; 3]) -> (f64, f64) {
        // μ = max −λ_min of the strain tensor on the padded grid
        let mu = self.max_neg_strain(uh);
        let gp = &self.gp;
        // velocity and vorticity on the padded grid
        for c in 0..3 {
            self.spec[c].iter_mut().for_each(|v| *v = Complex64::default());
            for &(pi, ni) in &self.pairs {
                self.spec[c][pi] = uh[c][ni];
            }
        }
        for c in 0..3 {
            self.fft.inverse(&self.spec[c], &mut self.tmp, &mut self.phys_u[c]);
        }
        for c in 0..3 {
            let (a, b) = ((c + 1) % 3, (c + 2) % 3);
            let (sa, sb) = (&self.spec[a], &self.spec[b]);
            self.tmp.par_iter_mut().enumerate().for_each(|(idx, w)| {
                let k = gp.k(idx);
                *w = I * (sb[idx] * k[a] - sa[idx] * k[b]);
            });
            self.fft.inverse_inplace(&mut self.tmp, &mut self.phys_w[c]);
        }
        // product v × ω exactly, then its discarded, projected modes
        let mut prod: [Vec<f64>; 3] = std::array::from_fn(|_| vec![0.0; gp.real_len()]);
        {
            let [u0, u1, u2] = &self.phys_u;
            let [w0, w1, w2] = &self.phys_w;
            let [p0, p1, p2] = &mut prod;
            p0.par_iter_mut().zip(p1.par_iter_mut()).zip(p2.par_iter_mut()).enumerate().for_each(
                |(q, ((a, b), c))| {
                    *a = u1[q] * w2[q] - u2[q] * w1[q];
                    *b = u2[q] * w0[q] - u0[q] * w2[q];
                    *c = u0[q] * w1[q] - u1[q] * w0[q];
                },
            );
        }
        for c in 0..3 {
            self.fft.forward(&prod[c], &mut self.spec[c]);
        }
        // Leray projection of the discarded part (the product has no aliasing
        // on this grid); `project` would also zero our high modes via the
        // padded grid's own 2/3 mask, so project by hand.
        let from_n = &self.from_n;
        let [a, b, c] = &self.spec;
        let h1 = (0..gp.spec_len())
            .into_par_iter()
            .map(|idx| {
                if from_n[idx].is_some() {
                    return 0.0; // retained mode: part of the discrete equations
                }
                let k = gp.k(idx);
                let k2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
                if k2 == 0.0 {
                    return 0.0;
                }
                let d = (a[idx] * k[0] + b[idx] * k[1] + c[idx] * k[2]) / k2;
                let p = [a[idx] - d * k[0], b[idx] - d * k[1], c[idx] - d * k[2]];
                gp.weight(idx) * (p[0].norm_sqr() + p[1].norm_sqr() + p[2].norm_sqr()) / k2
            })
            .sum::<f64>();
        (h1, mu)
    }

    /// max over padded grid points of −λ_min(sym ∇v), λ_min in closed form.
    fn max_neg_strain(&mut self, uh: &[Spec; 3]) -> f64 {
        let gp = &self.gp;
        let nr = gp.real_len();
        // six independent strain components S_ij = (∂_i u_j + ∂_j u_i)/2
        let comps = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];
        let mut s: Vec<Vec<f64>> = Vec::with_capacity(6);
        for &(i, j) in &comps {
            self.tmp.iter_mut().for_each(|v| *v = Complex64::default());
            for &(pi, ni) in &self.pairs {
                let k = gp.k(pi);
                self.tmp[pi] = 0.5 * I * (uh[j][ni] * k[i] + uh[i][ni] * k[j]);
            }
            let mut out = vec![0.0; nr];
            self.fft.inverse_inplace(&mut self.tmp, &mut out);
            s.push(out);
        }
        (0..nr)
            .into_par_iter()
            .map(|q| {
                let (a, b, c) = (s[0][q], s[1][q], s[2][q]);
                let (d, e, f) = (s[3][q], s[4][q], s[5][q]);
                -sym3_min_eig(a, b, c, d, e, f)
            })
            .reduce(|| f64::NEG_INFINITY, f64::max)
            .max(0.0)
    }
}

/// Smallest eigenvalue of the symmetric matrix [[a,d,e],[d,b,f],[e,f,c]]
/// (trigonometric solution of the characteristic cubic).
pub fn sym3_min_eig(a: f64, b: f64, c: f64, d: f64, e: f64, f: f64) -> f64 {
    let p1 = d * d + e * e + f * f;
    let q = (a + b + c) / 3.0;
    let p2 = (a - q).powi(2) + (b - q).powi(2) + (c - q).powi(2) + 2.0 * p1;
    if p2 <= 1e-300 {
        return q;
    }
    let p = (p2 / 6.0).sqrt();
    let (ba, bb, bc) = ((a - q) / p, (b - q) / p, (c - q) / p);
    let (bd, be, bf) = (d / p, e / p, f / p);
    let det = ba * (bb * bc - bf * bf) - bd * (bd * bc - bf * be) + be * (bd * bf - bb * be);
    let r = (det / 2.0).clamp(-1.0, 1.0);
    let phi = r.acos() / 3.0;
    q + 2.0 * p * (phi + 2.0 * std::f64::consts::PI / 3.0).cos()
}

/// Running certificate state, attached to a solver.
pub struct Certifier {
    padded: Padded,
    prev_u: Option<[Spec; 3]>,
    prev_f: Option<[Spec; 3]>,
    prev_mu: f64,
    prev_r2: f64,
    /// y(t) = certified bound on ‖u − v‖²
    pub y: f64,
    pub history: Vec<CertStep>,
}

impl Certifier {
    /// `initial_error`: a bound on ‖u(0) − v(0)‖ (zero when the initial data
    /// is exactly representable, e.g. Taylor–Green or ABC flows).
    pub fn new(grid: &Grid, initial_error: f64) -> Self {
        Self {
            padded: Padded::new(grid),
            prev_u: None,
            prev_f: None,
            prev_mu: 0.0,
            prev_r2: 0.0,
            y: initial_error * initial_error,
            history: Vec::new(),
        }
    }
}

impl Solver {
    /// Discrete tendency F(v) = P[(v × ω)^] − ν|k|² v̂ (no forcing/penalization).
    fn tendency(&mut self, uh: &[Spec; 3]) -> [Spec; 3] {
        let mut f = zeros_spec(self.grid.spec_len());
        self.rhs(uh, &mut f);
        let nu = self.nu;
        for c in 0..3 {
            f[c].par_iter_mut()
                .zip(uh[c].par_iter().zip(self.k2.par_iter()))
                .for_each(|(fv, (u, &k2))| *fv -= *u * (nu * k2));
        }
        f
    }

    /// Advance one IF-RK4 step and extend the certificate over it.
    /// Requires no forcing and no penalization (the bound is for the
    /// unforced Navier–Stokes equations).
    pub fn step_certified(&mut self, h: f64, cert: &mut Certifier) -> CertStep {
        assert!(self.forcing.is_none() && self.penalization.is_none() && self.body_force.is_none());
        assert!(!self.grid.spherical, "the certificate assumes the exact 2/3-rule discretization");
        let nu = self.nu;
        if cert.prev_u.is_none() {
            let u0 = self.u.clone();
            let (r2, mu) = cert.padded.evaluate(&u0);
            cert.prev_f = Some(self.tendency(&u0));
            cert.prev_u = Some(u0);
            cert.prev_mu = mu;
            cert.prev_r2 = r2;
        }
        let (u0, f0) = (cert.prev_u.take().unwrap(), cert.prev_f.take().unwrap());
        self.step_ifrk4(h);
        let u1 = self.u.clone();
        let f1 = self.tendency(&u1);

        // Hermite midpoint: v = (v0+v1)/2 + h/8 (F0 − F1),  v̇ = 3/(2h)(v1 − v0) − (F0 + F1)/4
        let ns = self.grid.spec_len();
        let mut vm = zeros_spec(ns);
        let mut vdot = zeros_spec(ns);
        for c in 0..3 {
            vm[c].par_iter_mut()
                .zip(vdot[c].par_iter_mut())
                .enumerate()
                .for_each(|(i, (m, d))| {
                    *m = 0.5 * (u0[c][i] + u1[c][i]) + (h / 8.0) * (f0[c][i] - f1[c][i]);
                    *d = (1.5 / h) * (u1[c][i] - u0[c][i]) - 0.25 * (f0[c][i] + f1[c][i]);
                });
        }
        let fm = self.tendency(&vm);
        // R_time on retained modes (already solenoidal): ‖v̇ − F(v)‖²_{H⁻¹}
        let mut rt = zeros_spec(ns);
        for c in 0..3 {
            rt[c].par_iter_mut()
                .enumerate()
                .for_each(|(i, r)| *r = vdot[c][i] - fm[c][i]);
        }
        project(&self.grid, &mut rt);
        let g = &self.grid;
        let k2 = &self.k2;
        let rtime2: f64 = (0..ns)
            .into_par_iter()
            .map(|i| {
                if k2[i] == 0.0 {
                    return 0.0;
                }
                g.weight(i) * (rt[0][i].norm_sqr() + rt[1][i].norm_sqr() + rt[2][i].norm_sqr()) / k2[i]
            })
            .sum();
        let (rtr_m, mu_m) = cert.padded.evaluate(&vm);
        let (rtr_1, mu_1) = cert.padded.evaluate(&u1);

        // ∫‖R‖² by Simpson: endpoints carry only truncation residual
        let r2_mid = rtime2 + rtr_m;
        let int_r2 = h / 6.0 * (cert.prev_r2 + 4.0 * r2_mid + rtr_1);
        let mu = cert.prev_mu.max(mu_m).max(mu_1);
        // y' ≤ 2μ y + r/(2ν), μ held at its step maximum
        let grow = (2.0 * mu * h).exp();
        cert.y = grow * (cert.y + int_r2 / (2.0 * nu));

        cert.prev_u = Some(u1);
        cert.prev_f = Some(f1);
        cert.prev_mu = mu_1;
        cert.prev_r2 = rtr_1;
        let st = CertStep {
            t: self.t,
            resid_trunc: rtr_m.sqrt(),
            resid_time: rtime2.sqrt(),
            mu,
            bound: cert.y.sqrt(),
        };
        cert.history.push(st);
        st
    }
}

// ---------------------------------------------------------------------------
// Certified adaptive resolution ("fold/unfold" of the whole field)

impl Solver {
    /// The same flow on an `n_new`³ grid: retained modes are copied, modes the
    /// new grid cannot hold are dropped. Returns the new solver and the exact
    /// box-mean L² norm of what was dropped (zero when refining).
    pub fn regrid(&self, n_new: usize) -> (Solver, f64) {
        let gn = Grid::new(n_new);
        let mut s = Solver::new(gn.clone(), self.nu, std::array::from_fn(|_| vec![0.0; gn.real_len()]));
        s.t = self.t;
        let go = &self.grid;
        let km = gn.kmax() as i64;
        let mut dropped2 = 0.0;
        for idx in 0..go.spec_len() {
            if !go.kept(idx) {
                continue;
            }
            let (i, j, l) = go.spec_ijl(idx);
            let (kx, ky, kz) = (wavenumber(i, go.n), wavenumber(j, go.n), l as i64);
            if kx.abs() <= km && ky.abs() <= km && kz <= km {
                let w = |k: i64| k.rem_euclid(n_new as i64) as usize;
                let ni = (w(kx) * n_new + w(ky)) * gn.nh + l;
                for c in 0..3 {
                    s.u[c][ni] = self.u[c][idx];
                }
            } else {
                dropped2 += go.weight(idx) * (0..3).map(|c| self.u[c][idx].norm_sqr()).sum::<f64>();
            }
        }
        (s, dropped2.sqrt())
    }
}

impl Certifier {
    /// Continue the certificate on a new grid after `Solver::regrid`; the
    /// dropped norm enters the bound by the triangle inequality.
    pub fn rebind(&mut self, grid: &Grid, dropped: f64) {
        self.padded = Padded::new(grid);
        self.prev_u = None;
        self.prev_f = None;
        self.y = (self.y.sqrt() + dropped).powi(2);
    }
}

/// Resolution control from the certificate's residual budget.
#[derive(Clone, Debug)]
pub struct AdaptPolicy {
    /// allowed grid sizes, ascending (FFT-friendly, even)
    pub sizes: Vec<usize>,
    /// target certified L² error at `t_end`
    pub target: f64,
    pub t_end: f64,
    /// coarsen only after this many consecutive steps far below budget
    pub patience: usize,
}

impl AdaptPolicy {
    pub fn new(target: f64, t_end: f64) -> Self {
        Self { sizes: vec![16, 24, 32, 48, 64, 96, 128, 192, 256], target, t_end, patience: 20 }
    }
}

/// Adaptive driver state.
pub struct Adaptive {
    pub policy: AdaptPolicy,
    quiet_steps: usize,
    /// (time, new N) for every change
    pub changes: Vec<(f64, usize)>,
}

impl Adaptive {
    pub fn new(policy: AdaptPolicy) -> Self {
        Self { policy, quiet_steps: 0, changes: Vec::new() }
    }

    /// Decide after a certified step of size `h`.
    ///
    /// Residual added at time t is amplified by exp(2∫_t^T μ) before t_end,
    /// so the residual budget is discounted by that factor, predicted with
    /// the current μ: β_eff = (target²/t_end)·exp(−2μ(T − t)). The bound's
    /// growth rate from the residual, ‖R‖²/(2ν), is compared with β_eff/2
    /// (headroom); far over budget, the grid jumps two sizes. Returns a new
    /// grid size if the solver should be regridded.
    pub fn decide(&mut self, s: &Solver, cert: &Certifier, st: &CertStep, h: f64) -> Option<usize> {
        let p = &self.policy;
        let remaining_t = (p.t_end - s.t).max(0.0);
        let beta = p.target * p.target / p.t_end * (-2.0 * st.mu * remaining_t).exp();
        let rate = (st.resid_trunc.powi(2) + st.resid_time.powi(2)) / (2.0 * s.nu);
        let n = s.grid.n;
        let pos = p.sizes.iter().position(|&m| m == n)?;
        if rate > 0.5 * beta && pos + 1 < p.sizes.len() {
            self.quiet_steps = 0;
            let jump = if rate > 16.0 * beta { 2 } else { 1 };
            return Some(p.sizes[(pos + jump).min(p.sizes.len() - 1)]);
        }
        if rate < beta / 100.0 {
            self.quiet_steps += 1;
        } else {
            self.quiet_steps = 0;
        }
        if self.quiet_steps >= p.patience && pos > 0 {
            // coarsen only if the dropped modes fit comfortably in the budget
            let (_, dropped) = s.regrid(p.sizes[pos - 1]);
            let remaining = (p.target - cert.y.sqrt()).max(0.0);
            self.quiet_steps = 0;
            let _ = h;
            if dropped < 0.05 * remaining {
                return Some(p.sizes[pos - 1]);
            }
        }
        None
    }
}

// ---------------------------------------------------------------------------
// Rigorous certificate: exact time quadrature, rigorous sup bounds, rounding
// allowance, and an H¹ a posteriori *regularity* certificate.
//
// Exact quadrature: in each step v(t) is a cubic Hermite polynomial in t, so
// the residual R(t) = v̇ + P(v·∇v) − νΔv (with its truncation part) is a
// polynomial of degree 6 in t and ‖R(t)‖² one of degree 12; 7-point
// Gauss–Legendre integrates it exactly.
//
// Sup bounds: sampled maxima of −λ_min(S) and |∇v|_F over the padded grid
// at 9 times per step, plus Lipschitz corrections (Weyl: eigenvalues move by
// at most the Frobenius change) from coefficient-space derivative bounds.
//
// Rounding: residual norms carry an allowance δ = 10·ε·log₂(M)·(‖v̇‖ +
// ‖F(v)‖ + ‖v×ω‖), a standard FFT error model (Higham) — not interval
// arithmetic.
//
// Regularity (Chernyshenko, Constantin, Robinson & Titi 2007): with
// w = u − v, test the error equation with −Δw. The terms are bounded by
//   B(w,w): Agmon (C_A from `constants`) + Young → ν/4‖Δw‖² + (27/4)C_A⁴y³/ν³
//   B(w,v): integrate by parts, −∫(∂w·∇)v·∂w − ∫(w·∇)∂v·∂w ≤ (G + H)y
//           with Poincaré ‖w‖ ≤ ‖∇w‖ (a G²/ν rate from Young here would be
//           useless: ~800 per unit time at Re = 200)
//   B(v,w): ∫(v·∇)w·Δw = −∫∂v·∂w·∂w ≤ G y   (div v = 0)
//   R:      ≤ ν/4‖Δw‖² + ‖R‖²/ν
// so  y' ≤ 2(2G + H) y + (27/2) C_A⁴ y³/ν³ + 2‖R‖²/ν,  y = ‖∇w‖²,
// G = sup |∇v|, H = sup |∇²v|. While y stays finite the true solution is strong (no
// blowup) and ‖∇(u − v)‖ ≤ √y; a bootstrap per step either extends that
// guarantee or reports that blowup can no longer be excluded.

const GL7: [(f64, f64); 7] = [
    (-0.949_107_912_342_758_5, 0.129_484_966_168_869_7),
    (-0.741_531_185_599_394_4, 0.279_705_391_489_276_7),
    (-0.405_845_151_377_397_2, 0.381_830_050_505_118_9),
    (0.0, 0.417_959_183_673_469_4),
    (0.405_845_151_377_397_2, 0.381_830_050_505_118_9),
    (0.741_531_185_599_394_4, 0.279_705_391_489_276_7),
    (0.949_107_912_342_758_5, 0.129_484_966_168_869_7),
];

/// Padded-grid evaluation of one state.
#[derive(Clone, Copy, Debug, Default)]
struct PadEval {
    trunc_hm1: f64,
    trunc_l2: f64,
    neg_lmin: f64,
    grad_f: f64,
    prod_l2: f64,
}

impl Padded {
    fn evaluate_full(&mut self, uh: &[Spec; 3]) -> PadEval {
        let gp = self.gp.clone();
        let nr = gp.real_len();
        // strain components (6) and vorticity (3) on the padded grid
        let comps = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];
        let mut sfield: Vec<Vec<f64>> = Vec::with_capacity(6);
        for &(i, j) in &comps {
            self.tmp.iter_mut().for_each(|v| *v = Complex64::default());
            for &(pi, ni) in &self.pairs {
                let k = gp.k(pi);
                self.tmp[pi] = 0.5 * I * (uh[j][ni] * k[i] + uh[i][ni] * k[j]);
            }
            let mut out = vec![0.0; nr];
            self.fft.inverse_inplace(&mut self.tmp, &mut out);
            sfield.push(out);
        }
        // velocity and vorticity, then the exact product
        for c in 0..3 {
            self.spec[c].iter_mut().for_each(|v| *v = Complex64::default());
            for &(pi, ni) in &self.pairs {
                self.spec[c][pi] = uh[c][ni];
            }
            self.fft.inverse(&self.spec[c], &mut self.tmp, &mut self.phys_u[c]);
        }
        for c in 0..3 {
            let (a, b) = ((c + 1) % 3, (c + 2) % 3);
            let (sa, sb) = (&self.spec[a], &self.spec[b]);
            self.tmp.par_iter_mut().enumerate().for_each(|(idx, w)| {
                let k = gp.k(idx);
                *w = I * (sb[idx] * k[a] - sa[idx] * k[b]);
            });
            self.fft.inverse_inplace(&mut self.tmp, &mut self.phys_w[c]);
        }
        let (neg_lmin, grad_f) = (0..nr)
            .into_par_iter()
            .map(|q| {
                let s = [sfield[0][q], sfield[1][q], sfield[2][q], sfield[3][q], sfield[4][q], sfield[5][q]];
                let nl = -sym3_min_eig(s[0], s[1], s[2], s[3], s[4], s[5]);
                let s2 = s[0] * s[0] + s[1] * s[1] + s[2] * s[2] + 2.0 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5]);
                let w2 = self.phys_w[0][q].powi(2) + self.phys_w[1][q].powi(2) + self.phys_w[2][q].powi(2);
                (nl, (s2 + 0.5 * w2).sqrt())
            })
            .reduce(|| (f64::NEG_INFINITY, 0.0), |a, b| (a.0.max(b.0), a.1.max(b.1)));
        let mut prod: [Vec<f64>; 3] = std::array::from_fn(|_| vec![0.0; nr]);
        {
            let [u0, u1, u2] = &self.phys_u;
            let [w0, w1, w2] = &self.phys_w;
            let [p0, p1, p2] = &mut prod;
            p0.par_iter_mut().zip(p1.par_iter_mut()).zip(p2.par_iter_mut()).enumerate().for_each(|(q, ((a, b), c))| {
                *a = u1[q] * w2[q] - u2[q] * w1[q];
                *b = u2[q] * w0[q] - u0[q] * w2[q];
                *c = u0[q] * w1[q] - u1[q] * w0[q];
            });
        }
        for c in 0..3 {
            self.fft.forward(&prod[c], &mut self.spec[c]);
        }
        let from_n = &self.from_n;
        let [a, b, c] = &self.spec;
        let (hm1, l2, pl2) = (0..gp.spec_len())
            .into_par_iter()
            .map(|idx| {
                let k = gp.k(idx);
                let k2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
                let w = gp.weight(idx);
                let all = w * (a[idx].norm_sqr() + b[idx].norm_sqr() + c[idx].norm_sqr());
                if from_n[idx].is_some() || k2 == 0.0 {
                    return (0.0, 0.0, all);
                }
                let d = (a[idx] * k[0] + b[idx] * k[1] + c[idx] * k[2]) / k2;
                let p = [a[idx] - d * k[0], b[idx] - d * k[1], c[idx] - d * k[2]];
                let m2 = w * (p[0].norm_sqr() + p[1].norm_sqr() + p[2].norm_sqr());
                (m2 / k2, m2, all)
            })
            .reduce(|| (0.0, 0.0, 0.0), |x, y| (x.0 + y.0, x.1 + y.1, x.2 + y.2));
        PadEval { trunc_hm1: hm1, trunc_l2: l2, neg_lmin: neg_lmin.max(0.0), grad_f, prod_l2: pl2.sqrt() }
    }
}

/// Fine sampling grid (twice the padded grid) for the maxima of
/// −λ_min(S) and |∇v|_F, which halves the spatial Lipschitz correction.
struct Sampler {
    gs: Grid,
    fft: Fft3,
    pairs: Vec<(usize, usize)>,
    tmp: Spec,
}

impl Sampler {
    fn new(gn: &Grid, m_padded: usize) -> Self {
        let m = 2 * m_padded;
        let gs = Grid::new(m);
        let km = gn.kmax() as i64;
        let mut pairs = Vec::new();
        for idx in 0..gn.spec_len() {
            if !gn.kept(idx) {
                continue;
            }
            let (i, j, l) = gn.spec_ijl(idx);
            let (kx, ky) = (wavenumber(i, gn.n), wavenumber(j, gn.n));
            if kx.abs() <= km && ky.abs() <= km && (l as i64) <= km {
                let w = |k: i64| k.rem_euclid(m as i64) as usize;
                pairs.push(((w(kx) * m + w(ky)) * gs.nh + l, idx));
            }
        }
        Self { fft: Fft3::new(&gs), tmp: vec![Complex64::default(); gs.spec_len()], pairs, gs }
    }

    /// (max −λ_min(S), max |∇v|_F) over the fine grid.
    fn maxima(&mut self, uh: &[Spec; 3]) -> (f64, f64) {
        let gs = self.gs.clone();
        let nr = gs.real_len();
        let mut fields: Vec<Vec<f64>> = Vec::with_capacity(9);
        // S_ij (6) and ω (3)
        let comps = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];
        for &(i, j) in &comps {
            self.tmp.iter_mut().for_each(|v| *v = Complex64::default());
            for &(pi, ni) in &self.pairs {
                let k = gs.k(pi);
                self.tmp[pi] = 0.5 * I * (uh[j][ni] * k[i] + uh[i][ni] * k[j]);
            }
            let mut out = vec![0.0; nr];
            self.fft.inverse_inplace(&mut self.tmp, &mut out);
            fields.push(out);
        }
        for c in 0..3 {
            let (a, b) = ((c + 1) % 3, (c + 2) % 3);
            self.tmp.iter_mut().for_each(|v| *v = Complex64::default());
            for &(pi, ni) in &self.pairs {
                let k = gs.k(pi);
                self.tmp[pi] = I * (uh[b][ni] * k[a] - uh[a][ni] * k[b]);
            }
            let mut out = vec![0.0; nr];
            self.fft.inverse_inplace(&mut self.tmp, &mut out);
            fields.push(out);
        }
        (0..nr)
            .into_par_iter()
            .map(|q| {
                let f = |n: usize| fields[n][q];
                let nl = -sym3_min_eig(f(0), f(1), f(2), f(3), f(4), f(5));
                let s2 = f(0) * f(0) + f(1) * f(1) + f(2) * f(2) + 2.0 * (f(3) * f(3) + f(4) * f(4) + f(5) * f(5));
                let w2 = f(6) * f(6) + f(7) * f(7) + f(8) * f(8);
                (nl, (s2 + 0.5 * w2).sqrt())
            })
            .reduce(|| (0.0, 0.0), |a, b| (a.0.max(b.0), a.1.max(b.1)))
    }
}

/// One rigorously certified step.
#[derive(Clone, Copy, Debug, Default)]
pub struct RigorousStep {
    pub t: f64,
    /// rigorous bound on ‖u − v‖ (box-mean L²), any Leray–Hopf solution u
    pub l2_bound: f64,
    /// rigorous bound on ‖∇(u − v)‖; `None` once blowup can no longer be
    /// excluded (the regularity certificate has lapsed)
    pub h1_bound: Option<f64>,
    /// rigorous sup over the step of −λ_min(sym ∇v) and of |∇v|
    pub mu_sup: f64,
    pub grad_sup: f64,
    /// sampled (non-rigorous) counterparts, to show the cost of rigor
    pub mu_sampled: f64,
    pub grad_sampled: f64,
    /// rigorous sup |∇²v| over the step
    pub hess_sup: f64,
}

pub struct RigorousCertifier {
    padded: Padded,
    sampler: Sampler,
    prev_u: Option<[Spec; 3]>,
    prev_f: Option<[Spec; 3]>,
    prev_eval: PadEval,
    y_l2: f64,
    y_h1: Option<f64>,
    c_a: f64,
    pub history: Vec<RigorousStep>,
}

impl RigorousCertifier {
    pub fn new(grid: &Grid, initial_l2: f64, initial_h1: f64) -> Self {
        let padded = Padded::new(grid);
        let sampler = Sampler::new(grid, padded.gp.n);
        Self {
            padded,
            sampler,
            prev_u: None,
            prev_f: None,
            prev_eval: PadEval::default(),
            y_l2: initial_l2 * initial_l2,
            y_h1: Some(initial_h1 * initial_h1),
            c_a: crate::constants::agmon().2,
            history: Vec::new(),
        }
    }
}

/// Hermite basis values at τ ∈ [0, 1] (h00, h10, h01, h11) and derivatives.
fn hermite(tau: f64) -> ([f64; 4], [f64; 4]) {
    let (t2, t3) = (tau * tau, tau * tau * tau);
    (
        [2.0 * t3 - 3.0 * t2 + 1.0, t3 - 2.0 * t2 + tau, -2.0 * t3 + 3.0 * t2, t3 - t2],
        [6.0 * t2 - 6.0 * tau, 3.0 * t2 - 4.0 * tau + 1.0, -6.0 * t2 + 6.0 * tau, 3.0 * t2 - 2.0 * tau],
    )
}

impl Solver {
    /// Weighted H⁻¹ and L² norms² of a (projected) spectral field.
    fn hm1_l2(&self, r: &[Spec; 3]) -> (f64, f64) {
        let g = &self.grid;
        let k2 = &self.k2;
        (0..g.spec_len())
            .into_par_iter()
            .map(|i| {
                if k2[i] == 0.0 {
                    return (0.0, 0.0);
                }
                let m2 = g.weight(i) * (r[0][i].norm_sqr() + r[1][i].norm_sqr() + r[2][i].norm_sqr());
                (m2 / k2[i], m2)
            })
            .reduce(|| (0.0, 0.0), |a, b| (a.0 + b.0, a.1 + b.1))
    }

    fn l2_norm(&self, r: &[Spec; 3]) -> f64 {
        let g = &self.grid;
        (0..g.spec_len())
            .into_par_iter()
            .map(|i| g.weight(i) * (r[0][i].norm_sqr() + r[1][i].norm_sqr() + r[2][i].norm_sqr()))
            .sum::<f64>()
            .sqrt()
    }

    /// Advance one IF-RK4 step and extend the rigorous certificates.
    pub fn step_certified_rigorous(&mut self, h: f64, cert: &mut RigorousCertifier) -> RigorousStep {
        assert!(self.forcing.is_none() && self.penalization.is_none() && self.body_force.is_none() && self.surface.is_none());
        assert!(!self.grid.spherical, "the certificate assumes the exact 2/3-rule discretization");
        let nu = self.nu;
        let eps_r = 10.0 * f64::EPSILON * (cert.padded.gp.n as f64).log2();
        if cert.prev_u.is_none() {
            let u0 = self.u.clone();
            cert.prev_eval = cert.padded.evaluate_full(&u0);
            let (a, b) = cert.sampler.maxima(&u0);
            cert.prev_eval.neg_lmin = a;
            cert.prev_eval.grad_f = b;
            cert.prev_f = Some(self.tendency(&u0));
            cert.prev_u = Some(u0);
        }
        let (u0, f0) = (cert.prev_u.take().unwrap(), cert.prev_f.take().unwrap());
        self.step_ifrk4(h);
        let u1 = self.u.clone();
        let f1 = self.tendency(&u1);
        let ns = self.grid.spec_len();

        // residual at the 7 Gauss nodes, exact quadrature of ‖R‖² (H⁻¹, L²)
        let (mut int_hm1, mut int_l2) = (0.0, 0.0);
        let (mut mu_s, mut grad_s) = (cert.prev_eval.neg_lmin, cert.prev_eval.grad_f);
        for &(xg, wg) in &GL7 {
            let tau = 0.5 * (xg + 1.0);
            let (b, db) = hermite(tau);
            let mut v = zeros_spec(ns);
            let mut vdot = zeros_spec(ns);
            for c in 0..3 {
                v[c].par_iter_mut().zip(vdot[c].par_iter_mut()).enumerate().for_each(|(i, (x, d))| {
                    *x = u0[c][i] * b[0] + f0[c][i] * (h * b[1]) + u1[c][i] * b[2] + f1[c][i] * (h * b[3]);
                    *d = u0[c][i] * (db[0] / h) + f0[c][i] * db[1] + u1[c][i] * (db[2] / h) + f1[c][i] * db[3];
                });
            }
            let fv = self.tendency(&v);
            let mut rt = zeros_spec(ns);
            for c in 0..3 {
                rt[c].par_iter_mut().enumerate().for_each(|(i, r)| *r = vdot[c][i] - fv[c][i]);
            }
            project(&self.grid, &mut rt);
            let (t_hm1, t_l2) = self.hm1_l2(&rt);
            let pe = cert.padded.evaluate_full(&v);
            let (sm, sg) = cert.sampler.maxima(&v);
            mu_s = mu_s.max(sm);
            grad_s = grad_s.max(sg);
            // rigorous node values: computed + rounding allowance (‖·‖_{H⁻¹} ≤ ‖·‖_{L²})
            let delta = eps_r * (self.l2_norm(&vdot) + self.l2_norm(&fv) + pe.prod_l2);
            let r_l2 = (t_l2 + pe.trunc_l2).sqrt() + delta;
            let r_hm1 = (t_hm1 + pe.trunc_hm1).sqrt() + delta;
            int_hm1 += 0.5 * h * wg * r_hm1 * r_hm1;
            int_l2 += 0.5 * h * wg * r_l2 * r_l2;
        }
        let mut pe1 = cert.padded.evaluate_full(&u1);
        let (sm1, sg1) = cert.sampler.maxima(&u1);
        pe1.neg_lmin = sm1;
        pe1.grad_f = sg1;
        mu_s = mu_s.max(sm1);
        grad_s = grad_s.max(sg1);

        // Lipschitz corrections from coefficient bounds over the step
        let g = &self.grid;
        let (lx, lt) = (0..ns)
            .into_par_iter()
            .map(|i| {
                if !g.kept(i) {
                    return (0.0, 0.0);
                }
                let k = g.k(i);
                let kn = (k[0] * k[0] + k[1] * k[1] + k[2] * k[2]).sqrt();
                let w = g.weight(i);
                let mag = |s: &[Spec; 3]| (s[0][i].norm_sqr() + s[1][i].norm_sqr() + s[2][i].norm_sqr()).sqrt();
                // v(τ) = v0 + h01(τ)Δ + h(h10 F0 + h11 F1), Δ = v1 − v0 (keeps the
                // cancellation v1 ≈ v0; bounding |v0| + |v1| separately put a
                // 1/h-sized term in the time derivative)
                let delta = (0..3).map(|c| (u1[c][i] - u0[c][i]).norm_sqr()).sum::<f64>().sqrt();
                let (a0, b0, b1) = (mag(&u0), mag(&f0), mag(&f1));
                let m = a0 + delta + (4.0 / 27.0) * h * (b0 + b1); // sup_t |v̂_k(t)|
                let d = 1.5 * delta / h + b0 + b1; // sup_t |∂_t v̂_k(t)|
                (w * kn * kn * m, w * kn * d)
            })
            .reduce(|| (0.0, 0.0), |a, b| (a.0 + b.0, a.1 + b.1));
        let hs = cert.sampler.gs.dx();
        let mut times: Vec<f64> = vec![0.0, 1.0];
        times.extend(GL7.iter().map(|&(x, _)| 0.5 * (x + 1.0)));
        times.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let gap = times.windows(2).map(|w| w[1] - w[0]).fold(0.0, f64::max);
        // ∇v: |k v̂ᵀ|_F = |k||v̂|; strain of a solenoidal field: |Ŝ_k|_F = |k||v̂_k|/√2
        let corr = lx * 3f64.sqrt() / 2.0 * hs + lt * 0.5 * gap * h;
        let (mu_sup, grad_sup) = (mu_s + corr / 2f64.sqrt(), grad_s + corr);
        // sup |∇²v|_F over the step: Σ_k |k|² sup_t |v̂_k(t)| (= lx)
        let hess_sup = lx;

        // L² weak–strong bound
        cert.y_l2 = (2.0 * mu_sup * h).exp() * (cert.y_l2 + int_hm1 / (2.0 * nu));
        // H¹ regularity bootstrap: find Ȳ with e^{(A + BȲ²)h}(y + ∫r) ≤ Ȳ
        cert.y_h1 = cert.y_h1.and_then(|y| {
            let a_c = 2.0 * (2.0 * grad_sup + hess_sup);
            let b_c = 13.5 * cert.c_a.powi(4) / nu.powi(3);
            let base = y + 2.0 * int_l2 / nu;
            let mut ybar = 2.0 * base * (a_c * h).exp() + 1e-300;
            for _ in 0..60 {
                let val = ((a_c + b_c * ybar * ybar) * h).exp() * base;
                if val <= ybar {
                    return Some(val);
                }
                ybar *= 2.0;
                if !ybar.is_finite() || b_c * ybar * ybar * h > 50.0 {
                    return None;
                }
            }
            None
        });
        cert.prev_u = Some(u1);
        cert.prev_f = Some(f1);
        cert.prev_eval = pe1;
        let st = RigorousStep {
            t: self.t,
            l2_bound: cert.y_l2.sqrt(),
            h1_bound: cert.y_h1.map(f64::sqrt),
            mu_sup,
            grad_sup,
            mu_sampled: mu_s,
            grad_sampled: grad_s,
            hess_sup,
        };
        cert.history.push(st);
        st
    }
}
