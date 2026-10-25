//! IMEX time stepping with implicit viscosity and implicit Brinkman
//! penalization, for flows with solid objects.
//!
//! Splitting du/dt = N(u) + L(u) with
//!   N(u) = P[u × ω + sponge + f]                    (explicit)
//!   L(u) = P[ν Δu − (χ/η)(u − u_s)]                 (implicit)
//! and the ARS(3,4,3) scheme (Ascher, Ruuth & Spiteri 1997), third order,
//! L-stable and stiffly accurate. Its steady states are exactly the steady
//! states of the penalized equations, independent of dt, and η can be made
//! tiny (true no-slip) without shrinking dt.
//!
//! Each implicit stage solves (I − aL)U = r, a = hγ. With H = I + aν|k|²
//! (diagonal in Fourier) and the penalty force g = (a/η) χ (U − u_s), which
//! lives only on the body:
//!   Û = H⁻¹(r̂ − P ĝ),     (D + S) g = χ>0 · (F⁻¹H⁻¹r̂ − u_s),
//!   D = η / (a χ),          S = F⁻¹ H⁻¹ P F   (symmetric, 0 ≤ S ≤ 1).
//! D + S is symmetric positive definite on the body's support, so it is
//! solved with conjugate gradients, preconditioned by the exact Fourier
//! inverse of the full-space operator restricted to the body, and
//! warm-started from the previous solve of the same stage.
//!
//! Where the constraint lives matters. The velocity is a trigonometric
//! polynomial with only the 2/3-retained modes, |k_i| ≤ k_max; imposing the
//! penalization at all N³ grid points over-determines it (a band-limited
//! field cannot vanish at every solid grid point), which crushes the flow as
//! η → 0 (plane Poiseuille: −27% at N = 64). When the body's geometry is
//! known, the system can be solved on the odd "solution grid"
//! M = 2k_max + 1, whose Fourier modes are exactly the retained modes: one
//! constraint per degree of freedom, and ~3× fewer points for the CG
//! (`Penalization::on_solution_grid`). Measured, it is under-constrained
//! between collocation points and converges from above, no better than the
//! N grid with grid-matched η, so it is opt-in.

use rayon::prelude::*;

use crate::fft::Fft3;
use crate::nufft::Nufft;
use crate::grid::{wavenumber, Grid};
use crate::solver::{project_with, zeros_spec, Solver, Spec};

const GAMMA: f64 = 0.435_866_521_508_459;

/// ARS(3,4,3) tableaux. Stage 1 is explicit (U1 = uⁿ). Shared with the
/// CUDA backend.
pub struct Ars343 {
    pub ae: [[f64; 4]; 4],
    pub ai: [[f64; 4]; 4],
    pub b: [f64; 4],
}

pub fn ars343() -> Ars343 {
    let g = GAMMA;
    let b1 = -1.5 * g * g + 4.0 * g - 0.25;
    let b2 = 1.5 * g * g - 5.0 * g + 1.25;
    Ars343 {
        ae: [
            [0.0; 4],
            [g, 0.0, 0.0, 0.0],
            [0.321_278_886_0, 0.396_654_374_7, 0.0, 0.0],
            [-0.105_858_296, 0.552_929_147_9, 0.552_929_147_9, 0.0],
        ],
        ai: [
            [0.0; 4],
            [0.0, g, 0.0, 0.0],
            [0.0, 0.5 * (1.0 - g), g, 0.0],
            [0.0, b1, b2, g],
        ],
        b: [0.0, b1, b2, g],
    }
}

/// The grid on which the penalty-force system is posed and solved.
pub(crate) struct SolveSpace {
    grid: Grid,
    fft: Fft3,
    k2: Vec<f64>,
    /// apply the 2/3 mask (only when solving on the N grid itself)
    dealias: bool,
    /// spectral index in this space → spectral index in the N-grid spectrum
    /// (`None`: the space is the N grid)
    to_n: Option<Vec<usize>>,
    chi: Vec<f64>,
    tmp: Spec,
    r: [Vec<f64>; 3],
    z: [Vec<f64>; 3],
    p: [Vec<f64>; 3],
    ap: [Vec<f64>; 3],
    b: [Vec<f64>; 3],
    spec: [Spec; 3],
    spec2: [Spec; 3],
    /// warm starts for the penalty force, per implicit stage
    g_prev: [[Vec<f64>; 3]; 4],
}

impl SolveSpace {
    fn new(grid: Grid, dealias: bool, to_n: Option<Vec<usize>>, chi: Vec<f64>) -> Self {
        let (ns, nr) = (grid.spec_len(), grid.real_len());
        let k2 = (0..ns)
            .map(|i| {
                let k = grid.k(i);
                k[0] * k[0] + k[1] * k[1] + k[2] * k[2]
            })
            .collect();
        Self {
            fft: Fft3::new(&grid),
            k2,
            dealias,
            to_n,
            chi,
            tmp: vec![Default::default(); ns],
            r: zeros_phys(nr),
            z: zeros_phys(nr),
            p: zeros_phys(nr),
            ap: zeros_phys(nr),
            b: zeros_phys(nr),
            spec: zeros_spec(ns),
            spec2: zeros_spec(ns),
            g_prev: std::array::from_fn(|_| zeros_phys(nr)),
            grid,
        }
    }

    /// Solution grid M = 2 k_max + 1 for an N grid, sampling `shape` sharply.
    fn solution_grid(gn: &Grid, shape: &crate::penalization::Shape) -> Self {
        let m = 2 * gn.kmax() + 1;
        let gm = Grid::new_any(m);
        let to_n = (0..gm.spec_len())
            .map(|idx| {
                let (i, j, l) = gm.spec_ijl(idx);
                let wrap = |k: i64| k.rem_euclid(gn.n as i64) as usize;
                let (kx, ky) = (wavenumber(i, m), wavenumber(j, m));
                (wrap(kx) * gn.n + wrap(ky)) * gn.nh + l
            })
            .collect();
        Self::new(gm, false, Some(to_n), shape.sharp_mask(m))
    }

}

/// Leray projection on a solve space (k = 0 kept when `keep_mean`); the 2/3
/// mask is applied only when the space is the N grid itself.
fn leray(g: &Grid, k2: &[f64], dealias: bool, v: &mut [Spec; 3], keep_mean: bool) {
    {
        if dealias {
            return project_with(g, v, keep_mean);
        }
        let [a, b, c] = v;
        a.par_iter_mut()
            .zip(b.par_iter_mut())
            .zip(c.par_iter_mut())
            .enumerate()
            .for_each(|(idx, ((x, y), z))| {
                if k2[idx] == 0.0 {
                    if !keep_mean {
                        *x = Default::default();
                        *y = Default::default();
                        *z = Default::default();
                    }
                    return;
                }
                let k = g.k(idx);
                let d = (*x * k[0] + *y * k[1] + *z * k[2]) / k2[idx];
                *x -= d * k[0];
                *y -= d * k[1];
                *z -= d * k[2];
            });
    }
}

/// Workspace for the surface-constraint solve (marker forces F, 3 per marker).
struct SurfaceWork {
    nufft: Nufft,
    /// smooth forcing extension operator and its per-stage warm starts
    ext: Option<crate::surface::ExtensionOp>,
    a_prev: [Vec<f64>; 4],
    phys: [Vec<f64>; 3],
    f_prev: [[Vec<f64>; 3]; 4],
    r: [Vec<f64>; 3],
    p: [Vec<f64>; 3],
    ap: [Vec<f64>; 3],
    b: [Vec<f64>; 3],
    spec: [Spec; 3],
}

impl SurfaceWork {
    fn new(g: &Grid, nm: usize) -> Self {
        let z = || std::array::from_fn(|_| vec![0.0; nm]);
        Self {
            nufft: Nufft::new(g, 12),
            ext: None,
            a_prev: std::array::from_fn(|_| vec![0.0; 3 * nm]),
            phys: std::array::from_fn(|_| vec![0.0; g.real_len()]),
            f_prev: std::array::from_fn(|_| z()),
            r: z(),
            p: z(),
            ap: z(),
            b: z(),
            spec: zeros_spec(g.spec_len()),
        }
    }
}

fn dot_m(a: &[Vec<f64>; 3], b: &[Vec<f64>; 3]) -> f64 {
    (0..3).map(|c| a[c].iter().zip(&b[c]).map(|(x, y)| x * y).sum::<f64>()).sum()
}

/// Persistent IMEX workspace (allocated on first use).
pub struct ImexWork {
    nl: [[Spec; 3]; 4],
    lin: [[Spec; 3]; 4],
    stage: [Spec; 3],
    rhs: [Spec; 3],
    space: Option<SolveSpace>,
    surf: Option<SurfaceWork>,
    pub last_iters: [usize; 4],
    /// final relative CG residual per stage
    pub last_resid: [f64; 4],
    pub last_force: [f64; 3],
}

fn zeros_phys(n: usize) -> [Vec<f64>; 3] {
    std::array::from_fn(|_| vec![0.0; n])
}

impl ImexWork {
    fn new(ns: usize) -> Self {
        Self {
            nl: std::array::from_fn(|_| zeros_spec(ns)),
            lin: std::array::from_fn(|_| zeros_spec(ns)),
            stage: zeros_spec(ns),
            rhs: zeros_spec(ns),
            space: None,
            surf: None,
            last_iters: [0; 4],
            last_resid: [0.0; 4],
            last_force: [0.0; 3],
        }
    }
}

fn dot3(a: &[Vec<f64>; 3], b: &[Vec<f64>; 3]) -> f64 {
    (0..3).map(|c| a[c].par_iter().zip(&b[c]).map(|(x, y)| x * y).sum::<f64>()).sum()
}

impl Solver {
    /// Relative residual tolerance of the penalty-force CG solve.
    pub const CG_TOL: f64 = 1e-10;
    pub const CG_MAX_ITERS: usize = 500;
    /// GMRES tolerance / budget for the smooth forcing extension solve.
    pub const EXT_TOL: f64 = 1e-8;
    pub const EXT_MAX_MATVECS: usize = 400;

    /// One ARS(3,4,3) IMEX step of size h (used when the solid is implicit).
    pub fn step_imex(&mut self, h: f64) {
        let tab = ars343();
        let ns = self.grid.spec_len();
        let mut w = self.imex.take().unwrap_or_else(|| Box::new(ImexWork::new(ns)));
        let u0 = std::mem::take(&mut self.u);

        // Stage 1 (explicit): U1 = uⁿ
        self.rhs(&u0, &mut w.nl[0]);
        for i in 1..4 {
            // r = uⁿ + h Σ_j ae[i][j] N_j + h Σ_{1≤j<i} ai[i][j] L_j
            for c in 0..3 {
                let rc = &mut w.rhs[c];
                rc.par_iter_mut().enumerate().for_each(|(idx, v)| {
                    let mut acc = u0[c][idx];
                    for j in 0..i {
                        acc += w.nl[j][c][idx] * (h * tab.ae[i][j]);
                        if j >= 1 {
                            acc += w.lin[j][c][idx] * (h * tab.ai[i][j]);
                        }
                    }
                    *v = acc;
                });
            }
            let a = h * tab.ai[i][i];
            let rhs = std::mem::take(&mut w.rhs);
            let mut stage = std::mem::take(&mut w.stage);
            let (iters, resid, force) = self.implicit_solve(&rhs, a, &mut stage, &mut w, i);
            w.last_iters[i] = iters;
            w.last_resid[i] = resid;
            if i == 3 {
                w.last_force = force;
            }
            // L(U_i) = (U_i − r)/a, exact algebra (no division by η)
            for c in 0..3 {
                w.lin[i][c]
                    .par_iter_mut()
                    .zip(stage[c].par_iter().zip(rhs[c].par_iter()))
                    .for_each(|(l, (u, r))| *l = (*u - *r) / a);
            }
            self.rhs(&stage, &mut w.nl[i]);
            w.rhs = rhs;
            w.stage = stage;
        }
        // uⁿ⁺¹ = uⁿ + h Σ_j b_j (N_j + L_j)   (b for both tableaux)
        let mut u1 = u0;
        for c in 0..3 {
            u1[c].par_iter_mut().enumerate().for_each(|(idx, v)| {
                let mut acc = *v;
                for j in 0..4 {
                    if tab.b[j] != 0.0 {
                        acc += (w.nl[j][c][idx] + w.lin[j][c][idx]) * (h * tab.b[j]);
                    }
                }
                *v = acc;
            });
        }
        self.u = u1;
        self.imex = Some(w);
        self.t += h;
    }

    /// Solve (I − aL)U = r. Returns (CG iterations, residual, force on the solid).
    fn implicit_solve(
        &mut self,
        r: &[Spec; 3],
        a: f64,
        out: &mut [Spec; 3],
        w: &mut ImexWork,
        stage: usize,
    ) -> (usize, f64, [f64; 3]) {
        let nu = self.nu;
        let k2n = &self.k2;
        // Û = H⁻¹ r̂ (k = 0 untouched: H(0) = 1)
        for c in 0..3 {
            out[c]
                .par_iter_mut()
                .zip(r[c].par_iter().zip(k2n.par_iter()))
                .for_each(|(o, (x, &kk))| *o = *x / (1.0 + a * nu * kk));
        }
        if let Some(sf) = &self.surface {
            return if sf.extension.is_some() {
                self.extension_solve(a, out, w, stage)
            } else {
                self.surface_solve(a, out, w, stage)
            };
        }
        let Some(pen) = self.penalization.as_ref().filter(|p| p.implicit) else {
            return (0, 0.0, [0.0; 3]);
        };
        let (eta, us) = (pen.eta, pen.u_solid);
        if w.space.is_none() {
            w.space = Some(match (&pen.shape, pen.on_solution_grid) {
                (Some(shape), true) => SolveSpace::solution_grid(&self.grid, shape),
                _ => SolveSpace::new(self.grid.clone(), true, None, pen.chi.clone()),
            });
        }
        let keep_mean = self.body_force.is_some();
        let sp = w.space.as_mut().unwrap();
        let SolveSpace { grid: g, fft, k2, dealias, to_n, chi, tmp, r: cr, z, p, ap, b, spec, spec2, g_prev } = sp;
        let (g, fft, k2, chi, dealias, to_n) = (&*g, &*fft, &*k2, &*chi, *dealias, to_n.as_ref());
        let leray = |v: &mut [Spec; 3]| leray(g, k2, dealias, v, keep_mean);

        // b = χ>0 · (F⁻¹H⁻¹r̂ − u_s) on the solve space
        for c in 0..3 {
            match to_n {
                None => spec[c].copy_from_slice(&out[c]),
                Some(map) => spec[c].par_iter_mut().zip(map.par_iter()).for_each(|(s, &i)| *s = out[c][i]),
            }
            fft.inverse_inplace(&mut spec[c], &mut b[c]);
            b[c].par_iter_mut().zip(chi.par_iter()).for_each(|(bv, &x)| {
                *bv = if x > 0.0 { *bv - us[c] } else { 0.0 };
            });
        }
        let dvec = |x: f64| if x > 0.0 { eta / (a * x) } else { 0.0 };

        // A x = D x + χ>0 · S x
        let apply = |x: &[Vec<f64>; 3], ax: &mut [Vec<f64>; 3], spec: &mut [Spec; 3], tmp: &mut Spec| {
            for c in 0..3 {
                fft.forward(&x[c], &mut spec[c]);
            }
            leray(spec);
            for c in 0..3 {
                spec[c].par_iter_mut().zip(k2.par_iter()).for_each(|(v, &kk)| *v /= 1.0 + a * nu * kk);
                fft.inverse(&spec[c], tmp, &mut ax[c]);
                ax[c].par_iter_mut().enumerate().for_each(|(q, v)| {
                    let xq = chi[q];
                    *v = if xq > 0.0 { *v + dvec(xq) * x[c][q] } else { 0.0 };
                });
            }
        };
        // Preconditioner: restricted exact Fourier inverse of D̄ + S, D̄ = η/a:
        // H/(1 + D̄H) on solenoidal modes, 1/D̄ on gradient (and removed) modes.
        let dbar = eta / a;
        let precond = |rr: &[Vec<f64>; 3], zz: &mut [Vec<f64>; 3], spec: &mut [Spec; 3], sol: &mut [Spec; 3], tmp: &mut Spec| {
            for c in 0..3 {
                fft.forward(&rr[c], &mut spec[c]);
                sol[c].copy_from_slice(&spec[c]);
            }
            leray(sol);
            for c in 0..3 {
                spec[c].par_iter_mut().zip(sol[c].par_iter().zip(k2.par_iter())).for_each(|(v, (s, &kk))| {
                    let h = 1.0 + a * nu * kk;
                    *v = *s * (h / (1.0 + dbar * h)) + (*v - *s) / dbar;
                });
                fft.inverse(&spec[c], tmp, &mut zz[c]);
                zz[c].par_iter_mut().zip(chi.par_iter()).for_each(|(zv, &x)| {
                    if x <= 0.0 {
                        *zv = 0.0;
                    }
                });
            }
        };

        // PCG, warm-started
        let gvec = &mut g_prev[stage];
        apply(gvec, ap, spec, tmp);
        for c in 0..3 {
            cr[c].par_iter_mut().zip(b[c].par_iter().zip(ap[c].par_iter())).for_each(|(rv, (bv, av))| *rv = bv - av);
        }
        let bnorm = dot3(b, b).sqrt().max(1e-300);
        precond(cr, z, spec, spec2, tmp);
        for c in 0..3 {
            p[c].copy_from_slice(&z[c]);
        }
        let mut rz = dot3(cr, z);
        let mut iters = 0;
        while dot3(cr, cr).sqrt() > Self::CG_TOL * bnorm && iters < Self::CG_MAX_ITERS {
            apply(p, ap, spec, tmp);
            let alpha = rz / dot3(p, ap);
            for c in 0..3 {
                gvec[c].par_iter_mut().zip(p[c].par_iter()).for_each(|(x, pv)| *x += alpha * pv);
                cr[c].par_iter_mut().zip(ap[c].par_iter()).for_each(|(x, av)| *x -= alpha * av);
            }
            precond(cr, z, spec, spec2, tmp);
            let rz_new = dot3(cr, z);
            let beta = rz_new / rz;
            rz = rz_new;
            for c in 0..3 {
                p[c].par_iter_mut().zip(z[c].par_iter()).for_each(|(pv, zv)| *pv = zv + beta * *pv);
            }
            iters += 1;
        }
        let resid = dot3(cr, cr).sqrt() / bnorm;

        // Û = H⁻¹(r̂ − P ĝ): ĝ from the solve space, mapped to the N spectrum
        for c in 0..3 {
            fft.forward(&gvec[c], &mut spec[c]);
        }
        leray(spec);
        for c in 0..3 {
            match to_n {
                None => out[c]
                    .par_iter_mut()
                    .zip(spec[c].par_iter().zip(k2n.par_iter()))
                    .for_each(|(o, (gs, &kk))| *o -= *gs / (1.0 + a * nu * kk)),
                Some(map) => {
                    // distinct space indices map to distinct N indices
                    for (idx, &i) in map.iter().enumerate() {
                        out[c][i] -= spec[c][idx] / (1.0 + a * nu * k2n[i]);
                    }
                }
            }
        }
        // force density on the solid = g / a
        let dv = g.dx().powi(3);
        let force = [0, 1, 2].map(|c| gvec[c].par_iter().sum::<f64>() * dv / a);
        (iters, resid, force)
    }
}

impl Solver {
    /// Surface constraint in an implicit stage. On entry `out` = H⁻¹r̂; solves
    ///   a J H⁻¹ P S F = u_s − J H⁻¹ r̂
    /// for the marker forces F (CG; the operator is symmetric PSD because
    /// J = (2π)³ S*), then sets Û = H⁻¹ r̂ + a H⁻¹ P S F. The force on the
    /// body is −Σ_m F_m.
    fn surface_solve(&mut self, a: f64, out: &mut [Spec; 3], w: &mut ImexWork, stage: usize) -> (usize, f64, [f64; 3]) {
        let surface = self.surface.as_ref().unwrap();
        let pts = &surface.points;
        let nm = pts.len();
        if w.surf.is_none() {
            w.surf = Some(SurfaceWork::new(&self.grid, nm));
        }
        let sw = w.surf.as_mut().unwrap();
        let (g, k2, nu) = (&self.grid, &self.k2, self.nu);
        let keep_mean = self.body_force.is_some();
        let SurfaceWork { nufft, f_prev, r, p, ap, b, spec, .. } = sw;

        // A x = a J H⁻¹ P S x
        let apply = |x: &[Vec<f64>; 3], ax: &mut [Vec<f64>; 3], spec: &mut [Spec; 3], nufft: &mut Nufft| {
            for c in 0..3 {
                nufft.spread(&x[c], pts, &mut spec[c]);
            }
            project_with(g, spec, keep_mean);
            for c in 0..3 {
                spec[c].iter_mut().zip(k2.iter()).for_each(|(v, &kk)| *v *= a / (1.0 + a * nu * kk));
                ax[c] = nufft.interpolate(&spec[c], pts);
            }
        };
        for c in 0..3 {
            let v = nufft.interpolate(&out[c], pts);
            for (m, (bv, x)) in b[c].iter_mut().zip(v).enumerate() {
                *bv = surface.target(m, c) - x;
            }
        }
        let f = &mut f_prev[stage];
        apply(f, ap, spec, nufft);
        for c in 0..3 {
            for m in 0..nm {
                r[c][m] = b[c][m] - ap[c][m];
                p[c][m] = r[c][m];
            }
        }
        let bnorm = dot_m(b, b).sqrt().max(1e-300);
        let mut rr = dot_m(r, r);
        let mut iters = 0;
        while rr.sqrt() > Self::CG_TOL * bnorm && iters < Self::CG_MAX_ITERS {
            apply(p, ap, spec, nufft);
            let alpha = rr / dot_m(p, ap);
            for c in 0..3 {
                for m in 0..nm {
                    f[c][m] += alpha * p[c][m];
                    r[c][m] -= alpha * ap[c][m];
                }
            }
            let rr_new = dot_m(r, r);
            let beta = rr_new / rr;
            rr = rr_new;
            for c in 0..3 {
                for m in 0..nm {
                    p[c][m] = r[c][m] + beta * p[c][m];
                }
            }
            iters += 1;
        }
        // Û = H⁻¹r̂ + a H⁻¹ P S F
        for c in 0..3 {
            nufft.spread(&f[c], pts, &mut spec[c]);
        }
        project_with(g, spec, keep_mean);
        for c in 0..3 {
            out[c]
                .iter_mut()
                .zip(spec[c].iter().zip(k2.iter()))
                .for_each(|(o, (v, &kk))| *o += *v * (a / (1.0 + a * nu * kk)));
        }
        let force = [0, 1, 2].map(|c| -f[c].iter().sum::<f64>());
        (iters, rr.sqrt() / bnorm, force)
    }
}

impl Solver {
    /// Smooth forcing extension in an implicit stage (`out` = H⁻¹r̂ on entry):
    /// find extension amplitudes A (3 per marker) with J H⁻¹ P (a·E A) = u_s − J H⁻¹ r̂
    /// by GMRES, then Û = H⁻¹ r̂ + a H⁻¹ P E A. Force on the body: −∫ E A.
    fn extension_solve(&mut self, a: f64, out: &mut [Spec; 3], w: &mut ImexWork, stage: usize) -> (usize, f64, [f64; 3]) {
        let surface = self.surface.as_ref().unwrap();
        let pts = &surface.points;
        let nm = pts.len();
        let nb = surface.centres().len();
        let lsq = surface.basis.is_some();
        if w.surf.is_none() {
            let mut sw = SurfaceWork::new(&self.grid, nm);
            sw.a_prev = std::array::from_fn(|_| vec![0.0; 3 * nb]);
            w.surf = Some(sw);
        }
        let sw = w.surf.as_mut().unwrap();
        if sw.ext.is_none() {
            sw.ext = Some(surface.extension_op(&self.grid, surface.extension.unwrap(), surface.spacing));
        }
        let (g, k2, nu, fft) = (&self.grid, &self.k2, self.nu, &self.fft);
        let keep_mean = self.body_force.is_some();
        let SurfaceWork { nufft, ext, a_prev, phys, spec, .. } = sw;
        let ext = ext.as_ref().unwrap();

        // physical extension forcing E·A (component c uses A[c·nb .. (c+1)·nb])
        let build = |amp: &[f64], phys: &mut [Vec<f64>; 3]| {
            for c in 0..3 {
                phys[c].iter_mut().for_each(|v| *v = 0.0);
                let ac = &amp[c * nb..(c + 1) * nb];
                for (pidx, coup) in &ext.rows {
                    phys[c][*pidx] = coup.iter().map(|&(m, wgt)| wgt * ac[m]).sum();
                }
            }
        };
        // spectral a H⁻¹ P (E A) into spec
        let to_spec = |phys: &[Vec<f64>; 3], spec: &mut [Spec; 3]| {
            for c in 0..3 {
                fft.forward(&phys[c], &mut spec[c]);
            }
            project_with(g, spec, keep_mean);
            for c in 0..3 {
                spec[c].iter_mut().zip(k2.iter()).for_each(|(v, &kk)| *v *= a / (1.0 + a * nu * kk));
            }
        };
        let mut b = vec![0.0; 3 * nm];
        for c in 0..3 {
            let v = nufft.interpolate(&out[c], pts);
            for m in 0..nm {
                b[c * nm + m] = surface.target(m, c) - v[m];
            }
        }
        let x = &mut a_prev[stage];
        let n3 = g.real_len() as f64;
        let two_pi3 = (2.0 * std::f64::consts::PI).powi(3);
        let (matvecs, resid) = if lsq {
            // Over-determined: minimise the wall slip over all collocation
            // markers by CGLS. Aᵀ = Eᵀ · (1/N³)F⁻¹ · aH⁻¹P · (2π)³ S, the exact
            // adjoint of A = J · aH⁻¹P · F · E for Euclidean inner products.
            let mut fwd = |amp: &[f64], res: &mut [f64], phys: &mut [Vec<f64>; 3], spec: &mut [Spec; 3], nufft: &mut Nufft| {
                build(amp, phys);
                to_spec(phys, spec);
                for c in 0..3 {
                    let v = nufft.interpolate(&spec[c], pts);
                    res[c * nm..(c + 1) * nm].copy_from_slice(&v);
                }
            };
            let mut adj = |y: &[f64], res: &mut [f64], phys: &mut [Vec<f64>; 3], spec: &mut [Spec; 3], nufft: &mut Nufft| {
                for c in 0..3 {
                    nufft.spread(&y[c * nm..(c + 1) * nm], pts, &mut spec[c]);
                    spec[c].iter_mut().for_each(|v| *v *= two_pi3);
                }
                project_with(g, spec, keep_mean);
                for c in 0..3 {
                    spec[c].iter_mut().zip(k2.iter()).for_each(|(v, &kk)| *v *= a / (1.0 + a * nu * kk) / n3);
                    fft.inverse_inplace(&mut spec[c], &mut phys[c]);
                }
                res.iter_mut().for_each(|v| *v = 0.0);
                for c in 0..3 {
                    for (pidx, coup) in &ext.rows {
                        let gp = phys[c][*pidx];
                        for &(m, wgt) in coup {
                            res[c * nb + m] += wgt * gp;
                        }
                    }
                }
            };
            cgls(&mut fwd, &mut adj, &b, x, 3 * nb, phys, spec, nufft, Self::EXT_TOL, Self::EXT_MAX_MATVECS)
        } else {
            let mut apply = |amp: &[f64], res: &mut [f64]| {
                build(amp, phys);
                to_spec(phys, spec);
                for c in 0..3 {
                    let v = nufft.interpolate(&spec[c], pts);
                    res[c * nm..(c + 1) * nm].copy_from_slice(&v);
                }
            };
            // The boundary system has many near-null directions (fine
            // tangential detail deep in the shell barely reaches the wall),
            // so GMRES stagnates near 1e-9; 1e-8 is converged for practical
            // purposes and caps the cost.
            crate::krylov::gmres(&mut apply, &b, x, 60, Self::EXT_TOL, Self::EXT_MAX_MATVECS)
        };
        build(x, phys);
        let dv = g.dx().powi(3);
        let force = [0, 1, 2].map(|c| -phys[c].iter().sum::<f64>() * dv);
        to_spec(phys, spec);
        for c in 0..3 {
            out[c].iter_mut().zip(spec[c].iter()).for_each(|(o, v)| *o += *v);
        }
        (matvecs, resid, force)
    }
}

/// CGLS for min ‖A x − b‖ with A given by `fwd` and Aᵀ by `adj` (workspaces
/// threaded through). Returns (products with A or Aᵀ, final ‖b − Ax‖/‖b‖).
#[allow(clippy::too_many_arguments)]
fn cgls(
    fwd: &mut impl FnMut(&[f64], &mut [f64], &mut [Vec<f64>; 3], &mut [Spec; 3], &mut Nufft),
    adj: &mut impl FnMut(&[f64], &mut [f64], &mut [Vec<f64>; 3], &mut [Spec; 3], &mut Nufft),
    b: &[f64],
    x: &mut [f64],
    nx: usize,
    phys: &mut [Vec<f64>; 3],
    spec: &mut [Spec; 3],
    nufft: &mut Nufft,
    tol: f64,
    max_products: usize,
) -> (usize, f64) {
    let dot = |u: &[f64], v: &[f64]| u.iter().zip(v).map(|(p, q)| p * q).sum::<f64>();
    let ny = b.len();
    let bnorm = dot(b, b).sqrt().max(1e-300);
    let mut ax = vec![0.0; ny];
    fwd(x, &mut ax, phys, spec, nufft);
    let mut r: Vec<f64> = b.iter().zip(&ax).map(|(p, q)| p - q).collect();
    let mut s = vec![0.0; nx];
    adj(&r, &mut s, phys, spec, nufft);
    let mut p = s.clone();
    let mut gamma = dot(&s, &s);
    let gamma0 = gamma.max(1e-300);
    let mut products = 2;
    let mut q = vec![0.0; ny];
    let mut history: Vec<f64> = Vec::new();
    while products + 2 <= max_products {
        let rnorm = dot(&r, &r).sqrt() / bnorm;
        // stop on small slip, a converged normal equation, or stagnation
        // Stagnation: < 2% progress over 5 iterations. Warm starts usually
        // land at the reachable residual floor, and a 20-iteration window
        // spent 40 products per solve just confirming that.
        let stagnant = history.len() >= 5 && rnorm > 0.98 * history[history.len() - 5];
        if rnorm < tol || gamma < tol * tol * gamma0 || stagnant {
            break;
        }
        history.push(rnorm);
        fwd(&p, &mut q, phys, spec, nufft);
        let alpha = gamma / dot(&q, &q).max(1e-300);
        x.iter_mut().zip(&p).for_each(|(xv, pv)| *xv += alpha * pv);
        r.iter_mut().zip(&q).for_each(|(rv, qv)| *rv -= alpha * qv);
        adj(&r, &mut s, phys, spec, nufft);
        products += 2;
        let gamma_new = dot(&s, &s);
        let beta = gamma_new / gamma;
        gamma = gamma_new;
        p.iter_mut().zip(&s).for_each(|(pv, sv)| *pv = sv + beta * *pv);
    }
    (products, dot(&r, &r).sqrt() / bnorm)
}

impl Solver {
    /// Force on each part of the installed surface (last implicit stage).
    pub fn surface_part_forces(&self) -> Option<Vec<[f64; 3]>> {
        let surface = self.surface.as_ref()?;
        let sw = self.imex.as_ref()?.surf.as_ref()?;
        Some(match &sw.ext {
            Some(op) => surface.part_forces_extension(op, &sw.a_prev[3]),
            None => {
                let nm = surface.points.len();
                let mut flat = vec![0.0; 3 * nm];
                for c in 0..3 {
                    flat[c * nm..(c + 1) * nm].copy_from_slice(&sw.f_prev[3][c]);
                }
                surface.part_forces_point(&flat)
            }
        })
    }
}
