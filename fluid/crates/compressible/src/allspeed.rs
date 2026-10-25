//! Mach-adaptive all-speed time stepping (Kwatra, Su, Grétarsson & Fedkiw
//! 2009, "A method for avoiding the acoustic time step restriction in
//! compressible flow"), in the form used here for any `System` with a
//! pressure (Euler, the five-equation multiphase model, …).
//!
//! Each step of size Δt is
//!  1. advection without pressure, SSP-RK3: U* from ∂U/∂t + ∇·(U u) = 0
//!     with WENO5-Z reconstruction of the primitives and a Rusanov flux
//!     whose dissipation scales with |u| only (not |u| + c — the acoustic
//!     dissipation that ruins compressible schemes at low Mach); the
//!     pressure p_a and volume fractions α are advected along
//!     (∂q/∂t + u·∇q = 0);
//!  2. one implicit pressure equation from Dp/Dt = −ρc²∇·u, θ-weighted per
//!     face: with p_θ = θ p + (1 − θ) p_a, face velocities
//!     u_f = u*_f − Δt/ρ_f ∂p_θ/∂n and u_θ = θ u_f + (1 − θ) u*_f,
//!       p/(ρc²Δt²) − ∇·(θ²∇p/ρ_f) = p_a/(ρc²Δt²)
//!                                   − ∇·(u*_f − θ(1−θ)Δt ∇p_a/ρ_f)/Δt,
//!     symmetric positive definite, solved by Jacobi-preconditioned CG;
//!  3. conservative corrections with the face values:
//!     ρu −= Δt ∇p_θ, E −= Δt ∇·(p_θ u_θ), α += Δt K ∇·u_θ (Kapila).
//!
//! θ follows the face's acoustic Courant number C = cΔt/Δx: θ = ½ (Crank–
//! Nicolson: resolved sound waves propagate undamped, second order) for
//! C ≤ 1, θ = 1 − 1/(2C) → 1 (backward Euler: unresolved sound is damped,
//! the projection enforces ∇·u → 0) for C ≫ 1.
//!
//! The time step is set by the flow speed, |u|Δt/Δx ≤ CFL, not by the
//! sound speed. The local acoustic Courant number C = cΔt/Δx = CFL/M
//! decides the regime cell by cell: where C ≫ 1 (low Mach) the pressure
//! equation is dominated by its Laplacian and step 2 is an incompressible
//! projection; where C ≲ 1 it is diagonal and the scheme is an explicit
//! compressible one. So the switch between incompressible and compressible
//! physics follows the local Mach number inside one conservative scheme —
//! no interface between two solvers to leak mass or energy.
//!
//! Boundary faces: walls (u_f = 0), outflow (∂p/∂n = 0), Dirichlet (p_f =
//! far-field pressure, a pressure outlet) and characteristic (impedance
//! condition p_f − p_∞ = ±ρc u_f, non-reflecting for normal acoustic waves).

use crate::scheme::weno5z;
use crate::solver::{Bc, Mesh, Solver};
use crate::system::{System, V};
use rayon::prelude::*;

#[derive(Clone, Copy, Debug)]
pub struct AllSpeedOpts {
    /// advective Courant number |u|Δt/Δx (summed over axes)
    pub cfl: f64,
    /// cap on the acoustic Courant number (Σ(|u| + c)/Δx)Δt; ∞ = none,
    /// ≤ 1 resolves acoustics (then the scheme is an explicit one)
    pub max_acoustic_cfl: f64,
    /// relative residual of the pressure solve
    pub tol: f64,
    pub max_iter: usize,
    /// fixed implicitness θ ∈ [½, 1] for every face (None: Mach-adaptive;
    /// Some(1.0) is the backward-Euler scheme of Kwatra et al.)
    pub theta: Option<f64>,
    /// first-order acoustic (Rusanov-like) face-velocity dissipation where
    /// sound is resolved (off: resolved sound waves are not damped)
    pub acoustic_damping: bool,
    /// take explicit steps instead when the all-speed Δt is below this
    /// multiple of the explicit one (an all-speed step costs ~2 explicit;
    /// 0 forces the all-speed scheme)
    pub min_gain: f64,
    /// pressure-based energy, balanced-force momentum and face-velocity
    /// memory also in pure-phase cells with acoustic Courant > 1 (off: only
    /// diffuse-interface cells of the multiphase models use them)
    pub regime_blend: bool,
    /// acoustic Courant number from which pure-phase cells blend to the
    /// elliptic pressure (fully at twice it)
    pub incompressible_courant: f64,
}

impl Default for AllSpeedOpts {
    fn default() -> Self {
        AllSpeedOpts { cfl: 0.4, max_acoustic_cfl: f64::INFINITY, tol: 1e-10, max_iter: 5000, theta: None, acoustic_damping: false, min_gain: 3.0, regime_blend: false, incompressible_courant: 64.0 }
    }
}

/// Per-run statistics of the all-speed mode.
#[derive(Clone, Copy, Debug, Default)]
pub struct AllSpeedStats {
    pub solves: usize,
    pub cg_iterations: usize,
    /// cells × stages in the incompressible regime (acoustic Courant > 1)
    pub incompressible_cells: usize,
    pub cells: usize,
    pub fallback_faces: usize,
    pub max_acoustic_cfl: f64,
    /// largest |Δu|/Δt of the last step (acceleration for the next Δt)
    pub last_accel: f64,
    /// steps taken explicitly because the all-speed Δt gained too little
    pub explicit_steps: usize,
    /// steps redone with explicit sub-steps (all-speed step inadmissible)
    pub explicit_fallbacks: usize,
    pub explicit_substeps: usize,
    /// total energy changed by taking the elliptic pressure in the
    /// incompressible regime (∫ΔE dV, summed over steps; 0 where C ≤ 1)
    pub energy_change: f64,
}

/// Affine boundary face velocity u_f = a + b p_cell and face pressure
/// p_f = c + e p_cell for one boundary face.
#[derive(Clone, Copy)]
struct Affine {
    a: f64,
    b: f64,
    c: f64,
    e: f64,
}

impl<const N: usize, S: System<N>> Solver<N, S> {
    /// All-speed time step: advective CFL, acceleration limit
    /// Δt ≤ CFL √(Δx ρ / |∇p|), optional acoustic cap, multiphase dilatation
    /// stiffness, and the explicit viscous limit.
    pub(crate) fn allspeed_dt(&self, o: &AllSpeedOpts) -> f64 {
        let m = &self.mesh;
        let sys = &self.sys;
        let dims: Vec<usize> = (0..3).filter(|&d| m.active(d)).collect();
        let pi = pressure_index(sys);
        let dxm = dims.iter().map(|&d| m.dx[d]).fold(f64::INFINITY, f64::min);
        let need_gradient = !(self.allspeed_stats.last_accel > 0.0);
        let need_source = !sys.nc_vars().is_empty();
        let need_acoustic = o.max_acoustic_cfl.is_finite();
        let (adv, acc, ac, src, diffusion) = (0..m.n[0])
            .into_par_iter()
            .map(|i| {
                let (mut ra, mut rc, mut rac, mut rs, mut rd) =
                    (0.0f64, 0.0f64, 0.0f64, 0.0f64, crate::viscous::DiffusionBounds::default());
                for j in 0..m.n[1] {
                    for k in 0..m.n[2] {
                        let c = [i + m.g(0), j + m.g(1), k + m.g(2)];
                        let w = sys.to_prim(&self.u[m.idx(c[0], c[1], c[2])]);
                        let cs = if need_acoustic {
                            sys.sound_speed(&w)
                        } else {
                            0.0
                        };
                        let (mut rate, mut acoustic, mut divu, mut g2) = (0.0, 0.0, 0.0, 0.0);
                        for &d in &dims {
                            let vi = sys.vel_index() + d;
                            rate += w[vi].abs() / m.dx[d];
                            if need_acoustic {
                                acoustic += (w[vi].abs() + cs) / m.dx[d];
                            }
                            if need_gradient {
                                // Evaluate CURRENT boundary states, including the first query.
                                let (wa, wb) = (
                                    self.neighbor_prim(c, d, false),
                                    self.neighbor_prim(c, d, true),
                                );
                                if need_source {
                                    divu += (wb[vi] - wa[vi]).abs() / (2.0 * m.dx[d]);
                                }
                                let gp = (wb[pi] - wa[pi]) / (2.0 * m.dx[d]);
                                g2 += gp * gp;
                            } else if need_source {
                                let (va, vb) = self.neighbor_velocities(c, d);
                                divu += (vb - va).abs() / (2.0 * m.dx[d]);
                            }
                        }
                        ra = ra.max(rate);
                        rac = rac.max(acoustic);
                        if need_gradient {
                            rc = rc.max((g2.sqrt() / (sys.density(&w) * dxm)).sqrt());
                        }
                        if need_source {
                            rs = rs.max(sys.dilatation_rate(&w, divu));
                        }
                        if let Some(tr) = &self.viscous {
                            rd.observe(sys, tr, &w);
                        }
                    }
                }
                (ra, rc, rac, rs, rd)
            })
            .reduce(
                || (0.0, 0.0, 0.0, 0.0, crate::viscous::DiffusionBounds::default()),
                |a, b| {
                    (
                        a.0.max(b.0),
                        a.1.max(b.1),
                        a.2.max(b.2),
                        a.3.max(b.3),
                        a.4.merge(b.4),
                    )
                },
            );
        let acc = if need_gradient {
            acc
        } else {
            (self.allspeed_stats.last_accel / dxm).sqrt()
        };
        let mut dt = o.cfl / adv.max(1e-300);
        dt = dt.min(o.cfl / acc.max(1e-300));
        if need_acoustic {
            dt = dt.min(o.max_acoustic_cfl / ac);
        }
        if need_source && src > 0.0 {
            dt = dt.min(0.5 / src);
        }
        dt = dt.min(self.viscous_dt(diffusion));
        dt
    }

    /// One all-speed step: SSP-RK3 advection, then one pressure solve.
    pub(crate) fn allspeed_step(&mut self, dt: f64, o: &AllSpeedOpts) {
        let m = self.mesh.clone();
        let interior = interior(&m);
        let ncell = m.cells();
        let pi = pressure_index(&self.sys);
        let u0 = self.u.clone();
        let mut uu = self.u.clone();
        self.fill_ghosts(&mut uu);
        // pressure carried from the last solve where sound is unresolved
        // (the EOS pressure of a stiff phase at C ≫ 1 is round-off noise of
        // the conserved energy, which itself stays exact), the EOS pressure
        // where it is resolved
        let mut p0 = vec![0.0; m.len()];
        let carried = self.allspeed_p.len() == interior.len();
        for (l, &q) in interior.iter().enumerate() {
            let pe = self.sys.to_prim(&uu[q])[pi];
            p0[q] = if carried {
                let (pc, w) = self.allspeed_p[l];
                w * pe + (1.0 - w) * pc
            } else {
                pe
            };
        }
        let mut pp = p0.clone();
        let mut fallback = 0;
        let phi_old = if self.scalars.is_empty() { Vec::new() } else { self.scalar_phis(&uu) };
        if !self.scalars.is_empty() {
            self.reset_flux_acc();
        }
        // effective face velocity of the RK3 transport (weights 1/6, 1/6, 2/3)
        let mut ueff: [Vec<f64>; 3] = Default::default();
        // SSP-RK3 on (U, p_a): stage weights (α_old, α_new), transport weight
        // advection ghosts: characteristic (acoustic) boundaries become
        // zero-gradient — the advected stage states carry no meaningful
        // pressure (no pressure work yet), and the acoustic boundary
        // condition is the impedance law of the pressure equation
        let bc_adv: [[Bc<N>; 2]; 3] = self.bc.map(|p| p.map(|b| if let Bc::Characteristic(_) = b { Bc::Outflow } else { b }));
        for (a0, a1, bw) in [(0.0, 1.0, 1.0 / 6.0), (0.75, 0.25, 1.0 / 6.0), (1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0)] {
            self.fill_ghosts_with(&mut uu, &bc_adv);
            self.fill_scalar_ghosts(&mut pp);
            let (adv, padv, fb, ub, cf) = self.advection_rates(&uu, &pp);
            fallback += fb;
            for (ci, per) in cf.iter().enumerate() {
                for d in 0..3 {
                    for (a, v) in self.flux_acc[ci][d].iter_mut().zip(per[d].iter()) {
                        *a += bw * v;
                    }
                }
            }
            for d in 0..3 {
                if ueff[d].is_empty() {
                    ueff[d] = vec![0.0; ub[d].len()];
                }
                for (e, v) in ueff[d].iter_mut().zip(ub[d].iter()) {
                    *e += bw * v;
                }
            }
            for (l, &q) in interior.iter().enumerate() {
                for k in 0..N {
                    uu[q][k] = a0 * u0[q][k] + a1 * (uu[q][k] + dt * adv[l][k]);
                }
                pp[q] = a0 * p0[q] + a1 * (pp[q] + dt * padv[l]);
            }
        }
        self.fill_ghosts_with(&mut uu, &bc_adv);
        let pa: Vec<f64> = interior.iter().map(|&q| pp[q]).collect();
        // kinetics (phase change, dissolution) on the advected state, before
        // the pressure solve: its pressure change at fixed volume,
        // Δp = −ρc²ΔV, enters the pressure equation as a volume source, so
        // the projection drives the flow that supplies (or removes) the
        // volume the mass transfer frees (a dissolving bubble draws liquid
        // in instead of stretching it)
        let mut dpk = vec![0.0; interior.len()];
        let mut rc2_pre = vec![0.0; interior.len()];
        if self.kinetics.is_some() || self.conduction.is_some() {
            // ρc² at the advected pressure, and the EOS pressure before
            let sys = &self.sys;
            let pre: Vec<(f64, f64)> = (0..interior.len())
                .into_par_iter()
                .map(|l| {
                    let mut w0 = sys.to_prim(&uu[interior[l]]);
                    let pe = w0[pi];
                    w0[pi] = pa[l];
                    (sys.density(&w0) * sys.sound_speed(&w0).powi(2), pe)
                })
                .collect();
            let p_before: Vec<f64> = pre.iter().map(|x| x.1).collect();
            for (l, (rc2, _)) in pre.into_iter().enumerate() {
                rc2_pre[l] = rc2;
            }
            // kinetics first, on the state at the advected pressure p_a (its
            // rates depend on p; the cells it acts on are interface cells,
            // pressure-based anyway); its pressure change counts from p_a
            let mut dp_kin = vec![0.0; interior.len()];
            let mut kin_cell = vec![false; interior.len()];
            if let Some(kin) = &self.kinetics {
                let ns = self.scalars.len();
                let sys = &self.sys;
                let res: Vec<(usize, V<N>, Vec<f64>)> = (0..interior.len())
                    .into_par_iter()
                    .filter_map(|l| {
                        let q = interior[l];
                        let mut w0 = sys.to_prim(&uu[q]);
                        w0[pi] = pa[l];
                        let mut uc = sys.to_cons(&w0);
                        let mut sc: Vec<f64> = (0..ns).map(|si| self.scalars[si].u[q]).collect();
                        kin(sys, &mut uc, &mut sc, dt).then_some((l, uc, sc))
                    })
                    .collect();
                self.kinetics_cells = res.len();
                for (l, uc, sc) in res {
                    let q = interior[l];
                    uu[q] = uc;
                    for si in 0..ns {
                        self.scalars[si].u[q] = sc[si];
                    }
                    dp_kin[l] = sys.to_prim(&uc)[pi] - pa[l];
                    kin_cell[l] = true;
                }
            }
            // conduction on the conservative state (its temperature is the
            // energy's); its pressure change counts from before it
            let p_mid: Vec<f64> = interior.iter().map(|&q| self.sys.to_prim(&uu[q])[pi]).collect();
            if self.conduction.is_some() {
                self.conduct(&mut uu, dt);
            }
            // pressure change at fixed volume → volume source of the solve
            let sys = &self.sys;
            let _ = &p_before;
            for (l, &q) in interior.iter().enumerate() {
                let dp_cond = sys.to_prim(&uu[q])[pi] - p_mid[l];
                dpk[l] = dp_cond + if kin_cell[l] { dp_kin[l] } else { 0.0 };
                if !dpk[l].is_finite() {
                    dpk[l] = 0.0;
                }
            }
            self.fill_ghosts_with(&mut uu, &bc_adv);
        }
        let _ = ncell;
        self.pressure_stage(uu, pa, ueff, dt, o, fallback, &dpk, &rc2_pre);
        self.finish_scalars_with(&phi_old, dt, false);
        let sys = &self.sys;
        let vi = sys.vel_index();
        let dims: Vec<usize> = (0..3).filter(|&d| m.active(d)).collect();
        let accel = interior
            .par_iter()
            .map(|&q| {
                let (a, b) = (sys.to_prim(&u0[q]), sys.to_prim(&self.u[q]));
                dims.iter().map(|&d| (b[vi + d] - a[vi + d]).powi(2)).sum::<f64>().sqrt() / dt
            })
            .reduce(|| 0.0, f64::max);
        self.allspeed_stats.last_accel = accel;
    }

    /// Ghost cells of a scalar field (the advected pressure): periodic,
    /// zero-gradient (outflow, characteristic), mirror (wall), far value
    /// (Dirichlet).
    fn fill_scalar_ghosts(&self, f: &mut [f64]) {
        let m = &self.mesh;
        let pi = pressure_index(&self.sys);
        for d in 0..3 {
            if !m.active(d) {
                continue;
            }
            let (n, ng) = (m.n[d], m.ng);
            let others: Vec<usize> = (0..3).filter(|&e| e != d).collect();
            let (e0, e1) = (others[0], others[1]);
            for a in 0..m.m(e0) {
                for b in 0..m.m(e1) {
                    let at = |c: usize| {
                        let mut ijk = [0; 3];
                        ijk[d] = c;
                        ijk[e0] = a;
                        ijk[e1] = b;
                        m.idx(ijk[0], ijk[1], ijk[2])
                    };
                    for gi in 0..ng {
                        f[at(gi)] = match self.bc[d][0] {
                            Bc::Periodic => f[at(n + gi)],
                            Bc::Reflect => f[at(2 * ng - 1 - gi)],
                            Bc::Dirichlet(v) => self.sys.to_prim(&v)[pi],
                            Bc::Outflow | Bc::Characteristic(_) => f[at(ng)],
                        };
                        f[at(n + ng + gi)] = match self.bc[d][1] {
                            Bc::Periodic => f[at(ng + gi)],
                            Bc::Reflect => f[at(n + ng - 1 - gi)],
                            Bc::Dirichlet(v) => self.sys.to_prim(&v)[pi],
                            Bc::Outflow | Bc::Characteristic(_) => f[at(n + ng - 1)],
                        };
                    }
                }
            }
        }
    }

    /// Advection rates (interior, linear index) of U (−∇·(U u), plus α∇·ū
    /// for the non-conservative volume fractions) and of the pressure field
    /// p (−u·∇p); ghosts of both must be filled.
    #[allow(clippy::type_complexity)]
    fn advection_rates(&self, u: &[V<N>], pfield: &[f64]) -> (Vec<V<N>>, Vec<f64>, usize, [Vec<f64>; 3], Vec<[Vec<f64>; 3]>) {
        let carriers = self.carriers();
        let cidx: Vec<std::ops::Range<usize>> = carriers.iter().map(|c| crate::scalars::carrier_indices(&self.sys, *c)).collect();
        let m = &self.mesh;
        let sys = &self.sys;
        let dims: Vec<usize> = (0..3).filter(|&d| m.active(d)).collect();
        let vi = sys.vel_index();
        let nc = sys.nc_vars();
        let ncell = m.cells();
        let lin = |c: [usize; 3]| (c[0] * m.n[1] + c[1]) * m.n[2] + c[2];
        let w: Vec<V<N>> = u.par_iter().map(|x| sys.to_prim(x)).collect();
        let mut adv = vec![[0.0; N]; ncell];
        let mut padv = vec![0.0; ncell];
        let mut fallback = 0usize;
        let mut ubar: [Vec<f64>; 3] = Default::default();
        let mut cflux: Vec<[Vec<f64>; 3]> = carriers.iter().map(|_| Default::default()).collect();
        for &d in &dims {
            let others: Vec<usize> = (0..3).filter(|&e| e != d).collect();
            let (e0, e1) = (others[0], others[1]);
            let pencils: Vec<(usize, usize)> = (0..m.n[e0]).flat_map(|a| (0..m.n[e1]).map(move |b| (a, b))).collect();
            let inv = 1.0 / m.dx[d];
            let res: Vec<(Vec<(usize, V<N>, f64)>, usize, Vec<f64>, Vec<Vec<f64>>)> = pencils
                .par_iter()
                .map(|&(a, b)| {
                    let at = |cc: usize| {
                        let mut ijk = [0; 3];
                        ijk[d] = cc;
                        ijk[e0] = a + m.g(e0);
                        ijk[e1] = b + m.g(e1);
                        m.idx(ijk[0], ijk[1], ijk[2])
                    };
                    let ng = m.ng;
                    let n = m.n[d];
                    let line: Vec<V<N>> = (0..m.m(d)).map(|cc| w[at(cc)]).collect();
                    let pl: Vec<f64> = (0..m.m(d)).map(|cc| pfield[at(cc)]).collect();
                    // faces f = 0..n between line cells ng−1+f and ng+f
                    let mut fl = vec![[0.0; N]; n + 1];
                    let mut fp = vec![0.0; n + 1];
                    let mut ub = vec![0.0; n + 1];
                    let mut fb = 0;
                    for f in 0..=n {
                        let i = ng - 1 + f;
                        let mut wl = [0.0; N];
                        let mut wr = [0.0; N];
                        for q in 0..N {
                            wl[q] = weno5z(line[i - 2][q], line[i - 1][q], line[i][q], line[i + 1][q], line[i + 2][q]);
                            wr[q] = weno5z(line[i + 3][q], line[i + 2][q], line[i + 1][q], line[i][q], line[i - 1][q]);
                        }
                        let (mut p_lf, mut p_rf) = (weno5z(pl[i - 2], pl[i - 1], pl[i], pl[i + 1], pl[i + 2]), weno5z(pl[i + 3], pl[i + 2], pl[i + 1], pl[i], pl[i - 1]));
                        if !(sys.admissible(&wl) && sys.admissible(&wr)) {
                            wl = line[i];
                            wr = line[i + 1];
                            p_lf = pl[i];
                            p_rf = pl[i + 1];
                            fb += 1;
                        }
                        let (ul, ur) = (sys.to_cons(&wl), sys.to_cons(&wr));
                        let (vl, vr) = (wl[vi + d], wr[vi + d]);
                        let s = vl.abs().max(vr.abs());
                        for q in 0..N {
                            fl[f][q] = 0.5 * (vl * ul[q] + vr * ur[q]) - 0.5 * s * (ur[q] - ul[q]);
                        }
                        fp[f] = 0.5 * (vl * p_lf + vr * p_rf) - 0.5 * s * (p_rf - p_lf);
                        ub[f] = 0.5 * (vl + vr);
                    }
                    let mut out = Vec::with_capacity(n);
                    for cc in 0..n {
                        let wc = &line[cc + ng];
                        let mut r = [0.0; N];
                        for q in 0..N {
                            r[q] = -(fl[cc + 1][q] - fl[cc][q]) * inv;
                        }
                        let divu = (ub[cc + 1] - ub[cc]) * inv;
                        for q in nc.clone() {
                            r[q] += wc[q] * divu;
                        }
                        let rp = -(fp[cc + 1] - fp[cc]) * inv + pl[cc + ng] * divu;
                        let mut c = [0usize; 3];
                        c[d] = cc;
                        c[e0] = a;
                        c[e1] = b;
                        out.push((lin(c), r, rp));
                    }
                    let cf: Vec<Vec<f64>> = cidx.iter().map(|r| (0..=n).map(|f| r.clone().map(|k| fl[f][k]).sum()).collect()).collect();
                    (out, fb, ub, cf)
                })
                .collect();
            for (cells, fb, ub, cf) in res {
                fallback += fb;
                ubar[d].extend_from_slice(&ub);
                for (ci, v) in cf.into_iter().enumerate() {
                    cflux[ci][d].extend_from_slice(&v);
                }
                for (q, r, rp) in cells {
                    for k in 0..N {
                        adv[q][k] += r[k];
                    }
                    padv[q] += rp;
                }
            }
        }
        // explicit extra physics with the advection
        if self.viscous.is_some() || self.surface_tension.is_some() {
            let mut extra = vec![[0.0; N]; m.len()];
            if let Some(tr) = &self.viscous {
                crate::viscous::add_viscous(sys, m, tr, u, &mut extra);
            }
            if let Some(st) = &self.surface_tension {
                st.add(sys, m, u, &mut extra);
            }
            for (l, q) in interior(m).into_iter().enumerate() {
                for k in 0..N {
                    adv[l][k] += extra[q][k];
                }
            }
        }
        (adv, padv, fallback, ubar, cflux)
    }

    /// Implicit pressure step and conservative corrections from the
    /// advected state (interior, linear index); writes `self.u`.
    #[allow(clippy::too_many_arguments)]
    fn pressure_stage(&mut self, ufull: Vec<V<N>>, pa: Vec<f64>, ueff: [Vec<f64>; 3], dt: f64, o: &AllSpeedOpts, fallback: usize, dpk: &[f64], rc2_pre: &[f64]) {
        let m = self.mesh.clone();
        let sys = &self.sys;
        let dims: Vec<usize> = (0..3).filter(|&d| m.active(d)).collect();
        let pi = pressure_index(sys);
        let vi = sys.vel_index();
        let nc = sys.nc_vars();
        let ncell = m.cells();
        let lin = |c: [usize; 3]| (c[0] * m.n[1] + c[1]) * m.n[2] + c[2];
        let interior = interior(&m);
        let ustar: Vec<V<N>> = interior.iter().map(|&q| ufull[q]).collect();
        // primitives of the star state with the advected pressure
        let wstar: Vec<V<N>> = ustar
            .par_iter()
            .zip(pa.par_iter())
            .map(|(x, &p)| {
                let mut ws = sys.to_prim(x);
                ws[pi] = p;
                ws
            })
            .collect();
        let rho: Vec<f64> = wstar.iter().map(|x| sys.density(x)).collect();
        let rc2: Vec<f64> = wstar
            .iter()
            .enumerate()
            .map(|(l, x)| if dpk[l] != 0.0 && rc2_pre[l] > 0.0 { rc2_pre[l] } else { sys.density(x) * sys.sound_speed(x).powi(2) })
            .collect();
        // ---- 2. pressure equation
        let nb = |c: [usize; 3], d: usize, s: i64| -> Option<[usize; 3]> {
            let mut cc = c;
            let n = m.n[d] as i64;
            let v = c[d] as i64 + s;
            if v < 0 || v >= n {
                if self.bc[d][(v >= n) as usize] == Bc::Periodic {
                    cc[d] = v.rem_euclid(n) as usize;
                    return Some(cc);
                }
                return None;
            }
            cc[d] = v as usize;
            Some(cc)
        };
        let coords = |l: usize| [l / (m.n[1] * m.n[2]), (l / m.n[2]) % m.n[1], l % m.n[2]];
        // + face of each cell: u*_f, 1/(ρ_f Δx²), θ, θ²/(ρ_f Δx²) and
        // ũ_f = u*_f − θ(1−θ)Δt ∂p_a/∂n / ρ_f
        let mut ustar_f = vec![[0.0; 3]; ncell];
        let mut coef = vec![[0.0; 3]; ncell];
        let mut theta = vec![[1.0; 3]; ncell];
        let mut cface = vec![[0.0; 3]; ncell];
        let mut coefi = vec![[0.0; 3]; ncell];
        let mut ut = vec![[0.0; 3]; ncell];
        // boundary faces: (cell, axis, side) → affine law
        let far_p = |d: usize, s: usize| -> Option<f64> {
            match self.bc[d][s] {
                Bc::Dirichlet(v) | Bc::Characteristic(v) => Some(sys.to_prim(&v)[pi]),
                _ => None,
            }
        };
        let mut bface: Vec<(usize, usize, usize, Affine)> = Vec::new();
        for l in 0..ncell {
            let c = coords(l);
            for &d in &dims {
                if let Some(cn) = nb(c, d, 1) {
                    let ln = lin(cn);
                    let uc = wstar[l][vi + d];
                    let un = wstar[ln][vi + d];
                    // average of the cell velocities, plus (incompressible
                    // regime only) the memory of the last projected face
                    // velocity — the balanced-force (MAC-like) face velocity
                    // at density jumps
                    ustar_f[l][d] = 0.5 * (uc + un) + self.allspeed_delta.get(l).map_or(0.0, |v| v[d]);
                    coef[l][d] = 2.0 / ((rho[l] + rho[ln]) * m.dx[d] * m.dx[d]);
                    let th = o.theta.unwrap_or_else(|| {
                        let cf = 0.5 * ((rc2[l] / rho[l]).sqrt() + (rc2[ln] / rho[ln]).sqrt()) * dt / m.dx[d];
                        if cf <= 1.0 { 0.5 } else { 1.0 - 0.5 / cf }
                    });
                    theta[l][d] = th;
                    cface[l][d] = 0.5 * ((rc2[l] / rho[l]).sqrt() + (rc2[ln] / rho[ln]).sqrt()) * dt / m.dx[d];
                    // acoustic (HLLC-like) face-velocity dissipation
                    // −Δp_a/(ρ_L c_L + ρ_R c_R), faded by 1/C where the
                    // implicit projection supplies the pressure–velocity
                    // coupling: where sound is resolved (small C, e.g. the
                    // low Wood sound speed of diffuse multiphase interfaces)
                    // the projection's Δt/(ρΔx) coupling is C times weaker
                    // than a Riemann solver's and the Kapila mixture modes
                    // are otherwise undamped. O(M u) at low Mach.
                    if o.acoustic_damping {
                        let imp = (rc2[l] * rho[l]).sqrt() + (rc2[ln] * rho[ln]).sqrt();
                        ustar_f[l][d] -= (1.0f64).min(1.0 / cface[l][d].max(1e-300)) * (pa[ln] - pa[l]) / imp;
                    }
                    coefi[l][d] = th * th * coef[l][d];
                    ut[l][d] = ustar_f[l][d] - th * (1.0 - th) * dt * coef[l][d] * m.dx[d] * (pa[ln] - pa[l]);
                }
                for s in 0..2usize {
                    if nb(c, d, if s == 1 { 1 } else { -1 }).is_some() {
                        continue;
                    }
                    let uc = wstar[l][vi + d];
                    let beta = 2.0 * dt / (rho[l] * m.dx[d]);
                    let sg = if s == 1 { 1.0 } else { -1.0 };
                    // u_f = u* − sg β (p_f − p_c) along the outward axis sign
                    let law = match self.bc[d][s] {
                        Bc::Reflect => Affine { a: 0.0, b: 0.0, c: 0.0, e: 1.0 },
                        Bc::Outflow => Affine { a: uc, b: 0.0, c: 0.0, e: 1.0 },
                        Bc::Dirichlet(_) => {
                            let pf = far_p(d, s).unwrap();
                            Affine { a: uc - sg * beta * pf, b: sg * beta, c: pf, e: 0.0 }
                        }
                        Bc::Characteristic(_) => {
                            // p_f − p∞ = sg Z u_f, Z = ρc
                            let pf = far_p(d, s).unwrap();
                            let z = (rc2[l] * rho[l]).sqrt();
                            let den = 1.0 + z * beta;
                            let a = (uc - sg * beta * pf) / den;
                            let b = sg * beta / den;
                            Affine { a, b, c: pf + sg * z * a, e: sg * z * b }
                        }
                        Bc::Periodic => unreachable!(),
                    };
                    bface.push((l, d, s, law));
                }
            }
        }
        let mut diag: Vec<f64> = (0..ncell).map(|l| 1.0 / (rc2[l] * dt * dt)).collect();
        let mut rhs: Vec<f64> = (0..ncell).map(|l| pa[l] / (rc2[l] * dt * dt)).collect();
        for l in 0..ncell {
            let c = coords(l);
            for &d in &dims {
                let inv = 1.0 / (dt * m.dx[d]);
                if nb(c, d, 1).is_some() {
                    rhs[l] -= ut[l][d] * inv;
                    diag[l] += coefi[l][d];
                }
                if let Some(cm) = nb(c, d, -1) {
                    let lm = lin(cm);
                    rhs[l] += ut[lm][d] * inv;
                    diag[l] += coefi[lm][d];
                }
            }
        }
        for &(l, d, s, law) in &bface {
            let inv = 1.0 / (dt * m.dx[d]);
            let sg = if s == 1 { 1.0 } else { -1.0 };
            diag[l] += sg * law.b * inv;
            rhs[l] -= sg * law.a * inv;
        }
        // matrix-free operator A p
        let apply = |p: &[f64], out: &mut [f64]| {
            out.par_iter_mut().with_min_len(4096).enumerate().for_each(|(l, o)| {
                let c = coords(l);
                let mut v = diag[l] * p[l];
                for &d in &dims {
                    if let Some(cp) = nb(c, d, 1) {
                        v -= coefi[l][d] * p[lin(cp)];
                    }
                    if let Some(cm) = nb(c, d, -1) {
                        let lm = lin(cm);
                        v -= coefi[lm][d] * p[lm];
                    }
                }
                *o = v;
            });
        };
        // solve for the increment δ = p − p_a: at low Mach p ≈ p₀ + O(ρu²)
        // with p₀ ~ ρc², and a tolerance relative to p₀ would leave the
        // dynamic pressure unconverged
        // r₀ = b − A p_a in difference form: the diagonal p_a/(ρc²Δt²)
        // terms cancel exactly and the Laplacian acts on differences of p_a
        // (forming A p_a would leave ~1e-16 C² p_a of round-off at large
        // acoustic Courant numbers C — Pa-level pressure noise in water)
        let mut r0 = vec![0.0; ncell];
        for l in 0..ncell {
            let c = coords(l);
            let mut r = 0.0;
            for &d in &dims {
                let inv = 1.0 / (dt * m.dx[d]);
                if let Some(cp) = nb(c, d, 1) {
                    let lp = lin(cp);
                    r -= ut[l][d] * inv + coefi[l][d] * (pa[l] - pa[lp]);
                }
                if let Some(cm) = nb(c, d, -1) {
                    let lm = lin(cm);
                    r += ut[lm][d] * inv - coefi[lm][d] * (pa[l] - pa[lm]);
                }
            }
            // volume source of the mass transfer: Δp_kin/(ρc²Δt²)
            r0[l] = r + dpk[l] / (rc2[l] * dt * dt);
        }
        for &(l, d, s, law) in &bface {
            let inv = 1.0 / (dt * m.dx[d]);
            let sg = if s == 1 { 1.0 } else { -1.0 };
            r0[l] -= sg * (law.a + law.b * pa[l]) * inv;
        }
        let _ = &rhs;
        let mut dp = vec![0.0; ncell];
        let its = pcg(&apply, &diag, &r0, &mut dp, o.tol, 0.0, o.max_iter);
        let p: Vec<f64> = pa.iter().zip(dp.iter()).map(|(a, d)| a + d).collect();
        // ---- 3. corrections
        let mut fu = vec![[f64::NAN; 6]; ncell]; // face velocity u_θ [d][side]
        let mut fpr = vec![[f64::NAN; 6]; ncell]; // face pressure p_θ
        let mut fa = vec![[0.0; 6]; ncell]; // face acceleration −∂p_θ/∂n / ρ_f
        let mut unew_f = vec![[0.0; 3]; ncell];
        for l in 0..ncell {
            let c = coords(l);
            for &d in &dims {
                if let Some(cn) = nb(c, d, 1) {
                    let ln = lin(cn);
                    let th = theta[l][d];
                    let g = dt * coef[l][d] * m.dx[d];
                    let unew = ustar_f[l][d] - g * (th * (p[ln] - p[l]) + (1.0 - th) * (pa[ln] - pa[l]));
                    let uf = th * unew + (1.0 - th) * ustar_f[l][d];
                    let pf = th * 0.5 * (p[l] + p[ln]) + (1.0 - th) * 0.5 * (pa[l] + pa[ln]);
                    fu[l][2 * d + 1] = uf;
                    fpr[l][2 * d + 1] = pf;
                    fu[ln][2 * d] = uf;
                    fpr[ln][2 * d] = pf;
                    let acc = (unew - ustar_f[l][d]) / dt;
                    fa[l][2 * d + 1] = acc;
                    fa[ln][2 * d] = acc;
                    unew_f[l][d] = unew;
                }
            }
        }
        for &(l, d, s, law) in &bface {
            fu[l][2 * d + s] = law.a + law.b * p[l];
            fpr[l][2 * d + s] = law.c + law.e * p[l];
            let sg = if s == 1 { 1.0 } else { -1.0 };
            fa[l][2 * d + s] = -sg * (fpr[l][2 * d + s] - p[l]) * 2.0 / (rho[l] * m.dx[d]);
        }
        // re-transport every variable with δu = u_θ − ū (ū: the face velocity
        // the RK3 advection used), so that all dilatation terms — internal
        // energy, volume fractions, pressure work — see the same θ-weighted
        // face velocity as the pressure equation (otherwise the internal-
        // energy part is forward Euler in time and amplifies sound waves)
        let wfull: Vec<V<N>> = ufull.par_iter().map(|x| sys.to_prim(x)).collect();
        let mut tr = vec![[0.0; N]; ncell];
        let mut dflux: Vec<[Vec<f64>; 3]> = self.carriers().iter().map(|_| std::array::from_fn(|d| if m.active(d) { vec![0.0; crate::scalars::face_count(&m, d)] } else { Vec::new() })).collect();
        for &d in &dims {
            let others: Vec<usize> = (0..3).filter(|&e| e != d).collect();
            let (e0, e1) = (others[0], others[1]);
            let n = m.n[d];
            let inv = 1.0 / m.dx[d];
            let pencils: Vec<(usize, usize)> = (0..m.n[e0]).flat_map(|a| (0..m.n[e1]).map(move |b| (a, b))).collect();
            let carriers = self.carriers();
            let cidx: Vec<std::ops::Range<usize>> = carriers.iter().map(|c| crate::scalars::carrier_indices(sys, *c)).collect();
            let res: Vec<(Vec<(usize, V<N>)>, Vec<Vec<f64>>)> = pencils
                .par_iter()
                .enumerate()
                .map(|(pidx, &(a, b))| {
                    let at = |cc: usize| {
                        let mut ijk = [0; 3];
                        ijk[d] = cc;
                        ijk[e0] = a + m.g(e0);
                        ijk[e1] = b + m.g(e1);
                        m.idx(ijk[0], ijk[1], ijk[2])
                    };
                    let cell = |cc: usize| {
                        let mut c = [0usize; 3];
                        c[d] = cc;
                        c[e0] = a;
                        c[e1] = b;
                        lin(c)
                    };
                    let ng = m.ng;
                    let mut fl = vec![[0.0; N]; n + 1];
                    let mut du = vec![0.0; n + 1];
                    for f in 0..=n {
                        // θ-weighted face velocity of this face
                        let uth = if f < n { fu[cell(f)][2 * d] } else { fu[cell(n - 1)][2 * d + 1] };
                        let dlt = uth - ueff[d][pidx * (n + 1) + f];
                        du[f] = dlt;
                        let i = ng - 1 + f;
                        let mut wl = [0.0; N];
                        let mut wr = [0.0; N];
                        for q in 0..N {
                            let v = |c: usize| wfull[at(c)][q];
                            wl[q] = weno5z(v(i - 2), v(i - 1), v(i), v(i + 1), v(i + 2));
                            wr[q] = weno5z(v(i + 3), v(i + 2), v(i + 1), v(i), v(i - 1));
                        }
                        if !(sys.admissible(&wl) && sys.admissible(&wr)) {
                            wl = wfull[at(i)];
                            wr = wfull[at(i + 1)];
                        }
                        let (ul, ur) = (sys.to_cons(&wl), sys.to_cons(&wr));
                        for q in 0..N {
                            fl[f][q] = 0.5 * (ul[q] + ur[q]) * dlt - 0.5 * dlt.abs() * (ur[q] - ul[q]);
                        }
                    }
                    let cells: Vec<(usize, V<N>)> = (0..n)
                        .map(|cc| {
                            let wc = &wfull[at(cc + ng)];
                            let mut r = [0.0; N];
                            for q in 0..N {
                                r[q] = -(fl[cc + 1][q] - fl[cc][q]) * inv;
                            }
                            let ddiv = (du[cc + 1] - du[cc]) * inv;
                            for q in nc.clone() {
                                r[q] += wc[q] * ddiv;
                            }
                            (cell(cc), r)
                        })
                        .collect();
                    let cf: Vec<Vec<f64>> = cidx.iter().map(|r| (0..=n).map(|f| r.clone().map(|k| fl[f][k]).sum()).collect()).collect();
                    let _ = pidx;
                    (cells, cf)
                })
                .collect();
            let mut off = 0;
            for (cells, cf) in res {
                for (l, r) in cells {
                    for k in 0..N {
                        tr[l][k] += r[k];
                    }
                }
                for (ci, v) in cf.into_iter().enumerate() {
                    let acc = &mut dflux[ci][d];
                    for (f, x) in v.iter().enumerate() {
                        acc[off + f] += x;
                    }
                }
                off += n + 1;
            }
        }
        for (ci, per) in dflux.iter().enumerate() {
            for d in 0..3 {
                for (a, v) in self.flux_acc[ci][d].iter_mut().zip(per[d].iter()) {
                    *a += v;
                }
            }
        }
        let mut n_inc = 0usize;
        let mut max_ac: f64 = 0.0;
        let mut e_fix = 0.0;
        let mut wnext = vec![1.0; ncell];
        // acoustic Courant number per cell, and its maximum over ±2 cells
        // per axis: the regime of a cell is that of its neighbourhood (an
        // interface cell has a small Wood sound speed of its own but sits
        // next to stiff liquid)
        let acell: Vec<f64> = (0..ncell).map(|l| dims.iter().map(|&d| (rc2[l] / rho[l]).sqrt() * dt / m.dx[d]).sum()).collect();
        let mut anb = acell.clone();
        for &d in &dims {
            let prev = anb.clone();
            for l in 0..ncell {
                let c = coords(l);
                for s2 in [-2i64, -1, 1, 2] {
                    if let Some(cn) = nb(c, d, s2) {
                        anb[l] = anb[l].max(prev[lin(cn)]);
                    }
                }
            }
        }
        for l in 0..ncell {
            let q = interior[l];
            let mut v = ustar[l];
            for k in 0..N {
                v[k] += dt * tr[l][k];
            }
            let mut divu = 0.0;
            let ac = acell[l];
            // regime weight: 1 where sound is resolved (C ≤ 1), 0 for C ≥ 2
            // diffuse-interface cells of the multiphase models always take
            // the pressure-based update: their mixture (Wood) stiffness is
            // the small difference of two O(γp∞) terms in a conservative
            // energy update, which the projection cannot keep consistent
            let mixed = {
                // two or more phases of very different stiffness (gas and
                // liquid; a solid–liquid front is not such a cell)
                let ph = sys.phase_stiffness(&wstar[l]);
                let present: Vec<f64> = ph.iter().filter(|x| x.0 > 1e-4 && x.0 < 1.0 - 1e-4).map(|x| x.1).collect();
                present.len() >= 2 && present.iter().cloned().fold(0.0, f64::max) > 100.0 * present.iter().cloned().fold(f64::INFINITY, f64::min)
            };
            // deep in the incompressible regime (C ≫ 1) the conservative
            // energy cannot follow the elliptic pressure either: a velocity
            // divergence at the solver's precision moves p by ρc²Δt ∇·u,
            // i.e. ~1e-16 C² p of noise that grows step to step
            let c_inc = o.incompressible_courant;
            let wcons = if mixed {
                0.0
            } else if o.regime_blend {
                (2.0 - anb[l]).clamp(0.0, 1.0)
            } else {
                1.0
            };
            // weight of the EOS pressure in the next advected pressure
            wnext[l] = if mixed { 1.0 } else { ((2.0 * c_inc - anb[l]) / c_inc).clamp(0.0, 1.0) };
            let rho_new = sys.density(&sys.to_prim(&v));
            for &d in &dims {
                let inv = 1.0 / m.dx[d];
                let (ulo, uhi) = (fu[l][2 * d], fu[l][2 * d + 1]);
                let (plo, phi) = (fpr[l][2 * d], fpr[l][2 * d + 1]);
                // momentum: conservative pressure flux (compressible regime)
                // or the mean of the two face accelerations, each with its
                // face density (balanced force: no spurious acceleration of
                // the light phase at a density jump)
                let dm_cons = -dt * (phi - plo) * inv;
                let dm_bal = rho_new * dt * 0.5 * (fa[l][2 * d] + fa[l][2 * d + 1]);
                v[sys.momentum_index(d)] += wcons * dm_cons + (1.0 - wcons) * dm_bal;
                if let Some(e) = sys.energy_index() {
                    v[e] -= dt * (phi * uhi - plo * ulo) * inv;
                }
                divu += (uhi - ulo) * inv;
            }
            for k in nc.clone() {
                // Kapila part of the volume-fraction source (the α∇·u part
                // went with the advection)
                let kk = sys.nc_coef_k(&wstar[l], k) - wstar[l][k];
                v[k] += dt * kk * divu;
            }
            if ac > 1.0 {
                n_inc += 1;
            }
            // thermodynamic pressure: conservative energy where sound is
            // resolved (C ≤ 1), the elliptic pressure where it is not (C ≥ 2),
            // linear in between. With unresolved sound a conservative update
            // compresses ρe with face values while the pressure equation
            // assumes the cell's own ρc²; at a stiffened-gas interface (ρe ∋
            // γp∞/(γ−1)) the mismatch is O(100) and the projection, which
            // corrects it fully every step, amplifies it.
            if wcons < 1.0 {
                let mut wn = sys.to_prim(&v);
                let e_before = sys.energy_index().map(|e| v[e]).unwrap_or(0.0);
                wn[pi] = wcons * wn[pi] + (1.0 - wcons) * p[l];
                v = sys.to_cons(&wn);
                e_fix += sys.energy_index().map(|e| v[e]).unwrap_or(0.0) - e_before;
            }
            max_ac = max_ac.max(ac);
            self.u[q] = v;
        }
        self.allspeed_p = (0..ncell).map(|l| (p[l], wnext[l])).collect();
        // face-velocity memory, weighted toward the incompressible regime
        let mut delta = vec![[0.0; 3]; ncell];
        for l in 0..ncell {
            let c = coords(l);
            for &d in &dims {
                if let Some(cn) = nb(c, d, 1) {
                    let ln = lin(cn);
                    let wf = (2.0 - cface[l][d]).clamp(0.0, 1.0);
                    if o.regime_blend && wf < 1.0 {
                        let vel = |q: usize| sys.to_prim(&self.u[interior[q]])[vi + d];
                        delta[l][d] = (1.0 - wf) * (unew_f[l][d] - 0.5 * (vel(l) + vel(ln)));
                    }
                }
            }
        }
        self.allspeed_delta = delta;
        let st = &mut self.allspeed_stats;
        st.solves += 1;
        st.cg_iterations += its;
        st.incompressible_cells += n_inc;
        st.cells += ncell;
        st.fallback_faces += fallback;
        st.max_acoustic_cfl = st.max_acoustic_cfl.max(max_ac);
        let vol: f64 = dims.iter().map(|&d| m.dx[d]).product();
        st.energy_change += e_fix * vol;
    }
}

/// Index of the pressure among the primitives (`energy_index`, which holds
/// E in the conserved and p in the primitive vector for every system here).
fn pressure_index<const N: usize, S: System<N>>(sys: &S) -> usize {
    sys.energy_index().expect("the all-speed mode needs a system with an energy equation")
}

fn interior(m: &Mesh) -> Vec<usize> {
    let mut v = Vec::with_capacity(m.cells());
    for i in 0..m.n[0] {
        for j in 0..m.n[1] {
            for k in 0..m.n[2] {
                v.push(m.idx(i + m.g(0), j + m.g(1), k + m.g(2)));
            }
        }
    }
    v
}

/// Jacobi-preconditioned conjugate gradients; returns the iteration count.
fn pcg(apply: &(impl Fn(&[f64], &mut [f64]) + Sync), diag: &[f64], b: &[f64], x: &mut [f64], tol: f64, floor: f64, max_iter: usize) -> usize {
    let n = b.len();
    const CH: usize = 8192; // rayon chunk: small grids stay (mostly) serial
    let dot = |a: &[f64], c: &[f64]| -> f64 { a.par_iter().with_min_len(CH).zip(c.par_iter()).map(|(x, y)| x * y).sum() };
    let mut r = vec![0.0; n];
    apply(x, &mut r);
    r.par_iter_mut().zip(b.par_iter()).for_each(|(ri, bi)| *ri = bi - *ri);
    let bn = dot(b, b).sqrt().max(1e-300);
    let mut z: Vec<f64> = r.iter().zip(diag).map(|(ri, di)| ri / di).collect();
    let mut pdir = z.clone();
    let mut rz = dot(&r, &z);
    let mut ap = vec![0.0; n];
    for it in 0..max_iter {
        if dot(&r, &r).sqrt() <= (tol * bn).max(floor) {
            return it;
        }
        apply(&pdir, &mut ap);
        let alpha = rz / dot(&pdir, &ap);
        x.par_iter_mut().with_min_len(CH).zip(pdir.par_iter()).for_each(|(xi, pi)| *xi += alpha * pi);
        r.par_iter_mut().with_min_len(CH).zip(ap.par_iter()).for_each(|(ri, ai)| *ri -= alpha * ai);
        z.par_iter_mut().with_min_len(CH).zip(r.par_iter().zip(diag.par_iter())).for_each(|(zi, (ri, di))| *zi = ri / di);
        let rz2 = dot(&r, &z);
        let beta = rz2 / rz;
        rz = rz2;
        pdir.par_iter_mut().with_min_len(CH).zip(z.par_iter()).for_each(|(pi, zi)| *pi = zi + beta * *pi);
    }
    max_iter
}
