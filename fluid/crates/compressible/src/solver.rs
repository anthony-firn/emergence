//! Finite-volume solver on a uniform Cartesian mesh (1D, 2D or 3D: axes
//! with one cell are inactive), SSP-RK3 in time.

use crate::scheme::{line_fluxes, FluxMode, SchemeOpts, ShockScheme, Stats};
use crate::system::{System, V};
use rayon::prelude::*;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Bc<const N: usize> {
    Periodic,
    /// zero-gradient extrapolation
    Outflow,
    /// slip wall (mirror, normal momentum reversed)
    Reflect,
    /// fixed conserved state in the ghosts (supersonic inflow, far field)
    Dirichlet(V<N>),
    /// non-reflecting open boundary toward a far-field conserved state
    /// (`boundary::characteristic_ghost`)
    Characteristic(V<N>),
}

#[derive(Clone, Debug)]
pub struct Mesh {
    pub n: [usize; 3],
    pub lo: [f64; 3],
    pub dx: [f64; 3],
    pub ng: usize,
}

impl Mesh {
    pub fn new(n: [usize; 3], lo: [f64; 3], hi: [f64; 3]) -> Self {
        let dx = std::array::from_fn(|d| (hi[d] - lo[d]) / n[d] as f64);
        Mesh { n, lo, dx, ng: 6 }
    }
    /// Use `ng` ghost layers (6 for viscous terms).
    pub fn with_ghosts(mut self, ng: usize) -> Self {
        self.ng = ng;
        self
    }
    #[inline]
    pub fn active(&self, d: usize) -> bool {
        self.n[d] > 1
    }
    /// Allocated extent along d (ghosts only on active axes).
    #[inline]
    pub fn m(&self, d: usize) -> usize {
        if self.active(d) {
            self.n[d] + 2 * self.ng
        } else {
            1
        }
    }
    #[inline]
    pub fn g(&self, d: usize) -> usize {
        if self.active(d) {
            self.ng
        } else {
            0
        }
    }
    #[inline]
    pub(crate) fn periodic_ghost_sources(&self, d: usize, gi: usize) -> (usize, usize) {
        let (n, ng) = (self.n[d], self.ng);
        if n >= ng {
            (n + gi, ng + gi)
        } else {
            // Deep ghost stencils may wrap around a short axis repeatedly.
            (ng + (gi + n - ng % n) % n, ng + gi % n)
        }
    }
    #[inline]
    pub fn idx(&self, i: usize, j: usize, k: usize) -> usize {
        (i * self.m(1) + j) * self.m(2) + k
    }
    pub fn len(&self) -> usize {
        self.m(0) * self.m(1) * self.m(2)
    }
    /// Cell centre of interior cell (i, j, k) (0-based interior indices).
    pub fn x(&self, i: usize, j: usize, k: usize) -> [f64; 3] {
        let c = [i, j, k];
        std::array::from_fn(|d| self.lo[d] + (c[d] as f64 + 0.5) * self.dx[d])
    }
    pub fn cells(&self) -> usize {
        self.n[0] * self.n[1] * self.n[2]
    }
}

struct SyncPtr<const N: usize>(*mut V<N>);
unsafe impl<const N: usize> Sync for SyncPtr<N> {}
unsafe impl<const N: usize> Send for SyncPtr<N> {}
struct SyncF64(*mut f64);
unsafe impl Sync for SyncF64 {}
unsafe impl Send for SyncF64 {}

/// Independent RK cell updates; contiguous rows keep writes disjoint and
/// support any orientation of a 1D/2D/3D mesh. Preserve the scalar formulas.
fn rk_combine<const N: usize, const STAGE: usize>(
    m: &Mesh,
    destination: &mut [V<N>],
    base: &[V<N>],
    rhs: &[V<N>],
    dt: f64,
) {
    let row_len = m.m(2);
    let ny = m.m(1);
    let update = |q: usize, dest: &mut [V<N>], old: &[V<N>], residual: &[V<N>]| {
        let (i, j) = (q / ny, q % ny);
        if i < m.g(0) || i >= m.g(0) + m.n[0] || j < m.g(1) || j >= m.g(1) + m.n[1] {
            return;
        }
        for k in m.g(2)..m.g(2) + m.n[2] {
            for v in 0..N {
                dest[k][v] = match STAGE {
                    1 => old[k][v] + dt * residual[k][v],
                    2 => 0.75 * old[k][v] + 0.25 * (dest[k][v] + dt * residual[k][v]),
                    3 => dest[k][v] / 3.0 + 2.0 / 3.0 * (old[k][v] + dt * residual[k][v]),
                    _ => unreachable!(),
                };
            }
        }
    };
    if m.cells() >= 8192 && (0..3).filter(|&d| m.active(d)).count() >= 2 {
        destination
            .par_chunks_mut(row_len)
            .zip(base.par_chunks(row_len))
            .zip(rhs.par_chunks(row_len))
            .enumerate()
            .for_each(|(q, ((dest, old), rhs))| update(q, dest, old, rhs));
    } else {
        destination
            .chunks_mut(row_len)
            .zip(base.chunks(row_len))
            .zip(rhs.chunks(row_len))
            .enumerate()
            .for_each(|(q, ((dest, old), rhs))| update(q, dest, old, rhs));
    }
}

pub struct Solver<const N: usize, S: System<N>> {
    pub sys: S,
    pub mesh: Mesh,
    /// conserved state, ghost-padded
    pub u: Vec<V<N>>,
    pub bc: [[Bc<N>; 2]; 3],
    pub cfl: f64,
    pub mode: FluxMode,
    pub shock: ShockScheme,
    pub tol: f64,
    pub t: f64,
    pub stats: Stats,
    /// Hu–Adams–Shu positivity-preserving flux limiter (needs CFL ≤ ~0.5)
    pub positivity: bool,
    /// Local relaxation applied to every interior cell after each step
    /// (e.g. liquid–vapour equilibrium); returns true if it changed the cell.
    pub relax: Option<Box<dyn Fn(&S, &mut V<N>, f64) -> bool + Send + Sync>>,
    /// cells changed by `relax` in the last step (the active region)
    pub relaxed_last: usize,
    /// Strang splitting: relax Δt/2 before and after the transport step
    /// (second order in Δt/τ) instead of Δt after it.
    pub strang: bool,
    stage_dt: f64,
    /// surface tension (diffuse-interface multiphase systems)
    pub surface_tension: Option<crate::capillary::SurfaceTension>,
    /// viscosity and heat conduction (None: Euler)
    pub viscous: Option<crate::viscous::Transport>,
    /// Mach-adaptive all-speed time stepping (`allspeed.rs`): time step
    /// from the flow speed, implicit pressure. None: explicit compressible.
    pub allspeed: Option<crate::allspeed::AllSpeedOpts>,
    pub allspeed_stats: crate::allspeed::AllSpeedStats,
    /// all-speed face-velocity memory (+ face of each interior cell, axis)
    pub(crate) allspeed_delta: Vec<[f64; 3]>,
    /// all-speed carried pressure (interior, linear) and its weight against
    /// the EOS pressure (1: EOS) for the next step's advected pressure
    pub(crate) allspeed_p: Vec<(f64, f64)>,
    /// transported scalars (`scalars.rs`)
    pub scalars: Vec<crate::scalars::Scalar>,
    /// per-cell kinetics after each step: (system, conserved state,
    /// conserved scalars, Δt) → changed (phase change, dissolution)
    pub kinetics: Option<Box<dyn Fn(&S, &mut V<N>, &mut [f64], f64) -> bool + Send + Sync>>,
    /// cells changed by `kinetics` in the last step
    pub kinetics_cells: usize,
    /// heat conduction (`thermal.rs`), with phase relaxation per sub-step
    pub conduction: Option<crate::thermal::Conduction>,
    /// conduction sub-steps taken so far
    pub conduction_substeps: usize,
    /// effective carrier face mass fluxes of the current step [carrier][axis]
    pub(crate) flux_acc: Vec<[Vec<f64>; 3]>,
    /// weight of the next `rhs` evaluation in `flux_acc` (0: off)
    flux_w: f64,
    /// carrier densities at the start of the step [carrier][cell]
    pub(crate) carrier_old: Vec<Vec<f64>>,
    rhs_buf: Vec<V<N>>,
    primitive_buf: Vec<V<N>>,
    stage: Vec<V<N>>,
}

impl<const N: usize, S: System<N>> Solver<N, S> {
    /// Build from a primitive initial condition w(x).
    pub fn new(sys: S, mesh: Mesh, init: impl Fn([f64; 3]) -> V<N> + Sync) -> Self {
        let mut u = vec![[0.0; N]; mesh.len()];
        let (gx, gy, gz) = (mesh.g(0), mesh.g(1), mesh.g(2));
        for i in 0..mesh.n[0] {
            for j in 0..mesh.n[1] {
                for k in 0..mesh.n[2] {
                    u[mesh.idx(i + gx, j + gy, k + gz)] = sys.to_cons(&init(mesh.x(i, j, k)));
                }
            }
        }
        let len = mesh.len();
        Solver {
            sys,
            mesh,
            u,
            bc: [[Bc::Periodic; 2]; 3],
            cfl: 0.5,
            mode: FluxMode::Hybrid,
            shock: ShockScheme::Wcns,
            tol: 1e-3,
            t: 0.0,
            stats: Stats::default(),
            viscous: None,
            surface_tension: None,
            positivity: true,
            relax: None,
            relaxed_last: 0,
            strang: false,
            stage_dt: 0.0,
            allspeed: None,
            allspeed_stats: Default::default(),
            allspeed_delta: Vec::new(),
            allspeed_p: Vec::new(),
            scalars: Vec::new(),
            kinetics: None,
            kinetics_cells: 0,
            conduction: None,
            conduction_substeps: 0,
            flux_acc: Vec::new(),
            flux_w: 0.0,
            carrier_old: Vec::new(),
            rhs_buf: vec![[0.0; N]; len],
            primitive_buf: Vec::new(),
            stage: vec![[0.0; N]; len],
        }
    }

    /// Primitive state of interior cell (i, j, k).
    pub fn prim(&self, i: usize, j: usize, k: usize) -> V<N> {
        let m = &self.mesh;
        self.sys.to_prim(&self.u[m.idx(i + m.g(0), j + m.g(1), k + m.g(2))])
    }

    pub(crate) fn fill_ghosts(&self, u: &mut [V<N>]) {
        let bc = self.bc;
        self.fill_ghosts_with(u, &bc);
    }

    /// Ghost cells with the given boundary conditions.
    #[inline]
    pub(crate) fn fill_ghosts_with(&self, u: &mut [V<N>], bc: &[[Bc<N>; 2]; 3]) {
        let m = &self.mesh;
        let deep = (0..3).any(|d| m.active(d) && m.n[d] < m.ng
            && bc[d].iter().any(|b| *b == Bc::Reflect));
        if deep {
            self.fill_short_reflecting_ghosts(u, bc);
        } else {
            self.fill_ghosts_impl::<false>(u, bc);
        }
    }

    #[cold]
    #[inline(never)]
    fn fill_short_reflecting_ghosts(&self, u: &mut [V<N>], bc: &[[Bc<N>; 2]; 3]) {
        self.fill_ghosts_impl::<true>(u, bc);
    }

    // Sharing the implementation keeps both boundary paths consistent, while
    // the constant flag removes deep-reflection calls from the common kernel.
    fn fill_ghosts_impl<const DEEP: bool>(&self, u: &mut [V<N>], bc: &[[Bc<N>; 2]; 3]) {
        let m = &self.mesh;
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
                    let characteristic = [
                        if let Bc::Characteristic(far) = bc[d][0] {
                            Some(crate::boundary::CharacteristicStencil::new(&self.sys,
                                &u[at(ng)], &u[at(ng + 1)], &far, d, 0))
                        } else { None },
                        if let Bc::Characteristic(far) = bc[d][1] {
                            Some(crate::boundary::CharacteristicStencil::new(&self.sys,
                                &u[at(n + ng - 1)], &u[at(n + ng - 2)], &far, d, 1))
                        } else { None },
                    ];
                    for gi in 0..ng {
                        // low side ghost gi (cell index gi), high side n+ng+gi
                        let (lo_dst, hi_dst) = (at(gi), at(n + ng + gi));
                        let (lo_src, hi_src) = m.periodic_ghost_sources(d, gi);
                        let (lo_src, hi_src) = (at(lo_src), at(hi_src));
                        let (lo_ref, hi_ref) = (at(2 * ng - 1 - gi), at(n + ng - 1 - gi));
                        u[lo_dst] = match bc[d][0] {
                            Bc::Periodic => u[lo_src],
                            Bc::Outflow => u[at(ng)],
                            Bc::Dirichlet(v) => v,
                            Bc::Characteristic(_) => characteristic[0].as_ref().unwrap().conserved(&self.sys, (ng - gi) as f64),
                            Bc::Reflect => {
                                if DEEP && n < ng {
                                    crate::boundary::reflected_ghost(&self.sys, u, &at,
                                        gi as isize - ng as isize, n, ng, d, &bc[d], &characteristic)
                                } else {
                                    let mut v = u[lo_ref];
                                    v[self.sys.momentum_index(d)] *= -1.0;
                                    v
                                }
                            }
                        };
                        u[hi_dst] = match bc[d][1] {
                            Bc::Periodic => u[hi_src],
                            Bc::Outflow => u[at(n + ng - 1)],
                            Bc::Dirichlet(v) => v,
                            Bc::Characteristic(_) => characteristic[1].as_ref().unwrap().conserved(&self.sys, (gi + 1) as f64),
                            Bc::Reflect => {
                                if DEEP && n < ng {
                                    crate::boundary::reflected_ghost(&self.sys, u, &at,
                                        (n + gi) as isize, n, ng, d, &bc[d], &characteristic)
                                } else {
                                    let mut v = u[hi_ref];
                                    v[self.sys.momentum_index(d)] *= -1.0;
                                    v
                                }
                            }
                        };
                    }
                }
            }
        }
    }

    /// Spatial residual −∇·F (interior cells), ghosts of `u` refreshed.
    fn rhs(&mut self, which: Which) -> Stats {
        let mut u = match which {
            Which::U => std::mem::take(&mut self.u),
            Which::Stage => std::mem::take(&mut self.stage),
        };
        self.fill_ghosts(&mut u);
        let mut out = std::mem::take(&mut self.rhs_buf);
        out.par_iter_mut().for_each(|v| *v = [0.0; N]);
        let m = self.mesh.clone();
        let nv = N;
        // Reuse the same stage's primitive conversion along every direction.
        // Skip edge/corner ghosts, which no interior pencil needs and whose
        // states a custom EOS is not required to accept. Cache only when it
        // actually saves conversions relative to the pencils.
        let pencil_work: usize = (0..3).filter(|&d|m.active(d)).map(|d|{
            m.m(d)*(0..3).filter(|&e|e!=d).map(|e|m.n[e]).product::<usize>()
        }).sum();
        let cached_primitives = m.len() < pencil_work;
        let mut primitive = std::mem::take(&mut self.primitive_buf);
        if cached_primitives {
            primitive.resize(m.len(), [0.0; N]);
            let plane = m.m(1)*m.m(2);
            primitive.par_chunks_mut(plane).enumerate().for_each(|(i,row)| {
                let ox = i<m.g(0) || i>=m.g(0)+m.n[0];
                for j in 0..m.m(1) {
                    let oy = j<m.g(1) || j>=m.g(1)+m.n[1];
                    if ox && oy {continue;}
                    let kr = if ox || oy {m.g(2)..m.g(2)+m.n[2]}else{0..m.m(2)};
                    for k in kr {let q=j*m.m(2)+k;row[q]=self.sys.to_prim(&u[i*plane+q]);}
                }
            });
        }
        let mut total = Stats::default();
        for d in 0..3 {
            if !m.active(d) {
                continue;
            }
            let others: Vec<usize> = (0..3).filter(|&e| e != d).collect();
            let (e0, e1) = (others[0], others[1]);
            let pencils: Vec<(usize, usize)> = (0..m.n[e0]).flat_map(|a| (0..m.n[e1]).map(move |b| (a, b))).collect();
            let ptr = SyncPtr(out.as_mut_ptr());
            let ndim = (0..3).filter(|&e| m.active(e)).count() as f64;
            let opts = SchemeOpts { mode: self.mode, shock: self.shock, tol: self.tol, lambda: if self.positivity { ndim * self.stage_dt / m.dx[d] } else { 0.0 } };
            let (sys, uref) = (&self.sys, &u);
            let nc = sys.nc_vars();
            // carrier mass-flux accumulators of the scalars (this axis)
            let fw = self.flux_w;
            let carriers = if fw != 0.0 { self.carriers() } else { Vec::new() };
            let cidx: Vec<std::ops::Range<usize>> = carriers.iter().map(|c| crate::scalars::carrier_indices(sys, *c)).collect();
            let accp: Vec<SyncF64> = (0..carriers.len()).map(|ci| SyncF64(self.flux_acc[ci][d].as_mut_ptr())).collect();
            let st = pencils
                .par_iter()
                .map(|&(a, b)| {
                    let _ = &ptr;
                    let _ = &accp;
                    let at = |c: usize| {
                        let mut ijk = [0; 3];
                        ijk[d] = c;
                        ijk[e0] = a + m.g(e0);
                        ijk[e1] = b + m.g(e1);
                        m.idx(ijk[0], ijk[1], ijk[2])
                    };
                    let len = m.m(d);
                    let w: Vec<V<N>> = (0..len).map(|c| if cached_primitives { primitive[at(c)] } else { sys.to_prim(&uref[at(c)]) }).collect();
                    let mut fl = vec![[0.0; N]; m.n[d] + 1];
                    let mut uf = vec![0.0; m.n[d] + 1];
                    let mut wf = vec![[0.0; N]; m.n[d] + 1];
                    let mut s = Stats::default();
                    line_fluxes(sys, &w, d, m.ng, opts, &mut fl, &mut uf, &mut wf, &mut s);
                    if fw != 0.0 {
                        let pidx = a * m.n[e1] + b;
                        for (ci, r) in cidx.iter().enumerate() {
                            for f in 0..=m.n[d] {
                                let fm: f64 = r.clone().map(|k| fl[f][k]).sum();
                                // SAFETY: pencils own disjoint face ranges
                                unsafe { *accp[ci].0.add(pidx * (m.n[d] + 1) + f) += fw * fm };
                            }
                        }
                    }
                    let ncg = sys.nc_general();
                    let inv = 1.0 / m.dx[d];
                    for c in 0..m.n[d] {
                        let idx = at(c + m.ng);
                        // SAFETY: each pencil owns its cells; pencils of one
                        // direction are disjoint and directions run in turn.
                        let o = unsafe { &mut *ptr.0.add(idx) };
                        for k in 0..nv {
                            o[k] -= (fl[c + 1][k] - fl[c][k]) * inv;
                        }
                        for q in nc.clone() {
                            o[q] += sys.nc_coef_k(&w[c + m.ng], q) * (uf[c + 1] - uf[c]) * inv;
                        }
                        if ncg {
                            let mut dw = [0.0; N];
                            for k in 0..N {
                                dw[k] = (wf[c + 1][k] - wf[c][k]) * inv;
                            }
                            let add = sys.nc_rhs(&w[c + m.ng], &dw, d);
                            for k in 0..N {
                                o[k] += add[k];
                            }
                        }
                    }
                    s
                })
                .reduce(Stats::default, |mut a, b| {
                    a.add(&b);
                    a
                });
            total.add(&st);
        }
        if let Some(tr) = &self.viscous {
            crate::viscous::add_viscous(&self.sys, &m, tr, &u, &mut out);
        }
        if let Some(st) = &self.surface_tension {
            st.add(&self.sys, &m, &u, &mut out);
        }
        self.rhs_buf = out;
        self.primitive_buf = primitive;
        match which {
            Which::U => self.u = u,
            Which::Stage => self.stage = u,
        }
        total
    }

    /// Stable time step from the CFL condition.
    pub fn dt_cfl(&self) -> f64 {
        if let Some(o) = &self.allspeed {
            return self.allspeed_dt(o);
        }
        let m = &self.mesh;
        let sys = &self.sys;
        let need_source = !sys.nc_vars().is_empty();
        let (rate, source, diffusion) = (0..m.n[0])
            .into_par_iter()
            .map(|i| {
                let (mut rate, mut source, mut diffusion) = (0.0f64, 0.0f64, crate::viscous::DiffusionBounds::default());
                for j in 0..m.n[1] {
                    for k in 0..m.n[2] {
                        let c = [i + m.g(0), j + m.g(1), k + m.g(2)];
                        let w = sys.to_prim(&self.u[m.idx(c[0], c[1], c[2])]);
                        let acoustic: f64 = (0..3)
                            .filter(|&d| m.active(d))
                            .map(|d| sys.max_speed(&w, d) / m.dx[d])
                            .sum();
                        rate = rate.max(acoustic);
                        if need_source {
                            let mut divu = 0.0;
                            for d in (0..3).filter(|&d| m.active(d)) {
                                let (va, vb) = self.neighbor_velocities(c, d);
                                divu += (vb - va).abs() / (2.0 * m.dx[d]);
                            }
                            source = source.max(sys.dilatation_rate(&w, divu));
                        }
                        if let Some(tr) = &self.viscous {
                            diffusion.observe(sys, tr, &w);
                        }
                    }
                }
                (rate, source, diffusion)
            })
            .reduce(
                || (0.0, 0.0, crate::viscous::DiffusionBounds::default()),
                |a, b| (a.0.max(b.0), a.1.max(b.1), a.2.merge(b.2)),
            );
        let mut dt = self.cfl / rate;
        if source > 0.0 {
            dt = dt.min(0.5 / source);
        }
        dt = dt.min(self.viscous_dt(diffusion));
        dt
    }

    /// Include the three ghost flux nodes on each boundary pencil using the
    /// same current-state characteristic closure as the viscous RHS. Periodic,
    /// outflow and reflecting ghosts inherit interior transport coefficients.
    pub(crate) fn viscous_dt(&self, mut bounds: crate::viscous::DiffusionBounds) -> f64 {
        let Some(tr) = &self.viscous else { return f64::INFINITY; };
        let (m, sys) = (&self.mesh, &self.sys);
        let mut inv = 0.0;
        for d in (0..3).filter(|&d| m.active(d)) {
            inv += 1.0 / (m.dx[d] * m.dx[d]);
            for side in 0..2 {
                match self.bc[d][side] {
                    Bc::Dirichlet(far) => bounds.observe_transport(sys, tr, &sys.to_prim(&far)),
                    Bc::Characteristic(far) => {
                        let others: Vec<_> = (0..3).filter(|&e| e != d).collect();
                        for a in 0..m.n[others[0]] {
                            for b in 0..m.n[others[1]] {
                                let mut c = [0; 3];
                                c[others[0]] = a + m.g(others[0]);
                                c[others[1]] = b + m.g(others[1]);
                                c[d] = m.g(d) + if side == 0 { 0 } else { m.n[d] - 1 };
                                let center = self.u[m.idx(c[0], c[1], c[2])];
                                if side == 0 { c[d] += 1; } else { c[d] -= 1; }
                                let neighbor = self.u[m.idx(c[0], c[1], c[2])];
                                let stencil = crate::boundary::CharacteristicStencil::new(
                                    sys, &center, &neighbor, &far, d, side);
                                for layer in 1..=3 {
                                    let ghost = stencil.conserved(sys, layer as f64);
                                    bounds.observe_transport(sys, tr, &sys.to_prim(&ghost));
                                }
                            }
                        }
                    }
                    _ => {}
                }
            }
        }
        0.3 / (bounds.diffusivity(tr) * inv)
    }

    #[inline]
    pub(crate) fn neighbor_velocities(&self, c: [usize; 3], d: usize) -> (f64, f64) {
        let m = &self.mesh;
        if c[d] > m.g(d) && c[d] < m.g(d) + m.n[d] - 1 {
            let q = m.idx(c[0], c[1], c[2]);
            let stride = match d {
                0 => m.m(1) * m.m(2),
                1 => m.m(2),
                _ => 1,
            };
            (
                self.sys.velocity_conserved(&self.u[q - stride], d),
                self.sys.velocity_conserved(&self.u[q + stride], d),
            )
        } else {
            (
                self.sys.velocity_conserved(&self.neighbor_conserved(c, d, false), d),
                self.sys.velocity_conserved(&self.neighbor_conserved(c, d, true), d),
            )
        }
    }

    // CFL estimates must not read stale (or never initialized) ghost cells.
    // Evaluate the nearest boundary neighbor without copying the whole grid.
    pub(crate) fn neighbor_prim(&self, c: [usize; 3], d: usize, high: bool) -> V<N> {
        self.sys.to_prim(&self.neighbor_conserved(c, d, high))
    }

    pub(crate) fn neighbor_conserved(&self, c: [usize; 3], d: usize, high: bool) -> V<N> {
        let m = &self.mesh;
        let boundary = if high {
            c[d] == m.g(d) + m.n[d] - 1
        } else {
            c[d] == m.g(d)
        };
        let mut q = c;
        let u = if !boundary {
            if high {
                q[d] += 1;
            } else {
                q[d] -= 1;
            }
            self.u[m.idx(q[0], q[1], q[2])]
        } else {
            let center = self.u[m.idx(c[0], c[1], c[2])];
            match self.bc[d][high as usize] {
                Bc::Periodic => {
                    q[d] = if high { m.g(d) } else { m.g(d) + m.n[d] - 1 };
                    self.u[m.idx(q[0], q[1], q[2])]
                }
                Bc::Outflow => center,
                Bc::Dirichlet(u) => u,
                Bc::Reflect => {
                    let mut u = center;
                    u[self.sys.momentum_index(d)] *= -1.0;
                    u
                }
                Bc::Characteristic(far) => {
                    if high {
                        q[d] -= 1;
                    } else {
                        q[d] += 1;
                    }
                    crate::boundary::characteristic_ghost_layer(
                        &self.sys,
                        &center,
                        &self.u[m.idx(q[0], q[1], q[2])],
                        &far,
                        d,
                        high as usize,
                        1.0,
                    )
                }
            }
        };
        u
    }

    /// One SSP-RK3 step of size dt (Shu–Osher form).
    fn apply_relax(&mut self, interior: &[usize], dt: f64) {
        if let Some(rx) = &self.relax {
            let sys = &self.sys;
            let ptr = SyncPtr(self.u.as_mut_ptr());
            self.relaxed_last = interior
                .par_iter()
                .filter(|&&q| {
                    let _ = &ptr;
                    // SAFETY: each interior index is visited once
                    let c = unsafe { &mut *ptr.0.add(q) };
                    rx(sys, c, dt)
                })
                .count();
        }
    }

    pub fn step(&mut self, dt: f64) {
        self.stage_dt = dt;
        let interior = self.interior_indices();
        if self.strang {
            self.apply_relax(&interior, 0.5 * dt);
        }
        if let Some(o) = self.allspeed {
            // global switch: the all-speed step only where it pays (Δt at
            // least `min_gain` × the explicit one, i.e. a low-Mach phase)
            self.allspeed = None;
            let dt_explicit = self.dt_cfl();
            self.allspeed = Some(o);
            if dt <= o.min_gain * dt_explicit {
                let n = self.explicit_substeps(dt);
                self.allspeed_stats.explicit_steps += 1;
                self.allspeed_stats.explicit_substeps += n;
                self.allspeed_delta.clear();
                self.allspeed_p.clear();
                self.allspeed_stats.last_accel = 0.0;
                self.t += dt;
                return;
            }
            let backup = self.u.clone();
            self.allspeed_step(dt, &o);
            // accept only states the explicit scheme can continue from
            // (admissible, and its Δt not collapsed — e.g. near-vacuum gas
            // left by the elliptic pressure at a strong shock)
            let ok = interior.par_iter().all(|&q| self.sys.admissible(&self.sys.to_prim(&self.u[q]))) && {
                self.allspeed = None;
                let dte_new = self.dt_cfl();
                self.allspeed = Some(o);
                dte_new > 0.1 * dt_explicit
            };
            if !ok {
                // the low-Mach scheme lost admissibility (a strong shock):
                // redo the step with explicit compressible sub-steps
                self.u = backup;
                let n = self.explicit_substeps(dt);
                self.allspeed_stats.explicit_fallbacks += 1;
                self.allspeed_stats.explicit_substeps += n;
                self.allspeed_delta.clear();
                self.allspeed_p.clear();
                self.allspeed_stats.last_accel = 0.0;
            } else {
                self.apply_relax(&interior, if self.strang { 0.5 * dt } else { dt });
            }
            self.t += dt;
            return;
        }
        let phi_old = if self.scalars.is_empty() { Vec::new() } else { self.scalar_phis(&self.u.clone()) };
        if !self.scalars.is_empty() {
            self.reset_flux_acc();
        }
        let sc = !self.scalars.is_empty();
        // stage 1: u1 = u + dt L(u)
        self.flux_w = if sc { 1.0 / 6.0 } else { 0.0 };
        let s1 = self.rhs(Which::U);
        self.stage.clone_from(&self.u);
        rk_combine::<N, 1>(&self.mesh, &mut self.stage, &self.u, &self.rhs_buf, dt);
        // stage 2: u2 = 3/4 u + 1/4 (u1 + dt L(u1))
        self.flux_w = if sc { 1.0 / 6.0 } else { 0.0 };
        let s2 = self.rhs(Which::Stage);
        rk_combine::<N, 2>(&self.mesh, &mut self.stage, &self.u, &self.rhs_buf, dt);
        // stage 3: u = 1/3 u + 2/3 (u2 + dt L(u2))
        self.flux_w = if sc { 2.0 / 3.0 } else { 0.0 };
        let s3 = self.rhs(Which::Stage);
        self.flux_w = 0.0;
        rk_combine::<N, 3>(&self.mesh, &mut self.u, &self.stage, &self.rhs_buf, dt);
        for s in [s1, s2, s3] {
            self.stats.add(&s);
        }
        self.finish_scalars(&phi_old, dt);
        if self.conduction.is_some() {
            let mut u = std::mem::take(&mut self.u);
            self.conduct(&mut u, dt);
            self.u = u;
        }
        self.apply_relax(&interior, if self.strang { 0.5 * dt } else { dt });
        self.t += dt;
    }

    /// Cover `dt` with explicit compressible steps (all-speed mode off);
    /// `self.t` is left unchanged. Returns the number of sub-steps.
    fn explicit_substeps(&mut self, dt: f64) -> usize {
        let o = self.allspeed.take();
        let (t0, t_end) = (self.t, self.t + dt);
        let mut n = 0;
        while self.t < t_end - 1e-12 * dt {
            let dte = self.dt_cfl().min(t_end - self.t);
            assert!(n < 10_000 && dte > 1e-9 * dt, "explicit sub-steps collapsed (Δt {dte:e} of {dt:e}) at t = {}", self.t);
            self.step(dte);
            n += 1;
        }
        self.t = t0;
        self.allspeed = o;
        n
    }

    /// Advance to time `t_end` with CFL steps.
    pub fn run_to(&mut self, t_end: f64) -> usize {
        let mut steps = 0;
        // stop within 1e-10 relative (summed steps leave rounding residues)
        while self.t < t_end - 1e-10 * t_end.abs() {
            let dt = self.dt_cfl().min(t_end - self.t);
            assert!(dt.is_finite() && dt > 1e-12 * t_end.abs().max(1e-300), "time step collapsed (dt = {dt:e} at t = {}): the solution lost admissibility", self.t);
            self.step(dt);
            steps += 1;
        }
        steps
    }

    fn interior_indices(&self) -> Vec<usize> {
        let m = &self.mesh;
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

    /// Totals of the conserved variables over the interior (× cell volume).
    pub fn totals(&self) -> V<N> {
        let m = &self.mesh;
        let vol: f64 = (0..3).filter(|&d| m.active(d)).map(|d| m.dx[d]).product();
        let mut s = [0.0; N];
        for q in self.interior_indices() {
            for k in 0..N {
                s[k] += self.u[q][k] * vol;
            }
        }
        s
    }
}

#[derive(Clone, Copy)]
enum Which {
    U,
    Stage,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{eos::Nasg, Euler};

    #[test]
    fn deep_periodic_ghosts_match_the_torus_including_corners() {
        for n in [[2, 3, 4], [8, 2, 1]] {
            for ng in [6, 10] {
                let s = Solver::new(
                    Euler {
                        eos: Nasg::ideal(1.4),
                    },
                    Mesh::new(n, [0.0; 3], [1.0; 3]).with_ghosts(ng),
                    |x| [1.0 + x[0] + 2.0 * x[1] + 3.0 * x[2], 0.3, -0.2, 0.1, 1.0],
                );
                let m = &s.mesh;
                let mut u = vec![[f64::NAN; 5]; m.len()];
                for i in 0..n[0] {
                    for j in 0..n[1] {
                        for k in 0..n[2] {
                            let q = m.idx(i + m.g(0), j + m.g(1), k + m.g(2));
                            u[q] = s.u[q];
                        }
                    }
                }
                s.fill_ghosts(&mut u);
                for i in 0..m.m(0) {
                    for j in 0..m.m(1) {
                        for k in 0..m.m(2) {
                            let c = [i, j, k];
                            let wrapped: [usize; 3] = std::array::from_fn(|d| {
                                m.g(d)
                                    + (c[d] as isize - m.g(d) as isize).rem_euclid(n[d] as isize)
                                        as usize
                            });
                            assert_eq!(
                                u[m.idx(i, j, k)],
                                s.u[m.idx(wrapped[0], wrapped[1], wrapped[2])]
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn short_reflecting_ghosts_match_mirrored_images_including_corners() {
        for n in [[2, 3, 4], [8, 2, 1]] {
            for ng in [6, 10] {
                let mut s = Solver::new(Euler { eos: Nasg::ideal(1.4) },
                    Mesh::new(n, [0.0; 3], [1.0; 3]).with_ghosts(ng),
                    |x| [1.0 + x[0] + 2.0*x[1] + 3.0*x[2], 0.3, -0.2, 0.1, 1.0]);
                s.bc = [[Bc::Reflect; 2]; 3];
                let m = &s.mesh;
                let mut u = vec![[f64::NAN; 5]; m.len()];
                for i in 0..n[0] { for j in 0..n[1] { for k in 0..n[2] {
                    let q = m.idx(i + m.g(0), j + m.g(1), k + m.g(2));
                    u[q] = s.u[q];
                }}}
                s.fill_ghosts(&mut u);
                for i in 0..m.m(0) { for j in 0..m.m(1) { for k in 0..m.m(2) {
                    let c = [i,j,k];
                    let image: [usize; 3] = std::array::from_fn(|d| {
                        let r = (c[d] as isize - m.g(d) as isize).rem_euclid(2*n[d] as isize) as usize;
                        m.g(d) + if r < n[d] { r } else { 2*n[d] - 1 - r }
                    });
                    let mut expected = s.u[m.idx(image[0], image[1], image[2])];
                    for d in 0..3 {
                        let r = (c[d] as isize - m.g(d) as isize).rem_euclid(2*n[d] as isize) as usize;
                        if r >= n[d] { expected[1+d] *= -1.0; }
                    }
                    assert_eq!(u[m.idx(i,j,k)], expected);
                }}}
            }
        }
    }

    #[test]
    fn deep_reflections_use_the_opposite_current_boundary_condition() {
        for d in 0..3 {
            for reflecting_side in 0..2 {
                for opposite_kind in 0..4 {
                    let mut n = [1;3]; n[d]=2;
                    let sys = Euler { eos: Nasg::ideal(1.4) };
                    let mut s = Solver::new(sys, Mesh::new(n, [0.0;3], [1.0;3]), |x| {
                        let mut w = [1.0 + x[d], 0.0, 0.0, 0.0, 1.0 + x[d]];
                        w[1+d] = 3.0;
                        w
                    });
                    let mut far = [2.0, 0.0, 0.0, 0.0, 2.0]; far[1+d] = 3.0;
                    let far = sys.to_cons(&far);
                    let opposite = 1-reflecting_side;
                    let bc = [Bc::Outflow, Bc::Dirichlet(far), Bc::Periodic, Bc::Characteristic(far)][opposite_kind];
                    s.bc[d][reflecting_side] = Bc::Reflect;
                    s.bc[d][opposite] = bc;
                    let m = &s.mesh;
                    let at = |x: usize| { let mut c = [0;3]; c[d]=x; m.idx(c[0],c[1],c[2]) };
                    let mut u = vec![[f64::NAN;5];m.len()];
                    for i in 0..2 {u[at(m.ng+i)] = s.u[at(m.ng+i)];}
                    s.fill_ghosts(&mut u);
                    for layer in 1..=m.ng {
                        let image = if reflecting_side == 0 { layer as isize-1 } else { n[d] as isize-layer as isize };
                        let mut expected = if image >= 0 && image < 2 {
                            s.u[at(m.ng+image as usize)]
                        } else {
                            match bc {
                                Bc::Outflow => s.u[at(m.ng + if opposite==0 {0} else {1})],
                                Bc::Periodic => s.u[at(m.ng + image.rem_euclid(2) as usize)],
                                Bc::Dirichlet(far) => far,
                                Bc::Characteristic(far) => {
                                    let center = s.u[at(m.ng + if opposite==0 {0} else {1})];
                                    let neighbor = s.u[at(m.ng + if opposite==0 {1} else {0})];
                                    let depth = if opposite==0 {-image} else {image-2+1};
                                    crate::boundary::characteristic_ghost_layer(&sys,&center,&neighbor,&far,d,opposite,depth as f64)
                                }
                                _ => unreachable!(),
                            }
                        };
                        expected[1+d] *= -1.0;
                        let q = at(if reflecting_side==0 {m.ng-layer} else {m.ng+2+layer-1});
                        assert_eq!(u[q],expected,"axis {d}, side {reflecting_side}, opposite {bc:?}, layer {layer}");
                    }
                }
            }
        }
    }
}
