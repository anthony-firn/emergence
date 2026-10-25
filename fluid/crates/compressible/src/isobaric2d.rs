//! Slow phase change in two dimensions: the Eulerian (fixed-grid) form of the
//! isobaric enthalpy model, **with the expansion flow** driven by the density
//! change on melting/freezing.
//!
//! The one-dimensional model in [`crate::isobaric`] lives in mass coordinates:
//! a cell keeps its mass and its width follows the density of the phase it
//! holds, so expansion on melting moves material exactly. In two dimensions
//! that Lagrangian grid would tangle, and the natural rung is a fixed uniform
//! finite-volume grid carrying the **volumetric enthalpy** `U = ρ h` and the
//! **mass density** `ρ`.
//!
//! Field equations (fixed pressure `p₀`):
//! ```text
//!   ∂U/∂t + ∇·(U u) = ∇·(k ∇T) + q,      U = ρ h,
//!   ∂ρ/∂t + ∇·(ρ u) = 0,
//! ```
//! with `ρ` and the phase mass fractions `y_k` following the NASG equation of
//! state at `p₀` and the phase equilibrium of [`PhaseSet`] (`Ladder` on the
//! Gibbs energies). The temperature and the phase state are recovered from the
//! **specific** enthalpy `h = U/ρ` by an exact inversion (the
//! enthalpy–temperature iteration): within one phase `h_k(T)` is monotone in
//! `T` (Newton/bisection), and at a transition `T = T_tr` while the two-phase
//! mixture moves along the tie line, so `h(φ) = (1−φ) h_a + φ h_b` gives the
//! fraction `φ` in closed form. Latent heat therefore costs no regularisation.
//!
//! **Expansion flow.** At constant pressure the density is a function of the
//! local equilibrium state, so mass conservation with `ρ = ρ(h, y)` requires a
//! dilatation `∇·u = −(1/ρ) Dρ/Dt`. In the slow limit the flow has no vorticity
//! source, so it is written as the gradient of a potential, `u = ∇ψ`, with the
//! variable-coefficient Poisson problem
//! ```text
//!   ∇·(ρ ∇ψ) = −Dρ/Dt = (ρ_i − ρ_i^eq)/Δt  =: s,
//! ```
//! solved each sub-step. The right-hand side is the density change produced by
//! conduction at fixed material point: after a conduction update the EOS gives a
//! provisional density `ρ^eq = 1/v(h)`, and the difference from the carried `ρ`
//! is the volume source. Solving with the *same* face coefficients used to
//! update `ρ` makes the projection exact up to round-off (the cell ends the
//! step at `ρ^eq`), and the same conservative mass flux advects `U`, so
//! `∫ρ dV` and `∫U dV` are both conserved by telescoping. A fixed wall is
//! `∂ψ/∂n = 0` (`u = 0`); an open, pressure-fixed side is `ψ = 0` (the normal
//! velocity is free, like the free end of the 1D Lagrangian mesh). With no open
//! side the problem is pure Neumann/periodic and the source is made mean-zero
//! (constant-volume enclosure). `expansion = false` restores the old zero-flow
//! closure, for measuring the error it removes.
//!
//! The Poisson problem is solved matrix-free by conjugate gradients with a
//! Jacobi preconditioner, warm-started from the previous sub-step; 1D grids with
//! an open end use a direct tridiagonal (Thomas) solve. A relative density
//! change below `1e-6` per step produces a flow far below the conduction
//! velocity scale and applies the zero-flow projection instead, so the
//! equal-density verification material (thermal expansion only,
//! `Δρ/ρ ~ 10⁻⁷` per step) runs at the old speed.
//!
//! Conduction is explicit and stable, limited by each cell's *sensible*
//! volumetric heat capacity `ρ min_k c_{p,k}` against its face conductances, as
//! in [`crate::isobaric`]. The expansion flow adds a CFL limit
//! `Δt < 0.4 Δx / max|u|`.
//!
//! For a melt pool the next unfold is convection. [`rayleigh_benard`] and
//! [`rayleigh_internal`] are the cheap indicators (the critical values
//! `1707.762` and `2772` are those used by `photothermal::indicators`), and
//! [`Isobaric2D::buoyancy`] returns the Boussinesq acceleration
//! `g β (T − T_ref)` that a flow solver would take as its momentum source.

use crate::eos::Nasg;
use crate::isobaric::{Fr, Ladder, PhaseSet, MAXP};
use crate::phase_change::v_pt;

/// Boundary condition on one side of the box.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Bc2D {
    /// Fixed temperature (the face sits half a cell outside the last centre).
    Dirichlet(f64),
    /// Zero normal heat flux.
    Adiabatic,
}

/// The four sides and the periodic pairings.
#[derive(Clone, Copy, Debug)]
pub struct Boundaries2D {
    pub west: Bc2D,
    pub east: Bc2D,
    pub south: Bc2D,
    pub north: Bc2D,
    /// Wrap west↔east (west/east are then ignored).
    pub periodic_x: bool,
    /// Wrap south↔north (south/north are then ignored).
    pub periodic_y: bool,
}

impl Default for Boundaries2D {
    fn default() -> Self {
        Boundaries2D {
            west: Bc2D::Adiabatic,
            east: Bc2D::Adiabatic,
            south: Bc2D::Adiabatic,
            north: Bc2D::Adiabatic,
            periodic_x: false,
            periodic_y: false,
        }
    }
}

/// Flow boundary condition for the expansion potential on one side.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FlowBc {
    /// Fixed wall: `∂ψ/∂n = 0`, no normal flow (`u = 0`).
    Wall,
    /// Open, pressure-fixed side: `ψ = 0`, normal velocity free (outflow or
    /// inflow), like the free end of the 1D Lagrangian mesh.
    Open,
}

/// Expansion-flow boundary conditions. Periodic directions (set on
/// [`Boundaries2D`]) override these. Default: a closed enclosure (all walls).
#[derive(Clone, Copy, Debug)]
pub struct ExpansionBc {
    pub west: FlowBc,
    pub east: FlowBc,
    pub south: FlowBc,
    pub north: FlowBc,
}

impl Default for ExpansionBc {
    fn default() -> Self {
        ExpansionBc {
            west: FlowBc::Wall,
            east: FlowBc::Wall,
            south: FlowBc::Wall,
            north: FlowBc::Wall,
        }
    }
}

/// Two-dimensional isobaric phase-change grid.
pub struct Isobaric2D {
    pub model: PhaseSet,
    pub ladder: Ladder,
    /// thermal conductivity of each phase, W/(m K), in the phase-set order
    pub k: Vec<f64>,
    pub p: f64,
    pub nx: usize,
    pub ny: usize,
    pub dx: f64,
    pub dy: f64,
    /// conserved volumetric enthalpy `ρ h`, J/m³
    pub u: Vec<f64>,
    /// conserved mass density, kg/m³ (carried by the expansion flow)
    pub rho: Vec<f64>,
    /// equilibrium phase mass fractions, row-major `[j*nx + i]`
    pub y: Vec<Fr>,
    /// volumetric heat source per cell, W/m³
    pub source: Vec<f64>,
    pub bc: Boundaries2D,
    /// Optional per-face ghost temperatures overriding a `Dirichlet` side:
    /// `[west(ny), east(ny), south(nx), north(nx)]`. Used for manufactured
    /// solutions whose boundary data varies along the side.
    pub bc_row: [Option<Vec<f64>>; 4],
    /// physical time
    pub t: f64,
    /// temperature range over which the ladder was built (inversion bounds)
    pub t_lo: f64,
    pub t_hi: f64,
    /// solve the expansion flow (`true`) or use the zero-flow closure (`false`)
    pub expansion: bool,
    /// expansion-flow side conditions (ignored in periodic directions)
    pub flow_bc: ExpansionBc,
    /// cached expansion potential, used as the next solve's initial guess
    pub psi: Vec<f64>,
}

impl Isobaric2D {
    /// Uniform grid of `nx × ny` cells of size `dx × dy`, initially in pure
    /// phase `ph` at temperature `t0` and pressure `p`.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        model: PhaseSet,
        p: f64,
        k: Vec<f64>,
        nx: usize,
        ny: usize,
        dx: f64,
        dy: f64,
        ph: usize,
        t0: f64,
        t_range: (f64, f64),
    ) -> Self {
        assert!(ph < model.len());
        let ladder = Ladder::new(&model, p, t_range.0, t_range.1);
        let n = nx * ny;
        let v0 = v_pt(&model.eos[ph], p, t0);
        let u0 = model.enthalpy(ph, p, t0) / v0;
        let mut y = [0.0; MAXP];
        y[ph] = 1.0;
        Isobaric2D {
            model,
            ladder,
            k,
            p,
            nx,
            ny,
            dx,
            dy,
            u: vec![u0; n],
            rho: vec![1.0 / v0; n],
            y: vec![y; n],
            source: vec![0.0; n],
            bc: Boundaries2D::default(),
            bc_row: [None, None, None, None],
            t: 0.0,
            t_lo: t_range.0,
            t_hi: t_range.1,
            expansion: true,
            flow_bc: ExpansionBc::default(),
            psi: vec![0.0; n],
        }
    }

    #[inline]
    pub fn idx(&self, i: usize, j: usize) -> usize {
        j * self.nx + i
    }

    pub fn len(&self) -> usize {
        self.u.len()
    }
    pub fn is_empty(&self) -> bool {
        self.u.is_empty()
    }
    pub fn eos(&self) -> &[Nasg] {
        &self.model.eos
    }
    fn np(&self) -> usize {
        self.model.len()
    }

    /// Specific enthalpy `h_k(T)` of a pure phase at temperature `t`.
    pub fn h_pure(&self, k: usize, t: f64) -> f64 {
        self.model.enthalpy(k, self.p, t)
    }

    /// Volumetric enthalpy `h_k(T)/v_k(T)` of a pure phase at temperature `t`
    /// (used to build initial conditions; the model carries `U` and `ρ`
    /// separately).
    pub fn u_pure(&self, k: usize, t: f64) -> f64 {
        self.h_pure(k, t) / v_pt(&self.eos()[k], self.p, t)
    }

    /// Specific volume `Σ y_k v_k(p,T)` of a mixture.
    fn mix_v(&self, y: &Fr, t: f64) -> f64 {
        (0..self.np()).map(|k| y[k] * v_pt(&self.eos()[k], self.p, t)).sum()
    }

    /// Solve `h_k(T) = h` for a pure phase on the bracket `[lo, hi]`, where
    /// `h_k` is monotone. For a NASG phase without a magnetic or thermal term
    /// both `h` and `v` are linear in `T`, so `T = (h − c₀)/c_p` with
    /// `h = c_p T + c₀`; the magnetic term (iron bcc) and the tabulated
    /// thermal EOS make `h_k(T)` nonlinear, so those cases use a safeguarded
    /// Newton iteration.
    pub fn solve_h(&self, k: usize, h: f64, lo: f64, hi: f64) -> f64 {
        if self.model.magnetic[k].is_none() && self.model.thermal[k].is_none() {
            let e = &self.eos()[k];
            let cp = e.gamma * e.cv;
            let c0 = e.q + self.p * e.b;
            let t = (h - c0) / cp;
            if t.is_finite() && t > 0.0 {
                return t.clamp(lo.min(hi), lo.max(hi));
            }
        }
        let (mut a, mut b) = (lo, hi);
        let fa = self.h_pure(k, a) - h;
        let fb = self.h_pure(k, b) - h;
        if fa == 0.0 {
            return a;
        }
        if fb == 0.0 {
            return b;
        }
        // If the root is outside the bracket, clamp to the near end.
        if fa * fb > 0.0 {
            return if fa.abs() < fb.abs() { a } else { b };
        }
        let mut t = 0.5 * (a + b);
        for _ in 0..80 {
            let f = self.h_pure(k, t) - h;
            if f.abs() < 1e-13 * h.abs().max(1.0) {
                return t;
            }
            if f * fa < 0.0 {
                b = t;
            } else {
                a = t;
            }
            let cp = self.model.cp(k, self.p, t);
            let tn = if cp.abs() > 0.0 { t - f / cp } else { f64::NAN };
            t = if tn > a && tn < b { tn } else { 0.5 * (a + b) };
        }
        t
    }

    /// Equilibrium `(T, y)` of a **specific** enthalpy `h` at the ladder
    /// pressure. Exact piecewise inversion.
    pub fn state_h(&self, h: f64) -> (f64, Fr) {
        let mut y = [0.0; MAXP];
        let mut phase = self.ladder.lowest;
        let mut prev_t = self.t_lo;
        // Below the built range: stay in the lowest phase.
        if h <= self.h_pure(phase, prev_t) {
            y[phase] = 1.0;
            return (prev_t, y);
        }
        for &(tt, a, b) in &self.ladder.transitions {
            // pure phase a on [prev_t, tt]
            let ha = self.h_pure(a, tt);
            if h <= ha {
                let t = self.solve_h(a, h, prev_t, tt);
                y[a] = 1.0;
                return (t, y);
            }
            // two-phase window at tt
            let hb = self.h_pure(b, tt);
            if h <= hb {
                let phi = ((h - ha) / (hb - ha)).clamp(0.0, 1.0);
                y[b] = phi;
                y[a] = 1.0 - phi;
                return (tt, y);
            }
            prev_t = tt;
            phase = b;
        }
        // Above the last transition: pure top phase.
        let hi = self.t_hi.max(2.0 * prev_t);
        let t = if h <= self.h_pure(phase, prev_t) { prev_t } else { self.solve_h(phase, h, prev_t, hi) };
        y[phase] = 1.0;
        (t, y)
    }

    // --- the old volumetric inversion, retained for the zero-flow closure ---

    /// Solve `U_k(T) = u` for a pure phase on the bracket `[lo, hi]`, where
    /// `U_k = h_k/v_k` is monotone.
    pub fn solve_pure(&self, k: usize, u: f64, lo: f64, hi: f64) -> f64 {
        if self.model.magnetic[k].is_none() && self.model.thermal[k].is_none() {
            let e = &self.eos()[k];
            let cp = e.gamma * e.cv;
            let c0 = e.q + self.p * e.b;
            let a1 = (e.gamma - 1.0) * e.cv / (self.p + e.p_inf);
            let t = (c0 - u * e.b) / (u * a1 - cp);
            if t.is_finite() && t > 0.0 {
                return t.clamp(lo.min(hi), lo.max(hi));
            }
        }
        let (mut a, mut b) = (lo, hi);
        let fa = self.u_pure(k, a) - u;
        let fb = self.u_pure(k, b) - u;
        if fa == 0.0 {
            return a;
        }
        if fb == 0.0 {
            return b;
        }
        if fa * fb > 0.0 {
            return if fa.abs() < fb.abs() { a } else { b };
        }
        let mut t = 0.5 * (a + b);
        for _ in 0..50 {
            let f = self.u_pure(k, t) - u;
            if f.abs() < 1e-13 * u.abs().max(1.0) {
                return t;
            }
            if f * fa < 0.0 {
                b = t;
            } else {
                a = t;
            }
            let (hh, cp) = (self.model.enthalpy(k, self.p, t), self.model.cp(k, self.p, t));
            let v = v_pt(&self.eos()[k], self.p, t);
            let vp = (v - self.eos()[k].b) / t;
            let d = (cp * v - hh * vp) / (v * v);
            let tn = if d.abs() > 0.0 { t - f / d } else { f64::NAN };
            t = if tn > a && tn < b { tn } else { 0.5 * (a + b) };
        }
        t
    }

    #[inline]
    fn within(u: f64, a: f64, b: f64) -> bool {
        let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
        u >= lo && u <= hi
    }

    /// Equilibrium `(T, y)` of a volumetric enthalpy `u` at the ladder
    /// pressure (the zero-flow closure's inversion).
    pub fn state(&self, u: f64) -> (f64, Fr) {
        let p = self.p;
        let mut y = [0.0; MAXP];
        let mut phase = self.ladder.lowest;
        let mut prev_t = self.t_lo;
        if u <= self.u_pure(phase, prev_t) {
            y[phase] = 1.0;
            return (prev_t, y);
        }
        for &(tt, a, b) in &self.ladder.transitions {
            let ua = self.u_pure(a, prev_t);
            let ub = self.u_pure(a, tt);
            if Self::within(u, ua, ub) {
                let t = self.solve_pure(a, u, prev_t, tt);
                y[a] = 1.0;
                return (t, y);
            }
            let (ha, hb) = (self.model.enthalpy(a, p, tt), self.model.enthalpy(b, p, tt));
            let (va, vb) = (v_pt(&self.eos()[a], p, tt), v_pt(&self.eos()[b], p, tt));
            let (ua2, ub2) = (ha / va, hb / vb);
            if Self::within(u, ua2, ub2) {
                let denom = (hb - ha) - u * (vb - va);
                let phi = if denom.abs() > 0.0 { ((u * va - ha) / denom).clamp(0.0, 1.0) } else { 0.0 };
                y[b] = phi;
                y[a] = 1.0 - phi;
                return (tt, y);
            }
            prev_t = tt;
            phase = b;
        }
        let hi = self.t_hi.max(2.0 * prev_t);
        let t = if u <= self.u_pure(phase, prev_t) { prev_t } else { self.solve_pure(phase, u, prev_t, hi) };
        y[phase] = 1.0;
        (t, y)
    }

    /// Equilibrium state of one cell: `(T, y, ρ_eq)`.
    pub fn state_at(&self, i: usize) -> (f64, Fr, f64) {
        let (t, y) = if self.expansion {
            self.state_h(self.u[i] / self.rho[i])
        } else {
            self.state(self.u[i])
        };
        let v = self.mix_v(&y, t);
        (t, y, 1.0 / v)
    }

    /// Temperature of every cell.
    pub fn temperatures(&self) -> Vec<f64> {
        (0..self.len()).map(|i| self.state_at(i).0).collect()
    }

    /// Equilibrium state of every cell: `(T, y, ρ_eq)`.
    pub fn states(&self) -> Vec<(f64, Fr, f64)> {
        (0..self.len()).map(|i| self.state_at(i)).collect()
    }

    /// Volume-fraction-weighted conductivity of a cell.
    pub fn cell_conductivity(&self, i: usize, t: f64) -> f64 {
        let y = &self.y[i];
        let mut vtot = 0.0;
        let mut kv = 0.0;
        for k in 0..self.np() {
            if y[k] > 0.0 {
                let vk = y[k] * v_pt(&self.eos()[k], self.p, t);
                vtot += vk;
                kv += vk * self.k[k];
            }
        }
        if vtot > 0.0 {
            kv / vtot
        } else {
            0.0
        }
    }

    /// Sensible volumetric heat capacity `ρ min_k c_{p,k}` over the phases
    /// present, used only for the explicit stability limit.
    fn sensible_cvol(&self, i: usize, t: f64) -> f64 {
        let y = &self.y[i];
        let cp_min = (0..self.np())
            .filter(|&k| y[k] > 1e-12)
            .map(|k| self.model.cp(k, self.p, t))
            .fold(f64::INFINITY, f64::min);
        if cp_min.is_finite() && cp_min > 0.0 && self.rho[i] > 0.0 {
            self.rho[i] * cp_min
        } else {
            0.0
        }
    }

    /// Attach a per-cell ghost temperature to one side,
    /// `side ∈ {0=west, 1=east, 2=south, 3=north}`; `vals` has length `ny`
    /// for the west/east sides and `nx` for the south/north sides.
    pub fn set_dirichlet_row(&mut self, side: usize, vals: Vec<f64>) {
        assert!(side < 4);
        let want = if side < 2 { self.ny } else { self.nx };
        assert_eq!(vals.len(), want, "boundary row length");
        self.bc_row[side] = Some(vals);
    }

    /// Set one cell to a pure phase at temperature `t`.
    pub fn set_cell(&mut self, i: usize, j: usize, ph: usize, t: f64) {
        let idx = self.idx(i, j);
        self.u[idx] = self.u_pure(ph, t);
        self.rho[idx] = 1.0 / v_pt(&self.eos()[ph], self.p, t);
        let mut y = [0.0; MAXP];
        y[ph] = 1.0;
        self.y[idx] = y;
    }

    /// Set every cell from a temperature field `f(x, y)`, taking the stable
    /// phase from the ladder.
    pub fn set_temperatures(&mut self, f: impl Fn(f64, f64) -> f64) {
        for j in 0..self.ny {
            for i in 0..self.nx {
                let x = (i as f64 + 0.5) * self.dx;
                let yc = (j as f64 + 0.5) * self.dy;
                let t = f(x, yc);
                let ph = self.stable_phase(t);
                self.set_cell(i, j, ph, t);
            }
        }
    }

    /// Stable phase at temperature `t` on the ladder.
    pub fn stable_phase(&self, t: f64) -> usize {
        let mut cur = self.ladder.lowest;
        for &(tt, _, b) in &self.ladder.transitions {
            if t > tt {
                cur = b;
            }
        }
        cur
    }

    /// Total energy `∫ U dV` (per unit depth).
    pub fn energy(&self) -> f64 {
        self.u.iter().sum::<f64>() * self.dx * self.dy
    }

    /// Total mass `∫ ρ dV` (per unit depth).
    pub fn mass(&self) -> f64 {
        self.rho.iter().sum::<f64>() * self.dx * self.dy
    }

    /// Area occupied by phase `ph`, from the equilibrium volume fractions.
    pub fn phase_area(&self, ph: usize) -> f64 {
        let cell = self.dx * self.dy;
        let mut a = 0.0;
        for i in 0..self.len() {
            let (t, y) = if self.expansion {
                self.state_h(self.u[i] / self.rho[i])
            } else {
                self.state(self.u[i])
            };
            // volume of phase ph = (mass of ph) × v_ph
            a += self.rho[i] * cell * y[ph] * v_pt(&self.eos()[ph], self.p, t);
        }
        a
    }

    /// Boussinesq buoyancy acceleration `g β (T − T_ref)` per cell, the
    /// momentum source a flow solver would take when the conduction-only
    /// assumption fails.
    pub fn buoyancy(&self, g: f64, beta: f64, t_ref: f64) -> Vec<f64> {
        self.temperatures().iter().map(|&t| g * beta * (t - t_ref)).collect()
    }

    // ------------------------------------------------------ expansion Poisson

    /// Is there any open (Dirichlet) side? If not, the Poisson problem is
    /// singular (constant null space) and the source must be mean-zero.
    fn has_open(&self) -> bool {
        (!self.bc.periodic_x && (self.flow_bc.west == FlowBc::Open || self.flow_bc.east == FlowBc::Open))
            || (!self.bc.periodic_y && (self.flow_bc.south == FlowBc::Open || self.flow_bc.north == FlowBc::Open))
    }

    /// Apply the variable-coefficient operator `∇·(ρ ∇ψ)` into `out`, zeroing
    /// `out` first, with the same face coefficients used by [`Self::advect`].
    fn apply_l_into(&self, psi: &[f64], rho: &[f64], out: &mut [f64]) {
        for v in out.iter_mut() {
            *v = 0.0;
        }
        let (nx, ny, dx, dy) = (self.nx, self.ny, self.dx, self.dy);
        // x interior
        if nx > 1 {
            let inv2 = 1.0 / (dx * dx);
            for j in 0..ny {
                for i in 0..nx - 1 {
                    let a = j * nx + i;
                    let b = a + 1;
                    let c = harm(rho[a], rho[b]) * inv2;
                    out[a] += c * (psi[b] - psi[a]);
                    out[b] += c * (psi[a] - psi[b]);
                }
            }
        }
        // x periodic (last↔first)
        if self.bc.periodic_x && nx > 1 {
            let inv2 = 1.0 / (dx * dx);
            for j in 0..ny {
                let a = j * nx;
                let b = j * nx + nx - 1;
                let c = harm(rho[a], rho[b]) * inv2;
                out[a] += c * (psi[b] - psi[a]);
                out[b] += c * (psi[a] - psi[b]);
            }
        }
        if !self.bc.periodic_x {
            let inv2 = 1.0 / (dx * dx);
            if self.flow_bc.west == FlowBc::Open {
                for j in 0..ny {
                    let a = j * nx;
                    let c = 2.0 * rho[a] * inv2;
                    out[a] += c * (0.0 - psi[a]);
                }
            }
            if self.flow_bc.east == FlowBc::Open {
                for j in 0..ny {
                    let a = j * nx + nx - 1;
                    let c = 2.0 * rho[a] * inv2;
                    out[a] += c * (0.0 - psi[a]);
                }
            }
        }
        // y interior
        if ny > 1 {
            let inv2 = 1.0 / (dy * dy);
            for j in 0..ny - 1 {
                for i in 0..nx {
                    let a = j * nx + i;
                    let b = a + nx;
                    let c = harm(rho[a], rho[b]) * inv2;
                    out[a] += c * (psi[b] - psi[a]);
                    out[b] += c * (psi[a] - psi[b]);
                }
            }
        }
        // y periodic
        if self.bc.periodic_y && ny > 1 {
            let inv2 = 1.0 / (dy * dy);
            for i in 0..nx {
                let a = i;
                let b = (ny - 1) * nx + i;
                let c = harm(rho[a], rho[b]) * inv2;
                out[a] += c * (psi[b] - psi[a]);
                out[b] += c * (psi[a] - psi[b]);
            }
        }
        if !self.bc.periodic_y {
            let inv2 = 1.0 / (dy * dy);
            if self.flow_bc.south == FlowBc::Open {
                for i in 0..nx {
                    let a = i;
                    let c = 2.0 * rho[a] * inv2;
                    out[a] += c * (0.0 - psi[a]);
                }
            }
            if self.flow_bc.north == FlowBc::Open {
                for i in 0..nx {
                    let a = (ny - 1) * nx + i;
                    let c = 2.0 * rho[a] * inv2;
                    out[a] += c * (0.0 - psi[a]);
                }
            }
        }
    }

    /// Diagonal of [`Self::apply_l_into`]: the sum of the face coefficients
    /// incident on each cell (`L(1)` is zero in the interior, so this cannot be
    /// obtained by applying the operator to a constant).
    fn diag_l(&self, rho: &[f64]) -> Vec<f64> {
        let (nx, ny, dx, dy) = (self.nx, self.ny, self.dx, self.dy);
        let n = rho.len();
        let mut d = vec![0.0; n];
        if nx > 1 {
            let inv2 = 1.0 / (dx * dx);
            for j in 0..ny {
                for i in 0..nx - 1 {
                    let a = j * nx + i;
                    let b = a + 1;
                    let c = harm(rho[a], rho[b]) * inv2;
                    d[a] += c;
                    d[b] += c;
                }
            }
        }
        if self.bc.periodic_x && nx > 1 {
            let inv2 = 1.0 / (dx * dx);
            for j in 0..ny {
                let a = j * nx;
                let b = j * nx + nx - 1;
                let c = harm(rho[a], rho[b]) * inv2;
                d[a] += c;
                d[b] += c;
            }
        }
        if !self.bc.periodic_x {
            let inv2 = 1.0 / (dx * dx);
            if self.flow_bc.west == FlowBc::Open {
                for j in 0..ny {
                    d[j * nx] += 2.0 * rho[j * nx] * inv2;
                }
            }
            if self.flow_bc.east == FlowBc::Open {
                for j in 0..ny {
                    let a = j * nx + nx - 1;
                    d[a] += 2.0 * rho[a] * inv2;
                }
            }
        }
        if ny > 1 {
            let inv2 = 1.0 / (dy * dy);
            for j in 0..ny - 1 {
                for i in 0..nx {
                    let a = j * nx + i;
                    let b = a + nx;
                    let c = harm(rho[a], rho[b]) * inv2;
                    d[a] += c;
                    d[b] += c;
                }
            }
        }
        if self.bc.periodic_y && ny > 1 {
            let inv2 = 1.0 / (dy * dy);
            for i in 0..nx {
                let a = i;
                let b = (ny - 1) * nx + i;
                let c = harm(rho[a], rho[b]) * inv2;
                d[a] += c;
                d[b] += c;
            }
        }
        if !self.bc.periodic_y {
            let inv2 = 1.0 / (dy * dy);
            if self.flow_bc.south == FlowBc::Open {
                for i in 0..nx {
                    d[i] += 2.0 * rho[i] * inv2;
                }
            }
            if self.flow_bc.north == FlowBc::Open {
                for i in 0..nx {
                    let a = (ny - 1) * nx + i;
                    d[a] += 2.0 * rho[a] * inv2;
                }
            }
        }
        d
    }

    /// Solve a 1D variable-coefficient Poisson problem directly by the Thomas
    /// algorithm. Only called for a single row/column with an open end, so the
    /// tridiagonal matrix is non-singular.
    fn solve_tridiag(&self, rhs: &[f64], rho: &[f64], d: f64, lo: FlowBc, hi: FlowBc) -> Vec<f64> {
        let m = rho.len();
        let inv2 = 1.0 / (d * d);
        let mut sub = vec![0.0; m];
        let mut diag = vec![0.0; m];
        let mut sup = vec![0.0; m];
        let mut s = vec![0.0; m];
        for k in 0..m - 1 {
            let c = harm(rho[k], rho[k + 1]) * inv2;
            diag[k] += c;
            diag[k + 1] += c;
            sup[k] = -c;
            sub[k + 1] = -c;
        }
        if lo == FlowBc::Open {
            diag[0] += 2.0 * rho[0] * inv2;
        }
        if hi == FlowBc::Open {
            diag[m - 1] += 2.0 * rho[m - 1] * inv2;
        }
        for k in 0..m {
            // solve A ψ = −rhs with A = −L positive definite
            s[k] = -rhs[k];
        }
        for k in 1..m {
            let w = sub[k] / diag[k - 1];
            diag[k] -= w * sup[k - 1];
            s[k] -= w * s[k - 1];
        }
        let mut x = vec![0.0; m];
        x[m - 1] = s[m - 1] / diag[m - 1];
        for k in (0..m - 1).rev() {
            x[k] = (s[k] - sup[k] * x[k + 1]) / diag[k];
        }
        x
    }

    /// Solve `∇·(ρ ∇ψ) = rhs`, warm-started from `self.psi`, and cache the
    /// result. 1D grids use a direct tridiagonal solve; otherwise conjugate
    /// gradients with a Jacobi preconditioner. With no open side the constant
    /// null space is removed (the source is made mean-zero).
    fn solve_poisson(&mut self, rhs_in: &[f64]) -> Vec<f64> {
        let n = rhs_in.len();
        let rho = self.rho.clone();
        let singular = !self.has_open();
        let mut rhs = rhs_in.to_vec();
        if singular {
            let mean = rhs.iter().sum::<f64>() / n as f64;
            for r in rhs.iter_mut() {
                *r -= mean;
            }
        }
        // direct solve for a single row/column with an open end
        if !singular && self.ny == 1 && !self.bc.periodic_x {
            let x = self.solve_tridiag(&rhs, &rho, self.dx, self.flow_bc.west, self.flow_bc.east);
            self.psi = x.clone();
            return x;
        }
        if !singular && self.nx == 1 && !self.bc.periodic_y {
            let x = self.solve_tridiag(&rhs, &rho, self.dy, self.flow_bc.south, self.flow_bc.north);
            self.psi = x.clone();
            return x;
        }
        let diag = self.diag_l(&rho);
        let mut x = if self.psi.len() == n { self.psi.clone() } else { vec![0.0; n] };
        if singular {
            let mean = x.iter().sum::<f64>() / n as f64;
            for v in x.iter_mut() {
                *v -= mean;
            }
        }
        let mut r = vec![0.0; n];
        self.apply_l_into(&x, &rho, &mut r);
        for i in 0..n {
            r[i] = rhs[i] - r[i];
        }
        if singular {
            let mean = r.iter().sum::<f64>() / n as f64;
            for v in r.iter_mut() {
                *v -= mean;
            }
        }
        let mut z = vec![0.0; n];
        let mut p = vec![0.0; n];
        let mut ap = vec![0.0; n];
        for i in 0..n {
            z[i] = r[i] / diag[i].max(1e-300);
            p[i] = z[i];
        }
        let mut rz: f64 = r.iter().zip(&z).map(|(a, b)| a * b).sum();
        let rhs_norm = rhs.iter().map(|a| a * a).sum::<f64>().sqrt().max(1e-300);
        let max_iter = (4 * n).min(4000).max(40);
        for _ in 0..max_iter {
            self.apply_l_into(&p, &rho, &mut ap);
            let pap: f64 = p.iter().zip(&ap).map(|(a, b)| a * b).sum();
            if pap.abs() < 1e-300 {
                break;
            }
            let alpha = rz / pap;
            for i in 0..n {
                x[i] += alpha * p[i];
                r[i] -= alpha * ap[i];
            }
            if singular {
                let mean = r.iter().sum::<f64>() / n as f64;
                for v in r.iter_mut() {
                    *v -= mean;
                }
            }
            let rnorm = r.iter().map(|a| a * a).sum::<f64>().sqrt();
            if rnorm < 1e-10 * rhs_norm {
                break;
            }
            for i in 0..n {
                z[i] = r[i] / diag[i].max(1e-300);
            }
            let rz_new: f64 = r.iter().zip(&z).map(|(a, b)| a * b).sum();
            let beta = rz_new / rz;
            for i in 0..n {
                p[i] = z[i] + beta * p[i];
            }
            rz = rz_new;
        }
        self.psi = x.clone();
        x
    }

    /// Largest `|u|` of the expansion flow from a potential `ψ`.
    fn max_flow_speed(&self, psi: &[f64]) -> f64 {
        let (nx, ny, dx, dy) = (self.nx, self.ny, self.dx, self.dy);
        let mut m: f64 = 0.0;
        if nx > 1 {
            for j in 0..ny {
                for i in 0..nx - 1 {
                    let a = j * nx + i;
                    m = m.max((psi[a + 1] - psi[a]).abs() / dx);
                }
            }
        }
        if self.bc.periodic_x && nx > 1 {
            for j in 0..ny {
                let a = j * nx;
                let b = j * nx + nx - 1;
                m = m.max((psi[b] - psi[a]).abs() / dx);
            }
        }
        if !self.bc.periodic_x {
            if self.flow_bc.west == FlowBc::Open {
                for j in 0..ny {
                    let a = j * nx;
                    m = m.max(2.0 * psi[a].abs() / dx);
                }
            }
            if self.flow_bc.east == FlowBc::Open {
                for j in 0..ny {
                    let a = j * nx + nx - 1;
                    m = m.max(2.0 * psi[a].abs() / dx);
                }
            }
        }
        if ny > 1 {
            for j in 0..ny - 1 {
                for i in 0..nx {
                    let a = j * nx + i;
                    m = m.max((psi[a + nx] - psi[a]).abs() / dy);
                }
            }
        }
        if self.bc.periodic_y && ny > 1 {
            for i in 0..nx {
                let a = i;
                let b = (ny - 1) * nx + i;
                m = m.max((psi[b] - psi[a]).abs() / dy);
            }
        }
        if !self.bc.periodic_y {
            if self.flow_bc.south == FlowBc::Open {
                for i in 0..nx {
                    m = m.max(2.0 * psi[i].abs() / dy);
                }
            }
            if self.flow_bc.north == FlowBc::Open {
                for i in 0..nx {
                    let a = (ny - 1) * nx + i;
                    m = m.max(2.0 * psi[a].abs() / dy);
                }
            }
        }
        m
    }

    /// Conservative advection of mass and volumetric enthalpy with the mass
    /// flux `ρ ∇ψ`. Uses the post-conduction enthalpies `hpred` for upwinding.
    /// The flux coefficients are exactly those of [`Self::apply_l`], so the
    /// mass update lands on the EOS density of the conduction step.
    fn advect(&mut self, psi: &[f64], hpred: &[f64], hstep: f64) {
        let (nx, ny, dx, dy) = (self.nx, self.ny, self.dx, self.dy);
        let mut rho = self.rho.clone();
        let mut u = self.u.clone();
        let add_face = |a: usize, b: usize, f: f64, rho: &mut Vec<f64>, u: &mut Vec<f64>| {
            let hup = if f >= 0.0 { hpred[a] } else { hpred[b] };
            rho[a] -= hstep * f / dx;
            rho[b] += hstep * f / dx;
            u[a] -= hstep * hup * f / dx;
            u[b] += hstep * hup * f / dx;
        };
        // x interior
        if nx > 1 {
            for j in 0..ny {
                for i in 0..nx - 1 {
                    let a = j * nx + i;
                    let b = a + 1;
                    let f = harm(self.rho[a], self.rho[b]) * (psi[b] - psi[a]) / dx;
                    add_face(a, b, f, &mut rho, &mut u);
                }
            }
        }
        if self.bc.periodic_x && nx > 1 {
            for j in 0..ny {
                let a = j * nx;
                let b = j * nx + nx - 1;
                let f = harm(self.rho[a], self.rho[b]) * (psi[b] - psi[a]) / dx;
                add_face(a, b, f, &mut rho, &mut u);
            }
        }
        if !self.bc.periodic_x {
            if self.flow_bc.west == FlowBc::Open {
                for j in 0..ny {
                    let a = j * nx;
                    // flux from the ghost into a (positive +x)
                    let f = 2.0 * self.rho[a] * psi[a] / dx;
                    rho[a] += hstep * f / dx;
                    u[a] += hstep * hpred[a] * f / dx;
                }
            }
            if self.flow_bc.east == FlowBc::Open {
                for j in 0..ny {
                    let a = j * nx + nx - 1;
                    // flux from a into the ghost (positive +x)
                    let f = 2.0 * self.rho[a] * psi[a] / dx;
                    rho[a] -= hstep * f / dx;
                    u[a] -= hstep * hpred[a] * f / dx;
                }
            }
        }
        // y interior
        if ny > 1 {
            for j in 0..ny - 1 {
                for i in 0..nx {
                    let a = j * nx + i;
                    let b = a + nx;
                    let f = harm(self.rho[a], self.rho[b]) * (psi[b] - psi[a]) / dy;
                    let hup = if f >= 0.0 { hpred[a] } else { hpred[b] };
                    rho[a] -= hstep * f / dy;
                    rho[b] += hstep * f / dy;
                    u[a] -= hstep * hup * f / dy;
                    u[b] += hstep * hup * f / dy;
                }
            }
        }
        if self.bc.periodic_y && ny > 1 {
            for i in 0..nx {
                let a = i;
                let b = (ny - 1) * nx + i;
                let f = harm(self.rho[a], self.rho[b]) * (psi[b] - psi[a]) / dy;
                let hup = if f >= 0.0 { hpred[a] } else { hpred[b] };
                rho[a] -= hstep * f / dy;
                rho[b] += hstep * f / dy;
                u[a] -= hstep * hup * f / dy;
                u[b] += hstep * hup * f / dy;
            }
        }
        if !self.bc.periodic_y {
            if self.flow_bc.south == FlowBc::Open {
                for i in 0..nx {
                    let a = i;
                    let f = 2.0 * self.rho[a] * psi[a] / dy;
                    rho[a] += hstep * f / dy;
                    u[a] += hstep * hpred[a] * f / dy;
                }
            }
            if self.flow_bc.north == FlowBc::Open {
                for i in 0..nx {
                    let a = (ny - 1) * nx + i;
                    let f = 2.0 * self.rho[a] * psi[a] / dy;
                    rho[a] -= hstep * f / dy;
                    u[a] -= hstep * hpred[a] * f / dy;
                }
            }
        }
        self.rho = rho;
        self.u = u;
    }

    // ------------------------------------------------------------- time step

    /// Advance by `dt` in stable explicit sub-steps; returns their number.
    pub fn advance(&mut self, dt: f64) -> usize {
        let (nx, ny, dx, dy) = (self.nx, self.ny, self.dx, self.dy);
        let n = self.len();
        let mut done = 0.0;
        let mut steps = 0;
        while done < dt * (1.0 - 1e-12) {
            // 1. equilibrium state and phase-dependent conductivity
            let mut temps = vec![0.0; n];
            for i in 0..n {
                let (t, y) = if self.expansion {
                    self.state_h(self.u[i] / self.rho[i])
                } else {
                    let (t, y) = self.state(self.u[i]);
                    // zero-flow closure: density follows the EOS directly
                    self.rho[i] = 1.0 / self.mix_v(&y, t);
                    (t, y)
                };
                temps[i] = t;
                self.y[i] = y;
            }
            let kc: Vec<f64> = (0..n).map(|i| self.cell_conductivity(i, temps[i])).collect();

            let mut du = vec![0.0; n];
            let mut g = vec![0.0; n];

            // x faces
            let xpairs: Vec<(usize, usize, usize, usize)> = if nx > 1 {
                let mut v = Vec::with_capacity(nx * ny);
                for j in 0..ny {
                    for i in 0..nx.saturating_sub(1) {
                        v.push((i, j, i + 1, j));
                    }
                }
                if self.bc.periodic_x {
                    for j in 0..ny {
                        v.push((nx - 1, j, 0, j));
                    }
                }
                v
            } else {
                Vec::new()
            };
            for &(i, j, ii, jj) in &xpairs {
                let a = self.idx(i, j);
                let b = self.idx(ii, jj);
                let kf = harm(kc[a], kc[b]);
                let q = kf * (temps[a] - temps[b]) / dx;
                du[a] -= q / dx;
                du[b] += q / dx;
                g[a] += kf / (dx * dx);
                g[b] += kf / (dx * dx);
            }
            // y faces
            let ypairs: Vec<(usize, usize, usize, usize)> = if ny > 1 {
                let mut v = Vec::with_capacity(nx * ny);
                for j in 0..ny.saturating_sub(1) {
                    for i in 0..nx {
                        v.push((i, j, i, j + 1));
                    }
                }
                if self.bc.periodic_y {
                    for i in 0..nx {
                        v.push((i, ny - 1, i, 0));
                    }
                }
                v
            } else {
                Vec::new()
            };
            for &(i, j, ii, jj) in &ypairs {
                let a = self.idx(i, j);
                let b = self.idx(ii, jj);
                let kf = harm(kc[a], kc[b]);
                let q = kf * (temps[a] - temps[b]) / dy;
                du[a] -= q / dy;
                du[b] += q / dy;
                g[a] += kf / (dy * dy);
                g[b] += kf / (dy * dy);
            }
            // Dirichlet faces (half a cell out); skipped when periodic.
            if nx > 1 && !self.bc.periodic_x {
                if let Bc2D::Dirichlet(tw) = self.bc.west {
                    for j in 0..ny {
                        let a = self.idx(0, j);
                        let tg = self.bc_row[0].as_ref().map_or(tw, |v| v[j]);
                        let q = 2.0 * kc[a] * (temps[a] - tg) / dx;
                        du[a] -= q / dx;
                        g[a] += 2.0 * kc[a] / (dx * dx);
                    }
                }
                if let Bc2D::Dirichlet(tw) = self.bc.east {
                    for j in 0..ny {
                        let a = self.idx(nx - 1, j);
                        let tg = self.bc_row[1].as_ref().map_or(tw, |v| v[j]);
                        let q = 2.0 * kc[a] * (temps[a] - tg) / dx;
                        du[a] -= q / dx;
                        g[a] += 2.0 * kc[a] / (dx * dx);
                    }
                }
            }
            if nx == 1 {
                for bc in [self.bc.west, self.bc.east] {
                    if let Bc2D::Dirichlet(tw) = bc {
                        let a = self.idx(0, 0);
                        let q = 2.0 * kc[a] * (temps[a] - tw) / dx;
                        du[a] -= q / dx;
                        g[a] += 2.0 * kc[a] / (dx * dx);
                    }
                }
            }
            if ny > 1 && !self.bc.periodic_y {
                if let Bc2D::Dirichlet(tw) = self.bc.south {
                    for i in 0..nx {
                        let a = self.idx(i, 0);
                        let tg = self.bc_row[2].as_ref().map_or(tw, |v| v[i]);
                        let q = 2.0 * kc[a] * (temps[a] - tg) / dy;
                        du[a] -= q / dy;
                        g[a] += 2.0 * kc[a] / (dy * dy);
                    }
                }
                if let Bc2D::Dirichlet(tw) = self.bc.north {
                    for i in 0..nx {
                        let a = self.idx(i, ny - 1);
                        let tg = self.bc_row[3].as_ref().map_or(tw, |v| v[i]);
                        let q = 2.0 * kc[a] * (temps[a] - tg) / dy;
                        du[a] -= q / dy;
                        g[a] += 2.0 * kc[a] / (dy * dy);
                    }
                }
            }

            // stability limit (conduction)
            let mut lim = f64::INFINITY;
            for idx in 0..n {
                let cvol = self.sensible_cvol(idx, temps[idx]);
                if cvol > 0.0 && g[idx] > 0.0 {
                    lim = lim.min(0.4 * cvol / g[idx]);
                }
            }
            let mut hstep = lim.min(dt - done).max(f64::MIN_POSITIVE);

            // 2. expansion flow
            let mut hpred = vec![0.0; n];
            let mut psi = vec![0.0; n];
            let mut advect_now = false;
            if self.expansion && n > 0 {
                let dmin = dx.min(dy);
                let rho_scale = self.rho.iter().cloned().fold(f64::MIN_POSITIVE, f64::max);
                let mut rho_eq = vec![0.0; n];
                let mut rhs = vec![0.0; n];
                let mut dmax = 0.0f64;
                for i in 0..n {
                    hpred[i] = (self.u[i] + hstep * (du[i] + self.source[i])) / self.rho[i];
                    let (tp, yp) = self.state_h(hpred[i]);
                    rho_eq[i] = 1.0 / self.mix_v(&yp, tp);
                    dmax = dmax.max((self.rho[i] - rho_eq[i]).abs());
                    rhs[i] = (self.rho[i] - rho_eq[i]) / hstep;
                }
                // A relative density change below 1e-6 per step produces an
                // expansion velocity far below the conduction velocity scale;
                // solve the flow only when it is larger.
                if dmax < 1e-6 * rho_scale {
                    // the density change over this step is negligible for the
                    // flow: apply the zero-flow projection (keep the EOS
                    // density, do not advect).
                    for i in 0..n {
                        self.rho[i] = rho_eq[i];
                    }
                } else {
                    for _attempt in 0..3 {
                        psi = self.solve_poisson(&rhs);
                        let umax = self.max_flow_speed(&psi);
                        if umax * hstep > 0.4 * dmin {
                            // CFL: shrink the step and rebuild the source
                            hstep = 0.4 * dmin / umax.max(f64::MIN_POSITIVE);
                            for i in 0..n {
                                hpred[i] = (self.u[i] + hstep * (du[i] + self.source[i])) / self.rho[i];
                                let (tp, yp) = self.state_h(hpred[i]);
                                rho_eq[i] = 1.0 / self.mix_v(&yp, tp);
                                rhs[i] = (self.rho[i] - rho_eq[i]) / hstep;
                            }
                        } else {
                            advect_now = true;
                            break;
                        }
                    }
                }
            }

            // 3. apply conduction + source, then advect
            for i in 0..n {
                self.u[i] += hstep * (du[i] + self.source[i]);
            }
            if advect_now {
                self.advect(&psi, &hpred, hstep);
            }
            done += hstep;
            self.t += hstep;
            steps += 1;
        }
        steps
    }

    /// Liquid mass fraction of a cell (the fraction of the phase marked as
    /// `liquid` in the phase set). Used by the Boussinesq / Darcy coupling in
    /// [`crate::convection2d`]: `0` = rigid solid, `1` = free liquid.
    pub fn liquid_fraction(&self, i: usize) -> f64 {
        let (_, y) = if self.expansion {
            self.state_h(self.u[i] / self.rho[i])
        } else {
            self.state(self.u[i])
        };
        y[self.model.liquid]
    }

    /// Conservative advection of the volumetric enthalpy `U` by a prescribed
    /// staggered velocity field (the unfold from conduction to natural
    /// convection). `ux` holds the `x`-face velocities, row-major
    /// `(nx+1) × ny` at `x = i dx`, so `ux[j*(nx+1)+i]` is the flux between
    /// cells `i−1` and `i`; `uy` holds the `y`-face velocities, row-major
    /// `nx × (ny+1)` at `y = j dy`. Wall faces carry zero velocity, so nothing
    /// leaves the box and every interior flux telescopes: `∫U dV` is conserved
    /// to round-off. The face enthalpy is a minmod-limited second-order upwind
    /// value, which keeps the moving melt front monotone. The velocity is
    /// sub-stepped by CFL 0.4.
    pub fn advect_enthalpy(&mut self, ux: &[f64], uy: &[f64], dt: f64) {
        let (nx, ny, dx, dy) = (self.nx, self.ny, self.dx, self.dy);
        assert_eq!(ux.len(), (nx + 1) * ny, "ux shape");
        assert_eq!(uy.len(), nx * (ny + 1), "uy shape");
        let n = self.len();
        let mut umax: f64 = 0.0;
        for &a in ux {
            umax = umax.max(a.abs());
        }
        for &a in uy {
            umax = umax.max(a.abs());
        }
        let dmin = dx.min(dy);
        let nsub = if umax > 0.0 { (umax * dt / (0.4 * dmin)).ceil() as usize } else { 1 };
        let nsub = nsub.max(1);
        let h = dt / nsub as f64;
        let mut hc = vec![0.0; n];
        let mut dudt = vec![0.0; n];
        for _ in 0..nsub {
            for c in 0..n {
                hc[c] = self.u[c] / self.rho[c];
            }
            for v in dudt.iter_mut() {
                *v = 0.0;
            }
            // x faces: i = 1..nx-1 carry flux (0 and nx are walls)
            for j in 0..ny {
                for i in 1..nx {
                    let face = ux[j * (nx + 1) + i];
                    let l = j * nx + (i - 1);
                    let r = j * nx + i;
                    let sl = if i >= 2 { hc[l] - hc[l - 1] } else { 0.0 };
                    let sr = hc[r] - hc[l];
                    let hl = hc[l] + 0.5 * minmod(sl, sr);
                    let sr2 = if i + 1 < nx { hc[r + 1] - hc[r] } else { 0.0 };
                    let hr = hc[r] + 0.5 * minmod(sr, sr2);
                    let hup = if face >= 0.0 { hl } else { hr };
                    // enthalpy flux = mass flux (ρu) × specific enthalpy
                    let f = 0.5 * (self.rho[l] + self.rho[r]) * face * hup;
                    dudt[l] -= f / dx;
                    dudt[r] += f / dx;
                }
            }
            // y faces: j = 1..ny-1 carry flux
            for j in 1..ny {
                for i in 0..nx {
                    let face = uy[j * nx + i];
                    let b = (j - 1) * nx + i;
                    let t = j * nx + i;
                    let sb = if j >= 2 { hc[b] - hc[b - nx] } else { 0.0 };
                    let st = hc[t] - hc[b];
                    let hb = hc[b] + 0.5 * minmod(sb, st);
                    let st2 = if j + 1 < ny { hc[t + nx] - hc[t] } else { 0.0 };
                    let ht = hc[t] + 0.5 * minmod(st, st2);
                    let hup = if face >= 0.0 { hb } else { ht };
                    let f = 0.5 * (self.rho[b] + self.rho[t]) * face * hup;
                    dudt[b] -= f / dy;
                    dudt[t] += f / dy;
                }
            }
            for c in 0..n {
                self.u[c] += h * dudt[c];
            }
        }
    }
}

#[inline]
fn minmod(a: f64, b: f64) -> f64 {
    if a * b <= 0.0 {
        0.0
    } else if a.abs() < b.abs() {
        a
    } else {
        b
    }
}

#[inline]
fn harm(a: f64, b: f64) -> f64 {
    if a <= 0.0 || b <= 0.0 {
        0.0
    } else {
        2.0 * a * b / (a + b)
    }
}

// ---------------------------------------------------------------- indicators

/// Critical Rayleigh number for the onset of steady Rayleigh–Bénard
/// convection between rigid, isothermal boundaries. Drazin & Reid (1981),
/// Chandrasekhar (1961): `Ra_c = 1707.762`. Same constant as
/// `photothermal::indicators::RA_CRIT_RAYLEIGH_BENARD`.
pub const RA_CRIT_RAYLEIGH_BENARD: f64 = 1707.762;

/// Critical internal Rayleigh number for uniform volumetric heating of a
/// rigid, isothermal layer. Sparrow, Goldstein & Jonsson, J. Fluid Mech. 18,
/// 513 (1964): `Ra_I,c = 2772`.
pub const RA_CRIT_INTERNAL: f64 = 2772.0;

/// Classical Rayleigh number `Ra = g β ΔT L³/(ν α)`.
pub fn rayleigh_benard(g: f64, beta: f64, delta_t: f64, l: f64, nu: f64, alpha: f64) -> f64 {
    g * beta * delta_t * l * l * l / (nu * alpha)
}

/// Internal-heating Rayleigh number `Ra_I = g β q L⁵/(k ν α)`.
pub fn rayleigh_internal(g: f64, beta: f64, q: f64, l: f64, k: f64, nu: f64, alpha: f64) -> f64 {
    g * beta * q * l.powi(5) / (k * nu * alpha)
}

/// Does the conduction-only isobaric rung hold, or must it unfold to flow?
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ConvectionRegime {
    Conduction,
    Convection,
}

/// Decide from the surface-like Rayleigh number.
pub fn convection_regime(ra: f64) -> ConvectionRegime {
    if ra > RA_CRIT_RAYLEIGH_BENARD {
        ConvectionRegime::Convection
    } else {
        ConvectionRegime::Conduction
    }
}
