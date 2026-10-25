//! Certified resolution ladder on the emergence engine.
//!
//! Each rung runs the pseudo-spectral solver at one resolution `N` with the
//! rigorous weak–strong certificate ([`crate::certify::RigorousCertifier`]):
//! a proven bound on `‖u − v‖` (box-mean `L²`) between the computed field `v`
//! and *any* Leray–Hopf solution `u` with the same data. The claim is
//! therefore [`Rigor::Rigorous`] and known only after the run, so every rung
//! is [`Assessment::Deferred`]; the engine tries resolutions in order of
//! cost and unfolds to the next one exactly when a certificate does not meet
//! the tolerance.
//!
//! The observable is the whole final velocity field; [`FlowResult::energy`]
//! carries the kinetic energy and [`FlowResult::energy_bound`] its rigorous
//! error, `|½‖u‖² − ½‖v‖²| ≤ δ ‖v‖ + δ²/2` for `δ` the `L²` bound.

use emergence::{Assessment, Audit, Claim, Ladder, Rung};
use num_complex::Complex64;

use crate::certify::RigorousCertifier;
use crate::grid::wavenumber;
use crate::{init, Grid, Solver};

/// Initial data exactly representable on every grid, so the initial error
/// is zero.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Initial {
    TaylorGreen,
    /// Arnold–Beltrami–Childress flow `(a, b, c)` at wavenumber 1 (an exact
    /// steady Euler solution, decaying exactly under viscosity).
    Abc(f64, f64, f64),
}

/// A periodic `[0, 2π)³` flow problem.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Flow {
    pub initial: Initial,
    pub nu: f64,
    pub t_end: f64,
    pub dt: f64,
}

impl Flow {
    pub fn steps(&self) -> usize {
        (self.t_end / self.dt).round().max(1.0) as usize
    }
}

/// A computed final state.
#[derive(Clone, Debug)]
pub struct FlowResult {
    pub grid: Grid,
    pub u: [Vec<Complex64>; 3],
    /// Rigorous `L²` bound at `t_end` (∞ if the certificate was not run).
    pub l2_bound: f64,
}

impl FlowResult {
    /// `‖v‖` (box-mean `L²`).
    pub fn norm(&self) -> f64 {
        let g = &self.grid;
        (0..g.spec_len())
            .map(|i| g.weight(i) * (0..3).map(|c| self.u[c][i].norm_sqr()).sum::<f64>())
            .sum::<f64>()
            .sqrt()
    }
    /// Kinetic energy `½‖v‖²`.
    pub fn energy(&self) -> f64 {
        0.5 * self.norm().powi(2)
    }
    /// Rigorous bound on the energy error, `δ‖v‖ + δ²/2`.
    pub fn energy_bound(&self) -> f64 {
        let d = self.l2_bound;
        d * self.norm() + 0.5 * d * d
    }
}

/// `L²` distance between two fields on (possibly) different grids: the
/// coarser spectrum is embedded in the finer one by wavenumber.
pub fn l2_distance(a: &FlowResult, b: &FlowResult) -> f64 {
    let (c, f) = if a.grid.n <= b.grid.n { (a, b) } else { (b, a) };
    let (gc, gf) = (&c.grid, &f.grid);
    let mut map = vec![None; gf.spec_len()];
    for idx in 0..gc.spec_len() {
        let (i, j, l) = gc.spec_ijl(idx);
        let w = |k: i64| k.rem_euclid(gf.n as i64) as usize;
        map[(w(wavenumber(i, gc.n)) * gf.n + w(wavenumber(j, gc.n))) * gf.nh + l] = Some(idx);
    }
    let mut e2 = 0.0;
    for (fi, m) in map.iter().enumerate() {
        for comp in 0..3 {
            let coarse = m.map(|ci| c.u[comp][ci]).unwrap_or_default();
            e2 += gf.weight(fi) * (f.u[comp][fi] - coarse).norm_sqr();
        }
    }
    e2.sqrt()
}

fn initial_field(g: &Grid, init: Initial) -> [Vec<f64>; 3] {
    match init {
        Initial::TaylorGreen => init::taylor_green(g),
        Initial::Abc(a, b, c) => init::abc(g, a, b, c, 1.0),
    }
}

/// Run `flow` at resolution `n`, with or without the rigorous certificate.
pub fn run(flow: &Flow, n: usize, certify: bool) -> FlowResult {
    let g = Grid::new(n);
    let mut s = Solver::new(g.clone(), flow.nu, initial_field(&g, flow.initial));
    let mut l2_bound = f64::INFINITY;
    if certify {
        let mut rc = RigorousCertifier::new(&g, 0.0, 0.0);
        for _ in 0..flow.steps() {
            s.step_certified_rigorous(flow.dt, &mut rc);
        }
        l2_bound = rc.history.last().map_or(f64::INFINITY, |st| st.l2_bound);
        if !l2_bound.is_finite() {
            l2_bound = f64::INFINITY;
        }
    } else {
        for _ in 0..flow.steps() {
            s.step(flow.dt);
        }
    }
    FlowResult { grid: g, u: s.u, l2_bound }
}

/// One certified resolution.
pub struct Resolution {
    pub n: usize,
    name: String,
}

impl Resolution {
    pub fn new(n: usize) -> Self {
        Self {
            n,
            name: format!("spectral {n}³ certified"),
        }
    }
}

impl Rung<Flow, FlowResult> for Resolution {
    fn name(&self) -> &str {
        &self.name
    }
    /// `N³ log₂ N` per step (FFT work; the certificate multiplies every
    /// resolution by the same factor). Units are relative.
    fn cost(&self, flow: &Flow) -> f64 {
        let n = self.n as f64;
        n.powi(3) * n.log2() * flow.steps() as f64
    }
    fn assess(&self, _: &Flow) -> Assessment {
        Assessment::Deferred
    }
    fn eval(&self, flow: &Flow) -> FlowResult {
        run(flow, self.n, true)
    }
    fn certify(&self, _: &Flow, out: &FlowResult) -> Option<Claim> {
        Some(Claim::rigorous(out.l2_bound))
    }
}

/// The ladder over the given resolutions (the finest is the audit
/// reference).
pub fn ladder(resolutions: &[usize], audit_every: u64) -> Ladder<Flow, FlowResult> {
    let mut l = Ladder::new("certified spectral NS", l2_distance);
    let mut last = 0;
    for &n in resolutions {
        last = l.push(Resolution::new(n));
    }
    l.reference(last).audit(if audit_every == 0 { Audit::Off } else { Audit::Every(audit_every) })
}
