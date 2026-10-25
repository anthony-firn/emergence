//! Real-fluid equilibrium equation of state for the flow solver, built on
//! a reference Helmholtz formulation and tabulated for speed.
//!
//! The state (ρ, e) is either single-phase (liquid, vapour, supercritical:
//! straight from the Helmholtz formulation) or inside the vapour dome,
//! where the equilibrium fluid is a liquid–vapour mixture at the
//! saturation temperature T_s with quality x (the homogeneous equilibrium
//! model, so cavitation and flashing follow from the EOS alone). The
//! equilibrium sound speed is c² = (∂p/∂ρ)_s = p_ρ|_e + (p/ρ²) p_e|_ρ.
//!
//! Emergent efficiency: exact evaluation needs Newton iterations on the
//! Helmholtz energy (and a saturation solve inside the dome), far too slow
//! for every flux. Both regions are smooth on their own, so each gets a
//! bicubic table over (ln ρ, e), extended a little across the dome
//! boundary; the boundary itself is a 1-D table e_sat(ρ). At run time a
//! state picks its region from e_sat and interpolates — no exact
//! evaluations in the time loop. `RealFluid::error_report` measures the
//! interpolation error against the exact evaluation.

use crate::helmholtz::Helmholtz;
use rayon::prelude::*;

#[derive(Clone, Copy, Debug)]
pub struct FluidState {
    pub p: f64,
    pub t: f64,
    pub c: f64,
    /// vapour mass fraction inside the dome, NaN outside
    pub x: f64,
}

/// Saturation curve sampled in T (dense; cubic interpolation).
#[derive(Clone, Debug)]
struct SatTable {
    t0: f64,
    dt: f64,
    p: Vec<f64>,
    rl: Vec<f64>,
    rv: Vec<f64>,
    el: Vec<f64>,
    ev: Vec<f64>,
}

fn cubic(v: &[f64], s: f64) -> f64 {
    // Catmull–Rom on a uniform grid; s in index units
    let n = v.len();
    let i = (s.floor() as isize).clamp(1, n as isize - 3) as usize;
    let f = s - i as f64;
    let (p0, p1, p2, p3) = (v[i - 1], v[i], v[i + 1], v[i + 2]);
    p1 + 0.5 * f * (p2 - p0 + f * (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3 + f * (3.0 * (p1 - p2) + p3 - p0)))
}

impl SatTable {
    fn build(h: &Helmholtz, n: usize) -> Self {
        let t0 = h.t_triple;
        let t1 = h.t_c * (1.0 - 2e-4);
        let dt = (t1 - t0) / (n - 1) as f64;
        let rows: Vec<(f64, f64, f64, f64, f64)> = (0..n)
            .into_par_iter()
            .map(|i| {
                let t = t0 + dt * i as f64;
                let (p, rl, rv) = h.saturation(t).unwrap_or((f64::NAN, f64::NAN, f64::NAN));
                (p, rl, rv, h.props(rl, t).u, h.props(rv, t).u)
            })
            .collect();
        SatTable {
            t0,
            dt,
            p: rows.iter().map(|r| r.0).collect(),
            rl: rows.iter().map(|r| r.1).collect(),
            rv: rows.iter().map(|r| r.2).collect(),
            el: rows.iter().map(|r| r.3).collect(),
            ev: rows.iter().map(|r| r.4).collect(),
        }
    }
    fn tmax(&self) -> f64 {
        self.t0 + self.dt * (self.p.len() - 1) as f64
    }
    fn at(&self, t: f64) -> (f64, f64, f64, f64, f64) {
        let s = (t - self.t0) / self.dt;
        (cubic(&self.p, s), cubic(&self.rl, s), cubic(&self.rv, s), cubic(&self.el, s), cubic(&self.ev, s))
    }
}

/// Exact (slow) evaluations used to build and check the tables.
pub struct Exact<'a> {
    pub h: &'a Helmholtz,
    sat: &'a SatTable,
}

impl<'a> Exact<'a> {
    /// Single-phase state from the Helmholtz energy (metastable inside the
    /// dome near its boundary).
    pub fn single(&self, rho: f64, e: f64, t_guess: f64) -> Option<FluidState> {
        let _ = t_guess;
        let t = self.h.temperature_stable(rho, e)?;
        let pr = self.h.props(rho, t);
        // only stable or metastable states: near the critical point the
        // spinodal (∂p/∂ρ|_T = 0) nearly touches the dome boundary, and
        // beyond it the Helmholtz extension is unphysical
        if !(pr.dpdrho_t > 0.0 && pr.cv > 0.0 && pr.w.is_finite() && pr.w > 0.0 && pr.p.is_finite()) {
            return None;
        }
        Some(FluidState { p: pr.p, t, c: pr.w, x: f64::NAN })
    }

    /// Two-phase equilibrium mixture at (ρ, e) (quality may leave [0, 1]
    /// slightly for the table overlap): T_s solves e_l + x(e_v − e_l) = e,
    /// x = (1/ρ − 1/ρ_l)/(1/ρ_v − 1/ρ_l).
    pub fn two_phase_pt(&self, rho: f64, e: f64) -> Option<(f64, f64, f64)> {
        // the state lies on the tie line at T when the cross product
        // (v − v_l)(e_v − e_l) − (e − e_l)(v_v − v_l) vanishes; it changes
        // sign as the tie lines sweep the dome (no division, robust near
        // the critical point where the tie line shrinks)
        let v = 1.0 / rho;
        let h = |t: f64| {
            let (_, rl, rv, el, ev) = self.sat.at(t);
            let (vl, vv) = (1.0 / rl, 1.0 / rv);
            (v - vl) * (ev - el) - (e - el) * (vv - vl)
        };
        let (mut lo, mut hi) = (self.sat.t0 * (1.0 + 1e-9), self.sat.tmax());
        let hlo = h(lo);
        if !(hlo * h(hi) < 0.0) {
            return None;
        }
        for _ in 0..80 {
            let mid = 0.5 * (lo + hi);
            if (h(mid) < 0.0) == (hlo < 0.0) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let t = 0.5 * (lo + hi);
        let (p, rl, rv, ..) = self.sat.at(t);
        let x = (v - 1.0 / rl) / (1.0 / rv - 1.0 / rl);
        Some((p, t, x))
    }

    pub fn two_phase(&self, rho: f64, e: f64) -> Option<FluidState> {
        let (p, t, x) = self.two_phase_pt(rho, e)?;
        // equilibrium sound speed from p(ρ, e) of the mixture
        let (hr, he) = (1e-6 * rho, 1e-6 * e.abs().max(1e3));
        let pr = |r: f64, ee: f64| self.two_phase_pt(r, ee).map(|v| v.0);
        let (pa, pb, pc, pd) = (pr(rho + hr, e)?, pr(rho - hr, e)?, pr(rho, e + he)?, pr(rho, e - he)?);
        let c2 = (pa - pb) / (2.0 * hr) + p / (rho * rho) * (pc - pd) / (2.0 * he);
        Some(FluidState { p, t, c: c2.max(0.0).sqrt(), x })
    }
}

#[derive(Clone, Debug)]
struct Grid {
    p: Vec<f64>,
    t: Vec<f64>,
    c: Vec<f64>,
    x: Vec<f64>,
}

#[derive(Clone, Debug)]
pub struct RealFluid {
    pub h: Helmholtz,
    sat: SatTable,
    /// density axis: φ = ln ρ (vapour side) or φ = ρ (liquid side)
    log_axis: bool,
    lr0: f64,
    dlr: f64,
    nr: usize,
    e0: f64,
    de: f64,
    ne: usize,
    single: Grid,
    two: Grid,
    /// dome boundary e_sat(ln ρ) (NaN where no dome)
    esat: Vec<f64>,
}

impl RealFluid {
    /// Tabulate over ρ ∈ [rho_min, rho_max] (log-spaced) and the energies
    /// spanned by T ∈ [t_min, t_max].
    pub fn build(h: Helmholtz, rho_min: f64, rho_max: f64, t_min: f64, t_max: f64, nr: usize, ne: usize) -> Self {
        let er = Self::energy_range(&h, rho_min, rho_max, t_min, t_max);
        Self::build_axis(h, rho_min, rho_max, er, nr, ne, true)
    }

    /// Energy range from well-behaved states only (the single-phase
    /// Helmholtz energy deep inside the spinodal is meaningless): max over
    /// the t_max isotherm (supercritical), min of the triple-point liquid
    /// and the compressed liquid at (ρ_max, t_min).
    pub fn energy_range(h: &Helmholtz, rho_min: f64, rho_max: f64, t_min: f64, t_max: f64) -> (f64, f64) {
        let mut e1 = f64::NEG_INFINITY;
        for i in 0..400 {
            let rho = (rho_min.ln() + (rho_max / rho_min).ln() * i as f64 / 399.0).exp();
            let u = h.props(rho, t_max).u;
            if u.is_finite() {
                e1 = e1.max(u);
            }
        }
        let tt = h.t_triple;
        let el_t = h.saturation(tt * 1.0001).map(|(_, rl, _)| h.props(rl, tt * 1.0001).u).unwrap_or(f64::INFINITY);
        let e0 = el_t.min(h.props(rho_max, t_min).u);
        (e0, e1)
    }

    /// As `build`, with a log (`true`) or linear density axis.
    #[allow(clippy::too_many_arguments)]
    pub fn build_axis(h: Helmholtz, rho_min: f64, rho_max: f64, e_range: (f64, f64), nr: usize, ne: usize, log_axis: bool) -> Self {
        let sat = SatTable::build(&h, 4000);
        let phi = |r: f64| if log_axis { r.ln() } else { r };
        let (lr0, lr1) = (phi(rho_min), phi(rho_max));
        let dlr = (lr1 - lr0) / (nr - 1) as f64;
        let rho_at = move |i: usize| if log_axis { (lr0 + dlr * i as f64).exp() } else { lr0 + dlr * i as f64 };
        let (e0, e1) = e_range;
        let de = (e1 - e0) / (ne - 1) as f64;
        let ex = Exact { h: &h, sat: &sat };
        let nodes: Vec<(usize, usize)> = (0..nr).flat_map(|i| (0..ne).map(move |j| (i, j))).collect();
        let t_mid = 0.5 * (h.t_triple + h.t_c);
        let vals: Vec<(FluidState, FluidState)> = nodes
            .par_iter()
            .map(|&(i, j)| {
                let rho = rho_at(i);
                let e = e0 + de * j as f64;
                let nan = FluidState { p: f64::NAN, t: f64::NAN, c: f64::NAN, x: f64::NAN };
                let s = ex.single(rho, e, t_mid).unwrap_or(nan);
                let tp = ex.two_phase(rho, e).unwrap_or(nan);
                (s, tp)
            })
            .collect();
        let grid = |sel: fn(&(FluidState, FluidState)) -> FluidState| -> Grid {
            Grid {
                p: vals.iter().map(|v| sel(v).p).collect(),
                t: vals.iter().map(|v| sel(v).t).collect(),
                c: vals.iter().map(|v| sel(v).c).collect(),
                x: vals.iter().map(|v| sel(v).x).collect(),
            }
        };
        let single = grid(|v| v.0);
        let two = grid(|v| v.1);
        // dome boundary: e on the liquid branch (ρ > ρ_c) or vapour branch
        let esat: Vec<f64> = (0..nr)
            .into_par_iter()
            .map(|i| {
                let rho = rho_at(i);
                let n = sat.p.len();
                let (rl_t, rv_t) = (sat.rl[0], sat.rv[0]);
                if !(rho > rv_t && rho < rl_t) {
                    return f64::NAN;
                }
                // densities inside the top of the dome (between the vapour
                // and liquid ends of the last tabulated tie line, just below
                // T_c): the boundary is that tie line
                let (rv_top, rl_top) = (sat.rv[n - 1], sat.rl[n - 1]);
                if rho >= rv_top && rho <= rl_top {
                    let (vv, vl) = (1.0 / rv_top, 1.0 / rl_top);
                    let x = (1.0 / rho - vl) / (vv - vl);
                    return sat.el[n - 1] + x * (sat.ev[n - 1] - sat.el[n - 1]);
                }
                let branch = |k: usize| if rho > h.rho_c { sat.rl[k] } else { sat.rv[k] };
                // find T on the branch with ρ_branch(T) = ρ (monotone in T)
                let (mut lo, mut hi) = (0.0, (n - 2) as f64);
                let g = |s: f64| {
                    let k = s.floor() as usize;
                    let f = s - k as f64;
                    (1.0 - f) * branch(k) + f * branch(k + 1) - rho
                };
                let glo = g(lo);
                if !(glo * g(hi) <= 0.0) {
                    return f64::NAN;
                }
                for _ in 0..60 {
                    let mid = 0.5 * (lo + hi);
                    if (g(mid) < 0.0) == (glo < 0.0) {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                let t = sat.t0 + sat.dt * 0.5 * (lo + hi);
                let (_, rl, rv, el, ev) = sat.at(t);
                // refine: e at the branch density
                let _ = (rl, rv);
                if rho > h.rho_c {
                    el
                } else {
                    ev
                }
            })
            .collect();
        let mut rf = RealFluid { h, sat, log_axis, lr0, dlr, nr, e0, de, ne, single, two, esat };
        // keep each table only within a band across the dome boundary: the
        // single-phase nodes are already restricted to (meta)stable states;
        // the two-phase extension is kept for qualities in [−0.5, 1.5]. The
        // band is wide because the liquid branch of e_sat(ρ) is steep: a
        // bicubic stencil two density cells away must still find values.
        let band = 30.0;
        for i in 0..nr {
            let es = rf.esat[i];
            for j in 0..ne {
                let e = e0 + de * j as f64;
                let k = i * ne + j;
                let inside = es.is_finite() && e < es;
                let depth = if es.is_finite() { (e - es).abs() / de } else { f64::INFINITY };
                if inside && depth > band {
                    rf.single.p[k] = f64::NAN;
                    rf.single.t[k] = f64::NAN;
                    rf.single.c[k] = f64::NAN;
                }
                let xq = rf.two.x[k];
                if (!inside && depth > band) || !(xq > -0.5 && xq < 1.5) {
                    rf.two.p[k] = f64::NAN;
                    rf.two.t[k] = f64::NAN;
                    rf.two.c[k] = f64::NAN;
                    rf.two.x[k] = f64::NAN;
                }
            }
        }
        rf
    }

    #[inline]
    fn coord(&self, rho: f64) -> f64 {
        ((if self.log_axis { rho.ln() } else { rho }) - self.lr0) / self.dlr
    }
    #[inline]
    fn rho_of(&self, s: f64) -> f64 {
        let v = self.lr0 + self.dlr * s;
        if self.log_axis {
            v.exp()
        } else {
            v
        }
    }

    pub fn in_domain(&self, rho: f64, e: f64) -> bool {
        let s = self.coord(rho);
        let q = (e - self.e0) / self.de;
        s >= 1.0 && s <= (self.nr - 3) as f64 && q >= 1.0 && q <= (self.ne - 3) as f64
    }

    /// Is (ρ, e) inside the vapour dome?
    pub fn two_phase_region(&self, rho: f64, e: f64) -> bool {
        let s = self.coord(rho);
        if s < 0.0 || s > (self.nr - 1) as f64 {
            return false;
        }
        let k = (s.floor() as usize).min(self.nr - 2);
        let f = s - k as f64;
        let (a, b) = (self.esat[k], self.esat[k + 1]);
        if !(a.is_finite() && b.is_finite()) {
            return false;
        }
        e < (1.0 - f) * a + f * b
    }

    fn bicubic(&self, g: &[f64], s: f64, q: f64) -> f64 {
        let i = (s.floor() as isize).clamp(1, self.nr as isize - 3) as usize;
        let j = (q.floor() as isize).clamp(1, self.ne as isize - 3) as usize;
        let (fs, fq) = (s - i as f64, q - j as f64);
        let w = |f: f64| -> [f64; 4] {
            [
                0.5 * (-f + 2.0 * f * f - f * f * f),
                0.5 * (2.0 - 5.0 * f * f + 3.0 * f * f * f),
                0.5 * (f + 4.0 * f * f - 3.0 * f * f * f),
                0.5 * (-f * f + f * f * f),
            ]
        };
        let (ws, wq) = (w(fs), w(fq));
        let mut acc = 0.0;
        for a in 0..4 {
            for b in 0..4 {
                acc += ws[a] * wq[b] * g[(i + a - 1) * self.ne + (j + b - 1)];
            }
        }
        acc
    }

    /// Tabulated state at (ρ, e); falls back to the exact evaluation when
    /// the interpolation stencil touches nodes the table cannot vouch for
    /// (outside its band near the dome boundary, or the domain edge).
    pub fn eval(&self, rho: f64, e: f64) -> FluidState {
        let s = self.coord(rho);
        let q = (e - self.e0) / self.de;
        let two = self.two_phase_region(rho, e);
        let g = if two { &self.two } else { &self.single };
        let st = FluidState { p: self.bicubic(&g.p, s, q), t: self.bicubic(&g.t, s, q), c: self.bicubic(&g.c, s, q), x: if two { self.bicubic(&g.x, s, q) } else { f64::NAN } };
        if st.p.is_finite() && st.t.is_finite() && st.c.is_finite() && st.c > 0.0 {
            return st;
        }
        self.exact(rho, e).unwrap_or(st)
    }

    /// Exact evaluation (slow): single-phase or two-phase by region.
    pub fn exact(&self, rho: f64, e: f64) -> Option<FluidState> {
        let ex = Exact { h: &self.h, sat: &self.sat };
        if self.two_phase_region(rho, e) {
            ex.two_phase(rho, e)
        } else {
            ex.single(rho, e, 0.5 * (self.h.t_triple + self.h.t_c))
        }
    }

    /// Maximum relative interpolation error of p, T, c at random states,
    /// against the exact evaluation, split by region.
    pub fn error_report(&self, samples: usize, seed: u64) -> [(usize, [f64; 3]); 2] {
        let mut st = seed;
        let mut rnd = || {
            st ^= st << 13;
            st ^= st >> 7;
            st ^= st << 17;
            (st >> 11) as f64 / (1u64 << 53) as f64
        };
        let pts: Vec<(f64, f64)> = (0..samples)
            .map(|_| {
                let s = 1.0 + rnd() * (self.nr as f64 - 4.0);
                let q = 1.0 + rnd() * (self.ne as f64 - 4.0);
                (self.rho_of(s), self.e0 + self.de * q)
            })
            .collect();
        let res: Vec<Option<(bool, [f64; 3])>> = pts
            .par_iter()
            .map(|&(rho, e)| {
                let ex = self.exact(rho, e)?;
                if !(ex.p.is_finite() && ex.c > 0.0) {
                    return None;
                }
                let tb = self.eval(rho, e);
                let pscale = ex.p.abs().max(1e-3 * rho * ex.c * ex.c);
                if std::env::var("RF_DEBUG").is_ok() && ((tb.p - ex.p).abs() / pscale > 1e-2 || (tb.t / ex.t - 1.0).abs() > 1e-2) {
                    eprintln!("  bad: rho {rho:.4} e {e:.4e} two {}: exact {:?} table {:?}", self.two_phase_region(rho, e), ex, tb);
                }
                Some((self.two_phase_region(rho, e), [(tb.p - ex.p).abs() / pscale, (tb.t / ex.t - 1.0).abs(), (tb.c / ex.c - 1.0).abs()]))
            })
            .collect();
        let mut out = [(0usize, [0.0f64; 3]); 2];
        for (two, e) in res.into_iter().flatten() {
            let k = two as usize;
            out[k].0 += 1;
            for m in 0..3 {
                if e[m].is_finite() {
                    out[k].1[m] = out[k].1[m].max(e[m]);
                }
            }
        }
        out
    }

    /// Saturation pressure at T from the dense table.
    pub fn p_sat(&self, t: f64) -> f64 {
        self.sat.at(t).0
    }
}

/// Real fluid over its whole density range: a log-density table on the
/// vapour side and a linear-density table on the liquid side (a log axis
/// would put ~6% density steps — ~1 GPa — between liquid nodes), split at
/// the critical density with an overlap.
#[derive(Clone, Debug)]
pub struct RealFluidEos {
    pub low: RealFluid,
    pub high: RealFluid,
    pub rho_split: f64,
}

impl RealFluidEos {
    pub fn build(h: Helmholtz, rho_min: f64, rho_max: f64, t_min: f64, t_max: f64, n: usize) -> Self {
        let rs = h.rho_c;
        let er = RealFluid::energy_range(&h, rho_min, rho_max, t_min, t_max);
        let low = RealFluid::build_axis(h.clone(), rho_min, 1.3 * rs, er, n, n, true);
        let high = RealFluid::build_axis(h, 0.7 * rs, rho_max, er, n, n, false);
        RealFluidEos { low, high, rho_split: rs }
    }
    #[inline]
    fn pick(&self, rho: f64) -> &RealFluid {
        if rho < self.rho_split {
            &self.low
        } else {
            &self.high
        }
    }
    pub fn eval(&self, rho: f64, e: f64) -> FluidState {
        self.pick(rho).eval(rho, e)
    }
    /// Inside the tabulated domain (the EOS's validity range here).
    pub fn in_domain(&self, rho: f64, e: f64) -> bool {
        self.pick(rho).in_domain(rho, e)
    }
    pub fn exact(&self, rho: f64, e: f64) -> Option<FluidState> {
        self.pick(rho).exact(rho, e)
    }
    pub fn two_phase_region(&self, rho: f64, e: f64) -> bool {
        self.pick(rho).two_phase_region(rho, e)
    }
    pub fn p_sat(&self, t: f64) -> f64 {
        self.low.p_sat(t)
    }
    pub fn h(&self) -> &Helmholtz {
        &self.low.h
    }
    /// Interpolation errors of both halves: [low, high] × [single, two-phase].
    pub fn error_report(&self, samples: usize) -> [[(usize, [f64; 3]); 2]; 2] {
        [self.low.error_report(samples, 0x1234567), self.high.error_report(samples, 0x7654321)]
    }
}
