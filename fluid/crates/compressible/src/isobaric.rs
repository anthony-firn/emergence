//! Slow phase change at constant pressure: the low-Mach (isobaric) limit
//! of the multiphase flow model, for melting, freezing, boiling and
//! condensation driven by heat conduction.
//!
//! Why a separate rung. When heat diffusion is the slowest process, the
//! acoustic Courant number of a time step that follows it is 10⁵–10⁷, and a
//! stiff condensed phase stored as (ρ, ρu, E) cannot resolve its temperature:
//! at fixed energy T depends on density through the elastic energy
//! (∂T/∂ρ ~ ΓT/ρ), so the O(1/C²) density errors of any compressible update
//! become kelvins. In this limit acoustics and elasticity fold away; what
//! remains is the isobaric enthalpy balance
//!   ρ Dh/Dt = ∇·(k∇T),   p = p₀,
//! with density following the equation of state, ρ = 1/Σ y_k v_k(p₀, T),
//! and the phase mass fractions y_k given by phase equilibrium (or by
//! nucleation and growth kinetics). The phases are the same NASG phases as
//! everywhere else (`elements`, `substances`).
//!
//! Mesh: 1D Lagrangian (mass coordinates): cells keep their mass, their
//! width is m v(h, y), so expansion on melting or condensation of a vapour
//! moves material exactly; a wall at x = 0 is fixed. With the new phase
//! adjacent to the wall the exact Neumann solution holds even with a density
//! change, using the new phase's properties (`examples/phase_change`).

use crate::elements::ThreePhase;
use crate::eos::Nasg;
use crate::phase_change::{e_pt, v_pt};
use crate::thermal_eos::ThermalEos;

/// Most phases a substance may have (solids, liquid, vapour).
pub const MAXP: usize = 6;
pub type Fr = [f64; MAXP];

const R_GAS: f64 = 8.314_462_618;

/// Inden / Hillert–Jarl magnetic contribution to the Gibbs energy, the
/// CALPHAD model of Dinsdale (Calphad 15, 1991). For τ = T/T_C,
///   g_mag = R T ln(β+1) f(τ),
/// and the enthalpy and heat capacity follow from h = g − T ∂g/∂T. This is
/// the temperature-dependent c_p that a NASG phase cannot carry: the
/// magnetic ordering of bcc iron decays above T_C = 1043 K and its heat
/// capacity falls with T, unlike the constant c_p·γ of NASG.
#[derive(Clone, Copy, Debug)]
pub struct Magnetic {
    /// Curie temperature T_C, K.
    pub tc: f64,
    /// structure parameter p (0.4 bcc, 0.28 fcc).
    pub p: f64,
    /// mean magnetic moment per atom, Bohr magnetons.
    pub beta: f64,
    /// molar mass, kg/mol (the model is per mole; NASG carries J/kg).
    pub mm: f64,
}

impl Magnetic {
    /// f(τ), f′(τ), f″(τ): the A-normalised SGTE polynomial. The piecewise
    /// form is continuous in value and slope at τ = 1 (the small mismatch
    /// from the rounded CALPHAD constants is the usual one).
    fn poly(&self, tau: f64) -> (f64, f64, f64) {
        let p = self.p;
        let a = 518.0 / 1125.0 + (11692.0 / 15975.0) * (1.0 / p - 1.0);
        if tau <= 1.0 {
            let inv = 1.0 / tau;
            let t1 = 79.0 / (140.0 * p) * inv;
            let t1p = -79.0 / (140.0 * p) * inv * inv;
            let t1pp = 2.0 * 79.0 / (140.0 * p) * inv * inv * inv;
            let c = (474.0 / 497.0) * (1.0 / p - 1.0);
            let u = tau * tau * tau / 6.0 + tau.powi(9) / 135.0 + tau.powi(15) / 600.0;
            let up = tau * tau / 2.0 + tau.powi(8) / 15.0 + tau.powi(14) / 40.0;
            let upp = tau + 8.0 * tau.powi(7) / 15.0 + 7.0 * tau.powi(13) / 20.0;
            (1.0 - (t1 + c * u) / a, -(t1p + c * up) / a, -(t1pp + c * upp) / a)
        } else {
            let inv = 1.0 / tau;
            let t3 = inv.powi(5) / 10.0 + inv.powi(15) / 315.0 + inv.powi(25) / 1500.0;
            let t3p = -0.5 * inv.powi(6) - inv.powi(16) / 21.0 - inv.powi(26) / 60.0;
            let t3pp = 3.0 * inv.powi(7) + 16.0 * inv.powi(17) / 21.0 + 13.0 * inv.powi(27) / 30.0;
            (-t3 / a, -t3p / a, -t3pp / a)
        }
    }
    /// Magnetic Gibbs energy per unit mass, J/kg.
    pub fn gibbs(&self, t: f64) -> f64 {
        let tau = t / self.tc;
        R_GAS * t * (self.beta + 1.0).ln() * self.poly(tau).0 / self.mm
    }
    /// Magnetic enthalpy per unit mass, J/kg (h = g − T ∂g/∂T).
    pub fn enthalpy(&self, t: f64) -> f64 {
        let tau = t / self.tc;
        -R_GAS * (self.beta + 1.0).ln() * t * tau * self.poly(tau).1 / self.mm
    }
    /// Magnetic heat capacity per unit mass, J/(kg K).
    pub fn cp(&self, t: f64) -> f64 {
        let tau = t / self.tc;
        let (_, fp, fpp) = self.poly(tau);
        -R_GAS * (self.beta + 1.0).ln() * (2.0 * tau * fp + tau * tau * fpp) / self.mm
    }
}

/// The phases of one substance: several solids (allotropes) allowed.
#[derive(Clone, Debug)]
pub struct PhaseSet {
    pub names: Vec<String>,
    pub eos: Vec<Nasg>,
    /// index of the liquid (the freezing kinetics' parent phase)
    pub liquid: usize,
    /// index of the vapour
    pub vapor: usize,
    /// above the critical temperature the liquid is not a separate phase
    pub t_crit: Option<f64>,
    /// magnetic contribution of each phase (None: non-magnetic)
    pub magnetic: Vec<Option<Magnetic>>,
    /// optional temperature-dependent solid EOS overriding [`Self::eos`]
    /// (iron allotropes); None keeps the NASG phase bit-identical
    pub thermal: Vec<Option<ThermalEos>>,
}

impl From<&ThreePhase> for PhaseSet {
    fn from(m: &ThreePhase) -> Self {
        let n = 3;
        PhaseSet {
            names: vec!["solid".into(), "liquid".into(), "vapour".into()],
            eos: vec![m.solid, m.liquid, m.vapor],
            liquid: 1,
            vapor: 2,
            t_crit: m.t_crit,
            magnetic: vec![None; n],
            thermal: vec![None; n],
        }
    }
}

impl PhaseSet {
    pub fn len(&self) -> usize {
        self.eos.len()
    }
    pub fn is_empty(&self) -> bool {
        self.eos.is_empty()
    }
    /// Gibbs energy of phase k at (p, T), including any magnetic term. A
    /// `thermal` EOS, when present, replaces the NASG phase (iron); it is
    /// per mole and divided by its molar mass.
    pub fn gibbs(&self, k: usize, p: f64, t: f64) -> f64 {
        let base = match &self.thermal[k] {
            Some(e) => e.gibbs(p, t) / e.mm,
            None => self.eos[k].gibbs(p, t),
        };
        base + self.magnetic[k].map_or(0.0, |m| m.gibbs(t))
    }
    /// Specific enthalpy of phase k at (p, T), including any magnetic term.
    pub fn enthalpy(&self, k: usize, p: f64, t: f64) -> f64 {
        let base = match &self.thermal[k] {
            Some(e) => e.enthalpy(p, t) / e.mm,
            None => enthalpy(&self.eos[k], p, t),
        };
        base + self.magnetic[k].map_or(0.0, |m| m.enthalpy(t))
    }
    /// Specific heat capacity c_p of phase k at (p, T), including any
    /// magnetic term.
    pub fn cp(&self, k: usize, p: f64, t: f64) -> f64 {
        let base = match &self.thermal[k] {
            Some(e) => e.cp(p, t) / e.mm,
            None => self.eos[k].gamma * self.eos[k].cv,
        };
        base + self.magnetic[k].map_or(0.0, |m| m.cp(t))
    }
    /// Specific volume of phase k at (p, T): the thermal EOS when present,
    /// otherwise the NASG phase.
    pub fn specific_volume(&self, k: usize, p: f64, t: f64) -> f64 {
        match &self.thermal[k] {
            Some(e) => e.volume(p, t) / e.mm,
            None => v_pt(&self.eos[k], p, t),
        }
    }
    /// Stable phase at (p, T): least Gibbs energy (no separate liquid above
    /// the critical temperature).
    pub fn stable(&self, p: f64, t: f64) -> usize {
        let mut best = (f64::INFINITY, 0);
        for k in 0..self.eos.len() {
            if k == self.liquid && self.t_crit.is_some_and(|tc| t > tc) {
                continue;
            }
            // A thermal solid EOS is only valid up to `t_max` at 0.1 MPa;
            // above it the phase is omitted (the published EOS turns over).
            if self.thermal[k].is_some_and(|e| t > e.t_max) {
                continue;
            }
            let g = self.gibbs(k, p, t);
            if g < best.0 {
                best = (g, k);
            }
        }
        best.1
    }
    /// Is phase k a solid?
    pub fn is_solid(&self, k: usize) -> bool {
        k != self.liquid && k != self.vapor
    }
}

/// Specific enthalpy h(p, T) of a NASG phase.
pub fn enthalpy(e: &Nasg, p: f64, t: f64) -> f64 {
    e_pt(e, p, t) + p * v_pt(e, p, t)
}

/// Phase transitions at fixed pressure, in increasing temperature.
#[derive(Clone, Debug)]
pub struct Ladder {
    pub p: f64,
    /// (T, phase below, phase above)
    pub transitions: Vec<(f64, usize, usize)>,
    /// stable phase below the first transition
    pub lowest: usize,
}

impl Ladder {
    /// Scan T in [t_lo, t_hi] for changes of the stable phase and refine
    /// each by bisection on the Gibbs energies.
    pub fn new(model: &PhaseSet, p: f64, t_lo: f64, t_hi: f64) -> Self {
        let n = 4000;
        let stable = |t: f64| model.stable(p, t);
        let lowest = stable(t_lo);
        let mut transitions = Vec::new();
        let mut prev = (t_lo, lowest);
        for i in 1..=n {
            let t = t_lo * (t_hi / t_lo).powf(i as f64 / n as f64);
            let s = stable(t);
            if s != prev.1 {
                let (mut a, mut b) = (prev.0, t);
                for _ in 0..80 {
                    let mid = 0.5 * (a + b);
                    if stable(mid) == prev.1 {
                        a = mid;
                    } else {
                        b = mid;
                    }
                }
                transitions.push((0.5 * (a + b), prev.1, s));
                prev = (t, s);
            }
        }
        Ladder { p, transitions, lowest }
    }
}

/// Homogeneous freezing kinetics of a supercooled liquid: nucleation rate
/// J(T) (m⁻³ s⁻¹) and crystal growth velocity v(ΔT) (m/s).
pub struct Freezing {
    pub rate: Box<dyn Fn(f64) -> f64 + Send + Sync>,
    pub growth: Box<dyn Fn(f64) -> f64 + Send + Sync>,
}

pub struct Isobaric {
    pub model: PhaseSet,
    pub ladder: Ladder,
    /// conductivity per phase (order of the phase set), W/(m K)
    pub k: Vec<f64>,
    /// mass per unit area of each cell (kg/m²)
    pub m: Vec<f64>,
    /// specific enthalpy (J/kg) and phase mass fractions
    pub h: Vec<f64>,
    pub y: Vec<Fr>,
    /// wall temperatures at x = 0 and at the far end (None: adiabatic)
    pub wall: [Option<f64>; 2],
    /// heat removed per unit mass and time from every cell (J/(kg s));
    /// the 0D cooling of a droplet
    pub sink: f64,
    /// finite-rate freezing (None: phase equilibrium at all times)
    pub freezing: Option<Freezing>,
    /// nuclei per unit mass of each cell (freezing kinetics)
    pub nuclei: Vec<f64>,
    pub t: f64,
}

impl Isobaric {
    /// Cells of width `dx` (at the initial state) with temperature t0 in
    /// phase `ph` (pure).
    #[allow(clippy::too_many_arguments)]
    pub fn new(model: PhaseSet, p: f64, k: Vec<f64>, n: usize, dx: f64, ph: usize, t0: f64, t_range: (f64, f64)) -> Self {
        let ladder = Ladder::new(&model, p, t_range.0, t_range.1);
        let rho = 1.0 / model.specific_volume(ph, p, t0);
        let h0 = model.enthalpy(ph, p, t0);
        let mut y = [0.0; MAXP];
        y[ph] = 1.0;
        Isobaric {
            model,
            ladder,
            k,
            m: vec![rho * dx; n],
            h: vec![h0; n],
            y: vec![y; n],
            wall: [None, None],
            sink: 0.0,
            freezing: None,
            nuclei: vec![0.0; n],
            t: 0.0,
        }
    }

    pub fn eos(&self) -> &[Nasg] {
        &self.model.eos
    }
    fn np(&self) -> usize {
        self.model.len()
    }

    /// Temperature of enthalpy h with fixed mass fractions y (Newton on
    /// Σ y_k h_k(T) = h).
    pub fn temperature(&self, h: f64, y: &Fr) -> f64 {
        let p = self.ladder.p;
        let np = self.np();
        let hy = |t: f64| (0..np).map(|k| y[k] * self.model.enthalpy(k, p, t)).sum::<f64>();
        let cp = |t: f64| (0..np).map(|k| y[k] * (self.model.enthalpy(k, p, t * (1.0 + 1e-7)) - self.model.enthalpy(k, p, t)) / (1e-7 * t)).sum::<f64>();
        let mut t = self.ladder.transitions.first().map_or(300.0, |x| x.0);
        for _ in 0..100 {
            let dt = (hy(t) - h) / cp(t);
            let tn = (t - dt).clamp(0.2 * t, 5.0 * t);
            if (tn - t).abs() < 1e-12 * t {
                return tn;
            }
            t = tn;
        }
        t
    }

    /// Equilibrium state (T, y) of enthalpy h at the ladder pressure.
    pub fn equilibrium(&self, h: f64) -> (f64, Fr) {
        let p = self.ladder.p;
        let mut cur = self.ladder.lowest;
        for &(tt, a, b) in &self.ladder.transitions {
            let (ha, hb) = (self.model.enthalpy(a, p, tt), self.model.enthalpy(b, p, tt));
            if h < ha {
                break;
            }
            if h <= hb {
                let mut y = [0.0; MAXP];
                y[b] = (h - ha) / (hb - ha);
                y[a] = 1.0 - y[b];
                return (tt, y);
            }
            cur = b;
        }
        let mut y = [0.0; MAXP];
        y[cur] = 1.0;
        (self.temperature(h, &y), y)
    }

    /// Specific volume of a cell.
    pub fn volume(&self, i: usize) -> f64 {
        let t = self.temperature(self.h[i], &self.y[i]);
        (0..self.np()).map(|k| self.y[i][k] * self.model.specific_volume(k, self.ladder.p, t)).sum()
    }

    /// Cell widths and centres (wall at x = 0).
    pub fn geometry(&self) -> (Vec<f64>, Vec<f64>) {
        let w: Vec<f64> = (0..self.m.len()).map(|i| self.m[i] * self.volume(i)).collect();
        let mut x = Vec::with_capacity(w.len());
        let mut acc = 0.0;
        for &wi in &w {
            x.push(acc + 0.5 * wi);
            acc += wi;
        }
        (w, x)
    }

    fn conductivity(&self, i: usize, t: f64) -> f64 {
        let p = self.ladder.p;
        let np = self.np();
        let vk: Vec<f64> = (0..np).map(|k| self.y[i][k] * self.model.specific_volume(k, p, t)).collect();
        let v: f64 = vk.iter().sum();
        (0..np).map(|k| vk[k] / v * self.k[k]).sum()
    }

    /// Advance by `dt` in stable explicit sub-steps; returns their number.
    pub fn advance(&mut self, dt: f64) -> usize {
        let n = self.m.len();
        let mut done = 0.0;
        let mut steps = 0;
        while done < dt * (1.0 - 1e-12) {
            let (w, _) = self.geometry();
            let t: Vec<f64> = (0..n).map(|i| self.temperature(self.h[i], &self.y[i])).collect();
            let k: Vec<f64> = (0..n).map(|i| self.conductivity(i, t[i])).collect();
            let (rl, rr) = self.resistances(&w, &t);
            // stable step: each cell's heat capacity (single-phase c_p)
            // against the conductances of its two faces
            let mut lim = f64::INFINITY;
            for i in 0..n {
                let cp_min = (0..self.np()).filter(|&j| self.y[i][j] > 1e-12).map(|j| self.model.cp(j, self.ladder.p, t[i])).fold(f64::INFINITY, f64::min);
                let gl = if i > 0 { 1.0 / (rr[i - 1] + rl[i]) } else if self.wall[0].is_some() { 1.0 / rl[i] } else { 0.0 };
                let gr = if i + 1 < n { 1.0 / (rr[i] + rl[i + 1]) } else if self.wall[1].is_some() { 1.0 / rr[i] } else { 0.0 };
                if gl + gr > 0.0 {
                    lim = lim.min(0.45 * self.m[i] * cp_min / (gl + gr));
                }
            }
            let h = lim.min(dt - done);
            // thermal resistance from each cell's temperature point to its
            // left and right faces. Single phase: the centre, half a width
            // through the cell. Two phases (at the transition temperature):
            // the phases are layered and the interface — the temperature
            // point — sits between them; each side conducts through its own
            // phase only (a condensing film conducts through the liquid, not
            // through the isothermal saturated vapour beside it).
            let mut q = vec![0.0; n + 1];
            if let Some(tw) = self.wall[0] {
                q[0] = (tw - t[0]) / rl[0];
            }
            if let Some(tw) = self.wall[1] {
                q[n] = -(tw - t[n - 1]) / rr[n - 1];
            }
            for f in 1..n {
                q[f] = (t[f - 1] - t[f]) / (rr[f - 1] + rl[f]);
            }
            let _ = &k;
            for i in 0..n {
                self.h[i] += h * ((q[i] - q[i + 1]) / self.m[i] - self.sink);
            }
            self.relax(h);
            done += h;
            self.t += h;
            steps += 1;
        }
        steps
    }

    /// Left and right thermal resistances (per unit area) of every cell.
    fn resistances(&self, w: &[f64], t: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let p = self.ladder.p;
        let n = self.m.len();
        let mut rl = vec![0.0; n];
        let mut rr = vec![0.0; n];
        for i in 0..n {
            let present: Vec<usize> = (0..self.np()).filter(|&k| self.y[i][k] > 1e-12).collect();
            if present.len() == 2 {
                let (a, b) = (present[0], present[1]);
                let va = self.y[i][a] * self.model.specific_volume(a, p, t[i]);
                let vb = self.y[i][b] * self.model.specific_volume(b, p, t[i]);
                let (fa, fb) = (va / (va + vb), vb / (va + vb));
                // which phase lies on the left: the one the left neighbour
                // (or a wall) holds more of
                let share = |j: Option<usize>, k: usize| -> f64 { j.map_or(0.5, |j| self.y[j][k]) };
                let left = if i > 0 { Some(i - 1) } else { None };
                let right = if i + 1 < n { Some(i + 1) } else { None };
                let a_left = if left.is_none() && self.wall[0].is_some() {
                    // the wall side holds the phase stable at the wall temperature
                    self.equilibrium_phase_at(self.wall[0].unwrap()) == a
                } else {
                    share(left, a) - share(right, a) >= 0.0
                };
                let (l_phase, l_frac, r_phase, r_frac) = if a_left { (a, fa, b, fb) } else { (b, fb, a, fa) };
                rl[i] = l_frac * w[i] / self.k[l_phase];
                rr[i] = r_frac * w[i] / self.k[r_phase];
            } else {
                let kc = self.conductivity(i, t[i]);
                rl[i] = 0.5 * w[i] / kc;
                rr[i] = 0.5 * w[i] / kc;
            }
        }
        (rl, rr)
    }

    /// Stable phase at temperature t on the ladder.
    fn equilibrium_phase_at(&self, t: f64) -> usize {
        let mut cur = self.ladder.lowest;
        for &(tt, _, b) in &self.ladder.transitions {
            if t > tt {
                cur = b;
            }
        }
        cur
    }

    /// Phase relaxation: equilibrium, or finite-rate freezing of a
    /// supercooled liquid (nucleation, then growth until equilibrium).
    fn relax(&mut self, dt: f64) {
        let n = self.m.len();
        for i in 0..n {
            let (teq, yeq) = self.equilibrium(self.h[i]);
            let Some(fr) = &self.freezing else {
                self.y[i] = yeq;
                continue;
            };
            // freezing of the liquid is the only kinetically limited
            // transition here: everything else relaxes to equilibrium
            let lq = self.model.liquid;
            let Some(&(t_fr, sol, _)) = self.ladder.transitions.iter().find(|x| x.2 == lq && self.model.is_solid(x.1)) else {
                self.y[i] = yeq;
                continue;
            };
            let liquid_to_solid = yeq[sol] > self.y[i][sol] + 1e-15 && self.y[i][lq] > 0.0;
            if !liquid_to_solid {
                self.y[i] = yeq;
                continue;
            }
            let t = self.temperature(self.h[i], &self.y[i]);
            let _ = teq;
            let undercool = (t_fr - t).max(0.0);
            let v = self.volume(i);
            // nucleation (per unit mass)
            self.nuclei[i] += (fr.rate)(t) * v * dt;
            let nb = self.nuclei[i] / v; // crystals per unit volume
            if self.nuclei[i] * self.m[i] < 1.0 && self.y[i][sol] <= 0.0 {
                // fewer than one crystal in the cell: still liquid (the cell
                // mass per area times one unit area — 0D: the droplet mass)
                continue;
            }
            // Avrami growth of nb crystals per volume at velocity v
            let x = self.y[i][sol];
            let x_ext = -(1.0 - x).max(1e-300).ln();
            let gv = (fr.growth)(undercool);
            let dx = (1.0 - x) * (36.0 * std::f64::consts::PI * nb).cbrt() * x_ext.max(1e-12).powf(2.0 / 3.0) * gv * dt;
            // a fresh crystal population starts at the critical size scale
            let dx = if x <= 0.0 { (1e-9_f64).max(dx) } else { dx };
            let new_s = (x + dx).min(yeq[sol]);
            let mut ynew = self.y[i];
            let take = new_s - x;
            ynew[sol] = new_s;
            ynew[lq] = (ynew[lq] - take).max(0.0);
            self.y[i] = ynew;
        }
    }

    /// Mean temperature and phase fractions (mass-weighted).
    pub fn temperatures(&self) -> Vec<f64> {
        (0..self.m.len()).map(|i| self.temperature(self.h[i], &self.y[i])).collect()
    }
}
