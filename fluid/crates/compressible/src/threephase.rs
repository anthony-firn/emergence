//! Solid, liquid and vapour of one substance in one flow.
//!
//! Flow model: Kapila's mechanical-equilibrium model generalized to three
//! phases. Conserved u = (α₀ρ₀, α₁ρ₁, α₂ρ₂, ρu, ρv, ρw, E, α₀, α₁), with
//! α₂ = 1 − α₀ − α₁ and one NASG equation of state per phase:
//!   ∂(α_kρ_k)/∂t + ∇·(α_kρ_k u) = 0,
//!   ∂α_k/∂t + ∇·(α_k u) = ρc² z_k ∇·u,   z_k = α_k/(ρ_k c_k²),
//! where 1/(ρc²) = Σ z_k is Wood's mixture sound speed (Kapila's
//! K_k = ρc² z_k − α_k reduces to the two-phase formula). The phase
//! coefficients ρc² z_k sum to one, so Σα_k = 1 is kept by the dynamics.
//! The solid here has no shear strength ("hydrodynamic solid"); its
//! elasticity (GPR distortion per phase) is a separate rung.
//!
//! Thermodynamics (`Equilibrium3`): at given specific volume and energy
//! the stable state maximizes entropy. Candidates: each single phase, each
//! two-phase coexistence (equal p, T, Gibbs energy), and the triple point
//! (three phases at fixed T_t, p_t with mass fractions from v and e). One
//! rule covers melting, freezing, boiling, condensation, sublimation and
//! triple-point coexistence. Relaxation runs only where it can matter:
//! cells with two or more phases present, or a single phase that has
//! become metastable (another phase has lower Gibbs energy at its p, T).

use crate::elements::{Phase, ThreePhase};
use crate::eos::Nasg;
use crate::phase_change::{e_pt, pair_equilibrium, v_pt};
use crate::system::{Mat, System};

pub const N3: usize = 9;
type V = [f64; N3];
type M = Mat<N3>;
const IU: usize = 3;
const IP: usize = 6;
const IA: usize = 7;

#[derive(Clone, Copy, Debug)]
pub struct Kapila3 {
    pub eos: [Nasg; 3],
}

impl Kapila3 {
    fn phase_properties(&self, w: &V) -> [(f64, f64, f64); 3] {
        let al = Self::alphas(w);
        std::array::from_fn(|k| {
            let e = &self.eos[k];
            let rk = if al[k] > 1e-12 { w[k] / al[k] } else { 0.0 };
            let t = if rk > 0.0 { e.temperature(rk, w[IP]) } else { f64::NAN };
            (al[k], t, rk * e.cv)
        })
    }
    #[inline]
    pub fn alphas(w: &V) -> [f64; 3] {
        [w[IA], w[IA + 1], 1.0 - w[IA] - w[IA + 1]]
    }
    #[inline]
    pub fn rho(w: &V) -> f64 {
        w[0] + w[1] + w[2]
    }
    fn db(&self, w: &V) -> (f64, f64) {
        let al = Self::alphas(w);
        let (mut d, mut b) = (0.0, 0.0);
        for k in 0..3 {
            let e = &self.eos[k];
            let vol = al[k] - e.b * w[k];
            d += vol / (e.gamma - 1.0);
            b += e.gamma * e.p_inf * vol / (e.gamma - 1.0) + e.q * w[k];
        }
        (d, b)
    }
    fn rho_e(&self, w: &V) -> f64 {
        let (d, b) = self.db(w);
        w[IP] * d + b
    }
    /// z_k = α_k/(ρ_k c_k²) = (α_k − b_k α_kρ_k)/(γ_k (p + p∞k)).
    fn z(&self, w: &V) -> [f64; 3] {
        let al = Self::alphas(w);
        std::array::from_fn(|k| {
            let e = &self.eos[k];
            (al[k] - e.b * w[k]) / (e.gamma * (w[IP] + e.p_inf))
        })
    }
    /// Wood's ρc².
    pub fn rho_c2(&self, w: &V) -> f64 {
        let z = self.z(w);
        1.0 / (z[0] + z[1] + z[2])
    }
    fn kk(&self, w: &V) -> [f64; 3] {
        let rc2 = self.rho_c2(w);
        let z = self.z(w);
        let al = Self::alphas(w);
        std::array::from_fn(|k| rc2 * z[k] - al[k])
    }
    fn drhoe_dalpha(&self, p: f64, k: usize) -> f64 {
        let f = |e: &Nasg| (p + e.gamma * e.p_inf) / (e.gamma - 1.0);
        f(&self.eos[k]) - f(&self.eos[2])
    }
    fn p_inf_mix(&self, w: &V) -> f64 {
        let al = Self::alphas(w);
        (0..3).map(|k| al[k].clamp(0.0, 1.0) * self.eos[k].p_inf).sum()
    }
}

impl System<N3> for Kapila3 {
    fn to_prim(&self, u: &V) -> V {
        let rho = u[0] + u[1] + u[2];
        let mut w = *u;
        for c in 0..3 {
            w[IU + c] = u[IU + c] / rho;
        }
        let rhoe = u[IP] - 0.5 * rho * (w[3] * w[3] + w[4] * w[4] + w[5] * w[5]);
        let (d, b) = self.db(&w);
        w[IP] = (rhoe - b) / d;
        w
    }

    #[inline]
    fn velocity_conserved(&self, u: &V, d: usize) -> f64 {
        u[3 + d] / (u[0] + u[1] + u[2])
    }
    fn to_cons(&self, w: &V) -> V {
        let rho = Self::rho(w);
        let mut u = *w;
        for c in 0..3 {
            u[IU + c] = rho * w[IU + c];
        }
        u[IP] = self.rho_e(w) + 0.5 * rho * (w[3] * w[3] + w[4] * w[4] + w[5] * w[5]);
        u
    }
    fn sound_speed(&self, w: &V) -> f64 {
        (self.rho_c2(w).max(0.0) / Self::rho(w)).sqrt()
    }
    fn vel_index(&self) -> usize {
        IU
    }
    fn momentum_index(&self, d: usize) -> usize {
        IU + d
    }
    fn energy_index(&self) -> Option<usize> {
        Some(IP)
    }
    fn density(&self, w: &V) -> f64 {
        Self::rho(w)
    }
    fn density_indices(&self) -> std::ops::Range<usize> {
        0..3
    }
    fn phases(&self, w: &V) -> Vec<(f64, f64, f64)> {
        self.phase_properties(w).to_vec()
    }
    fn conduction_properties(&self, w: &V, conductivity: &[f64]) -> (f64, f64, f64) {
        let phases = self.phase_properties(w);
        crate::thermal::properties(self.temperature_from_phases(w, &phases), &phases, conductivity)
    }
    fn temperature(&self, w: &V) -> f64 {
        crate::thermal::mixture_temperature(&self.phase_properties(w), &[w[0], w[1], w[2]])
    }
    fn temperature_from_phases(&self, w: &V, phases: &[(f64, f64, f64)]) -> f64 {
        crate::thermal::mixture_temperature(phases, &[w[0], w[1], w[2]])
    }
    fn phase_stiffness(&self, w: &V) -> Vec<(f64, f64)> {
        let al = Self::alphas(w);
        (0..3)
            .map(|k| {
                let e = &self.eos[k];
                let rk = if al[k] > 1e-12 { w[k] / al[k] } else { 0.0 };
                (al[k], e.gamma * (w[IP] + e.p_inf) / (1.0 - e.b * rk).max(1e-12))
            })
            .collect()
    }
    fn eigen(&self, w: &V, d: usize, l: &mut M, r: &mut M) {
        let rho = Self::rho(w);
        let c = self.sound_speed(w);
        let c2 = c * c;
        let kk = self.kk(w);
        let un = IU + d;
        let (t1, t2) = (IU + (d + 1) % 3, IU + (d + 2) % 3);
        for row in l.iter_mut().chain(r.iter_mut()) {
            *row = [0.0; N3];
        }
        // fields: 0..3 densities, 3 u−c, 4 t1, 5 t2, 6 u+c, 7, 8 volume fractions
        for k in 0..3 {
            l[k][k] = 1.0;
            l[k][IP] = -w[k] / rho / c2;
            r[k][k] = 1.0;
        }
        l[3][un] = -rho / (2.0 * c);
        l[3][IP] = 1.0 / (2.0 * c2);
        l[4][t1] = 1.0;
        l[5][t2] = 1.0;
        l[6][un] = rho / (2.0 * c);
        l[6][IP] = 1.0 / (2.0 * c2);
        for j in 0..2 {
            l[7 + j][IA + j] = 1.0;
            l[7 + j][IP] = kk[j] / rho / c2;
            r[IA + j][7 + j] = 1.0;
        }
        r[t1][4] = 1.0;
        r[t2][5] = 1.0;
        for (col, sg) in [(3usize, -1.0), (6usize, 1.0)] {
            for k in 0..3 {
                r[k][col] = w[k] / rho;
            }
            r[un][col] = sg * c / rho;
            r[IP][col] = c2;
            r[IA][col] = -kk[0] / rho;
            r[IA + 1][col] = -kk[1] / rho;
        }
    }
    fn eigenvalues(&self, w: &V, d: usize) -> V {
        let (u, c) = (w[IU + d], self.sound_speed(w));
        [u, u, u, u - c, u, u, u + c, u, u]
    }
    fn riemann_flux(&self, wl: &V, wr: &V, d: usize) -> (V, f64) {
        let un = IU + d;
        let (cl, cr) = (self.sound_speed(wl), self.sound_speed(wr));
        let (ul, ur) = (wl[un], wr[un]);
        let sl = (ul - cl).min(ur - cr);
        let sr = (ul + cl).max(ur + cr);
        if sl >= 0.0 {
            return (self.flux(wl, d), ul);
        }
        if sr <= 0.0 {
            return (self.flux(wr, d), ur);
        }
        let (rl, rr) = (Self::rho(wl), Self::rho(wr));
        let ss = (wr[IP] - wl[IP] + rl * ul * (sl - ul) - rr * ur * (sr - ur)) / (rl * (sl - ul) - rr * (sr - ur));
        let (w, s) = if ss >= 0.0 { (wl, sl) } else { (wr, sr) };
        let q = self.to_cons(w);
        let fx = self.flux(w, d);
        let (rho, u) = (Self::rho(w), w[un]);
        let fac = (s - u) / (s - ss);
        let mut qs = [0.0; N3];
        for k in 0..3 {
            qs[k] = w[k] * fac;
        }
        for c in 0..3 {
            qs[IU + c] = rho * fac * w[IU + c];
        }
        qs[un] = rho * fac * ss;
        qs[IP] = rho * fac * (q[IP] / rho + (ss - u) * (ss + w[IP] / (rho * (s - u))));
        qs[IA] = w[IA] * fac;
        qs[IA + 1] = w[IA + 1] * fac;
        let mut f = [0.0; N3];
        for k in 0..N3 {
            f[k] = fx[k] + s * (qs[k] - q[k]);
        }
        (f, ss)
    }
    fn split_flux(&self, a: &V, b: &V, d: usize) -> V {
        let un = IU + d;
        let u = 0.5 * (a[un] + b[un]);
        let rho = 0.5 * (Self::rho(a) + Self::rho(b));
        let c = rho * u;
        let mut f = [0.0; N3];
        for k in 0..3 {
            f[k] = 0.5 * (a[k] + b[k]) * u;
        }
        for k in IU..IU + 3 {
            f[k] = c * 0.5 * (a[k] + b[k]);
        }
        f[un] += 0.5 * (a[IP] + b[IP]);
        let kin = c * 0.5 * (a[3] * b[3] + a[4] * b[4] + a[5] * b[5]);
        let rhoe = 0.5 * (self.rho_e(a) + self.rho_e(b));
        let pw = 0.5 * (a[un] * b[IP] + b[un] * a[IP]);
        f[IP] = kin + u * rhoe + pw;
        f[IA] = 0.5 * (a[IA] + b[IA]) * u;
        f[IA + 1] = 0.5 * (a[IA + 1] + b[IA + 1]) * u;
        f
    }
    fn flux(&self, w: &V, d: usize) -> V {
        let un = IU + d;
        let q = self.to_cons(w);
        let u = w[un];
        let mut f = [0.0; N3];
        for k in 0..6 {
            f[k] = q[k] * u;
        }
        f[un] += w[IP];
        f[IP] = (q[IP] + w[IP]) * u;
        f[IA] = w[IA] * u;
        f[IA + 1] = w[IA + 1] * u;
        f
    }
    fn dudw(&self, w: &V) -> M {
        let rho = Self::rho(w);
        let p = w[IP];
        let (d, _) = self.db(w);
        let ke = 0.5 * (w[3] * w[3] + w[4] * w[4] + w[5] * w[5]);
        let mut j = [[0.0; N3]; N3];
        for k in 0..3 {
            j[k][k] = 1.0;
            let e = &self.eos[k];
            j[IP][k] = -e.b * (p + e.gamma * e.p_inf) / (e.gamma - 1.0) + e.q + ke;
        }
        for c in 0..3 {
            let k = IU + c;
            for m in 0..3 {
                j[k][m] = w[k];
            }
            j[k][k] = rho;
            j[IP][k] = rho * w[k];
        }
        j[IP][IP] = d;
        j[IP][IA] = self.drhoe_dalpha(p, 0);
        j[IP][IA + 1] = self.drhoe_dalpha(p, 1);
        j[IA][IA] = 1.0;
        j[IA + 1][IA + 1] = 1.0;
        j
    }
    fn admissible(&self, w: &V) -> bool {
        w.iter().all(|v| v.is_finite()) && Self::rho(w) > 0.0 && Self::alphas(w).iter().all(|a| (-0.05..=1.05).contains(a)) && self.rho_c2(w) > 0.0
    }
    fn nc_vars(&self) -> std::ops::Range<usize> {
        IA..IA + 2
    }
    fn nc_coef_k(&self, w: &V, q: usize) -> f64 {
        self.rho_c2(w) * self.z(w)[q - IA]
    }
    fn dilatation_rate(&self, w: &V, divu: f64) -> f64 {
        let al = Self::alphas(w);
        let mut pm = f64::INFINITY;
        for k in 0..3 {
            if al[k] > 1e-6 {
                pm = pm.min(self.eos[k].p_inf);
            }
        }
        self.rho_c2(w).max(0.0) * divu.abs() / (w[IP] + pm).max(1e-300)
    }
    fn sensor(&self, w: &V, out: &mut [f64; 8]) -> usize {
        out[0] = Self::rho(w).max(1e-300).ln();
        out[1] = (w[IP] + self.p_inf_mix(w)).max(1e-300).ln();
        out[2] = w[IA];
        out[3] = w[IA + 1];
        4
    }
}

/// NASG specific entropy s(p, T) (log form).
pub fn entropy(e: &Nasg, p: f64, t: f64) -> f64 {
    e.cv * (e.gamma * t.ln() - (e.gamma - 1.0) * (p + e.p_inf).ln()) + e.q_prime
}

/// Equilibrium state of the three-phase substance.
#[derive(Clone, Copy, Debug)]
pub struct State3 {
    pub p: f64,
    pub t: f64,
    /// mass fractions (solid, liquid, vapour)
    pub y: [f64; 3],
    pub entropy: f64,
}

/// Maximum-entropy phase equilibrium of one substance (solid, liquid,
/// vapour NASG phases, e.g. from `elements::fit`).
#[derive(Clone, Debug)]
pub struct Equilibrium3 {
    pub model: ThreePhase,
    pub triple: (f64, f64),
    /// Gibbs-energy margin (J/kg) by which a single phase must exceed
    /// another before the new phase nucleates (0: no metastability)
    pub nucleation_margin: f64,
    /// a phase is "present" above this volume fraction
    pub present: f64,
}

impl Equilibrium3 {
    pub fn new(model: ThreePhase, t_guess: f64) -> Option<Self> {
        let triple = model.triple_point(t_guess)?;
        Some(Equilibrium3 { model, triple, nucleation_margin: 0.0, present: 1e-6 })
    }

    fn eos(&self) -> [Nasg; 3] {
        [self.model.solid, self.model.liquid, self.model.vapor]
    }

    /// Single-phase state of phase k at (v, e), if physical.
    fn single(&self, k: usize, v: f64, e: f64) -> Option<State3> {
        let n = self.eos()[k];
        if v <= n.b {
            return None;
        }
        let p = (n.gamma - 1.0) * (e - n.q) / (v - n.b) - n.gamma * n.p_inf;
        if !(p + n.p_inf > 0.0) {
            return None;
        }
        let t = (p + n.p_inf) * (v - n.b) / ((n.gamma - 1.0) * n.cv);
        if !(t > 0.0) {
            return None;
        }
        let mut y = [0.0; 3];
        y[k] = 1.0;
        Some(State3 { p, t, y, entropy: entropy(&n, p, t) })
    }

    /// Two-phase coexistence of phases a < b at (v, e).
    fn pair(&self, a: usize, b: usize, v: f64, e: f64) -> Option<State3> {
        let eos = self.eos();
        let eq = pair_equilibrium(&eos[a], &eos[b], v, e, self.triple.0)?;
        let mut y = [0.0; 3];
        y[a] = 1.0 - eq.y;
        y[b] = eq.y;
        let s = (1.0 - eq.y) * entropy(&eos[a], eq.p, eq.t) + eq.y * entropy(&eos[b], eq.p, eq.t);
        Some(State3 { p: eq.p, t: eq.t, y, entropy: s })
    }

    /// Three-phase coexistence at the triple point: mass fractions from
    /// Σy = 1, Σ y_k v_k = v, Σ y_k e_k = e.
    fn triple_state(&self, v: f64, e: f64) -> Option<State3> {
        let (tt, pt) = self.triple;
        let eos = self.eos();
        let vk: [f64; 3] = std::array::from_fn(|k| v_pt(&eos[k], pt, tt));
        let ek: [f64; 3] = std::array::from_fn(|k| e_pt(&eos[k], pt, tt));
        // solve [1 1 1; v0 v1 v2; e0 e1 e2] y = [1, v, e]
        let a = [[1.0, 1.0, 1.0], [vk[0], vk[1], vk[2]], [ek[0], ek[1], ek[2]]];
        let rhs = [1.0, v, e];
        let det = |m: &[[f64; 3]; 3]| m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
        let d0 = det(&a);
        if d0.abs() < 1e-300 {
            return None;
        }
        let y: [f64; 3] = std::array::from_fn(|c| {
            let mut m = a;
            for r in 0..3 {
                m[r][c] = rhs[r];
            }
            det(&m) / d0
        });
        if y.iter().any(|&x| x < 0.0 || x > 1.0) {
            return None;
        }
        let s = (0..3).map(|k| y[k] * entropy(&eos[k], pt, tt)).sum();
        Some(State3 { p: pt, t: tt, y, entropy: s })
    }

    /// Stable state at (v, e): the admissible candidate of largest entropy.
    pub fn equilibrate(&self, v: f64, e: f64) -> Option<State3> {
        let mut best: Option<State3> = None;
        let mut consider = |c: Option<State3>| {
            if let Some(c) = c {
                if c.entropy.is_finite() && best.map_or(true, |b| c.entropy > b.entropy) {
                    best = Some(c);
                }
            }
        };
        for k in 0..3 {
            consider(self.single(k, v, e));
        }
        for (a, b) in [(0, 1), (1, 2), (0, 2)] {
            consider(self.pair(a, b, v, e));
        }
        consider(self.triple_state(v, e));
        best
    }

    /// Relax a cell of the `Kapila3` flow to equilibrium where it can
    /// change: two or more phases present, or one phase metastable.
    pub fn relax(&self, sys: &Kapila3, u: &mut V) -> bool {
        let w = sys.to_prim(u);
        let al = Kapila3::alphas(&w);
        let present: Vec<usize> = (0..3).filter(|&k| al[k] > self.present).collect();
        if present.len() == 1 {
            let k = present[0];
            let n = sys.eos[k];
            let rho_k = w[k] / al[k];
            let (p, t) = (w[IP], n.temperature(rho_k, w[IP]));
            let gk = n.gibbs(p, t);
            let eos = self.eos();
            let metastable = (0..3).any(|j| j != k && eos[j].gibbs(p, t) + self.nucleation_margin < gk);
            if !metastable {
                return false;
            }
        }
        let rho = Kapila3::rho(&w);
        let ke = 0.5 * rho * (w[3] * w[3] + w[4] * w[4] + w[5] * w[5]);
        let (v, e) = (1.0 / rho, (u[IP] - ke) / rho);
        let Some(st) = self.equilibrate(v, e) else { return false };
        let eos = self.eos();
        // absent phases are exactly absent (zero volume and mass): a trace of
        // a soft phase (vapour) inside a condensed one would take a share of
        // every compression far beyond its volume (Kapila: ∝ α/(ρc²)) and go
        // negative under the large steps of the all-speed mode
        let trace = 0.0;
        let mut alpha = [0.0; 3];
        let mut part = [0.0; 3];
        for k in 0..3 {
            let rk = 1.0 / v_pt(&eos[k], st.p, st.t);
            let y = st.y[k].max(0.0);
            part[k] = y * rho;
            alpha[k] = part[k] / rk;
        }
        for k in 0..3 {
            if alpha[k] < trace {
                alpha[k] = trace;
            }
        }
        let s: f64 = alpha.iter().sum();
        for a in alpha.iter_mut() {
            *a /= s;
        }
        let mut wn = w;
        wn[0] = part[0];
        wn[1] = part[1];
        wn[2] = part[2];
        wn[IA] = alpha[0];
        wn[IA + 1] = alpha[1];
        wn[IP] = st.p;
        let un = sys.to_cons(&wn);
        for k in 0..3 {
            u[k] = un[k];
        }
        u[IA] = un[IA];
        u[IA + 1] = un[IA + 1];
        // total energy is kept exactly: the closure with the new fractions
        // reproduces p up to the trace-phase renormalization
        true
    }

    /// Phase label of a state (for reports).
    pub fn label(st: &State3) -> &'static str {
        let n = st.y.iter().filter(|&&y| y > 1e-9).count();
        match n {
            3 => "solid+liquid+vapour",
            2 => {
                if st.y[2] < 1e-9 {
                    "solid+liquid"
                } else if st.y[0] < 1e-9 {
                    "liquid+vapour"
                } else {
                    "solid+vapour"
                }
            }
            _ => {
                if st.y[0] > 0.5 {
                    "solid"
                } else if st.y[1] > 0.5 {
                    "liquid"
                } else {
                    "vapour"
                }
            }
        }
    }
}

#[allow(dead_code)]
fn _phase(_p: Phase) {}
