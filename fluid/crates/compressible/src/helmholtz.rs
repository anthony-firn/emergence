//! Reference equations of state: multiparameter Helmholtz free energy
//! formulations (IAPWS-95 for water, Span et al. 2000 for nitrogen),
//!   a(ρ, T)/(RT) = α(δ, τ) = α⁰(δ, τ) + αʳ(δ, τ),  δ = ρ/ρ_c, τ = T_c/T.
//! Every property follows from α and its first and second derivatives,
//! which are computed exactly by forward-mode automatic differentiation
//! (a second-order dual number in (δ, τ)), so each term type is written
//! once and its derivatives cannot drift from it. Coefficients:
//! `reference_data.rs`, generated from CoolProp's fluid files.
//!
//! These formulations are the experimental-accuracy rung of the EOS ladder
//! (water: density to ~0.0001–0.003%, sound speed ~0.03–0.3% in the
//! liquid), valid in the stable and metastable regions and through the
//! critical point, where the NASG fits are not.

/// Value and first/second partial derivatives of f(δ, τ).
#[derive(Clone, Copy, Debug, Default)]
pub struct D2 {
    pub v: f64,
    pub x: f64,
    pub y: f64,
    pub xx: f64,
    pub xy: f64,
    pub yy: f64,
}

impl D2 {
    pub const fn c(v: f64) -> Self {
        D2 { v, x: 0.0, y: 0.0, xx: 0.0, xy: 0.0, yy: 0.0 }
    }
    /// unary chain rule with f, f', f'' at self.v
    #[inline]
    fn chain(self, f: f64, f1: f64, f2: f64) -> Self {
        D2 {
            v: f,
            x: f1 * self.x,
            y: f1 * self.y,
            xx: f2 * self.x * self.x + f1 * self.xx,
            xy: f2 * self.x * self.y + f1 * self.xy,
            yy: f2 * self.y * self.y + f1 * self.yy,
        }
    }
    pub fn exp(self) -> Self {
        let e = self.v.exp();
        self.chain(e, e, e)
    }
    pub fn ln(self) -> Self {
        self.chain(self.v.ln(), 1.0 / self.v, -1.0 / (self.v * self.v))
    }
    pub fn powf(self, a: f64) -> Self {
        if a == 0.0 {
            return D2::c(1.0);
        }
        let v = self.v;
        self.chain(v.powf(a), a * v.powf(a - 1.0), a * (a - 1.0) * v.powf(a - 2.0))
    }
    pub fn recip(self) -> Self {
        let v = self.v;
        self.chain(1.0 / v, -1.0 / (v * v), 2.0 / (v * v * v))
    }
}

impl std::ops::Add for D2 {
    type Output = D2;
    fn add(self, o: D2) -> D2 {
        D2 { v: self.v + o.v, x: self.x + o.x, y: self.y + o.y, xx: self.xx + o.xx, xy: self.xy + o.xy, yy: self.yy + o.yy }
    }
}
impl std::ops::Sub for D2 {
    type Output = D2;
    fn sub(self, o: D2) -> D2 {
        D2 { v: self.v - o.v, x: self.x - o.x, y: self.y - o.y, xx: self.xx - o.xx, xy: self.xy - o.xy, yy: self.yy - o.yy }
    }
}
impl std::ops::Mul for D2 {
    type Output = D2;
    fn mul(self, o: D2) -> D2 {
        D2 {
            v: self.v * o.v,
            x: self.x * o.v + self.v * o.x,
            y: self.y * o.v + self.v * o.y,
            xx: self.xx * o.v + 2.0 * self.x * o.x + self.v * o.xx,
            xy: self.xy * o.v + self.x * o.y + self.y * o.x + self.v * o.xy,
            yy: self.yy * o.v + 2.0 * self.y * o.y + self.v * o.yy,
        }
    }
}
impl std::ops::Mul<f64> for D2 {
    type Output = D2;
    fn mul(self, k: f64) -> D2 {
        D2 { v: self.v * k, x: self.x * k, y: self.y * k, xx: self.xx * k, xy: self.xy * k, yy: self.yy * k }
    }
}
impl std::ops::Add<f64> for D2 {
    type Output = D2;
    fn add(self, k: f64) -> D2 {
        D2 { v: self.v + k, ..self }
    }
}

#[derive(Clone, Debug)]
pub enum Term {
    /// n δ^d τ^t [exp(−δ^l) if l > 0]
    Power { n: f64, d: f64, t: f64, l: f64 },
    /// n δ^d τ^t exp(−η(δ−ε)² − β(τ−γ)²)
    Gauss { n: f64, d: f64, t: f64, eta: f64, eps: f64, beta: f64, gamma: f64 },
    /// IAPWS-95 near-critical term n Δ^b δ ψ, ψ = exp(−C(δ−1)² − D(τ−1)²),
    /// θ = (1−τ) + A((δ−1)²)^{1/(2β)}, Δ = θ² + B((δ−1)²)^a
    NonAnalytic { n: f64, a: f64, b: f64, big_b: f64, c: f64, big_d: f64, big_a: f64, beta: f64 },
    /// ln δ + a₁ + a₂τ
    Lead { a1: f64, a2: f64 },
    /// a ln τ
    LogTau { a: f64 },
    /// n ln(1 − exp(−tτ))
    PlanckEinstein { n: f64, t: f64 },
    /// n τ^t
    PowerTau { n: f64, t: f64 },
}

impl Term {
    fn eval(&self, del: D2, tau: D2) -> D2 {
        match *self {
            Term::Power { n, d, t, l } => {
                let mut r = del.powf(d) * tau.powf(t) * n;
                if l > 0.0 {
                    r = r * (del.powf(l) * -1.0).exp();
                }
                r
            }
            Term::Gauss { n, d, t, eta, eps, beta, gamma } => {
                let a = del + (-eps);
                let b = tau + (-gamma);
                del.powf(d) * tau.powf(t) * ((a * a) * -eta - (b * b) * beta).exp() * n
            }
            Term::NonAnalytic { n, a, b, big_b, c, big_d, big_a, beta } => {
                let dm = del + (-1.0);
                let dm2 = dm * dm;
                let tm = tau + (-1.0);
                let psi = (dm2 * -c - (tm * tm) * big_d).exp();
                let theta = (tau * -1.0 + 1.0) + dm2.powf(1.0 / (2.0 * beta)) * big_a;
                let delta_ = theta * theta + dm2.powf(a) * big_b;
                delta_.powf(b) * del * psi * n
            }
            Term::Lead { a1, a2 } => del.ln() + tau * a2 + a1,
            Term::LogTau { a } => tau.ln() * a,
            Term::PlanckEinstein { n, t } => ((tau * -t).exp() * -1.0 + 1.0).ln() * n,
            Term::PowerTau { n, t } => tau.powf(t) * n,
        }
    }
}

/// Ancillary correlation (initial guesses only).
#[derive(Clone, Debug)]
pub struct Ancillary {
    pub kind: &'static str,
    pub t_r: f64,
    pub reducing: f64,
    pub tau_r: bool,
    pub n: &'static [f64],
    pub t: &'static [f64],
}

impl Ancillary {
    pub fn eval(&self, temp: f64) -> f64 {
        let th = 1.0 - temp / self.t_r;
        let s: f64 = self.n.iter().zip(self.t).map(|(n, t)| n * th.powf(*t)).sum();
        match self.kind {
            "rhoLnoexp" => self.reducing * (1.0 + s),
            _ => {
                let f = if self.tau_r { self.t_r / temp } else { 1.0 };
                self.reducing * (f * s).exp()
            }
        }
    }
}

#[derive(Clone, Debug)]
pub struct Helmholtz {
    pub name: &'static str,
    pub molar_mass: f64,
    pub gas_constant: f64,
    pub t_c: f64,
    pub rho_c: f64,
    pub p_c: f64,
    pub t_triple: f64,
    pub residual: Vec<Term>,
    pub ideal: Vec<Term>,
    pub anc_p: Ancillary,
    pub anc_rho_l: Ancillary,
    pub anc_rho_v: Ancillary,
}

/// α⁰ and αʳ with derivatives at a state.
#[derive(Clone, Copy, Debug)]
pub struct Alphas {
    pub delta: f64,
    pub tau: f64,
    pub a0: D2,
    pub ar: D2,
}

impl Helmholtz {
    /// Specific gas constant, J/(kg K).
    pub fn r(&self) -> f64 {
        self.gas_constant / self.molar_mass
    }

    pub fn alphas(&self, rho: f64, t: f64) -> Alphas {
        let mut delta = rho / self.rho_c;
        // the non-analytic terms have an integrable singularity exactly at δ = 1
        if (delta - 1.0).abs() < 1e-12 {
            delta = 1.0 + 1e-12;
        }
        let tau = self.t_c / t;
        let del = D2 { v: delta, x: 1.0, ..Default::default() };
        let ta = D2 { v: tau, y: 1.0, ..Default::default() };
        let sum = |ts: &[Term]| ts.iter().fold(D2::c(0.0), |acc, term| acc + term.eval(del, ta));
        Alphas { delta, tau, a0: sum(&self.ideal), ar: sum(&self.residual) }
    }

    pub fn pressure(&self, rho: f64, t: f64) -> f64 {
        let a = self.alphas(rho, t);
        rho * self.r() * t * (1.0 + a.delta * a.ar.x)
    }

    /// All the thermodynamic properties at (ρ, T).
    pub fn props(&self, rho: f64, t: f64) -> Props {
        let a = self.alphas(rho, t);
        let (d, ta, r) = (a.delta, a.tau, self.r());
        let (ar, a0) = (a.ar, a.a0);
        let p = rho * r * t * (1.0 + d * ar.x);
        let u = r * t * ta * (a0.y + ar.y);
        let s = r * (ta * (a0.y + ar.y) - a0.v - ar.v);
        let h = r * t * (1.0 + ta * (a0.y + ar.y) + d * ar.x);
        let g = r * t * (1.0 + a0.v + ar.v + d * ar.x);
        let cv = -r * ta * ta * (a0.yy + ar.yy);
        let num = 1.0 + d * ar.x - d * ta * ar.xy;
        let dpdrho = r * t * (1.0 + 2.0 * d * ar.x + d * d * ar.xx);
        let cp = cv + r * num * num / (1.0 + 2.0 * d * ar.x + d * d * ar.xx);
        let w2 = r * t * (1.0 + 2.0 * d * ar.x + d * d * ar.xx - num * num / (ta * ta * (a0.yy + ar.yy)));
        Props { rho, t, p, u, s, h, g, cv, cp, w: w2.max(0.0).sqrt(), dpdrho_t: dpdrho }
    }

    /// Temperature at density ρ and specific internal energy e (Newton on
    /// u(ρ, T) = e; ∂u/∂T|ρ = c_v), starting from `t0`.
    pub fn temperature(&self, rho: f64, e: f64, t0: f64) -> Option<f64> {
        let mut t = t0;
        for _ in 0..100 {
            let pr = self.props(rho, t);
            let dt = (pr.u - e) / pr.cv;
            let tn = (t - dt).clamp(0.5 * t, 2.0 * t);
            if !tn.is_finite() {
                return None;
            }
            if (tn - t).abs() < 1e-12 * t {
                return Some(tn);
            }
            t = tn;
        }
        Some(t)
    }

    /// Temperature at (ρ, e) choosing the mechanically stable root. Inside
    /// the dome u(ρ, T) is not monotone (c_v < 0 in the unstable region)
    /// and Newton can converge to an unphysical root; at fixed ρ at most one
    /// root is stable or metastable. Scan T on a log grid for sign changes
    /// of u − e, bisect each, keep the root with ∂p/∂ρ|_T > 0 and c_v > 0.
    pub fn temperature_stable(&self, rho: f64, e: f64) -> Option<f64> {
        let (t_lo, t_hi) = (0.5 * self.t_triple, 6.0 * self.t_c);
        let n = 160;
        let tt = |i: usize| t_lo * (t_hi / t_lo).powf(i as f64 / n as f64);
        let f = |t: f64| self.props(rho, t).u - e;
        let mut prev = (tt(0), f(tt(0)));
        let mut best: Option<f64> = None;
        for i in 1..=n {
            let t = tt(i);
            let ft = f(t);
            if prev.1.is_finite() && ft.is_finite() && prev.1 * ft <= 0.0 {
                let (mut lo, mut hi, flo) = (prev.0, t, prev.1);
                for _ in 0..80 {
                    let mid = 0.5 * (lo + hi);
                    if (f(mid) < 0.0) == (flo < 0.0) {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                let root = 0.5 * (lo + hi);
                let pr = self.props(rho, root);
                if pr.dpdrho_t > 0.0 && pr.cv > 0.0 {
                    best = Some(root);
                }
            }
            prev = (t, ft);
        }
        best
    }

    /// Density at pressure p and temperature T (Newton on p(ρ, T) = p),
    /// starting from `rho0` (chooses the branch).
    pub fn density(&self, p: f64, t: f64, rho0: f64) -> Option<f64> {
        let mut rho = rho0;
        for _ in 0..200 {
            let pr = self.props(rho, t);
            if !(pr.dpdrho_t > 0.0) {
                return None;
            }
            let rn = (rho - (pr.p - p) / pr.dpdrho_t).clamp(0.5 * rho, 2.0 * rho);
            if (rn - rho).abs() < 1e-13 * rho {
                return Some(rn);
            }
            rho = rn;
        }
        Some(rho)
    }

    /// Saturation at temperature T: (p, ρ_liquid, ρ_vapour) from the Maxwell
    /// conditions p_l = p_v, g_l = g_v (Newton on (ρ_l, ρ_v), ancillary
    /// starting values; ∂g/∂ρ|_T = (1/ρ) ∂p/∂ρ|_T).
    pub fn saturation(&self, t: f64) -> Option<(f64, f64, f64)> {
        if !(t < self.t_c && t >= 0.95 * self.t_triple) {
            return None;
        }
        // the density ancillaries are molar (mol/m³)
        let (mut rl, mut rv) = (self.anc_rho_l.eval(t) * self.molar_mass, self.anc_rho_v.eval(t) * self.molar_mass);
        for _ in 0..200 {
            let (l, v) = (self.props(rl, t), self.props(rv, t));
            let (f1, f2) = (l.p - v.p, l.g - v.g);
            // Jacobian w.r.t. (ρ_l, ρ_v)
            let (a11, a12) = (l.dpdrho_t, -v.dpdrho_t);
            let (a21, a22) = (l.dpdrho_t / rl, -v.dpdrho_t / rv);
            let det = a11 * a22 - a12 * a21;
            if det.abs() < 1e-300 {
                return None;
            }
            let dl = (f1 * a22 - f2 * a12) / det;
            let dv = (a11 * f2 - a21 * f1) / det;
            let (nl, nv) = ((rl - dl).clamp(0.8 * rl, 1.2 * rl), (rv - dv).clamp(0.5 * rv, 1.5 * rv));
            let conv = (nl - rl).abs() < 1e-13 * rl && (nv - rv).abs() < 1e-13 * rv;
            rl = nl;
            rv = nv;
            if conv {
                break;
            }
        }
        // the vapour side: on the liquid side p is a small difference of
        // large numbers at low temperature (1 + δα_δ ~ 1e-6 at 275 K)
        let p = self.props(rv, t).p;
        (p.is_finite() && rl > rv).then_some((p, rl, rv))
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Props {
    pub rho: f64,
    pub t: f64,
    pub p: f64,
    pub u: f64,
    pub s: f64,
    pub h: f64,
    pub g: f64,
    pub cv: f64,
    pub cp: f64,
    pub w: f64,
    pub dpdrho_t: f64,
}

#[cfg(test)]
mod tests {
    use crate::reference_data::water_iapws95;

    /// IAPWS-95 release, Table 7 (single-phase verification values).
    #[test]
    fn iapws95_table7() {
        let w = water_iapws95();
        // (T, ρ, p [MPa], c_v [kJ/kg K], w [m/s], s [kJ/kg K])
        let rows = [
            (300.0, 0.9965560e3, 0.992418352e-1, 0.413018112e1, 0.150151914e4, 0.393062643),
            (300.0, 0.1005308e4, 0.200022515e2, 0.406798347e1, 0.153492501e4, 0.387405401),
            (300.0, 0.1188202e4, 0.700004704e3, 0.346135580e1, 0.244357992e4, 0.132609616),
            (500.0, 0.4350000, 0.999679423e-1, 0.150817541e1, 0.548314253e3, 0.794488271e1),
            (500.0, 0.4532000e1, 0.999938125, 0.166991025e1, 0.535739001e3, 0.682502725e1),
            (500.0, 0.8380250e3, 0.100003858e2, 0.322106219e1, 0.127128441e4, 0.256690919e1),
            (647.0, 0.3580000e3, 0.220384756e2, 0.618315728e1, 0.252145078e3, 0.432092307e1),
            (900.0, 0.2410000, 0.100062559, 0.175890657e1, 0.724027147e3, 0.916653194e1),
            (900.0, 0.5261500e2, 0.200000690e2, 0.193510526e1, 0.698445674e3, 0.659070225e1),
            (900.0, 0.8707690e3, 0.700000006e3, 0.266422350e1, 0.201933608e4, 0.417223802e1),
        ];
        for (t, rho, p, cv, ws, s) in rows {
            let pr = w.props(rho, t);
            let rel = |a: f64, b: f64| (a / b - 1.0).abs();
            assert!(rel(pr.p * 1e-6, p) < 1e-8, "p at {t} {rho}: {} vs {p}", pr.p * 1e-6);
            assert!(rel(pr.cv * 1e-3, cv) < 1e-8, "cv at {t} {rho}: {} vs {cv}", pr.cv * 1e-3);
            assert!(rel(pr.w, ws) < 1e-8, "w at {t} {rho}: {} vs {ws}", pr.w);
            assert!(rel(pr.s * 1e-3, s) < 1e-7, "s at {t} {rho}: {} vs {s}", pr.s * 1e-3);
        }
    }

    /// IAPWS-95 release, Table 8 (saturation verification values).
    #[test]
    fn iapws95_table8_saturation() {
        let w = water_iapws95();
        for (t, p, rl, rv) in [(275.0, 0.698451167e-3, 0.999887406e3, 0.550664919e-2), (450.0, 0.932203564, 0.890341250e3, 0.481200360e1), (625.0, 0.169082693e2, 0.567090385e3, 0.118290280e3)] {
            let (ps, l, v) = w.saturation(t).unwrap();
            assert!((ps * 1e-6 / p - 1.0).abs() < 1e-8, "p_sat {t}: {} vs {p}", ps * 1e-6);
            assert!((l / rl - 1.0).abs() < 1e-8, "rho_l {t}: {l} vs {rl}");
            assert!((v / rv - 1.0).abs() < 1e-8, "rho_v {t}: {v} vs {rv}");
        }
    }

    /// Against CoolProp 8 (single-phase states and the saturation curve)
    /// for both fluids: data/reference_eos_check.csv.
    #[test]
    fn matches_coolprop() {
        use crate::reference_data::nitrogen_span2000;
        let (w, n) = (water_iapws95(), nitrogen_span2000());
        let csv = include_str!("../data/reference_eos_check.csv");
        let mut worst = [0.0f64; 6];
        let mut worst_sat = [0.0f64; 3];
        let mut count = (0, 0);
        for line in csv.lines().skip(1) {
            let f: Vec<&str> = line.split(',').collect();
            let fl = if f[1] == "Water" { &w } else { &n };
            let num = |i: usize| f[i].parse::<f64>().unwrap();
            let (t, rho) = (num(2), num(3));
            if f[0] == "state" {
                let pr = fl.props(rho, t);
                let ours = [pr.p, pr.u, pr.s, pr.cv, pr.cp, pr.w];
                for k in 0..6 {
                    let r = num(4 + k);
                    // u and s: compare with an absolute scale (they cross zero)
                    // p: liquid pressures are small differences of large numbers;
                    // the natural scale is the bulk modulus ρw²
                    let scale = match k {
                        0 => r.abs().max(1e-3 * rho * num(9) * num(9)),
                        1 => r.abs().max(1e5),
                        2 => r.abs().max(1e3),
                        _ => r.abs(),
                    };
                    worst[k] = worst[k].max((ours[k] - r).abs() / scale);
                }
                count.0 += 1;
            } else {
                let (p, rl, rv) = fl.saturation(t).unwrap_or_else(|| panic!("saturation {} {t}", f[1]));
                for (k, (a, b)) in [(p, num(4)), (rl, num(5)), (rv, num(6))].iter().enumerate() {
                    worst_sat[k] = worst_sat[k].max((a / b - 1.0).abs());
                }
                count.1 += 1;
            }
        }
        println!("{} states: worst relative p {:.1e} u {:.1e} s {:.1e} cv {:.1e} cp {:.1e} w {:.1e}", count.0, worst[0], worst[1], worst[2], worst[3], worst[4], worst[5]);
        println!("{} saturation points: worst relative p {:.1e} rho_l {:.1e} rho_v {:.1e}", count.1, worst_sat[0], worst_sat[1], worst_sat[2]);
        // liquid-nitrogen pressures differ by up to ~6e-6 relative: with
        // ρc² ~ 1e9 Pa that is a 2e-8 relative density offset (unit
        // conversion), far below the EOS uncertainty (~0.02% in density)
        assert!(worst[0] < 1e-5 && worst[1..].iter().all(|&e| e < 1e-6), "{worst:?}");
        assert!(worst_sat.iter().all(|&e| e < 1e-6), "{worst_sat:?}");
    }
}
