//! Two-phase compressible flow: Kapila's five-equation model (the
//! mechanical-equilibrium limit of the Baer–Nunziato seven-equation model
//! under stiff pressure and velocity relaxation):
//!   ∂(α₁ρ₁)/∂t + ∇·(α₁ρ₁u) = 0,   ∂(α₂ρ₂)/∂t + ∇·(α₂ρ₂u) = 0,
//!   ∂(ρu)/∂t + ∇·(ρuu + pI) = 0,  ∂E/∂t + ∇·((E + p)u) = 0,
//!   ∂α₁/∂t + ∇·(α₁u) = (α₁ + K)∇·u,
//!   K = α₁α₂(ρ₂c₂² − ρ₁c₁²)/(α₂ρ₁c₁² + α₁ρ₂c₂²).
//! Each phase has its own NASG equation of state; the mixture pressure
//! follows in closed form from ρe = p·D + B with D, B linear in α_k and
//! α_kρ_k. With K the mixture sound speed is Wood's,
//!   1/(ρc²) = Σ α_k/(ρ_k c_k²),
//! the physical low-frequency sound speed of a bubbly liquid (it drops to
//! ~20 m/s for 50% air in water). `kapila = false` drops K (Allaire et al.
//! 2002), whose sound speed is monotone but not Wood's.
//!
//! The material interface is a diffuse zone a few cells wide in which α
//! goes from 0 to 1; the sensor flags it and the WCNS/HLLC path (α treated
//! like a partial density, Johnsen & Colonius 2006) keeps pressure and
//! velocity exactly uniform across a moving interface.

use crate::eos::Nasg;
use crate::system::{Mat, System};

type V = [f64; 7];
type M7 = Mat<7>;

#[derive(Clone, Copy, Debug)]
pub struct FiveEq {
    pub eos: [Nasg; 2],
    pub kapila: bool,
}

// conserved: a1, a2, ρu, ρv, ρw, E, α1 ; primitive: a1, a2, u, v, w, p, α1
const IA: usize = 6;
const IP: usize = 5;
const IU: usize = 2;

impl FiveEq {
    fn phase_properties(&self, w: &V) -> [(f64, f64, f64); 2] {
        let al = [w[IA], 1.0 - w[IA]];
        std::array::from_fn(|k| {
            let e = &self.eos[k];
            let rk = if al[k] > 1e-12 { w[k] / al[k] } else { 0.0 };
            let t = if rk > 0.0 { e.temperature(rk, w[IP]) } else { f64::NAN };
            (al[k], t, rk * e.cv)
        })
    }
    /// Volume fractions, *not* clamped: ρe = pD + B must stay exactly linear
    /// in α for pressure equilibrium to survive tiny over/undershoots of the
    /// high-order schemes (clamping broke it in the advection test).
    #[inline]
    fn alphas(w: &V) -> [f64; 2] {
        [w[IA], 1.0 - w[IA]]
    }

    /// D, B of ρe = p D + B.
    #[inline]
    fn db(&self, w: &V) -> (f64, f64) {
        let al = Self::alphas(w);
        let (mut d, mut b) = (0.0, 0.0);
        for k in 0..2 {
            let e = &self.eos[k];
            let vol = al[k] - e.b * w[k];
            d += vol / (e.gamma - 1.0);
            b += e.gamma * e.p_inf * vol / (e.gamma - 1.0) + e.q * w[k];
        }
        (d, b)
    }

    #[inline]
    fn rho_e(&self, w: &V) -> f64 {
        let (d, b) = self.db(w);
        w[IP] * d + b
    }

    /// Phase ρ_k c_k² (well defined as α_k → 0).
    #[inline]
    fn rho_c2_phase(&self, w: &V, k: usize) -> f64 {
        let al = Self::alphas(w);
        let e = &self.eos[k];
        let rho_k = if al[k] > 1e-12 && w[k] > 0.0 { w[k] / al[k] } else { 0.0 };
        e.gamma * (w[IP] + e.p_inf) / (1.0 - e.b * rho_k).max(1e-12)
    }

    /// Kapila coefficient K.
    pub fn kapila_k(&self, w: &V) -> f64 {
        if !self.kapila {
            return 0.0;
        }
        let al = [w[IA].clamp(0.0, 1.0), 1.0 - w[IA].clamp(0.0, 1.0)];
        let (z1, z2) = (self.rho_c2_phase(w, 0), self.rho_c2_phase(w, 1));
        let den = al[1] * z1 + al[0] * z2;
        if den <= 0.0 {
            return 0.0;
        }
        al[0] * al[1] * (z2 - z1) / den
    }

    /// Mixture ρc² of the model (Wood's with K).
    pub fn rho_c2(&self, w: &V) -> f64 {
        let (d, _) = self.db(w);
        let p = w[IP];
        let mut num = self.rho_e(w) + p;
        for k in 0..2 {
            let e = &self.eos[k];
            num += w[k] * (e.b * (p + e.gamma * e.p_inf) / (e.gamma - 1.0) - e.q);
        }
        num += self.kapila_k(w) * self.drhoe_dalpha(p);
        num / d
    }

    #[inline]
    fn drhoe_dalpha(&self, p: f64) -> f64 {
        let (e1, e2) = (&self.eos[0], &self.eos[1]);
        (p + e1.gamma * e1.p_inf) / (e1.gamma - 1.0) - (p + e2.gamma * e2.p_inf) / (e2.gamma - 1.0)
    }

    /// Mixture density.
    #[inline]
    pub fn rho(w: &V) -> f64 {
        w[0] + w[1]
    }

    /// Effective p∞ of the mixture (B/D-like), for the sensor's log scale.
    fn p_inf_mix(&self, w: &V) -> f64 {
        let a = w[IA].clamp(0.0, 1.0);
        a * self.eos[0].p_inf + (1.0 - a) * self.eos[1].p_inf
    }
}

impl System<7> for FiveEq {
    fn to_prim(&self, u: &V) -> V {
        let rho = u[0] + u[1];
        let (vx, vy, vz) = (u[2] / rho, u[3] / rho, u[4] / rho);
        let mut w = [0.0; 7];
        w[0] = u[0];
        w[1] = u[1];
        w[2] = vx;
        w[3] = vy;
        w[4] = vz;
        w[IA] = u[6];
        let rhoe = u[5] - 0.5 * rho * (vx * vx + vy * vy + vz * vz);
        let (d, b) = self.db(&w);
        w[IP] = (rhoe - b) / d;
        w
    }

    #[inline]
    fn velocity_conserved(&self, u: &V, d: usize) -> f64 {
        u[2 + d] / (u[0] + u[1])
    }
    fn to_cons(&self, w: &V) -> V {
        let rho = w[0] + w[1];
        let mut u = [0.0; 7];
        u[0] = w[0];
        u[1] = w[1];
        u[2] = rho * w[2];
        u[3] = rho * w[3];
        u[4] = rho * w[4];
        u[5] = self.rho_e(w) + 0.5 * rho * (w[2] * w[2] + w[3] * w[3] + w[4] * w[4]);
        u[6] = w[IA];
        u
    }

    fn sound_speed(&self, w: &V) -> f64 {
        (self.rho_c2(w).max(0.0) / Self::rho(w)).sqrt()
    }

    fn vel_index(&self) -> usize {
        IU
    }

    fn momentum_index(&self, d: usize) -> usize {
        2 + d
    }

    fn energy_index(&self) -> Option<usize> {
        Some(5)
    }

    fn density(&self, w: &V) -> f64 {
        Self::rho(w)
    }
    fn density_indices(&self) -> std::ops::Range<usize> {
        0..2
    }
    fn phases(&self, w: &V) -> Vec<(f64, f64, f64)> {
        self.phase_properties(w).to_vec()
    }
    fn conduction_properties(&self, w: &V, conductivity: &[f64]) -> (f64, f64, f64) {
        let phases = self.phase_properties(w);
        crate::thermal::properties(self.temperature_from_phases(w, &phases), &phases, conductivity)
    }
    fn temperature(&self, w: &V) -> f64 {
        crate::thermal::mixture_temperature(&self.phase_properties(w), &[w[0], w[1]])
    }
    fn temperature_from_phases(&self, w: &V, phases: &[(f64, f64, f64)]) -> f64 {
        crate::thermal::mixture_temperature(phases, &[w[0], w[1]])
    }
    fn phase_stiffness(&self, w: &V) -> Vec<(f64, f64)> {
        vec![(w[IA], self.rho_c2_phase(w, 0)), (1.0 - w[IA], self.rho_c2_phase(w, 1))]
    }

    fn eigen(&self, w: &V, d: usize, l: &mut M7, r: &mut M7) {
        let rho = Self::rho(w);
        let c = self.sound_speed(w);
        let kk = self.kapila_k(w);
        let un = IU + d;
        let (t1, t2) = (IU + (d + 1) % 3, IU + (d + 2) % 3);
        for row in l.iter_mut().chain(r.iter_mut()) {
            *row = [0.0; 7];
        }
        let c2 = c * c;
        // fields: 0 a1, 1 a2, 2 u−c, 3 t1, 4 t2, 5 u+c, 6 α
        l[0][0] = 1.0;
        l[0][IP] = -w[0] / rho / c2;
        l[1][1] = 1.0;
        l[1][IP] = -w[1] / rho / c2;
        l[2][un] = -rho / (2.0 * c);
        l[2][IP] = 1.0 / (2.0 * c2);
        l[3][t1] = 1.0;
        l[4][t2] = 1.0;
        l[5][un] = rho / (2.0 * c);
        l[5][IP] = 1.0 / (2.0 * c2);
        l[6][IA] = 1.0;
        l[6][IP] = kk / rho / c2;
        r[0][0] = 1.0;
        r[1][1] = 1.0;
        r[t1][3] = 1.0;
        r[t2][4] = 1.0;
        r[IA][6] = 1.0;
        for (col, sgn) in [(2usize, -1.0), (5usize, 1.0)] {
            r[0][col] = w[0] / rho;
            r[1][col] = w[1] / rho;
            r[un][col] = sgn * c / rho;
            r[IP][col] = c2;
            r[IA][col] = -kk / rho;
        }
    }

    fn eigenvalues(&self, w: &V, d: usize) -> V {
        let (u, c) = (w[IU + d], self.sound_speed(w));
        let mut l = [0.0; 7];
        l[..7].copy_from_slice(&[u, u, u - c, u, u, u + c, u]);
        l
    }

    fn riemann_flux(&self, wl: &V, wr: &V, d: usize) -> (V, f64) {
        let un = IU + d;
        let (cl, cr) = (self.sound_speed(wl), self.sound_speed(wr));
        let (ul, ur) = (wl[un], wr[un]);
        let sl = (ul - cl).min(ur - cr);
        let sr = (ul + cl).max(ur + cr);
        let (ql, qr) = (self.to_cons(wl), self.to_cons(wr));
        if sl >= 0.0 {
            return (self.flux(wl, d), ul);
        }
        if sr <= 0.0 {
            return (self.flux(wr, d), ur);
        }
        let (rl, rr) = (Self::rho(wl), Self::rho(wr));
        let (pl, pr) = (wl[IP], wr[IP]);
        let ss = (pr - pl + rl * ul * (sl - ul) - rr * ur * (sr - ur)) / (rl * (sl - ul) - rr * (sr - ur));
        let star = |w: &V, q: &V, s: f64| -> V {
            let (rho, u) = (Self::rho(w), w[un]);
            let fac = (s - u) / (s - ss);
            let mut qs = [0.0; 7];
            qs[0] = w[0] * fac;
            qs[1] = w[1] * fac;
            for k in IU..IU + 3 {
                qs[k] = rho * fac * w[k];
            }
            qs[un] = rho * fac * ss;
            qs[5] = rho * fac * (q[5] / rho + (ss - u) * (ss + w[IP] / (rho * (s - u))));
            qs[6] = w[IA] * fac;
            qs
        };
        let mut f = [0.0; 7];
        if ss >= 0.0 {
            let (fl, qs) = (self.flux(wl, d), star(wl, &ql, sl));
            for k in 0..7 {
                f[k] = fl[k] + sl * (qs[k] - ql[k]);
            }
        } else {
            let (fr, qs) = (self.flux(wr, d), star(wr, &qr, sr));
            for k in 0..7 {
                f[k] = fr[k] + sr * (qs[k] - qr[k]);
            }
        }
        (f, ss)
    }

    fn split_flux(&self, a: &V, b: &V, d: usize) -> V {
        let un = IU + d;
        let u = 0.5 * (a[un] + b[un]);
        let rho = 0.5 * (Self::rho(a) + Self::rho(b));
        let c = rho * u;
        let mut f = [0.0; 7];
        f[0] = 0.5 * (a[0] + b[0]) * u;
        f[1] = 0.5 * (a[1] + b[1]) * u;
        for k in IU..IU + 3 {
            f[k] = c * 0.5 * (a[k] + b[k]);
        }
        f[un] += 0.5 * (a[IP] + b[IP]);
        let kin = c * 0.5 * (a[2] * b[2] + a[3] * b[3] + a[4] * b[4]);
        let rhoe = 0.5 * (self.rho_e(a) + self.rho_e(b));
        let pw = 0.5 * (a[un] * b[IP] + b[un] * a[IP]);
        f[5] = kin + u * rhoe + pw;
        f[6] = 0.5 * (a[IA] + b[IA]) * u;
        f
    }

    fn flux(&self, w: &V, d: usize) -> V {
        let un = IU + d;
        let q = self.to_cons(w);
        let u = w[un];
        let mut f = [0.0; 7];
        for k in 0..5 {
            f[k] = q[k] * u;
        }
        f[un] += w[IP];
        f[5] = (q[5] + w[IP]) * u;
        f[6] = w[IA] * u;
        f
    }

    fn dudw(&self, w: &V) -> M7 {
        let rho = Self::rho(w);
        let p = w[IP];
        let (d, _) = self.db(w);
        let ke = 0.5 * (w[2] * w[2] + w[3] * w[3] + w[4] * w[4]);
        let mut j = [[0.0; 7]; 7];
        j[0][0] = 1.0;
        j[1][1] = 1.0;
        j[6][IA] = 1.0;
        for c in 0..3 {
            let k = IU + c;
            j[k][0] = w[k];
            j[k][1] = w[k];
            j[k][k] = rho;
            j[5][k] = rho * w[k];
        }
        for k in 0..2 {
            let e = &self.eos[k];
            j[5][k] = -e.b * (p + e.gamma * e.p_inf) / (e.gamma - 1.0) + e.q + ke;
        }
        j[5][IP] = d;
        j[5][IA] = self.drhoe_dalpha(p);
        j
    }

    fn admissible(&self, w: &V) -> bool {
        // what the physics needs is ρ > 0 and ρc² > 0; small over/undershoots
        // of α and α_kρ_k beyond the pure-phase limits are harmless because
        // the mixture EOS is used in its exact linear extension
        w[..7].iter().all(|v| v.is_finite()) && Self::rho(w) > 0.0 && (-0.05..=1.05).contains(&w[IA]) && self.rho_c2(w) > 0.0
    }

    fn nonconservative(&self) -> Option<usize> {
        Some(6)
    }

    fn nc_coef(&self, w: &V) -> f64 {
        w[IA] + self.kapila_k(w)
    }

    /// Dp/Dt = −ρc²∇·u must not overshoot the admissibility limit of the
    /// softest phase present (p > −p∞ of that phase; p > 0 for a vapour or
    /// gas): with α_v small, Wood's ρc² ≈ γ_v p/α_v makes this stiff.
    fn dilatation_rate(&self, w: &V, divu: f64) -> f64 {
        let a = w[IA];
        let mut pmin_inf = f64::INFINITY;
        for (k, al) in [(0usize, a), (1usize, 1.0 - a)] {
            if al > 1e-6 {
                pmin_inf = pmin_inf.min(self.eos[k].p_inf);
            }
        }
        // floor: a vapour or gas phase just nucleated in liquid under
        // tension sits at a negative mixture pressure (its sub-grid
        // interior is at p_v; the single-pressure closure cannot say so)
        let floor = 1e-6 * self.eos[0].p_inf.max(self.eos[1].p_inf);
        let head = (w[IP] + pmin_inf).max(floor).max(1e-300);
        self.rho_c2(w).max(0.0) * divu.abs() / head
    }

    fn sensor(&self, w: &V, out: &mut [f64; 8]) -> usize {
        out[0] = Self::rho(w).max(1e-300).ln();
        out[1] = (w[IP] + self.p_inf_mix(w)).max(1e-300).ln();
        out[2] = w[IA];
        3
    }
}
