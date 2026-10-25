//! Single-fluid compressible flow with a real-fluid equilibrium EOS
//! (`realfluid::RealFluidEos`: reference Helmholtz formulation, liquid,
//! vapour, supercritical and the homogeneous-equilibrium two-phase dome).
//! Cavitation and flashing need no extra model: the EOS puts the fluid at
//! saturation when the pressure falls to p_sat(T).
//!
//! Primitive variables w = (ρ, u, v, w, e) with e the specific internal
//! energy, so conversions never invert e(ρ, p) (two-valued for water near
//! its density maximum). Characteristic vectors in these variables:
//!   acoustic r± = (1, ±c/ρ, 0, 0, p/ρ²), entropy r₀ = (p_e, 0, 0, 0, −p_ρ),
//! with c² = p_ρ|_e + (p/ρ²) p_e|_ρ; left vectors by inversion.

use crate::realfluid::RealFluidEos;
use crate::system::{invert, Mat, System};
use std::sync::Arc;

type V = [f64; 5];
type M = Mat<5>;

#[derive(Clone)]
pub struct EulerReal {
    pub eos: Arc<RealFluidEos>,
}

impl EulerReal {
    #[inline]
    pub fn state(&self, w: &V) -> crate::realfluid::FluidState {
        self.eos.eval(w[0], w[4])
    }
    /// (p_ρ|_e, p_e|_ρ) by central differences of the table.
    fn dp(&self, w: &V) -> (f64, f64) {
        let (rho, e) = (w[0], w[4]);
        let hr = 1e-5 * rho;
        let he = 1e-5 * e.abs().max(1e4);
        let pr = (self.eos.eval(rho + hr, e).p - self.eos.eval(rho - hr, e).p) / (2.0 * hr);
        let pe = (self.eos.eval(rho, e + he).p - self.eos.eval(rho, e - he).p) / (2.0 * he);
        (pr, pe)
    }
}

impl System<5> for EulerReal {
    fn to_prim(&self, u: &V) -> V {
        let rho = u[0];
        let (vx, vy, vz) = (u[1] / rho, u[2] / rho, u[3] / rho);
        [rho, vx, vy, vz, u[4] / rho - 0.5 * (vx * vx + vy * vy + vz * vz)]
    }
    fn to_cons(&self, w: &V) -> V {
        let rho = w[0];
        [rho, rho * w[1], rho * w[2], rho * w[3], rho * (w[4] + 0.5 * (w[1] * w[1] + w[2] * w[2] + w[3] * w[3]))]
    }
    fn sound_speed(&self, w: &V) -> f64 {
        self.state(w).c
    }
    fn vel_index(&self) -> usize {
        1
    }
    fn momentum_index(&self, d: usize) -> usize {
        1 + d
    }
    fn energy_index(&self) -> Option<usize> {
        Some(4)
    }
    fn temperature(&self, w: &V) -> f64 {
        self.state(w).t
    }
    fn eigen(&self, w: &V, d: usize, l: &mut M, r: &mut M) {
        let st = self.state(w);
        let (rho, c, p) = (w[0], st.c, st.p);
        let (pr, pe) = self.dp(w);
        let un = 1 + d;
        let (t1, t2) = (1 + (d + 1) % 3, 1 + (d + 2) % 3);
        *r = [[0.0; 5]; 5];
        // columns: 0 u−c, 1 entropy, 2 t1, 3 t2, 4 u+c
        r[0][0] = 1.0;
        r[un][0] = -c / rho;
        r[4][0] = p / (rho * rho);
        let nrm = (pe * pe + pr * pr).sqrt().max(1e-300);
        r[0][1] = pe / nrm;
        r[4][1] = -pr / nrm;
        r[t1][2] = 1.0;
        r[t2][3] = 1.0;
        r[0][4] = 1.0;
        r[un][4] = c / rho;
        r[4][4] = p / (rho * rho);
        *l = invert(r);
    }
    fn eigenvalues(&self, w: &V, d: usize) -> V {
        let (u, c) = (w[1 + d], self.sound_speed(w));
        [u - c, u, u, u, u + c]
    }
    fn riemann_flux(&self, wl: &V, wr: &V, d: usize) -> (V, f64) {
        let un = 1 + d;
        let (sl_, sr_) = (self.state(wl), self.state(wr));
        let (cl, cr, pl, pr) = (sl_.c, sr_.c, sl_.p, sr_.p);
        let (ul, ur) = (wl[un], wr[un]);
        let sl = (ul - cl).min(ur - cr);
        let sr = (ul + cl).max(ur + cr);
        let flux = |w: &V, p: f64| -> V {
            let q = self.to_cons(w);
            let u = w[un];
            let mut f = [q[0] * u, q[1] * u, q[2] * u, q[3] * u, (q[4] + p) * u];
            f[un] += p;
            f
        };
        if sl >= 0.0 {
            return (flux(wl, pl), ul);
        }
        if sr <= 0.0 {
            return (flux(wr, pr), ur);
        }
        let (rl, rr) = (wl[0], wr[0]);
        let ss = (pr - pl + rl * ul * (sl - ul) - rr * ur * (sr - ur)) / (rl * (sl - ul) - rr * (sr - ur));
        let (w, p, s) = if ss >= 0.0 { (wl, pl, sl) } else { (wr, pr, sr) };
        let q = self.to_cons(w);
        let fx = flux(w, p);
        let (rho, u) = (w[0], w[un]);
        let fac = rho * (s - u) / (s - ss);
        let mut qs = [fac, fac * w[1], fac * w[2], fac * w[3], 0.0];
        qs[un] = fac * ss;
        qs[4] = fac * (q[4] / rho + (ss - u) * (ss + p / (rho * (s - u))));
        let mut f = [0.0; 5];
        for k in 0..5 {
            f[k] = fx[k] + s * (qs[k] - q[k]);
        }
        (f, ss)
    }
    fn split_flux(&self, a: &V, b: &V, d: usize) -> V {
        let un = 1 + d;
        let (pa, pb) = (self.state(a).p, self.state(b).p);
        let rho = 0.5 * (a[0] + b[0]);
        let u = 0.5 * (a[un] + b[un]);
        let c = rho * u;
        let mut f = [c, c * 0.5 * (a[1] + b[1]), c * 0.5 * (a[2] + b[2]), c * 0.5 * (a[3] + b[3]), 0.0];
        f[un] += 0.5 * (pa + pb);
        let kin = c * 0.5 * (a[1] * b[1] + a[2] * b[2] + a[3] * b[3]);
        let rhoe = 0.5 * (a[0] * a[4] + b[0] * b[4]);
        let pw = 0.5 * (a[un] * pb + b[un] * pa);
        f[4] = kin + u * rhoe + pw;
        f
    }
    fn flux(&self, w: &V, d: usize) -> V {
        let un = 1 + d;
        let p = self.state(w).p;
        let q = self.to_cons(w);
        let u = w[un];
        let mut f = [q[0] * u, q[1] * u, q[2] * u, q[3] * u, (q[4] + p) * u];
        f[un] += p;
        f
    }
    fn dudw(&self, w: &V) -> M {
        let rho = w[0];
        let ke = 0.5 * (w[1] * w[1] + w[2] * w[2] + w[3] * w[3]);
        let mut j = [[0.0; 5]; 5];
        j[0][0] = 1.0;
        for k in 1..4 {
            j[k][0] = w[k];
            j[k][k] = rho;
            j[4][k] = rho * w[k];
        }
        j[4][0] = w[4] + ke;
        j[4][4] = rho;
        j
    }
    fn admissible(&self, w: &V) -> bool {
        // the admissible set includes the EOS's validity domain: outside the
        // table the bicubic stencil extrapolates without bound
        if !(w.iter().all(|v| v.is_finite()) && w[0] > 0.0 && self.eos.in_domain(w[0], w[4])) {
            return false;
        }
        let st = self.state(w);
        st.p.is_finite() && st.c.is_finite() && st.c > 0.0
    }
    fn sensor(&self, w: &V, out: &mut [f64; 8]) -> usize {
        let st = self.state(w);
        out[0] = w[0].max(1e-300).ln();
        out[1] = st.p / (w[0] * st.c * st.c).max(1e-300);
        out[2] = st.c.max(1e-300).ln();
        3
    }
}
