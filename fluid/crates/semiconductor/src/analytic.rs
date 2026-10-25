//! Exact and classical references for the device solver.
//!
//! - MOS on uniformly doped p-type silicon, Boltzmann statistics: the 1D
//!   Poisson–Boltzmann equation integrates once in closed form,
//!   E_s = √(2qp₀V_t/ε) F(u, v),
//!   F² = (e^{−u} + u − 1) + (n₀/p₀) e^{−v}(e^{u} − u − 1),
//!   u = band bending/V_t, v = electron quasi-Fermi split/V_t. This gives
//!   the gate charge and low-frequency C–V exactly, and with the gradual-
//!   channel approximation the Pao–Sah drain current (the exact long-
//!   channel limit of drift–diffusion at constant mobility).
//! - Shockley's ideal diode and the depletion capacitance.

use crate::silicon::Q;

/// Uniform p-type MOS structure.
#[derive(Clone, Copy, Debug)]
pub struct Mos {
    /// acceptor density (m⁻³)
    pub na: f64,
    pub ni: f64,
    pub vt: f64,
    pub eps_s: f64,
    /// oxide capacitance per area (F/m²)
    pub cox: f64,
    /// flat-band voltage (V)
    pub vfb: f64,
}

impl Mos {
    /// Bulk majority and minority densities (exact neutrality).
    pub fn p0n0(&self) -> (f64, f64) {
        let psi_b = -(0.5 * self.na / self.ni).asinh();
        (self.ni * (-psi_b).exp(), self.ni * psi_b.exp())
    }
    pub fn f(&self, u: f64, v: f64) -> f64 {
        let (p0, n0) = self.p0n0();
        let a = (-u).exp() + u - 1.0;
        let b = n0 / p0 * (-v).exp() * (u.exp() - u - 1.0);
        (a + b).max(0.0).sqrt()
    }
    /// Semiconductor charge per area (C/m²) at band bending u (units of V_t).
    pub fn qs(&self, u: f64, v: f64) -> f64 {
        let (p0, _) = self.p0n0();
        -u.signum() * (2.0 * Q * p0 * self.eps_s * self.vt).sqrt() * self.f(u, v)
    }
    /// Gate voltage that produces band bending u.
    pub fn vg_of(&self, u: f64, v: f64) -> f64 {
        self.vfb + self.vt * u - self.qs(u, v) / self.cox
    }
    /// Band bending at gate voltage vg (bisection; vg is monotone in u).
    pub fn u_of(&self, vg: f64, v: f64) -> f64 {
        let (mut lo, mut hi) = (-60.0, 80.0);
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if self.vg_of(mid, v) < vg {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }
    /// Deep depletion (no minority carriers): semiconductor charge with the
    /// electron term removed.
    pub fn qs_deep(&self, u: f64) -> f64 {
        let (p0, _) = self.p0n0();
        -u.signum() * (2.0 * Q * p0 * self.eps_s * self.vt).sqrt() * ((-u).exp() + u - 1.0).max(0.0).sqrt()
    }
    /// Gate charge per area at vg in deep depletion.
    pub fn qg_deep(&self, vg: f64) -> f64 {
        let vg_of = |u: f64| self.vfb + self.vt * u - self.qs_deep(u) / self.cox;
        let (mut lo, mut hi) = (-60.0, 400.0);
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if vg_of(mid) < vg {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        -self.qs_deep(0.5 * (lo + hi))
    }
    /// Gate charge per area at vg (equilibrium).
    pub fn qg(&self, vg: f64) -> f64 {
        -self.qs(self.u_of(vg, 0.0), 0.0)
    }
    /// Inversion (electron) charge per area at gate voltage vg and channel
    /// quasi-Fermi split v (V): q∫ n dx = q∫₀^{u_s} n(u)/|du/dx| du.
    pub fn qi(&self, vg: f64, v: f64) -> f64 {
        let vv = v / self.vt;
        let us = self.u_of(vg, vv);
        if us <= 0.0 {
            return 0.0;
        }
        let (p0, n0) = self.p0n0();
        let k = (2.0 * Q * p0 / (self.eps_s * self.vt)).sqrt(); // |du/dx| = k F
        // integrand peaks at u_s (e^u); Gauss–Legendre panels of width ≤ 0.25
        let u0 = 1e-6 * us;
        let panels = ((us - u0) / 0.25).ceil().max(4.0) as usize;
        let (gx, gw) = gauss8();
        let mut s = 0.0;
        for pnl in 0..panels {
            let a = u0 + (us - u0) * pnl as f64 / panels as f64;
            let b = u0 + (us - u0) * (pnl + 1) as f64 / panels as f64;
            for (x, w) in gx.iter().zip(gw.iter()) {
                let u = 0.5 * (a + b) + 0.5 * (b - a) * x;
                let n = n0 * (u - vv).exp();
                s += 0.5 * (b - a) * w * n / (k * self.f(u, vv));
            }
        }
        Q * s
    }
    /// Pao–Sah drain current per unit width (A/m) for channel length l and
    /// mobility mu, source and body at 0.
    pub fn pao_sah(&self, vg: f64, vd: f64, l: f64, mu: f64) -> f64 {
        let panels = 16;
        let (gx, gw) = gauss8();
        let mut s = 0.0;
        for pnl in 0..panels {
            let a = vd * pnl as f64 / panels as f64;
            let b = vd * (pnl + 1) as f64 / panels as f64;
            for (x, w) in gx.iter().zip(gw.iter()) {
                let v = 0.5 * (a + b) + 0.5 * (b - a) * x;
                s += 0.5 * (b - a) * w * self.qi(vg, v);
            }
        }
        mu / l * s
    }
}

/// 8-point Gauss–Legendre nodes and weights on [−1, 1].
pub fn gauss8() -> ([f64; 8], [f64; 8]) {
    let x = [-0.960_289_856_497_536_2, -0.796_666_477_413_626_7, -0.525_532_409_916_329, -0.183_434_642_495_649_8, 0.183_434_642_495_649_8, 0.525_532_409_916_329, 0.796_666_477_413_626_7, 0.960_289_856_497_536_2];
    let w = [0.101_228_536_290_376_3, 0.222_381_034_453_374_5, 0.313_706_645_877_887_3, 0.362_683_783_378_362, 0.362_683_783_378_362, 0.313_706_645_877_887_3, 0.222_381_034_453_374_5, 0.101_228_536_290_376_3];
    (x, w)
}

/// Shockley long-base diode saturation current density (A/m²):
/// J₀ = q n_i² (D_n/(L_n N_a) + D_p/(L_p N_d)).
#[allow(clippy::too_many_arguments)]
pub fn shockley_j0(ni: f64, vt: f64, na: f64, nd: f64, mu_n_p_side: f64, tau_n: f64, mu_p_n_side: f64, tau_p: f64) -> f64 {
    let (dn, dp) = (mu_n_p_side * vt, mu_p_n_side * vt);
    let (ln, lp) = ((dn * tau_n).sqrt(), (dp * tau_p).sqrt());
    Q * ni * ni * (dn / (ln * na) + dp / (lp * nd))
}

/// Abrupt-junction depletion capacitance per area including the majority-
/// carrier tail correction (V_bi − V − 2V_t).
pub fn depletion_capacitance(eps: f64, na: f64, nd: f64, vbi: f64, v: f64, vt: f64) -> f64 {
    (Q * eps * na * nd / (2.0 * (na + nd) * (vbi - v - 2.0 * vt))).sqrt()
}
