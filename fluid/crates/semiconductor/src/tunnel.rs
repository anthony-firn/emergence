//! Gate tunnelling through the oxide (Tsu–Esaki current, WKB transmission).
//!
//!   J = (q m_s kT / (2π² ħ³)) ∫ T(E) ln[(1 + e^{(E_F,s − E)/kT}) /
//!                                         (1 + e^{(E_F,g − E)/kT})] dE
//!   T(E) = exp(−2 ∫ κ dz),  κ = √(2 m_ox (E_ox(z) − E))/ħ  where E_ox > E
//!
//! with the oxide conduction-band edge E_ox(z) = Φ_B − qF z (energies from
//! the silicon conduction-band edge at the interface; trapezoidal for
//! direct tunnelling, triangular beyond qV_ox > Φ_B: Fowler–Nordheim). The
//! WKB integral over a linear barrier is closed-form. Quantization of the
//! inversion layer and image-force lowering are not included (the supply
//! is a 3D electron gas).
//!
//! Fold: the gate current is a post-processing of the drift–diffusion
//! solution. The indicator is I_g/I_d; above ~1e-3 the gate current would
//! have to enter the device equations (unfold to self-consistency).

use crate::silicon::{KB, Q};

pub const HBAR: f64 = 1.054_571_817e-34;
pub const M0: f64 = 9.109_383_7e-31;

/// WKB exponent 2∫κ dz for a linear barrier top U(z) = φ_b − qF z (J) over
/// [0, t] at energy E (J), mass m: closed form over the classically
/// forbidden part.
pub fn wkb_exponent(phi_b: f64, field: f64, t: f64, e: f64, m: f64) -> f64 {
    let qf = Q * field;
    let u0 = phi_b - e; // barrier height above E at z = 0
    if u0 <= 0.0 {
        return 0.0;
    }
    let u1 = phi_b - qf * t - e; // at z = t
    let k = 2.0 * (2.0 * m).sqrt() / HBAR;
    if qf.abs() < 1e-30 {
        return k * u0.sqrt() * t;
    }
    // ∫ √(u0 − qF z) dz from 0 to z_end = (2/(3qF)) (u0^{3/2} − u_end^{3/2})
    let u_end = u1.max(0.0);
    k * 2.0 / (3.0 * qf) * (u0.powf(1.5) - u_end.powf(1.5))
}

/// Tsu–Esaki current density (A/m², positive from the semiconductor into
/// the gate) through an oxide of thickness t with field F (V/m, positive
/// lowering the barrier toward the gate), barrier φ_b (eV), Fermi levels
/// (eV, relative to the semiconductor band edge at the interface),
/// temperature T, supply mass m_s and oxide mass m_ox (units of m₀).
#[allow(clippy::too_many_arguments)]
pub fn tsu_esaki(t_ox: f64, field: f64, phi_b_ev: f64, ef_s_ev: f64, ef_g_ev: f64, temp: f64, m_s: f64, m_ox: f64) -> f64 {
    let kt = KB * temp;
    let phi_b = phi_b_ev * Q;
    let (efs, efg) = (ef_s_ev * Q, ef_g_ev * Q);
    // integrate E from the band edge (0) to well above both Fermi levels
    let e_hi = (efs.max(efg) + 40.0 * kt).max(40.0 * kt);
    let n = 4000;
    let mut s = 0.0;
    let (gx, gw) = crate::analytic::gauss8();
    for p in 0..n {
        let a = e_hi * p as f64 / n as f64;
        let b = e_hi * (p + 1) as f64 / n as f64;
        for (x, w) in gx.iter().zip(gw.iter()) {
            let e = 0.5 * (a + b) + 0.5 * (b - a) * x;
            let tr = (-wkb_exponent(phi_b, field, t_ox, e, m_ox * M0)).exp();
            let l = |ef: f64| {
                let z = (ef - e) / kt;
                if z > 40.0 { z } else { z.exp().ln_1p() }
            };
            s += 0.5 * (b - a) * w * tr * (l(efs) - l(efg));
        }
    }
    Q * m_s * M0 * kt / (2.0 * std::f64::consts::PI.powi(2) * HBAR.powi(3)) * s
}

/// Fowler–Nordheim current density (A/m²) for barrier φ_b (eV), field F
/// (V/m), oxide mass m_ox (m₀): J = q³ F² m₀/(8π h m_ox φ_b)
/// · exp(−4 √(2 m_ox) φ_b^{3/2}/(3 ħ q F)).
pub fn fowler_nordheim(phi_b_ev: f64, field: f64, m_ox: f64) -> f64 {
    let h = 2.0 * std::f64::consts::PI * HBAR;
    let phi = phi_b_ev * Q;
    let a = Q.powi(3) / (8.0 * std::f64::consts::PI * h * phi * m_ox);
    let b = 4.0 * (2.0 * m_ox * M0).sqrt() * phi.powf(1.5) / (3.0 * HBAR * Q);
    a * field * field * (-b / field).exp()
}
