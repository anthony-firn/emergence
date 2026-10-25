//! Exact Riemann solver for the 1D Euler equations with a stiffened-gas
//! material on each side (γ_K, p∞_K; b = q = 0), after Toro ch. 4 with
//! p → p + p∞ per side. Ideal gas on both sides is the classic problem;
//! different materials give the exact gas–liquid shock tube (the contact
//! is the material interface). Used only for validation.

use crate::eos::Nasg;

#[derive(Clone, Copy, Debug)]
pub struct Side {
    pub rho: f64,
    pub u: f64,
    pub p: f64,
    pub eos: Nasg,
}

#[derive(Clone, Copy, Debug)]
pub struct Star {
    pub p: f64,
    pub u: f64,
    pub rho_l: f64,
    pub rho_r: f64,
}

/// f_K(p) and its derivative (pressure function across the left/right wave).
fn f_k(p: f64, s: &Side) -> (f64, f64) {
    let (g, pi) = (s.eos.gamma, s.eos.p_inf);
    let c = (g * (s.p + pi) / s.rho).sqrt();
    let (pk, ps) = (s.p + pi, p + pi);
    if p > s.p {
        // shock
        let a = 2.0 / ((g + 1.0) * s.rho);
        let b = (g - 1.0) / (g + 1.0) * pk;
        let q = (a / (ps + b)).sqrt();
        ((p - s.p) * q, q * (1.0 - 0.5 * (p - s.p) / (ps + b)))
    } else {
        // rarefaction
        let r = ps / pk;
        let f = 2.0 * c / (g - 1.0) * (r.powf((g - 1.0) / (2.0 * g)) - 1.0);
        (f, 1.0 / (s.rho * c) * r.powf(-(g + 1.0) / (2.0 * g)))
    }
}

/// Star-region state by Newton iteration on f_L(p) + f_R(p) + Δu = 0.
pub fn star(l: &Side, r: &Side) -> Star {
    let du = r.u - l.u;
    let pmin = -(l.eos.p_inf.min(r.eos.p_inf)) + 1e-12 * (1.0 + l.p.abs() + r.p.abs());
    let mut p = (0.5 * (l.p + r.p)).max(pmin + 1e-8);
    for _ in 0..200 {
        let (fl, dl) = f_k(p, l);
        let (fr, dr) = f_k(p, r);
        let mut pn = p - (fl + fr + du) / (dl + dr);
        if pn <= pmin {
            pn = 0.5 * (p + pmin);
        }
        if (pn - p).abs() <= 1e-14 * (pn.abs() + l.eos.p_inf.max(r.eos.p_inf)) {
            p = pn;
            break;
        }
        p = pn;
    }
    let (fl, _) = f_k(p, l);
    let (fr, _) = f_k(p, r);
    let u = 0.5 * (l.u + r.u) + 0.5 * (fr - fl);
    let rho_star = |s: &Side| {
        let (g, pi) = (s.eos.gamma, s.eos.p_inf);
        let ratio = (p + pi) / (s.p + pi);
        if p > s.p {
            let gr = (g - 1.0) / (g + 1.0);
            s.rho * (ratio + gr) / (gr * ratio + 1.0)
        } else {
            s.rho * ratio.powf(1.0 / g)
        }
    };
    Star { p, u, rho_l: rho_star(l), rho_r: rho_star(r) }
}

/// Exact solution (ρ, u, p, side) at similarity coordinate ξ = x/t; side
/// is 0 left of the contact, 1 right of it.
pub fn sample(l: &Side, r: &Side, st: &Star, xi: f64) -> (f64, f64, f64, usize) {
    let left = xi <= st.u;
    let (s, sign) = if left { (l, 1.0) } else { (r, -1.0) };
    let rho_s = if left { st.rho_l } else { st.rho_r };
    let (g, pi) = (s.eos.gamma, s.eos.p_inf);
    let c = (g * (s.p + pi) / s.rho).sqrt();
    // work in the frame mirrored so the wave of interest is a left wave
    let (u, x) = (sign * s.u, sign * xi);
    let us = sign * st.u;
    let side = if left { 0 } else { 1 };
    if st.p > s.p {
        let ratio = (st.p + pi) / (s.p + pi);
        let sh = u - c * ((g + 1.0) / (2.0 * g) * ratio + (g - 1.0) / (2.0 * g)).sqrt();
        if x <= sh { (s.rho, s.u, s.p, side) } else { (rho_s, st.u, st.p, side) }
    } else {
        let cs = c * ((st.p + pi) / (s.p + pi)).powf((g - 1.0) / (2.0 * g));
        let (head, tail) = (u - c, us - cs);
        if x <= head {
            (s.rho, s.u, s.p, side)
        } else if x >= tail {
            (rho_s, st.u, st.p, side)
        } else {
            let k = 2.0 / (g + 1.0) + (g - 1.0) / ((g + 1.0) * c) * (u - x);
            let rho = s.rho * k.powf(2.0 / (g - 1.0));
            let uu = 2.0 / (g + 1.0) * (c + (g - 1.0) / 2.0 * u + x);
            let p = (s.p + pi) * k.powf(2.0 * g / (g - 1.0)) - pi;
            (rho, sign * uu, p, side)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Toro's Table 4.2/4.3 star values for his Tests 1–5 (γ = 1.4).
    #[test]
    fn toro_star_values() {
        let gas = Nasg::ideal(1.4);
        let cases = [
            ((1.0, 0.0, 1.0), (0.125, 0.0, 0.1), 0.30313, 0.92745),
            ((1.0, -2.0, 0.4), (1.0, 2.0, 0.4), 0.00189, 0.0),
            ((1.0, 0.0, 1000.0), (1.0, 0.0, 0.01), 460.894, 19.5975),
            ((1.0, 0.0, 0.01), (1.0, 0.0, 100.0), 46.0950, -6.19633),
            ((5.99924, 19.5975, 460.894), (5.99242, -6.19633, 46.0950), 1691.64, 8.68975),
        ];
        for (a, b, ps, us) in cases {
            let l = Side { rho: a.0, u: a.1, p: a.2, eos: gas };
            let r = Side { rho: b.0, u: b.1, p: b.2, eos: gas };
            let s = star(&l, &r);
            assert!((s.p - ps).abs() <= 1e-4 * ps.abs().max(1.0) + 1e-5, "p* {} vs {ps}", s.p);
            assert!((s.u - us).abs() <= 1e-4 * us.abs().max(1.0), "u* {} vs {us}", s.u);
        }
    }
}
