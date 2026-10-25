//! Hertz–Mindlin soft-sphere contact law with viscous damping and Coulomb
//! friction (the standard DEM contact model).
//!
//! **Normal.** The elastic part is Hertzian, `F_n = (4/3) E* sqrt(R*)
//! delta^(3/2)`, with `1/E* = (1-nu_1^2)/E_1 + (1-nu_2^2)/E_2` and
//! `1/R* = 1/R_1 + 1/R_2`. The damping follows Tsuji, Tanaka & Ishida
//! (1992): `F_n^d = -eta_n v_n`, `eta_n = 2 sqrt(5/6) beta sqrt(S_n m*)`,
//! `S_n = 2 E* sqrt(R* delta)`, `beta = -ln(e)/sqrt(ln^2 e + pi^2)`.
//!
//! **Tangential.** A Cundall–Strack tangential spring
//! `F_t = -k_t xi - c_t v_t` with `k_t = S_t = 8 G* sqrt(R* delta)` and
//! `1/G* = 2(2-nu_1)(1+nu_1)/E_1 + 2(2-nu_2)(1+nu_2)/E_2`. Sliding is
//! regularised by the Coulomb limit `|F_t| <= mu F_n`; when it is exceeded
//! the spring is reset to the limit.
//!
//! **Rolling.** An optional rolling-resistance torque
//! `tau_r = -mu_r R* |F_n| (omega_i - omega_j)/|omega_i - omega_j|`, which
//! is what keeps a pile of spheres from flattening (Zhou et al. 1999;
//! Ai et al. 2011).
//!
//! Sources: P. A. Cundall & O. D. L. Strack (1979) *Géotechnique* **29**,
//! 47; Y. Tsuji, T. Tanaka & T. Ishida (1992) *Powder Technol.* **71**, 239;
//! R. D. Mindlin (1949) *J. Appl. Mech.* **16**, 259.

/// Elastic and frictional parameters shared by all contacts.
#[derive(Clone, Copy, Debug)]
pub struct ContactParams {
    /// Young's modulus `E`.
    pub young: f64,
    /// Poisson ratio `nu`.
    pub poisson: f64,
    /// Interparticle Coulomb friction coefficient `mu`.
    pub mu: f64,
    /// Rolling friction coefficient `mu_r` (0 disables rolling resistance).
    pub mu_rolling: f64,
    /// Target normal coefficient of restitution `e` in `(0, 1]`.
    pub restitution: f64,
}

impl ContactParams {
    /// Effective Young's modulus for two identical particles.
    pub fn e_star(&self) -> f64 {
        self.young / (2.0 * (1.0 - self.poisson * self.poisson))
    }

    /// Effective shear modulus for two identical particles.
    pub fn g_star(&self) -> f64 {
        self.young / (4.0 * (2.0 - self.poisson) * (1.0 + self.poisson))
    }

    /// Positive Tsuji damping parameter `beta = -ln e / sqrt(ln^2 e + pi^2)`.
    pub fn beta(&self) -> f64 {
        let e = self.restitution.clamp(1e-6, 1.0);
        if e >= 1.0 {
            return 0.0;
        }
        let le = e.ln();
        -le / (le * le + std::f64::consts::PI * std::f64::consts::PI).sqrt()
    }
}

/// Persistent tangential state of one contact.
#[derive(Clone, Copy, Debug, Default)]
pub struct ContactState {
    /// Tangential spring displacement.
    pub xi: [f64; 3],
}

/// Forces and torques produced by one contact, plus the updated spring.
#[derive(Clone, Copy, Debug)]
pub struct Contact {
    /// Force on particle `i` from particle `j`.
    pub force_i: [f64; 3],
    /// Torque on particle `i`.
    pub torque_i: [f64; 3],
    /// Torque on particle `j`.
    pub torque_j: [f64; 3],
    /// Normal force magnitude (positive in compression).
    pub fn_mag: f64,
    /// Tangential force magnitude.
    pub ft_mag: f64,
    /// Updated tangential spring displacement.
    pub state: ContactState,
}

/// Evaluate one Hertz–Mindlin contact. `xi`, `xj` are the (minimum-image,
/// Lees–Edwards-unwrapped) centres, `vi`, `vj` the velocities, `wi`, `wj`
/// the angular velocities. `dt` advances the tangential spring. Returns
/// `None` if the spheres are not overlapping. `xi_state` is updated in
/// place.
#[allow(clippy::too_many_arguments)]
pub fn hertz_mindlin(
    p: &ContactParams,
    ri: f64,
    rj: f64,
    mi: f64,
    mj: f64,
    xi_pos: [f64; 3],
    vi: [f64; 3],
    wi: [f64; 3],
    xj_pos: [f64; 3],
    vj: [f64; 3],
    wj: [f64; 3],
    xi_state: &mut ContactState,
    dt: f64,
) -> Option<Contact> {
    let dr = [xj_pos[0] - xi_pos[0], xj_pos[1] - xi_pos[1], xj_pos[2] - xi_pos[2]];
    let dist = (dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]).sqrt();
    let delta = ri + rj - dist;
    if delta <= 0.0 || dist < 1e-300 {
        // Contact broken: forget the tangential history.
        *xi_state = ContactState::default();
        return None;
    }
    let n = [dr[0] / dist, dr[1] / dist, dr[2] / dist];

    let r_star = ri * rj / (ri + rj);
    let m_star = mi * mj / (mi + mj);
    let e_star = p.e_star();
    let g_star = p.g_star();
    let sqrt_r = r_star.sqrt();
    let sqrt_d = delta.sqrt();

    // Normal.
    let f_elastic = (4.0 / 3.0) * e_star * sqrt_r * delta * sqrt_d;
    let s_n = 2.0 * e_star * sqrt_r * sqrt_d;
    let beta = p.beta();
    let eta_n = 2.0 * (5.0 / 6.0_f64).sqrt() * beta * (s_n * m_star).sqrt();

    let v_rel = [vj[0] - vi[0], vj[1] - vi[1], vj[2] - vi[2]];
    let v_n = v_rel[0] * n[0] + v_rel[1] * n[1] + v_rel[2] * n[2];

    let fn_mag = f_elastic - eta_n * v_n;

    // Tangential: surface relative velocity at the contact.
    // v_s = v_rel - (r_i w_i + r_j w_j) x n.
    let wx = [
        ri * wi[0] + rj * wj[0],
        ri * wi[1] + rj * wj[1],
        ri * wi[2] + rj * wj[2],
    ];
    let wxn = [
        wx[1] * n[2] - wx[2] * n[1],
        wx[2] * n[0] - wx[0] * n[2],
        wx[0] * n[1] - wx[1] * n[0],
    ];
    let vt = [v_rel[0] - wxn[0], v_rel[1] - wxn[1], v_rel[2] - wxn[2]];
    // Project out the normal component.
    let vt_n = vt[0] * n[0] + vt[1] * n[1] + vt[2] * n[2];
    let vt = [vt[0] - vt_n * n[0], vt[1] - vt_n * n[1], vt[2] - vt_n * n[2]];

    // Spring update.
    let mut xi = xi_state.xi;
    for k in 0..3 {
        xi[k] += vt[k] * dt;
    }
    let xi_n = xi[0] * n[0] + xi[1] * n[1] + xi[2] * n[2];
    for k in 0..3 {
        xi[k] -= xi_n * n[k];
    }

    let s_t = 8.0 * g_star * sqrt_r * sqrt_d;
    let k_t = s_t;
    let eta_t = 2.0 * (5.0 / 6.0_f64).sqrt() * beta * (s_t * m_star).sqrt();
    let mut ft = [0.0; 3];
    for k in 0..3 {
        ft[k] = -k_t * xi[k] - eta_t * vt[k];
    }
    let mut ft_mag = (ft[0] * ft[0] + ft[1] * ft[1] + ft[2] * ft[2]).sqrt();
    let ft_limit = p.mu * fn_mag.max(0.0);
    if ft_mag > ft_limit && ft_mag > 0.0 {
        let scale = ft_limit / ft_mag;
        for k in 0..3 {
            ft[k] *= scale;
        }
        // Reset the spring so the elastic force sits on the yield surface.
        if k_t > 0.0 {
            for k in 0..3 {
                xi[k] = -ft[k] / k_t;
            }
        }
        ft_mag = ft_limit;
    }

    // Force on i is the negative of the force on j.
    let force_i = [
        -(fn_mag * n[0] + ft[0]),
        -(fn_mag * n[1] + ft[1]),
        -(fn_mag * n[2] + ft[2]),
    ];

    // Torques: tau = r_contact x F.
    let mut torque_i = [
        -ri * (n[1] * ft[2] - n[2] * ft[1]),
        -ri * (n[2] * ft[0] - n[0] * ft[2]),
        -ri * (n[0] * ft[1] - n[1] * ft[0]),
    ];
    let mut torque_j = [
        -rj * (n[1] * ft[2] - n[2] * ft[1]),
        -rj * (n[2] * ft[0] - n[0] * ft[2]),
        -rj * (n[0] * ft[1] - n[1] * ft[0]),
    ];

    // Rolling resistance.
    if p.mu_rolling > 0.0 {
        let dw = [wi[0] - wj[0], wi[1] - wj[1], wi[2] - wj[2]];
        let dw_mag = (dw[0] * dw[0] + dw[1] * dw[1] + dw[2] * dw[2]).sqrt();
        if dw_mag > 1e-300 {
            let tmag = p.mu_rolling * r_star * fn_mag.abs();
            for k in 0..3 {
                torque_i[k] -= tmag * dw[k] / dw_mag;
                torque_j[k] += tmag * dw[k] / dw_mag;
            }
        }
    }

    *xi_state = ContactState { xi };
    Some(Contact {
        force_i,
        torque_i,
        torque_j,
        fn_mag,
        ft_mag,
        state: ContactState { xi },
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn no_contact_when_separated() {
        let p = ContactParams {
            young: 1000.0,
            poisson: 0.3,
            mu: 0.4,
            mu_rolling: 0.0,
            restitution: 0.5,
        };
        let mut st = ContactState::default();
        let c = hertz_mindlin(
            &p,
            0.5,
            0.5,
            1.0,
            1.0,
            [0.0; 3],
            [0.0; 3],
            [0.0; 3],
            [1.5, 0.0, 0.0],
            [0.0; 3],
            [0.0; 3],
            &mut st,
            1e-3,
        );
        assert!(c.is_none());
    }

    #[test]
    fn static_contact_force_matches_hertz() {
        let p = ContactParams {
            young: 1000.0,
            poisson: 0.3,
            mu: 0.4,
            mu_rolling: 0.0,
            restitution: 1.0, // no damping
        };
        let mut st = ContactState::default();
        let (ri, rj) = (0.5, 0.5);
        let delta = 0.01;
        let xj = [ri + rj - delta, 0.0, 0.0];
        let c = hertz_mindlin(
            &p,
            ri,
            rj,
            1.0,
            1.0,
            [0.0; 3],
            [0.0; 3],
            [0.0; 3],
            xj,
            [0.0; 3],
            [0.0; 3],
            &mut st,
            1e-4,
        )
        .unwrap();
        let r_star = 0.25_f64;
        let expected = (4.0 / 3.0) * p.e_star() * r_star.sqrt() * delta.powf(1.5);
        assert!((c.fn_mag - expected).abs() < 1e-12 * expected, "{} vs {expected}", c.fn_mag);
        // Force on i points in -x (away from j).
        assert!(c.force_i[0] < 0.0);
    }
}
