//! Hyperbolic systems solved by the finite-volume scheme. Each system works
//! in primitive variables w for reconstruction (its characteristic
//! decomposition needs only the sound speed, so any equation of state
//! works) and supplies two numerical fluxes: an upwind Riemann flux for
//! shocks and a symmetric two-point split flux for the non-dissipative
//! central scheme used where the flow is smooth.

use crate::eos::Nasg;

/// State vector and square matrix of an N-variable system.
pub type V<const N: usize> = [f64; N];
pub type Mat<const N: usize> = [[f64; N]; N];

/// Product of n×n blocks.
pub fn matmul<const N: usize>(a: &Mat<N>, b: &Mat<N>) -> Mat<N> {
    let n = N;
    let mut c = [[0.0; N]; N];
    for i in 0..n {
        for k in 0..n {
            let aik = a[i][k];
            if aik != 0.0 {
                for j in 0..n {
                    c[i][j] += aik * b[k][j];
                }
            }
        }
    }
    c
}

/// Inverse of an n×n block (Gauss–Jordan, partial pivoting).
pub fn invert<const N: usize>(a: &Mat<N>) -> Mat<N> {
    let n = N;
    let mut m = *a;
    let mut inv = [[0.0; N]; N];
    for i in 0..n {
        inv[i][i] = 1.0;
    }
    for c in 0..n {
        let piv = (c..n).max_by(|&x, &y| m[x][c].abs().partial_cmp(&m[y][c].abs()).unwrap()).unwrap();
        m.swap(c, piv);
        inv.swap(c, piv);
        let d = 1.0 / m[c][c];
        for j in 0..n {
            m[c][j] *= d;
            inv[c][j] *= d;
        }
        for r in 0..n {
            if r != c {
                let f = m[r][c];
                if f != 0.0 {
                    for j in 0..n {
                        m[r][j] -= f * m[c][j];
                        inv[r][j] -= f * inv[c][j];
                    }
                }
            }
        }
    }
    inv
}

pub trait System<const N: usize>: Sync + Send {
    fn to_prim(&self, u: &V<N>) -> V<N>;
    fn to_cons(&self, w: &V<N>) -> V<N>;
    fn sound_speed(&self, w: &V<N>) -> f64;
    /// One velocity component without reconstructing unused thermodynamics.
    /// Custom systems retain the primitive conversion unless they specialize it.
    fn velocity_conserved(&self, u: &V<N>, d: usize) -> f64 {
        self.to_prim(u)[self.vel_index() + d]
    }

    /// Index of the x-velocity in the primitive vector (y, z follow).
    fn vel_index(&self) -> usize;
    /// Index of the d-momentum in the conserved vector (for reflection).
    fn momentum_index(&self, d: usize) -> usize;
    /// Left (rows) and right (columns) eigenvectors of the primitive
    /// quasi-linear system in direction d at state w.
    fn eigen(&self, w: &V<N>, d: usize, l: &mut Mat<N>, r: &mut Mat<N>);
    /// Upwind (Riemann-solver) flux from primitive face states, and the
    /// face normal velocity (HLLC contact speed) used by non-conservative
    /// terms.
    fn riemann_flux(&self, wl: &V<N>, wr: &V<N>, d: usize) -> (V<N>, f64);
    /// Conserved index of a variable obeying ∂q/∂t + ∇·(q u) = c(w) ∇·u
    /// (e.g. a volume fraction), with c from `nc_coef`.
    fn nonconservative(&self) -> Option<usize> {
        None
    }
    fn nc_coef(&self, _w: &V<N>) -> f64 {
        0.0
    }
    /// All conserved indices obeying ∂q/∂t + ∇·(q u) = c_q(w) ∇·u
    /// (several volume fractions for more than two phases); defaults to
    /// `nonconservative()`.
    /// Conserved indices of the (partial) densities, summing to ρ.
    fn density_indices(&self) -> std::ops::Range<usize> {
        0..1
    }
    fn nc_vars(&self) -> std::ops::Range<usize> {
        match self.nonconservative() {
            Some(q) => q..q + 1,
            None => 0..0,
        }
    }
    /// Coefficient c_q of variable q in `nc_vars`.
    fn nc_coef_k(&self, w: &V<N>, _q: usize) -> f64 {
        self.nc_coef(w)
    }
    /// General non-conservative products: when true, the scheme also returns
    /// face values of the primitive variables (same combination as the
    /// fluxes) and the solver adds `nc_rhs(w_i, ∂_d w|_i, d)` to the
    /// conserved right-hand side.
    fn nc_general(&self) -> bool {
        false
    }
    fn nc_rhs(&self, _w: &V<N>, _dw: &V<N>, _d: usize) -> V<N> {
        [0.0; N]
    }
    /// Symmetric, consistent two-point flux F̃(a, b) used by the high-order
    /// split-form central scheme.
    fn split_flux(&self, wa: &V<N>, wb: &V<N>, d: usize) -> V<N>;
    fn admissible(&self, w: &V<N>) -> bool;
    /// Physical flux in direction d.
    fn flux(&self, w: &V<N>, d: usize) -> V<N>;
    /// Jacobian ∂u/∂w (conserved w.r.t. primitive).
    fn dudw(&self, w: &V<N>) -> Mat<N>;
    /// Values watched by the smoothness sensor (log-scaled where a relative
    /// measure is meant, e.g. ln ρ, ln(p + p∞)); returns how many.
    fn sensor(&self, w: &V<N>, out: &mut [f64; 8]) -> usize;
    /// Mixture density.
    fn density(&self, w: &V<N>) -> f64 {
        w[0]
    }
    /// Relative rate |Dp/Dt|/(p − p_min) produced by a dilatation rate
    /// `divu` (for the time-step limit of stiff non-conservative terms);
    /// 0 when the CFL condition already covers it.
    fn dilatation_rate(&self, _w: &V<N>, _divu: f64) -> f64 {
        0.0
    }
    /// Temperature (K or model units) from primitive state.
    fn temperature(&self, _w: &V<N>) -> f64 {
        f64::NAN
    }

    /// Temperature using phase properties already evaluated at this state.
    /// The default preserves custom temperature closures.
    fn temperature_from_phases(&self, w: &V<N>, _phases: &[(f64, f64, f64)]) -> f64 {
        self.temperature(w)
    }
    /// Temperature, conductivity, and volumetric heat capacity for conduction.
    /// The default preserves custom phase and temperature closures; models may
    /// override this to evaluate fixed-size phase data without allocations.
    fn conduction_properties(&self, w: &V<N>, conductivity: &[f64]) -> (f64, f64, f64) {
        let phases = self.phases(w);
        crate::thermal::properties(self.temperature_from_phases(w, &phases), &phases, conductivity)
    }
    /// Specific heat at constant pressure (for κ = μ c_p / Pr).
    fn cp(&self, _w: &V<N>) -> f64 {
        f64::NAN
    }
    /// Specific heat at fixed density, for the explicit viscous thermal limit.
    /// Custom systems can supply this directly; the default uses their phase
    /// heat capacities at fixed volume.
    fn cv(&self, w: &V<N>) -> f64 {
        self.phases(w).iter().map(|&(alpha, _, capacity)| alpha * capacity).sum::<f64>()
            / self.density(w)
    }
    /// Phases present in the cell: (volume fraction, temperature, ρ_k c_v,k)
    /// (heat conduction: heating at fixed volume raises T at 1/(ρc_v)).
    /// Single-phase systems: one entry.
    fn phases(&self, w: &V<N>) -> Vec<(f64, f64, f64)> {
        vec![(1.0, self.temperature(w), self.density(w) * self.cp(w) / 1.4)]
    }
    /// Phases present: (volume fraction, ρ_k c_k²) (stiffness contrast).
    fn phase_stiffness(&self, w: &V<N>) -> Vec<(f64, f64)> {
        vec![(1.0, self.density(w) * self.sound_speed(w).powi(2))]
    }
    /// Conserved index of total energy, if the system has one.
    fn energy_index(&self) -> Option<usize> {
        None
    }
    /// Characteristic speeds in direction d, in the order of `eigen`'s fields.
    fn eigenvalues(&self, w: &V<N>, d: usize) -> V<N>;
    /// Largest |u_d| + c.
    fn max_speed(&self, w: &V<N>, d: usize) -> f64 {
        w[self.vel_index() + d].abs() + self.sound_speed(w)
    }
}

/// Single-material Euler equations: w = (ρ, u, v, w, p), u = (ρ, ρu, ρv, ρw, E).
#[derive(Clone, Copy, Debug)]
pub struct Euler {
    pub eos: Nasg,
}

impl System<5> for Euler {
    #[inline]
    fn to_prim(&self, u: &V<5>) -> V<5> {
        let rho = u[0];
        let (vx, vy, vz) = (u[1] / rho, u[2] / rho, u[3] / rho);
        let e = u[4] / rho - 0.5 * (vx * vx + vy * vy + vz * vz);
        [rho, vx, vy, vz, self.eos.pressure(rho, e)]
    }

    #[inline]
    fn velocity_conserved(&self, u: &V<5>, d: usize) -> f64 {
        u[1 + d] / u[0]
    }

    #[inline]
    fn to_cons(&self, w: &V<5>) -> V<5> {
        let rho = w[0];
        let ke = 0.5 * rho * (w[1] * w[1] + w[2] * w[2] + w[3] * w[3]);
        let mut u = [0.0; 5];
        u[..5].copy_from_slice(&[rho, rho * w[1], rho * w[2], rho * w[3], self.eos.rho_e(rho, w[4]) + ke]);
        u
    }

    #[inline]
    fn sound_speed(&self, w: &V<5>) -> f64 {
        self.eos.sound_speed_sq(w[0], w[4]).max(0.0).sqrt()
    }

    fn vel_index(&self) -> usize {
        1
    }

    fn momentum_index(&self, d: usize) -> usize {
        1 + d
    }

    fn eigen(&self, w: &V<5>, d: usize, l: &mut Mat<5>, r: &mut Mat<5>) {
        let (rho, c) = (w[0], self.sound_speed(w));
        let un = 1 + d;
        let (t1, t2) = (1 + (d + 1) % 3, 1 + (d + 2) % 3);
        for row in l.iter_mut().chain(r.iter_mut()) {
            *row = [0.0; 5];
        }
        // waves: 0 = u − c, 1 = entropy, 2, 3 = shear, 4 = u + c
        l[0][un] = -rho / (2.0 * c);
        l[0][4] = 1.0 / (2.0 * c * c);
        l[1][0] = 1.0;
        l[1][4] = -1.0 / (c * c);
        l[2][t1] = 1.0;
        l[3][t2] = 1.0;
        l[4][un] = rho / (2.0 * c);
        l[4][4] = 1.0 / (2.0 * c * c);
        // right eigenvectors as columns
        r[0][0] = 1.0;
        r[un][0] = -c / rho;
        r[4][0] = c * c;
        r[0][1] = 1.0;
        r[t1][2] = 1.0;
        r[t2][3] = 1.0;
        r[0][4] = 1.0;
        r[un][4] = c / rho;
        r[4][4] = c * c;
    }

    /// HLLC (Toro, Spruce & Speares) with Davis wave-speed bounds.
    fn riemann_flux(&self, wl: &V<5>, wr: &V<5>, d: usize) -> (V<5>, f64) {
        let un = 1 + d;
        let (cl, cr) = (self.sound_speed(wl), self.sound_speed(wr));
        let (ul, ur) = (wl[un], wr[un]);
        let sl = (ul - cl).min(ur - cr);
        let sr = (ul + cl).max(ur + cr);
        let (ql, qr) = (self.to_cons(wl), self.to_cons(wr));
        let flux = |w: &V<5>, q: &V<5>| -> V<5> {
            let u = w[un];
            let mut f = [0.0; 5];
            f[0] = q[0] * u;
            f[1] = q[1] * u;
            f[2] = q[2] * u;
            f[3] = q[3] * u;
            f[un] += w[4];
            f[4] = (q[4] + w[4]) * u;
            f
        };
        if sl >= 0.0 {
            return (flux(wl, &ql), ul);
        }
        if sr <= 0.0 {
            return (flux(wr, &qr), ur);
        }
        let (rl, rr, pl, pr) = (wl[0], wr[0], wl[4], wr[4]);
        let ss = (pr - pl + rl * ul * (sl - ul) - rr * ur * (sr - ur)) / (rl * (sl - ul) - rr * (sr - ur));
        let star = |w: &V<5>, q: &V<5>, s: f64| -> V<5> {
            let (rho, u) = (w[0], w[un]);
            let fac = rho * (s - u) / (s - ss);
            let mut qs = [0.0; 5];
            qs[0] = fac;
            for k in 1..4 {
                qs[k] = fac * w[k];
            }
            qs[un] = fac * ss;
            qs[4] = fac * (q[4] / rho + (ss - u) * (ss + w[4] / (rho * (s - u))));
            qs
        };
        let mut f = [0.0; 5];
        if ss >= 0.0 {
            let (fl, qs) = (flux(wl, &ql), star(wl, &ql, sl));
            for k in 0..5 {
                f[k] = fl[k] + sl * (qs[k] - ql[k]);
            }
        } else {
            let (fr, qs) = (flux(wr, &qr), star(wr, &qr, sr));
            for k in 0..5 {
                f[k] = fr[k] + sr * (qs[k] - qr[k]);
            }
        }
        (f, ss)
    }

    /// KEEP-PE two-point flux (Kuya, Totani & Kawai 2018; Shima et al.
    /// 2021): preserves kinetic energy (no spurious production by the
    /// convective terms) and pressure equilibrium across contacts.
    #[inline]
    fn split_flux(&self, a: &V<5>, b: &V<5>, d: usize) -> V<5> {
        let un = 1 + d;
        let rho = 0.5 * (a[0] + b[0]);
        let u = 0.5 * (a[un] + b[un]);
        let c = rho * u;
        let mut f = [0.0; 5];
        f[0] = c;
        for k in 1..4 {
            f[k] = c * 0.5 * (a[k] + b[k]);
        }
        f[un] += 0.5 * (a[4] + b[4]);
        let kin = c * 0.5 * (a[1] * b[1] + a[2] * b[2] + a[3] * b[3]);
        let rhoe = 0.5 * (self.eos.rho_e(a[0], a[4]) + self.eos.rho_e(b[0], b[4]));
        let pw = 0.5 * (a[un] * b[4] + b[un] * a[4]);
        f[4] = kin + u * rhoe + pw;
        f
    }

    #[inline]
    fn flux(&self, w: &V<5>, d: usize) -> V<5> {
        let un = 1 + d;
        let q = self.to_cons(w);
        let u = w[un];
        let mut f = [0.0; 5];
        for k in 0..4 {
            f[k] = q[k] * u;
        }
        f[un] += w[4];
        f[4] = (q[4] + w[4]) * u;
        f
    }

    fn temperature(&self, w: &V<5>) -> f64 {
        self.eos.temperature(w[0], w[4])
    }

    fn temperature_from_phases(&self, _w: &V<5>, phases: &[(f64, f64, f64)]) -> f64 {
        phases[0].1
    }

    fn conduction_properties(&self, w: &V<5>, conductivity: &[f64]) -> (f64, f64, f64) {
        let t = self.temperature(w);
        crate::thermal::properties(t, &[(1.0, t, w[0] * self.eos.cv)], conductivity)
    }

    fn cp(&self, _w: &V<5>) -> f64 {
        self.eos.gamma * self.eos.cv
    }

    fn cv(&self, _w: &V<5>) -> f64 {
        self.eos.cv
    }

    fn phases(&self, w: &V<5>) -> Vec<(f64, f64, f64)> {
        vec![(1.0, self.temperature(w), w[0] * self.eos.cv)]
    }

    fn energy_index(&self) -> Option<usize> {
        Some(4)
    }

    fn eigenvalues(&self, w: &V<5>, d: usize) -> V<5> {
        let (u, c) = (w[1 + d], self.sound_speed(w));
        let mut l = [0.0; 5];
        l[..5].copy_from_slice(&[u - c, u, u, u, u + c]);
        l
    }

    fn dudw(&self, w: &V<5>) -> Mat<5> {
        let (rho, p) = (w[0], w[4]);
        let e = &self.eos;
        let mut j = [[0.0; 5]; 5];
        j[0][0] = 1.0;
        for k in 1..4 {
            j[k][0] = w[k];
            j[k][k] = rho;
            j[4][k] = rho * w[k];
        }
        // ρe = (p + γp∞)(1 − bρ)/(γ − 1) + qρ
        let drhoe_drho = -e.b * (p + e.gamma * e.p_inf) / (e.gamma - 1.0) + e.q;
        let drhoe_dp = (1.0 - e.b * rho) / (e.gamma - 1.0);
        j[4][0] = drhoe_drho + 0.5 * (w[1] * w[1] + w[2] * w[2] + w[3] * w[3]);
        j[4][4] = drhoe_dp;
        j
    }

    #[inline]
    fn admissible(&self, w: &V<5>) -> bool {
        w[0].is_finite() && w[4].is_finite() && self.eos.admissible(w[0], w[4])
    }

    fn sensor(&self, w: &V<5>, out: &mut [f64; 8]) -> usize {
        out[0] = w[0].max(1e-300).ln();
        out[1] = (w[4] + self.eos.p_inf).max(1e-300).ln();
        2
    }
}
