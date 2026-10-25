//! Leaky-dielectric (Taylor–Melcher) drop in a uniform electric field.
//!
//! Rung: quasi-static electrics (∇×E = 0, ∇·J = 0, no free charge in the
//! bulk) plus Stokes flow. The leaky-dielectric model (Taylor 1966; Melcher &
//! Taylor 1969; Saville 1997) gives the electric potential in closed form.
//! With
//!     R = σ_i/σ_o,  S = ε_i/ε_o,
//! and an applied field E₀ along z:
//!     φ_o = (−E₀ r + A/r²) cosθ,   φ_i = B r cosθ,
//!     A = E₀ a³ (R−1)/(R+2),      B = −3 E₀/(R+2).
//! The interface free charge q_s = [[ε E_n]] and the tangential field E_t
//! produce the local electric tractions
//!     f_t^E = q_s E_t = −T sinθ cosθ,      (tangential, along +θ̂)
//!     p_n^E = [[½ ε (E_n² − E_t²)]]        (normal),
//! with
//!     T = 9 ε_o ε₀ E₀² (R−S)/(R+2)².
//! The tangential traction shears the interface into a (Hadamard)
//! recirculation; the normal Maxwell-stress jump enters the normal-stress
//! balance. The two are coupled through the same P₂ Stokes mode: the
//! tangential-driven viscous normal stress is what the normal balance sees
//! together with p_n^E.
//!
//! ## Complete small-deformation solution
//!
//! The electric field is evaluated on the undeformed sphere (its O(D)
//! back-reaction on the traction is O(Ca²) and is dropped). The P₂ Stokes
//! solution driven by the tangential traction f_t^E, with [[u]] = 0 and the
//! steady kinematic condition u_r(a,θ) = 0, is
//!     u_r = U(r) P₂(cosθ),   u_θ = V(r) P₂'(cosθ),
//!     inside : U = ee (r/a⁵ − r³/a⁷),  V = ee (r/(2a⁵) − 5r³/(6a⁷)),
//!     outside: U = ee (r⁻⁴ − a⁻²r⁻²),  V = −ee r⁻⁴/3,
//!     ee = −T a⁵ / [5 (μ_o + μ_i)].
//! (The sign of ee makes the physical viscous tangential jump
//! [[σ_rθ]] = −f_t^E = +T sinθ cosθ.) The interfacial velocity is
//!     u_T(θ) ≡ u_θ(a,θ) = −(9/10) (ε_o ε₀ E₀² a)/(μ_o (1+λ)) · (R−S)/(R+2)²
//!                         · sin2θ,   λ = μ_i/μ_o,
//! i.e. equator→pole for R > S (prolate A) and pole→equator for R < S,
//! matching the standard Taylor circulation (Vlahovska 2019, sidebar
//! "Interfacial velocity"; note that reference prints an extra 1/S).
//!
//! The P₂ normal-stress balance 4γ s/a = p_n^E + [[σ_rr^visc]] fixes the
//! *shape amplitude* s, defined by r = a(1 + s P₂(cosθ)):
//!     s = (3/4) Ca_E Φ / (R+2)² ,
//!     Φ = (R² + 1 − 2S) + (3/5) λ_v (R − S),
//!     λ_v = (2 μ_o + 3 μ_i)/(μ_o + μ_i),
//! which is 4/3 of Taylor's deformation parameter
//!     D = (L − B)/(L + B) = (9/16) Ca_E Φ/(R+2)².
//! The factor 3/4 is pure geometry: for r = a(1 + s P₂),
//! L/a = 1 + s, B/a = 1 − s/2, so D = (3/4) s exactly
//! (Vlahovska 2019, "D = 3s/4"; Das & Saintillan 2020, Eq. 5.5).
//! The historical "4/3 discrepancy" between the P₂ stress balance and
//! Taylor's D is therefore *not* missing normal-stress-driven flow: the flow
//! is driven by the tangential traction and the normal Maxwell stress enters
//! the normal balance. Adding a normal-traction-driven mean flow would
//! violate the steady kinematic condition u_r(a) = 0 and would give the
//! wrong viscosity structure. D > 0 is prolate, D < 0 oblate; the conductor
//! limit R → ∞ gives D = (9/16) Ca_E.
//!
//! The finite-difference check [`LeakyDrop::verify`] verifies the interface
//! velocity continuity, the physical tangential traction jump
//! [[σ_rθ]] = −f_t^E and (via [`LeakyDrop::d_from_stress_balance`]) the
//! normal balance.

use crate::pnp::EPS0;

pub struct LeakyDrop {
    pub a: f64,
    pub eps_o: f64,
    pub eps_i: f64,
    pub sigma_o: f64,
    pub sigma_i: f64,
    pub mu_o: f64,
    pub mu_i: f64,
    pub gamma: f64,
    pub e0: f64,
}

impl LeakyDrop {
    pub fn r_ratio(&self) -> f64 {
        self.sigma_i / self.sigma_o
    }
    pub fn s_ratio(&self) -> f64 {
        self.eps_i / self.eps_o
    }
    pub fn ca_e(&self) -> f64 {
        self.eps_o * EPS0 * self.e0 * self.e0 * self.a / self.gamma
    }
    /// Tangential electric-traction amplitude T (the traction is
    /// f_t^E = −T sinθ cosθ on the interface).
    pub fn t_amp(&self) -> f64 {
        let (r, s) = (self.r_ratio(), self.s_ratio());
        9.0 * self.eps_o * EPS0 * self.e0 * self.e0 * (r - s) / ((r + 2.0) * (r + 2.0))
    }
    pub fn lambda_v(&self) -> f64 {
        (2.0 * self.mu_o + 3.0 * self.mu_i) / (self.mu_o + self.mu_i)
    }
    /// Discriminating function in the (9/16) Ca_E normalisation (Taylor 1966).
    pub fn discriminating(&self) -> f64 {
        let (r, s) = (self.r_ratio(), self.s_ratio());
        (r * r + 1.0 - 2.0 * s) + (3.0 / 5.0) * self.lambda_v() * (r - s)
    }
    /// Taylor's deformation parameter D = (L−B)/(L+B) = (9/16) Ca_E Φ/(R+2)².
    pub fn deformation(&self) -> f64 {
        let r = self.r_ratio();
        (9.0 / 16.0) * self.ca_e() * self.discriminating() / ((r + 2.0) * (r + 2.0))
    }
    /// P₂ shape amplitude s in r = a(1 + s P₂), from Taylor's D via the exact
    /// geometric relation D = 3s/4.
    pub fn taylor_shape_amplitude(&self) -> f64 {
        (4.0 / 3.0) * self.deformation()
    }
    pub fn prolate(&self) -> bool {
        self.discriminating() > 0.0
    }

    /// Electric potential (o: r>a, i: r<a).
    pub fn potential(&self, r: f64, theta: f64) -> (f64, f64) {
        let r_ratio = self.r_ratio();
        let a3 = self.a.powi(3);
        let aa = self.e0 * a3 * (r_ratio - 1.0) / (r_ratio + 2.0);
        let bb = -3.0 * self.e0 / (r_ratio + 2.0);
        let c = theta.cos();
        let phi_o = (-self.e0 * r + aa / (r * r)) * c;
        let phi_i = bb * r * c;
        (phi_o, phi_i)
    }

    /// Flow amplitude ee appearing in the P₂ Stokes solution. The sign is
    /// chosen so that the physical viscous tangential jump is
    /// [[σ_rθ]] = −f_t^E = +T sinθ cosθ.
    fn ee(&self) -> f64 {
        -self.t_amp() * self.a.powi(5) / (5.0 * (self.mu_o + self.mu_i))
    }

    /// Flow: returns (u_r, u_θ) at (r, θ).
    pub fn velocity(&self, r: f64, theta: f64) -> (f64, f64) {
        let a = self.a;
        let ee = self.ee();
        let p2 = 0.5 * (3.0 * theta.cos().powi(2) - 1.0);
        let dp2 = -3.0 * theta.sin() * theta.cos();
        if r >= a {
            let u = ee * (r.powi(-4) - a.powi(-2) * r.powi(-2));
            let v = -ee * r.powi(-4) / 3.0;
            (u * p2, v * dp2)
        } else {
            let u = ee * r / a.powi(5) - ee * r.powi(3) / a.powi(7);
            let v = ee * r / (2.0 * a.powi(5)) - 5.0 * ee * r.powi(3) / (6.0 * a.powi(7));
            (u * p2, v * dp2)
        }
    }

    /// Interfacial tangential velocity u_θ(a, θ).
    pub fn interfacial_velocity(&self, theta: f64) -> f64 {
        self.velocity(self.a, theta).1
    }

    /// Amplitude U_T of the interfacial circulation, u_θ(a,θ) = U_T sin2θ.
    /// Equals −(9/10) ε_oε₀E₀²a/(μ_o(1+λ)) · (R−S)/(R+2)².
    pub fn circulation_amplitude(&self) -> f64 {
        self.ee() / (2.0 * self.a.powi(4))
    }

    /// Pressure P₂ coefficient: p = P(r) P₂.
    pub fn pressure_mode(&self, r: f64) -> f64 {
        let a = self.a;
        let ee = self.ee();
        if r >= a {
            -2.0 * self.mu_o * ee * a.powi(-2) * r.powi(-3)
        } else {
            -7.0 * self.mu_i * ee * r.powi(2) / a.powi(7)
        }
    }

    /// P₂ coefficient of the normal Maxwell-stress jump [[σ_rr^E]] on r=a
    /// (this is the electric pressure jump [[p^E]]).
    pub fn electric_normal_p2(&self) -> f64 {
        let (r, _s) = (self.r_ratio(), self.s_ratio());
        let e0 = self.e0;
        let en_o = 3.0 * r * e0 / (r + 2.0);
        let en_i = 3.0 * e0 / (r + 2.0);
        let et = -3.0 * e0 / (r + 2.0);
        let a_c = 0.5 * (self.eps_o * en_o * en_o - self.eps_i * en_i * en_i) * EPS0;
        let a_t = -0.5 * (self.eps_o - self.eps_i) * et * et * EPS0;
        (2.0 / 3.0) * (a_c - a_t)
    }

    /// P₂ coefficient of the viscous normal-stress jump [[σ_rr^visc]] on r=a
    /// produced by the interfacial circulation.
    pub fn viscous_normal_p2(&self) -> f64 {
        let a = self.a;
        let ee = self.ee();
        -(2.0 * self.mu_o + 3.0 * self.mu_i) * ee * a.powi(-5)
    }

    /// P₂ shape amplitude from the exact normal-stress balance
    /// 4γ s/a = [[σ_rr^E]] + [[σ_rr^visc]].
    pub fn shape_amplitude(&self) -> f64 {
        self.a / (4.0 * self.gamma) * (self.electric_normal_p2() + self.viscous_normal_p2())
    }

    /// The tangential-only P₂ normal-stress balance (kept for continuity).
    /// It now returns the *shape amplitude* s and therefore equals
    /// (4/3)·[`deformation`](Self::deformation) exactly; it is not 4/3 of
    /// Taylor's D by error but by the geometry D = 3s/4.
    pub fn d_from_stress_balance(&self) -> f64 {
        self.shape_amplitude()
    }

    /// Finite-difference verification of the analytic solution. Checks
    /// (a) velocity continuity across r=a and (b) the physical tangential
    /// stress jump [[σ_rθ]] = −f_t^E = +T sinθ cosθ, and returns
    /// (vel_jump, tang_mismatch, s_balance, D_taylor). The tangential stress
    /// is computed by central differences of the analytic velocity:
    ///     σ_rθ = μ[ r ∂_r(u_θ/r) + (1/r) ∂_θ u_r ].
    pub fn verify(&self) -> (f64, f64, f64, f64) {
        let a = self.a;
        let h = 1e-3 * a;
        let vscale = (self.t_amp() * a / self.mu_o).abs().max(1e-300);
        // velocity continuity at a few angles
        let mut vjump = 0.0f64;
        for &th in &[0.3, 0.8, 1.2, 1.9, 2.6] {
            let (uri, uti) = self.velocity(a - h, th);
            let (uro, uto) = self.velocity(a + h, th);
            vjump = vjump.max(((uri - uro).abs() + (uti - uto).abs()) / vscale);
        }
        // physical tangential stress jump +T sinθ cosθ, scaled by T
        let tscale = self.t_amp().abs().max(1e-300);
        let mut tjump = 0.0f64;
        let r = a;
        let drut = |rr: f64, th: f64| {
            let f = |x: f64| self.velocity(x, th).1 / x;
            (f(rr + h) - f(rr - h)) / (2.0 * h)
        };
        let dur =
            |rr: f64, th: f64| (self.velocity(rr, th + h).0 - self.velocity(rr, th - h).0) / (2.0 * h);
        for &th in &[0.4, 1.0, 1.5, 2.0, 2.5] {
            let s_o = self.mu_o * (r * drut(r + h, th) + dur(r + h, th) / r);
            let s_i = self.mu_i * (r * drut(r - h, th) + dur(r - h, th) / r);
            // f_t^E = −T sinθ cosθ, so [[σ_rθ]] = −f_t^E = +T sinθ cosθ
            let tau_t = self.t_amp() * th.sin() * th.cos();
            tjump = tjump.max((s_o - s_i - tau_t).abs() / tscale);
        }
        (vjump, tjump, self.d_from_stress_balance(), self.deformation())
    }
}
