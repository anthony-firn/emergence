//! Temperature-dependent thermal equation of state for solid iron phases,
//! after Dorogokupets, Dymshits, Litasov & Sokolova, Sci. Rep. 7, 41863
//! (2017), Table 1 and eqs (1)–(21).
//!
//! The Helmholtz free energy of a solid phase is
//!   F = U₀ + E₀(V) + F_th(V,T) − F_th(V,T₀) + F_e(V,T) − F_e(V,T₀),
//! with T₀ = 298.15 K:
//! * a Vinet cold curve E₀(V), eqs (2), (5);
//! * an Einstein thermal term F_th = 3nRT ln(1 − e^{−Θ/T}) with a
//!   volume-dependent Einstein temperature Θ(V), eq. (6), and the
//!   Altshuler volume dependence of the Grüneisen parameter, eqs (13)–(15);
//! * a free-electron term F_e = −(3/2) n R e₀ (V/V₀)^g T², eq. (16).
//!
//! The magnetic contribution (Inden / Hillert–Jarl, eqs (19)–(21)) is not
//! duplicated here: it is the same heat-capacity function already carried by
//! [`crate::isobaric::Magnetic`], which the phase set adds on top. The two
//! conventions differ only by a term linear in T (the reference of the
//! entropy), which the gauge below absorbs.
//!
//! The published `U₀` fixes the zero of energy but not the phase-boundary
//! temperatures. As in the NASG allotropes fit, a gauge
//!   G += q − q′T + c_off [ (T−T₀) − T ln(T/T₀) ]
//! is carried by the low-temperature phase and fitted so that the two
//! transitions and ΔH(α→γ) come out right; `c_off` is a constant offset to
//! c_v, the direct analogue of the single fitted `cv` of the NASG model.
//! Nothing else is adjusted: volume, bulk modulus, expansion, sound speed
//! and the temperature dependence of c_p are those of the published EOS.
//!
//! Everything in this module is per mole. The [`crate::isobaric::PhaseSet`]
//! divides by the molar mass to return the per-mass quantities the NASG
//! phases use.

/// Molar gas constant, J/(mol K) (CODATA 2018).
pub const R_GAS: f64 = 8.314_462_618;

/// Reference temperature of the Dorogokupets fit, K.
pub const T_REF: f64 = 298.15;

/// One solid phase of the iron thermal EOS.
#[derive(Clone, Copy, Debug)]
pub struct ThermalEos {
    /// Molar volume at the reference isotherm, m³/mol.
    pub v0: f64,
    /// Isothermal bulk modulus at the reference isotherm, Pa.
    pub k0: f64,
    /// Pressure derivative of the bulk modulus at the reference isotherm.
    pub k0p: f64,
    /// Einstein temperature at the reference volume, K.
    pub theta0: f64,
    /// Grüneisen parameter at the reference volume.
    pub gamma0: f64,
    /// Altshuler exponent of the volume dependence of γ.
    pub beta: f64,
    /// Grüneisen parameter at infinite compression.
    pub ginf: f64,
    /// Electronic free-energy coefficient, 1/K.
    pub e0: f64,
    /// Electronic analogue of the Grüneisen parameter.
    pub ge: f64,
    /// Reference energy, J/mol.
    pub u0: f64,
    /// Gauge constant added to G, J/mol.
    pub q: f64,
    /// Gauge slope, J/(mol K): G += q − q′T.
    pub q_prime: f64,
    /// Constant c_v offset, J/(mol K), the analogue of the NASG fitted cv.
    pub cv_off: f64,
    /// Upper temperature of validity at 0.1 MPa, K. The published EOS cannot
    /// be extrapolated far above the melting point at ambient pressure (the
    /// paper says so explicitly; above ~2000 K the isotherm turns over and
    /// K_T < 0), so the phase is omitted from the stability search above it.
    pub t_max: f64,
    /// Molar mass, kg/mol (per-mass conversion for `PhaseSet`).
    pub mm: f64,
}

impl ThermalEos {
    /// η = 3K₀′/2 − 3/2 of eq (2).
    #[inline]
    fn eta(&self) -> f64 {
        1.5 * self.k0p - 1.5
    }

    #[inline]
    fn xx(&self, v: f64) -> f64 {
        (v / self.v0).cbrt()
    }

    /// Cold pressure P₀(V), eq (2).
    pub fn cold_pressure(&self, v: f64) -> f64 {
        let x = self.xx(v);
        let eta = self.eta();
        3.0 * self.k0 * x.powi(-2) * (1.0 - x) * (eta * (1.0 - x)).exp()
    }

    /// Cold energy E₀(V), eq (5).
    pub fn cold_energy(&self, v: f64) -> f64 {
        let x = self.xx(v);
        let eta = self.eta();
        9.0 * self.k0 * self.v0 * eta.powi(-2) * (1.0 - (1.0 - eta * (1.0 - x)) * ((1.0 - x) * eta).exp())
    }

    /// Grüneisen parameter γ(V), eq (13).
    #[inline]
    pub fn gruneisen(&self, v: f64) -> f64 {
        let x = v / self.v0;
        self.ginf + (self.gamma0 - self.ginf) * x.powf(self.beta)
    }

    /// Einstein temperature Θ(V), eq (15).
    #[inline]
    pub fn theta(&self, v: f64) -> f64 {
        let x = v / self.v0;
        self.theta0 * x.powf(-self.ginf) * (((self.gamma0 - self.ginf) / self.beta) * (1.0 - x.powf(self.beta))).exp()
    }

    /// Thermal pressure P_th(V,T), eq (10).
    #[inline]
    pub fn thermal_pressure(&self, v: f64, t: f64) -> f64 {
        let th = self.theta(v);
        3.0 * R_GAS * self.gruneisen(v) * th / (v * ((th / t).exp() - 1.0))
    }

    /// Electronic pressure P_e(V,T), eq (17).
    #[inline]
    pub fn electronic_pressure(&self, v: f64, t: f64) -> f64 {
        let x = v / self.v0;
        (self.ge / v) * 1.5 * R_GAS * self.e0 * x.powf(self.ge) * t * t
    }

    /// Total pressure P(V,T) at the reference isotherm offset, eq (1).
    pub fn pressure(&self, v: f64, t: f64) -> f64 {
        self.cold_pressure(v) + self.thermal_pressure(v, t) - self.thermal_pressure(v, T_REF)
            + self.electronic_pressure(v, t) - self.electronic_pressure(v, T_REF)
    }

    /// Thermal part of the Helmholtz free energy F_th(V,T), eq (6).
    #[inline]
    fn f_th(&self, v: f64, t: f64) -> f64 {
        let th = self.theta(v);
        3.0 * R_GAS * t * (1.0 - (-th / t).exp()).ln()
    }

    /// Thermal internal energy E_th(V,T), eq (8).
    #[inline]
    fn e_th(&self, v: f64, t: f64) -> f64 {
        let th = self.theta(v);
        3.0 * R_GAS * th / ((th / t).exp() - 1.0)
    }

    /// Thermal entropy S_th(V,T), eq (7), J/(mol K).
    #[inline]
    fn s_th(&self, v: f64, t: f64) -> f64 {
        let th = self.theta(v);
        3.0 * R_GAS * (-(1.0 - (-th / t).exp()).ln() + (th / t) / ((th / t).exp() - 1.0))
    }

    #[inline]
    fn xg(&self, v: f64) -> f64 {
        (v / self.v0).powf(self.ge)
    }

    /// Free-electron Helmholtz free energy F_e(V,T), eq (16).
    #[inline]
    fn f_e(&self, v: f64, t: f64) -> f64 {
        -1.5 * R_GAS * self.e0 * self.xg(v) * t * t
    }

    /// Free-electron internal energy E_e(V,T), eq (17).
    #[inline]
    fn e_e(&self, v: f64, t: f64) -> f64 {
        1.5 * R_GAS * self.e0 * self.xg(v) * t * t
    }

    /// Free-electron entropy S_e(V,T), eq (17), J/(mol K).
    #[inline]
    fn s_e(&self, v: f64, t: f64) -> f64 {
        3.0 * R_GAS * self.e0 * self.xg(v) * t
    }

    /// Helmholtz free energy F(V,T), eq (1) without the magnetic part.
    pub fn helmholtz(&self, v: f64, t: f64) -> f64 {
        self.u0 + self.cold_energy(v) + self.f_th(v, t) - self.f_th(v, T_REF) + self.f_e(v, t) - self.f_e(v, T_REF)
    }

    /// Internal energy E(V,T) without the magnetic part.
    pub fn internal_energy(&self, v: f64, t: f64) -> f64 {
        self.u0 + self.cold_energy(v) + self.e_th(v, t) - self.f_th(v, T_REF) + self.e_e(v, t) - self.f_e(v, T_REF)
    }

    /// Entropy S(V,T) without the magnetic part, J/(mol K).
    pub fn entropy(&self, v: f64, t: f64) -> f64 {
        self.s_th(v, t) + self.s_e(v, t)
    }

    /// Gauge contribution to the Gibbs energy, J/mol.
    #[inline]
    fn gauge_g(&self, t: f64) -> f64 {
        self.q - self.q_prime * t + self.cv_off * ((t - T_REF) - t * (t / T_REF).ln())
    }

    /// Gauge contribution to the enthalpy, J/mol.
    #[inline]
    fn gauge_h(&self, t: f64) -> f64 {
        self.q + self.cv_off * (t - T_REF)
    }

    /// Equilibrium molar volume at (p, T), by bisection on P(V,T) = p.
    pub fn volume(&self, p: f64, t: f64) -> f64 {
        let (mut lo, mut hi) = (0.5 * self.v0, 2.0 * self.v0);
        // P decreases with V; if the root is outside the bracket, clamp.
        if self.pressure(lo, t) < p {
            return lo;
        }
        if self.pressure(hi, t) > p {
            return hi;
        }
        for _ in 0..100 {
            let mid = 0.5 * (lo + hi);
            if self.pressure(mid, t) > p {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }

    /// Molar Gibbs energy G(p,T) including the gauge, J/mol.
    pub fn gibbs(&self, p: f64, t: f64) -> f64 {
        let v = self.volume(p, t);
        self.helmholtz(v, t) + p * v + self.gauge_g(t)
    }

    /// Molar enthalpy H(p,T) including the gauge, J/mol.
    pub fn enthalpy(&self, p: f64, t: f64) -> f64 {
        let v = self.volume(p, t);
        self.internal_energy(v, t) + p * v + self.gauge_h(t)
    }

    /// Molar isochoric heat capacity C_v(V,T), J/(mol K), without magnetism.
    pub fn cv(&self, v: f64, t: f64) -> f64 {
        let th = self.theta(v);
        let e = (th / t).exp();
        3.0 * R_GAS * (th / t).powi(2) * e / (e - 1.0).powi(2) + 3.0 * R_GAS * self.e0 * self.xg(v) * t
            + self.cv_off
    }

    /// Molar isothermal bulk modulus K_T(p,T), Pa.
    pub fn k_t(&self, p: f64, t: f64) -> f64 {
        let v = self.volume(p, t);
        let h = v * 1e-6;
        -v * (self.pressure(v + h, t) - self.pressure(v - h, t)) / (2.0 * h)
    }

    /// Volumetric thermal expansion α = (∂P/∂T)_V / K_T, 1/K.
    pub fn alpha(&self, p: f64, t: f64) -> f64 {
        let v = self.volume(p, t);
        let th = self.theta(v);
        let e = (th / t).exp();
        let cv_th = 3.0 * R_GAS * (th / t).powi(2) * e / (e - 1.0).powi(2);
        let cv_e = 3.0 * R_GAS * self.e0 * self.xg(v) * t;
        let dpdt = (self.gruneisen(v) * cv_th + self.ge * cv_e) / v;
        dpdt / self.k_t(p, t)
    }

    /// Molar isobaric heat capacity C_p(p,T), J/(mol K), without magnetism.
    pub fn cp(&self, p: f64, t: f64) -> f64 {
        let v = self.volume(p, t);
        let kt = self.k_t(p, t);
        let alpha = self.alpha(p, t);
        self.cv(v, t) + alpha * alpha * t * v * kt
    }

    /// Adiabatic bulk modulus K_S(p,T), Pa.
    pub fn k_s(&self, p: f64, t: f64) -> f64 {
        let v = self.volume(p, t);
        let kt = self.k_t(p, t);
        let alpha = self.alpha(p, t);
        kt + v * t * (alpha * kt).powi(2) / self.cv(v, t)
    }

    /// Bulk sound speed √(K_S/ρ), m/s.
    pub fn sound_speed(&self, p: f64, t: f64) -> f64 {
        let v = self.volume(p, t);
        (self.k_s(p, t) / (self.mm / v)).sqrt()
    }
}
