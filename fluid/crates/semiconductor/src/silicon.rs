//! Crystalline silicon: band structure, carrier statistics, mobility and
//! recombination (SI units unless a name says otherwise).
//!
//! Sources:
//! - band gap: Varshni form with Thurmond's (1975) constants,
//!   E_g(T) = 1.170 eV − 4.73e-4 T²/(T + 636 K);
//! - effective densities of states: Green (1990),
//!   N_c = 2.86e19 (T/300)^1.58, N_v = 3.10e19 (T/300)^1.85 cm⁻³;
//! - measured intrinsic density (for comparison and as the device default):
//!   Misiakos & Tsamakis (1993), n_i = 5.29e19 (T/300)^2.54 e^{−6726/T} cm⁻³
//!   (78–340 K);
//! - low-field mobility vs doping: Masetti, Severi & Solmi (1983),
//!   phosphorus for electrons, boron for holes; lattice term
//!   μ_max ∝ (T/300)^−2.5 (electrons), ^−2.2 (holes);
//! - velocity saturation: Canali et al. (1975);
//! - Auger coefficients: Dziewior & Schmid (1977);
//! - SRH lifetimes: Scharfetter doping dependence, τ_max 10 µs (electrons),
//!   3 µs (holes), N_ref 1e16 cm⁻³ (common TCAD defaults);
//! - permittivities: Si 11.7, SiO₂ 3.9.

pub const Q: f64 = 1.602_176_634e-19;
pub const KB: f64 = 1.380_649e-23;
pub const EPS0: f64 = 8.854_187_812_8e-12;
pub const EPS_SI: f64 = 11.7 * EPS0;
pub const EPS_OX: f64 = 3.9 * EPS0;

/// Masetti mobility parameters (cm²/(V s), cm⁻³).
#[derive(Clone, Copy, Debug)]
pub struct Masetti {
    pub mu_min1: f64,
    pub mu_min2: f64,
    pub mu1: f64,
    pub mu_max: f64,
    pub c_r: f64,
    pub c_s: f64,
    pub alpha: f64,
    pub beta: f64,
    pub p_c: f64,
    /// lattice-scattering temperature exponent of μ_max
    pub zeta: f64,
}

pub const MASETTI_N: Masetti = Masetti { mu_min1: 68.5, mu_min2: 68.5, mu1: 56.1, mu_max: 1414.0, c_r: 9.20e16, c_s: 3.41e20, alpha: 0.711, beta: 1.98, p_c: 0.0, zeta: 2.5 };
pub const MASETTI_P: Masetti = Masetti { mu_min1: 44.9, mu_min2: 0.0, mu1: 29.0, mu_max: 470.5, c_r: 2.23e17, c_s: 6.10e20, alpha: 0.719, beta: 2.0, p_c: 9.23e16, zeta: 2.2 };

impl Masetti {
    /// Low-field mobility (m²/(V s)) at total ionized doping `n_tot` (m⁻³).
    pub fn mobility(&self, n_tot: f64, t: f64) -> f64 {
        let n = (n_tot * 1e-6).max(1.0); // cm⁻³
        let mu_max = self.mu_max * (t / 300.0).powf(-self.zeta);
        let mu = self.mu_min1 * (-self.p_c / n).exp() + (mu_max - self.mu_min2) / (1.0 + (n / self.c_r).powf(self.alpha)) - self.mu1 / (1.0 + (self.c_s / n).powf(self.beta));
        mu * 1e-4
    }
}

#[derive(Clone, Debug)]
pub struct Silicon {
    pub t: f64,
    /// intrinsic density used by the device solver (m⁻³)
    pub ni: f64,
    pub tau_n_max: f64,
    pub tau_p_max: f64,
    pub n_ref_tau: f64,
    /// Auger coefficients (m⁶/s)
    pub c_n: f64,
    pub c_p: f64,
    pub mob_n: Masetti,
    pub mob_p: Masetti,
}

impl Silicon {
    /// Silicon at temperature t with the measured intrinsic density.
    pub fn new(t: f64) -> Self {
        Silicon {
            t,
            ni: ni_measured(t),
            tau_n_max: 1e-5,
            tau_p_max: 3e-6,
            n_ref_tau: 1e22,
            c_n: 2.8e-31 * 1e-12,
            c_p: 9.9e-32 * 1e-12,
            mob_n: MASETTI_N,
            mob_p: MASETTI_P,
        }
    }
    /// Thermal voltage kT/q (V).
    pub fn vt(&self) -> f64 {
        KB * self.t / Q
    }
    pub fn tau_n(&self, n_tot: f64) -> f64 {
        self.tau_n_max / (1.0 + n_tot / self.n_ref_tau)
    }
    pub fn tau_p(&self, n_tot: f64) -> f64 {
        self.tau_p_max / (1.0 + n_tot / self.n_ref_tau)
    }
}

/// Band gap (eV), Varshni/Thurmond.
pub fn band_gap(t: f64) -> f64 {
    1.170 - 4.73e-4 * t * t / (t + 636.0)
}

/// Effective densities of states (m⁻³), Green 1990.
pub fn nc(t: f64) -> f64 {
    2.86e25 * (t / 300.0).powf(1.58)
}
pub fn nv(t: f64) -> f64 {
    3.10e25 * (t / 300.0).powf(1.85)
}

/// Intrinsic density from the band structure: √(N_c N_v) e^{−E_g/2kT}.
pub fn ni_from_bands(t: f64) -> f64 {
    (nc(t) * nv(t)).sqrt() * (-band_gap(t) * Q / (2.0 * KB * t)).exp()
}

/// Measured intrinsic density (m⁻³), Misiakos & Tsamakis 1993.
pub fn ni_measured(t: f64) -> f64 {
    5.29e25 * (t / 300.0).powf(2.54) * (-6726.0 / t).exp()
}

/// Saturation velocity (m/s) and Caughey–Thomas exponent, Canali 1975.
pub fn vsat_n(t: f64) -> (f64, f64) {
    (1.07e5 * (300.0 / t).powf(0.87), 1.109 * (t / 300.0).powf(0.66))
}
pub fn vsat_p(t: f64) -> (f64, f64) {
    (8.37e4 * (300.0 / t).powf(0.52), 1.213 * (t / 300.0).powf(0.17))
}

// ---------------------------------------------------------------------
// Heavy doping: Fermi–Dirac integrals and band-gap narrowing

/// Fermi–Dirac integrals F_{1/2}(η) and F_{−1/2}(η) (normalized, → e^η for
/// η → −∞), by Gauss–Legendre panels in t = √x (relative error < 1e-12).
pub fn fermi_half(eta: f64) -> (f64, f64) {
    let tmax = (eta.max(0.0) + 50.0).sqrt();
    let panels = (tmax / 0.25).ceil() as usize;
    let (gx, gw) = crate::analytic::gauss8();
    let (mut f12, mut fm12) = (0.0, 0.0);
    for p in 0..panels {
        let a = tmax * p as f64 / panels as f64;
        let b = tmax * (p + 1) as f64 / panels as f64;
        for (x, w) in gx.iter().zip(gw.iter()) {
            let t = 0.5 * (a + b) + 0.5 * (b - a) * x;
            let e = t * t - eta;
            // 1/(1 + e^e) without overflow
            let occ = if e > 0.0 { (-e).exp() / (1.0 + (-e).exp()) } else { 1.0 / (1.0 + e.exp()) };
            let ww = 0.5 * (b - a) * w;
            f12 += ww * 2.0 * t * t * occ;
            fm12 += ww * 2.0 * occ;
        }
    }
    let s = 1.0 / std::f64::consts::PI.sqrt();
    (2.0 * s * f12, s * fm12)
}

/// Degeneracy correction θ = ln(F_{1/2}(η)/e^η) and dθ/dη. Folded to the
/// Boltzmann limit (θ = 0) below η = −12, where |θ| < 2.2e-6.
pub fn degeneracy(eta: f64) -> (f64, f64) {
    if eta < -12.0 {
        return (0.0, 0.0);
    }
    let (f, fm) = fermi_half(eta);
    (f.ln() - eta, fm / f - 1.0)
}

/// Band-gap narrowing (eV), Slotboom & de Graaff (1976):
/// ΔE_g = 6.92 meV [ln(N/1.3e17) + √(ln²(N/1.3e17) + 0.5)], N in cm⁻³
/// (zero below 1.3e17 by construction of the fit, ~1e-3 eV there).
pub fn bgn_slotboom(n_tot: f64) -> f64 {
    let n = (n_tot * 1e-6).max(1.0);
    let l = (n / 1.3e17).ln();
    6.92e-3 * (l + (l * l + 0.5).sqrt())
}

// ---------------------------------------------------------------------
// Inversion-layer mobility: Lombardi et al. (1988), parameters as
// distributed with common TCAD (electrons / holes)

#[derive(Clone, Copy, Debug)]
pub struct Lombardi {
    /// acoustic phonons: B/E⊥ + C N^λ/(E⊥^{1/3} (T/300)^k); B in cm/s,
    /// C in cm^{5/3} V^{−2/3} s^{−1} (N in cm⁻³)
    pub b: f64,
    pub c: f64,
    pub lambda: f64,
    pub k: f64,
    /// surface roughness δ/E⊥² (V/s)
    pub delta: f64,
}
pub const LOMBARDI_N: Lombardi = Lombardi { b: 4.75e7, c: 580.0, lambda: 0.125, k: 1.0, delta: 5.82e14 };
pub const LOMBARDI_P: Lombardi = Lombardi { b: 9.925e6, c: 2947.0, lambda: 0.0317, k: 1.0, delta: 2.0546e14 };

impl Lombardi {
    /// Combine with the bulk mobility (Matthiessen): field E⊥ in V/m,
    /// doping in m⁻³, mobilities in m²/(V s).
    pub fn combine(&self, mu_bulk: f64, e_perp: f64, n_tot: f64, t: f64) -> f64 {
        let e = (e_perp * 1e-2).max(1.0); // V/cm
        let n = (n_tot * 1e-6).max(1.0);
        let mu_ac = self.b / e + self.c * n.powf(self.lambda) / (e.cbrt() * (t / 300.0).powf(self.k));
        let mu_sr = self.delta / (e * e);
        1.0 / (1.0 / mu_bulk + 1.0 / (mu_ac * 1e-4) + 1.0 / (mu_sr * 1e-4))
    }
}

// ---------------------------------------------------------------------
// Impact ionization: van Overstraeten & de Man (1970), α = a e^{−b/E}

/// Electron ionization coefficient (1/m) at field E (V/m).
pub fn alpha_n(e: f64) -> f64 {
    if e < 5e6 {
        return 0.0; // fold: below 5e4 V/cm α_n < 1e-5 /cm
    }
    7.03e7 * (-1.231e8 / e).exp()
}
/// Hole ionization coefficient (1/m) at field E (V/m).
pub fn alpha_p(e: f64) -> f64 {
    if e < 5e6 {
        return 0.0;
    }
    if e < 4e7 {
        1.582e8 * (-2.036e8 / e).exp()
    } else {
        6.71e7 * (-1.693e8 / e).exp()
    }
}

// ---------------------------------------------------------------------
// Thermal

/// Thermal conductivity of silicon (W/(m K)): 148 (300/T)^1.3 (fit to
/// Glassbrenner & Slack 1964 above 200 K).
pub fn k_si(t: f64) -> f64 {
    148.0 * (300.0 / t).powf(1.3)
}
/// Thermal conductivity of SiO₂ (W/(m K)).
pub const K_OX: f64 = 1.4;



// ---------------------------------------------------------------------
// Carrier screening of impurity scattering (Brooks–Herring)

/// Brooks–Herring screening function F(b) = ln(1 + b) − b/(1 + b),
/// b = 24 m* ε (kT)²/(ħ² q² n'), with n' the screening carrier density
/// (m⁻³) and m* the conductivity mass (electrons 0.26 m₀, holes 0.37 m₀).
pub fn bh_screening(n_screen: f64, t: f64, electrons: bool) -> f64 {
    let m0 = 9.109_383_7e-31;
    let hbar = 1.054_571_817e-34;
    let m = if electrons { 0.26 } else { 0.37 } * m0;
    let kt = KB * t;
    let b = 24.0 * m * EPS_SI * kt * kt / (hbar * hbar * Q * Q * n_screen.max(1e10));
    (1.0 + b).ln() - b / (1.0 + b)
}

impl Masetti {
    /// Mobility with the impurity-scattering part re-screened by the local
    /// carrier density: Masetti is fitted to neutral silicon, where the
    /// screening density equals the doping; the impurity-limited mobility
    /// μ_I = (1/μ_M − 1/μ_L)⁻¹ scales as F(b(N))/F(b(n')) (Brooks–Herring).
    /// In neutral bulk (n' = N) this is Masetti exactly.
    pub fn mobility_screened(&self, n_tot: f64, n_screen: f64, t: f64, electrons: bool) -> f64 {
        let mu_m = self.mobility(n_tot, t);
        let mu_l = self.mu_max * (t / 300.0).powf(-self.zeta) * 1e-4;
        if mu_m >= mu_l || n_tot <= 0.0 {
            return mu_m;
        }
        let mu_i = 1.0 / (1.0 / mu_m - 1.0 / mu_l);
        let ratio = bh_screening(n_tot, t, electrons) / bh_screening(n_screen.max(n_tot), t, electrons);
        1.0 / (1.0 / mu_l + 1.0 / (mu_i * ratio))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn fermi_integrals_limits() {
        // non-degenerate: F_{1/2} → e^η (1 − e^η/2^{3/2})
        for eta in [-20.0f64, -8.0] {
            let (f, fm) = fermi_half(eta);
            let ser = eta.exp() * (1.0 - eta.exp() / 8f64.sqrt() + (2.0 * eta).exp() / 27f64.sqrt());
            assert!((f / ser - 1.0).abs() < 1e-9, "eta {eta}: {f} vs {ser}");
            let serm = eta.exp() * (1.0 - eta.exp() / 2f64.sqrt());
            assert!((fm / serm - 1.0).abs() < 1e-6);
        }
        // degenerate: Sommerfeld expansion
        for eta in [30.0f64, 60.0] {
            let (f, fm) = fermi_half(eta);
            let som = 4.0 / (3.0 * std::f64::consts::PI.sqrt()) * eta.powf(1.5) * (1.0 + std::f64::consts::PI.powi(2) / (8.0 * eta * eta) + 7.0 * std::f64::consts::PI.powi(4) / (640.0 * eta.powi(4)));
            assert!((f / som - 1.0).abs() < 5e-8, "eta {eta}: {f} vs {som}");
            let somm = 2.0 / std::f64::consts::PI.sqrt() * eta.sqrt() * (1.0 - std::f64::consts::PI.powi(2) / (24.0 * eta * eta) - 7.0 * std::f64::consts::PI.powi(4) / (384.0 * eta.powi(4)));
            assert!((fm / somm - 1.0).abs() < 1e-6);
        }
        // F_{-1/2} = dF_{1/2}/dη
        for eta in [-3.0f64, 0.0, 2.5, 7.0] {
            let h = 1e-4;
            let d = (fermi_half(eta + h).0 - fermi_half(eta - h).0) / (2.0 * h);
            assert!((d / fermi_half(eta).1 - 1.0).abs() < 1e-7);
        }
    }
}
