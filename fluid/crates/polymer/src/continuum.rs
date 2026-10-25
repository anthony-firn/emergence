//! Continuum folds: Maxwell, standard linear solid, multi-mode Maxwell
//! (Prony series) and the Oldroyd-B fluid.
//!
//! The Oldroyd-B fluid is the *exact* continuum fold of a dilute solution
//! of Hookean dumbbells: a single normal mode gives one Maxwell element
//! with modulus `G = n k_B T`, relaxation time `λ`, and polymer viscosity
//! `η_p = G λ`; adding the solvent viscosity `η_s` gives Oldroyd-B.
//!
//! Benchmarks implemented here (exact, no fitting):
//! * start-up of steady shear `σ₁₂(t) = η_s γ̇ + η_p γ̇ (1 - e^{-t/λ})`;
//! * first normal-stress difference growth
//!   `N₁(t) = 2 η_p λ γ̇² [1 - (1 + t/λ) e^{-t/λ}]` and the steady value
//!   `N₁ = 2 η_p λ γ̇²`;
//! * the Oldroyd-B identity `N₁ = 2 λ σ₁₂,p² / η_p` at steady state;
//! * oscillatory `G'(ω) = G (ωλ)²/(1+(ωλ)²)`,
//!   `G''(ω) = η_s ω + G ωλ/(1+(ωλ)²)`, with terminal
//!   `G' ∝ ω², G'' ∝ ω`.
//!
//! References: J. G. Oldroyd, *Proc. R. Soc. A* **200**, 523 (1950);
//! R. B. Bird, R. C. Armstrong & O. Hassager, *Dynamics of Polymeric
//! Liquids* Vol. 1 (2nd ed. 1987), ch. 8; R. G. Larson, *Constitutive
//! Equations for Polymer Melts and Solutions* (1988).

/// A single Maxwell mode `G(t) = G e^{-t/λ}`.
#[derive(Clone, Copy, Debug)]
pub struct Maxwell {
    /// plateau / elastic modulus `G`
    pub g: f64,
    /// relaxation time `λ`
    pub lambda: f64,
}

impl Maxwell {
    pub fn g_t(&self, t: f64) -> f64 {
        self.g * (-t / self.lambda).exp()
    }
    /// zero-shear viscosity `η₀ = G λ`
    pub fn eta0(&self) -> f64 {
        self.g * self.lambda
    }
    /// `(G'(ω), G''(ω))`.
    pub fn complex(&self, omega: f64) -> (f64, f64) {
        let x = omega * self.lambda;
        let d = 1.0 + x * x;
        (self.g * x * x / d, self.g * x / d)
    }
}

/// Standard linear solid: a spring `G∞` in parallel with one Maxwell mode
/// `G₁`, `G(t) = G∞ + G₁ e^{-t/λ}`.
#[derive(Clone, Copy, Debug)]
pub struct StandardLinearSolid {
    pub g_inf: f64,
    pub g1: f64,
    pub lambda: f64,
}

impl StandardLinearSolid {
    pub fn g_t(&self, t: f64) -> f64 {
        self.g_inf + self.g1 * (-t / self.lambda).exp()
    }
    pub fn eta0(&self) -> f64 {
        self.g1 * self.lambda
    }
    pub fn complex(&self, omega: f64) -> (f64, f64) {
        let x = omega * self.lambda;
        let d = 1.0 + x * x;
        (self.g_inf + self.g1 * x * x / d, self.g1 * x / d)
    }
}

/// Multi-mode Maxwell (Prony) series, `G(t) = Σ_i g_i e^{-t/λ_i}`.
#[derive(Clone, Debug, Default)]
pub struct PronySeries {
    pub modes: Vec<(f64, f64)>,
}

impl PronySeries {
    pub fn new(modes: Vec<(f64, f64)>) -> Self {
        Self { modes }
    }
    pub fn g_t(&self, t: f64) -> f64 {
        self.modes.iter().map(|&(g, l)| g * (-t / l).exp()).sum()
    }
    pub fn eta0(&self) -> f64 {
        self.modes.iter().map(|&(g, l)| g * l).sum()
    }
    pub fn complex(&self, omega: f64) -> (f64, f64) {
        let (mut gp, mut gpp) = (0.0, 0.0);
        for &(g, l) in &self.modes {
            let x = omega * l;
            let d = 1.0 + x * x;
            gp += g * x * x / d;
            gpp += g * x / d;
        }
        (gp, gpp)
    }
}

/// Oldroyd-B fluid: solvent `η_s` plus one Hookean-dumbbell mode with
/// polymer viscosity `η_p = G λ`.
#[derive(Clone, Copy, Debug)]
pub struct OldroydB {
    pub eta_s: f64,
    pub eta_p: f64,
    pub lambda: f64,
}

impl OldroydB {
    /// A solvent-free upper-convected Maxwell fluid.
    pub fn ucm(g: f64, lambda: f64) -> Self {
        Self { eta_s: 0.0, eta_p: g * lambda, lambda }
    }

    /// `G = n k_B T = η_p / λ`.
    pub fn g_modulus(&self) -> f64 {
        self.eta_p / self.lambda
    }

    /// Polymer stress relaxation modulus `G_p(t) = (η_p/λ) e^{-t/λ}`.
    pub fn g_polymer_t(&self, t: f64) -> f64 {
        self.eta_p / self.lambda * (-t / self.lambda).exp()
    }

    /// Total zero-shear viscosity.
    pub fn eta0(&self) -> f64 {
        self.eta_s + self.eta_p
    }

    /// Start-up of steady shear from rest at `t = 0`.
    pub fn sigma12_startup(&self, t: f64, gdot: f64) -> f64 {
        self.eta_s * gdot + self.eta_p * gdot * (1.0 - (-t / self.lambda).exp())
    }

    /// Steady shear stress.
    pub fn sigma12_steady(&self, gdot: f64) -> f64 {
        self.eta0() * gdot
    }

    /// First normal-stress difference during start-up,
    /// `N₁ = 2 η_p λ γ̇² [1 - (1 + t/λ) e^{-t/λ}]`.
    pub fn n1_startup(&self, t: f64, gdot: f64) -> f64 {
        let x = t / self.lambda;
        2.0 * self.eta_p * self.lambda * gdot * gdot * (1.0 - (1.0 + x) * (-x).exp())
    }

    /// Steady first normal-stress difference `N₁ = 2 η_p λ γ̇²`.
    pub fn n1_steady(&self, gdot: f64) -> f64 {
        2.0 * self.eta_p * self.lambda * gdot * gdot
    }

    /// Dynamic moduli including the solvent: `G'(ω)` and `G''(ω)`.
    pub fn complex(&self, omega: f64) -> (f64, f64) {
        let g = self.g_modulus();
        let (gp, gpp) = Maxwell { g, lambda: self.lambda }.complex(omega);
        (gp, gpp + self.eta_s * omega)
    }

    /// Steady extensional viscosity in uniaxial extension,
    /// `η_E = 3η_s + 3η_p / [(1 - 2Wi)(1 + Wi)]` (Wi = λ ε̇), so the
    /// Trouton ratio is 3 as Wi → 0. It diverges at Wi = ½, the
    /// coil–stretch transition.
    pub fn eta_elongational(&self, ext_rate: f64) -> f64 {
        let wi = self.lambda * ext_rate;
        3.0 * self.eta_s + 3.0 * self.eta_p / ((1.0 - 2.0 * wi) * (1.0 + wi))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn maxwell_terminal_and_crossover() {
        let m = Maxwell { g: 2.0, lambda: 0.5 };
        assert!((m.eta0() - 1.0).abs() < 1e-14);
        // Terminal: G' ∝ ω², G'' ∝ ω, G''/ω → η₀.
        let w = 1e-5;
        let (gp, gpp) = m.complex(w);
        assert!((gpp / w - m.eta0()).abs() / m.eta0() < 1e-8);
        let (gp2, _) = m.complex(2.0 * w);
        assert!((gp2 / gp - 4.0).abs() < 1e-6);
        // At ωλ = 1, G' = G'' = G/2.
        let (gp, gpp) = m.complex(1.0 / m.lambda);
        assert!((gp - m.g / 2.0).abs() < 1e-12);
        assert!((gpp - m.g / 2.0).abs() < 1e-12);
    }

    #[test]
    fn oldroyd_b_startup_normal_stress() {
        let ob = OldroydB { eta_s: 0.1, eta_p: 0.9, lambda: 2.0 };
        let gdot = 0.3;
        // Steady limits.
        assert!((ob.sigma12_steady(gdot) - ob.eta0() * gdot).abs() < 1e-14);
        assert!((ob.n1_steady(gdot) - 2.0 * ob.eta_p * ob.lambda * gdot * gdot).abs() < 1e-14);
        // Short-time N₁ ∝ t².
        let t = 1e-4;
        let exact = ob.n1_startup(t, gdot);
        let series = ob.eta_p * gdot * gdot * t * t / ob.lambda;
        assert!((exact / series - 1.0).abs() < 1e-3, "exact {exact} series {series}");
        // Approach to steady state.
        let big = ob.n1_startup(200.0, gdot);
        assert!((big - ob.n1_steady(gdot)).abs() / ob.n1_steady(gdot) < 1e-6);

        // Oldroyd-B identity N₁ = 2 λ σ₁₂,p² / η_p at steady state.
        let s12p = ob.eta_p * gdot;
        let ident = 2.0 * ob.lambda * s12p * s12p / ob.eta_p;
        assert!((ident - ob.n1_steady(gdot)).abs() < 1e-12);
    }

    #[test]
    fn oldroyd_b_extension_diverges_at_half_weissenberg() {
        let ob = OldroydB { eta_s: 0.0, eta_p: 1.0, lambda: 1.0 };
        // Wi → ½ from below: η_E grows without bound.
        let near = ob.eta_elongational(0.4999);
        let far = ob.eta_elongational(0.499);
        assert!(near > far && near > 1e4, "near {near} far {far}");
        // Wi just above ½ the formula changes sign — the Hookean dumbbell
        // cannot reach a steady stretched state.
        assert!(ob.eta_elongational(0.5001) < 0.0);
    }

    #[test]
    fn standard_linear_solid_has_elastic_plateau() {
        let s = StandardLinearSolid { g_inf: 3.0, g1: 2.0, lambda: 1.0 };
        assert!((s.g_t(0.0) - 5.0).abs() < 1e-14);
        assert!((s.g_t(100.0) - 3.0).abs() < 1e-10);
        // High frequency: G' → G∞ + G₁.
        let (gp, _) = s.complex(1e6);
        assert!((gp - 5.0).abs() < 1e-6);
    }
}
