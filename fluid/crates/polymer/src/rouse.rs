//! Rouse model: normal modes, stress relaxation and viscosity.
//!
//! A Rouse chain is `N` beads joined by `N-1` entropic springs of stiffness
//! `k = 3 k_B T / b²` (`b` the Kuhn length), each bead with friction `ζ`,
//! with no hydrodynamic interaction. The bead positions obey
//! `ζ dR_i/dt = -k (A R)_i + f_i`, with `A` the Rouse matrix (2 on the
//! diagonal, -1 for neighbours, 1 at the free ends). Its eigenvalues are
//!
//! ```text
//! λ_p = 4 sin²(p π / (2N)),   p = 1..N-1,
//! τ_p = ζ / (k λ_p) = ζ b² / (12 k_B T sin²(p π / (2N))).
//! ```
//!
//! For `p ≪ N` this is the familiar `τ_p ≈ ζ N² b² / (3 π² k_B T p²)`. The
//! stress relaxation modulus of the chain is the sum of the normal-mode
//! contributions (`G(t) = n k_B T Σ_p e^{-2t/τ_p}`, `n` the number of
//! chains per unit volume; the factor 2 is because the stress is bilinear
//! in the Gaussian mode amplitudes, whose autocorrelation decays as
//! `e^{-t/τ_p}`). Integrating gives the Rouse zero-shear viscosity
//! `η₀ = n k_B T Σ_p τ_p / 2 = n ζ N² b² / 36` for large `N`.
//!
//! References: P. E. Rouse, *J. Chem. Phys.* **21**, 1272 (1953);
//! M. Doi & S. F. Edwards, *The Theory of Polymer Dynamics* (1986), ch. 4;
//! M. Rubinstein & R. H. Colby, *Polymer Physics* (2003), ch. 8.

use std::f64::consts::PI;

/// A Rouse chain and its environment.
#[derive(Clone, Copy, Debug)]
pub struct Rouse {
    /// number of beads `N` (so `N-1` springs)
    pub n_beads: usize,
    /// Kuhn length `b`
    pub b: f64,
    /// bead friction coefficient `ζ`
    pub zeta: f64,
    /// thermal energy `k_B T`
    pub kbt: f64,
    /// number density of chains `n`
    pub n_chains: f64,
}

impl Rouse {
    pub fn new(n_beads: usize, b: f64, zeta: f64, kbt: f64, n_chains: f64) -> Self {
        assert!(n_beads >= 2, "a Rouse chain needs at least two beads");
        Self { n_beads, b, zeta, kbt, n_chains }
    }

    /// Number of springs `N-1`.
    #[inline]
    pub fn n_seg(&self) -> usize {
        self.n_beads - 1
    }

    /// Entropic spring constant `k = 3 k_B T / b²`.
    #[inline]
    pub fn spring_k(&self) -> f64 {
        3.0 * self.kbt / (self.b * self.b)
    }

    /// Exact discrete Rouse relaxation time of mode `p` (1-based, `p ≤ N-1`).
    pub fn tau(&self, p: usize) -> f64 {
        assert!(p >= 1 && p < self.n_beads);
        let x = p as f64 * PI / (2.0 * self.n_beads as f64);
        let lam = 4.0 * x.sin() * x.sin();
        self.zeta / (self.spring_k() * lam)
    }

    /// Continuum (large-`N`) approximation `τ_p ≈ ζ N² b² / (3 π² k_B T p²)`
    /// with `N` the number of beads (the standard Rouse expression).
    pub fn tau_continuum(&self, p: usize) -> f64 {
        let n = self.n_beads as f64;
        self.zeta * n * n * self.b * self.b / (3.0 * PI * PI * self.kbt * (p * p) as f64)
    }

    /// Longest Rouse time `τ_1`.
    pub fn tau1(&self) -> f64 {
        self.tau(1)
    }

    /// Mean-square end-to-end distance `⟨R²⟩ = (N-1) b²`.
    pub fn mean_sq_end_to_end(&self) -> f64 {
        self.n_seg() as f64 * self.b * self.b
    }

    /// Stress relaxation modulus `G(t) = n k_B T Σ_p e^{-2t/τ_p}`.
    pub fn g_modulus(&self, t: f64) -> f64 {
        (1..self.n_beads)
            .map(|p| (self.n_chains * self.kbt) * (-2.0 * t / self.tau(p)).exp())
            .sum()
    }

    /// Zero-shear viscosity `η₀ = n k_B T Σ_p τ_p / 2 = ∫G(t)dt`.
    pub fn eta0(&self) -> f64 {
        (1..self.n_beads)
            .map(|p| self.n_chains * self.kbt * self.tau(p) / 2.0)
            .sum()
    }

    /// Dynamic moduli `(G'(ω), G''(ω))`; each Rouse term is a Maxwell mode
    /// with modulus `n k_B T` and time constant `τ_p/2`.
    pub fn complex(&self, omega: f64) -> (f64, f64) {
        let (mut gp, mut gpp) = (0.0, 0.0);
        for p in 1..self.n_beads {
            let x = omega * self.tau(p) / 2.0;
            let d = 1.0 + x * x;
            gp += self.n_chains * self.kbt * x * x / d;
            gpp += self.n_chains * self.kbt * x / d;
        }
        (gp, gpp)
    }

    /// The Rouse spectrum as a Maxwell/Prony series `(g_i, λ_i)`.
    ///
    /// Rouse mode `p` contributes `g_p = n k_B T` with `λ_p = τ_p/2`.
    pub fn prony(&self) -> Vec<(f64, f64)> {
        (1..self.n_beads)
            .map(|p| (self.n_chains * self.kbt, self.tau(p) / 2.0))
            .collect()
    }

    /// Fold the Rouse spectrum onto `m` Maxwell modes: keep the `m` slowest
    /// Rouse modes exactly and lump all the faster ones into a single mode
    /// that preserves the viscosity `η = Σ g λ` and the first moment
    /// `Σ g λ²`. The returned series is ordered slowest first.
    ///
    /// The lumping error in the terminal region is tiny (the slow modes are
    /// exact); the error is concentrated at high frequency. See
    /// `analysis/polymer.md` for numbers.
    pub fn prony_lumped(&self, m: usize) -> Vec<(f64, f64)> {
        let all = self.prony();
        let m = m.clamp(1, all.len());
        let mut kept: Vec<(f64, f64)> = all[..m].to_vec();
        let tail = &all[m..];
        if !tail.is_empty() {
            let eta: f64 = tail.iter().map(|(g, l)| g * l).sum();
            let m1: f64 = tail.iter().map(|(g, l)| g * l * l).sum();
            let lam = m1 / eta;
            kept.push((eta / lam, lam));
        }
        kept
    }
}

/// Dynamic moduli of a Prony (multi-mode Maxwell) series
/// `G(t) = Σ_i g_i e^{-t/λ_i}`.
pub fn prony_complex(series: &[(f64, f64)], omega: f64) -> (f64, f64) {
    let (mut gp, mut gpp) = (0.0, 0.0);
    for &(g, lam) in series {
        let x = omega * lam;
        let d = 1.0 + x * x;
        gp += g * x * x / d;
        gpp += g * x / d;
    }
    (gp, gpp)
}

/// Time-domain `G(t)` of a Prony series.
pub fn prony_g(series: &[(f64, f64)], t: f64) -> f64 {
    series.iter().map(|&(g, lam)| g * (-t / lam).exp()).sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rouse8() -> Rouse {
        // Reduced units: k_B T = 1, b = 1, ζ = 1, one chain per unit volume.
        Rouse::new(8, 1.0, 1.0, 1.0, 1.0)
    }

    #[test]
    fn tau_scaling_and_viscosity_match_rouse() {
        let r = rouse8();
        // τ_1 discrete vs continuum approximation.
        let exact = r.tau(1);
        let approx = r.tau_continuum(1);
        assert!((exact - approx).abs() / approx < 0.05, "exact {exact} approx {approx}");
        // Rouse viscosity η₀ = n ζ N² b² / 36 for large N; finite-N exact sum.
        let n = r.n_beads as f64;
        let rouse_limit = r.n_chains * r.zeta * n * n * r.b * r.b / 36.0;
        let finite = r.eta0();
        // The discrete sum is slightly below the continuum sum.
        assert!((finite - rouse_limit).abs() / rouse_limit < 0.05);

        // Terminal behaviour: G''/ω → η₀, G' ∝ ω².
        let w = 1e-4;
        let (gp, gpp) = r.complex(w);
        assert!((gpp / w - r.eta0()).abs() / r.eta0() < 1e-3);
        assert!(gp / (w * w) > 0.0);
    }

    #[test]
    fn rouse_high_frequency_half_power_law() {
        // Between the slowest and fastest mode, G' and G'' ∝ ω^{1/2}.
        let r = Rouse::new(200, 1.0, 1.0, 1.0, 1.0);
        let w1 = 10.0 / r.tau_continuum(1); // above slowest mode
        let w2 = w1 * 4.0;
        let (_, gpp1) = r.complex(w1);
        let (_, gpp2) = r.complex(w2);
        let slope = (gpp2 / gpp1).ln() / (w2 / w1).ln();
        assert!((slope - 0.5).abs() < 0.05, "Rouse slope {slope}");
    }

    #[test]
    fn lumped_fold_reproduces_terminal_viscosity_exactly() {
        let r = Rouse::new(50, 1.0, 1.0, 1.0, 1.0);
        let full = r.prony();
        let lump = r.prony_lumped(3);
        let eta_full: f64 = full.iter().map(|(g, l)| g * l).sum();
        let eta_lump: f64 = lump.iter().map(|(g, l)| g * l).sum();
        assert!((eta_full - eta_lump).abs() / eta_full < 1e-12);

        // Terminal G'' / ω is exact for any ω → 0.
        let w = 1e-6;
        let (_, gpp_l) = prony_complex(&lump, w);
        let (_, gpp_f) = prony_complex(&full, w);
        assert!((gpp_l - gpp_f).abs() / gpp_f < 1e-9);
    }
}
