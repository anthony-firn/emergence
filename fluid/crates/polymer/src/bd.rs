//! Brownian dynamics of a Rouse chain, and its Green–Kubo stress modulus.
//!
//! Overdamped bead-spring chain with no hydrodynamic interaction. The
//! Ito SDE for bead `i` is
//!
//! ```text
//! ζ dR_i = F_i dt + √(2 k_B T ζ) dW_i,
//! F_i = k (R_{i+1} - 2 R_i + R_{i-1})   (free ends),
//! ```
//!
//! integrated with Euler–Maruyama. The Kramers shear stress is
//! `S_xy = Σ_bonds k r_x r_y`, and the Green–Kubo estimator
//! `G(t)/(n k_B T) = ⟨S_xy(t) S_xy(0)⟩ / (k_B T)²` should equal the
//! analytic Rouse sum `Σ_p e^{-2t/τ_p}`. The end-to-end vector
//! `R = R_{N-1} - R_0` decays on the longest Rouse time `τ_1`.
//!
//! At equilibrium `⟨S_xy²⟩ = (N-1)(k_B T)²` (one `k_B T` per spring), which
//! fixes the normalisation of `G(0)`.

use crate::rng::Rng;
use crate::rouse::Rouse;

/// A Rouse chain in a Brownian-dynamics box (positions in length units).
pub struct RouseChain {
    pub n_beads: usize,
    pub k: f64,
    pub zeta: f64,
    pub kbt: f64,
    pub pos: Vec<[f64; 3]>,
}

impl RouseChain {
    /// Chain with `n_beads` beads, spring constant `k`, friction `ζ` and
    /// thermal energy `k_B T`, started at Gaussian equilibrium positions.
    pub fn new(n_beads: usize, k: f64, zeta: f64, kbt: f64, rng: &mut Rng) -> Self {
        // Equilibrium bond variance per component: k_B T / k.
        let s = (kbt / k).sqrt();
        let mut pos = vec![[0.0; 3]];
        for _ in 1..n_beads {
            let prev = pos[pos.len() - 1];
            pos.push([
                prev[0] + s * rng.normal(),
                prev[1] + s * rng.normal(),
                prev[2] + s * rng.normal(),
            ]);
        }
        Self { n_beads, k, zeta, kbt, pos }
    }

    /// One Euler–Maruyama step of size `dt`.
    pub fn step(&mut self, dt: f64, rng: &mut Rng) {
        let n = self.n_beads;
        let mut f = vec![[0.0; 3]; n];
        for i in 0..n {
            if i > 0 {
                for a in 0..3 {
                    f[i][a] += self.k * (self.pos[i - 1][a] - self.pos[i][a]);
                }
            }
            if i + 1 < n {
                for a in 0..3 {
                    f[i][a] += self.k * (self.pos[i + 1][a] - self.pos[i][a]);
                }
            }
        }
        let amp = (2.0 * self.kbt * dt / self.zeta).sqrt();
        for i in 0..n {
            for a in 0..3 {
                self.pos[i][a] += dt / self.zeta * f[i][a] + amp * rng.normal();
            }
        }
    }

    /// Kramers shear stress `S_xy = Σ k r_x r_y` (energy units).
    pub fn stress_xy(&self) -> f64 {
        let mut s = 0.0;
        for i in 0..self.n_beads - 1 {
            s += self.k * (self.pos[i + 1][0] - self.pos[i][0]) * (self.pos[i + 1][1] - self.pos[i][1]);
        }
        s
    }

    /// End-to-end vector.
    pub fn end_to_end(&self) -> [f64; 3] {
        let a = self.pos[self.n_beads - 1];
        let b = self.pos[0];
        [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
    }
}

/// Result of a Brownian-dynamics Green–Kubo run.
pub struct BdResult {
    /// lag times
    pub times: Vec<f64>,
    /// normalised stress modulus `G(t)/(n k_B T)`
    pub g: Vec<f64>,
    /// normalised end-to-end autocorrelation `⟨R(t)·R(0)⟩/⟨R²⟩`
    pub ee: Vec<f64>,
}

/// Run `n_traj` independent Rouse chains for `n_steps` steps of size `dt`,
/// sampling every `sample_every` steps, and return the time-averaged stress
/// and end-to-end autocorrelations up to `max_lag` samples.
pub fn green_kubo(
    rouse: &Rouse,
    n_traj: usize,
    n_steps: usize,
    n_eq: usize,
    dt: f64,
    sample_every: usize,
    max_lag: usize,
    seed: u64,
) -> BdResult {
    let k = rouse.spring_k();
    let n_beads = rouse.n_beads;
    let n_samples = (n_steps - n_eq) / sample_every;
    let mut stress = vec![vec![0.0; n_samples]; n_traj];
    let mut ee = vec![vec![[0.0; 3]; n_samples]; n_traj];

    for (tr, _) in (0..n_traj).enumerate() {
        let mut rng = Rng::new(seed.wrapping_add(tr as u64 * 0x9E37_79B9));
        let mut chain = RouseChain::new(n_beads, k, rouse.zeta, rouse.kbt, &mut rng);
        for _ in 0..n_eq {
            chain.step(dt, &mut rng);
        }
        for s in 0..n_samples {
            for _ in 0..sample_every {
                chain.step(dt, &mut rng);
            }
            stress[tr][s] = chain.stress_xy();
            ee[tr][s] = chain.end_to_end();
        }
    }

    // Autocorrelation, averaged over trajectories and time origins. The
    // true equilibrium mean of the off-diagonal stress (and of R) is zero,
    // so no sample mean is subtracted: subtracting the per-trajectory mean
    // of a slowly decaying bilinear stress biases the long-time tail
    // strongly negative, while the residual mean bias without subtraction
    // is only O(τ/T) and is small for the long runs used here.
    let mut g = vec![0.0; max_lag + 1];
    let mut ee_dot = vec![0.0; max_lag + 1];
    let mut count = 0usize;
    for tr in 0..n_traj {
        let m = stress[tr].len();
        if m <= max_lag {
            continue;
        }
        for lag in 0..=max_lag {
            let mut acc = 0.0;
            let mut acc_ee = 0.0;
            for i in 0..m - lag {
                acc += stress[tr][i] * stress[tr][i + lag];
                for a in 0..3 {
                    acc_ee += ee[tr][i][a] * ee[tr][i + lag][a];
                }
            }
            g[lag] += acc / (m - lag) as f64;
            ee_dot[lag] += acc_ee / (m - lag) as f64;
        }
        count += 1;
    }
    let inv = 1.0 / count as f64;
    let kbt2 = rouse.kbt * rouse.kbt;
    let ee0 = rouse.mean_sq_end_to_end();
    for x in g.iter_mut() {
        *x *= inv / kbt2;
    }
    for x in ee_dot.iter_mut() {
        *x *= inv / ee0;
    }
    let times = (0..=max_lag)
        .map(|l| (l * sample_every) as f64 * dt)
        .collect();
    BdResult { times, g, ee: ee_dot }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rouse::Rouse;

    #[test]
    fn green_kubo_matches_analytic_rouse() {
        // N = 8 beads, reduced units. dt resolves the fastest mode
        // (τ_7 ≈ 0.087). Each trajectory is ~80 τ_1 long so that removing
        // the per-trajectory mean does not bias the slowest mode.
        let r = Rouse::new(8, 1.0, 1.0, 1.0, 1.0);
        let res = green_kubo(&r, 200, 40_000, 4_000, 0.005, 10, 120, 12345);
        let analytic = |t: f64| r.g_modulus(t) / r.n_chains / r.kbt;

        // G(0) = N-1 exactly in expectation.
        assert!((res.g[0] - 7.0).abs() / 7.0 < 0.1, "G(0)={}", res.g[0]);

        // Compare the curve where it is well above both the noise and the
        // finite-window mean bias O(Στ_p/T) ≈ 0.02 that affects the tail.
        let mut max_rel = 0.0_f64;
        for (i, &t) in res.times.iter().enumerate() {
            let a = analytic(t);
            if a > 0.3 {
                max_rel = max_rel.max((res.g[i] - a).abs() / a);
            }
        }
        assert!(max_rel < 0.1, "max relative error {max_rel}");

        // End-to-end autocorrelation: starts at ⟨R²⟩, decays on τ_1.
        assert!((res.ee[0] - 1.0).abs() < 0.1, "ee(0)={}", res.ee[0]);
        let tau1 = r.tau1();
        // Long-time log slope ≈ 1/τ_1.
        let i1 = res.times.iter().position(|&t| t >= tau1).unwrap();
        let i2 = res.times.iter().position(|&t| t >= 1.6 * tau1).unwrap();
        if res.ee[i2] > 1e-3 {
            let slope = (res.ee[i1].ln() - res.ee[i2].ln()) / (res.times[i2] - res.times[i1]);
            assert!((slope - 1.0 / tau1).abs() < 0.35 / tau1, "ee slope {slope} vs {}", 1.0 / tau1);
        }
    }
}




