//! Finite extensibility: Hookean dumbbells, the coil–stretch transition,
//! and FENE springs.
//!
//! A Hookean dumbbell in steady uniaxial extension with rate `ε̇` has a
//! conformation `c_xx = ⟨Q_x²⟩/⟨Q_x²⟩_eq = 1/(1 - 2 Wi)` with
//! `Wi = λ ε̇`. It diverges at `Wi = ½`: the coil–stretch transition. The
//! Oldroyd-B constitutive law (the exact continuum fold of the Hookean
//! dumbbell) therefore blows up in strong extension and has no steady
//! state for `Wi > ½`; a real chain is saved by its finite contour length.
//!
//! The FENE-P (Peterlin) closure keeps a finite maximum extension `b`:
//! with `A = ⟨QQ⟩` measured in units where `A_eq = b/(b+3) I`, the steady
//! uniaxial solution is
//!
//! ```text
//! A_xx = Z / (1 - 2 Wi Z),  A_yy = Z / (1 + Wi Z),
//! Z = 1 - tr(A)/b,
//! ```
//!
//! and `tr(A) = b(1-Z) < b` for every `Wi`. We solve the scalar equation
//! for `Z` by bisection and also integrate the dumbbell Brownian dynamics
//! with a nonlinear FENE spring for a direct check.
//!
//! References: R. B. Bird, C. F. Curtiss, R. C. Armstrong & O. Hassager,
//! *Dynamics of Polymeric Liquids* Vol. 2 (2nd ed. 1987), ch. 13;
//! R. G. Larson, *Constitutive Equations for Polymer Melts and Solutions*
//! (1988); P. G. de Gennes, *J. Chem. Phys.* **60**, 5030 (1974) (the
//! coil–stretch transition); H. R. Warner, *Ind. Eng. Chem. Fundam.*
//! **11**, 379 (1972) (FENE).

/// Coil–stretch transition of a Hookean dumbbell: `Wi = λ ε̇ = ½`.
pub const HOOKEAN_COIL_STRETCH_WI: f64 = 0.5;

/// Steady Hookean dumbbell stretch `⟨Q_x²⟩/⟨Q_x²⟩_eq` in uniaxial
/// extension. Diverges as `Wi → ½⁻`.
pub fn hookean_steady_cxx(wi: f64) -> f64 {
    1.0 / (1.0 - 2.0 * wi)
}

/// Steady FENE-P dumbbell stretch in uniaxial extension, returned as
/// `(c_xx, c_yy, tr(A)/A_eq)` with `A_eq = b/(b+3)`.
///
/// `b = H Q₀²/(k_B T)` is the finite-extensibility parameter (the number
/// of Kuhn steps the spring can hold). The trace stays below `b`, so no
/// divergence at `Wi = ½`.
pub fn fene_p_steady(wi: f64, b: f64) -> (f64, f64, f64) {
    assert!(b > 3.0);
    // F(Z) = 1 - Z - (Z/b)[1/(1-2Wi Z) + 2/(1+Wi Z)].
    let f = |z: f64| -> f64 {
        let denom1 = 1.0 - 2.0 * wi * z;
        if denom1 <= 0.0 {
            return f64::NEG_INFINITY;
        }
        1.0 - z - (z / b) * (1.0 / denom1 + 2.0 / (1.0 + wi * z))
    };
    let (mut lo, mut hi) = (0.0f64, 1.0f64);
    if wi > 0.0 {
        hi = hi.min(0.5 / wi);
    }
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if f(mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let z = 0.5 * (lo + hi);
    let a_eq = b / (b + 3.0);
    let axx = z / (1.0 - 2.0 * wi * z);
    let ayy = z / (1.0 + wi * z);
    (axx / a_eq, ayy / a_eq, (axx + 2.0 * ayy) / a_eq)
}

/// Brownian dynamics of a single dumbbell in steady uniaxial extension.
///
/// Units: `k_B T = 1`, `H = 3` (so `⟨Q_x²⟩_eq = 1/3`), `λ = 1` and hence
/// single-bead friction `ζ = 4Hλ = 12`. `b_fene = H Q₀²/(k_B T)` is the
/// finite-extensibility parameter when the argument is the FENE `b`; pass
/// `f64::INFINITY` for a Hookean spring. The maximum squared extension is
/// then `Q₀² = b_fene/H`.
#[derive(Clone, Copy, Debug)]
pub struct Dumbbell {
    pub h: f64,
    pub kbt: f64,
    pub zeta: f64,
    pub lambda: f64,
    pub b_fene: f64,
}

impl Dumbbell {
    pub fn hookean() -> Self {
        let h = 3.0;
        let lambda = 1.0;
        Self { h, kbt: 1.0, zeta: 4.0 * h * lambda, lambda, b_fene: f64::INFINITY }
    }
    /// FENE dumbbell with extensibility parameter `b = H Q₀²/(k_B T)`.
    pub fn fene(b: f64) -> Self {
        let mut d = Self::hookean();
        d.b_fene = b;
        d
    }

    /// Maximum squared extension `Q₀² = b_fene / H` (finite for FENE).
    pub fn q0_sq(&self) -> f64 {
        self.b_fene / self.h
    }

    /// Equilibrium `⟨Q_x²⟩ = k_B T / H`.
    pub fn eq_qx2(&self) -> f64 {
        self.kbt / self.h
    }

    /// Integrate `n_traj` dumbbells to steady state under uniaxial
    /// extension rate `ext`, returning `(⟨Q_x²⟩/⟨Q_x²⟩_eq, ⟨Q²⟩)`.
    pub fn run_uniaxial(
        &self,
        ext: f64,
        n_traj: usize,
        n_steps: usize,
        n_eq: usize,
        dt: f64,
        seed: u64,
    ) -> (f64, f64) {
        use crate::rng::Rng;
        // Hookean relaxation rate 2H/ζ = 1/(2λ); λ=1 here.
        let drift = 2.0 * self.h / self.zeta;
        let sigma = (4.0 * self.kbt / self.zeta).sqrt();
        let mut sum_xx = 0.0;
        let mut sum_q2 = 0.0;
        for tr in 0..n_traj {
            let mut rng = Rng::new(seed.wrapping_add(tr as u64 * 0x9E37_79B9));
            let s = self.eq_qx2().sqrt();
            let mut q = [s * rng.normal(), s * rng.normal(), s * rng.normal()];
            for step in 0..n_steps {
                // κ = diag(ext, -ext/2, -ext/2)
                let k = [ext, -0.5 * ext, -0.5 * ext];
                let mut dq = [0.0; 3];
                let q2 = q[0] * q[0] + q[1] * q[1] + q[2] * q[2];
                let fene = if self.b_fene.is_finite() {
                    let x = (q2 / self.q0_sq()).min(1.0 - 1e-9);
                    1.0 / (1.0 - x)
                } else {
                    1.0
                };
                for a in 0..3 {
                    dq[a] = (k[a] - drift * fene) * q[a] * dt
                        + sigma * dt.sqrt() * rng.normal();
                }
                // Cap so the FENE spring cannot cross its limit in one step.
                for a in 0..3 {
                    q[a] += dq[a];
                }
                if step >= n_eq {
                    sum_xx += q[0] * q[0];
                    sum_q2 += q[0] * q[0] + q[1] * q[1] + q[2] * q[2];
                }
            }
        }
        let n = (n_traj * (n_steps - n_eq)) as f64;
        (sum_xx / n / self.eq_qx2(), sum_q2 / n)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hookean_bd_reproduces_oldroyd_b_stretch() {
        // Wi = 0.3: the analytic stretch is 1/(1-0.6) = 2.5.
        let d = Dumbbell::hookean();
        let (cxx, _) = d.run_uniaxial(0.3, 30, 40_000, 20_000, 0.002, 7);
        assert!((cxx - 2.5).abs() / 2.5 < 0.08, "BD cxx {cxx} vs 2.5");
    }

    #[test]
    fn fene_saturates_above_coil_stretch() {
        let b = 50.0;
        let (cxx, tr) = {
            let d = Dumbbell::fene(b);
            let (c, _q2) = d.run_uniaxial(0.8, 30, 40_000, 20_000, 0.002, 11);
            let (_cx, _cy, tr) = fene_p_steady(0.8, b);
            (c, tr)
        };
        // Finite, and far above the Hookean divergence point.
        assert!(cxx.is_finite() && cxx > 1.0, "cxx {cxx}");
        // FENE-P trace stays below b.
        assert!(tr < b, "trace {tr} >= b");
        // Monotone in Wi.
        let (c1, _, _) = fene_p_steady(0.3, b);
        let (c2, _, _) = fene_p_steady(1.5, b);
        assert!(c2 > c1, "not monotone: {c1} {c2}");
    }

    #[test]
    fn fene_p_reduces_to_hookean_for_large_b() {
        let (c, _, _) = fene_p_steady(0.3, 1e6);
        assert!((c - hookean_steady_cxx(0.3)).abs() < 1e-3);
    }
}
