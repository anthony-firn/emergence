//! The fold to continuum rheology: `mu(I)` for dry granular flow and
//! `mu(J)` / `phi(J)` for the fluid-immersed (viscous-number) case, plus
//! the cheap indicators that bound them.
//!
//! # Dry granular flow: Jop, Forterre & Pouliquen (2006)
//!
//! Local rheology `mu(I) = tau/P` with the inertial number
//!
//! ```text
//! I = gdot d / sqrt(P / rho_p),
//! ```
//!
//! and the closure
//!
//! ```text
//! mu(I) = mu_s + (mu_2 - mu_s) / (I0 / I + 1).
//! ```
//!
//! For glass beads Jop et al. measured `mu_s = tan 20.9°`,
//! `mu_2 = tan 32.8°` and `I0 = 0.279` (Nature **441**, 727, 2006).
//! Regimes: `I <~ 1e-3` quasi-static (Coulomb friction, no rate
//! dependence); `1e-3 <~ I <~ 0.1` dense flow (`mu(I)`); `I >~ 0.1–0.3`
//! collisional (kinetic theory).
//!
//! # Immersed grains: Boyer, Guazzelli & Pouliquen (2011)
//!
//! Replacing inertia by viscous stress gives the viscous number
//! `J = eta_f gdot / P` and the closures (PRL **107**, 188301, 2011):
//!
//! ```text
//! mu(J)  = mu1 + (mu2 - mu1)/(1 + I0/J) + J + (5/2) phi_m sqrt(J),
//! phi(J) = phi_m / (1 + sqrt(J)).
//! ```
//!
//! Frictional particles: `phi_m = 0.575`, `mu1 = 0.32`, `mu2 = 0.7`,
//! `I0 = 0.005`. Frictionless: `phi_m = 0.644`, `mu1 = 0.1`. The constants
//! are those quoted for Boyer et al. (2011) in the review of Guazzelli &
//! Pouliquen (2018, arXiv:1801.03805).
//!
//! [`Jop2006::fit`] least-squares-fits the three parameters to a table of
//! `(I, mu)` points; the DEM's choices of interparticle friction and contact
//! stiffness shift the fitted numbers away from the glass-bead values, which
//! is expected and reported in `analysis/granular.md`.

/// An inertial-number rheology `mu(I)`.
#[derive(Clone, Copy, Debug)]
pub struct Jop2006 {
    /// Quasi-static friction `mu_s`.
    pub mu_s: f64,
    /// Large-`I` friction `mu_2`.
    pub mu_2: f64,
    /// Cross-over inertial number `I0`.
    pub i0: f64,
}

/// Alias for the closure parameters.
pub type MuI = Jop2006;

impl Jop2006 {
    /// The glass-bead values of Jop, Forterre & Pouliquen (2006):
    /// `mu_s = tan 20.9° = 0.3819`, `mu_2 = tan 32.8° = 0.6447`,
    /// `I0 = 0.279`.
    pub fn glass_beads() -> Self {
        Self {
            mu_s: (20.9_f64).to_radians().tan(),
            mu_2: (32.8_f64).to_radians().tan(),
            i0: 0.279,
        }
    }

    /// `mu(I)`.
    pub fn mu(&self, i: f64) -> f64 {
        self.mu_s + (self.mu_2 - self.mu_s) / (self.i0 / i + 1.0)
    }

    /// Inverse: the inertial number at which `mu` is reached (valid for
    /// `mu_s < mu < mu_2`).
    pub fn inertial_number_at(&self, mu: f64) -> f64 {
        assert!(mu > self.mu_s && mu < self.mu_2, "mu outside (mu_s, mu_2)");
        self.i0 / ((self.mu_2 - self.mu_s) / (mu - self.mu_s) - 1.0)
    }
}

/// `mu(I)` with the glass-bead constants, as a free function.
pub fn mu_i_glass_beads(i: f64) -> f64 {
    Jop2006::glass_beads().mu(i)
}

/// Inertial number `I = gdot d / sqrt(P/rho_p)`. `d` is the particle
/// diameter, `rho_p` the particle density and `P` the pressure.
pub fn inertial_number(gdot: f64, d: f64, p: f64, rho_p: f64) -> f64 {
    gdot * d / (p / rho_p).sqrt()
}

/// Regimes of the `I`-rheology.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GranularRegime {
    /// `I < 1e-3`: quasi-static, Coulomb friction.
    QuasiStatic,
    /// `1e-3 <= I < 0.1`: dense flow, `mu(I)`.
    DenseFlow,
    /// `0.1 <= I < 0.3`: transition to collisional.
    Transition,
    /// `I >= 0.3`: collisional / kinetic theory.
    Collisional,
}

impl GranularRegime {
    pub fn label(self) -> &'static str {
        match self {
            GranularRegime::QuasiStatic => "quasi-static (Coulomb)",
            GranularRegime::DenseFlow => "dense flow mu(I)",
            GranularRegime::Transition => "transition",
            GranularRegime::Collisional => "collisional (kinetic theory)",
        }
    }
}

/// Classify by the inertial number.
pub fn granular_regime(i: f64) -> GranularRegime {
    if i < 1.0e-3 {
        GranularRegime::QuasiStatic
    } else if i < 0.1 {
        GranularRegime::DenseFlow
    } else if i < 0.3 {
        GranularRegime::Transition
    } else {
        GranularRegime::Collisional
    }
}

impl Jop2006 {
    /// Least-squares fit of `(mu_s, mu_2, I0)` to points `(I, mu)`. Because
    /// `mu(I) = a + b I/(I + I0)` is linear in `a = mu_s` and
    /// `b = mu_2 - mu_s` for fixed `I0`, the fit scans `I0` and solves the
    /// two-parameter linear problem for each value; returns the parameters
    /// and the RMS residual.
    pub fn fit(data: &[(f64, f64)]) -> (Self, f64) {
        assert!(data.len() >= 3, "need at least 3 points to fit");
        let mut best: Option<(Self, f64)> = None;
        // Log-spaced scan, then two refinement decades.
        let mut lo = 1e-4_f64;
        let mut hi = 1e2_f64;
        for _stage in 0..3 {
            let n = 600;
            for k in 0..=n {
                let t = k as f64 / n as f64;
                let i0 = lo * (hi / lo).powf(t);
                let rows: Vec<Vec<f64>> = data
                    .iter()
                    .map(|&(i, _)| vec![1.0, i / (i + i0)])
                    .collect();
                let rhs: Vec<f64> = data.iter().map(|&(_, m)| m).collect();
                let sol = crate::linalg::least_squares(&rows, &rhs, 2);
                let a = sol[0];
                let b = sol[1];
                // A physical fit must have b > 0 and i0 > 0.
                if b <= 0.0 || a < 0.0 {
                    continue;
                }
                let cand = Jop2006 {
                    mu_s: a,
                    mu_2: a + b,
                    i0,
                };
                let mut sse = 0.0;
                for &(i, m) in data {
                    let r = cand.mu(i) - m;
                    sse += r * r;
                }
                let rmse = (sse / data.len() as f64).sqrt();
                if best.is_none() || rmse < best.as_ref().unwrap().1 {
                    best = Some((cand, rmse));
                }
            }
            if let Some((c, _)) = best {
                // Refine around the best I0.
                lo = (c.i0 / 30.0).max(1e-6);
                hi = c.i0 * 30.0;
            }
        }
        best.unwrap()
    }
}

/// The fluid-immersed (viscous-number) closure of Boyer, Guazzelli &
/// Pouliquen (2011).
#[derive(Clone, Copy, Debug)]
pub struct Boyer2011 {
    /// Critical volume fraction at jamming.
    pub phi_m: f64,
    /// Low-`J` friction.
    pub mu1: f64,
    /// High-`J` friction.
    pub mu2: f64,
    /// Cross-over viscous number.
    pub i0: f64,
}

impl Boyer2011 {
    /// Frictional particles (glass beads): `phi_m = 0.575`, `mu1 = 0.32`,
    /// `mu2 = 0.7`, `I0 = 0.005` (Boyer et al. 2011).
    pub fn frictional() -> Self {
        Self {
            phi_m: 0.575,
            mu1: 0.32,
            mu2: 0.7,
            i0: 0.005,
        }
    }

    /// Frictionless particles: `phi_m = 0.644`, `mu1 = 0.1` (Boyer et al.
    /// 2011).
    pub fn frictionless() -> Self {
        Self {
            phi_m: 0.644,
            mu1: 0.1,
            mu2: 0.7,
            i0: 0.005,
        }
    }

    /// `mu(J) = mu1 + (mu2-mu1)/(1 + I0/J) + J + (5/2) phi_m sqrt(J)`.
    pub fn mu(&self, j: f64) -> f64 {
        self.mu1
            + (self.mu2 - self.mu1) / (1.0 + self.i0 / j)
            + j
            + 2.5 * self.phi_m * j.sqrt()
    }

    /// `phi(J) = phi_m / (1 + sqrt(J))`.
    pub fn phi(&self, j: f64) -> f64 {
        self.phi_m / (1.0 + j.sqrt())
    }
}

/// Viscous number `J = eta_f gdot / P` (Boyer, Guazzelli & Pouliquen 2011).
pub fn viscous_number(eta_f: f64, gdot: f64, p: f64) -> f64 {
    eta_f * gdot / p
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn glass_bead_constants() {
        let j = Jop2006::glass_beads();
        assert!((j.mu_s - 0.381_86).abs() < 1e-4, "{}", j.mu_s);
        assert!((j.mu_2 - 0.644_46).abs() < 1e-4, "{}", j.mu_2);
        assert!((j.i0 - 0.279).abs() < 1e-12);
        // Limits.
        assert!((j.mu(1e-12) - j.mu_s).abs() < 1e-6);
        assert!(j.mu(1e6) < j.mu_2 && (j.mu(1e6) - j.mu_2).abs() < 1e-4);
    }

    #[test]
    fn fit_recovers_its_own_parameters() {
        let truth = Jop2006 {
            mu_s: 0.38,
            mu_2: 0.65,
            i0: 0.3,
        };
        let data: Vec<(f64, f64)> = (0..40)
            .map(|k| {
                let i = 1e-3 * (300.0_f64).powf(k as f64 / 39.0);
                (i, truth.mu(i))
            })
            .collect();
        let (fit, rmse) = Jop2006::fit(&data);
        assert!(rmse < 1e-3, "rmse {rmse}");
        assert!((fit.mu_s - truth.mu_s).abs() < 2e-3, "{fit:?}");
        assert!((fit.mu_2 - truth.mu_2).abs() < 5e-3, "{fit:?}");
        assert!((fit.i0 - truth.i0).abs() / truth.i0 < 0.05, "{fit:?}");
    }

    #[test]
    fn boyer_limits() {
        let b = Boyer2011::frictional();
        // Small J: phi -> phi_m, mu -> mu1.
        assert!((b.phi(1e-8) - b.phi_m).abs() < 1e-3);
        assert!((b.mu(1e-8) - b.mu1).abs() < 1e-3);
        // phi decreases with J.
        assert!(b.phi(0.1) < b.phi(0.001));
    }
}
