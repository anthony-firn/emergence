//! Texture homogenisation: folding a pillar array into an effective contact
//! angle.
//!
//! Two coarse models, both valid for a drop much larger than the texture
//! period:
//!
//! * **Wenzel** (liquid penetrates the texture):
//!   `cos θ* = r cos θ_Y`, with roughness
//!   `r = (true wetted area)/(projected area) ≥ 1` (Wenzel 1936).
//! * **Cassie–Baxter** (liquid rests on the post tops, air in the gaps):
//!   `cos θ* = φ_s (cos θ_Y + 1) − 1`, with solid fraction `φ_s`
//!   (Cassie & Baxter 1944).
//!
//! The Cassie state has the lower surface energy when `cos θ_CB > cos θ_W`,
//! i.e.
//!
//! ```text
//!   cos θ_Y < (φ_s − 1)/(r − φ_s).
//! ```
//!
//! The Cassie state also fails by **impalement** when the liquid pressure
//! exceeds the Laplace pressure a meniscus can hold across a gap of width `w`:
//! `p_imp ≈ −2 γ cos θ_Y / w` (Bico, Thiele & Quéré, *Europhys. Lett.* **61**,
//! 348 (2003); Lafuma & Quéré, *Nature Mater.* **2**, 457 (2003)). A drop of
//! radius `R` has base pressure `≈ 2γ/R`; impalement is expected below
//! `R_c ≈ w/|cos θ_Y|` (the drop/sag criterion).
//!
//! Geometry is exact; the hysteresis estimates at the end are **empirical**
//! scaling models and are labelled as such.

use std::f64::consts::PI;

/// Square array of square posts: centre-to-centre pitch `p`, post width `a`,
/// post height `h`.
#[derive(Clone, Copy, Debug)]
pub struct SquarePosts {
    pub pitch: f64,
    pub post_width: f64,
    pub height: f64,
}

impl SquarePosts {
    pub fn new(pitch: f64, post_width: f64, height: f64) -> Self {
        Self { pitch, post_width, height }
    }

    /// Projected solid fraction `φ_s = a²/p²`.
    pub fn solid_fraction(&self) -> f64 {
        (self.post_width / self.pitch).powi(2)
    }

    /// Fully wetted roughness `r = 1 + 4 a h / p²`.
    pub fn roughness(&self) -> f64 {
        1.0 + 4.0 * self.post_width * self.height / self.pitch.powi(2)
    }

    /// Smallest horizontal gap between posts, `w = p − a`.
    pub fn gap(&self) -> f64 {
        self.pitch - self.post_width
    }

    pub fn wenzel_cos(&self, theta_y: f64) -> f64 {
        self.roughness() * theta_y.cos()
    }

    pub fn cassie_cos(&self, theta_y: f64) -> f64 {
        cassie_cos_general(self.solid_fraction(), theta_y)
    }

    /// Does the Cassie state have the lower energy?
    pub fn cassie_is_stable(&self, theta_y: f64) -> bool {
        self.cassie_cos(theta_y) > self.wenzel_cos(theta_y)
    }

    /// Young contact angle below which the Cassie state is favoured,
    /// `cos θ_Y = (φ_s − 1)/(r − φ_s)`.
    pub fn cassie_transition_cos(&self) -> f64 {
        (self.solid_fraction() - 1.0) / (self.roughness() - self.solid_fraction())
    }

    /// Signed Laplace pressure of a bridging meniscus, `p = −2 γ cos θ_Y / w`
    /// (Bico et al. 2003). For `θ_Y > 90°` this is the positive applied
    /// pressure required to impale the Cassie state; for `θ_Y < 90°` it is
    /// negative, i.e. the grooves fill spontaneously (Wenzel).
    pub fn impalement_pressure(&self, gamma: f64, theta_y: f64) -> f64 {
        -2.0 * gamma * theta_y.cos() / self.gap()
    }

    /// Critical drop radius (m) below which the drop's own Laplace pressure
    /// `2γ/R` reaches the impalement pressure.
    pub fn critical_drop_radius(&self, theta_y: f64) -> f64 {
        let c = theta_y.cos().abs();
        if c < 1e-12 {
            f64::INFINITY
        } else {
            self.gap() / c
        }
    }

    /// Maximum sag of a bridging meniscus across the gap at the impalement
    /// pressure. With curvature radius `R_c = γ/|p_imp|` and pinned contact at
    /// the post edges, the sagitta is `s = R_c − √(R_c² − w²/4)`. The Cassie
    /// state fails when `s ≥ h` (the meniscus reaches the groove floor).
    pub fn sag_before_impalement(&self, gamma: f64, theta_y: f64) -> f64 {
        let p = self.impalement_pressure(gamma, theta_y).abs();
        if p < 1e-12 {
            return f64::INFINITY;
        }
        let rc = gamma / p;
        let half = self.gap() / 2.0;
        rc - (rc * rc - half * half).max(0.0).sqrt()
    }
}

/// Hexagonal (triangular-lattice) array of regular hexagonal posts: lattice
/// constant `pitch` (centre-to-centre), hexagon side `side`, height `h`.
#[derive(Clone, Copy, Debug)]
pub struct HexPosts {
    pub pitch: f64,
    pub side: f64,
    pub height: f64,
}

impl HexPosts {
    pub fn new(pitch: f64, side: f64, height: f64) -> Self {
        Self { pitch, side, height }
    }

    /// Unit cell area `A_cell = (√3/2) p²`.
    pub fn cell_area(&self) -> f64 {
        (3.0f64).sqrt() / 2.0 * self.pitch.powi(2)
    }

    /// `φ_s = (3√3/2) s² / A_cell = 3 s²/p²`.
    pub fn solid_fraction(&self) -> f64 {
        3.0 * self.side.powi(2) / self.pitch.powi(2)
    }

    /// `r = 1 + 6 s h / A_cell = 1 + 4√3 s h / p²`.
    pub fn roughness(&self) -> f64 {
        1.0 + 6.0 * self.side * self.height / self.cell_area()
    }

    /// Gap between opposite sides of neighbouring hexagons.
    pub fn gap(&self) -> f64 {
        self.pitch - (3.0f64).sqrt() * self.side
    }

    pub fn wenzel_cos(&self, theta_y: f64) -> f64 {
        self.roughness() * theta_y.cos()
    }

    pub fn cassie_cos(&self, theta_y: f64) -> f64 {
        cassie_cos_general(self.solid_fraction(), theta_y)
    }

    pub fn cassie_is_stable(&self, theta_y: f64) -> bool {
        self.cassie_cos(theta_y) > self.wenzel_cos(theta_y)
    }

    pub fn cassie_transition_cos(&self) -> f64 {
        (self.solid_fraction() - 1.0) / (self.roughness() - self.solid_fraction())
    }

    pub fn impalement_pressure(&self, gamma: f64, theta_y: f64) -> f64 {
        -2.0 * gamma * theta_y.cos() / self.gap()
    }

    pub fn critical_drop_radius(&self, theta_y: f64) -> f64 {
        let c = theta_y.cos().abs();
        if c < 1e-12 {
            f64::INFINITY
        } else {
            self.gap() / c
        }
    }
}

/// `cos θ_CB = φ_s (cos θ_Y + 1) − 1`.
pub fn cassie_cos_general(phi_s: f64, theta_y: f64) -> f64 {
    phi_s * (theta_y.cos() + 1.0) - 1.0
}

/// Wenzel / Cassie apparent contact angle (radians) for a given effective
/// cosine, clamped to `[0, π]`.
pub fn apparent_angle(cos_theta_eff: f64) -> f64 {
    cos_theta_eff.clamp(-1.0, 1.0).acos()
}

/// **Empirical** hysteresis width from strong dilute defects.
///
/// Joanny & de Gennes (*J. Chem. Phys.* **81**, 552 (1984)) show that a dilute
/// array of strong pinning defects gives a hysteresis proportional to the
/// defect density; Reyssat & Quéré (*J. Phys. Chem. B* **113**, 7996 (2009))
/// confirm the scaling on fakir surfaces and find it controlled by the solid
/// fraction `φ_s`. The dimensionless model used here,
/// `Δcos θ = −C φ_s` (with `C` a supplied pinning strength of order unity),
/// reproduces the scaling only: **the coefficient is not predicted**.
/// The result is `(θ_A − θ_R)` in radians around the equilibrium angle.
pub fn dilute_defect_hysteresis_delta_cos(phi_s: f64, pinning_strength: f64) -> f64 {
    -pinning_strength * phi_s
}

/// Areal density of posts, `φ_s/(post area)` — a convenience for the defect
/// models.
pub fn defect_density(phi_s: f64, pitch: f64) -> f64 {
    phi_s / (pitch * pitch)
}

/// Number of texture periods across a drop of base radius `a`: the size-ratio
/// indicator `a/λ`. The fold is valid when this is large.
pub fn size_ratio(base_radius: f64, period: f64) -> f64 {
    base_radius / period
}

/// Laplace pressure of a drop of radius `R` with surface tension `γ`,
/// `Δp = 2γ/R`.
pub fn drop_laplace_pressure(gamma: f64, radius: f64) -> f64 {
    2.0 * gamma / radius
}

/// Gravity pressure at the base of a puddle of height `h`, `ρ g h`.
pub fn puddle_gravity_pressure(rho: f64, g: f64, height: f64) -> f64 {
    rho * g * height
}

/// Perimeter/area helper: area of a regular hexagon of side `s`.
pub fn hexagon_area(side: f64) -> f64 {
    3.0 * (3.0f64).sqrt() / 2.0 * side * side
}

/// Circumference helper reserved for documentation/tests.
pub fn circle_area(radius: f64) -> f64 {
    PI * radius * radius
}
