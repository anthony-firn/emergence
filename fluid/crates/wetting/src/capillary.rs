//! Exact planar capillary results used as benchmarks.
//!
//! * **Jurin rise** in a circular capillary of radius `r`:
//!   `h = 2 γ cos θ / (ρ g r)` (Jurin 1718; de Gennes, Brochard-Wyart &
//!   Quéré, *Capillarity and Wetting Phenomena*, 2004, ch. 2).
//! * **Meniscus at a vertical wall.** With `z` the height above the flat
//!   far-field level and `φ` the tangent angle to the horizontal, the planar
//!   Young–Laplace equation `γ dφ/ds = −ρ g z` has the first integral
//!   `z = 2 l_c sin(φ/2)`. The tangent at a vertical wall is
//!   `φ_w = π/2 − θ`, so the wall height is
//!   `z_w = l_c √(2 (1 − sin θ))` and the profile is
//!   `x(z) = l_c [F(u_w) − F(u)]`, `F(u) = 2√(1−u²) − artanh√(1−u²)`,
//!   `u = z/(2 l_c)`.
//! * **2D (planar) puddle height.** A flat-topped puddle has zero curvature at
//!   the top and the exact first integral `z = 2 l_c sin(φ/2)`, so its height is
//!   `h = 2 l_c sin(θ/2)`.

use crate::constants::capillary_length;

/// Jurin capillary rise (m). `radius` is the tube radius.
pub fn jurin_rise(gamma: f64, theta: f64, rho: f64, g: f64, radius: f64) -> f64 {
    2.0 * gamma * theta.cos() / (rho * g * radius)
}

/// Height of the planar meniscus on a vertical wall, `l_c √(2(1−sin θ))` (m).
pub fn planar_meniscus_wall_height(gamma: f64, theta: f64, rho: f64, g: f64) -> f64 {
    let lc = capillary_length(gamma, rho, g);
    lc * (2.0 * (1.0 - theta.sin())).sqrt()
}

/// Height of a large 2D (planar) puddle, `2 l_c sin(θ/2)` (m).
pub fn planar_puddle_height(gamma: f64, theta: f64, rho: f64, g: f64) -> f64 {
    let lc = capillary_length(gamma, rho, g);
    2.0 * lc * (theta / 2.0).sin()
}

/// Antiderivative of the planar meniscus slope, `F(u) = 2√(1−u²) − artanh√(1−u²)`.
fn meniscus_f(u: f64) -> f64 {
    let s = (1.0 - u * u).max(0.0).sqrt();
    // artanh(s) = 0.5 ln((1+s)/(1-s)); handle s → 1 (u → 0) which sends x → ∞.
    if s >= 1.0 - 1e-15 {
        return f64::NEG_INFINITY;
    }
    2.0 * s - 0.5 * ((1.0 + s) / (1.0 - s)).ln()
}

/// Horizontal distance `x(z)` from a vertical wall to the planar meniscus at
/// height `z` above the flat far field. `theta` is the contact angle through
/// the liquid. Returns `f64::INFINITY` for `z <= 0`.
pub fn planar_meniscus_x(z: f64, gamma: f64, theta: f64, rho: f64, g: f64) -> f64 {
    if z <= 0.0 {
        return f64::INFINITY;
    }
    let lc = capillary_length(gamma, rho, g);
    let zw = planar_meniscus_wall_height(gamma, theta, rho, g);
    let u = z / (2.0 * lc);
    let uw = zw / (2.0 * lc);
    lc * (meniscus_f(uw) - meniscus_f(u))
}

/// One sample of the planar meniscus profile.
#[derive(Clone, Copy, Debug)]
pub struct MeniscusSample {
    pub x: f64,
    pub z: f64,
    pub phi: f64,
}

/// Integrate the planar meniscus ODE away from a vertical wall, for validation
/// against [`planar_meniscus_x`]. With `s` the arclength outward from the wall
/// and `φ` the tangent angle to the horizontal,
/// `dφ/ds = −(2/l_c) sin(φ/2)`, `dz/ds = −sin φ`, `dx/ds = cos φ`,
/// whose first integral is `z = 2 l_c sin(φ/2)`. Returns samples with
/// increasing `x`.
pub fn planar_meniscus_profile(gamma: f64, theta: f64, rho: f64, g: f64, s_max: f64, ds: f64) -> Vec<MeniscusSample> {
    let lc = capillary_length(gamma, rho, g);
    let mut phi = std::f64::consts::FRAC_PI_2 - theta;
    let mut z = 2.0 * lc * (phi / 2.0).sin();
    let mut x = 0.0;
    let mut s = 0.0;
    let mut out = vec![MeniscusSample { x, z, phi }];
    let dphi = |p: f64| -(2.0 / lc) * (p / 2.0).sin();
    while s < s_max && phi > 1e-6 {
        let k1 = ((phi).cos(), -(phi).sin(), dphi(phi));
        let k2 = ((phi + 0.5 * ds * k1.2).cos(), -(phi + 0.5 * ds * k1.2).sin(), dphi(phi + 0.5 * ds * k1.2));
        let k3 = ((phi + 0.5 * ds * k2.2).cos(), -(phi + 0.5 * ds * k2.2).sin(), dphi(phi + 0.5 * ds * k2.2));
        let k4 = ((phi + ds * k3.2).cos(), -(phi + ds * k3.2).sin(), dphi(phi + ds * k3.2));
        x += ds / 6.0 * (k1.0 + 2.0 * k2.0 + 2.0 * k3.0 + k4.0);
        z += ds / 6.0 * (k1.1 + 2.0 * k2.1 + 2.0 * k3.1 + k4.1);
        phi += ds / 6.0 * (k1.2 + 2.0 * k2.2 + 2.0 * k3.2 + k4.2);
        s += ds;
        out.push(MeniscusSample { x, z, phi });
    }
    out
}

/// A 2D (planar) sessile puddle integrated in `φ` from the flat apex; the
/// exact solution is `z = 2 l_c sin(φ/2)`, returned as `(φ, z)`. The final
/// step is shortened so that the integration lands exactly on `theta`.
pub fn planar_puddle_profile(gamma: f64, theta: f64, rho: f64, g: f64, ds: f64) -> Vec<(f64, f64)> {
    let lc = capillary_length(gamma, rho, g);
    let z = |p: f64| 2.0 * lc * (p / 2.0).sin();
    let mut phi = 1e-6;
    let mut out = vec![(phi, z(phi))];
    let f = |p: f64| 2.0 / lc * (p / 2.0).sin();
    while phi < theta {
        // cap the step so the integration ends on `theta`
        let remaining = theta - phi;
        let h = ds.min(remaining / f(phi).max(1e-12)).max(1e-12);
        let k1 = f(phi);
        let k2 = f(phi + 0.5 * h * k1);
        let k3 = f(phi + 0.5 * h * k2);
        let k4 = f(phi + h * k3);
        phi += h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        out.push((phi.min(theta), z(phi.min(theta))));
    }
    out
}
