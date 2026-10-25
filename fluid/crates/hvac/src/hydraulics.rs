//! Ducts, pipes, fans and pumps.
//!
//! Pressure drop by Darcy–Weisbach, Δp = (f L/D + ΣK) ρV²/2, with the
//! friction factor exact in laminar flow (64/Re) and from the Colebrook
//! equation in turbulent flow, solved to machine precision (ASHRAE
//! Fundamentals 2017 Ch. 3, eq. 31). Between Re 2300 and 4000 the flow is
//! transitional and neither law holds; the turbulent value, the larger, is
//! used and flagged. Rectangular ducts use the circular equivalent for
//! friction, D_e = 1.30 (ab)^0.625 / (a + b)^0.25 (ASHRAE Ch. 21, eq. 25).
//! Fans and pumps: a quadratic through three catalogue points, intersected
//! with the system's own curve.

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Regime {
    Laminar,
    Transitional,
    Turbulent,
}

/// Darcy friction factor for Reynolds number `re` and relative roughness
/// `rr` = ε/D.
pub fn friction(re: f64, rr: f64) -> (f64, Regime) {
    if re < 2300.0 {
        return (64.0 / re, Regime::Laminar);
    }
    // Colebrook: 1/√f = −2 log10(ε/(3.7D) + 2.51/(Re √f)), Newton on x = 1/√f.
    let mut x = 1.0 / swamee_jain(re.max(4000.0), rr).sqrt();
    for _ in 0..50 {
        let a = rr / 3.7 + 2.51 * x / re;
        let g = x + 2.0 * a.log10();
        let dg = 1.0 + 2.0 / std::f64::consts::LN_10 * (2.51 / re) / a;
        let step = g / dg;
        x -= step;
        if step.abs() < 1e-15 * x {
            break;
        }
    }
    (1.0 / (x * x), if re < 4000.0 { Regime::Transitional } else { Regime::Turbulent })
}

/// Swamee & Jain (1976) explicit approximation, used as the starting value
/// (and, in tests, as an independent check: published within ~1 %).
pub fn swamee_jain(re: f64, rr: f64) -> f64 {
    0.25 / (rr / 3.7 + 5.74 / re.powf(0.9)).log10().powi(2)
}

/// Haaland (1983) explicit approximation (published within ~1.5 %).
pub fn haaland(re: f64, rr: f64) -> f64 {
    let x = -1.8 * ((rr / 3.7).powf(1.11) + 6.9 / re).log10();
    1.0 / (x * x)
}

/// Friction-equivalent circular diameter of an `a` × `b` rectangular duct.
pub fn equivalent_diameter(a: f64, b: f64) -> f64 {
    1.30 * (a * b).powf(0.625) / (a + b).powf(0.25)
}

#[derive(Clone, Copy, Debug)]
pub struct Run {
    /// Hydraulic diameter for friction (m) and flow area (m²).
    pub diameter: f64,
    pub area: f64,
    pub length: f64,
    pub roughness: f64,
    /// Sum of fitting loss coefficients (on the run's velocity pressure).
    pub k_sum: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct Drop {
    pub dp: f64,
    pub velocity: f64,
    pub reynolds: f64,
    pub friction: f64,
    pub regime: Regime,
}

impl Run {
    pub fn round(diameter: f64, length: f64, roughness: f64, k_sum: f64) -> Run {
        Run { diameter, area: std::f64::consts::PI * diameter * diameter / 4.0, length, roughness, k_sum }
    }
    pub fn rectangular(a: f64, b: f64, length: f64, roughness: f64, k_sum: f64) -> Run {
        Run { diameter: equivalent_diameter(a, b), area: a * b, length, roughness, k_sum }
    }

    /// Pressure drop at volume flow `q` (m³/s) of a fluid of density `rho`
    /// and dynamic viscosity `mu`.
    pub fn drop(&self, q: f64, rho: f64, mu: f64) -> Drop {
        let v = q / self.area;
        // For the circular equivalent the velocity is that of the same flow
        // in a circular duct of D_e (ASHRAE's convention for friction).
        let v_f = q / (std::f64::consts::PI * self.diameter * self.diameter / 4.0);
        let re = rho * v_f * self.diameter / mu;
        let (f, regime) = if re > 0.0 { friction(re, self.roughness / self.diameter) } else { (0.0, Regime::Laminar) };
        let dp = f * self.length / self.diameter * 0.5 * rho * v_f * v_f + self.k_sum * 0.5 * rho * v * v;
        Drop { dp, velocity: v, reynolds: re, friction: f, regime }
    }
}

/// A fan or pump curve Δp(Q) = c₀ + c₁Q + c₂Q² through three catalogue
/// points, with efficiency at the same points (interpolated the same way).
#[derive(Clone, Copy, Debug)]
pub struct Curve {
    pub c: [f64; 3],
    pub eta: [f64; 3],
}

fn quadratic(p: [(f64, f64); 3]) -> [f64; 3] {
    // Lagrange to monomial.
    let mut c = [0.0; 3];
    for i in 0..3 {
        let (xi, yi) = p[i];
        let others: Vec<f64> = (0..3).filter(|&j| j != i).map(|j| p[j].0).collect();
        let den = (xi - others[0]) * (xi - others[1]);
        let s = yi / den;
        c[0] += s * others[0] * others[1];
        c[1] -= s * (others[0] + others[1]);
        c[2] += s;
    }
    c
}

impl Curve {
    /// From (flow m³/s, pressure Pa, efficiency) at three points.
    pub fn through(points: [(f64, f64, f64); 3]) -> Curve {
        Curve { c: quadratic(points.map(|p| (p.0, p.1))), eta: quadratic(points.map(|p| (p.0, p.2))) }
    }
    pub fn pressure(&self, q: f64) -> f64 {
        self.c[0] + q * (self.c[1] + q * self.c[2])
    }
    pub fn efficiency(&self, q: f64) -> f64 {
        self.eta[0] + q * (self.eta[1] + q * self.eta[2])
    }
}

#[derive(Clone, Copy, Debug)]
pub struct OperatingPoint {
    pub flow: f64,
    pub pressure: f64,
    pub efficiency: f64,
    /// Shaft power Q·Δp/η (W).
    pub power: f64,
}

/// Where the curve meets the system: the flow at which the machine's rise
/// equals the system's drop `system(q)` (monotone increasing in q), found by
/// bisection on [0, q_max]. `None` if they do not cross there.
pub fn operating_point(curve: &Curve, system: &dyn Fn(f64) -> f64, q_max: f64) -> Option<OperatingPoint> {
    let f = |q: f64| curve.pressure(q) - system(q);
    if f(0.0) <= 0.0 || f(q_max) >= 0.0 {
        return None;
    }
    let (mut lo, mut hi) = (0.0, q_max);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if f(mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
        if hi - lo < 1e-15 * q_max {
            break;
        }
    }
    let q = 0.5 * (lo + hi);
    let p = curve.pressure(q);
    let eta = curve.efficiency(q);
    Some(OperatingPoint { flow: q, pressure: p, efficiency: eta, power: q * p / eta })
}
