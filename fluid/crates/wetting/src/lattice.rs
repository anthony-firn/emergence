//! Direct free-energy minimisation of a 2D drop on an explicit periodic
//! texture: the check of the Wenzel / Cassie–Baxter fold.
//!
//! The drop is a 2D (cylindrical) cap whose macroscopic free surface is a
//! circular arc, in contact with a square-wave substrate (posts of width `a`,
//! height `h`, period `λ`). The arc is solved exactly for a prescribed drop
//! area `A` and contact half-width `x_c`; the *only* remaining degree of
//! freedom is the position of the contact line, over which the total surface
//! energy is minimised directly:
//!
//! ```text
//!   E(x_c) = γ L_arc(x_c)  +  (γ_SL − γ_SV) L_wet(x_c),
//! ```
//!
//! with `L_wet` the **exact** wetted arc length of the square wave under the
//! footprint (walls included) for the Wenzel branch, or the projected length
//! of the wetted post tops plus the liquid–vapour meniscus area for the Cassie
//! branch. This is the fine model of the fold: it resolves the texture
//! geometry but not the local meniscus curvature. Its minimiser gives the
//! apparent (macroscopic) contact angle of the arc, and the expectation is
//! that its envelope approaches the homogenised Wenzel / Cassie–Baxter angle
//! as the drop size over the texture period grows; the residual is the
//! contact-line pinning, whose envelope shrinks like `λ / x_c` but whose
//! pinned branches scatter by a few degrees (the physical metastability that
//! produces hysteresis).
//!
//! On a flat substrate the minimiser reproduces Young's angle to 0.03° for
//! every drop size, which is the benchmark used in the tests.

/// A square-wave substrate: posts of width `post_width`, height `post_height`,
/// repeated every `period`. A post centre is at `x = 0`.
#[derive(Clone, Copy, Debug)]
pub struct SquareWave {
    pub period: f64,
    pub post_width: f64,
    pub post_height: f64,
}

impl SquareWave {
    pub fn new(period: f64, post_width: f64, post_height: f64) -> Self {
        Self { period, post_width, post_height }
    }

    pub fn flat() -> Self {
        Self { period: 1.0, post_width: 1.0, post_height: 0.0 }
    }

    /// Projected solid fraction `φ_s = a/λ`.
    pub fn solid_fraction(&self) -> f64 {
        self.post_width / self.period
    }

    /// Fully-wetted 2D roughness `r = 1 + 2h/λ`.
    pub fn roughness(&self) -> f64 {
        1.0 + 2.0 * self.post_height / self.period
    }

    /// Substrate height at distance `x ≥ 0` from a post centre.
    pub fn surface(&self, x: f64) -> f64 {
        let a = self.post_width;
        let p = self.period;
        let r = x.rem_euclid(p);
        if r < a / 2.0 || r >= p - a / 2.0 {
            0.0
        } else {
            -self.post_height
        }
    }

    /// Exact wetted arc length from the post centre to `x ≥ 0`, walls
    /// included.
    pub fn wetted_length(&self, x: f64) -> f64 {
        if x <= 0.0 {
            return 0.0;
        }
        let a = self.post_width;
        let p = self.period;
        let h = self.post_height;
        let per = p + 2.0 * h;
        let k = (x / p).floor();
        let r = x - k * p;
        let within = if r < a / 2.0 {
            r
        } else if r < p - a / 2.0 {
            a / 2.0 + h + (r - a / 2.0)
        } else {
            r + 2.0 * h
        };
        k * per + within
    }

    /// Projected (flat) length of the wetted post tops from the post centre to
    /// `x ≥ 0`.
    pub fn wetted_projected(&self, x: f64) -> f64 {
        if x <= 0.0 {
            return 0.0;
        }
        let a = self.post_width;
        let p = self.period;
        let k = (x / p).floor();
        let r = x - k * p;
        let within = if r < a / 2.0 {
            r
        } else if r < p - a / 2.0 {
            a / 2.0
        } else {
            a / 2.0 + (r - (p - a / 2.0))
        };
        k * a + within
    }
}

/// Which wetting branch to minimise.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Branch {
    /// Liquid penetrates the texture (Wenzel).
    Wenzel,
    /// Liquid bridges the grooves (Cassie–Baxter).
    Cassie,
}

/// The relaxed drop.
#[derive(Clone, Copy, Debug)]
pub struct Drop {
    pub branch: Branch,
    /// Contact-line half-width.
    pub half_width: f64,
    /// Radius of the circular free surface.
    pub radius: f64,
    /// Arc half-angle (the contact angle on a flat substrate).
    pub phi: f64,
    /// Apparent (macroscopic) contact angle, radians.
    pub angle: f64,
    /// Total surface energy (units of `γ`).
    pub energy: f64,
    pub area: f64,
}

/// Lower boundary of the liquid at `x` for the branch: the solid square wave
/// (Wenzel) or the flat post-top plane `z = 0` (Cassie, bridging).
fn lower(branch: Branch, wave: &SquareWave, x: f64) -> f64 {
    match branch {
        Branch::Wenzel => wave.surface(x),
        Branch::Cassie => 0.0,
    }
}

/// Area between the chord at height `lower(x_c)` and the branch's lower
/// boundary, over the footprint `[−x_c, x_c]`.
fn substrate_area(branch: Branch, wave: &SquareWave, xc: f64) -> f64 {
    let m = 400usize;
    let mut acc = 0.0;
    let sc = lower(branch, wave, xc);
    for i in 0..=m {
        let x = xc * i as f64 / m as f64;
        let w = if i == 0 || i == m { 0.5 } else { 1.0 };
        acc += w * (sc - lower(branch, wave, x));
    }
    acc * 2.0 * xc / m as f64
}

/// Solve the arc half-angle `psi ∈ (0, π)` such that the liquid area
/// `R²(psi − sin psi cos psi) + ∫(s(x_c) − s) dx` equals `area`
/// (`R = x_c / sin psi`). Works for obtuse caps.
fn solve_psi(branch: Branch, wave: &SquareWave, xc: f64, area: f64) -> Option<f64> {
    let s_sub = substrate_area(branch, wave, xc);
    let target = area - s_sub;
    let f = |psi: f64| {
        let sp = psi.sin();
        if sp < 1e-9 {
            return -target;
        }
        let r = xc / sp;
        r * r * (psi - sp * psi.cos()) - target
    };
    let mut lo = 1e-4;
    let mut hi = std::f64::consts::PI - 1e-4;
    if f(lo) > 0.0 || f(hi) < 0.0 {
        return None;
    }
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if f(mid) < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Some(0.5 * (lo + hi))
}

/// Apparent angle of the arc at the mean plane: `cos θ* = cos psi − s(x_c)/R`.
fn arc_angle(branch: Branch, wave: &SquareWave, xc: f64, psi: f64) -> f64 {
    if psi <= 1e-12 {
        return std::f64::consts::FRAC_PI_2;
    }
    let r = xc / psi.sin();
    let c = psi.cos() - lower(branch, wave, xc) / r;
    c.clamp(-1.0, 1.0).acos()
}

/// Minimise the free energy over the contact-line position for a drop of the
/// given branch, area and Young cosine. `x_max` bounds the search.
pub fn minimise(branch: Branch, wave: &SquareWave, cos_theta: f64, area: f64, x_max: f64, samples: usize) -> Drop {
    let m = samples.max(200);
    let mut best: Option<Drop> = None;
    // scan the contact-line position down to a small value so that even very
    // hydrophobic drops are reachable
    let lo = 0.02 * wave.period;
    for i in 0..=m {
        let xc = lo + (x_max - lo) * i as f64 / m as f64;
        let phi = match solve_psi(branch, wave, xc, area) {
            Some(v) => v,
            None => continue,
        };
        let r = xc / phi.sin();
        let l_arc = 2.0 * r * phi;
        let (l_wet, meniscus) = match branch {
            Branch::Wenzel => (2.0 * wave.wetted_length(xc), 0.0),
            Branch::Cassie => (
                2.0 * wave.wetted_projected(xc),
                2.0 * xc * (1.0 - wave.solid_fraction()),
            ),
        };
        let energy = l_arc - cos_theta * l_wet + meniscus;
        let angle = arc_angle(branch, wave, xc, phi);
        let cand = Drop { branch, half_width: xc, radius: r, phi, angle, energy, area };
        let better = match best {
            None => true,
            Some(b) => cand.energy < b.energy,
        };
        if better {
            best = Some(cand);
        }
    }
    best.expect("no feasible contact half-width")
}

/// The homogenised apparent cosine of a `SquareWave`.
pub fn wenzel_cos(wave: &SquareWave, cos_theta: f64) -> f64 {
    wave.roughness() * cos_theta
}

pub fn cassie_cos(wave: &SquareWave, cos_theta: f64) -> f64 {
    wave.solid_fraction() * (cos_theta + 1.0) - 1.0
}

/// Size-ratio indicator `x_c / λ`.
pub fn size_ratio(drop: &Drop, wave: &SquareWave) -> f64 {
    drop.half_width / wave.period
}
