//! Sphere tracing the implicit kernel.
//!
//! The renderer draws `cad::Solid` directly, with no tessellation: a ray is
//! advanced by the field's own value until it reaches the surface. That is
//! only sound if the field never overstates the clearance, and the CAD
//! kernel guarantees exactly that.
//!
//! > **Why this is correct.** Every node of a [`Solid`] keeps `sdf`
//! > 1-Lipschitz and zero on the boundary, so
//! > `|sdf(p)| <= dist(p, boundary)`. A step of `|sdf(p)| ` therefore
//! > cannot reach the boundary, let alone cross it. The march converges on
//! > the *first* crossing along the ray and cannot skip a surface, however
//! > thin — the usual failure of a hand-written distance estimator, which
//! > tunnels through thin walls, cannot happen here.
//!
//! What the guarantee does **not** cover is cost. A ray that grazes a
//! surface tangentially takes arbitrarily many steps, because the field
//! value falls without the ray ever crossing. So the march has a step
//! budget, and running out is reported as [`Trace::Exhausted`] rather than
//! silently as a miss: a renderer that confuses the two produces holes in
//! grazing geometry and no way to tell.
//!
//! Sphere tracing is Hart's (*The Visual Computer* 12(10), 1996); the
//! guarantee above is what he calls a Lipschitz bound on the implicit
//! function, which this kernel provides by construction rather than by
//! assumption.

use cad::solid::Solid;
use geometry::{norm, P3};
use radiation::math::{Bounds, Ray, Vec3};

/// Where a ray met a surface.
#[derive(Clone, Copy, Debug)]
pub struct Hit {
    /// Distance along the ray.
    pub t: f64,
    pub point: P3,
    /// Outward unit normal, from the field gradient.
    pub normal: P3,
    /// Field steps taken, which is the cost of the pixel.
    pub steps: usize,
    /// Whether the ray started inside the solid.
    pub inside: bool,
}

/// The three outcomes of a march, kept distinct on purpose.
#[derive(Clone, Copy, Debug)]
pub enum Trace {
    Hit(Hit),
    /// The ray left the solid's bounding box without reaching the surface.
    Miss,
    /// The step budget ran out while still approaching. The partial state
    /// is carried so a caller can report where, and how close.
    Exhausted { t: f64, distance: f64, steps: usize },
}

impl Trace {
    pub fn hit(self) -> Option<Hit> {
        match self {
            Trace::Hit(h) => Some(h),
            _ => None,
        }
    }

    pub fn is_exhausted(&self) -> bool {
        matches!(self, Trace::Exhausted { .. })
    }
}

/// March settings.
#[derive(Clone, Copy, Debug)]
pub struct MarchOptions {
    /// Surface tolerance: the march stops when `|sdf|` falls below this.
    /// It sets the geometric accuracy of the image.
    pub tolerance: f64,
    pub max_steps: usize,
    pub t_min: f64,
    pub t_max: f64,
    /// Step size for the central-difference normal.
    pub normal_epsilon: f64,
}

impl Default for MarchOptions {
    fn default() -> Self {
        Self {
            tolerance: 1e-6,
            max_steps: 512,
            t_min: 1e-9,
            t_max: f64::INFINITY,
            normal_epsilon: 1e-7,
        }
    }
}

impl MarchOptions {
    /// Settings scaled to a scene of the given size: a tolerance one part
    /// in a million of it, and a normal step ten times smaller again.
    pub fn for_scale(scale: f64) -> Self {
        Self {
            tolerance: 1e-6 * scale,
            max_steps: 512,
            t_min: 1e-9 * scale,
            t_max: f64::INFINITY,
            normal_epsilon: 1e-5 * scale,
        }
    }
}

/// Clip a ray to a solid's bounding box, returning the interval of `t`
/// worth marching. `None` when the ray misses the box entirely.
pub fn clip_to_bounds(solid: &Solid, ray: Ray, opts: &MarchOptions) -> Option<(f64, f64)> {
    let (lo, hi) = solid.bbox();
    if (0..3).any(|k| !lo[k].is_finite() || !hi[k].is_finite()) {
        // Unbounded: march the whole ray.
        return Some((opts.t_min, opts.t_max));
    }
    let bounds = Bounds {
        lo: Vec3(lo),
        hi: Vec3(hi),
    };
    let limit = if opts.t_max.is_finite() { opts.t_max } else { 1e30 };
    let (a, b) = bounds.interval(ray, limit)?;
    let a = a.max(opts.t_min);
    if a > b {
        None
    } else {
        Some((a, b))
    }
}

/// The outward unit normal of the field at `p`, by central differences.
///
/// For an exact solid this is the true surface normal. For a boolean tree
/// the field is only piecewise smooth, and the gradient is still the
/// correct normal of whichever surface is nearest — except exactly on a
/// seam, where the one-sided normals differ and the difference quotient
/// returns their average. That is a measure-zero set in the image and
/// shows as a single shaded pixel line along a crease.
pub fn normal_at(solid: &Solid, p: P3, epsilon: f64) -> P3 {
    let g = solid.gradient(p, epsilon);
    let n = norm(g);
    if n > 0.0 {
        [g[0] / n, g[1] / n, g[2] / n]
    } else {
        [0.0, 0.0, 1.0]
    }
}

/// Sphere-trace `ray` against `solid`.
pub fn march(solid: &Solid, ray: Ray, opts: &MarchOptions) -> Trace {
    let Some((t_enter, t_exit)) = clip_to_bounds(solid, ray, opts) else {
        return Trace::Miss;
    };
    let origin = ray.at(t_enter);
    let first = solid.sdf(origin.0);
    // Already on or inside the surface at the entry point.
    let inside = first < 0.0;
    let sign = if inside { -1.0 } else { 1.0 };
    if first.abs() <= opts.tolerance {
        return Trace::Hit(Hit {
            t: t_enter,
            point: origin.0,
            normal: normal_at(solid, origin.0, opts.normal_epsilon),
            steps: 1,
            inside,
        });
    }

    let mut t = t_enter;
    let mut steps = 0usize;
    let mut distance = first.abs();
    while steps < opts.max_steps {
        let p = ray.at(t);
        let d = sign * solid.sdf(p.0);
        steps += 1;
        distance = d;
        if d <= opts.tolerance {
            return Trace::Hit(Hit {
                t,
                point: p.0,
                normal: normal_at(solid, p.0, opts.normal_epsilon),
                steps,
                inside,
            });
        }
        // The step the guarantee allows. Nothing larger is safe; nothing
        // smaller is needed.
        t += d;
        if t > t_exit {
            return Trace::Miss;
        }
    }
    Trace::Exhausted { t, distance, steps }
}

/// Is the segment from `from` towards `dir` of length `limit` blocked?
///
/// Shadow rays stop at the first hit, so they are cheaper than a full
/// march: the tolerance can be loose and the budget small.
pub fn occluded(solid: &Solid, from: P3, dir: P3, limit: f64, opts: &MarchOptions) -> bool {
    let ray = Ray::new(Vec3(from), Vec3(dir));
    let mut local = *opts;
    local.t_max = limit;
    match march(solid, ray, &local) {
        Trace::Hit(h) => h.t < limit,
        // An exhausted shadow ray is treated as grazing, which is the
        // conservative choice for a shadow: it does not darken.
        Trace::Miss | Trace::Exhausted { .. } => false,
    }
}
