//! Implicit solids: the modelling kernel.
//!
//! A [`Solid`] is a tree of primitives, placements, booleans and sweeps
//! evaluated as a scalar field. Every node in the tree keeps one invariant:
//!
//! > **Kernel invariant.** `sdf` is 1-Lipschitz on R^3, is negative exactly
//! > on the interior of the solid it denotes, positive exactly outside, and
//! > zero exactly on its boundary.
//!
//! The invariant is what makes the kernel usable without a B-rep. A
//! 1-Lipschitz function that vanishes on the boundary satisfies
//! `|sdf(p)| <= dist(p, boundary)`, so
//!
//! * the **sign is always exact** — membership, meshing and the
//!   inside/outside test a solver needs are never approximate;
//! * the **magnitude is a safe lower bound** on the true clearance, so
//!   `sdf` can be used to march rays, to size a grid, or to offset a
//!   surface, and it can only ever be conservative;
//! * [`Solid::exact`] reports when the bound is tight, i.e. when `sdf` *is*
//!   the signed distance. Primitives and sweeps are exact; booleans and
//!   arrays keep the sign but may fall short of the true distance near a
//!   concave seam.
//!
//! This is the standard function-representation algebra (Pasko, Adzhiev,
//! Sourin & Savchenko, *The Visual Computer* 11(8), 1995). The exact
//! primitive distances for the capped cylinder and cone follow Quilez's
//! derivations (iquilezles.org/articles/distfunctions); the smooth union is
//! the quadratic polynomial blend from the same source. The triply periodic
//! level sets are the standard nodal approximations (Schoen, *NASA TN
//! D-5541*, 1970; von Schnering & Nesper, *Z. Phys. B* 83, 1991), each
//! divided here by an analytic bound on its gradient so that the kernel
//! invariant holds for them too.

use std::sync::Arc;

use geometry::{add, cross, dot, norm, normalize, scale, sub, P3};

use crate::sketch::{Region2, V2};

/// Rotation (or reflection), uniform scale and translation.
///
/// `rot` must be orthogonal; `scale` must be positive. A reflection
/// (`det = -1`) is allowed and is how mirroring is done.
#[derive(Clone, Copy, Debug)]
pub struct Placement {
    pub rot: [[f64; 3]; 3],
    pub translate: P3,
    pub scale: f64,
}

impl Default for Placement {
    fn default() -> Self {
        Self::IDENTITY
    }
}

impl Placement {
    pub const IDENTITY: Placement = Placement {
        rot: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        translate: [0.0; 3],
        scale: 1.0,
    };

    pub fn translation(t: P3) -> Self {
        Placement { translate: t, ..Self::IDENTITY }
    }

    pub fn scaling(s: f64) -> Self {
        Placement { scale: s, ..Self::IDENTITY }
    }

    /// Right-handed rotation of `angle` radians about `axis`
    /// (Rodrigues' formula).
    pub fn rotation(axis: P3, angle: f64) -> Self {
        let a = normalize(axis);
        let (s, c) = angle.sin_cos();
        let t = 1.0 - c;
        Placement {
            rot: [
                [c + a[0] * a[0] * t, a[0] * a[1] * t - a[2] * s, a[0] * a[2] * t + a[1] * s],
                [a[1] * a[0] * t + a[2] * s, c + a[1] * a[1] * t, a[1] * a[2] * t - a[0] * s],
                [a[2] * a[0] * t - a[1] * s, a[2] * a[1] * t + a[0] * s, c + a[2] * a[2] * t],
            ],
            ..Self::IDENTITY
        }
    }

    /// Reflection in the plane through the origin with unit normal `n`.
    pub fn mirror(n: P3) -> Self {
        let n = normalize(n);
        let mut rot = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                rot[i][j] = if i == j { 1.0 } else { 0.0 } - 2.0 * n[i] * n[j];
            }
        }
        Placement { rot, ..Self::IDENTITY }
    }

    /// The rotation taking `+z` to `dir`.
    pub fn align_z(dir: P3) -> Self {
        let d = normalize(dir);
        let axis = cross([0.0, 0.0, 1.0], d);
        let s = norm(axis);
        if s < 1e-15 {
            if d[2] > 0.0 {
                Self::IDENTITY
            } else {
                Self::rotation([1.0, 0.0, 0.0], std::f64::consts::PI)
            }
        } else {
            Self::rotation(axis, d[2].clamp(-1.0, 1.0).acos())
        }
    }

    pub fn with_translation(mut self, t: P3) -> Self {
        self.translate = t;
        self
    }

    pub fn apply(&self, p: P3) -> P3 {
        let mut q = [0.0; 3];
        for i in 0..3 {
            q[i] = self.scale * (self.rot[i][0] * p[0] + self.rot[i][1] * p[1] + self.rot[i][2] * p[2]);
        }
        add(q, self.translate)
    }

    /// The inverse map, used to pull a query point into local coordinates.
    pub fn inverse_apply(&self, p: P3) -> P3 {
        let d = sub(p, self.translate);
        let mut q = [0.0; 3];
        for i in 0..3 {
            q[i] = (self.rot[0][i] * d[0] + self.rot[1][i] * d[1] + self.rot[2][i] * d[2]) / self.scale;
        }
        q
    }
}

/// A field supplied by a generator that has no closed-form node.
pub struct CustomField {
    pub name: String,
    #[allow(clippy::type_complexity)]
    pub f: Box<dyn Fn(P3) -> f64 + Send + Sync>,
    pub bbox: (P3, P3),
    /// Whether `f` is the exact signed distance (it must always be
    /// 1-Lipschitz with the right sign).
    pub exact: bool,
}

impl std::fmt::Debug for CustomField {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "CustomField({})", self.name)
    }
}

/// A triply periodic level set.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Tpms {
    /// `sin x cos y + sin y cos z + sin z cos x` (Schoen's gyroid).
    Gyroid,
    /// `cos x + cos y + cos z` (Schwarz primitive).
    SchwarzP,
    /// Schwarz diamond.
    Diamond,
}

/// An implicit solid.
#[derive(Clone, Debug)]
pub enum Solid {
    Empty,
    /// All of space.
    Full,
    /// `{ p : dot(unit(normal), p) <= offset / |normal| }`.
    HalfSpace { normal: P3, offset: f64 },
    Sphere { center: P3, radius: f64 },
    /// Box of half-extents `half` with edges rounded by `round`.
    Cuboid { center: P3, half: P3, round: f64 },
    /// Capped cylinder between `a` and `b`.
    Cylinder { a: P3, b: P3, radius: f64 },
    /// Capped cone (frustum) between `a` and `b`.
    Cone { a: P3, b: P3, radius_a: f64, radius_b: f64 },
    /// Points within `radius` of the segment `a`-`b`.
    Capsule { a: P3, b: P3, radius: f64 },
    /// Torus about `+z` through `center`.
    Torus { center: P3, major: f64, minor: f64 },
    /// A sketch given thickness along `+z`.
    Extrude { region: Region2, z0: f64, z1: f64 },
    /// A sketch in the `(r, z)` meridian plane turned about the `z` axis.
    ///
    /// Exact, because for a full revolution the distance in the meridian
    /// plane *is* the distance in space: the nearest boundary point of the
    /// solid is always at the same azimuth.
    ///
    /// **A solid that contains the axis needs a profile that crosses it.**
    /// The field is evaluated at `r = sqrt(x^2 + y^2) >= 0`, so a profile
    /// bounded by the line `r = 0` has a boundary there — and the axis is
    /// *interior* to the revolved solid, not on its surface. The field then
    /// reads exactly zero along the axis, which is neither inside nor out.
    /// Draw the whole meridian section, from negative `r` to positive, and
    /// the axis becomes interior as it should.
    Revolve { region: Region2 },
    /// An imported or swept triangle surface, through its exact signed
    /// distance.
    Mesh(Arc<geometry::Sdf>),
    Union(Vec<Solid>),
    Intersection(Vec<Solid>),
    Difference(Box<Solid>, Box<Solid>),
    /// Quadratic polynomial blend of width `k` (a different solid from the
    /// sharp union: this is a fillet, not an approximation of one).
    SmoothUnion { a: Box<Solid>, b: Box<Solid>, k: f64 },
    /// Grow (`distance > 0`) or shrink the solid. Exact offsets of exact
    /// solids are the kernel's fillet and clearance operation.
    Offset { solid: Box<Solid>, distance: f64 },
    /// Wall of thickness `thickness` centred on the boundary.
    Shell { solid: Box<Solid>, thickness: f64 },
    Placed { solid: Box<Solid>, placement: Placement },
    /// Non-uniform scaling. Conservative, never exact: the field is divided
    /// by the largest factor.
    Stretch { solid: Box<Solid>, factors: P3 },
    /// Rotation about `+z` by `rate` radians per unit `z`, valid out to
    /// `radius` (the field is divided by `1 + rate * radius`).
    Twist { solid: Box<Solid>, rate: f64, radius: f64 },
    /// `counts` copies spaced by `period`, centred on the original.
    Array { solid: Box<Solid>, period: P3, counts: [usize; 3] },
    /// A triply periodic level set: the solid `f <= level`, or, with
    /// `thickness`, the sheet `|f - level| <= thickness/2`.
    ///
    /// The nodal function is divided by an analytic bound on its gradient,
    /// which keeps the kernel invariant but is conservative: the field
    /// grows more slowly than true distance, so the sheet's geometric wall
    /// is **thicker** than `thickness` by `1/|grad|`, a factor of about 1.6
    /// for the gyroid. `thickness` is therefore a lower bound on the wall;
    /// `kits::lattice::TpmsSpec::wall_thickness` measures the achieved
    /// value.
    Periodic { kind: Tpms, cell: P3, level: f64, thickness: Option<f64> },
    Custom(Arc<CustomField>),
}

const INF_BOX: (P3, P3) = ([f64::NEG_INFINITY; 3], [f64::INFINITY; 3]);

impl Solid {
    /// The kernel field: negative inside, 1-Lipschitz, zero on the
    /// boundary, and never larger in magnitude than the true distance.
    pub fn sdf(&self, p: P3) -> f64 {
        match self {
            Solid::Empty => f64::INFINITY,
            Solid::Full => f64::NEG_INFINITY,
            Solid::HalfSpace { normal, offset } => {
                let s = norm(*normal);
                if s == 0.0 {
                    f64::INFINITY
                } else {
                    dot(scale(*normal, 1.0 / s), p) - offset / s
                }
            }
            Solid::Sphere { center, radius } => norm(sub(p, *center)) - radius,
            Solid::Cuboid { center, half, round } => {
                let d = sub(p, *center);
                let q = [
                    d[0].abs() - (half[0] - round),
                    d[1].abs() - (half[1] - round),
                    d[2].abs() - (half[2] - round),
                ];
                let out = [q[0].max(0.0), q[1].max(0.0), q[2].max(0.0)];
                norm(out) + q[0].max(q[1]).max(q[2]).min(0.0) - round
            }
            Solid::Cylinder { a, b, radius } => {
                // Exact capped cylinder.
                let ba = sub(*b, *a);
                let pa = sub(p, *a);
                let baba = dot(ba, ba);
                let paba = dot(pa, ba);
                if baba == 0.0 {
                    return norm(pa) - radius;
                }
                let x = norm(sub(scale(pa, baba), scale(ba, paba))) - radius * baba;
                let y = (paba - 0.5 * baba).abs() - 0.5 * baba;
                let (x2, y2) = (x * x, y * y * baba);
                let d = if x.max(y) < 0.0 {
                    -x2.min(y2)
                } else {
                    (if x > 0.0 { x2 } else { 0.0 }) + (if y > 0.0 { y2 } else { 0.0 })
                };
                d.signum() * d.abs().sqrt() / baba
            }
            Solid::Cone { a, b, radius_a, radius_b } => {
                // Exact capped cone.
                let (ra, rb) = (*radius_a, *radius_b);
                let ba = sub(*b, *a);
                let baba = dot(ba, ba);
                if baba == 0.0 {
                    return norm(sub(p, *a)) - ra.max(rb);
                }
                let pa = sub(p, *a);
                let papa = dot(pa, pa);
                let paba = dot(pa, ba) / baba;
                let x = (papa - paba * paba * baba).max(0.0).sqrt();
                let cax = 0.0f64.max(x - if paba < 0.5 { ra } else { rb });
                let cay = (paba - 0.5).abs() - 0.5;
                let rba = rb - ra;
                let k = rba * rba + baba;
                let f = ((rba * (x - ra) + paba * baba) / k).clamp(0.0, 1.0);
                let cbx = x - ra - f * rba;
                let cby = paba - f;
                let s = if cbx < 0.0 && cay < 0.0 { -1.0 } else { 1.0 };
                s * (cax * cax + cay * cay * baba).min(cbx * cbx + cby * cby * baba).sqrt()
            }
            Solid::Capsule { a, b, radius } => {
                let ba = sub(*b, *a);
                let ee = dot(ba, ba);
                let t = if ee > 0.0 { (dot(sub(p, *a), ba) / ee).clamp(0.0, 1.0) } else { 0.0 };
                norm(sub(p, add(*a, scale(ba, t)))) - radius
            }
            Solid::Torus { center, major, minor } => {
                let d = sub(p, *center);
                let q = [(d[0] * d[0] + d[1] * d[1]).sqrt() - major, d[2]];
                (q[0] * q[0] + q[1] * q[1]).sqrt() - minor
            }
            Solid::Extrude { region, z0, z1 } => {
                let zc = 0.5 * (z0 + z1);
                let hh = 0.5 * (z1 - z0).abs();
                let d = region.sdf([p[0], p[1]]);
                let w = [d, (p[2] - zc).abs() - hh];
                let out = [w[0].max(0.0), w[1].max(0.0)];
                (out[0] * out[0] + out[1] * out[1]).sqrt() + w[0].max(w[1]).min(0.0)
            }
            Solid::Revolve { region } => region.sdf([(p[0] * p[0] + p[1] * p[1]).sqrt(), p[2]]),
            Solid::Mesh(sdf) => sdf.distance(p),
            Solid::Union(v) => v.iter().fold(f64::INFINITY, |d, s| d.min(s.sdf(p))),
            Solid::Intersection(v) => v.iter().fold(f64::NEG_INFINITY, |d, s| d.max(s.sdf(p))),
            Solid::Difference(a, b) => a.sdf(p).max(-b.sdf(p)),
            Solid::SmoothUnion { a, b, k } => {
                let (da, db) = (a.sdf(p), b.sdf(p));
                if *k <= 0.0 {
                    return da.min(db);
                }
                let h = (0.5 + 0.5 * (db - da) / k).clamp(0.0, 1.0);
                db + (da - db) * h - k * h * (1.0 - h)
            }
            Solid::Offset { solid, distance } => solid.sdf(p) - distance,
            Solid::Shell { solid, thickness } => solid.sdf(p).abs() - 0.5 * thickness,
            Solid::Placed { solid, placement } => placement.scale * solid.sdf(placement.inverse_apply(p)),
            Solid::Stretch { solid, factors } => {
                let m = factors[0].abs().min(factors[1].abs()).min(factors[2].abs());
                let q = [p[0] / factors[0], p[1] / factors[1], p[2] / factors[2]];
                m * solid.sdf(q)
            }
            Solid::Twist { solid, rate, radius } => {
                let a = -rate * p[2];
                let (s, c) = a.sin_cos();
                let q = [c * p[0] - s * p[1], s * p[0] + c * p[1], p[2]];
                solid.sdf(q) / (1.0 + rate.abs() * radius.abs())
            }
            Solid::Array { solid, period, counts } => {
                let idx = |k: usize| -> i64 {
                    if counts[k] <= 1 || period[k] == 0.0 {
                        0
                    } else {
                        (p[k] / period[k] + 0.5 * (counts[k] as f64 - 1.0)).round() as i64
                    }
                };
                let base = [idx(0), idx(1), idx(2)];
                let mut best = f64::INFINITY;
                for di in -1..=1i64 {
                    for dj in -1..=1i64 {
                        for dk in -1..=1i64 {
                            let c = [base[0] + di, base[1] + dj, base[2] + dk];
                            if (0..3).any(|k| c[k] < 0 || c[k] >= counts[k] as i64) {
                                continue;
                            }
                            let t: P3 = std::array::from_fn(|k| {
                                (c[k] as f64 - 0.5 * (counts[k] as f64 - 1.0)) * period[k]
                            });
                            best = best.min(solid.sdf(sub(p, t)));
                        }
                    }
                }
                best
            }
            Solid::Periodic { kind, cell, level, thickness } => {
                let (f, lip) = periodic_field(*kind, *cell, p);
                let g = (f - level) / lip;
                match thickness {
                    None => g,
                    Some(t) => g.abs() - 0.5 * t,
                }
            }
            Solid::Custom(c) => (c.f)(p),
        }
    }

    /// Whether [`Solid::sdf`] is the exact signed distance.
    pub fn exact(&self) -> bool {
        match self {
            Solid::Empty
            | Solid::Full
            | Solid::HalfSpace { .. }
            | Solid::Sphere { .. }
            | Solid::Cuboid { .. }
            | Solid::Cylinder { .. }
            | Solid::Cone { .. }
            | Solid::Capsule { .. }
            | Solid::Torus { .. }
            | Solid::Mesh(_) => true,
            Solid::Extrude { region, .. } | Solid::Revolve { region } => region.exact(),
            Solid::Offset { solid, .. } | Solid::Shell { solid, .. } => solid.exact(),
            Solid::Placed { solid, .. } => solid.exact(),
            Solid::Custom(c) => c.exact,
            Solid::Union(_)
            | Solid::Intersection(_)
            | Solid::Difference(_, _)
            | Solid::SmoothUnion { .. }
            | Solid::Stretch { .. }
            | Solid::Twist { .. }
            | Solid::Array { .. }
            | Solid::Periodic { .. } => false,
        }
    }

    /// A conservative axis-aligned bounding box. Infinite components mean
    /// the solid is unbounded and must be intersected with a box before it
    /// can be meshed.
    pub fn bbox(&self) -> (P3, P3) {
        match self {
            Solid::Empty => ([0.0; 3], [0.0; 3]),
            Solid::Full | Solid::HalfSpace { .. } | Solid::Periodic { .. } => INF_BOX,
            Solid::Sphere { center, radius } => ball(*center, *radius),
            Solid::Cuboid { center, half, .. } => (sub(*center, *half), add(*center, *half)),
            // A capped cylinder or cone ends in flat discs, so its box is
            // tight: the disc at `a` extends r sqrt(1 - d_k^2) along axis
            // k. Using a ball at each end would inflate the box along the
            // axis by the radius, which for a short wide part is most of
            // the grid a mesh would then have to cover.
            Solid::Cylinder { a, b, radius } => disc_box(*a, *b, *radius, *radius),
            Solid::Cone { a, b, radius_a, radius_b } => disc_box(*a, *b, *radius_a, *radius_b),
            Solid::Capsule { a, b, radius } => {
                let (l1, h1) = ball(*a, *radius);
                let (l2, h2) = ball(*b, *radius);
                (mins(l1, l2), maxs(h1, h2))
            }
            Solid::Torus { center, major, minor } => (
                [center[0] - major - minor, center[1] - major - minor, center[2] - minor],
                [center[0] + major + minor, center[1] + major + minor, center[2] + minor],
            ),
            Solid::Extrude { region, z0, z1 } => {
                let (lo, hi) = region.bbox();
                ([lo[0], lo[1], z0.min(*z1)], [hi[0], hi[1], z0.max(*z1)])
            }
            Solid::Revolve { region } => {
                let (lo, hi) = region.bbox();
                let r = hi[0].abs().max(lo[0].abs());
                ([-r, -r, lo[1]], [r, r, hi[1]])
            }
            Solid::Mesh(sdf) => sdf.mesh.bbox(),
            Solid::Union(v) => v
                .iter()
                .fold(([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]), |b, s| {
                    let (lo, hi) = s.bbox();
                    (mins(b.0, lo), maxs(b.1, hi))
                }),
            Solid::Intersection(v) => v
                .iter()
                .fold(INF_BOX, |b, s| {
                    let (lo, hi) = s.bbox();
                    (maxs(b.0, lo), mins(b.1, hi))
                }),
            Solid::Difference(a, _) => a.bbox(),
            Solid::SmoothUnion { a, b, k } => {
                let (l1, h1) = a.bbox();
                let (l2, h2) = b.bbox();
                let g = k.abs();
                (
                    [mins(l1, l2)[0] - g, mins(l1, l2)[1] - g, mins(l1, l2)[2] - g],
                    [maxs(h1, h2)[0] + g, maxs(h1, h2)[1] + g, maxs(h1, h2)[2] + g],
                )
            }
            Solid::Offset { solid, distance } => {
                let (lo, hi) = solid.bbox();
                let d = *distance;
                ([lo[0] - d, lo[1] - d, lo[2] - d], [hi[0] + d, hi[1] + d, hi[2] + d])
            }
            Solid::Shell { solid, thickness } => {
                let (lo, hi) = solid.bbox();
                let t = 0.5 * thickness.abs();
                ([lo[0] - t, lo[1] - t, lo[2] - t], [hi[0] + t, hi[1] + t, hi[2] + t])
            }
            Solid::Placed { solid, placement } => {
                let (lo, hi) = solid.bbox();
                let mut l = [f64::INFINITY; 3];
                let mut h = [f64::NEG_INFINITY; 3];
                for i in 0..8 {
                    let c = [
                        if i & 1 == 0 { lo[0] } else { hi[0] },
                        if i & 2 == 0 { lo[1] } else { hi[1] },
                        if i & 4 == 0 { lo[2] } else { hi[2] },
                    ];
                    let q = placement.apply(c);
                    l = mins(l, q);
                    h = maxs(h, q);
                }
                (l, h)
            }
            Solid::Stretch { solid, factors } => {
                let (lo, hi) = solid.bbox();
                let a: P3 = std::array::from_fn(|k| lo[k] * factors[k]);
                let b: P3 = std::array::from_fn(|k| hi[k] * factors[k]);
                (mins(a, b), maxs(a, b))
            }
            Solid::Twist { solid, radius, .. } => {
                let (lo, hi) = solid.bbox();
                let r = radius.abs();
                ([-r, -r, lo[2]], [r, r, hi[2]])
            }
            Solid::Array { solid, period, counts } => {
                let (lo, hi) = solid.bbox();
                let ext = |k: usize| 0.5 * (counts[k] as f64 - 1.0) * period[k].abs();
                (
                    [lo[0] - ext(0), lo[1] - ext(1), lo[2] - ext(2)],
                    [hi[0] + ext(0), hi[1] + ext(1), hi[2] + ext(2)],
                )
            }
            Solid::Custom(c) => c.bbox,
        }
    }

    /// Volume by midpoint quadrature of step `h`; second order in `h`.
    /// `bounds` must contain the solid.
    pub fn volume_quadrature(&self, bounds: (P3, P3), h: f64) -> f64 {
        use rayon::prelude::*;
        let n: [usize; 3] = std::array::from_fn(|k| (((bounds.1[k] - bounds.0[k]) / h).ceil() as usize).max(1));
        let inside: usize = (0..n[0])
            .into_par_iter()
            .map(|i| {
                let mut c = 0usize;
                for j in 0..n[1] {
                    for k in 0..n[2] {
                        let p = [
                            bounds.0[0] + h * (i as f64 + 0.5),
                            bounds.0[1] + h * (j as f64 + 0.5),
                            bounds.0[2] + h * (k as f64 + 0.5),
                        ];
                        if self.sdf(p) < 0.0 {
                            c += 1;
                        }
                    }
                }
                c
            })
            .sum();
        inside as f64 * h * h * h
    }

    /// Central-difference gradient of the field (the outward normal for an
    /// exact solid).
    pub fn gradient(&self, p: P3, h: f64) -> P3 {
        let mut g = [0.0; 3];
        for k in 0..3 {
            let mut a = p;
            let mut b = p;
            a[k] -= h;
            b[k] += h;
            g[k] = (self.sdf(b) - self.sdf(a)) / (2.0 * h);
        }
        g
    }

    // --- builders -------------------------------------------------------

    pub fn sphere(center: P3, radius: f64) -> Self {
        Solid::Sphere { center, radius }
    }
    pub fn cuboid(center: P3, half: P3) -> Self {
        Solid::Cuboid { center, half, round: 0.0 }
    }
    /// Axis-aligned box from two opposite corners.
    pub fn box_between(lo: P3, hi: P3) -> Self {
        Solid::Cuboid {
            center: scale(add(lo, hi), 0.5),
            half: scale(sub(hi, lo), 0.5),
            round: 0.0,
        }
    }
    pub fn rounded_box(center: P3, half: P3, round: f64) -> Self {
        Solid::Cuboid { center, half, round }
    }
    pub fn cylinder(a: P3, b: P3, radius: f64) -> Self {
        Solid::Cylinder { a, b, radius }
    }
    pub fn cylinder_z(center: P3, radius: f64, height: f64) -> Self {
        Solid::Cylinder {
            a: [center[0], center[1], center[2] - 0.5 * height],
            b: [center[0], center[1], center[2] + 0.5 * height],
            radius,
        }
    }
    pub fn cone(a: P3, b: P3, radius_a: f64, radius_b: f64) -> Self {
        Solid::Cone { a, b, radius_a, radius_b }
    }
    pub fn capsule(a: P3, b: P3, radius: f64) -> Self {
        Solid::Capsule { a, b, radius }
    }
    pub fn torus(center: P3, major: f64, minor: f64) -> Self {
        Solid::Torus { center, major, minor }
    }
    pub fn half_space(normal: P3, offset: f64) -> Self {
        Solid::HalfSpace { normal, offset }
    }
    pub fn extrude(region: Region2, z0: f64, z1: f64) -> Self {
        Solid::Extrude { region, z0, z1 }
    }
    pub fn revolve(region: Region2) -> Self {
        Solid::Revolve { region }
    }
    /// A sketch revolved through `angle` radians about `+z`, starting at the
    /// `+x` half-plane.
    pub fn revolve_arc(region: Region2, angle: f64) -> Self {
        let full = Solid::Revolve { region };
        if angle >= 2.0 * std::f64::consts::PI {
            return full;
        }
        let wedge = if angle <= std::f64::consts::PI {
            Solid::Intersection(vec![
                Solid::half_space([0.0, -1.0, 0.0], 0.0),
                Solid::half_space([angle.sin(), -angle.cos(), 0.0], 0.0),
            ])
        } else {
            Solid::Union(vec![
                Solid::half_space([0.0, -1.0, 0.0], 0.0),
                Solid::half_space([angle.sin(), -angle.cos(), 0.0], 0.0),
            ])
        };
        Solid::Intersection(vec![full, wedge])
    }
    pub fn mesh(sdf: geometry::Sdf) -> Self {
        Solid::Mesh(Arc::new(sdf))
    }
    pub fn from_trimesh(mesh: geometry::TriMesh) -> Self {
        Solid::Mesh(Arc::new(geometry::Sdf::new(mesh)))
    }
    pub fn custom(name: &str, bbox: (P3, P3), exact: bool, f: impl Fn(P3) -> f64 + Send + Sync + 'static) -> Self {
        Solid::Custom(Arc::new(CustomField {
            name: name.to_string(),
            f: Box::new(f),
            bbox,
            exact,
        }))
    }
    pub fn union(parts: Vec<Solid>) -> Self {
        Solid::Union(parts)
    }
    pub fn intersection(parts: Vec<Solid>) -> Self {
        Solid::Intersection(parts)
    }
    pub fn plus(self, other: Solid) -> Self {
        Solid::Union(vec![self, other])
    }
    pub fn minus(self, other: Solid) -> Self {
        Solid::Difference(Box::new(self), Box::new(other))
    }
    pub fn and(self, other: Solid) -> Self {
        Solid::Intersection(vec![self, other])
    }
    /// Union with a fillet of radius ~`k` along the seam.
    pub fn blend(self, other: Solid, k: f64) -> Self {
        Solid::SmoothUnion { a: Box::new(self), b: Box::new(other), k }
    }
    pub fn offset(self, distance: f64) -> Self {
        Solid::Offset { solid: Box::new(self), distance }
    }
    /// Round convex edges by `r`: a morphological opening.
    ///
    /// **This cannot be done by offsets alone.** In a distance algebra,
    /// eroding and dilating cancel exactly — `(f + r) - r = f` — so the
    /// eroded solid has to be re-distanced before it is grown back. That is
    /// what `cell` is for: the erosion is meshed at that resolution and its
    /// exact distance rebuilt. The result therefore carries the mesh's
    /// `O(cell^2)` error and costs one remesh, and `cell` must be a few
    /// times smaller than `r` for the rounding to be resolved at all.
    pub fn round_convex(self, r: f64, cell: f64) -> Result<Self, String> {
        let eroded = crate::meshing::mesh_solid(&self.offset(-r), cell)?;
        if !eroded.report().watertight {
            return Err(format!(
                "round_convex: the erosion by {r} did not mesh watertight at cell {cell}; \
                 the feature is probably thinner than 2*{r}"
            ));
        }
        Ok(Solid::from_trimesh(eroded).offset(r))
    }

    /// Fill concave edges with a fillet of radius `r`: a morphological
    /// closing, with the same re-distancing cost as
    /// [`Solid::round_convex`].
    pub fn round_concave(self, r: f64, cell: f64) -> Result<Self, String> {
        let dilated = crate::meshing::mesh_solid(&self.offset(r), cell)?;
        if !dilated.report().watertight {
            return Err(format!("round_concave: the dilation by {r} did not mesh watertight at cell {cell}"));
        }
        Ok(Solid::from_trimesh(dilated).offset(-r))
    }
    /// Hollow the solid, leaving a wall of `thickness`.
    pub fn shell(self, thickness: f64) -> Self {
        Solid::Shell { solid: Box::new(self), thickness }
    }
    pub fn placed(self, placement: Placement) -> Self {
        Solid::Placed { solid: Box::new(self), placement }
    }
    pub fn translated(self, t: P3) -> Self {
        self.placed(Placement::translation(t))
    }
    pub fn rotated(self, axis: P3, angle: f64) -> Self {
        self.placed(Placement::rotation(axis, angle))
    }
    pub fn scaled(self, s: f64) -> Self {
        self.placed(Placement::scaling(s))
    }
    pub fn mirrored(self, normal: P3) -> Self {
        self.placed(Placement::mirror(normal))
    }
    pub fn stretched(self, factors: P3) -> Self {
        Solid::Stretch { solid: Box::new(self), factors }
    }
    pub fn twisted(self, rate: f64, radius: f64) -> Self {
        Solid::Twist { solid: Box::new(self), rate, radius }
    }
    pub fn array(self, period: P3, counts: [usize; 3]) -> Self {
        Solid::Array { solid: Box::new(self), period, counts }
    }
    /// `count` copies rotated about `axis` through the origin.
    pub fn polar_array(self, axis: P3, count: usize) -> Self {
        Solid::Union(
            (0..count)
                .map(|i| {
                    let a = 2.0 * std::f64::consts::PI * i as f64 / count as f64;
                    self.clone().rotated(axis, a)
                })
                .collect(),
        )
    }
    /// Clip to a box, which is how an unbounded solid is made meshable.
    pub fn clipped(self, lo: P3, hi: P3) -> Self {
        self.and(Solid::box_between(lo, hi))
    }
}

/// A triply periodic nodal function and a bound on its gradient.
fn periodic_field(kind: Tpms, cell: P3, p: P3) -> (f64, f64) {
    let w: P3 = std::array::from_fn(|k| 2.0 * std::f64::consts::PI / cell[k]);
    let (x, y, z) = (w[0] * p[0], w[1] * p[1], w[2] * p[2]);
    let wn = (w[0] * w[0] + w[1] * w[1] + w[2] * w[2]).sqrt();
    match kind {
        // |d/dx| <= 2 w_x, so |grad| <= 2 |w|.
        Tpms::Gyroid => (
            x.sin() * y.cos() + y.sin() * z.cos() + z.sin() * x.cos(),
            2.0 * wn,
        ),
        // |d/dx| <= w_x.
        Tpms::SchwarzP => (x.cos() + y.cos() + z.cos(), wn),
        // |d/dx| <= 4 w_x.
        Tpms::Diamond => (
            x.sin() * y.sin() * z.sin()
                + x.sin() * y.cos() * z.cos()
                + x.cos() * y.sin() * z.cos()
                + x.cos() * y.cos() * z.sin(),
            4.0 * wn,
        ),
    }
}

/// Tight box of a capped cylinder or cone between `a` and `b`.
fn disc_box(a: P3, b: P3, ra: f64, rb: f64) -> (P3, P3) {
    let axis = sub(b, a);
    let len = norm(axis);
    if len < 1e-300 {
        return ball(a, ra.max(rb));
    }
    let d = scale(axis, 1.0 / len);
    let spread = |r: f64, k: usize| r * (1.0 - d[k] * d[k]).max(0.0).sqrt();
    let lo = std::array::from_fn(|k| (a[k] - spread(ra, k)).min(b[k] - spread(rb, k)));
    let hi = std::array::from_fn(|k| (a[k] + spread(ra, k)).max(b[k] + spread(rb, k)));
    (lo, hi)
}

fn ball(c: P3, r: f64) -> (P3, P3) {
    ([c[0] - r, c[1] - r, c[2] - r], [c[0] + r, c[1] + r, c[2] + r])
}
fn mins(a: P3, b: P3) -> P3 {
    std::array::from_fn(|k| a[k].min(b[k]))
}
fn maxs(a: P3, b: P3) -> P3 {
    std::array::from_fn(|k| a[k].max(b[k]))
}

/// Lift a 2D point to the sketch plane `z = 0`.
pub fn lift(p: V2) -> P3 {
    [p[0], p[1], 0.0]
}
