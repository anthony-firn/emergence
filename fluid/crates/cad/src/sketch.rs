//! 2D sketches: closed regions with an exact signed distance, exact area
//! and exact second moments.
//!
//! A [`Region2`] is the 2D half of the kernel. Everything that is drawn
//! flat and then given thickness — a plate outline, a copper pour, a gear
//! tooth, a mask layer, a yarn cross-section — starts here, and is lifted
//! to 3D by [`crate::solid::Solid::extrude`] or
//! [`crate::solid::Solid::revolve`], both of which preserve exactness.
//!
//! Every variant keeps the kernel invariant of [`crate::solid`]: `sdf` is
//! 1-Lipschitz, vanishes exactly on the region boundary, is negative
//! inside, and so satisfies `|sdf(p)| <= dist(p, boundary)`.
//! [`Region2::exact`] says when it *is* that distance.
//!
//! Polygon distance and the even-odd inside test are the textbook forms
//! (Ericson, *Real-Time Collision Detection*, 2005, §5.1); the area and
//! second moments come from Green's theorem and are exact for a polygon.

use serde::{Deserialize, Serialize};

/// A point or vector in the sketch plane.
pub type V2 = [f64; 2];

#[inline]
pub fn add2(a: V2, b: V2) -> V2 {
    [a[0] + b[0], a[1] + b[1]]
}
#[inline]
pub fn sub2(a: V2, b: V2) -> V2 {
    [a[0] - b[0], a[1] - b[1]]
}
#[inline]
pub fn scale2(a: V2, s: f64) -> V2 {
    [a[0] * s, a[1] * s]
}
#[inline]
pub fn dot2(a: V2, b: V2) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}
/// The z component of `a x b`.
#[inline]
pub fn cross2(a: V2, b: V2) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}
#[inline]
pub fn len2(a: V2) -> f64 {
    dot2(a, a).sqrt()
}
#[inline]
pub fn normalize2(a: V2) -> V2 {
    let n = len2(a);
    if n > 0.0 {
        scale2(a, 1.0 / n)
    } else {
        a
    }
}
/// `a` rotated by `angle` radians about the origin.
#[inline]
pub fn rot2(a: V2, angle: f64) -> V2 {
    let (s, c) = angle.sin_cos();
    [c * a[0] - s * a[1], s * a[0] + c * a[1]]
}

/// A closed polygon, optionally with holes.
///
/// Loops may be given in any order and any winding: inside/outside and
/// area both follow the even-odd rule, so a hole is simply a loop nested
/// inside another. Loops must not intersect one another.
/// [`Polygon::normalize_orientation`] imposes the outer-counter-clockwise,
/// holes-clockwise convention that export formats expect.
#[derive(Clone, Debug, Default, Deserialize, Serialize)]
pub struct Polygon {
    pub loops: Vec<Vec<V2>>,
}

impl Polygon {
    /// A single closed loop (the last point need not repeat the first).
    pub fn new(outer: Vec<V2>) -> Self {
        Self { loops: vec![outer] }
    }

    pub fn with_hole(mut self, hole: Vec<V2>) -> Self {
        self.loops.push(hole);
        self
    }

    /// Signed area of one loop (positive counter-clockwise), exact.
    pub fn loop_area(points: &[V2]) -> f64 {
        let n = points.len();
        let mut a = 0.0;
        for i in 0..n {
            a += cross2(points[i], points[(i + 1) % n]);
        }
        0.5 * a
    }

    /// Does one loop enclose `p`? (Even-odd, exact for a loop that does
    /// not self-intersect.)
    pub fn loop_contains(points: &[V2], p: V2) -> bool {
        let n = points.len();
        let mut inside = false;
        for i in 0..n {
            let (a, b) = (points[i], points[(i + 1) % n]);
            if (a[1] > p[1]) != (b[1] > p[1]) {
                let x = a[0] + (p[1] - a[1]) / (b[1] - a[1]) * (b[0] - a[0]);
                if p[0] < x {
                    inside = !inside;
                }
            }
        }
        inside
    }

    /// Enclosed area, exact and independent of how the loops are wound:
    /// each loop contributes `|signed area|` with the sign of its nesting
    /// depth, which matches the even-odd rule [`Polygon::contains`] uses.
    /// Loops must not intersect one another.
    pub fn area(&self) -> f64 {
        let mut total = 0.0;
        for (i, li) in self.loops.iter().enumerate() {
            if li.is_empty() {
                continue;
            }
            let depth = self
                .loops
                .iter()
                .enumerate()
                .filter(|&(j, lj)| j != i && Self::loop_contains(lj, li[0]))
                .count();
            let sign = if depth % 2 == 0 { 1.0 } else { -1.0 };
            total += sign * Self::loop_area(li).abs();
        }
        total
    }

    /// Set outer loops counter-clockwise and holes clockwise, in place.
    pub fn normalize_orientation(&mut self) {
        *self = self.oriented();
    }

    /// A copy with outer loops counter-clockwise and holes clockwise, so
    /// the signed sums in [`Polygon::centroid`] and
    /// [`Polygon::second_moments`] are consistent.
    pub fn oriented(&self) -> Polygon {
        let mut out = self.clone();
        for (i, li) in self.loops.iter().enumerate() {
            let depth = self
                .loops
                .iter()
                .enumerate()
                .filter(|&(j, lj)| j != i && !li.is_empty() && Self::loop_contains(lj, li[0]))
                .count();
            let want_ccw = depth % 2 == 0;
            if (Self::loop_area(li) > 0.0) != want_ccw {
                out.loops[i].reverse();
            }
        }
        out
    }

    /// Centroid of the enclosed area, exact.
    pub fn centroid(&self) -> V2 {
        let poly = self.oriented();
        let mut c = [0.0; 2];
        let mut a = 0.0;
        for l in &poly.loops {
            let n = l.len();
            for i in 0..n {
                let (p, q) = (l[i], l[(i + 1) % n]);
                let w = cross2(p, q);
                a += w;
                c = add2(c, scale2(add2(p, q), w));
            }
        }
        scale2(c, 1.0 / (3.0 * a))
    }

    /// Second moments of area `(Ixx, Iyy, Ixy)` about the centroid, exact.
    ///
    /// These are the section properties a beam model needs: `Ixx` bends
    /// about the x axis. Green's theorem on each loop
    /// (Steger, *CVGIP* 60(3), 1996, for the general moment form).
    pub fn second_moments(&self) -> (f64, f64, f64) {
        let poly = self.oriented();
        let c = poly.centroid();
        let (mut ixx, mut iyy, mut ixy) = (0.0, 0.0, 0.0);
        for l in &poly.loops {
            let n = l.len();
            for i in 0..n {
                let p = sub2(l[i], c);
                let q = sub2(l[(i + 1) % n], c);
                let w = cross2(p, q);
                ixx += w * (p[1] * p[1] + p[1] * q[1] + q[1] * q[1]);
                iyy += w * (p[0] * p[0] + p[0] * q[0] + q[0] * q[0]);
                ixy += w * (2.0 * p[0] * p[1] + p[0] * q[1] + q[0] * p[1] + 2.0 * q[0] * q[1]);
            }
        }
        (ixx / 12.0, iyy / 12.0, ixy / 24.0)
    }

    /// Total length of every loop, exact.
    pub fn perimeter(&self) -> f64 {
        self.loops
            .iter()
            .map(|l| {
                let n = l.len();
                (0..n).map(|i| len2(sub2(l[(i + 1) % n], l[i]))).sum::<f64>()
            })
            .sum()
    }

    /// Exact unsigned distance to the boundary.
    pub fn distance(&self, p: V2) -> f64 {
        let mut best = f64::INFINITY;
        for l in &self.loops {
            let n = l.len();
            for i in 0..n {
                let (a, b) = (l[i], l[(i + 1) % n]);
                let e = sub2(b, a);
                let ee = dot2(e, e);
                let t = if ee > 0.0 { (dot2(sub2(p, a), e) / ee).clamp(0.0, 1.0) } else { 0.0 };
                let d = len2(sub2(p, add2(a, scale2(e, t))));
                if d < best {
                    best = d;
                }
            }
        }
        best
    }

    /// Even-odd inside test (exact for loops that do not self-intersect).
    pub fn contains(&self, p: V2) -> bool {
        let mut inside = false;
        for l in &self.loops {
            let n = l.len();
            for i in 0..n {
                let (a, b) = (l[i], l[(i + 1) % n]);
                if (a[1] > p[1]) != (b[1] > p[1]) {
                    let x = a[0] + (p[1] - a[1]) / (b[1] - a[1]) * (b[0] - a[0]);
                    if p[0] < x {
                        inside = !inside;
                    }
                }
            }
        }
        inside
    }

    /// Exact signed distance (negative inside).
    pub fn sdf(&self, p: V2) -> f64 {
        let d = self.distance(p);
        if self.contains(p) {
            -d
        } else {
            d
        }
    }

    pub fn bbox(&self) -> (V2, V2) {
        let mut lo = [f64::INFINITY; 2];
        let mut hi = [f64::NEG_INFINITY; 2];
        for l in &self.loops {
            for p in l {
                for k in 0..2 {
                    lo[k] = lo[k].min(p[k]);
                    hi[k] = hi[k].max(p[k]);
                }
            }
        }
        (lo, hi)
    }
}

/// A closed 2D region.
#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(tag = "op", rename_all = "snake_case")]
pub enum Region2 {
    /// Nothing (the empty region).
    Empty,
    Circle { center: V2, radius: f64 },
    /// Rectangle with half-extents `half` and corner radius `round`.
    Rect {
        center: V2,
        half: V2,
        #[serde(default)]
        round: f64,
    },
    /// Regular `n`-gon of circumradius `r`, first vertex at angle `phase`.
    RegularPolygon { center: V2, radius: f64, sides: usize, #[serde(default)] phase: f64 },
    Polygon(Polygon),
    /// `{ p : dot(normal, p) - offset <= 0 }`; `normal` need not be unit.
    HalfPlane { normal: V2, offset: f64 },
    /// Capsule: points within `radius` of the segment `a`-`b`.
    Slot { a: V2, b: V2, radius: f64 },
    /// Annulus between `inner` and `outer`.
    Ring { center: V2, inner: f64, outer: f64 },
    Union { parts: Vec<Region2> },
    Intersection { parts: Vec<Region2> },
    /// `a` minus `b`.
    Difference { a: Box<Region2>, b: Box<Region2> },
    /// Grow (`distance > 0`) or shrink the region.
    Offset { region: Box<Region2>, distance: f64 },
    /// Shell of wall thickness `thickness` centred on the boundary.
    Shell { region: Box<Region2>, thickness: f64 },
    /// Rotate by `angle` about the origin, then translate.
    Placed {
        region: Box<Region2>,
        #[serde(default)]
        angle: f64,
        #[serde(default)]
        translate: V2,
    },
    /// `counts` copies spaced by `period`, centred on the original.
    Array { region: Box<Region2>, period: V2, counts: [usize; 2] },
    /// `count` copies equally spaced about the origin, evaluated by
    /// folding the query's angle into one sector instead of unioning every
    /// copy.
    ///
    /// This is how a gear's teeth, a bolt circle or an impeller's blades
    /// are drawn: three evaluations instead of `count`, exact whenever a
    /// copy stays within one sector of its own, which is the case for any
    /// feature that does not overlap its neighbour but one.
    PolarRepeat { region: Box<Region2>, count: usize },
}

impl Region2 {
    /// Signed distance: negative inside, 1-Lipschitz, zero exactly on the
    /// boundary, and `|sdf(p)| <= dist(p, boundary)` always.
    pub fn sdf(&self, p: V2) -> f64 {
        match self {
            Region2::Empty => f64::INFINITY,
            Region2::Circle { center, radius } => len2(sub2(p, *center)) - radius,
            Region2::Rect { center, half, round } => {
                let q = [
                    (p[0] - center[0]).abs() - (half[0] - round),
                    (p[1] - center[1]).abs() - (half[1] - round),
                ];
                let outside = [q[0].max(0.0), q[1].max(0.0)];
                len2(outside) + q[0].max(q[1]).min(0.0) - round
            }
            Region2::RegularPolygon { center, radius, sides, phase } => {
                Self::regular_polygon(*center, *radius, *sides, *phase).sdf(p)
            }
            Region2::Polygon(poly) => poly.sdf(p),
            Region2::HalfPlane { normal, offset } => {
                let n = normalize2(*normal);
                let scale = len2(*normal);
                if scale == 0.0 {
                    f64::INFINITY
                } else {
                    dot2(n, p) - offset / scale
                }
            }
            Region2::Slot { a, b, radius } => {
                let e = sub2(*b, *a);
                let ee = dot2(e, e);
                let t = if ee > 0.0 { (dot2(sub2(p, *a), e) / ee).clamp(0.0, 1.0) } else { 0.0 };
                len2(sub2(p, add2(*a, scale2(e, t)))) - radius
            }
            Region2::Ring { center, inner, outer } => {
                let r = len2(sub2(p, *center));
                (r - outer).max(inner - r)
            }
            Region2::Union { parts } => parts.iter().fold(f64::INFINITY, |d, r| d.min(r.sdf(p))),
            Region2::Intersection { parts } => parts.iter().fold(f64::NEG_INFINITY, |d, r| d.max(r.sdf(p))),
            Region2::Difference { a, b } => a.sdf(p).max(-b.sdf(p)),
            Region2::Offset { region, distance } => region.sdf(p) - distance,
            Region2::Shell { region, thickness } => region.sdf(p).abs() - 0.5 * thickness,
            Region2::Placed { region, angle, translate } => region.sdf(rot2(sub2(p, *translate), -angle)),
            Region2::Array { region, period, counts } => {
                let mut best = f64::INFINITY;
                // Nearest cell in each axis, plus one neighbour either side:
                // exact whenever the shape is contained in its own cell.
                let idx = |k: usize| -> i64 {
                    if counts[k] <= 1 || period[k] == 0.0 {
                        0
                    } else {
                        let off = 0.5 * (counts[k] as f64 - 1.0);
                        (p[k] / period[k] + off).round() as i64
                    }
                };
                let (i0, j0) = (idx(0), idx(1));
                for di in -1..=1 {
                    for dj in -1..=1 {
                        let (i, j) = (i0 + di, j0 + dj);
                        if i < 0 || j < 0 || i >= counts[0] as i64 || j >= counts[1] as i64 {
                            continue;
                        }
                        let t = [
                            (i as f64 - 0.5 * (counts[0] as f64 - 1.0)) * period[0],
                            (j as f64 - 0.5 * (counts[1] as f64 - 1.0)) * period[1],
                        ];
                        best = best.min(region.sdf(sub2(p, t)));
                    }
                }
                best
            }
            Region2::PolarRepeat { region, count } => {
                let count = (*count).max(1);
                let step = 2.0 * std::f64::consts::PI / count as f64;
                let k = (p[1].atan2(p[0]) / step).round();
                let mut best = f64::INFINITY;
                for d in -1..=1 {
                    best = best.min(region.sdf(rot2(p, -(k + d as f64) * step)));
                }
                best
            }
        }
    }

    /// Whether [`Region2::sdf`] is the exact signed distance rather than a
    /// conservative bound. Boolean and array nodes return `false`: their
    /// sign is exact, their magnitude can be short near a concave seam.
    pub fn exact(&self) -> bool {
        match self {
            Region2::Empty
            | Region2::Circle { .. }
            | Region2::Rect { .. }
            | Region2::RegularPolygon { .. }
            | Region2::Polygon(_)
            | Region2::HalfPlane { .. }
            | Region2::Slot { .. } => true,
            Region2::Ring { .. } => true,
            Region2::Offset { region, .. } | Region2::Shell { region, .. } => region.exact(),
            Region2::Placed { region, .. } => region.exact(),
            Region2::Union { .. }
            | Region2::Intersection { .. }
            | Region2::Difference { .. }
            | Region2::Array { .. }
            | Region2::PolarRepeat { .. } => false,
        }
    }

    /// Exact area where a closed form exists.
    pub fn area_exact(&self) -> Option<f64> {
        match self {
            Region2::Empty => Some(0.0),
            Region2::Circle { radius, .. } => Some(std::f64::consts::PI * radius * radius),
            Region2::Rect { half, round, .. } => {
                let (w, h) = (2.0 * half[0], 2.0 * half[1]);
                Some(w * h - (4.0 - std::f64::consts::PI) * round * round)
            }
            Region2::RegularPolygon { center, radius, sides, phase } => {
                Some(Self::regular_polygon(*center, *radius, *sides, *phase).area())
            }
            Region2::Polygon(p) => Some(p.area()),
            Region2::Ring { inner, outer, .. } => Some(std::f64::consts::PI * (outer * outer - inner * inner)),
            Region2::Slot { a, b, radius } => {
                Some(std::f64::consts::PI * radius * radius + 2.0 * radius * len2(sub2(*b, *a)))
            }
            Region2::Placed { region, .. } => region.area_exact(),
            Region2::Array { region, counts, .. } => {
                region.area_exact().map(|a| a * (counts[0] * counts[1]) as f64)
            }
            Region2::PolarRepeat { region, count } => {
                region.area_exact().map(|a| a * (*count).max(1) as f64)
            }
            _ => None,
        }
    }

    /// Area by midpoint quadrature on a grid of step `h` over the bounding
    /// box. Second order in `h` for a smooth boundary.
    pub fn area_quadrature(&self, h: f64) -> f64 {
        let (lo, hi) = self.bbox();
        let n: Vec<usize> = (0..2).map(|k| (((hi[k] - lo[k]) / h).ceil() as usize + 2).max(1)).collect();
        let mut a = 0.0;
        for i in 0..n[0] {
            for j in 0..n[1] {
                let p = [lo[0] + h * (i as f64 + 0.5) - 0.5 * h, lo[1] + h * (j as f64 + 0.5) - 0.5 * h];
                if self.sdf(p) < 0.0 {
                    a += h * h;
                }
            }
        }
        a
    }

    /// A conservative axis-aligned bounding box.
    pub fn bbox(&self) -> (V2, V2) {
        match self {
            Region2::Empty => ([0.0; 2], [0.0; 2]),
            Region2::Circle { center, radius } => (
                [center[0] - radius, center[1] - radius],
                [center[0] + radius, center[1] + radius],
            ),
            Region2::Rect { center, half, .. } => (sub2(*center, *half), add2(*center, *half)),
            Region2::RegularPolygon { center, radius, .. } => (
                [center[0] - radius, center[1] - radius],
                [center[0] + radius, center[1] + radius],
            ),
            Region2::Polygon(p) => p.bbox(),
            Region2::HalfPlane { .. } => ([f64::NEG_INFINITY; 2], [f64::INFINITY; 2]),
            Region2::Slot { a, b, radius } => (
                [a[0].min(b[0]) - radius, a[1].min(b[1]) - radius],
                [a[0].max(b[0]) + radius, a[1].max(b[1]) + radius],
            ),
            Region2::Ring { center, outer, .. } => (
                [center[0] - outer, center[1] - outer],
                [center[0] + outer, center[1] + outer],
            ),
            Region2::Union { parts } => parts.iter().fold(([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]), |b, r| {
                let (lo, hi) = r.bbox();
                ([b.0[0].min(lo[0]), b.0[1].min(lo[1])], [b.1[0].max(hi[0]), b.1[1].max(hi[1])])
            }),
            Region2::Intersection { parts } => parts.iter().fold(([f64::NEG_INFINITY; 2], [f64::INFINITY; 2]), |b, r| {
                let (lo, hi) = r.bbox();
                ([b.0[0].max(lo[0]), b.0[1].max(lo[1])], [b.1[0].min(hi[0]), b.1[1].min(hi[1])])
            }),
            Region2::Difference { a, .. } => a.bbox(),
            Region2::Offset { region, distance } => {
                let (lo, hi) = region.bbox();
                ([lo[0] - distance, lo[1] - distance], [hi[0] + distance, hi[1] + distance])
            }
            Region2::Shell { region, thickness } => {
                let (lo, hi) = region.bbox();
                let t = 0.5 * thickness;
                ([lo[0] - t, lo[1] - t], [hi[0] + t, hi[1] + t])
            }
            Region2::Placed { region, angle, translate } => {
                let (lo, hi) = region.bbox();
                let mut l = [f64::INFINITY; 2];
                let mut h = [f64::NEG_INFINITY; 2];
                for &x in &[lo[0], hi[0]] {
                    for &y in &[lo[1], hi[1]] {
                        let q = add2(rot2([x, y], *angle), *translate);
                        for k in 0..2 {
                            l[k] = l[k].min(q[k]);
                            h[k] = h[k].max(q[k]);
                        }
                    }
                }
                (l, h)
            }
            Region2::Array { region, period, counts } => {
                let (lo, hi) = region.bbox();
                let ext = |k: usize| 0.5 * (counts[k] as f64 - 1.0) * period[k].abs();
                ([lo[0] - ext(0), lo[1] - ext(1)], [hi[0] + ext(0), hi[1] + ext(1)])
            }
            Region2::PolarRepeat { region, .. } => {
                // A disc containing every rotated copy.
                let (lo, hi) = region.bbox();
                let r = [lo[0].abs(), hi[0].abs()]
                    .into_iter()
                    .flat_map(|x| [lo[1].abs(), hi[1].abs()].map(move |y| len2([x, y])))
                    .fold(0.0, f64::max);
                ([-r, -r], [r, r])
            }
        }
    }

    /// Vertices of a regular polygon as a [`Polygon`].
    pub fn regular_polygon(center: V2, radius: f64, sides: usize, phase: f64) -> Polygon {
        let pts = (0..sides.max(3))
            .map(|i| {
                let a = phase + 2.0 * std::f64::consts::PI * i as f64 / sides.max(3) as f64;
                [center[0] + radius * a.cos(), center[1] + radius * a.sin()]
            })
            .collect();
        Polygon::new(pts)
    }

    // --- builders -------------------------------------------------------

    pub fn circle(center: V2, radius: f64) -> Self {
        Region2::Circle { center, radius }
    }
    pub fn rect(center: V2, half: V2) -> Self {
        Region2::Rect { center, half, round: 0.0 }
    }
    pub fn rounded_rect(center: V2, half: V2, round: f64) -> Self {
        Region2::Rect { center, half, round }
    }
    pub fn polygon(points: Vec<V2>) -> Self {
        Region2::Polygon(Polygon::new(points))
    }
    pub fn slot(a: V2, b: V2, radius: f64) -> Self {
        Region2::Slot { a, b, radius }
    }
    pub fn ring(center: V2, inner: f64, outer: f64) -> Self {
        Region2::Ring { center, inner, outer }
    }
    pub fn union(parts: Vec<Region2>) -> Self {
        Region2::Union { parts }
    }
    pub fn intersection(parts: Vec<Region2>) -> Self {
        Region2::Intersection { parts }
    }
    pub fn minus(self, other: Region2) -> Self {
        Region2::Difference { a: Box::new(self), b: Box::new(other) }
    }
    pub fn plus(self, other: Region2) -> Self {
        Region2::Union { parts: vec![self, other] }
    }
    pub fn and(self, other: Region2) -> Self {
        Region2::Intersection { parts: vec![self, other] }
    }
    pub fn offset(self, distance: f64) -> Self {
        Region2::Offset { region: Box::new(self), distance }
    }
    pub fn shell(self, thickness: f64) -> Self {
        Region2::Shell { region: Box::new(self), thickness }
    }
    pub fn placed(self, angle: f64, translate: V2) -> Self {
        Region2::Placed { region: Box::new(self), angle, translate }
    }
    pub fn translated(self, t: V2) -> Self {
        self.placed(0.0, t)
    }
    pub fn rotated(self, angle: f64) -> Self {
        self.placed(angle, [0.0; 2])
    }
    pub fn array(self, period: V2, counts: [usize; 2]) -> Self {
        Region2::Array { region: Box::new(self), period, counts }
    }
    /// `count` copies of this region rotated about the origin.
    pub fn polar_repeat(self, count: usize) -> Self {
        Region2::PolarRepeat { region: Box::new(self), count }
    }
    /// `count` copies on a circle of radius `r` about the origin.
    pub fn polar_array(self, radius: f64, count: usize) -> Self {
        Region2::Union { parts: (0..count)
                .map(|i| {
                    let a = 2.0 * std::f64::consts::PI * i as f64 / count as f64;
                    self.clone().placed(a, [radius * a.cos(), radius * a.sin()])
                })
                .collect() }
    }
}
