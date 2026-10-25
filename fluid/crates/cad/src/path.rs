//! Paths and swept solids: the boundary-representation half of the kernel.
//!
//! Some shapes are defined by a curve, not by a field: a yarn in a weave, a
//! screw thread, a bent wire, a pipe run, a coil. [`sweep`] carries a closed
//! 2D [`Profile`] along a [`Path3`] and produces triangles directly;
//! [`crate::Solid::from_trimesh`] turns the result back into an exact signed
//! distance, so swept solids compose with the implicit algebra.
//!
//! The frames are rotation-minimising (Bishop) frames, computed by
//! transporting one normal along the path with the minimal rotation that
//! takes each tangent to the next (Bishop, *Amer. Math. Monthly* 82(3),
//! 1975). A Frenet frame would flip at inflection points; this one does
//! not. For a closed path the residual twist on return is distributed
//! evenly over the stations, which is the usual closure convention.
//!
//! **Validity.** A sweep is self-intersection-free only while the section
//! stays clear of the path's centre of curvature: at each station, how far
//! the section reaches toward that centre must be less than the radius of
//! curvature there. [`Path3::curvature`] gives both, [`sweep_checked`]
//! takes the worst ratio over the path, and the generators that sweep
//! (threads, yarns) all use the checked form. This is the one place in the
//! crate where geometry could be silently wrong, so the indicator is not
//! optional.
//!
//! The test is the *directional* reach, not the section's largest radius:
//! a flattened section bending about its thin axis is perfectly valid at a
//! radius its width would forbid, which is exactly the case for a yarn in
//! a close-set fabric.

use std::f64::consts::{PI, TAU};

use geometry::{add, cross, dot, norm, normalize, scale, sub, TriMesh, P3};

use crate::sketch::V2;

/// A polyline path through space.
#[derive(Clone, Debug)]
pub struct Path3 {
    pub points: Vec<P3>,
    pub closed: bool,
}

impl Path3 {
    pub fn new(points: Vec<P3>, closed: bool) -> Self {
        Self { points, closed }
    }

    /// `n` segments from `a` to `b`.
    pub fn line(a: P3, b: P3, n: usize) -> Self {
        let n = n.max(1);
        Self::new(
            (0..=n)
                .map(|i| {
                    let t = i as f64 / n as f64;
                    add(scale(a, 1.0 - t), scale(b, t))
                })
                .collect(),
            false,
        )
    }

    /// Sample a parametric curve on `t in [0, 1]` with `n` segments.
    pub fn from_fn(n: usize, closed: bool, f: impl Fn(f64) -> P3) -> Self {
        let n = n.max(1);
        let last = if closed { n } else { n + 1 };
        Self::new((0..last).map(|i| f(i as f64 / n as f64)).collect(), closed)
    }

    /// Right-handed helix about `+z`: `turns` turns of `pitch` rise each,
    /// starting at angle `phase` and height `z0`.
    pub fn helix(radius: f64, pitch: f64, turns: f64, phase: f64, z0: f64, per_turn: usize) -> Self {
        let n = ((per_turn as f64 * turns).ceil() as usize).max(2);
        Self::from_fn(n, false, |t| {
            let a = phase + TAU * turns * t;
            [radius * a.cos(), radius * a.sin(), z0 + pitch * turns * t]
        })
    }

    /// Circular arc in the plane spanned by `u` and `v` about `center`.
    pub fn arc(center: P3, u: P3, v: P3, radius: f64, a0: f64, a1: f64, n: usize) -> Self {
        let (u, v) = (normalize(u), normalize(v));
        Self::from_fn(n, false, |t| {
            let a = a0 + (a1 - a0) * t;
            add(center, add(scale(u, radius * a.cos()), scale(v, radius * a.sin())))
        })
    }

    pub fn length(&self) -> f64 {
        let n = self.points.len();
        if n < 2 {
            return 0.0;
        }
        let last = if self.closed { n } else { n - 1 };
        (0..last).map(|i| norm(sub(self.points[(i + 1) % n], self.points[i]))).sum()
    }

    /// Smallest radius of curvature over the path, from the circumradius of
    /// each consecutive triple. `inf` for a straight path.
    pub fn min_curvature_radius(&self) -> f64 {
        let n = self.points.len();
        if n < 3 {
            return f64::INFINITY;
        }
        let mut best = f64::INFINITY;
        let range: Vec<usize> = if self.closed { (0..n).collect() } else { (1..n - 1).collect() };
        for i in range {
            let (a, b, c) = (
                self.points[(i + n - 1) % n],
                self.points[i],
                self.points[(i + 1) % n],
            );
            let (ab, bc, ca) = (norm(sub(b, a)), norm(sub(c, b)), norm(sub(a, c)));
            let area = 0.5 * norm(cross(sub(b, a), sub(c, a)));
            if area > 1e-30 {
                best = best.min(ab * bc * ca / (4.0 * area));
            }
        }
        best
    }

    /// Curvature at each station: the radius and the unit direction
    /// toward the centre of curvature. A straight stretch gives an
    /// infinite radius and an arbitrary direction.
    ///
    /// This is what decides whether a sweep self-intersects: the section
    /// must not reach past the centre of curvature, so its extent along
    /// the inward direction must stay below the radius. The section's
    /// largest radius is the wrong test for a flattened section bending
    /// about its thin axis, which is exactly the case for a yarn in a
    /// fabric.
    pub fn curvature(&self) -> Vec<(f64, P3)> {
        let n = self.points.len();
        (0..n)
            .map(|i| {
                if n < 3 {
                    return (f64::INFINITY, [0.0, 0.0, 1.0]);
                }
                let (im, ip) = if self.closed {
                    ((i + n - 1) % n, (i + 1) % n)
                } else if i == 0 {
                    (0, 1)
                } else if i == n - 1 {
                    (n - 2, n - 1)
                } else {
                    (i - 1, i + 1)
                };
                if im == i || ip == i {
                    return (f64::INFINITY, [0.0, 0.0, 1.0]);
                }
                let (a, b, c) = (self.points[im], self.points[i], self.points[ip]);
                let (ab, bc, ca) = (norm(sub(b, a)), norm(sub(c, b)), norm(sub(a, c)));
                let area = 0.5 * norm(cross(sub(b, a), sub(c, a)));
                if area <= 1e-30 {
                    return (f64::INFINITY, [0.0, 0.0, 1.0]);
                }
                let radius = ab * bc * ca / (4.0 * area);
                // Inward: the component of (a - b) + (c - b) perpendicular
                // to the tangent.
                let t = normalize(sub(c, a));
                let m = add(sub(a, b), sub(c, b));
                let inward = normalize(sub(m, scale(t, dot(m, t))));
                (radius, inward)
            })
            .collect()
    }

    /// Unit tangents at each station (one-sided at the ends of an open
    /// path, central inside).
    pub fn tangents(&self) -> Vec<P3> {
        let n = self.points.len();
        (0..n)
            .map(|i| {
                if self.closed {
                    normalize(sub(self.points[(i + 1) % n], self.points[(i + n - 1) % n]))
                } else if i == 0 {
                    normalize(sub(self.points[1], self.points[0]))
                } else if i == n - 1 {
                    normalize(sub(self.points[n - 1], self.points[n - 2]))
                } else {
                    normalize(sub(self.points[i + 1], self.points[i - 1]))
                }
            })
            .collect()
    }

    /// Rotation-minimising frames `(tangent, u, v)` at every station, with
    /// the starting `u` taken from the least-aligned coordinate axis.
    ///
    /// Which axis that is depends on the path, so a profile that is not
    /// symmetric about both of its own axes should be swept with
    /// [`Path3::frames_with_up`] instead, where the orientation is stated
    /// rather than inferred.
    pub fn frames(&self) -> Vec<[P3; 3]> {
        let t = self.tangents();
        if t.is_empty() {
            return Vec::new();
        }
        let k = (0..3).min_by(|&i, &j| t[0][i].abs().partial_cmp(&t[0][j].abs()).unwrap()).unwrap();
        let mut e = [0.0; 3];
        e[k] = 1.0;
        self.frames_with_up(e)
    }

    /// Frames whose `u` axis starts along `up` (projected perpendicular to
    /// the tangent) and is then transported along the path. Use this
    /// whenever the profile's two axes mean different things — a flattened
    /// yarn, a rectangular bar, a keyed shaft.
    pub fn frames_with_up(&self, up: P3) -> Vec<[P3; 3]> {
        let t = self.tangents();
        let n = t.len();
        if n == 0 {
            return Vec::new();
        }
        let mut a = up;
        if norm(sub(a, scale(t[0], dot(a, t[0])))) < 1e-12 * norm(a).max(1e-30) {
            // `up` is along the tangent: fall back to an axis that is not.
            let k = (0..3).min_by(|&i, &j| t[0][i].abs().partial_cmp(&t[0][j].abs()).unwrap()).unwrap();
            a = [0.0; 3];
            a[k] = 1.0;
        }
        let mut u = normalize(sub(a, scale(t[0], dot(a, t[0]))));
        let mut frames = Vec::with_capacity(n);
        frames.push([t[0], u, cross(t[0], u)]);
        for i in 1..n {
            u = transport(t[i - 1], t[i], u);
            frames.push([t[i], u, cross(t[i], u)]);
        }
        if self.closed && n > 2 {
            // Close the frame: distribute the residual twist.
            let u_wrapped = transport(t[n - 1], t[0], u);
            let f0 = &frames[0];
            let residual = dot(u_wrapped, f0[2]).atan2(dot(u_wrapped, f0[1]));
            for (i, f) in frames.iter_mut().enumerate() {
                let a = -residual * i as f64 / n as f64;
                let (s, c) = a.sin_cos();
                let (uu, vv) = (f[1], f[2]);
                f[1] = add(scale(uu, c), scale(vv, s));
                f[2] = cross(f[0], f[1]);
            }
        }
        frames
    }
}

/// Minimal rotation taking `t0` to `t1`, applied to `u`.
fn transport(t0: P3, t1: P3, u: P3) -> P3 {
    let axis = cross(t0, t1);
    let s = norm(axis);
    if s < 1e-14 {
        return if dot(t0, t1) > 0.0 { u } else { scale(u, -1.0) };
    }
    let axis = scale(axis, 1.0 / s);
    let angle = dot(t0, t1).clamp(-1.0, 1.0).acos();
    let (sa, ca) = angle.sin_cos();
    let parallel = scale(axis, dot(axis, u));
    let perp = sub(u, parallel);
    normalize(add(parallel, add(scale(perp, ca), scale(cross(axis, perp), sa))))
}

/// A closed cross-section in the frame's `(u, v)` coordinates, listed
/// counter-clockwise.
#[derive(Clone, Debug)]
pub struct Profile {
    pub points: Vec<V2>,
}

impl Profile {
    pub fn new(points: Vec<V2>) -> Self {
        Self { points }
    }

    pub fn circle(radius: f64, n: usize) -> Self {
        Self::ellipse(radius, radius, n)
    }

    pub fn ellipse(a: f64, b: f64, n: usize) -> Self {
        Self::new(
            (0..n.max(3))
                .map(|i| {
                    let t = TAU * i as f64 / n.max(3) as f64;
                    [a * t.cos(), b * t.sin()]
                })
                .collect(),
        )
    }

    pub fn rect(half_u: f64, half_v: f64) -> Self {
        Self::new(vec![
            [-half_u, -half_v],
            [half_u, -half_v],
            [half_u, half_v],
            [-half_u, half_v],
        ])
    }

    /// Lenticular (two circular arcs) cross-section of width `2a` and
    /// height `2b`: the standard flattened-yarn shape (Kemp, *J. Text.
    /// Inst.* 49, 1958, uses the race-track form; the lens is the limiting
    /// shape of a jammed yarn).
    pub fn lens(a: f64, b: f64, n: usize) -> Self {
        let n = n.max(4) / 2 * 2;
        // Arc through (-a,0), (0,b), (a,0): radius r, centre at (0, b - r).
        let r = (a * a + b * b) / (2.0 * b);
        let half = (a / r).asin();
        let mut pts = Vec::with_capacity(n);
        for i in 0..n / 2 {
            let t = -half + 2.0 * half * i as f64 / (n / 2) as f64;
            pts.push([r * t.sin(), b - r + r * t.cos()]);
        }
        for i in 0..n / 2 {
            let t = -half + 2.0 * half * i as f64 / (n / 2) as f64;
            pts.push([-r * t.sin(), r - b - r * t.cos()]);
        }
        Self::new(pts)
    }

    /// Race-track: a rectangle capped by half circles, width `2a`, height `2b`
    /// (`a >= b`).
    pub fn racetrack(a: f64, b: f64, n: usize) -> Self {
        let flat = (a - b).max(0.0);
        let n = n.max(8) / 4 * 4;
        let mut pts = Vec::new();
        for i in 0..n / 2 {
            let t = -PI / 2.0 + PI * i as f64 / (n / 2) as f64;
            pts.push([flat + b * t.cos(), b * t.sin()]);
        }
        for i in 0..n / 2 {
            let t = PI / 2.0 + PI * i as f64 / (n / 2) as f64;
            pts.push([-flat + b * t.cos(), b * t.sin()]);
        }
        Self::new(pts)
    }

    /// Signed area (positive counter-clockwise), exact.
    pub fn area(&self) -> f64 {
        crate::sketch::Polygon::loop_area(&self.points)
    }

    /// Largest distance from the profile origin to its boundary.
    pub fn outer_radius(&self) -> f64 {
        self.points.iter().map(|p| (p[0] * p[0] + p[1] * p[1]).sqrt()).fold(0.0, f64::max)
    }
}

/// Carry `profile` along `path`. An open path is capped flat at both ends;
/// a closed path yields a torus-like tube.
pub fn sweep(path: &Path3, profile: &Profile) -> TriMesh {
    sweep_framed(path, profile, &path.frames())
}

/// [`sweep`] with the profile's `u` axis along `up`.
pub fn sweep_oriented(path: &Path3, profile: &Profile, up: P3) -> TriMesh {
    sweep_framed(path, profile, &path.frames_with_up(up))
}

fn sweep_framed(path: &Path3, profile: &Profile, frames: &[[P3; 3]]) -> TriMesh {
    let n = frames.len();
    let m = profile.points.len();
    assert!(n >= 2 && m >= 3, "sweep needs at least 2 stations and 3 profile points");
    let mut vertices = Vec::with_capacity(n * m + 2);
    for (i, f) in frames.iter().enumerate() {
        for q in &profile.points {
            vertices.push(add(path.points[i], add(scale(f[1], q[0]), scale(f[2], q[1]))));
        }
    }
    let mut triangles = Vec::with_capacity(2 * n * m);
    let stations = if path.closed { n } else { n - 1 };
    for i in 0..stations {
        let (a, b) = (i, (i + 1) % n);
        for j in 0..m {
            let j1 = (j + 1) % m;
            let (v00, v01) = ((a * m + j) as u32, (a * m + j1) as u32);
            let (v10, v11) = ((b * m + j) as u32, (b * m + j1) as u32);
            triangles.push([v00, v10, v11]);
            triangles.push([v00, v11, v01]);
        }
    }
    if !path.closed {
        let c0 = vertices.len() as u32;
        vertices.push(path.points[0]);
        let c1 = vertices.len() as u32;
        vertices.push(path.points[n - 1]);
        for j in 0..m {
            let j1 = (j + 1) % m;
            // The side quads traverse the first ring as j1 -> j and the
            // last as j -> j1, so the caps must traverse each the other
            // way or the surface is not consistently oriented.
            triangles.push([c0, j as u32, j1 as u32]);
            triangles.push([c1, ((n - 1) * m + j1) as u32, ((n - 1) * m + j) as u32]);
        }
    }
    let mut mesh = TriMesh::new(vertices, triangles);
    mesh.orient_outward();
    mesh
}

/// What a sweep's validity check found.
#[derive(Clone, Copy, Debug)]
pub struct SweepCheck {
    /// The worst ratio, over the path's stations, of how far the section
    /// reaches toward the centre of curvature to the radius of curvature
    /// there. Below 1 the sweep cannot self-intersect.
    pub curvature_ratio: f64,
    /// Which station was worst.
    pub worst_station: usize,
    pub watertight: bool,
    pub volume: f64,
}

/// [`sweep`] plus the self-intersection indicator and a watertightness
/// check.
pub fn sweep_checked(path: &Path3, profile: &Profile) -> (TriMesh, SweepCheck) {
    sweep_checked_framed(path, profile, path.frames())
}

/// [`sweep_checked`] with the profile's `u` axis along `up`.
pub fn sweep_checked_oriented(path: &Path3, profile: &Profile, up: P3) -> (TriMesh, SweepCheck) {
    sweep_checked_framed(path, profile, path.frames_with_up(up))
}

fn sweep_checked_framed(path: &Path3, profile: &Profile, frames: Vec<[P3; 3]>) -> (TriMesh, SweepCheck) {
    let mesh = sweep_framed(path, profile, &frames);
    let report = mesh.report();
    let curvature = path.curvature();
    let mut ratio: f64 = 0.0;
    let mut worst_station = 0usize;
    for (i, (radius, inward)) in curvature.iter().enumerate() {
        if !radius.is_finite() || i >= frames.len() {
            continue;
        }
        let (u, v) = (frames[i][1], frames[i][2]);
        let (du, dv) = (dot(u, *inward), dot(v, *inward));
        // How far the section reaches along the inward direction.
        let reach = profile
            .points
            .iter()
            .map(|q| q[0] * du + q[1] * dv)
            .fold(f64::NEG_INFINITY, f64::max);
        let r = reach / radius;
        if r > ratio {
            ratio = r;
            worst_station = i;
        }
    }
    (
        mesh,
        SweepCheck {
            curvature_ratio: ratio,
            worst_station,
            watertight: report.watertight,
            volume: report.volume,
        },
    )
}

/// Loft through a list of sections, each given as a path station index and
/// its own profile (all profiles must have the same point count).
pub fn loft(path: &Path3, profiles: &[Profile]) -> TriMesh {
    let frames = path.frames();
    let n = frames.len();
    assert_eq!(profiles.len(), n, "loft needs one profile per station");
    let m = profiles[0].points.len();
    assert!(profiles.iter().all(|p| p.points.len() == m), "profiles must match in point count");
    let mut vertices = Vec::with_capacity(n * m + 2);
    for (i, f) in frames.iter().enumerate() {
        for q in &profiles[i].points {
            vertices.push(add(path.points[i], add(scale(f[1], q[0]), scale(f[2], q[1]))));
        }
    }
    let mut triangles = Vec::new();
    let stations = if path.closed { n } else { n - 1 };
    for i in 0..stations {
        let (a, b) = (i, (i + 1) % n);
        for j in 0..m {
            let j1 = (j + 1) % m;
            triangles.push([(a * m + j) as u32, (b * m + j) as u32, (b * m + j1) as u32]);
            triangles.push([(a * m + j) as u32, (b * m + j1) as u32, (a * m + j1) as u32]);
        }
    }
    if !path.closed {
        let c0 = vertices.len() as u32;
        vertices.push(path.points[0]);
        let c1 = vertices.len() as u32;
        vertices.push(path.points[n - 1]);
        for j in 0..m {
            let j1 = (j + 1) % m;
            // The side quads traverse the first ring as j1 -> j and the
            // last as j -> j1, so the caps must traverse each the other
            // way or the surface is not consistently oriented.
            triangles.push([c0, j as u32, j1 as u32]);
            triangles.push([c1, ((n - 1) * m + j1) as u32, ((n - 1) * m + j) as u32]);
        }
    }
    let mut mesh = TriMesh::new(vertices, triangles);
    mesh.orient_outward();
    mesh
}
