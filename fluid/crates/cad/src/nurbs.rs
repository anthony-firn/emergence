//! NURBS curves and surfaces: the exact smooth geometry the rest of the
//! kernel lacks.
//!
//! Everything else in this crate is either an implicit field or a
//! polyline. Real CAD is neither: an airfoil, a cam profile, a duct
//! transition, a lens edge or a hull section is a rational spline, and
//! sampling one into a polygon early throws away the exactness that makes
//! it useful. This module keeps the spline, evaluates it exactly, and
//! samples only at the end — into a [`crate::Path3`] to sweep, a
//! [`crate::Polygon`] to extrude, or a [`TriMesh`] to mesh.
//!
//! The algorithms are the standard ones (Piegl & Tiller, *The NURBS Book*,
//! 2nd ed., 1997): the Cox–de Boor recurrence for the basis (§2.2), de
//! Boor's algorithm in homogeneous coordinates for evaluation (§4.4), the
//! derivative as a B-spline of one degree lower (§3.3), the exact rational
//! quadratic circle (§7.5) and the exact surface of revolution (§8.5,
//! A8.1).
//!
//! # What "exact" means here
//!
//! A rational quadratic with weights `1, cos(θ/2), 1` on a 90° span is a
//! circular arc **exactly**, not to within a tolerance. So
//! [`Nurbs::arc`] evaluated anywhere satisfies `|p − c| = r` to rounding,
//! and [`Surface::revolve`] of an arc is an exact sphere or torus. The
//! tests assert this to 1e-14, which is the property that separates a
//! rational spline from a polynomial one and is the reason CAD uses them.
//!
//! # What is not here
//!
//! Trimmed surfaces, surface-surface intersection, exact fillet blends and
//! STEP export all need a topological layer (faces bounded by trimmed
//! edges with tolerances) that this module does not have. A surface here
//! is a single untrimmed patch. Booleans still happen in the implicit
//! kernel, via [`crate::Solid::from_trimesh`].

use geometry::{add, cross, dot, norm, normalize, scale, sub, TriMesh, P3};

use crate::path::Path3;
use crate::sketch::{Polygon, V2};

/// A point in homogeneous coordinates, `(w x, w y, w z, w)`.
type H4 = [f64; 4];

#[inline]
fn h_add(a: H4, b: H4) -> H4 {
    std::array::from_fn(|k| a[k] + b[k])
}
#[inline]
fn h_sub(a: H4, b: H4) -> H4 {
    std::array::from_fn(|k| a[k] - b[k])
}
#[inline]
fn h_scale(a: H4, s: f64) -> H4 {
    std::array::from_fn(|k| a[k] * s)
}
#[inline]
fn h_point(p: P3, w: f64) -> H4 {
    [p[0] * w, p[1] * w, p[2] * w, w]
}
#[inline]
fn h_project(a: H4) -> P3 {
    [a[0] / a[3], a[1] / a[3], a[2] / a[3]]
}

/// A non-uniform rational B-spline curve.
///
/// `knots` has `control.len() + degree + 1` entries and must be
/// non-decreasing. A clamped curve (the usual kind) repeats its first and
/// last knot `degree + 1` times, which makes it pass through its first and
/// last control point.
#[derive(Clone, Debug)]
pub struct Nurbs {
    pub degree: usize,
    pub control: Vec<P3>,
    pub weights: Vec<f64>,
    pub knots: Vec<f64>,
}

impl Nurbs {
    /// Build and check a curve.
    pub fn new(degree: usize, control: Vec<P3>, weights: Vec<f64>, knots: Vec<f64>) -> Result<Self, String> {
        let n = control.len();
        if n == 0 {
            return Err("a curve needs control points".into());
        }
        if degree == 0 || degree >= n {
            return Err(format!("degree {degree} needs at least {} control points", degree + 1));
        }
        if weights.len() != n {
            return Err(format!("{} weights for {n} control points", weights.len()));
        }
        if weights.iter().any(|&w| !(w > 0.0)) {
            return Err("weights must be positive".into());
        }
        if knots.len() != n + degree + 1 {
            return Err(format!(
                "{} knots for {n} control points of degree {degree}: expected {}",
                knots.len(),
                n + degree + 1
            ));
        }
        if knots.windows(2).any(|w| w[1] < w[0]) {
            return Err("knots must be non-decreasing".into());
        }
        Ok(Self { degree, control, weights, knots })
    }

    /// A polynomial (all weights 1) clamped curve with uniform interior
    /// knots.
    pub fn clamped(degree: usize, control: Vec<P3>) -> Result<Self, String> {
        let n = control.len();
        let knots = clamped_knots(degree, n);
        let weights = vec![1.0; n];
        Self::new(degree, control, weights, knots)
    }

    /// A single Bézier segment through its control points.
    pub fn bezier(control: Vec<P3>) -> Result<Self, String> {
        let n = control.len();
        Self::clamped(n - 1, control)
    }

    /// A straight line as a degree-1 curve.
    pub fn line(a: P3, b: P3) -> Self {
        Self {
            degree: 1,
            control: vec![a, b],
            weights: vec![1.0, 1.0],
            knots: vec![0.0, 0.0, 1.0, 1.0],
        }
    }

    /// An exact circular arc of `quarters` quadrants, starting on the
    /// `u` axis and turning toward `v` (The NURBS Book §7.5). `v` need
    /// only be non-parallel to `u`; it is orthogonalised.
    ///
    /// Each quadrant is a rational quadratic Bézier with weights
    /// `1, √2/2, 1`, which represents a 90° arc exactly. `quarters` of 4
    /// closes the circle.
    pub fn arc(center: P3, u: P3, v: P3, radius: f64, quarters: usize) -> Result<Self, String> {
        if !(1..=4).contains(&quarters) {
            return Err("arc takes 1 to 4 quarters".into());
        }
        // `v` is orthogonalised against `u`: it names the direction of
        // turning, and a non-perpendicular pair would silently draw an
        // ellipse rather than a circle.
        let u = normalize(u);
        let v_perp = sub(v, scale(u, dot(u, v)));
        if norm(v_perp) < 1e-12 * norm(v).max(1e-30) {
            return Err("the arc's two directions are parallel".into());
        }
        let v = normalize(v_perp);
        // Corner and mid points of each quadrant on the unit circle.
        let ring: [V2; 9] = [
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
            [-1.0, 1.0],
            [-1.0, 0.0],
            [-1.0, -1.0],
            [0.0, -1.0],
            [1.0, -1.0],
            [1.0, 0.0],
        ];
        let count = 2 * quarters + 1;
        let control: Vec<P3> = ring[..count]
            .iter()
            .map(|c| add(center, add(scale(u, radius * c[0]), scale(v, radius * c[1]))))
            .collect();
        let s = 0.5f64.sqrt();
        let weights: Vec<f64> = (0..count).map(|i| if i % 2 == 0 { 1.0 } else { s }).collect();
        let mut knots = vec![0.0, 0.0, 0.0];
        for q in 1..quarters {
            let t = q as f64 / quarters as f64;
            knots.push(t);
            knots.push(t);
        }
        knots.extend([1.0, 1.0, 1.0]);
        Self::new(2, control, weights, knots)
    }

    /// A full exact circle.
    pub fn circle(center: P3, u: P3, v: P3, radius: f64) -> Self {
        Self::arc(center, u, v, radius, 4).expect("a full circle is well formed")
    }

    /// Interpolate `points` with a clamped curve of the given degree,
    /// using chord-length parameterisation (The NURBS Book §9.2.1).
    ///
    /// The curve passes through every point exactly; the tests check that
    /// to 1e-10.
    pub fn interpolate(points: &[P3], degree: usize) -> Result<Self, String> {
        let n = points.len();
        if n < degree + 1 {
            return Err(format!("interpolating degree {degree} needs {} points", degree + 1));
        }
        // Chord-length parameters.
        let mut chords = vec![0.0; n];
        for i in 1..n {
            chords[i] = chords[i - 1] + norm(sub(points[i], points[i - 1]));
        }
        let total = chords[n - 1];
        if total <= 0.0 {
            return Err("interpolation points are all coincident".into());
        }
        let t: Vec<f64> = chords.iter().map(|c| c / total).collect();
        // Averaged knot vector.
        let mut knots = vec![0.0; degree + 1];
        for j in 1..n - degree {
            knots.push(t[j..j + degree].iter().sum::<f64>() / degree as f64);
        }
        knots.extend(vec![1.0; degree + 1]);
        // Collocation matrix, solved by Gaussian elimination with partial
        // pivoting (the system is banded and well conditioned, and n is
        // small for a sketch curve).
        let basis = |i: usize, row: usize| -> f64 {
            let span = find_span(degree, &knots, t[row], n);
            let vals = basis_values(degree, &knots, span, t[row]);
            if i + degree >= span && i <= span {
                vals[i + degree - span]
            } else {
                0.0
            }
        };
        let mut a = vec![vec![0.0; n + 3]; n];
        for row in 0..n {
            for i in 0..n {
                a[row][i] = basis(i, row);
            }
            for k in 0..3 {
                a[row][n + k] = points[row][k];
            }
        }
        solve_in_place(&mut a, n, 3)?;
        let control: Vec<P3> = (0..n).map(|i| [a[i][n], a[i][n + 1], a[i][n + 2]]).collect();
        Self::new(degree, control, vec![1.0; n], knots)
    }

    pub fn domain(&self) -> (f64, f64) {
        (self.knots[self.degree], self.knots[self.knots.len() - 1 - self.degree])
    }

    /// Value of basis function `i` at `u`, for checking the basis itself.
    pub fn basis(&self, i: usize, u: f64) -> f64 {
        let n = self.control.len();
        let span = find_span(self.degree, &self.knots, u, n);
        let vals = basis_values(self.degree, &self.knots, span, u);
        if i + self.degree >= span && i <= span {
            vals[i + self.degree - span]
        } else {
            0.0
        }
    }

    /// The curve at parameter `u`, by de Boor's algorithm in homogeneous
    /// coordinates.
    pub fn eval(&self, u: f64) -> P3 {
        h_project(self.eval_homogeneous(u))
    }

    fn eval_homogeneous(&self, u: f64) -> H4 {
        self.homogeneous().eval(u)
    }

    /// This curve's control points in homogeneous coordinates.
    fn homogeneous(&self) -> HSpline {
        HSpline {
            degree: self.degree,
            control: (0..self.control.len())
                .map(|i| h_point(self.control[i], self.weights[i]))
                .collect(),
            knots: self.knots.clone(),
        }
    }

    /// First derivative of the rational curve, by the quotient rule
    /// `C' = (A' - w' C) / w`.
    ///
    /// `A'` and `w'` come from the homogeneous derivative spline, which is
    /// carried as raw 4-vectors rather than as (point, weight) pairs: for
    /// a polynomial curve every derivative weight is zero, and projecting
    /// through it would destroy the whole derivative.
    pub fn derivative(&self, u: f64) -> P3 {
        let h = self.homogeneous();
        let Some(dh) = h.derivative() else {
            return [0.0; 3];
        };
        let aw = h.eval(u);
        let w = aw[3];
        let c = h_project(aw);
        let d = dh.eval(u);
        std::array::from_fn(|k| (d[k] - d[3] * c[k]) / w)
    }

    /// Unit tangent at `u`.
    pub fn tangent(&self, u: f64) -> P3 {
        normalize(self.derivative(u))
    }

    /// Arc length by composite 5-point Gauss–Legendre quadrature,
    /// `per_span` intervals in each knot span.
    ///
    /// The subdivision follows the knots rather than the parameter range:
    /// the speed `|C'(u)|` has a kink at every knot, and a quadrature
    /// interval straddling one loses its order. Splitting at the knots
    /// first takes a circle from 5e-6 to 1e-14.
    pub fn length(&self, per_span: usize) -> f64 {
        const X: [f64; 5] = [
            -0.906_179_845_938_664,
            -0.538_469_310_105_683,
            0.0,
            0.538_469_310_105_683,
            0.906_179_845_938_664,
        ];
        const W: [f64; 5] = [
            0.236_926_885_056_189,
            0.478_628_670_499_366,
            0.568_888_888_888_889,
            0.478_628_670_499_366,
            0.236_926_885_056_189,
        ];
        let (lo, hi) = self.domain();
        let per_span = per_span.max(1);
        // Distinct knots inside the domain, which bound the smooth spans.
        let mut breaks = vec![lo];
        for &k in &self.knots {
            if k > *breaks.last().unwrap() + 1e-15 && k < hi - 1e-15 {
                breaks.push(k);
            }
        }
        breaks.push(hi);
        let mut total = 0.0;
        for w in breaks.windows(2) {
            let h = (w[1] - w[0]) / per_span as f64;
            for i in 0..per_span {
                let (a, b) = (w[0] + h * i as f64, w[0] + h * (i + 1) as f64);
                let (mid, half) = (0.5 * (a + b), 0.5 * (b - a));
                for k in 0..5 {
                    total += W[k] * norm(self.derivative(mid + half * X[k])) * half;
                }
            }
        }
        total
    }

    /// Insert knot `u` `times` times without changing the curve
    /// (The NURBS Book §5.2). Used to split a curve and to refine it, and
    /// a good check of the evaluator: the curve before and after must
    /// agree everywhere.
    pub fn insert_knot(&self, u: f64, times: usize) -> Result<Nurbs, String> {
        let (lo, hi) = self.domain();
        if u <= lo || u >= hi {
            return Err("can only insert an interior knot".into());
        }
        let mut out = self.clone();
        for _ in 0..times {
            let p = out.degree;
            let n = out.control.len();
            let span = find_span(p, &out.knots, u, n);
            let pw: Vec<H4> = (0..n).map(|i| h_point(out.control[i], out.weights[i])).collect();
            let mut new_pw = Vec::with_capacity(n + 1);
            new_pw.extend_from_slice(&pw[..span - p + 1]);
            for i in span - p + 1..=span {
                let denom = out.knots[i + p] - out.knots[i];
                let alpha = if denom.abs() > 0.0 { (u - out.knots[i]) / denom } else { 0.0 };
                new_pw.push(h_add(h_scale(pw[i - 1], 1.0 - alpha), h_scale(pw[i], alpha)));
            }
            new_pw.extend_from_slice(&pw[span..]);
            let mut knots = out.knots.clone();
            knots.insert(span + 1, u);
            out = Nurbs {
                degree: p,
                control: new_pw.iter().map(|&q| h_project(q)).collect(),
                weights: new_pw.iter().map(|q| q[3]).collect(),
                knots,
            };
        }
        Ok(out)
    }

    /// Sample into a path, uniformly in parameter.
    pub fn to_path(&self, samples: usize, closed: bool) -> Path3 {
        let (lo, hi) = self.domain();
        let n = samples.max(2);
        if closed {
            Path3::new((0..n).map(|i| self.eval(lo + (hi - lo) * i as f64 / n as f64)).collect(), true)
        } else {
            Path3::new(
                (0..=n).map(|i| self.eval(lo + (hi - lo) * i as f64 / n as f64)).collect(),
                false,
            )
        }
    }

    /// Sample into a path with points spaced equally in arc length, which
    /// is what a sweep wants: a uniform-parameter sampling of a rational
    /// curve bunches up where the weights are large.
    pub fn to_path_by_length(&self, samples: usize, closed: bool) -> Path3 {
        let n = samples.max(2);
        let (lo, hi) = self.domain();
        // Tabulate arc length on a fine grid, then invert by interpolation.
        let fine = (20 * n).max(400);
        let mut s = vec![0.0; fine + 1];
        let mut prev = self.eval(lo);
        for i in 1..=fine {
            let p = self.eval(lo + (hi - lo) * i as f64 / fine as f64);
            s[i] = s[i - 1] + norm(sub(p, prev));
            prev = p;
        }
        let total = s[fine];
        let at = |target: f64| -> P3 {
            let mut j = 0;
            while j < fine && s[j + 1] < target {
                j += 1;
            }
            let span = s[j + 1] - s[j];
            let f = if span > 0.0 { (target - s[j]) / span } else { 0.0 };
            let u = lo + (hi - lo) * (j as f64 + f) / fine as f64;
            self.eval(u)
        };
        let count = if closed { n } else { n + 1 };
        let divisor = if closed { n } else { n };
        Path3::new(
            (0..count).map(|i| at(total * i as f64 / divisor as f64)).collect(),
            closed,
        )
    }

    /// Sample a planar curve into a closed polygon in the `xy` plane, for
    /// extruding or revolving. The `z` component is dropped.
    ///
    /// Sampling is uniform in **parameter**, which for a rational curve is
    /// not uniform in arc length or in angle: the exact circle's parameter
    /// runs faster near each quadrant's corner control point. Use
    /// [`Nurbs::to_polygon_by_length`] where the spacing matters.
    pub fn to_polygon(&self, samples: usize) -> Polygon {
        let (lo, hi) = self.domain();
        let n = samples.max(3);
        Polygon::new(
            (0..n)
                .map(|i| {
                    let p = self.eval(lo + (hi - lo) * i as f64 / n as f64);
                    [p[0], p[1]]
                })
                .collect(),
        )
    }

    /// As [`Nurbs::to_polygon`], with the points spaced equally in arc
    /// length.
    pub fn to_polygon_by_length(&self, samples: usize) -> Polygon {
        let path = self.to_path_by_length(samples, true);
        Polygon::new(path.points.iter().map(|p| [p[0], p[1]]).collect())
    }

    /// Translate every control point.
    pub fn translated(&self, t: P3) -> Nurbs {
        Nurbs {
            control: self.control.iter().map(|&p| add(p, t)).collect(),
            ..self.clone()
        }
    }

    /// Control-polygon bounding box, which contains the curve (the convex
    /// hull property).
    pub fn control_box(&self) -> (P3, P3) {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for p in &self.control {
            for k in 0..3 {
                lo[k] = lo[k].min(p[k]);
                hi[k] = hi[k].max(p[k]);
            }
        }
        (lo, hi)
    }
}

/// A B-spline whose control points are homogeneous 4-vectors.
///
/// Keeping the 4-vectors raw is what makes the derivative work: the
/// derivative of a polynomial curve has every weight component zero, and
/// a (point, weight) representation cannot hold that.
#[derive(Clone, Debug)]
struct HSpline {
    degree: usize,
    control: Vec<H4>,
    knots: Vec<f64>,
}

impl HSpline {
    fn domain(&self) -> (f64, f64) {
        (self.knots[self.degree], self.knots[self.knots.len() - 1 - self.degree])
    }

    /// de Boor's algorithm.
    fn eval(&self, u: f64) -> H4 {
        let n = self.control.len();
        let p = self.degree;
        let (lo, hi) = self.domain();
        let u = u.clamp(lo, hi);
        if p == 0 {
            return self.control[find_span(0, &self.knots, u, n)];
        }
        let span = find_span(p, &self.knots, u, n);
        let mut d: Vec<H4> = (0..=p).map(|j| self.control[span - p + j]).collect();
        for r in 1..=p {
            for j in (r..=p).rev() {
                let i = span - p + j;
                let denom = self.knots[i + p + 1 - r] - self.knots[i];
                let alpha = if denom.abs() > 0.0 { (u - self.knots[i]) / denom } else { 0.0 };
                d[j] = h_add(h_scale(d[j - 1], 1.0 - alpha), h_scale(d[j], alpha));
            }
        }
        d[p]
    }

    /// The derivative spline: degree `p - 1` with control points
    /// `p (Q_{i+1} - Q_i) / (U_{i+p+1} - U_{i+1})` (The NURBS Book §3.3).
    fn derivative(&self) -> Option<HSpline> {
        let p = self.degree;
        if p < 1 || self.control.len() < 2 {
            return None;
        }
        let n = self.control.len();
        let control = (0..n - 1)
            .map(|i| {
                let denom = self.knots[i + p + 1] - self.knots[i + 1];
                let f = if denom.abs() > 0.0 { p as f64 / denom } else { 0.0 };
                h_scale(h_sub(self.control[i + 1], self.control[i]), f)
            })
            .collect();
        Some(HSpline {
            degree: p - 1,
            control,
            knots: self.knots[1..self.knots.len() - 1].to_vec(),
        })
    }
}

/// Clamped uniform knot vector for `n` control points of degree `p`.
pub fn clamped_knots(p: usize, n: usize) -> Vec<f64> {
    let mut knots = vec![0.0; p + 1];
    let interior = n.saturating_sub(p + 1);
    for i in 1..=interior {
        knots.push(i as f64 / (interior + 1) as f64);
    }
    knots.extend(vec![1.0; p + 1]);
    knots
}

/// The knot span containing `u` (The NURBS Book A2.1).
fn find_span(p: usize, knots: &[f64], u: f64, n: usize) -> usize {
    if u >= knots[n] {
        return n - 1;
    }
    if u <= knots[p] {
        return p;
    }
    let (mut lo, mut hi) = (p, n);
    let mut mid = (lo + hi) / 2;
    while u < knots[mid] || u >= knots[mid + 1] {
        if u < knots[mid] {
            hi = mid;
        } else {
            lo = mid;
        }
        mid = (lo + hi) / 2;
    }
    mid
}

/// The `p + 1` non-zero basis functions at `u` (The NURBS Book A2.2).
fn basis_values(p: usize, knots: &[f64], span: usize, u: f64) -> Vec<f64> {
    let mut n = vec![0.0; p + 1];
    let mut left = vec![0.0; p + 1];
    let mut right = vec![0.0; p + 1];
    n[0] = 1.0;
    for j in 1..=p {
        left[j] = u - knots[span + 1 - j];
        right[j] = knots[span + j] - u;
        let mut saved = 0.0;
        for r in 0..j {
            let denom = right[r + 1] + left[j - r];
            let temp = if denom.abs() > 0.0 { n[r] / denom } else { 0.0 };
            n[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        n[j] = saved;
    }
    n
}

/// Gaussian elimination with partial pivoting on an `n x (n + rhs)`
/// augmented matrix, leaving the solutions in the right-hand columns.
fn solve_in_place(a: &mut [Vec<f64>], n: usize, rhs: usize) -> Result<(), String> {
    for col in 0..n {
        let pivot = (col..n).max_by(|&i, &j| a[i][col].abs().partial_cmp(&a[j][col].abs()).unwrap()).unwrap();
        if a[pivot][col].abs() < 1e-300 {
            return Err("singular interpolation system".into());
        }
        a.swap(col, pivot);
        let d = a[col][col];
        for k in col..n + rhs {
            a[col][k] /= d;
        }
        for row in 0..n {
            if row == col {
                continue;
            }
            let f = a[row][col];
            if f == 0.0 {
                continue;
            }
            for k in col..n + rhs {
                a[row][k] -= f * a[col][k];
            }
        }
    }
    Ok(())
}

/// A NURBS surface patch: a tensor product of two curve bases.
///
/// `control[i][j]` is indexed `u` first. `weights` matches. The patch is
/// untrimmed: its boundary is its parameter rectangle.
#[derive(Clone, Debug)]
pub struct Surface {
    pub degree: [usize; 2],
    pub control: Vec<Vec<P3>>,
    pub weights: Vec<Vec<f64>>,
    pub knots_u: Vec<f64>,
    pub knots_v: Vec<f64>,
}

impl Surface {
    pub fn new(
        degree: [usize; 2],
        control: Vec<Vec<P3>>,
        weights: Vec<Vec<f64>>,
        knots_u: Vec<f64>,
        knots_v: Vec<f64>,
    ) -> Result<Self, String> {
        let nu = control.len();
        if nu == 0 {
            return Err("a surface needs control points".into());
        }
        let nv = control[0].len();
        if control.iter().any(|r| r.len() != nv) {
            return Err("the control net must be rectangular".into());
        }
        if weights.len() != nu || weights.iter().any(|r| r.len() != nv) {
            return Err("weights must match the control net".into());
        }
        if knots_u.len() != nu + degree[0] + 1 {
            return Err(format!("{} u knots, expected {}", knots_u.len(), nu + degree[0] + 1));
        }
        if knots_v.len() != nv + degree[1] + 1 {
            return Err(format!("{} v knots, expected {}", knots_v.len(), nv + degree[1] + 1));
        }
        Ok(Self { degree, control, weights, knots_u, knots_v })
    }

    pub fn domain_u(&self) -> (f64, f64) {
        (self.knots_u[self.degree[0]], self.knots_u[self.knots_u.len() - 1 - self.degree[0]])
    }
    pub fn domain_v(&self) -> (f64, f64) {
        (self.knots_v[self.degree[1]], self.knots_v[self.knots_v.len() - 1 - self.degree[1]])
    }

    /// The isoparametric curve at fixed `u`: a NURBS curve in `v`.
    pub fn iso_v(&self, u: f64) -> Nurbs {
        let nu = self.control.len();
        let nv = self.control[0].len();
        let p = self.degree[0];
        let span = find_span(p, &self.knots_u, u.clamp(self.domain_u().0, self.domain_u().1), nu);
        let vals = basis_values(p, &self.knots_u, span, u.clamp(self.domain_u().0, self.domain_u().1));
        let mut control = Vec::with_capacity(nv);
        let mut weights = Vec::with_capacity(nv);
        for j in 0..nv {
            let mut acc: H4 = [0.0; 4];
            for k in 0..=p {
                let i = span - p + k;
                acc = h_add(acc, h_scale(h_point(self.control[i][j], self.weights[i][j]), vals[k]));
            }
            control.push(h_project(acc));
            weights.push(acc[3]);
        }
        Nurbs { degree: self.degree[1], control, weights, knots: self.knots_v.clone() }
    }

    /// The surface point at `(u, v)`.
    pub fn eval(&self, u: f64, v: f64) -> P3 {
        self.iso_v(u).eval(v)
    }

    /// Unit normal at `(u, v)`, from the two parametric derivatives.
    /// Returns `None` at a degenerate point (a pole of a sphere).
    pub fn normal(&self, u: f64, v: f64) -> Option<P3> {
        let (lo_u, hi_u) = self.domain_u();
        let h = 1e-6 * (hi_u - lo_u).max(1e-12);
        let du = {
            let a = self.eval((u - h).max(lo_u), v);
            let b = self.eval((u + h).min(hi_u), v);
            sub(b, a)
        };
        let dv = self.iso_v(u).derivative(v);
        let n = cross(du, dv);
        if norm(n) < 1e-18 {
            None
        } else {
            Some(normalize(n))
        }
    }

    /// An exact surface of revolution about the `z` axis
    /// (The NURBS Book A8.1).
    ///
    /// `profile` must lie in the `xz` plane with `x >= 0`; `x` is the
    /// distance from the axis. Revolving an exact arc gives an exact
    /// sphere, torus or cone — which is what the tests check.
    pub fn revolve(profile: &Nurbs, quarters: usize) -> Result<Self, String> {
        if !(1..=4).contains(&quarters) {
            return Err("revolve takes 1 to 4 quarters".into());
        }
        if profile.control.iter().any(|p| p[0] < -1e-12 || p[1].abs() > 1e-12) {
            return Err("the profile must lie in the xz plane with x >= 0".into());
        }
        let ring: [V2; 9] = [
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
            [-1.0, 1.0],
            [-1.0, 0.0],
            [-1.0, -1.0],
            [0.0, -1.0],
            [1.0, -1.0],
            [1.0, 0.0],
        ];
        let count = 2 * quarters + 1;
        let s = 0.5f64.sqrt();
        let mut knots_u = vec![0.0, 0.0, 0.0];
        for q in 1..quarters {
            let t = q as f64 / quarters as f64;
            knots_u.push(t);
            knots_u.push(t);
        }
        knots_u.extend([1.0, 1.0, 1.0]);

        let nv = profile.control.len();
        let mut control = Vec::with_capacity(count);
        let mut weights = Vec::with_capacity(count);
        for i in 0..count {
            let cw = if i % 2 == 0 { 1.0 } else { s };
            let mut row = Vec::with_capacity(nv);
            let mut wrow = Vec::with_capacity(nv);
            for j in 0..nv {
                let r = profile.control[j][0];
                let z = profile.control[j][2];
                row.push([r * ring[i][0], r * ring[i][1], z]);
                wrow.push(cw * profile.weights[j]);
            }
            control.push(row);
            weights.push(wrow);
        }
        Self::new([2, profile.degree], control, weights, knots_u, profile.knots.clone())
    }

    /// Sweep `profile` along a straight line to `offset` (a general
    /// cylinder).
    pub fn extrude(profile: &Nurbs, offset: P3) -> Result<Self, String> {
        let nv = profile.control.len();
        let control = vec![
            profile.control.clone(),
            profile.control.iter().map(|&p| add(p, offset)).collect(),
        ];
        let weights = vec![profile.weights.clone(); 2];
        let _ = nv;
        Self::new([1, profile.degree], control, weights, vec![0.0, 0.0, 1.0, 1.0], profile.knots.clone())
    }

    /// A ruled surface between two curves (which must share a degree and
    /// knot vector; use [`Nurbs::insert_knot`] to make them compatible).
    pub fn ruled(a: &Nurbs, b: &Nurbs) -> Result<Self, String> {
        if a.degree != b.degree || a.knots.len() != b.knots.len() {
            return Err("ruled surfaces need two curves of the same degree and knot count".into());
        }
        let control = vec![a.control.clone(), b.control.clone()];
        let weights = vec![a.weights.clone(), b.weights.clone()];
        Self::new([1, a.degree], control, weights, vec![0.0, 0.0, 1.0, 1.0], a.knots.clone())
    }

    /// Tessellate into triangles, `nu` by `nv` quads. `wrap_u` joins the
    /// last column of stations to the first, which closes a full
    /// revolution without duplicated vertices.
    pub fn tessellate(&self, nu: usize, nv: usize, wrap_u: bool) -> TriMesh {
        let (lo_u, hi_u) = self.domain_u();
        let (lo_v, hi_v) = self.domain_v();
        let (nu, nv) = (nu.max(2), nv.max(2));
        let cols = if wrap_u { nu } else { nu + 1 };
        let mut vertices = Vec::with_capacity(cols * (nv + 1));
        for i in 0..cols {
            let u = lo_u + (hi_u - lo_u) * i as f64 / nu as f64;
            let iso = self.iso_v(u);
            for j in 0..=nv {
                let v = lo_v + (hi_v - lo_v) * j as f64 / nv as f64;
                vertices.push(iso.eval(v));
            }
        }
        let id = |i: usize, j: usize| ((i % cols) * (nv + 1) + j) as u32;
        let mut triangles = Vec::with_capacity(2 * nu * nv);
        let last = if wrap_u { cols } else { cols - 1 };
        for i in 0..last {
            for j in 0..nv {
                let (a, b, c, d) = (id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1));
                // Skip degenerate quads at a pole.
                if a != b && c != d {
                    triangles.push([a, b, c]);
                    triangles.push([a, c, d]);
                } else if a != b {
                    triangles.push([a, b, c]);
                } else if c != d {
                    triangles.push([a, c, d]);
                }
            }
        }
        let mut mesh = TriMesh::new(vertices, triangles);
        mesh.weld_relative(1e-10);
        mesh
    }

    /// Surface area by tessellation, which converges at second order.
    pub fn area(&self, nu: usize, nv: usize, wrap_u: bool) -> f64 {
        self.tessellate(nu, nv, wrap_u).area()
    }
}

