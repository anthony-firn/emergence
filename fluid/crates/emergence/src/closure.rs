//! Closure tables: a fine model folded, once, into a cheap validated surrogate.
//!
//! A closure is the law a finer model obeys, computed once and consumed by a
//! coarser one (an equation of state, a cross section, a rate). Here it is a
//! piecewise Chebyshev interpolant of a vector-valued function of one
//! variable, refined adaptively until a held-out check passes.
//!
//! # Construction and the error claim
//!
//! On each patch the fine model is sampled at the `2n + 1` Chebyshev–Lobatto
//! points of degree `2n`. The even-indexed points are the Lobatto points of
//! degree `n`; the degree-`n` interpolant through them is evaluated at the
//! `n` odd (held-out) points and compared with the fine model there. If the
//! largest held-out error is within tolerance the patch is accepted and the
//! degree-`2n` interpolant through all points is stored; otherwise the patch
//! is halved. The stored claim is the held-out error of the *degree-n*
//! interpolant, which for a function resolved by the patch overestimates the
//! stored degree-`2n` interpolant's error by orders of magnitude (spectral
//! convergence). This is the standard nested (Clenshaw–Curtis-style) error
//! estimate.
//!
//! The claim is [`Rigor::Validated`](crate::Rigor::Validated), not rigorous:
//! a feature narrower than the sample spacing (an undamped Mie resonance, a
//! cusp) can hide between samples. Patches that still fail at the minimum
//! width are kept as *unconverged* and report out-of-domain, so the ladder
//! unfolds to the fine model there.

use rayon::prelude::*;

use crate::claim::{Assessment, Claim};
use crate::rung::Rung;

/// The interpolation variable.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Axis {
    Linear,
    /// Interpolate in `ln x` (for quantities spanning decades).
    Log,
}

impl Axis {
    fn to_u(self, x: f64) -> f64 {
        match self {
            Axis::Linear => x,
            Axis::Log => x.ln(),
        }
    }
    fn from_u(self, u: f64) -> f64 {
        match self {
            Axis::Linear => u,
            Axis::Log => u.exp(),
        }
    }
}

/// Construction options.
#[derive(Clone, Copy, Debug)]
pub struct ClosureOptions {
    /// Degree `n` of the checked interpolant; the stored one has degree `2n`.
    pub degree: usize,
    /// Smallest patch width, in axis units, before a patch is given up.
    pub min_width: f64,
    /// Stop refining (leaving the rest unconverged) beyond this many patches.
    pub max_patches: usize,
}

impl Default for ClosureOptions {
    fn default() -> Self {
        Self {
            degree: 8,
            min_width: 1e-6,
            max_patches: 1 << 16,
        }
    }
}

#[derive(Clone, Debug)]
struct Patch<const K: usize> {
    a: f64,
    b: f64,
    values: Vec<[f64; K]>,
    error: f64,
    converged: bool,
}

/// A validated piecewise-Chebyshev closure of `f: R → R^K`.
///
/// Stored flat for lookup speed: patch ends in one contiguous array (the
/// binary search touches nothing else) and each patch's interpolant as
/// Chebyshev coefficients (DCT-I of the Lobatto values, the same polynomial),
/// evaluated by Clenshaw's recurrence without divisions.
#[derive(Clone, Debug)]
pub struct Closure<const K: usize> {
    axis: Axis,
    lo: f64,
    hi: f64,
    nodes: Vec<f64>,
    starts: Vec<f64>,
    ends: Vec<f64>,
    /// Claimed error per patch; `NaN` marks an unconverged patch.
    errors: Vec<f64>,
    /// Chebyshev coefficients, `nodes.len()` per patch (`c_0`, `c_N` halved).
    coeffs: Vec<[f64; K]>,
    /// Uniform bucket index: `first[b]` is the first patch reaching bucket `b`.
    first: Vec<u32>,
    bucket: f64,
    /// Interpolation tolerance the closure was built to (absolute, per
    /// component); claims add the samples' own error.
    pub tol: f64,
    /// Fine-model evaluations spent building it.
    pub samples: usize,
}

/// Chebyshev–Lobatto points `cos(π j / N)`, `j = 0..=N`.
fn lobatto(n: usize) -> Vec<f64> {
    (0..=n).map(|j| (std::f64::consts::PI * j as f64 / n as f64).cos()).collect()
}

/// Barycentric interpolation (second kind) at `t ∈ [-1, 1]`; `stride` picks
/// every `stride`-th node (stride 2 = the nested coarse interpolant).
fn barycentric<const K: usize>(t: f64, nodes: &[f64], values: &[[f64; K]], stride: usize) -> [f64; K] {
    let n = (nodes.len() - 1) / stride;
    let mut num = [0.0; K];
    let mut den = 0.0;
    for k in 0..=n {
        let j = k * stride;
        let d = t - nodes[j];
        if d == 0.0 {
            return values[j];
        }
        let mut w = if k % 2 == 0 { 1.0 } else { -1.0 };
        if k == 0 || k == n {
            w *= 0.5;
        }
        let c = w / d;
        den += c;
        for (acc, v) in num.iter_mut().zip(values[j].iter()) {
            *acc += c * v;
        }
    }
    num.map(|v| v / den)
}

impl<const K: usize> Closure<K> {
    /// Build a closure of `f` on `[lo, hi]` to absolute tolerance `tol` in
    /// every component. Sampling of each refinement level runs in parallel.
    pub fn build<F>(f: F, lo: f64, hi: f64, axis: Axis, tol: f64, opt: ClosureOptions) -> Self
    where
        F: Fn(f64) -> [f64; K] + Sync,
    {
        Self::build_resolving(f, lo, hi, axis, tol, opt, |_, _| f64::INFINITY)
    }

    /// As [`Closure::build`], but a patch `[x_a, x_b]` is accepted only if,
    /// besides the held-out check, its largest node spacing (in `x`) is at
    /// most `resolution(x_a, x_b)`. `resolution` encodes physics the samples
    /// cannot reveal: the narrowest feature the function has in (or next to)
    /// the interval, e.g. a resonance linewidth or a boundary-layer width.
    /// Without it, a feature narrower than the spacing can pass the held-out
    /// check unseen.
    pub fn build_resolving<F, R>(f: F, lo: f64, hi: f64, axis: Axis, tol: f64, opt: ClosureOptions, resolution: R) -> Self
    where
        F: Fn(f64) -> [f64; K] + Sync,
        R: Fn(f64, f64) -> f64,
    {
        Self::build_inexact(|x| (f(x), 0.0), lo, hi, axis, tol, opt, resolution)
    }

    /// The general builder: `f` returns each sample *and its own claimed
    /// error* (the fine model is rarely exact). A patch's claim is its
    /// held-out interpolation error plus the largest sample error on it, and
    /// only the interpolation error is held to `tol`.
    pub fn build_inexact<F, R>(f: F, lo: f64, hi: f64, axis: Axis, tol: f64, opt: ClosureOptions, resolution: R) -> Self
    where
        F: Fn(f64) -> ([f64; K], f64) + Sync,
        R: Fn(f64, f64) -> f64,
    {
        assert!(hi > lo && opt.degree >= 2);
        let n_fine = 2 * opt.degree;
        let nodes = lobatto(n_fine);
        let (ulo, uhi) = (axis.to_u(lo), axis.to_u(hi));
        let mut pending = vec![(ulo, uhi)];
        let mut patches = Vec::new();
        let mut samples = 0;
        while !pending.is_empty() {
            let sampled: Vec<(f64, f64, Vec<[f64; K]>, f64)> = pending
                .par_iter()
                .map(|&(a, b)| {
                    let (mid, half) = (0.5 * (a + b), 0.5 * (b - a));
                    let (v, e): (Vec<[f64; K]>, Vec<f64>) = nodes.iter().map(|t| f(axis.from_u(mid + half * t))).unzip();
                    (a, b, v, e.into_iter().fold(0.0, f64::max))
                })
                .collect();
            samples += sampled.len() * nodes.len();
            let mut next = Vec::new();
            let budget_left = opt.max_patches.saturating_sub(patches.len() + sampled.len());
            for (a, b, values, sample_error) in sampled {
                let error = held_out_error(&nodes, &values);
                let ok = error <= tol && resolves(axis, &nodes, a, b, &resolution);
                if ok || b - a <= opt.min_width || next.len() + 2 > budget_left {
                    patches.push(Patch {
                        a,
                        b,
                        values,
                        error: error + sample_error,
                        converged: ok,
                    });
                } else {
                    let m = 0.5 * (a + b);
                    next.push((a, m));
                    next.push((m, b));
                }
            }
            pending = next;
        }
        patches.sort_by(|p, q| p.a.total_cmp(&q.a));
        let ends: Vec<f64> = patches.iter().map(|p| p.b).collect();
        let n_buckets = 2 * patches.len().max(1);
        let bucket = (uhi - ulo) / n_buckets as f64;
        let first = (0..=n_buckets)
            .map(|b| ends.partition_point(|&e| e < ulo + b as f64 * bucket).min(ends.len().saturating_sub(1)) as u32)
            .collect();
        Self {
            axis,
            lo: ulo,
            hi: uhi,
            first,
            bucket,
            starts: patches.iter().map(|p| p.a).collect(),
            ends,
            errors: patches.iter().map(|p| if p.converged { p.error } else { f64::NAN }).collect(),
            coeffs: patches.iter().flat_map(|p| chebyshev_coefficients(&p.values)).collect(),
            nodes,
            tol,
            samples,
        }
    }

    /// The interpolated value and its claimed error, or `None` outside the
    /// domain or on an unconverged patch.
    #[inline]
    pub fn eval(&self, x: f64) -> Option<([f64; K], f64)> {
        let (i, u) = self.locate(x)?;
        let (a, b) = (self.starts[i], self.ends[i]);
        let t = ((2.0 * u - a - b) / (b - a)).clamp(-1.0, 1.0);
        let stride = self.nodes.len();
        let c = &self.coeffs[i * stride..(i + 1) * stride];
        // Clenshaw: b_k = c_k + 2t b_{k+1} − b_{k+2}; p = c_0 + t b_1 − b_2.
        let (mut b1, mut b2) = ([0.0; K], [0.0; K]);
        let t2 = 2.0 * t;
        for ck in c[1..].iter().rev() {
            for k in 0..K {
                let b0 = ck[k] + t2 * b1[k] - b2[k];
                b2[k] = b1[k];
                b1[k] = b0;
            }
        }
        let v = std::array::from_fn(|k| c[0][k] + t * b1[k] - b2[k]);
        Some((v, self.errors[i]))
    }

    /// The claimed error at `x` without interpolating (the indicator).
    #[inline]
    pub fn claim_at(&self, x: f64) -> Option<f64> {
        self.locate(x).map(|(i, _)| self.errors[i])
    }

    /// Patch index and axis coordinate of `x`, if on a converged patch.
    #[inline]
    fn locate(&self, x: f64) -> Option<(usize, f64)> {
        let u = self.axis.to_u(x);
        if !(u >= self.lo && u <= self.hi) {
            return None;
        }
        let b = (((u - self.lo) / self.bucket) as usize).min(self.first.len() - 2);
        let (i0, i1) = (self.first[b] as usize, self.first[b + 1] as usize + 1);
        let mut i = i0 + self.ends[i0..i1.min(self.ends.len())].partition_point(|&e| e < u);
        // Rounding in the bucket index can land one bucket off at a boundary:
        // fall back to the full search if the patch does not contain u.
        if i >= self.ends.len() || u < self.starts[i] || u > self.ends[i] {
            i = self.ends.partition_point(|&e| e < u);
        }
        let e = *self.errors.get(i)?;
        (!e.is_nan()).then_some((i, u))
    }

    pub fn patches(&self) -> usize {
        self.ends.len()
    }

    /// Fraction of the domain (in axis units) that converged.
    pub fn coverage(&self) -> f64 {
        let ok: f64 = (0..self.patches())
            .filter(|&i| !self.errors[i].is_nan())
            .map(|i| self.ends[i] - self.starts[i])
            .sum();
        ok / (self.hi - self.lo)
    }

    /// The unconverged intervals, in `x`.
    pub fn unconverged(&self) -> Vec<(f64, f64)> {
        (0..self.patches())
            .filter(|&i| self.errors[i].is_nan())
            .map(|i| (self.axis.from_u(self.starts[i]), self.axis.from_u(self.ends[i])))
            .collect()
    }

    /// Nodes per patch (interpolation work per evaluation is linear in this).
    pub fn nodes_per_patch(&self) -> usize {
        self.nodes.len()
    }
}

/// Chebyshev coefficients of the interpolant through values at the Lobatto
/// points `cos(π j/N)`: `c_k = (2/N) Σ''_j f_j cos(π j k/N)` (DCT-I), with
/// `c_0` and `c_N` halved so that `p(t) = Σ_k c_k T_k(t)`.
fn chebyshev_coefficients<const K: usize>(values: &[[f64; K]]) -> Vec<[f64; K]> {
    let n = values.len() - 1;
    (0..=n)
        .map(|k| {
            let mut c = [0.0; K];
            for (j, v) in values.iter().enumerate() {
                let w = if j == 0 || j == n { 0.5 } else { 1.0 };
                let cos = (std::f64::consts::PI * (j * k) as f64 / n as f64).cos();
                for (ck, vk) in c.iter_mut().zip(v) {
                    *ck += w * vk * cos;
                }
            }
            let scale = if k == 0 || k == n { 1.0 / n as f64 } else { 2.0 / n as f64 };
            c.map(|v| v * scale)
        })
        .collect()
}

/// Whether the patch `[a, b]` (axis units) has node spacing within
/// `resolution(x_a, x_b)`.
fn resolves<R: Fn(f64, f64) -> f64>(axis: Axis, nodes: &[f64], a: f64, b: f64, resolution: &R) -> bool {
    let (mid, half) = (0.5 * (a + b), 0.5 * (b - a));
    let xs: Vec<f64> = nodes.iter().map(|t| axis.from_u(mid + half * t)).collect();
    let spacing = xs.windows(2).map(|w| (w[0] - w[1]).abs()).fold(0.0, f64::max);
    spacing <= resolution(axis.from_u(a), axis.from_u(b))
}

/// Largest held-out error of the nested degree-`n` interpolant, with a
/// round-off floor relative to the patch's magnitude.
fn held_out_error<const K: usize>(nodes: &[f64], values: &[[f64; K]]) -> f64 {
    let mut err: f64 = 0.0;
    let mut scale: f64 = 0.0;
    for v in values {
        for c in v {
            scale = scale.max(c.abs());
        }
    }
    for j in (1..nodes.len()).step_by(2) {
        let p = barycentric(nodes[j], nodes, values, 2);
        for (a, b) in p.iter().zip(values[j].iter()) {
            let e = (a - b).abs();
            err = err.max(if e.is_nan() { f64::INFINITY } else { e });
        }
    }
    // Round-off in the barycentric sum grows like N ε · scale.
    err.max(4.0 * nodes.len() as f64 * f64::EPSILON * scale)
}

type KeyFn<I> = Box<dyn Fn(&I) -> Option<f64> + Send + Sync>;
type WrapFn<I, O, const K: usize> = Box<dyn Fn(&I, [f64; K]) -> O + Send + Sync>;

/// A closure on a ladder. `key` maps an input to the closure's variable, or
/// `None` when the input is not in the closure's family (another material,
/// another wavelength); `wrap` turns the interpolated vector into an output.
pub struct ClosureRung<I, O, const K: usize> {
    name: String,
    closure: Closure<K>,
    key: KeyFn<I>,
    wrap: WrapFn<I, O, K>,
    cost: f64,
}

impl<I, O, const K: usize> ClosureRung<I, O, K> {
    /// `cost` is the per-evaluation cost in the ladder's units.
    pub fn new(
        name: impl Into<String>,
        closure: Closure<K>,
        cost: f64,
        key: impl Fn(&I) -> Option<f64> + Send + Sync + 'static,
        wrap: impl Fn(&I, [f64; K]) -> O + Send + Sync + 'static,
    ) -> Self {
        Self {
            name: name.into(),
            closure,
            key: Box::new(key),
            wrap: Box::new(wrap),
            cost,
        }
    }

    pub fn closure(&self) -> &Closure<K> {
        &self.closure
    }
}

impl<I, O, const K: usize> Rung<I, O> for ClosureRung<I, O, K> {
    fn name(&self) -> &str {
        &self.name
    }
    fn cost(&self, _input: &I) -> f64 {
        self.cost
    }
    fn assess(&self, input: &I) -> Assessment {
        match (self.key)(input).and_then(|x| self.closure.claim_at(x)) {
            Some(e) => Assessment::Claim(Claim::validated(e)),
            None => Assessment::OutOfDomain,
        }
    }
    fn eval(&self, input: &I) -> O {
        let x = (self.key)(input).expect("closure evaluated outside its family");
        let (v, _) = self.closure.eval(x).expect("closure evaluated outside its domain");
        (self.wrap)(input, v)
    }
}
