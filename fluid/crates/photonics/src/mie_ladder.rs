//! The Mie ladder on the emergence engine.
//!
//! Rungs for the efficiencies of a homogeneous sphere, each with its error
//! claim in the ladder metric (largest absolute error in `Q_ext`, `Q_sca` and
//! `g Q_sca`):
//!
//! | rung | domain | claim |
//! |---|---|---|
//! | [`Rayleigh2`] | `max(x, abs(m) x) ≤ 1` | asymptotic: `abs(R2 − R1)` |
//! | closure table ([`closure_rung`]) | one `m`, a range of `x`, every resonance resolved | validated: held-out error |
//! | [`ExactMie`] | all `x` | validated: Wiscombe truncation `2e-8` (`x ≤ 2000`), `1e-6` beyond |
//!
//! `R1` is the classical Rayleigh result; `R2` adds the next order of the
//! small-`x` expansion of `a_1`, `b_1`, `a_2` (Bohren & Huffman 1983 eq. 5.7,
//! complex-conjugated to this crate's `n + iκ` convention). The rung returns
//! `R2` and claims `abs(R2 − R1)`, the embedded estimate of `R1`'s error, which
//! bounds `R2`'s error with a wide margin while `x` is small (local
//! extrapolation, as in embedded Runge–Kutta pairs).
//!
//! Costs are nanoseconds per evaluation measured on an i7-13620H
//! (`analysis/emergence_engine.md`; estimates on other machines).

use std::f64::consts::PI;

use emergence::{Assessment, Audit, Axis, Claim, Closure, ClosureOptions, ClosureRung, Functional, Ladder, Request, Rung};

use crate::complex::C;
use crate::mie::{wiscombe_terms, GammaSizeDistribution, MieScattering};

/// A sphere: size parameter `x = 2π n_medium a / λ` and relative index `m`
/// (`n + iκ`, absorption `κ > 0`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MieInput {
    pub x: f64,
    pub m: C,
}

/// The efficiencies the ladder computes. `gq_sca = g Q_sca` is used instead of
/// `g`, which is ill-defined as `Q_sca → 0` and is not additive over sizes.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MieQ {
    pub q_ext: f64,
    pub q_sca: f64,
    pub gq_sca: f64,
}

impl MieQ {
    pub fn q_abs(&self) -> f64 {
        self.q_ext - self.q_sca
    }
    pub fn g(&self) -> f64 {
        if self.q_sca > 0.0 {
            self.gq_sca / self.q_sca
        } else {
            0.0
        }
    }
    pub fn to_vec(self) -> [f64; 3] {
        [self.q_ext, self.q_sca, self.gq_sca]
    }
    pub fn from_slice(v: [f64; 3]) -> Self {
        Self {
            q_ext: v[0],
            q_sca: v[1],
            gq_sca: v[2],
        }
    }
}

/// The ladder metric: largest absolute difference of the three components.
pub fn metric(a: &MieQ, b: &MieQ) -> f64 {
    (a.q_ext - b.q_ext)
        .abs()
        .max((a.q_sca - b.q_sca).abs())
        .max((a.gq_sca - b.gq_sca).abs())
}

/// The exact series, Wiscombe-truncated.
pub fn exact(x: f64, m: C) -> MieQ {
    let e = MieScattering::new(x, m).efficiencies();
    MieQ {
        q_ext: e.q_ext,
        q_sca: e.q_sca,
        gq_sca: e.g * e.q_sca,
    }
}

/// First-order (classical) Rayleigh efficiencies `R1`.
pub fn rayleigh1(x: f64, m: C) -> MieQ {
    let m2 = m * m;
    let k = (m2 - 1.0) / (m2 + 2.0);
    let q_sca = 8.0 / 3.0 * x.powi(4) * k.abs2();
    MieQ {
        q_ext: 4.0 * x * k.im + q_sca,
        q_sca,
        gq_sca: 0.0,
    }
}

/// Second-order Rayleigh efficiencies `R2` from the small-`x` expansions of
/// `a_1`, `b_1`, `a_2` (Bohren & Huffman 1983 eq. 5.7, conjugated).
pub fn rayleigh2(x: f64, m: C) -> MieQ {
    if x == 0.0 {
        return MieQ::default();
    }
    let m2 = m * m;
    let k = (m2 - 1.0) / (m2 + 2.0);
    let i = C::i();
    // Bohren & Huffman's expansions, then conjugated (checked against the
    // exact a_1, b_1, a_2: relative residual O(x^2)).
    let a1 = -(i * k) * (2.0 * x.powi(3) / 3.0)
        - (i * ((m2 - 2.0) * (m2 - 1.0) / ((m2 + 2.0) * (m2 + 2.0)))) * (2.0 * x.powi(5) / 5.0)
        + k * k * (4.0 * x.powi(6) / 9.0);
    let b1 = -(i * (m2 - 1.0)) * (x.powi(5) / 45.0);
    let a2 = -(i * ((m2 - 1.0) / (m2 * 2.0 + 3.0))) * (x.powi(5) / 15.0);
    let (a1, b1, a2) = (a1.conj(), b1.conj(), a2.conj());
    let x2 = x * x;
    MieQ {
        q_ext: 2.0 / x2 * (3.0 * (a1.re + b1.re) + 5.0 * a2.re),
        q_sca: 2.0 / x2 * (3.0 * (a1.abs2() + b1.abs2()) + 5.0 * a2.abs2()),
        gq_sca: 4.0 / x2 * 1.5 * ((a1 * a2.conj()).re + (a1 * b1.conj()).re),
    }
}

/// Measured cost of the exact series (ns): `250 + 28 n_terms`
/// (i7-13620H, `x = 0.1 … 1000`; estimate elsewhere).
pub fn exact_cost(x: f64) -> f64 {
    250.0 + 28.0 * wiscombe_terms(x) as f64
}

/// Second-order Rayleigh with the embedded error estimate.
pub struct Rayleigh2;

impl Rung<MieInput, MieQ> for Rayleigh2 {
    fn name(&self) -> &str {
        "rayleigh-2"
    }
    fn cost(&self, _: &MieInput) -> f64 {
        // ~2× the R2 evaluation for the claim (estimate).
        60.0
    }
    fn assess(&self, i: &MieInput) -> Assessment {
        if i.x.max(i.m.abs() * i.x) > 1.0 {
            return Assessment::OutOfDomain;
        }
        Assessment::Claim(Claim::asymptotic(metric(&rayleigh2(i.x, i.m), &rayleigh1(i.x, i.m))))
    }
    fn eval(&self, i: &MieInput) -> MieQ {
        rayleigh2(i.x, i.m)
    }
}

/// The exact series.
pub struct ExactMie;

impl Rung<MieInput, MieQ> for ExactMie {
    fn name(&self) -> &str {
        "exact-mie"
    }
    fn cost(&self, i: &MieInput) -> f64 {
        exact_cost(i.x)
    }
    fn assess(&self, i: &MieInput) -> Assessment {
        // Truncation measured against n_stop + 30 + n_stop/5 terms over
        // x ≤ 2000 and seven indices (max 7.3e-9); Wiscombe (1980) beyond.
        Assessment::Claim(Claim::validated(if i.x <= 2000.0 { 2e-8 } else { 1e-6 }))
    }
    fn eval(&self, i: &MieInput) -> MieQ {
        exact(i.x, i.m)
    }
}

/// One multipole resonance of a sphere.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Resonance {
    /// Size parameter of the peak.
    pub x: f64,
    /// Full width at half maximum, in `x`.
    pub width: f64,
    /// Multipole order `n` (1-based) and type (`true`: electric `a_n`).
    pub order: usize,
    pub electric: bool,
    /// The mode's contribution to `Q_ext` at the peak,
    /// `(2/x²)(2n+1) Re c_n`.
    pub peak: f64,
}

/// All multipole resonances of the (Wiscombe-truncated) series with peak
/// contribution above `min_peak`, for `x_lo ≤ x ≤ x_hi`.
///
/// Each coefficient obeys `1/c_n = 1 + D_n + i C_n(x)` with `C_n` real and
/// smooth between poles and `D_n ≥ 0` the absorption, so `abs(c_n)² =
/// 1/((1 + D_n)² + C_n²)` is a Lorentzian in `C_n`. A resonance is a zero of
/// `C_n = Im(1/c_n)`: a sign scan cannot miss it however narrow it is (unlike
/// sampling `Q` itself). Each sign change is refined by Illinois regula falsi;
/// a pole of `C_n` (`c_n → 0`, lossless spheres) is recognised by `abs(C_n)`
/// growing past its bracket values and dropped, and anti-resonances of
/// absorbing spheres fall below `min_peak`. The full width at half
/// maximum is `Γ = 2 (1 + D_n)/abs(C_n')`. `step` must be smaller than the
/// distance between a zero and the neighbouring pole of the same `C_n`
/// (`~1` in `x`; the tests use `0.01` and check the count against `0.005`).
pub fn resonances(m: C, x_lo: f64, x_hi: f64, step: f64, min_peak: f64) -> Vec<Resonance> {
    use rayon::prelude::*;
    let steps = ((x_hi - x_lo) / step).ceil() as usize;
    let xs: Vec<f64> = (0..=steps).map(|i| x_lo + (x_hi - x_lo) * i as f64 / steps as f64).collect();
    // Im(1/c_n) for c = a (electric) or b, up to the series length at x.
    let coef = |x: f64, n: usize, electric: bool| -> Option<C> {
        let s = MieScattering::with_terms(x, m, n + 1);
        let c = if electric { s.a[n] } else { s.b[n] };
        c.is_finite().then_some(c)
    };
    let chunk = 256;
    let mut found: Vec<Resonance> = (0..steps.div_ceil(chunk))
        .into_par_iter()
        .flat_map_iter(|k| {
            let i0 = k * chunk;
            let i1 = ((k + 1) * chunk).min(steps);
            let mut out = Vec::new();
            let mut prev = MieScattering::new(xs[i0], m);
            for i in i0..i1 {
                let next = MieScattering::new(xs[i + 1], m);
                let n_common = prev.a.len().min(next.a.len());
                for n in 0..n_common {
                    for electric in [true, false] {
                        let (c0, c1) = if electric { (prev.a[n], next.a[n]) } else { (prev.b[n], next.b[n]) };
                        let (f0, f1) = (c0.recip().im, c1.recip().im);
                        if !(f0 * f1 < 0.0) {
                            continue;
                        }
                        // Illinois regula falsi: C_n is linear near a zero,
                        // so this converges in a few steps; a pole of C_n
                        // shows up as |C| exceeding the bracket values.
                        let (mut lo, mut hi, mut flo, mut fhi) = (xs[i], xs[i + 1], f0, f1);
                        let bound = 1e3 * f0.abs().max(f1.abs());
                        let mut side = 0i8;
                        let mut pole = false;
                        for _ in 0..100 {
                            let mid = (lo * fhi - hi * flo) / (fhi - flo);
                            let mid = if mid > lo && mid < hi { mid } else { 0.5 * (lo + hi) };
                            let Some(c) = coef(mid, n, electric) else {
                                pole = true;
                                break;
                            };
                            let fm = c.recip().im;
                            if fm.abs() > bound {
                                pole = true;
                                break;
                            }
                            if fm == 0.0 || hi - lo <= 4.0 * f64::EPSILON * hi {
                                lo = mid;
                                hi = mid;
                                break;
                            }
                            if fm * flo > 0.0 {
                                lo = mid;
                                flo = fm;
                                if side == -1 {
                                    fhi *= 0.5;
                                }
                                side = -1;
                            } else {
                                hi = mid;
                                fhi = fm;
                                if side == 1 {
                                    flo *= 0.5;
                                }
                                side = 1;
                            }
                            if (hi - lo) <= 1e-13 * hi {
                                break;
                            }
                        }
                        if pole {
                            continue;
                        }
                        let x = 0.5 * (lo + hi);
                        let Some(c) = coef(x, n, electric) else { continue };
                        let inv = c.recip();
                        if !(inv.abs() < 1e6) {
                            continue; // a pole of C_n, not a resonance
                        }
                        let order = n + 1;
                        let peak = 2.0 / (x * x) * (2 * order + 1) as f64 * c.re;
                        if peak < min_peak {
                            continue;
                        }
                        let d = 1e-7 * x.max(1.0);
                        let slope = match (coef(x + d, n, electric), coef(x - d, n, electric)) {
                            (Some(p), Some(q)) => (p.recip().im - q.recip().im) / (2.0 * d),
                            _ => continue,
                        };
                        let width = 2.0 * inv.re / slope.abs();
                        out.push(Resonance {
                            x,
                            width,
                            order,
                            electric,
                            peak,
                        });
                    }
                }
                prev = next;
            }
            out
        })
        .collect();
    found.sort_by(|a, b| a.x.total_cmp(&b.x));
    found
}

/// Measured cost of a closure evaluation (ns, degree-16 patch, bucket index
/// and Clenshaw; 55 ns serial on an i7-13620H, plus the lookup in `assess`).
pub const CLOSURE_COST: f64 = 70.0;

/// Fold the exact series into a closure table for one index `m` on
/// `[x_lo, x_hi]`, to absolute tolerance `tol`.
///
/// Every resonance with a peak above `tol/4` ([`resonances`]) must be sampled
/// at a spacing of at most half its width; a patch around a resonance too
/// narrow for `max_patches` stays out of domain, so the ladder unfolds to the
/// series exactly there.
pub fn closure_rung(m: C, x_lo: f64, x_hi: f64, tol: f64, max_patches: usize) -> ClosureRung<MieInput, MieQ, 3> {
    let res = resonances(m, x_lo, x_hi, 0.01, 0.25 * tol);
    let opt = ClosureOptions {
        degree: 8,
        min_width: 1e-9 * (x_hi - x_lo),
        max_patches,
    };
    let closure = Closure::build_resolving(
        |x| exact(x, m).to_vec(),
        x_lo,
        x_hi,
        Axis::Linear,
        tol,
        opt,
        |a, b| narrowest(&res, a, b),
    );
    ClosureRung::new(
        format!("closure m={:.4}+{:.2e}i", m.re, m.im),
        closure,
        CLOSURE_COST,
        move |i: &MieInput| (i.m == m).then_some(i.x),
        |_, v| MieQ::from_slice(v),
    )
}

/// Half the narrowest width among resonances whose peak lies within
/// `[a − 10Γ, b + 10Γ]` (the Lorentzian tails count too).
fn narrowest(res: &[Resonance], a: f64, b: f64) -> f64 {
    // Resonances are sorted by x; widths are at most ~1, so scan a window.
    let start = res.partition_point(|r| r.x < a - 10.0);
    res[start..]
        .iter()
        .take_while(|r| r.x <= b + 10.0)
        .filter(|r| r.x + 10.0 * r.width >= a && r.x - 10.0 * r.width <= b)
        .map(|r| 0.5 * r.width)
        .fold(f64::INFINITY, f64::min)
}

/// The ladder: Rayleigh-2, any closures given, and the exact series as the
/// reference, audited on every `audit_every`-th accepted result per rung.
pub fn ladder(closures: Vec<ClosureRung<MieInput, MieQ, 3>>, audit_every: u64) -> Ladder<MieInput, MieQ> {
    let mut l = Ladder::new("mie", metric);
    l.push(Rayleigh2);
    for c in closures {
        l.push(c);
    }
    let r = l.push(ExactMie);
    l.reference(r).audit(if audit_every == 0 { Audit::Off } else { Audit::Every(audit_every) })
}

/// Size-averaged properties with their claimed errors.
#[derive(Clone, Copy, Debug)]
pub struct AveragedQ {
    /// Area-weighted mean efficiencies `<Q_ext>`, `<Q_sca>`, `<g Q_sca>`.
    pub q: MieQ,
    /// Mean geometric cross section `<π r²>` (m²).
    pub mean_area: f64,
    /// Claim on each area-weighted mean efficiency (model error only).
    pub claim: Claim,
    /// Engine cost (ns-equivalent, excluding audits).
    pub cost: f64,
}

impl AveragedQ {
    /// Asymmetry parameter and its claimed error:
    /// `abs(δg) ≤ e (1 + abs(g)) / (<Q_sca> − e)`.
    pub fn g(&self) -> (f64, f64) {
        let g = self.q.g();
        let e = self.claim.error;
        (g, e * (1.0 + g.abs()) / (self.q.q_sca - e).max(0.0))
    }
}

/// Simpson nodes and weights on `[0, r_max]` (the same rule as
/// [`crate::mie::gamma_averaged`]).
fn simpson_nodes(r_max: f64, n: usize) -> (Vec<f64>, Vec<f64>) {
    let n = n + n % 2;
    let h = r_max / n as f64;
    let r = (0..=n).map(|i| i as f64 * h).collect();
    let w = (0..=n)
        .map(|i| {
            h / 3.0
                * if i == 0 || i == n {
                    1.0
                } else if i % 2 == 1 {
                    4.0
                } else {
                    2.0
                }
        })
        .collect();
    (r, w)
}

/// Gamma-distribution average of the efficiencies through the ladder, with a
/// goal-oriented budget: `req.tol` bounds the error of each area-weighted mean
/// efficiency, and each node gets only the accuracy its weight requires.
/// Quadrature error is not included (same rule as `gamma_averaged`).
pub fn gamma_averaged(
    ladder: &Ladder<MieInput, MieQ>,
    dist: &GammaSizeDistribution,
    m: C,
    wavelength: f64,
    n_medium: f64,
    n_quad: usize,
    req: Request,
) -> (AveragedQ, Functional<MieQ>) {
    let r_max = dist.r_eff * (1.0 + 8.0 * dist.v_eff.sqrt());
    let k = 2.0 * PI * n_medium / wavelength;
    let (r, w) = simpson_nodes(r_max, n_quad);
    let area_w: Vec<f64> = r.iter().zip(&w).map(|(&r, &w)| w * PI * r * r * dist.pdf(r)).collect();
    let mean_area: f64 = area_w.iter().sum();
    let weights: Vec<f64> = area_w.iter().map(|a| a / mean_area).collect();
    let inputs: Vec<MieInput> = r.iter().map(|&r| MieInput { x: k * r, m }).collect();
    let f = ladder.solve_functional(&inputs, &weights, req);
    let q = MieQ {
        q_ext: f.sum(&weights, |q| q.q_ext),
        q_sca: f.sum(&weights, |q| q.q_sca),
        gq_sca: f.sum(&weights, |q| q.gq_sca),
    };
    (
        AveragedQ {
            q,
            mean_area,
            claim: f.claim,
            cost: f.cost,
        },
        f,
    )
}
