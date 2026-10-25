//! The diffusive → ballistic ladder of the grey phonon slab
//! ([`crate::bte1d`]) on the emergence engine.
//!
//! The quantity is the effective-conductivity ratio `κ_eff/κ₀` of a slab
//! between diffuse walls as a function of the phonon Knudsen number
//! `Kn = Λ/L`. The metric is the absolute error of the ratio.
//!
//! | rung | cost | claim |
//! |---|---|---|
//! | Fourier, `1` | ~1 ns | `abs(1 − slip) + slip claim` (asymptotic) |
//! | Fourier + Milne jump, `1/(1 + 2ζ Kn)` | ~2 ns | `0.05 E₂(1/Kn)` (asymptotic) |
//! | ballistic, `(3/4)/Kn` | ~1 ns | `(1 + 0.1 ln(1 + Kn))/Kn²` (asymptotic) |
//! | closure table in `ln Kn` | ~100 ns | held-out + reference error (validated) |
//! | integral equation, Richardson, refinement `s = 1, 2, 4` | `~5 µs m²`, `m = s max(40, 20/Kn)` | `0.1 abs(r_2m − r_m)` a posteriori |
//!
//! **Why the jump rung is so good.** Away from the walls the linear
//! profile `i₀ ∝ x` solves the BTE exactly, and each wall's boundary layer is
//! the half-space Milne problem, whose extrapolation length `ζ = 0.710446 Λ`
//! (Case & Zweifel 1967) is exact. The only error is the interaction of the
//! two layers through the slab, which decays like the transmission kernel
//! `E₂(L/Λ) = E₂(1/Kn)`. The constant `0.05` is about 4× the measured
//! prefactor (0.009–0.013 for `Kn = 0.3 … 10`, `analysis/emergence_engine.md`).
//!
//! **Ballistic.** The measured error is `c/Kn²` with `c = 0.63, 0.70, 0.73`
//! at `Kn = 10, 30, 100`; the claim allows a slow logarithmic growth of `c`.
//!
//! **The integral equation** (the exact solution, discretised) converges at
//! order `p ≈ 1.77` (measured 3.4× per doubling; the `E₁` kernel is
//! log-singular). Two grids and Richardson extrapolation with `2^p − 1 = 2.41`
//! agree with the next-finer extrapolation to 0.5–0.6 % of `abs(r_2m − r_m)`,
//! so the a posteriori claim `0.1 abs(r_2m − r_m)` has a ~20× margin.

use emergence::{Assessment, Audit, Axis, Claim, Closure, ClosureOptions, ClosureRung, Ladder, Rung};

use crate::bte1d::solve_gray_bte_integral;
use crate::linalg::expn;

/// Milne extrapolation length in mean free paths (Case & Zweifel 1967).
pub const MILNE_ZETA: f64 = 0.710446;

/// `2^p − 1` for the measured convergence order `p ≈ 1.77` of the
/// integral-equation solver.
const RICHARDSON: f64 = 2.41;

/// The ladder's output: the conductivity ratio, plus the a posteriori
/// estimate when the integral equation produced it (zero otherwise).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct KappaRatio {
    pub ratio: f64,
    pub estimate: f64,
}

pub fn metric(a: &KappaRatio, b: &KappaRatio) -> f64 {
    (a.ratio - b.ratio).abs()
}

/// Fourier with the Milne temperature jump at both walls.
pub fn slip(kn: f64) -> f64 {
    1.0 / (1.0 + 2.0 * MILNE_ZETA * kn)
}

/// Claimed error of [`slip`] (see module docs).
pub fn slip_claim(kn: f64) -> f64 {
    0.05 * expn(2, 1.0 / kn)
}

/// Claimed error of the ballistic limit `(3/4)/Kn` (see module docs).
pub fn ballistic_claim(kn: f64) -> f64 {
    (1.0 + 0.1 * (1.0 + kn).ln()) / (kn * kn)
}

/// Base grid of the integral-equation rung: resolve the `O(Kn)` wall layer.
pub fn ie_cells(kn: f64) -> usize {
    ((20.0 / kn).ceil() as usize).max(40)
}

/// Measured cost of the two integral-equation solves (ns): assembly
/// `~5 µs` per matrix entry plus `m³/3` LU, for `m` and `2m`
/// (i7-13620H; estimate elsewhere).
pub fn ie_cost(kn: f64) -> f64 {
    ie_cost_refined(kn, 1)
}

/// [`ie_cost`] on a grid refined `s` times.
pub fn ie_cost_refined(kn: f64, s: usize) -> f64 {
    let m = (s * ie_cells(kn)) as f64;
    5000.0 * (m * m + 4.0 * m * m) + (m.powi(3) + 8.0 * m.powi(3)) / 3.0
}

/// The integral equation on `m` and `2m` cells, Richardson-extrapolated.
pub fn integral_equation(kn: f64) -> KappaRatio {
    integral_equation_refined(kn, 1)
}

/// [`integral_equation`] on grids refined `s` times.
pub fn integral_equation_refined(kn: f64, s: usize) -> KappaRatio {
    let m = s * ie_cells(kn);
    let r1 = solve_gray_bte_integral(kn, m).kappa_ratio;
    let r2 = solve_gray_bte_integral(kn, 2 * m).kappa_ratio;
    KappaRatio {
        ratio: r2 + (r2 - r1) / RICHARDSON,
        estimate: 0.1 * (r2 - r1).abs() + 1e-13,
    }
}

fn analytic(
    name: &'static str,
    cost: f64,
    value: fn(f64) -> f64,
    claim: impl Fn(f64) -> Option<f64> + Send + Sync + 'static,
) -> emergence::FnRung<f64, KappaRatio> {
    emergence::FnRung::new(
        name,
        move |_| cost,
        move |&kn: &f64| match claim(kn) {
            Some(e) => Assessment::Claim(Claim::asymptotic(e)),
            None => Assessment::OutOfDomain,
        },
        move |&kn: &f64| KappaRatio {
            ratio: value(kn),
            estimate: 0.0,
        },
    )
}

/// The exact rung at refinement `s` (grid `s · max(40, 20/Kn)`). Several
/// refinements on one ladder let the engine escalate only when the cheaper
/// grid's a posteriori certificate fails.
pub struct IntegralEquation {
    pub refine: usize,
    name: String,
}

impl IntegralEquation {
    pub fn new(refine: usize) -> Self {
        Self {
            refine,
            name: format!("integral-equation ×{refine}"),
        }
    }
}

impl Rung<f64, KappaRatio> for IntegralEquation {
    fn name(&self) -> &str {
        &self.name
    }
    fn cost(&self, &kn: &f64) -> f64 {
        ie_cost_refined(kn, self.refine)
    }
    fn assess(&self, _: &f64) -> Assessment {
        Assessment::Deferred
    }
    fn eval(&self, &kn: &f64) -> KappaRatio {
        integral_equation_refined(kn, self.refine)
    }
    fn certify(&self, _: &f64, out: &KappaRatio) -> Option<Claim> {
        Some(Claim::asymptotic(out.estimate))
    }
}

/// Fold the integral equation into a closure over `ln Kn` on
/// `[kn_lo, kn_hi]`, interpolation tolerance `tol`; each patch's claim adds
/// the largest a posteriori error of its samples.
pub fn closure_rung(kn_lo: f64, kn_hi: f64, tol: f64) -> ClosureRung<f64, KappaRatio, 1> {
    let closure = Closure::build_inexact(
        |kn| {
            let r = integral_equation(kn);
            ([r.ratio], r.estimate)
        },
        kn_lo,
        kn_hi,
        Axis::Log,
        tol,
        ClosureOptions {
            degree: 6,
            ..ClosureOptions::default()
        },
        |_, _| f64::INFINITY,
    );
    ClosureRung::new(
        "closure ln Kn",
        closure,
        100.0,
        |&kn: &f64| Some(kn),
        |_, v| KappaRatio {
            ratio: v[0],
            estimate: 0.0,
        },
    )
}

/// The ladder: Fourier, jump, ballistic, an optional closure, and the
/// integral equation at refinements 1, 2 and 4 (the base grid is the audit
/// reference).
pub fn ladder(closure: Option<ClosureRung<f64, KappaRatio, 1>>, audit_every: u64) -> Ladder<f64, KappaRatio> {
    let mut l = Ladder::new("grey-slab κ_eff/κ₀", metric);
    l.push(analytic("fourier", 1.0, |_| 1.0, |kn| Some((1.0 - slip(kn)) + slip_claim(kn))));
    l.push(analytic("fourier+milne-jump", 2.0, slip, |kn| Some(slip_claim(kn))));
    l.push(analytic("ballistic", 1.0, |kn| 0.75 / kn, |kn| (kn >= 1.0).then(|| ballistic_claim(kn))));
    if let Some(c) = closure {
        l.push(c);
    }
    let r = l.push(IntegralEquation::new(1));
    l.push(IntegralEquation::new(2));
    l.push(IntegralEquation::new(4));
    l.reference(r).audit(if audit_every == 0 { Audit::Off } else { Audit::Every(audit_every) })
}

/// One phonon mode of the kinetic quadrature: its contribution to the bulk
/// RTA conductivity (W/(m K)) and its mean free path `v τ` (m).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Mode {
    pub kappa: f64,
    pub mfp: f64,
}

/// The bulk mode spectrum on the same Gauss–Legendre quadrature as
/// [`crate::kinetic::kappa_rta`] (so `Σ kappa` equals it), with the
/// boundary length of `p` ignored (bulk lifetimes).
pub fn mode_spectrum(branches: &[crate::phonon::Branch], p: &crate::scattering::ScatteringParams, t: f64, nq: usize) -> Vec<Mode> {
    use crate::phonon::{mode_heat_capacity, HBAR, KB};
    let bulk = crate::scattering::ScatteringParams {
        length: f64::INFINITY,
        ..*p
    };
    let (xg, wg) = crate::linalg::gauss_legendre_01(nq);
    let mut out = Vec::with_capacity(branches.len() * nq);
    for b in branches {
        let qmax = b.q_max();
        for i in 0..nq {
            let q = xg[i] * qmax;
            let omega = b.omega(q);
            let vg = b.group_velocity(q);
            let c = mode_heat_capacity(HBAR * omega / (KB * t));
            let g = q * q / (2.0 * std::f64::consts::PI * std::f64::consts::PI);
            let rate = bulk.total(b, omega, t);
            let tau = if rate > 0.0 { 1.0 / rate } else { f64::INFINITY };
            out.push(Mode {
                kappa: b.degeneracy * wg[i] * qmax * c * vg * vg * g * tau / 3.0,
                mfp: vg * tau,
            });
        }
    }
    out
}

/// A certified cross-plane film conductivity.
#[derive(Clone, Debug)]
pub struct FilmKappa {
    pub kappa: f64,
    pub kappa_bulk: f64,
    /// Claim on `kappa` (W/(m K)): `Σ_i κ_i e_i` over the per-mode claims.
    pub claim: Claim,
    pub cost: f64,
    /// Modes solved by each rung.
    pub histogram: Vec<usize>,
}

/// Cross-plane conductivity of a film of thickness `l` by per-mode grey
/// suppression, `κ_film = Σ_i κ_i S(Λ_i/l)` with `S` the slab ratio of this
/// ladder, and a goal-oriented budget: `req.tol` (W/(m K)) bounds the error
/// of `κ_film`, so modes carrying little heat get cheap rungs.
///
/// Model assumption: each mode is suppressed as an independent grey slab
/// (the standard suppression-function approximation, e.g. Minnich, PRL 109,
/// 205901 (2012)); it neglects energy exchange between modes inside the
/// film. The claim covers the ladder's error within that model only.
pub fn film_kappa(ladder: &Ladder<f64, KappaRatio>, modes: &[Mode], l: f64, req: emergence::Request) -> FilmKappa {
    let kns: Vec<f64> = modes.iter().map(|m| m.mfp / l).collect();
    let weights: Vec<f64> = modes.iter().map(|m| m.kappa).collect();
    let f = ladder.solve_functional(&kns, &weights, req);
    FilmKappa {
        kappa: f.sum(&weights, |r| r.ratio),
        kappa_bulk: weights.iter().sum(),
        claim: f.claim,
        cost: f.cost,
        histogram: f.histogram(ladder.len()),
    }
}
