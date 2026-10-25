//! Two-phase validation against linear theory.
//!
//!   capillary_probe cap N ratio      capillary wave, zero gravity:
//!                                    ω² = σk³/(ρ₁ + ρ₂)
//!   capillary_probe rt  N ratio      Rayleigh–Taylor (heavy over light):
//!                                    n² = A g k − σk³/(ρ₁ + ρ₂)
//! ratio = ρ₁/ρ₂ (857 for water / nitrogen gas). A slab of fluid 1 occupies
//! z ∈ (π/2, 3π/2) with both interfaces perturbed by a₀ cos(kx), k = 2.
use spectral_core::twophase::{TwoPhaseParams, TwoPhaseSolver};
use spectral_core::{init, Grid};

const PI: f64 = std::f64::consts::PI;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let case = a[1].clone();
    let n: usize = a[2].parse().unwrap();
    let ratio: f64 = a[3].parse().unwrap();
    let g = Grid::new(n);
    let eta = std::env::var("ETA").ok().map(|v| v.parse::<f64>().unwrap()).unwrap_or(1.2) * g.dx();
    let (k, a0) = (2.0, 0.03);
    let rho = [1.0, 1.0 / ratio];
    let sigma = if case == "rt" { 0.02 } else { 1.0 };
    let gravity = if case == "rt" { [0.0, 0.0, -1.0] } else { [0.0; 3] };
    // physical kinematic-viscosity ratio for water / nitrogen gas (ν_g ≈ 15 ν_w);
    // small enough that the inviscid formulas apply (env MU1 overrides)
    let mu1: f64 = std::env::var("MU1").ok().map(|v| v.parse().unwrap()).unwrap_or(2e-4);
    let nu_ratio = if ratio > 1.0 { 15.0 } else { 1.0 };
    let prm = TwoPhaseParams {
        rho,
        mu: [mu1, mu1 * nu_ratio / ratio],
        sigma,
        gravity,
        eta,
        mobility: std::env::var("MOB").ok().map(|v| v.parse().unwrap()).unwrap_or(2e-4),
        stab: 2.0,
    };
    let w = 2f64.sqrt() * eta;
    let (zlo, zhi) = (0.5 * PI, 1.5 * PI);
    let phi = init::from_fn(&g, |x, _, z| {
        let zeta = a0 * (k * x).cos();
        let v = ((z - zlo - zeta) / w).tanh() - ((z - zhi - zeta) / w).tanh() - 1.0;
        [v, 0.0, 0.0]
    })[0]
        .clone();
    let dt_scale: f64 = std::env::var("DTS").ok().map(|v| v.parse().unwrap()).unwrap_or(1.0);
    let mut s = TwoPhaseSolver::new(g.clone(), prm, init::from_fn(&g, |_, _, _| [0.0; 3]), phi);
    let rsum = rho[0] + rho[1];
    let atwood = (rho[0] - rho[1]) / rsum;
    let (dt, t_end, theory) = if case == "rt" {
        // the lower interface (heavy above light) is unstable
        let n2 = atwood * 1.0 * k - sigma * k.powi(3) / rsum;
        (2e-3, 3.0 / n2.sqrt(), n2.sqrt())
    } else {
        let om = (sigma * k.powi(3) / rsum).sqrt();
        let dt_cap = 0.4 * (0.5 * rsum * g.dx().powi(3) / (2.0 * PI * sigma)).sqrt();
        (dt_cap.min(0.005), 2.5 * 2.0 * PI / om, om)
    };
    let dt = dt * dt_scale;
    let mass0 = s.phase_mean();
    let col_a = 0usize; // x = 0
    let col_b = n / 4; // x = π/2, where cos(2x) = −1
    let (zi, sign) = if case == "rt" { (zlo, 1.0) } else { (zhi, 1.0) };
    let amp = |s: &TwoPhaseSolver| sign * 0.5 * (s.interface_height(col_a, 0, zi) - s.interface_height(col_b, 0, zi));
    let mut samples = vec![(0.0, amp(&s))];
    let every = ((t_end / dt) / 120.0).ceil() as usize;
    let mut step = 0;
    let t0 = std::time::Instant::now();
    while s.t < t_end {
        s.step(dt);
        step += 1;
        if step % every == 0 {
            samples.push((s.t, amp(&s)));
        }
    }
    if case == "rt" {
        // The interface starts displaced at rest, so linear theory gives
        // a(t) = a₀ cosh(nt) (a log-slope fit is biased low early on). Fit n
        // by least squares on ln a over the linear range k·a < 0.2.
        let pts: Vec<(f64, f64)> = samples.iter().filter(|(_, a)| a.abs() * k < 0.2).map(|&(t, a)| (t, a.abs().ln())).collect();
        let a_init = samples[0].1.abs();
        let cost = |nn: f64| pts.iter().map(|&(t, la)| (la - (a_init * (nn * t).cosh()).ln()).powi(2)).sum::<f64>();
        let (mut lo, mut hi) = (0.1 * theory, 2.0 * theory);
        for _ in 0..100 {
            let (m1, m2) = (lo + (hi - lo) / 3.0, hi - (hi - lo) / 3.0);
            if cost(m1) < cost(m2) { hi = m2 } else { lo = m1 }
        }
        let slope = 0.5 * (lo + hi);
        println!("RT N={n} ratio={ratio}: growth rate {slope:.4} vs inviscid theory {theory:.4} ({:+.2}%), {} fit points, {:.0}s", 100.0 * (slope / theory - 1.0), pts.len(), t0.elapsed().as_secs_f64());
    } else {
        // period from zero crossings of the amplitude
        let mut zc = Vec::new();
        for w2 in samples.windows(2) {
            let ((t1, a1), (t2, a2)) = (w2[0], w2[1]);
            if a1 * a2 < 0.0 {
                zc.push(t1 + (t2 - t1) * a1 / (a1 - a2));
            }
        }
        let period = 2.0 * (zc[zc.len() - 1] - zc[0]) / (zc.len() - 1) as f64;
        let om = 2.0 * PI / period;
        println!(
            "capillary N={n} ratio={ratio}: omega {om:.4} vs inviscid theory {theory:.4} ({:+.2}%), amplitude {:.4} -> {:.4}, mass drift {:.1e}, {:.0}s",
            100.0 * (om / theory - 1.0),
            samples[0].1,
            samples.last().unwrap().1,
            s.phase_mean() - mass0,
            t0.elapsed().as_secs_f64()
        );
    }
}
