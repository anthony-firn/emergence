//! Surface tension in the compressible five-equation model.
//!
//!   capillary_compressible laplace N        static 2D droplet: Δp = σ/R
//!   capillary_compressible wave N ratio     planar capillary wave dispersion
//!
//! Wave: liquid (bottom, depth π) under gas (top, depth π), interface
//! perturbed by a₀ cos x, slip walls top/bottom, periodic in x:
//!   ω² = σk³ / (ρ_l coth(k h_l) + ρ_g coth(k h_g)).
//! Low Mach number by construction (wave speed ~1, sound speeds ≥ 14).
use compressible::capillary::SurfaceTension;
use compressible::eos::Nasg;
use compressible::multiphase::FiveEq;
use compressible::{Bc, Mesh, Solver};

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let case = a.get(1).cloned().unwrap_or("laplace".into());
    let n: usize = a.get(2).map(|v| v.parse().unwrap()).unwrap_or(64);
    let liquid = Nasg::stiffened(4.0, 50.0);
    let gas = Nasg::ideal(1.4);
    if case == "laplace" {
        let (rho_l, rho_g, p0, sigma, r0) = (1.0, 0.1, 1.0, 0.05, 0.2);
        let sys = FiveEq { eos: [liquid, gas], kapila: true };
        let mesh = Mesh::new([n, n, 1], [0.0; 3], [1.0, 1.0, 1.0]);
        let dx = 1.0 / n as f64;
        let mut s = Solver::new(sys, mesh, |x| {
            let r = ((x[0] - 0.5).powi(2) + (x[1] - 0.5).powi(2)).sqrt();
            let al = 0.5 * (1.0 - ((r - r0) / (1.0 * dx)).tanh());
            let al = al.clamp(1e-8, 1.0 - 1e-8);
            let p = if r < r0 { p0 + sigma / r0 } else { p0 };
            [al * rho_l, (1.0 - al) * rho_g, 0.0, 0.0, 0.0, p, al]
        });
        s.surface_tension = Some(SurfaceTension { sigma, alpha_index: 6, smoothing: 2 });
        s.cfl = 0.4;
        for t_end in [0.5, 1.0, 2.0] {
            s.run_to(t_end);
            let (mut pin, mut nin, mut pout, mut nout, mut umax) = (0.0, 0.0, 0.0, 0.0, 0.0f64);
            for i in 0..n {
                for j in 0..n {
                    let x = s.mesh.x(i, j, 0);
                    let r = ((x[0] - 0.5).powi(2) + (x[1] - 0.5).powi(2)).sqrt();
                    let w = s.prim(i, j, 0);
                    if r < 0.5 * r0 {
                        pin += w[5];
                        nin += 1.0;
                    }
                    if r > 1.5 * r0 {
                        pout += w[5];
                        nout += 1.0;
                    }
                    umax = umax.max((w[2] * w[2] + w[3] * w[3]).sqrt());
                }
            }
            let dp = pin / nin - pout / nout;
            println!(
                "laplace N={n} t={t_end}: Δp {dp:.5} vs σ/R {:.5} ({:+.2}%), parasitic |u|max {umax:.2e} (capillary velocity σ/μ-free scale √(σ/(ρR)) = {:.3})",
                sigma / r0,
                100.0 * (dp / (sigma / r0) - 1.0),
                (sigma / (rho_l * r0)).sqrt()
            );
        }
        return;
    }
    let ratio: f64 = a.get(3).map(|v| v.parse().unwrap()).unwrap_or(1.0);
    let pi = std::f64::consts::PI;
    // ambient pressure high enough that the gas sound speed is ≥ 10 (the
    // dispersion relation is incompressible; at p0 = 1 the ratio-1 gas has
    // c = 1.2 and the wave runs at Mach 0.6)
    let (rho_l, rho_g, sigma, a0) = (1.0, 1.0 / ratio, 1.0, 0.05);
    let p0 = (100.0 * rho_g / 1.4f64).max(1.0);
    let (h, k) = (pi, 1.0f64);
    let omega = (sigma * k.powi(3) / (rho_l / (k * h).tanh() + rho_g / (k * h).tanh())).sqrt();
    let sys = FiveEq { eos: [liquid, gas], kapila: true };
    let mesh = Mesh::new([n, n, 1], [0.0; 3], [2.0 * pi, 2.0 * pi, 1.0]);
    let dx = 2.0 * pi / n as f64;
    let mut s = Solver::new(sys, mesh, |x| {
        let eta = pi + a0 * x[0].cos();
        let al = (0.5 * (1.0 - ((x[1] - eta) / dx).tanh())).clamp(1e-8, 1.0 - 1e-8);
        [al * rho_l, (1.0 - al) * rho_g, 0.0, 0.0, 0.0, p0, al]
    });
    s.bc[1] = [Bc::Reflect; 2];
    s.surface_tension = Some(SurfaceTension { sigma, alpha_index: 6, smoothing: 2 });
    s.cfl = 0.4;
    let height = |s: &Solver<7, FiveEq>, i: usize| -> f64 {
        // α = 0.5 crossing in column i
        for j in 0..n - 1 {
            let (a1, a2) = (s.prim(i, j, 0)[6], s.prim(i, j + 1, 0)[6]);
            if (a1 - 0.5) * (a2 - 0.5) <= 0.0 && a1 != a2 {
                return s.mesh.x(i, j, 0)[1] + (a1 - 0.5) / (a1 - a2) * dx;
            }
        }
        f64::NAN
    };
    let amp = |s: &Solver<7, FiveEq>| 0.5 * (height(s, 0) - height(s, n / 2));
    let period = 2.0 * pi / omega;
    let t_end = 2.5 * period;
    let mut samples = vec![(0.0, amp(&s))];
    let t0 = std::time::Instant::now();
    while s.t < t_end {
        let dt = s.dt_cfl().min(t_end - s.t);
        s.step(dt);
        samples.push((s.t, amp(&s)));
        if std::env::var("TRACE").is_ok() && samples.len() % 200 == 0 {
            let umax = (0..n).flat_map(|i| (0..n).map(move |j| (i, j))).map(|(i, j)| { let w = s.prim(i, j, 0); (w[2] * w[2] + w[3] * w[3]).sqrt() }).fold(0.0, f64::max);
            println!("  t {:.3} amp {:+.5} |u|max {umax:.3e}", s.t, samples.last().unwrap().1);
        }
    }
    // zero crossings with hysteresis: a crossing counts only after the
    // amplitude has passed ±0.2 a₀ on the other side (jitter near zero
    // otherwise produces spurious crossings)
    let mut zc = Vec::new();
    let mut armed_sign = if samples[0].1 > 0.0 { 1.0 } else { -1.0 };
    let mut last = samples[0];
    for &(t, a) in samples.iter().skip(1) {
        if a * armed_sign < 0.0 && last.1 * armed_sign >= 0.0 {
            zc.push(last.0 + (t - last.0) * last.1 / (last.1 - a));
            armed_sign = 0.0; // wait for the threshold on the new side
        }
        if armed_sign == 0.0 && a.abs() > 0.2 * a0 {
            armed_sign = a.signum();
        }
        last = (t, a);
    }
    if std::env::var("TRACE").is_ok() {
        println!("  crossings: {:?}", zc.iter().map(|t| (t * 1000.0).round() / 1000.0).collect::<Vec<_>>());
        let n_nan = samples.iter().filter(|x| !x.1.is_finite()).count();
        println!("  NaN amplitude samples: {n_nan} of {}", samples.len());
    }
    // first three half-periods (late in the run parasitic-current noise can
    // exceed the hysteresis band at high density ratio)
    let m = zc.len().min(4);
    let om = if m >= 2 { pi * (m - 1) as f64 / (zc[m - 1] - zc[0]) } else { f64::NAN };
    println!(
        "capillary wave N={n} ratio={ratio}: ω {om:.4} vs theory {omega:.4} ({:+.2}%), amplitude {:.4} → {:.4}, {} steps, {:.0}s",
        100.0 * (om / omega - 1.0),
        samples[0].1,
        samples.last().unwrap().1,
        samples.len(),
        t0.elapsed().as_secs_f64()
    );
}
