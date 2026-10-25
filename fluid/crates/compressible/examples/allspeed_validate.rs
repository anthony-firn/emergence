//! Mach-adaptive all-speed scheme vs the explicit compressible scheme.
//!
//!   allspeed_validate gresho N        Gresho vortex at Mach 0.1, 0.01, 0.001
//!                                     (steady incompressible solution; low-Mach accuracy and cost)
//!   allspeed_validate sod N           Sod shock tube vs the exact solution (compressible end)
//!   allspeed_validate pulse N         acoustic pulse, acoustic Courant ≤ 0.5 (acoustics resolved)
//!
//! Env EXPLICIT=0 skips the explicit runs (Gresho at Mach 0.001 takes
//! ~1000× more explicit steps than all-speed ones).
use compressible::allspeed::AllSpeedOpts;
use compressible::eos::Nasg;
use compressible::exact::{sample, star, Side};
use compressible::{Bc, Euler, Mesh, Solver};
use std::time::Instant;

fn gresho(n: usize, mach: f64, allspeed: bool, t_end: f64) -> (f64, f64, usize, f64, f64) {
    let eos = Nasg::ideal(1.4);
    let p0 = 1.0 / (1.4 * mach * mach);
    let vel = |r: f64| if r < 0.2 { 5.0 * r } else if r < 0.4 { 2.0 - 5.0 * r } else { 0.0 };
    let pres = |r: f64| {
        if r < 0.2 {
            p0 + 12.5 * r * r
        } else if r < 0.4 {
            p0 + 12.5 * r * r + 4.0 - 20.0 * r + 4.0 * (5.0 * r).ln()
        } else {
            p0 - 2.0 + 4.0 * 2f64.ln()
        }
    };
    let init = move |x: [f64; 3]| {
        let (dx, dy) = (x[0] - 0.5, x[1] - 0.5);
        let r = (dx * dx + dy * dy).sqrt();
        let v = vel(r);
        let (ux, uy) = if r > 0.0 { (-v * dy / r, v * dx / r) } else { (0.0, 0.0) };
        [1.0, ux, uy, 0.0, pres(r)]
    };
    let mut s = Solver::new(Euler { eos }, Mesh::new([n, n, 1], [0.0; 3], [1.0, 1.0, 1.0]), init);
    s.cfl = 0.4;
    if allspeed {
        s.allspeed = Some(AllSpeedOpts { cfl: 0.4, ..Default::default() });
    }
    let ke = |s: &Solver<5, Euler>| -> f64 {
        let mut e = 0.0;
        for i in 0..n {
            for j in 0..n {
                let w = s.prim(i, j, 0);
                e += 0.5 * w[0] * (w[1] * w[1] + w[2] * w[2]);
            }
        }
        e
    };
    let ke0 = ke(&s);
    let t0 = Instant::now();
    let steps = s.run_to(t_end);
    let el = t0.elapsed().as_secs_f64();
    // L1 velocity error vs the steady solution
    let mut err = 0.0;
    for i in 0..n {
        for j in 0..n {
            let x = s.mesh.x(i, j, 0);
            let w = s.prim(i, j, 0);
            let e = init(x);
            err += ((w[1] - e[1]).powi(2) + (w[2] - e[2]).powi(2)).sqrt();
        }
    }
    let st = s.allspeed_stats;
    let cg = if st.solves > 0 { st.cg_iterations as f64 / st.solves as f64 } else { 0.0 };
    (ke(&s) / ke0, err / (n * n) as f64, steps, el, cg)
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let case = a.get(1).cloned().unwrap_or("gresho".into());
    let n: usize = a.get(2).map(|v| v.parse().unwrap()).unwrap_or(64);
    let explicit = std::env::var("EXPLICIT").map(|v| v != "0").unwrap_or(true);
    let eos = Nasg::ideal(1.4);
    match case.as_str() {
        "gresho" => {
            let t_end = 1.0;
            println!("Gresho vortex {n}², t = {t_end} (0.8 turns of the vortex core): KE(t)/KE(0), L1 |u − u_exact|, steps, wall time");
            for mach in [0.1, 0.01, 0.001] {
                let (k, e, st, el, cg) = gresho(n, mach, true, t_end);
                println!("  Mach {mach:<6} all-speed: KE {k:.6}, L1 err {e:.3e}, {st} steps, {el:.2}s, {cg:.0} CG iterations/solve");
                // explicit at Mach 0.001 needs ~10⁵ steps: only with EXPLICIT=all
                if explicit && (mach > 0.002 || std::env::var("EXPLICIT").as_deref() == Ok("all")) {
                    let (k, e, st, el, _) = gresho(n, mach, false, t_end);
                    println!("  Mach {mach:<6} explicit : KE {k:.6}, L1 err {e:.3e}, {st} steps, {el:.2}s");
                }
            }
        }
        "sod" => {
            let (l, r) = (Side { rho: 1.0, u: 0.0, p: 1.0, eos }, Side { rho: 0.125, u: 0.0, p: 0.1, eos });
            let st = star(&l, &r);
            let t_end = 0.2;
            for allspeed in [false, true] {
                let mut s = Solver::new(Euler { eos }, Mesh::new([n, 1, 1], [0.0; 3], [1.0, 1.0, 1.0]), |x| if x[0] < 0.5 { [1.0, 0.0, 0.0, 0.0, 1.0] } else { [0.125, 0.0, 0.0, 0.0, 0.1] });
                s.bc[0] = [Bc::Outflow; 2];
                s.cfl = 0.4;
                if allspeed {
                    s.allspeed = Some(AllSpeedOpts { cfl: 0.4, ..Default::default() });
                }
                let steps = s.run_to(t_end);
                let mut err = [0.0; 3];
                for i in 0..n {
                    let x = s.mesh.x(i, 0, 0)[0];
                    let e = sample(&l, &r, &st, (x - 0.5) / t_end);
                    let w = s.prim(i, 0, 0);
                    err[0] += (w[0] - e.0).abs() / n as f64;
                    err[1] += (w[1] - e.1).abs() / n as f64;
                    err[2] += (w[4] - e.2).abs() / n as f64;
                }
                let tot = s.totals();
                println!(
                    "Sod {n} {}: L1 error ρ {:.3e}, u {:.3e}, p {:.3e}; {steps} steps, max acoustic Courant {:.2}, mass {:.12}",
                    if allspeed { "all-speed" } else { "explicit " },
                    err[0],
                    err[1],
                    err[2],
                    s.allspeed_stats.max_acoustic_cfl,
                    tot[0]
                );
            }
        }
        "pulse" => {
            // Gaussian pressure pulse in air at rest: two halves at ±c
            let (p0, amp, sig) = (1.0, 1e-3, 0.03);
            let c = (1.4f64 * p0).sqrt();
            let t_end = 0.25 / c;
            let init = move |x: [f64; 3]| {
                let dp = amp * (-((x[0] - 0.5) / sig).powi(2)).exp();
                [1.0 + dp / (c * c), 0.0, 0.0, 0.0, p0 + dp]
            };
            for allspeed in [false, true] {
                let mut s = Solver::new(Euler { eos }, Mesh::new([n, 1, 1], [0.0; 3], [1.0, 1.0, 1.0]), init);
                s.cfl = 0.4;
                if allspeed {
                    // min_gain 0: force the all-speed scheme (resolved acoustics)
                    s.allspeed = Some(AllSpeedOpts { cfl: 0.4, max_acoustic_cfl: 0.5, min_gain: 0.0, ..Default::default() });
                }
                let steps = s.run_to(t_end);
                // linear exact: two half pulses at 0.5 ± c t
                let mut err = 0.0;
                let mut peak: (f64, f64) = (0.0, 0.0);
                for i in 0..n {
                    let x = s.mesh.x(i, 0, 0)[0];
                    let ex = 0.5 * amp * ((-((x - 0.5 - c * t_end) / sig).powi(2)).exp() + (-((x - 0.5 + c * t_end) / sig).powi(2)).exp());
                    let dp = s.prim(i, 0, 0)[4] - p0;
                    err = f64::max(err, (dp - ex).abs());
                    if x > 0.5 && dp > peak.1 {
                        peak = (x, dp);
                    }
                }
                println!(
                    "pulse {n} {}: max |Δp − exact| / amplitude {:.3e}, right peak at {:.4} (exact {:.4}), height {:.4} of {:.4}; {steps} steps",
                    if allspeed { "all-speed" } else { "explicit " },
                    err / amp,
                    peak.0,
                    0.5 + c * t_end,
                    peak.1 / amp,
                    0.5
                );
            }
        }
        "mix1d" => {
            // five-equation model: (a) pure gas in a converging flow, all-speed
            // vs explicit; (b) water slab advected at 10 m/s through air at
            // uniform pressure (Abgrall: p and u must stay uniform); (c) the
            // water | air shock tube (1e9 | 1e5 Pa) vs the exact solution
            use compressible::eos::materials::{AIR, WATER_SG};
            use compressible::multiphase::FiveEq;
            let sys = FiveEq { eos: [WATER_SG, AIR], kapila: true };
            let run = |init: &(dyn Fn([f64; 3]) -> [f64; 7] + Sync), allspeed: bool, bc: Bc<7>, t_end: f64| {
                let mut s = Solver::new(sys, Mesh::new([n, 1, 1], [0.0; 3], [1.0, 1.0, 1.0]), init);
                s.bc[0] = [bc; 2];
                s.cfl = 0.4;
                if allspeed {
                    s.allspeed = Some(AllSpeedOpts { cfl: 0.4, min_gain: std::env::var("GAIN").ok().map(|v| v.parse().unwrap()).unwrap_or(3.0), ..Default::default() });
                }
                let steps = s.run_to(t_end);
                (s, steps)
            };
            let gas = |x: [f64; 3]| [1e-6 * 1000.0, (1.0 - 1e-6) * 1.19, -20.0 * (2.0 * std::f64::consts::PI * x[0]).sin(), 0.0, 0.0, 1e5, 1e-6];
            let (a, na) = run(&gas, false, Bc::Periodic, 2e-3);
            let (b, nb) = run(&gas, true, Bc::Periodic, 2e-3);
            let d = (0..n).map(|i| (a.prim(i, 0, 0)[5] - b.prim(i, 0, 0)[5]).abs()).fold(0.0, f64::max);
            println!("five-eq pure gas, converging flow: max |p_allspeed − p_explicit| {d:.3e} Pa of 1e5 ({na} explicit, {nb} all-speed steps)");
            let slab = |x: [f64; 3]| {
                let al = 1e-6 + (1.0 - 2e-6) * 0.5 * (((x[0] - 0.3) / 0.01).tanh() - ((x[0] - 0.7) / 0.01).tanh());
                [al * 1000.0, (1.0 - al) * 1.19, 10.0, 0.0, 0.0, 1e5, al]
            };
            let t_end = 0.1; // one period
            let (b, nb) = run(&slab, true, Bc::Periodic, t_end);
            let (mut dp, mut du) = (0.0f64, 0.0f64);
            for i in 0..n {
                let w = b.prim(i, 0, 0);
                dp = dp.max((w[5] / 1e5 - 1.0).abs());
                du = du.max((w[2] / 10.0 - 1.0).abs());
            }
            println!("water slab advected one period (all-speed, {nb} steps, max acoustic Courant {:.0}): max |p/p0 − 1| {dp:.2e}, |u/u0 − 1| {du:.2e}", b.allspeed_stats.max_acoustic_cfl);
            let p_hi: f64 = std::env::var("PHI").ok().map(|v| v.parse().unwrap()).unwrap_or(1e9);
            let (l, r) = (Side { rho: 1000.0, u: 0.0, p: p_hi, eos: WATER_SG }, Side { rho: 1.19, u: 0.0, p: 1e5, eos: AIR });
            let st = star(&l, &r);
            let tube = |x: [f64; 3]| {
                let al = if x[0] < 0.7 { 1.0 - 1e-6 } else { 1e-6 };
                [al * 1000.0, (1.0 - al) * 1.19, 0.0, 0.0, 0.0, if x[0] < 0.7 { p_hi } else { 1e5 }, al]
            };
            let t_end = std::env::var("TT").ok().map(|v| v.parse().unwrap()).unwrap_or(2.4e-4);
            for allspeed in [false, true] {
                let (s, ns) = run(&tube, allspeed, Bc::Outflow, t_end);
                let mut err = [0.0; 2];
                for i in 0..n {
                    let x = s.mesh.x(i, 0, 0)[0];
                    let e = sample(&l, &r, &st, (x - 0.7) / t_end);
                    let w = s.prim(i, 0, 0);
                    err[0] += ((w[0] + w[1]) - e.0).abs() / n as f64;
                    err[1] += (w[5] - e.2).abs() / n as f64;
                }
                println!("water|air tube {}: L1 error ρ {:.3e}, p {:.3e}; {ns} steps", if allspeed { "all-speed" } else { "explicit " }, err[0], err[1]);
            }
        }
        "contact" => {
            // Euler: dense slab (ρ 1000) in ρ = 1 at uniform p, u = 10 (contact)
            let rho_s: f64 = std::env::var("RHO").ok().map(|v| v.parse().unwrap()).unwrap_or(1000.0);
            let init = move |x: [f64; 3]| {
                let h = 0.5 * (((x[0] - 0.3) / 0.01).tanh() - ((x[0] - 0.7) / 0.01).tanh());
                [1.0 + (rho_s - 1.0) * h, 10.0, 0.0, 0.0, 1e5]
            };
            let mut s = Solver::new(Euler { eos }, Mesh::new([n, 1, 1], [0.0; 3], [1.0, 1.0, 1.0]), init);
            let cap: f64 = std::env::var("CAP").ok().map(|v| v.parse().unwrap()).unwrap_or(f64::INFINITY);
            s.allspeed = Some(AllSpeedOpts { cfl: 0.4, max_acoustic_cfl: cap, ..Default::default() });
            for st in 0..40 {
                let dt = s.dt_cfl();
                s.step(dt);
                let mut dp = 0.0f64;
                for i in 0..n {
                    let w = s.prim(i, 0, 0);
                    let e = (w[4] / 1e5 - 1.0).abs();
                    dp = if e > dp || e.is_nan() { e } else { dp };
                }
                if st % 10 == 9 {
                    println!("contact step {st}: dt {dt:.3e}, max |p/p0 − 1| {dp:.3e}");
                }
            }
        }
        "slab" => {
            use compressible::eos::materials::{AIR, WATER_SG};
            use compressible::multiphase::FiveEq;
            let liquid = match std::env::var("LIQ").as_deref() {
                Ok("ideal") => Nasg::ideal(1.4),
                Ok("ideal44") => Nasg::ideal(4.4),
                _ => WATER_SG,
            };
            let kapila = std::env::var("NOKAPILA").is_err();
            let sys = FiveEq { eos: [liquid, AIR], kapila };
            let slab = |x: [f64; 3]| {
                let al = 1e-6 + (1.0 - 2e-6) * 0.5 * (((x[0] - 0.3) / 0.01).tanh() - ((x[0] - 0.7) / 0.01).tanh());
                [al * 1000.0, (1.0 - al) * 1.19, 10.0, 0.0, 0.0, 1e5, al]
            };
            let mut s = Solver::new(sys, Mesh::new([n, 1, 1], [0.0; 3], [1.0, 1.0, 1.0]), slab);
            s.cfl = 0.4;
            let cap: f64 = std::env::var("CAP").ok().map(|v| v.parse().unwrap()).unwrap_or(f64::INFINITY);
            if std::env::var("EXPL").is_err() {
                let theta = std::env::var("THETA").ok().map(|v| v.parse().unwrap());
                s.allspeed = Some(AllSpeedOpts { cfl: 0.4, max_acoustic_cfl: cap, theta, ..Default::default() });
            }
            let nsteps: usize = std::env::var("NS").ok().map(|v| v.parse().unwrap()).unwrap_or(20);
            for st in 0..nsteps {
                let dt = s.dt_cfl();
                s.step(dt);
                let (mut dp, mut du, mut at) = (0.0f64, 0.0f64, 0);
                for i in 0..n {
                    let w = s.prim(i, 0, 0);
                    let e = (w[5] / 1e5 - 1.0).abs();
                    if !(e <= dp) {
                        dp = e;
                        at = i;
                    }
                    du = du.max((w[2] / 10.0 - 1.0).abs());
                }
                println!("step {st}: dt {dt:.3e}, max |p/p0 − 1| {dp:.3e} at cell {at} (α {:.4}), |u/u0 − 1| {du:.3e}, CG its total {}", s.prim(at, 0, 0)[6], s.allspeed_stats.cg_iterations);
                if std::env::var("PROFILE_AT").ok().and_then(|v| v.parse::<usize>().ok()) == Some(st) {
                    for i in at.saturating_sub(8)..(at + 8).min(n) {
                        let w = s.prim(i, 0, 0);
                        println!("    cell {i}: α {:.6} ρ {:.4} u−10 {:+.3e} p−1e5 {:+.3e}", w[6], w[0] + w[1], w[2] - 10.0, w[5] - 1e5);
                    }
                }
            }
        }
        _ => panic!("unknown case"),
    }
}
