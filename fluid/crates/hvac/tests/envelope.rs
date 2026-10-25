//! Assemblies: steady U, the exact periodic response, and the
//! finite-volume wall against it.

use hvac::envelope::{Assembly, Flow, Layer, Wall};
use std::f64::consts::PI;

fn masonry() -> Assembly {
    // Inside to outside: plaster, brick, insulation, brick.
    Assembly::new(
        vec![
            Layer { thickness: 0.013, k: 0.57, rho: 1300.0, cp: 1000.0 },
            Layer { thickness: 0.100, k: 0.77, rho: 1700.0, cp: 1000.0 },
            Layer { thickness: 0.080, k: 0.035, rho: 30.0, cp: 1450.0 },
            Layer { thickness: 0.100, k: 0.77, rho: 1700.0, cp: 1000.0 },
        ],
        Flow::Horizontal,
    )
}

#[test]
fn steady_and_periodic_limits() {
    let a = masonry();
    let r = 0.13 + 0.013 / 0.57 + 0.1 / 0.77 + 0.08 / 0.035 + 0.1 / 0.77 + 0.04;
    assert!((a.u_value() - 1.0 / r).abs() < 1e-15);
    // Very slow swings pass straight through (decrement → 1); the lag tends
    // to the wall's finite thermal delay, not to zero.
    let slow = a.periodic(1e9);
    assert!((slow.decrement - 1.0).abs() < 1e-6 && slow.lag < 86_400.0, "{slow:?}");
    // A day: heavy masonry damps and delays.
    let day = a.periodic(86_400.0);
    println!("masonry wall: U {:.3} W/m²K, decrement {:.3}, lag {:.2} h, internal admittance {:.2} W/m²K", a.u_value(), day.decrement, day.lag / 3600.0, day.admittance);
    assert!(day.decrement < 0.5 && day.lag > 4.0 * 3600.0);
}

/// Periodic steady state of the finite-volume wall under an outdoor swing
/// of 10 K amplitude: amplitude and lag of the heat flux into a room held
/// at 0 °C.
fn simulated(a: &Assembly, nodes: usize, steps_per_day: usize) -> (f64, f64) {
    let p = 86_400.0;
    let dt = p / steps_per_day as f64;
    let mut wall = Wall::new(a, nodes, 0.0);
    let tout = |t: f64| 10.0 * (2.0 * PI * t / p).sin();
    // Spin up 12 days, then record one.
    let mut t = 0.0;
    for _ in 0..12 * steps_per_day {
        t += dt;
        wall.step(dt, 0.0, tout(t), 0.0);
    }
    let (mut s, mut c) = (0.0, 0.0);
    for _ in 0..steps_per_day {
        t += dt;
        let q = wall.step(dt, 0.0, tout(t), 0.0);
        let ph = 2.0 * PI * t / p;
        s += q * ph.sin();
        c += q * ph.cos();
    }
    let (s, c) = (2.0 * s / steps_per_day as f64, 2.0 * c / steps_per_day as f64);
    let amp = s.hypot(c);
    // q = amp·sin(ωt − φ): φ = −atan2(c, s)
    let mut lag = -c.atan2(s) / (2.0 * PI) * p;
    while lag < 0.0 {
        lag += p;
    }
    (amp, lag)
}

#[test]
fn finite_volume_wall_converges_to_the_exact_periodic_response() {
    let a = masonry();
    let exact = a.periodic(86_400.0);
    let (amp_exact, lag_exact) = (10.0 * exact.transmittance, exact.lag);
    let mut runs = Vec::new();
    for (nodes, steps) in [(2, 96), (4, 192), (8, 384)] {
        let (amp, lag) = simulated(&a, nodes, steps);
        println!(
            "{nodes} nodes/layer, {steps} steps/day: amplitude {amp:.5} W/m² (exact {amp_exact:.5}, {:+.2e}), lag {:.3} h (exact {:.3} h)",
            amp / amp_exact - 1.0,
            lag / 3600.0,
            lag_exact / 3600.0
        );
        runs.push((amp, lag));
    }
    let e: Vec<f64> = runs.iter().map(|r| (r.0 / amp_exact - 1.0).abs()).collect();
    println!("amplitude error ratios {:.2}, {:.2} (second order: 4)", e[0] / e[1], e[1] / e[2]);
    assert!(e[0] / e[1] > 3.0 && e[1] / e[2] > 3.0, "not second order");
    // Richardson on the two finest runs must land on the exact solution
    // far closer than either run.
    let rich = (4.0 * runs[2].0 - runs[1].0) / 3.0;
    let rich_lag = (4.0 * runs[2].1 - runs[1].1) / 3.0;
    println!("Richardson: amplitude {rich:.5} ({:+.1e}), lag {:.3} h", rich / amp_exact - 1.0, rich_lag / 3600.0);
    assert!((rich / amp_exact - 1.0).abs() < 0.1 * e[2]);
    assert!((rich_lag - lag_exact).abs() < 0.1 * (runs[2].1 - lag_exact).abs().max(1.0));
}

#[test]
fn steady_state_flux_is_u_delta_t() {
    let a = masonry();
    let mut wall = Wall::new(&a, 4, 0.0);
    let mut q = 0.0;
    for _ in 0..2000 {
        q = wall.step(3600.0, 20.0, -5.0, 0.0);
    }
    // Heat flows out of the room: −U·ΔT.
    assert!((q + a.u_value() * 25.0).abs() < 1e-9 * a.u_value() * 25.0, "{q}");
}
