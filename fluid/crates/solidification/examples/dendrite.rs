//! 2D free dendrite in Karma–Rappel units, the Table II benchmark of
//! Karma & Rappel, Phys. Rev. E **57**, 4323 (1998): `Delta = 0.55`,
//! `eps4 = 0.05`, and the KR parameter point `W0 = 1`, `D = 2`,
//! `d0 = a1 a2 W0/D = 0.277` (which makes `tau0 = 1`). Run in the moving frame
//! of KR Sec. III so the tip stays in a small quarter box.
//!
//! The corrected operator `PhaseField::step_kr` (variational anisotropy term
//! of KR Eq. 88) is used; the plain `PhaseField::step` passes through the
//! right tip speed and then decays without ever selecting a parabolic tip.
//! The matching GPU convergence study lives in `solidification-cuda`
//! (`examples/dendrite_convergence.rs`); on this point it gives
//! `V d0/D ~ 0.0164-0.0172` versus KR's 0.0168 (Green-function 0.0170).
//!
//! Run with `cargo run --release -p solidification --example dendrite`.
//! Environment overrides: `D`, `W0`, `D0`, `DX`, `NX`, `NY`, `STEPS`,
//! `REPORT`, `SHIFT` (moving frame, default 1), `DT`, `SEED_R`.

use solidification::phase_field::{Bc, PhaseField};
use std::time::Instant;

const EPS4: f64 = 0.05;
const DELTA: f64 = 0.55;

fn env_f64(name: &str, default: f64) -> f64 {
    std::env::var(name).ok().and_then(|s| s.parse().ok()).unwrap_or(default)
}

fn parse<T: std::str::FromStr>(args: &[String], i: usize, default: T) -> T {
    args.get(i).and_then(|s| s.parse().ok()).unwrap_or(default)
}

/// Translate the fields `sc` cells in `-x`, filling the trailing columns with
/// liquid at the far-field temperature (KR moving-frame technique). Ghosts are
/// refreshed by the next `step_kr`.
fn translate(pf: &mut PhaseField, sc: usize) {
    let m = pf.nx + 2;
    for j in 0..(pf.ny + 2) {
        let base = j * m;
        for i in 1..=(pf.nx - sc) {
            pf.p[base + i] = pf.p[base + i + sc];
            pf.u[base + i] = pf.u[base + i + sc];
        }
        for i in (pf.nx - sc + 1)..=pf.nx {
            pf.p[base + i] = -1.0;
            pf.u[base + i] = pf.u_far;
        }
    }
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let d = parse(&a, 1, env_f64("D", 2.0));
    let w0 = parse(&a, 2, env_f64("W0", 1.0));
    // KR Table II: d0 = a1 a2 W0 / D.
    let d0 = parse(
        &a,
        3,
        env_f64(
            "D0",
            solidification::phase_field::A1 * solidification::phase_field::A2 * w0 / d,
        ),
    );
    let dx = parse(&a, 4, env_f64("DX", 0.4));
    let nx = parse(&a, 5, env_f64("NX", 400.0) as usize);
    let ny = parse(&a, 6, env_f64("NY", 200.0) as usize);
    let steps = parse(&a, 7, env_f64("STEPS", 100_000.0) as usize);
    let report = parse(&a, 8, env_f64("REPORT", 2000.0) as usize);
    let shift_on = env_f64("SHIFT", 1.0) != 0.0;

    let mut pf = PhaseField::new(nx, ny, dx, d0, d, w0, EPS4);
    pf.bc = Bc::quadrant();
    pf.u_far = -DELTA;
    let seed_r = env_f64("SEED_R", 3.0 * w0);
    pf.seed_corner(seed_r);
    let dt = env_f64("DT", 0.0);
    let dt = if dt > 0.0 { dt } else { pf.stable_dt() };
    let lx = nx as f64 * dx;
    println!(
        "KR dendrite: D={d} W0={w0} d0={d0:.4} W0/d0={:.3} dx={dx} W0/dx={:.2} N={nx}x{ny} dt={dt:.5e} shift={shift_on}",
        w0 / d0,
        w0 / dx
    );

    let mut last_tip = pf.crossing_x(1, 0.0).unwrap_or(0.0);
    let mut last_t = 0.0;
    let mut shift = 0.0;
    let start = Instant::now();
    for k in 0..steps {
        pf.step_kr(dt);
        let t = (k + 1) as f64 * dt;
        if (k + 1) % report != 0 && k + 1 != steps {
            continue;
        }
        let tip = pf.crossing_x(1, 0.0).unwrap_or(f64::NAN);
        let r = pf.tip_radius();
        let v = (tip - last_tip) / (t - last_t);
        let vd = v * d0 / d;
        let sel = v * r * r / (2.0 * d * d0);
        println!(
            "step={:>8} t={:9.2} tip={:8.3} v={v:9.5} Vd0/D={vd:9.5} R={r:8.3} sel={sel:7.3}",
            k + 1,
            t,
            tip + shift
        );
        last_tip = tip;
        last_t = t;
        if shift_on && tip > 0.55 * lx {
            let sc = (((tip - 0.25 * lx) / dx).round() as usize).clamp(1, nx / 2);
            translate(&mut pf, sc);
            shift += sc as f64 * dx;
            last_tip -= sc as f64 * dx;
        }
        if start.elapsed().as_secs_f64() > 170.0 {
            println!("wall budget reached after {} steps", k + 1);
            break;
        }
    }
    println!("published (Karma & Rappel 1998, Table II): V d0/D = 0.0168 (Green function 0.0170)");
}
