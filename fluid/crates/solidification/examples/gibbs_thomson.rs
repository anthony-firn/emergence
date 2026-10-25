//! Gibbs–Thomson fold: the critical radius of a circular seed, and its
//! convergence as the interface width `W/d0 -> 0`.
//!
//! A seed of radius `R` in a melt held at `u_far = -Delta` grows when
//! `R > R_c = d0/Delta` and shrinks below it (the Gibbs–Thomson condition
//! `u_i = -d0/R` balances the far field). With the Karma–Rappel `lambda`–`W`
//! choice the interface kinetics vanish, so `R_c -> d0/Delta` as `W/d0 -> 0`,
//! with a finite-`W` offset that this example measures.
//!
//! Run with `cargo run --release -p solidification --example gibbs_thomson`.

use solidification::phase_field::{Bc, PhaseField};

const D0: f64 = 1.0;
const DIFF: f64 = 2.0;
const DELTA: f64 = 0.3;

/// Signed front displacement for a seed of radius `r0` after a fixed time.
fn displacement(r0: f64, w0: f64, n: usize, dx: f64, steps: usize) -> f64 {
    let mut pf = PhaseField::new(n, n, dx, D0, DIFF, w0, 0.0);
    pf.bc = Bc::far_field();
    pf.u_far = -DELTA;
    let xc = 0.5 * pf.length_x();
    pf.seed_circle(xc, xc, r0);
    let jmid = n / 2 + 1;
    let r_start = pf.crossing_x(jmid, xc).unwrap_or(0.0);
    let dt = pf.stable_dt();
    for _ in 0..steps {
        pf.step(dt);
    }
    let r_end = pf.crossing_x(jmid, xc).unwrap_or(0.0);
    r_end - r_start
}

/// Critical radius by bisection on the sign of the displacement.
fn critical_radius(w0: f64, n: usize, dx: f64, steps: usize) -> f64 {
    let mut lo = 1.0; // shrinks
    let mut hi = 6.0; // grows
    for _ in 0..14 {
        let mid = 0.5 * (lo + hi);
        if displacement(mid, w0, n, dx, steps) > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    0.5 * (lo + hi)
}

fn main() {
    println!(
        "d0={D0} D={DIFF} Delta={DELTA} -> d0/Delta = {}",
        D0 / DELTA
    );
    let n = 96;
    let dx = 0.25;
    let steps = 4000;
    for &w0 in &[0.8, 0.4, 0.2, 0.1] {
        let rc = critical_radius(w0, n, dx, steps);
        println!(
            "W0/d0={w0:.2}: Rc={rc:.4}  Rc/(d0/Delta)={:.4}  offset={:+.2}%",
            rc / (D0 / DELTA),
            (rc - D0 / DELTA) / (D0 / DELTA) * 100.0
        );
    }
    // Kinetic check: with a deliberately large tau0 the interface is no longer
    // kinetic-free and the critical radius must shift. This is not run here;
    // see the `kinetic_choice` test for the unit-scale version.
}
