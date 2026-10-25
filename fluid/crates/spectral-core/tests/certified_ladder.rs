//! The certified resolution ladder: the engine escalates the grid exactly
//! until a *rigorous* bound meets the tolerance, and the bound is checked
//! against a much finer reference run.

use emergence::{Request, Rigor};
use spectral_core::ladder::*;

#[test]
fn engine_escalates_resolution_until_the_rigorous_bound_meets_the_tolerance() {
    // Taylor–Green, Re = 200, to t = 0.5 (the 16³ certificate gives ~3.5e-4
    // there, 32³ ~3e-7; analysis/certificate_results.md).
    let flow = Flow {
        initial: Initial::TaylorGreen,
        nu: 1.0 / 200.0,
        t_end: 0.5,
        dt: 0.02,
    };
    let reference = run(&flow, 48, false);
    let ladder = ladder(&[8, 12, 16, 24], 0);
    let mut chosen = Vec::new();
    for &tol in &[1e-2, 1e-3, 3e-5] {
        let s = ladder.solve(&flow, Request::with_rigor(tol, Rigor::Rigorous));
        let err = l2_distance(&s.value, &reference);
        println!(
            "tol {tol:e}: {} (bound {:.3e}, distance to 48³ {:.3e}, energy {:.8} ± {:.2e}), cost incl. failed tries {:.2e}",
            ladder.rung_name(s.rung),
            s.claim.error,
            err,
            s.value.energy(),
            s.value.energy_bound(),
            s.cost
        );
        assert!(s.met, "tolerance {tol} not certified");
        assert_eq!(s.claim.rigor, Rigor::Rigorous);
        // The 48³ reference is itself not exact; its error is far below these
        // bounds (spectral convergence), so the distance must sit under them.
        assert!(err <= s.claim.error, "distance {err:e} above the rigorous bound {:e}", s.claim.error);
        chosen.push(s.value.grid.n);
    }
    assert!(chosen.windows(2).all(|w| w[0] <= w[1]), "tighter tolerances must not pick coarser grids: {chosen:?}");
    assert!(chosen[0] < chosen[2], "the tolerances should select different grids: {chosen:?}");
    println!("{}", ladder.report());
}
