//! Refrigerant equations of state against CoolProp's own evaluation of the
//! same formulations (values recorded by `data/build_refrigerants.py`).
//!
//! Tolerances by derivative order: pressure is α's first δ-derivative,
//! enthalpy and entropy first τ-derivatives (1e-9: two implementations
//! summing the same terms in a different order); c_p and the sound speed
//! need second derivatives and cancel more (1e-8).

use compressible::helmholtz::Helmholtz;
use compressible::refrigerant_data::*;

fn single_phase(eos: &Helmholtz, states: &[State]) -> [f64; 5] {
    let mut worst = [0.0f64; 5];
    for s in states {
        let pr = eos.props(s.rho, s.t);
        let errs = [(pr.p / s.p - 1.0).abs(), ((pr.h - s.h) / s.h.abs().max(1e5)).abs(), ((pr.s - s.s) / s.s.abs().max(1e3)).abs(), (pr.cp / s.cp - 1.0).abs(), (pr.w / s.w - 1.0).abs()];
        for (w, e) in worst.iter_mut().zip(errs) {
            *w = w.max(e);
        }
    }
    worst
}

fn check(eos: Helmholtz, sat: &[Sat], states: &[State]) {
    let w = single_phase(&eos, states);
    let mut ws = 0.0f64;
    for s in sat {
        let (p, rl, rv) = eos.saturation(s.t).expect("saturation");
        let (l, v) = (eos.props(rl, s.t), eos.props(rv, s.t));
        for e in [p / s.p - 1.0, rl / s.rho_l - 1.0, rv / s.rho_v - 1.0, (l.h - s.h_l) / s.h_l.abs().max(1e5), (v.h - s.h_v) / s.h_v.abs().max(1e5), (l.s - s.s_l) / s.s_l.abs().max(1e3), (v.s - s.s_v) / s.s_v.abs().max(1e3)] {
            ws = ws.max(e.abs());
        }
    }
    println!("{}: worst relative error p {:.1e}, h {:.1e}, s {:.1e}, cp {:.1e}, w {:.1e}; saturation {:.1e}", eos.name, w[0], w[1], w[2], w[3], w[4], ws);
    assert!(w[0] < 1e-9 && w[1] < 1e-9 && w[2] < 1e-9 && w[3] < 1e-8 && w[4] < 1e-8, "{} single-phase {w:?}", eos.name);
    assert!(ws < 1e-8, "{} saturation {ws:.2e}", eos.name);
}

#[test]
fn refrigerants_match_coolprop() {
    check(r32(), R32_SAT, R32_STATES);
    check(r134a(), R134A_SAT, R134A_STATES);
    check(r290(), R290_SAT, R290_STATES);
    check(r1234yf(), R1234YF_SAT, R1234YF_STATES);
    check(r744(), R744_SAT, R744_STATES);
}

#[test]
fn r410a_pseudo_pure() {
    // A blend treated as one fluid. Single-phase states must match exactly;
    // at saturation CoolProp takes bubble and dew pressures from ancillary
    // fits, and the one Maxwell pressure of the pseudo-pure equation must
    // lie between them.
    let eos = r410a();
    let w = single_phase(&eos, R410A_STATES);
    assert!(w[0] < 1e-9 && w[1] < 1e-9 && w[2] < 1e-9 && w[3] < 1e-8 && w[4] < 1e-8, "single-phase {w:?}");
    for s in R410A_SAT {
        let (p, _, _) = eos.saturation(s.t).unwrap();
        let (lo, hi) = (s.p.min(s.p_dew), s.p.max(s.p_dew));
        println!("R410A T {:.1} K: bubble {:.1} kPa, Maxwell {:.1} kPa, dew {:.1} kPa", s.t, s.p / 1e3, p / 1e3, s.p_dew / 1e3);
        assert!(p >= lo * (1.0 - 1e-6) && p <= hi * (1.0 + 1e-6), "Maxwell pressure outside the bubble–dew interval");
    }
}
