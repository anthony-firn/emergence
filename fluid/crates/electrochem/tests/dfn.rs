//! Tests and fold benchmarks for the pseudo-2D DFN model
//! ([`electrochem::dfn`]).
//!
//! The DFN is the reference rung above SPM/SPMe. These tests check the laws
//! the DFN must satisfy (lithium conservation, charge balance), that it
//! reduces to SPM/SPMe as the C-rate vanishes with the expected scaling, that
//! it converges under refinement, and that the coarse folds (SPM, SPMe, ECM)
//! lose accuracy as the indicator grows.

use electrochem::dfn::{solve_equilibrated, Dfn, DfnParams};
use electrochem::ecm::{extract_rc, TheveninRc};
use electrochem::spm::Spm;
use electrochem::spme::Spme;

/// Build a DFN cell, run it at `c_rate` for `t_end` seconds, and return the
/// terminal voltage and reaction non-uniformity.
fn run_dfn(
    n: (usize, usize, usize),
    m: usize,
    c_rate: f64,
    dt: f64,
    t_end: f64,
) -> (f64, f64) {
    let mut cell = Dfn::new(DfnParams::lg_m50(), n.0, n.1, n.2, m);
    let i_app = c_rate * cell.params.cell.i_ref();
    let steps = (t_end / dt).round() as usize;
    for _ in 0..steps {
        cell.step(dt, i_app);
    }
    (cell.voltage(), cell.reaction_nonuniformity())
}

#[test]
fn dense_solver_round_trip() {
    let a = vec![1.0, 2.0, 3.0, 4.0];
    let mut b = vec![5.0, 6.0];
    let x = solve_equilibrated(&a, 2, &mut b);
    assert!((x[0] - (-4.0)).abs() < 1e-9);
    assert!((x[1] - 4.5).abs() < 1e-9);
}

#[test]
fn lithium_is_conserved_to_round_off() {
    let p = DfnParams::lg_m50();
    let mut cell = Dfn::new(p, 5, 2, 5, 6);
    let i_app = cell.params.cell.i_ref();
    let sl0 = cell.solid_lithium();
    let el0 = cell.electrolyte_lithium();
    // Solid and electrolyte inventories are each conserved: the reaction moves
    // Li between the two electrodes, and the (1−t₊) electrolyte source/sink
    // balances across the cell.
    for _ in 0..30 {
        cell.step(10.0, i_app);
    }
    let sdrift = (cell.solid_lithium() - sl0).abs() / sl0;
    let edrift = (cell.electrolyte_lithium() - el0).abs() / el0;
    assert!(sdrift < 1e-9, "solid Li drift {sdrift}");
    assert!(edrift < 1e-9, "electrolyte Li drift {edrift}");
}

#[test]
fn integrated_reaction_current_balances_the_applied_current() {
    let p = DfnParams::lg_m50();
    let mut cell = Dfn::new(p, 5, 2, 5, 6);
    let i_app = 3.0 * cell.params.cell.i_ref();
    cell.step(10.0, i_app);
    let (neg, pos) = cell.integrated_reaction_current();
    // Discharge: anodic (positive) in the negative, cathodic in the positive.
    assert!((neg - i_app).abs() / i_app < 1e-8, "neg {neg} vs {i_app}");
    assert!((pos + i_app).abs() / i_app < 1e-8, "pos {pos} vs {}", -i_app);
    assert!((neg + pos).abs() / i_app < 1e-12, "net reaction current {}", neg + pos);
}

#[test]
fn reduces_to_spm_and_spme_at_low_rate() {
    // The DFN resolves the solid/electrolyte ohmic drop internally, so the
    // SPM/SPMe comparison sets the crate's lumped ESTIMATE r_ohmic to zero to
    // avoid double counting; what is left is the reaction/electrolyte physics.
    let p = DfnParams::lg_m50();
    let mut cp = p.cell.clone();
    cp.r_ohmic = 0.0;
    let dt = 20.0;
    let t_end = 600.0;
    let rate = |cr: f64| -> (f64, f64) {
        let mut dfn = Dfn::new(DfnParams::lg_m50(), 6, 2, 6, 6);
        let i_app = cr * p.cell.i_ref();
        let mut spm = Spm::new(cp.clone(), 40);
        let mut spme = Spme::new(cp.clone(), 40, 80, p.t_plus);
        let steps = (t_end / dt) as usize;
        let mut d_spm = 0.0f64;
        let mut d_spme = 0.0f64;
        for _ in 0..steps {
            dfn.step(dt, i_app);
            spm.step(dt, i_app);
            spme.step(dt, i_app);
            d_spm = d_spm.max((dfn.voltage() - spm.voltage(i_app)).abs());
            d_spme = d_spme.max((dfn.voltage() - spme.voltage(i_app)).abs());
        }
        (d_spm, d_spme)
    };
    let (d_spm_lo, d_spme_lo) = rate(0.01);
    let (d_spm_hi, d_spme_hi) = rate(0.02);
    println!(
        "low-rate: C/100 SPM {:.3} mV SPMe {:.3} mV; C/50 SPM {:.3} mV SPMe {:.3} mV",
        d_spm_lo * 1e3,
        d_spme_lo * 1e3,
        d_spm_hi * 1e3,
        d_spme_hi * 1e3
    );
    assert!(d_spm_lo < 3e-3, "C/100 SPM error {d_spm_lo}");
    assert!(d_spme_lo < 3e-3, "C/100 SPMe error {d_spme_lo}");
    // Errors scale linearly with C-rate (the expected fold scaling).
    let ratio = d_spm_hi / d_spm_lo;
    assert!((1.5..2.5).contains(&ratio), "SPM error ratio {ratio}");
    let ratio_e = d_spme_hi / d_spme_lo;
    assert!((1.5..2.5).contains(&ratio_e), "SPMe error ratio {ratio_e}");
    // SPMe is never worse than SPM at these low rates.
    assert!(d_spme_lo <= d_spm_lo * 1.05, "SPMe {} vs SPM {}", d_spme_lo, d_spm_lo);
}

#[test]
fn spatial_convergence() {
    // Voltage after 300 s at 1 C on a sequence of grids, versus the finest.
    let dt = 5.0;
    let t_end = 300.0;
    let reference = run_dfn((12, 3, 12), 10, 1.0, dt, t_end).0;
    let seq = [((4, 2, 4), 4), ((6, 2, 6), 5), ((8, 2, 8), 6)];
    let mut errs = Vec::new();
    for &(b, m) in &seq {
        let v = run_dfn(b, m, 1.0, dt, t_end).0;
        errs.push((v - reference).abs());
    }
    for w in errs.windows(2) {
        assert!(w[1] < w[0], "spatial error not decreasing: {errs:?}");
    }
    println!("spatial errors (V): {errs:?}");
    assert!(errs[2] < 2e-2, "finest spatial error {} V", errs[2]);
}

#[test]
fn temporal_convergence() {
    // Backward Euler is first order: halving dt should roughly halve the error.
    let grid = (6, 2, 6);
    let m = 6;
    let t_end = 400.0;
    let reference = run_dfn(grid, m, 1.0, 2.5, t_end).0;
    let dts = [40.0, 20.0, 10.0];
    let mut errs = Vec::new();
    for &dt in &dts {
        let v = run_dfn(grid, m, 1.0, dt, t_end).0;
        errs.push((v - reference).abs());
    }
    for w in errs.windows(2) {
        assert!(w[1] < w[0], "temporal error not decreasing: {errs:?}");
    }
    println!("temporal errors (V): {errs:?}");
    // Rough first-order check between the last two.
    let ratio = errs[1] / errs[2];
    assert!((1.3..3.5).contains(&ratio), "temporal order ratio {ratio}");
}

#[test]
fn reaction_nonuniformity_grows_with_c_rate() {
    // The position dependence SPMe cannot represent: evaluated early, before
    // the low-rate profile has had time to build.
    let dt = 20.0;
    let t_end = 100.0;
    let (_, u_lo) = run_dfn((8, 3, 8), 8, 0.05, dt, t_end);
    let (_, u_hi) = run_dfn((8, 3, 8), 8, 5.0, dt, t_end);
    assert!(u_hi > 1.5 * u_lo, "nonuniformity {u_lo} -> {u_hi}");
}

/// The fold table: voltage error of SPM and SPMe relative to the DFN over a
/// 600 s constant-current window (with the double-counted lumped R_Ω
/// removed), plus the ECM 60 s pulse error. Also checks that the errors grow
/// with the C-rate.
#[test]
fn fold_table() {
    let p = DfnParams::lg_m50();
    let mut cp = p.cell.clone();
    cp.r_ohmic = 0.0;
    let dt: f64 = 20.0;
    // A fixed 600 s constant-current window isolates the polarisation response
    // (at the same instant the same charge has passed in every model).
    let t_end: f64 = 600.0;
    let rates = [0.1f64, 0.5, 1.0, 2.0, 3.0, 5.0];
    println!("C-rate | V_DFN | SPM err mV | SPMe err mV | unif | ECM 60s err mV");
    let mut last_spm = 0.0f64;
    let mut ecm_lo = 0.0f64;
    let mut ecm_hi = 0.0f64;
    for &cr in &rates {
        let i_app = cr * p.cell.i_ref();
        let nsteps = (t_end / dt).round() as usize;
        let mut dfn = Dfn::new(DfnParams::lg_m50(), 8, 3, 8, 8);
        let mut spm = Spm::new(cp.clone(), 60);
        let mut spme = Spme::new(cp.clone(), 60, 120, p.t_plus);
        let mut d_spm = 0.0f64;
        let mut d_spme = 0.0f64;
        let mut unif_early = 0.0f64;
        for k in 0..nsteps {
            dfn.step(dt, i_app);
            spm.step(dt, i_app);
            spme.step(dt, i_app);
            if k == 5 {
                unif_early = dfn.reaction_nonuniformity();
            }
            d_spm = d_spm.max((dfn.voltage() - spm.voltage(i_app)).abs());
            d_spme = d_spme.max((dfn.voltage() - spme.voltage(i_app)).abs());
        }
        let v_dfn = dfn.voltage();
        let unif = unif_early;
        // ECM fold: fit one RC branch to a 60 s DFN pulse (tiny first sample
        // gives R0) and measure the pulse fit error.
        let pulse_dt = 2.0;
        let ocv0 = Dfn::new(DfnParams::lg_m50(), 8, 3, 8, 8).voltage();
        let mut pulse = Dfn::new(DfnParams::lg_m50(), 8, 3, 8, 8);
        let mut times = vec![0.01];
        pulse.step(0.01, i_app);
        let mut dv = vec![ocv0 - pulse.voltage()];
        let fit_steps = (60.0 / pulse_dt) as usize;
        let mut samples = vec![(0.01, pulse.voltage())];
        for k in 0..fit_steps {
            pulse.step(pulse_dt, i_app);
            let t = (k as f64 + 1.0) * pulse_dt + 0.01;
            times.push(t);
            dv.push(ocv0 - pulse.voltage());
            samples.push((t, pulse.voltage()));
        }
        let (r0, r1, tau) = extract_rc(&times, &dv, i_app);
        let mut ecm = TheveninRc::new(r0, r1, tau / r1, i_app * 3600.0, 1.0);
        let mut d_ecm = 0.0f64;
        for (i, &(t, v)) in samples.iter().enumerate() {
            let v_ecm = ocv0 - i_app * r0 - ecm.v_rc;
            d_ecm = d_ecm.max((v_ecm - v).abs());
            if i + 1 < samples.len() {
                ecm.step(samples[i + 1].0 - t, i_app);
            }
        }
        println!(
            "{cr:5.1}  | {v_dfn:.4} | {:9.2} | {:10.2} | {unif:.3} | {:9.2}",
            d_spm * 1e3,
            d_spme * 1e3,
            d_ecm * 1e3
        );
        if cr == 0.1 {
            ecm_lo = d_ecm;
        }
        if cr == 5.0 {
            ecm_hi = d_ecm;
        }
        last_spm = d_spm;
    }
    assert!(last_spm > 5.0 * 1e-3, "SPM error should reach >5 mV at 5C");
    assert!(ecm_hi > 3.0 * ecm_lo, "ECM fold error should grow: {ecm_lo} -> {ecm_hi}");
}

/// Fine-grid 1C discharge to 2.5 V, and the SPM capacity for comparison.
/// Ignored by default (dense Newton on a fine grid); run with
/// `cargo test -p electrochem --test dfn -- --ignored --nocapture`.
#[test]
#[ignore]
fn fine_grid_full_discharge() {
    let p = DfnParams::lg_m50();
    let mut cell = Dfn::new(p.clone(), 12, 4, 12, 10);
    let i_app = p.cell.i_ref();
    let dt = 10.0;
    let mut t = 0.0;
    let mut n = 0;
    while cell.voltage() > 2.5 && n < 3000 {
        cell.step(dt, i_app);
        t += dt;
        n += 1;
    }
    let cap = i_app * t / 3600.0 * 1000.0; // mAh m⁻²
    println!(
        "fine DFN 1C: capacity = {cap:.0} mAh m^-2, t_end = {t:.0} s, V_end = {:.3}, unif = {:.3}",
        cell.voltage(),
        cell.reaction_nonuniformity()
    );
    let spm = Spm::new(p.cell.clone(), 60);
    let (cap_spm, _) = spm.discharge(1.0, 2.5, 10.0, 2000);
    println!("SPM 1C capacity = {cap_spm:.0} mAh m^-2 (analysis/electrochem.md)");
    assert!(cap > 0.0 && cap.is_finite());
}
