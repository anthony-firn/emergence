//! Tests for the 2D Boussinesq convection solver and the conduction→convection
//! unfold (`compressible::convection2d`).
//!
//! Benchmarks:
//!   * de Vahl Davis (1983) differentially heated square cavity, `Pr = 0.71`:
//!     average Nusselt number `1.118 / 2.243 / 4.519` at `Ra = 10³/10⁴/10⁵` and
//!     the mid-plane velocity maxima, with grid convergence;
//!   * Rayleigh–Bénard onset between rigid plates: the perturbation growth
//!     rate changes sign at `Ra_c = 1707.762` (`k_c = 3.116`);
//!   * coupling to the isobaric enthalpy rung with a Carman–Kozeny Darcy drag
//!     (gallium melting, Gau & Viskanta 1986 / Brent, Voller & Reid 1988);
//!   * the unfold indicator: the flow speed and the melt front for a
//!     sub-critical and a super-critical gallium cell;
//!   * energy conservation of the coupled step to round-off.

use compressible::convection2d::{
    benard_growth, Boussinesq2D, EnthalpyConvection2D, K_CRIT_RIGID, RA_CRIT_RIGID,
};
use compressible::elements::{by_symbol, fit, P_ATM};
use compressible::isobaric::{Ladder, PhaseSet};
use compressible::isobaric2d::{convection_regime, rayleigh_benard, Bc2D, ConvectionRegime, Isobaric2D};

// ------------------------------------------------------------- de Vahl Davis

/// Run the cavity at `(Ra, n)` to steady state and return
/// `(Nu, u_max, v_max, steps)`. `u_max`/`v_max` are the mid-plane maxima.
fn run_cavity(ra: f64, n: usize, tol: f64, max_steps: usize) -> (f64, f64, f64, usize) {
    let mut s = Boussinesq2D::cavity(ra, n);
    let (steps, nu) = s.run_cavity_steady(tol, max_steps);
    let (um, vm) = s.midline_maxima();
    (nu, um, vm, steps)
}

/// de Vahl Davis (1983) benchmark table (`Pr = 0.71`, square cavity, heated
/// from the left): average Nusselt number and mid-plane velocity maxima.
/// `Ra = 10³`: Nu = 1.118, u = 3.649, v = 3.697;
/// `Ra = 10⁴`: Nu = 2.243, u = 16.178, v = 19.617;
/// `Ra = 10⁵`: Nu = 4.519, u = 34.73, v = 68.59.
#[test]
fn cavity_ra_1e3_matches_de_vahl_davis() {
    let (nu, um, vm, steps) = run_cavity(1e3, 24, 1e-8, 500);
    println!("cavity Ra=1e3 n=24: Nu={nu:.4} (ref 1.118), u_max={um:.3} (3.649), v_max={vm:.3} (3.697), {steps} outer steps");
    assert!((nu / 1.118 - 1.0).abs() < 0.02, "Nu {nu}");
    assert!((um / 3.649 - 1.0).abs() < 0.03, "u_max {um}");
    assert!((vm / 3.697 - 1.0).abs() < 0.03, "v_max {vm}");
}

#[test]
fn cavity_ra_1e4_matches_de_vahl_davis() {
    let (nu, um, vm, steps) = run_cavity(1e4, 40, 1e-7, 1000);
    println!("cavity Ra=1e4 n=40: Nu={nu:.4} (ref 2.243), u_max={um:.3} (16.178), v_max={vm:.3} (19.617), {steps} outer steps");
    assert!((nu / 2.243 - 1.0).abs() < 0.03, "Nu {nu}");
    assert!((um / 16.178 - 1.0).abs() < 0.04, "u_max {um}");
    assert!((vm / 19.617 - 1.0).abs() < 0.04, "v_max {vm}");
}

#[test]
fn cavity_ra_1e5_matches_de_vahl_davis() {
    let (nu, um, vm, steps) = run_cavity(1e5, 48, 1e-6, 1500);
    println!("cavity Ra=1e5 n=48: Nu={nu:.4} (ref 4.519), u_max={um:.3} (34.73), v_max={vm:.3} (68.59), {steps} outer steps");
    assert!((nu / 4.519 - 1.0).abs() < 0.05, "Nu {nu}");
    assert!((um / 34.73 - 1.0).abs() < 0.06, "u_max {um}");
    assert!((vm / 68.59 - 1.0).abs() < 0.06, "v_max {vm}");
}

/// Grid convergence of the Nusselt number at `Ra = 10⁴` (de Vahl Davis'
/// extrapolated value is 2.243). The error must fall towards the benchmark,
/// not away from it.
#[test]
fn cavity_grid_convergence() {
    let (nu24, _, _, _) = run_cavity(1e4, 24, 1e-7, 1500);
    let (nu48, _, _, _) = run_cavity(1e4, 48, 1e-7, 1500);
    println!("cavity grid convergence Ra=1e4: Nu(24)={nu24:.4}, Nu(48)={nu48:.4} (ref 2.243)");
    assert!((nu48 / 2.243 - 1.0).abs() < (nu24 / 2.243 - 1.0).abs(), "not converging: {nu24} -> {nu48}");
    assert!((nu48 / 2.243 - 1.0).abs() < 0.02, "Nu48 {nu48}");
}

/// Converged grid study (slow): `Ra = 10⁵` at `n = 64` and `n = 96`. Run with
/// `cargo test -p compressible --test convection2d -- --ignored --nocapture`.
#[test]
#[ignore = "converged grid study"]
fn cavity_ra_1e5_converged_grid_study() {
    for n in [64usize, 96] {
        let (nu, um, vm, steps) = run_cavity(1e5, n, 1e-6, 4000);
        println!("cavity converged Ra=1e5 n={n}: Nu={nu:.4} (ref 4.519), u_max={um:.3} (34.73), v_max={vm:.3} (68.59), {steps} outer steps");
    }
}

// ------------------------------------------------------------- Rayleigh-Bénard

/// The onset of Rayleigh–Bénard convection between rigid, isothermal plates:
/// at the critical wavenumber `k_c = 3.116` a small perturbation must decay
/// below `Ra_c = 1707.762` and grow above it. This is the sign change the
/// unfold indicator predicts.
#[test]
fn rayleigh_benard_growth_changes_sign_at_ra_c() {
    let pr = 0.71;
    let k = K_CRIT_RIGID;
    let t_end = 6.0;
    let (rate_hi, _) = benard_growth(2000.0, pr, k, 24, t_end);
    let (rate_lo, _) = benard_growth(1400.0, pr, k, 24, t_end);
    println!(
        "Bénard growth at k={k:.3}: Ra=2000 rate={rate_hi:+.4}, Ra=1400 rate={rate_lo:+.4}; Ra_c={RA_CRIT_RIGID}"
    );
    assert!(rate_hi > 0.0, "supercritical perturbation did not grow ({rate_hi})");
    assert!(rate_lo < 0.0, "subcritical perturbation did not decay ({rate_lo})");
    assert!(1400.0 < RA_CRIT_RIGID && RA_CRIT_RIGID < 2000.0);
}

// ------------------------------------------------------------- gallium melting

/// Gallium (element fit, `substances`/`elements.csv`): the metal is
/// ~32 W/(m K) in both condensed phases (Brent, Voller & Reid 1988).
fn gallium_model() -> (PhaseSet, [f64; 3]) {
    let m = fit(by_symbol("Ga").unwrap()).unwrap();
    (PhaseSet::from(&m), [32.0, 32.0, 0.05])
}

/// Melting point, lowest (solid) phase index and the liquid index of the
/// fitted gallium.
fn gallium_phases() -> (PhaseSet, [f64; 3], f64, usize, usize) {
    let (set, k) = gallium_model();
    let ladder = Ladder::new(&set, P_ATM, 250.0, 350.0);
    let t_m = ladder.transitions.first().map_or(302.9, |x| x.0);
    let solid = ladder.lowest;
    let liquid = set.liquid;
    (set, k, t_m, solid, liquid)
}

/// Physical gallium properties for the Boussinesq momentum:
/// `ρ = 6093 kg/m³`, `ν = 3.0e-7 m²/s`, `c_p = 381.5 J/(kg K)`,
/// `k = 32 W/(m K)`, `β = 1.2e-4 /K` (Brent, Voller & Reid 1988, from Gau &
/// Viskanta 1986). `α = k/(ρ c_p) = 1.38e-5 m²/s`, `Pr = ν/α = 0.0216`.
fn gallium_flow(nx: usize, ny: usize, w: f64, hm: f64, t_ref: f64) -> Boussinesq2D {
    let rho = 6093.0;
    let cp = 381.5;
    let k = 32.0;
    let nu = 3.0e-7;
    let alpha = k / (rho * cp);
    let beta = 1.2e-4;
    let mut flow = Boussinesq2D::new(nx, ny, w / nx as f64, hm / ny as f64, nu, alpha, beta, 9.81, t_ref);
    // liquid metal: cell Reynolds number ~10², central advection is unstable
    flow.upwind_momentum = true;
    flow
}

/// Gau & Viskanta (1986) rectangular cell: 6.35 cm wide, 8.94 cm high (aspect
/// ratio 0.71); the heated vertical wall is on the left.
const GALLIUM_W: f64 = 0.0635;
const GALLIUM_H: f64 = 0.0894;

/// Rayleigh number per kelvin of wall superheat, from the physical properties
/// above and the cavity height.
fn gallium_ra_per_k() -> f64 {
    let alpha = 32.0 / (6093.0 * 381.5);
    rayleigh_benard(9.81, 1.2e-4, 1.0, GALLIUM_H, 3.0e-7, alpha)
}

/// Build a gallium melting cell: solid at `T_m`, left wall at `T_m + ΔT`,
/// right wall at `T_m`, adiabatic top/bottom, coupled to a Carman–Kozeny
/// Darcy drag (`C = 1e6`, Voller & Prakash 1987). The Boussinesq flow is
/// constant-density (`pc.expansion = false`).
fn gallium_melt(nx: usize, ny: usize, delta_t: f64) -> (EnthalpyConvection2D, f64) {
    let (set, k, t_m, solid, _) = gallium_phases();
    let mut pc = Isobaric2D::new(
        set,
        P_ATM,
        k.to_vec(),
        nx,
        ny,
        GALLIUM_W / nx as f64,
        GALLIUM_H / ny as f64,
        solid,
        t_m,
        (200.0, 500.0),
    );
    pc.expansion = false;
    pc.bc.west = Bc2D::Dirichlet(t_m + delta_t);
    pc.bc.east = Bc2D::Dirichlet(t_m);
    pc.bc.south = Bc2D::Adiabatic;
    pc.bc.north = Bc2D::Adiabatic;
    for j in 0..ny {
        for i in 0..nx {
            pc.set_cell(i, j, solid, t_m);
        }
    }
    let mut flow = gallium_flow(nx, ny, GALLIUM_W, GALLIUM_H, t_m);
    flow.t_wall = [Some(t_m + delta_t), Some(t_m), None, None];
    let sys = EnthalpyConvection2D::new(pc, flow, 1e6);
    let ra = gallium_ra_per_k() * delta_t;
    (sys, ra)
}

/// The coupled gallium cell runs, develops a convection roll (buoyant melt
/// rises along the heated wall), conserves energy, and advances the melt
/// front. The computed interface is nearly vertical because the low-Prandtl
/// liquid is well mixed; this is reported honestly in `analysis/convection2d.md`
/// rather than tuned.
#[test]
fn gallium_melting_develops_convection() {
    let (nx, ny) = (17, 24);
    let dt = 0.05;
    let t_end = 300.0;
    let delta_t = 8.91; // Ste = c_p ΔT / L = 0.0424 (Brent, Voller & Reid 1988)
    let (mut sys, ra) = gallium_melt(nx, ny, delta_t);
    sys.run(t_end, dt);
    let front = sys.front_x(ny / 2);
    let lean = sys.front_x(ny - 2) - sys.front_x(1);
    let flow_speed = sys.flow.max_speed();
    let melt = sys.melt_fraction();
    // the liquid is warmer at the top than at the bottom (the roll)
    let i = nx / 2;
    let t_bot = sys.flow.t[sys.pc.idx(i, 1)];
    let t_top = sys.flow.t[sys.pc.idx(i, ny - 2)];
    println!(
        "gallium melting: Ra={ra:.3e}, Pr={:.4}, ΔT={delta_t} K, t={:.0}s, front(mid)={front:.4} m, melt={melt:.3}, |u|max={flow_speed:.3e}, lean(top−bottom)={lean:.2e} m, T(top−bottom)={:.3} K",
        3.0e-7 / (32.0 / (6093.0 * 381.5)),
        sys.pc.t,
        t_top - t_bot
    );
    assert!(flow_speed > 1e-3, "no convection developed ({flow_speed})");
    assert!(melt > 0.05, "too little melt {melt}");
    assert!(t_top >= t_bot - 1e-6, "top of the liquid must not be colder");
    // convection makes the interface lean: the top advances past the bottom
    assert!(lean > 0.005, "interface did not lean with convection ({lean} m)");
}

// ------------------------------------------------------------- indicator

/// The cheap indicator (the existing `isobaric2d` functions): the same water
/// layer is conduction-only below `Ra_c` and unfolds to convection above it.
#[test]
fn indicator_selects_conduction_and_convection() {
    let (g, beta, nu, alpha) = (9.81, 2.07e-4, 1e-6, 1.4e-7);
    let l = 1e-2; // 1 cm water layer
    let ra_cold = rayleigh_benard(g, beta, 0.01, l, nu, alpha);
    let ra_hot = rayleigh_benard(g, beta, 1.0, l, nu, alpha);
    println!("indicator: Ra(0.01 K)={ra_cold:.1} -> {:?}, Ra(1 K)={ra_hot:.1} -> {:?}", convection_regime(ra_cold), convection_regime(ra_hot));
    assert_eq!(convection_regime(ra_cold), ConvectionRegime::Conduction);
    assert_eq!(convection_regime(ra_hot), ConvectionRegime::Convection);
}

/// The indicator in action on the gallium cell: below `Ra_c` the coupled flow
/// is negligible and conduction-only is exact; above it the flow is orders of
/// magnitude larger and changes the melt front. We report the two fronts.
#[test]
fn supercritical_unfolds_subcritical_does_not() {
    let (nx, ny) = (17, 24);
    let t_end = 300.0;
    // super-critical physical gallium
    let (mut hot, ra_hot) = gallium_melt(nx, ny, 8.91);
    hot.run(t_end, 0.05);
    let speed_hot = hot.flow.max_speed();
    let melt_hot = hot.melt_fraction();
    let lean_hot = hot.front_x(ny - 2) - hot.front_x(1);
    // sub-critical: ΔT so that Ra ≈ 0.5 Ra_c
    let dt_sub = 0.5 * RA_CRIT_RIGID / gallium_ra_per_k();
    let (mut cold, ra_cold) = gallium_melt(nx, ny, dt_sub);
    cold.run(t_end, 0.05);
    let speed_cold = cold.flow.max_speed();
    let melt_cold = cold.melt_fraction();
    let lean_cold = cold.front_x(ny - 2) - cold.front_x(1);
    println!(
        "indicator unfold: super-critical Ra={ra_hot:.2e} |u|max={speed_hot:.3e} melt={melt_hot:.4} lean={lean_hot:.3e}; sub-critical Ra={ra_cold:.1} |u|max={speed_cold:.3e} melt={melt_cold:.4} lean={lean_cold:.2e}"
    );
    assert_eq!(convection_regime(ra_hot), ConvectionRegime::Convection);
    assert_eq!(convection_regime(ra_cold), ConvectionRegime::Conduction);
    assert!(speed_hot > 20.0 * speed_cold, "super-critical flow should dominate: {speed_hot} vs {speed_cold}");
    assert!(melt_hot > melt_cold + 0.05, "super-critical melt should be ahead: {melt_hot} vs {melt_cold}");
    assert!(lean_hot > lean_cold, "super-critical interface should lean more");
}

// ------------------------------------------------------------- conservation

/// The coupled step conserves the total enthalpy to round-off in a closed
/// (adiabatic) box: conduction, phase change and the conservative advection
/// all telescope. A hot spot melts outward and drives a flow.
#[test]
fn coupled_energy_is_conserved_to_round_off() {
    let (set, k, t_m, solid, _) = gallium_phases();
    let nx = 32;
    let ny = 32;
    let w = 0.02;
    let mut pc = Isobaric2D::new(
        set,
        P_ATM,
        k.to_vec(),
        nx,
        ny,
        w / nx as f64,
        w / ny as f64,
        solid,
        t_m,
        (200.0, 500.0),
    );
    pc.expansion = false;
    for j in 0..ny {
        for i in 0..nx {
            let x = (i as f64 + 0.5) * w / nx as f64;
            let y = (j as f64 + 0.5) * w / ny as f64;
            let t = t_m + 4.0 * (-((x - 0.5 * w).powi(2) + (y - 0.5 * w).powi(2)) / (2.0 * (0.2 * w).powi(2))).exp();
            pc.set_cell(i, j, solid, t);
        }
    }
    let mut flow = gallium_flow(nx, ny, w, w, t_m);
    flow.t_wall = [None, None, None, None];
    let mut sys = EnthalpyConvection2D::new(pc, flow, 1e6);
    let e0 = sys.pc.energy();
    sys.run(20.0, 0.02);
    let e1 = sys.pc.energy();
    println!("coupled energy drift {:.3e} (t={:.2}s)", e1 / e0 - 1.0, sys.pc.t);
    assert!((e1 / e0 - 1.0).abs() < 1e-11, "energy drift {:e}", e1 / e0 - 1.0);
}
