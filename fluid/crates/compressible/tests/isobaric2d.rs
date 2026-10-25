//! Tests for the two-dimensional Eulerian isobaric phase-change model
//! (`compressible::isobaric2d`).
//!
//! Benchmarks:
//!   * the one-phase Stefan/Neumann solution in a 2D domain, wall-driven
//!     (aligned) and as a manufactured solution rotated 45° (grid-orientation
//!     error), with the expansion flow on and off (before/after);
//!   * the exact cylindrical similarity solution for melting outward from a
//!     line source, with grid convergence;
//!   * a radially symmetric **density-change** case: cylindrical freezing from
//!     a cold line, where the solid properties set the exact front and the
//!     expansion flow is the only density effect;
//!   * energy and mass conservation to round-off with adiabatic walls;
//!   * consistency with the one-dimensional Lagrangian model (`isobaric.rs`);
//!   * the Rayleigh-number unfold indicator.

use compressible::elements::{by_symbol, fit, fit_inputs, FitInputs, Grade, P_ATM};
use compressible::isobaric::{enthalpy, Isobaric, PhaseSet};
use compressible::isobaric2d::{
    convection_regime, rayleigh_benard, rayleigh_internal, Bc2D, ConvectionRegime, FlowBc, Isobaric2D,
    RA_CRIT_RAYLEIGH_BENARD,
};
use compressible::phase_change::v_pt;

// ------------------------------------------------------------------ exact math

/// `erf` (Maclaurin + continued-fraction `erfc`), accurate to ~1e-15.
fn erf(x: f64) -> f64 {
    if x.abs() >= 1.5 {
        return if x >= 0.0 { 1.0 - erfc_pos(x) } else { erfc_pos(-x) - 1.0 };
    }
    let (mut term, mut sum, mut n) = (x, x, 0u32);
    loop {
        n += 1;
        term *= -x * x * (2.0 * n as f64 - 1.0) / (n as f64 * (2.0 * n as f64 + 1.0));
        sum += term;
        if term.abs() < 1e-18 * sum.abs() || n > 400 {
            break;
        }
    }
    2.0 / std::f64::consts::PI.sqrt() * sum
}

fn erfc_pos(x: f64) -> f64 {
    const TINY: f64 = 1e-300;
    let mut f = if x == 0.0 { TINY } else { x };
    let (mut c, mut d, mut n) = (f, 0.0, 0u32);
    loop {
        n += 1;
        let an = n as f64 * 0.5;
        d = x + an * d;
        if d == 0.0 {
            d = TINY;
        }
        c = x + an / c;
        if c == 0.0 {
            c = TINY;
        }
        d = 1.0 / d;
        let delta = c * d;
        f *= delta;
        if (delta - 1.0).abs() < 1e-16 || n > 300 {
            break;
        }
    }
    (-x * x).exp() / std::f64::consts::PI.sqrt() / f
}

/// Exponential integral `E1(x) = ∫_x^∞ e^{-t}/t dt` (series below 1,
/// continued fraction above; Numerical Recipes §6.3).
fn e1(x: f64) -> f64 {
    const EULER: f64 = 0.577_215_664_901_532_9;
    assert!(x > 0.0);
    if x < 1.0 {
        let mut sum = -EULER - x.ln();
        let mut term = 1.0;
        for n in 1..=200 {
            term *= -x / n as f64;
            sum -= term / n as f64;
            if term.abs() < 1e-18 {
                break;
            }
        }
        sum
    } else {
        let (mut b, mut c, mut d) = (x + 1.0, 1e30, 1.0 / (x + 1.0));
        let mut h = d;
        for i in 1..=300 {
            let a = -(i as f64) * (i as f64);
            b += 2.0;
            d = 1.0 / (a * d + b);
            c = b + a / c;
            let del = c * d;
            h *= del;
            if (del - 1.0).abs() < 1e-16 {
                break;
            }
        }
        h * (-x).exp()
    }
}

/// Planar one-phase Stefan root `λ e^{λ²} erf λ = St/√π`.
fn neumann_lambda_planar(st: f64) -> f64 {
    let f = |l: f64| l * (l * l).exp() * erf(l) - st / std::f64::consts::PI.sqrt();
    let (mut lo, mut hi) = (1e-9, 0.5);
    while f(hi) < 0.0 {
        hi *= 1.5;
    }
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if f(mid) < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// Line-source Stefan root `λ² e^{λ²} = St_Q`.
///
/// A line source of power per length `Q` melts outward: the interface is
/// `s = 2λ√(αt)` and the liquid temperature is
/// `T = T_m + (Q/4πk)[E1(η²) − E1(λ²)]`, `η = r/(2√(αt))`. The Stefan
/// condition at the cylinder front gives
/// `λ² e^{λ²} = Q/(4π ρ L α) = Q c/(4π k L)` (the source Stefan number).
/// The same transcendental equation holds for cylindrical freezing from a cold
/// line, with the **solid** properties and `Q` the magnitude of the heat
/// removed. Carslaw & Jaeger, *Conduction of Heat in Solids*, 2nd ed., §10.4–11.7.
fn line_source_lambda(st: f64) -> f64 {
    assert!(st > 0.0, "line-source Stefan number must be positive");
    let f = |l: f64| l * l * (l * l).exp() - st;
    let (mut lo, mut hi) = (1e-9, 0.5);
    while f(hi) < 0.0 {
        hi *= 1.5;
    }
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if f(mid) < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

// ------------------------------------------------------------------- materials

/// Water/ice/vapour, fitted to 1-atm data (`substances.rs`); conductivities
/// [ice, liquid, vapour] W/(m K) as used by `examples/phase_change`.
fn water_model() -> (PhaseSet, [f64; 3]) {
    let m = compressible::substances::by_name("water").unwrap().model();
    (PhaseSet::from(&m), [2.22, 0.561, 0.025])
}

/// Iron from the periodic-table fit; conductivities [solid, liquid, vapour].
fn iron_model() -> (PhaseSet, [f64; 3]) {
    let m = fit(by_symbol("Fe").unwrap()).unwrap();
    (PhaseSet::from(&m), [33.0, 35.0, 0.05])
}

/// Iron-like model material with equal densities and negligible expansion, so
/// the constant-property Neumann solution is exact.
fn equal_density_model() -> (PhaseSet, [f64; 3]) {
    let m = 55.845e-3;
    let model = fit_inputs(&FitInputs {
        molar_mass: m,
        t_melt: 1811.0,
        t_boil: 3134.0,
        sublimes: false,
        p_triple: None,
        rho_s: 7000.0,
        t_ref_s: 1811.0,
        alpha_s: 1e-7,
        c_bulk_s: 4000.0,
        shear_speed: 2000.0,
        rho_l: 7000.0,
        t_ref_l: 1811.0,
        alpha_l: 1e-7,
        c_l: 4000.0,
        cp_s: 800.0,
        cp_l: 800.0,
        gamma_v: 5.0 / 3.0,
        cv_v: 1.5 * 8.314 / m,
        h_fus: 247e3,
        h_vap: 6.1e6,
        t_crit: None,
        grade: Grade::C,
        estimated: vec![],
    });
    (PhaseSet::from(&model), [35.0, 35.0, 0.05])
}

/// One-phase Stefan parameters for a transition `old → new` at `dt_wall`.
struct StefanParams {
    tt: f64,
    st: f64,
    lambda: f64,
    alpha: f64,
}

fn stefan_params(set: &PhaseSet, k: &[f64; 3], old: usize, new: usize, dt_wall: f64, t_range: (f64, f64)) -> StefanParams {
    let lad = compressible::isobaric::Ladder::new(set, P_ATM, t_range.0, t_range.1);
    let tt = lad
        .transitions
        .iter()
        .find(|x| (x.1 == old && x.2 == new) || (x.1 == new && x.2 == old))
        .expect("transition")
        .0;
    let latent = (enthalpy(&set.eos[new], P_ATM, tt) - enthalpy(&set.eos[old], P_ATM, tt)).abs();
    let cp = set.cp(new, P_ATM, tt);
    let rho = 1.0 / v_pt(&set.eos[new], P_ATM, tt);
    let alpha = k[new] / (rho * cp);
    let st = cp * dt_wall.abs() / latent;
    StefanParams { tt, st, lambda: neumann_lambda_planar(st), alpha }
}

/// Exact one-phase temperature at depth `x` (from the wall) and time `t`.
fn stefan_exact(p: &StefanParams, x: f64, t: f64, tw: f64) -> f64 {
    if t <= 0.0 {
        return p.tt;
    }
    let eta = x / (2.0 * (p.alpha * t).sqrt());
    tw + (p.tt - tw) * erf(eta) / erf(p.lambda)
}

// ------------------------------------------------------------------- inversion

#[test]
fn volumetric_state_inverts_pure_and_mixture() {
    let (set, k) = water_model();
    let mut s = Isobaric2D::new(set, P_ATM, k.to_vec(), 1, 1, 1.0, 1.0, 0, 260.0, (150.0, 1000.0));
    s.set_cell(0, 0, 0, 250.0);
    let (t, y, _) = s.state_at(0);
    assert!((t - 250.0).abs() < 1e-9, "ice T {t}");
    assert!(y[0] > 0.999);
    s.set_cell(0, 0, 1, 300.0);
    let (t, y, _) = s.state_at(0);
    assert!((t - 300.0).abs() < 1e-9, "water T {t}");
    assert!(y[1] > 0.999);
    // a mixture halfway through the latent heat
    let tt = s.ladder.transitions.iter().find(|x| x.2 == 1).unwrap().0;
    let set2 = {
        let m = compressible::substances::by_name("water").unwrap().model();
        PhaseSet::from(&m)
    };
    let (ha, hb) = (enthalpy(&set2.eos[0], P_ATM, tt), enthalpy(&set2.eos[1], P_ATM, tt));
    let h = 0.5 * (ha + hb);
    let (t, y) = s.state_h(h);
    assert!((t - tt).abs() < 1e-9, "mixture T {t} vs {tt}");
    assert!((y[0] + y[1] - 1.0).abs() < 1e-12);
    assert!(y[0] > 0.4 && y[0] < 0.6, "ice fraction {}", y[0]);
}

/// `h` (and the volumetric `U`) are strictly monotone along the equilibrium
/// path for the condensed phases (this is what makes the inversion well posed).
#[test]
fn volumetric_enthalpy_is_monotone() {
    for (set, k, name, t_range, temps) in [
        (water_model().0, water_model().1, "water", (150.0, 350.0), [200.0, 300.0]),
        (iron_model().0, iron_model().1, "iron", (300.0, 2500.0), [400.0, 1900.0]),
        (equal_density_model().0, equal_density_model().1, "model", (500.0, 3000.0), [800.0, 1900.0]),
    ] {
        let mut s = Isobaric2D::new(set.clone(), P_ATM, k.to_vec(), 1, 1, 1.0, 1.0, 0, t_range.0, t_range);
        let (lo, hi) = (t_range.0 * 1.02, t_range.1 * 0.98);
        let n = 4000;
        let mut prev = f64::NEG_INFINITY;
        for i in 0..=n {
            let t = lo * (hi / lo).powf(i as f64 / n as f64);
            let ph = s.stable_phase(t);
            let u = s.u_pure(ph, t);
            assert!(u > prev, "{name}: U not monotone at T={t} ({u} <= {prev})");
            prev = u;
        }
        for (kk, &tk) in temps.iter().enumerate() {
            s.set_cell(0, 0, kk, tk);
            let (tr, y, _) = s.state_at(0);
            assert!((tr / tk - 1.0).abs() < 1e-9, "{name} phase {kk}: {tr} vs {tk}");
            assert!(y[kk] > 0.999);
        }
    }
}

// --------------------------------------------------------------- conservation

/// With adiabatic walls the total energy `∫U dV` is conserved exactly by the
/// flux update (telescoping) and the total mass `∫ρ dV` is conserved exactly by
/// the conservative expansion mass flux. Iron has a 4.9 % density jump on
/// melting, so the expansion flow is active.
#[test]
fn energy_and_mass_are_conserved_to_round_off() {
    let (set, k) = iron_model();
    let nx = 32;
    let dx = 1e-3 / nx as f64;
    let mut s = Isobaric2D::new(set, P_ATM, k.to_vec(), nx, nx, dx, dx, 0, 1500.0, (300.0, 5000.0));
    for j in 0..nx {
        for i in 0..nx {
            let x = (i as f64 + 0.5) * dx;
            let y = (j as f64 + 0.5) * dx;
            let t = 1500.0 + 600.0 * (-((x - 0.5e-3).powi(2) + (y - 0.5e-3).powi(2)) / (2.0 * (1.5e-4f64).powi(2))).exp();
            s.set_cell(i, j, s.stable_phase(t), t);
        }
    }
    let e0 = s.energy();
    let m0 = s.mass();
    s.advance(0.2);
    let e1 = s.energy();
    let m1 = s.mass();
    assert!(s.t > 0.0);
    println!("energy drift {:+.3e}, mass drift {:+.3e}", e1 / e0 - 1.0, m1 / m0 - 1.0);
    assert!((e1 / e0 - 1.0).abs() < 1e-13, "energy drift {:e}", e1 / e0 - 1.0);
    assert!((m1 / m0 - 1.0).abs() < 1e-13, "mass drift {:e}", m1 / m0 - 1.0);
}

// ------------------------------------------------------------------- aligned

/// Wall-driven one-phase Stefan in a 2D domain, front parallel to the grid.
/// `expansion` selects the expansion-flow closure over the zero-flow one.
#[allow(clippy::too_many_arguments)]
fn run_aligned(
    tag: &str,
    set: PhaseSet,
    k: [f64; 3],
    old: usize,
    new: usize,
    dt_wall: f64,
    t_range: (f64, f64),
    nx: usize,
    ny: usize,
    expansion: bool,
) -> f64 {
    let p = stefan_params(&set, &k, old, new, dt_wall, t_range);
    let t_end = 1.0;
    let s_end = 2.0 * p.lambda * (p.alpha * t_end).sqrt();
    let t0 = p.tt * (1.0 - 1e-9 * (new as f64 - old as f64).signum());
    let lx = 3.0 * s_end;
    let dx = lx / nx as f64;
    let dy = lx / ny as f64;
    let mut s = Isobaric2D::new(set, P_ATM, k.to_vec(), nx, ny, dx, dy, old, t0, t_range);
    s.bc.west = Bc2D::Dirichlet(p.tt + dt_wall);
    s.expansion = expansion;
    s.advance(t_end);
    let area = s.phase_area(new);
    let front = area / (ny as f64 * dy);
    let exact = 2.0 * p.lambda * (p.alpha * s.t).sqrt();
    let err = front / exact - 1.0;
    println!(
        "{tag}: St {:.4}, λ {:.5}, front {:.6} mm (exact {:.6}, {:+.3}%)",
        p.st,
        p.lambda,
        front * 1e3,
        exact * 1e3,
        100.0 * err
    );
    err
}

#[test]
fn stefan_aligned_matches_neumann() {
    // equal-density model material: the constant-property Neumann solution is
    // exact; the expansion flow is negligible (purely thermal) and the model
    // reproduces the solution either way
    let e = run_aligned("model aligned      ", equal_density_model().0, [35.0; 3], 0, 1, 189.0, (500.0, 5000.0), 64, 1, true);
    assert!(e.abs() < 0.02, "aligned model front error {e}");
    // iron melting and water freezing: the density jump drives an expansion
    // flow; with it the fixed grid reproduces the Neumann solution, without it
    // the zero-flow closure was +10.6 % / −14.5 % off.
    let iron_on = run_aligned("iron aligned (on)  ", iron_model().0, [33.0, 35.0, 0.05], 0, 1, 189.0, (500.0, 5000.0), 64, 1, true);
    let iron_off = run_aligned("iron aligned (off) ", iron_model().0, [33.0, 35.0, 0.05], 0, 1, 189.0, (500.0, 5000.0), 64, 1, false);
    let water_on = run_aligned("water freeze (on)  ", water_model().0, [2.22, 0.561, 0.025], 1, 0, -10.0, (150.0, 350.0), 64, 1, true);
    let water_off = run_aligned("water freeze (off) ", water_model().0, [2.22, 0.561, 0.025], 1, 0, -10.0, (150.0, 350.0), 64, 1, false);
    assert!(iron_on.abs() < 0.05, "iron aligned (expansion) {iron_on}");
    assert!(water_on.abs() < 0.05, "water aligned (expansion) {water_on}");
    // the expansion flow removes the zero-flow bias
    assert!(iron_on.abs() < 0.25 * iron_off.abs(), "iron expansion did not close the gap");
    assert!(water_on.abs() < 0.25 * water_off.abs(), "water expansion did not close the gap");
}

#[test]
fn stefan_aligned_with_outflow_matches() {
    // Physically motivated expansion BC: fixed wall on the cold side and an
    // open (pressure-fixed) outflow on the far side, like the free end of the
    // 1D Lagrangian mesh. Also exercises the direct tridiagonal Poisson path.
    let (set, k) = water_model();
    let p = stefan_params(&set, &k, 1, 0, -10.0, (150.0, 350.0));
    let t_end = 1.0;
    let n = 64;
    let s_end = 2.0 * p.lambda * (p.alpha * t_end).sqrt();
    let dx = 3.0 * s_end / n as f64;
    let t0 = p.tt * (1.0 - 1e-9 * (0.0 - 1.0f64).signum());
    let mut s = Isobaric2D::new(set, P_ATM, k.to_vec(), n, 1, dx, dx, 1, t0, (150.0, 350.0));
    s.bc.west = Bc2D::Dirichlet(p.tt - 10.0);
    s.flow_bc.west = FlowBc::Wall;
    s.flow_bc.east = FlowBc::Open;
    s.advance(t_end);
    let front = s.phase_area(0) / dx;
    let exact = 2.0 * p.lambda * (p.alpha * s.t).sqrt();
    println!(
        "water freeze open east: front {:.6} mm (exact {:.6}, {:+.3}%)",
        front * 1e3,
        exact * 1e3,
        100.0 * (front / exact - 1.0)
    );
    assert!((front / exact - 1.0).abs() < 0.05, "open-outflow front error");
}

// ------------------------------------------------------------------- rotated

/// Manufactured one-phase Stefan solution rotated 45°: the exact temperature
/// depends only on `ξ = (x+y)/√2`, and all four boundaries are driven by the
/// exact solution (per-cell `bc_row`). The front is a diagonal staircase, so
/// comparing its error with the aligned case isolates the grid-orientation
/// error. `aligned = true` drives the same machinery with `ξ = x`.
fn run_manufactured(nx: usize, aligned: bool) -> f64 {
    let (set0, k) = equal_density_model();
    let p = stefan_params(&set0, &k, 0, 1, 189.0, (500.0, 5000.0));
    let tw = p.tt + 189.0;
    let l = 1.0;
    let dx = l / nx as f64;
    // start with a front resolved over ~10 cells and let it advance ~3 more
    let (s0, s1) = (0.30, 0.40);
    let t_init = (s0 / (2.0 * p.lambda)).powi(2) / p.alpha;
    let t_end = (s1 / (2.0 * p.lambda)).powi(2) / p.alpha;
    let mut s = Isobaric2D::new(set0, P_ATM, k.to_vec(), nx, nx, dx, dx, 0, tw, (500.0, 5000.0));
    let xi = |x: f64, y: f64| if aligned { x } else { (x + y) / std::f64::consts::SQRT_2 };
    for j in 0..nx {
        for i in 0..nx {
            let (x, y) = ((i as f64 + 0.5) * dx, (j as f64 + 0.5) * dx);
            if xi(x, y) < 2.0 * p.lambda * (p.alpha * t_init).sqrt() {
                s.set_cell(i, j, 0, stefan_exact(&p, xi(x, y), t_init, tw));
            } else {
                s.set_cell(i, j, 1, p.tt);
            }
        }
    }
    let nstep = 40;
    let dtau = (t_end - t_init) / nstep as f64;
    for _ in 0..nstep {
        let tn = s.t + dtau;
        let edge = |x: f64, y: f64| stefan_exact(&p, xi(x, y), tn, tw);
        let west: Vec<f64> = (0..nx).map(|j| edge(0.0, (j as f64 + 0.5) * dx)).collect();
        let east: Vec<f64> = (0..nx).map(|j| edge(l, (j as f64 + 0.5) * dx)).collect();
        let south: Vec<f64> = (0..nx).map(|i| edge((i as f64 + 0.5) * dx, 0.0)).collect();
        let north: Vec<f64> = (0..nx).map(|i| edge((i as f64 + 0.5) * dx, l)).collect();
        s.set_dirichlet_row(0, west);
        s.set_dirichlet_row(1, east);
        s.set_dirichlet_row(2, south);
        s.set_dirichlet_row(3, north);
        s.bc.west = Bc2D::Dirichlet(tw);
        s.bc.east = Bc2D::Dirichlet(tw);
        s.bc.south = Bc2D::Dirichlet(tw);
        s.bc.north = Bc2D::Dirichlet(tw);
        s.advance(dtau);
    }
    // front crossing along the normal
    let temps = s.temperatures();
    let mut front = f64::NAN;
    for i in 0..nx.saturating_sub(1) {
        let (ia, ib) = if aligned { (s.idx(i, 0), s.idx(i + 1, 0)) } else { (s.idx(i, i), s.idx(i + 1, i + 1)) };
        let (a, b) = (temps[ia], temps[ib]);
        if (a - p.tt) * (b - p.tt) <= 0.0 && a != b {
            let f = (p.tt - a) / (b - a);
            let coord = |k: usize| {
                if aligned {
                    (k as f64 + 0.5) * dx
                } else {
                    (k as f64 + 0.5) * dx * std::f64::consts::SQRT_2
                }
            };
            front = coord(i) + f * (coord(i + 1) - coord(i));
            break;
        }
    }
    let exact = 2.0 * p.lambda * (p.alpha * s.t).sqrt();
    println!(
        "manufactured {} nx={nx}: front {:.5} (exact {:.5}, {:+.3}%), t={:.4}",
        if aligned { "aligned " } else { "rotated45" },
        front,
        exact,
        100.0 * (front / exact - 1.0),
        s.t
    );
    front / exact - 1.0
}

#[test]
fn rotated_45_grid_orientation_error() {
    let ea = run_manufactured(32, true);
    let er = run_manufactured(32, false);
    println!("grid-orientation error (rotated − aligned) at nx=32: {:+.3}%", 100.0 * (er - ea));
    assert!(ea.abs() < 0.05, "aligned manufactured {ea}");
    assert!(er.abs() < 0.10, "rotated manufactured {er}");
}

// ------------------------------------------------------------------- cylinder

/// Exact line-source melting: solid at `tm`, liquid
/// `T = tm + (Q/4πk)[E1(η²) − E1(λ²)]`, interface `s = 2λ√(αt)`.
fn line_source_exact(q: f64, k: f64, alpha: f64, tm: f64, lambda: f64, r: f64, t: f64) -> f64 {
    if t <= 0.0 {
        return tm;
    }
    let eta = r / (2.0 * (alpha * t).sqrt());
    tm + q / (4.0 * std::f64::consts::PI * k) * (e1(eta * eta) - e1(lambda * lambda))
}

/// Melt outward from a line source at the origin (equal-density material, so
/// the constant-property similarity solution is exact). Returns front error.
fn run_cylinder(nx: usize) -> f64 {
    let (set, kmodel) = equal_density_model();
    let k = kmodel[1];
    let tm = 1811.0;
    let cp = set.cp(1, P_ATM, tm);
    let latent = (set.enthalpy(1, P_ATM, tm) - set.enthalpy(0, P_ATM, tm)).abs();
    let rho = 1.0 / v_pt(&set.eos[1], P_ATM, tm);
    let alpha = k / (rho * cp);
    let st = 0.3;
    let lambda = line_source_lambda(st);
    let q = st * 4.0 * std::f64::consts::PI * k * latent / cp;
    let t_end = 1.0;
    let s_end = 2.0 * lambda * (alpha * t_end).sqrt();
    let l = 3.0 * s_end;
    let dx = l / nx as f64;
    let t_init = {
        let want = 4.0 * dx;
        (want / (2.0 * lambda)).powi(2) / alpha
    };
    let mut s = Isobaric2D::new(set, P_ATM, kmodel.to_vec(), nx, nx, dx, dx, 0, tm, (500.0, 5000.0));
    let s0 = 2.0 * lambda * (alpha * t_init).sqrt();
    for j in 0..nx {
        for i in 0..nx {
            let x = (i as f64 + 0.5) * dx - 0.5 * l;
            let y = (j as f64 + 0.5) * dx - 0.5 * l;
            let r = (x * x + y * y).sqrt();
            if r < s0 {
                s.set_cell(i, j, 1, line_source_exact(q, k, alpha, tm, lambda, r, t_init));
            } else {
                s.set_cell(i, j, 0, tm);
            }
        }
    }
    // regularised line source in the four cells around the centre
    let c = nx / 2;
    let cell_area = dx * dx;
    for (di, dj) in [(0usize, 0usize), (0, 1), (1, 0), (1, 1)] {
        let i = c - 1 + di;
        let j = c - 1 + dj;
        let idx = s.idx(i, j);
        s.source[idx] = q / (4.0 * cell_area);
    }
    s.advance(t_end - t_init);
    let tmax = s.temperatures().iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let area = s.phase_area(1);
    let front = (area / std::f64::consts::PI).sqrt();
    let exact = 2.0 * lambda * (alpha * s.t).sqrt();
    println!("cylinder nx={nx}: St_Q {st}, λ {lambda:.4}, front {:.6} mm (exact {:.6}, {:+.3}%), T_max {:.0} K", front * 1e3, exact * 1e3, 100.0 * (front / exact - 1.0), tmax);
    front / exact - 1.0
}

#[test]
fn cylindrical_line_source_converges() {
    let e24 = run_cylinder(24);
    let e48 = run_cylinder(48);
    let e96 = run_cylinder(96);
    println!("cylinder |err| 24/48/96: {:.3}% / {:.3}% / {:.3}%", 100.0 * e24.abs(), 100.0 * e48.abs(), 100.0 * e96.abs());
    assert!(e24.abs() < 0.25, "coarse cylinder {e24}");
    assert!(e48.abs() < 0.15, "mid cylinder {e48}");
    assert!(e96.abs() < 0.08, "fine cylinder {e96}");
}

/// Cylindrical freezing from a cold line at the origin: water at `T_m` freezes
/// outward, the solid (ice) inside carries the similarity profile
/// `T = T_m − (Q/4πk_s)[E1(η²) − E1(λ²)]`, and the front is
/// `s = 2λ√(α_s t)`. The liquid stays at `T_m`, so it is isothermal and the
/// **only** density effect is the expansion flow (ice is 9 % lighter); the
/// front is the exact one plus the expansion. `expansion` selects the closure.
/// Returns the front error vs the exact solution.
fn run_cylinder_freezing(nx: usize, expansion: bool) -> f64 {
    let (set, kmodel) = water_model();
    let ks = kmodel[0]; // ice
    let tm = 273.15;
    let cs = set.cp(0, P_ATM, tm);
    let latent = (set.enthalpy(1, P_ATM, tm) - set.enthalpy(0, P_ATM, tm)).abs();
    let rho_s = 1.0 / v_pt(&set.eos[0], P_ATM, tm);
    let alpha = ks / (rho_s * cs);
    let st_q = 0.1;
    let lambda = line_source_lambda(st_q);
    let q = st_q * 4.0 * std::f64::consts::PI * ks * latent / cs; // W/m removed
    let t_end = 1.0;
    let s_end = 2.0 * lambda * (alpha * t_end).sqrt();
    let l = 3.0 * s_end;
    let dx = l / nx as f64;
    let t_init = {
        let want = 4.0 * dx;
        (want / (2.0 * lambda)).powi(2) / alpha
    };
    let mut s = Isobaric2D::new(set, P_ATM, kmodel.to_vec(), nx, nx, dx, dx, 1, tm, (50.0, 1000.0));
    s.expansion = expansion;
    // walls are adiabatic; the cold line is the only sink
    let s0 = 2.0 * lambda * (alpha * t_init).sqrt();
    for j in 0..nx {
        for i in 0..nx {
            let x = (i as f64 + 0.5) * dx - 0.5 * l;
            let y = (j as f64 + 0.5) * dx - 0.5 * l;
            let r = (x * x + y * y).sqrt();
            if r < s0 {
                // exact solid profile; q negative makes the E1 bracket cool
                s.set_cell(i, j, 0, line_source_exact(-q, ks, alpha, tm, lambda, r, t_init));
            } else {
                s.set_cell(i, j, 1, tm);
            }
        }
    }
    let c = nx / 2;
    let cell_area = dx * dx;
    for (di, dj) in [(0usize, 0usize), (0, 1), (1, 0), (1, 1)] {
        let i = c - 1 + di;
        let j = c - 1 + dj;
        let idx = s.idx(i, j);
        s.source[idx] = -q / (4.0 * cell_area);
    }
    s.advance(t_end - t_init);
    let area = s.phase_area(0); // ice
    let front = (area / std::f64::consts::PI).sqrt();
    let exact = 2.0 * lambda * (alpha * s.t).sqrt();
    let err = front / exact - 1.0;
    println!(
        "cylinder freeze nx={nx} expansion={expansion}: front {:.6} mm (exact {:.6}, {:+.3}%), T_min {:.0} K",
        front * 1e3,
        exact * 1e3,
        100.0 * err,
        s.temperatures().iter().cloned().fold(f64::INFINITY, f64::min)
    );
    err
}

#[test]
fn cylindrical_freezing_with_density_change() {
    let on24 = run_cylinder_freezing(24, true);
    let off24 = run_cylinder_freezing(24, false);
    let on48 = run_cylinder_freezing(48, true);
    let off48 = run_cylinder_freezing(48, false);
    println!(
        "cylinder freezing |err| with/without expansion: nx=24 {:.3}%/{:.3}%, nx=48 {:.3}%/{:.3}%",
        100.0 * on24.abs(),
        100.0 * off24.abs(),
        100.0 * on48.abs(),
        100.0 * off48.abs()
    );
    // With the expansion flow the error is first order in the grid (the line
    // sink is regularised over one cell) and converges; the zero-flow closure
    // is biased and does not. The coarse nx=24 expansion overshoots only
    // because the source-regularisation error dominates there.
    assert!(on48.abs() < on24.abs(), "expansion not converging");
    assert!(on48.abs() < 0.5 * off48.abs(), "expansion did not improve the radial density-change case");
    assert!(on48.abs() < 0.05, "radial freezing with expansion {on48}");
}

// ------------------------------------------------------------- 1D consistency

/// Compare the 2D fixed-grid model in a 1D configuration with the 1D
/// Lagrangian model (`isobaric.rs`), freezing water/ice from a cold wall. The
/// Lagrangian mesh carries the expansion flow exactly; the Eulerian model with
/// the expansion flow must now agree with it, unlike the zero-flow closure.
#[test]
fn two_d_matches_one_d_lagrangian() {
    let (set, k) = water_model();
    let (old, new) = (1usize, 0usize);
    let t_range = (150.0, 1000.0);
    let dt_wall = -10.0;
    let p1 = stefan_params(&set, &k, old, new, dt_wall, t_range);
    let t_end = 1.0;
    let s_end = 2.0 * p1.lambda * (p1.alpha * t_end).sqrt();
    let l = 3.0 * s_end;
    let n = 96;
    let dx = l / n as f64;
    let t0 = p1.tt * (1.0 - 1e-9 * (new as f64 - old as f64).signum());

    let mut one = Isobaric::new(set.clone(), P_ATM, k.to_vec(), n, dx, old, t0, t_range);
    one.wall[0] = Some(p1.tt + dt_wall);
    one.advance(t_end);
    let temps = one.temperatures();
    let front1: f64 = (0..n).map(|c| one.m[c] * one.y[c][new] * v_pt(&one.eos()[new], P_ATM, temps[c])).sum();

    let mut two = Isobaric2D::new(set, P_ATM, k.to_vec(), n, 1, dx, dx, old, t0, t_range);
    two.bc.west = Bc2D::Dirichlet(p1.tt + dt_wall);
    two.expansion = true;
    two.advance(t_end);
    let front2 = two.phase_area(new) / dx;

    // zero-flow closure, the previous rung
    let mut off = Isobaric2D::new(water_model().0, P_ATM, k.to_vec(), n, 1, dx, dx, old, t0, t_range);
    off.bc.west = Bc2D::Dirichlet(p1.tt + dt_wall);
    off.expansion = false;
    off.advance(t_end);
    let front_off = off.phase_area(new) / dx;

    let exact = 2.0 * p1.lambda * (p1.alpha * t_end).sqrt();
    println!(
        "1D vs 2D (water freezing): 1D {:.6} mm, 2D-EX {:.6} mm, 2D-zero {:.6} mm, exact {:.6} mm; 2D/1D {:+.3}%, zero/1D {:+.3}%",
        front1 * 1e3,
        front2 * 1e3,
        front_off * 1e3,
        exact * 1e3,
        100.0 * (front2 / front1 - 1.0),
        100.0 * (front_off / front1 - 1.0)
    );
    assert!((front1 / exact - 1.0).abs() < 0.02, "1D front vs exact");
    assert!((front2 / front1 - 1.0).abs() < 0.01, "2D expansion vs 1D");
    assert!((front_off / front1 - 1.0).abs() > 0.05, "zero-flow gap unexpectedly small");
}

/// Equal-density material: the Eulerian 2D and Lagrangian 1D rungs must agree
/// to the discretisation, because the expansion flow vanishes.
#[test]
fn equal_density_2d_matches_1d() {
    let (set, k) = equal_density_model();
    let (old, new) = (0usize, 1usize);
    let t_range = (500.0, 5000.0);
    let dt_wall = 189.0;
    let p1 = stefan_params(&set, &k, old, new, dt_wall, t_range);
    let t_end = 1.0;
    let s_end = 2.0 * p1.lambda * (p1.alpha * t_end).sqrt();
    let l = 3.0 * s_end;
    let n = 96;
    let dx = l / n as f64;
    let t0 = p1.tt * (1.0 - 1e-9 * (new as f64 - old as f64).signum());

    let mut one = Isobaric::new(set.clone(), P_ATM, k.to_vec(), n, dx, old, t0, t_range);
    one.wall[0] = Some(p1.tt + dt_wall);
    one.advance(t_end);
    let temps = one.temperatures();
    let front1: f64 = (0..n).map(|c| one.m[c] * one.y[c][new] * v_pt(&one.eos()[new], P_ATM, temps[c])).sum();

    let mut two = Isobaric2D::new(set, P_ATM, k.to_vec(), n, 1, dx, dx, old, t0, t_range);
    two.bc.west = Bc2D::Dirichlet(p1.tt + dt_wall);
    two.advance(t_end);
    let front2 = two.phase_area(new) / dx;
    println!("equal-density 1D {:.6} mm, 2D {:.6} mm, ratio {:+.4}%", front1 * 1e3, front2 * 1e3, 100.0 * (front2 / front1 - 1.0));
    assert!((front2 / front1 - 1.0).abs() < 0.02);
}

// ------------------------------------------------------------------ indicator

#[test]
fn rayleigh_indicator_switches_at_drazin_reid() {
    let (g, beta, nu, alpha) = (9.81, 2.07e-4, 1e-6, 1.4e-7);
    let l = 1e-2;
    let ra_cold = rayleigh_benard(g, beta, 0.01, l, nu, alpha);
    let ra_hot = rayleigh_benard(g, beta, 1.0, l, nu, alpha);
    println!("Ra(0.01 K) = {ra_cold:.1}, Ra(1 K) = {ra_hot:.1}, Ra_c = {RA_CRIT_RAYLEIGH_BENARD}");
    assert!(ra_cold < RA_CRIT_RAYLEIGH_BENARD);
    assert_eq!(convection_regime(ra_cold), ConvectionRegime::Conduction);
    assert!(ra_hot > RA_CRIT_RAYLEIGH_BENARD);
    assert_eq!(convection_regime(ra_hot), ConvectionRegime::Convection);
    let ra_i = rayleigh_internal(g, beta, 1e4, l, 0.6, nu, alpha);
    println!("Ra_I = {ra_i:.3e} (critical 2772)");
    assert!(ra_i > 2772.0);
    let (set, k) = water_model();
    let mut s = Isobaric2D::new(set, P_ATM, k.to_vec(), 2, 1, 1e-3, 1e-3, 1, 300.0, (150.0, 1000.0));
    s.set_cell(1, 0, 1, 310.0);
    let b = s.buoyancy(9.81, 2.07e-4, 300.0);
    assert!(b[0].abs() < 1e-12);
    assert!(b[1] > 0.0);
}

/// The bcc phase of iron carries the CALPHAD magnetic Gibbs energy, so its
/// `h(T)` is nonlinear and the inversion takes the safeguarded Newton path
/// (not the NASG closed form). Check the round trip and monotonicity there.
#[test]
fn magnetic_phase_inversion_round_trips() {
    let (set, _) = compressible::allotropes::iron();
    let k = vec![33.0, 34.0, 35.0, 0.05];
    let mut s = Isobaric2D::new(set, P_ATM, k, 1, 1, 1.0, 1.0, 0, 300.0, (300.0, 5000.0));
    assert!(s.model.magnetic[0].is_some(), "bcc should carry magnetism");
    let mut prev = f64::NEG_INFINITY;
    for t in [300.0, 500.0, 900.0, 1040.0, 1043.0, 1100.0, 1180.0] {
        s.set_cell(0, 0, 0, t);
        let (tr, y) = s.state_h(s.u[0] / s.rho[0]);
        assert!((tr / t - 1.0).abs() < 1e-9, "magnetic bcc {tr} vs {t}");
        assert!(y[0] > 0.999);
        assert!(s.u[0] > prev, "U not monotone at {t}");
        prev = s.u[0];
    }
}
