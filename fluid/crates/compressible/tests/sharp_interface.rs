//! Validation of the sharp solid–fluid interface module.
//!
//! Every assertion states a measured bound; the numbers and the convergence
//! tables are in `analysis/sharp_interface.md`.

use compressible::eos::Nasg;
use compressible::exact;
use compressible::sharp_interface::*;
use phononic::layered::{Layer, Stack};

fn rel(a: f64, b: f64) -> f64 {
    (a - b).abs() / b.abs().max(1e-300)
}

// ---------------------------------------------------------------------
// Exact two-material Riemann solver
// ---------------------------------------------------------------------

/// The exact solid–fluid interface state must reduce to the acoustic
/// (impedance) result as the wave strength vanishes.
#[test]
fn exact_star_reduces_to_acoustic() {
    let w = water();
    let s = steel();
    let wr = State { rho: 7850.0, u: 0.0, p: P0 };
    for amp in [1.0, 1e2, 1e4, 1e6] {
        let l = State { rho: 998.2, u: 0.0, p: P0 + amp };
        let ex = exact_star(&w, &l, &s, &wr);
        let ac = acoustic_star(&w, &l, &s, &wr);
        assert!((ex.u - ac.u).abs() < 1e-3 * ac.u.abs().max(1e-9), "u: {} vs {} at amp {amp}", ex.u, ac.u);
        assert!((ex.tau - ac.tau).abs() < 2e-3 * amp.max(1.0), "tau: {} vs {} at amp {amp}", ex.tau, ac.tau);
    }
}

/// A stiffened gas on both sides: the module's exact solver must agree with
/// the crate's independent `exact::star`.
#[test]
fn exact_star_matches_crate_for_two_gases() {
    let g = Nasg::ideal(1.4);
    let m = Material { kind: Kind::Fluid, eos: g, rho0: 1.0, cs: 0.0 };
    let cases = [
        ((1.0, 0.0, 1.0), (0.125, 0.0, 0.1)),
        ((1.0, 0.0, 1000.0), (1.0, 0.0, 0.01)),
        ((1.0, -2.0, 0.4), (1.0, 2.0, 0.4)),
        ((1.0, 0.0, 0.01), (1.0, 0.0, 100.0)),
        ((5.99924, 19.5975, 460.894), (5.99242, -6.19633, 46.0950)),
    ];
    let mut max_tau = 0.0f64;
    let mut max_u = 0.0f64;
    let mut max_rho = 0.0f64;
    for ((rl, ul, pl), (rr, ur, pr)) in cases {
        let wl = State { rho: rl, u: ul, p: pl };
        let wr = State { rho: rr, u: ur, p: pr };
        let mine = exact_star(&m, &wl, &m, &wr);
        let sl = exact::Side { rho: rl, u: ul, p: pl, eos: g };
        let sr = exact::Side { rho: rr, u: ur, p: pr, eos: g };
        let theirs = exact::star(&sl, &sr);
        max_tau = max_tau.max(rel(mine.tau, theirs.p));
        max_u = max_u.max((mine.u - theirs.u).abs());
        max_rho = max_rho.max(rel(mine.rho_l, theirs.rho_l)).max(rel(mine.rho_r, theirs.rho_r));
        assert!(rel(mine.tau, theirs.p) < 1e-9, "tau {} vs {}", mine.tau, theirs.p);
        assert!((mine.u - theirs.u).abs() < 1e-12, "u {} vs {}", mine.u, theirs.u);
        assert!(rel(mine.rho_l, theirs.rho_l) < 1e-8, "rho_l {} vs {}", mine.rho_l, theirs.rho_l);
        assert!(rel(mine.rho_r, theirs.rho_r) < 1e-8, "rho_r {} vs {}", mine.rho_r, theirs.rho_r);
    }
    if std::env::var("SHARP_REPORT").is_ok() {
        eprintln!("riemann_vs_crate max_rel_tau={max_tau:.3e} max_abs_u={max_u:.3e} max_rel_rho={max_rho:.3e}");
    }
}

/// The GPR 1D solid reproduces the handbook longitudinal speed and has zero
/// elastic stress at the reference density.
#[test]
fn solid_longitudinal_speed() {
    let s = steel();
    assert!(rel(s.c_longitudinal(), 5960.0) < 1e-12, "{}", s.c_longitudinal());
    assert!((s.sound_speed(7850.0, P0) / 5960.0 - 1.0).abs() < 1e-12);
    assert!(s.elastic_stress(7850.0).abs() < 1e-9);
}

// ---------------------------------------------------------------------
// Eulerian high-order solver
// ---------------------------------------------------------------------

/// A stationary water/steel interface at uniform pressure must stay uniform.
#[test]
fn stationary_contact_stays_uniform() {
    let w = water();
    let s = steel();
    let n = 200;
    let mut sol = SharpSolver::new(n, 3, 0.0, 1.0, [w, s], |x| {
        let m = if x < 0.5 { w } else { s };
        (m.to_cons(State { rho: m.rho0, u: 0.0, p: P0 }), x - 0.5)
    });
    sol.cfl = 0.4;
    let dt = sol.dt_cfl();
    for _ in 0..400 {
        sol.step(dt);
    }
    let mut dev = 0.0f64;
    for i in 0..n {
        dev = dev.max(((sol.prim(i).p - P0) / P0).abs());
    }
    // stiff EOS round-off (~γ p∞/p₀ ~ 2.6e4) accumulates in the conservative
    // pressure inversion over 400 steps.
    assert!(dev < 1e-6, "max |dp|/p = {dev:.3e}");
}

/// The finite-volume update is a telescoping sum of single-valued face fluxes.
/// A velocity pulse at uniform pressure has no boundary force and stays
/// inside the domain, so all three conserved totals are unchanged.
#[test]
fn eulerian_conservation_is_exact() {
    let w = water();
    let s = steel();
    let n = 300;
    let mut sol = SharpSolver::new(n, 3, 0.0, 1.0, [w, s], |x| {
        let m = if x < 0.5 { w } else { s };
        let u = 5.0 * (-((x - 0.2) / 0.03).powi(2)).exp();
        (m.to_cons(State { rho: m.rho0, u, p: P0 }), x - 0.5)
    });
    sol.cfl = 0.3;
    let before = sol.totals();
    let dt = sol.dt_cfl();
    for _ in 0..300 {
        sol.step(dt);
    }
    let after = sol.totals();
    for k in 0..3 {
        let scale = before[k].abs().max(1e-300);
        assert!((after[k] - before[k]).abs() / scale < 1e-8, "component {k}: {} vs {}", after[k], before[k]);
    }
}

// ---------------------------------------------------------------------
// Lagrangian sharp-interface solver
// ---------------------------------------------------------------------

/// Lagrangian stationary water/steel contact at uniform pressure.
#[test]
fn lagrangian_stationary_contact() {
    let w = water();
    let s = steel();
    let mut sol = LagrangianSolver::new(200, 0.0, 1.0, [w, s], |x| {
        let m = if x < 0.5 { w } else { s };
        (State { rho: m.rho0, u: 0.0, p: P0 }, if x < 0.5 { 0 } else { 1 })
    });
    sol.cfl = 0.4;
    let dt = sol.dt_cfl();
    for _ in 0..400 {
        sol.step(dt);
    }
    let mut dev = 0.0f64;
    for j in 0..sol.n {
        dev = dev.max(((sol.prim(j).p - P0) / P0).abs());
    }
    assert!(dev < 1e-6, "max |dp|/p = {dev:.3e}");
}

/// With reflecting walls the mass is exactly constant and the energy is
/// exactly constant (the walls do no work); the momentum change equals the
/// integrated wall traction.
#[test]
fn lagrangian_conservation() {
    let w = water();
    let s = steel();
    let mut sol = LagrangianSolver::new(300, 0.0, 1.0, [w, s], |x| {
        let m = if x < 0.5 { w } else { s };
        let u = if x < 0.25 { 20.0 } else { 0.0 };
        (State { rho: m.rho0, u, p: P0 }, if x < 0.5 { 0 } else { 1 })
    });
    sol.cfl = 0.3;
    let before = sol.totals();
    for _ in 0..400 {
        let dt = sol.dt_cfl();
        sol.step(dt);
    }
    let after = sol.totals();
    assert!(rel(after[0], before[0]) < 1e-13, "mass {} vs {}", after[0], before[0]);
    assert!(rel(after[2], before[2]) < 1e-11, "energy {} vs {}", after[2], before[2]);
    // The solver records the Heun-averaged wall tractions; the momentum change
    // must equal the left minus right boundary impulse exactly.
    let expected = before[1] + sol.wall_impulse[0] - sol.wall_impulse[1];
    assert!(rel(after[1], expected) < 1e-10, "momentum {} vs accounting {}", after[1], expected);
}

/// Lagrangian strong water shock onto steel against the exact two-material
/// Riemann solution.
#[test]
fn lagrangian_strong_shock_matches_exact() {
    let w = water();
    let s = steel();
    let wl = State { rho: w.rho0, u: 0.0, p: 1.0e9 };
    let wr = State { rho: s.rho0, u: 0.0, p: P0 };
    let star = exact_star(&w, &wl, &s, &wr);
    let n = 1600;
    let mut sol = LagrangianSolver::new(n, -1.0, 1.0, [w, s], |x| {
        if x < 0.0 { (wl, 0) } else { (wr, 1) }
    });
    sol.cfl = 0.3;
    let t_end = 5.0e-5;
    sol.run_to(t_end);

    let (mut e_p, mut e_u, mut den) = (0.0, 0.0, 0.0);
    for j in (n / 4)..(3 * n / 4) {
        let x = 0.5 * (sol.x[j] + sol.x[j + 1]);
        let num = sol.prim(j);
        let ex = sample_exact(&w, &wl, &s, &wr, &star, x / sol.t);
        e_p += (num.p - ex.p).abs();
        e_u += (num.u - ex.u).abs();
        den += ex.p.abs().max(P0);
    }
    let l1_p = e_p / den;
    let l1_u = e_u / (n as f64 / 2.0);
    if std::env::var("SHARP_REPORT").is_ok() {
        eprintln!("shock star tau={:.5e} u={:.5e} rho_l={:.5e} rho_r={:.5e} L1p={l1_p:.3e} L1u={l1_u:.3e}", star.tau, star.u, star.rho_l, star.rho_r);
    }
    assert!(l1_p < 2e-2, "relative L1 pressure error {l1_p:.3e}");
    assert!(l1_u < 2.0, "mean velocity error {l1_u:.3e} m/s");

    // exact interface stress at the star region
    let j = sol.cell_at(0.5 * star.u * t_end);
    let num = sol.prim(j);
    let tau = sol.mats[sol.mat[j]].tau(num.rho, num.p);
    assert!((tau / star.tau - 1.0).abs() < 1e-2, "star stress {:.4e} vs {:.4e}", tau, star.tau);
}

// ---------------------------------------------------------------------
// Acoustic plate benchmark against the exact transfer matrix
// ---------------------------------------------------------------------

const NU: usize = 8192;

fn fft(re: &mut [f64], im: &mut [f64]) {
    let n = re.len();
    let mut j = 0usize;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if i < j {
            re.swap(i, j);
            im.swap(i, j);
        }
    }
    let mut len = 2;
    while len <= n {
        let ang = -2.0 * std::f64::consts::PI / len as f64;
        let (wr, wi) = (ang.cos(), ang.sin());
        let mut i = 0;
        while i < n {
            let (mut wlr, mut wli) = (1.0f64, 0.0f64);
            for k in 0..len / 2 {
                let ur = re[i + k];
                let ui = im[i + k];
                let vr = re[i + k + len / 2] * wlr - im[i + k + len / 2] * wli;
                let vi = re[i + k + len / 2] * wli + im[i + k + len / 2] * wlr;
                re[i + k] = ur + vr;
                im[i + k] = ui + vi;
                re[i + k + len / 2] = ur - vr;
                im[i + k + len / 2] = ui - vi;
                let nwr = wlr * wr - wli * wi;
                wli = wlr * wi + wli * wr;
                wlr = nwr;
            }
            i += len;
        }
        len <<= 1;
    }
}

fn spectrum(x: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let mut re = x.to_vec();
    let mut im = vec![0.0; x.len()];
    fft(&mut re, &mut im);
    (re, im)
}

fn resample(ts: &[f64], vals: &[f64], nu: usize, t_end: f64) -> Vec<f64> {
    let dtu = t_end / nu as f64;
    let mut out = vec![0.0; nu];
    let mut j = 0usize;
    for m in 0..nu {
        let t = m as f64 * dtu;
        while j + 1 < ts.len() && ts[j + 1] < t {
            j += 1;
        }
        out[m] = if ts.is_empty() {
            0.0
        } else if j + 1 >= ts.len() {
            *vals.last().unwrap()
        } else {
            let f = (t - ts[j]) / (ts[j + 1] - ts[j]);
            vals[j] * (1.0 - f) + vals[j + 1] * f
        };
    }
    out
}

/// Geometry of the plate run: water/steel plate in water, pulse launched in
/// water, probes just outside the plate.
struct Plate {
    n: usize,
    lo: f64,
    hi: f64,
    x1: f64,
    x2: f64,
    xc: f64,
    width: f64,
    f0: f64,
    amp: f64,
    t_end: f64,
}

impl Plate {
    fn standard(n: usize) -> Self {
        Plate { n, lo: 0.0, hi: 0.25, x1: 0.10, x2: 0.12, xc: 0.02, width: 0.006, f0: 80_000.0, amp: 100.0, t_end: 1.6e-4 }
    }
    fn dx(&self) -> f64 {
        (self.hi - self.lo) / self.n as f64
    }
    fn front_probe_x(&self) -> f64 {
        self.x1 - self.dx()
    }
    fn back_probe_x(&self) -> f64 {
        self.x2 + self.dx()
    }
    /// Actual cell-centre positions of the two probe cells. This mirrors
    /// `LagrangianSolver::cell_at` exactly (same tie-breaking), so the exact
    /// stack uses the true probe distances rather than the nominal `dx`
    /// spacers; a half-cell placement mismatch would otherwise show up as a
    /// non-converging phase error.
    fn probe_cells(&self) -> (f64, f64) {
        let n = self.n;
        let dx = self.dx();
        let center = |j: usize| self.lo + (j as f64 + 0.5) * dx;
        let nearest = |x: f64| {
            let mut best = 0usize;
            let mut d = f64::INFINITY;
            for j in 0..n {
                if (center(j) - x).abs() < d {
                    d = (center(j) - x).abs();
                    best = j;
                }
            }
            center(best)
        };
        (nearest(self.front_probe_x()), nearest(self.back_probe_x()))
    }
}

/// Run the Lagrangian plate experiment and record pressure at the two probes.
fn plate_run(p: &Plate, with_plate: bool, dt_override: Option<f64>) -> (Vec<f64>, Vec<f64>, Vec<f64>, f64, usize) {
    let w = water();
    let s = steel();
    let plx = p.front_probe_x();
    let prx = p.back_probe_x();
    let plate = (p.x1, p.x2);
    let mut sol = LagrangianSolver::new(p.n, p.lo, p.hi, [w, s], |x| {
        let mat = if with_plate && x >= plate.0 && x < plate.1 { 1 } else { 0 };
        let m = if mat == 1 { s } else { w };
        let c = w.c_longitudinal();
        let r = x - p.xc;
        let env = (-(r / p.width).powi(2)).exp();
        let dp = p.amp * env * (2.0 * std::f64::consts::PI * p.f0 * r / c).sin();
        // One-way (right-going) acoustic initial condition: p' = Z u' and
        // rho' = p'/c^2, so no spurious left-going half-pulse is launched.
        (State { rho: m.rho0 + dp / (c * c), u: dp / (w.rho0 * c), p: P0 + dp }, mat)
    });
    sol.cfl = 0.4;
    sol.acoustic_interface = true;
    let jl = sol.cell_at(plx);
    let jr = sol.cell_at(prx);
    let dt = dt_override.unwrap_or_else(|| sol.dt_cfl());
    let mut t = Vec::new();
    let mut pl = Vec::new();
    let mut pr = Vec::new();
    let mut steps = 0;
    while sol.t < p.t_end - 1e-15 {
        let d = dt.min(p.t_end - sol.t);
        sol.step(d);
        t.push(sol.t);
        pl.push(sol.prim(jl).p - P0);
        pr.push(sol.prim(jr).p - P0);
        steps += 1;
    }
    (t, pl, pr, dt, steps)
}

/// Exact pressure reflection/transmission of the probe-to-probe stack.
fn exact_rt(p: &Plate, f: f64) -> (f64, f64) {
    let z = water().rho0 * water().c_longitudinal();
    let (pfl, pbl) = p.probe_cells();
    let df = p.x1 - pfl; // probe - plate front (positive)
    let db = pbl - p.x2; // plate back - probe (positive)
    let layers = vec![
        Layer::lossless(water().rho0, water().c_longitudinal(), df.abs().max(1e-9)),
        Layer::lossless(steel().rho0, steel().c_longitudinal(), p.x2 - p.x1),
        Layer::lossless(water().rho0, water().c_longitudinal(), db.abs().max(1e-9)),
    ];
    let sol = Stack::new(z, layers, z).solve(2.0 * std::f64::consts::PI * f);
    ((sol.r.re.powi(2) + sol.r.im.powi(2)).sqrt(), (sol.t.re.powi(2) + sol.t.im.powi(2)).sqrt())
}

/// Measured vs exact plate transfer: (|T| abs error, |R| abs error,
/// energy-identity residual, steps, dt).
fn plate_errors(p: &Plate, f_lo: f64, f_hi: f64) -> (f64, f64, f64, usize, f64) {
    if std::env::var("SHARP_REPORT").is_ok() {
        let (a, b) = p.probe_cells();
        eprintln!(
            "probe n={} pfl={:.7} pbl={:.7} df={:.6e} db={:.6e}",
            p.n,
            a,
            b,
            p.x1 - a,
            b - p.x2
        );
    }
    let with = plate_run(p, true, None);
    let without = plate_run(p, false, Some(with.3));
    let t_end = p.t_end;
    let dtu = t_end / NU as f64;
    let inc = spectrum(&resample(&without.0, &without.1, NU, t_end));
    let tot = spectrum(&resample(&with.0, &with.1, NU, t_end));
    let tra = spectrum(&resample(&with.0, &with.2, NU, t_end));
    let df = 1.0 / (NU as f64 * dtu);
    let k0 = (f_lo / df).round() as usize;
    let k1 = (f_hi / df).round() as usize;
    let (mut num_t, mut num_r, mut ident, mut den): (f64, f64, f64, f64) = (0.0, 0.0, 0.0, 0.0);
    for k in k0..=k1 {
        let f = k as f64 * df;
        let (ir, ii) = (inc.0[k], inc.1[k]);
        let mag2 = ir * ir + ii * ii;
        if mag2 < 1e-24 {
            continue;
        }
        let tr = ((tra.0[k].powi(2) + tra.1[k].powi(2)) / mag2).sqrt();
        let rr = (((tot.0[k] - ir).powi(2) + (tot.1[k] - ii).powi(2)) / mag2).sqrt();
        let (re, te) = exact_rt(p, f);
        // Weight by the incident power so frequencies the pulse does not
        // excite (where the ratio is numerical noise) do not dominate.
        num_t += mag2 * (tr - te).powi(2);
        num_r += mag2 * (rr - re).powi(2);
        ident += mag2 * (rr * rr + tr * tr - 1.0).abs();
        den += mag2;
    }
    let n = den.max(1e-300);
    ((num_t / n).sqrt(), (num_r / n).sqrt(), ident / n, with.4, with.3)
}

/// The sharp plate transfer must beat the diffuse-interface absolute |T|
/// error (6.2e-2 at 1024 cells in `analysis/fsi.md`) and converge with
/// resolution.
#[test]
fn acoustic_plate_transfer_converges() {
    let mut prev = f64::INFINITY;
    let mut errs = Vec::new();
    for n in [200usize, 400, 800] {
        let p = Plate::standard(n);
        let cells = (p.x2 - p.x1) / p.dx();
        let (et, er, ident, steps, dt) = plate_errors(&p, 40_000.0, 120_000.0);
        if std::env::var("SHARP_REPORT").is_ok() {
            eprintln!("plate n={n} cells/plate={cells:.0} |T|abs={et:.4e} |R|abs={er:.4e} identity={ident:.4e} steps={steps} dt={dt:.3e}");
        }
        errs.push((n, cells, et, er, ident));
        assert!(et < prev, "|T| absolute error did not fall: {errs:?}");
        prev = et;
    }
    let finest = errs.last().unwrap().2;
    // Absolute incident-weighted |T| error; the plate is 16/32/64 cells thick.
    assert!(finest < 1.0e-2, "|T| absolute error {finest:.3e} at the finest grid");
    assert!(errs.last().unwrap().3 < 0.1, "|R| absolute error {:.3e}", errs.last().unwrap().3);
    assert!(errs.last().unwrap().4 < 0.1, "energy identity residual {:.3e}", errs.last().unwrap().4);
}
