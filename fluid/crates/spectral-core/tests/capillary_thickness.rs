//! Finite-interface-thickness correction to the capillary-wave dispersion.
//!
//! The residual capillary-frequency error is measured from the *initial
//! restoring acceleration* of a single mode (a few BDF1 steps at u = 0,
//! Richardson-extrapolated in dt).  This avoids the transient-fitting
//! ambiguity of a long time series and is a deterministic functional of the
//! initial condition and the discrete operators.  With ρ₁ = ρ₂ the exact
//! sharp-interface value is ω² = σk³/(ρ₁+ρ₂).
//!
//! Linearising the Cahn–Hilliard capillary force `μ_c∇φ` about the tanh
//! profile gives a force density ∝ φ₀'(z)², so the restoring stiffness is
//! reduced by the flow-weighted overlap
//!
//!   ω²_diffuse / ω²_sharp = ∫ e^{−k|z|} φ₀'(z)² dz / ∫ φ₀'(z)² dz
//!                         = 1 − 0.6267 (ηk) + O((ηk)²),
//!
//! the thin-diffuse-interface correction (Jacqmin 1999, JCP 155, 96; derived
//! in `analysis/twophase_and_mhd.md`).  The tests below check that the
//! measured correction is the physical finite-thickness one — it is invariant
//! under refinement at fixed physical η and k, and its magnitude matches the
//! overlap integral — rather than a projection/dealiasing inconsistency.
use spectral_core::twophase::{Projection, TwoPhaseParams, TwoPhaseSolver};
use spectral_core::{init, Grid};

const PI: f64 = std::f64::consts::PI;

/// tanh interface of physical thickness `eta`, single mode `a0 cos(k x)`.
fn build(n: usize, eta: f64, k: f64, a0: f64) -> TwoPhaseSolver {
    let g = Grid::new(n);
    let prm = TwoPhaseParams {
        rho: [1.0, 1.0],
        mu: [0.0, 0.0],
        sigma: 1.0,
        gravity: [0.0; 3],
        eta,
        mobility: 1e-8,
        stab: 2.0,
    };
    let w = 2f64.sqrt() * eta;
    let (zlo, zhi) = (0.5 * PI, 1.5 * PI);
    let phi = init::from_fn(&g, |x, _, z| {
        let zeta = a0 * (k * x).cos();
        [((z - zlo - zeta) / w).tanh() - ((z - zhi - zeta) / w).tanh() - 1.0, 0.0, 0.0]
    })[0]
        .clone();
    let mut s = TwoPhaseSolver::new(g.clone(), prm, init::from_fn(&g, |_, _, _| [0.0; 3]), phi);
    s.set_projection(Projection::Constant);
    s
}

/// Reference flat profile and its z-derivative.
fn reference(g: &Grid, eta: f64, zlo: f64, zhi: f64) -> (Vec<f64>, Vec<f64>) {
    let w = 2f64.sqrt() * eta;
    let n = g.n;
    let mut phi0 = vec![0.0; g.real_len()];
    let mut g0 = vec![0.0; g.real_len()];
    for i in 0..n {
        for j in 0..n {
            for l in 0..n {
                let z = g.x(l);
                let a = ((z - zlo) / w).tanh();
                let b = ((z - zhi) / w).tanh();
                phi0[(i * n + j) * n + l] = a - b - 1.0;
                let sa = 1.0 / ((z - zlo) / w).cosh().powi(2);
                let sb = 1.0 / ((z - zhi) / w).cosh().powi(2);
                g0[(i * n + j) * n + l] = (sa - sb) / w;
            }
        }
    }
    (phi0, g0)
}

/// Interface-shift projection ζ = −∫δφ ∂zφ₀ dz / ∫(∂zφ₀)² dz, at x = 0 and
/// x = π/2, combined so the k=0 bulk drift cancels.
fn proj_amp(s: &TwoPhaseSolver, phi0: &[f64], g0: &[f64]) -> f64 {
    let g = &s.grid;
    let n = g.n;
    let phi = s.phi_phys();
    let norm: f64 = (0..n).map(|l| g0[l] * g0[l]).sum();
    let zeta = |i: usize| {
        let base = (i * n) * n;
        let mut acc = 0.0;
        for l in 0..n {
            acc += (phi[base + l] - phi0[base + l]) * g0[base + l];
        }
        -acc / norm
    };
    0.5 * (zeta(0) - zeta(n / 4))
}

/// ω² = −ζ̈(0)/ζ(0) from a 3-parameter least-squares fit over `steps` BDF1
/// steps (still in the linear regime at a₀ = 1e-4).
fn accel_omega2(n: usize, eta: f64, k: f64, dt: f64, steps: usize) -> f64 {
    let a0 = 1e-4;
    let mut s = build(n, eta, k, a0);
    let (zlo, zhi) = (0.5 * PI, 1.5 * PI);
    let (phi0, g0) = reference(&s.grid, eta, zlo, zhi);
    let z0 = proj_amp(&s, &phi0, &g0);
    let mut ts = vec![0.0];
    let mut zs = vec![z0];
    for i in 0..steps {
        s.step(dt);
        ts.push((i + 1) as f64 * dt);
        zs.push(proj_amp(&s, &phi0, &g0));
    }
    let mut m = [[0.0f64; 3]; 3];
    let mut rhs = [0.0f64; 3];
    for (&t, &z) in ts.iter().zip(zs.iter()) {
        let phi = [1.0, t, 0.5 * t * t];
        for i in 0..3 {
            for j in 0..3 {
                m[i][j] += phi[i] * phi[j];
            }
            rhs[i] += phi[i] * (z - z0);
        }
    }
    let (mut a, mut b) = (m, rhs);
    for i in 0..3 {
        let p = (i..3).max_by(|&x, &y| a[x][i].abs().partial_cmp(&a[y][i].abs()).unwrap()).unwrap();
        a.swap(i, p);
        b.swap(i, p);
        for j in i + 1..3 {
            let f = a[j][i] / a[i][i];
            for kk in i..3 {
                a[j][kk] -= f * a[i][kk];
            }
            b[j] -= f * b[i];
        }
    }
    let mut x = [0.0; 3];
    for i in (0..3).rev() {
        let mut sum = b[i];
        for j in i + 1..3 {
            sum -= a[i][j] * x[j];
        }
        x[i] = sum / a[i][i];
    }
    -x[2] / z0
}

/// Richardson (first order in dt) acceleration estimate.
fn omega2(n: usize, eta: f64, k: f64) -> f64 {
    let coarse = accel_omega2(n, eta, k, 2e-3, 8);
    let fine = accel_omega2(n, eta, k, 1e-3, 8);
    2.0 * fine - coarse
}

/// ∫_0^∞ e^{−a u} sech⁴u du by Simpson, used for the thin-interface overlap.
fn overlap(a: f64) -> f64 {
    let n = 200_000usize;
    let h = 40.0 / n as f64;
    let f = |u: f64| (-a * u).exp() / u.cosh().powi(4);
    let mut s = f(0.0) + f(40.0);
    for i in 1..n {
        let w = if i % 2 == 0 { 2.0 } else { 4.0 };
        s += w * f(i as f64 * h);
    }
    s * h / 3.0
}

/// Analytic diffuse-interface correction to ω² at finite ηk:
/// (3/2)∫ e^{−√2(ηk)u} sech⁴u du, which reduces to the sharp value 1 at ηk=0
/// and has slope −0.6267.
fn theory_ratio(eta_k: f64) -> f64 {
    1.5 * overlap(2f64.sqrt() * eta_k)
}

#[test]
fn restoring_stiffness_is_finite_thickness_not_grid() {
    // Fixed physical η and k (ηk = 0.6283), so only η/Δx changes: N = 24 has
    // 1.2 points per η, N = 48 has 2.4.  A discretisation inconsistency tied
    // to the interface resolution would move the result; a physical
    // finite-thickness correction does not.
    let k = 2.0f64;
    let eta = 2.0 * PI / 20.0; // = 0.31416
    let exact2 = k.powi(3) / 2.0;
    let w2_24 = omega2(24, eta, k);
    let w2_48 = omega2(48, eta, k);
    let e24 = w2_24 / exact2 - 1.0;
    let e48 = w2_48 / exact2 - 1.0;
    eprintln!(
        "fixed eta k={:.4}: N=24 eta/dx={:.2} errW2={:+.3}%, N=48 eta/dx={:.2} errW2={:+.3}%",
        eta * k,
        eta / (2.0 * PI / 24.0),
        100.0 * e24,
        eta / (2.0 * PI / 48.0),
        100.0 * e48
    );
    assert!(
        (e24 - e48).abs() < 0.02,
        "restoring stiffness changed by {:.2}% under refinement at fixed physical eta and k",
        100.0 * (e24 - e48).abs()
    );
}

#[test]
fn diffuse_correction_matches_thin_interface_overlap() {
    // A resolved interface (η = 1.2Δx at N = 64), ηk = 0.2356.  The measured
    // reduction of ω² must equal the force-density overlap integral to a few
    // percent; this pins the correction as the physical one.
    let n = 64;
    let etas = 1.2;
    let k = 2.0f64;
    let dx = 2.0 * PI / n as f64;
    let eta = etas * dx;
    let exact2 = k.powi(3) / 2.0;
    let w2 = omega2(n, eta, k);
    let meas = w2 / exact2 - 1.0;
    let theory = theory_ratio(eta * k) - 1.0;
    eprintln!(
        "N={n} eta/dx={etas} eta k={:.4}: measured errW2={:+.3}% theory={:+.3}%, difference {:+.3}%",
        eta * k,
        100.0 * meas,
        100.0 * theory,
        100.0 * (meas - theory)
    );
    assert!(
        (meas - theory).abs() < 0.04,
        "measured finite-thickness correction {:.2}% disagrees with the overlap integral {:.2}%",
        100.0 * meas,
        100.0 * theory
    );
}
