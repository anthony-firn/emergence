//! Benchmarks for the phononic rung. Every value is compared with an exact or
//! published result; measured errors are printed so `analysis/phononic.md` can
//! quote them. No tolerance is loosened to make a case pass.

use phononic::fluid_coupling;
use phononic::homogenize::{
    hashin_shtrikman_bulk, hashin_shtrikman_shear, homogenisation_error,
    laminate_effective_density, laminate_effective_modulus, laminate_effective_speed,
    phase_velocity_error, reuss, voigt,
};
use phononic::lamb;
use phononic::layered::{quarter_wave_band_edges, Layer, Periodic, Stack};
use phononic::materials::{Solid, ALUMINIUM, EPOXY, PMMA, STEEL, WATER};
use phononic::pwe::{self, SquareLattice};
use phononic::resonant::MassInMass;
use std::f64::consts::PI;

fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:.12e} vs {b:.12e}, tol {tol:.1e}");
}

/// Steel plate in water, normal incidence. T = 1 at thickness resonances
/// d = n lambda / 2, i.e. f = n c_p / (2 d).
#[test]
fn single_layer_thickness_resonances() {
    let d = 5.0e-3;
    let mut worst = 0.0f64;
    for n in 1..=5u32 {
        let f = fluid_coupling::thickness_resonance(&STEEL, d, n);
        let omega = 2.0 * PI * f;
        let layer = Layer::from_solid(&STEEL, d);
        let s = Stack::from_media(&WATER, vec![layer], &WATER).solve(omega);
        println!(
            "resonance n={n}: f={:.6} MHz, T={:.12}, R={:.3e}",
            f / 1e6,
            s.transmittance,
            s.reflectance
        );
        close(s.transmittance, 1.0, 1e-8);
        worst = worst.max((s.transmittance - 1.0).abs());
    }
    // Off resonance the plate reflects: well below the first resonance.
    let omega = 2.0 * PI * 0.25 * fluid_coupling::thickness_resonance(&STEEL, d, 1);
    let s = Stack::single(
        WATER.impedance(),
        Layer::from_solid(&STEEL, d),
        WATER.impedance(),
    )
    .solve(omega);
    println!("off resonance: T={:.6}", s.transmittance);
    assert!(s.transmittance < 0.5);
    println!("worst resonance error = {worst:.3e}");
}

/// Lossless energy conservation in the symmetric convention
/// `|r|^2 + (Z_out/Z_in) |t_v|^2 = 1`.
#[test]
fn energy_conservation_lossless_stacks() {
    let z0 = WATER.impedance();
    let zs = WATER.impedance();
    let layers = vec![
        Layer::lossless(STEEL.rho, STEEL.c_p, 0.8e-3),
        Layer::lossless(ALUMINIUM.rho, ALUMINIUM.c_p, 1.7e-3),
        Layer::lossless(PMMA.rho, PMMA.c_p, 0.3e-3),
        Layer::lossless(STEEL.rho, STEEL.c_p, 0.4e-3),
    ];
    let stack = Stack::new(z0, layers, zs);
    let mut worst = 0.0f64;
    for k in 1..=200 {
        let omega = 2.0 * PI * (0.1e6 + 9.9e6 * k as f64 / 200.0);
        let s = stack.solve(omega);
        let tv = s.velocity_transmission();
        let e = s.reflectance + (s.exit_z.re / s.incident_z.re) * tv.abs2();
        worst = worst.max((e - 1.0).abs());
        assert!(s.absorptance.abs() < 1e-9);
    }
    println!("worst |R + (Zout/Zin)|tv|^2 - 1| = {worst:.3e}");
    assert!(worst < 1e-9);
}

/// Lossy stack: `R + T + A = 1`, positive absorptance, and the lossless stack
/// recovered as loss -> 0.
#[test]
fn lossy_stack_absorbs() {
    let stack = Stack::new(
        WATER.impedance(),
        vec![Layer::with_loss(STEEL.rho, STEEL.c_p, 5e-3, 0.01)],
        WATER.impedance(),
    );
    let omega = 2.0 * PI * 1.0e6;
    let s = stack.solve(omega);
    println!(
        "lossy steel: R={:.6}, T={:.6}, A={:.6}",
        s.reflectance, s.transmittance, s.absorptance
    );
    close(s.energy() + s.absorptance, 1.0, 1e-12);
    assert!(s.absorptance > 0.0);
    let lossless = Stack::new(
        WATER.impedance(),
        vec![Layer::lossless(STEEL.rho, STEEL.c_p, 5e-3)],
        WATER.impedance(),
    )
    .solve(omega);
    assert!((s.transmittance - lossless.transmittance).abs() < 0.05);
}

/// Exact quarter-wave band edges, from the Rytov relation and from the
/// closed form `(4/pi) asin((Z1-Z2)/(Z1+Z2))`.
#[test]
fn quarter_wave_band_edges_analytic() {
    let (z1, z2, f0) = (2.0, 1.0, 1.0e6);
    let per = Periodic::new(
        Layer::lossless(z1, 1.0, 1.0 / (4.0 * f0)),
        Layer::lossless(z2, 1.0, 1.0 / (4.0 * f0)),
    );
    let (flo, fhi) = quarter_wave_band_edges(z1, z2, f0);
    println!("analytic edges: {:.6} .. {:.6} MHz", flo / 1e6, fhi / 1e6);
    // At the edges, cos(q Lambda) must be exactly -1.
    for f in [flo, fhi] {
        let h = per.rytov_cos(2.0 * PI * f);
        close(h.re, -1.0, 1e-10);
        assert!(h.im.abs() < 1e-12);
    }
    // Numerical scan of the exact Rytov gap matches the closed form.
    let steps = 2_000_000;
    let (lo, hi) = (0.5 * f0, 1.5 * f0);
    let (mut a, mut b) = (f64::NAN, f64::NAN);
    for i in 0..=steps {
        let f = lo + (hi - lo) * i as f64 / steps as f64;
        if per.in_band_gap(2.0 * PI * f) {
            if a.is_nan() {
                a = f;
            }
            b = f;
        }
    }
    println!("scanned edges:  {:.6} .. {:.6} MHz", a / 1e6, b / 1e6);
    assert!((a - flo).abs() / flo < 1e-3, "lower edge {a} vs {flo}");
    assert!((b - fhi).abs() / fhi < 1e-3, "upper edge {b} vs {fhi}");
}

/// The Bloch gap of the infinite quarter-wave stack coincides with the
/// transmission dip of a long finite stack.
#[test]
fn bloch_gap_matches_finite_stack_dip() {
    let (z1, z2, f0) = (46.8e6, 17.3e6, 1.0e6); // steel / aluminium scale
    let (c1, c2) = (5960.0, 6420.0);
    let per = Periodic::quarter_wave(z1 / c1, c1, z2 / c2, c2, f0);
    let (flo, fhi) = quarter_wave_band_edges(z1, z2, f0);
    let periods = 20usize;
    let mut layers = Vec::new();
    for _ in 0..periods {
        layers.push(per.a);
        layers.push(per.b);
    }
    let stack = Stack::new(z1, layers, z1); // half-spaces matched to layer a
    let steps = 200_000;
    let (lo, hi) = (0.4 * f0, 1.6 * f0);
    let (mut a, mut b) = (f64::NAN, f64::NAN);
    for i in 0..=steps {
        let f = lo + (hi - lo) * i as f64 / steps as f64;
        let t = stack.solve(2.0 * PI * f).transmittance;
        if t < 1e-6 {
            if a.is_nan() {
                a = f;
            }
            b = f;
        }
    }
    println!(
        "Bloch gap {:.6}..{:.6} MHz  finite dip {:.6}..{:.6} MHz",
        flo / 1e6,
        fhi / 1e6,
        a / 1e6,
        b / 1e6
    );
    let gc = 0.5 * (flo + fhi);
    let fc = 0.5 * (a + b);
    assert!((fc - gc).abs() / gc < 0.02, "centres differ");
    assert!((b - a).abs() / (fhi - flo) > 0.8, "finite dip too narrow");
}

/// Bloch dispersion converges to the Reuss effective speed as q Lambda -> 0.
#[test]
fn bloch_converges_to_reuss_speed() {
    // Two layers with a moderate contrast.
    let (rho1, c1) = (3000.0, 4000.0);
    let (rho2, c2) = (2000.0, 2000.0);
    let f0 = 1.0e6;
    let per = Periodic::quarter_wave(rho1, c1, rho2, c2, f0);
    let c_eff = laminate_effective_speed(
        &[rho1 * c1 * c1, rho2 * c2 * c2],
        &[rho1, rho2],
        &[per.a.thickness, per.b.thickness],
    );
    println!("c_eff = {c_eff:.6} m/s");
    let mut prev_err = f64::NAN;
    println!("  x=omega*L/c_eff     rel_err        err/x^2");
    for &x in &[0.4, 0.2, 0.1, 0.05, 0.025, 0.0125] {
        let omega = x * c_eff / per.period();
        let err = phase_velocity_error(&per, omega, c_eff);
        println!("  {x:>10.4}   {err:.6e}   {:.6e}", err / (x * x));
        if !prev_err.is_nan() {
            // Clean quadratic convergence: error falls 4x when x halves.
            assert!(
                (prev_err / err - 4.0).abs() < 1.5,
                "not quadratic: {prev_err} -> {err}"
            );
        }
        prev_err = err;
    }
    // The empirical envelope bounds the measured error at small x.
    for &x in &[0.05, 0.025, 0.0125] {
        let omega = x * c_eff / per.period();
        let err = phase_velocity_error(&per, omega, c_eff);
        let env = homogenisation_error(omega, per.period(), c_eff);
        println!(
            "  x={x:.4} err={err:.3e} envelope={env:.3e} ratio={:.3}",
            err / env
        );
        assert!(err <= env, "envelope too small at x={x}: {err} > {env}");
    }
}

/// Hashin–Shtrikman bounds bracket each other and sit inside Voigt/Reuss.
#[test]
fn hashin_shtrikman_ordering() {
    let (k1, g1) = (100.0e9, 80.0e9);
    let (k2, g2) = (10.0e9, 5.0e9);
    for f1 in [0.1, 0.25, 0.5, 0.75, 0.9] {
        let (k_lo, k_hi) = hashin_shtrikman_bulk(k1, g1, k2, g2, f1);
        let (g_lo, g_hi) = hashin_shtrikman_shear(k1, g1, k2, g2, f1);
        let k_reuss = reuss(&[k1, k2], &[f1, 1.0 - f1]);
        let k_voigt = voigt(&[k1, k2], &[f1, 1.0 - f1]);
        let g_reuss = reuss(&[g1, g2], &[f1, 1.0 - f1]);
        let g_voigt = voigt(&[g1, g2], &[f1, 1.0 - f1]);
        println!(
            "f1={f1:.2}: K {k_reuss:.3e} <= {k_lo:.3e} <= {k_hi:.3e} <= {k_voigt:.3e} | G {g_reuss:.3e} <= {g_lo:.3e} <= {g_hi:.3e} <= {g_voigt:.3e}"
        );
        assert!(k_reuss <= k_lo + 1e-6);
        assert!(k_lo <= k_hi + 1e-6);
        assert!(k_hi <= k_voigt + 1e-6);
        assert!(g_reuss <= g_lo + 1e-6);
        assert!(g_lo <= g_hi + 1e-6);
        assert!(g_hi <= g_voigt + 1e-6);
    }
    // End members are exact.
    let (k_lo, k_hi) = hashin_shtrikman_bulk(k1, g1, k2, g2, 1.0);
    close(k_lo, k1, 1e-6);
    close(k_hi, k1, 1e-6);
    let (k_lo, k_hi) = hashin_shtrikman_bulk(k1, g1, k2, g2, 0.0);
    close(k_lo, k2, 1e-6);
    close(k_hi, k2, 1e-6);
}

/// The effective-mass dispersion must reproduce the exact 2-DOF unit-cell
/// eigenfrequencies, and the negative-mass interval sits inside the gap.
#[test]
fn mass_in_mass_effective_mass_matches_eigen() {
    let m = MassInMass::new(1.0, 1.0, 1.0, 1.0, 1.0);
    let mut worst = 0.0f64;
    for k in 0..=200 {
        let q = PI / m.a * k as f64 / 200.0;
        let (w_lo, w_hi) = m.eigen_frequencies(q);
        for w in [w_lo, w_hi] {
            if w <= 0.0 {
                continue;
            }
            let arg = 1.0 - w * w * m.effective_mass(w) / (2.0 * m.k1);
            let err = (arg - (q * m.a).cos()).abs();
            worst = worst.max(err);
        }
    }
    println!("max |effective-mass cos - eigen cos| = {worst:.3e}");
    assert!(worst < 1e-12, "effective mass and eigenproblem disagree");

    let (w0, wtop) = m.negative_mass_interval();
    let (glo, ghi) = m.band_gap();
    println!("omega0={w0:.6}, negative mass up to {wtop:.6}; gap {glo:.6}..{ghi:.6}");
    assert!(m.effective_mass(0.5 * (w0 + wtop)) < 0.0);
    assert!(m.effective_mass(0.5 * w0) > 0.0);
    assert!(m.effective_mass(1.05 * wtop) > 0.0);
    // The negative-mass interval is inside the true lattice gap.
    assert!(glo <= w0 + 1e-12);
    assert!(ghi >= wtop - 1e-12);
    // And the true gap is where the Bloch wavevector is complex.
    let wm = 0.5 * (glo + ghi);
    assert!(m.wavevector(wm).im.abs() > 0.0);
    assert!(m.wavevector(0.5 * glo).im.abs() < 1e-12);
}

/// A soft inner resonator gives a gap whose wavelength is much larger than
/// the lattice constant: the sub-wavelength gap.
#[test]
fn mass_in_mass_subwavelength_gap() {
    let m = MassInMass::new(1.0, 1.0, 1.0, 1.0, 0.01);
    let ratio = m.gap_wavelength_over_lattice();
    let (glo, ghi) = m.band_gap();
    println!(
        "sub-wavelength: gap {:.6}..{:.6} rad/s, lambda/a={:.3}, Bragg/bottom={:.2}",
        glo,
        ghi,
        ratio,
        m.bragg_frequency() / glo
    );
    assert!(ratio > 20.0, "gap is not sub-wavelength: {ratio}");
    assert!(ghi < m.bragg_frequency());
    // Static limit: the effective mass tends to m1 + m2.
    close(m.effective_mass(1e-9), m.m1 + m.m2, 1e-6);
}

/// At normal incidence the fluid–solid–fluid solver must reduce exactly to the
/// longitudinal layered stack, for both energy and transmittance.
#[test]
fn fluid_plate_reduces_to_layered_at_normal_incidence() {
    let d = 5.0e-3;
    let mut worst = 0.0f64;
    for k in 0..=100 {
        let f = 0.1e6 + 9.9e6 * k as f64 / 100.0;
        let omega = 2.0 * PI * f;
        let a = fluid_coupling::solve(&WATER, &STEEL, d, omega, 0.0);
        let b = fluid_coupling::normal_incidence_reference(&WATER, &STEEL, d, omega);
        worst = worst.max((a.transmittance - b).abs());
        close(a.energy(), 1.0, 1e-8);
    }
    println!("max |oblique(theta=0) - layered| = {worst:.3e}");
    assert!(worst < 1e-8);
}

/// Energy conservation and the coincidence/Lamb transmission peak at oblique
/// incidence.
#[test]
fn fluid_plate_oblique_energy_and_coincidence() {
    let d = 5.0e-3;
    let f = 2.0e6;
    let omega = 2.0 * PI * f;
    let theta_p = fluid_coupling::coincidence_angle(&WATER, &STEEL).unwrap();
    let c_s0 = fluid_coupling::plate_speed_s0(&STEEL);
    let theta_s0 = fluid_coupling::coincidence_angle_for_speed(&WATER, c_s0).unwrap();
    println!(
        "P coincidence = {:.4} deg; S0 plate speed = {:.1} m/s, coincidence = {:.4} deg",
        theta_p.to_degrees(),
        c_s0,
        theta_s0.to_degrees()
    );

    let mut worst = 0.0f64;
    let steps = 8000;
    let samples: Vec<(f64, f64)> = (0..=steps)
        .map(|i| {
            let theta = (40.0_f64).to_radians() * i as f64 / steps as f64;
            let s = fluid_coupling::solve(&WATER, &STEEL, d, omega, theta);
            worst = worst.max((s.energy() - 1.0).abs());
            (theta, s.transmittance)
        })
        .collect();
    let mut peaks: Vec<(f64, f64)> = Vec::new();
    for i in 1..steps {
        let (_, t) = samples[i];
        if t > 0.99 && t >= samples[i - 1].1 && t > samples[i + 1].1 {
            peaks.push(samples[i]);
        }
    }
    print!("near-unity peaks (deg, T):");
    for &(a, t) in &peaks {
        print!(" ({:.2}, {t:.4})", a.to_degrees());
    }
    println!();
    println!("max |R+T-1| = {worst:.3e}");
    close(worst, 0.0, 1e-6);
    // Coincidence shows up as near-unity transmission peaks. The compressional
    // (trace = c_p) and extensional S0 (trace = c_plate) coincidences are both
    // present; fluid loading shifts each by a fraction of a degree.
    let near = |target: f64, tol_deg: f64| {
        peaks
            .iter()
            .any(|&(a, _)| (a - target).abs() < tol_deg.to_radians())
    };
    assert!(near(theta_p, 0.6), "no P coincidence peak");
    assert!(near(theta_s0, 0.8), "no S0 coincidence peak");

    // A coincidence angle exists whenever the fluid is slower than the plate
    // speed; water is slower than all three solids here.
    let alu = fluid_coupling::coincidence_angle(&WATER, &ALUMINIUM).unwrap();
    let pmma = fluid_coupling::coincidence_angle(&WATER, &PMMA).unwrap();
    println!(
        "P coincidence: Al {:.3} deg, PMMA {:.3} deg",
        alu.to_degrees(),
        pmma.to_degrees()
    );
    // A faster plate has an earlier coincidence angle.
    assert!(alu < theta_p);
    assert!(pmma > theta_p);
}

/// The fluid–plate solver conserves energy for a lossless plate at oblique
/// incidence across a range of angles, including beyond the shear critical
/// angle where the SV wave is evanescent. The thickness is kept small enough
/// that the evanescent decay exponent stays moderate; the transfer matrix is
/// ill-conditioned when `|Im k_zS| d` becomes very large (documented limit).
#[test]
fn fluid_plate_energy_beyond_critical() {
    let d = 5.0e-3;
    let f = 1.0e6;
    let omega = 2.0 * PI * f;
    let theta_s = (WATER.c / STEEL.c_s).asin();
    println!("shear critical angle = {:.4} deg", theta_s.to_degrees());
    let mut worst = 0.0f64;
    for i in 0..=800 {
        let theta = (0.2_f64 + 39.6 * i as f64 / 800.0).to_radians();
        let s = fluid_coupling::solve(&WATER, &STEEL, d, omega, theta);
        worst = worst.max((s.energy() - 1.0).abs());
    }
    println!("max |R+T-1| beyond critical = {worst:.3e}");
    assert!(theta_s.to_degrees() < 39.0);
    close(worst, 0.0, 1e-8);
}

/// Laminate effective modulus/density are the exact Reuss/volume averages.
#[test]
fn laminate_reuss_exact() {
    let layers = [
        (STEEL.rho * STEEL.c_p * STEEL.c_p, 1.0e-3),
        (PMMA.rho * PMMA.c_p * PMMA.c_p, 2.0e-3),
    ];
    let density_layers = [(STEEL.rho, 1.0e-3), (PMMA.rho, 2.0e-3)];
    let f1 = 1.0 / 3.0;
    let m_eff = laminate_effective_modulus(&layers);
    let want = 1.0 / (f1 / layers[0].0 + (1.0 - f1) / layers[1].0);
    close(m_eff / want, 1.0, 1e-12);
    let rho_eff = laminate_effective_density(&density_layers);
    let want_rho = STEEL.rho * f1 + PMMA.rho * (1.0 - f1);
    close(rho_eff / want_rho, 1.0, 1e-12);
    println!(
        "M_eff={m_eff:.6e}, rho_eff={rho_eff:.6}, c_eff={:.6}",
        (m_eff / rho_eff).sqrt()
    );
}

// ---------------------------------------------------------------------------
// Global (stiffness) matrix fluid-plate solver
// ---------------------------------------------------------------------------

/// The unconditionally stable global-matrix solver must conserve energy to
/// round-off in the regime where the historical 4x4 transfer matrix is
/// ill-conditioned: past the shear critical angle at large `f d`. The old
/// propagator mixes growing and decaying evanescent partial waves and its
/// energy error grows without bound, reaching O(1) at 3 MHz / 5 mm steel; the
/// stable solver stays at ~1e-14. Both are printed.
#[test]
fn fluid_plate_stable_energy_high_fd() {
    let d = 5.0e-3;
    for &f in &[1.0e6, 3.0e6, 5.0e6, 10.0e6] {
        let omega = 2.0 * PI * f;
        let mut worst_stable = 0.0f64;
        let mut worst_transfer = 0.0f64;
        let steps = 800;
        for i in 0..=steps {
            let theta = (0.2 + 39.6 * i as f64 / steps as f64).to_radians();
            let s = fluid_coupling::solve_stable(&WATER, &STEEL, d, omega, theta);
            worst_stable = worst_stable.max((s.energy() - 1.0).abs());
            let o = fluid_coupling::solve_transfer_matrix(&WATER, &STEEL, d, omega, theta);
            worst_transfer = worst_transfer.max((o.energy() - 1.0).abs());
        }
        println!(
            "f={:.1} MHz, d={} mm: stable |R+T-1| <= {:.3e}, transfer |R+T-1| <= {:.3e}",
            f / 1e6,
            d * 1e3,
            worst_stable,
            worst_transfer
        );
        assert!(
            worst_stable < 1e-9,
            "stable solver energy error {worst_stable:.3e}"
        );
    }
    // The documented failure: the old method is genuinely bad at 3 MHz / 5 mm.
    let omega = 2.0 * PI * 3.0e6;
    let mut old = 0.0f64;
    for i in 0..=800 {
        let theta = (0.2 + 39.6 * i as f64 / 800.0).to_radians();
        let o = fluid_coupling::solve_transfer_matrix(&WATER, &STEEL, d, omega, theta);
        old = old.max((o.energy() - 1.0).abs());
    }
    assert!(
        old > 1.0,
        "expected the transfer matrix to fail at 3 MHz / 5 mm, got {old:.3e}"
    );
}

/// Where the transfer matrix is well conditioned the stable solver must
/// reproduce it, so the replacement does not change the physics there.
#[test]
fn fluid_plate_stable_matches_transfer_matrix() {
    let d = 5.0e-3;
    let omega = 2.0 * PI * 0.5e6; // small f d: evanescent decay is modest
    let mut worst = 0.0f64;
    for i in 0..=200 {
        let theta = (40.0 * i as f64 / 200.0).to_radians();
        let a = fluid_coupling::solve_stable(&WATER, &STEEL, d, omega, theta);
        let b = fluid_coupling::solve_transfer_matrix(&WATER, &STEEL, d, omega, theta);
        worst = worst.max((a.r.re - b.r.re).abs().max((a.t.re - b.t.re).abs()));
    }
    println!("max |stable - transfer| coefficient at 0.5 MHz = {worst:.3e}");
    assert!(worst < 1e-9);
}

// ---------------------------------------------------------------------------
// Free-plate Rayleigh-Lamb dispersion
// ---------------------------------------------------------------------------

/// The pole-free Rayleigh-Lamb residual and the (independently scaled) global
/// matrix determinant vanish on the same speeds, and the cross-multiplied
/// spurious roots at the S and P body-wave speeds are excluded. Thin-plate
/// limits: `S0 -> sqrt(E/(rho(1-nu^2)))` and `A0 -> 0`.
#[test]
fn rayleigh_lamb_dispersion_free_plate() {
    let d = 5.0e-3;
    // S0 approaches the thin-plate extensional speed and A0 approaches zero
    // as f d -> 0; both are checked at two small f d.
    let s0_ref = fluid_coupling::plate_speed_s0(&STEEL);
    let mut s0_prev = f64::NAN;
    let mut a0_prev = f64::NAN;
    for &f in &[0.05e6, 0.02e6] {
        let omega = 2.0 * PI * f;
        let modes = lamb::modes(&STEEL, d, omega, 1.0, 1.2 * STEEL.c_p, 400_000);
        let s0 = modes
            .iter()
            .filter(|m| m.kind == lamb::LambKind::Symmetric)
            .map(|m| m.speed)
            .fold(f64::INFINITY, f64::min);
        let a0 = modes
            .iter()
            .filter(|m| m.kind == lamb::LambKind::Antisymmetric)
            .map(|m| m.speed)
            .fold(f64::INFINITY, f64::min);
        println!(
            "fd={:.3} MHz.mm: S0={s0:.1} m/s (ref {s0_ref:.1}), A0={a0:.1} m/s",
            (f / 1e6) * (d * 1e3)
        );
        assert!((s0 - s0_ref).abs() / s0_ref < 0.01, "S0 {s0} vs {s0_ref}");
        if s0_prev.is_finite() {
            assert!(
                s0 > s0_prev,
                "S0 should approach the plate speed as fd -> 0"
            );
            assert!(a0 < a0_prev, "A0 should vanish as fd -> 0");
        }
        s0_prev = s0;
        a0_prev = a0;
    }

    // Cross-check the residual and the 4x4 global determinant at 3 MHz, and
    // confirm that the body-wave speeds are not reported as modes.
    let f = 3.0e6;
    let omega = 2.0 * PI * f;
    let modes = lamb::modes(&STEEL, d, omega, WATER.c + 1.0, 1.2 * STEEL.c_p, 400_000);
    let mut worst_det_ratio = 0.0f64;
    let mut n_checked = 0;
    for m in &modes {
        // Skip speeds at which the split residual is only the removable root.
        if (m.speed - STEEL.c_s).abs() < 1e-4 * STEEL.c_s
            || (m.speed - STEEL.c_p).abs() < 1e-4 * STEEL.c_p
        {
            continue;
        }
        for factor in [0.995, 1.005] {
            let c = m.speed * factor;
            let det = fluid_coupling::free_plate_determinant(&STEEL, d, omega, c).abs();
            let at_mode = fluid_coupling::free_plate_determinant(&STEEL, d, omega, m.speed).abs();
            worst_det_ratio = worst_det_ratio.max(at_mode / det);
        }
        n_checked += 1;
    }
    println!(
        "checked {n_checked} Lamb modes: worst |det(mode)|/|det(+-0.5%)| = {worst_det_ratio:.3e}"
    );
    assert!(n_checked > 5);
    assert!(
        worst_det_ratio < 1e-4,
        "determinant is not minimal on a Lamb mode"
    );
    // The spurious roots are filtered: no reported speed sits on c_s or c_p.
    for m in &modes {
        assert!(
            (m.speed - STEEL.c_s).abs() > 1e-6 * STEEL.c_s,
            "spurious c_s root reported"
        );
        assert!(
            (m.speed - STEEL.c_p).abs() > 1e-6 * STEEL.c_p,
            "spurious c_p root reported"
        );
    }
}

/// Transmission peaks of the water-loaded steel plate sit at the coincidence
/// angles of the free plate's leaky Lamb modes: for every mode faster than the
/// shear speed (coincidence below the shear critical angle) there is a
/// near-unity transmission maximum within 1.5 degrees of `asin(c_f/c)`. The
/// fluid loading shifts each peak, so the window is not the order of a scan
/// step.
#[test]
fn immersed_plate_peaks_at_lamb_coincidence() {
    let d = 5.0e-3;
    for &f in &[2.0e6, 3.0e6] {
        let omega = 2.0 * PI * f;
        let modes = lamb::modes(&STEEL, d, omega, WATER.c + 1.0, 1.2 * STEEL.c_p, 400_000);
        let mut checked = 0usize;
        for m in &modes {
            if m.speed <= STEEL.c_s {
                continue; // coincidence lies beyond the shear critical angle
            }
            let Some(theta_c) = lamb::coincidence_angle(WATER.c, m.speed) else {
                continue;
            };
            if theta_c > 40.0_f64.to_radians() {
                continue;
            }
            let n = 4000;
            let (lo, hi) = (
                theta_c - 1.5_f64.to_radians(),
                theta_c + 1.5_f64.to_radians(),
            );
            let mut best = (f64::NAN, 0.0f64);
            let mut prev2 = f64::NAN;
            let mut prev = f64::NAN;
            for i in 0..=n {
                let a = lo + (hi - lo) * i as f64 / n as f64;
                let t = fluid_coupling::solve(&WATER, &STEEL, d, omega, a).transmittance;
                if prev2.is_finite() && prev.is_finite() && prev >= prev2 && prev > t {
                    let a_peak = a - (hi - lo) / n as f64;
                    if (a_peak - theta_c).abs() <= 1.5_f64.to_radians() && prev > best.1 {
                        best = (a_peak, prev);
                    }
                }
                prev2 = prev;
                prev = t;
            }
            println!(
                "f={:.0} MHz: Lamb c={:.1} m/s ({:?}) -> coincidence {:.3} deg, peak {:.3} deg T={:.3}",
                f / 1e6,
                m.speed,
                m.kind,
                theta_c.to_degrees(),
                best.0.to_degrees(),
                best.1
            );
            assert!(
                best.1 > 0.9,
                "no near-unity peak near Lamb coincidence for c={}",
                m.speed
            );
            checked += 1;
        }
        assert!(checked >= 5, "too few Lamb coincidence peaks checked");
    }
}

// ---------------------------------------------------------------------------
// 2D plane-wave expansion
// ---------------------------------------------------------------------------

/// Build a `Solid` with a prescribed shear modulus and density (the P speed is
/// irrelevant to the scalar and in-plane operators in the tests that use this).
fn solid_from_mu(name: &'static str, mu: f64, rho: f64) -> Solid {
    Solid {
        name,
        rho,
        c_p: 2.0 * (mu / rho).sqrt(),
        c_s: (mu / rho).sqrt(),
        source: "test fixture",
    }
}

/// The anti-plane (SH) scalar problem `div(mu grad w) + rho omega^2 w = 0` is
/// mathematically identical to TM photonics
/// `div((1/eps) grad H) + (omega^2/c^2) H = 0` under `mu <-> 1/eps`,
/// `rho <-> 1/c^2`. Setting `mu = 1/eps` and `rho = 1` must make the SH band
/// frequencies equal the TM frequencies evaluated with `c = 1`.
#[test]
fn sh_reduces_to_tm_photonics() {
    let a = 1.0;
    let f = 0.35;
    let r = (f * a * a / PI).sqrt();
    let (eps_inc, eps_mat) = (12.0, 1.0);
    // mu = 1/eps, rho = 1 (so c_s^2 = mu).
    let inc = solid_from_mu("inc", 1.0 / eps_inc, 1.0);
    let mat = solid_from_mu("mat", 1.0 / eps_mat, 1.0);
    let lat = SquareLattice::new(a, r, inc, mat);
    let m = 3;
    let mut worst = 0.0f64;
    for &(kx, ky) in &[
        (0.0, 0.0),
        (0.3, 0.0),
        (PI / a, 0.0),
        (PI / a, PI / (2.0 * a)),
    ] {
        let sh = pwe::sh_frequencies(&lat, m, kx, ky);
        let tm = pwe::tm_frequencies(eps_inc, eps_mat, a, r, m, kx, ky);
        for i in 0..8.min(sh.len()) {
            worst = worst.max((sh[i] - tm[i]).abs());
        }
    }
    println!("max |SH - TM| frequency (c=1, mu=1/eps) = {worst:.3e}");
    assert!(worst < 1e-9);
}

/// An empty (homogeneous) lattice is the folded free-phonon spectrum
/// `omega = c |k + G|`, which is an exact check of the scalar PWE, the disk
/// Fourier coefficients and the eigensolver.
#[test]
fn pwe_homogeneous_limit() {
    let a = 1.0;
    let f = 0.3;
    let r = (f * a * a / PI).sqrt();
    let homo = SquareLattice::new(a, r, EPOXY, EPOXY);
    let c = EPOXY.c_s;
    let m = 3;
    let mut worst = 0.0f64;
    for &(kx, ky) in &[(0.0, 0.0), (1.0, 0.3), (PI / a, 0.0)] {
        let w = pwe::sh_frequencies(&homo, m, kx, ky);
        // Compare every computed eigenvalue with the nearest exact folded value.
        let mut exact: Vec<f64> = Vec::new();
        for nx in -m..=m {
            for ny in -m..=m {
                let gx = kx + 2.0 * PI * nx as f64 / a;
                let gy = ky + 2.0 * PI * ny as f64 / a;
                exact.push(c * (gx * gx + gy * gy).sqrt());
            }
        }
        exact.sort_by(|x, y| x.partial_cmp(y).unwrap());
        for i in 0..w.len() {
            worst = worst.max((w[i] - exact[i]).abs());
        }
    }
    println!("max |PWE homogeneous - exact free phonon| = {worst:.3e} rad/s");
    assert!(worst < 1e-6 * c / a);
}

/// The in-plane P-SV operator reduces exactly to the folded isotropic Navier
/// dispersion in a homogeneous medium: for every plane wave `G` the 2x2 block
/// has eigenvalues `c_s^2 |k+G|^2` and `c_p^2 |k+G|^2`. This catches any
/// lambda/mu ordering error in the mixed operator.
#[test]
fn pwe_inplane_homogeneous_limit() {
    let a = 1.0;
    let f = 0.3;
    let r = (f * a * a / PI).sqrt();
    let homo = SquareLattice::new(a, r, EPOXY, EPOXY);
    let m = 3;
    let (cp, cs) = (EPOXY.c_p, EPOXY.c_s);
    let mut worst = 0.0f64;
    for &(kx, ky) in &[(0.0, 0.0), (1.0, 0.3), (PI / a, 0.0)] {
        let w = pwe::inplane_frequencies(&homo, m, kx, ky);
        let mut exact: Vec<f64> = Vec::new();
        for nx in -m..=m {
            for ny in -m..=m {
                let gx = kx + 2.0 * PI * nx as f64 / a;
                let gy = ky + 2.0 * PI * ny as f64 / a;
                let mag = (gx * gx + gy * gy).sqrt();
                exact.push(cp * mag);
                exact.push(cs * mag);
            }
        }
        exact.sort_by(|x, y| x.partial_cmp(y).unwrap());
        for i in 0..w.len() {
            worst = worst.max((w[i] - exact[i]).abs());
        }
    }
    println!("max |in-plane PWE homogeneous - exact Navier| = {worst:.3e} rad/s");
    assert!(worst < 1e-6 * cp / a);
}

/// Anti-plane effective medium. The long-wavelength slope of the lowest SH
/// band must lie between the rigorous Reuss (harmonic) and Voigt (arithmetic)
/// modulus bounds from `homogenize`, and approach the published 2D
/// Maxwell-Garnett / Hashin effective shear modulus as the plane-wave count
/// grows. The 3D Hashin-Shtrikman bounds are also reported.
#[test]
fn sh_low_frequency_effective_medium() {
    let a = 1.0;
    let f = 0.4;
    let r = (f * a * a / PI).sqrt();
    let lat = SquareLattice::new(a, r, STEEL, EPOXY);
    let rho_eff = f * STEEL.rho + (1.0 - f) * EPOXY.rho;
    let (mu_i, mu_m) = (STEEL.mu(), EPOXY.mu());
    let mu_voigt = voigt(&[mu_i, mu_m], &[f, 1.0 - f]);
    let mu_reuss = reuss(&[mu_i, mu_m], &[f, 1.0 - f]);
    // Hashin (1965) exact dilute-limit / 2D Hashin-Shtrikman lower bound for
    // anti-plane shear of circular cylinders in a matrix.
    let mu_mg = mu_m + f / (1.0 / (mu_i - mu_m) + (1.0 - f) / (2.0 * mu_m));
    let (g_lo3, g_hi3) =
        hashin_shtrikman_shear(STEEL.bulk(), STEEL.mu(), EPOXY.bulk(), EPOXY.mu(), f);
    println!(
        "rho_eff={rho_eff:.1}; mu Reuss={mu_reuss:.4e}, Voigt={mu_voigt:.4e}, 2D-MG={mu_mg:.4e}"
    );
    println!("3D HS shear G in [{g_lo3:.4e}, {g_hi3:.4e}] (not the 2D anti-plane bracket)");
    let k = 1.0e-3; // k a << 1
    let mut prev_err = f64::NAN;
    for m in [5, 7, 9, 11] {
        let w = pwe::sh_frequencies(&lat, m, k, 0.0);
        let c = w[0] / k;
        let mu_eff = c * c * rho_eff;
        let c_lo = (mu_reuss / rho_eff).sqrt();
        let c_hi = (mu_voigt / rho_eff).sqrt();
        let err = (c - (mu_mg / rho_eff).sqrt()).abs() / (mu_mg / rho_eff).sqrt();
        println!(
            "m={m:2}: c={c:.4} m/s, mu_eff={mu_eff:.4e} in [{c_lo:.4},{c_hi:.4}], |c-c_MG|/c_MG={err:.3e}"
        );
        assert!(
            c > c_lo && c < c_hi,
            "SH speed outside the Reuss-Voigt bracket"
        );
        if prev_err.is_finite() {
            assert!(err < prev_err, "effective-medium error should fall with m");
        }
        prev_err = err;
    }
}

/// Convergence of the first SH gap of a square lattice of steel cylinders in
/// epoxy at filling fraction 0.4, alongside the published result. See
/// `analysis/phononic.md`: this is the high-contrast case for which the
/// plane-wave expansion is known to converge slowly (the displacement is
/// discontinuous), so the gap edges are extrapolated as `m` is refined and
/// every measured edge is printed. The test asserts a complete gap exists and
/// is broad, which is the qualitative published benchmark (Vasseur et al.,
/// J. Phys.: Condens. Matter 6, 8759 (1994), square steel/epoxy).
#[test]
fn sh_band_gap_steel_epoxy() {
    let a = 1.0;
    let f = 0.4;
    let r = (f * a * a / PI).sqrt();
    let lat = SquareLattice::new(a, r, STEEL, EPOXY);
    let scale = EPOXY.c_s / a; // Omega = omega a / c_T(epoxy)
                               // Complete-gap extrema for the lowest gap sit at the zone boundary.
    let x = PI / a;
    let pts = [(0.0, 0.0), (x, 0.0), (x, x)];
    println!("first SH gap, steel/epoxy f=0.4 (Omega = omega a / c_T):");
    let mut edges = Vec::new();
    for m in [3, 5, 7, 9, 11] {
        let (mut lo, mut hi) = (0.0f64, f64::INFINITY);
        for (kx, ky) in pts {
            let w = pwe::sh_frequencies(&lat, m, kx, ky);
            lo = lo.max(w[0]);
            hi = hi.min(w[1]);
        }
        let (lo, hi) = (lo / scale, hi / scale);
        println!(
            "  m={m:2}: {:.4} .. {:.4}  (reduced freq f a/c_T {:.4} .. {:.4})",
            lo,
            hi,
            lo / (2.0 * PI),
            hi / (2.0 * PI)
        );
        edges.push((lo, hi, m as f64));
        assert!(hi > lo, "no complete SH gap at m={m}");
    }
    // Extrapolate the monotone edges with a 1/m fit using the two finest runs.
    let fit = |which: usize| -> f64 {
        let (a1, b1, m1) = (
            edges[edges.len() - 1].0,
            edges[edges.len() - 1].1,
            edges[edges.len() - 1].2,
        );
        let (a0, b0, m0) = (
            edges[edges.len() - 2].0,
            edges[edges.len() - 2].1,
            edges[edges.len() - 2].2,
        );
        let (v1, v0) = if which == 0 { (a1, a0) } else { (b1, b0) };
        let slope = (v1 - v0) / (1.0 / m1 - 1.0 / m0);
        v1 - slope / m1
    };
    let (lo_inf, hi_inf) = (fit(0), fit(1));
    let width = hi_inf - lo_inf;
    let mid = 0.5 * (hi_inf + lo_inf);
    println!(
        "  extrapolated (1/m): {lo_inf:.4} .. {hi_inf:.4}, width/mid = {:.3}",
        width / mid
    );
    // The gap is broad, as published for this high-contrast system.
    assert!(width / mid > 0.3, "gap is not broad");
    // Refinement changes the edges monotonically and by a shrinking amount.
    for w in edges.windows(2) {
        assert!(w[1].0 < w[0].0, "lower edge not converging");
        assert!(w[1].1 < w[0].1, "upper edge not converging");
    }

    // Confirm the zone-boundary extrema are the complete-gap extrema by
    // scanning a coarse Gamma-X-M-Gamma path at m=5.
    let m_path = 5;
    let steps = 10;
    let (mut lo_path, mut hi_path) = (0.0f64, f64::INFINITY);
    for i in 0..=steps {
        let t = i as f64 / steps as f64;
        let ks = [(x * t, 0.0), (x, x * t), (x * (1.0 - t), x * (1.0 - t))];
        for (kx, ky) in ks {
            let w = pwe::sh_frequencies(&lat, m_path, kx, ky);
            lo_path = lo_path.max(w[0]);
            hi_path = hi_path.min(w[1]);
        }
    }
    let (lo_path, hi_path) = (lo_path / scale, hi_path / scale);
    let (lo_x, hi_x) = (edges[1].0, edges[1].1);
    println!(
        "  full-path check at m={m_path}: {lo_path:.4} .. {hi_path:.4} vs Gamma/X/M {lo_x:.4} .. {hi_x:.4}"
    );
    assert!((lo_path - lo_x).abs() < 1e-6 && (hi_path - hi_x).abs() < 1e-6);
}

/// The in-plane mixed P-SV operator: a square steel/epoxy crystal at f=0.4 has
/// a gap between bands 2 and 3 over the high-symmetry points Gamma, X, M. The
/// same slow convergence applies; the test checks the gap exists for two
/// plane-wave counts and prints the edges. A full-BZ scan is not run (the
/// in-plane solver is the expensive one), so this is not a completeness proof.
#[test]
fn inplane_band_gap_steel_epoxy() {
    let a = 1.0;
    let f = 0.4;
    let r = (f * a * a / PI).sqrt();
    let lat = SquareLattice::new(a, r, STEEL, EPOXY);
    let scale = EPOXY.c_s / a;
    let x = PI / a;
    let pts = [(0.0, 0.0), (x, 0.0), (x, x)]; // Gamma, X, M
    for m in [3, 5] {
        let mut lo = 0.0f64;
        let mut hi = f64::INFINITY;
        for (kx, ky) in pts {
            let w = pwe::inplane_frequencies(&lat, m, kx, ky);
            lo = lo.max(w[2]); // top of band 2
            hi = hi.min(w[3]); // bottom of band 3
        }
        println!(
            "in-plane m={m}: complete gap {:.4} .. {:.4} (Omega = omega a / c_T)",
            lo / scale,
            hi / scale
        );
        assert!(hi > lo, "no complete in-plane gap at m={m}");
    }
}
