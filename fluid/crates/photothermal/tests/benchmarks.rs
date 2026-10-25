//! Benchmarks for the photothermal crate. Every physical result is compared
//! with an exact or published value; measured errors are printed so
//! `analysis/photothermal.md` can quote them. No tolerance is loosened to make
//! a case pass.

use photothermal::fft::{fft_real, ifft_real};
use photothermal::indicators::*;
use photothermal::photoacoustic::*;
use photothermal::thermal::*;
use std::f64::consts::PI;

fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:.12e} vs {b:.12e}, tol {tol:.1e}");
}

fn rel(a: f64, b: f64) -> f64 {
    (a - b).abs() / b.abs().max(1e-300)
}

// ---------------------------------------------------------------------------
// Numerical helper
// ---------------------------------------------------------------------------

#[test]
fn erfc_reference_values() {
    // Known values of erfc at half-integer arguments.
    let cases = [
        (0.0, 1.0),
        (0.5, 0.479_500_122_186_953_5),
        (1.0, 0.157_299_207_050_285_13),
        (2.0, 0.004_677_734_981_047_266),
        (3.0, 2.209_049_699_858_544e-5),
    ];
    let mut worst = 0.0f64;
    for (x, want) in cases {
        let got = erfc(x);
        let e = (got - want).abs();
        println!("erfc({x}) = {got:.16e} (want {want:.16e}, err {e:.2e})");
        worst = worst.max(e);
    }
    assert!(worst < 5e-13, "erfc worst error {worst:e}");
}

// ---------------------------------------------------------------------------
// 1. Thermal
// ---------------------------------------------------------------------------

/// Semi-infinite absorber with a constant surface flux: compare the cell-centre
/// temperature with the exact erfc solution.
#[test]
fn semi_infinite_surface_flux_matches_erfc() {
    let q = 1.0e4; // W/m^2 absorbed
    let cells = 400;
    let thickness = 0.02; // 2 cm, far beyond the diffusive length at these times
    let model = ThermalModel::new(
        vec![ThermalLayer {
            material: WATER,
            thickness,
        }],
        cells,
        Boundary::Flux(q),
        Boundary::Adiabatic,
    );
    let n = model.n();
    let dt = 5.0e-4;
    let steps = 2000; // t = 1 s
    let mut t = vec![0.0; n];
    model.run(&mut t, &vec![0.0; n], dt, steps);
    let tf = dt * steps as f64;

    let mut worst = 0.0f64;
    let mut sum2 = 0.0;
    let mut ref2 = 0.0;
    for i in 0..n {
        let z = model.z[i];
        // Only compare where the finite back boundary and the smeared source
        // cell do not dominate: z greater than one cell, inside the diffusive
        // length.
        if z < model.dz[i] {
            continue;
        }
        let want = semi_infinite_surface_flux(&WATER, z, tf, q);
        if want.abs() < 1e-6 {
            continue;
        }
        worst = worst.max((t[i] - want).abs());
        sum2 += (t[i] - want).powi(2);
        ref2 += want * want;
    }
    let l2 = (sum2 / ref2).sqrt();
    println!("surface flux erfc: worst abs err {worst:.3e} K, relative L2 {l2:.3e} at t={tf} s");
    assert!(l2 < 2.0e-2, "relative L2 error {l2:e} too large");

    // Surface temperature rise at z ~ dz/2 against the closed form.
    let zc = model.z[0];
    let want_c = semi_infinite_surface_flux(&WATER, zc, tf, q);
    println!(
        "  cell-centre T({zc:.3e} m) = {:.6} K, exact {want_c:.6} K, err {:.2e}",
        t[0],
        (t[0] - want_c).abs()
    );
}

/// Steady Beer–Lambert heating in a slab with fixed-temperature ends, and
/// second-order spatial convergence.
#[test]
fn beer_lambert_steady_state_second_order() {
    let mu_a = 1200.0; // 1/m
    let i0 = 1.0e4; // W/m^2
    let thickness = 1.0e-3;
    let del = thickness;
    let k = WATER.k;
    // Exact steady state with T(0)=T(L)=0.
    let a = -i0 / (k * mu_a);
    let c1 = a * (1.0 - (-mu_a * thickness).exp()) / thickness;
    let exact = |z: f64| a * ((-mu_a * z).exp() - 1.0) + c1 * z;

    let mut prev: Option<(usize, f64)> = None;
    for cells in [50usize, 100, 200, 400] {
        let model = ThermalModel::new(
            vec![ThermalLayer {
                material: WATER,
                thickness,
            }],
            cells,
            Boundary::Fixed(0.0),
            Boundary::Fixed(0.0),
        );
        let source: Vec<f64> = (0..model.n())
            .map(|i| {
                let z0 = model.z[i] - 0.5 * model.dz[i];
                beer_lambert_cell_average(mu_a, i0, z0, model.dz[i])
            })
            .collect();
        let t = model.steady_state(&source);
        let mut worst = 0.0f64;
        for i in 0..model.n() {
            let e = (t[i] - exact(model.z[i])).abs();
            worst = worst.max(e);
        }
        if let Some((pc, pe)) = prev {
            let order = (pe / worst).ln() / ((cells as f64 / pc as f64).ln());
            println!(
                "steady Beer-Lambert: cells={cells}, max err {worst:.3e} K, observed order {order:.3}"
            );
            assert!(order > 1.7, "spatial order {order} too low");
        } else {
            println!("steady Beer-Lambert: cells={cells}, max err {worst:.3e} K");
        }
        prev = Some((cells, worst));
    }
    let _ = del;
}

/// Transient Beer–Lambert heating of a slab with Dirichlet ends against the
/// exact Fourier series, and first-order time convergence.
#[test]
fn beer_lambert_transient_series_and_time_order() {
    let mu_a = 1000.0;
    let i0 = 1.0e4;
    let thickness = 1.0e-3;
    let k = WATER.k;
    let alpha = WATER.alpha();
    let a = -i0 / (k * mu_a);
    let c1 = a * (1.0 - (-mu_a * thickness).exp()) / thickness;
    let tss = |z: f64| a * ((-mu_a * z).exp() - 1.0) + c1 * z;

    // Series coefficients B_n = -(2/L) int_0^L Tss sin(q_n z) dz.
    let modes = 300;
    let series = |z: f64, t: f64| -> f64 {
        let mut v = tss(z);
        for n in 1..=modes {
            let qn = n as f64 * PI / thickness;
            let sign = if n % 2 == 0 { 1.0 } else { -1.0 };
            let i1 = (1.0 - sign) / qn;
            let i2 = -thickness * sign / qn;
            let i3 = qn * (1.0 - sign * (-mu_a * thickness).exp()) / (mu_a * mu_a + qn * qn);
            let bn = -(2.0 / thickness) * (a * (i3 - i1) + c1 * i2);
            v += bn * (qn * z).sin() * (-alpha * qn * qn * t).exp();
        }
        v
    };

    let cells = 800;
    let model = ThermalModel::new(
        vec![ThermalLayer {
            material: WATER,
            thickness,
        }],
        cells,
        Boundary::Fixed(0.0),
        Boundary::Fixed(0.0),
    );
    let source: Vec<f64> = (0..model.n())
        .map(|i| {
            let z0 = model.z[i] - 0.5 * model.dz[i];
            beer_lambert_cell_average(mu_a, i0, z0, model.dz[i])
        })
        .collect();
    let tf = 5.0;

    let mut prev: Option<(f64, f64)> = None;
    for dt in [0.02f64, 0.01, 0.005] {
        let mut t = vec![0.0; model.n()];
        let steps = (tf / dt).round() as usize;
        model.run(&mut t, &source, dt, steps);
        // Compare near mid-depth, away from the boundary layers.
        let mut sum2 = 0.0;
        let mut ref2 = 0.0;
        for i in 0..model.n() {
            let want = series(model.z[i], tf);
            sum2 += (t[i] - want).powi(2);
            ref2 += want * want;
        }
        let l2 = (sum2 / ref2).sqrt();
        println!("transient series: dt={dt}, relative L2 {l2:.3e}");
        if let Some((pdt, perr)) = prev {
            let order = (perr / l2).ln() / (pdt / dt).ln();
            println!("  observed time order {order:.3}");
            assert!(order > 0.7, "time order {order} too low");
        }
        prev = Some((dt, l2));
    }
}

/// Energy conservation with insulating boundaries: every deposited joule is
/// stored, to round-off.
#[test]
fn energy_conservation_layered_adiabatic() {
    let mu_a = 500.0;
    let i0 = 2.0e4;
    let layers = vec![
        ThermalLayer {
            material: GOLD,
            thickness: 100e-9,
        },
        ThermalLayer {
            material: PYREX,
            thickness: 1e-3,
        },
        ThermalLayer {
            material: WATER,
            thickness: 2e-3,
        },
    ];
    let model = ThermalModel::new(layers, 40, Boundary::Adiabatic, Boundary::Adiabatic);
    let z_off = model.layers[0].thickness; // absorption starts in the glass
    let source: Vec<f64> = (0..model.n())
        .map(|i| {
            let z0 = (model.z[i] - 0.5 * model.dz[i] - z_off).max(0.0);
            beer_lambert_cell_average(mu_a, i0, z0, model.dz[i])
        })
        .collect();
    let dt = 1e-5;
    let steps = 5000;
    let tf = dt * steps as f64;
    let mut t = vec![0.0; model.n()];
    model.run(&mut t, &source, dt, steps);

    // Storage is measured relative to 0 K, so subtract the initial stored value
    // (which is zero here).
    let stored = model.stored_energy(&t);
    let deposited = model.source_energy(&source, tf);
    println!(
        "energy: deposited {deposited:.6e} J/m^2, stored {stored:.6e} J/m^2, rel err {:.2e}",
        rel(stored, deposited)
    );
    close(stored / deposited, 1.0, 1e-7);
}

/// At steady state with no internal source and an insulating back, the heat
/// flux is the same through every face, including material interfaces.
#[test]
fn steady_flux_continuous_across_interfaces() {
    let q = 5000.0;
    let layers = vec![
        ThermalLayer {
            material: GOLD,
            thickness: 200e-9,
        },
        ThermalLayer {
            material: FUSED_SILICA,
            thickness: 1e-3,
        },
        ThermalLayer {
            material: WATER,
            thickness: 2e-3,
        },
    ];
    let model = ThermalModel::new(layers, 30, Boundary::Flux(q), Boundary::Fixed(0.0));
    let source = vec![0.0; model.n()];
    let t = model.steady_state(&source);
    let mut worst = 0.0f64;
    // Interior face fluxes g[i](T[i-1]-T[i]); front flux is q.
    for i in 1..model.n() {
        let flux = model.g[i] * (t[i - 1] - t[i]);
        worst = worst.max((flux - q).abs() / q);
    }
    println!("steady flux continuity: worst relative deviation {worst:.3e}");
    assert!(worst < 1e-6, "flux not continuous: {worst:e}");
}

/// A photonic absorbed-power profile feeds the thermal solver, and the
/// integrated source recovers `I0 * A`.
#[test]
fn photonic_source_integrates_to_absorptance() {
    use photonics::dispersion::{Dispersion, DrudeLorentz};
    use photonics::multilayer::{Layer, Polarization, Stack};

    let lambda = 633e-9;
    let gold = DrudeLorentz::gold().index(lambda);
    let film = Layer::new(gold.re, gold.im.abs(), 50e-9);
    let stack = Stack::new(1.33, vec![film], 1.5); // water -> gold -> glass
    let sol = stack.solve(lambda, 0.0, Polarization::TE);
    let i0 = 1.0e5; // W/m^2

    // Thermal grid through the gold film plus a glass substrate.
    let cells = 256;
    let model = ThermalModel::new(
        vec![
            ThermalLayer {
                material: GOLD,
                thickness: 50e-9,
            },
            ThermalLayer {
                material: PYREX,
                thickness: 1e-3,
            },
        ],
        cells,
        Boundary::Adiabatic,
        Boundary::Fixed(293.15),
    );
    let z_gold: Vec<f64> = model.z.clone();
    let src = source_from_optics(&sol, i0, &z_gold);
    let absorbed_thermal: f64 = src.iter().zip(&model.dz).map(|(&p, &d)| p * d).sum();
    let absorbed_optics = i0 * sol.layer_absorptance[0];
    println!(
        "photonics->thermal: absorbed {absorbed_thermal:.6e} W/m^2 vs I0*A {absorbed_optics:.6e}, rel {:.3e}",
        (absorbed_thermal - absorbed_optics).abs() / absorbed_optics
    );
    // Sampling at cell centres across the gold film plus a coarse substrate is
    // second-order; require a few percent.
    close(absorbed_thermal / absorbed_optics, 1.0, 0.05);

    let t = model.steady_state(&src);
    let peak = t.iter().cloned().fold(f64::MIN, f64::max);
    println!("  steady peak temperature rise {:.3e} K", peak - 293.15);
    assert!(peak.is_finite());
}

// ---------------------------------------------------------------------------
// 2. Photoacoustics
// ---------------------------------------------------------------------------

/// The Grüneisen parameter of water from its definition.
#[test]
fn gruneisen_water_value() {
    let beta = 2.07e-4; // volumetric thermal expansion, 1/K (CRC, 20 C)
    let c = 1482.0;
    let cp = 4182.0;
    let g = gruneisen(beta, c, cp);
    println!("Gamma_water = {g:.6} (quoted round value 0.11)");
    close(g, GRUNEISEN_WATER, 2e-4);
    assert!((g - 0.11).abs() < 0.01);
}

/// Free surface gives a bipolar signal, rigid surface a monopolar one.
#[test]
fn free_bipolar_rigid_monopolar() {
    let mu_a = 1.0e4;
    let f = 0.02;
    let c = 1482.0;
    let p0 = |z: f64| initial_pressure(mu_a, f, GRUNEISEN_WATER, z);
    let z_obs = 3.0 / mu_a;
    let mut min_free = f64::MAX;
    let mut min_rigid = f64::MAX;
    let mut max_free = f64::MIN;
    for k in 0..=2000 {
        let t = 3.0e-6 * k as f64 / 2000.0;
        let pf = dalembert(p0, z_obs, t, c, Surface::Free);
        let pr = dalembert(p0, z_obs, t, c, Surface::Rigid);
        min_free = min_free.min(pf);
        max_free = max_free.max(pf);
        min_rigid = min_rigid.min(pr);
    }
    println!(
        "free surface: range [{min_free:.3e}, {max_free:.3e}] Pa (sign change: {})",
        min_free < 0.0
    );
    println!(
        "rigid surface: minimum {min_rigid:.3e} Pa (monopolar: {})",
        min_rigid >= 0.0
    );
    assert!(
        min_free < 0.0 && max_free > 0.0,
        "free surface should be bipolar"
    );
    assert!(min_rigid >= -1e-12, "rigid surface should be monopolar");

    // Rigid surface pressure is the exact exponential decay.
    for &t in &[1.0e-7f64, 5.0e-7, 1.0e-6] {
        let want = rigid_surface_pressure(mu_a, f, GRUNEISEN_WATER, c, t);
        let got = dalembert(p0, 0.0, t, c, Surface::Rigid);
        close(got, want, 1e-9);
    }
}

/// Staggered FDTD against the exact d'Alembert field, for a smooth Gaussian
/// initial pressure (rigid surface) and the physical exponential (both
/// surfaces).
#[test]
fn fdtd_matches_dalembert() {
    let c = 1482.0;
    let rho = 998.0;
    // --- smooth Gaussian, rigid surface: clean convergence
    let sigma = 3.0e-5;
    let z0 = 20.0 * sigma;
    let total = 80.0 * sigma;
    let dz = sigma / 15.0;
    let model = AcousticModel::from_layers(
        &[AcousticLayer::new(rho, c, total)],
        dz,
        Surface::Rigid,
        FarBoundary::Absorbing,
    );
    let p0 = |z: f64| (-((z - z0) / sigma).powi(2)).exp();
    let dt = model.max_dt(0.9);
    let (p, v) = model.initial_pressure_profile(p0, dt);
    // Observe before either travelling wave reaches a boundary.
    let t_obs = 10.0 * sigma / c;
    let steps = (t_obs / dt).round() as usize;
    let (p, _) = model.run(p, v, dt, steps);
    let tf = steps as f64 * dt;
    let (mut s2, mut r2) = (0.0, 0.0);
    for i in 0..model.n {
        let want = dalembert(p0, model.z[i], tf, c, Surface::Rigid);
        s2 += (p[i] - want).powi(2);
        r2 += want * want;
    }
    let l2 = (s2 / r2).sqrt();
    println!("FDTD Gaussian vs d'Alembert (rigid): relative L2 {l2:.3e} at t={tf:.3e} s");
    assert!(l2 < 5.0e-3, "Gaussian FDTD L2 {l2:e}");

    // --- physical exponential, rigid surface (even extension is smooth)
    let mu_a = 1.0e4;
    let f = 0.02;
    let delta = 1.0 / mu_a;
    let total = 20.0 * delta;
    let model = AcousticModel::from_layers(
        &[AcousticLayer::new(rho, c, total)],
        delta / 80.0,
        Surface::Rigid,
        FarBoundary::Absorbing,
    );
    let p0 = |z: f64| initial_pressure(mu_a, f, GRUNEISEN_WATER, z);
    let dt = model.max_dt(0.9);
    let (p, v) = model.initial_pressure_profile(p0, dt);
    let t_obs = 0.4 * total / c;
    let steps = (t_obs / dt).round() as usize;
    let (p, _) = model.run(p, v, dt, steps);
    let tf = steps as f64 * dt;
    let (mut s2, mut r2) = (0.0, 0.0);
    for i in 0..model.n {
        let want = dalembert(p0, model.z[i], tf, c, Surface::Rigid);
        s2 += (p[i] - want).powi(2);
        r2 += want * want;
    }
    let l2 = (s2 / r2).sqrt();
    println!("FDTD exponential vs d'Alembert (rigid): relative L2 {l2:.3e} at t={tf:.3e} s");
    // The even extension of the exponential has a kink in its third derivative
    // at the origin, so the exact field carries a moving kink and the FDTD
    // error falls only first order in `dz`; the value is small at 80 cells per
    // absorption depth.
    assert!(l2 < 1.0e-2, "exponential rigid-surface FDTD L2 {l2:e}");

    // --- smooth profile compatible with the free surface, p0(0)=0
    // The physical exponential has a pressure jump at a free surface, so the
    // exact solution is discontinuous and the FDTD smears it; use a smooth
    // odd-compatible profile to test the free boundary itself.
    let model = AcousticModel::from_layers(
        &[AcousticLayer::new(rho, c, total)],
        delta / 20.0,
        Surface::Free,
        FarBoundary::Absorbing,
    );
    let p0 = |z: f64| {
        let x = z / delta;
        x * (-x).exp()
    };
    let dt = model.max_dt(0.9);
    let (p, v) = model.initial_pressure_profile(p0, dt);
    let t_obs = 0.4 * total / c;
    let steps = (t_obs / dt).round() as usize;
    let (p, _) = model.run(p, v, dt, steps);
    let tf = steps as f64 * dt;
    let (mut s2, mut r2) = (0.0, 0.0);
    for i in 0..model.n {
        let want = dalembert(p0, model.z[i], tf, c, Surface::Free);
        s2 += (p[i] - want).powi(2);
        r2 += want * want;
    }
    let l2 = (s2 / r2).sqrt();
    println!("FDTD smooth vs d'Alembert (free): relative L2 {l2:.3e} at t={tf:.3e} s");
    assert!(l2 < 1.0e-2, "smooth free-surface FDTD L2 {l2:e}");
}

/// Frequency-domain launch through a phononic transfer matrix. A layer with the
/// same impedance as its surroundings is a pure delay `d/c`.
#[test]
fn spectral_transmit_matches_pure_delay() {
    use phononic::layered::{Layer, Stack};
    use phononic::materials::WATER as PWATER;

    let c = PWATER.c;
    let n = 4096;
    let f0 = 1.0e6; // 1 MHz
    let dt = 1.0 / (f0 * 40.0); // several samples per period
    let burst = tone_burst(n, dt, f0, 0.4);
    // Choose the layer thickness so the delay is exactly an integer number of
    // samples; then the pure-delay transfer function must reproduce a circular
    // shift to round-off.
    let shift: isize = 135;
    let d = shift as f64 * dt * c;
    let stack = Stack::new(
        PWATER.impedance(),
        vec![Layer::lossless(PWATER.rho, c, d)],
        PWATER.impedance(),
    );
    let out = spectral_transmit(&burst, dt, Some(&stack));
    let delay = d / c;
    let shift = (delay / dt).round() as isize;
    // Compare the transmitted pulse with the input shifted by `shift` samples.
    let mut s2 = 0.0;
    let mut r2 = 0.0;
    for i in 0..n {
        let j = i as isize - shift;
        let want = if j >= 0 && (j as usize) < n {
            burst[j as usize]
        } else {
            0.0
        };
        s2 += (out[i] - want).powi(2);
        r2 += want * want;
    }
    let l2 = (s2 / r2).sqrt();
    println!("spectral water layer: delay {delay:.3e} s = {shift} samples, relative L2 {l2:.3e}");
    assert!(l2 < 1.0e-6, "spectral delay L2 {l2:e}");

    // FFT round trip is exact.
    let back = ifft_real(&fft_real(&burst));
    let e = burst
        .iter()
        .zip(&back)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0f64, f64::max);
    assert!(e < 1e-10);
}

// ---------------------------------------------------------------------------
// 3. Emergent-boundary indicators
// ---------------------------------------------------------------------------

#[test]
fn confinement_clear_cases() {
    let alpha = 1.5e-7; // water
    let c = 1482.0;
    let delta = 1.0e-5; // 10 micron absorption depth
                        // Thermal relaxation of the absorption depth: delta^2/alpha = 0.67 ms.
    let thermal_time = delta * delta / alpha;
    let acoustic_time = delta / c;
    println!(
        "delta={delta} m: thermal time {thermal_time:.3e} s, acoustic time {acoustic_time:.3e} s"
    );
    // 1 ns pulse is confined both ways; 1 s is not.
    let tau_short = 1.0e-9;
    let tau_long = 1.0;
    assert!(is_thermally_confined(tau_short, alpha, delta));
    assert!(is_stress_confined(tau_short, c, delta));
    assert!(!is_thermally_confined(tau_long, alpha, delta));
    assert!(!is_stress_confined(tau_long, c, delta));
    println!(
        "tau=1 ns: thermal ratio {:.3e}, stress ratio {:.3e}",
        thermal_confinement_ratio(tau_short, alpha, delta),
        stress_confinement_ratio(tau_short, c, delta)
    );
}

#[test]
fn coherence_boundary() {
    // Optical thickness / coherence length.
    assert_eq!(optical_regime(1.5, 1e-6, 1e-3), OpticalRegime::Coherent);
    assert_eq!(optical_regime(1.5, 1e-3, 1e-3), OpticalRegime::Incoherent);
    let ind = coherence_indicator(1.5, 1e-3, 1e-3);
    println!("coherence indicator n d/L_c = {ind}");
    assert!((ind - 1.5).abs() < 1e-12);
}

#[test]
fn rayleigh_clear_cases() {
    let g = 9.81;
    let beta = 2.07e-4;
    let nu = 1.0e-6;
    let alpha = 1.4e-7;
    let l = 1.0e-2;
    let ra_cold = rayleigh_benard(g, beta, 0.01, l, nu, alpha);
    let ra_hot = rayleigh_benard(g, beta, 1.0, l, nu, alpha);
    println!("Ra(dT=0.01 K) = {ra_cold:.1}, Ra(dT=1 K) = {ra_hot:.1}, Ra_c = 1708");
    assert!(!convects_rayleigh_benard(ra_cold));
    assert!(convects_rayleigh_benard(ra_hot));

    // Internal heating, both branches of the absorbed-heating selector.
    let q = 1.0e4; // W/m^3
    let k = 0.613;
    let ra_i = rayleigh_internal(g, beta, q, l, k, nu, alpha);
    println!("Ra_I(q=1e4) = {ra_i:.3}, Ra_I,c = 2772");
    assert!(ra_i > 0.0);
    // Surface-like branch when mu_a L >> 1.
    let (ra_s, crit_s) = absorbed_heating_rayleigh(g, beta, 1.0, q, l, k, nu, alpha, 1.0e5);
    assert_eq!(crit_s, RA_CRIT_RAYLEIGH_BENARD);
    assert!(ra_s > 0.0);
    // Volumetric branch when mu_a L << 1.
    let (ra_v, crit_v) = absorbed_heating_rayleigh(g, beta, 1.0, q, l, k, nu, alpha, 1.0e-1);
    assert_eq!(crit_v, RA_CRIT_INTERNAL);
    close(ra_v, ra_i, 1e-12);
}

/// Frequency-domain launch across a single interface: the transmitted pulse is
/// the incident pulse scaled by the exact normal-incidence pressure
/// transmission coefficient `t = 2 Z2 / (Z1 + Z2)` (frequency independent).
#[test]
fn spectral_fresnel_interface_transmission() {
    use phononic::layered::Stack;
    use phononic::materials::{STEEL, WATER as PWATER};

    let z1 = PWATER.impedance();
    let z2 = STEEL.impedance_p();
    let stack = Stack::new(z1, vec![], z2);
    let n = 4096;
    let f0 = 1.0e6;
    let dt = 1.0 / (f0 * 40.0);
    let burst = tone_burst(n, dt, f0, 0.3);
    let out = spectral_transmit(&burst, dt, Some(&stack));
    let t = 2.0 * z2 / (z1 + z2);
    let mut worst = 0.0f64;
    for i in 0..n {
        worst = worst.max((out[i] - t * burst[i]).abs());
    }
    println!(
        "Fresnel interface: t = {t:.9}, worst |out - t*in| = {worst:.3e} (peak {:.3e})",
        burst.iter().cloned().fold(f64::MIN, f64::max)
    );
    assert!(
        worst < 1.0e-9,
        "Fresnel spectral transmission error {worst:e}"
    );
}

/// Strong water absorption at 10.6 micron: the absorption coefficient from the
/// complex index, the absorption depth, and Beer-Lambert volumetric heating
/// with exact energy conservation.
#[test]
fn water_absorption_10um() {
    // Hale & Querry (1973), Appl. Opt. 12, 555: water at 10.6 micron has
    // n = 1.18, kappa = 0.050 (the CO2-laser line).
    let lambda = 10.6e-6;
    let n_w = 1.18;
    let kappa = 0.050;
    let mu_a = 4.0 * PI * kappa / lambda;
    let delta = absorption_depth(mu_a);
    println!("water 10.6 um: n={n_w}, kappa={kappa}, mu_a={mu_a:.4e} 1/m, delta={delta:.4e} m");
    assert!((mu_a - 5.93e4).abs() / 5.93e4 < 0.01);
    assert!((delta * 1e6 - 16.9).abs() < 0.2);
    let _ = n_w;

    // Beer-Lambert heating of a water slab with adiabatic ends stores exactly
    // the absorbed energy.
    let thickness = 5.0 * delta;
    let i0 = 1.0e5;
    let model = ThermalModel::new(
        vec![ThermalLayer {
            material: WATER,
            thickness,
        }],
        100,
        Boundary::Adiabatic,
        Boundary::Adiabatic,
    );
    let source: Vec<f64> = (0..model.n())
        .map(|i| {
            let z0 = model.z[i] - 0.5 * model.dz[i];
            beer_lambert_cell_average(mu_a, i0, z0, model.dz[i])
        })
        .collect();
    let absorbed: f64 = source.iter().zip(&model.dz).map(|(&p, &d)| p * d).sum();
    let want_absorbed = i0 * (1.0 - (-mu_a * thickness).exp());
    println!(
        "  absorbed {absorbed:.6e} W/m^2 vs exact {want_absorbed:.6e} W/m^2, rel {:.3e}",
        rel(absorbed, want_absorbed)
    );
    close(absorbed / want_absorbed, 1.0, 1e-9);

    let dt = 1e-5;
    let steps = 2000;
    let mut t = vec![0.0; model.n()];
    model.run(&mut t, &source, dt, steps);
    let stored = model.stored_energy(&t);
    let deposited = model.source_energy(&source, dt * steps as f64);
    close(stored / deposited, 1.0, 1e-7);
    println!("  stored/deposited = {:.9}", stored / deposited);

    // Photoacoustic initial pressure from the same absorption.
    let fluence = 0.01; // J/m^2
    let p_peak = initial_peak_pressure(mu_a, fluence, GRUNEISEN_WATER);
    println!("  photoacoustic peak p0 = {p_peak:.3e} Pa");
    assert!(p_peak > 0.0);
}

/// Antireflection-coated window: air | MgF2 quarter wave | glass (faintly
/// absorbing) | air. The absorbed profile integrates to `I0 A`, and the ARC
/// keeps the absorptance small.
#[test]
fn arc_coated_window_source() {
    use photonics::multilayer::{Layer, Polarization, Stack};

    let lambda = 550e-9;
    let n_s = 1.52;
    let n_c = 1.38; // MgF2
    let arc = Layer::quarter_wave(n_c, lambda);
    // Keep the coherent transfer matrix in its range of validity: a 2 micron
    // glass pane (a few wavelengths), with a small but nonzero absorption.
    let glass = Layer::new(n_s, 1.0e-4, 2.0e-6);
    let stack = Stack::new(1.0, vec![arc, glass], 1.0);
    let sol = stack.solve(lambda, 0.0, Polarization::TE);
    let i0 = 1.0e3;
    println!(
        "ARC window: R={:.3e}, T={:.3e}, A={:.3e}",
        sol.reflectance, sol.transmittance, sol.absorptance
    );
    assert!(sol.absorptance > 0.0 && sol.absorptance < 0.05);
    assert!(
        sol.reflectance < 0.15,
        "ARC should suppress the front reflection"
    );

    // Thermal grid coincident with the optical layers; the integrated source
    // recovers I0*A to sampling accuracy.
    let model = ThermalModel::new(
        vec![
            ThermalLayer {
                material: MGF2,
                thickness: arc.thickness,
            },
            ThermalLayer {
                material: PYREX,
                thickness: glass.thickness,
            },
        ],
        64,
        Boundary::Fixed(293.15),
        Boundary::Fixed(293.15),
    );
    let src = source_from_optics(&sol, i0, &model.z);
    let absorbed: f64 = src.iter().zip(&model.dz).map(|(&p, &d)| p * d).sum();
    let want = i0 * sol.absorptance;
    println!(
        "  integrated source {absorbed:.6e} W/m^2 vs I0*A {want:.6e}, rel {:.3e}",
        rel(absorbed, want)
    );
    close(absorbed / want, 1.0, 0.05);
}

/// Full energy ledger with outflow: deposited = stored + conducted out the
/// cold back face.
#[test]
fn energy_ledger_with_outflow() {
    let q = 8.0e3;
    let model = ThermalModel::new(
        vec![ThermalLayer {
            material: WATER,
            thickness: 4e-3,
        }],
        200,
        Boundary::Flux(q),
        Boundary::Fixed(293.15),
    );
    let n = model.n();
    let source = vec![0.0; n];
    let mut t = vec![293.15; n];
    let stored0 = model.stored_energy(&t);
    let dt = 1e-3;
    let steps = 4000;
    let tf = dt * steps as f64;
    let mut conducted_out = 0.0;
    for _ in 0..steps {
        model.step(&mut t, &source, dt);
        conducted_out += -model.back_flux(&t) * dt;
    }
    let deposited = q * tf;
    let stored = model.stored_energy(&t) - stored0;
    let balance = (deposited - stored - conducted_out).abs() / deposited;
    println!(
        "energy ledger: deposited {deposited:.6e}, stored {stored:.6e}, conducted {conducted_out:.6e}, residual {balance:.2e}"
    );
    assert!(balance < 1e-6, "energy ledger residual {balance:e}");
}
