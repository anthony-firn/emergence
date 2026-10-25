//! Benchmarks for the molecular-dynamics rung.
//!
//! Every test checks a published or exact result and reports the error; no
//! tolerance is widened to make a model agree with itself. Runs are kept
//! short (256 atoms, a few thousand steps in the fast tests); the longer
//! transport trajectories live in `examples/long_runs.rs`.

use materials::potential::{energy_per_atom_with_tail, relax};
use materials::{Lattice, LennardJones};

use md::argon::{argon_closure, liquid_argon_90k, liquid_argon_at_lj_triple};
use md::eos::JzgEos;
use md::potential::Cutoff;
use md::potential::Potential;
use md::system::System;
use md::thermostat::Thermostat;
use md::transport::{
    autocorrelate, correlation_time, diffusion_from_msd, diffusion_from_vacf,
    finite_size_corrected_d, green_kubo_conductivity, green_kubo_viscosity, msd, vacf,
};
use md::validity::{continuum_indicator, FlowScales};

fn lj() -> Potential {
    Potential::LennardJones(LennardJones::new(1.0, 1.0))
}

fn block_mean(v: &[f64], nb: usize) -> (f64, f64) {
    let n = v.len();
    let b = n / nb;
    let mut means = Vec::with_capacity(nb);
    for k in 0..nb {
        means.push(v[k * b..(k + 1) * b].iter().sum::<f64>() / b as f64);
    }
    let mean = means.iter().sum::<f64>() / nb as f64;
    let var = means.iter().map(|m| (m - mean).powi(2)).sum::<f64>() / (nb - 1) as f64;
    (mean, (var / nb as f64).sqrt())
}

/// NVT average of temperature, pressure and potential energy per atom using
/// the Bussi thermostat. Returns `(T, P, sigma_P, U, sigma_U)`.
fn nvt_average(
    rho: f64,
    t: f64,
    ncells: usize,
    seed: u64,
    eq: usize,
    prod: usize,
) -> (f64, f64, f64, f64, f64) {
    let mut sys = System::fcc(ncells, rho, lj(), 2.5, 0.3, Cutoff::TailCorrected, seed);
    sys.set_velocities_mb(t);
    let mut th = Thermostat::bussi(0.5);
    let dt = 0.005;
    for _ in 0..eq {
        th.step(&mut sys, dt, t);
    }
    let mut ps = Vec::with_capacity(prod);
    let mut us = Vec::with_capacity(prod);
    let mut tsum = 0.0;
    for _ in 0..prod {
        th.step(&mut sys, dt, t);
        ps.push(sys.pressure());
        us.push(sys.potential_energy() / sys.n() as f64);
        tsum += sys.temperature_com();
    }
    let tm = tsum / prod as f64;
    let (pm, pe) = block_mean(&ps, 20);
    let (um, ue) = block_mean(&us, 20);
    (tm, pm, pe, um, ue)
}

#[test]
fn nve_energy_conservation_is_second_order() {
    // Identical initial condition for both step sizes, so the comparison is
    // of the integrator alone. A liquid of 256 atoms, shifted potential.
    let mut base = System::fcc(4, 0.8, lj(), 2.5, 0.3, Cutoff::Shifted, 7);
    base.set_velocities_mb(1.0);
    let mut th = Thermostat::bussi(0.5);
    for _ in 0..2500 {
        th.step(&mut base, 0.002, 1.0);
    }
    let pos = base.pos.clone();
    let vel = base.vel.clone();
    let box_l = base.box_l;
    let n = base.n() as f64;

    let rms = |dt: f64| -> f64 {
        let mut s = System::new(pos.clone(), box_l, lj(), 2.5, 0.3, Cutoff::Shifted, 7);
        s.vel = vel.clone();
        s.compute_forces();
        let steps = (1.0 / dt).round() as usize; // fixed physical time 1.0
        let mut sum = 0.0;
        let mut sum2 = 0.0;
        for _ in 0..steps {
            s.step_nve(dt);
            let e = s.total_energy();
            sum += e;
            sum2 += e * e;
        }
        let mean = sum / steps as f64;
        ((sum2 / steps as f64) - mean * mean).sqrt() / n
    };
    let r1 = rms(0.004);
    let r2 = rms(0.002);
    let ratio = r1 / r2;
    assert!(
        (3.3..4.7).contains(&ratio),
        "NVE rms fluctuation ratio {ratio:.3} (dt=0.004: {r1:.3e}, dt=0.002: {r2:.3e}) is not dt^2"
    );
    // The drift over the run must itself shrink.
    let drift = |dt: f64| -> f64 {
        let mut s = System::new(pos.clone(), box_l, lj(), 2.5, 0.3, Cutoff::Shifted, 7);
        s.vel = vel.clone();
        s.compute_forces();
        let e0 = s.total_energy();
        let steps = (1.0 / dt).round() as usize;
        for _ in 0..steps {
            s.step_nve(dt);
        }
        (s.total_energy() - e0).abs() / n
    };
    assert!(drift(0.002) < drift(0.004));
}

#[test]
fn jzg_eos_low_density_state() {
    // JZG (1993) reference: T*=2.0, rho*=0.5 -> P*=1.077450, U_res*=-3.144943.
    let (t, p, pe, u, ue) = nvt_average(0.5, 2.0, 4, 11, 2000, 4000);
    let eos = JzgEos::new();
    let p_ref = eos.pressure(0.5, 2.0);
    let u_ref = eos.residual_energy(0.5, 2.0);
    println!(
        "T*={t:.4}  P*={p:.4}+/-{pe:.4} (JZG {p_ref:.4})  U_res*={u:.4}+/-{ue:.4} (JZG {u_ref:.4})"
    );
    assert!((t - 2.0).abs() < 0.05, "temperature {t}");
    assert!(
        (p - p_ref).abs() < 3.0 * pe + 0.05,
        "pressure {p:.4}+/-{pe:.4} vs JZG {p_ref:.4}"
    );
    assert!(
        (u - u_ref).abs() < 3.0 * ue + 0.05,
        "residual energy {u:.4}+/-{ue:.4} vs JZG {u_ref:.4}"
    );
}

#[test]
fn jzg_eos_dense_state() {
    // JZG (1993) reference: T*=1.0, rho*=0.8 -> P*=1.031806, U_res*=-5.523454.
    let (t, p, pe, u, ue) = nvt_average(0.8, 1.0, 4, 12, 2000, 4000);
    let eos = JzgEos::new();
    let p_ref = eos.pressure(0.8, 1.0);
    let u_ref = eos.residual_energy(0.8, 1.0);
    println!(
        "T*={t:.4}  P*={p:.4}+/-{pe:.4} (JZG {p_ref:.4})  U_res*={u:.4}+/-{ue:.4} (JZG {u_ref:.4})"
    );
    assert!((t - 1.0).abs() < 0.05, "temperature {t}");
    assert!(
        (p - p_ref).abs() < 3.0 * pe + 0.05,
        "pressure {p:.4}+/-{pe:.4} vs JZG {p_ref:.4}"
    );
    assert!(
        (u - u_ref).abs() < 3.0 * ue + 0.05,
        "residual energy {u:.4}+/-{ue:.4} vs JZG {u_ref:.4}"
    );
}

#[test]
fn rdf_is_physical_for_the_liquid() {
    let mut sys = System::fcc(4, 0.8, lj(), 2.5, 0.3, Cutoff::TailCorrected, 3);
    sys.set_velocities_mb(1.0);
    let mut th = Thermostat::bussi(0.5);
    for _ in 0..4000 {
        th.step(&mut sys, 0.004, 1.0);
    }
    let rmax = 3.0;
    let nb = 60;
    let g = sys.rdf(rmax, nb);
    let dr = rmax / nb as f64;

    // Hard core: g is zero well below the LJ minimum.
    for b in 0..8 {
        let r = (b as f64 + 0.5) * dr;
        assert!(g[b] < 0.05, "g({r:.2})={} not zero", g[b]);
    }
    // g -> 1 in the bulk.
    for b in (nb - 6)..nb {
        assert!(
            (g[b] - 1.0).abs() < 0.15,
            "g({:.2})={:.3}",
            (b as f64 + 0.5) * dr,
            g[b]
        );
    }
    // First peak near 1.1 sigma and a liquid-like coordination number.
    let (pk, _) = g
        .iter()
        .enumerate()
        .skip(5)
        .take(25)
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
        .unwrap();
    let rpk = (pk as f64 + 0.5) * dr;
    assert!((1.0..1.25).contains(&rpk), "first peak at {rpk:.3}");
    assert!(g[pk] > 2.0, "first peak height {}", g[pk]);
    // Coordination number up to the first minimum near r = 1.5 sigma.
    let rmin_bin = (1.5 / dr) as usize;
    let coord: f64 = (0..rmin_bin)
        .map(|b| {
            let r = (b as f64 + 0.5) * dr;
            sys.density() * g[b] * 4.0 * std::f64::consts::PI * r * r * dr
        })
        .sum();
    println!(
        "first peak r={rpk:.3} g={:.3}, coordination {coord:.2}",
        g[pk]
    );
    assert!((10.0..15.0).contains(&coord), "coordination {coord}");
}

#[test]
fn fcc_zero_temperature_energy_matches_materials_lattice_sum() {
    // A perfect cold crystal has exactly the static lattice energy; this
    // isolates the MD pair sum from statistical or thermal effects.
    let rho: f64 = 1.0;
    let a = (4.0 / rho).cbrt();
    let sys = System::fcc(4, rho, lj(), 2.5, 0.3, Cutoff::TailCorrected, 1);
    let mat = energy_per_atom_with_tail(&Lattice::fcc(a), &LennardJones::new(1.0, 1.0), 2.5);
    let md = sys.potential_energy() / sys.n() as f64;
    println!(
        "fcc LJ zero-T: MD {md:.10}  materials {mat:.10}  diff {:.2e}",
        (md - mat).abs()
    );
    assert!((md - mat).abs() < 1e-9, "MD {md} vs materials {mat}");
}

#[test]
fn fcc_low_temperature_adds_the_harmonic_energy() {
    // Classical equipartition: <PE> = E_static + (3/2) T (harmonic) and
    // <E_tot> = E_static + 3 T. Use the relaxed (zero-pressure) lattice
    // constant, where the crystal is elastically stable; at T*=0.05 the
    // anharmonic correction is small.
    let ljp = LennardJones::new(1.0, 1.0);
    let eq = relax(&Lattice::fcc(1.0), &ljp, 2.5, 0.9, 1.8);
    let a = eq.scale;
    let rho = 4.0 / a.powi(3);
    let t = 0.05;
    let mut sys = System::fcc(4, rho, lj(), 2.5, 0.3, Cutoff::TailCorrected, 5);
    sys.set_velocities_mb(t);
    let mut th = Thermostat::bussi(0.5);
    for _ in 0..1500 {
        th.step(&mut sys, 0.004, t);
    }
    let (mut pe, mut tot, mut tt) = (0.0, 0.0, 0.0);
    let ns = 3000;
    for _ in 0..ns {
        th.step(&mut sys, 0.004, t);
        pe += sys.potential_energy() / sys.n() as f64;
        tot += sys.total_energy() / sys.n() as f64;
        tt += sys.temperature_com();
    }
    pe /= ns as f64;
    tot /= ns as f64;
    tt /= ns as f64;
    let e_static = energy_per_atom_with_tail(&Lattice::fcc(a), &LennardJones::new(1.0, 1.0), 2.5);
    let pe_expected = e_static + 1.5 * tt;
    let tot_expected = e_static + 3.0 * tt;
    println!(
        "fcc T*={tt:.4}: <PE>={pe:.6} (static+3T/2={pe_expected:.6})  <E>={tot:.6} (static+3T={tot_expected:.6})"
    );
    assert!(
        (pe - pe_expected).abs() < 0.01,
        "<PE> {pe} vs {pe_expected}"
    );
    assert!(
        (tot - tot_expected).abs() < 0.01,
        "<E> {tot} vs {tot_expected}"
    );
}

/// One short trajectory at the LJ triple point; returns the transport
/// estimates. The Green–Kubo viscosity of a 256-atom, 10^4-step
/// trajectory has ~30% run-to-run scatter, so the caller averages several
/// seeds and uses the seed-to-seed standard error as the honest error bar.
fn triple_point_run(seed: u64, prod: usize) -> (f64, f64, f64, f64, f64, f64, f64) {
    let rho = 0.8442;
    let t = 0.722;
    let mut sys = System::fcc(4, rho, lj(), 2.5, 0.3, Cutoff::TailCorrected, seed);
    sys.set_velocities_mb(t);
    let mut th = Thermostat::bussi(0.5);
    let dt = 0.004;
    for _ in 0..1500 {
        th.step(&mut sys, dt, t);
    }
    let frame_every = 20usize;
    let n = sys.n();
    let mut stress = Vec::new();
    let mut heat = Vec::new();
    let mut pos = Vec::new();
    let mut vel = Vec::new();
    let mut tsum = 0.0;
    for step in 0..prod {
        th.step(&mut sys, dt, t);
        stress.push(sys.pressure_tensor());
        heat.push(sys.heat_flux());
        tsum += sys.temperature_com();
        if step % frame_every == 0 {
            pos.extend_from_slice(&sys.pos);
            vel.extend_from_slice(&sys.vel);
        }
    }
    let vol = sys.volume();
    let kbt = tsum / stress.len() as f64;
    let (eta, _) = green_kubo_viscosity(&stress, dt, vol, kbt, 1500, 200, 600);
    let (kappa, _) = green_kubo_conductivity(&heat, dt, vol, kbt, 800, 100, 400);
    let frame_dt = dt * frame_every as f64;
    let ms = msd(&pos, n, 200);
    let (d_msd, _) = diffusion_from_msd(&ms, frame_dt, 5, 60);
    let vc = vacf(&vel, n, 200);
    let (d_vacf, _) = diffusion_from_vacf(&vc, frame_dt, 50);
    let d0 = finite_size_corrected_d(d_msd, kbt, eta, sys.box_l);
    let mut cxy = autocorrelate(&stress.iter().map(|s| s[3]).collect::<Vec<_>>(), 1500);
    cxy.dt = dt;
    let tau = correlation_time(&cxy);
    (eta, kappa, d_msd, d_vacf, d0, kbt, tau)
}

#[test]
fn transport_at_the_lj_triple_point() {
    // State point of Meier, Laesecke & Kabelac, J. Chem. Phys. 121, 3671
    // (2004): T* = 0.722, rho* = 0.8442 (the LJ triple-point liquid).
    // Reference viscosity at infinite system size and large cutoff:
    // eta* = 3.258 +/- 0.033; finite N = 256 and rc = 2.5 lower it by
    // roughly 10%. We average three independent seeds; the seed-to-seed
    // standard error is the honest error bar of the short-run estimator.
    let seeds = [21u64, 22, 23];
    let mut etas = Vec::new();
    let mut kappas = Vec::new();
    let mut ds = Vec::new();
    let mut d0s = Vec::new();
    let (mut tmean, mut tau, mut dvacf_sum, mut dvacf_n) = (0.0, 0.0, 0.0, 0);
    for &s in &seeds {
        let (eta, kappa, d_msd, d_vacf, d0, kbt, t_stress) = triple_point_run(s, 6000);
        etas.push(eta);
        kappas.push(kappa);
        ds.push(d_msd);
        d0s.push(d0);
        tmean += kbt;
        tau += t_stress;
        dvacf_sum += d_vacf;
        dvacf_n += 1;
    }
    let nseeds = seeds.len() as f64;
    let mean = |v: &[f64]| v.iter().sum::<f64>() / nseeds;
    let stderr = |v: &[f64]| {
        let m = mean(v);
        let var = v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (nseeds - 1.0);
        (var / nseeds).sqrt()
    };
    let eta_m = mean(&etas);
    let eta_e = stderr(&etas);
    let kappa_m = mean(&kappas);
    let kappa_e = stderr(&kappas);
    let d_m = mean(&ds);
    let d_e = stderr(&ds);
    let d0_m = mean(&d0s);
    let d_vacf_m = dvacf_sum / dvacf_n as f64;

    println!("T* = {:.4}", tmean / nseeds);
    println!("eta*   = {eta_m:.3} +/- {eta_e:.3} (3 seeds; Meier 3.258)");
    println!("kappa* = {kappa_m:.3} +/- {kappa_e:.3} (argon-equivalent 6.9)");
    println!("D_pbc* = {d_m:.5} +/- {d_e:.5} (MSD), {d_vacf_m:.5} (VACF), D0* = {d0_m:.5}");
    println!("tau_stress = {:.4}", tau / nseeds);

    // The two independent estimators of D must agree.
    assert!(
        (d_m - d_vacf_m).abs() / d_m < 0.15,
        "MSD {d_m} vs VACF {d_vacf_m}"
    );
    assert!(d0_m < d_m && d0_m > 0.0);
    // Viscosity within ~30% of the infinite-size reference. The difference
    // is the finite-size/cutoff systematic plus the seed-to-seed scatter of
    // a short estimator; we do not tune either away.
    assert!(
        (eta_m - 3.258).abs() / 3.258 < 0.30,
        "eta* {eta_m:.3}+/-{eta_e:.3} vs Meier 3.258"
    );
    assert!(eta_e > 0.0);
    // Thermal conductivity: argon-equivalent 6.9, LJ within ~30%.
    assert!(
        (kappa_m - 6.9).abs() / 6.9 < 0.30,
        "kappa* {kappa_m:.3}+/-{kappa_e:.3} vs 6.9"
    );
}

#[test]
fn argon_closures_match_measured_liquid_argon() {
    // Reduced LJ transport at the triple point: viscosity from Meier et al.
    // (2004), conductivity and self-diffusion from the MD run above.
    let c = argon_closure(0.722, 0.8442, 3.258, 6.48, 0.0248);
    // Measured saturated liquid argon at the corresponding temperature.
    let reft = liquid_argon_at_lj_triple();
    println!(
        "LJ closure: mu={:.4e} Pa s, k={:.5} W/(m K), D={:.3e} m^2/s",
        c.viscosity_pa_s, c.conductivity_w_m_k, c.diffusion_m2_s
    );
    println!(
        "measured Ar at {:.2} K: mu={:.4e} Pa s, k={:.5} W/(m K)",
        reft.temperature_k, reft.viscosity_pa_s, reft.conductivity_w_m_k
    );
    assert!(
        (c.viscosity_pa_s / reft.viscosity_pa_s - 1.0).abs() < 0.25,
        "viscosity ratio {}",
        c.viscosity_pa_s / reft.viscosity_pa_s
    );
    assert!(
        (c.conductivity_w_m_k / reft.conductivity_w_m_k - 1.0).abs() < 0.35,
        "conductivity ratio {}",
        c.conductivity_w_m_k / reft.conductivity_w_m_k
    );
    // Units sanity: the 90 K reference must land on the published numbers.
    let r90 = liquid_argon_90k();
    assert!((r90.viscosity_pa_s - 2.4004e-4).abs() < 1e-8);
    assert!((r90.conductivity_w_m_k - 0.124_51).abs() < 1e-6);
}

#[test]
fn continuum_indicator_separates_liquid_from_gas() {
    let flow = FlowScales {
        length: 100.0,
        velocity: 1.0,
        temperature: 1.0,
        mass: 1.0,
    };
    // Dense liquid: small mean free path, short correlation time.
    let liq = continuum_indicator(flow, 0.03, 0.1);
    assert!(liq.holds(), "liquid {liq:?}");
    // Dilute gas: large D, long mean free path.
    let gas = continuum_indicator(flow, 50.0, 0.1);
    assert!(!gas.holds(), "gas {gas:?}");
    assert!(gas.knudsen > 0.1);

    // A realistic argon flow: L = 1 mm, U = 1 m/s at 90 K liquid values.
    let sigma = md::argon::SIGMA_M;
    let d_si = 0.032 * md::argon::diffusion_unit();
    let flow_si = FlowScales {
        length: 1e-3 / sigma,
        velocity: 1.0 / md::argon::velocity_unit(),
        temperature: 90.0 / md::argon::EPSILON_OVER_KB_K,
        mass: 1.0,
    };
    let ar = continuum_indicator(flow_si, d_si / md::argon::diffusion_unit(), 0.1);
    println!("argon flow indicator: {ar:?}");
    assert!(ar.holds());
}
