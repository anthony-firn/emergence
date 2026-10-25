//! Longer production runs for the molecular-dynamics rung.
//!
//! `cargo test -p md` uses short trajectories (a few thousand steps). This
//! example runs a longer single trajectory at the Lennard-Jones triple
//! point and prints the Green–Kubo viscosity and thermal conductivity and
//! the self-diffusion coefficient, so the estimates and their scatter can
//! be inspected without slowing the test suite down.
//!
//! Run with e.g. `cargo run --release -p md --example long_runs`.
//! Optional first argument: number of production steps (default 60000).

use materials::LennardJones;

use md::transport::{
    autocorrelate, correlation_time, diffusion_from_msd, diffusion_from_msd_blocks,
    diffusion_from_vacf, finite_size_corrected_d, green_kubo_conductivity, green_kubo_viscosity,
    msd, vacf,
};
use md::validity::{continuum_indicator, FlowScales};
use md::{Cutoff, Potential, System, Thermostat};

fn main() {
    let prod: usize = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(60_000);
    let rho = 0.8442;
    let t = 0.722;
    let dt = 0.003; // Meier et al. (2004) used dt* = 0.003
    let pot = Potential::LennardJones(LennardJones::new(1.0, 1.0));
    let mut sys = System::fcc(5, rho, pot, 3.0, 0.3, Cutoff::TailCorrected, 1234);
    sys.set_velocities_mb(t);
    let n = sys.n();
    let mut thermo = Thermostat::bussi(0.5);
    let eq = 5000;
    for _ in 0..eq {
        thermo.step(&mut sys, dt, t);
    }
    let every = 2usize;
    let frame_every = 10usize;
    let mut stress = Vec::new();
    let mut heat = Vec::new();
    let mut pos = Vec::new();
    let mut vel = Vec::new();
    let mut tsum = 0.0;
    for step in 0..prod {
        thermo.step(&mut sys, dt, t);
        if step % every == 0 {
            stress.push(sys.pressure_tensor());
            heat.push(sys.heat_flux());
            tsum += sys.temperature_com();
        }
        if step % frame_every == 0 {
            pos.extend_from_slice(&sys.pos);
            vel.extend_from_slice(&sys.vel);
        }
    }
    let vol = sys.volume();
    let kbt = tsum / stress.len() as f64;
    let dt_s = dt * every as f64;
    let frame_dt = dt * frame_every as f64;

    let (eta, eta_scatter) = green_kubo_viscosity(&stress, dt_s, vol, kbt, 4000, 300, 1200);
    let (kappa, kappa_scatter) = green_kubo_conductivity(&heat, dt_s, vol, kbt, 2000, 100, 500);
    let ms = msd(&pos, n, 600);
    let (d_msd, _) = diffusion_from_msd(&ms, frame_dt, 10, 200);
    let (d_blocks, d_blocks_err) = diffusion_from_msd_blocks(&pos, n, frame_dt, 10, 200, 8);
    let vc = vacf(&vel, n, 600);
    let (d_vacf, _) = diffusion_from_vacf(&vc, frame_dt, 200);
    let d0 = finite_size_corrected_d(d_msd, kbt, eta, sys.box_l);
    let mut cxy = autocorrelate(&stress.iter().map(|s| s[3]).collect::<Vec<_>>(), 4000);
    cxy.dt = dt_s;
    let tau = correlation_time(&cxy);

    println!("N = {n}, rho* = {rho}, T* = {kbt:.4} (target {t}), {prod} steps");
    println!("eta*     = {eta:.3} +/- {eta_scatter:.3}   (Meier 2004: 3.258 +/- 0.033)");
    println!("kappa*   = {kappa:.3} +/- {kappa_scatter:.3}   (argon-equivalent 6.9)");
    println!("D_pbc*   = {d_msd:.5} (MSD) vs {d_vacf:.5} (VACF); block error {d_blocks_err:.5}");
    println!("D0*      = {d0:.5} (Yeh-Hummer finite-size corrected)");
    println!("D blocks = {d_blocks:.5}");
    println!("tau_stress = {tau:.4} reduced time units");

    // Continuum validity for a millimetre-scale flow at 1 m/s in argon.
    let sigma = md::argon::SIGMA_M;
    let flow = FlowScales {
        length: 1e-3 / sigma,
        velocity: 1.0 / md::argon::velocity_unit(),
        temperature: kbt,
        mass: 1.0,
    };
    let ind = continuum_indicator(flow, d_msd, tau);
    println!("continuum indicator (L=1 mm, U=1 m/s): {ind:?}");
    println!("Navier-Stokes with these closures holds: {}", ind.holds());
}
