//! GPU Green–Kubo size and cutoff study at the Lennard-Jones triple point.
//!
//! Usage:
//!   cargo run -p md-cuda --release --example gk_study -- <k> <rc> <steps> <seeds> [tag]
//!
//! `N = 4 k^3` atoms on an fcc lattice at `rho* = 0.8442`, `T* = 0.722`
//! (Meier, Laesecke & Kabelac 2004). `<seeds>` is a comma-separated list
//! of RNG seeds; each seed is an independent trajectory and the quoted
//! error is the standard error across seeds (plus, for a single long
//! trajectory, the block error from `running_gk_viscosity`). The scheme is
//! `Cutoff::TailCorrected` (bare truncation plus the analytic long-range
//! corrections; the tail is isotropic and does not enter the off-diagonal
//! stress that sets the shear viscosity). Each system is melted and
//! equilibrated with periodic velocity rescaling, then the production
//! trajectory is NVE. The f64-accumulated pressure tensor is sampled every
//! `sample_every` steps on the device, and unwrapped positions are copied
//! back every `frame_every` steps for the self-diffusion coefficient.
//!
//! Results go to `<tag>.csv` (running GK curve, seed-averaged) and
//! `<tag>.txt` (summary) under `$MDVISC_OUT` (default `/tmp/mdvisc`).

use md::potential::{Cutoff, Potential};
use md::system::System;
use md::transport::{diffusion_from_msd, msd};
use md_cuda::gk::{
    average_curves, einstein_helfand_viscosity, mean_sem, plateau_peak_lag,
    running_gk_viscosity, window_mean, BlockCurve,
};
use md_cuda::GpuMd;
use std::fs;
use std::time::Instant;

fn lj() -> Potential {
    Potential::LennardJones(materials::LennardJones::new(1.0, 1.0))
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.len() < 4 {
        eprintln!("usage: gk_study <k> <rc> <steps> <seeds> [tag]");
        std::process::exit(2);
    }
    let k: usize = args[0].parse().unwrap();
    let rc: f64 = args[1].parse().unwrap();
    let steps: usize = args[2].parse().unwrap();
    let seeds: Vec<u64> = args[3]
        .split(',')
        .map(|s| s.parse().unwrap())
        .collect();
    let tag = args.get(4).cloned().unwrap_or_else(|| {
        format!("k{k}_rc{rc}_n{steps}_nseed{}", seeds.len())
    });

    if !md_cuda::gpu_available() {
        eprintln!("no CUDA device; skipping {tag}");
        return;
    }

    let rho = 0.8442;
    let t_target = 0.722;
    let dt = 0.004;
    let sample_every = 2usize;
    let frame_every = (steps / 400).max(1);
    let max_lag = 1500usize; // 1500 * 2 * dt = 12 reduced time units
    let nsamp = steps.div_ceil(sample_every) + 2;

    let t0 = Instant::now();
    let mut seed_press: Vec<Vec<[f64; 6]>> = Vec::new();
    let mut seed_kbt: Vec<f64> = Vec::new();
    let mut seed_t: Vec<f64> = Vec::new();
    let mut seed_eta_eh: Vec<f64> = Vec::new();
    let mut seed_d: Vec<f64> = Vec::new();
    let mut seed_d_short: Vec<f64> = Vec::new();
    let mut n = 0usize;
    let mut box_l = 0.0f64;
    let mut nframes = 0usize;

    for &seed in &seeds {
        let cpu = System::fcc(k, rho, lj(), rc, 0.3, Cutoff::TailCorrected, seed);
        n = cpu.n();
        box_l = cpu.box_l;
        let mut gpu =
            GpuMd::<f32>::new(&cpu.pos, &cpu.vel, box_l, rc, &cpu.potential, cpu.scheme).unwrap();
        // Melt the metastable fcc, then settle at the target state.
        for s in 0..1500usize {
            gpu.step(dt).unwrap();
            if s % 20 == 0 {
                gpu.rescale_to_temperature(1.5).unwrap();
            }
        }
        for s in 0..5000usize {
            gpu.step(dt).unwrap();
            if s % 20 == 0 {
                gpu.rescale_to_temperature(t_target).unwrap();
            }
        }
        let t_eq = gpu.temperature().unwrap();

        gpu.reserve_pressure_samples(nsamp).unwrap();
        nframes = steps / frame_every + 1;
        let mut frames: Vec<[f64; 3]> = Vec::with_capacity(nframes * n);
        let mut tsum = 0.0;
        let mut tcount = 0usize;
        for s in 0..steps {
            gpu.step(dt).unwrap();
            if s % sample_every == 0 {
                gpu.sample_pressure_tensor().unwrap();
            }
            if s % frame_every == 0 {
                tsum += gpu.temperature().unwrap();
                tcount += 1;
                frames.extend_from_slice(&gpu.positions().unwrap());
            }
        }
        let press = gpu.take_pressure_samples().unwrap();
        seed_press.push(press);
        let kbt = tsum / tcount as f64;
        seed_kbt.push(kbt);
        seed_t.push(t_eq);

        // Per-seed Einstein–Helfand estimate (slope over the plateau).
        let vol = box_l.powi(3);
        let dt_s = dt * sample_every as f64;
        let (eta_eh, _) =
            einstein_helfand_viscosity(seed_press.last().unwrap(), dt_s, vol, kbt, max_lag, 1, 100, max_lag / 2);
        seed_eta_eh.push(eta_eh);

        // Per-seed self-diffusion. Two fit windows are reported: the whole
        // diffusive window, and the short, CPU-`md`-compatible window
        // t in [0.3, 6.0] (the `long_runs` example uses [0.3, 6.0]).
        let frame_dt = dt * frame_every as f64;
        let ms = msd(&frames, n, nframes.saturating_sub(1));
        let (d_msd, _) = diffusion_from_msd(&ms, frame_dt, 5, nframes / 2);
        let lo_s = ((0.3 / frame_dt).round() as usize).max(1);
        let hi_s = ((6.0 / frame_dt).round() as usize)
            .min(nframes.saturating_sub(1))
            .max(lo_s + 1);
        let (d_short, _) = diffusion_from_msd(&ms, frame_dt, lo_s, hi_s);
        seed_d.push(d_msd);
        seed_d_short.push(d_short);
        println!(
            "  seed {seed}: T_eq={t_eq:.4} T_prod={kbt:.4} eta_EH={eta_eh:.3} D_long={d_msd:.5} D_short={d_short:.5}"
        );
    }
    let elapsed = t0.elapsed().as_secs_f64();
    let vol = box_l.powi(3);
    let dt_s = dt * sample_every as f64;
    let _frame_dt = dt * frame_every as f64;

    // Seed-averaged running GK curves. Each seed is one full trajectory
    // (nblocks = 1); the spread across seeds is the honest error bar.
    let curves: Vec<BlockCurve> = seed_press
        .iter()
        .zip(&seed_kbt)
        .map(|(p, &kbt)| running_gk_viscosity(p, dt_s, vol, kbt, max_lag, 1))
        .collect();
    let avg = average_curves(&curves, dt_s);
    // Locate the plateau as the peak of the running integral smoothed over
    // ~0.4 reduced time units (the stress correlation time is ~0.14), then
    // average each seed's curve over a +/-1 reduced-time window around it.
    // The per-seed window means give an honest seed-averaged error bar.
    let peak_lag = plateau_peak_lag(&avg, ((0.4 / dt_s).round() as usize).max(1));
    let half = ((1.0 / dt_s).round() as usize).max(1);
    let w_lo = peak_lag.saturating_sub(half);
    let w_hi = (peak_lag + half).min(max_lag);
    let per_seed_plateau: Vec<f64> = curves.iter().map(|c| window_mean(c, w_lo, w_hi)).collect();
    let (eta_gk, eta_gk_err) = mean_sem(&per_seed_plateau);
    let plateau_time = peak_lag as f64 * dt_s;
    let plateau_lo = w_lo as f64 * dt_s;
    let plateau_hi = w_hi as f64 * dt_s;
    // Cross-check: a fixed physical integration window t in [1, 4] reduced
    // time, averaged per seed with the same seed-SEM error.
    let f_lo = ((1.0 / dt_s).round() as usize).min(max_lag);
    let f_hi = ((4.0 / dt_s).round() as usize).min(max_lag);
    let per_seed_fixed: Vec<f64> = curves.iter().map(|c| window_mean(c, f_lo, f_hi)).collect();
    let (eta_gk_fixed, eta_gk_fixed_err) = mean_sem(&per_seed_fixed);
    // Value of the averaged curve at a fixed 2-time-unit cutoff, for
    // comparison with the earlier md-crate short runs.
    let fixed_lag = (2.0 / dt_s).round() as usize;
    let fixed_lag = fixed_lag.min(avg.mean.len() - 1);
    let (eta_gk_2, eta_gk_2_err) = (avg.mean[fixed_lag], avg.err[fixed_lag]);

    let (eta_eh, eta_eh_err) = mean_sem(&seed_eta_eh);
    let (d_m, d_e) = mean_sem(&seed_d);
    let (d_s, d_se) = mean_sem(&seed_d_short);
    let kbt_mean = seed_kbt.iter().sum::<f64>() / seed_kbt.len() as f64;
    let d0 = md::transport::finite_size_corrected_d(d_m, kbt_mean, eta_gk.max(1e-6), box_l);
    let d0_short = md::transport::finite_size_corrected_d(d_s, kbt_mean, eta_gk.max(1e-6), box_l);

    let out_dir = std::env::var("MDVISC_OUT").unwrap_or_else(|_| "/tmp/mdvisc".into());
    fs::create_dir_all(&out_dir).ok();
    let mut csv = String::from("lag,time,eta,err\n");
    for (i, (&m, &e)) in avg.mean.iter().zip(&avg.err).enumerate() {
        csv.push_str(&format!("{i},{:.6},{m:.8},{e:.8}\n", i as f64 * dt_s));
    }
    fs::write(format!("{out_dir}/{tag}.csv"), csv).ok();

    let seeds_str = seeds
        .iter()
        .map(|s| s.to_string())
        .collect::<Vec<_>>()
        .join(",");
    let t_eq_str = seed_t
        .iter()
        .map(|t| format!("{t:.4}"))
        .collect::<Vec<_>>()
        .join(",");
    let summary = format!(
        "tag={tag} N={n} k={k} rc={rc} steps={steps} seeds={seeds_str} dt={dt} sample_every={sample_every} frame_every={frame_every}\n\
         T_eq=[{t_eq_str}] T_prod={kbt_mean:.5} volume={vol:.6} elapsed_s={elapsed:.1}\n\
         eta_GK={eta_gk:.5} eta_GK_err={eta_gk_err:.5} plateau_lag={peak_lag} plateau_time={plateau_time:.5} plateau_window=[{plateau_lo:.4},{plateau_hi:.4}]\n\
         eta_GK_fixed_1_4={eta_gk_fixed:.5} eta_GK_fixed_1_4_err={eta_gk_fixed_err:.5}\n\
         eta_GK_at_2={eta_gk_2:.5} eta_GK_at_2_err={eta_gk_2_err:.5}\n\
         eta_EH={eta_eh:.5} eta_EH_err={eta_eh_err:.5}\n\
         D_msd={d_m:.6} D_msd_err={d_e:.6} D0={d0:.6}\n\
         D_short={d_s:.6} D_short_err={d_se:.6} D0_short={d0_short:.6}\n\
         nframes={nframes}\n"
    );
    fs::write(format!("{out_dir}/{tag}.txt"), &summary).ok();
    print!("{summary}");
    println!("wrote {out_dir}/{tag}.csv and {out_dir}/{tag}.txt");
}
