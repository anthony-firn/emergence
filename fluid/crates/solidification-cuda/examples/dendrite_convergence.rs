//! 2D free-dendrite convergence example (GPU).
//!
//! Benchmark: Karma-Rappel anisotropic phase field at Delta = 0.55,
//! eps4 = 0.05 (Karma & Rappel, Phys. Rev. E 57, 4323 (1998), Table II).
//! The solver is a quarter domain with symmetry on x=0, y=0 and far field on
//! the outer edges, progressively translated in -x so the tip stays in a
//! small box (the moving-frame technique of KR, Sec. III).
//!
//! Usage:
//!   dendrite_convergence D w0 d0 dx nx ny steps state_file
//! Environment:
//!   SEED_R   initial quarter-disk radius (default 3*w0)
//!   REPORT   report interval in steps (default 2000)
//!   SHIFT    if non-zero, translate when tip passes SHIFT_AT*Lx (default 1)
//!   SHIFT_AT tip position (fraction of Lx) that triggers a translation (default 0.55)
//!   SHIFT_TO tip position (fraction of Lx) after the translation (default 0.25)
//!   DIAG     print far-field diagnostics at each report (default 0)
//!   LABEL    label printed in the report lines
//!   PREC     "f32" (default, fast) or "f64"
//!   DT       override the time step (default = stable_dt)
//!   TAU_FACTOR multiply tau0 after construction (kinetic-free sensitivity, default 1)
//!   MAXSEC   wall-clock budget; stops early and saves (default none)
//!
//! State is written after every report; re-running with the same `state_file`
//! resumes and accumulates. This keeps every invocation under a few minutes.

use solidification::phase_field::{Bc, PhaseField};
use solidification_cuda::{HostView, PhaseFieldGpu, Real};
use std::io::{Read, Write};
use std::time::Instant;

const EPS4: f64 = 0.05;
const DELTA: f64 = 0.55;

#[derive(Clone)]
struct Header {
    nx: usize,
    ny: usize,
    dx: f64,
    d0: f64,
    d: f64,
    w0: f64,
    u_far: f64,
    shift: f64,
    t: f64,
    /// Diagnostic only: multiplies the constructed `tau0` (not persisted).
    tau_factor: f64,
}

fn save(path: &str, h: &Header, p: &[f64], u: &[f64]) -> std::io::Result<()> {
    let tmp = format!("{path}.tmp");
    let mut f = std::fs::File::create(&tmp)?;
    f.write_all(b"PFD1")?;
    for v in [h.nx as f64, h.ny as f64, h.dx, h.d0, h.d, h.w0, h.u_far, h.shift, h.t] {
        // tau_factor is diagnostic and deliberately not persisted
        f.write_all(&v.to_le_bytes())?;
    }
    for v in p {
        f.write_all(&v.to_le_bytes())?;
    }
    for v in u {
        f.write_all(&v.to_le_bytes())?;
    }
    drop(f);
    std::fs::rename(&tmp, path)
}

fn load(path: &str) -> Option<(Header, Vec<f64>, Vec<f64>)> {
    let mut f = std::fs::File::open(path).ok()?;
    let mut magic = [0u8; 4];
    f.read_exact(&mut magic).ok()?;
    if &magic != b"PFD1" {
        return None;
    }
    let mut vals = [0.0f64; 9];
    for v in vals.iter_mut() {
        let mut b = [0u8; 8];
        f.read_exact(&mut b).ok()?;
        *v = f64::from_le_bytes(b);
    }
    let h = Header {
        nx: vals[0] as usize,
        ny: vals[1] as usize,
        dx: vals[2],
        d0: vals[3],
        d: vals[4],
        w0: vals[5],
        u_far: vals[6],
        shift: vals[7],
        t: vals[8],
        tau_factor: 1.0,
    };
    let n = (h.nx + 2) * (h.ny + 2);
    let mut p = vec![0.0; n];
    let mut u = vec![0.0; n];
    for v in p.iter_mut() {
        let mut b = [0u8; 8];
        f.read_exact(&mut b).ok()?;
        *v = f64::from_le_bytes(b);
    }
    for v in u.iter_mut() {
        let mut b = [0u8; 8];
        f.read_exact(&mut b).ok()?;
        *v = f64::from_le_bytes(b);
    }
    Some((h, p, u))
}

fn build(h: &Header) -> PhaseField {
    let mut pf = PhaseField::new(h.nx, h.ny, h.dx, h.d0, h.d, h.w0, EPS4);
    pf.bc = Bc::quadrant();
    pf.u_far = h.u_far;
    pf.tau0 *= h.tau_factor;
    pf
}

/// Translate the phase and temperature fields `sc` cells in -x, filling the
/// trailing columns with liquid at the far-field temperature. Ghosts are
/// refreshed so the device restart is clean.
fn translate(h: &Header, p: &mut [f64], u: &mut [f64], sc: usize) {
    let m = h.nx + 2;
    for j in 0..(h.ny + 2) {
        let base = j * m;
        for i in 1..=(h.nx - sc) {
            p[base + i] = p[base + i + sc];
            u[base + i] = u[base + i + sc];
        }
        for i in (h.nx - sc + 1)..=h.nx {
            p[base + i] = -1.0;
            u[base + i] = h.u_far;
        }
        // left edge (symmetry)
        p[base] = p[base + 1];
        u[base] = u[base + 1];
        p[base + h.nx + 1] = -1.0;
        u[base + h.nx + 1] = h.u_far;
    }
}

/// Crossing of `p = level` going outward on row `j`, or `None`.
fn crossing_level(hw: &HostView, j: usize, level: f64) -> Option<f64> {
    let nx = hw.nx;
    let at = |i: usize| hw.p[j * (nx + 2) + i];
    let mut prev = at(1);
    for i in 2..=nx {
        let cur = at(i);
        if prev >= level && cur < level {
            let x0 = (i as f64 - 1.5) * hw.dx;
            let x1 = (i as f64 - 0.5) * hw.dx;
            let t = (prev - level) / (prev - cur);
            return Some(x0 + t * (x1 - x0));
        }
        prev = cur;
    }
    None
}

/// Tip radius from a least-squares parabola fit `x(y) = x0 - y^2/(2R)` to the
/// `p = 0` needle contour.
///
/// The diffuse interface rounds the very tip on the scale `W0`, so the first
/// ~`1.5 W0` of the contour (where the local radius of curvature is only a few
/// `W0`) is excluded; the fit then runs over the single-valued monotone needle
/// out to `40 W0`. The cap matters: the selected dendrite is only
/// approximately parabolic and fitting the flanks of a very long needle
/// overshoots `R` by ~25 %. `HostView::tip_radius` fits the `p = 0.5` contour
/// over a fixed 8 rows, which at `dx/W0 = 0.4` is only `2.8 W0` and sits inside
/// the rounded cap: it underestimates `R` by ~2.4x and was the source of the
/// reported `sigma* ~ 0.5-0.7`. See `analysis/solidification_cuda.md` §5.
fn robust_tip_radius(hw: &HostView, w0: f64) -> f64 {
    let y0 = 1.5 * w0;
    let y1 = 40.0 * w0;
    let mut ys = Vec::new();
    let mut xs = Vec::new();
    let mut prev = f64::INFINITY;
    for j in 1..=hw.ny {
        let Some(x) = crossing_level(hw, j, 0.0) else {
            break;
        };
        let y = (j as f64 - 0.5) * hw.dx;
        if y < y0 {
            prev = x;
            continue;
        }
        if y > y1 {
            break;
        }
        if x > prev + 0.25 * hw.dx {
            break;
        }
        ys.push(y * y);
        xs.push(x);
        prev = x;
    }
    if ys.len() < 4 {
        return f64::NAN;
    }
    let n = ys.len() as f64;
    let sy: f64 = ys.iter().sum();
    let sx: f64 = xs.iter().sum();
    let syy: f64 = ys.iter().map(|v| v * v).sum();
    let sxy: f64 = ys.iter().zip(xs.iter()).map(|(a, b)| a * b).sum();
    let den = n * syy - sy * sy;
    if den.abs() < 1e-18 {
        return f64::NAN;
    }
    let b = (n * sxy - sy * sx) / den;
    if b >= -1e-12 {
        return f64::NAN;
    }
    -1.0 / (2.0 * b)
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let mut h = Header {
        d: a.get(1).map(|v| v.parse().unwrap()).unwrap_or(2.0),
        w0: a.get(2).map(|v| v.parse().unwrap()).unwrap_or(1.0),
        d0: a.get(3).map(|v| v.parse().unwrap()).unwrap_or(1.0),
        dx: a.get(4).map(|v| v.parse().unwrap()).unwrap_or(0.4),
        nx: a.get(5).map(|v| v.parse().unwrap()).unwrap_or(300),
        ny: a.get(6).map(|v| v.parse().unwrap()).unwrap_or(200),
        u_far: -DELTA,
        shift: 0.0,
        t: 0.0,
        tau_factor: 1.0,
    };
    let steps: usize = a.get(7).map(|v| v.parse().unwrap()).unwrap_or(20000);
    let state = a.get(8).cloned().unwrap_or_default();
    let seed_r: f64 = std::env::var("SEED_R").ok().and_then(|s| s.parse().ok()).unwrap_or(3.0 * h.w0);
    let report: usize = std::env::var("REPORT").ok().and_then(|s| s.parse().ok()).unwrap_or(2000);
    let shift_on = std::env::var("SHIFT").map(|v| v != "0").unwrap_or(true);
    let shift_at: f64 = std::env::var("SHIFT_AT").ok().and_then(|s| s.parse().ok()).unwrap_or(0.55);
    let shift_to: f64 = std::env::var("SHIFT_TO").ok().and_then(|s| s.parse().ok()).unwrap_or(0.25);
    let diag = std::env::var("DIAG").map(|v| v != "0").unwrap_or(false);
    let label = std::env::var("LABEL").unwrap_or_else(|_| "run".into());
    let f32run = std::env::var("PREC").map(|v| v == "f32").unwrap_or(true);
    let mode: i32 = std::env::var("STENCIL").ok().and_then(|v| v.parse().ok()).unwrap_or(3);
    let maxsec: f64 = std::env::var("MAXSEC").ok().and_then(|s| s.parse().ok()).unwrap_or(1e18);
    let dt_override: Option<f64> = std::env::var("DT").ok().and_then(|s| s.parse().ok());
    let tau_factor: f64 = std::env::var("TAU_FACTOR").ok().and_then(|s| s.parse().ok()).unwrap_or(1.0);

    if f32run {
        run::<f32>(&mut h, steps, &state, seed_r, report, shift_on, shift_at, shift_to, diag, &label, maxsec, dt_override, mode, tau_factor);
    } else {
        run::<f64>(&mut h, steps, &state, seed_r, report, shift_on, shift_at, shift_to, diag, &label, maxsec, dt_override, mode, tau_factor);
    }
}

#[allow(clippy::too_many_arguments)]
fn run<R: Real>(
    h: &mut Header,
    steps: usize,
    state: &str,
    seed_r: f64,
    report: usize,
    shift_on: bool,
    shift_at: f64,
    shift_to: f64,
    diag: bool,
    label: &str,
    maxsec: f64,
    dt_override: Option<f64>,
    mode: i32,
    tau_factor: f64,
) {
    let loaded = if state.is_empty() { None } else { load(state) };
    h.tau_factor = tau_factor;
    let (pf, mut p_host, mut u_host) = match loaded {
        Some((hh, pp, uu)) => {
            *h = hh;
            h.tau_factor = tau_factor;
            (build(h), pp, uu)
        }
        None => {
            let mut pf = build(h);
            match std::env::var("SEED_PARAB").ok().and_then(|s| s.parse::<f64>().ok()) {
                Some(rtip) => {
                    let xtip: f64 = std::env::var("SEED_X").ok().and_then(|s| s.parse().ok()).unwrap_or(15.0);
                    pf.seed_parabola(xtip, rtip);
                }
                None => pf.seed_corner(seed_r),
            }
            (pf, Vec::new(), Vec::new())
        }
    };
    let mut gpu = PhaseFieldGpu::<R>::from_cpu(&pf).unwrap();
    if !p_host.is_empty() {
        gpu.upload_p(&p_host).unwrap();
        gpu.upload_u(&u_host).unwrap();
    }
    let dt = dt_override.unwrap_or_else(|| gpu.stable_dt());
    println!(
        "[{label}] D={} w0={} d0={} W0/d0={:.4} dx={} W0/dx={:.3} Nx={} Ny={} dt={dt:.5e} prec={} stencil={mode} t0={:.2} shift={:.2}",
        h.d, h.w0, h.d0, h.w0 / h.d0, h.dx, h.w0 / h.dx, h.nx, h.ny, R::SUFFIX, h.t, h.shift
    );

    let lx = h.nx as f64 * h.dx;
    let mut last_tip = {
        let p = gpu.download_p().unwrap();
        let hw = HostView { nx: h.nx, ny: h.ny, dx: h.dx, p: &p };
        hw.crossing_x(1).unwrap_or(0.0)
    };
    let mut last_t = h.t;
    let start = Instant::now();
    for s in 0..steps {
        gpu.step_kr_mode(dt, mode).unwrap();
        h.t += dt;
        let done = s + 1;
        if done % report != 0 && done != steps {
            continue;
        }
        let p = gpu.download_p().unwrap();
        let u = gpu.download_u().unwrap();
        let hw = HostView { nx: h.nx, ny: h.ny, dx: h.dx, p: &p };
        let tip = hw.crossing_x(1).unwrap_or(f64::NAN);
        let r = robust_tip_radius(&hw, h.w0);
        let v = (tip - last_tip) / (h.t - last_t);
        let vd = v * h.d0 / h.d;
        let sel = v * r * r / (2.0 * h.d * h.d0);
        let mut extra = String::new();
        if diag {
            // Far-field diagnostics: the temperature on the last interior column
            // and row (should stay at u_far) and the measured diffusion length.
            let m = h.nx + 2;
            let u_right = u[(1) * m + h.nx];
            let u_top = u[h.ny * m + 1];
            let l_diff = if v.abs() > 1e-9 { 2.0 * h.d / v } else { f64::INFINITY };
            extra = format!(
                " uR={u_right:+.5} uT={u_top:+.5} (ufar={:+.5}) l=2D/V={l_diff:6.2} Ly={:6.1} ahead={:6.1}",
                h.u_far,
                h.ny as f64 * h.dx,
                (h.nx as f64 * h.dx) - tip
            );
        }
        println!(
            "[{label}] step={done:>8} t={:9.2} tip={tip:8.3} lab={:9.3} v={v:9.5} Vd0/D={vd:9.5} R={r:8.3} sel={sel:7.3}{extra}",
            h.t,
            tip + h.shift
        );
        last_tip = tip;
        last_t = h.t;
        p_host = p;
        u_host = u;
        // Persist for restart.
        if !state.is_empty() {
            save(state, h, &p_host, &u_host).unwrap();
        }
        // Move-frame translation: keep the tip near 0.25 Lx.
        if shift_on && tip > shift_at * lx {
            let sc = ((tip - shift_to * lx) / h.dx).round() as usize;
            let sc = sc.clamp(1, h.nx / 2);
            translate(h, &mut p_host, &mut u_host, sc);
            h.shift += sc as f64 * h.dx;
            last_tip -= sc as f64 * h.dx;
            gpu.upload_p(&p_host).unwrap();
            gpu.upload_u(&u_host).unwrap();
            if !state.is_empty() {
                save(state, h, &p_host, &u_host).unwrap();
            }
        }
        if start.elapsed().as_secs_f64() > maxsec {
            println!("[{label}] wall budget reached after {done} steps");
            break;
        }
    }
    if state.is_empty() {
        return;
    }
    // Final persist (the last report already saved, but do it for steps<report).
    let p = gpu.download_p().unwrap();
    let u = gpu.download_u().unwrap();
    save(state, h, &p, &u).unwrap();
    gpu.synchronize().unwrap();
}
