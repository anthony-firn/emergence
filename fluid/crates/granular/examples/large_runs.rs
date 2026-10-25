//! Longer production runs for the granular / suspension rung.
//!
//! `cargo test -p granular` uses a 108-particle shear cell for 12 000 steps
//! and two-sphere collisions. This example runs
//!
//! 1. a larger constant-pressure shear cell (256 particles, more rates) and
//!    fits the `mu(I)` closure; and
//! 2. quasi-2D granular heaps at two interparticle friction coefficients to
//!    measure the angle of repose.
//!
//! Run with `cargo run --release -p granular --example large_runs`. The
//! command is bounded by `timeout` in normal use; the defaults finish in
//! well under two minutes.

use granular::contact::{hertz_mindlin, ContactParams, ContactState};
use granular::dem::{Particle, ShearCell};
use granular::rheology::Jop2006;
use granular::rng::Rng;

fn contact_params(mu: f64, restitution: f64) -> ContactParams {
    ContactParams {
        young: 1000.0,
        poisson: 0.3,
        mu,
        mu_rolling: 0.15,
        restitution,
    }
}

fn shear_cell_sweep() {
    println!("=== constant-pressure shear cell: mu(I) ===");
    let rates = [0.005_f64, 0.01, 0.02, 0.04, 0.08, 0.16];
    let mut data = Vec::new();
    for (k, &g) in rates.iter().enumerate() {
        let mut cell = ShearCell::new_cubic(
            150,
            0.58,
            g,
            contact_params(0.4, 0.5),
            0.5,
            1.0,
            0.005,
            500 + k as u64,
        );
        cell.set_pressure(0.05, 2.0);
        let (p, tau, mu, i, tg) = cell.run(8000);
        println!(
            "  gdot={g:.4} P={p:.4} tau={tau:.4} mu={mu:.4} I={i:.4} Tg={tg:.2e} phi={:.4}",
            cell.volume_fraction()
        );
        data.push((i, mu));
    }
    let (fit, rmse) = Jop2006::fit(&data);
    let glass = Jop2006::glass_beads();
    println!(
        "  fit: mu_s={:.4} mu_2={:.4} I0={:.4} rmse={rmse:.4}",
        fit.mu_s, fit.mu_2, fit.i0
    );
    println!(
        "  Jop (2006) glass: mu_s={:.4} mu_2={:.4} I0={:.4}",
        glass.mu_s, glass.mu_2, glass.i0
    );
}

/// One flat-wall contact. `n` points from the wall into the fluid, `gap` is
/// the signed normal distance from the particle centre to the wall plane
/// (`gap < r` for contact).
#[allow(clippy::too_many_arguments)]
fn wall_force(
    p: &ContactParams,
    r: f64,
    m: f64,
    v: [f64; 3],
    w: [f64; 3],
    n: [f64; 3],
    gap: f64,
    state: &mut ContactState,
    dt: f64,
) -> ([f64; 3], [f64; 3]) {
    let delta = r - gap;
    if delta <= 0.0 {
        *state = ContactState::default();
        return ([0.0; 3], [0.0; 3]);
    }
    let e_star = p.e_star();
    let g_star = p.young / (2.0 * (2.0 - p.poisson) * (1.0 + p.poisson));
    let beta = p.beta();
    let sqrt_r = r.sqrt();
    let sqrt_d = delta.sqrt();
    let f_elastic = (4.0 / 3.0) * e_star * sqrt_r * delta * sqrt_d;
    let s_n = 2.0 * e_star * sqrt_r * sqrt_d;
    let eta_n = 2.0 * (5.0 / 6.0_f64).sqrt() * beta * (s_n * m).sqrt();
    let v_n = v[0] * n[0] + v[1] * n[1] + v[2] * n[2];
    let fn_mag = f_elastic - eta_n * v_n;

    // Contact point relative to the particle centre: c = -r n.
    let c = [-r * n[0], -r * n[1], -r * n[2]];
    let wxc = [
        w[1] * c[2] - w[2] * c[1],
        w[2] * c[0] - w[0] * c[2],
        w[0] * c[1] - w[1] * c[0],
    ];
    let vs = [v[0] + wxc[0], v[1] + wxc[1], v[2] + wxc[2]];
    let vs_n = vs[0] * n[0] + vs[1] * n[1] + vs[2] * n[2];
    let vt = [
        vs[0] - vs_n * n[0],
        vs[1] - vs_n * n[1],
        vs[2] - vs_n * n[2],
    ];
    let mut xi = state.xi;
    for k in 0..3 {
        xi[k] += vt[k] * dt;
    }
    let xin = xi[0] * n[0] + xi[1] * n[1] + xi[2] * n[2];
    for k in 0..3 {
        xi[k] -= xin * n[k];
    }
    let s_t = 8.0 * g_star * sqrt_r * sqrt_d;
    let eta_t = 2.0 * (5.0 / 6.0_f64).sqrt() * beta * (s_t * m).sqrt();
    let mut ft = [0.0; 3];
    for k in 0..3 {
        ft[k] = -s_t * xi[k] - eta_t * vt[k];
    }
    let ft_mag = (ft[0] * ft[0] + ft[1] * ft[1] + ft[2] * ft[2]).sqrt();
    let limit = p.mu * fn_mag.max(0.0);
    if ft_mag > limit && ft_mag > 0.0 {
        let sc = limit / ft_mag;
        for k in 0..3 {
            ft[k] *= sc;
            xi[k] = -ft[k] / s_t;
        }
    }
    state.xi = xi;
    let force = [
        fn_mag * n[0] + ft[0],
        fn_mag * n[1] + ft[1],
        fn_mag * n[2] + ft[2],
    ];
    let torque = [
        c[1] * ft[2] - c[2] * ft[1],
        c[2] * ft[0] - c[0] * ft[2],
        c[0] * ft[1] - c[1] * ft[0],
    ];
    (force, torque)
}

/// Quasi-2D heap: gravity, a rough floor and two side walls, periodic in
/// `z`. Returns the angle of repose estimated from the two free surfaces.
fn heap_angle(mu: f64, mu_roll: f64, n: usize, steps: usize, seed: u64) -> f64 {
    let r = 0.5_f64;
    let density = 1.0;
    let (lx, lz) = (40.0_f64, 2.0_f64);
    let dt = 0.003_f64;
    let g = 0.5_f64;
    let mut p = contact_params(mu, 0.4);
    p.mu_rolling = mu_roll;
    let mut rng = Rng::new(seed);
    // Dense rectangular block that collapses into a pile.
    let cols = 22usize;
    let rows = (n / cols).max(1);
    let mut parts: Vec<Particle> = Vec::with_capacity(cols * rows);
    let x0 = lx * 0.5 - 0.5 * cols as f64 * 1.05;
    for i in 0..cols {
        for j in 0..rows {
            let x = [
                x0 + i as f64 * 1.05 + rng.uniform(-0.01, 0.01),
                0.5 + j as f64 * 1.0 + rng.uniform(-0.01, 0.01),
                rng.uniform(0.0, lz),
            ];
            parts.push(Particle::new(x, r, density));
        }
    }
    let n = parts.len();
    let mut contacts: std::collections::HashMap<(u32, u32), ContactState> =
        std::collections::HashMap::new();
    let mut wall_state: Vec<ContactState> = vec![ContactState::default(); n * 3];

    for _ in 0..steps {
        // Pair forces (z periodic).
        let mut f = vec![[0.0; 3]; n];
        let mut t = vec![[0.0; 3]; n];
        let mut new_contacts = std::collections::HashMap::with_capacity(contacts.len() + 16);
        for i in 0..n {
            for j in (i + 1)..n {
                let xi = parts[i].x;
                let xj = parts[j].x;
                let mut dr = [xj[0] - xi[0], xj[1] - xi[1], xj[2] - xi[2]];
                let nz = (dr[2] / lz).round();
                dr[2] -= nz * lz;
                let dist2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
                if dist2 >= 4.0 * r * r {
                    continue;
                }
                let xju = [xi[0] + dr[0], xi[1] + dr[1], xi[2] + dr[2]];
                let key = (i as u32, j as u32);
                let mut st = *contacts.get(&key).unwrap_or(&ContactState::default());
                if let Some(c) = hertz_mindlin(
                    &p,
                    r,
                    r,
                    parts[i].m,
                    parts[j].m,
                    xi,
                    parts[i].v,
                    parts[i].w,
                    xju,
                    parts[j].v,
                    parts[j].w,
                    &mut st,
                    dt,
                ) {
                    for k in 0..3 {
                        f[i][k] += c.force_i[k];
                        f[j][k] -= c.force_i[k];
                        t[i][k] += c.torque_i[k];
                        t[j][k] += c.torque_j[k];
                    }
                    new_contacts.insert(key, st);
                }
            }
        }
        contacts = new_contacts;

        // Gravity + walls + integration (semi-implicit Euler).
        for i in 0..n {
            let mut force = f[i];
            let mut torque = t[i];
            force[1] -= parts[i].m * g;
            let x = parts[i].x;
            let v = parts[i].v;
            let w = parts[i].w;
            // Floor.
            let (wf, wt) = wall_force(
                &p,
                r,
                parts[i].m,
                v,
                w,
                [0.0, 1.0, 0.0],
                x[1],
                &mut wall_state[i * 3],
                dt,
            );
            for k in 0..3 {
                force[k] += wf[k];
                torque[k] += wt[k];
            }
            // Left wall.
            let (wf, wt) = wall_force(
                &p,
                r,
                parts[i].m,
                v,
                w,
                [1.0, 0.0, 0.0],
                x[0],
                &mut wall_state[i * 3 + 1],
                dt,
            );
            for k in 0..3 {
                force[k] += wf[k];
                torque[k] += wt[k];
            }
            // Right wall.
            let (wf, wt) = wall_force(
                &p,
                r,
                parts[i].m,
                v,
                w,
                [-1.0, 0.0, 0.0],
                lx - x[0],
                &mut wall_state[i * 3 + 2],
                dt,
            );
            for k in 0..3 {
                force[k] += wf[k];
                torque[k] += wt[k];
            }
            let m = parts[i].m;
            let inertia = parts[i].inertia();
            for k in 0..3 {
                parts[i].v[k] += force[k] / m * dt;
                parts[i].w[k] += torque[k] / inertia * dt;
                parts[i].x[k] += parts[i].v[k] * dt;
            }
            // Keep z in the box.
            if parts[i].x[2] >= lz {
                parts[i].x[2] -= lz;
            } else if parts[i].x[2] < 0.0 {
                parts[i].x[2] += lz;
            }
        }
    }

    // Surface profile: maximum height per x-bin.
    let nb = 80usize;
    let dxb = lx / nb as f64;
    let mut h = vec![0.0_f64; nb];
    for part in &parts {
        let b = ((part.x[0] / dxb) as usize).min(nb - 1);
        h[b] = h[b].max(part.x[1]);
    }
    // Find the apex and fit each flank from slightly below the apex down to
    // the point where the surface reaches one particle diameter.
    let (mut apex, mut hmax) = (0usize, 0.0);
    for (b, &v) in h.iter().enumerate() {
        if v > hmax {
            hmax = v;
            apex = b;
        }
    }
    let mut angle_sum = 0.0;
    let mut count = 0;
    for &(lo, hi, rev) in &[(0usize, apex, true), (apex + 1, nb, false)] {
        let idxs: Vec<usize> = if rev {
            (lo..hi).rev().collect()
        } else {
            (lo..hi).collect()
        };
        let mut xs = Vec::new();
        let mut ys = Vec::new();
        for b in idxs {
            if h[b] < 1.0 {
                break;
            }
            if h[b] < 0.8 * hmax {
                xs.push((b as f64 + 0.5) * dxb);
                ys.push(h[b]);
            }
        }
        if xs.len() >= 3 {
            let n = xs.len() as f64;
            let sx: f64 = xs.iter().sum();
            let sy: f64 = ys.iter().sum();
            let sxx: f64 = xs.iter().map(|x| x * x).sum();
            let sxy: f64 = xs.iter().zip(&ys).map(|(x, y)| x * y).sum();
            let slope = (n * sxy - sx * sy) / (n * sxx - sx * sx);
            angle_sum += slope.abs().atan().to_degrees();
            count += 1;
        }
    }
    angle_sum / count.max(1) as f64
}

fn main() {
    let mode = std::env::args().nth(1).unwrap_or_else(|| "all".into());
    if mode != "heap" {
        shear_cell_sweep();
    }
    if mode != "shear" {
        println!("\n=== angle of repose (quasi-2D heap, sliding friction 0.5) ===");
        for &mr in &[0.0_f64, 0.2, 0.4] {
            let angle = heap_angle(0.5, mr, 176, 40000, 11 + (mr * 100.0) as u64);
            println!("  mu_rolling={mr:.2} -> angle of repose {angle:.1} deg");
        }
    }
}
