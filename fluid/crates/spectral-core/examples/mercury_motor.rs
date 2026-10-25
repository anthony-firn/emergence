//! Homopolar (Faraday-type) liquid-metal motor against its exact solution.
//!
//!   mercury_motor N [current] [penal|point|ext]
//!
//! Liquid metal fills the annulus r_i < r < r_o around a z-axis (periodic
//! in z, walls by implicit Brinkman penalization). Current I per box height
//! 2π enters at the inner electrode and leaves at the outer one, so in the
//! liquid J = I/(4π² r) r̂; with B = B ẑ the Lorentz force −J B/ρ θ̂ spins
//! the liquid. The motional EMF u_θ B r̂ is a pure gradient (no current), so
//! the steady laminar profile solves ν(u'' + u'/r − u/r²) = K/r,
//! K = I B/(4π²ρ):  u = (K/2ν) r ln r + A r + C/r,  u(r_i) = u(r_o) = 0.
use spectral_core::mhd::{coaxial_source, Mhd};
use spectral_core::surface::{resolution, Extension};
use spectral_core::{fft::Fft3, init, Grid, Penalization, Solver, Surface};

const PI: f64 = std::f64::consts::PI;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let n: usize = a[1].parse().unwrap();
    let current: f64 = a.get(2).map(|v| v.parse().unwrap()).unwrap_or(2.0);
    let (nu, rho, sigma, bz) = (0.05, 1.0, 1.0, 1.0);
    let walls = a.get(3).cloned().unwrap_or("penal".into());
    let (ri, ro) = (0.9, 2.3);
    let g = Grid::new(n);
    let c = [PI, PI];
    let eta = 0.25 * g.dx() * g.dx() / nu;
    let radius = |p: usize| {
        let (i, j) = (p / (n * n), (p / n) % n);
        ((g.x(i) - c[0]).powi(2) + (g.x(j) - c[1]).powi(2)).sqrt()
    };
    let mut s = Solver::new(g.clone(), nu, init::from_fn(&g, |_, _, _| [0.0; 3]));
    match walls.as_str() {
        "penal" => {
            let mut pen = Penalization::sphere_sharp(&g, [PI; 3], 0.1, eta).implicit();
            pen.chi = (0..g.real_len()).map(|p| if radius(p) < ri || radius(p) > ro { 1.0 } else { 0.0 }).collect();
            s.penalization = Some(pen);
        }
        "point" => {
            let h = resolution(&g);
            s.surface = Some(Surface::union(vec![Surface::cylinder_z(c, ri, h), Surface::cylinder_z_inverse(c, ro, h)]));
        }
        _ => {
            let h = resolution(&g);
            let width = 0.75;
            let part = |mut sf: Surface| {
                sf.extension = Some(Extension { width, sigma_rel: 1.0 });
                sf.with_basis(h)
            };
            s.surface = Some(Surface::union(vec![part(Surface::cylinder_z(c, ri, 0.5 * h)), part(Surface::cylinder_z_inverse(c, ro, 0.5 * h))]));
        }
    }
    let mut m = Mhd::uniform(sigma / rho, 1.0 / rho, [0.0, 0.0, bz]);
    m.source = Some(coaxial_source(&g, c, ri, ro, 1.5 * g.dx(), current));
    if walls != "penal" {
        // solids are conductors at rest: mask the fictitious interior flow's EMF
        m.fluid_mask = Some((0..g.real_len()).map(|p| if radius(p) < ri || radius(p) > ro { 0.0 } else { 1.0 }).collect());
    }
    s.mhd = Some(m);

    // exact profile
    let k = current * bz / (4.0 * PI * PI * rho);
    let p = |r: f64| k / (2.0 * nu) * r * r.ln();
    // A r + C/r = −p(r) at r_i, r_o
    let det = ri / ro - ro / ri;
    let aa = (-p(ri) / ro + p(ro) / ri) / det;
    let cc = (-ri * p(ro) + ro * p(ri)) / det;
    let exact = |r: f64| p(r) + aa * r + cc / r;
    let umax = (0..400).map(|q| exact(ri + (ro - ri) * q as f64 / 399.0).abs()).fold(0.0, f64::max);
    let re = umax * (ro - ri) / nu;
    let dt = (0.3 * g.dx() / umax).min(0.5);
    let t_end = 6.0 * (ro - ri).powi(2) / nu / PI.powi(2) * 2.0;
    let t0 = std::time::Instant::now();
    let fft = Fft3::new(&g);
    let uth = |s: &Solver| -> Vec<(f64, f64)> {
        let u: Vec<Vec<f64>> = (0..2)
            .map(|cc| {
                let mut w = s.u[cc].clone();
                let mut o = vec![0.0; g.real_len()];
                fft.inverse_inplace(&mut w, &mut o);
                o
            })
            .collect();
        (0..n * n)
            .map(|ij| {
                let p = ij * n; // l = 0 plane
                let (i, j) = (p / (n * n), (p / n) % n);
                let (x, y) = (g.x(i) - c[0], g.x(j) - c[1]);
                let r = (x * x + y * y).sqrt();
                (r, if r > 0.0 { (-y * u[0][p] + x * u[1][p]) / r } else { 0.0 })
            })
            .collect()
    };
    while s.t < t_end {
        s.step(dt);
    }
    let prof = uth(&s);
    let (mut e2, mut n2, mut emax) = (0.0, 0.0, 0.0f64);
    for &(r, u) in prof.iter().filter(|&&(r, _)| r > ri && r < ro) {
        let ex = exact(r);
        e2 += (u - ex).powi(2);
        n2 += ex * ex;
        emax = emax.max((u - ex).abs());
    }
    if std::env::var("PROFILE").is_ok() {
        let mut pr: Vec<_> = prof.iter().filter(|&&(r, _)| r < 3.0).cloned().collect();
        pr.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        pr.dedup_by(|a, b| (a.0 - b.0).abs() < 1e-9);
        for (r, u) in pr.iter().step_by(3) {
            println!("r={r:.3} u={u:+.5} exact={:+.5}", if *r > ri && *r < ro { exact(*r) } else { 0.0 });
        }
    }
    let mhd = s.mhd.as_ref().unwrap();
    let joule = mhd.joule_dissipation(&g, &fft, &s.u);
    println!(
        "N={n} I={current} walls={walls} Re={re:.0}: rel L2 error {:.3e}, max {:.3e} (u_max {umax:.4}), t={:.0} ({} steps, {:.0}s); Joule dissipation {joule:.3e}",
        (e2 / n2).sqrt(),
        emax / umax,
        s.t,
        (s.t / dt).round(),
        t0.elapsed().as_secs_f64()
    );
}
