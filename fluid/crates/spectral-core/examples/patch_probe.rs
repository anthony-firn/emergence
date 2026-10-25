//! Stage 1 of overset refinement: can a patch reproduce an exact flow from
//! boundary data alone?
//!
//! Physical flow: steady forced ABC (Beltrami) field, forcing f = νk²u, an
//! exact Navier–Stokes solution. Patch: a periodic box of side Lp centred
//! at xc, coordinates x = xc + (ξ − π)·s with s = Lp/(2π). In ξ-space the
//! equations are standard Navier–Stokes with ν_p = ν/s, time τ = t/s and
//! forcing s·f. The patch's fluid region is the sphere |ξ − π| < R; its
//! boundary velocity comes from the exact field at the markers, imposed by
//! the smooth forcing extension in the buffer outside the sphere. Forcing
//! is windowed to zero before the patch edges (the buffer is fictitious).
//!
//! usage: patch_probe N_patch [ext_width]
use spectral_core::nufft::Nufft;
use spectral_core::surface::{resolution, Extension};
use spectral_core::{fft::Fft3, init, Grid, Solver, Surface};

const PI: f64 = std::f64::consts::PI;

fn abc(x: [f64; 3]) -> [f64; 3] {
    let (a, b, c) = (1.0, 0.7, 0.4);
    [a * x[2].sin() + c * x[1].cos(), b * x[0].sin() + a * x[2].cos(), c * x[1].sin() + b * x[0].cos()]
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let n: usize = args[1].parse().unwrap();
    let wext: f64 = args.get(2).map(|v| v.parse().unwrap()).unwrap_or(0.9);
    let (nu, lp, xc) = (1.0, 2.5, [1.0, 2.0, 3.0]);
    let s = lp / (2.0 * PI);
    let to_phys = |xi: [f64; 3]| [xc[0] + (xi[0] - PI) * s, xc[1] + (xi[1] - PI) * s, xc[2] + (xi[2] - PI) * s];
    let r_fluid = 2.0;
    let g = Grid::new(n);

    // windowed forcing s·νk²·u (k = 1): 1 inside r < R + 0.6, C^∞ to 0 by R + 1.0
    let window = |r: f64| {
        let (r0, r1) = (r_fluid + 0.6, r_fluid + 1.0);
        if r <= r0 {
            1.0
        } else if r >= r1 {
            0.0
        } else {
            // C^∞ step from 1 to 0 on (r0, r1): any kink in the total forcing,
            // even in the fictitious buffer, would cap the convergence order
            let t = (r - r0) / (r1 - r0);
            let e = |y: f64| if y > 0.0 { (-1.0 / y).exp() } else { 0.0 };
            e(1.0 - t) / (e(1.0 - t) + e(t))
        }
    };
    let f_phys = init::from_fn(&g, |a, b, c| {
        let xi = [a, b, c];
        let r = ((a - PI).powi(2) + (b - PI).powi(2) + (c - PI).powi(2)).sqrt();
        let u = abc(to_phys(xi));
        let wgt = s * nu * window(r);
        [wgt * u[0], wgt * u[1], wgt * u[2]]
    });
    let fft = Fft3::new(&g);
    let mut fh: [Vec<num_complex::Complex64>; 3] = std::array::from_fn(|_| vec![Default::default(); g.spec_len()]);
    for c in 0..3 {
        fft.forward(&f_phys[c], &mut fh[c]);
    }

    let mut surf = Surface::sphere_inverse([PI; 3], r_fluid, 0.5 * resolution(&g));
    surf.u_markers = Some(surf.points.iter().map(|&xm| abc(to_phys(xm))).collect());
    surf.extension = Some(Extension { width: wext, sigma_rel: 1.0 });
    let surf = surf.with_basis(resolution(&g));
    let mut sol = Solver::new(g.clone(), nu / s, init::from_fn(&g, |_, _, _| [0.0; 3]));
    sol.forcing = Some(fh);
    sol.surface = Some(surf.clone());
    let t0 = std::time::Instant::now();
    // explicit nonlinear term: respect CFL (|u| ~ 1.5); run ~5 viscous times
    let t_end = 5.0 * r_fluid * r_fluid / (nu / s);
    while sol.t < t_end {
        let dt = sol.cfl_dt(0.5).min(0.05).min(t_end - sol.t).max(1e-9);
        sol.step(dt);
    }

    // compare with the exact field at interior points
    let mut seed = 0x2545_F491_4F6C_DD1Du64;
    let mut rnd = || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        (seed >> 11) as f64 / (1u64 << 53) as f64
    };
    let pts: Vec<[f64; 3]> = (0..400)
        .map(|_| loop {
            let p = [PI + (2.0 * rnd() - 1.0) * r_fluid, PI + (2.0 * rnd() - 1.0) * r_fluid, PI + (2.0 * rnd() - 1.0) * r_fluid];
            let r = ((p[0] - PI).powi(2) + (p[1] - PI).powi(2) + (p[2] - PI).powi(2)).sqrt();
            if r < r_fluid {
                break p;
            }
        })
        .collect();
    let mut nf = Nufft::new(&g, 12);
    let vals: Vec<Vec<f64>> = (0..3).map(|c| nf.interpolate(&sol.u[c], &pts)).collect();
    let mut err: f64 = 0.0;
    let mut scale: f64 = 0.0;
    for (q, &p) in pts.iter().enumerate() {
        let ue = abc(to_phys(p));
        for c in 0..3 {
            err = err.max((vals[c][q] - ue[c]).abs());
            scale = scale.max(ue[c].abs());
        }
    }
    println!(
        "patch N={n:3} (resolution {:.4} physical) ext w={wext}: max |u − u_exact| / max|u| inside = {:.2e}  krylov {:?} resid {:?}  {:.1}s",
        resolution(&g) * s,
        err / scale,
        sol.last_cg_iters().unwrap(),
        sol.last_cg_resid().unwrap().map(|r| format!("{r:.0e}")),
        t0.elapsed().as_secs_f64()
    );
}
