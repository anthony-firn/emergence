//! Time a few steps of uniform flow past a sphere with surface walls and
//! report the wall-solve iteration counts (diagnostic).
//! usage: sphere_flow_gpu N point|ext steps
use spectral_core::surface::{resolution, Extension};
use spectral_core::{init, Backend, Grid, Penalization, Surface};
use spectral_cuda::CudaSolver;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let n: usize = a[1].parse().unwrap();
    let ext = a[2] == "ext";
    let steps: usize = a[3].parse().unwrap();
    let two_pi = 2.0 * std::f64::consts::PI;
    let (d, re) = (two_pi / 8.0, 100.0);
    let nu = d / re;
    let g = Grid::new(n);
    let center = [0.3 * two_pi, 0.5 * two_pi, 0.5 * two_pi];
    let mut sf = Surface::sphere(center, 0.5 * d, if ext { 0.5 } else { 1.0 } * resolution(&g));
    if ext {
        sf.extension = Some(Extension { width: 0.3 * d, sigma_rel: 1.0 });
        sf = sf.with_basis(resolution(&g));
    }
    let pen = Penalization::sphere(&g, center, 0.0, g.dx(), 1.0).with_sponge(&g, 0.8 * two_pi, 0.15 * two_pi, 0.05, [1.0, 0.0, 0.0]);
    let mut s = CudaSolver::new(g.clone(), nu, init::from_fn(&g, |_, _, _| [0.0; 3])).unwrap();
    s.set_penalization(&pen).unwrap();
    s.set_surface(&sf).unwrap();
    s.set_mean_flow([1.0, 0.0, 0.0]).unwrap();
    for k in 0..steps {
        let t0 = std::time::Instant::now();
        let dt = s.cfl_dt(0.5);
        s.try_step(dt).unwrap();
        let (it, res) = s.last_cg().unwrap();
        let f = s.try_solid_force().unwrap().unwrap();
        println!("step {k}: dt={dt:.4} iters {it:?} resid {res:?} Fx={:.4e} {:.2}s", f[0], t0.elapsed().as_secs_f64());
    }
}
