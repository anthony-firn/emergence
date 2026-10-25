//! DFG benchmark 2D-1 (Schäfer & Turek 1996; reference values from the
//! FeatFlow benchmark page): steady channel flow past a cylinder, Re = 20.
//! The cylinder enters as a surface file (STL written by the geometry
//! crate, or a STEP file tessellated by gmsh with `--step`).
//!
//! `cargo run --release -p immersed --example schafer_turek_2d -- 82 [--step] [--vtk out_dir]`

use std::time::Instant;

use geometry::primitives::cylinder_z;
use geometry::{Sdf, TriMesh};
use immersed::{Bc, Body, Grid, Solver};

const CD_REF: f64 = 5.57953523384;
const CL_REF: f64 = 0.010618948146;
const DP_REF: f64 = 0.11752016697;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let ny: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(82);
    let step = args.iter().any(|a| a == "--step");
    let vtk_dir = args.iter().position(|a| a == "--vtk").and_then(|i| args.get(i + 1)).cloned();
    let (hgt, len, d, umax, nu) = (0.41, 2.2, 0.1, 0.3, 1e-3);
    let h = hgt / ny as f64;
    let nx = (len / h).round() as usize;
    let dir = std::env::temp_dir();
    let mesh = if step {
        let path = dir.join("st_cylinder.step");
        let py = geometry::cad::gmsh_python().expect("gmsh: run fluid/scripts/setup-cad.sh");
        let make = format!(
            "import gmsh\ngmsh.initialize()\ngmsh.option.setNumber('General.Terminal',0)\ngmsh.model.occ.addCylinder(0.2,0.2,-0.5,0,0,1.5,0.05)\ngmsh.model.occ.synchronize()\ngmsh.write('{}')\ngmsh.finalize()\n",
            path.display()
        );
        assert!(std::process::Command::new(py).args(["-c", &make]).status().unwrap().success());
        let opt = geometry::cad::Tessellation { elements_per_2pi: 2048, max_size: 0.05 };
        geometry::cad::tessellate(&path, &opt).unwrap()
    } else {
        let path = dir.join("st_cylinder.stl");
        cylinder_z([0.2, 0.2, -0.5], 0.5 * d, 1.5, 4096).write_stl(&path).unwrap();
        TriMesh::read(&path).unwrap()
    };
    let r = mesh.report();
    println!("surface: {} triangles, watertight {}, {}", r.triangles, r.watertight, if step { "STEP via gmsh" } else { "STL" });
    let inlet = move |x: [f64; 3], _t: f64| [4.0 * umax * x[1] * (hgt - x[1]) / (hgt * hgt), 0.0, 0.0];
    let bcs = [[Bc::inlet(inlet), Bc::Outlet], [Bc::wall(), Bc::wall()], [Bc::Periodic, Bc::Periodic]];
    let grid = Grid::new([nx, ny, 1], h, [0.0; 3]);
    let t0 = Instant::now();
    let mut s = Solver::new(grid, nu, 1.0, bcs, vec![Body::fixed(Sdf::new(mesh))]);
    s.initialize(|x| inlet(x, 0.0));
    let setup = t0.elapsed();
    let umean = 2.0 / 3.0 * umax;
    let coef = |f: f64| 2.0 * f / (umean * umean * d * h); // per unit depth (depth = h)
    let mut last = 0.0;
    for block in 1..=200 {
        let (steps, rate) = s.run_until_steady(0.8, 1.0 * block as f64, 1e-7);
        let fm = s.body_force_momentum(0);
        let cd = coef(fm[0]);
        let change = (cd - last).abs();
        last = cd;
        if block % 5 == 0 || rate < 1e-7 {
            println!("t {:5.1}: Cd(momentum) {cd:.5}, rate {rate:.2e}, ΔCd/1t {change:.1e} ({steps} steps)", s.t);
        }
        if rate < 1e-7 {
            break;
        }
    }
    let solve = t0.elapsed() - setup;
    let fm = s.body_force_momentum(0);
    let fs = s.surface_force(0);
    let ft = fs.total();
    let dp = s.wall_pressure(0, [0.15, 0.2, 0.5 * h]) - s.wall_pressure(0, [0.25, 0.2, 0.5 * h]);
    let rel = |x: f64, r: f64| 100.0 * (x - r) / r;
    println!("ny {ny} (D/h = {:.1}), {nx}×{ny} cells, setup {setup:.2?}, solve {solve:.2?}", d / h);
    println!("  Cd momentum {:.5} ({:+.2}%), surface {:.5} ({:+.2}%) [pressure {:.4} + viscous {:.4}]", coef(fm[0]), rel(coef(fm[0]), CD_REF), coef(ft[0]), rel(coef(ft[0]), CD_REF), coef(fs.pressure[0]), coef(fs.viscous[0]));
    println!("  Cl momentum {:.5} ({:+.1}%), surface {:.5} ({:+.1}%)", coef(fm[1]), rel(coef(fm[1]), CL_REF), coef(ft[1]), rel(coef(ft[1]), CL_REF));
    println!("  Δp {dp:.5} ({:+.2}%); wetted length {:.5} (πD = {:.5})", rel(dp, DP_REF), fs.area / h, std::f64::consts::PI * d);
    if let Some(dir) = vtk_dir {
        std::fs::create_dir_all(&dir).unwrap();
        s.write_vtk(format!("{dir}/st2d_fields_{ny}.vtk")).unwrap();
        s.write_surface_vtk(0, format!("{dir}/st2d_cylinder_{ny}.vtk"), 0.0, umean).unwrap();
        println!("  wrote VTK to {dir}");
    }
}
