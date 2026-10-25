//! Does the local truncation-residual density predict *where* an
//! under-resolved run is wrong? Taylor–Green Re = 400 at 32³ vs a 96³
//! reference: correlate indicators with the local velocity error, and check
//! the continuum regime for water and nitrogen gas at lab scales.
use spectral_core::regime::{Material, RegimeMap};
use spectral_core::{init, Grid, Solver};

fn pearson(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len() as f64;
    let (ma, mb) = (a.iter().sum::<f64>() / n, b.iter().sum::<f64>() / n);
    let (mut sab, mut saa, mut sbb) = (0.0, 0.0, 0.0);
    for (x, y) in a.iter().zip(b) {
        sab += (x - ma) * (y - mb);
        saa += (x - ma).powi(2);
        sbb += (y - mb).powi(2);
    }
    sab / (saa * sbb).sqrt()
}

fn rank(v: &[f64]) -> Vec<f64> {
    let mut idx: Vec<usize> = (0..v.len()).collect();
    idx.sort_by(|&i, &j| v[i].partial_cmp(&v[j]).unwrap());
    let mut r = vec![0.0; v.len()];
    for (k, &i) in idx.iter().enumerate() {
        r[i] = k as f64;
    }
    r
}

fn main() {
    let (re, dt) = (400.0, 0.01);
    let t_end: f64 = std::env::args().nth(1).map(|v| v.parse().unwrap()).unwrap_or(3.0);
    let (gc, gf) = (Grid::new(32), Grid::new(96));
    let mut c = Solver::new(gc.clone(), 1.0 / re, init::taylor_green(&gc));
    let mut f = Solver::new(gf.clone(), 1.0 / re, init::taylor_green(&gf));
    let steps = (t_end / dt).round() as usize;
    for _ in 0..steps {
        c.step(dt);
        f.step(dt);
    }
    let uc: Vec<Vec<f64>> = (0..3).map(|k| c.velocity_phys(k)).collect();
    let uf: Vec<Vec<f64>> = (0..3).map(|k| f.velocity_phys(k)).collect();
    let n = 32;
    let err: Vec<f64> = (0..n * n * n)
        .map(|p| {
            let (i, j, l) = (p / (n * n), (p / n) % n, p % n);
            let q = (3 * i * 96 + 3 * j) * 96 + 3 * l;
            (0..3).map(|k| (uc[k][p] - uf[k][q]).powi(2)).sum::<f64>().sqrt()
        })
        .collect();
    let map = c.regime_map();
    let inv_k: Vec<f64> = map.kolmogorov.iter().map(|x| 1.0 / x).collect();
    println!("t = {t_end}: correlation of indicator with the local velocity error |u_32 - u_96|");
    for (name, ind) in [("truncation density", &map.truncation), ("1 / (k_max eta)", &inv_k), ("|vorticity|", &map.vorticity)] {
        println!("  {name:20}: Pearson {:+.3}  Spearman {:+.3}", pearson(ind, &err), pearson(&rank(ind), &rank(&err)));
    }
    let thr = 0.5 * map.truncation.iter().cloned().fold(0.0, f64::max);
    let regs = RegimeMap::regions(&map.truncation, n, thr);
    println!("  regions with truncation density > 50% of max: {}", regs.len());
    for r in regs.iter().take(4) {
        println!("    {} points, centroid ({:.2}, {:.2}, {:.2}), max {:.3e}", r.points, r.centroid[0], r.centroid[1], r.centroid[2], r.max);
    }
    let kmin = map.kolmogorov.iter().cloned().fold(f64::INFINITY, f64::min);
    println!("  min k_max*eta = {kmin:.2} (DNS practice wants >= 1.5)");
    // continuum regime: the same flow as water / nitrogen in a 1 cm box at 1 m/s
    for (mat, len, vel) in [(Material::WATER, 0.01, 1.0), (Material::NITROGEN_GAS, 0.01, 1.0), (Material::NITROGEN_GAS, 1e-6, 100.0)] {
        let cc = c.continuum_check(&mat, len, vel);
        println!(
            "  {} in a {len:.0e} m box at {vel} m/s: lambda|grad u|/c {:.1e}, Mach {:.1e}, lambda/dx {:.1e} -> continuum {}",
            mat.name, cc.non_equilibrium, cc.mach, cc.mfp_over_dx, if cc.valid() { "valid" } else { "NOT valid" }
        );
    }
}
