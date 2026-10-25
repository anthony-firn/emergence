//! Build the tabulated real-fluid EOS and measure its interpolation error.
//!   realfluid_table water|nitrogen N
use compressible::realfluid::RealFluidEos;
use compressible::reference_data::{nitrogen_span2000, water_iapws95};
use std::time::Instant;
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let fl = a.get(1).cloned().unwrap_or("water".into());
    let n: usize = a.get(2).map(|v| v.parse().unwrap()).unwrap_or(200);
    let t0 = Instant::now();
    let rf = if fl == "water" { RealFluidEos::build(water_iapws95(), 0.01, 1100.0, 275.0, 900.0, n) } else { RealFluidEos::build(nitrogen_span2000(), 0.05, 900.0, 64.0, 400.0, n) };
    println!("{fl}: two {n}×{n} tables built in {:.1}s", t0.elapsed().as_secs_f64());
    let rep = rf.error_report(4000);
    for (h, half) in ["vapour side (log ρ)", "liquid side (linear ρ)"].iter().enumerate() {
        for (k, name) in ["single-phase", "two-phase"].iter().enumerate() {
            println!("  {half:22} {name:12}: {:4} samples, max relative error p {:.2e}, T {:.2e}, c {:.2e}", rep[h][k].0, rep[h][k].1[0], rep[h][k].1[1], rep[h][k].1[2]);
        }
    }
    let t0 = Instant::now();
    let mut acc = 0.0;
    let m = 1_000_000;
    for i in 0..m {
        let rho = 900.0 + 90.0 * (i as f64 / m as f64);
        acc += rf.eval(rho, 3e5).p;
    }
    println!("  table lookup (liquid): {:.0} ns per evaluation (checksum {acc:.3e})", t0.elapsed().as_secs_f64() * 1e9 / m as f64);
}
