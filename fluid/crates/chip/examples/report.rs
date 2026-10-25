//! Synthesize a Verilog file and report timing and power at the nominal
//! corner: `cargo run -p chip --release --example report -- design.v [MHz]`.
use chip::library::{self, Corner};
use chip::{equiv, power, rtl, sta, synth};
fn main() {
    let path = std::env::args().nth(1).expect("Verilog file");
    let mhz: f64 = std::env::args().nth(2).map(|s| s.parse().unwrap()).unwrap_or(500.0);
    let m = rtl::parse(&std::fs::read_to_string(&path).unwrap()).unwrap_or_else(|e| panic!("{e}"));
    let n = synth::synthesize(&m).unwrap_or_else(|e| panic!("{e}"));
    let eq = equiv::check(&m, &n, 20_000).unwrap();
    println!("{}: {:?}, {} transistors; equivalence on {} points{}: {}", n.name, n.cells(), n.transistors(), eq.points, if eq.exhaustive { " (all)" } else { "" }, if eq.mismatch.is_none() { "ok" } else { "MISMATCH" });
    let t0 = std::time::Instant::now();
    let lib = library::at(Corner::NOMINAL).unwrap();
    println!("library ready in {:.1} s", t0.elapsed().as_secs_f64());
    let c = sta::Constraints { period: 1e-6 / mhz, input_delay: 100e-12, wire_cap: 1e-15, ..Default::default() };
    let t = sta::analyse(&n, &lib, &c);
    println!("fmax {:.0} MHz (period {:.1} ± {:.1} ps), hold slack {:.1} ps, endpoint {}", 1e-6 / t.min_period, t.min_period * 1e12, t.min_period_error * 1e12, t.hold_slack * 1e12, t.endpoint);
    for s in &t.critical {
        println!("   {:>10} {:>4} {:7.1} ps", s.cell, s.edge, s.arrival * 1e12);
    }
    let a = power::activity(&n, 0.25, 4000, 1);
    let p = power::analyse(&n, &lib, &c, &a);
    println!("at {mhz} MHz: {:.3} mW switching+clock ({:.2} pJ/cycle), {:.3} µW leakage", (p.dynamic + p.clock) * 1e3, p.energy_per_cycle * 1e12, p.leakage * 1e6);
}
