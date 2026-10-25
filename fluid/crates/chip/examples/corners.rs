//! Check the supply × temperature interpolation of `chip::corners` against
//! libraries characterised directly at held-out corners:
//! `cargo run -p chip --release --example corners -- design.v`.
use chip::corners::{at, exact};
use chip::library::Corner;
use chip::{power, rtl, sta, synth};
fn main() {
    let path = std::env::args().nth(1).expect("Verilog file");
    let m = rtl::parse(&std::fs::read_to_string(&path).unwrap()).unwrap();
    let n = synth::synthesize(&m).unwrap();
    let c = sta::Constraints { period: 4e-9, input_delay: 100e-12, wire_cap: 1e-15, ..Default::default() };
    let a = power::activity(&n, 0.25, 4000, 1);
    let t0 = std::time::Instant::now();
    for (vdd, temp) in [(1.14, 293.75), (1.26, 381.25), (1.2, 300.0), (1.10, 400.0)] {
        let (mi, ind) = at(&n, &c, &a, vdd, temp).unwrap();
        let (me, _) = exact(&n, &c, &a, Corner { vdd, temp }).unwrap();
        println!(
            "{vdd:.2} V {temp:6.2} K: period interp {:.2} ps exact {:.2} ps (error {:+.2} ps, indicator {:.2}); energy/cycle error {:+.2} % (indicator {:.2} %); leakage error {:+.2} % (indicator {:.2} %)",
            mi.min_period * 1e12,
            me.min_period * 1e12,
            (mi.min_period - me.min_period) * 1e12,
            ind.min_period * 1e12,
            100.0 * (mi.energy_per_cycle / me.energy_per_cycle - 1.0),
            100.0 * ind.energy_per_cycle / me.energy_per_cycle,
            100.0 * (mi.leakage / me.leakage - 1.0),
            100.0 * ind.leakage
        );
    }
    println!("({:.0} s including characterisation)", t0.elapsed().as_secs_f64());
}
