//! Gate-level power against the supply current of every transistor.

use chip::flat::run_cycles;
use chip::library::{self, Corner};
use chip::power::{analyse, random_vectors, toggle_counts};
use chip::rtl::parse;
use chip::sta::Constraints;
use chip::synth::synthesize;
use semiconductor::gates::{DeviceModel, Ekv};

#[test]
fn switching_energy_against_transistors() {
    let m = parse(
        "module pipe(input clk, input [1:0] a, input [1:0] b, output reg [2:0] s);
           reg [1:0] ra, rb;
           always @(posedge clk) begin ra <= a; rb <= b; s <= ra + rb; end
         endmodule",
    )
    .unwrap();
    let n = synthesize(&m).unwrap();
    let lib = library::at(Corner::NOMINAL).unwrap();
    let (en, ep) = (Ekv::n(), Ekv::p());
    let model = DeviceModel::new(&en, &ep);
    let cycles = 10;
    let vectors = random_vectors(&n, 0.5, cycles, 7);
    let period = 1e-9;
    let tr = run_cycles(&model, &n, &vectors, period, 100e-12, 80e-12, 2e-12);
    let c = Constraints { period, wire_cap: 0.0, output_load: 0.0, ..Default::default() };
    let counts = toggle_counts(&n, &vectors);
    let per_cycle: Vec<f64> = counts.iter().map(|t| t / cycles as f64).collect();
    let p = analyse(&n, &lib, &c, &per_cycle);
    // Gate-level energy over the same window: switching + clock + leakage.
    let gate = p.energy_per_cycle * cycles as f64 + p.leakage * tr.duration;
    let toggles: f64 = counts.iter().sum();
    println!(
        "{cycles} cycles, {toggles} net toggles: transistor-level {:.1} fJ, gate-level {:.1} fJ ({:+.1} %)",
        tr.energy * 1e15,
        gate * 1e15,
        100.0 * (gate / tr.energy - 1.0)
    );
    // The functional check rides along: outputs sampled each cycle.
    assert_eq!(tr.outputs.len(), cycles);
    let err = (gate / tr.energy - 1.0).abs();
    assert!(err < 0.25, "gate-level energy off by {:.1} %", 100.0 * err);
}
