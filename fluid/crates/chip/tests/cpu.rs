//! A small processor through the whole flow: it must run its program, in
//! the RTL and in the gate netlist, cycle for cycle.

use std::collections::BTreeMap;

use chip::equiv::{check, netlist_step, Interp};
use chip::rtl::parse;
use chip::synth::synthesize;

#[test]
fn mini8_runs_fibonacci_in_rtl_and_gates() {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../projects/mini-cpu/mini8.v");
    let m = parse(&std::fs::read_to_string(path).unwrap()).unwrap_or_else(|e| panic!("{e}"));
    let n = synthesize(&m).unwrap_or_else(|e| panic!("{e}"));
    println!("mini8: {:?}, {} transistors, {} flip-flops", n.cells(), n.transistors(), n.flops.len());
    let mut state: BTreeMap<String, u128> = m.regs.iter().map(|r| (r.clone(), 0)).collect();
    let mut gstate = state.clone();
    let mut seen = Vec::new();
    for cycle in 0..120 {
        let rst = (cycle == 0) as u128;
        let mut vals = state.clone();
        vals.insert("rst".into(), rst);
        let (outs, next) = Interp::new(&m, vals.clone()).step().unwrap();
        let mut gvals = gstate.clone();
        gvals.insert("rst".into(), rst);
        let (gouts, gnext) = netlist_step(&n, &m, &gvals);
        assert_eq!(outs, gouts, "cycle {cycle}: outputs differ");
        assert_eq!(next, gnext, "cycle {cycle}: next state differs");
        if seen.last() != Some(&outs["out"]) {
            seen.push(outs["out"]);
        }
        state = next;
        gstate = gnext;
    }
    // 0 after reset, then the Fibonacci numbers mod 256.
    let want: Vec<u128> = vec![0, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89];
    assert_eq!(&seen[..want.len()], &want[..], "port sequence {seen:?}");
    // Equivalence over the whole (input, state) space is too large to
    // enumerate (54 state bits): it is sampled.
    let eq = check(&m, &n, 20_000).unwrap();
    assert!(eq.mismatch.is_none() && !eq.exhaustive);
}
