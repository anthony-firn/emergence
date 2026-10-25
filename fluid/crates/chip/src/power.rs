//! Power from switching activity and the characterised cell energies.
//!
//! Activity comes from simulating the netlist cycle by cycle: random input
//! vectors in which each input bit toggles with probability `activity` per
//! cycle, from the all-zero register state, counting how often every net's
//! settled value changes. Each output transition of a cell costs the supply
//! energy characterised at transistor level for that cell at the net's
//! slew and load (it includes the short-circuit current that ½CV² misses).
//! Every flip-flop pays its clock energy every cycle plus its clock-to-Q
//! energy when Q changes. Leakage is each cell's worst-state leakage.
//!
//! Not included, and so the stated error: glitches (settled values only),
//! internal nodes that switch without the output switching, and the input
//! pins' load, which the external drivers charge.

use crate::library::CellLib;
use crate::sta::{self, Constraints};
use crate::synth::Netlist;

#[derive(Clone, Copy, Debug)]
pub struct Power {
    /// Supply energy per clock cycle (J).
    pub energy_per_cycle: f64,
    pub dynamic: f64,
    pub clock: f64,
    pub leakage: f64,
    pub total: f64,
}

/// Random input vectors (module input bits, in `Netlist::inputs` order) in
/// which each bit toggles with probability `input_activity` per cycle.
pub fn random_vectors(n: &Netlist, input_activity: f64, cycles: usize, seed: u64) -> Vec<Vec<bool>> {
    let bits: usize = n.inputs.iter().map(|(_, b)| b.len()).sum();
    let mut rng = seed | 1;
    let mut rand = || {
        rng ^= rng << 13;
        rng ^= rng >> 7;
        rng ^= rng << 17;
        (rng >> 11) as f64 / (1u64 << 53) as f64
    };
    let mut cur = vec![false; bits];
    (0..cycles)
        .map(|_| {
            for b in cur.iter_mut() {
                if rand() < input_activity {
                    *b = !*b;
                }
            }
            cur.clone()
        })
        .collect()
}

/// Total toggles of every net over the cycles of `vectors` (zero-delay:
/// settled values), from the all-zero register state with all-zero inputs.
pub fn toggle_counts(n: &Netlist, vectors: &[Vec<bool>]) -> Vec<f64> {
    let pos: std::collections::BTreeMap<usize, usize> = n.comb.inputs.iter().enumerate().map(|(i, &x)| (x, i)).collect();
    let in_nets: Vec<usize> = n.inputs.iter().flat_map(|(_, b)| b.iter().copied()).collect();
    let mut iv = vec![false; n.comb.inputs.len()];
    iv[pos[&n.one]] = true;
    let mut prev = n.comb.eval(&iv);
    let mut toggles = vec![0.0; n.comb.n_nets];
    for v in vectors {
        // Clock edge: registers take their D values; inputs move.
        for f in &n.flops {
            iv[pos[&f.q]] = prev[f.d];
        }
        for (k, &net) in in_nets.iter().enumerate() {
            iv[pos[&net]] = v[k];
        }
        let now = n.comb.eval(&iv);
        for (t, (a, b)) in toggles.iter_mut().zip(prev.iter().zip(&now)) {
            *t += (a != b) as u8 as f64;
        }
        prev = now;
    }
    toggles
}

/// Toggles per cycle of every net, from `cycles` random cycles.
pub fn activity(n: &Netlist, input_activity: f64, cycles: usize, seed: u64) -> Vec<f64> {
    let v = random_vectors(n, input_activity, cycles, seed);
    toggle_counts(n, &v).into_iter().map(|t| t / cycles as f64).collect()
}

pub fn analyse(n: &Netlist, lib: &CellLib, c: &Constraints, toggles: &[f64]) -> Power {
    let load = sta::loads(n, lib, c);
    let slew = c.input_slew;
    let mut e_cycle = 0.0;
    for inst in &n.comb.insts {
        let ct = lib.comb.cell(inst.kind);
        for (o, &out) in inst.outs.iter().enumerate() {
            // Mean over the pins of the rise and fall energies: half of the
            // toggles rise, half fall.
            let pins = &ct.outs[o];
            let e: f64 = pins
                .iter()
                .map(|p| 0.5 * (p.e_rise.interp(slew, load[out]).0 + p.e_fall.interp(slew, load[out]).0))
                .sum::<f64>()
                / pins.len() as f64;
            e_cycle += toggles[out] * e;
        }
    }
    let ff = &lib.dff;
    let mut clock = 0.0;
    for f in &n.flops {
        clock += ff.e_clock;
        let e = 0.5 * (ff.e_rise.interp(c.clock_slew, load[f.q]).0 + ff.e_fall.interp(c.clock_slew, load[f.q]).0);
        e_cycle += toggles[f.q] * e;
    }
    let freq = 1.0 / c.period;
    let leak_i: f64 = n.comb.insts.iter().map(|i| lib.comb.cell(i.kind).i_leak).sum::<f64>() + ff.i_leak * n.flops.len() as f64;
    let leakage = leak_i * lib.corner.vdd;
    let dynamic = e_cycle * freq;
    let clock_p = clock * freq;
    Power {
        energy_per_cycle: e_cycle + clock,
        dynamic,
        clock: clock_p,
        leakage,
        total: dynamic + clock_p + leakage,
    }
}
