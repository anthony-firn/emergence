# Gate-level logic with certified timing: results

Crate: `fluid/crates/semiconductor`. Every number below comes from
`examples/adder.rs` (`cargo run -p semiconductor --example adder`) and from
the tests in `tests/gates.rs`.

This is the rung above the circuit level in `analysis/emergence_ladder.md`:
standard cells built from the compact MOSFET law, characterised into
non-linear delay model (NLDM) timing/energy tables, a gate-level
event-driven logic simulator on those tables, static timing analysis, and
the fold back to transistor level with an error certificate.

The transistor rung here is the same all-region EKV law that
`compact.rs` uses as the held-out reference for its PCHIP table
(Enz–Krummenacher–Vittoz). It is calibrated to the scaled 150 nm
technology of `analysis/semiconductor.md` §6–7: `V_dd = 1.2 V`,
`I_on,n = 330 µA/µm`, `I_on,p = 230 µA/µm`, gate capacitance
`1.031 fF/µm`, drain diffusion `1.4 fF/µm`, source `0.7 fF/µm`,
`W_p/W_n = 1.5`. The same lumped-capacitance model as the ring
oscillator (§7) is used, and is the largest approximation of this rung.

## 1. Standard-cell library

Cells are transistor templates (`CellKind::template`):

| cell | transistors | topology |
|---|---|---|
| INV | 2 | 1 nMOS, 1 pMOS |
| NAND2 | 4 | 2 series nMOS, 2 parallel pMOS |
| NOR2 | 4 | 2 parallel nMOS, 2 series pMOS |
| XOR2 | 16 | four NAND2 (classic 4-NAND XOR) |
| FA | 44 | two XOR2 + three NAND2, sum and carry |

Characterised at input slew 100 ps and load 4 fF
(`E_sw` is the energy delivered by `V_dd` over one switching transition;
`P_leak` is the worst-case leakage power over all input states):

| cell | c_in (fF) | c_out (fF) | d_rise (ps) | d_fall (ps) | E_sw (fJ) | P_leak (nW) |
|---|---|---|---|---|---|---|
| INV | 2.577 | 3.500 | 39.89 | 40.96 | 11.29 | 0.175 |
| NAND2 | 2.577 | 5.600 | 45.20 | 58.48 | 14.25 | 0.351 |
| NOR2 | 2.577 | 4.900 | 54.12 | 44.56 | 13.26 | 0.175 |
| XOR2 | 5.155 | 5.600 | 89.94 | 83.76 | 15.29 | 0.701 |
| FA | 7.732 | 5.600 | 176.20 | 178.89 | 37.39 | 1.402 |

- `c_in` is the sum of the gate capacitances of the transistors on an
  input; the pMOS contributes `W_p/W_n = 1.5` times the nMOS. All inputs
  of a cell have the same `c_in`.
- Leakage is a **switch-level** estimate, not a direct DC Newton. The
  floating mid-node of a series stack has no DC path to a rail, so a direct
  Newton on the resulting near-singular system ran to ±30 V. Instead each
  transistor is first treated as a switch (on 10 mS, off 1e-14 S), the
  linear network is solved, switches are re-classified from the solved
  gate–source voltages to a fixed point, and the non-linear EKV current is
  evaluated at the resulting node voltages. The values scale with the
  number of stacked off devices, as expected.

## 2. NLDM characterisation and the certificate

Each (cell, input pin, output, edge) is characterised by a transistor-level
transient:

- the input is a linear ramp of the requested 10–90 % slew through a
  settled state;
- delay is the input-50 % to output-50 % time;
- output slew is the output 10–90 % (rise) or 90–10 % (fall) time;
- switching energy is `V_dd ∫ I_supply dt` over the transition, so it
  includes the short-circuit term the `½CV²` estimate misses.

The result is a bilinear NLDM table on

`SLEWS = {15, 40, 120, 350} ps`, `LOADS = {1, 4, 12} fF`.

The slew grid deliberately reaches fast edges: an inverter driven by
another inverter settles to an output slew of ~26–43 ps, so a grid starting
at 30 ps would clamp most of a chain. Extending it below 30 ps cut the
inverter-chain fold error from +9.5 % to +4.4 % and brought every chain
stage inside the characterised box.

The certificate is the maximum relative delay error of the interpolated
table against **held-out** transistor-level points not on the grid
(`(55 ps, 2 fF), (55, 8), (180, 2), (180, 8)`), together with the absolute
output-slew error:

| cell | worst delay certificate |
|---|---|
| INV | 5.45 % |
| NAND2 | 5.22 % |
| NOR2 | 5.27 % |
| XOR2 | 2.82 % |
| FA | 2.85 % |

An independent dense scan of the INV box (every 10 ps × 0.5 fF) found a
worst interpolation error of 5.97 %, so the four held-out points are a
fair (slightly conservative) sample of the sup for that cell.

## 3. Gate-level event simulator

`timing::simulate` is an event-driven logic simulator. Each transition
carries (arrival time, new logic value, slew); a changing input schedules
the affected gate outputs, whose delay and output slew are bilinearly
interpolated from the NLDM table. The load of a net is the sum of the
fan-out input capacitances plus a 0.2 fF wire estimate; the driver's own
drain-diffusion capacitance is already inside the characterisation and is
not added. Energy is accumulated one `e_rise`/`e_fall` table entry per
output transition, and leakage energy is `i_leak · V_dd · Δt`.

`timing::sta` propagates late arrivals and slews through the netlist and
returns the critical path. `timing::critical_path_bound` accumulates the
per-cell certificates along the worst path, choosing the later of the rise
and fall tables at each gate (so the bound is conservative over the edge
direction).

## 4. The fold: gate level against transistor level

All transistor-level references flatten the gate netlist into one
transistor network and integrate it with backward Euler, a scattered
analytic Jacobian and damped Newton (`gates::simulate_flat`). The supply
energy is `V_dd ∫ I_supply dt` on the same waveform.

**Inverter chain** (8 stages, terminated by an inverter load so the whole
measured path is inside the characterised box), 100 ps input edge:

| edge | delay gate | delay transistor | error | interpolation bound | energy error |
|---|---|---|---|---|---|
| low→high | 176.7 ps | 169.3 ps | +4.4 % | 10.5 ps (within) | +4.0 % |
| high→low | 176.1 ps | 168.7 ps | +4.4 % | 10.5 ps (within) | +7.8 % |

**4-bit ripple-carry adder**, 100 ps input edge (carry-out and sum MSB):

| transition | carry gate | carry transistor | error | sumMSB gate | sumMSB transistor | error | energy error |
|---|---|---|---|---|---|---|---|
| 0+0 → 15+0 | — | — | — | 168.8 ps | 164.2 ps | +2.8 % | +0.8 % |
| 15+0 → 15+1 | 364.3 ps | 353.7 ps | +3.0 % | 352.5 ps | 402.0 ps | −12.3 % | +33.2 % |
| 10+5 → 10+6 | 270.2 ps | 262.2 ps | +3.1 % | 258.4 ps | 310.5 ps | −16.8 % | +30.1 % |
| 0+15 → 1+15 | 364.3 ps | 353.6 ps | +3.0 % | 352.5 ps | 402.0 ps | −12.3 % | +33.6 % |

**8-bit ripple-carry adder**, worst-case carry 255+0 → 255+1:

| signal | gate | transistor | error | energy error |
|---|---|---|---|---|
| carry-out | 740.5 ps | 716.1 ps | +3.4 % | +29.6 % |
| sum MSB | 728.7 ps | 764.2 ps | −4.6 % | |

STA with the rise-edge table matches the event simulator to about 1 %:
368.8 ps vs 364.3 ps (4-bit) and 745.0 ps vs 740.5 ps (8-bit),
confirming the event engine. The conservative both-edge walk used for the
certificate gives 399.0 ps and 815.1 ps; it is ~10 % above the event result
where the actual transition uses the faster edge.

### Energy

The gate-level model is a conservative upper bound on the adder energy
(+30 %), but only 4–8 % high on the chain. The reason is the same
state-independence that limits any `½CV²`-style table: the switching energy
is characterised for a cell whose internal nodes are fully settled and
whose input edge has the table's slew. In the ripple adder many cells have
partial internal-node swings and short input edges, so the transistor
supply charge is below the sum of the isolated-cell table entries. The
chain, with one clean transition per stage, is much closer. This is a
model-form error, not a table interpolation error.

## 5. The error certificate

The accumulated interpolation certificate along the inverter chain is
`Σ ε_i d_i` over the cells on the path. For the 8-stage chain the test
`certificate_covers_critical_path_error` reports

| path bound (sum) | RSS bound | gate delay | transistor delay | absolute error |
|---|---|---|---|---|
| 9.5 ps | 3.5 ps | 176.7 ps | 169.9 ps | 6.8 ps |

so the transistor-level result **lies inside the accumulated
interpolation bound** for the inverter chain.

For the adder the interpolation bound is *marginally* below the total
fold error:

| circuit | interpolation bound | carry error | in bound? |
|---|---|---|---|
| 4-bit adder | 10.4 ps | 10.6 ps | no (by 0.2 ps) |
| 8-bit adder | 21.4 ps | 24.4 ps | no (by 3.0 ps) |

This is the honest boundary of the fold, and it is the second kind of
error, not interpolation: the NLDM delay depends on the *shape* of the
input edge, and a cell driven by another gate sees an RC-like edge rather
than the linear ramp used in characterisation. The table error itself is
bounded by the certificate; the extra term is the single-edge model form.
It is small (~0.4 ps per full-adder stage) and it grows with path length,
so on long combinational paths the certificate stops covering the total
error before it stops covering the interpolation error. That is exactly
where the gate-level model should unfold to transistor level.

## 6. Unfold indicators and limits

The gate-level model is used only while its validity indicators hold:

- **input slew or load leaves the characterised box** — `TimingTable::contains`
  returns false and the caller unfolds (`unfold_indicator_fires_outside_box`).
  This fired for the un-terminated chain: its last stage drove a wire-only
  0.2 fF load, below the 1 fF grid edge, and the clamped table overestimated
  that stage's delay by 16 %. Terminating the chain removed the
  discrepancy.
- **the accumulated interpolation bound is compared with the observed fold
  error** — where the model-form term pushes the total error past the bound
  (long adder paths), transistor level is required.
- **signal integrity** — glitches shorter than a gate delay, crosstalk and
  coupling are not represented at all: the event simulator evaluates
  settled Boolean values and the tables were characterised one cell at a
  time. These need transistor level (or a coupling-aware extension).

Other limitations, in order of size:

1. Lumped, bias-averaged capacitances (the same approximation as
   `analysis/semiconductor.md` §7); a charge-conserving terminal-charge
   model is the next fold.
2. State-independent switching energy (the +30 % adder overestimate).
3. Single-edge NLDM (the model-form term above).
4. STA takes the later edge at every gate, so it does not track the
   transition direction along a path; it is used only as a conservative
   check of the event simulator.
5. Leakage is a switch-level estimate, not a solved off-state.

## 7. Run

```
cd fluid
cargo test -p semiconductor --test gates          # ~13 s
cargo run  -p semiconductor --example adder        # full tables, ~2 min
```

The tests are small (2–9 cells, coarse time steps where accuracy allows):
the whole gate suite runs in about 13 s, well under the 30 s budget.
