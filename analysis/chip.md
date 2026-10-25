# Chip design: Verilog to signed-off numbers on transistor-simulated cells

Crates: `fluid/crates/chip` (flow), `fluid/crates/semiconductor` (devices,
cells, flip-flop). Project analysis: `chip_design`. Example:
`fluid/projects/mini-cpu/`.

## The ladder

```
TCAD (drift–diffusion, semiconductor.md)
  └─ fold: EKV compact law, temperature ratios fitted to TCAD     (temperature)
       └─ cells: INV NAND2 NOR2 XOR2 FA + D flip-flop, transistor-level transients
            └─ fold: NLDM tables with held-out certificates, per supply × temperature
                 └─ Verilog → synthesis → gate netlist (proven equal to the RTL)
                      └─ static timing, power → project graph (package, heat sink, budgets)
```

Every fold is checked against the rung below it. The checks are the rest
of this document.

## 1. Temperature of the device law

The cells use the all-region EKV law calibrated to the 150 nm technology
at 300 K (`analysis/gates.md`). Its temperature dependence is folded from
this repository's drift–diffusion transistors (full mobility model, |V_ds|
= 1.2 V, `semiconductor/examples/ekv_temperature.rs`).

A plain 3-parameter EKV law cannot follow the TCAD curve shape: it has no
velocity saturation or mobility degradation. A direct fit at each
temperature leaves 10–16 % rms error and misses I_on by 20–28 %. What the
cells depend on is how currents *change*. So the fold matches the ratios
I_on(T)/I_on(300 K) and I_off(T)/I_off(300 K), through two physical
parameters per polarity:

* the threshold's temperature slope;
* the exponent `a` in I_s ∝ (T/300)^a (mobility and V_T² together).

They are fitted on 250, 300, 350 and 400 K, with 275, 325, 375 and 425 K
held out.

| | dV_th/dT | a | held-out I_on | held-out I_off (≤ 375 K) | I_off at 425 K |
|---|---|---|---|---|---|
| nMOS | −0.586 mV/K | 0.98 | within 3.1 % | within 2.2 % | −8.9 % |
| pMOS | −0.620 mV/K | 0.94 | within 1.8 % | within 2.9 % | −10.9 % |

Over the range, TCAD leakage grows 2370× (nMOS, 250 → 425 K) and on-current
falls 28 %. `semiconductor/tests/temperature.rs` records the TCAD currents
and checks the coefficients against them; the TCAD run itself takes about
15 minutes. Outside 250–425 K the law refuses.

## 2. The D flip-flop

An 18-transistor static master–slave flip-flop (four transmission gates,
five inverters; Weste & Harris Fig. 10.19), simulated with the same law and
backward-Euler engine as the combinational cells. It is characterised the
industry way:

* clock-to-Q NLDM tables with energy;
* setup and hold by bisection on 10 % clock-to-Q push-out;
* pin capacitances, energy per clock cycle and leakage.

| at 1.2 V, 300 K | value |
|---|---|
| clock-to-Q rise / fall (80 ps clock, 4 fF) | 92.0 / 103.3 ps |
| setup rise / fall | 57.3 / 30.8 ps |
| hold rise / fall | 2.8 / −36.2 ps |
| D / CLK pin capacitance | 4.33 / 7.73 fF |
| energy per clock cycle | 13.0 fJ |
| leakage | 0.71 nA |

A transistor-level run over six clock cycles captures D on every rising
edge and holds it otherwise. A full characterisation takes 2.3 s.

## 3. Verilog and synthesis, proven equivalent

`chip::rtl` reads a synthesizable Verilog-2001 subset. It covers ANSI
ports, wire/reg, assign (including concatenated targets), localparam,
`always @(posedge clk)` with if/else and case, sized literals, selects,
concatenation and replication, and the usual operators. It refuses
latches, asynchronous resets, multiple clocks, `initial`, x/z and
combinational loops, each with a reason.

`chip::synth` bit-blasts into a hash-consed graph with constant folding
and sharing, mapped one-to-one onto the characterised cells:

* full-adder ripple chains for `+`, `−`, `*` and comparisons;
* barrel shifters for shifts;
* NAND multiplexers for selection;
* priority chains for `case`.

**Cone resynthesis** then takes every node whose function depends on at
most 6 inputs, computes its truth table bit-parallel, and rebuilds it by
Shannon expansion. Complementary cofactors become XOR, constant ones
AND/OR, and identical subfunctions are shared. The rebuilt cone is kept
where it is shallower by the characterised cell delays.

`chip::equiv` checks the result. An independent word-level interpreter
runs the RTL, executing only the true branch and using the CPU's adder, and
is compared with the gate netlist on outputs and next state. Registers map
one-to-one, so agreement on every (input, state) point is a proof. Up to
2^22 points are enumerated; beyond that the check is sampled and says so.

| design | cells (after resynthesis) | transistors | equivalence |
|---|---|---|---|
| 8-bit accumulator ALU (add/and/xor/sub, sync reset) | 14 FA, 116 NAND2, … 8 DFF | 1540 | **proven**, 524 288 points |
| 4×4 multiplier | 7 FA, 43 NAND2, … | 624 | proven, 256 |
| 8-bit barrel shifter | 126 NAND2, 18 INV | 540 | proven, 4096 |
| 6-bit comparators and max | 9 FA, … | 668 | proven, 4096 |
| traffic-light FSM (case, timer) | 6 DFF, … | 588 | proven, 256 |
| mini8 CPU | 53 DFF, 7 FA, 403 NAND2, … | 3696 | sampled, 20 000 (54 state bits) |

Changing one gate's kind in a proven netlist is caught.

**mini8** (`projects/mini-cpu/mini8.v`) is an 8-bit accumulator CPU with
four registers, an ALU, jumps, an output port and a 16-word ROM. The ROM
holds a Fibonacci program, which the RTL and the gate netlist run
identically for 120 cycles (0, 1, 2, 3, 5, 8, 13, …).

Resynthesis turns the ROM's 13-deep `?:` decode chain into a 4-level tree:

| mini8 | before | after |
|---|---|---|
| transistors | 5642 | 3696 (−34 %) |
| fmax (1.2 V, 300 K, 1 fF wires) | 321 MHz | 489 MHz (+52 %) |
| energy per cycle (25 % activity) | 1.32 pJ | 1.16 pJ |

The choice is depth-driven, so arithmetic designs can grow (the
accumulator ALU gains 15 %).

## 4. Clocked static timing, against transistors

`chip::sta`:

* launches paths at flip-flop clock-to-Q and at input delays;
* captures at setup and hold and at output delays;
* carries rise and fall separately through each cell's unateness
  (inverting gates, a positive-unate carry, non-unate XOR and sum);
* computes late and early arrivals, and the critical path with its
  accumulated held-out interpolation certificate.

**Check against the transistor level.** A register → 2-bit adder →
register pipeline is flattened to one transistor network with a driven
clock, about 230 transistors including the flip-flops. The clock period is
bisected for the worst carry ripple (00+00, then 11+01):

| | shortest period |
|---|---|
| STA | 340.4 ± 2.2 ps |
| every transistor simulated | 337.5 ps |

STA is safe and within 0.9 %. With zero input delay STA reports a 2.8 ps
hold violation, which is real: an input changing exactly at the clock edge
would break the flip-flop's hold requirement.

## 5. Power, against transistors

Activity comes from cycle simulation of the netlist. Each output
transition is charged the supply energy characterised at transistor level
for that cell at the net's slew and load, which includes short-circuit
current. Each flip-flop pays its clock energy every cycle, and leakage is
each cell's worst input state.

Check: the same pipeline over 10 random cycles at 1 GHz, with the supply
current of every transistor integrated. The gate-level estimate is
**14.8 % low**: it leaves out glitches and internal nodes that switch
without the output switching. The project analysis therefore claims ±25 %
on power. A second design, out of sample, would make that claim firmer.

## 6. Corners: supply × temperature

Characterising a library takes about 2.2 s for the combinational cells
(after parallelising: it was 8.9 s) plus 2.3 s for the flip-flop. A project
asks for many temperatures, because the junction temperature is solved
together with the leakage it causes.

So the design is analysed exactly on a 3 × 3 grid, supply 1.08/1.2/1.32 V
by temperature 250/337.5/425 K, and interpolated quadratically, leakage in
its logarithm. The indicator is |quadratic − linear|. Libraries are
cached on disk, keyed by a hash of every constant they depend on.

Check on mini8 at four corners characterised directly:

| corner | period error | energy error | leakage error |
|---|---|---|---|
| 1.14 V, 293.75 K | +3.0 ps (indicator 4.9) | −0.08 % | −10.7 % |
| 1.26 V, 381.25 K | −3.2 ps (8.6) | +0.09 % | +8.9 % |
| 1.20 V, 300 K | −1.1 ps (16.2) | −0.04 % | −9.8 % |
| 1.10 V, 400 K | +7.4 ps (**5.1**) | +0.06 % | +7.9 % |

The period indicator was optimistic once, by 1.43×, so the project
analysis claims twice the indicator. The energy and leakage indicators
covered every case.

## 7. In the project graph

`chip_design` (config `rtl`) takes supply, junction temperature, clock,
activity, wire load, input delay and instance count. It gives fmax, setup
and hold slack, switching and leakage power, and energy per cycle, each
with its claim, and names the critical path.

`projects/mini-cpu`: 20 000 mini8 cores at 400 MHz. Silicon owns supply
(±5 %), activity and wire load (±50 %, pre-layout). Packaging owns
junction-to-case and case-to-ambient. Systems owns ambient, clock and
budget. The junction temperature closes the leakage loop through the
package:

| requirement | result |
|---|---|
| junction ≤ 105 °C | 77.1 °C [62.9, 97.8], met; loop settled in 10 iterations, existence proven |
| power ≤ 15 W | 9.45 W [6.16, 13.58], met; leakage 51 mW |
| no hold violations | 191 ps [166, 222], met |
| fmax ≥ 1.1 × clock | 1.12 [0.93, 1.35], **undecided**: regulator tolerance 51 %, wire estimate 19 %, packaging (via temperature) 3 % |

The critical path at 1.2 V, 337.5 K has 34 stages: instruction decode,
the register-file multiplexer, a 7-stage ripple carry, the ALU select and
the zero flag.

## Honest limits

* **One technology.** The 150 nm scaled process of `semiconductor.md`,
  with lumped capacitances. There is no layout: no placement, routing,
  parasitic extraction or DRC/LVS. Wires are a capacitance per
  connection.
* **Ideal clock.** No skew and no clock tree, and one clock domain.
* **Synthesis is basic.** No carry-lookahead or carry-select adders, no
  retiming, no technology-independent optimisation beyond 6-input cones,
  and no sizing (one drive strength per cell).
* **Sequential equivalence assumes the same register encoding.**
  Equivalence is combinational over (input, state), which proves the
  synthesis but would not prove a re-encoded state machine.
* **Power omits glitches and internal nodes** (−15 % against
  transistors on one design), and leakage is worst-state, so it is
  conservative.
* **Temperature fold.** Fitted on TCAD ratios; leakage is 11 % off at
  425 K. The 300 K absolute leakage is the cell calibration's (EKV
  I_off 1.4e-10 A/µm against TCAD 3.1e-10), a known 2.2× gap inherited
  from `gates.md`.
