//! Chip design from Verilog to signed-off numbers, on cells that were
//! simulated transistor by transistor.
//!
//! ```text
//!  Verilog (rtl) ─► synthesis (synth) ─► gate netlist ─► timing, power (sta, power)
//!        └──────── equivalence (equiv): the netlist is the design ───┘
//!  cells: semiconductor::gates (INV NAND2 NOR2 XOR2 FA), semiconductor::seq (DFF)
//! ```

pub mod equiv;
pub mod rtl;
pub mod synth;
pub mod flat;
pub mod library;
pub mod sta;
pub mod power;
pub mod corners;
