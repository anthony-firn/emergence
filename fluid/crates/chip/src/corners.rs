//! A design's timing and power anywhere in the supply × temperature box,
//! from libraries characterised once on a 3 × 3 grid.
//!
//! Characterising a library from the transistor level takes seconds, and a
//! project asks for many temperatures (the junction temperature is solved
//! together with the leakage it causes). So the transistor level is folded
//! once: the design is analysed exactly at the grid corners, and the
//! results interpolated quadratically in each direction, leakage in its
//! logarithm (it is exponential in temperature). The indicator of the
//! interpolation error is the difference between the quadratic and the
//! linear interpolant; held-out corners check it in `tests/corners.rs`.

use crate::library::{self, Corner};
use crate::power;
use crate::sta::{self, Constraints, Timing};
use crate::synth::Netlist;

/// Supply grid (V) around the 1.2 V nominal, and junction-temperature grid
/// (K) spanning the device law's calibrated range.
pub const VDDS: [f64; 3] = [1.08, 1.2, 1.32];
pub const TEMPS: [f64; 3] = [250.0, 337.5, 425.0];

#[derive(Clone, Copy, Debug, Default)]
pub struct Metrics {
    pub min_period: f64,
    pub min_period_error: f64,
    pub hold_slack: f64,
    pub energy_per_cycle: f64,
    pub leakage: f64,
}

impl Metrics {
    fn of(t: &Timing, p: &power::Power) -> Metrics {
        Metrics { min_period: t.min_period, min_period_error: t.min_period_error, hold_slack: t.hold_slack, energy_per_cycle: p.energy_per_cycle, leakage: p.leakage }
    }
    fn fields(&self) -> [f64; 5] {
        [self.min_period, self.min_period_error, self.hold_slack, self.energy_per_cycle, self.leakage.ln()]
    }
    fn from_fields(f: [f64; 5]) -> Metrics {
        Metrics { min_period: f[0], min_period_error: f[1], hold_slack: f[2], energy_per_cycle: f[3], leakage: f[4].exp() }
    }
}

/// Exact analysis at one characterised corner.
pub fn exact(n: &Netlist, c: &Constraints, toggles: &[f64], corner: Corner) -> Result<(Metrics, Timing), String> {
    let lib = library::at(corner)?;
    let t = sta::analyse(n, &lib, c);
    let p = power::analyse(n, &lib, c, toggles);
    Ok((Metrics::of(&t, &p), t))
}

fn lagrange(x: [f64; 3], f: [f64; 3], t: f64) -> f64 {
    let l = |i: usize, j: usize, k: usize| (t - x[j]) * (t - x[k]) / ((x[i] - x[j]) * (x[i] - x[k]));
    f[0] * l(0, 1, 2) + f[1] * l(1, 0, 2) + f[2] * l(2, 0, 1)
}

fn linear(x: [f64; 3], f: [f64; 3], t: f64) -> f64 {
    let i = if t <= x[1] { 0 } else { 1 };
    f[i] + (f[i + 1] - f[i]) * (t - x[i]) / (x[i + 1] - x[i])
}

/// Metrics at (vdd, temp) by interpolation on the grid, and the indicator
/// of the interpolation error (|quadratic − linear| per field, leakage as a
/// relative error).
pub fn at(n: &Netlist, c: &Constraints, toggles: &[f64], vdd: f64, temp: f64) -> Result<(Metrics, Metrics), String> {
    if vdd < VDDS[0] - 1e-12 || vdd > VDDS[2] + 1e-12 {
        return Err(format!("supply {vdd:.3} V is outside the characterised {:.2}–{:.2} V", VDDS[0], VDDS[2]));
    }
    if temp < TEMPS[0] - 1e-9 || temp > TEMPS[2] + 1e-9 {
        return Err(format!("junction temperature {temp:.1} K is outside the characterised {:.0}–{:.0} K", TEMPS[0], TEMPS[2]));
    }
    // The nine corners, in parallel (each characterises once, then caches).
    use rayon::prelude::*;
    let corners: Vec<(usize, usize)> = (0..3).flat_map(|i| (0..3).map(move |j| (i, j))).collect();
    let grid: Vec<Result<[f64; 5], String>> = corners
        .par_iter()
        .map(|&(i, j)| exact(n, c, toggles, Corner { vdd: VDDS[i], temp: TEMPS[j] }).map(|(m, _)| m.fields()))
        .collect();
    let mut g = [[[0.0; 5]; 3]; 3];
    for (&(i, j), r) in corners.iter().zip(grid) {
        g[i][j] = r?;
    }
    let mut quad = [0.0; 5];
    let mut lin = [0.0; 5];
    for f in 0..5 {
        let qt: [f64; 3] = std::array::from_fn(|i| lagrange(TEMPS, [g[i][0][f], g[i][1][f], g[i][2][f]], temp));
        let lt: [f64; 3] = std::array::from_fn(|i| linear(TEMPS, [g[i][0][f], g[i][1][f], g[i][2][f]], temp));
        quad[f] = lagrange(VDDS, qt, vdd);
        lin[f] = linear(VDDS, lt, vdd);
    }
    let value = Metrics::from_fields(quad);
    let mut err = [0.0; 5];
    for f in 0..5 {
        err[f] = (quad[f] - lin[f]).abs();
    }
    // Leakage's indicator is in ln: report it as a relative error.
    let indicator = Metrics {
        min_period: err[0],
        min_period_error: err[1],
        hold_slack: err[2],
        energy_per_cycle: err[3],
        leakage: err[4].exp() - 1.0,
    };
    Ok((value, indicator))
}
