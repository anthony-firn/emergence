//! µMAG standard problem #4 (NIST) as a reusable run.
//!
//! A 500 x 125 x 3 nm permalloy film is saturated along `[1,1,1]`, relaxed
//! to zero field to obtain the equilibrium s-state, and then reversed by
//! the standard field 1, `mu0 H = (-24.6, 4.3, 0) mT`.  The reference
//! solutions (OOMMF and Albuquerque et al.) give
//! `first <m_y> zero crossing ~ 0.197 ns` from `s-state (0.9666, 0.1260, 0)`;
//! see `analysis/magnetism.md`.
//!
//! Everything here is the same physics as the example; keeping it in the
//! library lets the fast test and the full example share it.
//!
//! Because a 200 x 50 x 1 grid needs tens of minutes, the reversal is split
//! into resumable segments: [`Sp4State`] carries the accumulated time,
//! crossings and samples and can be written to and read from a binary
//! checkpoint together with the grid magnetisation
//! ([`save_checkpoint`], [`load_checkpoint`]).

use crate::demag::DemagFft;
use crate::geometry::V3;
use crate::llg3d::{Llg3d, LlgScratch};
use crate::micromag::{GAMMA, MU0};
use std::fs::File;
use std::io::{self, Read, Write};

/// NIST µMAG problem #4 material and geometry (permalloy).
pub const A_EX: f64 = 1.3e-11;
pub const MS: f64 = 8.0e5;
pub const ALPHA: f64 = 0.02;
pub const LX: f64 = 500e-9;
pub const LY: f64 = 125e-9;
pub const LZ: f64 = 3e-9;

/// Reference initial s-state `(<m_x>, <m_y>, <m_z>)` and the first zero
/// crossings (ns) of the OOMMF reference solution.
pub const REF_S_STATE: V3 = [0.9666, 0.1260, 0.0];
pub const REF_MX_CROSS_NS: f64 = 0.1384;
pub const REF_MY_CROSS_NS: f64 = 0.1967;

const CHECKPOINT_MAGIC: &[u8; 8] = b"MUMAGSP4";

/// Relaxation tolerance on the maximum torque `|m x H_eff|` (A/m).
pub const RELAX_TOL: f64 = 1.0e2;

/// Outcome of one standard-problem-4 run.
pub struct Sp4Run {
    pub s_state: V3,
    pub m_final: V3,
    pub mx_cross_ns: f64,
    pub my_cross_ns: f64,
    pub samples: Vec<(f64, V3)>,
    pub grid: [usize; 3],
}

/// Accumulated state of a standard-problem-4 reversal, so a long run can be
/// split into resumable segments.
#[derive(Clone)]
pub struct Sp4State {
    /// Time (s) already integrated in the reversal.
    pub t: f64,
    /// Last averaged magnetisation components, for crossing detection.
    pub mx_prev: f64,
    pub my_prev: f64,
    pub mx_cross_ns: f64,
    pub my_cross_ns: f64,
    /// Whether the heavy-damped s-state relaxation has converged.
    pub relaxed: bool,
    pub s_state: V3,
    /// Samples `(t_ns, <m>)` accumulated so far.
    pub samples: Vec<(f64, V3)>,
}

impl Sp4State {
    pub fn new(s_state: V3) -> Self {
        Self {
            t: 0.0,
            mx_prev: s_state[0],
            my_prev: s_state[1],
            mx_cross_ns: f64::NAN,
            my_cross_ns: f64::NAN,
            relaxed: true,
            s_state,
            samples: Vec::new(),
        }
    }

    /// A state whose s-state relaxation has not finished; only `mx_prev`,
    /// `my_prev` and `s_state` are filled in when the relaxation converges.
    pub fn new_unrelaxed() -> Self {
        Self {
            t: 0.0,
            mx_prev: 1.0,
            my_prev: 0.0,
            mx_cross_ns: f64::NAN,
            my_cross_ns: f64::NAN,
            relaxed: false,
            s_state: [0.0; 3],
            samples: Vec::new(),
        }
    }
}

/// Build the solver for an `nx x ny x nz` grid (constant material and field
/// parameters, no relaxation yet).
pub fn build(nx: usize, ny: usize, nz: usize) -> Llg3d {
    let dx = LX / nx as f64;
    let dy = LY / ny as f64;
    let dz = LZ / nz as f64;
    let mut s = Llg3d::new(nx, ny, nz, dx, dy, dz, A_EX, 0.0, MS, ALPHA);
    s.demag = Some(DemagFft::new(nx, ny, nz, [dx, dy, dz]));
    s.easy = [1.0, 0.0, 0.0];
    s
}

/// Heavy-damped relaxation time step (s) on an `nx` grid.
pub fn relaxation_dt(nx: usize) -> f64 {
    let dx = LX / nx as f64;
    0.1 * dx * dx * MU0 * MS / (2.0 * A_EX * GAMMA)
}

/// Relax to the equilibrium s-state (saturate along `[1,1,1]`, then damp).
pub fn relax_s_state(s: &mut Llg3d, ws: &mut LlgScratch) -> V3 {
    s.set_uniform([1.0, 1.0, 1.0]);
    s.relax_ws(relaxation_dt(s.nx), RELAX_TOL, 200_000, ws);
    s.average()
}

/// One resumable chunk of the heavy-damped s-state relaxation.  Runs at
/// most `max_steps` steps and returns `(steps used, final max torque)`.  A
/// returned torque below [`RELAX_TOL`] means the s-state has converged.
pub fn relax_chunk(s: &mut Llg3d, ws: &mut LlgScratch, max_steps: usize) -> (usize, f64) {
    let dt = relaxation_dt(s.nx);
    let used = s.relax_ws(dt, RELAX_TOL, max_steps, ws);
    (used, s.max_torque_ws(ws))
}

/// Explicit RK4 time step (s) used for the reversal on an `nx` grid.
pub fn reversal_dt(nx: usize) -> f64 {
    let dx = LX / nx as f64;
    let c_ex = 2.0 * A_EX / (MU0 * MS);
    0.05 * dx * dx / (c_ex * GAMMA)
}

/// Integrate the reversal from `st.t` to `t_end`, appending to `st`.
pub fn integrate(s: &mut Llg3d, ws: &mut LlgScratch, t_end: f64, st: &mut Sp4State) {
    let dt = reversal_dt(s.nx);
    s.h_ext = [-24.6e-3 / MU0, 4.3e-3 / MU0, 0.0];
    s.alpha = ALPHA;
    if t_end <= st.t {
        return;
    }
    let steps = ((t_end - st.t) / dt).ceil() as usize;
    let save_every = (steps / 12).max(1);
    for step in 0..steps {
        s.step_ws(dt, ws);
        let t = st.t + (step + 1) as f64 * dt;
        let avg = s.average();
        if st.mx_cross_ns.is_nan() && st.mx_prev >= 0.0 && avg[0] < 0.0 {
            st.mx_cross_ns = t / 1e-9;
        }
        if st.my_cross_ns.is_nan() && st.my_prev >= 0.0 && avg[1] < 0.0 {
            st.my_cross_ns = t / 1e-9;
        }
        st.mx_prev = avg[0];
        st.my_prev = avg[1];
        if step % save_every == 0 {
            st.samples.push((t / 1e-9, avg));
        }
    }
    st.t += steps as f64 * dt;
}

/// Relax to the s-state and integrate the reversal for `t_end` seconds.
pub fn run(nx: usize, ny: usize, nz: usize, t_end: f64) -> Sp4Run {
    let mut s = build(nx, ny, nz);
    let mut ws = s.make_scratch();
    let s_state = relax_s_state(&mut s, &mut ws);
    let mut st = Sp4State::new(s_state);
    integrate(&mut s, &mut ws, t_end, &mut st);
    Sp4Run {
        s_state,
        m_final: s.average(),
        mx_cross_ns: st.mx_cross_ns,
        my_cross_ns: st.my_cross_ns,
        samples: st.samples,
        grid: [nx, ny, nz],
    }
}

/// Write the grid magnetisation and accumulated reversal state to `path`.
pub fn save_checkpoint(path: &str, s: &Llg3d, st: &Sp4State) -> io::Result<()> {
    let mut f = File::create(path)?;
    f.write_all(CHECKPOINT_MAGIC)?;
    for v in [s.nx as u64, s.ny as u64, s.nz as u64] {
        f.write_all(&v.to_le_bytes())?;
    }
    for v in [
        st.t,
        st.mx_prev,
        st.my_prev,
        st.mx_cross_ns,
        st.my_cross_ns,
        st.s_state[0],
        st.s_state[1],
        st.s_state[2],
    ] {
        f.write_all(&v.to_le_bytes())?;
    }
    f.write_all(&(st.relaxed as u64).to_le_bytes())?;
    f.write_all(&(st.samples.len() as u64).to_le_bytes())?;
    for (t, m) in &st.samples {
        f.write_all(&t.to_le_bytes())?;
        for c in 0..3 {
            f.write_all(&m[c].to_le_bytes())?;
        }
    }
    for m in &s.m {
        for c in 0..3 {
            f.write_all(&m[c].to_le_bytes())?;
        }
    }
    Ok(())
}

fn read_u64(f: &mut File, buf: &mut [u8; 8]) -> io::Result<u64> {
    f.read_exact(buf)?;
    Ok(u64::from_le_bytes(*buf))
}

fn read_f64(f: &mut File, buf: &mut [u8; 8]) -> io::Result<f64> {
    f.read_exact(buf)?;
    Ok(f64::from_le_bytes(*buf))
}

/// Read a checkpoint written by [`save_checkpoint`] and rebuild the solver.
pub fn load_checkpoint(path: &str) -> io::Result<(Llg3d, Sp4State)> {
    let mut f = File::open(path)?;
    let mut magic = [0u8; 8];
    f.read_exact(&mut magic)?;
    if &magic != CHECKPOINT_MAGIC {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "not an SP4 checkpoint",
        ));
    }
    let mut u = [0u8; 8];
    let nx = read_u64(&mut f, &mut u)? as usize;
    let ny = read_u64(&mut f, &mut u)? as usize;
    let nz = read_u64(&mut f, &mut u)? as usize;
    let t = read_f64(&mut f, &mut u)?;
    let mx_prev = read_f64(&mut f, &mut u)?;
    let my_prev = read_f64(&mut f, &mut u)?;
    let mx_cross_ns = read_f64(&mut f, &mut u)?;
    let my_cross_ns = read_f64(&mut f, &mut u)?;
    let s_state = [
        read_f64(&mut f, &mut u)?,
        read_f64(&mut f, &mut u)?,
        read_f64(&mut f, &mut u)?,
    ];
    let relaxed = read_u64(&mut f, &mut u)? != 0;
    let nsamp = read_u64(&mut f, &mut u)? as usize;
    let mut samples = Vec::with_capacity(nsamp);
    for _ in 0..nsamp {
        let ts = read_f64(&mut f, &mut u)?;
        let m = [
            read_f64(&mut f, &mut u)?,
            read_f64(&mut f, &mut u)?,
            read_f64(&mut f, &mut u)?,
        ];
        samples.push((ts, m));
    }
    let mut s = build(nx, ny, nz);
    for cell in s.m.iter_mut() {
        for c in 0..3 {
            cell[c] = read_f64(&mut f, &mut u)?;
        }
    }
    let st = Sp4State {
        t,
        mx_prev,
        my_prev,
        mx_cross_ns,
        my_cross_ns,
        relaxed,
        s_state,
        samples,
    };
    Ok((s, st))
}
