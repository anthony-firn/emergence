//! µMAG standard problem #4 (NIST): field reversal of a 500 x 125 x 3 nm
//! permalloy film.  See [`magnetism::sp4`] and `analysis/magnetism.md`.
//!
//! Run with an optional grid, end time, checkpoint and relaxation-chunk
//! size:
//!
//! ```text
//! cargo run --release --example mumag_sp4 -- [nx ny nz] [t_end_ns] [checkpoint] [relax_chunk_steps]
//! ```
//!
//! The default `64 16 1` (7.8 nm cells) finishes in seconds.  A resolved
//! grid such as `128 32 1` (3.9 nm cells) takes under a minute.  The finest
//! `200 50 1` (2.5 nm cells) is run as resumable segments: give a
//! checkpoint file and the run relaxes the s-state in bounded chunks and
//! integrates to `t_end_ns`, saving the grid state and the accumulated
//! crossing times after each invocation.  Re-run with a larger `t_end_ns`
//! (or simply again, if the relaxation has not converged) to continue.
//!
//! Reference solution (OOMMF, from https://www.ctcms.nist.gov/~rdm/std4/):
//! s-state `(0.9666, 0.1260, 0)`, first `<m_x>` zero crossing 0.1384 ns,
//! first `<m_y>` zero crossing 0.1967 ns.

use magnetism::sp4;
use std::env;
use std::path::Path;
use std::time::Instant;

fn main() {
    let args: Vec<String> = env::args().collect();
    let (nx, ny, nz) = if args.len() >= 4 {
        (
            args[1].parse().unwrap(),
            args[2].parse().unwrap(),
            args[3].parse().unwrap(),
        )
    } else {
        (64, 16, 1)
    };
    let t_end = if args.len() >= 5 {
        args[4].parse::<f64>().unwrap() * 1e-9
    } else {
        0.6e-9
    };
    let checkpoint: Option<&str> = args.get(5).map(|s| s.as_str());
    let relax_chunk_steps: usize = args.get(6).map(|s| s.parse().unwrap()).unwrap_or(5_000);

    let t0 = Instant::now();
    let (mut s, mut st) = if let Some(path) = checkpoint {
        if Path::new(path).exists() {
            let (s, st) = sp4::load_checkpoint(path).expect("failed to read checkpoint");
            println!(
                "resuming {} x {} x {} at t = {:.4} ns, relaxed = {}, {} samples",
                s.nx,
                s.ny,
                s.nz,
                st.t / 1e-9,
                st.relaxed,
                st.samples.len()
            );
            (s, st)
        } else {
            let mut s = sp4::build(nx, ny, nz);
            s.set_uniform([1.0, 1.0, 1.0]);
            (s, sp4::Sp4State::new_unrelaxed())
        }
    } else {
        let mut s = sp4::build(nx, ny, nz);
        let mut ws = s.make_scratch();
        let s_state = sp4::relax_s_state(&mut s, &mut ws);
        (s, sp4::Sp4State::new(s_state))
    };

    let mut ws = s.make_scratch();
    if !st.relaxed {
        // Bounded relaxation chunk so an invocation stays short.
        let (used, torque) = sp4::relax_chunk(&mut s, &mut ws, relax_chunk_steps);
        println!(
            "relax chunk: {} steps, max torque {:.4e} (tol {:.1e})",
            used,
            torque,
            sp4::RELAX_TOL
        );
        if torque >= sp4::RELAX_TOL {
            if let Some(path) = checkpoint {
                sp4::save_checkpoint(path, &s, &st).expect("failed to write checkpoint");
                println!("wrote checkpoint {path} (still relaxing)");
            } else {
                println!("relaxation not converged; re-run with a checkpoint to continue");
            }
            return;
        }
        st.relaxed = true;
        st.s_state = s.average();
        st.mx_prev = st.s_state[0];
        st.my_prev = st.s_state[1];
        println!(
            "relaxed s-state <m> = ({:.4}, {:.4}, {:.4})",
            st.s_state[0], st.s_state[1], st.s_state[2]
        );
        if let Some(path) = checkpoint {
            sp4::save_checkpoint(path, &s, &st).expect("failed to write checkpoint");
        }
    }

    sp4::integrate(&mut s, &mut ws, t_end, &mut st);
    let elapsed = t0.elapsed().as_secs_f64();
    if let Some(path) = checkpoint {
        sp4::save_checkpoint(path, &s, &st).expect("failed to write checkpoint");
        println!("wrote checkpoint {path} at t = {:.4} ns", st.t / 1e-9);
    }

    let (gx, gy, gz) = (s.nx, s.ny, s.nz);
    println!(
        "grid {} x {} x {} ({} nm cells)",
        gx,
        gy,
        gz,
        sp4::LX / gx as f64 * 1e9
    );
    println!(
        "s-state <m> = ({:.4}, {:.4}, {:.4})   reference ({:.4}, {:.4}, {:.4})",
        st.s_state[0],
        st.s_state[1],
        st.s_state[2],
        sp4::REF_S_STATE[0],
        sp4::REF_S_STATE[1],
        sp4::REF_S_STATE[2]
    );
    println!(
        "first <m_x> zero crossing {:.4} ns   reference {:.4} ns",
        st.mx_cross_ns,
        sp4::REF_MX_CROSS_NS
    );
    println!(
        "first <m_y> zero crossing {:.4} ns   reference {:.4} ns",
        st.my_cross_ns,
        sp4::REF_MY_CROSS_NS
    );
    let m_final = s.average();
    println!(
        "final <m> = ({:.4}, {:.4}, {:.4})   [{:.1} s in this invocation]",
        m_final[0], m_final[1], m_final[2], elapsed
    );
    println!("   t(ns)      <mx>      <my>      <mz>");
    for (t, m) in &st.samples {
        println!("  {:.3}    {:+.4}   {:+.4}   {:+.4}", t, m[0], m[1], m[2]);
    }
}
