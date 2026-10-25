//! µMAG standard problem #4: coarse field-reversal test against the
//! published reference time series.
//!
//! The reference curves (OOMMF, Donahue & Porter; Albuquerque, Miltat &
//! Thiaville) were downloaded from https://www.ctcms.nist.gov/~rdm/std4/.
//! For the s-state and field 1 they give a first `<m_x>` zero crossing at
//! 0.1384 ns and a first `<m_y>` zero crossing at 0.1967 ns.  This test is
//! deliberately coarse (32 x 8 x 1 cells, 15.6 nm in-plane) and checks the
//! crossing times to 15 %.

use magnetism::sp4;

fn rel(a: f64, b: f64) -> f64 {
    (a - b).abs() / b.abs()
}

#[test]
fn standard_problem_4_field_reversal() {
    let r = sp4::run(32, 8, 1, 0.6e-9);

    // The relaxed s-state must be close to the reference s-state.  On this
    // coarse grid the end-domain canting is under-resolved, so only the
    // dominant x component is checked tightly.
    assert!(
        r.s_state[0] > 0.94,
        "s-state mx = {} (reference {})",
        r.s_state[0],
        sp4::REF_S_STATE[0]
    );
    assert!(
        rel(r.s_state[1], sp4::REF_S_STATE[1]) < 0.8,
        "s-state my = {} (reference {})",
        r.s_state[1],
        sp4::REF_S_STATE[1]
    );

    // First zero crossings of the averaged magnetisation.
    let ex = rel(r.mx_cross_ns, sp4::REF_MX_CROSS_NS);
    let ey = rel(r.my_cross_ns, sp4::REF_MY_CROSS_NS);
    assert!(
        ex < 0.15,
        "first <mx> crossing {} ns vs reference {} ns ({:.1} %)",
        r.mx_cross_ns,
        sp4::REF_MX_CROSS_NS,
        ex * 100.0
    );
    assert!(
        ey < 0.15,
        "first <my> crossing {} ns vs reference {} ns ({:.1} %)",
        r.my_cross_ns,
        sp4::REF_MY_CROSS_NS,
        ey * 100.0
    );

    // Qualitative reversal features: the layer ends up close to -x and the
    // in-plane magnetisation grows during the first tenth of a nanosecond.
    assert!(r.m_final[0] < -0.6, "final mx = {}", r.m_final[0]);
    let early = r
        .samples
        .iter()
        .find(|(t, _)| *t >= 0.08)
        .map(|(_, m)| m[1])
        .unwrap();
    assert!(early > 0.3, "early <my> = {early}");
}

/// A checkpoint must contain everything needed to resume: the grid
/// magnetisation and the accumulated reversal state.
#[test]
fn checkpoint_round_trips() {
    let mut s = sp4::build(32, 8, 1);
    let mut ws = s.make_scratch();
    let s_state = sp4::relax_s_state(&mut s, &mut ws);
    let mut st = sp4::Sp4State::new(s_state);
    // A few reversal steps, enough to move <m> off the s-state.
    sp4::integrate(&mut s, &mut ws, 5.0e-12, &mut st);
    assert!(st.t > 0.0);

    let path = std::env::temp_dir().join(format!("mumag_ckpt_{}.bin", std::process::id()));
    let path = path.to_str().unwrap().to_string();
    sp4::save_checkpoint(&path, &s, &st).unwrap();
    let (mut s2, mut st2) = sp4::load_checkpoint(&path).unwrap();
    std::fs::remove_file(&path).ok();

    assert_eq!(s2.nx, s.nx);
    assert_eq!(s2.ny, s.ny);
    assert_eq!(s2.nz, s.nz);
    assert!((st2.t - st.t).abs() < 1e-18);
    for (a, b) in s.m.iter().zip(&s2.m) {
        for c in 0..3 {
            assert!((a[c] - b[c]).abs() < 1e-15, "{a:?} vs {b:?}");
        }
    }

    // Continuing each copy for the same interval must agree.
    let end = st.t + 5.0e-12;
    sp4::integrate(&mut s, &mut ws, end, &mut st);
    let mut ws2 = s2.make_scratch();
    sp4::integrate(&mut s2, &mut ws2, end, &mut st2);
    for (a, b) in s.m.iter().zip(&s2.m) {
        for c in 0..3 {
            assert!((a[c] - b[c]).abs() < 1e-13, "{a:?} vs {b:?}");
        }
    }
}
