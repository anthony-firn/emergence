//! The developing-channel solver against its two exact limits.

use project::thermal::{march, similarity};

#[test]
fn blasius_and_pohlhausen() {
    // f''(0) = 0.332057336215... (Blasius; e.g. Schlichting, Boundary-Layer
    // Theory, 8th ed., Table 7.1).
    let s = similarity::blasius();
    println!("f''(0) = {s:.10}, theta'(0, Pr=0.7) = {:.6}", similarity::theta_wall(0.7));
    assert!((s - 0.332_057_336).abs() < 1e-8, "{s}");
    // Pr = 1: the energy equation is the momentum equation (θ = f'), so
    // θ'(0) = f''(0) exactly.
    assert!((similarity::theta_wall(1.0) - s).abs() < 1e-8);
    // Pr = 0.7: Pohlhausen's 0.2927 (Schlichting Table 12.1 quotes 0.293);
    // the correlation 0.332 Pr^(1/3) = 0.2948 is 0.7 % high.
    let t = similarity::theta_wall(0.7);
    assert!((t - 0.2927).abs() < 1e-4, "{t}");
}

#[test]
fn fully_developed_limit() {
    // Far downstream: Poiseuille flow (u = 1.5 at the mid-plane,
    // dP/dX = −3 in these units, i.e. dp/dx = −12 μU/b²) and Nu_Dh = 7.5407
    // for uniform wall temperature (Shah & London 1978).
    let pr = 0.7;
    let st = march(pr, 16.0 * pr * 0.5, 80, 800);
    let end = st.last().unwrap();
    let before = st[st.len() - 11];
    let dpdx = (end.p - before.p) / (end.x - before.x);
    println!("dP/dX = {dpdx:.6}, u_c = {:.6}, Nu = {:.6}", end.u_center, end.nu_local);
    assert!((dpdx + 3.0).abs() < 2e-3, "dP/dX = {dpdx}");
    assert!((end.u_center - 1.5).abs() < 1e-3, "u_c = {}", end.u_center);
    assert!((end.nu_local / 7.5407 - 1.0).abs() < 2e-3, "Nu = {}", end.nu_local);
}

#[test]
fn leading_edge_limit() {
    // Near the inlet each wall is a flat plate: Nu_Dh,x = 4 θ'(0) / √X.
    // The channel departs from the plate at O(√X) (displacement of the
    // core), so the ratio, extrapolated linearly in √X to X = 0, must return
    // the similarity value.
    let pr = 0.7;
    let exact = 4.0 * similarity::theta_wall(pr);
    let st = march(pr, 4e-3, 160, 1600);
    let at = |x: f64| {
        let s = st.iter().min_by(|a, b| (a.x - x).abs().total_cmp(&(b.x - x).abs())).unwrap();
        (s.x.sqrt(), s.nu_local * s.x.sqrt())
    };
    let (s1, r1) = at(1e-3);
    let (s2, r2) = at(4e-3);
    let r0 = r1 - (r2 - r1) / (s2 - s1) * s1;
    println!("ratio at X=1e-3: {:.5}, X=4e-3: {:.5}, extrapolated {:.5}, exact {:.5}", r1 / exact, r2 / exact, r0 / exact, 1.0);
    assert!((r0 / exact - 1.0).abs() < 5e-3, "extrapolated {r0} vs {exact}");
    // and the departure is upward (core acceleration thins the layer).
    assert!(r2 > r1 && r1 > exact);
}
