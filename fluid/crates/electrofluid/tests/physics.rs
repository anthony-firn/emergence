//! Physics benchmarks for `electrofluid`. Each test is small and fast.

use electrofluid::eo::Channel;
use electrofluid::er::{dipole_force, ErSim};
use electrofluid::leaky::LeakyDrop;
use electrofluid::pnp::*;
use std::f64::consts::PI;

#[test]
fn debye_length_matches_published() {
    // Published 1:1 electrolyte Debye lengths in water at 25 °C:
    // ~9.61 nm at 1 mM and ~0.962 nm at 100 mM (see analysis/electrofluid.md).
    let l1 = debye_length(EPS_WATER, T0, 1.0e-3, 1.0);
    let l100 = debye_length(EPS_WATER, T0, 1.0e-1, 1.0);
    assert!((l1 * 1e9 - 9.62).abs() < 0.15, "1 mM: {} nm", l1 * 1e9);
    assert!((l100 * 1e9 - 0.962).abs() < 0.02, "100 mM: {} nm", l100 * 1e9);
}

#[test]
fn gouy_chapman_matches_poisson_boltzmann() {
    let c = 1e-3;
    let lam = debye_length(EPS_WATER, T0, c, 1.0);
    let psi0 = -0.100;
    let (x, psi) = pb_solve(EPS_WATER, T0, c, 1.0, psi0, 8.0 * lam, 4000);
    let mut err = 0.0f64;
    for (i, (&xi, &pi)) in x.iter().zip(psi.iter()).enumerate() {
        if i % 500 == 0 {
            err = err.max((pi - gouy_chapman(xi, lam, psi0, T0, 1.0)).abs());
        }
    }
    assert!(err < 1e-4, "GC error {err} V");
}

#[test]
fn debye_huckel_is_the_linear_limit() {
    let c = 1e-3;
    let lam = debye_length(EPS_WATER, T0, c, 1.0);
    let psi0 = -0.005; // 5 mV, |u| = 0.2
    let (x, psi) = pb_solve(EPS_WATER, T0, c, 1.0, psi0, 8.0 * lam, 2000);
    let mut err = 0.0f64;
    for (i, (&xi, &pi)) in x.iter().zip(psi.iter()).enumerate() {
        if i % 200 == 0 {
            err = err.max((pi - debye_huckel(xi, lam, psi0)).abs());
        }
    }
    assert!(err / psi0.abs() < 5e-3, "DH relative error {}", err / psi0.abs());
}

#[test]
fn grahame_matches_profile_charge() {
    let c = 1e-3;
    let lam = debye_length(EPS_WATER, T0, c, 1.0);
    let psi0 = -0.100;
    let (x, psi) = pb_solve(EPS_WATER, T0, c, 1.0, psi0, 8.0 * lam, 4000);
    let sig = grahame_sigma(EPS_WATER, T0, c, 1.0, psi0);
    let sig_num = sigma_from_profile(&psi, &x, EPS_WATER);
    assert!(((sig - sig_num) / sig).abs() < 1e-3, "Grahame {sig} vs {sig_num}");
}

#[test]
fn electroosmotic_slip_fold() {
    // Helmholtz–Smoluchowski is the λ_D/h → 0 limit of the resolved channel.
    let ch = |c: f64| Channel {
        eps_r: EPS_WATER,
        t: T0,
        c_molar: c,
        z: 1.0,
        zeta: -0.025,
        e_field: 1.0e4,
        mu: MU_WATER,
        h: 60e-9,
    };
    let (uc10, us) = ch(1e-2).centerline(600);
    assert!(((uc10 - us) / us).abs() < 5e-3, "10 mM rel err {}", (uc10 - us).abs() / us.abs());
    let (uc100, us2) = ch(1e-1).centerline(600);
    assert!(((uc100 - us2) / us2).abs() < 1e-6, "100 mM rel err {}", (uc100 - us2).abs() / us2.abs());
    // and it fails for overlapping double layers
    let (uc01, us3) = ch(1e-4).centerline(600);
    assert!((uc01 - us3).abs() / us3.abs() > 0.1);
}

#[test]
fn leaky_drop_prolate_oblate_and_stress_balance() {
    let base = |r: f64, s: f64| LeakyDrop {
        a: 1e-6,
        eps_o: 2.0,
        eps_i: 2.0 * s,
        sigma_o: 1e-9,
        sigma_i: 1e-9 * r,
        mu_o: 1e-3,
        mu_i: 1e-3,
        gamma: 4e-3,
        e0: 1e6,
    };
    // R > S → prolate; R < S → oblate
    assert!(base(10.0, 2.0).prolate());
    assert!(!base(0.1, 2.0).prolate());
    // The P₂ normal-stress balance gives the shape amplitude s in
    // r = a(1 + s P₂). Taylor's deformation parameter is D = (L−B)/(L+B) =
    // 3s/4 exactly, so s = 4D/3. The historical "4/3" is this geometry, not
    // a missing flow: the normal Maxwell stress enters the balance while the
    // flow is driven by the tangential electric traction.
    for r in [0.1, 1.0, 10.0] {
        let d = base(r, 1.5);
        let dan = d.deformation();
        let dbal = d.shape_amplitude();
        assert!((dbal / dan - 4.0 / 3.0).abs() < 1e-10, "s/D = {}", dbal / dan);
        assert!((d.taylor_shape_amplitude() - dbal).abs() < 1e-12);
    }
    // D scales linearly with Ca_E
    let mut d0 = base(10.0, 2.0);
    let d1 = base(10.0, 2.0);
    d0.gamma *= 2.0;
    assert!((d0.deformation() / d1.deformation() - 0.5).abs() < 1e-12);
}

#[test]
fn leaky_drop_circulation_matches_published() {
    // Interfacial tangential velocity u_θ(a,θ) = U_T sin2θ with
    // U_T = −(9/10) ε_oε₀E₀²a/(μ_o(1+λ)) (R−S)/(R+2)², λ = μ_i/μ_o
    // (Vlahovska 2019 review, sidebar "Interfacial velocity", without the
    // spurious extra 1/S). Equator→pole for R>S (prolate A).
    let cases = [
        // a, eps_o, eps_i, sigma_o, sigma_i, mu_o, mu_i, gamma, e0, R, S
        (1e-6, 2.0, 8.0, 1e-9, 1e-9, 1e-3, 1e-3, 4e-3, 1e6, 1.0, 4.0),
        (1e-6, 2.0, 1.0, 1e-9, 1e-9, 1e-3, 1e-3, 4e-3, 1e6, 1.0, 0.5),
        (1e-6, 2.0, 3.0, 1e-9, 1e-10, 1e-3, 1e-3, 4e-3, 1e6, 0.1, 1.5),
        (1e-6, 2.0, 3.0, 1e-9, 1e-8, 1e-3, 5e-3, 4e-3, 1e6, 10.0, 1.5),
    ];
    for (a, eo, ei, so, si, mo, mi, g, e0, r, _s) in cases {
        let d = LeakyDrop { a, eps_o: eo, eps_i: ei, sigma_o: so, sigma_i: si, mu_o: mo, mu_i: mi, gamma: g, e0 };
        let lambda = mi / mo;
        let u_pub = -(9.0 / 10.0) * eo * EPS0 * e0 * e0 * a / (mo * (1.0 + lambda))
            * (d.r_ratio() - d.s_ratio())
            / ((d.r_ratio() + 2.0) * (d.r_ratio() + 2.0));
        assert!(
            (d.circulation_amplitude() - u_pub).abs() / u_pub.abs() < 1e-12,
            "U_T code {} vs pub {}",
            d.circulation_amplitude(),
            u_pub
        );
        // analytic interfacial velocity equals U_T sin2θ
        for &th in &[0.4, 1.1, 2.0] {
            let lhs = d.interfacial_velocity(th);
            let rhs = d.circulation_amplitude() * (2.0 * th).sin();
            let scale = d.circulation_amplitude().abs().max(1e-300);
            assert!((lhs - rhs).abs() / scale < 1e-12);
        }
        // flow direction: equator→pole for R>S (prolate A), pole→equator otherwise
        let ug_eq = d.interfacial_velocity(std::f64::consts::FRAC_PI_4);
        if r > d.s_ratio() {
            assert!(ug_eq < 0.0);
        } else if r < d.s_ratio() {
            assert!(ug_eq > 0.0);
        }
    }
}

#[test]
fn dipole_pair_force_is_attractive_along_field() {
    // aligned spheres (j above i): force on i must point +y (toward j)
    let (fx, fy) = dipole_force(0.0, 2.0);
    assert!(fx.abs() < 1e-12);
    assert!(fy > 0.0);
    // contact value: 12π/r⁴ · 2 = 1.5π in units of ε₀ε_fβ²E²a²
    assert!((fy - 1.5 * PI).abs() < 1e-9, "fy = {fy}");
}

#[test]
fn er_chains_form_and_shear_breaks_them() {
    let n = 36usize;
    let (lx, ly) = (15.0, 15.0);
    let dt = 1e-3;
    let mut sim = ErSim::new(n, lx, ly);
    for _ in 0..2500 {
        sim.step(dt);
    }
    // chains: strongly negative interaction energy and alignment
    assert!(sim.energy_per_particle() < -1.0, "U/N = {}", sim.energy_per_particle());
    assert!(sim.chain_order() > 0.7, "order = {}", sim.chain_order());
    assert!(sim.min_sep() > 1.99);

    // a strong shear destroys alignment: compare energy at small vs large Mn
    let mut low = ErSim::new(n, lx, ly);
    for _ in 0..2500 {
        low.step(dt);
    }
    let mut high = low.clone();
    high.shear = 6.0 * PI * 2.0; // Mn = 2
    for _ in 0..2500 {
        high.step(dt);
    }
    assert!(high.energy_per_particle() > low.energy_per_particle() + 0.5);
    assert!(high.min_sep() > 1.99);
}
