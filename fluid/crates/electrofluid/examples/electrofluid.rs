//! Benchmarks for the electrofluid crate. Run:
//!     cargo run -p electrofluid --example electrofluid
//!
//! Prints the numbers quoted in `analysis/electrofluid.md`.

use electrofluid::eo::Channel;
use electrofluid::er::ErSim;
use electrofluid::leaky::LeakyDrop;
use electrofluid::pnp::*;
use std::f64::consts::PI;

fn main() {
    let t = T0;
    println!("=== 1. Debye length (water, 298 K) ===");
    for &c in &[1.0e-3, 1.0e-2, 1.0e-1] {
        let lam = debye_length(EPS_WATER, t, c, 1.0);
        println!("c = {:>6.1} mM : lambda_D = {:8.3} nm", c * 1e3, lam * 1e9);
    }

    println!("\n=== 1b. Gouy-Chapman, Debye-Huckel, Grahame ===");
    let c = 1e-3;
    let lam = debye_length(EPS_WATER, t, c, 1.0);
    let psi0 = -0.100; // 100 mV
    println!("lambda_D = {:.4} nm, psi0 = {:.0} mV", lam * 1e9, psi0 * 1e3);
    let (x, psi) = pb_solve(EPS_WATER, t, c, 1.0, psi0, 8.0 * lam, 4000);
    let mut max_gc = 0.0f64;
    let mut max_dh = 0.0f64;
    for (i, (&xi, &pi)) in x.iter().zip(psi.iter()).enumerate() {
        if i % 400 == 0 {
            let gc = gouy_chapman(xi, lam, psi0, t, 1.0);
            let dh = debye_huckel(xi, lam, psi0);
            max_gc = max_gc.max((pi - gc).abs());
            max_dh = max_dh.max((pi - dh).abs());
        }
    }
    println!(
        "max |PB_num - GouyChapman| = {:.3e} V ; max |PB_num - DebyeHuckel| = {:.3e} V",
        max_gc, max_dh
    );
    let sig = grahame_sigma(EPS_WATER, t, c, 1.0, psi0);
    let sig_num = sigma_from_profile(&psi, &x, EPS_WATER);
    println!(
        "Grahame sigma = {:.6e} C/m^2 ; from numerical profile = {:.6e} C/m^2 (rel {:.2e})",
        sig,
        sig_num,
        (sig - sig_num).abs() / sig.abs()
    );
    // Debye-Huckel limit: small potential
    let psi_small = -0.005;
    let (xs, psis) = pb_solve(EPS_WATER, t, c, 1.0, psi_small, 8.0 * lam, 2000);
    let mut max_dh_small = 0.0f64;
    for (i, (&xi, &pi)) in xs.iter().zip(psis.iter()).enumerate() {
        if i % 200 == 0 {
            max_dh_small = max_dh_small.max((pi - debye_huckel(xi, lam, psi_small)).abs());
        }
    }
    println!(
        "small psi0 = {:.0} mV: max |PB - DH| = {:.3e} V (rel {:.2e})",
        psi_small * 1e3,
        max_dh_small,
        max_dh_small / psi_small.abs()
    );

    println!("\n=== 2. Electro-osmotic slip fold (Helmholtz-Smoluchowski) ===");
    let e_field = 1.0e4;
    let zeta = -0.025;
    let h = 60e-9;
    for &c in &[1e-4, 1e-3, 1e-2, 1e-1] {
        let ch = Channel {
            eps_r: EPS_WATER,
            t,
            c_molar: c,
            z: 1.0,
            zeta,
            e_field,
            mu: MU_WATER,
            h,
        };
        let n = 600;
        let (uc, us) = ch.centerline(n);
        let lam = ch.lambda();
        println!(
            "c={:>6.1} mM  lam_D/h={:.4}  u_center={:+.6e}  u_HS={:+.6e}  rel err={:.3e}",
            c * 1e3,
            lam / h,
            uc,
            us,
            (uc - us).abs() / us.abs()
        );
    }

    println!("\n=== 3. Leaky-dielectric (Taylor-Melcher) drop ===");
    let cases = [
        ("water in oil (R>>1, S>1)", 6.0e-6, 2.0, 80.0, 1.0e-12, 1.0e-6, 2.0e-3, 5.0e-3, 4.0e5),
        ("R=1, more polarisable (S>1)", 1.0e-6, 2.0, 8.0, 1.0e-9, 1.0e-9, 1.0e-3, 1.0e-3, 1.0e6),
        ("R=1, less polarisable (S<1)", 1.0e-6, 2.0, 1.0, 1.0e-9, 1.0e-9, 1.0e-3, 1.0e-3, 1.0e6),
        ("R=0.1, S=1.5", 1.0e-6, 2.0, 3.0, 1.0e-9, 1.0e-10, 1.0e-3, 1.0e-3, 1.0e6),
    ];
    for (name, a, eo, ei, so, si, mo, mi, e0) in cases {
        let d = LeakyDrop {
            a,
            eps_o: eo,
            eps_i: ei,
            sigma_o: so,
            sigma_i: si,
            mu_o: mo,
            mu_i: mi,
            gamma: 4.0e-3,
            e0,
        };
        let (res, bc, dnum, _dana) = d.verify();
        println!(
            "{name}: R={:.3} S={:.3} Ca_E={:.3e} D={:+.4e} ({}) Phi={:+.4}",
            d.r_ratio(),
            d.s_ratio(),
            d.ca_e(),
            d.deformation(),
            if d.prolate() { "prolate" } else { "oblate" },
            d.discriminating()
        );
        println!(
            "    shape amplitude s={:+.4e} (= 4/3 D = {:+.4e}); vel jump/scale={res:.2e}, tang-stress mismatch/T={bc:.2e}",
            dnum,
            d.shape_amplitude(),
        );
        println!(
            "    circulation U_T (sin2θ coeff)={:+.4e} m/s",
            d.circulation_amplitude()
        );
    }
    println!(
        "Note: the P2 balance returns the shape amplitude s in r=a(1+s P2);\n\
         Taylor's D=(L-B)/(L+B) equals 3s/4 exactly, so s/D=4/3 is geometry,\n\
         not missing physics (Vlahovska 2019; Das & Saintillan 2020)."
    );

    // Published finite-Ca comparison. Sorgentone & Vlahovska (JFM 2021) BEM
    // for an isolated leaky-dielectric drop: R = σ_i/σ_o = 0.1, S = 1,
    // λ = μ_i/μ_o = 1, Ca = 0.1. Their figure 4 shows D ≈ -0.030; the
    // small-deformation Taylor value is
    let e0 = 1.0e6_f64;
    let a = 1.0e-3_f64;
    let eps_o = 1.0_f64;
    let gamma = eps_o * EPS0 * e0 * e0 * a / 0.1; // gives Ca_E = 0.1
    let cmp = LeakyDrop {
        a,
        eps_o,
        eps_i: 1.0,
        sigma_o: 1e-9,
        sigma_i: 1e-10,
        mu_o: 1e-3,
        mu_i: 1e-3,
        gamma,
        e0,
    };
    println!(
        "BEM comparison (Sorgentone & Vlahovska): R={:.2} S={:.2} lambda={:.2} Ca={:.2}: D_Taylor={:+.4} (BEM fig4 ~ -0.030)",
        cmp.r_ratio(),
        cmp.s_ratio(),
        cmp.mu_i / cmp.mu_o,
        cmp.ca_e(),
        cmp.deformation()
    );

    println!("\n=== 4. Electrorheology: chains, yield stress, Mason number ===");
    let phi = 0.311;
    let n = 48usize;
    let lx = 22.0;
    let ly = 22.0;
    let mut sim = ErSim::new(n, lx, ly);
    // form chains
    let dt = 1e-3;
    for _ in 0..4000 {
        sim.step(dt);
    }
    println!(
        "zero shear after 4000 steps: chain order = {:.3}, U/N = {:+.4}",
        sim.chain_order(),
        sim.energy_per_particle()
    );
    println!("Mason number Mn = shear/(6*pi); tau in units eps0 eps_f beta^2 E^2, phi = {phi}");
    for &mn in &[0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0] {
        let shear = mn * 6.0 * PI;
        let mut s = ErSim::new(n, lx, ly);
        // form chains at zero shear
        for _ in 0..4000 {
            s.step(dt);
        }
        s.shear = shear;
        for _ in 0..4000 {
            s.step(dt);
        }
        // average over a window
        let mut tau = 0.0;
        let mut ord = 0.0;
        let mut energy = 0.0;
        let navg = 600;
        for _ in 0..navg {
            s.step(dt);
            tau += s.shear_stress();
            ord += s.chain_order();
            energy += s.energy_per_particle();
        }
        println!(
            "Mn={:>5.2}  tau_E={:+.4}  U/N={:+.4}  chain_order={:.3}  min_sep={:.3}  tau_E/phi={:+.4}",
            mn,
            tau / navg as f64,
            energy / navg as f64,
            ord / navg as f64,
            s.min_sep(),
            tau / navg as f64 / phi
        );
    }
}
