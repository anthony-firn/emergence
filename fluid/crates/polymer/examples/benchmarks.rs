//! Prints the benchmark numbers quoted in `analysis/polymer.md`.
//!
//! Run with `cargo run -p polymer --release --example benchmarks`.

use polymer::continuum::{Maxwell, OldroydB, PronySeries};
use polymer::fene::{Dumbbell, fene_p_steady, hookean_steady_cxx};
use polymer::gpr_fold::{MaxwellMode, maxwell_to_gpr};
use polymer::indicator::shear_indicator;
use polymer::rheology::{Rheometer, Wlf, polystyrene_melt, shift_prony};
use polymer::rouse::{Rouse, prony_complex};
use polymer::zimm::Zimm;

fn main() {
    println!("=== Rouse, N = 50 ===");
    let r = Rouse::new(50, 1.0, 1.0, 1.0, 1.0);
    let n = r.n_beads as f64;
    println!(
        "tau1 exact {:.6}, continuum {:.6}, ratio {:.4}",
        r.tau(1),
        r.tau_continuum(1),
        r.tau(1) / r.tau_continuum(1)
    );
    println!(
        "eta0 exact {:.6}, Rouse limit n zeta N^2 b^2/36 {:.6}, ratio {:.4}",
        r.eta0(),
        n * n / 36.0,
        r.eta0() / (n * n / 36.0)
    );
    println!("G(0) = n kT (N-1) = {}", r.g_modulus(0.0));

    // Fold onto m slow modes + one lumped fast mode.
    let full = r.prony();
    let tau1 = r.tau(1);
    for m in [1usize, 2, 3, 5] {
        let fold = r.prony_lumped(m);
        // terminal viscosity
        let eta_full: f64 = full.iter().map(|(g, l)| g * l).sum();
        let eta_fold: f64 = fold.iter().map(|(g, l)| g * l).sum();
        // G(t) relative error over [0.1 tau1, 10 tau1]
        let mut max_gt = 0.0f64;
        for k in 0..200 {
            let t = tau1 * 10f64.powf(-1.0 + 2.0 * k as f64 / 199.0);
            let a = r.g_modulus(t);
            let b = polymer::rouse::prony_g(&fold, t);
            max_gt = max_gt.max((b - a).abs() / a);
        }
        // G'' relative error over omega tau1 in [1e-3, 1e3] and in the
        // terminal region [1e-3, 0.1].
        let mut max_gpp = 0.0f64;
        let mut max_gpp_term = 0.0f64;
        for k in 0..200 {
            let w = 10f64.powf(-3.0 + 6.0 * k as f64 / 199.0) / tau1;
            let (_, a) = r.complex(w);
            let (_, b) = prony_complex(&fold, w);
            max_gpp = max_gpp.max((b - a).abs() / a);
        }
        for k in 0..200 {
            let w = 10f64.powf(-3.0 - 1.0 * k as f64 / 199.0) / tau1;
            let (_, a) = r.complex(w);
            let (_, b) = prony_complex(&fold, w);
            max_gpp_term = max_gpp_term.max((b - a).abs() / a);
        }
        println!(
            "m={m}: modes={} eta={:.6} (rel {:.3e}), max|dG|/G over [0.1,10]tau1={:.3e}, max|dG''|/G'' over w tau1 in [1e-3,1e3]={:.3e}, terminal={:.3e}",
            fold.len(),
            eta_fold,
            (eta_fold - eta_full).abs() / eta_full,
            max_gt,
            max_gpp,
            max_gpp_term
        );
    }
    // high-frequency Rouse slope
    let big = Rouse::new(200, 1.0, 1.0, 1.0, 1.0);
    let w1 = 10.0 / big.tau1();
    let (_, g1) = big.complex(w1);
    let (_, g2) = big.complex(4.0 * w1);
    println!(
        "Rouse high-frequency G'' slope d ln G''/d ln w = {:.4} (expect 0.5)",
        (g2 / g1).ln() / 4f64.ln()
    );

    println!("\n=== Brownian dynamics Green-Kubo, N = 8 ===");
    let rb = Rouse::new(8, 1.0, 1.0, 1.0, 1.0);
    let res = polymer::bd::green_kubo(&rb, 200, 40_000, 4_000, 0.005, 10, 120, 12345);
    let analytic = |t: f64| rb.g_modulus(t) / rb.n_chains / rb.kbt;
    let mut max_rel = 0.0f64;
    let mut max_t = 0.0;
    for (i, &t) in res.times.iter().enumerate() {
        let a = analytic(t);
        if a > 0.3 {
            let e = (res.g[i] - a).abs() / a;
            if e > max_rel {
                max_rel = e;
                max_t = t;
            }
        }
    }
    println!("G(0) BD {:.4} vs {} (rel {:.3})", res.g[0], rb.n_beads - 1, (res.g[0] - 7.0).abs() / 7.0);
    println!("max relative error in G(t) for analytic>0.3: {max_rel:.4} at t={max_t:.3} (tau1={:.3})", rb.tau1());
    let tau1b = rb.tau1();
    let i1 = res.times.iter().position(|&t| t >= tau1b).unwrap();
    let i2 = res.times.iter().position(|&t| t >= 1.6 * tau1b).unwrap();
    let slope = (res.ee[i1].ln() - res.ee[i2].ln()) / (res.times[i2] - res.times[i1]);
    println!("end-to-end ee(0) {:.4}, long-time slope {:.4} vs 1/tau1 {:.4}", res.ee[0], slope, 1.0 / tau1b);

    println!("\n=== Zimm ===");
    for (a, b) in [(20usize, 40usize), (40, 80), (80, 160)] {
        let ta = Zimm { n_beads: a, b: 1.0, zeta: 1.0, kbt: 1.0, h_star: 0.25 }.tau1();
        let tb = Zimm { n_beads: b, b: 1.0, zeta: 1.0, kbt: 1.0, h_star: 0.25 }.tau1();
        println!("h*=0.25 tau1({a})={ta:.4} tau1({b})={tb:.4} exponent {:.4}", (tb / ta).ln() / 2f64.ln());
    }
    let r0a = Zimm { n_beads: 40, b: 1.0, zeta: 1.0, kbt: 1.0, h_star: 0.0 }.tau1();
    let r0b = Zimm { n_beads: 80, b: 1.0, zeta: 1.0, kbt: 1.0, h_star: 0.0 }.tau1();
    println!("h*=0 (Rouse) 40->80 exponent {:.4}", (r0b / r0a).ln() / 2f64.ln());

    println!("\n=== Continuum benchmarks ===");
    let ob = OldroydB { eta_s: 0.1, eta_p: 0.9, lambda: 2.0 };
    let gdot = 0.3;
    println!("Oldroyd-B N1 steady {:.6}, 2 eta_p lambda gdot^2 {:.6}", ob.n1_steady(gdot), 2.0 * ob.eta_p * ob.lambda * gdot * gdot);
    println!("sigma12 steady {:.6}, eta0 gdot {:.6}", ob.sigma12_steady(gdot), ob.eta0() * gdot);
    let m = Maxwell { g: 2.0, lambda: 0.5 };
    let (gp, gpp) = m.complex(1e-5);
    println!("Maxwell terminal G'/w^2 {:.6}, G''/w {:.6} (eta0 {:.6})", gp / 1e-10, gpp / 1e-5, m.eta0());

    println!("\n=== GPR fold ===");
    let mm = MaxwellMode { g: 2.5e5, lambda: 1e-3 };
    let r = maxwell_to_gpr(mm, 1000.0);
    println!(
        "Maxwell G={} lambda={} -> rho c_s^2={} c_s={:.6} tau_GPR={} ; mu={} eta={}",
        mm.g, mm.lambda, r.shear_modulus(), r.cs, r.tau, r.viscosity(), mm.eta()
    );

    println!("\n=== FENE ===");
    println!(
        "coil-stretch Wi = 0.5; Hookean cxx(0.3) = {:.3}, cxx(0.49) = {:.3}",
        hookean_steady_cxx(0.3),
        hookean_steady_cxx(0.49)
    );
    let (cx, cy, tr) = fene_p_steady(0.8, 50.0);
    println!("FENE-P b=50 Wi=0.8: cxx={:.4} cyy={:.4} tr/Aeq={:.4} (< b=50)", cx, cy, tr);
    let d = Dumbbell::fene(50.0);
    let (bd_cxx, q2) = d.run_uniaxial(0.8, 30, 40_000, 20_000, 0.002, 11);
    println!("FENE BD b=50 Wi=0.8: cxx={bd_cxx:.4}, <Q^2>={q2:.4}");

    println!("\n=== WF / TTS / DMA ===");
    let wlf = Wlf::universal(373.15);
    println!("WLF log10 a at Tg+50 = {:.5} (a={:.3e})", wlf.log10_a(423.15), wlf.a(423.15));
    let base = polystyrene_melt();
    let mut max_tts = 0.0f64;
    for &t in &[400.0, 430.0, 470.0] {
        let a = wlf.a(t);
        let shifted = PronySeries::new(shift_prony(&base.modes, a));
        for &w in &[0.1, 1.0, 10.0, 100.0] {
            let (g1, _) = shifted.complex(w);
            let (g2, _) = base.complex(a * w);
            max_tts = max_tts.max((g1 - g2).abs() / g2);
        }
    }
    println!("TTS collapse max relative error = {max_tts:.3e}");
    let rh = Rheometer::new(base.clone(), 0.0);
    let mut max_dma = 0.0f64;
    for &w in &[0.1, 1.0, 10.0] {
        let (gm, gpm) = rh.oscillatory_from_ode(1e-3, w, 0.002, 20);
        let (ga, gpa) = rh.complex(w);
        max_dma = max_dma.max((gm - ga).abs() / ga).max((gpm - gpa).abs() / gpa);
    }
    println!("virtual DMA max relative error = {max_dma:.3e}");
    println!("PS melt G_N0 = {} Pa, eta0 = {:.4} Pa s (estimate)", base.g_t(0.0), base.eta0());

    println!("\n=== Indicator ===");
    for gamma in [1e-3, 1.0, 1e3] {
        let ind = shear_indicator(1.0, gamma);
        println!("lambda=1 gamma={gamma}: De={:.3} Wi={:.3} regime={:?}", ind.de, ind.wi, ind.regime);
    }
}
