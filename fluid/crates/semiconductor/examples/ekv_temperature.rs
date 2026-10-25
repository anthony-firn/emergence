//! Fold the drift–diffusion transistors' temperature dependence into the
//! EKV compact law used by the standard cells.
//!
//! Run: `cargo run -p semiconductor --release --example ekv_temperature`.
//!
//! For the 150 nm n- and pMOS of the scaled technology (full mobility
//! model, as in `analysis/semiconductor.md` §13), at 250–425 K, the TCAD
//! on-current (|V_gs| = |V_ds| = 1.2 V) and off-current (V_gs = 0) are
//! sampled. The plain EKV law cannot follow the TCAD curve's shape (no
//! velocity saturation or mobility degradation: a direct fit misses I_on by
//! 20 %), so what is folded is what the cells depend on: the *ratios*
//! I_on(T)/I_on(300 K) and I_off(T)/I_off(300 K), through two physical
//! parameters per polarity, the threshold's temperature slope and the
//! exponent of I_s ∝ (T/300)^a (mobility and V_T² together). They are fitted
//! on 250, 300, 350 and 400 K; 275, 325, 375 and 425 K are held out. The
//! printed coefficients are what `gates::Ekv::at_temperature` uses.

use semiconductor::compact::sample;
use semiconductor::mosfet::MosfetSpec;
use semiconductor::silicon::Silicon;

const VDD: f64 = 1.2;
const KB_Q: f64 = 8.617333262e-5; // Boltzmann constant / e, V/K (CODATA 2018)

fn phi(x: f64) -> f64 {
    (1.0 + (0.5 * x).exp()).ln().powi(2)
}

fn ekv(p: [f64; 3], vgs: f64, t: f64) -> f64 {
    let vt = KB_Q * t;
    let vp = (vgs - p[1]) / (p[2] * vt);
    p[0] * (phi(vp) - phi(vp - VDD / vt))
}

#[allow(dead_code)]
/// Nelder–Mead on (ln I_s, V_th, n).
fn fit(vgs: &[f64], id: &[f64], t: f64) -> ([f64; 3], f64) {
    let cost = |q: &[f64; 3]| -> f64 {
        let p = [q[0].exp(), q[1], q[2]];
        vgs.iter().zip(id).map(|(&v, &i)| (ekv(p, v, t).ln() - i.ln()).powi(2)).sum::<f64>()
    };
    let mut s: Vec<[f64; 3]> = vec![[(-12.0f64), 0.3, 1.3], [-11.0, 0.3, 1.3], [-12.0, 0.4, 1.3], [-12.0, 0.3, 1.5]];
    for _ in 0..4000 {
        s.sort_by(|a, b| cost(a).total_cmp(&cost(b)));
        let c: [f64; 3] = std::array::from_fn(|k| (s[0][k] + s[1][k] + s[2][k]) / 3.0);
        let w = s[3];
        let r: [f64; 3] = std::array::from_fn(|k| c[k] + (c[k] - w[k]));
        if cost(&r) < cost(&s[0]) {
            let e: [f64; 3] = std::array::from_fn(|k| c[k] + 2.0 * (c[k] - w[k]));
            s[3] = if cost(&e) < cost(&r) { e } else { r };
        } else if cost(&r) < cost(&s[2]) {
            s[3] = r;
        } else {
            let k: [f64; 3] = std::array::from_fn(|k| c[k] + 0.5 * (w[k] - c[k]));
            if cost(&k) < cost(&w) {
                s[3] = k;
            } else {
                let b = s[0];
                for v in s.iter_mut().skip(1) {
                    *v = std::array::from_fn(|k| b[k] + 0.5 * (v[k] - b[k]));
                }
            }
        }
    }
    s.sort_by(|a, b| cost(a).total_cmp(&cost(b)));
    let q = s[0];
    let rms = (cost(&q) / vgs.len() as f64).sqrt();
    ([q[0].exp(), q[1], q[2]], rms)
}

fn main() {
    let temps = [250.0, 275.0, 300.0, 325.0, 350.0, 375.0, 400.0, 425.0];
    let train = [0usize, 2, 4, 6];
    for pmos in [false, true] {
        let pol = if pmos { -1.0 } else { 1.0 };
        let name = if pmos { "pMOS" } else { "nMOS" };
        let spec = MosfetSpec::scaled(150e-9, pmos);
        let mut dev = spec.build(Silicon::new(300.0));
        dev.velocity_saturation = true;
        dev.lombardi = true;
        dev.screening = true;
        dev.solve().expect("equilibrium");
        let mut on = Vec::new();
        let mut off = Vec::new();
        // Walk up from 300 K, then back down through every point (the
        // Newton solve needs small temperature steps) to 275 and 250 K.
        let mut got = vec![(0.0, 0.0); temps.len()];
        let path = [2usize, 3, 4, 5, 6, 7, 6, 5, 4, 3, 2, 1, 0];
        for &i in &path {
            dev.set_temperature(temps[i]);
            dev.solve().expect("temperature step");
            if got[i].1 == 0.0 {
                let id = sample(&mut dev, &[0.0, pol * VDD], &[pol * VDD]);
                got[i] = ((id[0][0] * pol).abs(), (id[1][0] * pol).abs());
            }
        }
        for (i, &t) in temps.iter().enumerate() {
            off.push(got[i].0);
            on.push(got[i].1);
            println!("{name} T {t:.0} K: I_off {:.4e} A/µm, I_on {:.5e} A/µm", got[i].0, got[i].1);
        }
        // The cell-level law at 300 K (gates::Ekv::n / p) and its ratios.
        let (vth0, n, ion300) = (0.32, 1.3, if pmos { 230e-6 } else { 330e-6 });
        let vt0 = KB_Q * 300.0;
        let p0 = [ion300, vth0, n];
        let isat0 = {
            let vp = (VDD - vth0) / (n * vt0);
            ion300 / (phi(vp) - phi(vp - VDD / vt0))
        };
        let law = |dv: f64, a: f64, t: f64, vgs: f64| -> f64 {
            ekv([isat0 * (t / 300.0f64).powf(a), vth0 + dv * (t - 300.0), n], vgs, t)
        };
        let _ = p0;
        let ratio_err = |dv: f64, a: f64, i: usize| -> (f64, f64) {
            let r_on = law(dv, a, temps[i], VDD) / law(dv, a, 300.0, VDD);
            let r_off = law(dv, a, temps[i], 0.0) / law(dv, a, 300.0, 0.0);
            (r_on / (on[i] / on[2]), r_off / (off[i] / off[2]))
        };
        let mut best = (f64::INFINITY, 0.0, 0.0);
        for k in 0..=2000 {
            let dv = -2e-3 + 1e-6 * k as f64;
            for j in 0..=300 {
                let a = -1.5 + 0.01 * j as f64;
                let c: f64 = train.iter().map(|&i| {
                    let (x, y) = ratio_err(dv, a, i);
                    x.ln().powi(2) + y.ln().powi(2)
                }).sum();
                if c < best.0 {
                    best = (c, dv, a);
                }
            }
        }
        let (_, dv, a) = best;
        println!("{name}: dV_th/dT = {dv:.4e} V/K, I_s ∝ (T/300)^{a:.2} (fit on 250/300/350/400 K)");
        for (i, &t) in temps.iter().enumerate() {
            let (x, y) = ratio_err(dv, a, i);
            let tag = if train.contains(&i) { "fit " } else { "held" };
            println!("  {tag} T {t:.0} K: I_on ratio model/TCAD {x:.4}, I_off ratio {y:.4}");
        }
    }
}
