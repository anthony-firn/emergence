//! Quasi-static MHD against exact solutions.

use spectral_core::mhd::Mhd;
use spectral_core::{init, Grid, Solver};

fn noise(len: usize, seed: u64) -> Vec<f64> {
    let mut s = seed;
    (0..len)
        .map(|_| {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            (s as f64 / u64::MAX as f64) - 0.5
        })
        .collect()
}

/// Moffatt (1967): in a uniform field each Fourier mode of a small
/// disturbance decays at exactly ν|k|² + (σ/ρ)(k̂·B)². Checked mode by mode,
/// with the Joule term explicit and with it inside the integrating factor.
#[test]
fn joule_damping_is_anisotropic_and_exact() {
    let g = Grid::new(16);
    let amp = 1e-9; // nonlinearity negligible
    let u0: [Vec<f64>; 3] = std::array::from_fn(|c| noise(g.real_len(), 0xABCDEF + c as u64).iter().map(|v| v * amp).collect());
    let (nu, sor) = (0.02, 3.0);
    let b = [0.3, -0.5, 0.8];
    let t_end = 0.4;
    for exact_if in [false, true] {
        let mut s = Solver::new(g.clone(), nu, u0.clone());
        let mut m = Mhd::uniform(sor, 1.0, b);
        m.joule_in_if = exact_if;
        let mut s0 = s.u.clone();
        s.mhd = Some(m.clone());
        let h: f64 = 0.01;
        let steps = (t_end / h).round() as usize;
        for _ in 0..steps {
            s.step(h);
        }
        let t = s.t;
        let (mut err, mut scale) = (0.0f64, 0.0f64);
        for c in 0..3 {
            for i in 0..g.spec_len() {
                let k = g.k(i);
                let k2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
                s0[c][i] *= (-(nu * k2 + m.joule_rate(k)) * t).exp();
                err = err.max((s.u[c][i] - s0[c][i]).norm());
                scale = scale.max(s0[c][i].norm());
            }
        }
        let rel = err / scale;
        println!("exact_if={exact_if}: max relative mode error {rel:.2e}");
        assert!(rel < 1e-6, "rel {rel:e}");
    }
}

/// Electrode-driven flow: current source s = S cos x in a field B ẑ drives
/// J = S sin x x̂ and the force −S B sin x ŷ /ρ. The motional EMF of the
/// resulting u_y(x) ŷ is a pure gradient (no current), and the nonlinear
/// term vanishes, so u_y = −(S B/ρν)(1 − e^{−νt}) sin x exactly.
#[test]
fn electrode_current_drives_exact_flow() {
    let g = Grid::new(16);
    let (nu, rho, sigma, bz, s_amp) = (0.1, 2.0, 5.0, 0.7, 1.3);
    let mut s = Solver::new(g.clone(), nu, init::from_fn(&g, |_, _, _| [0.0; 3]));
    let mut m = Mhd::uniform(sigma / rho, 1.0 / rho, [0.0, 0.0, bz]);
    m.source = Some(init::from_fn(&g, |x, _, _| [s_amp * x.cos(), 0.0, 0.0])[0].clone());
    s.mhd = Some(m);
    let h = 0.02;
    for _ in 0..250 {
        s.step(h);
    }
    let t = s.t;
    let uy = {
        let fft = spectral_core::fft::Fft3::new(&g);
        let mut w = s.u[1].clone();
        let mut out = vec![0.0; g.real_len()];
        fft.inverse_inplace(&mut w, &mut out);
        out
    };
    let n = g.n;
    let mut err = 0.0f64;
    for p in 0..g.real_len() {
        let i = p / (n * n);
        let exact = -(s_amp * bz / (rho * nu)) * (1.0 - (-nu * t).exp()) * g.x(i).sin();
        err = err.max((uy[p] - exact).abs());
    }
    let scale = s_amp * bz / (rho * nu) * (1.0 - (-nu * t).exp());
    println!("electrode flow: max error {:.2e} relative", err / scale);
    assert!(err / scale < 1e-9);
}
