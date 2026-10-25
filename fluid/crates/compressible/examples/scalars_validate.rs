//! Transported scalars (`scalars.rs`).
//!
//!   scalars_validate
//!
//! 1. Gaussian advected one period in uniform flow: exact conservation,
//!    boundedness, TVD smearing error.
//! 2. φ = 1 carried by the liquid of a water slab moving through air
//!    (five-equation model), explicit and all-speed (acoustic Courant 64):
//!    φ must stay exactly 1 wherever there is liquid.
//! 3. Diffusion of a Gaussian at rest vs σ² = σ₀² + 2Dt.
use compressible::allspeed::AllSpeedOpts;
use compressible::eos::materials::{AIR, WATER_SG};
use compressible::eos::Nasg;
use compressible::multiphase::FiveEq;
use compressible::scalars::{Carrier, Scalar};
use compressible::{Euler, Mesh, Solver};

fn main() {
    let eos = Nasg::ideal(1.4);
    let n = 200;
    // 1. advection
    {
        let mut s = Solver::new(Euler { eos }, Mesh::new([n, 1, 1], [0.0; 3], [1.0, 1.0, 1.0]), |_| [1.0, 1.0, 0.0, 0.0, 1.0]);
        s.cfl = 0.4;
        let g = |x: f64| (-((x - 0.5) / 0.05).powi(2)).exp();
        let sc = Scalar::new(&s, "dye", Carrier::Mixture, 0.0, |x| g(x[0]));
        s.scalars.push(sc);
        let m0 = s.scalar_total(0);
        s.run_to(1.0);
        let (mut mn, mut mx, mut err) = (f64::INFINITY, f64::NEG_INFINITY, 0.0);
        for i in 0..n {
            let phi = s.scalars[0].phi(&s, i, 0, 0);
            mn = mn.min(phi);
            mx = mx.max(phi);
            err += (phi - g(s.mesh.x(i, 0, 0)[0])).abs() / n as f64;
        }
        println!("advection, one period: total {:+.2e} relative change, φ in [{mn:.3e}, {mx:.4}], L1 error {err:.3e} (TVD, MC limiter)", s.scalar_total(0) / m0 - 1.0);
    }
    // 2. liquid-carried φ = 1 through a moving interface
    for allspeed in [false, true] {
        let sys = FiveEq { eos: [WATER_SG, AIR], kapila: true };
        let slab = |x: [f64; 3]| {
            let al = 1e-6 + (1.0 - 2e-6) * 0.5 * (((x[0] - 0.3) / 0.01).tanh() - ((x[0] - 0.7) / 0.01).tanh());
            [al * 1000.0, (1.0 - al) * 1.19, 10.0, 0.0, 0.0, 1e5, al]
        };
        let mut s = Solver::new(sys, Mesh::new([n, 1, 1], [0.0; 3], [1.0, 1.0, 1.0]), slab);
        s.cfl = 0.4;
        if allspeed {
            s.allspeed = Some(AllSpeedOpts { min_gain: 0.0, ..Default::default() });
        }
        let sc = Scalar::new(&s, "dissolved", Carrier::Phase(0), 0.0, |_| 1.0);
        s.scalars.push(sc);
        let m0 = s.scalar_total(0);
        let steps = s.run_to(0.02);
        let mut dev: f64 = 0.0;
        for i in 0..n {
            dev = dev.max((s.scalars[0].phi(&s, i, 0, 0) - 1.0).abs());
        }
        println!(
            "liquid-carried φ = 1, slab at 10 m/s, {} ({steps} steps): max |φ − 1| {dev:.2e}, total {:+.2e}",
            if allspeed { "all-speed" } else { "explicit " },
            s.scalar_total(0) / m0 - 1.0
        );
    }
    // 3. diffusion at rest
    {
        let mut s = Solver::new(Euler { eos }, Mesh::new([n, 1, 1], [0.0; 3], [1.0, 1.0, 1.0]), |_| [1.0, 0.0, 0.0, 0.0, 1.0]);
        let (sig0, d, t) = (0.03f64, 1e-3, 0.5);
        let sc = Scalar::new(&s, "heat-like", Carrier::Mixture, d, |x| (-((x[0] - 0.5) / sig0).powi(2) / 2.0).exp());
        s.scalars.push(sc);
        s.allspeed = Some(AllSpeedOpts { min_gain: 0.0, ..Default::default() });
        let m0 = s.scalar_total(0);
        let steps = s.run_to(t);
        let sig = (sig0 * sig0 + 2.0 * d * t).sqrt();
        let mut err: f64 = 0.0;
        for i in 0..n {
            let x = s.mesh.x(i, 0, 0)[0];
            let ex = sig0 / sig * (-((x - 0.5) / sig).powi(2) / 2.0).exp();
            err = err.max((s.scalars[0].phi(&s, i, 0, 0) - ex).abs());
        }
        println!("diffusion at rest (all-speed, {steps} steps): max |φ − exact| {err:.2e} (peak {:.4}), total {:+.2e}", sig0 / sig, s.scalar_total(0) / m0 - 1.0);
    }
}
