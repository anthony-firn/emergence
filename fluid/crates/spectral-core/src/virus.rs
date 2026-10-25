//! Virus-particle (virion) transport in a fluid: advection by the simulated
//! flow + Brownian diffusion + gravitational settling, in the overdamped
//! (Stokes number → 0) limit.
//!
//! A virion is 20–300 nm across, far below any fluid grid cell, so it is a
//! point particle. Its momentum relaxation time τ = (2/9)ρr²C_c/η is ~1e-7 s
//! in air and less in water, versus fluid timescales ≥ 1e-3 s, so inertia is
//! negligible and the exact model is the overdamped Langevin equation
//!
//!   dx = u(x, t) dt − v_s ẑ dt + √(2D) dW,
//!
//! with Stokes–Einstein diffusivity D = k_B T C_c / (6π η r) and Stokes
//! settling speed v_s = (2/9)(ρ_p − ρ_f) g r² C_c / η. C_c is the Cunningham
//! slip correction (≈ 2.5 for a 120 nm virion in air; 1 in water). With
//! additive noise and constant D, Euler–Maruyama integrates the Brownian part
//! exactly in distribution; the drift part is exact to O(dt²) locally, and
//! the flow time step is already CFL-limited well below that.
//!
//! Velocities are sampled from the physical N³ grid by trilinear
//! interpolation (exact for resolved smooth fields; the virion is always a
//! sub-grid point). Particles are independent, so stepping is parallel over
//! particles with one deterministic xorshift stream per particle.
//!
//! Exact checks (see tests): quiescent fluid ⇒ ⟨|Δx|²⟩ = 6Dt (Einstein);
//! uniform flow ⇒ drift = U t; still fluid ⇒ ⟨z⟩ = −v_s t.

use rayon::prelude::*;

/// Boltzmann constant, J/K (exact SI value).
pub const K_B: f64 = 1.380_649e-23;
/// Standard gravity, m/s² (exact SI value).
pub const GRAVITY: f64 = 9.806_65;

const TWO_PI: f64 = 2.0 * std::f64::consts::PI;

/// Suspending medium properties at a fixed temperature.
#[derive(Clone, Debug)]
pub struct Medium {
    pub name: &'static str,
    /// Temperature, K.
    pub temperature: f64,
    /// Dynamic viscosity, Pa·s.
    pub eta: f64,
    /// Density, kg/m³.
    pub rho: f64,
    /// Gas mean free path, m; 0 means continuum (no slip correction).
    pub mean_free_path: f64,
}

impl Medium {
    /// Dry air at 20 °C, 1 atm.
    pub fn air() -> Self {
        Self {
            name: "air",
            temperature: 293.15,
            eta: 1.81e-5,
            rho: 1.204,
            mean_free_path: 66e-9,
        }
    }

    /// Water at 20 °C.
    pub fn water() -> Self {
        Self {
            name: "water",
            temperature: 293.15,
            eta: 1.002e-3,
            rho: 998.2,
            mean_free_path: 0.0,
        }
    }

    /// Kinematic viscosity, m²/s.
    pub fn nu(&self) -> f64 {
        self.eta / self.rho
    }
}

/// A spherical virus particle ("simple virus": no shape, no surface chemistry).
#[derive(Clone, Debug)]
pub struct Virion {
    /// Capsid radius, m.
    pub radius: f64,
    /// Mean density, kg/m³ (virions are ≈ 1.1 g/cm³: protein + RNA + water).
    pub rho: f64,
}

impl Virion {
    /// Influenza-like virion: 120 nm diameter.
    pub fn influenza() -> Self {
        Self { radius: 60e-9, rho: 1100.0 }
    }

    /// Cunningham slip correction C_c = 1 + (λ/r)(1.257 + 0.4 e^{−1.1 r/λ}).
    /// 1 for a continuum medium (water); ≈ 2.5 for a 120 nm virion in air.
    pub fn cunningham(&self, m: &Medium) -> f64 {
        if m.mean_free_path == 0.0 {
            return 1.0;
        }
        let l_over_r = m.mean_free_path / self.radius;
        1.0 + l_over_r * (1.257 + 0.4 * (-1.1 / l_over_r).exp())
    }

    /// Stokes–Einstein diffusivity D = k_B T C_c / (6π η r), m²/s.
    pub fn diffusivity(&self, m: &Medium) -> f64 {
        K_B * m.temperature * self.cunningham(m) / (6.0 * std::f64::consts::PI * m.eta * self.radius)
    }

    /// Stokes settling speed v_s = (2/9)(ρ_p − ρ_f) g r² C_c / η, m/s.
    /// ≈ 1.2 µm/s for a 120 nm virion in air; ≈ 0.8 nm/s in water.
    pub fn settling_speed(&self, m: &Medium) -> f64 {
        2.0 / 9.0 * (self.rho - m.rho) * GRAVITY * self.radius * self.radius * self.cunningham(m)
            / m.eta
    }

    /// Momentum relaxation time τ = v_s / g (per unit mass, Stokes regime).
    /// ≈ 1.2e-7 s in air: inertia is negligible on any resolved fluid
    /// timescale, which is what makes the overdamped model exact here.
    pub fn relaxation_time(&self, m: &Medium) -> f64 {
        2.0 / 9.0 * self.rho * self.radius * self.radius * self.cunningham(m) / m.eta
    }
}

/// Map between the solver's periodic box [0, 2π)³ and SI units: the box side
/// is `box_len` metres and unit solver velocity is `velocity` m/s.
#[derive(Clone, Copy, Debug)]
pub struct Scaling {
    pub box_len: f64,
    pub velocity: f64,
}

impl Scaling {
    /// Metres per unit solver length.
    pub fn l0(&self) -> f64 {
        self.box_len / TWO_PI
    }

    /// Seconds per unit solver time.
    pub fn t0(&self) -> f64 {
        self.l0() / self.velocity
    }

    /// Medium viscosity in solver units (ν̃ = ν t0 / l0²).
    pub fn nu_solver(&self, m: &Medium) -> f64 {
        m.nu() * self.t0() / (self.l0() * self.l0())
    }

    /// Diffusivity in solver units (D̃ = D t0 / l0² = 1/Pe).
    pub fn diffusivity_solver(&self, d_si: f64) -> f64 {
        d_si * self.t0() / (self.l0() * self.l0())
    }

    /// Settling speed in solver units (ṽ_s = v_s / U).
    pub fn settling_solver(&self, vs_si: f64) -> f64 {
        vs_si / self.velocity
    }
}

/// Deterministic per-particle xorshift64 stream (same schedule as the
/// solver's grid-shift RNG).
struct Rng(u64);

impl Rng {
    fn new(seed: u64, index: usize) -> Self {
        let mut s = Self(seed ^ (index as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1);
        for _ in 0..8 {
            s.next_u64(); // warm up: decorrelate neighbouring particles
        }
        s
    }

    #[inline]
    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }

    /// Uniform on (0, 1].
    #[inline]
    fn uniform(&mut self) -> f64 {
        1.0 - (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Standard normal, Box–Muller.
    #[inline]
    fn gauss(&mut self) -> f64 {
        let u1 = self.uniform();
        let u2 = self.uniform();
        (-2.0 * u1.ln()).sqrt() * (TWO_PI * u2).cos()
    }
}

/// Trilinear interpolation of a component field (layout `(i*N + j)*N + l`,
/// nodes at x_i = i·dx, periodic) at a point in [0, 2π)³.
pub fn trilinear(u: &[Vec<f64>; 3], n: usize, p: [f64; 3]) -> [f64; 3] {
    let dx = TWO_PI / n as f64;
    let mut out = [0.0; 3];
    let mut cell = [0usize; 3];
    let mut frac = [0.0; 3];
    for c in 0..3 {
        let s = p[c] / dx;
        let i = s.floor() as isize;
        frac[c] = s - i as f64;
        cell[c] = i.rem_euclid(n as isize) as usize;
    }
    let next = [
        (cell[0] + 1) % n,
        (cell[1] + 1) % n,
        (cell[2] + 1) % n,
    ];
    for c in 0..3 {
        let f = &u[c];
        let mut v = 0.0;
        for a in 0..2 {
            let i = if a == 0 { cell[0] } else { next[0] };
            let wa = if a == 0 { 1.0 - frac[0] } else { frac[0] };
            for b in 0..2 {
                let j = if b == 0 { cell[1] } else { next[1] };
                let wb = if b == 0 { 1.0 - frac[1] } else { frac[1] };
                for d in 0..2 {
                    let l = if d == 0 { cell[2] } else { next[2] };
                    let wd = if d == 0 { 1.0 - frac[2] } else { frac[2] };
                    v += wa * wb * wd * f[(i * n + j) * n + l];
                }
            }
        }
        out[c] = v;
    }
    out
}

/// A cloud of identical virions advected/diffusing in the periodic box.
pub struct Cloud {
    /// Wrapped positions in [0, 2π)³ (solver units) — used for flow sampling.
    pub pos: Vec<[f64; 3]>,
    /// Unwrapped displacements since t = 0 (solver units) — exact statistics
    /// even when particles cross the periodic boundary.
    pub disp: Vec<[f64; 3]>,
    /// Diffusivity in solver units.
    pub d: f64,
    /// Settling speed in solver units (gravity along −z).
    pub vs: f64,
    rng: Vec<Rng>,
}

/// Ensemble statistics in solver units.
#[derive(Clone, Copy, Debug)]
pub struct Stats {
    /// Mean displacement ⟨Δx⟩ per axis (drift: U t in x, −v_s t in z).
    pub mean: [f64; 3],
    /// ⟨|Δx − ⟨Δx⟩|²⟩: the Einstein term, = 6 D t exactly in still/uniform fluid.
    pub msd_about_mean: f64,
    /// ⟨|Δx|²⟩ = 6 D t + |drift|² t² in still/uniform fluid.
    pub msd: f64,
}

impl Cloud {
    /// `np` virions uniformly distributed over the box; `d` and `vs` in
    /// solver units (see [`Scaling`]).
    pub fn new(np: usize, d: f64, vs: f64, seed: u64) -> Self {
        let mut c = Self {
            pos: vec![[0.0; 3]; np],
            disp: vec![[0.0; 3]; np],
            d,
            vs,
            rng: (0..np).map(|i| Rng::new(seed, i)).collect(),
        };
        for i in 0..np {
            for k in 0..3 {
                c.pos[i][k] = c.rng[i].uniform() * TWO_PI;
            }
        }
        c
    }

    /// One Euler–Maruyama step: Δx = (u(x) − v_s ẑ) dt + √(2 D dt) ξ.
    /// `u` is the physical velocity on the `n`³ grid, or `None` for a
    /// quiescent fluid. Parallel over particles; deterministic for fixed seed.
    pub fn step(&mut self, u: Option<&[Vec<f64>; 3]>, n: usize, dt: f64) {
        let (d, vs) = (self.d, self.vs);
        let kick = (2.0 * d * dt).sqrt();
        self.pos
            .par_iter_mut()
            .zip(self.disp.par_iter_mut())
            .zip(self.rng.par_iter_mut())
            .for_each(|((p, disp), rng)| {
                let mut vel = match u {
                    Some(u) => trilinear(u, n, *p),
                    None => [0.0; 3],
                };
                vel[2] -= vs;
                for c in 0..3 {
                    let delta = vel[c] * dt + kick * rng.gauss();
                    disp[c] += delta;
                    p[c] = (p[c] + delta).rem_euclid(TWO_PI);
                }
            });
    }

    /// Ensemble displacement statistics (solver units).
    pub fn stats(&self) -> Stats {
        let np = self.pos.len() as f64;
        let mut mean = [0.0; 3];
        for d in &self.disp {
            for c in 0..3 {
                mean[c] += d[c];
            }
        }
        for m in &mut mean {
            *m /= np;
        }
        let (mut msd, mut var) = (0.0, 0.0);
        for d in &self.disp {
            let r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
            msd += r2;
            var += r2
                - 2.0 * (d[0] * mean[0] + d[1] * mean[1] + d[2] * mean[2])
                + mean[0] * mean[0]
                + mean[1] * mean[1]
                + mean[2] * mean[2];
        }
        Stats { mean, msd: msd / np, msd_about_mean: var / np }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 120 nm virion in 20 °C air: C_c ≈ 2.54, D ≈ 5.0e-10 m²/s,
    /// v_s ≈ 1.2 µm/s. Values cross-checked against the aerosol literature.
    #[test]
    fn stokes_einstein_air() {
        let m = Medium::air();
        let v = Virion::influenza();
        assert!((v.cunningham(&m) - 2.545).abs() < 0.01);
        let d = v.diffusivity(&m);
        assert!((d - 5.02e-10).abs() / 5.02e-10 < 0.01, "D = {d:e}");
        let vs = v.settling_speed(&m);
        assert!((vs - 1.21e-6).abs() / 1.21e-6 < 0.01, "v_s = {vs:e}");
        assert!(v.relaxation_time(&m) < 2e-7, "overdamped: τ ≪ fluid timescales");
    }

    /// In water: no slip correction, D ≈ 3.6e-12 m²/s, settling negligible.
    #[test]
    fn stokes_einstein_water() {
        let m = Medium::water();
        let v = Virion::influenza();
        assert_eq!(v.cunningham(&m), 1.0);
        let d = v.diffusivity(&m);
        assert!((d - 3.57e-12).abs() / 3.57e-12 < 0.01, "D = {d:e}");
        assert!(v.settling_speed(&m) < 1e-9);
    }

    /// Trilinear sampling is exact on a linear field (u = x · x̂, sampled
    /// away from the wrap seam).
    #[test]
    fn trilinear_exact_on_linear() {
        let n = 16;
        let dx = TWO_PI / n as f64;
        let mut u = [vec![0.0; n * n * n], vec![0.0; n * n * n], vec![0.0; n * n * n]];
        for i in 0..n {
            for j in 0..n {
                for l in 0..n {
                    u[0][(i * n + j) * n + l] = i as f64 * dx;
                }
            }
        }
        let p = [1.3 * dx + 0.31 * dx, 5.2 * dx, 9.7 * dx];
        let v = trilinear(&u, n, p);
        assert!((v[0] - p[0]).abs() < 1e-15);
        assert!(v[1] == 0.0 && v[2] == 0.0);
    }

    /// Einstein relation: quiescent fluid, ⟨|Δx|²⟩ = 6 D t.
    /// Statistical error σ/MSD = √(2/(3 np)) ≈ 0.4%; tolerance is 5σ.
    #[test]
    fn einstein_msd_quiescent() {
        let (m, v) = (Medium::air(), Virion::influenza());
        let d_si = v.diffusivity(&m);
        // Identity scaling: 1 solver length = 1 m, 1 solver time = 1 s.
        let sc = Scaling { box_len: TWO_PI, velocity: 1.0 };
        let d = sc.diffusivity_solver(d_si);
        let np = 50_000;
        let mut c = Cloud::new(np, d, 0.0, 0xC0FFEE);
        let (t_end, steps) = (1.0, 20);
        for _ in 0..steps {
            c.step(None, 0, t_end / steps as f64);
        }
        let st = c.stats();
        let exact = 6.0 * d * t_end;
        let rel = (st.msd_about_mean - exact).abs() / exact;
        assert!(rel < 0.02, "MSD {:.3}% off 6Dt", rel * 100.0);
        // No drift: |mean| ≤ 5σ, σ = √(2Dt/np).
        let sigma = (2.0 * d * t_end / np as f64).sqrt();
        for k in 0..3 {
            assert!(st.mean[k].abs() < 5.0 * sigma, "spurious drift {}", st.mean[k]);
        }
    }

    /// Uniform flow: drift = U t exactly (D = 0 isolates the advection).
    #[test]
    fn uniform_advection_exact() {
        let n = 8;
        let u_inf = [0.7, -0.2, 0.1];
        let u: [Vec<f64>; 3] =
            std::array::from_fn(|c| vec![u_inf[c]; n * n * n]);
        let mut c = Cloud::new(1_000, 0.0, 0.0, 7);
        let (t_end, steps) = (3.0, 30);
        for _ in 0..steps {
            c.step(Some(&u), n, t_end / steps as f64);
        }
        let st = c.stats();
        for k in 0..3 {
            assert!((st.mean[k] - u_inf[k] * t_end).abs() < 1e-12);
        }
        assert!(st.msd_about_mean < 1e-20, "deterministic: no spread");
    }

    /// Still air: mean descent = v_s t exactly (D = 0 isolates settling).
    #[test]
    fn settling_exact() {
        let (m, v) = (Medium::air(), Virion::influenza());
        let sc = Scaling { box_len: TWO_PI, velocity: 1.0 };
        let vs = sc.settling_solver(v.settling_speed(&m));
        let mut c = Cloud::new(1_000, 0.0, vs, 11);
        let (t_end, steps) = (2.0, 40);
        for _ in 0..steps {
            c.step(None, 0, t_end / steps as f64);
        }
        let st = c.stats();
        assert!((st.mean[2] + vs * t_end).abs() / (vs * t_end) < 1e-12);
    }

    /// All three at once in solver units with a nontrivial scaling: the
    /// Einstein term is independent of the drift terms.
    #[test]
    fn advection_settling_diffusion_together() {
        let (m, v) = (Medium::air(), Virion::influenza());
        let sc = Scaling { box_len: 2.0e-3 * TWO_PI, velocity: 0.05 };
        let d = sc.diffusivity_solver(v.diffusivity(&m));
        let vs = sc.settling_solver(v.settling_speed(&m));
        let n = 8;
        let u_inf = [1.0, 0.0, 0.0]; // 5 cm/s in SI
        let u: [Vec<f64>; 3] = std::array::from_fn(|c| vec![u_inf[c]; n * n * n]);
        let np = 100_000;
        let mut c = Cloud::new(np, d, vs, 99);
        let (t_end, steps) = (10.0, 100); // solver units: 10 · t0 = 0.4 s
        for _ in 0..steps {
            c.step(Some(&u), n, t_end / steps as f64);
        }
        let st = c.stats();
        let sigma = (2.0 * d * t_end / np as f64).sqrt();
        assert!((st.mean[0] - t_end).abs() < 5.0 * sigma + 1e-12, "drift x");
        assert!(st.mean[1].abs() < 5.0 * sigma, "drift y");
        assert!((st.mean[2] + vs * t_end).abs() < 5.0 * sigma, "settling");
        let exact = 6.0 * d * t_end;
        assert!((st.msd_about_mean - exact).abs() / exact < 0.015, "6Dt");
    }
}
