//! Photoacoustics: light absorbed in a short pulse launches a pressure wave.
//!
//! Under **stress confinement** the absorbed energy cannot relax during the
//! pulse, so the initial pressure is
//!
//! \[
//!   p_0(z) = \Gamma\, \mu_a F\, e^{-\mu_a z},
//!   \qquad \Gamma = \frac{\beta c^2}{c_p},
//! \]
//!
//! with `F` the fluence [J m^-2], `mu_a` the absorption coefficient [m^-1]
//! and `Gamma` the Grüneisen parameter. The 1D acoustic wave
//! `p_tt = c^2 p_zz` with `p(z,0)=p_0`, `p_t(z,0)=0` has the exact
//! d'Alembert solution
//!
//! \[
//!   p(z,t) = \tfrac12\bigl[P_0(z-ct) + P_0(z+ct)\bigr],
//! \]
//!
//! where `P_0` is the odd extension of `p_0` for a **free** (pressure-release)
//! surface and the even extension for a **rigid** surface. The free surface
//! therefore gives a bipolar (N-shaped) signal and the rigid surface a
//! monopolar one: the classic photoacoustic signature difference.
//!
//! For a layered stack, a time-domain finite-difference solver (pressure at
//! nodes, particle velocity at staggered faces, exact interface averaging)
//! propagates the initial condition. As an independent route, a pulse can be
//! launched in the frequency domain: FFT the incident pressure, multiply by
//! the layer stack's pressure transmission coefficient from
//! [`phononic::layered`], then inverse FFT. Both routes are validated against
//! the d'Alembert solution.
//!
//! Grüneisen parameter of water at 20 C: `beta = 2.07e-4 K^-1` (volumetric
//! thermal expansion, CRC Handbook), `c = 1482 m/s`, `c_p = 4182 J kg^-1 K^-1`,
//! giving `Gamma = beta c^2 / c_p = 0.1088`; the commonly quoted round value
//! is `0.11`.

use crate::fft::{fft, Z};
use phononic::layered;
use std::f64::consts::PI;

/// Grüneisen parameter of water at 20 C (see the module note).
pub const GRUNEISEN_WATER: f64 = 0.1088;

/// Grüneisen parameter from its definition `Gamma = beta c^2 / c_p`.
pub fn gruneisen(beta: f64, c: f64, cp: f64) -> f64 {
    beta * c * c / cp
}

/// Initial pressure [Pa] of an exponential absorber, `Gamma mu_a F e^{-mu_a z}`.
pub fn initial_pressure(mu_a: f64, fluence: f64, gamma: f64, z: f64) -> f64 {
    gamma * mu_a * fluence * (-mu_a * z).exp()
}

/// Peak initial pressure `Gamma mu_a F` at the surface.
pub fn initial_peak_pressure(mu_a: f64, fluence: f64, gamma: f64) -> f64 {
    gamma * mu_a * fluence
}

/// Surface condition for the semi-infinite d'Alembert solution.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Surface {
    /// Pressure-release (`p = 0`): odd image, bipolar signal.
    Free,
    /// Rigid (`v = 0`): even image, monopolar signal.
    Rigid,
}

/// Exact d'Alembert pressure on `z > 0` for an initial pressure profile `p0`
/// defined there, with the image convention of `surface`.
pub fn dalembert<F: Fn(f64) -> f64>(p0: F, z: f64, t: f64, c: f64, surface: Surface) -> f64 {
    let zm = z - c * t;
    let zp = z + c * t;
    let left = if zm >= 0.0 {
        p0(zm)
    } else {
        match surface {
            Surface::Free => -p0(-zm),
            Surface::Rigid => p0(-zm),
        }
    };
    0.5 * (left + p0(zp))
}

/// Pressure on a rigid semi-infinite surface for an exponential initial
/// profile: `p(0,t) = Gamma mu_a F exp(-mu_a c t)` (monopolar).
pub fn rigid_surface_pressure(mu_a: f64, fluence: f64, gamma: f64, c: f64, t: f64) -> f64 {
    initial_peak_pressure(mu_a, fluence, gamma) * (-mu_a * c * t).exp()
}

/// One lossless acoustic layer.
#[derive(Clone, Copy, Debug)]
pub struct AcousticLayer {
    pub rho: f64,
    pub c: f64,
    pub thickness: f64,
}

impl AcousticLayer {
    pub fn new(rho: f64, c: f64, thickness: f64) -> Self {
        assert!(rho > 0.0 && c > 0.0 && thickness > 0.0);
        Self { rho, c, thickness }
    }
    pub fn impedance(&self) -> f64 {
        self.rho * self.c
    }
    pub fn from_fluid(f: &phononic::materials::Fluid, thickness: f64) -> Self {
        Self::new(f.rho, f.c, thickness)
    }
}

/// Far-field boundary of the finite acoustic grid.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FarBoundary {
    Rigid,
    Free,
    /// First-order characteristic (impedance) outflow.
    Absorbing,
}

/// Staggered 1D acoustic grid: pressure at nodes, velocity at the faces
/// between them. Layer boundaries are snapped to the nearest node.
pub struct AcousticModel {
    pub n: usize,
    pub dz: f64,
    pub z: Vec<f64>,
    /// Density on each face (length `n-1`).
    pub rho_face: Vec<f64>,
    /// Effective bulk modulus at each node, `2 K_L K_R/(K_L+K_R)` at an
    /// interface so that `p` and `v` stay continuous.
    pub k_node: Vec<f64>,
    pub c_max: f64,
    pub surface: Surface,
    pub far: FarBoundary,
}

impl AcousticModel {
    /// Build a grid from a layer sequence. `dz` should divide every layer
    /// thickness; the boundary is snapped to the closest node.
    pub fn from_layers(
        layers: &[AcousticLayer],
        dz: f64,
        surface: Surface,
        far: FarBoundary,
    ) -> Self {
        assert!(!layers.is_empty() && dz > 0.0);
        let total: f64 = layers.iter().map(|l| l.thickness).sum();
        let n = (total / dz).round() as usize + 1;
        let dz = total / (n - 1) as f64;
        // Interval material for interval i = [z_i, z_{i+1}].
        let mut rho_i = vec![0.0; n - 1];
        let mut k_i = vec![0.0; n - 1];
        let mut edges = Vec::with_capacity(layers.len() + 1);
        let mut acc = 0.0;
        edges.push(0.0);
        for l in layers {
            acc += l.thickness;
            edges.push(acc);
        }
        for i in 0..n - 1 {
            let zm = dz * (i as f64 + 0.5);
            let mut idx = layers.len() - 1;
            for (j, _) in layers.iter().enumerate() {
                if zm <= edges[j + 1] + 1e-12 {
                    idx = j;
                    break;
                }
            }
            let c = layers[idx].c;
            let rho = layers[idx].rho;
            rho_i[i] = rho;
            k_i[i] = rho * c * c;
        }
        let mut k_node = vec![0.0; n];
        k_node[0] = k_i[0];
        k_node[n - 1] = k_i[n - 2];
        for i in 1..n - 1 {
            k_node[i] = 2.0 * k_i[i - 1] * k_i[i] / (k_i[i - 1] + k_i[i]);
        }
        let c_max = layers.iter().map(|l| l.c).fold(0.0f64, f64::max);
        let z = (0..n).map(|i| i as f64 * dz).collect();
        Self {
            n,
            dz,
            z,
            rho_face: rho_i,
            k_node,
            c_max,
            surface,
            far,
        }
    }

    /// Stability-limited step, `cfl * dz / c_max`.
    pub fn max_dt(&self, cfl: f64) -> f64 {
        cfl * self.dz / self.c_max
    }

    /// Advance one step. `p` has length `n`, `v` length `n-1`.
    pub fn step(&self, p: &mut [f64], v: &mut [f64], dt: f64) {
        let n = self.n;
        // Momentum: rho dv/dt = -dp/dz at the faces.
        for i in 0..n - 1 {
            v[i] -= dt / (self.rho_face[i] * self.dz) * (p[i + 1] - p[i]);
        }
        // Continuity: (1/K) dp/dt = -dv/dz at interior nodes.
        for i in 1..n - 1 {
            p[i] -= dt * self.k_node[i] / self.dz * (v[i] - v[i - 1]);
        }
        match self.surface {
            Surface::Rigid => {
                // Ghost face velocity zero at the wall.
                p[0] -= dt * self.k_node[0] / self.dz * (v[0] - 0.0);
            }
            Surface::Free => {
                p[0] = 0.0;
            }
        }
        match self.far {
            FarBoundary::Rigid => {
                p[n - 1] -= dt * self.k_node[n - 1] / self.dz * (0.0 - v[n - 2]);
            }
            FarBoundary::Free => {
                p[n - 1] = 0.0;
            }
            FarBoundary::Absorbing => {
                let zl = self.rho_face[n - 2] * (self.k_node[n - 1] / self.rho_face[n - 2]).sqrt();
                let v_ghost = p[n - 1] / zl;
                p[n - 1] -= dt * self.k_node[n - 1] / self.dz * (v_ghost - v[n - 2]);
            }
        }
    }

    /// Initialise a pressure profile `p0(z)` with the particle velocity at the
    /// half step `-dt/2`, so the first [`step`](Self::step) lands at `+dt/2`
    /// and the scheme is second-order in time from the start. For a free
    /// surface the surface node is set to zero, consistent with `p(0,t)=0`.
    pub fn initial_pressure_profile<F: Fn(f64) -> f64>(
        &self,
        p0: F,
        dt: f64,
    ) -> (Vec<f64>, Vec<f64>) {
        let mut p: Vec<f64> = self.z.iter().map(|&z| p0(z)).collect();
        if self.surface == Surface::Free {
            p[0] = 0.0;
        }
        // v^{-1/2} = (dt/2rho) dp/dz so that v^{+1/2} = -(dt/2rho) dp/dz,
        // the exact half-step velocity of the d'Alembert split.
        let v: Vec<f64> = (0..self.n - 1)
            .map(|i| 0.5 * dt / (self.rho_face[i] * self.dz) * (p[i + 1] - p[i]))
            .collect();
        (p, v)
    }

    /// Run `steps` steps, returning the final `(p, v)`.
    pub fn run(
        &self,
        mut p: Vec<f64>,
        mut v: Vec<f64>,
        dt: f64,
        steps: usize,
    ) -> (Vec<f64>, Vec<f64>) {
        for _ in 0..steps {
            self.step(&mut p, &mut v, dt);
        }
        (p, v)
    }
}

/// Frequency-domain transmission of a sampled incident pressure pulse through
/// a [`phononic::layered::Stack`].
///
/// `p_inc` is the incident pressure at the front face of the stack, sampled
/// uniformly at `dt`; the returned vector is the transmitted pressure at the
/// back face over the same time base. The transfer function is the pressure
/// transmission coefficient `t(omega)` from the phononic transfer matrix, with
/// Hermitian symmetry `H(-omega) = conj H(omega)`. Use `None` for a
/// homogeneous medium (`H = 1`). The length of `p_inc` must be a power of two.
pub fn spectral_transmit(p_inc: &[f64], dt: f64, stack: Option<&layered::Stack>) -> Vec<f64> {
    let n = p_inc.len();
    assert!(n.is_power_of_two());
    let mut spec = crate::fft::fft_real(p_inc);
    for (k, s) in spec.iter_mut().enumerate() {
        let kk = if k <= n / 2 {
            k as isize
        } else {
            k as isize - n as isize
        };
        // Hermitian transfer function: evaluate H at |omega| and conjugate it
        // for negative frequencies, so the real input maps to a real output.
        // The zero-frequency (DC) limit of a finite stack is the static
        // pressure transmission `2 Z_exit / (Z_in + Z_exit)`, exact because
        // every layer matrix tends to the identity as omega -> 0.
        let h = if kk == 0 {
            match stack {
                Some(st) => Z::real(2.0 * st.exit_z.re / (st.incident_z.re + st.exit_z.re)),
                None => Z::real(1.0),
            }
        } else {
            let omega = 2.0 * PI * kk.unsigned_abs() as f64 / (n as f64 * dt);
            let h_pos = match stack {
                Some(st) => {
                    let tc = st.solve(omega).t;
                    Z::new(tc.re, tc.im)
                }
                None => Z::real(1.0),
            };
            if kk < 0 {
                h_pos.conj()
            } else {
                h_pos
            }
        };
        *s = s.mul(h);
    }
    // The transfer function is Hermitian, so the output is real.
    let _ = fft(&mut spec, true);
    spec.iter().map(|z| z.re).collect()
}

/// A Gaussian-modulated tone burst with centre frequency `f0` and fractional
/// bandwidth `bw`, sampled at `dt` for `n` points. Useful as a photoacoustic
/// pulse for the spectral launch. Peak amplitude is 1.
pub fn tone_burst(n: usize, dt: f64, f0: f64, bw: f64) -> Vec<f64> {
    assert!(n.is_power_of_two());
    let t0 = 0.5 * n as f64 * dt;
    let sigma = bw / (2.0 * PI * f0); // envelope width from the fractional bandwidth
    (0..n)
        .map(|i| {
            let t = i as f64 * dt;
            let env = (-((t - t0) / sigma).powi(2) / 2.0).exp();
            env * (2.0 * PI * f0 * (t - t0)).cos()
        })
        .collect()
}
