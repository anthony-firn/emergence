//! Surface radial magnetic flux transport on equal-area mu=sin(latitude) cells.
//! Prescribed flows/effective diffusion; source bipoles must be supplied.
use crate::R_SUN;
use num_complex::Complex64;
use rustfft::{Fft, FftPlanner};
use std::{
    f64::consts::{PI, TAU},
    sync::Arc,
};
#[derive(Clone, Copy, Debug)]
pub struct Parameters {
    pub radius: f64,
    /// Poleward v_lat=v_peak sin(2 latitude), m/s.
    pub meridional_peak: f64,
    /// Effective horizontal diffusivity, m²/s; not microscopic resistivity.
    pub diffusivity: f64,
    /// Omega(mu)=a+b mu²+c mu⁴, rad/s in caller's reference frame.
    pub rotation: [f64; 3],
    pub max_step: f64,
}
impl Default for Parameters {
    fn default() -> Self {
        let c = PI / (180.0 * 86400.0);
        Self {
            radius: R_SUN,
            meridional_peak: 15.0,
            diffusivity: 250e6,
            rotation: [(13.38 - 13.2) * c, -2.30 * c, -1.62 * c],
            max_step: 86400.0,
        }
    }
}
impl Parameters {
    fn validate(self) -> Result<(), &'static str> {
        if !self.radius.is_finite()
            || self.radius <= 0.0
            || !self.meridional_peak.is_finite()
            || !self.diffusivity.is_finite()
            || self.diffusivity < 0.0
            || self.rotation.iter().any(|v| !v.is_finite())
            || !self.max_step.is_finite()
            || self.max_step <= 0.0
        {
            return Err("invalid solar transport coefficients");
        }
        Ok(())
    }
}
const DIRK_GAMMA: f64 = 1.0 - std::f64::consts::FRAC_1_SQRT_2;
struct Tridiagonal {
    elimination: Vec<f64>,
    inverse_diagonal: Vec<f64>,
    upper: Vec<f64>,
}
impl Tridiagonal {
    fn solve(&self, rhs: &mut [Complex64]) {
        for i in 1..rhs.len() {
            let previous = rhs[i - 1];
            rhs[i] -= previous * self.elimination[i];
        }
        let last = rhs.len() - 1;
        rhs[last] *= self.inverse_diagonal[last];
        for i in (0..last).rev() {
            let next = rhs[i + 1];
            rhs[i] = (rhs[i] - next * self.upper[i]) * self.inverse_diagonal[i];
        }
    }
}
struct DiffusionPlan {
    factors: Vec<Tridiagonal>,
}
impl DiffusionPlan {
    fn new(s: &Surface, dt: f64) -> Result<Option<Self>, &'static str> {
        if s.parameters.diffusivity == 0.0 {
            return Ok(None);
        }
        let n = s.n_mu;
        let h = 2.0 / n as f64;
        let rate = s.parameters.diffusivity / s.parameters.radius.powi(2);
        let faces: Vec<_> = (0..=n)
            .map(|i| {
                let mu = -1.0 + i as f64 * h;
                rate * (1.0 - mu * mu).max(0.0) / (h * h)
            })
            .collect();
        let mut factors = Vec::new();
        // Positive and negative longitude modes share the same real matrix.
        for mode in 0..=s.n_phi / 2 {
            let mut diagonal = vec![0.0; n];
            let mut upper = vec![0.0; n];
            let mut elimination = vec![0.0; n];
            for i in 0..n {
                let a = DIRK_GAMMA * dt * faces[i];
                let b = DIRK_GAMMA * dt * faces[i + 1];
                let lambda =
                    DIRK_GAMMA * dt * rate * (mode as f64).powi(2) / (1.0 - s.mu(i).powi(2));
                diagonal[i] = 1.0 + a + b + lambda;
                upper[i] = -b;
                if i > 0 {
                    elimination[i] = -a / diagonal[i - 1];
                    diagonal[i] -= elimination[i] * upper[i - 1];
                }
                if !diagonal[i].is_finite() || diagonal[i] <= 0.0 {
                    return Err("unrepresentable solar diffusion matrix");
                }
            }
            factors.push(Tridiagonal {
                elimination,
                upper,
                inverse_diagonal: diagonal.iter().map(|v| 1.0 / v).collect(),
            });
        }
        Ok(Some(Self { factors }))
    }
}
pub struct Surface {
    n_mu: usize,
    n_phi: usize,
    parameters: Parameters,
    field: Vec<f64>,
    time: f64,
    forward: Arc<dyn Fft<f64>>,
    inverse: Arc<dyn Fft<f64>>,
}
impl Surface {
    /// Equal-area cell averages in tesla, south to north and increasing longitude.
    pub fn new(
        n_mu: usize,
        n_phi: usize,
        parameters: Parameters,
        field: Vec<f64>,
    ) -> Result<Self, &'static str> {
        parameters.validate()?;
        if n_mu < 4
            || n_phi == 0
            || field.len() != n_mu.checked_mul(n_phi).ok_or("solar grid overflow")?
            || field.iter().any(|v| !v.is_finite())
        {
            return Err("invalid solar magnetogram grid");
        }
        let mut planner = FftPlanner::new();
        let forward = planner.plan_fft_forward(n_phi);
        let inverse = planner.plan_fft_inverse(n_phi);
        let mut field = field;
        // Remove the even-grid Nyquist mode whose rotating sine is unresolvable.
        if n_phi > 1 && n_phi.is_multiple_of(2) {
            for row in field.chunks_mut(n_phi) {
                let a = row
                    .iter()
                    .enumerate()
                    .map(|(j, b)| if j % 2 == 0 { *b } else { -*b })
                    .sum::<f64>()
                    / n_phi as f64;
                for (j, b) in row.iter_mut().enumerate() {
                    *b -= if j % 2 == 0 { a } else { -a };
                }
            }
        }
        Ok(Self {
            n_mu,
            n_phi,
            parameters,
            field,
            time: 0.0,
            forward,
            inverse,
        })
    }
    pub fn dimensions(&self) -> (usize, usize) {
        (self.n_mu, self.n_phi)
    }
    pub fn field(&self) -> &[f64] {
        &self.field
    }
    pub fn time(&self) -> f64 {
        self.time
    }
    pub fn parameters(&self) -> Parameters {
        self.parameters
    }
    pub fn mu(&self, i: usize) -> f64 {
        -1.0 + (i as f64 + 0.5) * 2.0 / self.n_mu as f64
    }
    pub fn longitude(&self, j: usize) -> f64 {
        (j as f64 + 0.5) * TAU / self.n_phi as f64
    }
    pub fn mean_field(&self) -> Vec<f64> {
        self.field
            .chunks(self.n_phi)
            .map(|r| r.iter().sum::<f64>() / self.n_phi as f64)
            .collect()
    }
    /// Exact discrete fold for longitude-independent linear transport; retains
    /// the future mean/dipole, discarding local magnetic-region information.
    pub fn axisymmetric(&self) -> Result<Self, &'static str> {
        let mut s = Self::new(self.n_mu, 1, self.parameters, self.mean_field())?;
        s.time = self.time;
        Ok(s)
    }
    pub fn axial_dipole(&self) -> f64 {
        1.5 * 2.0 / self.n_mu as f64
            * self
                .mean_field()
                .iter()
                .enumerate()
                .map(|(i, b)| b * self.mu(i))
                .sum::<f64>()
    }
    fn area(&self) -> f64 {
        4.0 * PI * self.parameters.radius.powi(2) / (self.n_mu * self.n_phi) as f64
    }
    pub fn net_flux(&self) -> f64 {
        self.field.iter().sum::<f64>() * self.area()
    }
    pub fn unsigned_flux(&self) -> f64 {
        self.field.iter().map(|b| b.abs()).sum::<f64>() * self.area()
    }
    /// Explicit observation preprocessing; returns uniform B removed in tesla.
    pub fn balance_flux(&mut self) -> f64 {
        let mean = self.field.iter().sum::<f64>() / self.field.len() as f64;
        for b in &mut self.field {
            *b -= mean;
        }
        mean
    }
    fn spectrum(&self, field: &[f64]) -> Vec<Complex64> {
        let mut s: Vec<_> = field.iter().map(|v| Complex64::new(*v, 0.0)).collect();
        let mut scratch = vec![Complex64::default(); self.forward.get_inplace_scratch_len()];
        for row in s.chunks_mut(self.n_phi) {
            self.forward.process_with_scratch(row, &mut scratch);
            for v in row {
                *v /= self.n_phi as f64;
            }
        }
        s
    }
    fn physical(&self, s: &mut [Complex64]) -> Vec<f64> {
        let mut scratch = vec![Complex64::default(); self.inverse.get_inplace_scratch_len()];
        for row in s.chunks_mut(self.n_phi) {
            self.inverse.process_with_scratch(row, &mut scratch);
        }
        s.iter().map(|v| v.re).collect()
    }
    fn rotate(&self, s: &mut [Complex64], dt: f64) {
        for i in 0..self.n_mu {
            let mu = self.mu(i);
            let r = self.parameters.rotation;
            let omega = r[0] + r[1] * mu * mu + r[2] * mu.powi(4);
            for j in 0..self.n_phi {
                let m = if j <= self.n_phi / 2 {
                    j as f64
                } else {
                    j as f64 - self.n_phi as f64
                };
                s[i * self.n_phi + j] *= Complex64::from_polar(1.0, -m * omega * dt);
            }
        }
    }
    fn diffuse(&self, spectrum: &mut [Complex64], plan: &Option<DiffusionPlan>) {
        let Some(plan) = plan else {
            return;
        };
        let mut original = vec![Complex64::default(); self.n_mu];
        let mut rhs = original.clone();
        let weight = (1.0 - DIRK_GAMMA) / DIRK_GAMMA;
        for mode in 0..self.n_phi {
            let factor = &plan.factors[mode.min(self.n_phi - mode)];
            for i in 0..self.n_mu {
                original[i] = spectrum[i * self.n_phi + mode];
                rhs[i] = original[i];
            }
            factor.solve(&mut rhs);
            // Stiffly accurate second stage: (I-gamma*h*L)Y2 =
            // Y0 + (1-gamma)/gamma*(Y1-Y0). Reuse the same factorization.
            for (value, initial) in rhs.iter_mut().zip(&original) {
                *value = *initial + (*value - *initial) * weight;
            }
            factor.solve(&mut rhs);
            for i in 0..self.n_mu {
                spectrum[i * self.n_phi + mode] = rhs[i];
            }
        }
    }
    fn advection_rhs(&self, b: &[f64]) -> Vec<f64> {
        let h = 2.0 / self.n_mu as f64;
        let mut out = vec![0.0; b.len()];
        for face in 1..self.n_mu {
            let mu = -1.0 + face as f64 * h;
            let speed = 2.0 * self.parameters.meridional_peak / self.parameters.radius
                * mu
                * (1.0 - mu * mu);
            let donor = if speed >= 0.0 { face - 1 } else { face };
            let direction = if speed >= 0.0 { 1.0 } else { -1.0 };
            for j in 0..self.n_phi {
                let q = donor * self.n_phi + j;
                // Unlimited centered MUSCL reconstruction keeps the operator
                // linear, permitting the exact longitude-average fold.
                let slope = if donor > 0 && donor + 1 < self.n_mu {
                    0.25 * (b[q + self.n_phi] - b[q - self.n_phi])
                } else {
                    0.0
                };
                let flux = speed * (b[q] + direction * slope) / h;
                out[(face - 1) * self.n_phi + j] -= flux;
                out[face * self.n_phi + j] += flux;
            }
        }
        out
    }
    fn timestep(&self) -> f64 {
        let h = 2.0 / self.n_mu as f64;
        let rate = (1..self.n_mu)
            .map(|i| {
                let mu = -1.0 + i as f64 * h;
                (2.0 * self.parameters.meridional_peak / self.parameters.radius
                    * mu
                    * (1.0 - mu * mu))
                    .abs()
                    / h
            })
            .fold(0.0, f64::max);
        self.parameters.max_step.min(if rate > 0.0 {
            0.4 / rate
        } else {
            f64::INFINITY
        })
    }
    fn step(&self, field: &[f64], dt: f64, diffusion: &Option<DiffusionPlan>) -> Vec<f64> {
        let mut s = self.spectrum(field);
        self.rotate(&mut s, dt * 0.5);
        self.diffuse(&mut s, diffusion);
        let b = self.physical(&mut s);
        let a = self.advection_rhs(&b);
        let intermediate: Vec<_> = b.iter().zip(&a).map(|(b, a)| b + dt * a).collect();
        let c = self.advection_rhs(&intermediate);
        let advected: Vec<_> = b
            .iter()
            .zip(a.iter().zip(&c))
            .map(|(b, (a, c))| b + dt * 0.5 * (a + c))
            .collect();
        let mut s = self.spectrum(&advected);
        self.diffuse(&mut s, diffusion);
        self.rotate(&mut s, dt * 0.5);
        self.physical(&mut s)
    }
    /// Second-order split evolution, exact Fourier rotation and implicit
    /// L-stable diffusion. Signed fields can ring on underresolved gradients; there is
    /// no nonlinear slope limiter. dt accuracy must be checked by refinement.
    pub fn advance(&mut self, dt: f64) -> Result<usize, &'static str> {
        if !dt.is_finite() || dt < 0.0 || !(self.time + dt).is_finite() {
            return Err("invalid surface forecast interval");
        }
        if dt == 0.0 {
            return Ok(0);
        }
        let count = (dt / self.timestep()).ceil().max(1.0);
        if !count.is_finite() || count > 1_000_000.0 {
            return Err("too many solar transport steps");
        }
        let steps = count as usize;
        let h = dt / steps as f64;
        let diffusion = DiffusionPlan::new(self, 0.5 * h)?;
        let mut field = self.field.clone();
        for _ in 0..steps {
            field = self.step(&field, h, &diffusion);
            if field.iter().any(|b| !b.is_finite()) {
                return Err("unrepresentable solar magnetic field");
            }
        }
        self.field = field;
        self.time += dt;
        Ok(steps)
    }
    /// Supplied source event. flux is per polarity [Wb], angles are radians;
    /// positive tilt rotates the eastward bipole axis northward.
    pub fn add_bipole(&mut self, e: Bipole) -> Result<(), &'static str> {
        e.validate()?;
        if self.n_phi == 1 {
            return Err("resolve the source on a longitude grid before folding it");
        }
        let (sin, cos) = e.latitude.sin_cos();
        let (ps, pc) = e.longitude.sin_cos();
        let center = [cos * pc, cos * ps, sin];
        let north = [-sin * pc, -sin * ps, cos];
        let east = [-ps, pc, 0.0];
        let axis: [f64; 3] =
            std::array::from_fn(|d| east[d] * e.tilt.cos() + north[d] * e.tilt.sin());
        let positive: [f64; 3] = std::array::from_fn(|d| {
            center[d] * (e.separation * 0.5).cos() + axis[d] * (e.separation * 0.5).sin()
        });
        let negative: [f64; 3] = std::array::from_fn(|d| {
            center[d] * (e.separation * 0.5).cos() - axis[d] * (e.separation * 0.5).sin()
        });
        let mut plus = vec![0.0; self.field.len()];
        let mut minus = plus.clone();
        for i in 0..self.n_mu {
            let mu = self.mu(i);
            let cosine = (1.0 - mu * mu).sqrt();
            for j in 0..self.n_phi {
                let phi = self.longitude(j);
                let v = [cosine * phi.cos(), cosine * phi.sin(), mu];
                let q = i * self.n_phi + j;
                let dot = |a: [f64; 3]| (0..3).map(|d| v[d] * a[d]).sum::<f64>();
                plus[q] = ((dot(positive) - 1.0) / e.width.powi(2)).exp();
                minus[q] = ((dot(negative) - 1.0) / e.width.powi(2)).exp();
            }
        }
        let sp = plus.iter().sum::<f64>();
        let sm = minus.iter().sum::<f64>();
        if sp == 0.0 || sm == 0.0 {
            return Err("bipole below magnetogram resolution");
        }
        let mut next = self.field.clone();
        let factor = e.flux / self.area();
        for q in 0..next.len() {
            next[q] += factor * (plus[q] / sp - minus[q] / sm);
        }
        let normalized = Self::new(self.n_mu, self.n_phi, self.parameters, next)?;
        let discarded = normalized
            .field
            .iter()
            .zip(&self.field)
            .zip(plus.iter().zip(&minus))
            .map(|((new, old), (p, m))| (new - old - factor * (p / sp - m / sm)).abs())
            .sum::<f64>()
            * self.area();
        if discarded > 0.001 * 2.0 * e.flux {
            return Err(
                "bipole contains unresolved longitudinal power; refine grid or widen source",
            );
        }
        self.field = normalized.field;
        Ok(())
    }
}
#[derive(Clone, Copy, Debug)]
pub struct Bipole {
    pub latitude: f64,
    pub longitude: f64,
    pub separation: f64,
    pub tilt: f64,
    pub width: f64,
    pub flux: f64,
}
impl Bipole {
    fn validate(self) -> Result<(), &'static str> {
        let v = [
            self.latitude,
            self.longitude,
            self.separation,
            self.tilt,
            self.width,
            self.flux,
        ];
        if v.iter().any(|x| !x.is_finite())
            || self.latitude.abs() > PI * 0.5
            || self.separation <= 0.0
            || self.separation >= PI
            || self.width <= 0.0
            || self.width >= PI
            || self.flux < 0.0
        {
            Err("invalid magnetic bipole")
        } else {
            Ok(())
        }
    }
}
#[derive(Clone, Copy, Debug)]
pub struct Skill {
    pub rmse: f64,
    pub persistence_rmse: f64,
    pub mse_skill: Option<f64>,
}
/// Area-weighted held-out skill. Never fit parameters or future source terms
/// using the target map that is being evaluated.
pub fn score(
    observed: &[f64],
    forecast: &[f64],
    persistence: &[f64],
) -> Result<Skill, &'static str> {
    if observed.is_empty()
        || forecast.len() != observed.len()
        || persistence.len() != observed.len()
        || observed
            .iter()
            .chain(forecast)
            .chain(persistence)
            .any(|v| !v.is_finite())
    {
        return Err("invalid held-out magnetogram comparison");
    }
    let mse = observed
        .iter()
        .zip(forecast)
        .map(|(a, b)| (a - b).powi(2))
        .sum::<f64>()
        / observed.len() as f64;
    let base = observed
        .iter()
        .zip(persistence)
        .map(|(a, b)| (a - b).powi(2))
        .sum::<f64>()
        / observed.len() as f64;
    if !mse.is_finite() || !base.is_finite() {
        return Err("unrepresentable forecast score");
    }
    Ok(Skill {
        rmse: mse.sqrt(),
        persistence_rmse: base.sqrt(),
        mse_skill: if base > 0.0 {
            Some(1.0 - mse / base)
        } else {
            None
        },
    })
}
