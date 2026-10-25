//! Massive and null geodesics in static isotropic metrics. Four-momenta are
//! local orthonormal-frame quantities; x,t are coordinates. Local decay laws
//! are applied at the proper-time event. No self-consistent Einstein evolution.
use crate::{decay, transport::Flight, Domain, FourMomentum, Particle, Rng, Simulation, Vec3, C};
use gravity::{
    metric::{Geometry, Metric},
    norm,
};
use rayon::prelude::*;
#[derive(Clone, Copy, Debug)]
pub struct Options {
    pub relative: f64,
    pub position_absolute: f64,
    pub max_step: f64,
    pub max_attempts: usize,
}
impl Default for Options {
    fn default() -> Self {
        Self {
            relative: 1e-10,
            position_absolute: 1e-3,
            max_step: 1e4,
            max_attempts: 100_000,
        }
    }
}
impl Options {
    fn validate(self) -> Result<(), &'static str> {
        if !self.relative.is_finite()
            || self.relative <= 0.0
            || self.relative > 1e-2
            || !self.position_absolute.is_finite()
            || self.position_absolute <= 0.0
            || !self.max_step.is_finite()
            || self.max_step <= 0.0
            || self.max_attempts == 0
        {
            Err("invalid geodesic tolerances")
        } else {
            Ok(())
        }
    }
}
#[derive(Clone, Copy)]
struct State {
    x: Vec3,
    p: Vec3,
    proper: f64,
}
impl State {
    fn add(self, b: Self, h: f64) -> Self {
        Self {
            x: self.x + b.x * h,
            p: self.p + b.p * h,
            proper: self.proper + b.proper * h,
        }
    }
}
fn derivative(y: State, m: f64, g: Geometry) -> Result<State, &'static str> {
    let p = y.p / g.scale;
    let e = m.hypot(norm(p));
    if !e.is_finite() || e == 0.0 {
        return Err("invalid geodesic energy");
    }
    let velocity = p * (C * g.lapse / (g.scale * e));
    let force =
        (g.grad_lapse * (-e) + g.grad_scale * (g.lapse * (norm(p) / e) * norm(p) / g.scale)) * C;
    Ok(State {
        x: velocity,
        p: force,
        proper: g.lapse * (m / e),
    })
}
fn rk4<M: Metric>(field: &M, y: State, m: f64, h: f64) -> Result<State, &'static str> {
    let a = derivative(y, m, field.geometry(y.x)?.validate()?)?;
    let z = y.add(a, h * 0.5);
    let b = derivative(z, m, field.geometry(z.x)?.validate()?)?;
    let z = y.add(b, h * 0.5);
    let c = derivative(z, m, field.geometry(z.x)?.validate()?)?;
    let z = y.add(c, h);
    let d = derivative(z, m, field.geometry(z.x)?.validate()?)?;
    Ok(y.add(a, h / 6.0)
        .add(b, h / 3.0)
        .add(c, h / 3.0)
        .add(d, h / 6.0))
}
#[derive(Clone, Copy, Debug)]
pub struct Orbit {
    pub position: Vec3,
    pub four: FourMomentum,
    pub coordinate_time: f64,
    pub proper_time: f64,
    pub decayed: bool,
    pub attempts: usize,
}
/// Adaptive RK4 with step doubling for the nonseparable geodesic Hamiltonian
/// H=alpha sqrt(m²+P²/beta²). This is an accuracy-controlled test-particle
/// integrator, not a long-term symplectic N-body evolution method.
#[allow(clippy::too_many_arguments)] // physical initial data, both clocks, and tolerances
pub fn propagate<M: Metric>(
    metric: &M,
    position: Vec3,
    four: FourMomentum,
    mass: f64,
    dt: f64,
    proper_limit: f64,
    opts: Options,
) -> Result<Orbit, &'static str> {
    opts.validate()?;
    four.validate_on_shell(mass)?;
    if !dt.is_finite() || dt < 0.0 || proper_limit.is_nan() || proper_limit < 0.0 {
        return Err("invalid gravity flight interval");
    }
    let initial = metric.geometry(position)?.validate()?;
    let mut y = State {
        x: position,
        p: four.momentum * initial.scale,
        proper: 0.0,
    };
    let mut time = 0.0;
    let mut h = dt.min(opts.max_step);
    let mut attempts = 0;
    let mut event = proper_limit == 0.0;
    while time < dt && !event {
        attempts += 1;
        if attempts > opts.max_attempts {
            return Err("gravity flight exceeded integration budget");
        }
        let g = metric.geometry(y.x)?.validate()?;
        let rate = derivative(y, mass, g)?.proper;
        h = h.min(dt - time).min(opts.max_step);
        if proper_limit.is_finite() && rate > 0.0 {
            h = h.min(1.02 * (proper_limit - y.proper) / rate);
        }
        if h <= 0.0 || time + h == time {
            return Err("unrepresentable geodesic step");
        }
        let full = rk4(metric, y, mass, h)?;
        let half = rk4(metric, y, mass, h * 0.5)?;
        let refined = rk4(metric, half, mass, h * 0.5)?;
        let increment = refined.proper - y.proper;
        if refined.proper > proper_limit * (1.0 + 2e-13) {
            h *= ((proper_limit - y.proper) / increment).clamp(0.01, 0.999999);
            continue;
        }
        let ex = norm(refined.x - full.x)
            / (15.0 * (opts.position_absolute + opts.relative * norm(refined.x - position)));
        let ep = norm(refined.p - full.p) / (15.0 * opts.relative * (mass + norm(refined.p)));
        let et = if mass > 0.0 {
            (refined.proper - full.proper).abs()
                / (15.0 * opts.relative * refined.proper.max(h * rate))
        } else {
            0.0
        };
        let error = ex.max(ep).max(et);
        if !error.is_finite() {
            return Err("unrepresentable geodesic integration error");
        }
        if error > 1.0 {
            h *= (0.9 * error.powf(-0.2)).clamp(0.1, 0.8);
            continue;
        }
        y = refined;
        time += h;
        if proper_limit.is_finite() && proper_limit - y.proper <= 4e-12 * proper_limit {
            event = true;
            y.proper = proper_limit;
        }
        h *= if error > 1e-20 {
            (0.9 * error.powf(-0.2)).clamp(0.5, 4.0)
        } else {
            4.0
        };
    }
    let g = metric.geometry(y.x)?.validate()?;
    let out = FourMomentum::on_shell(mass, y.p / g.scale)?;
    if !time.is_finite() || !y.proper.is_finite() {
        return Err("unrepresentable gravity flight");
    }
    Ok(Orbit {
        position: y.x,
        four: out,
        coordinate_time: time,
        proper_time: y.proper,
        decayed: event,
        attempts,
    })
}
pub fn fly<M: Metric>(
    particle: Particle,
    dt: f64,
    domain: Domain,
    metric: &M,
    opts: Options,
) -> Result<Flight, &'static str> {
    let mut queue = vec![(particle, dt)];
    let mut flight = Flight::default();
    while let Some((mut p, remaining)) = queue.pop() {
        if !p.species.allowed(domain) {
            return Err("gravity particle/domain mismatch");
        }
        let orbit = propagate(
            metric,
            p.position,
            p.four,
            p.species.mass(),
            remaining,
            p.proper_left,
            opts,
        )?;
        flight.field_exchange = flight.field_exchange + (orbit.four - p.four);
        p.position = orbit.position;
        p.four = orbit.four;
        p.time += orbit.coordinate_time;
        p.proper_elapsed += orbit.proper_time;
        if !p.time.is_finite() {
            return Err("unrepresentable gravity coordinate time");
        }
        if p.proper_left.is_finite() {
            p.proper_left = (p.proper_left - orbit.proper_time).max(0.0);
        }
        if !orbit.decayed {
            flight.live.push(p);
            continue;
        }
        let mut rng = Rng(p.rng_state);
        let products = decay::sample(p.species, p.four, domain, &mut rng)?;
        flight.decays += 1;
        for (origin, reason, four, charges) in products.unresolved {
            flight.unresolved.push(decay::Unresolved {
                origin,
                reason,
                four,
                charges,
                position: p.position,
                time: p.time,
            });
        }
        for (species, four) in products.particles {
            let child = Particle::from_four(species, p.position, four, p.time, rng.next_u64())?;
            queue.push((child, (remaining - orbit.coordinate_time).max(0.0)));
        }
    }
    Ok(flight)
}
impl Simulation {
    /// Gravity interacts with massive particles AND photons/massless neutrinos.
    /// Source geometry is frozen during this call. Returns changes in summed
    /// LOCAL tetrad four-momenta for bookkeeping, not a conserved global GR
    /// vector. Use killing_energy() for the stationary energy invariant.
    pub fn advance_gravity<M: Metric>(
        &mut self,
        dt: f64,
        metric: &M,
        opts: Options,
    ) -> Result<FourMomentum, &'static str> {
        opts.validate()?;
        if !dt.is_finite() || dt < 0.0 || !(self.time + dt).is_finite() {
            return Err("invalid gravity simulation interval");
        }
        if dt == 0.0 {
            return Ok(FourMomentum::default());
        }
        let run = |p: &Particle| fly(p.clone(), dt, self.domain, metric, opts);
        let flights: Result<Vec<_>, _> = if self.particles.len() >= 1024 {
            self.particles.par_iter().map(run).collect()
        } else {
            self.particles.iter().map(run).collect()
        };
        let flights = flights?;
        let mut exchange = FourMomentum::default();
        self.particles.clear();
        for f in flights {
            exchange = exchange + f.field_exchange;
            self.particles.extend(f.live);
            self.unresolved.extend(f.unresolved);
            self.decays += f.decays;
        }
        self.time += dt;
        Ok(exchange)
    }
    /// Sum alpha(x) E_local [GeV], including retained unresolved records.
    /// Conserved for local interactions/decays in this stationary metric.
    pub fn killing_energy<M: Metric>(&self, metric: &M) -> Result<f64, &'static str> {
        let mut sum = 0.0;
        for p in &self.particles {
            sum += metric.geometry(p.position)?.validate()?.lapse * p.four.energy;
        }
        for p in &self.unresolved {
            sum += metric.geometry(p.position)?.validate()?.lapse * p.four.energy;
        }
        if !sum.is_finite() {
            return Err("unrepresentable Killing energy");
        }
        Ok(sum)
    }
}
