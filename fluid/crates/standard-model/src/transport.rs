//! Ballistic transport with event-time decay cascades (no c/dx step limit).
//! Independent per-particle random streams make Rayon scheduling reproducible.
use crate::{
    decay::{self, Unresolved},
    species::{Charges, Domain},
    FourMomentum, Rng, Species, UniformField, Vec3, C, ELEMENTARY_CHARGE, GEV_J,
};
use rayon::prelude::*;
#[derive(Clone, Debug)]
pub struct Particle {
    pub(crate) species: Species,
    pub(crate) four: FourMomentum,
    pub(crate) position: Vec3,
    pub(crate) time: f64,
    pub(crate) proper_left: f64,
    pub(crate) proper_elapsed: f64,
    pub(crate) rng_state: u64,
}
impl Particle {
    pub fn new(
        species: Species,
        position: Vec3,
        momentum: Vec3,
        seed: u64,
        domain: Domain,
    ) -> Result<Self, &'static str> {
        if !species.allowed(domain) {
            return Err("colored particle requires Partonic domain or a QCD closure");
        }
        if position.0.iter().any(|v| !v.is_finite()) {
            return Err("invalid particle position");
        }
        Self::from_four(
            species,
            position,
            FourMomentum::on_shell(species.mass(), momentum)?,
            0.0,
            seed,
        )
    }
    pub(crate) fn from_four(
        species: Species,
        position: Vec3,
        four: FourMomentum,
        time: f64,
        seed: u64,
    ) -> Result<Self, &'static str> {
        let mut rng = Rng(seed);
        let proper_left = species
            .lifetime()
            .map(|tau| -tau * rng.uniform().ln())
            .unwrap_or(f64::INFINITY);
        Ok(Self {
            species,
            four,
            position,
            time,
            proper_left,
            proper_elapsed: 0.0,
            rng_state: rng.0,
        })
    }
    pub fn species(&self) -> Species {
        self.species
    }
    pub fn four(&self) -> FourMomentum {
        self.four
    }
    pub fn position(&self) -> Vec3 {
        self.position
    }
    pub fn time(&self) -> f64 {
        self.time
    }
    /// Massive-particle proper age [s]; massless species have no proper clock.
    pub fn proper_age(&self) -> f64 {
        self.proper_elapsed
    }
    pub fn proper_time_to_decay(&self) -> f64 {
        self.proper_left
    }
    /// Apply an externally scheduled impulse while preserving the sampled
    /// proper lifetime. Energy is recomputed from this species' mass shell.
    pub fn set_momentum(&mut self, momentum: Vec3) -> Result<(), &'static str> {
        self.four = FourMomentum::on_shell(self.species.mass(), momentum)?;
        Ok(())
    }
    fn lab_time_to_decay(&self) -> f64 {
        if self.proper_left.is_finite() {
            self.proper_left * (self.four.energy / self.species.mass())
        } else {
            f64::INFINITY
        }
    }
    fn drift(&mut self, dt: f64) -> Result<(), &'static str> {
        let position = self.position + (self.four.momentum / self.four.energy) * C * dt;
        let time = self.time + dt;
        if !time.is_finite() || position.0.iter().any(|v| !v.is_finite()) {
            return Err("flight exceeds representable position/time");
        }
        self.position = position;
        self.time = time;
        self.proper_elapsed += dt * (self.species.mass() / self.four.energy);
        if self.proper_left.is_finite() {
            self.proper_left =
                (self.proper_left - dt * self.species.mass() / self.four.energy).max(0.0);
        }
        Ok(())
    }
    /// Relativistic Boris momentum push with uniform SI E [V/m], B [T].
    /// Operator-split force only: advance() performs drift/decay. Caller must
    /// resolve gyro motion and force/decay splitting. Returns field work [GeV].
    pub fn push_uniform(
        &mut self,
        electric: Vec3,
        magnetic: Vec3,
        dt: f64,
    ) -> Result<f64, &'static str> {
        if !dt.is_finite()
            || dt < 0.0
            || electric
                .0
                .iter()
                .chain(magnetic.0.iter())
                .any(|x| !x.is_finite())
        {
            return Err("invalid field/time step");
        }
        let q = self.species.charges().electric3 as f64 / 3.0 * ELEMENTARY_CHARGE;
        if q == 0.0 || dt == 0.0 {
            return Ok(0.0);
        }
        let impulse = electric * (q * C * dt / (2.0 * GEV_J));
        let minus = self.four.momentum + impulse;
        let energy = self.species.mass().hypot(minus.length());
        let t = magnetic * (q * C * C * dt / (2.0 * GEV_J * energy));
        let s = t * (2.0 / (1.0 + t.dot(t)));
        let prime = minus + minus.cross(t);
        let plus = minus + prime.cross(s);
        let updated = FourMomentum::on_shell(self.species.mass(), plus + impulse)?;
        let work = updated.energy - self.four.energy;
        self.four = updated;
        Ok(work)
    }
}
#[derive(Default)]
pub struct Flight {
    pub live: Vec<Particle>,
    pub unresolved: Vec<Unresolved>,
    pub decays: usize,
    /// Four-momentum supplied by prescribed external fields, in GeV.
    pub field_exchange: FourMomentum,
}
/// A cascade has finite depth in this closure (top->W->tau->mu->e).
/// Daughters move only for the remaining interval after their birth.
pub fn fly(particle: Particle, dt: f64, domain: Domain) -> Result<Flight, &'static str> {
    if !dt.is_finite() || dt < 0.0 {
        return Err("invalid flight interval");
    }
    if !particle.species.allowed(domain) {
        return Err("particle/domain mismatch");
    }
    let mut queue = vec![(particle, dt)];
    let mut flight = Flight::default();
    while let Some((mut p, remaining)) = queue.pop() {
        let lab_left = p.lab_time_to_decay();
        if lab_left > remaining || !lab_left.is_finite() {
            p.drift(remaining)?;
            flight.live.push(p);
            continue;
        }
        p.drift(lab_left)?;
        let mut rng = Rng(p.rng_state);
        let result = decay::sample(p.species, p.four, domain, &mut rng)?;
        flight.decays += 1;
        for (s, reason, four, charges) in result.unresolved {
            flight.unresolved.push(Unresolved {
                origin: s,
                reason,
                four,
                charges,
                position: p.position,
                time: p.time,
            });
        }
        for (s, four) in result.particles {
            let child = Particle::from_four(s, p.position, four, p.time, rng.next_u64())?;
            queue.push((child, (remaining - lab_left).max(0.0)));
        }
    }
    Ok(flight)
}
/// Covariant field flight with proper-time decay events and daughter orbits.
/// Radiation reaction/emission is absent. Uniform fields are externally supplied;
/// `field_exchange` records both their energy work and momentum impulse.
pub fn fly_uniform_field(
    particle: Particle,
    dt: f64,
    domain: Domain,
    field: UniformField,
) -> Result<Flight, &'static str> {
    field.validate()?;
    if !dt.is_finite() || dt < 0.0 || !particle.species.allowed(domain) {
        return Err("invalid field flight");
    }
    if field.electric == Vec3::default() && field.magnetic == Vec3::default() {
        return fly(particle, dt, domain);
    }
    let mut queue = vec![(particle, dt)];
    let mut flight = Flight::default();
    while let Some((mut p, remaining)) = queue.pop() {
        let q = p.species.charges().electric3 as f64 / 3.0;
        let decay;
        let used;
        if q == 0.0 {
            let left = p.lab_time_to_decay();
            decay = left <= remaining;
            used = left.min(remaining);
            p.drift(used)?;
        } else {
            let (orbit, event) =
                field.until(p.four, p.species.mass(), q, remaining, p.proper_left)?;
            decay = event;
            used = orbit.lab_time;
            let position = p.position + orbit.displacement;
            if position.0.iter().any(|x| !x.is_finite()) || !(p.time + used).is_finite() {
                return Err("unrepresentable field position/time");
            }
            flight.field_exchange = flight.field_exchange + (orbit.four - p.four);
            p.four = orbit.four;
            p.position = position;
            p.time += used;
            p.proper_elapsed += orbit.proper_time;
            if p.proper_left.is_finite() {
                p.proper_left = (p.proper_left - orbit.proper_time).max(0.0);
            }
        }
        if !decay {
            flight.live.push(p);
            continue;
        }
        let mut rng = Rng(p.rng_state);
        let result = decay::sample(p.species, p.four, domain, &mut rng)?;
        flight.decays += 1;
        for (s, reason, four, charges) in result.unresolved {
            flight.unresolved.push(Unresolved {
                origin: s,
                reason,
                four,
                charges,
                position: p.position,
                time: p.time,
            });
        }
        for (s, four) in result.particles {
            let child = Particle::from_four(s, p.position, four, p.time, rng.next_u64())?;
            queue.push((child, (remaining - used).max(0.0)));
        }
    }
    Ok(flight)
}
/// No all-pairs collision search. Call QED kernels at scheduled collisions;
/// medium collision rates and QCD require separate closures.
pub struct Simulation {
    pub(crate) domain: Domain,
    pub(crate) particles: Vec<Particle>,
    pub(crate) unresolved: Vec<Unresolved>,
    pub(crate) time: f64,
    pub(crate) decays: usize,
}
enum History {
    /// Most species/steps only drift: no queue or per-particle heap allocation.
    Ballistic(Particle),
    Cascade(Flight),
}
impl Simulation {
    pub fn new(domain: Domain) -> Self {
        Self {
            domain,
            particles: vec![],
            unresolved: vec![],
            time: 0.0,
            decays: 0,
        }
    }
    pub fn insert(&mut self, mut p: Particle) -> Result<(), &'static str> {
        if !p.species.allowed(self.domain) {
            return Err("particle/domain mismatch");
        }
        // Inject at current simulation time, retaining sampled proper lifetime.
        p.time = self.time;
        self.particles.push(p);
        Ok(())
    }
    pub fn particles(&self) -> &[Particle] {
        &self.particles
    }
    /// Mutable particle access for scheduled field/interaction sources.
    pub fn particles_mut(&mut self) -> &mut [Particle] {
        &mut self.particles
    }
    pub fn unresolved(&self) -> &[Unresolved] {
        &self.unresolved
    }
    /// Hand retained events to an external closure. The caller then owns their
    /// energy/charge ledger; totals() includes only events still retained here.
    pub fn take_unresolved(&mut self) -> Vec<Unresolved> {
        std::mem::take(&mut self.unresolved)
    }
    pub fn time(&self) -> f64 {
        self.time
    }
    pub fn decays(&self) -> usize {
        self.decays
    }
    /// Apply one caller-scheduled photon/lepton encounter at the current time.
    /// The caller supplies the spatial/rate closure. Proper decay clocks survive
    /// recoil, so time dilation follows the new momentum on the next flight.
    pub fn compton(&mut self, photon: usize, lepton: usize) -> Result<(), &'static str> {
        if photon == lepton || photon >= self.particles.len() || lepton >= self.particles.len() {
            return Err("invalid collision indices");
        }
        let a = &self.particles[photon];
        let b = &self.particles[lepton];
        if a.species != Species::Photon {
            return Err("first collision particle must be a photon");
        }
        let mut rng = Rng(a.rng_state);
        let out = crate::qed::compton(a.four, b.four, b.species, &mut rng)?;
        self.particles[photon].four = out[0];
        self.particles[photon].rng_state = rng.0;
        self.particles[lepton].four = out[1];
        Ok(())
    }
    /// Replace a caller-scheduled comoving e-/e+ encounter with two photons.
    /// Uses the electron position as the scheduler's collision point. This
    /// threshold kernel rejects pairs with relative centre-of-momentum motion.
    pub fn annihilate_at_threshold(
        &mut self,
        electron: usize,
        positron: usize,
    ) -> Result<(), &'static str> {
        if electron == positron
            || electron >= self.particles.len()
            || positron >= self.particles.len()
        {
            return Err("invalid collision indices");
        }
        let a = &self.particles[electron];
        let b = &self.particles[positron];
        if a.species != Species::Electron || b.species != Species::Positron {
            return Err("annihilation requires electron and positron");
        }
        let position = a.position;
        let mut rng = Rng(a.rng_state);
        let out = crate::qed::annihilate_at_threshold(a.four, b.four, &mut rng)?;
        let mut children = Vec::with_capacity(2);
        for four in out {
            children.push(Particle::from_four(
                Species::Photon,
                position,
                four,
                self.time,
                rng.next_u64(),
            )?);
        }
        self.particles.remove(electron.max(positron));
        self.particles.remove(electron.min(positron));
        self.particles.extend(children);
        Ok(())
    }
    fn collision_indices(&self, a: usize, b: usize) -> Result<(), &'static str> {
        if a == b || a >= self.particles.len() || b >= self.particles.len() {
            return Err("invalid collision indices");
        }
        Ok(())
    }
    fn replace_pair(
        &mut self,
        a: usize,
        b: usize,
        species: [Species; 2],
        out: [FourMomentum; 2],
        rng: &mut Rng,
    ) -> Result<(), &'static str> {
        let position = self.particles[a].position;
        let children = [
            Particle::from_four(species[0], position, out[0], self.time, rng.next_u64())?,
            Particle::from_four(species[1], position, out[1], self.time, rng.next_u64())?,
        ];
        self.particles.remove(a.max(b));
        self.particles.remove(a.min(b));
        self.particles.extend(children);
        Ok(())
    }
    /// Scheduled in-flight e-/e+ annihilation with the finite-mass QED kernel.
    /// New photons are born at the electron's scheduler-selected position/time.
    pub fn annihilate(&mut self, electron: usize, positron: usize) -> Result<(), &'static str> {
        self.collision_indices(electron, positron)?;
        let a = &self.particles[electron];
        let b = &self.particles[positron];
        if a.species != Species::Electron || b.species != Species::Positron {
            return Err("annihilation requires electron and positron");
        }
        let mut rng = Rng(a.rng_state);
        let out = crate::qed::annihilate(a.four, b.four, &mut rng)?;
        self.replace_pair(electron, positron, [Species::Photon; 2], out, &mut rng)
    }
    /// Scheduled real-photon pair production into e, mu or tau and its anti.
    /// The flavor input must be the negatively charged lepton. Its proper
    /// lifetime is sampled at birth; the first photon supplies event position.
    pub fn pair_produce(
        &mut self,
        first: usize,
        second: usize,
        flavor: Species,
    ) -> Result<(), &'static str> {
        self.collision_indices(first, second)?;
        if !matches!(flavor, Species::Electron | Species::Muon | Species::Tau) {
            return Err("pair flavor must be a negative charged lepton");
        }
        let a = &self.particles[first];
        let b = &self.particles[second];
        if a.species != Species::Photon || b.species != Species::Photon {
            return Err("pair production requires two photons");
        }
        let mut rng = Rng(a.rng_state);
        let out = crate::qed::breit_wheeler(a.four, b.four, flavor, &mut rng)?;
        self.replace_pair(first, second, [flavor, flavor.conjugate()], out, &mut rng)
    }
    /// Scheduled elastic neutrino-electron scattering in the Fermi domain.
    /// Caller supplies flavor-appropriate chiral coefficients and event rates.
    /// Both identities and proper decay clocks survive recoil.
    pub fn scatter_neutrino(
        &mut self,
        neutrino: usize,
        electron: usize,
        couplings: crate::weak::Couplings,
    ) -> Result<(), &'static str> {
        self.collision_indices(neutrino, electron)?;
        let a = &self.particles[neutrino];
        let b = &self.particles[electron];
        if !a.species.neutrino() || b.species != Species::Electron {
            return Err("weak elastic encounter requires neutrino and electron");
        }
        let mut rng = Rng(a.rng_state);
        let out = crate::weak::scatter(a.four, b.four, couplings, &mut rng)?;
        self.particles[neutrino].four = out[0];
        self.particles[neutrino].rng_state = rng.0;
        self.particles[electron].four = out[1];
        Ok(())
    }
    pub fn totals(&self) -> (FourMomentum, Charges) {
        let mut p = FourMomentum::default();
        let mut q = Charges::default();
        for v in &self.particles {
            p = p + v.four;
            q = q + v.species.charges();
        }
        for v in &self.unresolved {
            p = p + v.four;
            q = q + v.charges;
        }
        (p, q)
    }
    /// Integrate trajectories and decay clocks together in a prescribed uniform
    /// SI field. Returns the external field's supplied total four-momentum.
    /// No kick/drift splitting and no gyro-step truncation for constant fields.
    pub fn advance_uniform_field(
        &mut self,
        dt: f64,
        field: UniformField,
    ) -> Result<FourMomentum, &'static str> {
        field.validate()?;
        if !dt.is_finite() || dt < 0.0 || !(self.time + dt).is_finite() {
            return Err("invalid field simulation interval");
        }
        if field.electric == Vec3::default() && field.magnetic == Vec3::default() {
            self.advance(dt)?;
            return Ok(FourMomentum::default());
        }
        if dt == 0.0 {
            return Ok(FourMomentum::default());
        }
        let history = |p: &Particle| fly_uniform_field(p.clone(), dt, self.domain, field);
        let flights: Result<Vec<_>, _> = if self.particles.len() >= 1024 {
            self.particles.par_iter().map(history).collect()
        } else {
            self.particles.iter().map(history).collect()
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
    /// Compute histories before committing; rejected input leaves state intact.
    /// Cloning each small particle is O(n), with no pairwise search or tiny
    /// lifetime substeps. Results reduce in input order across thread counts.
    pub fn advance(&mut self, dt: f64) -> Result<(), &'static str> {
        if !dt.is_finite() || dt < 0.0 || !(self.time + dt).is_finite() {
            return Err("invalid simulation interval");
        }
        if dt == 0.0 {
            return Ok(());
        }
        let history = |p: &Particle| -> Result<History, &'static str> {
            if !p.species.allowed(self.domain) {
                return Err("particle/domain mismatch");
            }
            if p.lab_time_to_decay() > dt {
                let mut live = p.clone();
                live.drift(dt)?;
                Ok(History::Ballistic(live))
            } else {
                Ok(History::Cascade(fly(p.clone(), dt, self.domain)?))
            }
        };
        let flights: Result<Vec<_>, _> = if self.particles.len() >= 1024 {
            self.particles.par_iter().map(history).collect()
        } else {
            self.particles.iter().map(history).collect()
        };
        let flights = flights?;
        self.particles.clear();
        for f in flights {
            match f {
                History::Ballistic(p) => self.particles.push(p),
                History::Cascade(f) => {
                    self.particles.extend(f.live);
                    self.unresolved.extend(f.unresolved);
                    self.decays += f.decays;
                }
            }
        }
        self.time += dt;
        Ok(())
    }
}
