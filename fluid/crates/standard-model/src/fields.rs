//! Covariant constant-field Lorentz force integrated in proper time.
//! Test-particle external fields: no radiation reaction or quantum emission.
use crate::{FourMomentum, Vec3, C, ELEMENTARY_CHARGE, GEV_J, HBAR_GEV_S};

#[derive(Clone, Copy, Debug, Default)]
pub struct UniformField {
    pub electric: Vec3,
    pub magnetic: Vec3,
}
#[derive(Clone, Copy, Debug)]
pub struct Orbit {
    pub four: FourMomentum,
    pub displacement: Vec3,
    pub lab_time: f64,
    pub proper_time: f64,
}
impl UniformField {
    pub fn validate(self) -> Result<(), &'static str> {
        if self
            .electric
            .0
            .iter()
            .chain(self.magnetic.0.iter())
            .any(|x| !x.is_finite())
        {
            return Err("nonfinite electromagnetic field");
        }
        Ok(())
    }
    /// Fields in the rest frame of the timelike reference (SI units).
    pub fn to_rest(self, reference: FourMomentum, mass: f64) -> Result<Self, &'static str> {
        self.validate()?;
        reference.validate_on_shell(mass)?;
        if mass <= 0.0 {
            return Err("field frame must be timelike");
        }
        let u = reference.momentum / mass;
        let gamma = reference.energy / mass;
        let speed = crate::kinematics::norm(reference.momentum);
        if speed == 0.0 {
            return Ok(self);
        }
        let axis = reference.momentum / speed;
        let parallel_e = axis * axis.dot(self.electric);
        let parallel_b = axis * axis.dot(self.magnetic);
        // The parallel field is invariant. Project before multiplying by
        // gamma, avoiding cancellation of two huge longitudinal terms.
        Ok(Self {
            electric: parallel_e
                + (self.electric - parallel_e) * gamma
                + u.cross(self.magnetic) * C,
            magnetic: parallel_b + (self.magnetic - parallel_b) * gamma
                - u.cross(self.electric) / C,
        })
    }
    /// Quantum nonlinearity chi using the rest-frame electric field. Classical
    /// orbits omit emission/reaction even at chi<<1; chi is a domain indicator,
    /// not a guarantee that radiative energy loss is negligible.
    pub fn quantum_chi(
        self,
        p: FourMomentum,
        mass: f64,
        charge_e: f64,
    ) -> Result<f64, &'static str> {
        if !charge_e.is_finite() {
            return Err("invalid charge");
        }
        let rest = self.to_rest(p, mass)?;
        Ok(charge_e.abs()
            * ELEMENTARY_CHARGE
            * HBAR_GEV_S
            * C
            * crate::kinematics::norm(rest.electric)
            / (mass * mass * GEV_J))
    }
    /// Exact constant-field solution up to floating-point matrix exponential
    /// error. Integrates the four-momentum and its integral in the same update.
    pub fn propagate(
        self,
        p: FourMomentum,
        mass: f64,
        charge_e: f64,
        tau: f64,
    ) -> Result<Orbit, &'static str> {
        self.validate()?;
        p.validate_on_shell(mass)?;
        if mass <= 0.0 || !charge_e.is_finite() || !tau.is_finite() || tau < 0.0 {
            return Err("invalid massive proper-time orbit inputs");
        }
        if charge_e == 0.0 || (self.electric == Vec3::default() && self.magnetic == Vec3::default())
        {
            return checked(p, p.momentum * (C * tau / mass), p.energy * tau / mass, tau);
        }
        let q = charge_e * ELEMENTARY_CHARGE;
        let en = crate::kinematics::norm(self.electric);
        let bn = crate::kinematics::norm(self.magnetic);
        if en == 0.0 {
            let axis = self.magnetic / bn;
            let parallel = axis * axis.dot(p.momentum);
            let perp = p.momentum - parallel;
            let theta = q * C * C * bn * tau / (mass * GEV_J);
            let (sin, cos) = theta.sin_cos();
            let cross = perp.cross(axis);
            let momentum = parallel + perp * cos + cross * sin;
            let four = FourMomentum::on_shell(mass, momentum)?;
            let displacement =
                (parallel + perp * sinc(theta) + cross * cosc(theta)) * (C * tau / mass);
            return checked(four, displacement, p.energy * tau / mass, tau);
        }
        if bn == 0.0 {
            let axis = self.electric / en;
            let longitudinal = axis.dot(p.momentum);
            let perp = p.momentum - axis * longitudinal;
            let theta = q * C * en * tau / (mass * GEV_J);
            if theta.abs() > 700.0 {
                return Err("unrepresentable electric rapidity");
            }
            let transverse_mass2 = mass * mass + perp.dot(perp);
            let (plus, minus) = if longitudinal >= 0.0 {
                let plus = p.energy + longitudinal;
                (plus, transverse_mass2 / plus)
            } else {
                let minus = p.energy - longitudinal;
                (transverse_mass2 / minus, minus)
            };
            let (boost_plus, boost_minus) = (plus * theta.exp(), minus * (-theta).exp());
            let energy = 0.5 * (boost_plus + boost_minus);
            let momentum = perp + axis * (0.5 * (boost_plus - boost_minus));
            let integral_e = 0.5 * tau * (plus * exprel(theta) + minus * exprel(-theta));
            let integral_p =
                perp * tau + axis * (0.5 * tau * (plus * exprel(theta) - minus * exprel(-theta)));
            return checked(
                FourMomentum { energy, momentum },
                integral_p * (C / mass),
                integral_e / mass,
                tau,
            );
        }
        let e = self.electric * (q * C / (mass * GEV_J));
        let b = self.magnetic * (q * C * C / (mass * GEV_J));
        let k = [
            [0.0, e.0[0], e.0[1], e.0[2]],
            [e.0[0], 0.0, b.0[2], -b.0[1]],
            [e.0[1], -b.0[2], 0.0, b.0[0]],
            [e.0[2], b.0[1], -b.0[0], 0.0],
        ];
        let (exponential, integral) = exp_integral(k, tau)?;
        let v = [p.energy, p.momentum.0[0], p.momentum.0[1], p.momentum.0[2]];
        let out = apply(exponential, v);
        let ix = apply(integral, v);
        checked(
            FourMomentum {
                energy: out[0],
                momentum: Vec3([out[1], out[2], out[3]]),
            },
            Vec3([ix[1], ix[2], ix[3]]) * (C / mass),
            ix[0] / mass,
            tau,
        )
    }
    /// Find the proper-time endpoint for a laboratory interval, or stop at the
    /// sampled proper decay clock. t(tau) is monotone since gamma>=1. Bounded
    /// growth bracketing avoids evaluating exp(K*dt) at an enormous rapidity.
    pub fn until(
        self,
        p: FourMomentum,
        mass: f64,
        charge: f64,
        dt: f64,
        proper_limit: f64,
    ) -> Result<(Orbit, bool), &'static str> {
        self.validate()?;
        p.validate_on_shell(mass)?;
        if mass <= 0.0
            || !dt.is_finite()
            || dt < 0.0
            || proper_limit.is_nan()
            || proper_limit < 0.0
            || !charge.is_finite()
        {
            return Err("invalid laboratory orbit interval");
        }
        if dt == 0.0 {
            return Ok((self.propagate(p, mass, charge, 0.0)?, proper_limit == 0.0));
        }
        // Constant gamma: solve directly before any growth bracketing, even
        // for steps spanning many gyro periods.
        if self.electric == Vec3::default() || charge == 0.0 {
            let tau = dt * mass / p.energy;
            let event = proper_limit <= tau;
            return Ok((
                self.propagate(p, mass, charge, tau.min(proper_limit))?,
                event,
            ));
        }
        let max_tau = dt.min(proper_limit);
        let rate = charge.abs() * ELEMENTARY_CHARGE * C / (mass * GEV_J)
            * (crate::kinematics::norm(self.electric) + C * crate::kinematics::norm(self.magnetic));
        let seed = (dt * mass / p.energy).min(if rate > 0.0 { 1.0 / rate } else { dt });
        let mut hi = max_tau.min(seed);
        let mut lo = 0.0;
        if hi == 0.0 && max_tau > 0.0 {
            return Err("unrepresentable proper-time interval");
        }
        let mut upper = self.propagate(p, mass, charge, hi)?;
        while upper.lab_time < dt && hi < max_tau {
            lo = hi;
            hi = (2.0 * hi).min(max_tau);
            upper = self.propagate(p, mass, charge, hi)?;
        }
        if hi == proper_limit && upper.lab_time <= dt {
            return Ok((upper, true));
        }
        let mut guess = (lo + hi) * 0.5;
        for _ in 0..80 {
            let trial = self.propagate(p, mass, charge, guess)?;
            let residual = trial.lab_time - dt;
            if residual.abs() <= 4e-13 * dt {
                return Ok((trial, false));
            }
            if residual > 0.0 {
                hi = guess;
            } else {
                lo = guess;
            }
            // Safeguarded Newton: derivative dt/dtau=E/m, from the same orbit.
            let next = guess - residual * mass / trial.four.energy;
            guess = if next > lo && next < hi {
                next
            } else {
                (lo + hi) * 0.5
            };
        }
        Err("proper-time orbit endpoint did not converge")
    }
}
fn checked(
    four: FourMomentum,
    displacement: Vec3,
    lab_time: f64,
    proper_time: f64,
) -> Result<Orbit, &'static str> {
    if !four.energy.is_finite()
        || four.energy <= 0.0
        || !lab_time.is_finite()
        || lab_time < 0.0
        || four
            .momentum
            .0
            .iter()
            .chain(displacement.0.iter())
            .any(|x| !x.is_finite())
    {
        return Err("unrepresentable field orbit");
    }
    Ok(Orbit {
        four,
        displacement,
        lab_time,
        proper_time,
    })
}
fn sinc(x: f64) -> f64 {
    if x.abs() < 1e-4 {
        1.0 - x * x / 6.0 + x.powi(4) / 120.0
    } else {
        x.sin() / x
    }
}
fn cosc(x: f64) -> f64 {
    if x.abs() < 1e-4 {
        x / 2.0 - x.powi(3) / 24.0 + x.powi(5) / 720.0
    } else {
        2.0 * (x / 2.0).sin().powi(2) / x
    }
}
fn exprel(x: f64) -> f64 {
    if x.abs() < 1e-5 {
        1.0 + x / 2.0 + x * x / 6.0 + x.powi(3) / 24.0 + x.powi(4) / 120.0
    } else {
        x.exp_m1() / x
    }
}
type Matrix = [[f64; 4]; 4];
fn identity() -> Matrix {
    std::array::from_fn(|i| std::array::from_fn(|j| if i == j { 1.0 } else { 0.0 }))
}
fn product(a: Matrix, b: Matrix) -> Matrix {
    std::array::from_fn(|i| std::array::from_fn(|j| (0..4).map(|k| a[i][k] * b[k][j]).sum()))
}
fn apply(a: Matrix, v: [f64; 4]) -> [f64; 4] {
    std::array::from_fn(|i| (0..4).map(|j| a[i][j] * v[j]).sum())
}
fn exp_integral(k: Matrix, tau: f64) -> Result<(Matrix, Matrix), &'static str> {
    let size = k
        .iter()
        .map(|r| r.iter().map(|v| v.abs() * tau).sum::<f64>())
        .fold(0.0, f64::max);
    if !size.is_finite() {
        return Err("unrepresentable field generator");
    }
    let scale = if size > 0.25 {
        (size / 0.25).log2().ceil() as i32
    } else {
        0
    };
    if scale > 60 {
        return Err("field orbit requires excessive matrix scaling");
    }
    let h = tau / 2f64.powi(scale);
    let a = k.map(|r| r.map(|v| v * h));
    let mut exponential = identity();
    let mut integral = identity();
    let mut term = identity();
    for n in 1..=20 {
        term = product(term, a).map(|r| r.map(|v| v / n as f64));
        for i in 0..4 {
            for j in 0..4 {
                exponential[i][j] += term[i][j];
                integral[i][j] += term[i][j] / (n + 1) as f64;
            }
        }
    }
    integral = integral.map(|r| r.map(|v| v * h));
    for _ in 0..scale {
        let extra = product(exponential, integral);
        for i in 0..4 {
            for j in 0..4 {
                integral[i][j] += extra[i][j];
            }
        }
        exponential = product(exponential, exponential);
    }
    if exponential
        .iter()
        .chain(integral.iter())
        .flatten()
        .any(|v| !v.is_finite())
    {
        return Err("unrepresentable field exponential");
    }
    Ok((exponential, integral))
}
