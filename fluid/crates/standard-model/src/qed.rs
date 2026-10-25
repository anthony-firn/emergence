//! Unpolarized tree-level free-lepton Compton and photon/pair reactions. No bound electron,
//! radiative, strong-field or medium corrections. Cross sections [m²].
use crate::{
    kinematics::{isotropic, two_body},
    radiation_math, FourMomentum, Rng, Species, ALPHA, HBAR_C,
};

fn lepton(s: Species) -> bool {
    matches!(s.pdg_id().abs(), 11 | 13 | 15)
}
/// Total Klein-Nishina cross section for a photon energy in the lepton rest frame.
pub fn compton_cross_section(species: Species, rest_energy: f64) -> Result<f64, &'static str> {
    if !lepton(species) || rest_energy < 0.0 || !rest_energy.is_finite() {
        return Err("Compton target must be a charged lepton and energy finite/nonnegative");
    }
    let x = rest_energy / species.mass();
    if !x.is_finite() {
        return Err("unrepresentable dimensionless photon energy");
    }
    let r = ALPHA * HBAR_C / species.mass();
    let thomson = 8.0 * std::f64::consts::PI * r * r / 3.0;
    if x < 1e-3 {
        return Ok(thomson * (1.0 - 2.0 * x + 5.2 * x * x - 13.3 * x * x * x));
    }
    if x > 1e3 {
        // Same exact expression in inverse x; no overflowing x³ or 2x.
        let inv = 1.0 / x;
        let log = (x + 0.5).ln() + std::f64::consts::LN_2;
        let value = inv
            * ((1.0 + inv) * ((1.0 + inv) / (1.0 + 0.5 * inv) - inv * log) + 0.5 * log
                - (3.0 + inv) / (2.0 + inv).powi(2));
        return Ok(2.0 * std::f64::consts::PI * r * r * value);
    }
    let log = (2.0 * x).ln_1p();
    let value = (1.0 + x) / x.powi(3) * (2.0 * x * (1.0 + x) / (1.0 + 2.0 * x) - log)
        + log / (2.0 * x)
        - (1.0 + 3.0 * x) / (1.0 + 2.0 * x).powi(2);
    Ok(2.0 * std::f64::consts::PI * r * r * value)
}
/// Returns outgoing photon and lepton. Input target may be relativistically
/// moving; the angular kernel is sampled in its rest frame, then boosted.
pub fn compton(
    photon: FourMomentum,
    target: FourMomentum,
    species: Species,
    rng: &mut Rng,
) -> Result<[FourMomentum; 2], &'static str> {
    if !lepton(species) {
        return Err("Compton scattering requires a charged lepton");
    }
    validate(photon, 0.0)?;
    validate(target, species.mass())?;
    let m = species.mass();
    let incoming = photon.to_rest_on_shell(0.0, target, m);
    let x = incoming.energy / m;
    if !x.is_finite() || x <= 0.0 {
        return Err("nonpositive rest-frame photon energy");
    }
    let (delta, fraction) = if x < 0.1 {
        loop {
            let c = 2.0 * rng.uniform() - 1.0;
            let y = 1.0 / (1.0 + x * (1.0 - c));
            let w = 0.5 * y * y * (y + 1.0 / y - 1.0 + c * c);
            if rng.uniform() <= w {
                break (1.0 - c, y);
            }
        }
    } else {
        // Mixture envelope y + 1/y in y=E_out/E_in. Unlike uniform angle
        // rejection this remains efficient for high-energy forward scattering.
        let ymin = 0.5 / (x + 0.5);
        let a = 0.5 * (1.0 - ymin * ymin);
        let b = (x + 0.5).ln() + std::f64::consts::LN_2;
        loop {
            let y = if rng.uniform() * (a + b) < a {
                (ymin * ymin + (1.0 - ymin * ymin) * rng.uniform()).sqrt()
            } else {
                (-b * rng.uniform()).exp()
            };
            let delta = ((1.0 - y) / (x * y)).clamp(0.0, 2.0);
            if rng.uniform() <= 1.0 - delta * (2.0 - delta) / (y + 1.0 / y) {
                break (delta, y);
            }
        }
    };
    // Preserve the sampled energy fraction and the small scattering angle.
    // Reconstructing y from a rounded cos(theta) erases recoil at large x.
    let energy = incoming.energy * fraction;
    let axis = incoming.momentum / crate::kinematics::norm(incoming.momentum);
    let basis = if axis.0[2].abs() < 0.9 {
        crate::Vec3([0.0, 0.0, 1.0])
    } else {
        crate::Vec3([1.0, 0.0, 0.0])
    };
    let u = axis.cross(basis).unit();
    let v = axis.cross(u);
    let phi = 2.0 * std::f64::consts::PI * rng.uniform();
    let dir =
        axis * (1.0 - delta) + (u * phi.cos() + v * phi.sin()) * (delta * (2.0 - delta)).sqrt();
    let outgoing = FourMomentum {
        energy,
        momentum: dir * energy,
    };
    let recoil_p = incoming.momentum - outgoing.momentum;
    // Recompute on shell to retain a tiny recoil energy without cancellation.
    let recoil = FourMomentum::on_shell(m, recoil_p)?;
    Ok([
        outgoing.from_rest_on_shell(0.0, target, m),
        recoil.from_rest_on_shell(m, target, m),
    ])
}
fn validate(p: FourMomentum, mass: f64) -> Result<(), &'static str> {
    p.validate_on_shell(mass)
}
/// Two-photon annihilation of a comoving e-/e+ pair at threshold (possibly
/// boosted). Above threshold use a differential annihilation closure, not an
/// isotropic replacement. Returns photons; caller schedules the encounter.
pub fn annihilate_at_threshold(
    electron: FourMomentum,
    positron: FourMomentum,
    rng: &mut Rng,
) -> Result<[FourMomentum; 2], &'static str> {
    let m = Species::Electron.mass();
    validate(electron, m)?;
    validate(positron, m)?;
    let total = electron + positron;
    let rest_e = electron.to_rest_on_shell(m, total, 2.0 * m);
    if rest_e.momentum.length() > 1e-7 * m || (rest_e.energy - m).abs() > 1e-10 * m {
        return Err("annihilation kernel is valid only for comoving threshold pairs");
    }
    two_body(total, 2.0 * m, 0.0, 0.0, rng)
}
/// Generic elastic phase-space kinematics, with an externally calibrated
/// angular sample cos(theta). Does not prescribe a rate or a strong interaction.
pub fn elastic(
    a: FourMomentum,
    b: FourMomentum,
    ma: f64,
    mb: f64,
    cos_theta: f64,
    rng: &mut Rng,
) -> Result<[FourMomentum; 2], &'static str> {
    validate(a, ma)?;
    validate(b, mb)?;
    if !cos_theta.is_finite() || !(-1.0..=1.0).contains(&cos_theta) {
        return Err("invalid scattering angle");
    }
    let (total, mass) = collision_frame(a, ma, b, mb)?;
    let mut pair = two_body(total, mass, ma, mb, rng)?;
    let ar = a.to_rest_on_shell(ma, total, mass);
    let out = pair[0].to_rest_on_shell(ma, total, mass);
    let axis = if ar.momentum.length() > 0.0 {
        ar.momentum.unit()
    } else {
        isotropic(rng)
    };
    let p = radiation_math::around(axis, cos_theta, 2.0 * std::f64::consts::PI * rng.uniform())
        * out.momentum.length();
    pair[0] = FourMomentum {
        energy: ma.hypot(p.length()),
        momentum: p,
    }
    .from_rest_on_shell(ma, total, mass);
    pair[1] = FourMomentum {
        energy: mb.hypot(p.length()),
        momentum: -p,
    }
    .from_rest_on_shell(mb, total, mass);
    Ok(pair)
}

/// Unpolarized finite-mass Born kernel for gamma gamma <-> l- l+.
/// mu is the negative lepton angle relative to the first photon in the CM.
/// The two-photon annihilation event distribution includes an identical-final
/// factor of 1/2 over the full sphere; pair production has distinct final states.
pub fn pair_angular_kernel(s: f64, mass: f64, mu: f64) -> Result<f64, &'static str> {
    if !s.is_finite()
        || !mass.is_finite()
        || mass <= 0.0
        || s < 4.0 * mass * mass
        || !mu.is_finite()
        || mu.abs() > 1.0
    {
        return Err("invalid pair invariant/angle");
    }
    let r = 4.0 * mass * mass / s;
    let beta2 = (1.0 - r).max(0.0);
    let d = r + beta2 * (1.0 - mu) * (1.0 + mu);
    let ratio = r / d;
    Ok(2.0 * (1.0 + beta2 * mu * mu) / d + 4.0 * ratio - 4.0 * ratio * ratio)
}
/// Breit-Wheeler total event cross section [m²] for one charged-lepton flavor.
/// s [GeV²]. Below and at threshold the production rate is zero.
pub fn pair_cross_section(species: Species, s: f64) -> Result<f64, &'static str> {
    if !lepton(species) || !s.is_finite() || s < 0.0 {
        return Err("invalid pair-production species/invariant");
    }
    let m = species.mass();
    if s <= 4.0 * m * m {
        return Ok(0.0);
    }
    let r = 4.0 * m * m / s;
    let beta = (1.0 - r).sqrt();
    let bracket = pair_bracket(beta, r);
    Ok(4.0 * std::f64::consts::PI * ALPHA * ALPHA * HBAR_C * HBAR_C / s * bracket)
}
fn pair_bracket(beta: f64, r: f64) -> f64 {
    if beta < 1e-3 {
        // Avoid cancellation in the threshold logarithm.
        beta * (1.0 + 2.0 * beta * beta - 0.4 * beta.powi(4) + 2.0 / 21.0 * beta.powi(6))
    } else {
        let log = 2.0 * (1.0 + beta).ln() - r.ln();
        0.5 * (3.0 - beta.powi(4)) * log - beta * (2.0 - beta * beta)
    }
}
/// Finite-mass two-photon annihilation cross section [m²]. It diverges as
/// 1/beta at threshold; annihilation_sigma_v has the finite physical limit.
pub fn annihilation_cross_section(s: f64) -> Result<f64, &'static str> {
    let m = Species::Electron.mass();
    if !s.is_finite() || s < 4.0 * m * m {
        return Err("invalid annihilation invariant");
    }
    if s == 4.0 * m * m {
        return Ok(f64::INFINITY);
    }
    let beta2 = 1.0 - 4.0 * m * m / s;
    Ok(pair_cross_section(Species::Electron, s)? / (2.0 * beta2))
}
/// sigma times CM Møller flux speed 2 beta c [m³/s], finite at threshold.
/// The flux speed is a collision-rate factor, not a particle velocity.
pub fn annihilation_sigma_v(s: f64) -> Result<f64, &'static str> {
    let m = Species::Electron.mass();
    if !s.is_finite() || s < 4.0 * m * m {
        return Err("invalid annihilation invariant");
    }
    let r = 4.0 * m * m / s;
    let beta = (1.0 - r).sqrt();
    let ratio = if beta == 0.0 {
        1.0
    } else {
        pair_bracket(beta, r) / beta
    };
    Ok(4.0 * std::f64::consts::PI * ALPHA * ALPHA * HBAR_C * HBAR_C * crate::C / s * ratio)
}
fn sample_pair_angle(s: f64, mass: f64, rng: &mut Rng) -> Result<f64, &'static str> {
    let r = 4.0 * mass * mass / s;
    let beta = (1.0 - r).max(0.0).sqrt();
    if r < 1e-14 {
        return Err("CM collinear peak exceeds double-precision angular resolution");
    }
    let rapidity = (1.0 + beta).ln() - 0.5 * r.ln();
    loop {
        // Inverse CDF of 1/(1-beta² mu²); no inefficient uniform-angle rejection.
        let mu = if beta < 1e-5 {
            2.0 * rng.uniform() - 1.0
        } else {
            (((2.0 * rng.uniform() - 1.0) * rapidity).tanh() / beta).clamp(-1.0, 1.0)
        };
        let d = r + (1.0 - r) * (1.0 - mu) * (1.0 + mu);
        let weight = pair_angular_kernel(s, mass, mu)? * d / (4.0 + 4.0 * r);
        if rng.uniform() <= weight {
            return Ok(mu);
        }
    }
}
fn collision_frame(
    a: FourMomentum,
    ma: f64,
    b: FourMomentum,
    mb: f64,
) -> Result<(FourMomentum, f64), &'static str> {
    let mass = a.pair_mass(ma, b, mb)?;
    if mass <= 0.0 {
        return Err("collision has no timelike CM frame");
    }
    let total = FourMomentum::on_shell(mass, a.momentum + b.momentum)?;
    Ok((total, mass))
}
fn along(
    total: FourMomentum,
    mass: f64,
    ma: f64,
    mb: f64,
    dir: crate::Vec3,
) -> Result<[FourMomentum; 2], &'static str> {
    if mass < ma + mb {
        return Err("closed final-state threshold");
    }
    let p = 0.5
        * mass
        * ((1.0 - ((ma + mb) / mass).powi(2)) * (1.0 - ((ma - mb) / mass).powi(2)))
            .max(0.0)
            .sqrt();
    let a = FourMomentum {
        energy: ma.hypot(p),
        momentum: dir * p,
    };
    let b = FourMomentum {
        energy: mb.hypot(p),
        momentum: -dir * p,
    };
    Ok([
        a.from_rest_on_shell(ma, total, mass),
        b.from_rest_on_shell(mb, total, mass),
    ])
}
/// In-flight e-/e+ -> two photons with the finite-electron-mass angular kernel.
/// Input may have arbitrary relative motion and a common Lorentz boost.
pub fn annihilate(
    electron: FourMomentum,
    positron: FourMomentum,
    rng: &mut Rng,
) -> Result<[FourMomentum; 2], &'static str> {
    let m = Species::Electron.mass();
    let (total, mass) = collision_frame(electron, m, positron, m)?;
    let incoming = electron.to_rest_on_shell(m, total, mass);
    let axis = if crate::kinematics::norm(incoming.momentum) > 1e-12 * m {
        incoming.momentum.unit()
    } else {
        isotropic(rng)
    };
    let mu = sample_pair_angle(mass * mass, m, rng)?;
    let dir = radiation_math::around(axis, mu, 2.0 * std::f64::consts::PI * rng.uniform());
    along(total, mass, 0.0, 0.0, dir)
}
/// Two real photons -> charged-lepton pair with finite-mass Born angular law.
/// Output order is negative lepton then positive lepton; threshold is invariant.
pub fn breit_wheeler(
    first: FourMomentum,
    second: FourMomentum,
    species: Species,
    rng: &mut Rng,
) -> Result<[FourMomentum; 2], &'static str> {
    if !lepton(species) {
        return Err("pair kernel requires a charged lepton flavor");
    }
    let m = species.mass();
    let (total, mass) = collision_frame(first, 0.0, second, 0.0)?;
    if mass <= 2.0 * m {
        return Err("photon pair is at or below production threshold");
    }
    let incoming = first.to_rest_on_shell(0.0, total, mass);
    let mu = sample_pair_angle(mass * mass, m, rng)?;
    let dir = radiation_math::around(
        incoming.momentum.unit(),
        mu,
        2.0 * std::f64::consts::PI * rng.uniform(),
    );
    along(total, mass, m, m, dir)
}
