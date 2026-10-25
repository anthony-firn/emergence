use crate::{Rng, Vec3};
/// Scaled norm avoids intermediate squares overflowing or underflowing.
pub fn norm(v: Vec3) -> f64 {
    v.0[0].hypot(v.0[1]).hypot(v.0[2])
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FourMomentum {
    pub energy: f64,
    pub momentum: Vec3,
}
impl std::ops::Add for FourMomentum {
    type Output = Self;
    fn add(self, b: Self) -> Self {
        Self {
            energy: self.energy + b.energy,
            momentum: self.momentum + b.momentum,
        }
    }
}
impl std::ops::Sub for FourMomentum {
    type Output = Self;
    fn sub(self, b: Self) -> Self {
        Self {
            energy: self.energy - b.energy,
            momentum: self.momentum - b.momentum,
        }
    }
}
impl FourMomentum {
    pub fn on_shell(mass: f64, momentum: Vec3) -> Result<Self, &'static str> {
        if !mass.is_finite() || mass < 0.0 || momentum.0.iter().any(|x| !x.is_finite()) {
            return Err("invalid mass or momentum");
        }
        let energy = mass.hypot(norm(momentum));
        if !energy.is_finite() || energy <= 0.0 {
            return Err("zero or nonfinite particle energy");
        }
        Ok(Self { energy, momentum })
    }
    pub fn validate_on_shell(self, mass: f64) -> Result<(), &'static str> {
        let expected = Self::on_shell(mass, self.momentum)?;
        if !self.energy.is_finite()
            || self.energy <= 0.0
            || (self.energy - expected.energy).abs() > 1e-10 * self.energy
        {
            return Err("invalid or off-shell four-momentum");
        }
        Ok(())
    }
    /// Stable on-shell Minkowski product for nearly parallel ultrarelativistic
    /// momenta. Known masses retain information that E²-|p|² cannot represent.
    pub fn dot_on_shell(self, mass: f64, b: Self, bmass: f64) -> Result<f64, &'static str> {
        self.validate_on_shell(mass)?;
        b.validate_on_shell(bmass)?;
        let pa = norm(self.momentum);
        let pb = norm(b.momentum);
        if pa == 0.0 || pb == 0.0 {
            return Ok(self.energy * b.energy);
        }
        let sa = self.energy + pa;
        let sb = b.energy + pb;
        let na = self.momentum / pa;
        let nb = b.momentum / pb;
        let cosine = na.dot(nb);
        let angle = if cosine >= 0.0 {
            // sin²(theta)/(1+cos(theta)) avoids interpreting normalization
            // roundoff as a tiny longitudinal angular separation.
            norm(na.cross(nb)).powi(2) / (1.0 + cosine)
        } else {
            norm(na - nb).powi(2) * 0.5
        };
        let dot = 0.5 * (mass * mass * (sb / sa) + bmass * bmass * (sa / sb)) + pa * pb * angle;
        if !dot.is_finite() {
            return Err("unrepresentable invariant product");
        }
        Ok(dot)
    }
    /// Center-of-momentum energy with known incoming masses (GeV).
    pub fn pair_mass(self, mass: f64, b: Self, bmass: f64) -> Result<f64, &'static str> {
        Ok((mass * mass + bmass * bmass + 2.0 * self.dot_on_shell(mass, b, bmass)?).sqrt())
    }
    /// On-shell boosts use light-cone components and avoid subtracting almost
    /// equal E and p when returning a fast, nearly collinear particle to rest.
    /// Reference must be timelike/on shell with the supplied reference mass.
    pub fn from_rest_on_shell(self, input_mass: f64, reference: Self, reference_mass: f64) -> Self {
        self.boost_shell(input_mass, reference, reference_mass, false)
    }
    pub fn to_rest_on_shell(self, input_mass: f64, reference: Self, reference_mass: f64) -> Self {
        self.boost_shell(input_mass, reference, reference_mass, true)
    }
    fn boost_shell(
        self,
        input_mass: f64,
        reference: Self,
        reference_mass: f64,
        inverse: bool,
    ) -> Self {
        let rp = norm(reference.momentum);
        if rp == 0.0 {
            return self;
        }
        let axis = reference.momentum / rp;
        let longitudinal = axis.dot(self.momentum);
        let perpendicular = self.momentum - axis * longitudinal;
        let transverse_mass2 = input_mass * input_mass + perpendicular.dot(perpendicular);
        let (mut plus, mut minus) = if longitudinal >= 0.0 {
            let plus = self.energy + longitudinal;
            (plus, transverse_mass2 / plus)
        } else {
            let minus = self.energy - longitudinal;
            (transverse_mass2 / minus, minus)
        };
        let factor = (reference.energy + rp) / reference_mass;
        if inverse {
            plus /= factor;
            minus *= factor;
        } else {
            plus *= factor;
            minus /= factor;
        }
        Self {
            energy: 0.5 * (plus + minus),
            momentum: perpendicular + axis * (0.5 * (plus - minus)),
        }
    }
    pub fn dot(self, b: Self) -> f64 {
        self.energy * b.energy - self.momentum.dot(b.momentum)
    }
    pub fn mass2(self) -> f64 {
        (self.energy - norm(self.momentum)) * (self.energy + norm(self.momentum))
    }
    /// Boost from the rest frame of a timelike reference, using its known mass.
    /// This avoids 1/sqrt(1-beta^2) cancellation at large gamma.
    pub fn from_rest(self, reference: Self, mass: f64) -> Self {
        let u = reference.momentum / mass;
        let gamma = reference.energy / mass;
        let up = u.dot(self.momentum);
        Self {
            energy: gamma * self.energy + up,
            momentum: self.momentum + u * (self.energy + up / (gamma + 1.0)),
        }
    }
    pub fn to_rest(self, reference: Self, mass: f64) -> Self {
        self.from_rest(
            Self {
                energy: reference.energy,
                momentum: -reference.momentum,
            },
            mass,
        )
    }
}
pub fn isotropic(rng: &mut Rng) -> Vec3 {
    let z = 2.0 * rng.uniform() - 1.0;
    let phi = 2.0 * std::f64::consts::PI * rng.uniform();
    let r = (1.0 - z * z).sqrt();
    Vec3([r * phi.cos(), r * phi.sin(), z])
}
/// Exact two-body kinematics; isotropic for an unpolarized parent.
/// This is a phase-space primitive, not a collision cross section.
pub fn two_body(
    parent: FourMomentum,
    mass: f64,
    a: f64,
    b: f64,
    rng: &mut Rng,
) -> Result<[FourMomentum; 2], &'static str> {
    if !mass.is_finite()
        || mass <= 0.0
        || !a.is_finite()
        || !b.is_finite()
        || a < 0.0
        || b < 0.0
        || mass < a + b
    {
        return Err("closed or invalid two-body phase space");
    }
    parent.validate_on_shell(mass)?;
    let p = 0.5
        * mass
        * ((1.0 - ((a + b) / mass).powi(2)) * (1.0 - ((a - b) / mass).powi(2)))
            .max(0.0)
            .sqrt();
    let v = isotropic(rng) * p;
    let ra = FourMomentum {
        energy: a.hypot(p),
        momentum: v,
    };
    let rb = FourMomentum {
        energy: b.hypot(p),
        momentum: -v,
    };
    Ok([
        ra.from_rest_on_shell(a, parent, mass),
        rb.from_rest_on_shell(b, parent, mass),
    ])
}
/// Unpolarized tree-level V-A charged-lepton decay with finite daughter mass,
/// massless neutrinos. Rejection from uniform Lorentz-invariant phase space.
/// Order: charged daughter, antineutrino of daughter flavor, parent neutrino.
pub fn leptonic(
    parent: FourMomentum,
    mass: f64,
    daughter: f64,
    rng: &mut Rng,
) -> Result<[FourMomentum; 3], &'static str> {
    if mass <= daughter || daughter < 0.0 || !mass.is_finite() || !daughter.is_finite() {
        return Err("invalid leptonic decay masses");
    }
    parent.validate_on_shell(mass)?;
    let rest = FourMomentum {
        energy: mass,
        momentum: Vec3::default(),
    };
    let smax = (mass - daughter).powi(2);
    loop {
        let s = smax * rng.uniform();
        // dPhi3/ds proportional to lambda^(1/2)(M²,m²,s)/M².
        let lam = ((1.0 - ((daughter + s.sqrt()) / mass).powi(2))
            * (1.0 - ((daughter - s.sqrt()) / mass).powi(2)))
        .max(0.0)
        .sqrt();
        if rng.uniform() > lam {
            continue;
        }
        let [charged, pair] = two_body(rest, mass, daughter, s.sqrt(), rng)?;
        let [anti, nu] = two_body(pair, s.sqrt(), 0.0, 0.0, rng)?;
        // |M|² proportional to (p_parent.p_anti)(p_charged.p_nu),
        // each factor <= M²/2; global bound M⁴/4 is safe.
        let weight = 4.0 * (mass * anti.energy) * charged.dot(nu) / mass.powi(4);
        if rng.uniform() <= weight {
            return Ok([
                charged.from_rest_on_shell(daughter, parent, mass),
                anti.from_rest_on_shell(0.0, parent, mass),
                nu.from_rest_on_shell(0.0, parent, mass),
            ]);
        }
    }
}
