//! Values with bounds: a nominal number, an interval that contains the true
//! value, its dimension, and how far that interval can be trusted.
//!
//! The interval is what crosses team boundaries. A thermal engineer does not
//! receive "2.1 W" from electrical but "2.1 W in [1.8, 2.5] W, rigorous given
//! the stated tolerances", so a requirement downstream can be called met, failed
//! or *undecided* honestly. Interval arithmetic is rigorous up to floating
//! point round-off (no directed rounding), which is what [`Rigor::Rigorous`]
//! means throughout the workspace.

use emergence::Rigor;

use crate::units::Dim;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Val {
    pub nom: f64,
    pub lo: f64,
    pub hi: f64,
    pub dim: Dim,
    /// How far `[lo, hi]` can be trusted to contain the true value.
    pub rigor: Rigor,
}

fn same(a: &Val, b: &Val, op: &str) -> Result<Dim, String> {
    if a.dim != b.dim {
        return Err(format!("cannot {op} {} and {}", a.dim, b.dim));
    }
    Ok(a.dim)
}

impl Val {
    /// An exact number.
    pub fn exact(v: f64, dim: Dim) -> Val {
        Val { nom: v, lo: v, hi: v, dim, rigor: Rigor::Rigorous }
    }

    pub fn new(nom: f64, lo: f64, hi: f64, dim: Dim, rigor: Rigor) -> Result<Val, String> {
        if !(lo <= nom && nom <= hi) || !nom.is_finite() || !lo.is_finite() || !hi.is_finite() {
            return Err(format!("bad interval: {nom} not in [{lo}, {hi}]"));
        }
        Ok(Val { nom, lo, hi, dim, rigor })
    }

    pub fn width(&self) -> f64 {
        self.hi - self.lo
    }

    /// Half-width relative to the nominal value.
    pub fn rel(&self) -> f64 {
        if self.nom == 0.0 {
            if self.width() == 0.0 { 0.0 } else { f64::INFINITY }
        } else {
            0.5 * self.width() / self.nom.abs()
        }
    }

    pub fn is_exact(&self) -> bool {
        self.lo == self.hi
    }

    pub fn add(&self, o: &Val) -> Result<Val, String> {
        let dim = same(self, o, "add")?;
        Ok(Val { nom: self.nom + o.nom, lo: self.lo + o.lo, hi: self.hi + o.hi, dim, rigor: self.rigor.min(o.rigor) })
    }

    pub fn sub(&self, o: &Val) -> Result<Val, String> {
        let dim = same(self, o, "subtract")?;
        Ok(Val { nom: self.nom - o.nom, lo: self.lo - o.hi, hi: self.hi - o.lo, dim, rigor: self.rigor.min(o.rigor) })
    }

    pub fn neg(&self) -> Val {
        Val { nom: -self.nom, lo: -self.hi, hi: -self.lo, ..*self }
    }

    pub fn mul(&self, o: &Val) -> Val {
        let p = [self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi];
        Val {
            nom: self.nom * o.nom,
            lo: p.iter().copied().fold(f64::INFINITY, f64::min),
            hi: p.iter().copied().fold(f64::NEG_INFINITY, f64::max),
            dim: self.dim.mul(o.dim),
            rigor: self.rigor.min(o.rigor),
        }
    }

    pub fn div(&self, o: &Val) -> Result<Val, String> {
        if o.lo <= 0.0 && o.hi >= 0.0 {
            return Err(format!("division by an interval containing zero [{}, {}]", o.lo, o.hi));
        }
        let inv = Val { nom: 1.0 / o.nom, lo: 1.0 / o.hi, hi: 1.0 / o.lo, dim: Dim::NONE.div(o.dim), rigor: o.rigor };
        Ok(self.mul(&inv))
    }

    pub fn powi(&self, n: i32) -> Result<Val, String> {
        if n == 0 {
            return Ok(Val::exact(1.0, Dim::NONE));
        }
        if n < 0 {
            return Val::exact(1.0, Dim::NONE).div(&self.powi(-n)?);
        }
        let (a, b) = (self.lo.powi(n), self.hi.powi(n));
        let (lo, hi) = if n % 2 == 0 && self.lo < 0.0 && self.hi > 0.0 {
            (0.0, a.max(b))
        } else if n % 2 == 0 && self.hi <= 0.0 {
            (b, a)
        } else {
            (a, b)
        };
        Ok(Val { nom: self.nom.powi(n), lo, hi, dim: self.dim.powi(n), rigor: self.rigor })
    }

    /// Apply a monotone increasing function to a dimensionless value.
    fn increasing(&self, name: &str, f: fn(f64) -> f64, domain_lo: f64) -> Result<Val, String> {
        if !self.dim.is_none() {
            return Err(format!("{name}() needs a dimensionless argument, got {}", self.dim));
        }
        if self.lo < domain_lo {
            return Err(format!("{name}() of an interval reaching {}", self.lo));
        }
        Ok(Val { nom: f(self.nom), lo: f(self.lo), hi: f(self.hi), dim: Dim::NONE, rigor: self.rigor })
    }

    pub fn sqrt(&self) -> Result<Val, String> {
        let dim = self.dim.sqrt().ok_or_else(|| format!("sqrt of {} has no dimension", self.dim))?;
        if self.lo < 0.0 {
            return Err(format!("sqrt of an interval reaching {}", self.lo));
        }
        Ok(Val { nom: self.nom.sqrt(), lo: self.lo.sqrt(), hi: self.hi.sqrt(), dim, rigor: self.rigor })
    }
    pub fn exp(&self) -> Result<Val, String> {
        self.increasing("exp", f64::exp, f64::NEG_INFINITY)
    }
    pub fn ln(&self) -> Result<Val, String> {
        if self.lo <= 0.0 {
            return Err(format!("ln of an interval reaching {}", self.lo));
        }
        self.increasing("ln", f64::ln, 0.0)
    }
    pub fn tanh(&self) -> Result<Val, String> {
        self.increasing("tanh", f64::tanh, f64::NEG_INFINITY)
    }
    pub fn abs(&self) -> Val {
        let (lo, hi) = if self.lo >= 0.0 {
            (self.lo, self.hi)
        } else if self.hi <= 0.0 {
            (-self.hi, -self.lo)
        } else {
            (0.0, (-self.lo).max(self.hi))
        };
        Val { nom: self.nom.abs(), lo, hi, ..*self }
    }
    pub fn min(&self, o: &Val) -> Result<Val, String> {
        let dim = same(self, o, "compare")?;
        Ok(Val { nom: self.nom.min(o.nom), lo: self.lo.min(o.lo), hi: self.hi.min(o.hi), dim, rigor: self.rigor.min(o.rigor) })
    }
    pub fn max(&self, o: &Val) -> Result<Val, String> {
        let dim = same(self, o, "compare")?;
        Ok(Val { nom: self.nom.max(o.nom), lo: self.lo.max(o.lo), hi: self.hi.max(o.hi), dim, rigor: self.rigor.min(o.rigor) })
    }

    /// Widen by an absolute model error `err` (same units) of rigor `rigor`.
    pub fn widened(&self, err: f64, rigor: Rigor) -> Val {
        Val { lo: self.lo - err, hi: self.hi + err, rigor: self.rigor.min(rigor), ..*self }
    }

    /// The same value pinned at its nominal (used to measure how much one
    /// input contributes to a downstream interval).
    pub fn pinned(&self) -> Val {
        Val { lo: self.nom, hi: self.nom, ..*self }
    }
}
