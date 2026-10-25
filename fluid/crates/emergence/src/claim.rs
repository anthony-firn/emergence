//! Error claims and requests.
//!
//! Every rung states how far its answer may be from the exact one, in the
//! units of the ladder's error metric, and how much that statement can be
//! trusted. The engine never compares a number without its rigor.

/// How much an error claim can be trusted, weakest first.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Rigor {
    /// A leading-order estimate: `C · indicator^p`, or the difference to an
    /// embedded higher-order model. It can be exceeded outside the asymptotic
    /// regime (near a resonance, at high contrast), which is what audits catch.
    Asymptotic,
    /// An estimate checked against the finer model on held-out points of the
    /// same domain (closure tables), or an empirically established bound with
    /// a published validation (e.g. Wiscombe's Mie truncation).
    Validated,
    /// A proven upper bound, up to floating-point round-off.
    Rigorous,
}

impl Rigor {
    pub fn name(self) -> &'static str {
        match self {
            Rigor::Asymptotic => "asymptotic",
            Rigor::Validated => "validated",
            Rigor::Rigorous => "rigorous",
        }
    }
}

/// A rung's statement about its own error.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Claim {
    /// Error bound or estimate, in the ladder's metric (non-negative).
    pub error: f64,
    pub rigor: Rigor,
}

impl Claim {
    pub fn new(error: f64, rigor: Rigor) -> Self {
        Self { error, rigor }
    }
    pub fn rigorous(error: f64) -> Self {
        Self::new(error, Rigor::Rigorous)
    }
    pub fn validated(error: f64) -> Self {
        Self::new(error, Rigor::Validated)
    }
    pub fn asymptotic(error: f64) -> Self {
        Self::new(error, Rigor::Asymptotic)
    }

    /// The claim on a sum: errors add (triangle inequality) and the rigor is
    /// the weakest of the two.
    pub fn plus(self, other: Claim) -> Claim {
        Claim::new(self.error + other.error, self.rigor.min(other.rigor))
    }

    /// The claim after a linear map with gain `s` (`|s|` times the error).
    pub fn scaled(self, s: f64) -> Claim {
        Claim::new(self.error * s.abs(), self.rigor)
    }

    /// True when this claim satisfies `req`.
    pub fn meets(&self, req: &Request) -> bool {
        self.error.is_finite() && self.error <= req.tol && self.rigor >= req.min_rigor
    }
}

/// What the caller needs: an error tolerance in the ladder's metric and the
/// weakest rigor it will accept.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Request {
    pub tol: f64,
    pub min_rigor: Rigor,
}

impl Request {
    /// Tolerance `tol`, any rigor.
    pub fn new(tol: f64) -> Self {
        Self {
            tol,
            min_rigor: Rigor::Asymptotic,
        }
    }
    pub fn with_rigor(tol: f64, min_rigor: Rigor) -> Self {
        Self { tol, min_rigor }
    }
}

/// A rung's a priori answer to "does your model apply here, and how well?".
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Assessment {
    /// The cheap indicator gives an error claim before evaluation.
    Claim(Claim),
    /// The model applies, but its error is known only after evaluation
    /// (a posteriori certificate, see [`crate::Rung::certify`]).
    Deferred,
    /// The indicator says the model does not apply (outside its domain).
    OutOfDomain,
}
