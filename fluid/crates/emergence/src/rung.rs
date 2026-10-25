//! A rung: one model of a ladder, with its cost, its validity indicator and
//! (optionally) an a posteriori certificate.

use crate::claim::{Assessment, Claim};

/// One model on a ladder from input `I` to output `O`.
///
/// `assess` must be cheap compared with `eval`: it is the indicator, and the
/// engine calls it for every rung before choosing one.
pub trait Rung<I, O>: Send + Sync {
    /// Short name for reports.
    fn name(&self) -> &str;

    /// Estimated cost of [`Rung::eval`] at `input`, in units shared by all the
    /// rungs of a ladder (the engine only compares and sums them).
    fn cost(&self, input: &I) -> f64;

    /// A priori validity and error claim (the indicator).
    fn assess(&self, input: &I) -> Assessment;

    /// Run the model.
    fn eval(&self, input: &I) -> O;

    /// A posteriori certificate computed from the output. When it returns
    /// `Some`, it replaces the a priori claim. The default has none.
    fn certify(&self, _input: &I, _out: &O) -> Option<Claim> {
        None
    }
}

type CostFn<I> = Box<dyn Fn(&I) -> f64 + Send + Sync>;
type AssessFn<I> = Box<dyn Fn(&I) -> Assessment + Send + Sync>;
type EvalFn<I, O> = Box<dyn Fn(&I) -> O + Send + Sync>;
type CertifyFn<I, O> = Box<dyn Fn(&I, &O) -> Option<Claim> + Send + Sync>;

/// A rung assembled from closures: the quickest way to put an existing
/// function and its indicator on a ladder.
pub struct FnRung<I, O> {
    name: String,
    cost: CostFn<I>,
    assess: AssessFn<I>,
    eval: EvalFn<I, O>,
    certify: Option<CertifyFn<I, O>>,
}

impl<I, O> FnRung<I, O> {
    pub fn new(
        name: impl Into<String>,
        cost: impl Fn(&I) -> f64 + Send + Sync + 'static,
        assess: impl Fn(&I) -> Assessment + Send + Sync + 'static,
        eval: impl Fn(&I) -> O + Send + Sync + 'static,
    ) -> Self {
        Self {
            name: name.into(),
            cost: Box::new(cost),
            assess: Box::new(assess),
            eval: Box::new(eval),
            certify: None,
        }
    }

    /// Attach an a posteriori certificate.
    pub fn with_certificate(
        mut self,
        certify: impl Fn(&I, &O) -> Option<Claim> + Send + Sync + 'static,
    ) -> Self {
        self.certify = Some(Box::new(certify));
        self
    }
}

impl<I, O> Rung<I, O> for FnRung<I, O> {
    fn name(&self) -> &str {
        &self.name
    }
    fn cost(&self, input: &I) -> f64 {
        (self.cost)(input)
    }
    fn assess(&self, input: &I) -> Assessment {
        (self.assess)(input)
    }
    fn eval(&self, input: &I) -> O {
        (self.eval)(input)
    }
    fn certify(&self, input: &I, out: &O) -> Option<Claim> {
        self.certify.as_ref().and_then(|c| c(input, out))
    }
}
