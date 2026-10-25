//! The ladder: choose the cheapest rung whose claim meets the request, unfold
//! when it does not, and audit the choices against the reference rung.
//!
//! # Selection
//!
//! For an input the engine asks every rung for its a priori [`Assessment`]
//! (the cheap indicator), keeps those whose (calibrated) claim meets the
//! [`Request`] plus those that certify only after running, and evaluates them
//! in order of estimated cost until one is certified. If none is, it returns
//! the result with the smallest claimed error and `met = false`: the engine
//! never reports a tolerance it did not reach.
//!
//! # Audits and calibration
//!
//! Non-rigorous claims are estimates. With [`Audit::Every`]`(n)` the engine
//! re-runs the reference rung on the first and then every `n`-th accepted
//! result of each other rung, and measures the real error. An asymptotic or
//! validated rung that exceeds its claim gets a safety factor
//! (`1.5 ×` the observed excess) that multiplies its later claims, so the
//! ladder corrects itself where an indicator was too optimistic. A rigorous
//! rung that exceeds its claim is a defect: it is counted and never
//! calibrated away.

use std::fmt::Write;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Mutex;

use rayon::prelude::*;

use crate::claim::{Assessment, Claim, Request, Rigor};
use crate::rung::Rung;

/// Audit policy.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Audit {
    Off,
    /// Audit the first and then every `n`-th accepted result of each rung.
    Every(u64),
}

/// Safety margin applied on top of an observed excess (dimensionless; an
/// engineering choice, not a physical constant).
const CALIBRATION_MARGIN: f64 = 1.5;

/// Running statistics for one rung.
#[derive(Clone, Debug, Default)]
pub struct RungStats {
    pub name: String,
    /// Evaluations, including failed certificate attempts and audits.
    pub evaluations: u64,
    /// Estimated cost of those evaluations.
    pub cost: f64,
    /// Results returned to the caller from this rung.
    pub accepted: u64,
    /// Accepted results re-checked against the reference rung.
    pub audits: u64,
    /// Audits where the measured error exceeded the claim.
    pub violations: u64,
    /// Largest measured error / raw (uncalibrated) claim.
    pub worst_ratio: f64,
    /// Multiplier on non-rigorous claims (starts at 1).
    pub safety: f64,
    /// Audits where a rigorous claim was exceeded.
    pub rigorous_defects: u64,
}

/// Audit outcomes for one rung (updated only on audits, which are rare).
#[derive(Clone, Debug, Default)]
struct AuditTally {
    audits: u64,
    violations: u64,
    worst_ratio: f64,
    rigorous_defects: u64,
}

/// Lock-free live counters for one rung: the hot path only does relaxed
/// atomic adds, and batch paths merge thread-local tallies once.
#[derive(Debug)]
pub(crate) struct Live {
    name: String,
    evaluations: AtomicU64,
    cost: AtomicU64,
    accepted: AtomicU64,
    safety: AtomicU64,
    audit: Mutex<AuditTally>,
}

impl Live {
    fn new(name: String) -> Self {
        Self {
            name,
            evaluations: AtomicU64::new(0),
            cost: AtomicU64::new(0f64.to_bits()),
            accepted: AtomicU64::new(0),
            safety: AtomicU64::new(1f64.to_bits()),
            audit: Mutex::new(AuditTally::default()),
        }
    }
    fn add_cost(&self, c: f64) {
        let _ = self.cost.fetch_update(Ordering::Relaxed, Ordering::Relaxed, |b| Some((f64::from_bits(b) + c).to_bits()));
    }
    fn safety(&self) -> f64 {
        f64::from_bits(self.safety.load(Ordering::Relaxed))
    }
}

/// Per-rung counts accumulated locally by a batch and merged once.
#[derive(Clone, Debug)]
pub(crate) struct Tally {
    pub(crate) evaluations: Vec<u64>,
    pub(crate) cost: Vec<f64>,
    pub(crate) accepted: Vec<u64>,
}

impl Tally {
    pub(crate) fn new(n: usize) -> Self {
        Self {
            evaluations: vec![0; n],
            cost: vec![0.0; n],
            accepted: vec![0; n],
        }
    }
    pub(crate) fn merge(mut self, o: Tally) -> Self {
        for r in 0..self.cost.len() {
            self.evaluations[r] += o.evaluations[r];
            self.cost[r] += o.cost[r];
            self.accepted[r] += o.accepted[r];
        }
        self
    }
}

/// Deterministic audit sampling for batch item `i` (SplitMix64 hash), so
/// batches need no shared counter.
pub(crate) fn sampled(i: usize, every: u64) -> bool {
    let mut z = (i as u64).wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    (z ^ (z >> 31)) % every.max(1) == 0
}

/// The result of one solve.
#[derive(Clone, Debug)]
pub struct Solved<O> {
    pub value: O,
    /// Index of the rung that produced `value`.
    pub rung: usize,
    /// Its (calibrated) error claim.
    pub claim: Claim,
    /// Whether `claim` meets the request.
    pub met: bool,
    /// Estimated cost spent on this answer (excluding audits).
    pub cost: f64,
    /// Measured error against the reference, if this result was audited.
    pub audited: Option<f64>,
}

type Metric<O> = Box<dyn Fn(&O, &O) -> f64 + Send + Sync>;

/// A ladder of models for the same quantity.
pub struct Ladder<I, O> {
    name: String,
    pub(crate) rungs: Vec<Box<dyn Rung<I, O>>>,
    pub(crate) metric: Metric<O>,
    reference: Option<usize>,
    pub(crate) audit: Audit,
    pub(crate) live: Vec<Live>,
}

impl<I: Sync, O: Send> Ladder<I, O> {
    /// An empty ladder. `metric(a, b)` is the error of `a` against the exact
    /// `b`, in the units every rung's claims use.
    pub fn new(name: impl Into<String>, metric: impl Fn(&O, &O) -> f64 + Send + Sync + 'static) -> Self {
        Self {
            name: name.into(),
            rungs: Vec::new(),
            metric: Box::new(metric),
            reference: None,
            audit: Audit::Off,
            live: Vec::new(),
        }
    }

    /// Add a rung; returns its index.
    pub fn push(&mut self, rung: impl Rung<I, O> + 'static) -> usize {
        self.push_boxed(Box::new(rung))
    }

    pub fn push_boxed(&mut self, rung: Box<dyn Rung<I, O>>) -> usize {
        self.live.push(Live::new(rung.name().to_string()));
        self.rungs.push(rung);
        self.rungs.len() - 1
    }

    /// Builder form of [`Ladder::push`].
    pub fn with(mut self, rung: impl Rung<I, O> + 'static) -> Self {
        self.push(rung);
        self
    }

    /// Set the reference rung (the finest model) used by audits.
    pub fn reference(mut self, index: usize) -> Self {
        assert!(index < self.rungs.len(), "reference rung {index} does not exist");
        self.reference = Some(index);
        self
    }

    pub fn audit(mut self, policy: Audit) -> Self {
        self.audit = policy;
        self
    }

    pub fn set_audit(&mut self, policy: Audit) {
        self.audit = policy;
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn len(&self) -> usize {
        self.rungs.len()
    }

    pub fn is_empty(&self) -> bool {
        self.rungs.is_empty()
    }

    pub fn rung_name(&self, index: usize) -> &str {
        self.rungs[index].name()
    }

    pub fn rung(&self, index: usize) -> &dyn Rung<I, O> {
        self.rungs[index].as_ref()
    }

    /// The error metric.
    pub fn error(&self, a: &O, exact: &O) -> f64 {
        (self.metric)(a, exact)
    }

    /// A snapshot of the per-rung statistics.
    pub fn stats(&self) -> Vec<RungStats> {
        self.live
            .iter()
            .map(|l| {
                let a = l.audit.lock().unwrap().clone();
                RungStats {
                    name: l.name.clone(),
                    evaluations: l.evaluations.load(Ordering::Relaxed),
                    cost: f64::from_bits(l.cost.load(Ordering::Relaxed)),
                    accepted: l.accepted.load(Ordering::Relaxed),
                    audits: a.audits,
                    violations: a.violations,
                    worst_ratio: a.worst_ratio,
                    safety: l.safety(),
                    rigorous_defects: a.rigorous_defects,
                }
            })
            .collect()
    }

    /// Clear counters but keep the learned safety factors.
    pub fn reset_counters(&self) {
        for l in &self.live {
            l.evaluations.store(0, Ordering::Relaxed);
            l.cost.store(0f64.to_bits(), Ordering::Relaxed);
            l.accepted.store(0, Ordering::Relaxed);
            *l.audit.lock().unwrap() = AuditTally::default();
        }
    }

    /// Apply the learned safety factor to a non-rigorous claim.
    pub(crate) fn calibrate(&self, rung: usize, c: Claim) -> Claim {
        if c.rigor == Rigor::Rigorous {
            return c;
        }
        Claim::new(c.error * self.live[rung].safety(), c.rigor)
    }

    /// The calibrated a priori assessment of one rung.
    pub fn assess(&self, rung: usize, input: &I) -> Assessment {
        match self.rungs[rung].assess(input) {
            Assessment::Claim(c) => Assessment::Claim(self.calibrate(rung, c)),
            a => a,
        }
    }

    pub(crate) fn record_eval(&self, rung: usize, cost: f64) {
        self.live[rung].evaluations.fetch_add(1, Ordering::Relaxed);
        self.live[rung].add_cost(cost);
    }

    pub(crate) fn merge_tally(&self, t: &Tally) {
        for (r, l) in self.live.iter().enumerate() {
            l.evaluations.fetch_add(t.evaluations[r], Ordering::Relaxed);
            l.accepted.fetch_add(t.accepted[r], Ordering::Relaxed);
            l.add_cost(t.cost[r]);
        }
    }

    /// Whether the audit policy samples this result of rung `r`.
    pub(crate) fn wants_audit(&self, r: usize, nth: u64) -> bool {
        match (self.audit, self.reference) {
            (Audit::Every(n), Some(reference)) if r != reference => nth % n.max(1) == 0,
            _ => false,
        }
    }

    /// Solve one input to the request.
    pub fn solve(&self, input: &I, req: Request) -> Solved<O> {
        // (cost, rung, a priori claim if any)
        let mut options: Vec<(f64, usize, Option<Claim>)> = Vec::new();
        // the in-domain rung with the smallest a priori claim, if none meets
        let mut fallback: Option<(Claim, usize)> = None;
        for r in 0..self.rungs.len() {
            match self.assess(r, input) {
                Assessment::Claim(c) => {
                    if c.meets(&req) {
                        options.push((self.rungs[r].cost(input), r, Some(c)));
                    } else if fallback.is_none_or(|(f, _)| better(c, f)) {
                        fallback = Some((c, r));
                    }
                }
                Assessment::Deferred => options.push((self.rungs[r].cost(input), r, None)),
                Assessment::OutOfDomain => {}
            }
        }
        options.sort_by(|a, b| a.0.total_cmp(&b.0));

        let mut spent = 0.0;
        let mut best: Option<(O, usize, Claim)> = None;
        for (cost, r, prior) in options {
            let out = self.rungs[r].eval(input);
            spent += cost;
            self.record_eval(r, cost);
            let post = self.rungs[r].certify(input, &out).map(|c| self.calibrate(r, c));
            let Some(claim) = post.or(prior) else { continue };
            if claim.meets(&req) {
                return self.accept(input, r, out, claim, spent, true);
            }
            if best.as_ref().is_none_or(|(_, _, b)| better(claim, *b)) {
                best = Some((out, r, claim));
            }
        }

        // Nothing met the request: return the most accurate answer available.
        match (best, fallback) {
            (Some((out, r, claim)), Some((fc, _))) if !better(fc, claim) => {
                self.accept(input, r, out, claim, spent, false)
            }
            (_, Some((fc, fr))) => {
                let cost = self.rungs[fr].cost(input);
                let out = self.rungs[fr].eval(input);
                self.record_eval(fr, cost);
                self.accept(input, fr, out, fc, spent + cost, false)
            }
            (Some((out, r, claim)), None) => self.accept(input, r, out, claim, spent, false),
            (None, None) => panic!("ladder `{}`: no rung applies to this input", self.name),
        }
    }

    /// Solve many independent inputs in parallel.
    pub fn solve_many(&self, inputs: &[I], req: Request) -> Vec<Solved<O>> {
        inputs.par_iter().map(|i| self.solve(i, req)).collect()
    }

    /// Record an accepted result and audit it if the policy says so.
    pub(crate) fn accept(&self, input: &I, r: usize, out: O, claim: Claim, cost: f64, met: bool) -> Solved<O> {
        let nth = self.live[r].accepted.fetch_add(1, Ordering::Relaxed);
        let audit = self.wants_audit(r, nth);
        self.finish(input, r, out, claim, cost, met, audit)
    }

    /// Build the result, auditing it when `audit` is set (no counting).
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn finish(&self, input: &I, r: usize, out: O, claim: Claim, cost: f64, met: bool, audit: bool) -> Solved<O> {
        let audited = if audit { Some(self.audit_one(input, r, &out, claim)) } else { None };
        Solved {
            value: out,
            rung: r,
            claim,
            met,
            cost,
            audited,
        }
    }

    /// Measure the real error of an accepted result against the reference
    /// rung, and calibrate. Returns the measured error.
    fn audit_one(&self, input: &I, r: usize, out: &O, claim: Claim) -> f64 {
        let reference = self.reference.expect("audit without a reference rung");
        let rr = &self.rungs[reference];
        let cost = rr.cost(input);
        let exact = rr.eval(input);
        self.record_eval(reference, cost);
        let ref_error = rr
            .certify(input, &exact)
            .or(match rr.assess(input) {
                Assessment::Claim(c) => Some(c),
                _ => None,
            })
            .map_or(0.0, |c| c.error);
        let measured = (self.metric)(out, &exact);
        // Only the part of the difference the reference cannot explain counts.
        let excess = (measured - ref_error).max(0.0);
        let live = &self.live[r];
        let safety = live.safety();
        let raw = claim.error / if claim.rigor == Rigor::Rigorous { 1.0 } else { safety };
        let ratio = if raw > 0.0 {
            excess / raw
        } else if excess > 0.0 {
            f64::INFINITY
        } else {
            0.0
        };
        let mut a = live.audit.lock().unwrap();
        a.audits += 1;
        a.worst_ratio = a.worst_ratio.max(ratio);
        if excess > claim.error {
            a.violations += 1;
            if claim.rigor == Rigor::Rigorous {
                a.rigorous_defects += 1;
            } else if ratio.is_finite() {
                let new = safety.max(CALIBRATION_MARGIN * ratio);
                live.safety.store(new.to_bits(), Ordering::Relaxed);
            }
        }
        measured
    }

    /// A plain-text table of the per-rung statistics.
    pub fn report(&self) -> String {
        let st = self.stats();
        let total: f64 = st.iter().map(|s| s.cost).sum();
        let mut out = String::new();
        let _ = writeln!(out, "ladder `{}`", self.name);
        let _ = writeln!(
            out,
            "  {:<26} {:>9} {:>9} {:>7} {:>6} {:>5} {:>10} {:>7}",
            "rung", "accepted", "evals", "cost%", "audits", "viol", "worst", "safety"
        );
        for s in &st {
            let _ = writeln!(
                out,
                "  {:<26} {:>9} {:>9} {:>6.1}% {:>6} {:>5} {:>10.3e} {:>7.2}",
                s.name,
                s.accepted,
                s.evaluations,
                if total > 0.0 { 100.0 * s.cost / total } else { 0.0 },
                s.audits,
                s.violations,
                s.worst_ratio,
                s.safety
            );
        }
        out
    }
}

/// `a` is a better claim than `b`: smaller error, then higher rigor.
fn better(a: Claim, b: Claim) -> bool {
    match a.error.total_cmp(&b.error) {
        std::cmp::Ordering::Less => true,
        std::cmp::Ordering::Equal => a.rigor > b.rigor,
        std::cmp::Ordering::Greater => false,
    }
}
