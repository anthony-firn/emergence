//! Goal-oriented error budgets.
//!
//! Many models are evaluated to compute one number: a quadrature over a size
//! distribution, a sum over the cells of a field, or anything downstream that
//! depends linearly (or after linearisation) on the local results,
//!
//! ```text
//! J = Σ_i w_i f(x_i),     |J_computed − J| ≤ Σ_i |w_i| e_i ,
//! ```
//!
//! where `e_i` is the claimed error of item `i` and `w_i` is its weight
//! (quadrature weight, cell volume, or the adjoint sensitivity `∂J/∂f_i`).
//! A uniform per-item tolerance wastes work: items with small weight can use
//! a crude rung. [`Ladder::plan`] picks one rung per item so that
//! `Σ |w_i| e_i ≤ tol` at the least total cost. This is a multiple-choice
//! knapsack; the greedy walk along each item's lower convex hull of
//! `(cost, |w| e)` solves its linear relaxation exactly and is within one
//! upgrade of the integer optimum.
//!
//! Items with no a priori claim (deferred certificates only) are solved one by
//! one with an equal share of the budget.

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use rayon::prelude::*;

use crate::claim::{Assessment, Claim, Request, Rigor};
use crate::ladder::{sampled, Audit, Ladder, Solved, Tally};

/// The rung chosen for every item, and the predicted bound.
#[derive(Clone, Debug)]
pub struct Plan {
    /// Chosen rung per item (`None`: solved individually, see module docs).
    pub choice: Vec<Option<usize>>,
    /// The a priori (calibrated) claim of each planned item.
    pub claims: Vec<Option<Claim>>,
    /// Predicted cost of the planned items.
    pub cost: f64,
    /// Predicted `Σ |w_i| e_i` over the planned items.
    pub bound: f64,
    /// Budget reserved for the individually solved items.
    pub reserved: f64,
}

/// The evaluated functional inputs.
#[derive(Clone, Debug)]
pub struct Functional<O> {
    pub items: Vec<Solved<O>>,
    pub plan: Plan,
    /// The claim on `J`: `Σ |w_i| e_i`, with the weakest rigor used.
    pub claim: Claim,
    pub met: bool,
    /// Total estimated cost (excluding audits).
    pub cost: f64,
}

impl<O> Functional<O> {
    /// `Σ w_i g(value_i)` for a scalar observable `g`.
    pub fn sum(&self, weights: &[f64], g: impl Fn(&O) -> f64) -> f64 {
        self.items.iter().zip(weights).map(|(s, w)| w * g(&s.value)).sum()
    }

    /// Number of items solved by each rung.
    pub fn histogram(&self, n_rungs: usize) -> Vec<usize> {
        let mut h = vec![0; n_rungs];
        for s in &self.items {
            h[s.rung] += 1;
        }
        h
    }
}

/// One option of one item: `(cost, weighted error, rung, claim)`.
#[derive(Clone, Copy, Debug)]
struct Option_ {
    cost: f64,
    err: f64,
    rung: usize,
    claim: Claim,
}

/// Heap entry: the marginal error reduction per unit cost of upgrading
/// `item` to hull position `next`.
struct Upgrade {
    ratio: f64,
    item: usize,
    next: usize,
}
impl PartialEq for Upgrade {
    fn eq(&self, o: &Self) -> bool {
        self.ratio == o.ratio
    }
}
impl Eq for Upgrade {}
impl PartialOrd for Upgrade {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}
impl Ord for Upgrade {
    fn cmp(&self, o: &Self) -> Ordering {
        self.ratio.total_cmp(&o.ratio).then(o.item.cmp(&self.item))
    }
}

/// Lower convex hull of options sorted by cost: strictly decreasing error,
/// and decreasing marginal gain `Δerr/Δcost` along the hull.
fn convex_hull(mut opts: Vec<Option_>) -> Vec<Option_> {
    opts.sort_by(|a, b| a.cost.total_cmp(&b.cost).then(a.err.total_cmp(&b.err)));
    let mut hull: Vec<Option_> = Vec::new();
    for o in opts {
        if hull.last().is_some_and(|l| o.err >= l.err) {
            continue;
        }
        while hull.len() >= 2 {
            let (p, q) = (hull[hull.len() - 2], hull[hull.len() - 1]);
            // gain per cost of p→q versus q→o; q is redundant if not steeper.
            let g1 = (p.err - q.err) * (o.cost - q.cost);
            let g2 = (q.err - o.err) * (q.cost - p.cost);
            if g1 <= g2 {
                hull.pop();
            } else {
                break;
            }
        }
        hull.push(o);
    }
    hull
}

fn ratio(a: &Option_, b: &Option_) -> f64 {
    let dc = b.cost - a.cost;
    if dc > 0.0 {
        (a.err - b.err) / dc
    } else {
        f64::INFINITY
    }
}

impl<I: Sync, O: Send + Sync> Ladder<I, O> {
    /// The options of one item that satisfy the rigor requirement.
    fn options(&self, input: &I, w: f64, req: &Request, mut each: impl FnMut(Option_)) {
        for r in 0..self.len() {
            if let Assessment::Claim(c) = self.assess(r, input) {
                if c.rigor >= req.min_rigor && c.error.is_finite() {
                    each(Option_ {
                        cost: self.rungs[r].cost(input),
                        err: w.abs() * c.error,
                        rung: r,
                        claim: c,
                    });
                }
            }
        }
    }

    /// Choose a rung per item so that `Σ |w_i| e_i ≤ req.tol` at least cost.
    ///
    /// Fast path: the cheapest option of every item is found without
    /// allocating; the convex hulls and the upgrade heap are built only when
    /// those choices exceed the budget.
    pub fn plan(&self, inputs: &[I], weights: &[f64], req: Request) -> Plan {
        assert_eq!(inputs.len(), weights.len());
        let cheapest: Vec<Option<Option_>> = inputs
            .par_iter()
            .zip(weights.par_iter())
            .map(|(input, &w)| {
                let mut best: Option<Option_> = None;
                self.options(input, w, &req, |o| {
                    let better = best.is_none_or(|b| o.cost < b.cost || (o.cost == b.cost && o.err < b.err));
                    if better {
                        best = Some(o);
                    }
                });
                best
            })
            .collect();

        // (Σ|w|, Σ|w| of unplannable items, Σ|w| e of the cheapest choices)
        let (w_total, w_deferred, mut bound) = cheapest
            .par_iter()
            .zip(weights.par_iter())
            .map(|(h, w)| match h {
                Some(o) => (w.abs(), 0.0, o.err),
                None => (w.abs(), w.abs(), 0.0),
            })
            .reduce(|| (0.0, 0.0, 0.0), |a, b| (a.0 + b.0, a.1 + b.1, a.2 + b.2));
        let reserved = if w_total > 0.0 { req.tol * w_deferred / w_total } else { 0.0 };
        let budget = req.tol - reserved;

        let mut chosen = cheapest;
        if bound > budget {
            let hulls: Vec<Vec<Option_>> = inputs
                .par_iter()
                .zip(weights.par_iter())
                .map(|(input, &w)| {
                    let mut opts = Vec::new();
                    self.options(input, w, &req, |o| opts.push(o));
                    convex_hull(opts)
                })
                .collect();
            let mut pos = vec![0usize; inputs.len()];
            let mut heap = BinaryHeap::new();
            for (i, h) in hulls.iter().enumerate() {
                if h.len() > 1 {
                    heap.push(Upgrade {
                        ratio: ratio(&h[0], &h[1]),
                        item: i,
                        next: 1,
                    });
                }
            }
            while bound > budget {
                let Some(u) = heap.pop() else { break };
                let h = &hulls[u.item];
                bound += h[u.next].err - h[pos[u.item]].err;
                pos[u.item] = u.next;
                if u.next + 1 < h.len() {
                    heap.push(Upgrade {
                        ratio: ratio(&h[u.next], &h[u.next + 1]),
                        item: u.item,
                        next: u.next + 1,
                    });
                }
            }
            chosen = hulls.iter().zip(&pos).map(|(h, &p)| h.get(p).copied()).collect();
        }

        let (choice, claims): (Vec<Option<usize>>, Vec<Option<Claim>>) =
            chosen.par_iter().map(|o| (o.map(|o| o.rung), o.map(|o| o.claim))).unzip();
        let cost = chosen.par_iter().map(|o| o.map_or(0.0, |o| o.cost)).sum();
        Plan {
            choice,
            claims,
            cost,
            bound,
            reserved,
        }
    }

    /// Plan, then evaluate every item in parallel. The returned claim is the
    /// bound on `Σ w_i f(x_i)` (in the metric's units per unit weight).
    pub fn solve_functional(&self, inputs: &[I], weights: &[f64], req: Request) -> Functional<O> {
        let plan = self.plan(inputs, weights, req);
        self.execute(plan, inputs, weights, req, None)
    }

    /// Re-plan with new weights (typically after a pilot solve fixed the
    /// linearisation point of a nonlinear observable), reusing every item of
    /// `previous` whose chosen rung is unchanged; only changed items are
    /// evaluated and counted in `cost`.
    pub fn refine_functional(&self, previous: Functional<O>, inputs: &[I], weights: &[f64], req: Request) -> Functional<O> {
        assert_eq!(previous.items.len(), inputs.len());
        let plan = self.plan(inputs, weights, req);
        self.execute(plan, inputs, weights, req, Some(previous.items))
    }

    fn execute(&self, plan: Plan, inputs: &[I], weights: &[f64], req: Request, reuse: Option<Vec<Solved<O>>>) -> Functional<O> {
        let w_deferred: f64 = plan
            .choice
            .iter()
            .zip(weights)
            .filter(|(c, _)| c.is_none())
            .map(|(_, w)| w.abs())
            .sum();
        let every = match self.audit {
            Audit::Every(n) => n,
            Audit::Off => 0,
        };
        // Planned items touch no shared state: counts are tallied afterwards
        // and audits are sampled by a hash of the item index.
        let fresh = |i: usize| -> (Solved<O>, bool) {
            let input = &inputs[i];
            match (plan.choice[i], plan.claims[i]) {
                (Some(r), Some(prior)) => {
                    let cost = self.rungs[r].cost(input);
                    let out = self.rungs[r].eval(input);
                    let claim = self.rungs[r].certify(input, &out).map(|c| self.calibrate(r, c)).unwrap_or(prior);
                    let audit = every > 0 && self.wants_audit(r, 0) && sampled(i, every);
                    (self.finish(input, r, out, claim, cost, true, audit), true)
                }
                _ => {
                    // Equal share of the reserved budget per unit weight.
                    let tol = if w_deferred > 0.0 { plan.reserved / w_deferred } else { req.tol };
                    (self.solve(input, Request::with_rigor(tol, req.min_rigor)), false)
                }
            }
        };
        let (items, planned): (Vec<Solved<O>>, Vec<bool>) = match reuse {
            None => (0..inputs.len()).into_par_iter().map(fresh).unzip(),
            Some(prev) => prev
                .into_par_iter()
                .enumerate()
                .map(|(i, p)| {
                    if plan.choice[i] == Some(p.rung) {
                        (Solved { cost: 0.0, ..p }, false)
                    } else {
                        fresh(i)
                    }
                })
                .unzip(),
        };
        let n_rungs = self.len();
        let tally = items
            .par_iter()
            .zip(planned.par_iter())
            .filter(|(_, &p)| p)
            .fold(
                || Tally::new(n_rungs),
                |mut t, (s, _)| {
                    t.evaluations[s.rung] += 1;
                    t.accepted[s.rung] += 1;
                    t.cost[s.rung] += s.cost;
                    t
                },
            )
            .reduce(|| Tally::new(n_rungs), Tally::merge);
        self.merge_tally(&tally);
        let zero = (Claim::new(0.0, Rigor::Rigorous), 0.0);
        let (claim, cost) = items
            .par_iter()
            .zip(weights.par_iter())
            .map(|(s, w)| (if *w != 0.0 { s.claim.scaled(*w) } else { zero.0 }, s.cost))
            .reduce(|| zero, |a, b| (a.0.plus(b.0), a.1 + b.1));
        let met = claim.meets(&req);
        Functional {
            items,
            plan,
            claim,
            met,
            cost,
        }
    }
}
