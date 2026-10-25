# The emergence engine: fold/unfold as shared, tested code

Implementation: `fluid/crates/emergence` (engine), with ladders in
`photonics::mie_ladder`, `heat::bte_ladder`, `spectral_core::ladder`, and the
cross-crate test `testlab/tests/emergence_links.rs`. The organising idea is
[`emergence_ladder.md`](emergence_ladder.md); until this crate it was a pattern
re-implemented by hand in every crate. Now it is one engine that every crate
plugs rungs into.

All timings: i7-13620H (16 threads), `cargo test` profile (opt-level 3).

## What the engine does

**Typed claims.** Every rung states its error in the ladder's metric *and*
how far that statement can be trusted:

| rigor | meaning | examples |
|---|---|---|
| `Asymptotic` | leading-order estimate; can fail outside its regime | Rayleigh `abs(R2 − R1)`, Milne-jump `0.05 E₂(1/Kn)` |
| `Validated` | checked against the finer model on held-out points, or an empirically established bound | closure tables; Wiscombe truncation (measured ≤ 7.3e-9) |
| `Rigorous` | proven bound (up to round-off) | weak–strong `L²` certificate of the spectral solver; Lagrange remainder |

A request is `(tol, min_rigor)`. The engine never reports a tolerance it did
not reach: unmet requests return the best available answer with `met = false`.

**Selection.** For each input the engine asks every rung for its a priori
claim (the cheap indicator), evaluates the cheapest rung whose claim meets the
request, and unfolds to the next when an a posteriori certificate fails
(`Assessment::Deferred` rungs, e.g. a solver at a given resolution).

**Audits and calibration.** With `Audit::Every(n)` the engine re-runs the
reference rung on the first and every `n`-th accepted result of each rung and
measures the real error. A non-rigorous rung that exceeded its claim gets a
safety factor (1.5 × the observed excess) on later claims; a rigorous rung
that exceeds its claim is counted as a defect and never calibrated away.

**Closures.** `Closure` folds a fine model once into a piecewise-Chebyshev
table: each patch samples `2n + 1` Lobatto points, the nested degree-`n`
interpolant is checked on the `n` held-out points, and the patch is halved
until the check passes. Two physics hooks make the claim honest:
`build_resolving` forces the node spacing below the narrowest feature the
samples cannot reveal (a resonance linewidth), and `build_inexact` adds the
fine model's own claimed error to every patch claim. Unconverged patches are
out of domain, so the ladder unfolds exactly there.

**Goal-oriented budgets.** For an observable `J = Σ w_i f(x_i)` (quadrature,
field sum, or linearised downstream sensitivity) `solve_functional` picks one
rung per item so that `Σ |w_i| e_i ≤ tol` at least cost: a multiple-choice
knapsack, solved greedily on each item's lower convex hull of
`(cost, |w| e)` (exact for the LP relaxation, within one upgrade of the
integer optimum; checked against brute force on 200 random instances).
`refine_functional` re-plans with new weights (pilot → final) and reuses every
item whose rung did not change.

## Results

### Engine, against exact references (`emergence/tests/engine.rs`)

| check | result |
|---|---|
| Taylor ladder of `exp` with the Lagrange bound, 6000 solves at three tolerances | every answer within its rigorous claim; the chosen degree is always the cheapest whose bound meets the tolerance |
| over-optimistic indicator (claims 1/10 of the true remainder) | first audit measures 10.8× the claim, safety becomes 14.4; zero violations after the first 1000 inputs |
| closure of `sin(10x) eˣ`, `1/(1+25x²)`, 20 000 random checks | claims hold; worst error/claim `5.7e-7` (tol 1e-6) and `9.9e-3` (tol 1e-10), i.e. the nested estimate is conservative |
| power law over 12 decades (log axis); a jump | claims hold; the jump's patch stays unconverged and reports out of domain |
| knapsack plan vs brute force | always feasible, never more than one upgrade above the optimum |
| pilot → refine reuse | identical values and rungs to a fresh solve, at lower cost |

### Mie scattering (`photonics/tests/mie_ladder.rs`)

Rungs: second-order Rayleigh (B&H eq. 5.7, conjugated to the crate's `n + iκ`
convention, checked against the exact `a₁, b₁, a₂`), resonance-gated closure
tables, and the Wiscombe series.

- **Per point**, 3000 random `(x, m)` per tolerance for soot, dust and water
  at 2.13 µm, every result audited against the series: zero violations; the
  ladder is 11.7× / 9.4× / 6.9× cheaper than the series at tolerance
  1e-3 / 1e-5 / 1e-7 (worst error/claim 0.79 / 0.67 / 0.02). A tolerance below
  the series' own accuracy (1e-10) is reported unmet.
- **Resonances.** Weakly absorbing spheres have resonances far narrower than
  any table spacing. `resonances()` finds every multipole resonance as a zero
  of `Im(1/c_n)` (a sign scan cannot miss it however narrow), refines it by
  Illinois regula falsi (identical set to brute-force bisection, 12 870
  resonances, 4× faster), and gets the FWHM from the Lorentzian slope (checked:
  half maximum at `x* ± Γ/2` within 2 %). For visible water at `60 ≤ x ≤ 64`
  (111 resonances, widths `3.7e-6 … 40`) a table built without this list is
  wrong by **3.9e6×** its own claim on the resonances; the gated table
  resolves all 111 (431 patches), or with a 48-patch budget resolves 22 and
  excludes 89, and the ladder returns every one within its claim either way.
  For water at 2.13 µm the gated closure holds its claim on all 472 resonance
  peaks in `5 ≤ x ≤ 40` (worst error/claim 0.04).
- **Cloud field**, 200 gamma distributions (`r_eff` 4–20 µm) × 2001 Simpson
  nodes at 2.13 µm: every cell matches exact Mie at every node to `4.7e-12`
  (claim ≤ 1e-5); **16 ms vs 174 ms** for the series (10.7×), plus a one-off
  0.2 s closure build (including the resonance scan).

### Grey phonon slab, `κ_eff/κ₀(Kn)` (`heat/tests/bte_ladder.rs`)

Rungs: Fourier; Fourier + Milne jump `1/(1 + 2ζ Kn)`, `ζ = 0.710446`
(claim `0.05 E₂(1/Kn)`: the bulk linear profile and each wall's Milne layer
are exact, only their interaction through the slab is not); ballistic
`(3/4)/Kn`; a closure in `ln Kn`; the integral equation, Richardson-extrapolated
(measured order 1.77, from the log-singular `E₁` kernel) at refinements 1, 2, 4.

- Claims on an independent 48-point sweep (`0.2 ≤ Kn ≤ 300`): worst
  error/claim 0.27 (jump), 0.52 (ballistic). The reference's a posteriori claim
  holds against a finer extrapolation with ≥ 9× margin (e.g. `Kn = 0.3`:
  `2.7e-6` vs claim `2.5e-5`).
- 1500 `Kn` over `[1e-6, 1e6]`: all met. Tolerance 1e-3 uses 11 exact solves;
  1e-5 uses 238, of which 75 escalated to the ×2 grid and 3 to ×4 because the
  coarser certificate failed. The integral equation everywhere would cost
  `4e23` ns (`Kn = 1e-6` needs `m = 2e7`); the ladder spends `8e10`.
- **Film conductivity** (heat kinetic spectrum → slab ladder, per-mode grey
  suppression): silicon at 300 K, `L = 0.1 / 1 / 10 µm` gives
  `21.8 / 66.4 / 103.2 W/(m K)` ± 0.013 (tolerance `1e-4 κ_bulk`). The
  goal-oriented budget is **3.6e5× cheaper** than a uniform per-mode tolerance
  with the same bound, because a uniform tolerance forces integral-equation
  solves on transitional modes that carry almost no heat.

### Rigorous certified resolution (`spectral-core/tests/certified_ladder.rs`)

Each rung runs the spectral solver at one `N` with the rigorous weak–strong
certificate; the engine escalates `N` until the proof meets the tolerance.
Taylor–Green, Re = 200, `t = 0.5`, distances against a 48³ run:

| tolerance | chosen | rigorous bound | true distance |
|---|---|---|---|
| 1e-2 | 12³ | 9.6e-3 | 1.9e-3 |
| 1e-3 | 16³ | 3.5e-4 | 8.7e-5 |
| 3e-5 | 24³ | 1.8e-5 | 5.1e-6 |

A tolerance of 1e-5 is reported unmet with the 24³ answer. The kinetic
energy comes with its rigorous bound `δ‖v‖ + δ²/2`.

### Cross-crate: Mie → radiation (`testlab/tests/emergence_links.rs`)

A 12³ cloud (200 m cube, `N` ≈ 1e7 m⁻³, `r_eff` 5–15 µm) at 2.13 µm, 80 rays.
Every (cell, size) item is budgeted by its adjoint weight
`max_p T_p ℓ_pc N_c w π r² pdf` from a pilot solve; the per-ray bound
`T_p (e^{δτ_p} − 1)` is then computed from the final per-item claims and the
reference's own truncation error. Over 362 504 items (3626 Rayleigh, 358 878
closure, 0 series), the worst transmittance error against exact Mie through
the same `radiation::Volume` is **2.3e-14** with every bound ≤ 1e-6;
**38 ms vs 127 ms** for the series.

## Engine performance (measured, not assumed)

The first version made the cloud field *3× slower* than calling the series
directly (641 ms vs 222 ms): one statistics mutex locked by every assess and
evaluation across 16 threads. Fixes, each measured:

| change | effect |
|---|---|
| per-rung atomics; batches tally locally and hash-sample audits | 641 → 164 ms |
| parallelise over cells (nested rayon) instead of tiny per-cell batches | 164 → 18 ms |
| planner fast path (no hulls/heap when the cheapest choices meet the budget); parallel reductions; no intermediate copies | 362k-item functional 18.8 → 8.9 ms |
| closures stored flat, bucket index, Chebyshev coefficients + Clenshaw | lookup 115 → 55 ns |
| resonance refinement by Illinois instead of 64-step bisection | 371 → 92 ms |

The remaining engine overhead is about 35 % over the bare rung evaluations on
362k-item batches.

## Honest misses and limits

- **Closures are one-dimensional.** `Q(x)` at fixed `m`, `S(Kn)`. Multi-
  parameter closures (`Q(x, m)`, `η(ρ, T)` from MD) need sparse/tensor tables.
- **Visible-water integrals.** A per-point closure cannot resolve resonances
  of width ~1e-6 at large `x` within a budget, so it is (correctly) out of
  domain there and the series is used. A size integral does not need them
  resolved point by point; a resonance-averaged closure would be the right
  fold and is not built.
- **Deferred rungs restart from scratch.** Escalating 8³ → 12³ → 16³ → 24³
  costs 38 % more than running 24³ directly at the finest tolerance. Predicting
  the needed resolution from the failed certificates' convergence is not
  implemented.
- **The radiation chain is not rigorous end to end.** The per-item claims are
  validated (closure) or asymptotic (Rayleigh); the propagation through
  `T = e^{−τ}` is exact in form.
- **Film conductivity assumes independent grey suppression per mode** and is
  not validated against measured films; it differs from the Casimir
  (Matthiessen) model by 10–20 % for thin films. The claim covers the ladder's
  error within the model only.
- **Asymptotic constants are fitted, then checked out of sample**: `0.05` in
  the Milne-jump claim (4× the fitted prefactor), `(1 + 0.1 ln(1+Kn))/Kn²` for
  ballistic. Audits catch a regime where they fail; they are not proofs.
- **Effective-medium indicators in `photonics`/`phononic` are still switches**
  and are not on the engine yet; neither are `compressible`, `materials`, `md`,
  `electrochem`, `semiconductor`, or the GPU crates.

## Adding a rung

Implement `emergence::Rung<I, O>` (or wrap functions with `FnRung`): a cost in
the ladder's units, an `assess` that is cheap compared with `eval`, and either
an a priori claim or `Deferred` plus `certify`. Give the ladder a metric that
matches the claims, mark the finest rung as the audit reference, and validate
every claim against the reference out of sample before relying on it (the
test files above are the template).
