# Two-phase validation status, and liquid-metal MHD

## Two-phase (Cahn–Hilliard + variable-density NS), `capillary_probe`

Slab of liquid between gas layers, both interfaces perturbed by
a₀ cos 2x, a₀ = 0.03; interface width η = 1.2Δx.

### Rayleigh–Taylor (density ratio 857, physical ν_g = 15 ν_w)

The interface starts displaced at rest, so linear theory gives
a(t) = a₀ cosh(nt). The first version fitted the log-slope of a(t),
which is biased low (slope n·tanh(nt)); the cosh fit removes that bias.

| N | measured n | vs sharp theory 1.3548 | diffuse-profile theory |
|---|---|---|---|
| 32 | 1.2655 | −6.6% | −2.95% |
| 64 | (log-slope fit: 1.2200) | — | −0.72% |

"Diffuse-profile theory" is the exact inviscid Rayleigh eigenvalue for
the solver's own tanh density profile (`rt_diffuse.py`): about half the remaining
N = 32 error is the finite interface thickness, a model error that
vanishes as η → 0.

### Capillary waves: variable-density projection implemented

A planar capillary wave of the two-fluid slab (a₀ cos 2x, a₀ = 0.03) is
compared with ω² = σk³/(ρ₁+ρ₂) (for the sinuous slab mode this differs
from the single-interface value by < 0.5% at kh = 2π).

**Getting the mode right matters.** Zero-crossing timing of the raw
interface height and a single-mode damped fit both pick up the internal
interface relaxation and a non-oscillatory transient, and this is what
produced most of the earlier scatter. The frequency must be measured from
the interface-shift projection, `ζ = −∫δφ ∂zφ₀ dz / ∫(∂zφ₀)² dz`, and fit
with the transient absorbed (`h = e^{−γt}(a cos ωt + b sin ωt) + c₀ + c₁t`;
`robust_fit` in `tests/capillary.rs`). That estimator removes the
worst of the scatter, but its absolute value still drifts by several
percent with the fit window (the projected signal has an
order-parameter-slip transient); the deterministic initial-acceleration
estimator of “The residual is the finite interface thickness” below is
the one used for the error budget.

**The variable-coefficient projection.** `Projection::Variable` in
`twophase.rs` replaces the constant-coefficient Dong & Shen (2012)
splitting with the exact variable-density projection. Given the BDF2
momentum numerator `u*` (including the implicit ν₀ solve, `D = γ₀/Δt +
ν₀|k|²`), it solves

    A p = −∇·u*,   A p = −∇·(ρ⁻¹∇p),   x₀ = 2pⁿ − pⁿ⁻¹,

by preconditioned conjugate gradients, then sets `u^{n+1} = (u* −
ρ⁻¹∇p)/D`. `A` is applied pseudo-spectrally with the solver's usual 2/3
(dealised) mask: `∂_d p` in Fourier → physical → multiply `1/ρ` →
forward → divergence, six real FFTs per matvec. The preconditioner is
the constant-coefficient FFT Poisson solve `ρ_ref/|k|²` with `ρ_ref` the
volume harmonic mean of ρ. The conjugate-gradient inner product is the
weighted half-spectrum (real-field) inner product, so `A` is symmetric
positive definite. The path is selected with `s.set_projection(...)` or
`VCPROJ=1`; the default remains `Constant`.

**The projection weight is *not* the dominant dispersion error.** The
two standard capillary force forms differ by
`μ_c∇φ − (−φ∇μ_c) = ∇(φ μ_c)`, a pressure gradient. The previous version
of this note attributed the constant-projection error to the wrong
removal weight `1/ρ₀`. Direct measurement contradicts that attribution:
the discrete variable projection absorbs that gradient to **0.67 %** of
its norm at η = 1.2Δx and to 0.046 % once η ≥ 2Δx, and the *solenoidal*
part of `μ_c∇φ` (the part that actually drives the flow) is essentially
the same under both projections — 64.5 % (constant) vs 64.9 % (variable)
at ratio 10. Consequently the dispersion is nearly projection-independent.
The remaining error is **not** a force-discretisation/aliasing effect:
the projected force forms differ by only a few percent, and the residual
is independent of the interface resolution η/Δx at fixed physical η and k.
It is the physical finite-interface-thickness correction, measured below.
The variable projection is nonetheless the correct projection and does
help when the density contrast is large.

Errors (N = 24, η = 1.2Δx, ν small, M = 1e-8 so Cahn–Hilliard is frozen,
2-period damped-sinusoid fit, single thread; exact ω from σk³/(ρ₁+ρ₂)).
The fit window is short and the displaced tanh profile is not an eigenmode,
so the absolute values scatter by a few percent with run length (the
ratio-1 match of the two force forms is the meaningful statement); the
residual itself is the finite-interface-thickness effect quantified
below.

| ρ₁/ρ₂ | exact ω | μ_c∇φ / Constant | μ_c∇φ / Variable | −φ∇μ_c / Constant | −φ∇μ_c / Variable |
|---|---|---|---|---|---|
| 1 | 2.000 | −2.8% | −2.8% | −2.5% | −2.5% |
| 10 | 2.697 | +6.0% | +6.0% | +1.3% | +1.3% |
| 100 | 2.814 | +5.8% | +7.1% | +31.5% (unstable) | +25.7% |
| 1000 | 2.827 | −16.8% | **+11.4%** | not tested | unstable |

At ratio 1 the variable path is *identical* to the constant path (one CG
iteration, residual 5e-16), as it must be. At ratio 1000 it improves the
default form from −16.8 % to +11.4 %. At ratios 10–100 it neither fixes
nor worsens the default form, because the projection weight is not the
leading error there (the finite-thickness correction and the force-form
difference dominate); the balanced `−φ∇μ_c` form is better at ratio 10
but unstable by ratio 100–1000. A grid study at ratio 10 (N = 32, longer
run) gives constant/default −31.9 %, variable/default −32.3 %,
variable/balanced −1.9 %, confirming both the projection-independence and
the poor convergence of the default form.

### The residual is the finite interface thickness

The clean, deterministic way to measure the restoring stiffness without
fitting a time series is the **initial acceleration** of the projected
interface displacement: with u = 0 and a small mode a₀ cos kx, a few BDF1
steps give ζ̈(0), and ω² = −ζ̈(0)/ζ(0), Richardson-extrapolated in Δt.
This is a functional of the initial condition and the discrete operators
only, so it is reproducible to the last digit (unlike a damped-sinusoid
fit, whose residual here is ~25 % of the amplitude because the displaced
tanh profile is not an eigenmode — the order parameter must first slip —
and the long-time fit is correspondingly biased).

Two experiments separate physics from discretisation:

* **Fixed physical η and k, varying η/Δx.** With η = 0.31416 and k = 2
  (ηk = 0.6283) at N = 24 (η/Δx = 1.2) through N = 64 (η/Δx = 3.2, a
  2.7× finer interface), the restoring-stiffness error is constant:
  −35.28 %, −34.83 %, −34.79 %, −34.79 %. A dealiasing or projection
  inconsistency tied to the interface resolution would move; this does not.
* **Fixed η/Δx = 1.2, varying ηk (N = 24, 32, 48, 64 at k = 2).** The
  error scales roughly linearly with ηk and vanishes as ηk → 0:

| N | ηk | measured (ω² error) | thin-interface overlap |
|---|---|---|---|
| 24 | 0.6283 | −35.3 % | −29.5 % |
| 32 | 0.4712 | −28.2 % | −23.7 % |
| 48 | 0.3142 | −19.5 % | −16.9 % |
| 64 | 0.2356 | −14.8 % | −13.2 % |

The measured asymptotic slope is d(ω²/ω²_sharp − 1)/d(ηk) → **−0.627**
as ηk → 0 (0.5616, 0.5981, 0.6211, 0.6272 over the rows above).

**Analytic finite-thickness correction.** Linearising the Cahn–Hilliard
capillary force `μ_c∇φ` about the equilibrium tanh profile (thickness
w = √2 η) gives φ₁ = −ζφ₀′ and, using Lφ₀′ = 0 for
L = −∂z² + (3φ₀²−1)/η²,

    μ_c,1 = λk²φ₁ = −λk²ζ φ₀′,   δf = μ_c,1 ∇φ₀ = −λk²ζ φ₀′² ẑ,

so the force density is ∝ φ₀′(z)². The initial acceleration is the
incompressible response to that normal force; the streamfunction Green's
function e^{−k|z|}/(2k) makes the restoring stiffness

    ω²/ω²_sharp = ∫ e^{−k|z|} φ₀′(z)² dz / ∫ φ₀′(z)² dz
                = (3/2)∫₀^∞ e^{−√2(ηk)u} sech⁴u du
                = 1 − 0.6267 (ηk) + O((ηk)²),

exactly the thin-diffuse-interface result of Jacqmin (1999), *J. Comput.
Phys.* **155**, 96. The predicted coefficient −0.6267 matches the
measured −0.627. This is model physics: a sharper interface is the only
remedy, and it is why the leading-order discrete scheme is *consistent*
— the error is O(ηk) in the physical interface thickness, not in Δx.

The discrete energy budget confirms the pairing. For the inviscid,
M = 0 solver the total energy (kinetic + Cahn–Hilliard free energy) is
non-increasing, with a relative drift of −1.3×10⁻⁶ (N = 24) and
−1.9×10⁻⁶ (N = 32) over 300 BDF2 steps: the capillary force is the
variational derivative of the discrete free energy, so there is no
operator/dealiasing mismatch between μ_c in the Cahn–Hilliard update and
in the momentum force. Together these rule out the remaining suspects
listed at the top of this note (projection weight, force pairing, energy
inconsistency); the earlier “grid dependence” of the default form was a
short-run damped-fit artifact.

Tests `tests/capillary_thickness.rs` (< 5 s) pin the two facts above:
resolution invariance at fixed physical η and k, and agreement (within
4 %) of the measured correction with the overlap integral at
ηk = 0.236. `tests/energy_consistency.rs` (< 11 s) asserts the discrete
energy law. The old `scale_probe`/`force_response`/`cap_series` scratch
probes are removed.

**Cost.** CG iterations per step (N = 24, η = 1.2Δx, tolerance 1e-8,
mean final residual ~1e-9): ratio 1 → 1, ratio 10 → 13–15, ratio 100 →
38, ratio 1000 → 72. Wall time on the same N = 24 ratio-10 trajectory is
≈ 2.6× the constant path at ratio 10 (single thread), rising with the
iteration count at higher contrast. The constant-coefficient
preconditioner does not capture the 1/ρ contrast, so a better
(variable-coefficient) preconditioner is the natural next step.

**Static droplet.** A cylindrical (2D-in-3D) drop, R = π/2, ratio 10,
η = 1.2Δx, 120 steps of Δt = 0.01 gives Δp = +2.7 % of σ/R and a
parasitic-current magnitude |u| ≈ 1.7e-2 with *either* projection (the
two agree to four digits); the radial capillary force is an almost pure
gradient, so this benchmark is — correctly — projection-independent. The
+2.7 % offset is the diffuse-interface/curvature discretisation.

**Should it be the default?** Not yet, not universally: at ρ₁/ρ₂ ≲ 10
the constant path is exact to the projection order and ≈ 3× cheaper, and
the variable path does not improve the (separately limited) dispersion.
It should be enabled by default for ρ₁/ρ₂ ≳ 10, where it removes the
constant-coefficient weight error (ratio 1000 row) and is the physically
correct operator. The capillary dispersion is not a projection or
force-pairing defect (see “The residual is the finite interface
thickness” above); it needs a thinner interface, and the emerging
compressible multiphase engine will revisit the interface treatment
anyway.

Tests. Unit tests next to the implementation
(`twophase::vc_absorption_tests`, < 1 s) check: ratio-1 equivalence and
one CG iteration; discrete-gradient absorption (0.67 % at η = 1.2Δx,
0.046 % at η ≥ 2Δx); non-triviality vs the constant path; projected force
forms agree to < 15 %. `tests/vc_projection.rs` (< 10 s) covers the
static droplet, the one-step force-form difference, and the CG cost vs
density ratio. The capillary dispersion is measured by the env-driven
`diag_capillary` probe, which prints all estimator variants and both
projections (`VCPROJ=1`); it is `#[ignore]`d because a full robust-fit run
is ~30 s.

## Liquid-metal MHD (quasi-static, low Rm), `mhd.rs`

On the torus with uniform σ the current has a closed form,
J = σ P[u×B] + ∇Δ⁻¹s + J̄. For uniform B the Joule term is diagonal,
−(σ/ρ)(k̂·B)²û, and can be put into the integrating factor exactly.

Exact tests (`tests/mhd.rs`):

| test | error |
|---|---|
| Moffatt anisotropic Joule decay, explicit | 1.6e-9 |
| same, Joule term inside the integrating factor | 1.0e-10 |
| electrode-driven flow u_y = −(SB/ρν)(1 − e^{−νt}) sin x | 2.6e-14 |

### Homopolar mercury motor (`mercury_motor`)

Annulus between coaxial electrodes, axial B, radial current: exact
laminar profile u = (K/2ν) r ln r + A r + C/r.

| walls | N | rel. L2 error |
|---|---|---|
| Brinkman (implicit, sharp mask) | 32 | 43% |
| Brinkman | 64 | 16% |

That is the expected first-order Brinkman convergence, with the liquid
slipping about one cell into the electrodes. Coaxial-cylinder marker
walls with the smooth extension were added for this (`Shape::CylinderZ`,
`Surface::cylinder_z`); results are below once the runs finish.

## Which model holds: `flow_numbers`

| flow | Ma | cavitation σ_c | Rm | Ha | τ_e U/L |
|---|---|---|---|---|---|
| mercury motor, 5 cm, 5 cm/s, 0.3 T | 3e-5 | 6e3 | 3e-3 | 392 | 9e-18 |
| tap-water stream, 1 cm, 1 m/s | 7e-4 | 2e2 | 7e-14 | 0 | 1.3e-2 |
| ship propeller tip, 0.5 m, 30 m/s | 2e-2 | **0.24** | 1e-10 | 0 | 8e-3 |
| nitrogen jet, 1 cm, 500 m/s | **1.43** | ∞ | 6e-21 | 0 | 4e8 |
| liquid nitrogen at boiling point | 6e-4 | **0** | 6e-24 | 0 | 6e5 |

What follows:
- pure water sits in the leaky-dielectric regime (charge relaxation is
  comparable to the flow time), so its electrical behaviour needs the
  Taylor–Melcher model, not a conductor or insulator limit;
- the mercury motor is safely low-Rm, but its Hartmann layers
  (L/Ha ≈ 0.13 mm) are sub-grid, so walls need a Hartmann-layer model;
- the supersonic jet and the cavitating propeller break the
  incompressible model and need the compressible solver.
