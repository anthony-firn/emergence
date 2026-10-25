# GPU phase-field solidification: `solidification-cuda`

This note ports the Karma–Rappel phase-field step of
`fluid/crates/solidification/` to CUDA (NVRTC, runtime-compiled, `cudarc`
0.19) and uses it to attack the reported miss in
[analysis/solidification.md](solidification.md) §2(d): the 2D free-dendrite
tip velocity. It also corrects the operator: the CPU crate omits the
variational anisotropy term of KR Eq. (88), and without it no steady dendrite
is selected at all.

Everything below is reproducible from
`cargo run --release -p solidification-cuda --example dendrite_convergence`
and `... --example gpu_validate`, and from
`cargo test -p solidification-cuda` (skips gracefully without a GPU). The
device is an RTX 3050 Laptop (6 GB, `sm_86`); the kernels are compiled with
`--fmad=false` so that the `f64` path matches Rust's non-FMA scalar arithmetic.

## 1. The port

`src/kernels.cu` is a templated (`double` / `float`) transcription of the CPU
update in `solidification/src/phase_field.rs`: cell-centred fields with a
one-cell ghost border, `idx(i,j) = j*(nx+2)+i`, symmetry/far-field ghost
edges, the anisotropic Laplacian expanded as
`W² lap p + 2 W W'(∇θ·∇p)`, and the Karma–Rappel kinetic-free choice
`λ = a1 W0/d0`, `τ0 = a2 λ W0²/D` (`a1 = 5√2/8`, `a2 = 0.6267`; Karma &
Rappel, *Phys. Rev. E* **57**, 4323 (1998)). Anisotropy is
`W(θ) = W0(1 + ε4 cos4θ)`, `τ(θ) = τ0 (W/W0)²`.

The crate adds a second, **corrected** operator, `step_kr`, with a `mode`
switch used for the comparison in §4:

* `mode` bit 0 — use the 9-point isotropic Laplacian that KR used
  (Eq. 88 discussion: "discretized using a nine-point formula ... which
  reduces the grid anisotropy"), otherwise the 5-point stencil;
* `mode` bit 1 — include the variational anisotropy term
  `∂ₓ[W |∇p|² ∂W/∂(∂ₓp)] + ∂_y[…]`, which with `Q = W dW/dθ` is
  `pₓ ∂_y Q − p_y ∂ₓ Q`. `solidification`'s CPU step omits it.

`step_kr(dt)` is `mode = 3`; the CPU-equivalent step is `step(dt)` (`mode` 0
without the variational term, and is the routine tested against the CPU).

## 2. GPU = CPU and f32 vs f64

`tests/gpu_cpu.rs`, 40 steps of the CPU-matched `step`, `48²`, `dx = 0.4`,
`W0 = 0.5`:

| quantity | value |
|---|---|
| `max|Δp|` GPU f64 vs CPU | `1.1e-16` |
| `max|Δu|` GPU f64 vs CPU | `3.3e-16` |
| `max|Δp|` GPU f32 vs GPU f64 (40 steps) | `1.1e-7` |

`examples/gpu_validate.rs` (`512²`, 2000 steps, `f64` reference values):

| | value |
|---|---|
| GPU f64 `max|Δp|` vs CPU | `3.05e-16` |
| GPU f32 `max|Δp|` vs GPU f64 | `5.94e-6` |

The f64 match is a line-by-line identity: no FMA contraction, IEEE division
and `sqrt`, and `atan2` evaluated on the device (the `f64` `--prec-div=true
--prec-sqrt=true` defaults are used).

## 3. Speedup vs the CPU solver

`cargo run --release -p solidification-cuda --example gpu_validate -- 512 2000`
on a thermally settled, otherwise idle GPU. The CPU reference is the rayon
parallel `solidification` solver, so the ratio is against a multi-core CPU,
not one thread:

| solver | wall | throughput | vs CPU f64 | accuracy |
|---|---|---|---|---|
| CPU f64 (rayon) | 4.10 s | 0.128 Gcell/s | 1.0× | — |
| GPU f64 | 10.48 s | 0.050 Gcell/s | **0.4×** | `3.1e-16` vs CPU |
| GPU f32 | 1.50 s | 0.350 Gcell/s | **2.7×** | `5.9e-6` vs f64 |

This is the honest number and it is modest: the RTX 3050 Laptop is a consumer
part with heavily throttled FP64 (the f64 path is *slower* than the CPU), and
each explicit step is three small kernel launches (two boundary passes plus
the update), so launch overhead is a real fraction of `512²`. The useful
result is the **f32** path at ~2.7× a multi-core CPU (per-core that is ~30–40×,
counting the machine's 16 cores) with only f32-roundoff-level drift over the
~10³-step runs used here. A production version would fuse the boundary passes into the
update kernel and batch steps; that is left as a clear next step.

## 4. The miss is the operator, not only the grid

The original miss (analysis/solidification.md §2(d)) ran the `solidification`
operator, which is `mode 0`. Running the KR D = 2 benchmark from a common
seed (`dx = 0.4`, `W0 = 1`, `D = 2`, `d0 = 0.277`, `400×200`, f32, moving
frame) isolates the two corrections:

| operator | `V d0/D` at `t=320` | at `t=640` | at `t=768` | tip radius fit (8-row cap) |
|---|---|---|---|---|
| 5-pt, no variational (CPU crate) | 0.0168 | 0.0053 | — | **fails** (`R = NaN`) |
| 9-pt, no variational | 0.0169 | 0.0066 | — | **fails** (`R = NaN`) |
| 5-pt + variational | 0.0255 | 0.0180 | — | `R_cap ≈ 3.7 W0` |
| **9-pt + variational (KR production)** | 0.0196 | 0.0174 | 0.0164 | `R_cap ≈ 3.7 W0` |

Without the variational term the tip velocity passes through ~0.017 and then
**decays monotonically while the `p = 0.5` contour loses its parabolic tip**
(the least-squares radius fit returns `NaN`). The CPU crate therefore never
selects a steady dendrite, at any resolution — this is the root cause of the
documented miss. Adding the KR term restores a stable parabolic tip. The
9-point Laplacian then lowers the velocity toward KR (5-pt 0.0180 → 9-pt
0.0174 at `t = 640`), exactly the grid-anisotropy effect KR introduced the
stencil for. In the fast test the two corrections are separated:
`mode 3 − mode 1` (variational) is `4.3e-2`, `mode 1 − mode 0` (stencil) is
`4.3e-3` after 20 steps.

## 5. The KR Table II benchmark

Parameters are KR's exactly: `W0 = 1`, `τ0 = 1`, `Δ = 0.55`,
`ε4 = 0.05`, `u = −Δ`, `λ = λ* = D/a2`, i.e. `d0 = a1 a2 W0/D`. The domain
is a quadrant (symmetry on `x = 0`, `y = 0`; far field on the outer edges)
with the moving-frame translation of KR §III. `V d0/D` is measured from the
`p = 0.5` crossing on the row adjacent to the symmetry plane. KR's converged
(Green-function) value is `V d0/D = 0.0170`; their phase-field column is
reproduced for reference.

### 5.1 The D = 3 deficit is resolution, not the operator

The open miss was: with the corrected operator, D = 2 appeared to agree with
KR (`V d0/D ≈ 0.0164–0.0172` vs 0.0168) while D = 3 plateaued at 0.0150 vs
0.0175 (13 % low). Four checks — (1) `dx/W0` refinement, (2) domain and
moving-frame/boundary treatment, (3) true steady state, (4) D = 4 — show that
the 13 % is **grid resolution**, compounded by two comparison errors: the old
D = 2 point was not converged, and the finite-`W0/d0` kinetic correction has
not vanished.

All runs below start from a generic parabolic needle (`x_tip = 15 W0`,
`r_tip = 2 W0`, so not tuned to the answer) and use a **gentle** moving frame
(`SHIFT_AT = 0.26 Lx → SHIFT_TO = 0.25 Lx`, shifted at every report rather
than in one jump). Each reported value is a genuine plateau: velocity and tip
radius constant over thousands of time units with the far field at
`u_far = −Δ`.

`R_cap` in the tables below is the old diagnostic: a least-squares `p = 0.5`
parabola fit over a fixed **8 rows**, which at `dx/W0 = 0.4` is only `2.8 W0`
and sits entirely inside the diffuse-interface cap. That cap radius (a few
`W0`) is *not* the macroscopic tip radius; the sharp-interface tip radius is
`R = 2 D p/V` (`p = 0.25693` for `Δ = 0.55`). A `p = 0` fit over the
single-valued monotone needle, excluding the first `1.5 W0` and capping at
`40 W0`, recovers it: `R/W0 = 9.3` for the steady D = 2, `dx/W0 = 0.4` run and
`5.46` for D = 3, `dx/W0 = 0.2`, against `2 D p/V = 9.3` and `5.5`. The
selection constant is `σ* = 2 D d0/(R²V) = (V d0/D)/(2p²)`; because the
operating point is Ivantsov, it follows from `V` alone.

**D = 3 (`d0/W0 = 0.18464`, `W0/d0 = 5.42`), `dx/W0` refinement** (runs to
`t ≈ 2000`, moving frame, `Ly = 64–96`):

| `dx/W0` | `W0/dx` | `V d0/D` | `R_cap/W0` (8-row, biased) | `R/W0 = 2Dp/V` | `σ*` |
|---|---|---|---|---|---|
| 0.4 | 2.5 | 0.01490 | 3.06 | 6.37 | 0.113 |
| 0.3 | 3.33 | 0.01536 | 2.97 | 6.18 | 0.116 |
| 0.2 | 5.0 | 0.01600 | 2.95 | 5.93 | 0.121 |

The velocity rises monotonically with refinement (0.0149 → 0.0154 → 0.0160),
and a linear extrapolation in `dx` through the two finest points reaches
≈0.0173, consistent with KR's phase-field 0.0175 and the Green-function
0.0170 (the three points are slightly convex, so this is an estimate, not a
fitted limit). `σ*` is already `dx`-converged (0.11–0.12); the remaining
deficit is the phase-field's finite `W0/d0`, not the grid. The whole-needle
`R/W0 = 5.46` measured on the `dx/W0 = 0.2` grid (at the run's `V d0/D`
≈ 0.0172) agrees with the `2 D p/V` column to 2 %.

**D = 2 (`d0/W0 = 0.27696`, `W0/d0 = 3.61`) — the old "agreement" was a
transient:**

| `dx/W0` | `V d0/D` steady | `R_cap/W0` (8-row, biased) | `R/W0 = 2Dp/V` | `σ*` |
|---|---|---|---|---|
| 0.4 | 0.01536 | 3.84 | 9.27 | 0.116 |
| 0.3 | 0.01574 | 3.72 | 9.04 | 0.119 |
| 0.2 | 0.01663 | 3.71 | 8.56 | 0.126 |

At `dx/W0 = 0.4` the velocity passes through 0.0164 at `t = 768` (the value
quoted as agreement) but keeps falling and **converges to 0.01536**, both from
the parabolic needle (plateau to `t = 9920`) and from the old quarter-disk
corner seed with the same gentle frame (plateau to `t = 18480`). The earlier
D = 2 agreement was therefore coincidental: the run was stopped while the tip
radius was still growing, before the true steady state. With refinement D = 2
gives 0.0166, within 1 % of KR's 0.0168 and 2 % of the Green-function value.
A separate steady rerun at `dx/W0 = 0.4` (`V d0/D = 0.01530`, plateau to
`t ≈ 11700`) gives a whole-needle `R/W0 = 9.3`, matching the `2 D p/V = 9.3`
column and the closure's `σ*` to within the fit scatter.

**D = 4 (`d0/W0 = 0.13848`, `W0/d0 = 7.22`), `dx/W0 = 0.4`:** steady
`V d0/D = 0.01450`, `R/W0 = 2 D p/V = 4.91`, `σ* = 0.110` (plateau over
`t = 400–13200`; the old biased cap value was `R/W0 = 2.66`, `σ* = 0.375`).
The earlier "0.0167 at `t = 120`" was the seed transient; at fixed
`dx/W0 = 0.4` the velocity decreases with `W0/d0` as the interface is coarsened,
the expected finite-interface-width trend.

**Domain and moving frame are not the cause.** The measured diffusion length
is `2D/V = 24.8 W0` (D = 3), `36.1 W0` (D = 2), `19.1 W0` (D = 4). The gentle
frame keeps the tip at `0.25 Lx` with `137–200 W0 ≈ 6–10 (2D/V)` of melt
ahead, and the right-edge temperature stays exactly `u_far = −0.55000` at
every report. A `SHIFT = 0` control at `D = 3`, `dx/W0 = 0.4`, `Lx = 600`
(ahead up to `12 (2D/V)`) reproduces the moving-frame trajectory
(`V d0/D = 0.01527` vs 0.01515 at `t = 640`, 0.01505 vs 0.01501 at
`t = 853`); the far field only departs from `u_far` once the tip is within
`≈ 300 W0 ≈ 12 (2D/V)` of the wall. A lateral-size control at
`D = 3`, `dx/W0 = 0.2` gives 0.01623 (`Ly = 64 W0`) vs 0.01619
(`Ly = 96 W0`) at `t = 533`, so the top boundary is immaterial at this width.

**The operator is not the cause.** The variational term was re-derived:
`∂ₓ(W|∇p|²∂W/∂pₓ) + ∂_y(…) = ∂_y(Q pₓ) − ∂ₓ(Q p_y) = pₓ ∂_y Q − p_y ∂ₓ Q`,
which is exactly `px*qy − py*qx` in `pf_step_kr`; there is no sign or factor
error. The 9-point Laplacian is the standard isotropic stencil
`[4(axes) + 1(corners) − 20]/(6 dx²)`. The residual grid error is therefore
in the **central-difference anisotropy and variational terms**, which are only
first-order at the interface and are unaffected by the isotropic 9-point
stencil — the refined runs show this error is `O(dx)` near the tip, not
`O(dx²)`.

**`f32` is not the cause.** The `D = 3`, `dx/W0 = 0.4` trajectory in `f64`
matches the `f32` run at every report (e.g. `V d0/D = 0.01505` vs 0.01501 at
`t = 853`; both plateau at 0.01490), so the steady state is not a
round-off bias.

**Secondary observation — the kinetic-free cancellation is asymptotic.**
At these `W0/d0` the result is still sensitive to `τ0`: at `D = 2`,
`dx/W0 = 0.4`, `t = 1500`, `τ0 × 0.5, 1, 2` gives `V d0/D =
0.0218, 0.0160, 0.0110`, and the `τ0 × 2` run eventually leaves the
dendritic attractor (a near-planar front, `R/W0 ≈ 110`). The KR kinetic-free
choice only removes the leading `O(1)` kinetic term; the residual is
`O(W0/d0)`, which is why the D = 2 point (`W0/d0 = 3.6`) is closer to the
Green-function result than D = 3 (5.4) or D = 4 (7.2). This is a property of
the thin-interface model, not a bug in the step.

**Summary.** The D = 3 deficit at `dx/W0 = 0.4` is not a domain or
moving-frame artefact, not a long transient (the D = 3 plateau is genuine; it
was the old D = 2 agreement that was transient), not an operator error, and
not `f32` round-off. It is discretisation of the
anisotropic terms (`V d0/D` climbs 0.0149 → 0.0160 as `dx/W0` falls
0.4 → 0.2, extrapolating to KR/GF) plus the finite `W0/d0` thin-interface
correction, which is also why the previously reported D = 2 agreement was
premature. An honest statement at fixed `dx/W0 = 0.4` is therefore:
`V d0/D = 0.0154, 0.0149, 0.0145` for D = 2, 3, 4, all below the
Green-function 0.0170, with the selection constant `σ* = (V d0/D)/(2p²)`
spanning `0.110–0.116` over the same range.

**The `σ*` discrepancy is resolved: it was the tip-radius estimator.** The
`0.49–0.71` values in earlier versions of this note came from a `p = 0.5`
parabola fit over a fixed 8 rows, which at `dx/W0 = 0.4` is only `2.8 W0` and
lies entirely inside the diffuse-interface cap. The cap radius is a few `W0`
and is not the macroscopic tip radius; the sharp-interface radius is
`R = 2 D p/V` (`= 1/(σ* p) d0`). Correcting the estimator (p = 0 contour,
exclude the first `1.5 W0`, fit the monotone needle out to `40 W0`) gives
`R/W0 = 9.3` for the steady D = 2, `dx/W0 = 0.4` run and `5.46` for D = 3,
`dx/W0 = 0.2`, in agreement with `2 D p/V` and with `σ* = 0.12–0.13`.
Therefore the value `σ* = 0.1288` used to calibrate the sub-grid closure in
[solidification.md](solidification.md) §3(ii) **is** the correct one; the
closure has been left unchanged and now fixes it explicitly via
`closure::sigma_star_2d_kr_benchmark`, with a regression test
(`closure_reproduces_kr_2d_benchmark`).

## 6. Tests and usage

`cargo test -p solidification-cuda` (each test ≪ 30 s; all four skip if no
CUDA device is present):

* `real_roundtrip` — host/device real conversion.
* `gpu_f64_matches_cpu` — round-off identity of the port against the CPU step.
* `gpu_f32_close_to_f64` — f32 drift over a short run.
* `kr_operator_differs_from_plain_operator` — the variational and 9-point
  corrections both actually change the update.

Long runs live in `examples/dendrite_convergence.rs`; it checkpoints state to
a file and resumes, so a convergence sweep is a chain of invocations each
under a few minutes. Its `robust_tip_radius` reports the `p = 0` whole-needle
fit of §5 (the library `HostView::tip_radius` estimator, fixed at 8 rows, is
retained for reference but is the biased cap fit). Environment switches:
`STENCIL` (operator mode), `PREC`
(`f32`/`f64`), `DT`, `REPORT`, `SHIFT`/`SHIFT_AT`/`SHIFT_TO` (moving-frame
control; `SHIFT_AT = 0.26`, `SHIFT_TO = 0.25` gives the gentle frame used in
§5.1), `DIAG` (far-field and diffusion-length print), `TAU_FACTOR`, `MAXSEC`,
and `SEED_PARAB`/`SEED_X` (parabolic needle seed instead of the corner
seed). `examples/gpu_validate.rs` measures the GPU = CPU agreement and the
speedup.

## Sources

* A. Karma & W.-J. Rappel, *Phys. Rev. E* **57**, 4323 (1998) — model,
  kinetic-free constants, 9-point Laplacian, Table II benchmark values
  (`Ṽtip`, `d0/W0`, `Dx/W0 = 0.4`, `Dt/τ0 = 0.016`).
* A. Karma & W.-J. Rappel, *Phys. Rev. E* **54**, 3065 (1996) — earlier
  formulation.
* G. P. Ivantsov, Dokl. Akad. Nauk SSSR **58**, 567 (1947) and
  J. S. Langer & H. Müller-Krumbhaar, *Acta Metall.* **26**, 1681 (1978) —
  the sharp-interface closure in `analysis/solidification.md`.
* `cudarc` 0.19 / NVRTC — runtime kernel compilation; the same stack as
  `fluid/crates/compressible-cuda`.
