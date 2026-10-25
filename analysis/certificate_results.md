# A posteriori error certificate and certified adaptivity

Implementation: `fluid/crates/spectral-core/src/certify.rs` (CPU, f64).
Bound: for the Hermite interpolant v of the numerical solution and any
Leray–Hopf solution u, d/dt‖u−v‖² ≤ 2μ‖u−v‖² + ‖R‖²_{H⁻¹}/(2ν),
μ = max(−λ_min(sym ∇v)). Residual = exact truncation part (3N/2 padded
grid) + Hermite time part (Simpson). Caveats: μ sampled on the padded grid,
Simpson quadrature, floating point.

## Validation (tests/certify.rs)

* ABC flow (exact solution): bound 1.6e-14.
* Under-resolved Taylor–Green, 16³, Re = 200, vs 64³ reference:
  bound ≥ true error at every checkpoint, **3.3–5.3× the true error**.

| t | true error | bound | ratio |
|---|---|---|---|
| 0.25 | 3.0e-6 | 1.6e-5 | 5.3 |
| 0.50 | 8.7e-5 | 3.4e-4 | 3.9 |
| 1.00 | 2.0e-3 | 6.5e-3 | 3.3 |
| 1.50 | 9.6e-3 | 3.3e-2 | 3.4 |

## Horizon of usefulness (Taylor–Green, Re = 400, 64³)

The time residual stays ~1e-8 (CFL 0.5 RK4 is ample); truncation grows to
2e-3 by t = 4. The bound is 3.5e-3 at t = 2.4 (< 1% of ‖u‖ ≈ 0.49) but 1.1
at t = 4: μ rises to 4.7 and exp(2∫μ) amplifies ~100× over t ∈ [3.5, 4].
This is the expected limit of any stability-based bound for flows that
develop strong strain.

## Certified adaptive resolution (Re = 400, t ∈ [0, 2.5], start 16³)

Refine when the residual's contribution to the bound, discounted by the
predicted amplification exp(2μ(T−t)), exceeds the budget; coarsening adds
the exact dropped norm to the bound.

| run | bound at t = 2.5 | point-steps |
|---|---|---|
| fixed 32³ | 6.3e-2 | 8.2e5 |
| fixed 48³ | 1.7e-2 | 4.1e6 |
| fixed 64³ | 5.0e-3 | 1.29e7 |
| fixed 96³ | 8.5e-4 | 6.46e7 |
| adaptive, target 5e-3 | 6.6e-3 (16→24→32→48→64→96) | 1.95e7 |
| adaptive, target 1e-3 | 5.3e-3 (16→32→…→128) | 5.49e7 |

Findings: for Taylor–Green, global adaptivity does **not** beat the best
fixed grid — the resolution need grows monotonically, so the fine late
phase dominates cost either way — and the controller undershoots its
target because μ keeps growing after the prediction. The certificate's
value here is diagnostic: it rigorously says whether a run was resolved
enough, within a factor of 3–5. Gains from adaptivity require flows whose
resolution need varies non-monotonically in time or (with local
refinement) in space.

## GPU certificate (`CudaCertifier`)

Same algorithm on the GPU (padded-grid truncation residual via cuFFT,
Hermite time residual, closed-form strain eigenvalues, f64 reductions).
Cross-check, Taylor–Green 24³, Re = 200, t = 1.5: truncation residual and
μ agree with the f64 CPU certificate to 4+ digits; bound 9.2909e-3 (GPU)
vs 9.2906e-3 (CPU).

Taylor–Green Re = 400 to t = 2.5 (`fluid tg --certify --backend cuda`):

| run | certified bound | time |
|---|---|---|
| CPU 96³ | 8.53e-4 | 33.2 s |
| GPU 96³ | 8.67e-4 | 5.0 s |
| GPU 192³ | 3.28e-4 | 76.7 s |

The f32 time residual has a roundoff floor (~1e-6 vs 1e-11 in f64); it
only inflates the bound (conservative), but it limits how small a GPU
certificate can get: 192³ improves on 96³ by only 2.6×, where resolution
alone would give far more. A certificate below ~1e-4 needs the residual
evaluated in f64.
