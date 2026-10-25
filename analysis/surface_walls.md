# Surface-constrained walls (markers on the true geometry)

Implementation: `fluid/crates/spectral-core/src/{nufft,surface,imex,krylov}.rs`;
probe: `fluid/crates/spectral-core/examples/wall_probe.rs`. CPU, f64, driven
mean flow, steady state (IMEX steady states are dt-independent; checked).

## Plane Poiseuille (fluid-region mean; exact G h²/(3ν))

| N | point forces | smooth extension (w = 1.0, square GMRES) |
|---|---|---|
| 16 | −24.0% | +4.26% |
| 24 | — | +0.96% |
| 32 | −11.0% | +0.35% |
| 48 | — | +0.021% |
| 64 | −5.8% | −0.0048% |
| 96 | −3.9% | −0.0047% |

Point forces: first order (a δ-force puts a kink in u at the wall).
Smooth extension: ~4th order down to a floor ~3–5e-5 from leakage of the
truncated C^∞ bump (Fourier decay ~exp(−c√(k w)); w = 1.25 → −2.8e-5).
Wall slip 1e-7 of the centre-line speed; the error is a smooth profile.

## Sphere array (Hasimoto / Sangani–Acrivos K(c))

R = π/8 (K = 1.2139), shell 0.3 = 1.5–3 cells at N ≤ 64: the extension is
erratic (square GMRES stalls at 3e-4–2e-3; least squares +0.14/−4.8/−1.4%).

R = π/3 (K = 1.8314):

| N | point forces | smooth extension (w = 0.6, least squares) |
|---|---|---|
| 32 | +3.05% | −2.97% (slip resid 7e-4) |
| 48 | +1.27% | −2.21% (2e-4) |
| 64 | +0.29% | −1.99% (1e-4) |

On the sphere, point forces converge well (faster than first order here)
while the extension stalls near −2% — under investigation.

### Why the extension stalls on the sphere (resolved)

The least-squares extension reproduces the channel exactly (−0.0048% at
64³, slip 1e-8), so the machinery and adjoint are correct. Different
extension bases (σ = 1 vs 2) give *identical* K (1.79094), as they must:
any interior forcing meeting the same wall condition yields the same fluid
flow. Measured wall slip at N = 32:

| method | max slip / flow speed | K error |
|---|---|---|
| smooth extension (lsq) | 3.8e-3 (CGLS stalls at 7e-4 after 42 products) | −2.97% |
| point forces | 8.9e-9 | +3.05% |

On a curved body the wall velocity has tangential structure down to the
grid scale; a smooth shell basis at that resolution cannot control it, the
boundary solve stalls and the flow slips ~0.4%, lowering the drag.

**Recommendation:** point-force surface walls for curved bodies (exact wall
position, no-slip to 1e-12, monotone convergence), smooth extension for
planar / tangentially simple walls (spectral down to ~5e-5). Higher order on
curved walls would need singularity subtraction of the δ-force kink.

## Correction (supersedes "Why the extension stalls on the sphere" above)

Two errors in the sphere analysis above:

1. **Drag convention.** The surface probes computed K from G·(V − V_s).
   Sangani & Acrivos / Zick & Homsy define F = τ|∇⟨p⟩| over the whole cell
   (Darcy k = τ/(6πaK)), which includes the mean-pressure-gradient force
   G·V_s on the sphere — the convention the penalization probe already used.
   Correct K = old K / (1 − c), c = 0.0194 for R = π/3 (factor 1.0198).
2. **Slip ratio.** The "slip / flow speed" ratio divided by the velocity at
   (π, π, π), which is the sphere's centre — fictitious interior flow — so
   the 3.8e-3 figure was meaningless.

Corrected R = π/3 results (K exact 1.8314):

| N | point forces (GPU) | smooth extension, lsq (CPU) |
|---|---|---|
| 32 | +5.08% | −1.05% |
| 48 | +3.27% | −0.27% |
| 64 | +2.28% | **−0.05%** |
| 80 | +1.81% | |
| 96 | +1.51% | |
| 128 | +1.10% | |
| 160 | +0.87% | |
| 192 | +0.73% | |
| 256 | +0.54% | |

* Point forces: monotone, first order (128 → 256 halves the error);
  first-order Richardson from 192/256 gives **−0.03%**.
* Smooth extension: ~4th order on the sphere as well (21× from 32³ to 64³).
  It needs a resolved shell: on the small R = π/8 sphere (shell 1.5–3
  cells) it was erratic, which stands.

Revised recommendation: smooth extension wherever the body is thick enough
for a shell of ≳ 5 cells; point forces (+ Richardson) for thin bodies.

## GPU convergence, both wall methods (R = π/3, K exact 1.8314)

| N | point forces | smooth extension (w = 0.6, lsq) |
|---|---|---|
| 32 | +5.08% (0.6 s) | −1.05% (5 s) |
| 48 | +3.27% | −0.275% |
| 64 | +2.28% (1.8 s) | −0.054% (12 s) |
| 96 | +1.51% | **+0.0042%** (41 s) |
| 128 | +1.10% | +0.0067% (89 s) |
| 192 | +0.73% | +0.0027% (8.3 min) |

The extension reaches the 1e-5 level by 96³ (360× more accurate than
point forces at the same grid); beyond that it sits at ~3–7e-5, consistent
with the bump-leakage floor seen on the channel and f32 arithmetic.
Plot: `surface_convergence.png`.

## Cutoff (bump) shape study

Radial bump ψ_p(s) = exp(1 − 1/(1 − s²)^p), Fourier tail ~exp(−c k^{p/(p+1)}).
Plane Poiseuille, shell 1.0, fluid-mean error:

| p | 64³ | 96³ |
|---|---|---|
| 0.5 | +6.1e-4 | +2.2e-4 |
| 1 | −4.8e-5 | −4.7e-5 |
| 2 | −9.3e-5 | **+1.95e-5** |
| 3 | +2.2e-4 | +3.4e-5 |

Sharper bumps win once resolved: p = 2 lowers the 96³ floor 2.4×. Rule now
built in: p = 2 when the shell spans ≥ 12 cells, else p = 1 (`BUMP_P`
overrides). Auto: −4.83e-5 at 64³ (p = 1), +1.95e-5 at 96³ (p = 2).
