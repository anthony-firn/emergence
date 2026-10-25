# Overset spectral refinement patch — stage 1

A refinement patch is a separate periodic Fourier box (side Lp, scale
s = Lp/2π) whose fluid region is a sphere; in patch coordinates the
equations are standard Navier–Stokes with ν/s, τ = t/s and forcing s·f.
The sphere's boundary velocity is prescribed (here: the exact field; in
the coupled scheme: the coarse grid), imposed by the least-squares smooth
forcing extension in the fictitious buffer outside the sphere. Physical
forcing is windowed by a C^∞ step so that the total forcing stays smooth.

Test: steady forced ABC flow (exact NS solution), patch side 2.5, fluid
sphere radius 2 (patch units), extension shell 0.9, GPU
(`patch_probe_gpu`):

| patch N | physical resolution | max error inside / max|u| | wall solve | time |
|---|---|---|---|---|
| 48 | 0.081 | 9.6e-4 | 42–62 products | 2 min |
| 64 | 0.058 | 2.3e-4 | 42–62 | 6 min |
| 96 | 0.040 | 6.4e-5 | 42–58 | 35 min |

~4th-order convergence: the patch engine reproduces an exact flow from
boundary data alone. Cost is dominated by the host-side CGLS (40–60
products per stage, each with host↔device transfers of marker data);
moving the Krylov loop onto the device is the obvious speed-up.

Lessons on the way: dt = 10 with |u| ~ 1.5 was a CFL of ~150 (blew up);
the GPU IMEX path was not tracking max|u| for CFL (fixed); the extension
shell must be ≳ 5 cells (N = 24, 3.4 cells, diverged).

## Two-way coupling, attempt 1: relaxation (failed — kept for the record)

Coarse grid relaxed (explicit penalization) toward the patch velocity in a
ball r < 0.5 around the sphere (zero inside the body), patch outer
boundary fed by the coarse grid. Result (`runs/overset_relaxation_48_64.txt`):
the coarse mean flow keeps accelerating (K(momentum) 0.92 → 0.25 by
t = 15) while the patch body force just scales with it (K ≈ 15).
Momentum is not transferred: the coarse grid's own body is a 3-cell
penalized ball (rigid core ~1 cell after mask smoothing), and relaxing it
toward a patch solution that itself takes its boundary data from the
coarse grid is circular — neither side pushes back.

Fix: overlapping Schwarz (Dirichlet–Dirichlet, as in overset/Chimera
methods). The coarse grid gets a hole r < R_h whose surface velocity is
the patch solution (a surface with prescribed marker velocities; its wall
force carries the body's momentum sink), and the patch outer sphere
R_o > R_h takes the coarse velocity.

## Two-way coupling, attempt 2: overlapping Schwarz (works)

Coarse grid with a hole r < R_h (walls with prescribed velocity = patch
solution); patch outer sphere R_o > R_h fed by the coarse grid; exchange
every step. Small sphere a = π/8, ν = 10, exact K = 1.2140.

**Under-resolved overlap** (coarse 48³, R_h = 0.55 ≈ 4 cells with point
walls, overlap 1.5 coarse cells; patch 64³): converges to a *wrong* fixed
point — at t = 100, K(momentum) −14.7%, K(patch body force) +38.3%, and
F_body/(G·V) stuck at 1.62 instead of 1. The coarse grid must resolve the
overlap.

**Resolved overlap** (coarse 64³, R_h = 0.7 with a smooth-extension hole
wall; patch 96³, R_o = 1.0, 3× finer than the coarse grid):

| t | K(momentum) | K(patch body force) |
|---|---|---|
| 5 | −4.08% | +0.93% |
| 10 | −1.00% | +0.33% |
| 15 | −0.32% | **+0.16%** |

(momentum estimate still converging with the mean flow)

### Is it worth it here?

| method | K error | points | time |
|---|---|---|---|
| coarse 64³ alone (point walls) | +4.73% | 0.26 M | 2 s |
| uniform 192³, point walls | +1.38% | 7.1 M | 41 s |
| uniform 192³, smooth extension | −0.03% | 7.1 M | 53 s |
| overset 64³ + 96³ patch | +0.16% | 1.15 M | 25 min |

The coupling is consistent and 30× more accurate than the coarse grid
alone with 6× fewer points than the equivalent uniform grid, but much
slower: the lagged exchange needs ~300 small coupled steps where the
uniform steady solve needs 16 large ones, and a 6 GB GPU holds 192³
anyway. Overset refinement pays off where a uniform fine grid would not
fit (≥ 8× local refinement in a large domain), and would then need
sub-iterated exchange (large steps) and GPU-side transfers.
