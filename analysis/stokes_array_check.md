# Penalized sphere vs. exact Stokes drag in a periodic array

Exact result (Hasimoto 1959; Sangani & Acrivos 1982, simple cubic array):
drag = 6πμaU·K(c), K = 1/(1 − 1.7601 c^{1/3} + c − 1.5593 c² + …).
Sphere D = 2π/8 in the (2π)³ box: c = 1.022e-3, K = 1.2139, so at
Re = UD/ν = 7.854e-4 (U = 1e-3, ν = 1): **Cd = 24/Re·K = 37,094**.

GPU runs, N = 64 (8 cells per diameter), no sponge, steady state (t = 8):

| mask | η = 1e-3 | η = 3e-4 | η = 1e-4 | trend as η → 0 |
|---|---|---|---|---|
| tanh (infinite tails) | 47,311 | 55,412 | 63,120 | diverges: tails χ/η inflate the body |
| compact C^∞, layer on [R−δ, R+δ] | 36,257 | 39,143 | 41,396 | converges to sphere of radius R+δ (overshoot) |
| compact C^∞, layer on [R−2δ, R] | 25,345 | — | 30,000 | converges from below to the true R |

Findings:
1. Mask tails matter: with an infinitely supported mask the η → 0 limit
   does not exist. Fixed by a compactly supported C^∞ transition.
2. Where the smoothing layer sits sets the η → 0 geometry. Placing it inside
   the body makes the limit the true sphere.
3. With the correct mask the error falls like η^0.22, close to the η^{1/4}
   rate proved by Angot–Bruneau–Fabrie. That is too slow to buy accuracy
   by shrinking η with explicit time stepping, because dt ≤ 2.5η. Next step:
   implicit penalization (solve (I + hχ/η)u + h∇p = r, ∇·u = 0 with
   FFT-preconditioned CG), so η can go to ~1e-8 at the CFL time step.

## Implicit (IMEX) penalization, driven mean flow (2026-10-03)

Driven by a uniform pressure gradient G (drag = G·V exactly at steady
state, mean velocity measured), sharp mask, CPU f64:

| η | N = 32 | N = 64 | N = 96 |
|---|---|---|---|
| 1e-8 (→ 0) | 1.302 | 1.369 | 1.390 |
| 0.05 dx²/ν | 0.977 | 1.096 | 1.170 |
| 0.10 dx²/ν | 0.879 | 1.056 | 1.140 |
| 0.25 dx²/ν | 0.713 | 0.977 | 1.086 |
| 1.00 dx²/ν | 0.404 | 0.770 | 0.942 |

(K = F/(6πμaU), exact 1.2139.) Grid-matched families rise monotonically
toward the exact value but are pre-asymptotic at N ≤ 96 (fitted orders
0.4–0.6, see `stokes_convergence.out`). η → 0 converges to ~1.4, a wrong limit.

### Plane Poiseuille check (`channel_probe`, exact answer elementary)

| N | η = 1e-10 | η = 0.25 dx²/ν |
|---|---|---|
| 16 | −61.9% | +13.9% |
| 32 | −54.9% | +5.5% |
| 64 | −27.4% | +2.1% |

Grid-matched penalization is consistent (order ≈ 1.35). The η → 0 failure is
an over-determination: the 2/3-dealiased velocity is a trigonometric
polynomial with ~2N/3 modes per direction, which cannot vanish at all solid
grid points (more constraints than degrees of freedom) without being
crushed. Principled fix: impose the penalization on the M = 2N/3 solution
grid, where constraints and degrees of freedom match.

### Re = 100 sphere (GPU, implicit, η = 0.25 dx²/ν, D/L = 1/8, sponge)

Cd = 1.115 (64³), 1.149 (128³). Not directly comparable to the isolated
sphere value 1.087 because of blockage and periodic images at D/L = 1/8.

### Penalization on the solution grid (negative result)

Idea: impose the constraint on the odd grid M = 2k_max + 1, whose modes are
exactly the retained ones (one constraint per degree of freedom). Plane
Poiseuille, η = 1e-6, error against the exact flow *for the staircase
geometry each grid samples* (walls midway between fluid and solid points):

| N | solution grid M | full grid N |
|---|---|---|
| 32 | +23.5% (M = 21) | −19.6% |
| 64 | +10.5% (M = 43) | −4.6% |
| 96 | +7.2% (M = 63) | — |

The solution grid is under-constrained between collocation points and
converges from above; the full grid converges from below and is more
accurate at N = 64. It is kept opt-in (`Penalization::on_solution_grid`).
With η → 1e-10 the solution-grid system is numerically singular (CG stalls).

Re = 100 sphere, GPU, η = 0.25 dx²/ν: Cd = 1.115 / 1.149 / 1.119 at
64³ / 128³ / 192³ — not monotone (staircase sampling noise).

**Conclusion:** Fourier volume penalization is a first-order wall model in
every variant tried. Grid-matched η on the full grid converges cleanly for
plane walls (+2.1% at 64³); curved bodies show staircase noise of a few
percent. Higher wall accuracy needs a different boundary treatment
(smooth immersed-boundary kernels, or a body-fitted discretization).
