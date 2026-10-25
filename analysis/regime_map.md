# Regime map: where each model holds ("emergent boundaries")

`Solver::regime_map` / `continuum_check` (`regime.rs`). Example:
`regime_probe` (output in `regime_probe.out`).

## Which local indicator finds the error?

Taylor–Green Re = 400 at 32³ (under-resolved, min k_max·η = 0.22),
t = 3, against a 96³ reference; correlation of each indicator with the
local velocity error |u₃₂ − u₉₆| over all grid points:

| indicator | Pearson | Spearman |
|---|---|---|
| truncation-residual density (local certificate residual) | **+0.56** | **+0.56** |
| local Kolmogorov resolution 1/(k_max η) | +0.34 | +0.29 |
| vorticity magnitude (common refinement trigger) | −0.04 | −0.02 |

The certificate-derived indicator locates under-resolution; the common
vorticity heuristic carries no information here.

## Where the continuum itself fails

Non-equilibrium parameter λ|∇u|/c (collision time × strain rate), Mach
number, and mean free path / grid spacing, for the same flow in physical
units:

| case | λ|∇u|/c | Mach | λ/Δx | continuum |
|---|---|---|---|---|
| water, 1 cm box, 1 m/s | 6e-10 | 6e-4 | 1e-6 | valid |
| nitrogen gas, 1 cm box, 1 m/s | 6e-7 | 3e-3 | 2e-4 | valid |
| nitrogen gas, 1 µm box, 100 m/s | 0.59 | 0.26 | 2.1 | not valid |

(A velocity-based gradient-length Knudsen number was tried first; it blows
up at stagnation points and is not used.)
