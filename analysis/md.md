# Molecular dynamics: the statistical fold from atoms to the continuum

Implemented in `fluid/crates/md` (new crate; it depends on `materials` for
the Lennard-Jones potential but does not modify it). This is the
statistical rung of [emergence_ladder.md](emergence_ladder.md): the fine
model is a system of classical atoms with a pair potential in a periodic
box; the coarse model is the Navier–Stokes fluid (or the elastic solid)
with an equation of state and transport coefficients; the fold between
them is

- an **equation of state** `P(ρ,T)`, `U(ρ,T)`, taken from the
  Johnson–Zollweg–Gubbins (1993) fit and reproduced by simulation;
- **constitutive closures** — shear viscosity `η`, thermal conductivity
  `κ` and self-diffusion `D` — computed from equilibrium time-correlation
  functions and returned in SI units for argon;
- a **validity indicator** — the Knudsen number and the ratio of the flow
  scale to the stress/heat-flux correlation length and time — that says
  when the continuum with those closures may be used, and where the
  molecular description has to be unfolded.

## The engine

Reduced Lennard-Jones units are used throughout: length `σ`, energy `ε`,
mass `m`, so temperature is `ε/k_B`, pressure `ε/σ³`, time
`σ√(m/ε)`, viscosity `√(mε)/σ²`, thermal conductivity
`k_B√(ε/m)/σ²` and self-diffusion `σ√(ε/m)`. The argon map
(`argon.rs`) uses `σ = 3.405 Å`, `ε/k_B = 119.8 K`, `m = 39.948 u`
(Rowley & Nicholson 1972, via the standard LJ tables used by `materials`).

| module | content |
|---|---|
| `potential.rs` | Lennard-Jones and Morse; three cutoff schemes: truncated-and-shifted (energy-continuous), shifted-force (energy- and force-continuous), and bare truncation with the standard homogeneous tail corrections to energy and pressure |
| `neighbor.rs` | cell-linked, minimum-image neighbour lists (full lists, so the force loop parallelises over atoms without races); rayon |
| `system.rs` | periodic cubic box, force / per-atom energy / stress / heat-flux kernels, velocity-Verlet integration; face-centred-cubic lattice constructor |
| `thermostat.rs` | Bussi stochastic velocity rescaling (exact canonical kinetic-energy distribution), Berendsen weak coupling, single-chain Nosé–Hoover |
| `barostat.rs` | Berendsen isotropic barostat (equilibration only) |
| `eos.rs` | the JZG (1993) modified Benedict–Webb–Rubin equation of state |
| `transport.rs` | Green–Kubo shear viscosity and thermal conductivity, the generalized Einstein (Helfand) viscosity, MSD/VACF self-diffusion, the Yeh–Hummer finite-size correction |
| `argon.rs` | reduced↔SI closure for argon and the measured liquid-argon reference values |
| `validity.rs` | Knudsen and scale-separation indicators |

The force on atom `i` from `j` is `u'(r) (r_j − r_i)/r`, so the virial
tensor is `Π = Σ_i m v_i⊗v_i − Σ_{i<j}(u'/r) r_ij⊗r_ij` and the pressure is
`tr Π/(3V)`. The heat flux follows the Irving–Kirkwood/LAMMPS convention,
`J V = Σ_i e_i v_i + (1/2) Σ_{i<j} (F_ij·(v_i+v_j)) r_ij`, with
`e_i = ½ m v_i² + ½ Σ_j u(r_ij)` and the pair term carrying the `1/2` of
the symmetrised sum; the implementation was checked against the LAMMPS
`compute heat/flux` source, where the per-atom virial already contains
`1/2`. Getting this factor right changed `κ` by a factor of four, from a
wrong `κ* ≈ 26` to the correct `κ* ≈ 6.5`, and is documented here because
it is the easiest place to be off by two.

## Benchmarks

All runs use 256 atoms (fcc), reduced units, a 2.5 σ cutoff, and the Bussi
thermostat unless stated. Errors are one standard error; where the
estimator's own scatter understates the true run-to-run error this is
noted. The numbers below are from `cargo test -p md`; the machine was
heavily loaded by other work, so the wall times are not representative.

| benchmark | measured | reference | note |
|---|---|---|---|
| NVE energy fluctuation, `dt = 0.004` / `0.002` | ratio `3.9` | `4` | velocity-Verlet is second order; same initial condition, fixed physical time |
| JZG EOS, `T* = 2.0064`, `ρ* = 0.5` | `P* = 1.0434 ± 0.0370`, `U_res* = −3.1447 ± 0.0091` | `P* = 1.0775`, `U_res* = −3.1449` | JZG 1993 Eq. (7), (9) |
| JZG EOS, `T* = 1.0064`, `ρ* = 0.8` | `P* = 1.0165 ± 0.0427`, `U_res* = −5.5410 ± 0.0081` | `P* = 1.0318`, `U_res* = −5.5235` | same |
| fcc zero-temperature energy, `ρ* = 1` | `E/N = −8.2970212419` | `−8.2970212419` (materials) | exact to `3.2e−14` |
| fcc low-temperature, `T* = 0.0513` | `<PE> = −8.600637`, `<E> = −8.523964` | `−8.600758`, `−8.523783` | static energy + `(3/2)T` (harmonic equipartition) |
| radial distribution function, `ρ* = 0.8`, `T* = 1` | first peak `r = 1.125 σ`, `g = 2.76`, coordination `11.25` | LJ liquid values | `g → 1` in the bulk |
| shear viscosity, `T* = 0.7205`, `ρ* = 0.8442` | `η* = 2.45 ± 0.41` | `η* = 3.258 ± 0.033` (Meier 2004 I) | mean of 3 seeds; `N = 256`, `rc = 2.5` |
| thermal conductivity, same state | `κ* = 6.48 ± 0.26` | `κ* ≈ 6.9` (argon-equivalent) | see below |
| self-diffusion, same state | `D_pbc* = 0.0248 ± 0.0046` (MSD), `0.0243` (VACF) | see below | two estimators agree to 2% |
| self-diffusion, Yeh–Hummer corrected | `D_0* = 0.0178` | Stokes–Einstein `0.024–0.035` | correction depends on the measured `η` |

Notes on the transport errors:

- **Viscosity.** The Green–Kubo integral of the stress autocorrelation is
  the hardest of the three because the LJ triple-point liquid has a long
  stress tail. A single 12 000-step, 256-atom trajectory gives
  `η*` anywhere in `2.0–4.3` depending on the seed; averaging three seeds
  gives `2.45 ± 0.41`. The seed-to-seed standard error is quoted. The
  reference `3.258` is the `N → ∞`, `rc → 5.5` limit; at `N = 256`,
  `rc = 2.5` Meier's own simulations give `η* ≈ 3.0–3.1`, so part of the
  25% gap is the systematic finite-size/cutoff effect, which we do not
  tune away.
- **Thermal conductivity.** The heat-flux autocorrelation decays on the
  scale of `τ ≈ 0.1` reduced time units, so `κ* = 6.48 ± 0.26` is much
  better determined. The argon-equivalent reduced value is
  `κ* = κ/(k_B √(ε/m)/σ²) = 0.12964/0.01880 ≈ 6.9` at 86.5 K, so the LJ
  model is about 6% low, within the LJ finite-size error.
- **Self-diffusion.** The MSD slope and the VACF integral agree to a few
  per cent, which is the internal consistency check. The Yeh–Hummer
  correction `D_0 = D_pbc − T ξ/(6πηL)`, `ξ = 2.837297`, subtracts a
  sizeable `0.006–0.007` at this box size; because it uses the measured
  (noisy) `η`, `D_0*` carries that uncertainty. The uncorrected
  `D_pbc* = 0.025` sits between the stick and slip Stokes–Einstein values
  computed from the Meier viscosity (`0.023` and `0.035`).

For the viscosity, a longer trajectory is the only fix; the
`examples/long_runs.rs` example runs 60 000+ steps so the plateau can be
inspected without slowing the test suite.

## The fold: closures for argon

`argon.rs` converts reduced transport coefficients to SI with the argon
constants and compares with measured saturated-liquid argon. The measured
values are those of B. A. Younglove and H. J. M. Hanley, *J. Phys. Chem.
Ref. Data* **15**, 1323 (1986) (viscosity and thermal conductivity), as
evaluated by CoolProp together with the argon reference equation of state
of Tegeler, Span & Wagner (1999). At the LJ triple-point temperature
`T* = 0.722` (86.50 K), saturated liquid argon has `ρ = 1400.3 kg/m³`,
`μ = 2.6676e−4 Pa s`, `k = 0.12964 W/(m K)`.

| closure | reduced input | SI output | measured Ar (86.5 K) | ratio |
|---|---|---|---|---|
| shear viscosity | `η* = 3.258` (Meier 2004 I) | `2.94e−4 Pa s` | `2.67e−4 Pa s` | `1.10` |
| thermal conductivity | `κ* = 6.48` (this work) | `0.122 W/(m K)` | `0.130 W/(m K)` | `0.94` |
| self-diffusion (PBC) | `D* = 0.0248` | `1.33e−9 m²/s` | — | — |
| self-diffusion (Yeh–Hummer) | `D_0* = 0.0178` | `9.6e−10 m²/s` | — | — |

The LJ viscosity is 10% high and the thermal conductivity 6% low relative
to real argon at the same reduced temperature. That is the expected size
of the corresponding-states error for the 12-6 model; the point of the
exercise is that the fold produces correctly dimensioned closures whose
errors are stated, not that the model replaces argon data.

## Emergent boundary: when does the continuum hold?

`validity.rs` returns a `ContinuumIndicator` from a self-diffusion
coefficient, a correlation time and the flow scales `(L, U, T, m)`. Two
independent numbers decide:

1. the **Knudsen number** `Kn = λ/L`, with the mean free path estimated
   from self-diffusion, `λ = 3D/v_mean`, `v_mean = √(8T/πm)`. The
   Navier–Stokes constitutive law with no-slip boundaries needs
   `Kn ≲ 0.1`;
2. the **scale separation** `L/(c_s τ)`, the ratio of the flow length to
   the correlation length `ℓ = c_s τ`, with `c_s = √(5T/3m)` the
   (ideal-gas) sound-speed estimate. Continuum behaviour needs
   `L ≫ c_s τ`; `τ` is the integrated stress or heat-flux correlation
   time, measured directly from the autocorrelation decay
   (`transport::correlation_time`).

For the LJ liquid above, `τ_stress = 0.139` reduced time units and
`c_s ≈ 1.1`, so `ℓ ≈ 0.15 σ`. A millimetre-scale flow (`L/σ ≈ 3e6`) has
`Kn ≈ 2e−8` and `L/(c_s τ) ≈ 2e7`: the continuum with these closures is
comfortably valid. A gas at the same temperature with `D* ≈ 50` has
`λ ≈ 111 σ`, so `Kn ≈ 0.1` already at `L/σ ≈ 1.1e3`, and the indicator
switches off. This is the
cheap check the ladder calls for: it uses only the already-computed
closure `D` and the already-computed correlation time, so it costs
nothing at run time and says when the atomic rung must be unfolded.

## Honest limitations

- The transport benchmark at the LJ triple point is short. The viscosity
  in particular has ~30% seed-to-seed scatter at 256 atoms and 10⁴ steps;
  the quoted error bar is the seed standard error, and the 25% offset from
  the infinite-size reference is partly physical finite-size/cutoff and
  partly statistics. Longer runs (`examples/long_runs.rs`) are needed for
  a 5% number.
- The finite-size correction to `D` uses the measured viscosity, so its
  uncertainty propagates.
- The barostat is a weak-coupling Berendsen scheme: correct averages, but
  wrong NPT fluctuations; production benchmarks are NVT.
- No long-range electrostatics, no bond constraints, no many-body
  potentials: the rung is pair potentials only, which is enough for the
  LJ closures and the argon comparison.

## Sources

- J. K. Johnson, W. Zollweg, K. E. Gubbins, *Mol. Phys.* **78**, 591
  (1993) — LJ equation of state (coefficients Table 10, pressure Eq. 7,
  internal energy Eq. 9).
- K. Meier, A. Laesecke, S. Kabelac, *J. Chem. Phys.* **121**, 3671
  (2004) — LJ viscosity; `η* = 3.258 ± 0.033` at `T* = 0.722`,
  `ρ* = 0.8442`.
- K. Meier, A. Laesecke, S. Kabelac, *J. Chem. Phys.* **121**, 9526
  (2004) — LJ self-diffusion (reference state).
- B. A. Younglove, H. J. M. Hanley, *J. Phys. Chem. Ref. Data* **15**,
  1323 (1986) — argon viscosity and thermal conductivity.
- R. Tegeler, R. Span, W. Wagner, *J. Phys. Chem. Ref. Data* **28**, 779
  (1999) — argon reference equation of state; via CoolProp.
- I.-C. Yeh, G. Hummer, *J. Phys. Chem. B* **108**, 15873 (2004) —
  finite-size correction `ξ = 2.837297`.
- G. Bussi, D. Donadio, M. Parrinello, *J. Chem. Phys.* **126**, 014101
  (2007) — stochastic velocity rescaling.
- H. J. C. Berendsen et al., *J. Chem. Phys.* **81**, 3684 (1984) —
  weak-coupling barostat.
- D. Frenkel, B. Smit, *Understanding Molecular Simulation* — LJ tail
  corrections.
- M. P. Allen, D. J. Tildesley, *Computer Simulation of Liquids* —
  velocity Verlet and Nosé–Hoover integration.
- LAMMPS `compute heat/flux` (S. J. Plimpton; Samolyuk & Pinto) — the
  heat-flux convention used here.

## Addendum: the triple-point viscosity miss is statistics-limited (GPU study)

The `η* = 2.45 ± 0.41` entry above was three 6 000-step, 256-atom seeds.
A GPU Green–Kubo study with per-step f64 reduction of the pressure-tensor
off-diagonals and seed/block error bars (see
[md_cuda.md](md_cuda.md#greenkubo-transport-on-the-device-size-and-cutoff-study))
reproduces that estimator and runs 40 000-step trajectories at
`N = 500, 2 048, 8 788, 32 000` and `rc = 2.5`, plus an `rc = 4.0` case.
At `N = 500`, `rc = 2.5` the plateau is `η* = 3.18 ± 0.09` (16 seeds),
within 3 % of Meier's `3.258`; the four sizes agree within the seed errors,
so the earlier 25 % gap was **statistics, not finite size or cutoff**. The
`rc = 4.0` running integral has a longer unresolved stress tail rather than
a shifted plateau. Self-diffusion with the CPU-compatible window gives
`D_pbc* ≈ 0.0316–0.0335` over the same sizes; the Yeh–Hummer correction is
the right size but does not flatten the residual `L`-dependence (the
measured `D_pbc` rises weakly with `L`, opposite to the correction's sign),
reported as a negative result.
