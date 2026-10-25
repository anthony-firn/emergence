# Compressible, multiphase, solid and per-element models: results

Crate: `fluid/crates/compressible`. Every number below comes from an
example in that crate and can be rerun.

## 1. Shocks and supersonic flow (`examples/validate`, `becker`)

Hybrid scheme: a 6th-order KEEP-PE central flux where the local
truncation estimate δ⁴(ln ρ, ln p) is below its tolerance, and WCNS
elsewhere (characteristic WENO5-Z interpolation, HLLC, conservative
6th-order midpoint-to-node difference). A Hu–Adams–Shu positivity limiter
runs in every stage.

| test | result |
|---|---|
| Sod, Lax, Toro-3 (10⁵ pressure ratio), 1 GPa water tube | ~first-order convergence (expected at discontinuities); Hybrid matches full-WCNS accuracy, ~2× faster, shock path on 1.4–2.5% of faces at N = 1600 |
| isentropic vortex (smooth) | central 6th order; WCNS more accurate than FD-WENO at every N |
| Becker viscous shock, Mach 2 | L1(ρ) 4.2e-3 → 5.7e-8 as Δx/δ goes 0.4 → 0.05; once the shock is resolved the sensor hands it to the central scheme (0% shock path) and the error is 12× below pure WENO |

## 2. Multiphase and cavitation (`multiphase_validate`, `bubble_collapse`, `cavitation_tube`)

Kapila five-equation model with a NASG equation of state per phase;
liquid–vapour phase change by relaxation to thermodynamic equilibrium.

| test | result |
|---|---|
| mixture sound speed | equals Wood's formula at every fraction (1625 m/s → 119 m/s at 1% air → 23.6 m/s at 50%) |
| water/air interface advected one period (833:1) | \|Δp\|/p ≤ 2e-10, \|Δu\|/u ≤ 1e-11 |
| 1 GPa water / 1 bar air tube vs exact | L1(ρ) 1.0e-3, L1(p) 3.8e-4 at N = 1600 |
| 3D Rayleigh collapse, p∞/p_g = 100 | follows Keller–Miksis, not Rayleigh–Plesset (R/R₀ 0.41 vs 0.44 at 1.3 t_c; RP gives 0.75); R_min 0.18 vs 0.10 at 12.8 cells/R₀, under-resolved at minimum |
| expansion tube, water at 355 K | without mass transfer 5.4 kPa (tension); with it p = p_sat(T_l) exactly, α_v 0.72, latent cooling 0.54 K |

Fixed on the way:
- clamping α broke pressure equilibrium;
- Kapila's K∇·u stiffness drove the vapour to p < 0, so there is now a
  dilatation time-step limit.

Known limit: the published NASG water fit is +20% in p_sat at 300–373 K
(+2% at 473 K).

## 3. Solids (`gpr_validate`)

GPR unified continuum model (distortion tensor A, strain-relaxation time
τ). Relaxation is exponential on the singular values of A and applied with
Strang splitting.

| test | result |
|---|---|
| elastic shear wave | speed = c_s to 1e-8 |
| elastic longitudinal wave | √(c₀² + 4c_s²/3) in the linear limit (shift ∝ amplitude) |
| Stokes' first problem | GPR → Navier–Stokes as De → 0: 4.2e-5 at De = 0.06; second-order convergence at De = 9e-4, limited by the (Δt/τ)² splitting error |
| elastic–plastic impact | precursor at c_L − u; plastic shock at its Rankine–Hugoniot speed; plateau stress equals the elastic relation to 0.03–0.07%. The power law yields at 0.93 σ_Y (n = 400), so the precursor is 7% below ideal plasticity |

## 4. Periodic table (`periodic_table`, `periodic_flow`)

118-element dataset with per-value sources (`fluid/data/build_elements.py`).
Each element with data gets solid, liquid and vapour NASG phases fitted at
1 atm.

Grade P (11 elements: Pa, Bk, Cf, Es, Rf, Cn, Nh, Mc, Lv, Ts, Og): no
measured melting or boiling point exists, so the missing inputs come from
published predictions and estimates (the Wikipedia element infoboxes,
`fluid/data/build_predicted.py` → `elements_predicted.csv`). Measured
values are never overwritten. 104 of 118 elements are now modelled
(B 65, C 28, P 11; 14 not modelled).

| check (never fitted) | result |
|---|---|
| standard state at 298 K | 98/98 correct |
| triple points (13 found) | median \|ΔT\|/T 0.02%, pressure within a factor 1.06 |
| vapour curves 1 Pa – 100 kPa | median rms \|ΔT\|/T 1.03% (63 grade-B elements) |
| equilibrium round trips (melting, boiling) | 203/204 to ~1e-15 |
| cavitation tube per element in the flow solver | 101/102 settle on p_sat(T_l) |

Source-data errors found by the sweep:
- potassium's latent heats are swapped in mendeleev;
- the tin "melting point" is its allotropic transition;
- Se, Te and Cd vaporization heats are 35–55% low;
- fluorine's triple-point pressure is wrong in the dataset.

Not modelled or limited:
- 14 elements without melting/boiling data or predictions (superheavy:
  Fm–Lr, Db–Rg, Fl);
- He, a quantum liquid;
- S, whose vapour speciation needs a chemistry model;
- critical points, which need a reference equation of state.

## 5. Adaptive multiresolution grids (`amr.rs`, `examples/amr_validate`)

The adaptive grid is a block-structured octree of 16^D blocks. Refinement
is driven by interpolating-wavelet details on the sensor values: refine
above ε, coarsen below ε/128. Ancestor blocks are kept and refreshed by
6th-order restriction; ghost sampling recurses only to coarser levels;
coarse faces are refluxed from the fine fluxes.

| test (adaptive vs uniform finest grid) | result |
|---|---|
| fully refined | identical to 2e-13 |
| Shu–Osher, WCNS mode, ε = 1e-3 … 1e-5 | 2.9e-5 → 1.2e-6 → 1.6e-7 (faster than linear in ε), vs discretization error 4.4e-3 |
| Shu–Osher, hybrid mode | floor ~1.2e-4: the hybrid switch's sensitivity in the under-resolved wave train (2.6% of the discretization error) |
| 2D four-shock Riemann, ε = 1e-2 … 1e-4 | 2.9e-5 → 1.5e-6 → 7.9e-7, vs 2.8e-3 |

Performance: in the 2D test 65% of cells must be fine, so there is no gain
(0.7× after removing overheads: dense block index, edge-only restriction).
3D bubble collapse at full refinement costs the same as the uniform grid
(497 vs 498 s). Effective-128³ run with 26% of the cells: see the end of
this file.

## 6. GPU (`compressible-cuda`)

Kernels: per-face blend weight θ; midpoint WENO5-Z/HLLC fluxes computed
once and only where needed; face combine with the positivity limiter;
SSP-RK3.

- GPU f64 equals the CPU to round-off: 2.6e-13 on the 2D Riemann problem
  at 128², 2.7e-13 on the water/air tube.
- 512², Performance profile, CPU shared with a background run:
  - GPU f64: 3.3× the CPU;
  - GPU f32: 11×, at 19.6 M cell-steps/s at 256².
- Under the system76 Battery profile the GPU ran at 210 MHz and gave
  1.2 M cell-steps/s.
- f32 with water's p∞ = 6e8 has 3e-3 relative pressure error.

## 7. Solid + liquid + vapour in one flow (`threephase.rs`)

Three-phase Kapila model, plus maximum-entropy phase equilibrium over
single phases, phase pairs and the triple point. Nitrogen:

- triple point 63.144 K / 13.3 kPa (data 63.151 K / 12.5 kPa);
- eight equilibrium round trips to ≤ 1e-12;
- triple-point tube:
  - expansion freezes and evaporates the liquid (solid 0.30 → 0.84, then
    onto the sublimation line: the slush mechanism);
  - compression melts and condenses.

The solid has no shear strength yet.

## 8. Reference equations of state (`helmholtz.rs`, `realfluid.rs`, `realflow.rs`)

- **IAPWS-95 water and Span et al. (2000) nitrogen.** Exact derivatives by
  automatic differentiation; coefficients from CoolProp's files.
  - IAPWS-95 verification Tables 7 and 8 pass to < 1e-8.
  - vs CoolProp: ~1e-7 (liquid-N₂ p 6e-6 = a 2e-8 density offset).
- **Tabulated equilibrium EOS for the flow solver.** Single-phase plus a
  homogeneous-equilibrium dome, bicubic in (φ(ρ), e), a dome-boundary
  table, and exact fallback where the table cannot vouch for a value.
  - Water, 2 × 200²: max error p 2.5e-4, T 4e-5, c 1.9e-3; 105 ns per
    lookup (exact ~100 µs).
  - Nitrogen: same quality except ~10 of ~7000 samples within 1 K of the
    critical point.
- **Cavitation on the reference EOS.** Pure liquid water at 355 K:
  p_centre / p_sat,IAPWS(T) = 1.000031. The NASG model was ~20% high.

Fixed on the way:
- molar ancillary densities;
- an energy range polluted by spinodal states;
- the dome cap near T_c;
- unphysical Newton roots inside the dome (stable-root scan);
- out-of-table states counted as admissible.

### 3D Rayleigh collapse on the adaptive grid (effective 128³, ε = 1e-2)

Output: `analysis/bubble_collapse_amr128.out`.

| | uniform 32³ | uniform 64³ | adaptive, effective 128³ | Keller–Miksis |
|---|---|---|---|---|
| R_min / R₀ | 0.224 | 0.177 | **0.167** | 0.101 |
| t(R_min) / t_c | 1.10 | 1.10 | 1.10 | 1.05 |
| cells (% of uniform 128³) | — | — | **18%** | — |

- Early collapse (t ≤ 0.4 t_c) tracks Keller–Miksis to 0.1%.
- Later the simulation lags. The 64³ and 128³ results agree closely
  through 0.9 t_c (0.583 vs 0.591), so that part of the gap is the setup
  (5 R₀ box, outflow boundaries), not resolution.
- R_min is still converging (0.224 → 0.177 → 0.167).

Wall time: 3221 s, against an estimated ~7950 s for uniform 128³ (2.5×).
The profile shows the overhead is the remaining target:

| phase | seconds |
|---|---|
| ghost filling | 1112 |
| ancestor restriction | 580 |
| fluxes | 590 |
| regrid | 339 |

With ε = 1e-2 the hierarchy never changed after initialization (60
coarse + 32 fine blocks), so here every regrid cost time without
changing the grid. Next steps:
- block-to-block ghost copies instead of per-cell sampling;
- regridding triggered by detail growth rather than a fixed interval;
- subcycling of coarse levels (multirate).

## 9. Open boundaries and surface tension (`boundary.rs`, `capillary.rs`)

**Non-reflecting boundaries.** Ghost states come from the characteristic
decomposition at the boundary cell. Outgoing wave amplitudes are
extrapolated linearly from the interior; incoming amplitudes come from the
far field.

| residual after a pulse has left (`boundary_validate`) | zero-gradient | characteristic |
|---|---|---|
| 1D air pulse | 1.01 | 1.6e-3 (constant-amplitude ghosts: 3.4e-2) |
| 2D cylindrical pulse, all incidence angles | 0.226 | 7.5e-3 |
| water pulse, five-equation model | 1.01 | 1.6e-3 |

With the open far field, the 3D bubble collapse reaches R_min at the
Keller–Miksis time (1.05 t_c; zero-gradient boundaries lagged to 1.10).
R(0.8 t_c) = 0.655, against 0.673 for Keller–Miksis and 0.706 with
zero-gradient boundaries.

**Surface tension** uses the continuum surface force σκ∇α (plus its work
in the energy equation), with the curvature computed from a lightly
smoothed volume fraction.
- Laplace law, 2D droplet: Δp within 0.2–3% of σ/R.
- Planar capillary wave at N = 64, low Mach, frequency error by density
  ratio:

  | density ratio | 1 | 10 | 100 | 857 |
  |---|---|---|---|---|
  | ω error | −4.8% | −2.5% | +0.45% | −0.9% |

  At 128² the error is −2.5% (ratio 1) and −1.2% (ratio 10), so it
  converges with resolution. The earlier incompressible Cahn–Hilliard
  solver was 20–30% off whenever the density ratio was not 1.
- Weakness: parasitic currents reach about 15% of √(σ/ρR). A
  balanced-force discretization would fix this.

## 10. Adaptive grids: overheads and local time stepping (`amr.rs`)

- **Ghost exchange.** Ghost cells are filled block-to-block:
  - neighbouring blocks at the same level are slab copies;
  - coarser regions get one gathered coarse patch plus separable
    6th-order prediction;
  - physical boundaries are applied inside the gathered box, axis by axis.
- **Restriction** is separable, restricts only the needed edge slabs, and
  shifts its stencils per line.
- **Regrid interval** doubles while the grid doesn't change, capped at 2×.

  Overhead on the 2D Riemann problem fell from 290% to 27% of the flux
  work, so the AMR runs are 1.1–1.2× faster than uniform even with 65% of
  cells refined.
- **Characteristic boundaries in the AMR** first cost 7×, because
  boundary cells were re-sampled through the hierarchy. Filling them
  inside the gathered box fixed that.
  - 3D bubble, effective 128³, to 0.2 t_c: **954 s → 308 s**, same R.
  - Profile now: fluxes 186 s, ghosts 40 s, restriction 16 s, regrid 7 s.
- **Subcycling** (Berger–Oliger/Berger–Colella): each level steps with
  Δt/2ˡ, coarse ghosts are interpolated linearly in time, and flux
  registers keep conservation exact.

  | `amr_validate` (ε = 1e-3), error vs uniform grid | global step | subcycled |
  |---|---|---|
  | 2D isentropic vortex | 3.1e-5 (6.8× faster than uniform) | 2.8e-5 (7.7×) |
  | Shu–Osher | 1.6e-4 | 3.0e-4 |
  | 2D Riemann | 1.4e-6 | 4.6e-5 |

  On the 2D Riemann problem, linear-in-time coarse ghosts cost accuracy at
  the shocks, but the error stays 60× below the discretization error.
  Subcycling pays off when the refined region is small. On the 3D bubble
  (26% of cells fine, a regrid every coarse step) it did not finish
  0.2 t_c within 3000 s on a shared CPU, against 308 s for the global step.
- Shu–Osher (WCNS) currently gives 3.8e-5 / 6.6e-6 / 4.4e-7 at
  ε = 1e-3 / 1e-4 / 1e-5. The adaptive regrid interval raised the middle
  value from the 1.2e-6 in section 5. It still decreases with ε and stays
  far below the discretization error.

## 11. Adaptive grids on the GPU (`compressible-cuda/src/amr.rs`)

The CPU keeps the hierarchy (refine and coarsen decisions). The solution
lives on the GPU, with all blocks (leaves and ancestors) in one pool with
face ghost slabs. Each stage runs batched kernels:
1. restrict the edge ancestors;
2. fill ghost slabs: same-level copy, separable 6th-order prediction in
   coarser regions, boundary conditions including nested characteristic
   corners;
3. per level and axis: sensor, Riemann fluxes, face fluxes, and storage of
   the boundary-face fluxes;
4. reflux (per axis and side, in the CPU's order);
5. SSP-RK3 update.

The wavelet details are computed on the GPU too; the CPU only edits the
tree.

**Accuracy, GPU f64 vs the CPU AMR:**
- 1.5e-14 on the 2D Riemann problem, regrids included;
- 3e-14 after 8 steps, growing smoothly (round-off) to 7e-11 after 1157
  steps on a 2D bubble collapse with characteristic far field.

The uniform GPU solver still matches the CPU to 1.5e-13 and now supports
characteristic boundaries.

**Speed**, 3D bubble, effective 128³, 26% of cells, RTX 3050 laptop,
f32: to 0.2 t_c in 579 steps, **62 s against 308 s** for the CPU AMR on 11
threads (5.0×), or 0.11 s per step against 0.53 s. The f32 result differs
from f64 by 1.2e-4 (relative, volume fraction).

| phase (58 steps, f32) | first version | final |
|---|---|---|
| ghost filling (direct 6³ prediction → separable) | 12.0 s | 1.3 s |
| regrid (CPU details → GPU details; separable parallel refine) | 14.9 s | 1.2 s |
| fluxes | 3.2 s | 3.1 s |
| total | 31.7 s | 7.3 s |

## 12. Incompressible ↔ compressible by local Mach number (`allspeed.rs`)

The obvious version of this idea — couple the periodic spectral
incompressible solver to the compressible one — fails on two counts: a
periodic incompressible field cannot absorb a net volume change (a
collapsing bubble is one), and every interface between two solvers leaks
mass and energy. So the switch is built into one conservative scheme,
after Kwatra et al. (2009):
1. **Explicit advection without pressure.** SSP-RK3 with WENO5-Z, and a
   numerical dissipation that scales with |u|, not |u| + c (the acoustic
   dissipation is what ruins compressible schemes at low Mach).
2. **One implicit pressure solve per step.** Based on Dp/Dt = −ρc²∇·u, the
   system is symmetric positive definite and solved by CG for the pressure
   increment, with its right-hand side formed in difference form.
3. **Conservative corrections**, plus re-transport with the projected face
   velocity so that every dilatation term sees the same velocity.

The time step is set by the flow speed, not the sound speed. The local
acoustic Courant number C = cΔt/Δx = CFL/M sets the regime, face by face
and cell by cell:

| | C ≤ 1 (sound resolved) | C ≫ 1 (sound unresolved) |
|---|---|---|
| time weighting θ | ½ (Crank–Nicolson, undamped sound) | 1 − 1/(2C) → 1 (backward Euler; projection enforces ∇·u → 0) |
| energy | conservative | elliptic pressure from C ≥ 64 (fully at 128) |
| diffuse multiphase interfaces | elliptic pressure | elliptic pressure |
| whole step | explicit when the all-speed Δt is < 3× the explicit Δt, or when an all-speed step leaves a state the explicit scheme cannot continue from | all-speed |

| test (`allspeed_validate`, `bubble_collapse ALLSPEED=inf`) | all-speed | explicit |
|---|---|---|
| Gresho vortex 64², Mach 0.1 | 215 steps, velocity error 6.5e-3 | 3430 steps, 6.8e-3 |
| Gresho, Mach 0.01 | 214 steps, 8.7e-3, 3.6 s | 32223 steps, 5.1e-3, 72 s |
| Gresho, Mach 0.001 | 214 steps, 8.9e-3 | ~3×10⁵ steps (not run) |
| Sod shock tube | identical to explicit (switch picks explicit) | 2.4e-3 |
| acoustic pulse, C = 0.5, all-speed forced | 2.0e-2 of amplitude | 1.1e-3 |
| water slab advected one period at C = 64 | p uniform to 3e-11 | — |
| 3D bubble collapse, 32³, R(0.95 t_c) | 0.439 (KM 0.456), 430 steps | 0.404, 971 steps |

- The step count is independent of Mach number, at the explicit scheme's
  accuracy.
- Water/air shock tube at 1e7 | 1e5 Pa: runs, with pressure error
  6.5e4 vs 4.8e4 explicit.
- At 1e9 | 1e5 Pa it is outside the scheme's range; use the explicit
  scheme for strong water/air shocks.

**Instabilities found and fixed on the way.**
- *Collocated pressure–velocity mismatch.* Fixed by keeping the pressure
  out of the predicted velocity, which is otherwise averaged.
- *Forward-Euler internal-energy compression amplified sound.* Fixed by
  re-transporting with u_θ.
- *Kapila/Wood interface stiffness.* In a stiffened-gas mixture the
  stiffness is the small difference of two O(γp∞) terms, which a
  conservative energy update cannot track once C ≫ 1. Fixed by taking the
  elliptic pressure in interface cells.
- *Round-off amplified by C².* Fixed with the difference-form residual
  plus the elliptic pressure in pure cells at C ≥ 64.
- *Characteristic ghosts during advection.* They turned meaningless stage
  pressures into velocities. Advection stages now use zero-gradient
  ghosts, and the acoustic condition becomes an impedance law in the
  pressure equation.

## 13. Transported scalars, finite-rate cavitation, gas dissolution (`scalars.rs`, `kinetics.rs`)

**Scalars** are carried by the mixture or by one phase; dissolved nitrogen
rides on the liquid's mass flux. Each step uses the carrier's effective
face fluxes from that step (RK3-weighted), so a uniform φ stays exactly
uniform. TVD face values keep φ bounded. Diffusion acts within the
carrier, and a per-cell kinetics hook can change the state and its scalars
together.

| test (`scalars_validate`) | result |
|---|---|
| Gaussian, one period | L1 3.6e-3, bounded, total conserved to 3e-16 |
| liquid-carried φ = 1 through a moving water/air interface | 1 to 6e-13 (explicit), 9e-15 (all-speed, C = 64) |
| diffusion of a Gaussian | 8e-4 from the exact solution |

**Cavitation with nucleation.**
- Activation:
  - heterogeneous nuclei (gas pockets) at the Blake threshold;
  - homogeneous nucleation from classical nucleation theory.
- Growth and collapse at the Schnerr–Sauer rate (nuclei density,
  Rayleigh–Plesset velocity), never beyond the thermodynamic equilibrium
  state.

| check (`kinetics_validate cnt`) | value |
|---|---|
| CNT cavitation pressure of water at 300 K (J·1 mm³·1 s = 1) | −143 MPa (water in quartz inclusions: about −140 MPa, Azouzi et al. 2013) |
| Blake threshold, 1 / 10 / 100 µm nucleus at 1 bar | p_v − 42.7 kPa / − 2.0 kPa / − 67 Pa |

Expansion tube from pure liquid at 355 K (`kinetics_validate tube`):
liquid pulled apart at ±2 m/s, 200 cells, 10 µm nuclei.

| model | centre p | α_v | lowest p during the run |
|---|---|---|---|
| no phase change | −4.46 MPa (tension) | 0 | −4.50 MPa |
| equilibrium (HEM) | 60 325 Pa (= p_sat) | 0.5156 | 60 325 Pa |
| finite rate, n = 10⁸ m⁻³ | 60 390 Pa | 0.4809 | 3.0 kPa |
| finite rate, n = 10¹⁰ m⁻³ | 60 327 Pa | 0.5147 | 3.0 kPa |
| finite rate, n = 10¹² m⁻³ | 60 327 Pa | 0.5144 | 3.0 kPa |
| finite rate, n = 10¹⁴ m⁻³ | 60 326 Pa | 0.5150 | 58.4 kPa |

- The finite-rate model converges to the equilibrium model as the nuclei
  density grows: for n ≥ 10¹⁰ m⁻³ the vapour fraction agrees to 0.2%.
- With few nuclei, the pressure dips below saturation before the vapour
  catches up.
- One limit: a single-pressure model cannot hold vapour beside liquid in
  tension (the bubble interiors sit at p_v). So once vapour has
  nucleated, the transfer is at least enough to lift the mixture pressure
  to 5% of p_v; that 3.0 kPa floor shows up in the rows above. Holding
  tension after nucleation would need pressure non-equilibrium, for
  example the 6-equation model or Rayleigh–Plesset sub-grid bubbles.

**Gas dissolution.**
- Saturation from Henry's law: c_s = k_H p (N₂: 17.9 mg/L at 1 bar).
- Transfer through the sub-grid bubble area, with k_L = Sh·D/2R_b,
  integrated exactly over each step.
- In the all-speed mode the transfer acts before the pressure solve and
  enters it as a volume source, so a dissolving bubble draws liquid in
  rather than stretching it.

N₂ bubbles (3%, 10⁹ m⁻³) dissolving in degassed water, open ends at 1 bar
(`kinetics_validate henry`): 1000 s in 500 steps of 2 s. That is an
acoustic Courant number of 6.5×10⁶; an explicit run would need about
6.5×10⁹ steps.

| t (s) | 10 | 30 | 100 | 300 | 1000 |
|---|---|---|---|---|---|
| centre c (mg/L) | 0.870 | 2.469 | 6.870 | 13.43 | 17.66 |
| independent ODE of the transfer law | 0.868 | 2.461 | 6.821 | 13.31 | 17.62 |
| gas fraction | 0.0294 | 0.0284 | 0.0254 | 0.0206 | 0.0173 |

N₂ (gas + dissolved) is conserved to 5e-9. The dissolved gas's own volume
is neglected (2e-5 of the water volume at saturation).

**Still open:**
- water vapour inside the nitrogen bubbles (a third phase or a
  gas-mixture EOS);
- exsolution without existing bubbles (a degassing nucleation threshold);
- scalars on the adaptive grid and the GPU.

## 14. Molecular substances (`substances.rs`, `examples/substances`)

Water, CO₂, ammonia and methane as solid + liquid + vapour models. Each is
fitted with the same routine as the elements (`elements::fit_inputs`) to
1-atm data only: the transition temperatures, latent heats, densities,
heat capacities and sound speeds. Inputs that are estimates are listed per
substance (`Substance::inputs.estimated`). Everything in the table is
outside the fit.

| check (never fitted) | water | CO₂ | NH₃ | CH₄ |
|---|---|---|---|---|
| boiling curve, \|ΔT\|/T | ≤ 0.29% (3–476 kPa) | 1.3–8.2% (0.9–5.7 MPa, toward T_c) | 0.6–2.0% | 0.3–7.2% |
| sublimation curve | ≤ 0.30% (13–260 Pa) | — | — | — |
| triple point | 273.16 K, 642 Pa (611.7) | 216.57 K, 485 kPa (518.5) | 195.42 K, 6.25 kPa (6.06) | 90.67 K, 12.4 kPa (11.7) |

Water's ice melting curve (IAPWS R14-08):

| p (MPa) | 15 | 39 | 92 | 138 | 179 |
|---|---|---|---|---|---|
| model T (K) | 272.00 | 270.02 | 265.32 | 260.97 | 256.91 |
| IAPWS (K) | 272.0 | 270.0 | 265.0 | 260.0 | 255.0 |

- Melting-point slope at 1 atm: −0.0744 K/MPa (measured −0.074). It falls
  out of the two densities, because ice is lighter than water.
- CO₂ sublimes at 1 atm (194.7 K) and has its triple point above 1 atm.
  Both come out right.
- The errors grow toward the critical point. The NASG vapour is a
  stiffened ideal gas, and near T_c only a reference EOS (§8) is accurate.
- Water's density maximum at 4 °C is not representable (constant
  expansion coefficient).

## 15. State changes driven by heat (`isobaric.rs`, `allotropes.rs`, `examples/phase_change`)

Freezing a lake, melting iron or condensing nitrogen on a cold wall is
slow. The fronts move at mm/s, while sound crosses a cell in nanoseconds.

**Why not the compressible solver.**
- The compressible solver can conduct heat (`thermal.rs`: explicit
  sub-steps, phase relaxation after each, harmonic face conductivity).
  The all-speed mode can also take steps at acoustic Courant numbers near
  10⁵. To make that work, it now carries a conservative energy plus a
  solved pressure in pure cells, and resets the energy only at interfaces
  with stiffness contrast above 100×.
- But in stiff matter a relative density error δρ/ρ gives a temperature
  error of about Γ·T·δρ/ρ. At these Courant numbers that is kelvins: more
  than the whole driving temperature difference of a Stefan problem.
- The low-Mach limit is the right model: pressure fixed in time, and
  enthalpy carried by conduction. That is the isobaric model below. It
  is the emergent fold of the full equations for this regime, not a
  separate physics.

**The isobaric model.**
- One-dimensional, in mass coordinates (cells keep their mass, and their
  widths follow the density of whatever phase they hold).
- State: enthalpy plus phase mass fractions at fixed p.
- The phase ladder (which phase is stable at each T) is found from the
  Gibbs energies of the fitted phases. Up to 6 phases, so allotropes are
  allowed.
- A cell at a transition holds the two phases as layers. The interface is
  its temperature point, and each side conducts through its own phase only:
  a condensing film conducts through the liquid, not the saturated vapour.
  Without this layering, N₂ condensation ran 86% slow.
- Explicit steps, with the stability limit set by each cell's heat
  capacity against its face conductances.

**Stefan problems vs the exact Neumann solution.** The new phase grows
from a wall held ΔT past the transition, and is compared with
s = 2λ√(αt), λe^{λ²}erf λ = St/√π. This solution holds even with a density
change, because the new phase sits at rest against the wall. 400 cells,
t = 1 s.

| case | transition | St | front error at t = 0.2 / 1 s | T halfway into the new phase |
|---|---|---|---|---|
| model material (equal densities) | melting 1811 K, wall +189 K | 0.612 | +0.00% / +0.00% | 1898.295 K (exact 1898.295) |
| iron (periodic-table fit) | melting 1811 K, wall +189 K | 0.395 | +0.31% / +0.31% | 1901.21 K (exact 1900.90) |
| lake ice (water at 0 °C) | freezing 273.15 K, wall −10 K | 0.063 | −0.04% / −0.04% | 268.247 K (exact 268.245) |
| liquid N₂ on a cold wall | condensation 77.35 K, wall −10 K | 0.105 | −1.92% / −1.34% | 72.52 K (exact 72.45) |
| silicon (periodic-table fit) | melting 1687 K, wall +100 K | 0.046 | +0.13% / +0.13% | 1736.86 K (exact 1736.80) |

- Liquid silicon is denser than the solid (2570 vs 2305 kg/m³), so its
  melting point falls with pressure, like ice. From the two fits:
  −0.042 K/MPa at 1 atm. The average slope from 1 atm to the
  Si-I/Si-II/liquid triple point (≈ 10.5 GPa, ≈ 1000 K) is about
  −0.065 K/MPa. Silicon as a semiconductor is in `semiconductor.md`.
- The offsets come from thermal expansion of the new phase, which the
  constant-property Neumann solution leaves out. The model material has
  negligible expansion and matches to 0.00%.
- Measuring the front by mass instead of thickness (each cell's new-phase
  mass at the transition density) gives:

  | case | front by thickness | front by mass |
  |---|---|---|
  | iron | +0.31% | −0.29% |
  | lake ice | −0.04% | +0.04% |
  | liquid N₂ | −1.34% | +1.04% |

- Iron: the hot liquid is about 0.6% thicker than at 1811 K, and its lower
  density lengthens the conduction path.
- Liquid N₂ expands by about 0.5%/K. Its early-time offset also includes
  the first cell's layering while the film is thinner than a cell.

**Supercooled droplets: water to ice by nucleation.**
- Finite-rate freezing: nuclei accumulate at the homogeneous rate J(T),
  then the ice grows by Avrami kinetics up to equilibrium.
- J(T) is from Koop et al. (2000), log₁₀J = −906.7 + 8502Δa − 26924Δa² +
  29180Δa³ (cm⁻³s⁻¹), with Δa = 1 − p_ice/p_liq.
- Both vapour pressures come from the fitted water model, not from
  freezing data. At 238 K the model's p_ice/p_liq is 0.6955, the same as
  Murphy & Koop (2005).
- Droplets cooled at 1 K/min:

| diameter | freezes at | same law, Murphy–Koop vapour pressures | ice fraction after recalescence (energy balance) |
|---|---|---|---|
| 1 µm | 235.28 K | 234.73 K | 0.4787 (0.4787) |
| 10 µm | 237.21 K | 236.78 K | 0.4544 (0.4544) |
| 100 µm | 238.91 K | 238.58 K | 0.4328 (0.4328) |

- Measured homogeneous freezing of pure water droplets at about 1 K/min is
  235–237 K.
- On nucleation, the latent heat warms the droplet back to 273.15 K
  (recalescence). The ice formed in that flash matches the energy
  balance to every printed digit.
- Ice-crystal shapes (dendrites, habits) are not resolved; Avrami growth
  is a sub-grid model.

**Allotropes: iron and tin** (`phase_change allotropes`).

Iron at 1 atm goes bcc (α) → fcc (γ) at 1185 K → bcc (δ) at 1667 K →
liquid at 1811 K.
- α and δ are one structure, so a single bcc phase must win both below
  and above fcc.
- Both solid phases now carry the temperature-dependent thermal EOS of
  Dorogokupets, Dymshits, Litasov & Sokolova, Sci. Rep. 7, 41863 (2017),
  Table 1 (`thermal_eos.rs`): a Vinet cold curve, an Einstein thermal term
  with an Altshuler Grüneisen parameter and a free-electron term. One
  `ThermalEos` is held per solid phase in `PhaseSet`; the NASG fields are
  kept as a fallback for code that reads `eos` directly. The phase-set
  methods (`gibbs`, `enthalpy`, `cp`, `specific_volume`) dispatch to the
  thermal EOS when it is present, so every NASG phase is bit-identical.
- Volume, bulk modulus, thermal expansion, the temperature dependence of
  c_p and the electronic term are taken from the published EOS. The bcc
  phase also carries the Inden / Hillert–Jarl magnetic Gibbs energy
  (CALPHAD; Dinsdale 1991), T_C = 1043 K, bcc p = 0.4, β = 2.22 μ_B. (The
  paper's magnetic convention differs from Dinsdale's only by a term linear
  in T, which the reference gauge absorbs; the c_p is the paper's λ-shaped
  one either way.)
- Only the reference of the two phases is fitted, exactly as the NASG
  model fitted bcc c_v and q, q′: a gauge q − q′T plus a constant c_v
  offset on bcc (the analogue of the single fitted c_v) is chosen so that
  g_bcc = g_fcc at 1185 K and at 1667 K and ΔH(α→γ) = 0.90 kJ/mol. The
  offset is +2.345 J/(mol K); everything else is the published EOS. The
  volume jumps and ΔH(γ→δ) are predictions.
- The published EOS is not extrapolated above the melting point at
  0.1 MPa: the paper says it is unreliable there, and the bcc isotherm
  indeed turns over (K_T < 0 above ≈ 2000 K), which would spuriously
  re-stabilise bcc near 2890 K. Above 1811 K the solid phases are omitted
  from the stability search.

| iron | NASG (previous) | thermal EOS (this work) | measured |
|---|---|---|---|
| ΔV/V at α→γ | −0.79% | −0.77% | ≈ −0.98% (lattice parameters) |
| ΔH(γ→δ) | +0.630 kJ/mol | +0.554 kJ/mol | 0.84 kJ/mol |
| ΔV/V at γ→δ | +0.06% | −0.03% | ≈ +0.53% |
| ΔV/V on melting (from δ) | +4.94% | +4.3% | — |

- The α→γ contraction is essentially unchanged (−0.77% vs −0.79%): at
  1185 K the published EOS gives the same measured bcc and fcc volumes
  the NASG fit had used (7.377 and 7.320 cm³/mol).
- The γ→δ volume **does not improve**: −0.03%, against the measured
  +0.53%. The earlier attribution to a single NASG expansion coefficient
  is refuted — the full temperature-dependent bcc and fcc EOSs themselves
  put V_δ ≈ V_γ at 1667 K. The fcc EOS expands faster (α_v = 70.7e-6 /K
  at 1185 K, 83.1e-6 at 1667 K) and catches the bcc volume.
- ΔH(γ→δ) is +0.554 kJ/mol, below the measured 0.84 (the NASG fit gave
  +0.630). The α→γ latent heat is now carried mostly by the fitted c_v
  offset rather than by magnetism.
- What does improve is the temperature dependence: α(T) and c_p(T) are
  now the published functions, not constants. bcc α_v = 35.5, 53.4 and
  67.4e-6 /K at 300, 1185 and 1667 K, fcc α_v = 53.3, 70.7 and
  83.1e-6 /K (Basinski et al. 1955 for bcc, Kohlhaas et al. 1967 for
  fcc). The published bcc EOS with its magnetic term (before the fitted
  c_v offset) gives c_p = 25.3 J/(mol K) at 300 K (measured 25.1, Desai
  1986, +0.8%), 54.0 at 1000 K, 60.1 at the Curie peak, 41.5 at 1185 K
  and 41.1 at 1667 K (measured δ ≈ 41, +0.2%); fcc c_p = 36.2 at 1185 K
  (measured γ ≈ 34, +6%) and 42.8 at 1667 K. The fitted offset shifts
  every bcc c_p up by 2.345 J/(mol K), so the model gives 27.6 at 300 K
  (≈ 10% high) and 43.8 / 43.5 at 1185 / 1667 K (≈ 6% high), and the
  sound speed at 1185 K drops from the published 4314 to 4305 m/s
  (−0.2%). That is the price of reproducing ΔH(α→γ) = 0.90 kJ/mol with
  this EOS.

Heating one cell at 1 kW/kg gives the classical thermal-arrest curve (the
way these transitions were discovered). The temperature halts at each
transition for latent heat ÷ power:

| transition | 1185 K | 1667 K | 1811 K | 3134 K |
|---|---|---|---|---|
| arrest | 16.1 s | 9.9 s | 247.3 s | 6088 s |

Cooling gives the same arrests in reverse.

Tin: grey (α, diamond cubic) → white (β) at 286.35 K (2.09 kJ/mol,
ΔV/V −20.6%), then melting at 505.08 K.
- The tin densities are inputs.
- In reality white → grey on cooling ("tin pest") takes nucleation and
  can be delayed for years. The model relaxes to equilibrium instantly;
  only liquid → solid has kinetics so far.

**Still open:**
- the isobaric model in 2D/3D, and coupled to flow (convection in melt
  pools and lakes);
- temperature-dependent heat capacities beyond the iron magnetic term
  (water's density maximum);
- solid–solid nucleation kinetics (tin pest, martensite);
- dendritic ice and solidification microstructure.
