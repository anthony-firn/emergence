# The emergence ladder: one pattern behind every speed-up

Every accurate shortcut in this repository has the same shape:

1. a **coarse model** (cheap, valid in a regime);
2. a **validity indicator**: a local, cheap quantity that bounds or
   estimates the coarse model's error;
3. where the indicator fails, either **unfold** to a finer model locally,
   or replace the unresolved structure by a **closure**: a law the finer
   model obeys, computed once, which the coarse model can use (a jump
   condition, an equation of state, a rate).

The speed comes from step 2: the expensive model runs only where the
cheap one is certified to fail.

This pattern is now executable code shared by the crates: the
`emergence` engine (typed error claims, cost-optimal rung choice, audits,
validated closure tables and goal-oriented error budgets). Its design and
every measured result are in [emergence_engine.md](emergence_engine.md);
the Mie, grey-phonon-slab and certified spectral-resolution ladders below run
on it.

## Where the pattern already appears

| coarse model | indicator | unfold / closure | where |
|---|---|---|---|
| 6th-order central (non-dissipative) | local truncation estimate δ⁴ ln q | WCNS + HLLC Riemann solver; Rankine–Hugoniot is the closure for the sub-grid shock | `compressible/scheme.rs` |
| resolved viscous shock | same sensor turns *off* once the shock spans cells | Navier–Stokes interior (Becker profile reproduced to 5.7e-8) | `examples/becker.rs` |
| incompressible NS | Mach number, λ\|∇u\|/c, λ/Δx | compressible solver; kinetic theory beyond Kn ~ 0.01 | `regime.rs` |
| low-Rm quasi-static MHD | Rm = μ₀σUL | induction equation (not built) | `mhd.rs` |
| Joule damping (stiff when the Hartmann number is high) | exact for a uniform field | integrated exactly in the integrating factor | `mhd.rs` |
| mechanical-equilibrium mixture (Kapila) | — | Wood sound speed: the closure of stiff pressure relaxation | `multiphase.rs` |
| no phase change | p < p_sat(T) − tension | thermodynamic equilibrium (equal p, T, g): fast relaxation collapsed onto its equilibrium manifold | `phase_change.rs` |
| a posteriori L² / H¹ certificates | residual of the computed solution | refine where the certificate is weakest | `certify.rs` |

The rows above are the original solver core. The crates added since (materials,
waves, light, heat, interfaces, soft matter, electrochemistry) each carry the
same three-part shape, so they are listed by domain below. A **✓** marks a rung
with a checked-in benchmark and a stated error against an exact or published
result; **◇** marks a rung that is documented or analytic-only, or whose check
is order-of-magnitude or a mapping rather than a solver test.

### Solids, crystals and mechanical metamaterials

| coarse model | indicator | unfold / closure | where |
|---|---|---|---|
| homogeneous elasticity: `C_ij` and Christoffel sound speeds | relative dispersion error `(ω_lat − v k)/ω_lat`, `~ (k a)²` | full harmonic dynamical matrix `D(k)` | `materials/elastic.rs`, `continuum.rs`, `phonon.rs` ✓ |
| isotropic polycrystal (Voigt/Reuss/Hill, Hashin–Shtrikman) | grain size / wavelength `d/λ` (Rayleigh, stochastic, geometric) | single-crystal anisotropy, grain by grain | `materials/polycrystal.rs` ✓ |
| temperature-dependent EOS by quasi-harmonic free energy | `T/θ_D` (quantum/classical) and Lindemann `√⟨u²⟩/d_nn` | molecular dynamics / anharmonic rung | `materials/quasiharmonic.rs` ✓ |
| Cauchy–Born finite-element solid (quasicontinuum) | per-element Cauchy–Born energy error; strain jump | promote the coarse cell to fully atomistic (QNL) | `materials/quasicontinuum.rs` ✓ |
| GPR elasto-plastic solid (`TauModel::Plastic`) | Peierls stress `σ_P`; sample / dislocation spacing `L√ρ` | discrete cores and source-limited “smaller is stronger” | `materials/plasticity.rs` ◇ (PN crude) |
| homogenised lattice stiffness `C*` (architected) | cell count `N` (`~1.13/N` clamped boundary layer); `h_cell/λ`; strut stress / Euler load | discrete frame lattice; Bloch bands per cell | `architected/cell.rs`, `fold.rs`, `buckling.rs` ✓ |
| isobaric macroscopic cell with no interface | tip diameter / spacing vs `dx`; `W0/d0` | anisotropic Karma–Rappel phase field | `solidification/closure.rs`, `phase_field.rs` ✓ |
| Lennard-Jones molecular dynamics (EOS + Green–Kubo `η, κ, D`) | block/seed statistics; Knudsen `Λ/L`; cutoff `rc` tail | continuum EOS/transport closures; GPU sampling for 10⁵ atoms | `md`, `md-cuda` ✓ |
| 5-pt Karma–Rappel phase field without the variational term | tip-radius fit (`NaN`) and drifting `V d0/D` | 9-pt Laplacian + variational anisotropy (KR production operator) on GPU | `solidification-cuda` ✓ |

### Waves and metamaterials

| coarse model | indicator | unfold / closure | where |
|---|---|---|---|
| homogenised laminate (Reuss) | `(ω Λ/c_eff)²` | Rytov transfer matrix / Bloch dispersion | `phononic/homogenize.rs` ✓ |
| infinite-crystal band structure (one unit cell) | `N · abs(Im q) Λ` | finite stack | `phononic/layered.rs` ✓ |
| lumped Helmholtz `M_a, C_a` | `max(L_eff, L_c)/λ` | exact duct transfer matrix | `phononic/helmholtz.rs` ✓ |
| Fang effective bulk modulus `B_eff` | `k a = ω a/c` | exact cell matrix / Bloch wavevector | `phononic/helmholtz.rs` ✓ |
| rigid-wall cavity compliance | wall/fluid compliance `3B0/(4G)` | elastic-wall (Lamé) cavity compliance | `phononic/helmholtz.rs` ✓ |
| 2D phononic band structure (plane-wave expansion) | plane-wave count `m` (gap-edge drift) | refine `m` or extrapolate in `1/m` | `phononic/pwe.rs` ✓ |
| lossless 4×4 fluid–plate propagator | evanescent exponent `abs(Im k_zS) d` (conditioning) | global (stiffness) matrix | `phononic/fluid_coupling.rs` ✓ |
| continuum heat equation with `κ(T, L)` | phonon Knudsen number `Kn = Λ/L` | grey Boltzmann equation (exact integral equation) | `heat/effective.rs`, `bte1d.rs` ✓ |

### Light

| coarse model | indicator | unfold / closure | where |
|---|---|---|---|
| homogeneous uniaxial / 2D effective medium | `Λ/λ` (error `~ 2(Λ/λ)²(1+c)`); `homogenisation_error_2d` | transfer matrix / plane-wave expansion | `photonics/effective_medium.rs` ✓ |
| incoherent intensity sum (geometrical optics) | optical thickness `n d/L_c` | coherent transfer matrix | `photonics/incoherent.rs` ✓ |
| infinite 1D Bloch crystal (one unit cell) | `N · abs(Im K) Λ` | finite mirror | `photonics/bloch.rs` ✓ |
| 2D photonic-crystal band structure (PWE) | plane-wave count / eigenvalue drift | refine the plane-wave basis | `photonics/pwe.rs` ✓ |
| effective `ε(ω), μ(ω)` negative-index medium | loss `Im ε, Im μ`; `k_max d ~ ln(1/γ)` | exact multilayer solve | `photonics/metamaterial.rs` ✓ |
| local frequency-only `ε(ω)` polariton | `q a` (spatial dispersion); film/lattice `a/d` | lattice `D(k)` and the Born charge / LST | `photonics/polariton.rs` ✓ / ◇ |
| isotropic 2×2 TE/TM transfer matrix | anisotropy `abs(ε_∥−ε_⊥)/abs(ε)` | 4×4 Berreman matrix | `photonics/berreman.rs` ✓ |
| cholesteric homogenised to an effective medium | pitch/wavelength `p/λ` | fine Berreman stack | `photonics/berreman.rs` ✓ |
| adiabatic (Mauguin) director waveguiding | inverse Mauguin number `1/u` | exact 4×4 stack | `photonics/berreman.rs` ✓ |
| Rayleigh dipole scattering | `max(x, (abs m) x)` (relative error `~ indicator²`) | exact Mie series | `photonics/mie.rs` ✓ |
| anomalous diffraction (van de Hulst) | phase delay `ρ = 2x abs(m−1)` | exact Mie series | `photonics/mie.rs` ✓ |
| large-particle extinction `Q_ext → 2` | `1/x` | exact Mie series | `photonics/mie.rs` ✓ |
| perfect conductor `abs(m) → ∞` | `1/abs(m)` | finite-index Mie series | `photonics/mie.rs` ✓ |
| participating medium from particle optics | size distribution + scattering coefficients | exact Mie `Q_ext`, `g`, `p(θ)` | `photonics/mie.rs` → `radiation/` ✓ |

### Heat

| coarse model | indicator | unfold / closure | where |
|---|---|---|---|
| instantaneous heat deposition (no heat flow during the pulse) | thermal confinement `τ_p α/δ²` | solve conduction during the pulse | `photothermal/thermal.rs`, `indicators.rs` ✓ |
| instantaneous pressure `p0 = Γ μ_a F` | stress confinement `τ_p c/δ` | coupled thermo-acoustic solve during the pulse | `photothermal/photoacoustic.rs` ✓ |
| conduction-only transport | Rayleigh number `Ra` (`1708` surface, `2772` volumetric) | flow solver (convection) | `photothermal/indicators.rs` ✓ |
| Beer–Lambert source from `I0 A` | layer thickness vs absorption depth `δ = 1/μ_a` | `photonics` power-density profile | `photothermal/source_from_optics` ✓ |
| equilibrium diffusion → flux-limited diffusion → M1 → Monte Carlo | optical depth, gradients, angular structure | next transport rung | `radiation.md` (existing) ◇ |

### Fluids and interfaces

| coarse model | indicator | unfold / closure | where |
|---|---|---|---|
| diffuse-interface five-equation mixture | cells across each material layer | transfer-matrix homogenisation; tracked sharp interface | `compressible/multiphase.rs` ✓ |
| tracked sharp solid–fluid interface (1D) | cells per plate | exact solid–fluid Riemann flux | `compressible/sharp_interface.rs` ✓ |
| tracked sharp solid–fluid interface (2D, plane strain) | cells per wavelength; first-order transverse closure | P/SV transverse characteristic source in the flux; exact half-space `abs(R)(θ)` reference | `compressible/sharp_interface2d.rs` ✓ |
| effective contact angle `θ*` (Wenzel, Cassie–Baxter) | drop base radius / texture period `a/λ` | explicit textured-substrate free-energy minimisation | `wetting/texture.rs`, `lattice.rs` ✓ |
| no phase change | `p < p_sat(T) − tension` | equal-`p,T,g` equilibrium collapse | `compressible/phase_change.rs` ✓ (existing) |
| macroscopic isobaric cell with no interface | tip diameter / spacing vs `dx_macro`; `W0/d0` | Karma–Rappel phase field; Ivantsov + solvability closure | `solidification/closure.rs` ✓ |
| 1D Lagrangian isobaric enthalpy (planar front) | grid orientation `Δξ/Δx`; density jump `Δρ/ρ`; Rayleigh `Ra` | 2D Eulerian enthalpy with exact inversion; Boussinesq buoyancy source | `compressible/isobaric2d.rs` ✓ |
| conduction-only melt (isobaric 2D) | Rayleigh number `Ra` vs `1708` (layer) | 2D Boussinesq natural convection with Carman–Kozeny mushy drag; de Vahl Davis Nu within ~1 % | `compressible/convection2d.rs` |
| constant-coefficient variable-density projection | density contrast `ρ_1/ρ_2`; CG residual | exact variable-coefficient projection (PCG) | `twophase.rs` ✓ |

### Soft matter

| coarse model | indicator | unfold / closure | where |
|---|---|---|---|
| few-mode Maxwell / Prony (Rouse fold) | mode separation; Deborah / Weissenberg number | full Rouse spectrum | `polymer/rouse.rs` ✓ |
| Hookean dumbbell / Oldroyd-B | Weissenberg `Wi = λγ̇` (coil–stretch at ½) | FENE / FENE-P | `polymer/fene.rs`, `continuum.rs` ✓ |
| point-dipole chain / Bingham electrorheology | Mason number `Mn` (chains valid `Mn ≲ 0.1`) | particle simulation | `electrofluid/er.rs` ✓ |
| Taylor–Melcher leaky-dielectric drop | `Ca_E`; small deformation `D ≪ 1` | free-boundary Stokes / full Maxwell | `electrofluid/leaky.rs` ◇ |
| Newtonian suspending fluid with `η(φ)` | `φ/φ_m`; particle Reynolds `Re_p`; Péclet `Pe` | exact Stokes stresslet → Batchelor–Green → Krieger–Dougherty; `μ(J)`, `μ(I)` granular rheology | `granular/` ✓ |

### Electromagnetism of matter

| coarse model | indicator | unfold / closure | where |
|---|---|---|---|
| quasi-static potential (`∇×E = 0`, `∇·J = 0`) | charge-relaxation time `τ_e = ε/σ` vs flow time | full Maxwell | `electrofluid` electric solve ✓ |
| Debye–Hückel linear double layer | `u = zeψ/kT`, `abs(u) ≪ 1` | Gouy–Chapman / Poisson–Nernst–Planck | `electrofluid/pnp.rs` ✓ |
| Helmholtz–Smoluchowski slip | `λ_D/h` | resolved Stokes + PNP body force | `electrofluid/eo.rs` ✓ |
| micromagnetics (`A, K, M_s`, LLG) | wall width / lattice `δ/a`; `l_ex/L`; `k a` | atomistic Heisenberg spins | `magnetism/micromag.rs` ✓ |
| Stoner–Wohlfarth macrospin | particle size vs single-domain `D_c`; wall fits `d > πδ` | micromagnetics | `magnetism/hysteresis.rs` ✓ |
| local-shape demagnetising field `H_d = −N M_s (m·n)n` | internal-field uniformity; size vs `l_ex` | FFT (Newell–Williams–Dunlop) demagnetising field | `magnetism/demag.rs`, `llg3d.rs` ✓ |
| Thevenin RC equivalent circuit | pulse length / `τ_RC` | single-particle model / DFN | `electrochem/ecm.rs` ✓ |
| single-particle model (SPM) | electrolyte time / discharge time `t_e/t_d`; reaction non-uniformity `std(a_s j)/mean` | electrolyte-enhanced SPM; pseudo-2D Doyle–Fuller–Newman `dfn` | `electrochem/spm.rs`, `dfn.rs` ✓ |
| Nernst / Butler–Volmer kinetics | `nF abs(η)/RT` | full Butler–Volmer; Koutecký–Levich | `electrochem/interface.rs` ✓ |

### Electronics

| coarse model | indicator | unfold / closure | where |
|---|---|---|---|
| gate-level NLDM logic on timing/energy tables | input slew / load leaves the characterised box; accumulated interpolation bound | transistor-level EKV backward-Euler | `semiconductor/gates.rs` ✓ |
| compact MOSFET table (PCHIP) | held-out EKV points | all-region EKV reference | `semiconductor/compact.rs` ✓ |
| architecture / RTL / circuit / device / quantum | timing slack, noise margin, model range, mean free path | next finer rung (CPU design ladder above) | `emergence_ladder.md`, CPU section ◇ |

### Materials-ladder integration

| coarse model | indicator | unfold / closure | where |
|---|---|---|---|
| one material through every rung (atomistic → single crystal → polycrystal → continuum) | each rung's own indicator, reported per measurement | the crate-appropriate finer rung, with a `fold_valid` flag | `testlab` ✓ (see [materials_overview.md](materials_overview.md)) |

There are three kinds of emergent boundary:

- **Spatial:** a thin structure the grid cannot resolve (shock, interface,
  Hartmann layer, boundary layer, bubble smaller than a cell). The closure
  is a jump or matching condition.
- **Temporal:** a fast process (phase change, pressure relaxation, Joule
  damping, chemistry, radioactive decay). The closure is the slow manifold
  it relaxes onto, or its exact solution.
- **Statistical:** many unresolved degrees of freedom (molecules →
  continuum; small eddies → turbulence closure; nuclei → nucleation rate).
  The closure is a constitutive law with a known validity range.

## Cross-rung couplings that now exist

The rungs are not independent: each fold produces a closure that a cruder model
consumes. The couplings that exist in the crates today, and whether a
test exercises the hand-off, are:

- **phonons → polaritons → optics.** `materials::phonon` gives the zone-centre
  TO frequency; `photonics::polariton` identifies it with the pole of `ε(ω)`
  and the Born charge supplies the LO shift and the LST relation. Validated by
  the diatomic-chain test (pole matches `ω_TO` exactly); the Born-charge
  prediction itself is ◇.
- **polycrystal → GPR solid.** `materials::polycrystal::isotropic_parameters_stiffness`
  returns `ρ, c_s, c_b, v_p`, which the GPR adapter turns into a lossless elastic
  solid without a new dependency. Validated end-to-end through `testlab`'s
  plate impact (Cu/Al shock speeds within 0.6%).
- **MD → transport closures → Navier–Stokes.** `md::transport` produces `η, κ, D`
  and `md::validity` produces the Knudsen number that says when the continuum
  with those closures is legal. The CPU closures are compared with argon and
  JZG references, and `md-cuda` reproduces the CPU forces to a few `10⁻¹⁵` and
  the triple-point `η* = 3.18 ± 0.09` against Meier's `3.258`; the coupling to
  the flow solver is documented ◇ (no shared test).
- **photonics absorbed power → heat → photoacoustics → phononic.**
  `photonics::multilayer::power_density` is integrated by
  `photothermal::thermal`, the stress-confined source `p0 = Γ μ_a F` launches
  the acoustic field, and `photonics`/`phononic` transfer matrices give the
  layered waveform. The source integral `∫p dz = I0 A` (5.5e-6) and the
  pure-delay spectral layer (3.0e-14) are cross-crate tests. ✓
- **Mie → radiation participating medium.** `photonics::mie` supplies
  `σ_abs, σ_sca, g, p(θ)` for the droplet distribution, mapped onto
  `radiation::Medium`. Through the emergence engine the hand-off is now
  tested end to end (`testlab/tests/emergence_links.rs`): per-cell size
  distributions → Mie ladder → `radiation::Volume` → ray transmittance, with
  each droplet evaluation budgeted by its adjoint weight; worst transmittance
  error vs exact Mie everywhere 2.3e-14 under a bound ≤ 1e-6. ✓ (`p(θ)` is
  not yet passed; the medium uses `g`.)
- **phonon spectrum → film conductivity.** `heat::bte_ladder::film_kappa`
  suppresses every mode of the kinetic spectrum by the grey-slab ladder
  (Fourier / Milne jump / ballistic / closure / integral equation) under one
  goal-oriented error budget; claims checked against the integral equation.
  ✓ for the ladder error; the per-mode grey-suppression model itself is not
  validated against measured films.
- **defects → plasticity → GPR.** `materials::plasticity::GprPlastic` maps
  one-to-one onto `compressible::gpr::TauModel::Plastic` (the dislocation-drag
  coefficient `B` is an external closure). Read-only mapping; ◇.
- **quasicontinuum → Cauchy–Born elasticity.** The finite-element continuum
  energy is the Cauchy–Born `W(F)` of `materials::elastic`, so the same module
  supplies both the atomistic fold and its continuum limit. ✓
- **phonon heat → continuum/semiconductor thermal.** `heat::effective::effective_kappa(T, L)`
  returns `κ(T,L)` and a Fourier-validity flag for `∇·(κ∇T)`; the existing
  `semiconductor/src/thermal.rs` can call it unchanged. Documented ◇.
- **electrochem heat → heat.** `electrochem::thermal` supplies the Bernardi
  source density for `heat`'s `ρ c_p ∂T/∂t = ∇·(κ∇T) + q_vol`. Documented ◇.
- **wetting → capillary wall condition.** `wetting::texture` supplies `θ*`
  for the phase-field wall condition `κ ∂_n α = ±γ cosθ*` of
  `compressible/capillary.rs`. Documented ◇.
- **solidification closure → isobaric macro cell.** `solidification::closure`
  supplies the sub-grid tip kinetics and Avrami rate for the isobaric enthalpy
  model of `compressible` §15. Documented ◇.
- **testlab ties the chain together.** `testlab` runs one material (Ar, Cu, Al)
  through every rung — LJ/Cauchy–Born, single crystal, polycrystal and the GPR
  continuum — and prints each fold's error and a `fold_valid` flag, which is
  the end-to-end integration test of the materials ladder. ✓

## Validated, documented, and the honest misses that remain open

The tables mark validated rungs (✓) separately from documented or analytic-only
ones (◇). The following are the honest misses, gathered from the analysis notes;
numbers here are quoted from those notes and nowhere else.

- **Dendrite tip velocity — converging, not yet resolved.** After the operator
  fix (Karma–Rappel variational term), steady tip speeds converge toward KR
  with grid refinement: `D = 2` `0.01536 → 0.01663` and `D = 3`
  `0.01490 → 0.01600` for `dx/W0 = 0.4 → 0.2` (KR `0.0168` / `0.0175`; a linear
  extrapolation of `D = 3` gives `≈0.0173`). An earlier `D = 2` "agreement" at
  `0.0164` was a transient. With a robust tip-radius fit the selection
  constant is `σ* = 0.116` (`D = 2`) and `0.11–0.13` (`D = 3`), consistent with
  the sub-grid closure's `0.1288` (an earlier `0.5–0.7` was a radius fitted
  inside the diffuse cap)
  (`analysis/solidification_cuda.md`).
- **Peierls–Nabarro is a crude closure.** The atomistic `γ_us` is 6.5× smaller
  than the elastic/sinusoidal amplitude that would give the classical core
  width, and `σ_P/G = 0.081` against `10⁻³–10⁻⁴` for real fcc metals; the
  pair-potential stacking-fault energy is negative (`analysis/defects.md`).
- **Lattice trapping is not established.** The discrete crack shows only a
  metastable local bond-rupture window (`K_+/K_- ≈ 1.17`), not the true
  load-controlled trapping bounds; the saddle search did not converge
  (`analysis/fracture.md`).
- **Pair potentials lack directional bonding.** They overestimate vacancy
  energies, give a negative `γ_ISF`, and model diamond as central-force floppy
  (`analysis/materials.md`, `defects.md`).
- **Quasi-harmonic limits.** The Grüneisen parameter is +4 % to +20 % above the
  rare-gas band; the harmonic spinodal precedes melting by ~25 K in LJ argon;
  the direct low-`T` Debye extraction is grid-limited
  (`analysis/quasiharmonic.md`).
- **Silicon heat conductivity is not predicted.** The published parameters give
  `127 W/(m K)` from the RTA (−14 % vs 148) and `218` with the Callaway
  correction (+47 %); the isotope effect is 3–4× too small
  (`analysis/heat.md`).
- **Iron γ→δ volume.** The temperature-dependent thermal EOS gives `−0.03 %`
  against the measured `+0.53 %`, and `ΔH(γ→δ) = 0.554` against
  `0.84 kJ/mol`; the earlier attribution to a single NASG expansion
  coefficient is refuted — the published bcc and fcc EOSs themselves put
  `V_δ ≈ V_γ` at 1667 K (`analysis/compressible.md` §15).
- **Capillary-wave dispersion is a finite-interface-thickness effect.** The
  variable-coefficient projection is now implemented and helps at ratio 1000
  (`−16.8 % → +11.4 %`), but the restored stiffness is
  `ω²/ω²_sharp = 1 − 0.6267 ηk` (Jacqmin 1999), measured slope `−0.627`,
  independent of `η/Δx` at fixed physical `η`; the remedy is a thinner
  interface (`analysis/twophase_and_mhd.md`).
- **Isobaric 2D: expansion flow now resolved; grid orientation remains.** The
  volume change on melting/freezing now drives an irrotational expansion flow
  (Stefan fronts: iron `+10.65 % → +1.23 %`, water `−14.53 % → −0.71 %`; water vs
  the 1D Lagrangian model `−0.083 %`); a 45° front still shows a `+9.9`
  percentage-point grid-orientation error (`analysis/isobaric2d.md`).
- **FSI at full solid/fluid contrast.** The diffuse-interface shear-bearing
  two-material GPR extension is linearly unstable at water/steel and
  over-reflects (`abs(r)` up to `4.05` vs Fresnel `0.94`). The sharp 2D
  solid–fluid interface now converges obliquely — orders `1.37–2.04` across
  `0–32°` — but its transverse closure is only first-order and coarse-grid
  errors reach 38 % at 16 cells/λ; a tangential-wavenumber interface operator
  is the next rung (`analysis/fsi.md`, `analysis/sharp_interface2d.md`).
- **Battery absolute voltage not validated.** The pseudo-2D DFN is now
  implemented and satisfies lithium conservation, charge balance, the linear
  low-rate reduction to SPM/SPMe and first-order convergence; but no published
  DFN voltage curve was available offline, the positive-electrode OCP is a
  surrogate sigmoid, and several rate/transport constants are estimates
  (`analysis/electrochem.md`).
- **Taylor drop is analytic only.** The flow, circulation and small-deformation
  shape are closed-form and only numerically differentiated, not solved on a
  free boundary; the earlier `4/3` discrepancy from Taylor's `D` is the exact
  geometry conversion `D = 3s/4`, not missing physics
  (`analysis/electrofluid.md`).
- **Dense Krieger–Dougherty misses the Brownian divergence.** KD at
  `φ_m = 0.63` gives `5.24` and `11.19` against published low-shear `11` and
  `49` (−52 %, −77 %) at the dynamic crossover and melting
  (`analysis/granular.md`).
- **Coated-sphere Mie is not implemented** (`analysis/mie.md`).
- **Rayleigh trap stiffness is an overestimate.** The point-dipole gradient
  force gives `0.875 pN/nm` for a 1 µm bead, 10–100× above measured
  `0.01–0.1 pN/nm` per 100 mW; only the exact Mie `Q_pr` is reliable
  (`analysis/mie.md`).
- **Magnetism: µMAG standard problem 4 now converged.** With the FFT
  demagnetising field and cells below the exchange length, the 200×50×1 run
  gives the s-state `(0.9664, 0.1259)` vs NIST `(0.9666, 0.1260)` and first
  `⟨m_x⟩/⟨m_y⟩` zero crossings within `0.5 %`; the Curie-point Binder crossing is
  still a small-lattice `1.7 %` result (`analysis/magnetism.md`).
- **Phononic PWE converges slowly.** At steel/epoxy `f = 0.4` the SH gap edges
  still move ~1 % between `m = 9` and `m = 11`, so the reported gap is a bounded
  extrapolation, not a converged number (`analysis/phononic.md`).
- **Effective-medium indicators are switches, not certificates.** The photonics
  and phononic quadratic error envelopes can be exceeded at high contrast or
  near a constituent resonance (`analysis/photonics.md`, `phononic.md`).
- **Wetting hysteresis is not predicted.** `Δcosθ = −C φ_s` carries a supplied
  pinning strength, and the fold error is non-monotone because the contact line
  pins metastably on a sharp texture (`analysis/wetting.md`).
- **Electrorheology Mason-number conventions differ by factors of two**, so the
  quoted onset (`Mn ≈ 0.1–0.2`) carries that ambiguity (`analysis/electrofluid.md`).
- **Architected lattices are linear-only.** The octet truss is cubic, not
  isotropic; there is no post-buckling path or joint compliance
  (`analysis/architected.md`).
- **Gate-level energy is a loose bound.** The adder energy is +30 % high and the
  accumulated interpolation certificate does not cover the single-edge
  model-form error on long paths (`analysis/gates.md`).
- **Low-Rm MHD induction is not built, and Hartmann layers are sub-grid** so
  walls need a Hartmann-layer model (`emergence_ladder.md`, `twophase_and_mhd.md`).
- **Standard-Model hadronisation and off-shell decays remain explicit closures**
  rather than resolved physics (`standard_model.md`).

## Emergent optimizations not yet exploited (by expected payoff)

1. **Multiresolution grids with a guaranteed error bound** (Harten, and
   adaptive wavelet collocation). Wavelet detail coefficients below a
   threshold ε are dropped, and the perturbation error is bounded by C·ε.
   This is spatial adaptivity with a proof, the grid analogue of the
   shock sensor. The bubble collapse and every shock problem spend most of
   their cells on uniform regions today.
2. **Slow-manifold projection as a general rule.** The phase-change
   relaxation and the exact Joule factor are both instances: when a
   process is much faster than the step, replace its dynamics with its
   exact solution or its equilibrium manifold. The same applies to
   chemistry, decay chains, ionization equilibrium (Saha) and thermal
   relaxation between species. Implicit or exponential treatment removes
   the stiffness limits we hit (for example the dilatation limit in a
   bubbly mixture).
3. **Local time stepping / multirate.** The Wood-speed region of a bubbly
   mixture and the liquid around it differ by 10–100× in wave speed; a
   global time step wastes that.
4. **Sub-grid models matched to the resolved field:**
   - Lagrangian Rayleigh–Plesset bubbles that hand over to the resolved
     five-equation model when they grow past a few cells;
   - Hartmann-layer and turbulent wall models;
   - DSMC patches where a shock is thinner than about 20 mean free paths.
5. **Learned closures with certificates (the AlphaFold lesson).** AlphaFold
   works because the structure emerges from rules too expensive to run
   directly, yet is learnable from data. The analogues here:
   - equations of state and transport coefficients for every element and
     mixture, learned from ab initio (DFT / molecular dynamics) and
     experimental data;
   - nucleation rates;
   - turbulence closures.

   The rule that keeps it "provably correct": a learned closure is used
   only inside a domain where a cheap check certifies it (residual against
   conservation laws, thermodynamic consistency, conformal error bounds
   from held-out data). Out-of-distribution means an emergent boundary, so
   fall back or unfold.
6. **Emergent coordinates.** Coherent structures have few degrees of
   freedom: vortex filaments, wave packets (the planned WKB fold),
   interfaces as surfaces, bubbles as (R, Ṙ). Evolving them in their own
   coordinates is the largest possible compression when the indicator says
   they stay coherent.
7. **Hardware.** GPU port of the compressible solver (structure-of-arrays
   layout, the per-face hybrid switch as warp-uniform branches).

## Interactions as emergent systems

Each fundamental interaction reaches the continuum in three forms:

| interaction | force in the momentum/energy equations | constitutive closure | source terms |
|---|---|---|---|
| electromagnetic | Lorentz J×B, electrostatic ρ_e E, Maxwell stress | conductivity, permittivity, viscosity and EOS (electrons and phonons), surface tension | Joule heating, radiation |
| gravity | ρg (self-gravity: Poisson) | — | — |
| strong | — | nuclear EOS (dense matter) | fusion, fission |
| weak | — | — | β decay, neutrino emission and absorption |

Coupling all of them is operator splitting with per-term error estimates,
which is the same machinery.

## Nuclear and particle physics: what is physically sensible

- **Radioactive decay.** Bateman equations per nuclide: linear, with an
  exact solution (matrix exponential of the decay-chain matrix), so they
  are integrated exactly at any time step, the temporal case above.
  Half-lives and branching ratios come from evaluated nuclear data
  (ENSDF/NUBASE). Energy deposition: α and β locally as heat; γ and
  neutrinos by transport or escape.
- **Fusion.** Thermonuclear reaction rates ⟨σv⟩(T_i) (Bosch–Hale fits for
  D–T, D–D, D–³He) as species sources in a two-temperature (ion/electron)
  MHD plasma; α-particle heating; bremsstrahlung and synchrotron losses.
  Validation: ignition and Lawson-criterion conditions. Plasma indicators
  are the Debye length, plasma parameter and collisionality; when the
  plasma is collisionless the fluid model fails and kinetic (Vlasov/PIC)
  patches are the unfold.
- **Radiation.** Photon transport as a ladder: equilibrium diffusion →
  flux-limited diffusion → M1 moments → discrete ordinates / Monte Carlo,
  switched by optical depth.
  The first implementation now has conservative implicit grey exchange,
  closed 1D thick diffusion, 3D angular packets and a spectral camera;
  see [radiation.md](radiation.md). Eligibility also checks gradients and
  angular structure. FLD/M1 and automatic 3D handoffs remain future work.
- **Standard Model particles.** In a continuum simulator they appear as
  species (e⁻, e⁺, p, n, nuclei, γ, ν, …) with masses, charges and
  reactions taken from PDG data, and as closure data: cross sections,
  decay rates, and the equation of state of extreme matter (quark–gluon
  plasma and neutron-star matter from lattice QCD and chiral effective
  field theory tables). Simulating the quantum field theory itself costs
  exponentially more (a lattice-QCD calculation of one hadron property
  takes a supercomputer allocation), so QFT is an offline closure
  generator, never part of the fluid's time loop. That is not a
  compromise: it is the emergence ladder applied to the bottom rungs.
  The first kinetic implementation now provides every elementary species and
  antiparticle, relativistic flight, a conservative decay graph, selected QED
  collisions and an ideal quantum-gas closure. Unresolved hadronization and
  off-shell decays retain an explicit conservation ledger rather than being
  converted to heat. See [standard_model.md](standard_model.md) and the
  [module's domain limits](../fluid/crates/standard-model/README.md).

## Simulating a CPU design: the same ladder

A processor is the clearest engineered example of emergence: about 10¹⁰
transistors whose collective behaviour is designed to be describable by
Boolean logic and clocked state machines. Simulating it with the
coarsest valid model per region, and certified hand-offs, is exactly
this framework.

| rung | model | valid when (indicator) | unfold when |
|---|---|---|---|
| architecture | instruction / cycle-level state machine | timing closure holds (all path slacks > margin) | a path's slack approaches 0 |
| RTL / gate level | Boolean gates, delays from library tables | noise margins hold, signals settle within the clock | glitches, metastability, crosstalk |
| circuit (SPICE) | compact transistor models (BSIM-CMG / FinFET), RC(L) interconnect | device inside the model's calibrated range | new device physics, extreme bias or temperature |
| device (TCAD) | drift–diffusion: Poisson + electron/hole continuity (+ hydrodynamic / energy transport) | channel ≫ mean free path, quasi-equilibrium carriers | ballistic transport, quantum confinement |
| quantum transport | non-equilibrium Green's functions, Schrödinger–Poisson | — | (bottom rung) |
| electromagnetic | quasi-static RC → RLC → full-wave Maxwell | wire length ≪ wavelength/10 | multi-GHz long lines, packaging |
| thermal and cooling | heat conduction in the die and package, coupled to the liquid/air coolant (this repository's fluid solvers) | — | — |

Drift–diffusion is the same machinery as the solvers here (an elliptic
Poisson solve plus advection–diffusion with stiff source terms,
Scharfetter–Gummel fluxes playing the role of the Riemann solver).
Electro-thermal coupling closes the loop: leakage grows exponentially with
temperature, temperature comes from the coolant flow, the coolant flow
from the fluid solver.

The certified hand-off is what is new. Commercial flows already use the
ladder (static timing, then SPICE on critical paths), but the switch is
set by hand. Here it would be driven by indicators with error bounds:
slack and noise margin per path decide where gate-level simulation
unfolds into circuit simulation, and device-model range checks decide
where circuit simulation unfolds into TCAD.

## Order of work

1. Compressible core, multiphase, cavitation (done: `fluid/crates/compressible`).
2. Elastic and plastic solids in the same hyperbolic framework (GPR unified
   continuum model: a solid is a fluid with a long strain-relaxation time).
3. Per-element property tables (periodic table) with provenance and
   uncertainty, and reference equations of state where they exist
   (IAPWS-95 water, Span–Wagner N₂).
4. Multiresolution grids and multirate time stepping (items 1–3 above).
5. Reactive species: decay chains (exact), chemistry; then fusion plasma
   and radiation transport.
6. Learned, certified closures where tables run out.
7. Electronics ladder: drift–diffusion device solver → compact models →
   circuit → gate level, with electro-thermal and coolant coupling.
