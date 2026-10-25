# Materials ladder: one-page map from atoms to continuum

This is the index for the solid-state side of the repository. The organising
idea is [`emergence_ladder.md`](emergence_ladder.md): a cheap model, a cheap
indicator, and an unfold to a finer model (or a closure) where the indicator
fails. Every materials crate is one rung of that ladder, and `testlab` runs one
material through all of them at once.

```
  Heisenberg spins  --exchange, small-angle-->  micromagnetics  --coarse-grain-->  macrospin
        |                                              |
        |  electromagnetic rung (magnetism.md)         v
  atoms (LJ / pair)  --Cauchy-Born-->  elasticity C_ij  --orientation average-->  isotropic polycrystal
        |                 \                    |                        \
        |  statistical     \  QHA             |  Peierls-Nabarro        \  Voigt/Reuss/Hill
        v  (md.md)          v                 v                          v
  molecular dynamics  -->  P(ρ,T), η, κ, D   defects / γ-surface  -->  GPR elasto-plastic solid
        |                        |                                          |
        |  heat rung (heat.md)   |  Navier-Stokes closures                  |  dynamic tests
        v                        v                                          v
  κ(T,L), phonon BTE      fluid / spectral solver                  plate impact, HEL, Hugoniot
```

The dashed/unfold path is equally part of the map: where the continuum is not
valid, the crates promote the local region back to atoms.

- **atomistic** — `materials` (lattice sums, Cauchy–Born, phonons),
  `quasiharmonic` (QHA EOS), `defects` (point/planar defects, γ-surface),
  `quasicontinuum` (adaptive atomistic/continuum), `md` (statistical closures);
- **continuum** — `materials::anisotropic`/`polycrystal` (measured
  single-crystal constants and the isotropic fold), the GPR solid inside
  `compressible`, and `testlab` (the standard test programme);
- **adjacent rungs** — `heat` (phonon `κ(T,L)`), `polymer` (viscoelastic
  closure), `magnetism` (spins → micromagnetics), `solidification` (phase field
  and mushy-zone closure), `architected` (beam-lattice homogenisation),
  `wetting` (textured-wall closure), `electrochem` (cell-level rungs).

## Crate map

| crate | what it does | key benchmark | error |
|---|---|---|---|
| `materials` | lattices, pair potentials, Cauchy–Born elasticity, phonons, anisotropic single crystals, polycrystal fold | Christoffel vs phonon speeds at `q a = 1e-3`; hcp–fcc energy difference | max `4.0e-8`; `9.0e-4 ε` (0.0104 %) |
| `materials::quasiharmonic` | harmonic free energy minimised over volume → `α, B_T, C_V, γ`, QHA EOS | solid-argon `a(4.2 K)`; `θ_D` | `+0.06 %`; `−0.2 %` |
| `materials::defects` | vacancy, surface, stacking-fault and γ-surface energies; Peierls–Nabarro and GPR plasticity mapping | rigid broken-bond vs relaxed slab surface energy; unrelaxed vacancy identity | `1e-6` relative; `8.6050` vs `8.6102 ε` |
| `materials::quasicontinuum` | adaptive atomistic/continuum (QNL) coupling | 1D patch test; adaptive vs uniform refinement to `10⁻⁴` | QNL ghost force `< 10⁻¹⁰`; 151 vs ≥ 283 ndof |
| `materials::fracture` | Griffith and anisotropic (Stroh) toughness, atomistic K-field cracks, Rice emission criterion, cohesive zone | Stroh vs isotropic `K_IC` on {100}; cohesive work = `2γ_s` | `−5.6 %` (anisotropy); `1e-6` (lattice-trapping range not yet measured validly, see [fracture.md](fracture.md)) |
| `md` | Lennard-Jones MD: EOS, Green–Kubo `η, κ, D`, continuum-validity indicators | JZG EOS `P*, U_res*`; zero-`T` fcc energy | `−3.2 %`/`0.006 %`; `3.2e-14` |
| `md-cuda` | GPU Lennard-Jones MD (cell list, f32/f64) and block/seed Green–Kubo `η`, `D` | GPU vs CPU forces; triple-point `η*` | few `1e-15`; `3.18 ± 0.09` vs `3.258 ± 0.033` |
| `architected` | periodic homogenisation of beam lattices to `C*`; buckling and Bloch | regular honeycomb `E*/Es`, octet truss | `< 1e-13`; `< 1e-9` |
| `testlab` | runs Ar/Cu/Al through every rung and reports value, reference, error, `fold_valid` | Cu Hill `E`; Cu dynamic particle velocity | `−0.51 %`; `+0.015 %` |
| `heat` | Callaway/Holland kinetic `κ(T,L)`, size effects, grey BTE fold | natural-Si `Γ`; grey-BTE flux vs exact integral equation | `0.09 %`; `< 0.6 %` |
| `polymer` | Rouse/Zimm chains → Maxwell/Prony/Oldroyd-B; FENE; GPR mapping | Rouse `τ₁`; Oldroyd-B steady `N₁` | `0.03 %`; exact |
| `magnetism` | Heisenberg spins → micromagnetics → Stoner–Wohlfarth; local-shape and FFT (Newell–Williams–Dunlop) demagnetising fields | Bloch wall `σ`; macrospin astroid; uniform-cube/Lorentz-sphere demag | `6.5e-4`; `5e-5` at intermediate angles; demag `< 1e-12` |
| `solidification` | Karma–Rappel phase field, sharp-interface Stefan, Ivantsov closure | 1D Neumann front; `R_c = d0/Δ` | `0.51 %`; `< 1 %` for `W0/dx ≥ 1.6` |
| `solidification-cuda` | GPU (NVRTC) phase field with the KR 9-point + variational operator | KR Table II steady tip speed vs grid | `D = 2`: `0.01663` at `dx/W0 = 0.2` vs `0.0168`; `D = 3`: `0.0160` (extrapolated `≈0.0173` vs `0.0175`) |
| `electrochem` | interface kinetics, 1D/spherical transport, SPM/SPMe/DFN and Thevenin ECM | Randles–Ševčík CV peak; Cottrell transient; DFN low-rate limit | `−0.0018 %`; `0.0046 %`; `0.25 mV` at C/100 |
| `wetting` | capillary statics and Wenzel/Cassie–Baxter wall closure | flat-wall apparent angle; pendant drop vs OpenDrop | `0.031°`; `3.8e-5` in `z` |
| `granular` | Stokes sphere → suspension viscosity (Einstein, Batchelor, Krieger–Dougherty); Hertz–Mindlin DEM → `μ(I)` | Hertz contact time; restitution targets | `< 1e-3`; 4 digits (dense KD −52…−77 %, see [granular.md](granular.md)) |
| `photonics` / `phononic` | wave rungs used by the materials chain (dispersion, homogenisation) | Mie vs Wiscombe; layered/Bloch identities | see [photonics.md](photonics.md), [phononic.md](phononic.md) |

Benchmarks and the honest misses of each rung are in the analysis notes named
after the crate; every number above is quoted from those notes.

## Running a material through the ladder

`testlab` is the single entry point. It is SI throughout and every measurement
carries the material, the test, the level, the value, a reference, the error and
a `fold_valid` flag.

```sh
cd fluid
cargo test -p testlab                              # 10 tests, ~11 s total
cargo run -p testlab --example report              # static + ultrasonic tables
cargo run -p testlab --example report -- --full    # plus MD and plate impact
```

The three materials are Cu and Al (measured single-crystal constants as the
bottom rung) and argon (a fully self-consistent LJ chain). For argon all three
folds close: Cauchy–Born `C_ij` to `0.11 %`, phonon vs Christoffel to `4.0e-8`,
and MD vs the JZG EOS to `+2 %` temperature / `+1 %` pressure. For Cu and Al the
atomistic rung is **not** closed — a central pair potential gives the wrong
Cauchy relation and a negative stacking-fault energy — so the measured constants
are the bottom rung, and the report says so.

A quick single-rung run:

```sh
cd fluid
cargo test -p materials                            # Cauchy–Born, phonons, polycrystal
cargo test -p materials --test quasiharmonic       # QHA EOS (~10 s)
cargo test -p materials --test defects             # defects and γ-surface (~1.7 s)
cargo test -p materials --test quasicontinuum      # QC coupling (~20 s)
cargo test -p md                                   # EOS and transport closures
cargo test -p heat                                 # κ(T,L) and the BTE fold
```

## Where the folds do not close

The materials chain is honest about its edges. The main open misses (all
detailed in the per-crate notes and gathered in
[`emergence_ladder.md`](emergence_ladder.md)):

- **Pair potentials are not metals.** They miss directional bonding, give a
  negative `γ_ISF`, and overestimate vacancy energies; Cu and Al therefore have
  no atomistic rung.
- **Peierls–Nabarro is a crude closure.** The atomistic `γ_us` is 6.5× smaller
  than the elastic amplitude that would give the classical core width, and
  `σ_P/G = 0.081` against `10⁻³–10⁻⁴` for real fcc metals.
- **Lattice trapping is not established.** The discrete crack shows only a
  metastable local bond-rupture window (`K_+/K_- ≈ 1.17`); the true
  barrier-defined bounds did not converge.
- **Quasi-harmonic limits.** The Grüneisen parameter is +4 % to +20 % above the
  rare-gas band; the harmonic spinodal precedes melting by ~25 K in LJ argon.
- **Silicon heat conductivity is not predicted.** `127 W/(m K)` from the RTA
  (−14 % vs 148) and `218` with the Callaway correction (+47 %); the isotope
  effect is 3–4× too small.
- **Dendrite tip velocity** converges toward Karma–Rappel with grid refinement
  (`D = 2`: `0.01663` at `dx/W0 = 0.2` vs `0.0168`; `D = 3`: `0.0160`,
  extrapolated `≈0.0173` vs `0.0175`); the selection constant `σ*` (`0.116` / `0.11–0.13`) agrees with the
  sub-grid closure's `0.1288`.
- **Iron γ→δ volume is off** (`−0.03 %` vs `+0.53 %`, `ΔH = 0.554` vs
  `0.84 kJ/mol`); the temperature-dependent EOS refutes the earlier
  single-NASG-coefficient explanation, and the architected octet truss is
  cubic, not isotropic.
- **Dense Krieger–Dougherty misses the Brownian divergence.** KD at
  `φ_m = 0.63` gives `5.24` and `11.19` against published low-shear `11` and
  `49` (−52 %, −77 %) at the dynamic crossover and melting.
- **µMAG standard problem 4** is now converged: at 200×50×1 the s-state is
  `(0.9664, 0.1259)` vs `(0.9666, 0.1260)` and the crossings are within `0.5 %`.
- **The battery's absolute voltage is not validated.** The DFN is implemented
  and its laws hold, but no published voltage curve was available offline and
  several inputs are estimates.
- **Rayleigh trap stiffness is an overestimate.** The point-dipole gradient
  force gives `0.875 pN/nm` for a 1 µm bead, 10–100× above measured
  `0.01–0.1 pN/nm` per 100 mW.

The sharp 2D interface, FFT demagnetisation and pseudo-2D DFN belong to their
own modules and are listed in [`emergence_ladder.md`](emergence_ladder.md).
