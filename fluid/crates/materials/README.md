# Materials: atomistic rung to continuum elasticity

The first rungs of the materials ladder in
[`analysis/emergence_ladder.md`](../../../analysis/emergence_ladder.md):
crystal lattices + pair potentials (fine) folded to continuum elasticity
(coarse) by the Cauchy–Born rule, then real single crystals folded to
isotropic polycrystals, with the harmonic phonon dispersion (and, at the
polycrystal rung, the grain-size/wavelength ratio) as the validity
indicator and the unfold.

- `lattice`: primitive vectors + fractional basis for sc, bcc, fcc, hcp
  (ideal or given `c/a`), diamond and a custom cell, plus true spherical
  neighbour shells and bond lists with periodic images.
- `potential`: Lennard-Jones and Morse with exact first and second
  derivatives, a harmonic test spring, exact lattice sums, a homogeneous
  LJ tail correction and equilibrium by minimising the lattice constant.
  The tail-consistent relaxation scales the cutoff with the trial lattice,
  so it reproduces the exact infinite-sum lattice sums (`r0`) rather than
  a cutoff-biased value.
- `elastic`: Cauchy–Born energy density, analytic Hessians, internal
  basis relaxation (needed for diamond/hcp), the relaxed 6×6 Voigt
  stiffness, an independent finite-difference check, and an
  extrapolated (infinite-cutoff) stiffness.
- `phonon`: the harmonic dynamical matrix `D(k)` for any lattice and pair
  potential, phonon frequencies, the exact 1D-chain dispersion and the fcc
  Γ–X–W–L–Γ path.
- `continuum`: the Christoffel matrix `C_ijkl n_j n_l / ρ`, its three
  sound speeds, and the relative dispersion error that decides when
  elasticity may replace the phonons.
- `anisotropic`: measured room-temperature single-crystal elastic
  constants and densities (cubic Si, Ge, diamond, Cu, Al, Ni, α-Fe, W,
  NaCl; hexagonal Mg, Ti, Zn), symmetry classes, Bond-matrix orientation
  rotation, the compliance `S = C⁻¹`, directional Young's modulus and the
  Zener ratio, and the Christoffel speeds and polarisations.
- `polycrystal`: Voigt/Reuss/Hill orientation averages, the
  Hashin–Shtrikman shear bounds for cubic polycrystals, isotropic moduli
  and sound speeds, the parameters the GPR solid model needs, and the
  grain-size/wavelength validity indicator with the Rayleigh attenuation.

```sh
cd fluid
cargo test -p materials
```

Validation is in `tests/atomistic.rs`, `tests/anisotropic.rs`,
`tests/polycrystal.rs` and `analysis/materials.md`. Exact checks: fcc/bcc/
diamond/hcp shell counts; fcc LJ `A6 = 14.45392`, `A12 = 12.13188`,
`r0 = 1.09017297 σ` (now to 1e-8), cohesive energy; converged LJ fcc
constants `C11 = 105.182`, `C12 = C44 = 60.187 ε/σ³` against the
infinite-sum reference `105.19/60.19 ± 0.02`; the central-force Cauchy
relation and cubic symmetry; the diamond floppy shear cured by internal
relaxation; the 1D chain `ω = 2√(K/m)|sin(ka/2)|` to machine precision;
long-wavelength phonon velocities equal to the Christoffel speeds to 4e-8
with a `(k a)²` error. Real-crystal benchmarks: silicon `v_L[100] =
8434.8`, `v_T[100] = 5846.2`, `v_L[110] = 9136.2 m/s` and `E = 130.1 /
169.1 / 187.9 GPa` for `[100]/[110]/[111]`; tungsten `A = 1.011`,
copper `A = 3.209`; Cu polycrystal `G = 47.3`, `E = 127.4 GPa`, Al `G =
26.1`, `E = 70.4 GPa`. Limits: central pair potentials only (no
angular/covalent forces), Born tensor without stress (Huang) terms, a
mean-field tail, a leading-order Hashin–Shtrikman shear construction and
a leading-order Rayleigh anisotropy factor (see the write-up for what is
approximate).
