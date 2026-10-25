# Wave optics: transfer matrices, photonic bands, and where to fold

Implemented in `fluid/crates/photonics`. This rung sits below the geometrical
optics of `fluid/crates/radiation`. Geometrical optics is exact only once two
things are unresolved: the phase of multiply-reflected partial waves, and the
sub-wavelength structure of the matter. This crate keeps both, then provides
the cheap indicator and the unfold/fold for each. Nothing in `radiation` is
changed; the two crates share conventions (complex `N = n + i kappa`, SI
wavelengths, normal-incidence power normalised to the surface irradiance).

## The models

A plane stack is described by layer characteristic matrices. With the signed
normal wavevector `q = N cos(theta) = N sqrt(1 - (s/N)^2)`
(`s = n0 sin(theta0)`, `Im q >= 0` for a passive positive-index medium, and the
sign of `N` kept for a negative-index medium), tilted admittance `eta = q/mu`
(TE) or `eps/q = N^2/(mu q)` (TM), and `delta = 2 pi d q / lambda`:

\[
  M_j = \begin{pmatrix}\cos\delta & i\sin\delta/\eta\\
  i\eta\sin\delta & \cos\delta\end{pmatrix},
  \qquad
  r=\frac{\eta_0-Y}{\eta_0+Y},\quad Y=\frac{m_{22}\eta_s-m_{21}}{m_{11}-m_{12}\eta_s}.
\]

`r`, `t`, `R`, `T`, `A`, per-layer absorptance and the internal field all come
from one matrix product. Non-magnetic layers (`mu = 1`) recover the usual
`q = sqrt(N^2-s^2)`, `eta = q` / `N^2/q`. The 1D photonic-crystal band
structure is the same matrix on one unit cell: `cos(K Lambda) = (1/2) tr M`.
Material dispersion uses causal forms only (Sellmeier, Drude, Lorentz,
split-ring), so the Kramers-Kronig relations hold by construction and can be
checked numerically.

## Benchmarks (measured)

Run `cd fluid && cargo test -p photonics`. Every value below is printed by the
tests; "error" is the absolute deviation from the exact or published result.

| Benchmark | Exact / published | Computed | Error |
|---|---|---|---|
| Fresnel single interface, TE/TM, 0–80° | closed form | transfer matrix | \(<10^{-13}\) |
| Brewster `R_p = 0` at `tan θ_B = n2/n1` (3 pairs) | 0 | \(<10^{-24}\) | round-off |
| Total internal reflection, `θ > θ_c` | 1 | 1 | \(<10^{-14}\) |
| Quarter-wave ARC `n_c = sqrt(n_s)`, design λ | 0 | 1.6e-34 | round-off |
| Bragg mirror, N = 1,2,4,8,16 pairs | \((1-\frac{n_s}{n_0}(\frac{n_L}{n_H})^{2N})^2/(\ldots)^2\) | transfer matrix | \(<3\times10^{-16}\) |
| Bragg stop-band, relative **frequency** width `(4/pi) asin(r)` | 0.291131190 | 0.291131190 (`tr M` band edges) | round-off |
| Same gap, relative **wavelength** width `Delta/(1-Delta^2/4)` | 0.297433612 | 0.297433612 | round-off |
| Energy `R + T + A = 1`, lossless and lossy, 0–50° | 1 | 1 | \(<10^{-12}\) |
| Per-layer absorptance sums to `A` | `A` | `A` | \(<10^{-12}\) |
| Reciprocity (lossless), `R_front = R_back` normal incidence | equal | equal | \(<10^{-12}\) |
| Bloch gap vs long finite mirror (N = 40) | gap `tr M` = 163.6 nm, centre 561.9 nm | finite `R>0.9` = 170.9 nm, centre 563.0 nm | centre 0.2%, width 4.5% |
| Silica `n(589 nm)` (Malitson 1965) | 1.4585 | 1.4584132 | \(8.7\times10^{-5}\) |
| Water `n(589 nm)` (Daimon–Masumura 2007) | ≈1.3334 | 1.333358 | — |
| Drude high-ω `eps -> 1 - ω_p²/ω²` | asymptotic | matches to \(10^{-6}\) | — |
| Kramers–Kronig, Lorentz oscillator | direct Re ε(ω) | KK integral | \(2.0\times10^{-4}\) |
| Gold/silver 633 nm from Rakić LD fit | ε_Au ≈ -9.8+2.0i | -9.811+1.964i | tabulated comparison |
| Heat source `∫ p dz` vs `I0 · A` | equal | equal | \(2.3\times10^{-9}\) |
| SPR (gold 50 nm on n=1.5 prism) analytic angle | 44.5734° | TMM dip 44.9492° | 0.3758° (loss + finite film) |

Fresnel, Brewster, TIR, ARC, Bragg and energy conservation are exact
identities, so their errors are round-off.

**The quarter-wave stop-band formula is exact, not leading order.** The
standard result `(4/pi) asin(r)`, with `r = (n_H-n_L)/(n_H+n_L)`, is the exact
relative **frequency** width `Delta_omega/omega0` of the infinite quarter-wave
stack to all orders in the index contrast. Writing `x = omega/omega0`, the
normal-incidence unit-cell Bloch trace is
`(1/2) tr M = cos^2(pi x/2) - ((rho+rho^-1)/2) sin^2(pi x/2)` with
`rho = n_H/n_L`; the edges `(1/2) tr M = -1` are at
`x = 1 +/- (2/pi) asin(r)`, hence `Delta_omega/omega0 = (4/pi) asin(r)`
exactly. At `n_H/n_L = 2.30/1.45` the formula gives `0.291131190`, and the
band edges found independently from `|tr M/2| = 1` give the same value to
round-off (measured difference `0`). The `0.297434` quoted earlier is the same
gap measured on the **wavelength** axis: `omega = 2 pi c/lambda` gives
`Delta_lambda/lambda0 = Delta/(1 - Delta^2/4)` (the Jacobian between the two
axes), which is `0.297433612` here. The previous 2.1% "miss'' was the
axis change, not model error. The two widths are exposed separately as
`quarter_wave_stopband_ratio` and `quarter_wave_stopband_wavelength_ratio`.

The water and gold/silver lines are consistency checks against the fits, not
independent measurements.

## Fold/unfold map

The organising idea is `emergence_ladder.md`. The wave-optics instance:

| Coarse model (fold) | Cheap indicator | Unfold when | Where |
|---|---|---|---|
| Homogeneous (uniaxial) effective medium | `Λ/λ`; `homogenisation_error = 2(Λ/λ)²(1+contrast)` | error exceeds target reflectance accuracy | `effective_medium` |
| Homogeneous 2D medium (area average / Maxwell Garnett) | `a/λ`; `homogenisation_error_2d` | error exceeds target band accuracy | `effective_medium` + `pwe` |
| Incoherent intensity sum (geometrical optics) | `n d / L_c` (optical thickness / coherence length) | `n d / L_c <~ 1`: fringes resolved | `incoherent` |
| Infinite-crystal band structure (one unit cell) | `N \|Im K\| Λ` (periods × decay per period) | `N \|Im K\| Λ <~ 1`: finite size dominates | `bloch` |
| Infinite 2D crystal (one unit cell, PWE) | plane-wave count / residual eigenvalue drift | gap edge moves with cutoff | `pwe` |
| Ideal lossless negative-index slab (perfect lens) | loss `Im ε`, `Im μ`; `k_max d ~ ln(1/γ)` | evanescent gain saturates at `k_max` | `metamaterial`, `multilayer` |
| Lumped heat source `I0 A / d` | layer thickness vs absorption length | absorption length < layer thickness | `multilayer::power_density_profile` |

### Effective medium: measured convergence

A laminate of `n1 = 2.0`, `n2 = 1.5`, fill `f = 0.4` (so
`eps_parallel = 2.95`, `eps_perp = 2.7273`), total thickness `2 µm`, at
`λ = 550 nm` and 34.4°, compared with the exact uniaxial slab
(`anisotropic_slab_reflectance`). Reflectance errors:

| `Λ/λ` | err TE | err TM | indicator `2(Λ/λ)²(1+c)` |
|---|---|---|---|
| 0.1136 | \(1.90\times10^{-2}\) | \(2.77\times10^{-3}\) | \(3.31\times10^{-2}\) |
| 0.0568 | \(4.24\times10^{-3}\) | \(6.37\times10^{-4}\) | \(8.26\times10^{-3}\) |
| 0.0284 | \(1.03\times10^{-3}\) | \(1.56\times10^{-4}\) | \(2.07\times10^{-3}\) |
| 0.0142 | \(2.56\times10^{-4}\) | \(3.89\times10^{-5}\) | \(5.17\times10^{-4}\) |
| 0.0071 | \(6.39\times10^{-5}\) | \(9.71\times10^{-6}\) | \(1.29\times10^{-4}\) |
| 0.0036 | \(1.60\times10^{-5}\) | \(2.43\times10^{-6}\) | \(3.23\times10^{-5}\) |

The error falls by 4× each time the period halves: clean `O((Λ/λ)²)`. The
indicator is an empirical envelope with the measured quadratic prefactor (the
leading coefficient is ≈1.25 for TE, ≈0.19 for TM at this contrast); it is a
switch, not a certificate, and is documented as such. Unfold to the transfer
matrix when it exceeds the required accuracy.

## Two-dimensional photonic crystals (plane-wave expansion)

`pwe` computes the band structure of a 2D crystal of circular rods or holes by
the plane-wave expansion method of Ho, Chan & Soukoulis (PRL **65**, 3152
(1990)). Both polarisations are discretised with the **inverse-`ε` (Ho) rule**,
which inverts the truncated dielectric matrix `E_{GG'} = ε_{G-G'}`:

```text
  TM (E along the rods):  A_{GG'} = |k+G| |k+G'| [E^-1]_{GG'}
  TE (H along the rods):  A_{GG'} = (k+G).(k+G') [E^-1]_{GG'}
```

The inverse rule is used for both because the dielectric function jumps at the
rod surface: the direct Fourier coefficients of `1/ε` decay only as `1/|G|`,
so a direct truncation converges slowly (the “staircase'' problem). The matrix
inverse reproduces the exact answer with a handful of plane waves. The
lowest bands are found by Lanczos with full reorthogonalisation; the start
vector is a deterministic pseudo-random vector so it cannot be orthogonal to
the ground state.

Benchmarks against Joannopoulos et al., *Photonic Crystals: Molding the Flow
of Light*, 2nd ed. (2008):

| Benchmark | Book | Computed | Error |
|---|---|---|---|
| Square alumina rods `ε=8.9`, `r=0.2a`: TM gap | 31.4% gap–midgap | `0.3224–0.4426`, 31.4% | <0.1% |
| Same structure: TE gap | none | none (bands 0–3) | — |
| Triangular air holes `ε=13`, `r=0.48a`: complete gap | 18.6% | `0.4299–0.5205`, 19.1% at `n=6` | 0.5% |
| Long-wavelength TM effective index | area average | matches | <1e-4 |
| Long-wavelength TE effective index | Maxwell Garnett | above MG, below arithmetic mean | — |

The square-rod TM gap and its 31.4% ratio are reproduced to better than 0.2%
from `n=3` (49 plane waves) onward, and the TE spectrum has no gap. The
triangular air-hole complete gap (the intersection of the TE and TM gaps,
which sit between **different** band indices) converges to `19.0%` as the
plane-wave count grows, within 0.5% of the book's 18.6%. The book's numbers
were produced by MPB with a smoothed interface; the small residual is that
smoothing.

**Note on the originally requested square-rod numbers.** The task quoted a TM
gap of `0.302–0.443` (38%). Those figures do not appear in the book: the book
states a 31.4% gap–midgap ratio for this structure, and the computed gap is
`0.3224–0.4426` (31.4%). The numerical result matches the published 31.4% and
is reported as measured rather than tuned to the quoted 38%.

### The 2D emergent boundary: long-wavelength homogenisation

The long-wavelength limit of the lowest band gives the effective-medium
index. For **TM** (E along the rods) the exact result is the area average
`ε_eff = f ε_rod + (1-f) ε_bg`; the PWE slope at Γ reproduces it to `1e-4`.
For **TE** the dilute closure is the 2D Maxwell Garnett formula for circular
cylinders (`depolarisation 1/2`); the PWE value lies just above it and below
the arithmetic mean, the expected lattice correction. The validity indicator
is the same quadratic envelope as in 1D, `homogenisation_error_2d =
2 (a/λ)^2 (1+contrast)`: at `k=0.25 (2π/a)` (`λ=4a`) the TM band departs from
the effective-medium line by `1.1e-2`, well inside the indicator `2.2e-1`.
Unfold to the full PWE when the indicator exceeds the required band accuracy.

## Metamaterials with a magnetic response

The multilayer code carries a complex relative permeability `μ`. The tilted
admittance generalises to `η = q/μ` (TE) and `ε/q = N^2/(μ q)` (TM), and the
normal wavevector uses the signed branch `q = N cos θ = N sqrt(1-(s/N)^2)`.
The sign matters only for a negative-index medium: for `N<0` an evanescent
wave then grows in `+z`, which is exactly Pendry's perfect lens.

A wire medium (`Drude ε`) plus a split-ring resonator (`Lorentz μ`, Pendry
1999) gives a double-negative window. The model `Metamaterial` exposes `ε`,
`μ`, the signed index `n = sqrt(ε) sqrt(μ)` and the impedance, and locates the
negative-index band. For `ω_p = 1.2`, `ω_0 = 1`, `F = 0.8`, `γ = 0` the band
is `ω = 1.00–1.20`, where `Re ε < 0`, `Re μ < 0` and `Re n < 0`.

| Benchmark | Exact | Computed | Error |
|---|---|---|---|
| Lossless `ε=μ=-1` slab: `R`, `T` at 0–75° | 0, 1 | 0, 1 | <1e-12 |
| Transmission phase vs equal vacuum gap | `-k0 d` | `-k0 d` | round-off |
| Evanescent `\|t\|`, `s=1.2,1.5,2,3` | `exp(sqrt(s^2-1) k0 d)` | same | <1e-9 |
| Passive lossy slab: `R+T+A` | 1 | 1 | <1e-12 |
| Resolution limit `k_max d` vs `ln(1/γ)` | 9.2,6.9,4.6,2.3 | 10.0,7.8,5.6,3.5 | — |

The `ε=μ=-1` slab is impedance matched at every angle, so `T=1`, `R=0`, and
its transmission phase is the negative of an equal vacuum gap. In the
evanescent regime it amplifies by exactly `exp(sqrt(s^2-1) k0 d)`. Adding a
loss `Im ε = Im μ = γ` caps that amplification at a finite transverse
wavevector; the measured peak is `k_max d ≈ ln(1/γ) + (0.8–1.1)`, so the ideal
resolution is `~2 pi/k_max = 2 pi d / (k_max d)`. Energy conservation still
holds for the passive lossy slabs.

### Coherence: measured crossover

A free-standing `n = 1.5` slab, `r1 = -0.2`, `r2 = +0.2`, with the phase fixed
at `δ = 1.234` and round-trip optical path `2 n d`. The partial-coherence sum
weights each pair of round trips by `gamma = exp(-(path/L_c)²)`. Coherent Airy
and incoherent sums bracket the result:

| indicator `n d / L_c` | R |
|---|---|
| 10 (incoherent) | 0.076923 (= incoherent sum `2r²/(1+r²)`) |
| 2.5 | 0.076923 |
| 1.0 | 0.078065 |
| 0.5 | 0.099868 |
| 0.25 | 0.125724 |
| 0.05 (coherent) | 0.139262 (= coherent Airy) |

The coherent and incoherent limits are both exact checks: the `L_c -> 0` value
equals `2r²/(1+r²)` to \(10^{-6}\) and the `L_c -> ∞` value equals the Airy
reflectance to \(10^{-9}\). The crossover is at
`n d / L_c ~ 1`: below it the phase must be kept.

## Heat-source coupling

`Solution::power_density(z, I0)` returns W m⁻³, using
`p = (2 pi/lambda)(Im N² / q0) I0 |E|²` with `|E|` normalised to unit incident
amplitude. Integrating it across a lossy film recovers `I0 · A` to
\(2.3\times10^{-9}\) relative, and `layer_heat_source(j, I0)` gives the layer
average. `power_density_profile` samples it for a thermal grid. It is a
volumetric source for a split thermal step; it does not include thermal
emission back into the optics (that loop is `radiation`).

## Limits and honest misses

* The transfer matrix is exact only for planar, isotropic, static layers. It
  does not include roughness, anisotropy, magneto-optics or moving media.
* The quarter-wave stop-band formula is exact **in frequency**; the crates
  exposes both the frequency and the wavelength width, which differ by the
  Jacobian `1/(1-Delta^2/4)` (2.1% at `n_H/n_L = 1.586`).
* The effective-medium error indicator is an empirical quadratic envelope from
  the benchmark, not a proven bound; it can be exceeded for very high contrast
  or near a resonance of the constituents. The same caveat applies to
  `homogenisation_error_2d`.
* The incoherent closure and the partial-coherence sum assume lossless
  interfaces; the absorptive incoherent two-surface formula is provided
  separately (`incoherent_film_reflectance`).
* The Kretschmann SPR dip is 0.376° from the analytic surface-plasmon angle at
  50 nm gold on `n = 1.5` glass; the shift is the physical finite-film/loss
  effect, reported rather than tuned away.
* The plane-wave expansion is 2D only (rods or holes invariant along the axis),
  and uses the sharp-interface Fourier coefficients. The TM (inverse-rule) gap
  is converged at tens of plane waves; the TE inverse rule needs more plane
  waves for the same accuracy, which is why the complete-gap ratio is quoted at
  `n=6` and lands 0.5% high relative to the smoothed MPB result. No 3D
  crystal, no out-of-plane/slab modes, no nonlinear response.
* The negative-index model is a homogeneous effective medium (`ε(ω)`, `μ(ω)`
  with no spatial dispersion). The perfect-lens resolution law
  `k_max d ~ ln(1/γ)` is measured from the Lorentz model here, not derived; the
  branch selection is the signed one `q = N cos θ`, and a medium with only one
  of `Re ε`, `Re μ` negative has no negative-index band.
* No surface-plasmon propagation length and no waveguiding solver.

## References

* I. H. Malitson, "Interspecimen comparison of the refractive index of fused
  silica", JOSA **55**, 1205 (1965), doi:10.1364/JOSA.55.001205.
* M. Daimon and A. Masumura, "Measurement of the refractive index of distilled
  water...", Appl. Opt. **46**, 3811 (2007), doi:10.1364/AO.46.003811.
* A. D. Rakić, A. B. Djurišić, J. M. Elazar, M. L. Majewski, "Optical
  properties of metallic films for vertical-cavity optoelectronic devices",
  Appl. Opt. **37**, 5271 (1998), doi:10.1364/AO.37.005271.
* H. A. Macleod, *Thin-Film Optical Filters*, 4th ed. (transfer-matrix
  admittance conventions).
* A. Yariv and P. Yeh, *Optical Waves in Crystals* (Bloch waves, form
  birefringence).
* K. M. Ho, C. T. Chan, C. M. Soukoulis, "Existence of a photonic gap in
  periodic dielectric structures", Phys. Rev. Lett. **65**, 3152 (1990),
  doi:10.1103/PhysRevLett.65.3152 (plane-wave expansion and the inverse-`ε`
  rule).
* J. D. Joannopoulos, S. G. Johnson, J. N. Winn, R. D. Meade, *Photonic
  Crystals: Molding the Flow of Light*, 2nd ed., Princeton University Press
  (2008), Figs. 5.2 and 5.10 (square-rod 31.4% TM gap and triangular-hole
  18.6% complete gap).
* V. G. Veselago, "The electrodynamics of substances with simultaneously
  negative values of `ε` and `μ`", Sov. Phys. Usp. **10**, 509 (1968),
  doi:10.1070/PU1968v010n04ABEH003699.
* J. B. Pendry, A. J. Holden, D. J. Robbins, W. J. Stewart, "Magnetism from
  conductors and enhanced nonlinear phenomena", IEEE Trans. Microwave Theory
  Tech. **47**, 2075 (1999), doi:10.1109/22.798002 (split-ring `μ`).
* J. B. Pendry, "Negative refraction makes a perfect lens", Phys. Rev. Lett.
  **85**, 3966 (2000), doi:10.1103/PhysRevLett.85.3966 (perfect lens and the
  resolution limit).
* L. Li, "Use of Fourier series in the analysis of discontinuous periodic
  structures", J. Opt. Soc. Am. A **13**, 1870 (1996),
  doi:10.1364/JOSAA.13.001870 (inverse Fourier-factorisation rule).
