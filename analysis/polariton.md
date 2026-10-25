# Phonon polaritons: where crystal vibrations meet light

Implemented in `fluid/crates/photonics/src/polariton.rs` with tests in
`fluid/crates/photonics/tests/polariton.rs`. This rung sits on top of the
transfer-matrix and dispersion models of `analysis/photonics.md` and below the
lattice dynamics of `analysis/materials.md`. The pattern is the one in
[`emergence_ladder.md`](emergence_ladder.md):

* **coarse model** — a local, frequency-only permittivity `eps(omega)` built
  from the zone-centre optical phonons;
* **validity indicator** — `q a` (spatial dispersion, wavevector against the
  Brillouin zone) and `a/d` (film thickness against the lattice spacing);
* **fold** — the pole of `eps(omega)` *is* the transverse-optical (TO) phonon
  the dynamical matrix predicts; the Born effective charge supplies the
  longitudinal-optical (LO) shift through the non-analytic (depolarisation)
  term.

Everything below is printed by `cd fluid && cargo test -p photonics`. The full
polariton suite runs in about 0.01 s.

## 1. The infrared dielectric function

The factorised (Gervais–Piriou) permittivity of one polar mode is

```
eps(omega) = eps_inf (omega_LO^2 - omega^2 - i gamma_LO omega)
                     / (omega_TO^2 - omega^2 - i gamma_TO omega).
```

With `gamma_LO = gamma_TO = gamma` this is exactly the single Lorentz
oscillator `eps_inf + eps_inf(omega_LO^2-omega_TO^2)/(omega_TO^2-omega^2-i gamma omega)`;
the test checks the two forms agree to `1e-12` at seven frequencies. Two
limits are exact:

* `eps(0) = eps_inf (omega_LO/omega_TO)^2` — the **Lyddane–Sachs–Teller (LST)**
  relation. For a single oscillator the code returns it as an identity, with
  zero residual; the test asserts `eps_static = eps_inf * (LO/TO)^2` to
  `1e-12`.
* `eps(omega -> infinity) = eps_inf`; measured to `1e-4` at
  `omega/2 pi c = 10^6 cm^-1`.

The denominator pole sits at `omega_TO` (Re eps changes sign across it) and the
numerator zero at `omega_LO` (lossless `|eps(omega_LO)| = 0`).

### Sourced parameters

| material | `omega_TO` (cm^-1) | `omega_LO` (cm^-1) | `eps_inf` | `eps0` (handbook) | source |
|---|---|---|---|---|---|
| 3C-SiC | 796.2 | 972.7 | 6.52 | 9.72 | Palik (1985); Spitzer et al. (1959) |
| NaCl | 164 | 264 | 2.34 | 5.90 | Palik II (1991); Born & Huang (1954) |
| GaAs | 268.3 | 292.3 | 10.9 | 12.9 | Blakemore (1982); Palik (1985) |
| hBN (in-plane) | 1367 | 1580 | 4.95 | 6.70 | Geick, Perry & Rupprecht (1966) |

The damping (`gamma ~ 2–10 cm^-1`) is a typical room-temperature linewidth and
is flagged as an estimate; all exact identities use a lossless copy.

### LST and static permittivity vs handbook

`eps0` from the model is `eps_inf (omega_LO/omega_TO)^2`, compared with the
independently measured handbook value:

| material | `eps0` model | `eps0` handbook | relative error |
|---|---|---|---|
| 3C-SiC | 9.7311 | 9.720 | +0.11 % |
| NaCl | 6.0637 | 5.900 | +2.77 % |
| GaAs | 12.9373 | 12.900 | +0.29 % |
| hBN (in-plane) | 6.6128 | 6.700 | −1.30 % |

The residuals are the scatter between independently measured `omega_LO`,
`omega_TO`, `eps_inf` and `eps0`, not model error: LST is exact for the form.
NaCl is the largest because its TO/LO are the broad, strongly damped classic
values.

## 2. Reststrahlen bands

Between `omega_TO` and `omega_LO` the permittivity is **negative**, so a bulk
crystal is a perfect mirror: at normal incidence
`R = |(1-sqrt(eps))/(1+sqrt(eps))|^2 -> 1`. The lossless test finds `R = 1` to
`1e-12` at five interior frequencies for all four materials, and `R(TO+) =
R(LO-) = 1.0000` for each. Just outside the band `R` drops sharply
(e.g. SiC: `R(0.98 TO) = 0.65`, `R(1.02 LO) = 0.007`). The band edges are
exactly the TO and LO frequencies.

**SiC vs the published band.** The task quotes SiC ≈ 790–970 cm^-1. The model
band is

```
lossless Re eps < 0 :  796.2 - 972.7 cm^-1   (errors +6.2 / +2.7 cm^-1)
damped half-power    :  748.2 - 971.6 cm^-1
damped peak R        :  0.9695 at 844.5 cm^-1
```

The edge definition matters: the *Re eps < 0* edges are TO and LO and agree
with the published 790–970 cm^-1 to a few wavenumbers; half-power points of the
damped spectrum extend below TO because the resonance has a long low-frequency
tail. This is reported, not tuned.

**Thin films** go through the existing transfer-matrix code
(`Stack`/`Layer`): a lossless 2 µm SiC film on an `n = 1.5` substrate has
`R = 0.9943` mid-band, a 20 µm film `R = 1.0000`, and far below the band
(`0.5 omega_TO`, where `eps ~ 10.8`) `R = 0.5704`. The 2 µm film is thinner
than the evanescent decay length, so a little light tunnels through; the
in-band mirror behaviour and the out-of-band contrast are the benchmarks.

## 3. Bulk and surface polaritons

**Bulk polariton.** `k = (omega/c) sqrt(eps(omega))`. For lossless SiC the two
propagating branches are separated by the Reststrahlen gap:

* lower branch: `k` real from `omega -> 0` up to the `omega_TO` asymptote;
* upper branch: `k -> 0` as `omega -> omega_LO+` and photon-like for
  `omega >> omega_LO`;
* between TO and LO, `eps < 0`, `k` is purely imaginary — no propagating
  transverse mode (the test asserts `Re k = 0`, `Im k > 0` at three interior
  frequencies).

**Surface phonon polariton (SPhP).** At a crystal/vacuum interface

```
k_sp = (omega/c) sqrt(eps/(eps+1)),   exists where eps < -1,
```

so the surface mode lives between `omega_TO` and the surface-mode frequency
`omega_s` where `eps(omega_s) = -1`. For the single lossless oscillator this is
closed-form:

```
omega_s^2 = (eps_inf omega_LO^2 + omega_TO^2)/(1 + eps_inf).
```

| material | `omega_s` (cm^-1) | interval |
|---|---|---|
| 3C-SiC | 951.1 | TO 796.2 – LO 972.7 |
| NaCl | 238.5 | 164 – 264 |
| GaAs | 290.4 | 268.3 – 292.3 |
| hBN (in-plane) | 1546.3 | 1367 – 1580 |

The test confirms `eps(omega_s) = -1` to `1e-9` for both the analytic formula
and an independent bisection of `Re eps = -1`, that `k_sp` diverges as
`omega -> omega_s^-`, that `k_sp` is real for `omega < omega_s`, and that it is
purely imaginary on `(omega_s, omega_LO)` (no surface mode there).

### Otto/Kretschmann ATR verification

The SPhP is detected as an attenuated-total-reflection dip in the Otto
configuration prism / air gap / bulk crystal, using the transfer-matrix code
(`Stack::solve_s`, p-polarisation). The transfer-matrix reflection itself is
used; no separate polariton solver is involved.

* **Angle at fixed frequency.** SiC, prism `n_p = 2.4`, gap 4 µm,
  `nu = 900 cm^-1`. The analytic phase-matching sine is
  `sin theta = sqrt(eps/(eps+1))/n_p = 0.4652` (`27.7244°`); the TMM dip is at
  `28.3931°`, a shift of **+0.67°**, with `R_min = 0.410`.
* **Frequency at fixed angle.** Prism `n_p = 2.4`, gap 6 µm,
  `theta = 30°` (`s = 1.2`). Solving `Re eps = s^2/(1-s^2)` gives
  `917.44 cm^-1`; the TMM reflectance dip is at `916.84 cm^-1`, an error of
  **−0.60 cm^-1** (`R_min = 0.183`).

A lossless crystal gives `R = 1` at every angle (all media lossless, no
absorption channel), so the dip only exists because of the crystal damping —
the same reason the gold surface-plasmon dip needs the metal's loss. For small
gaps the dip shifts strongly with gap thickness (a finite-cavity effect, up to
several degrees); the quoted operating points (4–6 µm) are where the coupling
is weak enough that the shift is sub-degree and sub-wavenumber. This is
reported rather than tuned away.

## 4. The fold to lattice dynamics

`DiatomicChain` is a 1D ionic crystal: masses `m1`, `m2` at `0` and `a/2`,
nearest-neighbour harmonic spring `k`, cell length `a`. Its 2×2 dynamical
matrix gives

```
omega_TO^2 = 2 k (1/m1 + 1/m2)     at q = 0,
omega^2(q)  = k(1/m1+1/m2) +/- k sqrt((1/m1-1/m2)^2 + 4 cos^2(qa/2)/(m1 m2)).
```

The acoustic sum rule `omega_acoustic(0) = 0` and the closed form of the
product `omega_ac^2 omega_op^2 = 4 k^2 sin^2(qa/2)/(m1 m2)` are checked (the
product to `1e-9` relative for `q > 0`).

**The Born charge turns the lattice mode into light.** With Born effective
charge `Z*`, ion-pair density `N` and reduced mass `mu`, the oscillator
strength is

```
S = Z*^2 e^2 N / (eps0 mu)   [rad^2/s^2],
```

so the long-wavelength response is `eps(omega) = eps_inf + S/(omega_TO^2 -
omega^2 - i gamma omega)`, whose **denominator pole is the chain's
`omega_TO`**, and whose static value and LO shift are

```
eps0   = eps_inf + S/omega_TO^2,
omega_LO^2 = omega_TO^2 + S/eps_inf     =>    eps0/eps_inf = (omega_LO/omega_TO)^2.
```

The test uses NaCl masses (`22.99` and `35.45` u) and a spring chosen so the
zone-centre optical mode is the measured `164 cm^-1`. The dynamical matrix
gives `omega_TO = 164.000 cm^-1`, and the permittivity pole built from it is
`164.000 cm^-1`: the electromagnetic closure reproduces the lattice mode
exactly, and LST holds by construction. Using `Z* = 1` and the NaCl pair
density `N = 4/a^3` (`a = 5.64 Å`) gives `omega_LO = 246.0 cm^-1` and
`eps0 = 5.265`, against the measured `264 cm^-1` and `5.90`. That miss is the
Born charge: the *effective* (transverse) charge of NaCl is larger than 1
(≈1.1–1.2, and the simple relation omits the local-field/Lorentz factor). The
fold's content is the pole matching and the LST identity, not a parameter-free
prediction of `Z*`; the miss is reported honestly.

### The indicator: when the local `eps(omega)` stops holding

The local permittivity assumes `q -> 0`. The chain's optical frequency is
`q`-dependent, so the relative departure from the zone-centre value is

```
spatial_dispersion_indicator(q, a) = |q a|,   relative shift ~ O((q a)^2).
```

Measured for the same chain:

| `q` | relative shift `|omega_op(q) - omega_TO|/omega_TO` |
|---|---|
| `0.005 pi/a` | 7.36e-6 |
| `0.01 pi/a` | 2.94e-5 |
| `pi/a` (zone boundary) | 2.21e-1 |

The factor between `0.01 pi/a` and `0.005 pi/a` is `4.000`, confirming the
quadratic approach to the continuum; at the zone boundary the shift is 22 %,
so the local `eps(omega)` closure is simply wrong there. The second indicator
is `film_thickness_indicator(a, d) = a/d`: a film of 10 lattice spacings is
homogeneous (`a/d = 0.1`), a film of 2 spacings is not (`a/d = 0.5`). For a
film only a few lattice spacings thick the discrete atomic layers and the
thickness-dependent image modes matter; the homogeneous-slab transfer matrix
cannot see them.

## Limits and honest misses

* **The local model is isotropic and homogeneous.** A single `eps(omega)`
  cannot describe anisotropy (hBN is uniaxial: only the in-plane response is
  modelled here), spatial dispersion, or the surface reconstruction at a real
  interface.
* **The factorised form with a single set of `gamma` is not
  Gervais–Piriou's full multi-mode fit.** Independent `gamma_LO`, `gamma_TO`
  per mode are supported and differ from the Lorentz equal-damping form; the
  materials table uses one mode each. Materials with many overlapping reststrahlen
  bands (α-quartz, perovskites) need several factorised modes and a numeric
  `omega_s`, which the crate does not yet provide.
* **The surface-mode formula `omega_s` is the lossless single-oscillator
  result.** With damping the `Re eps = -1` frequency shifts; the code provides
  `surface_mode_wavenumber_numeric` for that case.
* **ATR dip shifts with the gap.** The dip is verified to 0.7° in angle and
  0.6 cm^-1 in frequency at weak coupling; at small gaps the finite cavity
  shifts it by several degrees. This is the physics of the Otto geometry, not a
  numerical error, and the tests state the operating point.
* **The lattice-dynamics fold is a model.** The Born charge needed to match
  NaCl is not predicted; the demo uses `Z* = 1` and reports the resulting
  `omega_LO`/`eps0` miss. The fold establishes the *identification* of the
  dielectric pole with the zone-centre lattice mode and the LST identity.
* **Neither SPhP propagation length nor a full polariton band structure is
  computed.** The rung is dispersion relations and reflectance spectra, not a
  waveguide solver.
* **The `DiatomicChain` is a separate, minimal reimplementation** of the
  `materials` crate's `Harmonic`/`chain_lattice` construction (the photonics
  crate takes no dependency on `materials`, matching the repository's
  module-boundary convention). It reproduces the exact diatomic-chain
  dispersion identities.

## References

* E. D. Palik, *Handbook of Optical Constants of Solids* (Academic Press,
  1985) and Vol. II (1991): SiC, GaAs, NaCl lattice-band parameters.
* W. G. Spitzer, D. Kleinman, D. Walsh, "Infrared properties of hexagonal
  silicon carbide", Phys. Rev. **113**, 127 (1959),
  doi:10.1103/PhysRev.113.127.
* J. S. Blakemore, "Semiconducting and other major properties of gallium
  arsenide", J. Appl. Phys. **53**, R123 (1982),
  doi:10.1063/1.331665.
* R. Geick, C. H. Perry, G. Rupprecht, "Normal modes in hexagonal boron
  nitride", Phys. Rev. **146**, 543 (1966),
  doi:10.1103/PhysRev.146.543.
* M. Born and K. Huang, *Dynamical Theory of Crystal Lattices* (Oxford, 1954):
  the TO/LO splitting, the non-analytic term, and the LST relation.
* R. Lyddane, R. Sachs, E. Teller, "On the polar vibrations of alkali
  halides", Phys. Rev. **59**, 673 (1941),
  doi:10.1103/PhysRev.59.673.
* F. Gervais and B. Piriou, "Anharmonicity in several-polar-mode crystals:
  adjusting phonon self-energy of LO and TO modes in Al2O3 and TiO2 to fit
  infrared reflectivity", J. Phys. C **7**, 2374 (1974),
  doi:10.1088/0022-3719/7/13/017 (factorised form).
* H. Raether, *Surface Plasmons on Smooth and Rough Surfaces and on Gratings*
  (Springer, 1988); A. Otto, "Excitation of nonradiative surface plasma waves
  in silver by the method of frustrated total reflection", Z. Phys. **216**,
  398 (1968) — attenuated total reflection geometry (also used for surface
  phonon polaritons).
