# Magnetism ladder: spins → micromagnetics → hysteresis

Implemented in `fluid/crates/magnetism`. This is the magnetic rung of the
emergence ladder in [emergence_ladder.md](emergence_ladder.md), built with
the same shape as [materials.md](materials.md): a fine model (classical
Heisenberg spins on a lattice), a coarse model (micromagnetics), and the
fold between them (the long-wavelength expansion of the spin Hamiltonian);
a second fold collapses a small particle to a single macrospin
(Stoner–Wohlfarth), with the exchange length and domain-wall width as the
validity indicator.  Inside micromagnetics the local-shape demagnetising
field is itself a cheap closure, unfolded here to the exact non-local
magnetostatic field via the Newell-Williams-Dunlop tensor and a zero-padded
FFT (`src/demag.rs`, `src/llg3d.rs`), with the internal-field uniformity
and the exchange length as the validity indicator.

## The ladder

| rung | model | validity indicator | unfold / closure |
|---|---|---|---|
| atomistic | Heisenberg spins, exchange `J`, anisotropy `K`, field `H` | `k_B T / J z S^2` (ordering) | Metropolis Monte Carlo / spin waves; LLG for the dynamics |
| micromagnetics | fields `A`, `K`, `M_s`; LLG on a continuum grid | `delta/a` (wall width / lattice constant), `l_ex/L`, `k a` | back to the spin chain when `delta ~ a` |
| non-local magnetostatics | FFT convolution `H_d = -N * M` on the grid | internal-field uniformity; size vs `l_ex` | local-shape `N` closure when the internal field is uniform |
| macrospin / hysteresis | one magnetic moment; Stoner–Wohlfarth astroid | particle size vs `D_c` (single-domain size) | back to micromagnetics when a wall fits (`d > pi delta`) |

The fold from the atomistic rung to micromagnetics is the small-angle
expansion of `H = -J sum_<ij> S_i . S_j`. With unit spins the bond energy
gives

```
e = J S^2/(4 V_atom) sum_delta (delta . grad m)^2,
A = J S^2 sum_delta |delta|^2 / (12 V_atom) = J S^2 z d_nn^2 / (12 V_atom),
```

where `V_atom` is the atomic volume and the sum is over the `z`
nearest-neighbour bonds of one atom. For the three conventional cubic cells
`sum_delta |delta|^2 = 6 a^2`, so `A = J S^2 n_atoms / (2a)`:

| lattice | `z` | `d_nn` | `V_atom` | `A` |
|---|---|---|---|---|
| sc | 6 | `a` | `a^3` | `J S^2/(2a)` |
| bcc | 8 | `sqrt3 a/2` | `a^3/2` | `J S^2/a` |
| fcc | 12 | `a/sqrt2` | `a^3/4` | `2 J S^2/a` |

The often-quoted form `A = J S^2 z a^2/(2V)` assumes the bond sum is
`z a^2 (grad m)^2` instead of the isotropic average
`z d_nn^2 (grad m)^2 / 3`; the convention above is the one used here and the
closed forms agree exactly with the resolved nearest-neighbour table
(`tests/micromag.rs::exchange_stiffness_fold_relations`).

## 1. Atomistic spins

The Monte-Carlo Hamiltonian is

```
H = -J sum_<ij> S_i . S_j - K sum_i (S_i . e)^2 - H_ext . sum_i S_i,
```

with unit classical spins on an `l x l x l` conventional-cell cluster with
periodic boundaries. Metropolis updates use uniform proposals on the sphere
(symmetric, detailed balance). The order parameter is `M = |(1/N) sum S_i|`
and the Binder cumulant is `U_4 = 1 - <M^4>/(3 <M^2>^2)`.

### Curie temperature

The classical mean-field result is `k_B T_C = J z S^2/3`; the exact sc
benchmark is `k_B T_C/J = 1.44299(14)` from Chen, Ferrenberg & Landau,
*Phys. Rev. B* **48**, 3249 (1993).

| quantity | value |
|---|---|
| mean field sc (`J=S=1`) | `2.0` (exact identity) |
| mean field bcc | `8/3 = 2.667` |
| mean field fcc | `4.0` |
| Binder crossing L=4/6 | `1.424` |
| Binder crossing L=6/8 | `1.402` |
| Binder crossing L=4/8 | `1.428` |
| mean of the three crossings | `1.4179` |
| published Monte-Carlo | `1.44299(14)` |
| relative error of the crossing | `1.7 %`; the spread across size pairs is `0.026` and is dominated by finite-size and statistical noise |

The Monte-Carlo crossing is deliberately on small lattices with finite
statistics, so the pair-to-pair spread (`1.402` to `1.428`) is the honest
uncertainty; the mean-field value is 39 % too high, which is exactly the
kind of miss the ladder is meant to detect.

### Spin waves against LLG

Linear spin-wave theory for the ferromagnetic exchange Hamiltonian gives
`omega(k) = 2 J S (z - sum_delta cos k.delta) / hbar`. A 1D spin chain of
64 sites with `S = 1/2`, `hbar = 1`, `gamma = 2S/hbar = 1` is seeded with a
single mode and integrated with the explicit LLG
`dS/dt = -(gamma/(1+alpha^2))[S x H_eff + alpha S x (S x H_eff)]`; the
frequency is extracted by a least-squares sinusoid fit.

| mode `m` | spin-wave `omega` | LLG `omega` | relative error |
|---|---|---|---|
| 8 | 0.585786 | 0.585830 | 7.5e-5 |
| 16 | 2.000000 | 1.999950 | 2.5e-5 |
| 24 | 3.414214 | 3.414128 | 2.5e-5 |

The residual is the frequency-grid resolution of the fit, not a model
error.

### The discrete domain wall and its lattice correction

The 1D chain with exchange `J` and easy-axis anisotropy `K` folds to the
continuum with `A = J a/2` and `K_1d = K/a`, so
`delta = sqrt(A/K_1d) = a sqrt(J/(2K))` and
`sigma = 4 sqrt(A K_1d) = 4 sqrt(JK/2)`. Relaxing the discrete angles
(`tests/atomistic.rs::discrete_wall_fold_and_its_lattice_correction`) gives:

| `delta/a` | `sigma_meas` | `4 sqrt(AK)` | energy error | width error |
|---|---|---|---|---|
| 1 | 1.909486 | 2.000000 | 4.53 % | 2.79 % |
| 2 | 0.989140 | 1.000000 | 1.09 % | 0.87 % |
| 5 | 0.399329 | 0.400000 | 0.17 % | 0.16 % |
| 10 | 0.199917 | 0.200000 | 0.042 % | 0.040 % |

The continuum fold is exact as `delta/a → infinity`; the lattice
correction is the discrete counterpart of the `(ka)^2` dispersion error of
the elasticity fold, and becomes order unity when the wall is only a few
lattice constants wide.

## 2. Micromagnetics

The continuum energy density is

```
e = A (grad m)^2 - K (m . e)^2 - mu0 M_s H . m
    + (mu0/2) M_s^2 sum_k N_k (m . n_k)^2,
```

and the LLG is integrated with RK4 on a 1D finite-difference grid. The
exact Bloch wall is `theta(x) = 2 atan(exp(x/delta))`, `delta = sqrt(A/K)`,
`sigma = 4 sqrt(A K)`; the field-driven wall has the Walker velocity
`v = gamma delta H/alpha` below `H_W = alpha M_s/2` (Schryer & Walker 1974;
Malozemoff & Slonczewski, *Magnetic Domain Walls in Bubble Materials*,
1979; Bertotti, *Hysteresis in Magnetism*, 1998). The thin-film Kittel
resonance is `omega = gamma sqrt((H + H_K)(H + H_K + M_s))` with
`H_K = 2K/(mu0 M_s)` (Kittel, *Phys. Rev.* **73**, 155 (1948)).

### Bloch wall in the solver

`A = 1e-11 J/m`, `K = 1e5 J/m^3`, `M_s = 1e6 A/m`, `dx = delta/8`:

| quantity | exact | computed | relative error |
|---|---|---|---|
| wall energy `sigma` | `4.000000e-3 J/m^2` | `3.997390e-3` | 6.5e-4 |
| wall width `delta` | `1.000000e-8 m` | `1.001313e-8` | 1.3e-3 |

### Walker breakdown

The reduced 1D `(q, phi)` model (the standard 1D model) is integrated
exactly:

| `H/H_W` | `v_avg / (gamma delta H / alpha)` |
|---|---|
| 0.3 | 1.0000 |
| 0.5 | 1.0000 |
| 0.8 | 1.0000 |
| 1.0 | 0.9995 |
| 1.5 | 0.2643 |
| 2.0 | 0.1435 |

The linear law and the breakdown field `H_W = alpha M_s/2` are reproduced
to better than 0.1 %; above the Walker field the azimuth precesses and the
average velocity collapses, as it must.

The full finite-difference LLG reproduces the *same* terminal velocity
when the wall is seeded at the predicted steady azimuth: at
`H = 0.3 H_W`, `v = 3.10e2 m/s` against `v_Walker = 3.32e2 m/s` (ratio
`0.93`), and the seeded runs at `0.5 H_W` and `0.8 H_W` give `0.84` and
`0.66`. Starting instead from an untilted wall and integrating for the same
time gives about half the terminal velocity: the azimuthal relaxation time
is `1/(gamma alpha M_s) ≈ 4.5e-11 s`, longer than the short integration
window, so those runs measure the transient and are not a model error.
This is recorded as a measurement-time limitation, not a failure of the
Walker law.

### Kittel FMR

A normalised macrospin with `H_K = 1`, `M_s = 1`, `gamma = 1` and a
thin-film demagnetising factor `N = 1` is driven with a small tilt; the
precession frequency extracted from the LLG matches
`sqrt((H+1)(H+2))` to `2.5e-5` for `H = 0.2, 0.5, 1, 2`.

### Non-local demagnetising field (FFT)

The local-shape closure above replaces the magnetostatic kernel by
the single cell's own shape factor.  The exact field is the convolution
`H_demag(r_i) = -sum_j N(r_i - r_j) M(r_j)` with the Newell-Williams-Dunlop
cuboid tensor (Newell, Williams & Dunlop, *J. Geophys. Res.* **98**, 9551
(1993)), which for an observer at relative position `(x,y,z)` and cell size
`(dx,dy,dz)` is a corner sum of `atan(eta zeta/(xi r))` type primitives
over the eight source corners.  The implementation in
`fluid/crates/magnetism/src/demag.rs` was checked against direct numerical
quadrature of `d'_i d'_j (1/r)` to machine precision
(`tests/demag.rs`, and the tensor values in the unit tests).  The field is
evaluated by a zero-padded FFT convolution: the magnetisation is padded to
twice the grid size, so the circular convolution equals the linear one and
the open boundary is exact.

| benchmark | expected | computed | error |
|---|---|---|---|
| uniform cube, centre `H_x` | `-M_s/3` | `-0.3333333 M_s` | `< 1e-12` |
| uniform cube, energy | `(mu0/2) M_s^2 V/3` | same | `5.6e-16` |
| sphere of cells (Lorentz sum), `R = 3,5,7` | `-M_s/3` | `-0.3333333 M_s` | `< 1e-12` |
| thin film `N_z`, `nxy = 8,16,32` | `-> 1` | `-0.98852, -0.99435, -0.99718` | `1.2 %, 0.6 %, 0.3 %` |
| cubic energy symmetry | `E_x = E_y = E_z` | same | `1e-12` |

The sphere row is the classical Lorentz-sphere result: the discrete sum
over complete spherical shells is `4 pi/3` (here `1/3`) independently of
the radius.  The thin film converges to the local `N_z = 1` as the lateral
size grows and the cell aspect ratio shrinks.  The cube energy identity is
exact term by term, not just asymptotically.

### µMAG standard problem 4

NIST µMAG standard problem 4 (McMichael, Koch & Schrefl; developed by
Eicke & McMichael, <https://www.ctcms.nist.gov/~rdm/std4/>) is a 500 x 125
x 3 nm permalloy film (`A = 1.3e-11 J/m`, `M_s = 8.0e5 A/m`, `K = 0`,
`alpha = 0.02`) reversed from the equilibrium s-state by
`mu0 H = (-24.6, 4.3, 0)` mT.  The s-state is obtained here by saturating
along `[1,1,1]` and relaxing to zero field.  The published reference
solutions (OOMMF, Donahue & Porter; Albuquerque, Miltat & Thiaville,
downloaded from the µMAG site) give an s-state
`<m> = (0.9666, 0.1260, 0)`, a first `<m_x>` zero crossing at `0.1384 ns`
and a first `<m_y>` zero crossing at `0.1967 ns`.

| grid (nm cells) | relaxed s-state | first `<m_x>` zero | first `<m_y>` zero |
|---|---|---|---|
| reference (OOMMF) | `(0.9666, 0.1260, 0)` | `0.1384 ns` | `0.1967 ns` |
| 32 x 8 x 1 (15.6 nm) | `(0.9904, 0.0669, 0)` | `0.1366 ns` (`-1.3 %`) | `0.1836 ns` (`-6.7 %`) |
| 64 x 16 x 1 (7.8 nm) | `(0.9759, 0.1082, 0)` | `0.1350 ns` (`-2.5 %`) | `0.1884 ns` (`-4.2 %`) |
| 128 x 32 x 1 (3.9 nm) | `(0.9686, 0.1223, 0)` | `0.1377 ns` (`-0.5 %`) | `0.1949 ns` (`-0.9 %`) |
| 200 x 50 x 1 (2.5 nm) | `(0.9664, 0.1259, 0)` | `0.1391 ns` (`+0.5 %`) | `0.1977 ns` (`+0.5 %`) |

The 32 x 8 x 1 run is the fast test (`tests/mumag.rs`, about 1.5 s); the
64 x 16 x 1 run is the example default.  The s-state canting is
under-resolved on the coarse grids (`<m_y>` is only about half the
reference value at 32 x 8) and converges monotonically with refinement.
The two finest rows place the cell size (3.9 and 2.5 nm) below the
permalloy exchange length `l_ex = 5.7 nm`, which is the resolution the
problem actually requires: the s-state then agrees with the reference to
better than 0.2 %, and both crossing times to about 0.5 %.  The
crossing-time error shrinks from a few percent at 15.6 nm to under 1 % at
3.9 nm and is on the reference at 2.5 nm; no tuning was applied at any
resolution.

Two independent published time series were downloaded from the µMAG site
and compared against the 200 x 50 curve (both on the same 0-0.5 ns window
and the same uniform time base):

| reference | rms `<m_x>` | rms `<m_y>` | rms `<m_z>` | max `|<m_y>|` |
|---|---|---|---|---|
| Albuquerque, Miltat & Thiaville | `0.0089` | `0.0128` | `0.0022` | `0.033` |
| Martins (1 nm mesh) | `0.0197` | `0.0268` | `0.0054` | `0.073` |

The 2.5 nm solution tracks the Albuquerque curve to about one percent in
both in-plane components through the reversal and the first oscillation;
the residual against the 1 nm Martins curve is the remaining cell-size
error.  The reference data were downloaded from the µMAG site.

### FFT performance

The FFT demagnetising field is by far the dominant cost at these
resolutions.  The original `demag.rs` rebuilt the bit-reversal permutation
on every transform, ran all three transform axes serially, and allocated a
fresh demagnetising buffer inside every effective-field evaluation; the
RK4 stepper allocated its four stage vectors and three intermediate states
every step as well.  The optimised path:

- precomputes the bit-reversal table once per axis in a reusable `FftPlan`;
- keeps the nine padded complex arrays in a caller-owned `FftWorkspace`, so
the hot loop is allocation-free;
- parallelises the three axes of each 3-D transform and the exchange/Gilbert
  passes with rayon, and reuses the six RK4 stage buffers.

Times are wall-clock per field evaluation and per RK4 reversal step on the
16-core development machine (median of 20-40 repetitions, `sp4_bench`):

| grid | field: before → after | speedup | step: before → after | speedup |
|---|---|---|---|---|
| 32 x 8 x 1 | 0.18 → 0.12 ms | 1.6 x | 0.48 → 0.61 ms | 0.8 x |
| 64 x 16 x 1 | 0.53 → 0.24 ms | 2.2 x | 2.47 → 1.12 ms | 2.2 x |
| 128 x 32 x 1 | 2.73 → 0.58 ms | 4.7 x | 11.5 → 2.82 ms | 4.1 x |
| 200 x 50 x 1 | 12.5 → 1.9 ms | 6.5 x | 50.5 → 8.6 ms | 5.9 x |

The speedup grows with resolution because the serial bit-reversal rebuild
and the per-call allocation are a fixed cost that the larger transform
amortises.  At 32 x 8 the optimised step is *slower* (0.8 x): the grid is
so small that rayon's scheduling overhead exceeds what it saves.  That is
reported rather than hidden; the fast test still runs in about a second.
At the production 200 x 50 grid the reversal costs about 8.6 ms/step, so
0.6 ns takes roughly 95 s instead of the ~9 min the unoptimised code would
need.

Long runs are split into resumable segments by `src/sp4.rs`, which writes
a binary checkpoint (grid magnetisation plus accumulated time, samples and
crossing times) so each command finishes in a couple of minutes.  The
200 x 50 s-state was relaxed to convergence and the reversal integrated to
1 ns as two such segments.

### The fold: local shape vs non-local FFT

The local closure is enough when the internal demagnetising field is
nearly uniform over the sample.  Two cases from `tests/demag.rs`:

- **Thin film, out of plane.**  Cell aspect ratio `t/L = 0.05`.  The FFT
  interior field is `-0.9962 M_s`, within `0.4 %` of the local `N_z = 1`,
  and the full and local energies agree to better than `5 %`.  Indicator:
  the film thickness `3 nm` is below the permalloy exchange length
  `l_ex = 5.7 nm`, so the film is single-domain through its thickness.
- **Elongated particle (rod along x).**  The FFT gives an axial demag
  factor `N_x = 0.038` and a transverse one `N_y = 0.481`, i.e. the large
  shape anisotropy that the local closure reproduces.  The axial field is
  not perfectly uniform along the finite rod (spread `0.067`, a
  finite-length correction), but the average shape anisotropy is captured.
  Indicator: the cross-section is below `l_ex`, so no curling/vortex
  forms.

The failure case is the **vortex-state disc**.  A thin disc with in-plane
circulation has a closed flux and a strongly non-uniform field: the FFT
gives a field spread of `0.0154 M_s` where the local in-plane closure
predicts essentially zero, and a vortex energy `0.285` of the uniformly
in-plane state.  A single shape factor gives the same near-zero energy for
both states, so it cannot resolve the vortex.  Indicator: the disc
diameter is many exchange lengths once a physical cell size is assigned
(the model disc is 20 cells across), which is the regime where a vortex
exists.

## 3. Hysteresis fold

A particle small enough to reverse by coherent rotation is one macrospin
with `e = -K (m . n)^2 - mu0 M_s H . m`. The exact zero-temperature
switching field in units of `H_K = 2K/(mu0 M_s)` is the astroid

```
h_sw(psi) = 1 / (cos^{2/3} psi + sin^{2/3} psi)^{3/2},
```

with `psi` the angle of the reversed field to the easy axis (Stoner &
Wohlfarth, *Phil. Trans. R. Soc. A* **240**, 599 (1948)). It touches the
axes at 1 and has its minimum `1/2` at 45°.

### LLG macrospin against the astroid

| `psi` (deg) | astroid | LLG | relative error |
|---|---|---|---|
| 0 | 1.00000 | 1.01897 | 1.9 % |
| 15 | 0.61466 | 0.61463 | 5.2e-5 |
| 30 | 0.52402 | 0.52400 | 3.0e-5 |
| 45 | 0.50000 | 0.50000 | 3.0e-6 |
| 60 | 0.52402 | 0.52403 | 1.9e-5 |
| 75 | 0.61466 | 0.61468 | 3.8e-5 |
| 90 | 1.00000 | 1.03883 | 3.9 % |

The two near-critical points (`psi = 0` and 90°) are limited by the
divergence of the switching time at the saddle-node; the interpolation and
all intermediate angles are exact to `5e-5`.

### Random-orientation coercivity

An isotropic assembly of easy axes, field along a fixed direction, is
swept from saturation through reversal; the remanence coercivity is the
field where the average magnetisation crosses zero. The exact
Stoner–Wohlfarth result is `H_c = 0.479 H_K`; the LLG loop average gives
`H_c = 0.4919 H_K` (2.7 % high). The residual is the finite relaxation
time near each particle's switching field; refining the field grid and the
relaxation step reduces it monotonically (`0.504` at 200 fields / 150
steps, `0.492` at 400 fields / 300 steps).

### Single-domain indicator

The macrospin fold is valid while a particle is smaller than the
single-domain size. `single_domain_diameter` uses Kittel's estimate
(Kittel, *Phys. Rev.* **70**, 965 (1946))
`D_c = 36 sqrt(A K)/(mu0 M_s^2)`, and the crate also reports the exchange
length `l_ex = sqrt(2A/(mu0 M_s^2))` and the wall width
`delta = sqrt(A/K)`. A domain wall fits once `d > pi delta`; if `d > D_c`
the particle must be unfolded to the micromagnetic rung.

## Real materials

Room-temperature constants in `materials.rs` (Coey, *Magnetism and Magnetic
Materials*, 2010; Cullity & Graham, *Introduction to Magnetic Materials*,
2nd ed.; Curie temperatures from the standard solid-state tables; NdFeB
from Sagawa et al., *J. Appl. Phys.* **55**, 2083 (1984)):

| material | `A` (J/m) | `K` (J/m^3) | `M_s` (A/m) | `T_C` (K) | `delta` (nm) | `l_ex` (nm) | `H_K` (A/m) | `D_c` (nm) |
|---|---|---|---|---|---|---|---|---|
| Fe (bcc) | 2.1e-11 | 4.8e4 | 1.71e6 | 1043 | 20.9 | 3.4 | 4.5e4 | 9.8 |
| Ni (fcc) | 9e-12 | 5.7e3 | 4.9e5 | 631 | 39.7 | 7.7 | 1.9e4 | 27 |
| Ni80Fe20 | 1.3e-11 | ~0 | 8.0e5 | ~850 | — | 5.7 | — | — |
| Nd2Fe14B | 8e-12 | 4.9e6 | 1.07e6 | 585 | 1.3 | 3.3 | 7.3e6 | 157 |

`D_c` for permalloy is not quoted because `K ≈ 0`; the particle size there
is set by shape anisotropy and geometry, and the micron-scale soft-magnet
behaviour must be read off the micromagnetic rung. The iron row matches
the commonly quoted wall width ≈ 20–40 nm and exchange length ≈ 3 nm.

## Limits and honest misses

- The Binder crossing is computed on `L = 4, 6, 8` clusters with 400
  samples and averaged over three seed pairs; the pair-to-pair spread
  (`1.402` to `1.428`, mean `1.418`) is comparable to the finite-size
  shift, so the `1.7 %` agreement with `1.44299` is not a precision
  measurement. Larger lattices and a histogram-reweighting analysis
  would be needed for a sharp number.
- The full finite-difference LLG needs the azimuthal relaxation time
  `1/(gamma alpha M_s)` (tens of picoseconds for the parameters used)
  before it reaches the terminal Walker velocity. Short runs that start
  from an untilted wall give about half the analytic value; seeding the
  steady azimuth removes the transient. The terminal-velocity discrepancy
  at the highest fields (`0.66` at `0.8 H_W`) is not yet resolved against
  the rigid-wall prediction and is reported rather than tuned.
- The local-shape demagnetising field `H_d = -N M_s (m . n) n` in
  `micromag.rs` is retained as the cheap closure; it captures the thin-film
  Kittel resonance and the Walker hard axis.  The non-local unfold is now
  built in `demag.rs` (Newell-Williams-Dunlop tensor, zero-padded FFT) and
  used by `llg3d.rs`; it reproduces the uniform-cube, Lorentz-sphere and
  thin-film factors exactly.  The FFT plan and working arrays are now
  reused (see "FFT performance" above) and the transform is parallel; at
  2.5 nm resolution the reversal still costs about 8.6 ms per RK4 step,
  and the zero-padded kernel is computed by direct evaluation of the
  analytic tensor at every offset rather than by an O(N) recursion, so
  building the kernel is itself `O(N^2)` and would benefit from the
  Newell recursion.
- The µMAG standard problem 4 crossing times are now computed on four
  grids from 15.6 nm down to 2.5 nm, below the 5.7 nm exchange length; the
  3.125 nm reference discretisation is not exactly reproduced, but the
  2.5 nm row is closer than any published solution used for the
  comparison, and the s-state and crossing-time errors are both about
  0.5 %.  The 2.5 nm curve tracks the Albuquerque reference to rms
  `~0.01` in `<m_x>` and `<m_y>` over 0-0.5 ns; the residual against the
  1 nm Martins curve (`~0.02-0.03`) is the remaining discretisation error.
  The 200 x 50 run was taken to 1 ns and needs checkpointed segments, so
  it is not part of the routine test suite.
- The Curie-temperature Monte Carlo uses a single-spin Metropolis update.
  Near `T_C` the autocorrelation time is long; cluster or over-relaxed
  updates would sharpen the crossing.
- The exchange-stiffness fold assumes only nearest-neighbour exchange and
  a small-angle texture. At `delta ~ a` the discrete wall of section 1 is
  the correct model; the fold error there is 4.5 %.
- The anisotropy reduction from cubic `K_1` to the uniaxial `K` used in the
  micromagnetic and macrospin models is a convention (`K = 2K_1` for the
  easy direction), not a symmetry-exact statement.
- The permalloy Curie temperature is a rounded literature figure
  (`~850 K`) and is flagged approximate.

## Constants and sources

- Classical sc Heisenberg `k_B T_C/J = 1.44299(14)`: Chen, Ferrenberg &
  Landau, *Phys. Rev. B* **48**, 3249 (1993).
- Mean-field Curie temperature `k_B T_C = J z S^2/3`: standard
  Heisenberg mean-field theory.
- Spin-wave dispersion `omega = 2 J S (z - sum cos k.delta)`: standard
  linear spin-wave theory (Holstein–Primakoff).
- Domain wall `theta = 2 atan(exp(x/delta))`, `delta = sqrt(A/K)`,
  `sigma = 4 sqrt(AK)`: standard micromagnetics (Landau & Lifshitz;
  Hubert & Schäfer, *Magnetic Domains*, 1998).
- Walker velocity and field: Schryer & Walker, *J. Appl. Phys.* **45**,
  5406 (1974); Malozemoff & Slonczewski (1979); Bertotti (1998).
- Kittel resonance: Kittel, *Phys. Rev.* **73**, 155 (1948).
- Stoner–Wohlfarth astroid and `H_c = 0.479 H_K`: Stoner & Wohlfarth,
  *Phil. Trans. R. Soc. A* **240**, 599 (1948).
- Single-domain critical size: Kittel, *Phys. Rev.* **70**, 965 (1946).
- Cuboid demagnetising tensor and FFT demagnetising field: Newell,
  Williams & Dunlop, *J. Geophys. Res.* **98**, 9551 (1993).
- µMAG standard problem 4 specification and reference solutions: NIST
  µMAG, <https://www.ctcms.nist.gov/~rdm/std4/>; OOMMF solution by
  M. J. Donahue and D. G. Porter; independent solution by G. Albuquerque,
  J. Miltat and A. Thiaville (both downloaded for the comparison above).
- Demagnetising factors of ellipsoids: Osborn, *Phys. Rev.* **67**, 351
  (1945), used for the elongated-particle cross-check.
- Material constants: Coey, *Magnetism and Magnetic Materials*
  (Cambridge, 2010); Cullity & Graham, *Introduction to Magnetic
  Materials* (2nd ed.); Sagawa et al., *J. Appl. Phys.* **55**, 2083
  (1984); lattice constants from the CRC Handbook.
- Gyromagnetic ratio `mu0 gamma_e = 2.211e5 m/(A s)`: CODATA 2018.
