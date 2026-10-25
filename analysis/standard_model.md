# Standard Model CPU kinetic closures

This records the initial implementation and its baseline runs. Subsequent
[relativistic field and interaction improvements](standard_model_relativistic.md)
extend its capabilities and validation; use that report and the module README
for current behavior.

Implemented October 5, 2026. The `standard-model` crate introduces all elementary
fermion flavors, antiparticles and boson components in a usable relativistic
kinetic module, with explicit model domains and conservation bookkeeping. It
adds 37 species entries and selected decay/QED processes, rather than a complete
Standard Model field or event generator. See the
[usage and physical limits](../fluid/crates/standard-model/README.md).

## Implemented behavior

- On-shell GeV four-momenta; SI positions, times and fields. Massless species
  travel at c and massive ones at c p/E. Uniform-field impulses use relativistic
  Boris rotation, preserve magnetic energy and report electric field work.
- Proper-time exponential lifetimes; finite event-time decay cascades without
  boson-lifetime substeps. Daughters travel only after birth. Independent RNG
  streams and ordered reductions preserve results across thread counts.
- Measured inclusive W/Z/tau weights, dominant top bW, approximate pinned Higgs
  weights. Charged-lepton decays use finite daughter mass and tree-level V−A
  weights; other resolved two-body decays are unpolarized and isotropic.
- Scheduled free-lepton Klein–Nishina Compton events, comoving threshold e+/e−
  annihilation and externally parameterized elastic two-body kinematics. No
  automatic pair search, medium collision rates or generic strong/weak kernel.
- Ideal zero-chemical-potential quantum-gas integrals, exact massless limits,
  explicit degeneracy factors and SI conversions.
- Single physical photons can become existing radiation packets. Scheduled
  transferred four-momenta can become conservative SI fluid sources, with
  rejection of inadmissible updates before any cell commits.
- Vacuum rejects independent light partons; colored decay products are retained
  for joint hadronization. Partonic mode permits their noninteracting transport.
  Hadronic tau and off-shell Higgs decays remain unresolved events. Their full
  energy, momentum, charges, origin, time and position survive and can be handed
  to an external closure. They are not silently absorbed or discarded.

Confinement, showers, neutrino mass/mixing/scattering, spin correlations,
electroweak production, light-quark weak decay closures and interacting thermal
EOS tables remain unimplemented. PDG light-quark running masses are kinetic
parameters, not measured free-quark masses. Formal ideal-gas output below QCD
confinement or across electroweak restoration is not a physical plasma EOS.

## Verification

```sh
cd fluid
cargo test -p standard-model -p compressible -p radiation -p semiconductor --release --lib --tests
cargo test -p standard-model --release --doc
cargo clippy -p standard-model --all-targets --no-deps -- -D warnings
cargo fmt --check -p standard-model
RAYON_NUM_THREADS=8 cargo run -p standard-model --release --example particles -- 1000 partonic
RAYON_NUM_THREADS=8 cargo run -p standard-model --release --example particles -- 1000 vacuum
```

All **88 library/integration tests** passed: the previous 68 tests plus 19
kinetic physics tests and one fluid source test. The new suite independently
checks the Michel electron mean/second moment and both neutrino means, lifetime
statistics and Lorentz dilation, open phase space, daughter charges, conservation
through every supported decay, thread reproducibility, partitioned time steps,
magnetic energy, Compton recoil and differential/total cross sections, ideal-gas
limits and thermodynamic identity, source conversion, and atomic rejection.
The thread test uses 1221 injected particles, exercising parallel histories.
The README usage example is also tested as Rust documentation.

## Complete-species smoke run

One thousand copies of each allowed species, evolved for 100 microseconds:

| Domain | Initial | Live after flight | Unresolved records | Decays | Relative energy error | Momentum error / initial energy |
|---|---:|---:|---:|---:|---:|---:|
| Partonic | 37,000 | 49,533 | 1,964 | 13,086 | 3.428e-14 | 1.625e-14 |
| Vacuum | 19,000 | 21,542 | 11,981 | 13,120 | 1.038e-13 | 4.697e-14 |

Additive electric, baryon and flavor lepton charges match exactly, including the
retained ledger. One- and eight-thread partonic runs produced identical printed
physics results and event counts. Recorded single-run flight times were 25.6 ms
(one thread), 11.7 ms (eight threads) and 8.0 ms (vacuum, eight threads) on the
existing CPU power profile. These are smoke timings, not a repeat benchmark or
a general speedup claim; setup and thermal diagnostics are outside this timing.
Ballistic simulation histories avoid per-particle queues/allocations; only
histories with an actual decay allocate their cascade buffers.

The formal high-temperature ideal-gas check gives g*=106.749997. The exact
massless degeneracy count is 28 bosonic + (7/8) × 90 fermionic = 106.75.
Photon conversion returns photon_count=1; a 511 keV rest-target Compton encounter
has cross section 2.865399e-29 m². Raw stdout is retained in
[standard_model_results.json](standard_model_results.json). Data provenance is
pinned in the [module README](../fluid/crates/standard-model/README.md#data-provenance).
