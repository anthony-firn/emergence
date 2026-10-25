# Relativistic particle physics improvements

Static gravitational metrics and particle geodesics were added after this
report. See [gravity and solar validation](gravity_solar.md); the results and
counts below describe the earlier flat-spacetime implementation.

Implemented October 5, 2026, extending the
[initial Standard Model module](standard_model.md). The new APIs and physical
limits are documented in the [crate README](../fluid/crates/standard-model/README.md).

The particle solver integrates constant-field relativistic worldlines, momentum
and proper clocks together. For mixed fields it computes exp(K tau) and the
integral of that exponential. Pure E and pure B use analytic light-cone and
gyro solutions. A safeguarded proper-time root finds a requested laboratory
endpoint or the sampled decay event; daughters then move in the field for the
remaining interval. This implements the proper-time constant-field formulation
described by [Chin](https://arxiv.org/html/0809.0859). Prescribed-field work and
impulse appear in a separate four-momentum ledger. Magnetic-only laboratory
endpoints use constant gamma directly, avoiding bracketing and gyro substeps.

Known species masses anchor light-cone Lorentz boosts and invariant pair energies
without subtracting nearly equal E and |p|. Parallel fields survive boosts
without subtracting large longitudinal terms. Small-angle invariants use the
transverse cross product, avoiding false angular separation from vector
normalization roundoff. Compton sampling preserves energy fraction and small
angle separately, so forward recoil survives cos(theta) rounding to one.

Scheduled moving e-/e+ annihilation and real-photon pair production use the
finite-mass, unpolarized Born QED angular law, with electron, muon and tau
production thresholds. Total annihilation rates agree with the independent
[Geant4/Heitler reference](https://geant4.web.cern.ch/documentation/pipelines/master/prm_html/PhysicsReferenceManual/electromagnetic/electron_incident/annil.html).
The electron-channel cross sections satisfy sigma_pair = 2 beta_CM² sigma_ann,
including identical final photons. Their threshold flux product is pi r_e² c.
The total Breit–Wheeler law and differential law also agree with the
[WarpX reference](https://warpx.readthedocs.io/en/26.09/_static/doxyhtml/_linear_breit_wheeler_cross_section_8_h_source.html).
An inverse hyperbolic angular proposal resolves the forward/backward peaks
efficiently; unsupported peaks below the declared angular resolution are rejected.

All three neutrino/antineutrino flavors have a selected elastic electron
scattering channel using the finite-recoil chiral Fermi law. Tree coefficients
include electron-flavor charged current and left/right interchange for anti.
G_F is pinned; sin²(theta_W) and its scheme are caller inputs. Calibrated
effective coefficients can instead be supplied. The recoil interference term,
endpoint and integrated rate follow the leading-order formulas of
[Tomalak and Hill](https://arxiv.org/html/1907.03379). Their radiative correction
calculation is not implemented. The kernel rejects 2 m_e E_nu/M_W² > 1e-3
in the target rest frame and supplies only this elastic channel.

## Validation

From `fluid/`:

```sh
cargo test -p standard-model -p compressible -p radiation -p semiconductor --release --lib --tests
cargo test -p standard-model --release --doc
cargo clippy -p standard-model --all-targets --no-deps -- -D warnings
cargo fmt --check -p standard-model
cargo run -p standard-model --release --example relativistic
```

The suites contain 105 library/integration tests, including 36 Standard Model
tests (19 original, nine relativity and eight interaction tests). Validation
covers Lorentz-transformed fields, momenta and spacetime displacement;
composition/partitioning of trajectories; 1,000 gyro turns in one update;
deceleration from gamma=1e8 to rest; axial boost roundtrips and parallel-field
invariance through gamma=1e12; accelerated proper-time decay and its daughter
birth events; conservation including external fields; analytic cross sections,
detailed balance, sampled angular moments and recoil spectra; and atomic event
rejection. Compton numerical stress tests use target-rest E/m from 1e8 through
1e100, checking the transfer invariant as well as four-momentum conservation.
Those extreme tests concern arithmetic within the Born model, not its physical
validity at arbitrarily high energy. All 105 library/integration tests and the
README doctest passed. Clippy with warnings denied, formatting and diff checks
also passed.

The runnable example prints:

| Quantity | Result |
|---|---:|
| Recovered electron rest energy at axial gamma=1e12 | 5.1099895e-4 GeV |
| Moving electron pair invariant mass, same velocity | 1.0219979e-3 GeV |
| Accelerated muon sampled proper decay time | 2.529252540e-6 s |
| Accelerated muon laboratory decay time | 4.586622949e-6 s |
| Muon daughter birth position along x | 1.047218596e3 m |
| Field-inclusive relative energy ledger error | 1.175e-16 |
| Field-inclusive momentum error / final energy | 3.525e-16 |
| Annihilation sigma at beta_CM=0.8 | 1.233350838e-29 m² |
| Pair production sigma at the same invariant | 1.578689073e-29 m² |
| Threshold annihilation sigma times CM flux speed | 7.478838581e-21 m³/s |
| Tree nu_e/e elastic sigma at E=1 MeV, sin²theta_W=0.23112 | 7.169164202e-49 m² |

Raw deterministic output is retained in
[standard_model_relativistic_results.json](standard_model_relativistic_results.json).
No general speedup claim is made: magnetic motion removes gyro substeps for
constant fields, but mixed-field exponentials and endpoint solves cost more
than ballistic flight. The original particle-history parallelism is retained.

These are flat-spacetime test-particle and selected collision closures. Spatially
varying self-consistent fields, radiation reaction/emission, polarization and
spin correlations, confinement/hadronization, showers, neutrino mass/oscillations,
inelastic weak channels, gravitational geodesics and general relativity are
still absent. Caller-scheduled encounters require an external rate and spatial
eligibility closure. The ordinary fluid/material solvers remain classical.
