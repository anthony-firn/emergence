# Source provenance

The initial simulation/CAD snapshot was extracted from local commit
`cdd946b` of the Navier–Stokes working repository. This public repository starts
a new history; it does not include that repository's proof audit or downloads.
Machine-specific checkout paths in retained benchmark metadata were replaced
with `<checkout>`; numerical results were left unchanged.

## References and data

* `fluid/reference/tgv_re1600_spectral512.gdiag`: published 512³ spectral
  Taylor–Green reference, 1st International Workshop on High-Order CFD Methods,
  [original source](https://cfd.ku.edu/hiocfd/spectral_Re1600_512.gdiag).
* `fluid/data/` and crate-local CSV tables: generated reference, thermophysical,
  refrigerant and element data. Sources and generation scripts are retained
  alongside the tables; estimated quantities are identified in code.
* Solar reference downloads are not vendored; their published sources and
  fetch instructions are in `fluid/scripts/fetch_solar_reference.py` and
  `analysis/gravity_solar.md`.
* `analysis/photographic_scene/assets/`: generated material appearance maps.
  Exact prompts, provenance and their limits are retained in that directory's
  README and `texture-prompts.json`; they are not measured material data.
* Retained screenshots, plots and render movies document the numerical
  experiments described alongside them. They are not substitutes for their
  stated exact/published validation checks.

Mie coefficients follow Bohren & Huffman (1983) and Wiscombe (1980); optical
constants in the browser demo follow the source comments in
`fluid/crates/emergence-web/src/lib.rs`. Fluid and engineering material
references remain in the original model comments and `analysis/` write-ups.

Dependency licenses are declared by their upstream Cargo packages. Optional
STEP/IGES import uses separately installed gmsh/OpenCASCADE; these tools are
not bundled here.
