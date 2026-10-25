# Browser adapters: models, validation and limitations

`fluid/crates/emergence-web` exposes a small C ABI to JavaScript; the model
implementations in `emergence`, `photonics`, `spectral-core` and `project` are
used directly. `web/worker.js` owns one WebAssembly instance. Model state never
leaves the browser. Rust computation uses f64; the flow display is an f32 slice.

## Exact and reference checks

* Planar Taylor–Green, 16³, Re = 40, ten steps to t = 0.4: the native and
  WebAssembly builds both gave mean energy `0.240197359788081`, relative error
  **8.09 × 10⁻¹⁶** against `¼ exp(−4νt)`. The native unit test requires relative
  energy error below 10⁻¹² and normalized divergence below 10⁻¹⁴. The display
  curl is also checked against `2 sin x sin y exp(−2νt)` to 10⁻⁶ (f32 buffer).
* The optical adapter calls the real resonance-aware closure builder and
  ladder, audits every answer against the series and plots 100 independently
  evaluated sizes. It inherits the series' validated absolute truncation
  claim of 2 × 10⁻⁸. The wrapper tests selection and an unmet 10⁻¹⁰ request;
  the domain's original tests check claims against the series separately.
  An additional WebAssembly sweep checked 1,212 plotted/selected results across
  all three materials and four tolerances; each discrepancy was below its
  claim plus the reference's 2 × 10⁻⁸ claim (worst ratio 0.556).
* The browser motor-driver graph keeps the original thermal/electrical
  definitions and three relevant requirements. Default nine fins and 2.5 m/s
  gives met junction and pressure checks and an undecided trace drop. At six
  fins and 3.5 m/s the laminar range fails and the junction check is unavailable.
  This preserves the original model's response; it is not an experimental
  hardware validation. The native CAD, printing and microstrip portions are
  omitted from the browser calculation and retained in the full example.

These checks are in `cargo test -p emergence-web`; all four tests together
initially took 0.25 s after compilation on the development machine. Numerical
budgets are test criteria chosen before measurement, not fitted constants.

## Scope

The optical view is one homogeneous sphere at a fixed monochromatic refractive
index. It is a building block for participating-media rendering, not a sky
forecast or a complete cloud model. Every point carries a rigor label;
asymptotic and validated claims are not promoted to rigorous bounds.

The fluid view uses CFL-limited IF-RK4 on an explicitly selected coarse grid.
The spectral-tail indicator is a warning, not an error bound and not automatic
emergence selection. Native rigorous grid-selection results are documented
separately. The low browser resolutions cannot validate turbulent dissipation
at Re = 1600. The planar preset is an exact solution, independent of the
grid's suitability for the separate 3D flow.

The design view preserves model assumptions and the rigor of the propagated
intervals. In particular the developing-channel discretization error is an
estimate inside a laminar model. Meeting a requirement in that model is not a
validated certification of the physical device.

The browser uses bounded Workers, keyboard-accessible tabs, explicit start/
pause controls, and stops fluid playback in a hidden tab. A 45-second model
timeout terminates the owned Worker and reports a reloadable error. Build
metadata records the source commit and SHA-256 of the deployed binary.
