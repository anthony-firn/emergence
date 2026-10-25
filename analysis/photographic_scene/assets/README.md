# Material maps

`oak-albedo.png` and `linen-albedo.png` were created with the built-in imagegen
tool as flat material textures. The exact prompts are in `texture-prompts.json`.
They are appearance estimates, not measured spectral reflectance or height.
Any residual shading in a generated map is an appearance limitation.

The Rust renderer decodes sRGB to linear light, builds a filtered image pyramid,
maps each texture onto geometry, and estimates micro-relief from local luminance
gradients. Relief scales and material roughness are explicitly estimates in the
example. Folds, liquid and bubble interfaces, reflected window structure,
occlusion, indirect illumination, camera blur and final pixels are rendered by
Rust; the generated full-scene calibration image is not sampled by the renderer.
