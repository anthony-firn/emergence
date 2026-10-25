//! Rendering and visual effects for the CAD kernel and the solvers.
//!
//! The point of this crate is that it draws the *same* geometry the
//! solvers use, straight from the kernel, with no export step and no
//! second representation to keep in sync. A `cad::Solid` is sphere-traced
//! directly; a simulation field is sampled where the ray lands; a
//! parametric model is re-evaluated per frame to animate it.
//!
//! ```text
//!   cad::Solid ----sphere trace (trace.rs)----+
//!   simulation field (field.rs) --------------+--> Render (image.rs)
//!   camera (camera.rs) -----------------------+     beauty + passes
//!        |                                          depth, normal, id,
//!   scene: materials, lights (scene.rs)              position, steps,
//!        |                                           field
//!   integrator (shade.rs): direct lighting,
//!        spectral path tracing via radiation::optics
//!        |
//!   anim.rs: parameter curves driving cad::model, camera paths
//! ```
//!
//! # Why sphere tracing is sound here
//!
//! The CAD kernel guarantees `|sdf(p)| <= dist(p, boundary)` at every node
//! of a model. A ray step of `|sdf(p)|` therefore cannot cross a surface,
//! so the march converges on the first crossing and **cannot tunnel
//! through thin geometry** — the standard failure of hand-written distance
//! estimators. [`trace`] states the argument and the one thing it does not
//! cover, which is cost near a tangency.
//!
//! # Rays are the coarse rung
//!
//! [`wave`] puts the ray tracer where it belongs: rung zero of an optical
//! ladder whose finer rungs are Fraunhofer, Fresnel and
//! Rayleigh-Sommerfeld diffraction, with the **Fresnel number** as the
//! indicator that decides between them. Ask for a tolerance and
//! `emergence::Ladder` evaluates the cheapest rung that can meet it —
//! which is what makes accuracy and efficiency the same question rather
//! than opposed ones. Thin-film colour, partial coherence and photon
//! counting live there too.
//!
//! # What is borrowed rather than rewritten
//!
//! * `radiation::optics` for Planck emission, Fresnel reflectance,
//!   refraction and the CIE observer — so the spectral and colour physics
//!   is the implementation already validated in that crate;
//! * `radiation::math` for rays, bounds and the deterministic generator;
//! * `waveoptics::color` for the sRGB transfer function, gamut handling
//!   and PNG output;
//! * `waveoptics::diffraction` for the wave propagators *and the rigorous
//!   error bounds* they carry, `waveoptics::source` for the coherence of a
//!   real spectrum, `waveoptics::detector` for photon statistics;
//! * `photonics::multilayer` for the complex amplitude reflectance of a
//!   layered surface;
//! * `emergence` for the ladder that chooses among them.
//!
//! What is new here is the geometry path (implicit, exact, untessellated),
//! the render passes, the field visualisation and the animation layer.
//!
//! # Determinism
//!
//! Every sample comes from a counter-based seed built from the pixel and
//! the sample index, so a render is reproducible to the bit regardless of
//! thread scheduling, and two renders of the same scene differ only where
//! the scene differs. The tests assert it.

pub mod anim;
pub mod camera;
pub mod caustic;
pub mod field;
pub mod image;
pub mod pbr;
pub mod scene;
pub mod shade;
pub mod trace;
pub mod wave;

pub use camera::Camera;
pub use image::{Channel, Image, Pass, Render};
pub use scene::{Light, Material, Object, Scene};
pub use shade::{render, RenderOptions};
pub use trace::{march, Hit, MarchOptions, Trace};
