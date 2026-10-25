//! Images and render passes.
//!
//! A render produces several buffers at once, not one: the beauty pass
//! plus whatever auxiliary channels were asked for. Compositing, denoising
//! and any kind of measurement off a render all need those channels, and
//! computing them in the same traversal is free — the ray has already been
//! marched.
//!
//! Colour output goes through `waveoptics::color`, so the transfer
//! function and the out-of-gamut handling are the ones already validated
//! there rather than a second implementation.

use std::collections::BTreeMap;
use std::path::Path;

use serde::{Deserialize, Serialize};

/// A linear-light RGB image.
#[derive(Clone, Debug)]
pub struct Image {
    pub width: usize,
    pub height: usize,
    pub pixels: Vec<[f64; 3]>,
}

impl Image {
    pub fn new(width: usize, height: usize) -> Self {
        Self {
            width,
            height,
            pixels: vec![[0.0; 3]; width * height],
        }
    }

    pub fn aspect(&self) -> f64 {
        self.width as f64 / self.height as f64
    }

    #[inline]
    pub fn at(&self, x: usize, y: usize) -> [f64; 3] {
        self.pixels[y * self.width + x]
    }

    /// The brightest channel anywhere, which is the natural exposure for
    /// an image with no other reference.
    pub fn peak(&self) -> f64 {
        self.pixels
            .iter()
            .flat_map(|p| p.iter().copied())
            .fold(0.0, f64::max)
    }

    /// Mean of each channel.
    pub fn mean(&self) -> [f64; 3] {
        let n = self.pixels.len().max(1) as f64;
        let mut s = [0.0; 3];
        for p in &self.pixels {
            for k in 0..3 {
                s[k] += p[k];
            }
        }
        s.map(|v| v / n)
    }

    /// Write an 8-bit sRGB PNG, mapping `white` to full scale.
    pub fn write_png(&self, path: impl AsRef<Path>, white: f64) -> std::io::Result<()> {
        waveoptics::color::write_srgb_png(path, self.width, self.height, &self.pixels, white)
    }

    /// Write a linear-light PFM, which keeps the full range for later
    /// grading.
    pub fn write_pfm(&self, path: impl AsRef<Path>) -> std::io::Result<()> {
        use std::io::Write;
        let mut out = Vec::new();
        write!(out, "PF\n{} {}\n-1.0\n", self.width, self.height)?;
        // PFM rows run bottom to top.
        for y in (0..self.height).rev() {
            for x in 0..self.width {
                for v in self.at(x, y) {
                    out.extend_from_slice(&(v as f32).to_le_bytes());
                }
            }
        }
        std::fs::write(path, out)
    }
}

/// A single-channel buffer: depth, a field value, a step count.
#[derive(Clone, Debug)]
pub struct Channel {
    pub width: usize,
    pub height: usize,
    pub values: Vec<f64>,
}

impl Channel {
    pub fn new(width: usize, height: usize, fill: f64) -> Self {
        Self {
            width,
            height,
            values: vec![fill; width * height],
        }
    }

    /// Smallest and largest finite value.
    pub fn range(&self) -> (f64, f64) {
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for &v in &self.values {
            if v.is_finite() {
                lo = lo.min(v);
                hi = hi.max(v);
            }
        }
        (lo, hi)
    }

    /// Write as a 16-bit greyscale PNG, scaled over `(lo, hi)`.
    pub fn write_png(&self, path: impl AsRef<Path>, lo: f64, hi: f64) -> std::io::Result<()> {
        let span = if (hi - lo).abs() > 0.0 { hi - lo } else { 1.0 };
        let v: Vec<u32> = self
            .values
            .iter()
            .map(|&x| {
                if x.is_finite() {
                    (((x - lo) / span).clamp(0.0, 1.0) * 65535.0).round() as u32
                } else {
                    0
                }
            })
            .collect();
        waveoptics::color::write_gray16_png(path, self.width, self.height, &v, 65535)
    }
}

/// Which auxiliary channels to produce alongside the beauty pass.
#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq, Eq, PartialOrd, Ord)]
#[serde(rename_all = "snake_case")]
pub enum Pass {
    /// Distance along the ray to the first hit.
    Depth,
    /// Surface normal, as an RGB image in `[-1, 1]`.
    Normal,
    /// Index of the object hit, plus one; zero for the background.
    ObjectId,
    /// World position of the hit.
    Position,
    /// March steps taken, which is where the render spent its time.
    Steps,
    /// A scalar field sampled at the hit point.
    Field,
}

/// Everything a render produced.
#[derive(Clone, Debug)]
pub struct Render {
    pub beauty: Image,
    pub channels: BTreeMap<Pass, Channel>,
    /// Per-pass RGB buffers, for the vector-valued passes.
    pub vectors: BTreeMap<Pass, Image>,
    /// Rays whose march ran out of steps: if this is not zero, the image
    /// has grazing geometry that was given up on, and the budget should go
    /// up.
    pub exhausted: usize,
    pub rays: usize,
}

impl Render {
    pub fn new(width: usize, height: usize) -> Self {
        Self {
            beauty: Image::new(width, height),
            channels: BTreeMap::new(),
            vectors: BTreeMap::new(),
            exhausted: 0,
            rays: 0,
        }
    }

    /// Fraction of rays that ran out of march steps.
    pub fn exhausted_fraction(&self) -> f64 {
        self.exhausted as f64 / self.rays.max(1) as f64
    }
}
