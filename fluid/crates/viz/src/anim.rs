//! Animation: parameters over time.
//!
//! The models in this system are parametric, so animating one is sweeping
//! its parameters — there is no separate rig to build. A [`Timeline`] maps
//! frame numbers to parameter values, a [`cad::model::Model`] is rebuilt
//! per frame, and the camera follows a NURBS path. That makes a sequence
//! reproducible from a few numbers rather than from a scene file full of
//! baked keys.
//!
//! Interpolation is by [`Curve`]: constant, linear, or a Catmull-Rom
//! spline through the keys. Catmull-Rom passes through every key exactly
//! (the tests check that) and is `C^1`, so a parameter does not jerk at a
//! key; it can overshoot between keys, which is sometimes what you want
//! and sometimes not, so linear is there too.

use std::collections::BTreeMap;

use cad::nurbs::Nurbs;
use cad::Params;
use geometry::P3;
use serde::{Deserialize, Serialize};

/// How a parameter moves between its keys.
#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum Curve {
    /// Hold each key until the next.
    Constant,
    Linear,
    /// Catmull-Rom through the keys: passes through every one, `C^1`, and
    /// may overshoot between them.
    Smooth,
}

/// One animated parameter: `(time, value)` keys and how to join them.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Track {
    pub keys: Vec<(f64, f64)>,
    #[serde(default = "default_curve")]
    pub curve: Curve,
}

fn default_curve() -> Curve {
    Curve::Linear
}

impl Track {
    pub fn new(keys: Vec<(f64, f64)>, curve: Curve) -> Self {
        let mut keys = keys;
        keys.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        Self { keys, curve }
    }

    pub fn constant(value: f64) -> Self {
        Self::new(vec![(0.0, value)], Curve::Constant)
    }

    /// Value at time `t`. Before the first key and after the last, the
    /// nearest key's value is held.
    pub fn at(&self, t: f64) -> f64 {
        if self.keys.is_empty() {
            return 0.0;
        }
        if t <= self.keys[0].0 {
            return self.keys[0].1;
        }
        let last = self.keys.len() - 1;
        if t >= self.keys[last].0 {
            return self.keys[last].1;
        }
        let i = self
            .keys
            .windows(2)
            .position(|w| t >= w[0].0 && t <= w[1].0)
            .unwrap_or(0);
        let (t0, v0) = self.keys[i];
        let (t1, v1) = self.keys[i + 1];
        let span = t1 - t0;
        let f = if span > 0.0 { (t - t0) / span } else { 0.0 };
        match self.curve {
            Curve::Constant => v0,
            Curve::Linear => v0 + (v1 - v0) * f,
            Curve::Smooth => {
                // Catmull-Rom with the end tangents clamped.
                let vm = if i == 0 { v0 } else { self.keys[i - 1].1 };
                let vp = if i + 2 <= last { self.keys[i + 2].1 } else { v1 };
                let (f2, f3) = (f * f, f * f * f);
                0.5 * ((2.0 * v0)
                    + (-vm + v1) * f
                    + (2.0 * vm - 5.0 * v0 + 4.0 * v1 - vp) * f2
                    + (-vm + 3.0 * v0 - 3.0 * v1 + vp) * f3)
            }
        }
    }

    pub fn duration(&self) -> f64 {
        self.keys.last().map_or(0.0, |k| k.0)
    }
}

/// A set of named tracks, plus the frame rate they are played at.
#[derive(Clone, Debug, Default, Deserialize, Serialize)]
pub struct Timeline {
    pub tracks: BTreeMap<String, Track>,
    #[serde(default = "default_fps")]
    pub fps: f64,
    /// Number of frames; `None` means from the longest track.
    #[serde(default)]
    pub frames: Option<usize>,
}

fn default_fps() -> f64 {
    25.0
}

impl Timeline {
    pub fn new(fps: f64) -> Self {
        Self {
            tracks: BTreeMap::new(),
            fps,
            frames: None,
        }
    }

    pub fn with(mut self, name: &str, track: Track) -> Self {
        self.tracks.insert(name.to_string(), track);
        self
    }

    pub fn duration(&self) -> f64 {
        self.tracks.values().map(|t| t.duration()).fold(0.0, f64::max)
    }

    pub fn frame_count(&self) -> usize {
        self.frames
            .unwrap_or_else(|| (self.duration() * self.fps).round() as usize + 1)
    }

    pub fn time_of(&self, frame: usize) -> f64 {
        frame as f64 / self.fps
    }

    /// Parameter values at a frame, as overrides for a model build.
    pub fn params_at(&self, frame: usize) -> Params {
        let t = self.time_of(frame);
        let mut p = Params::new();
        for (name, track) in &self.tracks {
            p.set(name, track.at(t));
        }
        p
    }
}

/// A camera move along a NURBS path, aiming at another.
#[derive(Clone, Debug)]
pub struct CameraPath {
    pub eye: Nurbs,
    pub target: Nurbs,
    /// Sample the paths by arc length rather than parameter, so the move
    /// has constant speed.
    pub by_length: bool,
}

impl CameraPath {
    pub fn new(eye: Nurbs, target: Nurbs) -> Self {
        Self {
            eye,
            target,
            by_length: true,
        }
    }

    /// Eye and target at a normalised time in `[0, 1]`.
    pub fn at(&self, t: f64) -> (P3, P3) {
        let pick = |c: &Nurbs| -> P3 {
            let (lo, hi) = c.domain();
            if !self.by_length {
                return c.eval(lo + (hi - lo) * t.clamp(0.0, 1.0));
            }
            // Invert arc length on a fine table.
            let n = 512;
            let mut s = vec![0.0; n + 1];
            let mut prev = c.eval(lo);
            for i in 1..=n {
                let p = c.eval(lo + (hi - lo) * i as f64 / n as f64);
                s[i] = s[i - 1] + geometry::norm(geometry::sub(p, prev));
                prev = p;
            }
            let target = s[n] * t.clamp(0.0, 1.0);
            let mut j = 0;
            while j < n && s[j + 1] < target {
                j += 1;
            }
            let span = s[j + 1] - s[j];
            let f = if span > 0.0 { (target - s[j]) / span } else { 0.0 };
            c.eval(lo + (hi - lo) * (j as f64 + f) / n as f64)
        };
        (pick(&self.eye), pick(&self.target))
    }

    /// The camera at a normalised time, keeping the given field of view.
    pub fn camera(&self, t: f64, fov_y: f64) -> crate::Camera {
        let (eye, target) = self.at(t);
        crate::Camera::looking_at(eye, target, fov_y)
    }

    /// An orbit: the eye on a circle of `radius` about `centre` in the
    /// plane normal to `up`, looking at the centre.
    pub fn orbit(centre: P3, radius: f64, up: P3, start: P3) -> Self {
        let u = geometry::normalize(geometry::sub(start, centre));
        let v = geometry::cross(geometry::normalize(up), u);
        let eye = Nurbs::circle(centre, u, v, radius);
        let target = Nurbs::line(centre, centre);
        Self {
            eye,
            target,
            by_length: true,
        }
    }
}
