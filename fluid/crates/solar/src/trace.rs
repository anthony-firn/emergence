//! Adaptive magnetic field geometry, dx/ds = +/- B/|B|. These are spatial
//! curves, not photon paths or time-evolving plasma trajectories.
use crate::corona::Pfss;
use gravity::{norm, Vec3};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum End {
    Photosphere,
    SourceSurface,
    Null,
    LengthLimit,
    StepLimit,
    StepUnderflow,
}
impl End {
    pub fn label(self) -> &'static str {
        match self {
            Self::Photosphere => "photosphere",
            Self::SourceSurface => "source_surface",
            Self::Null => "near_null",
            Self::LengthLimit => "length_limit",
            Self::StepLimit => "step_limit",
            Self::StepUnderflow => "step_underflow",
        }
    }
}
#[derive(Clone, Copy, Debug)]
pub struct Options {
    /// Local position-error estimate and boundary location tolerance, in
    /// photospheric radii. This is not a rigorous global path error bound.
    pub tolerance: f64,
    pub max_step: f64,
    pub min_step: f64,
    pub max_length: f64,
    pub max_attempts: usize,
    /// Stop near a magnetic null; threshold is in tesla.
    pub null_field: f64,
}
impl Default for Options {
    fn default() -> Self {
        Self {
            tolerance: 1e-6,
            max_step: 0.025,
            min_step: 1e-10,
            max_length: 20.0,
            max_attempts: 20000,
            null_field: 1e-12,
        }
    }
}
impl Options {
    fn validate(self) -> Result<(), &'static str> {
        if !self.tolerance.is_finite()
            || self.tolerance <= 0.0
            || self.tolerance >= 0.01
            || !self.max_step.is_finite()
            || self.max_step <= 0.0
            || !self.min_step.is_finite()
            || self.min_step <= 0.0
            || self.min_step >= self.max_step
            || !self.max_length.is_finite()
            || self.max_length <= 0.0
            || self.max_attempts == 0
            || !self.null_field.is_finite()
            || self.null_field < 0.0
        {
            return Err("invalid magnetic trace options");
        }
        Ok(())
    }
}
pub struct Line {
    /// Cartesian positions in meters, ordered from one end to the other.
    pub points: Vec<Vec3>,
    pub ends: [End; 2],
    pub attempts: usize,
}
impl Line {
    /// Only successfully located boundaries determine magnetic connectivity.
    pub fn closed(&self) -> bool {
        self.ends == [End::Photosphere; 2]
    }
    pub fn open(&self) -> bool {
        self.ends.contains(&End::Photosphere) && self.ends.contains(&End::SourceSurface)
    }
}
enum Trial {
    Point(Vec3),
    Outside,
    Null,
}
fn direction(field: &Pfss, p: Vec3, sign: f64, o: Options) -> Result<Trial, &'static str> {
    let r = norm(p);
    let outer = field.source_radius() / field.photospheric_radius();
    if r < 1.0 - 8.0 * f64::EPSILON || r > outer * (1.0 + 8.0 * f64::EPSILON) {
        return Ok(Trial::Outside);
    }
    let b = field.magnetic_field(p * field.photospheric_radius())?;
    let magnitude = norm(b);
    if magnitude <= o.null_field {
        return Ok(Trial::Null);
    }
    Ok(Trial::Point(b * (sign / magnitude)))
}
fn rk4(field: &Pfss, p: Vec3, h: f64, sign: f64, o: Options) -> Result<Trial, &'static str> {
    let k1 = match direction(field, p, sign, o)? {
        Trial::Point(v) => v,
        other => return Ok(other),
    };
    let k2 = match direction(field, p + k1 * (h * 0.5), sign, o)? {
        Trial::Point(v) => v,
        other => return Ok(other),
    };
    let k3 = match direction(field, p + k2 * (h * 0.5), sign, o)? {
        Trial::Point(v) => v,
        other => return Ok(other),
    };
    let k4 = match direction(field, p + k3 * h, sign, o)? {
        Trial::Point(v) => v,
        other => return Ok(other),
    };
    let result = p + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (h / 6.0);
    let r = norm(result);
    if r < 1.0 || r > field.source_radius() / field.photospheric_radius() {
        return Ok(Trial::Outside);
    }
    Ok(Trial::Point(result))
}
struct Half {
    points: Vec<Vec3>,
    end: End,
    attempts: usize,
}
fn half(field: &Pfss, seed: Vec3, sign: f64, o: Options) -> Result<Half, &'static str> {
    let mut points = vec![seed];
    let mut p = seed;
    let mut length = 0.0;
    let mut h = o.max_step;
    let outer = field.source_radius() / field.photospheric_radius();
    for attempt in 0..o.max_attempts {
        let finish = |points, end| Half {
            points,
            end,
            attempts: attempt + 1,
        };
        let v = match direction(field, p, sign, o)? {
            Trial::Point(v) => v,
            Trial::Null => return Ok(finish(points, End::Null)),
            Trial::Outside => return Err("magnetic seed outside coronal domain"),
        };
        let r = norm(p);
        let radial = v.dot(p) / r;
        let boundary = if r - 1.0 <= o.tolerance && radial < 0.0 {
            Some((1.0, End::Photosphere))
        } else if outer - r <= o.tolerance && radial > 0.0 {
            Some((outer, End::SourceSurface))
        } else {
            None
        };
        if let Some((radius, end)) = boundary {
            // Projection is bounded by the requested boundary tolerance.
            points.push(p * (radius / r));
            return Ok(finish(points, end));
        }
        if length >= o.max_length {
            return Ok(finish(points, End::LengthLimit));
        }
        h = h.min(o.max_length - length);
        if h < o.min_step {
            return Ok(finish(points, End::StepUnderflow));
        }
        let full = rk4(field, p, h, sign, o)?;
        let first = rk4(field, p, h * 0.5, sign, o)?;
        let second = match first {
            Trial::Point(mid) => rk4(field, mid, h * 0.5, sign, o)?,
            other => other,
        };
        match (full, second) {
            (Trial::Point(a), Trial::Point(b)) => {
                let error = norm(a - b) / 15.0;
                let factor = if error > 0.0 {
                    (0.9 * (o.tolerance / error).powf(0.2)).clamp(0.2, 2.0)
                } else {
                    2.0
                };
                if error <= o.tolerance {
                    p = b;
                    points.push(p);
                    length += h;
                }
                h = (h * factor).min(o.max_step);
            }
            (Trial::Null, _) | (_, Trial::Null) => return Ok(finish(points, End::Null)),
            _ => h *= 0.5,
        }
    }
    Ok(Half {
        points,
        end: End::StepLimit,
        attempts: o.max_attempts,
    })
}
/// Trace both signs with a frozen PFSS field. Limit/null terminations remain
/// explicit and are never labelled open or closed. Input and output use SI.
pub fn trace(field: &Pfss, seed: Vec3, options: Options) -> Result<Line, &'static str> {
    options.validate()?;
    field.magnetic_field(seed)?;
    let r = field.photospheric_radius();
    let negative = half(field, seed / r, -1.0, options)?;
    let positive = half(field, seed / r, 1.0, options)?;
    let attempts = negative.attempts + positive.attempts;
    let ends = [negative.end, positive.end];
    let mut points = negative.points;
    points.reverse();
    points.extend(positive.points.into_iter().skip(1));
    for p in &mut points {
        *p = *p * r;
    }
    Ok(Line {
        points,
        ends,
        attempts,
    })
}
