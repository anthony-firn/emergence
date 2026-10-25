use crate::math::{Bounds, Ray, Vec3};
use crate::optics::Spectrum;

#[derive(Clone, Debug)]
pub struct Medium {
    pub absorption: Spectrum,
    pub scattering: Spectrum,
    pub g: f64,
    pub n: Spectrum,
}
impl Medium {
    pub fn grey(absorption: f64, scattering: f64, g: f64) -> Self {
        assert!(g.abs() < 1.0);
        Self {
            absorption: Spectrum::constant(absorption),
            scattering: Spectrum::constant(scattering),
            g,
            n: Spectrum::constant(1.0),
        }
    }
    pub fn coefficients(&self, l: f64, scale: f64) -> (f64, f64) {
        (self.absorption.at(l) * scale, self.scattering.at(l) * scale)
    }
}
#[derive(Clone, Copy, Debug)]
pub struct Voxel {
    pub medium: usize,
    pub scale: f64,
    pub temperature: f64,
}
#[derive(Clone, Copy, Default)]
pub(crate) struct Coefficients {
    pub a: f64,
    pub s: f64,
    pub n2: f64,
    pub g: f64,
}
/// Inline storage avoids allocating per spectral path in ordinary scenes.
pub(crate) struct SpectralMedia {
    inline: [Coefficients; 8],
    extra: Vec<Coefficients>,
}
impl SpectralMedia {
    pub fn new(media: &[Medium], lambda: f64) -> Self {
        let mut result = Self {
            inline: [Coefficients::default(); 8],
            extra: Vec::new(),
        };
        for (i, m) in media.iter().enumerate() {
            let c = Coefficients {
                a: m.absorption.at(lambda),
                s: m.scattering.at(lambda),
                n2: m.n.at(lambda).powi(2),
                g: m.g,
            };
            if i < 8 {
                result.inline[i] = c;
            } else {
                result.extra.push(c);
            }
        }
        result
    }
    pub fn get(&self, i: usize) -> Coefficients {
        if i < 8 {
            self.inline[i]
        } else {
            self.extra[i - 8]
        }
    }
}
/// Cartesian finite-volume snapshot; k varies fastest, same as compressible.
#[derive(Clone, Debug)]
pub struct Volume {
    pub n: [usize; 3],
    pub bounds: Bounds,
    pub media: Vec<Medium>,
    pub cells: Vec<Voxel>,
}
impl Volume {
    pub fn new(n: [usize; 3], bounds: Bounds, media: Vec<Medium>, cells: Vec<Voxel>) -> Self {
        assert!(n.iter().all(|&n| n > 0));
        assert!(bounds
            .lo
            .0
            .iter()
            .chain(&bounds.hi.0)
            .all(|v| v.is_finite()));
        assert_eq!(n.iter().product::<usize>(), cells.len());
        assert!((0..3).all(|d| bounds.hi.0[d] > bounds.lo.0[d]));
        assert!(cells.iter().all(|v| v.medium < media.len()
            && v.scale >= 0.0
            && v.scale.is_finite()
            && v.temperature >= 0.0
            && v.temperature.is_finite()));
        assert!(media.iter().all(|m| m.g.abs() < 1.0));
        Self {
            n,
            bounds,
            media,
            cells,
        }
    }
    pub fn dx(&self) -> Vec3 {
        Vec3(std::array::from_fn(|d| {
            (self.bounds.hi.0[d] - self.bounds.lo.0[d]) / self.n[d] as f64
        }))
    }
    pub fn cell_volume(&self) -> f64 {
        self.dx().0.iter().product()
    }
    pub fn index(&self, i: [usize; 3]) -> usize {
        (i[0] * self.n[1] + i[1]) * self.n[2] + i[2]
    }
    pub fn coefficients(&self, q: usize, l: f64) -> (f64, f64) {
        let v = self.cells[q];
        self.media[v.medium].coefficients(l, v.scale)
    }
    pub fn segments(&self, r: Ray, limit: f64) -> Segments<'_> {
        Segments::new(self, r, limit)
    }
    pub fn transmittance(&self, r: Ray, limit: f64, l: f64) -> f64 {
        self.transmittance_cached(r, limit, &SpectralMedia::new(&self.media, l))
    }
    pub(crate) fn transmittance_cached(&self, r: Ray, limit: f64, media: &SpectralMedia) -> f64 {
        let tau: f64 = self
            .segments(r, limit)
            .map(|s| {
                let v = self.cells[s.cell];
                let c = media.get(v.medium);
                (c.a * v.scale + c.s * v.scale) * (s.end - s.start)
            })
            .sum();
        (-tau).exp()
    }
}
#[derive(Clone, Copy, Debug)]
pub struct Segment {
    pub cell: usize,
    pub start: f64,
    pub end: f64,
}
/// Allocation-free voxel traversal; simultaneous face/edge/corner crossings.
pub struct Segments<'a> {
    v: &'a Volume,
    r: Ray,
    i: [isize; 3],
    next: [f64; 3],
    t: f64,
    end: f64,
    done: bool,
    delta: [f64; 3],
}
impl<'a> Segments<'a> {
    fn new(v: &'a Volume, r: Ray, limit: f64) -> Self {
        let interval = v.bounds.interval(r, limit);
        let (t, end) = interval.unwrap_or((0.0, 0.0));
        let dx = v.dx();
        // Nudge only cell selection, never the path length.
        let p =
            r.at(t + (end - t).min(dx.0.into_iter().fold(f64::INFINITY, f64::min) * 1e-9) * 0.5);
        let i = std::array::from_fn(|d| {
            (((p.0[d] - v.bounds.lo.0[d]) / dx.0[d]).floor() as isize).clamp(0, v.n[d] as isize - 1)
        });
        let next = std::array::from_fn(|d| {
            if r.dir.0[d] == 0.0 {
                f64::INFINITY
            } else {
                let plane = i[d] + if r.dir.0[d] > 0.0 { 1 } else { 0 };
                ((v.bounds.lo.0[d] + plane as f64 * dx.0[d] - r.origin.0[d]) / r.dir.0[d]).max(t)
            }
        });
        Self {
            v,
            r,
            i,
            next,
            t,
            end,
            done: interval.is_none(),
            delta: std::array::from_fn(|d| dx.0[d] / r.dir.0[d].abs()),
        }
    }
}
impl Iterator for Segments<'_> {
    type Item = Segment;
    fn next(&mut self) -> Option<Segment> {
        if self.done {
            return None;
        }
        let finish = self.next.into_iter().fold(self.end, f64::min);
        let result = Segment {
            cell: self.v.index(self.i.map(|q| q as usize)),
            start: self.t,
            end: finish,
        };
        self.t = finish;
        if finish >= self.end {
            self.done = true;
        } else {
            for d in 0..3 {
                if self.next[d] <= finish {
                    self.i[d] += if self.r.dir.0[d] > 0.0 { 1 } else { -1 };
                    self.next[d] += self.delta[d];
                    if self.i[d] < 0 || self.i[d] >= self.v.n[d] as isize {
                        self.done = true;
                    }
                }
            }
        }
        Some(result)
    }
}
