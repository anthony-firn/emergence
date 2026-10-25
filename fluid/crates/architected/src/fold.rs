//! Fold validation: a finite `N`-cell lattice against the homogenised beam.
//!
//! The effective continuum is the \(N \to \infty\) limit. For a finite beam
//! made of cells, two corrections survive:
//!
//! * a **boundary layer** at the clamped end, a few cells long, where the
//!   free-surface kinematics of the lattice do not match the continuum
//!   plane sections; and
//! * **shear** deflection, which for the bending-dominated honeycomb is not
//!   negligible because \(G^* \sim E^*\).
//!
//! Both shrink as the beam length in cells grows, so the tip deflection
//! relative to the Euler-Bernoulli continuum using \(E^* I\) is the natural
//! size-effect indicator. The functions here build a finite honeycomb
//! cantilever with the ordinary frame solver and measure it.

use crate::frame::{Frame, Member};
use std::collections::HashMap;

/// A finite patch of a 2D lattice.
#[derive(Clone, Debug)]
pub struct Patch {
    pub nodes: Vec<[f64; 2]>,
    pub members: Vec<Member>,
}

fn key(x: f64, y: f64) -> (i64, i64) {
    ((x * 1e7).round() as i64, (y * 1e7).round() as i64)
}

/// Rectangular patch of `nx` by `ny` hexagonal cells (vertical-walled
/// orientation). `theta` is the wall inclination, `h_over_l` the vertical
/// wall length ratio; the regular honeycomb is `theta = 30 deg`, `h = l`.
pub fn honeycomb_patch(nx: usize, ny: usize, t: f64, l: f64, h_over_l: f64, theta: f64) -> Patch {
    let h = h_over_l * l;
    let lc = l * theta.cos();
    let ls = l * theta.sin();
    let mut nodes: Vec<[f64; 2]> = Vec::new();
    let mut index: HashMap<(i64, i64), usize> = HashMap::new();
    let mut members = Vec::new();
    let mut add = |p: [f64; 2], nodes: &mut Vec<[f64; 2]>| -> usize {
        let k = key(p[0], p[1]);
        *index.entry(k).or_insert_with(|| {
            nodes.push(p);
            nodes.len() - 1
        })
    };
    for j in 0..=ny {
        for i in 0..=nx {
            let cx = 2.0 * lc * i as f64 + if j % 2 == 1 { lc } else { 0.0 };
            let cy = (h + ls) * j as f64;
            let tr = [cx + lc, cy + h / 2.0];
            let top = [cx, cy + h / 2.0 + ls];
            let tl = [cx - lc, cy + h / 2.0];
            let bl = [cx - lc, cy - h / 2.0];
            let bot = [cx, cy - h / 2.0 - ls];
            let br = [cx + lc, cy - h / 2.0];
            let verts = [tr, top, tl, bl, bot, br];
            for k in 0..6 {
                let a = add(verts[k], &mut nodes);
                let b = add(verts[(k + 1) % 6], &mut nodes);
                if a != b {
                    members.push(Member {
                        a,
                        b,
                        e: 1.0,
                        area: t,
                        inertia: t * t * t / 12.0,
                    });
                }
            }
        }
    }
    // De-duplicate shared members: the hexagon tiling reaches each interior
    // wall twice.
    let mut seen = std::collections::HashSet::new();
    members.retain(|m| {
        let k = if m.a < m.b { (m.a, m.b) } else { (m.b, m.a) };
        seen.insert(k)
    });
    Patch { nodes, members }
}

impl Patch {
    pub fn bounding_box(&self) -> ([f64; 2], [f64; 2]) {
        let mut lo = [f64::INFINITY; 2];
        let mut hi = [f64::NEG_INFINITY; 2];
        for p in &self.nodes {
            for d in 0..2 {
                lo[d] = lo[d].min(p[d]);
                hi[d] = hi[d].max(p[d]);
            }
        }
        (lo, hi)
    }

    /// Clamp every node within `tol` of the minimum `x`, and push down with
    /// total load `p` at the nodes within `tol` of the maximum `x`.
    pub fn cantilever(&self, p: f64, tol: f64) -> (Frame, Vec<[f64; 3]>) {
        let (lo, hi) = self.bounding_box();
        let mut frame = Frame::new(self.nodes.clone(), self.members.clone());
        for (i, node) in self.nodes.iter().enumerate() {
            if node[0] - lo[0] < tol {
                frame.supports[i] = [true, true, true];
            }
        }
        let tip_nodes: Vec<usize> = self
            .nodes
            .iter()
            .enumerate()
            .filter(|(_, n)| hi[0] - n[0] < tol)
            .map(|(i, _)| i)
            .collect();
        let mut loads = vec![[0.0f64; 3]; self.nodes.len()];
        // Distribute the load over the tip nodes weighted by proximity to the
        // mid-height, which approximates a point load without introducing a
        // spurious moment.
        let mid = 0.5 * (lo[1] + hi[1]);
        let mut wsum = 0.0;
        let mut w = vec![0.0; tip_nodes.len()];
        for (k, &i) in tip_nodes.iter().enumerate() {
            // tip nodes come in a vertical pair; weight the one nearer mid.
            w[k] = 1.0 / (1.0 + (self.nodes[i][1] - mid).abs());
            wsum += w[k];
        }
        for (k, &i) in tip_nodes.iter().enumerate() {
            loads[i][1] = -p * w[k] / wsum;
        }
        (frame, loads)
    }

    /// Tip deflection of a honeycomb cantilever. Returns
    /// `(delta_fem, length, height)`.
    pub fn cantilever_tip(&self, p: f64) -> (f64, f64, f64) {
        let (lo, hi) = self.bounding_box();
        let tol = 1e-6 * (hi[0] - lo[0]).max(1.0);
        let (frame, loads) = self.cantilever(p, tol);
        let (u, _) = frame.solve(&loads);
        let mid = 0.5 * (lo[1] + hi[1]);
        let mut best_w = 0.0;
        let mut delta = 0.0;
        for (i, node) in self.nodes.iter().enumerate() {
            if hi[0] - node[0] < tol {
                let w = 1.0 / (1.0 + (node[1] - mid).abs());
                if w > best_w {
                    best_w = w;
                    delta = u[i][1];
                }
            }
        }
        (delta, hi[0] - lo[0], hi[1] - lo[1])
    }
}
