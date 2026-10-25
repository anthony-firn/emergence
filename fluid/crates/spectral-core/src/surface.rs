//! Surface-constrained walls: the velocity is constrained only at marker
//! points on the true body surface (u(X_m) = u_s), with Lagrange-multiplier
//! point forces spread onto the retained modes by the non-uniform FFT.
//!
//! Volume penalization constrains the band-limited velocity at every solid
//! grid point — more constraints than retained modes — which over-determines
//! it and makes walls first order with staircase noise (analysis/
//! stokes_array_check.md). Here the wall position is exact (no staircase),
//! the number of constraints is set by the marker spacing relative to the
//! resolution 2π/(2k_max+1), and the region inside the body carries a
//! fictitious flow that is excluded from fluid averages analytically.

use num_complex::Complex64;

use crate::grid::Grid;
use crate::penalization::Shape;
use crate::solver::Spec;

#[derive(Clone, Debug)]
pub struct Surface {
    pub points: Vec<[f64; 3]>,
    pub u_solid: [f64; 3],
    /// Per-marker prescribed velocity (overrides `u_solid`), e.g. boundary
    /// data interpolated from a coarse grid for a refinement patch.
    pub u_markers: Option<Vec<[f64; 3]>>,
    /// the analytic body (for solid-region integrals)
    pub shape: Shape,
    /// marker spacing used to build `points`
    pub spacing: f64,
    /// smooth forcing extension (spectral accuracy); `None` = point forces
    pub extension: Option<Extension>,
    /// centres of the extension's tangential Gaussians. `None`: the markers
    /// themselves (square system, GMRES); `Some`: fewer centres than
    /// collocation markers (over-determined least squares, CGLS), which is
    /// robust when the shell is thin relative to the grid.
    pub basis: Option<Vec<[f64; 3]>>,
    /// Further parts of a compound surface (see `Surface::union`); `self`'s
    /// own shape/extension/spacing describe part 0.
    pub extra: Vec<Surface>,
    /// (marker range, centre range) of each part, including part 0
    pub parts: Vec<(std::ops::Range<usize>, std::ops::Range<usize>)>,
}

const TWO_PI: f64 = 2.0 * std::f64::consts::PI;

/// Resolution length of a grid's retained modes, 2π/(2k_max + 1).
pub fn resolution(g: &Grid) -> f64 {
    TWO_PI / (2 * g.kmax() + 1) as f64
}

impl Surface {
    /// Sphere with quasi-uniform (Fibonacci) markers at mean spacing `spacing`.
    pub fn sphere(center: [f64; 3], radius: f64, spacing: f64) -> Self {
        let n = ((4.0 * std::f64::consts::PI * radius * radius) / (spacing * spacing)).ceil().max(12.0) as usize;
        let golden = std::f64::consts::PI * (3.0 - 5f64.sqrt());
        let points = (0..n)
            .map(|i| {
                let z = 1.0 - (2 * i + 1) as f64 / n as f64;
                let r = (1.0 - z * z).sqrt();
                let phi = golden * i as f64;
                [center[0] + radius * r * phi.cos(), center[1] + radius * r * phi.sin(), center[2] + radius * z]
            })
            .collect();
        Self { points, u_solid: [0.0; 3], u_markers: None, shape: Shape::Sphere { center, radius }, spacing, extension: None, basis: None, extra: Vec::new(), parts: Vec::new() }
    }

    /// Two plane walls y = center ± half_width (solid outside), markers on a
    /// regular lattice in x, z at spacing ≈ `spacing`.
    pub fn channel_y(center: f64, half_width: f64, spacing: f64) -> Self {
        let m = (TWO_PI / spacing).ceil() as usize;
        let d = TWO_PI / m as f64;
        let mut points = Vec::with_capacity(2 * m * m);
        for y in [center - half_width, center + half_width] {
            for i in 0..m {
                for k in 0..m {
                    points.push([(i as f64 + 0.5) * d, y, (k as f64 + 0.5) * d]);
                }
            }
        }
        Self { points, u_solid: [0.0; 3], u_markers: None, shape: Shape::ChannelY { center, half_width }, spacing, extension: None, basis: None, extra: Vec::new(), parts: Vec::new() }
    }

    /// Solid z-cylinder (axis through `center` in x, y), markers on rings
    /// in z at spacing ≈ `spacing`, staggered between rings.
    pub fn cylinder_z(center: [f64; 2], radius: f64, spacing: f64) -> Self {
        let nz = (TWO_PI / spacing).ceil() as usize;
        let nt = ((TWO_PI * radius) / spacing).ceil().max(6.0) as usize;
        let (dz, dt) = (TWO_PI / nz as f64, TWO_PI / nt as f64);
        let mut points = Vec::with_capacity(nz * nt);
        for k in 0..nz {
            let off = if k % 2 == 0 { 0.0 } else { 0.5 };
            for i in 0..nt {
                let th = (i as f64 + off) * dt;
                points.push([center[0] + radius * th.cos(), center[1] + radius * th.sin(), (k as f64 + 0.5) * dz]);
            }
        }
        Self { points, u_solid: [0.0; 3], u_markers: None, shape: Shape::CylinderZ { center, radius }, spacing, extension: None, basis: None, extra: Vec::new(), parts: Vec::new() }
    }

    /// Fluid inside a z-cylinder (pipe / container wall), markers as `cylinder_z`.
    pub fn cylinder_z_inverse(center: [f64; 2], radius: f64, spacing: f64) -> Self {
        let mut s = Self::cylinder_z(center, radius, spacing);
        s.shape = Shape::CylinderZInverse { center, radius };
        s
    }

    /// Forces on each part from extension amplitudes `amp` (3 blocks of the
    /// total centre count): −dV Σ_m colsum[m]·A[m] over the part's centres.
    pub fn part_forces_extension(&self, op: &ExtensionOp, amp: &[f64]) -> Vec<[f64; 3]> {
        let nb = op.colsum.len();
        let ranges: Vec<std::ops::Range<usize>> = if self.parts.is_empty() { vec![0..nb] } else { self.parts.iter().map(|p| p.1.clone()).collect() };
        ranges
            .iter()
            .map(|r| [0, 1, 2].map(|c| -r.clone().map(|m| op.colsum[m] * amp[c * nb + m]).sum::<f64>()))
            .collect()
    }

    /// Forces on each part from point forces `f` (3 blocks of the marker count).
    pub fn part_forces_point(&self, f: &[f64]) -> Vec<[f64; 3]> {
        let nm = self.points.len();
        let ranges: Vec<std::ops::Range<usize>> = if self.parts.is_empty() { vec![0..nm] } else { self.parts.iter().map(|p| p.0.clone()).collect() };
        ranges.iter().map(|r| [0, 1, 2].map(|c| -r.clone().map(|m| f[c * nm + m]).sum::<f64>())).collect()
    }

    /// Fluid inside a sphere (refinement-patch region), markers as `sphere`.
    pub fn sphere_inverse(center: [f64; 3], radius: f64, spacing: f64) -> Self {
        let mut s = Self::sphere(center, radius, spacing);
        s.shape = Shape::SphereInverse { center, radius };
        s
    }

    /// Prescribed velocity component `c` at marker `m`.
    #[inline]
    pub fn target(&self, m: usize, c: usize) -> f64 {
        match &self.u_markers {
            Some(v) => v[m][c],
            None => self.u_solid[c],
        }
    }

    /// Exact box integral over the solid region of the field with
    /// coefficients `uh` (one component), using the analytic shape.
    pub fn solid_integral(&self, g: &Grid, uh: &Spec) -> f64 {
        let mut acc = 0.0;
        for idx in 0..g.spec_len() {
            if !g.kept(idx) {
                continue;
            }
            let k = g.k(idx);
            let w: Complex64 = match self.shape {
                Shape::SphereInverse { center: c, radius: r } => {
                    let kn = (k[0] * k[0] + k[1] * k[1] + k[2] * k[2]).sqrt();
                    let ball = if kn == 0.0 {
                        4.0 / 3.0 * std::f64::consts::PI * r.powi(3)
                    } else {
                        let x = kn * r;
                        4.0 * std::f64::consts::PI * r.powi(3) * (x.sin() - x * x.cos()) / x.powi(3)
                    };
                    let total = if kn == 0.0 { TWO_PI.powi(3) } else { 0.0 };
                    Complex64::from_polar(1.0, k[0] * c[0] + k[1] * c[1] + k[2] * c[2]) * (total - ball)
                }
                Shape::Sphere { center: c, radius: r } => {
                    let kn = (k[0] * k[0] + k[1] * k[1] + k[2] * k[2]).sqrt();
                    let vol = if kn == 0.0 {
                        4.0 / 3.0 * std::f64::consts::PI * r.powi(3)
                    } else {
                        let x = kn * r;
                        // 4πR³ j₁(x)/x, j₁(x) = (sin x − x cos x)/x²
                        4.0 * std::f64::consts::PI * r.powi(3) * (x.sin() - x * x.cos()) / x.powi(3)
                    };
                    Complex64::from_polar(vol, k[0] * c[0] + k[1] * c[1] + k[2] * c[2])
                }
                Shape::CylinderZ { center: c, radius: r } | Shape::CylinderZInverse { center: c, radius: r } => {
                    if k[2] != 0.0 {
                        Complex64::default()
                    } else {
                        // disk transform 2πR² J₁(κR)/(κR), times the height 2π
                        let kp = (k[0] * k[0] + k[1] * k[1]).sqrt();
                        let disk = if kp == 0.0 { std::f64::consts::PI * r * r } else { 2.0 * std::f64::consts::PI * r * r * bessel_j1(kp * r) / (kp * r) };
                        let v = match self.shape {
                            Shape::CylinderZ { .. } => disk,
                            _ => (if kp == 0.0 { TWO_PI * TWO_PI } else { 0.0 }) - disk,
                        };
                        Complex64::from_polar(TWO_PI * v, k[0] * c[0] + k[1] * c[1])
                    }
                }
                Shape::ChannelY { center, half_width: h } => {
                    if k[0] != 0.0 || k[2] != 0.0 {
                        Complex64::default()
                    } else {
                        // solid = box − fluid slab |y − center| < h
                        let slab = if k[1] == 0.0 { 2.0 * h } else { 2.0 * (k[1] * h).sin() / k[1] };
                        let total = if k[1] == 0.0 { TWO_PI } else { 0.0 };
                        Complex64::from_polar(TWO_PI * TWO_PI, k[1] * center) * (total - slab)
                    }
                }
            };
            acc += g.weight(idx) * (uh[idx] * w).re;
        }
        acc
    }

    /// Volume of the solid region.
    pub fn solid_volume(&self) -> f64 {
        match self.shape {
            Shape::Sphere { radius, .. } => 4.0 / 3.0 * std::f64::consts::PI * radius.powi(3),
            Shape::ChannelY { half_width, .. } => TWO_PI * TWO_PI * (TWO_PI - 2.0 * half_width),
            Shape::SphereInverse { radius, .. } => TWO_PI.powi(3) - 4.0 / 3.0 * std::f64::consts::PI * radius.powi(3),
            Shape::CylinderZ { radius, .. } => TWO_PI * std::f64::consts::PI * radius * radius,
            Shape::CylinderZInverse { radius, .. } => TWO_PI * (TWO_PI * TWO_PI - std::f64::consts::PI * radius * radius),
        }
    }
}

// ---------------------------------------------------------------------------
// Smooth forcing extension

/// Smooth forcing extension: instead of point forces on the surface, the
/// wall force is a C^∞ field in a shell of thickness `width` just inside
/// the body, ψ(d)·Σ_j a_j exp(−|π(x) − X_j|²/(2σ²)), with d the depth into
/// the solid, π(x) the closest surface point and ψ a C^∞ bump on (0, width)
/// that vanishes to all orders at the wall. Point forces put a kink in the
/// velocity at the wall, so its Fourier series converges only to first
/// order (plane channel: −24/−11/−5.8/−3.9% at N = 16/32/64/96); with a
/// smooth total forcing the solution is smooth and convergence is spectral.
#[derive(Clone, Copy, Debug)]
pub struct Extension {
    pub width: f64,
    /// tangential Gaussian width, in units of the marker spacing
    pub sigma_rel: f64,
}

/// Exponent p of the radial bump ψ_p(s) = exp(1 − 1/(1 − s²)^p); its
/// Fourier tail decays like exp(−c k^{p/(p+1)}), so a larger p wins once
/// the shell is well resolved but loses on coarse grids. Measured on plane
/// Poiseuille (shell 1.0): p = 1 → −4.8e-5 / −4.7e-5 at 64³ / 96³,
/// p = 2 → −9.3e-5 / +1.95e-5. Rule: p = 2 when the shell spans ≥ 12 cells,
/// else 1; `BUMP_P` in the environment overrides.
pub fn bump_exponent_for(cells_across_shell: f64) -> f64 {
    if let Some(p) = std::env::var("BUMP_P").ok().and_then(|v| v.parse().ok()) {
        return p;
    }
    if cells_across_shell >= 12.0 {
        2.0
    } else {
        1.0
    }
}

/// C^∞ bump on (0, w), 1 at the centre, all derivatives zero at 0 and w.
pub fn bump(d: f64, w: f64) -> f64 {
    bump_p(d, w, 1.0)
}

/// C^∞ bump of Gevrey type p on (0, w): exp(1 − 1/(1 − s²)^p), s = 2d/w − 1.
pub fn bump_p(d: f64, w: f64, p: f64) -> f64 {
    if d <= 0.0 || d >= w {
        return 0.0;
    }
    let s = 2.0 * d / w - 1.0;
    (1.0 - 1.0 / (1.0 - s * s).powf(p)).exp()
}

/// Sparse extension operator on the N grid: for each shell grid point,
/// the markers it couples to and their weights.
pub struct ExtensionOp {
    /// (grid point index, [(marker, weight)])
    pub rows: Vec<(usize, Vec<(usize, f64)>)>,
    /// Σ of weights per centre (∫ E·A dx = dV Σ_m colsum[m]·A[m])
    pub colsum: Vec<f64>,
}

/// Bessel J₁ by the trapezoid rule on J₁(x) = (1/π)∫₀^π cos(τ − x sin τ) dτ
/// (a periodic integrand, so the rule converges geometrically once the
/// point count exceeds ~|x|).
pub fn bessel_j1(x: f64) -> f64 {
    let m = (x.abs() as usize + 40) * 2;
    let h = std::f64::consts::PI / m as f64;
    let mut s = 0.0;
    for i in 0..=m {
        let t = i as f64 * h;
        let w = if i == 0 || i == m { 0.5 } else { 1.0 };
        s += w * (t - x * t.sin()).cos();
    }
    s * h / std::f64::consts::PI
}

fn periodic_delta(a: f64, b: f64) -> f64 {
    let mut d = (a - b).rem_euclid(TWO_PI);
    if d > std::f64::consts::PI {
        d -= TWO_PI;
    }
    d
}

impl Surface {
    /// Depth into the solid (positive inside) and closest surface point.
    fn depth_and_projection(&self, x: [f64; 3]) -> (f64, [f64; 3]) {
        match self.shape {
            Shape::Sphere { center: c, radius: r } => {
                let v = [periodic_delta(x[0], c[0]), periodic_delta(x[1], c[1]), periodic_delta(x[2], c[2])];
                let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt().max(1e-300);
                (r - len, [c[0] + r * v[0] / len, c[1] + r * v[1] / len, c[2] + r * v[2] / len])
            }
            Shape::ChannelY { center, half_width: h } => {
                let dy = periodic_delta(x[1], center);
                (dy.abs() - h, [x[0], center + h * dy.signum(), x[2]])
            }
            Shape::SphereInverse { center: c, radius: r } => {
                let v = [periodic_delta(x[0], c[0]), periodic_delta(x[1], c[1]), periodic_delta(x[2], c[2])];
                let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt().max(1e-300);
                (len - r, [c[0] + r * v[0] / len, c[1] + r * v[1] / len, c[2] + r * v[2] / len])
            }
            Shape::CylinderZ { center: c, radius: r } | Shape::CylinderZInverse { center: c, radius: r } => {
                let v = [periodic_delta(x[0], c[0]), periodic_delta(x[1], c[1])];
                let len = (v[0] * v[0] + v[1] * v[1]).sqrt().max(1e-300);
                let proj = [c[0] + r * v[0] / len, c[1] + r * v[1] / len, x[2]];
                let depth = if matches!(self.shape, Shape::CylinderZ { .. }) { r - len } else { len - r };
                (depth, proj)
            }
        }
    }

    /// Use extension centres at `basis_spacing` (same generator as the
    /// markers), making the boundary system over-determined.
    pub fn with_basis(mut self, basis_spacing: f64) -> Self {
        let b = match self.shape {
            Shape::Sphere { center, radius } | Shape::SphereInverse { center, radius } => Surface::sphere(center, radius, basis_spacing),
            Shape::ChannelY { center, half_width } => Surface::channel_y(center, half_width, basis_spacing),
            Shape::CylinderZ { center, radius } | Shape::CylinderZInverse { center, radius } => Surface::cylinder_z(center, radius, basis_spacing),
        };
        self.basis = Some(b.points);
        self
    }

    /// Centres of the extension basis.
    pub fn centres(&self) -> &[[f64; 3]] {
        self.basis.as_deref().unwrap_or(&self.points)
    }

    /// Spacing of the extension centres.
    fn centre_spacing(&self) -> f64 {
        let own = self.part_centres().len();
        match &self.basis {
            None => self.spacing,
            Some(_) => match self.shape {
                Shape::Sphere { radius, .. } | Shape::SphereInverse { radius, .. } => {
                    (4.0 * std::f64::consts::PI * radius * radius / own as f64).sqrt()
                }
                Shape::ChannelY { .. } => TWO_PI / ((own / 2) as f64).sqrt(),
                Shape::CylinderZ { radius, .. } | Shape::CylinderZInverse { radius, .. } => (TWO_PI * TWO_PI * radius / own as f64).sqrt(),
            },
        }
    }

    /// Build the extension operator on grid `g` (couplings to the centres);
    /// compound surfaces merge their parts' rows.
    pub fn extension_op(&self, g: &Grid, ext: Extension, _spacing: f64) -> ExtensionOp {
        let mut merged: std::collections::BTreeMap<usize, Vec<(usize, f64)>> = Default::default();
        let mut offset = 0;
        for (k, part) in std::iter::once(self).chain(self.extra.iter()).enumerate() {
            let pe = if k == 0 { ext } else { part.extension.expect("every part needs an extension") };
            for (p, coup) in part.single_extension_rows(g, pe) {
                merged.entry(p).or_default().extend(coup.into_iter().map(|(m, w)| (m + offset, w)));
            }
            offset += part.part_centres().len();
        }
        let dv = g.dx().powi(3);
        let mut colsum = vec![0.0; offset];
        for coup in merged.values() {
            for &(m, w) in coup {
                colsum[m] += w * dv;
            }
        }
        ExtensionOp { rows: merged.into_iter().collect(), colsum }
    }

    fn single_extension_rows(&self, g: &Grid, ext: Extension) -> Vec<(usize, Vec<(usize, f64)>)> {
        let spacing = self.centre_spacing();
        let sigma = ext.sigma_rel * spacing;
        let cutoff2 = (4.0 * sigma).powi(2);
        let n = g.n;
        let centres = self.part_centres();
        let gevrey = bump_exponent_for(ext.width / g.dx());
        (0..g.real_len())
            .filter_map(|p| {
                let (i, j, l) = (p / (n * n), (p / n) % n, p % n);
                let x = [g.x(i), g.x(j), g.x(l)];
                let (d, proj) = self.depth_and_projection(x);
                let psi = bump_p(d, ext.width, gevrey);
                if psi == 0.0 {
                    return None;
                }
                let couplings: Vec<(usize, f64)> = centres
                    .iter()
                    .enumerate()
                    .filter_map(|(m, xm)| {
                        let r2: f64 = (0..3).map(|q| periodic_delta(proj[q], xm[q]).powi(2)).sum();
                        (r2 < cutoff2).then(|| (m, psi * (-r2 / (2.0 * sigma * sigma)).exp()))
                    })
                    .collect();
                (!couplings.is_empty()).then_some((p, couplings))
            })
            .collect()
    }

    /// Union of surfaces (e.g. a body and a patch's outer boundary). Markers,
    /// targets and extension centres are concatenated; each part keeps its
    /// own shape and extension.
    pub fn union(parts: Vec<Surface>) -> Surface {
        let mut it = parts.into_iter();
        let first = it.next().expect("at least one part");
        let rest: Vec<Surface> = it.collect();
        let all: Vec<&Surface> = std::iter::once(&first).chain(rest.iter()).collect();
        let mut points = Vec::new();
        let mut targets = Vec::new();
        let mut centres = Vec::new();
        let mut ranges = Vec::new();
        for p in &all {
            let (m0, c0) = (points.len(), centres.len());
            for m in 0..p.points.len() {
                targets.push([p.target(m, 0), p.target(m, 1), p.target(m, 2)]);
            }
            points.extend_from_slice(&p.points);
            centres.extend_from_slice(p.centres());
            ranges.push((m0..points.len(), c0..centres.len()));
        }
        let mut s = first.clone();
        s.points = points;
        s.u_markers = Some(targets);
        s.basis = Some(centres);
        s.extra = rest.into_iter().map(|mut p| {
            // keep the part's own centres so its rows can be rebuilt
            p.basis = Some(p.centres().to_vec());
            p.extra.clear();
            p
        }).collect();
        s.parts = ranges;
        s
    }

    /// Centres belonging to this part alone (for compound surfaces, part 0).
    fn part_centres(&self) -> &[[f64; 3]] {
        match self.parts.first() {
            Some((_, c)) => &self.centres()[c.clone()],
            None => self.centres(),
        }
    }
}
