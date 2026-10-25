//! Hysteresis fold: Stoner-Wohlfarth single-domain particles.
//!
//! A particle small enough to reverse by coherent rotation is one
//! macrospin with energy density
//! `e = -K (m . n)^2 - mu0 M_s H . m`.  The exact zero-temperature
//! switching field, in units of the anisotropy field `H_K = 2K/(mu0 M_s)`,
//! is the astroid
//!
//! ```text
//! h_sw(psi) = 1 / (cos^{2/3} psi + sin^{2/3} psi)^{3/2},
//! ```
//!
//! with `psi` the angle of the reversed applied field to the easy axis
//! (Stoner & Wohlfarth, *Phil. Trans. R. Soc. A* **240**, 599 (1948)).  For
//! an isotropic assembly of easy axes the remanence coercivity is
//! `H_c = 0.479 H_K` (exact; Stoner-Wohlfarth 1948, see also Cullity &
//! Graham, *Introduction to Magnetic Materials*).
//!
//! The fold is valid while one macrospin represents the particle.  The
//! single-domain indicator is the particle size compared with the exchange
//! length `l_ex` and the domain-wall width `delta`: once `d > D_c` (Kittel
//! 1946) a domain wall fits and the macrospin picture must be unfolded to
//! micromagnetics.

use crate::geometry::{dot, V3};
use crate::micromag::normalize;
use crate::rng::Rng;

/// Exact Stoner-Wohlfarth astroid field in units of `H_K`, for reversed
/// field at angle `psi` (radians) to the easy axis.
pub fn astroid_hsw(psi: f64) -> f64 {
    let c = psi.cos().abs().powf(2.0 / 3.0);
    let s = psi.sin().abs().powf(2.0 / 3.0);
    1.0 / (c + s).powf(1.5)
}

/// A normalised macrospin: fields in units of `H_K = 2K/(mu0 M_s)`, energy
/// `e = -(1/2)(m . n)^2 - h . m`.
pub struct Macrospin {
    pub alpha: f64,
    pub easy: V3,
    pub m: V3,
    /// Normalised saturation magnetisation (default 1).
    pub ms: f64,
    /// Local-shape demagnetising terms `(axis, N)`.
    pub demag: Vec<(V3, f64)>,
}

impl Macrospin {
    pub fn new(alpha: f64) -> Self {
        Self {
            alpha,
            easy: [1.0, 0.0, 0.0],
            m: [1.0, 0.0, 0.0],
            ms: 1.0,
            demag: Vec::new(),
        }
    }

    pub fn h_eff(&self, h: V3, m: V3) -> V3 {
        let a = dot(m, self.easy);
        let mut he = [
            a * self.easy[0] + h[0],
            a * self.easy[1] + h[1],
            a * self.easy[2] + h[2],
        ];
        for (axis, nfac) in &self.demag {
            let an = dot(m, *axis);
            for c in 0..3 {
                he[c] -= nfac * self.ms * an * axis[c];
            }
        }
        he
    }

    /// One RK4 step of the normalised LLG.
    pub fn step(&mut self, h: V3, dt: f64) {
        let rhs = |m: V3| -> V3 {
            let g = 1.0 / (1.0 + self.alpha * self.alpha);
            let he = self.h_eff(h, m);
            let a = cross(m, he);
            let b = cross(m, a);
            [
                -g * (a[0] + self.alpha * b[0]),
                -g * (a[1] + self.alpha * b[1]),
                -g * (a[2] + self.alpha * b[2]),
            ]
        };
        let k1 = rhs(self.m);
        let s2 = add(self.m, scale(k1, 0.5 * dt));
        let k2 = rhs(s2);
        let s3 = add(self.m, scale(k2, 0.5 * dt));
        let k3 = rhs(s3);
        let s4 = add(self.m, scale(k3, dt));
        let k4 = rhs(s4);
        let mut s = self.m;
        for c in 0..3 {
            s[c] += dt / 6.0 * (k1[c] + 2.0 * k2[c] + 2.0 * k3[c] + k4[c]);
        }
        self.m = normalize(s);
    }

    pub fn relax(&mut self, h: V3, dt: f64, steps: usize) {
        for _ in 0..steps {
            self.step(h, dt);
        }
    }

    /// LLG-measured switching field (units of `H_K`) for a reversed field at
    /// angle `psi` to the easy axis, found by bisection.
    pub fn switching_field_llg(psi: f64, alpha: f64, dt: f64, steps: usize) -> f64 {
        let dir = [-psi.cos(), psi.sin(), 0.0];
        let switched = |h: f64| -> bool {
            let mut ms = Macrospin::new(alpha);
            // A tiny transverse perturbation so the collinear saddle does not
            // trap the integrator.
            ms.m = normalize([ms.easy[0], ms.easy[1] + 1.0e-4, ms.easy[2]]);
            ms.relax([h * dir[0], h * dir[1], h * dir[2]], dt, steps);
            dot(ms.m, ms.easy) < 0.0
        };
        let mut lo = 0.0;
        let mut hi = 1.5;
        // Ensure the bracket contains the transition.
        while !switched(hi) && hi < 8.0 {
            hi *= 1.5;
        }
        for _ in 0..40 {
            let mid = 0.5 * (lo + hi);
            if switched(mid) {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        0.5 * (lo + hi)
    }
}

/// Remanence coercivity of an isotropic Stoner-Wohlfarth assembly, in units
/// of `H_K`, measured by averaging macrospin hysteresis loops over random
/// easy-axis orientations with the field along z.
pub fn random_assembly_coercivity(
    n_orient: usize,
    n_fields: usize,
    alpha: f64,
    dt: f64,
    relax_steps: usize,
    seed: u64,
) -> f64 {
    let mut rng = Rng::new(seed);
    let mut orientations = Vec::with_capacity(n_orient);
    for _ in 0..n_orient {
        orientations.push(rng.unit_vector());
    }
    let h_max = 3.0_f64;
    let fields: Vec<f64> = (0..n_fields)
        .map(|f| h_max * (1.0 - 2.0 * f as f64 / (n_fields - 1) as f64))
        .collect();
    let mut mbar = vec![0.0; n_fields];
    // Carry each macrospin through the field sweep so the hysteresis branch
    // is followed: at T=0 the state stays in the metastable minimum until
    // it disappears at the Stoner-Wohlfarth switching field.
    for &easy in &orientations {
        let mut ms = Macrospin::new(alpha);
        ms.easy = normalize(easy);
        ms.m = [0.0, 0.0, 1.0];
        for (fi, &h) in fields.iter().enumerate() {
            ms.relax([0.0, 0.0, h], dt, relax_steps);
            mbar[fi] += ms.m[2];
        }
    }
    for v in mbar.iter_mut() {
        *v /= n_orient as f64;
    }
    // Find the zero crossing of the descending loop by linear interpolation.
    for i in 0..n_fields - 1 {
        if mbar[i] >= 0.0 && mbar[i + 1] < 0.0 {
            let f = mbar[i] / (mbar[i] - mbar[i + 1]);
            return fields[i] + f * (fields[i + 1] - fields[i]);
        }
    }
    f64::NAN
}

/// Kittel single-domain critical diameter for a spherical particle
/// `D_c = 36 sqrt(A K)/(mu0 M_s^2)` (Kittel, *Phys. Rev.* **70**, 965
/// (1946)); an order-of-magnitude estimate.
pub fn single_domain_diameter(a: f64, k: f64, ms: f64) -> f64 {
    36.0 * (a * k).sqrt() / (crate::micromag::MU0 * ms * ms)
}

fn cross(a: V3, b: V3) -> V3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn add(a: V3, b: V3) -> V3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
fn scale(a: V3, s: f64) -> V3 {
    [a[0] * s, a[1] * s, a[2] * s]
}
