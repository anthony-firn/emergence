//! Metamaterials with a magnetic response: negative refraction and the perfect
//! lens.
//!
//! An ordinary dielectric has `mu = 1`, so its refractive index is
//! `n = sqrt(eps)`. A metamaterial prepared from conducting scatterers can have
//! both `Re eps < 0` and `Re mu < 0`. The physical branch of
//! `n = sqrt(eps mu)` then has `Re n < 0` (Veselago 1968; Pendry 2000): the
//! phase velocity is opposite the Poynting vector, so a slab focuses light and
//! reverses the phase of a transmitted wave. The split-ring resonator gives a
//! causal Lorentz resonance in `mu` (Pendry 1999); a wire medium gives a Drude
//! `eps`. Where the two negative regions overlap, the medium supports a
//! negative-index band.
//!
//! The transfer matrix in [`crate::multilayer`] is exact for such a slab once
//! the layer carries a complex `mu` and the normal wavevector uses the signed
//! branch `q = n cos theta` (`multilayer::normal_wavevector`). The checks here
//! are the Veselago–Pendry ones:
//!
//! * a lossless `eps = mu = -1` slab in vacuum is impedance matched, so
//!   `T = 1` and `R = 0` at every angle, and its transmission phase is the
//!   negative of an equal vacuum gap;
//! * for evanescent waves (`s = k_x/k0 > 1`) the ideal slab amplifies them by
//!   exactly `|t| = exp(sqrt(s^2 - 1) k0 d)`, Pendry's perfect lens;
//! * a small loss caps that amplification at a finite transverse wavevector,
//!   which is the resolution limit.
//!
//! Sources: V. G. Veselago, Sov. Phys. Usp. **10**, 509 (1968);
//! J. B. Pendry, Phys. Rev. Lett. **85**, 3966 (2000),
//! doi:10.1103/PhysRevLett.85.3966; J. B. Pendry et al., IEEE Trans. Microwave
//! Theory Tech. **47**, 2075 (1999), doi:10.1109/22.798002.
use crate::complex::C;
use crate::dispersion::{drude, SplitRing};
use crate::multilayer::{Layer, Polarization, Stack};
use std::f64::consts::PI;

/// A wire-medium (Drude) electric response together with a split-ring
/// (Lorentz) magnetic response. Angular frequencies are in rad/s.
#[derive(Clone, Copy, Debug)]
pub struct Metamaterial {
    /// High-frequency permittivity (usually 1).
    pub eps_inf: f64,
    /// Wire-medium plasma frequency.
    pub omega_p: f64,
    /// Electric damping.
    pub gamma_e: f64,
    /// Split-ring permeability.
    pub srr: SplitRing,
}

impl Metamaterial {
    pub fn new(eps_inf: f64, omega_p: f64, gamma_e: f64, srr: SplitRing) -> Self {
        assert!(eps_inf > 0.0 && omega_p > 0.0 && gamma_e >= 0.0);
        Self {
            eps_inf,
            omega_p,
            gamma_e,
            srr,
        }
    }

    /// Complex permittivity `eps(omega)`; the wire medium is the Drude form
    /// `eps = eps_inf - omega_p^2 / (omega (omega + i gamma_e))`.
    pub fn eps(&self, omega: f64) -> C {
        assert!(omega > 0.0);
        let d = drude(omega, self.omega_p, self.gamma_e);
        // Scale the free-electron part by (eps_inf - 1) if eps_inf != 1.
        if (self.eps_inf - 1.0).abs() < 1e-15 {
            d
        } else {
            C::real(self.eps_inf) + (d - C::real(1.0))
        }
    }

    /// Complex permeability `mu(omega)` from the split-ring resonance.
    pub fn mu(&self, omega: f64) -> C {
        self.srr.mu(omega)
    }

    /// Signed refractive index `n = sqrt(eps) sqrt(mu)` on the physical branch
    /// (`Re n < 0` when `Re eps < 0` and `Re mu < 0`, `Im n >= 0`).
    pub fn index(&self, omega: f64) -> C {
        self.eps(omega).sqrt() * self.mu(omega).sqrt()
    }

    /// Relative impedance `Z/Z0 = sqrt(mu/eps)`.
    pub fn impedance(&self, omega: f64) -> C {
        (self.mu(omega) / self.eps(omega)).sqrt()
    }

    pub fn at_lambda(&self, lambda_m: f64) -> (C, C, C) {
        let omega = 2.0 * PI * crate::C0 / lambda_m;
        (self.eps(omega), self.mu(omega), self.index(omega))
    }

    /// True where the real parts of both `eps` and `mu` are negative.
    pub fn is_negative_index(&self, omega: f64) -> bool {
        self.eps(omega).re < 0.0 && self.mu(omega).re < 0.0
    }

    /// Interval `(omega_lo, omega_hi)` on which the medium is double negative,
    /// found by a coarse scan followed by bisection on `Re eps` and `Re mu`.
    /// Returns `None` if the two negative regions do not overlap.
    pub fn negative_index_band(&self, omega_lo: f64, omega_hi: f64) -> Option<(f64, f64)> {
        assert!(omega_lo > 0.0 && omega_hi > omega_lo);
        let n = 400usize;
        let mut edges: Vec<f64> = Vec::new();
        let mut inside = self.is_negative_index(omega_lo);
        if inside {
            edges.push(omega_lo);
        }
        for i in 1..=n {
            let w = omega_lo + (omega_hi - omega_lo) * i as f64 / n as f64;
            let now = self.is_negative_index(w);
            if now != inside {
                // Bisect the crossing of the combined predicate.
                let (mut a, mut b) = (w - (omega_hi - omega_lo) / n as f64, w);
                for _ in 0..80 {
                    let m = 0.5 * (a + b);
                    if self.is_negative_index(m) == inside {
                        a = m;
                    } else {
                        b = m;
                    }
                }
                edges.push(0.5 * (a + b));
                inside = now;
            }
        }
        if inside {
            edges.push(omega_hi);
        }
        if edges.len() >= 2 {
            // Longest double-negative interval.
            let mut best = (edges[0], edges[1]);
            for c in edges.chunks(2) {
                if c.len() == 2 && c[1] - c[0] > best.1 - best.0 {
                    best = (c[0], c[1]);
                }
            }
            Some(best)
        } else {
            None
        }
    }

    /// A metamaterial layer at angular frequency `omega`.
    pub fn layer(&self, omega: f64, thickness: f64) -> Layer {
        Layer::from_eps_mu(self.eps(omega), self.mu(omega), thickness)
    }
}

/// Ideal lossless Veselago layer with `eps = mu = -1` and thickness `d`.
pub fn perfect_lens_layer(thickness: f64) -> Layer {
    Layer::from_eps_mu(C::real(-1.0), C::real(-1.0), thickness)
}

/// Analytic transmission amplitude of an ideal (`eps = mu = -1`) slab in
/// vacuum for the evanescent wave with tangential wavevector `s = k_x/k0`:
/// `t = exp(sqrt(s^2 - 1) k0 d)`. This is the perfect lens: the slab exactly
/// undoes the vacuum decay `exp(-sqrt(s^2-1) k0 d)` over the same distance.
/// For `s <= 1` no amplification occurs and this returns `1` (the slab is
/// impedance matched, `|t| = 1`).
pub fn perfect_lens_evanescent_t(lambda: f64, s: f64, thickness: f64) -> C {
    let kappa = (s * s - 1.0).max(0.0).sqrt() * 2.0 * PI / lambda;
    C::real((kappa * thickness).exp())
}

/// Transmission amplitude of a single `eps = mu = -1` slab in vacuum at
/// tangential wavevector `s`, from the transfer matrix. Convenience wrapper.
pub fn perfect_lens_slab_t(lambda: f64, s: f64, thickness: f64, pol: Polarization) -> C {
    Stack::single(1.0, perfect_lens_layer(thickness), 1.0)
        .solve_s(lambda, s, pol)
        .t
}

/// Peak evanescent amplification `(k_x, |t|_peak)` of a lossy negative-index
/// slab, where `k_x = s k0`. The scan runs `s` from 1 to `s_max` on `steps`
/// points. The return value is the resolution metric: features with transverse
/// wavevector beyond `k_x` are not amplified, so the ideal resolution is about
/// `2 pi / k_x`.
pub fn evanescent_amplification_peak(
    lambda: f64,
    thickness: f64,
    eps: C,
    mu: C,
    s_max: f64,
    steps: usize,
) -> (f64, f64) {
    let stack = Stack::single(1.0, Layer::from_eps_mu(eps, mu, thickness), 1.0);
    let mut best = (0.0, 0.0);
    for i in 0..=steps {
        let s = 1.0 + (s_max - 1.0) * i as f64 / steps as f64;
        let t = stack.solve_s(lambda, s, Polarization::TM).t.abs();
        if t > best.1 {
            best = (s, t);
        }
    }
    (best.0 * 2.0 * PI / lambda, best.1)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn physical_index_branch() {
        let e = C::real(-1.0);
        let n = e.sqrt() * e.sqrt();
        assert!((n.re + 1.0).abs() < 1e-14 && n.im.abs() < 1e-14);
        // Lossy: Re n < 0 and Im n >= 0.
        let e = C::new(-1.0, 0.01);
        let n = e.sqrt() * e.sqrt();
        assert!(n.re < 0.0 && n.im >= 0.0, "{n:?}");
    }
}
