//! One-dimensional face fluxes along a line of cells: the hybrid scheme.
//!
//! Emergent-boundary logic. In smooth regions the solution is a band-
//! limited field and the cheapest accurate method is a non-dissipative
//! high-order central scheme (6th-order split form, KEEP-PE two-point
//! flux: no Riemann solver, no characteristic projection, no nonlinear
//! weights). Near shocks, contacts and material interfaces — sub-grid
//! structures whose interior the grid cannot represent, but whose jump
//! conditions are exact — it switches to a shock-capturing path:
//!
//! * `ShockScheme::Wcns` (default): weighted compact nonlinear scheme —
//!   characteristic WENO5-Z *interpolation* of point values to the face
//!   midpoints, an HLLC Riemann flux there, and the conservative 6th-order
//!   midpoint-to-node difference (Nonomura & Fujii 2013)
//!     ∂f/∂x|_i ≈ Σ_s c_s (F̂_{i+½+s} − F̂_{i−½−s})/h,  c = (75/64, −25/384, 3/640),
//!   written in flux form F̃_{i+½} = Σ_s c_s Σ_{|t|≤s} F̂_{i+½+t}, i.e. weights
//!   (c₀+c₁+c₂, c₁+c₂, c₂) on F̂_{i+½}, F̂_{i+½±1}, F̂_{i+½±2}.
//!   Riemann-solver based, so material interfaces keep pressure/velocity
//!   equilibrium; point-value consistent with the central scheme.
//! * `ShockScheme::FdWeno`: characteristic finite-difference WENO5-Z with
//!   per-field Lax–Friedrichs splitting (Jiang & Shu 1996), single
//!   material only.
//!
//! Switch: the local truncation-error estimate δ⁴q (fourth difference of the
//! system's sensor variables relative to the local level). In a smooth
//! field δ⁴q ≈ (kh)⁴q, so a face is "smooth" when the scheme's own error is
//! below `tol`; at a discontinuity δ⁴q is O(jump). The test uses only the
//! face's six-cell stencil, so every jump that can reach the face is seen.
//! The two fluxes are blended continuously across r = δ⁴/tol ∈ [½, 1], so
//! the scheme depends Lipschitz-continuously on the data.
//!
//! Every face also returns a face normal velocity Ũ, built with exactly the
//! same combination as the flux, so non-conservative terms c(w)∇·u
//! discretize as c_i (Ũ_{i+½} − Ũ_{i−½})/h consistently with the flux of
//! the advected variable (a uniform volume fraction stays uniform).

use crate::system::{invert, matmul, Mat, System, V};

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum FluxMode {
    /// sensor-switched central / shock-capturing
    Hybrid,
    /// shock-capturing path everywhere
    Weno,
    /// 6th-order central everywhere (smooth problems only)
    Central,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ShockScheme {
    Wcns,
    FdWeno,
}

const EPS: f64 = 1e-40;

/// WENO5-Z value at the right face of cell 2 from cell averages v[0..5]
/// (finite-volume reconstruction; used by the flux-split FD-WENO).
#[inline]
pub fn weno5z(v0: f64, v1: f64, v2: f64, v3: f64, v4: f64) -> f64 {
    let q0 = (2.0 * v0 - 7.0 * v1 + 11.0 * v2) / 6.0;
    let q1 = (-v1 + 5.0 * v2 + 2.0 * v3) / 6.0;
    let q2 = (2.0 * v2 + 5.0 * v3 - v4) / 6.0;
    weno_z_combine([q0, q1, q2], [0.1, 0.6, 0.3], v0, v1, v2, v3, v4)
}

/// WENO5-Z *interpolation* of point values to the midpoint x_{2+½}.
#[inline]
pub fn weno5z_interp(v0: f64, v1: f64, v2: f64, v3: f64, v4: f64) -> f64 {
    let q0 = (3.0 * v0 - 10.0 * v1 + 15.0 * v2) / 8.0;
    let q1 = (-v1 + 6.0 * v2 + 3.0 * v3) / 8.0;
    let q2 = (3.0 * v2 + 6.0 * v3 - v4) / 8.0;
    weno_z_combine([q0, q1, q2], [1.0 / 16.0, 10.0 / 16.0, 5.0 / 16.0], v0, v1, v2, v3, v4)
}

#[inline]
fn weno_z_combine(q: [f64; 3], d: [f64; 3], v0: f64, v1: f64, v2: f64, v3: f64, v4: f64) -> f64 {
    let b0 = 13.0 / 12.0 * (v0 - 2.0 * v1 + v2).powi(2) + 0.25 * (v0 - 4.0 * v1 + 3.0 * v2).powi(2);
    let b1 = 13.0 / 12.0 * (v1 - 2.0 * v2 + v3).powi(2) + 0.25 * (v1 - v3).powi(2);
    let b2 = 13.0 / 12.0 * (v2 - 2.0 * v3 + v4).powi(2) + 0.25 * (3.0 * v2 - 4.0 * v3 + v4).powi(2);
    let tau = (b0 - b2).abs();
    let a0 = d[0] * (1.0 + (tau / (b0 + EPS)).powi(2));
    let a1 = d[1] * (1.0 + (tau / (b1 + EPS)).powi(2));
    let a2 = d[2] * (1.0 + (tau / (b2 + EPS)).powi(2));
    (a0 * q[0] + a1 * q[1] + a2 * q[2]) / (a0 + a1 + a2)
}

/// 6th-order split-form coefficients a_l (central first derivative).
const A6: [f64; 3] = [0.75, -0.15, 1.0 / 60.0];
/// Flux-form weights of the 6th-order staggered difference
/// (c₀+c₁+c₂, c₁+c₂, c₂) with c = (75/64, −25/384, 3/640).
const S6: [f64; 3] = [75.0 / 64.0 - 25.0 / 384.0 + 3.0 / 640.0, -25.0 / 384.0 + 3.0 / 640.0, 3.0 / 640.0];

/// Fourth difference of a 5-cell window (of log-scaled sensor values, so a
/// relative measure).
#[inline]
fn d4(q0: f64, q1: f64, q2: f64, q3: f64, q4: f64) -> f64 {
    (q0 - 4.0 * q1 + 6.0 * q2 - 4.0 * q3 + q4).abs()
}

#[derive(Default, Clone, Copy, Debug)]
pub struct Stats {
    pub faces: usize,
    pub weno_faces: usize,
    pub fallback_faces: usize,
    /// faces where the positivity limiter blended in the first-order flux
    pub limited_faces: usize,
}

impl Stats {
    pub fn add(&mut self, o: &Stats) {
        self.faces += o.faces;
        self.weno_faces += o.weno_faces;
        self.fallback_faces += o.fallback_faces;
        self.limited_faces += o.limited_faces;
    }
}

#[derive(Clone, Copy, Debug)]
pub struct SchemeOpts {
    pub mode: FluxMode,
    pub shock: ShockScheme,
    pub tol: f64,
    /// D·Δt/Δx of the current (forward-Euler) stage for the positivity
    /// limiter (D = number of active dimensions); 0 disables it.
    pub lambda: f64,
}

/// Face fluxes F_{i+½} and face velocities Ũ_{i+½} for the n+1 faces
/// bounding the n interior cells of a line of primitive states `w`
/// (length n + 2·ng; ng ≥ 3, ≥ 5 for WCNS).
pub fn line_fluxes<const N: usize, S: System<N>>(sys: &S, w: &[V<N>], d: usize, ng: usize, o: SchemeOpts, out: &mut [V<N>], uface: &mut [f64], wface: &mut [V<N>], stats: &mut Stats) {
    let nv = N;
    let n = w.len() - 2 * ng;
    let un = sys.vel_index() + d;
    let mut ns = 0;
    let sens: Vec<[f64; 8]> = if o.mode == FluxMode::Hybrid {
        w.iter()
            .map(|wi| {
                let mut v = [0.0; 8];
                ns = sys.sensor(wi, &mut v);
                v
            })
            .collect()
    } else {
        Vec::new()
    };
    // midpoint Riemann fluxes for WCNS, computed on demand (index j = face j+½)
    let mut mid: Vec<Option<(V<N>, f64, V<N>)>> =
        if o.mode != FluxMode::Central && o.shock == ShockScheme::Wcns {
            vec![None; w.len()]
        } else {
            Vec::new()
        };
    // A lag-L pair occurs in L adjacent faces. Keep the last L pairs in
    // a stack ring; tags also handle gaps where the shock path is used.
    // Reuse the flux value without changing its accumulation order.
    let mut pair_flux = [[[0.0; N]; 3]; 2];
    let mut pair_index = [[usize::MAX; 3]; 2];
    let ncg = sys.nc_general();
    let mut l: Mat<N> = [[0.0; N]; N];
    let mut r: Mat<N> = [[0.0; N]; N];
    for f in 0..=n {
        let i = ng - 1 + f; // face between cells i and i+1
        // blending weight θ ∈ [0, 1] of the shock-capturing flux: a smooth
        // function of the sensor ratio r = δ⁴/tol (θ = 0 for r ≤ ½, 1 for
        // r ≥ 1, smoothstep between). A hard switch makes the scheme
        // discontinuous in the data: tiny perturbations flip faces between
        // schemes, which defeats error bounds (measured: an ε-independent
        // floor in the adaptive-vs-uniform difference).
        let theta = match o.mode {
            FluxMode::Central => 0.0,
            FluxMode::Weno => 1.0,
            FluxMode::Hybrid => {
                let mut r: f64 = 0.0;
                for k in 0..ns {
                    let s = |j: usize| sens[j][k];
                    r = r.max(d4(s(i - 2), s(i - 1), s(i), s(i + 1), s(i + 2))).max(d4(s(i - 1), s(i), s(i + 1), s(i + 2), s(i + 3)));
                }
                let x = ((r / o.tol - 0.5) / 0.5).clamp(0.0, 1.0);
                if r.is_finite() { x * x * (3.0 - 2.0 * x) } else { 1.0 }
            }
        };
        stats.faces += 1;
        let mut central = [0.0; N];
        let mut uf_c = 0.0;
        let mut wf_c = [0.0; N];
        if theta < 1.0 {
            for (li, &a) in A6.iter().enumerate() {
                let lag = li + 1;
                for m in 0..lag {
                    let (wa, wb) = (&w[i - m], &w[i - m + lag]);
                    let ft = if lag == 1 {
                        sys.split_flux(wa, wb, d)
                    } else {
                        let left = i - m;
                        let slot = if lag == 2 { left & 1 } else { left % 3 };
                        if pair_index[lag - 2][slot] != left {
                            pair_flux[lag - 2][slot] = sys.split_flux(wa, wb, d);
                            pair_index[lag - 2][slot] = left;
                        }
                        pair_flux[lag - 2][slot]
                    };
                    for k in 0..nv {
                        central[k] += 2.0 * a * ft[k];
                    }
                    uf_c += a * (wa[un] + wb[un]);
                    if ncg {
                        for k in 0..nv {
                            wf_c[k] += a * (wa[k] + wb[k]);
                        }
                    }
                }
            }
            if theta == 0.0 {
                out[f] = central;
                uface[f] = uf_c;
                if ncg {
                    wface[f] = wf_c;
                }
                continue;
            }
        }
        stats.weno_faces += 1;
        match o.shock {
            ShockScheme::Wcns => {
                let mut fl = [0.0; N];
                let mut uf = 0.0;
                let mut wf = [0.0; N];
                // F̃ = Σ_s S6[s] (F̂_{i+½+s} + F̂_{i+½−s}), the s = 0 term once
                for (s, &c) in S6.iter().enumerate() {
                    let js: &[usize] = if s == 0 { &[i] } else { &[i + s, i - s] };
                    for &j in js {
                        if mid[j].is_none() {
                            let (fv, us, fb, wm) = midpoint_flux(sys, w, j, d, &mut l, &mut r);
                            stats.fallback_faces += fb as usize;
                            mid[j] = Some((fv, us, wm));
                        }
                        let (fv, us, wm) = mid[j].unwrap();
                        for k in 0..nv {
                            fl[k] += c * fv[k];
                            wf[k] += c * wm[k];
                        }
                        uf += c * us;
                    }
                }
                out[f] = fl;
                uface[f] = uf;
                wface[f] = wf;
            }
            ShockScheme::FdWeno => {
                out[f] = fd_weno_flux(sys, w, i, d, &mut l, &mut r);
                for k in 0..nv {
                    wface[f][k] = 0.0;
                    for (li, &a) in A6.iter().enumerate() {
                        for m in 0..=li {
                            wface[f][k] += a * (w[i - m][k] + w[i - m + li + 1][k]);
                        }
                    }
                }
                let mut uf = 0.0;
                for (li, &a) in A6.iter().enumerate() {
                    for m in 0..=li {
                        uf += a * (w[i - m][un] + w[i - m + li + 1][un]);
                    }
                }
                uface[f] = uf;
            }
        }
        if theta < 1.0 {
            for k in 0..nv {
                out[f][k] = theta * out[f][k] + (1.0 - theta) * central[k];
                wface[f][k] = theta * wface[f][k] + (1.0 - theta) * wf_c[k];
            }
            uface[f] = theta * uface[f] + (1.0 - theta) * uf_c;
        }
    }
    if o.lambda > 0.0 {
        positivity_limit(sys, w, d, ng, o.lambda, out, uface, wface, stats);
    }
}

/// Characteristic WENO5-Z interpolation of the primitive variables to the
/// midpoint between cells j and j+1, and the Riemann flux there. Returns
/// (flux, face velocity, fell back to first order).
fn midpoint_flux<const N: usize, S: System<N>>(sys: &S, w: &[V<N>], j: usize, d: usize, l: &mut Mat<N>, r: &mut Mat<N>) -> (V<N>, f64, bool, V<N>) {
    let nv = N;
    let mut avg = [0.0; N];
    for k in 0..nv {
        avg[k] = 0.5 * (w[j][k] + w[j + 1][k]);
    }
    sys.eigen(&avg, d, l, r);
    let mut ch = [[0.0; N]; 6];
    for (s, cs) in ch.iter_mut().enumerate() {
        let ws = &w[j - 2 + s];
        for a in 0..nv {
            let mut acc = 0.0;
            for k in 0..nv {
                acc += l[a][k] * ws[k];
            }
            cs[a] = acc;
        }
    }
    let mut cl = [0.0; N];
    let mut cr = [0.0; N];
    for a in 0..nv {
        cl[a] = weno5z_interp(ch[0][a], ch[1][a], ch[2][a], ch[3][a], ch[4][a]);
        cr[a] = weno5z_interp(ch[5][a], ch[4][a], ch[3][a], ch[2][a], ch[1][a]);
    }
    let mut wl = [0.0; N];
    let mut wr = [0.0; N];
    for k in 0..nv {
        let (mut sl, mut sr) = (0.0, 0.0);
        for a in 0..nv {
            sl += r[k][a] * cl[a];
            sr += r[k][a] * cr[a];
        }
        wl[k] = sl;
        wr[k] = sr;
    }
    let mut fb = false;
    if !sys.admissible(&wl) || !sys.admissible(&wr) {
        fb = true;
        wl = w[j];
        wr = w[j + 1];
    }
    let (fv, us) = sys.riemann_flux(&wl, &wr, d);
    let mut wm = [0.0; N];
    for k in 0..nv {
        wm[k] = 0.5 * (wl[k] + wr[k]);
    }
    (fv, us, fb, wm)
}

/// Characteristic finite-difference WENO5-Z flux with per-field local
/// Lax–Friedrichs splitting at the face between cells i and i+1.
fn fd_weno_flux<const N: usize, S: System<N>>(sys: &S, w: &[V<N>], i: usize, d: usize, l: &mut Mat<N>, r: &mut Mat<N>) -> V<N> {
    let nv = N;
    let mut avg = [0.0; N];
    for k in 0..nv {
        avg[k] = 0.5 * (w[i][k] + w[i + 1][k]);
    }
    sys.eigen(&avg, d, l, r);
    let jw = sys.dudw(&avg);
    let lc = matmul(l, &invert(&jw));
    let rc = matmul(&jw, r);
    let mut alpha = [0.0f64; N];
    let mut gs = [[0.0; N]; 6];
    let mut vs = [[0.0; N]; 6];
    for s in 0..6 {
        let ws = &w[i - 2 + s];
        let q = sys.to_cons(ws);
        let fx = sys.flux(ws, d);
        for a in 0..nv {
            let (mut g, mut v) = (0.0, 0.0);
            for k in 0..nv {
                g += lc[a][k] * fx[k];
                v += lc[a][k] * q[k];
            }
            gs[s][a] = g;
            vs[s][a] = v;
        }
        let lam = sys.eigenvalues(ws, d);
        for a in 0..nv {
            alpha[a] = alpha[a].max(lam[a].abs());
        }
    }
    let mut fc = [0.0; N];
    for a in 0..nv {
        let fp = |s: usize| 0.5 * (gs[s][a] + alpha[a] * vs[s][a]);
        let fm = |s: usize| 0.5 * (gs[s][a] - alpha[a] * vs[s][a]);
        fc[a] = weno5z(fp(0), fp(1), fp(2), fp(3), fp(4)) + weno5z(fm(5), fm(4), fm(3), fm(2), fm(1));
    }
    let mut fl = [0.0; N];
    for k in 0..nv {
        for a in 0..nv {
            fl[k] += rc[k][a] * fc[a];
        }
    }
    if fl[..nv].iter().any(|v| !v.is_finite()) {
        fl = sys.riemann_flux(&w[i], &w[i + 1], d).0;
    }
    fl
}

/// Positivity-preserving flux limiter (Hu, Adams & Shu 2013). The update of
/// cell i is the mean of two one-sided updates,
///   U_i − 2λ(F_{i+½} − f(U_i))   and   U_i + 2λ(F_{i−½} − f(U_i)),
/// which are admissible for the first-order Lax–Friedrichs flux F_L when
/// 2λ·max|λ_k| ≤ 1. At each face the flux is replaced by θF + (1 − θ)F_L
/// with the largest θ ∈ [0, 1] keeping both adjacent one-sided states
/// admissible (the admissible set is convex in conserved variables, so a
/// bisection on θ is exact up to its tolerance). Untouched where the
/// high-order flux is already safe.
fn positivity_limit<const N: usize, S: System<N>>(sys: &S, w: &[V<N>], d: usize, ng: usize, lambda: f64, out: &mut [V<N>], uface: &mut [f64], wface: &mut [V<N>], stats: &mut Stats) {
    let nv = N;
    let n = w.len() - 2 * ng;
    let mut previous_cons = sys.to_cons(&w[ng - 1]);
    let mut previous_flux = sys.flux(&w[ng - 1], d);
    let ncg = sys.nc_general();
    for f in 0..=n {
        let i = ng - 1 + f;
        let (ul, fl_phys) = (previous_cons, previous_flux);
        let ur = sys.to_cons(&w[i + 1]);
        let fr_phys = sys.flux(&w[i + 1], d);
        previous_cons = ur;
        previous_flux = fr_phys;
        let ok = |fx: &V<N>| -> bool {
            let mut a = [0.0; N];
            let mut b = [0.0; N];
            for k in 0..nv {
                a[k] = ul[k] - 2.0 * lambda * (fx[k] - fl_phys[k]);
                b[k] = ur[k] + 2.0 * lambda * (fx[k] - fr_phys[k]);
            }
            sys.admissible(&sys.to_prim(&a)) && sys.admissible(&sys.to_prim(&b))
        };
        let fh = out[f];
        if fh[..nv].iter().all(|v| v.is_finite()) && ok(&fh) {
            continue;
        }
        // first-order Lax–Friedrichs flux
        let alpha = sys.max_speed(&w[i], d).max(sys.max_speed(&w[i + 1], d));
        let mut flf = [0.0; N];
        for k in 0..nv {
            flf[k] = 0.5 * (fl_phys[k] + fr_phys[k]) - 0.5 * alpha * (ur[k] - ul[k]);
        }
        stats.limited_faces += 1;
        let blend = |t: f64| -> V<N> {
            let mut v = [0.0; N];
            for k in 0..nv {
                v[k] = if fh[k].is_finite() { t * fh[k] + (1.0 - t) * flf[k] } else { flf[k] };
            }
            v
        };
        let (mut lo, mut hi) = (0.0, 1.0);
        for _ in 0..30 {
            let mid = 0.5 * (lo + hi);
            if ok(&blend(mid)) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        out[f] = blend(lo);
        // the face velocity of the non-conservative term follows the same
        // blend (first-order partner: the mean velocity, which is what the
        // Lax–Friedrichs flux of an advected fraction uses)
        let un = sys.vel_index() + d;
        let um = 0.5 * (w[i][un] + w[i + 1][un]);
        uface[f] = if uface[f].is_finite() { lo * uface[f] + (1.0 - lo) * um } else { um };
        // General non-conservative products must use the same limited
        // face interpolation as the flux and normal velocity.
        if ncg {
            for k in 0..nv {
                let wm = 0.5 * (w[i][k] + w[i + 1][k]);
                wface[f][k] = if wface[f][k].is_finite() {
                    lo * wface[f][k] + (1.0 - lo) * wm
                } else {
                    wm
                };
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        eos::Nasg,
        gpr::{Gpr, TauModel, NG},
    };

    #[test]
    fn limited_solid_flux_also_limits_the_nonconservative_face_state() {
        let sys = Gpr {
            eos: Nasg::ideal(1.4),
            cs: 1.0,
            tau: f64::INFINITY,
            tau_model: TauModel::Constant,
        };
        let mut wi = [0.0; NG];
        wi[0] = 1.0;
        wi[13] = 1.0;
        for k in 0..3 {
            wi[4 + 3 * k + k] = 1.0;
        }
        let w = [wi; 3]; // one interior cell and its two neighbors
        let physical = sys.flux(&wi, 0);
        let mut flux = [physical; 2];
        flux[0][13] += 100.0; // would drive a one-sided state negative
        let mut velocity = [3.0, 0.0];
        let high_face = std::array::from_fn(|k| wi[k] + 2.0);
        let mut face = [high_face, wi];
        face[0][12] = f64::NAN; // a failed reconstruction must also recover
        let mut stats = Stats::default();
        positivity_limit(
            &sys,
            &w,
            0,
            1,
            0.1,
            &mut flux,
            &mut velocity,
            &mut face,
            &mut stats,
        );
        assert_eq!(stats.limited_faces, 1);
        let theta = (flux[0][13] - physical[13]) / 100.0;
        assert!(theta > 0.0 && theta < 1.0);
        assert!((velocity[0] - 3.0 * theta).abs() < 1e-14);
        for k in 0..NG {
            let expected = if k == 12 {
                wi[k]
            } else {
                theta * high_face[k] + (1.0 - theta) * wi[k]
            };
            assert!((face[0][k] - expected).abs() < 1e-14);
        }
        // The actual one-sided conservative updates are admissible too.
        let u = sys.to_cons(&wi);
        for sign in [-1.0, 1.0] {
            let next = std::array::from_fn(|k| u[k] + sign * 0.2 * (flux[0][k] - physical[k]));
            assert!(sys.admissible(&sys.to_prim(&next)));
        }
        assert_eq!(face[1], wi);
    }
}
