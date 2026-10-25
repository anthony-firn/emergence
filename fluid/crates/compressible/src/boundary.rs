//! Non-reflecting (characteristic) boundary states.
//!
//! At an open boundary the state is decomposed into characteristic waves
//! at the boundary cell (primitive eigenvectors of the system, so any EOS
//! and the multiphase systems work): waves travelling out of the domain
//! keep their interior amplitude, waves travelling in take the far-field
//! value,
//!   w_ghost = w_far + Σ_{k outgoing} r_k l_k·(w_int − w_far).
//! For linear waves at normal incidence this is exactly non-reflecting;
//! for oblique and nonlinear waves it is the standard first-order
//! approximation (Thompson 1987; Hedstrom), far better than zero-gradient
//! extrapolation, which reflects part of every acoustic wave and lets the
//! far-field pressure drift.

use crate::system::{System, V};

/// Conserved ghost state at side `side` (0 low, 1 high) of axis `d` from the
/// interior boundary cell `u_int` and the far-field conserved state `u_far`.
pub fn characteristic_ghost<const N: usize, S: System<N>>(
    sys: &S,
    u_int: &V<N>,
    u_far: &V<N>,
    d: usize,
    side: usize,
) -> V<N> {
    characteristic_ghost_layer(sys, u_int, u_int, u_far, d, side, 0.0)
}

/// Ghost state `m` cells beyond the boundary: the outgoing amplitudes are
/// extrapolated linearly from the boundary cell `u_b` and its interior
/// neighbour `u_b1` (constant amplitudes would put a kink into the
/// outgoing wave; measured 3.4% reflection of a resolved 1D pulse),
/// incoming amplitudes are those of the far field.
pub fn characteristic_ghost_layer<const N: usize, S: System<N>>(
    sys: &S,
    u_b: &V<N>,
    u_b1: &V<N>,
    u_far: &V<N>,
    d: usize,
    side: usize,
    m: f64,
) -> V<N> {
    CharacteristicStencil::new(sys, u_b, u_b1, u_far, d, side).conserved(sys, m)
}

/// Prepared once per boundary pencil and RK stage, shared by its ghost layers.
/// Retain the original projection and accumulation order at every layer.
pub(crate) struct CharacteristicStencil<const N: usize> {
    wi: V<N>,
    wf: V<N>,
    w1: V<N>,
    l: crate::system::Mat<N>,
    r: crate::system::Mat<N>,
    outgoing: [bool; N],
}
impl<const N: usize> CharacteristicStencil<N> {
    pub(crate) fn new<S: System<N>>(
        sys: &S,
        u_b: &V<N>,
        u_b1: &V<N>,
        u_far: &V<N>,
        d: usize,
        side: usize,
    ) -> Self {
        let (wi, wf, w1) = (sys.to_prim(u_b), sys.to_prim(u_far), sys.to_prim(u_b1));
        let mut l = [[0.0; N]; N];
        let mut r = [[0.0; N]; N];
        sys.eigen(&wi, d, &mut l, &mut r);
        let lam = sys.eigenvalues(&wi, d);
        let outward = if side == 1 { 1.0 } else { -1.0 };
        Self {
            wi,
            wf,
            w1,
            l,
            r,
            outgoing: lam.map(|v| v * outward > 0.0),
        }
    }
    pub(crate) fn conserved<S: System<N>>(&self, sys: &S, m: f64) -> V<N> {
        let mut dw = [0.0; N];
        for k in 0..N {
            dw[k] = self.wi[k] - self.wf[k] + m * (self.wi[k] - self.w1[k]);
        }
        let mut wg = self.wf;
        for kf in 0..N {
            if self.outgoing[kf] {
                let amp: f64 = (0..N).map(|m| self.l[kf][m] * dw[m]).sum();
                for q in 0..N {
                    wg[q] += self.r[q][kf] * amp;
                }
            }
        }
        sys.to_cons(&wg)
    }
}

/// Resolve deep reflected images on short axes. Keep this uncommon path out
/// of the ordinary ghost-copy loop's code and register working set.
#[cold]
#[inline(never)]
pub(crate) fn reflected_ghost<const N: usize, S: System<N>>(
    sys: &S,
    states: &[V<N>],
    at: &impl Fn(usize) -> usize,
    mut x: isize,
    n: usize,
    ng: usize,
    d: usize,
    bc: &[crate::solver::Bc<N>; 2],
    characteristic: &[Option<CharacteristicStencil<N>>; 2],
) -> V<N> {
    use crate::solver::Bc;
    let mut flip = false;
    let mut v = loop {
        if x >= 0 && x < n as isize {
            break states[at(ng + x as usize)];
        }
        let side = (x >= n as isize) as usize;
        match bc[side] {
            Bc::Reflect => {
                x = if side == 0 {
                    -1 - x
                } else {
                    2 * n as isize - 1 - x
                };
                flip = !flip;
            }
            Bc::Periodic => break states[at(ng + x.rem_euclid(n as isize) as usize)],
            Bc::Outflow => break states[at(ng + if side == 0 { 0 } else { n - 1 })],
            Bc::Dirichlet(far) => break far,
            Bc::Characteristic(_) => {
                let layer = if side == 0 { -x } else { x - n as isize + 1 };
                break characteristic[side]
                    .as_ref()
                    .unwrap()
                    .conserved(sys, layer as f64);
            }
        }
    };
    if flip {
        v[sys.momentum_index(d)] *= -1.0;
    }
    v
}
