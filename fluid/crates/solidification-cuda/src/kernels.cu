// CUDA (NVRTC) port of the Karma-Rappel anisotropic phase-field
// solidification step. Mirrors fluid/crates/solidification/src/phase_field.rs
// expression for expression so that f64 on the GPU reproduces the CPU to
// round-off. Templated on the real type R (float/double).
//
// Layout (same as the CPU): row-major with a one-cell ghost border,
//   idx(i,j) = j*(nx+2) + i,   interior i=1..=nx, j=1..=ny.
//
// BC edge codes: 0 = symmetry (zero normal derivative), 1 = far field
// (p=-1, u=u_far).

#ifndef PF_PI
#define PF_PI 3.14159265358979323846
#endif

// ---------------------------------------------------------- boundary conditions
// Left/right edges of every row (including the ghost rows). Uses the state's
// own ghosts, so it must run before the bottom/top pass that reads row 1.
template <typename R>
__global__ void pf_bc_lr_impl(int nx, int ny, int bc_left, int bc_right, R ufar, R* p, R* u) {
    int m = nx + 2;
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (j >= ny + 2) return;
    long base = (long)j * m;
    if (bc_left == 0) {
        p[base] = p[base + 1];
        u[base] = u[base + 1];
    } else {
        p[base] = R(-1);
        u[base] = ufar;
    }
    if (bc_right == 0) {
        p[base + nx + 1] = p[base + nx];
        u[base + nx + 1] = u[base + nx];
    } else {
        p[base + nx + 1] = R(-1);
        u[base + nx + 1] = ufar;
    }
}

// Bottom/top edges (overwrite the corners, as the CPU does).
template <typename R>
__global__ void pf_bc_bt_impl(int nx, int ny, int bc_bot, int bc_top, R ufar, R* p, R* u) {
    int m = nx + 2;
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= m) return;
    if (bc_bot == 0) {
        p[i] = p[m + i];
        u[i] = u[m + i];
    } else {
        p[i] = R(-1);
        u[i] = ufar;
    }
    long top = (long)ny * m + i;
    if (bc_top == 0) {
        p[top + m] = p[top];
        u[top + m] = u[top];
    } else {
        p[top + m] = R(-1);
        u[top + m] = ufar;
    }
}

// ---------------------------------------------------------------- phase step
template <typename R>
__global__ void pf_step_impl(int nx, int ny, R dx, R inv_twodx, R inv_fourd2, R w0, R eps4,
                             R lambda, R tau0, R d, R dt, const R* __restrict__ p,
                             const R* __restrict__ u, R* __restrict__ pn, R* __restrict__ un) {
    int m = nx + 2;
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int nc = nx * ny;
    if (idx >= nc) return;
    int i = idx % nx + 1;
    int j = idx / nx + 1;
    long k = (long)j * m + i;

    R dx2 = dx * dx;
    R px = (p[k + 1] - p[k - 1]) * inv_twodx;
    R py = (p[k + m] - p[k - m]) * inv_twodx;
    R pxx = (p[k + 1] - R(2) * p[k] + p[k - 1]) / dx2;
    R pyy = (p[k + m] - R(2) * p[k] + p[k - m]) / dx2;
    R pxy = (p[k + m + 1] - p[k + m - 1] - p[k - m + 1] + p[k - m - 1]) * inv_fourd2;
    R g2 = px * px + py * py;
    R one_minus = R(1) - p[k] * p[k];
    R drive = (p[k] - lambda * u[k] * one_minus) * one_minus;

    R lap_anis, tau;
    if (g2 > R(1e-18)) {
        R theta = atan2(py, px);
        R c4 = cos(R(4) * theta);
        R s4 = sin(R(4) * theta);
        R w = w0 * (R(1) + eps4 * c4);
        R dw = R(-4) * w0 * eps4 * s4;
        R theta_x = (px * pxy - py * pxx) / g2;
        R theta_y = (px * pyy - py * pxy) / g2;
        lap_anis = w * w * (pxx + pyy) + R(2) * w * dw * (theta_x * px + theta_y * py);
        R ratio = w / w0;
        tau = tau0 * ratio * ratio;
    } else {
        lap_anis = w0 * w0 * (pxx + pyy);
        tau = tau0;
    }
    R ds = (lap_anis + drive) / tau;
    R dp = dt * ds;
    if (!isfinite(dp)) dp = R(0);
    pn[k] = p[k] + dp;

    R lap_u = (u[k + 1] - R(2) * u[k] + u[k - 1]) / dx2
            + (u[k + m] - R(2) * u[k] + u[k - m]) / dx2;
    un[k] = u[k] + d * lap_u * dt + R(0.5) * dp;
}

// ------------------------------------------------- corrected Karma-Rappel
// KR 1998 Eq. (88) has the variational anisotropy term in addition to
// div(W^2 grad p):
//   + d_x( |grad p|^2 W dW/dp_x ) + d_y( |grad p|^2 W dW/dp_y )
// which, with Q = W dW/dtheta, reduces to  p_x d_y Q - p_y d_x Q.
// `solidification/src/phase_field.rs` omits this term, so it does not
// reproduce KR's selected dendrite; `pf_step_kr` restores it.
template <typename R>
__global__ void pf_q_impl(int nx, int ny, R inv_twodx, R w0, R eps4,
                          const R* __restrict__ p, R* __restrict__ q) {
    int m = nx + 2;
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int nc = nx * ny;
    if (idx >= nc) return;
    int i = idx % nx + 1;
    int j = idx / nx + 1;
    long k = (long)j * m + i;
    R px = (p[k + 1] - p[k - 1]) * inv_twodx;
    R py = (p[k + m] - p[k - m]) * inv_twodx;
    R g2 = px * px + py * py;
    R qv = R(0);
    if (g2 > R(1e-18)) {
        R g2i = R(1) / g2;
        R c4 = (px * px * px * px - R(6) * px * px * py * py + py * py * py * py) * g2i * g2i;
        R s4 = R(4) * px * py * (px * px - py * py) * g2i * g2i;
        R w = w0 * (R(1) + eps4 * c4);
        R dw = R(-4) * w0 * eps4 * s4;
        qv = w * dw;
    }
    q[k] = qv;
}

// One corrected step. `q` is dW/dtheta * W at every cell (from pf_q).
// `mode` selects the discretisation so the variants can be compared:
//   bit 0 -> 9-point isotropic Laplacian (else the 5-point one)
//   bit 1 -> include the variational anisotropy term p_x d_y Q - p_y d_x Q
// The KR production operator is mode = 3.
template <typename R>
__global__ void pf_step_kr_impl(int nx, int ny, R dx, R inv_twodx, R inv_fourd2, R w0, R eps4,
                                R lambda, R tau0, R d, R dt, int mode, const R* __restrict__ p,
                                const R* __restrict__ u, const R* __restrict__ q,
                                R* __restrict__ pn, R* __restrict__ un) {
    int m = nx + 2;
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int nc = nx * ny;
    if (idx >= nc) return;
    int i = idx % nx + 1;
    int j = idx / nx + 1;
    long k = (long)j * m + i;

    R dx2 = dx * dx;
    R px = (p[k + 1] - p[k - 1]) * inv_twodx;
    R py = (p[k + m] - p[k - m]) * inv_twodx;
    R pxx = (p[k + 1] - R(2) * p[k] + p[k - 1]) / dx2;
    R pyy = (p[k + m] - R(2) * p[k] + p[k - m]) / dx2;
    R pxy = (p[k + m + 1] - p[k + m - 1] - p[k - m + 1] + p[k - m - 1]) * inv_fourd2;
    R g2 = px * px + py * py;
    R one_minus = R(1) - p[k] * p[k];
    R drive = (p[k] - lambda * u[k] * one_minus) * one_minus;
    // Isotropic 9-point Laplacian (KR 1998 Sec. III), which suppresses the
    // grid anisotropy of the 5-point stencil.
    R inv_sixdx2 = R(1) / (R(6) * dx2);
    R lap9 = (R(4) * (p[k + 1] + p[k - 1] + p[k + m] + p[k - m])
              + (p[k + m + 1] + p[k + m - 1] + p[k - m + 1] + p[k - m - 1])
              - R(20) * p[k]) * inv_sixdx2;
    R lap = (mode & 1) ? lap9 : (pxx + pyy);

    R lap_anis, tau;
    if (g2 > R(1e-18)) {
        R g2i = R(1) / g2;
        R c4 = (px * px * px * px - R(6) * px * px * py * py + py * py * py * py) * g2i * g2i;
        R s4 = R(4) * px * py * (px * px - py * py) * g2i * g2i;
        R w = w0 * (R(1) + eps4 * c4);
        R dw = R(-4) * w0 * eps4 * s4;
        R theta_x = (px * pxy - py * pxx) / g2;
        R theta_y = (px * pyy - py * pxy) / g2;
        // base div(W^2 grad p)
        R base = w * w * lap + R(2) * w * dw * (theta_x * px + theta_y * py);
        lap_anis = base;
        // variational anisotropy: p_x d_y Q - p_y d_x Q, Q = W dW/dtheta
        if (mode & 2) {
            R qy = (q[k + m] - q[k - m]) * inv_twodx;
            R qx = (q[k + 1] - q[k - 1]) * inv_twodx;
            lap_anis += px * qy - py * qx;
        }
        R ratio = w / w0;
        tau = tau0 * ratio * ratio;
    } else {
        lap_anis = w0 * w0 * lap;
        tau = tau0;
    }
    R ds = (lap_anis + drive) / tau;
    R dp = dt * ds;
    if (!isfinite(dp)) dp = R(0);
    pn[k] = p[k] + dp;

    R lap_u = (u[k + 1] - R(2) * u[k] + u[k - 1]) / dx2
            + (u[k + m] - R(2) * u[k] + u[k - m]) / dx2;
    un[k] = u[k] + d * lap_u * dt + R(0.5) * dp;
}

// ------------------------------------------------------------------ entries
#define PF_ENTRIES(NAME, R)                                                                          \
    extern "C" __global__ void pf_bc_lr_##NAME(int nx, int ny, int bl, int br, R ufar, R* p, R* u) { \
        pf_bc_lr_impl<R>(nx, ny, bl, br, ufar, p, u);                                                \
    }                                                                                                \
    extern "C" __global__ void pf_bc_bt_##NAME(int nx, int ny, int bb, int bt, R ufar, R* p, R* u) { \
        pf_bc_bt_impl<R>(nx, ny, bb, bt, ufar, p, u);                                                \
    }                                                                                                \
    extern "C" __global__ void pf_step_##NAME(int nx, int ny, R dx, R inv_twodx, R inv_fourd2,      \
                                              R w0, R eps4, R lambda, R tau0, R d, R dt,             \
                                              const R* p, const R* u, R* pn, R* un) {                \
        pf_step_impl<R>(nx, ny, dx, inv_twodx, inv_fourd2, w0, eps4, lambda, tau0, d, dt, p, u, pn, un); \
    }                                                                                                \
    extern "C" __global__ void pf_step_kr_##NAME(int nx, int ny, R dx, R inv_twodx, R inv_fourd2,    \
                                                 R w0, R eps4, R lambda, R tau0, R d, R dt, int mode, \
                                                 const R* p, const R* u, const R* q, R* pn, R* un) { \
        pf_step_kr_impl<R>(nx, ny, dx, inv_twodx, inv_fourd2, w0, eps4, lambda, tau0, d, dt, mode, p, u, q, pn, un); \
    }                                                                                                \
    extern "C" __global__ void pf_q_##NAME(int nx, int ny, R inv_twodx, R w0, R eps4, const R* p, R* q) { \
        pf_q_impl<R>(nx, ny, inv_twodx, w0, eps4, p, q);                                             \
    }

PF_ENTRIES(f64, double)
PF_ENTRIES(f32, float)
