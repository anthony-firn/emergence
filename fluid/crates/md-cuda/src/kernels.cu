// GPU molecular dynamics for Lennard-Jones pair potentials (NVRTC).
//
// Mirrors md/crates/md/src/{neighbor,system,potential}.rs. Data is
// structure-of-arrays, atom-major: pos/vel/force are [n][3], the per-atom
// energy `pe` is [n] and the per-atom (unsymmetrised) virial `stress` is
// [n][6].
//
// The neighbour scheme is a plain cell list rebuilt every force evaluation
// (no Verlet skin): atoms are binned, a serial exclusive scan turns the
// per-cell counts into offsets, and the fill kernel scatters atom indices
// into the cells. The force kernel then visits the (up to) 27 neighbouring
// cells of each atom, applying the minimum-image convention. Like the CPU
// engine this is a *full* list (each pair is seen from both endpoints), so
// the per-atom energy is halved and the pair virial is halved when totalled.
//
// Templated on the real type R (float/double).

__device__ __forceinline__ float r_sqrt(float x) { return sqrtf(x); }
__device__ __forceinline__ double r_sqrt(double x) { return sqrt(x); }
__device__ __forceinline__ float r_floor(float x) { return floorf(x); }
__device__ __forceinline__ double r_floor(double x) { return floor(x); }
__device__ __forceinline__ float r_round(float x) { return roundf(x); }
__device__ __forceinline__ double r_round(double x) { return round(x); }

// Truncated-and-shifted / shifted-force / bare-truncated Lennard-Jones pair
// energy and radial force. scheme: 0 = shifted, 1 = shifted-force,
// 2 = bare (tail corrections added on the host).
template <typename R>
__device__ __forceinline__ R lj_energy(R r, R rc, R sigma, R eps, int scheme) {
    R sr = sigma / r;
    R s2 = sr * sr;
    R s6 = s2 * s2 * s2;
    R u = R(4) * eps * (s6 * s6 - s6);
    if (scheme != 2) {
        R src = sigma / rc;
        R c2 = src * src;
        R c6 = c2 * c2 * c2;
        R uc = R(4) * eps * (c6 * c6 - c6);
        u -= uc;
        if (scheme == 1) {
            R d1c = R(24) * eps * (c6 - R(2) * c6 * c6) / rc;
            u -= (r - rc) * d1c;
        }
    }
    return u;
}

template <typename R>
__device__ __forceinline__ R lj_force(R r, R rc, R sigma, R eps, int scheme) {
    R sr = sigma / r;
    R s2 = sr * sr;
    R s6 = s2 * s2 * s2;
    R d = R(24) * eps * (s6 - R(2) * s6 * s6) / r;
    if (scheme == 1) {
        R src = sigma / rc;
        R c2 = src * src;
        R c6 = c2 * c2 * c2;
        d -= R(24) * eps * (c6 - R(2) * c6 * c6) / rc;
    }
    return d;
}

__device__ __forceinline__ int wrap_cell(int v, int nc) {
    v %= nc;
    if (v < 0) v += nc;
    return v;
}

// Bin every atom: counts[c]++ and cell_of[i] = c.
#define CELL_COUNTS(NAME, R)                                                          \
    extern "C" __global__ void md_cell_counts_##NAME(const R* pos, int* counts,       \
                                                     int* cell_of, int n, int ncx,    \
                                                     int ncy, int ncz, R invbox) {    \
        int i = blockIdx.x * blockDim.x + threadIdx.x;                                \
        if (i >= n) return;                                                           \
        int ix = (int)r_floor(pos[3 * i + 0] * invbox * ncx);                   \
        int iy = (int)r_floor(pos[3 * i + 1] * invbox * ncy);                   \
        int iz = (int)r_floor(pos[3 * i + 2] * invbox * ncz);                   \
        ix = wrap_cell(ix, ncx);                                                      \
        iy = wrap_cell(iy, ncy);                                                      \
        iz = wrap_cell(iz, ncz);                                                      \
        int c = ix + ncx * (iy + ncy * iz);                                           \
        cell_of[i] = c;                                                               \
        atomicAdd(&counts[c], 1);                                                     \
    }

// Exclusive scan of the per-cell counts into starts. One block of 256
// threads loops over the cells in chunks: an inclusive Hillis-Steele scan
// in shared memory plus a running carry. (A serial scan is O(ncells) on a
// single thread and dominates at 256k atoms, where ncells ~ 20 000.)
#define CELL_SCAN(NAME, R)                                                            \
    extern "C" __global__ void md_cell_scan_##NAME(int* counts, int* starts, int nc) {\
        __shared__ int tmp[256];                                                      \
        int tid = threadIdx.x;                                                        \
        int carry = 0;                                                                \
        for (int base = 0; base < nc; base += 256) {                                  \
            int i = base + tid;                                                       \
            int v = (i < nc) ? counts[i] : 0;                                         \
            tmp[tid] = v;                                                             \
            __syncthreads();                                                          \
            for (int off = 1; off < 256; off <<= 1) {                                 \
                int t = (tid >= off) ? tmp[tid - off] : 0;                            \
                __syncthreads();                                                      \
                tmp[tid] += t;                                                        \
                __syncthreads();                                                      \
            }                                                                         \
            if (i < nc) starts[i] = carry + tmp[tid] - v;                             \
            carry += tmp[255];                                                        \
            __syncthreads();                                                          \
        }                                                                             \
    }

// Scatter atom indices into their cells using starts as a cursor.
#define CELL_FILL(NAME, R)                                                            \
    extern "C" __global__ void md_cell_fill_##NAME(const int* cell_of, int* starts,   \
                                                   int* order, int n) {               \
        int i = blockIdx.x * blockDim.x + threadIdx.x;                                \
        if (i >= n) return;                                                           \
        int c = cell_of[i];                                                           \
        int slot = atomicAdd(&starts[c], 1);                                          \
        order[slot] = i;                                                              \
    }

// Lennard-Jones forces, per-atom energy (half the pair sum) and per-atom
// pair virial, visiting the up-to-27 unique neighbouring cells.
#define FORCE(NAME, R)                                                                \
    extern "C" __global__ void md_force_##NAME(                                       \
        const R* pos, const int* starts, const int* counts, const int* order,        \
        R* force, R* pe, R* stress, int n, int ncx, int ncy, int ncz, R box, R rc,   \
        R rc2, R sigma, R eps, int scheme) {                                         \
        int i = blockIdx.x * blockDim.x + threadIdx.x;                                \
        if (i >= n) return;                                                           \
        R xi = pos[3 * i + 0], yi = pos[3 * i + 1], zi = pos[3 * i + 2];              \
        R invbox = R(1) / box;                                                        \
        int hx = wrap_cell((int)r_floor(xi * invbox * ncx), ncx);                      \
        int hy = wrap_cell((int)r_floor(yi * invbox * ncy), ncy);                      \
        int hz = wrap_cell((int)r_floor(zi * invbox * ncz), ncz);                      \
        R f0 = R(0), f1 = R(0), f2 = R(0), p = R(0);                                  \
        R sxx = R(0), syy = R(0), szz = R(0), sxy = R(0), sxz = R(0), syz = R(0);     \
        for (int dz = -1; dz <= 1; dz++) {                                            \
            if (ncz == 1 && dz != 0) continue;                                        \
            if (ncz == 2 && dz == 1) continue;                                        \
            int cz = wrap_cell(hz + dz, ncz);                                         \
            for (int dy = -1; dy <= 1; dy++) {                                        \
                if (ncy == 1 && dy != 0) continue;                                    \
                if (ncy == 2 && dy == 1) continue;                                    \
                int cy = wrap_cell(hy + dy, ncy);                                     \
                for (int dx = -1; dx <= 1; dx++) {                                    \
                    if (ncx == 1 && dx != 0) continue;                                \
                    if (ncx == 2 && dx == 1) continue;                                \
                    int cx = wrap_cell(hx + dx, ncx);                                 \
                    int c = cx + ncx * (cy + ncy * cz);                               \
                    int e = starts[c];                                                \
                    int b = e - counts[c];                                            \
                    for (int idx = b; idx < e; idx++) {                               \
                        int j = order[idx];                                           \
                        if (j == i) continue;                                         \
                        R d0 = pos[3 * j + 0] - xi;                                   \
                        R d1 = pos[3 * j + 1] - yi;                                   \
                        R d2 = pos[3 * j + 2] - zi;                                   \
                        d0 -= box * r_round(d0 * invbox);                             \
                        d1 -= box * r_round(d1 * invbox);                             \
                        d2 -= box * r_round(d2 * invbox);                             \
                        R r2 = d0 * d0 + d1 * d1 + d2 * d2;                           \
                        if (r2 > rc2) continue;                                       \
                        R r = r_sqrt(r2);                                             \
                        R u = lj_energy(r, rc, sigma, eps, scheme);                   \
                        R fmag = lj_force(r, rc, sigma, eps, scheme);                 \
                        R invr = R(1) / r;                                            \
                        R fj0 = fmag * d0 * invr;                                     \
                        R fj1 = fmag * d1 * invr;                                     \
                        R fj2 = fmag * d2 * invr;                                     \
                        f0 += fj0; f1 += fj1; f2 += fj2;                              \
                        p += u;                                                       \
                        R cc = -fmag * invr;                                          \
                        sxx += cc * d0 * d0;                                          \
                        syy += cc * d1 * d1;                                          \
                        szz += cc * d2 * d2;                                          \
                        sxy += cc * d0 * d1;                                          \
                        sxz += cc * d0 * d2;                                          \
                        syz += cc * d1 * d2;                                          \
                    }                                                                 \
                }                                                                     \
            }                                                                         \
        }                                                                             \
        force[3 * i + 0] = f0; force[3 * i + 1] = f1; force[3 * i + 2] = f2;          \
        pe[i] = R(0.5) * p;                                                           \
        stress[6 * i + 0] = sxx; stress[6 * i + 1] = syy; stress[6 * i + 2] = szz;    \
        stress[6 * i + 3] = sxy; stress[6 * i + 4] = sxz; stress[6 * i + 5] = syz;    \
    }

// Instantaneous pressure-tensor numerator, accumulated in *double* even
// when the force path is f32: out[off + k] += sum_i [m v_i v_j]
//                                               + 0.5 sum_i stress_i[k]
// for k = 0..5 (xx,yy,zz,xy,xz,yz). The per-atom pair virial `stress` is
// already summed over neighbours and is halved when totalled, exactly as
// in the host `pressure_tensor`. The host divides by V and adds the
// isotropic tail pressure on the diagonal.
#define PRESS(NAME, R)                                                                \
    extern "C" __global__ void md_press_##NAME(const R* vel, const R* stress,         \
                                               double* out, int n, int off) {         \
        int i = blockIdx.x * blockDim.x + threadIdx.x;                                \
        if (i >= n) return;                                                           \
        R vx = vel[3 * i + 0], vy = vel[3 * i + 1], vz = vel[3 * i + 2];              \
        double k0 = (double)vx * (double)vx;                                          \
        double k1 = (double)vy * (double)vy;                                          \
        double k2 = (double)vz * (double)vz;                                          \
        double k3 = (double)vx * (double)vy;                                          \
        double k4 = (double)vx * (double)vz;                                          \
        double k5 = (double)vy * (double)vz;                                          \
        atomicAdd(&out[off + 0], k0 + 0.5 * (double)stress[6 * i + 0]);               \
        atomicAdd(&out[off + 1], k1 + 0.5 * (double)stress[6 * i + 1]);               \
        atomicAdd(&out[off + 2], k2 + 0.5 * (double)stress[6 * i + 2]);               \
        atomicAdd(&out[off + 3], k3 + 0.5 * (double)stress[6 * i + 3]);               \
        atomicAdd(&out[off + 4], k4 + 0.5 * (double)stress[6 * i + 4]);               \
        atomicAdd(&out[off + 5], k5 + 0.5 * (double)stress[6 * i + 5]);               \
    }

// Multiply every velocity by a scalar `lambda` (global thermostat move).
#define SCALE(NAME, R)                                                                \
    extern "C" __global__ void md_scale_##NAME(R* vel, int n, R lambda) {             \
        int i = blockIdx.x * blockDim.x + threadIdx.x;                                \
        if (i >= n) return;                                                           \
        vel[3 * i + 0] *= lambda;                                                     \
        vel[3 * i + 1] *= lambda;                                                     \
        vel[3 * i + 2] *= lambda;                                                     \
    }

// Velocity half-kick: v += 0.5 dt f (unit mass).
#define KICK(NAME, R)                                                                 \
    extern "C" __global__ void md_kick_##NAME(R* vel, const R* force, int n, R c) {   \
        int i = blockIdx.x * blockDim.x + threadIdx.x;                                \
        if (i >= n) return;                                                           \
        vel[3 * i + 0] += c * force[3 * i + 0];                                       \
        vel[3 * i + 1] += c * force[3 * i + 1];                                       \
        vel[3 * i + 2] += c * force[3 * i + 2];                                       \
    }

// Position drift: x += dt v.
#define DRIFT(NAME, R)                                                                \
    extern "C" __global__ void md_drift_##NAME(R* pos, const R* vel, int n, R dt) {   \
        int i = blockIdx.x * blockDim.x + threadIdx.x;                                \
        if (i >= n) return;                                                           \
        pos[3 * i + 0] += dt * vel[3 * i + 0];                                        \
        pos[3 * i + 1] += dt * vel[3 * i + 1];                                        \
        pos[3 * i + 2] += dt * vel[3 * i + 2];                                        \
    }

// Global reductions needed for diagnostics. out = [pe, stress[6],
// kinetic_tensor[6]] (13 entries). Kinetic tensor is m v_i v_j with m = 1.
#define DIAG(NAME, R)                                                                 \
    extern "C" __global__ void md_diag_##NAME(const R* vel, const R* pe,              \
                                              const R* stress, R* out, int n) {       \
        int i = blockIdx.x * blockDim.x + threadIdx.x;                                \
        if (i >= n) return;                                                           \
        R vx = vel[3 * i + 0], vy = vel[3 * i + 1], vz = vel[3 * i + 2];              \
        atomicAdd(&out[0], pe[i]);                                                    \
        for (int k = 0; k < 6; k++) atomicAdd(&out[1 + k], stress[6 * i + k]);        \
        atomicAdd(&out[7], vx * vx);                                                  \
        atomicAdd(&out[8], vy * vy);                                                  \
        atomicAdd(&out[9], vz * vz);                                                  \
        atomicAdd(&out[10], vx * vy);                                                 \
        atomicAdd(&out[11], vx * vz);                                                 \
        atomicAdd(&out[12], vy * vz);                                                 \
    }

#define MD_DEF(NAME, R)                                                               \
    CELL_COUNTS(NAME, R)                                                              \
    CELL_SCAN(NAME, R)                                                                \
    CELL_FILL(NAME, R)                                                                \
    FORCE(NAME, R)                                                                    \
    KICK(NAME, R)                                                                     \
    DRIFT(NAME, R)                                                                    \
    DIAG(NAME, R)                                                                     \
    PRESS(NAME, R)                                                                    \
    SCALE(NAME, R)

MD_DEF(f64, double)
MD_DEF(f32, float)
