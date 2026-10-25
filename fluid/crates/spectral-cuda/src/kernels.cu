// Kernels for the pseudo-spectral Navier–Stokes step (f32 storage).
// Layouts match spectral-core/src/grid.rs:
//   physical  N*N*N,       index (i*N + j)*N + l
//   spectral  N*N*(N/2+1), index (i*N + j)*NH + l
// cuFFT transforms are unnormalized; the forward 1/N^3 factor is applied in
// project_scale.


__device__ __forceinline__ float wavenumber(int i, int n) {
    return (float)(i <= n / 2 ? i : i - n);
}

__device__ __forceinline__ void kvec(long idx, int n, int nh, float &kx, float &ky, float &kz) {
    int l = (int)(idx % nh);
    long ij = idx / nh;
    int j = (int)(ij % n);
    int i = (int)(ij / n);
    kx = wavenumber(i, n);
    ky = wavenumber(j, n);
    kz = (float)l;
}

// Retained modes: 2/3 cube (kc2 == 0) or spherical truncation |k|^2 < kc2
// (phase-shift dealiasing); Nyquist modes always removed.
__device__ __forceinline__ bool kept(float kx, float ky, float kz, int n, float kc2) {
    float nyq = (float)(n / 2);
    if (fabsf(kx) >= nyq || fabsf(ky) >= nyq || fabsf(kz) >= nyq) return false;
    if (kc2 > 0.0f) return kx * kx + ky * ky + kz * kz < kc2;
    float lim = (float)n / 3.0f;
    return fabsf(kx) < lim && fabsf(ky) < lim && fabsf(kz) < lim;
}

// e^{i k.s}
__device__ __forceinline__ float2 shift_phase(float kx, float ky, float kz,
                                              float sx, float sy, float sz) {
    float sn, cs;
    sincosf(kx * sx + ky * sy + kz * sz, &sn, &cs);
    return make_float2(cs, sn);
}

__device__ __forceinline__ float2 cmul(float2 a, float2 b) {
    return make_float2(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
}

// dst = src * e^{i k.s}  (velocity evaluated on a grid shifted by s)
extern "C" __global__ void phase_copy(float2 *dst, const float2 *src, float sx, float sy,
                                      float sz, int n, int nh) {
    long idx = blockIdx.x * (long)blockDim.x + threadIdx.x;
    long ns = (long)n * n * nh;
    if (idx >= ns) return;
    float kx, ky, kz;
    kvec(idx, n, nh, kx, ky, kz);
    dst[idx] = cmul(src[idx], shift_phase(kx, ky, kz, sx, sy, sz));
}

// out += c * exp(-nu |k|^2 h) * x, applied as c*fma(expm1, x, x) (unbiased)
extern "C" __global__ void accum_if(float2 *out, const float2 *x, float c, float h, float nu,
                                    int n, int nh) {
    long idx = blockIdx.x * (long)blockDim.x + threadIdx.x;
    long ns = (long)n * n * nh;
    if (idx >= ns) return;
    float kx, ky, kz;
    kvec(idx, n, nh, kx, ky, kz);
    float em1 = expm1f(-nu * (kx * kx + ky * ky + kz * kz) * h);
    float2 v = x[idx], o = out[idx];
    o.x = fmaf(c, fmaf(em1, v.x, v.x), o.x);
    o.y = fmaf(c, fmaf(em1, v.y, v.y), o.y);
    out[idx] = o;
}

// w_hat = i k x u_hat (all three components per thread).
extern "C" __global__ void curl(const float2 *u0, const float2 *u1, const float2 *u2,
                                float2 *w0, float2 *w1, float2 *w2, int n, int nh,
                                float sx, float sy, float sz) {
    long idx = blockIdx.x * (long)blockDim.x + threadIdx.x;
    long ns = (long)n * n * nh;
    if (idx >= ns) return;
    float kx, ky, kz;
    kvec(idx, n, nh, kx, ky, kz);
    float2 a = u0[idx], b = u1[idx], c = u2[idx];
    // i*(z) = (-z.y, z.x)
    float2 cx = make_float2(ky * c.x - kz * b.x, ky * c.y - kz * b.y);
    float2 cy = make_float2(kz * a.x - kx * c.x, kz * a.y - kx * c.y);
    float2 cz = make_float2(kx * b.x - ky * a.x, kx * b.y - ky * a.y);
    float2 r0 = make_float2(-cx.y, cx.x), r1 = make_float2(-cy.y, cy.x), r2 = make_float2(-cz.y, cz.x);
    if (sx != 0.0f || sy != 0.0f || sz != 0.0f) {
        float2 ph = shift_phase(kx, ky, kz, sx, sy, sz);
        r0 = cmul(r0, ph); r1 = cmul(r1, ph); r2 = cmul(r2, ph);
    }
    w0[idx] = r0;
    w1[idx] = r1;
    w2[idx] = r2;
}

// In physical space: w <- u x w (in place, all components per thread).
extern "C" __global__ void cross_inplace(const float *u0, const float *u1, const float *u2,
                                         float *w0, float *w1, float *w2, long nr) {
    long p = blockIdx.x * (long)blockDim.x + threadIdx.x;
    if (p >= nr) return;
    float a0 = u0[p], a1 = u1[p], a2 = u2[p];
    float b0 = w0[p], b1 = w1[p], b2 = w2[p];
    w0[p] = a1 * b2 - a2 * b1;
    w1[p] = a2 * b0 - a0 * b2;
    w2[p] = a0 * b1 - a1 * b0;
}

// v <- mask * P(scale * v + f), P = I - k k^T/|k|^2; k = 0 mode zeroed.
// `has_f` selects whether forcing f0..f2 is added.
// keep_mean: leave the k = 0 mode (scaled) instead of zeroing it, and add the
// uniform body force (bf0, bf1, bf2) to it (free, driven mean flow).
extern "C" __global__ void project_scale(float2 *v0, float2 *v1, float2 *v2,
                                         const float2 *f0, const float2 *f1, const float2 *f2,
                                         int has_f, float scale, int n, int nh,
                                         int keep_mean, float bf0, float bf1, float bf2,
                                         float kc2, float sx, float sy, float sz) {
    long idx = blockIdx.x * (long)blockDim.x + threadIdx.x;
    long ns = (long)n * n * nh;
    if (idx >= ns) return;
    float kx, ky, kz;
    kvec(idx, n, nh, kx, ky, kz);
    float k2 = kx * kx + ky * ky + kz * kz;
    if (k2 == 0.0f && keep_mean) {
        float2 a = v0[idx], b = v1[idx], c = v2[idx];
        v0[idx] = make_float2(a.x * scale + bf0, 0.0f);
        v1[idx] = make_float2(b.x * scale + bf1, 0.0f);
        v2[idx] = make_float2(c.x * scale + bf2, 0.0f);
        return;
    }
    if (!kept(kx, ky, kz, n, kc2) || k2 == 0.0f) {
        v0[idx] = v1[idx] = v2[idx] = make_float2(0.0f, 0.0f);
        return;
    }
    float2 a = v0[idx], b = v1[idx], c = v2[idx];
    a.x *= scale; a.y *= scale; b.x *= scale; b.y *= scale; c.x *= scale; c.y *= scale;
    if (sx != 0.0f || sy != 0.0f || sz != 0.0f) {
        // undo the grid shift of the product: multiply by e^{-i k.s}
        float2 ph = shift_phase(kx, ky, kz, -sx, -sy, -sz);
        a = cmul(a, ph); b = cmul(b, ph); c = cmul(c, ph);
    }
    if (has_f) {
        float2 p = f0[idx], q = f1[idx], r = f2[idx];
        a.x += p.x; a.y += p.y; b.x += q.x; b.y += q.y; c.x += r.x; c.y += r.y;
    }
    float inv = 1.0f / k2;
    float dr = (kx * a.x + ky * b.x + kz * c.x) * inv;
    float di = (kx * a.y + ky * b.y + kz * c.y) * inv;
    v0[idx] = make_float2(a.x - dr * kx, a.y - di * kx);
    v1[idx] = make_float2(b.x - dr * ky, b.y - di * ky);
    v2[idx] = make_float2(c.x - dr * kz, c.y - di * kz);
}

// IF-RK4 stage combination for one component:
//   out = sum_m c[m] * exp(-nu |k|^2 h[m]) * x_m,   m = 0..4
// Terms with c[m] == 0 are skipped (their pointer may alias anything).
// `out` may alias x0 (each element reads only its own index).
extern "C" __global__ void combine5(float2 *out,
                                    const float2 *x0, const float2 *x1, const float2 *x2,
                                    const float2 *x3, const float2 *x4,
                                    float c0, float c1, float c2, float c3, float c4,
                                    float h0, float h1, float h2, float h3, float h4,
                                    float nu, int n, int nh) {
    long idx = blockIdx.x * (long)blockDim.x + threadIdx.x;
    long ns = (long)n * n * nh;
    if (idx >= ns) return;
    float kx, ky, kz;
    kvec(idx, n, nh, kx, ky, kz);
    float nk2 = -nu * (kx * kx + ky * ky + kz * kz);
    float2 acc = make_float2(0.0f, 0.0f);
    const float2 *xs[5] = {x0, x1, x2, x3, x4};
    float cs[5] = {c0, c1, c2, c3, c4};
    float hs[5] = {h0, h1, h2, h3, h4};
    // E = 1 + expm1(x) applied as fma(expm1, v, v): at low viscosity the
    // per-step decrement is ~1e-5, which a rounded f32 factor e^x would carry
    // with ~1e-3 relative error, identically every step (a systematic
    // dissipation bias). This form rounds once per step and stays unbiased.
#pragma unroll
    for (int m = 0; m < 5; ++m) {
        if (cs[m] != 0.0f) {
            float em1 = expm1f(nk2 * hs[m]);
            float2 v = xs[m][idx];
            acc.x = fmaf(cs[m], fmaf(em1, v.x, v.x), acc.x);
            acc.y = fmaf(cs[m], fmaf(em1, v.y, v.y), acc.y);
        }
    }
    out[idx] = acc;
}

// Block partial sums (in double) of energy, enstrophy and max divergence.
// out has 3 doubles per block: [sumE, sumZ, maxDiv].
extern "C" __global__ void diag_partial(const float2 *u0, const float2 *u1, const float2 *u2,
                                        double *out, int n, int nh) {
    extern __shared__ double sh[];
    double *se = sh, *sz = sh + blockDim.x, *sd = sh + 2 * blockDim.x;
    long ns = (long)n * n * nh;
    double e = 0.0, z = 0.0, d = 0.0;
    for (long idx = blockIdx.x * (long)blockDim.x + threadIdx.x; idx < ns;
         idx += (long)gridDim.x * blockDim.x) {
        float kx, ky, kz;
        kvec(idx, n, nh, kx, ky, kz);
        int l = (int)(idx % nh);
        double w = (l == 0 || l == n / 2) ? 1.0 : 2.0;
        float2 a = u0[idx], b = u1[idx], c = u2[idx];
        double m2 = (double)a.x * a.x + (double)a.y * a.y + (double)b.x * b.x +
                    (double)b.y * b.y + (double)c.x * c.x + (double)c.y * c.y;
        double k2 = (double)kx * kx + (double)ky * ky + (double)kz * kz;
        e += 0.5 * w * m2;
        z += 0.5 * w * k2 * m2;
        if (k2 > 0.0) {
            double dr = kx * (double)a.x + ky * (double)b.x + kz * (double)c.x;
            double di = kx * (double)a.y + ky * (double)b.y + kz * (double)c.y;
            d = fmax(d, sqrt((dr * dr + di * di) / k2));
        }
    }
    se[threadIdx.x] = e;
    sz[threadIdx.x] = z;
    sd[threadIdx.x] = d;
    __syncthreads();
    for (unsigned s = blockDim.x / 2; s > 0; s >>= 1) {
        if (threadIdx.x < s) {
            se[threadIdx.x] += se[threadIdx.x + s];
            sz[threadIdx.x] += sz[threadIdx.x + s];
            sd[threadIdx.x] = fmax(sd[threadIdx.x], sd[threadIdx.x + s]);
        }
        __syncthreads();
    }
    if (threadIdx.x == 0) {
        out[3 * blockIdx.x] = se[0];
        out[3 * blockIdx.x + 1] = sz[0];
        out[3 * blockIdx.x + 2] = sd[0];
    }
}

// Block partial max of |u|^2 over physical points.
extern "C" __global__ void max_speed2_partial(const float *u0, const float *u1, const float *u2,
                                              float *out, long nr) {
    extern __shared__ float shf[];
    float m = 0.0f;
    for (long p = blockIdx.x * (long)blockDim.x + threadIdx.x; p < nr;
         p += (long)gridDim.x * blockDim.x) {
        float s = u0[p] * u0[p] + u1[p] * u1[p] + u2[p] * u2[p];
        m = fmaxf(m, s);
    }
    shf[threadIdx.x] = m;
    __syncthreads();
    for (unsigned s = blockDim.x / 2; s > 0; s >>= 1) {
        if (threadIdx.x < s) shf[threadIdx.x] = fmaxf(shf[threadIdx.x], shf[threadIdx.x + s]);
        __syncthreads();
    }
    if (threadIdx.x == 0) out[blockIdx.x] = shf[0];
}

// Brinkman penalization in physical space (see spectral-core/src/penalization.rs):
//   n -= chi/eta (u - us) + sponge/eta_s (u - uinf)
extern "C" __global__ void penalize(const float *u0, const float *u1, const float *u2,
                                    float *n0, float *n1, float *n2,
                                    const float *chi, const float *sponge,
                                    float inv_eta, float inv_eta_s,
                                    float us0, float us1, float us2,
                                    float ui0, float ui1, float ui2, long nr) {
    long p = blockIdx.x * (long)blockDim.x + threadIdx.x;
    if (p >= nr) return;
    float cs = chi[p] * inv_eta, cf = sponge[p] * inv_eta_s;
    float a = u0[p], b = u1[p], c = u2[p];
    n0[p] -= cs * (a - us0) + cf * (a - ui0);
    n1[p] -= cs * (b - us1) + cf * (b - ui1);
    n2[p] -= cs * (c - us2) + cf * (c - ui2);
}

// Block partial sums (double) of chi*(u - us): out has 3 doubles per block.
extern "C" __global__ void solid_force_partial(const float *u0, const float *u1, const float *u2,
                                               const float *chi, float us0, float us1, float us2,
                                               double *out, long nr) {
    extern __shared__ double shs[];
    double *s0 = shs, *s1 = shs + blockDim.x, *s2 = shs + 2 * blockDim.x;
    double a = 0.0, b = 0.0, c = 0.0;
    for (long p = blockIdx.x * (long)blockDim.x + threadIdx.x; p < nr;
         p += (long)gridDim.x * blockDim.x) {
        double w = chi[p];
        a += w * ((double)u0[p] - us0);
        b += w * ((double)u1[p] - us1);
        c += w * ((double)u2[p] - us2);
    }
    s0[threadIdx.x] = a; s1[threadIdx.x] = b; s2[threadIdx.x] = c;
    __syncthreads();
    for (unsigned s = blockDim.x / 2; s > 0; s >>= 1) {
        if (threadIdx.x < s) {
            s0[threadIdx.x] += s0[threadIdx.x + s];
            s1[threadIdx.x] += s1[threadIdx.x + s];
            s2[threadIdx.x] += s2[threadIdx.x + s];
        }
        __syncthreads();
    }
    if (threadIdx.x == 0) {
        out[3 * blockIdx.x] = s0[0];
        out[3 * blockIdx.x + 1] = s1[0];
        out[3 * blockIdx.x + 2] = s2[0];
    }
}

// ---------------------------------------------------------------------------
// IMEX / implicit penalization (see spectral-core/src/imex.rs for the math)

// out += sum_m c[m] * x_m  (spectral, one component). Terms with c == 0 skipped.
extern "C" __global__ void accum4(float2 *out, const float2 *x0, const float2 *x1,
                                  const float2 *x2, const float2 *x3,
                                  float c0, float c1, float c2, float c3, long ns) {
    long i = blockIdx.x * (long)blockDim.x + threadIdx.x;
    if (i >= ns) return;
    float2 acc = out[i];
    if (c0 != 0.0f) { float2 v = x0[i]; acc.x = fmaf(c0, v.x, acc.x); acc.y = fmaf(c0, v.y, acc.y); }
    if (c1 != 0.0f) { float2 v = x1[i]; acc.x = fmaf(c1, v.x, acc.x); acc.y = fmaf(c1, v.y, acc.y); }
    if (c2 != 0.0f) { float2 v = x2[i]; acc.x = fmaf(c2, v.x, acc.x); acc.y = fmaf(c2, v.y, acc.y); }
    if (c3 != 0.0f) { float2 v = x3[i]; acc.x = fmaf(c3, v.x, acc.x); acc.y = fmaf(c3, v.y, acc.y); }
    out[i] = acc;
}

// out = x / (1 + a nu |k|^2)   (out may alias x)
extern "C" __global__ void hinv(float2 *out, const float2 *x, float anu, int n, int nh) {
    long idx = blockIdx.x * (long)blockDim.x + threadIdx.x;
    long ns = (long)n * n * nh;
    if (idx >= ns) return;
    float kx, ky, kz;
    kvec(idx, n, nh, kx, ky, kz);
    float h = 1.0f + anu * (kx * kx + ky * ky + kz * kz);
    float2 v = x[idx];
    out[idx] = make_float2(v.x / h, v.y / h);
}

// out -= hinv(x):  U_hat = H^-1 (r_hat - P g_hat), applied after U_hat = H^-1 r_hat
extern "C" __global__ void sub_hinv(float2 *out, const float2 *x, float anu, int n, int nh) {
    long idx = blockIdx.x * (long)blockDim.x + threadIdx.x;
    long ns = (long)n * n * nh;
    if (idx >= ns) return;
    float kx, ky, kz;
    kvec(idx, n, nh, kx, ky, kz);
    float h = 1.0f + anu * (kx * kx + ky * ky + kz * kz);
    float2 v = x[idx], o = out[idx];
    out[idx] = make_float2(o.x - v.x / h, o.y - v.y / h);
}

// Spectral preconditioner: full <- sol * h/(1 + dbar h) + (scale*full - sol)/dbar
// (`full` is an unnormalized cuFFT forward transform, `sol` its scaled projection)
extern "C" __global__ void precond_combine(float2 *full, const float2 *sol, float anu,
                                           float dbar, float scale, int n, int nh) {
    long idx = blockIdx.x * (long)blockDim.x + threadIdx.x;
    long ns = (long)n * n * nh;
    if (idx >= ns) return;
    float kx, ky, kz;
    kvec(idx, n, nh, kx, ky, kz);
    float h = 1.0f + anu * (kx * kx + ky * ky + kz * kz);
    float ms = h / (1.0f + dbar * h), mg = 1.0f / dbar;
    float2 f = full[idx], s = sol[idx];
    f.x *= scale; f.y *= scale;
    full[idx] = make_float2(s.x * ms + (f.x - s.x) * mg, s.y * ms + (f.y - s.y) * mg);
}

// L_i = (U - r) * inv_a
extern "C" __global__ void stage_lin(float2 *l, const float2 *u, const float2 *r, float inv_a, long ns) {
    long i = blockIdx.x * (long)blockDim.x + threadIdx.x;
    if (i >= ns) return;
    float2 a = u[i], b = r[i];
    l[i] = make_float2((a.x - b.x) * inv_a, (a.y - b.y) * inv_a);
}

// b = chi > 0 ? b - us : 0   (one component)
extern "C" __global__ void mask_sub(float *b, const float *chi, float us, long nr) {
    long p = blockIdx.x * (long)blockDim.x + threadIdx.x;
    if (p >= nr) return;
    b[p] = chi[p] > 0.0f ? b[p] - us : 0.0f;
}

// ap = chi > 0 ? ap + (eta/(a chi)) p : 0
extern "C" __global__ void apply_finish(float *ap, const float *p, const float *chi,
                                        float eta_over_a, long nr) {
    long q = blockIdx.x * (long)blockDim.x + threadIdx.x;
    if (q >= nr) return;
    float x = chi[q];
    ap[q] = x > 0.0f ? ap[q] + (eta_over_a / x) * p[q] : 0.0f;
}

// z = chi > 0 ? z : 0
extern "C" __global__ void mask_zero(float *z, const float *chi, long nr) {
    long q = blockIdx.x * (long)blockDim.x + threadIdx.x;
    if (q >= nr) return;
    if (chi[q] <= 0.0f) z[q] = 0.0f;
}

// x += alpha * y
extern "C" __global__ void axpy_f(float *x, const float *y, float alpha, long nr) {
    long q = blockIdx.x * (long)blockDim.x + threadIdx.x;
    if (q >= nr) return;
    x[q] = fmaf(alpha, y[q], x[q]);
}

// p = z + beta * p
extern "C" __global__ void xpby_f(float *p, const float *z, float beta, long nr) {
    long q = blockIdx.x * (long)blockDim.x + threadIdx.x;
    if (q >= nr) return;
    p[q] = fmaf(beta, p[q], z[q]);
}

// r = b - ap
extern "C" __global__ void sub_f(float *r, const float *b, const float *ap, long nr) {
    long q = blockIdx.x * (long)blockDim.x + threadIdx.x;
    if (q >= nr) return;
    r[q] = b[q] - ap[q];
}

// Block partial sums (double) of sum_c a_c . b_c over three components.
extern "C" __global__ void dot3_partial(const float *a0, const float *a1, const float *a2,
                                        const float *b0, const float *b1, const float *b2,
                                        double *out, long nr) {
    extern __shared__ double shd[];
    double acc = 0.0;
    for (long q = blockIdx.x * (long)blockDim.x + threadIdx.x; q < nr;
         q += (long)gridDim.x * blockDim.x) {
        acc += (double)a0[q] * b0[q] + (double)a1[q] * b1[q] + (double)a2[q] * b2[q];
    }
    shd[threadIdx.x] = acc;
    __syncthreads();
    for (unsigned s = blockDim.x / 2; s > 0; s >>= 1) {
        if (threadIdx.x < s) shd[threadIdx.x] += shd[threadIdx.x + s];
        __syncthreads();
    }
    if (threadIdx.x == 0) out[blockIdx.x] = shd[0];
}

// Block partial sums (double) of each component of g: out has 3 doubles per block.
extern "C" __global__ void sum3_partial(const float *g0, const float *g1, const float *g2,
                                        double *out, long nr) {
    extern __shared__ double sh3[];
    double *s0 = sh3, *s1 = sh3 + blockDim.x, *s2 = sh3 + 2 * blockDim.x;
    double a = 0.0, b = 0.0, c = 0.0;
    for (long q = blockIdx.x * (long)blockDim.x + threadIdx.x; q < nr;
         q += (long)gridDim.x * blockDim.x) {
        a += g0[q]; b += g1[q]; c += g2[q];
    }
    s0[threadIdx.x] = a; s1[threadIdx.x] = b; s2[threadIdx.x] = c;
    __syncthreads();
    for (unsigned s = blockDim.x / 2; s > 0; s >>= 1) {
        if (threadIdx.x < s) {
            s0[threadIdx.x] += s0[threadIdx.x + s];
            s1[threadIdx.x] += s1[threadIdx.x + s];
            s2[threadIdx.x] += s2[threadIdx.x + s];
        }
        __syncthreads();
    }
    if (threadIdx.x == 0) {
        out[3 * blockIdx.x] = s0[0]; out[3 * blockIdx.x + 1] = s1[0]; out[3 * blockIdx.x + 2] = s2[0];
    }
}

// v /= (1 + a nu |k|^2) in place
extern "C" __global__ void hinv_inplace(float2 *v, float anu, int n, int nh) {
    long idx = blockIdx.x * (long)blockDim.x + threadIdx.x;
    long ns = (long)n * n * nh;
    if (idx >= ns) return;
    float kx, ky, kz;
    kvec(idx, n, nh, kx, ky, kz);
    float h = 1.0f + anu * (kx * kx + ky * ky + kz * kz);
    float2 x = v[idx];
    v[idx] = make_float2(x.x / h, x.y / h);
}

// ---------------------------------------------------------------------------
// Non-uniform FFT by Gaussian gridding (see spectral-core/src/nufft.rs)

// out_m[pm[q]] = in_n[pn[q]] * scale_q  (scatter between spectra; out_m pre-zeroed)
extern "C" __global__ void spec_scatter(float2 *out_m, const float2 *in_n, const int *pn,
                                        const int *pm, const float *scale_q, int npairs) {
    int q = blockIdx.x * blockDim.x + threadIdx.x;
    if (q >= npairs) return;
    float2 v = in_n[pn[q]];
    float s = scale_q[q];
    out_m[pm[q]] = make_float2(v.x * s, v.y * s);
}

// out_n[pn[q]] = in_m[pm[q]] * scale_q * extra  (gather; out_n pre-zeroed)
extern "C" __global__ void spec_gather(float2 *out_n, const float2 *in_m, const int *pn,
                                       const int *pm, const float *scale_q, float extra,
                                       int npairs) {
    int q = blockIdx.x * blockDim.x + threadIdx.x;
    if (q >= npairs) return;
    float2 v = in_m[pm[q]];
    float s = scale_q[q] * extra;
    out_n[pn[q]] = make_float2(v.x * s, v.y * s);
}

__device__ __forceinline__ void stencil1(float x, float h, int hw, float tau, int &base,
                                         float *w) {
    base = (int)floorf(x / h) - hw + 1;
    for (int q = 0; q < 2 * hw; ++q) {
        float d = x - (float)(base + q) * h;
        w[q] = expf(-d * d / (4.0f * tau));
    }
}

// One thread per point: value = h^3 sum_j grid[j] g(X - x_j)
extern "C" __global__ void nufft_interp(float *values, const float *grid, const float *pts,
                                        int npts, int m, int hw, float tau, float h) {
    int p = blockIdx.x * blockDim.x + threadIdx.x;
    if (p >= npts) return;
    float wx[32], wy[32], wz[32];
    int bx, by, bz;
    stencil1(pts[3 * p], h, hw, tau, bx, wx);
    stencil1(pts[3 * p + 1], h, hw, tau, by, wy);
    stencil1(pts[3 * p + 2], h, hw, tau, bz, wz);
    float acc = 0.0f;
    for (int a = 0; a < 2 * hw; ++a) {
        int i = ((bx + a) % m + m) % m;
        for (int b = 0; b < 2 * hw; ++b) {
            int j = ((by + b) % m + m) % m;
            long row = ((long)i * m + j) * m;
            float wab = wx[a] * wy[b];
            for (int c = 0; c < 2 * hw; ++c) {
                int l = ((bz + c) % m + m) % m;
                acc = fmaf(wab * wz[c], grid[row + l], acc);
            }
        }
    }
    values[p] = acc * h * h * h;
}

// One thread per point: grid[j] += F g(x_j - X)  (atomic)
extern "C" __global__ void nufft_spread(float *grid, const float *pts, const float *forces,
                                        int npts, int m, int hw, float tau, float h) {
    int p = blockIdx.x * blockDim.x + threadIdx.x;
    if (p >= npts) return;
    float wx[32], wy[32], wz[32];
    int bx, by, bz;
    stencil1(pts[3 * p], h, hw, tau, bx, wx);
    stencil1(pts[3 * p + 1], h, hw, tau, by, wy);
    stencil1(pts[3 * p + 2], h, hw, tau, bz, wz);
    float f = forces[p];
    for (int a = 0; a < 2 * hw; ++a) {
        int i = ((bx + a) % m + m) % m;
        for (int b = 0; b < 2 * hw; ++b) {
            int j = ((by + b) % m + m) % m;
            long row = ((long)i * m + j) * m;
            float wab = f * wx[a] * wy[b];
            for (int c = 0; c < 2 * hw; ++c) {
                int l = ((bz + c) % m + m) % m;
                atomicAdd(&grid[row + l], wab * wz[c]);
            }
        }
    }
}

// out[c] += coef/(1 + anu |k|^2) * x[c]   (one component)
extern "C" __global__ void add_hinv_scaled(float2 *out, const float2 *x, float coef, float anu,
                                           int n, int nh) {
    long idx = blockIdx.x * (long)blockDim.x + threadIdx.x;
    long ns = (long)n * n * nh;
    if (idx >= ns) return;
    float kx, ky, kz;
    kvec(idx, n, nh, kx, ky, kz);
    float s = coef / (1.0f + anu * (kx * kx + ky * ky + kz * kz));
    float2 v = x[idx], o = out[idx];
    out[idx] = make_float2(fmaf(s, v.x, o.x), fmaf(s, v.y, o.y));
}

// ---------------------------------------------------------------------------
// Smooth forcing extension (see spectral-core/src/surface.rs)

// phys[grid_idx[r]] = sum_{q in row r} val[q] * amp[col[q]]   (phys pre-zeroed)
extern "C" __global__ void ext_build(float *phys, const float *amp, const int *grid_idx,
                                     const int *row_ptr, const int *col, const float *val,
                                     int nrows) {
    int r = blockIdx.x * blockDim.x + threadIdx.x;
    if (r >= nrows) return;
    float acc = 0.0f;
    for (int q = row_ptr[r]; q < row_ptr[r + 1]; ++q) acc = fmaf(val[q], amp[col[q]], acc);
    phys[grid_idx[r]] = acc;
}

// res[col[q]] += val[q] * phys[grid_idx[r]]   (res pre-zeroed; E^T)
extern "C" __global__ void ext_adjoint(float *res, const float *phys, const int *grid_idx,
                                       const int *row_ptr, const int *col, const float *val,
                                       int nrows) {
    int r = blockIdx.x * blockDim.x + threadIdx.x;
    if (r >= nrows) return;
    float g = phys[grid_idx[r]];
    for (int q = row_ptr[r]; q < row_ptr[r + 1]; ++q) atomicAdd(&res[col[q]], val[q] * g);
}

// v *= s  (spectral)
extern "C" __global__ void spec_scale(float2 *v, float s, long ns) {
    long i = blockIdx.x * (long)blockDim.x + threadIdx.x;
    if (i >= ns) return;
    v[i].x *= s;
    v[i].y *= s;
}

// Penalization toward a target field: n -= chi/eta (u - t) + sponge/eta_s (u - uinf)
extern "C" __global__ void penalize_field(const float *u0, const float *u1, const float *u2,
                                          float *n0, float *n1, float *n2,
                                          const float *chi, const float *sponge,
                                          const float *t0, const float *t1, const float *t2,
                                          float inv_eta, float inv_eta_s,
                                          float ui0, float ui1, float ui2, long nr) {
    long p = blockIdx.x * (long)blockDim.x + threadIdx.x;
    if (p >= nr) return;
    float cs = chi[p] * inv_eta, cf = sponge[p] * inv_eta_s;
    float a = u0[p], b = u1[p], c = u2[p];
    n0[p] -= cs * (a - t0[p]) + cf * (a - ui0);
    n1[p] -= cs * (b - t1[p]) + cf * (b - ui1);
    n2[p] -= cs * (c - t2[p]) + cf * (c - ui2);
}

// ---------------------------------------------------------------------------
// A posteriori certificate (see spectral-core/src/certify.rs)

// Block partial sums (double) of sum over padded modes NOT retained
// (|kx|,|ky|,|kz| <= kmax) of weight * |P v|^2 / |k|^2.
extern "C" __global__ void trunc_resid_partial(const float2 *v0, const float2 *v1,
                                               const float2 *v2, double *out, int m, int mh,
                                               int kmax) {
    extern __shared__ double shr[];
    long ns = (long)m * m * mh;
    double acc = 0.0;
    for (long idx = blockIdx.x * (long)blockDim.x + threadIdx.x; idx < ns;
         idx += (long)gridDim.x * blockDim.x) {
        float kx, ky, kz;
        kvec(idx, m, mh, kx, ky, kz);
        if (fabsf(kx) <= kmax && fabsf(ky) <= kmax && kz <= kmax) continue;
        double k2 = (double)kx * kx + (double)ky * ky + (double)kz * kz;
        if (k2 == 0.0) continue;
        float2 a = v0[idx], b = v1[idx], c = v2[idx];
        double dr = (kx * (double)a.x + ky * (double)b.x + kz * (double)c.x) / k2;
        double di = (kx * (double)a.y + ky * (double)b.y + kz * (double)c.y) / k2;
        double px = a.x - dr * kx, py = b.x - dr * ky, pz = c.x - dr * kz;
        double qx = a.y - di * kx, qy = b.y - di * ky, qz = c.y - di * kz;
        int l = (int)(idx % mh);
        double w = (l == 0 || (m % 2 == 0 && l == m / 2)) ? 1.0 : 2.0;
        acc += w * (px * px + py * py + pz * pz + qx * qx + qy * qy + qz * qz) / k2;
    }
    shr[threadIdx.x] = acc;
    __syncthreads();
    for (unsigned s = blockDim.x / 2; s > 0; s >>= 1) {
        if (threadIdx.x < s) shr[threadIdx.x] += shr[threadIdx.x + s];
        __syncthreads();
    }
    if (threadIdx.x == 0) out[blockIdx.x] = shr[0];
}

// s_ij = (i/2)(k_i u_j + k_j u_i) on a spectral grid (n, nh)
extern "C" __global__ void strain_spec(const float2 *ui, const float2 *uj, float2 *s, int i,
                                       int j, int n, int nh) {
    long idx = blockIdx.x * (long)blockDim.x + threadIdx.x;
    long ns = (long)n * n * nh;
    if (idx >= ns) return;
    float k[3];
    kvec(idx, n, nh, k[0], k[1], k[2]);
    float2 a = ui[idx], b = uj[idx];
    // t = k_i u_j + k_j u_i ; s = (i/2) t = (-t.y/2, t.x/2)
    float tx = k[i] * b.x + k[j] * a.x, ty = k[i] * b.y + k[j] * a.y;
    s[idx] = make_float2(-0.5f * ty, 0.5f * tx);
}

__device__ double sym3_min(double a, double b, double c, double d, double e, double f) {
    double p1 = d * d + e * e + f * f;
    double q = (a + b + c) / 3.0;
    double p2 = (a - q) * (a - q) + (b - q) * (b - q) + (c - q) * (c - q) + 2.0 * p1;
    if (p2 <= 1e-300) return q;
    double p = sqrt(p2 / 6.0);
    double ba = (a - q) / p, bb = (b - q) / p, bc = (c - q) / p, bd = d / p, be = e / p, bf = f / p;
    double det = ba * (bb * bc - bf * bf) - bd * (bd * bc - bf * be) + be * (bd * bf - bb * be);
    double r = fmin(1.0, fmax(-1.0, det / 2.0));
    double phi = acos(r) / 3.0;
    return q + 2.0 * p * cos(phi + 2.0943951023931957);
}

// Block partial max of -lambda_min(S) over physical points (S as 6 fields:
// xx, yy, zz, xy, xz, yz).
extern "C" __global__ void min_eig_partial(const float *sxx, const float *syy, const float *szz,
                                           const float *sxy, const float *sxz, const float *syz,
                                           double *out, long nr) {
    extern __shared__ double shm[];
    double m = 0.0;
    for (long p = blockIdx.x * (long)blockDim.x + threadIdx.x; p < nr;
         p += (long)gridDim.x * blockDim.x) {
        double v = -sym3_min(sxx[p], syy[p], szz[p], sxy[p], sxz[p], syz[p]);
        m = fmax(m, v);
    }
    shm[threadIdx.x] = m;
    __syncthreads();
    for (unsigned s = blockDim.x / 2; s > 0; s >>= 1) {
        if (threadIdx.x < s) shm[threadIdx.x] = fmax(shm[threadIdx.x], shm[threadIdx.x + s]);
        __syncthreads();
    }
    if (threadIdx.x == 0) out[blockIdx.x] = shm[0];
}

// f -= nu |k|^2 u   (turns the nonlinear tendency into the full tendency)
extern "C" __global__ void sub_viscous(float2 *f, const float2 *u, float nu, int n, int nh) {
    long idx = blockIdx.x * (long)blockDim.x + threadIdx.x;
    long ns = (long)n * n * nh;
    if (idx >= ns) return;
    float kx, ky, kz;
    kvec(idx, n, nh, kx, ky, kz);
    float s = nu * (kx * kx + ky * ky + kz * kz);
    float2 a = f[idx], b = u[idx];
    f[idx] = make_float2(a.x - s * b.x, a.y - s * b.y);
}

// Hermite midpoint: vm = (u0+u1)/2 + h/8 (f0-f1); vdot = 3/(2h)(u1-u0) - (f0+f1)/4
extern "C" __global__ void hermite_mid(float2 *vm, float2 *vdot, const float2 *u0,
                                       const float2 *u1, const float2 *f0, const float2 *f1,
                                       float h, long ns) {
    long i = blockIdx.x * (long)blockDim.x + threadIdx.x;
    if (i >= ns) return;
    float2 a = u0[i], b = u1[i], c = f0[i], d = f1[i];
    vm[i] = make_float2(0.5f * (a.x + b.x) + 0.125f * h * (c.x - d.x),
                        0.5f * (a.y + b.y) + 0.125f * h * (c.y - d.y));
    float s = 1.5f / h;
    vdot[i] = make_float2(s * (b.x - a.x) - 0.25f * (c.x + d.x), s * (b.y - a.y) - 0.25f * (c.y + d.y));
}

// r -= x  (spectral)
extern "C" __global__ void spec_sub(float2 *r, const float2 *x, long ns) {
    long i = blockIdx.x * (long)blockDim.x + threadIdx.x;
    if (i >= ns) return;
    r[i].x -= x[i].x;
    r[i].y -= x[i].y;
}

// Block partial sums (double) of weight * |v|^2 / |k|^2 over k != 0 (3 comps)
extern "C" __global__ void hm1_partial(const float2 *v0, const float2 *v1, const float2 *v2,
                                       double *out, int n, int nh) {
    extern __shared__ double shh[];
    long ns = (long)n * n * nh;
    double acc = 0.0;
    for (long idx = blockIdx.x * (long)blockDim.x + threadIdx.x; idx < ns;
         idx += (long)gridDim.x * blockDim.x) {
        float kx, ky, kz;
        kvec(idx, n, nh, kx, ky, kz);
        double k2 = (double)kx * kx + (double)ky * ky + (double)kz * kz;
        if (k2 == 0.0) continue;
        int l = (int)(idx % nh);
        double w = (l == 0 || (n % 2 == 0 && l == n / 2)) ? 1.0 : 2.0;
        float2 a = v0[idx], b = v1[idx], c = v2[idx];
        acc += w * ((double)a.x * a.x + (double)a.y * a.y + (double)b.x * b.x + (double)b.y * b.y +
                    (double)c.x * c.x + (double)c.y * c.y) / k2;
    }
    shh[threadIdx.x] = acc;
    __syncthreads();
    for (unsigned s = blockDim.x / 2; s > 0; s >>= 1) {
        if (threadIdx.x < s) shh[threadIdx.x] += shh[threadIdx.x + s];
        __syncthreads();
    }
    if (threadIdx.x == 0) out[blockIdx.x] = shh[0];
}
