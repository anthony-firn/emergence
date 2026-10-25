// Compressible hybrid central / WCNS scheme on the GPU (NVRTC).
// Mirrors fluid/crates/compressible/src/{scheme,system,multiphase}.rs.
// Templated on the real type R (float/double) and the system (Euler: 5
// variables, FiveEq: 7). Structure-of-arrays storage: u[q * ncell + cell].

#define MAXN 7

template <typename R> __device__ __forceinline__ R sq(R x) { return x * x; }
__device__ __forceinline__ float rlog(float x) { return logf(x); }
__device__ __forceinline__ double rlog(double x) { return log(x); }
__device__ __forceinline__ float rsqrt_(float x) { return sqrtf(x); }
__device__ __forceinline__ double rsqrt_(double x) { return sqrt(x); }
template <typename R> __device__ __forceinline__ R rabs(R x) { return x < R(0) ? -x : x; }
template <typename R> __device__ __forceinline__ R rmax(R a, R b) { return a > b ? a : b; }
template <typename R> __device__ __forceinline__ R rmin(R a, R b) { return a < b ? a : b; }
template <typename R> __device__ __forceinline__ bool fin(R x) { return isfinite(x); }

template <typename R> struct Nasg {
    R gamma, pinf, b, q;
    // ρe = (p + γp∞)(1 − bρ)/(γ − 1) + qρ
    __device__ R rho_e(R rho, R p) const { return (p + gamma * pinf) * (R(1) - b * rho) / (gamma - R(1)) + q * rho; }
    __device__ R pressure(R rho, R e) const { R v = R(1) / rho; return (gamma - R(1)) * (e - q) / (v - b) - gamma * pinf; }
    __device__ R c2(R rho, R p) const { R v = R(1) / rho; return gamma * (p + pinf) * v * v / (v - b); }
};

// ---------------------------------------------------------------- Euler
template <typename R> struct Euler {
    static const int N = 5;
    static const int VI = 1;
    Nasg<R> e;
    __device__ void to_prim(const R* u, R* w) const {
        R rho = u[0];
        R vx = u[1] / rho, vy = u[2] / rho, vz = u[3] / rho;
        R ein = u[4] / rho - R(0.5) * (vx * vx + vy * vy + vz * vz);
        w[0] = rho; w[1] = vx; w[2] = vy; w[3] = vz; w[4] = e.pressure(rho, ein);
    }
    __device__ void to_cons(const R* w, R* u) const {
        R rho = w[0];
        u[0] = rho; u[1] = rho * w[1]; u[2] = rho * w[2]; u[3] = rho * w[3];
        u[4] = e.rho_e(rho, w[4]) + R(0.5) * rho * (w[1] * w[1] + w[2] * w[2] + w[3] * w[3]);
    }
    __device__ R c(const R* w) const { return rsqrt_(rmax(e.c2(w[0], w[4]), R(0))); }
    __device__ R density(const R* w) const { return w[0]; }
    __device__ bool admissible(const R* w) const {
        return fin(w[0]) && fin(w[4]) && w[0] > R(0) && R(1) / w[0] > e.b && w[4] + e.pinf > R(0);
    }
    __device__ int sensor(const R* w, R* s) const {
        s[0] = rlog(rmax(w[0], R(1e-300)));
        s[1] = rlog(rmax(w[4] + e.pinf, R(1e-300)));
        return 2;
    }
    __device__ void flux(const R* w, int d, R* f) const {
        R q[5]; to_cons(w, q);
        R un = w[1 + d];
        for (int k = 0; k < 4; k++) f[k] = q[k] * un;
        f[1 + d] += w[4];
        f[4] = (q[4] + w[4]) * un;
    }
    __device__ void split(const R* a, const R* b, int d, R* f) const {
        int un = 1 + d;
        R rho = R(0.5) * (a[0] + b[0]), u = R(0.5) * (a[un] + b[un]), cc = rho * u;
        f[0] = cc;
        for (int k = 1; k < 4; k++) f[k] = cc * R(0.5) * (a[k] + b[k]);
        f[un] += R(0.5) * (a[4] + b[4]);
        R kin = cc * R(0.5) * (a[1] * b[1] + a[2] * b[2] + a[3] * b[3]);
        R rhoe = R(0.5) * (e.rho_e(a[0], a[4]) + e.rho_e(b[0], b[4]));
        R pw = R(0.5) * (a[un] * b[4] + b[un] * a[4]);
        f[4] = kin + u * rhoe + pw;
    }
    // primitive characteristic eigenvectors (rows of L, columns of R)
    __device__ void eigen(const R* w, int d, R L[][MAXN], R Rm[][MAXN]) const {
        for (int i = 0; i < 5; i++) for (int j = 0; j < 5; j++) { L[i][j] = 0; Rm[i][j] = 0; }
        R rho = w[0], cs = c(w);
        int un = 1 + d, t1 = 1 + (d + 1) % 3, t2 = 1 + (d + 2) % 3;
        L[0][un] = -rho / (R(2) * cs); L[0][4] = R(1) / (R(2) * cs * cs);
        L[1][0] = 1; L[1][4] = -R(1) / (cs * cs);
        L[2][t1] = 1; L[3][t2] = 1;
        L[4][un] = rho / (R(2) * cs); L[4][4] = R(1) / (R(2) * cs * cs);
        Rm[0][0] = 1; Rm[un][0] = -cs / rho; Rm[4][0] = cs * cs;
        Rm[0][1] = 1; Rm[t1][2] = 1; Rm[t2][3] = 1;
        Rm[0][4] = 1; Rm[un][4] = cs / rho; Rm[4][4] = cs * cs;
    }
    __device__ void eigenvalues(const R* w, int d, R* l) const {
        R u = w[1 + d], cs = c(w);
        l[0] = u - cs; l[1] = u; l[2] = u; l[3] = u; l[4] = u + cs;
    }
    // HLLC; returns face velocity
    __device__ R riemann(const R* wl, const R* wr, int d, R* f) const {
        int un = 1 + d;
        R cl = c(wl), cr = c(wr), ul = wl[un], ur = wr[un];
        R sl = rmin(ul - cl, ur - cr), sr = rmax(ul + cl, ur + cr);
        R ql[5], qr[5], fl[5], fr[5];
        to_cons(wl, ql); to_cons(wr, qr); flux(wl, d, fl); flux(wr, d, fr);
        if (sl >= 0) { for (int k = 0; k < 5; k++) f[k] = fl[k]; return ul; }
        if (sr <= 0) { for (int k = 0; k < 5; k++) f[k] = fr[k]; return ur; }
        R rl = wl[0], rr = wr[0], pl = wl[4], pr = wr[4];
        R ss = (pr - pl + rl * ul * (sl - ul) - rr * ur * (sr - ur)) / (rl * (sl - ul) - rr * (sr - ur));
        const R* w = ss >= 0 ? wl : wr; const R* q = ss >= 0 ? ql : qr; const R* fx = ss >= 0 ? fl : fr;
        R s = ss >= 0 ? sl : sr;
        R rho = w[0], u = w[un], fac = rho * (s - u) / (s - ss), qs[5];
        qs[0] = fac;
        for (int k = 1; k < 4; k++) qs[k] = fac * w[k];
        qs[un] = fac * ss;
        qs[4] = fac * (q[4] / rho + (ss - u) * (ss + w[4] / (rho * (s - u))));
        for (int k = 0; k < 5; k++) f[k] = fx[k] + s * (qs[k] - q[k]);
        return ss;
    }
    __device__ bool has_nc() const { return false; }
    __device__ int nc_index() const { return -1; }
    __device__ R nc_coef(const R*) const { return 0; }
    __device__ R dilatation_rate(const R*, R) const { return 0; }
};

// ---------------------------------------------------------------- FiveEq
// conserved: a1, a2, ρu, ρv, ρw, E, α1 ; primitive: a1, a2, u, v, w, p, α1
template <typename R> struct FiveEq {
    static const int N = 7;
    static const int VI = 2;
    Nasg<R> e1, e2;
    int kapila;
    __device__ void db(const R* w, R& D, R& B) const {
        R al[2] = {w[6], R(1) - w[6]};
        const Nasg<R>* es[2] = {&e1, &e2};
        D = 0; B = 0;
        for (int k = 0; k < 2; k++) {
            R vol = al[k] - es[k]->b * w[k];
            D += vol / (es[k]->gamma - R(1));
            B += es[k]->gamma * es[k]->pinf * vol / (es[k]->gamma - R(1)) + es[k]->q * w[k];
        }
    }
    __device__ R rho_e(const R* w) const { R D, B; db(w, D, B); return w[5] * D + B; }
    __device__ void to_prim(const R* u, R* w) const {
        R rho = u[0] + u[1];
        w[0] = u[0]; w[1] = u[1]; w[2] = u[2] / rho; w[3] = u[3] / rho; w[4] = u[4] / rho; w[6] = u[6];
        R rhoe = u[5] - R(0.5) * rho * (w[2] * w[2] + w[3] * w[3] + w[4] * w[4]);
        R D, B; db(w, D, B);
        w[5] = (rhoe - B) / D;
    }
    __device__ void to_cons(const R* w, R* u) const {
        R rho = w[0] + w[1];
        u[0] = w[0]; u[1] = w[1]; u[2] = rho * w[2]; u[3] = rho * w[3]; u[4] = rho * w[4];
        u[5] = rho_e(w) + R(0.5) * rho * (w[2] * w[2] + w[3] * w[3] + w[4] * w[4]);
        u[6] = w[6];
    }
    __device__ R rc2_phase(const R* w, int k) const {
        R al = k == 0 ? w[6] : R(1) - w[6];
        const Nasg<R>& e = k == 0 ? e1 : e2;
        R rk = (al > R(1e-12) && w[k] > 0) ? w[k] / al : R(0);
        return e.gamma * (w[5] + e.pinf) / rmax(R(1) - e.b * rk, R(1e-12));
    }
    __device__ R kk(const R* w) const {
        if (!kapila) return 0;
        R a = rmin(rmax(w[6], R(0)), R(1)), al[2] = {a, R(1) - a};
        R z1 = rc2_phase(w, 0), z2 = rc2_phase(w, 1);
        R den = al[1] * z1 + al[0] * z2;
        if (den <= 0) return 0;
        return al[0] * al[1] * (z2 - z1) / den;
    }
    __device__ R drhoe_dalpha(R p) const {
        return (p + e1.gamma * e1.pinf) / (e1.gamma - R(1)) - (p + e2.gamma * e2.pinf) / (e2.gamma - R(1));
    }
    __device__ R rho_c2(const R* w) const {
        R D, B; db(w, D, B);
        R p = w[5], num = rho_e(w) + p;
        const Nasg<R>* es[2] = {&e1, &e2};
        for (int k = 0; k < 2; k++) num += w[k] * (es[k]->b * (p + es[k]->gamma * es[k]->pinf) / (es[k]->gamma - R(1)) - es[k]->q);
        num += kk(w) * drhoe_dalpha(p);
        return num / D;
    }
    __device__ R density(const R* w) const { return w[0] + w[1]; }
    __device__ R c(const R* w) const { return rsqrt_(rmax(rho_c2(w), R(0)) / density(w)); }
    __device__ R pinf_mix(const R* w) const { R a = rmin(rmax(w[6], R(0)), R(1)); return a * e1.pinf + (R(1) - a) * e2.pinf; }
    __device__ bool admissible(const R* w) const {
        for (int k = 0; k < 7; k++) if (!fin(w[k])) return false;
        return density(w) > 0 && w[6] >= R(-0.05) && w[6] <= R(1.05) && rho_c2(w) > 0;
    }
    __device__ int sensor(const R* w, R* s) const {
        s[0] = rlog(rmax(density(w), R(1e-300)));
        s[1] = rlog(rmax(w[5] + pinf_mix(w), R(1e-300)));
        s[2] = w[6];
        return 3;
    }
    __device__ void flux(const R* w, int d, R* f) const {
        R q[7]; to_cons(w, q);
        R un = w[2 + d];
        for (int k = 0; k < 5; k++) f[k] = q[k] * un;
        f[2 + d] += w[5];
        f[5] = (q[5] + w[5]) * un;
        f[6] = w[6] * un;
    }
    __device__ void split(const R* a, const R* b, int d, R* f) const {
        int un = 2 + d;
        R u = R(0.5) * (a[un] + b[un]), rho = R(0.5) * (density(a) + density(b)), cc = rho * u;
        f[0] = R(0.5) * (a[0] + b[0]) * u; f[1] = R(0.5) * (a[1] + b[1]) * u;
        for (int k = 2; k < 5; k++) f[k] = cc * R(0.5) * (a[k] + b[k]);
        f[un] += R(0.5) * (a[5] + b[5]);
        R kin = cc * R(0.5) * (a[2] * b[2] + a[3] * b[3] + a[4] * b[4]);
        R rhoe = R(0.5) * (rho_e(a) + rho_e(b));
        R pw = R(0.5) * (a[un] * b[5] + b[un] * a[5]);
        f[5] = kin + u * rhoe + pw;
        f[6] = R(0.5) * (a[6] + b[6]) * u;
    }
    __device__ void eigen(const R* w, int d, R L[][MAXN], R Rm[][MAXN]) const {
        for (int i = 0; i < 7; i++) for (int j = 0; j < 7; j++) { L[i][j] = 0; Rm[i][j] = 0; }
        R rho = density(w), cs = c(w), K = kk(w), c2 = cs * cs;
        int un = 2 + d, t1 = 2 + (d + 1) % 3, t2 = 2 + (d + 2) % 3;
        L[0][0] = 1; L[0][5] = -w[0] / rho / c2;
        L[1][1] = 1; L[1][5] = -w[1] / rho / c2;
        L[2][un] = -rho / (R(2) * cs); L[2][5] = R(1) / (R(2) * c2);
        L[3][t1] = 1; L[4][t2] = 1;
        L[5][un] = rho / (R(2) * cs); L[5][5] = R(1) / (R(2) * c2);
        L[6][6] = 1; L[6][5] = K / rho / c2;
        Rm[0][0] = 1; Rm[1][1] = 1; Rm[t1][3] = 1; Rm[t2][4] = 1; Rm[6][6] = 1;
        int cols[2] = {2, 5}; R sg[2] = {-1, 1};
        for (int m = 0; m < 2; m++) {
            int col = cols[m];
            Rm[0][col] = w[0] / rho; Rm[1][col] = w[1] / rho; Rm[un][col] = sg[m] * cs / rho;
            Rm[5][col] = c2; Rm[6][col] = -K / rho;
        }
    }
    __device__ void eigenvalues(const R* w, int d, R* l) const {
        R u = w[2 + d], cs = c(w);
        l[0] = u; l[1] = u; l[2] = u - cs; l[3] = u; l[4] = u; l[5] = u + cs; l[6] = u;
    }
    __device__ R riemann(const R* wl, const R* wr, int d, R* f) const {
        int un = 2 + d;
        R cl = c(wl), cr = c(wr), ul = wl[un], ur = wr[un];
        R sl = rmin(ul - cl, ur - cr), sr = rmax(ul + cl, ur + cr);
        R ql[7], qr[7], fl[7], fr[7];
        to_cons(wl, ql); to_cons(wr, qr); flux(wl, d, fl); flux(wr, d, fr);
        if (sl >= 0) { for (int k = 0; k < 7; k++) f[k] = fl[k]; return ul; }
        if (sr <= 0) { for (int k = 0; k < 7; k++) f[k] = fr[k]; return ur; }
        R rl = density(wl), rr = density(wr), pl = wl[5], pr = wr[5];
        R ss = (pr - pl + rl * ul * (sl - ul) - rr * ur * (sr - ur)) / (rl * (sl - ul) - rr * (sr - ur));
        const R* w = ss >= 0 ? wl : wr; const R* q = ss >= 0 ? ql : qr; const R* fx = ss >= 0 ? fl : fr;
        R s = ss >= 0 ? sl : sr;
        R rho = density(w), u = w[un], fac = (s - u) / (s - ss), qs[7];
        qs[0] = w[0] * fac; qs[1] = w[1] * fac;
        for (int k = 2; k < 5; k++) qs[k] = rho * fac * w[k];
        qs[un] = rho * fac * ss;
        qs[5] = rho * fac * (q[5] / rho + (ss - u) * (ss + w[5] / (rho * (s - u))));
        qs[6] = w[6] * fac;
        for (int k = 0; k < 7; k++) f[k] = fx[k] + s * (qs[k] - q[k]);
        return ss;
    }
    __device__ bool has_nc() const { return true; }
    __device__ int nc_index() const { return 6; }
    __device__ R nc_coef(const R* w) const { return w[6] + kk(w); }
    __device__ R dilatation_rate(const R* w, R divu) const {
        R a = w[6], pm = R(1e300);
        if (a > R(1e-6)) pm = rmin(pm, e1.pinf);
        if (R(1) - a > R(1e-6)) pm = rmin(pm, e2.pinf);
        R head = rmax(rmax(w[5] + pm, R(1e-6) * rmax(e1.pinf, e2.pinf)), R(1e-300));
        return rmax(rho_c2(w), R(0)) * rabs(divu) / head;
    }
};

// ---------------------------------------------------------------- scheme
template <typename R> __device__ R weno_z(const R* q, const R* dd, R v0, R v1, R v2, R v3, R v4) {
    R b0 = R(13) / R(12) * sq(v0 - R(2) * v1 + v2) + R(0.25) * sq(v0 - R(4) * v1 + R(3) * v2);
    R b1 = R(13) / R(12) * sq(v1 - R(2) * v2 + v3) + R(0.25) * sq(v1 - v3);
    R b2 = R(13) / R(12) * sq(v2 - R(2) * v3 + v4) + R(0.25) * sq(R(3) * v2 - R(4) * v3 + v4);
    R tau = rabs(b0 - b2), eps = sizeof(R) == 4 ? R(1e-30) : R(1e-40);
    R a0 = dd[0] * (R(1) + sq(tau / (b0 + eps)));
    R a1 = dd[1] * (R(1) + sq(tau / (b1 + eps)));
    R a2 = dd[2] * (R(1) + sq(tau / (b2 + eps)));
    return (a0 * q[0] + a1 * q[1] + a2 * q[2]) / (a0 + a1 + a2);
}
template <typename R> __device__ R weno_interp(R v0, R v1, R v2, R v3, R v4) {
    R q[3] = {(R(3) * v0 - R(10) * v1 + R(15) * v2) / R(8), (-v1 + R(6) * v2 + R(3) * v3) / R(8), (R(3) * v2 + R(6) * v3 - v4) / R(8)};
    R dd[3] = {R(1) / R(16), R(10) / R(16), R(5) / R(16)};
    return weno_z(q, dd, v0, v1, v2, v3, v4);
}

// midpoint Riemann flux between local stencil cells j and j+1 (w[j-2..j+3])
template <class S, typename R> __device__ R midpoint_flux(const S& sys, R W[][MAXN], int j, int d, R* f, int* fb) {
    const int N = S::N;
    R avg[MAXN], L[MAXN][MAXN], Rm[MAXN][MAXN];
    for (int k = 0; k < N; k++) avg[k] = R(0.5) * (W[j][k] + W[j + 1][k]);
    sys.eigen(avg, d, L, Rm);
    R ch[6][MAXN];
    for (int s = 0; s < 6; s++)
        for (int a = 0; a < N; a++) {
            R acc = 0;
            for (int k = 0; k < N; k++) acc += L[a][k] * W[j - 2 + s][k];
            ch[s][a] = acc;
        }
    R cl[MAXN], cr[MAXN], wl[MAXN], wr[MAXN];
    for (int a = 0; a < N; a++) {
        cl[a] = weno_interp(ch[0][a], ch[1][a], ch[2][a], ch[3][a], ch[4][a]);
        cr[a] = weno_interp(ch[5][a], ch[4][a], ch[3][a], ch[2][a], ch[1][a]);
    }
    for (int k = 0; k < N; k++) {
        R sl = 0, sr = 0;
        for (int a = 0; a < N; a++) { sl += Rm[k][a] * cl[a]; sr += Rm[k][a] * cr[a]; }
        wl[k] = sl; wr[k] = sr;
    }
    if (!sys.admissible(wl) || !sys.admissible(wr)) {
        *fb += 1;
        for (int k = 0; k < N; k++) { wl[k] = W[j][k]; wr[k] = W[j + 1][k]; }
    }
    return sys.riemann(wl, wr, d, f);
}

// Geometry of one block (the uniform solver is one block). Launches run
// over `nb` blocks at once: block bb of the launch is pool block bl[bb]
// (field arrays u, w: u[q * ntot + bl[bb] * ncell + cell]) and right-hand
// side slot s0 + bb (rhs[q * rtot + slot * nint + interior cell]). Face and
// midpoint scratch (TH, MF, MU, F, UF) is indexed by launch-local block.
struct Geo {
    int m[3];    // allocated extents
    int n[3];    // interior cells
    int g[3];    // ghost offset (0 on inactive axes)
    long stride[3];
    long ncell;  // cells per block (with ghosts)
    long ntot;   // variable stride of u / w
    long rtot;   // variable stride of rhs
    int nb;      // blocks in this launch
    int s0;      // first rhs slot
};

// Face / midpoint indexing along direction d. Faces of a pencil: f = 0..n
// (face between interior cells f−1 and f, i.e. stencil cells i = g−1+f and
// i+1). Midpoints m = 0..n+4 cover faces f−2..f+2: midpoint m sits at face
// f = m − 2.
__device__ __forceinline__ void pencil_axes(int d, int& e0, int& e1) {
    e0 = (d + 1) % 3; e1 = (d + 2) % 3;
    if (e0 > e1) { int t = e0; e0 = e1; e1 = t; }
}

// K1: blending weight θ per face (sensor). TH[fi].
template <class S, typename R>
__device__ void theta_impl(const S& sys, const R* __restrict__ w, Geo geo, const int* __restrict__ bl, int d, int mode, R tol,
                           R* __restrict__ TH) {
    int e0, e1; pencil_axes(d, e0, e1);
    long nfp = geo.n[d] + 1, nf = nfp * geo.n[e0] * geo.n[e1];
    long T = blockIdx.x * (long)blockDim.x + threadIdx.x;
    if (T >= nf * geo.nb) return;
    if (mode != 0) { TH[T] = mode == 1 ? R(1) : R(0); return; }
    long bb = T / nf, t = T % nf;
    long f = t % nfp, rest = t / nfp;
    long b = rest % geo.n[e1], a = rest / geo.n[e1];
    int c[3]; c[d] = geo.g[d] - 1 + (int)f; c[e0] = geo.g[e0] + (int)a; c[e1] = geo.g[e1] + (int)b;
    long base = bl[bb] * geo.ncell + c[0] * geo.stride[0] + c[1] * geo.stride[1] + c[2] * geo.stride[2], sd = geo.stride[d];
    const int N = S::N;
    R sv[6][4], wc[MAXN];
    int ns = 0;
    for (int s = 0; s < 6; s++) {
        for (int k = 0; k < N; k++) wc[k] = w[k * geo.ntot + base + (long)(s - 2) * sd];
        ns = sys.sensor(wc, sv[s]);
    }
    R r = 0;
    for (int k = 0; k < ns; k++) {
        R d1 = rabs(sv[0][k] - R(4) * sv[1][k] + R(6) * sv[2][k] - R(4) * sv[3][k] + sv[4][k]);
        R d2 = rabs(sv[1][k] - R(4) * sv[2][k] + R(6) * sv[3][k] - R(4) * sv[4][k] + sv[5][k]);
        r = rmax(r, rmax(d1, d2));
    }
    R x = rmin(rmax((r / tol - R(0.5)) / R(0.5), R(0)), R(1));
    TH[T] = fin(r) ? x * x * (R(3) - R(2) * x) : R(1);
}

// K2: midpoint Riemann fluxes, only where a face within ±2 has θ > 0.
// MF[q*NM + mi], MU[mi]; mi indexes (block, a, b, m) with m = 0..n+4.
template <class S, typename R>
__device__ void midflux_impl(const S& sys, const R* __restrict__ w, Geo geo, const int* __restrict__ bl, int d,
                             const R* __restrict__ TH, R* __restrict__ MF, R* __restrict__ MU, unsigned long long* counters) {
    const int N = S::N;
    int e0, e1; pencil_axes(d, e0, e1);
    long nfp = geo.n[d] + 1, nmp = geo.n[d] + 5;
    long nf = nfp * geo.n[e0] * geo.n[e1], nm = nmp * geo.n[e0] * geo.n[e1];
    long T = blockIdx.x * (long)blockDim.x + threadIdx.x;
    long NM = nm * geo.nb;
    if (T >= NM) return;
    long bb = T / nm, t = T % nm;
    long m = t % nmp, rest = t / nmp;
    long f0 = m - 2; // the face at this midpoint
    // needed by faces f0−2 .. f0+2
    bool need = false;
    for (long ff = f0 - 2; ff <= f0 + 2; ff++)
        if (ff >= 0 && ff < nfp && TH[bb * nf + rest * nfp + ff] > R(0)) need = true;
    if (!need) return;
    long b = rest % geo.n[e1], a = rest / geo.n[e1];
    int c[3]; c[d] = geo.g[d] - 1 + (int)f0; c[e0] = geo.g[e0] + (int)a; c[e1] = geo.g[e1] + (int)b;
    long base = bl[bb] * geo.ncell + c[0] * geo.stride[0] + c[1] * geo.stride[1] + c[2] * geo.stride[2], sd = geo.stride[d];
    R W[6][MAXN];
    for (int s = 0; s < 6; s++)
        for (int k = 0; k < N; k++) W[s][k] = w[k * geo.ntot + base + (long)(s - 2) * sd];
    R fv[MAXN];
    int fb = 0;
    R um = midpoint_flux(sys, W, 2, d, fv, &fb);
    for (int k = 0; k < N; k++) MF[k * NM + T] = fv[k];
    MU[T] = um;
    if (fb) atomicAdd(&counters[1], (unsigned long long)fb);
}

// K3: face flux = (1 − θ) central + θ WCNS combination, positivity limiter.
template <class S, typename R>
__device__ void flux_impl(const S& sys, const R* __restrict__ w, Geo geo, const int* __restrict__ bl, int d, R lambda,
                          const R* __restrict__ TH, const R* __restrict__ MF, const R* __restrict__ MU, R* __restrict__ F,
                          R* __restrict__ UF, unsigned long long* counters) {
    const int N = S::N;
    int e0, e1; pencil_axes(d, e0, e1);
    long nfp = geo.n[d] + 1, nmp = geo.n[d] + 5;
    long nf = nfp * geo.n[e0] * geo.n[e1], nm = nmp * geo.n[e0] * geo.n[e1];
    long NF = nf * geo.nb, NM = nm * geo.nb;
    long T = blockIdx.x * (long)blockDim.x + threadIdx.x;
    if (T >= NF) return;
    long bb = T / nf, t = T % nf;
    long f = t % nfp, rest = t / nfp;
    long b = rest % geo.n[e1], a = rest / geo.n[e1];
    int c[3]; c[d] = geo.g[d] - 1 + (int)f; c[e0] = geo.g[e0] + (int)a; c[e1] = geo.g[e1] + (int)b;
    long base = bl[bb] * geo.ncell + c[0] * geo.stride[0] + c[1] * geo.stride[1] + c[2] * geo.stride[2], sd = geo.stride[d];
    int un = S::VI + d;
    R theta = TH[T];
    R W[6][MAXN]; // cells i−2 .. i+3 (face between W[2] and W[3])
    for (int s = 0; s < 6; s++)
        for (int k = 0; k < N; k++) W[s][k] = w[k * geo.ntot + base + (long)(s - 2) * sd];
    R out[MAXN], uf = 0;
    for (int k = 0; k < N; k++) out[k] = 0;
    if (theta < R(1)) {
        const R A6[3] = {R(0.75), R(-0.15), R(1) / R(60)};
        R ft[MAXN];
        for (int li = 0; li < 3; li++) {
            int lag = li + 1;
            for (int mm = 0; mm < lag; mm++) {
                const R* wa = W[2 - mm];
                const R* wb = W[2 - mm + lag];
                sys.split(wa, wb, d, ft);
                for (int k = 0; k < N; k++) out[k] += R(2) * A6[li] * ft[k];
                uf += A6[li] * (wa[un] + wb[un]);
            }
        }
        for (int k = 0; k < N; k++) out[k] *= (R(1) - theta);
        uf *= (R(1) - theta);
    }
    if (theta > R(0)) {
        const R S6[3] = {R(75) / R(64) - R(25) / R(384) + R(3) / R(640), R(-25) / R(384) + R(3) / R(640), R(3) / R(640)};
        long m0 = bb * nm + rest * nmp + f + 2; // midpoint of this face
        for (int s = -2; s <= 2; s++) {
            R wgt = S6[s < 0 ? -s : s];
            for (int k = 0; k < N; k++) out[k] += theta * wgt * MF[k * NM + m0 + s];
            uf += theta * wgt * MU[m0 + s];
        }
        atomicAdd(&counters[0], 1ull);
    }
    if (lambda > R(0)) {
        R ul[MAXN], ur[MAXN], fl[MAXN], fr[MAXN], tmp[MAXN], tw[MAXN];
        sys.to_cons(W[2], ul); sys.to_cons(W[3], ur);
        sys.flux(W[2], d, fl); sys.flux(W[3], d, fr);
        bool finite = true;
        for (int k = 0; k < N; k++) finite = finite && fin(out[k]);
        auto ok = [&](const R* fx) -> bool {
            for (int k = 0; k < N; k++) tmp[k] = ul[k] - R(2) * lambda * (fx[k] - fl[k]);
            sys.to_prim(tmp, tw);
            if (!sys.admissible(tw)) return false;
            for (int k = 0; k < N; k++) tmp[k] = ur[k] + R(2) * lambda * (fx[k] - fr[k]);
            sys.to_prim(tmp, tw);
            return sys.admissible(tw);
        };
        if (!(finite && ok(out))) {
            R alpha = rmax(rabs(W[2][un]) + sys.c(W[2]), rabs(W[3][un]) + sys.c(W[3]));
            R flf[MAXN], blf[MAXN];
            for (int k = 0; k < N; k++) flf[k] = R(0.5) * (fl[k] + fr[k]) - R(0.5) * alpha * (ur[k] - ul[k]);
            R lo = 0, hi = 1;
            for (int it = 0; it < 30; it++) {
                R mid = R(0.5) * (lo + hi);
                for (int k = 0; k < N; k++) blf[k] = fin(out[k]) ? mid * out[k] + (R(1) - mid) * flf[k] : flf[k];
                if (ok(blf)) lo = mid; else hi = mid;
            }
            for (int k = 0; k < N; k++) out[k] = fin(out[k]) ? lo * out[k] + (R(1) - lo) * flf[k] : flf[k];
            R umean = R(0.5) * (W[2][un] + W[3][un]);
            uf = fin(uf) ? lo * uf + (R(1) - lo) * umean : umean;
            atomicAdd(&counters[2], 1ull);
        }
    }
    for (int k = 0; k < N; k++) F[k * NF + T] = out[k];
    UF[T] = uf;
}

// rhs += −(F[f+1] − F[f])/dx + nc_coef (UF[f+1] − UF[f])/dx along d
template <class S, typename R>
__device__ void accum_impl(const S& sys, const R* __restrict__ w, Geo geo, const int* __restrict__ bl, int d, R dx,
                           const R* __restrict__ F, const R* __restrict__ UF, R* __restrict__ rhs, int first) {
    const int N = S::N;
    long nint = (long)geo.n[0] * geo.n[1] * geo.n[2];
    long T = blockIdx.x * (long)blockDim.x + threadIdx.x;
    if (T >= nint * geo.nb) return;
    long bb = T / nint, t = T % nint;
    int ci[3];
    ci[2] = t % geo.n[2]; ci[1] = (t / geo.n[2]) % geo.n[1]; ci[0] = t / ((long)geo.n[1] * geo.n[2]);
    int e0 = (d + 1) % 3, e1 = (d + 2) % 3;
    if (e0 > e1) { int tt = e0; e0 = e1; e1 = tt; }
    long nfp = geo.n[d] + 1;
    long nf = nfp * geo.n[e0] * geo.n[e1], NF = nf * geo.nb;
    long fi = bb * nf + ((long)ci[e0] * geo.n[e1] + ci[e1]) * nfp + ci[d];
    long cell = bl[bb] * geo.ncell + (ci[0] + geo.g[0]) * geo.stride[0] + (ci[1] + geo.g[1]) * geo.stride[1] + (ci[2] + geo.g[2]) * geo.stride[2];
    R inv = R(1) / dx;
    R add[MAXN];
    for (int k = 0; k < N; k++) add[k] = -(F[k * NF + fi + 1] - F[k * NF + fi]) * inv;
    if (sys.has_nc()) {
        R wc[MAXN];
        for (int k = 0; k < N; k++) wc[k] = w[k * geo.ntot + cell];
        add[sys.nc_index()] += sys.nc_coef(wc) * (UF[fi + 1] - UF[fi]) * inv;
    }
    long slot = geo.s0 + bb;
    for (int k = 0; k < N; k++) {
        long o = k * geo.rtot + slot * nint + t;
        rhs[o] = first ? add[k] : rhs[o] + add[k];
    }
}

// u ← α u0 + β (u + dt rhs) on interior cells
template <int N, typename R>
__device__ void rk_impl(R* __restrict__ u, const R* __restrict__ u0, const R* __restrict__ rhs, Geo geo, const int* __restrict__ bl,
                        R al, R be, R dt) {
    long nint = (long)geo.n[0] * geo.n[1] * geo.n[2];
    long T = blockIdx.x * (long)blockDim.x + threadIdx.x;
    if (T >= nint * geo.nb) return;
    long bb = T / nint, t = T % nint;
    int ci[3];
    ci[2] = t % geo.n[2]; ci[1] = (t / geo.n[2]) % geo.n[1]; ci[0] = t / ((long)geo.n[1] * geo.n[2]);
    long cell = bl[bb] * geo.ncell + (ci[0] + geo.g[0]) * geo.stride[0] + (ci[1] + geo.g[1]) * geo.stride[1] + (ci[2] + geo.g[2]) * geo.stride[2];
    long slot = geo.s0 + bb;
    for (int k = 0; k < N; k++) {
        long o = k * geo.ntot + cell;
        u[o] = al * u0[o] + be * (u[o] + dt * rhs[k * geo.rtot + slot * nint + t]);
    }
}

template <class S, typename R>
__device__ void prim_impl(const S& sys, const R* __restrict__ u, R* __restrict__ w, Geo geo, const int* __restrict__ bl) {
    const int N = S::N;
    long T = blockIdx.x * (long)blockDim.x + threadIdx.x;
    if (T >= geo.ncell * geo.nb) return;
    long c = bl[T / geo.ncell] * geo.ncell + T % geo.ncell;
    R uu[MAXN], ww[MAXN];
    for (int k = 0; k < N; k++) uu[k] = u[k * geo.ntot + c];
    sys.to_prim(uu, ww);
    for (int k = 0; k < N; k++) w[k * geo.ntot + c] = ww[k];
}

// Non-reflecting ghost state (mirrors compressible::boundary::
// characteristic_ghost_layer): outgoing characteristic amplitudes from the
// boundary cell ub (linearly extrapolated with its neighbour ub1 to m
// cells beyond), incoming ones from the far field uf (all conserved).
template <class S, typename R>
__device__ void char_layer(const S& sys, const R* ub, const R* ub1, const R* uf, int d, int side, R m, R* out) {
    const int N = S::N;
    R wi[MAXN], wf[MAXN], w1[MAXN], L[MAXN][MAXN], Rm[MAXN][MAXN], lam[MAXN], dw[MAXN], wg[MAXN];
    sys.to_prim(ub, wi); sys.to_prim(uf, wf); sys.to_prim(ub1, w1);
    sys.eigen(wi, d, L, Rm);
    sys.eigenvalues(wi, d, lam);
    R outward = side == 1 ? R(1) : R(-1);
    for (int k = 0; k < N; k++) { dw[k] = wi[k] - wf[k] + m * (wi[k] - w1[k]); wg[k] = wf[k]; }
    for (int kf = 0; kf < N; kf++) {
        if (lam[kf] * outward > R(0)) {
            R amp = 0;
            for (int q = 0; q < N; q++) amp += L[kf][q] * dw[q];
            for (int q = 0; q < N; q++) wg[q] += Rm[q][kf] * amp;
        }
    }
    sys.to_cons(wg, out);
}

// ghost cells along axis d: bc 0 periodic, 1 outflow, 2 reflect,
// 3 dirichlet, 4 characteristic (far state in dir_*)
template <class S, typename R>
__device__ void ghost_impl(const S& sys, R* __restrict__ u, Geo geo, int d, int bc_lo, int bc_hi, int mom, const R* __restrict__ dir_lo,
                           const R* __restrict__ dir_hi) {
    const int N = S::N;
    int e0 = (d + 1) % 3, e1 = (d + 2) % 3;
    if (e0 > e1) { int t = e0; e0 = e1; e1 = t; }
    int ng = geo.g[d];
    long np = (long)geo.m[e0] * geo.m[e1];
    long t = blockIdx.x * (long)blockDim.x + threadIdx.x;
    if (t >= np * ng) return;
    int gi = t % ng;
    long p = t / ng;
    int b = p % geo.m[e1], a = p / geo.m[e1];
    int n = geo.n[d];
    auto at = [&](int cc) -> long { int c[3]; c[d] = cc; c[e0] = a; c[e1] = b; return c[0] * geo.stride[0] + c[1] * geo.stride[1] + c[2] * geo.stride[2]; };
    long lo_dst = at(gi), hi_dst = at(n + ng + gi);
    R cl[MAXN], ch[MAXN];
    if (bc_lo == 4) {
        R ub[MAXN], ub1[MAXN];
        for (int k = 0; k < N; k++) { ub[k] = u[k * geo.ntot + at(ng)]; ub1[k] = u[k * geo.ntot + at(ng + 1)]; }
        char_layer(sys, ub, ub1, dir_lo, d, 0, R(ng - gi), cl);
    }
    if (bc_hi == 4) {
        R ub[MAXN], ub1[MAXN];
        for (int k = 0; k < N; k++) { ub[k] = u[k * geo.ntot + at(n + ng - 1)]; ub1[k] = u[k * geo.ntot + at(n + ng - 2)]; }
        char_layer(sys, ub, ub1, dir_hi, d, 1, R(gi + 1), ch);
    }
    for (int k = 0; k < N; k++) {
        long o = (long)k * geo.ntot;
        R vl, vh;
        switch (bc_lo) {
            case 0: vl = u[o + at(n + gi)]; break;
            case 1: vl = u[o + at(ng)]; break;
            case 2: vl = u[o + at(2 * ng - 1 - gi)] * (k == mom ? R(-1) : R(1)); break;
            case 4: vl = cl[k]; break;
            default: vl = dir_lo[k];
        }
        switch (bc_hi) {
            case 0: vh = u[o + at(ng + gi)]; break;
            case 1: vh = u[o + at(n + ng - 1)]; break;
            case 2: vh = u[o + at(n + ng - 1 - gi)] * (k == mom ? R(-1) : R(1)); break;
            case 4: vh = ch[k]; break;
            default: vh = dir_hi[k];
        }
        u[o + lo_dst] = vl;
        u[o + hi_dst] = vh;
    }
}

// per-block maxima of the CFL rate Σ_d (|u_d|+c)/dx_d and of the
// dilatation stiffness rate
template <class S, typename R>
__device__ void rate_impl(const S& sys, const R* __restrict__ w, Geo geo, const int* __restrict__ bl, R dx0, R dx1, R dx2,
                          R* __restrict__ out) {
    const int N = S::N;
    __shared__ R s1[256], s2[256];
    long nint = (long)geo.n[0] * geo.n[1] * geo.n[2];
    R r1 = 0, r2 = 0;
    R dxs[3] = {dx0, dx1, dx2};
    for (long T = blockIdx.x * (long)blockDim.x + threadIdx.x; T < nint * geo.nb; T += (long)gridDim.x * blockDim.x) {
        long bb = T / nint, t = T % nint;
        int ci[3];
        ci[2] = t % geo.n[2]; ci[1] = (t / geo.n[2]) % geo.n[1]; ci[0] = t / ((long)geo.n[1] * geo.n[2]);
        long cell = bl[bb] * geo.ncell + (ci[0] + geo.g[0]) * geo.stride[0] + (ci[1] + geo.g[1]) * geo.stride[1] + (ci[2] + geo.g[2]) * geo.stride[2];
        R wc[MAXN];
        for (int k = 0; k < N; k++) wc[k] = w[k * geo.ntot + cell];
        R cs = sys.c(wc), s = 0, divu = 0;
        for (int d = 0; d < 3; d++) {
            if (geo.n[d] <= 1) continue;
            int vi = S::VI + d;
            s += (rabs(wc[vi]) + cs) / dxs[d];
            divu += rabs(w[vi * geo.ntot + cell + geo.stride[d]] - w[vi * geo.ntot + cell - geo.stride[d]]) / (R(2) * dxs[d]);
        }
        r1 = rmax(r1, s);
        r2 = rmax(r2, sys.dilatation_rate(wc, divu));
    }
    s1[threadIdx.x] = r1; s2[threadIdx.x] = r2;
    __syncthreads();
    for (int h = blockDim.x / 2; h > 0; h >>= 1) {
        if (threadIdx.x < h) { s1[threadIdx.x] = rmax(s1[threadIdx.x], s1[threadIdx.x + h]); s2[threadIdx.x] = rmax(s2[threadIdx.x], s2[threadIdx.x + h]); }
        __syncthreads();
    }
    if (threadIdx.x == 0) { out[2 * blockIdx.x] = s1[0]; out[2 * blockIdx.x + 1] = s2[0]; }
}

// ---------------------------------------------------------------- adaptive grid
// Mirrors fluid/crates/compressible/src/amr.rs. The block pool holds every
// block (leaves and ancestors) with face ghost slabs of depth ng; block p,
// interior cell c sits at u[q * ntot + p * bcell + Σ (c_d + ng) · stride_d].
// Block lookup is a dense per-level index (block coordinates → pool id).
struct AmrT {
    const int* blevel;  // pool id → level
    const int* bidx;    // pool id → block coordinates (3 each)
    const long* mapoff; // level → offset into mapd
    const int* mapd;    // dense block index, −1 = absent
    const void* far;    // [axis][side][N] Dirichlet / far-field conserved states (R)
    long bcell;         // cells per pool block (with ghosts)
    long ntot;          // variable stride of the pool
    int nb0[3];         // level-0 blocks per side
    int act[3];         // active axes
    int bd[3];          // block interior extents
    int ext[3];         // block extents with ghosts
    int mom[3];         // momentum index per axis
    int bc[3][2];       // 0 periodic 1 outflow 2 reflect 3 dirichlet 4 characteristic
    int ng;
    int levels;
    int pad_;
};

__device__ __forceinline__ long fdiv(long a, long b) { long q = a / b; return (a % b != 0 && ((a < 0) != (b < 0))) ? q - 1 : q; }
__device__ __forceinline__ long fmod_(long a, long b) { long r = a % b; return r < 0 ? r + b : r; }
__device__ __forceinline__ long amr_n(const AmrT& T, int l, int d) { return T.act[d] ? ((long)T.nb0[d] * T.bd[d]) << l : 1; }

// pool offset of interior cell c (block-local, may be negative / beyond for ghosts) of block p
__device__ __forceinline__ long amr_off(const AmrT& T, long p, const long* c) {
    long o = 0;
    for (int d = 0; d < 3; d++) o = o * T.ext[d] + (T.act[d] ? c[d] + T.ng : 0);
    return p * T.bcell + o;
}

// pool offset of in-domain (already wrapped, so non-negative) global cell
// g at level l, −1 if no block (32-bit division: 64-bit is slow on GPUs)
__device__ long amr_find(const AmrT& T, int l, const long* g) {
    long s = 0, c[3];
    for (int d = 0; d < 3; d++) {
        int nbk = T.act[d] ? (T.nb0[d] << l) : 1;
        int gi = (int)g[d];
        int k = T.act[d] ? gi / T.bd[d] : 0;
        if (k < 0 || k >= nbk) return -1;
        s = s * nbk + k;
        c[d] = T.act[d] ? gi - k * T.bd[d] : 0;
    }
    int p = T.mapd[T.mapoff[l] + s];
    return p < 0 ? -1 : amr_off(T, p, c);
}

// value at global cell g of level l with periodic wrap; false if no block there
template <int N, typename R> __device__ bool amr_fetch(const AmrT& T, const R* __restrict__ u, int l, long* g, R* v) {
    for (int d = 0; d < 3; d++) {
        if (!T.act[d]) { g[d] = 0; continue; }
        long n = amr_n(T, l, d);
        if (g[d] < 0 || g[d] >= n) {
            if (T.bc[d][g[d] >= n] != 0) return false;
            g[d] = g[d] < 0 ? g[d] + n : g[d] - n; // ghost depths are far below n
        }
    }
    long o = amr_find(T, l, g);
    if (o < 0) return false;
    for (int k = 0; k < N; k++) v[k] = u[k * T.ntot + o];
    return true;
}

// `Amr::sample` restricted to cells that need no prediction: boundary
// conditions applied outside the domain (characteristic axes nested with
// the last axis outermost, like the axis sweeps of `bc_box`)
template <class S, typename R> __device__ bool amr_sample(const S& sys, const AmrT& T, const R* __restrict__ u, int l, const long* g0, R* out) {
    const int N = S::N;
    {
        // fast path: inside the domain or across periodic boundaries only
        bool plain = true;
        for (int d = 0; d < 3; d++) {
            long n = amr_n(T, l, d);
            if (T.act[d] && (g0[d] < 0 || g0[d] >= n) && T.bc[d][g0[d] >= n] != 0) plain = false;
        }
        if (plain) {
            long gg[3] = {g0[0], g0[1], g0[2]};
            return amr_fetch<N>(T, u, l, gg, out);
        }
    }
    const R* far = (const R*)T.far;
    long g[3] = {g0[0], g0[1], g0[2]};
    bool flip[3] = {false, false, false};
    int cd[3], cs[3], ncd = 0;
    for (int d = 2; d >= 0; d--) {
        if (!T.act[d]) { g[d] = 0; continue; }
        long n = amr_n(T, l, d);
        if (g[d] >= 0 && g[d] < n) continue;
        int side = g[d] >= n;
        switch (T.bc[d][side]) {
            case 0: g[d] = fmod_(g[d], n); break;
            case 1: g[d] = g[d] < 0 ? 0 : n - 1; break;
            case 2: g[d] = side == 0 ? -1 - g[d] : 2 * n - 1 - g[d]; flip[d] = true; break;
            case 3: for (int k = 0; k < N; k++) out[k] = far[(d * 2 + side) * N + k]; return true;
            default: cd[ncd] = d; cs[ncd] = side; ncd++;
        }
    }
    // cd is in decreasing axis order; the innermost (combined first) is the smallest axis
    R leaf[8][MAXN];
    for (int s = 0; s < (1 << ncd); s++) {
        long gg[3] = {g[0], g[1], g[2]};
        for (int i = 0; i < ncd; i++) {
            int d = cd[i];
            long n = amr_n(T, l, d);
            long gi = g[d] < 0 ? 0 : n - 1;
            if ((s >> i) & 1) gi += cs[i] == 1 ? -1 : 1;
            gg[d] = gi;
        }
        if (!amr_fetch<N>(T, u, l, gg, leaf[s])) return false;
        for (int d = 0; d < 3; d++)
            if (flip[d]) leaf[s][T.mom[d]] = -leaf[s][T.mom[d]];
    }
    for (int i = ncd - 1; i >= 0; i--) {
        int d = cd[i];
        long n = amr_n(T, l, d);
        R m = R(g[d] < 0 ? -g[d] : g[d] - (n - 1));
        R tmp[MAXN];
        for (int s = 0; s < (1 << ncd); s++) {
            // bits above i (smaller axes) are already combined: pair the
            // entries that differ in bit i and have those bits clear
            if ((s >> i) & 1) continue;
            bool lower_clear = true;
            for (int j = i + 1; j < ncd; j++) if ((s >> j) & 1) lower_clear = false;
            if (!lower_clear) continue;
            char_layer(sys, leaf[s], leaf[s | (1 << i)], far + (d * 2 + cs[i]) * N, d, cs[i], m, tmp);
            for (int k = 0; k < N; k++) leaf[s][k] = tmp[k];
        }
    }
    for (int k = 0; k < N; k++) out[k] = leaf[0][k];
    return true;
}

// Face ghost slabs along axis d of the leaves bl[0..nb): same-level block →
// copy, coarser region → 6th-order prediction from level l − 1 (6^D
// points, weights wp[parity][k]), outside a non-periodic boundary → the
// boundary condition. Slab cells only (no edges / corners: the scheme is
// dimension-by-dimension and never reads them).
template <class S, typename R>
__device__ void amr_ghost_impl(const S& sys, AmrT T, R* __restrict__ u, const int* __restrict__ bl, int nb, int d, const R* __restrict__ wp,
                               unsigned long long* err, int sep) {
    const int N = S::N;
    int e0, e1; pencil_axes(d, e0, e1);
    long per = 2L * T.ng * T.bd[e0] * T.bd[e1];
    long t = blockIdx.x * (long)blockDim.x + threadIdx.x;
    if (t >= per * nb) return;
    long bb = t / per, r = t % per;
    long b = r % T.bd[e1]; r /= T.bd[e1];
    long a = r % T.bd[e0]; r /= T.bd[e0];
    int gi = r % T.ng, side = r / T.ng;
    int p = bl[bb], l = T.blevel[p];
    long lc[3];
    lc[d] = side ? T.bd[d] + gi : gi - T.ng;
    lc[e0] = a; lc[e1] = b;
    long g[3];
    for (int e = 0; e < 3; e++) g[e] = T.act[e] ? (long)T.bidx[3 * p + e] * T.bd[e] + lc[e] : 0;
    R v[MAXN];
    bool ok;
    long n = amr_n(T, l, d);
    if ((g[d] < 0 || g[d] >= n) && T.bc[d][side] != 0) {
        ok = amr_sample(sys, T, u, l, g, v);
    } else {
        long gw[3] = {g[0], g[1], g[2]};
        ok = amr_fetch<N>(T, u, l, gw, v);
        if (!ok && l > 0 && sep) return; // filled by the separable kernels
        if (!ok && l > 0) {
            // prediction from the coarser level (index of the wrapped cell)
            for (int k = 0; k < N; k++) v[k] = 0;
            int dims[3], nd = 0;
            for (int e = 0; e < 3; e++) if (T.act[e]) dims[nd++] = e;
            int total = 1;
            for (int i = 0; i < nd; i++) total *= 6;
            ok = true;
            for (int s = 0; s < total && ok; s++) {
                R wgt = 1;
                long c[3] = {0, 0, 0};
                int rem = s;
                for (int i = 0; i < nd; i++) {
                    int e = dims[i], k = rem % 6;
                    rem /= 6;
                    int odd = (int)(gw[e] & 1);
                    c[e] = (gw[e] >> 1) + (odd ? -2 : -3) + k;
                    wgt *= wp[odd * 6 + k];
                }
                R cv[MAXN];
                ok = amr_sample(sys, T, u, l - 1, c, cv);
                for (int k = 0; k < N; k++) v[k] += wgt * cv[k];
            }
        }
    }
    if (!ok) {
        atomicAdd(err, 1ull);
        for (int k = 0; k < N; k++) v[k] = R(0) / R(0);
    }
    long o = amr_off(T, p, lc);
    for (int k = 0; k < N; k++) u[k * T.ntot + o] = v[k];
}

// Separable prediction of the face ghost slabs that lie in a coarser
// region (as the coarse-patch path of `Amr::fill_scratch`): region r =
// (pool block, axis a = launch axis, side) with fine origin f0 and coarse
// patch origin cmin (8 ints: p, side, f0[3], cmin[3]). Gather the coarse
// patch (szc), then one 6-point pass per active axis in increasing order.
template <class S, typename R>
__device__ void amr_pgather_impl(const S& sys, AmrT T, const R* __restrict__ u, const int* __restrict__ reg, int nreg, int c0, int c1, int c2,
                                 R* __restrict__ A, unsigned long long* err) {
    const int N = S::N;
    long szc[3] = {c0, c1, c2}, P = szc[0] * szc[1] * szc[2];
    long t = blockIdx.x * (long)blockDim.x + threadIdx.x;
    if (t >= P * nreg) return;
    long r = t / P, q = t % P;
    long c[3];
    c[2] = q % szc[2]; c[1] = (q / szc[2]) % szc[1]; c[0] = q / (szc[1] * szc[2]);
    const int* rg = reg + 8 * r;
    int p = rg[0], l = T.blevel[p];
    long g[3];
    for (int d = 0; d < 3; d++) g[d] = T.act[d] ? rg[5 + d] + c[d] : 0;
    R v[MAXN];
    if (!amr_sample(sys, T, u, l - 1, g, v)) {
        atomicAdd(err, 1ull);
        for (int k = 0; k < N; k++) v[k] = R(0) / R(0);
    }
    long PT = P * nreg;
    for (int k = 0; k < N; k++) A[k * PT + t] = v[k];
}

template <int N, typename R>
__device__ void amr_ppass_impl(AmrT T, const int* __restrict__ reg, int nreg, int a, int i0, int i1, int i2, int o_a, const R* __restrict__ A,
                               R* __restrict__ B, const R* __restrict__ wp, int fin_, R* __restrict__ u) {
    long szi[3] = {i0, i1, i2}, szo[3] = {i0, i1, i2};
    szo[a] = o_a;
    long Pi = szi[0] * szi[1] * szi[2], Po = szo[0] * szo[1] * szo[2];
    long t = blockIdx.x * (long)blockDim.x + threadIdx.x;
    if (t >= Po * nreg) return;
    long r = t / Po, q = t % Po;
    long c[3];
    c[2] = q % szo[2]; c[1] = (q / szo[2]) % szo[1]; c[0] = q / (szo[1] * szo[2]);
    const int* rg = reg + 8 * r;
    long f = rg[2 + a] + c[a];
    int par = (int)(f & 1);
    long base = (f >> 1) + (par ? -2 : -3) - rg[5 + a];
    long PTi = Pi * nreg;
    R acc[MAXN];
    for (int k = 0; k < N; k++) acc[k] = 0;
    for (int kk = 0; kk < 6; kk++) {
        long s2[3] = {c[0], c[1], c[2]};
        s2[a] = base + kk;
        long sp = r * Pi + (s2[0] * szi[1] + s2[1]) * szi[2] + s2[2];
        R wgt = wp[par * 6 + kk];
        for (int k = 0; k < N; k++) acc[k] += wgt * A[k * PTi + sp];
    }
    if (fin_) {
        int p = rg[0];
        long lc[3];
        for (int d = 0; d < 3; d++) lc[d] = T.act[d] ? rg[2 + d] + c[d] - (long)T.bidx[3 * p + d] * T.bd[d] : 0;
        long o = amr_off(T, p, lc);
        for (int k = 0; k < N; k++) u[k * T.ntot + o] = acc[k];
    } else {
        long PTo = Po * nreg;
        for (int k = 0; k < N; k++) B[k * PTo + t] = acc[k];
    }
}

// Largest relative wavelet detail per leaf (as `Amr::block_detail`): in
// each active axis, sensor value minus its 6th-order midpoint prediction
// from the cells at ±1, ±3, ±5; one thread block per leaf.
template <class S, typename R>
__device__ void amr_detail_impl(const S& sys, Geo geo, const R* __restrict__ w, const int* __restrict__ bl, const R* __restrict__ wr,
                                R* __restrict__ det) {
    const int N = S::N;
    __shared__ R sm[256];
    long b0 = bl[blockIdx.x] * geo.ncell;
    long nint = (long)geo.n[0] * geo.n[1] * geo.n[2];
    const int offs[6] = {-5, -3, -1, 1, 3, 5};
    R mx = 0;
    for (long t = threadIdx.x; t < nint; t += blockDim.x) {
        int ci[3];
        ci[2] = t % geo.n[2]; ci[1] = (t / geo.n[2]) % geo.n[1]; ci[0] = t / ((long)geo.n[1] * geo.n[2]);
        long cell = b0 + (ci[0] + geo.g[0]) * geo.stride[0] + (ci[1] + geo.g[1]) * geo.stride[1] + (ci[2] + geo.g[2]) * geo.stride[2];
        R wc[MAXN], s0[4], sv[4];
        for (int k = 0; k < N; k++) wc[k] = w[k * geo.ntot + cell];
        int ns = sys.sensor(wc, s0);
        for (int d = 0; d < 3; d++) {
            if (geo.n[d] <= 1) continue;
            R pred[4] = {0, 0, 0, 0};
            for (int m = 0; m < 6; m++) {
                for (int k = 0; k < N; k++) wc[k] = w[k * geo.ntot + cell + offs[m] * geo.stride[d]];
                sys.sensor(wc, sv);
                for (int q = 0; q < ns; q++) pred[q] += wr[m] * sv[q];
            }
            for (int q = 0; q < ns; q++) {
                R dd = rabs(s0[q] - pred[q]);
                mx = fin(dd) ? rmax(mx, dd) : R(1) / R(0);
            }
        }
    }
    sm[threadIdx.x] = mx;
    __syncthreads();
    for (int h = blockDim.x / 2; h > 0; h >>= 1) {
        if (threadIdx.x < h) sm[threadIdx.x] = rmax(sm[threadIdx.x], sm[threadIdx.x + h]);
        __syncthreads();
    }
    if (threadIdx.x == 0) det[blockIdx.x] = sm[0];
}

// Restriction of ancestors tl[0..nt) at level l from level l + 1, as
// `Amr::restrict_box` over the whole block: gather the fine patch
// [2c₀ − 4, 2c₀ + 2B + 4) with a validity mask, then one pass per axis.
template <class S, typename R>
__device__ void amr_rgather_impl(const S& sys, AmrT T, const R* __restrict__ u, const int* __restrict__ tl, int nt, int l, R* __restrict__ A,
                                 unsigned char* __restrict__ okA) {
    const int N = S::N;
    long sz[3], P = 1;
    for (int d = 0; d < 3; d++) { sz[d] = T.act[d] ? 2 * T.bd[d] + 8 : 1; P *= sz[d]; }
    long t = blockIdx.x * (long)blockDim.x + threadIdx.x;
    if (t >= P * nt) return;
    long task = t / P, r = t % P;
    long c[3];
    c[2] = r % sz[2]; c[1] = (r / sz[2]) % sz[1]; c[0] = r / (sz[1] * sz[2]);
    int p = tl[task];
    long g[3];
    bool outside = false;
    for (int d = 0; d < 3; d++) {
        g[d] = T.act[d] ? 2L * T.bidx[3 * p + d] * T.bd[d] - 4 + c[d] : 0;
        long n = amr_n(T, l + 1, d);
        if (T.act[d] && (g[d] < 0 || g[d] >= n) && T.bc[d][g[d] >= n] != 0) outside = true;
    }
    R v[MAXN];
    bool ok = outside ? amr_sample(sys, T, u, l + 1, g, v) : amr_fetch<N>(T, u, l + 1, g, v);
    long PT = P * nt;
    for (int k = 0; k < N; k++) A[k * PT + t] = ok ? v[k] : R(0);
    okA[t] = ok;
}

// pass along axis d: input extents szi (fine along d), output szo
// (coarse along d, = B). wts: 5 shifted 6-point weight sets (shift 0, −1,
// 1, −2, 2). Final pass (fin != 0) writes the block interiors in the pool.
template <int N, typename R>
__device__ void amr_rpass_impl(AmrT T, const int* __restrict__ tl, int nt, int d, int si0, int si1, int si2, const R* __restrict__ A,
                               const unsigned char* __restrict__ okA, R* __restrict__ B, unsigned char* __restrict__ okB, const R* __restrict__ wts,
                               int fin_, R* __restrict__ u) {
    long szi[3] = {si0, si1, si2}, szo[3] = {si0, si1, si2};
    szo[d] = T.bd[d];
    long Pi = szi[0] * szi[1] * szi[2], Po = szo[0] * szo[1] * szo[2];
    long t = blockIdx.x * (long)blockDim.x + threadIdx.x;
    if (t >= Po * nt) return;
    long task = t / Po, r = t % Po;
    long c[3];
    c[2] = r % szo[2]; c[1] = (r / szo[2]) % szo[1]; c[0] = r / (szo[1] * szo[2]);
    auto src = [&](long pos) -> long { long s2[3] = {c[0], c[1], c[2]}; s2[d] = pos; return task * Pi + (s2[0] * szi[1] + s2[1]) * szi[2] + s2[2]; };
    const int shifts[5] = {0, -1, 1, -2, 2};
    long base = 2 * c[d] + 2;
    long PTi = Pi * nt;
    R acc[MAXN];
    bool done = false;
    for (int m = 0; m < 5 && !done; m++) {
        long st = base + shifts[m];
        if (st < 0 || st + 5 >= szi[d]) continue;
        bool all = true;
        for (int kk = 0; kk < 6; kk++) all = all && okA[src(st + kk)];
        if (!all) continue;
        for (int q = 0; q < N; q++) acc[q] = 0;
        for (int kk = 0; kk < 6; kk++) {
            long sp = src(st + kk);
            for (int q = 0; q < N; q++) acc[q] += wts[m * 6 + kk] * A[q * PTi + sp];
        }
        done = true;
    }
    if (!done) {
        long pa = src(base + 2), pb = src(base + 3);
        if (okA[pa] && okA[pb]) {
            for (int q = 0; q < N; q++) acc[q] = R(0.5) * (A[q * PTi + pa] + A[q * PTi + pb]);
            done = true;
        } else {
            for (int q = 0; q < N; q++) acc[q] = 0;
        }
    }
    if (fin_) {
        long cc[3] = {c[0], c[1], c[2]};
        long o = amr_off(T, tl[task], cc);
        for (int q = 0; q < N; q++) u[q * T.ntot + o] = acc[q];
    } else {
        long PTo = Po * nt;
        for (int q = 0; q < N; q++) B[q * PTo + t] = acc[q];
        okB[t] = done;
    }
}

// boundary face fluxes of a launch (one level, axis d) into the per-slot
// store BF[q * bft + ((slot * 3 + d) * 2 + side) * nfc + face]
template <int N, typename R>
__device__ void amr_faces_impl(Geo geo, int d, const R* __restrict__ F, R* __restrict__ BF, long bft, long nfc) {
    int e0, e1; pencil_axes(d, e0, e1);
    long nfp = geo.n[d] + 1, nf = nfp * geo.n[e0] * geo.n[e1], NF = nf * geo.nb;
    long per = 2L * geo.n[e0] * geo.n[e1];
    long t = blockIdx.x * (long)blockDim.x + threadIdx.x;
    if (t >= per * geo.nb) return;
    long bb = t / per, r = t % per;
    long fi = r % (geo.n[e0] * geo.n[e1]);
    int side = r / (geo.n[e0] * geo.n[e1]);
    long src = bb * nf + fi * nfp + (side ? geo.n[d] : 0);
    long dst = ((long)((geo.s0 + bb) * 3 + d) * 2 + side) * nfc + fi;
    for (int q = 0; q < N; q++) BF[q * bft + dst] = F[q * NF + src];
}

// refluxing (as `Amr::reflux`) for one (axis, side): task i corrects the
// coarse cell next to coarse face cf of slot cs by the mean of up to four
// fine faces. Launched per (axis, side) in the CPU order, so no cell is
// written twice in a launch.
template <int N, typename R>
__device__ void amr_reflux_impl(int ntask, const int* __restrict__ task, const R* __restrict__ invdx, int d, int side, const R* __restrict__ BF,
                                long bft, long nfc, R* __restrict__ rhs, long rtot, long nint) {
    long i = blockIdx.x * (long)blockDim.x + threadIdx.x;
    if (i >= ntask) return;
    const int* tk = task + 13 * i; // cs, cell, cf, level, nfine, (fine slot, fine face) × 4
    int cs = tk[0], cell = tk[1], cf = tk[2], lev = tk[3], nfn = tk[4];
    R sum[MAXN];
    for (int q = 0; q < N; q++) sum[q] = 0;
    for (int j = 0; j < nfn; j++) {
        int fs = tk[5 + 2 * j], ff = tk[6 + 2 * j];
        long o = ((long)(fs * 3 + d) * 2 + (1 - side)) * nfc + ff;
        for (int q = 0; q < N; q++) sum[q] += BF[q * bft + o];
    }
    if (nfn == 0) return;
    long oc = ((long)(cs * 3 + d) * 2 + side) * nfc + cf;
    R inv = invdx[lev * 3 + d];
    for (int q = 0; q < N; q++) {
        R df = sum[q] / R(nfn) - BF[q * bft + oc];
        rhs[q * rtot + (long)cs * nint + cell] += side == 1 ? -df * inv : df * inv;
    }
}

// ---------------------------------------------------------------- entry points
#define EULER(R) Euler<R> sys; sys.e = Nasg<R>{g1, pi1, b1, q1};
#define FIVE(R) FiveEq<R> sys; sys.e1 = Nasg<R>{g1, pi1, b1, q1}; sys.e2 = Nasg<R>{g2, pi2, b2, q2}; sys.kapila = kap;
#define EOSARGS(R) R g1, R pi1, R b1, R q1, R g2, R pi2, R b2, R q2, int kap

#define ENTRIES(NAME, R, MAKE, NV)                                                                                         \
    extern "C" __global__ void theta_##NAME(EOSARGS(R), const R* w, Geo geo, const int* bl, int d, int mode, R tol,       \
                                            R* TH) {                                                                       \
        MAKE(R) theta_impl(sys, w, geo, bl, d, mode, tol, TH);                                                             \
    }                                                                                                                       \
    extern "C" __global__ void midflux_##NAME(EOSARGS(R), const R* w, Geo geo, const int* bl, int d, const R* TH, R* MF,  \
                                              R* MU, unsigned long long* counters) {                                       \
        MAKE(R) midflux_impl(sys, w, geo, bl, d, TH, MF, MU, counters);                                                    \
    }                                                                                                                       \
    extern "C" __global__ void flux_##NAME(EOSARGS(R), const R* w, Geo geo, const int* bl, int d, R lambda, const R* TH,  \
                                           const R* MF, const R* MU, R* F, R* UF, unsigned long long* counters) {          \
        MAKE(R) flux_impl(sys, w, geo, bl, d, lambda, TH, MF, MU, F, UF, counters);                                        \
    }                                                                                                                       \
    extern "C" __global__ void accum_##NAME(EOSARGS(R), const R* w, Geo geo, const int* bl, int d, R dx, const R* F,      \
                                            const R* UF, R* rhs, int first) {                                              \
        MAKE(R) accum_impl(sys, w, geo, bl, d, dx, F, UF, rhs, first);                                                     \
    }                                                                                                                       \
    extern "C" __global__ void prim_##NAME(EOSARGS(R), const R* u, R* w, Geo geo, const int* bl) {                        \
        MAKE(R) prim_impl(sys, u, w, geo, bl);                                                                             \
    }                                                                                                                       \
    extern "C" __global__ void rate_##NAME(EOSARGS(R), const R* w, Geo geo, const int* bl, R dx0, R dx1, R dx2, R* out) { \
        MAKE(R) rate_impl(sys, w, geo, bl, dx0, dx1, dx2, out);                                                            \
    }                                                                                                                       \
    extern "C" __global__ void rk_##NAME(R* u, const R* u0, const R* rhs, Geo geo, const int* bl, R al, R be, R dt) {     \
        rk_impl<NV, R>(u, u0, rhs, geo, bl, al, be, dt);                                                                   \
    }                                                                                                                       \
    extern "C" __global__ void ghost_##NAME(EOSARGS(R), R* u, Geo geo, int d, int bc_lo, int bc_hi, int mom, const R* dlo, \
                                            const R* dhi) {                                                                \
        MAKE(R) ghost_impl(sys, u, geo, d, bc_lo, bc_hi, mom, dlo, dhi);                                                   \
    }                                                                                                                       \
    extern "C" __global__ void amr_ghost_##NAME(EOSARGS(R), AmrT T, R* u, const int* bl, int nb, int d, const R* wp,      \
                                                unsigned long long* err, int sep) {                                        \
        MAKE(R) amr_ghost_impl(sys, T, u, bl, nb, d, wp, err, sep);                                                        \
    }                                                                                                                       \
    extern "C" __global__ void amr_pgather_##NAME(EOSARGS(R), AmrT T, const R* u, const int* reg, int nreg, int c0, int c1, \
                                                  int c2, R* A, unsigned long long* err) {                                 \
        MAKE(R) amr_pgather_impl(sys, T, u, reg, nreg, c0, c1, c2, A, err);                                                \
    }                                                                                                                       \
    extern "C" __global__ void amr_ppass_##NAME(AmrT T, const int* reg, int nreg, int a, int i0, int i1, int i2, int o_a,  \
                                                const R* A, R* B, const R* wp, int fin_, R* u) {                           \
        amr_ppass_impl<NV, R>(T, reg, nreg, a, i0, i1, i2, o_a, A, B, wp, fin_, u);                                        \
    }                                                                                                                       \
    extern "C" __global__ void amr_rgather_##NAME(EOSARGS(R), AmrT T, const R* u, const int* tl, int nt, int l, R* A,     \
                                                  unsigned char* okA) {                                                    \
        MAKE(R) amr_rgather_impl(sys, T, u, tl, nt, l, A, okA);                                                            \
    }                                                                                                                       \
    extern "C" __global__ void amr_rpass_##NAME(AmrT T, const int* tl, int nt, int d, int s0, int s1, int s2, const R* A, \
                                                const unsigned char* okA, R* B, unsigned char* okB, const R* wts, int fin_, \
                                                R* u) {                                                                    \
        amr_rpass_impl<NV, R>(T, tl, nt, d, s0, s1, s2, A, okA, B, okB, wts, fin_, u);                                     \
    }                                                                                                                       \
    extern "C" __global__ void amr_detail_##NAME(EOSARGS(R), Geo geo, const R* w, const int* bl, const R* wr, R* det) {    \
        MAKE(R) amr_detail_impl(sys, geo, w, bl, wr, det);                                                                 \
    }                                                                                                                       \
    extern "C" __global__ void amr_faces_##NAME(Geo geo, int d, const R* F, R* BF, long bft, long nfc) {                  \
        amr_faces_impl<NV, R>(geo, d, F, BF, bft, nfc);                                                                    \
    }                                                                                                                       \
    extern "C" __global__ void amr_reflux_##NAME(int ntask, const int* task, const R* invdx, int d, int side, const R* BF, \
                                                 long bft, long nfc, R* rhs, long rtot, long nint) {                       \
        amr_reflux_impl<NV, R>(ntask, task, invdx, d, side, BF, bft, nfc, rhs, rtot, nint);                                \
    }

ENTRIES(euler_f64, double, EULER, 5)
ENTRIES(euler_f32, float, EULER, 5)
ENTRIES(five_f64, double, FIVE, 7)
ENTRIES(five_f32, float, FIVE, 7)
