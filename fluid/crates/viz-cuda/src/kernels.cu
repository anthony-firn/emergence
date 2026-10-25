// Sphere tracing on the device. The generated `scene_sdf` is prepended to
// this file before compilation.
//
// The march is the same algorithm as viz::trace::march, including its
// three outcomes: a hit, a clean miss, or a budget exhausted while still
// approaching. Exhaustion is written as a negative depth so the host can
// tell it from a miss, because a renderer that confuses the two puts holes
// in grazing geometry.

__device__ __forceinline__ bool slab(REAL3 o, REAL3 d, REAL3 lo, REAL3 hi,
                                     REAL tmin, REAL tmax, REAL *a, REAL *b) {
    REAL near = tmin, far = tmax;
    REAL oc[3] = {o.x, o.y, o.z};
    REAL dc[3] = {d.x, d.y, d.z};
    REAL lc[3] = {lo.x, lo.y, lo.z};
    REAL hc[3] = {hi.x, hi.y, hi.z};
    for (int k = 0; k < 3; ++k) {
        if (fabs(dc[k]) < (REAL)1e-30) {
            if (oc[k] < lc[k] || oc[k] > hc[k]) return false;
            continue;
        }
        REAL t0 = (lc[k] - oc[k]) / dc[k];
        REAL t1 = (hc[k] - oc[k]) / dc[k];
        if (t0 > t1) { REAL t = t0; t0 = t1; t1 = t; }
        near = fmax(near, t0);
        far = fmin(far, t1);
    }
    if (near > far) return false;
    *a = near;
    *b = far;
    return true;
}

extern "C" __global__ void trace(
    double3 eye_d, double3 right_d, double3 up_d, double3 forward_d,
    double half_h_d, double aspect_d, int width, int height,
    double3 bbox_lo_d, double3 bbox_hi_d, int bounded,
    double tolerance_d, int max_steps, double t_max_d, double normal_eps_d,
    double *depth, double *normals, int *steps) {
    // The camera and the march settings arrive in double and are cast
    // once; the hot loop runs at REAL.
    REAL3 eye = mk3((REAL)eye_d.x, (REAL)eye_d.y, (REAL)eye_d.z);
    REAL3 right = mk3((REAL)right_d.x, (REAL)right_d.y, (REAL)right_d.z);
    REAL3 up = mk3((REAL)up_d.x, (REAL)up_d.y, (REAL)up_d.z);
    REAL3 forward = mk3((REAL)forward_d.x, (REAL)forward_d.y, (REAL)forward_d.z);
    REAL3 bbox_lo = mk3((REAL)bbox_lo_d.x, (REAL)bbox_lo_d.y, (REAL)bbox_lo_d.z);
    REAL3 bbox_hi = mk3((REAL)bbox_hi_d.x, (REAL)bbox_hi_d.y, (REAL)bbox_hi_d.z);
    REAL half_h = (REAL)half_h_d, aspect = (REAL)aspect_d;
    REAL tolerance = (REAL)tolerance_d, t_max = (REAL)t_max_d;
    REAL normal_eps = (REAL)normal_eps_d;
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= width * height) return;
    int x = idx % width;
    int y = idx / width;

    REAL u = R(2.0) * ((REAL)x + R(0.5)) / (REAL)width - R(1.0);
    REAL v = R(1.0) - R(2.0) * ((REAL)y + R(0.5)) / (REAL)height;
    REAL3 dir = add3(forward,
                     add3(scale3(right, u * half_h * aspect), scale3(up, v * half_h)));
    dir = scale3(dir, R(1.0) / v3len(dir));

    REAL t_enter = R(1e-9), t_exit = t_max;
    if (bounded) {
        REAL a, b;
        if (!slab(eye, dir, bbox_lo, bbox_hi, R(1e-9), t_max, &a, &b)) {
            depth[idx] = (double)VINF;
            steps[idx] = 0;
            normals[3 * idx] = 0.0;
            normals[3 * idx + 1] = 0.0;
            normals[3 * idx + 2] = 0.0;
            return;
        }
        t_enter = a;
        t_exit = b;
    }

    REAL3 p = add3(eye, scale3(dir, t_enter));
    REAL first = scene_sdf(p);
    REAL sign = first < R(0.0) ? R(-1.0) : R(1.0);
    REAL t = t_enter;
    int n = 0;
    REAL d = fabs(first);
    bool hit = fabs(first) <= tolerance;
    if (!hit) {
        while (n < max_steps) {
            p = add3(eye, scale3(dir, t));
            d = sign * scene_sdf(p);
            ++n;
            if (d <= tolerance) { hit = true; break; }
            t += d;
            if (t > t_exit) {
                depth[idx] = (double)VINF;
                steps[idx] = n;
                normals[3 * idx] = 0.0;
                normals[3 * idx + 1] = 0.0;
                normals[3 * idx + 2] = 0.0;
                return;
            }
        }
    } else {
        n = 1;
    }

    if (!hit) {
        // Budget exhausted while still approaching: a negative depth marks
        // it, carrying how far the ray had travelled.
        depth[idx] = -t;
        steps[idx] = n;
        normals[3 * idx] = 0.0;
        normals[3 * idx + 1] = 0.0;
        normals[3 * idx + 2] = 0.0;
        return;
    }

    p = add3(eye, scale3(dir, t));
    REAL e = normal_eps;
    REAL3 g = mk3(
        (scene_sdf(mk3(p.x + e, p.y, p.z)) - scene_sdf(mk3(p.x - e, p.y, p.z))) / (R(2.0) * e),
        (scene_sdf(mk3(p.x, p.y + e, p.z)) - scene_sdf(mk3(p.x, p.y - e, p.z))) / (R(2.0) * e),
        (scene_sdf(mk3(p.x, p.y, p.z + e)) - scene_sdf(mk3(p.x, p.y, p.z - e))) / (R(2.0) * e));
    REAL gl = v3len(g);
    if (gl > R(0.0)) g = scale3(g, R(1.0) / gl); else g = mk3(R(0.0), R(0.0), R(1.0));
    depth[idx] = (double)t;
    steps[idx] = n;
    normals[3 * idx] = g.x;
    normals[3 * idx + 1] = g.y;
    normals[3 * idx + 2] = g.z;
}

// Evaluate the field at a list of points: the test that the generated
// device function is the same function as the host one.
extern "C" __global__ void sample_field(const double *points, double *out, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    out[i] = (double)scene_sdf(mk3((REAL)points[3 * i], (REAL)points[3 * i + 1],
                                   (REAL)points[3 * i + 2]));
}
