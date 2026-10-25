// Vector helpers and exact primitive distances, transcribed from
// cad/src/solid.rs. The arithmetic is written in the same order as the
// Rust so that the two agree bit for bit in double precision; the tests
// check that over tens of thousands of points.

// nvrtc compiles without the standard headers, so the infinity the field
// algebra needs is spelled out from its bit pattern.
#ifndef INFINITY
#define INFINITY __longlong_as_double(0x7ff0000000000000LL)
#endif

__device__ __forceinline__ double3 sub3(double3 a, double3 b) {
    return make_double3(a.x - b.x, a.y - b.y, a.z - b.z);
}
__device__ __forceinline__ double3 add3(double3 a, double3 b) {
    return make_double3(a.x + b.x, a.y + b.y, a.z + b.z);
}
__device__ __forceinline__ double3 scale3(double3 a, double s) {
    return make_double3(a.x * s, a.y * s, a.z * s);
}
__device__ __forceinline__ double dot3(double3 a, double3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
__device__ __forceinline__ double v3len(double3 a) { return sqrt(dot3(a, a)); }

__device__ __forceinline__ double2 sub2(double2 a, double2 b) {
    return make_double2(a.x - b.x, a.y - b.y);
}
__device__ __forceinline__ double v2len(double2 a) {
    return sqrt(a.x * a.x + a.y * a.y);
}

// Exact capped cylinder (cad::Solid::Cylinder).
__device__ double sd_cylinder(double3 p, double3 a, double3 b, double r) {
    double3 ba = sub3(b, a);
    double3 pa = sub3(p, a);
    double baba = dot3(ba, ba);
    double paba = dot3(pa, ba);
    if (baba == 0.0) return v3len(pa) - r;
    double x = v3len(sub3(scale3(pa, baba), scale3(ba, paba))) - r * baba;
    double y = fabs(paba - 0.5 * baba) - 0.5 * baba;
    double x2 = x * x;
    double y2 = y * y * baba;
    double d;
    if (fmax(x, y) < 0.0) {
        d = -fmin(x2, y2);
    } else {
        d = (x > 0.0 ? x2 : 0.0) + (y > 0.0 ? y2 : 0.0);
    }
    return copysign(sqrt(fabs(d)), d) / baba;
}

// Exact capped cone (cad::Solid::Cone).
__device__ double sd_cone(double3 p, double3 a, double3 b, double ra, double rb) {
    double3 ba = sub3(b, a);
    double baba = dot3(ba, ba);
    if (baba == 0.0) return v3len(sub3(p, a)) - fmax(ra, rb);
    double3 pa = sub3(p, a);
    double papa = dot3(pa, pa);
    double paba = dot3(pa, ba) / baba;
    double x = sqrt(fmax(papa - paba * paba * baba, 0.0));
    double cax = fmax(0.0, x - (paba < 0.5 ? ra : rb));
    double cay = fabs(paba - 0.5) - 0.5;
    double rba = rb - ra;
    double k = rba * rba + baba;
    double fr = fmin(fmax((rba * (x - ra) + paba * baba) / k, 0.0), 1.0);
    double cbx = x - ra - fr * rba;
    double cby = paba - fr;
    double s = (cbx < 0.0 && cay < 0.0) ? -1.0 : 1.0;
    return s * sqrt(fmin(cax * cax + cay * cay * baba, cbx * cbx + cby * cby * baba));
}

// Capsule (cad::Solid::Capsule).
__device__ double sd_capsule(double3 p, double3 a, double3 b, double r) {
    double3 ba = sub3(b, a);
    double ee = dot3(ba, ba);
    double t = ee > 0.0 ? fmin(fmax(dot3(sub3(p, a), ba) / ee, 0.0), 1.0) : 0.0;
    return v3len(sub3(p, add3(a, scale3(ba, t)))) - r;
}
