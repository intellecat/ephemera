#include "internal.h"

#include "data.h"

namespace eph {

namespace {

double chebyshev(const float* c, int n, double u) {
    double b1 = 0.0, b2 = 0.0;
    for (int k = n - 1; k > 0; --k) {
        const double b0 = 2.0 * u * b1 - b2 + c[k];
        b2 = b1;
        b1 = b0;
    }
    return u * b1 - b2 + c[0];
}

}  // namespace

bool asteroid_heliocentric(int idx, double jd_tt, Vec3& v) {
    if (idx < 0 || idx > 3 || !data().asteroid_ok[idx]) return false;
    const ast::Body& b = data().asteroid[idx];
    if (jd_tt < b.t[0] || jd_tt > b.t[b.nseg]) return false;

    int lo = 0, hi = b.nseg - 1;
    while (lo < hi) {
        const int mid = (lo + hi + 1) / 2;
        if (jd_tt < b.t[mid])
            hi = mid - 1;
        else
            lo = mid;
    }
    const double t0 = b.t[lo], t1 = b.t[lo + 1];
    const double u = (2.0 * jd_tt - (t0 + t1)) / (t1 - t0);

    double out[3];
    for (int c = 0; c < 3; ++c) {
        const float* q = b.c + (lo * 3 + c) * b.ncoef;
        out[c] = chebyshev(q, b.ncoef, u);
    }
    v = {out[0], out[1], out[2]};
    return true;
}

}  // namespace eph
