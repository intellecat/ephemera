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

bool chebyshev_body(const double* t, const float* c, int nseg, int ncoef,
                    double jd_tt, double anchor, Vec3& v) {
    if (anchor < t[0] || anchor > t[nseg]) return false;

    int lo = 0, hi = nseg - 1;
    while (lo < hi) {
        const int mid = (lo + hi + 1) / 2;
        if (anchor < t[mid])
            hi = mid - 1;
        else
            lo = mid;
    }
    const double t0 = t[lo], t1 = t[lo + 1];
    const double u = (2.0 * jd_tt - (t0 + t1)) / (t1 - t0);

    double out[3];
    for (int k = 0; k < 3; ++k) out[k] = chebyshev(c + (lo * 3 + k) * ncoef, ncoef, u);
    v = {out[0], out[1], out[2]};
    return true;
}

bool asteroid_heliocentric(int idx, double jd_tt, double anchor, Vec3& v) {
    if (idx < 0 || idx > 3 || !data().asteroid_ok[idx]) return false;
    const ast::Body& b = data().asteroid[idx];
    return chebyshev_body(b.t, b.c, b.nseg, b.ncoef, jd_tt, anchor, v);
}

}  // namespace eph
