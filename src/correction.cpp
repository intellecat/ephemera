#include "internal.h"

#include "data.h"

namespace eph {

namespace {

/* Evaluate sum_k c[k] * T_k(u) by Clenshaw's recurrence. */
double chebyshev(const float* c, int n, double u) {
    double b1 = 0.0, b2 = 0.0;
    for (int k = n - 1; k > 0; --k) {
        const double b0 = 2.0 * u * b1 - b2 + c[k];
        b2 = b1;
        b1 = b0;
    }
    return u * b1 - b2 + c[0];
}

/* Secular part plus amplitude-modulated harmonics of the body's
 * fundamental arguments (see tools/fit_correction.py).  The amplitudes
 * are clamped to the fitted window; the arguments keep running. */
void apply(const corr::Series& s, double jd_tt, double& lon, double& lat,
           double& r) {
    const double T = (jd_tt - J2000) / DAYS_PER_CENTURY;
    double u = (2.0 * T - (s.tmin + s.tmax)) / (s.tmax - s.tmin);
    if (u < -1.0) u = -1.0;
    if (u > 1.0) u = 1.0;
    double out[3];
    const int stride = s.n0 + 2 * s.nharm * s.nh;
    for (int c = 0; c < 3; ++c) out[c] = chebyshev(s.coef + c * stride, s.n0, u);
    for (int k = 0; k < s.nharm; ++k) {
        const double th = s.harm[k].phase + s.harm[k].rate * T;
        const double ct = std::cos(th), st = std::sin(th);
        for (int c = 0; c < 3; ++c) {
            const float* a = s.coef + c * stride + s.n0 + 2 * k * s.nh;
            out[c] += chebyshev(a, s.nh, u) * ct + chebyshev(a + s.nh, s.nh, u) * st;
        }
    }
    lon += out[0];
    lat += out[1];
    r += out[2];
}

}  // namespace

/*
 * DE440 corrections for the planets and the Earth, whose VSOP87 (DE200)
 * series drift from DE440 (up to ~10" for Uranus/Neptune over 1550-2650),
 * fitted by tools/fit_correction.py.  The amplitudes are clamped to the
 * fitted window, so outside it the correction freezes rather than
 * diverging.
 */
void apply_correction(int planet, double jd_tt, double& lon, double& lat,
                      double& r) {
    if (planet >= 0 && planet < 8 && data().corr_ok[planet])
        apply(data().corr[planet], jd_tt, lon, lat, r);
}

/* DE440 correction for the ELP2000-82B (DE200) lunar drift. */
void apply_moon_correction(double jd_tt, double& lon, double& lat, double& r) {
    if (data().corr_ok[10]) apply(data().corr[10], jd_tt, lon, lat, r);
}

}  // namespace eph
