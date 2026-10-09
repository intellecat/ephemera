#include "internal.h"

#include "data.h"

namespace eph {

namespace {

/* Evaluate sum_k c[k] * T_k(u) by Clenshaw's recurrence. */
double chebyshev(const double* c, int n, double u) {
    double b1 = 0.0, b2 = 0.0;
    for (int k = n - 1; k > 0; --k) {
        const double b0 = 2.0 * u * b1 - b2 + c[k];
        b2 = b1;
        b1 = b0;
    }
    return u * b1 - b2 + c[0];
}

void apply(const corr::Series& s, double jd_tt, double& lon, double& lat,
           double& r) {
    const double T = (jd_tt - J2000) / DAYS_PER_CENTURY;
    double u = (2.0 * T - (s.tmin + s.tmax)) / (s.tmax - s.tmin);
    if (u < -1.0) u = -1.0;
    if (u > 1.0) u = 1.0;
    lon += chebyshev(s.lon, s.n, u);
    lat += chebyshev(s.lat, s.n, u);
    r += chebyshev(s.rad, s.n, u);
}

}  // namespace

/*
 * Small DE440s corrections for the planets and the Earth, whose VSOP87
 * (DE200) series drift from DE440 (up to ~1-2 arcsec for Uranus/Neptune)
 * over 1850-2140.  The Chebyshev argument is
 * clamped to the fitted window, so outside it the correction freezes at
 * the boundary value rather than diverging.
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
