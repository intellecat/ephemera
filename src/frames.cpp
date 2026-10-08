#include "internal.h"

namespace eph {

Vec3 spherical_to_cartesian(double lon, double lat, double r) {
    const double cl = std::cos(lat);
    return {r * cl * std::cos(lon), r * cl * std::sin(lon), r * std::sin(lat)};
}

void cartesian_to_spherical(const Vec3& v, double& lon, double& lat, double& r) {
    r = std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    lon = normalize(std::atan2(v.y, v.x));
    lat = (r > 0.0) ? std::asin(v.z / r) : 0.0;
}

/*
 * Precession of ecliptic coordinates from J2000.0 to the mean equinox and
 * ecliptic of date.  Meeus, Astronomical Algorithms, eq. (21.5); the
 * accumulated precession angles refer to the origin J2000 so the epoch
 * term T is zero.
 */
void precess_ecliptic_j2000_to_date(double jd_tt, double& lon, double& lat) {
    const double t = (jd_tt - J2000) / DAYS_PER_CENTURY;

    const double eta = (47.0029 * t - 0.03302 * t * t + 0.000060 * t * t * t) *
                       ARCSEC;
    const double Pi = (174.876384 * 3600.0 - 869.8089 * t + 0.000598 * t * t) *
                      ARCSEC;
    const double p = (5029.0966 * t + 1.11113 * t * t - 0.000006 * t * t * t) *
                     ARCSEC;

    const double cb = std::cos(lat);
    const double sb = std::sin(lat);
    const double ce = std::cos(eta);
    const double se = std::sin(eta);
    const double d = Pi - lon;

    const double A = ce * cb * std::cos(d) - se * sb;
    const double B = cb * std::sin(d);
    const double C = ce * sb + se * cb * std::cos(d);

    lon = normalize(p + Pi - std::atan2(B, A));
    lat = std::asin(C);
}

double mean_obliquity(double jd_tt) {
    const double T = (jd_tt - J2000) / DAYS_PER_CENTURY;
    const double seconds =
        84381.448 + T * (-46.8150 + T * (-0.00059 + T * 0.001813));
    return seconds * ARCSEC;
}

namespace {

struct NutTerm {
    double d, m, n, f, omega;  /* multipliers */
    double s0, s1;             /* dpsi: (s0 + s1*T) * 1e-4 arcsec */
    double c0, c1;             /* deps: (c0 + c1*T) * 1e-4 arcsec */
};

/* IAU 1980 nutation series, 63 terms (terms < 0.0003" dropped). */
const NutTerm kNutation[] = {
    {0, 0, 0, 0, 1, -171996, -174.2, 92025, 8.9},
    {-2, 0, 0, 2, 2, -13187, -1.6, 5736, -3.1},
    {0, 0, 0, 2, 2, -2274, -0.2, 977, -0.5},
    {0, 0, 0, 0, 2, 2062, 0.2, -895, 0.5},
    {0, 1, 0, 0, 0, 1426, -3.4, 54, -0.1},
    {0, 0, 1, 0, 0, 712, 0.1, -7, 0},
    {-2, 1, 0, 2, 2, -517, 1.2, 224, -0.6},
    {0, 0, 0, 2, 1, -386, -0.4, 200, 0},
    {0, 0, 1, 2, 2, -301, 0, 129, -0.1},
    {-2, -1, 0, 2, 2, 217, -0.5, -95, 0.3},
    {-2, 0, 1, 0, 0, -158, 0, 0, 0},
    {-2, 0, 0, 2, 1, 129, 0.1, -70, 0},
    {0, 0, -1, 2, 2, 123, 0, -53, 0},
    {2, 0, 0, 0, 0, 63, 0, 0, 0},
    {0, 0, 1, 0, 1, 63, 0.1, -33, 0},
    {2, 0, -1, 2, 2, -59, 0, 26, 0},
    {0, 0, -1, 0, 1, -58, -0.1, 32, 0},
    {0, 0, 1, 2, 1, -51, 0, 27, 0},
    {-2, 0, 2, 0, 0, 48, 0, 0, 0},
    {0, 0, -2, 2, 1, 46, 0, -24, 0},
    {2, 0, 0, 2, 2, -38, 0, 16, 0},
    {0, 0, 2, 2, 2, -31, 0, 13, 0},
    {0, 0, 2, 0, 0, 29, 0, 0, 0},
    {-2, 0, 1, 2, 2, 29, 0, -12, 0},
    {0, 0, 0, 2, 0, 26, 0, 0, 0},
    {-2, 0, 0, 2, 0, -22, 0, 0, 0},
    {0, 0, -1, 2, 1, 21, 0, -10, 0},
    {0, 2, 0, 0, 0, 17, -0.1, 0, 0},
    {2, 0, -1, 0, 1, 16, 0, -8, 0},
    {-2, 2, 0, 2, 2, -16, 0.1, 7, 0},
    {0, 1, 0, 0, 1, -15, 0, 9, 0},
    {-2, 0, 1, 0, 1, -13, 0, 7, 0},
    {0, -1, 0, 0, 1, -12, 0, 6, 0},
    {0, 0, 2, -2, 0, 11, 0, 0, 0},
    {2, 0, -1, 2, 1, -10, 0, 5, 0},
    {2, 0, 1, 2, 2, -8, 0, 3, 0},
    {0, 1, 0, 2, 2, 7, 0, -3, 0},
    {-2, 1, 1, 0, 0, -7, 0, 0, 0},
    {0, -1, 0, 2, 2, -7, 0, 3, 0},
    {2, 0, 0, 2, 1, -7, 0, 3, 0},
    {2, 0, 1, 0, 0, 6, 0, 0, 0},
    {-2, 0, 2, 2, 2, 6, 0, -3, 0},
    {-2, 0, 1, 2, 1, 6, 0, -3, 0},
    {2, 0, -2, 0, 1, -6, 0, 3, 0},
    {2, 0, 0, 0, 1, -6, 0, 3, 0},
    {0, -1, 1, 0, 0, 5, 0, 0, 0},
    {-2, -1, 0, 2, 1, -5, 0, 3, 0},
    {-2, 0, 0, 0, 1, -5, 0, 3, 0},
    {0, 0, 2, 2, 1, -5, 0, 3, 0},
    {-2, 0, 2, 0, 1, 4, 0, 0, 0},
    {-2, 1, 0, 2, 1, 4, 0, 0, 0},
    {0, 0, 1, -2, 0, 4, 0, 0, 0},
    {-1, 0, 1, 0, 0, -4, 0, 0, 0},
    {-2, 1, 0, 0, 0, -4, 0, 0, 0},
    {1, 0, 0, 0, 0, -4, 0, 0, 0},
    {0, 0, 1, 2, 0, 3, 0, 0, 0},
    {0, 0, -2, 2, 2, -3, 0, 0, 0},
    {-1, -1, 1, 0, 0, -3, 0, 0, 0},
    {0, 1, 1, 0, 0, -3, 0, 0, 0},
    {0, -1, 1, 2, 2, -3, 0, 0, 0},
    {2, -1, -1, 2, 2, -3, 0, 0, 0},
    {0, 0, 3, 2, 2, -3, 0, 0, 0},
    {2, -1, 0, 2, 2, -3, 0, 0, 0},
};

}  // namespace

void nutation(double jd_tt, double& dpsi, double& deps) {
    const double T = (jd_tt - J2000) / DAYS_PER_CENTURY;

    const double D = (297.85036 + 445267.11148 * T - 0.0019142 * T * T +
                      T * T * T / 189474.0) * DEG;
    const double M = (357.52772 + 35999.050340 * T - 0.0001603 * T * T -
                      T * T * T / 300000.0) * DEG;
    const double N = (134.96298 + 477198.867398 * T + 0.0086972 * T * T +
                      T * T * T / 5620.0) * DEG;
    const double F = (93.27191 + 483202.017538 * T - 0.0036825 * T * T +
                      T * T * T / 327270.0) * DEG;
    const double Om = (125.04452 - 1934.136261 * T + 0.0020708 * T * T +
                       T * T * T / 450000.0) * DEG;

    double dpsi_s = 0.0, deps_s = 0.0;
    for (const NutTerm& row : kNutation) {
        const double arg = row.d * D + row.m * M + row.n * N + row.f * F +
                           row.omega * Om;
        const double s = std::sin(arg);
        const double c = std::cos(arg);
        dpsi_s += s * (row.s0 + row.s1 * T);
        deps_s += c * (row.c0 + row.c1 * T);
    }
    dpsi = dpsi_s * 1e-4 * ARCSEC;
    deps = deps_s * 1e-4 * ARCSEC;
}

void ecliptic_to_equatorial(double lon, double lat, double jd_tt, double* ra,
                            double* dec) {
    double dpsi, deps;
    nutation(jd_tt, dpsi, deps);
    const double eps = mean_obliquity(jd_tt) + deps;
    const double ce = std::cos(eps), se = std::sin(eps);
    const double sl = std::sin(lon), cl = std::cos(lon);
    const double sb = std::sin(lat), cb = std::cos(lat);
    const double x = cb * cl;
    const double y = ce * cb * sl - se * sb;
    const double z = se * cb * sl + ce * sb;
    *ra = normalize(std::atan2(y, x));
    *dec = std::asin(z);
}

}  // namespace eph
