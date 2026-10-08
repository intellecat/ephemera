#ifndef EPH_INTERNAL_H
#define EPH_INTERNAL_H

#include <cmath>

namespace eph {

constexpr double PI = 3.14159265358979323846;
constexpr double TWO_PI = 2.0 * PI;
constexpr double DEG = PI / 180.0;
constexpr double ARCSEC = DEG / 3600.0;

constexpr double J2000 = 2451545.0;
constexpr double DAYS_PER_CENTURY = 36525.0;
constexpr double DAYS_PER_MILLENNIUM = 365250.0;

constexpr double AU_KM = 149597870.7;
constexpr double AU_PER_LIGHT_DAY = 173.1446326846693; /* c in AU/day */

struct Vec3 {
    double x, y, z;
};

inline Vec3 operator-(const Vec3& a, const Vec3& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

/* Reduce an angle to [0, 2pi). */
inline double normalize(double a) {
    a = std::fmod(a, TWO_PI);
    if (a < 0.0) a += TWO_PI;
    return a;
}

/* VSOP87D heliocentric spherical coordinates, mean ecliptic of date. */
void vsop87_heliocentric(int planet, double jd_tt, double& lon, double& lat,
                         double& r);

/* Small fitted DE440 corrections (Uranus, Neptune); src/correction.cpp. */
void apply_correction(int planet, double jd_tt, double& lon, double& lat,
                      double& r);

/* Frame helpers (src/frames.cpp). */
Vec3 spherical_to_cartesian(double lon, double lat, double r);
void cartesian_to_spherical(const Vec3& v, double& lon, double& lat, double& r);
void precess_ecliptic_j2000_to_date(double jd_tt, double& lon, double& lat);
double mean_obliquity(double jd_tt);
void nutation(double jd_tt, double& dpsi, double& deps);
void ecliptic_to_equatorial(double lon, double lat, double jd_tt, double* ra,
                            double* dec);

/* ELP2000-82B geocentric rectangular, J2000 ecliptic, km (src/elp.cpp). */
void elp_moon_j2000(double jd_tt, double& x, double& y, double& z);

/* Lunar position, geocentric mean ecliptic of date (src/moon.cpp). */
void moon_position(double jd_tt, double& lon, double& lat, double& dist_au);

/* Pluto heliocentric J2000 ecliptic (src/pluto.cpp). */
void pluto_heliocentric_j2000(double jd_tt, double& lon, double& lat, double& r);

}  // namespace eph

#endif
