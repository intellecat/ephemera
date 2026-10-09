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

inline double length(const Vec3& v) {
    return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
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

/* VSOP87D without the terms whose amplitude within tmax_days of J2000 is
 * below threshold (radians; scaled by the mean distance for R), and
 * without the DE440 correction: the base of the residual table. */
void vsop87_truncate(int planet, double threshold, double tmax_days);
void vsop87_truncated(int planet, double jd_tt, double& lon, double& lat,
                      double& r);

/* Small fitted DE440 corrections (Mercury..Neptune, Earth); src/correction.cpp. */
void apply_correction(int planet, double jd_tt, double& lon, double& lat,
                      double& r);

/* DE440 correction for the Moon (geocentric of date, AU); src/correction.cpp. */
void apply_moon_correction(double jd_tt, double& lon, double& lat, double& r);

/* Optional residual table (src/compressed.cpp): DE440 minus the base
 * series, IAU-1976 ecliptic of date (AU), heliocentric for the planets and
 * the Earth, geocentric for the Moon, with its velocity and acceleration
 * (per day).  `base` tells whether the truncated series must be
 * added.  compressed_covers: the advertised date range.
 * compressed_generation changes on every load/unload. */
bool compressed_residual(int body, double jd, Vec3& p, Vec3& v, Vec3& a,
                         bool& base);
bool compressed_covers(double jd);
int compressed_generation();

/* Frame helpers (src/frames.cpp). */
Vec3 spherical_to_cartesian(double lon, double lat, double r);
void cartesian_to_spherical(const Vec3& v, double& lon, double& lat, double& r);
void precess_ecliptic_j2000_to_date(double jd_tt, double& lon, double& lat);
double obliquity_iau1980(double jd_tt); /* native frame of the series */
double mean_obliquity(double jd_tt);    /* IAU 2006, output frame       */
void nutation(double jd_tt, double& dpsi, double& deps);
void ecliptic_to_equatorial(double lon, double lat, double jd_tt, double* ra,
                            double* dec);
void equatorial_to_ecliptic(double ra, double dec, double jd_tt, double* lon,
                            double* lat);

/* Greenwich mean / apparent sidereal time in degrees, from UT. */
double gmst_degrees(double jd_ut);
double gast_degrees(double jd_ut);

/* Rotate a vector from the series' native IAU-1976 mean ecliptic of date
 * to the IAU-2006 mean ecliptic of date used for all output. */
Vec3 frame_1976_to_2006(double jd_tt, const Vec3& v);

/* ELP2000-82B geocentric rectangular, J2000 ecliptic, km (src/elp.cpp). */
void elp_moon_j2000(double jd_tt, double& x, double& y, double& z);

/* The same truncation for ELP2000-82B (src/elp.cpp). */
void elp_truncate(double threshold, double tmax_days);
void elp_truncated_j2000(double jd_tt, double& x, double& y, double& z);

/* Lunar position, geocentric mean ecliptic of date (src/moon.cpp). */
void moon_position(double jd_tt, double& lon, double& lat, double& dist_au);

/* Pluto heliocentric J2000 ecliptic (src/pluto.cpp). */
void pluto_heliocentric_j2000(double jd_tt, double& lon, double& lat, double& r);

/* Asteroids: heliocentric J2000 ecliptic rectangular (AU), idx 0..3
 * (Ceres, Pallas, Juno, Vesta).  src/asteroids.cpp. */
bool asteroid_heliocentric(int idx, double jd_tt, Vec3& v);

}  // namespace eph

#endif
