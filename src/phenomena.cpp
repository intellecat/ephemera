#include "internal.h"

#include "../include/ephemeris/ephemeris.h"

/* Phase, elongation, apparent diameter and visual magnitude. */

namespace {

eph::Vec3 vec(const eph_position& p) {
    return eph::spherical_to_cartesian(p.longitude, p.latitude, p.distance);
}

double dot(const eph::Vec3& a, const eph::Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

double angle(const eph::Vec3& a, const eph::Vec3& b) {
    const double c = dot(a, b) / (eph::length(a) * eph::length(b));
    return std::acos(std::fmax(-1.0, std::fmin(1.0, c)));
}

/* Mean radii (km): IAU 2015 (Sun: nominal), Moon IAU. */
double radius_km(eph_body body) {
    switch (body) {
        case EPH_SUN: return 695700.0;
        case EPH_MOON: return 1737.4;
        case EPH_MERCURY: return 2439.4;
        case EPH_VENUS: return 6051.8;
        case EPH_MARS: return 3389.5;
        case EPH_JUPITER: return 69911.0;
        case EPH_SATURN: return 58232.0;
        case EPH_URANUS: return 25362.0;
        case EPH_NEPTUNE: return 24622.0;
        case EPH_PLUTO: return 1188.3;
        case EPH_CERES: return 469.7;
        case EPH_PALLAS: return 256.0;
        case EPH_JUNO: return 123.3;
        case EPH_VESTA: return 262.7;
        case EPH_CHIRON: return 103.0;
        default: return 0.0;
    }
}

/* H, G magnitude system (IAU 1985) for the minor bodies. */
double hg_magnitude(double h, double g, double r, double delta, double i) {
    const double t = std::tan(i / 2);
    const double phi1 = std::exp(-3.33 * std::pow(t, 0.63));
    const double phi2 = std::exp(-1.87 * std::pow(t, 1.22));
    return h + 5 * std::log10(r * delta) -
           2.5 * std::log10((1 - g) * phi1 + g * phi2);
}

/* Visual magnitude: Mallama & Hilton (2018), as in the Astronomical
 * Almanac, for the planets (i: phase angle in degrees; r, delta in AU;
 * ring: Saturn's ring tilt seen from the Earth, radians); the Moon after
 * Allen; H, G values from the JPL small-body database. */
double magnitude(eph_body body, double r, double delta, double i, double ring,
                 double year) {
    const double d = 5 * std::log10(r * delta);
    switch (body) {
        case EPH_SUN: return -26.74 + 5 * std::log10(delta);
        case EPH_MOON: return 0.23 + d + 0.026 * i + 4e-9 * std::pow(i, 4);
        case EPH_MERCURY:
            return -0.613 + d +
                   i * (6.3280e-2 +
                        i * (-1.6336e-3 +
                             i * (3.3644e-5 +
                                  i * (-3.4265e-7 + i * (1.6893e-9 - i * 3.0334e-12)))));
        case EPH_VENUS:
            return i < 163.7
                       ? -4.384 + d + i * (-1.044e-3 + i * (3.687e-4 +
                                           i * (-2.814e-6 + i * 8.938e-9)))
                       : -4.384 + 240.44228 + d + i * (-2.81914 + i * 8.39034e-3);
        case EPH_MARS:
            return i <= 50.0 ? -1.601 + d + i * (2.267e-2 - i * 1.302e-4)
                             : -0.367 + d + i * (-0.02573 + i * 3.445e-4);
        case EPH_JUPITER: {
            if (i <= 12.0) return -9.395 + d + i * (-3.7e-4 + i * 6.16e-4);
            const double x = i / 180.0;
            return -9.428 + d -
                   2.5 * std::log10(1 - 1.507 * x - 0.363 * x * x -
                                    0.062 * x * x * x + 2.809 * std::pow(x, 4) -
                                    1.876 * std::pow(x, 5));
        }
        case EPH_SATURN: {
            const double sb = std::fabs(std::sin(ring));
            if (i <= 6.5)
                return -8.914 + d - 1.825 * sb + 0.026 * i -
                       0.378 * sb * std::exp(-2.25 * i);
            return -8.94 + d +
                   i * (2.446e-4 + i * (2.672e-4 + i * (-1.505e-6 + i * 4.767e-9)));
        }
        case EPH_URANUS: /* ring: |planetographic sub-Earth latitude| */
            return -7.110 + d - 8.4e-4 * std::fabs(ring) / eph::DEG +
                   i * (6.587e-3 + i * 1.045e-4);
        case EPH_NEPTUNE: {
            /* Neptune brightened between 1980 and 2000. */
            const double v10 = year <= 1980.0   ? -6.89
                               : year >= 2000.0 ? -7.00
                                                : -6.89 - 0.0055 * (year - 1980.0);
            return v10 + d + i * (7.944e-3 + i * 9.617e-5);
        }
        case EPH_PLUTO: return -1.01 + d + 0.041 * i;
        case EPH_CERES: return hg_magnitude(3.33, 0.12, r, delta, i * eph::DEG);
        case EPH_PALLAS: return hg_magnitude(4.11, 0.11, r, delta, i * eph::DEG);
        case EPH_JUNO: return hg_magnitude(5.19, 0.32, r, delta, i * eph::DEG);
        case EPH_VESTA: return hg_magnitude(3.25, 0.32, r, delta, i * eph::DEG);
        case EPH_CHIRON: return hg_magnitude(5.8, 0.15, r, delta, i * eph::DEG);
        default: return std::nan("");
    }
}

/* Latitude of the Earth over a planet's equator (Saturn: its ring plane)
 * as seen from the planet: the pole is given in the ICRF (IAU) and
 * carried to the ecliptic of date by the general precession in
 * longitude. */
double sub_earth_latitude(double jd_tt, const eph::Vec3& planet, double ra_deg,
                          double dec_deg) {
    const double ra = ra_deg * eph::DEG, dec = dec_deg * eph::DEG;
    const double e0 = 23.4392911 * eph::DEG;
    const double x = std::cos(dec) * std::cos(ra), y = std::cos(dec) * std::sin(ra),
                 z = std::sin(dec);
    double lon, lat, r;
    eph::cartesian_to_spherical({x, y * std::cos(e0) + z * std::sin(e0),
                                 -y * std::sin(e0) + z * std::cos(e0)},
                                lon, lat, r);
    const double T = (jd_tt - eph::J2000) / eph::DAYS_PER_CENTURY;
    lon += 5028.796195 * T * eph::ARCSEC;
    const eph::Vec3 pole = eph::spherical_to_cartesian(lon, lat, 1.0);
    return std::asin(-dot(pole, planet) / eph::length(planet));
}

}  // namespace

extern "C" {

int eph_phenomena_tt(eph_body body, double jd_tt, eph_phenomena* out) {
    if (!out || body < 0 || body >= EPH_MEAN_NODE || body == EPH_EARTH) return -1;
    /* Apparent positions (as Swiss Ephemeris; geometric ones change the
     * phase angle by ~20-150", the magnitude by <0.001). */
    eph_position b, s;
    if (eph_position_apparent_tt(body, jd_tt, &b) != 0 ||
        eph_position_apparent_tt(EPH_SUN, jd_tt, &s) != 0)
        return -1;
    const eph::Vec3 vb = vec(b), vs = vec(s);
    const double delta = b.distance;
    double i = 0.0, r = s.distance;
    if (body != EPH_SUN) {
        const eph::Vec3 to_sun = vs - vb, to_earth = {-vb.x, -vb.y, -vb.z};
        i = angle(to_sun, to_earth);
        r = eph::length(to_sun);
    }
    out->phase_angle = i;
    out->phase = (1.0 + std::cos(i)) / 2.0;
    out->elongation = body == EPH_SUN ? 0.0 : angle(vb, vs);
    out->diameter = 2.0 * std::asin(std::fmin(1.0, radius_km(body) / (delta * eph::AU_KM)));
    const double tilt = body == EPH_SATURN   ? sub_earth_latitude(jd_tt, vb, 40.589, 83.537)
                        : body == EPH_URANUS ? sub_earth_latitude(jd_tt, vb, 257.311, -15.175)
                                             : 0.0;
    const double year = 2000.0 + (jd_tt - eph::J2000) / 365.25;
    out->magnitude = magnitude(body, r, delta, i / eph::DEG, tilt, year);
    return 0;
}

int eph_phenomena_ut(eph_body body, double jd_ut, double year, eph_phenomena* out) {
    return eph_phenomena_tt(body, eph_ut_to_tt(jd_ut, year), out);
}

}  // extern "C"
