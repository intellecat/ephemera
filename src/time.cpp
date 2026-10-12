#include "internal.h"

#include "../include/ephemeris/ephemeris.h"
#include "../data/delta_t_data.h"

namespace {

/* Espenak & Meeus polynomial expressions for Delta-T = TT - UT (seconds),
 * used outside the tabulated range. */
double delta_t_polynomial(double y) {
    double u, t;
    if (y < -500.0) {
        u = (y - 1820.0) / 100.0;
        return -20.0 + 32.0 * u * u;
    }
    if (y < 500.0) {
        u = y / 100.0;
        return 10583.6 + u * (-1014.41 +
               u * (33.78311 +
               u * (-5.952053 +
               u * (-0.1798452 +
               u * (0.022174192 + u * 0.0090316521)))));
    }
    if (y < 1600.0) {
        u = (y - 1000.0) / 100.0;
        return 1574.2 + u * (-556.01 +
               u * (71.23472 +
               u * (0.319781 +
               u * (-0.8503463 +
               u * (-0.005050998 + u * 0.0083572073)))));
    }
    if (y < 1700.0) {
        t = y - 1600.0;
        return 120.0 + t * (-0.9808 + t * (-0.01532 + t / 7129.0));
    }
    if (y < 1800.0) {
        t = y - 1700.0;
        return 8.83 + t * (0.1603 +
               t * (-0.0059285 +
               t * (0.00013336 - t / 1174000.0)));
    }
    if (y < 1860.0) {
        t = y - 1800.0;
        return 13.72 + t * (-0.332447 +
               t * (0.0068612 +
               t * (0.0041116 +
               t * (-0.00037436 +
               t * (0.0000121272 +
               t * (-0.0000001699 + t * 0.000000000875))))));
    }
    if (y < 1900.0) {
        t = y - 1860.0;
        return 7.62 + t * (0.5737 +
               t * (-0.251754 +
               t * (0.01680668 +
               t * (-0.0004473624 + t / 233174.0))));
    }
    if (y < 1920.0) {
        t = y - 1900.0;
        return -2.79 + t * (1.494119 +
               t * (-0.0598939 +
               t * (0.0061966 - t * 0.000197)));
    }
    if (y < 1941.0) {
        t = y - 1920.0;
        return 21.20 + t * (0.84493 + t * (-0.076100 + t * 0.0020936));
    }
    if (y < 1961.0) {
        t = y - 1950.0;
        return 29.07 + t * (0.407 + t * (-1.0 / 233.0 + t / 2547.0));
    }
    if (y < 1986.0) {
        t = y - 1975.0;
        return 45.45 + t * (1.067 + t * (-1.0 / 260.0 - t / 718.0));
    }
    if (y < 2005.0) {
        t = y - 2000.0;
        return 63.86 + t * (0.3345 +
               t * (-0.060374 +
               t * (0.0017275 +
               t * (0.000651814 + t * 0.00002373599))));
    }
    if (y < 2050.0) {
        t = y - 2000.0;
        return 62.92 + t * (0.32217 + t * 0.005589);
    }
    if (y < 2150.0) {
        u = (y - 1820.0) / 100.0;
        return -20.0 + 32.0 * u * u - 0.5628 * (2150.0 - y);
    }
    u = (y - 1820.0) / 100.0;
    return -20.0 + 32.0 * u * u;
}

/* Delta-T: yearly table (IERS observations, Stephenson et al. 2016) with
 * linear interpolation; outside it the polynomials, offset to join the
 * table continuously. */
double delta_t_seconds(double y) {
    using namespace eph::delta_t;
    if (y < kFirstYear)
        return delta_t_polynomial(y) + kTable[0] -
               delta_t_polynomial(kFirstYear);
    if (y >= kLastYear)
        return delta_t_polynomial(y) + kTable[kLastYear - kFirstYear] -
               delta_t_polynomial(kLastYear);
    const int i = static_cast<int>(std::floor(y)) - kFirstYear;
    const double f = y - std::floor(y);
    return kTable[i] + f * (kTable[i + 1] - kTable[i]);
}

}  // namespace

extern "C" {

double eph_julian_date(int year, int month, double day) {
    if (month <= 2) {
        year -= 1;
        month += 12;
    }
    const bool gregorian =
        (year > 1582) || (year == 1582 && (month > 10 || (month == 10 && day >= 15.0)));
    const long A = year / 100;
    const long B = gregorian ? 2 - A + A / 4 : 0;
    return std::floor(365.25 * (year + 4716)) +
           std::floor(30.6001 * (month + 1)) + day + B - 1524.5;
}

void eph_calendar_date(double jd, int* year, int* month, double* day) {
    const double z0 = jd + 0.5;
    double Z = std::floor(z0);
    const double F = z0 - Z;
    double A;
    if (Z < 2299161.0) {
        A = Z;
    } else {
        const double alpha = std::floor((Z - 1867216.25) / 36524.25);
        A = Z + 1.0 + alpha - std::floor(alpha / 4.0);
    }
    const double B = A + 1524.0;
    const double C = std::floor((B - 122.1) / 365.25);
    const double D = std::floor(365.25 * C);
    const double E = std::floor((B - D) / 30.6001);
    *day = B - D - std::floor(30.6001 * E) + F;
    *month = static_cast<int>(E < 14.0 ? E - 1.0 : E - 13.0);
    *year = static_cast<int>(*month > 2 ? C - 4716.0 : C - 4715.0);
}

double eph_delta_t(double year) { return delta_t_seconds(year); }

double eph_tdb_minus_tt(double jd_tt) {
    /* Fairhead & Bretagnon (1990), as in USNO Circular 179; T in Julian
     * centuries from J2000.  Dominant terms, good to about 1 microsecond. */
    const double T = (jd_tt - 2451545.0) / 36525.0;
    return 0.001657 * std::sin(628.3076 * T + 6.2401) +
           0.000022 * std::sin(575.3385 * T + 4.2970) +
           0.000014 * std::sin(1256.6152 * T + 6.1969) +
           0.000005 * std::sin(606.9777 * T + 4.0212) +
           0.000005 * std::sin(52.9691 * T + 0.4444) +
           0.000002 * std::sin(21.3299 * T + 5.5431) +
           0.000010 * T * std::sin(628.3076 * T + 4.2490);
}

double eph_ut_to_tt(double jd_ut, double year) {
    return jd_ut + delta_t_seconds(year) / 86400.0;
}

double eph_tt_to_ut(double jd_tt, double year) {
    return jd_tt - delta_t_seconds(year) / 86400.0;
}

double eph_mean_obliquity(double jd_tt) { return eph::mean_obliquity(jd_tt); }

void eph_nutation(double jd_tt, double* dpsi, double* deps) {
    eph::nutation(jd_tt, *dpsi, *deps);
}

void eph_ecliptic_to_equatorial(double lon, double lat, double jd_tt,
                                double* ra, double* dec) {
    eph::ecliptic_to_equatorial(lon, lat, jd_tt, ra, dec);
}

const char* eph_body_name(eph_body body) {
    if (body >= EPH_MINOR_BASE) {
        const char* name = eph_minor_name(body - EPH_MINOR_BASE);
        return name ? name : "?";
    }
    switch (body) {
        case EPH_MERCURY: return "Mercury";
        case EPH_VENUS:   return "Venus";
        case EPH_EARTH:   return "Earth";
        case EPH_MARS:    return "Mars";
        case EPH_JUPITER: return "Jupiter";
        case EPH_SATURN:  return "Saturn";
        case EPH_URANUS:  return "Uranus";
        case EPH_NEPTUNE: return "Neptune";
        case EPH_PLUTO:   return "Pluto";
        case EPH_SUN:     return "Sun";
        case EPH_MOON:    return "Moon";
        case EPH_CERES:   return "Ceres";
        case EPH_PALLAS:  return "Pallas";
        case EPH_JUNO:    return "Juno";
        case EPH_VESTA:   return "Vesta";
        case EPH_CHIRON:  return "Chiron";
        case EPH_MEAN_NODE:   return "Mean Node";
        case EPH_TRUE_NODE:   return "True Node";
        case EPH_MEAN_LILITH: return "Mean Lilith";
        case EPH_OSCU_LILITH: return "Osc. Lilith";
        case EPH_CUPIDO:      return "Cupido";
        case EPH_HADES:       return "Hades";
        case EPH_ZEUS:        return "Zeus";
        case EPH_KRONOS:      return "Kronos";
        case EPH_APOLLON:     return "Apollon";
        case EPH_ADMETOS:     return "Admetos";
        case EPH_VULKANUS:    return "Vulkanus";
        case EPH_POSEIDON:    return "Poseidon";
        default:          return "?";
    }
}

}  // extern "C"
