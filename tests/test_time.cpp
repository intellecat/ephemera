/* Time helpers and utilities.  Usage: test_time <series data file> */
#include <cmath>
#include <cstdio>

#include "ephemeris/ephemeris.h"
#include "load_data.h"

static const double D2R = 3.14159265358979323846 / 180.0;
static int g_fails = 0;

static void check(bool ok, const char* what) {
    if (!ok) {
        std::printf("FAIL: %s\n", what);
        g_fails++;
    }
}

int main(int argc, char** argv) {
    if (!load_data_args(argc, argv)) return 1;
    double tt, ut, tt2;

    /* UTC: TT = UTC + 37 s + 32.184 s since 2017; a leap second at the
     * end of 2016 makes 23:59:60 and 00:00:00 one second apart. */
    check(eph_utc_to_jd(2020, 6, 1, 12, 0, 0.0, &tt, &ut) == 0 &&
              std::fabs((tt - eph_julian_date(2020, 6, 1.5)) * 86400.0 - 69.184) < 1e-4,
          "TT - UTC in 2020");
    eph_utc_to_jd(2016, 12, 31, 23, 59, 60.0, &tt, nullptr);
    eph_utc_to_jd(2017, 1, 1, 0, 0, 0.0, &tt2, nullptr);
    check(std::fabs((tt2 - tt) * 86400.0 - 1.0) < 1e-4, "leap second 2016-12-31");
    check(eph_utc_to_jd(2021, 4, 31, 0, 0, 0.0, &tt, &ut) == -1, "31 April rejected");
    check(eph_utc_to_jd(2021, 13, 1, 0, 0, 0.0, &tt, &ut) == -1, "month 13 rejected");

    /* Equation of time vs Swiss Ephemeris (days). */
    const struct {
        double jd, e;
    } eot[] = {
        {2448908.50, 0.009518443475},
        {2451545.00, -0.002281427250},
        {2460400.25, -0.002929058894},
    };
    for (const auto& r : eot) {
        double e;
        check(eph_equation_of_time(r.jd, &e) == 0 && std::fabs(e - r.e) * 86400.0 < 0.02,
              "equation of time");
    }

    /* Sidereal time at J2000 UT (ERFA gst06a: 280.457072 deg). */
    check(std::fabs(eph_sidereal_time(2451545.0) / D2R - 280.457072) < 1e-6,
          "sidereal time at J2000");

    /* Utilities. */
    int sign, deg, min;
    double sec;
    eph_split_longitude((359.0 + 59.0 / 60 + 59.99999 / 3600) * D2R, &sign, &deg,
                        &min, &sec);
    check(sign == 11 && deg == 29 && min == 59 && std::fabs(sec - 59.99999) < 1e-5,
          "split longitude");
    eph_split_longitude((359.0 + 59.0 / 60 + 59.9999999 / 3600) * D2R, &sign, &deg,
                        &min, &sec);
    check(sign == 0 && deg == 0 && min == 0 && sec == 0.0,
          "split rounds 59.9999999\" up to 0 Aries");
    eph_split_longitude(45.5 * D2R, &sign, &deg, &min, &sec);
    check(sign == 1 && deg == 15 && min == 30 && std::fabs(sec) < 1e-6,
          "split 45.5 deg");
    check(std::fabs(eph_midpoint(350 * D2R, 20 * D2R) / D2R - 5.0) < 1e-9 &&
              std::fabs(eph_midpoint(20 * D2R, 350 * D2R) / D2R - 5.0) < 1e-9,
          "midpoint across 0 Aries");

    if (g_fails == 0) std::printf("All time checks passed.\n");
    return g_fails == 0 ? 0 : 1;
}
