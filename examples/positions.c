/* positions.c -- print apparent geocentric positions of all bodies.
 *
 * Build:  cmake -S . -B build && cmake --build build
 * Run:    ./build/ephcliapp 2026 1 1 12
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ephemeris/ephemeris.h"

static const double R2D = 180.0 / 3.14159265358979323846;

static const char* kSigns[12] = {
    "Aries", "Taurus", "Gemini", "Cancer", "Leo", "Virgo",
    "Libra", "Scorpio", "Sagittarius", "Capricorn", "Aquarius", "Pisces"};

static void format_zodiac(double lon_rad, char* buf, size_t n) {
    double deg = fmod(lon_rad * R2D, 360.0);
    if (deg < 0.0) deg += 360.0;
    int sign = (int)(deg / 30.0);
    double within = deg - sign * 30.0;
    int d = (int)within;
    double m = (within - d) * 60.0;
    snprintf(buf, n, "%2d deg %02d' %s", d, (int)(m + 0.5), kSigns[sign]);
}

int main(int argc, char** argv) {
    int year = 2000, month = 1;
    double day = 1.5;
    if (argc >= 4) {
        year = atoi(argv[1]);
        month = atoi(argv[2]);
        day = atof(argv[3]);
    }
    double jd_ut = eph_julian_date(year, month, day);
    double y = year + (month - 0.5) / 12.0;

    printf("Apparent geocentric positions for %04d-%02d-%g UT\n", year, month,
           day);
    printf("JD(UT) = %.6f   Delta-T = %.1f s\n\n", jd_ut, eph_delta_t(y));
    printf("%-8s %-24s %12s %14s %10s\n", "Body", "Ecliptic longitude",
           "Latitude", "Distance(AU)", "Speed/day");

    for (int i = 0; i < EPH_BODY_COUNT; ++i) {
        eph_body body = (eph_body)i;
        if (body == EPH_EARTH) continue;
        eph_position p;
        if (eph_position_apparent_ut(body, jd_ut, y, &p) != 0) continue;
        char z[64];
        format_zodiac(p.longitude, z, sizeof z);
        printf("%-8s %-24s %11.5f deg %14.9f %8.4f deg%s\n",
               eph_body_name(body), z, p.latitude * R2D, p.distance,
               p.longitude_speed * R2D,
               p.longitude_speed < 0 ? "  (R)" : "");
    }
    return 0;
}
