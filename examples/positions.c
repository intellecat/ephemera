/* positions.c -- print apparent geocentric positions of all bodies.
 *
 * Build:  cmake -S . -B build && cmake --build build
 * Run:    ./build/ephcliapp 2026 1 1 12      (year month day [UT hours])
 *
 * Data files come from $EPH_DATA (paths separated by ':'), or by default
 * data/eph-full.bin, data/eph-asteroids.bin and the DE440 tables
 * data/de440/de440-YEAR.bin; a build with embedded data needs none.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

/* Load the data files; returns how many were loaded. */
static int load_data(void) {
    const char* env = getenv("EPH_DATA");
    int n = 0;
    if (env && *env) {
        char path[1024];
        while (*env) {
            size_t len = strcspn(env, ":");
            if (len > 0 && len < sizeof path) {
                memcpy(path, env, len);
                path[len] = '\0';
                if (eph_load_file(path) == 0)
                    n++;
                else
                    fprintf(stderr, "cannot load %s\n", path);
            }
            env += len + (env[len] == ':');
        }
        return n;
    }
    const char* defaults[] = {"data/eph-full.bin", "data/eph-asteroids.bin"};
    for (size_t i = 0; i < sizeof defaults / sizeof defaults[0]; ++i)
        if (eph_load_file(defaults[i]) == 0) n++;
    for (int year = 1550; year < 2650; year += 100) {
        char path[64];
        snprintf(path, sizeof path, "data/de440/de440-%d.bin", year);
        if (eph_load_file(path) == 0) n++;
    }
    return n;
}

int main(int argc, char** argv) {
    /* A build with embedded data needs no file. */
    eph_position probe;
    if (load_data() == 0 && eph_position_tt(EPH_SUN, 2451545.0, &probe) != 0) {
        fprintf(stderr, "no data file loaded: set EPH_DATA or run from the "
                        "source directory\n");
        return 1;
    }
    int year = 2000, month = 1;
    double day = 1.5;
    if (argc >= 4) {
        year = atoi(argv[1]);
        month = atoi(argv[2]);
        day = atof(argv[3]);
        if (argc >= 5) day += atof(argv[4]) / 24.0;
    }
    double jd_ut = eph_julian_date(year, month, day);
    double y = year + (jd_ut - eph_julian_date(year, 1, 1.0)) / 365.25;

    int d = (int)day;
    int minutes = (int)((day - d) * 1440.0 + 0.5);
    printf("Apparent geocentric positions for %04d-%02d-%02d %02d:%02d UT\n",
           year, month, d, minutes / 60, minutes % 60);
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
