/* Event searches.  Usage: test_search <series data file> */
#include <cmath>
#include <cstdio>
#include <initializer_list>

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

static double wrap(double a) { return std::remainder(a, 2 * 3.14159265358979323846); }

int main(int argc, char** argv) {
    if (!load_data_args(argc, argv)) return 1;
    eph_position p;

    /* Solar and lunar crossings (TT) vs JPL DE440 with Skyfield
     * (tools/make_search_reference.py). */
    const struct {
        double jd, lon, sun, moon;
    } ref[] = {
    {2461756.230099, 103.254195, 2461957.231758161, 2461782.119385233},
    {2455078.703751, 279.912270, 2455196.971830348, 2455101.262528104},
    {2457901.865740, 329.536843, 2458167.762648552, 2457919.892697950},
    {2461509.841230, 330.565546, 2461821.204449217, 2461526.208840692},
    {2440664.693362, 157.409280, 2440829.737409717, 2440666.556949975},
    {2452123.608174, 23.455511, 2452378.076380586, 2452131.735376221},
    {2440140.651267, 299.023717, 2440240.776474946, 2440157.584987245},
    {2464582.556107, 282.459547, 2464695.788852571, 2464609.616260689},
    };
    double worst_sun = 0.0, worst_moon = 0.0;
    for (const auto& r : ref) {
        double t;
        check(eph_longitude_crossing(EPH_SUN, r.lon * D2R, r.jd, 1, 400, &t) == 0,
              "solar crossing found");
        worst_sun = std::fmax(worst_sun, std::fabs(t - r.sun) * 86400.0);
        check(eph_longitude_crossing(EPH_MOON, r.lon * D2R, r.jd, 1, 40, &t) == 0,
              "lunar crossing found");
        worst_moon = std::fmax(worst_moon, std::fabs(t - r.moon) * 86400.0);
    }
    std::printf("Crossings vs DE440: Sun %.3f s, Moon %.3f s\n", worst_sun, worst_moon);
    /* Measured with eph-full.bin: 0.13 s and 0.015 s (~0.005" and 0.008"). */
    check(worst_sun < 0.3 && worst_moon < 0.05, "crossing times");

    /* Every found event is one, and searching back from just after it finds
     * it again. */
    for (int k = 0; k < 20; ++k) {
        const double jd = 2440000.0 + 1234.567 * k;
        const double target = (17.0 * k) * D2R;
        for (eph_body b : {EPH_SUN, EPH_MOON, EPH_MARS, EPH_JUPITER, EPH_TRUE_NODE}) {
            double t, back;
            if (eph_longitude_crossing(b, target, jd, 1, 7000, &t) != 0) {
                check(false, "crossing found");
                continue;
            }
            eph_position_apparent_tt(b, t, &p);
            check(t > jd && std::fabs(wrap(p.longitude - target)) < 1e-9,
                  "longitude at crossing");
            check(eph_longitude_crossing(b, target, t + 0.5, -1, 1.0, &back) == 0 &&
                      std::fabs(back - t) < 1e-6,
                  "backward search");
        }
        double t, lon;
        const int node = eph_node_crossing(EPH_MOON, jd, 1, 40, &t, &lon);
        eph_position_apparent_tt(EPH_MOON, t, &p);
        check(node >= 1 && std::fabs(p.latitude) < 1e-10 &&
                  (node == 1) == (p.latitude_speed > 0),
              "node crossing");
        const int st = eph_station(EPH_MERCURY, jd, 1, 200, &t);
        eph_position_apparent_tt(EPH_MERCURY, t, &p);
        eph_position after;
        eph_position_apparent_tt(EPH_MERCURY, t + 1.0, &after);
        check(st >= 1 && std::fabs(p.longitude_speed) < 1e-9 &&
                  (st == 1) == (after.longitude_speed < 0),
              "station");
    }

    /* Heliocentric crossing: Mars' heliocentric longitude. */
    {
        double t, lon;
        check(eph_heliocentric_crossing(EPH_MARS, 1.0, 2451545.0, 1, 700, &t) == 0 &&
                  eph_heliocentric_tt(EPH_MARS, t, &lon, nullptr, nullptr) == 0 &&
                  std::fabs(wrap(lon - 1.0)) < 1e-9,
              "heliocentric crossing");
    }

    /* No event in the window, bad arguments. */
    double t;
    check(eph_longitude_crossing(EPH_SUN, 0.0, 2451545.0, 1, 10, &t) == -1,
          "no crossing in 10 days");
    check(eph_longitude_crossing(EPH_SUN, 0.0, 2451545.0, 0, 400, &t) == -1,
          "direction 0 rejected");

    /* UT variant: the TT result shifted by Delta-T. */
    {
        double ut, tt;
        eph_longitude_crossing_ut(EPH_SUN, 0.0, 2460000.5, 1, 400, &ut);
        eph_longitude_crossing(EPH_SUN, 0.0, eph_ut_to_tt(2460000.5, 2023.2), 1, 400, &tt);
        check(std::fabs((tt - ut) * 86400.0 - eph_delta_t(2024.2)) < 0.5, "UT crossing");
    }

    if (g_fails == 0) std::printf("All search checks passed.\n");
    return g_fails == 0 ? 0 : 1;
}
