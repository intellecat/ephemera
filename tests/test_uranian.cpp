/* The Uranian points: the heliocentric mean longitudes reproduce Schlyter's
 * table of Witte's ephemeris at 1900.0 and 1950.0, the geocentric positions are the heliocentric ones seen from the Earth, and
 * they work with the searches like the bodies. */
#include <cmath>
#include <cstdio>
#include <string_view>

#include "ephemeris/ephemeris.h"
#include "load_data.h"

namespace {

const double R2D = 180.0 / M_PI;
int fails = 0;

void check(bool ok, const char* what) {
    if (!ok) {
        std::printf("FAIL: %s\n", what);
        fails++;
    }
}

}  // namespace

int main(int argc, char** argv) {
    if (!load_data_args(argc, argv)) return 1;
    /* Mean longitude (1900.0, 1950.0), e and perihelion from the table. */
    static const struct {
        const char* name;
        double lon1900, lon1950, e, peri;
    } kTable[8] = {{"Cupido", 105.026, 174.294, 0.0045, 300},
                   {"Hades", 337.180, 27.790, 0.0024, 305},
                   {"Zeus", 104.172, 144.375, 0.0011, 295},
                   {"Kronos", 17.298, 52.495, 0.0028, 205},
                   {"Apollon", 138.060, 169.297, 0, 0},
                   {"Admetos", 351.334, 20.521, 0, 0},
                   {"Vulkanus", 55.900, 83.104, 0, 0},
                   {"Poseidon", 165.521, 189.739, 0, 0}};
    for (int k = 0; k < 8; ++k) {
        const eph_body b = (eph_body)(EPH_CUPIDO + k);
        check(std::string_view(eph_body_name(b)) == kTable[k].name, "name");
        const double epochs[2] = {2415020.0, 2415020.0 + 18262.5};
        const double table[2] = {kTable[k].lon1900, kTable[k].lon1950};
        for (int t = 0; t < 2; ++t) {
            double lon, lat, r;
            if (eph_heliocentric_tt(b, epochs[t], &lon, &lat, &r) != 0) {
                check(false, "heliocentric");
                continue;
            }
            /* Remove the equation of centre to get the mean longitude. */
            const double M = (table[t] - kTable[k].peri) / R2D;
            const double mean = lon * R2D - 2 * kTable[k].e * std::sin(M) * R2D;
            if (std::fabs(std::remainder(mean - table[t], 360.0)) > 0.0015) {
                std::printf("%s at %.1f: %.4f vs %.3f\n", kTable[k].name, epochs[t], mean, table[t]);
                check(false, "mean longitude vs the table");
            }
        }
    }
    /* Geocentric = heliocentric - Earth (geometric). */
    for (int k = 0; k < 8; ++k) {
        const eph_body b = (eph_body)(EPH_CUPIDO + k);
        const double jd = 2461323.5;
        eph_position g, sun;
        double lon, lat, r;
        if (eph_position_tt(b, jd, &g) != 0 || eph_position_tt(EPH_SUN, jd, &sun) != 0 ||
            eph_heliocentric_tt(b, jd, &lon, &lat, &r) != 0) {
            check(false, "geocentric position");
            continue;
        }
        const double hx = r * std::cos(lat) * std::cos(lon), hy = r * std::cos(lat) * std::sin(lon),
                     hz = r * std::sin(lat);
        const double sx = sun.distance * std::cos(sun.latitude) * std::cos(sun.longitude),
                     sy = sun.distance * std::cos(sun.latitude) * std::sin(sun.longitude),
                     sz = sun.distance * std::sin(sun.latitude);
        const double gx = hx + sx, gy = hy + sy, gz = hz + sz;
        const double dl = std::remainder(std::atan2(gy, gx) - g.longitude, 2 * M_PI) * R2D * 3600;
        const double db = (std::asin(gz / std::sqrt(gx * gx + gy * gy + gz * gz)) - g.latitude) *
                          R2D * 3600;
        check(std::fabs(dl) < 0.01 && std::fabs(db) < 0.01, "geocentric = heliocentric - Earth");
        /* Apparent = geometric + aberration (< 21") + nutation (< 18"). */
        eph_position a;
        eph_position_apparent_tt(b, jd, &a);
        check(std::fabs(std::remainder(a.longitude - g.longitude, 2 * M_PI)) * R2D * 3600 < 40,
              "apparent vs geometric");
    }
    /* Searches: Cupido stations twice a year, from the Earth's motion. */
    double t1, t2;
    const int s1 = eph_station(EPH_CUPIDO, 2461323.5, 1, 400, &t1);
    const int s2 = eph_station(EPH_CUPIDO, t1 + 1, 1, 400, &t2);
    check(s1 > 0 && s2 > 0 && s1 != s2 && t2 - t1 > 100 && t2 - t1 < 260, "Cupido stations");
    double tc;
    check(eph_longitude_crossing(EPH_KRONOS, 0.0, 2461323.5, 1, 365.25 * 400, &tc) == 0,
          "Kronos reaches 0 Aries");

    if (fails == 0) std::printf("All Uranian checks passed.\n");
    return fails == 0 ? 0 : 1;
}
