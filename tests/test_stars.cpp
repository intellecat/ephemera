/* Fixed stars and the star-based ayanamsas.
 * Usage: test_stars <series data file> <star data file> */
#include <cmath>
#include <cstdio>
#include <cstring>

#include "ephemeris/ephemeris.h"
#include "load_data.h"
#include "stars_ref.h"

static const double R2D = 180.0 / 3.14159265358979323846;
static int g_fails = 0;

static void check(bool ok, const char* what) {
    if (!ok) {
        std::printf("FAIL: %s\n", what);
        g_fails++;
    }
}

int main(int argc, char** argv) {
    check(eph_star_count() == 0 && eph_star_find("Regulus") == -1,
          "no stars before loading");
    if (!load_data_args(argc, argv)) return 1;

    /* Lookup. */
    const int regulus = eph_star_find("Regulus");
    check(regulus >= 0 && eph_star_find("regulus") == regulus &&
              eph_star_find("alf Leo") == regulus &&
              eph_star_find("Alpha Leo") == regulus &&
              eph_star_find("HIP 49669") == regulus &&
              std::strcmp(eph_star_designation(regulus), "alf Leo") == 0,
          "lookup by name, designation, HIP");
    check(eph_star_find("Asellus Australis") == eph_star_find("del Cnc"),
          "two-word names");
    check(eph_star_find("No Such Star") == -1 && eph_star_name(-1) == nullptr,
          "unknown star");
    check(std::fabs(eph_star_magnitude(regulus) - 1.36) < 0.01, "magnitude");

    /* Apparent positions vs ERFA's space motion + Skyfield/DE440. */
    double worst = 0.0;
    for (const StarRef& r : kStarRef) {
        eph_position p;
        if (eph_star_position_apparent_tt(eph_star_find(r.name), r.jd_tt, &p) != 0) {
            std::printf("FAIL: %s\n", r.name);
            g_fails++;
            continue;
        }
        const double dl = std::remainder(p.longitude * R2D - r.lon, 360.0) *
                          std::cos(r.lat / R2D);
        worst = std::fmax(worst, std::hypot(dl, p.latitude * R2D - r.lat) * 3600.0);
    }
    std::printf("stars vs ERFA + Skyfield: max %.4f\"\n", worst);
    check(worst < 0.05, "star positions");

    /* Star-based ayanamsas (1600, 1900, 2000, 2050, 2500) from their
     * definitions: the star's apparent longitude (ERFA + Skyfield/DE440)
     * minus nutation minus its defined sidereal longitude. */
    double worst_ayan = 0.0;
    for (int m = 0; m < 7; ++m)
        for (int k = 0; k < 5; ++k) {
            const eph_ayanamsa_system s =
                (eph_ayanamsa_system)(EPH_AYANAMSA_TRUE_CITRA + m);
            const double d = std::fabs(eph_ayanamsa(s, kStarAyanamsaDates[k]) * R2D -
                                       kStarAyanamsaRef[m][k]) * 3600.0;
            worst_ayan = std::fmax(worst_ayan, d);
            if (d > 0.05) {
                std::printf("FAIL: %s off by %.3f\"\n", eph_ayanamsa_name(s), d);
                g_fails++;
            }
        }
    std::printf("star-based ayanamsas: max %.4f\"\n", worst_ayan);
    eph_unload_data();
    check(std::isnan(eph_ayanamsa(EPH_AYANAMSA_TRUE_CITRA, 2451545.0)),
          "star-based ayanamsa without data is NaN");

    if (g_fails == 0) std::printf("All star checks passed.\n");
    return g_fails == 0 ? 0 : 1;
}
