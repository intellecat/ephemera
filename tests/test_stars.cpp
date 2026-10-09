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

    /* Star-based ayanamsas vs Swiss Ephemeris (1600, 1900, 2000, 2050,
     * 2500): the stars within 0.02", the Galactic Center within 0.3". */
    const double dates[5] = {2305448.0, 2415020.5, 2451545.0, 2469807.5, 2634000.0};
    const double ref[7][5] = {
        {18.262795595, 22.444754617, 23.840014948, 24.537915331, 30.821438879},
        {14.449309485, 18.645296224, 20.045164650, 20.745308082, 27.045732220},
        {17.139975634, 21.329379683, 22.727091551, 23.426200381, 29.708319887},
        {18.998723986, 23.183675868, 24.579972884, 25.278381288, 31.575776475},
        {21.263884081, 25.449546030, 26.846045806, 27.544548218, 33.842483732},
        {16.886943064, 21.072605013, 22.469104789, 23.167607201, 29.465542715},
        {-8.736115919, -4.550453970, -3.153954194, -2.455451782, 3.842483732},
    };
    for (int m = 0; m < 7; ++m)
        for (int k = 0; k < 5; ++k) {
            const eph_ayanamsa_system s =
                (eph_ayanamsa_system)(EPH_AYANAMSA_TRUE_CITRA + m);
            const double d = std::fabs(eph_ayanamsa(s, dates[k]) * R2D - ref[m][k]) * 3600.0;
            if (d > (m < 4 ? 0.02 : 0.3)) {
                std::printf("FAIL: %s off by %.3f\"\n", eph_ayanamsa_name(s), d);
                g_fails++;
            }
        }
    eph_unload_data();
    check(std::isnan(eph_ayanamsa(EPH_AYANAMSA_TRUE_CITRA, 2451545.0)),
          "star-based ayanamsa without data is NaN");

    if (g_fails == 0) std::printf("All star checks passed.\n");
    return g_fails == 0 ? 0 : 1;
}
