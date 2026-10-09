/* Lunar nodes and apogee (Lilith) against tests/lunar_points_ref.h.
 * Usage: test_lunar_points <series data file> */
#include <cmath>
#include <cstdio>

#include "ephemeris/ephemeris.h"
#include "load_data.h"
#include "lunar_points_ref.h"

static const double R2D = 180.0 / 3.14159265358979323846;

static double wrap180(double d) {
    d = std::fmod(d + 180.0, 360.0);
    return (d < 0.0 ? d + 360.0 : d) - 180.0;
}

int main(int argc, char** argv) {
    int fails = 0;
    eph_position p;

    /* The mean points need no data file. */
    if (eph_position_apparent_tt(EPH_MEAN_NODE, 2451545.0, &p) != 0 ||
        eph_position_apparent_tt(EPH_MEAN_LILITH, 2451545.0, &p) != 0) {
        std::printf("FAIL: mean points without data\n");
        fails++;
    }
    if (eph_position_apparent_tt(EPH_TRUE_NODE, 2451545.0, &p) == 0) {
        std::printf("FAIL: true node without the Moon's data\n");
        fails++;
    }
    if (!load_data_args(argc, argv)) return 1;

    /* Limits (arcsec): the mean points against another implementation of
     * the same polynomials; the osculating ones against DE440, where the
     * node magnifies the Moon's latitude error ~11x and the apogee its
     * velocity error ~20x. */
    struct {
        eph_body body;
        const char* name;
        double limit;
    } points[] = {{EPH_MEAN_NODE, "mean node", 0.5},
                  {EPH_TRUE_NODE, "true node", 1.5},
                  {EPH_MEAN_LILITH, "mean Lilith", 2.0},
                  {EPH_OSCU_LILITH, "osc. Lilith", 5.0}};
    for (const auto& pt : points) {
        double worst = 0.0;
        for (const LunarRef& r : kLunarRef) {
            if (eph_position_apparent_tt(pt.body, r.jd_tt, &p) != 0) {
                std::printf("FAIL: %s at %.1f\n", pt.name, r.jd_tt);
                return 1;
            }
            double lon = r.mean_node, lat = 0.0;
            if (pt.body == EPH_TRUE_NODE) lon = r.true_node;
            if (pt.body == EPH_MEAN_LILITH) lon = r.mean_apogee_lon, lat = r.mean_apogee_lat;
            if (pt.body == EPH_OSCU_LILITH) lon = r.oscu_apogee_lon, lat = r.oscu_apogee_lat;
            const double e = std::fmax(
                std::fabs(wrap180(p.longitude * R2D - lon)) *
                    std::cos(lat / R2D),
                std::fabs(p.latitude * R2D - lat)) * 3600.0;
            worst = std::fmax(worst, e);
        }
        std::printf("%-12s max %.3f\"\n", pt.name, worst);
        if (worst > pt.limit) {
            std::printf("FAIL: %s exceeds %.1f\"\n", pt.name, pt.limit);
            fails++;
        }
    }

    /* Speeds against central differences of the positions (+-0.01 day:
     * the osculating apogee swings fast, and a double Julian date resolves
     * ~5e-10 day, which jitters it by ~0.005"). */
    for (const auto& pt : points) {
        double worst = 0.0;
        for (const LunarRef& r : kLunarRef) {
            eph_position a, b;
            eph_position_apparent_tt(pt.body, r.jd_tt, &p);
            eph_position_apparent_tt(pt.body, r.jd_tt - 0.01, &a);
            eph_position_apparent_tt(pt.body, r.jd_tt + 0.01, &b);
            const double fd = wrap180((b.longitude - a.longitude) * R2D) / 0.02;
            const double err = std::fabs(p.longitude_speed * R2D - fd) * 3600.0;
            worst = std::fmax(worst, err / (1.0 + 0.002 * std::fabs(fd) * 3600.0));
        }
        if (worst > 1.0) {
            std::printf("FAIL: %s speed off by %.3f (1\"/day + 0.2%%)\n", pt.name,
                        worst);
            fails++;
        }
    }

    if (fails == 0) std::printf("\nAll lunar point checks passed.\n");
    return fails == 0 ? 0 : 1;
}
