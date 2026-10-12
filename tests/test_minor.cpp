/* Minor planets from a minor-planet file: apparent geocentric RA/Dec
 * against JPL Horizons (tests/minor_ref.h, tools/make_minor_reference.py),
 * lookup, names and phenomena.  Usage: test_minor <series file> <minor file> */
#include <cmath>
#include <cstdio>
#include <string_view>

#include "ephemeris/ephemeris.h"
#include "load_data.h"
#include "minor_ref.h"

int main(int argc, char** argv) {
    if (!load_data_args(argc, argv)) return 1;
    int fails = 0;
    double worst = 0.0;
    int worst_number = 0;
    for (const MinorRef& r : kMinorRef) {
        char name[16];
        std::snprintf(name, sizeof name, "%d", r.number);
        const int i = eph_minor_find(name);
        eph_position p;
        if (i < 0 || eph_position_apparent_tt((eph_body)(EPH_MINOR_BASE + i), r.jd, &p) != 0) {
            std::printf("FAIL: %d at %.1f not available\n", r.number, r.jd);
            fails++;
            continue;
        }
        double ra, dec;
        eph_ecliptic_to_equatorial(p.longitude, p.latitude, r.jd, &ra, &dec);
        const double R = M_PI / 180;
        const double d = std::acos(std::fmin(1.0, std::sin(dec) * std::sin(r.dec * R) +
                                                      std::cos(dec) * std::cos(r.dec * R) *
                                                          std::cos(ra - r.ra * R))) /
                         R * 3600;
        if (d > worst) worst = d, worst_number = r.number;
        if (d > 0.5) {
            std::printf("FAIL: %d at %.1f off by %.2f\"\n", r.number, r.jd, d);
            fails++;
        }
    }
    std::printf("Minor planets vs Horizons: worst %.3f\" (%d), %zu positions\n", worst,
                worst_number, sizeof kMinorRef / sizeof *kMinorRef);

    /* Lookup and names. */
    const int eros = eph_minor_find("Eros");
    if (eros < 0 || eph_minor_find("433") != eros || eph_minor_find("433 eros") != eros ||
        eph_minor_find("(433) Eros") != eros || eph_minor_find("Nonexistent") != -1 ||
        eph_minor_number(eros) != 433 || std::string_view(eph_minor_name(eros)) != "Eros" ||
        std::string_view(eph_body_name((eph_body)(EPH_MINOR_BASE + eros))) != "Eros" ||
        eph_minor_name(eph_minor_count()) != nullptr) {
        std::printf("FAIL: lookup\n");
        fails++;
    }
    /* Phenomena use the file's H, G and radius; searches and the
     * heliocentric position work. */
    eph_phenomena ph;
    double t, lon;
    if (eph_phenomena_tt((eph_body)(EPH_MINOR_BASE + eros), 2461324.0, &ph) != 0 ||
        !(ph.magnitude > 5 && ph.magnitude < 20) || !(ph.diameter > 0) ||
        eph_heliocentric_tt((eph_body)(EPH_MINOR_BASE + eros), 2461324.0, &lon, nullptr,
                            nullptr) != 0 ||
        eph_longitude_crossing((eph_body)(EPH_MINOR_BASE + eros), 0.0, 2461324.0, 1, 3000,
                               &t) != 0) {
        std::printf("FAIL: phenomena / heliocentric / search\n");
        fails++;
    }
    /* Out of range. */
    eph_position p;
    if (eph_position_tt((eph_body)(EPH_MINOR_BASE + eros), 2378496.5, &p) == 0 ||
        eph_position_tt((eph_body)(EPH_MINOR_BASE + eph_minor_count()), 2451545.0, &p) == 0) {
        std::printf("FAIL: out of range\n");
        fails++;
    }

    if (fails == 0) std::printf("All minor planet checks passed.\n");
    return fails ? 1 : 0;
}
