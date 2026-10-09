/* Validation of the DE440 tables (data/de440/) on their own.
 * Usage: test_compressed <DE440 table files covering 1850-2140> */
#include <cmath>
#include <cstdio>

#include "ephemeris/ephemeris.h"
#include "reference_data.h"

static const double R2D = 180.0 / 3.14159265358979323846;
static double wrap180(double d) {
    d = std::fmod(d + 180.0, 360.0);
    if (d < 0.0) d += 360.0;
    return d - 180.0;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::printf("usage: %s <de440-YEAR.bin>...\n", argv[0]);
        return 1;
    }
    /* The tables are self-contained: no series file is loaded. */
    for (int i = 1; i < argc; ++i)
        if (eph_load_file(argv[i]) != 0) {
            std::printf("FAIL: could not load %s\n", argv[i]);
            return 1;
        }

    double gmax[EPH_BODY_COUNT] = {0}, amax[EPH_BODY_COUNT] = {0};
    for (const RefSample& r : kRefSamples) {
        eph_position g, a;
        if (eph_position_tt((eph_body)r.body, r.tt_jd, &g) != 0 ||
            eph_position_apparent_tt((eph_body)r.body, r.tt_jd, &a) != 0) {
            std::printf("FAIL: no position for %s at JD %.1f\n",
                        eph_body_name((eph_body)r.body), r.tt_jd);
            return 1;
        }
        const double e_glon =
            std::fabs(wrap180(g.longitude * R2D - r.geom_lon)) * 3600.0;
        const double e_glat = std::fabs(g.latitude * R2D - r.geom_lat) * 3600.0;
        const double e_alon =
            std::fabs(wrap180(a.longitude * R2D - r.app_lon)) * 3600.0;
        const double e_alat = std::fabs(a.latitude * R2D - r.app_lat) * 3600.0;
        gmax[r.body] = std::fmax(gmax[r.body], std::fmax(e_glon, e_glat));
        amax[r.body] = std::fmax(amax[r.body], std::fmax(e_alon, e_alat));
    }

    std::printf("\nCompressed backend, max residual vs DE440s (arcsec)\n");
    std::printf("%-8s %12s %12s\n", "Body", "geometric", "apparent");
    int failures = 0;
    for (int i = 0; i < EPH_BODY_COUNT; ++i) {
        if (i == EPH_EARTH || gmax[i] == 0.0) continue;
        std::printf("%-8s %12.4f %12.4f\n", eph_body_name((eph_body)i),
                    gmax[i], amax[i]);
        /* The table is fitted to 0.002"; the apparent place adds the
         * reduction's own error (IAU 2000B nutation, ~1-3 mas). */
        if (gmax[i] > 0.005 || amax[i] > 0.006) {
            std::printf("FAIL: %s exceeds 0.005\" / 0.006\"\n",
                        eph_body_name((eph_body)i));
            failures++;
        }
    }
    if (failures == 0) std::printf("\nAll compressed checks passed.\n");
    return failures == 0 ? 0 : 1;
}
