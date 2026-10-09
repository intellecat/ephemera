/* Validation of the optional Chebyshev compressed backend against DE440s.
 * Usage: test_compressed <path to ephem.bin> */
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
        std::printf("usage: %s <ephem.bin>\n", argv[0]);
        return 1;
    }
    if (eph_load_compressed(argv[1]) != 0) {
        std::printf("SKIP: could not load %s\n", argv[1]);
        return 0; /* not a failure: the backend is optional */
    }
    if (!eph_compressed_available()) {
        std::printf("SKIP: backend unavailable\n");
        return 0;
    }

    double gmax[EPH_BODY_COUNT] = {0}, amax[EPH_BODY_COUNT] = {0};
    for (const RefSample& r : kRefSamples) {
        eph_position g, a;
        eph_position_tt((eph_body)r.body, r.tt_jd, &g);
        eph_position_apparent_tt((eph_body)r.body, r.tt_jd, &a);
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
        if (gmax[i] > 0.08 || amax[i] > 0.08) {
            std::printf("FAIL: %s exceeds 0.08\"\n",
                        eph_body_name((eph_body)i));
            failures++;
        }
    }
    if (failures == 0) std::printf("\nAll compressed checks passed.\n");
    return failures == 0 ? 0 : 1;
}
