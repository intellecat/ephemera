/* Phenomena vs Swiss Ephemeris (tests/phenomena_ref.h).
 * Usage: test_phenomena <series data file> <asteroid data file> */
#include <cmath>
#include <cstdio>

#include "ephemeris/ephemeris.h"
#include "load_data.h"
#include "phenomena_ref.h"

static const double R2D = 180.0 / 3.14159265358979323846;

int main(int argc, char** argv) {
    if (!load_data_args(argc, argv)) return 1;
    int fails = 0;
    double worst[5] = {0};
    for (const PhenoRef& r : kPhenoRef) {
        eph_phenomena p;
        if (eph_phenomena_tt((eph_body)r.body, r.jd_tt, &p) != 0) {
            std::printf("FAIL: no phenomena for body %d\n", r.body);
            return 1;
        }
        const double e[5] = {
            std::fabs(p.phase_angle * R2D - r.phase_angle) * 3600.0,
            /* Swiss Ephemeris gives the Sun phase 0; ours is 1 (fully lit). */
            r.body == EPH_SUN ? 0.0 : std::fabs(p.phase - r.phase),
            std::fabs(p.elongation * R2D - r.elongation) * 3600.0,
            std::fabs(p.diameter * R2D - r.diameter) / r.diameter,
            /* The Moon's magnitude formula holds to phase angles ~150 deg. */
            (r.body == EPH_MOON && r.phase_angle > 150.0)
                ? 0.0 : std::fabs(p.magnitude - r.magnitude)};
        for (int i = 0; i < 5; ++i) worst[i] = std::fmax(worst[i], e[i]);
    }
    std::printf("vs Swiss Ephemeris: phase angle %.1f\", phase %.1e, elongation "
                "%.2f\", diameter %.1e (relative), magnitude %.3f\n",
                worst[0], worst[1], worst[2], worst[3], worst[4]);
    const double limit[5] = {30.0, 1e-4, 2.0, 1e-3, 0.15};
    for (int i = 0; i < 5; ++i)
        if (worst[i] > limit[i]) {
            std::printf("FAIL: quantity %d exceeds %g\n", i, limit[i]);
            fails++;
        }

    /* Minor bodies: plausible magnitudes; the lunar points have none. */
    eph_phenomena p;
    if (eph_phenomena_tt(EPH_CERES, 2451545.0, &p) != 0 || !(p.magnitude > 6.0) ||
        !(p.magnitude < 10.0) ||
        eph_phenomena_tt(EPH_CHIRON, 2451545.0, &p) != 0 || !(p.magnitude > 15.0) ||
        !(p.magnitude < 19.0)) {
        std::printf("FAIL: minor body magnitudes\n");
        fails++;
    }
    if (eph_phenomena_tt(EPH_TRUE_NODE, 2451545.0, &p) != -1 ||
        eph_phenomena_tt(EPH_EARTH, 2451545.0, &p) != -1) {
        std::printf("FAIL: points and the Earth have no phenomena\n");
        fails++;
    }
    if (fails == 0) std::printf("All phenomena checks passed.\n");
    return fails == 0 ? 0 : 1;
}
