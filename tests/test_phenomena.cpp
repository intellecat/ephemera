/* Phenomena vs JPL DE440 with Skyfield (tests/phenomena_ref.h, made by
 * tools/make_phenomena_reference.py).
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
    double worst[5] = {0}, worst_mag_outer = 0.0;
    for (const PhenoRef& r : kPhenoRef) {
        eph_phenomena p;
        if (eph_phenomena_tt((eph_body)r.body, r.jd_tt, &p) != 0) {
            std::printf("FAIL: no phenomena for body %d\n", r.body);
            return 1;
        }
        const double e[4] = {
            std::fabs(p.phase_angle * R2D - r.phase_angle) * 3600.0,
            std::fabs(p.phase - r.phase),
            std::fabs(p.elongation * R2D - r.elongation) * 3600.0,
            std::fabs(p.diameter * R2D - r.diameter) / r.diameter};
        for (int i = 0; i < 4; ++i) worst[i] = std::fmax(worst[i], e[i]);
        /* No reference magnitude for the Moon and Pluto (NaN). */
        if (std::isnan(r.magnitude)) continue;
        const double dm = std::fabs(p.magnitude - r.magnitude);
        if (r.body == EPH_SATURN || r.body == EPH_URANUS || r.body == EPH_NEPTUNE)
            worst_mag_outer = std::fmax(worst_mag_outer, dm);
        else
            worst[4] = std::fmax(worst[4], dm);
    }
    std::printf("vs DE440: phase angle %.3f\", phase %.1e, elongation %.3f\", "
                "diameter %.1e (relative), magnitude %.4f (Saturn-Neptune %.3f)\n",
                worst[0], worst[1], worst[2], worst[3], worst[4], worst_mag_outer);
    /* Measured with eph-full.bin: 0.03", 5e-8, 0.05", 1.5e-7, <1e-4 mag. */
    const double limit[5] = {0.2, 1e-6, 0.2, 1e-6, 1e-4};
    for (int i = 0; i < 5; ++i)
        if (worst[i] > limit[i]) {
            std::printf("FAIL: quantity %d exceeds %g\n", i, limit[i]);
            fails++;
        }
    /* Saturn, Uranus and Neptune: ours use the Earth's sub-latitude only
     * (Mallama & Hilton average it with the Sun's) and apply the phase
     * term at all phase angles; up to 0.04 mag. */
    if (worst_mag_outer > 0.05) {
        std::printf("FAIL: Saturn-Neptune magnitudes\n");
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
