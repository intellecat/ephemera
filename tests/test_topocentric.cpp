/* Topocentric (observer) validation against Skyfield/DE440s references. */
#include <cmath>
#include <cstdio>

#include "ephemeris/ephemeris.h"
#include "topo_ref.h"

static const double R2D = 180.0 / 3.14159265358979323846;
static const double D2R = 3.14159265358979323846 / 180.0;

static double wrap180(double d) {
    d = std::fmod(d + 180.0, 360.0);
    if (d < 0.0) d += 360.0;
    return d - 180.0;
}

int main() {
    double maxang = 0.0, maxdist = 0.0;
    for (const TopoRef& r : kTopoRef) {
        eph_observer obs = {r.lon * D2R, r.lat * D2R, r.elev};
        eph_position p;
        if (eph_position_topocentric_ut((eph_body)r.body, r.jd_ut, r.year, &obs,
                                        &p) != 0) {
            std::printf("FAIL: body %d returned error\n", r.body);
            return 1;
        }
        const double dlon =
            std::fabs(wrap180(p.longitude * R2D - r.lon_deg)) * 3600.0;
        const double dlat = std::fabs(p.latitude * R2D - r.lat_deg) * 3600.0;
        const double dd = std::fabs(p.distance - r.dist_au);
        maxang = std::fmax(maxang, std::fmax(dlon, dlat));
        maxdist = std::fmax(maxdist, dd);
    }
    std::printf("\nTopocentric vs Skyfield/DE440s: max angular %.3f\", "
                "max distance %.2e AU\n", maxang, maxdist);
    if (maxang > 2.0) {
        std::printf("FAIL: angular error %.3f\" > 2\"\n", maxang);
        return 1;
    }
    std::printf("All topocentric checks passed.\n");
    return 0;
}
