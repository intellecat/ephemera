/* Asteroid validation against JPL Horizons apparent geocentric positions. */
#include <cmath>
#include <cstdio>

#include "ephemeris/ephemeris.h"
#include "asteroid_ref.h"

static const double R2D = 180.0 / 3.14159265358979323846;
static double wrap180(double d) {
    d = std::fmod(d + 180.0, 360.0);
    if (d < 0.0) d += 360.0;
    return d - 180.0;
}

int main() {
    double maxlon = 0.0, maxlat = 0.0;
    for (const AstRef& r : kAstRef) {
        eph_position p;
        if (eph_position_apparent_ut((eph_body)r.body, r.jd_ut, r.year, &p) != 0) {
            std::printf("FAIL: body %d error\n", r.body);
            return 1;
        }
        maxlon = std::fmax(maxlon,
                           std::fabs(wrap180(p.longitude * R2D - r.lon_deg)) *
                               3600.0);
        maxlat = std::fmax(maxlat,
                           std::fabs(p.latitude * R2D - r.lat_deg) * 3600.0);
    }
    std::printf("\nAsteroids vs Horizons: max dlon %.3f\", dlat %.3f\"\n",
                maxlon, maxlat);
    if (maxlon > 0.5 || maxlat > 0.5) {
        std::printf("FAIL: exceeds 0.5\"\n");
        return 1;
    }
    std::printf("All asteroid checks passed.\n");
    return 0;
}
