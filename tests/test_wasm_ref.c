/* Native reference for tests/test_wasm.mjs: apparent positions at J2000 TT
 * (degrees, degrees/day) with the given data file. */
#include <stdio.h>

#include "ephemeris/ephemeris.h"

int main(int argc, char** argv) {
    const double deg = 180.0 / 3.14159265358979323846;
    if (argc != 2 || eph_load_file(argv[1]) != 0) return 1;
    for (int b = 0; b <= EPH_MOON; ++b) {
        eph_position p;
        if (b == EPH_EARTH || eph_position_apparent_tt((eph_body)b, 2451545.0, &p))
            continue;
        printf("%d %.15g %.15g %.15g\n", b, p.longitude * deg, p.latitude * deg,
               p.longitude_speed * deg);
    }
    return 0;
}
