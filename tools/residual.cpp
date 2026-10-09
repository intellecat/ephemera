/* Helper for tools/pack_ephemeris.py: residual of a reference ephemeris
 * against the base series of a residual table, evaluated exactly as
 * src/ephemeris.cpp does.
 *
 *     residual BODY [BASE_FILE] < "jd x y z" lines
 *
 * Input: mean ecliptic J2000 (IAU-1976) vectors in AU, heliocentric or
 * geocentric for the Moon (10).  BASE_FILE is a data file with the base
 * series (pack_series base); without it the residual is the input itself.
 * Output: the residual "x y z" in the ecliptic of date.  Built by
 * pack_ephemeris.py against build/libephemeris.a. */
#include <cstdio>
#include <cstdlib>

#include "ephemeris/ephemeris.h"
#include "internal.h"

using namespace eph;

static Vec3 j2000_to_date(double jd, const Vec3& v) {
    double l, b, r;
    cartesian_to_spherical(v, l, b, r);
    precess_ecliptic_j2000_to_date(jd, l, b);
    return spherical_to_cartesian(l, b, r);
}

int main(int argc, char** argv) {
    if (argc != 2 && argc != 3) {
        std::fprintf(stderr, "usage: residual BODY [BASE_FILE]\n");
        return 1;
    }
    const int body = std::atoi(argv[1]);
    const bool base = (argc == 3);
    if (base && eph_load_file(argv[2]) != 0) {
        std::fprintf(stderr, "residual: cannot load %s\n", argv[2]);
        return 1;
    }
    double jd;
    Vec3 x;
    while (std::scanf("%lf %lf %lf %lf", &jd, &x.x, &x.y, &x.z) == 4) {
        Vec3 r = j2000_to_date(jd, x);
        if (base) {
            Vec3 b;
            bool ok;
            if (body == 10) {
                double bx, by, bz;
                ok = elp_base_j2000(jd, bx, by, bz);
                b = j2000_to_date(jd, {bx / AU_KM, by / AU_KM, bz / AU_KM});
            } else {
                double l, la, rr;
                ok = vsop87_base(body, jd, l, la, rr);
                b = spherical_to_cartesian(l, la, rr);
            }
            if (!ok) {
                std::fprintf(stderr, "residual: no base series for %d\n", body);
                return 1;
            }
            r = r - b;
        }
        std::printf("%.17g %.17g %.17g\n", r.x, r.y, r.z);
    }
    return 0;
}
