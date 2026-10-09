/* Helper for tools/pack_ephemeris.py: residual of a reference ephemeris
 * against the library's truncated series, exactly as src/ephemeris.cpp
 * evaluates them.
 *
 *     residual BODY BASE THRESHOLD JD0 JD1 < "jd x y z" lines
 *
 * Input: mean ecliptic J2000 (IAU-1976) vectors in AU, heliocentric or
 * geocentric for the Moon (10).  BASE 0 means no base series; 1 is
 * supported for the Moon and Mercury..Neptune except the Earth.  JD0/JD1 is
 * the table's advertised range (it sets the truncation).  Output: the
 * residual "x y z" in the ecliptic of date.  Built by pack_ephemeris.py. */
#include <cstdio>
#include <cstdlib>

#include "internal.h"

using namespace eph;

static Vec3 j2000_to_date(double jd, const Vec3& v) {
    double l, b, r;
    cartesian_to_spherical(v, l, b, r);
    precess_ecliptic_j2000_to_date(jd, l, b);
    return spherical_to_cartesian(l, b, r);
}

int main(int argc, char** argv) {
    if (argc != 6) {
        std::fprintf(stderr, "usage: residual BODY BASE THRESHOLD JD0 JD1\n");
        return 1;
    }
    const int body = std::atoi(argv[1]);
    const bool base = std::atoi(argv[2]) != 0;
    const double thr = std::atof(argv[3]);
    const double jd0 = std::atof(argv[4]), jd1 = std::atof(argv[5]);
    const double tmax = std::fmax(std::fabs(jd0 - J2000), std::fabs(jd1 - J2000));
    if (base) {
        if (body == 10)
            elp_truncate(thr, tmax);
        else if (body <= 7)
            vsop87_truncate(body, thr, tmax);
    }
    double jd;
    Vec3 x;
    while (std::scanf("%lf %lf %lf %lf", &jd, &x.x, &x.y, &x.z) == 4) {
        Vec3 r = j2000_to_date(jd, x);
        if (base) {
            Vec3 b;
            if (body == 10) {
                double bx, by, bz;
                elp_truncated_j2000(jd, bx, by, bz);
                b = j2000_to_date(jd, {bx / AU_KM, by / AU_KM, bz / AU_KM});
            } else {
                double l, la, rr;
                vsop87_truncated(body, jd, l, la, rr);
                b = spherical_to_cartesian(l, la, rr);
            }
            r = r - b;
        }
        std::printf("%.17g %.17g %.17g\n", r.x, r.y, r.z);
    }
    return 0;
}
