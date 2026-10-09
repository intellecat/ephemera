/* Sanity checks for ayanamsa, Ascendant/Midheaven and houses. */
#include <cmath>
#include <cstdio>

#include "ephemeris/ephemeris.h"

static const double R2D = 180.0 / 3.14159265358979323846;
static const double D2R = 3.14159265358979323846 / 180.0;
static const double PI = 3.14159265358979323846;

static double wrap360(double d) {
    d = std::fmod(d, 360.0);
    return d < 0.0 ? d + 360.0 : d;
}

int main() {
    int fails = 0;
    const double jd = 2451545.0; /* 2000-01-01 12h UT */

    /* Ayanamsa (Lahiri) ~ 23.853 deg at J2000. */
    const double ayan = eph_ayanamsa_lahiri(jd) * R2D;
    if (std::fabs(ayan - 23.8532) > 0.01) {
        std::printf("FAIL: ayanamsa %.4f deg\n", ayan);
        fails++;
    }

    struct {
        double lat, lon;
    } obs[] = {{40.7128, -74.0060}, {-33.8688, 151.2093}, {51.5074, -0.1278}};

    for (const auto& o : obs) {
        eph_observer ob = {o.lon * D2R, o.lat * D2R, 0.0};
        double asc, mc;
        eph_asc_mc(jd, &ob, &asc, &mc);

        /* The MC's right ascension is the local sidereal time; the Asc must
         * lie on the horizon (altitude 0). */
        double ra_mc, dec_mc, ra_asc, dec_asc;
        eph_ecliptic_to_equatorial(mc, 0.0, jd, &ra_mc, &dec_mc);
        eph_ecliptic_to_equatorial(asc, 0.0, jd, &ra_asc, &dec_asc);
        const double H = ra_mc - ra_asc; /* LST - RA_asc */
        const double alt =
            std::asin(std::sin(ob.latitude) * std::sin(dec_asc) +
                      std::cos(ob.latitude) * std::cos(dec_asc) * std::cos(H)) *
            R2D;
        if (std::fabs(alt) > 0.01) {
            std::printf("FAIL: asc altitude %.4f deg (lat %.4f)\n", alt, o.lat);
            fails++;
        }

        /* Porphyry cusp relations. */
        double cusps[12];
        eph_houses(jd, &ob, EPH_HOUSES_PORPHYRY, cusps);
        if (std::fabs(wrap360((cusps[0] - asc) * R2D)) > 1e-6 ||
            std::fabs(wrap360((cusps[9] - mc) * R2D)) > 1e-6 ||
            std::fabs(wrap360((cusps[3] - mc - PI) * R2D)) > 1e-6) {
            std::printf("FAIL: porphyry angular cusps (lat %.4f)\n", o.lat);
            fails++;
        }
    }

    if (fails == 0) std::printf("All house checks passed.\n");
    return fails == 0 ? 0 : 1;
}
