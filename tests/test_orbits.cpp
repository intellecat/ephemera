/* Planetary nodes and apsides: osculating elements and orbit points against
 * DE440 (tests/orbit_ref.h, tools/make_orbit_reference.py), run with the
 * DE440 tables (with the series alone the outer planets' elements carry the
 * series' ~1e-6 AU noise: Neptune 3", Pluto 11"); the mean elements
 * against Meeus' worked example 31.a. */
#include <cmath>
#include <cstdio>
#include <initializer_list>

#include "ephemeris/ephemeris.h"
#include "load_data.h"
#include "orbit_ref.h"

namespace {

const double R2D = 180.0 / M_PI;
int fails = 0;

double wrap(double d) { return std::remainder(d, 360.0); }

}  // namespace

int main(int argc, char** argv) {
    if (!load_data_args(argc, argv)) return 1;

    /* Osculating: the worst angle (arcseconds) per body. */
    static const char* kNames[9] = {"Mercury", "Venus", "", "Mars", "Jupiter",
                                    "Saturn", "Uranus", "Neptune", "Pluto"};
    double worst_angle[9] = {0}, worst_a[9] = {0}, worst_e[9] = {0}, worst_geo[9] = {0};
    for (const OrbitRef& r : kOrbitRef) {
        eph_orbit o;
        if (eph_orbit_tt((eph_body)r.body, r.jd, EPH_ORBIT_OSCULATING, &o) != 0) {
            std::printf("FAIL: no osculating orbit for body %d\n", r.body);
            return 1;
        }
        const double d[3] = {o.inclination * R2D - r.i, wrap(o.node * R2D - r.node),
                             wrap(o.perihelion * R2D - r.peri)};
        /* The node and perihelion directions are weighted by sin i and e:
         * an error in them moves the orbit by that much less. */
        const double w[3] = {1.0, std::sin(r.i / R2D), r.e};
        for (int k = 0; k < 3; ++k)
            worst_angle[r.body] = std::fmax(worst_angle[r.body], std::fabs(d[k]) * w[k] * 3600);
        worst_a[r.body] = std::fmax(worst_a[r.body], std::fabs(o.semi_major_axis / r.a - 1));
        worst_e[r.body] = std::fmax(worst_e[r.body], std::fabs(o.eccentricity - r.e));
        for (int k = 0; k < 4; ++k)
            worst_geo[r.body] = std::fmax(
                worst_geo[r.body],
                std::fabs(wrap(o.geocentric[k].longitude * R2D - r.geo[k])) * 3600 *
                    (k < 2 ? std::sin(r.i / R2D) : r.e));
    }
    for (int b = 0; b < 9; ++b) {
        if (!*kNames[b]) continue;
        std::printf("%-8s angles %.3f\"  a %.1e  e %.1e  geocentric points %.3f\"\n",
                    kNames[b], worst_angle[b], worst_a[b], worst_e[b], worst_geo[b]);
        if (worst_angle[b] > 1.0 || worst_a[b] > 1e-5 || worst_e[b] > 1e-5 || worst_geo[b] > 1.0) {
            std::printf("FAIL: %s osculating orbit\n", kNames[b]);
            fails++;
        }
    }

    /* Mean elements: Meeus, Astronomical Algorithms, example 31.a (Mercury,
     * 2065-06-24 0h TD): e 0.20564510, i 7.006171, node 49.107650,
     * perihelion 78.475382 degrees. */
    {
        eph_orbit o;
        if (eph_orbit_tt(EPH_MERCURY, 2475460.5, EPH_ORBIT_MEAN, &o) != 0 ||
            std::fabs(o.eccentricity - 0.20564510) > 5e-9 ||
            std::fabs(o.inclination * R2D - 7.006171) > 5e-7 ||
            std::fabs(o.node * R2D - 49.107650) > 5e-7 ||
            std::fabs(o.perihelion * R2D - 78.475382) > 5e-7) {
            std::printf("FAIL: Mercury mean elements (Meeus example 31.a)\n");
            fails++;
        }
        /* The points: the ascending node at the node's longitude, on the
         * ecliptic; perihelion at distance a (1 - e). */
        const eph_position& n = o.heliocentric[EPH_ASCENDING_NODE];
        const eph_position& p = o.heliocentric[EPH_PERIHELION];
        if (std::fabs(wrap((n.longitude - o.node) * R2D)) > 1e-9 || std::fabs(n.latitude) > 1e-12 ||
            std::fabs(p.distance - o.semi_major_axis * (1 - o.eccentricity)) > 1e-12) {
            std::printf("FAIL: orbit points\n");
            fails++;
        }
    }
    /* Mean and osculating orbits agree to the size of the perturbations. */
    for (int b : {0, 1, 3, 4, 5, 6, 7}) {
        eph_orbit m, o;
        eph_orbit_tt((eph_body)b, 2451545.0, EPH_ORBIT_MEAN, &m);
        eph_orbit_tt((eph_body)b, 2451545.0, EPH_ORBIT_OSCULATING, &o);
        const double dn = std::fabs(wrap((m.node - o.node) * R2D)) * std::sin(m.inclination);
        if (dn > 0.1 || std::fabs(m.inclination - o.inclination) * R2D > 0.1) {
            std::printf("FAIL: %s mean vs osculating node %.3f deg\n", kNames[b], dn);
            fails++;
        }
    }
    /* Unsupported: the Earth, the Moon, the Sun; mean elements for Pluto. */
    eph_orbit o;
    if (eph_orbit_tt(EPH_EARTH, 2451545.0, EPH_ORBIT_OSCULATING, &o) != -1 ||
        eph_orbit_tt(EPH_MOON, 2451545.0, EPH_ORBIT_MEAN, &o) != -1 ||
        eph_orbit_tt(EPH_SUN, 2451545.0, EPH_ORBIT_OSCULATING, &o) != -1 ||
        eph_orbit_tt(EPH_PLUTO, 2451545.0, EPH_ORBIT_MEAN, &o) != -1 ||
        eph_orbit_tt(EPH_PLUTO, 2451545.0, EPH_ORBIT_OSCULATING, &o) != 0) {
        std::printf("FAIL: supported bodies\n");
        fails++;
    }

    if (fails == 0) std::printf("All orbit checks passed.\n");
    return fails == 0 ? 0 : 1;
}
