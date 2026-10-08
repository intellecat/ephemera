/* Validation of the analytical engine against a JPL DE440s reference.
 *
 * Reference values (tests/reference_data.h) were generated with Skyfield +
 * DE440s for dates spanning 1850-2140 and are compared here.  The reference
 * is reduced to the library's frame (mean ecliptic/equinox of date, IAU
 * 1976) so the test measures the dynamics rather than a frame difference. */
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

static int check_time() {
    int failures = 0;
    struct {
        int y, m;
        double d, jd;
    } cases[] = {
        {2000, 1, 1.5, 2451545.0},     /* J2000                          */
        {1957, 10, 4.81, 2436116.31},  /* Meeus, example 7.a (Gregorian) */
        {333, 1, 27.5, 1842713.0},     /* Meeus, example 7.b (Julian)    */
    };
    for (const auto& c : cases) {
        const double jd = eph_julian_date(c.y, c.m, c.d);
        if (std::fabs(jd - c.jd) > 1e-6) {
            std::printf("FAIL: JD(%d-%d-%.2f) = %.5f, expected %.5f\n", c.y,
                        c.m, c.d, jd, c.jd);
            failures++;
        }
    }
    int y, m;
    double d;
    eph_calendar_date(2436116.31, &y, &m, &d);
    if (y != 1957 || m != 10 || std::fabs(d - 4.81) > 1e-6) {
        std::printf("FAIL: calendar round-trip -> %d-%d-%.2f\n", y, m, d);
        failures++;
    }
    if (std::fabs(eph_delta_t(2000.0) - 63.9) > 2.0) {
        std::printf("FAIL: Delta-T(2000) = %.2f s\n", eph_delta_t(2000.0));
        failures++;
    }
    return failures;
}

int main() {
    int time_failures = check_time();

    double gmax[11] = {0}, amax[11] = {0};
    double gdist[11] = {0}, adist[11] = {0};
    int gcount[11] = {0}, acount[11] = {0};

    for (const RefSample& r : kRefSamples) {
        eph_position g, a;
        if (eph_position_tt((eph_body)r.body, r.tt_jd, &g) != 0) {
            std::printf("body %d failed (geometric)\n", r.body);
            return 1;
        }
        if (eph_position_apparent_tt((eph_body)r.body, r.tt_jd, &a) != 0) {
            std::printf("body %d failed (apparent)\n", r.body);
            return 1;
        }

        const double glon = g.longitude * R2D, glat = g.latitude * R2D;
        const double alon = a.longitude * R2D, alat = a.latitude * R2D;

        /* Reference geometric is the mean ecliptic/equinox of date under the
         * IAU-1976 reduction, matching VSOP87D and this library. */
        const double e_glon = std::fabs(wrap180(glon - r.geom_lon)) * 3600.0;
        const double e_glat = std::fabs(glat - r.geom_lat) * 3600.0;
        const double e_alon = std::fabs(wrap180(alon - r.app_lon)) * 3600.0;
        const double e_alat = std::fabs(alat - r.app_lat) * 3600.0;
        const double e_gdist = std::fabs(g.distance - r.geom_dist);
        const double e_adist = std::fabs(a.distance - r.app_dist);

        gmax[r.body] = std::fmax(gmax[r.body], std::fmax(e_glon, e_glat));
        amax[r.body] = std::fmax(amax[r.body], std::fmax(e_alon, e_alat));
        gdist[r.body] = std::fmax(gdist[r.body], e_gdist);
        adist[r.body] = std::fmax(adist[r.body], e_adist);
        gcount[r.body]++;
        acount[r.body]++;

        if (e_glon > 300.0 || e_glat > 300.0 || e_alon > 300.0 ||
            e_alat > 300.0) {
            std::printf(
                "  LARGE %-8s jd=%.5f geomLon=%.4f\" geomLat=%.4f\" "
                "appLon=%.4f\" appLat=%.4f\"\n",
                eph_body_name((eph_body)r.body), r.tt_jd, e_glon, e_glat,
                e_alon, e_alat);
        }
    }

    std::printf("\nMaximum angular residuals vs DE440s (arcsec)\n");
    std::printf("%-8s %14s %14s\n", "Body", "geometric", "apparent");
    for (int i = 0; i < EPH_BODY_COUNT; ++i) {
        if (gcount[i] == 0) continue;
        std::printf("%-8s %14.3f %14.3f\n", eph_body_name((eph_body)i),
                    gmax[i], amax[i]);
    }
    std::printf("\nMaximum distance residuals (AU)\n");
    for (int i = 0; i < EPH_BODY_COUNT; ++i) {
        if (gcount[i] == 0) continue;
        std::printf("%-8s geom %.3e  app %.3e\n", eph_body_name((eph_body)i),
                    gdist[i], adist[i]);
    }

    /* Accuracy gates.  Inner planets are tight; outer planets and the Moon
     * are looser because of the theories used. */
    static const double glimit[11] = {5, 5, 5, 10, 15, 20, 15, 15, 2, 5, 30};
    static const double alimit[11] = {5, 5, 5, 10, 15, 20, 15, 15, 2, 5, 30};
    int failures = 0;
    for (int i = 0; i < EPH_BODY_COUNT; ++i) {
        if (gcount[i] == 0) continue;
        if (gmax[i] > glimit[i]) {
            std::printf("FAIL: %s geometric residual %.3f\" > %.0f\"\n",
                        eph_body_name((eph_body)i), gmax[i], glimit[i]);
            failures++;
        }
        if (amax[i] > alimit[i]) {
            std::printf("FAIL: %s apparent residual %.3f\" > %.0f\"\n",
                        eph_body_name((eph_body)i), amax[i], alimit[i]);
            failures++;
        }
    }

    failures += time_failures;
    if (failures == 0) std::printf("\nAll checks passed.\n");
    return failures == 0 ? 0 : 1;
}
