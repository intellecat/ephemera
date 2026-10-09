/* Sanity checks for ayanamsa, Ascendant/Midheaven and houses. */
#include <cmath>
#include <cstdio>
#include <initializer_list>

#include "ephemeris/ephemeris.h"
#include "houses_ref.h"

static const double R2D = 180.0 / 3.14159265358979323846;
static const double D2R = 3.14159265358979323846 / 180.0;
static const double PI = 3.14159265358979323846;

/* Angle difference in degrees, wrapped to [-180, 180). */
static double wrap180(double d) {
    d = std::fmod(d + 180.0, 360.0);
    return (d < 0.0 ? d + 360.0 : d) - 180.0;
}

int main() {
    int fails = 0;
    const double jd = 2451545.0; /* 2000-01-01 12h UT */

    /* Ayanamsas (degrees) vs Swiss Ephemeris (swe_get_ayanamsa) for
     * Lahiri, Fagan/Bradley, Raman, Krishnamurti, Yukteshwar. */
    const struct {
        double jd;
        double deg[5];
    } ayan[] = {
        {2415020.5, {22.460530593, 23.343738216, 21.014229126, 22.363678126, 21.082241114}},
        {2451545.0, {23.857092325, 24.740299966, 22.410791012, 23.760240012, 22.478803000}},
        {2469807.5, {24.555613070, 25.438820717, 23.109311808, 24.458760808, 23.177323796}},
    };
    for (const auto& a : ayan)
        for (int m = 0; m < 5; ++m) {
            const double d =
                eph_ayanamsa((eph_ayanamsa_system)m, a.jd) * R2D - a.deg[m];
            if (std::fabs(d) * 3600.0 > 0.01) {
                std::printf("FAIL: ayanamsa %d at %.1f off by %.4f\"\n", m,
                            a.jd, d * 3600.0);
                fails++;
            }
        }
    if (std::fabs(eph_ayanamsa_lahiri(jd) -
                  eph_ayanamsa(EPH_AYANAMSA_LAHIRI, jd)) > 0 ||
        !std::isnan(eph_ayanamsa((eph_ayanamsa_system)99, jd))) {
        std::printf("FAIL: ayanamsa API\n");
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
        if (std::fabs(wrap180((cusps[0] - asc) * R2D)) > 1e-6 ||
            std::fabs(wrap180((cusps[9] - mc) * R2D)) > 1e-6 ||
            std::fabs(wrap180((cusps[3] - mc - PI) * R2D)) > 1e-6) {
            std::printf("FAIL: porphyry angular cusps (lat %.4f)\n", o.lat);
            fails++;
        }
    }

    /* Asc/MC vs Swiss Ephemeris 2.10 (swe_houses_ex), UT. */
    struct {
        double jd, lat, lon, asc, mc;
    } ref[] = {
        {2451545.0, 51.5, 0.0, 24.2873051, 279.6110878},
        {2460676.5, 40.7, -74.0, 130.3345811, 28.9401645},
        {2447892.75, -33.9, 151.2, 59.7829494, 340.3196221},
        {2469807.3, 64.1, -21.9, 133.0235458, 7.3529408},
    };
    for (const auto& r : ref) {
        eph_observer ob = {r.lon * D2R, r.lat * D2R, 0.0};
        double asc, mc;
        eph_asc_mc(r.jd, &ob, &asc, &mc);
        const double da = std::fabs(std::remainder(asc * R2D - r.asc, 360.0));
        const double dm = std::fabs(std::remainder(mc * R2D - r.mc, 360.0));
        if (da * 3600.0 > 1.0 || dm * 3600.0 > 1.0) {
            std::printf("FAIL: asc/mc off by %.2f\"/%.2f\" (jd %.2f)\n",
                        da * 3600.0, dm * 3600.0, r.jd);
            fails++;
        }
    }

    /* All systems vs Swiss Ephemeris: 1900-2030, |latitude| <= 65 degrees
     * and polar cases for the systems defined there. */
    double worst[14] = {0};
    for (const HouseRef& r : kHouseRef) {
        eph_observer ob = {r.lon * D2R, r.lat * D2R, 0.0};
        double c[12];
        if (eph_houses(r.jd_ut, &ob, (eph_house_system)r.system, c) != 0) {
            std::printf("FAIL: system %d undefined at latitude %.2f\n",
                        r.system, r.lat);
            fails++;
            continue;
        }
        for (int i = 0; i < 12; ++i)
            worst[r.system] = std::fmax(
                worst[r.system],
                std::fabs(wrap180(c[i] * R2D - r.cusp[i])) * 3600.0);
    }
    static const char* kNames[14] = {
        "Whole Sign", "Equal",    "Porphyry",    "Placidus",   "Koch",
        "Regiomontanus", "Campanus", "Alcabitius", "Topocentric", "Vehlow",
        "Equal (MC)", "Equal (Aries)", "Meridian", "Morinus"};
    for (int s = 1; s < 14; ++s) {
        std::printf("%-14s max %.3f\" vs Swiss Ephemeris\n", kNames[s], worst[s]);
        if (worst[s] > 0.5) {
            std::printf("FAIL: %s\n", kNames[s]);
            fails++;
        }
    }

    /* Angles and sensitive points vs Swiss Ephemeris. */
    double worst_angle = 0.0;
    for (const AngleRef& r : kAngleRef) {
        eph_observer ob = {r.lon * D2R, r.lat * D2R, 0.0};
        eph_angles a;
        eph_house_angles(r.jd_ut, &ob, &a);
        const double v[8] = {a.ascendant, a.mc, a.armc, a.vertex,
                             a.equatorial_ascendant, a.co_ascendant_koch,
                             a.co_ascendant_munkasey, a.polar_ascendant};
        for (int i = 0; i < 8; ++i)
            worst_angle = std::fmax(
                worst_angle, std::fabs(wrap180(v[i] * R2D - r.value[i])) * 3600.0);
    }
    std::printf("Angles         max %.3f\" vs Swiss Ephemeris\n", worst_angle);
    if (worst_angle > 0.5) {
        std::printf("FAIL: angles\n");
        fails++;
    }

    /* Placidus and Koch are undefined inside the polar circles: Porphyry
     * cusps and return code 1. */
    {
        eph_observer ob = {0.0, 70.0 * D2R, 0.0};
        double c[12], p[12];
        eph_houses(jd, &ob, EPH_HOUSES_PORPHYRY, p);
        for (eph_house_system s : {EPH_HOUSES_PLACIDUS, EPH_HOUSES_KOCH}) {
            const int rc = eph_houses(jd, &ob, s, c);
            if (rc != 1 || std::fabs(wrap180((c[1] - p[1]) * R2D)) > 1e-9) {
                std::printf("FAIL: polar fallback for system %d (rc %d)\n", s, rc);
                fails++;
            }
        }
    }

    /* House positions: the start of each house, and halfway in Equal. */
    {
        eph_observer ob = {0.3, 0.8, 0.0};
        double c[12];
        eph_houses(jd, &ob, EPH_HOUSES_PLACIDUS, c);
        for (int i = 0; i < 12; ++i) {
            const double h = eph_house_position(c, c[i] + 1e-9);
            if (std::fabs(h - (i + 1)) > 1e-6) {
                std::printf("FAIL: house position %.6f at cusp %d\n", h, i + 1);
                fails++;
            }
        }
        eph_houses(jd, &ob, EPH_HOUSES_EQUAL, c);
        const double h = eph_house_position(c, c[4] + 15.0 * D2R);
        if (std::fabs(h - 5.5) > 1e-9) {
            std::printf("FAIL: house position %.6f, expected 5.5\n", h);
            fails++;
        }
    }

    if (fails == 0) std::printf("All house checks passed.\n");
    return fails == 0 ? 0 : 1;
}
