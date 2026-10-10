/* Sanity checks for ayanamsa, Ascendant/Midheaven and houses. */
#include <cmath>
#include <cstdio>
#include <initializer_list>

#include "ephemeris/ephemeris.h"
#include "houses_ref.h"
#include "lahiri_ref.h"

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

    /* Table ayanamsas: each moves with the IAU 2006 general precession in
     * longitude (ERFA p06e) from its J2000 value, 1700-2400. */
    double worst_ayan = 0.0;
    const int ntable = EPH_AYANAMSA_TRUE_CITRA; /* star-based: test_stars */
    for (int m = 0; m < ntable; ++m) {
        const double a0 = eph_ayanamsa((eph_ayanamsa_system)m, 2451545.0) * R2D;
        for (int k = 0; k < 5; ++k) {
            const double d = std::fabs(
                eph_ayanamsa((eph_ayanamsa_system)m, kAyanamsaDates[k]) * R2D - a0 -
                kPrecession[k]) * 3600.0;
            worst_ayan = std::fmax(worst_ayan, d);
            if (d > 1e-4) {
                std::printf("FAIL: ayanamsa %s precession off by %.6f\"\n",
                            eph_ayanamsa_name((eph_ayanamsa_system)m), d);
                fails++;
            }
        }
    }
    std::printf("Ayanamsas      max %.6f\" vs IAU 2006 precession (%d systems)\n",
                worst_ayan, ntable);

    /* Definitions: J2000, J1900 (JD 2415020.0) and B1950 (JD 2433282.4235)
     * are zero at their epochs.  Lahiri: the Indian Astronomical Ephemeris
     * defines the mean ayanamsa as 23 deg 51' 25.53" at J2000.0 plus the
     * general precession (checked tightly below); the Calendar Reform
     * Committee's original 23 deg 15' 00.658" on 1956-03-21 0h (including
     * nutation) is 0.14" from it with IAU 2006 precession.
     * Fagan/Bradley's SVP at B1950.0 is 335 deg 57' 28.64" (ayanamsa
     * 24 deg 02' 31.36"); our value, the one in common use, is 0.41" from
     * that literal definition carried by IAU 2006 precession. */
    struct {
        eph_ayanamsa_system system;
        double jd, value, nutation, limit; /* degrees; limit in arcseconds */
    } anchors[] = {
        {EPH_AYANAMSA_J2000, 2451545.0, 0.0, 0.0, 1e-6},
        {EPH_AYANAMSA_LAHIRI, 2451545.0, 23.0 + 51.0 / 60 + 25.53 / 3600, 0.0, 0.005},
        {EPH_AYANAMSA_J1900, 2415020.0, 0.0, 0.0, 1e-3},
        {EPH_AYANAMSA_B1950, 2433282.4235, 0.0, 0.0, 1e-3},
        {EPH_AYANAMSA_LAHIRI, 2435553.5, 23.0 + 15.0 / 60 + 0.658 / 3600, kNutation1956,
         0.5},
        {EPH_AYANAMSA_FAGAN_BRADLEY, 2433282.4235, 24.0 + 2.0 / 60 + 31.36 / 3600, 0.0,
         0.5},
    };
    for (const auto& a : anchors) {
        const double d =
            (eph_ayanamsa(a.system, a.jd) * R2D + a.nutation - a.value) * 3600.0;
        std::printf("  %-14s %+.4f\" at its defining epoch\n", eph_ayanamsa_name(a.system),
                    d);
        if (std::fabs(d) > a.limit) {
            std::printf("FAIL: %s definition\n", eph_ayanamsa_name(a.system));
            fails++;
        }
    }
    /* Lahiri against the official tables of the Indian Astronomical
     * Ephemeris (true ayanamsa every three days, 0.1" rounding;
     * tests/lahiri_ref.h).  The 2026 edition uses IAU 2006 precession like
     * us; the 2019 edition still used the IAU 1976 rate, 0.30"/century
     * faster, which puts it 0.06" above ours in 2019. */
    {
        double worst[2] = {0.0, 0.0};
        for (const LahiriRef& r : kLahiriRef) {
            const int edition = r.jd_ut > 2460000.0; /* 0: 2019, 1: 2026 */
            const double tt = eph_ut_to_tt(r.jd_ut, edition ? 2026.5 : 2019.5);
            double dpsi, deps;
            eph_nutation(tt, &dpsi, &deps);
            const double d = std::fabs(
                (eph_ayanamsa(EPH_AYANAMSA_LAHIRI, tt) + dpsi) * R2D * 3600.0 - r.arcsec);
            worst[edition] = std::fmax(worst[edition], d);
        }
        std::printf("  Lahiri         max %.3f\" vs IAE 2026, %.3f\" vs IAE 2019\n",
                    worst[1], worst[0]);
        if (worst[1] > 0.06 || worst[0] > 0.13) {
            std::printf("FAIL: Lahiri vs the Indian Astronomical Ephemeris\n");
            fails++;
        }
    }
    if (std::fabs(eph_ayanamsa_lahiri(jd) - eph_ayanamsa(EPH_AYANAMSA_LAHIRI, jd)) > 0 ||
        !std::isnan(eph_ayanamsa((eph_ayanamsa_system)99, jd)) ||
        eph_ayanamsa_name((eph_ayanamsa_system)99) != nullptr ||
        std::fabs(eph_ayanamsa_custom(2451545.0, 0.4, 2451545.0) - 0.4) > 1e-15 ||
        std::fabs(eph_ayanamsa_custom(2415020.0, 0.0, 2451545.0) -
                  eph_ayanamsa(EPH_AYANAMSA_J1900, 2451545.0)) > 1e-8) {
        std::printf("FAIL: ayanamsa API\n");
        fails++;
    }

    /* eph_asc_mc against the reference angles. */
    for (const AngleRef& r : kAngleRef) {
        eph_observer ob = {r.lon * D2R, r.lat * D2R, 0.0};
        double asc, mc;
        eph_asc_mc(r.jd_ut, &ob, &asc, &mc);
        const double da = std::fabs(wrap180(asc * R2D - r.value[0])) * 3600.0;
        const double dm = std::fabs(wrap180(mc * R2D - r.value[1])) * 3600.0;
        if (da > 0.02 || dm > 0.02) {
            std::printf("FAIL: asc/mc off by %.4f\"/%.4f\" (jd %.2f)\n", da, dm, r.jd_ut);
            fails++;
        }
    }

    /* All systems vs the reference (tools/make_houses_reference.py, an
     * independent implementation): 1900-2030, |latitude| <= 65 degrees and
     * polar cases for the systems defined there. */
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
        std::printf("%-14s max %.4f\"\n", kNames[s], worst[s]);
        if (worst[s] > 0.02) {
            std::printf("FAIL: %s\n", kNames[s]);
            fails++;
        }
    }

    /* Angles and sensitive points vs the reference. */
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
    std::printf("Angles         max %.4f\"\n", worst_angle);
    if (worst_angle > 0.02) {
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
